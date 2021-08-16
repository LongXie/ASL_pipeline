#!/bin/bash
#$ -S /bin/bash
set -x -e

################################################
# Setup environment
# Software PATH
ANTSBIN=/home/longxie/pkg/antsbin/bin
C3DBIN=/home/longxie/pkg/c3d_tool/bin
FSLROOT=/share/apps/fsl/5.0.5
FSLBIN=$FSLROOT/bin
# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

###############################################
# Directories
ROOT=/home/longxie/WolkMCI/PMC_MCI/ASLPET
BRAINTEMPDIR=/home/longxie/WolkMCI/BrainTemplate/mni_icbm152_nlin_asym_09c
ANALYSISDIR=$ROOT/analysis_input
THRESHOLDMAPSDIR=/data/jag/fyang/PUBLIC/MCI

###############################################
# Parameters needs to be specify
# 0. Experiment information
expid=003
EXPDIR=$ROOT/exp/exp${expid}
#MAPSDIR=$EXPDIR/CBFmaps_includ
RECONDIR=$ROOT/recon/
OUTROOT=$EXPDIR/VoxelWise
DUMPDIR=$EXPDIR/dump_fsl

# parameters
ITER=10000
#RUNS=$(cat $ANALYSISDIR/run.txt)
#RUNS_ARRAY=($RUNS)
#ASLTYPES="abs"
#TIMEPOINTS=$(cat $ANALYSISDIR/timepoint.txt)
#MAPS=("NoThresh" "WithThresh")

###############################################
# main function
function main()
{

  # 0. clean up the directory
  reset_dir

  #############################################
  # 1. Superager cross group comparison
  if [[ 1 == 1 ]]; then
  ##################################
    #CompareGroupsStruct NC MCI
    CopyThicknessMaps

  fi
}

################################################
function reset_dir()
{
  mkdir -p $DUMPDIR
  rm -rf $DUMPDIR/*
}

################################################
function CompareGroupsStruct()
{
  grp1=$1
  grp2=$2
  grp1name=$grp1
  grp2name=$grp2

  # outdir
  OUTDIR=$OUTROOT/CompareGroupsStruct/${grp1name}_VS_${grp2name}/
  if [ -d $OUTDIR ]; then
    rm -rf $OUTDIR
  fi

  ##############################################
  # prepare the first group
  mkdir -p $OUTDIR/masks
  # load spreadsheet
  TXT=$ANALYSISDIR/Demog_included.csv
  IDS=($(cat $TXT | sed "s/,/ /g" | awk '{print $1}' | tail -n +2))
  GROUP=($(cat $TXT | sed "s/,/ /g" | awk '{print $4}' | tail -n +2))
  AGE=($(cat $TXT | sed "s/,/ /g" | awk '{print $13}' | tail -n +2))

  # copy data
  idx1=0
  idx2=0
  DEMOG1FN=$OUTDIR/Demog_group${grp1name}.txt
  DEMOG2FN=$OUTDIR/Demog_group${grp2name}.txt
  rm -rf $DEMOG1FN $DEMOG2FN
  for ((i=0;i<${#IDS[*]};i++)); do

    id=${IDS[i]}
    tp="Time1"
    testgrp=${GROUP[i]}
    age=${AGE[i]}

    if [[ $testgrp -ne $grp1 ]] && [[ $testgrp -ne $grp2 ]]; then
      echo "error $testgrp"
      exit
    fi

    mkdir -p $OUTDIR/Group${testgrp}

    if [[ $testgrp == $grp1 ]]; then
      prefix=$(printf %03d $idx1)
      idx1=$((idx1+1))
      echo "$age" >> $DEMOG1FN
    else
      prefix=$(printf %03d $idx2)
      idx2=$((idx2+1))
      echo "$age" >> $DEMOG2FN
    fi

    THKMAP=$RECONDIR/$id/$tp/MPRAGE/antsCorticalThickness/T1CorticalThicknessNormalizedToTemplate.nii.gz
    STHKMAP=$RECONDIR/$id/$tp/MPRAGE/antsCorticalThickness/T1CorticalThicknessNormalizedToTemplate_smooth3vox.nii.gz
    THKMASK=$RECONDIR/$id/$tp/MPRAGE/antsCorticalThickness/T1CorticalThicknessNormalizedToTemplateMask.nii.gz
    THKMASKLN=$OUTDIR/masks/${prefix}_${id}_${tp}_mask.nii.gz
    if [ ! -f $STHKMAP ]; then
      c3d $THKMAP -smooth 3vox \
        -o $STHKMAP
    fi
    ln -sf $STHKMAP \
      $OUTDIR/Group${testgrp}/${prefix}_${id}_thick_template.nii.gz
    if [ ! -f $THKMASK ]; then
      c3d $THKMAP -thresh 0.01 inf 1 0 \
        -smooth 3vox -thresh 0.5 inf 1 0 \
        -o $THKMASK
    fi
    ln -sf $THKMASK $THKMASKLN

  done

  # generate mask
  mask=$OUTDIR/tmpmask.nii.gz
  if [ ! -f $OUTDIR/NCmean.nii.gz ]; then
  c3d $OUTDIR/masks/*_1*_Time1*.nii.gz -mean -o $OUTDIR/NCmean.nii.gz
  fi
  if [ ! -f $OUTDIR/MCImean.nii.gz ]; then
  c3d $OUTDIR/masks/*_2*_Time1*.nii.gz -mean -o $OUTDIR/MCImean.nii.gz
  fi
  N_NC=$(ls $OUTDIR/masks/*_1*_Time1*.nii.gz | wc -w)
  N_MCI=$(ls $OUTDIR/masks/*_2*_Time1*.nii.gz | wc -w)
  c3d $OUTDIR/NCmean.nii.gz -scale $N_NC \
    $OUTDIR/MCImean.nii.gz -scale $N_MCI \
    -add -scale $(echo 1 / $((N_NC+N_MCI)) | bc -l) \
    -thresh 0.5 inf 1 0 \
    -o $mask

  # submit job to run
  #mask=$BRAINTEMPDIR/mni_icbm152_t1_tal_nlin_asym_09c_loresmask.nii
  OUTPRE=Struct_${grp1name}_${grp2name}
  qsubp2 -cwd -o $DUMPDIR -j y \
       -l h_vmem=25.1G,s_vmem=25G \
       -N "TwoTTest_${OUTPRE}" \
       $0 UnpairedTwoSampleTTest_sub \
       $OUTDIR/Group${grp1} $OUTDIR/Group${grp2} \
       $mask $OUTDIR $OUTPRE $DEMOG1FN $DEMOG2FN
}

function UnpairedTwoSampleTTest_sub()
{
  GP1DIR=$1
  GP2DIR=$2
  MASKFN=$3
  OUTROOT=$4
  OUTPRE=$5
  DEMOG1FN=$6
  DEMOG2FN=$7
  idx=0

  # count number of files for each group
  GP1_FILES=($(ls $GP1DIR/*.nii.gz))
  GP2_FILES=($(ls $GP2DIR/*.nii.gz))
  N_GP1=${#GP1_FILES[*]}
  N_GP2=${#GP2_FILES[*]}
  N=$((N_GP1+N_GP2))
  DEMOGINFO=$(cat $DEMOG1FN | head -n 1 | tail -n 1)
  N_DEMOG=${#DEMOGINFO[*]}

  # Combine to 4D image
  if [ ! -f $OUTROOT/cbf4D.nii.gz ]; then
  fslmerge -t \
    $TMPDIR/cbf4D.nii.gz \
    $GP1DIR/*.nii.gz \
    $GP2DIR/*.nii.gz
  fi

  # Generate design file
  DESIGNFNTMP=$TMPDIR/design.txt
  DESIGNFN=$OUTROOT/design.mat
  rm -rf $DESIGNFN
  #echo "/NumWaves 2" > $DESIGNFN
  #echo "/NumPoints $N" >> $DESIGNFN
  #echo "/PPheights 1 1" >> $DESIGNFN
  #echo "/Matrix" >> $DESIGNFN
  for ((iter=0;iter<${N_GP1};iter++)); do
    idx=$((idx+1))
    echo "1 0 $(cat $DEMOG1FN | head -n $((iter+1)) | tail -n 1)" >> $DESIGNFNTMP
    #echo "1 0" >> $DESIGNFNTMP
  done
  for ((iter=0;iter<${N_GP2};iter++)); do
    idx=$((idx+1))
    echo "0 1 $(cat $DEMOG2FN | head -n $((iter+1)) | tail -n 1)">> $DESIGNFNTMP
    #echo "0 1">> $DESIGNFNTMP
  done
  Text2Vest $DESIGNFNTMP $DESIGNFN

  # Generate contrast file
  CONTRASTFNTMP=$TMPDIR/constrast.txt
  CONTRASTFN=$OUTROOT/design.con
  rm -rf $CONTRASTFN
  #echo "/NumWaves 2" > $CONTRASTFN
  #echo "/NumContrasts 2" >> $CONTRASTFN
  #echo "/PPheights 1 1" >> $CONTRASTFN
  #echo "/Matrix" >> $CONTRASTFN
  ZEROS=""
  for ((i=0;i<$N_DEMOG;i++)); do
    ZEROS="$ZEROS 0"
  done
  echo "1 -1 $ZEROS" >> $CONTRASTFNTMP
  echo "-1 1 $ZEROS" >> $CONTRASTFNTMP
  Text2Vest $CONTRASTFNTMP $CONTRASTFN

  # Perform statistical test
  randomise \
    -i $TMPDIR/cbf4D.nii.gz \
    -o $OUTROOT/$OUTPRE \
    -d $DESIGNFN \
    -t $CONTRASTFN \
    -m $MASKFN \
    -x -n $ITER -T
}

################################################
function CopyThicknessMaps()
{
  # outdir
  OUTDIR=$OUTROOT/CompareGroupsStruct/ThicknessMaps/
  if [ -d $OUTDIR ]; then
    rm -rf $OUTDIR
  fi
  mkdir -p $OUTDIR/subject $OUTDIR/template

  ##############################################
  # load spreadsheet
  TXT=$ANALYSISDIR/Demog_included.csv
  IDS=($(cat $TXT | sed "s/,/ /g" | awk '{print $1}' | tail -n +2))
  GROUP=($(cat $TXT | sed "s/,/ /g" | awk '{print $4}' | tail -n +2))
  AGE=($(cat $TXT | sed "s/,/ /g" | awk '{print $13}' | tail -n +2))

  # copy data
  idx1=0
  idx2=0
  for ((i=0;i<${#IDS[*]};i++)); do

    id=${IDS[i]}
    tp="Time1"
    testgrp=${GROUP[i]}
    age=${AGE[i]}

    THKMAP=$RECONDIR/$id/$tp/MPRAGE/antsCorticalThickness/T1CorticalThickness.nii.gz
    THKMAPTEMPLATE=$RECONDIR/$id/$tp/MPRAGE/antsCorticalThickness/T1CorticalThicknessNormalizedToTemplate.nii.gz
    cp $THKMAP $OUTDIR/subject/${id}_${tp}_thickmap_subject.nii.gz
    cp $THKMAPTEMPLATE $OUTDIR/template/${id}_${tp}_thickmap_template.nii.gz

  done
}
   
##################################################
# Main entrypoint
cmd=$0
if [[ $# -lt 2 ]]; then

  main $@

else
  cmd=$1
  shift
  $cmd $@
fi
