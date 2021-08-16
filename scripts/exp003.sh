#/bin/bash 
#$ -S /bin/bash
set -e

#set -x -e

##############################################
# Setup environment
# Software PATH
export ANTSPATH=/home/longxie/pkg/bin/antsbin/bin
ANTSCorticalThicknessPATH=/home/longxie/pkg/bin/antsbin/bin
#ANTSPATH=/data/picsl/longxie/pkg/antsbin/bin
#C3DPATH=/data/picsl/longxie/pkg/c3d_tool/bin
C3DPATH=/data/picsl/pauly/bin/c3d_affine_tool
FSLROOT=/data/picsl/longxie/pkg/fsl
FSLPATH=$FSLROOT/bin
JAVAEXEBIN=/data/picsl/longxie/CardiacMR/JavaExecutable
JAVABIN=/data/picsl/longxie/pkg/jdk1.7.0_71/bin
export ASHS_ROOT=/data/picsl/pauly/wolk/ashs
export PATH=$PATH:$ASHS_ROOT/bin
MATLABCODEDIR=/home/longxie/matlab
MATLAB_BIN=/share/apps/matlab/R2013a/bin/matlab
SRPATH=/home/longxie/pkg/PatchSuperResolution/build_release
export ASHS_ROOT=/data/picsl/pauly/wolk/ashs-fast
export PATH=$PATH:$ASHS_ROOT/bin
export ALOHA_ROOT=/data/picsl/longxie/pkg/aloha_LatestC3D


#################################################################
# INIT

# Directories
ALLDATADIR=/home/longxie/WolkMCI/PMC_MCI/recon
ROOT=/home/longxie/WolkMCI/PMC_MCI/ASLPET
ANALYSISDIR=$ROOT/analysis_input
INFODIR=$ROOT/info
DATADIR=$ROOT/recon
CODEDIR=$ROOT/code
ASLFUNCTIONDIR=/home/longxie/WolkMCI/code/ASL_function
OASISTEMPLATEDIR=/data/picsl/longxie/WolkMCI/BrainTemplate/OASIS

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

DEMOG=$ANALYSISDIR/Demographics_Combined_04192019.csv
TPTXT=$ANALYSISDIR/timepoint.txt
RUNTXT=$ANALYSISDIR/run.txt
SUBJTXT=$INFODIR/subj.txt
OUTDEMOG=$INFODIR/demog.csv


##################################################################
# Parameters needs to be specify
# 0. Experiment information
expid=003
EXPDIR=$ROOT/exp/exp${expid}
DUMPDIR=$EXPDIR/dump
LOGDIR=$EXPDIR/errorlog

##########################################
# Parameters
NOTES="-q all.q,basic.q"
QSUB4GOPT="$NOTES -l h_vmem=4.1G,s_vmem=4G"
QSUB5GOPT="$NOTES -l h_vmem=5.1G,s_vmem=5G"
QSUB6GOPT="$NOTES -l h_vmem=6.1G,s_vmem=6G"
QSUB8GOPT="$NOTES -l h_vmem=8.1G,s_vmem=8G"
QSUB10GOPT="$NOTES -l h_vmem=10.1G,s_vmem=10G"
QSUB12GOPT="$NOTES -l h_vmem=14.1G,s_vmem=12G"
QSUB14GOPT="$NOTES -l h_vmem=14.1G,s_vmem=14G"
QSUBWAITOPT="$NOTES -l h_vmem=2.1G,s_vmem=2G"
SKIP=0


##################################################################
function main()
{
  reset_dir

  # 0. make link to data
  #LinkData_Longi

  # 0. restore functional directory
  RestoreFunc

  # 1. motion correction using SPM8
  #MotionCorrection

  # 2. coregistration using SPM8
  #CoRegistration

  # 3. register to MNI using DARTEL
  #DARTEL_SPM8

  # 4. segment anatomical image
  #GenerateTissueMaps

  # 5. normalize structural images
  #NormMNI_structural_DARTEL

  # 6. smooth realigned functional image
  #SmoothFuncImage

  # 7. CBF quantification
  #CBFQuantification

  # 8. SCORE outlier cleaning
  #SCORECleaning

  # 9. compute global CBF
  #GlobalCBF

  # 10. compute relative CBF maps in subject space
  #RelativeSubjCBF 

  # 11. normalize to MNI space
  #NormMNI_DARTEL

  # 12. compute relative CBF maps in template space
  #RelativeTempCBF

  # 13. extract AAL ROI CBF
  #ExtractGMAALCBF


  ####################################
  # get hippocampus measurements
  ####################################
  # trim T1 image
  #TrimNeck

  # SR
  #SR

  # ASHST1
  #ASHST1
  #CleanASHST1

  # ALOHA
  #AlohaMTL
  

  # ASHSICV
  #ASHSICV
  #CleanASHSICV


  ####################################
  # Summarize ASL and hippo information
  ####################################
  #SummarizeHippo



  # 14. ants T1
  #antsT1

  exit
}

#############################################
function ReFormateDate()
{
  indate=$1

  if [[ $indate == "" || $indate == " " ]]; then

    outdate=$indate

  else

    DD=$(date -d "$indate" '+%d')
    MM=$(date -d "$indate" '+%m')
    YYYY=$(date -d "$indate" '+%Y')
    outdate="${YYYY}-${MM}-${DD}"

  fi

  echo $outdate
}

##################################################################
function LinkData_chooseASLPET()
{
  # go through all cases and find the qualify ones
  N=$(cat $DEMOG | wc -l)
  RUN=$(cat $RUNTXT)
  mkdir -p $INFODIR
  rm -f $SUBJTXT

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $DEMOG Study)
  TIMECol=$(csvcol.sh $DEMOG Time)
   
  # header
  ROW=$(cat -A $DEMOG | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  echo $ROW > $OUTDEMOG

  # go through each subjects
  PREFIX=LD
  for ((i=2;i<${N};i++)); do

    ROW=$(cat -A $DEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # select time points   
    study=$(echo $ROW | cut -f $STUDYCol -d ",")
    if [[ $study != "ASL-PET" ]]; then
      continue
    fi

    # link MPRAGE
    mkdir -p $DATADIR/$id/$tp/MPRAGE
    ln -sf $ALLDATADIR/$id/$tp/MPRAGE/* \
           $DATADIR/$id/$tp/MPRAGE/.

    # link ASL
    for run in $RUN; do
      if [[ -d $ALLDATADIR/$id/$tp/$run ]]; then
        mkdir -p $DATADIR/$id/$tp/$run
        ln -sf $ALLDATADIR/$id/$tp/$run/* \
               $DATADIR/$id/$tp/$run/.
      fi
    done

    echo $ROW >> $OUTDEMOG

  done
}

function LinkData_Longi()
{
  # go through all cases and find the qualify ones
  IDs=($(cat $DEMOG | awk -F, '{print $1}' | uniq))
  N=${#IDs[*]}
  RUN=$(cat $RUNTXT)
  mkdir -p $INFODIR

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $DEMOG Study)
  TIMECol=$(csvcol.sh $DEMOG Time)
  SCANDATECol=$(csvcol.sh $DEMOG MRI_Date)

  # header
  ROW=$(cat -A $DEMOG | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  echo "$ROW,AnalysisType,BaselineTimePoint,BaselineScanDate,DateDiffFromBaseline" > $OUTDEMOG

  # go through each subjects
  PREFIX=LD
  for ((i=0;i<${N};i++)); do

    # get id
    id=${IDs[$i]}

    # get cases and continue if less than 2
    set +e
    Ncases=$(cat $DEMOG | grep "^$id," | grep "ASL-PET" | wc -l)
    set -e
    if [[ $Ncases -lt 2 ]]; then
      continue
    fi

    # process each case
    TYPE=""
    BLTP=""
    BLSD=""
    cat -A $DEMOG | grep "^$id," | grep "ASL-PET" | \
      while read line; do

        line=$(echo $line | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
        tporig=$(echo $line | cut -f $TIMECol -d ",")
        tp=Time${tporig}
        scandate=$(echo $line | cut -f $SCANDATECol -d ",")
        scandate=$(ReFormateDate $scandate)

        # link ASL
        EXIST=0
        for run in $RUN; do
          if [[ -d $ALLDATADIR/$id/$tp/$run ]]; then
            EXIST=1
            mkdir -p $DATADIR/$id/$tp/$run
            ln -sf $ALLDATADIR/$id/$tp/$run/* \
                   $DATADIR/$id/$tp/$run/.
          fi
        done

        if [[ $EXIST == 1 ]]; then

          # link MPRAGE
          mkdir -p $DATADIR/$id/$tp/MPRAGE
          ln -sf $ALLDATADIR/$id/$tp/MPRAGE/* \
             $DATADIR/$id/$tp/MPRAGE/.

          if [[ $TYPE == "" ]]; then
            TYPE="Baseline"
            BLTP=$tporig
            BLSD=$scandate
            date_diff=0
          else
            TYPE="Longitudinal"
            date_diff=$(( \
              ($(date -d $scandate +%s) - \
              $(date -d $BLSD +%s) )/(60*60*24) ))
            date_diff=$(echo ${date_diff#-})
          fi

          echo "$line,$TYPE,$BLTP,$BLSD,$date_diff" >> $OUTDEMOG
        fi

      done
  done
}

##################################################################
function RestoreFunc()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=RF
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 RestoreFunc_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function RestoreFunc_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    restore_session('$DATADIR/$id/$tp','func','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function MotionCorrection()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=MC
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 MotionCorrection_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function MotionCorrection_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    mc_SPM8_session('$DATADIR/$id/$tp','RAWPCASL.nii.gz','mc','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function CoRegistration()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=CR
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -cwd -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 CoRegistration_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function CoRegistration_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    coregister_SPM8_session('$DATADIR/$id/$tp','$DATADIR/$id/$tp/MPRAGE/MPRAGE.nii.gz','mc/meanRAWPCASL.nii.gz','mc/rRAWPCASL.nii.gz','coreg','c','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function GenerateTissueMaps()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=GTM
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 GenerateTissueMaps_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function GenerateTissueMaps_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    generate_tissuemaps_session('$DATADIR/$id/$tp','$DATADIR/$id/$tp/MPRAGE/NewSegment','coreg/cmeanRAWPCASL.nii.gz','$DATADIR/$id/$tp/MPRAGE/segmentation','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function SmoothFuncImage()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=SFI
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 SmoothFuncImage_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function SmoothFuncImage_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    smooth_SPM8_session('$DATADIR/$id/$tp','coreg/crRAWPCASL.nii.gz','smooth','s','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function CBFQuantification()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=CQ
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 CBFQuantification_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function CBFQuantification_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    perfASL_SPM8_session('$DATADIR/$id/$tp','smooth/scrRAWPCASL.nii.gz','perf','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function SCORECleaning()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=SC
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 SCORECleaning_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function SCORECleaning_sub()
{
  id=$1
  tp=$2

  # Run matlab
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    OutlierCleaning_updatedSCORE_session('$DATADIR/$id/$tp','$DATADIR/$id/$tp/MPRAGE/segmentation','perf/cbf_0_scrRAWPCASL_4D.nii.gz','$DATADIR/$id/$tp/MPRAGE/segmentation/SmallBrainMask.nii.gz','clean_SCORE','cleaned_meanCBF.nii.gz','stat.mat','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function GlobalCBF()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=GC
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 GlobalCBF_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function GlobalCBF_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  #mkdir -p $TMPDIR/${id}
  #source $MATLABCODEDIR/matlab_batch.sh \
  #   $TMPDIR/${id}/ \
  #   $ASLFUNCTIONDIR \
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    globalCBF_session('$DATADIR/$id/$tp','$DATADIR/$id/$tp/MPRAGE/segmentation','clean_SCORE/cleaned_meanCBF.nii.gz','clean_SCORE','globalCBF.mat','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function RelativeSubjCBF()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=RSC
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 RelativeSubjCBF_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function RelativeSubjCBF_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    relativeCBF_session('$DATADIR/$id/$tp','clean_SCORE/cleaned_meanCBF.nii.gz','clean_SCORE/globalCBF.mat','clean_SCORE','rel_','GMWM','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function DARTEL_SPM8()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=NS
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 NewSegment_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
 
  # Build DARTEL template
  PREFIX=DARTEL
  qsubp3 -cwd -o $DUMPDIR -j y \
       $QSUB12GOPT \
       -N "${PREFIX}_" \
       $0 DARTEL_sub notmain
  sleep 0.1

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function NewSegment_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  mkdir -p $TMPDIR/${id}
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    newsegment_SPM8_session('$DATADIR/$id/$tp','$DATADIR/$id/$tp/MPRAGE/MPRAGE.nii.gz','NewSegment','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

function DARTEL_sub()
{
  # generate a list of unique subjects 
  cat $OUTDEMOG | grep -v ID \
    | awk -F, '{print $1}' | uniq  \
    > $SUBJTXT

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    SPM8_DARTEL_Template('$DATADIR','$SUBJTXT','$TPTXT','NewSegment','$LOGDIR/DARTEL.txt');
MATCODE
}

##################################################################
function NormMNI_DARTEL()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=NMD
  for ((i=2;i<=2;i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    TEMPLATE=$DATADIR/$id/$tp/MPRAGE/NewSegment/Template_6.nii
    TEMP2MNI=$DATADIR/$id/$tp/MPRAGE/NewSegment/Template_6_2mni.mat

    if [ ! -f $TEMP2MNI ]; then
    
      qsub -o $DUMPDIR -j y \
           $QSUB4GOPT \
           -N "${PREFIX}_start" \
           $0 NormMNI_DARTEL_sub \
           $id $tp $TEMPLATE
      sleep 0.1

    fi

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
  
  # go through each subjects
  PREFIX=NMD
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub  -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 NormMNI_DARTEL_sub $id $tp $TEMPLATE
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function NormMNI_DARTEL_sub()
{
  id=$1
  tp=$2
  template=$3

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    DARTEL_Norm_SPM8_session('$DATADIR/$id/$tp','$template','$DATADIR/$id/$tp/MPRAGE/NewSegment','clean_SCORE/cleaned_meanCBF.nii.gz','$ROOT/MNI/reslice_reference.nii.gz','normalized','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function NormMNI_structural_DARTEL()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=NSD
  for ((i=2;i<=2;i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    TEMPLATE=$DATADIR/$id/$tp/MPRAGE/NewSegment/Template_6.nii
    TEMP2MNI=$DATADIR/$id/$tp/MPRAGE/NewSegment/Template_6_2mni.mat

    if [ ! -f $TEMP2MNI ]; then

      qsub -o $DUMPDIR -j y \
           $QSUB4GOPT \
           -N "${PREFIX}_start" \
           $0 NormMNI_structural_DARTEL_sub \
           $id $tp $TEMPLATE
      sleep 0.1

    fi

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
     $QSUBWAITOPT \
     -hold_jid "${PREFIX}_*" -sync y -b y \
     sleep 1

  # go through each subjects
  PREFIX=NSD
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 NormMNI_structural_DARTEL_sub $id $tp $TEMPLATE
      sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function NormMNI_structural_DARTEL_sub()
{
  id=$1
  tp=$2
  template=$3

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/WolkMCI/code/ASL_function')

    DARTEL_Norm_SPM8_structural(...
      '$template', ...
      '$DATADIR/$id/$tp/MPRAGE/NewSegment/u_rc1MPRAGE_Template.nii', ...
      '$DATADIR/$id/$tp/MPRAGE/MPRAGE.nii.gz', ...
      '$ROOT/MNI/reslice_reference.nii.gz', ...
      '$DATADIR/$id/$tp/MPRAGE/normalized', ...
      '$LOGDIR/${id}_${tp}.txt');

    DARTEL_Norm_SPM8_structural(...
      '$template', ...
      '$DATADIR/$id/$tp/MPRAGE/NewSegment/u_rc1MPRAGE_Template.nii', ...
      '$DATADIR/$id/$tp/MPRAGE/NewSegment/c1MPRAGE.nii.gz', ...
      '$ROOT/MNI/reslice_reference.nii.gz', ...
      '$DATADIR/$id/$tp/MPRAGE/normalized', ...
      '$LOGDIR/${id}_${tp}.txt');

    DARTEL_Norm_SPM8_structural(...
      '$template', ...
      '$DATADIR/$id/$tp/MPRAGE/NewSegment/u_rc1MPRAGE_Template.nii', ...
      '$DATADIR/$id/$tp/MPRAGE/NewSegment/c2MPRAGE.nii.gz', ...
      '$ROOT/MNI/reslice_reference.nii.gz', ...
      '$DATADIR/$id/$tp/MPRAGE/normalized', ...
      '$LOGDIR/${id}_${tp}.txt');

    DARTEL_Norm_SPM8_structural(...
      '$template', ...
      '$DATADIR/$id/$tp/MPRAGE/NewSegment/u_rc1MPRAGE_Template.nii', ...
      '$DATADIR/$id/$tp/MPRAGE/NewSegment/c3MPRAGE.nii.gz', ...
      '$ROOT/MNI/reslice_reference.nii.gz', ...
      '$DATADIR/$id/$tp/MPRAGE/normalized', ...
      '$LOGDIR/${id}_${tp}.txt');

    DARTEL_Norm_SPM8_structural(...
      '$template', ...
      '$DATADIR/$id/$tp/MPRAGE/NewSegment/u_rc1MPRAGE_Template.nii', ...
      '$DATADIR/$id/$tp/MPRAGE/segmentation/SmallBrainMask.nii.gz', ...
      '$ROOT/MNI/reslice_reference.nii.gz', ...
      '$DATADIR/$id/$tp/MPRAGE/normalized', ...
      '$LOGDIR/${id}_${tp}.txt');

    exit;
MATCODE

  c3d $DATADIR/$id/$tp/MPRAGE/normalized/swSmallBrainMask.nii.gz -thresh 0.5 1 1 0 -o $DATADIR/$id/$tp/MPRAGE/normalized/swSmallBrainMask.nii.gz
}

##################################################################
function RelativeTempCBF()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=RTC
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 RelativeTempCBF_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function RelativeTempCBF_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    relativeCBF_session('$DATADIR/$id/$tp','normalized/swcleaned_meanCBF.nii.gz','clean_SCORE/globalCBF.mat','normalized','rel_','GMWM','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function ExtractGMAALCBF()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=EAC
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 ExtractGMAALCBF_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function ExtractGMAALCBF_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    extractGMAALROI_dataset_session('$DATADIR/$id/$tp','/data/picsl/longxie/WolkMCI/ROIs/AAL_ROI/reslice_ROI_MNI_V4.nii','normalized/swcleaned_meanCBF.nii.gz','normalized/rel_swcleaned_meanCBF.nii.gz','$DATADIR/$id/$tp/MPRAGE/normalized/swc1MPRAGE.nii.gz','AALROI','AALROI.mat','$LOGDIR/${id}_${tp}.txt');
MATCODE
}

##################################################################
function ExtractGMFrontalAALCBF()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=EAC
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 ExtractGMFrontalAALCBF_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # if there is error, stop
  if [ "$(ls $LOGDIR)" ]; then
    echo "Error occured. Check errorlog directory."
    exit
  fi
}

function ExtractGMFrontalAALCBF_sub()
{
  id=$1
  tp=$2

  # Run matlab to denoise T1 image
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('$ASLFUNCTIONDIR');
    extractGMFrontalAALROI_dataset_session('$DATADIR/$id/$tp','/data/picsl/longxie/WolkMCI/ROIs/AAL_ROI/reslice_ROI_MNI_V4.nii','normalized/swcleaned_meanCBF.nii.gz','normalized/rel_swcleaned_meanCBF.nii.gz','$DATADIR/$id/$tp/MPRAGE/normalized/swc1MPRAGE.nii.gz','FrontalAALROI','FrontalAALROI.mat','$LOGDIR/${id}_${tp}.txt');
MATCODE


}

########################################################
function antsT1()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=NS
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")

    if [[ $tp == 1 ]]; then

      # Submit jobs
      tp=Time${tp}
      qsub -o $DUMPDIR -j y \
           $QSUB12GOPT \
           -N "${PREFIX}_${id}_${tp}" \
           $0 antsT1_sub $id $tp
      sleep 0.1

    fi

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function antsT1_sub()
{
  id=$1
  tp=$2
  FILE=$DATADIR/$id/$tp/MPRAGE/MPRAGE.nii.gz
  OUT_DIR=$DATADIR/$id/$tp/MPRAGE/antsCorticalThickness
  OUT_PREFIX=$OUT_DIR/T1
  if [ -f $OUT_DIR/T1ACTStage5Complete.txt ]; then
    echo "already exist."
  else
    rm -rf $OUT_DIR
    mkdir -p $OUT_DIR

    # perform ANTS cortical thickness
    bash ${ANTSCorticalThicknessPATH}/antsCorticalThickness.sh -d 3 \
      -a $FILE \
      -e $OASISTEMPLATEDIR/T_template0.nii.gz \
      -m $OASISTEMPLATEDIR/T_template0_BrainCerebellumProbabilityMask.nii.gz \
      -f $OASISTEMPLATEDIR/T_template0_BrainCerebellumRegistrationMask.nii.gz \
      -p $OASISTEMPLATEDIR/Priors/priors%d.nii.gz \
      -t $OASISTEMPLATEDIR/T_template0_BrainCerebellum.nii.gz \
      -o $OUT_PREFIX \
      -k 1 \
      -q 1

  fi
}

##################################################################
function TrimNeck()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)
  #N=2

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=TN
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 TrimNeck_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function TrimNeck_sub()
{
  id=$1
  tp=$2
  ALLMPRAGEDIR=$ALLDATADIR/$id/$tp/MPRAGE
  ALLT1=$ALLMPRAGEDIR/MPRAGE.nii.gz
  ALLT1TRIM=$ALLMPRAGEDIR/MPRAGE_trim.nii.gz
  MPRAGEDIR=$DATADIR/$id/$tp/MPRAGE
  T1=$MPRAGEDIR/MPRAGE.nii.gz
  T1TRIM=$MPRAGEDIR/MPRAGE_trim.nii.gz

  if [[ ! -f $ALLT1TRIM ]]; then

    /data/picsl/pauly/bin/trim_neck_rf.sh \
      $ALLT1 $ALLT1TRIM

  fi

  rm -rf $T1TRIM
  if [[ -f $ALLT1TRIM ]]; then
    ln -sf $ALLT1TRIM $T1TRIM
  fi
}

##################################################################
function SR()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)
  #N=2

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=SR
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    # Submit jobs
    qsub -o $DUMPDIR -j y \
         $QSUB4GOPT \
         -N "${PREFIX}_${id}_${tp}" \
         $0 SR_sub $id $tp
    sleep 0.1

  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function SR_sub()
{
  id=$1
  tp=$2
  ALLMPRAGEDIR=$ALLDATADIR/$id/$tp/MPRAGE
  ALLT1=$ALLMPRAGEDIR/MPRAGE.nii.gz
  ALLT1TRIM=$ALLMPRAGEDIR/MPRAGE_trim.nii.gz
  ALLT1SR=$ALLMPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
  MPRAGEDIR=$DATADIR/$id/$tp/MPRAGE
  T1=$MPRAGEDIR/MPRAGE.nii.gz
  T1TRIM=$MPRAGEDIR/MPRAGE_trim.nii.gz
  T1SR=$MPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz

  #if [[ -f $ALLMPRAGEDIR/MPRAGE_trim_denosied_SR.nii.gz ]]; then
  #  mv $ALLMPRAGEDIR/MPRAGE_trim_denosied_SR.nii.gz \
  #     $ALLT1SR
  #fi

  if [[ ! -f $ALLT1SR ]]; then

    $SRPATH/NLMDenoise \
      -i $ALLT1TRIM \
      -o $TMPDIR/T1w_trim_denoised.nii.gz

    orient_code=$(c3d $T1TRIM -info | cut -d ';' -f 5 | cut -d ' ' -f 5)
    if [[ $orient_code == "Oblique," ]]; then
      orient_code=$(c3d $T1TRIM -info | cut -d ';' -f 5 | cut -d ' ' -f 8)
    fi

    c3d $TMPDIR/T1w_trim_denoised.nii.gz \
      -swapdim RPI \
      -o $TMPDIR/T1w_trim_denoised.nii.gz

    $SRPATH/NLMUpsample \
      -i $TMPDIR/T1w_trim_denoised.nii.gz \
      -o $TMPDIR/T1w_trim_denoised_SR.nii.gz \
      -lf 2 1 2

    c3d $TMPDIR/T1w_trim_denoised_SR.nii.gz\
      -swapdim $orient_code \
      -clip 0 inf \
      -o $ALLT1SR
  fi

  rm -f $T1SR
  if [[ -f $ALLT1SR ]]; then
    ln -sf $ALLT1SR $T1SR
  fi
}

##################################################################
function ASHST1()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)
  #N=2
  ATLAS=/home/longxie/ASHS_T1/ASHSexp/exp201/atlas/final
  Nrun=20

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=ASHST1
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    ALLMPRAGEDIR=$ALLDATADIR/$id/$tp/MPRAGE
    ALLT1=$ALLMPRAGEDIR/MPRAGE.nii.gz
    ALLT1TRIM=$ALLMPRAGEDIR/MPRAGE_trim.nii.gz
    ALLT1SR=$ALLMPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ALLASHSOUTDIR=$ALLMPRAGEDIR/ASHST1
    MPRAGEDIR=$DATADIR/$id/$tp/MPRAGE
    T1=$MPRAGEDIR/MPRAGE.nii.gz
    T1TRIM=$MPRAGEDIR/MPRAGE_trim.nii.gz
    T1SR=$MPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ASHSOUTDIR=$MPRAGEDIR/ASHST1

    if [[ -f $ALLT1 && -f $ALLT1SR && \
       ! -f $ALLASHSOUTDIR/final/${id}_${tp}_left_heur_volumes.txt ]]; then

      mkdir -p $ALLASHSOUTDIR
      $ASHS_ROOT/bin/ashs_main.sh \
        -a $ATLAS -d -T -I ${id}_${tp} \
        -g $ALLT1 -f $ALLT1SR \
        -Q -s 1-7 \
        -z $CODEDIR/scripts/ashs-fast-z.sh \
        -m $CODEDIR/scripts/identity.mat -M \
        -w $ALLASHSOUTDIR &
      sleep 0.1

      ALLJOBS=($(jobs -p))
      while [ ${#ALLJOBS[*]} -ge $Nrun ]; do
        sleep 60
        ALLJOBS=($(jobs -p))
      done
    fi

  done
}

##################################################################
function CleanASHST1()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)
  #N=2

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)

  # go through each subjects
  PREFIX=ASHST1
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}

    ALLMPRAGEDIR=$ALLDATADIR/$id/$tp/MPRAGE
    ALLT1=$ALLMPRAGEDIR/MPRAGE.nii.gz
    ALLT1TRIM=$ALLMPRAGEDIR/MPRAGE_trim.nii.gz
    ALLT1SR=$ALLMPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ALLASHSOUTDIR=$ALLMPRAGEDIR/ASHST1
    MPRAGEDIR=$DATADIR/$id/$tp/MPRAGE
    T1=$MPRAGEDIR/MPRAGE.nii.gz
    T1TRIM=$MPRAGEDIR/MPRAGE_trim.nii.gz
    T1SR=$MPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ASHSOUTDIR=$MPRAGEDIR/ASHST1

    if [[ -f $ALLASHSOUTDIR/final/${id}_${tp}_left_heur_volumes.txt && -f $ALLASHSOUTDIR/tse.nii.gz ]]; then

      rm -rf $ALLASHSOUTDIR/affine_t1_to_template \
        $ALLASHSOUTDIR/ants_t1_to_temp \
        $ALLASHSOUTDIR/bootstrap \
        $ALLASHSOUTDIR/dump \
        $ALLASHSOUTDIR/multiatlas \
        $ALLASHSOUTDIR/mprage* \
        $ALLASHSOUTDIR/tse.nii.gz \
        $ALLASHSOUTDIR/tse_to_chunktemp*.nii.gz

      if [[ -d $ASHSOUTDIR ]]; then
        rm -rf $ASHSOUTDIR
      fi
      ln -sf $ALLASHSOUTDIR $ASHSOUTDIR

    fi
  done
}

##################################################################
function ASHSICV()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)
  #N=2
  ICVATLAS=/data/jet/vpiskin/ASHS/build_ICV_atlas_trimmed/ICV_atlas/final/
  Nrun=20

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)
  TYPECol=$(csvcol.sh $OUTDEMOG AnalysisType)

  # go through each subjects
  PREFIX=ASHST1
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    if [[ $Type != "Baseline" ]]; then
      continue
    fi

    ALLMPRAGEDIR=$ALLDATADIR/$id/$tp/MPRAGE
    ALLT1=$ALLMPRAGEDIR/MPRAGE.nii.gz
    ALLT1TRIM=$ALLMPRAGEDIR/MPRAGE_trim.nii.gz
    ALLT1SR=$ALLMPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ALLASHSOUTDIR=$ALLMPRAGEDIR/ASHSICV
    MPRAGEDIR=$DATADIR/$id/$tp/MPRAGE
    T1=$MPRAGEDIR/MPRAGE.nii.gz
    T1TRIM=$MPRAGEDIR/MPRAGE_trim.nii.gz
    T1SR=$MPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ASHSOUTDIR=$MPRAGEDIR/ASHSICV

    if [[ -f $ALLT1 &&  \
       ! -f $ALLASHSOUTDIR/final/${id}_${tp}_left_corr_nogray_volumes.txt ]]; then

      mkdir -p $ALLASHSOUTDIR
      $ASHS_ROOT/bin/ashs_main.sh \
        -a $ICVATLAS -d -T -I ${id}_${tp} \
        -g $ALLT1 -f $ALLT1 \
        -Q -s 1-7 \
        -z $CODEDIR/scripts/ashs-fast-z.sh \
        -m $CODEDIR/scripts/identity.mat -M \
        -B \
        -w $ALLASHSOUTDIR &
      sleep 0.1

      ALLJOBS=($(jobs -p))
      while [ ${#ALLJOBS[*]} -ge $Nrun ]; do
        sleep 60
        ALLJOBS=($(jobs -p))
      done
    fi

  done
}

##################################################################
function CleanASHSICV()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)
  #N=2

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)
  TYPECol=$(csvcol.sh $OUTDEMOG AnalysisType)

  # go through each subjects
  PREFIX=ASHST1
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    if [[ $Type != "Baseline" ]]; then
      continue
    fi

    ALLMPRAGEDIR=$ALLDATADIR/$id/$tp/MPRAGE
    ALLT1=$ALLMPRAGEDIR/MPRAGE.nii.gz
    ALLT1TRIM=$ALLMPRAGEDIR/MPRAGE_trim.nii.gz
    ALLT1SR=$ALLMPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ALLASHSOUTDIR=$ALLMPRAGEDIR/ASHSICV
    MPRAGEDIR=$DATADIR/$id/$tp/MPRAGE
    T1=$MPRAGEDIR/MPRAGE.nii.gz
    T1TRIM=$MPRAGEDIR/MPRAGE_trim.nii.gz
    T1SR=$MPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ASHSOUTDIR=$MPRAGEDIR/ASHSICV

    if [[ -f $ALLASHSOUTDIR/final/${id}_${tp}_left_corr_nogray_volumes.txt && -f $ALLASHSOUTDIR/tse.nii.gz ]]; then

      rm -rf $ALLASHSOUTDIR/affine_t1_to_template \
        $ALLASHSOUTDIR/ants_t1_to_temp \
        $ALLASHSOUTDIR/bootstrap \
        $ALLASHSOUTDIR/dump \
        $ALLASHSOUTDIR/multiatlas \
        $ALLASHSOUTDIR/mprage* \
        $ALLASHSOUTDIR/tse.nii.gz \
        $ALLASHSOUTDIR/tse_to_chunktemp*.nii.gz \
        $ALLASHSOUTDIR/flirt_t2_to_t1 \
        $ALLASHSOUTDIR/tse_native_chunk*.nii.gz

      if [[ -d $ASHSOUTDIR ]]; then
        rm -rf $ASHSOUTDIR
      fi
      ln -sf $ALLASHSOUTDIR $ASHSOUTDIR

    fi
  done
}


##################################################################
function AlohaMTL()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)
  #N=3

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)
  TYPECol=$(csvcol.sh $OUTDEMOG AnalysisType)
  BLTPCol=$(csvcol.sh $OUTDEMOG BaselineTimePoint)

  # go through each subjects
  PREFIX=ALOHA
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}
    bltp=$(echo $ROW | cut -f $BLTPCol -d ",")
    bltp=Time${bltp}
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    if [[ $Type != "Longitudinal" ]]; then
      continue
    fi

    # baseline files 
    ALLBLMPRAGEDIR=$ALLDATADIR/$id/$bltp/MPRAGE
    ALLBLT1=$ALLBLMPRAGEDIR/MPRAGE.nii.gz
    ALLBLT1TRIM=$ALLBLMPRAGEDIR/MPRAGE_trim.nii.gz
    ALLBLT1SR=$ALLBLMPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ALLBLASHSOUTDIR=$ALLBLMPRAGEDIR/ASHST1
    ALLBLASHST1LSEG=$ALLBLASHSOUTDIR/final/${id}_${bltp}_left_lfseg_heur.nii.gz
    ALLBLASHST1RSEG=$ALLBLASHSOUTDIR/final/${id}_${bltp}_right_lfseg_heur.nii.gz

    # follow up files
    ALLMPRAGEDIR=$ALLDATADIR/$id/$tp/MPRAGE
    ALLFUT1=$ALLMPRAGEDIR/MPRAGE.nii.gz
    ALLFUT1TRIM=$ALLMPRAGEDIR/MPRAGE_trim.nii.gz
    ALLFUT1SR=$ALLMPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ALLALOHADIR=$ALLMPRAGEDIR/AlohaMTL_ASHST1_${bltp}
    MPRAGEDIR=$DATADIR/$id/$tp/MPRAGE
    ALOHADIR=$MPRAGEDIR/AlohaMTL_ASHST1_${bltp}

    if [[ -f $ALLBLT1SR && -f $ALLFUT1SR && -f $ALLBLASHST1LSEG && -f $ALLBLASHST1RSEG && ! -d $ALLALOHADIR/results_PHC ]]; then

      rm -rf $ALLALOHADIR
      qsub -o $DUMPDIR -j y \
           $QSUB5GOPT \
           -N "${PREFIX}_${id}_${tp}" \
           $0 AlohaMTL_sub \
                $ALLBLT1SR \
                $ALLFUT1SR \
                $ALLBLASHST1LSEG \
                $ALLBLASHST1RSEG \
                $ALLALOHADIR \
                $ALOHADIR
      sleep 0.1

    fi
  done

  # Wait for completion
  qsub -o $DUMPDIR -j y \
       $QSUBWAITOPT \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function AlohaMTL_sub()
{
  BLT1SRIMG=$1
  YT1SRIMG=$2
  BLASHST1LSEG=$3
  BLASHST1RSEG=$4
  SUBJALOHADIR=$5
  FINALALOHADIR=$6
  mkdir -p $SUBJALOHADIR

  if [[ ! -f $SUBJALOHADIR/results/volumes_left.txt || ! -f $SUBJALOHADIR/results/volumes_right.txt ]]; then

  c3d $BLT1SRIMG -swapdim RPI \
     -o $TMPDIR/baseline_RPI.nii.gz
  BLT1SRIMG=$TMPDIR/baseline_RPI.nii.gz

  c3d $YT1SRIMG -swapdim RPI \
     -o $TMPDIR/followup_RPI.nii.gz
  YT1SRIMG=$TMPDIR/followup_RPI.nii.gz

  c3d $BLT1SRIMG \
     $BLASHST1LSEG \
     -replace 1 999 2 999 10 999 11 999 12 999 13 999 \
     -thresh 999 999 1 0 \
     -int 0 -reslice-identity \
     -o $SUBJALOHADIR/bl_seg_left.nii.gz

  c3d $BLT1SRIMG \
     $BLASHST1RSEG \
     -replace 1 999 2 999 10 999 11 999 12 999 13 999 \
     -thresh 999 999 1 0 \
     -int 0 -reslice-identity \
     -o $SUBJALOHADIR/bl_seg_right.nii.gz

  $ALOHA_ROOT/scripts/aloha_main.sh \
    -b $BLT1SRIMG \
    -f $YT1SRIMG \
    -r $SUBJALOHADIR/bl_seg_left.nii.gz \
    -s $SUBJALOHADIR/bl_seg_right.nii.gz \
    -w $SUBJALOHADIR \
    -t 1-4

  fi

  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do

    if [[ ! -f $SUBJALOHADIR/results_${sub}/volumes_left.txt || ! -f $SUBJALOHADIR/results_${sub}/volumes_right.txt ]]; then

      if [[ $sub == "ERC" ]]; then
        label=(10 10)
      elif [[ $sub == "BA35" ]]; then
        label=(11 11)
      elif [[ $sub == "BA36" ]]; then
        label=(12 12)
      elif [[ $sub == "PHC" ]]; then
        label=(13 13)
      elif [[ $sub == "AHippo" ]]; then
        label=(1 1)
      elif [[ $sub == "PHippo" ]]; then
        label=(2 2)
      elif [[ $sub == "Hippo" ]]; then
        label=(1 2)
      fi

      SUBJALOHASUBDIR=$TMPDIR/aloha_${sub}
      mkdir -p $SUBJALOHASUBDIR/results
      ln -sf $SUBJALOHADIR/deformable $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/init $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/global $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/dump $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/final $SUBJALOHASUBDIR

      c3d $BLT1SRIMG \
        $BLASHST1LSEG \
        -thresh ${label[0]} ${label[1]} 1 0 \
        -int 0 -reslice-identity \
        -o $SUBJALOHASUBDIR/bl_seg_left.nii.gz

      c3d $BLT1SRIMG \
        $BLASHST1RSEG \
        -thresh ${label[0]} ${label[1]} 1 0 \
        -int 0 -reslice-identity \
        -o $SUBJALOHASUBDIR/bl_seg_right.nii.gz

      $ALOHA_ROOT/scripts/aloha_main.sh \
        -b $BLT1SRIMG \
        -f $YT1SRIMG \
        -r $SUBJALOHASUBDIR/bl_seg_left.nii.gz \
        -s $SUBJALOHASUBDIR/bl_seg_right.nii.gz \
        -w $SUBJALOHASUBDIR \
        -t 4

      mv $SUBJALOHASUBDIR/results \
         $SUBJALOHADIR/results_${sub}

    fi

    # measure thickness
    for side in left right; do

      if [[ ! -f $SUBJALOHADIR/results_${sub}/mean_thickness_${side}.txt || ! -f $SUBJALOHADIR/results_${sub}/median_thickness_${side}.txt ]]; then

        # measure thickness of bl mesh
        cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
          -T $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_thickmap.vtk \
          -p 1.2 -e 6 \
          $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}.vtk \
          $SUBJALOHADIR/results_${sub}/skel_blmptrim_seg_${side}.vtk

        # measure thickness of fu mesh
        cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
          -T $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om_thickmap.vtk \
          -p 1.2 -e 6 \
          $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om.vtk \
          $SUBJALOHADIR/results_${sub}/skel_blmptrim_seg_${side}_warped_to_futrim_om.vtk

        # get thickness
        $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
        addpath('/home/longxie/ASHS_T1/Application/TAUPET/longitudinal/code');
        MeasureMeanMedianThickness('$SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_thickmap.vtk','$SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om_thickmap.vtk','$SUBJALOHADIR/results_${sub}/mean_thickness_${side}.txt','$SUBJALOHADIR/results_${sub}/median_thickness_${side}.txt');
MATCODE

     fi
   done
  done

  rm -rf $SUBJALOHADIR/final $SUBJALOHADIR/global \
         $SUBJALOHADIR/init \
         $SUBJALOHADIR/dump  $SUBJALOHADIR/bl_seg_*.nii.gz

  rm -rf $FINALALOHADIR
  ln -sf $SUBJALOHADIR $FINALALOHADIR

}   

##################################################################
function SummarizeHippo()
{
  # go through all cases and find the qualify ones
  N=$(cat $OUTDEMOG | wc -l)

  # columns
  IDCol=1
  STUDYCol=$(csvcol.sh $OUTDEMOG Study)
  TIMECol=$(csvcol.sh $OUTDEMOG Time)
  TYPECol=$(csvcol.sh $OUTDEMOG AnalysisType)
  BLTPCol=$(csvcol.sh $OUTDEMOG BaselineTimePoint)
  DATEDIFFCol=$(csvcol.sh $OUTDEMOG DateDiffFromBaseline)

  # header
  HIPPOCSV=$INFODIR/structural_longitudinal.csv
  header=$(cat -A $OUTDEMOG | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header="$header,ALOHA_Success,ICV_BL"
  for side in L R; do
  for sub in AHippo PHippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_ASHST1VOL_BL"
  done
  done
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_PercChangeAnnualized"
  done
  done
  done
  header="$header$longheader,LONGIEND"
  echo "$header" > $HIPPOCSV


  # go through each subjects
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $OUTDEMOG | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
    id=$(echo $ROW | cut -f $IDCol -d ",")
    tp=$(echo $ROW | cut -f $TIMECol -d ",")
    tp=Time${tp}
    bltp=$(echo $ROW | cut -f $BLTPCol -d ",")
    bltp=Time${bltp}
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    date_diff=$(echo $ROW | cut -f $DATEDIFFCol -d ",")

    MPRAGEDIR=$DATADIR/$id/$tp/MPRAGE
    FUT1=$MPRAGEDIR/MPRAGE.nii.gz
    FUT1TRIM=$MPRAGEDIR/MPRAGE_trim.nii.gz
    FUT1SR=$MPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    ASHST1OUTDIR=$MPRAGEDIR/ASHST1
    ASHSICVOUTDIR=$MPRAGEDIR/ASHSICV
    ALOHADIR=$MPRAGEDIR/AlohaMTL_ASHST1_${bltp}

    # baseline info
    BLMPRAGEDIR=$DATADIR/$id/$bltp/MPRAGE
    BLT1=$BLMPRAGEDIR/MPRAGE.nii.gz
    BLT1TRIM=$BLMPRAGEDIR/MPRAGE_trim.nii.gz
    BLT1SR=$BLMPRAGEDIR/MPRAGE_trim_denoised_SR.nii.gz
    BLASHSICVOUTDIR=$BLMPRAGEDIR/ASHSICV
    BLASHST1OUTDIR=$BLMPRAGEDIR/ASHST1
    BLASHST1LSEG=$BLASHST1OUTDIR/final/${id}_${bltp}_left_lfseg_heur.nii.gz
    BLASHST1RSEG=$BLASHST1OUTDIR/final/${id}_${bltp}_right_lfseg_heur.nii.gz

    # check if exist
    EXIST=0
    if [[ $Type != "Longitudinal" ]]; then
      EXIST="-1"
    elif [[ -f $FUT1SR && -f $BLT1SR && -f $BLASHST1LSEG && -f $BLASHST1RSEG ]]; then
      EXIST=1
    else
      EXIST=0
    fi

    # get baseline ICV
    ICV=$(cat $BLASHSICVOUTDIR/final/${id}_${bltp}_left_corr_nogray_volumes.txt | awk -F' ' '{ print $5 }')

    # get baseline ASHST1 measurements
    ASHST1VOL=""
    for side in left right; do
      VOLFILE=$BLASHST1OUTDIR/final/${id}_${bltp}_${side}_heur_volumes.txt
      ASHST1VOL="$ASHST1VOL,$(cat $VOLFILE | grep Anterior_hippocampus | cut -d ' ' -f 5)"
      ASHST1VOL="$ASHST1VOL,$(cat $VOLFILE | grep Posterior_hippocampus | cut -d ' ' -f 5)"
      ASHST1VOL="$ASHST1VOL,$(cat $VOLFILE | grep ERC | cut -d ' ' -f 5)"
      ASHST1VOL="$ASHST1VOL,$(cat $VOLFILE | grep Br35 | cut -d ' ' -f 5)"
      ASHST1VOL="$ASHST1VOL,$(cat $VOLFILE | grep Br36 | cut -d ' ' -f 5)"
      ASHST1VOL="$ASHST1VOL,$(cat $VOLFILE | grep " PHC" | cut -d ' ' -f 5)"
    done

    # go through all the fields
    set +e
    RAWMEASURE=""
    LONGIMEASURE=""
    for type in volumes mean_thickness median_thickness; do
    for side in left right; do
      for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
      SUBJALOHASUBDIR=$ALOHADIR/results_${sub}
        if [[ $EXIST == 1 ]]; then
          # get volume and thickness
          if [[ -f $SUBJALOHASUBDIR/${type}_${side}.txt ]]; then
            MEA=$(cat $SUBJALOHASUBDIR/${type}_${side}.txt)
            MEA="${MEA// /,}"
            BL=$(echo $MEA | cut -d , -f 1)
            FU=$(echo $MEA | cut -d , -f 2)
            CHANGE=$(echo "scale=7;($BL-$FU)*100*365/$BL/$date_diff" | bc)
            LONGIMEASURE="$LONGIMEASURE,$CHANGE"

          else
            MEA=","
            EXIST=0
            LONGIMEASURE="$LONGIMEASURE,"
          fi
          RAWMEASURE="$RAWMEASURE,$MEA"
        else
          RAWMEASURE="$RAWMEASURE,,"
          LONGIMEASURE="$LONGIMEASURE,"
        fi
      done
    done
    done
     
    echo "$ROW,$EXIST,$ICV$ASHST1VOL$RAWMEASURE$LONGIMEASURE,-1" >> $HIPPOCSV

    set -e

  done
}


########################################################
function reset_dir()
{
  rm -rf $DUMPDIR 
  mkdir -p $DUMPDIR
  if [ -d $LOGDIR ]; then
    rm -r $LOGDIR
  fi
  mkdir -p $LOGDIR
}

########################################################
if [[ $# -lt 2 ]]; then

  main

else

  cmd=$1
  shift
  $cmd $@

fi
