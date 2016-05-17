#!/usr/bin/env sh

SVM_TRAIN=../../svm-train
SVM_PREDICT=../../svm-predict
SVM_PREDICT_EXTRA=../../libsvm-3.21/svm-predict

DATA=./data
MODEL=./model
LOG=./Log
OUTPUTS=./outputs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo "- - - - - - train"
$SVM_TRAIN $DATA/heart_scale $MODEL/model_default \
2>&1 | tee -a $LOG/train-default.log

echo "- - - - - - testing"
$SVM_PREDICT $DATA/heart_scale $MODEL/model_default $OUTPUTS/output_default \
2>&1 | tee -a $LOG/pred-default.log

echo "- - - - - - testing by libsvm"
$SVM_PREDICT_EXTRA $DATA/heart_scale $MODEL/model_default $OUTPUTS/output_default_extra \
2>&1 | tee -a $LOG/pred-default.log


