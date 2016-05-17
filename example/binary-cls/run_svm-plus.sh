#!/usr/bin/env sh


SVM_TRAIN=../../svm-train
SVM_PREDICT=../../svm-predict
SVM_PREDICT_EXTRA=../../libsvm-3.21/svm-predict

DATA=./data
MODEL=./model
LOG=./Log
OUTPUTS=./outputs

echo "- - - - - - training"
$SVM_TRAIN -s 5 -C 1 -f $DATA/heart_scale $DATA/heart_scale $MODEL/model_svm-plus \
2>&1 | tee -a ./Log/train-svm-plus.log

echo "- - - - - - testing"
$SVM_PREDICT $DATA/heart_scale $MODEL/model_default $OUTPUTS/output_svm-plus \
2>&1 | tee -a ./Log/pred-svm-plus.log

echo "- - - - - - testing by libsvm"
$SVM_PREDICT_EXTRA $DATA/heart_scale $MODEL/model_default $OUTPUTS/output_svm-plus \
2>&1 | tee -a ./Log/pred-svm-plus-extra.log
