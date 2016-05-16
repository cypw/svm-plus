#!/usr/bin/env sh

SVM_TRAIN=../svm-train
SVM_PREDICT=../svm-predict

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$SVM_TRAIN heart_scale model_default \
2>&1 | tee -a ./Log/train-default.log

$SVM_PREDICT heart_scale model_default output_default \
2>&1 | tee -a ./Log/pred-default.log


