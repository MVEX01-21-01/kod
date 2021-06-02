#!/bin/bash
zip -u9X "mvex01-21-01.zip" \
  README \
  *.r \
  batch.sh \
  report_out/* \
  experimental/hier2.R experimental/agnenv.r experimental/hier2.env.r
  #DATA_ENFS/* \
  #envelopes/envs499_K_REPR2.rds envelopes/envs499_single_REPR2.rds \
  #envelopes/csrenv* \
