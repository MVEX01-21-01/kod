export R_PROGRESSR_ENABLE=TRUE
Rscript --no-save --no-restore "batch_${1}_envelopes.r" "$2"
