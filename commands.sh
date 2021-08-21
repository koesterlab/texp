target/release/bexp preprocess --kallisto-quants /vol/nano/pietsch-rna-rela-ependymomas/results/kallisto/{D1112-0,D1644-0}/abundance.h5 --sample-ids D1112-0 D1644-0 > /vol/nano/tmp/bexp-norm.mpk
target/release/bexp sample-expression --threads 32 --sample-id D1112-0 /vol/nano/tmp/bexp-norm.mpk --output /vol/nano/tmp/bexp-sample-exp-D1112-0
