using OnlinePCA
using DelimitedFiles

infile = ARGS[1]
out_cov = ARGS[2]
out_mean = ARGS[3]

@info "ooc_cov" input=infile
cov_mat, colmeanvec = OnlinePCA.ooc_cov(infile, "raw", 1.0f0, 300000, "sparse_bincoo")

mkpath(dirname(out_cov))
writedlm(out_cov, cov_mat, ',')
writedlm(out_mean, colmeanvec, ',')
@info "Done" cov=out_cov mean=out_mean
