module DeconvolutionTests

import DeconvolveDistribution
import DeconvolveDistribution: 
        decon, DeconvolutionMethod, 
        FourierDeconv, EfronDeconv
using HypothesisTests
using Distributions
using StatsBase: sample, midpoints

# include("deconvolve.jl")
include("deconv_boot_test.jl")
include("test_statistics.jl")
include("simulate.jl")

export deconv_boot_test, distance_test_statistic, KS_test_statistic

end # module
