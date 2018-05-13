module DeconvolutionTests

import DeconvolveDistribution
using HypothesisTests
using Distributions
using StatsBase: sample, midpoints

include("deconvolve.jl")
include("deconv_boot_test.jl")
include("test_statistics.jl")
include("simulate.jl")

end # module
