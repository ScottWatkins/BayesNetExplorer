module BayesNetExplorer

println("Welcome to BayesNetExplorer v0.3.1\n")
println("This package is under active development.\nFeedback and bug reports are appreciated!\n[scott.watkins@utah.edu]\n")
println("Loading the package and all dependencies ...")

using CSV, DataFrames, DelimitedFiles
using StatsBase, Statistics, Random, Distributions, Combinatorics
using HypothesisTests
using Plots, GraphRecipes, RCall, Suppressor, Glob
using LinearAlgebra, Impute
using Printf
using OhMyREPL

println("Exporting user functions ...")

export CSV
export DataFrames
export format_file
export impute_dataframe
export bne
export RRcalculator
export plot_network
export feature_selector
export bnemle
export querywriter
export bnenbc
export bootstrapRRtable
export bnescan
export pairwise_correlations
export correlation_analyzer
export getRRtablerows

include("impute_dataframe.jl")
include("format_file.jl")
include("feature_selector.jl")
include("bne.jl")
include("run_bnlearn.jl")
include("pcor.jl")
include("ConProb.jl")
include("recoder.jl")
include("clean_dfvars_by_frequency.jl")
include("create_powerset.jl")
include("recompose_powerset.jl")
include("getcounts.jl")
include("get_markov_blanket.jl")
include("RRcalculator.jl")
include("bne_bootstrap.jl")
include("plot_network.jl")
include("bnemle.jl")
include("queryparser.jl")
include("printformat.jl")
include("moralizeDAG.jl")
include("querywriter.jl")
include("bnenbc.jl")
include("bootstrapRRtable.jl")
include("correlation_analyzer.jl")
include("pairwise_correlations.jl")
include("StatsFunctions.jl")
include("bnescan.jl")
include("getRRtablerows.jl")

end

