module VoronoiFVMDiffEq

using VoronoiFVM
using ExtendableGrids
using ExtendableSparse
using DocStringExtensions
using DiffResults
using ForwardDiff
using LinearAlgebra

import DifferentialEquations
using RecursiveArrayTools

import VoronoiFVM: eval_and_assemble, Node, BNode,
    _complete!, _fill!, _spec,
    _firstnodedof, _lastnodedof,
    AbstractSystem, num_species,isdata,nofunc2

include("vfvm_diffeq_interface.jl")

end
