using RCall JLD, Random, LinearAlgebra

include("Init2.jl")

R"source('processImages.r')"

@rget hold_images hold_labels
images = reshape(hold_images,(28,28,length(hold_labels)))

DSSModel = EM_Main_init(hold_images,1,1,3,10 ,1000,hold_labels)



