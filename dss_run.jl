using RCall, JLD, Random, LinearAlgebra

include("Init2.jl")

R"source('processImages.r')"

@rget hold_images hold_labels
images = reshape(hold_images,(28,28,length(hold_labels)))

DSSModel = EM_Main_init(images,1,1,3,10 ,1000,hold_labels)

# Map Labels function 



# Run model 
# DSSModel = EM_Main_init(resizedImages,1,1,3,10 ,1000,uniqLabels)


# Get class predictions
# predProbs = DSSModel[end-4]
# @rput predProbs

