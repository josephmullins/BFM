include("../src/estimation/estimation.jl")

M = CSV.read("data/MarriagePanel.csv",DataFrame,missingstring="NA")
K = CSV.read("data/KidPanel.csv",DataFrame,missingstring="NA")
F = CSV.read("data/MarriageFile.csv",DataFrame,missingstring="NA")

cprobs = CSV.read("data/CustodyMomentsSimple.csv",DataFrame).p[:]
τgrid = [0.1,0.3,0.5,0.7,0.9]

stage1(M)
stage2(K)
stage3(K,τgrid,cprobs)


data_moms(F,M)

model_moms(rand(1:10,100),rand(1:10,100),rand(0:1,100).==1,rand(0:2,100))
