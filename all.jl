# Import Packages
using Plots
using Random
using Distributions
using SpecialFunctions
using QuadGK
using StatsBase
using LinearAlgebra
# Define Functions & Constants
mu = 0
sigma = 3

function gaussian_unnorm(x,mu = 0, sigma = 3)
    return exp(-(x-mu)^2/2/sigma^2)
end

function gaussian_used(x,mu = 0, sigma = 3)
    return 1/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/2/sigma^2)
end

function mse(x,y)
    nx = length(x)
    ny = length(y)
    if (nx!=ny) 
        println(nx,ny)
        return false
    end
    sq = sum((x-y).*(x-y))
    return sqrt(sq/nx)
end

function getHistWeights(x,edge)
    hist = fit(Histogram,x,edge)
    return hist.weights/sum(hist.weights)
end

function getNormWrong(x,center)
    norm = fit(Normal,x)
    mumle = norm.μ
    sigmamle = norm.σ 
    cenpfit = [gaussian_unnorm(cen,mumle,sigmamle) for cen in center] 
    cenpfit./=sum(cenpfit)
    return cenpfit
end

function getNormRight(x,center)
    norm = fit(Normal,x)
    mumle = norm.μ
    sigmamle = norm.σ 
    cenpfit = [gaussian_used(cen,mumle,sigmamle) for cen in center] 
    return cenpfit
end

function analysis(nlog,ntimes=200)
    n = 10^nlog
    dis = Normal(mu,sigma)

    center = [i for i in -10:10]
    edge = [i for i in -10:(20/21):10] 
    centerCorrect = (edge[2:end]+edge[1:end-1])/2

    probs = [quadgk(x -> gaussian_used(x), edge[i], edge[i+1], rtol=1e-8)[1] for i in 1:21]
    
    cenpdfsW = [gaussian_unnorm(cen) for cen in center]
    cenpdfsW./=sum(cenpdfsW)

    cenpdfsR = [gaussian_used(cen) for cen in center]
    
    mse1 = [0.0 for i = 1:ntimes] 
    mse2 = [0.0 for i = 1:ntimes]
    mse3 = [0.0 for i = 1:ntimes] 
    mse4 = [0.0 for i = 1:ntimes]
    for i = 1:ntimes
        x = rand(dis,n)
        
        histW = getHistWeights(x,edge)    
        cenpfitW = getNormWrong(x,center)
        cenpfitR = getNormRight(x,center)

        mse1[i] += mse(cenpdfsW,cenpfitW)
        mse2[i] += mse(cenpdfsW,histW)
        mse3[i] += mse(cenpdfsR,cenpfitR)
        mse4[i] += mse(probs,histW)
    end

    return mean(mse1),mean(mse2),mean(mse3),mean(mse4)
end
# Analyzing & Ploting
nlogs = 1:8

mses = [analysis(i) for i in nlogs]
mse1 = [mse[1] for mse in mses]
mse2 = [mse[2] for mse in mses]
mse3 = [mse[3] for mse in mses]
mse4 = [mse[4] for mse in mses]

plot(nlogs,[mse1,mse2,mse3,mse4],yaxis=:log,marker=:auto,label = ["Fitted density: Ziming","Filled histogram:Ziming","Fitted density: Corrected","Filled histogram: Corrected"])