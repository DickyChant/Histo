# Loading Packages
using Plots
using Random
using Distributions
using SpecialFunctions
using QuadGK
using StatsBase
using LinearAlgebra

# Function used
function gaussian_unnorm(x,mu = 0, sigma = 3)
    return exp(-(x-mu)^2/2/sigma^2)
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

function analysis(nlog)
    mu = 0
    sigma = 3
    n = 10^nlog
    dis = Normal(mu,sigma)
    x = rand(dis,n)
    
    center = [i for i in -10:10]
    edge = [i for i in -10:(20/21):10] 
    
    cenpdfs = [gaussian_unnorm(cen) for cen in center]
    cenpdfs./=sum(cenpdfs)
    
    hist = fit(Histogram,x,edge)
    
    norm = fit(Normal,x)
    mumle = norm.μ
    sigmamle = norm.σ
    
    cenpfit = [gaussian_unnorm(cen,mumle,sigmamle) for cen in center] 
    cenpfit./=sum(cenpfit)
    
    return mse(cenpdfs,cenpfit),mse(cenpdfs,hist.weights/n)
end
# Main analysis

nlogs = 1:8

mses = [analysis(i) for i in nlogs]

mse1 = [mse[1] for mse in mses]
mse2 = [mse[2] for mse in mses]

plot(nlogs,[mse1,mse2],yaxis=:log,label = ["Fitted density" "Filled histogram"])






