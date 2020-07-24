# Loading Packages
using Plots
using Random
using Distributions
using SpecialFunctions
using QuadGK
using StatsBase
using LinearAlgebra

# Function used
function gaussian_used(x,mu = 0, sigma = 3)
    return 1/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/2/sigma^2)
end

function mse(x,y)
    nx = length(x)
    ny = length(y)
    if (nx!=ny) 
        return false
    end
    sq = sum((x-y).*(x-y))
    return sqrt(sq/n)
end

function analysis(nlog)
    mu = 0
    sigma = 3
    n = 10^nlog
    dis = Normal(mu,sigma)
    x = rand(dis,n)
    
    edge = [i for i in -10:10]
    center = (edge[2:end]+edge[1:end-1])/2
    
    probs = [quadgk(x -> gaussian_used(x), edge[i], edge[i+1], rtol=1e-8)[1] for i in 1:20]
    cenpdfs = [gaussian_used(cen) for cen in center]
    
    hist = fit(Histogram,x,edge)
    
    norm = fit(Normal,x)
    mumle = norm.μ
    sigmamle = norm.σ
    
    cenpfit = [gaussian_used(cen,mumle,sigmamle) for cen in center] 
    
    return mse(cenpdfs,cenpfit),mse(probs,hist.weights/n)
end
# Main analysis

nlogs = 1:8

mses = [analysis(i) for i in nlogs]

mse1 = [mse[1] for mse in mses]
mse2 = [mse[2] for mse in mses]

plot(nlogs,[mse1,mse2],yaxis=:log,label = ["Fitted density" "Filled histogram"])






