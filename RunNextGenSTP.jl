using DifferentialEquations,Plots,Parameters,LinearAlgebra,JLD,Statistics
include("functions.jl")
include("GlobalFunctions/_includes.jl")


const SC,dist,lags,N,FC,missingROIs = networksetup(13,0.1;digits=10,nSC=2,nFC=1,N=140)
if ~@isdefined meanFC
meanFC,missingROI = get_mean_all_functional_data(;ROI=140,type="control")
end
NGp = NextGen2PopParams()

NGp.η_0I = -13.2
NGp.κ=0.0203

tspan = (0.,30000.)
const AMP=20
const fr = 10.
const W = SC

u0 = rand(N*12)
#u0[39] = 100.

p = NGp
prob = ODEProblem(nextgen_stp_de,u0,tspan,p)
sol = solve(prob,BS3(),maxiters=10e20,saveat=0:0.1:tspan[2])

R_c = zeros(N,N)
R_s = zeros(N,N)

r_c = sol[1:N,20000:90000]
r_s = sol[1:N,150000:300000]
for i = 1:N
    for j = i+1:N
        R_c[i,j] = cor(r_c[i,:],r_c[j,:])
        R_c[j,i] = R_c[i,j]

        R_s[i,j] = cor(r_s[i,:],r_s[j,:])
        R_s[j,i] = R_s[i,j]
    end
end
print("fit is: ", fit_r(R_c.^2,meanFC.^2))

pc = heatmap(R_c.^2)
ps = heatmap(R_s.^2)

plot(pc,ps)


