function dW(du,u,p,t)
    @unpack τE,τI = p

    Threads.@threads for i = 1:N
          du[i+2N] =(1. /τE)*(0.001)
    
          du[i+3N] =(1. /τI)*(0.001)
    end
end


function F(z,η,Δ)
    return -0.5*im*((z-1)^2) + (im*η-Δ)*0.5*((z+1)^2)
end

function G(z,g,Vsyn)
    return im*Vsyn*g*0.5*((z+1)^2) - 0.5*g*(z^2-1)
end

function FR(z,τ)

    cz = conj(z)

    return (1. /(τ*π))*real((1. - cz)/(1. + cz))

end

function nextgen_stp_de(du,u,p,t)
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = p
        s=0.

        
        zE=u[1]
        zI=u[2]
        gEE=u[3]
        pEE=u[4]
        gIE=u[5]
        pIE=u[6]
        gEI=u[7]
        pEI=u[8]
        gII=u[9]
        pII=u[10]

        

        if  (10000 < t < 10100)
            s = AMP*sin(fr*t)
        else
            s = 0.0;
        end
        
        #zE
        du[1] = (1/τE)*(F(zE,η_0E+s,ΔE) + G(zE,gEE,VsynEE) + G(zE,gEI,VsynEI))
        #zI
        du[2] = (1/τI)*(F(zI,η_0I+s,ΔI) + G(zI,gIE,VsynIE) + G(zI,gII,VsynII))
        #gEE
        du[3] = αEE * (-gEE + pEE)
        #pEE 
        du[4] = αEE * (-pEE + κSEE*FR(zE,τE))
        #gIE
        du[5] = αIE * (-gIE + pIE)
        #pIE
        du[6] = αIE * (-pIE + κSIE * FR(zE,τE))
        #gEI
        du[7] = αEI * (-gEI + pEI)
        #pEI 
        du[8] = αEI * (-pEI + κSEI * FR(zI,τI))
        #gII
        du[9] = αII * (-gII + pII)
        #pII
        du[10] = αII * (-pII + κSII * FR(zI,τI))
end


@with_kw mutable struct NextGen2PopParams{R}
    ΔE::R = 0.5
    ΔI::R = 0.5
    η_0E::R = 20.
    η_0I::R =5.
    τE::R = 5.0
    τI::R = 5.0
    αEE::R = 1.
    αIE::R = 0.5
    αEI::R = 0.6
    αII::R = 2.0
    κSEE::R = 1.5
    κSIE::R = 1.5
    κSEI::R = 4.5
    κSII::R = 2.4
    κVEE::R = 0.
    κVIE::R = 0.
    κVEI::R = 0.
    κVII::R = 0.
    VsynEE::R = 12.
    VsynIE::R = 10.
    VsynEI::R = -10.
    VsynII::R = -10.
    κ::R = 0.3
end