function dW(du,u,p,t)
    @unpack τE,τI = p

    Threads.@threads for i = 1:N
          du[i+2N] =(1. /τE)*(0.001)
    
          du[i+3N] =(1. /τI)*(0.001)
    end
end



function nextgen_stp_de(du,u,p,t)
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = p
 
    Threads.@threads for i = 1:N



        d=0.
        s=0.

        
        rE=u[i]
        rI=u[i+N]
        vE=u[i+2N]
        vI=u[i+3N]
        gEE=u[i+4N]
        pEE=u[i+5N]
        gIE=u[i+6N]
        pIE=u[i+7N]
        gEI=u[i+8N]
        pEI=u[i+9N]
        gII=u[i+10N]
        pII=u[i+11N]

        if i == 39 && (10000 < t < 10100)
            s = AMP*sin(fr*t)
        else
            s = 0.0;
        end


        for j = 1:N
            d += W[i,j]*u[j]
        end
        
        #rE
        du[i] =(1. /τE)*(-gEE*rE -gEI*rE - κVEE*rE - κVEI*rE +2. * rE * vE + (ΔE / (τE*pi)))
        #rI
        du[i+N] =(1. /τI)*(-gIE*rI - gII*rI -κVIE*rI - κVII*rI + 2. * rI * vI + (ΔI / (τI*pi)))
        #vE
        du[i+2N] =(1. /τE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vE) + κVEI*(vI - vE) - (τE^2)*(pi^2) * (rE^2.) +  vE^2. + η_0E + s )
        #vI
        du[i+3N] =(1. /τI)*(gIE*(VsynIE - vI) + gII*(VsynII - vI) + κVIE*(vE - vI) - (τI^2)*(pi^2)* (rI^2.) + vI^2. + η_0I + s )
        #gEE
        du[i+4N] = αEE * (-gEE + pEE)
        #pEE 
        du[i+5N] = αEE * (-pEE + κ*d + κSEE*rE)
        #gIE
        du[i+6N] = αIE * (-gIE + pIE)
        #pIE
        du[i+7N] = αIE * (-pIE + κSIE * rE)
        #gEI
        du[i+8N] = αEI * (-gEI + pEI)
        #pEI 
        du[i+9N] = αEI * (-pEI + κSEI * rI)
        #gII
        du[i+10N] = αII * (-gII + pII)
        #pII
        du[i+11N] = αII * (-pII + κSII * rI)

    end
end


@with_kw mutable struct NextGen2PopParams{R}
    ΔE::R = 0.5
    ΔI::R = 0.5
    η_0E::R = 20.
    η_0I::R =-20.
    τE::R = 1.0
    τI::R = 1.0
    αEE::R = 1.
    αIE::R = 1.4
    αEI::R = 0.7
    αII::R = 0.4
    κSEE::R = 1.5
    κSIE::R = 1.
    κSEI::R = 2.
    κSII::R = 3.
    κVEE::R = 0.
    κVIE::R = 0.
    κVEI::R = 0.
    κVII::R = 0.
    VsynEE::R = 10.
    VsynIE::R = 8.
    VsynEI::R = -8.
    VsynII::R = -12.
    κ::R = 0.3
end

