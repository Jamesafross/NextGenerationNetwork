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
        s=0.

        
        rE=u[1]
        rI=u[2]
        vE=u[3]
        vI=u[4]
        gEE=u[5]
        pEE=u[6]
        gIE=u[7]
        pIE=u[8]
        gEI=u[9]
        pEI=u[10]
        gII=u[11]
        pII=u[12]

        if  (10000 < t < 10100)
            s = AMP*sin(fr*t)
        else
            s = 0.0;
        end
        
        if t < 9000
            τisp=0.0
        else
            τisp=0.00
        end 

        #rE
        du[1] =(1. /τE)*(-gEE*rE -gEI*rE + 2. * rE * vE + (ΔE / (τE*pi)))
        #rI
        du[2] =(1. /τI)*(-gIE*rI - gII*rI + 2. * rI * vI + (ΔI / (τI*pi)))
        #vE
        du[3] =(1. /τE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vE) - (τE^2)*(pi^2) * (rE^2.) +  vE^2. + η_0E +s)
        #vI
        du[4] =(1. /τI)*(gIE*(VsynIE - vI) + gII*(VsynII - vI) - (τI^2)*(pi^2) * (rI^2.) + vI^2. + η_0I)
        #gEE
        du[5] = αEE * (-gEE + pEE)
        #pEE 
        du[6] = αEE * (-pEE + κSEE*rE)
        #gIE
        du[7] = αIE * (-gIE + pIE)
        #pIE
        du[8] = αIE * (-pIE + κSIE * rE)
        #gEI
        du[9] = αEI * (-gEI + pEI)
        #pEI 
        du[10] = αEI * (-pEI + u[13]* rI)
        #gII
        du[11] = αII * (-gII + pII)
        #pII
        du[12] = αII * (-pII + κSII * rI)

        du[13] = τisp*((rI-0.7496555947801044)*(rE - rI))


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
    κSEE::R = 3.0
    κSIE::R = 2.
    κSEI::R = 1.
    κSII::R = 1.5
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