function normalise(W,N)
    for ii = 1:N
     W[W .< 0.0] .= 0.0
         if sum(W[:,ii]) != 0.0
         @views W[:,ii] = W[:,ii]./sum(W[:,ii])
         end
 
     end
      @inbounds for k=1:N #Maintaining symmetry in the weights between regions
         @inbounds for l = k:N
                  W[k,l] = (W[k,l]+W[l,k])/2
                  W[l,k] = W[k,l]
         end
     end
     return W
 end

 function fit_r(modelFC,realFC)
    N = size(modelFC,1)
    modelFCLT = zeros(Int((N^2-N)/2))
    realFCLT = zeros(Int((N^2-N)/2))
    c = 1
    for i = 2:N
            for j = 1:i-1
                    modelFCLT[c] = modelFC[i,j]
                    realFCLT[c] = realFC[i,j]
                    c+=1
            end
    end
    return cor(modelFCLT,realFCLT)
end

function make_uhist(tgrid,u)
    sizeT = size(tgrid,1)
    t = LinRange(tgrid[1],tgrid[end],size(u,2))
    interp = []
    for i = 1:size(u,1)
        if i == 1
            interp = [CubicSplineInterpolation(t,u[i,:])]
        else
            interp = cat(interp,[CubicSplineInterpolation(t,u[i,:])],dims=1) 
        end
    end
    return interp
end

function stim(t,i,stimParams,current_window)
    @unpack stim,stimWindow,stimNodes,stimStr,Tstim = stimParams

    if i ∈ stimNodes && (Tstim[1] <t < Tstim[2]) && (stim == true) && current_window == stimWindow
     
        return stimStr

    else 
        return 0.
    end
 end

 function getFCwindows(SIG,N)
    FC = zeros(N,N)
        for ii = 1:N
            for jj = ii+1:N
                FC[ii,jj] = cor(SIG[ii,:],SIG[jj,:])
                FC[jj,ii] =  FC[ii,jj]
            end
        end
    return FC
end


function get_FC(BOLD,tstart,step_i,step_j)
    N = size(BOLD,1)
    counterT = 1
    for i = tstart:step_i:size(BOLD,2)
        j = i + step_j
        if j < size(BOLD,2)
            counterT +=1
        else
            break
        end
    end
    MODEL_FC = zeros(N,N,counterT) 
    counter = 1
    for i = tstart:step_i:size(BOLD,2)
        j = i + step_j
        if j < size(BOLD,2)
            MODEL_FC[:,:,counter] = getFCwindows(BOLD[:,i:j],N)
            counter += 1
            #println(counter)
          
        else
            break
        end
    end
    return MODEL_FC[:,:,1:end-1]
end




