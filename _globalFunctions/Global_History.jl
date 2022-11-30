function h1(hparams,t;idxs = nothing)
    #history function used on first window
        if t < 0
        return hparams[idxs]
    end
end

function h2(hparams,t;idxs = nothing)
    #history function used on windows > 1
        if t < 0
            return hparams[idxs](t)
        end
end

function make_hist_mat2_threads!(h,W::Matrix{Float64},u::Vector{Float64},hparams,N::Real,lags::Array{Float64},t::Float64,WHistMat::Array{Float64})
    @inbounds Threads.@threads for i = 1:N
            for j = 1:N
                if lags[j,i] > 0
                  
                    WHistMat[j,i] = W[j,i]*h(hparams,t-lags[j,i]; idxs=i)
                else
                    WHistMat[j,i] = W[j,i]*u[i]
                end
            end
        end
end

function make_hist_mat2_unthreads!(h,W::Matrix{Float64},u::Vector{Float64},hparams,N::Real,lags::Array{Float64},t::Float64,WHistMat::Array{Float64})
    for i = 1:N
            for j = 1:N
                if lags[j,i] > 0
                    WHistMat[j,i] = W[j,i]*h(hparams,t-lags[j,i]; idxs=j)
                else
                    WHistMat[j,i] = W[j,i]*u[j]
                end
            end
        end
end


        
function make_hist_mat!(h,W::Matrix{Float64},u::Vector{Float64},hparams,N::Real,lags::Array{Float64},t::Float64,WHistMat::Array{Float64},non_zero_weights)
        @inbounds Threads.@threads for i in non_zero_weights
                    WHistMat[i] = W[i]*h(hparams,t-lags[i]; idxs=i[2])
            end
end

function make_d!(W::Matrix{Float64},HistMat::Matrix{Float64},d::Vector{Float64})
    d .= sum(W.*HistMat,dims=2)
end