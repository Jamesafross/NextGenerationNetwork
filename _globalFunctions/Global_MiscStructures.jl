mutable struct TimerStruct
    meanIntegrationTime
    meanBalloonTime
    meanAllTime
 end

 mutable struct times_save
    time
    network_size
    network_density
 end

mutable struct init
   u0
end

mutable struct adaptParameters
   LearningRate::Float64
   windowStart::Int
   tP::Float64
   HIST::Array{Float64}
end



 