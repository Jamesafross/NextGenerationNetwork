function networksetup(c,constant_delay;digits=10,nSC=2,nFC=1,N=140)

    SC,dist,FC,missingROIs= get_paul_data_all(;nSC=nSC,nFC=nFC,type="control",ROI=N)
    N = size(SC,1)
    lags = dist./c
    lags = round.(lags,digits=digits)
    lags[lags .> 0.0] = lags[lags .> 0.0] .+ constant_delay
    return SC,dist,lags,N,FC,missingROIs

    
end


function get_paul_data_all(;nSC=1,nFC=1,type="control",ROI=140)

    #FUNCTDIR="$DATADIR/Functional/$ROIDIR"
    
    
    SC = get_stuct_data(;n=nSC,ROI=ROI)
    dist = get_dist_data(;n=nSC,ROI=ROI)
    FC,missingROIs = get_functonal_data(;n=nFC,type=type,ROI=ROI)
    
    return SC,dist,FC,missingROIs
end
    
    
function get_stuct_data(;n=1,ROI=140,mean=false)

    if ROI in [18,64,140,246,503,673] 
        HOMEDIR = homedir()
        ROIDIR = "ROI"*"$ROI"
        DATADIR = "$HOMEDIR/PaulJuliaData_ALL"
        STRUCTDIR="$DATADIR/Structural/$ROIDIR"
        StructMats = readdir(STRUCTDIR)
        numMats = size(StructMats,1)
        if mean == false
            if n > numMats
                @error "There are only $numMats structural matrices in this directory. Please choose n < 5"
            else
                print()
                SC = load("$STRUCTDIR/$(StructMats[n])","$(split(StructMats[n],".")[1])")

                SC = log.(SC)
                SC[SC.==-Inf] .= 0

                SC .= SC./maximum(SC)

                return SC
            end
        else
            SC=zeros(ROI,ROI)
            for i = 1:3
                SC[:,:] += load("$STRUCTDIR/$(StructMats[i])","$(split(StructMats[i],".")[1])")/(3)
                SC = log.(SC)
                SC[SC.==-Inf] .= 0

                SC .= SC./maximum(SC)

                return SC
            end
        end
    else
        @error "no such ROI size, choose an ROI in [18,64,140,246,503,673]"
    end

end

function get_dist_data(;n=1,ROI=140)
    if ROI in [18,64,140,246,503,673] 
        HOMEDIR = homedir()
        ROIDIR = "ROI"*"$ROI"
        DATADIR = "$HOMEDIR/PaulJuliaData_ALL"
        DISTDIR="$DATADIR/Distance/$ROIDIR"
        DistMats = readdir(DISTDIR)
        numMats = size(DistMats,1)
        if n > numMats
            @error "There are only $numMats distance matrices in this directory. Please choose n < 5"
        else

            dist = load("$DISTDIR/$(DistMats[n])","$(split(DistMats[n],".")[1])")

            return dist
        end
    else
        @error "no such ROI size, choose an ROI in [18,64,140,246,503,673]"
    end
end

function get_functonal_data(;n=1,type="control",ROI=140)
    if ROI in [18,64,140,246,503,673] 
        HOMEDIR = homedir()
        ROIDIR = "ROI"*"$ROI"
        DATADIR = "$HOMEDIR/PaulJuliaData_ALL/Functional"
        if lowercase(type) == "control"
            FUNCTIONALDIR="$DATADIR/$ROIDIR/Control_Groups"
        elseif lowercase(type) == "stimulated"
            FUNCTIONALDIR="$DATADIR/$ROIDIR/Stimulated_Groups"
        end

        FunctionalDIR = readdir(FUNCTIONALDIR)[n]
        
        
        #println(FUNCTIONALDIR*"/"*FunctionalDIR)
        
        FunctionalMats = readdir(FUNCTIONALDIR*"/"*FunctionalDIR)
        FunctionalMats = FunctionalMats[FunctionalMats .!== "MISSING_ROIs.jld"]
        missingROIs = load(FUNCTIONALDIR*"/"*FunctionalDIR*"/"*"MISSING_ROIs.jld","MISSING_ROIs")

        sizeFC = ROI - size(missingROIs,1)
        FC = zeros(sizeFC,sizeFC)

        numMats = size(FunctionalMats,1)
        
        for i in FunctionalMats
            FC += load("$FUNCTIONALDIR/$FunctionalDIR/$i","$(split(i,".")[1])")
        
        end


        return FC./numMats,missingROIs
    else
        @error "no such ROI size, choose an ROI in [18,64,140,246,503,673]"
    end
end

function get_mean_all_functional_data(;ROI=140,type="control")
    HOMEDIR = homedir()
    ROIDIR = "ROI"*"$ROI"
    DATADIR = "$HOMEDIR/PaulJuliaData_ALL/Functional"
    if lowercase(type) == "control"
        FUNCTIONALDIR="$DATADIR/$ROIDIR/Control_Groups"
    elseif lowercase(type) == "stimulated"
        FUNCTIONALDIR="$DATADIR/$ROIDIR/Stimulated_Groups"
    end
    A = []
    missingROI = []
    ngroups = length(readdir(FUNCTIONALDIR))
    
    for i = 1:ngroups
 
        FC_GROUP,missingROIs = get_functonal_data(;n=i,type=type,ROI=ROI)
        missingROI = cat(missingROI,missingROIs,dims=2)

        A = cat(A,[FC_GROUP],dims=1)
    end
    
    missingROI = unique(missingROI)
    newROI = ROI-size(missingROI,1)
    meanFC = zeros(newROI,newROI)


    
    if size(missingROI,1) > 0
        keepElements = Int.(ones(N))
        for i in missingROIs
            keepElements .= collect(1:N) != i
        end
        for i = 1:ngroups
            meanFC += A[i][keepElements,keepElements]/ngroups
        end
    else
        for i = 1:ngroups
            meanFC += A[i]/ngroups
        end
        
    end

    

    
    return meanFC.-diagm(ones(ROI)) , missingROI


end


function get_HCP_mat()
    file = matopen("/Users/james/abeysuriya_wc_isp/SC.mat")
    SC = read(file, "SC") # note that this does NOT introduce a variable ``varname`` into scope
    close(file)

    N = size(SC,1)
    return SC,N
end