

using Distributed
Distributed.addprocs(active_cpu_cores-1); 
@everywhere begin 
    using ParallelDataTransfer
    using GeoArrays
    using Random
    using Statistics
    
    @noinline function generate_List(offset::Int64)::Vector{Int64}    
        """Generate a list of possible movements based on a user specified offset"""
        
        return push!(vcat([[n, -n] for n in 1:offset]...),0)  
    end
    
    @noinline function convert_to_offsets(x::String)::Vector{Int64}  
        """Convert a user specified number of neighboring pixels to an offset (x,y)"""
        
        return Dict("0" => [0,0], "8" => [1,1], "14" => [1,2], "20" => [1,3], "24" => [2,2], "34" => [2,3])[x]
    end
    
    @noinline function selectUP_minNx(idxy::Tuple{Int64, Int64, Int64}, 
                                    L1::Vector{Int64}, 
                                    L2::Vector{Int64}, 
                                    minNx::Int64,                                    
                                    targetImgBand::GeoArray{Float32, Array{Float32, 3}})
        """From a list of coordinates, select those who satisfy a user specified minNx criterion"""
        
        idx = idxy[1] 
        idy = idxy[2] 
        Nx = []
        @simd for k1 in L1
            m1 = (idx+k1)
            if (m1>=1) && (m1<=size(targetImgBand)[1])      
                @simd for k2 in L2
                    if (k1 != 0) || (k2 != 0)           
                        m2 = (idy+k2)
                        if (m2>=1) && (m2<=size(targetImgBand)[2]) 
                            @inbounds Zx = targetImgBand[m1,m2][1]
                            @inbounds if (Zx>=0) && (Zx<1e+20)
                                @inbounds push!(Nx, Zx)                            
                            end
                        end
                    end
                end
            end
        end
        if length(Nx)>=minNx
            return idx,idy
        end
    end
    
    @noinline function getZx(L1::Vector{Int64},
                            L2::Vector{Int64},
                            idx::Int64,
                            idy::Int64,
                            ImgBand::GeoArray{Float32, Array{Float32, 3}})
        """generate valid lagv, lagh, Zx values"""
        l1=Int64[]
        l2=Int64[]
        l3=Float64[]
        @simd for k1 in L1 
            m1 = idx+k1            
            if (m1>=1) && (m1<=size(ImgBand)[1])               
                @simd for k2 in L2    
                    if (k1!=0) || (k2!=0)                        
                        m2 = idy+k2                        
                        if (m2>=1) && (m2<=size(ImgBand)[2])   
                            @inbounds Zx = ImgBand[m1, m2][1]
                            @inbounds if (Zx>=0) && (Zx<1e+20)        
                                push!(l1,k1)
                                push!(l2,k2)
                                push!(l3,Zx)
                            end
                        end
                    end
                end
            end
        end
        return l1,l2,l3
    end 
    
    @noinline function getZy(L1::Vector{Int64},
                            L2::Vector{Int64},
                            idx::Int64,
                            idy::Int64,
                            ImgBand::GeoArray{Float32, Array{Float32, 3}})
        """generate valid Zy values"""
        
        l1=Float64[]
        @simd for p in 1:length(L1)
            b1 = (idx+L1[p])
            b2 = (idy+L2[p])
            if (b1>0) && (b2>0) && (b1<=size(ImgBand)[1]) && (b2<=size(ImgBand)[2]) 
                @inbounds Zy = ImgBand[b1,b2][1]
                @inbounds if (Zy>=0) && (Zy<1e+20)
                    push!(l1,Zy)  
                end
            end
        end 
        return l1
    end
        
    @noinline function find_subset1(Nx::Vector{Float64},
                                list_boundary_values::Vector{Float32},
                                coordsKP::Vector{Tuple{Int64, Int64, Int64}},
                                valuesKP::Vector{Float32},
                                f::Float64,
                                seed1::Int64)
        """A function that narrows down the search area to known pixels whose range corresponds to Nx pixel values"""
        
                subset_coordsKP_temp = []
                detectDuplicate = Int[]            
                for value in Nx
                    for l in 1:length(list_boundary_values)-1
                        if (value>=list_boundary_values[l]) && (value<list_boundary_values[l+1]) && !(l in detectDuplicate)
                            indices_temp = findall(X -> (X.>=list_boundary_values[l]) && (X.<list_boundary_values[l+1]), valuesKP) 
                            append!(subset_coordsKP_temp, [coordsKP[elt] for elt in indices_temp]) 
                            push!(detectDuplicate,l)
                        end
                    end
                end        
                Random.seed!(seed1)
                subset_list_KP=[subset_coordsKP_temp[tempID] for tempID in randperm(length(subset_coordsKP_temp))[1:trunc(Int, length(subset_coordsKP_temp)*f)]] # get subset based on specified fraction # Random selection without repetition
            
                return subset_list_KP                
            end
    
    
    @noinline function find_best_replicate(w)     
        """The core function for finding the best replicate value of a selected unknown pixel"""
        
        idx = re_selectedUP[w][1] 
        idy = re_selectedUP[w][2] 
        IDx = 999999                                            
        IDy = 999999 
        replicateValue = 999999 
        ned_temp = 999999
        ned_final = 999999  
        lagv1,lagh1,Nx1 = getZx(L1,L2,idx,idy,targetImgBand) 
        if (DS == "_UV_")
            subset_list_KP = find_subset1(Nx1, list_boundary_values, coordsKP, valuesKP, f, seed1) 
            for item in subset_list_KP                             
                Ny = Float64[]
                id1 = item[1]
                id2 = item[2]
                Ny = getZy(lagv1,lagh1,id1,id2,trainingImgBand)                    
                if length(Ny)==length(Nx1) 
                    @inbounds ned_temp = sqrt(sum([(Nx1[o]-Ny[o])^2 for o in 1:length(Ny)]))/dmax
                    if (ned_temp<t1)
                        IDx = id1
                        IDy = id2
                        ned_final = ned_temp                        
                        break
                    else
                        if (ned_temp<ned_final)
                            IDx = id1
                            IDy = id2
                            ned_final = ned_temp          
                        end
                    end
                end
            end
            if (ned_final!=999999)                                                   
                @inbounds replicateValue = trainingImgBand[IDx, IDy][1]                     
            else                                                    
                Random.seed!()
                replicateValue = [valuesKP[tempID] for tempID in randperm(length(valuesKP))[1]][1] 
            end
        elseif (DS == "_BV_")
            lagv2,lagh2,Nx2 = getZx(L1,L2,idx,idy,auxImgBand)
            subset_list_KP = find_subset1(Nx1, list_boundary_values, coordsKP1, valuesKP1, f, seed1) 
            for item in subset_list_KP                                   
                id1 = item[1]
                id2 = item[2]
                Ny1 = getZy(lagv1,lagh1,id1,id2,targetImgBand)                    
                Ny2 = getZy(lagv2,lagh2,id1,id2,auxImgBand)  

                if (length(Ny1)==length(Nx1)) && (length(Ny2)==length(Nx2))     
                    @inbounds temp1 = w1*sqrt(sum([(Nx1[o]-Ny1[o])^2 for o in 1:length(Ny1)]))/dmax1
                    @inbounds temp2 = w2*sqrt(sum([(Nx2[o]-Ny2[o])^2 for o in 1:length(Ny2)]))/dmax2
                    ned_temp = temp1+temp2 
                    if (ned_temp<t1)
                        IDx = id1
                        IDy = id2
                        ned_final = ned_temp;                        
                        break
                    else
                        if (ned_temp<ned_final)
                            IDx = id1
                            IDy = id2
                            ned_final = ned_temp  
                        end
                    end
                end
            end
            if (ned_final!=999999)                                                     
                @inbounds replicateValue = targetImgBand[IDx, IDy][1]
            else                                                    
                Random.seed!()
                replicateValue = [valuesKP1[tempID] for tempID in randperm(length(valuesKP1))[1]][1] 
            end
        end 
        return replicateValue, idx, idy 
    end    
end;

function analyze(N1, N2, T1, F, minNxs, W1, nb_repitions, bands, active_cpu_cores, 
    imgType, write_img, DS, pixels_to_be_filled, missingDataImg, folders, path, param)

    dir = joinpath(path, imgType)
    for folder in folders 
        input_dir1 = joinpath(dir, "Preprocessed", "Step5", folder, pixels_to_be_filled) 
        df1 = CSV.read(joinpath(input_dir1,"pathsData.csv"), DataFrame)

        output_dir = joinpath(dir, "Analysis", folder, pixels_to_be_filled)
        output_dir_img = joinpath(dir, "NewImg", folder, pixels_to_be_filled)
        output_dir_bias = joinpath(dir, "BiasImg",folder, pixels_to_be_filled)
        mkpath(output_dir)  
        mkpath(output_dir_img) 
        mkpath(output_dir_bias) 

        combinations = []
        totalLength = nb_repitions*length(df1.Ref_Img)*length(N1)*length(bands)*length(T1)*length(F)*length(N2)*length(minNxs)
        if (DS == "_BV_")
            totalLength = totalLength*length(W1)
        end
        p = Progress(totalLength, dt=1, barlen=50, color=:black)
        statusFile = [file for file in readdir(output_dir; join=true) if occursin("errorResults"*DS*imgType*pixels_to_be_filled*"_rep"*string(nb_repitions)*"_"*param*".csv",file)]
        for i in 1:length(df1.Ref_Img)
            for band in bands
                for seed1 in 1:nb_repitions
                    for n1 in N1
                        for n2 in N2
                            for t1 in T1
                                for f in F
                                    for minNx in minNxs
                                        if (DS == "_BV_")
                                            for w1 in W1
                                                push!(combinations, [i,band,seed1,n1,n2,t1,f,minNx,w1,1-w1])
                                                next!(p)
                                            end
                                        elseif (DS == "_UV_")
                                            push!(combinations, [i,band,seed1,n1,n2,t1,f,minNx])
                                            next!(p)
                                        end
                                    end 
                                end
                            end
                        end
                    end
                end
            end
        end
        Factor1 = String[]
        Factor2 = Int64[]
        Factor3 = Float64[]
        Factor4 = Float64[]
        Factor5 = Int64[]
        Factor6 = Float64[]
        Factor7 = Float64[]

        ITERATION = Int64[]
        DATE_REF_IMG = String[]
        DATE_TARGET_IMG = String[]
        LOCATION = String[]
        MASKP = Float64[]
        TOTAL_ELAPSED_TIME_SEC_PER_BAND = Float64[]
        ELAPSED_TIME_SEC_PER_BAND = Float64[]
        ELAPSED_TIME_MIN_PER_BAND = Float64[]
        ELAPSED_TIME_HOURS_PER_BAND = Float64[]
        ELAPSED_TIME_DAYS_PER_BAND = Float64[]
        BAND = Int64[]
        ACCURACY_MSE = Float64[]
        ACCURACY_MSLE = Float64[]
        ACCURACY_R2 = Float64[]
        N_missingPixels = Float64[]
        N_filledPixels = Float64[]
        if length(statusFile)!=0
            df2 = CSV.read(statusFile, DataFrame)
            for index in 1:length(df2.Location)           
                push!(ITERATION,df2.Iteration[index])
                push!(DATE_REF_IMG,string(df2.Date_Ref[index]))
                push!(DATE_TARGET_IMG,string(df2.Date_Target[index]))
                push!(LOCATION,df2.Location[index])
                push!(MASKP,df2.MaskedP[index])
                push!(TOTAL_ELAPSED_TIME_SEC_PER_BAND,df2.Total_Time_Sec[index])
                push!(ELAPSED_TIME_SEC_PER_BAND,df2.Time_Sec[index])
                push!(ELAPSED_TIME_MIN_PER_BAND,df2.Time_Min[index])
                push!(ELAPSED_TIME_HOURS_PER_BAND,df2.Time_Hours[index])
                push!(ELAPSED_TIME_DAYS_PER_BAND,df2.Time_Days[index])
                push!(BAND,df2.Band[index])
                push!(ACCURACY_MSE,df2.MSE[index])
                push!(ACCURACY_MSLE,df2.MSLE[index])
                push!(ACCURACY_R2,df2.R2[index])
                push!(N_missingPixels,df2.UP[index])
                push!(N_filledPixels,df2.KP[index])

                push!(Factor1,string(df2.n1[index]))
                push!(Factor2,df2.n2[index])
                push!(Factor3,df2.t1[index])
                push!(Factor4,df2.f[index])
                push!(Factor5,df2.minNx[index]) 
                if (DS == "_BV_")
                    push!(Factor6,df2.w1[index])
                    push!(Factor7,df2.w2[index]) 
                end  
            end
            combinations = combinations[(length(df2.Location)+1):end]
        end
        println(" >> N Combinations:"*string(size(combinations)[1]))
        
        p = Progress(length(combinations), dt=1, barlen=50, color=:black)
        for combination in combinations
            elapsed_sec = @elapsed begin

                i = combination[1]
                band = combination[2]
                seed1 = combination[3]
                n1 = combination[4]
                n2 = combination[5]
                t1 = combination[6]
                f = combination[7]
                minNx = combination[8] 
                w1 = 1.0
                w2 = 1.0
                if (DS == "_BV_")
                    w1 = combination[9]
                    w2 = combination[10]
                end

                state_r2 = "positive"
                statusFile = [file for file in readdir(output_dir; join=true) if occursin("errorResults"*DS*imgType*pixels_to_be_filled*"_rep"*string(nb_repitions)*"_"*param*".csv",file)]
                if length(statusFile) !==0
                    df2 = CSV.read(statusFile, DataFrame)
                    row = length(df2.Location)                 
                    if row !== 0
                        if ((df2[row,"R2"])<0) && (reverse(reverse(df1.Target_Img[i])[5:12])==string(df2[row,"Date_Target"])) && (df1.missingPixelsTarget[i]==df2[row,"MaskedP"]) && (df2[row,"minNx"]>0) && (n1==string(df2[row,"n1"])) && (band==df2[row,"Band"]) && (t1==df2[row,"t1"])
                            state_r2 = "negative"                     
                        end                
                    end
                end 

                if state_r2 == "positive"
                    L1 = generate_List(convert_to_offsets(n1)[1])   
                    L2 = generate_List(convert_to_offsets(n1)[2])    

                    @inbounds trainingImg = df1.Target_Img[i]
                    @inbounds auxImg = df1.Img_Before1[i]             
                    @inbounds refImgBand = GeoArrays.read(df1.Ref_Img[i], band=band)
                    @inbounds targetImgBand = GeoArrays.read(df1.Target_Img[i], band=band)
                    @inbounds targetImgBand_bias = GeoArrays.read(df1.Target_Img[i], band=band)         
                    if (DS == "_UV_")
                        case = "targetImg"
                        cloudValue = df1.missingPixelsTarget[i]
                        if (cloudValue>missingDataImg)          
                            case = "auxImg"
                            @inbounds if (df1.missingPixelsBefore1[i]!=100) || (df1.missingPixelsAfter1[i]!=100)
                                @inbounds if (cloudValue>df1.missingPixelsBefore1[i]) || (cloudValue>df1.missingPixelsAfter1[i]) 
                                    @inbounds if df1.missingPixelsBefore1[i]>df1.missingPixelsAfter1[i]  
                                        @inbounds trainingImg = df1.Img_After1[i]
                                    else
                                        @inbounds trainingImg = df1.Img_Before1[i]
                                    end
                                end
                            else
                                @inbounds if (cloudValue>df1.missingPixelsBefore2[i]) || (cloudValue>df1.missingPixelsAfter2[i])                                                                 
                                    @inbounds if df1.missingPixelsBefore2[i]>df1.missingPixelsAfter2[i]  
                                        @inbounds trainingImg = df1.Img_After2[i]
                                    else
                                        @inbounds trainingImg = df1.Img_Before2[i]
                                    end
                                end
                            end 
                        end  

                        @inbounds trainingImgBand = GeoArrays.read(trainingImg, band=band)  
                        coordsUP_temp = Tuple.(findall(x->(x==-99), refImgBand)) 
                        for item in coordsUP_temp
                            @inbounds targetImgBand[item[1],item[2],1] = 1e+21
                            @inbounds trainingImgBand[item[1],item[2],1] = 1e+21                
                        end
                        coordsKP = Tuple.(findall(x->(x>=0)&&(x<1e+20), trainingImgBand))
                        valuesKP = [trainingImgBand[coord[1],coord[2],coord[3]] for coord in coordsKP]
                        minKP = minimum(valuesKP) 
                        maxKP = maximum(valuesKP)
                        dmax = maxKP - minKP
                        list_boundary_values = [b_value for b_value in minKP:(dmax/n2):(maxKP+(dmax/n2))] 
                        actual_coordsUP = Tuple.(findall(x->(x==-99), targetImgBand))
                        NumberUP = 1
                        while (NumberUP != 0)
                            coordsUP = Tuple.(findall(x->(x==-99), targetImgBand)) 
                            selectedUP = [selectUP_minNx(coordsUP[v], L1, L2, minNx, targetImgBand) 
                                for v in 1:length(coordsUP) if selectUP_minNx(coordsUP[v], L1, L2, minNx, targetImgBand) !== nothing]
                            re_selectedUP = selectedUP
                            if(length(selectedUP) != 0)                                     
                                Random.seed!(seed1)
                                re_selectedUP = [selectedUP[tempID] for tempID in randperm(length(selectedUP))]                          

                            else                                                       
                                if length(coordsUP)<5
                                    re_selectedUP = coordsUP
                                else
                                    Random.seed!(seed1)   
                                    re_selectedUP = [coordsUP[tempID] for tempID in randperm(length(coordsUP))[1:5]] 
                                end
                            end   
                            sendto(workers(), re_selectedUP=re_selectedUP,
                                                L1=L1::Vector{Int64},
                                                L2=L2::Vector{Int64},
                                                targetImgBand=targetImgBand::GeoArray{Float32, Array{Float32, 3}},
                                                trainingImgBand=trainingImgBand::GeoArray{Float32, Array{Float32, 3}},
                                                f=f::Float64,
                                                coordsKP=coordsKP::Vector{Tuple{Int64, Int64, Int64}},
                                                valuesKP=valuesKP::Vector{Float32},
                                                list_boundary_values=list_boundary_values::Vector{Float32},
                                                t1=t1::Float64,
                                                dmax=dmax::Float32,
                                                seed1=seed1::Int64,
                                                DS=DS::String)
                            results = Distributed.pmap(find_best_replicate, 1:length(re_selectedUP))
                            for index in 1:length(results)
                                @inbounds targetImgBand[re_selectedUP[index][1],re_selectedUP[index][2],1] = results[index][1]
                            end
                            if case == "targetImg"
                                trainingImgBand = copy(targetImgBand)
                            end                
                            NumberUP = length(Tuple.(findall(x->(x==-99), targetImgBand)))  
                        end
                    elseif (DS == "_BV_")
                        if (df1.missingPixelsBefore1[i]>df1.missingPixelsAfter1[i]) && (df1.missingPixelsAfter1[i]<100)  
                            @inbounds auxImg = df1.Img_After1[i]
                        elseif (df1.missingPixelsBefore1[i]<df1.missingPixelsAfter1[i]) && (df1.missingPixelsBefore1[i]<100)
                            @inbounds auxImg = df1.Img_Before1[i]
                        elseif (df1.missingPixelsBefore2[i]>df1.missingPixelsAfter2[i]) && (df1.missingPixelsAfter2[i]<100)
                            @inbounds auxImg = df1.Img_After2[i]
                        elseif (df1.missingPixelsBefore2[i]<df1.missingPixelsAfter2[i]) && (df1.missingPixelsBefore2[i]<100)
                            @inbounds auxImg = df1.Img_Before2[i]
                        else
                            println("No valid auxImg")
                            break
                        end
                        @inbounds auxImgBand = GeoArrays.read(auxImg, band=band) 
                        coordsUP_temp = Tuple.(findall(x->(x==-99), refImgBand)) 
                        for item in coordsUP_temp
                            @inbounds targetImgBand[item[1],item[2],1] = 1e+21              
                        end
                        coordsKP1 = Tuple.(findall(x->(x>=0)&&(x<1e+20), targetImgBand)) 
                        valuesKP1 = [targetImgBand[coord[1],coord[2],coord[3]] for coord in coordsKP1]
                        minKP1 = minimum(valuesKP1) 
                        maxKP1 = maximum(valuesKP1)
                        dmax1 = maxKP1 - minKP1
                        coordsKP2 = Tuple.(findall(x->(x>=0)&&(x<1e+20), auxImgBand))
                        valuesKP2 = [auxImgBand[coord[1],coord[2],coord[3]] for coord in coordsKP2]
                        minKP2 = minimum(valuesKP2) 
                        maxKP2 = maximum(valuesKP2)
                        dmax2 = maxKP2 - minKP2
                        list_boundary_values = [b_value for b_value in minKP1:(dmax1/n2):(maxKP1+(dmax1/n2))] 
                        actual_coordsUP = Tuple.(findall(x->(x==-99), targetImgBand))
                        NumberUP = length(actual_coordsUP)
                        while (NumberUP != 0)
                            oordsUP = Tuple.(findall(x->(x==-99), targetImgBand)) 
                            selectedUP = [selectUP_minNx(coordsUP[v], L1, L2, minNx, targetImgBand) 
                                for v in 1:length(coordsUP) if selectUP_minNx(coordsUP[v], L1, L2, minNx, targetImgBand) !== nothing]
                            re_selectedUP = selectedUP
                            if(length(selectedUP) != 0)                                     
                                Random.seed!(seed1)                        
                                re_selectedUP = [selectedUP[tempID] for tempID in randperm(length(selectedUP))] 
                            else                                                       
                                if length(coordsUP)<5
                                    re_selectedUP = coordsUP
                                else
                                    Random.seed!(seed1)   
                                    re_selectedUP = [coordsUP[tempID] for tempID in randperm(length(coordsUP))[1:5]] 
                                end
                            end                    
                            sendto(workers(), re_selectedUP=re_selectedUP,
                                                L1=L1::Vector{Int64},
                                                L2=L2::Vector{Int64},
                                                targetImgBand=targetImgBand::GeoArray{Float32, Array{Float32, 3}},
                                                auxImgBand=auxImgBand::GeoArray{Float32, Array{Float32, 3}},
                                                f=f::Float64,
                                                coordsKP1=coordsKP1::Vector{Tuple{Int64, Int64, Int64}},
                                                valuesKP1=valuesKP1::Vector{Float32},
                                                list_boundary_values=list_boundary_values::Vector{Float32},
                                                t1=t1::Float64,
                                                dmax1=dmax1::Float32,
                                                dmax2=dmax2::Float32,
                                                w1=w1::Float64, 
                                                w2=w2::Float64,
                                                seed1=seed1::Int64,
                                                DS=DS::String)
                            results = Distributed.pmap(find_best_replicate, 1:length(re_selectedUP))
                            for index in 1:length(results)
                                @inbounds targetImgBand[re_selectedUP[index][1],re_selectedUP[index][2],1] = results[index][1]
                            end
                            NumberUP = length(Tuple.(findall(x->(x==-99), targetImgBand)))                     
                        end
                    end
                    acualKPvalues = Float64[]
                    newKPvalues = Float64[] 
                    for coordUP in actual_coordsUP 
                        push!(acualKPvalues, refImgBand[coordUP[1],coordUP[2],coordUP[3]])
                        push!(newKPvalues, targetImgBand[coordUP[1],coordUP[2],coordUP[3]])
                    end 
                    n_mPixels = length(Tuple.(findall(x->(x==-99) && (x==-101), targetImgBand))) 
                    n_nanPixels = length(Tuple.(findall(x->(isnan(x)), targetImgBand)))        
                    percent_mPixels = 100*round(n_mPixels/((size(targetImgBand)[1]*size(targetImgBand)[2])); digits=4) 
                    percent_nanPixels = 100*round(n_nanPixels/(size(targetImgBand)[1]*size(targetImgBand)[2]); digits=4)
                    percent_fPixels = 100 - percent_mPixels - percent_nanPixels
                    push!(BAND, band)
                    push!(ACCURACY_MSE, mse(newKPvalues, acualKPvalues))
                    push!(ACCURACY_MSLE, msle(newKPvalues, acualKPvalues))
                    push!(ACCURACY_R2, r2_score(newKPvalues, acualKPvalues))
                    push!(N_missingPixels,percent_mPixels)
                    push!(N_filledPixels,percent_fPixels)
                else
                    push!(BAND, band)
                    push!(ACCURACY_MSE, -99)
                    push!(ACCURACY_MSLE, -99)
                    push!(ACCURACY_R2, -99)
                    push!(N_missingPixels,-99)
                    push!(N_filledPixels,-99)
                end
            end
            push!(ITERATION,seed1)
            push!(DATE_REF_IMG,reverse(reverse(df1.Ref_Img[i])[12:19]))
            push!(DATE_TARGET_IMG,reverse(reverse(df1.Target_Img[i])[12:19]))    
            push!(LOCATION,df1.Location[i])
            push!(MASKP,df1.missingPixelsTarget[i])
            push!(Factor1,n1)
            push!(Factor2,n2)
            push!(Factor3,t1)
            push!(Factor4,f)
            push!(Factor5,minNx)
            if (DS == "_BV_")
                push!(Factor6,w1)
                push!(Factor7,w2) 
            end
            t_d = floor.(Int,elapsed_sec/60/60/24)
            t_hr = floor.(Int,abs(elapsed_sec-t_d*60*60*24)/60/60)
            t_min = floor.(Int,abs(abs(elapsed_sec-t_d*60*60*24)-t_hr*60*60)/60)
            t_sec = abs(abs(abs(elapsed_sec-t_d*60*60*24)-t_hr*60*60)-t_min*60)  
            push!(TOTAL_ELAPSED_TIME_SEC_PER_BAND, round(elapsed_sec))
            push!(ELAPSED_TIME_SEC_PER_BAND, t_sec)
            push!(ELAPSED_TIME_MIN_PER_BAND, t_min)
            push!(ELAPSED_TIME_HOURS_PER_BAND, t_hr)
            push!(ELAPSED_TIME_DAYS_PER_BAND, t_d)   

            if (DS == "_UV_") & (write_img != "y")
                df = DataFrame(Iteration=ITERATION, 
                            Location=LOCATION, 
                            Date_Ref=DATE_REF_IMG, 
                            Date_Target=DATE_TARGET_IMG, 
                            MaskedP=MASKP, 
                            Band=BAND, 
                            n1=Factor1,
                            n2=Factor2,
                            t1=Factor3,
                            f=Factor4,
                            minNx=Factor5,
                            MSLE=ACCURACY_MSLE, 
                            MSE=ACCURACY_MSE, 
                            R2=ACCURACY_R2, 
                            UP=N_missingPixels,
                            KP=N_filledPixels,
                            Total_Time_Sec=TOTAL_ELAPSED_TIME_SEC_PER_BAND, 
                            Time_Sec=ELAPSED_TIME_SEC_PER_BAND, 
                            Time_Min=ELAPSED_TIME_MIN_PER_BAND, 
                            Time_Hours=ELAPSED_TIME_HOURS_PER_BAND, 
                            Time_Days=ELAPSED_TIME_DAYS_PER_BAND)
                CSV.write(joinpath(output_dir, "errorResults"*DS*imgType*pixels_to_be_filled*"_rep"*string(nb_repitions)*"_"*param*".csv"), df)
            elseif (DS == "_BV_") & (write_img != "y")
                df = DataFrame(Iteration=ITERATION, 
                            Location=LOCATION, 
                            Date_Ref=DATE_REF_IMG, 
                            Date_Target=DATE_TARGET_IMG, 
                            MaskedP=MASKP, 
                            Band=BAND, 
                            n1=Factor1,
                            n2=Factor2,
                            t1=Factor3,
                            f=Factor4,
                            minNx=Factor5,
                            w1=Factor6,
                            w2=Factor7,
                            MSLE=ACCURACY_MSLE, 
                            MSE=ACCURACY_MSE, 
                            R2=ACCURACY_R2, 
                            UP=N_missingPixels,
                            KP=N_filledPixels,
                            Total_Time_Sec=TOTAL_ELAPSED_TIME_SEC_PER_BAND, 
                            Time_Sec=ELAPSED_TIME_SEC_PER_BAND, 
                            Time_Min=ELAPSED_TIME_MIN_PER_BAND, 
                            Time_Hours=ELAPSED_TIME_HOURS_PER_BAND, 
                            Time_Days=ELAPSED_TIME_DAYS_PER_BAND)
                CSV.write(joinpath(output_dir, "errorResults"*DS*imgType*pixels_to_be_filled*"_rep"*string(nb_repitions)*"_"*param*".csv"), df)
            end
            if write_img == "y"  
                for i1 in 1:size(targetImgBand_bias)[1]
                    for i2 in 1:size(targetImgBand_bias)[2]
                        @inbounds targetImgBand_bias[i1,i2,1] = 100*((targetImgBand_bias[i1,i2,1]-targetImgBand[i1,i2,1])/targetImgBand[i1,i2,1])
                    end
                end 
                GeoArrays.write!(joinpath(output_dir_img,"B"*string(band)*DS*reverse(reverse(df1.Target_Img[i])[begin:19])), targetImgBand)
                GeoArrays.write!(joinpath(output_dir_bias,"Bias_B"*string(band)*DS*reverse(reverse(df1.Target_Img[i])[begin:19])), targetImgBand_bias)
            end
            next!(p)    
        end
    end
end;


