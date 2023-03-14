

function differenceP(value1,value2)
    return 100*((value1-value2)/value1)
end

function postprocess(nb_rep1, nb_rep2, folders, params)
    for param in params
        for folder in folders
            files1 = [] 
            for imgType in ["L7", "S2"]
                for pixels_to_be_filled in ["_nan","_cloud"]
                    if isdir(joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled)) == true 
                        for file in readdir(joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled))
                            for DS in ["_UV_", "_BV_"]
                                if occursin("errorResults"*DS*imgType*pixels_to_be_filled*"_rep"*string(nb_rep1)*"_"*param*".csv",file)
                                    push!(files1, joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled,file))
                                end                        
                            end
                        end
                    end
                end
            end
            
            selectedParam = []
            for file in files1
                df1 = CSV.read(file, DataFrame) 
                for p in unique(df1[(df1.Location.==folder), :].MaskedP)
                    for band in unique(df1[(df1.Location.==folder), :].Band)
                        df1_filtered = df1[(df1.Location.==folder) .& (df1.MaskedP.==p) .& (df1.Band.==band), :]
                        df1_sorted = sort(df1_filtered, (:R2))                    
                        if length(df1_sorted[!,param])-1 == 0
                            push!(selectedParam, df1_sorted[!,param][length(df1_sorted[!,param])])
                        else
                            if (df1_sorted.Total_Time_Sec[length(df1_sorted.Total_Time_Sec)-1] < df1_sorted.Total_Time_Sec[length(df1_sorted.Total_Time_Sec)]) & (abs(df1_sorted.R2[length(df1_sorted.R2)]-df1_sorted.R2[length(df1_sorted.R2)-1])<=0.015)
                                push!(selectedParam, df1_sorted[!,param][length(df1_sorted[!,param])-1])
                            else
                                push!(selectedParam, df1_sorted[!,param][length(df1_sorted[!,param])])
                            end
                        end 
                    end
                end 
            end     
            df_temp = DataFrame(unique_values=unique(selectedParam), frequency=[sum(selectedParam .== item) for item in unique(selectedParam)])
            df_temp_sorted = sort(df_temp, [:frequency])  
            println(folder,"_",param,": ",df_temp_sorted.unique_values[length(df_temp_sorted.unique_values)]) 
        end
    end



    for param in ["n2", "minNx"]
        for folder in folders
            files1 = [] 
            output_dirs = []
            for imgType in ["L7", "S2"]
                for pixels_to_be_filled in ["_nan","_cloud"]
                    if isdir(joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled)) == true 
                        for file in readdir(joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled))
                            for DS in ["_UV_", "_BV_"]
                                if occursin("errorResults"*DS*imgType*pixels_to_be_filled*"_rep"*string(nb_rep1)*"_"*param*".csv",file)
                                    push!(files1, joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled,file))
                                    push!(output_dirs, joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled))
                                end                        
                            end
                        end
                    end
                end
            end
            
            for i in 1:length(files1)
                GAP = []
                BAND = []
                v1 = []
                v2 = []
                v3 = []
                v4 = []
                v5 = []
                t1 = []
                t2 = []
                t3 = []
                t4 = []
                file = files1[i]
                output_dir = output_dirs[i]
                df1 = CSV.read(file, DataFrame) 
                for p in unique(df1[(df1.Location.==folder), :].MaskedP)
                    for band in unique(df1[(df1.Location.==folder), :].Band)
                        df1_filtered = df1[(df1.Location.==folder) .& (df1.MaskedP.==p) .& (df1.Band.==band), :]
                        push!(GAP, df1_filtered.MaskedP[1])
                        push!(BAND, df1_filtered.Band[1])                                
                        push!(v1, df1_filtered.R2[1])
                        push!(v2, df1_filtered.R2[2])
                        push!(v3, df1_filtered.R2[3])
                        push!(v4, df1_filtered.R2[4]) 
                        push!(t1, differenceP(df1_filtered.Total_Time_Sec[1], df1_filtered.Total_Time_Sec[2]))
                        push!(t2, differenceP(df1_filtered.Total_Time_Sec[1], df1_filtered.Total_Time_Sec[3]))
                        push!(t3, differenceP(df1_filtered.Total_Time_Sec[1], df1_filtered.Total_Time_Sec[4]))
                        if param == "n2"
                            push!(v5, df1_filtered.R2[5])
                            push!(t4, differenceP(df1_filtered.Total_Time_Sec[1], df1_filtered.Total_Time_Sec[5]))
                        end                    
                    end
                end
                
                if param == "n2"
                    df_accuracy_n2 = DataFrame(MaskedP=GAP,Band=BAND,l1=v1,l50=v2,l100=v3,l150=v4,l200=v5)
                    df_time_n2 = DataFrame(MaskedP=GAP,Band=BAND,t1=t1,t2=t2,t3=t3,t4=t4)
                    
                    XLSX.writetable(joinpath(output_dir, param*"_Summary.xlsx"), overwrite=true, 
                        df_accuracy_n2 = (collect(DataFrames.eachcol(df_accuracy_n2)), DataFrames.names(df_accuracy_n2)),
                        df_time_n2 = (collect(DataFrames.eachcol(df_time_n2)), DataFrames.names(df_time_n2)),   
                    )
                else
                    df_accuracy_minNx = DataFrame(MaskedP=GAP,Band=BAND,l0=v1,l2=v2,l4=v3,l6=v4)
                    df_time_minNx = DataFrame(MaskedP=GAP,Band=BAND,t1=t1,t2=t2,t3=t3)
                    
                    XLSX.writetable(joinpath(output_dir, param*"_Summary.xlsx"), overwrite=true, 
                        df_accuracy_minNx = (collect(DataFrames.eachcol(df_accuracy_minNx)), DataFrames.names(df_accuracy_minNx)),
                        df_time_minNx = (collect(DataFrames.eachcol(df_time_minNx)), DataFrames.names(df_time_minNx)),   
                    )
                end                      
            end   
        end
    end     

    for folder in folders
        files2 = [] 
        for imgType in ["L7", "S2"]
            for pixels_to_be_filled in ["_nan","_cloud"]
                if isdir(joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled)) == true 
                    for file in readdir(joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled))
                        for DS in ["_UV_", "_BV_"]
                            if occursin("errorResults"*DS*imgType*pixels_to_be_filled*"_rep"*string(nb_rep2)*"_.csv",file)
                                push!(files2, joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled,file))
                            end                        
                        end
                    end
                end
            end
        end
        
        for file in files2
            df1 = CSV.read(file, DataFrame)

            LOCATION = String[]
            DATE_REF_IMG = String[]
            DATE_TARGET_IMG = String[]    
            MASKP = Float64[]
            BAND = Int64[]

            Factor1 = String[]
            Factor2 = Int64[]
            Factor3 = Float64[]
            Factor4 = Float64[]
            Factor5 = Int64[]
            Factor6 = Float64[]
            Factor7 = Float64[]

            ACCURACY_MSE = Float64[]
            ACCURACY_MSLE = Float64[]
            ACCURACY_R2 = Float64[]
            N_missingPixels = Float64[]
            N_filledPixels = Float64[]

            TOTAL_ELAPSED_TIME_SEC_PER_BAND = Float64[]
            ELAPSED_TIME_SEC_PER_BAND = Float64[]
            ELAPSED_TIME_MIN_PER_BAND = Float64[]
            ELAPSED_TIME_HOURS_PER_BAND = Float64[]
            ELAPSED_TIME_DAYS_PER_BAND = Float64[]
            for row in 1:nb_rep2:length(df1.R2)
                push!(LOCATION,df1.Location[row])
                push!(DATE_REF_IMG,string(df1.Date_Ref[row]))
                push!(DATE_TARGET_IMG,string(df1.Date_Target[row]))       
                push!(MASKP,df1.MaskedP[row])
                push!(BAND,df1.Band[row])               

                push!(Factor1,string(df1.n1[row]))
                push!(Factor2,df1.n2[row])
                push!(Factor3,df1.t1[row])
                push!(Factor4,df1.f[row])
                push!(Factor5,df1.minNx[row])
                if (occursin("_BV_",file))
                    push!(Factor6,df1.w1[row])
                    push!(Factor7,df1.w2[row]) 
                end

                push!(ACCURACY_MSLE, mean(df1.MSLE[row:row+nb_rep2-1]))
                push!(ACCURACY_MSE, mean(df1.MSE[row:row+nb_rep2-1]))
                push!(ACCURACY_R2, mean(df1.R2[row:row+nb_rep2-1]))
                push!(N_missingPixels, mean(df1.UP[row:row+nb_rep2-1]))
                push!(N_filledPixels, mean(df1.KP[row:row+nb_rep2-1]))

                push!(TOTAL_ELAPSED_TIME_SEC_PER_BAND, mean(df1.Total_Time_Sec[row:row+nb_rep2-1]))
                push!(ELAPSED_TIME_SEC_PER_BAND, mean(df1.Time_Sec[row:row+nb_rep2-1]))
                push!(ELAPSED_TIME_MIN_PER_BAND, mean(df1.Time_Min[row:row+nb_rep2-1]))
                push!(ELAPSED_TIME_HOURS_PER_BAND, mean(df1.Time_Hours[row:row+nb_rep2-1]))
                push!(ELAPSED_TIME_DAYS_PER_BAND, mean(df1.Time_Days[row:row+nb_rep2-1]))
            end

            output_dir = replace(file, "_rep"*string(nb_rep2)*"_" => "_rep"*"AVG")    
            if (occursin("_UV_",file))
                df2 = DataFrame(Location=LOCATION, 
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
                CSV.write(output_dir, df2)
            elseif (occursin("_BV_",file))
                df2 = DataFrame(Location=LOCATION, 
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
                CSV.write(output_dir, df2)
            end
        end   
    end

    for folder in folders
        files3 = [] 
        files4 = [] 
        for imgType in ["L7", "S2"]
            for pixels_to_be_filled in ["_nan","_cloud"]
                if isdir(joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled)) == true 
                    for file in readdir(joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled))
                        if occursin("errorResults_UV_"*imgType*pixels_to_be_filled*"_repAVG.csv",file)
                            push!(files3, joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled,file))
                        elseif occursin("errorResults_BV_"*imgType*pixels_to_be_filled*"_repAVG.csv",file)
                            push!(files4, joinpath(raw"C:\GF_Algorithm", imgType, "Analysis",folder, pixels_to_be_filled,file))
                        end
                    end
                end
            end
        end    
            
        for i in 1:length(files3)
            df_UV = CSV.read(files3[i], DataFrame)
            df_BV = CSV.read(files4[i], DataFrame)       
            
            DIFF_TIME = [100*((df_BV.Total_Time_Sec[row]-df_UV.Total_Time_Sec[row])/df_UV.Total_Time_Sec[row])  for row in 1:length(df_UV.R2)]
            output_dir = replace(files3[i], "repAVG" => "UV_BV_Comparison")    
            df3 = DataFrame(MaskedP=df_UV.MaskedP, 
                        Band=df_UV.Band, 
                        R2_UV=df_UV.R2,
                        R2_BV=df_BV.R2,
                        Time_UV=df_UV.Total_Time_Sec, 
                        Time_BV=df_BV.Total_Time_Sec,
                        Diff_Time=DIFF_TIME
                        )
            CSV.write(output_dir, df3)   
        end   
    end

    for folder in folders
        files3 = [] 
        files4 = [] 
        for imgType in ["L7", "S2"]
            for pixels_to_be_filled in ["_nan","_cloud"]
                input_dir3 = joinpath(raw"C:\GF_Algorithm", imgType, "BiasImg",folder, pixels_to_be_filled)
                outputHist = joinpath(raw"C:\GF_Algorithm", imgType, "Histograms",folder, pixels_to_be_filled)
                mkpath(outputHist)

                if isdir(input_dir3) == true               
                    xmin = 999
                    xmax = -999
                    
                    for sub_input_dir3 in readdir(input_dir3; join=true) 
                        fileName = reverse(reverse(sub_input_dir3)[5:22])*".png"
                        arrBias = GeoArrays.read(joinpath(outputHist, sub_input_dir3), band=1)
                        v1 = filter(x -> x !=0.0, vcat(arrBias...)) 
                        if minimum(v1)<xmin
                            xmin = minimum(v1)
                        end
                        if maximum(v1)>xmax
                            xmax = maximum(v1)
                        end
                    end
                    for sub_input_dir3 in readdir(input_dir3; join=true) 
                        fileName = reverse(reverse(sub_input_dir3)[5:22])*".png"
                        arrBias = GeoArrays.read(joinpath(outputHist, sub_input_dir3), band=1)
                        v1 = filter(x -> x !=0.0, vcat(arrBias...)) 
                        plotHist = histogram(v1, ylabel = "Number of Pixels", xlabel = "Difference %", color=:gray, legend=false, size = (600, 400))
                        xlims!(xmin, xmax)
                        ylims!(0, 10000)
                        savefig(plotHist,joinpath(outputHist, fileName))
                    end
                end
            end
        end
    end
end;


