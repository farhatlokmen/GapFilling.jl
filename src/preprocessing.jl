

function set_boundary(value::Int64, Bmin::Int64, Bmax::Int64)
    """Force any value to be within the boundaries of an array""" 
    
    if value<=Bmin
        return Bmin    
    elseif value>=Bmax
        return Bmax
    else
        return value
    end
end
function rename_resize_image(begin1::Int64, end1::Int64, begin2::Int64, end2::Int64, file::String, imgType::String,sub_output_dir::String)
    """rename an image file name & resize array to user specified range""" 
    
    fileName = reverse(reverse(file)[5:12])*"_"*reverse(reverse(file)[14:19])*".tif"
    if imgType == "S2"
        fileName = reverse(reverse(file)[35:42])*"_"*reverse(reverse(file)[5:10])*".tif"
    end  
    arr = GeoArrays.read(file)
    begin1 = set_boundary(begin1,1,size(arr)[1])
    begin2 = set_boundary(begin2,1,size(arr)[2])
    end1 = set_boundary(end1,begin1,size(arr)[1])
    end2 = set_boundary(end2,begin2,size(arr)[2])
    subset = arr[begin1:end1,begin2:end2,begin:end]    
    return GeoArrays.write!(joinpath(sub_output_dir,fileName), subset) 
end

function get_cloudP(file::String, imgType::String)
    """get summary info of imagery"""  
    
    if imgType == "L7"
        # 1, 9 - DDV (dark dense vegetation)
        # 2, 34 - cloud
        # 4, 12, 20, 36, 52 - cloud shadow
        # 8, 12, 24, 40, 56 - adjacent to cloud
        # 16, 20, 24, 48, 52, 56 - snow
        # 32, 34, 36, 40, 48, 52, 56 - water
        arr6 = GeoArrays.read(file,band=6)
        arr7 = GeoArrays.read(file,band=7)
        n_cloudPixels = length(Tuple.(findall(x->(x==2)||(x==34)||(x==4)||(x==12)||(x==20)||(x==36)||(x==52), arr7))) 
        nb_nan_pixels = length(Tuple.(findall(x->(isnan(x)), arr6)))
        percent_cloudPixels = 100*round(n_cloudPixels/(size(arr6)[1]*size(arr6)[2]); digits=4) # count based on total pixels in image valid or not
        percent_nan_pixels = 100*round(nb_nan_pixels/(size(arr6)[1]*size(arr6)[2]); digits=4)
        return percent_cloudPixels, percent_nan_pixels, file
    else        
        # 0 - No data
        # 1 - Saturated / Defective
        # 2 - Dark Area Pixels
        # 3 - Cloud Shadows
        # 4 - Vegetation
        # 5 - Bare Soils
        # 6 - Water
        # 7 - Clouds low probability / Unclassified
        # 8 - Clouds medium probability
        # 9 - Clouds high probability
        # 10 - Cirrus
        # 11 - Snow / Ice
        arr10 = GeoArrays.read(file,band=10)        
        arr11 = GeoArrays.read(file,band=11)
        n_cloudPixels = length(Tuple.(findall(x->(x==3)||(x==8)||(x==9)||(x==10), arr11))) # cloud_pixels = [3, 8, 9, 10]
        nb_nan_pixels = length(Tuple.(findall(x->(isnan(x)), arr10)))        
        percent_cloudPixels = 100*round(n_cloudPixels/((size(arr11)[1]*size(arr11)[2])); digits=4) # count based on total pixels in image valid or not
        percent_nan_pixels = 100*round(nb_nan_pixels/(size(arr11)[1]*size(arr11)[2]); digits=4)
        return percent_cloudPixels, percent_nan_pixels, file
    end
end
function get_path(date,files)
    """extract dir file path based on image acquisition date"""
    
    P = "No_Image"
    for file in files
        if occursin(string(date), file) && occursin("tif", file)
            P = file
        end
    end
    return P
end

function plotHistogram(values::Vector{Float64}, output_dir::String, pixels_to_be_filled::String) 
    plotHist1 = histogram(values, ylabel = "Number of Days", xlabel = "Gap %", color=:gray, legend=false, size = (600, 400))
    xlims!(minimum(values), maximum(values))
    ylims!(0, 100)
    savefig(plotHist1,joinpath(output_dir, "gapTemporalVariability"*pixels_to_be_filled*".png"))
end

function select_Img(df1::DataFrame, location::String, values::Vector{Int64})
    """extract info regarding cloud or nan pixels number"""
    
    LOCATION = String[]
    CLOUD_P = Float64[]
    NAN_P = Float64[]
    DATE = Int64[]
    tempDF = sort(df1, [:nanP,:CloudP]) 
    push!(LOCATION, location)
    push!(CLOUD_P, tempDF.CloudP[1])
    push!(NAN_P, tempDF.nanP[1])
    push!(DATE, tempDF.Date[1])
    
    for index in 1:(length(values)-1)                    
        filtered_df = df1[(df1.Location.==location) .& (df1.CloudP.>values[index]) .& (df1.CloudP.<=values[index+1]), :]
        if size(filtered_df)[1] >= 1
            Random.seed!(seed) 
            idx = rand(1:length(filtered_df.CloudP))
            push!(LOCATION, location)
            push!(CLOUD_P, filtered_df.CloudP[idx])
            push!(NAN_P, filtered_df.nanP[idx])
            push!(DATE, filtered_df.Date[idx])                
        end 
    end    
    return LOCATION,CLOUD_P,NAN_P,DATE
end

function getInfo(case,index,filtered_filePath,filtered_missingP)
    if case == "before"
        condition = (index>0)
        if condition
            return replace(filtered_filePath[index], "Step1" => "Step4"), filtered_missingP[index][1]
        else
            return "No_Image", 404.0
        end 
    else
        condition = (index<(length(filtered_filePath)+1))
        if condition
            return replace(filtered_filePath[index], "Step1" => "Step4"), filtered_missingP[index][1]
        else
            return "No_Image", 404.0
        end
    end
end;

function preprocess(path, imgType, folders, pixels_to_be_filled, values, seed, x, w)
    dir = joinpath(path, imgType)
    y = x+w
    for folder in folders    
        ####################################
        ### Step1: Rename & Resize image ###
        ####################################
        input_dir1 = joinpath(dir, "Donloaded", folder)    
        output_dir1 = joinpath(dir, "Preprocessed", "Step1", folder, pixels_to_be_filled)
        rm(output_dir1, force=true, recursive=true)  
        mkpath(output_dir1)     
        for file in readdir(input_dir1; join=true) 
            rename_resize_image(x, y, x, y, file, imgType,output_dir1) 
        end
        
        ##############################################################################
        ### Step2: Count UP per image & Select reference and SCL images to be used ###
        ##############################################################################
        input_dir2 = joinpath(dir, "Preprocessed", "Step1", folder, pixels_to_be_filled)
        output_dir2 = joinpath(dir, "Preprocessed", "Step2", folder, pixels_to_be_filled)
        rm(output_dir2, force=true, recursive=true)      
        mkpath(output_dir2) 
        DATE=Int[]
        LOCATION=String[]
        CLOUD_PIXELS=Float64[]
        NAN_PIXELS=Float64[]
        PATH = String[]
        for file in readdir(input_dir2; join=true) 
            cloudP, nanP, path = get_cloudP(file, imgType)  
            push!(DATE, parse(Int,reverse(reverse(file)[12:19])))
            push!(LOCATION, folder)
            push!(CLOUD_PIXELS, cloudP)
            push!(NAN_PIXELS, nanP)
            push!(PATH, path)
        end
        id_rows_delete = []
        checked_rows = []
        for i in 1:length(DATE)
            for t in [j for j in 1:length(DATE) if (j!=i) & !(j in checked_rows)]
                if (DATE[i] == DATE[t]) & (LOCATION[i] == LOCATION[t])            
                    if NAN_PIXELS[i] > NAN_PIXELS[t] 
                        push!(id_rows_delete,i)
                    else
                        push!(id_rows_delete,t)
                    end
                end        
            end
            push!(checked_rows, i)
        end
        df1 = DataFrame(Date=DATE, Location=LOCATION, CloudP=CLOUD_PIXELS, nanP=NAN_PIXELS, filePath=PATH)
        delete!(df1, id_rows_delete)
        CSV.write(joinpath(output_dir2, "gapP.csv"), df1);
        cP = df1[(df1.Location.==folder), :].CloudP
        nP = df1[(df1.Location.==folder), :].nanP
        plotHistogram(cP, output_dir2, "_cloud")
        plotHistogram(nP, output_dir2, "_nan")

        
        LOCATION = String[]
        CLOUD_P = Float64[]
        NAN_P = Float64[]
        DATE = Int64[]
        if imgType == "L7"
            df1_before_failure = df1[(df1.Location.==folder) .& (df1.Date.<20030531) .& (df1.nanP.<0.1), :]              
            if pixels_to_be_filled == "_nan"
                if size(df1_before_failure)[1] >= 1
                    tempDF_before = sort(df1_before_failure , [:nanP, :CloudP]) 
                    push!(LOCATION, folder)
                    push!(CLOUD_P, tempDF_before.CloudP[1])
                    push!(NAN_P, tempDF_before.nanP[1])
                    push!(DATE, tempDF_before.Date[1]) 
                    for index in 1:(length(values)-1)
                        df1_after_failure = df1[(df1.Location.==folder) .& (df1.Date.>=20030531) .& (df1.nanP.>values[index]) .& (df1.nanP.<=values[index+1]), :]
                        if size(df1_after_failure)[1] >= 1
                            tempDF_after = sort(df1_after_failure, [:CloudP]) 
                            push!(LOCATION, folder)
                            Random.seed!(seed) 
                            r = rand(1:length(tempDF_after.CloudP))                        
                            push!(CLOUD_P, tempDF_after.CloudP[r])
                            push!(NAN_P, tempDF_after.nanP[r])
                            push!(DATE, tempDF_after.Date[r]) 
                        end
                    end                
                end
            else
                v1,v2,v3,v4 = select_Img(df1_before_failure,folder,values)
                append!(LOCATION, v1)
                append!(CLOUD_P, v2)
                append!(NAN_P, v3)
                append!(DATE, v4) 
            end
        else        
            v1,v2,v3,v4 = select_Img(df1[(df1.Location.==folder), :],folder,values)
            append!(LOCATION, v1)
            append!(CLOUD_P, v2)
            append!(NAN_P, v3)
            append!(DATE, v4) 
        end
        df2 = DataFrame(Location=LOCATION, CloudP=CLOUD_P, nanP=NAN_P, Date=DATE)
        CSV.write(joinpath(output_dir2, "selectedImages.csv"), df2);    
        
        ###################################################
        ### Step3: Generate artificially masked images ###
        ###################################################
        input_dir3 = joinpath(dir, "Preprocessed", "Step1", folder, pixels_to_be_filled)
        output_dir3 = joinpath(dir, "Preprocessed", "Step3", folder, pixels_to_be_filled)
        rm(output_dir3, force=true, recursive=true)     
        mkpath(output_dir3)  
            
        date_ref = string(df2[(df2.Location.==folder), :].Date[1])
        for file1 in readdir(input_dir3; join=true)
            if occursin(date_ref, file1)
                selected_ref = GeoArrays.read(file1)
                for k in 2:(size(df2[(df2.Location.==folder), :].Date)[1])
                    date_selectedImg = string(df2[(df2.Location.==folder) , :].Date[k])                
                    for file2 in readdir(input_dir3; join=true)
                        if occursin(date_selectedImg, file2)
                            if imgType == "L7"
                                selected_img = GeoArrays.read(file2, band=6)
                                coord_UP = Tuple.(findall(x->(isnan(x)), selected_img))
                                if pixels_to_be_filled == "_cloud"  
                                    selected_img = GeoArrays.read(file2, band=7)                           
                                    coord_UP = Tuple.(findall(x->(x==2)||(x==34)||(x==4)||(x==12)||(x==20)||(x==36)||(x==52), selected_img)) 
                                end                            
                                temp_arr = deepcopy(selected_ref)
                                for coord in coord_UP
                                    temp_arr[coord[1],coord[2],1:6].=-99
                                end                            
                                fileName = reverse(reverse(file2)[begin:19])
                                GeoArrays.write!(joinpath(output_dir3, fileName), temp_arr)
                            else
                                selected_img = GeoArrays.read(file2, band=11)
                                coord_UP = Tuple.(findall(x->(x==3)||(x==8)||(x==9)||(x==10), selected_img))
                                temp_arr = deepcopy(selected_ref)
                                for coord in coord_UP
                                    temp_arr[coord[1],coord[2],1:10].=-99
                                end
                                fileName = reverse(reverse(file2)[begin:19])
                                GeoArrays.write!(joinpath(output_dir3, fileName), temp_arr)
                            end
                        end
                    end                
                end            
            end
        end;
        
        ######################################################
        ### Step4: Apply SCL cloud masks on all S2 imagery ###
        ######################################################
        input_dir4 = joinpath(dir, "Preprocessed", "Step1", folder, pixels_to_be_filled)
        output_dir4 = joinpath(dir, "Preprocessed", "Step4", folder, pixels_to_be_filled)
        rm(output_dir4, force=true, recursive=true)      
        mkpath(output_dir4)  
        
        for file in readdir(input_dir4; join=true)
            if imgType == "L7"
                img = GeoArrays.read(file)  
                scl = GeoArrays.read(file, band=6)
                coord_UP = Tuple.(findall(x->(isnan(x)), scl))
                if pixels_to_be_filled == "_cloud"  
                    scl = GeoArrays.read(file, band=7)                           
                    coord_UP = Tuple.(findall(x->(x==2)||(x==34)||(x==4)||(x==12)||(x==20)||(x==36)||(x==52), scl)) 
                end           
                for coord in coord_UP
                    img[coord[1],coord[2],1:6].=-99
                end        
                fileName = reverse(reverse(file)[begin:19])
                GeoArrays.write!(joinpath(output_dir4, fileName), img)
            else
                img = GeoArrays.read(file)  
                scl = GeoArrays.read(file, band=11)
                coord_UP = Tuple.(findall(x->(x==3)||(x==8)||(x==9)||(x==10), scl))
                for coord in coord_UP
                    img[coord[1],coord[2],1:10].=-99
                end        
                fileName = reverse(reverse(file)[begin:19])
                GeoArrays.write!(joinpath(output_dir4, fileName), img)
            end
        end
        println("Step4: Apply SCL cloud masks on all S2 imagery >> Done!")
        
        input_dir5_1 = joinpath(dir, "Preprocessed", "Step3", folder, pixels_to_be_filled)
        input_dir5_2 = joinpath(dir, "Preprocessed", "Step4", folder, pixels_to_be_filled)
        output_dir5 = joinpath(dir, "Preprocessed", "Step5", folder, pixels_to_be_filled)
        rm(output_dir5, force=true, recursive=true)     
        mkpath(output_dir5) 
        
        LOCATION = String[]
        REFERENCE = String[]
        pathBEFORE1 = String[]
        pathBEFORE2 = String[]
        pathTARGET = String[]
        pathAFTER1 = String[]
        pathAFTER2 = String[]
        imgBEFORE1 = Float64[]
        imgBEFORE2 = Float64[]
        imgTARGET = Float64[]
        imgAFTER1 = Float64[]
        imgAFTER2 = Float64[]
        
        files1 = readdir(input_dir5_2; join=true)  
        files2 = readdir(input_dir5_1; join=true)  
        date_cloudFree = df2[(df2.Location.==folder), :].Date[1]
        refPath = get_path(date_cloudFree,files1)    
        masked_dates = df2[(df2.Location.==folder), :].Date[2:end]  
        for masked_date in masked_dates
            df1 =sort(df1, [:Date])
            filtered_dates = df1[(df1.Location.==folder), :].Date
            filtered_filePath = df1[(df1.Location.==folder), :].filePath
            filtered_missingP = []
            for i in 1:length(filtered_dates) 
                push!(filtered_missingP,df1[(df1.Location.==folder), :].CloudP[i]+df1[(df1.Location.==folder), :].nanP[i])
            end

            idx = 0
            for counter in 1:length(filtered_dates)
                if masked_date==filtered_dates[counter]
                    idx = counter
                end
            end
            push!(LOCATION, folder)
            push!(REFERENCE, refPath)        
            push!(pathTARGET, get_path(masked_date,files2))                
            push!(imgTARGET, filtered_missingP[idx][1])

            p,m = getInfo("before",idx-1,filtered_filePath,filtered_missingP)
            push!(pathBEFORE1, p) 
            push!(imgBEFORE1, m) 
            p,m = getInfo("before",idx-2,filtered_filePath,filtered_missingP)
            push!(pathBEFORE2, p) 
            push!(imgBEFORE2, m)
            p,m = getInfo("after",idx+1,filtered_filePath,filtered_missingP)
            push!(pathAFTER1, p) 
            push!(imgAFTER1, m)
            p,m = getInfo("after",idx+2,filtered_filePath,filtered_missingP)
            push!(pathAFTER2, p) 
            push!(imgAFTER2, m)
        end 

        df3 = DataFrame(Location=LOCATION, 
                        Ref_Img=REFERENCE, 
                        Img_Before2=pathBEFORE2, 
                        Img_Before1=pathBEFORE1, 
                        Target_Img=pathTARGET, 
                        Img_After1=pathAFTER1, 
                        Img_After2=pathAFTER2,
                        missingPixelsBefore2=imgBEFORE2, 
                        missingPixelsBefore1=imgBEFORE1, 
                        missingPixelsTarget=imgTARGET, 
                        missingPixelsAfter1=imgAFTER1, 
                        missingPixelsAfter2=imgAFTER2)
        CSV.write(joinpath(output_dir5, "pathsData.csv"), df3);    
    end
end

