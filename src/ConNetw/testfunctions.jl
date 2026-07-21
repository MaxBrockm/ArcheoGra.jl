"""
    ArcheoGra_graph_comparison()

Test function creating Haversine graphs and testing graph comparison functions. 
"""
function ArcheoGra_graph_comparison()
    file_magda_long = file_props(
        filePath = joinpath("src", "DataFolder", "Rohmaterial_Datenbank_Magdalenien_copy.xlsx"),
        file_type="xlsx", sheetname ="Tabelle1", material_cols=8:115, 
        material_age = "Lang", material_distance = 250, long_ind = 3, lat_ind = 4,
        name_ind = 2, class_ind = 5
    )
    file_magda_short = file_props(
        filePath = joinpath("src", "DataFolder", "Rohmaterial_Datenbank_Magdalenien_copy.xlsx"),
        file_type="xlsx", sheetname ="Tabelle1", material_cols=8:115, 
        material_age = "Kurz", material_distance = 150, long_ind = 3, lat_ind = 4,
        name_ind = 2, class_ind = 5
    )

    data_test = Arch_DataSet(Vector{ArchDataEntry}(undef, 0)); MaterialData_test=MaterialDataList(Dict{String, MaterialEntry}())
    data_test, MaterialData_test = read_data(file_magda_long; data=data_test, MaterialData=MaterialData_test)
    data_test, MaterialData_test = read_data(file_magda_short; data=data_test, MaterialData=MaterialData_test)

    graphen_test = Haversine_connect_nodes(data_test, MaterialData_test);
   
    adj_1 = graphen_test[1][2]; w_adj_1 = graphen_test[1][3]
    adj_2 = graphen_test[2][2]; w_adj_2 = graphen_test[2][3]

    long = [entry.coords[1] for entry in graphen_test[1][1]]
    lat = [entry.coords[2] for entry in graphen_test[1][1]]
    cluster_comp_plot(w_adj_1, w_adj_2, (long, lat))
    @info "Start produce results"
    produce_results(w_adj_1)
    w_lap_mat1, _ = get_matrices(w_adj_1)
    w_lap_mat2, _ = get_matrices(w_adj_2)

    @show effective_resistance_long = effective_resistance(w_lap_mat1)
    @show effective_resistance_short = effective_resistance(w_lap_mat2)

    produce_results(w_adj_1)
end


"""
    Output_function_1104(;compare=false)

Test function for Haversine-networks. Parameter compare decides wether plots include comparison of material- to non-material-graphs-
"""
function Output_function_1104(;compare=false)
    try 
        mkdir(joinpath("src", "OutputFolder", string(Dates.today()))) 
    catch 
        @info "Overwriting existing files"
    end
    file_fullLBK = file_props(
        filePath = joinpath("src", "DataFolder", "fullLBK_ds02preuss_GrdPts_belegtklein.txt"),
        file_type="txt", material_cols=6:9, delim = ',',
        material_age = "fullLBK", material_distance = 350*ones(length(6:9)), long_ind = 2, lat_ind = 3
    )

    file_oldLBK = file_props(
        filePath = joinpath("src", "DataFolder", "oLBK_ds02_belegt.txt"),
        file_type="txt",  material_cols=6:9, 
        material_age = "oldLBK", material_distance = 350*ones(length(6:9)), long_ind = 2, lat_ind = 3,
        delim = '\t'
    )

    file_Magdafinal = file_props(
        filePath = joinpath("src", "DataFolder", "final_Magda_new.csv"), file_type="csv",  material_cols= 4:39, 
        material_age = "Magdafinal", material_distance = 350*ones(length(4:39)), long_ind = 2, lat_ind = 3,
        delim = ','
    )
    files = [file_fullLBK, file_oldLBK, file_Magdafinal]#, file_magda, file_Magdafinal, file_MagdaGrid]
    for file in files
        try 
            mkdir(joinpath("src", "OutputFolder", string(Dates.today()), file.material_age)) 
        catch 
            @info "Overwriting existing files"
        end
    end
    
    data_test = Arch_DataSet(Vector{ArchDataEntry}(undef, 0)); MaterialData_test=MaterialDataList(Dict{String, MaterialEntry}())
    data_test, MaterialData_test = read_data(file_Magdafinal; data=data_test, MaterialData=MaterialData_test)
    data_test, MaterialData_test = read_data(file_fullLBK; data=data_test, MaterialData=MaterialData_test)
    data_test, MaterialData_test = read_data(file_oldLBK; data=data_test, MaterialData=MaterialData_test)
    if compare
        graphen_test_compare = Haversine_connect_nodes(data_test, MaterialData_test);
    else
        graphen_test = Haversine_connect_nodes(data_test, MaterialData_test, Plot_compare=false);
    end
end

"""
    test_LCPAdistance()

Test function for LCPA networks.
"""
function test_LCPAdistance()
    try 
        mkdir(joinpath("src", "OutputFolder", string(Dates.today()))) 
    catch 
        @info "Overwriting existing files"
    end
    file_fullLBK = file_props(
        filePath = joinpath("src", "DataFolder", "fullLBK_ds02preuss_GrdPts_belegtklein.txt"),
        file_type="txt", material_cols=6:9, delim = ',',
        material_age = "fullLBK", material_distance = 350*ones(length(6:9)), long_ind = 2, lat_ind = 3
    )

    file_oldLBK = file_props(
        filePath = joinpath("src", "DataFolder", "oLBK_ds02_belegt.txt"),
        file_type="txt", material_cols=6:9, 
        material_age = "oldLBK", material_distance = 350*ones(length(6:9)), long_ind = 2, lat_ind = 3,
        delim = '\t'
    )
    files = [file_fullLBK, file_oldLBK]#, file_magda, file_Magdafinal, file_MagdaGrid]
    for file in files
        try 
            mkdir(joinpath("src", "OutputFolder", string(Dates.today()), file.material_age)) 
        catch 
            @info "Overwriting existing files"
        end
    end
    
    data_test = Arch_DataSet(Vector{ArchDataEntry}(undef, 0)); MaterialData_test=MaterialDataList(Dict{String, MaterialEntry}())
    # data_test, MaterialData_test = read_data(file; data=data_test, MaterialData=MaterialData_test)
    data_test, MaterialData_test = read_data(file_oldLBK; data=data_test, MaterialData=MaterialData_test)
    data_test, MaterialData_test = read_data(file_fullLBK; data=data_test, MaterialData=MaterialData_test)
    k = 16; ℓ = 16
    connect_nodes_LCPA(data_test, MaterialData_test, k=k, ℓ=ℓ, PlotOut=false, conn_func=Tobler_Hiking_mod)
    return nothing
end


"""
    test_LCPAdistance_nonsymm

Function to test (legacy) non-symmetric LCPA network creation.
"""
function test_LCPAdistance_nonsymm()
    try 
        mkdir(joinpath("src", "OutputFolder", string(Dates.today()))) 
    catch 
        @info "Overwriting existing files"
    end
    file_fullLBK = file_props(
        filePath = joinpath("src", "DataFolder", "fullLBK_ds02preuss_GrdPts_belegtklein.txt"),
        file_type="txt", material_cols=6:9, delim = ',',
        material_age = "fullLBK", material_distance = 350*ones(length(6:9)), long_ind = 2, lat_ind = 3
    )

    file_oldLBK = file_props(
        filePath = joinpath("src", "DataFolder", "oLBK_ds02_belegt.txt"),
        file_type="txt",  material_cols=6:9, 
        material_age = "oldLBK", material_distance = 350*ones(length(6:9)), long_ind = 2, lat_ind = 3,
        delim = '\t'
    )
    files = [file_fullLBK, file_oldLBK]#, file_magda, file_Magdafinal, file_MagdaGrid]
    for file in files
        try 
            mkdir(joinpath("src", "OutputFolder", string(Dates.today()), file.material_age, "nonsymm")) 
        catch 
            @info "Overwriting existing files"
        end
    end
    
    data_test = Arch_DataSet(Vector{ArchDataEntry}(undef, 0)); MaterialData_test=MaterialDataList(Dict{String, MaterialEntry}())
    # data_test, MaterialData_test = read_data(file; data=data_test, MaterialData=MaterialData_test)
    data_test, MaterialData_test = read_data(file_oldLBK; data=data_test, MaterialData=MaterialData_test)
    # data_test, MaterialData_test = read_data(file_fullLBK; data=data_test, MaterialData=MaterialData_test)
    k = 16; ℓ = 16
    connect_nodes_LCPA_non_symm(data_test, MaterialData_test, k=k, ℓ=ℓ, PlotOut=false, conn_func=mod_mod_T_Hiking)
    return nothing
end

"""
    test_crosstime()

Test function for LCPA networks.
"""
function test_crosstime()
    try 
        mkdir(joinpath("src", "OutputFolder", string(Dates.today()))) 
    catch 
        @info "Overwriting existing files"
    end
    file_fullLBK = file_props(
        filePath = joinpath("src", "DataFolder", "fullLBK_ds02preuss_GrdPts_belegtklein.txt"),
        file_type="txt", material_cols=6:9, delim = ',',
        material_age = "fullLBK", material_distance = 350*ones(length(6:9)), long_ind = 2, lat_ind = 3
    )

    file_oldLBK = file_props(
        filePath = joinpath("src", "DataFolder", "oLBK_ds02_belegt.txt"),
        file_type="txt",  material_cols=6:9, 
        material_age = "oldLBK", material_distance = 350*ones(length(6:9)), long_ind = 2, lat_ind = 3,
        delim = '\t'
    )
        try 
            mkdir(joinpath("src", "OutputFolder", string(Dates.today()), "cross_time")) 
        catch 
            @info "Overwriting existing files"
        end
    
    data_test = Arch_DataSet(Vector{ArchDataEntry}(undef, 0)); MaterialData_test=MaterialDataList(Dict{String, MaterialEntry}())
    data_test, MaterialData_test = read_data(file_oldLBK; data=data_test, MaterialData=MaterialData_test)
    data_test, MaterialData_test = read_data(file_fullLBK; data=data_test, MaterialData=MaterialData_test)
    k = 16; ℓ = 16
    create_connections(data_test, ["oldLBK", "fullLBK"], k=k, ℓ=ℓ, conn_func=Tobler_Hiking_mod)
    return nothing
end