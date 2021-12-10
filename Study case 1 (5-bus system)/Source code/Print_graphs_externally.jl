#Project: A lossy DC Optimal Power Flow with an accuratemodeling of Li-Ion batteries
#Student: P. Serna Torre
#Couse: MAE207 / Professor: M. Davidson

#Warning: Paths used here have Linux format (i.e. /..../), which is different to Windows format ("\\.....\\")

println(string("file.jl to print graphs (in case the main program is run in the server)"))
println(string("  "))
println(string("      .................. #  # ...................      "))
println(string("      Lossy DCOPF with an accurate BESS modelling      "))
println(string("      Student: Paul Serna-Torre"))
println(string("      Course: MAE207 UCSD, Professor: M. Davidson"))
println(string("      Note: The formulation was implemented in Julia 1.6.3 under Ubuntu 20.04"))
println(string("      .................. #  # ...................      "))
println(string("  "))
println(string("  "))


#Packages
using DataFrames, CSV, CPLEX, JuMP, LinearAlgebra, Dates, Plots, VegaLite
println(string("Julia packages were loaded successfully"))

println(string("Reading input files (csv) and creating dataframes and arrays. . ."))

#Configure the fontfamily of the plots
default(; # Plots defaults
fontfamily="Computer modern",
label="" # only explicit legend entries
)

#Load directory path
paths=CSV.read("directory_paths.csv", DataFrame);
main_path=paths[1,2];
input_folder=paths[2,2];
output_folder=paths[3,2];
renewable_folder=paths[4,2];

#Main folder
folder_name=main_path;

#Input data folder
path_inputdata=string(folder_name, "/",input_folder,"/");

#Output data folder
path_output=string(folder_name, "/",output_folder, "/");

#Load Configuration.csv
config=CSV.read(string(path_inputdata,"configuration.csv"), DataFrame);
start_case=config[1,2];
final_case=config[2,2];
date_format=config[3,2];
resolution=parse(Int,config[4,2]);
price_p_shed=parse(Float64,config[5,2]);
tolerance_solver=parse(Float64,config[6,2]);
cores_solver=parse(Float64, config[7,2]);
MVAbase=parse(Float64, config[8,2]);
digits_round=parse(Int, config[9,2]);
iter_losses_max=parse(Int, config[10,2]);
tol_losses=parse(Float64, config[11,2]);
bess_modeling=parse(Int, config[12,2]);

#Set main parameters
start_case=Date(start_case,date_format);
final_case=Date(final_case,date_format);
ndays=Dates.value(final_case-start_case)+1;
deltaT=resolution/60;
nper=Int(24/deltaT);

#Load bus.csv
bus=CSV.read(string(path_inputdata,"bus.csv"), DataFrame)[:,1:3];
nbuses=size(bus,1);
insertcols!(bus, 1, :order => 1:nbuses);

#Create a time chart.csv to know the time horizon and duration of each period
time_chart=Array{Any}(nothing,ndays*nper+1,8);
time_chart[1,1:8]=["nday" "nper" "year" "month" "day" "hour" "min" "date"];

    for d in 1:ndays, t in 1:nper    
            time_chart[t+1+(d-1)*nper,1]=d;
            time_chart[t+1+(d-1)*nper,2]=t;
            time_chart[t+1+(d-1)*nper,3]=Dates.year(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,4]=Dates.month(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,5]=Dates.day(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,6]=Dates.hour(DateTime(start_case)+Dates.Minute(resolution*(t-1)));
            time_chart[t+1+(d-1)*nper,7]=Dates.minute(DateTime(start_case)+Dates.Minute(resolution*(t-1)));        
            time_chart[t+1+(d-1)*nper,8]=DateTime(start_case+Dates.Day(d-1))+Dates.Minute(resolution*(t-1));
    end
    time_chart_df = DataFrame(time_chart[2:end,:], :auto);
    rename!(time_chart_df, Symbol.(time_chart[1,:]));
    CSV.write(string(path_output,"time_chart.csv"), time_chart_df);

    time_chart_df[!,:nday] = convert.(Int, time_chart_df[!,:nday]);
    time_chart_df[!,:nper] = convert.(Int, time_chart_df[!,:nper]);
    time_chart_df[!,:year] = convert.(Int, time_chart_df[!,:year]);
    time_chart_df[!,:month] = convert.(Int, time_chart_df[!,:month]);
    time_chart_df[!,:day] = convert.(Int, time_chart_df[!,:day]);
    time_chart_df[!,:hour] = convert.(Int, time_chart_df[!,:hour] );
    time_chart_df[!,:min] = convert.(Int, time_chart_df[!,:min]);

#Upload plant.csv
plant=CSV.read(string(path_inputdata,"plant.csv"), DataFrame)[:,[1,2,9,10,11,27,29,30]];
ngen=size(plant,1);


#Non-dispatchable units
dis_units=plant[ ((plant.type .!= "wind") .& (plant.type .!= "solar")) .& (plant.status .==1), :];
ngendis=size(dis_units,1);
insertcols!(dis_units, 1, :order => 1:ngendis);
gencost(g)=dis_units[g, :GenFuelCost]*dis_units[g, :GenIOB]; 

nondis_units=plant[ (plant.type .== "wind") .| (plant.type .== "solar") .& (plant.status .==1), :];
ngennodis=size(nondis_units,1);
insertcols!(nondis_units, 1, :order => 1:ngennodis);

#Set the path to pull out renewable energy data
path_reprofile=string(path_inputdata, renewable_folder);

#Function that finds the maximum renewable power availability of the generator g on day d at period t
function prnw(g,d,t)
    zone= bus[ bus.bus_agg .== nondis_units[g,:bus_agg], :zone_id][1,1];
    type= nondis_units[g,:type];
    profile=CSV.read(string(path_reprofile,"/",type,"_zone",zone,".csv"),DataFrame,header=4)[:,[2,3]];
    profile.time = DateTime.(profile.local_time, "yyyy-mm-dd H:M");
    value=profile[ (profile.time .== time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1]),:electricity][1,1]*nondis_units[g,:Pmax];
        if ismissing(value)
                return 0
        else
                return value
        end
end

#Load load.csv
load=CSV.read(string(path_inputdata,"load.csv"), DataFrame);
load.Date = Date.(load.Date, "mm/dd/yyyy")

#Function that finds the demand (MW) of the bus n on day d at period t
function pc(n,d,t)  
    global column=0;
    for i in 1:size(names(load),1)
        if names(load)[i]==string(bus[n, :bus_agg])
           global column = i;
           break
        end
   end
   if column !=0 
    value=load[ (load.Date .== Date(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])) .& 
    (load.Time .== Time(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])),column][1,1];    
   else
    value=0
   end

    if ismissing(value)
            return 0
    else
            return value
    end
end 

#Some functions to arrange branch data

#Function that find the branch if (m,n) belongs to branches matrix
function map(m,n,branches)
    a=[m n]; 
    found=0;
    rowfound=0;
    Searchi=Int[ a == [branches[i,1] branches[i,2]] for i=1:size(branches,1) ]
    for i in 1:size(Searchi,1)
        if Searchi[i,1] == 1
            found=1
            rowfound=i
        end
    end
    if found==0
        a=[n m];
        Searchi=Int[ a == [branches[i,1] branches[i,2]] for i=1:size(branches,1) ]
        for i in 1:size(Searchi,1)
            if Searchi[i,1] == 1
                found=1
                rowfound=i
            end
        end
    end
    return rowfound
end

#Function that verify if the branch (m,n) exists
function exist(m,n,branches)
    a=[m n]; 
    found=0;
    rowfound=0;
    Searchi=Int[ a == [branches[i,1] branches[i,2]] for i=1:size(branches,1) ]
    for i in 1:size(Searchi,1)
        if Searchi[i,1] == 1
            found=1
            rowfound=i
        end
    end
  
    return found
end

#Load branch.csv
branch=CSV.read(string(path_inputdata,"branch.csv"), DataFrame)[:,1:8];
branch=branch[branch.status.==1,:];
nbranches=size(branch,1);
insertcols!(branch, 1, :order => 1:nbranches);
mapbranch=zeros(nbranches,2);
for br in 1:nbranches
    mapbranch[br,1:2]=[bus[bus.bus_agg.==branch[br,:from_bus_agg],:order][1] bus[bus.bus_agg.==branch[br,:to_bus_agg],:order][1]];
end
mapbranch=convert.(Int, mapbranch);
r(br)=branch[br,:r];
x(br)=branch[br,:x];
Capacity(br)=branch[br,:rateA];

#Load battery.csv
bess=CSV.read(string(path_inputdata,"battery.csv"), DataFrame);
nbess=size(bess,1);
insertcols!(bess, 1, :order => 1:nbess);

#Main sets of the optimization problem
N = 1:nbuses;
BR = 1:nbranches;
NPER = 1:nper;
TH = 1:ngendis;
RE = 1:ngennodis;
BESS = 1:nbess;

println(string("All input files (csv) were loaded successfully"))

# Comparison of inverter outputs
global pdiss_results = zeros(nper,nbess, ndays);
global pchh_results = zeros(nper,nbess, ndays);

for d in 1:1
    println(string("    "))
    println(string("Day ",d))
    #Measure the total execution time of the program
    #total_timing=@elapsed begin

    #Creat a folder for OPF results of each day
    global daterun=start_case+Dates.Day(d-1);
    global opffolder=string(path_output,"Results ","Day ", d, " Date ", 
                            string(daterun)," ", 
                            string(Dates.dayname(daterun))," ",
                            bess_modeling==1 ? "(Considering accurate BESS model)" : "(Considering traditional BESS model)"
                            );
       
    #Create folder
    mkpath(opffolder);
        
    global pcc=zeros(nbuses,nper);
    global premax=zeros(ngennodis,nper);

    maxre_results=Array{Any}(nothing,nper+3,ngennodis+1);
    maxre_results[1:3,1]=["plant_id" "bus_agg" "type"]; 

    loadp_results=Array{Any}(nothing,nper+2,nbuses+1);
    loadp_results[1:2,1]=["order" "bus_agg"];

        for t in NPER
            maxre_results[t+3,1]=t;
            for g in RE
                maxre_results[1:3,g+1]=[nondis_units[g, :plant_id][1] nondis_units[g, :bus_agg][1] nondis_units[g, :type][1]] ;
                premax[g,t]=prnw(g,d,t);
                maxre_results[t+3,g+1]=premax[g,t];
            end

            loadp_results[t+2,1]=t;
            for n in N
                loadp_results[1:2,n+1]=[bus[n, :order][1] bus[n, :bus_agg][1]];
                pcc[n,t]=pc(n,d,t);
                loadp_results[t+2,n+1]=pcc[n,t];              
             end
        end

        global pathopf = string(opffolder, "/", "bess discharging power(MW).csv");
        global dfopf = CSV.read(pathopf, DataFrame);
        global pdiss_results[:,:,d] = Matrix(dfopf[1:end,2:end]); 

        global pathopf = string(opffolder, "/", "bess charging power(MW).csv");
        global dfopf = CSV.read(pathopf, DataFrame);
        global pchh_results[:,:,d] = Matrix(dfopf[1:end,2:end]); 


 #File for printing graphs

 total_demand=[sum(pcc[n,t] for n in N) for t in NPER];
 min_y=ceil(minimum(total_demand)-50);
 max_y=ceil(maximum(total_demand)+20);
 load = bar(1:1:24, 
             total_demand,
             title="Total demand of the system", 
             xticks=1:1:nper,
             yticks=min_y:20:max_y,
             ylims=(min_y,max_y),
             xlabel="Hour (h)", 
             ylabel="Power (MW)",
             linewidth=2,
             color=:white,
             titlefontsize=11,
             labelfontsize=11,
             xtickfontsize=9,
             ytickfontsize=9,
             legendfontsize=7)
 namegraphic = string(opffolder,"/","Total demand(MW).pdf") 
 savefig(load, namegraphic);

 sample_battery=2;
 sample_bess = bar(1:1:24, 
             [pdiss_results[:,sample_battery,d] pchh_results[:,sample_battery,d] ],
             label=["Discharging power" "Charging power"], 
             title=string("Operation of the battery in bus ", bess[sample_battery, :bus_agg]), 
             xticks=1:1:nper,
             xlabel="Hour (h)", 
             ylabel="Power (MW)",
             legend=:bottomright, 
             linewidth=1,
             left_margin = 0.5Plots.mm, right_margin = 12Plots.mm,
             bottom_margin=3Plots.mm,
             titlefontsize=11,
             labelfontsize=11,
             xtickfontsize=9,
             ytickfontsize=9,
             legendfontsize=7);

 sample_bess = plot!(twinx(),
             ee.data[sample_battery,:],
             color=:green,
             xticks = :false,   
             left_margin = 0.5Plots.mm, right_margin = 14Plots.mm,
             bottom_margin=3Plots.mm,
             label="SOC",
             ylabel="Energy (MWh)",
             grid=true,
             grid_linewidth=1,      
             linewidth=2,
             legend=:topright,
             titlefontsize=11,
             labelfontsize=11,
             xtickfontsize=7,
             ytickfontsize=9,
             legendfontsize=6);

 namegraphic = string(opffolder,"/","Operation of a battery.pdf") 
 savefig(sample_bess, namegraphic);

 #Print the generation by resource

 global  pathopf= string(opffolder, "/", "generation by resource(MW).csv");
 global  genres_results_df = CSV.read(pathopf, DataFrame);
 

 # Rescale hours (to start from zero)
 genres_results_df.hour = genres_results_df.hour .-1;

 gen_res=genres_results_df |>
 @vlplot(:area, 
 x=:hour, y={:generation, stack=:zero}, 
 color={"resource:n", scale={scheme="category10"}})
 save(string(opffolder,"/","generation by resource(MW).pdf"),  gen_res)

 println(string("All plots were printed successfully"))

    end