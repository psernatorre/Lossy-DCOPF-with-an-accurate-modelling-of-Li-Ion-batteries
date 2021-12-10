#Packages
using DataFrames, CSV, CPLEX, JuMP, LinearAlgebra, Dates, Plots, VegaLite

using Plots.PlotMeasures

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

path_output_comparison = string(folder_name, "/",output_folder, "/", "Comparison/");

mkpath(path_output_comparison);

# Enter the path of the traditional battery model
path_output_accurate =string(folder_name, "/",output_folder,"/", "accurate");

# Enter the path of the accurate battery model
path_output_traditional =string(folder_name, "/",output_folder,"/", "traditional");

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
plant=CSV.read(string(path_inputdata,"plant.csv"), DataFrame)[:,[1,2,10,11,27,29,30]];
ngen=size(plant,1);


#Non-dispatchable units
dis_units=plant[ (plant.type .!= "wind") .& (plant.type .!= "solar") , :];
ngendis=size(dis_units,1);
insertcols!(dis_units, 1, :order => 1:ngendis);
gencost(g)=dis_units[g, :GenFuelCost]*dis_units[g, :GenIOB]; 

nondis_units=plant[ (plant.type .== "wind") .| (plant.type .== "solar") , :];
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
branch=CSV.read(string(path_inputdata,"branch.csv"), DataFrame)[:,1:7];
nbranches=size(branch,1);
insertcols!(branch, 1, :order => 1:nbranches);
mapbranch=zeros(nbranches,nbranches);
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

# Plot Objective Function
objtrad_results = zeros(nper, ndays);
objtraddis_results = zeros(nper, ndays);
objtradch_results = zeros(nper, ndays);
objacc_results = zeros(nper, ndays);
objaccdis_results = zeros(nper, ndays);
objaccch_results = zeros(nper, ndays);

for d in 1:ndays
    global daterun=start_case+Dates.Day(d-1);
    global opffolder=string(path_output,"Results ","Day ", d, " Date ", 
                            string(daterun)," ", 
                            string(Dates.dayname(daterun))," ",
                            "(Considering traditional BESS model)");
    global pathopf = string(opffolder,"/", "hourly objective function.csv");
    global objtrad = CSV.read(pathopf, DataFrame);
    global objtrad_results[:,d] = Matrix(objtrad[1:end,2:end])[:,5]./1000;
    global objtraddis_results[:,d] = Matrix(objtrad[1:end,2:end])[:,2]./1000 ;  
    global objtradch_results[:,d] = Matrix(objtrad[1:end,2:end])[:,3]./1000 ;  
end

for d in 1:ndays
    global daterun=start_case+Dates.Day(d-1);
    global opffolder=string(path_output,"Results ","Day ", d, " Date ", 
                            string(daterun)," ", 
                            string(Dates.dayname(daterun))," ",
                            "(Considering accurate BESS model)");
    global pathopf = string(opffolder,"/", "hourly objective function.csv");
    global objacc = CSV.read(pathopf, DataFrame);
    global objacc_results[:,d] = Matrix(objacc[1:end,2:end])[:,5]./1000 ; 
    global objaccdis_results[:,d] = Matrix(objacc[1:end,2:end])[:,2]./1000 ;  
    global objaccch_results[:,d] = Matrix(objacc[1:end,2:end])[:,3]./1000 ;  
end



for d in 1:ndays
#multiply by 1000 if the system is 5 bus.
min=round(minimum(objtraddis_results[:,d]+objtradch_results[:,d]),digits=0)-3;
max=round(maximum(objtraddis_results[:,d]+objtradch_results[:,d]),digits=0)+3;

            fig = plot(range(1, nper, length=nper), 
                [(objtraddis_results[:,d]+objtradch_results[:,d]) (objaccdis_results[:,d]+objaccch_results[:,d])], 
                label=["Net BESS cost (traditional)" "Net BESS cost (accurate)"],
                left_margin = 0.3Plots.mm, right_margin = 13Plots.mm,
                xticks=1:1:nper;
                xlabel="Hour (h)", 
                yticks=min:5:max,
                ylim=(min,max),
                ylabel="Net cost of BESS operation (kUSD)",
                color=[:black :green],
                linewidth=2,
                legend=:topleft,
                titlefontsize=11,
                labelfontsize=11,
                xtickfontsize=9,
                ytickfontsize=9,
                legendfontsize=7)
      
delta=(objtraddis_results[:,d]+objtradch_results[:,d])-(objaccdis_results[:,d]+objaccch_results[:,d])
        fig = plot!(twinx(),
                [objtrad_results[:,d] objacc_results[:,d]],
                label=["O.F. (traditional)" "O.F. (accurate)"],
                ylabel="Total objective function (kUSD)",
                xticks = :false,   
                left_margin = 0.2Plots.mm, right_margin = 0.2Plots.mm,
                bottom_margin=3Plots.mm,
                grid=true,
                grid_linewidth=1,      
                linewidth=2,
                legend=:bottomright,
                titlefontsize=11,
                labelfontsize=11,
                xtickfontsize=9,
                ytickfontsize=9,
                legendfontsize=7);
                global daterun=start_case+Dates.Day(d-1);
                global namegraphic = string(path_output_comparison, "Objective Function (kUSD)"," Day ", d," ",string(daterun),".pdf");     
                savefig(fig, namegraphic);
end