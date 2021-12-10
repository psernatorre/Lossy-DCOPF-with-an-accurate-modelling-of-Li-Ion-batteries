#Project: A lossy DC Optimal Power Flow with an accuratemodeling of Li-Ion batteries
#Student: P. Serna Torre
#Couse: MAE207 / Professor: M. Davidson

#Warning: Paths used here have Linux format (i.e. /..../), which is different to Windows format ("\\.....\\")
println(string("  "))
println(string("      .................. #  # ...................      "))
println(string("      Lossy DCOPF with an accurate BESS modelling      "))
println(string("      Student: Paul Serna-Torre"))
println(string("      Course: MAE207 UCSD, Professor: M. Davidson"))
println(string("      Language: Julia 1.6.3 - Ubuntu 20.04"))
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
println(string(" "))
println(string("General characteristics of the network under analysis: "))
println(string("Buses: ", nbuses))
println(string("Transmission lines: ", nbranches))
println(string("Resolution (periods per day): ", nper))
println(string("Dispatchable generators: ", ngendis))
println(string("Non-dispatchable generators: ", ngennodis))
println(string("Batteries (BESS): ", nbess))

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


        maxre_results_df = DataFrame(maxre_results[2:end,:], :auto);
        rename!(maxre_results_df, Symbol.(maxre_results[1,:]));
        CSV.write(string(opffolder,"/","renewable energy availability(MW).csv"), maxre_results_df);
 
        loadp_results_df = DataFrame(loadp_results[2:end,:], :auto);
        rename!(loadp_results_df, Symbol.(loadp_results[1,:]));
        CSV.write(string(opffolder,"/","load(MW).csv"), loadp_results_df);

        println(string("The files renewable generation availability.csv and load.csv were printed successfully"))

        #Power Flow Model
        global flow=zeros(nbranches,nper);
        global pthh=zeros(ngendis,nper);
        global pree=zeros(ngennodis,nper);
        global pshedd=zeros(nbuses,nper);
        global thetaa=zeros(nbuses,nper);
        global pdiss=zeros(nbess,nper);
        global pchh=zeros(nbess,nper);
        global udd=zeros(nbess,nper);
        global ucc=zeros(nbess,nper);
        global ee=zeros(nbess,nper);
        global lmp=zeros(nbuses,nper);
    
        println(string("Power Flow"))
        println(string("Day ", d, " Running..."))
        iter=1;
        error_losses=ones(iter_losses_max+1);
        error_losses[1]=1000;
        global losses_branch=zeros(iter_losses_max+1,nbranches,nper);
        global pthh_iter=ones(iter_losses_max+1,ngendis,nper);

        while error_losses[iter] >= tol_losses && iter<=iter_losses_max;
            
            println(" ")
            println("Iteration ", iter, " Refining Losses")
            
            PowerFlow = JuMP.Model(CPLEX.Optimizer);

            set_optimizer_attribute(PowerFlow, "CPX_PARAM_BARQCPEPCOMP",tolerance_solver)
            set_optimizer_attribute(PowerFlow, "CPX_PARAM_THREADS", cores_solver)
            set_silent(PowerFlow)

            #Decision variables
            @variable(PowerFlow, Cost);
            @variable(PowerFlow, f[BR, NPER]);
            @variable(PowerFlow, 0 <= pth[TH, NPER]);
            @variable(PowerFlow, 0 <= pre[RE, NPER]);
            @variable(PowerFlow, 0 <= pshed[N, NPER]);
            @variable(PowerFlow, -pi/2 <= theta[N, NPER]  <= pi/2);
            @variable(PowerFlow, 0 <= e[BESS, NPER]);
            @variable(PowerFlow, 0 <= pch[BESS, NPER]);
            @variable(PowerFlow, 0 <= pdis[BESS, NPER]);

            if bess_modeling==1
                @variable(PowerFlow, ud[BESS, NPER], Bin);
                @variable(PowerFlow, uc[BESS, NPER], Bin);
            end
           
            #Objective Function
            @objective(PowerFlow, Min, Cost);
        
            #Active power balance constraint (losses are considered as fictitious loads)
            @constraint(PowerFlow, pbalance[n in N, t in NPER], 
                        deltaT*(sum( dis_units[g, :bus_agg]==bus[n,:bus_agg] ? pth[g,t] : 0 for g in TH) 
                        + sum( nondis_units[g, :bus_agg]==bus[n,:bus_agg] ? pre[g,t] : 0 for g in RE) 
                        + sum( bess[b, :bus_agg]==bus[n,:bus_agg] ? (pdis[b,t] - pch[b,t]) : 0 for b in BESS) 
                        + (bus[n, :Pd]!=0 ? pshed[n,t] : 0) - pcc[n,t])
                        == deltaT*sum( (exist(n,mapbranch[br,2],mapbranch[[br],:])==1) ? (f[br,t]+losses_branch[iter,br,t]/2) : 0 for br in BR )  
                           -deltaT*sum( (exist(mapbranch[br,1],n,mapbranch[[br],:])==1) ? (f[br,t]-losses_branch[iter,br,t]/2) : 0 for br in BR )                   
                        );
                     #   sum( (exist(n,mapbranch[br,2],mapbranch[[br],:])==1) ? 1 : 0 for br in BR ) 
                     #   deltaT*sum( (exist(mapbranch[br,1],n,mapbranch[[br],:])==1) ? 1 : 0 for br in BR)                  
            #DC power flow equation
            @constraint(PowerFlow, pflow[br in BR, t in NPER],
                        f[br,t] ==  MVAbase*x(br)/(x(br)^2+r(br)^2)*(theta[mapbranch[br,1],t] -  theta[mapbranch[br,2],t]) );
            
            #Capacity of the transmission line (in one direction)     
            @constraint(PowerFlow, flowcap1[br in BR, t in NPER],
                        f[br,t]+losses_branch[iter,br,t]/2<= branch[br,:rateA]);
    
            #Capacity of the transmission line (in the opposite direction)        
            @constraint(PowerFlow, flowcap2[br in BR, t in NPER],
                        f[br,t]-losses_branch[iter,br,t]/2>= -branch[br,:rateA]);
    
            #Maximum bound of power generation
            @constraint(PowerFlow, maxpowerth[g in TH, t in NPER], 
                        pth[g,t] <=  dis_units[g, :Pmax]);
    
            #Minimum bound of power generation
            @constraint(PowerFlow, minpowerth[g in TH, t in NPER], 
                        pth[g,t] >=  dis_units[g, :Pmin]);
    
            #Maximum renewable power generation
            @constraint(PowerFlow, maxpowere[g in RE, t in NPER], 
                        pre[g,t] <=   premax[g,t]);
    
            #Maximum power shedding
            @constraint(PowerFlow, maxshed[n=1:nbuses, t=1:nper], 
                        pshed[n,t] <=  (bus[n, :Pd]!=0 ? pcc[n,t] : 0) );
    
            #Fix the reference bus
            @constraint(PowerFlow, thref[t in NPER], 
                        theta[1,t] == 0);
    
            #Bess state equation and initial state
            @constraint(PowerFlow, BessState[b in BESS, t in NPER],
                        e[b,t]== (t==1 ? (bess[b, :Initial_state]*bess[b, :Emax]) : e[b,t-1])
                                 +(pch[b,t]*bess[b, :ChargEff] - pdis[b,t]/bess[b, :DischEff])*deltaT);

            #Bess final state
            @constraint(PowerFlow, BessFin[b in BESS], 
                        e[b,nper]==bess[b, :Final_state]*bess[b, :Emax]);
    
            #If the user enter 1 in the row "BESS modeling" of configuration.csv, we add the following constraints (accurate Bess modeling)
            if bess_modeling==1
         
            #Maximum energy capacity of a battery (considering the relation between the max. energy capacity and the charging power)
            @constraint(PowerFlow, BessEmax[b in BESS, t in NPER],
                        e[b,t] <= bess[b, :Emax]*(1-pch[b, t]/bess[b, :Pmax]*(1-bess[b, :a2_charg])));
            
            #Minimum energy capacity of a battery (considering the relation between the min. energy capacity and the discharging power)
            @constraint(PowerFlow, BessEmin[b in BESS, t in NPER],
                        e[b,t] >= bess[b, :Emax]*bess[b, :a1_charg]/bess[b, :Pmax]*pdis[b, t]);
    
            #Complementary (charging and discharging)
            @constraint(PowerFlow, BessComplementary[b in BESS, t in NPER], 
                        uc[b,t]+ud[b,t] <= 1);
            
            #Maximum bound bess discharging power
            @constraint(PowerFlow, BessDmax[b in BESS, t in NPER],
                        pdis[b,t] <= ud[b,t]*bess[b, :Pmax]); 

            #Minimum bound bess discharging power
            @constraint(PowerFlow, BessDmin[b in BESS, t in NPER],
                        pdis[b,t] >= ud[b,t]*bess[b, :PercMin]*bess[b, :Pmax]); 

            #Maximum bound bess charging power
            @constraint(PowerFlow, BessCmax[b in BESS, t in NPER],
                        pch[b,t] <= uc[b,t]*bess[b, :Pmax]);   

            #Minimum bound bess charging power
            @constraint(PowerFlow, BessCmin[b in BESS, t in NPER],
                        pch[b,t] >= uc[b,t]*bess[b, :PercMin]*bess[b, :Pmax]); 
            
            else #Otherwise, we add the traditional constraints of a BESS modeling

            @constraint(PowerFlow, BessEmax[b in BESS, t in NPER],
                    e[b,t] <= bess[b, :Emax]);
            
            @constraint(PowerFlow, BessDmax[b in BESS, t in NPER],
                    pdis[b,t] <= bess[b, :Pmax]); 
   
            @constraint(PowerFlow, BessCmax[b in BESS, t in NPER],
                    pch[b,t] <= bess[b, :Pmax]);   

            end
        
            #Objective function
            @constraint(PowerFlow, z, 
                        Cost==
                        price_p_shed*deltaT*sum(pshed[n,t] for n in N, t in NPER if bus[n, :Pd]!=0) +
                        deltaT*sum(gencost(g)*pth[g,t] for g in TH,t in NPER) +
                        deltaT*sum(bess[b, :DischCost]*pdis[b,t] for b in BESS,t in NPER) +
                        (-1)*deltaT*sum(bess[b, :ChargCost]*pch[b,t] for b in BESS,t in NPER))
                        ;                              
    
            JuMP.optimize!(PowerFlow);
            global status1=string(primal_status(PowerFlow));
            global status2=string(termination_status(PowerFlow));
    
            println(string("Day ", d," - ", primal_status(PowerFlow)) ,"-", termination_status(PowerFlow));
    
            if status1=="FEASIBLE_POINT"  
                   global  flow=round.(JuMP.value.(f),digits=digits_round);
                   global  pthh=round.(JuMP.value.(pth),digits=digits_round);
                   global  pthh_iter[iter,:,:]=round.(JuMP.value.(pth),digits=digits_round);
                   global  pree=round.(JuMP.value.(pre),digits=digits_round);
                   global  pshedd=round.(JuMP.value.(pshed),digits=digits_round);
                   global  thetaa=round.(JuMP.value.(theta),digits=digits_round);
                   global  pdiss=round.(JuMP.value.(pdis),digits=digits_round);
                   global  pchh=round.(JuMP.value.(pch), digits=digits_round);
                   global  ee=round.(JuMP.value.(e), digits=digits_round);

            #If the modeling of BESS is accurate, we save the solution of the binary variables
                    if bess_modeling==1
                        global  udd=round.(JuMP.value.(ud), digits=digits_round);
                        global  ucc=round.(JuMP.value.(uc), digits=digits_round);
                    end
            end

            #Calculate the losses by using a non-linear expression
           # global losses_branch[iter+1,:,:]=[MVAbase*2*r(br)/(r(br)^2+x(br)^2)*(1-cos(thetaa[mapbranch[br,1],t]-thetaa[mapbranch[br,2],t])) for br in BR, t in NPER];
            global losses_branch[iter+1,:,:]=[1/MVAbase*r(br)*flow[br,t]*flow[br,t] for br in BR, t in NPER];

            #Calculate error of losses between the present iteration and the following.
            if iter>1
             global delta=pthh_iter[iter,:,:] - pthh_iter[iter-1,:,:];
             global error_losses[iter+1]=abs(findmax(delta)[1])/pthh_iter[iter,findmax(delta)[2][1],findmax(delta)[2][2]];
            end
            #global error_losses[iter+1]=abs(sum(pthh_iter[iter,g,t] - pthh_iter[iter-1,g,t] for g in TH, t in NPER))/sum(pthh_iter[iter,g,t] for g in TH,t in NPER);
            #global error_losses[iter+1]=abs(sum(abs(losses_branch[iter+1,br,t] - losses_branch[iter,br,t]) for br in BR, t in NPER))/sum(losses_branch[iter+1,br,t] for br in BR,t in NPER);
            println("Iteration gap (refining losses): ", round(error_losses[iter+1]*100, digits=3),"%")
            global iter=iter+1;
        end     

        #Since the last iteration was i+1, we reset the iteration
        global iter=iter-1;

        #The next model is to calculate LMPs (We fix the binary variables, and copy the same constraints of the previous model)
        println(string(" "))
        println(string("LMP calculation"))

        PowerFlow_LMP = JuMP.Model(CPLEX.Optimizer);

        set_optimizer_attribute(PowerFlow_LMP, "CPX_PARAM_BARQCPEPCOMP",tolerance_solver)
        set_optimizer_attribute(PowerFlow_LMP, "CPX_PARAM_THREADS", cores_solver)
        set_silent(PowerFlow_LMP)

        #Decision variables
        @variable(PowerFlow_LMP, Cost);
        @variable(PowerFlow_LMP, f[BR, NPER]);
        @variable(PowerFlow_LMP, 0 <= pth[TH, NPER]);
        @variable(PowerFlow_LMP, 0 <= pre[RE, NPER]);
        @variable(PowerFlow_LMP, 0 <= pshed[N, NPER]);
        @variable(PowerFlow_LMP, -pi/2 <= theta[N, NPER]  <= pi/2);
        @variable(PowerFlow_LMP, 0 <= e[BESS, NPER]);
        @variable(PowerFlow_LMP, 0 <= pch[BESS, NPER]);
        @variable(PowerFlow_LMP, 0 <= pdis[BESS, NPER]);

        #Objective Function
        @objective(PowerFlow_LMP, Min, Cost);
    
        #Active power balance constraint
        @constraint(PowerFlow_LMP, pbalance[n in N, t in NPER], 
        deltaT*(sum( dis_units[g, :bus_agg]==bus[n,:bus_agg] ? pth[g,t] : 0 for g in TH) 
        + sum( nondis_units[g, :bus_agg]==bus[n,:bus_agg] ? pre[g,t] : 0 for g in RE) 
        + sum( bess[b, :bus_agg]==bus[n,:bus_agg] ? (pdis[b,t] - pch[b,t]) : 0 for b in BESS) 
        + (bus[n, :Pd]!=0 ? pshed[n,t] : 0) - pcc[n,t])
        == deltaT*sum( (exist(n,mapbranch[br,2],mapbranch[[br],:])==1) ? (f[br,t]+losses_branch[iter,br,t]/2) : 0 for br in BR )  
            -deltaT*sum( (exist(mapbranch[br,1],n,mapbranch[[br],:])==1) ? (f[br,t]-losses_branch[iter,br,t]/2) : 0 for br in BR )                
        );
        
        #DC power flow equation
        @constraint(PowerFlow_LMP, pflow[br in BR, t in NPER],
        f[br,t] ==  MVAbase*x(br)/(x(br)^2+r(br)^2)*(theta[mapbranch[br,1],t] -  theta[mapbranch[br,2],t]) );

        #Capacity of the transmission line (in one direction)     
        @constraint(PowerFlow_LMP, flowcap1[br in BR, t in NPER],
        f[br,t]+losses_branch[iter,br,t]/2<= branch[br,:rateA]);

        #Capacity of the transmission line (in the opposite direction)        
        @constraint(PowerFlow_LMP, flowcap2[br in BR, t in NPER],
        f[br,t]-losses_branch[iter,br,t]/2>= -branch[br,:rateA]);


        @constraint(PowerFlow_LMP, maxpowerth[g in TH, t in NPER],
                    pth[g,t] <=   dis_units[g, :Pmax]);

        @constraint(PowerFlow_LMP, minpowerth[g in TH, t in NPER],
                    pth[g,t] >=   dis_units[g, :Pmin]);

        @constraint(PowerFlow_LMP, maxpowere[g in RE, t in NPER],
                    pre[g,t] <=   premax[g,t]);

        @constraint(PowerFlow_LMP, maxshed[n=1:nbuses, t=1:nper],
                    pshed[n,t] <=  (bus[n, :Pd]!=0 ? pcc[n,t] : 0) );

        @constraint(PowerFlow_LMP, thref[t in NPER],
                    theta[1,t] == 0);

        @constraint(PowerFlow_LMP, BessState[b in BESS, t in NPER],
                    e[b,t]== (t==1 ? (bess[b, :Initial_state]*bess[b, :Emax]) : e[b,t-1])
                             +(pch[b,t]*bess[b, :ChargEff] - pdis[b,t]/bess[b, :DischEff])*deltaT);
     
        @constraint(PowerFlow_LMP, BessFin[b in BESS],
                    e[b,nper]==bess[b, :Final_state]*bess[b, :Emax]);
       
        if bess_modeling==1
         
        @constraint(PowerFlow_LMP, BessEmax[b in BESS, t in NPER],
                    e[b,t] <= bess[b, :Emax]*(1-pch[b, t]/bess[b, :Pmax]*(1-bess[b, :a2_charg])));
            
        @constraint(PowerFlow_LMP, BessEmin[b in BESS, t in NPER],
                    e[b,t] >= bess[b, :Emax]*bess[b, :a1_charg]/bess[b, :Pmax]*pdis[b, t]);
                                      
        @constraint(PowerFlow_LMP, BessDmax[b in BESS, t in NPER],
                    pdis[b,t] <= udd[b,t]*bess[b, :Pmax]);  #Here udd is not a decision variable, but a input data that comes from the previous model "Power Flow"
            
        @constraint(PowerFlow_LMP, BessDmin[b in BESS, t in NPER],
                    pdis[b,t] >= udd[b,t]*bess[b, :PercMin]*bess[b, :Pmax]); #Here udd is not a decision variable, but a input data that comes from the previous model "Power Flow"
                        
        @constraint(PowerFlow_LMP, BessCmax[b in BESS, t in NPER],
                    pch[b,t] <= ucc[b,t]*bess[b, :Pmax]);   #Here ucc is not a decision variable, but a input data that comes from the previous model "Power Flow"

        @constraint(PowerFlow_LMP, BessCmin[b in BESS, t in NPER],
                    pch[b,t] >= ucc[b,t]*bess[b, :PercMin]*bess[b, :Pmax]);  #Here ucc is not a decision variable, but a input data that comes from the previous model "Power Flow"
                        
        else
            
        @constraint(PowerFlow_LMP, BessEmax[b in BESS, t in NPER],
                    e[b,t] <= bess[b, :Emax]);
               
        @constraint(PowerFlow_LMP, BessDmax[b in BESS, t in NPER],
                    pdis[b,t] <= bess[b, :Pmax]); 
               
        @constraint(PowerFlow_LMP, BessCmax[b in BESS, t in NPER],
                    pch[b,t] <= bess[b, :Pmax]);   
           
        end

        @constraint(PowerFlow_LMP, z, 
                    Cost==
                    price_p_shed*deltaT*sum(pshed[n,t] for n in N, t in NPER if bus[n, :Pd]!=0) +
                    deltaT*sum(gencost(g)*pth[g,t] for g in TH,t in NPER) +
                    deltaT*sum(bess[b, :DischCost]*pdis[b,t] for b in BESS,t in NPER) +
                    (-1)*deltaT*sum(bess[b, :ChargCost]*pch[b,t] for b in BESS,t in NPER))
                    ;                              

        JuMP.optimize!(PowerFlow_LMP);
        global status1=string(primal_status(PowerFlow_LMP));
        global status2=string(termination_status(PowerFlow_LMP));

        println(string("Day ", d," - ", primal_status(PowerFlow_LMP)) ,"-", termination_status(PowerFlow_LMP));

        if status1=="FEASIBLE_POINT"  
               global  flow=round.(JuMP.value.(f),digits=digits_round);
               global  pthh=round.(JuMP.value.(pth),digits=digits_round);
               global  pree=round.(JuMP.value.(pre),digits=digits_round);
               global  pshedd=round.(JuMP.value.(pshed),digits=digits_round);
               global  thetaa=round.(JuMP.value.(theta),digits=digits_round);
               global  pdiss=round.(JuMP.value.(pdis),digits=digits_round);
               global  pchh=round.(JuMP.value.(pch), digits=digits_round);
               global  ee=round.(JuMP.value.(e), digits=digits_round);
               global  lmp=round.(JuMP.dual.(pbalance), digits=digits_round);

        end

        #Call the file Print_CSVfiles.jl and Print_graphs.jl
        println(string(" "))
        println(string("Creating output files . . ."))

        include(string(main_path,"/Source code/","Print_CSVfiles.jl"))     

        include(string(main_path,"/Source code/","Print_graphs.jl"))     

    end
    
    println(string(" "))
    println(string("All the process was completed successfully"))