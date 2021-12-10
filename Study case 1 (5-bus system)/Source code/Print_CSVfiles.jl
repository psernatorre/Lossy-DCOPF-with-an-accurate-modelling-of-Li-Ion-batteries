        #File for printing csv files

        #Results of thermal power
        pth_results=Array{Any}(nothing,nper+1,ngendis+1);
        pth_results[1,1]=string("periods");
        for g in TH, t in NPER
                pth_results[1,g+1]=string("gen ", dis_units[g, :plant_id]);
                pth_results[t+1,[1 g+1]]=[t pthh[g,t]];
        end
        pth_results_df = DataFrame(pth_results[2:end,:], :auto);
        rename!( pth_results_df, Symbol.(pth_results[1,:]));
        CSV.write(string(opffolder,"/","thermal power generation(MW).csv"), pth_results_df);
        

        #Results of renewable power
        pre_results=Array{Any}(nothing,nper+1,ngennodis+1);
        pre_results[1,1]=string("periods");
        for g in RE, t in NPER
                pre_results[1,g+1]=string("gen ", nondis_units[g, :plant_id]);
                pre_results[t+1,[1 g+1]]=[t pree[g,t]];
        end
        pre_results_df = DataFrame(pre_results[2:end,:], :auto);
        rename!( pre_results_df, Symbol.(pre_results[1,:]));
        CSV.write(string(opffolder,"/","renewable power generation(MW).csv"), pre_results_df);

        #Results of bess discharging power
        dis_results=Array{Any}(nothing,nper+1,nbess+1);
        dis_results[1,1]=string("periods");
        for b in BESS, t in NPER
                 dis_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                 dis_results[t+1,[1 b+1]]=[t pdiss[b,t]];
        end
        dis_results_df = DataFrame(dis_results[2:end,:], :auto);
        rename!( dis_results_df, Symbol.(dis_results[1,:]));
        CSV.write(string(opffolder,"/","bess discharging power(MW).csv"), dis_results_df);

        #Results of bess charging power
        cha_results=Array{Any}(nothing,nper+1,nbess+1);
        cha_results[1,1]=string("periods");
        for b in BESS, t in NPER
                cha_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                cha_results[t+1,[1 b+1]]=[t pchh[b,t]];
        end
        cha_results_df = DataFrame(cha_results[2:end,:], :auto);
        rename!( cha_results_df, Symbol.(cha_results[1,:]));
        CSV.write(string(opffolder,"/","bess charging power(MW).csv"), cha_results_df);

        #Results of bess charging status
        ucc_results=Array{Any}(nothing,nper+1,nbess+1);
        ucc_results[1,1]=string("periods");
        for b in BESS, t in NPER
                ucc_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                ucc_results[t+1,[1 b+1]]=[t (ucc[b,t]>=0.99 ? 1 : 0)];
        end
        ucc_results_df = DataFrame(ucc_results[2:end,:], :auto);
        rename!( ucc_results_df, Symbol.(ucc_results[1,:]));
        CSV.write(string(opffolder,"/","bess charging status.csv"), ucc_results_df);

        #Results of bess discharging status
        udd_results=Array{Any}(nothing,nper+1,nbess+1);
        udd_results[1,1]=string("periods");
        for b in BESS, t in NPER
                udd_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                udd_results[t+1,[1 b+1]]=[t (udd[b,t]>=0.99 ? 1 : 0)]   ;
        end
        udd_results_df = DataFrame(udd_results[2:end,:], :auto);
        rename!( udd_results_df, Symbol.(udd_results[1,:]));
        CSV.write(string(opffolder,"/","bess discharging status.csv"), udd_results_df);
        
        #Results of bess energy
        ee_results=Array{Any}(nothing,nper+1,nbess+1);
        ee_results[1,1]=string("periods");
        for b in BESS, t in NPER
                ee_results[1,b+1]=string("bess ", bess[b, :bess_id]);
                ee_results[t+1,[1 b+1]]=[t ee[b,t]];
        end
        ee_results_df = DataFrame(ee_results[2:end,:], :auto);
        rename!( ee_results_df, Symbol.(ee_results[1,:]));
        CSV.write(string(opffolder,"/","bess energy (MWh).csv"), ee_results_df);

        #Results of buses angle
        theta_results=Array{Any}(nothing,nper+1,nbuses+1);
        theta_results[1,1]=string("periods");
        for n in N, t in NPER
                theta_results[1,n+1]=string("bus ", n);
                theta_results[t+1,[1 n+1]]=[t thetaa[n,t]*180/pi];
        end
        theta_results_df = DataFrame(theta_results[2:end,:], :auto);
        rename!(theta_results_df, Symbol.(theta_results[1,:]));
        CSV.write(string(opffolder,"/","theta(degree).csv"), theta_results_df);

        #Results of LMP
        lmp_results=Array{Any}(nothing,nper+1,nbuses+1);
        lmp_results[1,1]=string("periods");
        for n in N, t in NPER
                lmp_results[1,n+1]=string("bus ", n);
                lmp_results[t+1,[1 n+1]]=[t lmp[n,t]];
        end
        lmp_results_df = DataFrame(lmp_results[2:end,:], :auto);
        rename!(lmp_results_df, Symbol.(lmp_results[1,:]));
        CSV.write(string(opffolder,"/","LMP(Dolar per MWh).csv"), lmp_results_df);

        flow_results=Array{Any}(nothing,nper+1,2*nbranches+1);
        flow_results[1,1]=string("periods");
        for br in BR
            flow_results[1,2*br]=string("id",branch[br,:branch_id]," ",mapbranch[br,1],"->",mapbranch[br,2]);
            flow_results[1,2*br+1]=string("id",branch[br,:branch_id]," ",mapbranch[br,2],"->",mapbranch[br,1]);
            for t in NPER
                flow_results[t+1,1]=t;
                flow_results[t+1,2*br]=flow[br,t]+losses_branch[iter,br,t]/2;
                flow_results[t+1,2*br+1]=-flow[br,t]+losses_branch[iter,br,t]/2;
            end
        end
        flow_results_df = DataFrame(flow_results[2:end,:], :auto);
        rename!(flow_results_df, Symbol.(flow_results[1,:]));
        CSV.write(string(opffolder,"/","power flow branches(MW).csv"), flow_results_df);

        flowdc_results=Array{Any}(nothing,nper+1,2*nbranches+1);
        flowdc_results[1,1]=string("periods");
        for br in BR
            flowdc_results[1,2*br]=string("id",branch[br,:branch_id]," ",mapbranch[br,1],"->",mapbranch[br,2]);
            flowdc_results[1,2*br+1]=string("id",branch[br,:branch_id]," ",mapbranch[br,2],"->",mapbranch[br,1]);
            for t in NPER
                flowdc_results[t+1,1]=t;
                flowdc_results[t+1,2*br]=flow[br,t];
                flowdc_results[t+1,2*br+1]=-flow[br,t];
            end
        end
        flowdc_results_df = DataFrame(flowdc_results[2:end,:], :auto);
        rename!(flowdc_results_df, Symbol.(flowdc_results[1,:]));
        CSV.write(string(opffolder,"/","dc power flow branches(MW).csv"), flowdc_results_df);

        losses_results=Array{Any}(nothing,nper+1,nbranches+1);
        losses_results[1,1]=string("periods");
        for br in BR, t in NPER
            losses_results[1,br+1]=string("id",branch[br,:branch_id]," ",mapbranch[br,1],"-",mapbranch[br,2]);
            losses_results[t+1,1]=t;
            losses_results[t+1,br+1]=losses_branch[iter,br,t];
        end
        losses_results_df = DataFrame(losses_results[2:end,:], :auto);
        rename!(losses_results_df, Symbol.(losses_results[1,:]));
        CSV.write(string(opffolder,"/","losses(MW).csv"), losses_results_df);

        #Results of p shedding 
        pshed_results=Array{Any}(nothing,nper+1,nbuses+1);
        pshed_results[1,1]=string("periods");
        for n in N, t in NPER
            pshed_results[1,n+1]=string("bus ",n);
            pshed_results[t+1,[1 n+1]]=[t pshedd[n,t]];
        end
        pshed_results_df = DataFrame(pshed_results[2:end,:], :auto);
        rename!(pshed_results_df, Symbol.(pshed_results[1,:]));
        CSV.write(string(opffolder,"/","demand response(MW).csv"), pshed_results_df);

        #Results of curtailment
        curt_results=Array{Any}(nothing,nper+1,ngennodis+1);
        curt_results[1,1]=string("periods");
        for g in RE, t in NPER
                 curt_results[1,g+1]=string("gen ", nondis_units[g, :plant_id]);
                 curt_results[t+1,[1 g+1]]=[t premax[g,t]-pree[g,t]];
        end
        curt_results_df = DataFrame(curt_results[2:end,:], :auto);
        rename!( curt_results_df, Symbol.(curt_results[1,:]));
        CSV.write(string(opffolder,"/","curtailment(MW).csv"), curt_results_df );
         
         #Results of Energy balance
        activbal_results=Array{Any}(nothing,nper+1,9);
        activbal_results[1,1]="Periods";
        activbal_results[1,2]="Thermal generation(MW)";
        activbal_results[1,3]="Renewable generation(MW)";
        activbal_results[1,4]="BESS discharging(MW)";
        activbal_results[1,5]="Demand response(MW)";
        activbal_results[1,6]="Load(MW)";
        activbal_results[1,7]="BESS charging(MW)";
        activbal_results[1,8]="Losses(MW)";
        activbal_results[1,9]="Error";
        for t in 1:nper 
            activbal_results[t+1,1]=t;
            activbal_results[t+1,2]=sum(pthh[g,t] for g in TH);
            activbal_results[t+1,3]=sum(pree[g,t] for g in RE);
            activbal_results[t+1,4]=sum(pdiss[b,t] for b in BESS);
            activbal_results[t+1,5]=sum(pshedd[n,t] for n in N);
            activbal_results[t+1,6]=sum(pcc[n,t] for n in N);
            activbal_results[t+1,7]=sum(pchh[b,t] for b in BESS);
            activbal_results[t+1,8]=sum(losses_branch[iter,br,t] for br in BR);
            activbal_results[t+1,9]=activbal_results[t+1,2]+activbal_results[t+1,3]+activbal_results[t+1,4]+activbal_results[t+1,5]+
                                    -activbal_results[t+1,6]-activbal_results[t+1,7]-activbal_results[t+1,8];
        end
        activbal_results_df = DataFrame(activbal_results[2:end,:], :auto);
        rename!(activbal_results_df, Symbol.(activbal_results[1,:]));
        CSV.write(string(opffolder,"/","active power balance(MW).csv"), activbal_results_df);

        #Results of Objective Function

        OF_results=Array{Any}(nothing,9,2);
        OF_results[1,1]="Parameter";
        OF_results[1,2]="Value";
        OF_results[2,1]="Power shedding cost (USD)";
        OF_results[2,2]=price_p_shed*deltaT*sum(pshedd[n,t] for n in N, t in NPER);
        OF_results[3,1]="Thermal variable cost (USD)";
        OF_results[3,2]=deltaT*sum(gencost(g)*pthh[g,t] for g in TH,t in NPER);   
        OF_results[4,1]="Discharging power cost (USD)";
        OF_results[4,2]=deltaT*sum(bess[b, :DischCost]*pdiss[b,t] for b in BESS,t in NPER);
        OF_results[5,1]="Charging power cost (USD)";
        OF_results[5,2]=(-1)*deltaT*sum(bess[b, :ChargCost]*pchh[b,t] for b in BESS,t in NPER);
        OF_results[6,1]="Total cost (USD)";
        OF_results[6,2]=OF_results[2,2]+OF_results[3,2]+OF_results[4,2]+OF_results[5,2];
        OF_results[7,1]="Status";
        OF_results[7,2]=string(status1," ",status2);
        OF_results[8,1]="Type of BESS modeling";
        OF_results[8,2]=bess_modeling==1 ? "Accurate" : "Traditional";
        OF_results[9,1]="Capacity factors and load data from";
        OF_results[9,2]=string(daterun);
        OF_results_df = DataFrame(OF_results[2:end,:], :auto);
        rename!(OF_results_df, Symbol.(OF_results[1,:]));
        CSV.write(string(opffolder,"/","objective function.csv"), OF_results_df);

        resource=unique(plant[:,:type]);
        resource=vcat(resource,"_discharge_bess");
        resource=vcat(resource,"_charge_bess");
        resource=vcat(resource, "curtailment")
        nresource=size(resource,1);

        genres_results=Array{Any}(nothing,nper*nresource+1,3);
        genres_results[1,1]=string("resource");
        genres_results[1,2]=string("hour");
        genres_results[1,3]=string("generation");
        for t in NPER
                for r in 1:nresource
                        genres_results[(t-1)*nresource+r+1,1]=resource[r];
                        genres_results[(t-1)*nresource+r+1,2]=t;
                        genres_results[(t-1)*nresource+r+1,3]=sum(dis_units[g,:type]==resource[r] ? pthh[g,t] : 0 for g in TH)+
                                                              sum(nondis_units[g,:type]==resource[r] ? pree[g,t] : 0 for g in RE)+
                                                              sum(resource[r]=="_discharge_bess" ? pdiss[b,t] : 0 for b in BESS)+
                                                              sum(resource[r]=="_charge_bess" ? pchh[b,t] : 0 for b in BESS)   +
                                                              sum(resource[r]=="curtailment" ? (premax[g,t]-pree[g,t]) : 0 for g in RE);
                end
        end
        genres_results_df = DataFrame(genres_results[2:end,:], :auto);
        rename!(genres_results_df, Symbol.(genres_results[1,:]));
        CSV.write(string(opffolder,"/","generation by resource(MW).csv"), genres_results_df);

        genres_results_df[!,:resource] = convert.(String, genres_results_df[!,:resource]);
        genres_results_df[!,:hour] = convert.(Int, genres_results_df[!,:hour]);
        genres_results_df[!,:generation] = convert.(Float64,  genres_results_df[!,:generation]);        

        #Results OF in details
        OFd_results=Array{Any}(nothing,nper+1,6);
        OFd_results[1,1]=string("Periods");
        OFd_results[1,2]=string("Thermal variable cost");
        OFd_results[1,3]=string("Discharging power cost");
        OFd_results[1,4]=string("Charging power cost");
        OFd_results[1,5]=string("Power shedding cost");
        OFd_results[1,6]=string("Total Cost");

        for t in NPER
                OFd_results[t+1,1]=t;
                OFd_results[t+1,2]=deltaT*sum(gencost(g)*pthh[g,t] for g in TH);
                OFd_results[t+1,3]=deltaT*sum(bess[b, :DischCost]*pdiss[b,t] for b in BESS);
                OFd_results[t+1,4]=(-1)*deltaT*sum(bess[b, :ChargCost]*pchh[b,t] for b in BESS);
                OFd_results[t+1,5]=price_p_shed*deltaT*sum(pshedd[n,t] for n in N if bus[n, :Pd]!=0);
                OFd_results[t+1,6]=OFd_results[t+1,2]+OFd_results[t+1,3]+OFd_results[t+1,4]+OFd_results[t+1,5];   
        end
    
        OFd_results_df = DataFrame(OFd_results[2:end,:], :auto);
        rename!(OFd_results_df, Symbol.(OFd_results[1,:]));
        CSV.write(string(opffolder,"/","hourly objective function.csv"), OFd_results_df);
  
        println(string("All CSV files were printed successfully"))

      