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
                    [pdiss.data[sample_battery,:] pchh.data[sample_battery,:] ],
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

        # Rescale hours (to start from zero)
        genres_results_df.hour = genres_results_df.hour .-1;

        gen_res=genres_results_df |>
        @vlplot(:area, 
        x=:hour, y={:generation, stack=:zero}, 
        color={"resource:n", scale={scheme="category10"}})
        save(string(opffolder,"/","generation by resource(MW).pdf"),  gen_res)
      
        println(string("All plots were printed successfully"))
       