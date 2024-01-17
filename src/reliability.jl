mutable struct SystemRunData
    objective::Array{Float64, 3}
    ENS::Array{Float64, 3}
    LOL::Array{Float64, 3}
    curtailment::Array{Float64, 3}
end

@enum State B=1 ST=2 LT=3

SystemRunData(runs::Integer, datasize::Integer) = SystemRunData([zeros(runs, 3, datasize) for _ in 1:4]...)

get_objective(val::SystemRunData, run::Integer, i::Integer) = sum(val.objective[run, :, i])

function gather_run_data!(vals::SystemRunData, run::Integer, i::Integer, model::Model, opfs::OPFsystem, 
    Pc::Dict{<:Integer,ExprC}, Pcc::Dict{<:Integer,ExprCC}, atol::Real=1e-6
)
    MOI.get(model, MOI.ResultCount()) < 1 && return

    vals.objective[run, 1, i] = calc_objective(model, opfs)
    vals.objective[run, 2, i] = calc_objective(model, opfs, Pc)
    vals.objective[run, 3, i] = calc_objective(model, opfs, Pcc)
    # @assert sum(vals.objective[run, :, i]) ≈ JuMP.objective_value(model)

    ens = sum_value_property(model, :ls0)
    vals.ENS[run, 1, i] = ens
    vals.LOL[run, 1, i] = abs(ens)>atol ? 1 : 0
    ens = sum_value_property(model, Pc, :lsc)
    vals.ENS[run, 2, i] = sum(ens)
    vals.LOL[run, 2, i] = count(x->abs(x)>atol, ens)
    ens = sum_value_property(model, Pcc, :lscc)
    vals.ENS[run, 3, i] = sum(ens)
    vals.LOL[run, 3, i] = count(x->abs(x)>atol, ens)

    vals.curtailment[run, 1, i] = sum(sum_value_property(model, :pr0))
    vals.curtailment[run, 2, i] = sum(sum_value_property(model, Pc, :prc))
    vals.curtailment[run, 3, i] = sum(sum_value_property(model, Pcc, :prcc))

    return vals
end

function print_data(val::SystemRunData)
    for sn_n in propertynames(val)
        println(sn_n)
        sn = getproperty(val, sn_n)
        for run in 1:size(val.ENS, 1)
            println(string(run))
            for col in axes(sn, 2)
                @printf("%12s ", string(State(col)))
                for row in axes(sn, 3)
                    @printf("%12.4f ", sn[run,col,row])
                end
                println()
            end
        end
    end
end

sum_value_property(model::Model, symb::Symbol) = sum(get_value(model, symb))
sum_value_property(model::Model, P::Dict, symb::Symbol) = [sum(get_value(model, getproperty(c, symb))) for (_,c) in P]

function run_case!(result, i, c, case, goal, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    model, opf, pf, oplim, Pc, Pcc, Pccx = opf_base(case, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec);
    solve_model!(model);
    MOI.get(model, MOI.ResultCount()) < 1 && return
    if case != goal
        fix_base_case!(model)
        case.C1 && fix_contingencies!(model, Pc)
        case.C2 && fix_contingencies!(model, Pcc)
        model, opf, pf, oplim, Pc, Pcc, Pccx = add_all_contingencies!(goal - case, opf, oplim, model, pf, Pc, Pcc, Pccx)
        solve_model!(model);
    end
    gather_run_data!(result, c, i, model, opf, Pc, Pcc)
    return
end

function run_cases!(result, i, cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    for (c, case) in enumerate(cases)
        run_case!(result, i, c, case, cases[end], optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)
    end
    return
end

function run_contingency_select_case!(result, i, c, case, goal, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    model, opf, pf, oplim, Pc, Pcc, Pccx = opf_base(case, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec);
    solve_model!(model);
    MOI.get(model, MOI.ResultCount()) < 1 && return
    if case != goal
        fix_base_case!(model)
        case.C1 && fix_contingencies!(model, Pc)
        case.C2 && fix_contingencies!(model, Pcc)
        solve_model!(model);
        model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_contingency_select!(goal - case, model, opf, pf, oplim, Pc, Pcc, Pccx)
    end
    gather_run_data!(result, c, i, model, opf, Pc, Pcc)
    return
end

function run_contingency_select_cases!(result, i, cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    for (c, case) in enumerate(cases)
        run_contingency_select_case!(result, i, c, case, cases[end], optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)
    end
    return
end

function run_benders_case!(result, i, c, case, goal, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    model, opf, pf, oplim, Pc, Pcc, Pccx = opf_base(case, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec);
    solve_model!(model);
    MOI.get(model, MOI.ResultCount()) < 1 && return
    if case != goal
        fix_base_case!(model)
        case.C1 && fix_contingencies!(model, Pc)
        case.C2 && fix_contingencies!(model, Pcc)
        solve_model!(model);
        model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_benders!(goal - case, model, opf, pf, oplim, Pc, Pcc, Pccx)
    end
    gather_run_data!(result, c, i, model, opf, Pc, Pcc)
    return
end

function run_benders_cases!(result, i, cases, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    for (c, case) in enumerate(cases)
        run_benders_case!(result, i, c, case, cases[end], optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)
    end
    return
end

function run_reliability_calculation_benders(cases, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    # demands = sort_components!(get_demands(system))
    # pd = demands .|> get_active_power
    # hours = read_x_data("data\\ieee_std_load_profile.txt")
    # hours = 0.5:0.1:1.3
    hours = 0.8:0.1:1.1

    result = SystemRunData(4, length(hours))
    Logging.disable_logging(Logging.Info)
    # buffer_sys = [deepcopy(system) for _ in 1:Threads.nthreads()]
    # Threads.@threads for i in eachindex(hours)
    for i in eachindex(hours)
        h = hours[i]
        # sys = buffer_sys[Threads.threadid()]
        set_active_power!.(get_components(StaticLoad, system), (get_components(StaticLoad, system) .|> get_max_active_power) * h)
        # cont = sort_components!(get_branches(sys))

        run_benders_cases!(result, i, cases, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)

        print(i, " ")
    end
    Logging.disable_logging(Logging.Debug)
    # set_active_power!.(get_components(StaticLoad, system), (get_components(StaticLoad, system) .|> get_max_active_power))
    return result
end

read_x_data(fname) = vec(DelimitedFiles.readdlm(fname,Float64))

get_objective(mod::Model) = termination_status(mod) != MOI.OPTIMAL ? -1.0 : objective_value(mod)

get_ens(mod::Model, cont::Dict, symb::Symbol) = 
    sum(sum(get_value(mod, getfield(c, symb))) for (i,c) in cont)
    
get_lol(mod::Model, cont::Dict, symb::Symbol, atol::Real=1e-6) = 
    sum(count(x->(abs(x)>atol), get_value(mod, getfield(c, symb))) for (i,c) in cont)
