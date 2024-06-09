mutable struct SystemRunData
    objective::Array{Float64, 3}
    ENS::Array{Float64, 3}
    LOL::Array{Float64, 3}
    curtailment::Array{Float64, 3}
end

@enum State B=1 ST=2 LT=3

SystemRunData(runs::Integer, datasize::Integer) = SystemRunData([zeros(runs, 3, datasize) for _ in 1:4]...)

get_objective(val::SystemRunData, run::Integer, i::Integer) = sum(val.objective[run, :, i])

function gather_run_data!(vals::SystemRunData, run::Integer, i::Integer, case::Case, atol::Real=1e-14
)
    model = case.model
    MOI.get(model, MOI.ResultCount()) < 1 && return

    vals.objective[run, 1, i] = calc_objective(model, case.opf)
    vals.objective[run, 2, i] = calc_objective(model, case.opf, case.Pc)
    vals.objective[run, 3, i] = calc_objective(model, case.opf, case.Pcc)
    # @assert sum(vals.objective[run, :, i]) ≈ JuMP.objective_value(model)

    ens = sum_value_property(model, :ls0)
    vals.ENS[run, 1, i] = ens
    vals.LOL[run, 1, i] = abs(ens)>atol ? 1 : 0
    ens = sum_value_property(model, case.Pc, :lsc)
    vals.ENS[run, 2, i] = sum(ens)
    vals.LOL[run, 2, i] = count(x->abs(x)>atol, ens)
    ens = sum_value_property(model, case.Pcc, :lscc)
    vals.ENS[run, 3, i] = sum(ens)
    vals.LOL[run, 3, i] = count(x->abs(x)>atol, ens)

    vals.curtailment[run, 1, i] = sum(sum_value_property(model, :pr0))
    vals.curtailment[run, 2, i] = sum(sum_value_property(model, case.Pc, :prc))
    vals.curtailment[run, 3, i] = sum(sum_value_property(model, case.Pcc, :prcc))

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

function run_type!(result, i, c, type, goal, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    case = Case(opf_base(type, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec)...);
    add_branch_constraints!(case.model, case.pf.ϕ, case.model[:p0], case.brc_up, case.brc_down, case.oplim.branch_rating)
    solve_model!(case.model)
    MOI.get(case.model, MOI.ResultCount()) < 1 && return
    if type != goal
        fix_base_case!(case.model)
        type.C1 && fix_contingencies!(case.model, case.Pc)
        type.C2 && fix_contingencies!(case.model, case.Pcc)
        add_all_contingencies!(goal - type, case)
        solve_model!(case.model);
    end
    gather_run_data!(result, c, i, case)
    return
end

function run_types!(result, i, types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    for (c, type) in enumerate(types)
        run_type!(result, i, c, type, types[end], optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)
    end
    return
end

function run_contingency_select_type!(result, i, c, type, goal, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    case = Case(opf_base(type, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec)...);
    tot_t = constrain_branches!(case.model, case.pf, case.oplim, case.brc_up, case.brc_down, 0.0)
    MOI.get(case.model, MOI.ResultCount()) < 1 && return
    if type != goal
        fix_base_case!(case.model)
        type.C1 && fix_contingencies!(case.model, case.Pc)
        type.C2 && fix_contingencies!(case.model, case.Pcc)
        solve_model!(case.model);
        case, tot_t = run_contingency_select!(goal - type, case)
    end
    gather_run_data!(result, c, i, case)
    return
end

function run_contingency_select_types!(result, i, types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    for (c, type) in enumerate(types)
        run_contingency_select_type!(result, i, c, type, types[end], optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)
    end
    return
end

function run_benders_type!(result, i, c, type, goal, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    case = Case(opf_base(type, sys, optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, time_limit_sec=time_limit_sec)...);
    tot_t = constrain_branches!(case.model, case.pf, case.oplim, case.brc_up, case.brc_down, 0.0)
    MOI.get(case.model, MOI.ResultCount()) < 1 && return
    if type != goal
        fix_base_case!(case.model)
        type.C1 && fix_contingencies!(case.model, case.Pc)
        type.C2 && fix_contingencies!(case.model, case.Pcc)
        solve_model!(case.model);
        case, tot_t = run_benders!(goal - type, case)
    end
    gather_run_data!(result, c, i, case)
    return
end

function run_benders_types!(result, i, types, optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
    for (c, type) in enumerate(types)
        run_benders_type!(result, i, c, type, types[end], optimizer, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)
    end
    return
end

function run_reliability_calculation_benders(types, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long; p_failure=0.00, time_limit_sec=600)
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

        run_benders_types!(result, i, types, optimizer, system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=p_failure, time_limit_sec=time_limit_sec)

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
    
get_lol(mod::Model, cont::Dict, symb::Symbol, atol::Real=1e-14) = 
    sum(count(x->(abs(x)>atol), get_value(mod, getfield(c, symb))) for (i,c) in cont)
