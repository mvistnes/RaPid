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
    Pc::Dict{<:Integer,Main.SCOPF.ExprC}, Pcc::Dict{<:Integer,Main.SCOPF.ExprCC}, ramp_mult::Real
)
    vals.objective[run, 1, i] = calc_objective(model, opfs)
    vals.objective[run, 2, i] = calc_objective(model, opfs, Pc, ramp_mult)
    vals.objective[run, 3, i] = calc_objective(model, opfs, Pcc, ramp_mult)
    # @assert sum(vals.objective[run, :, i]) ≈ JuMP.objective_value(model)

    ens = sum_value_property(model, :ls0)
    vals.ENS[run, 1, i] = ens
    vals.LOL[run, 1, i] = abs(ens)>1e-6 ? 1 : 0
    ens = sum_value_property(model, Pc, :lsc)
    vals.ENS[run, 2, i] = sum(ens)
    vals.LOL[run, 2, i] = count(x->abs(x)>1e-6, ens)
    ens = sum_value_property(model, Pcc, :lscc)
    vals.ENS[run, 3, i] = sum(ens)
    vals.LOL[run, 3, i] = count(x->abs(x)>1e-6, ens)

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

function run_cases!(result, i, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long)
    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    MOI.get(model, MOI.ResultCount()) < 1 && return
    idx = SCOPF.get_nodes_idx(opf.nodes);
    list = SCOPF.make_list(opf, idx, opf.nodes);
    SCOPF.fix_base_case(model)
    model, opf, pf, oplim, Pc, Pcc, Pccx = add_all_contingencies!(SCOPF.SC, SCOPF.PCSC, opf, oplim, model, list, pf, idx, Pc, Pcc, Pccx)
    SCOPF.solve_model!(model);
    gather_run_data!(result, 1, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    model, opf, pf, oplim, Pc, Pcc, Pccx = add_all_contingencies!(SCOPF.PSC, SCOPF.PCSC, opf, oplim, model, list, pf, idx, Pc, Pcc, Pccx)
    SCOPF.solve_model!(model);
    SCOPF.fix_base_case(model)
    SCOPF.fix_short_term(model, Pc)
    model, opf, pf, oplim, Pc, Pcc, Pccx = add_all_contingencies!(SCOPF.SC, SCOPF.PSC, opf, oplim, model, list, pf, idx, Pc, Pcc, Pccx)
    SCOPF.solve_model!(model);
    gather_run_data!(result, 2, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.PSC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    SCOPF.fix_base_case(model)
    SCOPF.fix_short_term(model, Pc)
    model, opf, pf, oplim, Pc, Pcc, Pccx = add_all_contingencies!(SCOPF.PSC, SCOPF.PCSC, opf, oplim, model, list, pf, idx, Pc, Pcc, Pccx)
    SCOPF.solve_model!(model);
    gather_run_data!(result, 3, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.PCSC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    gather_run_data!(result, 4, i, model, opf, Pc, Pcc, ramp_mult)
end

function run_contingency_select_cases!(result, i, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=0.00)
    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    MOI.get(model, MOI.ResultCount()) < 1 && return
    SCOPF.fix_base_case(model)
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_contingency_select!(SCOPF.SC, SCOPF.PCSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    gather_run_data!(result, 1, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_contingency_select!(SCOPF.PSC, SCOPF.PCSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    SCOPF.fix_base_case(model)
    SCOPF.fix_long_term(model, Pcc)
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_contingency_select!(SCOPF.SC, SCOPF.PSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    gather_run_data!(result, 2, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_contingency_select(SCOPF.PSC, sys, Gurobi.Optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
    SCOPF.fix_base_case(model)
    SCOPF.fix_short_term(model, Pc)
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_contingency_select!(SCOPF.PSC, SCOPF.PCSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    gather_run_data!(result, 3, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_contingency_select(SCOPF.PCSC, sys, Gurobi.Optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
    gather_run_data!(result, 4, i, model, opf, Pc, Pcc, ramp_mult)
    return result
end

function run_benders_cases!(result, i, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=0.00)
    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    SCOPF.solve_model!(model);
    MOI.get(model, MOI.ResultCount()) < 1 && return
    SCOPF.fix_base_case(model)
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_benders!(SCOPF.SC, SCOPF.PCSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    gather_run_data!(result, 1, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx = SCOPF.opf(SCOPF.SC, sys, Gurobi.Optimizer, voll=voll, contingencies=cont, prob=prob, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long);
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_benders!(SCOPF.PSC, SCOPF.PCSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    SCOPF.fix_base_case(model)
    SCOPF.fix_long_term(model, Pcc)
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_benders!(SCOPF.SC, SCOPF.PSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    gather_run_data!(result, 2, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(SCOPF.PSC, sys, Gurobi.Optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
    SCOPF.fix_base_case(model)
    SCOPF.fix_short_term(model, Pc)
    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = run_benders!(SCOPF.PSC, SCOPF.PCSC, model, opf, pf, oplim, Pc, Pcc, Pccx)
    gather_run_data!(result, 3, i, model, opf, Pc, Pcc, ramp_mult)

    model, opf, pf, oplim, Pc, Pcc, Pccx, tot_t = SCOPF.run_benders(SCOPF.PCSC, sys, Gurobi.Optimizer, voll, prob, cont, max_shed=max_shed,
        ramp_mult=ramp_mult, ramp_minutes=ramp_minutes, short_term_multi=short, long_term_multi=long, p_failure=p_failure);
    gather_run_data!(result, 4, i, model, opf, Pc, Pcc, ramp_mult)
    return result
end

function run_reliability_calculation_benders(system, voll, prob, contingencies, max_shed, ramp_mult, ramp_minutes, short, long, p_failure=0.00)
    # demands = sort_components!(get_demands(system))
    # pd = demands .|> get_active_power
    hours = read_x_data("data\\ieee_std_load_profile.txt")
    hours = 0.5:0.1:1.3

    result = SystemRunData(4, length(hours))
    Logging.disable_logging(Logging.Info)
    buffer_sys = [deepcopy(system) for _ in 1:Threads.nthreads()]
    # Threads.@threads for i in eachindex(hours)
    for i in eachindex(hours)
        h = hours[i]
        sys = buffer_sys[Threads.threadid()]
        set_active_power!.(get_components(StaticLoad, sys), (get_components(StaticLoad, sys) .|> get_max_active_power) * h)
        cont = sort_components!(get_branches(sys))

        run_benders_cases!(result, i, sys, voll, prob, cont, max_shed, ramp_mult, ramp_minutes, short, long, p_failure)

        print(i, " ")
    end
    Logging.disable_logging(Logging.Debug)
    # set_active_power!.(get_components(StaticLoad, system), (get_components(StaticLoad, system) .|> get_max_active_power))
    return result
end

read_x_data(fname) = vec(readdlm(fname,Float64))

get_objective(mod::Model) = termination_status(mod) != MOI.OPTIMAL ? -1.0 : objective_value(mod)

get_ens(mod::Model, cont::Dict, symb::Symbol) = 
    sum(sum(SCOPF.get_value(mod, getfield(c, symb))) for (i,c) in cont)
    
get_lol(mod::Model, cont::Dict, symb::Symbol) = 
    sum(count(x->(abs(x)>1e-5), SCOPF.get_value(mod, getfield(c, symb))) for (i,c) in cont)