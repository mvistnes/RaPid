using .SCOPF
using PowerSystems

system = SCOPF.get_system("ELK14.json")
voll, prob, contingencies = setup(system)
opfm, imml = SCOPF.run_benders2(system, voll, prob)