# Matias Vistnes 2022
# Based on formulae in "Post-contingency corrective control failure: a risk to neglect or a risk to control?"

import psutil
import pyomo.opt as po
import pyomo.environ as pe


class gencoOPT:
    
    def __init__(self):
        self.VOLL          = 100
        self.penalty       = 1E5
                           # prob : contingency 
        self.C             = {0.01:('L1'), 0.02:('L2'), 0.001:('L3'), 0.002:('L4'), 0.02:('L5'), 0.0003:('L1', 'B2'), 0.0004:('L1', 'B3')}
        self.N             = ['B1', 'B2', 'B3', 'B4', 'B5']
        self.L             = ['L1', 'L2', 'L3', 'L4', 'L5']
        self.D             = {'B1':0, 'B2':10,'B3':5, 'B4':10, 'B5':60}
        self.generators    = {'B1':'G1', 'B2':'G2', 'B3':'G3'}
        self.G             = {'G1':40, 'G2':20,'G3':40}
        self.cost          = {'G1': 1, 'G2': 4, 'G3': 6}
                           # line : fBus, tBus, reactance, long-term rating, short-term ratio
        self.L_para        = {'L1':('B1', 'B2', 0.7, 40, 0.7), 
                              'L2':('B2', 'B3', 0.7, 40, 0.7), 
                              'L3':('B3', 'B4', 0.6, 50, 0.7), 
                              'L4':('B3', 'B5', 0.3, 70, 0.7), 
                              'L5':('B1', 'B3', 0.8, 30, 0.7)}
    
    def createOPTModel(self):
        '''
        Establish a model with constraints and objective.
        '''
        m = pe.ConcreteModel()
        
        # Sets
        m.N               = pe.Set(initialize = self.N)
        m.L               = pe.Set(initialize = self.L)
        m.C               = pe.Set(initialize = self.C)
        m.generators      = pe.Set(m.N, initialize = self.generators)
        
        # Create variables and parameters
        m.a             = pe.Param(m.L, self.C, initialize = lambda m,l,c: 0, within = pe.Boolean) #Binary
        m.marginal_cost = pe.Param(m.generators, bounds = lambda m,g: (0, m.cost[g]))
        m.P_max_g       = pe.Param(m.generators, initialize = self.G)  
        m.P_d           = pe.Param(m.N, initialize = self.D)
        m.cost          = pe.Param(m.generators, initialize = self.cost)
        m.VOLL          = pe.Param(initialize = self.VOLL)
        m.X             = pe.Param(m.L, initialize = lambda m,l: self.L_para[l][2])
        m.L_max         = pe.Param(m.L, bounds = lambda m,l: (0, self.L_para[l][3]))
        m.L_max_short   = pe.Param(m.L, bounds = lambda m,l: (0, self.L_para[l][4]))

        m.prod_level      = pe.Var(m.generators, bounds=lambda m,g: (0, m.G[g]), within = pe.NonNegativeReals) # C6
        m.P_flow          = pe.Var(m.L, initialize = lambda m,l: 0, within = pe.Reals)
        m.P_flow_pc       = pe.Var(m.L, m.C, initialize = lambda m,l,c: 0, within = pe.Reals)
        m.angle           = pe.Var(m.L, initialize = lambda m,l: 0, within = pe.Reals)
        m.angle_pc        = pe.Var(m.L, m.C, initialize = lambda m,l,c: 0, within = pe.Reals)
        m.spot_price      = pe.Var(bounds = (0,self.VOLL), within = pe.NonNegativeReals)
        m.ls_d            = pe.Var(m.N, bounds = lambda m,n: (0, m.D[n]), within = pe.NonNegativeReals) # C7
        
        # m.xib             = pe.Var(m.generators, initialize = lambda m,g: 0)
        # m.yib             = pe.Var(m.generators, initialize = lambda m,g: 0)
        # m.gen_dual_min    = pe.Var(m.generators, initialize = lambda m,g: 0)
        # m.gen_dual_max    = pe.Var(m.generators, initialize = lambda m,g: 0)
        # m.binding_min     = pe.Var(m.generators, initialize = lambda m,g: 0, within = pe.Boolean)
        # m.binding_max     = pe.Var(m.generators, initialize = lambda m,g: 0, within = pe.Boolean)
        
        # Objective
        m.OBJ = pe.Objective(expr = sum(m.cost[g]*m.prod_level[g] for g in m.generators) + m.penalty * sum(m.ls_d[n] for n in m.N), sense = pe.minimize)
        
        # Constraints
        m.C2  = pe.ConstraintList()
        for n in m.N:
                m.C2.add(expr = lambda m: sum(m.prod_level[i] if self.L_para[i] == n else 0 for i in m.L) - 
                                sum(m.P_flow if self.L_para[i] == n else 0 for i in m.L) == 
                                sum(m.P_d[i] - m.ls_d[i] if self.L_para[i] == n else 0 for i in m.L))
        m.C3  = pe.Constraint(m.L, expr = lambda m,l: (m.flow[l] - 1/m.X[l]*sum(m.angle[i] if self.L_para[i] == n else 0 for i in m.L) == 0))
        m.C4  = pe.Constraint(m.L, expr = lambda m,l: (m.P_flow[l] <= m.L_max[l]))
        m.C5  = pe.Constraint(m.L, expr = lambda m,l: (-m.P_flow[l] <= m.L_max[l]))
        m.C8  = pe.ConstraintList()
        for n in m.N:
                for c in self.C:
                        m.C8.add(expr = lambda m: sum(m.prod_level[i] if self.L_para[i] == n else 0 for i in m.L) - 
                                        sum(m.P_flow_pc[i,c] if self.L_para[i] == n else 0 for i in m.L) == 
                                        sum(m.P_d[i] - m.ls_d[i] if self.L_para[i] == n else 0 for i in m.L))
        m.C9  = pe.ConstraintList()
        for l in m.L:
                for c in self.C:
                        m.C9.add(expr = lambda m: m.flow_pc[l,c] - m.a[l,c]*1/m.X[l]*sum(m.angle_pc[i,c] if self.L_para[i] == n else 0 for i in m.L) == 0)
        m.C10 = pe.Constraint(m.L, m.C, expr = lambda m,l: (m.P_flow_pc[l,c] <= m.a[l,c] * m.L_max[l]))
        m.C11 = pe.Constraint(m.L, m.C, expr = lambda m,l: (-m.P_flow_pc[l,c] <= m.a[l,c] * m.L_max[l]))

        self.m = m
    
    def solve(self, solver_name, multi_thread = False):
        # Specify solver and solve the problem
        Solver = po.SolverFactory(solver_name)
        if multi_thread: Solver.options['Threads'] = int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2)
        results = Solver.solve(self.m, load_solutions = True)
        if (results.solver.status == po.SolverStatus.ok) and (results.solver.termination_condition == po.TerminationCondition.optimal):
            mod.print_market(disp_all= False)
        else:
            print('No optimal solution found. Try again')
    
    def print_market(self, disp_all):
        if disp_all: self.m.display()
        for c in self.m.gencos:
            if c[0] == 'C':
                print(f'Genco {c} with profit {self.m.genco_profit[c].value:.1f} $/h', end='')
                for g in self.m.generators:
                    print(f'    Gen {g}: {self.m.prod_level[g].value:5.1f} MW, {self.m.price_bid[g].value:7.1f} $/MWh')
            else:
                g_sum = 0
                for g in self.m.generators:
                    g_sum += self.m.prod_level[g].value
                print(f'Load shedding: {g_sum:5.1f} MW')
                
        print(f'Optimal price: {self.m.spot_price.value:.2f} $/MWh')

if __name__ == '__main__':
        mod = gencoOPT()
        mod.createOPTModel()
        mod.solve('gurobi')
        del mod
        print('')