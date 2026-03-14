"""
Created on Jul 19, 2013
Modified on Jul 13, 2023

This class represents the OptiType model for MHC typing based on NGS data

It is dependent on Coopr and uses an external ILP solver such as GLPK or CPLEX

@original author: Benjamin Schubert
@follow-up editors: JJY ZYX

                
"""

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
from pyomo.environ import ConcreteModel, Set, Param, Var, Binary, Objective, Constraint, ConstraintList, maximize
from pyomo.opt import SolverFactory, TerminationCondition
from collections import defaultdict
import pandas as pd
import itertools


class OptiType(object):
    """
    classdocs

    """

    def __init__(self, cov, occ, allele_groups, beta, t_max_allele=2, solver="glpk", threads=1,
                 verbosity=0):
        """
        Constructor
        """

        self.__beta = float(beta)
        self.__t_max_allele = t_max_allele
        self.__solver = SolverFactory(solver)
        self.__threads = threads
        self.__opts = {"threads": threads} if threads > 1 else {}
        self.__verbosity = verbosity
        self.__changed = True  # model needs to know if it changed from last run or not
        self.__ks = 1
        self.__allele_groups = allele_groups


        #loci = loci_alleles

        self.__allele_to_type = {allele: allele_type for allele_type, group in allele_groups.items() for allele in
                                   group}

        '''
            generates the basic ILP model
        '''

        model = ConcreteModel()

        # init Sets
        model.LociNames = Set(initialize=list(allele_groups.keys()))
        model.Loci = Set(model.LociNames, initialize=lambda m, l: allele_groups[l])

        L = list(itertools.chain(*list(allele_groups.values()))) ## L: allele ids
        R = set([r for (r, _) in list(cov.keys())])  ## R: read ids
        model.L = Set(initialize=L)
        #reconst = {allele_id: 0.01 for allele_id in L if '_' in allele_id}
        model.R = Set(initialize=R)

        # init Params
        model.cov = Param(model.R, model.L, initialize=lambda model, r, a: cov.get((r, a), 0))
        #model.reconst = Param(model.L, initialize=lambda model, a: reconst.get(a, 0))

        model.occ = Param(model.R, initialize=occ)
        model.t_allele = Param(initialize=self.__t_max_allele, mutable=True)

        model.beta = Param(initialize=self.__beta,
                           validate=lambda val, model: 0.0 <= float(self.__beta) <= 0.999,
                           mutable=True)
        model.nof_loci = Param(initialize=len(allele_groups))

        # init variables
        model.x = Var(model.L, domain=Binary)
        model.y = Var(model.R, domain=Binary)

        model.re = Var(model.R, bounds=(0.0, None))
        model.hetero = Var(bounds=(0.0, model.nof_loci))

        # init objective
        #model.read_cov = Objective(
            #rule=lambda model: sum(model.occ[r] * (model.y[r] - model.beta * (model.re[r])) for r in model.R) - sum(
                #model.reconst[a] * model.x[a] for a in model.L), sense=maximize)
        model.read_cov = Objective(
            rule=lambda model: sum(model.occ[r] * (model.y[r] - model.beta * (model.re[r])) for r in model.R), sense=maximize)

                

        # init Constraints
        model.max_allel_selection = Constraint(model.LociNames, rule=lambda model, l: sum(
            model.x[a] for a in model.Loci[l]) <= model.t_allele)
        model.min_allel_selection = Constraint(model.LociNames,
                                               rule=lambda model, l: sum(model.x[a] for a in model.Loci[l]) >= 1)
        model.is_read_cov = Constraint(model.R,
                                       rule=lambda model, r: sum(model.cov[r, a] * model.x[a] for a in model.L) >=
                                                             model.y[r])
        model.heterozygot_count = Constraint(
            rule=lambda model: model.hetero >= sum(model.x[a] for a in model.L) - model.nof_loci)

        # regularization constraints
        model.reg1 = Constraint(model.R, rule=lambda model, r: model.re[r] <= model.nof_loci * model.y[r])
        model.reg2 = Constraint(model.R, rule=lambda model, r: model.re[r] <= model.hetero)
        model.reg3 = Constraint(model.R,
                                rule=lambda model, r: model.re[r] >= model.hetero - model.nof_loci * (1 - model.y[r]))

        # generate constraint list for solution enumeration
        model.c = ConstraintList()
        # Generate instance. Used to be .create() but deprecated since,
        # as ConcreteModels are instances on their own now.
        self.__instance = model

    def set_beta(self, beta):
        """
            Sets the parameter beta
        """
        self.__changed = True
        getattr(self.__instance, str(self.__instance.beta)).set_value(float(beta))

    def set_t_max_allele(self, t_max_allele):
        """
            Sets the upper bound of alleles selected per loci
        """
        self.__changed = True
        getattr(self.__instance, str(self.__instance.t_allele)).set_value(t_max_allele)

    def solve(self, ks):
        """
            solves the problem k times and discards the found solutions in the next run.
        """
        d = defaultdict(list)  # in there we store the typing +objective and generate afterwards a DatarFrame with it

        if self.__changed or self.__ks != ks:
            self.__ks = ks
            for k in range(ks):
                #expr = 0

                self.__instance.preprocess()
                try:
                    res = self.__solver.solve(self.__instance, options=self.__opts, tee=self.__verbosity)
                except:
                    print ("WARNING: Solver does not support multi-threading. Please change the config"
                          " file accordingly. Falling back to single-threading.")
                    res = self.__solver.solve(self.__instance, options={}, tee=self.__verbosity)
                self.__instance.solutions.load_from(res)  # solution loading changed recently.

                # if self.__verbosity > 0:
                #     res.write(num=1)

                if res.solver.termination_condition != TerminationCondition.optimal:
                    print("Optimal solution hasn't been obtained. This is a terminal problem.")  # TODO message, exit
                    break

                selected = []
                #indices = []
                #encountered_4digit = []
                for j in self.__instance.x:
                    #if self.__allele_to_4digit[j][0] in 'HJG':
                        #if 0.99 <= self.__instance.x[j].value <= 1.01:
                            #selected.append(j)
                        #indices.append(j)
                        #continue
                    if 0.99 <= self.__instance.x[j].value <= 1.01:
                        selected.append(j)
                        #exp_i = 0
                        #exp_i += self.__instance.x[j]

                        #if self.__allele_to_4digit[j] in encountered_4digit:
                            #continue
                        #encountered_4digit.append(self.__allele_to_4digit[j])


                        #for i_allele in self.__groups_4digit[self.__allele_to_4digit[j]]:
                            #if self.__instance.x[i_allele].value <= 0:
                                #exp_i += self.__instance.x[i_allele]
                            #indices.append(i_allele)
                        #expr += (1.0 - exp_i)
                #zero_indices = set([j for j in self.__instance.x]).difference(set(indices))
                #for j in zero_indices:
                    #expr += self.__instance.x[j]

                #self.__instance.c.add(expr >= 1)

                # if self.__verbosity > 0:
                #     print selected
                #     self.__instance.c.pprint()

                print(selected)

                aas = [self.__allele_to_type[x] for x in selected]
                c = dict.fromkeys(aas, 1)
                for i in range(len(aas)):
                    if aas.count(aas[i]) < 2:  # homo
                        d[aas[i] + "1"].append(selected[i])
                        d[aas[i] + "2"].append(selected[i])
                    else: # hetero
                        d[aas[i] + str(c[aas[i]])].append(selected[i])
                        c[aas[i]] += 1

                nof_reads = sum((self.__instance.occ[j] * self.__instance.y[j].value for j in self.__instance.y))
                # if self.__verbosity > 0:
                #     print "Obj", res.Solution.Objective.__default_objective__.Value
                d['obj'].append(self.__instance.read_cov())
                d['nof_reads'].append(nof_reads)

            self.__instance.c.clear()
            self.__changed = False
            self.__enumeration = pd.DataFrame(d)

            # self.__rank()
            return self.__enumeration
        else:
            return self.__enumeration

   
