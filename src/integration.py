import cobra
import mewpy.solvers
import pandas as pd
from mewpy.omics import ExpressionSet
from reframed.io.sbml import load_cbmodel
from mewpy.omics.integration.eflux import eFlux
from mewpy.omics.integration.gimme import GIMME
import numpy as np
import os
from cobra.flux_analysis import flux_variability_analysis


class Integration:
    def __init__(self, model, condition=None):
        self.solver = 'gurobi'
        model = cobra.io.read_sbml_model(model)
        self.model = model
        if model:
            print("Model loaded")
            print(model.summary())
        else:
            print("Model not loaded")

    def set_solver(self, solver):
        self.solver = solver
        self.model.solver = solver
        mewpy.solvers.set_default_solver(solver)
        mewpy.solvers.get_default_ode_solver()

    def __del__(self):
        del self.model
        print("Model deleted")

    def gimme(self, expr, conditions, **kwargs):
        """

        :param conditions: list of conditions
        :param expr:
        :return:
        """

        if conditions is None:
            conditions = ['']
        identifiers = expr['Geneid'].tolist()
        expression = expr['tpm'].to_numpy()[:, np.newaxis]
        set_expression = ExpressionSet(identifiers, conditions, expression)
        gimme = GIMME(self.model, set_expression, "e_Biomass__cytop", condition="tpm", parsimonious=True)
        # gimme = GIMME(self.model, expr, "e_Biomass__cytop", condition, parsimonious=True)
        # gimme.save_results(local_path, file_name)
        return gimme

    def eflux(self, expr, conditions, scale_rxn=None, scale_value=1, constraints=None,
              parsimonious=False,
              max_exp=None, **kwargs):
        """

        :param conditions:
        :param file_name:
        :param local_path:
        :param expr:
        :param scale_rxn:
        :param scale_value:
        :param constraints:
        :param parsimonious:
        :param max_exp:
        :param kwargs:
        :return:
        """
        identifiers = expr['Geneid'].tolist()
        expression = expr['tpm'].to_numpy()[:, np.newaxis]
        set_expression = ExpressionSet(identifiers, conditions, expression)
        eflux = eFlux(self.model, set_expression, conditions[0], scale_rxn, scale_value, constraints, parsimonious,
                      max_exp)
        # eflux.save_results(local_path, file_name)
        return eflux

    def pFBA(self):
        """

        :return:
        """
        fba_solution = self.model.optimize()
        pfba_solution = cobra.flux_analysis.pfba(self.model)

        return abs(fba_solution.fluxes["e_Biomass__cytop"] - pfba_solution.fluxes[
            "e_Biomass__cytop"])

    def fva(self, fva=0.95, fraction_of_optimum=0.9, loopless=False, **kwargs):
        """

        :param fva:
        :param loopless: Use the loopless argument to avoid loops (high absolute flux values) that only can be high if
        they are allowed
        to participate in loops
        :param fraction_of_optimum:
        :return:
        """
        float(fraction_of_optimum)
        print('FVA finds the ranges of each metabolic flux at the optimum.')

        print(cobra.flux_analysis.flux_variability_analysis(
            self.model, self.model.reactions[:10], fraction_of_optimum=fraction_of_optimum))

        loop_reactions = [self.model.reactions.FRD7, self.model.reactions.SUCDi]
        flux_variability_analysis(self.model, reaction_list=loop_reactions, loopless=loopless)
        print(f'Loop reactions: {loop_reactions}')
        print('\n')
        print('Running FVA in summary methods')
        print('Model summary. default fva=0.95')
        print(self.model.summary(fva=fva))
        print('Variability in metabolite mass balances with flux variability analysis:')
        print(self.model.metabolites.pyr_c.summary(fva=fva))
        print('Variability in reaction fluxes with flux variability analysis:')
        print(self.model.reactions.EX_pyr_c.summary(fva=fva))


    # def simulation(self):
    #    # build a phenotype simulator
    #    simul = get_simulator(self.model)
    #    return simul.genes[:10]
    #

