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
import re


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

        return abs(fba_solution.fluxes["e_Biomass__cytop"] - pfba_solution.fluxes["e_Biomass__cytop"])

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

        # loop_reactions = [self.model.reactions.FRD7, self.model.reactions.SUCDi]
        reaction_list = self.model.reactions
        flux_variability_analysis(self.model, reaction_list=reaction_list, loopless=loopless)
        # print(f'Loop reactions: {loop_reactions}')
        print('\n')
        print('Running FVA in summary methods')
        print('Model summary. default fva=0.95')
        print(self.model.summary(fva=fva))
        print('Variability in metabolite mass balances with flux variability analysis:')
        # search for omega-3 fatty acids
        print(self.model.metabolites.get_by_id('omega3_fatty_acids').summary(fva=fva))
        print(self.model.metabolites.omega.summary(fva=fva))
        print('Variability in reaction fluxes with flux variability analysis:')
        print(self.model.reactions.EX_pyr_c.summary(fva=fva))

    def biodiesel_fatty_acids(self):
        """ Method to grab the biodiesel fatty acids from the model. Usually fatty acids used in biodisel have between
        14 and 20 carbons """
        for metabolite in self.model.metabolites:
            if 'C14' or 'C15' or 'C16' or 'C17' or 'C18' or 'C19' or 'C20' in metabolite.formula:
                print(f"Potential biodiesel fatty acid: {metabolite.name}, id: {metabolite.id}, formula: "
                      f"{metabolite.formula}, price: {metabolite.shadow_price}")

    def pufas(self):
        """Method to grab potential PUFAs from the model"""
        pufa = ['Linoleic acid', 'alpha-Linolenic acid', 'gamma-Linolenic acid', 'Columbinic acid', 'Stearidonic acid',
                'Mead acid', 'Dihomo-γ-linolenic acid', 'Arachidonic acid', 'Eicosapentaenoic acid', 'Docosapentaenoic '
                                                                                                     'acid',
                'Docosahexaenoic acid', 'Methylene']

        real_pufa = []
        
        for metabolite in self.model.metabolites:
            beg = metabolite.formula[0:3]
            
            if any(x in metabolite.formula for x in pufa) and beg[0] =='C' and beg[1] in '23456789' and beg[2] in \
                    '0123456789':
                real_pufa.append(metabolite)
                print(f"Very Potential PUFA: {metabolite.name}, id: {metabolite.id}, formula: {metabolite.formula},"
                      f"shadow price: {metabolite.shadow_price}")

        for metabolite in self.model.metabolites:
            # search for C20, C21, C22, C23, C24, C25, C26, C27, C28, C29, C30 until C100 in the formula
            beg = metabolite.formula[0:3]
            if beg[0] == 'C' and beg[1] in '23456789' and beg[2] in '0123456789':
                print(f"Potential PUFA: {metabolite.name}, id: {metabolite.id}, formula: {metabolite.formula}, "
                      f"shadow price: {metabolite.shadow_price}")
                            
            #if re.search(r'\bC20\b', metabolite.formula):
            #    print(f"Potential PUFA: {metabolite.name}, id: {metabolite.id}, formula: {metabolite.formula}, "
            #          f" shadow price: {metabolite.shadow_price}")

    def search_for_metabolites(self, search_term):
        """Method to search for metabolites in the model"""
        for metabolite in self.model.metabolites:
            if search_term in metabolite.name:
                print(f"Potential PUFA: {metabolite.name}, id: {metabolite.id}, formula: {metabolite.formula}, "
                      f" shadow price: {metabolite.shadow_price}")

    def most_valuable_metabolites(self):
        """Method to find the most valuable metabolites in the model"""
        valuable = []
        for metabolite in self.model.metabolites:
            if metabolite.shadow_price > 0:
                valuable.append(metabolite)
        valuable.sort(key=lambda x: x.shadow_price, reverse=True)

    def pigments(self):
        """ Method to grab potential pigments from the model """
        list_pigments = ['zeaxanthin', 'Zeaxanthin', 'canthaxanthin', 'Canthaxanthin', 'Astaxanthin', 'astaxanthin',
                         'Chlorophyll a', 'chlorophyll a', 'Chlorophyll b', 'chlorophyll b', 'Chlorophyll c',
                         'chlorophyll c', 'beta-Carotene', 'Violaxanthin', 'violaxanthin', 'Vaucheriaxanthin',
                         'vaucheriaxanthin''canthaxanthin, Canthaxanthin', 'carotenoid', 'Carotenoid']

        for metabolite in self.model.metabolites:
            if any(x in metabolite.name for x in list_pigments):
                print(f"Potential pigment: {metabolite.name}, id: {metabolite.id}, formula: {metabolite.formula}, "
                      f" shadow price: {metabolite.shadow_price}")

    def from_metabolite_get_reaction(self, metabolite):
        """ Method to get the reactions that are associated with a metabolite """
        for reaction in self.model.reactions:
            if metabolite in reaction.metabolites:

                print(f"{metabolite.name} is associated with {reaction.name} and has a id: {reaction.id}")
                print(f"stoichometry of {reaction.metabolites[metabolite]}, lower bound: {reaction.lower_bound}, "
                      f"upper bound: {reaction.upper_bound}")

