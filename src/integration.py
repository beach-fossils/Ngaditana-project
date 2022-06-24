import cobra
import mewpy.solvers
import pandas as pd
from mewpy.omics import ExpressionSet
from reframed.io.sbml import load_cbmodel
from mewpy.omics.integration.eflux import eFlux
from mewpy.omics.integration.gimme import GIMME
import numpy as np
import os


class Integration:
    def __init__(self, model, condition=None):
        model = cobra.io.read_sbml_model(model)
        self.model = model
        if model:
            print("Model loaded")
            print(model.summary())
        else:
            print("Model not loaded")

    def __del__(self):
        del self.model
        print("Model deleted")

    def gimme(self, expr, condition, local_path, file_name, **kwargs):
        """

        :param file_name:
        :param local_path:
        :param condition: list of conditions
        :param expr:
        :return:
        """

        if condition is None:
            condition = ''
        identifiers = expr['Geneid'].tolist()
        expression = expr['tpm'].to_numpy()[:, np.newaxis]
        set_expression = ExpressionSet(identifiers, condition, expression)
        gimme = GIMME(self.model, set_expression, "e_Biomass__cytop", condition="tpm", parsimonious=True)
        # gimme = GIMME(self.model, expr, "e_Biomass__cytop", condition, parsimonious=True)
        # gimme.save_results(local_path, file_name)
        return gimme

    def eflux(self, expr, local_path, file_name, conditions, scale_rxn=None, scale_value=1, constraints=None,
              parsimonious=False,
              max_exp=None, **kwargs):
        """

        :param file_name:
        :param local_path:
        :param expr:
        :param condition:
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
        eflux = eFlux(self.model, set_expression, conditions[0], scale_rxn, scale_value, constraints, parsimonious, max_exp)
        # eflux.save_results(local_path, file_name)
        return eflux

    # def simulation(self):
    #    # build a phenotype simulator
    #    simul = get_simulator(self.model)
    #    return simul.genes[:10]
    #
