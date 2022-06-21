import pandas as pd
from reframed.io.sbml import load_cbmodel
from mewpy.omics.integration.eflux import eFlux
from mewpy.omics.integration.gimme import GIMME
import numpy as np
import os


class Integration:
    def __init__(self, model):
        model = load_cbmodel(model)
        self.model = model

        if model:
            print("Model loaded")
            print(model.summary())
        else:
            print("Model not loaded")

    def __del__(self):
        del self.model
        print("Model deleted")

    def gimme(self, expr, biomass=None, condition=0, cutoff=0.25, growth_frac=0.9, constraints=None, parsimonious=False,
              **kwargs):
        """

        :param expr:
        :param biomass:
        :param condition:
        :param cutoff:
        :param growth_frac:
        :param constraints:
        :param parsimonious:
        :param kwargs:
        :return:
        """
        expr = pd.read_csv(expr)
        expr = pd.DataFrame(expr)
        print(expr)
        return GIMME(self.model, expr, biomass, condition, cutoff, growth_frac, constraints, parsimonious, **kwargs)

    def eflux(self, expr, condition=0, scale_rxn=None, scale_value=1, constraints=None, parsimonious=False,
              max_exp=None, **kwargs):
        """

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
        expr = pd.read_csv(expr, sep="\t")
        return eFlux(self.model, expr, condition, scale_rxn, scale_value, constraints, parsimonious, max_exp, **kwargs)

    # def simulation(self):
    #    # build a phenotype simulator
    #    simul = get_simulator(self.model)
    #    return simul.genes[:10]
    #
