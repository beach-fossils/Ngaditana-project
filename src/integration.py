#from reframed.io.sbml import load_cbmodel
from mewpy.simulation import get_simulator

from mewpy.omics import eFlux, GIMME


class Integration:
    def __init__(self, model):
        self.model = model
        #load_cbmodel(self.model, flavor='cobra')

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
        return eFlux(self.model, expr, condition, scale_rxn, scale_value, constraints, parsimonious, max_exp, **kwargs)

    # def simulation(self):
    #    # build a phenotype simulator
    #    simul = get_simulator(self.model)
    #    return simul.genes[:10]
    #
