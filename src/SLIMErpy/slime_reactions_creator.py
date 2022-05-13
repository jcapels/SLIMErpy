import cobra


class SLIMEReactionsGenerator:

    def __init__(self, model: cobra.Model):
        self.model = model

    @staticmethod
    def get_backbone_name(metabolite: cobra.Metabolite) -> str:
        """
        Taken from getBackboneName.m
        :param metabolite:
        :return:
        """
        group1 = {'1-phosphatidyl-1D-myo-inositol': 'cytoplasm',
                  'sn-2-acyl-1-lysophosphatidylinositol': 'endoplasmic reticulum membrane',
                  'phosphatidyl-L-serine': 'endoplasmic reticulum membrane',
                  'phosphatidylcholine': 'endoplasmic reticulum membrane',
                  'phosphatidylethanolamine': 'endoplasmic reticulum membrane',
                  'phosphatidate': 'endoplasmic reticulum membrane',
                  'diglyceride': 'endoplasmic reticulum membrane',
                  'triglyceride': 'endoplasmic reticulum membrane',
                  'phosphatidylglycerol': 'mitochondrial membrane',
                  'cardiolipin': 'mitochondrial membrane',
                  'ceramide': 'Golgi',
                  'inositol-P-ceramide': 'Golgi',
                  'inositol phosphomannosylinositol phosphoceramide': 'Golgi',
                  'mannosylinositol phosphorylceramide': 'Golgi'}

        group2 = {'palmitate': ('fatty acid', 'cytoplasm'),
                  'palmitoleate': ('fatty acid', 'cytoplasm'),
                  'stearate': ('fatty acid', 'cytoplasm'),
                  'oleate': ('fatty acid', 'cytoplasm'),
                  'ergosterol': ('ergosterol', 'cytoplasm'),
                  'ergosteryl palmitoleate': ('ergosterol ester', 'endoplasmic reticulum membrane'),
                  'ergosteryl oleate': ('ergosterol ester', 'endoplasmic reticulum membrane'),
                  'sphinganine': ('long-chain base', 'endoplasmic reticulum'),
                  'phytosphingosine': ('long-chain base', 'endoplasmic reticulum'),
                  'sphinganine 1-phosphate': ('long-chain base phosphate', 'endoplasmic reticulum'),
                  'phytosphingosine 1-phosphate': ('long-chain base phosphate', 'endoplasmic reticulum')
                  }
        pass

    def add_lipid_species(self, backbone_name, formula, exchange):
        """
        Implementation of addLipidSpecies.m

        :param backbone_name:
        :param formula:
        :param exchange:
        :return:
        """
        pass

    def __main(self):
        metabolite_names = [metabolite.name for metabolite in self.model.metabolites]
        for metabolite in self.model.metabolites:
            backbone_name = self.get_backbone_name(metabolite)
            if backbone_name:
                if backbone_name in metabolite_names:
                    if not backbone_name.startswith('ergosterol ['):
                        metabolite.formula = ""
                    else:
                        self.add_lipid_species(backbone_name, "", False)

                if "cytoplasm" in backbone_name:
                    pass
