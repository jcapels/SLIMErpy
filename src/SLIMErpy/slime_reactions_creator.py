import re

import cobra
import pandas as pd
from cobra import Metabolite, Reaction
from pandas import Series


class Table:  # fulldata2009
    def __init__(self, metabolite_names_field, abundances_fields, sd_deviations_fields):
        self.metabolite_names_field = metabolite_names_field
        self.abundances_fields = abundances_fields
        self.sd_deviations_fields = sd_deviations_fields

    def from_csv(self, csv_file):
        self.data = pd.read_csv(csv_file)
        self.metabolites_names = list(self.data.iloc[:, self.metabolite_names_field])
        self.abundances = self.data.iloc[:, self.abundances_fields]  # dataframe por isso não clocamos list()
        self.sd_deviations = self.data.iloc[:, self.sd_deviations_fields]


class LipidTable:  # lipidData2017
    def __init__(self, short_name_field, name_field, ID_in_model_field, abundancies_fields):
        self.short_name_field = short_name_field
        self.name_field = name_field
        self.ID_in_model_field = ID_in_model_field
        self.abundancies_fields = abundancies_fields

    def from_csv(self, csv_file):
        self.data = pd.read_csv(csv_file)
        self.short_name = list(self.data.iloc[:, self.short_name_field])
        self.name = list(self.data.iloc[:, self.name_field])
        self.ID_in_model = list(self.data.iloc[:, self.ID_in_model_field])
        self.abundancies = self.data.iloc[:, self.abundancies_fields]


class ChainData:  # ChainData2017
    def __init__(self, name_field, formula_field, abundances_fields, sd_deviations_fields):
        self.name_field = name_field
        self.formula_field = formula_field
        self.abundances_fields = abundances_fields
        self.sd_deviations_fields = sd_deviations_fields

    def from_csv(self, csv_file):
        self.data = pd.read_csv(csv_file)
        self.name = list(self.data.iloc[:, self.name_field])
        self.formula = list(self.data.iloc[:, self.formula_field])
        self.abundances = self.data.iloc[:, self.abundances_fields]
        self.sd_deviations = self.data.iloc[:, self.sd_deviations_fields]


class CompData:  # compData2017
    def __init__(self, names_field, ID_in_model_field, abundances_fields):
        self.names_field = names_field
        self.ID_in_model_field = ID_in_model_field
        self.abundances_fields = abundances_fields

    def from_csv(self, csv_file):
        self.data = pd.read_csv(csv_file)
        self.names = list(self.data.iloc[:, self.names_field])
        self.ID_in_model = list(self.data.iloc[:, self.ID_in_model_field])
        self.abundances = self.data.iloc[:, self.abundances_fields]


class FluxData:  # fluxData2017
    def __init__(self, name_field, ID_in_model_field, abundances_fields, sd_deviations_fields):
        self.name_field = name_field
        self.ID_in_model_field = ID_in_model_field
        self.abundances_fields = abundances_fields
        self.sd_deviations_fields = sd_deviations_fields

    def from_csv(self, csv_file):
        self.data = pd.read_csv(csv_file)
        self.name = list(self.data.iloc[:, self.name_field])
        self.abundances = self.data.iloc[:, self.abundances_fields]
        self.sd_deviations = self.data.iloc[:, self.sd_deviations_fields]


class SLIMEReactionsGenerator:

    def __init__(self, model: cobra.Model, backbone_table, chains_table,
                 generic_lipid_id,
                 cytoplasm_id: str = "c"):

        self.tail_metabolites = []
        self.model = model
        self.cytoplasm_id = cytoplasm_id
        self.backbone_table = backbone_table
        self.chains_table = chains_table
        self.generic_lipid_id = generic_lipid_id

    @staticmethod
    def get_backbone_name(metabolite: cobra.Metabolite) -> str:
        """
        Taken from SLIMEr/models/getBackboneName.m
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
        metabolite_name = metabolite.name
        # sugestão: usar expressões regulares: https://regexr.com/

        x = re.search("\[[a-zA-Z0-9_ ]*\]", metabolite_name)
        compartment_name = metabolite_name[x.start():x.end()]
        compartment_name = compartment_name.replace('[', '').replace(']', '')
        new_metabolite_name = metabolite_name[:x.start() - 1]
        model_backbone_name = None
        for backbone_name in group1.keys():
            if new_metabolite_name.startswith(backbone_name) and new_metabolite_name != backbone_name and \
                    compartment_name == group1[backbone_name] and "phosphate" not in new_metabolite_name:
                model_backbone_name = backbone_name + " [" + compartment_name + "]"
                break

        # vamos verificar se o nome do metabolito é igual ao nome do fatty acid name (que é a chave do grupo)
        # e depois vamos ver se o compartimento é o mesmo
        for fatty_acid_name in group2:
            if fatty_acid_name == new_metabolite_name and group2[fatty_acid_name][1] == compartment_name:
                model_backbone_name = group2[fatty_acid_name][0] + " [" + compartment_name + "]"
                break
        return model_backbone_name

    def get_new_metabolite_index(self):
        metabolites_list = self.model.metabolites
        list_of_identifiers = []
        for metabolite in metabolites_list:
            identifier = re.sub('[^(\d*)]', '', metabolite.id)
            identifier = int(identifier)
            list_of_identifiers.append(identifier)
        max_id = max(list_of_identifiers)
        max_id += 1
        max_id = str(max_id)
        return max_id

    def get_new_reactions_index(self):
        reactions_list = self.model.reactions
        list_of_identifiers = []
        for reactions in reactions_list:
            identifier = re.sub('[^(\d*)]', '', reactions.id)
            identifier = int(identifier)
            list_of_identifiers.append(identifier)
        max_id = max(list_of_identifiers)
        max_id += 1
        max_id = str(max_id)
        return max_id

    def add_lipid_species(self, metabolite_name, compartment_name, formula, exchange):
        last_index = self.get_new_metabolite_index()
        x = re.search("\[[a-zA-Z0-9_ ]*\]", metabolite_name)
        # new_metabolite_name = metabolite_name[:x.start() - 1]
        new_id = 's_' + last_index + '[' + compartment_name + ']'

        new_metabolite = Metabolite(
            new_id,
            formula=formula,
            name=metabolite_name,
            compartment=compartment_name
        )
        if exchange:
            new_reaction_index = self.get_new_reactions_index()
            reaction_name = metabolite_name + ' exchange'
            exchange_reaction = Reaction('r_' + new_reaction_index,
                                         name=reaction_name,
                                         lower_bound=0,
                                         upper_bound=1000)
            exchange_reaction.add_metabolites({new_metabolite: -1})  # -1 pq éum reagente. 1 se fosse produto

            # defini a reação e o metabolito, falta adicionar ao modelo:

            self.model.add_reactions([exchange_reaction])

        else:
            self.model.add_metabolites([new_metabolite])

        return new_metabolite

    def get_metabolite_by_name(self, name):
        for metabolite in self.model.metabolites:
            if metabolite.name == name:
                return metabolite

    def add_backbones(self, include_tails):
        met_names = self.backbone_table.name  # vai receber a tabela, para ir buscar o nome dos metabolitos
        new_met_names = met_names  # ...
        met_ids = []  # vai ao modelo ver qual é o id do metabolito com aquele nome, e dps adiciona-o a uma lista
        for i in range(len(new_met_names)):
            new_met_names[i] = new_met_names[i] + ' [cytoplasm]'
            # get_metabolite_by_name é um método da classe, por isso para o chamarmos usamos o self.
            metabolite = self.get_metabolite_by_name(new_met_names[i])
            met_ids.append(metabolite.id)

        self.backbone_table.ID_in_model = Series(met_ids)

        self.add_lipid_species('lipid - backbones [cytoplasm]', self.cytoplasm_id, '', include_tails)

        self.change_lipid_composition(self.backbone_table)

    def add_metabolites_and_abundacies_to_dictionary(self, lipid_table, pseudo_backbone_metabolite):
        dic_add_metabolites = {}
        for i in range(len(lipid_table.ID_in_model)):
            metabolite_ID = lipid_table.ID_in_model[i]
            metabolite = self.model.metabolites.get_by_id(metabolite_ID)
            abundance = lipid_table.abundancies.iloc[i, 0]  # depois ver melhor maneira de definir condições ambientais
            dic_add_metabolites[metabolite] = - abundance

        dic_add_metabolites[pseudo_backbone_metabolite] = 1

        return dic_add_metabolites

    def change_lipid_composition(self, lipid_table):
        pseudo_backbone_metabolite = self.get_metabolite_by_name('lipid - backbones [cytoplasm]')
        if pseudo_backbone_metabolite is not None:
            new_reaction_index = self.get_new_reactions_index()
            new_reaction_ID = 'r_' + new_reaction_index
            reaction_name = 'lipid pseudoreaction - backbone'

            pseudo_reaction = Reaction(new_reaction_ID,
                                       name=reaction_name,
                                       lower_bound=0,
                                       upper_bound=1000)

            dic_add_metabolites = self.add_metabolites_and_abundacies_to_dictionary(lipid_table,
                                                                                    pseudo_backbone_metabolite)

            pseudo_reaction.add_metabolites(dic_add_metabolites)
            self.model.add_reactions([pseudo_reaction])
        else:
            new_reaction_ID = 'r_2108'
            pseudo_backbone_metabolite = self.model.metabolites.get_by_id('s_1096[c]')
            lipid_pseudoreaction = self.model.reactions.get_by_id(new_reaction_ID)

            dictionary_with_metabolites_and_abundances = self.add_metabolites_and_abundacies_to_dictionary(lipid_table,
                                                                                                           pseudo_backbone_metabolite)
            lipid_pseudoreaction.add_metabolites(dictionary_with_metabolites_and_abundances)

    def add_chains(self, include_tails):
        new_metabolite_names = self.chains_table.name
        for i in range(len(new_metabolite_names)):
            new_metabolite_names[i] = new_metabolite_names[i] + ' [cytoplasm]'
            self.add_lipid_species(new_metabolite_names[i], self.cytoplasm_id, self.chains_table.formula[i],
                                   not include_tails)

        if include_tails:
            self.add_lipid_species('lipid - tails [cytoplasm]', self.cytoplasm_id, '', include_tails)
            self.change_chain_composition(self.chains_table)

    def change_chain_composition(self, chain_data):
        tail_name = chain_data.name
        metabolites_for_reaction = {}
        for i in range(len(tail_name)):
            tail_name[i] = tail_name[i]
            tail = self.get_metabolite_by_name(tail_name[i])
            metabolites_for_reaction[tail] = - chain_data.abundances.iloc[i, 0]  # definir depois condiçoes ambientais
            self.tail_metabolites.append(tail)

        new_id = self.get_new_reactions_index()
        reaction_id = 'r_' + new_id
        reaction_name = 'lipid pseudoreaction - tail'
        tail_pseudo_metabolite = self.get_metabolite_by_name('lipid - tails [cytoplasm]')
        # tail_id = self.

        pseudo_reaction_tails = Reaction(reaction_id,
                                         name=reaction_name,
                                         lower_bound=0,
                                         upper_bound=1000)

        metabolites_for_reaction[tail_pseudo_metabolite] = 1

        pseudo_reaction_tails.add_metabolites(metabolites_for_reaction)
        self.model.add_reactions([pseudo_reaction_tails])

    @staticmethod
    def get_reaction_reactants_and_products(reaction):
        reactants = []
        products = []
        for metabolite in reaction.metabolites:
            if reaction.metabolites[metabolite] < 0:
                reactants.append(metabolite)
            else:
                products.append(metabolite)

        return reactants, products

    @staticmethod
    def get_atoms_in_formula(formula):
        # one uppercase letter followed by optional lowercase letters
        # followed by zero or more digits
        rx = re.compile(r'([A-Z][a-z]*)(\d*)')

        return {k: int(count) if count else 1 for k, count in rx.findall(formula)}

    def get_stoichiometry_from_formula(self, formulas, atom):
        stoichiometries = []
        for formula in formulas:
            atom_counts = self.get_atoms_in_formula(formula)
            if atom in atom_counts:
                stoichiometries.append(atom_counts[atom])

            else:
                stoichiometries.append(0)

        return stoichiometries

    def get_metabolite_by_startswith_name(self, metabolite_name):
        for metabolite in self.model.metabolites:
            if metabolite.name.startswith(metabolite_name):
                return metabolite

    def add_slimer_reaction(self, reaction, structurally_defined_lipid):

        if structurally_defined_lipid is None:
            reactants, products = self.get_reaction_reactants_and_products(reaction)
            structurally_defined_lipid = reactants[0]
            reaction_id = reaction.id

        if reaction is None:
            reaction_id = "r_" + self.get_new_reactions_index()

        backbone_name = self.get_backbone_name(structurally_defined_lipid)

        fatty_acids = {
                       'palmitate': 'C16:0',
                       'palmitoleate': 'C16:1',
                       'stearate': 'C18:0',
                       'oleate': 'C18:1',
                       'ergosteryl palmitoleate': 'C16:1',
                       'ergosteryl oleate': 'C18:1',
                       'sphinganine': 'C18:0',
                       'phytosphingosine': 'C18:0',
                       'sphinganine 1-phosphate': 'C18:0',
                       'phytosphingosine 1-phosphate': 'C18:0'}

        to_delete = False
        if backbone_name is None:
            # i = structurally_defined_lipid.name.index("[")
            # structurally_defined_lipid_name = structurally_defined_lipid.name[:i - 1]
            #
            # if structurally_defined_lipid_name in fatty_acids:
            #     backbone_name = structurally_defined_lipid.name
            # else:
            to_delete = True
            return to_delete

        i = backbone_name.index("[")
        backbone_name = backbone_name[:i - 1]

        i = structurally_defined_lipid.name.index("[")
        structurally_defined_lipid_name = structurally_defined_lipid.name[:i - 1]
        case_1 = ['1-phosphatidyl-1D-myo-inositol',
                  'sn-2-acyl-1-lysophosphatidylinositol',
                  'phosphatidyl-L-serine',
                  'phosphatidylcholine',
                  'phosphatidylethanolamine',
                  'phosphatidate',
                  'phosphatidylglycerol',
                  'cardiolipin',
                  'diglyceride',
                  'triglyceride']

        case_2 = ['ceramide',
                  'inositol-P-ceramide',
                  'inositol phosphomannosylinositol phosphoceramide',
                  'mannosylinositol phosphorylceramide']

        case_3 = [
            'fatty acid',
            'ergosterol ester',
            'long-chain base',
            'long-chain base phosphate'
        ]

        tail_reactions = []
        # if reaction is not None:
        #     if any([case_name in reaction.name for case_name in case_1]):
        #         matches = re.finditer(":", structurally_defined_lipid_name)
        #
        #         for match in matches:
        #             tail_reactions.append("C" + structurally_defined_lipid_name[match.start() - 2: match.start() + 2])
        #
        #     elif any([case_name in reaction.name for case_name in case_2]):
        #         if "(C24)" in backbone_name:
        #             tail_reactions.extend(['C18:0', 'C24:0'])
        #         elif "(C26)" in backbone_name:
        #             tail_reactions.extend(['C18:0', 'C26:0'])
        #
        #     elif any([case_name in reaction.name for case_name in case_3]):
        #         tail_reactions.append(fatty_acids[structurally_defined_lipid_name])
        #
        # else:
        if backbone_name in case_1:
            matches = re.finditer(":", structurally_defined_lipid_name)

            for match in matches:
                tail_reactions.append("C" + structurally_defined_lipid_name[match.start() - 2: match.start() + 2])

        elif backbone_name in case_2:
            if "(C24)" in backbone_name:
                tail_reactions.extend(['C18:0', 'C24:0'])
            elif "(C26)" in backbone_name:
                tail_reactions.extend(['C18:0', 'C26:0'])

        elif backbone_name in case_3:
            tail_reactions.append(fatty_acids[structurally_defined_lipid_name])

        tails_coefficients = {}
        prod_formulas = []
        for rxn in tail_reactions:
            tail_name = rxn + " chain [cytoplasm]"
            tail_metabolite = self.get_metabolite_by_name(tail_name)
            tail_coef = tail_metabolite.formula_weight
            tails_coefficients[tail_metabolite] = tail_coef
            prod_formulas.append(tail_metabolite.formula)

        formula = ""
        atoms = ['C', 'H', 'N', 'O', 'P', 'S']

        for atom in atoms:
            Nin = self.get_stoichiometry_from_formula([structurally_defined_lipid.formula], atom)
            Nout = self.get_stoichiometry_from_formula(prod_formulas, atom)
            diff = sum(Nin) - sum(Nout)
            if diff == 1:
                formula = formula + atom
            elif diff > 1:
                formula = formula + atom + str(diff)

        backbone_name = backbone_name + " ["
        backbone_metabolite = self.get_metabolite_by_startswith_name(backbone_name)
        backbone_metabolite.formula = formula

        reaction_name = structurally_defined_lipid.name + ' SLIME rxn'

        slimer_reaction = Reaction(id=reaction_id,
                                   name=reaction_name,
                                   lower_bound=0,
                                   upper_bound=1000,
                                   )
        if reaction is not None:
            reaction.remove_from_model()

        formula_weight = structurally_defined_lipid.formula_weight
        formula_weight = formula_weight / 1000
        metabolites_to_add = {structurally_defined_lipid: -1,
                              backbone_metabolite: formula_weight}

        metabolites_to_add.update(tails_coefficients)
        slimer_reaction.add_metabolites(metabolites_to_add)
        self.model.add_reactions([slimer_reaction])
        return to_delete

    def get_reaction_by_name(self, name) -> Reaction:
        for rxn in self.model.reactions:
            if rxn.name == name:
                return rxn

    def get_fraction(self, mw_data, compound_type):
        reaction_name = compound_type + " pseudoreaction"
        if compound_type == "lipid":
            reaction = self.get_reaction_by_name(reaction_name + " - backbone")
            reactants, products = self.get_reaction_reactants_and_products(reaction)
            coefficients = reaction.get_coefficients(reactants)
            F = -sum(coefficients)
        else:
            reaction = self.get_reaction_by_name(reaction_name)
            mw_data_compound_type = mw_data[mw_data['compound_type'] == compound_type]
            F = 0
            for i, row in mw_data_compound_type.iterrows():
                metabolite = self.model.metabolites.get_by_id(row['id_in_model'])
                if metabolite in reaction.metabolites:
                    coefficient = reaction.get_coefficient(metabolite)
                    abundance = -coefficient * (row['mw'] - 18) / 1000
                    F += abundance

        return F

    def sum_biomass(self, mw_data):
        X = 0
        P = self.get_fraction(mw_data, "protein")
        X += P
        L = self.get_fraction(mw_data, "lipid")
        X += L
        R = self.get_fraction(mw_data, "RNA")
        X += R
        D = self.get_fraction(mw_data, "DNA")
        X += D
        C = self.get_fraction(mw_data, "carbohydrate")
        X += C

        biomass_reaction = self.model.reactions.get_by_id('r_4041')
        for i, row in mw_data.iterrows():
            metabolite = self.model.metabolites.get_by_id(row['id_in_model'])
            if metabolite in biomass_reaction.metabolites:
                coefficient = biomass_reaction.get_coefficient(metabolite)
                abundance = - coefficient * row['mw'] / 1000
                X += abundance

        return X, P, C, R, D, L

    def get_metabolites_contain_name(self, name):
        metabolites = []
        for metabolite in self.model.metabolites:
            if name in metabolite.name:
                metabolites.append(metabolite)
        return metabolites

    def change_metabolite_coefficients(self, coefficients, fraction, metabolites, reaction):
        j = 0
        for coefficient in coefficients:
            metabolite = metabolites[j]
            if metabolite in reaction.metabolites:
                coefficient_fraction = coefficient * fraction
                reaction.subtract_metabolites({metabolite: coefficient})
                reaction.add_metabolites({metabolite: coefficient_fraction})
            j += 1

    def change_other_metabolites_composition(self, data: CompData, mw_data: pd.DataFrame):
        X, P, C, R, D, L = self.sum_biomass(mw_data)
        protein_pseudo_reaction = self.get_reaction_by_name('protein pseudoreaction')
        rna_pseudo_reaction = self.get_reaction_by_name('RNA pseudoreaction')

        for i in range(len(data.ID_in_model)):
            metabolite_id = data.ID_in_model[i]
            if metabolite_id == "protein":
                protein_fraction = data.abundances.iloc[i, 0] / P
                amino_acids = self.get_metabolites_contain_name('tRNA')
                amino_acids_in_reaction = []
                for amino_acid in amino_acids:
                    if amino_acid in protein_pseudo_reaction.metabolites:
                        amino_acids_in_reaction.append(amino_acid)
                coefficients = protein_pseudo_reaction.get_coefficients(amino_acids_in_reaction)
                self.change_metabolite_coefficients(coefficients, protein_fraction, amino_acids,
                                                    protein_pseudo_reaction)
            elif metabolite_id == "RNA":
                rna_fraction = data.abundances.iloc[i, 0] / R

                nucs = mw_data[mw_data['compound_type'] == 'RNA'].loc[:, "id_in_model"].tolist()
                nucs = [self.model.metabolites.get_by_id(nuc) for nuc in nucs]

                nucs_in_reaction = []
                for nuc in nucs:
                    if nuc in protein_pseudo_reaction.metabolites:
                        nucs_in_reaction.append(nuc)
                coefficients = rna_pseudo_reaction.get_coefficients(nucs_in_reaction)

                self.change_metabolite_coefficients(coefficients, rna_fraction, nucs,
                                                    rna_pseudo_reaction)

            else:
                metabolite = self.model.metabolites.get_by_id(metabolite_id)
                abundance = data.abundances.iloc[i, 0]
                compound_type = mw_data.loc[mw_data['id_in_model'] == metabolite_id, "compound_type"].values[0]
                if compound_type == "carbohydrate":
                    reaction = self.get_reaction_by_name('carbohydrate pseudoreaction')
                    mw = mw_data.loc[mw_data['id_in_model'] == metabolite_id, "mw"].values[0] - 18
                elif compound_type == "DNA":
                    reaction = self.get_reaction_by_name('DNA pseudoreaction')
                    mw = mw_data.loc[mw_data['id_in_model'] == metabolite_id, "mw"].values[0] - 18
                elif compound_type == "N":
                    reaction = self.get_reaction_by_name('biomass pseudoreaction')
                    mw = mw_data.loc[mw_data['id_in_model'] == metabolite_id, "mw"].values[0]

                if metabolite in reaction.metabolites:
                    coefficient = reaction.get_coefficient(metabolite)
                    coefficient_altered = - abundance / (mw * 1000)
                    reaction.subtract_metabolites({metabolite: coefficient})
                    reaction.add_metabolites({metabolite: coefficient_altered})

        X, P, C, R, D, L = self.sum_biomass(mw_data)
        delta = X - 1
        metabolite_ids = mw_data[mw_data["compound_type"] == "carbohydrate"].iloc[:, 0].values
        metabolites = [self.model.metabolites.get_by_id(metabolite_id) for metabolite_id in metabolite_ids]
        mass_pre = C
        mass_post = mass_pre - delta
        carbohydrate_fraction = mass_post / mass_pre
        carbohydrate_reaction = self.get_reaction_by_name('carbohydrate pseudoreaction')

        carbohydrate_in_reaction = []
        for metabolite in metabolites:
            if metabolite in carbohydrate_reaction.metabolites:
                carbohydrate_in_reaction.append(metabolite)

        coefficients = carbohydrate_reaction.get_coefficients(carbohydrate_in_reaction)

        self.change_metabolite_coefficients(coefficients,
                                            carbohydrate_fraction,
                                            metabolites,
                                            carbohydrate_reaction
                                            )
        X, P, C, R, D, L = self.sum_biomass(mw_data)
        GAMpol = P * 37.7 + C * 12.8 + R * 26.0 + D * 26.0
        return GAMpol

    def run(self, include_tails):
        metabolite_names = [metabolite.name for metabolite in self.model.metabolites]

        for metabolite in self.model.metabolites:
            backbone_name = self.get_backbone_name(metabolite)
            if backbone_name:
                if backbone_name in metabolite_names:
                    if not backbone_name.startswith('ergosterol ['):
                        backbone_metabolite = self.get_metabolite_by_name(backbone_name)
                        backbone_metabolite.formula = ""  # verificar se não falta nada dentro de parêntesis
                else:
                    self.add_lipid_species(backbone_name, metabolite.compartment, "", False)
                    metabolite_names.append(backbone_name)

                if "cytoplasm" not in backbone_name:
                    x = re.search("\[[a-zA-Z0-9_ ]*\]", backbone_name)
                    backbone_without_compartiment = backbone_name[:x.start()]
                    backbone_with_cytoplasm = backbone_without_compartiment + '[cytoplasm]'

                    if backbone_with_cytoplasm in metabolite_names:
                        backbone_metabolite = self.get_metabolite_by_name(backbone_with_cytoplasm)
                        backbone_metabolite.formula = ""
                    else:
                        self.add_lipid_species(backbone_with_cytoplasm, self.cytoplasm_id, "", False)

                    backbone_metabolite = self.get_metabolite_by_name(backbone_name)  # objeto metabolito
                    cytoplasm_metabolite = self.get_metabolite_by_name(backbone_with_cytoplasm)  # objeto metabolito
                    transport_id = 'r_' + self.get_new_reactions_index()
                    transport_name = backbone_without_compartiment + 'transport'

                    transport_reaction = Reaction(transport_id,
                                                  name=transport_name,
                                                  lower_bound=0,
                                                  upper_bound=1000)
                    metabolite_dictionary = {backbone_metabolite: -1, cytoplasm_metabolite: 1}
                    transport_reaction.add_metabolites(metabolite_dictionary)
                    reaction_list = [transport_reaction]
                    self.model.add_reactions(reaction_list)

        self.add_backbones(include_tails=include_tails)
        self.add_chains(include_tails=include_tails)

        reactions_to_remove = []
        print(self.model.reactions.get_by_id("r_3982"))
        # Add SLIME reactions replacing existing ISA rxns
        reactions_in_model = self.model.reactions.copy()
        for reaction in reactions_in_model:
            to_delete = False
            if "isa " in reaction.name:
                to_delete = self.add_slimer_reaction(reaction, None)

            elif "complex sphingolipid transport" == reaction.name:
                to_delete = True

            if to_delete:
                reactions_to_remove.append(reaction)

        self.model.remove_reactions(reactions_to_remove)

        # Add new SLIME reactions (with no corresponding ISA rxns):
        for metabolite in self.model.metabolites:
            backbone_name = self.get_backbone_name(metabolite)
            if backbone_name is not None:
                self.add_slimer_reaction(None, metabolite)

        backbone_metabolite = self.get_metabolite_by_name('lipid - backbones [cytoplasm]')
        generic_lipid = self.model.metabolites.get_by_id(self.generic_lipid_id)
        if include_tails:
            tail_metabolite = self.get_metabolite_by_name('lipid - tails [cytoplasm]')

            metabolite_coefficient = {

                backbone_metabolite: -1,
                tail_metabolite: -1,
                generic_lipid: 1

            }

        else:

            metabolite_coefficient = {

                backbone_metabolite: -1,
                generic_lipid: 1

            }

        lipid_pseudoreaction_merge = Reaction(

            id='r_2108',
            name='lipid pseudoreaction - merge',
            lower_bound=0,
            upper_bound=1000

        )

        lipid_pseudoreaction_merge.add_metabolites(metabolite_coefficient)
        lipid_pseudo_reaction = self.model.reactions.get_by_id('r_2108')
        lipid_pseudo_reaction.remove_from_model()
        self.model.add_reactions([lipid_pseudoreaction_merge])
