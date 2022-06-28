import re

import cobra
import pandas as pd
from cobra import Metabolite, Reaction


class SLIMEReactionsGenerator:

    def __init__(self, model: cobra.Model, cytoplasm_id: str = "c"):
        self.model = model
        self.cytoplasm_id = cytoplasm_id

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
        new_metabolite_name = metabolite_name[:x.start() - 1]
        model_backbone_name = None
        for backbone_name in group1.keys():
            if new_metabolite_name.startswith(backbone_name):
                model_backbone_name = backbone_name + " " + compartment_name

        # vamos verificar se o nome do metabolito é igual ao nome do fatty acid name (que é a chave do grupo)
        # e depois vamos ver se o compartimento é o mesmo
        for fatty_acid_name in group2:
            if fatty_acid_name == new_metabolite_name and group2[fatty_acid_name][1] == compartment_name:
                model_backbone_name = group2[fatty_acid_name][0] + " " + compartment_name

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
        new_metabolite_name = metabolite_name[:x.start() - 1]
        new_id = 's_' + last_index + '[' + compartment_name + ']'

        new_metabolite = Metabolite(
            new_id,
            formula=formula,
            name=new_metabolite_name,
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

    def __main(self):
        metabolite_names = [metabolite.name for metabolite in self.model.metabolites]

        for metabolite in self.model.metabolites:
            backbone_name = self.get_backbone_name(metabolite)
            if backbone_name:
                if backbone_name in metabolite_names:
                    if not backbone_name.startswith('ergosterol ['):
                        metabolite.formula = ""  # verificar se não falta nada dentro de parêntesis
                    else:
                        self.add_lipid_species(backbone_name, metabolite.compartment, "", False)

                if "cytoplasm" not in backbone_name:
                    x = re.search("\[[a-zA-Z0-9_ ]*\]", backbone_name)
                    backbone_without_compartiment = backbone_name[:x.start()]
                    backbone_with_cytoplasm = backbone_without_compartiment + '[cytoplasm]'

                    if backbone_with_cytoplasm in metabolite_names:
                        metabolite.formula = ""
                    else:
                        self.add_lipid_species(backbone_with_cytoplasm, metabolite.compartment, "", False)

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

    def add_backbones(self, table, include_tails):
        met_names = table.metabolites_names  # vai receber a tabela, para ir buscar o nome dos metabolitos
        new_met_names = met_names  # ...
        met_ids = [] #vai ao modelo ver qual é o id do metabolito com aquele nome, e dps adiciona-o a uma lista
        metabolite = None
        for i in range(len(new_met_names)):
            new_met_names[i] = new_met_names[i] + ' [cytoplasm]'
#get_metabolite_by_name é um método da classe, por isso para o chamarmos usamos o self.
            metabolite = self.get_metabolite_by_name(new_met_names[i])
            met_ids.append(metabolite.id)

        if metabolite is not None:
            self.add_lipid_species('lipid - backbones [cytoplasm]', metabolite.compartment, '', include_tails)

        self.change_lipid_composition(table)

    def add_metabolites_and_abundacies_to_dictionary(self,lipid_table, pseudo_backbone_metabolite):
        dic_add_metabolites = {}
        for i in range(len(lipid_table.ID_in_model)):
            metabolite_ID = lipid_table.ID_in_model[i]
            metabolite = self.model.metabolites.get_by_id(metabolite_ID)
            abundance = lipid_table.abundancies.iloc[i, 0]  #depois ver melhor maneira de definir condições ambientais
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

            dic_add_metabolites = self.add_metabolites_and_abundacies_to_dictionary(lipid_table, pseudo_backbone_metabolite)

            pseudo_reaction.add_metabolites(dic_add_metabolites)
            self.model.add_reactions([pseudo_reaction])
        else:
            new_reaction_ID = 'r_2108'
            pseudo_backbone_metabolite = self.model.metabolites.get_by_id('s_1096[c]')
            lipid_pseudoreaction = self.model.reactions.get_by_id(new_reaction_ID)

            dictionary_with_metabolites_and_abundances = self.add_metabolites_and_abundacies_to_dictionary(lipid_table, pseudo_backbone_metabolite)
            lipid_pseudoreaction.add_metabolites(dictionary_with_metabolites_and_abundances)

    def add_chains(self, chain_data, include_tails):
        new_metabolite_names = chain_data.name
        mets = None
        for i in range(len(new_metabolite_names)):
            new_metabolite_names[i] = new_metabolite_names[i] + ' [cytoplasm]'
            mets = self.add_lipid_species(new_metabolite_names[i], self.cytoplasm_id, chain_data.formula[i],
                                          not include_tails)

        if include_tails:
            self.add_lipid_species('lipid - tails [cytoplasm]', self.cytoplasm_id,'', include_tails)
            self.change_chain_composition(chain_data)

    def change_chain_composition(self, chain_data):
        tail_name = chain_data.name
        metabolites_for_reaction = {}
        for i in range(len(tail_name)):
            tail_name[i] = tail_name[i] + ' [cytoplasm]'
            tail = self.get_metabolite_by_name(tail_name[i])
            metabolites_for_reaction[tail] = - chain_data.abundances.iloc[i, 0] # definir depois condiçoes ambientais

        new_id = self.get_new_reactions_index()
        reaction_id = 'r_' + new_id
        reaction_name = 'lipid pseudoreaction - tail'
        tail_pseudo_metabolite = self.get_metabolite_by_name('lipid - tails [cytoplasm]')
        #tail_id = self.

        pseudo_reaction_tails = Reaction(reaction_id,
                                         name=reaction_name,
                                         lower_bond=0,
                                         upper_bound=1000)

        metabolites_for_reaction[tail_pseudo_metabolite] = 1

        pseudo_reaction_tails.add_metabolites(metabolites_for_reaction)


lipid_id = 's_1096[c]'

reaction_name = 'lipid pseudoreaction - merge'
add_reaction = Reaction(add_reaction, 'r_2108',
                        'reactionName', reaction_name,


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


