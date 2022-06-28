import os
from unittest import TestCase
from unittest.mock import Mock

import pandas as pd
from cobra import Metabolite

from SLIMErpy.io_slimer.cobra import read_mat_model, write_sbml_model, read_sbml_model
from SLIMErpy.slime_reactions_creator import SLIMEReactionsGenerator, Table, ChainData, LipidTable, CompData
from tests import TEST_DIR


class TestSlimerReactionGenerator(TestCase):

    def test_get_new_index(self):
        model = read_mat_model("../data/yeast780.mat")
        generator = SLIMEReactionsGenerator(model, backbone_table=None, chains_table=None, generic_lipid_id="s_1096[c]",
                                            cytoplasm_id="c")

        last_index_plus_one = generator.get_new_metabolite_index()
        self.assertEqual("3721", last_index_plus_one)

    def test_add_backbone(self):
        csv_file = os.path.join(TEST_DIR, "data/lipidData_Lahtvee2017.csv")

        table = LipidTable(short_name_field=0,
                           name_field=1,
                           ID_in_model_field=2,
                           abundancies_fields=slice(3, -1, 1))
        table.from_csv(csv_file)

        model = read_mat_model(os.path.join(TEST_DIR, "data/yeast780.mat"))
        generator = SLIMEReactionsGenerator(model, backbone_table=table, chains_table=None,
                                            generic_lipid_id="s_1096[c]",
                                            cytoplasm_id="c")

        generator.add_backbones(True)

        metabolite = generator.get_metabolite_by_name("triglyceride [cytoplasm]")
        self.assertIsInstance(metabolite, Metabolite)
        self.assertTrue("[c]" in generator.backbone_table.ID_in_model[0])

    def test_add_chains(self):
        csv_file = os.path.join(TEST_DIR, "data/chainData_Lahtvee2017.csv")

        chaindata = ChainData(name_field=0,
                              formula_field=1,
                              abundances_fields=slice(2, -1, 2),
                              sd_deviations_fields=slice(3, -1, 2))

        chaindata.from_csv(csv_file)

        model = read_mat_model(os.path.join(TEST_DIR, "data/yeast780.mat"))
        generator = SLIMEReactionsGenerator(model, backbone_table=None, chains_table=chaindata,
                                            generic_lipid_id="s_1096[c]",
                                            cytoplasm_id="c")

        generator.add_chains(True)

        metabolite = generator.get_metabolite_by_name("lipid - tails [cytoplasm]")
        self.assertIsInstance(metabolite, Metabolite)

    def test_chaindata(self):
        csv_file = "\\Users\\HP-PC\\projeto\\SLIMErpy\\SLIMEr\\data\\chainData_Lahtvee2017.csv"

        chaindata = ChainData(name_field=0,
                              formula_field=1,
                              abundances_fields=slice(2, -1, 2),
                              sd_deviations_fields=slice(3, -1, 2))

        chaindata.from_csv(csv_file)

        self.assertEqual(chaindata.name[0], "C16:0 chain")
        self.assertEqual(chaindata.abundances[0, 0], 0.00808584)

    def test_comp_data(self):
        csv_file = os.path.join(TEST_DIR, "data/compData_Lahtvee2017.csv")

        comp_data = CompData(names_field=0,
                             ID_in_model_field=1,
                             abundances_fields=slice(2, -1)
                             )

        comp_data.from_csv(csv_file)
        self.assertEqual(comp_data.abundances.iloc[0, 0], 0.050958)

    def test_table(self):
        csv_file = "\\Users\\HP-PC\\projeto\\SLIMErpy\\SLIMEr\\data\\fullData_Ejsing2009.csv"

        table = Table(metabolite_names_field=0,
                      abundances_fields=slice(1, -1, 2),
                      sd_deviations_fields=slice(2, -1, 2))

        table.from_csv(csv_file)

        self.assertEqual(table.abundances[0, 0], 0.030568242)
        self.assertEqual(table.sd_deviations[0, 0], 0.009517803)

    def test_lipidtable(self):
        csv_file = "\\Users\\HP-PC\\projeto\\SLIMErpy\\SLIMEr\\data\\lipidData_Lahtvee2017.csv"

        lipidtable = LipidTable(short_name_field=0,
                                name_field=1,
                                ID_in_model_field=2,
                                abundancies_fields=slice(3, -1, 1))

        lipidtable.from_csv(csv_file)

        self.assertEqual(lipidtable.name[0], "1-phosphatidyl-1D-myo-inositol")
        self.assertEqual(lipidtable.abundancies[0, 0], 0.003637)

    def test_main(self):
        csv_file = os.path.join(TEST_DIR, "data/lipidData_Lahtvee2017.csv")

        table = LipidTable(short_name_field=0,
                           name_field=1,
                           ID_in_model_field=2,
                           abundancies_fields=slice(3, -1, 1))
        table.from_csv(csv_file)

        csv_file = os.path.join(TEST_DIR, "data/chainData_Lahtvee2017.csv")

        chaindata = ChainData(name_field=0,
                              formula_field=1,
                              abundances_fields=slice(2, -1, 2),
                              sd_deviations_fields=slice(3, -1, 2))

        chaindata.from_csv(csv_file)

        model = read_mat_model(os.path.join(TEST_DIR, "data/yeast780.mat"))
        generator = SLIMEReactionsGenerator(model, chains_table=chaindata, backbone_table=table,
                                            cytoplasm_id="c", generic_lipid_id="s_1096[c]")

        csv_file = os.path.join(TEST_DIR, "data/lipidData_Lahtvee2017.csv")

        table = LipidTable(short_name_field=0,
                           name_field=1,
                           ID_in_model_field=2,
                           abundancies_fields=slice(3, -1, 1))
        table.from_csv(csv_file)

        generator.run(True)

        csv_file = os.path.join(TEST_DIR, "data/compData_Lahtvee2017.csv")

        comp_data = CompData(names_field=0,
                             ID_in_model_field=1,
                             abundances_fields=slice(2, -1)
                             )

        comp_data.from_csv(csv_file)

        write_sbml_model(generator.model, "yeast_lipidomics.xml")

        mw_data = pd.read_csv(os.path.join(TEST_DIR, "data/mw.csv"))
        generator.change_other_metabolites_composition(comp_data, mw_data)
        write_sbml_model(generator.model, "yeast_lipidomics2.xml")

    def test_model_optimization(self):
        model = read_mat_model(os.path.join(TEST_DIR, "data/yeast780.mat"))
        print(model.optimize())
        model = read_sbml_model(os.path.join(TEST_DIR, "yeast_lipidomics.xml"))
        print(model.objective)
        print(model.optimize())
