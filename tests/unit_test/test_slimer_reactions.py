from unittest import TestCase
from unittest.mock import Mock

from SLIMErpy.io_slimer.cobra import read_mat_model
from SLIMErpy.slime_reactions_creator import SLIMEReactionsGenerator, Table, ChainData, LipidTable


class TestSlimerReactionGenerator(TestCase):

    def test_get_new_index(self):
        model = read_mat_model("yeast780.mat")
        generator = SLIMEReactionsGenerator(model, cytoplasm_id="c")

        last_index_plus_one = generator.get_new_metabolite_index()
        self.assertEqual("3721", last_index_plus_one)


    def test_add_backbone(self):

        csv_file = "\\Users\\HP-PC\\projeto\\SLIMErpy\\SLIMEr\\data\\fullData_Ejsing2009.csv"

        table = Table(metabolite_names_field=0,
                      abundances_fields=slice(1, -1, 2),
                      sd_deviations_fields=slice(2, -1, 2))
        table.from_csv(csv_file)

        model = Mock()
        generator = SLIMEReactionsGenerator(model)

        generator.add_backbones(table, True)


    def test_chaindata(self):

        csv_file = "\\Users\\HP-PC\\projeto\\SLIMErpy\\SLIMEr\\data\\chainData_Lahtvee2017.csv"

        chaindata = ChainData(name_field=0,
                              formula_field=1,
                              abundances_fields=slice(2, -1, 2),
                              sd_deviations_fields=slice(3, -1, 2) )

        chaindata.from_csv(csv_file)

        self.assertEqual(chaindata.name[0], "C16:0 chain")
        self.assertEqual(chaindata.abundances[0, 0], 0.00808584)

    def test_table(self):

        csv_file = "\\Users\\HP-PC\\projeto\\SLIMErpy\\SLIMEr\\data\\fullData_Ejsing2009.csv"

        table =  Table(metabolite_names_field=0,
                       abundances_fields=slice(1, -1, 2),
                       sd_deviations_fields=slice(2, -1, 2))

        table.from_csv(csv_file)

        self.assertEqual(table.abundances[0, 0], 0.030568242)
        self.assertEqual(table.sd_deviations[0,0], 0.009517803)



    def test_lipidtable(self):

        csv_file = "\\Users\\HP-PC\\projeto\\SLIMErpy\\SLIMEr\\data\\lipidData_Lahtvee2017.csv"

        lipidtable = LipidTable(short_name_field=0,
                                name_field=1,
                                ID_in_model_field=2,
                                abundancies_fields=slice(3, -1, 1))

        lipidtable.from_csv(csv_file)

        self.assertEqual(lipidtable.name[0], "1-phosphatidyl-1D-myo-inositol")
        self.assertEqual(lipidtable.abundancies[0, 0], 0.003637)



