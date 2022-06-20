import cobra


def read_mat_model(file_path: str) -> cobra.Model:
    """
    Read matlab model

    :param str file_path: path for the model file
    :return cobra.Model: the model object from cobra
    """
    model = cobra.io.load_matlab_model(file_path)

    return model


def write_mat_model(model: cobra.Model, file_path: str):
    """

    :param cobra.Model model: model object to be written
    :param str file_path: path of the file to be written
    :return:
    """

    cobra.io.save_matlab_model(model, file_name=file_path)
    
    
def read_sbml_model(file_path: str) -> cobra.Model:
    """

    :param str file_path: path for the model file
    :return cobra.Model: the model object from cobra
    """
    model = cobra.io.read_sbml_model(file_path)
    return model


def write_sbml_model(model: cobra.Model, file_path: str):
    """

    :param cobra.Model model: model object to be written
    :param str file_path: path of the file to be written
    :return:
    """
    cobra.io.write_sbml_model(model, filename=file_path)
