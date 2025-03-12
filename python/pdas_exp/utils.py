import os

def check_meshdir(dirpath):
    assert os.path.isdir(dirpath), f"No mesh directory found at {dirpath}"
    assert os.path.isfile(os.path.join(dirpath, "connectivity.dat"))
    assert os.path.isfile(os.path.join(dirpath, "coordinates.dat"))
    assert os.path.isfile(os.path.join(dirpath, "info.dat"))


def mkdir(dirpath):
    if not os.path.isdir(dirpath):
        os.mkdir(dirpath)


def catchlist(var, dtype, count):
    if not isinstance(var, list):
        assert isinstance(var, dtype)
        var = [var] * count
    assert len(var) == count
    return var

def catchinput(inputdict, inputstr, default_type=None, default_val="nodefault", listval=False):

    # read, apply default if necessary
    try:
        inputval = inputdict[inputstr]
    except KeyError:
        if default_val != "nodefault":
            inputval = default_val
        else:
            raise KeyError("Could not find required input " + inputstr)

    # just break for None input
    if inputval is None:
        return inputval

    # cast to default type
    if (default_type is None) and (default_val != "nodefault"):
        default_type = type(default_val)
    if (default_type is not None) and (default_type != list) and (not isinstance(inputval, list)):
        inputval = default_type(inputval)

    # make into list if it's a list-able input
    if listval:
        if (default_type == list):
            if not isinstance(inputval[0], list):
                inputval = [inputval]
        else:
            if not isinstance(inputval, list):
                inputval = [inputval]

    return inputval