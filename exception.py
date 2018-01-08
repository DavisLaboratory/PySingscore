class error(Exception):
    """
    Base class for other exceptions
    """
    pass

class InvalidNormalisation(error):
    """
    Raised when an invalid normalisation string is entered
    """
    pass

class InvalidIDType(error):
    """
    Raised when the gene identifier types don't match
    """
    pass

class InvalidGrid(error):
    """
    Raised when the grid for graphs is not of sufficient dimensions
    """
    pass