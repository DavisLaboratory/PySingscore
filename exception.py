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