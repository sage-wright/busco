class BatchFatalError(Exception):
    """
    Error that prevents batch mode from running.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


class BuscoError(Exception):
    """
    Module-specific exception
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value
