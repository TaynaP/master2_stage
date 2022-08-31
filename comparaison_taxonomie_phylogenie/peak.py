class Peak(object):
    """Classe pour représenter les pics des spectres de masse
    """

    def __init__(self, mass, intensity):
        self.__mass = mass
        self.__intensity = intensity

    # Getter
    def get_mass(self):
        return self.__mass
    
    def get_intensity(self):
        return self.__intensity
    

    def __str__(self):
        return f"mass : {self.__mass}; intensity : {self.__intensity} "

    def __repr__(self):
        return f'Peak(mass={self.__mass}, intensity={self.__intensity})'