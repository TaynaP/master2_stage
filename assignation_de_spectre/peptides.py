class Peptides(object):
    """Classe pour représenter les peptides donnés par RPG
    """

    def __init__(self, seq, mass, nodes):
        self.__seq = seq
        self.__mass = mass
        self.__nodes = nodes

    # Getter
    def get_seq(self):
        return self.__seq

    def get_mass(self):
        return self.__mass
    
    def get_nodes(self):
        return self.__nodes
    
    # Setter
    def set_nodes(self, liste):
        self.__nodes = liste
    
    def set_mass(self, mass):
        self.__mass = mass
    
    def append_mass(self, mass):
        self.__mass.append(mass)
    

    def __str__(self):
        return f"seq : {self.__seq}; mass : {self.__mass}; nodes : {self.__nodes} "

    def __repr__(self):
        return f'Peptide(seq={self.__seq}, mass={self.__mass}, nodes={self.__nodes})'