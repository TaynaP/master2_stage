class Peptides(object):
    """Classe pour représenter les peptides donnés par RPG
    """

    def __init__(self, seq, mass, nodes, label = None, positions = None, target_seq = None, target_mass = None):
        self.__seq = seq
        self.__mass = mass
        self.__nodes = nodes #les noms des séquences où il apparait => species
        self.__label = label
        self.__positions = positions
        self.__target_seq = target_seq
        self.__target_mass = target_mass

    # Getter
    def get_seq(self):
        return self.__seq

    def get_mass(self):
        return self.__mass
    
    def get_nodes(self):
        return self.__nodes
    
    def get_label(self):
        return self.__label
    
    def get_positions(self):
        return self.__positions
    
    def get_target_seq(self):
        return self.__target_seq
    
    def get_target_mass(self):
        return self.__target_mass
    
    # Setter
    def set_nodes(self, liste):
        self.__nodes = liste
    
    def set_positions(self, liste):
        self.__positions = liste
    
    def set_target_seq(self, liste):
        self.__target_seq = liste
    
    def set_target_mass(self, liste):
        self.__target_mass = liste
    
    def set_mass(self, mass):
        self.__mass = mass
    
    def append_mass(self, mass):
        self.__mass.append(mass) #dans le cas des PTM => utile en présence de PTM
    

    def __str__(self):
        return f"seq : {self.__seq}; mass : {self.__mass}; nodes : {self.__nodes}, label : {self.__label}, positions : {self.__positions}, target_seq : {self.__target_seq}, target_mass : {self.__target_mass}"

    def __repr__(self):
        return f'Peptide(seq={self.__seq}, mass={self.__mass}, nodes={self.__nodes}, label={self.__label}, positions={self.__positions}, target_seq={self.__target_seq}, target_mass={self.target_mass})'