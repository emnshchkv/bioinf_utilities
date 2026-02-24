from abc import ABC, abstractmethod

class BiologicalSequence(ABC):
    
    def __init__(self, sequence):
        self.sequence = sequence
                
    def __len__(self):
        return len(self.sequence)
        
    @property
    @abstractmethod
    def ALPHABET(self):
        pass
    
    def check_alphabet(self):
        return set(self.sequence).issubset(self.ALPHABET)
    
    @abstractmethod
    def __str__(self):
        pass
    
    @abstractmethod
    def __getitem__(self, index):
        pass
    
class NucleicAcidSequence(BiologicalSequence):
    
    def __str__(self):
        ...
        
    def __getitem__(self, index):
        ...
    
    def complement(self):
        ...
        
    def reverse(self):
        return self.sequence[::-1]
        
    def reverse_complement(self):
        ...
    
    
class DNASequence(NucleicAcidSequence):
    ALPHABET = set('ATGCatgc')
    
    def transcribe(self):
        ...
    
class RNASequence(NucleicAcidSequence):
    ALPHABET = set('AUGCaugc')
    ...
    
class AminoAcidSequence(BiologicalSequence):
    ...