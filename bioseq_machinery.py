from abc import ABC, abstractmethod

class BiologicalSequence(ABC):
    
    def __init__(self, sequence):
        self.sequence = sequence
                
    def __len__(self):
        return len(self.sequence)
        
    @property
    @abstractmethod
    def ALPHABET(self) -> set[str]:
        pass
    
    def check_alphabet(self):
        return set(self.sequence).issubset(self.ALPHABET)
    
    @abstractmethod
    def __str__(self) -> str:
        pass
    
    @abstractmethod
    def __getitem__(self, index):
        pass
    
class NucleicAcidSequence(BiologicalSequence):
    
    def __init__(self, sequence):
        if type(self) is NucleicAcidSequence:
            raise NotImplementedError(
                "Direct instantiation of the Nucleic Acid Sequence class is not allowed."
            )
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError(f"Sequence contains invalid symbols. Allowed symbols: {self.ALPHABET}")
    
    def __str__(self):
        return "Oligonucleotide : " +self.sequence
        
    def __getitem__(self, index):
        return self.sequence[index]
    
    def complement(self):
        replication_dict = {
        "A": "T",
        "a": "t",
        "T": "A",
        "t": "a",
        "G": "C",
        "g": "c",
        "C": "G",
        "c": "g",
        "U": "A",
        "u": "a",
        }
        return self.__class__(''.join(replication_dict[nucleotide] for nucleotide in self.sequence))
        
    def reverse(self):
        return self.__class__(self.sequence[::-1])
        
    def reverse_complement(self):
        return self.complement().reverse()
    
    
class DNASequence(NucleicAcidSequence):
    @property
    def ALPHABET(self):
        return set('ATGCatgc')
    
    def transcribe(self):
        transcription_dict = {
        "A": "U",
        "a": "u",
        "T": "A",
        "t": "a",
        "G": "C",
        "g": "c",
        "C": "G",
        "c": "g",
        }
        return RNASequence(''.join(transcription_dict[nucleotide] for nucleotide in self.sequence))
    
class RNASequence(NucleicAcidSequence):
    @property
    def ALPHABET(self):
        return set('AUGCaugc')
    
class AminoAcidSequence(BiologicalSequence):
    AMINO_ACID_CATEGORIES = {
    'nonpolar': {'A', 'V', 'L', 'I', 'M', 'F', 'W'},
    'polar': {'S', 'T', 'C', 'Y', 'N', 'Q'},
    '+ charged': {'K', 'R', 'H'},
    '- charged': {'D', 'E'},
    'special': {'P', 'G'},
    'unconventional': {'O', 'U'}
    }
    
    @property
    def ALPHABET(self):
        return set('ACDEFGHIKLMNOPQRSTVWUXYacdefghiklmnopqrstvwuxy')
    
    def __init__(self, sequence):
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError("Sequence contains unknown aminoacid.")
        
    def __str__(self):
        return "Oligoprotein : " + self.sequence
        
    def __getitem__(self, index):
        return self.sequence[index]
    
    def categorize(self):
        counts = {
            'nonpolar': 0,
            'polar': 0,
            '+ charged': 0,
            '- charged': 0,
            'special': 0,
            'unconventional': 0
        }
        for aa in self.sequence.upper():
            for key, value in self.AMINO_ACID_CATEGORIES.items():
                if aa in value:
                    counts[key] += 1
                    
        return counts