from abc import ABC, abstractmethod


class BiologicalSequence(ABC):
    """Abstract base class for biological sequences (DNA, RNA, proteins)."""

    def __init__(self, sequence: str) -> None:
        """
        Initialize the sequence.

        Args:
            sequence (str): The biological sequence.
        """
        self.sequence = sequence

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)

    @property
    @abstractmethod
    def ALPHABET(self) -> set[str]:
        """Return the set of valid symbols for this type of sequence."""
        pass

    def check_alphabet(self) -> bool:
        """
        Check if the sequence contains only valid symbols.

        Returns:
            bool: True if all symbols are in the ALPHABET, False otherwise.
        """
        return set(self.sequence).issubset(self.ALPHABET)

    @abstractmethod
    def __str__(self) -> str:
        """Return a human-readable string representation of the sequence."""
        pass

    @abstractmethod
    def __getitem__(self, index: int | slice) -> "str|BiologicalSequence":
        """
        Allow indexing and slicing of the sequence.

        Args:
            index (int or slice): Position or slice of the sequence.

        Returns:
            str: Single symbol if index is int.
            Same class object: subsequence if index is slice.
        """
        pass


class NucleicAcidSequence(BiologicalSequence):
    """Base class for nucleic acid sequences (DNA/RNA). Cannot be instantiated directly."""

    def __init__(self, sequence: str) -> None:
        if type(self) is NucleicAcidSequence:
            raise NotImplementedError(
                "Direct instantiation of the Nucleic Acid Sequence class is not allowed."
            )
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError(
                f"Sequence contains invalid symbols. Allowed symbols: {self.ALPHABET}"
            )

    def __str__(self) -> str:
        """Return a human-readable string representation of the nucleic acid sequence."""
        return "Oligonucleotide : " + self.sequence

    def __getitem__(self, index: int | slice) -> "str|NucleicAcidSequence":
        """
        Return a single symbol or a subsequence object.

        Args:
            index (int | slice): Position or slice of the sequence.

        Returns:
            str: Single symbol if index is int.
            Same class object: subsequence if index is slice.
        """
        if isinstance(index, int):
            return self.sequence[index]
        return self.__class__(self.sequence[index])

    def complement(self) -> "NucleicAcidSequence":
        """
        Return the complementary sequence.

        Returns:
            NucleicAcidSequence: Complement sequence of the same type as self.
        """
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
        return self.__class__(
            "".join(replication_dict[nucleotide] for nucleotide in self.sequence)
        )

    def reverse(self) -> "NucleicAcidSequence":
        """
        Return the reversed sequence.

        Returns:
            NucleicAcidSequence: Reversed sequence of the same type as self.
        """
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        """
        Return the reverse complement of the sequence.

        Returns:
            NucleicAcidSequence: Reverse complement of the sequence.
        """
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """Class representing a DNA sequence."""

    @property
    def ALPHABET(self) -> set[str]:
        """Return the set of valid DNA nucleotides."""
        return set("ATGCatgc")

    def transcribe(self) -> "RNASequence":
        """
        Transcribe DNA to RNA.

        Returns:
            RNASequence: Transcribed RNA sequence.
        """
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
        return RNASequence(
            "".join(transcription_dict[nucleotide] for nucleotide in self.sequence)
        )


class RNASequence(NucleicAcidSequence):
    """Class representing an RNA sequence."""

    @property
    def ALPHABET(self) -> set[str]:
        """Return the set of valid RNA nucleotides."""
        return set("AUGCaugc")


class AminoAcidSequence(BiologicalSequence):
    """Class representing an amino acid (protein) sequence."""

    AMINO_ACID_CATEGORIES = {
        "nonpolar": {"A", "V", "L", "I", "M", "F", "W"},
        "polar": {"S", "T", "C", "Y", "N", "Q"},
        "+ charged": {"K", "R", "H"},
        "- charged": {"D", "E"},
        "special": {"P", "G"},
        "unconventional": {"O", "U"},
    }

    @property
    def ALPHABET(self) -> set[str]:
        """Return the set of valid amino acid symbols (including lowercase)."""
        return set("ACDEFGHIKLMNOPQRSTVWUXYacdefghiklmnopqrstvwuxy")

    def __init__(self, sequence: str) -> None:
        """
        Initialize amino acid sequence and check validity.

        Args:
            sequence (str): Amino acid sequence.

        Raises:
            ValueError: If sequence contains unknown amino acids.
        """
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError("Sequence contains unknown aminoacid.")

    def __str__(self) -> str:
        """Return a human-readable string representation of the protein."""
        return "Oligoprotein : " + self.sequence

    def __getitem__(self, index: int | slice) -> "str|AminoAcidSequence":
        """
        Return a single symbol or a subsequence object.

        Args:
            index (int | slice): Position or slice of the sequence.

        Returns:
            str: Single symbol if index is int.
            Same class object: subsequence if index is slice.
        """
        if isinstance(index, int):
            return self.sequence[index]
        return self.__class__(self.sequence[index])

    def categorize(self) -> dict[str, int]:
        """
        Count amino acids in each chemical category.

        Returns:
            dict: Dictionary mapping categories to counts.
        """
        counts = {
            "nonpolar": 0,
            "polar": 0,
            "+ charged": 0,
            "- charged": 0,
            "special": 0,
            "unconventional": 0,
        }
        for aa in self.sequence.upper():
            for key, value in self.AMINO_ACID_CATEGORIES.items():
                if aa in value:
                    counts[key] += 1

        return counts
