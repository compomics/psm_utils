"""SQLAlchemy models for Mascot MSF files."""

from sqlalchemy import (
    CHAR,
    BigInteger,
    Boolean,
    Column,
    DateTime,
    Float,
    Index,
    Integer,
    LargeBinary,
    SmallInteger,
    String,
    Table,
    Text,
    UniqueConstraint,
    text,
)

try:
    from sqlalchemy.orm import declarative_base
except ImportError:
    from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql.sqltypes import NullType

Base = declarative_base()
metadata = Base.metadata


class AminoAcidModification(Base):
    __tablename__ = "AminoAcidModifications"

    AminoAcidModificationID = Column(Integer, primary_key=True)
    ModificationName = Column(String, nullable=False)
    DeltaMass = Column(Float)
    Substitution = Column(String)
    LeavingGroup = Column(String)
    Abbreviation = Column(String, nullable=False)
    PositionType = Column(Integer, nullable=False)
    IsActive = Column(Boolean)
    DeltaAverageMass = Column(Float)
    UnimodAccession = Column(String)
    IsSubstitution = Column(Boolean, nullable=False, server_default=text("0"))


class AminoAcidModificationsAminoAcid(Base):
    __tablename__ = "AminoAcidModificationsAminoAcids"

    AminoAcidModificationID = Column(Integer, primary_key=True, nullable=False)
    AminoAcidID = Column(Integer, primary_key=True, nullable=False)
    Classification = Column(Integer, nullable=False)


class AminoAcidModificationsAminoAcidsNL(Base):
    __tablename__ = "AminoAcidModificationsAminoAcidsNL"

    AminoAcidModificationID = Column(Integer, primary_key=True, nullable=False)
    AminoAcidID = Column(Integer, primary_key=True, nullable=False)
    NeutralLossID = Column(Integer, primary_key=True, nullable=False)


class AminoAcidModificationsNeutralLoss(Base):
    __tablename__ = "AminoAcidModificationsNeutralLosses"

    NeutralLossID = Column(Integer, primary_key=True)
    Name = Column(String, nullable=False)
    MonoisotopicMass = Column(Float, nullable=False)
    AverageMass = Column(Float, nullable=False)


class AminoAcid(Base):
    __tablename__ = "AminoAcids"

    AminoAcidID = Column(Integer, primary_key=True)
    AminoAcidName = Column(String, nullable=False)
    OneLetterCode = Column(CHAR)
    ThreeLetterCode = Column(CHAR)
    MonoisotopicMass = Column(Float, nullable=False)
    AverageMass = Column(Float, nullable=False)
    SumFormula = Column(String)


class AnnotationDataVersion(Base):
    __tablename__ = "AnnotationDataVersion"

    PcDataVersion = Column(Integer, primary_key=True)
    PcDataRelease = Column(BigInteger, nullable=False)


class AnnotationDataset(Base):
    __tablename__ = "AnnotationDataset"

    DatasetId = Column(Integer, primary_key=True)
    Name = Column(String, nullable=False)
    DisplayName = Column(String, nullable=False)
    Guid = Column(String, nullable=False)
    Description = Column(Text)


class AnnotationGroup(Base):
    __tablename__ = "AnnotationGroups"

    AnnotationGroupId = Column(Integer, primary_key=True, nullable=False)
    Description = Column(Text)
    DatasetId = Column(Integer, primary_key=True, nullable=False)
    Position = Column(Integer, nullable=False)
    ColorR = Column(Integer, nullable=False)
    ColorG = Column(Integer, nullable=False)
    ColorB = Column(Integer, nullable=False)
    GroupDefinition = Column(LargeBinary)


class AnnotationType(Base):
    __tablename__ = "AnnotationTypes"

    AnnotationTypeId = Column(Integer, primary_key=True)
    Name = Column(String, nullable=False)
    Description = Column(Text)


class Annotation(Base):
    __tablename__ = "Annotations"

    AnnotationId = Column(Integer, primary_key=True)
    Accession = Column(String, nullable=False)
    Description = Column(Text)
    type = Column(Integer)


class AnnotationsAnnotationGroup(Base):
    __tablename__ = "AnnotationsAnnotationGroups"

    AnnotationId = Column(Integer, primary_key=True, nullable=False)
    AnnotationGroupId = Column(Integer, primary_key=True, nullable=False)


class AnnotationsProtein(Base):
    __tablename__ = "AnnotationsProtein"

    proteinID = Column(Integer, primary_key=True, nullable=False)
    AnnotationId = Column(Integer, primary_key=True, nullable=False)
    Evidence = Column(Integer, primary_key=True)
    PositionBegin = Column(Integer, primary_key=True)
    PositionEnd = Column(Integer)
    ProteinAccession = Column(String, primary_key=True, nullable=False)


class Chromatogram(Base):
    __tablename__ = "Chromatograms"

    FileID = Column(Integer, primary_key=True, nullable=False)
    TraceType = Column(Integer, primary_key=True, nullable=False)
    Chromatogram = Column(String, nullable=False)


class CustomDataField(Base):
    __tablename__ = "CustomDataFields"

    FieldID = Column(Integer, primary_key=True)
    Guid = Column(String, nullable=False)
    DisplayName = Column(String, nullable=False)
    SourceNodeNumber = Column(Integer, nullable=False)
    TargetNodeNumber = Column(Integer, nullable=False)
    DataType = Column(Integer, nullable=False)
    DataTarget = Column(Integer, nullable=False)
    Version = Column(Float, nullable=False)
    AccessMode = Column(Integer, server_default=text("0"))
    Visibility = Column(Integer, server_default=text("0"))
    GroupVisibility = Column(Integer, server_default=text("0"))
    Format = Column(String)
    PlotType = Column(Integer, nullable=False)
    DataPurpose = Column(String)


class CustomDataPeptide(Base):
    __tablename__ = "CustomDataPeptides"

    FieldID = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue = Column(String)


class CustomDataPeptidesDecoy(Base):
    __tablename__ = "CustomDataPeptides_decoy"

    FieldID = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue = Column(String)


class CustomDataProcessingNode(Base):
    __tablename__ = "CustomDataProcessingNodes"

    FieldID = Column(Integer, primary_key=True, nullable=False)
    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue = Column(String)


class CustomDataProtein(Base):
    __tablename__ = "CustomDataProteins"

    FieldID = Column(Integer, primary_key=True, nullable=False)
    ProteinID = Column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue = Column(String)


class CustomDataProteinsDecoy(Base):
    __tablename__ = "CustomDataProteins_decoy"

    FieldID = Column(Integer, primary_key=True, nullable=False)
    ProteinID = Column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue = Column(String)


class CustomDataSpectra(Base):
    __tablename__ = "CustomDataSpectra"

    FieldID = Column(Integer, primary_key=True, nullable=False)
    SpectrumID = Column(Integer, primary_key=True, nullable=False, index=True)
    FieldValue = Column(String)


class Enzyme(Base):
    __tablename__ = "Enzymes"

    EnzymeID = Column(Integer, primary_key=True)
    Name = Column(String, nullable=False)
    Abbreviation = Column(String, nullable=False)
    Seperator = Column(String, nullable=False)
    NonSeperator = Column(String, nullable=False)
    Offset = Column(Integer, nullable=False)


class EnzymesCleavageSpecificity(Base):
    __tablename__ = "EnzymesCleavageSpecificities"

    EnzymeID = Column(Integer, primary_key=True, nullable=False)
    Specificity = Column(Integer, primary_key=True, nullable=False)


class EventAnnotation(Base):
    __tablename__ = "EventAnnotations"
    __table_args__ = (
        Index(
            "IX_EventAnnotations_IsotopePatternID_QuanResultID", "IsotopePatternID", "QuanResultID"
        ),
        Index("IX_EventAnnotations_QuanResultID_QuanChannelID", "QuanResultID", "QuanChannelID"),
    )

    EventID = Column(Integer, primary_key=True)
    Charge = Column(SmallInteger, nullable=False)
    IsotopePatternID = Column(Integer, nullable=False)
    QuanResultID = Column(Integer, nullable=False)
    QuanChannelID = Column(Integer, nullable=False)


class EventAreaAnnotation(Base):
    __tablename__ = "EventAreaAnnotations"

    EventID = Column(Integer, primary_key=True)
    Charge = Column(SmallInteger, nullable=False)
    IsotopePatternID = Column(Integer, nullable=False, index=True)
    QuanResultID = Column(Integer, nullable=False)


class Event(Base):
    __tablename__ = "Events"
    __table_args__ = (
        Index("IX_Events_FileID_LeftRT_RightRT", "FileID", "LeftRT", "RightRT"),
        Index("IX_Events_FileID_RT", "FileID", "RT"),
    )

    EventID = Column(Integer, primary_key=True)
    Mass = Column(Float, nullable=False)
    MassAvg = Column(Float, nullable=False)
    Area = Column(Float, nullable=False)
    Intensity = Column(Float, nullable=False)
    PeakWidth = Column(Float, nullable=False)
    RT = Column(Float, nullable=False)
    LeftRT = Column(Float, nullable=False)
    RightRT = Column(Float, nullable=False)
    SN = Column(Float, nullable=False, server_default=text("0.0"))
    FileID = Column(Integer, nullable=False)


class FastaFile(Base):
    __tablename__ = "FastaFiles"

    FastaFileID = Column(Integer, primary_key=True)
    FileName = Column(String, nullable=False)
    State = Column(Integer, nullable=False)
    VirtualFileName = Column(String, nullable=False)
    FileSize = Column(BigInteger, nullable=False)
    FileTime = Column(BigInteger, nullable=False)
    NumberOfProteins = Column(BigInteger)
    NumberOfAminoAcids = Column(BigInteger)
    FileHashCode = Column(BigInteger)
    Hidden = Column(Boolean, nullable=False)
    IsSrfImport = Column(Boolean, nullable=False)
    IsScheduledForDeletion = Column(Boolean, nullable=False, server_default=text("0"))


class FastaFilesProteinAnnotation(Base):
    __tablename__ = "FastaFilesProteinAnnotations"

    FastaFileID = Column(Integer, primary_key=True, nullable=False)
    ProteinAnnotationID = Column(Integer, primary_key=True, nullable=False, index=True)


class FileInfo(Base):
    __tablename__ = "FileInfos"

    FileID = Column(Integer, primary_key=True)
    FileName = Column(String, nullable=False)
    FileTime = Column(String, nullable=False)
    FileSize = Column(BigInteger, nullable=False)
    PhysicalFileName = Column(String, nullable=False)
    FileType = Column(SmallInteger, nullable=False)


class MassPeakRelation(Base):
    __tablename__ = "MassPeakRelations"

    MassPeakID = Column(Integer, primary_key=True, nullable=False)
    RelatedMassPeakID = Column(Integer, primary_key=True, nullable=False)


class MassPeak(Base):
    __tablename__ = "MassPeaks"

    MassPeakID = Column(Integer, primary_key=True)
    Charge = Column(SmallInteger)
    Intensity = Column(Float)
    Mass = Column(Float)
    ScanNumbers = Column(String)
    FileID = Column(Integer)
    PercentIsolationInterference = Column(Float)
    IonInjectTime = Column(Integer)


class PeptideScore(Base):
    __tablename__ = "PeptideScores"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False)
    ScoreID = Column(Integer, primary_key=True, nullable=False)
    ProcessingNodeID = Column(Integer)
    ScoreValue = Column(Float, nullable=False)


class PeptideScoreDecoy(Base):
    __tablename__ = "PeptideScores_decoy"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False)
    ScoreID = Column(Integer, primary_key=True, nullable=False)
    ProcessingNodeID = Column(Integer)
    ScoreValue = Column(Float, nullable=False)


class Peptide(Base):
    __tablename__ = "Peptides"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False, index=True)
    SpectrumID = Column(Integer, nullable=False, index=True)
    TotalIonsCount = Column(SmallInteger, nullable=False)
    MatchedIonsCount = Column(SmallInteger, nullable=False)
    ConfidenceLevel = Column(SmallInteger, nullable=False)
    SearchEngineRank = Column(Integer, nullable=False)
    Hidden = Column(Boolean, nullable=False, server_default=text("0"))
    Sequence = Column(String)
    Annotation = Column(String)
    UniquePeptideSequenceID = Column(Integer, nullable=False, server_default=text("1"))
    MissedCleavages = Column(SmallInteger, nullable=False)


class PeptidesAminoAcidModification(Base):
    __tablename__ = "PeptidesAminoAcidModifications"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False)
    AminoAcidModificationID = Column(Integer, primary_key=True, nullable=False)
    Position = Column(Integer, primary_key=True, nullable=False)


class PeptidesAminoAcidModificationsDecoy(Base):
    __tablename__ = "PeptidesAminoAcidModifications_decoy"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False)
    AminoAcidModificationID = Column(Integer, primary_key=True, nullable=False)
    Position = Column(Integer, primary_key=True, nullable=False)


class PeptidesProtein(Base):
    __tablename__ = "PeptidesProteins"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False, index=True)
    ProteinID = Column(Integer, primary_key=True, nullable=False)


class PeptidesProteinDecoy(Base):
    __tablename__ = "PeptidesProteins_decoy"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False, index=True)
    ProteinID = Column(Integer, primary_key=True, nullable=False)


class PeptidesReferenceSpectra(Base):
    __tablename__ = "PeptidesReferenceSpectra"

    PeptideID = Column(Integer, primary_key=True)
    ReferenceSpectrumID = Column(Integer, nullable=False)


class PeptidesTerminalModification(Base):
    __tablename__ = "PeptidesTerminalModifications"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False)
    TerminalModificationID = Column(Integer, primary_key=True, nullable=False)


class PeptidesTerminalModificationDecoy(Base):
    __tablename__ = "PeptidesTerminalModifications_decoy"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False)
    TerminalModificationID = Column(Integer, primary_key=True, nullable=False)


class PeptideDecoy(Base):
    __tablename__ = "Peptides_decoy"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    PeptideID = Column(Integer, primary_key=True, nullable=False, index=True)
    SpectrumID = Column(Integer, nullable=False, index=True)
    TotalIonsCount = Column(SmallInteger, nullable=False)
    MatchedIonsCount = Column(SmallInteger, nullable=False)
    ConfidenceLevel = Column(SmallInteger, nullable=False)
    SearchEngineRank = Column(Integer, nullable=False)
    Sequence = Column(String)
    Annotation = Column(String)
    UniquePeptideSequenceID = Column(Integer, nullable=False, server_default=text("1"))
    MissedCleavages = Column(SmallInteger, nullable=False)


t_PrecursorIonAreaSearchSpectra = Table(
    "PrecursorIonAreaSearchSpectra",
    metadata,
    Column("QuanResultID", Integer, nullable=False, index=True),
    Column("SearchSpectrumID", Integer),
)


t_PrecursorIonQuanResults = Table(
    "PrecursorIonQuanResults",
    metadata,
    Column("QuanChannelID", Integer, nullable=False),
    Column("QuanResultID", Integer, nullable=False),
    Column("Mass", Float, nullable=False),
    Column("Charge", Integer, nullable=False),
    Column("Area", Float),
    Column("RetentionTime", Float),
    Index(
        "IX_PrecursorIonQuanResults_QuanResultID_QuanChannelID", "QuanResultID", "QuanChannelID"
    ),
)


t_PrecursorIonQuanResultsSearchSpectra = Table(
    "PrecursorIonQuanResultsSearchSpectra",
    metadata,
    Column("ProcessingNodeNumber", Integer, nullable=False),
    Column("QuanResultID", Integer, nullable=False, index=True),
    Column("SearchSpectrumID", Integer, index=True),
)


t_ProcessingNodeConnectionPoints = Table(
    "ProcessingNodeConnectionPoints",
    metadata,
    Column("ProcessingNodeID", Integer, nullable=False),
    Column("Interface", String, nullable=False),
    Column("ConnectionDirection", Integer, nullable=False),
    Column("ConnectionMode", Integer, nullable=False),
    Column("ConnectionMultiplicity", Integer, nullable=False),
    Column("ConnectionRequirement", Integer, nullable=False),
    Column("DataTypeSpecialization", String, nullable=False),
    Column("ConnectionDisplayName", String, nullable=False),
)


class ProcessingNodeExtension(Base):
    __tablename__ = "ProcessingNodeExtensions"

    ExtensionID = Column(Integer, primary_key=True)
    ProcessingNodeNumber = Column(Integer, nullable=False)
    Guid = Column(String, nullable=False)
    Purpose = Column(String, nullable=False)
    PurposeDetail = Column(String)
    MajorVersion = Column(Integer, nullable=False)
    MinorVersion = Column(Integer, nullable=False)
    Settings = Column(Text)


class ProcessingNodeFilterParameter(Base):
    __tablename__ = "ProcessingNodeFilterParameters"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    FilterParameterName = Column(String, primary_key=True, nullable=False)
    FilterModuleTypeID = Column(Integer, nullable=False)
    FilterModuleNumber = Column(Integer, nullable=False)
    ProcessingNodeID = Column(Integer, nullable=False)
    FilterParameterValue = Column(Float, nullable=False)


t_ProcessingNodeInterfaces = Table(
    "ProcessingNodeInterfaces",
    metadata,
    Column("ProcessingNodeID", Integer, nullable=False),
    Column("InterfaceKind", Integer, nullable=False),
    Column("InterfaceName", String, nullable=False),
)


class ProcessingNodeParameter(Base):
    __tablename__ = "ProcessingNodeParameters"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    ParameterName = Column(String, primary_key=True, nullable=False)
    FriendlyName = Column(String, nullable=False)
    ProcessingNodeID = Column(Integer, nullable=False)
    IntendedPurpose = Column(Integer, nullable=False)
    PurposeDetails = Column(String, nullable=False)
    Hidden = Column(Boolean, nullable=False)
    Advanced = Column(Boolean, nullable=False)
    Category = Column(String, nullable=False)
    Position = Column(Integer, nullable=False)
    ParameterValue = Column(String, nullable=False)
    ValueDisplayString = Column(String, nullable=False)


class ProcessingNodeScore(Base):
    __tablename__ = "ProcessingNodeScores"
    __table_args__ = (UniqueConstraint("ProcessingNodeID", "ScoreName"),)

    ProcessingNodeID = Column(Integer, nullable=False)
    ScoreID = Column(Integer, primary_key=True)
    ScoreName = Column(String, nullable=False)
    FriendlyName = Column(String, nullable=False)
    Description = Column(String, nullable=False)
    FormatString = Column(String, nullable=False)
    ScoreCategory = Column(Integer, nullable=False)
    Hidden = Column(Boolean, nullable=False)
    IsMainScore = Column(Boolean, nullable=False)
    ScoreGUID = Column(String, nullable=False)


class ProcessingNode(Base):
    __tablename__ = "ProcessingNodes"

    ProcessingNodeNumber = Column(Integer, primary_key=True)
    ProcessingNodeID = Column(Integer, nullable=False)
    ProcessingNodeParentNumber = Column(String, nullable=False)
    NodeName = Column(String)
    FriendlyName = Column(String, nullable=False)
    MajorVersion = Column(Integer, nullable=False)
    MinorVersion = Column(Integer, nullable=False)
    NodeComment = Column(String)
    NodeGUID = Column(String, nullable=False)
    ProcessingNodeState = Column(Integer, nullable=False, server_default=text("0"))


class ProcessingNodesSpectra(Base):
    __tablename__ = "ProcessingNodesSpectra"

    SendingProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    SpectrumID = Column(Integer, primary_key=True, nullable=False, index=True)


class ProteinAnnotation(Base):
    __tablename__ = "ProteinAnnotations"
    __table_args__ = (
        Index(
            "IX_ProteinAnnotations_ProteinID_DescriptionHashCode",
            "ProteinID",
            "DescriptionHashCode",
        ),
    )

    ProteinAnnotationID = Column(Integer, primary_key=True)
    ProteinID = Column(Integer, nullable=False)
    DescriptionHashCode = Column(BigInteger, nullable=False)
    Description = Column(Text, nullable=False)
    TaxonomyID = Column(Integer, nullable=False, index=True)


class ProteinIdentificationGroup(Base):
    __tablename__ = "ProteinIdentificationGroups"

    ProteinIdentificationGroupId = Column(Integer, primary_key=True, nullable=False)
    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)


class ProteinScore(Base):
    __tablename__ = "ProteinScores"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    ProteinID = Column(Integer, primary_key=True, nullable=False)
    ProteinIdentificationGroupID = Column(Integer, nullable=False)
    ProteinScore = Column(Float, nullable=False)
    Coverage = Column(Float, nullable=False, server_default=text("0"))


class ProteinScoresDecoy(Base):
    __tablename__ = "ProteinScores_decoy"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    ProteinID = Column(Integer, primary_key=True, nullable=False)
    ProteinIdentificationGroupID = Column(Integer, nullable=False)
    ProteinScore = Column(Float, nullable=False)
    Coverage = Column(Float, nullable=False, server_default=text("0"))


class Protein(Base):
    __tablename__ = "Proteins"

    ProteinID = Column(Integer, primary_key=True)
    Sequence = Column(Text, nullable=False)
    SequenceHashCode = Column(BigInteger, nullable=False, index=True)
    IsMasterProtein = Column(Boolean, nullable=False, server_default=text("0"))


t_ProteinsProteinGroups = Table(
    "ProteinsProteinGroups",
    metadata,
    Column("ProteinID", Integer, nullable=False),
    Column("ProteinGroupID", Integer, nullable=False),
)


class PtmAnnotationDatum(Base):
    __tablename__ = "PtmAnnotationData"

    AnnotationType = Column(Integer, primary_key=True, nullable=False)
    ProteinId = Column(Integer, primary_key=True, nullable=False)
    AnnotationId = Column(Integer, primary_key=True, nullable=False)
    Position = Column(Integer, primary_key=True, nullable=False)
    Annotation = Column(String)


class ReferenceSpectra(Base):
    __tablename__ = "ReferenceSpectra"

    ReferenceSpectrumId = Column(Integer, primary_key=True)
    Sequence = Column(String, nullable=False)
    SequenceHashCode = Column(BigInteger, nullable=False)
    Spectrum = Column(String, nullable=False)
    SpectrumHashCode = Column(BigInteger, nullable=False)
    Comment = Column(Text)
    CommentHashCode = Column(BigInteger, nullable=False)


class ReporterIonQuanResult(Base):
    __tablename__ = "ReporterIonQuanResults"
    __table_args__ = (
        Index(
            "IX_ReporterIonQuanResults_ProcessingNodeNumber_SpectrumID",
            "ProcessingNodeNumber",
            "SpectrumID",
        ),
    )

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    QuanChannelID = Column(Integer, primary_key=True, nullable=False)
    SpectrumID = Column(Integer, primary_key=True, nullable=False)
    Mass = Column(Float, nullable=False)
    Height = Column(Float)


t_ReporterIonQuanResultsSearchSpectra = Table(
    "ReporterIonQuanResultsSearchSpectra",
    metadata,
    Column("ProcessingNodeNumber", Integer, nullable=False),
    Column("SpectrumID", Integer, nullable=False),
    Column("SearchSpectrumID", Integer, index=True),
)


class ScanEvent(Base):
    __tablename__ = "ScanEvents"

    ScanEventID = Column(Integer, primary_key=True)
    MSLevel = Column(Integer, nullable=False)
    Polarity = Column(Integer, nullable=False)
    ScanType = Column(Integer, nullable=False)
    Ionization = Column(Integer, nullable=False)
    MassAnalyzer = Column(Integer, nullable=False)
    ActivationType = Column(Integer, nullable=False)


class SchemaInfo(Base):
    __tablename__ = "SchemaInfo"

    Version = Column(Integer, primary_key=True)
    Kind = Column(String, nullable=False)
    Date = Column(DateTime, nullable=False)
    SoftwareVersion = Column(String, nullable=False)
    Comment = Column(Text, nullable=False)


class Spectrum(Base):
    __tablename__ = "Spectra"

    UniqueSpectrumID = Column(Integer, primary_key=True)
    Spectrum = Column(String, nullable=False)
    SpectrumHashCode = Column(BigInteger)


class SpectrumHeader(Base):
    __tablename__ = "SpectrumHeaders"

    SpectrumID = Column(Integer, primary_key=True)
    MassPeakID = Column(Integer)
    ScanEventID = Column(Integer)
    LastScan = Column(Integer)
    FirstScan = Column(Integer)
    RetentionTime = Column(Float)
    Hidden = Column(Boolean, nullable=False, server_default=text("0"))
    ScanNumbers = Column(String)
    Charge = Column(SmallInteger)
    Mass = Column(Float)
    CreatingProcessingNodeNumber = Column(Integer, nullable=False)
    UniqueSpectrumID = Column(Integer, nullable=False, server_default=text("0"))


class SpectrumScore(Base):
    __tablename__ = "SpectrumScores"

    ProcessingNodeNumber = Column(Integer, primary_key=True, nullable=False)
    SpectrumID = Column(Integer, primary_key=True, nullable=False)
    Score = Column(Float, nullable=False)


t_TaxonomyNames = Table(
    "TaxonomyNames",
    metadata,
    Column("TaxonomyID", Integer, nullable=False, index=True),
    Column("Name", String),
    Column("NameCategory", Integer, nullable=False),
)


class TaxonomyNode(Base):
    __tablename__ = "TaxonomyNodes"
    __table_args__ = (
        Index("IX_TaxonomyNodes_LeftNodeIndex_RightNodeIndex", "LeftNodeIndex", "RightNodeIndex"),
    )

    TaxonomyID = Column(Integer, primary_key=True, unique=True)
    ParentTaxonomyID = Column(Integer, nullable=False)
    TaxonomyRank = Column(Integer, nullable=False)
    LeftNodeIndex = Column(Integer, nullable=False)
    RightNodeIndex = Column(Integer, nullable=False)


t_WorkflowInfo = Table(
    "WorkflowInfo",
    metadata,
    Column("WorkflowName", String, nullable=False),
    Column("WorkflowDescription", String, nullable=False),
    Column("WorkflowState", Integer, nullable=False, server_default=text("0")),
    Column("WorkflowStartDate", DateTime, nullable=False),
    Column("WorkflowTemplate", String, nullable=False),
    Column("User", String, nullable=False),
    Column("WorkflowGUID", String, nullable=False),
    Column("MachineGUID", String, nullable=False),
    Column("MachineName", String, nullable=False),
    Column("MergeSimilarIdentificationResults", Boolean, nullable=False),
    Column("IsValid", Boolean, nullable=False),
    Column("Version", Integer, nullable=False),
)


class WorkflowMessage(Base):
    __tablename__ = "WorkflowMessages"

    MessageID = Column(Integer, primary_key=True)
    ProcessingNodeID = Column(Integer, nullable=False)
    ProcessingNodeNumber = Column(Integer, nullable=False)
    Time = Column(BigInteger, nullable=False)
    MessageKind = Column(Integer, nullable=False)
    Message = Column(String, nullable=False)


t_sqlite_sequence = Table(
    "sqlite_sequence", metadata, Column("name", NullType), Column("seq", NullType)
)
