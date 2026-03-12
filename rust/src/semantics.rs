use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Semantics {
    #[serde(rename = "legacy-v1")]
    LegacyV1,
    #[serde(rename = "science-v1")]
    ScienceV1,
}

impl Semantics {
    pub fn parse_cli(s: &str) -> Result<Self, String> {
        match s {
            "legacy-v1" => Ok(Self::LegacyV1),
            "science-v1" => Ok(Self::ScienceV1),
            other => Err(format!(
                "invalid --semantics: {other:?} (expected legacy-v1|science-v1)"
            )),
        }
    }

    pub const fn as_str(self) -> &'static str {
        match self {
            Self::LegacyV1 => "legacy-v1",
            Self::ScienceV1 => "science-v1",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum HydrogenPolicy {
    #[serde(rename = "legacy-mmcif-bug")]
    LegacyMmcifBug,
    #[serde(rename = "discard-all")]
    DiscardAll,
    #[serde(rename = "keep-all")]
    KeepAll,
}

impl HydrogenPolicy {
    pub fn parse_cli(s: &str) -> Result<Self, String> {
        match s {
            "legacy-mmcif-bug" => Ok(Self::LegacyMmcifBug),
            "discard-all" => Ok(Self::DiscardAll),
            "keep-all" => Ok(Self::KeepAll),
            other => Err(format!(
                "invalid --hydrogen-policy: {other:?} (expected legacy-mmcif-bug|discard-all|keep-all)"
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MissingInsertionCodePolicy {
    #[serde(rename = "legacy-question-mark")]
    LegacyQuestionMark,
    #[serde(rename = "none")]
    None,
}

impl MissingInsertionCodePolicy {
    pub fn parse_cli(s: &str) -> Result<Self, String> {
        match s {
            "legacy-question-mark" => Ok(Self::LegacyQuestionMark),
            "none" => Ok(Self::None),
            other => Err(format!(
                "invalid --missing-insertion-code: {other:?} (expected legacy-question-mark|none)"
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ChainIdPolicy {
    #[serde(rename = "legacy-1char")]
    Legacy1Char,
    #[serde(rename = "unique-1char")]
    Unique1Char,
}

impl ChainIdPolicy {
    pub fn parse_cli(s: &str) -> Result<Self, String> {
        match s {
            "legacy-1char" => Ok(Self::Legacy1Char),
            "unique-1char" => Ok(Self::Unique1Char),
            other => Err(format!(
                "invalid --chain-id-policy: {other:?} (expected legacy-1char|unique-1char)"
            )),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StructurePolicies {
    pub hydrogen_policy: HydrogenPolicy,
    pub missing_insertion_code_policy: MissingInsertionCodePolicy,
    pub chain_id_policy: ChainIdPolicy,
}

impl StructurePolicies {
    pub const fn legacy_v1_defaults() -> Self {
        Self {
            hydrogen_policy: HydrogenPolicy::LegacyMmcifBug,
            missing_insertion_code_policy: MissingInsertionCodePolicy::LegacyQuestionMark,
            chain_id_policy: ChainIdPolicy::Legacy1Char,
        }
    }

    pub const fn science_v1_defaults() -> Self {
        Self {
            hydrogen_policy: HydrogenPolicy::DiscardAll,
            // Keep non-hydrogen behaviors aligned with legacy-v1 by default, so Gate C focuses on
            // scientifically-motivated changes (starting with the mmCIF hydrogen handling bug).
            // Users can still override via explicit policy flags.
            missing_insertion_code_policy: MissingInsertionCodePolicy::LegacyQuestionMark,
            chain_id_policy: ChainIdPolicy::Legacy1Char,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChemistryPolicies {
    pub residue_identity_policy: String,
    pub edge_catalog_policy: String,
}

impl ChemistryPolicies {
    pub fn science_v1_defaults() -> Self {
        Self {
            residue_identity_policy: "explicit-polymer-v1".to_string(),
            edge_catalog_policy: "dna-aware-lw-v1".to_string(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Policies {
    pub structure: StructurePolicies,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chemistry: Option<ChemistryPolicies>,
}

impl Policies {
    pub const fn legacy_v1_defaults() -> Self {
        Self {
            structure: StructurePolicies::legacy_v1_defaults(),
            chemistry: None,
        }
    }

    pub fn science_v1_defaults() -> Self {
        Self {
            structure: StructurePolicies::science_v1_defaults(),
            chemistry: Some(ChemistryPolicies::science_v1_defaults()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SemanticsConfig {
    pub semantics: Semantics,
    pub policies: Policies,
}

impl SemanticsConfig {
    pub fn defaults(semantics: Semantics) -> Self {
        match semantics {
            Semantics::LegacyV1 => Self {
                semantics,
                policies: Policies::legacy_v1_defaults(),
            },
            Semantics::ScienceV1 => Self {
                semantics,
                policies: Policies::science_v1_defaults(),
            },
        }
    }
}
