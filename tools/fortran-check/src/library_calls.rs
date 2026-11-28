use std::collections::{BTreeMap, BTreeSet};

use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};
use tree_sitter::Node;

use crate::{
    config::{Config, matches},
    flow::{Edge, unwrap},
    semantic::{Project, Scope, Symbol},
    source::{Unit, child, children},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Family {
    Mpi,
    Petsc,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArgumentSelector {
    #[serde(default = "enabled")]
    pub enabled: bool,
    pub index: i32,
    #[serde(default)]
    pub keywords: Vec<String>,
}

impl ArgumentSelector {
    pub fn new(index: i32, keywords: &[&str]) -> Self {
        Self {
            enabled: true,
            index,
            keywords: strings(keywords),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResourcePair {
    #[serde(default = "enabled")]
    pub enabled: bool,
    pub family: Family,
    pub release: String,
    pub acquire_handle: ArgumentSelector,
    pub release_handle: ArgumentSelector,
    pub acquire_status: Option<ArgumentSelector>,
    pub release_status: Option<ArgumentSelector>,
    #[serde(default)]
    pub acquire_context: Vec<ArgumentSelector>,
    #[serde(default)]
    pub release_context: Vec<ArgumentSelector>,
    #[serde(default)]
    pub borrowed: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct MpiSignature {
    pub enabled: bool,
    pub status: Option<ArgumentSelector>,
    pub outputs: Vec<ArgumentSelector>,
    pub rank: Option<RankSignature>,
    pub request: Option<RequestSignature>,
    pub collective: Option<CollectiveSignature>,
}

impl Default for MpiSignature {
    fn default() -> Self {
        Self {
            enabled: true,
            status: None,
            outputs: Vec::new(),
            rank: None,
            request: None,
            collective: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RankSignature {
    #[serde(default = "enabled")]
    pub enabled: bool,
    pub rank: ArgumentSelector,
    pub communicator: ArgumentSelector,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RequestAction {
    Post,
    Initialize,
    Start,
    Wait,
    Test,
    Cancel,
    Free,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RequestSignature {
    #[serde(default = "enabled")]
    pub enabled: bool,
    pub action: RequestAction,
    pub handle: ArgumentSelector,
    pub count: Option<ArgumentSelector>,
    pub completion_flag: Option<ArgumentSelector>,
    #[serde(default)]
    pub collective: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CollectiveSignature {
    #[serde(default = "enabled")]
    pub enabled: bool,
    pub communicator: ArgumentSelector,
    #[serde(default)]
    pub compare_arguments: Vec<ArgumentSelector>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SplitPair {
    #[serde(default = "enabled")]
    pub enabled: bool,
    pub end: String,
    #[serde(default)]
    pub context: Vec<ArgumentSelector>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SafetyPolicy {
    pub mpi_modules: Vec<String>,
    pub petsc_modules: Vec<String>,
    pub mpi_prefixes: Vec<String>,
    pub petsc_prefixes: Vec<String>,
    pub call_families: BTreeMap<String, Family>,
    pub resources: BTreeMap<String, ResourcePair>,
    pub mpi_calls: BTreeMap<String, MpiSignature>,
    pub petsc_status: ArgumentSelector,
    pub petsc_call_status: BTreeMap<String, ArgumentSelector>,
    pub petsc_no_status: Vec<String>,
    pub petsc_splits: BTreeMap<String, SplitPair>,
    pub status_handlers: BTreeMap<String, ArgumentSelector>,
    pub terminators: Vec<String>,
    pub finalizers: Vec<String>,
    pub success_constants: Vec<String>,
    pub self_communicators: Vec<String>,
    pub request_nulls: Vec<String>,
    pub pointer_outputs: BTreeMap<String, ArgumentSelector>,
    pub collective_path_limit: usize,
    pub check_pure: bool,
    pub check_do_concurrent: bool,
    pub check_openmp: bool,
}

fn enabled() -> bool {
    true
}

fn strings(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| (*item).to_owned()).collect()
}

impl Default for SafetyPolicy {
    fn default() -> Self {
        let resources = resources();
        let mut mpi_calls = mpi_signatures();
        for (name, pair) in &resources {
            if pair.family == Family::Mpi {
                mpi_calls
                    .entry(name.clone())
                    .or_insert_with(|| MpiSignature {
                        status: pair.acquire_status.clone(),
                        outputs: vec![pair.acquire_handle.clone()],
                        ..MpiSignature::default()
                    });
                mpi_calls
                    .entry(pair.release.clone())
                    .or_insert_with(|| MpiSignature {
                        status: pair.release_status.clone(),
                        ..MpiSignature::default()
                    });
            }
        }
        let petsc_splits = [
            "dmglobaltolocal",
            "dmlocaltoglobal",
            "dmlocaltolocal",
            "vecscatter",
            "vecassembly",
            "matassembly",
            "vecghostupdate",
            "petscsfbcast",
            "petscsfreduce",
        ]
        .into_iter()
        .map(|base| {
            (
                format!("{base}begin"),
                SplitPair {
                    enabled: true,
                    end: format!("{base}end"),
                    context: Vec::new(),
                },
            )
        })
        .collect();
        let status_handlers = [
            ("mpi_abort", 1, "errorcode"),
            ("mpiu_abort", 1, "errorcode"),
            ("mpi_error_class", 0, "errorcode"),
            ("mpi_error_string", 0, "errorcode"),
            ("petscerrorf", 0, "ierr"),
            ("petscerror", 0, "ierr"),
        ]
        .into_iter()
        .map(|(name, index, key)| (name.to_owned(), ArgumentSelector::new(index, &[key])))
        .collect();
        let petsc_call_status = [
            ("vecsetfromoptions", 1),
            ("vecassemblybegin", 1),
            ("vecassemblyend", 1),
            ("matassemblybegin", 2),
            ("matassemblyend", 2),
            ("vecscatterbegin", 5),
            ("vecscatterend", 5),
        ]
        .into_iter()
        .map(|(name, index)| {
            (
                name.to_owned(),
                ArgumentSelector::new(index, &["ierr", "ierror"]),
            )
        })
        .collect();
        Self {
            mpi_modules: strings(&["mpi*"]),
            petsc_modules: strings(&["petsc*"]),
            mpi_prefixes: strings(&["mpi_", "mpiu_"]),
            petsc_prefixes: strings(&[
                "petsc", "ao", "dm", "draw", "is", "ksp", "mat", "pc", "pf", "snes", "tao", "ts",
                "vec",
            ]),
            call_families: BTreeMap::new(),
            resources,
            mpi_calls,
            petsc_status: ArgumentSelector::new(-1, &["ierr", "ierror"]),
            petsc_call_status,
            petsc_no_status: strings(&["petscerror", "petscerrorf", "petscobjectnullify"]),
            petsc_splits,
            status_handlers,
            terminators: strings(&["mpi_abort", "mpiu_abort"]),
            finalizers: strings(&["mpi_finalize", "petscfinalize"]),
            success_constants: strings(&["mpi_success", "petsc_success"]),
            self_communicators: strings(&["mpi_comm_self", "petsc_comm_self"]),
            request_nulls: strings(&["mpi_request_null"]),
            pointer_outputs: BTreeMap::from([(
                "c_f_pointer".into(),
                ArgumentSelector::new(1, &["fptr"]),
            )]),
            collective_path_limit: 128,
            check_pure: true,
            check_do_concurrent: true,
            check_openmp: true,
        }
    }
}

impl SafetyPolicy {
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.collective_path_limit > 0,
            "safety.collective_path_limit must be greater than zero"
        );
        let mut releases = BTreeMap::new();
        for (name, pair) in &self.resources {
            ensure!(
                !name.is_empty() && !pair.release.is_empty(),
                "safety.resources requires named acquisition and release calls"
            );
            ensure!(
                pair.acquire_context.len() == pair.release_context.len(),
                "safety.resources.{name} must select the same number of acquisition and release context arguments"
            );
            if pair.enabled {
                let signature = (
                    &pair.release_handle,
                    &pair.release_status,
                    &pair.release_context,
                    pair.borrowed,
                );
                if let Some(previous) = releases.insert((pair.family, &pair.release), signature) {
                    ensure!(
                        signature == previous,
                        "safety.resources.{name} conflicts with another release signature for {}",
                        pair.release
                    );
                }
            }
        }
        let mut endings = BTreeSet::new();
        for (name, pair) in &self.petsc_splits {
            ensure!(
                !name.is_empty() && !pair.end.is_empty() && *name != pair.end,
                "safety.petsc_splits requires distinct begin and end calls"
            );
            ensure!(
                !pair.enabled || endings.insert(&pair.end),
                "safety.petsc_splits.{name} duplicates an end call"
            );
        }
        Ok(())
    }

    pub fn mpi(&self, call: &Call<'_>) -> Option<&MpiSignature> {
        if call.family != Family::Mpi {
            return None;
        }
        self.mpi_calls
            .get(&call.name)
            .filter(|signature| signature.enabled)
    }

    pub fn expected_status(&self, call: &Call<'_>) -> Option<&ArgumentSelector> {
        let selector = match call.family {
            Family::Mpi => self.mpi(call)?.status.as_ref(),
            Family::Petsc if self.petsc_no_status.iter().any(|name| name == &call.name) => None,
            Family::Petsc => Some(
                self.petsc_call_status
                    .get(&call.name)
                    .unwrap_or(&self.petsc_status),
            ),
        };
        selector.filter(|selector| selector.enabled)
    }

    pub fn terminates(&self, call: Option<&Call<'_>>) -> bool {
        call.is_some_and(|call| self.terminators.contains(&call.name))
    }
}

fn resources() -> BTreeMap<String, ResourcePair> {
    let mut pairs = BTreeMap::new();
    let petsc_created = [
        (
            "vecdestroy",
            "veccreate veccreateseq veccreatempi veccreatenest veccreateseqwitharray veccreatempiwitharray vecduplicate dmcreateglobalvector dmcreatelocalvector",
        ),
        (
            "matdestroy",
            "matcreate matcreateaij matcreateseqaij matcreatempiaij matcreatedense matduplicate dmcreatematrix matcreatesnesmf",
        ),
        (
            "dmdestroy",
            "dmcreate dmclone dmplexcreate dmplexcreatefromfile dmplexcreatefromcelllistpetsc dmplexcreateboxmesh",
        ),
        (
            "isdestroy",
            "iscreategeneral iscreatestride iscreateblock isduplicate dmgetlabelidis dmlabelgetstratumis",
        ),
        (
            "petscviewerdestroy",
            "petscviewercreate petscviewerasciiopen petscviewerbinaryopen petscviewerhdf5open petscviewervtkopen petscviewerdrawopen",
        ),
        ("petscsectiondestroy", "petscsectioncreate"),
        (
            "petscfedestroy",
            "petscfecreatedefault petscfecreatelagrange",
        ),
        ("petscdsdestroy", "petscdscreate"),
        ("petscsfdestroy", "petscsfcreate"),
        ("petscquadraturedestroy", "petscquadraturecreate"),
        ("snesdestroy", "snescreate"),
        ("kspdestroy", "kspcreate"),
        ("pcdestroy", "pccreate"),
        ("tsdestroy", "tscreate"),
    ];
    for (release, acquisitions) in petsc_created {
        for acquire in acquisitions.split_whitespace() {
            pairs.insert(
                acquire.into(),
                ResourcePair {
                    enabled: true,
                    family: Family::Petsc,
                    release: release.into(),
                    borrowed: false,
                    acquire_handle: ArgumentSelector::new(-2, &[]),
                    release_handle: ArgumentSelector::new(0, &[]),
                    acquire_status: Some(ArgumentSelector::new(-1, &["ierr"])),
                    release_status: Some(ArgumentSelector::new(-1, &["ierr"])),
                    acquire_context: Vec::new(),
                    release_context: Vec::new(),
                },
            );
        }
    }
    let petsc_borrowed: [(&str, &str, i32, &[i32]); 11] = [
        ("vecgetarray", "vecrestorearray", 1, &[0]),
        ("vecgetarrayread", "vecrestorearrayread", 1, &[0]),
        ("vecgetarraywrite", "vecrestorearraywrite", 1, &[0]),
        ("dmgetlocalvector", "dmrestorelocalvector", 1, &[0]),
        ("dmgetglobalvector", "dmrestoreglobalvector", 1, &[0]),
        ("isgetindices", "isrestoreindices", 1, &[0]),
        ("dmplexgetcone", "dmplexrestorecone", 2, &[0, 1]),
        ("dmplexgetsupport", "dmplexrestoresupport", 2, &[0, 1]),
        ("dmplexgetjoin", "dmplexrestorejoin", 4, &[0, 1, 2]),
        (
            "dmplexgettransitiveclosure",
            "dmplexrestoretransitiveclosure",
            4,
            &[0, 1, 2],
        ),
        (
            "dmplexvecgetclosure",
            "dmplexvecrestoreclosure",
            5,
            &[0, 1, 2, 3],
        ),
    ];
    for (acquire, release, handle, contexts) in petsc_borrowed {
        let context: Vec<_> = contexts
            .iter()
            .map(|index| ArgumentSelector::new(*index, &[]))
            .collect();
        pairs.insert(
            acquire.into(),
            ResourcePair {
                enabled: true,
                family: Family::Petsc,
                release: release.into(),
                borrowed: true,
                acquire_handle: ArgumentSelector::new(handle, &[]),
                release_handle: ArgumentSelector::new(handle, &[]),
                acquire_status: Some(ArgumentSelector::new(-1, &["ierr"])),
                release_status: Some(ArgumentSelector::new(-1, &["ierr"])),
                acquire_context: context.clone(),
                release_context: context,
            },
        );
    }
    type MpiResourceGroup = (
        &'static str,
        &'static str,
        &'static str,
        &'static [(&'static str, i32)],
    );
    let mpi_created: [MpiResourceGroup; 8] = [
        (
            "comm_free",
            "comm",
            "newcomm",
            &[
                ("comm_dup", 1),
                ("comm_dup_with_info", 2),
                ("comm_create", 2),
                ("comm_create_group", 3),
                ("comm_split", 3),
                ("comm_split_type", 4),
            ],
        ),
        (
            "group_free",
            "group",
            "newgroup",
            &[
                ("comm_group", 1),
                ("comm_remote_group", 1),
                ("group_incl", 3),
                ("group_excl", 3),
                ("group_range_incl", 3),
                ("group_range_excl", 3),
                ("group_union", 2),
                ("group_intersection", 2),
                ("group_difference", 2),
            ],
        ),
        (
            "type_free",
            "datatype",
            "newtype",
            &[
                ("type_dup", 1),
                ("type_contiguous", 2),
                ("type_vector", 4),
                ("type_create_hvector", 4),
                ("type_indexed", 4),
                ("type_create_hindexed", 4),
                ("type_create_indexed_block", 4),
                ("type_create_hindexed_block", 4),
                ("type_create_struct", 4),
                ("type_create_resized", 3),
                ("type_create_subarray", 6),
            ],
        ),
        (
            "win_free",
            "win",
            "win",
            &[
                ("win_create", 5),
                ("win_allocate", 5),
                ("win_allocate_shared", 5),
                ("win_create_dynamic", 2),
            ],
        ),
        ("op_free", "op", "op", &[("op_create", 2)]),
        (
            "info_free",
            "info",
            "newinfo",
            &[("info_create", 0), ("info_dup", 1)],
        ),
        ("file_close", "fh", "fh", &[("file_open", 4)]),
        (
            "session_finalize",
            "session",
            "session",
            &[("session_init", 2)],
        ),
    ];
    for (release, old, new, acquisitions) in mpi_created {
        for (acquire, index) in acquisitions {
            let keyword = if matches!(*acquire, "comm_group" | "comm_remote_group" | "info_create")
            {
                old
            } else {
                new
            };
            pairs.insert(
                format!("mpi_{acquire}"),
                ResourcePair {
                    enabled: true,
                    family: Family::Mpi,
                    release: format!("mpi_{release}"),
                    borrowed: false,
                    acquire_handle: ArgumentSelector::new(*index, &[keyword]),
                    release_handle: ArgumentSelector::new(0, &[old]),
                    acquire_status: Some(ArgumentSelector::new(index + 1, &["ierror", "ierr"])),
                    release_status: Some(ArgumentSelector::new(1, &["ierror", "ierr"])),
                    acquire_context: Vec::new(),
                    release_context: Vec::new(),
                },
            );
        }
    }
    pairs
}

fn mpi_signatures() -> BTreeMap<String, MpiSignature> {
    let statuses = [
        ("abort", 2),
        ("init", 0),
        ("init_thread", 2),
        ("initialized", 1),
        ("finalize", 0),
        ("finalized", 1),
        ("comm_rank", 2),
        ("comm_size", 2),
        ("comm_compare", 3),
        ("comm_dup", 2),
        ("comm_split", 4),
        ("comm_free", 1),
        ("comm_group", 2),
        ("group_free", 1),
        ("type_free", 1),
        ("type_commit", 1),
        ("win_free", 1),
        ("send", 6),
        ("bsend", 6),
        ("ssend", 6),
        ("rsend", 6),
        ("recv", 7),
        ("sendrecv", 11),
        ("start", 1),
        ("startall", 2),
        ("request_free", 1),
        ("cancel", 1),
        ("wait", 2),
        ("test", 3),
        ("waitall", 3),
        ("testall", 4),
        ("waitany", 4),
        ("testany", 5),
        ("waitsome", 5),
        ("testsome", 5),
        ("error_class", 2),
        ("error_string", 3),
        ("query_thread", 1),
        ("comm_set_errhandler", 2),
    ];
    let status = |index| Some(ArgumentSelector::new(index, &["ierror", "ierr"]));
    let mut table: BTreeMap<_, _> = statuses
        .into_iter()
        .map(|(name, index)| {
            (
                format!("mpi_{name}"),
                MpiSignature {
                    status: status(index),
                    ..MpiSignature::default()
                },
            )
        })
        .collect();
    table.get_mut("mpi_comm_rank").unwrap().rank = Some(RankSignature {
        enabled: true,
        rank: ArgumentSelector::new(1, &["rank"]),
        communicator: ArgumentSelector::new(0, &["comm"]),
    });
    let operations = [
        ("start", RequestAction::Start, 0, false, None),
        ("startall", RequestAction::Start, 1, true, None),
        ("wait", RequestAction::Wait, 0, false, None),
        ("waitall", RequestAction::Wait, 1, true, None),
        ("test", RequestAction::Test, 0, false, Some(1)),
        ("testall", RequestAction::Test, 1, true, Some(2)),
        ("request_free", RequestAction::Free, 0, false, None),
        ("cancel", RequestAction::Cancel, 0, false, None),
    ];
    for (name, action, handle, array, flag) in operations {
        table.get_mut(&format!("mpi_{name}")).unwrap().request = Some(RequestSignature {
            enabled: true,
            action,
            handle: ArgumentSelector::new(
                handle,
                &[if array {
                    "array_of_requests"
                } else {
                    "request"
                }],
            ),
            count: array.then(|| ArgumentSelector::new(0, &["count"])),
            completion_flag: flag.map(|index| ArgumentSelector::new(index, &["flag"])),
            collective: false,
        });
    }
    for routine in ["send", "bsend", "ssend", "rsend", "recv"] {
        for (name, action) in [
            (format!("mpi_i{routine}"), RequestAction::Post),
            (format!("mpi_{routine}_init"), RequestAction::Initialize),
        ] {
            table.insert(
                name,
                MpiSignature {
                    status: status(7),
                    request: Some(RequestSignature {
                        enabled: true,
                        action,
                        handle: ArgumentSelector::new(6, &["request"]),
                        count: None,
                        completion_flag: None,
                        collective: false,
                    }),
                    ..MpiSignature::default()
                },
            );
        }
    }
    let collectives = [
        ("barrier", 0),
        ("bcast", 4),
        ("reduce", 6),
        ("allreduce", 5),
        ("gather", 7),
        ("gatherv", 8),
        ("scatter", 7),
        ("scatterv", 8),
        ("allgather", 6),
        ("allgatherv", 7),
        ("alltoall", 6),
        ("alltoallv", 8),
        ("alltoallw", 8),
        ("scan", 5),
        ("exscan", 5),
        ("reduce_scatter", 5),
        ("reduce_scatter_block", 5),
        ("neighbor_allgather", 6),
        ("neighbor_allgatherv", 7),
        ("neighbor_alltoall", 6),
        ("neighbor_alltoallv", 8),
        ("neighbor_alltoallw", 8),
    ];
    for (routine, index) in collectives {
        let mut compare = Vec::new();
        if matches!(
            routine,
            "bcast" | "reduce" | "gather" | "gatherv" | "scatter" | "scatterv"
        ) {
            compare.push(ArgumentSelector::new(index - 1, &["root"]));
        }
        if routine == "bcast" {
            compare.extend([
                ArgumentSelector::new(1, &["count"]),
                ArgumentSelector::new(2, &["datatype"]),
            ]);
        }
        if matches!(routine, "allreduce" | "reduce" | "scan" | "exscan") {
            compare.extend([
                ArgumentSelector::new(2, &["count"]),
                ArgumentSelector::new(3, &["datatype"]),
                ArgumentSelector::new(4, &["op"]),
            ]);
        }
        let collective = CollectiveSignature {
            enabled: true,
            communicator: ArgumentSelector::new(index, &["comm"]),
            compare_arguments: compare,
        };
        table.insert(
            format!("mpi_{routine}"),
            MpiSignature {
                status: status(index + 1),
                collective: Some(collective.clone()),
                ..MpiSignature::default()
            },
        );
        for (name, action, handle) in [
            (format!("mpi_i{routine}"), RequestAction::Post, index + 1),
            (
                format!("mpi_{routine}_init"),
                RequestAction::Initialize,
                index + 2,
            ),
        ] {
            table.insert(
                name,
                MpiSignature {
                    status: status(handle + 1),
                    collective: Some(collective.clone()),
                    request: Some(RequestSignature {
                        enabled: true,
                        action,
                        handle: ArgumentSelector::new(handle, &["request"]),
                        count: None,
                        completion_flag: None,
                        collective: true,
                    }),
                    ..MpiSignature::default()
                },
            );
        }
    }
    let outputs: [(&str, &[(i32, &str)]); 11] = [
        ("comm_rank", &[(1, "rank")]),
        ("comm_size", &[(1, "size")]),
        ("initialized", &[(0, "flag")]),
        ("finalized", &[(0, "flag")]),
        ("test", &[(1, "flag"), (2, "status")]),
        ("testall", &[(2, "flag"), (3, "array_of_statuses")]),
        ("testany", &[(3, "flag"), (2, "index"), (4, "status")]),
        ("wait", &[(1, "status")]),
        ("waitall", &[(2, "array_of_statuses")]),
        ("bcast", &[(0, "buffer")]),
        ("allreduce", &[(1, "recvbuf")]),
    ];
    for (routine, slots) in outputs {
        table.get_mut(&format!("mpi_{routine}")).unwrap().outputs = slots
            .iter()
            .map(|(index, keyword)| ArgumentSelector::new(*index, &[*keyword]))
            .collect();
    }
    table
}

#[derive(Debug, Clone)]
pub struct Argument<'a> {
    pub node: Node<'a>,
    pub keyword: Option<String>,
}

pub fn arguments<'a>(unit: &Unit, node: Node<'a>) -> Vec<Argument<'a>> {
    child(node, "argument_list").map_or_else(Vec::new, |items| {
        children(items)
            .into_iter()
            .filter_map(|item| {
                if item.kind() == "keyword_argument" {
                    Some(Argument {
                        node: item.child_by_field_name("value")?,
                        keyword: Some(unit.name(item.child_by_field_name("name")?)),
                    })
                } else {
                    Some(Argument {
                        node: item,
                        keyword: None,
                    })
                }
            })
            .collect()
    })
}

#[derive(Debug, Clone)]
pub struct Call<'a> {
    pub family: Family,
    pub name: String,
    pub arguments: Vec<Argument<'a>>,
}

impl<'a> Call<'a> {
    pub fn argument(&self, selector: &ArgumentSelector) -> Option<Node<'a>> {
        if !selector.enabled {
            return None;
        }
        self.arguments
            .iter()
            .find(|argument| {
                argument
                    .keyword
                    .as_ref()
                    .is_some_and(|name| selector.keywords.contains(name))
            })
            .map(|argument| argument.node)
            .or_else(|| {
                let count = self
                    .arguments
                    .iter()
                    .take_while(|argument| argument.keyword.is_none())
                    .count();
                let index = if selector.index < 0 {
                    count.checked_add_signed(selector.index as isize)?
                } else {
                    selector.index as usize
                };
                (index < count).then(|| self.arguments[index].node)
            })
    }
}

fn family(policy: &SafetyPolicy, module: &str, name: &str) -> Option<Family> {
    if let Some(family) = policy.call_families.get(&format!("{module}::{name}")) {
        return Some(*family);
    }
    let petsc = policy
        .petsc_modules
        .iter()
        .any(|pattern| matches(pattern, module));
    let mpi = petsc
        || policy
            .mpi_modules
            .iter()
            .any(|pattern| matches(pattern, module));
    if mpi
        && policy
            .mpi_prefixes
            .iter()
            .any(|prefix| name.starts_with(prefix))
    {
        Some(Family::Mpi)
    } else if petsc
        && policy
            .petsc_prefixes
            .iter()
            .any(|prefix| name.starts_with(prefix))
        && !name.contains('_')
    {
        Some(Family::Petsc)
    } else {
        None
    }
}

pub fn call_name(project: &Project<'_>, scope: usize, node: Node<'_>) -> Option<String> {
    let info = &project.scopes[scope];
    let unit = &project.units[info.unit];
    if let Some(symbol) = project.expression(info.unit, node) {
        return Some(project.symbols[symbol].name.clone());
    }
    let target = node.child_by_field_name("subroutine")?;
    let local = unit.name(target);
    Some(
        project
            .lookup(scope, &local)
            .map_or(local, |symbol| project.symbols[symbol].name.clone()),
    )
}

pub fn read_call<'a>(
    project: &Project<'a>,
    config: &Config,
    scope: usize,
    node: Node<'a>,
) -> Option<Call<'a>> {
    if node.kind() != "subroutine_call" {
        return None;
    }
    let target = node.child_by_field_name("subroutine")?;
    if target.kind() != "identifier" {
        return None;
    }
    let info = &project.scopes[scope];
    let unit = &project.units[info.unit];
    let local = unit.name(target);
    let symbol = project
        .expression(info.unit, node)
        .or_else(|| project.lookup(scope, &local));
    let (name, family) = if let Some(id) = symbol {
        let symbol = &project.symbols[id];
        let recognized = if let Some(module) = symbol.external.as_deref() {
            family(&config.safety, module, &symbol.name)
        } else {
            config.safety.call_families.get(&symbol.qualified).copied()
        };
        (symbol.name.clone(), recognized?)
    } else {
        let mut candidates = BTreeSet::new();
        let mut parent = Some(scope);
        let mut seen = BTreeSet::new();
        while let Some(current) = parent {
            if !seen.insert(current) {
                break;
            }
            for import in &project.scopes[current].uses {
                if project.modules.contains_key(&import.module) {
                    continue;
                }
                let remote = import.imports.get(&local).map(String::as_str);
                if import.only && remote.is_none() {
                    continue;
                }
                if remote.is_none() && import.imports.values().any(|name| name == &local) {
                    continue;
                }
                let remote = remote.unwrap_or(&local);
                if let Some(family) = family(&config.safety, &import.module, remote) {
                    candidates.insert((remote.to_owned(), family));
                }
            }
            parent = project.host(current);
        }
        if candidates.len() != 1 {
            return None;
        }
        candidates.into_iter().next()?
    };
    Some(Call {
        family,
        name,
        arguments: arguments(unit, node),
    })
}

#[derive(Debug, Clone)]
pub struct ResourceAccess<'a> {
    pub handle: Node<'a>,
    pub acquire: bool,
    pub borrowed: bool,
    pub release: String,
    pub context: Vec<Option<Node<'a>>>,
    pub status: Option<Node<'a>>,
}

pub fn resource_access<'a>(policy: &SafetyPolicy, call: &Call<'a>) -> Option<ResourceAccess<'a>> {
    let acquire = policy
        .resources
        .get(&call.name)
        .filter(|pair| pair.enabled && pair.family == call.family);
    let release = || {
        policy
            .resources
            .values()
            .find(|pair| pair.enabled && pair.family == call.family && pair.release == call.name)
    };
    let (pair, is_acquire) = if let Some(pair) = acquire {
        (pair, true)
    } else {
        (release()?, false)
    };
    let (handle, context, status) = if is_acquire {
        (
            &pair.acquire_handle,
            &pair.acquire_context,
            &pair.acquire_status,
        )
    } else {
        (
            &pair.release_handle,
            &pair.release_context,
            &pair.release_status,
        )
    };
    let explicit_status = status
        .as_ref()
        .filter(|selector| selector.enabled)
        .is_some_and(|selector| {
            call.arguments.iter().any(|arg| {
                arg.keyword
                    .as_ref()
                    .is_some_and(|keyword| selector.keywords.contains(keyword))
            })
        });
    let mut handle = handle.clone();
    if pair.family == Family::Petsc && is_acquire && handle.index == -2 && explicit_status {
        handle.index = -1;
    }
    Some(ResourceAccess {
        handle: call.argument(&handle)?,
        acquire: is_acquire,
        borrowed: pair.borrowed,
        release: pair.release.clone(),
        context: context
            .iter()
            .map(|selector| call.argument(selector))
            .collect(),
        status: status.as_ref().and_then(|selector| call.argument(selector)),
    })
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Part {
    Component(String),
    Indices(Vec<Option<i64>>),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Storage {
    pub root: usize,
    pub parts: Vec<Part>,
}

pub fn int(text: &str) -> Option<i64> {
    let compact: String = text.chars().filter(|ch| !ch.is_whitespace()).collect();
    compact.split('_').next()?.parse().ok()
}

pub fn storage(project: &Project<'_>, scope: usize, node: Node<'_>) -> Option<Storage> {
    if node.parent().is_some_and(|parent| {
        parent.kind() == "keyword_argument" && parent.child_by_field_name("name") == Some(node)
    }) {
        return None;
    }
    let info = &project.scopes[scope];
    let unit = &project.units[info.unit];
    fn designator(unit: &Unit, node: Node<'_>, parts: &mut Vec<Part>) -> Option<String> {
        let node = unwrap(node);
        match node.kind() {
            "identifier" | "name" => Some(unit.name(node)),
            "derived_type_member_expression" => {
                let members = children(node);
                let (root, component) = (members.first()?, members.get(1)?);
                let name = designator(unit, *root, parts)?;
                parts.push(Part::Component(unit.name(*component)));
                Some(name)
            }
            "call_expression" | "sized_allocation" => {
                let name = designator(unit, node.named_child(0)?, parts)?;
                let nodes = child(node, "argument_list").or_else(|| child(node, "size"));
                let indices_nodes = nodes.map(children).unwrap_or_default();
                if indices_nodes.len() == 1 && unit.name(indices_nodes[0]) == ":" {
                    return Some(name);
                }
                let indices: Vec<_> = indices_nodes
                    .into_iter()
                    .map(|arg| int(unit.text(arg)))
                    .collect();
                if !indices.is_empty() {
                    parts.push(Part::Indices(indices));
                }
                Some(name)
            }
            _ => None,
        }
    }

    let mut parts = Vec::new();
    let name = designator(unit, node, &mut parts)?;
    let root_node = base(node);
    let root = project
        .expression(info.unit, root_node)
        .or_else(|| project.lookup(scope, &name))?;
    let mut parent = root_node.parent();
    while let Some(construct) = parent {
        if construct == info.node {
            break;
        }
        if construct.kind() == "associate_statement"
            && let Some(list) = child(construct, "association_list")
                .filter(|list| list.end_byte() <= root_node.start_byte())
            && let Some(selector) = children(list).into_iter().find_map(|item| {
                let alias = item.child_by_field_name("name")?;
                let selector = item.child_by_field_name("selector")?;
                (unit.name(alias) == name && project.expression(info.unit, selector) == Some(root))
                    .then_some(selector)
            })
            && let Some(mut selected) = storage(project, scope, selector)
        {
            selected.parts.extend(parts);
            return Some(selected);
        }
        parent = construct.parent();
    }
    // A parameter may serve as a comparison constant but cannot carry mutable storage.
    if !matches!(
        project.symbols[root].kind,
        "variable" | "association" | "associate"
    ) || project.symbols[root].attrs.contains("parameter")
    {
        return None;
    }
    Some(Storage { root, parts })
}

pub fn base(mut node: Node<'_>) -> Node<'_> {
    node = unwrap(node);
    while matches!(
        node.kind(),
        "derived_type_member_expression" | "call_expression" | "sized_allocation"
    ) {
        let Some(first) = node.named_child(0) else {
            break;
        };
        node = first;
    }
    node
}

pub fn assigned(node: Node<'_>) -> Option<(Node<'_>, Node<'_>)> {
    if !matches!(
        node.kind(),
        "assignment_statement" | "pointer_association_statement"
    ) {
        return None;
    }
    Some((
        node.child_by_field_name("left")?,
        node.child_by_field_name("right")?,
    ))
}

pub fn executable(scope: &Scope<'_>) -> bool {
    !scope.interface_body
        && matches!(
            scope.kind,
            "program" | "subroutine" | "function" | "module_procedure"
        )
}

pub fn body_nodes<'a>(scope: &Scope<'a>) -> Vec<Node<'a>> {
    fn body_children(node: Node<'_>) -> Vec<Node<'_>> {
        let mut cursor = node.walk();
        node.named_children(&mut cursor).collect()
    }
    let mut pending: Vec<_> = body_children(scope.node).into_iter().rev().collect();
    let mut nodes = Vec::new();
    while let Some(node) = pending.pop() {
        if node == scope.header
            || matches!(
                node.kind(),
                "internal_procedures" | "interface" | "derived_type_definition"
            )
        {
            continue;
        }
        nodes.push(node);
        pending.extend(body_children(node).into_iter().rev());
    }
    nodes
}

pub fn declarator<'a>(symbol: &Symbol<'a>) -> Node<'a> {
    let mut node = symbol.node;
    while let Some(parent) = node.parent() {
        if !matches!(
            parent.kind(),
            "init_declarator"
                | "pointer_init_declarator"
                | "sized_declarator"
                | "coarray_declarator"
        ) {
            break;
        }
        node = parent;
    }
    node
}

pub fn escapes(project: &Project<'_>, scope: usize, place: &Storage) -> bool {
    let symbol = &project.symbols[place.root];
    let owner = &project.scopes[scope];
    symbol.scope != scope
        || owner.parameters.contains(&symbol.name)
        || owner.result.as_ref() == Some(&symbol.name)
        || (owner.kind != "program"
            && (symbol.attrs.contains("save")
                || matches!(
                    declarator(symbol).kind(),
                    "init_declarator" | "pointer_init_declarator"
                )))
}

fn constant_name(project: &Project<'_>, scope: usize, node: Node<'_>) -> String {
    let info = &project.scopes[scope];
    project
        .expression(info.unit, node)
        .filter(|id| project.symbols[*id].external.is_some())
        .map_or_else(
            || project.units[info.unit].name(node),
            |id| project.symbols[id].name.clone(),
        )
}

pub fn identity(project: &Project<'_>, scope: usize, node: Node<'_>) -> String {
    let info = &project.scopes[scope];
    let unit = &project.units[info.unit];
    let node = unwrap(node);
    if let Some(place) = storage(project, scope, node)
        && place
            .parts
            .iter()
            .all(|part| !matches!(part, Part::Indices(items) if items.contains(&None)))
    {
        return format!("#{}{:?}", place.root, place.parts);
    }
    if matches!(node.kind(), "identifier" | "name") {
        return constant_name(project, scope, node);
    }
    unit.name(node)
}

pub fn operator(unit: &Unit, left: Node<'_>, right: Node<'_>) -> String {
    unit.source
        .text
        .get(left.end_byte()..right.start_byte())
        .unwrap_or_default()
        .chars()
        .filter(|item| !item.is_whitespace() && *item != '&')
        .collect::<String>()
        .to_ascii_lowercase()
}

pub fn status_tests(
    project: &Project<'_>,
    config: &Config,
    scope: usize,
    edge: &Edge<'_>,
) -> Vec<(Storage, bool)> {
    fn test(
        project: &Project<'_>,
        policy: &SafetyPolicy,
        scope: usize,
        node: Node<'_>,
        truth: bool,
        out: &mut Vec<(Storage, bool)>,
    ) {
        let node = unwrap(node);
        let unit = &project.units[project.scopes[scope].unit];
        if unit.name(node).starts_with(".not.") {
            if let Some(part) = children(node).last() {
                test(project, policy, scope, *part, !truth, out);
            }
            return;
        }
        let (Some(left), Some(right)) = (
            node.child_by_field_name("left"),
            node.child_by_field_name("right"),
        ) else {
            return;
        };
        let op = operator(unit, left, right);
        if (op == ".and." && truth) || (op == ".or." && !truth) {
            test(project, policy, scope, left, truth, out);
            test(project, policy, scope, right, truth, out);
            return;
        }
        if !matches!(op.as_str(), "==" | ".eq." | "/=" | ".ne.") {
            return;
        }
        let success = |item| {
            let name = constant_name(project, scope, unwrap(item));
            int(&name) == Some(0) || policy.success_constants.contains(&name)
        };
        let place = if success(left) {
            storage(project, scope, right)
        } else if success(right) {
            storage(project, scope, left)
        } else {
            None
        };
        if let Some(place) = place {
            out.push((place, truth == matches!(op.as_str(), "==" | ".eq.")));
        }
    }
    let mut tests = Vec::new();
    if let Some(condition) = edge.condition {
        test(
            project,
            &config.safety,
            scope,
            condition,
            edge.truth,
            &mut tests,
        );
    }
    tests
}
