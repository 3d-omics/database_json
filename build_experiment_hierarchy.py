"""Utility to build a simplified experiment hierarchy JSON file."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List


def generate_experiment_hierarchy(
    input_path: str | Path, output_path: str | Path
) -> Dict[str, Any]:
    """
    Build a hierarchical JSON structure for the 3D'omics project experiments.

    The input is expected to be the Airtable export stored in
    ``data/animaltrialexperiment.json`` plus specimens in
    ``data/animalspecimen.json``, macrosamples in
    ``data/intestinalsectionsample.json``, cryosections in
    ``data/cryosection.json``, and microsamples in ``data/microsample.json``.
    The resulting structure keeps only the experiment ``ID``, ``Name``,
    ``Bioproject accession``, scheduling metadata, linked specimens (ID,
    Biosample accession, treatment code, weight), linked macrosamples (sample
    type, ENA accession, Metabolights accession), linked cryosections (date),
    linked microsamples (date plus linked cryosections), and adds static project
    metadata.

    Parameters
    ----------
    input_path:
        Path to the raw ``animaltrialexperiment.json`` file.
    output_path:
        Destination path for the generated JSON file.

    Returns
    -------
    dict
        The hierarchy that was written to ``output_path``.
    """

    input_file = Path(input_path)
    specimens_file = input_file.parent / "animalspecimen.json"
    macrosamples_file = input_file.parent / "intestinalsectionsample.json"
    cryosections_file = input_file.parent / "cryosection.json"
    microsamples_file = input_file.parent / "microsample.json"
    output_file = Path(output_path)

    with input_file.open("r", encoding="utf-8") as handle:
        records: List[Dict[str, Any]] = json.load(handle)

    experiments: Dict[str, Dict[str, Any]] = {}
    experiment_ids: List[str] = []
    experiment_record_to_id: Dict[str, str] = {}
    for record in records:
        fields = record.get("fields", {})
        experiment_id = fields.get("ID")
        name = fields.get("Name")
        bioproject = fields.get("Bioproject accession")
        start_date = fields.get("StartDate")
        end_date = fields.get("EndDate")
        experiment_record_id = record.get("id")

        if not (experiment_id and name and bioproject):
            # Skip entries without the required metadata.
            continue

        experiment_data: Dict[str, str] = {
            "ID": experiment_id,
            "Name": name,
            "Bioproject accession": bioproject,
        }

        # Add optional scheduling metadata when present.
        if start_date:
            experiment_data["Start date"] = start_date
        if end_date:
            experiment_data["End date"] = end_date

        if experiment_record_id:
            experiment_record_to_id[experiment_record_id] = experiment_id
        # Preserve insertion order of IDs.
        experiment_ids.append(experiment_id)
        # Store experiment keyed by its ID, omitting the redundant top-level ID inside.
        experiments[experiment_id] = {
            key: value for key, value in experiment_data.items() if key != "ID"
        }

    # Attach individuals (specimens) to their experiments.
    specimens_by_experiment: Dict[str, List[str]] = defaultdict(list)
    individuals: Dict[str, Dict[str, Any]] = {}
    with specimens_file.open("r", encoding="utf-8") as handle:
        specimen_records: List[Dict[str, Any]] = json.load(handle)

    for specimen in specimen_records:
        fields = specimen.get("fields", {})
        specimen_id = fields.get("ID")
        biosample = fields.get("Biosample accession")
        linked_experiments: List[str] = fields.get("Experiment", []) or []
        treatment_code = fields.get("Treatment_flat")
        weight = fields.get("Weight")
        sex = fields.get("Sex")

        if not (specimen_id and biosample and linked_experiments):
            continue

        individual_data: Dict[str, Any] = {"Biosample accession": biosample}

        if treatment_code:
            individual_data["Treatment"] = (
                treatment_code[0] if isinstance(treatment_code, list) else treatment_code
            )
        if weight is not None:
            individual_data["Weight"] = weight
        if sex:
            individual_data["Sex"] = sex[0] if isinstance(sex, list) else sex

        individuals[specimen_id] = individual_data

        for experiment_record_id in linked_experiments:
            experiment_id = experiment_record_to_id.get(experiment_record_id)
            if not experiment_id:
                continue

            specimens_by_experiment[experiment_id].append(specimen_id)

    for experiment_id, individual_ids in specimens_by_experiment.items():
        experiments[experiment_id]["Individual IDs"] = individual_ids

    # Attach macrosamples to individuals.
    macrosamples_by_individual: Dict[str, List[str]] = defaultdict(list)
    macrosamples: Dict[str, Dict[str, Any]] = {}
    with macrosamples_file.open("r", encoding="utf-8") as handle:
        macrosample_records: List[Dict[str, Any]] = json.load(handle)

    for sample in macrosample_records:
        fields = sample.get("fields", {})
        sample_id = fields.get("ID")
        individual_id = fields.get("Individual")
        sample_type = fields.get("Sample type")
        ena_accession = fields.get("ENA accession")
        metabolights_accession = fields.get("Metabolights accession")

        if not (sample_id and individual_id):
            continue

        sample_data: Dict[str, Any] = {}
        if sample_type:
            sample_data["Sample type"] = sample_type
        if ena_accession:
            ena_value = (
                ena_accession[0] if isinstance(ena_accession, list) else ena_accession
            )
            if ena_value:
                sample_data["ENA accession"] = ena_value
        if metabolights_accession:
            sample_data["Metabolights accession"] = metabolights_accession

        macrosamples[sample_id] = sample_data
        macrosamples_by_individual[individual_id].append(sample_id)

    for individual_id, sample_ids in macrosamples_by_individual.items():
        if individual_id not in individuals:
            continue
        individuals[individual_id]["Macrosample IDs"] = sample_ids

    # Cryosections keyed by ID with prefix map for microsample linkage.
    cryosections_by_prefix: Dict[str, List[str]] = defaultdict(list)
    cryosections: Dict[str, Dict[str, Any]] = {}
    cryosection_microsamples: Dict[str, List[str]] = defaultdict(list)
    with cryosections_file.open("r", encoding="utf-8") as handle:
        cryosection_records: List[Dict[str, Any]] = json.load(handle)

    for record in cryosection_records:
        fields = record.get("fields", {})
        cryo_id = fields.get("ID")
        date_value = fields.get("SlideDate")

        if not cryo_id:
            continue

        cryo_data: Dict[str, Any] = {}
        if date_value:
            cryo_data["Date"] = (
                date_value[0] if isinstance(date_value, list) else date_value
            )

        cryosections[cryo_id] = cryo_data
        prefix = cryo_id[:6] if len(cryo_id) >= 6 else cryo_id
        cryosections_by_prefix[prefix].append(cryo_id)

    # Microsamples keyed by ID with links to cryosections and macrosamples.
    microsamples_by_macrosample: Dict[str, List[str]] = defaultdict(list)
    microsamples: Dict[str, Dict[str, Any]] = {}
    with microsamples_file.open("r", encoding="utf-8") as handle:
        microsample_records: List[Dict[str, Any]] = json.load(handle)

    for record in microsample_records:
        fields = record.get("fields", {})
        code = fields.get("Code")
        date_value = fields.get("Date")
        x_coord = fields.get("Xcoord")
        y_coord = fields.get("Ycoord")
        size = fields.get("Size")
        batch = fields.get("LMBatch_flat")
        ena_accession = fields.get("ENA accession")
        sample_type = fields.get("Sample_type")

        if not code:
            continue

        microsample_data: Dict[str, Any] = {}
        if date_value:
            microsample_data["Date"] = date_value
        if x_coord is not None:
            microsample_data["Xcoord"] = x_coord
        if y_coord is not None:
            microsample_data["Ycoord"] = y_coord
        if size is not None:
            microsample_data["Size"] = size
        if batch:
            microsample_data["Batch"] = (
                batch[0] if isinstance(batch, list) else batch
            )
        if ena_accession:
            microsample_data["ENA accession"] = (
                ena_accession[0]
                if isinstance(ena_accession, list)
                else ena_accession
            )
        if sample_type:
            microsample_data["Sample type"] = (
                sample_type[0] if isinstance(sample_type, list) else sample_type
            )

        macrosample_id = code[:6] if len(code) >= 6 else code
        if macrosample_id in macrosamples:
            microsample_data["Macrosample ID"] = macrosample_id
            microsamples_by_macrosample[macrosample_id].append(code)

        cryo_ids = cryosections_by_prefix.get(macrosample_id)
        if cryo_ids:
            for cryo_id in cryo_ids:
                cryosection_microsamples[cryo_id].append(code)

        microsamples[code] = microsample_data

    for cryo_id, micro_ids in cryosection_microsamples.items():
        if cryo_id in cryosections:
            cryosections[cryo_id]["Microsample IDs"] = micro_ids

    for macrosample_id, micro_ids in microsamples_by_macrosample.items():
        if macrosample_id in macrosamples:
            macrosamples[macrosample_id]["Microsample IDs"] = micro_ids

    project_metadata = {
        "Name": "3D'omics",
        "Bioproject accession": "PRJEB86267",
        "Start date": "2021-09-01",
        "End date": "2025-12-31",
    }

    # Store experiments under the project to make the parent/child hierarchy explicit.
    hierarchy = {
        "Projects": {
            project_metadata["Name"]: {
                **project_metadata,
                "Experiment IDs": experiment_ids,
            }
        },
        "Experiments": experiments,
        "Individuals": individuals,
        "Macrosamples": macrosamples,
        "Cryosections": cryosections,
        "Microsamples": microsamples,
    }

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8") as handle:
        json.dump(hierarchy, handle, indent=2)

    return hierarchy


if __name__ == "__main__":
    generate_experiment_hierarchy(
        input_path="data/animaltrialexperiment.json",
        output_path="output/experiment_hierarchy.json",
    )
