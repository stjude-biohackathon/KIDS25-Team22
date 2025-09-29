def extract_coordinates(file_path, limit=None):
    """Parse a VCF and return variant coordinates as dicts."""
    coordinates = []

    with open(file_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            columns = line.rstrip().split("\t")
            if len(columns) < 5:
                continue

            chrom, pos, _vid, ref, alt = columns[:5]

            if not pos.isdigit():
                continue

            for alt_allele in alt.split(","):
                entry = {
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt_allele,
                }
                coordinates.append(entry)

                if limit is not None and len(coordinates) >= limit:
                    return coordinates

    return coordinates


def format_biomcp_variant_output(raw_output):
    """Convert BioMCP stdout into table-formatted sections for readability."""
    lines = raw_output.splitlines()
    tables = []
    heading_stack = []
    rows = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        if stripped.startswith("#"):
            if rows:
                title = " > ".join(heading_stack) or "Details"
                tables.append((title, rows))
                rows = []

            level = len(stripped) - len(stripped.lstrip("#"))
            heading = stripped[level:].strip()
            heading_stack = heading_stack[: level - 1]
            heading_stack.append(heading)
        else:
            if ":" in line:
                key, value = line.split(":", 1)
                rows.append((key.strip(), value.strip()))
            else:
                rows.append((stripped, ""))

    if rows:
        title = " > ".join(heading_stack) or "Details"
        tables.append((title, rows))

    formatted_sections = []
    for title, section_rows in tables:
        header = title or "Details"
        key_width = max(
            len("Field"),
            max((len(field) for field, _ in section_rows), default=0),
        )
        value_width = max(
            len("Value"),
            max((len(value) for _, value in section_rows), default=0),
        )

        formatted_rows = [
            f"{header}",
            f"{'Field':<{key_width}} | {'Value':<{value_width}}",
            f"{'-' * key_width}-+-{'-' * value_width}",
        ]

        for field, value in section_rows:
            formatted_rows.append(f"{field:<{key_width}} | {value:<{value_width}}")

        formatted_sections.append("\n".join(formatted_rows))

    return "\n\n".join(formatted_sections)
