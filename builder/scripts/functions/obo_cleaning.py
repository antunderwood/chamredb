"""Lossless clean-up of the CARD ontology before parsing (ticket #12). Standard library only."""
import re

# `property_value: <relation>` with no value is invalid OBO; CARD 4.0.2 writes 1,402 of them and pronto refuses the file.
VALUELESS = re.compile(r"^property_value: \S+[ \t]*$")


def strip_valueless_property_values(text, max_removed=5000):
    """Return (text without valueless property_value lines, number removed). Other lines are untouched; more than
    max_removed removals means the file is not the problem we know and the build stops."""
    kept, removed = [], 0
    for line in text.splitlines(keepends=True):
        if VALUELESS.match(line.rstrip("\r\n")):
            removed += 1
        else:
            kept.append(line)
    if removed > max_removed:
        raise ValueError(f"removed more than {max_removed} valueless property_value lines ({removed}); not the known CARD defect")
    return "".join(kept), removed
