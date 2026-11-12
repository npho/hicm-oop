"""Rewrite the author's historical absolute paths in source text, loudly.

The launchers copy scripts into a run and point each historical root at the
run. ``re.sub`` does nothing when its target is absent, so a root that stops
matching (an edited ``root_map.json``, an edited script) would leave the old
path in place without any error. Each rewrite therefore counts its matches,
and the caller requires the roots it depends on to have matched.
"""

import re

# A user-home path that no historical root accounts for.
HOME_PATH = re.compile(r"/(?:Users|home)/[^\s\"'`)\]]+")


def rewrite_roots(text, mapping, counts):
    """Replace each historical root in ``text`` and tally the matches.

    Args:
        text: Source text to rewrite.
        mapping: ``{historical_root: new_root}``.
        counts: ``collections.Counter`` that receives one count per root.

    Returns:
        The rewritten text.
    """
    pattern = "|".join(re.escape(old) for old in sorted(mapping, key=len, reverse=True))

    def replace(match):
        counts[match.group(0)] += 1
        return mapping[match.group(0)]

    return re.sub(pattern, replace, text)


def require_roots(counts, required, where):
    """Raise if a root that ``where`` depends on matched nothing.

    Args:
        counts: Tally filled by ``rewrite_roots``.
        required: Historical roots that must each match at least once.
        where: Name of the rewritten tree, for the error message.

    Raises:
        RuntimeError: If any required root was never found.
    """
    missing = [root for root in required if counts[root] == 0]
    if missing:
        raise RuntimeError(f"No file in {where} contains the root(s): {missing}")


def unmapped_paths(text, package):
    """List user-home paths in ``text`` that lie outside ``package``.

    These are historical locations that ``root_map.json`` does not cover.

    Args:
        text: Rewritten source text.
        package: Package directory; paths inside it are expected.

    Returns:
        Sorted list of distinct paths.
    """
    return sorted(
        {
            found
            for found in HOME_PATH.findall(text)
            if not found.startswith(str(package))
        }
    )
