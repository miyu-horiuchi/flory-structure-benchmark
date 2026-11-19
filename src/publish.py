"""Build the hosted page from an untouched compiled report and static metadata."""

import argparse
import re
from pathlib import Path


def build_html(compiled: str, metadata: str) -> str:
    head, separator, body = compiled.partition("</head>")
    if not separator:
        raise ValueError("Compiled report has no closing head tag")
    head, count = re.subn(
        r"<title>[^<]*</title>", lambda _: metadata.strip(), head, count=1
    )
    if not count:
        raise ValueError("Compiled report has no static title")
    return head + separator + body


def main() -> None:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("compiled", type=Path)
    parser.add_argument("--output", type=Path, default=root.parent / "index.html")
    args = parser.parse_args()
    if args.compiled.resolve() == args.output.resolve():
        parser.error("Keep the compiled artifact separate from the hosted output")
    metadata = (root / "preview.head.html").read_text(encoding="utf-8")
    output = build_html(args.compiled.read_text(encoding="utf-8"), metadata)
    args.output.write_text(output, encoding="utf-8")


if __name__ == "__main__":
    main()
