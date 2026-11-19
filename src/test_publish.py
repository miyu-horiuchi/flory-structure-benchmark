import unittest

from publish import build_html


class PublishingTests(unittest.TestCase):
    def test_replaces_static_title_without_changing_report(self) -> None:
        head = '<html><head><meta charset="utf-8"><title>Artifact</title>'
        body = '<body><script>const report = "unchanged";</script></body></html>'
        metadata = (
            '<title>Benchmark</title>\n<meta property="og:title" content="Benchmark">'
        )
        result = build_html(head + "</head>" + body, metadata)
        self.assertNotIn("<title>Artifact</title>", result)
        self.assertIn(metadata, result.split("</head>")[0])
        self.assertIn('<meta charset="utf-8">', result)
        self.assertEqual(result.split("</head>", 1)[1], body)

    def test_rejects_missing_head(self) -> None:
        with self.assertRaisesRegex(ValueError, "head"):
            build_html("<title>Artifact</title>", "<title>Benchmark</title>")

    def test_rejects_missing_title(self) -> None:
        with self.assertRaisesRegex(ValueError, "title"):
            build_html("<head></head><body></body>", "<title>Benchmark</title>")


if __name__ == "__main__":
    unittest.main()
