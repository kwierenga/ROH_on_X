"""
md_to_pdf.py — render a Markdown file to a clean PDF via headless Edge.
======================================================================
No LaTeX / pandoc / reportlab. Markdown -> styled HTML -> PDF (Chromium print).

Usage:  python md_to_pdf.py "Union typing paper (autosomes + X).md" [out.pdf]
        python md_to_pdf.py file1.md file2.md     # batch (each -> same basename .pdf)
"""
import os, sys, subprocess, tempfile
import markdown

EDGE = next((p for p in [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
] if os.path.exists(p)), None)

CSS = """
@page { size: A4; margin: 18mm 17mm; }
* { -webkit-print-color-adjust: exact; }
body { font-family: Georgia, 'Times New Roman', serif; font-size: 10.7pt;
       line-height: 1.5; color: #1a1a1a; max-width: 100%; }
h1 { font-family: 'Segoe UI', Arial, sans-serif; font-size: 19pt; line-height:1.25;
     margin: 0 0 .3em; color:#111; }
h2 { font-family: 'Segoe UI', Arial, sans-serif; font-size: 14pt; margin: 1.4em 0 .4em;
     border-bottom: 1.5px solid #ccc; padding-bottom: 3px; color:#111; }
h3 { font-family: 'Segoe UI', Arial, sans-serif; font-size: 11.5pt; margin: 1.1em 0 .3em; color:#222; }
p, li { orphans: 2; widows: 2; }
em { color:#333; }
strong { color:#000; }
code { font-family: 'Cascadia Mono', Consolas, monospace; font-size: 9.3pt;
       background:#f4f4f4; padding:1px 4px; border-radius:3px; }
pre { background:#f6f8fa; padding:10px 12px; border-radius:6px; overflow-x:auto;
      font-size:9pt; border:1px solid #e5e5e5; }
pre code { background:none; padding:0; }
table { border-collapse: collapse; margin: .8em 0; font-size: 9.8pt; width:auto; }
th, td { border: 1px solid #c9c9c9; padding: 4px 9px; text-align: left; }
th { background:#f0f0f0; font-family:'Segoe UI',Arial,sans-serif; }
tr:nth-child(even) td { background:#fafafa; }
hr { border:none; border-top:1px solid #ddd; margin:1.5em 0; }
blockquote { border-left:3px solid #ccc; margin:.8em 0; padding:.2em 1em; color:#444; }
h2, h3 { page-break-after: avoid; }
table, pre { page-break-inside: avoid; }
"""

def convert(md_path, pdf_path=None):
    if pdf_path is None:
        pdf_path = os.path.splitext(md_path)[0] + ".pdf"
    with open(md_path, encoding="utf-8") as f:
        body = markdown.markdown(f.read(),
                                 extensions=["tables", "fenced_code", "sane_lists", "smarty"])
    html = (f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
            f"<style>{CSS}</style></head><body>{body}</body></html>")
    tmp = tempfile.NamedTemporaryFile("w", suffix=".html", delete=False, encoding="utf-8")
    tmp.write(html); tmp.close()
    url = "file:///" + tmp.name.replace("\\", "/")
    out_abs = os.path.abspath(pdf_path)
    subprocess.run([EDGE, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={out_abs}", url],
                   check=True, capture_output=True, timeout=120)
    os.unlink(tmp.name)
    print(f"wrote {pdf_path}  ({os.path.getsize(out_abs)//1024} KB)")

if __name__ == "__main__":
    if EDGE is None:
        sys.exit("Microsoft Edge not found.")
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    for a in args:
        convert(a)
