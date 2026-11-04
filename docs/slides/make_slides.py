"""Build the WPG 2026.11 slide deck (16:9) for the repository front page and talks.

Run:  python docs/slides/make_slides.py
Needs python-pptx; PyMuPDF (fitz) and LibreOffice for the PDF and the title PNG.
Facts and numbers come from README.md and CHANGELOG.md, figures from docs/images/.
"""
import re
import subprocess
import sys
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, MSO_AUTO_SIZE, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

HERE = Path(__file__).resolve().parent
IMAGES = HERE.parent / "images"
OUT = HERE / "WPG_2026.11.pptx"
SOFFICE = r"C:\Program Files\LibreOffice\program\soffice.exe"
KEYNOTE_FIX_DIR = r"C:\Users\Alexey\home\src\ai_stuff\stuff\tools"

ACCENT = RGBColor(0x0B, 0x5C, 0xAD)
DARK = RGBColor(0x1E, 0x23, 0x2B)
GREY = RGBColor(0x5F, 0x66, 0x70)
LIGHT = RGBColor(0xEE, 0xF3, 0xF9)
CODE_BG = RGBColor(0xF5, 0xF6, 0xF8)
COMMENT = RGBColor(0x6A, 0x8A, 0x4A)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

FONT = "Calibri"
MONO = "Consolas"

SLIDE_W, SLIDE_H = 13.333, 7.5
MARGIN = 0.6

REPO_URL = "https://github.com/samoylv/WPG"
DOCS_URL = "https://wpg.readthedocs.io"
DOI_URL = "https://doi.org/10.1107/S160057671600995X"

prs = Presentation()
prs.slide_width = Emu(12192000)
prs.slide_height = Emu(6858000)
BLANK = prs.slide_layouts[6]


# ---------------------------------------------------------------- helpers
def no_style(shp):
    """Drop the theme style reference of an autoshape: no shadow, no theme colours."""
    st = shp._element.find(qn("p:style"))
    if st is not None:
        shp._element.remove(st)
    return shp


def new_slide():
    s = prs.slides.add_slide(BLANK)
    bg = s.background.fill
    bg.solid()
    bg.fore_color.rgb = WHITE
    return s


def textbox(slide, left, top, width, height, anchor=MSO_ANCHOR.TOP):
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.auto_size = MSO_AUTO_SIZE.NONE
    tf.vertical_anchor = anchor
    for side in ("left", "right", "top", "bottom"):
        setattr(tf, f"margin_{side}", Inches(0.05))
    return tb, tf


def add_runs(par, text, size, color=DARK, bold=False, italic=False, font=FONT):
    """Text with `inline code` spans set in the monospace font."""
    for i, part in enumerate(re.split(r"`", text)):
        if not part:
            continue
        r = par.add_run()
        r.text = part
        f = r.font
        f.size = Pt(size)
        f.bold = bold
        f.italic = italic
        f.color.rgb = ACCENT if (i % 2 and color == DARK) else color
        f.name = MONO if i % 2 else font
    return par


def title(slide, text, size=36):
    _, tf = textbox(slide, MARGIN, 0.35, SLIDE_W - 2 * MARGIN, 0.8, MSO_ANCHOR.BOTTOM)
    add_runs(tf.paragraphs[0], text, size, DARK, bold=True)
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(MARGIN + 0.05), Inches(1.2),
                                 Inches(1.2), Inches(0.07))
    bar.fill.solid()
    bar.fill.fore_color.rgb = ACCENT
    bar.line.fill.background()
    no_style(bar)


def takeaway(slide, text, size=20):
    _, tf = textbox(slide, MARGIN, 1.35, SLIDE_W - 2 * MARGIN, 0.5)
    add_runs(tf.paragraphs[0], text, size, GREY, italic=True)


def set_bullet(par, indent_in=0.3):
    pPr = par._p.get_or_add_pPr()
    pPr.set("marL", str(int(Inches(indent_in))))
    pPr.set("indent", str(-int(Inches(indent_in))))
    for tag in ("a:buNone", "a:buChar", "a:buAutoNum"):
        for el in pPr.findall(qn(tag)):
            pPr.remove(el)
    buClr = pPr.makeelement(qn("a:buClr"), {})
    clr = buClr.makeelement(qn("a:srgbClr"), {"val": str(ACCENT)})
    buClr.append(clr)
    pPr.append(buClr)
    buFont = pPr.makeelement(qn("a:buFont"), {"typeface": "Arial"})
    pPr.append(buFont)
    buChar = pPr.makeelement(qn("a:buChar"), {"char": "\u2022"})
    pPr.append(buChar)


def bullets(slide, items, left, top, width, height, size=20, space=10):
    _, tf = textbox(slide, left, top, width, height)
    for i, item in enumerate(items):
        par = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        set_bullet(par)
        par.space_after = Pt(space)
        add_runs(par, item, size)
    return tf


def picture(slide, name, left, top, width=None, height=None):
    """Place a figure keeping its pixel aspect ratio. Give width or height, not both."""
    path = IMAGES / name
    w_px, h_px = Image.open(path).size
    ratio = w_px / h_px
    if width is None:
        width = height * ratio
    height = width / ratio
    pic = slide.shapes.add_picture(str(path), Inches(left), Inches(top), Inches(width), Inches(height))
    assert abs(pic.width / pic.height - ratio) / ratio < 0.01, name
    return left, top, width, height


def caption(slide, text, left, top, width, height=0.6, size=13, align=PP_ALIGN.LEFT):
    _, tf = textbox(slide, left, top, width, height)
    par = tf.paragraphs[0]
    par.alignment = align
    add_runs(par, text, size, GREY)


def slide_number(slide, n):
    _, tf = textbox(slide, SLIDE_W - MARGIN - 2.5, SLIDE_H - 0.45, 2.5, 0.3)
    par = tf.paragraphs[0]
    par.alignment = PP_ALIGN.RIGHT
    add_runs(par, f"WPG 2026.11  |  {n}", 11, GREY)


def link(run, url):
    run.hyperlink.address = url


def code_block(slide, lines, left, top, width, height, size=12, bg=CODE_BG, leading=None):
    box = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top), Inches(width), Inches(height))
    box.fill.solid()
    box.fill.fore_color.rgb = bg
    box.line.fill.background()
    no_style(box)
    tf = box.text_frame
    tf.word_wrap = False
    tf.auto_size = MSO_AUTO_SIZE.NONE
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    for side in ("left", "right"):
        setattr(tf, f"margin_{side}", Inches(0.2))
    for side in ("top", "bottom"):
        setattr(tf, f"margin_{side}", Inches(0.1))
    for i, line in enumerate(lines):
        par = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        par.alignment = PP_ALIGN.LEFT
        par.line_spacing = (Pt(leading / 2) if not line else Pt(leading)) if leading else 1.0
        code, hash_, comment = line.partition("#")
        for text, color in ((code, DARK), (hash_ + comment, COMMENT)):
            if not text:
                continue
            r = par.add_run()
            r.text = text
            r.font.name = MONO
            r.font.size = Pt(size)
            r.font.color.rgb = color
        if not line:
            r = par.add_run()
            r.text = " "
            r.font.name = MONO
            r.font.size = Pt(size / 2 if leading else size)
    return box


def readme_quick_example():
    """Code of the 'Quick example' section of README.md and the result printed in its last comment."""
    text = (HERE.parents[1] / "README.md").read_text(encoding="utf-8")
    section = text.split("## Quick example", 1)[1]
    code = section.split("```python", 1)[1].split("```", 1)[0].strip("\n").split("\n")
    last, _, printed = code[-1].partition("#")
    code[-1] = last.rstrip()
    return code, printed.strip()


def rounded_box(slide, left, top, width, height, fill, line=None):
    shp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(left), Inches(top),
                                 Inches(width), Inches(height))
    shp.adjustments[0] = 0.12
    shp.fill.solid()
    shp.fill.fore_color.rgb = fill
    if line is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line
        shp.line.width = Pt(1.5)
    no_style(shp)
    tf = shp.text_frame
    tf.word_wrap = True
    tf.auto_size = MSO_AUTO_SIZE.NONE
    for side in ("left", "right", "top", "bottom"):
        setattr(tf, f"margin_{side}", Inches(0.12))
    return shp


# ---------------------------------------------------------------- 1. title
s = new_slide()
_, tf = textbox(s, MARGIN, 0.7, SLIDE_W - 2 * MARGIN, 1.1)
add_runs(tf.paragraphs[0], "WPG \u2014 WavePropaGator 2026.11", 48, DARK, bold=True)
bar = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(MARGIN + 0.05), Inches(1.85), Inches(1.6), Inches(0.09))
bar.fill.solid()
bar.fill.fore_color.rgb = ACCENT
bar.line.fill.background()
no_style(bar)
_, tf = textbox(s, MARGIN, 2.05, SLIDE_W - 2 * MARGIN, 0.6)
add_runs(tf.paragraphs[0], "Wave optics of X-ray beamlines in Python, from the XFEL source to the sample", 26, GREY)
_, tf = textbox(s, MARGIN, 2.65, SLIDE_W - 2 * MARGIN, 0.5)
add_runs(tf.paragraphs[0], "Built on SRW by O. Chubar. Developed for the European XFEL beamlines.", 20, GREY)

fw = 9.6
l, t, w, h = picture(s, "kb_focus.png", (SLIDE_W - fw) / 2, 3.2, width=fw)
caption(s, "A 9 keV beam focused by a pair of KB mirrors (elliptical, B4C coating, 3 mrad): (a) at the mirrors, "
           "50 m from the source, (b) focus with ideal mirrors, (c) focus with 1.5 nm RMS height errors. "
           "Each case takes 0.2 s on a laptop.",
        l, t + h + 0.03, w, 0.5, size=12)
_, tf = textbox(s, MARGIN, SLIDE_H - 0.62, SLIDE_W - 2 * MARGIN, 0.45)
par = tf.paragraphs[0]
par.alignment = PP_ALIGN.CENTER
add_runs(par, REPO_URL, 22, ACCENT, bold=True)
link(par.runs[0], REPO_URL)

# ---------------------------------------------------------------- 2. what WPG is
s = new_slide()
title(s, "What WPG is")
takeaway(s, "A beamline is a Python list of optical elements, a wavefront is a Python object stored in HDF5")

boxes = [
    ("Source / pulse", "Gaussian beams, GENESIS output, XFEL Photon Pulses Database"),
    ("Beamline", "list of elements: drifts, lenses, CRL, mirrors, gratings"),
    ("Wavefront", "parameters, field and history in one HDF5 file"),
    ("Analysis", "beam sizes, plots, the next code of a start-to-end chain (SIMEX)"),
]
bw, bh, gap, top = 2.55, 2.3, 0.58, 2.35
x0 = (SLIDE_W - (4 * bw + 3 * gap)) / 2
for i, (head, body) in enumerate(boxes):
    x = x0 + i * (bw + gap)
    fill = ACCENT if i in (1, 2) else LIGHT
    shp = rounded_box(s, x, top, bw, bh, fill, None if i in (1, 2) else ACCENT)
    tf = shp.text_frame
    tf.vertical_anchor = MSO_ANCHOR.TOP
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    p.space_after = Pt(8)
    add_runs(p, head, 22, WHITE if i in (1, 2) else ACCENT, bold=True)
    p = tf.add_paragraph()
    p.alignment = PP_ALIGN.CENTER
    add_runs(p, body, 18, WHITE if i in (1, 2) else DARK)
    if i < 3:
        arr = s.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, Inches(x + bw + 0.08), Inches(top + bh / 2 - 0.22),
                                 Inches(gap - 0.16), Inches(0.44))
        arr.fill.solid()
        arr.fill.fore_color.rgb = GREY
        arr.line.fill.background()
        no_style(arr)

# SRW engine bar under the source, beamline and wavefront boxes
ex, ew = x0, 3 * bw + 2 * gap
eng = rounded_box(s, ex, top + bh + 0.3, ew, 0.75, DARK)
tf = eng.text_frame
tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p = tf.paragraphs[0]
p.alignment = PP_ALIGN.CENTER
add_runs(p, "SRW engine (O. Chubar) propagates the field: serial, OpenMP or CUDA", 20, WHITE)
for r in p.runs:
    r.font.color.rgb = WHITE

_, tf = textbox(s, MARGIN, 6.35, SLIDE_W - 2 * MARGIN, 0.5)
add_runs(tf.paragraphs[0], "Wave optics: diffraction, mirror surface errors and the finite size of optical elements "
         "are in the result",
         20, GREY)
slide_number(s, 2)

# ---------------------------------------------------------------- 3. quick example
s = new_slide()
title(s, "A beamline in a few lines")
code, printed = readme_quick_example()
assert "1.17e-06" in printed and "1.03e-06" in printed, printed  # the numbers in the line below
n_full = sum(1 for c in code if c)
code_top, code_h = 1.35, 5.15
leading = min(13.5, (code_h - 0.2) * 72 / (n_full + (len(code) - n_full) / 2))
assert leading >= 12.5, f"README example too long for one slide at 12 pt: leading {leading:.1f} pt"
code_block(s, code, MARGIN - 0.2, code_top, SLIDE_W - 2 * MARGIN + 0.4, code_h, size=12, leading=leading)
res = rounded_box(s, MARGIN - 0.2, 6.55, SLIDE_W - 2 * MARGIN + 0.4, 0.5, LIGHT)
tf = res.text_frame
tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p = tf.paragraphs[0]
add_runs(p, "Output:  ", 18, ACCENT, bold=True)
add_runs(p, "`" + printed + "`", 16)
add_runs(p, "   a 9 keV focus of 1.17 × 1.03 µm FWHM after a KB pair", 18, DARK)
slide_number(s, 3)

# ---------------------------------------------------------------- 4. real mirrors
s = new_slide()
title(s, "Real mirrors: coating reflectivity")
takeaway(s, "Reflectivity from `wpg.materials` differs from the CXRO calculator by less than 0.007 at 1\u20136 mrad")
l, t, w, h = picture(s, "mirror_reflectivity.png", MARGIN, 2.0, width=6.6)
caption(s, "Reflectivity of thick mirrors versus grazing angle: B4C and Pt at 9 keV, Si at 12 keV. "
           "Lines: `wpg.materials`, circles: the CXRO calculator.",
        l, t + h + 0.05, w, 0.6)
bullets(s, [
    "Fresnel reflectivity of a thick coating from xraylib optical constants",
    "N\u00e9vot\u2013Croce roughness, \u03c3 and \u03c0 polarizations",
    "`reflectivity=` takes a constant R or a table R(energy, angle)",
    "New mirror shapes: hyperbolic, parabolic, spherical, toroidal",
], l + w + 0.35, 2.15, SLIDE_W - MARGIN - (l + w + 0.35), 4.6, size=20, space=14)
slide_number(s, 4)

# ---------------------------------------------------------------- 5. surface errors
s = new_slide()
title(s, "Surface height errors and misalignment")
takeaway(s, "Height errors of 1.5 nm RMS on the KB mirrors spread part of the focus into a cross-shaped halo")
fw = 9.6
l, t, w, h = picture(s, "kb_focus.png", (SLIDE_W - fw) / 2, 1.95, width=fw)
caption(s, "A 9 keV beam: (a) at the mirrors, 50 m from the source, linear scale. (b), (c) at the focus, "
           "log scale normalized to the ideal peak: (b) ideal mirrors, (c) 1.5 nm RMS height errors on both mirrors.",
        l, t + h + 0.02, w, 0.45, size=12)
cols = [
    ("`Surface_errors`", "1D profile or 2D map of heights, placed right after its mirror"),
    ("`random_profile`", "Gaussian random surface, power-law spectrum, given RMS height"),
    ("`fixed_axis=True`", "a pitch error deflects the beam by 2\u00b7(`theta0` \u2212 `thetaE`) off the design axis"),
]
cw, cg, ct = 3.85, 0.25, t + h + 0.52
cx0 = (SLIDE_W - (3 * cw + 2 * cg)) / 2
for i, (head, body) in enumerate(cols):
    shp = rounded_box(s, cx0 + i * (cw + cg), ct, cw, SLIDE_H - ct - 0.5, LIGHT)
    tf = shp.text_frame
    tf.vertical_anchor = MSO_ANCHOR.TOP
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    p.space_after = Pt(2)
    add_runs(p, head, 18, bold=True)
    p = tf.add_paragraph()
    p.alignment = PP_ALIGN.LEFT
    add_runs(p, body, 18)
slide_number(s, 5)

# ---------------------------------------------------------------- 6. numerics
s = new_slide()
title(s, "Numerics under control")
takeaway(s, "By default the beam wraps around the window, the suggested window gives the analytic profile")
l, t, w, h = picture(s, "sampling_suggest.png", MARGIN, 2.0, width=6.6)
caption(s, "A 9 keV Gaussian beam after a 30 m drift, intensity summed over y. Default parameters (blue): "
           "RMS size 11.5 % below the analytic one (dashed). Suggested window (orange): 3.38 times wider.",
        l, t + h + 0.05, w, 0.6)
bullets(s, [
    "`check_sampling` measures the phase step between pixels and the intensity at the window border",
    "`suggest_propagation_parameters` chooses `Use_PP` zoom and sampling by trial propagations",
    "`set_radius_from_field` corrects the radius of curvature near the waist",
], l + w + 0.35, 2.15, SLIDE_W - MARGIN - (l + w + 0.35), 4.6, size=20, space=14)
slide_number(s, 6)

# ---------------------------------------------------------------- 7. speed
s = new_slide()
title(s, "Speed: SRW with OpenMP or CUDA")
takeaway(s, "OpenMP on 32 cores or CUDA on an RTX 3090 run the drift about 3 times faster than the serial wheel")
l, t, w, h = picture(s, "srw_builds.png", MARGIN, 2.0, width=6.6)
caption(s, "A 5 m drift of a 5 keV pulse of 1024 \u00d7 1024 points and 128 slices. Grey: serial builds, "
           "blue: OpenMP on 32 cores, green: CUDA on an RTX 3090. Linux server, median of 3 runs.",
        l, t + h + 0.05, w, 0.6)
rx = l + w + 0.35
rw = SLIDE_W - MARGIN - rx
_, tf = textbox(s, rx, 2.1, rw, 0.5)
add_runs(tf.paragraphs[0], "Build from the latest SRW sources, Linux:", 20, DARK, bold=True)
code_block(s, ["make install SRW_MODE=omp", "make install SRW_MODE=cuda"], rx, 2.7, rw, 0.95, size=18)
bullets(s, [
    "Any SRW branch, tag or commit: `SRW_REF=...`",
    "GPU on request with `gpu=1`, experimental",
    "Windows uses the PyPI wheel",
], rx, 3.95, rw, 2.8, size=20, space=12)
slide_number(s, 7)

# ---------------------------------------------------------------- 8. new in 2026.11
s = new_slide()
title(s, "New in 2026.11")
takeaway(s, "A modernized WPG with new mirror physics and tools that check the numerics")
bullets(s, [
    "SRW from the official `srwpy` package instead of the copy of 2019 bundled with WPG",
    "Current Python 3.10+, NumPy 2, h5py 3 and Matplotlib 3, installed from `pyproject.toml`",
    "Mirror reflectivity, new mirror shapes, surface height errors, pitch errors that deflect the beam",
    "Sampling reports, automatic choice of the mesh, SRW builds with OpenMP or CUDA",
    "About 150 unit tests and 7 sample notebooks run by GitHub Actions on Linux and Windows",
], MARGIN, 2.1, SLIDE_W - 2 * MARGIN, 3.6, size=22, space=20)
note = rounded_box(s, MARGIN, 5.75, SLIDE_W - 2 * MARGIN, 0.95, LIGHT, ACCENT)
tf = note.text_frame
tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p = tf.paragraphs[0]
add_runs(p, "Fixes of beam parameters, converters and optics that gave wrong results. ", 20, ACCENT, bold=True)
add_runs(p, "Results of old scripts can change: see the upgrade notes in `CHANGELOG.md`.", 20)
slide_number(s, 8)

# ---------------------------------------------------------------- 9. get started + cite
s = new_slide()
title(s, "Get started and cite")
_, tf = textbox(s, MARGIN, 1.45, 6, 0.5)
add_runs(tf.paragraphs[0], "Install (Python 3.10 or newer)", 22, ACCENT, bold=True)
code_block(s, [
    "git clone https://github.com/samoylv/WPG.git",
    "cd WPG",
    'python -m pip install -e ".[test,materials]" --config-settings editable_mode=compat',
    "python -m pytest",
], MARGIN, 2.0, SLIDE_W - 2 * MARGIN, 1.45, size=16)

_, tf = textbox(s, MARGIN, 3.6, SLIDE_W - 2 * MARGIN, 1.0)
p = tf.paragraphs[0]
add_runs(p, "Documentation:  ", 22, ACCENT, bold=True)
add_runs(p, DOCS_URL, 22, DARK)
link(p.runs[-1], DOCS_URL)
p = tf.add_paragraph()
add_runs(p, "Code and tutorials:  ", 22, ACCENT, bold=True)
add_runs(p, REPO_URL, 22, DARK)
link(p.runs[-1], REPO_URL)

cite = rounded_box(s, MARGIN, 4.85, SLIDE_W - 2 * MARGIN, 1.85, LIGHT)
tf = cite.text_frame
tf.vertical_anchor = MSO_ANCHOR.MIDDLE
for side in ("left", "right"):
    setattr(tf, f"margin_{side}", Inches(0.25))
p = tf.paragraphs[0]
p.space_after = Pt(6)
p.alignment = PP_ALIGN.LEFT
add_runs(p, "If you use WPG, please cite", 20, ACCENT, bold=True)
p = tf.add_paragraph()
p.alignment = PP_ALIGN.LEFT
add_runs(p, "L. Samoylova, A. Buzmakov, O. Chubar, H. Sinn. ", 18)
add_runs(p, "WavePropaGator: interactive framework for X-ray free-electron laser optics design and simulations. ",
         18, italic=True)
add_runs(p, "J. Appl. Cryst. ", 18)
add_runs(p, "49", 18, bold=True)
add_runs(p, ", 1347\u20131355 (2016). ", 18)
add_runs(p, "doi:10.1107/S160057671600995X", 18, ACCENT)
link(p.runs[-1], DOI_URL)
slide_number(s, 9)

# ---------------------------------------------------------------- save and export
prs.save(OUT)
sys.path.insert(0, KEYNOTE_FIX_DIR)
from keynote_fix import fix  # noqa: E402

fix(str(OUT))
print("saved", OUT, len(prs.slides), "slides")

if "--no-pdf" not in sys.argv:
    subprocess.run([SOFFICE, "--headless", "--convert-to", "pdf", "--outdir", str(HERE), str(OUT)],
                   check=True, timeout=180)
    pdf = OUT.with_suffix(".pdf")
    import fitz  # PyMuPDF

    doc = fitz.open(pdf)
    for page in doc:
        r = page.rect
        print(f"page {page.number + 1}: {r.width:.0f} x {r.height:.0f} pt, ratio {r.width / r.height:.4f}")
    doc[0].get_pixmap(dpi=150).save(HERE / "WPG_2026.11_title.png")
    print("saved", pdf, "and", HERE / "WPG_2026.11_title.png")
