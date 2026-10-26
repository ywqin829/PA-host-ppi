from __future__ import annotations

from pathlib import Path
import shutil

import fitz


root = Path(__file__).resolve().parents[2]
original_dir = root / "figures" / "original"
panel_dir = root / "figures"
final_dir = root / "figures" / "final"
final_dir.mkdir(parents=True, exist_ok=True)


def place_pdf(
    page: fitz.Page,
    panel_path: Path,
    target: fitz.Rect,
    padding: float = 3.0,
    panel_label: str | None = None,
) -> None:
    page.add_redact_annot(
        target,
        fill=(1, 1, 1),
        cross_out=False,
    )
    page.apply_redactions(
        images=fitz.PDF_REDACT_IMAGE_NONE,
        graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
        text=fitz.PDF_REDACT_TEXT_REMOVE,
    )
    inner = fitz.Rect(
        target.x0 + padding,
        target.y0 + padding,
        target.x1 - padding,
        target.y1 - padding,
    )
    with fitz.open(panel_path) as panel:
        page.show_pdf_page(
            inner,
            panel,
            0,
            keep_proportion=True,
        )
    if panel_label:
        page.insert_text(
            fitz.Point(target.x0 + 3, target.y0 + 11),
            panel_label,
            fontsize=9,
            fontname="Helvetica-Bold",
            color=(0, 0, 0),
        )


# Figure 4: keep the original Adobe composition except the two panels whose
# underlying statistics changed (descriptive single-cell panel and corrected
# TOST panel).
fig4 = fitz.open(original_dir / "Fig4.pdf")
page4 = fig4[0]
place_pdf(
    page4,
    panel_dir / "Fig4D_single_cell_descriptive.pdf",
    fitz.Rect(292, 161, 594, 334),
    panel_label="D",
)
place_pdf(
    page4,
    panel_dir / "Fig4F_updated_v9.pdf",
    fitz.Rect(292, 334, 594, 502),
)
fig4.set_metadata(
    {
        "title": "PA-HIR Figure 4",
        "author": "Jianhao Xu; Cuiqi Yang; Liangqiong Chen; Jian Zhang; Jiahui Tang; Yongwei Qin",
        "subject": "Mechanistic interrogation and descriptive cellular localization",
        "keywords": "Pseudomonas aeruginosa; TLR5; single-cell RNA-seq; TOST",
    }
)
fig4.save(final_dir / "Fig4.pdf", garbage=4, deflate=True)
fig4.close()


# Figure 5: preserve the original panel E and replace A-D with the corrected
# nine-dataset meta-analysis, panel analysis, leave-one-out and sensitivity
# panels.
fig5 = fitz.open(original_dir / "Fig5.pdf")
page5 = fig5[0]
place_pdf(page5, panel_dir / "Fig5A_updated_v9.pdf", fitz.Rect(0, 0, 595.3, 259))
place_pdf(page5, panel_dir / "Fig5B_updated_v9.pdf", fitz.Rect(0, 260, 257, 447))
place_pdf(page5, panel_dir / "Fig5C_updated_v9.pdf", fitz.Rect(258, 260, 595.3, 447))
place_pdf(page5, panel_dir / "Fig5D_updated_v9.pdf", fitz.Rect(0, 448, 257, 655))
fig5.set_metadata(
    {
        "title": "PA-HIR Figure 5",
        "author": "Jianhao Xu; Cuiqi Yang; Liangqiong Chen; Jian Zhang; Jiahui Tang; Yongwei Qin",
        "subject": "Corrected TLR5 meta-analysis and receptor/adaptor panel",
        "keywords": "Pseudomonas aeruginosa; TLR5; meta-analysis; receptor panel",
    }
)
fig5.save(final_dir / "Fig5.pdf", garbage=4, deflate=True)
fig5.close()


for source, destination in [
    (original_dir / "Fig1.pdf", final_dir / "Fig1.pdf"),
    (original_dir / "Fig2.pdf", final_dir / "Fig2.pdf"),
    (original_dir / "Fig3.pdf", final_dir / "Fig3.pdf"),
    (
        panel_dir / "Fig6_state_decoupling.pdf",
        final_dir / "Fig6.pdf",
    ),
]:
    shutil.copy2(source, destination)


for number in range(1, 7):
    document = fitz.open(final_dir / f"Fig{number}.pdf")
    pixmap = document[0].get_pixmap(
        matrix=fitz.Matrix(150 / 72, 150 / 72),
        alpha=False,
    )
    pixmap.save(final_dir / f"Fig{number}_preview.png")
    document.close()

print(f"Wrote final figure PDFs and previews to {final_dir}")
