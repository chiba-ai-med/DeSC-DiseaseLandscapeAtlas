"""Build composite figures from individual panels using Pillow."""
import sys
sys.path.insert(0, '')
from PIL import Image, ImageDraw, ImageFont

base = "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"

def make_composite(panels, title, outpath, layout="2x2"):
    """Compose 4 panels into a single figure with a title and panel labels."""
    imgs = []
    for path in panels:
        im = Image.open(path).convert("RGB")
        imgs.append(im)

    # Target: 4800x3600 at 300 dpi
    W, H = 4800, 3600
    title_h = 140
    panel_w = W // 2
    panel_h = (H - title_h) // 2

    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)

    # Title
    try:
        font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 48)
        font_label = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 42)
    except:
        font_title = ImageFont.load_default()
        font_label = font_title

    bbox = draw.textbbox((0, 0), title, font=font_title)
    tw = bbox[2] - bbox[0]
    draw.text(((W - tw) // 2, 40), title, fill="black", font=font_title)

    labels = ["a", "b", "c", "d"]
    positions = [
        (0, title_h),
        (panel_w, title_h),
        (0, title_h + panel_h),
        (panel_w, title_h + panel_h),
    ]

    for i, (img, (px, py)) in enumerate(zip(imgs, positions)):
        resized = img.resize((panel_w, panel_h), Image.LANCZOS)
        canvas.paste(resized, (px, py))
        draw.text((px + 20, py + 10), labels[i], fill="black", font=font_label)

    canvas.save(outpath, dpi=(300, 300))
    print(f"  Written: {outpath}")

# Thyroid composite
print("=== Fig06: Thyroid composite ===")
make_composite(
    [
        f"{base}/plot/Figures/Fig06/Fig06a_thyroid_basin3_vs_4.png",
        f"{base}/plot/Figures/Fig06/Fig06b_thyroid_km_inter_basin.png",
        f"{base}/plot/Figures/Fig06/Fig06c_thyroid_disorder_discon.png",
        f"{base}/plot/Figures/Fig06/Fig06d_thyroid_disorder_sankey.png",
    ],
    "Thyroid disorder (Structured reversed)",
    f"{base}/plot/Figures/Fig06/Fig06_thyroid.png",
)

# Diabetes composite
print("=== FigS7: Diabetes composite ===")
make_composite(
    [
        f"{base}/plot/Figures/Supplementary/Fig06a_diabetes_basin1_vs_3.png",
        f"{base}/plot/Figures/Supplementary/Fig06b_diabetes_km_inter_basin.png",
        f"{base}/plot/Figures/Supplementary/Fig06c_diabetes_discon.png",
        f"{base}/plot/Figures/Supplementary/Fig06d_diabetes_sankey.png",
    ],
    "Diabetes (Progressive)",
    f"{base}/plot/Figures/Supplementary/Fig06_diabetes.png",
)

print("Done.")
