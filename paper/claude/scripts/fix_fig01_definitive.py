"""Fix Fig01 — definitive version based on full pixel scan of bottom half."""
from PIL import Image, ImageDraw, ImageFont
import numpy as np

base = "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
src = f"{base}/plot/Figures/Fig01/Fig01_workflow.png"
im = Image.open(src).convert("RGB")
arr = np.array(im)
W, H = im.size  # 1536x1024
draw = ImageDraw.Draw(im)

def F(bold=False, size=12):
    p = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
    try: return ImageFont.truetype(p, size)
    except: return ImageFont.load_default()

def bg(y, x, r=5):
    region = arr[max(0,y-r):y+r+1, max(0,x-r):x+r+1]
    light = region[region.sum(axis=2) > 650]
    return tuple(light.mean(axis=0).astype(int)) if len(light) > 0 else (253,253,253)

# ============================================================
# Fix 1: Top annotation — "(2-5 basins per disease)" → "(1-4)"
# The text "Energy landscape with basins (2-5 basins per disease)"
# is ONE line at y=415-424, x=858-1243
# "(2-5" starts at ~x=1068 (after "Energy landscape with basins ")
# ============================================================
print("Fix 1: Top annotation")
bg1 = bg(420, 1150)
# Cover from x=1063 to end of text
draw.rectangle([1063, 413, 1246, 425], fill=bg1)
draw.text((1064, 415), "(1–4 basins per disease)", fill=(40, 40, 40), font=F(size=8))

# ============================================================
# Fix 2: "Basin structure (2-5 per disease)"
# y=510-521, x=580-742 — one text block in the green summary box
# Need to change "2-5" to "1-4"
# "(2-5" starts at ~x=660 (after "Basin structure ")
# ============================================================
print("Fix 2: Basin structure box")
bg2 = bg(516, 680)
# Cover from "(" onwards: x=658 to end
draw.rectangle([656, 508, 745, 523], fill=bg2)
draw.text((658, 510), "(1–4 per disease)", fill=(40, 40, 40), font=F(size=8))

# ============================================================
# Fix 3: "Mortality gap (up to 20.4%)" → "(up to 16.4%)"
# Line 1 "Mortality gap" at y=618-630 — keep
# Line 2 "(up to 20.4%)" at y=640-653, x=595-730
# ============================================================
print("Fix 3: Mortality gap")
bg3 = bg(648, 660)
draw.rectangle([593, 638, 732, 656], fill=bg3)
draw.text((596, 640), "(up to 16.4%)", fill=(80, 80, 80), font=F(size=9))

# ============================================================
# Fix 4-5: Cox validation text
# The Cox text spans a LARGE area:
#   "Per-basin Mortality:" at y=750-758, x=571-752
#   "HR < 1 in 33 / 51" at y=770-782, x=573-751
#   [gap]
#   "Per-basin Hospitalization:" at y=812-826, x=567-754
# Cover the entire Cox text area and redraw
# ============================================================
print("Fix 4-5: Cox validation text")
# Sample bg from the Cox box area
cox_bg = bg(760, 650)
# Cover ALL Cox text: y=744-830, x=560-755
draw.rectangle([558, 744, 756, 830], fill=cox_bg)

# Redraw Per-basin Mortality
draw.text((563, 746), "Per-basin Mortality:", fill=(40, 40, 40), font=F(bold=True, size=9))
draw.text((563, 760), "HR < 1 in", fill=(40, 40, 40), font=F(bold=True, size=11))
draw.text((653, 760), "27 / 51", fill=(20, 60, 140), font=F(bold=True, size=13))

# Redraw Per-basin Hospitalization
draw.text((563, 790), "Per-basin Hospitalization:", fill=(40, 40, 40), font=F(bold=True, size=9))
draw.text((563, 804), "HR < 1 in", fill=(40, 40, 40), font=F(bold=True, size=11))
draw.text((653, 804), "35 / 51", fill=(20, 60, 140), font=F(bold=True, size=13))

# ============================================================
# Fix 6: Forest plot — cover ENTIRE right-side area
# "Cross-disease analysis (51 disease groups)" header at y=486-498
# Header text: "BasinNA bidirectionality" at y=506-525, x=739-947
# Forest plot content: y=594-970, x=680-1080
# Frailty codes: y=586-940, x=1080-1536
# Cover everything from y=504 down to y=975, x=680-1534
# ============================================================
print("Fix 6-7: Forest plot + frailty codes")

# --- FOREST PLOT (left half of right section) ---
forest_bg = (252, 253, 254)
draw.rectangle([680, 504, 1080, 975], fill=forest_bg)
draw.rectangle([680, 504, 1080, 975], outline=(190, 200, 215))

y = 510
draw.text((694, y), "Per-basin majority criterion", fill=(40, 40, 40), font=F(bold=True, size=11))
y += 16
draw.text((694, y), "(overall HR, landmark design)", fill=(110, 110, 110), font=F(size=8))
y += 16

draw.line([(694, y), (1066, y)], fill=(200, 200, 210), width=1)
y += 5
draw.text((698, y), "Disease", fill=(80, 80, 80), font=F(bold=True, size=8))
draw.text((845, y), "HR (95% CI)", fill=(80, 80, 80), font=F(bold=True, size=8))
y += 12
draw.line([(694, y), (1066, y)], fill=(220, 220, 225), width=1)
y += 4

pass_list = [
    ("Hypertension", "0.91", "0.90–0.92"),
    ("CKD", "0.88", "0.86–0.89"),
    ("Diabetes", "0.95", "0.94–0.96"),
    ("Heart failure", "0.93", "0.91–0.94"),
    ("Stroke", "0.94", "0.93–0.95"),
    ("AF", "0.94", "0.92–0.95"),
    ("Stomach ca.", "0.87", "0.84–0.89"),
    ("Sleep dis.", "0.90", "0.88–0.91"),
    ("Asp. pneumonia", "0.98", "0.97–0.99"),
]
fail_list = [
    ("Depression", "1.17", "1.15–1.19"),
    ("RA", "1.22", "1.19–1.24"),
    ("Thyroid dis.", "1.05", "1.04–1.06"),
]

for name, hr, ci in pass_list:
    draw.ellipse([698, y+3, 703, y+8], fill=(33, 100, 200))
    draw.text((708, y), name, fill=(40, 40, 40), font=F(size=8))
    draw.text((835, y), hr, fill=(33, 100, 200), font=F(bold=True, size=8))
    draw.text((870, y), f"({ci})", fill=(120, 120, 120), font=F(size=7))
    y += 12

y += 3
draw.line([(694, y), (1066, y)], fill=(220, 220, 225), width=1)
y += 4

for name, hr, ci in fail_list:
    draw.ellipse([698, y+3, 703, y+8], fill=(200, 50, 50))
    draw.text((708, y), name, fill=(40, 40, 40), font=F(size=8))
    draw.text((835, y), hr, fill=(200, 50, 50), font=F(bold=True, size=8))
    draw.text((870, y), f"({ci})", fill=(120, 120, 120), font=F(size=7))
    y += 12

y += 6
draw.line([(694, y), (1066, y)], fill=(200, 200, 210), width=1)
y += 8

# Summary totals
draw.text((694, y), "Mortality:", fill=(40, 40, 40), font=F(bold=True, size=9))
draw.text((775, y), "27/51 PASS (53%)", fill=(33, 100, 200), font=F(bold=True, size=9))
y += 14
draw.text((694, y), "Hospitalization:", fill=(40, 40, 40), font=F(bold=True, size=9))
draw.text((810, y), "35/51 PASS (69%)", fill=(33, 100, 200), font=F(bold=True, size=9))
y += 18

# Legend
draw.ellipse([694, y+2, 699, y+7], fill=(33, 100, 200))
draw.text((704, y), "PASS (HR < 1)", fill=(33, 100, 200), font=F(size=7))
draw.ellipse([810, y+2, 815, y+7], fill=(200, 50, 50))
draw.text((820, y), "FAIL (HR ≥ 1)", fill=(200, 50, 50), font=F(size=7))
y += 12
draw.text((694, y), "Hazard ratio (log scale)", fill=(120, 120, 120), font=F(size=7))

# --- FRAILTY CODES (right half of right section) ---
frailty_bg = (253, 249, 252)
draw.rectangle([1080, 504, 1534, 975], fill=frailty_bg)
draw.rectangle([1080, 504, 1534, 975], outline=(200, 185, 195))

y = 510
draw.text((1094, y), "Basin comorbidity examples", fill=(40, 40, 40), font=F(bold=True, size=11))
y += 16
draw.text((1094, y), "(landmark design)", fill=(110, 110, 110), font=F(size=8))
y += 20

# Thyroid B3
draw.text((1094, y), "Thyroid B3", fill=(0, 70, 150), font=F(bold=True, size=9))
draw.text((1185, y), "(death 13.0%)", fill=(100, 100, 100), font=F(size=8))
y += 14
draw.text((1104, y), "CKD stage 5, hyperkalemia", fill=(60, 60, 60), font=F(size=8))
y += 12
draw.text((1104, y), "Renal failure, anemia", fill=(60, 60, 60), font=F(size=8))
y += 18

# Thyroid B4
draw.text((1094, y), "Thyroid B4", fill=(0, 70, 150), font=F(bold=True, size=9))
draw.text((1185, y), "(death 12.8%)", fill=(100, 100, 100), font=F(size=8))
y += 14
draw.text((1104, y), "Periodontitis, dental caries", fill=(60, 60, 60), font=F(size=8))
y += 12
draw.text((1104, y), "Hypertension, constipation", fill=(60, 60, 60), font=F(size=8))
y += 20

# Dementia B1
draw.text((1094, y), "Dementia B1", fill=(180, 40, 40), font=F(bold=True, size=9))
draw.text((1200, y), "(death 33.6%)", fill=(100, 100, 100), font=F(size=8))
y += 14
draw.text((1104, y), "Constipation, xerosis cutis", fill=(60, 60, 60), font=F(size=8))
y += 12
draw.text((1104, y), "Heart failure, dermatitis", fill=(60, 60, 60), font=F(size=8))
y += 18

# Dementia B2
draw.text((1094, y), "Dementia B2", fill=(180, 40, 40), font=F(bold=True, size=9))
draw.text((1200, y), "(death 17.2%)", fill=(100, 100, 100), font=F(size=8))
y += 14
draw.text((1104, y), "Hypertension, diabetes", fill=(60, 60, 60), font=F(size=8))
y += 12
draw.text((1104, y), "Hyperlipidemia", fill=(60, 60, 60), font=F(size=8))
y += 20

# Sepsis B3
draw.text((1094, y), "Sepsis B3", fill=(160, 80, 0), font=F(bold=True, size=9))
draw.text((1175, y), "(death 45.3%)", fill=(100, 100, 100), font=F(size=8))
y += 14
draw.text((1104, y), "Highest mortality basin", fill=(60, 60, 60), font=F(size=8))
y += 12
draw.text((1104, y), "B1: 31.4%, B2: 31.0%", fill=(60, 60, 60), font=F(size=8))
y += 20

# Malnutrition B1
draw.text((1094, y), "Malnutrition B1", fill=(100, 120, 0), font=F(bold=True, size=9))
draw.text((1220, y), "(death 44.3%)", fill=(100, 100, 100), font=F(size=8))
y += 14
draw.text((1104, y), "Entropion, presbyopia", fill=(60, 60, 60), font=F(size=8))
y += 12
draw.text((1104, y), "Macular degeneration", fill=(60, 60, 60), font=F(size=8))
y += 20

# Diabetes B1
draw.text((1094, y), "Diabetes B1", fill=(150, 0, 100), font=F(bold=True, size=9))
draw.text((1200, y), "(death 9.1%)", fill=(100, 100, 100), font=F(size=8))
y += 14
draw.text((1104, y), "Tooth loss, dental caries", fill=(60, 60, 60), font=F(size=8))
y += 12
draw.text((1104, y), "Periodontitis, keratitis", fill=(60, 60, 60), font=F(size=8))

print("All fixes applied")
im.save(src)
print(f"Saved: {src}")
