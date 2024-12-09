from PIL import Image
import astroalign as aa
source = Image.open(r"D:\space_photos2024\pk111_Ha\haoiimaybe.tif")
target = Image.open(r"D:\space_photos2024\pk111_Ha\rgbdone.tif")
registered, footprint = aa.register(source, target)
# Convert back to pillow image if necessary:
registered = Image.fromarray(registered.astype("unit8"))
