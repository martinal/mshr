import Image

input_image = Image.open("Icons.png")

for i, name in enumerate(["polygon", "rectangle", "circle", "ellipse", "disk", "cylinder", "box", "cone", "sphere", "tetrahedron"]) :
    x = (i%5)*268 + 5
    y = (i/5)*268
    print (x, y)
    cropped_image = input_image.crop((x, y, x+250, y+250))
    cropped_image.save(name+"-large.png")
    small_image = cropped_image.resize( (32, 32), Image.ANTIALIAS)
    small_image.save(name+"-small.png")
    
