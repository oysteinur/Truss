from tkinter import *
import matplotlib as plt
from Plot import *

window = Tk()

Heading_Label = Label(window, text="Stavmodell")

Heading_Label.pack()

load = Image.open('Konstruksjon.png')
render = ImageTk.photoImage(load)

img = Label(self, image=render)
img.image = render
img.place(x=0, y=0)

window.mainloop()