# -*- coding:utf-8 -*-
from Tkinter import *

def menu(action, config):
    frame = Tk()
    frame.geometry("550x400")


    #CHOICE OF THE SEQUENCE
    seq = IntVar()
    """
    choose_seq = Listbox(frame, width=15)
    choose_seq.grid(row=10)
    sequences = ["AGCT", "ATGC", "CGCG", "CGGC", "CGTA", "TAAT", "TATA", "TCGA"]
    for element in sequences:
      choose_seq.insert(END, element)
    def clic(inutile): #inutile is a useless but necessary variable
      seq.set(choose_seq.curselection()[0])
      print seq.get()
    choose_seq.bind('<ButtonRelease-1>', clic)
    ###
    """

    #CHOICE OF EIGENVALUE 1
    texte_vp1 = Label(frame, text='Eigenvalue 1 (0 for equilibrium) :', fg="black")
    vp1 = IntVar()
    vp1.set(206)
    get_vp1 = Entry(textvariable=vp1, width=8)
    texte_vp1.grid(row=20, sticky=W, padx=10, pady=5)
    get_vp1.grid(row=20,column=1)

    #CHOICE OF TRUNCATION 1
    texte_m1 = Label(frame, text='Truncation filter : ', fg="black")
    m1 =IntVar()
    m1.set(0)
    get_m11 = Radiobutton(frame, text="No",  variable=m1, value=0)
    get_m12 = Radiobutton(frame, text="M66", variable=m1, value=1)
    get_m13 = Radiobutton(frame, text="M55", variable=m1, value=2)
    get_m14 = Radiobutton(frame, text="M53", variable=m1, value=3)
    texte_m1.grid(row=30, sticky=W, padx=50, pady=5)
    get_m11.grid(row=30, column=1, sticky=W, padx=10)
    get_m12.grid(row=30, column=2, sticky=W, padx=10)
    get_m13.grid(row=30, column=3, sticky=W, padx=10)
    get_m14.grid(row=30, column=4, sticky=W, padx=10)

    #CHOICE OF EIGENVALUE 2
    texte_vp2 = Label(frame, text='Eigenvalue 2 (optional) : ', fg="black")
    vp2 = IntVar()
    vp2.set(0)
    get_vp2 = Entry(textvariable=vp2, width=8)
    texte_vp2.grid(row=40, sticky=W, padx=10, pady=5)
    get_vp2.grid(row=40,column=1)

    #CHOICE OF TRUNCATION 2
    texte_m2 = Label(frame, text='Truncation filter : ', fg="black")
    m2 =IntVar()
    m2.set(0)
    get_m21 = Radiobutton(frame, text="No",  variable=m2, value=0)
    get_m22 = Radiobutton(frame, text="M66", variable=m2, value=1)
    get_m23 = Radiobutton(frame, text="M55", variable=m2, value=2)
    get_m24 = Radiobutton(frame, text="M53", variable=m2, value=3)
    texte_m2.grid(row=50, sticky=W, padx=50, pady=5)
    get_m21.grid(row=50, column=1, sticky=W, padx=10)
    get_m22.grid(row=50, column=2, sticky=W, padx=10)
    get_m23.grid(row=50, column=3, sticky=W, padx=10)
    get_m24.grid(row=50, column=4, sticky=W, padx=10)

    #SPEED OF THE ANIMATION
    texte_precision = Label(frame, text='Precision (lowers speed) :', fg="black")
    precision = IntVar()
    precision.set(10)
    get_precision = Entry(textvariable=precision, width=8)
    texte_precision.grid(row=60, sticky=W, padx=10, pady=5)
    get_precision.grid(row=60,column=1)

    #SCALAR MULTIPLE OF THE EIGENVECTOR
    texte_mult = Label(frame, text='Lambda (mult. of the eigenvector) :', fg="black")
    mult = DoubleVar()
    mult.set(2)
    get_mult = Entry(textvariable=mult, width=8)
    texte_mult.grid(row=65, sticky=W, padx=10, pady=5)
    get_mult.grid(row=65,column=1)

    #SHOW BASES/CENTERLINE OR NOT
    texte_affichage = Label(frame, text='Show ')
    centerline = IntVar(); showbases = IntVar(); ghost = IntVar();
    centerline.set(0); showbases.set(1); ghost.set(1);
    get_showbases = Checkbutton(frame, text="bases", variable=showbases)
    get_centerline = Checkbutton(frame, text="center line", variable=centerline)
    get_ghost = Checkbutton(frame, text="ghost", variable=ghost)
    texte_affichage.grid(row=70, sticky=W, padx=10, pady=5)
    get_showbases.grid(row=70, column=1, sticky=W, padx=10)
    get_centerline.grid(row=70, column=2, sticky=W, padx=10)
    get_ghost.grid(row=70, column=3, sticky=W, padx=10)

    #MAIN BUTTONS
    quit_button  = Button(frame,text="Quit",width=15,command=frame.quit)
    start_button = Button(frame,text="Draw",width=15,command = lambda: \
                          action(config,vp1.get(),vp2.get(),precision.get(), \
                                 m1.get(),m2.get(),centerline.get(),showbases.get(), \
                                 seq.get(),mult.get(),ghost.get() ))
    start_button.grid(row=90, columnspan=5, pady=15)
    quit_button.grid(row=100, columnspan=5)


    frame.mainloop()

