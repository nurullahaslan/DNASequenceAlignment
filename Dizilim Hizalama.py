from tkinter import *
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from Bio.SubsMat import MatrixInfo
def maximum(a, b, c): 
    list = [a, b, c]
    return 0 if a<0 and b<0 and c<0 else max(list) 
        
def findAlignment(m,s,seq):
    for i in range(len(s)):
        if s[i] == m:
            seq.configure(bg="yellow")
            
def drawMatrix(window2 , s , numMatrix):
    print(numMatrix)
    for i in range(len(s)):
        for j in range(len(s[0])):
            if i==0 and j == 0:
                seq = Entry(window2)
                seq.insert(END , str(""))
                seq.grid(row=0 , column=0)
            elif i == 0 and j > 0:
                seq = Entry(window2)
                seq.insert(END , str(s[i][j]))
                seq.grid(row=i , column=j)
            elif i > 0:
                if j >= 0:
                    seq = Entry(window2)
                    seq.insert(END , str(s[i][j]))
                    seq.grid(row=i , column=j)
            findAlignment([i,j],numMatrix,seq)
def score_2():
    blosum62 = MatrixInfo.blosum62

    seq1 = firstSeq.get()
    seq1 = seq1.upper()
    seq2 = secondSeq.get()
    seq2 = seq2.upper()
    gap_pen = int(Gap.get())
    u1 = len(seq1)
    u2 = len(seq2)
    score = 0
    matrix = []
    numMatrix = [[]]
    
    for i in range(u2+1):
        matrix.append([])
        if i == 0:
            matrix[i].append('+')
            for j in range(u1):
                matrix[i].append(seq1[j])
        else:
            matrix[i].append(seq2[i-1])
            for j in range(u1):
                matrix[i].append(0)
    
    for r,row in enumerate(matrix):
        for c,col in enumerate(row):
            if r>0 and c>0:
                A = matrix[r][0] 
                B = matrix[0][c]
                matrix[r][c] = blosum62[(A,B)] if (A,B) in blosum62 else blosum62[(B,A)]
                
        
    d = maximum(u1,u2,0)
    r = u2
    c = u1
    g = 0
    for i in range(d):
        if matrix[0][c] == matrix[r][0] or r==c:    
            score +=matrix[r][c]
            numMatrix.append([r , c])
            c-=1
            r-=1
        else:
            g +=1
            c-=1
    score += g*gap_pen
    window2 = Tk()
    drawMatrix(window2 , matrix , numMatrix)
    return score

def score_3(seq1,seq2,gap_pen,mismatch,match):
    seq1 = firstSeq.get()
    seq1 = seq1.upper()
    seq2 = secondSeq.get()
    seq2 = seq2.upper()
    gap_pen = int(Gap.get())
    mismatch = int(Mismatch.get())
    match = int(Match.get())
    u1 = len(seq1)
    u2 = len(seq2)
    matrix = []
    i = 2
    biggest_number = 0
    r_of_bg = 0
    c_of_bg = 0
    
    for i in range(u2+2):
        matrix.append([])
        if i == 0:
            matrix[i].append('+')
            for j in range(u1+1):
                if j == 0:
                    matrix[i].append('-')
                else:
                    matrix[i].append(seq1[j-1])
        elif i== 1:
            matrix[i].append('-')
            for j in range(u1+1):
                matrix[i].append(0)
        else:
            matrix[i].append(seq2[i-2])
            for j in range(u1+1):
                matrix[i].append(0)
    
    for r,row in enumerate(matrix):
        for c,col in enumerate(row):
            first_row = matrix[r][0] 
            first_col = matrix[0][c]
            
            if c>1 and r>1:
                up = matrix[r-1][c]
                left = matrix[r][c-1]
                upleft =matrix[r-1][c-1]
                
                if first_row == first_col:
                    matrix[r][c] = upleft + match
                    if matrix[r][c] > biggest_number:
                        biggest_number = matrix[r][c]
                        r_of_bg = r
                        c_of_bg = c
                else:
                    u = up+gap_pen
                    l = left+gap_pen
                    ul = upleft + mismatch
                    
                    matrix[r][c] = maximum(u,l,ul)
    
    return matrix,biggest_number,r_of_bg,c_of_bg
def dot_matrix():
    """Create a dot matrix of the two DNA sequences specified by the user"""
    window = 1
    seq_one = firstSeq.get()
    seq_two = secondSeq.get()
    data = [[(seq_one[i:i + window] != seq_two[j:j + window]) \
        for j in range(len(seq_one) - window + 1)] \
        for i in range(len(seq_two) - window + 1)]
    
    plt.figure()
    plt.imshow(data, cmap='gray', vmin=0, vmax=1)

def Local_Alignment():
    numMatrix =[[]]
    seq1 = firstSeq.get()
    seq1 = seq1.upper()
    seq2 = secondSeq.get()
    seq2 = seq2.upper()
    gap_pen = int(Gap.get())
    mismatch = int(Mismatch.get())
    match = int(Match.get())
    matrix, bg,r,c = score_3(seq1,seq2,gap_pen,mismatch,match)
    window2 = Tk()
    
    numMatrix.append([r , c])
    listseq1 = []
    listseq2 = []
    
    m=0
    mm = 0
    g=0
    
    for i in range(r):
        last = matrix[r][c]
        up = matrix[r-1][c]
        left = matrix[r][c-1]
        upleft =matrix[r-1][c-1]
        
        if matrix[0][c] == matrix[r][0]:
            listseq1.append(matrix[0][c])
            listseq2.append(matrix[r][0])
            numMatrix.append([r , c])
            r -=1
            c -=1
            m +=1
        else:
            if last-upleft == mismatch:
                listseq1.append(matrix[0][c])
                listseq2.append(matrix[r][0])
                numMatrix.append([r , c])
                r -=1
                c -=1
                mm +=1
            elif last-up == gap_pen:
                listseq1.append('-')
                listseq2.append(matrix[r][0])
                numMatrix.append([r , c])
                r -=1
                g +=1
            elif last-left == gap_pen:
                listseq1.append(matrix[0][c])
                listseq2.append('-')
                numMatrix.append([r , c])
                c -=1
                g +=1
    drawMatrix(window2 , matrix , numMatrix)






root = Tk()
root.title ('Dizilim Hizalama')
frame = ttk.Frame(root, padding = 10)


# Widgets in the control frame
label1 = Label(frame, text = 'Dizilimleri giriniz.')
label2 = Label(frame, text = 'First sequence')
label3 = Label(frame, text = 'Second sequence')
label4 = Label(frame, text = 'Match')
label5 = Label(frame, text = 'Mismatch')
label6 = Label(frame, text = 'Gap')

firstSeq = Entry(frame)
secondSeq = Entry(frame)
Match = Entry(frame)
Mismatch = Entry(frame)
Gap = Entry(frame)

button1 = Button(frame, text = 'Global', command = score_2)
button2 = Button(frame, text = 'Local', command = Local_Alignment)
button3 = Button(frame, text = 'Dot Matrix', command = dot_matrix)

frame.grid(column = 0, row = 0)
label1.grid(column = 0, row = 0, columnspan = 2)
label2.grid(column = 0, row = 1, sticky = W)
label3.grid(column = 0, row = 2, sticky = W)
label4.grid(column = 0, row = 3, sticky = W)
label5.grid(column = 0, row = 4, sticky = W)
label6.grid(column = 0, row = 5, sticky = W)

firstSeq.grid(column = 1, row = 1)
secondSeq.grid(column = 1, row = 2)
Match.grid(column = 1, row = 3)
Mismatch.grid(column = 1, row = 4)
Gap.grid(column = 1, row = 5)

button1.grid(column = 2, row = 2)
button2.grid(column = 2, row = 3)
button3.grid(column = 2, row = 4)

for child in frame.winfo_children(): 
    child.grid_configure(padx = 3, pady = 3)

root.mainloop()