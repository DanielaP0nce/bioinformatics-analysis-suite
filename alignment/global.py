import sys

def scoring_matrix(file):
	scoring_dict = {}
	with open(file, 'r') as text:
		text = text.readlines()
		nucleotides = text[0].strip().split(" ")[1:]
		for lines in text[1:]:
			lines = lines.strip().split()
			scoring_dict[lines[0]] = {nucleotides[j-1]: lines[j] for j in range(1, len(nucleotides)+1)}
	return scoring_dict



scoring_matrix('standard.m')

# txt_sequences_filepath=sys.argv[1]
# txt_scoring_filepath=sys.argv[2]
# gap_penalty=sys.argv[3]

# def global_alignment(txt_sequences_filepath,txt_scoring_filepath,gap_penalty):
#     gap_penalty=int(gap_penalty)
#     with open(txt_sequences_filepath,"r") as text,open(txt_scoring_filepath,"r") as file:
#         string_line=text.readlines()
#         n_str=string_line[1].strip().upper() #rows
#         m_str=string_line[0].strip().upper() #columns
#         m_len = len(m_str)
#         n_len=len(n_str)
#         seq_dir_matrix = [[["0"]] * (m_len+2) for _ in range(n_len+2)]

#         scoring_matrix=[]
#         for line in file:
#           scoring_matrix.append([x for x in line.split()])
#         match_dict={}
#         mismatch_dict={}

#         for j in range(1,len(scoring_matrix[0])):
#            for i in range(1,len(scoring_matrix)):
#                if scoring_matrix[0][j]== scoring_matrix[i][0]:
#                    match_dict[scoring_matrix[0][j]+scoring_matrix[i][0]]=scoring_matrix[i][j]
#                else:
#                    mismatch_dict[scoring_matrix[0][j]+scoring_matrix[i][0]]=scoring_matrix[i][j]

#         seq_dir_matrix[1][0]="-"
#         seq_dir_matrix[0][0]="-"
#         seq_dir_matrix[0][1]="-"
        
#         for row in range(1,n_len+1):
#             seq_dir_matrix[row+1][0]=n_str[row-1]
#             x=int(seq_dir_matrix[row][1][0])+gap_penalty
#             seq_dir_matrix[row+1][1]=[x,"up"]
#         for col in range(1, m_len+1):
#            seq_dir_matrix[0][col+1]=m_str[col-1]
#            y = int(seq_dir_matrix[1][col][0]) - 1
#            seq_dir_matrix[1][col+1] = [y, "left"]

#            for row in range(1,n_len+1):
#                 if n_str[row-1]==m_str[col-1]:
#                     temp_str=''
#                     temp_str=temp_str+n_str[row-1]+m_str[col-1]
                    
#                     match=int(seq_dir_matrix[row][col][0]) + int(match_dict[temp_str])
#                     indel_up = int(seq_dir_matrix[row][col+1][0]) + gap_penalty
#                     indel_left = int(seq_dir_matrix[row+1][col][0]) + gap_penalty
#                     max_score=max(match,indel_left,indel_up)
#                     if (max_score==match and max_score==indel_left and max_score==indel_up):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"diagonal,left,up"]
#                     elif(max_score==match and max_score==indel_left):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"diagonal,left"]
#                     elif(max_score==match and max_score==indel_up):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"diagonal,up"]
#                     elif (max_score==indel_left and max_score==indel_up):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"left,up"]
#                     elif(max_score==match):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"diagonal"]
#                     elif(max_score==indel_left):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"left"]
#                     else:
#                         seq_dir_matrix[row+1][col+1]=[max_score,"up"]
#                 else:
#                     temp_str=""
#                     temp_str=temp_str+ n_str[row-1]+m_str[col-1]
                    
#                     mismatch=int(seq_dir_matrix[row][col][0])+int(mismatch_dict[temp_str])
#                     indel_up = int(seq_dir_matrix[row][col+1][0]) + gap_penalty
#                     indel_left = int(seq_dir_matrix[row+1][col][0]) + gap_penalty
#                     max_score=max(mismatch,indel_up,indel_left)
#                     if (max_score==mismatch and max_score==indel_left and max_score==indel_up):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"diagonal,left,up"]
#                     elif(max_score==mismatch and max_score==indel_left):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"diagonal,left"]
#                     elif(max_score==mismatch and max_score==indel_up):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"diagonal,up"]
#                     elif (max_score==indel_left and max_score==indel_up):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"left,up"]
#                     elif(max_score==mismatch):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"diagonal"]
#                     elif(max_score==indel_left):
#                         seq_dir_matrix[row+1][col+1]=[max_score,"left"]
#                     else:
#                         seq_dir_matrix[row+1][col+1]=[max_score,"up"]
            
#         backtracking_path=[]
#         index_matrix=[]
#         max_score=seq_dir_matrix[n_len+1][m_len+1][0]
#         row = n_len+1 #x
#         col = m_len+1 #y
#         x = 0
#         y = 0
#         index_matrix=[]


#         while seq_dir_matrix[row][col] !=  ["0"]:
#             backtracking_path.append(seq_dir_matrix[row][col])
#             index_matrix.append((row,col))
#             diagonal = seq_dir_matrix[row-1][col-1][0]
#             up = seq_dir_matrix[row-1][col][0]
#             left= seq_dir_matrix[row][col-1][0]

#             if ("diagonal" in seq_dir_matrix[row][col][1]) and ("left" in seq_dir_matrix[row][col][1]) and ("up" in seq_dir_matrix[row][col][1]):
#                    if diagonal >  up and diagonal > left:   
#                        seq_dir_matrix[row][col][1]="diagonal" 
#                        x = row - 1
#                        row = x
#                        y = col - 1
#                        col = y
#                    elif(left> diagonal and left>up): 
#                        seq_dir_matrix[row][col][1]="left"
#                        y = col - 1
#                        col = y
#                    else: 
#                        seq_dir_matrix[row][col][1]="up"
#                        x=row-1
#                        row=x
#             elif("diagonal" in seq_dir_matrix[row][col][1]) and ("left" in seq_dir_matrix[row][col][1]):
#                    if (diagonal>left):
#                        seq_dir_matrix[row][col][1]="diagonal"     
#                        x = row - 1
#                        row = x
#                        y = col - 1
#                        col = y
#                    else:     
#                        seq_dir_matrix[row][col][1]="left"  
#                        y = col - 1
#                        col = y
#             elif(("diagonal" in seq_dir_matrix[row][col][1]) and ("up" in seq_dir_matrix[row][col][1])):
#                    if (diagonal>up):
#                        seq_dir_matrix[row][col][1]="diagonal"
#                        x = row-1
#                        row=x
#                        y=col-1
#                        col=y           
#                    else:
#                        seq_dir_matrix[row][col][1]="up" 
#                        x=row-1
#                        row=x 
#             elif(("left" in seq_dir_matrix[row][col][1]) and ("up" in seq_dir_matrix[row][col][1])):
#                    if left>up: 
#                        seq_dir_matrix[row][col][1]="left"
#                        y=col-1
#                        col=y
#                    else:
#                        seq_dir_matrix[row][col][1]="up" 

#                        x=row-1
#                        row=x
#             elif(("diagonal" in seq_dir_matrix[row][col][1]) and ("up" not in seq_dir_matrix[row][col][1]) and ("left" not in seq_dir_matrix[row][col][1])):
#                x = row-1
#                row=x
#                y=col-1
#                col=y
#             elif(("left" in seq_dir_matrix[row][col][1]) and ("up" not in seq_dir_matrix[row][col][1]) and ("diagonal" not in seq_dir_matrix[row][col][1])):
#                y=col-1
#                col=y
#             else:
#                x = row-1
#                row=x

#     backtracking_path.insert(0,"")
#     index_matrix.insert(0,"")
#     b_len=len(backtracking_path) 
#     b_index=b_len-1
#     z=0

#     top_strand=["0"]*b_index
#     bottom_strand=["0"]*b_index
 
#     while backtracking_path[b_index]!="":
#         if backtracking_path[b_index][1]=="diagonal":
#             top_strand[b_index-1] = m_str[index_matrix[b_index][1]-2]
#             bottom_strand[b_index-1]=n_str[index_matrix[b_index][0]-2]
#             z=b_index-1
#             b_index=z
          
#         elif backtracking_path[b_index][1]=="left" :
#             top_strand[b_index-1] = m_str[index_matrix[b_index][1]-2]
#             bottom_strand[b_index-1]="-"
#             z=b_index-1
#             b_index=z
#         elif backtracking_path[b_index][1]=="up" :
#             top_strand[b_index-1]="-"
#             bottom_strand[b_index-1]=n_str[index_matrix[b_index][0]-2]
#             z=b_index-1
#             b_index=z
 
#     print("".join(top_strand[::-1]))
#     print("".join(bottom_strand[::-1]))
#     print(max_score)

# global_alignment(txt_sequences_filepath,txt_scoring_filepath,gap_penalty)