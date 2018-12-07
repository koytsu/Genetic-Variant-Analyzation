"""
Organizes a spreadsheet of data based on the gene and paitent to variants with gain, loss, and neutral

Author: Adam Cyze
2018
"""

import pandas as pd

data = pd.read_csv("~/Downloads/CNV_challenge_5.csv") #reads the file from the directory the data is in,into the program
FinalData = data.drop_duplicates(subset=['orderid', 'gene','call']).reset_index(drop=True) #Drops all Duplicates, and Re-aranges index
ListOfData = list() #Creates an emply list, This will be a list of lists that we can add to our final Dataframe FinalOutput

for Value in range(len(FinalData.index)): #Here we will itterate through all of the values in FinalData, and the corolating values to ListOfData 
    if FinalData.call[Value] == 'neutral':
       ListOfData.append([FinalData.orderid[Value], FinalData.gene[Value],0,1,0])
    elif FinalData.call[Value] == 'loss':
       ListOfData.append([FinalData.orderid[Value], FinalData.gene[Value],1,0,0])
    elif FinalData.call[Value] == 'gain':
       ListOfData.append([FinalData.orderid[Value], FinalData.gene[Value],0,0,1]) 
       
   #Since there are alot of neutral's its best to have that be the first if statment.   
FinalOutput = pd.DataFrame(ListOfData, columns = ['Orderid','Gene', 'Loss', 'neutral', 'gain'])#Here we add our ListOfData to the FinalOutput DataFrame
ToBreak = len(FinalOutput.index) #used to stop nullpointer error, and allows for next conditonal checking
#The Pizza time dataframe will have 2 entries for loss and neutral even if they are the same orderid and gene we need to merge them and delete the duplicate
for Value in range(len(FinalOutput.index)):
    if Value == ToBreak-1: #will stop us once we reach the end of the list so we don't try to compare the next value to something that does not exist.
        break
    elif FinalOutput.Orderid[Value] == FinalOutput.Orderid[Value+1]:#Here we compare the gene and order id if they match we add them together, then delete the duplicate
        if FinalOutput.Gene[Value] == FinalOutput.Gene[Value+1]:
            FinalOutput.loc[Value+1, 'Loss'] = FinalOutput.Loss[Value] + FinalOutput.Loss[Value+1]
            FinalOutput.loc[Value+1, 'neutral'] = FinalOutput.neutral[Value] + FinalOutput.neutral[Value+1]
            FinalOutput.loc[Value+1, 'gain'] = FinalOutput.gain[Value] + FinalOutput.gain[Value+1]
            FinalOutput = FinalOutput.drop(Value)
            #loc vs ix as it was discontinued#
            
FinalOutput = FinalOutput.reset_index(drop=True)#Here we have to rest the index as its out of order from us deleting duplicate rows. Better to do this out side the for loop.
print(FinalOutput)

#Finds the individual values for the gene on Loss, Neutral and Gain for TP53, EGFR, and CDKN2A
TP53_Loss = FinalOutput.loc[(FinalOutput['Gene'] == 'TP53') & (FinalOutput['Loss'] == 1)]
TP53_Neutral = FinalOutput.loc[(FinalOutput['Gene'] == 'TP53') & (FinalOutput['neutral'] == 1)]
TP53_Gain = FinalOutput.loc[(FinalOutput['Gene'] == 'TP53') & (FinalOutput['gain'] == 1)]
TP_Loss_count = TP53_Loss.count()
TP_Neutral_count = TP53_Neutral.count()
TP_Gain_count = TP53_Gain.count()
TP_Total = TP_Loss_count + TP_Neutral_count + TP_Gain_count
TP_Loss_percentage = TP_Loss_count/TP_Total
TP_Neutral_percentage = TP_Neutral_count/TP_Total
TP_Gain_percentage = TP_Gain_count/TP_Total
print(TP_Loss_count)
print(TP_Neutral_count)
print(TP_Gain_count)
print(TP_Loss_percentage)
print(TP_Neutral_percentage)
print(TP_Gain_percentage)

EGFR_Loss = FinalOutput.loc[(FinalOutput['Gene'] == 'EGFR') & (FinalOutput['Loss'] == 1)]
EGFR_Neutral = FinalOutput.loc[(FinalOutput['Gene'] == 'EGFR') & (FinalOutput['neutral'] == 1)]
EGFR_Gain = FinalOutput.loc[(FinalOutput['Gene'] == 'EGFR') & (FinalOutput['gain'] == 1)]
EGFR_Loss_count = EGFR_Loss.count()
EGFR_Neutral_count = EGFR_Neutral.count()
EGFR_Gain_count = EGFR_Gain.count()
EGFR_Total = EGFR_Loss_count + EGFR_Neutral_count + EGFR_Gain_count
EGFR_Loss_percentage = EGFR_Loss_count/EGFR_Total
EGFR_Neutral_percentage = EGFR_Neutral_count/EGFR_Total
EGFR_Gain_percentage = EGFR_Gain_count/EGFR_Total
print(EGFR_Loss_count)
print(EGFR_Neutral_count)
print(EGFR_Gain_count)
print(EGFR_Loss_percentage)
print(EGFR_Neutral_percentage)
print(EGFR_Gain_percentage)

CDKN2A_Loss = FinalOutput.loc[(FinalOutput['Gene'] == 'CDKN2A') & (FinalOutput['Loss'] == 1)]
CDKN2A_Neutral = FinalOutput.loc[(FinalOutput['Gene'] == 'CDKN2A') & (FinalOutput['neutral'] == 1)]
CDKN2A_Gain = FinalOutput.loc[(FinalOutput['Gene'] == 'CDKN2A') & (FinalOutput['gain'] == 1)]
CDKN2A_Loss_count = CDKN2A_Loss.count()
CDKN2A_Neutral_count = CDKN2A_Neutral.count()
CDKN2A_Gain_count = CDKN2A_Gain.count()
CDKN2A_Total = CDKN2A_Loss_count + CDKN2A_Neutral_count + CDKN2A_Gain_count
CDKN2A_Loss_percentage = CDKN2A_Loss_count/CDKN2A_Total
CDKN2A_Neutral_percentage = CDKN2A_Neutral_count/CDKN2A_Total
CDKN2A_Gain_percentage = CDKN2A_Gain_count/CDKN2A_Total
print(CDKN2A_Loss_count)
print(CDKN2A_Neutral_count)
print(CDKN2A_Gain_count)
print(CDKN2A_Loss_percentage)
print(CDKN2A_Neutral_percentage)
print(CDKN2A_Gain_percentage)

Loss_And_Neutral = FinalOutput.loc[(FinalOutput['Loss'] == 1) & (FinalOutput['neutral'] == 1) & (FinalOutput['gain'] == 0)] #checks to see what values and patients have both neutral and loss 
Gain_And_Loss = FinalOutput.loc[(FinalOutput['Loss'] == 1) & (FinalOutput['gain'] == 1) & (FinalOutput['neutral'] == 0)] #checks to see what values and patients have both gain and loss
Gain_And_Neutral = FinalOutput.loc[(FinalOutput['Loss'] == 0) & (FinalOutput['gain'] == 1) & (FinalOutput['neutral'] == 1)] #checks to see what values and patients have both neutral and gain
print(Loss_And_Neutral)
print(Gain_And_Loss)
print(Gain_And_Neutral)

LossNeutral = Loss_And_Neutral.count() #counts the total amount of entries in Loss_And Neutral
GainLoss = Gain_And_Loss.count()#counts the total amount of entries in Gain_And_Loss
GainNeutral = Gain_And_Neutral.count() #counts the total amount of entries in Gain_And_Neutral
print(LossNeutral)
print(GainLoss)
print(GainNeutral)

FinalOutput[['Loss', 'gain', 'neutral']].count() #counts the amount of times the gene is counted for each column in those 3
FinalOutput[['Loss', 'gain', 'neutral']].sum() #adds up each column separately
percentage = FinalOutput[['Loss', 'gain', 'neutral']].sum()/len(FinalOutput) #divides the sum by the length of the dataset
print(percentage)