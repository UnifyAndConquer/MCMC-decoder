# an implementation of the Markov Chain Monte Carlo Revolution paper

library(Matrix)

# permute numbers subroutine: makes random permutation of a string of numbers
randPermute <- function(codeMapping, n)
{
  for(i in 1:n)
  {
    flipPair1 = sample(alphabetIndices, 2, replace = FALSE)
    temp = codeMapping[flipPair1[1]]
    codeMapping[flipPair1[1]] = codeMapping[flipPair1[2]]
    codeMapping[flipPair1[2]] = temp 
  }
  return(codeMapping)
}

# generate matrix of co-occurrence coefficients from a training text sample (The Hitchhiker's Guide to the Galaxy)
n = 30
text = readLines("Hitchhiker Sample Text.txt")
mat = matrix(data = 0, nrow = n, ncol = n)
numLetters = 0   # number of letters in training text
for(i in 1:length(text))  # loop through lines
{
  for(j in 1:(nchar(text[i])-1))  # loop across each line
  {
    firstLetter = substr(text[i], j, j)
    secondLetter = substr(text[i], j+1, j+1)
    
    if(firstLetter == "") {firstLetter = " "}
    if(secondLetter == "") {secondLetter = " "}
    
    fLetterPos = utf8ToInt(firstLetter) - 64
    sLetterPos = utf8ToInt(secondLetter) - 64

    if(firstLetter == " ") { fLetterPos = 27 } # space
    if(firstLetter == "'") { fLetterPos = 28 } # apostrophe
    if(firstLetter == ",") { fLetterPos = 29 } # comma
    if(firstLetter == ".") { fLetterPos = 30 } # period
    
    if(secondLetter == " ") { sLetterPos = 27 } # space
    if(secondLetter == "'") { sLetterPos = 28 } # apostrophe
    if(secondLetter == ",") { sLetterPos = 29 } # comma
    if(secondLetter == ".") { sLetterPos = 30 } # period
    
    if(fLetterPos <= 30 && sLetterPos <= 30 && fLetterPos > 0 && sLetterPos > 0) 
    {
      numLetters = numLetters + 1
      mat[fLetterPos, sLetterPos] = mat[fLetterPos, sLetterPos] + 1
    }
  }
}
mat = mat + 0.00001
rowSum = apply(mat, 1, sum)
mat = mat[1:30,] / rowSum
# mat = abs(log(mat))


# generate fake coded text using a given mapping
set.seed(42)
sourceMapping = sample(alphabetIndices, n, replace = FALSE)
# 28 30  9 23 17 13 18  4 15 22 10 14 26  5  8 21 19  2  6  7 20 25 16 11  1  3 29 24 12 27 ### KEY

trialOrigText = readLines("Short Sample.txt")
trialCodedText = rep("", length(trialOrigText))

applyMapping <- function(origText, codedText, mapping)
{
  for(i in 1:length(origText))
  {
    newString = rep("", nchar(origText[i]))
    for(j in 1:length(newString))
    {
      letter = substr(origText[i], j, j)
      
      oldAscii = utf8ToInt(letter)
      oldIndex = oldAscii - 64
      
      # old letter -> old ascii -> old index -> map -> new index -> new ascii -> new letter
      
      if(letter == " ") # position 27 in the mapping sequence means "space": map[27] = "D" means "space" maps to "D"
      {
          oldIndex = 27 
      }
      else if(letter == "'") # apostrophe
      { 
          oldIndex = 28 
      }
      else if(letter == ",") # comma
      { 
          oldIndex = 29 
      }
      else if(letter == ".") # period
      {   
          oldIndex = 30
      }
      
      newIndex = mapping[oldIndex]
      newAscii = newIndex + 64
      
      if(newIndex == 27) # position 27 in the mapping sequence means "space": map[27] = "D" means "space" maps to "D"
      {
        newAscii = 32 
      }
      else if(newIndex == 28) # apostrophe
      { 
        newAscii = 39
      }
      else if(newIndex == 29) # comma
      { 
        newAscii = 44
      }
      else if(newIndex == 30) # period
      {   
        newAscii = 46
      }

      newString[j] = intToUtf8(newAscii)
      
    }
    codedText[i] = paste(newString, collapse = "")
  }
  return(codedText)
}

trialCodedText = applyMapping(trialOrigText, trialCodedText, sourceMapping)

# compute inverse mapping  !!!!!! DOESN'T WORK YET, TRIED TO DO IT IN CONSTANT SPACE. VERY EASY TO DO IN LINEAR SPACE
# inverseMap <- function(mapping)
# {
#   last = 1
#   current = mapping[last]
#   for(i in 1:length(mapping))
#   {
#     temp = mapping[current]
#     mapping[current] = last
#     last = current
#     current = temp
#   }
#   return(mapping)
# }
# 5 4 3 2 1 ---> 5 4 3 2 1  |  temp:1  current:5  last:1


# initialize random mapping
n = 30    # alphabet length, letters + space, comma, period, apostrophe
alphabetIndices = 1:n  # output alphabet indices
codeMapping = sample(alphabetIndices, n, replace = FALSE)  # input alphabet indices (mapping)

# compute plausibility of a mapping
getPlaus <- function(codedText, matrix)
{
  plaus = 1
  for(i in 1:length(codedText))
  {
    for(j in 1:(nchar(codedText[i])-1))
    {
      firstLetter = substr(codedText[i], j, j)
      secondLetter = substr(codedText[i], j+1, j+1)
      
      fLetterPos = utf8ToInt(firstLetter) - 64
      sLetterPos = utf8ToInt(secondLetter) - 64
      
      if(firstLetter == " ") { fLetterPos = 27 } # space
      if(firstLetter == "'") { fLetterPos = 28 } # apostrophe
      if(firstLetter == ",") { fLetterPos = 29 } # comma
      if(firstLetter == ".") { fLetterPos = 30 } # period
      
      if(secondLetter == " ") { sLetterPos = 27 } # space
      if(secondLetter == "'") { sLetterPos = 28 } # apostrophe
      if(secondLetter == ",") { sLetterPos = 29 } # comma
      if(secondLetter == ".") { sLetterPos = 30 } # period

      # plaus = (plaus^(j) * mat[fLetterPos, sLetterPos])^(1/(j+1))
      # plaus = plaus * mat[fLetterPos, sLetterPos]
      plaus = plaus + log(mat[fLetterPos, sLetterPos])
    }
  }
  return(plaus)
}


set.seed(42)
sourceMapping = sample(alphabetIndices, n, replace = FALSE)
# 28 30  9 23 17 13 18  4 15 22 10 14 26  5  8 21 19  2  6  7 20 25 16 11  1  3 29 24 12 27 ### KEY



# MAIN LOOP
set.seed(344)
sourceMapping = sample(alphabetIndices, n, replace = FALSE)

trialOrigText = readLines("Short Sample.txt")
trialCodedText = rep("", length(trialOrigText))
trialCodedText = applyMapping(trialOrigText, trialCodedText, sourceMapping)
trialDecodedText = rep("", length(trialCodedText))

bestText = rep("", length(trialCodedText))
rejections = 0
rejectionRate = 0
pmax = -10000

for(i in 1:100000)
{
  trialDecodedText = applyMapping(trialCodedText, trialDecodedText, sourceMapping)
  p1 = getPlaus(trialDecodedText, mat)
  newMapping = randPermute(sourceMapping, 1)
  trialDecodedText = applyMapping(trialCodedText, trialDecodedText, newMapping)
  p2 = getPlaus(trialDecodedText, mat)
  
  if(p2 > p1)
  {
    sourceMapping = newMapping  # found a better mapping
  }
  else if(rbinom(1, 1, exp(p2-p1)) == 1)  # accept with probability min(1, p2/p1)
  {
    sourceMapping = newMapping
  }
  else # reject
  {
    rejections = rejections + 1
    rejectionRate = rejections/i
  }
  
  if(p2 > pmax)
  {
    bestText = trialDecodedText
    pmax = p2
  }

  if(i%%10 == 0)
  {
  print(bestText)
  print(pmax)
  print(rejectionRate)}
}
 
# got into a situation where the optimal plausibility isn't the one where the sentence is in english.
# that's because the trial sentence is too short. decode larger texts to reduce the variance of the plausibility.
# normalize the plausibility so you can decode arbitrarily large texts.

