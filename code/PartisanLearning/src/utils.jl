import DataFrames as DF


function DictionariesToDataFrame(dictlist)
  ret = Dict()                 #Holds dataframe's columns while we build it
  #Get all unique keys from dictlist and make them entries in ret
  for x in unique([y for x in [collect(keys(x)) for x in dictlist] for y in x])
    ret[x] = []
  end
  for row in dictlist          #Loop through each row
    for (key,value) in ret     #Use ret to check all possible keys in row
      if haskey(row,key)       #Is key present in row?
        push!(value, row[key]) #Yes
      else                     #Nope
        push!(value, nothing)  #So add nothing. Keeps columns same length.
      end
    end
  end
  #Fix the data types of the columns
  for (k,v) in ret                             #Consider each column
    row_type = unique([typeof(x) for x in v])  #Get datatypes of each row
    if length(row_type)==1                     #All rows had same datatype
      row_type = row_type[1]                   #Fetch datatype
      ret[k]   = convert(Array{row_type,1}, v) #Convert column to that type
    end
  end
  #DataFrame is ready to go!
  return DF.DataFrame(ret)
end



function dictmap(l,d)
    Dict(Pair(k,l(v)) for (k,v) in d)
end

function kvdictmap(l,d)
    Dict(Pair(k,l(k,v)) for (k,v) in d)
end

function kdictmap(l,d)
    Dict(Pair(k,l(k)) for (k,_) in d)
end
