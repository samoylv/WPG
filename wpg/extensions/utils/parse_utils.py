# -*- coding: utf-8 -*-
import optparse

def var_merge(_arOpt1, _arOpt2):
    """Merge arrays of options specified in _arOpt1 and _arOpt2, eliminating duplicates
    :param _arOpt1: list providing compact description of options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    :param _arOpt2: list providing compact description of extra options
    """
    nOpt1 = len(_arOpt1)
    nOpt2 = len(_arOpt2)
    arOptRes = copy(_arOpt1)
    
    for i in range(nOpt2):
        curOpt2 = _arOpt2[i]
        curOpt2nm = curOpt2[0]
        iOpt1 = nOpt1
        for j in range(nOpt1):
            curOpt1 = arOptRes[j]
            curOpt1nm = curOpt1[0]
            if(curOpt2nm == curOpt1nm):
                iOpt1 = j
                break
        curOpt2c = copy(curOpt2)
        if((iOpt1 >= 0) and (iOpt1 < nOpt1)):
            arOptRes[iOpt1] = curOpt2c
        else:
            arOptRes.append(curOpt2c)
    return arOptRes


#****************************************************************************
def var_parse(_descr):
    """Set and parse command-prompt options from a compact description provided in _descr
    :param _descr: list providing compact description of all options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    """

    p = optparse.OptionParser()
    nOpt = len(_descr)

    listOptNamesPostParse = []
    for i in range(nOpt):
        curOpt = _descr[i]
        
        sTypeShort = curOpt[1]
        sType = 'string'
        if(sTypeShort == 'f'): sType = 'float'
        elif(sTypeShort == 'i'): sType = 'int'        
        #elif(sTypeShort == 's'): sType = 'string'

        sAct = 'store'
        if(len(curOpt) > 4): sAct = curOpt[4]

        defVal = curOpt[2]
        
        optIsList = False
        if(isinstance(defVal, list)): optIsList = True

        if(optIsList):
            sType = 'string'
            listOptNamesPostParse.append(curOpt[0])

        if(len(sTypeShort) <= 0):
            p.add_option('--' + curOpt[0], default=defVal, help=curOpt[3], action=sAct)
        else:
            p.add_option('--' + curOpt[0], type=sType, default=defVal, help=curOpt[3], action=sAct)

    v, args = p.parse_args()

    #"post-parsing" list-type options
    for i in range(len(listOptNamesPostParse)):
        curOptName = listOptNamesPostParse[i]
        valCurOpt = getattr(v, curOptName)
        
    return v
