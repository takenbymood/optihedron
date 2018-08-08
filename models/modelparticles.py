def modelpatchy36():
    #4 x 6 strips (a bit bent) + 2 x 6 star
    ind = [0]*72
    express = [6, 10, 46, 30, 20, 41, 43, 63, 61, 38, 35, 50, 53, 55, 29, 70, 12, 44, 68, 42, 40, 52, 11, 47, 49, 34, 27, 37, 25, 0, 36, 28, 17, 57, 60, 2]
    

    for express_i in express:
        ind[express_i] = 1
    print(len(set(express)))
    if len(set(express)) != len(express):
        raise ValueError
    return ind

def modelpatchy42():
    #5 x 6 strips (clean) + 2 x 6 star
    ind = [0]*72
    express = [7, 11, 5, 4, 6, 9, 10, 1, 3, 8, 12, 2, 34, 49, 45, 19, 25, 0, 69, 22, 14, 17, 60, 57, 55, 53, 29, 70, 68, 40, 33, 47, 46, 30, 20, 41, 43, 63, 61, 38, 35, 50]
    

    for express_i in express:
        ind[express_i] = 1
    print(len(set(express)))
    if len(set(express)) != len(express):
        raise ValueError
    return ind

def modelpatchy33():
    #10 x 3 blobs + 1 x 3 (1 blob decorator to spiral)
    ind = [0]*72
    express = [7, 11, 5, 4, 6, 9, 10, 1, 3, 8, 12, 2, 22, 25, 56, 14, 46, 30, 49, 45, 24, 48, 21, 36, 13, 39, 59, 55, 68, 33, 35, 66, 58]
    

    for express_i in express:
        ind[express_i] = 1
    print(len(set(express)))
    if len(set(express)) != len(express):
        raise ValueError
    return ind

def modelpatchy32():
    #10 x 3 blobs + 1 x 2 (1 blob decorator to christian cross)
    ind = [0]*72
    express = [23, 26, 7, 40, 33, 47, 51, 20, 59, 53, 12, 70, 19, 71, 24, 35, 10, 50, 16, 21, 1, 27, 37, 25, 18, 31, 43, 2, 60, 17, 14, 32]
    

    for express_i in express:
        ind[express_i] = 1
    print(len(set(express)))
    if len(set(express)) != len(express):
        raise ValueError
    return ind

def modelpatchy34():
    #8 x 4 blobs + 1 x 2 (1 free small line)
    ind = [0]*72
    express = [56, 43, 41, 6, 26, 23, 7, 42, 47, 11, 66, 71, 34, 27, 5, 37, 22, 28, 36, 69, 61, 10, 16, 63, 9, 44, 51, 55, 60, 32, 2, 39, 65, 12]
    

    for express_i in express:
        ind[express_i] = 1
    print(len(set(express)))
    if len(set(express)) != len(express):
        raise ValueError
    return ind

def modelpatchy30():
    #3 x 10
    ind = [0]*72
    express = [70, 24, 30, 71, 19, 40, 47, 33, 12, 67, 23, 60, 57, 2, 53, 16, 0, 1, 10, 35, 50, 64, 59, 8, 27, 49, 15, 3, 28, 31]
    

    for express_i in express:
        ind[express_i] = 1
    print(len(set(express)))
    if len(set(express)) != len(express):
        raise ValueError
    return ind