def getRGB(dWave,maxPix=1,gamma=1):
    waveArea = [380,440,490,510,580,645,780]
    minusWave = [0,440,440,510,510,645,780]
    deltWave = [1,60,50,20,70,65,35]
    for p in range(len(waveArea)):
        if dWave<waveArea[p]:
            break

    pVar = abs(minusWave[p]-dWave)/deltWave[p]
    rgbs = [[0,0,0],[pVar,0,1],[0,pVar,1],[0,1,pVar],
            [pVar,1,0],[1,pVar,0],[1,0,0],[0,0,0]]
    
    #在光谱边缘处颜色变暗
    if (dWave>=380) & (dWave<420):
        alpha = 0.3+0.7*(dWave-380)/(420-380)
    elif (dWave>=420) & (dWave<701):
        alpha = 1.0
    elif (dWave>=701) & (dWave<780):
        alpha = 0.3+0.7*(780-dWave)/(780-700)
    else:
        alpha = 0       #非可见区

    return [maxPix*(c*alpha)**gamma for c in rgbs[p]]
getRGB(550)