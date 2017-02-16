import matplotlib.pyplot as plt
import seaborn as sns

def plot_realignment(realigned_signal,t_model):
    fig = plt.figure()
    sns.set_style('white')
    plt.plot(range(len(realigned_signal)),[t_model.model['A'.join(x[1].split('M'))][0] for x in realigned_signal],color='red',label='reference')
    plt.plot(range(len(realigned_signal)),[x[2] for x in realigned_signal],color='blue',label='aligned signal')
    inds = [index for index in range(len(realigned_signal)) if len(realigned_signal[index][1].split('M')) == 2]
    plt.plot([inds[0]]*2,[40,70],'y--')
    plt.plot([inds[-1]]*2,[40,70],'y--')
    plt.xticks(range(len(realigned_signal)),[x[1] for x in realigned_signal])
    plt.ylabel('adjusted current (pA)')
    plt.legend()
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90)
    sns.despine()
    plt.savefig('realigned_signal'+str(realigned_signal[-1][-1])+'.png',dpi=500)
    #plt.show()
