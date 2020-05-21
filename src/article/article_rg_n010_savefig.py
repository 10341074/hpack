from src import *
import dill
''' 
    This script export all problems from first level of directory and save in backup folder 
    This script receives the folder as input parameter
'''
# ----------------------------------------------------------------------------------------
article_options = True
# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
    folder=sys.argv[1]
    for filename in os.listdir(folder):
        dill.load_session(folder + '/' + filename)
        print(np.all(p.z == p.isolver.save.ratio[-1, :]))

        if article_options:
            p.plot(pltype='c', colorbar=False)
        else:
            p.plot(pltype='cf')
        
        # ----- parsing folder name -----
        idx1 = folder.rfind('/')
        idx2 = folder[:idx1].rfind('/')
        shortfolder = folder[idx2+1:idx1] + '+' + folder[idx1+1:]
        print(shortfolder)
        # -----
        plt.title(shortfolder + '+' + filename)
        plt.show(block=False)
        os.system('mkdir %s_fig' %folder)
        if article_options:
            plt.savefig('%s_fig/%s+%s.eps' %(folder, filename, shortfolder), bbox_inches='tight')
        else:
            plt.savefig('%s_fig/%s+%s.svg' %(folder, filename, shortfolder), bbox_inches='tight')
    pass