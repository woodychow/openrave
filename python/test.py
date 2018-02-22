from openravepy import *

import traceback
def ikfast_print_stack(): 
    tb = traceback.extract_stack() 
    pattern = '%-30s %5s %24s'  
    print( '\n'+pattern % ('        FUNCTION','LINE', 'FILE      ')) 
    keyword_of_interest = [ 'ikfast_IKFastSolver.py', 'ikfast_AST.py', 'ikfast.py', 'inversekinematics.py'] 
    print('--------------------------------------------------------------') 
    for function_call in tb: 
        for keyword in keyword_of_interest: 
            if (keyword in function_call[0]) and (function_call[2] not in 'ikfast_print_stack'): 
                print(pattern % (function_call[2], function_call[1], keyword)) 
                break 
ipython_str = 'ikfast_print_stack(); ' + \
              'from IPython.terminal import embed; ' + \
              'ipshell = embed.InteractiveShellEmbed(banner1="", config=embed.load_default_config())(local_ns=locals())' 

env = Environment()
env.Load('sanyo.mujin.dae')
robot = env.GetRobots()[0]
manip = robot.GetManipulator('flange')
iktype = IkParameterizationType.TranslationDirection5D
ikmodel = databases.inversekinematics.InverseKinematicsModel(manip = manip, \
                                                             iktype = iktype, \
                                                             freeindices = [manip.GetArmIndices()[0]])
ikmodel.autogenerate()
exec(ipython_str, globals(), locals())

ikmodel.testik('1000')
