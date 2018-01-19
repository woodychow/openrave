
    

        # MAIN FUNCTION
        chaintree = solvefn(self, LinksRaw, jointvars, isolvejointvars)
        if self.useleftmultiply:
            chaintree.leftmultiply(Tleft=self.multiplyMatrix(LinksLeft), Tleftinv=self.multiplyMatrix(LinksLeftInv[::-1]))
        chaintree.dictequations += self.globalsymbols.items()
        return chaintree

    def solveFullIK_Direction3D(self, LinksRaw, jointvars, isolvejointvars, \
                                rawmanipdir = Matrix(3,1,[S.Zero,S.Zero,S.One])):
        """
        manipdir is a 3x1 column vector that prescribes the initial direction to control
        """
        self._iktype = 'direction3d'
        manipdir = Matrix(3,1,[Float(x,30) for x in rawmanipdir])
        manipdir /= sqrt(manipdir[0]*manipdir[0]+manipdir[1]*manipdir[1]+manipdir[2]*manipdir[2])
        for i in range(3):
            manipdir[i] = self.convertRealToRational(manipdir[i])

        Links = LinksRaw[:] # copy
        LinksInv = [self.affineInverse(link) for link in Links]
        
        T = self.multiplyMatrix(Links)
        self.Tfinal = zeros((4,4))
        self.Tfinal[0,0:3] = (T[0:3,0:3]*manipdir).transpose()
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)
        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        
        solvejointvars = [jointvars[i] for i in isolvejointvars]
        if len(solvejointvars) != 2:
            raise self.CannotSolveError('Need 2 joints; now there are %i' % len(solvejointvars))
        log.info('ikfast direction3d: %s', solvejointvars)

        Daccum = self.Tee[0,0:3].transpose() # 3x1
        numvarsdone = 2
        Ds = []
        Dsee = []
        for i in range(len(Links)-1):
            T = self.multiplyMatrix(Links[i:])
            D = T[0:3,0:3]*manipdir
            hasvars = [self.has(D, v) for v in solvejointvars]
            if __builtin__.sum(hasvars) == numvarsdone:
                Ds.append(D)
                Dsee.append(Daccum)
                numvarsdone -= 1
            # Tinv = self.affineInverse(Links[i])
            # Daccum = Tinv[0:3,0:3]*Daccum
            Daccum = LinksInv[i][0:3,0:3]*Daccum # still 3x1
            
        AllEquations = self.buildEquationsFromTwoSides(Ds, Dsee, jointvars, \
                                                       uselength = False)

        # check, solve, verify
        self.checkSolvability(AllEquations, solvejointvars, self.freejointvars)
        tree = self.SolveAllEquations(AllEquations, \
                                      curvars = solvejointvars, \
                                      othersolvedvars = self.freejointvars[:], \
                                      solsubs = self.freevarsubs[:], \
                                      endbranchtree = endbranchtree)
        tree = self.verifyAllEquations(AllEquations, \
                                       solvejointvars, \
                                       self.freevarsubs, \
                                       tree)
        # call AST
        chaintree = AST.SolverIKChainDirection3D([(jointvars[ijoint], ijoint) for ijoint in isolvejointvars], \
                                                 [(v,i) for v,i in izip(self.freejointvars, self.ifreejointvars)], \
                                                 Dee = self.Tee[0,0:3].transpose().subs(self.freevarsubs), \
                                                 jointtree = tree, \
                                                 Dfk = self.Tfinal[0,0:3].transpose())
        return chaintree

    def solveFullIK_Lookat3D(self, LinksRaw, jointvars, isolvejointvars, \
                             rawmanipdir = Matrix(3,1,[S.Zero,S.Zero,S.One]), \
                             rawmanippos = Matrix(3,1,[S.Zero,S.Zero,S.Zero])):
        """
        manipdir, manippos needs to be filled with a direction and position of the ray to control the lookat
        """
        self._iktype = 'lookat3d'
        
        manipdir = Matrix(3,1,[Float(x,30) for x in rawmanipdir])
        manipdir /= sqrt(manipdir[0]*manipdir[0]+manipdir[1]*manipdir[1]+manipdir[2]*manipdir[2])
        for i in range(3):
            manipdir[i] = self.convertRealToRational(manipdir[i])

        manippos = Matrix(3,1,[self.convertRealToRational(x) for x in rawmanippos])
        manippos = manippos-manipdir*manipdir.dot(manippos)
        
        Links = LinksRaw[:] # copy
        LinksInv = [self.affineInverse(link) for link in Links]
        T = self.multiplyMatrix(Links)
        self.Tfinal = zeros((4,4))
        self.Tfinal[0,0:3] = (T[0:3,0:3]*manipdir).transpose()
        self.Tfinal[0:3,3] = T[0:3,0:3]*manippos + T[0:3,3]
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)
        solvejointvars = [jointvars[i] for i in isolvejointvars]
        if len(solvejointvars) != 2:
            raise self.CannotSolveError('Need 2 joints; now there are %i' % len(solvejointvars))

        log.info('ikfast lookat3d: %s', solvejointvars)
        
        Paccum = self.Tee[0:3,3]
        numvarsdone = 2
        Positions   = []
        Positionsee = []
        
        for i in range(len(Links)-1):
            T = self.multiplyMatrix(Links[i:])
            P = T[0:3,0:3]*manippos+T[0:3,3]
            D = T[0:3,0:3]*manipdir
            hasvars = [self.has(P,v) or self.has(D,v) for v in solvejointvars]
            
            if __builtin__.sum(hasvars) == numvarsdone:
                Positions.append(P.cross(D))
                Positionsee.append(Paccum.cross(D))
                numvarsdone -= 1
                
            Tinv = self.affineInverse(Links[i])
            Paccum = Tinv[0:3,0:3]*Paccum + Tinv[0:3,3]

        frontcond = (Links[-1][0:3,0:3]*manipdir).dot(Paccum - \
                                                      (Links[-1][0:3,0:3]*manippos+Links[-1][0:3,3]))
        
        for v in jointvars:
            frontcond = frontcond.subs(self.getVariable(v).subs)

        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 checkgreaterzero = [frontcond], \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        AllEquations = self.buildEquationsFromTwoSides(Positions, Positionsee, jointvars)

        # check, solve, verify
        self.checkSolvability(AllEquations, solvejointvars, self.freejointvars)
        tree = self.SolveAllEquations(AllEquations, \
                                      curvars = solvejointvars, \
                                      othersolvedvars = self.freejointvars[:], \
                                      solsubs = self.freevarsubs[:], \
                                      endbranchtree = endbranchtree)
        tree = self.verifyAllEquations(AllEquations, \
                                       solvejointvars, \
                                       self.freevarsubs, \
                                       tree)

        # call AST
        chaintree = AST.SolverIKChainLookat3D([(jointvars[ijoint],ijoint) for ijoint in isolvejointvars], \
                                              [(v,i) for v,i in izip(self.freejointvars,self.ifreejointvars)], \
                                              Pee = self.Tee[0:3,3].subs(self.freevarsubs), \
                                              jointtree = tree, \
                                              Dfk = self.Tfinal[0,0:3].transpose(), \
                                              Pfk = self.Tfinal[0:3,3])
        chaintree.dictequations += self.ppsubs
        return chaintree

    def solveFullIK_Rotation3D(self,LinksRaw, jointvars, isolvejointvars, \
                               Rbaseraw = eye(3)):
        self._iktype = 'rotation3d'
        
        Rbase = eye(4)
        for i in range(3):
            for j in range(3):
                Rbase[i,j] = self.convertRealToRational(Rbaseraw[i,j])
                
        Tfirstright = LinksRaw[-1]*Rbase
        Links = LinksRaw[:-1]
        LinksInv = [self.affineInverse(link) for link in Links]
        self.Tfinal = self.multiplyMatrix(Links)
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)

        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        solvejointvars = [jointvars[i] for i in isolvejointvars]
        if len(solvejointvars) != 3:
            raise self.CannotSolveError('Need 3 joints; now there are %i' % len(solvejointvars))
        
        log.info('ikfast rotation3d: %s',solvejointvars)

        AllEquations = self.buildEquationsFromRotation(Links, self.Tee[0:3,0:3], solvejointvars, \
                                                       self.freejointvars)

        # check, solve, verify
        self.checkSolvability(AllEquations, solvejointvars, self.freejointvars)
        tree = self.SolveAllEquations(AllEquations, \
                                      curvars = solvejointvars[:], \
                                      othersolvedvars = self.freejointvars, \
                                      solsubs = self.freevarsubs[:], \
                                      endbranchtree = endbranchtree)
        tree = self.verifyAllEquations(AllEquations, \
                                       solvejointvars, \
                                       self.freevarsubs, \
                                       tree)

        # call AST
        chaintree =  AST.SolverIKChainRotation3D([(jointvars[ijoint],ijoint) for ijoint in isolvejointvars], \
                                                 [(v,i) for v,i in izip(self.freejointvars, self.ifreejointvars)], \
                                                 (self.Tee[0:3,0:3] * \
                                                  self.affineInverse(Tfirstright)[0:3,0:3]).subs(self.freevarsubs), \
                                                 jointtree = tree, \
                                                 Rfk = self.Tfinal[0:3,0:3] * Tfirstright[0:3,0:3])
        return chaintree

    def solveFullIK_TranslationLocalGlobal6D(self, LinksRaw, jointvars, isolvejointvars, \
                                             Tmanipraw = eye(4)):
        self._iktype = 'translation3d'
        Tgripper = eye(4)
        for i in range(4):
            for j in range(4):
                Tgripper[i,j] = self.convertRealToRational(Tmanipraw[i,j])
                
        localpos = Matrix(3,1, [self.Tee[0,0], self.Tee[1,1], self.Tee[2,2]])
        chaintree = self._solveFullIK_Translation3D(LinksRaw, \
                                                    jointvars, \
                                                    isolvejointvars, \
                                                    Tgripper[0:3,3]+Tgripper[0:3,0:3]*localpos, \
                                                    False)
        chaintree.uselocaltrans = True
        return chaintree
    
    def solveFullIK_Translation3D(self, LinksRaw, jointvars, isolvejointvars, \
                                  rawmanippos = Matrix(3,1,[S.Zero,S.Zero,S.Zero])):
        self._iktype = 'translation3d'
        manippos = Matrix(3,1,[self.convertRealToRational(x) for x in rawmanippos])
        chaintree = self._solveFullIK_Translation3D(LinksRaw, jointvars, isolvejointvars, manippos)
        return chaintree
    
    def _solveFullIK_Translation3D(self, LinksRaw, jointvars, isolvejointvars, manippos, \
                                   check = True):
        """
        Called by solveFullIK_TranslationLocalGlobal6D and solveFullIK_Translation3D
        """
        Links = LinksRaw[:]
        LinksInv = [self.affineInverse(link) for link in Links]
        self.Tfinal = self.multiplyMatrix(Links)
        
        self.Tfinal[0:3,3] = self.Tfinal[0:3,0:3]*manippos+self.Tfinal[0:3,3]
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)
        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        solvejointvars = [jointvars[i] for i in isolvejointvars]
        if len(solvejointvars) != 3:
            raise self.CannotSolveError('Need 3 joints; now there are %i' % len(solvejointvars))
        
        log.info('ikfast translation3d: %s',solvejointvars)
        Tmanipposinv = eye(4)
        Tmanipposinv[0:3,3] = -manippos
        T1links = [Tmanipposinv] + LinksInv[::-1] + [self.Tee]
        T1linksinv = [self.affineInverse(Tmanipposinv)] + Links[::-1] + [self.Teeinv]
        AllEquations = self.buildEquationsFromPositions(T1links, \
                                                        T1linksinv, \
                                                        solvejointvars, \
                                                        self.freejointvars, \
                                                        uselength = True)
        # check, solve, verify
        if check:
            self.checkSolvability(AllEquations, \
                                  solvejointvars, \
                                  self.freejointvars)
        tree = self.SolveAllEquations(AllEquations, \
                                      curvars = solvejointvars[:], \
                                      othersolvedvars = self.freejointvars, \
                                      solsubs = self.freevarsubs[:], \
                                      endbranchtree = endbranchtree)
        tree = self.verifyAllEquations(AllEquations, \
                                       solvejointvars, \
                                       self.freevarsubs, \
                                       tree)
        # call AST
        chaintree = AST.SolverIKChainTranslation3D([(jointvars[ijoint],ijoint) for ijoint in isolvejointvars], \
                                                   [(v,i) for v,i in izip(self.freejointvars, self.ifreejointvars)], \
                                                   Pee = self.Tee[0:3,3], \
                                                   jointtree = tree, \
                                                   Pfk = self.Tfinal[0:3,3])
        chaintree.dictequations += self.ppsubs
        return chaintree

    def solveFullIK_TranslationXY2D(self, LinksRaw, jointvars, isolvejointvars, \
                                    rawmanippos = Matrix(2,1,[S.Zero,S.Zero])):
        self._iktype = 'translationxy2d'
        self.ppsubs = [] # disable since pz is not valid
        self.pp = None
        manippos = Matrix(2,1,[self.convertRealToRational(x) for x in rawmanippos])
        
        Links = LinksRaw[:] # copy
        LinksInv = [self.affineInverse(link) for link in Links]
        self.Tfinal = self.multiplyMatrix(Links)
        
        self.Tfinal[0:2,3] = self.Tfinal[0:2,0:2]*manippos + self.Tfinal[0:2,3]
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)
        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        solvejointvars = [jointvars[i] for i in isolvejointvars]
        if len(solvejointvars) != 2:
            raise self.CannotSolveError('Need 2 joints; now there are %i' % len(solvejointvars))

        log.info('ikfast translationxy2d: %s',solvejointvars)
        Tmanipposinv = eye(4)
        Tmanipposinv[2,2] = S.Zero
        Tmanipposinv[0:2,3] = -manippos
        
        Tmanippos = eye(4)
        Tmanippos[2,2] = S.Zero
        Tmanippos[0:2,3] = manippos
        
        T1links = [Tmanipposinv] + LinksInv[::-1] + [self.Tee]
        T1linksinv = [Tmanippos] + Links[::-1] + [self.Teeinv]
        Taccum = eye(4)
        numvarsdone = 1
        Positions = []
        Positionsee = []
        for i in range(len(T1links)-1):
            Taccum = T1linksinv[i]*Taccum
            hasvars = [self.has(Taccum,v) for v in solvejointvars]
            if __builtin__.sum(hasvars) == numvarsdone:
                Positions.append(Taccum[0:2,3])
                Positionsee.append(self.multiplyMatrix(T1links[(i+1):])[0:2,3])
                numvarsdone += 1
            if numvarsdone > 2:
                # more than 2 variables is almost always useless
                break
            
        if len(Positions) == 0:
            Positions.append(zeros((2,1)))
            Positionsee.append(self.multiplyMatrix(T1links)[0:2,3])
            
        AllEquations = self.buildEquationsFromTwoSides(Positions, \
                                                       Positionsee, \
                                                       solvejointvars + self.freejointvars)

        self.checkSolvability(AllEquations, \
                              solvejointvars, \
                              self.freejointvars)
        
        transtree = self.SolveAllEquations(AllEquations, \
                                           curvars = solvejointvars[:], \
                                           othersolvedvars = self.freejointvars, \
                                           solsubs = self.freevarsubs[:], \
                                           endbranchtree = endbranchtree)
        transtree = self.verifyAllEquations(AllEquations, \
                                            solvejointvars, \
                                            self.freevarsubs, \
                                            transtree)
        
        chaintree = AST.SolverIKChainTranslationXY2D([(jointvars[ijoint],ijoint) for ijoint in isolvejointvars], \
                                                     [(v,i) for v,i in izip(self.freejointvars,self.ifreejointvars)], \
                                                     Pee = self.Tee[0:2,3], \
                                                     jointtree = transtree, \
                                                     Pfk = self.Tfinal[0:2,3])
        chaintree.dictequations += self.ppsubs
        return chaintree

    def solveFullIK_TranslationXYOrientation3D(self, LinksRaw, jointvars, isolvejointvars, \
                                               rawmanippos = Matrix(2,1,[S.Zero,S.Zero]), \
                                               rawangle = S.Zero):
        self._iktype = 'translationxyorientation3d'
        raise self.CannotSolveError('TranslationXYOrientation3D not implemented yet')

    def solveFullIK_Ray4D(self, LinksRaw, jointvars, isolvejointvars, \
                          rawmanipdir = Matrix(3,1,[S.Zero,S.Zero,S.One]), \
                          rawmanippos = Matrix(3,1,[S.Zero,S.Zero,S.Zero])):
        """
        manipdir, manippos needs to be filled with a direction and position of the ray to control
        """
        self._iktype = 'ray4d'
        
        manipdir = Matrix(3,1,[Float(x,30) for x in rawmanipdir])
        manipdir /= sqrt(manipdir[0]*manipdir[0]+manipdir[1]*manipdir[1]+manipdir[2]*manipdir[2])
        for i in range(3):
            manipdir[i] = self.convertRealToRational(manipdir[i])

        manippos = Matrix(3,1,[self.convertRealToRational(x) for x in rawmanippos])
        manippos = manippos-manipdir*manipdir.dot(manippos)
        
        Links = LinksRaw[:] # copy
        LinksInv = [self.affineInverse(link) for link in Links]
        T = self.multiplyMatrix(Links)
        
        self.Tfinal = zeros((4,4))
        self.Tfinal[0,0:3] = (T[0:3,0:3]*manipdir).transpose()
        self.Tfinal[0:3,3] =  T[0:3,0:3]*manippos + T[0:3,3]
        
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)
        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        solvejointvars = [jointvars[i] for i in isolvejointvars]
        if len(solvejointvars) != 4:
            raise self.CannotSolveError('Need 4 joints; now there are %i' % len(solvejointvars))

        log.info('ikfast ray4d: %s',solvejointvars)
        
        Pee = self.Tee[0:3,3]
        Dee = self.Tee[0,0:3].transpose()
        numvarsdone = 2
        Positions = []
        Positionsee = []
        for i in range(len(Links)-1):
            T = self.multiplyMatrix(Links[i:])
            P = T[0:3,0:3]*manippos + T[0:3,3]
            D = T[0:3,0:3]*manipdir
            hasvars = [self.has(P,v) or self.has(D,v) for v in solvejointvars]
            if __builtin__.sum(hasvars) == numvarsdone:
                Positions.append(P.cross(D))
                Positionsee.append(Pee.cross(Dee))
                Positions.append(D)
                Positionsee.append(Dee)
                break
            Tinv = LinksInv[i] #self.affineInverse(Links[i])
            Pee = Tinv[0:3,0:3]*Pee + Tinv[0:3,3]
            Dee = Tinv[0:3,0:3]*Dee
            
        AllEquations = self.buildEquationsFromTwoSides(Positions, \
                                                       Positionsee, \
                                                       jointvars)
        # check, solve, verify
        self.checkSolvability(AllEquations, solvejointvars, self.freejointvars)
        #try:
        tree = self.SolveAllEquations(AllEquations, \
                                      curvars = solvejointvars[:], \
                                      othersolvedvars = self.freejointvars[:], \
                                      solsubs = self.freevarsubs[:], \
                                      endbranchtree = endbranchtree)
        #except self.CannotSolveError:
            # build the raghavan/roth equations and solve with higher power methods
        #    pass
        tree = self.verifyAllEquations(AllEquations, \
                                       solvejointvars, \
                                       self.freevarsubs, \
                                       tree)

        # call AST
        chaintree = AST.SolverIKChainRay([(jointvars[ijoint], ijoint) for ijoint in isolvejointvars], \
                                         [(v,i) for v,i in izip(self.freejointvars, self.ifreejointvars)], \
                                         Pee = self.Tee[0:3,3].subs(self.freevarsubs), \
                                         Dee = self.Tee[0,0:3].transpose().subs(self.freevarsubs), \
                                         jointtree = tree, \
                                         Dfk = self.Tfinal[0,0:3].transpose(), \
                                         Pfk = self.Tfinal[0:3,3],
                                         is5dray = False)
        chaintree.dictequations += self.ppsubs
        return chaintree
    
    def solveFullIK_TranslationDirection5D(self, LinksRaw, jointvars, isolvejointvars, \
                                           rawmanipdir = Matrix(3,1,[S.Zero,S.Zero,S.One ]), \
                                           rawmanippos = Matrix(3,1,[S.Zero,S.Zero,S.Zero])):
        """
        Solves 3D translation + 3D direction
        """
        self._iktype = 'translationdirection5d'

        manippos = Matrix(3,1,[self.convertRealToRational(x) for x in rawmanippos])
        
        manipdir = Matrix(3,1,[Float(x,30) for x in rawmanipdir])
        manipdir /= sqrt(manipdir[0]*manipdir[0]+manipdir[1]*manipdir[1]+manipdir[2]*manipdir[2])
        
        # try to simplify manipdir based on possible angles
        for i in range(3):
            value = None
            manipdiri = manipdir[i]
            # TODO should restore 12 once we can capture stuff like pi/12+sqrt(12531342/5141414)
            for num in [3,4,5,6,7,8]:#,12]:
                angle = pi/num
                for v in [cos(angle), -cos(angle), sin(angle), -sin(angle)]:
                    if abs(manipdiri-v).evalf() <= (10**-self.precision):
                        value = v
                        break
                if value is not None:
                    break
            manipdir[i] = self.convertRealToRational(manipdiri, 5) if value is None else value 

        # unfortunately have to do it again...
        manipdir /= sqrt(trigsimp(manipdir[0]*manipdir[0]+manipdir[1]*manipdir[1]+manipdir[2]*manipdir[2]))
        
        offsetdist = manipdir.dot(manippos)
        manippos = manippos - manipdir*offsetdist
        Links = LinksRaw[:] # copy

        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
#        numzeros = len([manipdiri for manipdiri in manipdir if manipdiri==S.Zero])
#        # int(manipdir[0]==S.Zero) + int(manipdir[1]==S.Zero) + int(manipdir[2]==S.Zero)
        
#         if numzeros < 2:
#             try:
#                 log.info('try to rotate the last joint so that numzeros increases')
#                 assert(not self.has(Links[-1],*solvejointvars))
#                 localdir = Links[-1][0:3,0:3]*manipdir
#                 localpos = Links[-1][0:3,0:3]*manippos+Links[-1][0:3,3]
#                 AllEquations = Links[-2][0:3,0:3]*localdir
#                 tree=self.SolveAllEquations(AllEquations,curvars=solvejointvars[-1:],othersolvedvars = [],solsubs = [],endbranchtree=[])
#                 offset = tree[0].jointeval[0]
#                 endbranchtree[0].offsetvalues = [S.Zero]*len(solvejointvars)
#                 endbranchtree[0].offsetvalues[-1] = offset
#                 Toffset = Links[-2].subs(solvejointvars[-1],offset).evalf()
#                 localdir2 = Toffset[0:3,0:3]*localdir
#                 localpos2 = Toffset[0:3,0:3]*localpos+Toffset[0:3,3]
#                 Links[-1]=eye(4)
#                 for i in range(3):
#                     manipdir[i] = self.convertRealToRational(localdir2[i])
#                 manipdir /= sqrt(manipdir[0]*manipdir[0]+manipdir[1]*manipdir[1]+manipdir[2]*manipdir[2]) # unfortunately have to do it again...
#                 manippos = Matrix(3,1,[self.convertRealToRational(x) for x in localpos2])
#             except Exception, e:
#                 print 'failed to rotate joint correctly',e

        LinksInv = [self.affineInverse(link) for link in Links]
        T = self.multiplyMatrix(Links)
        
        self.Tfinal = zeros((4,4))
        self.Tfinal[0,0:3] = (T[0:3,0:3]*manipdir).transpose()
        self.Tfinal[0:3,3] =  T[0:3,0:3]*manippos + T[0:3,3]
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)

        solvejointvars = [jointvars[i] for i in isolvejointvars]
        if len(solvejointvars) != 5:
            raise self.CannotSolveError('Need 5 joints; now there are %i' % len(solvejointvars))
        
        log.info('ikfast translation direction 5d: %r, direction = %r', solvejointvars, manipdir)
        
        # if last two axes are intersecting, can divide computing position and direction
        ilinks = [i for i, Tlink in enumerate(Links) if self.has(Tlink, *solvejointvars)]
        T = self.multiplyMatrix(Links[ilinks[-2]:])
        P = T[0:3,0:3]*manippos + T[0:3,3]
        D = T[0:3,0:3]*manipdir
        tree = None
        if not self.has(P, *solvejointvars):
            Tposinv = eye(4)
            Tposinv[0:3,3] = -P
            T0links = [Tposinv] + Links[:ilinks[-2]]
            try:
                log.info('Last 2 axes are intersecting')
                tree = self.solve5DIntersectingAxes(T0links, manippos, D, solvejointvars, endbranchtree)
            except self.CannotSolveError, e:
                log.warn('%s', e)

        if tree is None:
            rawpolyeqs2_dict = {}
            coupledsolutions = None
            endbranchtree2 = []

            # Try Li&Woernle&Hiller, Kohli&Osvatic, and (commented out) Manocha&Canny solvers
            # TGN: should swap the inner INDEX for-loop with the outer SOLVEMETHOD for-loop
            #      so that the same set of equations are not to be set up twice
            
            for solvemethod in [self.solveLiWoernleHiller, self.solveKohliOsvatic]:#, self.solveManochaCanny]:

                # inner index for-loop
                for index in [2, 3]:

                    if index not in rawpolyeqs2_dict:
                        # set up equations only once, add to dictionary, and reuse it for a different solvemethod
                        T0links = LinksInv[:ilinks[index]][::-1]
                        T0 = self.multiplyMatrix(T0links)
                        T1links = Links[ilinks[index]:]
                        T1 = self.multiplyMatrix(T1links)
                        
                        p0 = T0[0:3,0:3]*self.Tee[0:3,3] + T0[0:3,3]
                        p1 = T1[0:3,0:3]*manippos + T1[0:3,3]
                        l0 = T0[0:3,0:3]*self.Tee[0,0:3].transpose()
                        l1 = T1[0:3,0:3]*manipdir

                        AllEquations = []
                        for i in range(3):
                            AllEquations.append(self.SimplifyTransform(p0[i]-p1[i]))#.expand())
                            AllEquations.append(self.SimplifyTransform(l0[i]-l1[i]))#.expand())
                        
                        # check if all joints in solvejointvars[index:] are revolute and oriented in the same way
                        # initialize
                        Tcheckorientation = None
                        checkorientationjoints = None
                        leftside = None

                        tempcheckorientationjoints = solvejointvars[:index] # up to index-1
                        if len(tempcheckorientationjoints) == 3 and \
                           all([self.IsHinge(j.name) for j in tempcheckorientationjoints]):
                            Taccums = None
                            for T in T0links:
                                if self.has(T, solvejointvars[0]):
                                    Taccums = [T]
                                elif Taccums is not None:
                                    Taccums.append(T)
                                if self.has(T, solvejointvars[index-1]):
                                    break
                            if Taccums is not None:
                                # update
                                Tcheckorientation = self.multiplyMatrix(Taccums)
                                checkorientationjoints = tempcheckorientationjoints
                                leftside = True

                        else: # TGN: if above case passes, should we still work on the other case???
                              # below there is a "if not leftside" block
                            tempcheckorientationjoints = solvejointvars[index:] # from index and on
                            if len(tempcheckorientationjoints) == 3 and \
                               all([self.IsHinge(j.name) for j in tempcheckorientationjoints]):
                                Taccums = None
                                for T in T1links:
                                    if self.has(T, solvejointvars[index]):
                                        Taccums = [T]
                                    elif Taccums is not None:
                                        Taccums.append(T)
                                    if self.has(T, solvejointvars[-1]):
                                        break
                                if Taccums is not None:
                                    # update
                                    Tcheckorientation = self.multiplyMatrix(Taccums)
                                    checkorientationjoints = tempcheckorientationjoints
                                    leftside = False
                                       
                        newsolvejointvars = solvejointvars
                        if checkorientationjoints is not None:
                            assert(Tcheckorientation is not None and \
                                   len(checkorientationjoints)==3 and \
                                   leftside is not None)
                            # TODO: consider different signs of the joints, now +, +, +
                            sumj = sum(checkorientationjoints)
                            cvar3 = cos(sumj).expand(trig = True)
                            svar3 = sin(sumj).expand(trig = True)
                            # check if T is a rotation matrix of angle sumj, hence composed of 0, 1, +/-cvar3, +/-svar3 only
                            sameorientation = True
                            for i in range(3):
                                for j in range(3):
                                    if not any([self.equal(Tcheckorientation[i,j], value) \
                                                for value in [S.Zero, cvar3, -cvar3, svar3, -svar3, S.One]]):
                                        sameorientation = False
                                        break # j for-loop
                                if sameorientation == False:
                                    break # i for-loop
                              
                            if sameorientation:
                                log.info('Joints %r have same orientation; introduce j100 and add more equations', \
                                         checkorientationjoints)
                            
                                sumjoint = Symbol('j100')
                                self.gen_trigsubs([sumjoint]) # add to trigsubs
                            
                                Tdict = {cvar3: cos(sumjoint), -cvar3: -cos(sumjoint), \
                                         svar3: sin(sumjoint), -svar3: -sin(sumjoint)  }
                                for i in range(3):
                                    for j in range(3):
                                        for value in Tdict:
                                            if self.equal(Tcheckorientation[i,j], value):
                                                Tcheckorientation[i,j] = Tdict[value]
                                                break # value for-loop
                            
                                if not leftside:
                                    # latter case where checkorientationjoints is solvejointvars[index:]
                                    newT1links = [Tcheckorientation] + Links[ilinks[-1]+1:]
                                    newT1 = self.multiplyMatrix(newT1links)
                                    newp1 = newT1[0:3,0:3]*manippos + newT1[0:3,3]
                                    newl1 = newT1[0:3,0:3]*manipdir

                                    add01 = checkorientationjoints[0] + checkorientationjoints[1]
                                    add12 = checkorientationjoints[1] + checkorientationjoints[2]
                                    add20 = checkorientationjoints[2] + checkorientationjoints[0]
                                    sumsub2 = sumjoint - checkorientationjoints[2]
                                    sumsub0 = sumjoint - checkorientationjoints[0]
                                    sumsub1 = sumjoint - checkorientationjoints[1]

                                    tosubs = (sin(checkorientationjoints[2]), \
                                              sin(sumjoint - add01).expand(trig = True))
                                    newp1 = newp1.subs(tosubs).expand()
                                    newl1 = newl1.subs(tosubs).expand()

                                    # TGN: ensure a subset of self.trigvars_subs
                                    assert(all([z in self.trigvars_subs \
                                                for z in [sumjoint, checkorientationjoints[0], checkorientationjoints[1]]]))
                                
                                    for i in range(3):
                                        newp1[i] = self.trigsimp_new(newp1[i])
                                        newl1[i] = self.trigsimp_new(newl1[i])

                                    for i in range(3):
                                        AllEquations.append(self.SimplifyTransform(p0[i]-newp1[i]).expand())
                                        AllEquations.append(self.SimplifyTransform(l0[i]-newl1[i]).expand())
                                    
                                    AllEquations.append(sumj - sumjoint)
                                    toappend = [sin(add01) - sin(sumsub2), \
                                                cos(add01) - cos(sumsub2), \
                                                sin(add12) - sin(sumsub0), \
                                                cos(add12) - cos(sumsub0), \
                                                sin(add20) - sin(sumsub1), \
                                                cos(add20) - cos(sumsub1)  ]
                                    AllEquations += [eq.expand(trig = True) for eq in toappend]
                                
                                    for consistentvalues in self.testconsistentvalues:
                                        var = self.getVariable(sumjoint)
                                        consistentvalues += var.getsubs(sumj.subs(consistentvalues))
                                    newsolvejointvars = solvejointvars + [sumjoint]
                                
                        self.sortComplexity(AllEquations)
                        log.info('Finished setting up AllEquations for index %i' % index)
                        rawpolyeqs2 = self.buildRaghavanRothEquations(p0, p1, l0, l1, solvejointvars)
                        log.info('Finished setting up Raghavan-Roth equations for index %i' % index)
                        
                        rawpolyeqs2_dict[index] = (rawpolyeqs2, AllEquations, newsolvejointvars)
                    else:
                        # reuse equations set up in the first iteration of solvemethod
                        rawpolyeqs2, AllEquations, newsolvejointvars = rawpolyeqs2_dict[index]
                    try:
                        # solvemethod is Li&Woernle&Hiller, Kohli&Osvatic, or (commented out) Manocha&Canny
                        log.info('Calling %r' % solvemethod)
                        coupledsolutions, usedvars = solvemethod(rawpolyeqs2, \
                                                                 newsolvejointvars, \
                                                                 endbranchtree = [AST.SolverSequence([endbranchtree2])], \
                                                                 AllEquationsExtra = AllEquations)
                        assert(coupledsolutions is not None)
                        log.info('solvemethod has found a solution')
                        break # inner index for-loop
                    except self.CannotSolveError, e:
                        log.warn('%s', e)
                        log.info('solvemethod has NOT found a solution')
                        continue
                    
                if coupledsolutions is not None:
                    break # outer solvemethod for-loop
                
            if coupledsolutions is None:
                raise self.CannotSolveError('Raghavan-Roth equations too complex')
            
            log.info('Solved coupled variables: %s', usedvars)
            if len(usedvars) < len(solvejointvars):
                curvars = solvejointvars[:]
                solsubs = self.freevarsubs[:]
                for var in usedvars:
                    curvars.remove(var)
                    solsubs += self.getVariable(var).subs

                # check, solve, verify
                self.checkSolvability(AllEquations, curvars, self.freejointvars + usedvars)
                localtree = self.SolveAllEquations(AllEquations, \
                                                   curvars = curvars, \
                                                   othersolvedvars = self.freejointvars + usedvars, \
                                                   solsubs = solsubs, \
                                                   endbranchtree = endbranchtree)
                localtree = self.verifyAllEquations(AllEquations, curvars, solsubs, localtree)
                
                # make it a function so compiled code is smaller
                endbranchtree2.append(AST.SolverFunction('innerfn', localtree))
                tree = coupledsolutions
            else:
                endbranchtree2 += endbranchtree
                tree = coupledsolutions

        # call AST
        chaintree = AST.SolverIKChainRay([(jointvars[ijoint],ijoint) for ijoint in isolvejointvars], \
                                         [(v,i) for v,i in izip(self.freejointvars,self.ifreejointvars)], \
                                         Pee = (self.Tee[0:3,3]-self.Tee[0,0:3].transpose()*offsetdist).subs(self.freevarsubs), \
                                         Dee = self.Tee[0,0:3].transpose().subs(self.freevarsubs), \
                                         jointtree = tree, \
                                         Dfk = self.Tfinal[0,0:3].transpose(), \
                                         Pfk = self.Tfinal[0:3,3], \
                                         is5dray = True)
        chaintree.dictequations += self.ppsubs
        return chaintree

    def solve5DIntersectingAxes(self, T0links, manippos, D, solvejointvars, endbranchtree):
        """
        Called by solveFullIK_TranslationDirection5D only.
        """
        Tmanipposinv = eye(4)
        Tmanipposinv[0:3,3] = -manippos

        LinksInv = [self.affineInverse(T) for T in T0links]
        T0 = self.multiplyMatrix(T0links)
        T1links = [Tmanipposinv] + LinksInv[::-1] + [self.Tee]
        T1linksinv = [self.affineInverse(Tmanipposinv)] + T0links[::-1] + [self.Teeinv]
        
        AllEquations = self.buildEquationsFromPositions(T1links, \
                                                        T1linksinv, \
                                                        solvejointvars, \
                                                        self.freejointvars, \
                                                        uselength = True)
        transvars = [v for v in solvejointvars if self.has(T0, v)]
        
        dirtree = []
        # AST.SolverSequence
        newendbranchtree = [AST.SolverSequence([dirtree])]

        # check, solve, verify
        self.checkSolvability(AllEquations, transvars, self.freejointvars)
        transtree = self.SolveAllEquations(AllEquations, \
                                           curvars = transvars[:], \
                                           othersolvedvars = self.freejointvars, \
                                           solsubs = self.freevarsubs[:], \
                                           endbranchtree = newendbranchtree)
        transtree = self.verifyAllEquations(AllEquations, \
                                            solvejointvars, \
                                            self.freevarsubs, \
                                            transtree)
        
        rotvars = [v for v in solvejointvars if self.has(D,v)]
        solsubs = self.freevarsubs[:]
        for v in transvars:
            solsubs += self.getVariable(v).subs
        AllEquations = self.buildEquationsFromTwoSides([D], \
                                                       [T0[0:3,0:3].transpose()*self.Tee[0,0:3].transpose()], \
                                                       solvejointvars, \
                                                       uselength = False)

        # check, solve, verify
        self.checkSolvability(AllEquations, rotvars, self.freejointvars + transvars)
        localdirtree = self.SolveAllEquations(AllEquations, \
                                              curvars = rotvars[:], \
                                              othersolvedvars = self.freejointvars + transvars, \
                                              solsubs = solsubs, \
                                              endbranchtree = endbranchtree)
        localdirtree = self.verifyAllEquations(AllEquations, rotvars, solsubs, localdirtree)
        
        # make it a function so compiled code is smaller
        dirtree.append(AST.SolverFunction('innerfn', localdirtree))
        return transtree

    def solveFullIK_6D(self, LinksRaw, jointvars, isolvejointvars, \
                       Tmanipraw = eye(4)):
        """
        Default IK solver. Solves the full 6D translation + rotation IK.

        Methods to attempt:
        (1) Check if some set of 3 intersecting axes exists. 
        (2) Try sliding non-hinge variables to left/right and check intersecting sets again.
        (3) Li-Woernle-Hiller
        (4) Kohli-Osvatic & Manocha-Canny
        """
        self._iktype = 'transform6d'

        Tgripper = Matrix(4,4, \
                          [self.convertRealToRational(Tij) for Tij in list(Tmanipraw.flat)])
        Tfirstright = LinksRaw[-1]*Tgripper
        Links = LinksRaw[:-1]
        
        #         if Links[0][0:3,0:3] == eye(3):
        #             # first axis is prismatic, so zero out self.Tee
        #             for i in range(3):
        #                 if Links[0][i,3] != S.Zero:
        #                     self.Tee[i,3] = S.Zero
        #             self.Teeinv = self.affineInverse(self.Tee)
    
        LinksInv = [self.affineInverse(link) for link in Links]
        self.Tfinal = self.multiplyMatrix(Links)
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)
        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        
        solvejointvars = [jointvars[i] for i in isolvejointvars]
        if len(solvejointvars) != 6:
            raise self.CannotSolveError('Need 6 joints; now there are %i' % len(solvejointvars))
        log.info('ikfast 6d: %s',solvejointvars)

        # check if some set of three consecutive axes intersect at one point
        # if so, the IK solution will be easier to derive
        tree = self.TestIntersectingAxes(solvejointvars, Links, LinksInv, endbranchtree)

        # Try sliding strategy
        if tree is None:
            sliderjointvars = [var for var in solvejointvars if not self.IsHinge(var.name)]
            if len(sliderjointvars) > 0:
                ZeroMatrix = zeros(4)
                for i, Tlink in enumerate(Links):
                    if self.has(Tlink, *sliderjointvars):
                        # try sliding left
                        if i > 0:
                            ileftsplit = None
                            for isplit in range(i-1, -1, -1):
                                M = self.multiplyMatrix(Links[isplit:i])
                                # test if they can swap
                                if M*Tlink - Tlink*M != ZeroMatrix:
                                    break
                                if self.has(M, *solvejointvars):
                                    # surpassed a variable!
                                    ileftsplit = isplit
                                    
                            if ileftsplit is not None:
                                # try with the new order
                                log.info('Rearranging Links[%d : %d]', ileftsplit, i+1)
                                NewLinks = list(Links)
                                NewLinks[(ileftsplit+1):(i+1)] = Links[ileftsplit:i]
                                NewLinks[ileftsplit] = Links[i]
                                NewLinksInv = list(LinksInv)
                                NewLinksInv[(ileftsplit+1):(i+1)] = Links[ileftsplit:i]
                                NewLinksInv[ileftsplit] = LinksInv[i]
                                # check again intersecting sets
                                tree = self.TestIntersectingAxes(solvejointvars, \
                                                                 NewLinks, \
                                                                 NewLinksInv, \
                                                                 endbranchtree)
                                if tree is not None:
                                    break
                                
                        # try sliding right                            
                        if i+1 < len(Links):
                            irightsplit = None
                            for isplit in range(i+1,len(Links)):
                                M = self.multiplyMatrix(Links[i+1:(isplit+1)])
                                # test if they can swap
                                if M*Tlink - Tlink*M != ZeroMatrix:
                                    break
                                if self.has(M,*solvejointvars):
                                    # surpassed a variable!
                                    irightsplit = isplit
                                    
                            if irightsplit is not None:
                                log.info('Rearranging Links[%d : %d]', i, irightsplit+1)
                                # try with the new order
                                NewLinks = list(Links)
                                NewLinks[i:irightsplit] = Links[(i+1):(irightsplit+1)]
                                NewLinks[irightsplit] = Links[i]
                                NewLinksInv = list(LinksInv)
                                NewLinksInv[i:irightsplit] = LinksInv[(i+1):(irightsplit+1)]
                                NewLinksInv[irightsplit] = LinksInv[i]
                                # check again intersecting sets
                                tree = self.TestIntersectingAxes(solvejointvars, \
                                                                 NewLinks, \
                                                                 NewLinksInv, \
                                                                 endbranchtree)
                                if tree is not None:
                                    break
                                
        if tree is None:
            # find all sets of non-intersecting axes
            linklist = list(self.iterateThreeNonIntersectingAxes(solvejointvars, Links, LinksInv))
            
            for ilinklist, (T0links, T1links) in enumerate(linklist):
                # first try LiWoernleHiller since it is the most robust
                log.info('Try Li-Woernle-Hiller on %d/%d', ilinklist, len(linklist))
                try:
                    # If T1links[-1] has no solvejointvars, then we remove it in T1links and append its inverse to T0links.
                    # As the rest of T1links contains the position unknowns, doing so simplifies computations.
                    if not self.has(T1links[-1], *solvejointvars):
                        T0links.append(self.affineInverse(T1links.pop(-1)))
                        
                    tree = self.solveFullIK_6DGeneral(T0links, \
                                                      T1links, \
                                                      solvejointvars, \
                                                      endbranchtree, \
                                                      usesolvers = 0b1) # 0b1 for LiWoernleHiller
                    break
                except (self.CannotSolveError, self.IKFeasibilityError), e:
                    log.warn('%s', e)
                    log.info('Li-Woernle-Hiller cannot solve %d/%d', ilinklist, len(linklist))
            
            if tree is None:
                for ilinklist, (T0links, T1links) in enumerate(linklist):
                    log.info('Try Kohli-Osvatic & Manocha-Canny on %d/%d', ilinklist, len(linklist))
                    try:
                        # same as above
                        if not self.has(T1links[-1], *solvejointvars):
                            T0links.append(self.affineInverse(T1links.pop(-1)))
                        # 0b10 for KohliOsvatic, 0b100 for ManochaCanny
                        tree = self.solveFullIK_6DGeneral(T0links, \
                                                          T1links, \
                                                          solvejointvars, \
                                                          endbranchtree, \
                                                          usesolvers = 0b110)
                        break
                    except (self.CannotSolveError, self.IKFeasibilityError), e:
                        log.warn('%s', e)
                        log.info('Kohli-Osvatic & Manocha-Canny cannot solve %d/%d', ilinklist, len(linklist))
                
        if tree is None:
            raise self.CannotSolveError('Cannot solve 6D mechanism!')

        # call AST
        chaintree = AST.SolverIKChainTransform6D([(jointvars[ijoint],ijoint) for ijoint in isolvejointvars], \
                                                 [(v,i) for v,i in izip(self.freejointvars, self.ifreejointvars)], \
                                                 (self.Tee*self.affineInverse(Tfirstright)).subs(self.freevarsubs), \
                                                 tree, \
                                                 Tfk = self.Tfinal*Tfirstright)
        
        chaintree.dictequations += self.ppsubs + self.npxyzsubs + self.rxpsubs
        return chaintree
    


    def solve6DIntersectingAxes(self, T0links, T1links, transvars, rotvars, solveRotationFirst, endbranchtree):
        """
        Solve 6D equations where 3 axes intersect at a point.
        These axes correspond to T0links; we use them to compute the orientation.
        The remaining 3 axes correspond to T1links; we use them to compute the position first.

        Called by TestIntersectingAxes only.
        """
        
        self._iktype = 'transform6d'
        assert(len(transvars)==3 and len(rotvars) == 3)
        T0 = self.multiplyMatrix(T0links)
        T0posoffset = eye(4)
        T0posoffset[0:3,3] = -T0[0:3,3]
        T0links = [T0posoffset] + T0links
        T1links = [T0posoffset] + T1links
        T1 = self.multiplyMatrix(T1links)

        # TGN: getting into this function means solveRotationFirst is True?
        assert(solveRotationFirst)
        
        # othersolvedvars = rotvars + self.freejointvars if solveRotationFirst else self.freejointvars[:]
        # in original code, solveRotationFirst is either None or False
        othersolvedvars = self.freejointvars[:]
        T1linksinv = [self.affineInverse(T) for T in T1links]
        AllEquations = self.buildEquationsFromPositions(T1links, T1linksinv, \
                                                        transvars, othersolvedvars, \
                                                        uselength = True)

        # TGN: This function simply passes
        self.checkSolvability(AllEquations, transvars, self.freejointvars)
        
        rottree = []
        #if solveRotationFirst:
        #    # can even get here?? it is either None or False
        #    assert(0)
        #    newendbranchtree = endbranchtree
        #else:

        # call IKFastSolver.SolverSequence constructor        
        newendbranchtree = [AST.SolverSequence([rottree])]

        # current variables (translation)
        curvars = transvars[:]
        # known values we can plug in
        solsubs = self.freevarsubs[:]

        transtree = self.SolveAllEquations(AllEquations, \
                                           curvars = curvars, \
                                           othersolvedvars = othersolvedvars[:], \
                                           solsubs = solsubs, \
                                           endbranchtree = newendbranchtree)

        transtree = self.verifyAllEquations(AllEquations, \
                                            transvars + rotvars, \
                                            # rotvars if solveRotationFirst \
                                            # else transvars+rotvars, \
                                            self.freevarsubs[:], transtree)

        solvedvarsubs = self.freevarsubs[:]
        #if solveRotationFirst:
        #    # can even get here?? it is either None or False
        #    assert(0)
        #    storesolutiontree = transtree
        #else:
        solvertree = transtree
        storesolutiontree = endbranchtree
        for tvar in transvars:
            solvedvarsubs += self.getVariable(tvar).subs

        Ree = [ Symbol('new_r%d%d'%(i,j)) for i in range(3) for j in range(3) ]
        try:
            T1sub = T1.subs(solvedvarsubs)
            othersolvedvars = self.freejointvars if solveRotationFirst else transvars+self.freejointvars
            AllEquations = self.buildEquationsFromRotation(T0links, Matrix(3,3,Ree), rotvars, othersolvedvars)
            currotvars = rotvars[:]

            # check, solve
            self.checkSolvability(AllEquations, rotvars, othersolvedvars)
            rottree += self.SolveAllEquations(AllEquations, \
                                              curvars = currotvars, \
                                              othersolvedvars = othersolvedvars, \
                                              solsubs = self.freevarsubs[:], \
                                              endbranchtree = storesolutiontree)
            
            # has to be after SolveAllEquations...?
            T1rot =  [T1[i,j] for i in range(3) for j in range(3)]
            self.globalsymbols.update(dict(izip(Ree, T1rot)))

            if len(rottree) == 0:
                raise self.CannotSolveError('Cannot solve for all rotation variables: %s:%s' % \
                                            (str(freevar), str(freevalue)))

            #if solveRotationFirst:
            #    solvertree.append(AST.SolverRotation(T1sub, rottree))
            #else:
            rottree[:] = [AST.SolverRotation(T1sub, rottree[:])]
            return solvertree
        
        finally:
            # remove the Ree global symbols out of dictionary whether a solution is found or not
            for Rij in Ree:
                self.globalsymbols.pop(Rij, None)
            
    def solveFullIK_6DGeneral(self, T0links, T1links, solvejointvars, endbranchtree, usesolvers = 0b111):
        """
        Solve 6D equations of a general kinematics structure.

        Will try Li&Woernle&Hiller (0b1), Kohli&Osvatic (0b10), and Manocha&Canny (0b100) solvers. Default is 0b111.

        These methode work only if NO sets of 3 intersecting consecutive joints exist (i.e. general 6D).

        Called by solveFullIK_6D only.
        """
        
        self._iktype = 'transform6d'
        rawpolyeqs2 = [None, None]
        coupledsolutions = None
        leftovervarstree = []
        origendbranchtree = endbranchtree
        
        solvemethods = []
        if usesolvers & 0b1:
            solvemethods.append(self.solveLiWoernleHiller)
        if usesolvers & 0b10:
            solvemethods.append(self.solveKohliOsvatic)
        if usesolvers & 0b100:
            solvemethods.append(self.solveManochaCanny)
            
        for solvemethod in solvemethods:
            if coupledsolutions is not None:
                break
            complexities = [0,0]
            for splitindex in [0, 1]:
                if rawpolyeqs2[splitindex] is None:
                    if splitindex == 0:
                        # invert, this seems to always give simpler solutions, so prioritize it
                        T0 = self.affineSimplify(self.multiplyMatrix([self.affineInverse(T) for T in T0links][::-1]))
                        T1 = self.affineSimplify(self.multiplyMatrix([self.affineInverse(T) for T in T1links][::-1]))
                    else:
                        T0 = self.affineSimplify(self.multiplyMatrix(T0links))
                        T1 = self.affineSimplify(self.multiplyMatrix(T1links))
                        
                    rawpolyeqs, numminvars = self.buildRaghavanRothEquationsFromMatrix(T0, T1, solvejointvars, \
                                                                                       simplify = False)
                    if numminvars <= 5 or len(rawpolyeqs[0][1].gens) <= 6:
                        rawpolyeqs2[splitindex] = rawpolyeqs
                        
                complexities[splitindex] = sum(self.ComputePolyComplexity(peq0) + \
                                               self.ComputePolyComplexity(peq1) \
                                               for peq0, peq1 in rawpolyeqs2[splitindex])
                
            # try the lowest complexity first and then simplify!
            sortedindices = sorted(zip(complexities, [0,1]))
            
            for complexity, splitindex in sortedindices:
                for peqs in rawpolyeqs2[splitindex]:
                    c = sum(self.codeComplexity(eq) for eq in peqs[0].coeffs())
                    if c < 5000:
                        peqs[0] = self.SimplifyTransformPoly(peqs[0])
                    else:
                        log.info('skip simplification since complexity of peq[0] = %d...', c)
                    #self.codeComplexity(poly0.as_expr()) < 2000:
                    c = sum(self.codeComplexity(eq) for eq in peqs[1].coeffs())
                    if c < 5000:
                        peqs[1] = self.SimplifyTransformPoly(peqs[1])
                    else:
                        log.info('skip simplification since complexity of peq[1] = %d...', c)
                        
                try:
                    if rawpolyeqs2[splitindex] is not None:
                        rawpolyeqs = rawpolyeqs2[splitindex]
                        endbranchtree = [AST.SolverSequence([leftovervarstree])]
                        unusedsymbols = []
                        for solvejointvar in solvejointvars:
                            usedinequs = any([var in rawpolyeqs[0][0].gens or \
                                              var in rawpolyeqs[0][1].gens \
                                              for var in self.getVariable(solvejointvar).vars])
                            if not usedinequs:
                                unusedsymbols += self.getVariable(solvejointvar).vars
                                
                        AllEquationsExtra = []
                        # prune equations for variables that are not used in rawpolyeqs
                        AllEquationsExtraPruned = [] 
                        for i in range(3):
                            for j in range(4):
                                # ensure any variable not in rawpolyeqs[0][0].gens and rawpolyeqs[0][1].gens is not used
                                eq = self.SimplifyTransform(T0[i,j]-T1[i,j])
                                if not eq.has(*unusedsymbols):
                                    AllEquationsExtraPruned.append(eq)
                                AllEquationsExtra.append(eq)
                                
                        self.sortComplexity(AllEquationsExtraPruned)
                        self.sortComplexity(AllEquationsExtra)
                        coupledsolutions, usedvars = solvemethod(rawpolyeqs, \
                                                                 solvejointvars, \
                                                                 endbranchtree = endbranchtree, \
                                                                 AllEquationsExtra = AllEquationsExtraPruned)
                        break # (complexity, splitindex) for-loop
                    
                except self.CannotSolveError, e:
                    if rawpolyeqs2[splitindex] is not None and \
                       len(rawpolyeqs2[splitindex]) > 0:
                        log.warn(u'solving %s: %s', rawpolyeqs2[splitindex][0][0].gens, e)
                    else:
                        log.warn(e)
                    continue
                
        if coupledsolutions is None:
            raise self.CannotSolveError('6D general method failed, Raghavan-Roth equations might be too complex')
        
        log.info('Solved coupled variables: %s', usedvars)
        curvars = solvejointvars[:]
        solsubs = self.freevarsubs[:]
        
        for var in usedvars:
            curvars.remove(var)
            solsubs += self.getVariable(var).subs
            
        if len(curvars) > 0:
            self.sortComplexity(AllEquationsExtra)

            # check, solve
            self.checkSolvability(AllEquationsExtra, curvars, \
                                  self.freejointvars + usedvars)
            leftovertree = self.SolveAllEquations(AllEquationsExtra, \
                                                  curvars = curvars, \
                                                  othersolvedvars = self.freejointvars + usedvars, \
                                                  solsubs = solsubs, \
                                                  endbranchtree = origendbranchtree)
            leftovervarstree.append(AST.SolverFunction('innerfn', leftovertree))
        else:
            leftovervarstree += origendbranchtree
        return coupledsolutions
    
    def solveFullIK_TranslationAxisAngle4D(self, LinksRaw, jointvars, isolvejointvars, \
                                           rawmanipdir  = Matrix(3,1,[S.One, S.Zero,S.Zero]), \
                                           rawmanippos  = Matrix(3,1,[S.Zero,S.Zero,S.Zero]), \
                                           rawglobaldir = Matrix(3,1,[S.Zero,S.Zero,S.One ]), \
                                           rawglobalnormaldir = None, \
                                           ignoreaxis = None, \
                                           rawmanipnormaldir = None, \
                                           Tmanipraw = None):
        """
        Solves 3D translation + Angle with respect to an axis

        :param rawglobalnormaldir: the axis in the base coordinate system about which we will compute a rotation.
        :param rawglobaldir: the axis normal to rawglobalnormaldir that represents the 0 angle.
        :param rawmanipnormaldir: the normal direction in the manip coordinate system for measuring the 0 angle offset. Complements rawglobalnormaldir, which should be in the base coordinate system.
        :param rawmanipdir: the axis in the manip coordinate system measuring the in-plane angle with.
        :param rawmanippos: the position in manip effector coordinate system for measuring position.
        :param Tmanipraw: extra transform of the manip coordinate system with respect to the end effector.
        """
        
        self._iktype = 'translationaxisangle4d'
        globaldir = Matrix(3,1,[Float(x,30) for x in rawglobaldir])
        globaldir /= sqrt(globaldir[0]*globaldir[0]+globaldir[1]*globaldir[1]+globaldir[2]*globaldir[2])
        for i in range(3):
            globaldir[i] = self.convertRealToRational(globaldir[i], 5)
            
        iktype = None
        if rawglobalnormaldir is None:
            globalnormaldir = None
            if globaldir[0] == S.One:
                iktype = IkType.TranslationXAxisAngle4D
            elif globaldir[1] == S.One:
                iktype = IkType.TranslationYAxisAngle4D
            elif globaldir[2] == S.One:
                iktype = IkType.TranslationZAxisAngle4D

        else:        
            globalnormaldir = Matrix(3,1,[Float(x,30) for x in rawglobalnormaldir])
            binormaldir = globalnormaldir.cross(globaldir).transpose()
            
            if globaldir[0] == S.One and globalnormaldir[2] == S.One:
                iktype = IkType.TranslationXYOrientation3D if ignoreaxis == 2 else \
                         IkType.TranslationXAxisAngleZNorm4D
                
            elif globaldir[1] == S.One and globalnormaldir[0] == S.One:
                iktype = IkType.TranslationYAxisAngleXNorm4D
                
            elif globaldir[2] == S.One and globalnormaldir[1] == S.One:
                iktype = IkType.TranslationZAxisAngleYNorm4D

        manipnormaldir = globalnormaldir if rawmanipnormaldir is None \
                         else Matrix(3,1,[self.convertRealToRational(x) for x in rawmanipnormaldir])
        
        if iktype is None:
            raise ValueError('Currently globaldir can only by x-, y-, or z-axis.')
        
        manippos = Matrix(3,1,[self.convertRealToRational(x) for x in rawmanippos])
        manipdir = Matrix(3,1,[Float(x,30) for x in rawmanipdir])
        manipdir /= sqrt(manipdir[0]*manipdir[0]+manipdir[1]*manipdir[1]+manipdir[2]*manipdir[2])
        for i in range(3):
            manipdir[i] = self.convertRealToRational(manipdir[i], 5)
        # unfortunately have to do it again...            
        manipdir /= sqrt(manipdir[0]*manipdir[0]+manipdir[1]*manipdir[1]+manipdir[2]*manipdir[2])

        Links = LinksRaw[:]
        if Tmanipraw is not None:
            Links.append(self.RoundMatrix(self.GetMatrixFromNumpy(Tmanipraw)))
        
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        
        LinksInv = [self.affineInverse(link) for link in Links]
        Tallmult = self.multiplyMatrix(Links)
        self.Tfinal = zeros(4)

        # TGN: refactored. Should we normalize v?
        v = Tallmult[0:3,0:3]*manipdir
        dotprod = globaldir.dot(v)
        self.Tfinal[0,0] = acos(dotprod) if globalnormaldir is None else \
                           atan2(binormaldir.dot(v), dotprod)

        if self.Tfinal[0,0] == nan:
            raise self.CannotSolveError('Cannot solve 4D axis angle IK. ' + \
                                        'Most likely manipulator direction is aligned with the rotation axis')
        
        self.Tfinal[0:3,3] = Tallmult[0:3,0:3]*manippos+Tallmult[0:3,3]
        self.testconsistentvalues = self.ComputeConsistentValues(jointvars, self.Tfinal, \
                                                                 numsolutions = self._numsolutions)
        
        solvejointvars = [jointvars[i] for i in isolvejointvars]
        expecteddof = 4 if ignoreaxis is None else 3
        if len(solvejointvars) != expecteddof:
            raise self.CannotSolveError('Need %d joints; now there are %i' % \
                                        (expecteddof, len(solvejointvars)))
        
        log.info('ikfast translation axis %dd, globaldir = %s, manipdir = %s : %s', \
                 expecteddof, globaldir, manipdir, solvejointvars)
        
        # if last two axes are intersecting, can divide computing position and direction
        ilinks = [i for i, Tlink in enumerate(Links) if self.has(Tlink, *solvejointvars)]
        
        Tmanipposinv = eye(4)
        Tmanipposinv[0:3,3] = -manippos
        T1links = [Tmanipposinv] + LinksInv[::-1] + [self.Tee]
        T1linksinv = [self.affineInverse(Tmanipposinv)] + Links[::-1] + [self.Teeinv]
        AllEquations = self.buildEquationsFromPositions(T1links, \
                                                        T1linksinv,\
                                                        solvejointvars, \
                                                        self.freejointvars, \
                                                        uselength = True, \
                                                        ignoreaxis = ignoreaxis)
        # TGN: only check one set of consistent values?
        if not all([abs(eq.subs(self.testconsistentvalues[0]).evalf())<=1e-10 for eq in AllEquations]):
            raise self.CannotSolveError('Some equations are not consistent with the IK, double check if IK type is correct')
        
        for ilink in ilinks: # index in range(len(ilinks)):
            # inv(T0) * T1 * manipdir = globaldir
            # => T1 * manipdir = T0 * globaldir
            T0links = LinksInv[:ilink][::-1]
            T0 = self.multiplyMatrix(T0links)
            T1links = Links[ilink:]
            T1 = self.multiplyMatrix(T1links)
            globaldir2 = T0[0:3,0:3]*globaldir
            manipdir2  = T1[0:3,0:3]*manipdir
            
            for i in range(3):
                if globaldir2[i].is_number:
                    globaldir2[i] = self.convertRealToRational(globaldir2[i])
                if manipdir2[i].is_number:
                    manipdir2[i] = self.convertRealToRational(manipdir2[i])

            # TGN: ensure solvejointvars is a subset of self.trigvars_subs
            assert(all([z in self.trigvars_subs for z in solvejointvars]))
            
            eq = self.SimplifyTransform(self.trigsimp_new(globaldir2.dot(manipdir2))) - cos(self.Tee[0])
            if self.CheckExpressionUnique(AllEquations, eq):
                AllEquations.append(eq)
                
            if globalnormaldir is not None:
                binormaldir2 = T0[0:3,0:3]*binormaldir
                for i in range(3):
                    if binormaldir2[i].is_number:
                        binormaldir2[i] = self.convertRealToRational(binormaldir2[i])
                        
                eq = self.SimplifyTransform(self.trigsimp_new(binormaldir2.dot(manipdir2))) - sin(self.Tee[0])
                if self.CheckExpressionUnique(AllEquations, eq):
                    AllEquations.append(eq)
        
        # check if planar with respect to globalnormaldir
        extravar = None
        if globalnormaldir is not None:
            if Tallmult[0:3,0:3]*manipnormaldir == globalnormaldir:
                Tnormaltest = self.rodrigues(manipnormaldir, pi/2)
                # planar, so know that the sum of all hinge joints is equal to the final angle
                # can use this fact to substitute one angle with the other values
                angles = []
                isanglepositive = []
                
                for solvejoint in solvejointvars:
                    if self.IsHinge(solvejoint.name):
                        Tall0 = Tallmult[0:3,0:3].subs(solvejoint, S.Zero)
                        Tall1 = Tallmult[0:3,0:3].subs(solvejoint, pi/2)
                        angles.append(solvejoint)
                        isanglepositive.append(all([f==S.Zero for f in Tall0*Tnormaltest-Tall1]))
                            
                Tzero = Tallmult.subs([(a,S.Zero) for a in angles])
                for i in range(3):
                    if binormaldir[i].is_number:
                        binormaldir[i] = self.convertRealToRational(binormaldir[i])
                    if manipdir[i].is_number:
                        manipdir[i] = self.convertRealToRational(manipdir[i])

                v = Tzero[0:3,0:3]*manipdir
                zeroangle = atan2(binormaldir.dot(v), globaldir.dot(v))
                eqangles = self.Tee[0]-zeroangle
                for iangle, a in enumerate(angles[:-1]):
                    eqangles += -a if isanglepositive[iangle] else a

                if not isanglepositive[-1]:
                    eqangles = -eqangles
                extravar = (angles[-1],eqangles)
                coseq = cos(eqangles).expand(trig = True)
                sineq = sin(eqangles).expand(trig = True)
                AllEquationsOld = AllEquations
                AllEquations = [self.trigsimp(eq.subs([(cos(angles[-1]), coseq), \
                                                       (sin(angles[-1]), sineq)]).expand(), \
                                              solvejointvars) for eq in AllEquationsOld]
                solvejointvarsold = list(solvejointvars)
                for var in solvejointvars:
                    if angles[-1].has(var):
                        solvejointvars.remove(var)
                        break

        self.sortComplexity(AllEquations)

        # AST.SolverStoreSolution
        endbranchtree = [AST.SolverStoreSolution(jointvars, \
                                                 isHinge = [self.IsHinge(var.name) for var in jointvars])]
        if extravar is not None:
            solution = AST.SolverSolution(extravar[0].name, \
                                          jointeval = [extravar[1]], \
                                          isHinge = self.IsHinge(extravar[0].name))
            endbranchtree.insert(0, solution)
        
        try:
            # solve, verify
            tree = self.SolveAllEquations(AllEquations, \
                                          curvars = solvejointvars[:], \
                                          othersolvedvars = self.freejointvars, \
                                          solsubs = self.freevarsubs[:], \
                                          endbranchtree = endbranchtree)
            tree = self.verifyAllEquations(AllEquations, solvejointvars, self.freevarsubs, \
                                           tree)
            
        except self.CannotSolveError, e:
            log.debug('Failed to solve using SolveAllEquations: %s', e)
            if 0:
                solvejointvar0sols = solve(AllEquations[4], solvejointvars[0])
                NewEquations = [eq.subs(solvejointvars[0], solvejointvar0sols[0]) for eq in AllEquations]
                newsolution = AST.SolverSolution(solvejointvars[0].name, \
                                                 jointeval = solvejointvar0sols, \
                                                 isHinge = self.IsHinge(solvejointvars[0].name))
                endbranchtree.insert(0, newsolution)
                tree = self.SolveAllEquations(NewEquations, \
                                              curvars = solvejointvars[1:], \
                                              othersolvedvars = self.freejointvars, \
                                              solsubs = self.freevarsubs[:], \
                                              endbranchtree = endbranchtree)
            else:
                othersolvedvars = self.freejointvars[:]
                solsubs = self.freevarsubs[:]
                freevarinvsubs = [(f[1], f[0]) for f in self.freevarsubs]
                solinvsubs = [(f[1], f[0]) for f in solsubs]
                
                # single variable solutions
                solutions = []
                gatheredexceptions = []
                for curvar in solvejointvars:
                    othervars = [var for var in solvejointvars if var != curvar]
                    curvarsym = self.getVariable(curvar)
                    raweqns = []
                    for e in AllEquations:
                        if (len(othervars) == 0 or not e.has(*othervars)) \
                           and e.has(curvar, curvarsym.htvar, curvarsym.cvar, curvarsym.svar):
                            eq = e.subs(self.freevarsubs + solsubs)
                            if self.CheckExpressionUnique(raweqns, eq):
                                raweqns.append(eq)
                    if len(raweqns) > 0:
                        try:
                            rawsolutions = self.solveSingleVariable(self.sortComplexity(raweqns), \
                                                                    curvar, \
                                                                    othersolvedvars, \
                                                                    unknownvars = solvejointvars)
                            for solution in rawsolutions:
                                self.ComputeSolutionComplexity(solution, othersolvedvars, \
                                                               solvejointvars)
                                solutions.append((solution, curvar))
                        except self.CannotSolveError, e:
                            gatheredexceptions.append((curvar.name, e))
                    else:
                        gatheredexceptions.append((curvar.name, None))
                        
                if len(solutions) == 0:
                    raise self.CannotSolveError('Failed to solve for equations. Possible errors are %s' % \
                                                gatheredexceptions)
                
                firstsolution, firstvar = solutions[0]
                othersolvedvars.append(firstvar)
                solsubs += self.getVariable(firstvar).subs
                curvars = solvejointvars[:]
                curvars.remove(firstvar)

                trigsubs = []
                polysubs = []
                polyvars = []
                for v in curvars:
                    if self.IsHinge(v.name):
                        var = self.getVariable(v)
                        polysubs += [(cos(v), var.cvar), \
                                     (sin(v), var.svar)]
                        polyvars += [var.cvar, var.svar]
                        trigsubs.append((var.svar**2, 1-var.cvar**2))
                        trigsubs.append((var.svar**3, var.svar*(1-var.cvar**2)))
                    else:
                        polyvars.append(v)
                        
                polysubsinv = [(b, a) for a, b in polysubs]
                rawpolyeqs = [Poly(Poly(eq.subs(polysubs), *polyvars).subs(trigsubs), *polyvars) \
                              for eq in AllEquations if eq.has(*curvars)]

                dummys     = []
                dummysubs  = []
                dummysubs2 = []
                dummyvars  = []
                numvars = len(polyvars)/2
                dummydenoms = [0] * numvars
                dummynums   = [0] * (numvars*2)

                for i in range(numvars):
                    dummy = Symbol('ht%s' % polyvars[2*i].name[1:]) # e.g. ht1
                    dummys.append(dummy)

                    dummydenoms[i]   = 1+dummy**2
                    dummynums[2*i]   = 1-dummy**2
                    dummynums[2*i+1] = 2*dummy
                    
                    dummysubs += [(polyvars[2*i],   dummynums[2*i]  /dummydenoms[i]), \
                                  (polyvars[2*i+1], dummynums[2*i+1]/dummydenoms[i])  ]
                    
                    var = polyvars[2*i].subs(self.invsubs).args[0]
                    dummysubs2.append((var, 2*atan(dummy)))
                    dummyvars.append((dummy, tan(0.5*var)))

                newreducedeqs = []
                for peq in rawpolyeqs:
                    maxdenom = [ max([monoms[2*i]+monoms[2*i+1] \
                                      for monoms in peq.monoms()]) \
                                 for i in range(numvar)]
                    eqnew = S.Zero
                    for monoms, c in peq.terms():
                        term = c
                        for i in range(numvar):
                            term *= dummynums[2*i]**monoms[2*i]
                            term *= dummynums[2*i+1]**monoms[2*i+1]
                            term *= dummydenoms[i]**(maxdenom[i]-monoms[2*i]-monoms[2*i+1])
                        eqnew += term
                    newreducedeqs.append(Poly(eqnew, *dummys))

                newreducedeqs.sort(lambda x,y: len(x.monoms()) - len(y.monoms()))
                ileftvar = 0
                leftvar = dummys[ileftvar]
                exportcoeffeqs = None
                for ioffset in range(len(newreducedeqs)):
                    try:
                        exportcoeffeqs, exportmonoms = self.solveDialytically(newreducedeqs[ioffset:], \
                                                                              ileftvar)
                        log.info('ioffset %d/%d' % (ioffset, len(newreducedeqs)))
                        break
                    except self.CannotSolveError, e:
                        log.debug('solveDialytically errors: %s', e)

                if exportcoeffeqs is None:
                    raise self.CannotSolveError('Failed to solveDialytically')

                coupledsolution = AST.SolverCoeffFunction(jointnames = [v.name for v in curvars], \
                                                          jointeval = [v[1] for v in dummysubs2], \
                                                          jointevalcos = [dummysubs[2*i][1] \
                                                                          for i in range(len(curvars))], \
                                                          jointevalsin = [dummysubs[2*i+1][1] \
                                                                          for i in range(len(curvars))], \
                                                          isHinges = [self.IsHinge(v.name) \
                                                                      for v in curvars], \
                                                          exportvar = [v.name for v in dummys], \
                                                          exportcoeffeqs = exportcoeffeqs, \
                                                          exportfnname = 'solvedialyticpoly12qep', \
                                                          rootmaxdim = 16)
                self.usinglapack = True
                tree = [firstsolution, coupledsolution] + endbranchtree

        # call AST to package final solution
        chaintree = AST.SolverIKChainAxisAngle([(jointvars[ijoint], ijoint) \
                                                for ijoint in isolvejointvars], \
                                               [(v,i) for v,i in izip(self.freejointvars, self.ifreejointvars)], \
                                               Pee = self.Tee[0:3,3].subs(self.freevarsubs), \
                                               angleee = self.Tee[0,0].subs(self.freevarsubs), \
                                               jointtree = tree, \
                                               Pfk = self.Tfinal[0:3,3], \
                                               anglefk = self.Tfinal[0,0], \
                                               iktype = iktype)
        chaintree.dictequations += self.ppsubs
        return chaintree


    def buildRaghavanRothEquationsFromMatrix(self, T0, T1, solvejointvars, \
                                             simplify = True, \
                                             currentcasesubs = None):
        """
        Builds the 14 equations using only 5 unknowns. 

        Method explained in [Raghavan1993]_. Basically take the position and one column/row so that the least number of variables are used.

        .. [Raghavan1993] M Raghavan and B Roth, "Inverse Kinematics of the General 6R Manipulator and related Linkages",  Journal of Mechanical Design, Volume 115, Issue 3, 1993.

        Called by solveFullIK_6DGeneral only.
        """
        
        p0 = T0[0:3,3]
        p1 = T1[0:3,3]
        p = p0-p1
        T = T0-T1
        numminvars = 100000
        for irow in range(3):
            hasvar = [self.has(T[0:3,irow], var) or \
                      self.has(p, var) for var in solvejointvars]
            numcurvars = __builtin__.sum(hasvar)
            if numminvars > numcurvars and numcurvars > 0:
                numminvars = numcurvars
                l0 = T0[0:3, irow]
                l1 = T1[0:3, irow]
                
            hasvar = [self.has(T[irow,0:3], var) or \
                      self.has(p,var) for var in solvejointvars]
            numcurvars = __builtin__.sum(hasvar)
            if numminvars > numcurvars and numcurvars > 0:
                numminvars = numcurvars
                l0 = T0[irow, 0:3].transpose()
                l1 = T1[irow, 0:3].transpose()
                
        if currentcasesubs is not None:
            p0 = p0.subs(currentcasesubs)
            p1 = p1.subs(currentcasesubs)
            l0 = l0.subs(currentcasesubs)
            l1 = l1.subs(currentcasesubs)
        return self.buildRaghavanRothEquations(p0, p1, l0, l1, \
                                               solvejointvars, simplify, currentcasesubs), \
                                               numminvars

    def buildRaghavanRothEquations(self, p0, p1, l0, l1, solvejointvars, \
                                   simplify = True, currentcasesubs = None):
        """
        Called by buildRaghavanRothEquationsFromMatrix and solveFullIK_TranslationDirection5D.
        """
        trigsubs = []
        polysubs = []
        polyvars = []
        for v in solvejointvars:
            polyvars.append(v)
            if self.IsHinge(v.name):
                var = self.getVariable(v)
                polysubs += [(cos(v),var.cvar),(sin(v),var.svar)]
                polyvars += [var.cvar,var.svar]
                trigsubs.append((var.svar**2,1-var.cvar**2))
                trigsubs.append((var.svar**3,var.svar*(1-var.cvar**2)))
        for v in self.freejointvars:
            if self.IsHinge(v.name):
                trigsubs.append((sin(v)**2,1-cos(v)**2))
                trigsubs.append((sin(v)**3,sin(v)*(1-cos(v)**2)))
        if currentcasesubs is not None:
            trigsubs += currentcasesubs
        polysubsinv = [(b,a) for a,b in polysubs]
        polyeqs = []
        for i in range(14):
            polyeqs.append([None,None])
            
        eqs = []
        for i in range(3):
            eqs.append([l0[i],l1[i]])
        for i in range(3):
            eqs.append([p0[i],p1[i]])
        l0xp0 = l0.cross(p0)
        l1xp1 = l1.cross(p1)
        for i in range(3):
            eqs.append([l0xp0[i],l1xp1[i]])
        eqs.append([p0.dot(p0),p1.dot(p1)])
        eqs.append([l0.dot(p0),l1.dot(p1)])
        starttime = time.time()
        usedvars = []
        for j in range(2):
            usedvars.append([var for var in polyvars if any([eq[j].subs(polysubs).has(var) for eq in eqs])])
        for i in range(len(eqs)):
            self._CheckPreemptFn(progress = 0.05)
            if not self.CheckEquationForVarying(eqs[i][0]) and \
               not self.CheckEquationForVarying(eqs[i][1]):
                for j in range(2):
                    if polyeqs[i][j] is not None:
                        continue
                    poly0 = Poly(eqs[i][j].subs(polysubs),*usedvars[j]).subs(trigsubs)
                    if self.codeComplexity(poly0.as_expr()) < 5000:
                        poly1 = Poly(poly0.expand().subs(trigsubs),*usedvars[j])
                        if not simplify or poly1 == S.Zero:
                            polyeqs[i][j] = poly1
                        else:
                            polyeqs[i][j] = self.SimplifyTransformPoly(poly1)
                    else:
                        polyeqs[i][j] = Poly(poly0.expand().subs(trigsubs),*usedvars[j])
        #ppl0 = p0.dot(p0)*l0 - 2*l0.dot(p0)*p0
        #ppl1 = p1.dot(p1)*l1 - 2*l1.dot(p1)*p1        
        ppl0 = polyeqs[9][0].as_expr()*l0 - 2*polyeqs[10][0].as_expr()*p0 # p0.dot(p0)*l0 - 2*l0.dot(p0)*p0
        ppl1 = polyeqs[9][1].as_expr()*l1 - 2*polyeqs[10][1].as_expr()*p1 # p1.dot(p1)*l1 - 2*l1.dot(p1)*p1
        for i in range(3):
            eqs.append([ppl0[i],ppl1[i]])
        for i in range(11, len(eqs)):
            if not self.CheckEquationForVarying(eqs[i][0]) and \
               not self.CheckEquationForVarying(eqs[i][1]):
                for j in range(2):
                    if polyeqs[i][j] is not None:
                        continue
                    poly0 = Poly(eqs[i][j].subs(polysubs),*usedvars[j]).subs(trigsubs)
                    if self.codeComplexity(poly0.as_expr()) < 5000:
                        poly1 = Poly(poly0.expand().subs(trigsubs),*usedvars[j])
                        if not simplify or poly1 == S.Zero:
                            polyeqs[i][j] = poly1
                        else:
                            polyeqs[i][j] = self.SimplifyTransformPoly(poly1)
                    else:
                        log.warn('raghavan roth equation (%d,%d) too complex', i, j)
                        polyeqs[i][j] = Poly(poly0.expand().subs(trigsubs),*usedvars[j])
        log.info('computed in %fs', time.time()-starttime)
        # prune any that have varying symbols
        # remove all fractions? having big integers could blow things up...
        return [[peq0, peq1] for peq0, peq1 in polyeqs \
                if peq0 is not None and \
                peq1 is not None and \
                not self.CheckEquationForVarying(peq0) and \
                not self.CheckEquationForVarying(peq1)]

    def reduceBothSides(self, polyeqs):
        """
        Reduces a set of equations in 5 unknowns to a set of equations with 3 unknowns by 
        solving for one side with respect to another.
        
        The input is usually the output of buildRaghavanRothEquations.

        Called by solveManochaCanny only.
        """
        
        usedvars = [polyeqs[0][0].gens, polyeqs[0][1].gens]
        reducedelayed = []
        for j in (0, 1):
            if len(usedvars[j]) <= 4:
                leftsideeqs = [polyeq[j] for polyeq in polyeqs if sum(polyeq[j].degree_list()) > 0]
                rightsideeqs = [polyeq[1-j] for polyeq in polyeqs if sum(polyeq[j].degree_list()) > 0]
                if all([all(d <= 2 for d in eq.degree_list()) for eq in leftsideeqs]):
                    try:
                        numsymbolcoeffs, _computereducedequations = self.reduceBothSidesSymbolicallyDelayed(leftsideeqs, \
                                                                                                            rightsideeqs)
                        reducedelayed.append([j, leftsideeqs, \
                                              rightsideeqs,__builtin__.sum(numsymbolcoeffs), \
                                              _computereducedequations])
                    except self.CannotSolveError:
                        continue
        
        # sort with respect to least number of symbols
        reducedelayed.sort(lambda x, y: x[3]-y[3])

        reducedeqs = []
        tree = []
        for j, leftsideeqs, rightsideeqs, numsymbolcoeffs, _computereducedequations in reducedelayed:
            self._CheckPreemptFn(progress = 0.06)
            try:
                reducedeqs2 = _computereducedequations()
                if len(reducedeqs2) == 0:
                    log.info('forcing matrix inverse (might take some time)')
                    reducedeqs2, tree = self.reduceBothSidesInverseMatrix(leftsideeqs, rightsideeqs)
                if len(reducedeqs2) > 0:
                    # success, add all the reduced equations
                    reducedeqs += [[Poly(eq[0], *usedvars[j]), \
                                    Poly(eq[1], *usedvars[1-j])] for eq in reducedeqs2] + \
                                    [[Poly(S.Zero, *polyeq[j].gens), polyeq[1-j]-polyeq[j].as_expr()] \
                                     for polyeq in polyeqs if sum(polyeq[j].degree_list()) == 0]
                    if len(reducedeqs) > 0:
                        break;
            except self.CannotSolveError, e:
                log.warn(e)
                continue

        if len(reducedeqs) > 0:
            # check if any substitutions are needed
#             for eq in reducedeqs:
#                 for j in range(2):
#                     eq[j] = Poly(eq[j].subs(trigsubs).as_expr().expand(),*eq[j].gens)
            polyeqs = reducedeqs
        return [eq for eq in polyeqs if eq[0] != S.Zero or eq[1] != S.Zero], tree

    def reduceBothSidesInverseMatrix(self, leftsideeqs, rightsideeqs):
        """
        Solve a linear system inside the program since the matrix cannot be reduced so easily.

        Called by reduceBothSides only.
        """
        allmonomsleft = set()
        for peq in leftsideeqs:
            allmonomsleft = allmonomsleft.union(set(peq.monoms()))
        allmonomsleft = list(allmonomsleft)
        allmonomsleft.sort()
        if __builtin__.sum(allmonomsleft[0]) == 0:
            allmonomsleft.pop(0)
        if len(leftsideeqs) < len(allmonomsleft):
            raise self.CannotSolveError('Left side has too few equations for the number of variables %d < %d' % \
                                        (len(leftsideeqs), len(allmonomsleft)))
        
        systemcoeffs = []
        for ileft, left in enumerate(leftsideeqs):
            coeffs = [S.Zero]*len(allmonomsleft)
            rank = 0
            for m,c in left.terms():
                if __builtin__.sum(m) > 0:
                    if c != S.Zero:
                        rank += 1
                    coeffs[allmonomsleft.index(m)] = c
            systemcoeffs.append((rank, ileft, coeffs))
        # ideally we want to try all combinations of simple equations first until we arrive to linearly independent ones.
        # However, in practice most of the first equations are linearly dependent and it takes a lot of time to prune all of them,
        # so start at the most complex
        systemcoeffs.sort(lambda x, y: -x[0]+y[0])
        # sort left and right in the same way
        leftsideeqs  = [leftsideeqs [ileft] for rank,ileft,coeffs in systemcoeffs]
        rightsideeqs = [rightsideeqs[ileft] for rank,ileft,coeffs in systemcoeffs]

        A = zeros((len(allmonomsleft), len(allmonomsleft)))
        Asymbols = []
        for i in range(A.shape[0]):
            Asymbols.append([Symbol('gconst%d_%d'%(i,j)) for j in range(A.shape[1])])
        solution = None
        for eqindices in combinations(range(len(leftsideeqs)), len(allmonomsleft)):
            self._CheckPreemptFn(progress = 0.06)
            for i,index in enumerate(eqindices):
                for k in range(len(allmonomsleft)):
                    A[i,k] = systemcoeffs[index][2][k]
            nummatrixsymbols = __builtin__.sum(1 for a in A if not a.is_number)
            if nummatrixsymbols > 10:
                # if too many symbols, evaluate numerically
                if not self.IsDeterminantNonZeroByEval(A, evalfirst=nummatrixsymbols>60): # pi_robot has 55 symbols and still finishes ok
                    continue
                log.info('found non-zero determinant by evaluation')
            else:
                det = self.det_bareis(A,*self.pvars)
                if det == S.Zero:
                    continue
                solution.checkforzeros = [self.removecommonexprs(det)]
            solution = AST.SolverMatrixInverse(A = A,Asymbols = Asymbols)
            self.usinglapack = True
            Aadj = A.adjugate() # too big to be useful for now, but can be used to see if any symbols are always 0
            break
        if solution is None:
            raise self.CannotSolveError('failed to find %d linearly independent equations' % len(allmonomsleft))
        
        reducedeqs = []
        for i in range(len(allmonomsleft)):
            var = S.One
            for k, kpower in enumerate(allmonomsleft[i]):
                if kpower != 0:
                    var *= leftsideeqs[0].gens[k]**kpower
            pright = S.Zero
            for k in range(len(allmonomsleft)):
                if Aadj[i,k] != S.Zero:
                    pright += Asymbols[i][k] * (rightsideeqs[eqindices[k]].as_expr()-leftsideeqs[eqindices[k]].TC())
            reducedeqs.append([var,pright.expand()])
            
        othereqindices = set(range(len(leftsideeqs))).difference(set(eqindices))
        for i in othereqindices:
            # have to multiply just the constant by the determinant
            neweq = rightsideeqs[i].as_expr()
            for m,c in leftsideeqs[i].terms():
                if __builtin__.sum(m) > 0:
                    neweq -= c*reducedeqs[allmonomsleft.index(m)][1]
                else:
                    neweq -= c
            reducedeqs.append([S.Zero, neweq])
        return reducedeqs, [solution]

    def reduceBothSidesSymbolically(self, *args, **kwargs):
        """
        Called by solveKohliOsvatic only.
        """
        numsymbolcoeffs, _computereducedequations = self.reduceBothSidesSymbolicallyDelayed(*args, **kwargs)
        return _computereducedequations()

    def reduceBothSidesSymbolicallyDelayed(self, leftsideeqs, rightsideeqs, \
                                           maxsymbols = 10, \
                                           usesymbols = True):
        """
        The left and right side of the equations need to have different variables.

        Called by reduceBothSides and reduceBothSidesSymbolically.
        """
        assert(len(leftsideeqs)==len(rightsideeqs))
        # first count the number of different monomials, then try to solve for each of them
        symbolgen = cse_main.numbered_symbols('const')
        vargen = cse_main.numbered_symbols('tempvar')
        rightsidedummy = []
        localsymbols = []
        dividesymbols = []
        allmonoms = dict()
        for left, right in izip(leftsideeqs, rightsideeqs):
            if right != S.Zero:
                rightsidedummy.append(symbolgen.next())
                localsymbols.append((rightsidedummy[-1], right.as_expr().expand()))
            else:
                rightsidedummy.append(S.Zero)
                
            for m in left.monoms():
                if __builtin__.sum(m) > 0 and not m in allmonoms:
                    newvar = vargen.next()
                    localsymbols.append((newvar, Poly.from_dict({m:S.One}, *left.gens).as_expr()))
                    allmonoms[m] = newvar

        if len(leftsideeqs) < len(allmonoms):
            raise self.CannotSolveError('left side has too few equations for the number of variables %d < %d' % \
                                        (len(leftsideeqs), len(allmonoms)))
        
        if len(allmonoms) == 0:
            def _returnequations():
                return [[left, right] for left, right in izip(leftsideeqs, rightsideeqs)]
            
            return 0, _returnequations
        
        unknownvars = leftsideeqs[0].gens
        newleftsideeqs = []
        numsymbolcoeffs = []
        for left,right in izip(leftsideeqs,rightsidedummy):
            left = left - right
            newleft = Poly(S.Zero,*allmonoms.values())
            leftcoeffs = [c for m,c in left.terms() if __builtin__.sum(m) > 0]
            allnumbers = all([c.is_number for c in leftcoeffs])
            if usesymbols and not allnumbers:
                # check if all the equations are within a constant from each other
                # This is neceesary since the current linear system solver cannot handle too many symbols.
                reducedeq0,common0 = self.removecommonexprs(leftcoeffs[0], returncommon = True)
                commonmults = [S.One]
                for c in leftcoeffs[1:]:
                    reducedeq1,common1 = self.removecommonexprs(c, returncommon = True)
                    if self.equal(reducedeq1,reducedeq0):
                        commonmults.append(common1/common0)
                    elif self.equal(reducedeq1,-reducedeq0):
                        commonmults.append(-common1/common0)
                    else:
                        break
                if len(commonmults) == len(leftcoeffs):
                    # divide everything by reducedeq0
                    index = 0
                    for m,c in left.terms():
                        if __builtin__.sum(m) > 0:
                            newleft = newleft + commonmults[index]*allmonoms.get(m)
                            index += 1
                        else:
                            # look in the dividesymbols for something similar
                            gmult = None
                            for gsym,geq in dividesymbols:
                                greducedeq,gcommon = self.removecommonexprs(S.One/geq, returncommon = True)
                                if self.equal(greducedeq,reducedeq0):
                                    gmult = gsym*(gcommon/common0)
                                    break
                                elif self.equal(greducedeq,-reducedeq0):
                                    gmult = gsym*(-gcommon/common0)
                                    break
                            if gmult is None:
                                gmult = symbolgen.next()
                                dividesymbols.append((gmult,S.One/leftcoeffs[0]))
                            newc = (c*gmult).subs(localsymbols).expand()
                            sym = symbolgen.next()
                            localsymbols.append((sym,newc))
                            newleft = newleft + sym
                    numsymbolcoeffs.append(0)
                    newleftsideeqs.append(newleft)
                    continue
            numsymbols = 0
            for m,c in left.terms():
                polyvar = S.One
                if __builtin__.sum(m) > 0:
                    polyvar = allmonoms.get(m)
                    if not c.is_number:
                        numsymbols += 1
                newleft = newleft + c*polyvar
            numsymbolcoeffs.append(numsymbols)
            newleftsideeqs.append(newleft)

        def _computereducedequations():
            reducedeqs = []
            # order the equations based on the number of terms
            newleftsideeqs.sort(lambda x,y: len(x.monoms()) - len(y.monoms()))
            newunknowns = newleftsideeqs[0].gens
            log.info('solving for all pairwise variables in %s, number of symbol coeffs are %s', \
                     unknownvars, __builtin__.sum(numsymbolcoeffs))
            systemcoeffs = []
            for eq in newleftsideeqs:
                eqdict = eq.as_dict()
                coeffs = []
                for i,var in enumerate(newunknowns):
                    monom = [0]*len(newunknowns)
                    monom[i] = 1
                    coeffs.append(eqdict.get(tuple(monom),S.Zero))
                monom = [0]*len(newunknowns)
                coeffs.append(-eqdict.get(tuple(monom),S.Zero))
                systemcoeffs.append(coeffs)
            
            detvars = [s for s,v in localsymbols] + self.pvars
            for eqindices in combinations(range(len(newleftsideeqs)),len(newunknowns)):
                # very quick rejection
                numsymbols = __builtin__.sum(numsymbolcoeffs[i] for i in eqindices)
                if numsymbols > maxsymbols:
                    continue
                M = Matrix([systemcoeffs[i] for i in eqindices])
                det = self.det_bareis(M[:,:-1], *detvars)
                if det == S.Zero:
                    continue
                try:
                    eqused = [newleftsideeqs[i] for i in eqindices]
                    solution=solve(eqused,newunknowns)
                except IndexError:
                    # not enough equations?
                    continue                
                if solution is not None and \
                   all([self.isValidSolution(value.subs(localsymbols)) for key,value in solution.iteritems()]):
                    # substitute 
                    solsubs = []
                    allvalid = True
                    for key,value in solution.iteritems():
                        valuesub = value.subs(localsymbols)
                        solsubs.append((key,valuesub))
                        reducedeqs.append([key.subs(localsymbols),valuesub])
                    othereqindices = set(range(len(newleftsideeqs))).difference(set(eqindices))
                    for i in othereqindices:
                        reducedeqs.append([S.Zero,(newleftsideeqs[i].subs(solsubs).subs(localsymbols)).as_expr().expand()])
                    break
            
            # remove the dividesymbols from reducedeqs
            for sym, ivalue in dividesymbols:
                value=1/ivalue
                for i in range(len(reducedeqs)):
                    eq = reducedeqs[i][1]
                    if eq.has(sym):
                        neweq = S.Zero
                        peq = Poly(eq,sym)
                        for m,c in peq.terms():
                            neweq += c*value**(peq.degree(0) - m[0])
                        reducedeqs[i][1] = neweq.expand()
                        reducedeqs[i][0] = (reducedeqs[i][0]*value**peq.degree(0)).expand()
            if len(reducedeqs) > 0:
                log.info('finished with %d equations',len(reducedeqs))
            return reducedeqs
        
        return numsymbolcoeffs, _computereducedequations
    
    def solveManochaCanny(self, \
                          rawpolyeqs, \
                          solvejointvars, \
                          endbranchtree, \
                          AllEquationsExtra = None, \
                          currentcases = None, \
                          currentcasesubs = None):
        """
        Solves the IK equations by solving a 12x12 quadratic eigenvalue problem. 

        Method explained in
        
        Dinesh Manocha and J.F. Canny. "Efficient inverse kinematics for general 6R manipulators", 
        IEEE Transactions on Robotics and Automation, Volume 10, Issue 5, Oct 1994.

        Called by solveFullIK_6DGeneral and (commented out) solveFullIK_TranslationDirection5D
        """
        log.info('Attempt Manocha-Canny general IK method')
        PolyEquations, raghavansolutiontree = self.reduceBothSides(rawpolyeqs)
        # find all equations with zeros on the left side
        RightEquations = []
        for ipeq,peq in enumerate(PolyEquations):
            if peq[0] == S.Zero:
                if len(raghavansolutiontree) > 0 or peq[1] == S.Zero:
                    # give up on optimization
                    RightEquations.append(peq[1])
                else:
                    RightEquations.append(self.SimplifyTransformPoly(peq[1]))
        
        if len(RightEquations) < 6:
            raise self.CannotSolveError('number of equations %d less than 6'%(len(RightEquations)))
        
        # sort with respect to the number of monomials
        RightEquations.sort(lambda x, y: len(x.monoms())-len(y.monoms()))
        
        # substitute with dummy=tan(half angle)
        symbols = RightEquations[0].gens
        symbolsubs = [(symbols[i].subs(self.invsubs),symbols[i]) for i in range(len(symbols))]
        unsolvedsymbols = []
        for solvejointvar in solvejointvars:
            testvars = self.getVariable(solvejointvar).vars
            if not any([v in symbols for v in testvars]):
                unsolvedsymbols += testvars

        # check that the coefficients of the reduced equations do not contain any unsolved variables
        for peq in RightEquations:
            if peq.has(*unsolvedsymbols):
                raise self.CannotSolveError('found unsolved symbol being used so ignoring: %s'%peq)
        
        log.info('solving simultaneously for symbols: %s',symbols)

        dummys = []
        dummysubs = []
        dummysubs2 = []
        dummyvars = []
        usedvars = []
        singlevariables = []
        i = 0
        while i < len(symbols):
            dummy = Symbol('ht%s'%symbols[i].name[1:])
            var = symbols[i].subs(self.invsubs)
            if not isinstance(var,Symbol):
                # [0] - cos, [1] - sin
                var = var.args[0]
                dummys.append(dummy)
                dummysubs += [(symbols[i],(1-dummy**2)/(1+dummy**2)),(symbols[i+1],2*dummy/(1+dummy**2))]
                dummysubs2.append((var,2*atan(dummy)))
                dummyvars.append((dummy,tan(0.5*var)))
                if not var in usedvars:
                    usedvars.append(var)
                i += 2
            else:
                singlevariables.append(var)
                # most likely a single variable
                dummys.append(var)
                dummysubs += [(var,var)]
                dummysubs2.append((var,var))
                if not var in usedvars:
                    usedvars.append(var)                    
                i += 1

        newreducedeqs = []
        for peq in RightEquations:
            maxdenom = dict()
            for monoms in peq.monoms():
                i = 0
                while i < len(monoms):
                    if peq.gens[i].name[0] == 'j':
                        # single variable
                        maxdenom[peq.gens[i]] = max(maxdenom.get(peq.gens[i],0),monoms[i])
                        i += 1
                    else:
                        maxdenom[peq.gens[i]] = max(maxdenom.get(peq.gens[i],0),monoms[i]+monoms[i+1])
                        i += 2
            eqnew = S.Zero
            for monoms,c in peq.terms():
                term = c
                for i in range(len(dummysubs)):
                    num,denom = fraction(dummysubs[i][1])
                    term *= num**monoms[i]
                # the denoms for 0,1 and 2,3 are the same
                i = 0
                while i < len(monoms):
                    if peq.gens[i].name[0] == 'j':
                        denom = fraction(dummysubs[i][1])[1]
                        term *= denom**(maxdenom[peq.gens[i]]-monoms[i])
                        i += 1
                    else:
                        denom = fraction(dummysubs[i][1])[1]
                        term *= denom**(maxdenom[peq.gens[i]]-monoms[i]-monoms[i+1])
                        i += 2
                eqnew += term
            newreducedeqs.append(Poly(eqnew,*dummys))
            
        # check for equations with a single variable
        if len(singlevariables) > 0:
            try:
                AllEquations = [eq.subs(self.invsubs).as_expr() for eq in newreducedeqs]
                tree = self.SolveAllEquations(AllEquations, \
                                              curvars = dummys, \
                                              othersolvedvars = [], \
                                              solsubs = self.freevarsubs, \
                                              endbranchtree = endbranchtree, \
                                              currentcases = currentcases, \
                                              currentcasesubs = currentcasesubs)
                return raghavansolutiontree + tree, usedvars
            except self.CannotSolveError:
                pass

            if 0:
                # try solving for the single variable and substituting for the rest of the equations in order to get a set of equations without the single variable
                var = singlevariables[0]
                monomindex = symbols.index(var)
                singledegreeeqs = []
                AllEquations = []
                for peq in newreducedeqs:
                    if all([m[monomindex] <= 1 for m in peq.monoms()]):
                        newpeq = Poly(peq,var)
                        if sum(newpeq.degree_list()) > 0:
                            singledegreeeqs.append(newpeq)
                        else:
                            AllEquations.append(peq.subs(self.invsubs).as_expr())
                for peq0, peq1 in combinations(singledegreeeqs,2):
                    AllEquations.append(simplify((peq0.TC()*peq1.LC() - peq0.LC()*peq1.TC()).subs(self.invsubs)))

                log.info(str(AllEquations))
                #sol=self.SolvePairVariablesHalfAngle(AllEquations,usedvars[1],usedvars[2],[])

        # choose which leftvar can determine the singularity of the following equations!
        exportcoeffeqs = None
        getsubs = raghavansolutiontree[0].getsubs if len(raghavansolutiontree) > 0 else None
        for ileftvar in range(len(dummys)):
            leftvar = dummys[ileftvar]
            try:
                exportcoeffeqs,exportmonoms = self.solveDialytically(newreducedeqs,ileftvar,getsubs=getsubs)
                break
            except self.CannotSolveError,e:
                log.warn('failed with leftvar %s: %s',leftvar,e)

        if exportcoeffeqs is None:
            raise self.CannotSolveError('failed to solve dialytically')
        if ileftvar > 0:
            raise self.CannotSolveError('solving equations dialytically succeeded with var index %d, unfortunately code generation supports only index 0'%ileftvar)
    
        jointevalcos=[d[1] for d in dummysubs if d[0].name[0] == 'c']
        jointevalsin=[d[1] for d in dummysubs if d[0].name[0] == 's']
        #jointeval=[d[1] for d in dummysubs if d[0].name[0] == 'j']
        coupledsolution = AST.SolverCoeffFunction(jointnames = [v.name for v in usedvars], \
                                                  jointeval = [v[1] for v in dummysubs2], \
                                                  jointevalcos = jointevalcos, \
                                                  jointevalsin = jointevalsin, \
                                                  isHinges = [self.IsHinge(v.name) for v in usedvars], \
                                                  exportvar = [v.name for v in dummys], \
                                                  exportcoeffeqs = exportcoeffeqs, \
                                                  exportfnname = 'solvedialyticpoly12qep', \
                                                  rootmaxdim = 16)
        self.usinglapack = True
        return raghavansolutiontree + [coupledsolution] + endbranchtree, usedvars

    def solveLiWoernleHiller(self,rawpolyeqs, \
                             solvejointvars, \
                             endbranchtree, \
                             AllEquationsExtra = [], \
                             currentcases = set(), \
                             currentcasesubs = []):
        """
        Li-Woernle-Hiller procedure covered in 
        
        Jorge Angeles, "Fundamentals of Robotics Mechanical Systems", Springer, 2007.

        Called by solveFullIK_6DGeneral and solveFullIK_TranslationDirection5D
        """
        log.info('Attempt the Li-Woernle-Hiller general IK method')
        
        if len(rawpolyeqs[0][0].gens) <len(rawpolyeqs[0][1].gens):
            for peq in rawpolyeqs:
                peq[0], peq[1] = peq[1], peq[0]
        
        originalsymbols = list(rawpolyeqs[0][0].gens)
        symbolsubs = [(originalsymbols[i].subs(self.invsubs), originalsymbols[i]) for i in range(len(originalsymbols))]
        numsymbols = 0
        for solvejointvar in solvejointvars:
            for var in self.getVariable(solvejointvar).vars:
                if var in originalsymbols:
                    numsymbols += 1
                    break # var for-loop
        if numsymbols != 3:
            raise self.CannotSolveError('Li/Woernle/Hiller method requires 3 unknown variables; now there are %d' % numsymbols)
        
        if len(originalsymbols) != 6:
            log.warn('Symbols %r are not all rotational; is this really necessary?' % originalsymbols)
            raise self.CannotSolveError('Symbols %r are not all rotational; is this really necessary?'%originalsymbols)
            
        # choose which leftvar can determine the singularity of the following equations!
        allowedindices = []
        for i in range(len(originalsymbols)):
            # if first symbol is cjX, then next should be sjX
            if originalsymbols[i].name[0] == 'c':
                assert( originalsymbols[i+1].name == 's'+originalsymbols[i].name[1:])
                if 8 == __builtin__.sum(int(peq[0].has(originalsymbols[i],originalsymbols[i+1])) for peq in rawpolyeqs):
                    allowedindices.append(i)
                    
        if len(allowedindices) == 0:
            log.warn('Could not find any variable where number of equations is exacty 8, trying all possibilities')
            for i in range(len(originalsymbols)):
                # if first symbol is cjX, then next should be sjX
                if originalsymbols[i].name[0] == 'c':
                    assert( originalsymbols[i+1].name == 's'+originalsymbols[i].name[1:])
                    allowedindices.append(i)
            #raise self.CannotSolveError('need exactly 8 equations of one variable')
            
        log.info('allowedindices = %s', allowedindices)
        for allowedindex in allowedindices:
            solutiontree = []
            checkforzeros = []
            symbols = list(originalsymbols)
            cvar = symbols[allowedindex]
            svar = symbols[allowedindex+1]
            varname = cvar.name[1:]
            tvar = Symbol('ht'+varname)
            symbols.remove(cvar)
            symbols.remove(svar)
            symbols.append(tvar)
            othersymbols = list(rawpolyeqs[0][1].gens)
            othersymbols.append(tvar)
            polyeqs = [[peq[0].as_expr(),peq[1]] for peq in rawpolyeqs if peq[0].has(cvar,svar)]
            neweqs=[]
            unusedindices = set(range(len(polyeqs)))
            for i in range(len(polyeqs)):
                if not i in unusedindices:
                    continue
                p0 = Poly(polyeqs[i][0],cvar,svar)
                p0dict=p0.as_dict()
                for j in unusedindices:
                    if j == i:
                        continue
                    p1 = Poly(polyeqs[j][0],cvar,svar) # TODO can be too complex
                    p1dict=p1.as_dict()
                    r0 = polyeqs[i][1].as_expr()
                    r1 = polyeqs[j][1].as_expr()
                    if self.equal(p0dict.get((1,0),S.Zero),-p1dict.get((0,1),S.Zero)) and \
                       self.equal(p0dict.get((0,1),S.Zero),p1dict.get((1,0),S.Zero)):
                        p0,p1 = p1,p0
                        p0dict,p1dict=p1dict,p0dict
                        r0,r1 = r1,r0
                    if self.equal(p0dict.get((1,0),S.Zero),p1dict.get((0,1),S.Zero)) and \
                       self.equal(p0dict.get((0,1),S.Zero),-p1dict.get((1,0),S.Zero)):
                        # p0+tvar*p1, p1-tvar*p0
                        # subs: tvar*svar + cvar = 1, svar-tvar*cvar=tvar
                        neweqs.append([Poly(p0dict.get((1,0),S.Zero) + \
                                            p0dict.get((0,1), S.Zero)*tvar + p0.TC() + \
                                            tvar*p1.TC(),*symbols), \
                                       Poly(r0+tvar*r1,*othersymbols)])
                        
                        neweqs.append([Poly(p0dict.get((1,0),S.Zero)*tvar - \
                                            p0dict.get((0,1),S.Zero) - \
                                            p0.TC()*tvar + \
                                            p1.TC(),*symbols), \
                                       Poly(r1-tvar*r0,*othersymbols)])
                        
                        unusedindices.remove(i)
                        unusedindices.remove(j)
                        break
            if len(neweqs) >= 8:
                break
            log.warn('allowedindex %d found %d equations where coefficients of equations match', \
                     allowedindex, len(neweqs))
            
        if len(neweqs) < 8:
            raise self.CannotSolveError('found %d equations where coefficients of equations match! need at least 8' \
                                        % len(neweqs))

        mysubs = []
        badjointvars = []
        for solvejointvar in solvejointvars:
            varsubs = self.getVariable(solvejointvar).subs
            # only choose if varsubs has entry in originalsymbols or othersymbols
            if len([s for s in varsubs if s[1] in originalsymbols+othersymbols]) > 0:
                mysubs += varsubs
            else:
                badjointvars.append(solvejointvar)

        AllEquationsExtra = [eq for eq in AllEquationsExtra if not eq.has(*badjointvars)]
        AllPolyEquationsExtra = []
        for eq in AllEquationsExtra:
            mysubs = []
            for solvejointvar in solvejointvars:
                mysubs += self.getVariable(solvejointvar).subs
            peq = Poly(eq.subs(mysubs), rawpolyeqs[0][0].gens)
            mixed = False
            for monom, coeff in peq.terms():
                if sum(monom) > 0:
                    # make sure coeff doesn't have any symbols from
                    if coeff.has(*rawpolyeqs[0][1].gens):
                        mixed = True
                        break
            if not mixed:
                AllPolyEquationsExtra.append((peq - peq.TC(), Poly(-peq.TC(), rawpolyeqs[0][1].gens)))
            
        for polyeq in [polyeqs[ipeq] for ipeq in unusedindices] + AllPolyEquationsExtra:
            p0 = Poly(polyeq[0],cvar,svar)
            p1 = polyeq[1]
            # need to substitute cvar and svar with tvar
            maxdenom = 0
            for monoms in p0.monoms():
                maxdenom=max(maxdenom,monoms[0]+monoms[1])
            eqnew = S.Zero
            for monoms,c in p0.terms():
                term = c*((1-tvar**2)**monoms[0])*(2*tvar)**monoms[1]*(1+tvar**2)**(maxdenom-monoms[0]-monoms[1])
                eqnew += term
            neweqs.append([Poly(eqnew,*symbols),Poly(p1.as_expr()*(1+tvar**2)**maxdenom,*othersymbols)])
            neweqs.append([Poly(eqnew*tvar,*symbols),Poly(p1.as_expr()*tvar*(1+tvar**2)**maxdenom,*othersymbols)])
        for ipeq,peq in enumerate(rawpolyeqs):
            if not peq[0].has(cvar,svar):
                neweqs.append([Poly(peq[0],*symbols),Poly(peq[1],*othersymbols)])
                neweqs.append([Poly(peq[0].as_expr()*tvar,*symbols),Poly(peq[1].as_expr()*tvar,*othersymbols)])
                
        # according to theory, neweqs should have 20 equations, however this isn't always the case
        # one side should have only numbers, this makes the following inverse operations trivial
        for peq in neweqs:
            peq0dict = peq[0].as_dict()
            peq[1] = peq[1] - tvar*peq0dict.get((0,0,0,0,1),S.Zero)-peq[0].TC()
            peq[0] = peq[0] - tvar*peq0dict.get((0,0,0,0,1),S.Zero)-peq[0].TC()
            
        hasreducedeqs = True
        while hasreducedeqs:
            self._CheckPreemptFn(progress = 0.08)
            hasreducedeqs = False
            for ipeq,peq in enumerate(neweqs):
                peq0dict = peq[0].as_dict()
                if len(peq0dict) == 1:
                    monomkey = peq0dict.keys()[0]
                    monomcoeff = peq0dict[monomkey]
                    monomvalue = peq[1].as_expr()
                    if sympy_smaller_073:
                        monomexpr = Monomial(*monomkey).as_expr(*peq[0].gens)
                    else:
                        monomexpr = Monomial(monomkey).as_expr(*peq[0].gens)
                    # for every equation that has this monom, substitute it
                    for ipeq2, peq2 in enumerate(neweqs):
                        if ipeq == ipeq2:
                            continue
                        for monoms,c in peq2[0].terms():
                            if monoms == monomkey:
                                # have to remove any common expressions between c and monomcoeff, or else equation can get huge
                                num2, denom2 = fraction(cancel(c/monomcoeff))
                                denom3, num3 = fraction(cancel(monomcoeff/c))
                                if denom2.is_number and denom3.is_number and abs(denom2.evalf()) > abs(denom3.evalf()):
                                    # have to select one with least abs value, or else equation will skyrocket
                                    denom2 = denom3
                                    num2 = num3
                                # have to be careful when multiplying or equation magnitude can get really skewed
                                if denom2.is_number and abs(denom2.evalf()) > 100:
                                    peq2[0] = (peq2[0] - c*monomexpr)*monomcoeff
                                    peq2[1] = peq2[1]*monomcoeff - c*monomvalue
                                else:
                                    peq2[0] = (peq2[0] - c*monomexpr)*denom2
                                    peq2[1] = peq2[1]*denom2 - num2*monomvalue
                                hasreducedeqs = True
                                break

        neweqs_full = []
        reducedeqs = []
        # filled with equations where one variable is singled out
        reducedsinglevars = [None, None, None, None]
        for ipeq, peq in enumerate(neweqs):
            peqcomb = Poly(peq[1].as_expr()-peq[0].as_expr(), peq[0].gens[:-1] + peq[1].gens)
            minimummonom = None
            for monom in (peqcomb).monoms():
                if minimummonom is None:
                    minimummonom = monom
                else:
                    minimummonom = [min(minimummonom[i], monom[i]) for i in range(len(monom))]
            
            if minimummonom is None:
                continue
            
            diveq = None
            for i in range(len(minimummonom)):
                if minimummonom[i] > 0:
                    if diveq is None:
                        diveq = peqcomb.gens[i]**minimummonom[i]
                    else:
                        diveq *= peqcomb.gens[i]**minimummonom[i]
            
            if diveq is not None:
                log.info(u'assuming equation %r is non-zero, dividing by %r', diveq, peqcomb)
                peqcombnum, r = div(peqcomb, diveq)
                assert(r==S.Zero)
                peqcombold = peqcomb # save for debugging
                peqcomb = Poly(peqcombnum, peqcomb.gens)#Poly(peqcomb / diveq, peqcomb.gens)
                peq0norm, r = div(peq[0], diveq)
                assert(r==S.Zero)
                peq1norm, r = div(peq[1], diveq)
                assert(r==S.Zero)
                peq = (Poly(peq0norm, *peq[0].gens), Poly(peq1norm, *peq[1].gens))

            coeff, factors = peqcomb.factor_list()
            
            # check if peq[1] can factor out certain monoms                
            if len(factors) > 1:
                # if either of the factors evaluate to 0, then we are ok
                # look for trivial factors that evaluate to 0 or some constant expression and put those into the checkforzeros
                eq = S.One
                divisoreq = S.One
                newzeros = []
                for factor, fdegree in factors:
                    if sum(factor.degree_list()) == 1:
                        # actually causes fractions to blow up, so don't use
                        #if factor.as_expr().has(*(originalsymbols+othersymbols)):
                        #    eq *= factor.as_expr()
                        #    continue
                        log.info(u'assuming equation %r is non-zero', factor)
                        newzeros.append(factor.as_expr())
                        divisoreq *= factor.as_expr()
                    else:
                        eq *= factor.as_expr()
                eq = coeff*eq.expand() # have to multiply by the coeff, or otherwise the equation will be weighted different and will be difficult to determine epsilons
                if peq[0] != S.Zero:
                    peq0norm, r = div(peq[0], divisoreq)
                    assert(r==S.Zero)
                    peq1norm, r = div(peq[1], divisoreq)
                    assert(r==S.Zero)
                    peq0norm = Poly(peq0norm, *peq[0].gens)
                    peq1norm = Poly(peq1norm, *peq[1].gens)
                    peq0dict = peq0norm.as_dict()
                    monom, value = peq0dict.items()[0]
                    if len(peq0dict) == 1 and __builtin__.sum(monom) == 1:
                        indices = [index for index in range(4) if monom[index] == 1]
                        if len(indices) > 0 and indices[0] < 4:
                            reducedsinglevars[indices[0]] = (value, peq1norm.as_expr())
                    isunique = True
                    for test0, test1 in neweqs_full:
                        if (self.equal(test0,peq0norm) and self.equal(test1,peq1norm)) or (self.equal(test0,-peq0norm) and self.equal(test1,-peq1norm)):
                            isunique = False
                            break
                    if isunique:
                        neweqs_full.append((peq0norm, peq1norm))
                    else:
                        log.info('not unique: %r', eq)
                else:
                    eq = eq.subs(self.freevarsubs)
                    if self.CheckExpressionUnique(reducedeqs, eq):
                        reducedeqs.append(eq)
                    else:
                        log.info('factors %d not unique: %r', len(factors), eq)
            else:
                if peq[0] != S.Zero:
                    peq0dict = peq[0].as_dict()
                    monom, value = peq0dict.items()[0]
                    if len(peq0dict) == 1 and __builtin__.sum(monom) == 1:
                        indices = [index for index in range(4) if monom[index] == 1]
                        if len(indices) > 0 and indices[0] < 4:
                            reducedsinglevars[indices[0]] = (value,peq[1].as_expr())
                    
                    isunique = True
                    for test0, test1 in neweqs_full:
                        if (self.equal(test0,peq[0]) and self.equal(test1,peq[1])) or (self.equal(test0,-peq[0]) and self.equal(test1,-peq[1])):
                            isunique = False
                            break
                    if isunique:
                        neweqs_full.append(peq)
                    else:
                        log.info('not unique: %r', peq)
                else:
                    eq = peq[1].as_expr().subs(self.freevarsubs)
                    if self.CheckExpressionUnique(reducedeqs, eq):
                        reducedeqs.append(eq)
                    else:
                        log.info('factors %d reduced not unique: %r', len(factors), eq)
        for ivar in range(2):
            if reducedsinglevars[2*ivar+0] is not None and reducedsinglevars[2*ivar+1] is not None:
                # a0*cos = b0, a1*sin = b1
                a0,b0 = reducedsinglevars[2*ivar+0]
                a1,b1 = reducedsinglevars[2*ivar+1]
                reducedeqs.append((b0*a1)**2 + (a0*b1)**2 - (a0*a1)**2)
                        
        haszeroequations = len(reducedeqs)>0
        
        neweqs_simple = []
        neweqs_complex = []
        for peq in neweqs_full:
            hassquare = False
            for monom in peq[0].monoms():
                if any([m > 1 for m in monom]):
                    hassquare = True
            if not hassquare:
                neweqs_simple.append(peq)
            else:
                neweqs_complex.append(peq)

        # add more equations by multiplying tvar. this makes it possible to have a fuller matrix
        neweqs2 = neweqs_simple + [(Poly(peq[0]*tvar, peq[0].gens), Poly(peq[1]*tvar, peq[1].gens)) for peq in neweqs_simple if not peq[0].has(tvar)]
        
        # check hacks for 5dof komatsu ik
        if 1:
            for itest in range(0,len(neweqs_simple)-1,2):
                if neweqs_simple[itest][0]*tvar-neweqs_simple[itest+1][0] == S.Zero:
                    eq = (neweqs_simple[itest][1]*tvar-neweqs_simple[itest+1][1]).as_expr()
                    if eq != S.Zero and self.CheckExpressionUnique(reducedeqs, eq):
                        reducedeqs.append(eq)
                if neweqs_simple[itest+1][0]*tvar-neweqs_simple[itest][0] == S.Zero:
                    eq = (neweqs_simple[itest+1][1]*tvar-neweqs_simple[itest][1]).as_expr()
                    if eq != S.Zero and self.CheckExpressionUnique(reducedeqs, eq):
                        reducedeqs.append(eq)
            for testrational in [Rational(-103651, 500000), Rational(-413850340369, 2000000000000), Rational(151,500), Rational(301463, 1000000)]:
                if ((neweqs_simple[0][0]*tvar - neweqs_simple[1][0])*testrational + neweqs_simple[6][0]*tvar - neweqs_simple[7][0]).expand() == S.Zero:
                    if (neweqs_simple[0][0]*tvar - neweqs_simple[1][0]) == S.Zero:
                        eq = (neweqs_simple[6][1]*tvar - neweqs_simple[7][1]).expand().as_expr()
                    else:
                        eq = ((neweqs_simple[0][1]*tvar - neweqs_simple[1][1])*testrational + neweqs_simple[6][1]*tvar - neweqs_simple[7][1]).expand().as_expr()
                    if self.CheckExpressionUnique(reducedeqs, eq):
                        reducedeqs.append(eq)
                if ((neweqs_simple[0][0]*tvar - neweqs_simple[1][0])*testrational + (neweqs_simple[2][0]*tvar-neweqs_simple[3][0])*sqrt(2)) == S.Zero:
                    if (neweqs_simple[0][0]*tvar - neweqs_simple[1][0]) == S.Zero:
                        eq = ((neweqs_simple[2][1]*tvar-neweqs_simple[3][1])*sqrt(2)).as_expr()
                    else:
                        eq = ((neweqs_simple[0][1]*tvar - neweqs_simple[1][1])*testrational + (neweqs_simple[2][1]*tvar-neweqs_simple[3][1])*sqrt(2)).as_expr()
                    if self.CheckExpressionUnique(reducedeqs, eq):
                        reducedeqs.append(eq)
                    
        neweqs_test = neweqs2#neweqs_simple
        
        allmonoms = set()
        for ipeq, peq in enumerate(neweqs_test):
            allmonoms = allmonoms.union(set(peq[0].monoms()))
        allmonoms = list(allmonoms)
        allmonoms.sort()
        
        if len(allmonoms) > len(neweqs_full) and len(reducedeqs) < 3:
            raise self.CannotSolveError('new monoms is %d>%d, reducedeqs=%d'%(len(allmonoms), len(neweqs_full), len(reducedeqs)))
        
        # the equations are ginac objects
        getsubs = None
        dictequations = []
        preprocesssolutiontree = []
        localsymbolmap = {}
        AUinv = None
        if len(allmonoms) < len(neweqs_test):
            # order with respect to complexity of [0], this is to make the inverse of A faster
            complexity = [(self.codeComplexity(peq[0].as_expr()),peq) for peq in neweqs_test]
            complexity.sort(key=itemgetter(0))
            neweqs_test = [peq for c,peq in complexity]
            A = zeros((len(neweqs_test),len(allmonoms)))
            B = zeros((len(neweqs_test),1))
            for ipeq,peq in enumerate(neweqs_test):
                for m,c in peq[0].terms():
                    A[ipeq,allmonoms.index(m)] = c.subs(self.freevarsubs)
                B[ipeq] = peq[1].as_expr().subs(self.freevarsubs)
            AU = zeros((len(allmonoms),len(allmonoms)))
            AL = zeros((A.shape[0]-len(allmonoms),len(allmonoms)))
            BU = zeros((len(allmonoms),1))
            BL = zeros((A.shape[0]-len(allmonoms),1))
            AUadjugate = None
            AU = A[:A.shape[1],:]
            nummatrixsymbols = 0
            numcomplexpows = 0 # for non-symbols, how many non-integer pows there are
            for a in AU:
                if not a.is_number:
                    nummatrixsymbols += 1
                    continue
                
                hascomplexpow = False
                for poweq in a.find(Pow):
                    if poweq.exp != S.One and poweq.exp != -S.One:
                        hascomplexpow = True
                        break
                if hascomplexpow:
                    numcomplexpows += 1
            
            # the 150 threshold is a guess
            if nummatrixsymbols > 150:
                log.info('found a non-singular matrix with %d symbols, but most likely there is a better one', nummatrixsymbols)
                raise self.CannotSolveError('matrix has too many symbols (%d), giving up since most likely will freeze'%nummatrixsymbols)
                    
            log.info('matrix has %d symbols', nummatrixsymbols)
            if nummatrixsymbols > 10:
                # if matrix symbols are great, yield so that other combinations can be tested?
                pass
            
            AUdetmat = None
            if self.IsDeterminantNonZeroByEval(AU):
                rows = range(A.shape[1])
                AUdetmat = AU
            elif not self.IsDeterminantNonZeroByEval(A.transpose()*A):
                raise self.CannotSolveError('coefficient matrix is singular')
            
            else:
                # prune the dependent vectors
                AU = A[0:1,:]
                rows = [0]
                for i in range(1,A.shape[0]):
                    self._CheckPreemptFn(progress = 0.09)
                    AU2 = AU.col_join(A[i:(i+1),:])
                    if AU2.shape[0] == AU2.shape[1]:
                        AUdetmat = AU2
                    else:
                        AUdetmat = AU2*AU2.transpose()
                    # count number of fractions/symbols
                    numausymbols = 0
                    numaufractions = 0
                    for f in AUdetmat:
                        if not f.is_number:
                            numausymbols += 1
                        if f.is_rational and not f.is_integer:
                            numaufractions += 1
                            # if fraction is really huge, give it more counts (the bigger it is, the slower it takes to compute)
                            flength = len(str(f))
                            numaufractions += int(flength/20)
                    #d = AUdetmat.det().evalf()
                    #if d == S.Zero:
                    if not self.IsDeterminantNonZeroByEval(AUdetmat, len(rows)>9 and \
                                                           (numaufractions > 120 or \
                                                            numaufractions+numausymbols > 120)):
                        log.info('skip dependent index %d, numausymbols = %d, numausymbols = %d', \
                                 i, numaufractions, numausymbols)
                        continue
                    AU = AU2
                    rows.append(i)
                    if AU.shape[0] == AU.shape[1]:
                        break
                if AU.shape[0] != AU.shape[1]:
                    raise self.CannotSolveError('could not find non-singular matrix %r'%(AU.shape,))
                
            otherrows = range(A.shape[0])
            for i,row in enumerate(rows):
                BU[i] = B[row]
                otherrows.remove(row)
            for i,row in enumerate(otherrows):
                BL[i] = B[row]
                AL[i,:] = A[row,:]
                
            if 0:#self.has(A,*self.freevars):
                AUinv = AU.inv()
                AUdet = AUdetmat.det()
                log.info('AU has symbols, so working with inverse might take some time')
                AUdet = self.trigsimp(AUdet.subs(self.freevarsubsinv),self.freejointvars).subs(self.freevarsubs)
                # find the adjugate by simplifying from the inverse
                AUadjugate = zeros(AUinv.shape)
                sinsubs = []
                for freevar in self.freejointvars:
                    var=self.getVariable(freevar)
                    for ideg in range(2,40):
                        if ideg % 2:
                            sinsubs.append((var.cvar**ideg,var.cvar*(1-var.svar**2)**int((ideg-1)/2)))
                        else:
                            sinsubs.append((var.cvar**ideg,(1-var.svar**2)**(ideg/2)))
                for i in range(AUinv.shape[0]):
                    log.info('replacing row %d', i)
                    for j in range(AUinv.shape[1]):
                        numerator,denominator = self.recursiveFraction(AUinv[i,j])
                        numerator = self.trigsimp(numerator.subs(self.freevarsubsinv),self.freejointvars).subs(self.freevarsubs)
                        numerator, common = self.removecommonexprs(numerator, onlygcd = True, returncommon = True)
                        denominator = self.trigsimp((denominator/common).subs(self.freevarsubsinv),self.freejointvars).subs(self.freevarsubs)
                        try:
                            q,r=div(numerator*AUdet,denominator,self.freevars)
                        except PolynomialError, e:
                            # 1/(-9000000*cj16 - 9000000) contains an element of the generators set
                            raise self.CannotSolveError('cannot divide for matrix inversion: %s'%e)
                        
                        if r != S.Zero:
                            # sines and cosines can mix things up a lot, so converto to half-tan
                            numerator2, numerator2d, htvarsubsinv = self.ConvertSinCosEquationToHalfTan((AUdet*numerator).subs(sinsubs).expand().subs(sinsubs).expand().subs(sinsubs).expand(), self.freejointvars)
                            denominator2, denominator2d, htvarsubsinv = self.ConvertSinCosEquationToHalfTan(denominator.subs(sinsubs).expand().subs(sinsubs).expand().subs(sinsubs).expand(), self.freejointvars)
                            extranumerator, extradenominator = fraction(numerator2d/denominator2d)
                            htvars = [v for v,eq in htvarsubsinv]
                            q,r=div((numerator2*extradenominator).expand(),(denominator2).expand(),*htvars)
                            if r != S.Zero:
                                log.warn('cannot get rid of denominator for element (%d, %d) in (%s/%s)',i, j, numerator2,denominator2)
                                #raise self.CannotSolveError('cannot get rid of denominator')
                                
                            # convert back to cos/sin in order to get rid of the denominator term?
                            sym = self.gsymbolgen.next()
                            dictequations.append((sym, q / extranumerator))
                            q = sym
                            #div(q.subs(htvarsubsinv).expand(), extranumerator.subs(htvarsubsinv).expand(), *self.freevars)
                            #newsubs=[(Symbol('htj4'), sin(self.freejointvars[0])/(1+cos(self.freejointvars[0])))]
                            #div(q.extranumerator

                        AUadjugate[i,j] = self.trigsimp(q.subs(self.freevarsubsinv),self.freejointvars).subs(self.freevarsubs)
                checkforzeros.append(self.removecommonexprs(AUdet))
                # reason we're multiplying by adjugate instead of inverse is to get rid of the potential divides by (free) parameters
                BUresult = AUadjugate*BU
                C = AL*BUresult-BL*AUdet
                for c in C:
                    reducedeqs.append(c)
            else:
                # usually if nummatrixsymbols == 0, we would just solve the inverse of the matrix. however if non-integer powers get in the way, we have to resort to solving the matrix dynamically...
                if nummatrixsymbols + numcomplexpows/4 > 40:
                    Asymbols = []
                    for i in range(AU.shape[0]):
                        Asymbols.append([Symbol('gclwh%d_%d'%(i,j)) for j in range(AU.shape[1])])
                    matrixsolution = AST.SolverMatrixInverse(A=AU,Asymbols=Asymbols)
                    getsubs = matrixsolution.getsubs
                    preprocesssolutiontree.append(matrixsolution)
                    self.usinglapack = True
                    # evaluate the inverse at various solutions and see which entries are always zero
                    isnotzero = zeros((AU.shape[0],AU.shape[1]))
                    epsilon = 1e-15
                    epsilondet = 1e-30
                    hasOneNonSingular = False
                    for itest,subs in enumerate(self.testconsistentvalues):
                        AUvalue = AU.subs(subs)
                        isallnumbers = True
                        for f in AUvalue:
                            if not f.is_number:
                                isallnumbers = False
                                break
                        if isallnumbers:
                            # compute more precise determinant
                            AUdetvalue = AUvalue.det()
                        else:
                            AUdetvalue = AUvalue.evalf().det().evalf()
                        if abs(AUdetvalue) > epsilondet:# != S.Zero:
                            hasOneNonSingular = True
                            AUinvvalue = AUvalue.evalf().inv()
                            for i in range(AUinvvalue.shape[0]):
                                for j in range(AUinvvalue.shape[1]):
                                    # since making numerical approximations, need a good value for zero
                                    if abs(AUinvvalue[i,j]) > epsilon:#!= S.Zero:
                                        isnotzero[i,j] = 1
                    if not hasOneNonSingular:
                        raise self.CannotSolveError('inverse matrix is always singular')
                    
                    AUinv = zeros((AU.shape[0],AU.shape[1]))
                    for i in range(AUinv.shape[0]):
                        for j in range(AUinv.shape[1]):
                            if isnotzero[i,j] == 0:
                                Asymbols[i][j] = None
                            else:
                                AUinv[i,j] = Asymbols[i][j]
                    BUresult = AUinv*BU
                    C = AL*BUresult-BL
                    for c in C:
                        reducedeqs.append(c)                        
                elif 0:#nummatrixsymbols > 60:
                    # requires swiginac
                    getsubs = lambda valuesubs: self.SubstituteGinacEquations(dictequations, valuesubs, localsymbolmap)
                    # cannot compute inverse since too many symbols
                    log.info('lu decomposition')
                    # PA = L DD**-1 U
                    P, L, DD, U = self.LUdecompositionFF(AU,*self.pvars)
                    log.info('lower triangular solve')
                    res0 = L.lower_triangular_solve(P*BU)
                    # have to use ginac, since sympy is too slow
                    # there are divides in res0, so have to simplify
                    gres1 = swiginac.symbolic_matrix(len(res0),1,'gres1')
                    for i in range(len(res0)):
                        gres0i = GinacUtils.ConvertToGinac(res0[i],localsymbolmap)
                        gDDi = GinacUtils.ConvertToGinac(DD[i,i],localsymbolmap)
                        gres1[i,0] = gres0i*gDDi
                    gothersymbols = [localsymbolmap[s.name] for s in othersymbols if s.name in localsymbolmap]
                    res2 = []
                    gres2 = swiginac.symbolic_matrix(len(res0),1,'gres2')
                    for icol in range(len(gres1)):
                        log.info('extracting poly monoms from L solving: %d', icol)
                        polyterms = GinacUtils.GetPolyTermsFromGinac(gres1[icol],gothersymbols,othersymbols)
                        # create a new symbol for every term
                        eq = S.Zero
                        for monom, coeff in polyterms.iteritems():
                            sym = self.gsymbolgen.next()
                            dictequations.append((sym,coeff))
                            localsymbolmap[sym.name] = swiginac.symbol(sym.name)
                            if sympy_smaller_073:
                                eq += sym*Monomial(*monom).as_expr(*othersymbols)
                            else:
                                eq += sym*Monomial(monom).as_expr(*othersymbols)
                        res2.append(eq)
                        gres2[icol] = GinacUtils.ConvertToGinac(eq,localsymbolmap)
                        
                    gU = GinacUtils.ConvertMatrixToGinac(U,'U',localsymbolmap)
                    log.info('upper triangular solve')
                    gres3 = GinacUtils.SolveUpperTriangular(gU, gres2, 'gres3')
                    res3 = []
                    for icol in range(len(gres3)):
                        log.info('extracting poly monoms from U solving: %d', icol)
                        polyterms = GinacUtils.GetPolyTermsFromGinac(gres3[icol],gothersymbols,othersymbols)
                        # create a new symbol for every term
                        eq = S.Zero
                        for monom, coeff in polyterms.iteritems():
                            sym = self.gsymbolgen.next()
                            dictequations.append((sym,coeff))
                            localsymbolmap[sym.name] = swiginac.symbol(sym.name)
                            if sympy_smaller_073:
                                eq += sym*Monomial(*monom).as_expr(*othersymbols)
                            else:
                                eq += sym*Monomial(monom).as_expr(*othersymbols)
                        res3.append(eq)
                    BUresult = Matrix(gres3.rows(),gres3.cols(),res3)
                    C = AL*BUresult-BL
                    for c in C:
                        reducedeqs.append(c)
                else:
                    # if AU has too many fractions, it can prevent computation
                    allzeros = True
                    for b in BU:
                        if b != S.Zero:
                            allzeros = False
                    if not allzeros:
                        try:
                            AUinv = AU.inv()
                        except ValueError, e:
                            raise self.CannotSolveError(u'failed to invert matrix: %e'%e)
                        
                        BUresult = AUinv*BU
                        C = AL*BUresult-BL
                    else:
                        C = -BL
                    for c in C:
                        if c != S.Zero:
                            reducedeqs.append(c)
            log.info('computed non-singular AU matrix')
        
        if len(reducedeqs) == 0:
            raise self.CannotSolveError('reduced equations are zero')
        
        # is now a (len(neweqs)-len(allmonoms))x1 matrix, usually this is 4x1
        htvars = []
        htvarsubs = []
        htvarsubs2 = []
        usedvars = []
        htvarcossinoffsets = []
        nonhtvars = []
        for iothersymbol, othersymbol in enumerate(othersymbols):
            if othersymbol.name[0] == 'c':
                assert(othersymbols[iothersymbol+1].name[0] == 's')
                htvarcossinoffsets.append(iothersymbol)
                name = othersymbol.name[1:]
                htvar = Symbol('ht%s'%name)
                htvarsubs += [(othersymbol,(1-htvar**2)/(1+htvar**2)),(othersymbols[iothersymbol+1],2*htvar/(1+htvar**2))]
                htvars.append(htvar)
                htvarsubs2.append((Symbol(name),2*atan(htvar)))
                usedvars.append(Symbol(name))
            elif othersymbol.name[0] != 'h' and othersymbol.name[0] != 's':
                # not half-tan, sin, or cos
                nonhtvars.append(othersymbol)
                usedvars.append(othersymbol)
        htvarsubs += [(cvar,(1-tvar**2)/(1+tvar**2)),(svar,2*tvar/(1+tvar**2))]
        htvars.append(tvar)
        htvarsubs2.append((Symbol(varname),2*atan(tvar)))
        usedvars.append(Symbol(varname))
        
        if haszeroequations:
            log.info('special structure in equations detected, try to solve through elimination')
            AllEquations = [eq.subs(self.invsubs) for eq in reducedeqs if self.codeComplexity(eq) < 2000]
            for curvar in usedvars[:-1]:
                try:
                    unknownvars = usedvars[:]
                    unknownvars.remove(curvar)
                    jointtrees2 = []
                    curvarsubs = self.getVariable(curvar).subs
                    treefirst = self.SolveAllEquations(AllEquations, \
                                                       curvars = [curvar], \
                                                       othersolvedvars = self.freejointvars, \
                                                       solsubs = self.freevarsubs[:], \
                                                       endbranchtree = [AST.SolverSequence([jointtrees2])], \
                                                       unknownvars = unknownvars+[tvar], \
                                                       canguessvars = False, \
                                                       currentcases = currentcases, \
                                                       currentcasesubs = currentcasesubs)
                    
                    # solvable, which means we now have len(AllEquations)-1 with two variables, solve with half angles
                    halfanglesolution=self.SolvePairVariablesHalfAngle(raweqns=[eq.subs(curvarsubs) for eq in AllEquations],var0=unknownvars[0],var1=unknownvars[1],othersolvedvars=self.freejointvars+[curvar])[0]
                    # sometimes halfanglesolution can evaluate to all zeros (katana arm), need to catch this and go to a different branch
                    halfanglesolution.AddHalfTanValue = True
                    jointtrees2.append(halfanglesolution)
                    halfanglevar = unknownvars[0] if halfanglesolution.jointname==unknownvars[0].name else unknownvars[1]
                    unknownvars.remove(halfanglevar)
                    
                    try:
                        # give that two variables are solved, can most likely solve the rest. Solving with the original
                        # equations yields simpler solutions since reducedeqs hold half-tangents
                        curvars = solvejointvars[:]
                        curvars.remove(curvar)
                        curvars.remove(halfanglevar)
                        subsinv = []
                        for v in solvejointvars:
                            subsinv += self.getVariable(v).subsinv
                        AllEquationsOrig = [(peq[0].as_expr()-peq[1].as_expr()).subs(subsinv) for peq in rawpolyeqs]
                        self.sortComplexity(AllEquationsOrig)
                        jointtrees2 += self.SolveAllEquations(AllEquationsOrig,curvars=curvars,othersolvedvars=self.freejointvars+[curvar,halfanglevar],solsubs=self.freevarsubs+curvarsubs+self.getVariable(halfanglevar).subs,endbranchtree=endbranchtree, canguessvars=False, currentcases=currentcases, currentcasesubs=currentcasesubs)
                        return preprocesssolutiontree+solutiontree+treefirst,solvejointvars
                    
                    except self.CannotSolveError,e:
                        # try another strategy
                        log.debug(e)
                        
                    # solve all the unknowns now
                    jointtrees3=[]
                    treesecond = self.SolveAllEquations(AllEquations,curvars=unknownvars,othersolvedvars=self.freejointvars+[curvar,halfanglevar],solsubs=self.freevarsubs+curvarsubs+self.getVariable(halfanglevar).subs,endbranchtree=[AST.SolverSequence([jointtrees3])], canguessvars=False, currentcases=currentcases, currentcasesubs=currentcasesubs)
                    for t in treesecond:
                        # most likely t is a solution...
                        t.AddHalfTanValue = True
                        if isinstance(t,AST.SolverCheckZeros):
                            for t2 in t.zerobranch:
                                t2.AddHalfTanValue = True
                            for t2 in t.nonzerobranch:
                                t2.AddHalfTanValue = True
                            if len(t.zerobranch) == 0 or isinstance(t.zerobranch[0],AST.SolverBreak):
                                log.info('detected zerobranch with SolverBreak, trying to fix')
                                
                    jointtrees2 += treesecond
                    # using these solutions, can evaluate all monoms and check for consistency, this step is crucial since
                    # AllEquations might not constrain all degrees of freedom (check out katana)
                    indices = []
                    for i in range(4):
                        monom = [0]*len(symbols)
                        monom[i] = 1
                        indices.append(allmonoms.index(tuple(monom)))
                    if AUinv is not None:
                        X = AUinv*BU
                        for i in [0,2]:
                            jointname=symbols[i].name[1:]
                            try:
                                # atan2(0,0) produces an invalid solution
                                jointtrees3.append(AST.SolverSolution(jointname,jointeval=[atan2(X[indices[i+1]],X[indices[i]])],isHinge=self.IsHinge(jointname)))
                                usedvars.append(Symbol(jointname))
                            except Exception, e:
                                log.warn(e)
                                
                        jointcheckeqs = []
                        for i,monom in enumerate(allmonoms):
                            if not i in indices:
                                eq = S.One
                                for isymbol,ipower in enumerate(monom):
                                    eq *= symbols[isymbol]**ipower
                                jointcheckeqs.append(eq-X[i])
                        # threshold can be a little more loose since just a sanity check
                        jointtrees3.append(AST.SolverCheckZeros('sanitycheck',jointcheckeqs,zerobranch=endbranchtree,nonzerobranch=[AST.SolverBreak('sanitycheck for solveLiWoernleHiller')],anycondition=False,thresh=0.001))
                        return preprocesssolutiontree+solutiontree+treefirst,usedvars
                    else:
                        log.warn('AUinv not initialized, perhaps missing important equations')
                        
                except self.CannotSolveError,e:
                    log.info(e)
                
            try:
                log.info('try to solve first two variables pairwise')
                
                #solution = self.SolvePairVariables(AllEquations,usedvars[0],usedvars[1],self.freejointvars,maxcomplexity=50)
                jointtrees=[]
                unusedvars = [s for s in solvejointvars if not s in usedvars]
                raweqns=[eq for eq in AllEquations if not eq.has(tvar, *unusedvars)]
                if len(raweqns) > 1:
                    halfanglesolution = self.SolvePairVariablesHalfAngle(raweqns=raweqns,var0=usedvars[0],var1=usedvars[1],othersolvedvars=self.freejointvars)[0]
                    halfanglevar = usedvars[0] if halfanglesolution.jointname==usedvars[0].name else usedvars[1]
                    unknownvar = usedvars[1] if halfanglesolution.jointname==usedvars[0].name else usedvars[0]
                    nexttree = self.SolveAllEquations(raweqns,curvars=[unknownvar],othersolvedvars=self.freejointvars+[halfanglevar],solsubs=self.freevarsubs+self.getVariable(halfanglevar).subs,endbranchtree=[AST.SolverSequence([jointtrees])], canguessvars=False, currentcases=currentcases, currentcasesubs=currentcasesubs)
                    #finalsolution = self.solveSingleVariable(AllEquations,usedvars[2],othersolvedvars=self.freejointvars+usedvars[0:2],maxsolutions=4,maxdegree=4)
                    try:
                        finaltree = self.SolveAllEquations(AllEquations,curvars=usedvars[2:],othersolvedvars=self.freejointvars+usedvars[0:2],solsubs=self.freevarsubs+self.getVariable(usedvars[0]).subs+self.getVariable(usedvars[1]).subs,endbranchtree=endbranchtree, canguessvars=False, currentcases=currentcases, currentcasesubs=currentcasesubs)
                        jointtrees += finaltree
                        return preprocesssolutiontree+[halfanglesolution]+nexttree,usedvars
                    
                    except self.CannotSolveError,e:
                        log.debug('failed to solve for final variable %s, so returning just two: %s'%(usedvars[2],str(usedvars[0:2])))
                        jointtrees += endbranchtree
                        # sometimes the last variable cannot be solved, so returned the already solved variables and let the higher function take care of it 
                        return preprocesssolutiontree+[halfanglesolution]+nexttree,usedvars[0:2]
                
            except self.CannotSolveError,e:
                log.debug(u'failed solving first two variables pairwise: %s', e)
                
        if len(reducedeqs) < 3:
            raise self.CannotSolveError('have need at least 3 reducedeqs (%d)'%len(reducedeqs))
        
        log.info('reducing %d equations', len(reducedeqs))
        newreducedeqs = []
        hassinglevariable = False
        for eq in reducedeqs:
            self._CheckPreemptFn(progress = 0.10)
            complexity = self.codeComplexity(eq)
            if complexity > 1500:
                log.warn('equation way too complex (%d), looking for another solution', complexity)
                continue
            
            if complexity > 1500:
                log.info('equation way too complex (%d), so try breaking it down', complexity)
                # don't support this yet...
                eq2 = eq.expand()
                assert(eq2.is_Add)
                log.info('equation has %d additions', len(eq2.args))
                indices = list(range(0, len(eq2.args), 100))
                indices[-1] = len(eq2.args)
                testpolyeqs = []
                startvalue = 0
                for nextvalue in indices[1:]:
                    log.info('computing up to %d', nextvalue)
                    testadd = S.Zero
                    for i in range(startvalue,nextvalue):
                        testadd += eq2.args[i]
                    testpolyeqs.append(Poly(testadd,*othersymbols))
                    startvalue = nextvalue
                # convert each poly's coefficients to symbols
                peq = Poly(S.Zero, *othersymbols)
                for itest, testpolyeq in enumerate(testpolyeqs):
                    log.info('adding equation %d', itest)
                    newpeq = Poly(S.Zero, *othersymbols)
                    for monom, coeff in newpeq.terms():
                        sym = self.gsymbolgen.next()
                        dictequations.append((sym,coeff))
                        if sympy_smaller_073:
                            newpeq += sym*Monomial(*monom).as_expr(*othersymbols)
                        else:
                            newpeq += sym*Monomial(monom).as_expr(*othersymbols)
                    peq += newpeq
            else:
                peq = Poly(eq,*othersymbols)
            maxdenom = [0]*len(htvarcossinoffsets)
            for monoms in peq.monoms():
                for i,ioffset in enumerate(htvarcossinoffsets):
                    maxdenom[i] = max(maxdenom[i],monoms[ioffset]+monoms[ioffset+1])
            eqnew = S.Zero
            for monoms,c in peq.terms():
                term = c
                for i,ioffset in enumerate(htvarcossinoffsets):
                    # for cos
                    num, denom = fraction(htvarsubs[2*i][1])
                    term *= num**monoms[ioffset]
                    # for sin
                    num, denom = fraction(htvarsubs[2*i+1][1])
                    term *= num**monoms[ioffset+1]
                # the denoms for sin/cos of the same joint variable are the same
                for i,ioffset in enumerate(htvarcossinoffsets):
                    denom = fraction(htvarsubs[2*i][1])[1]
                    term *= denom**(maxdenom[i]-monoms[ioffset]-monoms[ioffset+1])
                # multiply the rest of the monoms
                for imonom, monom in enumerate(monoms):
                    if not imonom in htvarcossinoffsets and not imonom-1 in htvarcossinoffsets:
                        # handle non-sin/cos variables yet
                        term *= othersymbols[imonom]**monom
                eqnew += term
            newpeq = Poly(eqnew,htvars+nonhtvars)
            if newpeq != S.Zero:
                newreducedeqs.append(newpeq)
                hassinglevariable |= any([all([__builtin__.sum(monom)==monom[i] for monom in newpeq.monoms()]) for i in range(3)])
        
        if hassinglevariable:
            log.info('hassinglevariable, trying with raw equations')
            AllEquations = []
            for eq in reducedeqs:
                peq = Poly(eq,tvar)
                if sum(peq.degree_list()) == 0:
                    AllEquations.append(peq.TC().subs(self.invsubs).expand())
                elif sum(peq.degree_list()) == 1 and peq.TC() == S.Zero:
                    AllEquations.append(peq.LC().subs(self.invsubs).expand())
                else:
                    # two substitutions: sin/(1+cos), (1-cos)/sin
                    neweq0 = S.Zero
                    neweq1 = S.Zero
                    for monoms,c in peq.terms():
                        neweq0 += c*(svar**monoms[0])*((1+cvar)**(peq.degree(0)-monoms[0]))
                        neweq1 += c*((1-cvar)**monoms[0])*(svar**(peq.degree(0)-monoms[0]))
                    AllEquations.append(neweq0.subs(self.invsubs).expand())
                    AllEquations.append(neweq1.subs(self.invsubs).expand())
            unusedvars = [solvejointvar for solvejointvar in solvejointvars if not solvejointvar in usedvars]
            for eq in AllEquationsExtra:
                #if eq.has(*usedvars) and not eq.has(*unusedvars):
                AllEquations.append(eq)
            self.sortComplexity(AllEquations)

            # first try to solve all the variables at once
            try:
                solutiontree = self.SolveAllEquations(AllEquations,curvars=solvejointvars,othersolvedvars=self.freejointvars[:], solsubs=self.freevarsubs[:], endbranchtree=endbranchtree, canguessvars=False, currentcases=currentcases, currentcasesubs=currentcasesubs)
                return solutiontree, solvejointvars
            except self.CannotSolveError, e:
                log.debug(u'failed solving all variables: %s', e)
                
            try:
                solutiontree = self.SolveAllEquations(AllEquations,curvars=usedvars,othersolvedvars=self.freejointvars[:],solsubs=self.freevarsubs[:], unknownvars=unusedvars, endbranchtree=endbranchtree, canguessvars=False, currentcases=currentcases, currentcasesubs=currentcasesubs)
                return solutiontree, usedvars
            except self.CannotSolveError, e:
                log.debug(u'failed solving used variables: %s', e)
            
            for ivar in range(3):
                try:
                    unknownvars = usedvars[:]
                    unknownvars.pop(ivar)
                    endbranchtree2 = []
                    if 1:
                        solutiontree = self.SolveAllEquations(AllEquations, \
                                                              curvars = [usedvars[ivar]], \
                                                              othersolvedvars = self.freejointvars[:], \
                                                              solsubs = self.freevarsubs[:], \
                                                              endbranchtree = [AST.SolverSequence([endbranchtree2])], \
                                                              unknownvars = unknownvars+unusedvars, \
                                                              canguessvars = False, \
                                                              currentcases = currentcases, \
                                                              currentcasesubs = currentcasesubs)
                        
                        endbranchtree2 += self.SolveAllEquations(AllEquations, \
                                                                 curvars = unknownvars[0:2], \
                                                                 othersolvedvars = self.freejointvars[:]+[usedvars[ivar]], \
                                                                 solsubs = self.freevarsubs[:]+self.getVariable(usedvars[ivar]).subs, \
                                                                 unknownvars = unusedvars, \
                                                                 endbranchtree = endbranchtree, \
                                                                 canguessvars = False, \
                                                                 currentcases = currentcases, \
                                                                 currentcasesubs = currentcasesubs)
                        
                    return preprocesssolutiontree+solutiontree, usedvars#+unusedvars#[unknownvars[1], usedvars[ivar]]#
                
                except self.CannotSolveError, e:
                    log.debug(u'single variable %s failed: %s', usedvars[ivar], e)
                    
#         try:
#             testvars = [Symbol(othersymbols[0].name[1:]),Symbol(othersymbols[2].name[1:]),Symbol(varname)]
#             AllEquations = [(peq[0].as_expr()-peq[1].as_expr()).expand() for peq in polyeqs if not peq[0].has(*symbols)]
#             coupledsolutions = self.SolveAllEquations(AllEquations,curvars=testvars,othersolvedvars=self.freejointvars[:],solsubs=self.freevarsubs[:],endbranchtree=endbranchtree)
#             return coupledsolutions,testvars
#         except self.CannotSolveError:
#             pass
#
        exportcoeffeqs = None
        # only support ileftvar=0 for now
        for ileftvar in [0]:#range(len(htvars)):
            # always take the equations 4 at a time....?
            if len(newreducedeqs) == 3:
                try:
                    exportcoeffeqs,exportmonoms = self.solveDialytically(newreducedeqs, ileftvar, \
                                                                         getsubs = getsubs)
                    break
                except self.CannotSolveError,e:
                    log.warn('failed with leftvar %s: %s',newreducedeqs[0].gens[ileftvar],e)
            else:
                for dialyticeqs in combinations(newreducedeqs,3):
                    try:
                        exportcoeffeqs,exportmonoms = self.solveDialytically(dialyticeqs, ileftvar, \
                                                                             getsubs = getsubs)
                        break
                    except self.CannotSolveError,e:
                        log.warn('failed with leftvar %s: %s',newreducedeqs[0].gens[ileftvar],e)
                
                for dialyticeqs in combinations(newreducedeqs,4):
                    try:
                        exportcoeffeqs,exportmonoms = self.solveDialytically(dialyticeqs, ileftvar, \
                                                                             getsubs = getsubs)
                        break
                    except self.CannotSolveError,e:
                        log.warn('failed with leftvar %s: %s',newreducedeqs[0].gens[ileftvar],e)

                if exportcoeffeqs is None:
                    filteredeqs = [peq for peq in newreducedeqs if peq.degree() <= 2] # has never worked for higher degrees than 2
                    for dialyticeqs in combinations(filteredeqs,6):
                        try:
                            exportcoeffeqs,exportmonoms = self.solveDialytically(dialyticeqs, ileftvar, \
                                                                                 getsubs = getsubs)
                            break
                        except self.CannotSolveError,e:
                            log.warn('failed with leftvar %s: %s',newreducedeqs[0].gens[ileftvar], e)
            if exportcoeffeqs is not None:
                break
        
        self._CheckPreemptFn(progress = 0.11)
        if exportcoeffeqs is None:
            if len(nonhtvars) > 0 and newreducedeqs[0].degree:
                log.info('try to solve one variable in terms of the others')
                doloop = True
                while doloop:
                    doloop = False
                    
                    # check if there is an equation that can solve for nonhtvars[0] easily
                    solvenonhtvareq = None
                    for peq in newreducedeqs:
                        if peq.degree(len(htvars)) == 1:
                            solvenonhtvareq = peq
                            break
                    if solvenonhtvareq is None:
                        break
                    
                    # nonhtvars[0] index is len(htvars)
                    usedvar0solution = solve(newreducedeqs[0],nonhtvars[0])[0]
                    num,denom = fraction(usedvar0solution)
                    igenoffset = len(htvars)
                    # substitute all instances of the variable
                    processedequations = []
                    for peq in newreducedeqs[1:]:
                        if self.codeComplexity(peq.as_expr()) > 10000:
                            log.warn('equation too big')
                            continue
                        maxdegree = peq.degree(igenoffset)
                        eqnew = S.Zero
                        for monoms,c in peq.terms():
                            term = c
                            term *= denom**(maxdegree-monoms[igenoffset])
                            term *= num**(monoms[igenoffset])
                            for imonom, monom in enumerate(monoms):
                                if imonom != igenoffset and imonom < len(htvars):
                                    term *= htvars[imonom]**monom
                            eqnew += term
                        try:
                            newpeq = Poly(eqnew,htvars)
                        except PolynomialError, e:
                            # most likel uservar0solution was bad...
                            raise self.CannotSolveError('equation %s cannot be represented as a polynomial'%eqnew)
                        
                        if newpeq != S.Zero:
                            processedequations.append(newpeq)
                    
                    if len(processedequations) == 0:
                        break
                    
                    # check if any variables have degree <= 1 for all equations
                    for ihtvar,htvar in enumerate(htvars):
                        leftoverhtvars = list(htvars)
                        leftoverhtvars.pop(ihtvar)
                        freeequations = []
                        linearequations = []
                        higherequations = []
                        for peq in processedequations:
                            if peq.degree(ihtvar) == 0:
                                freeequations.append(peq)
                            elif peq.degree(ihtvar) == 1:
                                linearequations.append(peq)
                            else:
                                higherequations.append(peq)
                        if len(freeequations) > 0:
                            log.info('found a way to solve this! still need to implement it though...')
                        elif len(linearequations) > 0 and len(leftoverhtvars) == 1:
                            # try substituting one into the other equations Ax = B
                            A = S.Zero
                            B = S.Zero
                            for monoms,c in linearequations[0].terms():
                                term = c
                                for imonom, monom in enumerate(monoms):
                                    if imonom != ihtvar:
                                        term *= htvars[imonom]**monom
                                if monoms[ihtvar] > 0:
                                    A += term
                                else:
                                    B -= term
                            Apoly = Poly(A,leftoverhtvars)
                            Bpoly = Poly(B,leftoverhtvars)
                            singlepolyequations = []
                            useequations = linearequations[1:]
                            if len(useequations) == 0:
                                useequations += higherequations
                            for peq in useequations:
                                complexity = self.codeComplexity(peq.as_expr())
                                if complexity < 2000:
                                    peqnew = Poly(S.Zero,leftoverhtvars)
                                    maxhtvardegree = peq.degree(ihtvar)
                                    for monoms,c in peq.terms():
                                        term = c
                                        for imonom, monom in enumerate(monoms):
                                            if imonom != ihtvar:
                                                term *= htvars[imonom]**monom
                                        termpoly = Poly(term,leftoverhtvars)
                                        peqnew += termpoly * (Bpoly**(monoms[ihtvar]) * Apoly**(maxhtvardegree-monoms[ihtvar]))
                                    singlepolyequations.append(peqnew)
                            if len(singlepolyequations) > 0:
                                jointsol = 2*atan(leftoverhtvars[0])
                                jointname = leftoverhtvars[0].name[2:]
                                firstsolution = AST.SolverPolynomialRoots(jointname = jointname, \
                                                                          poly = singlepolyequations[0], \
                                                                          jointeval = [jointsol], \
                                                                          isHinge = self.IsHinge(jointname))
                                firstsolution.checkforzeros = []
                                firstsolution.postcheckforzeros = []
                                firstsolution.postcheckfornonzeros = []
                                firstsolution.postcheckforrange = []
                                # in Ax=B, if A is 0 and B is non-zero, then equation is invalid
                                # however if both A and B evaluate to 0, then equation is still valid
                                # therefore equation is invalid only if A==0&&B!=0
                                firstsolution.postcheckforNumDenom = [(A.as_expr(), B.as_expr())]
                                firstsolution.AddHalfTanValue = True

                                # actually both A and B can evaluate to zero, in which case we have to use a different method to solve them
                                AllEquations = []
                                for eq in reducedeqs:
                                    if self.codeComplexity(eq) > 500:
                                        continue
                                    peq = Poly(eq, tvar)
                                    if sum(peq.degree_list()) == 0:
                                        AllEquations.append(peq.TC().subs(self.invsubs).expand())
                                    elif sum(peq.degree_list()) == 1 and peq.TC() == S.Zero:
                                        AllEquations.append(peq.LC().subs(self.invsubs).expand())
                                    else:
                                        # two substitutions: sin/(1+cos), (1-cos)/sin
                                        neweq0 = S.Zero
                                        neweq1 = S.Zero
                                        for monoms,c in peq.terms():
                                            neweq0 += c*(svar**monoms[0])*((1+cvar)**(peq.degree(0)-monoms[0]))
                                            neweq1 += c*((1-cvar)**monoms[0])*(svar**(peq.degree(0)-monoms[0]))
                                        if self.codeComplexity(neweq0) > 1000 or self.codeComplexity(neweq1) > 1000:
                                            break
                                        AllEquations.append(neweq0.subs(self.invsubs).expand())
                                        AllEquations.append(neweq1.subs(self.invsubs).expand())

                                #oldmaxcasedepth = self.maxcasedepth                            
                                try:
                                    #self.maxcasedepth = min(self.maxcasedepth, 2)
                                    solvevar = Symbol(jointname)
                                    curvars = list(usedvars)
                                    curvars.remove(solvevar)
                                    unusedvars = [solvejointvar for solvejointvar in solvejointvars if not solvejointvar in usedvars]
                                    solutiontree = self.SolveAllEquations(AllEquations + AllEquationsExtra, \
                                                                          curvars = curvars+unusedvars, \
                                                                          othersolvedvars = self.freejointvars[:]+[solvevar], \
                                                                          solsubs = self.freevarsubs[:]+self.getVariable(solvevar).subs, \
                                                                          endbranchtree = endbranchtree, \
                                                                          canguessvars = False, \
                                                                          currentcases = currentcases, \
                                                                          currentcasesubs = currentcasesubs)
                                    #secondSolutionComplexity = self.codeComplexity(B) + self.codeComplexity(A)
                                    #if secondSolutionComplexity > 500:
                                    #    log.info('solution for %s is too complex, so delaying its solving')
                                    #solutiontree = self.SolveAllEquations(AllEquations,curvars=curvars,othersolvedvars=self.freejointvars[:]+[solvevar],solsubs=self.freevarsubs[:]+self.getVariable(solvevar).subs,endbranchtree=endbranchtree)
                                    return preprocesssolutiontree + [firstsolution] + solutiontree, usedvars + unusedvars

                                except self.CannotSolveError, e:
                                    log.debug('could not solve full variables from scratch, so use existing solution: %s', e)
                                    secondsolution = AST.SolverSolution(htvar.name[2:], \
                                                                        isHinge = self.IsHinge(htvar.name[2:]))
                                    secondsolution.jointeval = [2*atan2(B.as_expr(), A.as_expr())]
                                    secondsolution.AddHalfTanValue = True
                                    thirdsolution = AST.SolverSolution(nonhtvars[0].name, \
                                                                       isHinge = self.IsHinge(nonhtvars[0].name))
                                    thirdsolution.jointeval = [usedvar0solution]
                                    return preprocesssolutiontree + [firstsolution, \
                                                                     secondsolution, \
                                                                     thirdsolution] + endbranchtree, usedvars
#                               finally:
#                                   self.maxcasedepth = oldmaxcasedepth
            # try to factor the equations manually
            deg1index = None
            for i in range(len(newreducedeqs)):
                if newreducedeqs[i].degree(2) == 1:
                    if self.codeComplexity(newreducedeqs[i].as_expr()) <= 5000:
                        deg1index = i
                        break
                    else:
                        log.warn('found equation with linear DOF, but too complex so skip')
            if deg1index is not None:
                # try to solve one variable in terms of the others
                if len(htvars) > 2:
                    usedvar0solutions = [solve(newreducedeqs[deg1index], htvars[2])[0]]
                    # check which index in usedvars matches htvars[2]
                    for igenoffset in range(len(usedvars)):
                        if htvars[2].name.find(usedvars[igenoffset].name) >= 0:
                            break
                    polyvars = htvars[0:2]
                elif len(htvars) > 1:
                    usedvar0solutions = solve(newreducedeqs[deg1index], htvars[1])
                    igenoffset = 1
                    polyvars = htvars[0:1] + nonhtvars
                else:
                    usedvar0solutions = []
                processedequations = []
                if len(usedvar0solutions) > 0:
                    usedvar0solution = usedvar0solutions[0]
                    num,denom = fraction(usedvar0solution)
                    # substitute all instances of the variable
                    
                    for ipeq, peq in enumerate(newreducedeqs):
                        if ipeq == deg1index:
                            continue
                        newpeq = S.Zero
                        if peq.degree(igenoffset) > 1:
                            # ignore higher powers
                            continue
                        elif peq.degree(igenoffset) == 0:
                            newpeq = Poly(peq,*polyvars)
                        else:
                            maxdegree = peq.degree(igenoffset)
                            eqnew = S.Zero
                            for monoms, c in peq.terms():
                                term = c*denom**(maxdegree-monoms[igenoffset])
                                term *= num**(monoms[igenoffset])
                                for imonom, monom in enumerate(monoms):
                                    if imonom != igenoffset:
                                        term *= peq.gens[imonom]**monom
                                eqnew += term.expand()
                            try:
                                newpeq = Poly(eqnew, *polyvars)
                            except PolynomialError, e:
                                # most likel uservar0solution was bad
                                raise self.CannotSolveError('equation %s cannot be represented as a polynomial'%eqnew)

                        if newpeq != S.Zero:
                            # normalize by the greatest coefficient in LC, or otherwise determinant will never succeed
                            LC=newpeq.LC()
                            highestcoeff = None
                            if LC.is_Add:
                                for arg in LC.args:
                                    coeff = None
                                    if arg.is_Mul:
                                        coeff = S.One
                                        for subarg in arg.args:
                                            if subarg.is_number:
                                                coeff *= abs(subarg)
                                    elif arg.is_number:
                                        coeff = abs(arg)
                                    if coeff is not None:
                                        if coeff > S.One:
                                            # round to the nearest integer
                                            coeff = int(round(coeff.evalf()))
                                        if highestcoeff is None or coeff > highestcoeff:
                                            highestcoeff = coeff
                            if highestcoeff.has(oo) or highestcoeff.has(-oo):
                                log.warn('an equation has inifinity?!')
                            else:
                                if highestcoeff is not None:
                                    processedequations.append(newpeq*(S.One/highestcoeff))
                                else:
                                    processedequations.append(newpeq)
                        else:
                            log.info('equation is zero, so ignoring')
                            
                for dialyticeqs in combinations(processedequations,3):
                    Mall = None
                    leftvar = None
                    for ileftvar in range(2):
                        # TODO, sometimes this works and sometimes this doesn't
                        try:
                            Mall, allmonoms = self.solveDialytically(dialyticeqs,ileftvar,returnmatrix=True)
                            if Mall is not None:
                                leftvar=processedequations[0].gens[ileftvar]
                                break
                        except self.CannotSolveError, e:
                            log.debug(e)
                    if Mall is None:
                        continue
                    log.info('success in solving sub-coeff matrix!')
                    shape=Mall[0].shape
                    Malltemp = [None]*len(Mall)
                    M = zeros(shape)
                    dictequations2 = list(dictequations)
                    for idegree in range(len(Mall)):
                        Malltemp[idegree] = zeros(shape)
                        for i in range(shape[0]):
                            for j in range(shape[1]):
                                if Mall[idegree][i,j] != S.Zero:
                                    sym = self.gsymbolgen.next()
                                    Malltemp[idegree][i,j] = sym
                                    dictequations2.append((sym,Mall[idegree][i,j]))
                        M += Malltemp[idegree]*leftvar**idegree
                    tempsymbols = [self.gsymbolgen.next() for i in range(len(M))]
                    tempsubs = []
                    for i in range(len(tempsymbols)):
                        if M[i] != S.Zero:
                            tempsubs.append((tempsymbols[i],Poly(M[i],leftvar)))
                        else:
                            tempsymbols[i] = S.Zero
                    Mtemp = Matrix(M.shape[0],M.shape[1],tempsymbols)                    
                    dettemp=Mtemp.det()
                    log.info('multiplying all determinant coefficients for solving %s',leftvar)
                    eqadds = []
                    for arg in dettemp.args:
                        eqmuls = [Poly(arg2.subs(tempsubs),leftvar) for arg2 in arg.args]
                        if sum(eqmuls[0].degree_list()) == 0:
                            eq = eqmuls.pop(0)
                            eqmuls[0] = eqmuls[0]*eq
                        while len(eqmuls) > 1:
                            ioffset = 0
                            eqmuls2 = []
                            while ioffset < len(eqmuls)-1:
                                eqmuls2.append(eqmuls[ioffset]*eqmuls[ioffset+1])
                                ioffset += 2
                            eqmuls = eqmuls2
                        eqadds.append(eqmuls[0])
                    det = Poly(S.Zero,leftvar)
                    for eq in eqadds:
                        det += eq
                        
                    jointsol = 2*atan(leftvar)
                    firstsolution = AST.SolverPolynomialRoots(jointname = usedvars[ileftvar].name, \
                                                              poly = det, \
                                                              jointeval = [jointsol], \
                                                              isHinge = self.IsHinge(usedvars[ileftvar].name))
                    firstsolution.checkforzeros = []
                    firstsolution.postcheckforzeros = []
                    firstsolution.postcheckfornonzeros = []
                    firstsolution.postcheckforrange = []
                    firstsolution.dictequations = dictequations2
                    firstsolution.AddHalfTanValue = True
                    
                    # just solve the lowest degree one
                    complexity = [(eq.degree(1-ileftvar)*100000 + \
                                   self.codeComplexity(eq.as_expr()),eq) for eq in processedequations if eq.degree(1-ileftvar) > 0]
                    complexity.sort(key = itemgetter(0))
                    
                    orderedequations = [peq for c,peq in complexity]
                    jointsol = 2*atan(htvars[1-ileftvar])
                    secondsolution = AST.SolverPolynomialRoots(jointname = usedvars[1-ileftvar].name, \
                                                               poly = Poly(orderedequations[0], htvars[1-ileftvar]), \
                                                               jointeval = [jointsol], \
                                                               isHinge = self.IsHinge(usedvars[1-ileftvar].name))
                    secondsolution.checkforzeros = []
                    secondsolution.postcheckforzeros = []
                    secondsolution.postcheckfornonzeros = []
                    secondsolution.postcheckforrange = []
                    secondsolution.AddHalfTanValue = True
                    
                    thirdsolution = AST.SolverSolution(usedvars[2].name, \
                                                       isHinge = self.IsHinge(usedvars[2].name))
                    thirdsolution.jointeval = [usedvar0solution]
                    return preprocesssolutiontree + [firstsolution, \
                                                     secondsolution, \
                                                     thirdsolution] + endbranchtree, usedvars

                    
            raise self.CannotSolveError('failed to solve dialytically')

        if 0:
            # quadratic equations
            iquadvar = 1
            quadpoly0 = Poly(newreducedeqs[0].as_expr(), htvars[iquadvar])
            quadpoly1 = Poly(newreducedeqs[2].as_expr(), htvars[iquadvar])
            a0, b0, c0 = quadpoly0.coeffs()
            a1, b1, c1 = quadpoly1.coeffs()
            quadsolnum = (-a1*c0 + a0*c1).expand()
            quadsoldenom = (-a1*b0 + a0*b1).expand()
            
        if ileftvar > 0:
            raise self.CannotSolveError(('solving equations dialytically succeeded with var index %d, ' + \
                                         'unfortunately code generation supports only index 0') % ileftvar)
        
        exportvar = [htvars[ileftvar].name]
        exportvar += [v.name for i,v in enumerate(htvars) if i != ileftvar]
        exportfnname = 'solvedialyticpoly12qep' if len(exportmonoms) == 9 else 'solvedialyticpoly8qep'
        coupledsolution = AST.SolverCoeffFunction(jointnames = [v.name for v in usedvars], \
                                                  jointeval = [v[1] for v in htvarsubs2], \
                                                  jointevalcos = [htvarsubs[2*i][1] for i in range(len(htvars))], \
                                                  jointevalsin = [htvarsubs[2*i+1][1] for i in range(len(htvars))], \
                                                  isHinges = [self.IsHinge(v.name) for v in usedvars], \
                                                  exportvar = exportvar, \
                                                  exportcoeffeqs = exportcoeffeqs, \
                                                  exportfnname = exportfnname, \
                                                  rootmaxdim = 16)
        coupledsolution.presetcheckforzeros = checkforzeros
        coupledsolution.dictequations = dictequations
        solutiontree.append(coupledsolution)
        self.usinglapack = True

        if 0:
            if currentcases is None:
                currentcases = set()
            if currentcasesubs is None:
                currentcasesubs = list()
            
            rotsymbols = set(self.Tee[:3,:3])
            possiblesub = [(self.Tee[1,2], S.Zero)]
            possiblesub2 = [(self.Tee[2,2], S.Zero)]
            possiblevar,possiblevalue = possiblesub[0]
            possiblevar2,possiblevalue2 = possiblesub2[0]
            cond = Abs(possiblevar-possiblevalue.evalf(n=30))
            evalcond = Abs(fmod(possiblevar-possiblevalue+pi,2*pi)-pi)# + evalcond
            cond2 = Abs(possiblevar2-possiblevalue2.evalf(n=30))
            evalcond2 = Abs(fmod(possiblevar2-possiblevalue2+pi,2*pi)-pi)# + evalcond
            if self._iktype == 'transform6d' and \
               possiblevar in rotsymbols and \
               possiblevalue == S.Zero and \
               possiblevar2 in rotsymbols and \
               possiblevalue2 == S.Zero:
                
                checkexpr = [[cond+cond2], evalcond+evalcond2, possiblesub+possiblesub2, []]
                #flatzerosubstitutioneqs.append(checkexpr)
                #localsubstitutioneqs.append(checkexpr)
                #handledconds.append(cond+cond2)
                row1 = int(possiblevar.name[-2])
                col1 = int(possiblevar.name[-1])
                row2 = int(possiblevar2.name[-2])
                col2 = int(possiblevar2.name[-1])
                row3 = 3 - row1 - row2
                col3 = 3 - col1 - col2
                if row1 == row2:
                    # (row1, col3) is either 1 or -1, but don't know which.
                    # know that (row1+1,col3) and (row1+2,col3) are zero though...
                    checkexpr[2].append((Symbol('%s%d%d'%(possiblevar.name[:-2], (row2+1)%3, col3)), S.Zero))
                    checkexpr[2].append((Symbol('%s%d%d'%(possiblevar.name[:-2], (row1+2)%3, col3)), S.Zero))
                    checkexpr[2].append((Symbol('%s%d%d'%(possiblevar.name[:-2], row1, col3))**2, S.One)) # squared in the corner should always be 1
                    checkexpr[2].append((Symbol('%s%d%d'%(possiblevar.name[:-2], row1, col3))**3, Symbol('%s%d%d'%(possiblevar.name[:-2], row1, col3)))) # squared in the corner should always be 1
                    # furthermore can defer that the left over 4 values are [cos(ang), sin(ang), cos(ang), -sin(ang)] = abcd
                    if row1 == 1:
                        minrow = 0
                        maxrow = 2
                    else:
                        minrow = (row1+1)%3
                        maxrow = (row1+2)%3
                    ra = Symbol('%s%d%d'%(possiblevar.name[:-2], minrow, col1))
                    rb = Symbol('%s%d%d'%(possiblevar.name[:-2], minrow, col2))
                    rc = Symbol('%s%d%d'%(possiblevar.name[:-2], maxrow, col1))
                    rd = Symbol('%s%d%d'%(possiblevar.name[:-2], maxrow, col2))
                    checkexpr[2].append((rb**2, S.One-ra**2))
                    checkexpr[2].append((rb**3, rb-rb*ra**2)) # need 3rd power since sympy cannot divide out the square
                    checkexpr[2].append((rc**2, S.One-ra**2))
                    #checkexpr[2].append((rc, -rb)) # not true
                    #checkexpr[2].append((rd, ra)) # not true
                elif col1 == col2:
                    # (row3, col1) is either 1 or -1, but don't know which.
                    # know that (row3,col1+1) and (row3,col1+2) are zero though...
                    checkexpr[2].append((Symbol('%s%d%d'%(possiblevar.name[:-2], row3, (col1+1)%3)), S.Zero))
                    checkexpr[2].append((Symbol('%s%d%d'%(possiblevar.name[:-2], row3, (col1+2)%3)), S.Zero))
                    checkexpr[2].append((Symbol('%s%d%d'%(possiblevar.name[:-2], row3, col1))**2, S.One)) # squared in the corner should always be 1
                    checkexpr[2].append((Symbol('%s%d%d'%(possiblevar.name[:-2], row3, col1))**3, Symbol('%s%d%d'%(possiblevar.name[:-2], row3, col1)))) # squared in the corner should always be 1
                    # furthermore can defer that the left over 4 values are [cos(ang), sin(ang), cos(ang), -sin(ang)] = abcd
                    if col1 == 1:
                        mincol = 0
                        maxcol = 2
                    else:
                        mincol = (col1+1)%3
                        maxcol = (col1+2)%3
                    ra = Symbol('%s%d%d'%(possiblevar.name[:-2], row1, mincol))
                    rb = Symbol('%s%d%d'%(possiblevar.name[:-2], row2, mincol))
                    rc = Symbol('%s%d%d'%(possiblevar.name[:-2], row1, maxcol))
                    rd = Symbol('%s%d%d'%(possiblevar.name[:-2], row2, maxcol))
                    checkexpr[2].append((rb**2, S.One-ra**2))
                    checkexpr[2].append((rb**3, rb-rb*ra**2)) # need 3rd power since sympy cannot divide out the square
                    checkexpr[2].append((rc**2, S.One-ra**2))
                    #checkexpr[2].append((rc, -rb)) # not true
                    #checkexpr[2].append((rd, ra)) # not true

        return preprocesssolutiontree + solutiontree + endbranchtree,usedvars

    def ConvertSinCosEquationToHalfTan(self, eq, convertvars):
        """
        Converts all the sin/cos of variables to half-tangents. Returns two equations (poly, denominator)

        Called by solveLiWoernleHiller.
        """
        cossinvars = []
        htvarsubs = []
        htvars = []
        htvarsubsinv = []
        cossinsubs = []
        for varsym in convertvars:
            var = self.getVariable(varsym)
            cossinvars.append(var.cvar)
            cossinvars.append(var.svar)
            htvar = Symbol('ht%s'%varsym.name)
            htvarsubs += [(var.cvar,(1-htvar**2)/(1+htvar**2)),(var.svar,2*htvar/(1+htvar**2))]
            htvarsubsinv.append((htvar, (1-var.cvar)/var.svar))
            htvars.append(htvar)
            cossinsubs.append((cos(varsym), var.cvar))
            cossinsubs.append((sin(varsym), var.svar))
        peq = Poly(eq.subs(cossinsubs),*cossinvars)
        maxdenom = [0]*len(convertvars)
        for monoms in peq.monoms():
            for i in range(len(convertvars)):
                maxdenom[i] = max(maxdenom[i],monoms[2*i]+monoms[2*i+1])
        eqnew = S.Zero
        for monoms,c in peq.terms():
            term = c
            for i in range(len(convertvars)):
                # for cos
                num, denom = fraction(htvarsubs[2*i][1])
                term *= num**monoms[2*i]
                # for sin
                num, denom = fraction(htvarsubs[2*i+1][1])
                term *= num**monoms[2*i+1]
            # the denoms for sin/cos of the same joint variable are the same
            for i in range(len(convertvars)):
                denom = fraction(htvarsubs[2*i][1])[1]
                exp = maxdenom[i] - monoms[2*i] - monoms[2*i+1]
                if exp > 0:
                    term *= denom**exp
            eqnew += term
        #newpeq = Poly(eqnew,htvars)
        othereq = S.One
        for i in range(len(convertvars)):
            othereq *= (1+htvars[i]**2)**maxdenom[i]
        return eqnew, othereq, htvarsubsinv

    
    def solveKohliOsvatic(self, \
                          rawpolyeqs, \
                          solvejointvars, \
                          endbranchtree, \
                          AllEquationsExtra = None, \
                          currentcases = None, \
                          currentcasesubs = None):
        """
        Find a 16x16 matrix where the entries are linear in tan(tj/2) [Kohli1993]. Takes in the 14 Raghavan-Roth equations.
        
        [Kohli1993] Dilip Kohli and M. Osvatic, "Inverse Kinematics of General 6R and 5R,P Serial Manipulators", 
        Journal of Mechanical Design, Volume 115, Issue 4, Dec 1993.

        Called by solveFullIK_6DGeneral and solveFullIK_TranslationDirection5D
        """
        
        log.info('Attempt the Kohli-Osvatic general IK method')
        if len(rawpolyeqs[0][0].gens) < len(rawpolyeqs[0][1].gens):
            for peq in rawpolyeqs:
                peq[0], peq[1] = peq[1], peq[0]

        symbols = list(rawpolyeqs[0][0].gens)
        othersymbols = list(rawpolyeqs[0][1].gens)
        othersymbolsnames = []
        for s in othersymbols:
            testeq = s.subs(self.invsubs)
            for solvejointvar in solvejointvars:
                if testeq.has(solvejointvar):
                    othersymbolsnames.append(solvejointvar)
                    break
        assert(len(othersymbols)==len(othersymbolsnames))
        symbolsubs = [(symbols[i].subs(self.invsubs),symbols[i]) for i in range(len(symbols))]
        if len(symbols) != 6:
            raise self.CannotSolveError('Kohli/Osvatic method requires 3 unknown variables')
            
        # choose which leftvar can determine the singularity of the following equations!
        for i in range(0,6,2):
            eqs = [peq for peq in rawpolyeqs if peq[0].has(symbols[i],symbols[i+1])]
            if len(eqs) <= 8:
                break
        if len(eqs) > 8:
            raise self.CannotSolveError('Need <=8 equations in one variable; currently there are %d' % len(eqs))
        
        cvar = symbols[i]
        svar = symbols[i+1]
        tvar = Symbol('t'+cvar.name[1:])
        symbols.remove(cvar)
        symbols.remove(svar)
        othereqs = [peq for peq in rawpolyeqs if not peq[0].has(cvar,svar)]

        polyeqs = [[eq[0].as_expr(),eq[1]] for eq in eqs]
        if len(polyeqs) < 8:
            raise self.CannotSolveError('solveKohliOsvatic: need 8 or more polyeqs')

        # solve the othereqs for symbols without the standalone symbols[2] and symbols[3]
        reducedeqs = []
        othersymbolsnamesunique = list(set(othersymbolsnames)) # get the unique names
        for jother in range(len(othersymbolsnamesunique)):
            if not self.IsHinge(othersymbolsnamesunique[jother].name):
                continue
            othervar=self.getVariable(othersymbolsnamesunique[jother])
            cosmonom = [0]*len(othersymbols)
            cosmonom[othersymbols.index(othervar.cvar)] = 1
            cosmonom = tuple(cosmonom)
            sinmonom = [0]*len(othersymbols)
            sinmonom[othersymbols.index(othervar.svar)] = 1
            sinmonom = tuple(sinmonom)
            leftsideeqs = []
            rightsideeqs = []
            finaleqsymbols = symbols + [othervar.cvar,othervar.svar]
            for eq0,eq1 in othereqs:
                leftsideeq = Poly(eq1,*othersymbols)
                leftsideeqdict = leftsideeq.as_dict()
                rightsideeq = Poly(eq0,*finaleqsymbols)
                coscoeff = leftsideeqdict.get(cosmonom,S.Zero)
                if coscoeff != S.Zero:
                    rightsideeq = rightsideeq - othervar.cvar*coscoeff
                    leftsideeq = leftsideeq - othervar.cvar*coscoeff
                sincoeff = leftsideeqdict.get(sinmonom,S.Zero)
                if sincoeff != S.Zero:
                    rightsideeq = rightsideeq - othervar.svar*sincoeff
                    leftsideeq = leftsideeq - othervar.svar*sincoeff
                const = leftsideeq.TC()
                if const != S.Zero:
                    rightsideeq = rightsideeq - const
                    leftsideeq = leftsideeq - const
                # check that leftsideeq doesn't hold any terms with cosmonom and sinmonom?
                rightsideeqs.append(rightsideeq)
                leftsideeqs.append(leftsideeq)
            # number of symbols for kawada-hiro robot is 16
            if len(othersymbols) > 2:
                reducedeqs = self.reduceBothSidesSymbolically(leftsideeqs,rightsideeqs,usesymbols=False,maxsymbols=18)
                for peq in reducedeqs:
                    peq[0] = Poly(peq[0],*othersymbols)
            else:
                reducedeqs = [[left,right] for left,right in izip(leftsideeqs,rightsideeqs)]
            if len(reducedeqs) > 0:
                break
            
        if len(reducedeqs) == 0:
            raise self.CannotSolveError('KohliOsvatic method: could not reduce the equations')

        finaleqs = []
        for peq0,eq1 in reducedeqs:
            if peq0 == S.Zero:
                finaleqs.append(Poly(eq1,*finaleqsymbols))

        if len(finaleqs) >= 2:
            # perhaps can solve finaleqs as is?
            # transfer othersymbols[2*jother:(2+2*jother)] to the leftside
            try:
                leftsideeqs = []
                rightsideeqs = []
                for finaleq in finaleqs:
                    peq=Poly(finaleq,*othersymbols[2*jother:(2+2*jother)])
                    leftsideeqs.append(peq.sub(peq.TC()))
                    rightsideeqs.append(-peq.TC())
                reducedeqs2 = self.reduceBothSidesSymbolically(leftsideeqs,rightsideeqs,usesymbols=False,maxsymbols=18)
                # find all the equations with left side = to zero
                usedvars = set()
                for symbol in symbols:
                    usedvars.add(Symbol(symbol.name[1:]))
                AllEquations = []
                for eq0, eq1 in reducedeqs2:
                    if eq0 == S.Zero:
                        AllEquations.append(eq1.subs(self.invsubs))
                if len(AllEquations) > 0:
                    otherjointtrees = []
                    tree = self.SolveAllEquations(AllEquations,curvars=list(usedvars),othersolvedvars=[],solsubs=self.freevarsubs,endbranchtree=[AST.SolverSequence([otherjointtrees])], canguessvars=False, currentcases=currentcases, currentcasesubs=currentcasesubs)
                    log.info('first SolveAllEquations successful: %s',usedvars)
#                     try:
#                         # although things can be solved at this point, it yields a less optimal solution than if all variables were considered...
#                         solsubs=list(self.freevarsubs)
#                         for usedvar in usedvars:
#                             solsubs += self.getVariable(usedvar).subs
#                         # solved, so substitute back into reducedeqs and see if anything new can be solved
#                         otherusedvars = set()
#                         for symbol in othersymbols:
#                             otherusedvars.add(Symbol(symbol.name[1:]))
#                         OtherAllEquations = []
#                         for peq0,eq1 in reducedeqs:
#                             OtherAllEquations.append((peq0.as_expr()-eq1).subs(self.invsubs).expand())
#                         otherjointtrees += self.SolveAllEquations(OtherAllEquations,curvars=list(otherusedvars),othersolvedvars=list(usedvars),solsubs=solsubs,endbranchtree=endbranchtree)
#                         return tree, list(usedvars)+list(otherusedvars)
#                     except self.CannotSolveError:
                        # still have the initial solution
                    otherjointtrees += endbranchtree
                    return tree, list(usedvars)
                
            except self.CannotSolveError,e:
                pass
        
        log.info('Build final equations for symbols: %s', finaleqsymbols)
        neweqs=[]
        for i in range(0,8,2):
            p0 = Poly(polyeqs[i][0],cvar,svar)
            p0dict = p0.as_dict()
            p1 = Poly(polyeqs[i+1][0],cvar,svar)
            p1dict = p1.as_dict()
            r0 = polyeqs[i][1].as_expr()
            r1 = polyeqs[i+1][1].as_expr()
            if self.equal(p0dict.get((1,0),S.Zero),-p1dict.get((0,1),S.Zero)) and self.equal(p0dict.get((0,1),S.Zero),p1dict.get((1,0),S.Zero)):
                p0,p1 = p1,p0
                p0dict,p1dict=p1dict,p0dict
                r0,r1 = r1,r0
            if self.equal(p0dict.get((1,0),S.Zero),p1dict.get((0,1),S.Zero)) and self.equal(p0dict.get((0,1),S.Zero),-p1dict.get((1,0),S.Zero)):
                # p0+tvar*p1, p1-tvar*p0
                # subs: tvar*svar + cvar = 1, svar-tvar*cvar=tvar
                neweqs.append([Poly(p0dict.get((1,0),S.Zero) + p0dict.get((0,1),S.Zero)*tvar + p0.TC() + tvar*p1.TC(),*symbols), Poly(r0+tvar*r1,*othersymbols)])
                neweqs.append([Poly(p0dict.get((1,0),S.Zero)*tvar - p0dict.get((0,1),S.Zero) - p0.TC()*tvar + p1.TC(),*symbols), Poly(r1-tvar*r0,*othersymbols)])
        if len(neweqs) != 8:
            raise self.CannotSolveError('Coefficients of equations need to match! only got %d reduced equations' % \
                                        len(neweqs))
    
        for eq0,eq1 in neweqs:
            commondenom = Poly(S.One,*self.pvars)
            hasunknown = False
            for m,c in eq1.terms():
                foundreq = [req[1] for req in reducedeqs if req[0].monoms()[0] == m]
                if len(foundreq) > 0:
                    n,d = fraction(foundreq[0])
                    commondenom = Poly(lcm(commondenom,d),*self.pvars)
                else:
                    if m[2*(1-jother)] > 0 or m[2*(1-jother)+1] > 0:
                        # perhaps there's a way to combine what's in reducedeqs?
                        log.warn('unknown %s',m)
                        hasunknown = True
            if hasunknown:
                continue
            commondenom = self.removecommonexprs(commondenom.as_expr(), onlygcd = True)
            finaleq = eq0.as_expr()*commondenom
            for m,c in eq1.terms():
                foundreq = [req[1] for req in reducedeqs if req[0].monoms()[0] == m]
                if len(foundreq) > 0:
                    finaleq = finaleq - c*simplify(foundreq[0]*commondenom)
                else:
                    finaleq = finaleq - Poly.from_dict({m:c*commondenom},*eq1.gens).as_expr()
            finaleqs.append(Poly(finaleq.expand(),*finaleqsymbols))
                
        # finally do the half angle substitution with symbols
        # set:
        # j=othersymbols[2]*(1+dummys[0]**2)*(1+dummys[1]**2)
        # k=othersymbols[3]*(1+dummys[0]**2)*(1+dummys[1]**2)
        dummys = []
        dummysubs = []
        dummysubs2 = []
        dummyvars = []
        usedvars = []

        dummys.append(tvar)
        dummyvars.append((tvar,tan(0.5*Symbol(tvar.name[1:]))))
        usedvars.append(Symbol(cvar.name[1:]))
        dummysubs2.append((usedvars[-1],2*atan(tvar)))
        dummysubs += [(cvar,(1-tvar**2)/(1+tvar**2)),(svar,2*tvar/(1+tvar**2))]

        for i in range(0,len(symbols),2):
            dummy = Symbol('ht%s'%symbols[i].name[1:])
            # [0] - cos, [1] - sin
            dummys.append(dummy)
            dummysubs += [(symbols[i],(1-dummy**2)/(1+dummy**2)),(symbols[i+1],2*dummy/(1+dummy**2))]
            var = symbols[i].subs(self.invsubs).args[0]
            dummyvars.append((dummy,tan(0.5*var)))
            dummysubs2.append((var,2*atan(dummy)))
            if not var in usedvars:
                usedvars.append(var)
        commonmult = (1+dummys[1]**2)*(1+dummys[2]**2)

        usedvars.append(Symbol(othersymbols[2*jother].name[1:]))
        dummyj = Symbol('dummyj')
        dummyk = Symbol('dummyk')
        dummyjk = Symbol('dummyjk')

        dummys.append(dummyj)
        dummyvars.append((dummyj,othersymbols[2*jother]*(1+dummyvars[1][1]**2)*(1+dummyvars[2][1]**2)))
        dummysubs.append((othersymbols[2*jother],cos(dummyjk)))        
        dummys.append(dummyk)
        dummyvars.append((dummyk,othersymbols[1+2*jother]*(1+dummyvars[1][1]**2)*(1+dummyvars[2][1]**2)))
        dummysubs.append((othersymbols[1+2*jother],sin(dummyjk)))
        dummysubs2.append((usedvars[-1],dummyjk))

        newreducedeqs = []
        for peq in finaleqs:
            eqnew = S.Zero
            for monoms,c in peq.terms():
                term = S.One
                for i in range(4):
                    term *= dummysubs[i+2][1]**monoms[i]
                if monoms[4] == 1:
                    eqnew += c * dummyj
                elif monoms[5] == 1:
                    eqnew += c * dummyk
                else:
                    eqnew += c*simplify(term*commonmult)
            newreducedeqs.append(Poly(eqnew,*dummys))

        exportcoeffeqs = None
        for ileftvar in range(len(dummys)):
            leftvar = dummys[ileftvar]
            try:
                exportcoeffeqs,exportmonoms = self.solveDialytically(newreducedeqs,ileftvar,getsubs=None)
                break
            except self.CannotSolveError,e:
                log.warn('failed with leftvar %s: %s',leftvar,e)

        if exportcoeffeqs is None:
            raise self.CannotSolveError('failed to solve dialytically')
        if ileftvar > 0:
            raise self.CannotSolveError('solving equations dialytically succeeded with var index %d, unfortunately code generation supports only index 0'%ileftvar)
    
        coupledsolution = AST.SolverCoeffFunction(jointnames = [v.name for v in usedvars], \
                                                  jointeval = [v[1] for v in dummysubs2], \
                                                  jointevalcos = [dummysubs[2*i][1] for i in range(len(usedvars))], \
                                                  jointevalsin = [dummysubs[2*i+1][1] for i in range(len(usedvars))], \
                                                  isHinges = [self.IsHinge(v.name) for v in usedvars], \
                                                  exportvar = dummys[0:3]+[dummyjk], \
                                                  exportcoeffeqs = exportcoeffeqs, \
                                                  exportfnname = 'solvedialyticpoly16lep', \
                                                  rootmaxdim = 16)
        self.usinglapack = True
        return [coupledsolution]+endbranchtree, usedvars

    def solveDialytically(self, dialyticeqs, ileftvar, \
                          returnmatrix = False, \
                          getsubs = None):
        """ 
        Return the coefficients to solve equations dialytically (Salmon 1885) leaving out variable index ileftvar.

        Extract the coefficients of 1, leftvar**1, leftvar**2, ... of every equation

        Every len(dialyticeqs)*len(monoms) coefficients specify one degree of all the equations 
        (order of monoms is specified in exportmonomorder

        There should be len(dialyticeqs)*len(monoms)*maxdegree coefficients

        Method also checks if the equations are linearly dependent.

        Called by solveManochaCanny, solveLiWoernleHiller, solveKohliOsvatic, SolvePairVariablesHalfAngle
        """
        self._CheckPreemptFn(progress = 0.12)
        if len(dialyticeqs) == 0:
            raise self.CannotSolveError('solveDialytically given zero equations')
        
        allmonoms = set()
        origmonoms = set()
        maxdegree = 0
        leftvar = dialyticeqs[0].gens[ileftvar]
        extradialyticeqs = []
        for peq in dialyticeqs:
            if sum(peq.degree_list()) == 0:
                log.warn('solveDialytically: polynomial %s degree is 0',peq)
                continue
            for m in peq.monoms():
                mlist = list(m)
                maxdegree=max(maxdegree,mlist.pop(ileftvar))
                allmonoms.add(tuple(mlist))
                origmonoms.add(tuple(mlist))
                mlist[0] += 1
                allmonoms.add(tuple(mlist))
            
            # check if any monoms are not expressed in this poly
            # if so, add another poly with the monom multiplied, will this give bad solutions?
            for igen in range(len(peq.gens)):
                if all([m[igen]==0 for m in peq.monoms()]):
                    log.debug('adding extra equation multiplied by %s', peq.gens[igen])
                    extradialyticeqs.append(peq*peq.gens[igen])
                    # multiply by peq.gens[igen]
                    for m in peq.monoms():
                        mlist = list(m)
                        mlist[igen] += 1
                        maxdegree=max(maxdegree,mlist.pop(ileftvar))
                        allmonoms.add(tuple(mlist))
                        origmonoms.add(tuple(mlist))
                        mlist[0] += 1
                        allmonoms.add(tuple(mlist))
        
        dialyticeqs = list(dialyticeqs) + extradialyticeqs # dialyticeqs could be a tuple
        allmonoms = list(allmonoms)
        allmonoms.sort()
        origmonoms = list(origmonoms)
        origmonoms.sort()
        if len(allmonoms)<2*len(dialyticeqs):
            log.warn('solveDialytically equations %d > %d, should be equal...', 2*len(dialyticeqs),len(allmonoms))
            # TODO not sure how to select the equations
            N = len(allmonoms)/2
            dialyticeqs = dialyticeqs[:N]
        if len(allmonoms) == 0 or len(allmonoms)>2*len(dialyticeqs):
            raise self.CannotSolveError('solveDialytically: more unknowns than equations %d>%d'%(len(allmonoms), 2*len(dialyticeqs)))
        
        Mall = [zeros((2*len(dialyticeqs),len(allmonoms))) for i in range(maxdegree+1)]
        Mallindices = [-ones((2*len(dialyticeqs),len(allmonoms))) for i in range(maxdegree+1)]
        exportcoeffeqs = [S.Zero]*(len(dialyticeqs)*len(origmonoms)*(maxdegree+1))
        for ipeq,peq in enumerate(dialyticeqs):
            for m,c in peq.terms():
                mlist = list(m)
                degree=mlist.pop(ileftvar)
                exportindex = degree*len(origmonoms)*len(dialyticeqs) + len(origmonoms)*ipeq+origmonoms.index(tuple(mlist))
                assert(exportcoeffeqs[exportindex] == S.Zero)
                exportcoeffeqs[exportindex] = c
                Mall[degree][len(dialyticeqs)+ipeq,allmonoms.index(tuple(mlist))] = c
                Mallindices[degree][len(dialyticeqs)+ipeq,allmonoms.index(tuple(mlist))] = exportindex
                mlist[0] += 1
                Mall[degree][ipeq,allmonoms.index(tuple(mlist))] = c
                Mallindices[degree][ipeq,allmonoms.index(tuple(mlist))] = exportindex

            # check if any monoms are not expressed in this poly
            # if so, add another poly with the monom multiplied, will this give bad solutions?
            for igen in range(len(peq.gens)):
                if all([m[igen]==0 for m in peq.monoms()]):
                    for m,c in peq.terms():
                        mlist = list(m)
                        mlist[igen] += 1
                        degree=mlist.pop(ileftvar)
                        exportindex = degree*len(origmonoms)*len(dialyticeqs) + len(origmonoms)*ipeq+origmonoms.index(tuple(mlist))
                        assert(exportcoeffeqs[exportindex] == S.Zero)
                        exportcoeffeqs[exportindex] = c
                        Mall[degree][len(dialyticeqs)+ipeq,allmonoms.index(tuple(mlist))] = c
                        Mallindices[degree][len(dialyticeqs)+ipeq,allmonoms.index(tuple(mlist))] = exportindex
                        mlist[0] += 1
                        Mall[degree][ipeq,allmonoms.index(tuple(mlist))] = c
                        Mallindices[degree][ipeq,allmonoms.index(tuple(mlist))] = exportindex

        # have to check that the determinant is not zero for several values of ileftvar!
        # It is very common that some equations are linearly dependent and not solvable through this method.
        if self.testconsistentvalues is not None:
            linearlyindependent = False
            for itest,subs in enumerate(self.testconsistentvalues):
                if getsubs is not None:
                    # have to explicitly evaluate since testsubs can be very complex
                    subsvals = [(s,v.evalf()) for s,v in subs]
                    try:
                        subs = subsvals+getsubs(subsvals)
                    except self.CannotSolveError, e:
                        # getsubs failed (sometimes it requires solving inverse matrix), so go to next set
                        continue
                # have to sub at least twice with the global symbols
                A = Mall[maxdegree].subs(subs)
                for i in range(A.shape[0]):
                    for j in range(A.shape[1]):
                        A[i,j] = self._SubstituteGlobalSymbols(A[i,j]).subs(subs).evalf()
                eps = 10**-(self.precision-3)
                try:
                    Anumpy = numpy.array(numpy.array(A), numpy.float64)
                except ValueError, e:
                    log.warn(u'could not convert to numpy array: %s', e)
                    continue
                
                if numpy.isnan(numpy.sum(Anumpy)):
                    log.info('A has NaNs')
                    break
                eigenvals = numpy.linalg.eigvals(Anumpy)
                if all([Abs(f) > eps for f in eigenvals]):
                    try:
                        Ainv = A.inv(method='LU')
                    except ValueError, e:
                        log.error('error when taking inverse: %s', e)
                        continue
                    B = Ainv*Mall[1].subs(subs)
                    for i in range(B.shape[0]):
                        for j in range(B.shape[1]):
                            B[i,j] = self._SubstituteGlobalSymbols(B[i,j]).subs(subs).evalf()
                    C = Ainv*Mall[0].subs(subs).evalf()
                    for i in range(C.shape[0]):
                        for j in range(C.shape[1]):
                            C[i,j] = self._SubstituteGlobalSymbols(C[i,j]).subs(subs).evalf()
                    A2 = zeros((B.shape[0],B.shape[0]*2))
                    for i in range(B.shape[0]):
                        A2[i,B.shape[0]+i] = S.One
                    A2=A2.col_join((-C).row_join(-B))
                    eigenvals2,eigenvecs2 = numpy.linalg.eig(numpy.array(numpy.array(A2),numpy.float64))
                    # check if solutions can actually be extracted
                    # find all the zero eigenvalues
                    roots = []
                    numrepeating = 0
                    for ieig,eigenvalue in enumerate(eigenvals2):
                        if abs(numpy.imag(eigenvalue)) < 1e-12:
                            if abs(numpy.real(eigenvalue)) > 1:
                                ev = eigenvecs2[A.shape[0]:,ieig]
                            else:
                                ev = eigenvecs2[:A.shape[0],ieig]
                            if abs(ev[0]) < 1e-14:
                                continue
                            br = ev[1:] / ev[0]
                            dists = abs(numpy.array(roots) - numpy.real(eigenvalue))
                            if any(dists<1e-7):
                                numrepeating += 1
                            roots.append(numpy.real(eigenvalue))
                    if numrepeating > 0:
                        log.info('found %d repeating roots in solveDialytically matrix: %s',numrepeating,roots)
                        # should go on even if there's repeating roots?
                        continue
                    Atotal = None
                    for idegree in range(maxdegree+1):
                        Adegree = Mall[idegree].subs(subs)
                        for i in range(Adegree.shape[0]):
                            for j in range(Adegree.shape[1]):
                                Adegree[i,j] = self._SubstituteGlobalSymbols(Adegree[i,j]).subs(subs).evalf()
                        if Atotal is None:
                            Atotal = Adegree
                        else:
                            Atotal += Adegree*leftvar**idegree
                    # make sure the determinant of Atotal is not-zero for at least several solutions
                    leftvarvalue = leftvar.subs(subs).evalf()
                    hasnonzerodet = False
                    for testvalue in [-10*S.One, -S.One,-0.5*S.One, 0.5*S.One, S.One, 10*S.One]:
                        Atotal2 = Atotal.subs(leftvar,leftvarvalue+testvalue).evalf()
                        detvalue = Atotal2.det()
                        if abs(detvalue) > 1e-10:
                            hasnonzerodet = True
                    if not hasnonzerodet:
                        log.warn('has zero det, so failed')
                    else:
                        linearlyindependent = True
                    break
                else:
                    log.info('not all abs(eigenvalues) > %e. min is %e', eps, min([Abs(f) for f in eigenvals if Abs(f) < eps]))
            if not linearlyindependent:
                raise self.CannotSolveError('equations are not linearly independent')

        if returnmatrix:
            return Mall,allmonoms

        return exportcoeffeqs,origmonoms

    def PropagateSolvedConstants(self, AllEquations, othersolvedvars, unknownvars, \
                                 constantSymbols = None):
        """
        Sometimes equations can be like "npz" or "pp - 1", meaning npz = 0 and pp = 1. 
        
        We check these constraints and apply them to the rest of the equations.

        Returns a new set of equations.

        :param constantSymbols: the variables to try to propagage, if None will use self.pvars

        Called by AddSolution only.
        """
        constantSymbols = list(self.pvars if constantSymbols is None else constantSymbols)
            
        for othersolvedvar in othersolvedvars:
            constantSymbols.append(othersolvedvar)
            if self.IsHinge(othersolvedvar.name):
                constantSymbols.append(cos(othersolvedvar))
                constantSymbols.append(sin(othersolvedvar))
                
        newsubsdict = {}
        for eq in AllEquations:
            if not eq.has(*unknownvars) and eq.has(*constantSymbols):
                try:
                    # TGN: comment out this since it's unused
                    # reducedeq = self.SimplifyTransform(eq)
                    # eq = reduceeq ???
                    for constantSymbol in constantSymbols:
                        if eq.has(constantSymbol):
                            try:
                                peq = Poly(eq, constantSymbol)
                                if peq.degree(0) == 1:
                                    # equation is only degree 1 in the variable, and doesn't have any solvevars multiplied with it
                                    newsolution = solve(peq, constantSymbol)[0]
                                    if constantSymbol not in newsubsdict or \
                                       self.codeComplexity(newsolution) < self.codeComplexity(newsubsdict[constantSymbol]):
                                        newsubsdict[constantSymbol] = newsolution
                            except PolynomialError:
                                pass
                except PolynomialError, e:
                    # expected from simplifyTransform if eq is too complex
                    pass
                
        # first substitute everything that doesn't have othersolvedvar or unknownvars
        numberSubstitutions = []
        otherSubstitutions = []
        for var, value in newsubsdict.iteritems():
            if not value.has(*constantSymbols):
                numberSubstitutions.append((var,value))
            else:
                otherSubstitutions.append((var,value))
                
        NewEquations = []
        for ieq, eq in enumerate(AllEquations):
            if True: # not eq.has(*unknownvars):
                neweq = eq.subs(numberSubstitutions).expand()
                if neweq != S.Zero:
                    # don't expand here since otherSubstitutions could make it very complicated
                    neweq2 = neweq.subs(otherSubstitutions)
                    if self.codeComplexity(neweq2) < self.codeComplexity(neweq)*2:
                        neweq2 = neweq2.expand()
                        if self.codeComplexity(neweq2) < self.codeComplexity(neweq) and neweq2 != S.Zero:
                            NewEquations.append(neweq2)
                        else:
                            NewEquations.append(neweq)
                    else:
                        NewEquations.append(neweq)
            else:
                NewEquations.append(eq)
        return NewEquations
    
    def SolveAllEquations(self, AllEquations, \
                          curvars, othersolvedvars, \
                          solsubs, endbranchtree, \
                          currentcases    = set(), \
                          unknownvars     = [], \
                          currentcasesubs = [], \
                          canguessvars    = True):
        """
        If canguessvars is True, then we can guess variable values, prodived they satisfy required conditions

        Primitive calls are by GuessValuesAndSolveEquations, those solveFullIK_* functions.

        Recursive calls are by AddSolution.
        """

        # range of progress is [0.15, 0.45].
        # Usually scopecounters can go to several hundred
        progress = 0.45 - 0.3/(1+self._scopecounter/100)
        self._CheckPreemptFn(progress = progress)
        
        if len(curvars) == 0:
            return endbranchtree
        
        self._scopecounter += 1
        scopecounter = int(self._scopecounter)

        if currentcases is None:
            exec(ipython_str, globals(), locals())
            
        log.info('depth = %d, c = %d\n' + \
                 '        %s, %s\n' + \
                 '        cases = %s', \
                 len(currentcases), \
                 self._scopecounter, othersolvedvars, curvars, \
                 None if len(currentcases) is 0 else \
                 ("\n"+" "*16).join(str(x) for x in list(currentcases)))

        # solsubs = solsubs[:]

        # inverse substitutions
        # inv_freevarsubs = [(f[1],f[0]) for f in self.freevarsubs]
        # inv_solsubs     = [(f[1],f[0]) for f in solsubs         ]

        # single variable solutions
        solutions = []
        freevar_sol_subs = set().union(*[solsubs, self.freevarsubs])
        # equivalent to
        # self.freevarsubs + [solsub for solsub in solsubs if not solsub in self.freevarsubs]
        freevar = [f[0] for f in freevar_sol_subs]
        
        for curvar in curvars:
            othervars = unknownvars + [var for var in curvars if var != curvar]
            curvarsym = self.getVariable(curvar)
            raweqns = []
            for eq in AllEquations:

                if (len(othervars) == 0 or \
                    not eq.has(*othervars)) and \
                    eq.has(curvar, curvarsym.htvar, curvarsym.cvar, curvarsym.svar):

                    if eq.has(*freevar):
                        # neweq = eq.subs(freevar_sol_subs)
                        # log.info('\n        %s\n\n-->     %s', eq, neweq)
                        # eq = neweq
                        eq = eq.subs(freevar_sol_subs)

                    if self.CheckExpressionUnique(raweqns, eq):
                        raweqns.append(eq)

            if len(raweqns) > 0:
                try:
                    rawsolutions = self.solveSingleVariable(\
                                                            self.sortComplexity(raweqns), \
                                                            curvar, othersolvedvars, \
                                                            unknownvars = curvars + unknownvars)

                    for solution in rawsolutions:
                        self.ComputeSolutionComplexity(solution, othersolvedvars, curvars)
                        if solution.numsolutions() > 0:
                            solutions.append((solution, curvar))
                        else:
                            log.warn('solution did not have any equations')

                except self.CannotSolveError:
                    pass

        # Only return here if a solution was found that perfectly determines the unknown.
        # Otherwise, the pairwise solver could come up with something.
        #
        # There is still a problem with this: (bertold robot)
        # Sometimes an equation like atan2(y,x) evaluates to atan2(0,0) during runtime.
        # This cannot be known at compile time, so the equation is selected and any other possibilities are rejected.
        # In the bertold robot case, the next possibility is a pair-wise solution involving two variables
        #
        # TGN: don't we check Abs(y)+Abs(x) for atan2?
        
        if any([s[0].numsolutions() == 1 for s in solutions]):
            return self.AddSolution(solutions, \
                                    AllEquations, \
                                    curvars, \
                                    othersolvedvars, \
                                    solsubs, \
                                    endbranchtree, \
                                    currentcases = currentcases, \
                                    currentcasesubs = currentcasesubs, \
                                    unknownvars = unknownvars)
        
        curvarsubssol = []
        for var0, var1 in combinations(curvars,2):
            othervars = unknownvars + \
                        [var for var in curvars if var != var0 and var != var1]
            raweqns = []
            complexity = 0
            for eq in AllEquations:
                if (len(othervars) == 0 or not eq.has(*othervars)) \
                   and eq.has(var0, var1):
                    
                    eq = eq.subs(self.freevarsubs + solsubs)
                    if self.CheckExpressionUnique(raweqns, eq):
                        raweqns.append(eq)
                        complexity += self.codeComplexity(eq)
                        
            if len(raweqns) > 1:
                curvarsubssol.append((var0, var1, raweqns, complexity))
                
        curvarsubssol.sort(lambda x, y: x[3]-y[3])
        
        if len(curvars) == 2 and \
           self.IsHinge(curvars[0].name) and \
           self.IsHinge(curvars[1].name) and \
           len(curvarsubssol) > 0:
            # There are only two variables left, so two possibilities:
            #
            # EITHER two axes are aligning, OR these two variables depend on each other.
            #
            # Note that the axes' anchors also have to be along the direction!
            
            var0, var1, raweqns, complexity = curvarsubssol[0]
            dummyvar = Symbol('dummy')
            dummyvalue = var0 + var1
            NewEquations = []
            NewEquationsAll = []
            hasExtraConstraints = False
            for eq in raweqns:

                # TGN: ensure curvars is a subset of self.trigvars_subs
                assert(all([z in self.trigvars_subs for z in curvars]))

                # try dummyvar = var0 + var1
                neweq = self.trigsimp_new(eq.subs(var0, dummyvar-var1).expand(trig=True))
                eq = neweq.subs(self.freevarsubs+solsubs)
                if self.CheckExpressionUnique(NewEquationsAll, eq):
                    NewEquationsAll.append(eq)
                    
                if neweq.has(dummyvar):
                    if neweq.has(*(othervars+curvars)):
                        hasExtraConstraints = True
                        # break
                        # don't know why breaking here ...
                        # sometimes equations can be very complex but variables can still be dependent
                    else:
                        eq = neweq.subs(self.freevarsubs + solsubs)
                        if self.CheckExpressionUnique(NewEquations, eq):
                            NewEquations.append(eq)
                            
            if len(NewEquations) < 2 and hasExtraConstraints:
                # try dummyvar = var0 - var1
                NewEquations = []
                NewEquationsAll = []
                hasExtraConstraints = False
                dummyvalue = var0 - var1
                
                for eq in raweqns:
                    # TGN: ensure curvars is a subset of self.trigvars_subs
                    assert(all([z in self.trigvars_subs for z in curvars]))
                        
                    neweq = self.trigsimp_new(eq.subs(var0, dummyvar + var1).expand(trig = True))
                    eq = neweq.subs(self.freevarsubs + solsubs)
                    
                    if self.CheckExpressionUnique(NewEquationsAll, eq):
                        NewEquationsAll.append(eq)
                        
                    if neweq.has(dummyvar):
                        if neweq.has(*(othervars + curvars)):
                            hasExtraConstraints = True
                            # break
                            # don't know why breaking here ...
                            # sometimes equations can be too complex but variables can still be dependent
                        else:
                            eq = neweq.subs(self.freevarsubs + solsubs)
                            if self.CheckExpressionUnique(NewEquations, eq):
                                NewEquations.append(eq)
                                
            if len(NewEquations) >= 2:
                dummysolutions = []
                try:
                    rawsolutions = self.solveSingleVariable(NewEquations, dummyvar, othersolvedvars, \
                                                            unknownvars = curvars+unknownvars)
                    for solution in rawsolutions:
                        self.ComputeSolutionComplexity(solution, othersolvedvars, curvars)
                        dummysolutions.append(solution)
                        
                except self.CannotSolveError:
                    pass
                
                if any([s.numsolutions()==1 for s in dummysolutions]):
                    # two axes are aligning, so modify the solutions to reflect the original variables and add a free variable
                    log.info('found two aligning axes %s: %r', dummyvalue, NewEquations)
                    solutions = []
                    for dummysolution in dummysolutions:
                        if dummysolution.numsolutions() != 1:
                            continue
                        if dummysolution.jointevalsin is not None or \
                           dummysolution.jointevalcos is not None:
                            log.warn('dummy solution should not have sin/cos parts!')

                        sindummyvarsols = []
                        cosdummyvarsols = []
                        for eq in NewEquations:
                            sols = solve(eq, sin(dummyvar))
                            sindummyvarsols += sols
                            sols = solve(eq, cos(dummyvar))
                            cosdummyvarsols += sols
                        
                        # double check with NewEquationsAll that everything evaluates to 0
                        newsubs = [( value,  sin(dummyvar)) for value in sindummyvarsols] + \
                                  [( value,  cos(dummyvar)) for value in cosdummyvarsols] + \
                                  [(-value, -sin(dummyvar)) for value in sindummyvarsols] + \
                                  [(-value, -cos(dummyvar)) for value in cosdummyvarsols]
                        allzeros = True
                        for eq in NewEquationsAll:
                            if trigsimp(eq.subs(newsubs)) != S.Zero:
                                allzeros = False
                                break
                            
                        if allzeros:
                            solution = AST.SolverSolution(curvars[0].name, \
                                                          isHinge = self.IsHinge(curvars[0].name))
                            solution.jointeval = [dummysolution.jointeval[0] - dummyvalue + curvars[0]]
                            self.ComputeSolutionComplexity(solution, othersolvedvars, curvars)
                            solutions.append((solution, curvars[0]))
                        else:
                            log.warn('Not all equations evaluate to zero, so variables %s are not collinear', curvars)
                            
                    if len(solutions) > 0:
                        tree = self.AddSolution(solutions, raweqns, curvars[0:1], \
                                                othersolvedvars + curvars[1:2], \
                                                solsubs + self.getVariable(curvars[1]).subs, \
                                                endbranchtree, \
                                                currentcases = currentcases, \
                                                currentcasesubs = currentcasesubs,
                                                unknownvars = unknownvars)
                        if tree is not None:
                            return [AST.SolverFreeParameter(curvars[1].name, tree)]
                else:
                    log.warn('Almost found two axes, but number of solutions is %r', \
                             [s.numsolutions() == 1 for s in dummysolutions])
                    
        for var0, var1, raweqns, complexity in curvarsubssol:
            try:
                rawsolutions = self.SolvePrismaticHingePairVariables(raweqns, var0, var1, \
                                                                     othersolvedvars, \
                                                                     unknownvars = curvars + unknownvars)
                for solution in rawsolutions:
                    # solution.subs(inv_freevarsubs)
                    self.ComputeSolutionComplexity(solution, othersolvedvars, curvars)
                    solutions.append((solution, Symbol(solution.jointname)))
                    
                if len(rawsolutions) > 0: # solving a pair is rare, so any solution will do
                    # TGN: so we don't try others in the for-loop?
                    break
            except self.CannotSolveError:
                pass
            
        for var0, var1, raweqns, complexity in curvarsubssol:
            try:
                rawsolutions = self.SolvePairVariables(raweqns, var0, var1, \
                                                       othersolvedvars, \
                                                       unknownvars = curvars + unknownvars)
            except self.CannotSolveError, e:
                log.debug(e)
#                 try:
#                     rawsolutions=self.SolvePrismaticHingePairVariables(raweqns,var0,var1,othersolvedvars,unknownvars=curvars+unknownvars)
#                 except self.CannotSolveError, e:
#                     log.debug(e)
                rawsolutions = []
            for solution in rawsolutions:
                #solution.subs(inv_freevarsubs)
                try:
                    self.ComputeSolutionComplexity(solution, othersolvedvars, curvars)
                    solutions.append((solution, Symbol(solution.jointname)))
                except self.CannotSolveError, e:
                    log.warn(u'equation failed to compute solution complexity: %s', solution.jointeval)
            if len(rawsolutions) > 0: # solving a pair is rare, so any solution will do
                # TGN: so we don't try others in the for-loop?
                break
                        
        # take the least complex solution and go on
        if len(solutions) > 0:
            return self.AddSolution(solutions, AllEquations, \
                                    curvars, othersolvedvars, \
                                    solsubs, \
                                    endbranchtree, \
                                    currentcases = currentcases, \
                                    currentcasesubs = currentcasesubs, \
                                    unknownvars = unknownvars)
        
        # test with higher degrees, necessary?
        for curvar in curvars:
            othervars = unknownvars + [var for var in curvars if var != curvar]
            raweqns = []
            for eq in AllEquations:
                if (len(othervars) == 0 or not eq.has(*othervars)) and eq.has(curvar):
                    eq = eq.subs(self.freevarsubs + solsubs)
                    if self.CheckExpressionUnique(raweqns, eq):
                        raweqns.append(eq)
                        
            for raweqn in raweqns:
                try:
                    log.debug('testing with higher degrees')
                    solution = self.solveHighDegreeEquationsHalfAngle([raweqn], self.getVariable(curvar))
                    self.ComputeSolutionComplexity(solution, othersolvedvars, curvars)
                    solutions.append((solution, curvar))
                    
                except self.CannotSolveError:
                    pass
               
        if len(solutions) > 0:
            return self.AddSolution(solutions, AllEquations, \
                                    curvars, othersolvedvars, \
                                    solsubs, \
                                    endbranchtree, \
                                    currentcases = currentcases, \
                                    currentcasesubs = currentcasesubs, \
                                    unknownvars = unknownvars)
        
        # solve with all 3 variables together?
#         htvars = [self.getVariable(varsym).htvar for varsym in curvars]
#         reducedeqs = []
#         for eq in AllEquations:
#             if eq.has(*curvars):
#                 num, denom, htvarsubsinv = self.ConvertSinCosEquationToHalfTan(eq, curvars)
#                 reducedeqs.append(Poly(num, *htvars))
# 

        # only guess if final joint to be solved, or there exists current cases and at least one joint has been solved already.
        # don't want to start guessing when no joints have been solved yet, this indicates bad equations
        if canguessvars and \
           len(othersolvedvars) + len(curvars) == len(self.freejointvars) + len(self._solvejointvars) and \
           (len(curvars) == 1 or (len(curvars) < len(self._solvejointvars) and \
                                  currentcases is not None and \
                                  len(currentcases) > 0)):
            # only estimate when deep in the hierarchy, do not want the guess to be executed all the time
            # perhaps there's a degree of freedom that is not trivial to compute?
            # take the highest hinge variable and set it
            log.info('trying to guess variable from %r', curvars)
            return self.GuessValuesAndSolveEquations(AllEquations, \
                                                     curvars, othersolvedvars, solsubs, \
                                                     endbranchtree, \
                                                     currentcases, \
                                                     unknownvars, \
                                                     currentcasesubs)
        
        # have got this far, so perhaps two axes are aligned?
        #
        # TGN: so we cannot detect such aligning before?
        #
        raise self.CannotSolveError('SolveAllEquations failed to find a variable to solve')
        
    def AddSolution(self, solutions, AllEquations, \
                    curvars, othersolvedvars, \
                    solsubs, endbranchtree, \
                    currentcases    = set(), \
                    currentcasesubs = [], \
                    unknownvars     = []):
        """
        Take the least complex solution of a set of solutions and resume solving.

        Called by SolveAllEquations and also calls it recursively.
        """
        
        self._CheckPreemptFn()
        self._scopecounter += 1
        scopecounter = int(self._scopecounter)
        # remove solution that has infinite scores
        solutions = [s for s in solutions if s[0].score < oo and \
                     s[0].checkValidSolution()] 
        if len(solutions) == 0:
            raise self.CannotSolveError('No valid solutions')
            
        solutions.sort(lambda x, y: x[0].score-y[0].score)
        hasonesolution = False
        
        for solution in solutions:
            checkforzeros = solution[0].checkforzeros
            hasonesolution |= solution[0].numsolutions() == 1
            if len(checkforzeros) == 0 and solution[0].numsolutions() == 1:
                # did find a good solution, so take it. Make sure to check any zero branches
                var = solution[1]
                newvars = curvars[:]
                newvars.remove(var)
                return [solution[0].subs(solsubs)] + \
                    self.SolveAllEquations(AllEquations, \
                                           curvars = newvars, \
                                           othersolvedvars = othersolvedvars+[var], \
                                           solsubs = solsubs+self.getVariable(var).subs, \
                                           endbranchtree = endbranchtree, \
                                           currentcases = currentcases, \
                                           currentcasesubs = currentcasesubs, \
                                           unknownvars = unknownvars)
            
        if not hasonesolution:
            # check again except without the number of solutions requirement
            for solution in solutions:
                checkforzeros = solution[0].checkforzeros
                if len(checkforzeros) == 0:
                    # did find a good solution, so take it. Make sure to check any zero branches
                    var = solution[1]
                    newvars = curvars[:]
                    newvars.remove(var)
                    return [solution[0].subs(solsubs)] + \
                        self.SolveAllEquations(AllEquations, \
                                               curvars = newvars, \
                                               othersolvedvars = othersolvedvars+[var], \
                                               solsubs = solsubs+self.getVariable(var).subs, \
                                               endbranchtree = endbranchtree, \
                                               currentcases = currentcases, \
                                               currentcasesubs = currentcasesubs, \
                                               unknownvars = unknownvars)

        originalGlobalSymbols = self.globalsymbols
        # all solutions have check for zero equations
        # choose the variable with the shortest solution and compute (this is a conservative approach)
        usedsolutions = []
        # remove any solutions with similar checkforzero constraints (because they are essentially the same)
        for solution, var in solutions:
            solution.subs(solsubs)
            if len(usedsolutions) == 0:
                usedsolutions.append((solution, var))
            else:
                match = False
                for usedsolution, usedvar in usedsolutions:
                    if len(solution.checkforzeros) == len(usedsolution.checkforzeros):
                        if not any([self.CheckExpressionUnique(usedsolution.checkforzeros, eq) \
                                    for eq in solution.checkforzeros]):
                            match = True
                            break
                if not match:
                    usedsolutions.append((solution, var))
                    if len(usedsolutions) >= 3:
                        # don't need more than three alternatives (used to be two, but then lookat barrettwam4 proved that wrong)
                        break

        allvars            = list(chain.from_iterable([self.getVariable(v).vars for v in curvars]))
        allothersolvedvars = list(chain.from_iterable([self.getVariable(v).vars for v in othersolvedvars]))
        
        prevbranch = lastbranch = []
        nextsolutions = dict()
            
        if self.degeneratecases is None:
            self.degeneratecases = self.DegenerateCases()
        handledconds = self.degeneratecases.GetHandledConditions(currentcases)
        
        # one to one correspondence with usedsolutions and the SolverCheckZeros hierarchies
        # (used for cross product of equations later on)
        
        zerosubstitutioneqs = [] # indexed by reverse ordering of usedsolutions (len(usedsolutions)-solutionindex-1)
        # zerosubstitutioneqs equations flattened for easier checking
        flatzerosubstitutioneqs = []
        hascheckzeros = False
        
        addhandleddegeneratecases = [] # for bookkeeping/debugging
        
        # iterate in reverse order and put the most recently processed solution at the front.
        # There is a problem with this algorithm transferring the degenerate cases correctly.
        # Although the zeros of the first equation are checked, they are not added as conditions to the later equations,
        # so that the later equations will also use variables as unknowns
        # (even though they are determined to be specific constants). This is most apparent in rotations.
        
        for solution, var in usedsolutions[::-1]:
            # there are divide by zeros, so check if they can be explicitly solved for joint variables
            checkforzeros = []
            localsubstitutioneqs = []
            for checkzero in solution.checkforzeros:
                if checkzero.has(*allvars):
                    log.info('Ignore special check for zero since it has symbols %s: %s', \
                             str(allvars), str(checkzero))
                    continue
                
                # Don't bother trying to extract something if too complex
                # (takes a lot of time to check and most likely nothing will be extracted).
                # 120 is from heuristics
                checkzeroComplexity = self.codeComplexity(checkzero)
                if checkzeroComplexity > 120: 
                    log.warn('Checkforzero too big (%d): %s', checkzeroComplexity, checkzero)
                    # don't even add it if it is too big
                    if checkzeroComplexity < 500:
                        checkforzeros.append(checkzero)
                        #self.removecommonexprs(checkzero.evalf())
                else:
                    checkzero2 = self._SubstituteGlobalSymbols(checkzero, originalGlobalSymbols)
                    checkzero2Complexity = self.codeComplexity(checkzero2)
                    if checkzero2Complexity < 2*checkzeroComplexity: # check that with substitutions, things don't get too big
                        checkzero = checkzero2
                        # fractions could get big, so evaluate directly
                        checkzeroeval = checkzero.evalf()
                        if checkzero2Complexity < self.codeComplexity(checkzeroeval):
                            checkforzeros.append(checkzero)
                        else:
                            checkforzeros.append(checkzero.evalf())
                            #self.removecommonexprs(checkzero.evalf())
                            
                    checksimplezeroexprs = [checkzero]
                    if not checkzero.has(*allothersolvedvars):
                        
                        sumsquaresexprs = self._GetSumSquares(checkzero)
                        
                        if sumsquaresexprs is not None:
                            checksimplezeroexprs += sumsquaresexprs
                            sumsquaresexprstozero = []
                            #sumsquaresexprstozero = [sumsquaresexpr \
                            #                         for sumsquaresexpr in sumsquaresexprs \
                            #                         if sumsquaresexpr.is_Symbol]
                            
                            for sumsquaresexpr in sumsquaresexprs:
                                if sumsquaresexpr.is_Symbol:
                                    sumsquaresexprstozero.append(sumsquaresexpr)
                                elif sumsquaresexpr.is_Mul:
                                    for arg in sumsquaresexpr.args:
                                        if arg.is_Symbol:
                                            sumsquaresexprstozero.append(arg)
                            
                            if len(sumsquaresexprstozero) > 0:
                               log.info('%r', [sumsquaresexprstozero,checkzero, \
                                               [(sumsquaresexpr,S.Zero) \
                                                for sumsquaresexpr in sumsquaresexprstozero], []])

                               toappend = [sumsquaresexprstozero, \
                                           checkzero, \
                                           [(sumsquaresexpr, S.Zero) \
                                            for sumsquaresexpr in sumsquaresexprstozero], []]
                               log.info('%r', toappend)
                               localsubstitutioneqs.append(toappend)
                               handledconds += sumsquaresexprstozero
                                
                    for checksimplezeroexpr in checksimplezeroexprs:
                        #if checksimplezeroexpr.has(*othersolvedvars): # cannot do this check since sjX,cjX might be used
                        for othervar in othersolvedvars:
                            sothervar = self.getVariable(othervar).svar
                            cothervar = self.getVariable(othervar).cvar
                            if checksimplezeroexpr.has(othervar,sothervar,cothervar):
                                # the easiest thing to check first is if the equation evaluates to zero on boundaries
                                # 0, pi/2, pi, -pi/2
                                s = AST.SolverSolution(othervar.name, \
                                                       jointeval = [], \
                                                       isHinge = self.IsHinge(othervar.name))
                                for value in [S.Zero, pi/2, pi, -pi/2]:
                                    try:
                                        # doing (1/x).subs(x,0) produces a RuntimeError (infinite recursion...)
                                        # TGN: may not need that many digits (used to be n=30)
                                        checkzerosub = checksimplezeroexpr.subs([(othervar,  value), \
                                                                                 (sothervar, sin(value).evalf(n=2)), \
                                                                                 (cothervar, cos(value).evalf(n=2))])
                                        
                                        if self.isValidSolution(checkzerosub) and \
                                           checkzerosub.evalf(n=30) == S.Zero:
                                            if s.jointeval is None:
                                                s.jointeval = []
                                            s.jointeval.append(S.One*value)
                                    except (RuntimeError, AssertionError),e: # 
                                        log.warn('othervar %s = %f: %s', str(othervar), value, e)

                                
                                ss = [s] if s.jointeval is not None and len(s.jointeval) > 0 else []
                                    
                                try:
                                    # checksimplezeroexpr can be simple like -cj4*r21 - r20*sj4
                                    # in which case the solutions would be
                                    # [-atan2(-r21, -r20), -atan2(-r21, -r20) + 3.14159265358979]
                                    ss += self.solveSingleVariable([checksimplezeroexpr.subs([(sothervar, sin(othervar)), \
                                                                                              (cothervar, cos(othervar))])], \
                                                                   othervar, \
                                                                   othersolvedvars)
                                except PolynomialError:
                                    # checksimplezeroexpr was too complex
                                    pass
                                
                                except self.CannotSolveError,e:
                                    # this is actually a little tricky
                                    # sometimes really good solutions can have a divide that looks like:
                                    # ((0.405 + 0.331*cj2)**2 + 0.109561*sj2**2 (manusarm_left)
                                    # This will never be 0, but the solution cannot be solved.
                                    # Instead of rejecting, add a condition to check if checksimplezeroexpr itself is 0 or not
                                    pass
                                
                                for s in ss:
                                    # can actually simplify Positions and possibly get a new solution!
                                    if s.jointeval is not None:
                                        for eq in s.jointeval:
                                            eq = self._SubstituteGlobalSymbols(eq, originalGlobalSymbols)
                                            # why checking for just number?
                                            # ok to check if solution doesn't contain any other variables?
                                            # if the equation is non-numerical, make sure it isn't deep in the degenerate cases
                                            if eq.is_number \
                                               or (len(currentcases) <= 1 and \
                                                   not eq.has(*allothersolvedvars) and \
                                                   self.codeComplexity(eq) < 100):
                                                
                                                isimaginary = self.AreAllImaginaryByEval(eq) \
                                                              or eq.evalf().has(I)
                                                
                                                # TODO should use the fact that eq is imaginary
                                                if isimaginary:
                                                    log.warn('eq %s is imaginary, but currently do not support this', eq)
                                                    continue
                                                
                                                dictequations = []
                                                if not eq.is_number and not eq.has(*allothersolvedvars):
                                                    # not dependent on variables
                                                    # so it could be in the form of atan(px,py),
                                                    # so we convert to a global symbol since it never changes
                                                    sym = self.gsymbolgen.next()
                                                    dictequations.append((sym, eq))
                                                    #eq = sym
                                                    sineq = self.gsymbolgen.next()
                                                    dictequations.append((sineq,self.SimplifyAtan2(sin(eq))))
                                                    coseq = self.gsymbolgen.next()
                                                    dictequations.append((coseq,self.SimplifyAtan2(cos(eq))))
                                                else:
                                                    sineq = sin(eq).evalf(n=30)
                                                    coseq = cos(eq).evalf(n=30)
                                                    
                                                cond = Abs(othervar-eq.evalf(n=30))
                                                
                                                if self.CheckExpressionUnique(handledconds, cond):
                                                    evalcond = fmod(cond+pi,2*pi)-pi if self.IsHinge(othervar.name) else cond
                                                    toappend = [[cond], \
                                                                evalcond, \
                                                                [(sothervar, sineq), \
                                                                 (sin(othervar), sineq), \
                                                                 (cothervar,coseq), \
                                                                 (cos(othervar),coseq), \
                                                                 (othervar,eq)], \
                                                                dictequations]
                                                    log.info('%r', toappend)
                                                    # exec(ipython_str,globals(),locals())
                                                    localsubstitutioneqs.append(toappend)
                                                    handledconds.append(cond)
                                                    
                                    elif s.jointevalsin is not None:
                                        for eq in s.jointevalsin:
                                            eq = self.SimplifyAtan2(self._SubstituteGlobalSymbols(eq, originalGlobalSymbols))
                                            if eq.is_number or (len(currentcases) <= 1 and \
                                                                not eq.has(*allothersolvedvars) and \
                                                                self.codeComplexity(eq) < 100):
                                                dictequations = []
                                                # test when cos(othervar) > 0
                                                # DO NOT use asin(eq)
                                                # because eq = (-pz**2/py**2)**(1/2) would produce imaginary numbers
                                                #
                                                # cond = othervar-asin(eq).evalf(n=30)
                                                # test if eq is imaginary
                                                # If so, then only solution is when sothervar==0 and eq==0
                                                isimaginary = self.AreAllImaginaryByEval(eq) \
                                                              or eq.evalf().has(I)
                                                if isimaginary:
                                                    cond = abs(sothervar) + abs((eq**2).evalf(n=30)) + abs(sign(cothervar)-1)
                                                else:
                                                    if not eq.is_number and not eq.has(*allothersolvedvars):
                                                        # not dependent on variables
                                                        # so it could be in the form of atan(px,py)
                                                        # so we convert to a global symbol since it never changes
                                                        sym = self.gsymbolgen.next()
                                                        dictequations.append((sym,eq))
                                                        #eq = sym
                                                    cond = abs(sothervar-eq.evalf(n=30)) + abs(sign(cothervar)-1)
                                                    
                                                if self.CheckExpressionUnique(handledconds, cond):
                                                    evalcond = fmod(cond+pi, 2*pi) - pi if self.IsHinge(othervar.name) else cond
                                                        
                                                    if isimaginary:
                                                        toappend = [[cond], \
                                                                    evalcond, \
                                                                    [(sothervar,     S.Zero), \
                                                                     (sin(othervar), S.Zero), \
                                                                     (cothervar,     S.One),  \
                                                                     (cos(othervar), S.One),  \
                                                                     (othervar,      S.One)], \
                                                                    dictequations]
                                                    else:
                                                        # exec(ipython_str, globals(), locals())
                                                        toappend = [[cond], \
                                                                    evalcond, \
                                                                    [(sothervar,     eq), \
                                                                     (sin(othervar), eq), \
                                                                     (cothervar,    sqrt(1-eq*eq).evalf(n=30)), \
                                                                     (cos(othervar),sqrt(1-eq*eq).evalf(n=30)), \
                                                                     (othervar, asin(eq).evalf(n=30))], \
                                                                    dictequations]
                                                        
                                                    localsubstitutioneqs.append(toappend)
                                                    handledconds.append(cond)
                                                    
                                                # test when cos(othervar) < 0
                                                cond = abs(sign(cothervar)+1) + (abs(sothervar) + abs((eq**2).evalf(n=30)) \
                                                                                 if isimaginary else \
                                                                                 abs(sothervar-eq.evalf(n=30)))
                                                    
                                                #cond=othervar-(pi-asin(eq).evalf(n=30))
                                                if self.CheckExpressionUnique(handledconds, cond):
                                                    
                                                    evalcond = fmod(cond+pi,2*pi)-pi if self.IsHinge(othervar.name) else cond
                                                        
                                                    if isimaginary:
                                                        toappend = [[cond], \
                                                                    evalcond, \
                                                                    [(sothervar,      S.Zero), \
                                                                     (sin(othervar),  S.Zero), \
                                                                     (cothervar,     -S.One),  \
                                                                     (cos(othervar), -S.One),  \
                                                                     (othervar, pi.evalf(n=30))], \
                                                                    dictequations]

                                                    else:
                                                        toappend = [[cond], \
                                                                    evalcond, \
                                                                    [(sothervar,     eq), \
                                                                     (sin(othervar), eq), \
                                                                     (cothervar,    -sqrt(1-eq*eq).evalf(n=30)), \
                                                                     (cos(othervar),-sqrt(1-eq*eq).evalf(n=30)), \
                                                                     (othervar, \
                                                                      (pi-asin(eq)).evalf(n=30))], \
                                                                    dictequations]
                                                        
                                                    localsubstitutioneqs.append(toappend)
                                                    handledconds.append(cond)
                                                    
                                    elif s.jointevalcos is not None:
                                        for eq in s.jointevalcos:
                                            eq = self.SimplifyAtan2(self._SubstituteGlobalSymbols(eq, originalGlobalSymbols))
                                            if eq.is_number or (len(currentcases) <= 1 and \
                                                                not eq.has(*allothersolvedvars) and \
                                                                self.codeComplexity(eq) < 100):
                                                
                                                dictequations = []
                                                # test when sin(othervar) > 0
                                                # DO NOT use acos(eq) because eq = (-pz**2/px**2)**(1/2),
                                                # which would produce imaginary numbers
                                                # that's why check eq.evalf().has(I)
                                                # cond = othervar-acos(eq).evalf(n=30)
                                                isimaginary = self.AreAllImaginaryByEval(eq) or \
                                                              eq.evalf().has(I)
                                                
                                                if isimaginary:
                                                    cond = abs(cothervar) + abs((eq**2).evalf(n=30)) + abs(sign(sothervar)-1)
                                                else:
                                                    if not (eq.is_number or eq.has(*allothersolvedvars)):
                                                        # not dependent on variables
                                                        # so it could be in the form of atan(px,py),
                                                        # so convert to a global symbol since it never changes
                                                        sym = self.gsymbolgen.next()
                                                        dictequations.append((sym, eq))
                                                        eq = sym
                                                    cond = abs(cothervar-eq.evalf(n=30)) + abs(sign(sothervar)-1)
                                                    
                                                if self.CheckExpressionUnique(handledconds, cond):
                                                    
                                                    evalcond = fmod(cond+pi,2*pi)-pi if self.IsHinge(othervar.name) else cond
                                                        
                                                    if isimaginary:
                                                        toappend = [[cond], \
                                                                    evalcond, \
                                                                    [(sothervar,     S.One),  \
                                                                     (sin(othervar), S.One),  \
                                                                     (cothervar,     S.Zero), \
                                                                     (cos(othervar), S.Zero), \
                                                                     (othervar, (pi/2).evalf(n=30))], \
                                                                    dictequations]
                                                    else:
                                                        toappend = [[cond], \
                                                                    evalcond, \
                                                                    [(sothervar,     sqrt(1-eq*eq).evalf(n=30)), \
                                                                     (sin(othervar), sqrt(1-eq*eq).evalf(n=30)), \
                                                                     (cothervar,     eq), \
                                                                     (cos(othervar), eq), \
                                                                     (othervar, acos(eq).evalf(n=30))], \
                                                                    dictequations]
                                                        
                                                    log.info('%r', toappend)
                                                    # exec(ipython_str,globals(),locals())
                                                    localsubstitutioneqs.append(toappend)
                                                    handledconds.append(cond)
                                                    
                                                #cond=othervar+acos(eq).evalf(n=30)
                                                cond = abs(sign(sothervar)+1) + (abs(cothervar) + abs((eq**2).evalf(n=30)) \
                                                                                 if isimaginary else \
                                                                                 abs(cothervar-eq.evalf(n=30)))

                                                if self.CheckExpressionUnique(handledconds, cond):

                                                    evalcond = fmod(cond+pi,2*pi)-pi if self.IsHinge(othervar.name) else cond
                                                        
                                                    if isimaginary:
                                                        toappend = [[cond], \
                                                                    evalcond, \
                                                                    [(sothervar,    -S.One),  \
                                                                     (sin(othervar),-S.One),  \
                                                                     (cothervar,     S.Zero), \
                                                                     (cos(othervar), S.Zero), \
                                                                     (othervar, (-pi/2).evalf(n=30))], \
                                                                    dictequations]
                                                    else:
                                                        toappend = [[cond], \
                                                                    evalcond, \
                                                                    [(sothervar,     -sqrt(1-eq*eq).evalf(n=30)), \
                                                                     (sin(othervar), -sqrt(1-eq*eq).evalf(n=30)), \
                                                                     (cothervar,     eq), \
                                                                     (cos(othervar), eq), \
                                                                     (othervar, -acos(eq).evalf(n=30))], \
                                                                    dictequations]
                                                        
                                                    log.info('%r', toappend)
                                                    localsubstitutioneqs.append(toappend)
                                                    handledconds.append(cond)

            if len(localsubstitutioneqs)>0:
                log.info('%r', localsubstitutioneqs)
                #exec(ipython_str, globals(), locals())
            flatzerosubstitutioneqs += localsubstitutioneqs
            zerosubstitutioneqs.append(localsubstitutioneqs)
            
            if not var in nextsolutions:
                try:
                    newvars = curvars[:]
                    newvars.remove(var)
                    # degenreate cases should get restored here since once we go down a particular branch, there's no turning back
                    olddegeneratecases = self.degeneratecases
                    self.degeneratecases = olddegeneratecases.Clone()
                    nextsolutions[var] = self.SolveAllEquations(AllEquations, \
                                                                curvars = newvars, \
                                                                othersolvedvars = othersolvedvars + [var], \
                                                                solsubs = solsubs + self.getVariable(var).subs, \
                                                                endbranchtree = endbranchtree, \
                                                                currentcases = currentcases, \
                                                                currentcasesubs = currentcasesubs, \
                                                                unknownvars = unknownvars)
                finally:
                    addhandleddegeneratecases += olddegeneratecases.handleddegeneratecases
                    self.degeneratecases = olddegeneratecases
                    
            if len(checkforzeros) > 0:
                hascheckzeros = True                
                solvercheckzeros = AST.SolverCheckZeros(jointname = var.name, \
                                                        jointcheckeqs = checkforzeros, \
                                                        nonzerobranch = [solution]+nextsolutions[var], \
                                                        zerobranch = prevbranch,anycondition = True, \
                                                        thresh = solution.GetZeroThreshold())
                # have to transfer the dictionary!
                solvercheckzeros.dictequations = originalGlobalSymbols.items() + solution.dictequations                    
                solvercheckzeros.equationsused = AllEquations
                solution.dictequations = []
                prevbranch=[solvercheckzeros]
            else:
                prevbranch = [solution] + nextsolutions[var]
        
        if len(prevbranch) == 0:
            raise self.CannotSolveError('failed to add solution!')

        maxlevel2scopecounter = 300 # used to limit how deep the hierarchy goes or otherwise IK can get too big
        if len(currentcases) >= self.maxcasedepth or \
           (scopecounter > maxlevel2scopecounter and \
            len(currentcases) >= 2):
            
            log.warn('c = %d, %d levels deep in checking degenerate cases, skip\n' + \
                     '        curvars = %r\n' + \
                     '        AllEquations = %s', \
                     scopecounter, len(currentcases), curvars, \
                     ("\n"+" "*23).join(str(x) for x in list(AllEquations)))

            varlist =  [(var, eq if eq.is_Symbol or eq.is_number else \
                         self.SimplifyAtan2(self._SubstituteGlobalSymbols(eq, originalGlobalSymbols))) \
                        for var, eq in currentcasesubs]
            
            lastbranch.append(AST.SolverBreak('%d cases reached' % self.maxcasedepth, \
                                              varlist, \
                                              othersolvedvars, \
                                              solsubs, \
                                              endbranchtree))

            return prevbranch
        
        # fill the last branch with all the zero conditions
        if hascheckzeros:
            # count the number of rotation symbols seen in the current cases
            numRotSymbolsInCases = 0
            if self._iktype == 'transform6d' or \
               self._iktype == 'rotation3d':
                rotsymbols = set(self.Tee[:3,:3]).union([Symbol('new_r00'), Symbol('new_r01'), Symbol('new_r02'), \
                                                         Symbol('new_r10'), Symbol('new_r11'), Symbol('new_r12'), \
                                                         Symbol('new_r20'), Symbol('new_r21'), Symbol('new_r22')])
                for var, eq in currentcasesubs:
                    if var in rotsymbols:
                        numRotSymbolsInCases += 1
            else:
                rotsymbols = []
                
            # if not equations found, try setting two variables at once
            # also try setting px, py, or pz to 0 (barrettwam4 lookat)
            # sometimes can get the following: cj3**2*sj4**2 + cj4**2
            threshnumsolutions = 1 # number of solutions to take usedsolutions[:threshnumsolutions] for the dual values
            for isolution,(solution,var) in enumerate(usedsolutions[::-1]):
                if isolution < len(usedsolutions) - threshnumsolutions and \
                   len(flatzerosubstitutioneqs) > 0:
                    # have at least one zero condition...
                    continue
                localsubstitutioneqs = []
                for checkzero in solution.checkforzeros:
                    if checkzero.has(*allvars):
                        log.info('ignoring special check for zero 2 since it has symbols %s: %s', \
                                 str(allvars), str(checkzero))
                        continue
                    
                    # don't bother trying to extract something if too complex
                    # (takes a lot of time to check and most likely nothing will be extracted)
                    # 120 is from heuristics
                    if self.codeComplexity(checkzero) > 120:
                        continue
                    
                    possiblesubs = []
                    ishinge = []
                    for preal in self.Tee[:3,3]:
                        if checkzero.has(preal):
                            possiblesubs.append([(preal,S.Zero)])
                            ishinge.append(False)
                            
                    # have to be careful with the rotations since they are dependent on each other.
                    # E.g. If r00 = r01 = 0, then r02 = 1 or -1 and r12 = r22 = 0.
                    # Then r10, r12, r20, r21 is a 2D rotation matrix
                    
                    if numRotSymbolsInCases < 2:
                        for preal in rotsymbols:
                            if checkzero.has(preal):
                                possiblesubs.append([(preal,S.Zero)])
                                ishinge.append(False)
                                
                    for othervar in othersolvedvars:
                        othervarobj = self.getVariable(othervar)
                        if checkzero.has(*othervarobj.vars):
                            if not self.IsHinge(othervar.name):
                                possiblesubs.append([(othervar,S.Zero)])
                                ishinge.append(False)
                                continue
                            else:
                                sothervar = othervarobj.svar
                                cothervar = othervarobj.cvar
                                for value in [S.Zero, pi/2, pi, -pi/2]:
                                    # TGN: may not need that many digits (used to be n=30)
                                    possiblesubs.append([(othervar,      value), \
                                                         (sothervar,     sin(value).evalf(n=2)), \
                                                         (sin(othervar), sin(value).evalf(n=2)), \
                                                         (cothervar,     cos(value).evalf(n=2)), \
                                                         (cos(othervar), cos(value).evalf(n=2))])
                                    ishinge.append(True)
                                    
                    # all possiblesubs are present in checkzero
                    for ipossiblesub, possiblesub in enumerate(possiblesubs):
                        try:
                            eq = checkzero.subs(possiblesub).evalf(n=30)
                        except RuntimeError, e:
                            # most likely doing (1/x).subs(x,0) produces a RuntimeError (infinite recursion...)
                            log.warn(e)
                            continue
                        
                        if not self.isValidSolution(eq):
                            continue
                        
                        # only take the first index
                        possiblevar, possiblevalue = possiblesub[0]
                        cond = Abs(possiblevar - possiblevalue.evalf(n=30))
                        if not self.CheckExpressionUnique(handledconds, cond):
                            # already present, so don't use it for double expressions
                            continue
                        
                        evalcond = Abs(fmod(possiblevar-possiblevalue+pi, 2*pi)-pi) \
                                   if ishinge[ipossiblesub] else cond
                            
                        if eq == S.Zero:
                            log.info('c = %d, adding case %s = %s in %s', \
                                     scopecounter, possiblevar, possiblevalue, checkzero)
                            
                            # if the variable is 1 and part of the rotation matrix, can deduce other variables
                            if possiblevar in rotsymbols and (possiblevalue == S.One or possiblevalue == -S.One):
                                row1 = int(possiblevar.name[-2])
                                col1 = int(possiblevar.name[-1])
                                possiblevarname = possiblevar.name[:-2]
                                possiblesub.append((Symbol('%s%d%d'%(possiblevarname, row1, (col1+1)%3)), S.Zero))
                                possiblesub.append((Symbol('%s%d%d'%(possiblevarname, row1, (col1+2)%3)), S.Zero))
                                possiblesub.append((Symbol('%s%d%d'%(possiblevarname, (row1+1)%3, col1)), S.Zero))
                                possiblesub.append((Symbol('%s%d%d'%(possiblevarname, (row1+2)%3, col1)), S.Zero))
                                
                            checkexpr = [[cond], evalcond, possiblesub, []]
                            log.info('%r', flatzerosubstitutioneqs)
                            log.info('%r', checkexpr)
                            # exec(ipython_str, globals(), locals())
                            flatzerosubstitutioneqs.append(checkexpr)
                            localsubstitutioneqs.append(checkexpr)
                            handledconds.append(cond)
                            continue
                        
                        # try another possiblesub
                        for ipossiblesub2, possiblesub2 in enumerate(possiblesubs[ipossiblesub+1:]):
                            possiblevar2, possiblevalue2 = possiblesub2[0]
                            if possiblevar == possiblevar2:
                                # same var, so skip
                                continue
                            try:
                                eq2 = eq.subs(possiblesub2).evalf(n=30)
                            except RuntimeError, e:
                                # most likely doing (1/x).subs(x,0) produces a RuntimeError (infinite recursion...)
                                log.warn(e)
                                continue
                            
                            if not self.isValidSolution(eq2):
                                continue
                            
                            if eq2 == S.Zero:
                                possiblevar2, possiblevalue2 = possiblesub2[0]
                                cond2 = Abs(possiblevar2-possiblevalue2.evalf(n=30))
                                if not self.CheckExpressionUnique(handledconds, cond2):
                                    # already present, so don't use it for double expressions
                                    continue
                                
                                # don't combine the conditions like cond+cond2
                                # instead test them individually so we can reduce the solution tree
                                
                                evalcond2 = Abs(fmod(possiblevar2-possiblevalue2+pi,2*pi)-pi) \
                                            if ishinge[ipossiblesub+ipossiblesub2+1] else cond2 # + evalcond
                                #cond2 += cond
                                
                                if self.CheckExpressionUnique(handledconds, cond+cond2):
                                    # if the variables are both part of the rotation matrix and both zeros,
                                    # then we can deduce other rotation variables
                                    if self._iktype == 'transform6d' and \
                                       possiblevar in rotsymbols and \
                                       possiblevalue == S.Zero and \
                                       possiblevar2 in rotsymbols and \
                                       possiblevalue2 == S.Zero:
                                        
                                        checkexpr = [[cond + cond2], \
                                                     evalcond + evalcond2, \
                                                     possiblesub + possiblesub2, \
                                                     []]
                                        log.info('%r', flatzerosubstitutioneqs)
                                        log.info('%r', checkexpr)
                                        # exec(ipython_str, globals(), locals())
                                        flatzerosubstitutioneqs.append(checkexpr)
                                        localsubstitutioneqs.append(checkexpr)
                                        handledconds.append(cond+cond2)
                                        row1 = int(possiblevar.name[-2])
                                        col1 = int(possiblevar.name[-1])
                                        row2 = int(possiblevar2.name[-2])
                                        col2 = int(possiblevar2.name[-1])
                                        row3 = 3 - row1 - row2
                                        col3 = 3 - col1 - col2

                                        possiblevarname = possiblevar.name[:-2]
                                        
                                        if row1 == row2:
                                            # (row1, col3) is either 1 or -1, but don't know which.
                                            # know that (row1+1,col3) and (row1+2,col3) are zero though...
                                            checkexpr[2].append((Symbol('%s%d%d'%(possiblevarname, \
                                                                                  (row2+1)%3, col3)), S.Zero))
                                            checkexpr[2].append((Symbol('%s%d%d'%(possiblevarname, \
                                                                                  (row1+2)%3, col3)), S.Zero))
                                            # furthermore can defer that the left over 4 values are
                                            # [cos(ang), sin(ang), cos(ang), -sin(ang)] = abcd
                                            if row1 == 1:
                                                minrow = 0
                                                maxrow = 2
                                            else:
                                                minrow = (row1+1)%3
                                                maxrow = (row1+2)%3

                                            ra = Symbol('%s%d%d'%(possiblevarname, minrow, col1))
                                            rb = Symbol('%s%d%d'%(possiblevarname, minrow, col2))
                                            rc = Symbol('%s%d%d'%(possiblevarname, maxrow, col1))
                                            rd = Symbol('%s%d%d'%(possiblevarname, maxrow, col2))
                                            checkexpr[2].append((rb**2, S.One-ra**2))
                                            # need 3rd power since sympy cannot divide out the square
                                            checkexpr[2].append((rb**3, rb-rb*ra**2)) 
                                            checkexpr[2].append((rc**2, S.One-ra**2))
                                            #checkexpr[2].append((rc, -rb)) # not true
                                            #checkexpr[2].append((rd, ra))  # not true
                                            
                                        elif col1 == col2:
                                            # (row3, col1) is either 1 or -1, but don't know which.
                                            # know that (row3,col1+1) and (row3,col1+2) are zero though...
                                            checkexpr[2].append((Symbol('%s%d%d'%(possiblevarname, \
                                                                                  row3, (col1+1)%3)), S.Zero))
                                            checkexpr[2].append((Symbol('%s%d%d'%(possiblevarname, \
                                                                                  row3, (col1+2)%3)), S.Zero))
                                            # furthermore can defer that the left over 4 values are
                                            # [cos(ang), sin(ang), cos(ang), -sin(ang)] = abcd
                                            if col1 == 1:
                                                mincol = 0
                                                maxcol = 2
                                            else:
                                                mincol = (col1+1)%3
                                                maxcol = (col1+2)%3
                                                
                                            ra = Symbol('%s%d%d'%(possiblevarname, row1, mincol))
                                            rb = Symbol('%s%d%d'%(possiblevarname, row2, mincol))
                                            rc = Symbol('%s%d%d'%(possiblevarname, row1, maxcol))
                                            rd = Symbol('%s%d%d'%(possiblevarname, row2, maxcol))
                                            checkexpr[2].append((rb**2, S.One-ra**2))
                                            # need 3rd power since sympy cannot divide out the square
                                            checkexpr[2].append((rb**3, rb-rb*ra**2)) 
                                            checkexpr[2].append((rc**2, S.One-ra**2))
                                            #checkexpr[2].append((rc, -rb)) # not true
                                            #checkexpr[2].append((rd, ra))  # not true
                                            
                                        log.info('dual constraint %s\n' + \
                                                 '	in %s', \
                                                 ("\n"+" "*24).join(str(x) for x in list(checkexpr[2])), \
                                                 checkzero)
                                    else:
                                        # shouldn't have any rotation vars
                                        if not possiblevar in rotsymbols and not possiblevar2 in rotsymbols:
                                            checkexpr = [[cond + cond2], \
                                                         evalcond + evalcond2, \
                                                         possiblesub + possiblesub2, \
                                                         []]
                                            log.info('%r', flatzerosubstitutioneqs)
                                            log.info('%r', checkexpr)
                                            #exec(ipython_str, globals(), locals())
                                            flatzerosubstitutioneqs.append(checkexpr)
                                            localsubstitutioneqs.append(checkexpr)
                                            handledconds.append(cond + cond2)
                                            
                zerosubstitutioneqs[isolution] += localsubstitutioneqs
        # test the solutions
        
        # PREV: have to take the cross product of all the zerosubstitutioneqs in order to form stronger constraints on the equations because the following condition will be executed only if all SolverCheckZeros evalute to 0
        # NEW: not sure why cross product is necessary anymore....
        zerobranches = []
        accumequations = []
#         # since sequence_cross_product requires all lists to be non-empty, insert None for empty lists
#         for conditioneqs in zerosubstitutioneqs:
#             if len(conditioneqs) == 0:
#                 conditioneqs.append(None)
#         for conditioneqs in self.sequence_cross_product(*zerosubstitutioneqs):
#             validconditioneqs = [c for c in conditioneqs if c is not None]
#             if len(validconditioneqs) > 1:
#                 # merge the equations, be careful not to merge equations constraining the same variable
#                 cond = []
#                 evalcond = S.Zero
#                 othervarsubs = []
#                 dictequations = []
#                 duplicatesub = False
#                 for subcond, subevalcond, subothervarsubs, subdictequations in validconditioneqs:
#                     cond += subcond
#                     evalcond += abs(subevalcond)
#                     for subothervarsub in subothervarsubs:
#                         if subothervarsub[0] in [sym for sym,value in othervarsubs]:
#                             # variable is duplicated
#                             duplicatesub = True
#                             break
#                         othervarsubs.append(subothervarsub)
#                     if duplicatesub:
#                         break
#                     dictequations += subdictequations
#                 if not duplicatesub:
#                     flatzerosubstitutioneqs.append([cond,evalcond,othervarsubs,dictequations])

        trysubstitutions = self.ppsubs + self.npxyzsubs + self.rxpsubs \
                           if self._iktype == 'transform6d' or \
                              self._iktype == 'rotation3d' else \
                              self.ppsubs
            
        log.debug('c = %d, %d zero-substitution(s)', scopecounter, len(flatzerosubstitutioneqs))
        
        for iflatzerosubstitutioneqs, (cond, evalcond, othervarsubs, dictequations) in enumerate(flatzerosubstitutioneqs):
            # have to convert to fractions before substituting!
            if not all([self.isValidSolution(v) for s,v in othervarsubs]):
                continue
            
            othervarsubs = [(s,self.ConvertRealToRationalEquation(v)) for s,v in othervarsubs]
            #NewEquations = [eq.subs(self.npxyzsubs + self.rxpsubs).subs(othervarsubs) for eq in AllEquations]
            NewEquations = [eq.subs(othervarsubs) for eq in AllEquations]
            NewEquationsClean = self.PropagateSolvedConstants(NewEquations, othersolvedvars, curvars)
            
            try:
                # forcing a value, so have to check if all equations in NewEquations that do not contain
                # unknown variables are really 0
                extrazerochecks = []
                for i in range(len(NewEquations)):
                    expr = NewEquations[i]
                    if not self.isValidSolution(expr):
                        log.warn('not valid: %s',expr)
                        extrazerochecks = None
                        break
                    if not expr.has(*allvars) and \
                       self.CheckExpressionUnique(extrazerochecks,expr):
                        if expr.is_Symbol:
                            # can set that symbol to zero and create a new set of equations!
                            extrazerochecks.append(expr.subs(solsubs).evalf(n=30))
                            
                if extrazerochecks is not None:
                    newcases = set(currentcases)
                    
                    for singlecond in cond:
                        newcases.add(singlecond)
                        
                    if not self.degeneratecases.CheckCases(newcases):
                        log.info('depth = %d, c = %d, iter = %d/%d\n' +\
                                 '        start new cases: %s', \
                                 len(currentcases), scopecounter, iflatzerosubstitutioneqs, len(flatzerosubstitutioneqs), \
                                 ("\n"+" "*25).join(str(x) for x in list(newcases)))
                        
                        if len(NewEquationsClean) > 0:
                            newcasesubs = currentcasesubs + othervarsubs
                            # empty global symbols dictionary
                            self.globalsymbols = {}
                            for casesub in newcasesubs:
                                self._AddToGlobalSymbols(casesub[0], casesub[1])
                            extradictequations = []
                            
                            for s, v in trysubstitutions:
                                neweq = v.subs(newcasesubs)
                                if neweq != v:
                                    # should we make sure we're not adding it a second time?
                                    newcasesubs.append((s, neweq))
                                    extradictequations.append((s, neweq))
                                    self._AddToGlobalSymbols(s, neweq)
                                    
                            for var, eq in chain(originalGlobalSymbols.items(), dictequations):
                                neweq = eq.subs(othervarsubs)
                                if not self.isValidSolution(neweq):
                                    raise self.CannotSolveError(('equation %s is invalid ' + \
                                                                 'because of the following substitutions: ' + \
                                                                 '%s') % (eq, othervarsubs))
                                
                                if neweq == S.Zero:
                                    extradictequations.append((var, S.Zero))
                                self._AddToGlobalSymbols(var, neweq)
                                
                            if len(extradictequations) > 0:
                                # have to re-substitute since some equations evaluated to zero
                                NewEquationsClean = [eq.subs(extradictequations).expand() for eq in NewEquationsClean]
                                
                            newtree = self.SolveAllEquations(NewEquationsClean, \
                                                             curvars, \
                                                             othersolvedvars, \
                                                             solsubs, \
                                                             endbranchtree, \
                                                             currentcases = newcases, \
                                                             currentcasesubs = newcasesubs, \
                                                             unknownvars = unknownvars)
                            accumequations.append(NewEquationsClean) # store the equations for debugging purposes
                            
                        else:
                            log.info('no new equations! probably can freely determine %r', curvars)
                            # unfortunately cannot add curvars as a FreeVariable
                            # because all the remaining variables will have complex dependencies
                            # therefore, iterate a couple of jointevals
                            newtree = []
                            for curvar in curvars:
                                newtree.append(AST.SolverSolution(curvar.name, \
                                                                  jointeval = [S.Zero, pi/2, pi, -pi/2], \
                                                                  isHinge = self.IsHinge(curvar.name)))
                            newtree += endbranchtree
                            
                        zerobranches.append(([evalcond] + extrazerochecks, \
                                             newtree, \
                                             dictequations)) # what about extradictequations?

                        # print flatzerosubstitutioneqs
                        log.info('depth = %d, c = %d, iter = %d/%d\n' \
                                 + '        add new cases: %s', \
                                 len(currentcases), scopecounter, iflatzerosubstitutioneqs, len(flatzerosubstitutioneqs), \
                                 ("\n"+" "*23).join(str(x) for x in list(newcases)))
                        self.degeneratecases.AddCases(newcases)
                    else:
                        log.warn('already has handled cases %r', newcases)
                        
            except self.CannotSolveError, e:
                log.debug(e)
                continue
            finally:
                # restore the global symbols
                self.globalsymbols = originalGlobalSymbols

        if len(zerobranches) > 0:
            branchconds = AST.SolverBranchConds(zerobranches + \
                                                [(None, \
                                                  [AST.SolverBreak('branch miss %r'%curvars, \
                                                                   [(var,self._SubstituteGlobalSymbols(eq, originalGlobalSymbols)) \
                                                                    for var, eq in currentcasesubs], \
                                                                   othersolvedvars, \
                                                                   solsubs, \
                                                                   endbranchtree)], \
                                                [])])
            branchconds.accumequations = accumequations
            lastbranch.append(branchconds)
        else:            
            # add GuessValuesAndSolveEquations?
            lastbranch.append(AST.SolverBreak('no branches %r'%curvars, \
                                              [(var,self._SubstituteGlobalSymbols(eq, originalGlobalSymbols)) \
                                               for var, eq in currentcasesubs], \
                                              othersolvedvars, \
                                              solsubs, \
                                              endbranchtree))
            
        return prevbranch
    
    def GuessValuesAndSolveEquations(self, AllEquations, \
                                     curvars, \
                                     othersolvedvars, \
                                     solsubs, \
                                     endbranchtree, \
                                     currentcases = None, \
                                     unknownvars = None, \
                                     currentcasesubs = None):
        """
        There might be a degree of freedom but it is not trivial to compute.
        In this case we can try setting some "guess" value to a hinge variable.
        We try the highest hinge.

        Called by SolveAllEquations only.
        """
        
        scopecounter = int(self._scopecounter)
        hingevariables = [curvar for curvar in sorted(curvars,reverse=True) if self.IsHinge(curvar.name)]
        if len(hingevariables) > 0 and len(curvars) >= 2:
            curvar = hingevariables[0]
            leftovervars = list(curvars)
            leftovervars.remove(curvar)
            newtree = [AST.SolverConditionedSolution([])]
            zerovalues = []
            
            for jointeval in [S.Zero, pi/2, pi, -pi/2]:
                checkzeroequations = []
                NewEquations = []
                for eq in AllEquations:
                    neweq = eq.subs(curvar, jointeval)
                    neweqeval = neweq.evalf()
                    if neweq.is_number:
                        # if zero, then can ignore
                        if neweq == S.Zero:
                            continue
                        else:
                            # contradiciton! so jointeval is bad
                            NewEquations = None
                            break                        
                    if neweq.has(*leftovervars):
                        NewEquations.append(neweq)
                    else:
                        checkzeroequations.append(neweq)
                        
                if NewEquations is None: # neweq != S.Zero above and jointeval is bad
                    continue
                
                # check to make sure all leftover vars are in scope
                cansolve = True
                for leftovervar in leftovervars:
                    if not any([eq.has(leftovervar) for eq in NewEquations]):
                        cansolve = False
                        break
                if not cansolve:
                    continue
                if len(checkzeroequations) > 0:
                    solution = AST.SolverSolution(curvar.name, jointeval=[jointeval], isHinge=self.IsHinge(curvar.name))
                    solution.checkforzeros = checkzeroequations
                    solution.FeasibleIsZeros = True
                    newtree[0].solversolutions.append(solution)
                else:
                    # one value is enough
                    zerovalues.append(jointeval)
                    
            if len(zerovalues) > 0:
                # prioritize these solutions since they don't come with any extra checks
                solution = AST.SolverSolution(curvar.name, \
                                              jointeval = zerovalues, \
                                              isHinge = self.IsHinge(curvar.name))
                solution.FeasibleIsZeros = True
                newtree = [solution]
            elif len(newtree[0].solversolutions) == 0:
                # nothing found so remove the condition node
                newtree = []
                
            if len(newtree) > 0:
                log.warn('c = %d; ' + \
                         'we conjecture there is a free variable, but cannot figure it out by maths so set variable %s', \
                         scopecounter, curvar)
                newtree += self.SolveAllEquations(AllEquations, \
                                                  leftovervars, \
                                                  othersolvedvars+[curvar], \
                                                  solsubs + self.getVariable(curvar).subs, \
                                                  endbranchtree, \
                                                  currentcases = currentcases, \
                                                  currentcasesubs = currentcasesubs, \
                                                  unknownvars = unknownvars)
                return newtree

        if len(curvars) == 1:
            log.info('Have only one remaining variable %r and it is likely not in equations %r', curvars[0], AllEquations)
            solution = AST.SolverSolution(curvars[0].name, \
                                          jointeval = [S.Zero], \
                                          isHinge = self.IsHinge(curvars[0].name))
            solution.FeasibleIsZeros = True
            return [solution]+endbranchtree
        
        raise self.CannotSolveError('Cannot find a good variable to solve for')
    
    def SolvePairVariablesHalfAngle(self, raweqns, var0, var1, \
                                    othersolvedvars, tosubs = []):
        """
        Solves equations of two variables var0, var1 in sin, cos

        Called by SolvePairVariables, solveLiWoernleHiller.
        """
        import itertools

        varsym0 = self.getVariable(var0)
        varsym1 = self.getVariable(var1)
        varsyms = [varsym0, varsym1]
        
        unknownvars = [varsym0.cvar, varsym0.svar, \
                       varsym1.cvar, varsym1.svar  ]
        
        varsubs    = varsym0.subs    + varsym1.subs
        varsubsinv = varsym0.subsinv + varsym1.subsinv

        # numerators
        halftannum   = [1-varsym0.htvar**2, 2*varsym0.htvar, \
                        1-varsym1.htvar**2, 2*varsym1.htvar]
        # denominators
        halftandenom = [1+varsym0.htvar**2, 1+varsym1.htvar**2]

        # half tangent substitutions
        halftansubs = [(varsym0.cvar, halftannum[0] / halftandenom[0]), \
                       (varsym0.svar, halftannum[1] / halftandenom[0]), \
                       (varsym1.cvar, halftannum[2] / halftandenom[1]), \
                       (varsym1.svar, halftannum[3] / halftandenom[1])  ]
                       

        dummyvars = [self.getVariable(othervar) for othervar in othersolvedvars]
        dummyvars = list(chain.from_iterable([[v.cvar, v.svar, v.var, v.htvar] \
                                              for v in dummyvars]))

        trigsubs = [(varsym0.svar**2,               1-varsym0.cvar**2),  \
                    (varsym0.svar**3, varsym0.svar*(1-varsym0.cvar**2)), \
                    (varsym1.svar**2,               1-varsym1.cvar**2),  \
                    (varsym1.svar**3, varsym1.svar*(1-varsym1.cvar**2))]
            
        polyeqs = []
        for eq in raweqns:

            # unknownvars = [ cv1, sv1, cv2, sv2 ]
            peq = Poly(eq.subs(varsubs).subs(trigsubs).expand().subs(trigsubs), *unknownvars)
            
            if peq.has(varsym0.var) or peq.has(varsym1.var):
                raise self.CannotSolveError('Expecting only sin and cos! %s' % peq)
            
            maxdenom  = [ max([monoms[0]+monoms[1] for monoms in peq.monoms()]), \
                          max([monoms[2]+monoms[3] for monoms in peq.monoms()])  ]

            eqnew = S.Zero
            
            for monoms, c in peq.terms():
                """
                print '\nmonoms = ', monoms, '\nc = ', c
                term = c 
                for i in range(4): 
                    num, denom = fraction(halftansubs[i][1]) 
                    term *= num**monoms[i] 
 
                # the denoms for 0,1 and 2,3 are the same 
                for i in [0, 2]: 
                    denom = fraction(halftansubs[i][1])[1] 
                    term *= denom**(maxdenom[i/2] - monoms[i] - monoms[i+1])
                complexityvalue = self.codeComplexity(term.expand()) 
                if complexityvalue < 450: 
                    eqnew += simplify(term) 
                else: 
                    # too big, so don't simplify? 
                    eqnew += term    
                """

                term = c
                for i in range(4):
                    term *= halftannum[i]**monoms[i]
                    
                term *= halftandenom[0]**(maxdenom[0]-monoms[0]-monoms[1])
                term *= halftandenom[1]**(maxdenom[1]-monoms[2]-monoms[3])

                # too big, so don't simplify?
                eqnew += simplify(term) if self.codeComplexity(term.expand())<450 else term
                
            polyeq = Poly(eqnew, varsym0.htvar, varsym1.htvar)

            
            if polyeq.TC() == S.Zero:
                log.info('Trailing coefficient is 0, i.e. for power (0,0)')


                # might be able to divide out variables?
                minmonoms = None
                for monom in polyeq.monoms():
                    if minmonoms is None:
                        minmonoms = list(monom)
                    else:
                        for i in range(len(minmonoms)):
                            minmonoms[i] = min(minmonoms[i],monom[i])
                            
                newpolyeq = Poly(S.Zero, *polyeq.gens)
                for m, c in polyeq.terms():
                    newm = list(m)
                    for i in range(len(minmonoms)):
                        newm[i] -= minmonoms[i]
                    newpolyeq = newpolyeq.add(Poly.from_dict({tuple(newm):c}, \
                                                             *newpolyeq.gens))
                    
                log.warn('Converting "%s" to "%s"' % (polyeq, newpolyeq))
                
                # check if any equations are only in one variable
                polyeq = newpolyeq
                
            polyeqs.append(polyeq)

        eqtosolve = self.sortComplexity([e.as_expr() for e in polyeqs if not e.has(varsym1.htvar)])
        if len(eqtosolve)>0:
            try:
                return self.solveSingleVariable(eqtosolve, varsym0.var, \
                                                othersolvedvars, unknownvars = [])
            except self.CannotSolveError:
                log.info('solveSingleVariable fails to solve %r', varsym0.var)

        else:
            eqtosolve = self.sortComplexity([e.as_expr() for e in polyeqs if not e.has(varsym0.htvar)])
            if len(eqtosolve)>0:
                try:
                    return self.solveSingleVariable(eqtosolve, varsym1.var, \
                                            othersolvedvars, unknownvars = [])
                except self.CannotSolveError:
                    log.info('solveSingleVariable fails to solve %r', varsym1.var)

        #complexity = [(self.codeComplexity(peq.as_expr()), peq) for peq in polyeqs]
        #complexity.sort(key = itemgetter(0))
        #polyeqs = [peq[1] for peq in complexity]

        polyeqs.sort(key = lambda x: x.count_ops())
        
        solutions = [None, None]
        linearsolution = None
        for ileftvar in (0,1):
            if linearsolution is not None:
                log.info('Found linearsolution')
                break
            
            leftvar = varsyms[ileftvar].htvar
            newpolyeqs = [Poly(eq,varsyms[1-ileftvar].htvar) for eq in polyeqs]
            maxdeglist = [max(peq.degree_list()) for peq in newpolyeqs]
            mindegree = __builtin__.min(maxdeglist)
            maxdegree = __builtin__.max(maxdeglist)

            for peq in newpolyeqs:
                if len(peq.monoms()) == 1:
                    possiblefinaleq = self.checkFinalEquation(Poly(peq.LC(), leftvar), tosubs)
                    if possiblefinaleq is not None:
                        solutions[ileftvar] = [possiblefinaleq]
                        break
                    
            for degree in range(mindegree, maxdegree+1):
                if not (solutions[ileftvar] is None and linearsolution is None):
                    log.info('Found either solutions[ileftvar] or linearsolution')
                    break
                
                newpolyeqs2 = [peq for peq in newpolyeqs \
                               if max(peq.degree_list()) <= degree]
                
                if degree+1 <= len(newpolyeqs2):
                    # To avoid wrong solutions, we get resultants for all equations
                    possibilities = []
                    unusedindices = range(len(newpolyeqs2))
                    for eqsindices in combinations(range(len(newpolyeqs2)), degree+1):
                        Mall = zeros((degree+1, degree+1))
                        totalcomplexity = 0
                        for i, eqindex in enumerate(eqsindices):
                            eq = newpolyeqs2[eqindex]
                            for j, c in eq.terms():
                                totalcomplexity += self.codeComplexity(c.expand())
                                Mall[i, j[0]] = c
                        if degree >= 4 and totalcomplexity > 5000:
                            # the determinant will never finish otherwise
                            continue
                        # det_bareis freezes when there are huge fractions
                        # det=self.det_bareis(Mall, *(self.pvars+dummyvars+[leftvar]))
                        # for i in range(Mall.shape[0]):
                        #      for j in range(Mall.shape[1]):
                        #          Mall[i,j] = Poly(Mall[i,j], leftvar)

                        try:
                            Malldet = Mall.berkowitz_det()
                            # log.info('Try simplifying det(Mall)')
                            # Malldet = simplify(Malldet)
                            # print Malldet
                            # log.info('Finished simplifying of of det(Mall)')
                        except Exception, e:
                            log.warn('Failed to compute det(Mall): %s', e)
                            continue
                        
                        complexity = self.codeComplexity(Malldet)
                        if complexity > 1200:
                            log.warn('Complexity of det(Mall) is too big: %d', complexity)
                            continue

                        possiblefinaleq = self.checkFinalEquation(Poly(Malldet, leftvar), tosubs)
                        if possiblefinaleq is not None:
                            # sometimes +- I are solutions, so remove them
                            q, r = div(possiblefinaleq, leftvar+I)
                            if r == S.Zero:
                                possiblefinaleq = Poly(q, leftvar)
                                
                            q, r = div(possiblefinaleq, leftvar-I)
                            if r == S.Zero:
                                possiblefinaleq = Poly(q, leftvar)
                                
                            possibilities.append(possiblefinaleq)

                            unusedindices = [ind for ind in unusedindices if ind not in eqsindices]
                            if len(unusedindices) == 0:
                                break # break degree in range(mindegree, maxdegree+1)

                    if len(possibilities) > 1:
                        try:
                            linearsolutions = self.solveVariablesLinearly(possibilities, othersolvedvars)
                            log.info('solveVariablesLinearly has found some solution for %r', possibilities[0].gens)
                            
                            # if can solve for a unique solution linearly, then prioritize this over anything
                            prevsolution = AST.SolverBreak('SolvePairVariablesHalfAngle fail')
                            for divisor, linearsolution in linearsolutions:
                                assert(len(linearsolution)==1)
                                divisorsymbol = self.gsymbolgen.next()
                                
                                # Call AST SolverSolution constructor
                                solversolution = AST.SolverSolution(varsyms[ileftvar].name, \
                                                                    jointeval = [2*atan(linearsolution[0]/divisorsymbol)], \
                                                                    isHinge = self.IsHinge(varsyms[ileftvar].name))
                                
                                # Call AST SolverCheckZeros constructor
                                prevsolution = AST.SolverCheckZeros(varsyms[ileftvar].name, \
                                                                    [divisorsymbol], \
                                                                    zerobranch    = [prevsolution], \
                                                                    nonzerobranch = [solversolution], \
                                                                    thresh = 1e-6)
                                    
                                prevsolution.dictequations = [(divisorsymbol, divisor)]
                            linearsolution = prevsolution
                            break
                        
                        except self.CannotSolveError:
                            log.info('solveVariablesLinearly failed to find %r', possibilities[0].gens)
                            pass
                        
                    if len(possibilities) > 0:
                        # sort with respect to degree
                        # equationdegrees = [(max(peq.degree_list())*100000 + \
                        #                     self.codeComplexity(peq.as_expr()), peq) \
                        #                    for peq in possibilities]
                        # equationdegrees.sort(key=itemgetter(0))
                        possibilities.sort(key = lambda x: max(x.degree_list())*100000+x.count_ops())
                        solutions[ileftvar] = possibilities # [peq[1] for peq in equationdegrees]
                        break
                    
        if linearsolution is not None:
            log.info('SolvePairVariablessHalfAngle returns linear solution.')
            return [linearsolution]
        
        # take the solution with the smallest degree
        if solutions[0] is None:
            if solutions[1] is None:

                log.info('SolvePairVariablesHalfAngle has not found any solution yet.')
                
                pfinals  = None
                ileftvar = None
            else:
                pfinals  = solutions[1]
                ileftvar = 1
        elif solutions[1] is None:
            pfinals = solutions[0]
            ileftvar = 0
        else:
            maxsol0deglist = max(solutions[0][0].degree_list())
            maxsol1deglist = max(solutions[1][0].degree_list())
                
            if   maxsol1deglist < maxsol0deglist:
                pfinals = solutions[1]
                ileftvar = 1
            elif maxsol1deglist == maxsol1deglist and \
                 self.codeComplexity(solutions[1][0].as_expr()) < self.codeComplexity(solutions[0][0].as_expr()):
                pfinals = solutions[1]
                ileftvar = 1
            else:
                pfinals = solutions[0]
                ileftvar = 0
                
        dictequations = []
        if pfinals is None:
            #simplifyfn = self._createSimplifyFn(self.freevarsubs, self.freevarsubsinv)
            for newreducedeqs in combinations(polyeqs, 2):
                try:
                    Mall = None
                    numrepeating = None
                    for ileftvar in range(2):
                        # TODO, sometimes this works and sometimes this doesn't
                        try:
                            Mall, allmonoms = self.solveDialytically(newreducedeqs, ileftvar, \
                                                                     returnmatrix = True)
                            if Mall is not None:
                                leftvar = polyeqs[0].gens[ileftvar]
                                break
                        except self.CannotSolveError, e:
                            log.debug(e)
                        
                    if Mall is None:
                        continue
                    
                    shape=Mall[0].shape
                    assert(shape[0] == 4 and shape[1] == 4)
                    Malltemp = [None]*len(Mall)
                    M = zeros(shape)
                    for idegree in range(len(Mall)):
                        Malltemp[idegree] = zeros(shape)
                        for i in range(shape[0]):
                            for j in range(shape[1]):
                                if Mall[idegree][i,j] != S.Zero:
                                    if self.codeComplexity(Mall[idegree][i,j]) > 5:
                                        sym = self.gsymbolgen.next()
                                        Malltemp[idegree][i,j] = sym
                                        dictequations.append((sym, Mall[idegree][i, j]))
                                    else:
                                        Malltemp[idegree][i, j] = Mall[idegree][i, j]
                        M += Malltemp[idegree]*leftvar**idegree
                        
                    tempsymbols = [self.gsymbolgen.next() for i in range(16)]
                    tempsubs = []
                    for i in range(16):
                        if M[i] != S.Zero:
                            tempsubs.append((tempsymbols[i], Poly(M[i], leftvar)))
                        else:
                            tempsymbols[i] = S.Zero
                            
                    Mtemp = Matrix(4, 4,tempsymbols)                    
                    dettemp = Mtemp.det()
                    log.info('Multiplying all determinant coefficients to solve %s', leftvar)
                    eqadds = []
                    for arg in dettemp.args:
                        eqmuls = [Poly(arg2.subs(tempsubs), leftvar) for arg2 in arg.args]
                        
                        if sum(eqmuls[0].degree_list()) == 0:
                            eq = eqmuls.pop(0)
                            eqmuls[0] = eqmuls[0]*eq
                            
                        while len(eqmuls) > 1:
                            ioffset = 0
                            eqmuls2 = []
                            while ioffset < len(eqmuls)-1:
                                eqmuls2.append(eqmuls[ioffset]*eqmuls[ioffset+1])
                                ioffset += 2
                            eqmuls = eqmuls2
                        eqadds.append(eqmuls[0])
                        
                    log.info('Done multiplying all determinant. Now convert to Poly')
                    det = Poly(S.Zero, leftvar)
                    for ieq, eq in enumerate(eqadds):
                        # log.info('adding to det %d/%d', ieq, len(eqadds))
                        det += eq
                        
                    if len(Mall) <= 3:
                        # need to simplify further since self.globalsymbols can have important substitutions
                        # that can yield the entire determinant to zero
                        log.info('Attempt to simplify determinant')
                        newdet = Poly(S.Zero,leftvar)
                        for m, c in det.terms():
                            origComplexity = self.codeComplexity(c)
                            # 100 is a guess
                            if origComplexity < 100:
                                neweq = c.subs(dictequations)
                                if self.codeComplexity(neweq) < 100:
                                    neweq = self._SubstituteGlobalSymbols(neweq).expand()
                                    newComplexity = self.codeComplexity(neweq)
                                    if newComplexity < origComplexity:
                                        c = neweq
                            newdet += c*leftvar**m[0]
                        det = newdet
                    if det.degree(0) <= 0:
                        continue
                    pfinals = [det]
                    break
                except self.CannotSolveError,e:
                    log.debug(e)
        else:
            log.info('SolvePairVariableHalfAngle has found %d pfinals', len(pfinals))
                    
                    
        if pfinals is None:
            raise self.CannotSolveError('SolvePairVariablesHalfAngle: ' + \
                                        'failed to solve dialytically with %d equations'% \
                                        len(polyeqs))

        log.info('SolvePairVariableHalfAngle has found %d pfinal(s)', len(pfinals))
        jointsol = 2*atan(varsyms[ileftvar].htvar)

        # Call AST SolverPolynomialRoots constructor
        solution = AST.SolverPolynomialRoots(jointname = varsyms[ileftvar].name, \
                                             poly      = pfinals[0], \
                                             jointeval = [jointsol], \
                                             isHinge   = self.IsHinge(varsyms[ileftvar].name))
        solution.checkforzeros = []
        solution.postcheckforzeros = []
        if len(pfinals) > 1:
            # verify with at least one solution
            solution.postcheckfornonzeros = [peq.as_expr() for peq in pfinals[1:2]]
            solution.polybackup = pfinals[1]
        solution.postcheckforrange = []
        solution.dictequations = dictequations
        solution.postcheckfornonzerosThresh = 1e-7
        # make threshold a little loose since can be a lot of numbers compounding.
        # depending on the degree, can expect small coefficients to be still valid
        solution.AddHalfTanValue = True

        log.info('End of SolvePairVariableHalfAngle')
        
        return [solution]

                
    def solveVariablesLinearly(self, polyeqs, othersolvedvars, maxsolvabledegree = 4):
        """
        Called by SolvePairVariablesHalfAngle only.
        """

        log.debug('solveVariablesLinearly:\n' + \
                  '        solvevariables  = %r\n' + \
                  '        othersolvedvars = %r', \
                  polyeqs[0].gens, \
                  othersolvedvars)

        # number of monomials (excluding constants) in each equation
        nummonoms = [(len(peq.monoms()) - int(peq.TC()!=S.Zero)) for peq in polyeqs]
        mindegree = __builtin__.min(nummonoms)
        maxdegree = min(__builtin__.max(nummonoms), len(polyeqs))

        polyeqs.sort(key = lambda x: x.count_ops())
        #complexity = [(self.codeComplexity(peq.as_expr()),peq) for peq in polyeqs]
        #complexity.sort(key=itemgetter(0))
        #polyeqs = [peq[1] for peq in complexity]
        
        v = [self.getVariable(othervar) for othervar in othersolvedvars]
        trigsubs            = list(chain.from_iterable([var.subs    for var in v]))
        trigsubsinv         = list(chain.from_iterable([var.subsinv for var in v]))
        othersolvedvarssyms = list(chain.from_iterable([var.vars    for var in v]))

        symbolscheck = []
        for i,solvevar in enumerate(polyeqs[0].gens):
            monom = [0]*len(polyeqs[0].gens)
            monom[i] = 1
            symbolscheck.append(tuple(monom))
        solutions = []
        
        for degree in range(mindegree, maxdegree+1):
            allindices = [i for i, n in enumerate(nummonoms) if n <= degree]
            if len(allindices) >= degree:
                allmonoms = set()
                for index in allindices:
                    allmonoms = allmonoms.union(set(polyeqs[index].monoms()))
                allmonoms = list(allmonoms)
                allmonoms.sort()
                if __builtin__.sum(allmonoms[0]) == 0:
                    allmonoms.pop(0)
                # allmonoms has to have symbols as a single variable
                if not all([check in allmonoms for check in symbolscheck]):
                    continue
                
                if len(allmonoms) == degree:
                    if degree > maxsolvabledegree:
                        log.warn('Cannot handle linear solving for more than 4 equations')
                        continue
                    
                    systemequations = []
                    consts = []
                    for index in allindices:
                        pdict = polyeqs[index].as_dict()
                        systemequations.append([pdict.get(monom, S.Zero) for monom in allmonoms])
                        consts.append(-polyeqs[index].TC())
                        
                    # generate at least two solutions in case the first solution has determinant = 0
                    solutions = []
                    for startrow in range(len(systemequations)):
                        rows = [startrow]
                        M = Matrix(1,len(allmonoms), systemequations[rows[0]])
                        for i in range(startrow+1, len(systemequations)):
                            numequationsneeded = M.shape[1] - M.shape[0]
                            if i+numequationsneeded > len(systemequations):
                                # cannot do anything
                                break
                            mergedsystemequations = list(systemequations[i])
                            for j in range(1,numequationsneeded):
                                mergedsystemequations += systemequations[i+j]
                            M2 = M.col_join(Matrix(numequationsneeded,len(allmonoms),mergedsystemequations))
                            complexity = 0
                            for i2 in range(M2.rows):
                                for j2 in range(M2.cols):
                                    complexity += self.codeComplexity(M2[i2, j2])
                            if self.IsDeterminantNonZeroByEval(M2):
                                if complexity < 5000:
                                    Mdet = M2.det()
                                    if Mdet != S.Zero:
                                        M = M2
                                        for j in range(numequationsneeded):
                                            rows.append(i+j)
                                        break
                                else:
                                    log.warn('Found solution, but matrix is too complex and ' + \
                                             'determinant will most likely freeze (%d)', complexity)
                                
                        if M.shape[0] == M.shape[1]:
                            Mdet = self.trigsimp_new(Mdet.subs(trigsubsinv)).subs(trigsubs)
                            #Minv = M.inv()
                            B = Matrix(M.shape[0],1,[consts[i] for i in rows])
                            Madjugate = M.adjugate()
                            solution = []
                            for check in symbolscheck:
                                value = Madjugate[allmonoms.index(check),:]*B
                                solution.append(self.trigsimp_new(value[0].subs(trigsubsinv)).subs(trigsubs))
                            solutions.append([Mdet, solution])
                            if len(solutions) >= 2:
                                break
                    if len(solutions) > 0:
                        break

        if len(solutions) == 0:            
            raise self.CannotSolveError('solveVariablesLinearly failed')
        else:
            # exec(ipython_str)
            log.info('solveVariablesLinearly has found some solution.')
        
        return solutions

    def solveSingleVariableLinearly(self,raweqns,solvevar,othervars,maxnumeqs=2,douniquecheck=True):
        """
        Solves a linear system for one variable, assuming everything else is constant.

        Need >=3 equations.

        Called by SolvePairVariables only.
        """
        cvar = Symbol('c%s' % solvevar.name)
        svar = Symbol('s%s' % solvevar.name)
        varsubs = [(cos(solvevar), cvar), (sin(solvevar), svar)]
        othervarsubs = [(sin(v)**2, 1-cos(v)**2) for v in othervars]
        eqpolys = [Poly(eq.subs(varsubs), cvar, svar) for eq in raweqns]
        eqpolys = [eq for eq in eqpolys if sum(eq.degree_list()) == 1 and not eq.TC().has(solvevar)]
        #eqpolys.sort(lambda x,y: iksolver.codeComplexity(x) - iksolver.codeComplexity(y))
        partialsolutions = []
        neweqs = []
        for p0,p1 in combinations(eqpolys, 2):
            p0dict = p0.as_dict()
            p1dict = p1.as_dict()

            M = Matrix(2, 3, \
                       [p0dict.get((1,0), S.Zero), \
                        p0dict.get((0,1), S.Zero), p0.TC(), \
                        p1dict.get((1,0), S.Zero), \
                        p1dict.get((0,1), S.Zero), p1.TC()])
            M = M.subs(othervarsubs).expand()
            
            partialsolution = [-M[1,1]*M[0,2]+M[0,1]*M[1,2], \
                               M[1,0]*M[0,2]-M[0,0]*M[1,2] , \
                               M[0,0]*M[1,1]-M[0,1]*M[1,0]]
            
            partialsolution = [eq.expand().subs(othervarsubs).expand() for eq in partialsolution]
            rank = [self.codeComplexity(eq) for eq in partialsolution]
            partialsolutions.append([rank, partialsolution])
            # cos(A)**2 + sin(A)**2 - 1 = 0, useful equation but the squares introduce wrong solutions
            #neweqs.append(partialsolution[0]**2+partialsolution[1]**2-partialsolution[2]**2)
        # try to cross
        partialsolutions.sort(lambda x, y: int(min(x[0])-min(y[0])))
        for (rank0,ps0),(rank1,ps1) in combinations(partialsolutions,2):
            if self.equal(ps0[0]*ps1[2]-ps1[0]*ps0[2],S.Zero):
                continue
            neweqs.append(ps0[0]*ps1[2]-ps1[0]*ps0[2])
            neweqs.append(ps0[1]*ps1[2]-ps1[1]*ps0[2])
            # probably a linear combination of the first two
            #neweqs.append(ps0[0]*ps1[1]-ps1[0]*ps0[1])
            # too long
            #neweqs.append(ps0[0]*ps1[0]+ps0[1]*ps1[1]-ps0[2]*ps1[2])
            if len(neweqs) >= maxnumeqs:
                break;
            
        neweqs2 = [eq.expand().subs(othervarsubs).expand() for eq in neweqs]
        
        if douniquecheck:
            reducedeqs = []
            i = 0
            while i < len(neweqs2):
                reducedeq = self.removecommonexprs(neweqs2[i])
                if neweqs2[i] != S.Zero and \
                   self.CheckExpressionUnique(reducedeqs, reducedeq):
                    reducedeqs.append(reducedeq)
                    i += 1
                else:
                    eq = neweqs2.pop(i)
        return neweqs2

    def solveHighDegreeEquationsHalfAngle(self, lineareqs, varsym, tosubs = []):
        """
        Solve a set of equations in one variable with half-angle substitution.

        Called by SolvePairVariables, SolveAllEquations, solveSingleVariable.
        """
        
        dummysubs = [(varsym.cvar, (1-varsym.htvar**2)/(1+varsym.htvar**2)), \
                     (varsym.svar,      2*varsym.htvar/(1+varsym.htvar**2))]
        polyeqs = []
        
        for eq in lineareqs:
            trigsubs = [(varsym.svar**2, 1-varsym.cvar**2), \
                        (varsym.svar**3, varsym.svar*(1-varsym.cvar**2))]
            try:
                peq = Poly(eq.subs(varsym.subs).subs(trigsubs), varsym.cvar, varsym.svar)
                
            except PolynomialError, e:
                raise self.CannotSolveError('solveHighDegreeEquationsHalfAngle: poly error (%r)' % eq)
            
            if peq.has(varsym.var):
                raise self.CannotSolveError('solveHighDegreeEquationsHalfAngle: expecting only sin and cos! %s' % peq)
            
            if sum(peq.degree_list()) == 0:
                continue
            
            # check if all terms are multiples of cos/sin
            maxmonoms = [0, 0]
            maxdenom = 0
            for monoms in peq.monoms():
                for i in range(2):
                    maxmonoms[i] = max(maxmonoms[i], monoms[i])
                maxdenom = max(maxdenom,monoms[0] + monoms[1])
            eqnew = S.Zero
            for monoms,c in peq.terms():
                if c.evalf() != S.Zero: # big fractions might make this difficult to reduce to 0
                    term = c
                    for i in range(2):
                        num,denom = fraction(dummysubs[i][1])
                        term *= num**monoms[i]
                    # the denoms for 0,1 and 2,3 are the same
                    denom = fraction(dummysubs[0][1])[1]
                    term *= denom**(maxdenom-monoms[0]-monoms[1])
                    eqnew += simplify(term)
            polyeqs.append(Poly(eqnew, varsym.htvar))

        for peq in polyeqs:
            # do some type of resultants, for now just choose first polynomial
            finaleq = simplify(peq.as_expr()).expand()
            pfinal = Poly(self.removecommonexprs(finaleq), \
                          varsym.htvar)
            pfinal = self.checkFinalEquation(pfinal, tosubs)
            if pfinal is not None and pfinal.degree(0) > 0:
                jointsol = 2*atan(varsym.htvar)
                solution = AST.SolverPolynomialRoots(jointname = varsym.name, \
                                                     poly = pfinal, \
                                                     jointeval = [jointsol], \
                                                     isHinge = self.IsHinge(varsym.name))
                solution.AddHalfTanValue      = True
                solution.checkforzeros        = []
                solution.postcheckforzeros    = []
                solution.postcheckfornonzeros = []
                solution.postcheckforrange    = []
                return solution

        raise self.CannotSolveError(('half-angle substitution for joint %s failed, ' + \
                                    '%d equations examined') % (varsym.var, len(polyeqs)))

    def solveSingleVariable(self, raweqns, \
                            var, othersolvedvars, \
                            maxsolutions = 4, maxdegree = 2, \
                            subs = None, \
                            unknownvars = None):
        """
        Called by 

        SolvePrismaticHingePairVariables
        SolvePairVariables
        solveFullIK_TranslationAxisAngle4D
        SolveAllEquations
        AddSolution
        SolvePairVariablesHalfAngle
        """
        
        varsym = self.getVariable(var)
        vars = [varsym.cvar, varsym.svar, varsym.htvar, var]
        othersubs = []
        for othersolvedvar in othersolvedvars:
            othersubs += self.getVariable(othersolvedvar).subs

#         eqns = []
#         for eq in raweqns:
#             if eq.has(*vars):
#                 # for equations that are very complex, make sure at least one set of values yields a non zero equation
#                 testeq = eq.subs(varsym.subs+othersubs)
#                 if any([testeq.subs(testconsistentvalue).evalf()!=S.Zero for testconsistentvalue in self.testconsistentvalues]):
#                     eqns.append(eq)
        eqns = [eq.expand() for eq in raweqns if eq.has(*vars)]
        if len(eqns) == 0:
            raise self.CannotSolveError('not enough equations')

        # prioritize finding a solution when var is alone
        returnfirstsolutions = []
        
        for eq in eqns:
            symbolgen = cse_main.numbered_symbols('const')
            eqnew, symbols = self.groupTerms(eq.subs(varsym.subs), vars, symbolgen)
            try:
                ps = Poly(eqnew, varsym.svar)
                pc = Poly(eqnew, varsym.cvar)
                if sum(ps.degree_list()) > 0 or \
                   sum(pc.degree_list()) > 0 or \
                   ps.TC() == S.Zero or \
                   pc.TC() == S.Zero:
                    continue
                
            except PolynomialError:
                continue
            
            numvar = self.countVariables(eqnew, var)
            if numvar in [1, 2]:
                try:
                    tempsolutions  = solve(eqnew, var)

                    # TGN: ensure curvars is a subset of self.trigvars_subs
                    assert(all([z in self.trigvars_subs for z in othersolvedvars]))
                    
                    jointsolutions = [self.SimplifyTransform(self.trigsimp_new(s.subs(symbols))) \
                                      for s in tempsolutions]
                    if all([self.isValidSolution(s) and s != S.Zero \
                            for s in jointsolutions]) and \
                                len(jointsolutions) > 0:
                        # check if any solutions don't have divide by zero problems
                        returnfirstsolutions.append(AST.SolverSolution(var.name, \
                                                                       jointeval = jointsolutions,\
                                                                       isHinge = self.IsHinge(var.name)))
                        hasdividebyzero = any([len(self.checkForDivideByZero(self._SubstituteGlobalSymbols(s))) > 0 \
                                               for s in jointsolutions])
                        if not hasdividebyzero:
                            return returnfirstsolutions
                        
                except NotImplementedError, e:
                    # when solve cannot solve an equation
                    log.warn(e)
            
            numvar = self.countVariables(eqnew, varsym.htvar)
            if Poly(eqnew, varsym.htvar).TC() != S.Zero and numvar in [1, 2]:
                try:
                    tempsolutions = solve(eqnew,varsym.htvar)
                    jointsolutions = []
                    for s in tempsolutions:
                        s2 = s.subs(symbols)

                        # TGN: ensure curvars is a subset of self.trigvars_subs
                        assert(all([z in self.trigvars_subs for z in othersolvedvars]))
                    
                        s3 = self.trigsimp_new(s2)
                        s4 = self.SimplifyTransform(s3)
                        try:
                            jointsolutions.append(2*atan(s4, evaluate = False))
                            # set evaluate to False; otherwise it takes long time to evaluate when s4 is a number
                        except RuntimeError, e:
                            log.warn('got runtime error when taking atan: %s', e)
                            
                    if all([self.isValidSolution(s) \
                            and s != S.Zero \
                            for s in jointsolutions]) \
                                and len(jointsolutions) > 0:
                        returnfirstsolutions.append(AST.SolverSolution(var.name, \
                                                                       jointeval = jointsolutions, \
                                                                       isHinge = self.IsHinge(var.name)))
                        hasdividebyzero = any([len(self.checkForDivideByZero(self._SubstituteGlobalSymbols(s))) > 0 \
                                               for s in jointsolutions])
                        if not hasdividebyzero:
                            return returnfirstsolutions

                except NotImplementedError, e:
                    # when solve cannot solve an equation
                    log.warn(e)
                    
        if len(returnfirstsolutions) > 0:
            # already computed some solutions, so return them
            # note that this means that all solutions have a divide-by-zero condition
            return returnfirstsolutions
        
        solutions = []
        if len(eqns) > 1:
            neweqns = []
            listsymbols = []
            symbolgen = cse_main.numbered_symbols('const')
            for e in eqns:
                enew, symbols = self.groupTerms(e.subs(varsym.subs), \
                                                [varsym.cvar, varsym.svar, var], \
                                                symbolgen)
                try:
                    # remove coupled equations
                    if any([(m[0]>0) + (m[1]>0) + (m[2]>0) > 1 \
                            for m in Poly(enew, varsym.cvar, varsym.svar, var).monoms()]):
                        continue
                except PolynomialError:
                    continue
                
                try:
                    # ignore any equations with degree 3 or more
                    if max(Poly(enew, varsym.svar).degree_list()) > maxdegree or \
                       max(Poly(enew, varsym.cvar).degree_list()) > maxdegree:
                        log.debug('ignoring equation: ', enew)
                        continue
                except PolynomialError:
                    continue
                
                try:
                    if Poly(enew,varsym.svar).TC() == S.Zero or \
                       Poly(enew,varsym.cvar)      == S.Zero or \
                       Poly(enew,varsym.var)       == S.Zero:
                        log.debug('%s allows trivial solution for %s, ignore', e, varsym.name)
                        continue
                except PolynomialError:
                    continue
                
                rank = self.codeComplexity(enew)
                for s in symbols:
                    rank += self.codeComplexity(s[1])
                neweqns.append((rank, enew))
                listsymbols += symbols
                
            # We only need two equations for two variables, so we sort all equations and
            # start with the least complicated ones until we find a solution
            eqcombinations = []
            for eqs in combinations(neweqns,2):
                eqcombinations.append((eqs[0][0] + eqs[1][0], [Eq(e[1], 0) for e in eqs]))
            eqcombinations.sort(lambda x, y: x[0]-y[0])
            hasgoodsolution = False
            for icomb,comb in enumerate(eqcombinations):
                # skip if too complex
                if len(solutions) > 0 and comb[0] > 200:
                    break
                # try to solve for both sin and cos terms
                if not (self.has(comb[1], varsym.svar) and \
                        self.has(comb[1], varsym.cvar)):
                    continue
                
                try:
                    s = solve(comb[1], [varsym.svar, varsym.cvar])
                except (PolynomialError, CoercionFailed), e:
                    log.debug('solveSingleVariable: failed: %s', e)
                    continue
                
                if s is not None:
                    sollist = [(s[varsym.svar], s[varsym.cvar])] if \
                              s.has_key(varsym.svar) and \
                              s.has_key(varsym.cvar) else [] if \
                              hasattr(s, 'has_key') else s
                    
                    # sollist = None
                    # if hasattr(s, 'has_key'):
                    #     if s.has_key(varsym.svar) and \
                    #        s.has_key(varsym.cvar):
                    #         sollist = [(s[varsym.svar], s[varsym.cvar])]
                    #     else:
                    #         sollist = []
                    # else:
                    #     sollist = s
                        
                    solversolution = AST.SolverSolution(var.name,jointeval = [], \
                                                        isHinge = self.IsHinge(var.name))
                    goodsolution = 0
                    for svarsol, cvarsol in sollist:
                        # solutions cannot be trivial
                        soldiff = (svarsol-cvarsol).subs(listsymbols)
                        soldiffComplexity = self.codeComplexity(soldiff)
                        if soldiffComplexity < 1000 and soldiff.expand() == S.Zero:
                            break
                        
                        svarComplexity = self.codeComplexity(svarsol.subs(listsymbols))
                        cvarComplexity = self.codeComplexity(cvarsol.subs(listsymbols))
                        
                        if  svarComplexity < 600 and \
                            svarsol.subs(listsymbols).expand() == S.Zero and \
                            cvarComplexity < 600 and \
                            Abs(cvarsol.subs(listsymbols).expand()) != S.One:
                            # TGN: this used to be ... - S.One != S.Zero
                            break
                        
                        if cvarComplexity < 600 and \
                           cvarsol.subs(listsymbols).expand() == S.Zero and \
                           svarComplexity < 600 and \
                           Abs(svarsol.subs(listsymbols).expand()) != S.One:
                            # TGN: this used to be ... - S.One != S.Zero
                            break
                        
                        # check the numerator and denominator if solutions are the same or for possible divide by zeros
                        svarfrac = fraction(svarsol)
                        svarfrac = [svarfrac[0].subs(listsymbols), \
                                    svarfrac[1].subs(listsymbols)]
                        cvarfrac = fraction(cvarsol)
                        cvarfrac = [cvarfrac[0].subs(listsymbols), \
                                    cvarfrac[1].subs(listsymbols)]
                        
                        if self.equal(svarfrac[0], cvarfrac[0]) and \
                           self.equal(svarfrac[1], cvarfrac[1]):
                            break
                        
                        if not (\
                                self.isValidSolution(svarfrac[0]) and \
                                self.isValidSolution(svarfrac[1]) and \
                                self.isValidSolution(cvarfrac[0]) and \
                                self.isValidSolution(cvarfrac[1]) ):
                            continue
                        
                        # check if there exists at least one test solution with non-zero denominators
                        if subs is None:
                            testeqs = [svarfrac[1].subs(othersubs), \
                                       cvarfrac[1].subs(othersubs)]
                        else:
                            testeqs = [svarfrac[1].subs(subs).subs(othersubs), \
                                       cvarfrac[1].subs(subs).subs(othersubs)]
                        testsuccess = False
                        
                        for testconsistentvalue in self.testconsistentvalues:
                            if all([testeq.subs(self.globalsymbols).subs(testconsistentvalue).evalf() != S.Zero \
                                    for testeq in testeqs]):
                                testsuccess = True
                                break
                            
                        if not testsuccess:
                            continue
                        scomplexity = self.codeComplexity(svarfrac[0]) + \
                                      self.codeComplexity(svarfrac[1])
                        ccomplexity = self.codeComplexity(cvarfrac[0]) + \
                                      self.codeComplexity(cvarfrac[1])
                        
                        if scomplexity > 1200 or ccomplexity > 1200:
                            log.debug('equation too complex for single variable solution (%d, %d) ' + \
                                      '... (probably wrong?)', scomplexity, ccomplexity)
                            break
                        
                        if scomplexity < 500 and len(str(svarfrac[1])) < 600:
                            # long fractions can take long time to simplify, so we check the length of equation
                            svarfrac[1] = simplify(svarfrac[1])
                            
                        if self.chop(svarfrac[1])== 0:
                            break
                        
                        if ccomplexity < 500 and len(str(cvarfrac[1])) < 600:
                            cvarfrac[1] = simplify(cvarfrac[1])
                            
                        if self.chop(cvarfrac[1])== 0:
                            break
                        # sometimes the returned simplest solution makes really gross approximations

                        # TGN: ensure othersolvedvars is a subset of self.trigvars_subs
                        assert(all([z in self.trigvars_subs for z in othersolvedvars]))
                        
                        svarfracsimp_denom = self.SimplifyTransform(self.trigsimp_new(svarfrac[1]))
                        cvarfracsimp_denom = self.SimplifyTransform(self.trigsimp_new(cvarfrac[1]))
                        # self.SimplifyTransform could help reduce denoms further...
                        denomsequal = False
                        if self.equal(svarfracsimp_denom, cvarfracsimp_denom):
                            denomsequal = True
                        elif self.equal(svarfracsimp_denom, -cvarfracsimp_denom):
                            cvarfrac[0] = -cvarfrac[0]
                            cvarfracsimp_denom = -cvarfracsimp_denom
                            
                        if self.equal(svarfracsimp_denom,cvarfracsimp_denom) and \
                           not svarfracsimp_denom.is_number:
                            log.debug('denom of %s = %s\n' + \
                                      '        do global subs', \
                                      var.name, svarfracsimp_denom)
                            #denom = self.gsymbolgen.next()
                            #solversolution.dictequations.append((denom,sign(svarfracsimp_denom)))

                            # TGN: ensure othersolvedvars is a subset of self.trigvars_subs
                            assert(all([z in self.trigvars_subs for z in othersolvedvars]))
                            
                            svarsolsimp = self.SimplifyTransform(self.trigsimp_new(svarfrac[0]))#*denom)
                            cvarsolsimp = self.SimplifyTransform(self.trigsimp_new(cvarfrac[0]))#*denom)
                            solversolution.FeasibleIsZeros = False
                            solversolution.presetcheckforzeros.append(svarfracsimp_denom)
                            # instead of doing atan2(sign(dummy)*s, sign(dummy)*c)
                            # we do atan2(s,c) + pi/2*(1-1/sign(dummy)) so equations become simpler
                            #
                            # TGN: or just 1-sign(dummy)?
                            #
                            expandedsol = atan2(svarsolsimp,cvarsolsimp) + pi/2*(-S.One + sign(svarfracsimp_denom))
                        else:
                            
                            # TGN: ensure othersolvedvars is a subset of self.trigvars_subs
                            assert(all([z in self.trigvars_subs for z in othersolvedvars]))
                            
                            svarfracsimp_num = self.SimplifyTransform(self.trigsimp_new(svarfrac[0]))
                            cvarfracsimp_num = self.SimplifyTransform(self.trigsimp_new(cvarfrac[0]))
                            svarsolsimp = svarfracsimp_num/svarfracsimp_denom
                            cvarsolsimp = cvarfracsimp_num/cvarfracsimp_denom
                            
                            if svarsolsimp.is_number and cvarsolsimp.is_number:
                                if Abs(svarsolsimp**2+cvarsolsimp**2-S.One).evalf() > 1e-10:
                                    log.debug('%s solution: atan2(%s, %s), sin/cos not on circle; ignore', \
                                              var.name, svarsolsimp, cvarsolsimp)
                                    continue
                                
                            svarsolsimpcomplexity = self.codeComplexity(svarsolsimp)
                            cvarsolsimpcomplexity = self.codeComplexity(cvarsolsimp)
                            if svarsolsimpcomplexity > 3000 or cvarsolsimpcomplexity > 3000:
                                log.warn('new substituted solutions too complex: %d, %d', \
                                         svarsolsimpcomplexity, \
                                         cvarsolsimpcomplexity)
                                continue
                            
                            try:
                                expandedsol = atan2check(svarsolsimp, cvarsolsimp)
                            except RuntimeError, e:
                                log.warn(u'most likely got recursion error when calling atan2: %s', e)
                                continue
                            
                            solversolution.FeasibleIsZeros = False
                            log.debug('solution for %s: atan2 check for joint', var.name)
                        solversolution.jointeval.append(expandedsol)
                        
                        if unknownvars is not None:
                            unsolvedsymbols = []
                            for unknownvar in unknownvars:
                                if unknownvar != var:
                                    unsolvedsymbols += self.getVariable(unknownvar).vars
                            if len(unsolvedsymbols) > 0:
                                solversolution.equationsused = [eq for eq in eqns \
                                                                if not eq.has(*unsolvedsymbols)]
                            else:
                                solversolution.equationsused = eqns
                                
                            if len(solversolution.equationsused) > 0:
                                log.info('%s = atan2( %s,\n' + \
                                         '                  %s%s )', \
                                         var.name, \
                                         str(solversolution.equationsused[0]),
                                         ' '*len(var.name), 
                                         str(solversolution.equationsused[1]) )
                                
                        if len(self.checkForDivideByZero(expandedsol.subs(solversolution.dictequations))) == 0:
                            goodsolution += 1
                            
                    if len(solversolution.jointeval) == len(sollist) and len(sollist) > 0:
                        solutions.append(solversolution)
                        if goodsolution > 0:
                            hasgoodsolution = True
                        if len(sollist) == goodsolution and goodsolution == 1 and len(solutions) >= 2:
                            break
                        if len(solutions) >= maxsolutions:
                            # probably more than enough already?
                            break

            if len(solutions) > 0 or hasgoodsolution:
                # found a solution without any divides, necessary for pr2 head_torso lookat3d ik
                return solutions

        # solve one equation
        for ieq, eq in enumerate(eqns):
            symbolgen = cse_main.numbered_symbols('const')
            eqnew, symbols = self.groupTerms(eq.subs(varsym.subs), \
                                             [varsym.cvar, varsym.svar, varsym.var], \
                                             symbolgen)
            try:
                # ignore any equations with degree 3 or more 
                ps = Poly(eqnew, varsym.svar)
                pc = Poly(eqnew, varsym.cvar)
                if max(ps.degree_list()) > maxdegree or \
                   max(pc.degree_list()) > maxdegree:
                    log.debug('cannot solve equation with high degree: %s', str(eqnew))
                    continue
                
                if ps.TC() == S.Zero and len(ps.monoms()) > 0:
                    log.debug('%s has trivial solution, ignore', ps)
                    continue
                
                if pc.TC() == S.Zero and len(pc.monoms()) > 0:
                    log.debug('%s has trivial solution, ignore', pc)
                    continue
                
            except PolynomialError:
                # might not be a polynomial, so ignore
                continue

            equationsused = None
            if unknownvars is not None:
                unsolvedsymbols = []
                for unknownvar in unknownvars:
                    if unknownvar != var:
                        unsolvedsymbols += self.getVariable(unknownvar).vars
                if len(unsolvedsymbols) > 0:
                    equationsused = [eq2 for ieq2, eq2 in enumerate(eqns) \
                                     if ieq2 != ieq and not eq2.has(*unsolvedsymbols)]
                else:
                    equationsused = eqns[:]
                    equationsused.pop(ieq)

            numcvar = self.countVariables(eqnew, varsym.cvar)
            numsvar = self.countVariables(eqnew, varsym.svar)
            if numcvar == 1 and numsvar == 1:
                a = Wild('a', exclude = [varsym.svar, varsym.cvar])
                b = Wild('b', exclude = [varsym.svar, varsym.cvar])
                c = Wild('c', exclude = [varsym.svar, varsym.cvar])
                m = eqnew.match(a*varsym.cvar + b*varsym.svar + c)
                if m is not None:
                    symbols += [(varsym.svar, sin(var)), \
                                (varsym.cvar, cos(var))]
                    asinsol = trigsimp(asin(-m[c]/Abs(sqrt(m[a]*m[a]+m[b]*m[b]))).subs(symbols), \
                                       deep = True)
                    # can't use atan2().evalf()... maybe only when m[a] or m[b] is complex?
                    if m[a].has(I) or m[b].has(I):
                        continue
                    constsol = (-atan2(m[a], m[b]).subs(symbols)).evalf()
                    jointsolutions = [constsol + asinsol, \
                                      constsol + pi.evalf() - asinsol]
                    
                    if not constsol.has(I) and \
                       all([self.isValidSolution(s) for s in jointsolutions]) and \
                       len(jointsolutions) > 0:
                        #self.checkForDivideByZero(expandedsol)
                        solutions.append(AST.SolverSolution(var.name, \
                                                            jointeval = jointsolutions, \
                                                            isHinge = self.IsHinge(var.name)))
                        solutions[-1].equationsused = equationsused
                    continue
            if numcvar > 0:
                try:
                    # substitute cos

                    # TGN: the following condition seems weird to me
                    # if  self.countVariables(eqnew, varsym.svar) <= 1 or \
                    #    (self.countVariables(eqnew, varsym.cvar) <= 2 and \
                    #     self.countVariables(eqnew, varsym.svar) == 0):

                    if self.countVariables(eqnew, varsym.svar) <= 1:
                        # anything more than 1 implies quartic equation
                        tempsolutions = solve(eqnew.subs(varsym.svar, sqrt(1-varsym.cvar**2)).expand(), \
                                              varsym.cvar)
                        jointsolutions = []
                        
                        for s in tempsolutions:
                            # TGN: ensure othersolvedvars is a subset of self.trigvars_subs
                            assert(all([z in self.trigvars_subs for z in othersolvedvars]))
                            
                            s2 = self.trigsimp_new(s.subs(symbols+varsym.subsinv))
                            if self.isValidSolution(s2):
                                jointsolutions.append(self.SimplifyTransform(s2))
                        if len(jointsolutions) > 0 and \
                           all([self.isValidSolution(s) \
                                and self.isValidSolution(s) \
                                for s in jointsolutions]):
                            solutions.append(AST.SolverSolution(var.name, \
                                                                jointevalcos = jointsolutions, \
                                                                isHinge = self.IsHinge(var.name)))
                            solutions[-1].equationsused = equationsused
                        continue
                except self.CannotSolveError, e:
                    log.debug(e)
                except NotImplementedError, e:
                    # when solve cannot solve an equation
                    log.warn(e)
            if numsvar > 0:
                # substitute sin
                try:
                    # TGN: the following condition seems weird to me
                    # if  self.countVariables(eqnew, varsym.svar) <= 1 or \
                    #    (self.countVariables(eqnew, varsym.svar) <= 2 and \
                    #     self.countVariables(eqnew, varsym.cvar) == 0):
                    if  self.countVariables(eqnew, varsym.svar) <= 1 or \
                       (self.countVariables(eqnew, varsym.svar) == 2 and \
                        self.countVariables(eqnew, varsym.cvar) == 0):
                        # anything more than 1 implies quartic equation
                        tempsolutions = solve(eqnew.subs(varsym.cvar, \
                                                         sqrt(1-varsym.svar**2)).expand(), \
                                              varsym.svar)

                        # TGN: ensure othersolvedvars is a subset of self.trigvars_subs
                        assert(all([z in self.trigvars_subs for z in othersolvedvars]))
                        
                        jointsolutions = [self.SimplifyTransform(self.trigsimp_new(s.subs(symbols+varsym.subsinv), \
                                                                               )) \
                                          for s in tempsolutions]
                        
                        if all([self.isValidSolution(s) for s in jointsolutions]) and \
                           len(jointsolutions) > 0:
                            solutions.append(AST.SolverSolution(var.name,
                                                                jointevalsin = jointsolutions, \
                                                                isHinge = self.IsHinge(var.name)))
                            solutions[-1].equationsused = equationsused
                        continue
                    
                except self.CannotSolveError, e:
                    log.debug(e)
                    
                except NotImplementedError, e:
                    # when solve cannot solve an equation
                    log.warn(e)

            if numcvar == 0 and numsvar == 0:
                try:
                    tempsolutions = solve(eqnew, var)
                    jointsolutions = []
                    for s in tempsolutions:
                        eqsub = s.subs(symbols)
                        if self.codeComplexity(eqsub) < 2000:
                            
                            # TGN: ensure othersolvedvars is a subset of self.trigvars_subs
                            assert(all([z not in self.trigvars_subs for z in othersolvedvars]))
                            
                            eqsub = self.SimplifyTransform(self.trigsimp_new(eqsub))
                        jointsolutions.append(eqsub)
                        
                    if all([self.isValidSolution(s) and s != S.Zero \
                            for s in jointsolutions]) and \
                                len(jointsolutions) > 0:
                        solutions.append(AST.SolverSolution(var.name, \
                                                            jointeval = jointsolutions, \
                                                            isHinge = self.IsHinge(var.name)))
                        solutions[-1].equationsused = equationsused
                        
                except NotImplementedError, e:
                    # when solve cannot solve an equation
                    log.warn(e)
                continue
            
            try:
                solution = self.solveHighDegreeEquationsHalfAngle([eqnew], varsym, symbols)
                solutions.append(solution.subs(symbols))
                solutions[-1].equationsused = equationsused
            except self.CannotSolveError, e:
                log.debug(e)
                
        if len(solutions) > 0:                
            return solutions
        
        return [self.solveHighDegreeEquationsHalfAngle(eqns, varsym)]

    def SolvePrismaticHingePairVariables(self, raweqns, var0,var1,othersolvedvars,unknownvars=None):
        """solves one hinge and one prismatic variable together
        """
        if self.IsPrismatic(var0.name) and self.IsHinge(var1.name):
            prismaticSymbol = var0
            hingeSymbol = var1
        elif self.IsHinge(var0.name) and self.IsPrismatic(var1.name):
            hingeSymbol = var0
            prismaticSymbol = var1
        else:
            raise self.CannotSolveError('need to have one hinge and one prismatic variable')
        
        prismaticVariable = self.getVariable(prismaticSymbol)
        hingeVariable = self.getVariable(hingeSymbol)
        chingeSymbol,shingeSymbol = hingeVariable.cvar, hingeVariable.svar
        varsubs=prismaticVariable.subs+hingeVariable.subs
        varsubsinv = prismaticVariable.subsinv+hingeVariable.subsinv
        unknownvars=[chingeSymbol,shingeSymbol,prismaticSymbol]
        reducesubs = [(shingeSymbol**2,1-chingeSymbol**2)]
        polyeqs = [Poly(eq.subs(varsubs).subs(reducesubs).expand(),unknownvars) for eq in raweqns if eq.has(prismaticSymbol,hingeSymbol)]
        if len(polyeqs) <= 1:
            raise self.CannotSolveError('not enough equations')
        
        # try to solve one variable in terms of the others
        solvevariables = []
        for polyeq in polyeqs:
            if polyeq.degree(0) == 1 and polyeq.degree(1) == 0:
                chingeSolutions = solve(polyeq,chingeSymbol)
                solvevariables.append((prismaticSymbol,[(chingeSymbol,chingeSolutions[0])]))
            elif polyeq.degree(0) == 0 and polyeq.degree(1) == 1:
                shingeSolutions = solve(polyeq,shingeSymbol)
                solvevariables.append((prismaticSymbol,[(shingeSymbol,shingeSolutions[0])]))
            elif polyeq.degree(2) == 1:
                prismaticSolutions = solve(polyeq,prismaticSymbol)
                solvevariables.append((hingeSymbol,[(prismaticSymbol,prismaticSolutions[0])]))
        
        # prioritize solving the hingeSymbol out
        for solveSymbol in [hingeSymbol,prismaticSymbol]:
            for solveSymbol2, solvesubs in solvevariables:
                if solveSymbol == solveSymbol2:
                    # have a solution for one variable, so substitute it in and see if the equations become solvable with one variable
                    reducedeqs = []
                    for polyeq2 in polyeqs:
                        eqnew = simplify(polyeq2.as_expr().subs(solvesubs))
                        if eqnew != S.Zero:
                            reducedeqs.append(eqnew)
                    self.sortComplexity(reducedeqs)
                    try:
                        rawsolutions = self.solveSingleVariable(reducedeqs,solveSymbol,othersolvedvars, unknownvars=unknownvars)
                        if len(rawsolutions) > 0:
                            return rawsolutions

                    except self.CannotSolveError:
                        pass
                
        raise self.CannotSolveError(u'SolvePrismaticHingePairVariables: failed to find variable with degree 1')
        
    def SolvePairVariables(self, raweqns, var0, var1, \
                           othersolvedvars, \
                           maxcomplexity = 50, \
                           unknownvars = None):
        """
        Solves two hinge variables together.

        Called by SolveAllEquations only.
        """
        # make sure both variables are hinges
        if not (self.IsHinge(var0.name) and self.IsHinge(var1.name)):
            raise self.CannotSolveError('pairwise variables only supports hinge joints')
        
        varsym0 = self.getVariable(var0)
        varsym1 = self.getVariable(var1)
        cvar0,svar0 = varsym0.cvar, varsym0.svar
        cvar1,svar1 = varsym1.cvar, varsym1.svar
        varsubs=varsym0.subs+varsym1.subs
        varsubsinv = varsym0.subsinv+varsym1.subsinv
        unknownvars=[cvar0,svar0,cvar1,svar1]
        reducesubs = [(svar0**2,1-cvar0**2),(svar1**2,1-cvar1**2)]
        eqns = [eq.subs(varsubs).subs(reducesubs).expand() for eq in raweqns if eq.has(var0,var1)]
        if len(eqns) <= 1:
            raise self.CannotSolveError('not enough equation')
        
        # group equations with single variables
        symbolgen = cse_main.numbered_symbols('const')
        orgeqns = []
        allsymbols = []
        for eq in eqns:
            eqnew, symbols = self.groupTerms(eq, unknownvars, symbolgen)
            allsymbols += symbols
            orgeqns.append([self.codeComplexity(eq),Poly(eqnew,*unknownvars)])
        orgeqns.sort(lambda x, y: x[0]-y[0])
        neweqns = orgeqns[:]
        
        pairwisesubs = [(svar0*cvar1,Symbol('s0c1')),(svar0*svar1,Symbol('s0s1')),(cvar0*cvar1,Symbol('c0c1')),(cvar0*svar1,Symbol('c0s1')),(cvar0*svar0,Symbol('s0c0')),(cvar1*svar1,Symbol('c1s1'))]
        pairwiseinvsubs = [(f[1],f[0]) for f in pairwisesubs]
        pairwisevars = [f[1] for f in pairwisesubs]
        reduceeqns = [Poly(eq.as_expr().subs(pairwisesubs),*pairwisevars) for rank,eq in orgeqns if rank < 4*maxcomplexity]
        for i,eq in enumerate(reduceeqns):
            if eq.TC != S.Zero and not eq.TC().is_Symbol:
                n=symbolgen.next()
                allsymbols.append((n,eq.TC().subs(allsymbols)))
                reduceeqns[i] += n-eq.TC()
        
        # try to at least subtract as much paired variables out
        eqcombs = [c for c in combinations(reduceeqns,2)]
        while len(eqcombs) > 0 and len(neweqns) < 20:
            eq0,eq1 = eqcombs.pop()
            eq0dict = eq0.as_dict()
            eq1dict = eq1.as_dict()
            for i in range(6):
                monom = [0,0,0,0,0,0]
                monom[i] = 1
                eq0value = eq0dict.get(tuple(monom),S.Zero)
                eq1value = eq1dict.get(tuple(monom),S.Zero)
                if eq0value != 0 and eq1value != 0:
                    tempeq = (eq0.as_expr()*eq1value-eq0value*eq1.as_expr()).subs(allsymbols+pairwiseinvsubs).expand()
                    if self.codeComplexity(tempeq) > 200:
                        continue
                    eq = simplify(tempeq)
                    if eq == S.Zero:
                        continue
                    
                    peq = Poly(eq,*pairwisevars)
                    if max(peq.degree_list()) > 0 and self.codeComplexity(eq) > maxcomplexity:
                        # don't need such complex equations
                        continue
                    
                    if not self.CheckExpressionUnique(eqns,eq):
                        continue
                    
                    if eq.has(*unknownvars): # be a little strict about new candidates
                        eqns.append(eq)
                        eqnew, symbols = self.groupTerms(eq, unknownvars, symbolgen)
                        allsymbols += symbols
                        neweqns.append([self.codeComplexity(eq),Poly(eqnew,*unknownvars)])

        orgeqns = neweqns[:]
        # try to solve for all pairwise variables
        systemofequations = []
        for i in range(len(reduceeqns)):
            if reduceeqns[i].has(pairwisevars[4],pairwisevars[5]):
                continue
            if not all([__builtin__.sum(m) <= 1 for m in reduceeqns[i].monoms()]):
                continue
            arr = [S.Zero]*5
            for m,c in reduceeqns[i].terms():
                if __builtin__.sum(m) == 1:
                    arr[list(m).index(1)] = c
                else:
                    arr[4] = c
            systemofequations.append(arr)

        if len(systemofequations) >= 4:
            singleeqs = None
            for eqs in combinations(systemofequations,4):
                M = zeros((4,4))
                B = zeros((4,1))
                for i,arr in enumerate(eqs):
                    for j in range(4):
                        M[i,j] = arr[j]
                    B[i] = -arr[4]
                det = self.det_bareis(M,*(self.pvars+unknownvars)).subs(allsymbols)
                if det.evalf() != S.Zero:
                    X = M.adjugate()*B
                    singleeqs = []
                    for i in range(4):
                        eq = (pairwisesubs[i][0]*det - X[i]).subs(allsymbols)
                        eqnew, symbols = self.groupTerms(eq, unknownvars, symbolgen)
                        allsymbols += symbols
                        singleeqs.append([self.codeComplexity(eq),Poly(eqnew,*unknownvars)])
                    break
            if singleeqs is not None:
                neweqns += singleeqs
                neweqns.sort(lambda x, y: x[0]-y[0])

        # check if any equations are at least degree 1 (if not, try to compute some)
        for ivar in range(2):
            polyunknown = []
            for rank,eq in orgeqns:
                p = Poly(eq,unknownvars[2*ivar],unknownvars[2*ivar+1])
                if sum(p.degree_list()) == 1 and __builtin__.sum(p.LM()) == 1:
                    polyunknown.append((rank,p))
            if len(polyunknown) > 0:
                break
        if len(polyunknown) == 0:
            addedeqs = eqns[:]
            polyeqs = []
            for ivar in range(2):
                polyunknown = []
                for rank,eq in orgeqns:
                    p = Poly(eq,unknownvars[2*ivar],unknownvars[2*ivar+1])
                    polyunknown.append(Poly(p.subs(unknownvars[2*ivar+1]**2,1-unknownvars[2*ivar]**2),unknownvars[2*ivar],unknownvars[2*ivar+1]))
                if len(polyunknown) >= 2:
                    monomtoremove = [[polyunknown,(2,0)],[polyunknown,(1,1)]]
                    for curiter in range(2):
                        # remove the square
                        polyunknown,monom = monomtoremove[curiter]
                        pbase = [p for p in polyunknown if p.as_dict().get(monom,S.Zero) != S.Zero]
                        if len(pbase) == 0:
                            continue
                        pbase = pbase[0]
                        pbasedict = pbase.as_dict()
                        for i in range(len(polyunknown)):
                            eq = (polyunknown[i]*pbasedict.get(monom,S.Zero)-pbase*polyunknown[i].as_dict().get(monom,S.Zero)).as_expr().subs(allsymbols)
                            if self.codeComplexity(eq) > 4000:
                                # .. way too complex
                                continue
                            eq = eq.expand()
                            if self.codeComplexity(eq) > 10000:
                                # .. way too complex
                                continue
                            if len(addedeqs) > 10 and self.codeComplexity(eq) > 2000:
                                # .. already have enough...
                                continue
                            if eq != S.Zero and self.CheckExpressionUnique(addedeqs,eq):
                                eqnew, symbols = self.groupTerms(eq, unknownvars, symbolgen)
                                allsymbols += symbols
                                p = Poly(eqnew,*pbase.gens)
                                if p.as_dict().get((1,1),S.Zero) != S.Zero and curiter == 0:
                                    monomtoremove[1][0].insert(0,p)
                                polyeqs.append([self.codeComplexity(eqnew),Poly(eqnew,*unknownvars)])
                                addedeqs.append(eq)
            neweqns += polyeqs
        neweqns.sort(lambda x,y: x[0]-y[0])

        rawsolutions = []
        # try single variable solution, only return if a single solution has been found
        # returning multiple solutions when only one exists can lead to wrong results.
        try:
            rawsolutions += self.solveSingleVariable(self.sortComplexity([e.as_expr().subs(varsubsinv).expand() for score,e in neweqns if not e.has(cvar1,svar1,var1)]),var0,othersolvedvars,subs=allsymbols,unknownvars=unknownvars)
        except self.CannotSolveError:
            pass

        try:
            rawsolutions += self.solveSingleVariable(self.sortComplexity([e.as_expr().subs(varsubsinv).expand() for score,e in neweqns if not e.has(cvar0,svar0,var0)]),var1,othersolvedvars,subs=allsymbols,unknownvars=unknownvars)                    
        except self.CannotSolveError:
            pass

        if len(rawsolutions) > 0:
            solutions = []
            for s in rawsolutions:
                try:
                    solutions.append(s.subs(allsymbols))
                except self.CannotSolveError:
                    pass
                
            if len(solutions) > 0:
                return solutions
        
        groups=[]
        for i,unknownvar in enumerate(unknownvars):
            listeqs = []
            listeqscmp = []
            for rank,eq in neweqns:
                # if variable ever appears, it should be alone
                if all([m[i] == 0 or (__builtin__.sum(m) == m[i] and m[i]>0) for m in eq.monoms()]) and any([m[i] > 0 for m in eq.monoms()]):
                    # make sure there's only one monom that includes other variables
                    othervars = [__builtin__.sum(m) - m[i] > 0 for m in eq.monoms()]
                    if __builtin__.sum(othervars) <= 1:
                        eqcmp = self.removecommonexprs(eq.subs(allsymbols).as_expr(), \
                                                       onlygcd = True, \
                                                       onlynumbers = False)
                        if self.CheckExpressionUnique(listeqscmp,eqcmp):
                            listeqs.append(eq)
                            listeqscmp.append(eqcmp)
            groups.append(listeqs)
        # find a group that has two or more equations:
        useconic=False
        goodgroup = [(i,g) for i,g in enumerate(groups) if len(g) >= 2]
        if len(goodgroup) == 0:
            # might have a set of equations that can be solved with conics
            # look for equations where the variable and its complement are alone
            groups=[]
            for i in [0,2]:
                unknownvar = unknownvars[i]
                complementvar = unknownvars[i+1]
                listeqs = []
                listeqscmp = []
                for rank,eq in neweqns:
                    # if variable ever appears, it should be alone
                    addeq = False
                    if all([__builtin__.sum(m) == m[i]+m[i+1] for m in eq.monoms()]):
                        addeq = True
                    else:
                        # make sure there's only one monom that includes other variables
                        othervars = 0
                        for m in eq.monoms():
                            if __builtin__.sum(m) >  m[i]+m[i+1]:
                                if m[i] == 0 and m[i+1]==0:
                                    othervars += 1
                                else:
                                    othervars = 10000
                        if othervars <= 1:
                            addeq = True
                    if addeq:
                        eqcmp = self.removecommonexprs(eq.subs(allsymbols).as_expr(), \
                                                       onlygcd = True, \
                                                       onlynumbers = False)
                        if self.CheckExpressionUnique(listeqscmp,eqcmp):
                            listeqs.append(eq)
                            listeqscmp.append(eqcmp)
                groups.append(listeqs)
                groups.append([]) # necessary to get indices correct
            goodgroup = [(i,g) for i,g in enumerate(groups) if len(g) >= 2]
            useconic=True
            if len(goodgroup) == 0:
                try:
                    return self.SolvePairVariablesHalfAngle(raweqns,var0,var1,othersolvedvars)
                except self.CannotSolveError,e:
                    log.warn('%s',e)

                # try a separate approach where the two variables are divided on both sides
                neweqs = []
                for rank,eq in neweqns:
                    p = Poly(eq,unknownvars[0],unknownvars[1])
                    iscoupled = False
                    for m,c in p.terms():
                        if __builtin__.sum(m) > 0:
                            if c.has(unknownvars[2],unknownvars[3]):
                                iscoupled = True
                                break
                    if not iscoupled:
                        neweqs.append([p-p.TC(),Poly(-p.TC(),unknownvars[2],unknownvars[3])])
                if len(neweqs) > 0:
                    for ivar in range(2):
                        lineareqs = [eq for eq in neweqs if __builtin__.sum(eq[ivar].LM())==1]
                        for paireq0,paireq1 in combinations(lineareqs,2):
                            log.info('solving separated equations with linear terms')
                            eq0 = paireq0[ivar]
                            eq0dict = eq0.as_dict()
                            eq1 = paireq1[ivar]
                            eq1dict = eq1.as_dict()
                            disc = (eq0dict.get((1,0),S.Zero)*eq1dict.get((0,1),S.Zero) - eq0dict.get((0,1),S.Zero)*eq1dict.get((1,0),S.Zero)).subs(allsymbols).expand()
                            if disc == S.Zero:
                                continue
                            othereq0 = paireq0[1-ivar].as_expr() - eq0.TC()
                            othereq1 = paireq1[1-ivar].as_expr() - eq1.TC()
                            csol = - eq1dict.get((0,1),S.Zero) * othereq0 + eq0dict.get((0,1),S.Zero) * othereq1
                            ssol = eq1dict.get((1,0),S.Zero) * othereq0 - eq0dict.get((1,0),S.Zero) * othereq1
                            polysymbols = paireq0[1-ivar].gens
                            totaleq = (csol**2+ssol**2-disc**2).subs(allsymbols).expand()
                            if self.codeComplexity(totaleq) < 4000:
                                log.info('simplifying final equation to %d', self.codeComplexity(totaleq))
                                totaleq = simplify(totaleq)
                            ptotal_cos = Poly(Poly(totaleq,*polysymbols).subs(polysymbols[0]**2,1-polysymbols[1]**2).subs(polysymbols[1]**2,1-polysymbols[0]**2),*polysymbols)
                            ptotal_sin = Poly(S.Zero,*polysymbols)
                            for m,c in ptotal_cos.terms():
                                if m[1] > 0:
                                    assert(m[1] == 1)
                                    ptotal_sin = ptotal_sin.sub(Poly.from_dict({(m[0],0):c},*ptotal_sin.gens))
                                    ptotal_cos = ptotal_cos.sub(Poly.from_dict({m:c},*ptotal_cos.gens))

                            ptotalcomplexity = self.codeComplexity(ptotal_cos.as_expr()) + self.codeComplexity(ptotal_sin.as_expr())
                            if ptotalcomplexity < 50000:
                                #log.info('ptotal complexity is %d', ptotalcomplexity)
                                finaleq = (ptotal_cos.as_expr()**2 - (1-polysymbols[0]**2)*ptotal_sin.as_expr()**2).expand()
                                # sometimes denominators can accumulate
                                pfinal = Poly(self.removecommonexprs(finaleq),polysymbols[0])
                                pfinal = self.checkFinalEquation(pfinal)
                                if pfinal is not None:
                                    jointsol = atan2(ptotal_cos.as_expr()/ptotal_sin.as_expr(), polysymbols[0])
                                    var = var1 if ivar == 0 else var0
                                    solution = AST.SolverPolynomialRoots(jointname=var.name,poly=pfinal,jointeval=[jointsol],isHinge=self.IsHinge(var.name))
                                    solution.postcheckforzeros = [ptotal_sin.as_expr()]
                                    solution.postcheckfornonzeros = []
                                    solution.postcheckforrange = []
                                    return [solution]
                                
                # if maxnumeqs is any less, it will miss linearly independent equations
                lineareqs = self.solveSingleVariableLinearly(raweqns,var0,[var1],maxnumeqs=len(raweqns))
                if len(lineareqs) > 0:
                    try:
                        return [self.solveHighDegreeEquationsHalfAngle(lineareqs,varsym1)]
                    except self.CannotSolveError,e:
                        log.warn('%s',e)

                raise self.CannotSolveError('cannot cleanly separate pair equations')

        varindex=goodgroup[0][0]
        var = var0 if varindex < 2 else var1
        varsym = varsym0 if varindex < 2 else varsym1
        unknownvar=unknownvars[goodgroup[0][0]]
        eqs = goodgroup[0][1][0:2]
        simpleterms = []
        complexterms = []
        domagicsquare = False
        for i in range(2):
            if useconic:
                terms=[(c,m) for m,c in eqs[i].terms() if __builtin__.sum(m) - m[varindex] - m[varindex+1] > 0]
            else:
                terms=[(c,m) for m,c in eqs[i].terms() if __builtin__.sum(m) - m[varindex] > 0]
            if len(terms) > 0:
                simpleterms.append(eqs[i].sub(Poly.from_dict({terms[0][1]:terms[0][0]},*eqs[i].gens)).as_expr()/terms[0][0]) # divide by the coeff
                complexterms.append(Poly({terms[0][1]:S.One},*unknownvars).as_expr())
                domagicsquare = True
            else:
                simpleterms.append(eqs[i].as_expr())
                complexterms.append(S.Zero)
        finaleq = None
        checkforzeros = []
        if domagicsquare:

            # TGN: ensure othersolvedvars+[var0,var1] is a subset of self.trigvars_subs
            assert(all([z in self.trigvars_subs for z in othersolvedvars+[var0,var1]]))

            # here is the magic transformation:
            finaleq = self.trigsimp_new(expand(((complexterms[0]**2+complexterms[1]**2) \
                                                - simpleterms[0]**2 - simpleterms[1]**2).subs(varsubsinv))).subs(varsubs)
            
            denoms = [fraction(simpleterms[0])[1], \
                      fraction(simpleterms[1])[1], \
                      fraction(complexterms[0])[1], \
                      fraction(complexterms[1])[1]]
            
            lcmvars = self.pvars+unknownvars
            for othersolvedvar in othersolvedvars:
                lcmvars += self.getVariable(othersolvedvar).vars
            denomlcm = Poly(S.One,*lcmvars)
            for denom in denoms:
                if denom != S.One:
                    checkforzeros.append(self.removecommonexprs(denom))
                    denomlcm = Poly(lcm(denomlcm,denom),*lcmvars)
            finaleq = simplify(finaleq*denomlcm.as_expr()**2)
            complementvarindex = varindex-(varindex%2)+((varindex+1)%2)
            complementvar = unknownvars[complementvarindex]
            finaleq = simplify(finaleq.subs(complementvar**2,1-unknownvar**2)).subs(allsymbols).expand()
        else:
            # try to reduce finaleq
            p0 = Poly(simpleterms[0],unknownvars[varindex],unknownvars[varindex+1])
            p1 = Poly(simpleterms[1],unknownvars[varindex],unknownvars[varindex+1])
            if max(p0.degree_list()) > 1 \
               and max(p1.degree_list()) > 1 \
               and max(p0.degree_list()) == max(p1.degree_list()) \
               and p0.LM() == p1.LM():
                finaleq = (p0*p1.LC()-p1*p0.LC()).as_expr()
                finaleq = expand(simplify(finaleq.subs(allsymbols)))
                if finaleq == S.Zero:
                    finaleq = expand(p0.as_expr().subs(allsymbols))
        if finaleq is None:
            log.warn('SolvePairVariables: did not compute a final variable. This is a weird condition...')
            return self.SolvePairVariablesHalfAngle(raweqns,var0,var1,othersolvedvars)
        
        if not self.isValidSolution(finaleq):
            log.warn('failed to solve pairwise equation: %s'%str(finaleq))
            return self.SolvePairVariablesHalfAngle(raweqns,var0,var1,othersolvedvars)

        newunknownvars = unknownvars[:]
        newunknownvars.remove(unknownvar)
        if finaleq.has(*newunknownvars):
            log.warn('equation relies on unsolved variables(%s):\n' + \
                     '        %s',newunknownvars, finaleq)
            return self.SolvePairVariablesHalfAngle(raweqns,var0,var1,othersolvedvars)

        if not finaleq.has(unknownvar):
            # somehow removed all variables, so try the general method
            return self.SolvePairVariablesHalfAngle(raweqns,var0,var1,othersolvedvars)

        try:
            if self.codeComplexity(finaleq) > 100000:
                return self.SolvePairVariablesHalfAngle(raweqns,var0,var1,othersolvedvars)
            
        except self.CannotSolveError:
            pass

        if useconic:
            # conic roots solver not as robust as half-angle transform!
            #return [SolverConicRoots(var.name,[finaleq],isHinge=self.IsHinge(var.name))]
            solution = self.solveHighDegreeEquationsHalfAngle([finaleq],varsym)
            solution.checkforzeros += checkforzeros
            return [solution]

        # now that everything is with respect to one variable, simplify and solve the equation
        eqnew, symbols = self.groupTerms(finaleq, unknownvars, symbolgen)
        allsymbols += symbols
        solutions = solve(eqnew,unknownvar)
        log.info('pair solution: %s, %s', eqnew,solutions)
        if solutions:
            solversolution = AST.SolverSolution(var.name, isHinge=self.IsHinge(var.name))
            processedsolutions = []
            for s in solutions:
                processedsolution = s.subs(allsymbols+varsubsinv).subs(varsubs)
                # trigsimp probably won't work on long solutions
                if self.codeComplexity(processedsolution) < 2000: # complexity of 2032 for pi robot freezes
                    log.info('solution complexity: %d', self.codeComplexity(processedsolution))
                    processedsolution = self.SimplifyTransform(self.trigsimp(processedsolution,othersolvedvars))
                processedsolutions.append(processedsolution.subs(varsubs))
            if (varindex%2)==0:
                solversolution.jointevalcos=processedsolutions
            else:
                solversolution.jointevalsin=processedsolutions
            return [solversolution]
        
        return self.SolvePairVariablesHalfAngle(raweqns,var0,var1,othersolvedvars)
        #raise self.CannotSolveError('cannot solve pair equation')
        

