import numpy as np
import time
import pickle as pkl
import pysalvador as sal

class CellMonosome:

    def __init__(self, uid, gen, mother):

        self.uid = uid
        self.born = gen
        self.daughters = []
        self.age = 0
        
        if mother != None:
            self.monosome_fitness = mother.monosome_fitness
            self.monosome_rate = mother.monosome_rate
            self.revert_rate = mother.revert_rate
            self.ploidy = mother.ploidy
            self.genealogy = f'{mother.genealogy}.1'
            self.mother = mother.uid
            self.monosome = mother.monosome
            self.revertant = mother.revertant
            self.generation_monosome = mother.generation_monosome
            self.generation_revert = mother.generation_revert
            self.rng = mother.rng
        
        else:
            self.genealogy = '0'
            self.mother = None
            self.monosome = False
            self.revertant = False
            self.generation_monosome = None
            self.generation_revert = None
            self.rng = np.random.default_rng()

    def founder(self, monosome_fitness, monosome_rate, revert_rate, ploidy):
        self.monosome_fitness = monosome_fitness
        self.monosome_rate = monosome_rate
        self.revert_rate = revert_rate
        self.ploidy = ploidy
        
        return self
        
    def become_monosome(self, gen):
        self.monosome = True
        self.ploidy -= 1
        self.generation_monosome = gen
        
    def become_revert(self, gen):
        self.revertant = True
        self.monosome = False
        self.ploidy += 1
        self.generation_revert = gen

    def decide_divide(self):
        if self.monosome:
            return self.rng.random() < self.monosome_fitness
        else:
            return True

    def decide_monosome(self):
        return self.rng.random() < self.monosome_rate

    def decide_revert(self):
        return self.rng.random() < self.revert_rate

    def divide(self, daughter_uid, gen):
        daughter = CellMonosome(daughter_uid, gen, self)

        event_monosome = False
        event_revert = False
        
        self.age += 1
        self.daughters.append(daughter_uid)
        self.genealogy = f'{self.genealogy}.0'
        
        if not self.monosome:
            if self.decide_monosome():
                daughter.become_monosome(gen)
                event_monosome = True
        else:
            if self.decide_revert():
                daughter.become_revert(gen)
                event_revert = True
            
        return daughter, event_monosome, event_revert
        
    def summary(self):
        summary_table = ['--------CELL SUMMARY--------',
                        f'| ploidy:    | {str(self.ploidy).ljust(12, " ")}|',
                        f'| monosome:  | {str(self.monosome).ljust(12, " ")}|',
                        f'| revertant: | {str(self.revertant).ljust(12, " ")}|',
                        f'| age:       | {str(self.age).ljust(12, " ")}|',
                         '--------CELL SUMMARY--------'
                        ]
        for l in summary_table:
            print(l)

class Population:
    
    def __init__(self, founder, target_div):
        self.founder = founder
        self.target_div = target_div
        self.monosome_fitness = founder.monosome_fitness
        self.monosome_rate = founder.monosome_rate
        self.revert_rate = founder.revert_rate
        self.ploidy = founder.ploidy

    def expand(self, verbose=False, cleanup=True):

        pop = [self.founder]
        len_pop = 1
        uid = self.founder.uid
        gen = self.founder.born
        Events_monosome = []
        Events_revert = []
        time_init = time.time()
        
        #for _ in range(self.n_gen):
        while len_pop < self.target_div:
            gen += 1
            new_gen = []
            len_new_gen = 0
            for c in pop:
                if c.decide_divide():
                    uid += 1
                    daughter, new_monosome, new_revert = c.divide(uid, gen)
                    new_gen.append(daughter)
                    len_new_gen += 1
                    
                    if new_monosome:
                        Events_monosome.append((uid, gen))
                    elif new_revert:
                        Events_revert.append((uid, gen))
                    if len_pop + len_new_gen >= self.target_div:
                        break
                    
            if verbose:
                print(f'END OF GEN {gen}: daughters/mothers: {len_new_gen}/{len_pop} cumul: {len_pop+len_new_gen}')

            pop.extend(new_gen)
            len_pop += len_new_gen
            
        if verbose:
            print(f'final pop size: {len_pop} Events monosome: {len(Events_monosome)} Events revert: {len(Events_revert)}')
        
        Pop = {c.uid:c for c in pop}
        time_final = time.time()
        
        self.Population = Pop
        self.Events_monosome = Events_monosome
        self.Events_revert = Events_revert
        self.compute_time = time_final-time_init
        self.Report = self.report()

        if cleanup:
            self.clean_population()

    def report(self):
        Report = {}
        Report['monosome_fitness'] = self.monosome_fitness
        Report['monosome_rate'] = self.monosome_rate
        Report['revert_rate'] = self.revert_rate
        Report['ploidy'] = self.ploidy
        Report['m_monosome'] = len(self.Events_monosome)
        Report['m_revert'] = len(self.Events_revert)
        Report['n_monosome'] = sum([c.monosome for c in self.Population.values()])
        Report['n_revert'] = sum([c.revertant for c in self.Population.values()])
        Report['n_total'] = sum([c.monosome or c.revertant for c in self.Population.values()])
        Report['final_size'] = len(self.Population)

        return Report

    def clean_population(self):
        for uid, c in list(self.Population.items()):
            if not c.monosome and not c.revertant:
                del self.Population[uid]

def load_populations(list_of_fn):
    Populations = []
    for fn in list_of_fn:
        with open(fn, 'rb') as fi:
            Populations.append(pkl.load(fi))
    return Populations

def FluctuationAssayResult(m, m_ci, Nt, upper_bound, mutant, model, w):
    Results = {}
    
    mu = m/Nt
    m_ci = np.array(m_ci)
    mu_ci = m_ci/Nt
    
    Results['m'] = m
    Results['m_ci'] = m_ci
    Results['mu'] = mu
    Results['mu_ci'] = mu_ci
    Results['upper_bound'] = upper_bound
    Results['mutant'] = mutant
    Results['model'] = model
    Results['w'] = w

    return Results

class FluctuationAssay:
    
    def __init__(self, Populations):

        self.Populations = Populations
        self.replicates = len(Populations)
        self.monosome_fitness = [pop.monosome_fitness for pop in Populations][0]
        self.monosome_rate = [pop.monosome_rate for pop in Populations][0]
        self.revert_rate = [pop.revert_rate for pop in Populations][0]
        self.ploidy = [pop.ploidy for pop in Populations][0]
        
        self.m_monosome = [pop.Report['m_monosome'] for pop in Populations]
        self.m_revert = [pop.Report['m_revert'] for pop in Populations]
        self.n_monosome = [pop.Report['n_monosome'] for pop in Populations]
        self.n_revert = [pop.Report['n_revert'] for pop in Populations]
        self.n_total = [pop.Report['n_total'] for pop in Populations]
        self.final_size = [pop.Report['final_size'] for pop in Populations]

        self.Nt = np.mean(self.final_size)
        self.Results = []

    def fit_LD(self, mutant):

        if mutant == 'monosome':
            counts = self.m_total
        elif mutant == 'revert':
            counts = self.m_revert
        
        upper_bound = False
        if np.all(np.array(counts)==0):
            counts[0] = 1
            upper_bound = True

        try:
            m = sal.newtonLD(counts, max_iter=100)
            m_ci = np.array(sal.confintLD(counts, max_iter=100))
        except:
            m = np.nan
            m_ci = (np.nan, np.nan)
        mu = m/self.Nt
        mu_ci = m_ci/self.Nt
        fres = FluctuationAssayResult(m, m_ci, self.Nt, upper_bound, mutant, 'LD', 1)
        
        self.Results.append(fres)

    def fit_MK(self, mutant, w=None):

        if w == None:
            w = self.monosome_fitness

        if mutant == 'monosome':
            counts = self.m_total
        elif mutant == 'revert':
            counts = self.m_revert
            
        upper_bound = False
        if np.all(np.array(counts)==0):
            counts[0] = 1
            upper_bound = True

        try:
            m = sal.newtonMK(counts, w=w, max_iter=100)
            m_ci = np.array(sal.confintMK(counts, w=w, max_iter=100))
        except:
            m = np.nan
            m_ci = (np.nan, np.nan)
        mu = m/self.Nt
        mu_ci = m_ci/self.Nt
        fres = FluctuationAssayResult(m, m_ci, self.Nt, upper_bound, mutant, 'MK', w)
        
        self.Results.append(fres)

def ngen_monosome_rate(mr, x=20):
    n_gen = np.int64(np.ceil(np.log2(1/mr*x)))
    return min([n_gen, 25])
