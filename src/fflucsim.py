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
        self.dead = False
        self.history = []
 
        if mother != None:
            self.w_mono = mother.w_mono
            self.w_triso = mother.w_triso
            self.mu_A = mother.mu_A
            self.mu_R = mother.mu_R
            self.mu_T = mother.mu_T
            self.ploidy = mother.ploidy
            self.genealogy = f'{mother.genealogy}.1'
            self.mother = mother.uid
            self.disome = mother.disome
            self.monosome = mother.monosome
            self.revertant = mother.revertant
            self.revertant2 = mother.revertant2
            self.trisome = mother.trisome
            self.tetrasome = mother.tetrasome
            self.generation_monosome = mother.generation_monosome
            self.generation_revert = mother.generation_revert
            self.generation_trisome = mother.generation_trisome
            self.generation_revert2 = mother.generation_revert2
            self.generation_tetrasome = mother.generation_tetrasome
            self.rng = mother.rng
            self.homologs = {}
            for (h,i) in mother.homologs.items():
                self.homologs[h] = i
        
        else:
            self.ploidy = 2
            self.genealogy = '0'
            self.mother = None
            self.disome = True
            self.monosome = False
            self.revertant = False
            self.revertant2 = False
            self.trisome = False
            self.tetrasome = False
            self.generation_monosome = None
            self.generation_revert = None
            self.generation_trisome = None
            self.generation_revert2 = None
            self.generation_tetrasome = None
            self.rng = np.random.default_rng()
            self.homologs = {0:1, 1:1}

    def founder(self, w_mono, w_triso, mu_A, mu_R, mu_T):
        self.w_mono = w_mono
        self.w_triso = w_triso
        self.mu_A = mu_A
        self.mu_R = mu_R
        self.mu_T = mu_T
        
        return self
    
    def get_homologs(self):
        return np.concatenate([np.repeat(h,i) for (h,i) in self.homologs.items()])

    def selection(self):
        return self.get_homologs().sum == 0
        
    def become_monosome(self, gen, h=None):
        self.disome = False
        self.monosome = True
        self.ploidy -= 1
        self.generation_monosome = gen
        if h == None:
            h = self.rng.choice(self.get_homologs(), 1)[0]
        self.homologs[h] -= 1
        return h
        
    def become_trisome(self, gen, h=None):
        self.disome = False
        self.trisome = True
        self.ploidy += 1
        self.generation_trisome = gen
        if h == None:
            h = self.rng.choice(self.get_homologs(), 1)[0]
        self.homologs[h] += 1
        return h
        
    def become_revert(self, gen, h=None):
        self.monosome = False
        self.revertant = True
        self.disome = True
        self.ploidy += 1
        self.generation_revert = gen
        if h == None:
            h = self.rng.choice(self.get_homologs(), 1)[0]
        self.homologs[h] += 1
        return h
    
    def become_revert2(self, gen, h=None):
        self.trisome = False
        self.revertant2 = True
        self.disome = True
        self.ploidy -= 1
        self.generation_revert2 = gen
        if h == None:
            h = self.rng.choice(self.get_homologs(), 1)[0]
        self.homologs[h] -= 1
        return h

    def become_tetrasome(self, gen, h):
        self.trisome = False
        self.tetrasome = True
        self.revertant2 = False
        self.ploidy += 1
        self.generation_tetrasome = gen
        self.homologs[h] += 1

    def decide_divide(self):
        if self.dead == False:
            if self.monosome:
                return self.rng.random() < self.w_mono
            elif self.trisome or self.tetrasome:
                return self.rng.random() < self.w_triso
            else:
                return True
        else:
            return False

    def decide_monosome(self):
        return self.rng.random() < self.mu_A

    def decide_revert(self):
        return self.rng.random() < self.mu_R

    def decide_revert2(self):
        return self.rng.random() < self.mu_T

    def divide(self, daughter_uid, gen):
        daughter = CellMonosome(daughter_uid, gen, self)

        event_monosome = False
        event_revert = False
        event_revert2 = False
    
        self.age += 1
        self.daughters.append(daughter_uid)
        self.genealogy = f'{self.genealogy}.0'
        
        if self.disome:
            if self.decide_monosome():
                h = daughter.become_monosome(gen)
                self.become_trisome(gen, h)
                event_monosome = True

        elif self.monosome:
            if self.decide_revert():
                h = daughter.become_revert(gen)
                self.dead = True
                event_revert = True
        
        elif self.trisome:
            if self.decide_revert2():
                h = daughter.become_revert2(gen)
                self.become_tetrasome(gen, h)
                event_revert2 = True


        self.history.append({'gen':gen, 
                            'homologs':self.get_homologs(), 
                            'event_monosome':event_monosome, 
                            'event_revert':event_revert, 
                            'event_revert2':event_revert2})
        return daughter, event_monosome, event_revert, event_revert2
        
    def summary_print(self):
        summary_table = ['--------CELL SUMMARY---------',
                        f'| ploidy:     | {str(self.ploidy).ljust(12, " ")}|',
                        f'| monosome:   | {str(self.monosome).ljust(12, " ")}|',
                        f'| trisome:    | {str(self.monosome).ljust(12, " ")}|',
                        f'| revertant:  | {str(self.revertant).ljust(12, " ")}|',
                        f'| revertant2: | {str(self.revertant2).ljust(12, " ")}|',
                        f'| age:        | {str(self.age).ljust(12, " ")}|',
                         '--------CELL SUMMARY--------'
                        ]
        for l in summary_table:
            print(l)
    
    def summary(self):
        L = {'uid':self.uid,
            'ploidy':self.ploidy,
            'age':self.age,
            'monosome':self.monosome,
            'trisome':self.trisome,
            'revertant':self.revertant,
            'revertant2':self.revertant2,
            'mother':self.mother,
            'daughters':';'.join([str(d) for d in self.daughters]),
            'genealogy':self.genealogy}
        return L

    def get_binary_flag(self):
        flags = [self.monosome, self.disome, self.revertant, self.revertant2, self.trisome, self.tetrasome]
        return ''.join([str(int(i)) for i in flags])
        
    def select(self):
        return self.homologs[1] > 0

class Population:
    
    def __init__(self, founder, target_div):
        self.founder = founder
        self.target_div = target_div
        self.w_mono = founder.w_mono
        self.w_triso = founder.w_triso
        self.mu_A = founder.mu_A
        self.mu_R = founder.mu_R
        self.mu_T = founder.mu_T
        self.ploidy = founder.ploidy

    def expand(self, verbose=False, cleanup=False, select=True):

        pop = [self.founder]
        founder_flag = self.founder.get_binary_flag()
        len_pop = 1
        uid = self.founder.uid
        gen = self.founder.born
        Events_monosome = []
        Events_revert = []
        Events_revert2 = []
        time_init = time.time()
        
        while len_pop < self.target_div:
            gen += 1
            new_gen = []
            len_new_gen = 0
            for c in pop:
                if c.decide_divide():
                    uid += 1
                    daughter, new_monosome, new_revert, new_revert2 = c.divide(uid, gen)
                    new_gen.append(daughter)
                    len_new_gen += 1
                    
                    if new_monosome:
                        Events_monosome.append((uid, gen))
                    elif new_revert:
                        Events_revert.append((uid, gen))
                    elif new_revert2:
                        Events_revert2.append((uid, gen))
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
        self.Events_revert2 = Events_revert2
        self.compute_time = time_final-time_init
        self.Report = self.report()

        if cleanup:
            self.clean_population(founder_flag)
        
        if select:
            self.select_population()
            self.ReportSelect = self.report()

    def report(self):
        Report = {}
        Report['w_mono'] = self.w_mono
        Report['w_triso'] = self.w_triso
        Report['mu_A'] = self.mu_A
        Report['mu_R'] = self.mu_R
        Report['mu_T'] = self.mu_T
        Report['ploidy'] = self.ploidy
        Report['m_monosome'] = len(self.Events_monosome)
        Report['m_revert'] = len(self.Events_revert)
        Report['m_revert2'] = len(self.Events_revert2)
        bool_cts = np.array([[c.monosome, c.revertant, c.revertant2, c.dead] for c in self.Population.values()]).sum(axis=0)
        Report['n_monosome'] = bool_cts[0]
        Report['n_revert'] = bool_cts[1]
        Report['n_revert2'] = bool_cts[2]
        Report['n_total'] = bool_cts[:3].sum()
        Report['n_dead'] = bool_cts[3]
        Report['final_size'] = len(self.Population)

        return Report

    def clean_population(self, founder_flag):
        for uid, c in list(self.Population.items()):
            if c.get_binary_flag() == founder_flag:
                del self.Population[uid]

    def select_population(self):
        for uid, c in list(self.Population.items()):
            if c.select():
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
            counts = self.n_total
        elif mutant == 'revert':
            counts = self.n_revert
        
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
            counts = self.n_total
        elif mutant == 'revert':
            counts = self.n_revert
            
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

def target_div_monosome_rate(mr, x=20):
    target_div = np.int64(x/mr)
    return min([target_div, 2**25])
