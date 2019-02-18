
from pysb.core import MonomerPattern, ComplexPattern, RuleExpression, ReactionPattern, ANY, WILD
from pysb.builder import Builder
from pysb.export.pysb_flat import PysbFlatExporter
from copy import deepcopy
from itertools import combinations, product
import re
import csv
from collections import defaultdict


class M_Node:

    def __init__(self):
        self.incidentNodes = []
        self.labels = []
        self.reactions = []
        self.optional_rxns = []
        self.compete = []
        self.sequence = []
        self.exclusive = False
        self.relax = False
        self.depth = 1000
        self.initial = None
        self.motifs = []
        self.partials = []
        self.self = []
        self.objective = False
        self.boolean = None
        self.boolean_list = None
        self.table = None


def importModel(model_name):

    nodes = {}

    # read in the Boolean model in Booleannet format; this should be separate in the future
    model = open(model_name, 'r')
    bnodes = []
    brules = []
    brules_parced = []
    binputs = []

    within_model = False
    for line in model:

        if '\"\"\"' in line and within_model:
            break

        if within_model:

            if '=' in line and '*' not in line:
                s = line.find(' ')
                nodes[line[:s]] = M_Node()
                if line[line.rfind(' ')+1:-1] == 'False':
                    nodes[line[:s]].initial = 0
                if line[line.rfind(' ')+1:-1] == 'True':
                    nodes[line[:s]].initial = 1

            if '*' in line and ':' in line:
                bnodes.append(line[line.index(':')+2:line.index('*')])
                brules.append('('+line[line.index('=') + 2:-1]+')')
                bp = deepcopy(brules[-1])
                bp = bp.replace('(', ' ( ')
                bp = bp.replace(')', ' ) ')
                bp = bp.strip()
                bp = re.split(r'\s*', bp)
                brules_parced.append(bp)
                binputs.append([])

            if '*' in line and ':' not in line:
                bnodes.append(line[:line.index('*')])
                brules.append(line[line.index('=')+2:-1])
                bp = deepcopy(brules[-1])
                bp = bp.replace('(', ' ( ')
                bp = bp.replace(')', ' ) ')
                bp = bp.strip()
                bp = re.split(r'\s*', bp)
                brules_parced.append(bp)
                binputs.append([])

        if '\"\"\"' in line and not within_model:
            within_model = True

    btable = []

    for i, each in enumerate(brules_parced):
        temp = []
        for item in bnodes:
            if item in each:
                temp.append([each.index(item), item])
        temp.sort()
        for item in temp:
            binputs[i].append(item[1])

    for each in binputs:
        btable.append(list(product([True, False], repeat=len(each))))
    for i, each in enumerate(btable):
        for j, item in enumerate(each):
            btable[i][j] = list(btable[i][j])
    for i, each in enumerate(brules_parced):
        for j, item in enumerate(btable[i]):
            rule = deepcopy(each)
            b_list = deepcopy(item)
            inputs = deepcopy(binputs[i])
            for k, every in enumerate(rule):
                for l, thing in enumerate(inputs):
                    if every == thing:
                        rule[k] = str(b_list[l])
            rule = ' '.join(rule)
            btable[i][j].append(eval(rule))

    for i, each in enumerate(bnodes):
        head = deepcopy(binputs[i])
        head.append(each)
        btable[i].insert(0, head)
    for i, each in enumerate(bnodes):
        nodes[each].boolean = brules[i]
        nodes[each].boolean_list = brules_parced[i]
        nodes[each].table = btable[i]
        for item in binputs[i]:
            nodes[each].incidentNodes.append(item)

    return nodes


class Rxn:

    def __init__(self, mole, base_states, par, react, part, direct, reacts, targs, rxnTemps, instructs, params):
        self.molecule = mole
        self.base_states = base_states
        self.parent = par
        self.reaction = react
        self.partners = part
        self.direction = direct
        self.reactants = reacts
        self.targets = targs
        self.rxnTemplates = rxnTemps
        self.rxnsParsed = []
        self.instructions = instructs
        self.parameters = params


def importLibrary(library_name):

    # read in the library of molecules and their associated reactions
    # todo: add consistency tests
    # todo: make list of reactions a dictionary

    molecule_list = defaultdict(list)

    molecule = None
    base_states = []
    parent = None
    reaction = None
    partners = []
    direction = None
    reactants = []
    targets = []
    rxnTemplate = []
    instructions = []
    parameters = []

    mol_library = open(library_name)
    for line in mol_library:
        if 'molecule:' in line:
            molecule = line.split(':', 1)[1].strip()
        if 'base_state:' in line:
            base_states.append(line.split(':', 1)[1].strip())
        if 'parent:' in line:
            parent = line.split(':', 1)[1].strip()
        if 'reaction:' in line:
            reaction = line.split(':', 1)[1].strip()
        if 'partner:' in line:
            partners.append(line.split(':', 1)[1].strip())
        if 'direction:' in line:
            direction = line.split(':', 1)[1].strip()
        if 'reactant:' in line:
            reactants.append(line.split(':', 1)[1].strip())
        if 'target:' in line:
            targets.append(line.split(':', 1)[1].strip())
        if 'rxnTemplate:' in line:
            rxnTemplate.append(line.split(':', 1)[1].strip())
        if 'instruction:' in line:
            instructions.append(line.split(':', 1)[1].strip())
        if 'parameters:' in line:
            instructions.append(line.split(':', 1)[1].strip())
        if '$$$' in line:
            molecule_list[molecule].append(
                Rxn(molecule, base_states, parent, reaction, partners, direction, reactants, targets, rxnTemplate, instructions, parameters))

            reaction = None
            partners = []
            direction = None
            reactants = []
            targets = []
            rxnTemplate = []
            instructions = []
            parameters = []
        if '+++' in line:
            base_states = []

    return molecule_list


def importLabels(labels, library, nodes):

    label_info = False
    initial = False
    required = False
    optional = False
    exclusive = False
    depth = False
    relax = False
    compete = False
    sequence = False
    objective = False
    self = False
    saturated = False

    with open(labels) as label_file:
        reader = csv.reader(label_file)
        label_list = list(reader)
        for each in label_list:
            if each:
                if each[0].strip() == 'labels':
                    label_info = True
                    initial = False
                    required = False
                    optional = False
                    exclusive = False
                    depth = False
                    relax = False
                    compete = False
                    sequence = False
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'initial values':
                    label_info = False
                    initial = True
                    required = False
                    optional = False
                    exclusive = False
                    depth = False
                    relax = False
                    compete = False
                    sequence = False
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'required reactions':
                    label_info = False
                    initial = False
                    required = True
                    optional = False
                    exclusive = False
                    depth = False
                    relax = False
                    compete = False
                    sequence = False
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'optional reactions':
                    label_info = False
                    initial = False
                    required = False
                    optional = True
                    exclusive = False
                    depth = False
                    relax = False
                    compete = False
                    sequence = False
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'restrict reactions':
                    label_info = False
                    initial = False
                    required = False
                    optional = False
                    exclusive = True
                    depth = False
                    relax = False
                    compete = False
                    sequence = False
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'depth':
                    label_info = False
                    initial = False
                    required = False
                    optional = False
                    exclusive = False
                    depth = True
                    relax = False
                    compete = False
                    sequence = False
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'relax':
                    label_info = False
                    initial = False
                    required = False
                    optional = False
                    exclusive = False
                    depth = False
                    relax = True
                    compete = False
                    sequence = False
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'compete':
                    label_info = False
                    initial = False
                    required = False
                    optional = False
                    exclusive = False
                    depth = False
                    relax = False
                    compete = True
                    sequence = False
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'sequence':
                    label_info = False
                    initial = False
                    required = False
                    optional = False
                    exclusive = False
                    depth = False
                    relax = False
                    compete = False
                    sequence = True
                    objective = False
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'objective':
                    label_info = False
                    initial = False
                    required = False
                    optional = False
                    exclusive = False
                    depth = False
                    relax = False
                    compete = False
                    sequence = False
                    objective = True
                    self = False
                    saturated = False
                    continue
                if each[0].strip() == 'self reactions':
                    label_info = False
                    initial = False
                    required = False
                    optional = False
                    exclusive = False
                    depth = False
                    relax = False
                    compete = False
                    sequence = False
                    objective = False
                    self = True
                    saturated = False
                    continue
                if each[0].strip() == 'saturated':
                    label_info = False
                    initial = False
                    required = False
                    optional = False
                    exclusive = False
                    depth = False
                    relax = False
                    compete = False
                    sequence = False
                    objective = False
                    self = False
                    saturated = True
                    continue

                if label_info:
                    for lab in each[1:]:
                        if lab.strip() != 'all':
                            nodes[each[0].strip()].labels.append(lab.strip())

                        # todo: add exception here
                        if lab.strip() not in library and lab.strip() != 'all':
                            print lab, 'not in library'

                        if lab.strip() == 'all':
                            for item in library:
                                nodes[each[0].strip()].labels.append(item)

                if initial:
                    nodes[each[0].strip()].initial = each[1]
                if required and [each[1].strip(), [x.strip() for x in each[2].split(':')], [x.strip() for x in each[3].split(':')]] not in nodes[each[0].strip()].reactions:
                    nodes[each[0].strip()].reactions.append([each[1].strip(), [x.strip() for x in each[2].split(':')], [x.strip() for x in each[3].split(':')], [x.strip() for x in each[4:]]])
                if optional and [each[1].strip(), [x.strip() for x in each[2].split(':')], [x.strip() for x in each[3].split(':')]] not in nodes[each[0].strip()].optional_rxns:
                    nodes[each[0].strip()].optional_rxns.append([each[1].strip(), [x.strip() for x in each[2].split(':')], [x.strip() for x in each[3].split(':')], [x.strip() for x in each[4:]]])
                    nodes[each[0].strip()].exclusive = True
                if self and [each[1].strip(), [x.strip() for x in each[2].split(':')], [x.strip() for x in each[3].split(':')]] not in nodes[each[0].strip()].self:
                    nodes[each[0].strip()].self.append([each[1].strip(), [x.strip() for x in each[2:]]])
                if exclusive:
                    nodes[each[0].strip()].exclusive = True
                if depth:
                    nodes[each[0].strip()].depth = each[1]
                if relax:
                    nodes[each[0].strip()].relax = True
                if compete:
                    nodes[each[0].strip()].compete.append(each[1].strip())
                if sequence:
                    nodes[each[0].strip()].sequence.append([])
                    for item in each[1:]:
                        nodes[each[0].strip()].sequence[-1].append(item.strip())
                if objective:
                    nodes[each[0].strip()].objective = True
                if saturated:
                    nodes[each[0].strip()].saturated = each[1:]
                    for i, item in enumerate(nodes[each[0].strip()].saturated):
                        nodes[each[0].strip()].saturated[i] = item.strip()


class MotifBuilder:

    def __init__(self, nodes, library, hidden=0, max_depth=100):
        self.nodes = nodes
        self.library = library
        self.hidden = hidden
        self.max_depth = max_depth
        self.count = 0
        self.master = deepcopy(self.nodes)
        self._build_motifs()

    def _findSelfInteractions(self, n):

        interacts = []
        for this in self.master[n].labels:
            for that in self.library[this]:
                if that.direction == 'self':
                    interacts.append([[None], [n], [None], [that.molecule], that.reaction])

        return interacts

    def _findInteractions(self, target_node, affecting_nodes):

        targ_list = []
        aff_list = []
        interacts = []

        # build lists of reactions for target 'input', affecting 'output', target 'mutual', and affecting 'mutual' nodes
        for target_label in self.master[target_node].labels:
            for mol in self.library[target_label]:
                if mol.direction == 'input' and\
                        [target_node, mol.reaction, mol.molecule, mol.reactants, 'input'] not in targ_list:
                    targ_list.append([target_node, mol.reaction, mol.molecule, mol.reactants, 'input'])

        for aff_node in affecting_nodes:
            for aff_label in self.master[aff_node].labels:
                for mol in self.library[aff_label]:
                    if mol.direction == 'output' and\
                            [aff_node, mol.reaction, mol.molecule, mol.targets, 'output'] not in aff_list:
                        aff_list.append([aff_node, mol.reaction, mol.molecule, mol.targets, 'output'])

        # group reactions - target 'input', affecting 'output'
        targ_groups = []
        aff_groups = []

        aff_max = 0
        for rxn in aff_list:
            if len(rxn[3]) > aff_max:
                aff_max = len(rxn[3])
        for r in range(aff_max):
            for combos in combinations(targ_list, r+1):
                targ_groups.append(list(combos))

        targ_max = 0
        for rxn in targ_list:
            if len(rxn[3]) > targ_max:
                targ_max = len(rxn[3])
        for r in range(targ_max):
            for combos in combinations(aff_list, r+1):
                aff_groups.append(list(combos))

        # match input and output reactions

        for t_group in targ_groups:

            t_group_o = deepcopy(t_group)

            reaction_name = t_group[0][1]
            t2 = []
            t3 = []

            for t in t_group:
                t2.append(t[2])
                t3.append(t[3])

            for i, t in enumerate(t3):
                for j, elem in enumerate(t):
                    if elem.split('_')[-1].isdigit():
                        t3[i][j] = t3[i][j].rsplit('_')[0]

            t2 = sorted(t2)
            for i, t in enumerate(t3):
                t3[i] = sorted(t)

            for a_group in aff_groups:

                a_group_o = deepcopy(a_group)
                match = True

                # match reaction names
                for t in t_group:
                    if t[1] != reaction_name:
                        match = False
                for a in a_group:
                    if a[1] != reaction_name:
                        match = False

                a2 = []
                a3 = []
                for a in a_group:
                    a2.append(a[2])
                    a3.append(a[3])

                for i, a in enumerate(a3):
                    for j, elem in enumerate(a):
                        if elem.split('_')[-1].isdigit():
                            a3[i][j] = a3[i][j].rsplit('_')[0]

                a2 = sorted(a2)
                for i, a in enumerate(a3):
                    a3[i] = sorted(a)
                    if a3[i] != t2:
                        match = False
                for i, t in enumerate(t3):
                    if t != a2:
                        match = False

                if match:
                    interact = [[], [], [], [], reaction_name]
                    for t in t_group_o:
                        interact[1].append(t[0])
                    for a in a_group_o:
                        interact[0].append(a[0])

                    interact[2] = deepcopy(t_group_o[0][3])
                    interact[3] = deepcopy(a_group_o[0][3])

                    interacts.append(interact)

        # # add instructions and partners
        # for i, each in enumerate(interacts):
        #     interacts[i].extend([[], []])
        #     for item in self.library[each[3][0]]:
        #         if item.reaction == each[5]:
        #             for j, every in enumerate(item.instructions):
        #                 if every not in interacts[i][6]:
        #                     interacts[i][6].append(every)
        #             for j, every in enumerate(item.partners):
        #                 if every not in interacts[i][7]:
        #                     interacts[i][7].append(every)

        return interacts

    def _findMotifs(self, targs, motif_build, interacts, in_nodes, mots):

        self.count += 1
        # todo: apply filters to partial motifs

        relax = self.nodes[motif_build[0][0][0]].relax
        working_motif = deepcopy(motif_build)

        def max_path_finder(node_i, node_paths, m2, connects):

            for k in connects[node_i]:
                if k != 0:
                    if k not in node_paths:
                        npa = deepcopy(node_paths)
                        npa.append(k)
                        max_path_finder(k, npa, m2, connects)
                else:
                    if len(node_paths) > m2[0]:
                        m2[0] = len(node_paths)

        # motif filters

        # =======================================================================

        # depth check

        depth_g = True

        new_working_motif = []
        for each in working_motif:
            new_working_motif.append(each[:2])

        # simplify the paths for scoring
        simple_motif = deepcopy(new_working_motif)
        still_simplifying = True
        while still_simplifying:
            still_simplifying = False
            for i, item in enumerate(simple_motif):
                if item[1] and len(item[0]) > 1:
                    for every in simple_motif:
                        if every[1] and set(every[0] + every[1]) == set(item[0]):
                            simple_motif[i][0] = every[1]
                            still_simplifying = True
                if item[1] and len(item[1]) > 1:
                    for every in simple_motif:
                        if every[1] and set(every[0] + every[1]) == set(item[1]):
                            simple_motif[i][1] = every[0]
                            still_simplifying = True

        # list all species
        seed = None
        for each in simple_motif:
            if not each[1]:
                seed = each[0]
        incidents = [seed]
        for each in simple_motif:
            if each[0] not in incidents:
                incidents.append(each[0])
            if each[1] and each[1] not in incidents:
                incidents.append(each[1])

        # find all connections
        connections = [[] for _ in incidents]
        for i, each in enumerate(incidents):
            for item in simple_motif:
                if each == item[0]:
                    for j, every in enumerate(incidents):
                        if every == item[1]:
                            connections[i].append(j)

        # find max path
        max_path = [0]
        for i, each in enumerate(incidents):
            max_path_finder(i, [i], max_path, connections)

        # test depth
        target = None
        for each in working_motif:
            if not each[1]:
                target = each[0][0]
        if max_path[0] > min(self.nodes[target].depth, self.max_depth):
            depth_g = False
            # return

        # =======================================================================

        # # cycle_check; eliminates motifs with cycles
        #
        # still_adding = True
        # paths = [working_motif[0][0]]
        # while still_adding:
        #     still_adding = False
        #     new_paths = []
        #     for p in paths:
        #         for inter in working_motif[1:]:
        #             if p[-1] in inter[1]:
        #                 for each in inter[0]:
        #                     q = deepcopy(p)
        #                     q.append(each)
        #                     new_paths.append(q)
        #     for p in new_paths:
        #         if p not in paths:
        #             paths.append(p)
        #             still_adding = True
        #     for each in paths:
        #         if each[-1] in each[:-1]:
        #             return

        # =======================================================================

        # target path check

        path = True
        new_working_motif = []
        for inter in working_motif:
            if not inter[1]:
                new_working_motif.append(inter)
        still_searching = True
        while still_searching:
            still_searching = False
            working_temp = []
            for inter1 in working_motif:
                for inter2 in new_working_motif:
                    if inter1[1] and set(inter1[1]).issubset(set(inter2[0])) and inter1 not in new_working_motif \
                            and inter1 not in working_temp:
                        working_temp.append(inter1)
                        still_searching = True
            new_working_motif.extend(working_temp)
        if len(new_working_motif) != len(working_motif):
            path = False

        # =======================================================================

        # reaction check; guarantees the specified reactions are in the motif

        reaction = True
        target = working_motif[0][0][0]
        for each in self.nodes[target].reactions:
            rc = False
            for item in working_motif:
                if each[0] == item[5] and set(each[1]) == set(item[0]) and set(each[2]) == set(item[1]):
                    rc = True
            if not rc:
                reaction = False

        # =======================================================================

        # motif coverage check

        coverage = False
        node_coverage = set()
        for inter in working_motif[1:]:
            for species in inter[0]:
                node_coverage.add(species)
        if set.intersection(node_coverage, set(deepcopy(in_nodes))) == set(deepcopy(in_nodes)):
            coverage = True

        if relax:
            coverage = True

        # =======================================================================

        # complex coverage check; checks if complexes are complete motifs

        complex_coverage = True
        # for inter1 in working_motif[1:]:
        #     if len(inter1[0]) > 1:
        #         current_inter = False
        #         for inter2 in working_motif[1:]:
        #             if set(inter2[0] + inter2[1]) == set(inter1[0]):
        #                 current_inter = True
        #         if not current_inter:
        #             complex_coverage = False
        # if relax:
        #     complex_coverage = True

        # =======================================================================

        if depth_g and path and reaction and coverage and complex_coverage:

            if working_motif not in mots:
                mots.append(working_motif)

        # find interactions for the current targs
        current_interactions = []
        new_interacts = deepcopy(interacts)
        for targ in targs:
            for n, inter in reversed(list(enumerate(new_interacts))):
                inter_set = set(deepcopy(inter[1]))
                targ_set = set(deepcopy(targ[0]))
                if inter_set.issubset(targ_set):
                    current_interactions.append(new_interacts.pop(n))

        # add the new interactions to create new motifs
        for m in range(len(current_interactions)):
            for combin in combinations(current_interactions, m+1):
                new_targs = deepcopy(list(combin))
                new_motif_build = deepcopy(motif_build) + deepcopy(list(combin))

                self._findMotifs(new_targs, new_motif_build, new_interacts, in_nodes, mots)

    @staticmethod
    def _motif_size(motif):

        def path_finder(node_i, node_paths, s, connects):

            for k in connects[node_i]:
                if k != 0:
                    if k not in node_paths:
                        npa = deepcopy(node_paths)
                        npa.append(k)
                        path_finder(k, npa, s, connects)
                else:
                    s[0] += len(node_paths)

        motif2 = []
        for each in motif:
            motif2.append(each[:2])

        # simplify the motif for scoring
        simple_motif = deepcopy(motif2)
        still_simplifying = True
        while still_simplifying:
            still_simplifying = False
            for i, item in enumerate(simple_motif):
                if item[1] and len(item[0]) > 1:
                    for every in simple_motif:
                        if every[1] and set(every[0] + every[1]) == set(item[0]):
                            simple_motif[i][0] = every[1]
                            still_simplifying = True
                if item[1] and len(item[1]) > 1:
                    for every in simple_motif:
                        if every[1] and set(every[0] + every[1]) == set(item[1]):
                            simple_motif[i][1] = every[0]
                            still_simplifying = True

        # duplicate reaction elimination
        simple_temp = []
        for each in simple_motif:
            if each not in simple_temp:
                simple_temp.append(deepcopy(each))

        simple_motif = simple_temp

        # list all species
        seed = None
        for each in simple_motif:
            if not each[1]:
                seed = each[0]
        incidents = [seed]
        for each in simple_motif:
            if each[0] not in incidents:
                incidents.append(each[0])
            if each[1] and each[1] not in incidents:
                incidents.append(each[1])

        # find all connections
        connections = [[] for _ in incidents]
        for i, each in enumerate(incidents):
            for item in simple_motif:
                if each == item[0]:
                    for j, every in enumerate(incidents):
                        if every == item[1]:
                            connections[i].append(j)

        # score motif
        score = [len(motif)-1]
        for i, each in enumerate(incidents):
            path_finder(i, [i], score, connections)

        return score[0]

    def _filter_interactions(self, root, interacts):

        new_interacts = []
        if self.nodes[root].exclusive:
            reactions = deepcopy(self.nodes[root].reactions)
            reactions.extend(deepcopy(self.nodes[root].optional_rxns))
            for each in reactions:
                for item in interacts:
                    if each[0] == item[5] and set(each[1]) == set(item[0]) and set(each[2]) == set(item[1]):
                        new_interacts.append(item)
        else:
            new_interacts = interacts

        return new_interacts

    def _build_motifs(self):

        for node in self.nodes:

            self.count = 0
            self.nodes[node].self = self._findSelfInteractions(node)
            nodelist = deepcopy(self.nodes[node].incidentNodes)
            if not nodelist:
                self.nodes[node].motifs.extend([[0, [[[node], None, None, None, None]]]])
            else:
                seed = [[[node], None, None, None, None]]
                interactions = []
                interactions += self._findInteractions(node, nodelist)
                # quit()
                if int(self.nodes[node].depth) > 1:  # default: complex
                    for j, thing in enumerate(nodelist):
                        temp_nodelist = deepcopy(nodelist)
                        temp_nodelist.pop(j)
                        interactions += self._findInteractions(thing, temp_nodelist)

                motifs = []
                interactions = self._filter_interactions(node, interactions)
                self._findMotifs(seed, seed, interactions, self.nodes[node].incidentNodes, motifs)
                if motifs:
                    self.nodes[node].motifs.extend(motifs)
                else:
                    self.nodes[node].motifs.extend([[[[node], None, None, None, None]]])


class Combine_and_Build:

    def __init__(self, nodes, library, top=1):
        self.nodes = nodes
        self.library = library
        self.models = []
        self.reduced_models = []
        self.pysb_models = []
        self._combine(top)

    def _combine(self, top):

        node_list = []
        motif_list = []
        for node in self.nodes:
            node_list.append(node)
            motif_list.append(self.nodes[node].motifs)

        motif_combos = []
        num_motifs = []
        for each in motif_list:
            num_motifs.append(range(len(each)))

        for indices in product(*num_motifs):
            motif = []
            for i, each in enumerate(motif_list):
                motif.append(each[indices[i]])
            motif_combos.append(motif)

        # else:
        #     top_list = [[score, [0 for _ in self.nodes]]]
        #     for i in range(top-1):
        #         top_list.append([1000000+i, []])
        #     current_model_list = [[score, [0 for _ in self.nodes]]]
        #
        #     # todo: this is inefficient
        #     # while new potential models are still being generated
        #     while current_model_list:
        #
        #         # a round of bfs
        #         new_model_list = []
        #         for model in current_model_list:
        #             for i, node in enumerate(model[1]):
        #                 if model[1][i]+1 < len(motif_list[i]):
        #                     new_model = model[1][:]
        #                     new_model[i] += 1
        #                     score = 0
        #                     for j, each in enumerate(motif_list):
        #                         score += each[new_model[j]][0]
        #                     if score <= top_list[-1][0] and [score, new_model] not in new_model_list:
        #                         new_model_list.append([score, new_model])
        #         new_model_list.sort()
        #         current_model_list = []
        #
        #         # add models to top list, sort and filter
        #         for each in new_model_list:
        #             if each[0] <= top_list[-1]:
        #                 top_list.append(deepcopy(each))
        #                 top_list.sort()
        #                 top_values = []
        #                 for item in top_list:
        #                     if item[0] not in top_values:
        #                         top_values.append(item[0])
        #                 top_values = top_values[:top]
        #
        #                 new_top_list = []
        #                 for i, item in enumerate(top_list):
        #                     if item[0] in top_values:
        #                         new_top_list.append(deepcopy(item))
        #                 top_list = new_top_list
        #                 current_model_list.append(deepcopy(each))
        #
        #     tl2 = []
        #     for each in top_list:
        #         if each[1]:
        #             tl2.append(each)
        #     top_list = tl2
        #
        #     # find combinations of motifs
        #     for each in top_list:
        #         motifs = []
        #         for i, node in enumerate(node_list):
        #             motifs.append(self.nodes[node].motifs[each[1][i]][1])
        #         motif_combos.append(motifs)

        # reorder interactions
        for i, each in enumerate(motif_combos):
            for j, item in enumerate(each):
                if item[0][1]:
                    for k, every in enumerate(item):
                        if not every[1]:
                            motif_combos[i][j].insert(0, motif_combos[i][j].pop(k))

        for i, each in enumerate(motif_combos):

            name = 'model_' + str(i) + '_motifs'
            f = open('/home/mak/models/motifs/'+name, 'w+')

            for item in each:
                for every in item:
                    f.write(str(every) + '\n')
                f.write('\n')
            f.close()

        # generate list of models, expanding for proxy nodes
        for i, each in enumerate(motif_combos):
            copy_nodes = deepcopy(self.nodes)
            for item in each:
                c_node = None
                for em in item:
                    if not em[1]:
                        c_node = em[0][0]
                copy_nodes[c_node].motifs = item
            table_combos = []
            n = 0

            for item in copy_nodes:
                t = copy_nodes[item].table[1:]
                t0 = list(product([True, False], repeat=0))

                tables = []
                for every in t0:
                    t2 = [[]]
                    t2[-1].extend(copy_nodes[item].table[0])
                    for thing in t:
                        t2.append(list(every))
                        t2[-1].extend(thing)
                    tables.append(t2)
                table_combos.append(tables)
                n += 1
            table_list = list(product(*table_combos))

            for item in table_list:
                copy_copy_nodes = deepcopy(copy_nodes)
                for every in item:
                    copy_copy_nodes[every[0][-1]].table = every
                ModelBuilder(i, copy_copy_nodes, self.library, self.pysb_models)


class ModelBuilder(Builder):
    """
    Assemble a PySB model from a Boolean model.
    """
    def __init__(self, num, nodes, library, pysb_models):

        super(ModelBuilder, self).__init__()
        self.num = num
        self.nodes = nodes
        self.library = library
        self.pysb_models = pysb_models
        self.monomer_info = defaultdict(list)
        self.action_info = defaultdict(list)
        self.base_states = defaultdict(list)
        self.active_states = defaultdict(list)
        self.inactive_states = defaultdict(list)
        self.iv_inactive_states = defaultdict(list)
        self._build()
        self._export()

    def _build(self):

        self._parse()
        self._get_monomer_info()
        self._add_monomers()
        self._add_rules()
        self._add_initials()
        self._add_observables()

    def _parse(self):

        # parse library reactions
        for every in self.library:
            for thing in self.library[every]:
                for each in thing.rxnTemplates:
                    mols = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', each)
                    ops = re.findall(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', each)
                    parced = []
                    for m in mols:
                        parced.append([])
                        sites = []
                        states = []
                        if '(' in m:
                            parced[-1].append(m[:m.index('(')])
                            if '()' not in m:
                                ms = re.split(r'\s*=\s*|\s*,\s*', m[m.index('(') + 1:-1])
                                sites.extend(ms[::2])
                                for s in ms[1::2]:
                                    states.append(s)
                        else:
                            parced[-1].append(m)
                        parced[-1].append(sites)
                        parced[-1].append(states)
                    parced_rxn = [parced.pop(0)]
                    for i, m in enumerate(parced):
                        parced_rxn.append(ops[i])
                        parced_rxn.append(m)
                    thing.rxnsParsed.append(parced_rxn)

        # for node in self.nodes:
        #     for inter in self.nodes[node].motifs:
        #         if inter[1]:
        #             for every in self.library[inter[3][0]]:
        #                 if every.reaction == inter[4]:
        #                     pass

    def _export(self):

        # self.model.name = 'model_' + str(self.num)
        # self.pysb_models.append(self.model)

        # ------------------------------

        # f = open('/home/mak/models/models/model_' + str(self.num) + '.py', 'w+')
        # f.write(PysbFlatExporter(self.model).export())
        # f.close()
        #
        # f = open('/home/mak/models/models/model_' + str(self.num) + '.py', "r")
        # contents = f.readlines()
        # f.close()
        #
        # contents.insert(2, 'import numpy as np\nfrom pysb.integrate import Solver\nimport pylab as pl\nimport matplotlib.pyplot as plt\n')
        #
        # f = open('/home/mak/models/models/model_' + str(self.num) + '.py', "w")
        # contents = "".join(contents)
        # f.write(contents)
        # f.close()

        # ------------------------------

        f = open('model_' + str(self.num) + '.py', 'w+')
        f.write(PysbFlatExporter(self.model).export())
        f.close()

        f = open('model_' + str(self.num) + '.py', "r")
        contents = f.readlines()
        f.close()

        contents.insert(2, 'import numpy as np\nfrom pysb.integrate import Solver\nimport pylab as pl\nimport matplotlib.pyplot as plt\n')

        f = open('model_' + str(self.num) + '.py', "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()

    def _find_sites(self, interaction):

        if interaction[1]:

            monomer_names = []
            monomer_labels = []
            monomer_sites = []

            # get names and labels from interaction and initialize lists
            for i, each in enumerate(interaction[0]):
                monomer_names.append(each)
                monomer_labels.append(interaction[2][i])
                monomer_sites.append({})
            for i, each in enumerate(interaction[1]):
                monomer_names.append(each)
                monomer_labels.append(interaction[3][i])
                monomer_sites.append({})

            # break down reaction template and add sites and states
            for every in self.library[interaction[3][0]]:
                if every.reaction == interaction[4]:
                    for thing in every.rxnTemplates:
                        rxn_template = thing
                        rxnTemp = re.split(r'\s*:', rxn_template)[0]
                        mol_list = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', rxnTemp)
                        for m in mol_list:
                            if '(' in m and '()' not in m:
                                mol = m[:m.index('(')]
                                sites = re.split(r'\s*=\s*|\s*,\s*', m[m.index('(') + 1:-1])[::2]
                                states = re.split(r'\s*=\s*|\s*,\s*', m[m.index('(') + 1:-1])[1::2]
                                for i, each in enumerate(sites):
                                    if each not in monomer_sites[monomer_labels.index(mol)]:
                                        monomer_sites[monomer_labels.index(mol)][each] = []
                                        monomer_sites[monomer_labels.index(mol)][each].append(states[i])
                                    else:
                                        if states[i] not in monomer_sites[monomer_labels.index(mol)][each]:
                                            monomer_sites[monomer_labels.index(mol)][each].append(states[i])

            # print
            # print 'mn', monomer_names
            # print 'ml', monomer_labels
            # print 'ms', monomer_sites
            # print

            # rename the sites appropriately
            for i, each in enumerate(monomer_names):
                for item in monomer_sites[i]:
                    if item in monomer_labels:
                        monomer_sites[i][monomer_names[monomer_labels.index(item)]] = monomer_sites[i].pop(item)
                    # s1 = item[:item.rfind('_')]
                    # s2 = item[item.rfind('_')+1:]
                    # if s1 in monomer_labels and s2 == 's':
                    #     monomer_sites[i][monomer_names[monomer_labels.index(s1)]+'_s'] = monomer_sites[i].pop(item)
                    # if s1 in monomer_labels and s2.isdigit():
                    #     monomer_sites[i][monomer_names[monomer_labels.index(s1)]+'_'+s2] = monomer_sites[i].pop(item)

            # add sites to monomer_info
            for i, each in enumerate(monomer_names):
                for item in monomer_sites[i]:
                    if item not in self.monomer_info[each][0]:
                        self.monomer_info[each][0].append(item)

            # initiate state sites in monomer_info
            for i, each in enumerate(monomer_names):
                for j, item in enumerate(monomer_sites[i]):
                    if item in self.monomer_info or (item[:item.rfind('_')] and item[item.rfind('_')+1:].isdigit()):
                        pass
                    else:
                        if item in self.monomer_info[each][1]:
                            pass
                        else:
                            self.monomer_info[each][1][item] = []

            # add state sites to monomer_info
            for i, each in enumerate(monomer_names):
                for j, item in enumerate(monomer_sites[i]):
                    if item in self.monomer_info or (item[:item.rfind('_')] and item[item.rfind('_')+1:].isdigit()):
                        pass
                    else:
                        for every in monomer_sites[i][item]:
                            if every not in self.monomer_info[each][1][item]:
                                self.monomer_info[each][1][item].append(every)

    def _get_monomer_info(self):

        # initiate monomer information
        for motif in self.nodes:
            for interaction in self.nodes[motif].motifs:
                if interaction[1]:
                    for i, species in enumerate(interaction[0]):
                        if species not in self.monomer_info:
                            self.monomer_info[species] = [[], {}]
                    for i, species in enumerate(interaction[1]):
                        if species not in self.monomer_info:
                            self.monomer_info[species] = [[], {}]
                else:
                    for i, species in enumerate(interaction[0]):
                        if species not in self.monomer_info:
                            self.monomer_info[species] = [[], {}]

        # find binding and state sites; inferred from motifs/self interactions and library
        for node in self.nodes:
            for interaction in self.nodes[node].self:
                self._find_sites(interaction)

        for node in self.nodes:
            for interaction in self.nodes[node].motifs:
                self._find_sites(interaction)

        # print
        # print 'MONOMER INFO'
        # for each in self.monomer_info:
        #     print each, self.monomer_info[each]
        # print
        # quit()

    def _add_monomers(self):

        for each in self.monomer_info:
            if self.nodes[each].compete:
                sites = self.monomer_info[each][0]
                for item in self.nodes[each].compete:
                    comps = item.split(':')
                    comp_site = ''
                    for every in comps:
                        comp_site += every + '_'
                    comp_site = comp_site[:-1]
                    for i, every in reversed(list(enumerate(sites))):
                        if every in comps:
                            sites.pop(i)
                    sites.append(comp_site)
                self.monomer(each, sites, self.monomer_info[each][1])
            else:
                self.monomer(each, self.monomer_info[each][0], self.monomer_info[each][1])

    def _get_action_info(self):

        for motif in self.nodes:
            motif_action = []
            for interaction in self.nodes[motif].motifs:
                if interaction[1]:
                    reactants = interaction[0]
                    targets = interaction[1]
                    target_states = [[[], [], []] for _ in range(len(targets))]
                    for every in self.library[interaction[3][0]]:
                        if every.reaction == interaction[4]:
                            for thing in every.rxnTemplates:
                                mol_list = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', thing)
                                mol_list2 = []
                                for m in mol_list:
                                    mol_list2.append(re.split(r'\s*\(\s*|\s*,\s*|\s*=\s*', m[:-1]))
                                for i, each in enumerate(mol_list2):
                                    for j, item in enumerate(each):
                                        for k, stuff in enumerate(interaction[2]):
                                            if stuff == item:
                                                mol_list2[i][j] = interaction[0][k]
                                            if stuff + '_s' == item:
                                                mol_list2[i][j] = interaction[0][k] + '_s'
                                        for k, stuff in enumerate(interaction[3]):
                                            if stuff == item:
                                                mol_list2[i][j] = interaction[1][k]
                                            if stuff + '_s' == item:
                                                mol_list2[i][j] = interaction[1][k] + '_s'
                                reactant_present = [False for _ in range(len(reactants))]
                                for each in mol_list2:
                                    for i, item in enumerate(reactants):
                                        if each[0] == item:
                                            reactant_present[i] = True
                                if all(reactant_present):
                                    for i, each in enumerate(mol_list2):
                                        if len(each) == 2 and each[1] == '':
                                            mol_list2[i].append('')
                                    for each in mol_list2:
                                        for i, item in enumerate(targets):
                                            if each[0] == item and each[1]:
                                                if not target_states[i][0]:
                                                    target_states[i][0] = each[1::2]
                                                    target_states[i][1] = each[2::2]
                                                else:
                                                    target_states[i][2] = each[2::2]
                                    for i in range(len(target_states)):
                                        if not target_states[i][2]:
                                            target_states[i][2] = target_states[i][1]
                    motif_action.append([interaction[0], interaction[1], target_states])
            self.action_info[motif] = motif_action

        # print
        # print 'ACTION INFO'
        # for each in self.action_info:
        #     for item in self.action_info[each]:
        #         print each, item
        # print
        # quit()

    def _add_rules(self):

        self._get_action_info()

        # print
        # print 'ACTION INFO'
        # for each in self.action_info:
        #     for item in self.action_info[each]:
        #         print each, item
        # print

        # create dictionary of base states
        for each in self.monomer_info:
            mon_obj = self.model.monomers[each]
            states = {}
            for item in mon_obj.sites:
                if item in mon_obj.site_states:
                    for every in self.nodes[each].labels:
                        for rxn in self.library[every]:
                            for base in rxn.base_states:
                                if base.split('=')[0] == item:
                                    states[item] = base.split('=')[1]
                else:
                    states[item] = 'None'
            state_list1 = []
            state_list2 = []
            for item in states:
                state_list1.append(item)
                state_list2.append(states[item])
            state_list3 = [deepcopy(state_list1), deepcopy(state_list2)]
            self.base_states[each] = state_list3

        # print 'BASE STATES'
        # for each in self.base_states:
        #     print each, self.base_states[each]
        # print
        # quit()

        # create dictionary of active states based on the Boolean equations
        for each in self.monomer_info:
            a_states = []
            i_states = []
            base_sites = {}
            for i, item in enumerate(self.base_states[each][0]):
                base_sites[item] = self.base_states[each][1][i]

            # Takes the first state for state sites without a base state.
            # Is that warranted? Should we force a base state and leave those
            # without as dwdc?
            # Either way the base states should be definable in the label file
            # for better refinement.
            for i, item in enumerate(self.monomer_info[each][0]):
                if item not in base_sites:
                    if item not in self.monomer_info[each][1]:
                        base_sites[item] = None
                    else:
                        base_sites[item] = self.monomer_info[each][1][item][0]

            if each in self.nodes:  # if model node; else proxy node
                table = self.nodes[each].table

                # find altered states and determine dominance for conflicting states
                altered_states = {}
                for item in self.action_info[each]:
                    for i, every in enumerate(item[2][0][0]):
                        if item[2][0][1][i] != item[2][0][2][i] and item[1][0] == each:
                            altered_states[every] = []

                for item in self.action_info[each]:
                    for i, every in enumerate(item[2][0][0]):
                        if item[2][0][1][i] != item[2][0][2][i] and item[1][0] == each:
                            altered_states[every].append([item[0], [item[2][0][1][i], item[2][0][2][i]]])

                per_site_counts = {}
                for item in altered_states:
                    state_table = deepcopy(table)
                    state_table[0].append('states')
                    for i, st in enumerate(state_table):
                        if i > 0:
                            state_table[i].append(set())
                    for i, st in enumerate(state_table):
                        if i > 0:
                            for every in altered_states[item]:
                                all_true = [False for _ in every[0]]
                                for k, thing in enumerate(every[0]):
                                    for j, lotsa in enumerate(state_table[0][:-2]):
                                        if thing == lotsa and state_table[i][j]:
                                            all_true[k] = True
                                if all(all_true):
                                    state_table[i][-1].add(every[1][1])
                                else:
                                    state_table[i][-1].add(every[1][0])

                    active_state_count = {}
                    for every in state_table[1:]:
                        if every[-2]:
                            for thing in every[-1]:
                                if thing in active_state_count:
                                    active_state_count[thing] += 1
                                else:
                                    active_state_count[thing] = 1

                    per_site_counts[item] = active_state_count

                table2 = [[] for _ in range(len(table))]
                direct = []
                for item in self.action_info[each]:
                    if item[1][0] == each:
                        direct.extend(item[0])
                direct = set(direct)
                for i, item in enumerate(table[0]):
                    if item in direct:
                        for j, thing in enumerate(table):
                            table2[j].append(table[j][i])

                for i, item in enumerate(table):
                    table2[i].append(table[i][-1])

                table3 = [table2[0]]
                for item in table2[1:]:
                    if item[-1]:
                        table3.append(item)
                for item in table2[1:]:
                    if not item[-1]:
                        t = deepcopy(item)
                        t[-1] = True
                        if t not in table3 and item not in table3:
                            table3.append(item)

                for item in table[1:]:  # find all combinations of active incident nodes that activate the target
                    if item[-1]:
                        boole = defaultdict(bool)  # dictionary of incident node Boolean states on a given row
                        for i, every in enumerate(table[0][:-1]):
                            boole[every] = item[i]

                        # run through all actions for a given motif; consider only those acting on the target node
                        sites = deepcopy(base_sites)
                        changed_sites = {}

                        # add sites for active incident nodes
                        for every in self.action_info[each]:
                            if every[1][0] == each:
                                for thing in every[0]:
                                    if thing in boole and not boole[thing]:
                                        for j, lotsa in enumerate(every[2][0][0]):
                                            sites[lotsa] = every[2][0][1][j]
                                            if lotsa not in changed_sites:
                                                changed_sites[lotsa] = [every[2][0][1][j]]
                                            else:
                                                changed_sites[lotsa].append(every[2][0][1][j])
                        for every in self.action_info[each]:
                            if every[1][0] == each:
                                for thing in every[0]:
                                    if thing in boole and boole[thing]:
                                        for j, lotsa in enumerate(every[2][0][0]):
                                            sites[lotsa] = every[2][0][2][j]
                                            if lotsa not in changed_sites:
                                                changed_sites[lotsa] = [every[2][0][2][j]]
                                            else:
                                                changed_sites[lotsa].append(every[2][0][2][j])
                        cs2 = {}
                        for every in changed_sites:
                            if len(changed_sites[every]) > 1:
                                cs2[every] = changed_sites[every]
                        sites = deepcopy(base_sites)
                        for every in sites:
                            if every in changed_sites:
                                if len(changed_sites[every]) > 1:
                                    dom = None
                                    dom_count = 0
                                    for thing in per_site_counts[every]:
                                        if per_site_counts[every][thing] > dom_count:
                                            dom = thing
                                            dom_count = per_site_counts[every][thing]
                                    sites[every] = dom
                                else:
                                    sites[every] = changed_sites[every][0]
                        if sites not in a_states:
                            a_states.append(sites)

                for item in table3[1:]:
                    if not item[-1] and base_sites:
                        boole = defaultdict(bool)  # dictionary of incident node Boolean states on a given row
                        for i, every in enumerate(table3[0][:-1]):
                            boole[every] = item[i]

                        # run through all actions for a given motif; consider only those acting on the target node
                        sites = deepcopy(base_sites)

                        # add sites for active incident nodes
                        for every in self.action_info[each]:
                            if every[1][0] == each:
                                for thing in every[0]:
                                    if thing in boole and not boole[thing]:
                                        for j, lotsa in enumerate(every[2][0][0]):
                                            sites[lotsa] = every[2][0][1][j]
                        for every in self.action_info[each]:
                            if every[1][0] == each:
                                for thing in every[0]:
                                    if thing in boole and boole[thing]:
                                        for j, lotsa in enumerate(every[2][0][0]):
                                            sites[lotsa] = every[2][0][2][j]

                        if sites not in i_states and sites not in a_states:
                            i_states.append(sites)

                self.active_states[each] = a_states
                if i_states:
                    self.inactive_states[each] = i_states
                else:
                    self.inactive_states[each] = [{}]

        # print
        # print 'ACTIVE STATES'
        # for each in self.active_states:
        #     for item in self.active_states[each]:
        #         print each, item
        # print

        # print 'INACTIVE STATES'
        # for each in self.inactive_states:
        #     for item in self.inactive_states[each]:
        #         print each, item
        # print
        # quit()

        # create rules for self interactions
        # used_self_interactions = []
        for node in self.nodes:
            for interaction in self.nodes[node].self:
                # used_self_interactions.append(interaction)

                # retrieve current rxns and reaction templates
                current_rxn = None
                for rxn in self.library[interaction[3][0]]:
                    if rxn.reaction == interaction[4]:
                        current_rxn = rxn
                n = 0
                for temp in current_rxn.rxnTemplates:

                    rxnTemp1 = re.split(r'\s*:\s*', temp)
                    rxnTemp = rxnTemp1[0]

                    # split template
                    rxn_split = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', rxnTemp)
                    rxn_split_mols = []
                    for each in rxn_split:
                        if each == 'None':
                            rxn_split_mols.append(['None', [], []])
                        else:
                            mol = each.split('(')[0]
                            sites = each.split('(')[1][:-1]
                            site_names = []
                            site_values = []
                            if sites:
                                sites = sites.split(',')
                                for item in sites:
                                    site_names.append(item.split('=')[0].strip())
                                    site_values.append(item.split('=')[1].strip())
                            rxn_split_mols.append([mol, site_names, site_values])

                    ops = re.findall(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', rxnTemp)
                    for i, each in enumerate(ops):
                        ops[i] = ops[i].strip()

                    rxn_split_parsed = [rxn_split_mols[0]]
                    for i, each in enumerate(ops):
                        rxn_split_parsed.append(each)
                        rxn_split_parsed.append(rxn_split_mols[i + 1])
                    # print rxn_split_parsed

                    # Make substitutions from the interaction to the template.
                    for j, each in enumerate(interaction[3]):
                        for i, item in enumerate(rxn_split_parsed):
                            if item != '>>' and item != '<>' and item != '+' and item != '%':
                                if item[0] == each:
                                    rxn_split_parsed[i][0] = interaction[1][j]
                                for k, every in enumerate(item[1]):
                                    if every == each:
                                        rxn_split_parsed[i][1][k] = interaction[1][j]
                                    if every[:every.rfind('_')] == each and every[every.rfind('_')+1:].isdigit():
                                        rxn_split_parsed[i][1][k] = interaction[1][j] + every[every.rfind('_'):]

                    # print rxn_split_parsed

                    # add base states for bond and states sites undefined in the reaction template
                    for j, species in enumerate(rxn_split_parsed):
                        if isinstance(species, list) and species[0] != 'None':
                            for i, each in enumerate(self.base_states[species[0]][0]):
                                if each not in species[1]:
                                    if each in self.nodes or each[:each.rfind('_')] in self.nodes:
                                        rxn_split_parsed[j][1].append(each)
                                        rxn_split_parsed[j][2].append(self.base_states[species[0]][1][i])

                    # print node, rxn_split_parsed

                    # routine for synthesis rxns
                    if rxn_split_parsed[0][0] == 'None':
                        rxn_split_parsed[2][1] = self.base_states[rxn_split_parsed[2][0]][0]
                        rxn_split_parsed[2][2] = self.base_states[rxn_split_parsed[2][0]][1]

                    # routine for degradation rxns
                    if rxn_split_parsed[2][0] == 'None':
                        binding_states = [[], []]
                        for i, each in enumerate(self.base_states[rxn_split_parsed[0][0]][0]):
                            if each in self.nodes:
                                binding_states[0].append(deepcopy(self.base_states[rxn_split_parsed[0][0]][0][i]))
                                binding_states[1].append(deepcopy(self.base_states[rxn_split_parsed[0][0]][1][i]))
                        rxn_split_parsed[0][1].extend(binding_states[0])
                        rxn_split_parsed[0][2].extend(binding_states[1])

                    # substitute integer for 'integer' and None for 'None'
                    for i, item in enumerate(rxn_split_parsed):
                        if item != '+' and item != '%' and item != '>>' and item != '<>':
                            for j, every in enumerate(item[2]):
                                if item[1][j][:item[1][j].rfind('_')] in self.nodes and every != 'None' and every != 'ANY' and every != 'WILD':
                                    rxn_split_parsed[i][2][j] = int(rxn_split_parsed[i][2][j])
                                if item[1][j] in self.nodes and every != 'None' and every != 'ANY' and every != 'WILD':
                                    rxn_split_parsed[i][2][j] = int(rxn_split_parsed[i][2][j])
                                if every == 'None':
                                    rxn_split_parsed[i][2][j] = None
                                if every == 'ANY':
                                    rxn_split_parsed[i][2][j] = ANY
                                if every == 'WILD':
                                    rxn_split_parsed[i][2][j] = WILD

                    # define the rule rule_name
                    rule_name = ''
                    for item in interaction[1]:
                        if item:
                            rule_name += item + '_'
                    for item in interaction[0]:
                        if item:
                            rule_name += item + '_'
                        else:
                            rule_name = rule_name[:-1] + '_'

                    rule_name += interaction[4]
                    rule_name += '_' + str(n)

                    # define monomer patterns
                    mon_pats = []
                    for item in rxn_split_parsed:
                        if item == '+' or item == '%' or item == '>>' or item == '<>':
                            mon_pats.append(item)
                        else:
                            if item[0] == 'None':
                                mon_pats.append('None')
                            else:
                                mon_states = {}
                                for i, every in enumerate(item[1]):
                                    mon_states[every] = item[2][i]
                                mon_obj = self.model.monomers[item[0]]
                                mon_pats.append(MonomerPattern(mon_obj, mon_states, None))

                    # define complex patterns
                    com_pats_temp = [[]]
                    for item in mon_pats:
                        if item == '>>' or item == '<>':
                            com_pats_temp.extend([item, []])
                        elif item == '+':
                            com_pats_temp.append([])
                        elif item == '%':
                            pass
                        else:
                            com_pats_temp[-1].append(item)
                    com_pats = []
                    for item in com_pats_temp:
                        if item == '>>' or item == '<>':
                            com_pats.append(item)
                        elif item == ['None']:
                            pass
                        else:
                            com_pats.append(ComplexPattern(item, None))

                    # define reversibility and split patterns into reactants and products
                    react_com_pats = []
                    prod_com_pats = []
                    carrot = 0
                    reversible = None
                    for item in com_pats:
                        if item == '<>':
                            carrot = 1
                            reversible = True
                        elif item == '>>':
                            carrot = 1
                            reversible = False
                        else:
                            if carrot == 0:
                                react_com_pats.append(item)
                            if carrot == 1:
                                prod_com_pats.append(item)

                    # define rule expression
                    order = [len(react_com_pats), len(prod_com_pats)]
                    rule_exp = RuleExpression(ReactionPattern(react_com_pats), ReactionPattern(prod_com_pats), reversible)

                    if reversible:
                        forward = rule_name + '_' + str(order[0]) + 'kf'
                        self.parameter(forward, 1)
                        reverse = rule_name + '_' + str(order[1]) + 'kr'
                        self.parameter(reverse, 1)
                        self.rule(rule_name, rule_exp, self.model.parameters[forward], self.model.parameters[reverse])
                    else:
                        forward = rule_name + '_' + str(order[0]) + 'kc'
                        self.parameter(forward, 1)
                        self.rule(rule_name, rule_exp, self.model.parameters[forward])
                    n += 1

        # create rules for motif interactions

        # routine for mechanistic output
        # todo: add functionality for explicitly determining reactant and target states through instructions

        used_interactions = []

        for node in self.nodes:

            # initialize motif specific reactant states
            reactant_states = defaultdict(list)
            reactant_states[node] = []
            for each in self.nodes[node].incidentNodes:
                reactant_states[each] = []

            # loop through truth table
            # find combinations of incident nodes for which the final target is activated
            for truth in self.nodes[node].table[1:]:
                if truth[-1]:

                    # find the active incident nodes
                    active_incidents = []
                    for i, each in enumerate(truth[:-1]):
                        if each:
                            active_incidents.append(self.nodes[node].table[0][i])

                    # find the active interactions in the motif
                    active_motif = []
                    for inter in self.nodes[node].motifs:
                        keep = True
                        for each in inter[0]:
                            if each not in active_incidents:
                                keep = False
                        if keep:
                            active_motif.append(inter)

                    # seed the reactant states with the active states
                    reactant_states_temp = defaultdict(list)
                    reactant_states_temp[node] = deepcopy(self.active_states[node])
                    for each in self.nodes[node].incidentNodes:
                        reactant_states_temp[each] = deepcopy(self.active_states[each])

                    # search through action information
                    for action in self.action_info[node]:
                        comp = False
                        for item in self.action_info[node]:
                            if set(action[0] + action[1]) == set(item[0]):
                                comp = True

                        # ignore if a complex
                        # binding will be explicit in the rules
                        # update the sites
                        if not comp and action[2][0][0]:  # and action[1][0] != node:
                            for j, every in enumerate(reactant_states_temp[action[1][0]]):
                                for i, thing in enumerate(action[2][0][0]):
                                    for lotsa in every:
                                        if lotsa == thing:
                                            reactant_states_temp[action[1][0]][j][thing] = action[2][0][2][i]

                    # update the reactant states
                    for each in reactant_states_temp:
                        for item in reactant_states_temp[each]:
                            if item not in reactant_states[each]:
                                reactant_states[each].append(item)

            # turn base states into a dictionary
            # todo: this will go away

            base_states = defaultdict(list)
            for each in reactant_states:
                base_states[each] = [{}]
                for i, item in enumerate(self.base_states[each][0]):

                    # todo: this 'if' statement should not be needed
                    # 'None' is getting inadvertently converted to None somewhere above
                    if not self.base_states[each][1][i]:
                        base_states[each][0][item] = 'None'
                    else:
                        base_states[each][0][item] = self.base_states[each][1][i]

            # define the target states for each node
            # ignores active states but implicitly requires activation
            # todo: make multiple target states for multiple active states
            target_states = defaultdict(list)
            for each in base_states:
                target_states[each] = [{}]

                if each != node:

                    for item in base_states[each][0]:

                        same = True
                        first = reactant_states[each][0][item]
                        for every in reactant_states[each]:
                            if every[item] != first:
                                same = False

                        base = True
                        for every in reactant_states[each]:
                            if base_states[each][0][item] != every[item]:
                                base = False

                        active = True
                        for every in reactant_states[each]:
                            for thing in self.active_states[each]:
                                if every[item] != thing[item]:
                                    active = False

                        if same and base:
                            target_states[each][0][item] = base_states[each][0][item]

                        if same and active:
                            target_states[each][0][item] = self.active_states[each][0][item]
                else:
                    # defaulting to this will work for all nodes but
                    # will allow for meaningless reactions to take place
                    for item in base_states[each][0]:
                        add_to_target = True

                        for every in self.inactive_states[each]:
                            if item in every:
                                if base_states[each][0][item] != every[item]:
                                    add_to_target = False
                        for every in self.active_states[each]:
                            if base_states[each][0][item] != every[item]:
                                add_to_target = False
                        for every in reactant_states[each]:
                            if base_states[each][0][item] != every[item]:
                                add_to_target = False

                        if add_to_target:
                            target_states[each][0][item] = base_states[each][0][item]

            # print
            # if node == 'fyn':
            #     for each in reactant_states:
            #         print each, 'base    ', base_states[each]
            #         print each, 'active  ', self.active_states[each]
            #         print each, 'inactive', self.inactive_states[each]
            #         print each, 'reactant', reactant_states[each]
            #         print each, 'target  ', target_states[each]
            #         print
                # quit()

            for interaction in self.nodes[node].motifs:
                # print
                # print interaction
                if interaction[1] and interaction not in used_interactions:
                    # print 'u1', used_interactions
                    used_interactions.append(interaction)
                    # print 'u2', used_interactions

                    active_list = []
                    for reactant in interaction[0]:
                        active_list.append(deepcopy(reactant_states[reactant]))

                    active_combos = list(product(*active_list))
                    for i, each in enumerate(active_combos):
                        active_combos[i] = list(each)
                    if not active_combos:
                        active_combos.append([{}])

                    # retrieve rxn, reaction templates, and instructions
                    current_rxn = None
                    for rxn in self.library[interaction[3][0]]:
                        if rxn.reaction == interaction[4]:
                            current_rxn = rxn

                    # loop through reactions for the current interaction
                    n = 0
                    for temp in current_rxn.rxnTemplates:
                        rxnTemp1 = re.split(r'\s*:\s*', temp)
                        rxnTemp = rxnTemp1[0]
                        rxnInstruction = current_rxn.instructions

                        # split template
                        rxn_split = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', rxnTemp)
                        rxn_split_mols = []
                        for each in rxn_split:
                            if each == 'None':
                                rxn_split_mols.append(['None', [], []])
                            else:
                                mol = each.split('(')[0]
                                sites = each.split('(')[1][:-1]
                                site_names = []
                                site_values = []
                                if sites:
                                    sites = sites.split(',')
                                    for item in sites:
                                        site_names.append(item.split('=')[0].strip())
                                        site_values.append(item.split('=')[1].strip())
                                rxn_split_mols.append([mol, site_names, site_values])

                        ops = re.findall(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', rxnTemp)
                        for i, each in enumerate(ops):
                            ops[i] = ops[i].strip()

                        rxn_split_parsed = [rxn_split_mols[0]]
                        for i, each in enumerate(ops):
                            rxn_split_parsed.append(each)
                            rxn_split_parsed.append(rxn_split_mols[i + 1])
                        # print
                        # print rxn_split_parsed


                        # Here we make substitutions from the interaction to the template.
                        # Note that the motif generation step does not differentiate between identical
                        # species from different cells. Because the substitutions relies on information
                        # from that step we must account for those identical species here.
                        # !!!!!!!!!!!!!! NEEDS SIMPLIFICATION !!!!!!!!!!!!!!!!!!!

                        # print
                        # print interaction
                        # print rxn_split_parsed
                        # quit()

                        # used_indexes = []
                        # for i, item in enumerate(rxn_split_parsed)

                        for j, each in enumerate(interaction[2]):
                            for i, item in enumerate(rxn_split_parsed):
                                if item != '>>' and item != '<>' and item != '+' and item != '%':
                                    if item[0] == each:
                                        rxn_split_parsed[i][0] = interaction[0][j]
                                    for k, every in enumerate(item[1]):
                                        if every == each:
                                            rxn_split_parsed[i][1][k] = interaction[0][j]
                                        if every[:every.rfind('_')] == each and (every[every.rfind('_') + 1:].isdigit() or every[every.rfind('_') + 1:] == 's'):
                                            rxn_split_parsed[i][1][k] = interaction[0][j] + every[every.rfind('_'):]

                        for j, each in enumerate(interaction[3]):
                            for i, item in enumerate(rxn_split_parsed):
                                if item != '>>' and item != '<>' and item != '+' and item != '%':
                                    if item[0] == each:
                                        rxn_split_parsed[i][0] = interaction[1][j]
                                    for k, every in enumerate(item[1]):
                                        if every == each:
                                            rxn_split_parsed[i][1][k] = interaction[1][j]
                                        if every[:every.rfind('_')] == each and (every[every.rfind('_') + 1:].isdigit() or every[every.rfind('_') + 1:] == 's'):
                                            rxn_split_parsed[i][1][k] = interaction[0][j] + every[every.rfind('_'):]
                        # print
                        # print node
                        # print rxn_split_parsed
                        # quit()
                        # add additional sites
                        for combo in active_combos:
                            rps = deepcopy(rxn_split_parsed)
                            combo2 = defaultdict(list)
                            reactants = []
                            targets = []
                            for i, each in enumerate(interaction[0]):
                                combo2[each] = combo[i]
                                reactants.append(each)
                            for each in interaction[1]:
                                targets.append(each)
                            monos = []
                            for each in rps:
                                if isinstance(each, list) and each[0] not in monos:
                                    monos.append(each[0])

                            # don't think I need this
                            counts = [0 for _ in range(len(monos))]
                            groups = [deepcopy(counts)]
                            for each in rps:
                                if isinstance(each, list):
                                    for i, item in enumerate(monos):
                                        if item == each[0]:
                                            groups[-1][i] += 1
                                if each == '+' or each == '>>' or each == '<>':
                                    groups.append(deepcopy(counts))

                            # print rps

                            # todo: instructions for defining states as base states

                            # routine for synthesis rxns
                            if rps[2][0] == 'None':
                                rps[0][1] = self.base_states[rps[0][0]][0]
                                rps[0][2] = self.base_states[rps[0][0]][1]

                            # routine for synthesis rxns
                            if rps[0][0] == 'None':
                                rps[2][1] = self.base_states[rps[2][0]][0]
                                rps[2][2] = self.base_states[rps[2][0]][1]

                            # add base states for newly created monomers
                            balance = [[], []]
                            side = 0
                            for i, each in enumerate(rps):
                                if each == '<>' or each == '>>':
                                    side = 1
                                if isinstance(each, list):
                                    balance[side].append(each[0])
                            new_mon = []
                            for each in balance[1]:
                                if each not in balance[0]:
                                    new_mon.append(each)
                            for i, each in enumerate(rps):
                                if isinstance(each, list) and each[0] in new_mon and each[0] != 'None':  # and not each[1]:
                                    rps[i][1] = deepcopy(self.base_states[each[0]][0])
                                    rps[i][2] = deepcopy(self.base_states[each[0]][1])

                            # print rps

                            # todo: fix this
                            for i, each in enumerate(rps):
                                if isinstance(each, list):
                                    for j, item in enumerate(each[2]):
                                        if not item:
                                            rps[i][2][j] = 'None'

                            # find removed monomers
                            balance = [[], []]
                            side = 0
                            for i, each in enumerate(rps):
                                if each == '<>' or each == '>>':
                                    side = 1
                                if isinstance(each, list):
                                    balance[side].append(each[0])
                            removed = []
                            for each in balance[0]:
                                if each not in balance[1]:
                                    removed.append(deepcopy(each))
                            # print
                            # print rps

                            # add active states for reactant monomers
                            # add additional bound monomers
                            # balance with any necessary unbound monomers
                            bonds = [0]
                            for i, each in enumerate(rps):
                                if isinstance(each,list):
                                    for item in each[2]:
                                        if item.isdigit():
                                            bonds.append(int(item))
                            bond = max(bonds)
                            rps2 = []
                            sm = []
                            for i, each in enumerate(rps):
                                rps2.append(each)
                                nm = []
                                if isinstance(each, list):
                                    if each[0] in reactants:
                                        for item in combo2[each[0]]:
                                            if item not in each[1]:
                                                rps2[-1][1].append(item)
                                                if item in self.nodes and item not in monos:
                                                    if combo2[each[0]][item] == 'None':
                                                        rps2[-1][2].append(combo2[each[0]][item])
                                                    else:
                                                        bond += 1
                                                        rps2[-1][2].append(str(bond))
                                                        nm.append('%')
                                                        nm.append([deepcopy(item), [each[0]], [str(bond)]])
                                                        if each[0] in removed:
                                                            sm.append('+')
                                                            sm.append([deepcopy(item), [each[0]], ['None']])
                                                elif item[:item.rfind('_')] in self.nodes and item[:item.rfind('_')] not in monos and item[item.rfind('_') + 1:].isdigit():
                                                    if combo2[each[0]][item] == 'None':
                                                        rps2[-1][2].append(combo2[each[0]][item])
                                                    else:
                                                        bond += 1
                                                        rps2[-1][2].append(str(bond))
                                                        nm.append('%')
                                                        nm.append([deepcopy(item), [each[0]], [str(bond)]])
                                                        if each[0] in removed:
                                                            sm.append('+')
                                                            sm.append([deepcopy(item), [each[0]], ['None']])
                                                else:
                                                    rps2[-1][2].append(combo2[each[0]][item])
                                    if each[0] in targets:
                                        for item in target_states[each[0]][0]:
                                            if item not in each[1]:
                                                rps2[-1][1].append(item)
                                                if item in self.nodes:
                                                    if target_states[each[0]][0][item] == 'None':
                                                        rps2[-1][2].append(target_states[each[0]][0][item])
                                                    else:
                                                        bond += 1
                                                        rps2[-1][2].append(str(bond))
                                                        nm.append('%')
                                                        nm.append([deepcopy(item), [each[0]], [str(bond)]])
                                                        if each[0] in removed:
                                                            sm.append('+')
                                                            sm.append([deepcopy(item), [each[0]], ['None']])
                                                else:
                                                    rps2[-1][2].append(target_states[each[0]][0][item])
                                rps2.extend(nm)
                            rps2.extend(sm)
                            rps = deepcopy(rps2)
                            # print rps

                            # update for sequential reaction info
                            only_target = True
                            for each in rps:
                                if isinstance(each, list) and each[0] != interaction[1][0]:
                                    only_target = False

                            if not only_target:

                                # parse reaction sequence info from the label file
                                for sequence in self.nodes[node].sequence:
                                    seq = []
                                    for item in sequence:
                                        if '|' in item:
                                            seq.append(item.split('|'))
                                        else:
                                            seq.append([item])
                                    seqs = list(product(*seq))
                                    for i, item in enumerate(seqs):
                                        seqs[i] = list(item)
                                    for i, each in enumerate(seqs):
                                        for j, item in enumerate(each):
                                            seqs[i][j] = [item]
                                    for i, each in enumerate(seqs):
                                        for j, item in enumerate(each):
                                            if ':' in item[0]:
                                                seqs[i][j] = item[0].split(':')

                                    # incorporate reaction sequence info
                                    for each in seqs:

                                        # for each reactant in sequence info find its position in the sequence
                                        seq_ind = 0
                                        for i, item in enumerate(each):
                                            if set(item) == set(interaction[0]):
                                                seq_ind = i

                                        seq2 = each[:seq_ind]  # all reactants before the reactant in question
                                        for item in seq2:

                                            on_off = None
                                            if item[0] == 'not':
                                                on_off = 0
                                                item = item[1:]
                                            else:
                                                on_off = 1

                                            # find activity of previous reactants
                                            action = None
                                            for every in self.action_info[node]:
                                                if set(item) == set(every[0]) and interaction[1] == every[1]:
                                                    action = every

                                            for i, every in enumerate(rps2):
                                                if isinstance(every, list) and every[0] == interaction[1][0]:
                                                    for j, thing in enumerate(action[2][0][0]):
                                                        if on_off:
                                                            if thing in rps[i][1]:
                                                                ind = rps[i][1].index(thing)
                                                                if thing in self.monomer_info:
                                                                    if action[2][0][2][j] == 'None':
                                                                        rps2[i][2][ind] = action[2][0][2][j]
                                                                    else:
                                                                        bond += 1
                                                                        rps2[i][2][ind] = str(bond)
                                                                else:
                                                                    rps2[i][2][ind] = action[2][0][2][j]
                                                            else:
                                                                rps2[i][1].append(thing)
                                                                if thing in self.monomer_info:
                                                                    if action[2][0][2][j] == 'None':
                                                                        rps2[i][2].append(action[2][0][2][j])
                                                                    else:
                                                                        bond += 1
                                                                        rps2[i][2].append(str(bond))
                                                                else:
                                                                    rps2[i][2].append(action[2][0][2][j])
                                                        else:
                                                            if thing in rps[i][1]:
                                                                ind = rps[i][1].index(thing)
                                                                if thing in self.monomer_info:
                                                                    rps2[i][2][ind] = 'None'
                                                                else:
                                                                    rps2[i][2][ind] = action[2][0][1][j]
                                                            else:
                                                                rps2[i][1].append(thing)
                                                                if thing in self.monomer_info:
                                                                    rps2[i][2].append('None')
                                                                else:
                                                                    rps2[i][2].append(action[2][0][1][j])

                            rps3 = []
                            for i, each in enumerate(rps2):
                                rps3.append(each)
                                if isinstance(each, list) and len(each[1]) != len(rps[i][1]):
                                    diff = list(set(each[1])-set(rps[i][1]))
                                    for item in diff:
                                        if item in self.monomer_info and each[2][each[1].index(item)] != 'None':
                                            rps3.append('%')
                                            rps3.append([item, [each[0]], [str(each[2][each[1].index(item)])]])

                            rps = rps3

                            for item in self.nodes[node].compete:
                                competitors = item.split(':')
                                site = ''
                                for every in competitors:
                                    site += every + '_'
                                site = site[:-1]
                                for j, every in enumerate(rps):
                                    if isinstance(every, list) and every[0] == node:
                                        for k, thing in enumerate(every[1]):
                                            if thing in competitors:
                                                rps[j][1][k] = site

                            # print rps

                            # substitute in integers, None, ANY, and WILD
                            for i, item in enumerate(rps):
                                if item != '+' and item != '%' and item != '>>' and item != '<>':
                                    for j, every in enumerate(item[2]):
                                        if item[1][j][:item[1][j].rfind('_')] in self.nodes and every != 'None' and every != 'ANY' and every != 'WILD':
                                            rps[i][2][j] = int(rps[i][2][j])
                                        if item[1][j] in self.nodes and every != 'None' and every != 'ANY' and every != 'WILD':
                                            rps[i][2][j] = int(rps[i][2][j])
                                        if every == 'None':
                                            rps[i][2][j] = None
                                        if every == 'ANY':
                                            rps[i][2][j] = ANY
                                        if every == 'WILD':
                                            rps[i][2][j] = WILD

                            # print rps

                            # define the rule rule_name
                            rule_name = ''
                            for item in interaction[0]:
                                if item:
                                    rule_name += item + '_'
                                else:
                                    rule_name = rule_name[:-1]
                            rule_name += interaction[4] + '_'
                            for item in interaction[1]:
                                if item:
                                    rule_name += item + '_'
                            rule_name += str(n)

                            # define monomer patterns
                            mon_pats = []

                            for item in rps:
                                if item == '+' or item == '%' or item == '>>' or item == '<>':
                                    mon_pats.append(item)
                                else:
                                    if item[0] == 'None':
                                        mon_pats.append('None')
                                    else:
                                        mon_states = {}
                                        for i, every in enumerate(item[1]):
                                            mon_states[every] = item[2][i]
                                        # print item[0]
                                        mon_obj = self.model.monomers[item[0]]
                                        mon_pats.append(MonomerPattern(mon_obj, mon_states, None))

                            # print
                            # for each in mon_pats:
                            #     print each
                            # print mon_pats

                            # define complex patterns
                            com_pats_temp = [[]]
                            for item in mon_pats:
                                if item == '>>' or item == '<>':
                                    com_pats_temp.extend([item, []])
                                elif item == '+':
                                    com_pats_temp.append([])
                                elif item == '%':
                                    pass
                                else:
                                    com_pats_temp[-1].append(item)
                            com_pats = []
                            for item in com_pats_temp:
                                if item == '>>' or item == '<>':
                                    com_pats.append(item)
                                elif item == ['None']:
                                    pass
                                else:
                                    com_pats.append(ComplexPattern(item, None))

                            # define reversibility and split patterns into reactants and products
                            react_com_pats = []
                            prod_com_pats = []
                            carrot = 0
                            reversible = None
                            for item in com_pats:
                                if item == '<>':
                                    carrot = 1
                                    reversible = True
                                elif item == '>>':
                                    carrot = 1
                                    reversible = False
                                else:
                                    if carrot == 0:
                                        react_com_pats.append(item)
                                    if carrot == 1:
                                        prod_com_pats.append(item)
                            order = [len(react_com_pats), len(prod_com_pats)]

                            # define rule expression
                            rule_exp = RuleExpression(ReactionPattern(react_com_pats), ReactionPattern(prod_com_pats), reversible)
                            if reversible:
                                forward = rule_name + '_' + str(order[0]) + 'kf'
                                self.parameter(forward, 1)
                                reverse = rule_name + '_' + str(order[1]) + 'kr'
                                self.parameter(reverse, 1)
                                self.rule(rule_name, rule_exp, self.model.parameters[forward],
                                          self.model.parameters[reverse])
                            else:
                                forward = rule_name + '_' + str(order[0]) + 'kc'
                                self.parameter(forward, 1)
                                self.rule(rule_name, rule_exp, self.model.parameters[forward])

                            n += 1
                else:
                    pass
                    # print 'used------------------------------------'
    def _add_initials(self):

        for each in self.monomer_info:

            sites = []
            states = []

            for item in self.monomer_info[each][0]:
                sites.append(item)
            for item in sites:
                if item not in self.monomer_info[each][1]:
                    states.append([None])
                else:
                    states.append(self.monomer_info[each][1][item])

            combos = list(product(*states))
            for i, item in enumerate(combos):
                combos[i] = list(item)
            for i, every in enumerate(combos):
                init_name = each + '_' + str(i) + '_0'
                mon_obj = self.model.monomers[each]
                state = {}
                for j, item in enumerate(every):
                    state[sites[j]] = item
                self.parameter(init_name, 0)
                self.initial(MonomerPattern(mon_obj, state, None), self.model.parameters[init_name])

    def _add_observables(self):

        # add all active states combinations to observables
        for each in self.active_states:
            for i, item in enumerate(self.active_states[each]):

                obs_name = each + '_' + str(i) + '_obs'

                mon_states = {}
                for every in item:
                    for thing in self.action_info[each]:
                        for lotsa in thing[2][0][0]:
                            if lotsa == every:
                                mon_states[every] = item[every]
                for every in mon_states:
                    if mon_states[every].isdigit():
                        mon_states[every] = ANY
                    if mon_states[every] == 'None':
                        mon_states[every] = None
                mon_pat = MonomerPattern(self.model.monomers[each], mon_states, None)

                self.observable(obs_name, mon_pat)

        # adds all possible site combinations to observables
        # for each in self.monomer_info:
        #
        #     sites = []
        #     states = []
        #
        #     for item in self.monomer_info[each][0]:
        #         sites.append(item)
        #     for item in sites:
        #         if item not in self.monomer_info[each][1]:
        #             states.append([None, ANY])
        #         else:
        #             states.append(self.monomer_info[each][1][item])
        #
        #     combos = list(product(*states))
        #     for i,item in enumerate(combos):
        #         combos[i] = list(item)
        #     for i,every in enumerate(combos):
        #         obs_name = each + '_' + str(i) + '_obs'
        #         mon_states = {}
        #         for j, thing in enumerate(sites):
        #             mon_states[thing] = every[j]
        #         pat = MonomerPattern(self.model.monomers[each], mon_states, None)
        #         self.observable(obs_name, pat)

        # # adds only monomers to observables
        for each in self.monomer_info:
            obs_name = each + '_obs'
            self.observable(obs_name, self.model.monomers[each])
