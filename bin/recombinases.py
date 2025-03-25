# pylint: disable=R0916

""" Recombinase rules and aliases """

import itertools as it

from dataclasses import dataclass


MGE_ALIASES = {
    "c1_n1ser": "ser_tn",
    "c2_n1ser": "ser_ce",
    "c3_n1ser": "ser_lsr",
    "casposons": "cas1",
}


@dataclass
class MgeRule:
    '''The following class defines the set of rules used to determine MGE type.
    Type classification is based on two criteria:
    1. Recombinase subfamily i.e.
    2. Structural information

    MGE categories include
    - IS_Tn(tn)
    - Phage(ph)
    - CE(conjugative elements)
    - Integron(int)
    - Cellular(cell)
    As well as newly introduced categories:
    - Phage_like(pli)
    - Mobility island(mi)'''
    subfamily: str = None
    is_tn: bool = False
    phage: bool = False
    ce: bool = False
    integron: bool = False
    cellular: bool = False
    recombinase_scan: bool = False

    def __post_init__(self):
        """ Deal with special case for Tn3
        since it can carry conjugative system
        - ignored if the rule is used during recombinase-scans
        """
        if all(
                (
                        self.subfamily is not None,
                        "tn3" in self.subfamily.lower(),
                        not self.recombinase_scan,
                )
        ):
            self.ce = 1

    def get_signals(self):
        """ Returns MGE signals of rule. """
        return tuple(
            k
            for k, v in self.__dict__.items()
            if v and k not in ("subfamily", "recombinase_scan")
        )

    def c_tn_check(self, island):
        """ Tn check. """
        # c_tn, n_recombinases = island.c_tn, len(island.recombinases)
        c_tn, n_recombinases = island.c_tn, sum(island.recombinases.values())
        if self.is_tn and not self.cellular and not self.ce and not self.phage:
            # IS_Tn
            c_tn += 1
        elif (
                self.is_tn and self.ce and island.conj_man_count < 1 and
                n_recombinases == 2 and not island.tn3_found and not island.ser_found
        ):
            # c2_n1ser(considers solo c2_n1ser and Tn3)
            c_tn = 1
        elif self.is_tn and self.ce and island.conj_man_count < 1 and n_recombinases == 1:
            # c2_n1ser(considers c2_n1ser and Tn3 as one tn when not together)
            c_tn += 1

        # disentangles recombinase shared by tn and ph
        if (self.is_tn and self.phage and island.phage_count < 2):
            c_tn += 1

        return c_tn

    def patch_c_tn_check(self, island):
        """Deals with special case when 2 recombinases.

        c_tn = is_tn and ce and conj_man_count < 1 and |recombinases|=2 and !(tn3 or ser)
        """
        # old check was:
        # if rule.is_tn and rule.ce and self.conj_man_count < 1
        # and len(self.recombinases) == 2 and not self.tn3_found:
        # 	self.c_tn = True
        # elif rule.is_tn and rule.ce and self.conj_man_count < 1
        # and len(self.recombinases) == 2 and not self.ser_found:
        # 	self.c_tn = True
        # old:
        # if len(island.recombinases) == 2 and self.is_tn and self.ce and island.conj_man_count < 1:
        #     recombinase_types = list(island.recombinases)
        #     two_tn3 = "tn3" in recombinase_types[0] and "tn3" in recombinase_types[1]
        #     two_ser_ce = "ser_ce" in recombinase_types[0] and "ser_ce" in recombinase_types[1]

        #     mixed = "tn3" in recombinase_types[0] or "tn3" in recombinase_types[1]
        #     mixed |= "ser_ce" in recombinase_types[0] or "ser_ce" in recombinase_types[1]

        #     return two_tn3 != two_ser_ce or mixed

        if sum(island.recombinases.values()) == 2 and self.is_tn and self.ce and island.conj_man_count < 1:
            # recombinase_types = ",".join(list(island.recombinases))
            recombinase_types = ",".join(it.chain(*it.chain((r,) * c for r, c in island.recombinases.items())))

            mixed = "tn3" in recombinase_types and "ser_ce" in recombinase_types

            two_tn3 = recombinase_types.count("tn3") == 2
            two_ser_ce = recombinase_types.count("ser_ce") == 2

            return (two_tn3 != two_ser_ce) or mixed

        return False

    def phage_check(self, island):
        """Phage annotation based on recombinase presence
        and phage structural genes in the neighbourhood.
        """
        phage, c_mi, nov = island.phage, island.c_mi, island.nov
        phage |= (self.is_tn and self.phage)
        phage |= (self.phage and not self.ce)
        phage |= (self.phage and self.ce)

        c_mi |= (not self.phage and self.ce)
        nov = c_mi

        return phage, c_mi, nov

    def phage_like_check(self, island, is_brujita):
        """Annotate phage_like element
        (presence of phage specific recombinase and absence of phage structural genes
        in the neighbourhood)
        and mobility island
        (presence of recombinase common to phages and conjugative elements
        and absence of phage structural genes in the neighbourhood)
        """
        c_pli, c_mi = island.c_pli, island.c_mi
        if not self.is_tn:
            if self.phage and (not self.ce or is_brujita):
                c_pli = 1
            elif self.ce and (not self.phage or not is_brujita):
                c_mi = 1

        return c_pli, c_mi

    def conjug_element_check(self, island):
        """Conjugative element annotation based on presence of recombinase
        and presence of conjugative machinery genes in the neighbourhood
        """
        c_ce, nov = island.c_ce, island.nov
        if self.ce:
            c_ce = 1
        elif self.phage:
            c_ce = nov = 1
        elif all(
                (
                        bool(self.is_tn),
                        bool(self.ce),
                        # len(island.recombinases) >= 3,
                        sum(island.recombinases.values()) >= 3,
                        (island.tn3_found or island.ser_found)
                )
        ):
            c_ce = 1

        return c_ce, nov

    def mobility_island_check(self, island):
        """Annotate MI(Mobility island) presence of both
        phage structural genes and conjugation machinery genes in the neighbourhood
        """
        phage, c_mi, nov = island.phage, island.c_mi, island.nov
        if self.is_tn and self.phage:
            phage = 1
        else:
            c_mi = nov = 1

        return phage, c_mi, nov
