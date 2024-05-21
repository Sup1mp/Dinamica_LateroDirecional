

class MILF8587C:
    def __init__(self, Class: int, Category: str):
        '''
        Requisitos da norma MIL-F-8587C
            Class 1 :   Small, light airplanes such as\n
                        Light utility\n
                        Primary trainer\n
                        Light observation\n\n

            Class 2:   Medium weight, low-to-medium maneuverability airplanes such as\n
                        Heavy utility/search and rescue\n
                        Light or medium transport/cargo/tanker\n
                        Early warning/electronic countermeasures/airborne command,\
                        control, or communications relay\n
                        Antisubmarine\n
                        Assault transport\n
                        Reconnaissance\n
                        Tactical bomber\n
                        Heavy attack\n
                        Trainer for Class II\n\n

            Class 3:    Large, heavy, low-to-medium maneuverability airplanes such as\n
                        Heavy transport/cargo/tanker\n
                        Heavy bomber\n
                        Patrol/early warning/electronic countermeasures/airborne command,\
                        control, or communications relay\n
                        Trainer for Class III\n\n

            Class 4:    High-maneuverability airplanes such as\n
                        Fighter/interceptor\n
                        Attack\n
                        Tactical reconnaissance\n
                        Observation\n
                        Trainer for Class IV\n
            
            Category A: Those nonterminal Flight Phases that require rapid maneuvering, precision\
                        tracking, or precise flight-path control.\n

            Category B: Those nonterminal Flight Phases that are normally accomplished using gradual\
                        maneuvers and without precision tracking, although accurate flight-path control\
                        may be required.\n

            Category C: Terminal Flight Phases are normally accomplished using gradual maneuvers and\
                        usually require accurate flight-path control.
        '''

        self.Cl = Class
        self.Ca = Category.upper()

        return
    
    def dutch_roll (self, Wnd, Cd):
        '''
        Topic 3.3.1.1
            Wnd : natural frequency rate
            Cd : damping ratio
        '''
        lvl = 0
        match self.Ca:
            case "A":
                if Cd >= 0.19 and Cd*Wnd >= 0.35:
                    if self.Cl == 1 or self.Cl == 4:
                        if Wnd >= 1:
                            lvl = 1
                    else:
                        if Wnd >= 0.4:
                            lvl = 1

            case "B":
                if Cd >= 0.08 and Cd*Wnd >= 0.15 and Wnd >= 0.4:
                    lvl = 1

            case "C":
                if Cd >= 0.08:
                    if self.Cl == 1 or self.Cl == 4:
                        if Cd*Wnd >= 0.15 and Wnd >= 1:
                            lvl = 1
                    elif self.Cl == 2 or self.Cl == 3:
                        if Cd*Wnd >= 0.1 and Wnd >= 0.4:
                            lvl = 1
        if lvl == 0:
            if Cd >= 0.02 and Cd*Wnd >= 0.05 and Wnd >= 0.4:
                lvl = 2
            elif Cd >= 0 and Cd*Wnd >= 0 and Wnd >= 0.4:
                lvl = 3

        return lvl
    
    def roll_mode (self, t_r):
        '''
        Topic 3.3.1.2
            t_r : roll mode time constant
        '''
        lvl = 0
        match self.Ca:
            case "A":
                if self.Cl == 1 or self.Cl == 4:
                    if t_r <= 1:
                        lvl = 1
                    elif t_r <= 1.4:
                        lvl = 2

                if self.Cl == 2 or self.Cl == 3:
                    if t_r <= 1.4:
                        lvl = 1
                    elif t_r <= 3:
                        lvl = 2
            case "B":
                if t_r <= 1.4:
                    lvl = 1
                elif t_r <= 3:
                    lvl = 2
            
            case "C":
                if self.Cl == 1 or self.Cl == 4:
                    if t_r <= 1:
                        lvl = 1
                    elif t_r <= 1.4:
                        lvl = 2

                if self.Cl == 2 or self.Cl == 3:
                    if t_r <= 1.4:
                        lvl = 1
                    elif t_r <= 3:
                        lvl = 2
        if lvl == 0:
            if t_r <= 10:
                lvl = 3
        
        return lvl
    
    def spiral_stability (self, t):
        '''
        Topic 3.3.1.3
            t : time for the bank angle to double
        '''
        lvl = 0
        if (t >= 12 and self.Ca in "AC") or (t >= 20 and self.Ca == "B"):
            lvl = 1
        elif t >= 8:
            lvl = 2
        elif t >= 4:
            lvl = 3

        return lvl
    
    def roll_spiral (self, Crs, Wnrs):
        '''
        Topic 3.3.1.4
            Wnrs : natural frequency of roll spiral
            Crs : damping ration of roll spiral
        '''
        lvl = 0
        if self.Ca == "B" or self.Ca == "C":
            pa = Crs*Wnrs
            if pa >= 0.5:
                lvl = 1
            elif pa >= 0.3:
                lvl = 2
            elif pa >= 0.15:
                lvl = 3
        
        return lvl
    
    def roll_rate (self, psi, por):
        '''
        Topic 3.3.2.2
            psi : roll rate at the first minimum following the first peak
            por : percentage of roll rate at first peak
        '''
        lvl = 0
        if psi*por > 0:     # checks if is same sign
            p = round(100*(psi/por), 2)
            if self.Ca == "A" or self.Ca == "C":
                if p > 60:
                    lvl = 1
                elif p > 25:
                    lvl = 2
            else:
                if p > 25:
                    lvl = 1
                elif p > 0:
                    lvl = 2
        return lvl
    
    def sideslip (self, dbeta, k):
        '''
        Topic 3.3.2.4
            dbeta : sideslip increment
            k : parameter k (6.2.6)
        '''
        lvl = 0
        r = dbeta/k     # ratio of dbeta and k
        if self.Ca == "A":
            if r <= 6:
                lvl = 1
            elif r <= 15:
                lvl = 2
        else:
            if r <= 10:
                lvl = 1
            elif r <= 15:
                lvl = 2

        return lvl
    
    def roll_control (self, t, speed_range = "", lv = 1):
        '''
        Topic 3.3.4 + 3.3.4.1 + 3.3.4.2
            t : time to achieve the given banking angle
            speed_range : the speed range
            lv : the level for the speed range
        '''
        lvl = 0
        match self.Cl:
            case 1:
                match self.Ca:
                    case "A":
                        # to achieve 60°
                        if t <= 1.3:
                            lvl = 1
                        elif t <= 1.7:
                            lvl = 2
                        elif t <= 2.6:
                            lvl = 3

                    case "B":
                        # to achieve 60°
                        if t <= 1.7:
                            lvl = 1
                        elif t <= 2.5:
                            lvl = 2
                        elif t <= 3.4:
                            lvl = 3

                    case "C":
                        # to achieve 30°
                        if t <= 1.3:
                            lvl = 1
                        elif t <= 1.8:
                            lvl = 2
                        elif t <= 2.6:
                            lvl = 3

            case 2:
                match self.Ca:
                    case "A":
                        # to achieve 45°
                        if t <= 1.4:
                            lvl = 1
                        elif t <= 1.9:
                            lvl = 2
                        elif t <= 2.8:
                            lvl = 3

                    case "B":
                        # to achieve 45°
                        if t <= 1.9:
                            lvl = 1
                        elif t <= 2.8:
                            lvl = 2
                        elif t <= 3.8:
                            lvl = 3

                    case "C":
                        # to achieve 30°
                        if t <= 1.8:
                            lvl = 1
                        elif t <= 2.5:
                            lvl = 2
                        elif t <= 3.6:
                            lvl = 3

            case 3:
                # to achieve 30°
                match self.Ca:
                    case "A":
                        match speed_range:
                            case "L":
                                if t <= 1.8:
                                    lvl = 1
                                elif t <= 2.4 and lv > 1:
                                    lvl = 2

                            case "M":
                                if t <= 1.5:
                                    lvl = 1
                                elif t <= 2.0 and lv > 1:
                                    lvl = 2

                            case "H":
                                if t <= 2.0:
                                    lvl = 1
                                elif t <= 2.5 and lv > 1:
                                    lvl = 2

                        if lvl == 0 and lv > 1:
                            if t <= 3.0:
                                lvl = 3

                    case "B":
                        match speed_range:
                            case "L":
                                if t <= 2.3:
                                    lvl = 1
                                elif t <= 3.9 and lv > 1:
                                    lvl = 2

                            case "M":
                                if t <= 2.0:
                                    lvl = 1
                                elif t <= 3.3 and lv > 1:
                                    lvl = 2

                            case "H":
                                if t <= 2.3:
                                    lvl = 1
                                elif t <= 3.9 and lv > 1:
                                    lvl = 2

                        if lvl == 0 and lv > 1:
                            if t <= 5.0:
                                lvl = 3

                    case "C":
                        if t <= 2.5:
                            lvl = 1
                        elif t <= 4.0 and lv > 1:
                            lvl = 2
                        elif t <= 6.0 and lv > 1:
                            lvl = 3

            case 4:
                match self.Ca:
                    case "A":
                        match speed_range:
                            case "VL":
                                # to achieve 30°
                                if t <= 1.1:
                                    lvl = 1
                                elif t <= 1.6 and lv > 1:
                                    lvl = 2
                                elif t <= 2.6 and lv > 1:
                                    lvl = 3

                            case "L":
                                # to achieve 30°
                                if t <= 1.1:
                                    lvl = 1
                                elif t <= 1.5 and lv > 1:
                                    lvl = 2
                                elif t <= 2.0 and lv > 1:
                                    lvl = 3
                            
                            case "M":
                                # to achieve 90°
                                if t <= 1.3:
                                    lvl = 1
                                elif t <= 1.7 and lv > 1:
                                    lvl = 2
                                elif t <= 2.6 and lv > 1:
                                    lvl = 3
                            
                            case "H":
                                # to achieve 50°
                                if t <= 1.1:
                                    lvl = 1
                                elif t <= 1.3 and lv > 1:
                                    lvl = 2
                                elif t <= 2.6 and lv > 1:
                                    lvl = 3

                    case "B":
                        # to achieve 90°
                        match speed_range:
                            case "VL":
                                if t <= 2.0:
                                    lvl = 1
                                elif t <= 2.8 and lv > 1:
                                    lvl = 2
                                elif t <= 3.7 and lv > 1:
                                    lvl = 3

                            case "L":
                                if t <= 1.7:
                                    lvl = 1
                                elif t <= 2.5 and lv > 1:
                                    lvl = 2
                                elif t <= 3.4 and lv > 1:
                                    lvl = 3
                            
                            case "M":
                                if t <= 1.7:
                                    lvl = 1
                                elif t <= 2.5 and lv > 1:
                                    lvl = 2
                                elif t <= 3.4 and lv > 1:
                                    lvl = 3
                            
                            case "H":
                                if t <= 1.7:
                                    lvl = 1
                                elif t <= 2.5 and lv > 1:
                                    lvl = 2
                                elif t <= 3.4 and lv > 1:
                                    lvl = 3

                    case "C":
                        # to achieve 30°
                        if t <= 1.1:
                            lvl = 1
                        elif t <= 1.3 and lv > 1:
                            lvl = 2
                        elif t <= 2.0 and lv > 1:
                            lvl = 3

        return lvl
    
    def speed_range (self, Vomin, Vmin, Vomax, Vmax, V, lv):
        '''
        From topic 3.3.4.1 + 3.3.4.2
        '''
        sp = ""
        if self.Cl == 3:
            if lv == 1:
                if Vomin <= V < 1.8*Vmin:
                    sp = "L"
                elif max(1.8*Vmin, Vomin) <= V < min(0.7*Vmax, Vomax):
                    sp = "M"
                elif min(0.7*Vmax, Vomax) <= V <= Vomax:
                    sp = "H"
            
            else:
                if Vmin <= V < 1.8*Vmin:
                    sp = "L"
                elif 1.8*Vmin <= V < 0.7*Vmax:
                    sp = "M"
                elif 0.7*Vmax <= V <= Vmax:
                    sp = "H"

        elif self.Cl == 4:
            if lv == 1:
                if Vomin <= V <= Vmin + 10.28:
                    sp = "VL"
                elif max(Vomin, Vmin + 10.28) <= V < 1.4*Vmin:
                    sp = "L"
                elif 1.4*Vomin <= V < min(0.7*Vmax, Vomax):
                    sp = "M"
                elif min(0.7*Vmax, Vomax) <= V <= Vomax:
                    sp = "H"
            
            else:
                if Vmin <= V <= Vmin + 10.28:
                    sp = "VL"
                elif Vmin + 10.28 <= V < 1.4*Vmin:
                    sp = "L"
                elif 1.4*Vmin <= V < 0.7*Vmax:
                    sp = "M"
                elif 0.7*Vmax <= V <= Vmax:
                    sp = "H"

        return sp
               
    if __name__ == "__main__":
        import main
            