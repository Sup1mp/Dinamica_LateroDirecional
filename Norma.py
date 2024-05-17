

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
        3.3.1.1
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
        3.3.1.2
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
        3.3.1.3
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
        3.3.1.4
        '''
        lvl = 0
        pa = Crs*Wnrs
        if pa >= 0.5:
            lvl = 1
        elif pa >= 0.3:
            lvl = 2
        elif pa >= 0.15:
            lvl = 3
        
        return lvl
            