import lcapy as lc
import sympy as sym
 
class RFCircuit:
    
    def __init__(self):
        
        self._circuit = None
        self.Z = None
        self.Gamma = None
        #self.deltaGamma = None
        #self.Reff = None
        self.wres = None
        
    def set_circuit(self, lcapy_circ):
        # save input lcapy circuit e.g.: lc.R('R') + lc.C('C') | lc.L('L')
        self._circuit = lcapy_circ
        
        # Get Z
        self.Z = self._circuit.Z.expr #get sympy expr
        self.args = {str(elem):elem for elem in self.Z.free_symbols} #get symbols involved in expression
        #replace s with Iw
        w=sym.Symbol('w', real=True)
        self.Z = self.Z.subs({self.args['s']:sym.I*w}).simplify()
        self.args = {str(elem):elem for elem in self.Z.free_symbols} #get symbols involved in expression
        
        # Get gamma
        #Z0=sym.Symbol('Z0', real=True)
        self.Gamma = (self.Z-50)/(self.Z+50)
        
        # get resonant frequency by solving Im(Z)=0
        self.wres = sym.solve(sym.im(self.Z), self.args['w'])
       
    def draw_circuit(self, filename=None, style='british', spacing=1.5, scale=2):
        if self._circuit:
            self._circuit.draw(filename=filename, svg=True, style=style, node_spacing=spacing, scale=scale)
            
    def get_Gamma(self, args={}):
        # args are dict with symbols and values to replace
        Gamma=self.Gamma.subs(args)
        return sym.lambdify(Gamma.free_symbols,Gamma)
    
    def get_wres(self, wres_expr, args={}):
        # args are dict with symbols and values to replace
        wres=wres_expr.subs(args)
        return sym.lambdify(wres.free_symbols,wres)
    
    def get_DeltaGamma(self, dX, args={}):
        # args are dict with symbols and values to replace
        self.DeltaGamma = sym.diff(self.Gamma, dX).simplify()
        DeltaGamma = self.DeltaGamma.subs(args)
        return sym.lambdify(DeltaGamma.free_symbols,DeltaGamma)
    
    def get_Reff(self, wres_expr, args={}):
        # put resonance frequency into real part of Z
        self.Reff = sym.re(self.Z).subs(self.args['w'], wres_expr).simplify()
        Reff=self.Reff.subs(args)
        return sym.lambdify(Reff.free_symbols,Reff)
        