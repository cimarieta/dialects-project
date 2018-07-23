# UnitTest general_model
import random
import numpy as np
import unittest

from general_model import DialectDynamicsSimulator, integrate_peer_ode


class BasicSetup(object):
    def __init__(self, conformity_mode='peer', lpub=1., n=3):
        self.params = {
            'pi_A': 0.8,
            'pi_B': 0.2,
            'pa': 0.5,
            'pb': 0.2,
            'xi': 10,
            'lpub': lpub,
            'n': n
        }
        self.params['lpriv'] = 1 - self.params['lpub']
        self.simulator = DialectDynamicsSimulator(self.params, conformity_mode=conformity_mode)
        self.init_vec = [random.random() for k in range(self.params['n'] + 1)]
        self.init_vec = [elem/sum(self.init_vec) for elem in self.init_vec]
        self.init_vec = self.init_vec[:-1]


class TestDialectDynamicsN3(unittest.TestCase):
    """
        Test dynamics when n = 3
    """
    def common_setup(self, conformity_mode='peer', lpub=1.):
        setup = BasicSetup(conformity_mode=conformity_mode, lpub=lpub)
        simulator, init_vec, params = setup.simulator, setup.init_vec, setup.params
        self.dz_dt = simulator.public_private_dz(init_vec)
        self.dz_results = self.calc_results(params, init_vec)

    def calc_results(self, params, init_vec):
        raise NotImplementedError('calc_results is not implemented')

    def check_all(self):
        self.assertAlmostEqual(self.dz_dt[0], self.dz_results[0])
        self.assertAlmostEqual(self.dz_dt[1], self.dz_results[1])
        self.assertAlmostEqual(self.dz_dt[2], self.dz_results[2])
        #self.assertAlmostEqual(self.dz_dt[3], self.dz_results[3])


class TestPublicDynamicsN3(TestDialectDynamicsN3):
    """
        Test public dynamics when n = 3
    """
    def setUp(self):
        self.common_setup(lpub=1.)

    def calc_results(self, params, init_vec):
        pi_A, pi_B, xi = params['pi_A'], params['pi_B'], params['xi']
        zA = (init_vec[1] + init_vec[2]*2 + (1-sum(init_vec))*3)/3
        return [
            -3*(pi_A + xi*zA)*init_vec[0] + (pi_B + xi*(1-zA))*init_vec[1],
            -(2*pi_A + pi_B + xi*(1+zA))*init_vec[1] + 3*(pi_A +
                xi*zA)*init_vec[0] + 2*(pi_B + xi*(1-zA))*init_vec[2],
            -(pi_A + 2*pi_B + xi*(2-zA))*init_vec[2] + 2*(pi_A +
                xi*zA)*init_vec[1] + 3*(pi_B + xi*(1-zA))*(1-sum(init_vec)),
            -3*(pi_B + xi*(1-zA))*(1-sum(init_vec)) + (pi_A + xi*zA)*init_vec[2]
        ]

    def test_dynamics(self):
        self.check_all()


class TestPeerPrivateDynamicsN3(TestDialectDynamicsN3):
    """
        Test private dynamics when n = 3 in the peer pressure conformity mode
    """
    def setUp(self):
        self.common_setup(conformity_mode='peer', lpub=0.)

    def calc_results(self, params, init_vec):
        pa, pb, xi = params['pa'], params['pb'], params['xi']
        return [
            -3*pa*init_vec[0] + (pb + xi)*init_vec[1],
            -(2*pa + pb + 2*xi)*init_vec[1] + 3*pa*init_vec[0] + (2*pb + xi)*init_vec[2],
            -(pa + 2*pb + 2*xi)*init_vec[2] + (2*pa + xi)*init_vec[1] + 3*pb*(1-sum(init_vec)),
            -3*pb*(1-sum(init_vec)) + (pa + xi)*init_vec[2]
        ]

    def test_dynamics(self):
        self.check_all()


class TestUnanimityPrivateDynamicsN3(TestDialectDynamicsN3):
    """
        Test private dynamics when n = 3 in the peer pressure conformity mode
    """
    def setUp(self):
        self.common_setup(conformity_mode='unanimity', lpub=0.)

    def calc_results(self, params, init_vec):
        pa, pb, xi = params['pa'], params['pb'], params['xi']
        return [
            -3*pa*init_vec[0] + (pb + xi)*init_vec[1],
            -(2*pa + pb + xi)*init_vec[1] + 3*pa*init_vec[0] + 2*pb*init_vec[2],
            -(pa + 2*pb + xi)*init_vec[2] + 2*pa*init_vec[1] + 3*pb*(1-sum(init_vec)),
            -3*pb*(1-sum(init_vec)) + (pa + xi)*init_vec[2]
        ]

    def test_dynamics(self):
        self.check_all()


#class TestSanityDialectDynamics(unittest.TestCase):
#    """
#        Sanity check (only works when using all the dimensions of dz/dt)
#    """
#    def setUp(self):
#        lpub = random.random()
#        setup = BasicSetup(conformity_mode='peer', lpub=lpub, n=random.randint(2, 11))
#        simulator, init_vec = setup.simulator, setup.init_vec
#        self.dz_dt = simulator.public_private_dz(init_vec)
#
#    def test_sum_diffs(self):
#        self.assertAlmostEqual(sum(self.dz_dt), 0)


class TestIntegrationPeerDynamics(unittest.TestCase):
    """
        Test zA final fraction
    """
    def setUp(self):
        self.params = {
            'pi_A': 0.8,
            'pi_B': 0.2,
            'pa': 0.5,
            'pb': 0.6,
            'xi': 10,
            'lpub': 0.,
            'n': 3
        }
        self.params['lpriv'] = 1 - self.params['lpub']
        self.simulator = DialectDynamicsSimulator(self.params, conformity_mode='peer')
        self.params_order_fun = ['pi_A', 'pi_B', 'xi', 'lpub', 'lpriv', 'pa', 'pb', 'n']

    def test_peer_integrate(self):
        params_nums = tuple([self.params[lab] for lab in self.params_order_fun])
        curva, lista_t = integrate_peer_ode([0.99, 0., 0.], params_nums)
        curva_transp = np.asarray(curva).T.tolist()
        vec = []
        for j in range(self.params['n']):
            vec.append(curva_transp[j][-1])
        vec.append(1-np.sum(vec))
        zA = np.sum([i*vec[i] for i in range(1, self.params['n'] + 1)])/self.params['n']
        self.assertAlmostEqual(zA, self.params['pa']/(self.params['pa'] + self.params['pb']))
