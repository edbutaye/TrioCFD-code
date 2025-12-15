import unittest

class triocfd_icoco_test(unittest.TestCase):

    def test_small_run(self):
        """ Execute a small run """
        import triocfdicoco as ti
        import medcoupling as mc

        pbT = ti.ProblemTrio()
        pbT.name = "TRUST"
        pbT.setDataFile("Marche.data")
        pbT.initialize()


        def run(pb, tmax):
            stop = False
            t= 0
            while (t < tmax) and (not stop):
                dt, stop = pb.computeTimeStep()
                if stop:
                    return
                pb.initTimeStep(dt)

                # Main time loop:
                ok = pb.solveTimeStep()
                if ok:
                    pb.validateTimeStep()
                    t = pb.presentTime()
                    if pb.isStationary():
                        print("pb stationary -> stop")
                        stop = True
                else:
                    pb.abortTimeStep()
                    dt2, stop2 = pb.computeTimeStep()
                    print("pb abortTimeStep", dt, dt2)
                    assert dt != dt2

        run(pbT,8.e-5)
        pbT.terminate()

if __name__ == "__main__":
    unittest.main()
