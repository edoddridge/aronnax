Verification
************************

Aronnax includes a number of diagnostic test to verify that the numerical core is satisfactorily solving the equations of motion.


Conservation of volume
========================

The use of sponge regions affects volume conservation in both global and layerwise sense, depending on the configuration. The following table shows when global or layerwise volume conservation can be expected.


+------------------+-----------------+----------------------+
| Physics          | Sponge regions? | Volume conservation  |
+------------------+-----------------+----------+-----------+
|                  |                 |  Global  | Layerwise |
+==================+=================+==========+===========+
| n-layer          |       Yes       |  Yes     | No        |
+------------------+-----------------+----------+-----------+
| n-layer          |       No        |  Yes     | Yes       |
+------------------+-----------------+----------+-----------+
| n + 1/2 layer    |       Yes       |  No      | No        |
+------------------+-----------------+----------+-----------+
| n + 1/2 layer    |       No        |  Yes     | Yes       |
+------------------+-----------------+----------+-----------+

TODO:
- Add graphs of layerwise/global volume in simulations with(out) sponges


Momentum
==========================
To do

