% Calculate matrix of centrifugal and coriolis load on the joints for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [2x2]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S2RR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert( isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:44:20
% EndTime: 2018-11-16 16:44:20
% DurationCPUTime: 0.49s
% Computational Cost: add. (688->82), mult. (1660->134), div. (0->0), fcn. (1556->4), ass. (0->60)
t53 = sin(qJ(1));
t54 = cos(qJ(2));
t55 = cos(qJ(1));
t75 = t54 * t55;
t52 = sin(qJ(2));
t79 = t52 * t55;
t61 = -rSges(3,1) * t75 + rSges(3,2) * t79;
t82 = rSges(3,3) + pkin(1);
t18 = t82 * t53 + t61;
t80 = t52 * t53;
t81 = rSges(3,1) * t54;
t62 = -rSges(3,2) * t80 + t53 * t81;
t19 = t82 * t55 + t62;
t60 = rSges(3,1) * t52 + rSges(3,2) * t54;
t35 = t60 * t53;
t36 = t60 * t55;
t66 = Icges(3,2) * t54;
t70 = Icges(3,4) * t52;
t71 = Icges(3,1) * t54;
t91 = (t70 - t71 / 0.2e1 + t66 / 0.2e1) * t52 - m(3) * (t18 * t36 - t19 * t35);
t69 = Icges(3,4) * t54;
t39 = Icges(3,2) * t52 - t69;
t89 = t53 ^ 2;
t88 = t55 ^ 2;
t85 = t53 / 0.2e1;
t83 = -t55 / 0.2e1;
t78 = t53 * t55;
t48 = t53 * t70;
t26 = Icges(3,5) * t55 + t53 * t71 - t48;
t77 = t54 * t26;
t49 = t55 * t70;
t27 = Icges(3,5) * t53 - t55 * t71 + t49;
t76 = t54 * t27;
t65 = Icges(3,6) * t53;
t68 = Icges(3,5) * t54;
t22 = Icges(3,3) * t55 - t52 * t65 + t53 * t68;
t24 = Icges(3,6) * t55 - t39 * t53;
t73 = t53 * t22 + t24 * t79;
t57 = Icges(3,6) * t52 - t68;
t23 = Icges(3,3) * t53 + t57 * t55;
t25 = t39 * t55 + t65;
t72 = t53 * t23 + t25 * t79;
t15 = t25 * t80;
t64 = t55 * t23 - t15;
t63 = t22 + t76;
t58 = Icges(3,1) * t52 + t69;
t56 = Icges(3,5) * t52 + Icges(3,6) * t54;
t44 = rSges(3,2) * t52 - t81;
t30 = t56 * t55;
t29 = t56 * t53;
t11 = -t27 * t75 + t72;
t10 = t26 * t75 - t73;
t9 = -t53 * t76 - t64;
t6 = (t58 / 0.2e1 - t39 / 0.2e1) * t54 - t91;
t5 = -t10 * t55 + t11 * t53;
t4 = t9 * t53 - t55 * (-(t52 * t24 - t77) * t53 + t55 * t22);
t3 = (-t15 + t9 + (t23 - t77) * t55 + t73) * t55 + t72 * t53;
t2 = (t63 * t55 + t11 - t72) * t55 + (t63 * t53 + t10 + t64) * t53;
t1 = (-t3 / 0.2e1 + t5 / 0.2e1) * t55 + (t4 / 0.2e1 + t2 / 0.2e1) * t53;
t7 = [t6 * qJD(2), t6 * qJD(1) + (m(3) * ((-t18 * t53 - t19 * t55) * t44 - (t35 * t55 - t36 * t53) * t60) + ((-t55 * t66 - t27 - t49) * t54 + (-t58 * t55 + t25) * t52) * t85 + t55 * t3 / 0.2e1 + (t88 / 0.2e1 + t89 / 0.2e1) * t57 - (t4 + t2) * t53 / 0.2e1 + ((-t53 * t66 + t26 - t48) * t54 + (-t58 * t53 - t24) * t52 + t5) * t83) * qJD(2); ((t39 - t58) * t54 / 0.2e1 + t91) * qJD(1) + t1 * qJD(2), t1 * qJD(1) + (m(3) * ((-t53 * t62 + t55 * t61) * (t53 * t35 + t36 * t55) - (t88 + t89) * t44 * t60) + (t88 * t29 - t30 * t78) * t83 + (-t29 * t78 + t89 * t30) * t85) * qJD(2);];
Cq  = t7;
