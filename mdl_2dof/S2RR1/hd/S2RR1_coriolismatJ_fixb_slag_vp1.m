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

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_coriolismatJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:06
% EndTime: 2020-01-03 11:19:08
% DurationCPUTime: 0.56s
% Computational Cost: add. (688->82), mult. (1660->134), div. (0->0), fcn. (1556->4), ass. (0->61)
t51 = sin(qJ(1));
t52 = cos(qJ(2));
t53 = cos(qJ(1));
t74 = t52 * t53;
t50 = sin(qJ(2));
t77 = t50 * t53;
t62 = -rSges(3,1) * t74 + rSges(3,2) * t77;
t79 = rSges(3,3) + pkin(1);
t18 = t79 * t51 + t62;
t76 = t51 * t52;
t78 = t50 * t51;
t63 = -rSges(3,1) * t76 + rSges(3,2) * t78;
t19 = -t79 * t53 + t63;
t60 = rSges(3,1) * t50 + rSges(3,2) * t52;
t34 = t60 * t51;
t35 = t60 * t53;
t67 = Icges(3,2) * t52;
t71 = Icges(3,4) * t50;
t37 = -t67 - t71;
t72 = Icges(3,1) * t52;
t58 = -t71 + t72;
t88 = -(t58 / 0.2e1 + t37 / 0.2e1) * t50 - m(3) * (t18 * t35 + t19 * t34);
t70 = Icges(3,4) * t52;
t56 = -Icges(3,2) * t50 + t70;
t86 = t51 ^ 2;
t85 = t53 ^ 2;
t83 = -t51 / 0.2e1;
t80 = t53 / 0.2e1;
t75 = t51 * t53;
t69 = Icges(3,5) * t52;
t55 = -Icges(3,6) * t50 + t69;
t22 = -Icges(3,3) * t51 + t55 * t53;
t26 = -Icges(3,5) * t51 + t58 * t53;
t73 = t53 * t22 + t26 * t76;
t66 = Icges(3,6) * t51;
t17 = t26 * t74;
t65 = -t51 * t22 + t17;
t21 = Icges(3,3) * t53 - t50 * t66 + t51 * t69;
t24 = t56 * t53 - t66;
t64 = -t50 * t24 - t21;
t23 = Icges(3,6) * t53 + t56 * t51;
t47 = t51 * t71;
t25 = Icges(3,5) * t53 + t51 * t72 - t47;
t61 = -t23 * t77 + t25 * t74;
t59 = -t50 * t23 + t52 * t25;
t57 = Icges(3,1) * t50 + t70;
t54 = -Icges(3,5) * t50 - Icges(3,6) * t52;
t43 = -rSges(3,1) * t52 + rSges(3,2) * t50;
t29 = t53 * t54;
t28 = t54 * t51;
t11 = -t24 * t77 + t65;
t10 = -t51 * t21 + t61;
t9 = -t24 * t78 + t73;
t8 = t53 * t21 + t59 * t51;
t6 = (t57 / 0.2e1 + t56 / 0.2e1) * t52 - t88;
t5 = t10 * t53 - t11 * t51;
t4 = -t9 * t51 + t53 * t8;
t3 = (t61 - t9 + t73) * t53 + (-t17 - t8 + (t22 + t59) * t51) * t51;
t2 = (t64 * t51 - t10 + t73) * t51 + (t64 * t53 - t11 + t65) * t53;
t1 = (t3 / 0.2e1 - t5 / 0.2e1) * t53 + (-t4 / 0.2e1 - t2 / 0.2e1) * t51;
t7 = [t6 * qJD(2), t6 * qJD(1) + (m(3) * ((-t18 * t51 + t19 * t53) * t43 - (t34 * t53 - t35 * t51) * t60) + ((-t37 * t53 - t26) * t52 + (t57 * t53 + t24) * t50) * t83 - t53 * t3 / 0.2e1 - (t85 / 0.2e1 + t86 / 0.2e1) * t55 + (t4 + t2) * t51 / 0.2e1 + ((t51 * t67 - t25 + t47) * t52 + (t57 * t51 + t23) * t50 + t5) * t80) * qJD(2); ((-t56 - t57) * t52 / 0.2e1 + t88) * qJD(1) + t1 * qJD(2), t1 * qJD(1) + (m(3) * ((t51 * t63 + t53 * t62) * (t51 * t34 + t35 * t53) - (t85 + t86) * t43 * t60) + (t85 * t28 - t29 * t75) * t80 + (-t28 * t75 + t86 * t29) * t83) * qJD(2);];
Cq = t7;
