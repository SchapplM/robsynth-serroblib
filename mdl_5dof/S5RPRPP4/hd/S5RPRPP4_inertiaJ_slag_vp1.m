% Calculate joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:21
% DurationCPUTime: 1.00s
% Computational Cost: add. (828->132), mult. (1176->180), div. (0->0), fcn. (1000->6), ass. (0->68)
t74 = qJ(3) + pkin(7);
t67 = cos(t74);
t143 = t67 ^ 2;
t79 = sin(qJ(1));
t75 = t79 ^ 2;
t81 = cos(qJ(1));
t76 = t81 ^ 2;
t58 = t75 + t76;
t141 = Icges(6,4) + Icges(5,5);
t140 = Icges(5,6) - Icges(6,6);
t109 = rSges(6,3) + qJ(5);
t104 = t109 * t67;
t131 = rSges(6,1) + pkin(4);
t66 = sin(t74);
t137 = (t131 * t66 - t104) * t81;
t78 = sin(qJ(3));
t80 = cos(qJ(3));
t136 = Icges(4,5) * t78 + Icges(4,6) * t80 + t140 * t67 + t141 * t66;
t135 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t130 = (rSges(4,1) * t78 + rSges(4,2) * t80) * t81;
t129 = t135 * t79 - t136 * t81;
t128 = t135 * t81 + t136 * t79;
t125 = m(6) * t67;
t124 = pkin(3) * t80;
t123 = t66 * t79;
t122 = t67 * t79;
t121 = t78 * t79;
t120 = t79 * t80;
t70 = t81 * rSges(5,3);
t119 = t109 * t66 + t131 * t67;
t118 = (m(5) + m(6)) * t58;
t77 = -qJ(4) - pkin(6);
t117 = t81 * t78 * pkin(3) + t79 * t77;
t116 = t81 * pkin(1) + t79 * qJ(2);
t108 = -t81 * rSges(6,2) - t123 * t131;
t107 = -rSges(5,1) * t123 - rSges(5,2) * t122 - t70;
t106 = rSges(4,1) * t121 + rSges(4,2) * t120 + t81 * rSges(4,3);
t69 = t81 * qJ(2);
t105 = t69 + t117;
t3 = (-rSges(6,2) - pkin(1)) * t79 + t137 + t105;
t61 = pkin(3) * t121;
t83 = -t81 * t77 + t116 + t61;
t4 = -t79 * t104 - t108 + t83;
t102 = t79 * t3 - t81 * t4;
t62 = pkin(3) * t120;
t8 = t119 * t79 + t62;
t9 = (-t119 - t124) * t81;
t101 = t8 * t79 - t9 * t81;
t99 = rSges(5,1) * t66 + rSges(5,2) * t67;
t10 = t69 + t130 + (-rSges(4,3) - pkin(1) - pkin(6)) * t79;
t11 = t81 * pkin(6) + t106 + t116;
t82 = m(4) * (t79 * t10 - t81 * t11);
t53 = t81 * rSges(2,1) - t79 * rSges(2,2);
t52 = t80 * rSges(4,1) - t78 * rSges(4,2);
t51 = -t79 * rSges(2,1) - t81 * rSges(2,2);
t44 = t67 * rSges(5,1) - t66 * rSges(5,2);
t35 = t61 + (-pkin(6) - t77) * t81;
t34 = -t81 * rSges(3,2) + t79 * rSges(3,3) + t116;
t33 = t81 * rSges(3,3) + t69 + (rSges(3,2) - pkin(1)) * t79;
t26 = t81 * (-t79 * pkin(6) - t117);
t13 = (-t44 - t124) * t81;
t12 = t79 * t44 + t62;
t7 = t83 - t107;
t6 = t99 * t81 + (-rSges(5,3) - pkin(1)) * t79 + t105;
t5 = -t79 * t106 + (t79 * rSges(4,3) - t130) * t81;
t2 = t26 - t99 * t76 + (t107 - t35 + t70) * t79;
t1 = t26 + (t109 * t122 + t108 - t35) * t79 + (t79 * rSges(6,2) - t137) * t81;
t14 = [Icges(4,1) * t80 ^ 2 + Icges(3,1) + Icges(2,3) + m(6) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(3) * (t33 ^ 2 + t34 ^ 2) + m(2) * (t51 ^ 2 + t53 ^ 2) + (Icges(5,1) + Icges(6,1)) * t143 + (-0.2e1 * Icges(4,4) * t80 + Icges(4,2) * t78) * t78 + (0.2e1 * (Icges(6,5) - Icges(5,4)) * t67 + (Icges(5,2) + Icges(6,3)) * t66) * t66; m(6) * t102 + m(5) * (t79 * t6 - t81 * t7) + t82 + m(3) * (t79 * t33 - t81 * t34); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t58 + t118; m(6) * (t8 * t3 + t9 * t4) + m(5) * (t12 * t6 + t13 * t7) + t52 * t82 + t58 * (Icges(4,5) * t80 - Icges(4,6) * t78 - t140 * t66 + t141 * t67); m(5) * (t12 * t79 - t13 * t81) + m(6) * t101 + m(4) * t58 * t52; m(6) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t2 ^ 2) + m(4) * (t58 * t52 ^ 2 + t5 ^ 2) + t129 * t79 * t75 + (t128 * t76 + (t128 * t79 + t129 * t81) * t79) * t81; m(6) * (t81 * t3 + t79 * t4) + m(5) * (t81 * t6 + t79 * t7); 0; m(6) * (t79 * t9 + t81 * t8) + m(5) * (t81 * t12 + t79 * t13); t118; -t102 * t125; -t58 * t125; m(6) * (t66 * t1 - t101 * t67); 0; m(6) * (t58 * t143 + t66 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
