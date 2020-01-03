% Calculate joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:09
% EndTime: 2020-01-03 11:25:14
% DurationCPUTime: 0.88s
% Computational Cost: add. (1661->204), mult. (2045->292), div. (0->0), fcn. (2125->8), ass. (0->100)
t82 = sin(pkin(8));
t83 = cos(pkin(8));
t85 = sin(qJ(4));
t87 = cos(qJ(4));
t48 = -Icges(6,3) * t83 + (Icges(6,5) * t87 - Icges(6,6) * t85) * t82;
t49 = -Icges(5,3) * t83 + (Icges(5,5) * t87 - Icges(5,6) * t85) * t82;
t129 = -t48 - t49;
t50 = -Icges(6,6) * t83 + (Icges(6,4) * t87 - Icges(6,2) * t85) * t82;
t51 = -Icges(5,6) * t83 + (Icges(5,4) * t87 - Icges(5,2) * t85) * t82;
t128 = (-t50 - t51) * t85;
t81 = qJ(1) + pkin(7);
t76 = sin(t81);
t120 = t76 * t82;
t73 = t87 * pkin(4) + pkin(3);
t121 = t73 * t83;
t113 = t83 * t85;
t77 = cos(t81);
t56 = -t76 * t113 - t77 * t87;
t112 = t83 * t87;
t117 = t77 * t85;
t57 = t76 * t112 - t117;
t127 = t57 * rSges(6,1) + t56 * rSges(6,2) + rSges(6,3) * t120 + t76 * t121;
t52 = -Icges(6,5) * t83 + (Icges(6,1) * t87 - Icges(6,4) * t85) * t82;
t53 = -Icges(5,5) * t83 + (Icges(5,1) * t87 - Icges(5,4) * t85) * t82;
t126 = (t52 + t53) * t82 * t87;
t119 = t76 * t85;
t58 = t77 * t113 - t76 * t87;
t59 = -t77 * t112 - t119;
t125 = -t59 * rSges(6,1) - t58 * rSges(6,2) + pkin(4) * t119;
t124 = t83 ^ 2;
t123 = pkin(3) * t83;
t84 = -qJ(5) - pkin(6);
t122 = pkin(6) + t84;
t118 = t77 * t82;
t116 = t82 * t84;
t111 = t82 * t128 + t129 * t83 + t126;
t94 = Icges(6,3) * t82;
t21 = Icges(6,5) * t57 + Icges(6,6) * t56 + t76 * t94;
t95 = Icges(5,3) * t82;
t23 = Icges(5,5) * t57 + Icges(5,6) * t56 + t76 * t95;
t110 = t23 + t21;
t22 = Icges(6,5) * t59 + Icges(6,6) * t58 - t77 * t94;
t24 = Icges(5,5) * t59 + Icges(5,6) * t58 - t77 * t95;
t109 = -t24 - t22;
t96 = Icges(6,6) * t82;
t26 = Icges(6,4) * t59 + Icges(6,2) * t58 - t77 * t96;
t97 = Icges(5,6) * t82;
t28 = Icges(5,4) * t59 + Icges(5,2) * t58 - t77 * t97;
t108 = t26 + t28;
t25 = Icges(6,4) * t57 + Icges(6,2) * t56 + t76 * t96;
t27 = Icges(5,4) * t57 + Icges(5,2) * t56 + t76 * t97;
t107 = t27 + t25;
t98 = Icges(6,5) * t82;
t30 = Icges(6,1) * t59 + Icges(6,4) * t58 - t77 * t98;
t99 = Icges(5,5) * t82;
t32 = Icges(5,1) * t59 + Icges(5,4) * t58 - t77 * t99;
t106 = t30 + t32;
t29 = Icges(6,1) * t57 + Icges(6,4) * t56 + t76 * t98;
t31 = Icges(5,1) * t57 + Icges(5,4) * t56 + t76 * t99;
t105 = t31 + t29;
t104 = -pkin(4) * t117 + (-t122 * t82 - t123) * t76 + t127;
t102 = pkin(6) * t118 + t77 * t123;
t103 = (t116 - t121) * t77 + t102 - rSges(6,3) * t118 - t125;
t86 = sin(qJ(1));
t78 = t86 * pkin(1);
t101 = t76 * pkin(2) + t78;
t100 = t76 ^ 2 + t77 ^ 2;
t36 = t57 * rSges(5,1) + t56 * rSges(5,2) + rSges(5,3) * t120;
t88 = cos(qJ(1));
t79 = t88 * pkin(1);
t92 = t77 * pkin(2) + t76 * qJ(3) + t79;
t91 = ((-t122 + rSges(6,3)) * t83 + (-rSges(6,1) * t87 + rSges(6,2) * t85 + pkin(3) - t73) * t82) * t82;
t90 = rSges(4,1) * t83 - rSges(4,2) * t82;
t38 = t59 * rSges(5,1) + t58 * rSges(5,2) - rSges(5,3) * t118;
t67 = t88 * rSges(2,1) - t86 * rSges(2,2);
t66 = t86 * rSges(2,1) + t88 * rSges(2,2);
t61 = t77 * rSges(3,1) - t76 * rSges(3,2) + t79;
t60 = t76 * rSges(3,1) + t77 * rSges(3,2) + t78;
t55 = -t83 * rSges(5,3) + (rSges(5,1) * t87 - rSges(5,2) * t85) * t82;
t40 = t76 * rSges(4,3) + t90 * t77 + t92;
t39 = (-rSges(4,3) - qJ(3)) * t77 + t90 * t76 + t101;
t18 = -t55 * t118 + t83 * t38;
t17 = -t55 * t120 - t83 * t36;
t16 = -t38 + t92 + t102;
t15 = -t77 * qJ(3) + (pkin(6) * t82 + t123) * t76 + t36 + t101;
t14 = (t121 + (rSges(6,3) - t84) * t82) * t77 + t92 + t125;
t13 = -t76 * t116 + (-pkin(4) * t85 - qJ(3)) * t77 + t101 + t127;
t12 = -t49 * t118 + t58 * t51 + t59 * t53;
t11 = -t48 * t118 + t58 * t50 + t59 * t52;
t10 = t49 * t120 + t56 * t51 + t57 * t53;
t9 = t48 * t120 + t56 * t50 + t57 * t52;
t8 = (t36 * t77 + t38 * t76) * t82;
t7 = t103 * t83 + t77 * t91;
t6 = -t104 * t83 + t76 * t91;
t5 = -t83 * t24 + (-t28 * t85 + t32 * t87) * t82;
t4 = -t83 * t23 + (-t27 * t85 + t31 * t87) * t82;
t3 = -t83 * t22 + (-t26 * t85 + t30 * t87) * t82;
t2 = -t83 * t21 + (-t25 * t85 + t29 * t87) * t82;
t1 = (t103 * t76 + t104 * t77) * t82;
t19 = [Icges(2,3) + Icges(3,3) + (Icges(4,2) * t83 + t129) * t83 + (Icges(4,1) * t82 + 0.2e1 * Icges(4,4) * t83 + t128) * t82 + m(2) * (t66 ^ 2 + t67 ^ 2) + m(3) * (t60 ^ 2 + t61 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2) + t126; 0; m(3) + m(4) + m(5) + m(6); m(4) * (-t76 * t39 - t77 * t40) + m(5) * (-t76 * t15 - t77 * t16) + m(6) * (-t76 * t13 - t77 * t14); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t100; -t111 * t83 + m(5) * (t17 * t15 + t18 * t16) + m(6) * (t6 * t13 + t7 * t14) + ((-t5 / 0.2e1 - t3 / 0.2e1 - t12 / 0.2e1 - t11 / 0.2e1) * t77 + (t4 / 0.2e1 + t2 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t76) * t82; m(5) * t8 + m(6) * t1; m(5) * (-t17 * t76 - t18 * t77) + m(6) * (-t6 * t76 - t7 * t77); m(5) * (t17 ^ 2 + t18 ^ 2 + t8 ^ 2) + m(6) * (t1 ^ 2 + t6 ^ 2 + t7 ^ 2) + t111 * t124 + (((t106 * t59 + t108 * t58 + t109 * t118) * t118 + (t11 + t12 + t3 + t5) * t83) * t77 + ((t105 * t57 + t107 * t56 + t110 * t120) * t120 + (-t4 - t10 - t2 - t9) * t83 + ((t109 * t76 + t110 * t77) * t82 - t105 * t59 - t107 * t58 - t106 * t57 - t108 * t56) * t118) * t76) * t82; m(6) * (-t13 * t77 + t14 * t76) * t82; -m(6) * t83; 0; m(6) * (-t83 * t1 + (-t6 * t77 + t7 * t76) * t82); m(6) * (t100 * t82 ^ 2 + t124);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
