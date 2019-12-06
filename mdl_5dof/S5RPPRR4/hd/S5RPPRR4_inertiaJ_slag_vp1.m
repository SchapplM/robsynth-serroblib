% Calculate joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:05
% EndTime: 2019-12-05 17:44:11
% DurationCPUTime: 1.30s
% Computational Cost: add. (3379->261), mult. (3595->387), div. (0->0), fcn. (3787->10), ass. (0->121)
t104 = cos(pkin(8));
t103 = cos(pkin(9));
t90 = t103 * pkin(3) + pkin(2);
t98 = pkin(9) + qJ(4);
t92 = cos(t98);
t81 = pkin(4) * t92 + t90;
t128 = t81 - t90;
t139 = t128 * t104;
t106 = sin(qJ(1));
t99 = t106 ^ 2;
t107 = cos(qJ(1));
t100 = t107 ^ 2;
t138 = m(4) / 0.2e1;
t137 = m(5) / 0.2e1;
t136 = m(6) / 0.2e1;
t101 = sin(pkin(9));
t135 = t101 * pkin(3);
t102 = sin(pkin(8));
t93 = qJ(5) + t98;
t88 = sin(t93);
t89 = cos(t93);
t60 = -Icges(6,6) * t104 + (Icges(6,4) * t89 - Icges(6,2) * t88) * t102;
t134 = t88 * t60;
t91 = sin(t98);
t64 = -Icges(5,6) * t104 + (Icges(5,4) * t92 - Icges(5,2) * t91) * t102;
t133 = t91 * t64;
t105 = -pkin(6) - qJ(3);
t120 = t107 * t135;
t125 = t102 * t106;
t82 = pkin(4) * t91 + t135;
t96 = -pkin(7) + t105;
t129 = t107 * t82 + t96 * t125;
t123 = t104 * t106;
t69 = t107 * t89 + t88 * t123;
t70 = t107 * t88 - t89 * t123;
t131 = t70 * rSges(6,1) + t69 * rSges(6,2);
t42 = -rSges(6,3) * t125 + t131;
t132 = t120 - (-t102 * t105 - t139) * t106 - t129 - t42;
t124 = t102 * t107;
t122 = t104 * t107;
t71 = t106 * t89 - t88 * t122;
t72 = t106 * t88 + t89 * t122;
t112 = -t72 * rSges(6,1) - t71 * rSges(6,2);
t43 = rSges(6,3) * t124 - t112;
t62 = -t104 * rSges(6,3) + (rSges(6,1) * t89 - rSges(6,2) * t88) * t102;
t26 = t104 * t43 + t62 * t124;
t75 = t107 * t92 + t91 * t123;
t76 = t107 * t91 - t92 * t123;
t130 = t76 * rSges(5,1) + t75 * rSges(5,2);
t61 = -Icges(6,5) * t104 + (Icges(6,1) * t89 - Icges(6,4) * t88) * t102;
t53 = t102 * t89 * t61;
t59 = -Icges(6,3) * t104 + (Icges(6,5) * t89 - Icges(6,6) * t88) * t102;
t23 = -t102 * t134 - t104 * t59 + t53;
t127 = t23 * t104;
t126 = t100 + t99;
t121 = t106 * t107;
t37 = Icges(6,5) * t72 + Icges(6,6) * t71 + Icges(6,3) * t124;
t39 = Icges(6,4) * t72 + Icges(6,2) * t71 + Icges(6,6) * t124;
t41 = Icges(6,1) * t72 + Icges(6,4) * t71 + Icges(6,5) * t124;
t10 = -t104 * t37 + (-t39 * t88 + t41 * t89) * t102;
t15 = -t59 * t125 + t69 * t60 + t70 * t61;
t16 = t59 * t124 + t71 * t60 + t72 * t61;
t36 = Icges(6,5) * t70 + Icges(6,6) * t69 - Icges(6,3) * t125;
t38 = Icges(6,4) * t70 + Icges(6,2) * t69 - Icges(6,6) * t125;
t40 = Icges(6,1) * t70 + Icges(6,4) * t69 - Icges(6,5) * t125;
t9 = -t104 * t36 + (-t38 * t88 + t40 * t89) * t102;
t117 = -(t15 + t9) * t125 / 0.2e1 + (t10 + t16) * t124 / 0.2e1;
t116 = -t104 * t81 - pkin(1);
t115 = -t104 * t90 - pkin(1);
t114 = t138 + t137 + t136;
t77 = t106 * t92 - t91 * t122;
t78 = t106 * t91 + t92 * t122;
t113 = -t78 * rSges(5,1) - t77 * rSges(5,2);
t111 = t101 * rSges(4,1) + t103 * rSges(4,2);
t110 = -rSges(3,1) * t104 + rSges(3,2) * t102 - pkin(1);
t1 = (-t16 * t104 - (t36 * t124 + t71 * t38 + t72 * t40) * t125 + (t37 * t124 + t71 * t39 + t72 * t41) * t124) * t124;
t2 = -t15 * t104 - (-t36 * t125 + t69 * t38 + t70 * t40) * t125 + (-t37 * t125 + t69 * t39 + t70 * t41) * t124;
t3 = -t127 + (t10 * t107 - t106 * t9) * t102;
t109 = -t104 * t3 - t2 * t125 + t1;
t108 = -pkin(1) + (-rSges(4,3) - qJ(3)) * t102 + (-rSges(4,1) * t103 + rSges(4,2) * t101 - pkin(2)) * t104;
t95 = t107 * qJ(2);
t86 = t105 * t124;
t85 = -t107 * rSges(2,1) + t106 * rSges(2,2);
t84 = -t106 * rSges(2,1) - t107 * rSges(2,2);
t66 = -t104 * rSges(5,3) + (rSges(5,1) * t92 - rSges(5,2) * t91) * t102;
t65 = -Icges(5,5) * t104 + (Icges(5,1) * t92 - Icges(5,4) * t91) * t102;
t63 = -Icges(5,3) * t104 + (Icges(5,5) * t92 - Icges(5,6) * t91) * t102;
t58 = (-rSges(3,3) - qJ(2)) * t106 + t110 * t107;
t57 = t107 * rSges(3,3) + t110 * t106 + t95;
t55 = t62 * t125;
t54 = t102 * t92 * t65;
t52 = (-t105 + t96) * t104 + t128 * t102;
t51 = rSges(5,3) * t124 - t113;
t50 = -rSges(5,3) * t125 + t130;
t49 = Icges(5,1) * t78 + Icges(5,4) * t77 + Icges(5,5) * t124;
t48 = Icges(5,1) * t76 + Icges(5,4) * t75 - Icges(5,5) * t125;
t47 = Icges(5,4) * t78 + Icges(5,2) * t77 + Icges(5,6) * t124;
t46 = Icges(5,4) * t76 + Icges(5,2) * t75 - Icges(5,6) * t125;
t45 = Icges(5,5) * t78 + Icges(5,6) * t77 + Icges(5,3) * t124;
t44 = Icges(5,5) * t76 + Icges(5,6) * t75 - Icges(5,3) * t125;
t34 = (-qJ(2) - t111) * t106 + t108 * t107;
t33 = t108 * t106 + t111 * t107 + t95;
t32 = t86 + (t82 - t135) * t106 + (-t102 * t96 + t139) * t107;
t30 = t104 * t51 + t66 * t124;
t29 = -t104 * t50 + t66 * t125;
t28 = t86 + (-qJ(2) - t135) * t106 + (-rSges(5,3) * t102 + t115) * t107 + t113;
t27 = t120 + t95 + ((-rSges(5,3) + t105) * t102 + t115) * t106 + t130;
t25 = -t104 * t42 + t55;
t24 = -t102 * t133 - t104 * t63 + t54;
t22 = (-qJ(2) - t82) * t106 + ((-rSges(6,3) + t96) * t102 + t116) * t107 + t112;
t21 = t95 + (-rSges(6,3) * t102 + t116) * t106 + t129 + t131;
t20 = (-t106 * t51 - t107 * t50) * t102;
t19 = t63 * t124 + t77 * t64 + t78 * t65;
t18 = -t63 * t125 + t75 * t64 + t76 * t65;
t17 = (-t106 * t43 - t107 * t42) * t102;
t12 = -t104 * t45 + (-t47 * t91 + t49 * t92) * t102;
t11 = -t104 * t44 + (-t46 * t91 + t48 * t92) * t102;
t6 = t104 * t32 + t52 * t124 + t26;
t5 = t132 * t104 + t52 * t125 + t55;
t4 = (t132 * t107 + (-t32 - t43) * t106) * t102;
t7 = [Icges(2,3) + t53 + t54 + (-t63 - t59 + (Icges(3,2) + Icges(4,3)) * t104) * t104 + m(2) * (t84 ^ 2 + t85 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + (-t133 - t134 + (Icges(4,1) * t103 ^ 2 + Icges(3,1) + (-0.2e1 * Icges(4,4) * t103 + Icges(4,2) * t101) * t101) * t102 + 0.2e1 * (-Icges(4,5) * t103 + Icges(4,6) * t101 + Icges(3,4)) * t104) * t102; m(3) * (t106 * t57 + t107 * t58) + m(4) * (t106 * t33 + t107 * t34) + m(5) * (t106 * t27 + t107 * t28) + m(6) * (t106 * t21 + t107 * t22); 0.2e1 * (m(3) / 0.2e1 + t114) * t126; 0.2e1 * ((-t106 * t34 + t107 * t33) * t138 + (-t106 * t28 + t107 * t27) * t137 + (-t106 * t22 + t107 * t21) * t136) * t102; 0; 0.2e1 * t114 * (t126 * t102 ^ 2 + t104 ^ 2); (-t24 - t23) * t104 + m(5) * (t29 * t27 + t30 * t28) + m(6) * (t5 * t21 + t6 * t22) + ((t12 / 0.2e1 + t19 / 0.2e1) * t107 + (-t11 / 0.2e1 - t18 / 0.2e1) * t106) * t102 + t117; m(5) * (t29 * t106 + t30 * t107) + m(6) * (t5 * t106 + t6 * t107); m(5) * (-t20 * t104 + (-t106 * t30 + t107 * t29) * t102) + m(6) * (-t4 * t104 + (-t106 * t6 + t107 * t5) * t102); t1 + (t24 * t104 - t3) * t104 + (-t106 * t2 + (-t106 * (-(t75 * t46 + t76 * t48) * t106 + (t75 * t47 + t76 * t49) * t107) + t107 * (-(t77 * t46 + t78 * t48) * t106 + (t77 * t47 + t78 * t49) * t107) + (-t106 * (-t45 * t121 + t99 * t44) + t107 * (t100 * t45 - t44 * t121)) * t102) * t102 + ((-t12 - t19) * t107 + (t11 + t18) * t106) * t104) * t102 + m(5) * (t20 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2); -t127 + m(6) * (t25 * t21 + t26 * t22) + t117; m(6) * (t25 * t106 + t26 * t107); m(6) * (-t17 * t104 + (-t106 * t26 + t107 * t25) * t102); m(6) * (t17 * t4 + t25 * t5 + t26 * t6) + t109; m(6) * (t17 ^ 2 + t25 ^ 2 + t26 ^ 2) + t109;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
