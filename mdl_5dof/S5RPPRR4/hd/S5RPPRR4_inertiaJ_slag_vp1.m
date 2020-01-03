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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:30:27
% EndTime: 2020-01-03 11:30:34
% DurationCPUTime: 1.38s
% Computational Cost: add. (3379->260), mult. (3595->385), div. (0->0), fcn. (3787->10), ass. (0->120)
t106 = sin(pkin(8));
t109 = -pkin(6) - qJ(3);
t100 = -pkin(7) + t109;
t124 = t100 - t109;
t142 = t124 * t106;
t110 = sin(qJ(1));
t103 = t110 ^ 2;
t111 = cos(qJ(1));
t104 = t111 ^ 2;
t141 = m(4) / 0.2e1;
t140 = m(5) / 0.2e1;
t139 = m(6) / 0.2e1;
t107 = cos(pkin(9));
t92 = t107 * pkin(3) + pkin(2);
t105 = sin(pkin(9));
t138 = t105 * pkin(3);
t108 = cos(pkin(8));
t102 = pkin(9) + qJ(4);
t95 = qJ(5) + t102;
t90 = sin(t95);
t91 = cos(t95);
t60 = -Icges(6,6) * t108 + (Icges(6,4) * t91 - Icges(6,2) * t90) * t106;
t137 = t90 * t60;
t93 = sin(t102);
t94 = cos(t102);
t64 = -Icges(5,6) * t108 + (Icges(5,4) * t94 - Icges(5,2) * t93) * t106;
t136 = t93 * t64;
t128 = t106 * t111;
t129 = t106 * t110;
t127 = t108 * t110;
t69 = -t111 * t91 - t90 * t127;
t70 = -t111 * t90 + t91 * t127;
t44 = t70 * rSges(6,1) + t69 * rSges(6,2) + rSges(6,3) * t129;
t126 = t108 * t111;
t71 = -t110 * t91 + t90 * t126;
t72 = -t110 * t90 - t91 * t126;
t116 = -t72 * rSges(6,1) - t71 * rSges(6,2);
t45 = -rSges(6,3) * t128 - t116;
t17 = t44 * t128 + t45 * t129;
t135 = t110 * t138 + t92 * t126;
t134 = t111 * pkin(1) + t110 * qJ(2);
t81 = pkin(4) * t94 + t92;
t133 = t108 * t81;
t132 = t108 * t92;
t83 = pkin(4) * t93 + t138;
t131 = t110 * t83;
t61 = -Icges(6,5) * t108 + (Icges(6,1) * t91 - Icges(6,4) * t90) * t106;
t55 = t106 * t91 * t61;
t59 = -Icges(6,3) * t108 + (Icges(6,5) * t91 - Icges(6,6) * t90) * t106;
t23 = -t106 * t137 - t108 * t59 + t55;
t130 = t23 * t108;
t125 = t110 * t111;
t123 = t104 + t103;
t75 = -t111 * t94 - t93 * t127;
t76 = -t111 * t93 + t94 * t127;
t52 = t76 * rSges(5,1) + t75 * rSges(5,2) + rSges(5,3) * t129;
t39 = Icges(6,5) * t72 + Icges(6,6) * t71 - Icges(6,3) * t128;
t41 = Icges(6,4) * t72 + Icges(6,2) * t71 - Icges(6,6) * t128;
t43 = Icges(6,1) * t72 + Icges(6,4) * t71 - Icges(6,5) * t128;
t10 = -t108 * t39 + (-t41 * t90 + t43 * t91) * t106;
t15 = t59 * t129 + t69 * t60 + t70 * t61;
t16 = -t59 * t128 + t71 * t60 + t72 * t61;
t38 = Icges(6,5) * t70 + Icges(6,6) * t69 + Icges(6,3) * t129;
t40 = Icges(6,4) * t70 + Icges(6,2) * t69 + Icges(6,6) * t129;
t42 = Icges(6,1) * t70 + Icges(6,4) * t69 + Icges(6,5) * t129;
t9 = -t108 * t38 + (-t40 * t90 + t42 * t91) * t106;
t120 = (t15 + t9) * t129 / 0.2e1 - (t10 + t16) * t128 / 0.2e1;
t62 = -t108 * rSges(6,3) + (rSges(6,1) * t91 - rSges(6,2) * t90) * t106;
t119 = t106 * (-t124 * t108 - (t81 - t92) * t106 - t62);
t118 = t141 + t140 + t139;
t77 = -t110 * t94 + t93 * t126;
t78 = -t110 * t93 - t94 * t126;
t117 = -t78 * rSges(5,1) - t77 * rSges(5,2);
t115 = rSges(3,1) * t108 - rSges(3,2) * t106;
t114 = t105 * rSges(4,1) + t107 * rSges(4,2);
t1 = (-t15 * t108 + (t129 * t38 + t69 * t40 + t70 * t42) * t129 - (t129 * t39 + t69 * t41 + t70 * t43) * t128) * t129;
t2 = -t16 * t108 + (-t128 * t38 + t71 * t40 + t72 * t42) * t129 - (-t128 * t39 + t71 * t41 + t72 * t43) * t128;
t3 = -t130 + (-t10 * t111 + t110 * t9) * t106;
t113 = -t108 * t3 - t2 * t128 + t1;
t112 = (rSges(4,3) + qJ(3)) * t106 + (rSges(4,1) * t107 - rSges(4,2) * t105 + pkin(2)) * t108;
t98 = t110 * pkin(1);
t85 = t111 * rSges(2,1) - t110 * rSges(2,2);
t84 = t110 * rSges(2,1) + t111 * rSges(2,2);
t79 = t81 * t127;
t66 = -t108 * rSges(5,3) + (rSges(5,1) * t94 - rSges(5,2) * t93) * t106;
t65 = -Icges(5,5) * t108 + (Icges(5,1) * t94 - Icges(5,4) * t93) * t106;
t63 = -Icges(5,3) * t108 + (Icges(5,5) * t94 - Icges(5,6) * t93) * t106;
t58 = t110 * rSges(3,3) + t115 * t111 + t134;
t57 = t98 + (-rSges(3,3) - qJ(2)) * t111 + t115 * t110;
t56 = t106 * t94 * t65;
t53 = -rSges(5,3) * t128 - t117;
t51 = Icges(5,1) * t78 + Icges(5,4) * t77 - Icges(5,5) * t128;
t50 = Icges(5,1) * t76 + Icges(5,4) * t75 + Icges(5,5) * t129;
t49 = Icges(5,4) * t78 + Icges(5,2) * t77 - Icges(5,6) * t128;
t48 = Icges(5,4) * t76 + Icges(5,2) * t75 + Icges(5,6) * t129;
t47 = Icges(5,5) * t78 + Icges(5,6) * t77 - Icges(5,3) * t128;
t46 = Icges(5,5) * t76 + Icges(5,6) * t75 + Icges(5,3) * t129;
t37 = t108 * t45;
t34 = t114 * t110 + t112 * t111 + t134;
t33 = t98 + (-qJ(2) - t114) * t111 + t112 * t110;
t32 = -t131 + (-t133 + t142) * t111 + t135;
t31 = t79 + (-t83 + t138) * t111 + (-t132 - t142) * t110;
t30 = t108 * t53 - t66 * t128;
t29 = -t108 * t52 - t66 * t129;
t28 = (rSges(5,3) - t109) * t128 + t117 + t134 + t135;
t27 = t98 + (-qJ(2) - t138) * t111 + (-t106 * t109 + t132) * t110 + t52;
t26 = -t62 * t128 + t37;
t25 = -t108 * t44 - t62 * t129;
t24 = -t106 * t136 - t108 * t63 + t56;
t22 = t131 + (t133 + (rSges(6,3) - t100) * t106) * t111 + t116 + t134;
t21 = -t100 * t129 + t79 + t98 + (-qJ(2) - t83) * t111 + t44;
t20 = (t110 * t53 + t111 * t52) * t106;
t19 = -t63 * t128 + t77 * t64 + t78 * t65;
t18 = t63 * t129 + t75 * t64 + t76 * t65;
t12 = -t108 * t47 + (-t49 * t93 + t51 * t94) * t106;
t11 = -t108 * t46 + (-t48 * t93 + t50 * t94) * t106;
t6 = t108 * t32 + t111 * t119 + t37;
t5 = (-t31 - t44) * t108 + t110 * t119;
t4 = (t110 * t32 + t111 * t31) * t106 + t17;
t7 = [Icges(2,3) + t55 + t56 + (-t63 - t59 + (Icges(3,2) + Icges(4,3)) * t108) * t108 + m(2) * (t84 ^ 2 + t85 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + (-t136 - t137 + (Icges(4,1) * t107 ^ 2 + Icges(3,1) + (-0.2e1 * Icges(4,4) * t107 + Icges(4,2) * t105) * t105) * t106 + 0.2e1 * (-Icges(4,5) * t107 + Icges(4,6) * t105 + Icges(3,4)) * t108) * t106; m(3) * (-t110 * t57 - t111 * t58) + m(4) * (-t110 * t33 - t111 * t34) + m(5) * (-t110 * t27 - t111 * t28) + m(6) * (-t110 * t21 - t111 * t22); 0.2e1 * (m(3) / 0.2e1 + t118) * t123; 0.2e1 * ((t110 * t34 - t111 * t33) * t141 + (t110 * t28 - t111 * t27) * t140 + (t110 * t22 - t111 * t21) * t139) * t106; 0; 0.2e1 * t118 * (t106 ^ 2 * t123 + t108 ^ 2); (-t24 - t23) * t108 + m(5) * (t29 * t27 + t30 * t28) + m(6) * (t5 * t21 + t6 * t22) + ((-t12 / 0.2e1 - t19 / 0.2e1) * t111 + (t11 / 0.2e1 + t18 / 0.2e1) * t110) * t106 + t120; m(5) * (-t29 * t110 - t30 * t111) + m(6) * (-t5 * t110 - t6 * t111); m(5) * (-t20 * t108 + (t110 * t30 - t111 * t29) * t106) + m(6) * (-t4 * t108 + (t110 * t6 - t111 * t5) * t106); t1 + (t24 * t108 - t3) * t108 + (-t111 * t2 + (t110 * ((t75 * t48 + t76 * t50) * t110 - (t75 * t49 + t76 * t51) * t111) - t111 * ((t77 * t48 + t78 * t50) * t110 - (t77 * t49 + t78 * t51) * t111) + (t110 * (t103 * t46 - t125 * t47) - t111 * (t104 * t47 - t125 * t46)) * t106) * t106 + ((t12 + t19) * t111 + (-t11 - t18) * t110) * t108) * t106 + m(5) * (t20 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2); -t130 + m(6) * (t25 * t21 + t26 * t22) + t120; m(6) * (-t25 * t110 - t26 * t111); m(6) * (-t17 * t108 + (t110 * t26 - t111 * t25) * t106); m(6) * (t17 * t4 + t25 * t5 + t26 * t6) + t113; m(6) * (t17 ^ 2 + t25 ^ 2 + t26 ^ 2) + t113;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
