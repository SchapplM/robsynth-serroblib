% Calculate joint inertia matrix for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR7_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR7_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:04
% EndTime: 2019-12-31 16:36:08
% DurationCPUTime: 1.37s
% Computational Cost: add. (5870->286), mult. (15146->446), div. (0->0), fcn. (19328->10), ass. (0->135)
t119 = cos(pkin(4));
t122 = sin(qJ(2));
t117 = sin(pkin(4));
t121 = sin(qJ(3));
t132 = t117 * t121;
t138 = cos(qJ(3));
t110 = -t119 * t138 + t122 * t132;
t116 = sin(pkin(8));
t118 = cos(pkin(8));
t124 = cos(qJ(2));
t128 = t124 * t119;
t106 = t116 * t122 - t118 * t128;
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t129 = t122 * t119;
t107 = t116 * t124 + t118 * t129;
t130 = t118 * t117;
t97 = t107 * t138 - t121 * t130;
t74 = t106 * t123 - t97 * t120;
t75 = t106 * t120 + t97 * t123;
t127 = t117 * t138;
t96 = t107 * t121 + t118 * t127;
t49 = Icges(5,5) * t75 + Icges(5,6) * t74 + Icges(5,3) * t96;
t51 = Icges(5,4) * t75 + Icges(5,2) * t74 + Icges(5,6) * t96;
t53 = Icges(5,1) * t75 + Icges(5,4) * t74 + Icges(5,5) * t96;
t108 = t116 * t128 + t118 * t122;
t109 = -t116 * t129 + t118 * t124;
t99 = t109 * t138 + t116 * t132;
t76 = t108 * t123 - t99 * t120;
t77 = t108 * t120 + t99 * t123;
t98 = t109 * t121 - t116 * t127;
t18 = t98 * t49 + t76 * t51 + t77 * t53;
t50 = Icges(5,5) * t77 + Icges(5,6) * t76 + Icges(5,3) * t98;
t52 = Icges(5,4) * t77 + Icges(5,2) * t76 + Icges(5,6) * t98;
t54 = Icges(5,1) * t77 + Icges(5,4) * t76 + Icges(5,5) * t98;
t19 = t98 * t50 + t76 * t52 + t77 * t54;
t111 = t119 * t121 + t122 * t127;
t131 = t117 * t124;
t100 = -t111 * t120 - t123 * t131;
t101 = t111 * t123 - t120 * t131;
t66 = Icges(5,5) * t101 + Icges(5,6) * t100 + Icges(5,3) * t110;
t67 = Icges(5,4) * t101 + Icges(5,2) * t100 + Icges(5,6) * t110;
t68 = Icges(5,1) * t101 + Icges(5,4) * t100 + Icges(5,5) * t110;
t27 = t98 * t66 + t76 * t67 + t77 * t68;
t2 = t27 * t110 + t18 * t96 + t19 * t98;
t142 = t2 / 0.2e1;
t141 = t96 / 0.2e1;
t140 = t98 / 0.2e1;
t139 = t110 / 0.2e1;
t55 = t75 * rSges(5,1) + t74 * rSges(5,2) + t96 * rSges(5,3);
t137 = t97 * pkin(3) + t96 * pkin(7) + t55;
t56 = t77 * rSges(5,1) + t76 * rSges(5,2) + t98 * rSges(5,3);
t136 = t99 * pkin(3) + t98 * pkin(7) + t56;
t69 = t101 * rSges(5,1) + t100 * rSges(5,2) + t110 * rSges(5,3);
t135 = t111 * pkin(3) + t110 * pkin(7) + t69;
t133 = t116 * t117;
t93 = t107 * pkin(2) + t106 * pkin(6);
t94 = t109 * pkin(2) + t108 * pkin(6);
t134 = t94 * t130 + t93 * t133;
t112 = (pkin(2) * t122 - pkin(6) * t124) * t117;
t89 = t111 * rSges(4,1) - t110 * rSges(4,2) - rSges(4,3) * t131;
t126 = (-t112 - t89) * t117;
t125 = (-t112 - t135) * t117;
t105 = t119 * rSges(3,3) + (rSges(3,1) * t122 + rSges(3,2) * t124) * t117;
t104 = Icges(3,5) * t119 + (Icges(3,1) * t122 + Icges(3,4) * t124) * t117;
t103 = Icges(3,6) * t119 + (Icges(3,4) * t122 + Icges(3,2) * t124) * t117;
t102 = Icges(3,3) * t119 + (Icges(3,5) * t122 + Icges(3,6) * t124) * t117;
t92 = t119 * t94;
t88 = Icges(4,1) * t111 - Icges(4,4) * t110 - Icges(4,5) * t131;
t87 = Icges(4,4) * t111 - Icges(4,2) * t110 - Icges(4,6) * t131;
t86 = Icges(4,5) * t111 - Icges(4,6) * t110 - Icges(4,3) * t131;
t85 = t109 * rSges(3,1) - t108 * rSges(3,2) + rSges(3,3) * t133;
t84 = t107 * rSges(3,1) - t106 * rSges(3,2) - rSges(3,3) * t130;
t83 = Icges(3,1) * t109 - Icges(3,4) * t108 + Icges(3,5) * t133;
t82 = Icges(3,1) * t107 - Icges(3,4) * t106 - Icges(3,5) * t130;
t81 = Icges(3,4) * t109 - Icges(3,2) * t108 + Icges(3,6) * t133;
t80 = Icges(3,4) * t107 - Icges(3,2) * t106 - Icges(3,6) * t130;
t79 = Icges(3,5) * t109 - Icges(3,6) * t108 + Icges(3,3) * t133;
t78 = Icges(3,5) * t107 - Icges(3,6) * t106 - Icges(3,3) * t130;
t71 = -t105 * t130 - t119 * t84;
t70 = -t105 * t133 + t119 * t85;
t65 = t99 * rSges(4,1) - t98 * rSges(4,2) + t108 * rSges(4,3);
t64 = t97 * rSges(4,1) - t96 * rSges(4,2) + t106 * rSges(4,3);
t63 = Icges(4,1) * t99 - Icges(4,4) * t98 + Icges(4,5) * t108;
t62 = Icges(4,1) * t97 - Icges(4,4) * t96 + Icges(4,5) * t106;
t61 = Icges(4,4) * t99 - Icges(4,2) * t98 + Icges(4,6) * t108;
t60 = Icges(4,4) * t97 - Icges(4,2) * t96 + Icges(4,6) * t106;
t59 = Icges(4,5) * t99 - Icges(4,6) * t98 + Icges(4,3) * t108;
t58 = Icges(4,5) * t97 - Icges(4,6) * t96 + Icges(4,3) * t106;
t57 = (t116 * t84 + t118 * t85) * t117;
t48 = -t108 * t89 - t65 * t131;
t47 = t106 * t89 + t64 * t131;
t46 = -t110 * t87 + t111 * t88 - t86 * t131;
t45 = -t106 * t65 + t108 * t64;
t44 = (-t64 - t93) * t119 + t118 * t126;
t43 = t116 * t126 + t119 * t65 + t92;
t42 = t108 * t86 - t98 * t87 + t99 * t88;
t41 = t106 * t86 - t96 * t87 + t97 * t88;
t40 = (t116 * t64 + t118 * t65) * t117 + t134;
t39 = t110 * t56 - t98 * t69;
t38 = -t110 * t55 + t96 * t69;
t37 = -t110 * t61 + t111 * t63 - t59 * t131;
t36 = -t110 * t60 + t111 * t62 - t58 * t131;
t35 = t100 * t67 + t101 * t68 + t110 * t66;
t34 = t108 * t59 - t98 * t61 + t99 * t63;
t33 = t108 * t58 - t98 * t60 + t99 * t62;
t32 = t106 * t59 - t96 * t61 + t97 * t63;
t31 = t106 * t58 - t96 * t60 + t97 * t62;
t30 = t98 * t55 - t96 * t56;
t29 = -t135 * t108 - t136 * t131;
t28 = t135 * t106 + t137 * t131;
t26 = t96 * t66 + t74 * t67 + t75 * t68;
t25 = (-t93 - t137) * t119 + t118 * t125;
t24 = t116 * t125 + t136 * t119 + t92;
t23 = -t136 * t106 + t137 * t108;
t22 = t100 * t52 + t101 * t54 + t110 * t50;
t21 = t100 * t51 + t101 * t53 + t110 * t49;
t20 = (t137 * t116 + t136 * t118) * t117 + t134;
t17 = t96 * t50 + t74 * t52 + t75 * t54;
t16 = t96 * t49 + t74 * t51 + t75 * t53;
t15 = t46 * t119 + (t116 * t37 - t118 * t36) * t117;
t14 = t36 * t106 + t37 * t108 - t46 * t131;
t13 = t42 * t119 + (t116 * t34 - t118 * t33) * t117;
t12 = t41 * t119 + (t116 * t32 - t118 * t31) * t117;
t11 = t33 * t106 + t34 * t108 - t42 * t131;
t10 = t31 * t106 + t32 * t108 - t41 * t131;
t9 = t35 * t119 + (t116 * t22 - t118 * t21) * t117;
t8 = t21 * t106 + t22 * t108 - t35 * t131;
t7 = t35 * t110 + t21 * t96 + t22 * t98;
t6 = t27 * t119 + (t116 * t19 - t118 * t18) * t117;
t5 = t26 * t119 + (t116 * t17 - t118 * t16) * t117;
t4 = t18 * t106 + t19 * t108 - t27 * t131;
t3 = t16 * t106 + t17 * t108 - t26 * t131;
t1 = t26 * t110 + t16 * t96 + t17 * t98;
t72 = [m(2) + m(3) + m(4) + m(5); m(3) * t57 + m(4) * t40 + m(5) * t20; m(5) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(4) * (t40 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(3) * (t57 ^ 2 + t70 ^ 2 + t71 ^ 2) + (t6 + t13 + (-t108 * t81 + t109 * t83 + t79 * t133) * t133) * t133 + (-t5 - t12 + (-t106 * t80 + t107 * t82 - t78 * t130) * t130 + (t106 * t81 - t107 * t83 + t108 * t80 - t109 * t82 + t79 * t130 - t78 * t133) * t133) * t130 + ((t102 * t133 - t108 * t103 + t109 * t104) * t133 - (-t102 * t130 - t106 * t103 + t107 * t104) * t130 + t9 + t15 + ((t122 * t83 + t124 * t81) * t116 - (t122 * t82 + t124 * t80) * t118) * t117 ^ 2 + ((t103 * t124 + t104 * t122 + t116 * t79 - t118 * t78) * t117 + t119 * t102) * t119) * t119; m(4) * t45 + m(5) * t23; (t8 / 0.2e1 + t14 / 0.2e1) * t119 + (t6 / 0.2e1 + t13 / 0.2e1) * t108 + (t5 / 0.2e1 + t12 / 0.2e1) * t106 + m(5) * (t20 * t23 + t24 * t29 + t25 * t28) + m(4) * (t45 * t40 + t48 * t43 + t47 * t44) + ((-t9 / 0.2e1 - t15 / 0.2e1) * t124 + (-t3 / 0.2e1 - t10 / 0.2e1) * t118 + (t4 / 0.2e1 + t11 / 0.2e1) * t116) * t117; (-t14 - t8) * t131 + (t4 + t11) * t108 + (t3 + t10) * t106 + m(5) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(4) * (t45 ^ 2 + t47 ^ 2 + t48 ^ 2); m(5) * t30; t9 * t139 + t119 * t7 / 0.2e1 + t5 * t141 + t6 * t140 + m(5) * (t20 * t30 + t24 * t39 + t25 * t38) + (t116 * t142 - t118 * t1 / 0.2e1) * t117; -t7 * t131 / 0.2e1 + m(5) * (t23 * t30 + t28 * t38 + t29 * t39) + t8 * t139 + t4 * t140 + t108 * t142 + t106 * t1 / 0.2e1 + t3 * t141; m(5) * (t30 ^ 2 + t38 ^ 2 + t39 ^ 2) + t98 * t2 + t96 * t1 + t110 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t72(1), t72(2), t72(4), t72(7); t72(2), t72(3), t72(5), t72(8); t72(4), t72(5), t72(6), t72(9); t72(7), t72(8), t72(9), t72(10);];
Mq = res;
