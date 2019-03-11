% Calculate Gravitation load on the joints for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% rSges [8x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [7x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [8x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [8,3]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [8x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:29:22
% EndTime: 2019-03-10 06:29:32
% DurationCPUTime: 3.68s
% Computational Cost: add. (1120->316), mult. (3014->522), div. (0->0), fcn. (3749->14), ass. (0->134)
t84 = sin(qJ(2));
t85 = sin(qJ(1));
t127 = t84 * t85;
t90 = cos(qJ(3));
t91 = cos(qJ(2));
t119 = t90 * t91;
t83 = sin(qJ(3));
t92 = cos(qJ(1));
t66 = t119 * t85 + t83 * t92;
t82 = sin(qJ(4));
t89 = cos(qJ(4));
t46 = t127 * t82 + t66 * t89;
t116 = t92 * t90;
t123 = t85 * t91;
t65 = t123 * t83 - t116;
t81 = sin(qJ(5));
t88 = cos(qJ(5));
t23 = t46 * t81 + t65 * t88;
t24 = t46 * t88 - t65 * t81;
t126 = t84 * t89;
t45 = -t126 * t85 + t66 * t82;
t80 = sin(qJ(6));
t87 = cos(qJ(6));
t4 = t24 * t87 + t45 * t80;
t79 = sin(qJ(7));
t86 = cos(qJ(7));
t158 = t23 * t86 + t4 * t79;
t157 = t23 * t79 - t4 * t86;
t5 = t24 * t80 - t45 * t87;
t152 = g(1) * t92 + g(2) * t85;
t149 = g(3) * t84;
t148 = t46 * pkin(3);
t124 = t84 * t92;
t70 = t116 * t91 - t83 * t85;
t52 = t124 * t82 + t70 * t89;
t147 = t52 * pkin(3);
t118 = t91 * t82;
t125 = t84 * t90;
t64 = t125 * t89 - t118;
t146 = t64 * pkin(3);
t145 = t84 * pkin(2);
t144 = -rSges(4,3) - pkin(2);
t143 = rSges(6,3) + pkin(3);
t142 = pkin(4) + rSges(8,3);
t141 = rSges(7,3) * t81;
t139 = t65 * t82;
t117 = t91 * t92;
t69 = -t117 * t83 - t85 * t90;
t137 = t69 * t82;
t136 = t79 * t81;
t135 = t79 * t87;
t134 = t80 * t82;
t133 = t80 * t88;
t132 = t81 * t86;
t131 = t81 * t89;
t130 = t82 * t87;
t129 = t83 * t84;
t128 = t83 * t91;
t122 = t86 * t87;
t121 = t87 * t88;
t120 = t88 * t89;
t115 = pkin(2) * t124;
t114 = pkin(2) * t123;
t113 = pkin(2) * t117;
t112 = t81 * t129;
t111 = t82 * t129;
t110 = t88 * t129;
t109 = t83 * t124;
t108 = t144 * t84;
t74 = pkin(2) * t127;
t107 = -pkin(3) * t45 + t74;
t106 = pkin(3) * t111;
t67 = t118 * t90 - t126;
t105 = t67 * pkin(3) - t145;
t104 = rSges(3,1) * t91 - rSges(3,2) * t84;
t102 = rSges(4,1) * t90 - rSges(4,2) * t83;
t101 = rSges(5,1) * t89 - rSges(5,2) * t82;
t100 = -rSges(6,1) * t88 + rSges(6,2) * t81;
t99 = -rSges(7,1) * t87 + rSges(7,2) * t80;
t98 = rSges(8,1) * t86 - rSges(8,2) * t79;
t51 = -t124 * t89 + t70 * t82;
t97 = t51 * pkin(3) - t115;
t63 = t125 * t82 + t89 * t91;
t55 = t63 * t85;
t96 = -t55 * pkin(3) - t114;
t57 = t63 * t92;
t95 = -t57 * pkin(3) - t113;
t94 = rSges(5,3) * t129 - pkin(2) * t91;
t68 = t119 * t89 + t82 * t84;
t60 = pkin(3) * t137;
t59 = pkin(3) * t139;
t58 = t64 * t92;
t56 = t64 * t85;
t54 = (-t120 * t83 - t81 * t90) * t84;
t53 = t112 * t89 - t125 * t88;
t50 = -t128 * t81 + t68 * t88;
t49 = t128 * t88 + t68 * t81;
t44 = t64 * t88 - t112;
t43 = t64 * t81 + t110;
t42 = t109 * t81 - t58 * t88;
t41 = -t109 * t88 - t58 * t81;
t40 = t112 * t85 - t56 * t88;
t39 = -t110 * t85 - t56 * t81;
t38 = -t111 * t80 + t54 * t87;
t37 = -t111 * t87 - t54 * t80;
t36 = t120 * t69 - t70 * t81;
t35 = t131 * t69 + t70 * t88;
t34 = -t120 * t65 - t66 * t81;
t33 = -t131 * t65 + t66 * t88;
t32 = -t121 * t63 + t64 * t80;
t31 = t133 * t63 + t64 * t87;
t30 = t52 * t88 + t69 * t81;
t29 = t52 * t81 - t69 * t88;
t28 = t50 * t87 + t67 * t80;
t27 = -t50 * t80 + t67 * t87;
t22 = t44 * t87 + t63 * t80;
t21 = -t44 * t80 + t63 * t87;
t20 = t42 * t87 - t57 * t80;
t19 = -t42 * t80 - t57 * t87;
t18 = t40 * t87 - t55 * t80;
t17 = -t40 * t80 - t55 * t87;
t16 = t134 * t69 + t36 * t87;
t15 = t130 * t69 - t36 * t80;
t14 = -t134 * t65 + t34 * t87;
t13 = -t130 * t65 - t34 * t80;
t12 = -t121 * t51 + t52 * t80;
t11 = t133 * t51 + t52 * t87;
t10 = -t121 * t45 + t46 * t80;
t9 = t133 * t45 + t46 * t87;
t8 = t30 * t87 + t51 * t80;
t7 = -t30 * t80 + t51 * t87;
t2 = -t29 * t79 + t8 * t86;
t1 = -t29 * t86 - t79 * t8;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t85 - rSges(2,2) * t92) + g(2) * (rSges(2,1) * t92 - rSges(2,2) * t85)) - m(3) * (g(1) * (rSges(3,3) * t92 - t104 * t85) + g(2) * (rSges(3,3) * t85 + t104 * t92)) - m(4) * (g(1) * (-rSges(4,1) * t66 + rSges(4,2) * t65 + rSges(4,3) * t127 + t74) + g(2) * (rSges(4,1) * t70 + rSges(4,2) * t69 + t108 * t92)) - m(5) * (g(1) * (-rSges(5,1) * t46 + rSges(5,2) * t45 + rSges(5,3) * t65 + t74) + g(2) * (rSges(5,1) * t52 - rSges(5,2) * t51 + rSges(5,3) * t69 - t115)) - m(6) * (g(1) * (-rSges(6,1) * t24 + rSges(6,2) * t23 - t143 * t45 + t74) + g(2) * (rSges(6,1) * t30 - rSges(6,2) * t29 + t143 * t51 - t115)) - m(7) * (g(1) * (-rSges(7,1) * t4 + rSges(7,2) * t5 - rSges(7,3) * t23 + t107) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + rSges(7,3) * t29 + t97)) - m(8) * (g(1) * (rSges(8,1) * t157 + rSges(8,2) * t158 + t142 * t5 + t107) + g(2) * (t2 * rSges(8,1) + t1 * rSges(8,2) + t142 * t7 + t97)) -m(3) * (g(3) * t104 + t152 * (-rSges(3,1) * t84 - rSges(3,2) * t91)) - m(4) * (g(3) * (t102 * t91 + t108) + t152 * (-t102 * t84 + t144 * t91)) - m(5) * (g(1) * (-rSges(5,1) * t58 + rSges(5,2) * t57 + t92 * t94) + g(2) * (-rSges(5,1) * t56 + rSges(5,2) * t55 + t85 * t94) + g(3) * (rSges(5,1) * t68 - rSges(5,2) * t67 - rSges(5,3) * t128 - t145)) - m(6) * (g(1) * (rSges(6,1) * t42 - rSges(6,2) * t41 - t143 * t57 - t113) + g(2) * (rSges(6,1) * t40 - rSges(6,2) * t39 - t143 * t55 - t114) + g(3) * (rSges(6,1) * t50 - rSges(6,2) * t49 + t143 * t67 - t145)) - m(7) * (g(1) * (rSges(7,1) * t20 + rSges(7,2) * t19 + rSges(7,3) * t41 + t95) + g(2) * (rSges(7,1) * t18 + rSges(7,2) * t17 + rSges(7,3) * t39 + t96) + g(3) * (rSges(7,1) * t28 + rSges(7,2) * t27 + rSges(7,3) * t49 + t105)) - m(8) * (g(1) * ((t20 * t86 - t41 * t79) * rSges(8,1) + (-t20 * t79 - t41 * t86) * rSges(8,2) + t142 * t19 + t95) + g(2) * ((t18 * t86 - t39 * t79) * rSges(8,1) + (-t18 * t79 - t39 * t86) * rSges(8,2) + t142 * t17 + t96) + g(3) * ((t28 * t86 - t49 * t79) * rSges(8,1) + (-t28 * t79 - t49 * t86) * rSges(8,2) + t142 * t27 + t105)) -m(4) * (g(1) * (rSges(4,1) * t69 - rSges(4,2) * t70) + g(2) * (-rSges(4,1) * t65 - rSges(4,2) * t66) + (-rSges(4,1) * t83 - rSges(4,2) * t90) * t149) - m(5) * (g(1) * (-rSges(5,3) * t70 + t101 * t69) + g(2) * (-rSges(5,3) * t66 - t101 * t65) + (-rSges(5,3) * t90 - t101 * t83) * t149) - m(6) * (g(1) * (rSges(6,1) * t36 - rSges(6,2) * t35 + rSges(6,3) * t137 + t60) + g(2) * (rSges(6,1) * t34 - rSges(6,2) * t33 - rSges(6,3) * t139 - t59) + g(3) * (rSges(6,1) * t54 + rSges(6,2) * t53 - t111 * t143)) - m(7) * (g(1) * (rSges(7,1) * t16 + rSges(7,2) * t15 + rSges(7,3) * t35 + t60) + g(2) * (rSges(7,1) * t14 + rSges(7,2) * t13 + rSges(7,3) * t33 - t59) + g(3) * (rSges(7,1) * t38 + rSges(7,2) * t37 - rSges(7,3) * t53 - t106)) - m(8) * (g(1) * (t60 + (t16 * t86 - t35 * t79) * rSges(8,1) + (-t16 * t79 - t35 * t86) * rSges(8,2) + t142 * t15) + g(2) * (-t59 + (t14 * t86 - t33 * t79) * rSges(8,1) + (-t14 * t79 - t33 * t86) * rSges(8,2) + t142 * t13) + g(3) * (-t106 + (t38 * t86 + t53 * t79) * rSges(8,1) + (-t38 * t79 + t53 * t86) * rSges(8,2) + t142 * t37)) -m(5) * (g(1) * (-rSges(5,1) * t51 - rSges(5,2) * t52) + g(2) * (-rSges(5,1) * t45 - rSges(5,2) * t46) + g(3) * (-rSges(5,1) * t63 - rSges(5,2) * t64)) - m(6) * (g(1) * (t100 * t51 + t143 * t52) + g(2) * (t100 * t45 + t143 * t46) + g(3) * (t100 * t63 + t143 * t64)) - m(7) * (g(1) * (rSges(7,1) * t12 + rSges(7,2) * t11 - t141 * t51 + t147) + g(2) * (rSges(7,1) * t10 + rSges(7,2) * t9 - t141 * t45 + t148) + g(3) * (rSges(7,1) * t32 + rSges(7,2) * t31 - t141 * t63 + t146)) - m(8) * (g(1) * (t147 + (t12 * t86 + t136 * t51) * rSges(8,1) + (-t12 * t79 + t132 * t51) * rSges(8,2) + t142 * t11) + g(2) * (t148 + (t10 * t86 + t136 * t45) * rSges(8,1) + (-t10 * t79 + t132 * t45) * rSges(8,2) + t142 * t9) + g(3) * (t146 + (t136 * t63 + t32 * t86) * rSges(8,1) + (t132 * t63 - t32 * t79) * rSges(8,2) + t142 * t31)) -m(6) * (g(1) * (-rSges(6,1) * t29 - rSges(6,2) * t30) + g(2) * (-rSges(6,1) * t23 - rSges(6,2) * t24) + g(3) * (-rSges(6,1) * t43 - rSges(6,2) * t44)) - m(7) * (g(1) * (rSges(7,3) * t30 + t29 * t99) + g(2) * (rSges(7,3) * t24 + t23 * t99) + g(3) * (rSges(7,3) * t44 + t43 * t99)) + (-g(1) * ((-t122 * t29 - t30 * t79) * rSges(8,1) + (t135 * t29 - t30 * t86) * rSges(8,2)) - g(2) * ((-t122 * t23 - t24 * t79) * rSges(8,1) + (t135 * t23 - t24 * t86) * rSges(8,2)) - g(3) * ((-t122 * t43 - t44 * t79) * rSges(8,1) + (t135 * t43 - t44 * t86) * rSges(8,2)) - (g(1) * t29 + g(2) * t23 + g(3) * t43) * t80 * t142) * m(8), -m(7) * (g(1) * (rSges(7,1) * t7 - rSges(7,2) * t8) + g(2) * (-rSges(7,1) * t5 - rSges(7,2) * t4) + g(3) * (rSges(7,1) * t21 - rSges(7,2) * t22)) - m(8) * (g(1) * (-t142 * t8 + t7 * t98) + g(2) * (-t142 * t4 - t5 * t98) + g(3) * (-t142 * t22 + t21 * t98)) -m(8) * (g(1) * (rSges(8,1) * t1 - rSges(8,2) * t2) + g(2) * (-rSges(8,1) * t158 + rSges(8,2) * t157) + g(3) * ((-t22 * t79 - t43 * t86) * rSges(8,1) + (-t22 * t86 + t43 * t79) * rSges(8,2)))];
taug  = t3(:);
