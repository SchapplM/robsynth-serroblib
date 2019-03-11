% Calculate Gravitation load on the joints for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:29
% EndTime: 2019-03-08 22:02:33
% DurationCPUTime: 1.55s
% Computational Cost: add. (1154->220), mult. (3093->353), div. (0->0), fcn. (3967->16), ass. (0->105)
t68 = sin(pkin(6));
t69 = cos(pkin(12));
t122 = t68 * t69;
t67 = sin(pkin(7));
t108 = cos(pkin(13));
t66 = sin(pkin(13));
t73 = sin(qJ(3));
t77 = cos(qJ(3));
t82 = t108 * t77 - t73 * t66;
t46 = t82 * t67;
t70 = cos(pkin(7));
t48 = t82 * t70;
t107 = sin(pkin(12));
t74 = sin(qJ(2));
t78 = cos(qJ(2));
t109 = cos(pkin(6));
t96 = t69 * t109;
t52 = -t107 * t74 + t78 * t96;
t53 = t107 * t78 + t74 * t96;
t57 = -t108 * t73 - t77 * t66;
t12 = -t122 * t46 + t48 * t52 + t53 * t57;
t104 = t67 * t122;
t118 = t70 * t77;
t105 = pkin(3) * t118;
t126 = t53 * t73;
t79 = t52 * t105 + (-t77 * t104 - t126) * pkin(3);
t137 = t12 * pkin(4) + t79;
t112 = t77 * t78;
t116 = t73 * t74;
t136 = t70 * t112 - t116;
t128 = rSges(7,3) + pkin(11);
t135 = g(1) * t128;
t47 = t57 * t67;
t49 = t57 * t70;
t23 = t68 * (-t49 * t78 + t74 * t82) - t109 * t47;
t133 = pkin(3) * t77;
t132 = g(3) * t68;
t76 = cos(qJ(5));
t131 = t76 * pkin(5);
t129 = -rSges(6,3) - pkin(10);
t127 = rSges(5,3) * t67;
t88 = t109 * t107;
t55 = t69 * t78 - t74 * t88;
t125 = t55 * t73;
t72 = sin(qJ(5));
t124 = t67 * t72;
t123 = t67 * t76;
t121 = t68 * t74;
t120 = t68 * t78;
t119 = t70 * t73;
t71 = sin(qJ(6));
t117 = t71 * t76;
t115 = t73 * t78;
t114 = t74 * t77;
t75 = cos(qJ(6));
t113 = t75 * t76;
t50 = pkin(3) * t119 + (-pkin(9) - qJ(4)) * t67;
t65 = pkin(2) + t133;
t111 = -t53 * t50 + t52 * t65;
t54 = -t69 * t74 - t78 * t88;
t110 = -t55 * t50 + t54 * t65;
t106 = -m(5) - m(6) - m(7);
t103 = t67 * t121;
t27 = t49 * t53 + t52 * t82;
t101 = t27 * pkin(4) + t111;
t29 = t49 * t55 + t54 * t82;
t100 = t29 * pkin(4) + t110;
t99 = g(2) * t128;
t98 = g(3) * t128;
t97 = t68 * t107;
t95 = t77 * t109;
t92 = t67 * t97;
t91 = rSges(6,1) * t76 - rSges(6,2) * t72;
t89 = -pkin(3) * t125 + t54 * t105 + t92 * t133;
t34 = (t49 * t74 + t78 * t82) * t68;
t60 = t65 * t120;
t87 = t34 * pkin(4) - t50 * t121 + t60;
t85 = -t52 * t70 + t104;
t15 = t46 * t97 + t54 * t48 + t55 * t57;
t84 = t15 * pkin(4) + t89;
t83 = (t136 * t68 + t95 * t67) * pkin(3);
t22 = t109 * t46 + (t48 * t78 + t57 * t74) * t68;
t81 = t22 * pkin(4) + t83;
t80 = t54 * t70 + t92;
t11 = -t122 * t47 + t49 * t52 - t53 * t82;
t14 = t47 * t97 + t54 * t49 - t55 * t82;
t51 = t109 * t70 - t120 * t67;
t38 = -t54 * t67 + t70 * t97;
t37 = -t122 * t70 - t52 * t67;
t33 = t120 * t57 - t121 * t48;
t31 = t103 * t72 + t34 * t76;
t30 = -t103 * t76 + t34 * t72;
t28 = -t48 * t55 + t54 * t57;
t26 = -t48 * t53 + t52 * t57;
t18 = t23 * t76 + t51 * t72;
t17 = -t23 * t72 + t51 * t76;
t10 = t124 * t55 + t29 * t76;
t9 = -t123 * t55 + t29 * t72;
t8 = t124 * t53 + t27 * t76;
t7 = -t123 * t53 + t27 * t72;
t4 = -t14 * t76 + t38 * t72;
t3 = t14 * t72 + t38 * t76;
t2 = -t11 * t76 + t37 * t72;
t1 = t11 * t72 + t37 * t76;
t5 = [(-m(2) - m(3) - m(4) + t106) * g(3), -m(3) * (g(1) * (rSges(3,1) * t54 - rSges(3,2) * t55) + g(2) * (rSges(3,1) * t52 - rSges(3,2) * t53) + (rSges(3,1) * t78 - rSges(3,2) * t74) * t132) - m(4) * (g(1) * (t54 * pkin(2) + (-t119 * t55 + t54 * t77) * rSges(4,1) + (-t118 * t55 - t54 * t73) * rSges(4,2)) + g(2) * (t52 * pkin(2) + (-t119 * t53 + t52 * t77) * rSges(4,1) + (-t118 * t53 - t52 * t73) * rSges(4,2)) + (t78 * pkin(2) + (-t70 * t116 + t112) * rSges(4,1) + (-t70 * t114 - t115) * rSges(4,2)) * t132 + (g(1) * t55 + g(2) * t53 + t74 * t132) * t67 * (rSges(4,3) + pkin(9))) - m(5) * (g(1) * (rSges(5,1) * t29 + rSges(5,2) * t28 + t127 * t55 + t110) + g(2) * (rSges(5,1) * t27 + rSges(5,2) * t26 + t127 * t53 + t111) + g(3) * (rSges(5,1) * t34 + rSges(5,2) * t33 + t60 + (-t50 + t127) * t121)) - m(6) * (g(1) * (rSges(6,1) * t10 - rSges(6,2) * t9 + t129 * t28 + t100) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t129 * t26 + t101) + g(3) * (rSges(6,1) * t31 - rSges(6,2) * t30 + t129 * t33 + t87)) - m(7) * (g(1) * (t10 * pkin(5) - t28 * pkin(10) + (t10 * t75 - t28 * t71) * rSges(7,1) + (-t10 * t71 - t28 * t75) * rSges(7,2) + t128 * t9 + t100) + g(2) * (t8 * pkin(5) - t26 * pkin(10) + (-t26 * t71 + t75 * t8) * rSges(7,1) + (-t26 * t75 - t71 * t8) * rSges(7,2) + t128 * t7 + t101) + g(3) * (t31 * pkin(5) - t33 * pkin(10) + (t31 * t75 - t33 * t71) * rSges(7,1) + (-t31 * t71 - t33 * t75) * rSges(7,2) + t128 * t30 + t87)) -m(4) * (g(1) * ((t77 * t80 - t125) * rSges(4,1) + (-t55 * t77 - t73 * t80) * rSges(4,2)) + g(2) * ((-t77 * t85 - t126) * rSges(4,1) + (-t53 * t77 + t73 * t85) * rSges(4,2)) + g(3) * ((-rSges(4,2) * t109 * t73 + rSges(4,1) * t95) * t67 + (t136 * rSges(4,1) + (-t115 * t70 - t114) * rSges(4,2)) * t68)) - m(5) * (g(1) * (rSges(5,1) * t15 + rSges(5,2) * t14 + t89) + g(2) * (rSges(5,1) * t12 + rSges(5,2) * t11 + t79) + g(3) * (rSges(5,1) * t22 - rSges(5,2) * t23 + t83)) - m(6) * (g(1) * (t129 * t14 + t15 * t91 + t84) + g(2) * (t11 * t129 + t12 * t91 + t137) + g(3) * (-t129 * t23 + t22 * t91 + t81)) - m(7) * (g(1) * (t15 * t131 - t14 * pkin(10) + (t113 * t15 - t14 * t71) * rSges(7,1) + (-t117 * t15 - t14 * t75) * rSges(7,2) + t84) + g(2) * (t12 * t131 - t11 * pkin(10) + (-t11 * t71 + t113 * t12) * rSges(7,1) + (-t11 * t75 - t117 * t12) * rSges(7,2) + t137) + g(3) * (t22 * t131 + t23 * pkin(10) + (t113 * t22 + t23 * t71) * rSges(7,1) + (-t117 * t22 + t23 * t75) * rSges(7,2) + t81) + (t12 * t99 + t135 * t15 + t22 * t98) * t72) t106 * (g(1) * t38 + g(2) * t37 + g(3) * t51) -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (rSges(6,1) * t17 - rSges(6,2) * t18)) - m(7) * (t4 * t135 + t18 * t98 + t2 * t99 + (g(1) * t3 + g(2) * t1 + g(3) * t17) * (t75 * rSges(7,1) - t71 * rSges(7,2) + pkin(5))) -m(7) * (g(1) * ((-t15 * t75 - t4 * t71) * rSges(7,1) + (t15 * t71 - t4 * t75) * rSges(7,2)) + g(2) * ((-t12 * t75 - t2 * t71) * rSges(7,1) + (t12 * t71 - t2 * t75) * rSges(7,2)) + g(3) * ((-t18 * t71 - t22 * t75) * rSges(7,1) + (-t18 * t75 + t22 * t71) * rSges(7,2)))];
taug  = t5(:);
