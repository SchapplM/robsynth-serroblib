% Calculate Gravitation load on the joints for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:18
% EndTime: 2019-03-10 03:39:20
% DurationCPUTime: 1.03s
% Computational Cost: add. (684->180), mult. (645->242), div. (0->0), fcn. (603->12), ass. (0->94)
t137 = rSges(5,3) + pkin(9);
t63 = qJ(4) + qJ(5);
t57 = cos(t63);
t68 = cos(qJ(4));
t60 = t68 * pkin(4);
t31 = pkin(5) * t57 + t60;
t27 = pkin(3) + t31;
t64 = qJ(2) + qJ(3);
t56 = sin(t64);
t58 = cos(t64);
t136 = rSges(7,3) * t56 + t58 * t27;
t53 = t60 + pkin(3);
t135 = rSges(6,3) * t56 + t58 * t53;
t134 = rSges(4,1) * t58 - rSges(4,2) * t56;
t67 = sin(qJ(1));
t70 = cos(qJ(1));
t133 = g(1) * t70 + g(2) * t67;
t80 = t58 * pkin(3) + t137 * t56;
t71 = -pkin(10) - pkin(9);
t113 = t58 * t67;
t59 = qJ(6) + t63;
t51 = sin(t59);
t52 = cos(t59);
t5 = t113 * t51 + t52 * t70;
t6 = -t113 * t52 + t51 * t70;
t132 = -rSges(7,1) * t5 + rSges(7,2) * t6;
t112 = t58 * t70;
t7 = -t112 * t51 + t52 * t67;
t8 = t112 * t52 + t51 * t67;
t131 = rSges(7,1) * t7 - rSges(7,2) * t8;
t66 = sin(qJ(2));
t130 = pkin(2) * t66;
t65 = sin(qJ(4));
t129 = pkin(4) * t65;
t55 = sin(t63);
t128 = pkin(5) * t55;
t125 = g(3) * t56;
t124 = rSges(3,3) + pkin(7);
t13 = t113 * t55 + t57 * t70;
t14 = -t113 * t57 + t55 * t70;
t123 = -rSges(6,1) * t13 + rSges(6,2) * t14;
t122 = rSges(5,1) * t68;
t121 = rSges(6,1) * t57;
t120 = rSges(7,1) * t52;
t118 = rSges(5,2) * t65;
t117 = rSges(6,2) * t55;
t116 = rSges(7,2) * t51;
t62 = -pkin(11) + t71;
t115 = t56 * t62;
t114 = t56 * t71;
t111 = t58 * t71;
t110 = t65 * t67;
t109 = t65 * t70;
t108 = t67 * t68;
t107 = t68 * t70;
t72 = -pkin(8) - pkin(7);
t106 = rSges(4,3) - t72;
t15 = -t112 * t55 + t57 * t67;
t16 = t112 * t57 + t55 * t67;
t105 = rSges(6,1) * t15 - rSges(6,2) * t16;
t97 = t56 * t116;
t104 = rSges(7,3) * t113 + t67 * t97;
t103 = rSges(7,3) * t112 + t70 * t97;
t98 = t56 * t117;
t102 = rSges(6,3) * t113 + t67 * t98;
t101 = rSges(6,3) * t112 + t70 * t98;
t30 = t128 + t129;
t100 = t30 - t72;
t99 = t56 * t118;
t96 = t113 * t137 + t67 * t99;
t95 = t112 * t137 + t70 * t99;
t93 = -t72 + t129;
t92 = -t53 - t121;
t91 = -t27 - t120;
t69 = cos(qJ(2));
t89 = rSges(3,1) * t69 - rSges(3,2) * t66;
t86 = -rSges(4,1) * t56 - rSges(4,2) * t58;
t85 = -rSges(6,1) * t55 - rSges(6,2) * t57;
t84 = -rSges(7,1) * t51 - rSges(7,2) * t52;
t83 = pkin(1) + t89;
t19 = -t109 * t58 + t108;
t17 = t110 * t58 + t107;
t82 = t135 + (-t117 + t121) * t58;
t81 = t136 + (-t116 + t120) * t58;
t79 = t80 + (-t118 + t122) * t58;
t77 = -t114 + t135;
t76 = -t115 + t136;
t73 = t133 * (-pkin(3) - t122) * t56;
t61 = t69 * pkin(2);
t54 = t61 + pkin(1);
t34 = t70 * t54;
t20 = t107 * t58 + t110;
t18 = -t108 * t58 + t109;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t67 - rSges(2,2) * t70) + g(2) * (rSges(2,1) * t70 - rSges(2,2) * t67)) - m(3) * ((g(1) * t124 + g(2) * t83) * t70 + (-g(1) * t83 + g(2) * t124) * t67) - m(4) * (g(2) * t34 + (g(1) * t106 + g(2) * t134) * t70 + (g(1) * (-t54 - t134) + g(2) * t106) * t67) - m(5) * (g(1) * (t18 * rSges(5,1) + t17 * rSges(5,2)) + g(2) * (t20 * rSges(5,1) + t19 * rSges(5,2) + t34) + (-g(1) * t72 + g(2) * t80) * t70 + (g(1) * (-t54 - t80) - g(2) * t72) * t67) - m(6) * (g(1) * (t14 * rSges(6,1) + t13 * rSges(6,2)) + g(2) * (t16 * rSges(6,1) + t15 * rSges(6,2) + t34) + (g(1) * t93 + g(2) * t77) * t70 + (g(1) * (-t54 - t77) + g(2) * t93) * t67) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2)) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t34) + (g(1) * t100 + g(2) * t76) * t70 + (g(1) * (-t54 - t76) + g(2) * t100) * t67) -m(3) * (g(3) * t89 + t133 * (-rSges(3,1) * t66 - rSges(3,2) * t69)) - m(4) * (g(3) * (t61 + t134) + t133 * (t86 - t130)) - m(5) * (g(1) * (-t130 * t70 + t95) + g(2) * (-t130 * t67 + t96) + g(3) * (t61 + t79) + t73) - m(6) * (g(1) * t101 + g(2) * t102 + g(3) * (t61 + t82 - t114) + t133 * (t56 * t92 - t111 - t130)) - m(7) * (g(1) * t103 + g(2) * t104 + g(3) * (t61 + t81 - t115) + t133 * (t56 * t91 - t58 * t62 - t130)) -m(4) * (g(3) * t134 + t133 * t86) - m(5) * (g(1) * t95 + g(2) * t96 + g(3) * t79 + t73) - m(6) * (g(1) * (-t111 * t70 + t101) + g(2) * (-t111 * t67 + t102) + g(3) * t82 + (-g(3) * t71 + t133 * t92) * t56) - m(7) * (g(1) * (-t112 * t62 + t103) + g(2) * (-t113 * t62 + t104) + g(3) * t81 + (-g(3) * t62 + t133 * t91) * t56) -m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (-rSges(5,1) * t17 + rSges(5,2) * t18)) - m(6) * (g(1) * (pkin(4) * t19 + t105) + g(2) * (-pkin(4) * t17 + t123)) - m(7) * (g(1) * (-t112 * t30 + t31 * t67 + t131) + g(2) * (-t113 * t30 - t31 * t70 + t132)) + (-m(5) * (-rSges(5,1) * t65 - rSges(5,2) * t68) - m(6) * (t85 - t129) - m(7) * (-t30 + t84)) * t125, -m(6) * (g(1) * t105 + g(2) * t123) - m(7) * (g(1) * (pkin(5) * t15 + t131) + g(2) * (-pkin(5) * t13 + t132)) + (-m(6) * t85 - m(7) * (t84 - t128)) * t125, -m(7) * (g(1) * t131 + g(2) * t132 + t125 * t84)];
taug  = t1(:);
