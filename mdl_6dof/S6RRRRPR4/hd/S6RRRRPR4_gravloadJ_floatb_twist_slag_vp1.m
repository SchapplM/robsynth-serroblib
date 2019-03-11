% Calculate Gravitation load on the joints for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:13
% EndTime: 2019-03-09 22:06:16
% DurationCPUTime: 0.99s
% Computational Cost: add. (621->174), mult. (597->233), div. (0->0), fcn. (553->12), ass. (0->90)
t130 = rSges(5,3) + pkin(9);
t60 = qJ(4) + pkin(11);
t53 = cos(t60);
t66 = cos(qJ(4));
t57 = t66 * pkin(4);
t28 = pkin(5) * t53 + t57;
t24 = pkin(3) + t28;
t61 = qJ(2) + qJ(3);
t55 = sin(t61);
t56 = cos(t61);
t129 = t55 * rSges(7,3) + t56 * t24;
t50 = t57 + pkin(3);
t128 = t55 * rSges(6,3) + t56 * t50;
t127 = t56 * rSges(4,1) - rSges(4,2) * t55;
t65 = sin(qJ(1));
t68 = cos(qJ(1));
t126 = g(1) * t68 + g(2) * t65;
t77 = t56 * pkin(3) + t130 * t55;
t125 = t126 * t55;
t107 = t56 * t65;
t54 = qJ(6) + t60;
t46 = sin(t54);
t47 = cos(t54);
t6 = t107 * t46 + t47 * t68;
t7 = -t107 * t47 + t46 * t68;
t124 = -t6 * rSges(7,1) + t7 * rSges(7,2);
t106 = t56 * t68;
t8 = -t106 * t46 + t47 * t65;
t9 = t106 * t47 + t46 * t65;
t123 = t8 * rSges(7,1) - t9 * rSges(7,2);
t64 = sin(qJ(2));
t122 = pkin(2) * t64;
t63 = sin(qJ(4));
t121 = pkin(4) * t63;
t118 = g(3) * t55;
t117 = rSges(3,3) + pkin(7);
t116 = rSges(5,1) * t66;
t115 = rSges(6,1) * t53;
t114 = rSges(7,1) * t47;
t112 = rSges(5,2) * t63;
t52 = sin(t60);
t111 = rSges(6,2) * t52;
t110 = rSges(7,2) * t46;
t62 = -qJ(5) - pkin(9);
t59 = -pkin(10) + t62;
t109 = t55 * t59;
t108 = t55 * t62;
t105 = t63 * t65;
t104 = t63 * t68;
t103 = t65 * t66;
t102 = t66 * t68;
t69 = -pkin(8) - pkin(7);
t101 = rSges(4,3) - t69;
t93 = t55 * t110;
t100 = rSges(7,3) * t107 + t65 * t93;
t99 = rSges(7,3) * t106 + t68 * t93;
t94 = t55 * t111;
t98 = rSges(6,3) * t107 + t65 * t94;
t97 = rSges(6,3) * t106 + t68 * t94;
t27 = pkin(5) * t52 + t121;
t96 = t27 - t69;
t95 = t55 * t112;
t92 = t107 * t130 + t65 * t95;
t91 = t106 * t130 + t68 * t95;
t89 = -t69 + t121;
t88 = -t50 - t115;
t87 = -t24 - t114;
t67 = cos(qJ(2));
t85 = rSges(3,1) * t67 - rSges(3,2) * t64;
t82 = -rSges(4,1) * t55 - rSges(4,2) * t56;
t81 = -rSges(7,1) * t46 - rSges(7,2) * t47;
t80 = pkin(1) + t85;
t17 = -t104 * t56 + t103;
t15 = t105 * t56 + t102;
t79 = t128 + (-t111 + t115) * t56;
t78 = t129 + (-t110 + t114) * t56;
t76 = t77 + (-t112 + t116) * t56;
t74 = -t108 + t128;
t73 = -t109 + t129;
t70 = (-pkin(3) - t116) * t125;
t58 = t67 * pkin(2);
t51 = t58 + pkin(1);
t31 = t68 * t51;
t18 = t102 * t56 + t105;
t16 = -t103 * t56 + t104;
t13 = t106 * t53 + t52 * t65;
t12 = -t106 * t52 + t53 * t65;
t11 = -t107 * t53 + t52 * t68;
t10 = t107 * t52 + t53 * t68;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t65 - rSges(2,2) * t68) + g(2) * (rSges(2,1) * t68 - rSges(2,2) * t65)) - m(3) * ((g(1) * t117 + g(2) * t80) * t68 + (-g(1) * t80 + g(2) * t117) * t65) - m(4) * (g(2) * t31 + (g(1) * t101 + g(2) * t127) * t68 + (g(1) * (-t51 - t127) + g(2) * t101) * t65) - m(5) * (g(1) * (t16 * rSges(5,1) + t15 * rSges(5,2)) + g(2) * (t18 * rSges(5,1) + t17 * rSges(5,2) + t31) + (-g(1) * t69 + g(2) * t77) * t68 + (g(1) * (-t51 - t77) - g(2) * t69) * t65) - m(6) * (g(1) * (t11 * rSges(6,1) + t10 * rSges(6,2)) + g(2) * (t13 * rSges(6,1) + t12 * rSges(6,2) + t31) + (g(1) * t89 + g(2) * t74) * t68 + (g(1) * (-t51 - t74) + g(2) * t89) * t65) - m(7) * (g(1) * (t7 * rSges(7,1) + t6 * rSges(7,2)) + g(2) * (t9 * rSges(7,1) + t8 * rSges(7,2) + t31) + (g(1) * t96 + g(2) * t73) * t68 + (g(1) * (-t51 - t73) + g(2) * t96) * t65) -m(3) * (g(3) * t85 + t126 * (-rSges(3,1) * t64 - rSges(3,2) * t67)) - m(4) * (g(3) * (t58 + t127) + t126 * (t82 - t122)) - m(5) * (g(1) * (-t122 * t68 + t91) + g(2) * (-t122 * t65 + t92) + g(3) * (t58 + t76) + t70) - m(6) * (g(1) * t97 + g(2) * t98 + g(3) * (t58 + t79 - t108) + t126 * (t55 * t88 - t56 * t62 - t122)) - m(7) * (g(1) * t99 + g(2) * t100 + g(3) * (t58 + t78 - t109) + t126 * (t55 * t87 - t56 * t59 - t122)) -m(4) * (g(3) * t127 + t126 * t82) - m(5) * (g(1) * t91 + g(2) * t92 + g(3) * t76 + t70) - m(6) * (g(1) * (-t106 * t62 + t97) + g(2) * (-t107 * t62 + t98) + g(3) * t79 + (-g(3) * t62 + t126 * t88) * t55) - m(7) * (g(1) * (-t106 * t59 + t99) + g(2) * (-t107 * t59 + t100) + g(3) * t78 + (-g(3) * t59 + t126 * t87) * t55) -m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (-rSges(5,1) * t15 + rSges(5,2) * t16)) - m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t13 + pkin(4) * t17) + g(2) * (-rSges(6,1) * t10 + rSges(6,2) * t11 - pkin(4) * t15)) - m(7) * (g(1) * (-t106 * t27 + t28 * t65 + t123) + g(2) * (-t107 * t27 - t28 * t68 + t124)) + (-m(5) * (-rSges(5,1) * t63 - rSges(5,2) * t66) - m(6) * (-rSges(6,1) * t52 - rSges(6,2) * t53 - t121) - m(7) * (-t27 + t81)) * t118 (-m(6) - m(7)) * (-g(3) * t56 + t125) -m(7) * (g(1) * t123 + g(2) * t124 + t81 * t118)];
taug  = t1(:);
