% Calculate Gravitation load on the joints for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:43
% EndTime: 2019-03-09 17:51:47
% DurationCPUTime: 1.27s
% Computational Cost: add. (696->200), mult. (1683->279), div. (0->0), fcn. (2032->10), ass. (0->83)
t128 = pkin(4) + pkin(9);
t70 = sin(qJ(2));
t71 = sin(qJ(1));
t74 = cos(qJ(2));
t100 = cos(pkin(6));
t117 = cos(qJ(1));
t85 = t100 * t117;
t49 = t70 * t85 + t71 * t74;
t69 = sin(qJ(3));
t73 = cos(qJ(3));
t67 = sin(pkin(6));
t94 = t67 * t117;
t21 = t49 * t69 + t73 * t94;
t48 = t70 * t71 - t74 * t85;
t68 = sin(qJ(5));
t72 = cos(qJ(5));
t3 = t21 * t68 + t48 * t72;
t4 = t21 * t72 - t48 * t68;
t22 = t49 * t73 - t69 * t94;
t111 = t67 * t71;
t89 = t71 * t100;
t51 = t117 * t74 - t70 * t89;
t26 = t111 * t69 + t51 * t73;
t110 = t67 * t73;
t47 = t100 * t69 + t110 * t70;
t127 = g(1) * t26 + g(2) * t22 + g(3) * t47;
t123 = g(3) * t67;
t122 = rSges(5,1) + pkin(9);
t121 = rSges(7,1) + pkin(5);
t120 = rSges(7,2) + pkin(10);
t119 = rSges(4,3) + pkin(9);
t118 = rSges(6,3) + pkin(10);
t114 = t48 * t73;
t50 = t117 * t70 + t74 * t89;
t113 = t50 * t73;
t112 = t67 * t70;
t109 = t67 * t74;
t108 = t68 * t69;
t107 = t68 * t74;
t106 = t69 * t72;
t105 = t117 * pkin(1) + pkin(8) * t111;
t104 = pkin(2) * t109 + pkin(9) * t112;
t103 = qJ(4) * t69;
t102 = rSges(5,3) + qJ(4);
t101 = rSges(7,3) + qJ(6);
t99 = t72 * t109;
t98 = t73 * t109;
t39 = t48 * pkin(2);
t97 = -pkin(3) * t114 - t48 * t103 - t39;
t43 = t50 * pkin(2);
t96 = -pkin(3) * t113 - t50 * t103 - t43;
t95 = t51 * pkin(2) + t105;
t93 = -t71 * pkin(1) + pkin(8) * t94;
t15 = t21 * pkin(3);
t92 = -pkin(10) * t21 - t15;
t25 = -t110 * t71 + t51 * t69;
t17 = t25 * pkin(3);
t91 = -pkin(10) * t25 - t17;
t46 = -t100 * t73 + t112 * t69;
t37 = t46 * pkin(3);
t90 = -pkin(10) * t46 - t37;
t88 = t26 * pkin(3) + t95;
t87 = pkin(3) * t98 + t103 * t109 + t104;
t86 = -t49 * pkin(2) + t93;
t84 = rSges(4,1) * t73 - rSges(4,2) * t69;
t83 = rSges(5,2) * t73 - rSges(5,3) * t69;
t82 = -pkin(3) * t22 + t86;
t81 = pkin(4) * t112 + pkin(10) * t98 + t87;
t79 = -pkin(10) * t114 + t128 * t49 + t97;
t78 = -pkin(10) * t113 + t128 * t51 + t96;
t76 = qJ(4) * t25 + t128 * t50 + t88;
t75 = -qJ(4) * t21 - t128 * t48 + t82;
t30 = (t107 * t69 + t70 * t72) * t67;
t29 = t112 * t68 - t69 * t99;
t20 = -t46 * t68 + t99;
t19 = t107 * t67 + t46 * t72;
t11 = -t108 * t50 + t51 * t72;
t10 = t106 * t50 + t51 * t68;
t9 = -t108 * t48 + t49 * t72;
t8 = t106 * t48 + t49 * t68;
t7 = t25 * t68 + t50 * t72;
t6 = -t25 * t72 + t50 * t68;
t1 = [-m(2) * (g(1) * (-t71 * rSges(2,1) - rSges(2,2) * t117) + g(2) * (rSges(2,1) * t117 - t71 * rSges(2,2))) - m(3) * (g(1) * (-t49 * rSges(3,1) + t48 * rSges(3,2) + rSges(3,3) * t94 + t93) + g(2) * (rSges(3,1) * t51 - rSges(3,2) * t50 + rSges(3,3) * t111 + t105)) - m(4) * (g(1) * (-rSges(4,1) * t22 + rSges(4,2) * t21 - t119 * t48 + t86) + g(2) * (rSges(4,1) * t26 - rSges(4,2) * t25 + t119 * t50 + t95)) - m(5) * (g(1) * (rSges(5,2) * t22 - t102 * t21 - t122 * t48 + t82) + g(2) * (-rSges(5,2) * t26 + t102 * t25 + t122 * t50 + t88)) - m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4 - t118 * t22 + t75) + g(2) * (rSges(6,1) * t7 - rSges(6,2) * t6 + t118 * t26 + t76)) - m(7) * (g(1) * (t101 * t4 - t120 * t22 - t121 * t3 + t75) + g(2) * (t101 * t6 + t120 * t26 + t121 * t7 + t76)) -m(3) * (g(1) * (-rSges(3,1) * t50 - rSges(3,2) * t51) + g(2) * (-rSges(3,1) * t48 - rSges(3,2) * t49) + (rSges(3,1) * t74 - rSges(3,2) * t70) * t123) - m(4) * (g(1) * (t119 * t51 - t50 * t84 - t43) + g(2) * (t119 * t49 - t48 * t84 - t39) + g(3) * t104 + (rSges(4,3) * t70 + t74 * t84) * t123) - m(5) * (g(1) * (t122 * t51 + t83 * t50 + t96) + g(2) * (t122 * t49 + t83 * t48 + t97) + g(3) * t87 + (rSges(5,1) * t70 - t83 * t74) * t123) - m(6) * (g(1) * (rSges(6,1) * t11 - rSges(6,2) * t10 - rSges(6,3) * t113 + t78) + g(2) * (rSges(6,1) * t9 - rSges(6,2) * t8 - rSges(6,3) * t114 + t79) + g(3) * (t30 * rSges(6,1) - t29 * rSges(6,2) + rSges(6,3) * t98 + t81)) - m(7) * (g(1) * (-rSges(7,2) * t113 + t10 * t101 + t11 * t121 + t78) + g(2) * (-rSges(7,2) * t114 + t101 * t8 + t121 * t9 + t79) + g(3) * (rSges(7,2) * t98 + t101 * t29 + t121 * t30 + t81)) -m(4) * (g(1) * (-rSges(4,1) * t25 - rSges(4,2) * t26) + g(2) * (-rSges(4,1) * t21 - rSges(4,2) * t22) + g(3) * (-rSges(4,1) * t46 - rSges(4,2) * t47)) - m(5) * (g(1) * (rSges(5,2) * t25 + t102 * t26 - t17) + g(2) * (rSges(5,2) * t21 + t102 * t22 - t15) + g(3) * (rSges(5,2) * t46 + t102 * t47 - t37)) + (-g(1) * (-rSges(7,2) * t25 + t91) - g(2) * (-rSges(7,2) * t21 + t92) - g(3) * (-rSges(7,2) * t46 + t90) - t127 * (-t101 * t72 + t121 * t68 + qJ(4))) * m(7) + (-g(1) * (-rSges(6,3) * t25 + t91) - g(2) * (-rSges(6,3) * t21 + t92) - g(3) * (-rSges(6,3) * t46 + t90) - t127 * (rSges(6,1) * t68 + rSges(6,2) * t72 + qJ(4))) * m(6) (-m(5) - m(6) - m(7)) * (g(1) * t25 + g(2) * t21 + g(3) * t46) -m(6) * (g(1) * (-rSges(6,1) * t6 - rSges(6,2) * t7) + g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3) + g(3) * (rSges(6,1) * t19 + rSges(6,2) * t20)) - m(7) * (g(1) * (t101 * t7 - t121 * t6) + g(2) * (t101 * t3 + t121 * t4) + g(3) * (-t101 * t20 + t121 * t19)) -m(7) * (g(1) * t6 - g(2) * t4 - g(3) * t19)];
taug  = t1(:);
