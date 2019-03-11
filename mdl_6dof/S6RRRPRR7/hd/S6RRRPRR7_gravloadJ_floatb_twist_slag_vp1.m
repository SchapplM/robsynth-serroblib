% Calculate Gravitation load on the joints for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:49
% EndTime: 2019-03-09 18:32:52
% DurationCPUTime: 1.19s
% Computational Cost: add. (900->195), mult. (1234->275), div. (0->0), fcn. (1423->14), ass. (0->86)
t125 = pkin(11) + rSges(7,3);
t76 = sin(qJ(6));
t80 = cos(qJ(6));
t139 = rSges(7,1) * t80 - rSges(7,2) * t76 + pkin(5);
t123 = cos(qJ(1));
t74 = sin(pkin(6));
t103 = t74 * t123;
t78 = sin(qJ(2));
t79 = sin(qJ(1));
t82 = cos(qJ(2));
t105 = cos(pkin(6));
t96 = t105 * t123;
t47 = t78 * t96 + t79 * t82;
t73 = qJ(3) + pkin(12);
t69 = qJ(5) + t73;
t64 = sin(t69);
t65 = cos(t69);
t15 = -t64 * t103 + t47 * t65;
t46 = t78 * t79 - t82 * t96;
t138 = t15 * t76 - t46 * t80;
t137 = -t15 * t80 - t46 * t76;
t115 = t74 * t82;
t136 = g(3) * t115;
t135 = t76 * rSges(7,1) + t80 * rSges(7,2);
t130 = g(2) * t46;
t98 = t79 * t105;
t48 = t123 * t78 + t82 * t98;
t134 = g(1) * t48 + t130;
t133 = t125 * t64 + t139 * t65;
t101 = -t65 * t103 - t47 * t64;
t132 = rSges(6,1) * t101 - t15 * rSges(6,2);
t129 = g(2) * t47;
t68 = cos(t73);
t81 = cos(qJ(3));
t70 = t81 * pkin(3);
t56 = pkin(4) * t68 + t70;
t54 = pkin(2) + t56;
t128 = t54 * t136;
t127 = g(3) * t74;
t126 = -pkin(9) - rSges(4,3);
t117 = t74 * t79;
t49 = t123 * t82 - t78 * t98;
t18 = -t65 * t117 + t49 * t64;
t19 = t64 * t117 + t49 * t65;
t124 = -t18 * rSges(6,1) - t19 * rSges(6,2);
t118 = t74 * t78;
t116 = t74 * t81;
t75 = -qJ(4) - pkin(9);
t72 = -pkin(10) + t75;
t112 = -t46 * t54 - t47 * t72;
t111 = -t48 * t54 - t49 * t72;
t35 = t105 * t65 - t64 * t118;
t36 = t105 * t64 + t65 * t118;
t110 = t35 * rSges(6,1) - t36 * rSges(6,2);
t67 = sin(t73);
t77 = sin(qJ(3));
t55 = pkin(3) * t77 + pkin(4) * t67;
t109 = t56 * t117 - t49 * t55;
t108 = t105 * t56 - t55 * t118;
t107 = t123 * pkin(1) + pkin(8) * t117;
t106 = t75 - rSges(5,3);
t104 = t77 * t117;
t102 = -t79 * pkin(1) + pkin(8) * t103;
t100 = t67 * t103 - t47 * t68;
t59 = t77 * t103;
t99 = -t47 * t81 + t59;
t97 = t55 * t117 - t48 * t72 + t49 * t54 + t107;
t95 = rSges(6,1) * t65 - rSges(6,2) * t64;
t94 = -t56 * t103 - t47 * t55;
t93 = rSges(4,1) * t81 - rSges(4,2) * t77 + pkin(2);
t22 = t79 * t116 - t49 * t77;
t66 = t70 + pkin(2);
t91 = rSges(5,1) * t68 - rSges(5,2) * t67 + t66;
t90 = t55 * t103 + t46 * t72 - t47 * t54 + t102;
t89 = t68 * t103 + t47 * t67;
t88 = t81 * t103 + t47 * t77;
t87 = t105 * t81 - t77 * t118;
t86 = t139 * t101 + t125 * t15;
t85 = t125 * t19 - t139 * t18;
t84 = t125 * t36 + t139 * t35;
t23 = t49 * t81 + t104;
t21 = t67 * t117 + t49 * t68;
t20 = t68 * t117 - t49 * t67;
t2 = t19 * t80 + t48 * t76;
t1 = -t19 * t76 + t48 * t80;
t3 = [-m(2) * (g(1) * (-t79 * rSges(2,1) - t123 * rSges(2,2)) + g(2) * (t123 * rSges(2,1) - t79 * rSges(2,2))) - m(3) * (g(1) * (-t47 * rSges(3,1) + t46 * rSges(3,2) + rSges(3,3) * t103 + t102) + g(2) * (rSges(3,1) * t49 - rSges(3,2) * t48 + rSges(3,3) * t117 + t107)) - m(4) * (g(1) * (t99 * rSges(4,1) + t88 * rSges(4,2) - t47 * pkin(2) + t126 * t46 + t102) + g(2) * (rSges(4,1) * t23 + rSges(4,2) * t22 + pkin(2) * t49 - t126 * t48 + t107)) - m(5) * (g(1) * (t100 * rSges(5,1) + t89 * rSges(5,2) + pkin(3) * t59 + t106 * t46 - t47 * t66 + t102) + g(2) * (rSges(5,1) * t21 + rSges(5,2) * t20 + pkin(3) * t104 - t106 * t48 + t49 * t66 + t107)) - m(6) * (g(1) * (-rSges(6,1) * t15 - rSges(6,2) * t101 - rSges(6,3) * t46 + t90) + g(2) * (rSges(6,1) * t19 - rSges(6,2) * t18 + rSges(6,3) * t48 + t97)) - m(7) * (g(1) * (rSges(7,1) * t137 + rSges(7,2) * t138 - t15 * pkin(5) + t125 * t101 + t90) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t19 + t125 * t18 + t97)) -m(3) * (g(1) * (-t48 * rSges(3,1) - t49 * rSges(3,2)) + g(2) * (-t46 * rSges(3,1) - t47 * rSges(3,2)) + (rSges(3,1) * t82 - rSges(3,2) * t78) * t127) - m(4) * (g(1) * (-t126 * t49 - t93 * t48) - t126 * t129 - t93 * t130 + (-t126 * t78 + t93 * t82) * t127) - m(5) * (g(1) * (-t106 * t49 - t91 * t48) - t106 * t129 - t91 * t130 + (-t106 * t78 + t91 * t82) * t127) - m(6) * (g(1) * (t49 * rSges(6,3) - t95 * t48 + t111) + g(2) * (t47 * rSges(6,3) - t95 * t46 + t112) + t128 + (t95 * t82 + (rSges(6,3) - t72) * t78) * t127) - m(7) * (g(1) * (t135 * t49 + t111) + g(2) * (t135 * t47 + t112) + t128 + ((-t72 + t135) * t78 + t133 * t82) * t127 - t134 * t133) -m(4) * (g(1) * (rSges(4,1) * t22 - rSges(4,2) * t23) + g(2) * (-t88 * rSges(4,1) + t99 * rSges(4,2)) + g(3) * (t87 * rSges(4,1) + (-t105 * t77 - t78 * t116) * rSges(4,2))) - m(5) * (g(1) * (rSges(5,1) * t20 - rSges(5,2) * t21) + g(2) * (-t89 * rSges(5,1) + t100 * rSges(5,2)) + g(3) * ((t105 * t68 - t67 * t118) * rSges(5,1) + (-t105 * t67 - t68 * t118) * rSges(5,2)) + (g(1) * t22 - g(2) * t88 + g(3) * t87) * pkin(3)) - m(6) * (g(1) * (t109 + t124) + g(2) * (t94 + t132) + g(3) * (t108 + t110)) - m(7) * (g(1) * (t85 + t109) + g(2) * (t86 + t94) + g(3) * (t84 + t108)) (-m(5) - m(6) - m(7)) * (t134 - t136) -m(6) * (g(1) * t124 + g(2) * t132 + g(3) * t110) - m(7) * (g(1) * t85 + g(2) * t86 + g(3) * t84) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t138 * rSges(7,1) + rSges(7,2) * t137) + g(3) * ((-t80 * t115 - t36 * t76) * rSges(7,1) + (t76 * t115 - t36 * t80) * rSges(7,2)))];
taug  = t3(:);
