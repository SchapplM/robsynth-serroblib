% Calculate Gravitation load on the joints for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:30
% EndTime: 2019-03-09 17:29:34
% DurationCPUTime: 1.65s
% Computational Cost: add. (838->207), mult. (1698->296), div. (0->0), fcn. (2050->12), ass. (0->85)
t71 = sin(pkin(11));
t134 = pkin(4) * t71;
t138 = pkin(9) + t134;
t125 = cos(qJ(1));
t72 = sin(pkin(6));
t103 = t72 * t125;
t123 = sin(qJ(1));
t124 = cos(qJ(2));
t76 = sin(qJ(2));
t105 = cos(pkin(6));
t92 = t105 * t125;
t50 = t123 * t124 + t76 * t92;
t75 = sin(qJ(3));
t77 = cos(qJ(3));
t22 = t77 * t103 + t50 * t75;
t101 = t72 * t123;
t91 = t105 * t123;
t52 = t125 * t124 - t76 * t91;
t26 = -t77 * t101 + t52 * t75;
t115 = t72 * t76;
t47 = -t105 * t77 + t75 * t115;
t135 = g(1) * t26 + g(2) * t22 + g(3) * t47;
t23 = -t75 * t103 + t50 * t77;
t49 = t123 * t76 - t124 * t92;
t70 = pkin(11) + qJ(5);
t67 = sin(t70);
t68 = cos(t70);
t2 = t23 * t67 - t49 * t68;
t3 = t23 * t68 + t49 * t67;
t127 = rSges(7,1) + pkin(5);
t107 = rSges(7,3) + qJ(6);
t106 = qJ(4) + rSges(5,3);
t51 = t124 * t91 + t125 * t76;
t136 = g(1) * t51 + g(2) * t49;
t128 = g(3) * t72;
t126 = rSges(4,3) + pkin(9);
t120 = t49 * t71;
t119 = t51 * t71;
t73 = cos(pkin(11));
t66 = pkin(4) * t73 + pkin(3);
t118 = t66 * t77;
t117 = t67 * t77;
t116 = t68 * t77;
t74 = -pkin(10) - qJ(4);
t112 = -t22 * t66 - t23 * t74;
t27 = t75 * t101 + t52 * t77;
t111 = -t26 * t66 - t27 * t74;
t48 = t105 * t75 + t77 * t115;
t110 = -t47 * t66 - t48 * t74;
t109 = t125 * pkin(1) + pkin(8) * t101;
t102 = t72 * t124;
t108 = pkin(2) * t102 + pkin(9) * t115;
t104 = t52 * pkin(2) + t109;
t100 = t75 * t124;
t99 = t77 * t124;
t98 = t124 * t74;
t95 = t72 * t99;
t97 = t115 * t134 + t66 * t95 + t108;
t94 = t67 * t102;
t93 = -t123 * pkin(1) + pkin(8) * t103;
t90 = -rSges(4,1) * t77 + rSges(4,2) * t75;
t89 = -rSges(6,1) * t68 + rSges(6,2) * t67;
t43 = t49 * pkin(2);
t88 = -t49 * t118 + t138 * t50 - t43;
t45 = t51 * pkin(2);
t87 = -t51 * t118 + t138 * t52 - t45;
t86 = t51 * pkin(9) + t104;
t85 = -t50 * pkin(2) + t93;
t84 = -t73 * rSges(5,1) + t71 * rSges(5,2) - pkin(3);
t83 = t71 * rSges(5,1) + t73 * rSges(5,2) + pkin(9);
t82 = -t49 * pkin(9) + t85;
t81 = pkin(4) * t119 - t26 * t74 + t27 * t66 + t86;
t79 = -pkin(4) * t120 + t22 * t74 - t23 * t66 + t82;
t78 = -t106 * t75 + t84 * t77;
t29 = (t67 * t76 + t68 * t99) * t72;
t28 = -t68 * t115 + t77 * t94;
t17 = t48 * t68 - t94;
t16 = t68 * t102 + t48 * t67;
t11 = -t51 * t116 + t52 * t67;
t10 = -t51 * t117 - t52 * t68;
t9 = -t49 * t116 + t50 * t67;
t8 = -t49 * t117 - t50 * t68;
t7 = t27 * t68 + t51 * t67;
t6 = t27 * t67 - t51 * t68;
t1 = [-m(2) * (g(1) * (-t123 * rSges(2,1) - t125 * rSges(2,2)) + g(2) * (t125 * rSges(2,1) - t123 * rSges(2,2))) - m(3) * (g(1) * (-t50 * rSges(3,1) + t49 * rSges(3,2) + rSges(3,3) * t103 + t93) + g(2) * (t52 * rSges(3,1) - t51 * rSges(3,2) + rSges(3,3) * t101 + t109)) - m(4) * (g(1) * (-rSges(4,1) * t23 + rSges(4,2) * t22 - t126 * t49 + t85) + g(2) * (rSges(4,1) * t27 - rSges(4,2) * t26 + t126 * t51 + t104)) - m(5) * (g(1) * (-t23 * pkin(3) + (-t23 * t73 - t120) * rSges(5,1) + (t23 * t71 - t49 * t73) * rSges(5,2) - t106 * t22 + t82) + g(2) * (t27 * pkin(3) + (t27 * t73 + t119) * rSges(5,1) + (-t27 * t71 + t51 * t73) * rSges(5,2) + t106 * t26 + t86)) - m(6) * (g(1) * (-rSges(6,1) * t3 + rSges(6,2) * t2 - rSges(6,3) * t22 + t79) + g(2) * (rSges(6,1) * t7 - rSges(6,2) * t6 + rSges(6,3) * t26 + t81)) - m(7) * (g(1) * (-rSges(7,2) * t22 - t107 * t2 - t127 * t3 + t79) + g(2) * (rSges(7,2) * t26 + t107 * t6 + t127 * t7 + t81)) -m(3) * (g(1) * (-rSges(3,1) * t51 - rSges(3,2) * t52) + g(2) * (-rSges(3,1) * t49 - rSges(3,2) * t50) + (t124 * rSges(3,1) - rSges(3,2) * t76) * t128) - m(4) * (g(1) * (t126 * t52 + t90 * t51 - t45) + g(2) * (t126 * t50 + t90 * t49 - t43) + g(3) * t108 + (rSges(4,1) * t99 - rSges(4,2) * t100 + rSges(4,3) * t76) * t128) - m(5) * (g(1) * (t78 * t51 + t83 * t52 - t45) + g(2) * (t78 * t49 + t83 * t50 - t43) + g(3) * (pkin(3) * t95 + t108 + (t106 * t100 + (t71 * t76 + t73 * t99) * rSges(5,1) + (-t71 * t99 + t73 * t76) * rSges(5,2)) * t72)) - m(6) * (g(1) * (rSges(6,1) * t11 - rSges(6,2) * t10 + t87) + g(2) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t88) + g(3) * (t29 * rSges(6,1) - t28 * rSges(6,2) + t97) + ((t124 * rSges(6,3) - t98) * t128 + t136 * (-rSges(6,3) + t74)) * t75) - m(7) * (g(1) * (t107 * t10 + t127 * t11 + t87) + g(2) * (t107 * t8 + t127 * t9 + t88) + g(3) * (t107 * t28 + t127 * t29 + t97) + ((t124 * rSges(7,2) - t98) * t128 + t136 * (-rSges(7,2) + t74)) * t75) -m(4) * (g(1) * (-rSges(4,1) * t26 - rSges(4,2) * t27) + g(2) * (-rSges(4,1) * t22 - rSges(4,2) * t23) + g(3) * (-rSges(4,1) * t47 - rSges(4,2) * t48)) - m(5) * (t135 * t84 + (g(1) * t27 + g(2) * t23 + g(3) * t48) * t106) - m(6) * (g(1) * (rSges(6,3) * t27 + t89 * t26 + t111) + g(2) * (rSges(6,3) * t23 + t89 * t22 + t112) + g(3) * (rSges(6,3) * t48 + t89 * t47 + t110)) + (-g(1) * (rSges(7,2) * t27 + t111) - g(2) * (rSges(7,2) * t23 + t112) - g(3) * (rSges(7,2) * t48 + t110) - t135 * (-t107 * t67 - t127 * t68)) * m(7) (-m(5) - m(6) - m(7)) * t135, -m(6) * (g(1) * (-rSges(6,1) * t6 - rSges(6,2) * t7) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3) + g(3) * (-rSges(6,1) * t16 - rSges(6,2) * t17)) - m(7) * (g(1) * (t107 * t7 - t127 * t6) + g(2) * (t107 * t3 - t127 * t2) + g(3) * (t107 * t17 - t127 * t16)) -m(7) * (g(1) * t6 + g(2) * t2 + g(3) * t16)];
taug  = t1(:);
