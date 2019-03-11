% Calculate Gravitation load on the joints for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:17
% EndTime: 2019-03-09 04:01:21
% DurationCPUTime: 1.50s
% Computational Cost: add. (1199->193), mult. (3175->296), div. (0->0), fcn. (4079->16), ass. (0->93)
t70 = sin(pkin(6));
t82 = cos(qJ(1));
t119 = t70 * t82;
t67 = sin(pkin(13));
t71 = cos(pkin(13));
t77 = sin(qJ(3));
t81 = cos(qJ(3));
t55 = t67 * t77 - t81 * t71;
t69 = sin(pkin(7));
t44 = t55 * t69;
t73 = cos(pkin(7));
t46 = t55 * t73;
t74 = cos(pkin(6));
t116 = t74 * t82;
t68 = sin(pkin(12));
t78 = sin(qJ(1));
t122 = t68 * t78;
t72 = cos(pkin(12));
t51 = -t72 * t116 + t122;
t114 = t78 * t72;
t52 = t68 * t116 + t114;
t97 = t67 * t81 + t71 * t77;
t18 = -t44 * t119 - t51 * t46 + t52 * t97;
t34 = t73 * t119 - t51 * t69;
t76 = sin(qJ(5));
t80 = cos(qJ(5));
t45 = t97 * t69;
t47 = t97 * t73;
t84 = t45 * t119 + t51 * t47 + t52 * t55;
t6 = t34 * t76 + t80 * t84;
t75 = sin(qJ(6));
t79 = cos(qJ(6));
t140 = t18 * t79 + t6 * t75;
t139 = -t18 * t75 + t6 * t79;
t106 = t69 * t119;
t117 = t73 * t81;
t108 = pkin(3) * t117;
t125 = t52 * t77;
t83 = (-t81 * t106 - t125) * pkin(3) - t51 * t108;
t138 = -t18 * pkin(4) + t83;
t135 = t72 * t117 - t68 * t77;
t134 = -t34 * t80 + t76 * t84;
t127 = rSges(7,3) + pkin(11);
t133 = g(1) * t127;
t132 = (t51 * t73 + t106) * t81 + t125;
t27 = t45 * t74 + (t47 * t72 - t55 * t68) * t70;
t130 = pkin(3) * t81;
t129 = t80 * pkin(5);
t128 = -rSges(6,3) - pkin(10);
t54 = -t74 * t122 + t72 * t82;
t124 = t54 * t77;
t121 = t69 * t74;
t120 = t70 * t78;
t118 = t73 * t77;
t115 = t75 * t80;
t113 = t79 * t80;
t112 = pkin(9) + qJ(4);
t110 = qJ(2) * t70;
t111 = t82 * pkin(1) + t78 * t110;
t109 = -m(5) - m(6) - m(7);
t107 = t69 * t120;
t104 = g(2) * t127;
t103 = g(3) * t127;
t102 = -t78 * pkin(1) + t82 * t110;
t48 = pkin(3) * t69 * t77 + t112 * t73;
t49 = pkin(3) * t118 - t112 * t69;
t53 = -t74 * t114 - t68 * t82;
t65 = pkin(2) + t130;
t100 = t48 * t120 + t53 * t49 + t54 * t65 + t111;
t99 = rSges(6,1) * t80 - rSges(6,2) * t76;
t96 = -pkin(3) * t124 + t107 * t130 + t53 * t108;
t22 = t45 * t120 + t47 * t53 - t54 * t55;
t95 = t22 * pkin(4) + t100;
t92 = t73 * t120 - t53 * t69;
t91 = t53 * t73 + t107;
t90 = t135 * pkin(3) * t70 + t121 * t130;
t21 = -t44 * t120 - t46 * t53 - t54 * t97;
t89 = t21 * pkin(4) + t96;
t88 = t48 * t119 + t51 * t49 - t52 * t65 + t102;
t87 = t77 * t106 + t51 * t118 - t52 * t81;
t26 = -t44 * t74 + (-t46 * t72 - t68 * t97) * t70;
t86 = t26 * pkin(4) + t90;
t85 = pkin(4) * t84 + t88;
t50 = -t69 * t70 * t72 + t73 * t74;
t29 = t54 * t81 + t91 * t77;
t28 = t91 * t81 - t124;
t14 = t27 * t80 + t50 * t76;
t13 = -t27 * t76 + t50 * t80;
t8 = t22 * t80 + t92 * t76;
t7 = t22 * t76 - t92 * t80;
t2 = -t21 * t75 + t79 * t8;
t1 = -t21 * t79 - t75 * t8;
t3 = [-m(2) * (g(1) * (-t78 * rSges(2,1) - rSges(2,2) * t82) + g(2) * (rSges(2,1) * t82 - t78 * rSges(2,2))) - m(3) * (g(1) * (-t52 * rSges(3,1) + t51 * rSges(3,2) + rSges(3,3) * t119 + t102) + g(2) * (rSges(3,1) * t54 + rSges(3,2) * t53 + rSges(3,3) * t120 + t111)) - m(4) * (g(1) * (t87 * rSges(4,1) + t132 * rSges(4,2) + t34 * rSges(4,3) - t52 * pkin(2) + t102) + g(2) * (t29 * rSges(4,1) + t28 * rSges(4,2) + t92 * rSges(4,3) + t54 * pkin(2) + t111) + (g(1) * t34 + g(2) * t92) * pkin(9)) - m(5) * (g(1) * (rSges(5,1) * t84 + rSges(5,2) * t18 + rSges(5,3) * t34 + t88) + g(2) * (t22 * rSges(5,1) + t21 * rSges(5,2) + t92 * rSges(5,3) + t100)) - m(6) * (g(1) * (rSges(6,1) * t6 - rSges(6,2) * t134 + t128 * t18 + t85) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t128 * t21 + t95)) - m(7) * (g(1) * (t139 * rSges(7,1) - t140 * rSges(7,2) + t6 * pkin(5) - t18 * pkin(10) + t127 * t134 + t85) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t8 - pkin(10) * t21 + t127 * t7 + t95)) (-m(3) - m(4) + t109) * (g(3) * t74 + (g(1) * t78 - g(2) * t82) * t70) -m(4) * (g(1) * (rSges(4,1) * t28 - rSges(4,2) * t29) + g(2) * (-t132 * rSges(4,1) + t87 * rSges(4,2)) + g(3) * ((t81 * rSges(4,1) - t77 * rSges(4,2)) * t121 + (t135 * rSges(4,1) + (-t72 * t118 - t68 * t81) * rSges(4,2)) * t70)) - m(5) * (g(1) * (rSges(5,1) * t21 - rSges(5,2) * t22 + t96) + g(2) * (-rSges(5,1) * t18 + rSges(5,2) * t84 + t83) + g(3) * (rSges(5,1) * t26 - rSges(5,2) * t27 + t90)) - m(6) * (g(1) * (-t128 * t22 + t99 * t21 + t89) + g(2) * (t128 * t84 - t18 * t99 + t138) + g(3) * (-t128 * t27 + t99 * t26 + t86)) - m(7) * (g(1) * (t21 * t129 + t22 * pkin(10) + (t21 * t113 + t22 * t75) * rSges(7,1) + (-t21 * t115 + t22 * t79) * rSges(7,2) + t89) + g(2) * (-t18 * t129 - t84 * pkin(10) + (-t113 * t18 - t75 * t84) * rSges(7,1) + (t115 * t18 - t79 * t84) * rSges(7,2) + t138) + g(3) * (t26 * t129 + t27 * pkin(10) + (t26 * t113 + t27 * t75) * rSges(7,1) + (-t26 * t115 + t27 * t79) * rSges(7,2) + t86) + (t26 * t103 - t104 * t18 + t21 * t133) * t76) t109 * (g(1) * t92 - g(2) * t34 + g(3) * t50) -m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (rSges(6,1) * t134 + rSges(6,2) * t6) + g(3) * (rSges(6,1) * t13 - rSges(6,2) * t14)) - m(7) * (t8 * t133 + t14 * t103 - t6 * t104 + (-g(1) * t7 + g(2) * t134 + g(3) * t13) * (t79 * rSges(7,1) - t75 * rSges(7,2) + pkin(5))) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (rSges(7,1) * t140 + t139 * rSges(7,2)) + g(3) * ((-t14 * t75 - t26 * t79) * rSges(7,1) + (-t14 * t79 + t26 * t75) * rSges(7,2)))];
taug  = t3(:);
