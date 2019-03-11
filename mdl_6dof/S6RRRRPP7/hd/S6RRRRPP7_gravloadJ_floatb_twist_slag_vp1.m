% Calculate Gravitation load on the joints for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:10
% EndTime: 2019-03-09 21:17:12
% DurationCPUTime: 1.67s
% Computational Cost: add. (882->212), mult. (1813->298), div. (0->0), fcn. (2191->12), ass. (0->93)
t135 = cos(qJ(1));
t75 = sin(pkin(6));
t106 = t75 * t135;
t111 = cos(pkin(6));
t100 = t111 * t135;
t134 = sin(qJ(1));
t79 = sin(qJ(2));
t82 = cos(qJ(2));
t54 = t100 * t79 + t134 * t82;
t78 = sin(qJ(3));
t81 = cos(qJ(3));
t25 = -t106 * t78 + t54 * t81;
t53 = -t100 * t82 + t134 * t79;
t77 = sin(qJ(4));
t80 = cos(qJ(4));
t153 = -t25 * t77 + t53 * t80;
t145 = pkin(4) * t77;
t154 = pkin(9) + t145;
t105 = t75 * t134;
t99 = t111 * t134;
t56 = t135 * t82 - t79 * t99;
t29 = t105 * t78 + t56 * t81;
t55 = t135 * t79 + t82 * t99;
t8 = -t29 * t77 + t55 * t80;
t24 = t106 * t81 + t54 * t78;
t28 = -t105 * t81 + t56 * t78;
t122 = t75 * t79;
t51 = -t111 * t81 + t122 * t78;
t152 = -g(1) * t28 - g(2) * t24 - g(3) * t51;
t129 = t53 * t77;
t151 = -t25 * t80 - t129;
t74 = qJ(4) + pkin(11);
t71 = sin(t74);
t72 = cos(t74);
t2 = t25 * t71 - t53 * t72;
t3 = t25 * t72 + t53 * t71;
t138 = rSges(7,1) + pkin(5);
t112 = rSges(7,3) + qJ(6);
t149 = g(1) * t55 + g(2) * t53;
t121 = t75 * t82;
t147 = (g(3) * t121 - t149) * t78;
t136 = pkin(10) + rSges(5,3);
t91 = rSges(5,1) * t80 - rSges(5,2) * t77 + pkin(3);
t146 = t136 * t78 + t91 * t81;
t139 = g(3) * t75;
t137 = rSges(4,3) + pkin(9);
t127 = t55 * t77;
t70 = pkin(4) * t80 + pkin(3);
t125 = t70 * t81;
t124 = t71 * t81;
t123 = t72 * t81;
t120 = t81 * t82;
t76 = -qJ(5) - pkin(10);
t117 = -t24 * t70 - t25 * t76;
t116 = -t28 * t70 - t29 * t76;
t52 = t111 * t78 + t122 * t81;
t115 = -t51 * t70 - t52 * t76;
t114 = pkin(1) * t135 + pkin(8) * t105;
t113 = pkin(2) * t121 + pkin(9) * t122;
t109 = t71 * t121;
t108 = pkin(2) * t56 + t114;
t107 = g(3) * t113;
t104 = t120 * t70 * t75 + t122 * t145 + t113;
t103 = t153 * pkin(4);
t102 = t8 * pkin(4);
t101 = -pkin(1) * t134 + pkin(8) * t106;
t98 = rSges(4,1) * t81 - rSges(4,2) * t78;
t97 = rSges(5,1) * t77 + rSges(5,2) * t80;
t96 = -rSges(6,1) * t72 + rSges(6,2) * t71;
t47 = t53 * pkin(2);
t95 = -t53 * t125 + t154 * t54 - t47;
t49 = t55 * pkin(2);
t94 = -t55 * t125 + t154 * t56 - t49;
t93 = pkin(9) * t55 + t108;
t92 = -pkin(2) * t54 + t101;
t90 = pkin(9) + t97;
t89 = -t121 * t80 - t52 * t77;
t88 = t89 * pkin(4);
t87 = -t53 * pkin(9) + t92;
t86 = pkin(4) * t127 - t28 * t76 + t29 * t70 + t93;
t84 = -pkin(4) * t129 + t24 * t76 - t25 * t70 + t87;
t31 = (t120 * t72 + t71 * t79) * t75;
t30 = t109 * t81 - t122 * t72;
t19 = t52 * t72 - t109;
t18 = t121 * t72 + t52 * t71;
t13 = -t123 * t55 + t56 * t71;
t12 = -t124 * t55 - t56 * t72;
t11 = -t123 * t53 + t54 * t71;
t10 = -t124 * t53 - t54 * t72;
t9 = t29 * t80 + t127;
t7 = t29 * t72 + t55 * t71;
t6 = t29 * t71 - t55 * t72;
t1 = [-m(2) * (g(1) * (-t134 * rSges(2,1) - t135 * rSges(2,2)) + g(2) * (t135 * rSges(2,1) - t134 * rSges(2,2))) - m(3) * (g(1) * (-t54 * rSges(3,1) + t53 * rSges(3,2) + rSges(3,3) * t106 + t101) + g(2) * (t56 * rSges(3,1) - t55 * rSges(3,2) + rSges(3,3) * t105 + t114)) - m(4) * (g(1) * (-rSges(4,1) * t25 + rSges(4,2) * t24 - t137 * t53 + t92) + g(2) * (rSges(4,1) * t29 - rSges(4,2) * t28 + t137 * t55 + t108)) - m(5) * (g(1) * (t151 * rSges(5,1) - rSges(5,2) * t153 - t25 * pkin(3) - t136 * t24 + t87) + g(2) * (rSges(5,1) * t9 + rSges(5,2) * t8 + pkin(3) * t29 + t136 * t28 + t93)) - m(6) * (g(1) * (-rSges(6,1) * t3 + rSges(6,2) * t2 - rSges(6,3) * t24 + t84) + g(2) * (rSges(6,1) * t7 - rSges(6,2) * t6 + rSges(6,3) * t28 + t86)) - m(7) * (g(1) * (-rSges(7,2) * t24 - t112 * t2 - t138 * t3 + t84) + g(2) * (rSges(7,2) * t28 + t112 * t6 + t138 * t7 + t86)) -m(3) * (g(1) * (-rSges(3,1) * t55 - rSges(3,2) * t56) + g(2) * (-rSges(3,1) * t53 - rSges(3,2) * t54) + (rSges(3,1) * t82 - rSges(3,2) * t79) * t139) - m(4) * (g(1) * (t137 * t56 - t55 * t98 - t49) + g(2) * (t137 * t54 - t53 * t98 - t47) + t107 + (rSges(4,3) * t79 + t82 * t98) * t139) - m(5) * (t107 + (t146 * t82 + t97 * t79) * t139 - t149 * t146 + (t54 * t90 - t47) * g(2) + (t56 * t90 - t49) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t13 - rSges(6,2) * t12 + t94) + g(2) * (rSges(6,1) * t11 - rSges(6,2) * t10 + t95) + g(3) * (rSges(6,1) * t31 - rSges(6,2) * t30 + t104) + (rSges(6,3) - t76) * t147) - m(7) * (g(1) * (t112 * t12 + t13 * t138 + t94) + g(2) * (t112 * t10 + t138 * t11 + t95) + g(3) * (t112 * t30 + t138 * t31 + t104) + (rSges(7,2) - t76) * t147) -m(4) * (g(1) * (-rSges(4,1) * t28 - rSges(4,2) * t29) + g(2) * (-rSges(4,1) * t24 - rSges(4,2) * t25) + g(3) * (-rSges(4,1) * t51 - rSges(4,2) * t52)) - m(5) * (t152 * t91 + (g(1) * t29 + g(2) * t25 + g(3) * t52) * t136) - m(6) * (g(1) * (rSges(6,3) * t29 + t96 * t28 + t116) + g(2) * (rSges(6,3) * t25 + t96 * t24 + t117) + g(3) * (rSges(6,3) * t52 + t96 * t51 + t115)) + (-g(1) * (rSges(7,2) * t29 + t116) - g(2) * (rSges(7,2) * t25 + t117) - g(3) * (rSges(7,2) * t52 + t115) + t152 * (-t112 * t71 - t138 * t72)) * m(7), -m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (rSges(5,1) * t153 + rSges(5,2) * t151) + g(3) * (t89 * rSges(5,1) + (t121 * t77 - t52 * t80) * rSges(5,2))) - m(6) * (g(1) * (-rSges(6,1) * t6 - rSges(6,2) * t7 + t102) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3 + t103) + g(3) * (-t18 * rSges(6,1) - t19 * rSges(6,2) + t88)) - m(7) * (g(1) * (t112 * t7 - t138 * t6 + t102) + g(2) * (t112 * t3 - t138 * t2 + t103) + g(3) * (t112 * t19 - t138 * t18 + t88)) -(-m(6) - m(7)) * t152, -m(7) * (g(1) * t6 + g(2) * t2 + g(3) * t18)];
taug  = t1(:);
