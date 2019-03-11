% Calculate Gravitation load on the joints for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:08
% EndTime: 2019-03-09 23:51:12
% DurationCPUTime: 1.92s
% Computational Cost: add. (959->244), mult. (2388->345), div. (0->0), fcn. (2970->12), ass. (0->96)
t76 = sin(qJ(6));
t126 = sin(pkin(6));
t139 = cos(qJ(1));
t103 = t139 * t126;
t138 = cos(qJ(3));
t127 = cos(pkin(6));
t105 = t127 * t139;
t137 = sin(qJ(1));
t79 = sin(qJ(2));
t82 = cos(qJ(2));
t59 = t79 * t105 + t137 * t82;
t78 = sin(qJ(3));
t32 = -t78 * t103 + t59 * t138;
t58 = -t82 * t105 + t137 * t79;
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t8 = t32 * t77 - t58 * t81;
t80 = cos(qJ(6));
t9 = t32 * t81 + t58 * t77;
t150 = t76 * t9 - t8 * t80;
t149 = t76 * t8 + t80 * t9;
t148 = pkin(4) * t81 + qJ(5) * t77;
t140 = rSges(7,3) + pkin(11);
t115 = t82 * t126;
t102 = t138 * t126;
t57 = t79 * t102 + t127 * t78;
t29 = t81 * t115 + t57 * t77;
t30 = -t77 * t115 + t57 * t81;
t147 = (t29 * t80 - t30 * t76) * rSges(7,1) - (t29 * t76 + t30 * t80) * rSges(7,2);
t146 = -rSges(7,1) * t150 - t149 * rSges(7,2);
t83 = (t76 * t77 + t80 * t81) * rSges(7,1) - (t76 * t81 - t77 * t80) * rSges(7,2) + t81 * pkin(5);
t143 = rSges(6,2) + pkin(10);
t142 = rSges(4,3) + pkin(9);
t141 = rSges(5,3) + pkin(10);
t136 = rSges(4,2) * t78;
t134 = t58 * t78;
t104 = t127 * t137;
t60 = t82 * t104 + t139 * t79;
t132 = t60 * t78;
t101 = t126 * t137;
t131 = t139 * pkin(1) + pkin(8) * t101;
t116 = t79 * t126;
t130 = pkin(2) * t115 + pkin(9) * t116;
t128 = rSges(6,3) + qJ(5);
t125 = pkin(10) - t140;
t117 = -t139 * t102 - t59 * t78;
t25 = t117 * pkin(3);
t124 = t148 * t117 + t25;
t61 = -t79 * t104 + t139 * t82;
t35 = -t138 * t101 + t61 * t78;
t27 = t35 * pkin(3);
t123 = -t148 * t35 - t27;
t56 = -t78 * t116 + t127 * t138;
t51 = t56 * pkin(3);
t122 = t148 * t56 + t51;
t121 = t61 * pkin(2) + t131;
t120 = pkin(3) * t138;
t119 = t77 * t138;
t118 = t81 * t138;
t111 = t78 * t115;
t91 = t82 * t102;
t114 = pkin(3) * t91 + pkin(10) * t111 + t130;
t113 = g(1) * t125;
t112 = g(2) * t125;
t110 = -t137 * pkin(1) + pkin(8) * t103;
t39 = t77 * t116 + t81 * t91;
t109 = t39 * pkin(4) + t114;
t36 = t78 * t101 + t61 * t138;
t12 = t36 * t77 - t60 * t81;
t13 = t36 * t81 + t60 * t77;
t2 = t12 * t80 - t13 * t76;
t3 = t12 * t76 + t13 * t80;
t108 = rSges(7,1) * t2 - rSges(7,2) * t3;
t100 = rSges(5,1) * t81 - rSges(5,2) * t77;
t99 = rSges(6,1) * t81 + rSges(6,3) * t77;
t52 = t58 * pkin(2);
t94 = t59 * pkin(9) - pkin(10) * t134 - t58 * t120 - t52;
t54 = t60 * pkin(2);
t93 = t61 * pkin(9) - pkin(10) * t132 - t60 * t120 - t54;
t92 = -t59 * pkin(2) + t110;
t17 = -t58 * t118 + t59 * t77;
t90 = t17 * pkin(4) + t94;
t19 = -t60 * t118 + t61 * t77;
t89 = t19 * pkin(4) + t93;
t88 = t36 * pkin(3) + pkin(9) * t60 + t121;
t87 = -t138 * rSges(4,1) + t136;
t86 = t13 * pkin(4) + t88;
t85 = -pkin(3) * t32 - t58 * pkin(9) + t92;
t84 = -pkin(4) * t9 + t85;
t38 = -t81 * t116 + t77 * t91;
t24 = t29 * pkin(4);
t18 = -t60 * t119 - t61 * t81;
t16 = -t58 * t119 - t59 * t81;
t6 = t12 * pkin(4);
t4 = t8 * pkin(4);
t1 = [-m(2) * (g(1) * (-t137 * rSges(2,1) - t139 * rSges(2,2)) + g(2) * (t139 * rSges(2,1) - t137 * rSges(2,2))) - m(3) * (g(1) * (-t59 * rSges(3,1) + t58 * rSges(3,2) + rSges(3,3) * t103 + t110) + g(2) * (t61 * rSges(3,1) - t60 * rSges(3,2) + rSges(3,3) * t101 + t131)) - m(4) * (g(1) * (-rSges(4,1) * t32 - rSges(4,2) * t117 - t142 * t58 + t92) + g(2) * (rSges(4,1) * t36 - rSges(4,2) * t35 + t142 * t60 + t121)) - m(5) * (g(1) * (-rSges(5,1) * t9 + rSges(5,2) * t8 + t117 * t141 + t85) + g(2) * (rSges(5,1) * t13 - rSges(5,2) * t12 + t141 * t35 + t88)) - m(6) * (g(1) * (-rSges(6,1) * t9 + t117 * t143 - t128 * t8 + t84) + g(2) * (rSges(6,1) * t13 + t128 * t12 + t143 * t35 + t86)) - m(7) * (g(1) * (-t149 * rSges(7,1) + rSges(7,2) * t150 - t9 * pkin(5) - t8 * qJ(5) + t84) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t13 + qJ(5) * t12 + t86) + t35 * t112 + t117 * t113) -m(3) * (g(1) * (-rSges(3,1) * t60 - rSges(3,2) * t61) + g(2) * (-rSges(3,1) * t58 - rSges(3,2) * t59) + g(3) * (rSges(3,1) * t115 - rSges(3,2) * t116)) - m(4) * (g(1) * (t142 * t61 + t87 * t60 - t54) + g(2) * (t142 * t59 + t87 * t58 - t52) + g(3) * (rSges(4,3) * t116 + (rSges(4,1) * t102 - t126 * t136) * t82 + t130)) - m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t18 - rSges(5,3) * t132 + t93) + g(2) * (rSges(5,1) * t17 - rSges(5,2) * t16 - rSges(5,3) * t134 + t94) + g(3) * (t39 * rSges(5,1) - t38 * rSges(5,2) + rSges(5,3) * t111 + t114)) - m(6) * (g(1) * (rSges(6,1) * t19 - rSges(6,2) * t132 + t128 * t18 + t89) + g(2) * (rSges(6,1) * t17 - rSges(6,2) * t134 + t128 * t16 + t90) + g(3) * (t39 * rSges(6,1) + rSges(6,2) * t111 + t128 * t38 + t109)) - m(7) * (g(1) * (t19 * pkin(5) + t18 * qJ(5) + (t18 * t76 + t19 * t80) * rSges(7,1) + (t18 * t80 - t19 * t76) * rSges(7,2) + t89) + g(2) * (t17 * pkin(5) + t16 * qJ(5) + (t16 * t76 + t17 * t80) * rSges(7,1) + (t16 * t80 - t17 * t76) * rSges(7,2) + t90) + g(3) * (t39 * pkin(5) + t38 * qJ(5) + (t38 * t76 + t39 * t80) * rSges(7,1) + (t38 * t80 - t39 * t76) * rSges(7,2) + t109) + (g(1) * t60 + g(2) * t58 - g(3) * t115) * t140 * t78) -m(4) * (g(1) * (-rSges(4,1) * t35 - rSges(4,2) * t36) + g(2) * (rSges(4,1) * t117 - rSges(4,2) * t32)) - m(5) * (g(1) * (-t100 * t35 + t141 * t36 - t27) + g(2) * (t100 * t117 + t141 * t32 + t25)) - m(6) * (g(1) * (t143 * t36 - t99 * t35 + t123) + g(2) * (t117 * t99 + t143 * t32 + t124)) - m(7) * (t32 * t112 + t36 * t113 + (t117 * t83 + t124) * g(2) + (-t35 * t83 + t123) * g(1)) + (-m(4) * (rSges(4,1) * t56 - rSges(4,2) * t57) - m(5) * (t100 * t56 + t141 * t57 + t51) - m(6) * (t143 * t57 + t99 * t56 + t122) - m(7) * (t125 * t57 + t83 * t56 + t122)) * g(3), -m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,2) * t13) + g(2) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(3) * (-rSges(5,1) * t29 - rSges(5,2) * t30)) - m(6) * (g(1) * (-rSges(6,1) * t12 + t128 * t13 - t6) + g(2) * (-rSges(6,1) * t8 + t128 * t9 - t4) + g(3) * (-rSges(6,1) * t29 + t128 * t30 - t24)) - m(7) * (g(1) * (-t12 * pkin(5) + t13 * qJ(5) - t108 - t6) + g(2) * (-t8 * pkin(5) + t9 * qJ(5) - t146 - t4) + g(3) * (-t29 * pkin(5) + t30 * qJ(5) - t147 - t24)) (-m(6) - m(7)) * (g(1) * t12 + g(2) * t8 + g(3) * t29) -m(7) * (g(1) * t108 + g(2) * t146 + g(3) * t147)];
taug  = t1(:);
