% Calculate Gravitation load on the joints for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR15_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:35
% EndTime: 2019-03-09 20:24:41
% DurationCPUTime: 1.98s
% Computational Cost: add. (1220->239), mult. (3269->348), div. (0->0), fcn. (4125->14), ass. (0->94)
t118 = cos(pkin(7));
t81 = sin(qJ(3));
t111 = t81 * t118;
t78 = sin(pkin(6));
t86 = cos(qJ(1));
t123 = t78 * t86;
t131 = cos(qJ(3));
t119 = cos(pkin(6));
t110 = t86 * t119;
t130 = sin(qJ(1));
t82 = sin(qJ(2));
t85 = cos(qJ(2));
t58 = -t110 * t85 + t130 * t82;
t59 = t110 * t82 + t130 * t85;
t77 = sin(pkin(7));
t23 = -t77 * t81 * t123 - t111 * t58 + t131 * t59;
t112 = t78 * t118;
t149 = t86 * t112 - t58 * t77;
t103 = t118 * t131;
t115 = t78 * t131;
t22 = t115 * t77 * t86 + t103 * t58 + t59 * t81;
t80 = sin(qJ(5));
t84 = cos(qJ(5));
t7 = -t149 * t84 + t22 * t80;
t79 = sin(qJ(6));
t83 = cos(qJ(6));
t153 = -t23 * t83 + t7 * t79;
t152 = -t23 * t79 - t7 * t83;
t104 = t119 * t130;
t60 = -t104 * t85 - t86 * t82;
t43 = t130 * t112 - t60 * t77;
t148 = t149 * t80 + t22 * t84;
t133 = rSges(6,3) + pkin(11);
t132 = rSges(7,3) + pkin(12);
t61 = -t104 * t82 + t85 * t86;
t147 = g(1) * t61 + g(2) * t59;
t114 = t78 * t130;
t109 = t77 * t114;
t27 = t61 * t131 + (t118 * t60 + t109) * t81;
t113 = t77 * t119;
t38 = t81 * t113 + (t111 * t85 + t131 * t82) * t78;
t146 = g(1) * t27 + g(2) * t23 + g(3) * t38;
t26 = -t103 * t60 - t109 * t131 + t61 * t81;
t124 = t78 * t85;
t125 = t78 * t82;
t37 = -t103 * t124 - t113 * t131 + t125 * t81;
t145 = g(1) * t26 + g(2) * t22 + g(3) * t37;
t127 = t77 * t80;
t126 = t77 * t84;
t116 = t77 * t125;
t122 = pkin(2) * t124 + pkin(10) * t116;
t121 = t86 * pkin(1) + pkin(9) * t114;
t120 = rSges(5,3) + qJ(4);
t117 = g(3) * t125;
t30 = t103 * t59 - t58 * t81;
t31 = -t111 * t59 - t131 * t58;
t53 = t58 * pkin(2);
t107 = t31 * pkin(3) + t30 * qJ(4) - t53;
t32 = t103 * t61 + t60 * t81;
t33 = -t111 * t61 + t131 * t60;
t55 = t60 * pkin(2);
t106 = t33 * pkin(3) + t32 * qJ(4) + t55;
t105 = -pkin(1) * t130 + pkin(9) * t123;
t51 = (t103 * t82 + t81 * t85) * t78;
t52 = -t111 * t125 + t115 * t85;
t102 = t52 * pkin(3) + t51 * qJ(4) + t122;
t101 = t83 * rSges(7,1) - t79 * rSges(7,2) + pkin(5);
t97 = pkin(4) * t116 + t102;
t96 = t31 * pkin(11) + t107;
t95 = t33 * pkin(11) + t106;
t94 = t61 * pkin(2) + t43 * pkin(10) + t121;
t93 = t27 * pkin(3) + t94;
t92 = -t59 * pkin(2) + t149 * pkin(10) + t105;
t91 = t147 * (pkin(4) + pkin(10)) * t77;
t90 = -pkin(3) * t23 + t92;
t89 = t43 * pkin(4) + qJ(4) * t26 + t93;
t87 = pkin(4) * t149 - qJ(4) * t22 + t90;
t57 = t118 * t119 - t124 * t77;
t36 = t37 * pkin(3);
t35 = t116 * t84 + t51 * t80;
t34 = t116 * t80 - t51 * t84;
t21 = t37 * t80 + t57 * t84;
t20 = t37 * t84 - t57 * t80;
t18 = t26 * pkin(3);
t16 = t22 * pkin(3);
t13 = t126 * t61 + t32 * t80;
t12 = -t127 * t61 + t32 * t84;
t11 = t126 * t59 + t30 * t80;
t10 = -t127 * t59 + t30 * t84;
t9 = t26 * t80 + t43 * t84;
t8 = -t26 * t84 + t43 * t80;
t2 = t27 * t79 + t83 * t9;
t1 = t27 * t83 - t79 * t9;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t130 - t86 * rSges(2,2)) + g(2) * (t86 * rSges(2,1) - rSges(2,2) * t130)) - m(3) * (g(1) * (-t59 * rSges(3,1) + t58 * rSges(3,2) + rSges(3,3) * t123 + t105) + g(2) * (t61 * rSges(3,1) + t60 * rSges(3,2) + rSges(3,3) * t114 + t121)) - m(4) * (g(1) * (-rSges(4,1) * t23 + rSges(4,2) * t22 + rSges(4,3) * t149 + t92) + g(2) * (rSges(4,1) * t27 - rSges(4,2) * t26 + rSges(4,3) * t43 + t94)) - m(5) * (g(1) * (rSges(5,1) * t149 + rSges(5,2) * t23 - t120 * t22 + t90) + g(2) * (rSges(5,1) * t43 - rSges(5,2) * t27 + t120 * t26 + t93)) - m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t148 - t133 * t23 + t87) + g(2) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t133 * t27 + t89)) - m(7) * (g(1) * (t152 * rSges(7,1) + t153 * rSges(7,2) - t7 * pkin(5) - t23 * pkin(11) + t132 * t148 + t87) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t9 + pkin(11) * t27 + t132 * t8 + t89)) -m(3) * (g(1) * (rSges(3,1) * t60 - rSges(3,2) * t61) + g(2) * (-rSges(3,1) * t58 - rSges(3,2) * t59) + g(3) * (rSges(3,1) * t85 - rSges(3,2) * t82) * t78) - m(4) * (g(1) * (rSges(4,1) * t33 - rSges(4,2) * t32 + t55) + g(2) * (rSges(4,1) * t31 - rSges(4,2) * t30 - t53) + g(3) * (rSges(4,1) * t52 - rSges(4,2) * t51 + t122) + (rSges(4,3) * t117 + t147 * (rSges(4,3) + pkin(10))) * t77) - m(5) * (g(1) * (-rSges(5,2) * t33 + rSges(5,3) * t32 + t106) + g(2) * (-rSges(5,2) * t31 + rSges(5,3) * t30 + t107) + g(3) * (-rSges(5,2) * t52 + rSges(5,3) * t51 + t102) + (rSges(5,1) * t117 + t147 * (rSges(5,1) + pkin(10))) * t77) - m(6) * (g(1) * (rSges(6,1) * t13 + rSges(6,2) * t12 + rSges(6,3) * t33 + t95) + g(2) * (rSges(6,1) * t11 + rSges(6,2) * t10 + rSges(6,3) * t31 + t96) + g(3) * (rSges(6,1) * t35 - rSges(6,2) * t34 + t133 * t52 + t97) + t91) - m(7) * (g(1) * (t13 * pkin(5) + (t13 * t83 + t33 * t79) * rSges(7,1) + (-t13 * t79 + t33 * t83) * rSges(7,2) + t95 - t132 * t12) + g(2) * (t11 * pkin(5) + (t11 * t83 + t31 * t79) * rSges(7,1) + (-t11 * t79 + t31 * t83) * rSges(7,2) + t96 - t132 * t10) + g(3) * (t35 * pkin(5) + t52 * pkin(11) + (t35 * t83 + t52 * t79) * rSges(7,1) + (-t35 * t79 + t52 * t83) * rSges(7,2) + t132 * t34 + t97) + t91) -m(4) * (g(1) * (-rSges(4,1) * t26 - rSges(4,2) * t27) + g(2) * (-rSges(4,1) * t22 - rSges(4,2) * t23) + g(3) * (-rSges(4,1) * t37 - rSges(4,2) * t38)) - m(5) * (g(1) * (rSges(5,2) * t26 + t120 * t27 - t18) + g(2) * (rSges(5,2) * t22 + t120 * t23 - t16) + g(3) * (rSges(5,2) * t37 + t120 * t38 - t36)) + (-g(1) * (-t133 * t26 - t18) - g(2) * (-t133 * t22 - t16) - g(3) * (-t133 * t37 - t36) - t146 * (rSges(6,1) * t80 + rSges(6,2) * t84 + qJ(4))) * m(6) + (g(1) * t18 + g(2) * t16 + g(3) * t36 - t146 * (t101 * t80 - t132 * t84 + qJ(4)) - t145 * (-t79 * rSges(7,1) - t83 * rSges(7,2) - pkin(11))) * m(7) (-m(5) - m(6) - m(7)) * t145, -m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (rSges(6,1) * t148 - rSges(6,2) * t7) + g(3) * (rSges(6,1) * t20 - rSges(6,2) * t21)) - m(7) * (g(1) * (-t101 * t8 + t132 * t9) + (t101 * t20 + t132 * t21) * g(3) + (t101 * t148 + t132 * t7) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t153 * rSges(7,1) + t152 * rSges(7,2)) + g(3) * ((-t21 * t79 + t38 * t83) * rSges(7,1) + (-t21 * t83 - t38 * t79) * rSges(7,2)))];
taug  = t3(:);
