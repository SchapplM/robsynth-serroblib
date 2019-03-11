% Calculate Gravitation load on the joints for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:14
% EndTime: 2019-03-08 23:48:18
% DurationCPUTime: 1.18s
% Computational Cost: add. (1053->189), mult. (2885->280), div. (0->0), fcn. (3658->14), ass. (0->94)
t65 = sin(qJ(6));
t68 = cos(qJ(6));
t135 = -t65 * rSges(7,1) - t68 * rSges(7,2);
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t132 = -pkin(4) * t69 - qJ(5) * t66;
t123 = pkin(11) + rSges(7,3);
t112 = cos(pkin(12));
t113 = cos(pkin(7));
t109 = sin(pkin(12));
t120 = sin(qJ(2));
t122 = cos(qJ(2));
t114 = cos(pkin(6));
t94 = t114 * t112;
t75 = t109 * t120 - t122 * t94;
t110 = sin(pkin(7));
t111 = sin(pkin(6));
t91 = t111 * t110;
t131 = t112 * t91 + t75 * t113;
t93 = t114 * t109;
t76 = t112 * t120 + t122 * t93;
t90 = t111 * t109;
t130 = -t110 * t90 + t76 * t113;
t92 = t113 * t111;
t129 = t110 * t114 + t122 * t92;
t128 = t135 * t66;
t87 = t68 * rSges(7,1) - t65 * rSges(7,2) + pkin(5) + pkin(10);
t127 = -m(6) - m(7);
t125 = rSges(6,1) + pkin(10);
t124 = rSges(5,3) + pkin(10);
t121 = cos(qJ(3));
t82 = t120 * t91;
t98 = t122 * t111;
t117 = pkin(2) * t98 + pkin(9) * t82;
t115 = rSges(6,3) + qJ(5);
t55 = t109 * t122 + t120 * t94;
t67 = sin(qJ(3));
t21 = t121 * t131 + t55 * t67;
t18 = t21 * pkin(3);
t108 = t132 * t21 - t18;
t56 = t112 * t122 - t120 * t93;
t23 = t121 * t130 + t56 * t67;
t19 = t23 * pkin(3);
t107 = t132 * t23 - t19;
t97 = t111 * t120;
t41 = -t121 * t129 + t67 * t97;
t40 = t41 * pkin(3);
t106 = t132 * t41 - t40;
t83 = t120 * t92;
t51 = t121 * t98 - t67 * t83;
t105 = t51 * pkin(3) + t117;
t104 = t110 * pkin(9);
t103 = t66 * t110;
t102 = t67 * t113;
t101 = t69 * t110;
t35 = t51 * t69 + t66 * t82;
t100 = t35 * pkin(4) + t105;
t99 = t113 * t121;
t96 = -rSges(5,1) * t69 + rSges(5,2) * t66;
t95 = rSges(6,2) * t69 - rSges(6,3) * t66;
t88 = qJ(5) - t135;
t30 = -t102 * t55 - t121 * t75;
t53 = t75 * pkin(2);
t86 = t30 * pkin(3) + t104 * t55 - t53;
t32 = -t102 * t56 - t121 * t76;
t54 = t76 * pkin(2);
t85 = t32 * pkin(3) + t104 * t56 - t54;
t11 = t103 * t55 + t30 * t69;
t81 = t11 * pkin(4) + t86;
t13 = t103 * t56 + t32 * t69;
t80 = t13 * pkin(4) + t85;
t77 = rSges(4,3) * t110 + t104;
t74 = t113 * t114 - t122 * t91;
t71 = t110 * t76 + t113 * t90;
t70 = t110 * t75 - t112 * t92;
t50 = t121 * t83 + t67 * t98;
t42 = t121 * t97 + t129 * t67;
t34 = t51 * t66 - t69 * t82;
t31 = t56 * t99 - t67 * t76;
t29 = t55 * t99 - t67 * t75;
t26 = t42 * t69 + t66 * t74;
t25 = t42 * t66 - t69 * t74;
t24 = t56 * t121 - t130 * t67;
t22 = t121 * t55 - t131 * t67;
t20 = t25 * pkin(4);
t12 = -t101 * t56 + t32 * t66;
t10 = -t101 * t55 + t30 * t66;
t7 = t24 * t69 + t66 * t71;
t6 = t24 * t66 - t69 * t71;
t5 = t22 * t69 + t66 * t70;
t4 = t22 * t66 - t69 * t70;
t3 = t6 * pkin(4);
t2 = t4 * pkin(4);
t1 = [(-m(2) - m(3) - m(4) - m(5) + t127) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t76 - t56 * rSges(3,2)) + g(2) * (-rSges(3,1) * t75 - t55 * rSges(3,2)) + g(3) * (rSges(3,1) * t98 - rSges(3,2) * t97)) - m(4) * (g(1) * (t32 * rSges(4,1) - t31 * rSges(4,2) + t56 * t77 - t54) + g(2) * (t30 * rSges(4,1) - t29 * rSges(4,2) + t55 * t77 - t53) + g(3) * (t51 * rSges(4,1) - t50 * rSges(4,2) + rSges(4,3) * t82 + t117)) - m(5) * (g(1) * (t13 * rSges(5,1) - t12 * rSges(5,2) + t124 * t31 + t85) + g(2) * (t11 * rSges(5,1) - t10 * rSges(5,2) + t124 * t29 + t86) + g(3) * (rSges(5,1) * t35 - rSges(5,2) * t34 + t124 * t50 + t105)) - m(6) * (g(1) * (-t13 * rSges(6,2) + t115 * t12 + t125 * t31 + t80) + g(2) * (-t11 * rSges(6,2) + t10 * t115 + t125 * t29 + t81) + g(3) * (-rSges(6,2) * t35 + t115 * t34 + t125 * t50 + t100)) - m(7) * (g(1) * (t12 * t88 + t123 * t13 + t31 * t87 + t80) + g(2) * (t10 * t88 + t11 * t123 + t29 * t87 + t81) + g(3) * (t123 * t35 + t34 * t88 + t50 * t87 + t100)) -m(4) * (g(1) * (-rSges(4,1) * t23 - rSges(4,2) * t24) + g(2) * (-rSges(4,1) * t21 - rSges(4,2) * t22) + g(3) * (-rSges(4,1) * t41 - rSges(4,2) * t42)) - m(5) * (g(1) * (t124 * t24 + t23 * t96 - t19) + g(2) * (t124 * t22 + t21 * t96 - t18) + g(3) * (t124 * t42 + t41 * t96 - t40)) - m(6) * (g(1) * (t125 * t24 + t23 * t95 + t107) + g(2) * (t125 * t22 + t21 * t95 + t108) + g(3) * (t125 * t42 + t41 * t95 + t106)) + (-g(1) * (t128 * t23 + t24 * t87 + t107) - g(2) * (t128 * t21 + t22 * t87 + t108) - g(3) * (t128 * t41 + t42 * t87 + t106) - (-g(1) * t23 - g(2) * t21 - g(3) * t41) * t69 * t123) * m(7), -m(5) * (g(1) * (-rSges(5,1) * t6 - rSges(5,2) * t7) + g(2) * (-rSges(5,1) * t4 - rSges(5,2) * t5) + g(3) * (-rSges(5,1) * t25 - rSges(5,2) * t26)) - m(6) * (g(1) * (rSges(6,2) * t6 + t115 * t7 - t3) + g(2) * (rSges(6,2) * t4 + t115 * t5 - t2) + g(3) * (rSges(6,2) * t25 + t115 * t26 - t20)) + (-g(1) * (-t123 * t6 - t3) - g(2) * (-t123 * t4 - t2) - g(3) * (-t123 * t25 - t20) - (g(1) * t7 + g(2) * t5 + g(3) * t26) * t88) * m(7), t127 * (g(1) * t6 + g(2) * t4 + g(3) * t25) -m(7) * (g(1) * ((-t23 * t65 + t6 * t68) * rSges(7,1) + (-t23 * t68 - t6 * t65) * rSges(7,2)) + g(2) * ((-t21 * t65 + t4 * t68) * rSges(7,1) + (-t21 * t68 - t4 * t65) * rSges(7,2)) + g(3) * ((t25 * t68 - t41 * t65) * rSges(7,1) + (-t25 * t65 - t41 * t68) * rSges(7,2)))];
taug  = t1(:);
