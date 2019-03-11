% Calculate Gravitation load on the joints for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:52
% EndTime: 2019-03-08 22:42:55
% DurationCPUTime: 1.20s
% Computational Cost: add. (632->159), mult. (1305->231), div. (0->0), fcn. (1561->12), ass. (0->82)
t59 = sin(qJ(4));
t115 = pkin(4) * t59;
t126 = pkin(8) + t115;
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t86 = sin(pkin(10));
t88 = cos(pkin(6));
t71 = t88 * t86;
t87 = cos(pkin(10));
t39 = -t61 * t71 + t64 * t87;
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t57 = sin(pkin(6));
t81 = t57 * t86;
t19 = t39 * t63 + t60 * t81;
t38 = t61 * t87 + t64 * t71;
t62 = cos(qJ(4));
t125 = -t19 * t59 + t38 * t62;
t72 = t88 * t87;
t37 = t61 * t72 + t64 * t86;
t82 = t57 * t87;
t17 = t37 * t63 - t60 * t82;
t36 = t61 * t86 - t64 * t72;
t124 = -t17 * t59 + t36 * t62;
t16 = t37 * t60 + t63 * t82;
t18 = t39 * t60 - t63 * t81;
t98 = t57 * t61;
t40 = t60 * t98 - t63 * t88;
t123 = -g(1) * t18 - g(2) * t16 - g(3) * t40;
t122 = rSges(7,1) + pkin(5);
t89 = rSges(7,3) + qJ(6);
t120 = g(1) * t38 + g(2) * t36;
t97 = t57 * t64;
t118 = t60 * (g(3) * t97 - t120);
t106 = rSges(5,3) + pkin(9);
t70 = rSges(5,1) * t62 - rSges(5,2) * t59 + pkin(3);
t117 = t106 * t60 + t63 * t70;
t116 = -m(6) - m(7);
t109 = g(3) * t57;
t107 = rSges(4,3) + pkin(8);
t53 = pkin(4) * t62 + pkin(3);
t101 = t53 * t63;
t56 = qJ(4) + pkin(11);
t54 = sin(t56);
t100 = t54 * t63;
t55 = cos(t56);
t99 = t55 * t63;
t96 = t63 * t64;
t58 = -qJ(5) - pkin(9);
t93 = -t16 * t53 - t17 * t58;
t92 = -t18 * t53 - t19 * t58;
t41 = t60 * t88 + t63 * t98;
t91 = -t40 * t53 - t41 * t58;
t90 = pkin(2) * t97 + pkin(8) * t98;
t84 = t54 * t97;
t83 = g(3) * t90;
t80 = t57 * t53 * t96 + t98 * t115 + t90;
t79 = t124 * pkin(4);
t78 = t125 * pkin(4);
t77 = rSges(4,1) * t63 - rSges(4,2) * t60;
t76 = t59 * rSges(5,1) + t62 * rSges(5,2);
t75 = -rSges(6,1) * t55 + rSges(6,2) * t54;
t34 = t36 * pkin(2);
t74 = -t36 * t101 + t126 * t37 - t34;
t35 = t38 * pkin(2);
t73 = -t38 * t101 + t126 * t39 - t35;
t69 = pkin(8) + t76;
t68 = -t41 * t59 - t62 * t97;
t67 = t68 * pkin(4);
t21 = (t54 * t61 + t55 * t96) * t57;
t20 = -t55 * t98 + t63 * t84;
t13 = t41 * t55 - t84;
t12 = t41 * t54 + t55 * t97;
t9 = -t38 * t99 + t39 * t54;
t8 = -t100 * t38 - t39 * t55;
t7 = -t36 * t99 + t37 * t54;
t6 = -t100 * t36 - t37 * t55;
t5 = t19 * t55 + t38 * t54;
t4 = t19 * t54 - t38 * t55;
t3 = t17 * t55 + t36 * t54;
t2 = t17 * t54 - t36 * t55;
t1 = [(-m(2) - m(3) - m(4) - m(5) + t116) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t38 - rSges(3,2) * t39) + g(2) * (-rSges(3,1) * t36 - t37 * rSges(3,2)) + (rSges(3,1) * t64 - rSges(3,2) * t61) * t109) - m(4) * (g(1) * (t107 * t39 - t38 * t77 - t35) + g(2) * (t107 * t37 - t36 * t77 - t34) + t83 + (rSges(4,3) * t61 + t64 * t77) * t109) - m(5) * (t83 + (t117 * t64 + t76 * t61) * t109 - t120 * t117 + (t69 * t37 - t34) * g(2) + (t69 * t39 - t35) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t73) + g(2) * (rSges(6,1) * t7 - rSges(6,2) * t6 + t74) + g(3) * (t21 * rSges(6,1) - t20 * rSges(6,2) + t80) + (rSges(6,3) - t58) * t118) - m(7) * (g(1) * (t122 * t9 + t89 * t8 + t73) + g(2) * (t122 * t7 + t89 * t6 + t74) + g(3) * (t122 * t21 + t89 * t20 + t80) + (rSges(7,2) - t58) * t118) -m(4) * (g(1) * (-rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (-rSges(4,1) * t16 - rSges(4,2) * t17) + g(3) * (-rSges(4,1) * t40 - rSges(4,2) * t41)) - m(5) * (t123 * t70 + (g(1) * t19 + g(2) * t17 + g(3) * t41) * t106) - m(6) * (g(1) * (rSges(6,3) * t19 + t18 * t75 + t92) + g(2) * (rSges(6,3) * t17 + t16 * t75 + t93) + g(3) * (rSges(6,3) * t41 + t40 * t75 + t91)) + (-g(1) * (rSges(7,2) * t19 + t92) - g(2) * (rSges(7,2) * t17 + t93) - g(3) * (rSges(7,2) * t41 + t91) + t123 * (-t122 * t55 - t54 * t89)) * m(7), -m(5) * (g(1) * (t125 * rSges(5,1) + (-t19 * t62 - t38 * t59) * rSges(5,2)) + g(2) * (t124 * rSges(5,1) + (-t17 * t62 - t36 * t59) * rSges(5,2)) + g(3) * (t68 * rSges(5,1) + (-t41 * t62 + t59 * t97) * rSges(5,2))) - m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5 + t78) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3 + t79) + g(3) * (-t12 * rSges(6,1) - t13 * rSges(6,2) + t67)) - m(7) * (g(1) * (-t122 * t4 + t5 * t89 + t78) + g(2) * (-t122 * t2 + t3 * t89 + t79) + g(3) * (-t12 * t122 + t13 * t89 + t67)) -t116 * t123, -m(7) * (g(1) * t4 + g(2) * t2 + g(3) * t12)];
taug  = t1(:);
