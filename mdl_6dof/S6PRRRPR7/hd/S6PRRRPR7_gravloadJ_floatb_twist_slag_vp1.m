% Calculate Gravitation load on the joints for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:39:02
% EndTime: 2019-03-08 23:39:04
% DurationCPUTime: 1.28s
% Computational Cost: add. (1136->209), mult. (2978->319), div. (0->0), fcn. (3780->16), ass. (0->99)
t100 = cos(pkin(12));
t101 = cos(pkin(7));
t114 = sin(qJ(2));
t116 = cos(qJ(2));
t102 = cos(pkin(6));
t84 = t102 * t100;
t97 = sin(pkin(12));
t66 = t114 * t97 - t116 * t84;
t98 = sin(pkin(7));
t99 = sin(pkin(6));
t81 = t99 * t98;
t132 = t100 * t81 + t66 * t101;
t83 = t102 * t97;
t67 = t100 * t114 + t116 * t83;
t80 = t99 * t97;
t131 = t67 * t101 - t98 * t80;
t82 = t101 * t99;
t130 = t98 * t102 + t116 * t82;
t115 = cos(qJ(3));
t59 = sin(qJ(3));
t86 = t99 * t114;
t30 = t115 * t86 + t130 * t59;
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t65 = t101 * t102 - t116 * t81;
t17 = t30 * t60 + t58 * t65;
t41 = t114 * t84 + t116 * t97;
t13 = t115 * t41 - t132 * t59;
t61 = -t100 * t82 + t66 * t98;
t3 = t13 * t60 + t58 * t61;
t42 = t100 * t116 - t114 * t83;
t15 = t42 * t115 - t131 * t59;
t62 = t101 * t80 + t67 * t98;
t5 = t15 * t60 + t58 * t62;
t129 = g(1) * t5 + g(2) * t3 + g(3) * t17;
t16 = t30 * t58 - t60 * t65;
t2 = t13 * t58 - t60 * t61;
t4 = t15 * t58 - t60 * t62;
t128 = g(1) * t4 + g(2) * t2 + g(3) * t16;
t12 = t115 * t132 + t41 * t59;
t14 = t115 * t131 + t42 * t59;
t29 = -t115 * t130 + t59 * t86;
t127 = t58 * (-g(1) * t14 - g(2) * t12 - g(3) * t29);
t124 = -m(6) - m(7);
t118 = t60 * pkin(4);
t117 = rSges(5,3) + pkin(10);
t55 = sin(pkin(13));
t113 = t13 * t55;
t112 = t15 * t55;
t111 = t30 * t55;
t54 = pkin(13) + qJ(6);
t52 = sin(t54);
t110 = t52 * t60;
t53 = cos(t54);
t109 = t53 * t60;
t108 = t55 * t60;
t56 = cos(pkin(13));
t107 = t56 * t60;
t51 = pkin(5) * t56 + pkin(4);
t106 = t60 * t51;
t105 = rSges(7,3) + pkin(11) + qJ(5);
t72 = t114 * t81;
t87 = t116 * t99;
t104 = pkin(2) * t87 + pkin(9) * t72;
t103 = rSges(6,3) + qJ(5);
t73 = t114 * t82;
t37 = t115 * t87 - t59 * t73;
t96 = t37 * pkin(3) + t104;
t10 = t12 * pkin(3);
t95 = t13 * pkin(10) - t10;
t11 = t14 * pkin(3);
t94 = t15 * pkin(10) - t11;
t28 = t29 * pkin(3);
t93 = t30 * pkin(10) - t28;
t92 = t98 * pkin(9);
t91 = t58 * t98;
t90 = t59 * t101;
t89 = t60 * t98;
t88 = t101 * t115;
t85 = -rSges(5,1) * t60 + rSges(5,2) * t58;
t77 = t53 * rSges(7,1) - t52 * rSges(7,2) + t51;
t21 = -t115 * t66 - t41 * t90;
t39 = t66 * pkin(2);
t76 = t21 * pkin(3) + t41 * t92 - t39;
t23 = -t115 * t67 - t42 * t90;
t40 = t67 * pkin(2);
t75 = t23 * pkin(3) + t42 * t92 - t40;
t69 = t52 * rSges(7,1) + t53 * rSges(7,2) + t55 * pkin(5) + pkin(10);
t68 = rSges(4,3) * t98 + t92;
t36 = t115 * t73 + t59 * t87;
t25 = t37 * t60 + t58 * t72;
t24 = t37 * t58 - t60 * t72;
t22 = t42 * t88 - t59 * t67;
t20 = t41 * t88 - t59 * t66;
t9 = t23 * t60 + t42 * t91;
t8 = t23 * t58 - t42 * t89;
t7 = t21 * t60 + t41 * t91;
t6 = t21 * t58 - t41 * t89;
t1 = [(-m(2) - m(3) - m(4) - m(5) + t124) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t67 - t42 * rSges(3,2)) + g(2) * (-rSges(3,1) * t66 - t41 * rSges(3,2)) + g(3) * (rSges(3,1) * t87 - rSges(3,2) * t86)) - m(4) * (g(1) * (t23 * rSges(4,1) - t22 * rSges(4,2) + t42 * t68 - t40) + g(2) * (t21 * rSges(4,1) - t20 * rSges(4,2) + t41 * t68 - t39) + g(3) * (t37 * rSges(4,1) - t36 * rSges(4,2) + rSges(4,3) * t72 + t104)) - m(5) * (g(1) * (t9 * rSges(5,1) - t8 * rSges(5,2) + t117 * t22 + t75) + g(2) * (t7 * rSges(5,1) - t6 * rSges(5,2) + t117 * t20 + t76) + g(3) * (rSges(5,1) * t25 - rSges(5,2) * t24 + t117 * t36 + t96)) - m(6) * (g(1) * (t9 * pkin(4) + t22 * pkin(10) + (t22 * t55 + t56 * t9) * rSges(6,1) + (t22 * t56 - t55 * t9) * rSges(6,2) + t103 * t8 + t75) + g(2) * (t7 * pkin(4) + t20 * pkin(10) + (t20 * t55 + t56 * t7) * rSges(6,1) + (t20 * t56 - t55 * t7) * rSges(6,2) + t103 * t6 + t76) + g(3) * (t25 * pkin(4) + t36 * pkin(10) + (t25 * t56 + t36 * t55) * rSges(6,1) + (-t25 * t55 + t36 * t56) * rSges(6,2) + t103 * t24 + t96)) - m(7) * (g(1) * (t105 * t8 + t22 * t69 + t77 * t9 + t75) + g(2) * (t105 * t6 + t20 * t69 + t7 * t77 + t76) + g(3) * (t105 * t24 + t25 * t77 + t36 * t69 + t96)) -m(4) * (g(1) * (-rSges(4,1) * t14 - rSges(4,2) * t15) + g(2) * (-rSges(4,1) * t12 - rSges(4,2) * t13) + g(3) * (-rSges(4,1) * t29 - rSges(4,2) * t30)) - m(5) * (g(1) * (t117 * t15 + t14 * t85 - t11) + g(2) * (t117 * t13 + t12 * t85 - t10) + g(3) * (t117 * t30 + t29 * t85 - t28)) + (-g(1) * (-t14 * t106 + pkin(5) * t112 + (-t109 * t14 + t15 * t52) * rSges(7,1) + (t110 * t14 + t15 * t53) * rSges(7,2) + t94) - g(2) * (-t12 * t106 + pkin(5) * t113 + (-t109 * t12 + t13 * t52) * rSges(7,1) + (t110 * t12 + t13 * t53) * rSges(7,2) + t95) - g(3) * (-t29 * t106 + pkin(5) * t111 + (-t109 * t29 + t30 * t52) * rSges(7,1) + (t110 * t29 + t30 * t53) * rSges(7,2) + t93) - t105 * t127) * m(7) + (-g(1) * (-t14 * t118 + (-t107 * t14 + t112) * rSges(6,1) + (t108 * t14 + t15 * t56) * rSges(6,2) + t94) - g(2) * (-t12 * t118 + (-t107 * t12 + t113) * rSges(6,1) + (t108 * t12 + t13 * t56) * rSges(6,2) + t95) - g(3) * (-t29 * t118 + (-t107 * t29 + t111) * rSges(6,1) + (t108 * t29 + t30 * t56) * rSges(6,2) + t93) - t103 * t127) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t4 - rSges(5,2) * t5) + g(2) * (-rSges(5,1) * t2 - rSges(5,2) * t3) + g(3) * (-rSges(5,1) * t16 - rSges(5,2) * t17)) - m(6) * (t128 * (-t56 * rSges(6,1) + t55 * rSges(6,2) - pkin(4)) + t129 * t103) - m(7) * (t105 * t129 - t128 * t77) t124 * t128, -m(7) * (g(1) * ((t14 * t53 - t5 * t52) * rSges(7,1) + (-t14 * t52 - t5 * t53) * rSges(7,2)) + g(2) * ((t12 * t53 - t3 * t52) * rSges(7,1) + (-t12 * t52 - t3 * t53) * rSges(7,2)) + g(3) * ((-t17 * t52 + t29 * t53) * rSges(7,1) + (-t17 * t53 - t29 * t52) * rSges(7,2)))];
taug  = t1(:);
