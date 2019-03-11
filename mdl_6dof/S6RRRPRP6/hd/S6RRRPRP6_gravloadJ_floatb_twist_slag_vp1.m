% Calculate Gravitation load on the joints for
% S6RRRPRP6
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:59
% EndTime: 2019-03-09 16:54:03
% DurationCPUTime: 1.54s
% Computational Cost: add. (796->203), mult. (1387->293), div. (0->0), fcn. (1627->12), ass. (0->84)
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t89 = cos(pkin(6));
t57 = sin(pkin(6));
t62 = sin(qJ(2));
t99 = t57 * t62;
t127 = -t61 * t99 + t89 * t65;
t110 = cos(qJ(1));
t66 = cos(qJ(2));
t63 = sin(qJ(1));
t81 = t63 * t89;
t37 = t110 * t66 - t62 * t81;
t97 = t57 * t65;
t20 = -t37 * t61 + t63 * t97;
t76 = t89 * t110;
t35 = t62 * t76 + t63 * t66;
t56 = qJ(3) + pkin(11);
t53 = sin(t56);
t54 = cos(t56);
t84 = t57 * t110;
t15 = t35 * t54 - t53 * t84;
t34 = t62 * t63 - t66 * t76;
t60 = sin(qJ(5));
t64 = cos(qJ(5));
t3 = t15 * t60 - t34 * t64;
t109 = t34 * t60;
t126 = -t15 * t64 - t109;
t68 = t35 * t61 + t65 * t84;
t125 = pkin(3) * t68;
t112 = rSges(6,3) + pkin(10);
t93 = rSges(7,3) + qJ(6) + pkin(10);
t115 = g(2) * t34;
t36 = t110 * t62 + t66 * t81;
t124 = g(1) * t36 + t115;
t96 = t57 * t66;
t123 = -g(3) * t96 + t124;
t14 = t35 * t53 + t54 * t84;
t98 = t57 * t63;
t18 = t37 * t53 - t54 * t98;
t28 = t53 * t99 - t54 * t89;
t122 = g(1) * t18 + g(2) * t14 + g(3) * t28;
t51 = pkin(5) * t64 + pkin(4);
t121 = t51 * t54 + t93 * t53;
t120 = m(7) * t122;
t119 = pkin(4) * t54;
t113 = g(3) * t57;
t111 = -pkin(9) - rSges(4,3);
t107 = t35 * t60;
t105 = t36 * t60;
t104 = t37 * t60;
t101 = t54 * t60;
t100 = t54 * t64;
t95 = t60 * t66;
t94 = t64 * t66;
t52 = pkin(3) * t65 + pkin(2);
t59 = -qJ(4) - pkin(9);
t92 = -t34 * t52 - t35 * t59;
t91 = -t36 * t52 - t37 * t59;
t90 = t110 * pkin(1) + pkin(8) * t98;
t87 = t61 * t98;
t83 = -t63 * pkin(1) + pkin(8) * t84;
t45 = t61 * t84;
t82 = -t35 * t65 + t45;
t78 = t20 * pkin(3);
t77 = pkin(3) * t87 - t36 * t59 + t37 * t52 + t90;
t75 = rSges(5,1) * t54 - rSges(5,2) * t53;
t19 = t37 * t54 + t53 * t98;
t5 = -t19 * t60 + t36 * t64;
t74 = t127 * pkin(3);
t73 = rSges(4,1) * t65 - rSges(4,2) * t61 + pkin(2);
t29 = t53 * t89 + t54 * t99;
t12 = -t29 * t60 - t57 * t94;
t70 = pkin(3) * t45 + t34 * t59 - t35 * t52 + t83;
t38 = t52 * t96;
t23 = (t54 * t94 + t60 * t62) * t57;
t22 = (-t54 * t95 + t62 * t64) * t57;
t21 = t37 * t65 + t87;
t13 = -t29 * t64 + t57 * t95;
t10 = -t100 * t36 + t104;
t9 = t101 * t36 + t37 * t64;
t8 = -t100 * t34 + t107;
t7 = t101 * t34 + t35 * t64;
t6 = t19 * t64 + t105;
t1 = [-m(2) * (g(1) * (-t63 * rSges(2,1) - rSges(2,2) * t110) + g(2) * (rSges(2,1) * t110 - t63 * rSges(2,2))) - m(3) * (g(1) * (-t35 * rSges(3,1) + t34 * rSges(3,2) + rSges(3,3) * t84 + t83) + g(2) * (rSges(3,1) * t37 - rSges(3,2) * t36 + rSges(3,3) * t98 + t90)) - m(4) * (g(1) * (rSges(4,1) * t82 + rSges(4,2) * t68 - t35 * pkin(2) + t111 * t34 + t83) + g(2) * (rSges(4,1) * t21 + rSges(4,2) * t20 + pkin(2) * t37 - t111 * t36 + t90)) - m(5) * (g(1) * (-rSges(5,1) * t15 + rSges(5,2) * t14 - rSges(5,3) * t34 + t70) + g(2) * (rSges(5,1) * t19 - rSges(5,2) * t18 + rSges(5,3) * t36 + t77)) - m(6) * (g(1) * (rSges(6,1) * t126 + rSges(6,2) * t3 - pkin(4) * t15 - t112 * t14 + t70) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 + pkin(4) * t19 + t112 * t18 + t77)) - m(7) * (g(1) * (rSges(7,1) * t126 + rSges(7,2) * t3 - pkin(5) * t109 - t14 * t93 - t15 * t51 + t70) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + pkin(5) * t105 + t18 * t93 + t19 * t51 + t77)) -m(3) * (g(1) * (-rSges(3,1) * t36 - rSges(3,2) * t37) + g(2) * (-rSges(3,1) * t34 - rSges(3,2) * t35) + (rSges(3,1) * t66 - rSges(3,2) * t62) * t113) - m(4) * (g(1) * (-t111 * t37 - t73 * t36) - g(2) * t111 * t35 - t73 * t115 + (-t111 * t62 + t73 * t66) * t113) - m(5) * (g(1) * (rSges(5,3) * t37 - t75 * t36 + t91) + g(2) * (rSges(5,3) * t35 - t75 * t34 + t92) + g(3) * t38 + (t75 * t66 + (rSges(5,3) - t59) * t62) * t113) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9 - t36 * t119 + t91) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 - t34 * t119 + t92) + g(3) * (t23 * rSges(6,1) + t22 * rSges(6,2) + t96 * t119 - t59 * t99 + t38) - t123 * t53 * t112) - m(7) * (g(1) * (rSges(7,1) * t10 + rSges(7,2) * t9 + pkin(5) * t104 + t91) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + pkin(5) * t107 + t92) + g(3) * (t23 * rSges(7,1) + t22 * rSges(7,2) + t38) + ((pkin(5) * t60 - t59) * t62 + t121 * t66) * t113 - t124 * t121) -m(4) * (g(1) * (rSges(4,1) * t20 - rSges(4,2) * t21) + g(2) * (-rSges(4,1) * t68 + rSges(4,2) * t82) + g(3) * (t127 * rSges(4,1) + (-t61 * t89 - t62 * t97) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t18 - rSges(5,2) * t19 + t78) + g(2) * (-t14 * rSges(5,1) - t15 * rSges(5,2) - t125) + g(3) * (-rSges(5,1) * t28 - rSges(5,2) * t29 + t74)) - m(7) * (g(1) * (t19 * t93 + t78) + g(2) * (t15 * t93 - t125) + g(3) * (t29 * t93 + t74)) - (-rSges(7,1) * t64 + rSges(7,2) * t60 - t51) * t120 + (-g(1) * (t112 * t19 + t78) - g(2) * (t112 * t15 - t125) - g(3) * (t112 * t29 + t74) - t122 * (-rSges(6,1) * t64 + rSges(6,2) * t60 - pkin(4))) * m(6) (-m(5) - m(6) - m(7)) * t123, -m(6) * (g(1) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t3 + rSges(6,2) * t126) + g(3) * (rSges(6,1) * t12 + rSges(6,2) * t13)) + (-g(1) * (rSges(7,1) * t5 - rSges(7,2) * t6) - g(2) * (-rSges(7,1) * t3 + rSges(7,2) * t126) - g(3) * (t12 * rSges(7,1) + t13 * rSges(7,2)) - (g(1) * t5 - g(2) * t3 + g(3) * t12) * pkin(5)) * m(7), -t120];
taug  = t1(:);
