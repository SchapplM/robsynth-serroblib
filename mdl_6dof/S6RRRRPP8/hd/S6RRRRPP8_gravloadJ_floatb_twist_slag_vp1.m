% Calculate Gravitation load on the joints for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:20
% EndTime: 2019-03-09 21:30:23
% DurationCPUTime: 1.36s
% Computational Cost: add. (833->214), mult. (2047->302), div. (0->0), fcn. (2509->10), ass. (0->91)
t75 = sin(qJ(4));
t78 = cos(qJ(4));
t141 = pkin(4) * t78 + qJ(5) * t75;
t130 = cos(qJ(1));
t74 = sin(pkin(6));
t106 = t74 * t130;
t129 = cos(qJ(3));
t128 = sin(qJ(1));
t77 = sin(qJ(2));
t79 = cos(qJ(2));
t114 = cos(pkin(6));
t93 = t114 * t130;
t57 = t128 * t79 + t77 * t93;
t76 = sin(qJ(3));
t30 = -t76 * t106 + t57 * t129;
t56 = t128 * t77 - t79 * t93;
t6 = t30 * t75 - t56 * t78;
t7 = t30 * t78 + t56 * t75;
t134 = rSges(7,1) + pkin(5);
t117 = rSges(7,2) + qJ(5);
t105 = t74 * t129;
t100 = -t130 * t105 - t57 * t76;
t104 = t74 * t128;
t92 = t114 * t128;
t59 = t130 * t79 - t77 * t92;
t33 = -t129 * t104 + t59 * t76;
t122 = t74 * t77;
t54 = t114 * t129 - t76 * t122;
t140 = m(7) * (-g(1) * t33 + g(2) * t100 + g(3) * t54);
t135 = g(3) * t74;
t133 = rSges(6,2) + pkin(10);
t132 = rSges(4,3) + pkin(9);
t131 = rSges(5,3) + pkin(10);
t127 = rSges(7,2) * t75;
t125 = t56 * t76;
t58 = t130 * t77 + t79 * t92;
t123 = t58 * t76;
t121 = t74 * t79;
t120 = t130 * pkin(1) + pkin(8) * t104;
t119 = pkin(2) * t121 + pkin(9) * t122;
t116 = rSges(6,3) + qJ(5);
t115 = rSges(7,3) + qJ(6);
t113 = t76 * t121;
t112 = pkin(10) - t115;
t23 = t100 * pkin(3);
t111 = t100 * t141 + t23;
t25 = t33 * pkin(3);
t110 = -t141 * t33 - t25;
t49 = t54 * pkin(3);
t109 = t141 * t54 + t49;
t108 = t59 * pkin(2) + t120;
t107 = pkin(3) * t129;
t103 = t75 * t129;
t102 = t78 * t129;
t101 = t79 * t129;
t96 = t74 * t101;
t99 = pkin(3) * t96 + pkin(10) * t113 + t119;
t98 = g(1) * t112;
t97 = g(2) * t112;
t95 = -t128 * pkin(1) + pkin(8) * t106;
t37 = (t78 * t101 + t75 * t77) * t74;
t94 = t37 * pkin(4) + t99;
t91 = rSges(5,1) * t78 - rSges(5,2) * t75;
t90 = rSges(6,1) * t78 + rSges(6,3) * t75;
t50 = t56 * pkin(2);
t89 = t57 * pkin(9) - pkin(10) * t125 - t56 * t107 - t50;
t52 = t58 * pkin(2);
t88 = t59 * pkin(9) - pkin(10) * t123 - t58 * t107 - t52;
t87 = -t57 * pkin(2) + t95;
t15 = -t56 * t102 + t57 * t75;
t86 = t15 * pkin(4) + t89;
t17 = -t58 * t102 + t59 * t75;
t85 = t17 * pkin(4) + t88;
t34 = t76 * t104 + t59 * t129;
t84 = t34 * pkin(3) + pkin(9) * t58 + t108;
t83 = t129 * rSges(4,1) - rSges(4,2) * t76;
t11 = t34 * t78 + t58 * t75;
t82 = t11 * pkin(4) + t84;
t81 = -pkin(3) * t30 - t56 * pkin(9) + t87;
t80 = -pkin(4) * t7 + t81;
t55 = t77 * t105 + t114 * t76;
t36 = -t78 * t122 + t75 * t96;
t28 = -t75 * t121 + t55 * t78;
t27 = t78 * t121 + t55 * t75;
t22 = t27 * pkin(4);
t16 = -t58 * t103 - t59 * t78;
t14 = -t56 * t103 - t57 * t78;
t10 = t34 * t75 - t58 * t78;
t4 = t10 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-t128 * rSges(2,1) - t130 * rSges(2,2)) + g(2) * (t130 * rSges(2,1) - t128 * rSges(2,2))) - m(3) * (g(1) * (-t57 * rSges(3,1) + t56 * rSges(3,2) + rSges(3,3) * t106 + t95) + g(2) * (t59 * rSges(3,1) - t58 * rSges(3,2) + rSges(3,3) * t104 + t120)) - m(4) * (g(1) * (-rSges(4,1) * t30 - rSges(4,2) * t100 - t132 * t56 + t87) + g(2) * (rSges(4,1) * t34 - rSges(4,2) * t33 + t132 * t58 + t108)) - m(5) * (g(1) * (-rSges(5,1) * t7 + rSges(5,2) * t6 + t100 * t131 + t81) + g(2) * (rSges(5,1) * t11 - rSges(5,2) * t10 + t131 * t33 + t84)) - m(6) * (g(1) * (-rSges(6,1) * t7 + t100 * t133 - t116 * t6 + t80) + g(2) * (rSges(6,1) * t11 + t116 * t10 + t133 * t33 + t82)) - m(7) * (g(1) * (-t117 * t6 - t134 * t7 + t80) + g(2) * (t117 * t10 + t134 * t11 + t82) + t33 * t97 + t100 * t98) -m(3) * (g(1) * (-rSges(3,1) * t58 - rSges(3,2) * t59) + g(2) * (-rSges(3,1) * t56 - rSges(3,2) * t57) + (rSges(3,1) * t79 - rSges(3,2) * t77) * t135) - m(4) * (g(1) * (t132 * t59 - t83 * t58 - t52) + g(2) * (t132 * t57 - t83 * t56 - t50) + g(3) * t119 + (rSges(4,3) * t77 + t83 * t79) * t135) - m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t16 - rSges(5,3) * t123 + t88) + g(2) * (rSges(5,1) * t15 - rSges(5,2) * t14 - rSges(5,3) * t125 + t89) + g(3) * (t37 * rSges(5,1) - t36 * rSges(5,2) + rSges(5,3) * t113 + t99)) - m(6) * (g(1) * (rSges(6,1) * t17 - rSges(6,2) * t123 + t116 * t16 + t85) + g(2) * (rSges(6,1) * t15 - rSges(6,2) * t125 + t116 * t14 + t86) + g(3) * (t37 * rSges(6,1) + rSges(6,2) * t113 + t116 * t36 + t94)) - m(7) * (g(1) * (t117 * t16 + t134 * t17 + t85) + g(2) * (t117 * t14 + t134 * t15 + t86) + g(3) * (t117 * t36 + t134 * t37 + t94) + (g(1) * t58 + g(2) * t56 - g(3) * t121) * t76 * t115) -m(4) * (g(1) * (-rSges(4,1) * t33 - rSges(4,2) * t34) + g(2) * (rSges(4,1) * t100 - rSges(4,2) * t30) + g(3) * (rSges(4,1) * t54 - rSges(4,2) * t55)) - m(5) * (g(1) * (t131 * t34 - t91 * t33 - t25) + g(2) * (t100 * t91 + t131 * t30 + t23) + g(3) * (t131 * t55 + t91 * t54 + t49)) - m(6) * (g(1) * (t133 * t34 - t90 * t33 + t110) + g(2) * (t100 * t90 + t133 * t30 + t111) + g(3) * (t133 * t55 + t90 * t54 + t109)) - m(7) * (g(1) * (-t33 * t127 + t110) + g(2) * (t100 * t127 + t111) + t34 * t98 + t30 * t97 + (t112 * t55 + t54 * t127 + t109) * g(3)) - t78 * t134 * t140, -m(5) * (g(1) * (-rSges(5,1) * t10 - rSges(5,2) * t11) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7) + g(3) * (-rSges(5,1) * t27 - rSges(5,2) * t28)) - m(6) * (g(1) * (-rSges(6,1) * t10 + t116 * t11 - t4) + g(2) * (-rSges(6,1) * t6 + t116 * t7 - t2) + g(3) * (-rSges(6,1) * t27 + t116 * t28 - t22)) - m(7) * (g(1) * (-t134 * t10 + t117 * t11 - t4) + g(2) * (t117 * t7 - t134 * t6 - t2) + g(3) * (t117 * t28 - t134 * t27 - t22)) (-m(6) - m(7)) * (g(1) * t10 + g(2) * t6 + g(3) * t27) -t140];
taug  = t1(:);
