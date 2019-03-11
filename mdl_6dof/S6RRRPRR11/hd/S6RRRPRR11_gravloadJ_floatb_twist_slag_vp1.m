% Calculate Gravitation load on the joints for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:53
% EndTime: 2019-03-09 19:24:57
% DurationCPUTime: 1.47s
% Computational Cost: add. (846->211), mult. (2092->296), div. (0->0), fcn. (2585->12), ass. (0->95)
t116 = cos(pkin(6));
t74 = sin(pkin(6));
t78 = sin(qJ(2));
t123 = t74 * t78;
t77 = sin(qJ(3));
t81 = cos(qJ(3));
t54 = -t116 * t81 + t123 * t77;
t55 = t116 * t77 + t123 * t81;
t76 = sin(qJ(5));
t80 = cos(qJ(5));
t20 = t54 * t80 - t55 * t76;
t75 = sin(qJ(6));
t79 = cos(qJ(6));
t87 = t79 * rSges(7,1) - t75 * rSges(7,2) + pkin(5);
t145 = t20 * g(3) * t87;
t129 = rSges(7,3) + pkin(11);
t21 = t54 * t76 + t55 * t80;
t144 = t129 * t21;
t126 = sin(qJ(1));
t109 = t74 * t126;
t127 = cos(qJ(2));
t128 = cos(qJ(1));
t93 = t116 * t126;
t59 = t127 * t128 - t78 * t93;
t36 = -t109 * t81 + t59 * t77;
t37 = t109 * t77 + t59 * t81;
t13 = t36 * t76 + t37 * t80;
t90 = -t36 * t80 + t37 * t76;
t143 = -rSges(6,1) * t90 - rSges(6,2) * t13;
t142 = t129 * t13 - t87 * t90;
t111 = t74 * t128;
t94 = t116 * t128;
t57 = t126 * t127 + t78 * t94;
t32 = t111 * t81 + t57 * t77;
t33 = -t77 * t111 + t57 * t81;
t138 = t32 * t80 - t33 * t76;
t7 = t32 * t76 + t33 * t80;
t141 = rSges(6,1) * t138 - rSges(6,2) * t7;
t140 = t129 * t7 + t138 * t87;
t139 = rSges(6,1) * t20 - rSges(6,2) * t21;
t135 = -t76 * t81 + t77 * t80;
t134 = pkin(10) - pkin(9);
t133 = g(3) * t74;
t132 = rSges(5,2) + pkin(9);
t131 = rSges(4,3) + pkin(9);
t130 = -rSges(6,3) - pkin(10);
t56 = t126 * t78 - t127 * t94;
t125 = t56 * t81;
t58 = t127 * t93 + t128 * t78;
t124 = t58 * t81;
t120 = t128 * pkin(1) + pkin(8) * t109;
t110 = t74 * t127;
t119 = pkin(2) * t110 + pkin(9) * t123;
t118 = qJ(4) * t77;
t117 = rSges(5,3) + qJ(4);
t115 = -pkin(9) - t130;
t50 = t56 * pkin(2);
t114 = -pkin(3) * t125 - t56 * t118 - t50;
t52 = t58 * pkin(2);
t113 = -pkin(3) * t124 - t58 * t118 - t52;
t112 = t59 * pkin(2) + t120;
t108 = t77 * t127;
t107 = t81 * t127;
t106 = t37 * pkin(3) + t112;
t105 = -pkin(4) * t125 + t114;
t104 = -pkin(4) * t124 + t113;
t100 = t74 * t107;
t101 = t74 * t108;
t103 = pkin(3) * t100 + qJ(4) * t101 + t119;
t102 = g(2) * t115;
t25 = t32 * pkin(3);
t99 = -t32 * pkin(4) + qJ(4) * t33 - t25;
t29 = t36 * pkin(3);
t98 = -t36 * pkin(4) + qJ(4) * t37 - t29;
t49 = t54 * pkin(3);
t97 = -t54 * pkin(4) + qJ(4) * t55 - t49;
t96 = -pkin(1) * t126 + pkin(8) * t111;
t95 = pkin(4) * t100 + t103;
t92 = -rSges(4,1) * t81 + rSges(4,2) * t77;
t91 = -rSges(5,1) * t81 - rSges(5,3) * t77;
t89 = -t76 * t77 - t80 * t81;
t88 = -t57 * pkin(2) + t96;
t86 = -pkin(3) * t33 + t88;
t85 = t75 * rSges(7,1) + t79 * rSges(7,2) + t134;
t83 = t37 * pkin(4) + qJ(4) * t36 + t106;
t82 = -pkin(4) * t33 - qJ(4) * t32 + t86;
t39 = (t107 * t80 + t108 * t76) * t74;
t38 = t100 * t76 - t101 * t80;
t17 = t89 * t58;
t16 = t135 * t58;
t15 = t89 * t56;
t14 = t135 * t56;
t2 = t13 * t79 - t58 * t75;
t1 = -t13 * t75 - t58 * t79;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t126 - rSges(2,2) * t128) + g(2) * (rSges(2,1) * t128 - rSges(2,2) * t126)) - m(3) * (g(1) * (-t57 * rSges(3,1) + t56 * rSges(3,2) + rSges(3,3) * t111 + t96) + g(2) * (t59 * rSges(3,1) - t58 * rSges(3,2) + rSges(3,3) * t109 + t120)) - m(4) * (g(1) * (-rSges(4,1) * t33 + rSges(4,2) * t32 - t131 * t56 + t88) + g(2) * (rSges(4,1) * t37 - rSges(4,2) * t36 + t131 * t58 + t112)) - m(5) * (g(1) * (-rSges(5,1) * t33 - t117 * t32 - t132 * t56 + t86) + g(2) * (rSges(5,1) * t37 + t117 * t36 + t132 * t58 + t106)) - m(6) * (g(2) * (rSges(6,1) * t13 - rSges(6,2) * t90 + t83) - t58 * t102 + (-rSges(6,1) * t7 - rSges(6,2) * t138 + t115 * t56 + t82) * g(1)) - m(7) * (g(1) * (t129 * t138 + t56 * t85 - t7 * t87 + t82) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t13 + t129 * t90 - t134 * t58 + t83)) -m(3) * (g(1) * (-rSges(3,1) * t58 - rSges(3,2) * t59) + g(2) * (-rSges(3,1) * t56 - rSges(3,2) * t57) + (rSges(3,1) * t127 - rSges(3,2) * t78) * t133) - m(4) * (g(1) * (t131 * t59 + t58 * t92 - t52) + g(2) * (t131 * t57 + t56 * t92 - t50) + g(3) * t119 + (rSges(4,1) * t107 - rSges(4,2) * t108 + rSges(4,3) * t78) * t133) - m(5) * (g(1) * (t132 * t59 + t58 * t91 + t113) + g(2) * (t132 * t57 + t56 * t91 + t114) + g(3) * t103 + (rSges(5,1) * t107 + rSges(5,2) * t78 + rSges(5,3) * t108) * t133) - m(6) * (g(2) * (rSges(6,1) * t15 - rSges(6,2) * t14 + t105) + g(3) * (rSges(6,1) * t39 - rSges(6,2) * t38 + t123 * t130 + t95) - t57 * t102 + (rSges(6,1) * t17 - rSges(6,2) * t16 - t115 * t59 + t104) * g(1)) - m(7) * (g(1) * (t129 * t16 + t17 * t87 - t59 * t85 + t104) + g(2) * (t129 * t14 + t15 * t87 - t57 * t85 + t105) + g(3) * (t39 * pkin(5) - pkin(10) * t123 + (-t123 * t75 + t39 * t79) * rSges(7,1) + (-t123 * t79 - t39 * t75) * rSges(7,2) + t129 * t38 + t95)) -m(4) * (g(1) * (-rSges(4,1) * t36 - rSges(4,2) * t37) + g(2) * (-rSges(4,1) * t32 - rSges(4,2) * t33) + g(3) * (-rSges(4,1) * t54 - rSges(4,2) * t55)) - m(5) * (g(1) * (-rSges(5,1) * t36 + t117 * t37 - t29) + g(2) * (-rSges(5,1) * t32 + t117 * t33 - t25) + g(3) * (-rSges(5,1) * t54 + t117 * t55 - t49)) - m(6) * (g(1) * (-t143 + t98) + g(2) * (-t141 + t99) + g(3) * (-t139 + t97)) - m(7) * (g(3) * (t97 - t144) - t145 + (-t140 + t99) * g(2) + (-t142 + t98) * g(1)) (-m(5) - m(6) - m(7)) * (g(1) * t36 + g(2) * t32 + g(3) * t54) -m(6) * (g(1) * t143 + g(2) * t141 + g(3) * t139) - m(7) * (t142 * g(1) + g(2) * t140 + g(3) * t144 + t145) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * ((-t56 * t79 - t7 * t75) * rSges(7,1) + (t56 * t75 - t7 * t79) * rSges(7,2)) + g(3) * ((t110 * t79 - t21 * t75) * rSges(7,1) + (-t110 * t75 - t21 * t79) * rSges(7,2)))];
taug  = t3(:);
