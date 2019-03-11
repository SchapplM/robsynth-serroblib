% Calculate Gravitation load on the joints for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:44:09
% EndTime: 2019-03-09 06:44:14
% DurationCPUTime: 1.62s
% Computational Cost: add. (1373->187), mult. (3717->272), div. (0->0), fcn. (4786->14), ass. (0->87)
t132 = cos(qJ(3));
t119 = sin(pkin(7));
t120 = sin(pkin(6));
t105 = t120 * t119;
t122 = cos(pkin(7));
t133 = cos(qJ(1));
t121 = cos(pkin(12));
t123 = cos(pkin(6));
t108 = t123 * t121;
t118 = sin(pkin(12));
t131 = sin(qJ(1));
t96 = -t108 * t133 + t131 * t118;
t143 = t133 * t105 + t96 * t122;
t107 = t123 * t118;
t64 = t107 * t133 + t121 * t131;
t78 = sin(qJ(3));
t44 = -t132 * t64 + t143 * t78;
t106 = t122 * t120;
t56 = -t133 * t106 + t96 * t119;
t77 = sin(qJ(4));
t80 = cos(qJ(4));
t20 = t44 * t80 - t56 * t77;
t41 = t143 * t132 + t64 * t78;
t76 = sin(qJ(5));
t79 = cos(qJ(5));
t3 = t20 * t76 + t41 * t79;
t4 = t20 * t79 - t41 * t76;
t142 = t44 * t77 + t56 * t80;
t91 = t108 * t131 + t118 * t133;
t82 = t131 * t106 + t91 * t119;
t135 = pkin(11) + rSges(7,2);
t141 = -t131 * t105 + t91 * t122;
t139 = t121 * t106 + t119 * t123;
t124 = rSges(7,3) + qJ(6);
t137 = pkin(5) + rSges(7,1);
t97 = t124 * t76 + t137 * t79;
t138 = pkin(4) * t80;
t136 = pkin(10) + rSges(5,3);
t134 = pkin(11) + rSges(6,3);
t130 = t41 * t77;
t65 = -t107 * t131 + t121 * t133;
t45 = t141 * t132 + t65 * t78;
t129 = t45 * t77;
t104 = t120 * t118;
t54 = t104 * t78 - t139 * t132;
t128 = t54 * t77;
t127 = t76 * t80;
t126 = t79 * t80;
t114 = t120 * t131;
t125 = t133 * pkin(1) + qJ(2) * t114;
t115 = t133 * t120;
t117 = -pkin(1) * t131 + qJ(2) * t115;
t113 = -rSges(5,1) * t80 + rSges(5,2) * t77;
t112 = rSges(6,1) * t79 - rSges(6,2) * t76;
t35 = t41 * pkin(3);
t111 = -pkin(10) * t44 - pkin(11) * t130 - t41 * t138 - t35;
t37 = t45 * pkin(3);
t46 = t65 * t132 - t141 * t78;
t110 = t46 * pkin(10) - pkin(11) * t129 - t45 * t138 - t37;
t53 = t54 * pkin(3);
t55 = t132 * t104 + t139 * t78;
t109 = t55 * pkin(10) - pkin(11) * t128 - t54 * t138 - t53;
t88 = -t64 * pkin(2) - t56 * pkin(9) + t117;
t86 = t44 * pkin(3) + t88;
t85 = t65 * pkin(2) + t82 * pkin(9) + t125;
t84 = t46 * pkin(3) + t85;
t83 = t20 * pkin(4) - pkin(10) * t41 + t86;
t22 = t46 * t80 + t77 * t82;
t81 = t22 * pkin(4) + t45 * pkin(10) + t84;
t63 = -t105 * t121 + t122 * t123;
t40 = t55 * t80 + t63 * t77;
t39 = -t55 * t77 + t63 * t80;
t34 = t39 * pkin(4);
t24 = -t126 * t54 + t55 * t76;
t23 = -t127 * t54 - t55 * t79;
t21 = t46 * t77 - t80 * t82;
t15 = t21 * pkin(4);
t13 = t142 * pkin(4);
t12 = t40 * t79 + t54 * t76;
t11 = t40 * t76 - t54 * t79;
t10 = -t126 * t45 + t46 * t76;
t9 = -t127 * t45 - t46 * t79;
t8 = -t126 * t41 - t44 * t76;
t7 = -t127 * t41 + t44 * t79;
t6 = t22 * t79 + t45 * t76;
t5 = t22 * t76 - t45 * t79;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t131 - rSges(2,2) * t133) + g(2) * (rSges(2,1) * t133 - rSges(2,2) * t131)) - m(3) * (g(1) * (-t64 * rSges(3,1) + rSges(3,2) * t96 + rSges(3,3) * t115 + t117) + g(2) * (t65 * rSges(3,1) - rSges(3,2) * t91 + rSges(3,3) * t114 + t125)) - m(4) * (g(1) * (t44 * rSges(4,1) + rSges(4,2) * t41 - rSges(4,3) * t56 + t88) + g(2) * (t46 * rSges(4,1) - t45 * rSges(4,2) + rSges(4,3) * t82 + t85)) - m(5) * (g(1) * (t20 * rSges(5,1) - rSges(5,2) * t142 - t136 * t41 + t86) + g(2) * (t22 * rSges(5,1) - t21 * rSges(5,2) + t136 * t45 + t84)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t134 * t142 + t83) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t134 * t21 + t81)) - m(7) * (g(1) * (t124 * t3 + t135 * t142 + t137 * t4 + t83) + g(2) * (t124 * t5 + t135 * t21 + t137 * t6 + t81)) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t114 - g(2) * t115 + g(3) * t123) -m(4) * (g(1) * (-rSges(4,1) * t45 - rSges(4,2) * t46) + g(2) * (-rSges(4,1) * t41 + rSges(4,2) * t44) + g(3) * (-rSges(4,1) * t54 - rSges(4,2) * t55)) - m(5) * (g(1) * (t113 * t45 + t136 * t46 - t37) + g(2) * (t113 * t41 - t136 * t44 - t35) + g(3) * (t113 * t54 + t136 * t55 - t53)) - m(6) * (g(1) * (rSges(6,1) * t10 - rSges(6,2) * t9 - rSges(6,3) * t129 + t110) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t130 + t111) + g(3) * (rSges(6,1) * t24 - rSges(6,2) * t23 - rSges(6,3) * t128 + t109)) - m(7) * (g(1) * (-rSges(7,2) * t129 + t10 * t137 + t124 * t9 + t110) + g(2) * (-rSges(7,2) * t130 + t124 * t7 + t137 * t8 + t111) + g(3) * (-rSges(7,2) * t128 + t124 * t23 + t137 * t24 + t109)) -m(5) * (g(1) * (-rSges(5,1) * t21 - rSges(5,2) * t22) + g(2) * (rSges(5,1) * t142 + rSges(5,2) * t20) + g(3) * (rSges(5,1) * t39 - rSges(5,2) * t40)) - m(6) * (g(1) * (-t112 * t21 + t134 * t22 - t15) + g(2) * (t112 * t142 - t134 * t20 + t13) + g(3) * (t112 * t39 + t134 * t40 + t34)) - m(7) * ((t135 * t40 + t97 * t39 + t34) * g(3) + (-t135 * t20 + t142 * t97 + t13) * g(2) + (t135 * t22 - t97 * t21 - t15) * g(1)) -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(3) * (-rSges(6,1) * t11 - rSges(6,2) * t12)) - m(7) * (g(1) * (t124 * t6 - t137 * t5) + g(2) * (-t124 * t4 + t137 * t3) + g(3) * (-t137 * t11 + t12 * t124)) -m(7) * (g(1) * t5 - g(2) * t3 + g(3) * t11)];
taug  = t1(:);
