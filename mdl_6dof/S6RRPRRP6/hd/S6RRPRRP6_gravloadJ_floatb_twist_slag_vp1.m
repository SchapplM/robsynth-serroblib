% Calculate Gravitation load on the joints for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:48
% EndTime: 2019-03-09 12:07:50
% DurationCPUTime: 1.32s
% Computational Cost: add. (990->195), mult. (2437->282), div. (0->0), fcn. (3097->12), ass. (0->91)
t82 = cos(qJ(2));
t83 = cos(qJ(1));
t114 = t82 * t83;
t78 = sin(qJ(2));
t79 = sin(qJ(1));
t117 = t79 * t78;
t75 = cos(pkin(6));
t131 = t75 * t114 - t117;
t110 = sin(pkin(6));
t103 = t83 * t110;
t111 = cos(pkin(11));
t74 = sin(pkin(11));
t93 = -t78 * t111 - t82 * t74;
t57 = t93 * t75;
t92 = t82 * t111 - t78 * t74;
t40 = -t57 * t83 + t79 * t92;
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t18 = -t77 * t103 + t40 * t81;
t56 = t92 * t75;
t39 = t83 * t56 + t79 * t93;
t76 = sin(qJ(5));
t80 = cos(qJ(5));
t1 = t18 * t76 + t39 * t80;
t2 = t18 * t80 - t39 * t76;
t128 = rSges(7,2) + pkin(10);
t112 = rSges(7,3) + qJ(6);
t129 = rSges(7,1) + pkin(5);
t89 = t112 * t76 + t129 * t80;
t130 = pkin(4) * t81;
t127 = rSges(5,3) + pkin(9);
t126 = rSges(6,3) + pkin(10);
t124 = t39 * t77;
t42 = -t56 * t79 + t83 * t93;
t122 = t42 * t77;
t106 = t78 * t110;
t96 = t111 * t110;
t54 = t74 * t106 - t82 * t96;
t121 = t54 * t77;
t120 = t76 * t81;
t119 = t78 * t83;
t73 = pkin(2) * t82 + pkin(1);
t118 = t79 * t73;
t116 = t79 * t82;
t115 = t80 * t81;
t104 = t82 * t110;
t71 = pkin(2) * t104;
t113 = -t54 * pkin(3) + t71;
t108 = t110 * pkin(8);
t107 = -t81 * t103 - t40 * t77;
t105 = t79 * t110;
t102 = t131 * pkin(2);
t41 = -t79 * t57 - t83 * t92;
t58 = t75 * t78 * pkin(2) - t110 * qJ(3) - t108;
t67 = t83 * t73;
t101 = -pkin(3) * t41 - t79 * t58 + t67;
t100 = rSges(5,1) * t81 - rSges(5,2) * t77;
t99 = rSges(6,1) * t80 - rSges(6,2) * t76;
t98 = t110 * rSges(4,3) - t58;
t97 = t39 * pkin(3) + t102;
t62 = -t75 * t116 - t119;
t95 = -t40 * pkin(3) - t83 * t58 - t118;
t55 = t74 * t104 + t78 * t96;
t94 = pkin(9) * t55 - pkin(10) * t121 - t54 * t130 + t113;
t91 = t62 * pkin(2);
t22 = t77 * t105 - t41 * t81;
t90 = t22 * pkin(4) - pkin(9) * t42 + t101;
t88 = t110 * rSges(3,3) + t108;
t87 = t42 * pkin(3) + t91;
t86 = -pkin(4) * t18 + t39 * pkin(9) + t95;
t85 = pkin(9) * t40 + pkin(10) * t124 + t39 * t130 + t97;
t84 = -t41 * pkin(9) + pkin(10) * t122 + t42 * t130 + t87;
t63 = -t75 * t117 + t114;
t61 = -t75 * t119 - t116;
t46 = t55 * t81 + t75 * t77;
t45 = -t55 * t77 + t75 * t81;
t44 = t45 * pkin(4);
t24 = -t54 * t115 + t55 * t76;
t23 = -t54 * t120 - t55 * t80;
t21 = -t81 * t105 - t41 * t77;
t15 = t21 * pkin(4);
t13 = t107 * pkin(4);
t12 = t46 * t80 + t54 * t76;
t11 = t46 * t76 - t54 * t80;
t10 = t42 * t115 - t41 * t76;
t9 = t42 * t120 + t41 * t80;
t8 = t39 * t115 + t40 * t76;
t7 = t39 * t120 - t40 * t80;
t6 = t22 * t80 - t42 * t76;
t5 = t22 * t76 + t42 * t80;
t3 = [-m(2) * (g(1) * (-t79 * rSges(2,1) - rSges(2,2) * t83) + g(2) * (rSges(2,1) * t83 - t79 * rSges(2,2))) - m(3) * (g(1) * (t61 * rSges(3,1) - rSges(3,2) * t131 - t79 * pkin(1) + t88 * t83) + g(2) * (t63 * rSges(3,1) + t62 * rSges(3,2) + t83 * pkin(1) + t88 * t79)) - m(4) * (g(1) * (-t40 * rSges(4,1) - t39 * rSges(4,2) + t98 * t83 - t118) + g(2) * (-rSges(4,1) * t41 + t42 * rSges(4,2) + t98 * t79 + t67)) - m(5) * (g(1) * (-rSges(5,1) * t18 - rSges(5,2) * t107 + t127 * t39 + t95) + g(2) * (rSges(5,1) * t22 - rSges(5,2) * t21 - t127 * t42 + t101)) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t107 * t126 + t86) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t126 * t21 + t90)) - m(7) * (g(1) * (-t1 * t112 + t107 * t128 - t129 * t2 + t86) + g(2) * (t112 * t5 + t128 * t21 + t129 * t6 + t90)) -m(3) * (g(1) * (rSges(3,1) * t62 - rSges(3,2) * t63) + g(2) * (rSges(3,1) * t131 + rSges(3,2) * t61) + g(3) * (rSges(3,1) * t104 - rSges(3,2) * t106)) - m(4) * (g(1) * (t42 * rSges(4,1) + t41 * rSges(4,2) + t91) + g(2) * (rSges(4,1) * t39 - rSges(4,2) * t40 + t102) + g(3) * (-rSges(4,1) * t54 - rSges(4,2) * t55 + t71)) - m(5) * (g(1) * (t100 * t42 - t127 * t41 + t87) + g(2) * (t100 * t39 + t127 * t40 + t97) + g(3) * (-t100 * t54 + t127 * t55 + t113)) - m(6) * (g(1) * (t10 * rSges(6,1) - t9 * rSges(6,2) + rSges(6,3) * t122 + t84) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 + rSges(6,3) * t124 + t85) + g(3) * (rSges(6,1) * t24 - rSges(6,2) * t23 - rSges(6,3) * t121 + t94)) - m(7) * (g(1) * (rSges(7,2) * t122 + t129 * t10 + t112 * t9 + t84) + g(2) * (rSges(7,2) * t124 + t112 * t7 + t129 * t8 + t85) + g(3) * (-rSges(7,2) * t121 + t112 * t23 + t129 * t24 + t94)) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t105 - g(2) * t103 + g(3) * t75) -m(5) * (g(1) * (-rSges(5,1) * t21 - rSges(5,2) * t22) + g(2) * (rSges(5,1) * t107 - rSges(5,2) * t18) + g(3) * (rSges(5,1) * t45 - rSges(5,2) * t46)) - m(6) * (g(1) * (t126 * t22 - t99 * t21 - t15) + g(2) * (t107 * t99 + t126 * t18 + t13) + g(3) * (t126 * t46 + t99 * t45 + t44)) - m(7) * ((t128 * t46 + t89 * t45 + t44) * g(3) + (t89 * t107 + t128 * t18 + t13) * g(2) + (t128 * t22 - t89 * t21 - t15) * g(1)) -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t11 - rSges(6,2) * t12)) - m(7) * (g(1) * (t112 * t6 - t129 * t5) + g(2) * (-t129 * t1 + t112 * t2) + g(3) * (-t129 * t11 + t112 * t12)) -m(7) * (g(1) * t5 + g(2) * t1 + g(3) * t11)];
taug  = t3(:);
