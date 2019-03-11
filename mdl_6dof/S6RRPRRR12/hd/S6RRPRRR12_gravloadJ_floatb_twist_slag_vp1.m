% Calculate Gravitation load on the joints for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:16
% EndTime: 2019-03-09 14:35:19
% DurationCPUTime: 1.00s
% Computational Cost: add. (666->182), mult. (1209->249), div. (0->0), fcn. (1403->12), ass. (0->77)
t129 = rSges(7,3) + pkin(11);
t67 = sin(pkin(6));
t75 = cos(qJ(1));
t113 = t67 * t75;
t70 = sin(qJ(2));
t71 = sin(qJ(1));
t74 = cos(qJ(2));
t107 = cos(pkin(6));
t98 = t75 * t107;
t46 = t70 * t71 - t74 * t98;
t69 = sin(qJ(4));
t73 = cos(qJ(4));
t89 = t69 * t113 + t46 * t73;
t115 = t67 * t71;
t99 = t71 * t107;
t48 = t75 * t70 + t74 * t99;
t20 = -t69 * t115 + t48 * t73;
t68 = sin(qJ(6));
t72 = cos(qJ(6));
t139 = rSges(7,1) * t72 - rSges(7,2) * t68 + pkin(5);
t128 = -pkin(9) - rSges(5,3);
t138 = -rSges(7,1) * t68 - rSges(7,2) * t72;
t47 = t70 * t98 + t71 * t74;
t49 = -t70 * t99 + t74 * t75;
t137 = g(1) * t49 + g(2) * t47;
t66 = qJ(4) + qJ(5);
t63 = sin(t66);
t64 = cos(t66);
t91 = t64 * t113 - t46 * t63;
t136 = t47 * t72 + t68 * t91;
t135 = -t47 * t68 + t72 * t91;
t14 = t63 * t115 - t48 * t64;
t15 = t64 * t115 + t48 * t63;
t134 = -t14 * rSges(6,1) - t15 * rSges(6,2);
t133 = pkin(4) * t69;
t130 = g(3) * t67;
t90 = t63 * t113 + t46 * t64;
t127 = rSges(6,1) * t90 + rSges(6,2) * t91;
t122 = t46 * t69;
t118 = t48 * t69;
t116 = t67 * t70;
t114 = t67 * t74;
t30 = -t107 * t63 - t64 * t114;
t31 = t107 * t64 - t63 * t114;
t112 = t30 * rSges(6,1) - t31 * rSges(6,2);
t111 = t89 * pkin(4);
t110 = pkin(2) * t114 + qJ(3) * t116;
t109 = t75 * pkin(1) + pkin(8) * t115;
t108 = rSges(4,3) + qJ(3);
t42 = t46 * pkin(2);
t76 = -pkin(10) - pkin(9);
t104 = t47 * t133 + t46 * t76 - t42;
t44 = t48 * pkin(2);
t103 = t49 * t133 + t48 * t76 - t44;
t102 = t49 * pkin(2) + t109;
t101 = -t71 * pkin(1) + pkin(8) * t113;
t100 = t73 * t113 - t122;
t97 = g(3) * (t116 * t133 + t110);
t96 = -t47 * pkin(2) + t101;
t95 = rSges(5,1) * t69 + rSges(5,2) * t73;
t94 = rSges(6,1) * t63 + rSges(6,2) * t64;
t93 = t20 * pkin(4);
t92 = qJ(3) * t48 + t102;
t86 = -t46 * qJ(3) + t96;
t85 = -t107 * t69 - t73 * t114;
t84 = t129 * t15 - t139 * t14;
t62 = pkin(4) * t73 + pkin(3);
t83 = pkin(4) * t118 + t62 * t115 - t49 * t76 + t92;
t82 = -t129 * t91 + t139 * t90;
t81 = t129 * t31 + t139 * t30;
t80 = t85 * pkin(4);
t79 = -pkin(4) * t122 + t62 * t113 + t47 * t76 + t86;
t78 = -t129 * t64 + t139 * t63;
t21 = t73 * t115 + t118;
t2 = t15 * t72 + t49 * t68;
t1 = -t15 * t68 + t49 * t72;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t71 - rSges(2,2) * t75) + g(2) * (rSges(2,1) * t75 - rSges(2,2) * t71)) - m(3) * (g(1) * (-rSges(3,1) * t47 + rSges(3,2) * t46 + rSges(3,3) * t113 + t101) + g(2) * (rSges(3,1) * t49 - rSges(3,2) * t48 + rSges(3,3) * t115 + t109)) - m(4) * (g(1) * (rSges(4,1) * t113 + rSges(4,2) * t47 - t108 * t46 + t96) + g(2) * (rSges(4,1) * t115 - rSges(4,2) * t49 + t108 * t48 + t102)) - m(5) * (g(1) * (t100 * rSges(5,1) - t89 * rSges(5,2) + pkin(3) * t113 + t128 * t47 + t86) + g(2) * (rSges(5,1) * t21 + rSges(5,2) * t20 + pkin(3) * t115 - t128 * t49 + t92)) - m(6) * (g(1) * (rSges(6,1) * t91 - rSges(6,2) * t90 - rSges(6,3) * t47 + t79) + g(2) * (rSges(6,1) * t15 - rSges(6,2) * t14 + rSges(6,3) * t49 + t83)) - m(7) * (g(1) * (t135 * rSges(7,1) - t136 * rSges(7,2) + t91 * pkin(5) + t129 * t90 + t79) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t15 + t129 * t14 + t83)) -m(3) * (g(1) * (-rSges(3,1) * t48 - rSges(3,2) * t49) + g(2) * (-rSges(3,1) * t46 - rSges(3,2) * t47) + (rSges(3,1) * t74 - rSges(3,2) * t70) * t130) - m(4) * (g(1) * (rSges(4,2) * t48 + t108 * t49 - t44) + g(2) * (rSges(4,2) * t46 + t108 * t47 - t42) + g(3) * ((-rSges(4,2) * t74 + rSges(4,3) * t70) * t67 + t110)) - m(5) * (g(1) * (t128 * t48 - t44) + g(2) * (t128 * t46 - t42) + g(3) * t110 + (-t128 * t74 + t95 * t70) * t130 + t137 * (qJ(3) + t95)) - m(6) * (g(1) * (-rSges(6,3) * t48 + t103) + g(2) * (-rSges(6,3) * t46 + t104) + t97 + ((rSges(6,3) - t76) * t74 + t94 * t70) * t130 + t137 * (qJ(3) + t94)) - m(7) * (g(1) * (t138 * t48 + t103) + g(2) * (t138 * t46 + t104) + t97 + ((-t76 - t138) * t74 + t78 * t70) * t130 + t137 * (qJ(3) + t78)) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t48 + g(2) * t46 - g(3) * t114) -m(5) * (g(1) * (rSges(5,1) * t20 - rSges(5,2) * t21) + g(2) * (t89 * rSges(5,1) + t100 * rSges(5,2)) + g(3) * (t85 * rSges(5,1) + (-t107 * t73 + t69 * t114) * rSges(5,2))) - m(6) * (g(1) * (t93 + t134) + g(2) * (t111 + t127) + g(3) * (t80 + t112)) - m(7) * (g(1) * (t84 + t93) + g(2) * (t82 + t111) + g(3) * (t80 + t81)) -m(6) * (g(1) * t134 + g(2) * t127 + g(3) * t112) - m(7) * (g(1) * t84 + g(2) * t82 + g(3) * t81) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (t136 * rSges(7,1) + t135 * rSges(7,2)) + g(3) * ((t72 * t116 - t31 * t68) * rSges(7,1) + (-t68 * t116 - t31 * t72) * rSges(7,2)))];
taug  = t3(:);
