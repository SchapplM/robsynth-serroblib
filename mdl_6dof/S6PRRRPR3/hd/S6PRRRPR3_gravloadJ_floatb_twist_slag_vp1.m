% Calculate Gravitation load on the joints for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:11:03
% EndTime: 2019-03-08 23:11:05
% DurationCPUTime: 0.89s
% Computational Cost: add. (647->153), mult. (1058->218), div. (0->0), fcn. (1230->12), ass. (0->71)
t122 = rSges(7,3) + pkin(10);
t60 = sin(qJ(6));
t63 = cos(qJ(6));
t118 = rSges(7,1) * t60 + rSges(7,2) * t63;
t59 = sin(pkin(6));
t62 = sin(qJ(2));
t100 = t59 * t62;
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t92 = cos(pkin(6));
t121 = -t61 * t100 + t92 * t64;
t57 = qJ(3) + qJ(4);
t55 = sin(t57);
t56 = cos(t57);
t120 = pkin(4) * t56 + qJ(5) * t55;
t65 = cos(qJ(2));
t58 = sin(pkin(11));
t88 = t58 * t92;
t91 = cos(pkin(11));
t43 = -t62 * t88 + t91 * t65;
t99 = t59 * t64;
t119 = -t43 * t61 + t58 * t99;
t75 = t92 * t91;
t41 = t58 * t65 + t62 * t75;
t87 = t59 * t91;
t18 = t41 * t55 + t56 * t87;
t19 = t41 * t56 - t55 * t87;
t117 = t18 * rSges(6,2) + t19 * rSges(6,3);
t116 = t122 * t56;
t115 = t118 * t19 - t122 * t18;
t114 = -m(6) - m(7);
t40 = t58 * t62 - t65 * t75;
t112 = g(2) * t40;
t111 = g(3) * t59;
t109 = rSges(4,3) + pkin(8);
t107 = -t18 * rSges(5,1) - t19 * rSges(5,2);
t103 = t55 * t60;
t102 = t55 * t63;
t101 = t58 * t59;
t98 = t59 * t65;
t20 = -t56 * t101 + t43 * t55;
t21 = t55 * t101 + t43 * t56;
t97 = -t20 * rSges(5,1) - t21 * rSges(5,2);
t54 = pkin(3) * t64 + pkin(2);
t66 = -pkin(9) - pkin(8);
t96 = -t40 * t54 - t41 * t66;
t42 = t91 * t62 + t65 * t88;
t95 = -t42 * t54 - t43 * t66;
t36 = t55 * t100 - t92 * t56;
t37 = t56 * t100 + t92 * t55;
t94 = -t36 * rSges(5,1) - t37 * rSges(5,2);
t85 = -t18 * pkin(4) + t19 * qJ(5);
t84 = -t20 * pkin(4) + qJ(5) * t21;
t83 = -t36 * pkin(4) + qJ(5) * t37;
t82 = -t120 * t40 + t96;
t81 = -t120 * t42 + t95;
t45 = t54 * t98;
t80 = g(3) * (t120 * t98 + t45);
t79 = t119 * pkin(3);
t78 = rSges(5,1) * t56 - rSges(5,2) * t55;
t77 = rSges(6,2) * t56 - rSges(6,3) * t55;
t76 = t121 * pkin(3);
t74 = t20 * rSges(6,2) + t21 * rSges(6,3) + t84;
t73 = t36 * rSges(6,2) + t37 * rSges(6,3) + t83;
t72 = rSges(4,1) * t64 - rSges(4,2) * t61 + pkin(2);
t71 = -t41 * t61 - t64 * t87;
t70 = t118 * t21 - t122 * t20 + t84;
t69 = t118 * t37 - t122 * t36 + t83;
t68 = t71 * pkin(3);
t67 = t68 + t85;
t1 = [(-m(2) - m(3) - m(4) - m(5) + t114) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t42 - rSges(3,2) * t43) + g(2) * (-rSges(3,1) * t40 - rSges(3,2) * t41) + (rSges(3,1) * t65 - rSges(3,2) * t62) * t111) - m(4) * (g(1) * (t109 * t43 - t72 * t42) + g(2) * t109 * t41 - t72 * t112 + (t109 * t62 + t72 * t65) * t111) - m(5) * (g(1) * (rSges(5,3) * t43 - t78 * t42 + t95) + g(2) * (rSges(5,3) * t41 - t78 * t40 + t96) + g(3) * t45 + (t78 * t65 + (rSges(5,3) - t66) * t62) * t111) - m(6) * (g(1) * (rSges(6,1) * t43 + t77 * t42 + t81) + g(2) * (rSges(6,1) * t41 + t77 * t40 + t82) + t80 + (-t77 * t65 + (rSges(6,1) - t66) * t62) * t111) - m(7) * (g(1) * (t43 * pkin(5) + (-t42 * t103 + t43 * t63) * rSges(7,1) + (-t42 * t102 - t43 * t60) * rSges(7,2) + t81) + g(2) * (t41 * pkin(5) + (-t40 * t103 + t41 * t63) * rSges(7,1) + (-t40 * t102 - t41 * t60) * rSges(7,2) + t82) + t80 + (-g(1) * t42 - t112) * t116 + ((t63 * rSges(7,1) - t60 * rSges(7,2) + pkin(5) - t66) * t62 + (t118 * t55 + t116) * t65) * t111) -m(4) * (g(1) * (t119 * rSges(4,1) + (-t61 * t101 - t43 * t64) * rSges(4,2)) + g(2) * (t71 * rSges(4,1) + (-t41 * t64 + t61 * t87) * rSges(4,2)) + g(3) * (t121 * rSges(4,1) + (-t92 * t61 - t62 * t99) * rSges(4,2))) - m(5) * (g(1) * (t79 + t97) + g(2) * (t68 + t107) + g(3) * (t76 + t94)) - m(6) * (g(1) * (t74 + t79) + g(2) * (t67 + t117) + g(3) * (t73 + t76)) - m(7) * (g(1) * (t70 + t79) + g(2) * (t67 + t115) + g(3) * (t69 + t76)) -m(5) * (g(1) * t97 + g(2) * t107 + g(3) * t94) - m(6) * (g(1) * t74 + g(2) * (t85 + t117) + g(3) * t73) - m(7) * (g(1) * t70 + g(2) * (t85 + t115) + g(3) * t69) t114 * (g(1) * t20 + g(2) * t18 + g(3) * t36) -m(7) * (g(1) * ((t20 * t63 - t42 * t60) * rSges(7,1) + (-t20 * t60 - t42 * t63) * rSges(7,2)) + g(2) * ((t18 * t63 - t40 * t60) * rSges(7,1) + (-t18 * t60 - t40 * t63) * rSges(7,2)) + g(3) * ((t36 * t63 + t60 * t98) * rSges(7,1) + (-t36 * t60 + t63 * t98) * rSges(7,2)))];
taug  = t1(:);
