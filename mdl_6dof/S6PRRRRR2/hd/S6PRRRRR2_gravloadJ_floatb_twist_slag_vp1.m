% Calculate Gravitation load on the joints for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:27
% EndTime: 2019-03-09 00:42:29
% DurationCPUTime: 0.90s
% Computational Cost: add. (785->150), mult. (1226->219), div. (0->0), fcn. (1437->14), ass. (0->73)
t134 = rSges(6,3) + pkin(10);
t133 = rSges(7,3) + pkin(11) + pkin(10);
t100 = cos(pkin(6));
t63 = sin(pkin(6));
t66 = sin(qJ(2));
t111 = t63 * t66;
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t132 = t100 * t68 - t65 * t111;
t110 = t63 * t68;
t69 = cos(qJ(2));
t62 = sin(pkin(12));
t96 = t62 * t100;
t99 = cos(pkin(12));
t47 = -t66 * t96 + t99 * t69;
t131 = t62 * t110 - t47 * t65;
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t130 = rSges(6,1) * t67 - rSges(6,2) * t64 + pkin(4);
t60 = qJ(5) + qJ(6);
t56 = sin(t60);
t58 = cos(t60);
t129 = rSges(7,1) * t58 - rSges(7,2) * t56 + pkin(5) * t67 + pkin(4);
t128 = t64 * rSges(6,1) + t67 * rSges(6,2);
t88 = t100 * t99;
t44 = t62 * t66 - t69 * t88;
t123 = g(2) * t44;
t46 = t99 * t66 + t69 * t96;
t127 = g(1) * t46 + t123;
t61 = qJ(3) + qJ(4);
t57 = sin(t61);
t59 = cos(t61);
t126 = t129 * t59 + t133 * t57;
t125 = t130 * t59 + t134 * t57;
t45 = t62 * t69 + t66 * t88;
t122 = g(2) * t45;
t109 = t63 * t69;
t55 = pkin(3) * t68 + pkin(2);
t121 = g(3) * t55 * t109;
t120 = g(3) * t63;
t119 = rSges(4,3) + pkin(8);
t112 = t62 * t63;
t95 = t63 * t99;
t25 = -t45 * t57 - t59 * t95;
t26 = t45 * t59 - t57 * t95;
t105 = t25 * rSges(5,1) - t26 * rSges(5,2);
t27 = t59 * t112 - t47 * t57;
t28 = t57 * t112 + t47 * t59;
t104 = t27 * rSges(5,1) - t28 * rSges(5,2);
t71 = -pkin(9) - pkin(8);
t103 = -t44 * t55 - t45 * t71;
t102 = -t46 * t55 - t47 * t71;
t40 = t100 * t59 - t57 * t111;
t41 = t100 * t57 + t59 * t111;
t101 = t40 * rSges(5,1) - t41 * rSges(5,2);
t93 = t131 * pkin(3);
t92 = rSges(5,1) * t59 - rSges(5,2) * t57;
t91 = -t26 * t64 + t44 * t67;
t90 = -t28 * t64 + t46 * t67;
t89 = t132 * pkin(3);
t87 = rSges(4,1) * t68 - rSges(4,2) * t65 + pkin(2);
t85 = -t67 * t109 - t41 * t64;
t83 = t129 * t25 + t133 * t26;
t82 = t129 * t27 + t133 * t28;
t81 = t129 * t40 + t133 * t41;
t80 = t56 * rSges(7,1) + t58 * rSges(7,2) + t64 * pkin(5);
t79 = -t45 * t65 - t68 * t95;
t78 = t130 * t25 + t134 * t26;
t77 = t130 * t27 + t134 * t28;
t76 = t130 * t40 + t134 * t41;
t75 = t79 * pkin(3);
t72 = m(7) * (g(1) * ((-t28 * t56 + t46 * t58) * rSges(7,1) + (-t28 * t58 - t46 * t56) * rSges(7,2)) + g(2) * ((-t26 * t56 + t44 * t58) * rSges(7,1) + (-t26 * t58 - t44 * t56) * rSges(7,2)) + g(3) * ((-t58 * t109 - t41 * t56) * rSges(7,1) + (t56 * t109 - t41 * t58) * rSges(7,2)));
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t46 - rSges(3,2) * t47) + g(2) * (-rSges(3,1) * t44 - rSges(3,2) * t45) + (rSges(3,1) * t69 - rSges(3,2) * t66) * t120) - m(4) * (g(1) * (t119 * t47 - t87 * t46) + t119 * t122 - t87 * t123 + (t119 * t66 + t87 * t69) * t120) - m(5) * (g(1) * (rSges(5,3) * t47 - t92 * t46 + t102) + g(2) * (rSges(5,3) * t45 - t92 * t44 + t103) + t121 + (t92 * t69 + (rSges(5,3) - t71) * t66) * t120) - m(6) * (g(1) * (t128 * t47 + t102) + g(2) * (t128 * t45 + t103) + t121 + ((-t71 + t128) * t66 + t125 * t69) * t120 - t127 * t125) - m(7) * (g(2) * t103 + t121 + t80 * t122 + ((-t71 + t80) * t66 + t126 * t69) * t120 - t127 * t126 + (t80 * t47 + t102) * g(1)) -m(4) * (g(1) * (t131 * rSges(4,1) + (-t65 * t112 - t47 * t68) * rSges(4,2)) + g(2) * (t79 * rSges(4,1) + (-t45 * t68 + t65 * t95) * rSges(4,2)) + g(3) * (t132 * rSges(4,1) + (-t100 * t65 - t66 * t110) * rSges(4,2))) - m(5) * (g(1) * (t93 + t104) + g(2) * (t75 + t105) + g(3) * (t89 + t101)) - m(6) * (g(1) * (t77 + t93) + g(2) * (t75 + t78) + g(3) * (t76 + t89)) - m(7) * (g(1) * (t82 + t93) + g(2) * (t75 + t83) + g(3) * (t81 + t89)) -m(5) * (g(1) * t104 + g(2) * t105 + g(3) * t101) - m(6) * (g(1) * t77 + g(2) * t78 + g(3) * t76) - m(7) * (g(1) * t82 + g(2) * t83 + g(3) * t81) -m(6) * (g(1) * (t90 * rSges(6,1) + (-t28 * t67 - t46 * t64) * rSges(6,2)) + g(2) * (t91 * rSges(6,1) + (-t26 * t67 - t44 * t64) * rSges(6,2)) + g(3) * (t85 * rSges(6,1) + (t64 * t109 - t41 * t67) * rSges(6,2))) - t72 - m(7) * (g(1) * t90 + g(2) * t91 + g(3) * t85) * pkin(5), -t72];
taug  = t1(:);
