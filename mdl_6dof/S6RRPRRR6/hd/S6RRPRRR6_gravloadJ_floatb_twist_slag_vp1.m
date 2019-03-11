% Calculate Gravitation load on the joints for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:32
% EndTime: 2019-03-09 13:50:34
% DurationCPUTime: 1.01s
% Computational Cost: add. (528->167), mult. (806->211), div. (0->0), fcn. (847->10), ass. (0->84)
t130 = -pkin(10) - rSges(7,3);
t57 = sin(qJ(6));
t61 = cos(qJ(6));
t132 = t61 * rSges(7,1) - t57 * rSges(7,2) + pkin(5);
t131 = -rSges(5,3) - pkin(8);
t59 = sin(qJ(2));
t63 = cos(qJ(2));
t93 = qJ(4) + qJ(5);
t89 = sin(t93);
t90 = cos(t93);
t27 = t59 * t89 + t63 * t90;
t28 = t59 * t90 - t63 * t89;
t102 = -t27 * rSges(6,1) - t28 * rSges(6,2);
t58 = sin(qJ(4));
t109 = t58 * t59;
t47 = pkin(4) * t109;
t62 = cos(qJ(4));
t49 = pkin(4) * t62 + pkin(3);
t129 = t63 * t49 + t47;
t52 = t59 * qJ(3);
t96 = t63 * pkin(2) + t52;
t60 = sin(qJ(1));
t14 = t27 * t60;
t128 = t130 * t14;
t64 = cos(qJ(1));
t82 = t64 * t89;
t83 = t64 * t90;
t15 = -t59 * t82 - t63 * t83;
t127 = t130 * t15;
t126 = t130 * t28;
t125 = g(1) * t64 + g(2) * t60;
t124 = t125 * t59;
t13 = t28 * t60;
t122 = t13 * rSges(6,1) - t14 * rSges(6,2);
t16 = -t59 * t83 + t63 * t82;
t121 = -t16 * rSges(6,1) + t15 * rSges(6,2);
t120 = pkin(3) * t63;
t119 = pkin(4) * t58;
t112 = -pkin(2) - t49;
t111 = rSges(4,1) * t63;
t108 = t59 * t62;
t107 = t59 * t64;
t106 = t60 * t63;
t104 = t62 * t63;
t103 = t63 * t64;
t37 = t106 * t119;
t101 = t60 * pkin(4) * t108 - t37;
t94 = qJ(3) * t63;
t44 = t60 * t94;
t100 = t37 + t44;
t92 = t58 * t103;
t40 = pkin(4) * t92;
t91 = t62 * t107;
t99 = pkin(4) * t91 - t40;
t46 = t64 * t94;
t98 = t40 + t46;
t55 = t64 * pkin(7);
t65 = -pkin(9) - pkin(8);
t97 = t64 * t65 + t55;
t95 = t64 * pkin(1) + t60 * pkin(7);
t88 = t96 + t129;
t87 = pkin(2) * t103 + t64 * t52 + t95;
t86 = -pkin(4) * t104 - t47;
t81 = rSges(3,1) * t63 - rSges(3,2) * t59;
t74 = t58 * t63 - t108;
t23 = t74 * t60;
t73 = t104 + t109;
t24 = t73 * t60;
t79 = -rSges(5,1) * t23 - rSges(5,2) * t24;
t25 = -t91 + t92;
t26 = t73 * t64;
t78 = -rSges(5,1) * t25 - rSges(5,2) * t26;
t77 = -rSges(5,1) * t73 + rSges(5,2) * t74;
t76 = -t14 * t61 - t57 * t64;
t75 = t14 * t57 - t61 * t64;
t71 = -pkin(1) - t96;
t70 = t49 * t103 + t64 * t47 + t60 * t65 + t87;
t69 = t132 * t13 - t128;
t68 = -t132 * t16 + t127;
t67 = -t132 * t27 - t126;
t66 = t112 * t124;
t2 = -t15 * t61 - t57 * t60;
t1 = t15 * t57 - t60 * t61;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t60 - rSges(2,2) * t64) + g(2) * (rSges(2,1) * t64 - rSges(2,2) * t60)) - m(3) * (g(1) * (rSges(3,3) * t64 + t55) + g(2) * (rSges(3,1) * t103 - rSges(3,2) * t107 + t95) + (g(1) * (-pkin(1) - t81) + g(2) * rSges(3,3)) * t60) - m(4) * (g(1) * (rSges(4,2) * t64 + t55) + g(2) * (rSges(4,1) * t103 + rSges(4,3) * t107 + t87) + (g(1) * (-rSges(4,3) * t59 - t111 + t71) + g(2) * rSges(4,2)) * t60) - m(5) * (g(1) * (-rSges(5,1) * t24 + rSges(5,2) * t23 + t131 * t64 + t55) + g(2) * (rSges(5,1) * t26 - rSges(5,2) * t25 + pkin(3) * t103 + t87) + (g(1) * (t71 - t120) + g(2) * t131) * t60) - m(6) * (g(1) * (-rSges(6,1) * t14 - rSges(6,2) * t13 - rSges(6,3) * t64 + t97) + g(2) * (-rSges(6,1) * t15 - rSges(6,2) * t16 + t70) + (g(1) * (t71 - t129) - g(2) * rSges(6,3)) * t60) - m(7) * (g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 - pkin(5) * t15 - t130 * t16 + t70) + (t76 * rSges(7,1) + t75 * rSges(7,2) - t14 * pkin(5) + t97 + (-pkin(1) + t112 * t63 + (-qJ(3) - t119) * t59) * t60 - t130 * t13) * g(1)) -m(3) * (g(3) * t81 + t125 * (-rSges(3,1) * t59 - rSges(3,2) * t63)) - m(4) * (g(1) * (rSges(4,3) * t103 + t46) + g(2) * (rSges(4,3) * t106 + t44) + g(3) * (t96 + t111) + (g(3) * rSges(4,3) + t125 * (-rSges(4,1) - pkin(2))) * t59) - m(5) * (g(1) * (t46 - t78) + g(2) * (t44 - t79) + g(3) * (-t77 + t96 + t120) + (-pkin(2) - pkin(3)) * t124) - m(6) * (g(1) * (t98 - t121) + g(2) * (t100 - t122) + g(3) * (t88 - t102) + t66) - m(7) * (g(1) * (t98 - t127) + g(2) * (t100 + t128) + g(3) * (t88 + t126) + t66 + (g(1) * t16 - g(2) * t13 + g(3) * t27) * t132) (-m(4) - m(5) - m(6) - m(7)) * (-g(3) * t63 + t124) -m(5) * (g(1) * t78 + g(2) * t79 + g(3) * t77) - m(6) * (g(1) * (t99 + t121) + g(2) * (t101 + t122) + g(3) * (t86 + t102)) - m(7) * (g(1) * (t68 + t99) + g(2) * (t69 + t101) + g(3) * (t67 + t86)) -m(6) * (g(1) * t121 + g(2) * t122 + g(3) * t102) - m(7) * (g(1) * t68 + g(2) * t69 + g(3) * t67) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-rSges(7,1) * t75 + rSges(7,2) * t76) + g(3) * (-t57 * rSges(7,1) - t61 * rSges(7,2)) * t28)];
taug  = t3(:);
