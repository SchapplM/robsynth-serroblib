% Calculate Gravitation load on the joints for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:09
% EndTime: 2019-03-09 18:44:13
% DurationCPUTime: 1.51s
% Computational Cost: add. (870->197), mult. (1434->277), div. (0->0), fcn. (1686->14), ass. (0->83)
t56 = sin(pkin(6));
t60 = sin(qJ(2));
t103 = t56 * t60;
t59 = sin(qJ(3));
t63 = cos(qJ(3));
t93 = cos(pkin(6));
t133 = -t59 * t103 + t93 * t63;
t101 = t56 * t63;
t111 = cos(qJ(1));
t64 = cos(qJ(2));
t61 = sin(qJ(1));
t85 = t61 * t93;
t33 = t111 * t64 - t60 * t85;
t18 = t61 * t101 - t33 * t59;
t80 = t93 * t111;
t31 = t60 * t80 + t61 * t64;
t54 = qJ(3) + pkin(12);
t49 = sin(t54);
t50 = cos(t54);
t89 = t56 * t111;
t13 = t31 * t50 - t49 * t89;
t30 = t60 * t61 - t64 * t80;
t58 = sin(qJ(5));
t62 = cos(qJ(5));
t132 = t13 * t58 - t30 * t62;
t108 = t30 * t58;
t131 = -t13 * t62 - t108;
t55 = qJ(5) + qJ(6);
t51 = sin(t55);
t52 = cos(t55);
t130 = t13 * t51 - t30 * t52;
t129 = -t13 * t52 - t30 * t51;
t71 = t31 * t59 + t63 * t89;
t128 = pkin(3) * t71;
t112 = pkin(10) + rSges(6,3);
t100 = t56 * t64;
t127 = g(3) * t100;
t94 = pkin(11) + pkin(10) + rSges(7,3);
t126 = t58 * rSges(6,1) + t62 * rSges(6,2);
t118 = g(2) * t30;
t32 = t111 * t60 + t64 * t85;
t125 = g(1) * t32 + t118;
t102 = t56 * t61;
t16 = -t102 * t50 + t33 * t49;
t24 = -t103 * t49 + t50 * t93;
t87 = -t31 * t49 - t50 * t89;
t124 = -g(1) * t16 + g(2) * t87 + g(3) * t24;
t47 = pkin(5) * t62 + pkin(4);
t73 = rSges(7,1) * t52 - rSges(7,2) * t51 + t47;
t123 = t49 * t94 + t50 * t73;
t75 = rSges(6,1) * t62 - rSges(6,2) * t58 + pkin(4);
t122 = t112 * t49 + t50 * t75;
t117 = g(2) * t31;
t48 = pkin(3) * t63 + pkin(2);
t115 = t48 * t127;
t114 = g(3) * t56;
t113 = -pkin(9) - rSges(4,3);
t105 = t32 * t58;
t57 = -qJ(4) - pkin(9);
t97 = -t30 * t48 - t31 * t57;
t96 = -t32 * t48 - t33 * t57;
t95 = t111 * pkin(1) + pkin(8) * t102;
t92 = t59 * t102;
t88 = -t61 * pkin(1) + pkin(8) * t89;
t42 = t59 * t89;
t86 = -t31 * t63 + t42;
t82 = t18 * pkin(3);
t81 = pkin(3) * t92 - t32 * t57 + t33 * t48 + t95;
t79 = rSges(5,1) * t50 - rSges(5,2) * t49;
t17 = t102 * t49 + t33 * t50;
t7 = -t17 * t58 + t32 * t62;
t77 = t133 * pkin(3);
t76 = t63 * rSges(4,1) - t59 * rSges(4,2) + pkin(2);
t25 = t103 * t50 + t49 * t93;
t74 = -t100 * t62 - t25 * t58;
t72 = pkin(3) * t42 + t30 * t57 - t31 * t48 + t88;
t70 = t51 * rSges(7,1) + t52 * rSges(7,2) + t58 * pkin(5);
t5 = -t17 * t51 + t32 * t52;
t6 = t17 * t52 + t32 * t51;
t66 = m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * (-rSges(7,1) * t130 + t129 * rSges(7,2)) + g(3) * ((-t100 * t52 - t25 * t51) * rSges(7,1) + (t100 * t51 - t25 * t52) * rSges(7,2)));
t19 = t33 * t63 + t92;
t8 = t17 * t62 + t105;
t1 = [-m(2) * (g(1) * (-t61 * rSges(2,1) - rSges(2,2) * t111) + g(2) * (rSges(2,1) * t111 - t61 * rSges(2,2))) - m(3) * (g(1) * (-t31 * rSges(3,1) + t30 * rSges(3,2) + rSges(3,3) * t89 + t88) + g(2) * (rSges(3,1) * t33 - rSges(3,2) * t32 + rSges(3,3) * t102 + t95)) - m(4) * (g(1) * (rSges(4,1) * t86 + rSges(4,2) * t71 - t31 * pkin(2) + t113 * t30 + t88) + g(2) * (rSges(4,1) * t19 + rSges(4,2) * t18 + pkin(2) * t33 - t113 * t32 + t95)) - m(5) * (g(1) * (-rSges(5,1) * t13 - rSges(5,2) * t87 - rSges(5,3) * t30 + t72) + g(2) * (rSges(5,1) * t17 - rSges(5,2) * t16 + rSges(5,3) * t32 + t81)) - m(6) * (g(1) * (t131 * rSges(6,1) + t132 * rSges(6,2) - t13 * pkin(4) + t112 * t87 + t72) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 + pkin(4) * t17 + t112 * t16 + t81)) - m(7) * (g(1) * (t129 * rSges(7,1) + rSges(7,2) * t130 - pkin(5) * t108 - t13 * t47 + t94 * t87 + t72) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + pkin(5) * t105 + t16 * t94 + t17 * t47 + t81)) -m(3) * (g(1) * (-rSges(3,1) * t32 - rSges(3,2) * t33) + g(2) * (-rSges(3,1) * t30 - rSges(3,2) * t31) + (rSges(3,1) * t64 - rSges(3,2) * t60) * t114) - m(4) * (g(1) * (-t113 * t33 - t32 * t76) - t113 * t117 - t76 * t118 + (-t113 * t60 + t64 * t76) * t114) - m(5) * (g(1) * (rSges(5,3) * t33 - t32 * t79 + t96) + g(2) * (rSges(5,3) * t31 - t30 * t79 + t97) + t115 + (t79 * t64 + (rSges(5,3) - t57) * t60) * t114) - m(6) * (g(1) * (t126 * t33 + t96) + g(2) * (t126 * t31 + t97) + t115 + ((-t57 + t126) * t60 + t122 * t64) * t114 - t125 * t122) - m(7) * (g(2) * t97 + t115 + t70 * t117 + ((-t57 + t70) * t60 + t123 * t64) * t114 - t125 * t123 + (t70 * t33 + t96) * g(1)) -m(4) * (g(1) * (rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (-rSges(4,1) * t71 + rSges(4,2) * t86) + g(3) * (t133 * rSges(4,1) + (-t101 * t60 - t59 * t93) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t16 - rSges(5,2) * t17 + t82) + g(2) * (rSges(5,1) * t87 - t13 * rSges(5,2) - t128) + g(3) * (rSges(5,1) * t24 - rSges(5,2) * t25 + t77)) + (-g(1) * (t94 * t17 + t82) - g(2) * (t94 * t13 - t128) - g(3) * (t94 * t25 + t77) - t124 * t73) * m(7) + (-g(1) * (t112 * t17 + t82) - g(2) * (t112 * t13 - t128) - g(3) * (t112 * t25 + t77) - t124 * t75) * m(6) (-m(5) - m(6) - m(7)) * (t125 - t127) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (-rSges(6,1) * t132 + t131 * rSges(6,2)) + g(3) * (t74 * rSges(6,1) + (t100 * t58 - t25 * t62) * rSges(6,2))) - t66 - m(7) * (g(1) * t7 - g(2) * t132 + g(3) * t74) * pkin(5), -t66];
taug  = t1(:);
