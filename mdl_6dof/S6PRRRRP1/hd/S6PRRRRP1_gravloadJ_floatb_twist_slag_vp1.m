% Calculate Gravitation load on the joints for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:46
% EndTime: 2019-03-08 23:55:49
% DurationCPUTime: 0.99s
% Computational Cost: add. (713->167), mult. (1179->252), div. (0->0), fcn. (1378->12), ass. (0->85)
t132 = rSges(6,3) + pkin(10);
t131 = rSges(7,3) + qJ(6) + pkin(10);
t67 = sin(pkin(6));
t71 = sin(qJ(2));
t108 = t67 * t71;
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t98 = cos(pkin(6));
t130 = -t70 * t108 + t98 * t73;
t74 = cos(qJ(2));
t96 = sin(pkin(11));
t86 = t98 * t96;
t97 = cos(pkin(11));
t53 = -t71 * t86 + t97 * t74;
t93 = t67 * t96;
t129 = -t53 * t70 + t73 * t93;
t87 = t98 * t97;
t50 = t96 * t71 - t74 * t87;
t124 = g(2) * t50;
t52 = t97 * t71 + t74 * t86;
t128 = g(1) * t52 + t124;
t72 = cos(qJ(5));
t62 = pkin(5) * t72 + pkin(4);
t66 = qJ(3) + qJ(4);
t64 = sin(t66);
t65 = cos(t66);
t127 = t131 * t64 + t62 * t65;
t126 = pkin(4) * t65;
t123 = g(3) * t67;
t122 = rSges(4,3) + pkin(8);
t51 = t71 * t87 + t96 * t74;
t94 = t67 * t97;
t27 = t51 * t64 + t65 * t94;
t69 = sin(qJ(5));
t120 = t27 * t69;
t119 = t27 * t72;
t29 = t53 * t64 - t65 * t93;
t118 = t29 * t69;
t117 = t29 * t72;
t46 = t64 * t108 - t98 * t65;
t116 = t46 * t69;
t115 = t46 * t72;
t114 = t51 * t69;
t113 = t53 * t69;
t110 = t65 * t69;
t109 = t65 * t72;
t107 = t67 * t74;
t106 = t69 * t74;
t105 = t72 * t74;
t28 = t51 * t65 - t64 * t94;
t103 = -t27 * rSges(5,1) - t28 * rSges(5,2);
t30 = t53 * t65 + t64 * t93;
t102 = -t29 * rSges(5,1) - t30 * rSges(5,2);
t63 = pkin(3) * t73 + pkin(2);
t75 = -pkin(9) - pkin(8);
t101 = -t50 * t63 - t51 * t75;
t100 = -t52 * t63 - t53 * t75;
t47 = t65 * t108 + t98 * t64;
t99 = -t46 * rSges(5,1) - t47 * rSges(5,2);
t91 = t129 * pkin(3);
t89 = rSges(5,1) * t65 - rSges(5,2) * t64;
t1 = -t28 * t69 + t50 * t72;
t3 = -t30 * t69 + t52 * t72;
t88 = t130 * pkin(3);
t85 = rSges(4,1) * t73 - rSges(4,2) * t70 + pkin(2);
t31 = -t67 * t105 - t47 * t69;
t83 = -rSges(7,1) * t119 + rSges(7,2) * t120 + t131 * t28 - t27 * t62;
t82 = -rSges(7,1) * t117 + rSges(7,2) * t118 + t131 * t30 - t29 * t62;
t81 = -rSges(7,1) * t115 + rSges(7,2) * t116 + t131 * t47 - t46 * t62;
t80 = -t51 * t70 - t73 * t94;
t79 = -rSges(6,1) * t119 + rSges(6,2) * t120 - t27 * pkin(4) + t132 * t28;
t78 = -rSges(6,1) * t117 + rSges(6,2) * t118 - t29 * pkin(4) + t132 * t30;
t77 = -rSges(6,1) * t115 + rSges(6,2) * t116 - t46 * pkin(4) + t132 * t47;
t76 = t80 * pkin(3);
t54 = t63 * t107;
t38 = (t65 * t105 + t69 * t71) * t67;
t37 = (-t65 * t106 + t71 * t72) * t67;
t32 = t67 * t106 - t47 * t72;
t8 = -t52 * t109 + t113;
t7 = t52 * t110 + t53 * t72;
t6 = -t50 * t109 + t114;
t5 = t50 * t110 + t51 * t72;
t4 = -t30 * t72 - t52 * t69;
t2 = -t28 * t72 - t50 * t69;
t9 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t52 - rSges(3,2) * t53) + g(2) * (-rSges(3,1) * t50 - t51 * rSges(3,2)) + (rSges(3,1) * t74 - rSges(3,2) * t71) * t123) - m(4) * (g(1) * (t122 * t53 - t85 * t52) + g(2) * t122 * t51 - t85 * t124 + (t122 * t71 + t85 * t74) * t123) - m(5) * (g(1) * (rSges(5,3) * t53 - t89 * t52 + t100) + g(2) * (rSges(5,3) * t51 - t89 * t50 + t101) + g(3) * t54 + (t89 * t74 + (rSges(5,3) - t75) * t71) * t123) - m(6) * (g(1) * (rSges(6,1) * t8 + rSges(6,2) * t7 - t52 * t126 + t100) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 - t50 * t126 + t101) + g(3) * (t38 * rSges(6,1) + t37 * rSges(6,2) + t107 * t126 - t75 * t108 + t54) + (g(3) * t107 - t128) * t64 * t132) - m(7) * (g(1) * (rSges(7,1) * t8 + rSges(7,2) * t7 + pkin(5) * t113 + t100) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + pkin(5) * t114 + t101) + g(3) * (t38 * rSges(7,1) + t37 * rSges(7,2) + t54) + ((pkin(5) * t69 - t75) * t71 + t127 * t74) * t123 - t128 * t127) -m(4) * (g(1) * (t129 * rSges(4,1) + (-t53 * t73 - t70 * t93) * rSges(4,2)) + g(2) * (t80 * rSges(4,1) + (-t51 * t73 + t70 * t94) * rSges(4,2)) + g(3) * (t130 * rSges(4,1) + (-t73 * t108 - t98 * t70) * rSges(4,2))) - m(5) * (g(1) * (t91 + t102) + g(2) * (t76 + t103) + g(3) * (t88 + t99)) - m(6) * (g(1) * (t78 + t91) + g(2) * (t76 + t79) + g(3) * (t77 + t88)) - m(7) * (g(1) * (t82 + t91) + g(2) * (t76 + t83) + g(3) * (t81 + t88)) -m(5) * (g(1) * t102 + g(2) * t103 + g(3) * t99) - m(6) * (g(1) * t78 + g(2) * t79 + g(3) * t77) - m(7) * (g(1) * t82 + g(2) * t83 + g(3) * t81) -m(6) * (g(1) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (rSges(6,1) * t31 + rSges(6,2) * t32)) + (-g(1) * (rSges(7,1) * t3 + rSges(7,2) * t4) - g(2) * (rSges(7,1) * t1 + rSges(7,2) * t2) - g(3) * (rSges(7,1) * t31 + rSges(7,2) * t32) - (g(1) * t3 + g(2) * t1 + g(3) * t31) * pkin(5)) * m(7), -m(7) * (g(1) * t29 + g(2) * t27 + g(3) * t46)];
taug  = t9(:);
