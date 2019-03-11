% Calculate Gravitation load on the joints for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:56
% EndTime: 2019-03-09 17:04:00
% DurationCPUTime: 1.35s
% Computational Cost: add. (877->204), mult. (1539->292), div. (0->0), fcn. (1833->12), ass. (0->90)
t132 = cos(qJ(1));
t74 = sin(pkin(6));
t106 = t74 * t132;
t78 = sin(qJ(2));
t79 = sin(qJ(1));
t82 = cos(qJ(2));
t112 = cos(pkin(6));
t95 = t112 * t132;
t49 = t78 * t95 + t79 * t82;
t73 = qJ(3) + pkin(11);
t70 = sin(t73);
t71 = cos(t73);
t104 = -t71 * t106 - t49 * t70;
t77 = sin(qJ(3));
t81 = cos(qJ(3));
t85 = t81 * t106 + t49 * t77;
t83 = t85 * pkin(3);
t141 = t104 * pkin(4) - t83;
t122 = t74 * t78;
t140 = t112 * t81 - t77 * t122;
t120 = t74 * t81;
t102 = t79 * t112;
t51 = -t78 * t102 + t132 * t82;
t24 = t79 * t120 - t51 * t77;
t19 = -t70 * t106 + t49 * t71;
t48 = t78 * t79 - t82 * t95;
t76 = sin(qJ(5));
t80 = cos(qJ(5));
t1 = t19 * t76 - t48 * t80;
t2 = t19 * t80 + t48 * t76;
t135 = rSges(7,2) + pkin(10);
t113 = rSges(7,3) + qJ(6);
t136 = rSges(7,1) + pkin(5);
t84 = t113 * t76 + t136 * t80;
t139 = pkin(4) * t71;
t138 = g(2) * t48;
t137 = g(3) * t74;
t134 = rSges(6,3) + pkin(10);
t133 = -pkin(9) - rSges(4,3);
t131 = t48 * t70;
t50 = t82 * t102 + t132 * t78;
t127 = t50 * t70;
t125 = t70 * t82;
t124 = t71 * t76;
t123 = t71 * t80;
t121 = t74 * t79;
t119 = t74 * t82;
t75 = -qJ(4) - pkin(9);
t118 = t75 * t78;
t117 = t80 * t82;
t69 = pkin(3) * t81 + pkin(2);
t116 = -t48 * t69 - t49 * t75;
t115 = -t50 * t69 - t51 * t75;
t114 = t132 * pkin(1) + pkin(8) * t121;
t111 = t77 * t121;
t108 = t76 * t119;
t54 = t69 * t119;
t107 = t54 + (pkin(10) * t70 + t139) * t119;
t105 = -t79 * pkin(1) + pkin(8) * t106;
t64 = t77 * t106;
t103 = -t49 * t81 + t64;
t100 = -pkin(10) * t131 - t48 * t139 + t116;
t99 = -pkin(10) * t127 - t50 * t139 + t115;
t97 = t24 * pkin(3);
t96 = pkin(3) * t111 - t50 * t75 + t51 * t69 + t114;
t94 = rSges(5,1) * t71 - rSges(5,2) * t70;
t93 = rSges(6,1) * t80 - rSges(6,2) * t76;
t92 = t140 * pkin(3);
t22 = -t71 * t121 + t51 * t70;
t91 = -t22 * pkin(4) + t97;
t23 = t70 * t121 + t51 * t71;
t90 = t23 * pkin(4) + t96;
t89 = rSges(4,1) * t81 - rSges(4,2) * t77 + pkin(2);
t37 = t112 * t71 - t70 * t122;
t88 = t37 * pkin(4) + t92;
t87 = pkin(3) * t64 + t48 * t75 - t49 * t69 + t105;
t86 = -pkin(4) * t19 + t87;
t38 = t112 * t70 + t71 * t122;
t27 = (t71 * t117 + t76 * t78) * t74;
t26 = t71 * t108 - t80 * t122;
t25 = t51 * t81 + t111;
t17 = t38 * t80 - t108;
t16 = t74 * t117 + t38 * t76;
t10 = -t50 * t123 + t51 * t76;
t9 = -t50 * t124 - t51 * t80;
t8 = -t48 * t123 + t49 * t76;
t7 = -t48 * t124 - t49 * t80;
t6 = t23 * t80 + t50 * t76;
t5 = t23 * t76 - t50 * t80;
t3 = [-m(2) * (g(1) * (-t79 * rSges(2,1) - t132 * rSges(2,2)) + g(2) * (t132 * rSges(2,1) - t79 * rSges(2,2))) - m(3) * (g(1) * (-t49 * rSges(3,1) + t48 * rSges(3,2) + rSges(3,3) * t106 + t105) + g(2) * (rSges(3,1) * t51 - rSges(3,2) * t50 + rSges(3,3) * t121 + t114)) - m(4) * (g(1) * (t103 * rSges(4,1) + t85 * rSges(4,2) - t49 * pkin(2) + t133 * t48 + t105) + g(2) * (rSges(4,1) * t25 + rSges(4,2) * t24 + pkin(2) * t51 - t133 * t50 + t114)) - m(5) * (g(1) * (-rSges(5,1) * t19 - rSges(5,2) * t104 - rSges(5,3) * t48 + t87) + g(2) * (rSges(5,1) * t23 - rSges(5,2) * t22 + rSges(5,3) * t50 + t96)) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t104 * t134 + t86) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t134 * t22 + t90)) - m(7) * (g(1) * (-t1 * t113 + t104 * t135 - t136 * t2 + t86) + g(2) * (t113 * t5 + t135 * t22 + t136 * t6 + t90)) -m(3) * (g(1) * (-rSges(3,1) * t50 - rSges(3,2) * t51) + g(2) * (-rSges(3,1) * t48 - rSges(3,2) * t49) + (rSges(3,1) * t82 - rSges(3,2) * t78) * t137) - m(4) * (g(1) * (-t133 * t51 - t89 * t50) - g(2) * t133 * t49 - t89 * t138 + (-t133 * t78 + t89 * t82) * t137) - m(5) * (g(1) * (rSges(5,3) * t51 - t94 * t50 + t115) + g(2) * (rSges(5,3) * t49 - t94 * t48 + t116) + g(3) * t54 + (t94 * t82 + (rSges(5,3) - t75) * t78) * t137) - m(6) * (g(1) * (rSges(6,1) * t10 - rSges(6,2) * t9 - rSges(6,3) * t127 + t99) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t131 + t100) + g(3) * (t27 * rSges(6,1) - t26 * rSges(6,2) + (rSges(6,3) * t125 - t118) * t74 + t107)) - m(7) * (g(1) * (-rSges(7,2) * t127 + t136 * t10 + t113 * t9 + t99) + g(2) * (-rSges(7,2) * t131 + t113 * t7 + t136 * t8 + t100) + g(3) * ((rSges(7,2) * t125 - t118) * t74 + t136 * t27 + t113 * t26 + t107)) -m(4) * (g(1) * (rSges(4,1) * t24 - rSges(4,2) * t25) + g(2) * (-t85 * rSges(4,1) + t103 * rSges(4,2)) + g(3) * (t140 * rSges(4,1) + (-t112 * t77 - t78 * t120) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t22 - rSges(5,2) * t23 + t97) + g(2) * (rSges(5,1) * t104 - t19 * rSges(5,2) - t83) + g(3) * (rSges(5,1) * t37 - rSges(5,2) * t38 + t92)) - m(6) * (g(1) * (t134 * t23 - t93 * t22 + t91) + g(2) * (t104 * t93 + t134 * t19 + t141) + g(3) * (t134 * t38 + t93 * t37 + t88)) - m(7) * ((t135 * t38 + t84 * t37 + t88) * g(3) + (t84 * t104 + t135 * t19 + t141) * g(2) + (t135 * t23 - t84 * t22 + t91) * g(1)) (-m(5) - m(6) - m(7)) * (g(1) * t50 - g(3) * t119 + t138) -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t16 - rSges(6,2) * t17)) - m(7) * (g(1) * (t113 * t6 - t136 * t5) + g(2) * (-t136 * t1 + t113 * t2) + g(3) * (t113 * t17 - t136 * t16)) -m(7) * (g(1) * t5 + g(2) * t1 + g(3) * t16)];
taug  = t3(:);
