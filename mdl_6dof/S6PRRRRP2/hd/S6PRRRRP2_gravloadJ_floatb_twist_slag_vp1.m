% Calculate Gravitation load on the joints for
% S6PRRRRP2
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:11
% EndTime: 2019-03-09 00:01:14
% DurationCPUTime: 0.91s
% Computational Cost: add. (796->170), mult. (1334->250), div. (0->0), fcn. (1583->12), ass. (0->92)
t141 = rSges(7,1) + pkin(5);
t114 = cos(pkin(6));
t80 = sin(pkin(6));
t83 = sin(qJ(2));
t126 = t80 * t83;
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t147 = t114 * t85 - t82 * t126;
t125 = t80 * t85;
t79 = sin(pkin(11));
t104 = t79 * t114;
t113 = cos(pkin(11));
t86 = cos(qJ(2));
t64 = -t83 * t104 + t113 * t86;
t146 = t79 * t125 - t64 * t82;
t103 = t80 * t113;
t96 = t114 * t113;
t62 = t79 * t86 + t83 * t96;
t78 = qJ(3) + qJ(4);
t76 = sin(t78);
t77 = cos(t78);
t29 = -t77 * t103 - t62 * t76;
t84 = cos(qJ(5));
t138 = t29 * t84;
t81 = sin(qJ(5));
t139 = t29 * t81;
t30 = -t76 * t103 + t62 * t77;
t145 = rSges(6,1) * t138 - rSges(6,2) * t139 + t30 * rSges(6,3);
t116 = qJ(6) * t81;
t144 = t30 * rSges(7,2) + rSges(7,3) * t139 + t29 * t116 + t141 * t138;
t143 = pkin(4) * t77;
t142 = g(3) * t80;
t140 = rSges(4,3) + pkin(8);
t127 = t79 * t80;
t31 = t77 * t127 - t64 * t76;
t137 = t31 * t81;
t136 = t31 * t84;
t53 = t114 * t77 - t76 * t126;
t135 = t53 * t81;
t134 = t53 * t84;
t61 = t79 * t83 - t86 * t96;
t133 = t61 * t76;
t63 = t86 * t104 + t113 * t83;
t132 = t63 * t76;
t130 = t76 * t86;
t129 = t77 * t81;
t128 = t77 * t84;
t124 = t80 * t86;
t87 = -pkin(9) - pkin(8);
t123 = t83 * t87;
t122 = t84 * t86;
t121 = t29 * rSges(5,1) - t30 * rSges(5,2);
t32 = t76 * t127 + t64 * t77;
t120 = t31 * rSges(5,1) - t32 * rSges(5,2);
t75 = pkin(3) * t85 + pkin(2);
t119 = -t61 * t75 - t62 * t87;
t118 = -t63 * t75 - t64 * t87;
t54 = t114 * t76 + t77 * t126;
t117 = t53 * rSges(5,1) - t54 * rSges(5,2);
t115 = rSges(7,3) + qJ(6);
t109 = t81 * t124;
t65 = t75 * t124;
t108 = t65 + (pkin(10) * t76 + t143) * t124;
t107 = t29 * pkin(4) + t30 * pkin(10);
t106 = t31 * pkin(4) + pkin(10) * t32;
t105 = t53 * pkin(4) + pkin(10) * t54;
t101 = -pkin(10) * t133 - t61 * t143 + t119;
t100 = -pkin(10) * t132 - t63 * t143 + t118;
t99 = t146 * pkin(3);
t98 = rSges(5,1) * t77 - rSges(5,2) * t76;
t97 = t147 * pkin(3);
t95 = rSges(4,1) * t85 - rSges(4,2) * t82 + pkin(2);
t94 = -t85 * t103 - t62 * t82;
t93 = t32 * rSges(7,2) + rSges(7,3) * t137 + t31 * t116 + t141 * t136 + t106;
t92 = t54 * rSges(7,2) + rSges(7,3) * t135 + t53 * t116 + t141 * t134 + t105;
t91 = rSges(6,1) * t136 - rSges(6,2) * t137 + t32 * rSges(6,3) + t106;
t90 = rSges(6,1) * t134 - rSges(6,2) * t135 + t54 * rSges(6,3) + t105;
t89 = t94 * pkin(3);
t88 = t107 + t89;
t41 = (t77 * t122 + t81 * t83) * t80;
t40 = t77 * t109 - t84 * t126;
t34 = t54 * t84 - t109;
t33 = t80 * t122 + t54 * t81;
t8 = -t63 * t128 + t64 * t81;
t7 = -t63 * t129 - t64 * t84;
t6 = -t61 * t128 + t62 * t81;
t5 = -t61 * t129 - t62 * t84;
t4 = t32 * t84 + t63 * t81;
t3 = t32 * t81 - t63 * t84;
t2 = t30 * t84 + t61 * t81;
t1 = t30 * t81 - t61 * t84;
t9 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t63 - rSges(3,2) * t64) + g(2) * (-rSges(3,1) * t61 - rSges(3,2) * t62) + (rSges(3,1) * t86 - rSges(3,2) * t83) * t142) - m(4) * (g(1) * (t140 * t64 - t95 * t63) + (t140 * t83 + t95 * t86) * t142 + (t140 * t62 - t95 * t61) * g(2)) - m(5) * (g(1) * (t64 * rSges(5,3) - t98 * t63 + t118) + g(2) * (t62 * rSges(5,3) - t98 * t61 + t119) + g(3) * t65 + (t98 * t86 + (rSges(5,3) - t87) * t83) * t142) - m(6) * (g(1) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t132 + t100) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 - rSges(6,3) * t133 + t101) + g(3) * (t41 * rSges(6,1) - t40 * rSges(6,2) + (rSges(6,3) * t130 - t123) * t80 + t108)) - m(7) * (g(1) * (-rSges(7,2) * t132 + t115 * t7 + t141 * t8 + t100) + g(2) * (-rSges(7,2) * t133 + t115 * t5 + t141 * t6 + t101) + g(3) * ((rSges(7,2) * t130 - t123) * t80 + t141 * t41 + t115 * t40 + t108)) -m(4) * (g(1) * (t146 * rSges(4,1) + (-t82 * t127 - t64 * t85) * rSges(4,2)) + g(2) * (t94 * rSges(4,1) + (t82 * t103 - t62 * t85) * rSges(4,2)) + g(3) * (t147 * rSges(4,1) + (-t114 * t82 - t83 * t125) * rSges(4,2))) - m(5) * (g(1) * (t99 + t120) + g(2) * (t89 + t121) + g(3) * (t97 + t117)) - m(6) * (g(1) * (t91 + t99) + g(2) * (t88 + t145) + g(3) * (t90 + t97)) - m(7) * (g(1) * (t93 + t99) + g(2) * (t88 + t144) + g(3) * (t92 + t97)) -m(5) * (g(1) * t120 + g(2) * t121 + g(3) * t117) - m(6) * (g(1) * t91 + g(2) * (t107 + t145) + g(3) * t90) - m(7) * (g(1) * t93 + g(2) * (t107 + t144) + g(3) * t92) -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t33 - rSges(6,2) * t34)) - m(7) * (g(1) * (t115 * t4 - t141 * t3) + g(2) * (-t141 * t1 + t115 * t2) + g(3) * (t115 * t34 - t141 * t33)) -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t33)];
taug  = t9(:);
