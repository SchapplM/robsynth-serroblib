% Calculate Gravitation load on the joints for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:31
% EndTime: 2019-03-10 01:34:35
% DurationCPUTime: 1.54s
% Computational Cost: add. (955->221), mult. (1597->322), div. (0->0), fcn. (1882->12), ass. (0->96)
t148 = rSges(6,3) + pkin(11);
t124 = rSges(7,3) + qJ(6) + pkin(11);
t117 = cos(pkin(6));
t83 = sin(pkin(6));
t87 = sin(qJ(2));
t130 = t83 * t87;
t86 = sin(qJ(3));
t90 = cos(qJ(3));
t156 = t117 * t90 - t86 * t130;
t128 = t83 * t90;
t88 = sin(qJ(1));
t110 = t88 * t117;
t146 = cos(qJ(1));
t91 = cos(qJ(2));
t63 = -t110 * t87 + t146 * t91;
t37 = t88 * t128 - t63 * t86;
t113 = t83 * t146;
t106 = t117 * t146;
t61 = t106 * t87 + t88 * t91;
t82 = qJ(3) + qJ(4);
t79 = sin(t82);
t80 = cos(t82);
t32 = -t79 * t113 + t61 * t80;
t60 = -t106 * t91 + t87 * t88;
t85 = sin(qJ(5));
t89 = cos(qJ(5));
t3 = t32 * t85 - t60 * t89;
t139 = t60 * t85;
t155 = -t32 * t89 - t139;
t150 = g(2) * t60;
t62 = t110 * t91 + t146 * t87;
t154 = g(1) * t62 + t150;
t77 = pkin(5) * t89 + pkin(4);
t153 = t124 * t79 + t77 * t80;
t152 = pkin(4) * t80;
t149 = g(3) * t83;
t147 = -pkin(9) - rSges(4,3);
t31 = t113 * t80 + t61 * t79;
t145 = t31 * t85;
t144 = t31 * t89;
t129 = t83 * t88;
t35 = -t129 * t80 + t63 * t79;
t143 = t35 * t85;
t142 = t35 * t89;
t54 = -t117 * t80 + t130 * t79;
t141 = t54 * t85;
t140 = t54 * t89;
t137 = t61 * t85;
t136 = t62 * t85;
t135 = t63 * t85;
t132 = t80 * t85;
t131 = t80 * t89;
t127 = t83 * t91;
t126 = t85 * t91;
t125 = t89 * t91;
t123 = -t31 * rSges(5,1) - t32 * rSges(5,2);
t36 = t129 * t79 + t63 * t80;
t122 = -t35 * rSges(5,1) - t36 * rSges(5,2);
t78 = pkin(3) * t90 + pkin(2);
t92 = -pkin(10) - pkin(9);
t121 = -t60 * t78 - t61 * t92;
t120 = -t62 * t78 - t63 * t92;
t55 = t117 * t79 + t130 * t80;
t119 = -t54 * rSges(5,1) - t55 * rSges(5,2);
t118 = t146 * pkin(1) + pkin(8) * t129;
t116 = t86 * t129;
t112 = -t88 * pkin(1) + pkin(8) * t113;
t71 = t86 * t113;
t111 = -t61 * t90 + t71;
t108 = t37 * pkin(3);
t107 = pkin(3) * t116 - t62 * t92 + t63 * t78 + t118;
t105 = rSges(5,1) * t80 - rSges(5,2) * t79;
t5 = -t36 * t85 + t62 * t89;
t104 = t156 * pkin(3);
t103 = rSges(4,1) * t90 - rSges(4,2) * t86 + pkin(2);
t29 = -t125 * t83 - t55 * t85;
t102 = pkin(3) * t71 + t60 * t92 - t61 * t78 + t112;
t100 = -rSges(7,1) * t144 + rSges(7,2) * t145 + t124 * t32 - t31 * t77;
t99 = -rSges(7,1) * t142 + rSges(7,2) * t143 + t124 * t36 - t35 * t77;
t98 = -rSges(7,1) * t140 + rSges(7,2) * t141 + t124 * t55 - t54 * t77;
t97 = t113 * t90 + t61 * t86;
t96 = -rSges(6,1) * t144 + rSges(6,2) * t145 - t31 * pkin(4) + t148 * t32;
t95 = -rSges(6,1) * t142 + rSges(6,2) * t143 - t35 * pkin(4) + t148 * t36;
t94 = -rSges(6,1) * t140 + rSges(6,2) * t141 - t54 * pkin(4) + t148 * t55;
t93 = t97 * pkin(3);
t64 = t78 * t127;
t44 = (t125 * t80 + t85 * t87) * t83;
t43 = (-t126 * t80 + t87 * t89) * t83;
t38 = t63 * t90 + t116;
t30 = t126 * t83 - t55 * t89;
t10 = -t131 * t62 + t135;
t9 = t132 * t62 + t63 * t89;
t8 = -t131 * t60 + t137;
t7 = t132 * t60 + t61 * t89;
t6 = t36 * t89 + t136;
t1 = [-m(2) * (g(1) * (-t88 * rSges(2,1) - rSges(2,2) * t146) + g(2) * (rSges(2,1) * t146 - t88 * rSges(2,2))) - m(3) * (g(1) * (-t61 * rSges(3,1) + t60 * rSges(3,2) + rSges(3,3) * t113 + t112) + g(2) * (rSges(3,1) * t63 - rSges(3,2) * t62 + rSges(3,3) * t129 + t118)) - m(4) * (g(1) * (rSges(4,1) * t111 + rSges(4,2) * t97 - t61 * pkin(2) + t147 * t60 + t112) + g(2) * (rSges(4,1) * t38 + rSges(4,2) * t37 + pkin(2) * t63 - t147 * t62 + t118)) - m(5) * (g(1) * (-rSges(5,1) * t32 + rSges(5,2) * t31 - rSges(5,3) * t60 + t102) + g(2) * (rSges(5,1) * t36 - rSges(5,2) * t35 + rSges(5,3) * t62 + t107)) - m(6) * (g(1) * (rSges(6,1) * t155 + rSges(6,2) * t3 - pkin(4) * t32 - t148 * t31 + t102) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 + pkin(4) * t36 + t148 * t35 + t107)) - m(7) * (g(1) * (rSges(7,1) * t155 + rSges(7,2) * t3 - pkin(5) * t139 - t124 * t31 - t32 * t77 + t102) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + pkin(5) * t136 + t124 * t35 + t36 * t77 + t107)) -m(3) * (g(1) * (-rSges(3,1) * t62 - rSges(3,2) * t63) + g(2) * (-rSges(3,1) * t60 - rSges(3,2) * t61) + (rSges(3,1) * t91 - rSges(3,2) * t87) * t149) - m(4) * (g(1) * (-t103 * t62 - t147 * t63) - g(2) * t147 * t61 - t103 * t150 + (t103 * t91 - t147 * t87) * t149) - m(5) * (g(1) * (rSges(5,3) * t63 - t105 * t62 + t120) + g(2) * (rSges(5,3) * t61 - t105 * t60 + t121) + g(3) * t64 + (t105 * t91 + (rSges(5,3) - t92) * t87) * t149) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9 - t152 * t62 + t120) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 - t152 * t60 + t121) + g(3) * (t44 * rSges(6,1) + t43 * rSges(6,2) + t127 * t152 - t130 * t92 + t64) + (g(3) * t127 - t154) * t79 * t148) - m(7) * (g(1) * (rSges(7,1) * t10 + rSges(7,2) * t9 + pkin(5) * t135 + t120) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + pkin(5) * t137 + t121) + g(3) * (t44 * rSges(7,1) + t43 * rSges(7,2) + t64) + ((pkin(5) * t85 - t92) * t87 + t153 * t91) * t149 - t154 * t153) -m(4) * (g(1) * (rSges(4,1) * t37 - rSges(4,2) * t38) + g(2) * (-rSges(4,1) * t97 + rSges(4,2) * t111) + g(3) * (t156 * rSges(4,1) + (-t117 * t86 - t128 * t87) * rSges(4,2))) - m(5) * (g(1) * (t108 + t122) + g(2) * (-t93 + t123) + g(3) * (t104 + t119)) - m(6) * (g(1) * (t108 + t95) + g(2) * (-t93 + t96) + g(3) * (t104 + t94)) - m(7) * (g(1) * (t108 + t99) + g(2) * (t100 - t93) + g(3) * (t104 + t98)) -m(5) * (g(1) * t122 + g(2) * t123 + g(3) * t119) - m(6) * (g(1) * t95 + g(2) * t96 + g(3) * t94) - m(7) * (g(1) * t99 + g(2) * t100 + g(3) * t98) -m(6) * (g(1) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t3 + rSges(6,2) * t155) + g(3) * (rSges(6,1) * t29 + rSges(6,2) * t30)) + (-g(1) * (rSges(7,1) * t5 - rSges(7,2) * t6) - g(2) * (-rSges(7,1) * t3 + rSges(7,2) * t155) - g(3) * (rSges(7,1) * t29 + rSges(7,2) * t30) - (g(1) * t5 - g(2) * t3 + g(3) * t29) * pkin(5)) * m(7), -m(7) * (g(1) * t35 + g(2) * t31 + g(3) * t54)];
taug  = t1(:);
