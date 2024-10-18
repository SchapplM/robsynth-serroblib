% Calculate Gravitation load on the joints for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR15_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR15_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:24:02
% EndTime: 2024-09-27 22:24:04
% DurationCPUTime: 0.91s
% Computational Cost: add. (908->196), mult. (1108->303), div. (0->0), fcn. (1151->26), ass. (0->133)
t102 = cos(qJ(2));
t101 = cos(qJ(3));
t100 = cos(qJ(4));
t90 = sin(pkin(6));
t165 = pkin(11) * t90;
t95 = sin(qJ(4));
t66 = -t95 * pkin(4) + t100 * t165;
t59 = t66 * t101;
t65 = -pkin(4) * t100 - t95 * t165;
t62 = pkin(3) - t65;
t96 = sin(qJ(3));
t19 = -t96 * t62 + t59;
t16 = t19 * t102;
t158 = t66 * t96;
t18 = -t62 * t101 - t158;
t17 = pkin(2) - t18;
t97 = sin(qJ(2));
t6 = -t17 * t97 + t16;
t166 = pkin(8) + pkin(9);
t164 = t97 * pkin(2);
t103 = cos(qJ(1));
t163 = -t103 / 0.2e1;
t82 = t102 * pkin(2) + pkin(1);
t161 = rSges(6,3) * t90;
t159 = t19 * t97;
t91 = sin(pkin(5));
t157 = t90 * t91;
t93 = cos(pkin(5));
t156 = t90 * t93;
t92 = cos(pkin(6));
t94 = sin(qJ(5));
t155 = t92 * t94;
t99 = cos(qJ(5));
t154 = t92 * t99;
t153 = t96 * t97;
t89 = qJ(2) + qJ(3);
t87 = qJ(4) + t89;
t79 = sin(t87);
t98 = sin(qJ(1));
t152 = t98 * t79;
t80 = cos(t87);
t151 = t98 * t80;
t85 = sin(t89);
t150 = t98 * t85;
t149 = t98 * t94;
t148 = t98 * t97;
t147 = t99 * t98;
t83 = pkin(5) + t89;
t77 = qJ(4) + t83;
t71 = cos(t77) / 0.2e1;
t84 = pkin(5) - t89;
t116 = -qJ(4) + t84;
t72 = cos(t116);
t52 = t72 / 0.2e1 + t71;
t26 = -t103 * t52 + t152;
t109 = sin(t116);
t127 = sin(t77) / 0.2e1;
t51 = t127 - t109 / 0.2e1;
t146 = -t26 * rSges(5,1) + (-t103 * t51 - t151) * rSges(5,2);
t141 = t103 * t80;
t142 = t103 * t79;
t27 = -t98 * t52 - t142;
t145 = t27 * rSges(5,1) + (t98 * t51 - t141) * rSges(5,2);
t144 = (t127 + t109 / 0.2e1) * rSges(5,1) + (t71 - t72 / 0.2e1) * rSges(5,2);
t143 = t102 * t91;
t140 = t103 * t85;
t139 = t103 * t94;
t138 = t103 * t97;
t137 = t103 * t99;
t136 = t17 * t102;
t135 = t98 * t102;
t134 = t103 * t102;
t133 = pkin(10) + t166;
t132 = t79 * t161;
t131 = t94 * t157;
t130 = t99 * t157;
t129 = t93 * t149;
t128 = t93 * t147;
t41 = t92 * t129 - t137;
t42 = -t92 * t128 - t139;
t47 = t92 * t139 + t128;
t122 = t92 * t137;
t48 = -t122 + t129;
t126 = (t41 * t79 - t47 * t80) * rSges(6,1) + (-t42 * t79 + t48 * t80) * rSges(6,2) + (-t93 * t152 + t141) * t161;
t121 = t93 * t139;
t43 = t121 * t92 + t147;
t44 = t122 * t93 - t149;
t45 = -t93 * t137 + t92 * t149;
t46 = t92 * t147 + t121;
t125 = (-t43 * t79 - t45 * t80) * rSges(6,1) + (-t44 * t79 - t46 * t80) * rSges(6,2) + (t93 * t142 + t151) * t161;
t124 = ((-t79 * t154 - t80 * t94) * rSges(6,2) + (-t79 * t155 + t80 * t99) * rSges(6,1) + t132) * t91;
t123 = t103 * t157;
t86 = cos(t89);
t120 = t86 * rSges(4,1) + t82;
t119 = t80 * rSges(5,1) + pkin(3) * t86 + t82;
t118 = t91 * t159 + t124;
t113 = t136 + t159;
t117 = pkin(1) + t113 + t132;
t115 = -t41 * t80 - t47 * t79;
t114 = -t44 * t80 + t46 * t79;
t24 = t65 * t101 - t158;
t25 = t65 * t96 + t59;
t112 = t102 * t24 - t25 * t97;
t81 = t101 * pkin(3) + pkin(2);
t111 = -pkin(3) * t153 + t102 * t81;
t110 = -t51 * rSges(5,1) - (t102 * t96 * pkin(3) + t97 * t81) * t93 + (rSges(5,3) + t133) * t91;
t57 = -t93 * t135 - t138;
t108 = t93 * t134 - t148;
t73 = sin(t83);
t74 = sin(t84);
t63 = t73 - t74;
t107 = -t93 * t164 - t63 * rSges(4,1) / 0.2e1 + (t166 + rSges(4,3)) * t91;
t106 = pkin(3) * (t101 * t102 - t153);
t105 = t91 * (t92 * pkin(11) + t133) + t6 * t93 + (t80 * t156 + t91 * t92) * rSges(6,3);
t75 = cos(t83);
t76 = cos(t84);
t64 = t76 + t75;
t37 = t64 * t163 + t150;
t38 = -t140 - t98 * t64 / 0.2e1;
t104 = g(1) * (t38 * rSges(4,1) + (-t103 * t86 + t98 * t63 / 0.2e1) * rSges(4,2)) + g(2) * (-t37 * rSges(4,1) + (t63 * t163 - t98 * t86) * rSges(4,2)) + g(3) * ((t73 / 0.2e1 + t74 / 0.2e1) * rSges(4,1) + (t75 / 0.2e1 - t76 / 0.2e1) * rSges(4,2));
t68 = -pkin(3) * t85 - t164;
t58 = -t93 * t148 + t134;
t56 = -t93 * t138 - t135;
t40 = t93 * t106;
t39 = t111 * t93;
t10 = t94 * t123 - t43 * t80 + t45 * t79;
t9 = t98 * t130 + t42 * t80 + t48 * t79;
t8 = t25 * t102 + t24 * t97;
t7 = t18 * t97 + t16;
t4 = t112 * t93;
t3 = (t102 * t18 - t159) * t93;
t2 = t113 * t93;
t1 = [-m(2) * (g(1) * (-t98 * rSges(2,1) - t103 * rSges(2,2)) + g(2) * (t103 * rSges(2,1) - t98 * rSges(2,2))) - m(3) * (g(1) * (t56 * rSges(3,1) - rSges(3,2) * t108 - t98 * pkin(1)) + g(2) * (t58 * rSges(3,1) + t57 * rSges(3,2) + t103 * pkin(1)) + (g(1) * t103 + g(2) * t98) * t91 * (rSges(3,3) + pkin(8))) - m(4) * (g(1) * (t37 * rSges(4,2) + t107 * t103 - t120 * t98) + g(2) * (t38 * rSges(4,2) + t120 * t103 + t107 * t98)) - m(5) * (g(1) * (t26 * rSges(5,2) + t110 * t103 - t119 * t98) + g(2) * (t27 * rSges(5,2) + t119 * t103 + t110 * t98)) - m(6) * (g(1) * (t10 * rSges(6,1) + t114 * rSges(6,2) - t117 * t98 + (rSges(6,2) * t130 + t105) * t103) + g(2) * (t115 * rSges(6,1) + t9 * rSges(6,2) + t117 * t103 + (rSges(6,1) * t131 + t105) * t98)), -m(3) * (g(1) * (t57 * rSges(3,1) - t58 * rSges(3,2)) + g(2) * (rSges(3,1) * t108 + t56 * rSges(3,2)) + g(3) * (rSges(3,1) * t102 - rSges(3,2) * t97) * t91) - m(4) * ((g(1) * t57 + g(2) * t108 + g(3) * t143) * pkin(2) + t104) - m(5) * (g(1) * (t103 * t68 - t98 * t39 + t145) + g(2) * (t103 * t39 + t98 * t68 + t146) + g(3) * (t111 * t91 + t144)) - m(6) * (g(1) * (t6 * t103 - t2 * t98 + t126) + g(2) * (t2 * t103 + t6 * t98 + t125) + g(3) * (t91 * t136 + t118)), -m(4) * t104 - m(5) * (g(1) * (-pkin(3) * t140 - t98 * t40 + t145) + g(2) * (-pkin(3) * t150 + t103 * t40 + t146) + g(3) * (t106 * t91 + t144)) - m(6) * (g(1) * (t7 * t103 + t3 * t98 + t126) + g(2) * (-t3 * t103 + t7 * t98 + t125) + g(3) * (-t18 * t143 + t118)), -m(5) * (g(1) * t145 + g(2) * t146 + g(3) * t144) - m(6) * (g(1) * (t8 * t103 + t4 * t98 + t126) + g(2) * (-t4 * t103 + t8 * t98 + t125) + g(3) * (-t112 * t91 + t124)), -m(6) * (g(1) * (t9 * rSges(6,1) + (-t98 * t131 - t115) * rSges(6,2)) + g(2) * ((-t99 * t123 - t114) * rSges(6,1) + t10 * rSges(6,2)) + g(3) * ((rSges(6,1) * t99 - rSges(6,2) * t94) * t156 + ((t80 * t154 - t79 * t94) * rSges(6,1) + (-t80 * t155 - t79 * t99) * rSges(6,2)) * t91))];
taug = t1(:);
