% Calculate Gravitation load on the joints for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:51:59
% EndTime: 2019-03-09 19:52:05
% DurationCPUTime: 2.22s
% Computational Cost: add. (1375->261), mult. (3252->394), div. (0->0), fcn. (4099->16), ass. (0->102)
t129 = cos(pkin(7));
t148 = cos(qJ(3));
t111 = t129 * t148;
t88 = sin(pkin(6));
t94 = cos(qJ(1));
t137 = t88 * t94;
t87 = sin(pkin(7));
t128 = t87 * t137;
t130 = cos(pkin(6));
t149 = cos(qJ(2));
t113 = t130 * t149;
t146 = sin(qJ(2));
t147 = sin(qJ(1));
t64 = -t113 * t94 + t146 * t147;
t112 = t130 * t146;
t65 = t112 * t94 + t147 * t149;
t92 = sin(qJ(3));
t24 = t111 * t64 + t128 * t148 + t65 * t92;
t120 = t88 * t129;
t164 = t94 * t120 - t64 * t87;
t119 = t92 * t129;
t25 = -t119 * t64 - t92 * t128 + t148 * t65;
t85 = pkin(13) + qJ(5);
t82 = sin(t85);
t83 = cos(t85);
t5 = -t164 * t82 + t25 * t83;
t91 = sin(qJ(6));
t93 = cos(qJ(6));
t168 = -t24 * t93 + t5 * t91;
t167 = -t24 * t91 - t5 * t93;
t103 = t113 * t147 + t146 * t94;
t97 = t103 * t87 + t147 * t120;
t123 = t88 * t147;
t162 = -t103 * t129 + t87 * t123;
t66 = -t112 * t147 + t149 * t94;
t28 = -t148 * t162 + t66 * t92;
t121 = t87 * t130;
t122 = t88 * t146;
t124 = t88 * t149;
t46 = -t111 * t124 - t121 * t148 + t122 * t92;
t159 = g(1) * t28 + g(2) * t24 + g(3) * t46;
t163 = -t164 * t83 - t25 * t82;
t150 = rSges(7,3) + pkin(12);
t161 = g(1) * t66 + g(2) * t65;
t158 = pkin(5) * t83;
t152 = t87 * pkin(10);
t86 = sin(pkin(13));
t145 = t164 * t86;
t143 = t82 * t87;
t142 = t83 * t87;
t141 = t83 * t91;
t140 = t83 * t93;
t139 = t86 * t87;
t89 = cos(pkin(13));
t138 = t87 * t89;
t80 = pkin(4) * t89 + pkin(3);
t90 = -pkin(11) - qJ(4);
t136 = -t24 * t80 - t25 * t90;
t29 = t66 * t148 + t162 * t92;
t135 = -t28 * t80 - t29 * t90;
t47 = t92 * t121 + (t119 * t149 + t146 * t148) * t88;
t134 = -t46 * t80 - t47 * t90;
t118 = t87 * t122;
t133 = pkin(2) * t124 + pkin(10) * t118;
t132 = t94 * pkin(1) + pkin(9) * t123;
t131 = qJ(4) + rSges(5,3);
t36 = t111 * t65 - t64 * t92;
t37 = -t119 * t65 - t148 * t64;
t59 = t64 * pkin(2);
t127 = -t36 * t90 + t37 * t80 - t59;
t38 = -t103 * t92 + t111 * t66;
t39 = -t103 * t148 - t119 * t66;
t61 = t103 * pkin(2);
t126 = -t38 * t90 + t39 * t80 - t61;
t115 = -pkin(1) * t147 + pkin(9) * t137;
t108 = t86 * t118;
t110 = t129 * t146;
t57 = (t110 * t148 + t149 * t92) * t88;
t58 = (-t110 * t92 + t148 * t149) * t88;
t114 = pkin(4) * t108 - t57 * t90 + t58 * t80 + t133;
t109 = -rSges(6,1) * t83 + rSges(6,2) * t82;
t106 = rSges(7,1) * t93 - rSges(7,2) * t91 + pkin(5);
t104 = -t65 * pkin(2) + pkin(10) * t164 + t115;
t102 = pkin(4) * t145 + t24 * t90 - t25 * t80 + t104;
t100 = t161 * (pkin(4) * t86 + pkin(10)) * t87;
t98 = t66 * pkin(2) + pkin(10) * t97 + t132;
t95 = t97 * t86;
t96 = pkin(4) * t95 - t28 * t90 + t29 * t80 + t98;
t63 = -t124 * t87 + t129 * t130;
t35 = t118 * t82 + t58 * t83;
t34 = -t118 * t83 + t58 * t82;
t19 = t47 * t83 + t63 * t82;
t18 = -t47 * t82 + t63 * t83;
t13 = t143 * t66 + t39 * t83;
t12 = -t142 * t66 + t39 * t82;
t11 = t143 * t65 + t37 * t83;
t10 = -t142 * t65 + t37 * t82;
t9 = t29 * t83 + t82 * t97;
t8 = t29 * t82 - t83 * t97;
t2 = t28 * t91 + t9 * t93;
t1 = t28 * t93 - t9 * t91;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t147 - rSges(2,2) * t94) + g(2) * (rSges(2,1) * t94 - rSges(2,2) * t147)) - m(3) * (g(1) * (-rSges(3,1) * t65 + rSges(3,2) * t64 + rSges(3,3) * t137 + t115) + g(2) * (t66 * rSges(3,1) - rSges(3,2) * t103 + rSges(3,3) * t123 + t132)) - m(4) * (g(1) * (-rSges(4,1) * t25 + rSges(4,2) * t24 + rSges(4,3) * t164 + t104) + g(2) * (t29 * rSges(4,1) - t28 * rSges(4,2) + rSges(4,3) * t97 + t98)) - m(5) * (g(1) * (-t25 * pkin(3) + (-t25 * t89 + t145) * rSges(5,1) + (t164 * t89 + t25 * t86) * rSges(5,2) - t131 * t24 + t104) + g(2) * (t29 * pkin(3) + (t29 * t89 + t95) * rSges(5,1) + (-t29 * t86 + t89 * t97) * rSges(5,2) + t131 * t28 + t98)) - m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t163 - rSges(6,3) * t24 + t102) + g(2) * (t9 * rSges(6,1) - t8 * rSges(6,2) + t28 * rSges(6,3) + t96)) - m(7) * (g(1) * (rSges(7,1) * t167 + rSges(7,2) * t168 - pkin(5) * t5 + t150 * t163 + t102) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t9 * pkin(5) + t150 * t8 + t96)) -m(3) * (g(1) * (-rSges(3,1) * t103 - t66 * rSges(3,2)) + g(2) * (-rSges(3,1) * t64 - rSges(3,2) * t65) + g(3) * (rSges(3,1) * t149 - rSges(3,2) * t146) * t88) - m(4) * (g(1) * (rSges(4,1) * t39 - rSges(4,2) * t38 - t61) + g(2) * (rSges(4,1) * t37 - rSges(4,2) * t36 - t59) + g(3) * (rSges(4,1) * t58 - rSges(4,2) * t57 + t133) + (rSges(4,3) * g(3) * t122 + t161 * (rSges(4,3) + pkin(10))) * t87) - m(5) * (g(1) * (t39 * pkin(3) - t61 + t66 * t152 + (t139 * t66 + t39 * t89) * rSges(5,1) + (t138 * t66 - t39 * t86) * rSges(5,2) + t131 * t38) + g(2) * (t37 * pkin(3) - t59 + t65 * t152 + (t139 * t65 + t37 * t89) * rSges(5,1) + (t138 * t65 - t37 * t86) * rSges(5,2) + t131 * t36) + g(3) * (t58 * pkin(3) + (t58 * t89 + t108) * rSges(5,1) + (t118 * t89 - t58 * t86) * rSges(5,2) + t131 * t57 + t133)) - m(6) * (g(1) * (rSges(6,1) * t13 - rSges(6,2) * t12 + rSges(6,3) * t38 + t126) + g(2) * (rSges(6,1) * t11 - rSges(6,2) * t10 + rSges(6,3) * t36 + t127) + g(3) * (rSges(6,1) * t35 - rSges(6,2) * t34 + rSges(6,3) * t57 + t114) + t100) - m(7) * (g(1) * (t13 * pkin(5) + (t13 * t93 + t38 * t91) * rSges(7,1) + (-t13 * t91 + t38 * t93) * rSges(7,2) + t126 + t150 * t12) + g(2) * (t11 * pkin(5) + (t11 * t93 + t36 * t91) * rSges(7,1) + (-t11 * t91 + t36 * t93) * rSges(7,2) + t127 + t150 * t10) + g(3) * (t35 * pkin(5) + (t35 * t93 + t57 * t91) * rSges(7,1) + (-t35 * t91 + t57 * t93) * rSges(7,2) + t150 * t34 + t114) + t100) -m(4) * (g(1) * (-rSges(4,1) * t28 - rSges(4,2) * t29) + g(2) * (-rSges(4,1) * t24 - rSges(4,2) * t25) + g(3) * (-rSges(4,1) * t46 - rSges(4,2) * t47)) - m(5) * ((g(1) * t29 + g(2) * t25 + g(3) * t47) * t131 + t159 * (-rSges(5,1) * t89 + rSges(5,2) * t86 - pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t29 + t109 * t28 + t135) + g(2) * (rSges(6,3) * t25 + t109 * t24 + t136) + g(3) * (rSges(6,3) * t47 + t109 * t46 + t134)) + (-g(1) * (-t28 * t158 + (-t140 * t28 + t29 * t91) * rSges(7,1) + (t141 * t28 + t29 * t93) * rSges(7,2) + t135) - g(2) * (-t24 * t158 + (-t140 * t24 + t25 * t91) * rSges(7,1) + (t141 * t24 + t25 * t93) * rSges(7,2) + t136) - g(3) * (-t46 * t158 + (-t140 * t46 + t47 * t91) * rSges(7,1) + (t141 * t46 + t47 * t93) * rSges(7,2) + t134) + t159 * t82 * t150) * m(7) (-m(5) - m(6) - m(7)) * t159, -m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (rSges(6,1) * t163 - rSges(6,2) * t5) + g(3) * (rSges(6,1) * t18 - rSges(6,2) * t19)) - m(7) * (g(1) * (-t106 * t8 + t150 * t9) + (t106 * t18 + t150 * t19) * g(3) + (t106 * t163 + t150 * t5) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-rSges(7,1) * t168 + rSges(7,2) * t167) + g(3) * ((-t19 * t91 + t46 * t93) * rSges(7,1) + (-t19 * t93 - t46 * t91) * rSges(7,2)))];
taug  = t3(:);
