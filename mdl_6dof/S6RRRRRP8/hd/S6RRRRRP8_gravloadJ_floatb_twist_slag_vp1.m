% Calculate Gravitation load on the joints for
% S6RRRRRP8
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:18
% EndTime: 2019-03-10 01:47:21
% DurationCPUTime: 1.43s
% Computational Cost: add. (1058->223), mult. (1788->314), div. (0->0), fcn. (2137->12), ass. (0->106)
t176 = rSges(7,1) + pkin(5);
t102 = sin(qJ(3));
t106 = cos(qJ(3));
t100 = sin(pkin(6));
t103 = sin(qJ(2));
t147 = t100 * t103;
t149 = cos(pkin(6));
t182 = -t102 * t147 + t149 * t106;
t104 = sin(qJ(1));
t145 = t100 * t106;
t107 = cos(qJ(2));
t127 = t104 * t149;
t172 = cos(qJ(1));
t77 = -t103 * t127 + t172 * t107;
t41 = -t102 * t77 + t104 * t145;
t101 = sin(qJ(5));
t105 = cos(qJ(5));
t121 = t149 * t172;
t74 = t103 * t104 - t107 * t121;
t152 = t74 * t105;
t131 = t100 * t172;
t75 = t103 * t121 + t104 * t107;
t99 = qJ(3) + qJ(4);
t96 = sin(t99);
t97 = cos(t99);
t36 = -t96 * t131 + t75 * t97;
t1 = t101 * t36 - t152;
t2 = t101 * t74 + t105 * t36;
t129 = -t97 * t131 - t75 * t96;
t156 = t105 * t129;
t162 = t101 * t129;
t181 = rSges(6,1) * t156 - rSges(6,2) * t162 + t36 * rSges(6,3);
t63 = -t96 * t147 + t149 * t97;
t154 = t105 * t63;
t160 = t101 * t63;
t64 = t97 * t147 + t149 * t96;
t180 = rSges(6,1) * t154 - rSges(6,2) * t160 + t64 * rSges(6,3);
t148 = qJ(6) * t101;
t179 = t36 * rSges(7,2) + rSges(7,3) * t162 + t129 * t148 + t176 * t156;
t178 = t64 * rSges(7,2) + rSges(7,3) * t160 + t63 * t148 + t176 * t154;
t177 = pkin(4) * t97;
t175 = rSges(7,2) + pkin(11);
t174 = rSges(6,3) + pkin(11);
t173 = -pkin(9) - rSges(4,3);
t171 = g(3) * t100;
t170 = t74 * t96;
t76 = t172 * t103 + t107 * t127;
t169 = t76 * t96;
t168 = rSges(5,1) * t129 - t36 * rSges(5,2);
t146 = t100 * t104;
t39 = -t97 * t146 + t77 * t96;
t40 = t96 * t146 + t77 * t97;
t167 = -t39 * rSges(5,1) - t40 * rSges(5,2);
t108 = -pkin(10) - pkin(9);
t95 = pkin(3) * t106 + pkin(2);
t166 = -t75 * t108 - t74 * t95;
t165 = -t77 * t108 - t76 * t95;
t164 = t63 * rSges(5,1) - t64 * rSges(5,2);
t163 = t172 * pkin(1) + pkin(8) * t146;
t161 = t101 * t39;
t158 = t101 * t97;
t155 = t105 * t39;
t153 = t107 * t96;
t151 = t76 * t105;
t150 = rSges(7,3) + qJ(6);
t144 = t100 * t107;
t143 = t103 * t108;
t142 = t105 * t107;
t78 = t95 * t144;
t141 = t78 + (pkin(11) * t96 + t177) * t144;
t138 = t102 * t146;
t135 = t101 * t144;
t134 = pkin(4) * t129 + t36 * pkin(11);
t133 = -t39 * pkin(4) + pkin(11) * t40;
t132 = t63 * pkin(4) + pkin(11) * t64;
t130 = -t104 * pkin(1) + pkin(8) * t131;
t90 = t102 * t131;
t128 = -t106 * t75 + t90;
t125 = -pkin(11) * t170 - t74 * t177 + t166;
t124 = -pkin(11) * t169 - t76 * t177 + t165;
t123 = t41 * pkin(3);
t122 = pkin(3) * t138 - t76 * t108 + t77 * t95 + t163;
t120 = rSges(5,1) * t97 - rSges(5,2) * t96;
t119 = t40 * pkin(4) + t122;
t118 = t182 * pkin(3);
t117 = rSges(4,1) * t106 - rSges(4,2) * t102 + pkin(2);
t116 = pkin(3) * t90 + t74 * t108 - t75 * t95 + t130;
t115 = -pkin(4) * t36 + t116;
t114 = t75 * t102 + t106 * t131;
t113 = t40 * rSges(7,2) - rSges(7,3) * t161 - t39 * t148 - t155 * t176 + t133;
t112 = -rSges(6,1) * t155 + rSges(6,2) * t161 + t40 * rSges(6,3) + t133;
t111 = t118 + t132;
t110 = t114 * pkin(3);
t109 = -t110 + t134;
t49 = (t101 * t103 + t97 * t142) * t100;
t48 = -t105 * t147 + t97 * t135;
t42 = t106 * t77 + t138;
t34 = t105 * t64 - t135;
t33 = t100 * t142 + t101 * t64;
t10 = t101 * t77 - t97 * t151;
t9 = -t77 * t105 - t76 * t158;
t8 = t101 * t75 - t97 * t152;
t7 = -t75 * t105 - t74 * t158;
t6 = t101 * t76 + t105 * t40;
t5 = t101 * t40 - t151;
t3 = [-m(2) * (g(1) * (-t104 * rSges(2,1) - t172 * rSges(2,2)) + g(2) * (t172 * rSges(2,1) - t104 * rSges(2,2))) - m(3) * (g(1) * (-t75 * rSges(3,1) + t74 * rSges(3,2) + rSges(3,3) * t131 + t130) + g(2) * (rSges(3,1) * t77 - rSges(3,2) * t76 + rSges(3,3) * t146 + t163)) - m(4) * (g(1) * (t128 * rSges(4,1) + t114 * rSges(4,2) - t75 * pkin(2) + t173 * t74 + t130) + g(2) * (rSges(4,1) * t42 + rSges(4,2) * t41 + pkin(2) * t77 - t173 * t76 + t163)) - m(5) * (g(1) * (-rSges(5,1) * t36 - rSges(5,2) * t129 - rSges(5,3) * t74 + t116) + g(2) * (rSges(5,1) * t40 - rSges(5,2) * t39 + rSges(5,3) * t76 + t122)) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t129 * t174 + t115) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t174 * t39 + t119)) - m(7) * (g(1) * (-t1 * t150 + t129 * t175 - t176 * t2 + t115) + g(2) * (t150 * t5 + t175 * t39 + t176 * t6 + t119)) -m(3) * (g(1) * (-rSges(3,1) * t76 - rSges(3,2) * t77) + g(2) * (-rSges(3,1) * t74 - rSges(3,2) * t75) + (rSges(3,1) * t107 - rSges(3,2) * t103) * t171) - m(4) * (g(1) * (-t117 * t76 - t173 * t77) + (-t173 * t103 + t117 * t107) * t171 + (-t117 * t74 - t173 * t75) * g(2)) - m(5) * (g(1) * (rSges(5,3) * t77 - t120 * t76 + t165) + g(2) * (rSges(5,3) * t75 - t120 * t74 + t166) + g(3) * t78 + (t120 * t107 + (rSges(5,3) - t108) * t103) * t171) - m(6) * (g(1) * (rSges(6,1) * t10 - rSges(6,2) * t9 - rSges(6,3) * t169 + t124) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t170 + t125) + g(3) * (t49 * rSges(6,1) - t48 * rSges(6,2) + (rSges(6,3) * t153 - t143) * t100 + t141)) - m(7) * (g(1) * (-rSges(7,2) * t169 + t176 * t10 + t150 * t9 + t124) + g(2) * (-rSges(7,2) * t170 + t150 * t7 + t176 * t8 + t125) + g(3) * (t176 * t49 + t150 * t48 + (rSges(7,2) * t153 - t143) * t100 + t141)) -m(4) * (g(1) * (rSges(4,1) * t41 - rSges(4,2) * t42) + g(2) * (-t114 * rSges(4,1) + t128 * rSges(4,2)) + g(3) * (t182 * rSges(4,1) + (-t149 * t102 - t103 * t145) * rSges(4,2))) - m(5) * (g(1) * (t123 + t167) + g(2) * (-t110 + t168) + g(3) * (t118 + t164)) - m(6) * (g(1) * (t112 + t123) + g(2) * (t109 + t181) + g(3) * (t111 + t180)) - m(7) * (g(1) * (t113 + t123) + g(2) * (t109 + t179) + g(3) * (t111 + t178)) -m(5) * (g(1) * t167 + g(2) * t168 + g(3) * t164) - m(6) * (g(1) * t112 + g(2) * (t134 + t181) + g(3) * (t132 + t180)) - m(7) * (g(1) * t113 + g(2) * (t134 + t179) + g(3) * (t132 + t178)) -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t33 - rSges(6,2) * t34)) - m(7) * (g(1) * (t150 * t6 - t176 * t5) + g(2) * (-t176 * t1 + t150 * t2) + g(3) * (t150 * t34 - t176 * t33)) -m(7) * (g(1) * t5 + g(2) * t1 + g(3) * t33)];
taug  = t3(:);
