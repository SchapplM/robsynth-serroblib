% Calculate Gravitation load on the joints for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:14:07
% EndTime: 2019-03-09 01:14:11
% DurationCPUTime: 1.78s
% Computational Cost: add. (2402->285), mult. (6917->447), div. (0->0), fcn. (9003->18), ass. (0->135)
t164 = rSges(7,3) + pkin(13);
t100 = cos(qJ(3));
t148 = sin(pkin(14));
t153 = cos(pkin(6));
t131 = t153 * t148;
t150 = cos(pkin(14));
t161 = sin(qJ(2));
t163 = cos(qJ(2));
t116 = t131 * t163 + t150 * t161;
t149 = sin(pkin(6));
t128 = t149 * t148;
t152 = cos(pkin(7));
t93 = sin(pkin(7));
t108 = -t116 * t152 + t128 * t93;
t85 = -t131 * t161 + t150 * t163;
t97 = sin(qJ(3));
t103 = -t100 * t108 + t85 * t97;
t109 = t116 * t93 + t128 * t152;
t151 = cos(pkin(8));
t92 = sin(pkin(8));
t174 = t103 * t151 - t109 * t92;
t132 = t153 * t150;
t115 = -t132 * t163 + t148 * t161;
t129 = t150 * t149;
t171 = -t115 * t152 - t93 * t129;
t84 = t132 * t161 + t148 * t163;
t104 = -t100 * t171 + t84 * t97;
t107 = t115 * t93 - t129 * t152;
t173 = t104 * t151 - t107 * t92;
t130 = t152 * t149;
t118 = t130 * t163 + t153 * t93;
t139 = t149 * t161;
t112 = -t100 * t118 + t139 * t97;
t140 = t163 * t149;
t117 = -t140 * t93 + t152 * t153;
t172 = t112 * t151 - t117 * t92;
t170 = pkin(10) * t93;
t169 = pkin(11) * t92;
t99 = cos(qJ(5));
t168 = t99 * pkin(5);
t165 = rSges(6,3) + pkin(12);
t162 = cos(qJ(4));
t122 = t161 * t130;
t80 = -t100 * t122 - t140 * t97;
t160 = t80 * t92;
t159 = t92 * t93;
t95 = sin(qJ(5));
t158 = t92 * t95;
t157 = t92 * t99;
t94 = sin(qJ(6));
t156 = t94 * t99;
t98 = cos(qJ(6));
t155 = t98 * t99;
t126 = t93 * t139;
t154 = pkin(2) * t140 + pkin(10) * t126;
t147 = t93 * t151;
t96 = sin(qJ(4));
t146 = t96 * t151;
t145 = t97 * t152;
t144 = t100 * t152;
t119 = t151 * t126;
t81 = t100 * t140 - t122 * t97;
t143 = t81 * pkin(3) + pkin(11) * t119 + t154;
t142 = t162 * t159;
t141 = t151 * t162;
t138 = -rSges(6,1) * t99 + rSges(6,2) * t95;
t66 = -t100 * t115 - t145 * t84;
t82 = t115 * pkin(2);
t137 = t66 * pkin(3) + t170 * t84 - t82;
t68 = -t100 * t116 - t145 * t85;
t83 = t116 * pkin(2);
t136 = t68 * pkin(3) + t170 * t85 - t83;
t59 = t84 * t100 + t171 * t97;
t30 = -t104 * t162 - t146 * t59;
t56 = t104 * pkin(3);
t135 = t30 * pkin(4) + t169 * t59 - t56;
t60 = t85 * t100 + t108 * t97;
t32 = -t103 * t162 - t146 * t60;
t57 = t103 * pkin(3);
t134 = t32 * pkin(4) + t169 * t60 - t57;
t76 = t100 * t139 + t118 * t97;
t46 = -t112 * t162 - t146 * t76;
t75 = t112 * pkin(3);
t133 = t46 * pkin(4) + t169 * t76 - t75;
t127 = t98 * rSges(7,1) - t94 * rSges(7,2) + pkin(5);
t124 = t92 * t126;
t50 = -t124 * t162 - t141 * t80 + t81 * t96;
t51 = t81 * t162 + (t151 * t80 + t124) * t96;
t123 = t51 * pkin(4) + t50 * pkin(12) + t143;
t65 = t115 * t97 - t144 * t84;
t47 = t147 * t84 - t65 * t92;
t67 = t116 * t97 - t144 * t85;
t48 = t147 * t85 - t67 * t92;
t25 = -t141 * t65 - t142 * t84 + t66 * t96;
t26 = t66 * t162 + (t151 * t65 + t159 * t84) * t96;
t121 = t26 * pkin(4) + t25 * pkin(12) + t137;
t27 = -t141 * t67 - t142 * t85 + t68 * t96;
t28 = t68 * t162 + (t151 * t67 + t159 * t85) * t96;
t120 = t28 * pkin(4) + t27 * pkin(12) + t136;
t111 = (g(1) * t48 + g(2) * t47 - g(3) * t160) * pkin(11);
t70 = t119 - t160;
t58 = t112 * t92 + t117 * t151;
t45 = -t112 * t96 + t141 * t76;
t41 = t103 * t92 + t109 * t151;
t40 = t104 * t92 + t107 * t151;
t39 = t76 * t162 - t172 * t96;
t38 = t162 * t172 + t76 * t96;
t37 = t38 * pkin(4);
t36 = t158 * t76 + t46 * t99;
t35 = -t157 * t76 + t46 * t95;
t34 = t51 * t99 + t70 * t95;
t33 = t51 * t95 - t70 * t99;
t31 = -t103 * t96 + t141 * t60;
t29 = -t104 * t96 + t141 * t59;
t20 = t60 * t162 - t174 * t96;
t19 = t162 * t174 + t60 * t96;
t18 = t59 * t162 - t173 * t96;
t17 = t162 * t173 + t59 * t96;
t16 = t39 * t99 + t58 * t95;
t15 = -t39 * t95 + t58 * t99;
t14 = t19 * pkin(4);
t13 = t17 * pkin(4);
t12 = t158 * t60 + t32 * t99;
t11 = -t157 * t60 + t32 * t95;
t10 = t158 * t59 + t30 * t99;
t9 = -t157 * t59 + t30 * t95;
t8 = t28 * t99 + t48 * t95;
t7 = t28 * t95 - t48 * t99;
t6 = t26 * t99 + t47 * t95;
t5 = t26 * t95 - t47 * t99;
t4 = t20 * t99 + t41 * t95;
t3 = -t20 * t95 + t41 * t99;
t2 = t18 * t99 + t40 * t95;
t1 = -t18 * t95 + t40 * t99;
t21 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t116 - t85 * rSges(3,2)) + g(2) * (-rSges(3,1) * t115 - t84 * rSges(3,2)) + g(3) * (rSges(3,1) * t140 - rSges(3,2) * t139)) - m(4) * (g(1) * (rSges(4,1) * t68 + rSges(4,2) * t67 - t83) + g(2) * (rSges(4,1) * t66 + rSges(4,2) * t65 - t82) + g(3) * (t81 * rSges(4,1) + t80 * rSges(4,2) + t154) + (rSges(4,3) * g(3) * t139 + (g(1) * t85 + g(2) * t84) * (rSges(4,3) + pkin(10))) * t93) - m(5) * (g(1) * (t28 * rSges(5,1) - t27 * rSges(5,2) + t48 * rSges(5,3) + t136) + g(2) * (t26 * rSges(5,1) - t25 * rSges(5,2) + t47 * rSges(5,3) + t137) + g(3) * (rSges(5,1) * t51 - rSges(5,2) * t50 + rSges(5,3) * t70 + t143) + t111) - m(6) * (g(1) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t27 * rSges(6,3) + t120) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t25 * rSges(6,3) + t121) + g(3) * (rSges(6,1) * t34 - rSges(6,2) * t33 + rSges(6,3) * t50 + t123) + t111) - m(7) * (g(1) * (t8 * pkin(5) + (t27 * t94 + t8 * t98) * rSges(7,1) + (t27 * t98 - t8 * t94) * rSges(7,2) + t120 + t164 * t7) + g(2) * (t6 * pkin(5) + (t25 * t94 + t6 * t98) * rSges(7,1) + (t25 * t98 - t6 * t94) * rSges(7,2) + t121 + t164 * t5) + g(3) * (t34 * pkin(5) + (t34 * t98 + t50 * t94) * rSges(7,1) + (-t34 * t94 + t50 * t98) * rSges(7,2) + t123 + t164 * t33) + t111) -m(4) * (g(1) * (-rSges(4,1) * t103 - t60 * rSges(4,2)) + g(2) * (-rSges(4,1) * t104 - t59 * rSges(4,2)) + g(3) * (-rSges(4,1) * t112 - t76 * rSges(4,2))) - m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t11 + t165 * t31 + t134) + g(2) * (rSges(6,1) * t10 - rSges(6,2) * t9 + t165 * t29 + t135) + g(3) * (rSges(6,1) * t36 - rSges(6,2) * t35 + t165 * t45 + t133)) - m(7) * (g(1) * (t12 * pkin(5) + t31 * pkin(12) + (t12 * t98 + t31 * t94) * rSges(7,1) + (-t12 * t94 + t31 * t98) * rSges(7,2) + t164 * t11 + t134) + g(2) * (t10 * pkin(5) + t29 * pkin(12) + (t10 * t98 + t29 * t94) * rSges(7,1) + (-t10 * t94 + t29 * t98) * rSges(7,2) + t164 * t9 + t135) + g(3) * (t36 * pkin(5) + t45 * pkin(12) + (t36 * t98 + t45 * t94) * rSges(7,1) + (-t36 * t94 + t45 * t98) * rSges(7,2) + t164 * t35 + t133)) + (-g(1) * (rSges(5,1) * t32 - rSges(5,2) * t31 - t57) - g(2) * (rSges(5,1) * t30 - rSges(5,2) * t29 - t56) - g(3) * (rSges(5,1) * t46 - rSges(5,2) * t45 - t75) - (g(1) * t60 + g(2) * t59 + g(3) * t76) * t92 * (rSges(5,3) + pkin(11))) * m(5), -m(5) * (g(1) * (-rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (-rSges(5,1) * t17 - rSges(5,2) * t18) + g(3) * (-rSges(5,1) * t38 - rSges(5,2) * t39)) - m(6) * (g(1) * (t138 * t19 + t165 * t20 - t14) + g(2) * (t138 * t17 + t165 * t18 - t13) + g(3) * (t138 * t38 + t165 * t39 - t37)) + (-g(1) * (-t19 * t168 - t14 + t20 * pkin(12) + (-t155 * t19 + t20 * t94) * rSges(7,1) + (t156 * t19 + t20 * t98) * rSges(7,2)) - g(2) * (-t17 * t168 - t13 + t18 * pkin(12) + (-t155 * t17 + t18 * t94) * rSges(7,1) + (t156 * t17 + t18 * t98) * rSges(7,2)) - g(3) * (-t38 * t168 - t37 + t39 * pkin(12) + (-t155 * t38 + t39 * t94) * rSges(7,1) + (t156 * t38 + t39 * t98) * rSges(7,2)) - (-g(1) * t19 - g(2) * t17 - g(3) * t38) * t95 * t164) * m(7), -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (rSges(6,1) * t15 - rSges(6,2) * t16)) - m(7) * (g(1) * (t127 * t3 + t164 * t4) + (t127 * t15 + t164 * t16) * g(3) + (t127 * t1 + t164 * t2) * g(2)) -m(7) * (g(1) * ((t19 * t98 - t4 * t94) * rSges(7,1) + (-t19 * t94 - t4 * t98) * rSges(7,2)) + g(2) * ((t17 * t98 - t2 * t94) * rSges(7,1) + (-t17 * t94 - t2 * t98) * rSges(7,2)) + g(3) * ((-t16 * t94 + t38 * t98) * rSges(7,1) + (-t16 * t98 - t38 * t94) * rSges(7,2)))];
taug  = t21(:);
