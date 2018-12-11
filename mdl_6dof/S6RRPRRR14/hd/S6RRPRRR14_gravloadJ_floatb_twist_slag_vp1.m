% Calculate Gravitation load on the joints for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:09:57
% EndTime: 2018-12-10 18:10:10
% DurationCPUTime: 3.80s
% Computational Cost: add. (9138->299), mult. (9337->423), div. (0->0), fcn. (9298->30), ass. (0->146)
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t119 = sin(qJ(4));
t114 = sin(pkin(7));
t201 = sin(pkin(6));
t202 = cos(pkin(7));
t170 = t202 * t201;
t211 = cos(qJ(1));
t161 = t211 * t170;
t192 = pkin(6) + qJ(2);
t175 = cos(t192) / 0.2e1;
t193 = pkin(6) - qJ(2);
t187 = cos(t193);
t155 = t187 / 0.2e1 + t175;
t208 = sin(qJ(2));
t209 = sin(qJ(1));
t94 = -t155 * t211 + t208 * t209;
t147 = t94 * t114 - t161;
t190 = pkin(8) + qJ(4);
t171 = sin(t190) / 0.2e1;
t191 = pkin(8) - qJ(4);
t183 = sin(t191);
t151 = t171 + t183 / 0.2e1;
t174 = cos(t190) / 0.2e1;
t186 = cos(t191);
t154 = t186 / 0.2e1 + t174;
t184 = sin(t192);
t172 = t184 / 0.2e1;
t185 = sin(t193);
t153 = t172 - t185 / 0.2e1;
t210 = cos(qJ(2));
t179 = t209 * t210;
t139 = t153 * t211 + t179;
t188 = pkin(7) + pkin(14);
t168 = sin(t188) / 0.2e1;
t189 = pkin(7) - pkin(14);
t181 = sin(t189);
t149 = t168 + t181 / 0.2e1;
t145 = t149 * t201;
t169 = cos(t189) / 0.2e1;
t182 = cos(t188);
t150 = t169 + t182 / 0.2e1;
t200 = sin(pkin(14));
t214 = t139 * t200 + t145 * t211 + t150 * t94;
t115 = cos(pkin(14));
t177 = t211 * t201;
t98 = t168 - t181 / 0.2e1;
t99 = t169 - t182 / 0.2e1;
t59 = t115 * t139 - t177 * t99 - t94 * t98;
t19 = t119 * t59 - t147 * t151 + t154 * t214;
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t100 = t171 - t183 / 0.2e1;
t101 = t174 - t186 / 0.2e1;
t122 = cos(qJ(4));
t129 = -t100 * t214 - t101 * t147 + t59 * t122;
t113 = sin(pkin(8));
t116 = cos(pkin(8));
t44 = t113 * t214 + t116 * t147;
t4 = t118 * t44 + t121 * t129;
t222 = t117 * t4 - t120 * t19;
t221 = -t117 * t19 - t120 * t4;
t218 = -t118 * t129 + t121 * t44;
t140 = t155 * t209 + t208 * t211;
t180 = t211 * t210;
t141 = -t153 * t209 + t180;
t128 = t140 * t150 + t141 * t200 - t145 * t209;
t217 = t140 * t114 + t209 * t170;
t124 = t128 * t113 + t116 * t217;
t212 = rSges(7,3) + pkin(13);
t213 = rSges(6,3) + pkin(12);
t207 = t121 * pkin(5);
t203 = cos(pkin(6));
t199 = qJ(3) * t114;
t198 = t114 * t101;
t197 = t114 * t116;
t196 = t117 * t121;
t195 = t120 * t121;
t176 = t201 * t209;
t194 = t211 * pkin(1) + pkin(10) * t176;
t178 = -pkin(1) * t209 + pkin(10) * t177;
t173 = t185 / 0.2e1;
t167 = -rSges(6,1) * t121 + rSges(6,2) * t118;
t163 = t173 - t184 / 0.2e1;
t95 = t163 * t211 - t179;
t72 = -t115 * t94 + t95 * t98;
t90 = t94 * pkin(2);
t166 = t72 * pkin(3) - t199 * t95 - t90;
t96 = -t163 * t209 - t180;
t74 = -t115 * t140 + t96 * t98;
t92 = t140 * pkin(2);
t165 = t74 * pkin(3) - t199 * t96 - t92;
t102 = t175 - t187 / 0.2e1;
t152 = t172 + t173;
t84 = t102 * t98 + t115 * t152;
t97 = t152 * pkin(2);
t164 = t84 * pkin(3) - t102 * t199 + t97;
t162 = rSges(7,1) * t120 - rSges(7,2) * t117 + pkin(5);
t71 = t150 * t95 + t200 * t94;
t49 = -t113 * t71 - t197 * t95;
t73 = t140 * t200 + t150 * t96;
t50 = -t113 * t73 - t197 * t96;
t83 = t102 * t150 - t152 * t200;
t68 = -t102 * t197 - t113 * t83;
t146 = t114 * t151;
t31 = t72 * t119 + t146 * t95 - t154 * t71;
t32 = t71 * t100 + t122 * t72 + t198 * t95;
t159 = t32 * pkin(4) + t31 * pkin(12) + t166;
t33 = t74 * t119 + t146 * t96 - t154 * t73;
t34 = t73 * t100 + t122 * t74 + t198 * t96;
t158 = t34 * pkin(4) + t33 * pkin(12) + t165;
t42 = t102 * t146 + t84 * t119 - t154 * t83;
t43 = t83 * t100 + t102 * t198 + t122 * t84;
t156 = t43 * pkin(4) + t42 * pkin(12) + t164;
t148 = -t139 * pkin(2) + qJ(3) * t161 - t199 * t94 + t178;
t144 = -pkin(3) * t59 - pkin(11) * t44 + t148;
t143 = -pkin(4) * t129 + t144;
t142 = t114 * t152 - t202 * t203;
t137 = t141 * pkin(2) + t217 * qJ(3) + t194;
t136 = (g(1) * t50 + g(2) * t49 + g(3) * t68) * pkin(11);
t132 = t102 * t200 + t149 * t203 + t150 * t152;
t75 = -t102 * t115 + t152 * t98 + t203 * t99;
t130 = t100 * t132 + t101 * t142 + t75 * t122;
t62 = t115 * t141 - t140 * t98 + t176 * t99;
t126 = t62 * pkin(3) + t124 * pkin(11) + t137;
t123 = -t100 * t128 - t101 * t217 + t62 * t122;
t125 = pkin(4) * t123 + t126;
t56 = -t113 * t132 - t116 * t142;
t36 = t119 * t75 - t132 * t154 + t142 * t151;
t35 = t36 * pkin(4);
t28 = t118 * t68 + t121 * t43;
t27 = t118 * t43 - t68 * t121;
t24 = t119 * t62 + t128 * t154 - t151 * t217;
t17 = t24 * pkin(4);
t15 = t19 * pkin(4);
t14 = t118 * t56 + t121 * t130;
t13 = -t118 * t130 + t121 * t56;
t12 = t118 * t50 + t121 * t34;
t11 = t118 * t34 - t50 * t121;
t10 = t118 * t49 + t121 * t32;
t9 = t118 * t32 - t49 * t121;
t8 = t118 * t124 + t121 * t123;
t7 = t118 * t123 - t121 * t124;
t2 = t117 * t24 + t120 * t8;
t1 = -t117 * t8 + t120 * t24;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t209 - rSges(2,2) * t211) + g(2) * (rSges(2,1) * t211 - rSges(2,2) * t209)) - m(3) * (g(1) * (-rSges(3,1) * t139 + t94 * rSges(3,2) + rSges(3,3) * t177 + t178) + g(2) * (rSges(3,1) * t141 - rSges(3,2) * t140 + rSges(3,3) * t176 + t194)) - m(4) * (g(1) * (-rSges(4,1) * t59 + rSges(4,2) * t214 - rSges(4,3) * t147 + t148) + g(2) * (t62 * rSges(4,1) - rSges(4,2) * t128 + rSges(4,3) * t217 + t137)) - m(5) * (g(1) * (-rSges(5,1) * t129 + rSges(5,2) * t19 - rSges(5,3) * t44 + t144) + g(2) * (rSges(5,1) * t123 - t24 * rSges(5,2) + rSges(5,3) * t124 + t126)) - m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t218 - t19 * t213 + t143) + g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t213 * t24 + t125)) - m(7) * (g(1) * (t221 * rSges(7,1) + t222 * rSges(7,2) - t4 * pkin(5) - t19 * pkin(12) + t212 * t218 + t143) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t8 * pkin(5) + t24 * pkin(12) + t212 * t7 + t125)) -m(3) * (g(1) * (-rSges(3,1) * t140 + t96 * rSges(3,2)) + g(2) * (-rSges(3,1) * t94 + rSges(3,2) * t95) + g(3) * (rSges(3,1) * t152 + t102 * rSges(3,2))) - m(4) * (g(1) * (rSges(4,1) * t74 + rSges(4,2) * t73 - t92) + g(2) * (rSges(4,1) * t72 + rSges(4,2) * t71 - t90) + g(3) * (rSges(4,1) * t84 + rSges(4,2) * t83 + t97) + (g(1) * t96 + g(2) * t95 + g(3) * t102) * t114 * (-rSges(4,3) - qJ(3))) - m(5) * (g(1) * (rSges(5,1) * t34 - rSges(5,2) * t33 + rSges(5,3) * t50 + t165) + g(2) * (rSges(5,1) * t32 - rSges(5,2) * t31 + rSges(5,3) * t49 + t166) + g(3) * (rSges(5,1) * t43 - rSges(5,2) * t42 + rSges(5,3) * t68 + t164) + t136) - m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t11 + rSges(6,3) * t33 + t158) + g(2) * (rSges(6,1) * t10 - rSges(6,2) * t9 + rSges(6,3) * t31 + t159) + g(3) * (rSges(6,1) * t28 - rSges(6,2) * t27 + rSges(6,3) * t42 + t156) + t136) - m(7) * (g(1) * (t12 * pkin(5) + (t117 * t33 + t12 * t120) * rSges(7,1) + (-t117 * t12 + t120 * t33) * rSges(7,2) + t158 + t212 * t11) + g(2) * (t10 * pkin(5) + (t10 * t120 + t117 * t31) * rSges(7,1) + (-t10 * t117 + t120 * t31) * rSges(7,2) + t159 + t212 * t9) + g(3) * (t28 * pkin(5) + (t117 * t42 + t120 * t28) * rSges(7,1) + (-t117 * t28 + t120 * t42) * rSges(7,2) + t156 + t212 * t27) + t136) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t217 + g(2) * t147 - g(3) * t142) -m(5) * (g(1) * (-rSges(5,1) * t24 - rSges(5,2) * t123) + g(2) * (-rSges(5,1) * t19 - rSges(5,2) * t129) + g(3) * (-rSges(5,1) * t36 - rSges(5,2) * t130)) - m(6) * (g(1) * (t123 * t213 + t167 * t24 - t17) + g(2) * (t129 * t213 + t167 * t19 - t15) + g(3) * (t130 * t213 + t167 * t36 - t35)) + (-g(1) * (-t24 * t207 - t17 + t123 * pkin(12) + (t117 * t123 - t195 * t24) * rSges(7,1) + (t120 * t123 + t196 * t24) * rSges(7,2)) - g(2) * (-t19 * t207 - t15 + t129 * pkin(12) + (t117 * t129 - t19 * t195) * rSges(7,1) + (t120 * t129 + t19 * t196) * rSges(7,2)) - g(3) * (-t36 * t207 - t35 + t130 * pkin(12) + (t117 * t130 - t195 * t36) * rSges(7,1) + (t120 * t130 + t196 * t36) * rSges(7,2)) - (-g(1) * t24 - g(2) * t19 - g(3) * t36) * t118 * t212) * m(7), -m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (rSges(6,1) * t218 - rSges(6,2) * t4) + g(3) * (rSges(6,1) * t13 - rSges(6,2) * t14)) - m(7) * (g(1) * (-t162 * t7 + t212 * t8) + (t13 * t162 + t14 * t212) * g(3) + (t162 * t218 + t212 * t4) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t222 * rSges(7,1) + t221 * rSges(7,2)) + g(3) * ((-t117 * t14 + t120 * t36) * rSges(7,1) + (-t117 * t36 - t120 * t14) * rSges(7,2)))];
taug  = t3(:);
