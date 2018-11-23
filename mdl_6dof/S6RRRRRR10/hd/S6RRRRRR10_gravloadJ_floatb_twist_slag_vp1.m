% Calculate Gravitation load on the joints for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 10:29:22
% EndTime: 2018-11-23 10:29:37
% DurationCPUTime: 5.00s
% Computational Cost: add. (11566->362), mult. (11677->520), div. (0->0), fcn. (11549->30), ass. (0->159)
t142 = sin(qJ(6));
t145 = cos(qJ(6));
t144 = sin(qJ(4));
t210 = pkin(8) + qJ(4);
t195 = sin(t210) / 0.2e1;
t211 = pkin(8) - qJ(4);
t203 = sin(t211);
t173 = t195 + t203 / 0.2e1;
t198 = cos(t210) / 0.2e1;
t206 = cos(t211);
t176 = t206 / 0.2e1 + t198;
t150 = cos(qJ(1));
t214 = pkin(6) + qJ(2);
t200 = cos(t214) / 0.2e1;
t215 = pkin(6) - qJ(2);
t208 = cos(t215);
t178 = t208 / 0.2e1 + t200;
t234 = sin(qJ(2));
t235 = sin(qJ(1));
t117 = -t150 * t178 + t234 * t235;
t140 = sin(pkin(7));
t225 = sin(pkin(6));
t226 = cos(pkin(7));
t194 = t226 * t225;
t240 = -t117 * t140 + t150 * t194;
t197 = sin(t214) / 0.2e1;
t205 = sin(t215);
t124 = t197 - t205 / 0.2e1;
t149 = cos(qJ(2));
t118 = t150 * t124 + t149 * t235;
t212 = pkin(7) + qJ(3);
t196 = sin(t212) / 0.2e1;
t213 = pkin(7) - qJ(3);
t204 = sin(t213);
t123 = t196 - t204 / 0.2e1;
t199 = cos(t212) / 0.2e1;
t207 = cos(t213);
t126 = t199 - t207 / 0.2e1;
t148 = cos(qJ(3));
t209 = t150 * t225;
t77 = -t117 * t123 + t118 * t148 + t126 * t209;
t174 = t196 + t204 / 0.2e1;
t169 = t174 * t225;
t177 = t207 / 0.2e1 + t199;
t233 = sin(qJ(3));
t79 = t117 * t177 + t118 * t233 + t150 * t169;
t23 = t144 * t77 + t173 * t240 + t79 * t176;
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t122 = t195 - t203 / 0.2e1;
t125 = t198 - t206 / 0.2e1;
t147 = cos(qJ(4));
t27 = t122 * t79 - t125 * t240 - t77 * t147;
t139 = sin(pkin(8));
t141 = cos(pkin(8));
t60 = -t139 * t79 + t141 * t240;
t6 = t143 * t60 + t146 * t27;
t248 = t142 * t6 + t145 * t23;
t247 = -t142 * t23 + t145 * t6;
t246 = t143 * t27 - t60 * t146;
t166 = t150 * t234 + t178 * t235;
t184 = -t124 * t235 + t150 * t149;
t157 = t166 * t177 - t169 * t235 + t184 * t233;
t241 = t166 * t140 + t235 * t194;
t152 = t157 * t139 + t141 * t241;
t236 = rSges(7,3) + pkin(14);
t237 = rSges(6,3) + pkin(13);
t232 = pkin(11) * t140;
t231 = pkin(12) * t139;
t230 = t146 * pkin(5);
t227 = cos(pkin(6));
t222 = t139 * t143;
t221 = t139 * t146;
t220 = t140 * t125;
t219 = t140 * t141;
t218 = t142 * t146;
t217 = t145 * t146;
t201 = t225 * t235;
t216 = t150 * pkin(1) + pkin(10) * t201;
t202 = -pkin(1) * t235 + pkin(10) * t209;
t193 = -rSges(6,1) * t146 + rSges(6,2) * t143;
t44 = -t122 * t77 - t147 * t79;
t73 = t79 * pkin(3);
t192 = t44 * pkin(4) + t231 * t77 - t73;
t81 = -t123 * t166 - t126 * t201 + t148 * t184;
t46 = -t122 * t81 - t147 * t157;
t75 = t157 * pkin(3);
t191 = t46 * pkin(4) + t231 * t81 - t75;
t127 = t200 - t208 / 0.2e1;
t175 = t197 + t205 / 0.2e1;
t159 = t127 * t233 + t174 * t227 + t175 * t177;
t98 = -t123 * t175 + t126 * t227 + t127 * t148;
t53 = t98 * t122 + t147 * t159;
t96 = t159 * pkin(3);
t190 = t53 * pkin(4) - t231 * t98 + t96;
t113 = t117 * pkin(2);
t93 = -t117 * t148 - t118 * t123;
t189 = t93 * pkin(3) + t118 * t232 - t113;
t115 = t166 * pkin(2);
t95 = -t123 * t184 - t148 * t166;
t188 = t95 * pkin(3) + t184 * t232 - t115;
t107 = t127 * t123 + t148 * t175;
t121 = t175 * pkin(2);
t187 = t107 * pkin(3) - t127 * t232 + t121;
t186 = t145 * rSges(7,1) - t142 * rSges(7,2) + pkin(5);
t92 = t117 * t233 - t118 * t177;
t63 = t118 * t219 - t92 * t139;
t94 = t166 * t233 - t177 * t184;
t64 = -t94 * t139 + t184 * t219;
t106 = t127 * t177 - t175 * t233;
t89 = -t106 * t139 - t127 * t219;
t170 = t140 * t173;
t39 = -t118 * t170 + t93 * t144 - t176 * t92;
t40 = -t118 * t220 + t122 * t92 + t147 * t93;
t182 = t40 * pkin(4) + t39 * pkin(13) + t189;
t41 = t95 * t144 - t170 * t184 - t176 * t94;
t42 = t122 * t94 + t147 * t95 - t184 * t220;
t181 = t42 * pkin(4) + t41 * pkin(13) + t188;
t57 = -t106 * t176 + t107 * t144 + t127 * t170;
t58 = t106 * t122 + t107 * t147 + t127 * t220;
t180 = t58 * pkin(4) + t57 * pkin(13) + t187;
t172 = -t118 * pkin(2) + pkin(11) * t240 + t202;
t168 = -t77 * pkin(3) + pkin(12) * t60 + t172;
t167 = t27 * pkin(4) + t168;
t165 = t140 * t175 - t226 * t227;
t163 = t184 * pkin(2) + pkin(11) * t241 + t216;
t162 = (g(1) * t64 + g(2) * t63 + g(3) * t89) * pkin(12);
t155 = t122 * t159 + t125 * t165 - t147 * t98;
t154 = t81 * pkin(3) + pkin(12) * t152 + t163;
t151 = -t122 * t157 - t125 * t241 + t81 * t147;
t153 = pkin(4) * t151 + t154;
t72 = -t139 * t159 - t141 * t165;
t52 = t144 * t159 - t176 * t98;
t48 = -t144 * t98 - t159 * t176 + t165 * t173;
t47 = t48 * pkin(4);
t45 = -t144 * t157 + t176 * t81;
t43 = -t144 * t79 + t176 * t77;
t34 = t146 * t53 - t222 * t98;
t33 = t143 * t53 + t221 * t98;
t32 = t143 * t89 + t146 * t58;
t31 = t143 * t58 - t89 * t146;
t28 = t144 * t81 + t157 * t176 - t173 * t241;
t21 = t28 * pkin(4);
t19 = t23 * pkin(4);
t18 = t143 * t72 + t146 * t155;
t17 = -t143 * t155 + t146 * t72;
t16 = t146 * t46 + t222 * t81;
t15 = t143 * t46 - t221 * t81;
t14 = t146 * t44 + t222 * t77;
t13 = t143 * t44 - t221 * t77;
t12 = t143 * t64 + t146 * t42;
t11 = t143 * t42 - t64 * t146;
t10 = t143 * t63 + t146 * t40;
t9 = t143 * t40 - t63 * t146;
t8 = t143 * t152 + t146 * t151;
t7 = t143 * t151 - t146 * t152;
t2 = t142 * t28 + t145 * t8;
t1 = -t142 * t8 + t145 * t28;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t235 - t150 * rSges(2,2)) + g(2) * (t150 * rSges(2,1) - rSges(2,2) * t235)) - m(3) * (g(1) * (-t118 * rSges(3,1) + t117 * rSges(3,2) + rSges(3,3) * t209 + t202) + g(2) * (rSges(3,1) * t184 - rSges(3,2) * t166 + rSges(3,3) * t201 + t216)) - m(4) * (g(1) * (-rSges(4,1) * t77 + t79 * rSges(4,2) + rSges(4,3) * t240 + t172) + g(2) * (t81 * rSges(4,1) - rSges(4,2) * t157 + rSges(4,3) * t241 + t163)) - m(5) * (g(1) * (t27 * rSges(5,1) + rSges(5,2) * t23 + t60 * rSges(5,3) + t168) + g(2) * (rSges(5,1) * t151 - t28 * rSges(5,2) + rSges(5,3) * t152 + t154)) - m(6) * (g(1) * (t6 * rSges(6,1) - rSges(6,2) * t246 - t23 * t237 + t167) + g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t237 * t28 + t153)) - m(7) * (g(1) * (t247 * rSges(7,1) - t248 * rSges(7,2) + t6 * pkin(5) - t23 * pkin(13) + t236 * t246 + t167) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t8 * pkin(5) + t28 * pkin(13) + t236 * t7 + t153)) -m(3) * (g(1) * (-rSges(3,1) * t166 - rSges(3,2) * t184) + g(2) * (-rSges(3,1) * t117 - rSges(3,2) * t118) + g(3) * (rSges(3,1) * t175 + t127 * rSges(3,2))) - m(4) * (g(1) * (rSges(4,1) * t95 + rSges(4,2) * t94 - t115) + g(2) * (rSges(4,1) * t93 + rSges(4,2) * t92 - t113) + g(3) * (rSges(4,1) * t107 + rSges(4,2) * t106 + t121) + (-g(1) * t184 - g(2) * t118 + g(3) * t127) * t140 * (-rSges(4,3) - pkin(11))) - m(5) * (g(1) * (rSges(5,1) * t42 - rSges(5,2) * t41 + rSges(5,3) * t64 + t188) + g(2) * (rSges(5,1) * t40 - rSges(5,2) * t39 + rSges(5,3) * t63 + t189) + g(3) * (rSges(5,1) * t58 - rSges(5,2) * t57 + rSges(5,3) * t89 + t187) + t162) - m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t11 + rSges(6,3) * t41 + t181) + g(2) * (rSges(6,1) * t10 - rSges(6,2) * t9 + rSges(6,3) * t39 + t182) + g(3) * (rSges(6,1) * t32 - rSges(6,2) * t31 + rSges(6,3) * t57 + t180) + t162) - m(7) * (g(1) * (t12 * pkin(5) + (t12 * t145 + t142 * t41) * rSges(7,1) + (-t12 * t142 + t145 * t41) * rSges(7,2) + t181 + t236 * t11) + g(2) * (t10 * pkin(5) + (t10 * t145 + t142 * t39) * rSges(7,1) + (-t10 * t142 + t145 * t39) * rSges(7,2) + t182 + t236 * t9) + g(3) * (t32 * pkin(5) + (t142 * t57 + t145 * t32) * rSges(7,1) + (-t142 * t32 + t145 * t57) * rSges(7,2) + t180 + t236 * t31) + t162) -m(4) * (g(1) * (-rSges(4,1) * t157 - rSges(4,2) * t81) + g(2) * (-rSges(4,1) * t79 - rSges(4,2) * t77) + g(3) * (rSges(4,1) * t159 + t98 * rSges(4,2))) - m(6) * (g(1) * (rSges(6,1) * t16 - rSges(6,2) * t15 + t237 * t45 + t191) + g(2) * (rSges(6,1) * t14 - rSges(6,2) * t13 + t237 * t43 + t192) + g(3) * (rSges(6,1) * t34 - rSges(6,2) * t33 + t237 * t52 + t190)) - m(7) * (g(1) * (t16 * pkin(5) + t45 * pkin(13) + (t142 * t45 + t145 * t16) * rSges(7,1) + (-t142 * t16 + t145 * t45) * rSges(7,2) + t236 * t15 + t191) + g(2) * (t14 * pkin(5) + t43 * pkin(13) + (t14 * t145 + t142 * t43) * rSges(7,1) + (-t14 * t142 + t145 * t43) * rSges(7,2) + t236 * t13 + t192) + g(3) * (t34 * pkin(5) + t52 * pkin(13) + (t142 * t52 + t145 * t34) * rSges(7,1) + (-t142 * t34 + t145 * t52) * rSges(7,2) + t236 * t33 + t190)) + (-g(1) * (rSges(5,1) * t46 - rSges(5,2) * t45 - t75) - g(2) * (rSges(5,1) * t44 - rSges(5,2) * t43 - t73) - g(3) * (rSges(5,1) * t53 - rSges(5,2) * t52 + t96) - (-g(1) * t81 - g(2) * t77 + g(3) * t98) * t139 * (-rSges(5,3) - pkin(12))) * m(5), -m(5) * (g(1) * (-rSges(5,1) * t28 - rSges(5,2) * t151) + g(2) * (-rSges(5,1) * t23 + rSges(5,2) * t27) + g(3) * (-rSges(5,1) * t48 - rSges(5,2) * t155)) - m(6) * (g(1) * (t151 * t237 + t193 * t28 - t21) + g(2) * (t193 * t23 - t237 * t27 - t19) + g(3) * (t155 * t237 + t193 * t48 - t47)) + (-g(1) * (-t28 * t230 - t21 + t151 * pkin(13) + (t142 * t151 - t217 * t28) * rSges(7,1) + (t145 * t151 + t218 * t28) * rSges(7,2)) - g(2) * (-t23 * t230 - t19 - t27 * pkin(13) + (-t142 * t27 - t217 * t23) * rSges(7,1) + (-t145 * t27 + t218 * t23) * rSges(7,2)) - g(3) * (-t48 * t230 - t47 + t155 * pkin(13) + (t142 * t155 - t217 * t48) * rSges(7,1) + (t145 * t155 + t218 * t48) * rSges(7,2)) - (-g(1) * t28 - g(2) * t23 - g(3) * t48) * t143 * t236) * m(7), -m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (rSges(6,1) * t246 + rSges(6,2) * t6) + g(3) * (rSges(6,1) * t17 - rSges(6,2) * t18)) - m(7) * (g(1) * (-t186 * t7 + t236 * t8) + (t186 * t17 + t236 * t18) * g(3) + (t186 * t246 - t236 * t6) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (t248 * rSges(7,1) + t247 * rSges(7,2)) + g(3) * ((-t142 * t18 + t145 * t48) * rSges(7,1) + (-t142 * t48 - t145 * t18) * rSges(7,2)))];
taug  = t3(:);
