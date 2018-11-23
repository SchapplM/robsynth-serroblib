% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_gravloadJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_matlab.m
t144 = sin(qJ(6));
t148 = cos(qJ(6));
t146 = sin(qJ(4));
t232 = pkin(8) + qJ(4);
t206 = sin(t232) / 0.2e1;
t233 = pkin(8) - qJ(4);
t220 = sin(t233);
t175 = t206 + t220 / 0.2e1;
t209 = cos(t232) / 0.2e1;
t223 = cos(t233);
t177 = t223 / 0.2e1 + t209;
t147 = sin(qJ(1));
t153 = cos(qJ(1));
t236 = pkin(6) + qJ(2);
t211 = cos(t236) / 0.2e1;
t237 = pkin(6) - qJ(2);
t225 = cos(t237);
t178 = t225 / 0.2e1 + t211;
t255 = sin(qJ(2));
t118 = t147 * t255 - t153 * t178;
t142 = sin(pkin(7));
t247 = sin(pkin(6));
t248 = cos(pkin(7));
t204 = t248 * t247;
t260 = -t118 * t142 + t153 * t204;
t208 = sin(t236) / 0.2e1;
t222 = sin(t237);
t125 = t208 - t222 / 0.2e1;
t152 = cos(qJ(2));
t119 = t125 * t153 + t147 * t152;
t234 = pkin(7) + qJ(3);
t207 = sin(t234) / 0.2e1;
t235 = pkin(7) - qJ(3);
t221 = sin(t235);
t124 = t207 - t221 / 0.2e1;
t210 = cos(t234) / 0.2e1;
t224 = cos(t235);
t127 = t210 - t224 / 0.2e1;
t151 = cos(qJ(3));
t226 = t153 * t247;
t78 = -t118 * t124 + t119 * t151 + t127 * t226;
t176 = t207 + t221 / 0.2e1;
t172 = t176 * t247;
t200 = t224 / 0.2e1 + t210;
t254 = sin(qJ(3));
t80 = t118 * t200 + t119 * t254 + t153 * t172;
t24 = t78 * t146 + t175 * t260 + t80 * t177;
t145 = sin(qJ(5));
t149 = cos(qJ(5));
t123 = t206 - t220 / 0.2e1;
t126 = t209 - t223 / 0.2e1;
t150 = cos(qJ(4));
t28 = t123 * t80 - t126 * t260 - t78 * t150;
t141 = sin(pkin(8));
t143 = cos(pkin(8));
t61 = -t141 * t80 + t143 * t260;
t7 = t145 * t61 + t149 * t28;
t269 = t144 * t7 + t148 * t24;
t268 = -t144 * t24 + t148 * t7;
t267 = t145 * t28 - t61 * t149;
t121 = t147 * t125 - t152 * t153;
t169 = t147 * t178 + t153 * t255;
t160 = -t121 * t254 - t147 * t172 + t169 * t200;
t261 = t169 * t142 + t147 * t204;
t262 = t160 * t141 + t143 * t261;
t241 = t142 * t143;
t93 = t118 * t254 - t119 * t200;
t64 = t119 * t241 - t141 * t93;
t95 = t121 * t200 + t169 * t254;
t65 = -t121 * t241 - t141 * t95;
t128 = t211 - t225 / 0.2e1;
t199 = t208 + t222 / 0.2e1;
t107 = t128 * t200 - t199 * t254;
t90 = -t107 * t141 - t128 * t241;
t259 = -g(1) * t65 - g(2) * t64 - g(3) * t90;
t253 = pkin(11) * t142;
t252 = pkin(12) * t141;
t249 = cos(pkin(6));
t244 = t141 * t145;
t243 = t141 * t149;
t242 = t142 * t126;
t240 = t144 * t149;
t239 = t148 * t149;
t227 = t147 * t247;
t238 = t153 * pkin(1) + pkin(10) * t227;
t231 = -t24 * pkin(4) - pkin(13) * t28;
t82 = -t121 * t151 - t169 * t124 - t127 * t227;
t154 = -t160 * t123 - t126 * t261 + t82 * t150;
t29 = t82 * t146 + t160 * t177 - t175 * t261;
t230 = -t29 * pkin(4) + pkin(13) * t154;
t162 = t128 * t254 + t249 * t176 + t199 * t200;
t171 = t199 * t142 - t249 * t248;
t99 = -t199 * t124 + t249 * t127 + t128 * t151;
t158 = t162 * t123 + t171 * t126 - t150 * t99;
t49 = -t146 * t99 - t162 * t177 + t171 * t175;
t229 = -t49 * pkin(4) + pkin(13) * t158;
t228 = -pkin(1) * t147 + pkin(10) * t226;
t219 = -t80 * pkin(3) + t252 * t78;
t218 = -t160 * pkin(3) + t252 * t82;
t217 = t162 * pkin(3) - t99 * t252;
t8 = t145 * t154 - t149 * t262;
t216 = g(1) * t267 + g(2) * t8;
t215 = -t118 * pkin(2) + t119 * t253;
t214 = -t169 * pkin(2) - t121 * t253;
t213 = t199 * pkin(2) - t128 * t253;
t212 = -g(1) * t24 + g(2) * t29;
t205 = -pkin(5) * t149 - pkin(14) * t145;
t94 = -t118 * t151 - t119 * t124;
t203 = t94 * pkin(3) + t215;
t96 = t121 * t124 - t169 * t151;
t202 = t96 * pkin(3) + t214;
t108 = t128 * t124 + t199 * t151;
t201 = t108 * pkin(3) + t213;
t73 = -t162 * t141 - t171 * t143;
t18 = -t145 * t158 + t149 * t73;
t196 = g(1) * t8 - g(2) * t267 - g(3) * t18;
t19 = t145 * t73 + t149 * t158;
t9 = t145 * t262 + t149 * t154;
t195 = g(1) * t9 - g(2) * t7 + g(3) * t19;
t41 = -t119 * t242 + t123 * t93 + t150 * t94;
t10 = t145 * t41 - t64 * t149;
t43 = t121 * t242 + t123 * t95 + t150 * t96;
t12 = t145 * t43 - t65 * t149;
t59 = t107 * t123 + t108 * t150 + t128 * t242;
t32 = t145 * t59 - t90 * t149;
t194 = g(1) * t12 + g(2) * t10 + g(3) * t32;
t45 = -t123 * t78 - t80 * t150;
t14 = t145 * t45 - t243 * t78;
t47 = -t123 * t82 - t160 * t150;
t16 = t145 * t47 - t243 * t82;
t54 = t99 * t123 + t162 * t150;
t34 = t145 * t54 + t99 * t243;
t193 = g(1) * t16 + g(2) * t14 + g(3) * t34;
t192 = g(1) * t29 + g(2) * t24 + g(3) * t49;
t191 = g(1) * t154 - g(2) * t28 + g(3) * t158;
t173 = t142 * t175;
t40 = -t119 * t173 + t94 * t146 - t93 * t177;
t42 = t121 * t173 + t96 * t146 - t95 * t177;
t58 = -t107 * t177 + t108 * t146 + t128 * t173;
t190 = g(1) * t42 + g(2) * t40 + g(3) * t58;
t44 = -t80 * t146 + t177 * t78;
t46 = -t160 * t146 + t177 * t82;
t53 = t162 * t146 - t99 * t177;
t189 = g(1) * t46 + g(2) * t44 + g(3) * t53;
t188 = -g(1) * t82 - g(2) * t78 + g(3) * t99;
t187 = g(1) * t121 - g(2) * t119 + g(3) * t128;
t186 = t45 * pkin(4) + pkin(13) * t44 + t219;
t185 = t47 * pkin(4) + pkin(13) * t46 + t218;
t184 = t54 * pkin(4) + pkin(13) * t53 + t217;
t183 = t41 * pkin(4) + t40 * pkin(13) + t203;
t182 = t43 * pkin(4) + t42 * pkin(13) + t202;
t181 = t59 * pkin(4) + t58 * pkin(13) + t201;
t180 = -t119 * pkin(2) + pkin(11) * t260 + t228;
t170 = -t78 * pkin(3) + pkin(12) * t61 + t180;
t167 = t28 * pkin(4) - pkin(13) * t24 + t170;
t166 = -t121 * pkin(2) + pkin(11) * t261 + t238;
t164 = t259 * pkin(12);
t157 = t82 * pkin(3) + pkin(12) * t262 + t166;
t155 = pkin(4) * t154 + t29 * pkin(13) + t157;
t35 = t149 * t54 - t99 * t244;
t33 = t145 * t90 + t149 * t59;
t17 = t149 * t47 + t244 * t82;
t15 = t149 * t45 + t244 * t78;
t13 = t145 * t65 + t149 * t43;
t11 = t145 * t64 + t149 * t41;
t3 = t144 * t29 + t148 * t9;
t2 = -t144 * t9 + t148 * t29;
t1 = t192 * t145;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t147 - g(2) * t153, g(1) * t153 + g(2) * t147, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t119 + g(2) * t121, -g(1) * t118 + g(2) * t169, -g(1) * t226 - g(2) * t227, -g(1) * t228 - g(2) * t238, 0, 0, 0, 0, 0, 0, g(1) * t78 - g(2) * t82, -g(1) * t80 + g(2) * t160, -g(1) * t260 - g(2) * t261, -g(1) * t180 - g(2) * t166, 0, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t154, t212, -g(1) * t61 - g(2) * t262, -g(1) * t170 - g(2) * t157, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, t216, -t212, -g(1) * t167 - g(2) * t155, 0, 0, 0, 0, 0, 0, -g(1) * t268 - g(2) * t3, g(1) * t269 - g(2) * t2, -t216, -g(1) * (pkin(5) * t7 + pkin(14) * t267 + t167) - g(2) * (t9 * pkin(5) + t8 * pkin(14) + t155); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t169 + g(2) * t118 - g(3) * t199, -t187, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t108, -g(1) * t95 - g(2) * t93 - g(3) * t107, t187 * t142, -g(1) * t214 - g(2) * t215 - g(3) * t213, 0, 0, 0, 0, 0, 0, -g(1) * t43 - g(2) * t41 - g(3) * t59, t190, t259, -g(1) * t202 - g(2) * t203 - g(3) * t201 + t164, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t33, t194, -t190, -g(1) * t182 - g(2) * t183 - g(3) * t181 + t164, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t148 + t144 * t42) - g(2) * (t11 * t148 + t144 * t40) - g(3) * (t144 * t58 + t148 * t33) -g(1) * (-t13 * t144 + t148 * t42) - g(2) * (-t11 * t144 + t148 * t40) - g(3) * (-t144 * t33 + t148 * t58) -t194, -g(1) * (t13 * pkin(5) + t12 * pkin(14) + t182) - g(2) * (t11 * pkin(5) + t10 * pkin(14) + t183) - g(3) * (t33 * pkin(5) + t32 * pkin(14) + t181) + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t160 + g(2) * t80 - g(3) * t162, -t188, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t47 - g(2) * t45 - g(3) * t54, t189, t188 * t141, -g(1) * t218 - g(2) * t219 - g(3) * t217, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15 - g(3) * t35, t193, -t189, -g(1) * t185 - g(2) * t186 - g(3) * t184, 0, 0, 0, 0, 0, 0, -g(1) * (t144 * t46 + t148 * t17) - g(2) * (t144 * t44 + t148 * t15) - g(3) * (t144 * t53 + t148 * t35) -g(1) * (-t144 * t17 + t148 * t46) - g(2) * (-t144 * t15 + t148 * t44) - g(3) * (-t144 * t35 + t148 * t53) -t193, -g(1) * (pkin(5) * t17 + pkin(14) * t16 + t185) - g(2) * (pkin(5) * t15 + pkin(14) * t14 + t186) - g(3) * (pkin(5) * t35 + pkin(14) * t34 + t184); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t191, 0, 0, 0, 0, 0, 0, 0, 0, t192 * t149, -t1, -t191, -g(1) * t230 - g(2) * t231 - g(3) * t229, 0, 0, 0, 0, 0, 0, -g(1) * (t144 * t154 - t29 * t239) - g(2) * (-t144 * t28 - t24 * t239) - g(3) * (t144 * t158 - t49 * t239) -g(1) * (t148 * t154 + t29 * t240) - g(2) * (-t148 * t28 + t24 * t240) - g(3) * (t148 * t158 + t49 * t240) t1, -g(1) * (t205 * t29 + t230) - g(2) * (t205 * t24 + t231) - g(3) * (t205 * t49 + t229); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t195, 0, 0, 0, 0, 0, 0, 0, 0, t196 * t148, -t196 * t144, -t195, -g(1) * (-pkin(5) * t8 + pkin(14) * t9) - g(2) * (pkin(5) * t267 - pkin(14) * t7) - g(3) * (pkin(5) * t18 + pkin(14) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t269 - g(3) * (-t144 * t19 + t148 * t49) g(1) * t3 - g(2) * t268 - g(3) * (-t144 * t49 - t148 * t19) 0, 0;];
taug_reg  = t4;
