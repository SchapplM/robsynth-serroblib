% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:17
% EndTime: 2019-12-31 21:35:28
% DurationCPUTime: 4.21s
% Computational Cost: add. (2891->452), mult. (6407->587), div. (0->0), fcn. (4342->8), ass. (0->220)
t152 = sin(qJ(5));
t156 = cos(qJ(5));
t158 = cos(qJ(2));
t153 = sin(qJ(3));
t159 = cos(qJ(1));
t250 = t159 * t153;
t155 = sin(qJ(1));
t157 = cos(qJ(3));
t253 = t155 * t157;
t84 = t158 * t250 - t253;
t251 = t157 * t159;
t85 = t153 * t155 + t158 * t251;
t190 = t152 * t85 - t156 * t84;
t244 = qJD(1) * t158;
t144 = pkin(6) * t244;
t118 = qJD(2) * pkin(7) + t144;
t238 = qJD(3) * t157;
t240 = qJD(3) * t153;
t154 = sin(qJ(2));
t203 = pkin(2) * t154 - pkin(7) * t158;
t108 = t203 * qJD(2);
t183 = pkin(2) * t158 + pkin(7) * t154 + pkin(1);
t59 = qJD(1) * t108 - qJDD(1) * t183;
t147 = t158 * qJDD(1);
t230 = qJD(1) * qJD(2);
t301 = -t154 * t230 + t147;
t79 = pkin(6) * t301 + qJDD(2) * pkin(7);
t90 = t183 * qJD(1);
t211 = t118 * t238 + t153 * t79 - t157 * t59 - t240 * t90;
t197 = qJDD(4) + t211;
t291 = pkin(3) + pkin(4);
t215 = t158 * t230;
t229 = t154 * qJDD(1);
t232 = t157 * qJD(2);
t239 = qJD(3) * t154;
t308 = qJD(1) * t239 - qJDD(2);
t40 = -qJD(3) * t232 + (-t215 - t229) * t157 + t308 * t153;
t95 = qJDD(3) - t301;
t3 = pkin(8) * t40 - t291 * t95 + t197;
t309 = qJD(3) - t244;
t121 = t309 * qJD(4);
t173 = -t118 * t240 + t153 * t59 + t157 * t79 - t238 * t90;
t88 = t95 * qJ(4);
t11 = t121 + t173 + t88;
t41 = ((qJD(3) + t244) * qJD(2) + t229) * t153 + t308 * t157;
t5 = pkin(8) * t41 + t11;
t223 = -t152 * t5 + t156 * t3;
t245 = qJD(1) * t154;
t117 = -qJD(2) * pkin(2) + pkin(6) * t245;
t222 = t157 * t245;
t243 = qJD(2) * t153;
t99 = t222 + t243;
t179 = qJ(4) * t99 - t117;
t97 = t153 * t245 - t232;
t25 = -t291 * t97 + t179;
t256 = t153 * t158;
t82 = t155 * t256 + t251;
t252 = t157 * t158;
t83 = t155 * t252 - t250;
t303 = -t152 * t83 + t156 * t82;
t47 = t152 * t97 + t156 * t99;
t254 = t154 * t157;
t255 = t154 * t156;
t75 = t152 * t254 - t153 * t255;
t311 = -g(1) * t190 + g(2) * t303 - g(3) * t75 + t25 * t47 - t223;
t281 = g(3) * t158;
t202 = g(1) * t159 + g(2) * t155;
t296 = t154 * t202;
t310 = t281 - t296;
t299 = qJD(3) - qJD(5);
t189 = t152 * t99 - t156 * t97;
t307 = -t189 ^ 2 + t47 ^ 2;
t304 = t47 * t189;
t51 = -t118 * t153 - t157 * t90;
t248 = qJD(4) - t51;
t302 = qJD(4) * t153 + t144;
t231 = -qJD(5) + t309;
t300 = -t231 - qJD(5);
t234 = qJD(5) * t156;
t235 = qJD(5) * t152;
t9 = t152 * t41 - t156 * t40 + t234 * t97 - t235 * t99;
t298 = t189 * t231 - t9;
t288 = pkin(7) * t95;
t37 = pkin(3) * t97 - t179;
t297 = -t309 * t37 + t288;
t262 = qJ(4) * t153;
t295 = -t157 * t291 - t262;
t191 = t152 * t82 + t156 * t83;
t39 = t152 * t84 + t156 * t85;
t101 = t152 * t153 + t156 * t157;
t76 = t101 * t154;
t293 = -g(1) * t39 - g(2) * t191 - g(3) * t76 - t189 * t25;
t292 = t99 ^ 2;
t290 = pkin(7) - pkin(8);
t289 = pkin(3) * t95;
t285 = pkin(3) * t153;
t282 = g(3) * t154;
t280 = t99 * t97;
t171 = t158 * t101;
t279 = -qJD(1) * t171 + t101 * t299;
t225 = t156 * t256;
t258 = t152 * t157;
t278 = qJD(1) * t225 + t152 * t238 + t153 * t234 - t156 * t240 - t157 * t235 - t244 * t258;
t226 = t291 * t153;
t261 = qJ(4) * t157;
t177 = -t226 + t261;
t277 = t177 * t309 + t302;
t195 = -t261 + t285;
t276 = t195 * t309 - t302;
t275 = t108 * t153 - t183 * t238;
t274 = pkin(7) * qJD(3);
t273 = qJ(4) * t97;
t123 = t309 * qJ(4);
t52 = t118 * t157 - t153 * t90;
t32 = t123 + t52;
t272 = t309 * t32;
t271 = t309 * t52;
t270 = t309 * t97;
t269 = t309 * t99;
t28 = pkin(8) * t97 + t52;
t23 = t123 + t28;
t268 = t152 * t23;
t266 = t153 * t40;
t265 = t157 * t99;
t135 = pkin(6) * t252;
t264 = -t153 * t183 + t135;
t105 = t203 * qJD(1);
t86 = t153 * t105;
t263 = qJ(4) * t245 + t86;
t260 = qJD(2) * t99;
t259 = t105 * t157;
t257 = t153 * t154;
t249 = pkin(8) * t99 - t248;
t247 = (g(1) * t251 + g(2) * t253) * t154;
t150 = t154 ^ 2;
t246 = -t158 ^ 2 + t150;
t242 = qJD(2) * t154;
t241 = qJD(2) * t158;
t236 = qJD(4) * t157;
t233 = t117 * qJD(3);
t17 = -t291 * t309 - t249;
t228 = t152 * t3 + t156 * t5 + t17 * t234;
t227 = qJ(4) * t242 + t275;
t120 = t290 * t157;
t224 = -pkin(6) * t153 - pkin(3);
t221 = t309 * t243;
t220 = t309 * t232;
t219 = t158 * t232;
t218 = t309 * t240;
t217 = t153 * t239;
t213 = -t152 * t40 - t156 * t41;
t134 = pkin(6) * t256;
t210 = -t157 * t183 - t134;
t209 = t231 ^ 2;
t207 = -g(1) * t82 + g(2) * t84;
t206 = g(1) * t83 - g(2) * t85;
t61 = -qJ(4) * t158 + t264;
t205 = -qJD(3) * t135 + t108 * t157 + t183 * t240;
t204 = t224 * t154;
t201 = g(1) * t155 - g(2) * t159;
t142 = pkin(6) * t229;
t200 = qJDD(2) * pkin(2) - pkin(6) * t215 - t142;
t119 = t290 * t153;
t199 = -qJD(5) * t119 + t290 * t240 + (-pkin(6) * t254 + pkin(8) * t256) * qJD(1) + t263;
t167 = -pkin(8) * t252 + (-pkin(4) + t224) * t154;
t198 = t167 * qJD(1) - t120 * t299 - t259;
t196 = pkin(3) * t157 + t262;
t194 = t233 - t288;
t8 = t152 * t17 + t156 * t23;
t149 = t158 * pkin(3);
t42 = pkin(4) * t158 + t134 + t149 + (-pkin(8) * t154 + t183) * t157;
t50 = pkin(8) * t257 + t61;
t193 = -t152 * t50 + t156 * t42;
t192 = t152 * t42 + t156 * t50;
t31 = -pkin(3) * t309 + t248;
t188 = -t153 * t32 + t157 * t31;
t187 = qJ(4) * t156 - t152 * t291;
t186 = qJ(4) * t152 + t156 * t291;
t185 = -t153 * t156 + t258;
t182 = pkin(2) + t196;
t181 = t274 * t309 + t281;
t178 = t23 * t235 - t228;
t176 = -0.2e1 * pkin(1) * t230 - pkin(6) * qJDD(2);
t175 = t153 * t95 + t238 * t309;
t174 = t157 * t95 - t218;
t165 = -qJ(4) * t40 + qJD(4) * t99 + t200;
t12 = pkin(3) * t41 - t165;
t172 = -t12 - t181;
t162 = qJD(1) ^ 2;
t170 = pkin(1) * t162 + t202;
t169 = g(1) * t84 + g(2) * t82 + g(3) * t257 - t211;
t168 = -t158 * t202 - t282;
t161 = qJD(2) ^ 2;
t166 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t161 + t201;
t10 = qJD(5) * t47 + t213;
t164 = t37 * t99 + qJDD(4) - t169;
t163 = g(1) * t85 + g(2) * t83 + g(3) * t254 + t309 * t51 - t173;
t130 = qJ(4) * t254;
t93 = pkin(2) - t295;
t89 = -qJDD(5) + t95;
t71 = -t130 + (pkin(6) + t285) * t154;
t62 = t149 - t210;
t60 = t130 + (-pkin(6) - t226) * t154;
t58 = pkin(3) * t99 + t273;
t57 = qJD(1) * t204 - t259;
t56 = -pkin(6) * t222 + t263;
t29 = -t291 * t99 - t273;
t26 = (qJD(3) * t196 - t236) * t154 + (pkin(6) + t195) * t241;
t24 = qJD(2) * t204 - t205;
t22 = -t40 + t270;
t21 = -qJD(4) * t158 + (-t154 * t232 - t158 * t240) * pkin(6) + t227;
t20 = t154 * t185 * t299 + qJD(2) * t171;
t19 = -qJD(2) * t225 + qJD(5) * t76 - t238 * t255 + (-t217 + t219) * t152;
t18 = (qJD(3) * t295 + t236) * t154 + (-pkin(6) + t177) * t241;
t15 = (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t254 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t153) * t158 + t227;
t14 = pkin(8) * t217 + qJD(2) * t167 - t205;
t13 = t197 - t289;
t7 = t156 * t17 - t268;
t6 = -t291 * t41 + t165;
t1 = [qJDD(1), t201, t202, qJDD(1) * t150 + 0.2e1 * t154 * t215, 0.2e1 * t147 * t154 - 0.2e1 * t230 * t246, qJDD(2) * t154 + t158 * t161, qJDD(2) * t158 - t154 * t161, 0, t154 * t176 + t158 * t166, -t154 * t166 + t158 * t176, t99 * t219 + (-t157 * t40 - t240 * t99) * t154, (-t153 * t99 - t157 * t97) * t241 + (t266 - t157 * t41 + (t153 * t97 - t265) * qJD(3)) * t154, (t40 + t220) * t158 + (t174 + t260) * t154, (t41 - t221) * t158 + (-qJD(2) * t97 - t175) * t154, -t158 * t95 + t242 * t309, t205 * t309 + t210 * t95 + ((pkin(6) * t97 + t117 * t153) * qJD(2) + t211) * t158 + (t157 * t233 + t51 * qJD(2) - t200 * t153 + (t41 + t221) * pkin(6)) * t154 + t206, -t275 * t309 - t264 * t95 + (t117 * t232 + (t218 + t260) * pkin(6) + t173) * t158 + (-t153 * t233 - t52 * qJD(2) - t200 * t157 + (-t40 + t220) * pkin(6)) * t154 + t207, -t309 * t24 + t26 * t97 + t41 * t71 - t62 * t95 + (t243 * t37 + t13) * t158 + (-qJD(2) * t31 + t12 * t153 + t238 * t37) * t154 + t206, -t21 * t97 + t24 * t99 - t40 * t62 - t41 * t61 + t188 * t241 + (-t11 * t153 + t13 * t157 + (-t153 * t31 - t157 * t32) * qJD(3) + t201) * t154, t309 * t21 - t26 * t99 + t40 * t71 + t61 * t95 + (-t232 * t37 - t11) * t158 + (qJD(2) * t32 - t12 * t157 + t240 * t37) * t154 - t207, t11 * t61 + t32 * t21 + t12 * t71 + t37 * t26 + t13 * t62 + t31 * t24 - g(1) * (-pkin(3) * t83 - qJ(4) * t82) - g(2) * (pkin(3) * t85 + qJ(4) * t84) + (-pkin(6) * g(1) - g(2) * t183) * t159 + (-g(2) * pkin(6) + g(1) * t183) * t155, t20 * t47 + t76 * t9, -t10 * t76 - t189 * t20 - t19 * t47 - t75 * t9, t158 * t9 - t20 * t231 - t242 * t47 - t76 * t89, -t10 * t158 + t189 * t242 + t19 * t231 + t75 * t89, -t158 * t89 + t231 * t242, -(t14 * t156 - t15 * t152) * t231 - t193 * t89 + t223 * t158 - t7 * t242 + t18 * t189 + t60 * t10 + t6 * t75 + t25 * t19 + g(1) * t191 - g(2) * t39 + (-t158 * t8 + t192 * t231) * qJD(5), (qJD(5) * t193 + t14 * t152 + t15 * t156) * t231 + t192 * t89 + t178 * t158 + t8 * t242 + t18 * t47 + t60 * t9 + t6 * t76 + t25 * t20 + g(1) * t303 + g(2) * t190; 0, 0, 0, -t154 * t162 * t158, t246 * t162, t229, t147, qJDD(2), t154 * t170 - t142 - t281, t282 + (-pkin(6) * qJDD(1) + t170) * t158, t265 * t309 - t266, (-t40 - t270) * t157 + (-t269 - t41) * t153, (-t154 * t99 - t252 * t309) * qJD(1) + t175, (t154 * t97 + t256 * t309) * qJD(1) + t174, -t309 * t245, -pkin(2) * t41 + t194 * t153 + (-t281 + t200 - (t105 + t274) * t309) * t157 + (-t117 * t256 - t51 * t154 + (-t158 * t97 - t257 * t309) * pkin(6)) * qJD(1) + t247, pkin(2) * t40 + t86 * t309 + t194 * t157 + (-t117 * t252 + t52 * t154 + (-t158 * t99 - t254 * t309) * pkin(6)) * qJD(1) + (t181 - t200 - t296) * t153, -t153 * t297 + t172 * t157 - t182 * t41 + t31 * t245 + t276 * t97 + t309 * t57 + t247, t56 * t97 - t57 * t99 + (t11 + t309 * t31 + (qJD(3) * t99 - t41) * pkin(7)) * t157 + (t13 - t272 + (qJD(3) * t97 - t40) * pkin(7)) * t153 + t168, -t32 * t245 - t182 * t40 - t309 * t56 - t276 * t99 + t297 * t157 + (t172 + t296) * t153, -t31 * t57 - t32 * t56 + t276 * t37 + (qJD(3) * t188 + t11 * t157 + t13 * t153 + t168) * pkin(7) + (-t12 - t310) * t182, -t185 * t9 + t279 * t47, t10 * t185 - t101 * t9 - t189 * t279 - t278 * t47, t185 * t89 - t231 * t279 + t245 * t47, t101 * t89 - t189 * t245 + t231 * t278, -t231 * t245, -(t119 * t156 - t120 * t152) * t89 + t93 * t10 + t277 * t189 + t278 * t25 - g(3) * t171 - (t152 * t199 - t156 * t198) * t231 + t7 * t245 + (t6 + t296) * t101, (t119 * t152 + t120 * t156) * t89 + t93 * t9 + t277 * t47 + t279 * t25 - (t152 * t198 + t156 * t199) * t231 - t8 * t245 + (-t6 + t310) * t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280, -t97 ^ 2 + t292, t22, -t41 + t269, t95, -t117 * t99 + t169 + t271, t117 * t97 + t163, -t58 * t97 - t164 + t271 + 0.2e1 * t289, pkin(3) * t40 - qJ(4) * t41 + (t32 - t52) * t99 + (t31 - t248) * t97, -t37 * t97 + t58 * t99 + 0.2e1 * t121 - t163 + 0.2e1 * t88, t11 * qJ(4) - t13 * pkin(3) - t37 * t58 - t31 * t52 - g(1) * (-pkin(3) * t84 + qJ(4) * t85) - g(2) * (-pkin(3) * t82 + qJ(4) * t83) - g(3) * (-pkin(3) * t257 + t130) + t248 * t32, -t304, -t307, t298, t231 * t47 + t10, t89, t186 * t89 - t29 * t189 - (t152 * t249 - t156 * t28) * t231 + (t187 * t231 + t8) * qJD(5) + t311, t187 * t89 - t29 * t47 - (t152 * t28 + t156 * t249) * t231 + (-t186 * t231 - t268) * qJD(5) + t228 + t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95 + t280, t22, -t309 ^ 2 - t292, t164 - t272 - t289, 0, 0, 0, 0, 0, -t152 * t209 - t156 * t89 - t189 * t99, t152 * t89 - t156 * t209 - t47 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t304, t307, -t298, t300 * t47 - t213, -t89, t300 * t8 - t311, -t231 * t7 + t178 - t293;];
tau_reg = t1;
