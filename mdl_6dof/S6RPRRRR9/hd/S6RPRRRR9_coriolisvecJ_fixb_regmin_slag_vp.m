% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x34]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:19
% EndTime: 2019-03-09 07:25:35
% DurationCPUTime: 5.96s
% Computational Cost: add. (5562->428), mult. (12511->621), div. (0->0), fcn. (8913->8), ass. (0->219)
t191 = sin(qJ(5));
t193 = sin(qJ(3));
t270 = qJD(1) * t193;
t180 = qJD(4) + t270;
t196 = cos(qJ(4));
t197 = cos(qJ(3));
t269 = qJD(1) * t197;
t246 = t196 * t269;
t192 = sin(qJ(4));
t268 = qJD(3) * t192;
t153 = t246 + t268;
t164 = pkin(3) * t193 - pkin(8) * t197 + qJ(2);
t136 = t164 * qJD(1);
t198 = -pkin(1) - pkin(7);
t178 = t198 * qJD(1) + qJD(2);
t163 = t193 * t178;
t142 = qJD(3) * pkin(8) + t163;
t83 = t196 * t136 - t142 * t192;
t65 = -pkin(9) * t153 + t83;
t60 = pkin(4) * t180 + t65;
t195 = cos(qJ(5));
t258 = t196 * qJD(3);
t151 = t192 * t269 - t258;
t84 = t136 * t192 + t142 * t196;
t66 = -pkin(9) * t151 + t84;
t64 = t195 * t66;
t20 = t191 * t60 + t64;
t94 = t195 * t151 + t153 * t191;
t328 = pkin(10) * t94;
t15 = t20 - t328;
t190 = sin(qJ(6));
t260 = qJD(6) * t190;
t13 = t15 * t260;
t194 = cos(qJ(6));
t217 = t151 * t191 - t195 * t153;
t48 = t190 * t217 - t194 * t94;
t292 = t178 * t197;
t143 = -qJD(3) * pkin(3) - t292;
t105 = pkin(4) * t151 + t143;
t53 = pkin(5) * t94 + t105;
t333 = -t53 * t48 + t13;
t257 = qJD(1) * qJD(3);
t182 = t197 * t257;
t263 = qJD(4) * t197;
t242 = t192 * t263;
t208 = -t193 * t258 - t242;
t109 = qJD(1) * t208 + qJD(4) * t258;
t223 = pkin(3) * t197 + pkin(8) * t193;
t149 = qJD(3) * t223 + qJD(2);
t122 = t149 * qJD(1);
t114 = t196 * t122;
t206 = -t84 * qJD(4) + t114;
t266 = qJD(3) * t197;
t25 = -pkin(9) * t109 + (pkin(4) * qJD(1) - t178 * t192) * t266 + t206;
t267 = qJD(3) * t193;
t245 = t192 * t267;
t110 = -qJD(1) * t245 + t153 * qJD(4);
t243 = t197 * t258;
t264 = qJD(4) * t196;
t265 = qJD(4) * t192;
t210 = t192 * t122 + t136 * t264 - t142 * t265 + t178 * t243;
t27 = -pkin(9) * t110 + t210;
t236 = -t191 * t27 + t195 * t25;
t204 = -qJD(5) * t20 + t236;
t261 = qJD(5) * t195;
t262 = qJD(5) * t191;
t34 = t195 * t109 - t191 * t110 - t151 * t261 - t153 * t262;
t2 = pkin(5) * t182 - pkin(10) * t34 + t204;
t332 = -t190 * t2 + t333;
t219 = t190 * t94 + t194 * t217;
t324 = t219 * t48;
t319 = t219 ^ 2 - t48 ^ 2;
t176 = qJD(5) + t180;
t170 = qJD(6) + t176;
t202 = qJD(5) * t217 - t109 * t191 - t195 * t110;
t259 = qJD(6) * t194;
t8 = t190 * t202 + t194 * t34 + t217 * t260 - t94 * t259;
t318 = -t170 * t48 + t8;
t233 = -t191 * t25 - t195 * t27 - t60 * t261 + t66 * t262;
t3 = pkin(10) * t202 - t233;
t249 = -t190 * t3 + t194 * t2;
t331 = t53 * t219 + t249;
t203 = qJD(6) * t219 - t190 * t34 + t194 * t202;
t313 = -t170 * t219 + t203;
t247 = t192 * t270;
t307 = pkin(8) + pkin(9);
t248 = qJD(4) * t307;
t160 = t223 * qJD(1);
t282 = t196 * t197;
t275 = t192 * t160 + t178 * t282;
t330 = pkin(9) * t247 + t192 * t248 + t275;
t141 = t196 * t160;
t286 = t192 * t197;
t251 = t178 * t286;
t284 = t193 * t196;
t255 = pkin(9) * t284;
t329 = t196 * t248 - t251 + t141 + (pkin(4) * t197 + t255) * qJD(1);
t327 = pkin(10) * t217;
t325 = t217 * t94;
t156 = t191 * t196 + t192 * t195;
t308 = qJD(4) + qJD(5);
t322 = t308 * t156;
t155 = t191 * t192 - t195 * t196;
t310 = t155 * t193;
t277 = -qJD(1) * t310 - t308 * t155;
t134 = t156 * qJD(1);
t276 = t193 * t134 + t322;
t241 = t196 * t263;
t321 = t241 - t245;
t320 = t217 ^ 2 - t94 ^ 2;
t317 = t176 * t94 + t34;
t316 = t105 * t94 + t233;
t62 = t191 * t66;
t19 = t195 * t60 - t62;
t14 = t19 + t327;
t12 = pkin(5) * t176 + t14;
t299 = t194 * t15;
t5 = t190 * t12 + t299;
t315 = -qJD(6) * t5 + t331;
t314 = t105 * t217 + t204;
t312 = -t176 * t217 + t202;
t256 = 0.2e1 * qJD(1);
t311 = t329 * t195;
t224 = -t163 + (t247 + t265) * pkin(4);
t171 = t307 * t192;
t172 = t307 * t196;
t274 = -t191 * t171 + t195 * t172;
t283 = t193 * t198;
t175 = t196 * t283;
t273 = t192 * t164 + t175;
t309 = t171 * t261 + t172 * t262 + t329 * t191 + t195 * t330;
t100 = t194 * t155 + t156 * t190;
t306 = -qJD(6) * t100 - t276 * t190 + t277 * t194;
t101 = -t155 * t190 + t156 * t194;
t305 = qJD(6) * t101 + t277 * t190 + t276 * t194;
t304 = t195 * t65 - t62;
t302 = t276 * pkin(5) + t224;
t104 = -pkin(9) * t286 + t273;
t148 = t196 * t164;
t285 = t192 * t198;
t239 = pkin(4) - t285;
t91 = -pkin(9) * t282 + t193 * t239 + t148;
t301 = t195 * t104 + t191 * t91;
t300 = t194 * t12;
t298 = t155 * qJD(1) - t156 * t266 + t308 * t310;
t126 = t155 * t197;
t297 = qJD(3) * t126 + t193 * t322 + t134;
t296 = t109 * t192;
t295 = t143 * t192;
t294 = t151 * t180;
t293 = t153 * t180;
t291 = t180 * t192;
t290 = t180 * t193;
t289 = t180 * t196;
t288 = t190 * t191;
t287 = t191 * t194;
t281 = t197 * t198;
t199 = qJD(3) ^ 2;
t280 = t199 * t193;
t279 = t199 * t197;
t200 = qJD(1) ^ 2;
t278 = t200 * qJ(2);
t189 = t197 ^ 2;
t272 = t193 ^ 2 - t189;
t271 = -t199 - t200;
t254 = pkin(4) * qJD(5) * t170;
t252 = qJD(2) * t256;
t250 = t192 * t283;
t186 = -pkin(4) * t196 - pkin(3);
t79 = pkin(4) * t110 + t178 * t267;
t238 = qJD(6) * t12 + t3;
t129 = t196 * t149;
t42 = t129 + (-t175 + (pkin(9) * t197 - t164) * t192) * qJD(4) + (t197 * t239 + t255) * qJD(3);
t207 = -qJD(4) * t250 + t192 * t149 + t164 * t264 + t198 * t243;
t51 = -t321 * pkin(9) + t207;
t235 = -t191 * t51 + t195 * t42;
t234 = -t191 * t65 - t64;
t232 = -t104 * t191 + t195 * t91;
t230 = -t195 * t171 - t172 * t191;
t150 = pkin(4) * t286 - t281;
t229 = t151 + t258;
t228 = -t153 + t268;
t227 = qJD(4) * t193 + qJD(1);
t78 = -pkin(10) * t155 + t274;
t226 = pkin(5) * t269 + t277 * pkin(10) + t274 * qJD(5) + qJD(6) * t78 - t191 * t330 + t311;
t77 = -pkin(10) * t156 + t230;
t225 = t276 * pkin(10) - qJD(6) * t77 + t309;
t123 = t156 * t193;
t222 = qJD(6) * t123 + t297;
t221 = -qJD(6) * t310 - t298;
t30 = pkin(5) * t193 + pkin(10) * t126 + t232;
t124 = t156 * t197;
t38 = -pkin(10) * t124 + t301;
t220 = t190 * t30 + t194 * t38;
t73 = t194 * t124 - t126 * t190;
t74 = -t124 * t190 - t126 * t194;
t216 = qJD(1) * t189 - t290;
t185 = pkin(4) * t195 + pkin(5);
t215 = pkin(4) * t287 + t185 * t190;
t214 = -pkin(4) * t288 + t185 * t194;
t213 = -pkin(8) * t266 + t143 * t193;
t212 = -t104 * t262 + t191 * t42 + t195 * t51 + t91 * t261;
t111 = t321 * pkin(4) + t198 * t267;
t174 = t193 * t182;
t119 = pkin(5) * t155 + t186;
t99 = pkin(5) * t124 + t150;
t67 = pkin(4) * t153 - pkin(5) * t217;
t59 = -t262 * t286 + (t308 * t282 - t245) * t195 + t208 * t191;
t57 = qJD(3) * t310 - t197 * t322;
t39 = pkin(5) * t59 + t111;
t18 = -pkin(5) * t202 + t79;
t17 = t304 + t327;
t16 = t234 + t328;
t11 = qJD(6) * t74 + t190 * t57 + t194 * t59;
t10 = -qJD(6) * t73 - t190 * t59 + t194 * t57;
t7 = -pkin(10) * t59 + t212;
t6 = pkin(5) * t266 - pkin(10) * t57 - qJD(5) * t301 + t235;
t4 = -t15 * t190 + t300;
t1 = [0, 0, 0, 0, t252, qJ(2) * t252, -0.2e1 * t174, 0.2e1 * t272 * t257, -t280, -t279, 0, -t198 * t280 + (qJ(2) * t266 + qJD(2) * t193) * t256, -t198 * t279 + (-qJ(2) * t267 + qJD(2) * t197) * t256, t109 * t282 + t153 * t208 (t151 * t196 + t153 * t192) * t267 + (-t296 - t110 * t196 + (t151 * t192 - t153 * t196) * qJD(4)) * t197, -t180 * t242 + t109 * t193 + (t153 * t197 + t196 * t216) * qJD(3), -t180 * t241 - t110 * t193 + (-t151 * t197 - t192 * t216) * qJD(3), t180 * t266 + t174, -t110 * t281 + t114 * t193 + t129 * t180 + (t143 * t282 - t273 * t180 - t84 * t193) * qJD(4) + ((t151 * t198 - t295) * t193 + (-t180 * t285 + (t148 - t250) * qJD(1) + t83) * t197) * qJD(3), -t207 * t180 - t210 * t193 + (-t198 * t109 - t143 * t265) * t197 + ((-t273 * qJD(1) - t84) * t197 + (t198 * t153 + (-t143 + t292) * t196) * t193) * qJD(3), -t126 * t34 - t217 * t57, -t124 * t34 - t126 * t202 + t217 * t59 - t57 * t94, t176 * t57 + t193 * t34 + (-qJD(1) * t126 - t217) * t266, -t176 * t59 + t193 * t202 + (-qJD(1) * t124 - t94) * t266, t176 * t266 + t174, t235 * t176 + t236 * t193 + t111 * t94 - t150 * t202 + t79 * t124 + t105 * t59 + (-t176 * t301 - t193 * t20) * qJD(5) + (qJD(1) * t232 + t19) * t266, -t212 * t176 + t233 * t193 - t111 * t217 + t150 * t34 - t79 * t126 + t105 * t57 + (-t301 * qJD(1) - t20) * t266, -t10 * t219 + t74 * t8, t10 * t48 + t11 * t219 + t203 * t74 - t73 * t8, t10 * t170 + t193 * t8 + (qJD(1) * t74 - t219) * t266, -t11 * t170 + t193 * t203 + (-qJD(1) * t73 + t48) * t266, t170 * t266 + t174 (-t190 * t7 + t194 * t6) * t170 + t249 * t193 - t39 * t48 - t99 * t203 + t18 * t73 + t53 * t11 + (-t170 * t220 - t193 * t5) * qJD(6) + ((-t190 * t38 + t194 * t30) * qJD(1) + t4) * t266, t53 * t10 + t13 * t193 + t18 * t74 - t39 * t219 + t99 * t8 + (-(-qJD(6) * t38 + t6) * t170 - t2 * t193) * t190 + (-(qJD(6) * t30 + t7) * t170 - t238 * t193) * t194 + (-qJD(1) * t220 - t5) * t266; 0, 0, 0, 0, -t200, -t278, 0, 0, 0, 0, 0, t271 * t193, t271 * t197, 0, 0, 0, 0, 0, -t197 * t110 - t227 * t289 + (t151 * t193 + (-t180 - t270) * t286) * qJD(3), -t109 * t197 + t227 * t291 + (-t180 * t282 + (t153 - t246) * t193) * qJD(3), 0, 0, 0, 0, 0, t197 * t202 + t298 * t176 + (-t123 * t269 + t193 * t94) * qJD(3), -t197 * t34 + t297 * t176 + (-t193 * t217 + t269 * t310) * qJD(3), 0, 0, 0, 0, 0, t197 * t203 + (t190 * t222 - t194 * t221) * t170 + ((-t123 * t194 + t190 * t310) * t269 - t193 * t48) * qJD(3), -t197 * t8 + (t190 * t221 + t194 * t222) * t170 + (-(-t123 * t190 - t194 * t310) * t269 - t193 * t219) * qJD(3); 0, 0, 0, 0, 0, 0, t197 * t200 * t193, -t272 * t200, 0, 0, 0, -t197 * t278, t193 * t278, t153 * t289 + t296 (t109 - t294) * t196 + (-t110 - t293) * t192, t180 * t264 + (t180 * t284 + t197 * t228) * qJD(1), -t180 * t265 + (-t192 * t290 + t197 * t229) * qJD(1), -t180 * t269, -pkin(3) * t110 - t141 * t180 + (t180 * t286 - t193 * t229) * t178 + (-pkin(8) * t289 + t295) * qJD(4) + (t192 * t213 - t83 * t197) * qJD(1), -pkin(3) * t109 + t275 * t180 + t228 * t163 + (pkin(8) * t291 + t143 * t196) * qJD(4) + (t196 * t213 + t84 * t197) * qJD(1), t156 * t34 - t217 * t277, -t155 * t34 + t156 * t202 + t217 * t276 - t277 * t94, t277 * t176 + (qJD(3) * t156 + t217) * t269, -t276 * t176 + (-qJD(3) * t155 + t94) * t269, -t176 * t269, t79 * t155 - t186 * t202 + t224 * t94 + (-t172 * t261 + (qJD(5) * t171 + t330) * t191 - t311) * t176 + t276 * t105 + (qJD(3) * t230 - t19) * t269, t79 * t156 + t186 * t34 - t224 * t217 + t309 * t176 + t277 * t105 + (-qJD(3) * t274 + t20) * t269, t101 * t8 - t219 * t306, -t100 * t8 + t101 * t203 + t219 * t305 + t306 * t48, t306 * t170 + (qJD(3) * t101 + t219) * t269, -t305 * t170 + (-qJD(3) * t100 - t48) * t269, -t170 * t269, t18 * t100 - t119 * t203 + t305 * t53 - t302 * t48 + (t190 * t225 - t194 * t226) * t170 + ((-t190 * t78 + t194 * t77) * qJD(3) - t4) * t269, t18 * t101 + t119 * t8 + t306 * t53 - t302 * t219 + (t190 * t226 + t194 * t225) * t170 + (-(t190 * t77 + t194 * t78) * qJD(3) + t5) * t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153 * t151, -t151 ^ 2 + t153 ^ 2, t109 + t294, -t110 + t293, t182, -qJD(3) * t251 - t143 * t153 + t180 * t84 + t206, t143 * t151 + t180 * t83 - t210, -t325, t320, t317, t312, t182, -t234 * t176 + (-t153 * t94 - t176 * t262 + t182 * t195) * pkin(4) + t314, t304 * t176 + (t153 * t217 - t176 * t261 - t182 * t191) * pkin(4) + t316, t324, t319, t318, t313, t182, t214 * t182 - (t16 * t194 - t17 * t190) * t170 + t67 * t48 + (-t190 * t195 - t287) * t254 + (-t170 * t215 - t5) * qJD(6) + t331, -t215 * t182 - t194 * t3 + (t16 * t190 + t17 * t194) * t170 + t67 * t219 - (t194 * t195 - t288) * t254 + (-t170 * t214 - t300) * qJD(6) + t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t325, t320, t317, t312, t182, t176 * t20 + t314, t176 * t19 + t316, t324, t319, t318, t313, t182 -(-t14 * t190 - t299) * t170 + (-t170 * t260 + t182 * t194 - t217 * t48) * pkin(5) + t315 (-t15 * t170 - t2) * t190 + (t14 * t170 - t238) * t194 + (-t170 * t259 - t182 * t190 - t217 * t219) * pkin(5) + t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t319, t318, t313, t182, t170 * t5 + t315, t170 * t4 - t194 * t238 + t332;];
tauc_reg  = t1;
