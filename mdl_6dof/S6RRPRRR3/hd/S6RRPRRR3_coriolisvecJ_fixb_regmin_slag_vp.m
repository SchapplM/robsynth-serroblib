% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x33]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:19
% EndTime: 2019-03-09 13:25:35
% DurationCPUTime: 6.38s
% Computational Cost: add. (9531->450), mult. (24248->619), div. (0->0), fcn. (19121->10), ass. (0->238)
t227 = sin(pkin(11));
t232 = sin(qJ(2));
t291 = qJD(1) * t232;
t228 = cos(pkin(11));
t236 = cos(qJ(2));
t306 = t228 * t236;
t192 = qJD(1) * t306 - t227 * t291;
t337 = qJD(4) + qJD(5);
t361 = t192 - t337;
t229 = sin(qJ(6));
t233 = cos(qJ(6));
t205 = t227 * t236 + t228 * t232;
t194 = t205 * qJD(1);
t231 = sin(qJ(4));
t235 = cos(qJ(4));
t284 = t235 * qJD(2);
t159 = t194 * t231 - t284;
t161 = qJD(2) * t231 + t194 * t235;
t230 = sin(qJ(5));
t234 = cos(qJ(5));
t252 = t159 * t230 - t234 * t161;
t94 = t234 * t159 + t161 * t230;
t255 = t229 * t94 + t233 * t252;
t60 = t229 * t252 - t233 * t94;
t353 = t255 * t60;
t207 = t230 * t231 - t234 * t235;
t295 = t361 * t207;
t303 = t230 * t235;
t208 = t231 * t234 + t303;
t294 = t361 * t208;
t349 = t255 ^ 2 - t60 ^ 2;
t186 = qJD(4) - t192;
t280 = -pkin(2) * t236 - pkin(1);
t256 = t280 * qJD(1);
t211 = qJD(3) + t256;
t117 = -pkin(3) * t192 - pkin(8) * t194 + t211;
t333 = -qJ(3) - pkin(7);
t213 = t333 * t232;
t209 = qJD(1) * t213;
t327 = qJD(2) * pkin(2);
t201 = t209 + t327;
t214 = t333 * t236;
t210 = qJD(1) * t214;
t307 = t228 * t210;
t143 = t227 * t201 - t307;
t137 = qJD(2) * pkin(8) + t143;
t79 = t235 * t117 - t137 * t231;
t63 = -pkin(9) * t161 + t79;
t43 = pkin(4) * t186 + t63;
t80 = t117 * t231 + t137 * t235;
t64 = -pkin(9) * t159 + t80;
t54 = t234 * t64;
t25 = t230 * t43 + t54;
t356 = pkin(10) * t94;
t15 = t25 - t356;
t286 = qJD(6) * t229;
t13 = t15 * t286;
t197 = t227 * t210;
t142 = t201 * t228 + t197;
t136 = -qJD(2) * pkin(3) - t142;
t88 = pkin(4) * t159 + t136;
t47 = pkin(5) * t94 + t88;
t348 = -t47 * t60 + t13;
t178 = qJD(5) + t186;
t174 = qJD(6) + t178;
t283 = qJD(1) * qJD(2);
t275 = t236 * t283;
t276 = t232 * t283;
t180 = -t227 * t276 + t228 * t275;
t290 = qJD(4) * t231;
t111 = qJD(4) * t284 + t235 * t180 - t194 * t290;
t112 = qJD(4) * t161 + t180 * t231;
t239 = t252 * qJD(5) - t111 * t230 - t234 * t112;
t285 = qJD(6) * t233;
t287 = qJD(5) * t234;
t288 = qJD(5) * t230;
t39 = t234 * t111 - t230 * t112 - t159 * t287 - t161 * t288;
t8 = t229 * t239 + t233 * t39 + t252 * t286 - t94 * t285;
t347 = -t174 * t60 + t8;
t193 = t205 * qJD(2);
t179 = qJD(1) * t193;
t217 = pkin(2) * t276;
t115 = pkin(3) * t179 - pkin(8) * t180 + t217;
t105 = t235 * t115;
t268 = qJD(2) * t333;
t190 = qJD(3) * t236 + t232 * t268;
t166 = t190 * qJD(1);
t191 = -qJD(3) * t232 + t236 * t268;
t167 = t191 * qJD(1);
t108 = t166 * t228 + t167 * t227;
t240 = -t80 * qJD(4) - t108 * t231 + t105;
t22 = pkin(4) * t179 - pkin(9) * t111 + t240;
t289 = qJD(4) * t235;
t247 = t235 * t108 + t231 * t115 + t117 * t289 - t137 * t290;
t28 = -pkin(9) * t112 + t247;
t273 = t234 * t22 - t230 * t28;
t242 = -t25 * qJD(5) + t273;
t2 = pkin(5) * t179 - pkin(10) * t39 + t242;
t266 = -t230 * t22 - t234 * t28 - t43 * t287 + t64 * t288;
t3 = pkin(10) * t239 - t266;
t279 = t233 * t2 - t229 * t3;
t360 = t47 * t255 + t279;
t241 = t255 * qJD(6) - t229 * t39 + t233 * t239;
t342 = -t174 * t255 + t241;
t219 = pkin(2) * t227 + pkin(8);
t334 = pkin(9) + t219;
t267 = qJD(4) * t334;
t129 = pkin(2) * t291 + pkin(3) * t194 - pkin(8) * t192;
t148 = t209 * t228 + t197;
t297 = t231 * t129 + t235 * t148;
t312 = t192 * t231;
t359 = -pkin(9) * t312 + t231 * t267 + t297;
t121 = t235 * t129;
t358 = pkin(4) * t194 - t148 * t231 + t121 + (-pkin(9) * t192 + t267) * t235;
t357 = t290 - t312;
t355 = pkin(10) * t252;
t149 = -t207 * t229 + t208 * t233;
t329 = qJD(6) * t149 + t295 * t229 - t294 * t233;
t354 = t252 * t94;
t277 = t205 * t289;
t204 = t227 * t232 - t306;
t196 = t204 * qJD(2);
t311 = t196 * t231;
t351 = t277 - t311;
t350 = t252 ^ 2 - t94 ^ 2;
t346 = t178 * t94 + t39;
t345 = t88 * t94 + t266;
t344 = t88 * t252 + t242;
t52 = t230 * t64;
t24 = t234 * t43 - t52;
t14 = t24 + t355;
t12 = pkin(5) * t178 + t14;
t322 = t233 * t15;
t7 = t229 * t12 + t322;
t343 = -t7 * qJD(6) + t360;
t341 = -t178 * t252 + t239;
t340 = -0.2e1 * t283;
t339 = t358 * t234;
t132 = t208 * t205;
t146 = t209 * t227 - t307;
t260 = t357 * pkin(4) - t146;
t202 = t334 * t231;
t203 = t334 * t235;
t293 = -t230 * t202 + t234 * t203;
t338 = t202 * t287 + t203 * t288 + t358 * t230 + t359 * t234;
t147 = t233 * t207 + t208 * t229;
t330 = -t147 * qJD(6) + t294 * t229 + t295 * t233;
t336 = -t149 * t179 - t330 * t174;
t335 = -t295 * t178 - t179 * t208;
t332 = t234 * t63 - t52;
t141 = pkin(3) * t204 - pkin(8) * t205 + t280;
t135 = t235 * t141;
t158 = t213 * t227 - t214 * t228;
t308 = t205 * t235;
t73 = pkin(4) * t204 - pkin(9) * t308 - t158 * t231 + t135;
t151 = t235 * t158;
t296 = t231 * t141 + t151;
t309 = t205 * t231;
t82 = -pkin(9) * t309 + t296;
t328 = t230 * t73 + t234 * t82;
t326 = t194 * t60;
t325 = t194 * t255;
t324 = t194 * t94;
t323 = t194 * t252;
t321 = -t294 * pkin(5) + t260;
t320 = t111 * t231;
t318 = t159 * t186;
t317 = t159 * t194;
t316 = t161 * t186;
t315 = t161 * t194;
t310 = t196 * t235;
t305 = t230 * t179;
t304 = t230 * t233;
t302 = t231 * t179;
t301 = t233 * t179;
t164 = t235 * t179;
t238 = qJD(1) ^ 2;
t300 = t236 * t238;
t237 = qJD(2) ^ 2;
t299 = t237 * t232;
t298 = t237 * t236;
t292 = t232 ^ 2 - t236 ^ 2;
t282 = t232 * t327;
t220 = -pkin(2) * t228 - pkin(3);
t278 = t205 * t290;
t274 = qJD(6) * t12 + t3;
t130 = pkin(3) * t193 + pkin(8) * t196 + t282;
t122 = t235 * t130;
t128 = t190 * t228 + t191 * t227;
t32 = pkin(9) * t310 + pkin(4) * t193 - t128 * t231 + t122 + (-t151 + (pkin(9) * t205 - t141) * t231) * qJD(4);
t246 = t235 * t128 + t231 * t130 + t141 * t289 - t158 * t290;
t35 = -pkin(9) * t351 + t246;
t271 = -t230 * t35 + t234 * t32;
t270 = -t230 * t63 - t54;
t269 = -t230 * t82 + t234 * t73;
t265 = pkin(1) * t340;
t107 = t166 * t227 - t228 * t167;
t127 = t190 * t227 - t228 * t191;
t263 = -t234 * t202 - t203 * t230;
t157 = -t228 * t213 - t214 * t227;
t262 = t186 * t235;
t261 = -t147 * t179 - t329 * t174;
t109 = -pkin(10) * t208 + t263;
t259 = -t294 * pkin(10) - qJD(6) * t109 + t338;
t110 = -pkin(10) * t207 + t293;
t258 = pkin(5) * t194 + t295 * pkin(10) + t293 * qJD(5) + qJD(6) * t110 - t230 * t359 + t339;
t257 = t294 * t178 - t207 * t179;
t125 = pkin(4) * t309 + t157;
t253 = t107 * t205 - t158 * t179;
t133 = t207 * t205;
t83 = t233 * t132 - t133 * t229;
t84 = -t132 * t229 - t133 * t233;
t212 = -pkin(4) * t235 + t220;
t87 = t351 * pkin(4) + t127;
t251 = -t357 * t186 + t164;
t74 = pkin(4) * t112 + t107;
t250 = t230 * t32 + t234 * t35 + t73 * t287 - t82 * t288;
t248 = -t278 - t310;
t244 = t186 * t136 - t219 * t179;
t223 = pkin(4) * t234 + pkin(5);
t163 = pkin(5) * t207 + t212;
t145 = t179 * t204;
t86 = pkin(5) * t132 + t125;
t85 = pkin(4) * t161 - pkin(5) * t252;
t46 = -t196 * t303 - t230 * t278 - t288 * t309 + (t308 * t337 - t311) * t234;
t45 = -t132 * t337 + t207 * t196;
t33 = pkin(5) * t46 + t87;
t29 = -pkin(10) * t132 + t328;
t27 = pkin(5) * t204 + pkin(10) * t133 + t269;
t23 = -pkin(5) * t239 + t74;
t17 = t332 + t355;
t16 = t270 + t356;
t11 = t84 * qJD(6) + t229 * t45 + t233 * t46;
t10 = -t83 * qJD(6) - t229 * t46 + t233 * t45;
t6 = t233 * t12 - t15 * t229;
t5 = -pkin(10) * t46 + t250;
t4 = pkin(5) * t193 - pkin(10) * t45 - qJD(5) * t328 + t271;
t1 = [0, 0, 0, 0.2e1 * t232 * t275, t292 * t340, t298, -t299, 0, -pkin(7) * t298 + t232 * t265, pkin(7) * t299 + t236 * t265, -t108 * t204 + t127 * t194 + t128 * t192 + t142 * t196 - t143 * t193 + t157 * t180 + t253, t107 * t157 + t108 * t158 - t142 * t127 + t143 * t128 + (t211 + t256) * t282, t111 * t308 + t161 * t248 -(-t159 * t235 - t161 * t231) * t196 + (-t320 - t112 * t235 + (t159 * t231 - t161 * t235) * qJD(4)) * t205, t111 * t204 + t161 * t193 + t205 * t164 + t186 * t248, -t112 * t204 - t159 * t193 - t186 * t351 - t205 * t302, t186 * t193 + t145 (-t158 * t289 + t122) * t186 + t135 * t179 + (-t137 * t289 + t105) * t204 + t79 * t193 + t127 * t159 + t157 * t112 + t136 * t277 + ((-qJD(4) * t141 - t128) * t186 + (-qJD(4) * t117 - t108) * t204 - t136 * t196 + t253) * t231, t107 * t308 + t157 * t111 + t127 * t161 + t248 * t136 - t296 * t179 - t246 * t186 - t80 * t193 - t247 * t204, -t133 * t39 - t252 * t45, -t132 * t39 - t133 * t239 + t252 * t46 - t45 * t94, -t133 * t179 + t178 * t45 - t193 * t252 + t204 * t39, -t132 * t179 - t178 * t46 - t193 * t94 + t204 * t239, t178 * t193 + t145, t87 * t94 - t125 * t239 + t74 * t132 + t88 * t46 + t271 * t178 + t269 * t179 + t273 * t204 + t24 * t193 + (-t178 * t328 - t204 * t25) * qJD(5), t125 * t39 - t74 * t133 - t250 * t178 - t328 * t179 - t25 * t193 + t266 * t204 - t252 * t87 + t88 * t45, -t10 * t255 + t8 * t84, t10 * t60 + t11 * t255 + t241 * t84 - t8 * t83, t10 * t174 + t179 * t84 - t193 * t255 + t204 * t8, -t11 * t174 - t179 * t83 + t193 * t60 + t204 * t241, t174 * t193 + t145 (-t229 * t5 + t233 * t4) * t174 + (-t229 * t29 + t233 * t27) * t179 + t279 * t204 + t6 * t193 - t33 * t60 - t86 * t241 + t23 * t83 + t47 * t11 + ((-t229 * t27 - t233 * t29) * t174 - t7 * t204) * qJD(6), t47 * t10 + t13 * t204 - t7 * t193 + t23 * t84 - t33 * t255 + t86 * t8 + (-(-qJD(6) * t29 + t4) * t174 - t27 * t179 - t2 * t204) * t229 + (-(qJD(6) * t27 + t5) * t174 - t29 * t179 - t274 * t204) * t233; 0, 0, 0, -t232 * t300, t292 * t238, 0, 0, 0, t238 * pkin(1) * t232, pkin(1) * t300 (t143 - t146) * t194 + (t142 - t148) * t192 + (-t179 * t227 - t180 * t228) * pkin(2), t142 * t146 - t143 * t148 + (-t107 * t228 + t108 * t227 - t211 * t291) * pkin(2), t161 * t262 + t320 (t111 - t318) * t235 + (-t112 - t316) * t231, t186 * t262 + t302 - t315, t251 + t317, -t186 * t194, -t107 * t235 + t220 * t112 - t146 * t159 - t79 * t194 + (-t219 * t289 - t121) * t186 + (t148 * t186 + t244) * t231, t107 * t231 + t220 * t111 - t146 * t161 + t80 * t194 + (t219 * t290 + t297) * t186 + t244 * t235, t208 * t39 - t252 * t295, -t207 * t39 + t208 * t239 - t252 * t294 - t295 * t94, t323 - t335, t257 + t324, -t178 * t194, -t212 * t239 + t74 * t207 + t263 * t179 - t24 * t194 + t260 * t94 - t294 * t88 + (-t203 * t287 + (qJD(5) * t202 + t359) * t230 - t339) * t178, t178 * t338 - t293 * t179 + t25 * t194 + t74 * t208 + t212 * t39 - t252 * t260 + t295 * t88, t149 * t8 - t255 * t330, -t147 * t8 + t149 * t241 + t255 * t329 + t330 * t60, t325 - t336, t261 - t326, -t174 * t194 (t109 * t233 - t110 * t229) * t179 - t163 * t241 + t23 * t147 - t6 * t194 - t321 * t60 + t329 * t47 + (t229 * t259 - t233 * t258) * t174 -(t109 * t229 + t110 * t233) * t179 + t163 * t8 + t23 * t149 + t7 * t194 - t321 * t255 + t330 * t47 + (t229 * t258 + t233 * t259) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192 ^ 2 - t194 ^ 2, t142 * t194 - t143 * t192 + t217, 0, 0, 0, 0, 0, t251 - t317, -t186 ^ 2 * t235 - t302 - t315, 0, 0, 0, 0, 0, t257 - t324, t323 + t335, 0, 0, 0, 0, 0, t261 + t326, t325 + t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161 * t159, -t159 ^ 2 + t161 ^ 2, t111 + t318, -t112 + t316, t179, -t136 * t161 + t186 * t80 + t240, t136 * t159 + t186 * t79 - t247, -t354, t350, t346, t341, t179, -t270 * t178 + (-t161 * t94 - t178 * t288 + t234 * t179) * pkin(4) + t344, t332 * t178 + (t161 * t252 - t178 * t287 - t305) * pkin(4) + t345, t353, t349, t347, t342, t179, t223 * t301 - (t16 * t233 - t17 * t229) * t174 + t85 * t60 + (-t229 * t305 + (-t229 * t234 - t304) * t174 * qJD(5)) * pkin(4) + ((-pkin(4) * t304 - t223 * t229) * t174 - t7) * qJD(6) + t360, t85 * t255 + (-t223 * t179 - t2 + (t16 - (-qJD(5) - qJD(6)) * t230 * pkin(4)) * t174) * t229 + (-pkin(4) * t305 + (-pkin(4) * t287 - qJD(6) * t223 + t17) * t174 - t274) * t233 + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t354, t350, t346, t341, t179, t178 * t25 + t344, t178 * t24 + t345, t353, t349, t347, t342, t179 -(-t14 * t229 - t322) * t174 + (-t174 * t286 - t252 * t60 + t301) * pkin(5) + t343 (-t15 * t174 - t2) * t229 + (t14 * t174 - t274) * t233 + (-t174 * t285 - t229 * t179 - t252 * t255) * pkin(5) + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, t349, t347, t342, t179, t174 * t7 + t343, t174 * t6 - t229 * t2 - t233 * t274 + t348;];
tauc_reg  = t1;
