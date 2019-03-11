% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:00
% EndTime: 2019-03-09 15:36:15
% DurationCPUTime: 5.54s
% Computational Cost: add. (5883->460), mult. (14673->619), div. (0->0), fcn. (10274->8), ass. (0->235)
t224 = sin(qJ(6));
t227 = cos(qJ(6));
t225 = sin(qJ(3));
t226 = sin(qJ(2));
t301 = qJD(1) * t226;
t280 = t225 * t301;
t228 = cos(qJ(3));
t291 = t228 * qJD(2);
t173 = t280 - t291;
t299 = qJD(2) * t225;
t175 = t228 * t301 + t299;
t222 = sin(pkin(10));
t223 = cos(pkin(10));
t121 = t223 * t173 + t175 * t222;
t248 = -t173 * t222 + t223 * t175;
t69 = t121 * t224 + t227 * t248;
t229 = cos(qJ(2));
t289 = qJD(1) * qJD(2);
t272 = t229 * t289;
t296 = qJD(3) * t225;
t276 = t226 * t296;
t288 = qJD(2) * qJD(3);
t137 = -qJD(1) * t276 + (t272 + t288) * t228;
t297 = qJD(2) * t229;
t278 = t225 * t297;
t295 = qJD(3) * t228;
t239 = t226 * t295 + t278;
t138 = qJD(1) * t239 + t225 * t288;
t84 = t137 * t222 + t223 * t138;
t85 = t137 * t223 - t138 * t222;
t16 = qJD(6) * t69 + t224 * t85 - t227 * t84;
t300 = qJD(1) * t229;
t203 = -qJD(3) + t300;
t290 = -qJD(6) - t203;
t336 = t290 * t69;
t363 = t16 + t336;
t292 = qJD(6) * t227;
t293 = qJD(6) * t224;
t266 = -t121 * t292 - t224 * t84 - t227 * t85 + t248 * t293;
t353 = -t227 * t121 + t224 * t248;
t335 = t290 * t353;
t362 = t266 + t335;
t361 = -t353 ^ 2 + t69 ^ 2;
t190 = t203 * qJD(5);
t210 = t226 * t289;
t215 = pkin(7) * t300;
t188 = qJD(2) * pkin(8) + t215;
t183 = -pkin(2) * t229 - pkin(8) * t226 - pkin(1);
t163 = t183 * qJD(1);
t321 = t225 * t163;
t130 = t188 * t228 + t321;
t257 = pkin(2) * t226 - pkin(8) * t229;
t177 = t257 * qJD(2);
t164 = qJD(1) * t177;
t263 = pkin(7) * t210;
t309 = -t228 * t164 - t225 * t263;
t235 = -qJD(3) * t130 - t309;
t38 = pkin(3) * t210 - qJ(4) * t137 - qJD(4) * t175 + t235;
t243 = t163 * t295 + t225 * t164 - t188 * t296;
t232 = -t228 * t263 + t243;
t44 = -qJ(4) * t138 - qJD(4) * t173 + t232;
t12 = t222 * t38 + t223 * t44;
t286 = qJ(5) * t210 + t12;
t9 = -t190 + t286;
t6 = pkin(9) * t84 + t9;
t345 = pkin(4) + pkin(5);
t284 = t345 * t226;
t262 = qJD(2) * t284;
t342 = t222 * t44 - t223 * t38;
t7 = -pkin(9) * t85 - qJD(1) * t262 + t342;
t281 = -t224 * t6 + t227 * t7;
t187 = -qJD(2) * pkin(2) + pkin(7) * t301;
t240 = -pkin(3) * t173 - qJD(4) - t187;
t234 = qJ(5) * t248 + t240;
t30 = -t345 * t121 + t234;
t360 = -t30 * t69 + t281;
t358 = t30 * t353;
t357 = t69 * t353;
t273 = t222 * t296;
t279 = t225 * t300;
t322 = t223 * t228;
t310 = -t222 * t279 - t223 * t295 + t300 * t322 + t273;
t356 = pkin(9) * t121;
t355 = t121 * t203;
t167 = t222 * t228 + t223 * t225;
t157 = t167 * qJD(3);
t311 = t167 * t300 - t157;
t354 = -t215 + (-t279 + t296) * pkin(3);
t352 = t248 ^ 2;
t351 = -0.2e1 * t289;
t350 = pkin(9) * t248;
t48 = pkin(4) * t121 - t234;
t349 = t248 * t48;
t317 = t228 * t229;
t242 = pkin(3) * t226 - qJ(4) * t317;
t176 = t257 * qJD(1);
t306 = pkin(7) * t280 + t228 * t176;
t108 = qJD(1) * t242 + t306;
t159 = t225 * t176;
t318 = t226 * t228;
t319 = t225 * t229;
t125 = t159 + (-pkin(7) * t318 - qJ(4) * t319) * qJD(1);
t343 = -qJ(4) - pkin(8);
t269 = qJD(3) * t343;
t294 = qJD(4) * t228;
t152 = t225 * t269 + t294;
t153 = -qJD(4) * t225 + t228 * t269;
t330 = (t108 - t153) * t223 + (-t125 + t152) * t222;
t104 = t223 * t152 + t222 * t153;
t57 = t222 * t108 + t223 * t125;
t51 = qJ(5) * t301 + t57;
t329 = t104 - t51;
t97 = -qJ(4) * t173 + t130;
t334 = t222 * t97;
t129 = t228 * t163 - t188 * t225;
t96 = -qJ(4) * t175 + t129;
t46 = t223 * t96 - t334;
t312 = qJD(5) - t46;
t346 = -t310 * qJ(5) + qJD(5) * t167 - t354;
t344 = pkin(7) * t225;
t166 = t222 * t225 - t322;
t249 = t227 * t166 - t167 * t224;
t341 = qJD(6) * t249 - t311 * t224 - t310 * t227;
t118 = t166 * t224 + t167 * t227;
t340 = qJD(6) * t118 - t310 * t224 + t311 * t227;
t206 = pkin(7) * t317;
t298 = qJD(2) * t226;
t307 = t228 * t177 + t298 * t344;
t58 = -t226 * t294 + t242 * qJD(2) + (-t206 + (qJ(4) * t226 - t183) * t225) * qJD(3) + t307;
t308 = t225 * t177 + t183 * t295;
t73 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t318 + (-qJD(4) * t226 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t229) * t225 + t308;
t24 = t222 * t58 + t223 * t73;
t339 = t311 * t345 + t346;
t338 = t311 * pkin(4) + t346;
t333 = t223 * t97;
t86 = -pkin(3) * t203 + t96;
t40 = t222 * t86 + t333;
t45 = t222 * t96 + t333;
t337 = t248 * t45;
t34 = -t203 * qJ(5) + t40;
t25 = t34 + t356;
t332 = t224 * t25;
t331 = pkin(4) * t301 + t330;
t328 = t137 * t225;
t327 = t173 * t203;
t326 = t175 * t203;
t325 = t187 * t225;
t324 = t187 * t228;
t323 = t203 * t228;
t320 = t225 * t226;
t231 = qJD(1) ^ 2;
t316 = t229 * t231;
t230 = qJD(2) ^ 2;
t315 = t230 * t226;
t314 = t230 * t229;
t313 = -t312 + t350;
t169 = t228 * t183;
t126 = -qJ(4) * t318 + t169 + (-pkin(3) - t344) * t229;
t304 = t225 * t183 + t206;
t131 = -qJ(4) * t320 + t304;
t75 = t222 * t126 + t223 * t131;
t185 = t343 * t225;
t186 = t343 * t228;
t134 = t222 * t185 - t223 * t186;
t303 = pkin(3) * t320 + t226 * pkin(7);
t220 = t226 ^ 2;
t302 = -t229 ^ 2 + t220;
t39 = t223 * t86 - t334;
t260 = qJD(5) - t39;
t22 = t345 * t203 + t260 - t350;
t287 = t22 * t292 + t224 * t7 + t227 * t6;
t283 = t239 * pkin(3) + pkin(7) * t297;
t282 = -pkin(3) * t228 - pkin(2);
t209 = -pkin(3) * t223 - pkin(4);
t277 = t229 * t291;
t274 = t203 * t295;
t119 = pkin(3) * t138 + pkin(7) * t272;
t23 = -t222 * t73 + t223 * t58;
t268 = pkin(1) * t351;
t74 = t126 * t223 - t222 * t131;
t133 = -t223 * t185 - t186 * t222;
t267 = t290 ^ 2;
t265 = t173 + t291;
t264 = -t175 + t299;
t261 = -t104 * t121 + t133 * t85 - t134 * t84;
t106 = -pkin(9) * t167 + t133;
t259 = t311 * pkin(9) - qJD(6) * t106 - t329;
t107 = pkin(9) * t166 + t134;
t258 = -t310 * pkin(9) - qJD(1) * t284 + qJD(6) * t107 - t330;
t70 = -qJ(5) * t229 + t75;
t148 = -t222 * t320 + t223 * t318;
t255 = qJ(5) * t148 - t303;
t72 = t229 * pkin(4) - t74;
t254 = -pkin(3) * t175 - qJ(5) * t121;
t2 = t224 * t22 + t227 * t25;
t47 = pkin(5) * t229 - pkin(9) * t148 + t72;
t147 = t167 * t226;
t49 = pkin(9) * t147 + t70;
t253 = -t224 * t49 + t227 * t47;
t252 = t224 * t47 + t227 * t49;
t251 = -t121 ^ 2 - t352;
t250 = t227 * t147 - t148 * t224;
t100 = t147 * t224 + t148 * t227;
t201 = -pkin(5) + t209;
t207 = pkin(3) * t222 + qJ(5);
t247 = t201 * t227 - t207 * t224;
t246 = t201 * t224 + t207 * t227;
t245 = qJD(1) * t220 - t203 * t229;
t19 = qJ(5) * t298 - qJD(5) * t229 + t24;
t244 = qJ(5) * t167 - t282;
t241 = -t25 * t293 + t287;
t10 = -pkin(4) * t210 + t342;
t102 = t157 * t226 + t222 * t278 - t223 * t277;
t238 = -qJ(5) * t102 + qJD(5) * t148 - t283;
t237 = qJ(5) * t85 + qJD(5) * t248 - t119;
t17 = pkin(4) * t84 - t237;
t115 = pkin(4) * t166 - t244;
t101 = -t222 * t277 - t239 * t223 + t226 * t273;
t88 = -t345 * t166 + t244;
t87 = pkin(4) * t147 - t255;
t71 = -t345 * t147 + t255;
t50 = pkin(4) * t248 - t254;
t33 = pkin(4) * t203 + t260;
t32 = -t248 * t345 + t254;
t31 = -pkin(4) * t101 - t238;
t28 = t45 + t356;
t27 = qJD(6) * t100 + t227 * t101 - t102 * t224;
t26 = qJD(6) * t250 - t101 * t224 - t102 * t227;
t21 = -pkin(4) * t298 - t23;
t18 = t345 * t101 + t238;
t14 = -pkin(9) * t101 + t19;
t13 = pkin(9) * t102 - t23 - t262;
t8 = -t345 * t84 + t237;
t1 = t22 * t227 - t332;
t3 = [0, 0, 0, 0.2e1 * t229 * t210, t302 * t351, t314, -t315, 0, -pkin(7) * t314 + t226 * t268, pkin(7) * t315 + t229 * t268, t137 * t318 + (-t276 + t277) * t175 (-t173 * t228 - t175 * t225) * t297 + (-t328 - t138 * t228 + (t173 * t225 - t175 * t228) * qJD(3)) * t226, t203 * t276 - t137 * t229 + (t175 * t226 + t228 * t245) * qJD(2), t226 * t274 + t138 * t229 + (-t173 * t226 - t225 * t245) * qJD(2) (-t203 - t300) * t298 -(-t183 * t296 + t307) * t203 + (t187 * t295 + pkin(7) * t138 + (qJD(1) * t169 + t129) * qJD(2)) * t226 + ((pkin(7) * t173 + t325) * qJD(2) + (t321 + (pkin(7) * t203 + t188) * t228) * qJD(3) + t309) * t229 (-pkin(7) * t229 * t296 + t308) * t203 + t243 * t229 + (pkin(7) * t137 - t187 * t296) * t226 + ((pkin(7) * t175 + t324) * t229 + (-pkin(7) * t323 - qJD(1) * t304 - t130) * t226) * qJD(2), t101 * t40 + t102 * t39 - t12 * t147 - t121 * t24 + t148 * t342 - t23 * t248 - t74 * t85 - t75 * t84, t119 * t303 + t12 * t75 + t39 * t23 + t40 * t24 - t240 * t283 - t342 * t74, t10 * t229 - t101 * t48 + t121 * t31 + t147 * t17 + t203 * t21 + t84 * t87 + (-qJD(1) * t72 - t33) * t298, t10 * t148 + t101 * t34 - t102 * t33 - t121 * t19 - t147 * t9 + t21 * t248 - t70 * t84 + t72 * t85, t102 * t48 - t248 * t31 - t148 * t17 - t19 * t203 - t229 * t9 - t85 * t87 + (qJD(1) * t70 + t34) * t298, t10 * t72 + t17 * t87 + t19 * t34 + t21 * t33 + t31 * t48 + t70 * t9, -t100 * t266 + t26 * t69, -t100 * t16 - t250 * t266 - t26 * t353 - t27 * t69, -t266 * t229 - t290 * t26 + (-qJD(1) * t100 - t69) * t298, -t16 * t229 + t290 * t27 + (-qJD(1) * t250 + t353) * t298 (t290 - t300) * t298 -(t13 * t227 - t14 * t224) * t290 + t281 * t229 + t18 * t353 + t71 * t16 - t8 * t250 + t30 * t27 + (-t2 * t229 + t252 * t290) * qJD(6) + (-qJD(1) * t253 - t1) * t298 (qJD(6) * t253 + t13 * t224 + t14 * t227) * t290 - t241 * t229 + t18 * t69 - t71 * t266 + t8 * t100 + t30 * t26 + (qJD(1) * t252 + t2) * t298; 0, 0, 0, -t226 * t316, t302 * t231, 0, 0, 0, t231 * pkin(1) * t226, pkin(1) * t316, -t175 * t323 + t328 (t137 + t327) * t228 + (-t138 + t326) * t225, -t274 + (t203 * t317 + t226 * t264) * qJD(1), t203 * t296 + (-t203 * t319 + t226 * t265) * qJD(1), t203 * t301, -pkin(2) * t138 + t306 * t203 + (pkin(8) * t323 + t325) * qJD(3) + ((-pkin(8) * t299 - t129) * t226 + (-pkin(7) * t265 - t325) * t229) * qJD(1), -pkin(2) * t137 - t159 * t203 + (-pkin(8) * t203 * t225 + t324) * qJD(3) + (-t187 * t317 + (-pkin(8) * t291 + t130) * t226 + (t203 * t318 + t229 * t264) * pkin(7)) * qJD(1), -t12 * t166 + t121 * t57 + t167 * t342 + t248 * t330 + t310 * t39 + t311 * t40 + t261, t12 * t134 + t342 * t133 + t119 * t282 + (t104 - t57) * t40 - t330 * t39 - t354 * t240, t115 * t84 + t166 * t17 - t311 * t48 + t331 * t203 - t338 * t121 + (-qJD(2) * t133 + t33) * t301, t10 * t167 + t121 * t51 - t166 * t9 + t248 * t331 - t310 * t33 + t311 * t34 + t261, -t115 * t85 - t167 * t17 + t310 * t48 - t329 * t203 + t338 * t248 + (qJD(2) * t134 - t34) * t301, t10 * t133 + t115 * t17 + t134 * t9 + t329 * t34 + t33 * t331 - t338 * t48, -t118 * t266 + t341 * t69, -t118 * t16 - t249 * t266 - t340 * t69 - t341 * t353, -t341 * t290 + (-qJD(2) * t118 + t69) * t301, t340 * t290 + (-qJD(2) * t249 - t353) * t301, -t290 * t301, -t8 * t249 + t88 * t16 + t339 * t353 + t340 * t30 - (t224 * t259 - t227 * t258) * t290 + (-(t106 * t227 - t107 * t224) * qJD(2) + t1) * t301, t8 * t118 - t88 * t266 + t339 * t69 + t341 * t30 - (t224 * t258 + t227 * t259) * t290 + ((t106 * t224 + t107 * t227) * qJD(2) - t2) * t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175 * t173, -t173 ^ 2 + t175 ^ 2, t137 - t327, -t138 - t326, t210, -t130 * t203 - t175 * t187 + t235, -t129 * t203 + t173 * t187 - t232, t248 * t40 - t337 + (-t222 * t84 - t223 * t85) * pkin(3) + (t46 - t39) * t121, t39 * t45 - t40 * t46 + (t12 * t222 + t175 * t240 - t223 * t342) * pkin(3), -t349 - t121 * t50 - t203 * t45 + (pkin(4) - t209) * t210 - t342, -t207 * t84 + t209 * t85 + t248 * t34 - t337 + (-t312 + t33) * t121, -t121 * t48 + t203 * t46 + t207 * t210 + t248 * t50 - 0.2e1 * t190 + t286, t10 * t209 + t207 * t9 + t312 * t34 - t33 * t45 - t48 * t50, -t357, -t361, t362, t363, t210, -t247 * t210 - t32 * t353 - (t224 * t313 - t227 * t28) * t290 + (t246 * t290 + t2) * qJD(6) - t360, t246 * t210 - t32 * t69 - t358 - (t224 * t28 + t227 * t313) * t290 + (t247 * t290 - t332) * qJD(6) + t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t121 * t40 + t248 * t39 + t119, -t203 * t248 + t84, t251, -t85 - t355, t121 * t34 - t248 * t33 + t17, 0, 0, 0, 0, 0, -t16 + t336, t266 - t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t248 - t210, t85 - t355, -t203 ^ 2 - t352, t203 * t34 + t10 + t349, 0, 0, 0, 0, 0, -t210 * t227 - t224 * t267 - t248 * t353, t210 * t224 - t227 * t267 - t248 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t357, t361, -t362, -t363, -t210 (-qJD(6) - t290) * t2 + t360, -t1 * t290 - t241 + t358;];
tauc_reg  = t3;
