% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:10
% EndTime: 2019-03-09 18:23:30
% DurationCPUTime: 9.18s
% Computational Cost: add. (15132->541), mult. (35088->685), div. (0->0), fcn. (25088->8), ass. (0->290)
t411 = qJD(5) + qJD(6);
t261 = sin(qJ(6));
t262 = sin(qJ(5));
t265 = cos(qJ(5));
t403 = cos(qJ(6));
t412 = -t261 * t265 - t403 * t262;
t167 = t411 * t412;
t263 = sin(qJ(3));
t264 = sin(qJ(2));
t266 = cos(qJ(2));
t404 = cos(qJ(3));
t223 = t263 * t266 + t264 * t404;
t347 = qJD(1) * t223;
t435 = t412 * t347 + t167;
t329 = t404 * t266;
t307 = qJD(1) * t329;
t346 = qJD(1) * t264;
t325 = t263 * t346;
t203 = -t307 + t325;
t257 = qJD(2) + qJD(3);
t175 = -t265 * t203 + t257 * t262;
t177 = t203 * t262 + t257 * t265;
t106 = t403 * t175 + t177 * t261;
t292 = -t261 * t175 + t177 * t403;
t170 = t257 * t223;
t159 = t170 * qJD(1);
t343 = qJD(5) * t265;
t344 = qJD(5) * t262;
t351 = -t203 * t344 - t257 * t343;
t304 = t265 * t159 + t351;
t149 = t262 * t159;
t95 = -t203 * t343 + t257 * t344 - t149;
t35 = qJD(6) * t292 - t261 * t95 - t403 * t304;
t222 = -t261 * t262 + t403 * t265;
t321 = t403 * qJD(6);
t341 = qJD(6) * t261;
t423 = -t261 * t344 - t262 * t341 + (qJD(5) * t403 + t321) * t265 + t222 * t347;
t12 = t106 * t423 - t412 * t35;
t362 = t263 * t264;
t302 = t257 * t362;
t350 = t257 * t307;
t158 = qJD(1) * t302 - t350;
t422 = qJD(5) + t347;
t192 = qJD(6) + t422;
t434 = -t412 * t158 - t192 * t423;
t34 = t175 * t321 + t177 * t341 - t261 * t304 + t403 * t95;
t11 = -t222 * t34 + t292 * t435;
t405 = -pkin(8) - pkin(7);
t231 = t405 * t264;
t225 = qJD(1) * t231;
t393 = qJD(2) * pkin(2);
t211 = t225 + t393;
t193 = t404 * t211;
t232 = t405 * t266;
t226 = qJD(1) * t232;
t363 = t263 * t226;
t162 = -t193 - t363;
t339 = qJD(4) + t162;
t400 = t347 * pkin(4);
t340 = t400 + t339;
t406 = pkin(3) + pkin(9);
t101 = -t257 * t406 + t340;
t249 = -t266 * pkin(2) - pkin(1);
t230 = qJD(1) * t249;
t278 = -qJ(4) * t347 + t230;
t104 = t203 * t406 + t278;
t57 = t265 * t101 - t104 * t262;
t41 = -pkin(10) * t177 + t57;
t39 = pkin(5) * t422 + t41;
t58 = t101 * t262 + t104 * t265;
t42 = -pkin(10) * t175 + t58;
t390 = t261 * t42;
t19 = t39 * t403 - t390;
t335 = t403 * t42;
t20 = t261 * t39 + t335;
t337 = qJD(1) * qJD(2);
t320 = t264 * t337;
t288 = pkin(2) * t320 + t158 * qJ(4) - qJD(4) * t347;
t54 = t159 * t406 + t288;
t331 = qJD(2) * t405;
t308 = qJD(1) * t331;
t214 = t264 * t308;
t215 = t266 * t308;
t322 = qJD(3) * t404;
t345 = qJD(3) * t263;
t90 = t211 * t345 + t263 * t214 - t404 * t215 - t226 * t322;
t75 = -pkin(4) * t158 + t90;
t14 = t101 * t343 - t104 * t344 + t262 * t75 + t265 * t54;
t10 = pkin(10) * t304 + t14;
t15 = -qJD(5) * t58 - t262 * t54 + t265 * t75;
t9 = -t158 * pkin(5) + t95 * pkin(10) + t15;
t273 = -t10 * t403 - t261 * t9 - t321 * t39 + t42 * t341;
t4 = -qJD(6) * t20 - t261 * t10 + t403 * t9;
t276 = -t19 * t435 - t423 * t20 - t4 * t222 - t273 * t412;
t313 = t422 * t177;
t432 = t304 + t313;
t431 = -t222 * t158 + t435 * t192;
t312 = t422 * t175;
t430 = t95 - t312;
t429 = -t203 * pkin(5) - pkin(10) * t344;
t409 = t422 * t58 + t15;
t253 = t264 * t393;
t428 = 0.2e1 * t253;
t160 = pkin(3) * t347 + qJ(4) * t203;
t145 = pkin(2) * t346 + t160;
t199 = t347 * pkin(9);
t112 = t145 + t199;
t209 = t404 * t226;
t164 = t225 * t263 - t209;
t401 = t203 * pkin(4);
t134 = t164 - t401;
t125 = t265 * t134;
t248 = -pkin(2) * t404 - pkin(3);
t244 = -pkin(9) + t248;
t336 = pkin(2) * t345;
t275 = -t244 * t344 + t265 * t336;
t402 = pkin(10) * t347;
t427 = -t275 + t125 + (-t112 - t402) * t262 + t429;
t373 = t347 * t265;
t187 = pkin(10) * t373;
t399 = -pkin(10) + t244;
t213 = t399 * t265;
t71 = t265 * t112 + t262 * t134;
t426 = qJD(5) * t213 + t262 * t336 - t187 - t71;
t163 = t263 * t211 - t209;
t127 = t163 - t401;
t121 = t265 * t127;
t122 = t160 + t199;
t342 = qJD(5) * t406;
t326 = t262 * t342;
t425 = t326 - t121 - (-t122 - t402) * t262 - t429;
t398 = -pkin(10) - t406;
t228 = t398 * t265;
t74 = t265 * t122 + t262 * t127;
t424 = -qJD(5) * t228 + t187 + t74;
t383 = t106 * t292;
t221 = -t329 + t362;
t151 = t222 * t221;
t421 = -t106 ^ 2 + t292 ^ 2;
t420 = t106 * t192 - t34;
t255 = t257 * qJ(4);
t113 = t255 + t127;
t86 = t175 * pkin(5) + t113;
t419 = t86 * t106 + t273;
t417 = -0.2e1 * t337;
t295 = -t223 * qJ(4) + t249;
t131 = t221 * t406 + t295;
t180 = -t404 * t231 - t232 * t263;
t146 = pkin(4) * t223 + t180;
t142 = t262 * t146;
t80 = t265 * t131 + t142;
t128 = t158 * t223;
t323 = qJD(2) * t404;
t169 = t302 + (-t322 - t323) * t266;
t414 = -t169 * t347 - t128;
t330 = t404 * t225;
t165 = t330 + t363;
t236 = pkin(2) * t322 + qJD(4);
t352 = t165 - t236;
t332 = -pkin(5) * t265 - pkin(4);
t413 = pkin(5) * t343 - t332 * t347 - t363;
t410 = -t86 * t292 + t4;
t408 = t192 * t292 - t35;
t407 = t347 ^ 2;
t212 = t399 * t262;
t154 = t212 * t403 + t261 * t213;
t397 = qJD(6) * t154 + t426 * t261 + t427 * t403;
t153 = -t261 * t212 + t213 * t403;
t396 = -qJD(6) * t153 + t427 * t261 - t426 * t403;
t227 = t398 * t262;
t173 = -t261 * t227 + t228 * t403;
t395 = -qJD(6) * t173 - t425 * t261 + t424 * t403;
t389 = t262 * t57;
t55 = t57 * t344;
t394 = t347 * t389 + t55;
t388 = t265 * t95;
t387 = t90 * t180;
t174 = t227 * t403 + t261 * t228;
t386 = -qJD(6) * t174 + t424 * t261 + t425 * t403;
t385 = t236 - t330 + t413;
t384 = qJD(4) - t193 + t413;
t296 = t405 * t323;
t310 = t263 * t331;
t110 = -t231 * t322 - t232 * t345 - t264 * t296 - t266 * t310;
t382 = t110 * t257;
t181 = t263 * t231 - t232 * t404;
t111 = qJD(3) * t181 + t264 * t310 - t266 * t296;
t381 = t111 * t257;
t380 = t163 * t257;
t378 = t177 * t175;
t377 = t192 * t203;
t376 = t422 * t203;
t375 = t203 * t347;
t371 = t221 * t262;
t370 = t221 * t265;
t369 = t257 * t169;
t368 = t257 * t170;
t367 = t257 * t203;
t364 = t262 * t158;
t150 = t265 * t158;
t269 = qJD(1) ^ 2;
t361 = t266 * t269;
t360 = t406 * t158;
t268 = qJD(2) ^ 2;
t359 = t268 * t264;
t358 = t268 * t266;
t349 = t400 - t352;
t348 = t264 ^ 2 - t266 ^ 2;
t333 = t264 * t361;
t324 = -t203 ^ 2 + t407;
t245 = pkin(2) * t263 + qJ(4);
t319 = -pkin(10) * t221 - t131;
t311 = -t211 * t322 - t404 * t214 - t263 * t215 - t226 * t345;
t87 = -t257 * qJD(4) + t311;
t65 = -t159 * pkin(4) - t87;
t318 = -t58 * t203 + t65 * t265;
t316 = pkin(1) * t417;
t315 = t422 ^ 2;
t314 = t422 * t113;
t306 = t266 * t320;
t305 = -t164 + t336;
t303 = qJD(5) * t177 + t351;
t301 = t262 * t58 + t265 * t57;
t300 = t265 * t58 - t389;
t298 = t159 * t221 + t170 * t203;
t297 = t57 * t203 + t65 * t262 + (t343 + t373) * t113;
t143 = t265 * t146;
t60 = t223 * pkin(5) + t262 * t319 + t143;
t72 = pkin(10) * t370 + t80;
t30 = -t261 * t72 + t403 * t60;
t31 = t261 * t60 + t403 * t72;
t294 = t304 * t262;
t293 = -qJD(5) * t175 + t149 + t95;
t291 = -t262 * t315 - t150;
t290 = t170 * t262 + t221 * t343;
t289 = -t170 * t265 + t221 * t344;
t287 = t169 * qJ(4) - t223 * qJD(4) + t253;
t64 = t170 * t406 + t287;
t84 = -t169 * pkin(4) + t111;
t23 = -t131 * t344 + t146 * t343 + t262 * t84 + t265 * t64;
t38 = -pkin(5) * t304 + t65;
t286 = t19 * t203 - t38 * t412 + t423 * t86;
t285 = -t20 * t203 + t38 * t222 + t435 * t86;
t137 = t203 * pkin(3) + t278;
t284 = t137 * t347 + t90;
t283 = -t230 * t347 - t90;
t282 = t230 * t203 + t311;
t281 = -t137 * t203 - t87;
t277 = -t265 * t315 + t364;
t83 = -t170 * pkin(4) - t110;
t274 = t158 * t221 - t159 * t223 + t169 * t203 - t170 * t347;
t116 = -t257 * t325 + t350 + t367;
t271 = qJD(5) * t300 + t14 * t262 + t15 * t265;
t270 = t110 * t203 + t111 * t347 - t158 * t180 - t159 * t181 + t223 * t90;
t256 = t262 * pkin(5);
t246 = qJ(4) + t256;
t229 = t245 + t256;
t161 = pkin(3) * t221 + t295;
t152 = t412 * t221;
t148 = -t255 - t163;
t147 = -t221 * pkin(4) + t181;
t144 = -pkin(3) * t257 + t339;
t117 = -t257 * t347 + t159;
t115 = t158 - t367;
t114 = t221 * t332 + t181;
t85 = pkin(3) * t170 + t287;
t82 = t265 * t84;
t79 = -t131 * t262 + t143;
t76 = pkin(3) * t159 + t288;
t73 = -t122 * t262 + t121;
t70 = -t112 * t262 + t125;
t53 = -t175 * t203 + t277;
t52 = t177 * t203 + t291;
t50 = pkin(5) * t289 + t83;
t48 = t265 * t312 - t294;
t47 = -t262 * t313 - t388;
t44 = -t167 * t221 - t170 * t222;
t43 = -t151 * t411 + t170 * t412;
t33 = -t106 * t203 + t434;
t32 = t203 * t292 + t431;
t25 = (-t313 + t304) * t265 + (t95 + t312) * t262;
t24 = -qJD(5) * t80 - t262 * t64 + t82;
t22 = t403 * t41 - t390;
t21 = -t261 * t41 - t335;
t18 = -pkin(10) * t289 + t23;
t13 = -t169 * pkin(5) + t82 + (-pkin(10) * t170 - t64) * t262 + (t265 * t319 - t142) * qJD(5);
t7 = -t106 * t435 - t222 * t35 - t292 * t423 - t34 * t412;
t6 = -qJD(6) * t31 + t13 * t403 - t261 * t18;
t5 = qJD(6) * t30 + t261 * t13 + t18 * t403;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t306, t348 * t417, t358, -0.2e1 * t306, -t359, 0, -pkin(7) * t358 + t264 * t316, pkin(7) * t359 + t266 * t316, 0, 0, t414, t274, -t369, t298, -t368, 0, -t381 + t249 * t159 + t230 * t170 + (qJD(1) * t221 + t203) * t253, -t249 * t158 - t230 * t169 + t347 * t428 + t382, -t162 * t169 - t163 * t170 + t221 * t311 + t270, -t163 * t110 + t162 * t111 - t181 * t311 + t230 * t428 + t387, 0, t369, t368, t414, t274, t298, -t144 * t169 + t148 * t170 + t221 * t87 + t270, -t137 * t170 - t159 * t161 - t203 * t85 - t221 * t76 + t381, t137 * t169 + t158 * t161 - t223 * t76 - t347 * t85 - t382, t110 * t148 + t111 * t144 + t137 * t85 + t161 * t76 - t181 * t87 + t387, t177 * t290 - t371 * t95 (-t175 * t262 + t177 * t265) * t170 + (t294 - t388 + (-t175 * t265 - t177 * t262) * qJD(5)) * t221, -t177 * t169 - t221 * t364 - t95 * t223 + t290 * t422, t175 * t289 + t304 * t370, -t150 * t221 + t175 * t169 + t223 * t304 - t289 * t422, -t169 * t422 - t128, t113 * t289 - t147 * t304 + t15 * t223 - t79 * t158 - t57 * t169 + t83 * t175 + t24 * t422 - t370 * t65, t113 * t290 - t14 * t223 - t147 * t95 + t80 * t158 + t58 * t169 + t83 * t177 - t23 * t422 + t371 * t65, -t23 * t175 + t80 * t304 - t24 * t177 + t79 * t95 + t300 * t170 + (-qJD(5) * t301 + t14 * t265 - t15 * t262) * t221, t113 * t83 + t14 * t80 + t147 * t65 + t15 * t79 + t23 * t58 + t24 * t57, t152 * t34 - t292 * t43, t106 * t43 - t151 * t34 + t152 * t35 - t292 * t44, t152 * t158 - t169 * t292 - t192 * t43 - t223 * t34, t106 * t44 - t151 * t35, t106 * t169 - t151 * t158 - t192 * t44 - t223 * t35, -t169 * t192 - t128, t106 * t50 + t114 * t35 - t151 * t38 - t158 * t30 - t169 * t19 + t192 * t6 + t223 * t4 + t44 * t86, -t114 * t34 - t152 * t38 + t158 * t31 + t169 * t20 - t192 * t5 + t223 * t273 + t292 * t50 - t43 * t86, -t106 * t5 - t151 * t273 + t152 * t4 + t19 * t43 - t20 * t44 - t292 * t6 + t30 * t34 - t31 * t35, t114 * t38 + t19 * t6 + t20 * t5 - t273 * t31 + t30 * t4 + t50 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t333, t348 * t269, 0, t333, 0, 0, t269 * pkin(1) * t264, pkin(1) * t361, 0, 0, t375, t324, t116, -t375, 0, 0, t164 * t257 + (-t203 * t346 - t257 * t345) * pkin(2) + t283, t165 * t257 + (-t257 * t322 - t346 * t347) * pkin(2) + t282 (t163 - t164) * t347 + (t162 + t165) * t203 + (t404 * t158 - t159 * t263 + (-t203 * t404 + t263 * t347) * qJD(3)) * pkin(2), -t162 * t164 - t163 * t165 + (-t230 * t346 - t404 * t90 - t263 * t311 + (t162 * t263 + t163 * t404) * qJD(3)) * pkin(2), 0, t115, t117, t375, t324, -t375, -t248 * t158 - t245 * t159 + (-t148 + t305) * t347 + (t144 + t352) * t203, t145 * t203 + t257 * t305 + t284, t145 * t347 - t257 * t352 + t281, -t137 * t145 + t144 * t305 + t148 * t352 - t87 * t245 + t90 * t248, t47, t25, t52, t48, t53, t376, -t244 * t150 - t245 * t304 + t349 * t175 + (t275 - t70) * t422 + t297, -t245 * t95 + (-t244 * t343 + t71) * t422 + t349 * t177 + (t244 * t158 - t336 * t422 - t314) * t262 + t318, t71 * t175 + t70 * t177 + (-t175 * t336 + t244 * t303 - t14) * t262 + (-t177 * t336 + t244 * t293 - t409) * t265 + t394, t113 * t349 + t244 * t271 + t65 * t245 + t301 * t336 - t57 * t70 - t58 * t71, t11, t7, t32, t12, t33, t377, t106 * t385 - t153 * t158 - t192 * t397 + t229 * t35 + t286, t154 * t158 + t192 * t396 - t229 * t34 + t292 * t385 + t285, t106 * t396 + t153 * t34 - t154 * t35 + t292 * t397 + t276, t4 * t153 - t154 * t273 - t19 * t397 - t20 * t396 + t38 * t229 + t385 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, t324, t116, -t375, 0, 0, t283 + t380, -t162 * t257 + t282, 0, 0, 0, t115, t117, t375, t324, -t375, pkin(3) * t158 - qJ(4) * t159 + (-t148 - t163) * t347 + (t144 - t339) * t203, t160 * t203 + t284 - t380, t160 * t347 + t257 * t339 + t281, -t90 * pkin(3) - t87 * qJ(4) - t137 * t160 - t144 * t163 - t148 * t339, t47, t25, t52, t48, t53, t376, t265 * t360 - qJ(4) * t304 + (-t73 + t326) * t422 + t340 * t175 + t297, -qJ(4) * t95 + (t265 * t342 + t74) * t422 + t340 * t177 + (-t314 - t360) * t262 + t318, t74 * t175 + t73 * t177 + (-t303 * t406 - t14) * t262 + (-t293 * t406 - t409) * t265 + t394, t65 * qJ(4) + t113 * t340 - t271 * t406 - t57 * t73 - t58 * t74, t11, t7, t32, t12, t33, t377, t106 * t384 - t173 * t158 + t192 * t386 + t246 * t35 + t286, t174 * t158 + t192 * t395 - t246 * t34 + t292 * t384 + t285, t106 * t395 + t173 * t34 - t174 * t35 - t292 * t386 + t276, t4 * t173 - t174 * t273 + t19 * t386 - t20 * t395 + t38 * t246 + t384 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, -t375, -t257 ^ 2 - t407, t148 * t257 + t284, 0, 0, 0, 0, 0, 0, -t257 * t175 + t291, -t257 * t177 + t277, t432 * t262 + t430 * t265, -t113 * t257 - t55 + (-t347 * t57 + t14) * t262 + t409 * t265, 0, 0, 0, 0, 0, 0, -t106 * t257 + t431, -t257 * t292 + t434, -t11 - t12, -t86 * t257 - t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t378, -t175 ^ 2 + t177 ^ 2, -t430, -t378, t432, -t158, -t113 * t177 + t409, t113 * t175 + t422 * t57 - t14, 0, 0, t383, t421, t420, -t383, t408, -t158, -t21 * t192 + (-t106 * t177 - t158 * t403 - t192 * t341) * pkin(5) + t410, t22 * t192 + (t158 * t261 - t177 * t292 - t192 * t321) * pkin(5) + t419, t20 * t292 + t22 * t106 - t19 * t106 + t21 * t292 + (t403 * t34 - t261 * t35 + (-t106 * t403 + t261 * t292) * qJD(6)) * pkin(5), -t19 * t21 - t20 * t22 + (t403 * t4 - t177 * t86 - t261 * t273 + (-t19 * t261 + t20 * t403) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t383, t421, t420, -t383, t408, -t158, t20 * t192 + t410, t19 * t192 + t419, 0, 0;];
tauc_reg  = t1;