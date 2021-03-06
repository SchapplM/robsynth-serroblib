% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:44:38
% EndTime: 2019-03-08 23:45:10
% DurationCPUTime: 15.73s
% Computational Cost: add. (16097->659), mult. (43761->947), div. (0->0), fcn. (36198->14), ass. (0->288)
t276 = cos(qJ(3));
t393 = cos(pkin(7));
t330 = t276 * t393;
t267 = sin(pkin(7));
t273 = sin(qJ(3));
t378 = t267 * t273;
t225 = pkin(2) * t330 - pkin(9) * t378;
t216 = qJD(3) * t225;
t274 = sin(qJ(2));
t277 = cos(qJ(2));
t332 = t273 * t393;
t290 = -t274 * t332 + t276 * t277;
t268 = sin(pkin(6));
t361 = qJD(1) * t268;
t363 = t290 * t361 - t216;
t309 = pkin(3) * t273 - pkin(10) * t276;
t298 = t309 * qJD(3);
t343 = t274 * t361;
t440 = (t298 - t343) * t267;
t272 = sin(qJ(4));
t275 = cos(qJ(4));
t327 = t393 * qJD(2);
t303 = t327 + qJD(3);
t359 = qJD(2) * t267;
t342 = t273 * t359;
t201 = t272 * t303 + t275 * t342;
t357 = qJD(2) * t276;
t341 = t267 * t357;
t249 = -qJD(4) + t341;
t266 = sin(pkin(13));
t269 = cos(pkin(13));
t161 = t201 * t269 - t249 * t266;
t271 = sin(qJ(6));
t184 = t266 * t201;
t325 = t249 * t269 + t184;
t416 = cos(qJ(6));
t297 = t416 * t325;
t86 = t161 * t271 + t297;
t439 = t86 ^ 2;
t242 = t275 * t303;
t199 = t272 * t342 - t242;
t193 = qJD(6) + t199;
t438 = t193 * t86;
t377 = t267 * t276;
t226 = pkin(2) * t332 + pkin(9) * t377;
t210 = pkin(10) * t393 + t226;
t310 = -pkin(3) * t276 - pkin(10) * t273;
t211 = (-pkin(2) + t310) * t267;
t351 = qJD(4) * t275;
t352 = qJD(4) * t272;
t399 = t210 * t352 - t211 * t351 - t440 * t272 + t275 * t363;
t217 = qJD(3) * t226;
t288 = t273 * t277 + t274 * t330;
t364 = t288 * t361 - t217;
t88 = -t161 * t416 + t271 * t325;
t437 = t88 ^ 2;
t355 = qJD(3) * t273;
t436 = -(qJ(5) * t355 - qJD(5) * t276) * t267 + t399;
t346 = t272 * t378;
t318 = qJD(4) * t346;
t331 = t275 * t393;
t356 = qJD(3) * t267;
t339 = t276 * t356;
t178 = -qJD(4) * t331 - t275 * t339 + t318;
t224 = t272 * t393 + t275 * t378;
t179 = qJD(4) * t224 + t272 * t339;
t435 = -t179 * pkin(4) - t178 * qJ(5) + t224 * qJD(5) + t364;
t409 = qJD(2) * pkin(2);
t239 = t277 * t361 + t409;
t328 = t393 * t239;
t228 = pkin(9) * t359 + t343;
t370 = t276 * t228;
t291 = t273 * t328 + t370;
t307 = pkin(4) * t272 - qJ(5) * t275;
t270 = cos(pkin(6));
t360 = qJD(1) * t270;
t434 = -(t273 * t360 + t307 * t357) * t267 - t291 + qJD(4) * t307 - qJD(5) * t272;
t410 = t266 * t436 - t269 * t435;
t347 = pkin(10) * t352;
t218 = t273 * t228;
t287 = t267 * t360 + t328;
t141 = t276 * t287 - t218;
t214 = t309 * t359;
t104 = t275 * t141 + t272 * t214;
t93 = qJ(5) * t342 + t104;
t398 = t434 * t269 + (t347 + t93) * t266;
t433 = t266 * t434 - t269 * t93;
t395 = t266 * t435 + t269 * t436;
t340 = t267 * t355;
t151 = -t178 * t269 + t266 * t340;
t432 = pkin(5) * t179 - pkin(11) * t151 + t410;
t150 = -t178 * t266 - t269 * t340;
t431 = pkin(11) * t150 + t395;
t371 = t275 * t276;
t186 = (t266 * t273 + t269 * t371) * t359;
t373 = t269 * t275;
t415 = pkin(5) * t272;
t430 = -pkin(11) * t186 + t341 * t415 - (-pkin(11) * t373 + t415) * qJD(4) - t398;
t400 = -t210 * t351 - t211 * t352 + t363 * t272 + t275 * t440;
t379 = t266 * t275;
t185 = -t269 * t342 + t341 * t379;
t374 = t269 * t272;
t429 = pkin(11) * t185 + (-pkin(10) * t374 - pkin(11) * t379) * qJD(4) + t433;
t426 = -t266 * t351 + t185;
t425 = -t269 * t351 + t186;
t142 = t273 * t287 + t370;
t125 = pkin(10) * t303 + t142;
t333 = t270 * t393;
t252 = qJD(1) * t333;
t158 = t252 + (qJD(2) * t310 - t239) * t267;
t171 = (t298 + t343) * t359;
t283 = t290 * qJD(2);
t316 = t270 * t339;
t97 = (t276 * t328 - t218) * qJD(3) + (t268 * t283 + t316) * qJD(1);
t324 = t125 * t351 + t158 * t352 - t275 * t171 + t272 * t97;
t76 = t125 * t275 + t158 * t272;
t424 = t249 * t76 + t324;
t293 = t125 * t352 - t158 * t351 - t272 * t171 - t275 * t97;
t75 = -t272 * t125 + t275 * t158;
t423 = t249 * t75 - t293;
t397 = -t269 * t347 + t433;
t394 = -pkin(4) * t340 - t400;
t349 = qJD(2) * qJD(3);
t334 = t276 * t349;
t313 = t267 * t334;
t353 = qJD(4) * t201;
t165 = t272 * t313 + t353;
t223 = -t331 + t346;
t131 = t165 * t223;
t422 = t179 * t199 + t131;
t329 = t277 * t393;
t421 = -t273 * t274 + t276 * t329;
t335 = qJD(6) * t416;
t350 = qJD(6) * t271;
t420 = -t266 * t350 + t269 * t335;
t314 = qJD(2) * t340;
t284 = qJ(5) * t314 - qJD(5) * t249;
t26 = t284 - t293;
t164 = qJD(2) * t318 - qJD(4) * t242 - t275 * t313;
t282 = t288 * qJD(2);
t317 = t270 * t340;
t98 = t291 * qJD(3) + (t268 * t282 + t317) * qJD(1);
t46 = t165 * pkin(4) + t164 * qJ(5) - t201 * qJD(5) + t98;
t12 = -t26 * t266 + t269 * t46;
t132 = -t164 * t269 + t266 * t314;
t10 = pkin(5) * t165 - pkin(11) * t132 + t12;
t13 = t269 * t26 + t266 * t46;
t159 = t266 * t164;
t292 = t269 * t314 + t159;
t11 = pkin(11) * t292 + t13;
t68 = -qJ(5) * t249 + t76;
t124 = -pkin(3) * t303 - t141;
t77 = t199 * pkin(4) - t201 * qJ(5) + t124;
t28 = -t266 * t68 + t269 * t77;
t21 = pkin(5) * t199 - pkin(11) * t161 + t28;
t29 = t266 * t77 + t269 * t68;
t23 = -pkin(11) * t325 + t29;
t300 = -t21 * t416 + t271 * t23;
t1 = -qJD(6) * t300 + t271 * t10 + t11 * t416;
t419 = t199 ^ 2;
t278 = qJD(2) ^ 2;
t177 = t224 * t269 - t266 * t377;
t209 = -pkin(3) * t393 - t225;
t126 = t223 * pkin(4) - t224 * qJ(5) + t209;
t144 = t275 * t210 + t272 * t211;
t129 = -qJ(5) * t377 + t144;
t65 = t269 * t126 - t129 * t266;
t47 = pkin(5) * t223 - pkin(11) * t177 + t65;
t176 = t224 * t266 + t269 * t377;
t66 = t266 * t126 + t269 * t129;
t59 = -pkin(11) * t176 + t66;
t18 = -t271 * t59 + t416 * t47;
t418 = qJD(6) * t18 + t432 * t271 - t431 * t416;
t19 = t271 * t47 + t416 * t59;
t417 = -qJD(6) * t19 + t431 * t271 + t432 * t416;
t414 = t88 * t86;
t413 = pkin(11) + qJ(5);
t243 = -pkin(4) * t275 - qJ(5) * t272 - pkin(3);
t230 = t269 * t243;
t169 = -pkin(11) * t374 + t230 + (-pkin(10) * t266 - pkin(5)) * t275;
t198 = pkin(10) * t373 + t266 * t243;
t180 = -pkin(11) * t266 * t272 + t198;
t106 = t271 * t169 + t180 * t416;
t412 = qJD(6) * t106 + t429 * t271 + t430 * t416;
t105 = t169 * t416 - t271 * t180;
t411 = -qJD(6) * t105 + t430 * t271 - t429 * t416;
t172 = -t268 * t421 - t270 * t377;
t408 = t172 * t98;
t27 = -pkin(4) * t314 + t324;
t405 = t27 * t266;
t404 = t27 * t269;
t247 = t413 * t266;
t248 = t413 * t269;
t183 = -t271 * t247 + t248 * t416;
t232 = t266 * t416 + t271 * t269;
t385 = t199 * t269;
t140 = pkin(4) * t201 + qJ(5) * t199;
t52 = t269 * t140 - t266 * t75;
t34 = pkin(5) * t201 + pkin(11) * t385 + t52;
t386 = t199 * t266;
t53 = t266 * t140 + t269 * t75;
t42 = pkin(11) * t386 + t53;
t403 = qJD(5) * t232 + qJD(6) * t183 - t271 * t42 + t34 * t416;
t182 = -t247 * t416 - t271 * t248;
t296 = -t271 * t266 + t269 * t416;
t402 = qJD(5) * t296 + qJD(6) * t182 - t271 * t34 - t416 * t42;
t401 = pkin(5) * t150 + t394;
t344 = pkin(5) * t266 + pkin(10);
t103 = -t272 * t141 + t214 * t275;
t96 = -pkin(4) * t342 - t103;
t396 = -pkin(5) * t185 + t344 * t351 - t96;
t392 = t132 * t266;
t391 = t165 * t275;
t192 = -t239 * t267 + t252;
t389 = t192 * t267;
t388 = t199 * t201;
t387 = t199 * t249;
t384 = t201 * t249;
t295 = t249 * t272;
t382 = t249 * t275;
t263 = t267 ^ 2;
t381 = t263 * t278;
t380 = t266 * t165;
t376 = t268 * t278;
t375 = t269 * t165;
t67 = t249 * pkin(4) + qJD(5) - t75;
t369 = -qJD(5) + t67;
t368 = t416 * t185 + t186 * t271 - t232 * t351 - t272 * t420;
t222 = t232 * qJD(6);
t367 = t222 * t272 - t426 * t271 + t425 * t416;
t366 = -t232 * t199 - t222;
t365 = -t296 * t199 - t420;
t362 = t273 ^ 2 - t276 ^ 2;
t358 = qJD(2) * t268;
t354 = qJD(3) * t275;
t348 = t263 * t409;
t345 = t274 * t376;
t326 = t271 * t132 - t416 * t292;
t143 = -t272 * t210 + t211 * t275;
t323 = t263 * t345;
t322 = t273 * t276 * t381;
t319 = t267 * t274 * t358;
t311 = t267 * t278 * t393;
t130 = pkin(4) * t377 - t143;
t289 = t273 * t329 + t274 * t276;
t173 = t268 * t289 + t270 * t378;
t220 = -t268 * t277 * t267 + t333;
t128 = t173 * t275 + t220 * t272;
t127 = t173 * t272 - t220 * t275;
t306 = t325 * qJD(4);
t305 = t263 * t273 * t334;
t302 = 0.2e1 * t327 + qJD(3);
t82 = -t128 * t266 + t172 * t269;
t83 = t128 * t269 + t172 * t266;
t40 = -t271 * t83 + t416 * t82;
t8 = t271 * t21 + t23 * t416;
t41 = t271 * t82 + t416 * t83;
t299 = -t348 + t389;
t109 = -t271 * t176 + t177 * t416;
t35 = qJD(6) * t297 - t416 * t132 + t161 * t350 - t271 * t292;
t286 = -t161 * t266 - t269 * t325;
t285 = t292 * t269;
t281 = qJD(3) * t228 + t327 * t343;
t280 = -t199 * t295 - t391;
t2 = -qJD(6) * t8 + t416 * t10 - t271 * t11;
t36 = -qJD(6) * t88 + t326;
t22 = -pkin(5) * t292 + t27;
t279 = -qJD(3) * t328 - t192 * t359 + (-t270 * t356 - t277 * t358) * qJD(1);
t262 = -pkin(5) * t269 - pkin(4);
t235 = t344 * t272;
t213 = t296 * t272;
t212 = t232 * t272;
t197 = -pkin(10) * t379 + t230;
t113 = t316 + (qJD(3) * t421 + t283) * t268;
t112 = t317 + (qJD(3) * t289 + t282) * t268;
t108 = t176 * t416 + t177 * t271;
t91 = pkin(5) * t176 + t130;
t63 = -pkin(5) * t386 + t76;
t62 = -qJD(4) * t127 + t113 * t275 + t272 * t319;
t61 = qJD(4) * t128 + t113 * t272 - t275 * t319;
t56 = pkin(5) * t325 + t67;
t49 = qJD(6) * t109 + t150 * t416 + t271 * t151;
t48 = t271 * t150 - t151 * t416 + t176 * t335 + t177 * t350;
t38 = t112 * t266 + t269 * t62;
t37 = t112 * t269 - t266 * t62;
t6 = -qJD(6) * t41 - t271 * t38 + t37 * t416;
t5 = qJD(6) * t40 + t271 * t37 + t38 * t416;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345, -t277 * t376, 0, 0, 0, 0, 0, 0, 0, 0, -t112 * t303 + t220 * t314 - t276 * t323, -t113 * t303 + t220 * t313 + t273 * t323 (t112 * t273 + t113 * t276 + (t172 * t276 - t173 * t273) * qJD(3)) * t359, -t112 * t141 + t113 * t142 + t408 + t173 * t97 + (qJD(1) * t220 + t192) * t319, 0, 0, 0, 0, 0, 0, t112 * t199 - t127 * t314 + t165 * t172 + t249 * t61, t112 * t201 - t128 * t314 - t164 * t172 + t249 * t62, -t127 * t164 - t128 * t165 - t199 * t62 + t201 * t61, t112 * t124 + t127 * t324 - t128 * t293 - t61 * t75 + t62 * t76 + t408, 0, 0, 0, 0, 0, 0, -t127 * t292 + t82 * t165 + t37 * t199 + t325 * t61, t127 * t132 + t161 * t61 - t165 * t83 - t199 * t38, -t82 * t132 - t37 * t161 + t292 * t83 - t325 * t38, t12 * t82 + t127 * t27 + t13 * t83 + t28 * t37 + t29 * t38 + t61 * t67, 0, 0, 0, 0, 0, 0, t127 * t36 + t165 * t40 + t193 * t6 + t61 * t86, -t127 * t35 - t165 * t41 - t193 * t5 - t61 * t88, t35 * t40 - t36 * t41 - t5 * t86 + t6 * t88, t1 * t41 + t127 * t22 + t2 * t40 - t300 * t6 + t5 * t8 + t56 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t305, -0.2e1 * t362 * t263 * t349, t302 * t339, -0.2e1 * t305, -t302 * t340, 0 (t273 * t299 + t364) * qJD(3) + (qJD(2) * t364 - t98) * t393 (t276 * t299 + t363) * qJD(3) + (qJD(2) * t363 - t97) * t393 (t273 * t98 + t276 * t97 + (-t141 * t276 - t142 * t273) * qJD(3) + ((-t363 - t216) * t276 + (-t364 - t217) * t273) * qJD(2)) * t267, -t225 * t98 + t226 * t97 - t363 * t142 + t364 * t141 + (-t348 - t389) * t343, -t164 * t224 - t178 * t201, t164 * t223 - t165 * t224 + t178 * t199 - t179 * t201, t178 * t249 + (t164 * t276 + (qJD(2) * t224 + t201) * t355) * t267, t422, t249 * t179 + (t165 * t276 + (-qJD(2) * t223 - t199) * t355) * t267 (-t249 * t267 - t263 * t357) * t355, t124 * t179 + t165 * t209 + t223 * t98 - t400 * t249 - t364 * t199 + (t276 * t324 + (qJD(2) * t143 + t75) * t355) * t267, -t124 * t178 - t164 * t209 + t224 * t98 - t399 * t249 - t364 * t201 + (-t276 * t293 + (-qJD(2) * t144 - t76) * t355) * t267, t143 * t164 - t144 * t165 + t178 * t75 - t179 * t76 + t199 * t399 - t201 * t400 + t223 * t293 + t224 * t324, -t124 * t364 - t143 * t324 - t144 * t293 + t209 * t98 - t399 * t76 + t400 * t75, t132 * t177 + t151 * t161, -t132 * t176 - t161 * t150 - t151 * t325 + t177 * t292, t132 * t223 + t151 * t199 + t161 * t179 + t165 * t177, t150 * t325 - t176 * t292, -t150 * t199 - t176 * t165 - t179 * t325 + t223 * t292, t422, t12 * t223 - t130 * t292 + t67 * t150 + t65 * t165 + t27 * t176 + t28 * t179 + t410 * t199 + t325 * t394, -t13 * t223 + t130 * t132 + t151 * t67 + t161 * t394 - t165 * t66 + t177 * t27 - t179 * t29 + t199 * t395, -t12 * t177 - t13 * t176 - t65 * t132 - t29 * t150 - t28 * t151 - t410 * t161 + t66 * t292 + t325 * t395, t12 * t65 + t13 * t66 + t130 * t27 + t28 * t410 - t29 * t395 + t394 * t67, -t109 * t35 + t48 * t88, t108 * t35 - t109 * t36 + t48 * t86 + t49 * t88, t109 * t165 - t179 * t88 - t193 * t48 - t223 * t35, t108 * t36 + t49 * t86, -t108 * t165 - t179 * t86 - t193 * t49 - t223 * t36, t179 * t193 + t131, t108 * t22 + t165 * t18 - t179 * t300 + t193 * t417 + t2 * t223 + t36 * t91 + t401 * t86 + t49 * t56, -t1 * t223 + t109 * t22 - t165 * t19 - t179 * t8 - t193 * t418 - t35 * t91 - t401 * t88 - t48 * t56, -t1 * t108 - t109 * t2 + t18 * t35 - t19 * t36 - t300 * t48 + t417 * t88 - t418 * t86 - t49 * t8, t1 * t19 + t18 * t2 + t22 * t91 - t300 * t417 + t401 * t56 + t418 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, t362 * t381, -t276 * t311, t322, t273 * t311, 0, t142 * t303 + t273 * t279 - t276 * t281, t141 * t303 + t273 * t281 + t276 * t279, 0, 0, -t164 * t272 - t201 * t382 (-t164 + t387) * t275 + (-t165 + t384) * t272, -t249 * t351 + (t249 * t371 + (qJD(3) * t272 - t201) * t273) * t359, t280, t249 * t352 + (-t276 * t295 + (t199 + t354) * t273) * t359, t249 * t342, -pkin(3) * t165 + t103 * t249 - t142 * t199 - t275 * t98 + (pkin(10) * t382 + t124 * t272) * qJD(4) + (-t273 * t75 + (-pkin(10) * t355 - t124 * t276) * t272) * t359, pkin(3) * t164 - t104 * t249 - t142 * t201 + t272 * t98 + (-pkin(10) * t295 + t124 * t275) * qJD(4) + (-t124 * t371 + (-pkin(10) * t354 + t76) * t273) * t359, t103 * t201 + t104 * t199 + ((-t165 + t353) * pkin(10) + t423) * t275 + ((qJD(4) * t199 - t164) * pkin(10) + t424) * t272, -pkin(3) * t98 - t103 * t75 - t104 * t76 - t124 * t142 + (t272 * t324 - t275 * t293 + (-t272 * t76 - t275 * t75) * qJD(4)) * pkin(10), t132 * t374 - t161 * t425, t186 * t325 + t161 * t185 + (t285 - t392) * t272 + t286 * t351, -t132 * t275 - t425 * t199 + (-t161 * t249 + t375) * t272, -t325 * t185 + (-t272 * t292 + t275 * t306) * t266, -t292 * t275 + t426 * t199 + (t325 * t341 - t306 - t380) * t272, t280, t197 * t165 - t96 * t325 - t67 * t185 + t398 * t199 + (-t12 + (pkin(10) * t325 + t67 * t266) * qJD(4)) * t275 + (-pkin(10) * t292 - t249 * t28 + t405) * t272, -t161 * t96 - t165 * t198 - t186 * t67 - t397 * t199 + (t13 + (pkin(10) * t161 + t269 * t67) * qJD(4)) * t275 + (pkin(10) * t132 + t249 * t29 + t404) * t272, t198 * t292 - t197 * t132 + t29 * t185 + t28 * t186 + (-t12 * t269 - t13 * t266) * t272 - t398 * t161 + (-t266 * t29 - t269 * t28) * t351 - t397 * t325, t12 * t197 + t13 * t198 - t67 * t96 + t397 * t29 + t398 * t28 + (t27 * t272 + t351 * t67) * pkin(10), -t213 * t35 + t367 * t88, t212 * t35 - t213 * t36 + t367 * t86 - t368 * t88, t165 * t213 - t193 * t367 + t275 * t35 + t295 * t88, t212 * t36 - t368 * t86, -t165 * t212 + t193 * t368 + t275 * t36 + t295 * t86, -t193 * t295 - t391, t105 * t165 - t193 * t412 - t2 * t275 + t212 * t22 + t235 * t36 + t295 * t300 - t368 * t56 + t396 * t86, t1 * t275 - t106 * t165 + t193 * t411 + t213 * t22 - t235 * t35 + t295 * t8 - t367 * t56 - t396 * t88, -t1 * t212 + t105 * t35 - t106 * t36 - t2 * t213 - t300 * t367 + t368 * t8 + t411 * t86 - t412 * t88, t1 * t106 + t105 * t2 + t22 * t235 + t300 * t412 + t396 * t56 - t411 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t388, t201 ^ 2 - t419, -t164 - t387, -t388, -t384 - t165, t314, -t124 * t201 - t424, t124 * t199 - t423, 0, 0, t161 * t385 + t392, t132 * t269 + t199 * t286 + t266 * t292, -t161 * t201 + t269 * t419 + t380, t325 * t386 + t285, t201 * t325 - t266 * t419 + t375, -t388, -qJ(5) * t380 + pkin(4) * t292 - t404 - t28 * t201 - t76 * t325 + (t266 * t369 - t52) * t199, -qJ(5) * t375 - pkin(4) * t132 - t161 * t76 + t201 * t29 + t405 + (t269 * t369 + t53) * t199, t52 * t161 + t53 * t184 + (qJ(5) * t159 - qJD(5) * t184 - t28 * t199 + t53 * t249 + t269 * t284 + t13) * t269 + (qJ(5) * t132 + qJD(5) * t161 - t199 * t29 - t12) * t266, -pkin(4) * t27 - t28 * t52 - t29 * t53 - t67 * t76 + (-t266 * t28 + t269 * t29) * qJD(5) + (-t12 * t266 + t13 * t269) * qJ(5), -t232 * t35 + t365 * t88, -t232 * t36 - t296 * t35 + t365 * t86 - t366 * t88, t165 * t232 - t193 * t365 + t201 * t88, -t296 * t36 - t366 * t86, t165 * t296 + t193 * t366 + t201 * t86, -t193 * t201, t165 * t182 - t193 * t403 + t201 * t300 - t22 * t296 + t262 * t36 - t366 * t56 - t63 * t86, -t165 * t183 - t193 * t402 + t201 * t8 + t22 * t232 - t262 * t35 - t365 * t56 + t63 * t88, t1 * t296 + t182 * t35 - t183 * t36 - t2 * t232 - t300 * t365 + t366 * t8 - t402 * t86 - t403 * t88, t1 * t183 + t182 * t2 + t22 * t262 + t300 * t403 + t402 * t8 - t56 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161 * t199 - t292, -t199 * t325 + t132, -t161 ^ 2 - t325 ^ 2, t161 * t28 + t29 * t325 + t27, 0, 0, 0, 0, 0, 0, -t88 * t193 + t36, -t35 - t438, -t437 - t439, t300 * t88 + t8 * t86 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t414, t437 - t439, -t35 + t438, t414, -t326 - (-qJD(6) + t193) * t88, t165, t8 * t193 + t56 * t88 + t2, -t193 * t300 + t56 * t86 - t1, 0, 0;];
tauc_reg  = t3;
