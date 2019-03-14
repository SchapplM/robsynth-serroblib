% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x38]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:38
% EndTime: 2019-03-10 03:43:01
% DurationCPUTime: 8.89s
% Computational Cost: add. (13773->517), mult. (32329->701), div. (0->0), fcn. (25170->10), ass. (0->291)
t305 = cos(qJ(2));
t431 = pkin(7) + pkin(8);
t276 = t431 * t305;
t266 = qJD(1) * t276;
t300 = sin(qJ(3));
t243 = t300 * t266;
t301 = sin(qJ(2));
t274 = t431 * t301;
t264 = qJD(1) * t274;
t429 = cos(qJ(3));
t195 = -t264 * t429 - t243;
t360 = qJD(3) * t429;
t470 = -pkin(2) * t360 + t195;
t432 = qJD(4) + qJD(5);
t361 = qJD(1) * t429;
t379 = qJD(1) * t301;
t462 = -t300 * t379 + t305 * t361;
t469 = t432 - t462;
t294 = qJD(2) + qJD(3);
t182 = t462 * t294;
t394 = t300 * t305;
t239 = -qJD(1) * t394 - t301 * t361;
t184 = -pkin(3) * t239 - pkin(9) * t462;
t167 = pkin(2) * t379 + t184;
t299 = sin(qJ(4));
t304 = cos(qJ(4));
t468 = -t304 * t167 + t299 * t470;
t207 = -t239 * t299 - t304 * t294;
t298 = sin(qJ(5));
t303 = cos(qJ(5));
t328 = t239 * t304 - t294 * t299;
t139 = t303 * t207 - t298 * t328;
t297 = sin(qJ(6));
t302 = cos(qJ(6));
t329 = t207 * t298 + t303 * t328;
t331 = t139 * t297 + t302 * t329;
t71 = t302 * t139 - t297 * t329;
t459 = t331 * t71;
t397 = t298 * t304;
t258 = t299 * t303 + t397;
t455 = t469 * t258;
t256 = t298 * t299 - t303 * t304;
t387 = t469 * t256;
t467 = t299 * t167 + t304 * t470;
t453 = t331 ^ 2 - t71 ^ 2;
t232 = qJD(4) - t462;
t290 = -pkin(2) * t305 - pkin(1);
t272 = t290 * qJD(1);
t164 = -pkin(3) * t462 + pkin(9) * t239 + t272;
t246 = t429 * t266;
t419 = qJD(2) * pkin(2);
t247 = -t264 + t419;
t189 = t300 * t247 + t246;
t171 = pkin(9) * t294 + t189;
t110 = t304 * t164 - t171 * t299;
t84 = pkin(10) * t328 + t110;
t67 = pkin(4) * t232 + t84;
t111 = t164 * t299 + t171 * t304;
t85 = -pkin(10) * t207 + t111;
t81 = t303 * t85;
t41 = t298 * t67 + t81;
t461 = pkin(11) * t139;
t29 = t41 - t461;
t373 = qJD(6) * t297;
t26 = t29 * t373;
t188 = t247 * t429 - t243;
t170 = -t294 * pkin(3) - t188;
t132 = t207 * pkin(4) + t170;
t68 = t139 * pkin(5) + t132;
t451 = t68 * t71 + t26;
t376 = qJD(4) * t304;
t377 = qJD(4) * t299;
t130 = t304 * t182 + t239 * t377 + t294 * t376;
t131 = -qJD(4) * t328 + t182 * t299;
t309 = qJD(5) * t329 - t130 * t298 - t303 * t131;
t372 = qJD(6) * t302;
t374 = qJD(5) * t303;
t375 = qJD(5) * t298;
t54 = t303 * t130 - t298 * t131 - t207 * t374 + t328 * t375;
t11 = -t139 * t372 + t297 * t309 + t302 * t54 + t329 * t373;
t222 = qJD(5) + t232;
t218 = qJD(6) + t222;
t450 = t218 * t71 + t11;
t259 = t301 * t429 + t394;
t200 = t294 * t259;
t183 = t200 * qJD(1);
t371 = qJD(1) * qJD(2);
t359 = t301 * t371;
t118 = pkin(2) * t359 + pkin(3) * t183 - pkin(9) * t182;
t114 = t304 * t118;
t367 = qJD(2) * t431;
t342 = qJD(1) * t367;
t248 = t301 * t342;
t249 = t305 * t342;
t378 = qJD(3) * t300;
t122 = t247 * t360 - t248 * t429 - t300 * t249 - t266 * t378;
t310 = -qJD(4) * t111 - t122 * t299 + t114;
t22 = pkin(4) * t183 - pkin(10) * t130 + t310;
t319 = t299 * t118 + t304 * t122 + t164 * t376 - t171 * t377;
t27 = -pkin(10) * t131 + t319;
t356 = t303 * t22 - t298 * t27;
t312 = -qJD(5) * t41 + t356;
t3 = pkin(5) * t183 - pkin(11) * t54 + t312;
t351 = -t298 * t22 - t303 * t27 - t67 * t374 + t85 * t375;
t4 = pkin(11) * t309 - t351;
t368 = -t297 * t4 + t302 * t3;
t466 = t68 * t331 + t368;
t311 = qJD(6) * t331 - t297 * t54 + t302 * t309;
t445 = -t218 * t331 + t311;
t286 = pkin(2) * t300 + pkin(9);
t423 = -pkin(10) - t286;
t352 = qJD(4) * t423;
t404 = t462 * t299;
t370 = pkin(10) * t404;
t458 = t299 * t352 + t370 - t467;
t403 = t462 * t304;
t341 = -t239 * pkin(4) - pkin(10) * t403;
t465 = t304 * t352 - t341 + t468;
t347 = t304 * t184 - t188 * t299;
t430 = -pkin(9) - pkin(10);
t366 = qJD(4) * t430;
t464 = t304 * t366 - t341 - t347;
t386 = t299 * t184 + t304 * t188;
t463 = -t299 * t366 - t370 + t386;
t460 = pkin(11) * t329;
t192 = -t256 * t297 + t258 * t302;
t417 = -qJD(6) * t192 + t297 * t387 - t302 * t455;
t191 = t302 * t256 + t258 * t297;
t416 = qJD(6) * t191 + t297 * t455 + t302 * t387;
t457 = t329 * t139;
t194 = -t300 * t264 + t246;
t340 = pkin(2) * t378 - t194;
t326 = -t300 * t301 + t305 * t429;
t199 = t294 * t326;
t395 = t299 * t199;
t454 = t259 * t376 + t395;
t452 = -t139 ^ 2 + t329 ^ 2;
t449 = t139 * t222 + t54;
t79 = t298 * t85;
t40 = t303 * t67 - t79;
t28 = t40 + t460;
t24 = pkin(5) * t222 + t28;
t418 = t29 * t302;
t10 = t24 * t297 + t418;
t448 = -qJD(6) * t10 + t466;
t447 = t132 * t139 + t351;
t446 = t132 * t329 + t312;
t444 = -t222 * t329 + t309;
t443 = -0.2e1 * t371;
t442 = t465 * t303;
t441 = t464 * t303;
t177 = t258 * t259;
t291 = pkin(4) * t377;
t440 = pkin(5) * t455 + t291;
t439 = t455 * pkin(11);
t251 = t423 * t299;
t293 = t304 * pkin(10);
t252 = t286 * t304 + t293;
t383 = t298 * t251 + t303 * t252;
t273 = t430 * t299;
t275 = pkin(9) * t304 + t293;
t382 = t298 * t273 + t303 * t275;
t219 = pkin(4) * t404;
t438 = -t219 + t340;
t437 = -t273 * t374 + t275 * t375 - t298 * t464 + t303 * t463;
t436 = -t251 * t374 + t252 * t375 - t465 * t298 - t303 * t458;
t435 = -t239 * pkin(5) - pkin(11) * t387;
t434 = -t429 * t274 - t300 * t276;
t433 = qJD(1) * t259;
t427 = pkin(11) * t258;
t425 = t256 * pkin(5);
t424 = t304 * pkin(4);
t422 = t303 * t84 - t79;
t187 = -pkin(3) * t326 - pkin(9) * t259 + t290;
t180 = t304 * t187;
t211 = -t300 * t274 + t276 * t429;
t400 = t259 * t304;
t100 = -pkin(4) * t326 - pkin(10) * t400 - t211 * t299 + t180;
t203 = t304 * t211;
t384 = t299 * t187 + t203;
t401 = t259 * t299;
t115 = -pkin(10) * t401 + t384;
t415 = t298 * t100 + t303 * t115;
t147 = t219 + t189;
t414 = -t147 + t440;
t413 = t130 * t299;
t411 = t183 * t298;
t410 = t199 * t304;
t409 = t207 * t232;
t408 = t328 * t232;
t407 = t218 * t239;
t406 = t222 * t239;
t405 = t232 * t239;
t402 = t239 * t462;
t287 = pkin(4) * t303 + pkin(5);
t399 = t287 * t183;
t398 = t298 * t302;
t396 = t299 * t183;
t393 = t304 * t183;
t307 = qJD(1) ^ 2;
t392 = t305 * t307;
t306 = qJD(2) ^ 2;
t391 = t306 * t301;
t390 = t306 * t305;
t389 = t438 + t440;
t381 = t291 + t438;
t380 = t301 ^ 2 - t305 ^ 2;
t369 = t301 * t419;
t289 = -pkin(3) - t424;
t363 = t259 * t377;
t158 = t170 * t376;
t357 = qJD(6) * t24 + t4;
t129 = pkin(3) * t200 - pkin(9) * t199 + t369;
t125 = t304 * t129;
t265 = t301 * t367;
t267 = t305 * t367;
t143 = qJD(3) * t434 - t429 * t265 - t300 * t267;
t35 = -pkin(10) * t410 + pkin(4) * t200 - t143 * t299 + t125 + (-t203 + (pkin(10) * t259 - t187) * t299) * qJD(4);
t318 = t299 * t129 + t304 * t143 + t187 * t376 - t211 * t377;
t44 = -pkin(10) * t454 + t318;
t354 = -t298 * t44 + t303 * t35;
t353 = -t298 * t84 - t81;
t350 = t303 * t100 - t115 * t298;
t349 = pkin(1) * t443;
t346 = t303 * t251 - t252 * t298;
t345 = t303 * t273 - t275 * t298;
t344 = t232 * t304;
t123 = t247 * t378 - t300 * t248 + t429 * t249 + t266 * t360;
t288 = -pkin(2) * t429 - pkin(3);
t339 = -t147 + t291;
t149 = t346 - t427;
t338 = -qJD(6) * t149 + t436 + t439;
t250 = t256 * pkin(11);
t150 = -t250 + t383;
t337 = t383 * qJD(5) + qJD(6) * t150 + t298 * t458 + t435 - t442;
t168 = t345 - t427;
t336 = -qJD(6) * t168 + t437 + t439;
t169 = -t250 + t382;
t335 = t382 * qJD(5) + qJD(6) * t169 - t298 * t463 + t435 - t441;
t334 = -t111 * t239 + t123 * t299 + t158;
t330 = -t170 * t462 - t183 * t286;
t178 = t256 * t259;
t116 = t302 * t177 - t178 * t297;
t117 = -t177 * t297 - t178 * t302;
t165 = pkin(4) * t401 - t434;
t327 = t110 * t239 - t123 * t304 + t170 * t377;
t324 = -t363 + t410;
t323 = t100 * t374 - t115 * t375 + t298 * t35 + t303 * t44;
t62 = pkin(4) * t131 + t123;
t271 = t288 - t424;
t23 = -pkin(5) * t309 + t62;
t9 = t302 * t24 - t29 * t297;
t322 = t23 * t191 + t9 * t239 - t417 * t68;
t321 = -t10 * t239 + t23 * t192 - t416 * t68;
t320 = t239 * t272 - t123;
t144 = -t300 * t265 + t267 * t429 - t274 * t378 + t276 * t360;
t315 = -t132 * t387 - t41 * t239 + t62 * t258;
t314 = t132 * t455 + t40 * t239 + t62 * t256;
t90 = pkin(4) * t454 + t144;
t308 = -t272 * t462 - t122;
t215 = t289 + t425;
t214 = t271 + t425;
t152 = t239 ^ 2 - t462 ^ 2;
t151 = t183 * t326;
t146 = (-t239 - t433) * t294;
t119 = t177 * pkin(5) + t165;
t109 = -pkin(4) * t328 - pkin(5) * t329;
t66 = t232 * t344 - t239 * t328 + t396;
t65 = -t232 ^ 2 * t299 - t207 * t239 + t393;
t61 = -t328 * t344 + t413;
t60 = t199 * t397 - t298 * t363 - t375 * t401 + (t400 * t432 + t395) * t303;
t59 = -t177 * t432 - t256 * t199;
t50 = -pkin(11) * t177 + t415;
t47 = -t139 * t239 - t183 * t256 - t222 * t455;
t46 = t183 * t258 - t222 * t387 - t239 * t329;
t45 = -pkin(5) * t326 + pkin(11) * t178 + t350;
t39 = t60 * pkin(5) + t90;
t32 = (t130 - t409) * t304 + (-t131 + t408) * t299;
t31 = t422 + t460;
t30 = t353 + t461;
t17 = qJD(6) * t117 + t297 * t59 + t302 * t60;
t16 = -qJD(6) * t116 - t297 * t60 + t302 * t59;
t15 = t258 * t54 + t329 * t387;
t14 = -t183 * t191 + t218 * t417 - t239 * t71;
t13 = t183 * t192 - t218 * t416 - t239 * t331;
t8 = t139 * t387 - t256 * t54 + t258 * t309 + t329 * t455;
t7 = -pkin(11) * t60 + t323;
t6 = pkin(5) * t200 - pkin(11) * t59 - qJD(5) * t415 + t354;
t5 = t11 * t192 + t331 * t416;
t1 = -t11 * t191 + t192 * t311 - t331 * t417 + t416 * t71;
t2 = [0, 0, 0, 0.2e1 * t305 * t359, t380 * t443, t390, -t391, 0, -pkin(7) * t390 + t301 * t349, pkin(7) * t391 + t305 * t349, t182 * t259 - t199 * t239, t182 * t326 - t183 * t259 + t199 * t462 + t200 * t239, t199 * t294, -t200 * t294, 0, -t144 * t294 + t183 * t290 + t200 * t272 + (-qJD(1) * t326 - t462) * t369, -t143 * t294 + t182 * t290 + t199 * t272 + (-t239 + t433) * t369, t130 * t400 - t324 * t328 (-t207 * t304 + t299 * t328) * t199 + (-t413 - t131 * t304 + (t207 * t299 + t304 * t328) * qJD(4)) * t259, -t130 * t326 - t200 * t328 + t232 * t324 + t259 * t393, t131 * t326 - t200 * t207 - t232 * t454 - t259 * t396, t200 * t232 - t151 (-t211 * t376 + t125) * t232 + t180 * t183 - (-t171 * t376 + t114) * t326 + t110 * t200 + t144 * t207 - t434 * t131 + t259 * t158 + ((-qJD(4) * t187 - t143) * t232 - t211 * t183 - (-qJD(4) * t164 - t122) * t326 + t123 * t259 + t170 * t199) * t299, -t111 * t200 + t123 * t400 - t130 * t434 - t144 * t328 + t170 * t324 - t183 * t384 - t232 * t318 + t319 * t326, -t178 * t54 - t329 * t59, -t139 * t59 - t177 * t54 - t178 * t309 + t329 * t60, -t178 * t183 - t200 * t329 + t222 * t59 - t326 * t54, -t139 * t200 - t177 * t183 - t222 * t60 - t309 * t326, t200 * t222 - t151, t354 * t222 + t350 * t183 - t356 * t326 + t40 * t200 + t90 * t139 - t165 * t309 + t62 * t177 + t132 * t60 + (-t222 * t415 + t326 * t41) * qJD(5), t132 * t59 + t165 * t54 - t62 * t178 - t183 * t415 - t41 * t200 - t222 * t323 - t326 * t351 - t329 * t90, t11 * t117 - t16 * t331, -t11 * t116 + t117 * t311 - t16 * t71 + t17 * t331, -t11 * t326 + t117 * t183 + t16 * t218 - t200 * t331, -t116 * t183 - t17 * t218 - t200 * t71 - t311 * t326, t200 * t218 - t151 (-t297 * t7 + t302 * t6) * t218 + (-t297 * t50 + t302 * t45) * t183 - t368 * t326 + t9 * t200 + t39 * t71 - t119 * t311 + t23 * t116 + t68 * t17 + ((-t297 * t45 - t302 * t50) * t218 + t10 * t326) * qJD(6), -t10 * t200 + t119 * t11 + t23 * t117 + t68 * t16 - t26 * t326 - t39 * t331 + (-(-qJD(6) * t50 + t6) * t218 - t45 * t183 + t3 * t326) * t297 + (-(qJD(6) * t45 + t7) * t218 - t50 * t183 + t357 * t326) * t302; 0, 0, 0, -t301 * t392, t380 * t307, 0, 0, 0, t307 * pkin(1) * t301, pkin(1) * t392, t402, t152, 0, t146, 0, t194 * t294 + (-t294 * t378 + t379 * t462) * pkin(2) + t320, t195 * t294 + (t239 * t379 - t294 * t360) * pkin(2) + t308, t61, t32, t66, t65, t405, t288 * t131 + t330 * t299 + t340 * t207 + (-t286 * t376 + t468) * t232 + t327, t288 * t130 + t330 * t304 - t340 * t328 + (t286 * t377 + t467) * t232 + t334, t15, t8, t46, t47, t406, t346 * t183 - t271 * t309 + (-t252 * t374 + (-qJD(5) * t251 - t458) * t298 + t442) * t222 + t381 * t139 + t314, -t383 * t183 + t222 * t436 + t271 * t54 - t329 * t381 + t315, t5, t1, t13, t14, t407 (t149 * t302 - t150 * t297) * t183 - t214 * t311 + t389 * t71 + (t297 * t338 - t302 * t337) * t218 + t322 -(t149 * t297 + t150 * t302) * t183 + t214 * t11 - t389 * t331 + (t297 * t337 + t302 * t338) * t218 + t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t402, t152, 0, t146, 0, t189 * t294 + t320, t188 * t294 + t308, t61, t32, t66, t65, t405, -pkin(3) * t131 - t347 * t232 - t189 * t207 - t170 * t404 + (-t232 * t376 - t396) * pkin(9) + t327, -pkin(3) * t130 + t386 * t232 + t189 * t328 - t170 * t403 + (t232 * t377 - t393) * pkin(9) + t334, t15, t8, t46, t47, t406, t345 * t183 - t289 * t309 + (-t275 * t374 + (-qJD(5) * t273 + t463) * t298 + t441) * t222 + t339 * t139 + t314, -t382 * t183 + t222 * t437 + t289 * t54 - t329 * t339 + t315, t5, t1, t13, t14, t407 (t168 * t302 - t169 * t297) * t183 - t215 * t311 + t414 * t71 + (t297 * t336 - t302 * t335) * t218 + t322 -(t168 * t297 + t169 * t302) * t183 + t215 * t11 - t414 * t331 + (t297 * t335 + t302 * t336) * t218 + t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t328 * t207, -t207 ^ 2 + t328 ^ 2, t130 + t409, -t131 - t408, t183, t111 * t232 + t170 * t328 + t310, t110 * t232 + t170 * t207 - t319, -t457, t452, t449, t444, t183, -t353 * t222 + (t139 * t328 + t183 * t303 - t222 * t375) * pkin(4) + t446, t422 * t222 + (-t222 * t374 - t328 * t329 - t411) * pkin(4) + t447, -t459, t453, t450, t445, t183, t302 * t399 - (-t297 * t31 + t30 * t302) * t218 - t109 * t71 + (-t297 * t411 + (-t297 * t303 - t398) * t218 * qJD(5)) * pkin(4) + ((-pkin(4) * t398 - t287 * t297) * t218 - t10) * qJD(6) + t466, t109 * t331 + (-t399 - t3 + (t30 - (-qJD(5) - qJD(6)) * t298 * pkin(4)) * t218) * t297 + (-pkin(4) * t411 + (-pkin(4) * t374 - qJD(6) * t287 + t31) * t218 - t357) * t302 + t451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t457, t452, t449, t444, t183, t222 * t41 + t446, t222 * t40 + t447, -t459, t453, t450, t445, t183 -(-t28 * t297 - t418) * t218 + (t183 * t302 - t218 * t373 + t329 * t71) * pkin(5) + t448 (-t218 * t29 - t3) * t297 + (t218 * t28 - t357) * t302 + (-t183 * t297 - t218 * t372 - t329 * t331) * pkin(5) + t451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t459, t453, t450, t445, t183, t10 * t218 + t448, t218 * t9 - t297 * t3 - t302 * t357 + t451;];
tauc_reg  = t2;