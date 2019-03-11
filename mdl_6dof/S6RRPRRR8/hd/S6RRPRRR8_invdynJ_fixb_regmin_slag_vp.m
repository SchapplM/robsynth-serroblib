% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:34
% EndTime: 2019-03-09 14:05:57
% DurationCPUTime: 10.34s
% Computational Cost: add. (13492->594), mult. (31233->803), div. (0->0), fcn. (24783->18), ass. (0->282)
t300 = cos(qJ(2));
t381 = qJD(1) * t300;
t269 = -qJD(4) + t381;
t262 = -qJD(5) + t269;
t256 = -qJD(6) + t262;
t290 = sin(pkin(11));
t295 = sin(qJ(2));
t382 = qJD(1) * t295;
t361 = t290 * t382;
t291 = cos(pkin(11));
t369 = t291 * qJD(2);
t234 = t361 - t369;
t360 = t291 * t382;
t380 = qJD(2) * t290;
t236 = t360 + t380;
t294 = sin(qJ(4));
t299 = cos(qJ(4));
t165 = t234 * t294 - t236 * t299;
t166 = t299 * t234 + t236 * t294;
t293 = sin(qJ(5));
t298 = cos(qJ(5));
t109 = t165 * t293 - t298 * t166;
t292 = sin(qJ(6));
t297 = cos(qJ(6));
t420 = t165 * t298 + t166 * t293;
t442 = t109 * t297 + t292 * t420;
t449 = t256 * t442;
t337 = pkin(2) * t295 - qJ(3) * t300;
t245 = t337 * qJD(1);
t225 = t290 * t245;
t396 = t291 * t295;
t397 = t290 * t300;
t325 = -pkin(7) * t396 - pkin(8) * t397;
t172 = qJD(1) * t325 + t225;
t448 = qJD(3) * t291 - t172;
t399 = t290 * t294;
t242 = -t299 * t291 + t399;
t323 = t242 * t300;
t447 = -qJD(1) * t323 + t242 * qJD(4);
t243 = t290 * t299 + t291 * t294;
t324 = t243 * t300;
t385 = -qJD(1) * t324 + t243 * qJD(4);
t335 = -t109 * t292 + t297 * t420;
t435 = t335 * t442;
t446 = t109 * t262;
t445 = t256 * t335;
t185 = pkin(7) * t361 + t291 * t245;
t395 = t291 * t300;
t331 = pkin(3) * t295 - pkin(8) * t395;
t153 = qJD(1) * t331 + t185;
t409 = pkin(8) + qJ(3);
t254 = t409 * t290;
t255 = t409 * t291;
t332 = qJD(3) * t290 + qJD(4) * t255;
t374 = qJD(4) * t299;
t444 = -t254 * t374 + t448 * t299 + (-t153 - t332) * t294;
t430 = t335 ^ 2 - t442 ^ 2;
t280 = t291 * qJDD(2);
t367 = qJD(1) * qJD(2);
t357 = t300 * t367;
t365 = t295 * qJDD(1);
t322 = t357 + t365;
t193 = t290 * t322 - t280;
t366 = t290 * qJDD(2);
t194 = t291 * t322 + t366;
t376 = qJD(4) * t294;
t94 = -t294 * t193 + t299 * t194 - t234 * t374 - t236 * t376;
t95 = -qJD(4) * t165 + t299 * t193 + t294 * t194;
t306 = qJD(5) * t420 - t293 * t94 - t298 * t95;
t372 = qJD(5) * t298;
t373 = qJD(5) * t293;
t33 = t165 * t373 - t166 * t372 - t293 * t95 + t298 * t94;
t370 = qJD(6) * t297;
t371 = qJD(6) * t292;
t8 = t109 * t370 + t292 * t306 + t297 * t33 + t371 * t420;
t428 = t8 + t449;
t330 = pkin(2) * t300 + qJ(3) * t295 + pkin(1);
t227 = t330 * qJD(1);
t278 = pkin(7) * t381;
t257 = qJD(2) * qJ(3) + t278;
t174 = -t291 * t227 - t257 * t290;
t364 = pkin(3) * t381;
t127 = -pkin(8) * t236 + t174 - t364;
t175 = -t290 * t227 + t291 * t257;
t131 = -pkin(8) * t234 + t175;
t72 = t299 * t127 - t131 * t294;
t59 = pkin(9) * t165 + t72;
t54 = -pkin(4) * t269 + t59;
t73 = t127 * t294 + t131 * t299;
t60 = -pkin(9) * t166 + t73;
t58 = t298 * t60;
t27 = t293 * t54 + t58;
t437 = pkin(10) * t109;
t19 = t27 + t437;
t17 = t19 * t371;
t287 = pkin(11) + qJ(4);
t285 = qJ(5) + t287;
t274 = qJ(6) + t285;
t264 = sin(t274);
t265 = cos(t274);
t301 = cos(qJ(1));
t296 = sin(qJ(1));
t393 = t296 * t300;
t188 = t264 * t301 - t265 * t393;
t392 = t300 * t301;
t190 = t264 * t296 + t265 * t392;
t411 = g(3) * t295;
t249 = -qJD(2) * pkin(2) + pkin(7) * t382 + qJD(3);
t184 = pkin(3) * t234 + t249;
t122 = pkin(4) * t166 + t184;
t62 = -pkin(5) * t109 + t122;
t427 = g(1) * t190 - g(2) * t188 + t265 * t411 - t442 * t62 + t17;
t187 = t264 * t393 + t265 * t301;
t189 = -t264 * t392 + t265 * t296;
t284 = t300 * qJDD(1);
t321 = t295 * t367 - t284;
t241 = qJDD(4) + t321;
t228 = qJDD(5) + t241;
t219 = qJD(2) * t337 - t295 * qJD(3);
t164 = qJD(1) * t219 - qJDD(1) * t330;
t207 = -pkin(7) * t321 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t118 = t291 * t164 - t290 * t207;
t80 = pkin(3) * t321 - t194 * pkin(8) + t118;
t119 = t290 * t164 + t291 * t207;
t97 = -pkin(8) * t193 + t119;
t348 = -t294 * t97 + t299 * t80;
t308 = -qJD(4) * t73 + t348;
t12 = t241 * pkin(4) - t94 * pkin(9) + t308;
t326 = t127 * t374 - t131 * t376 + t294 * t80 + t299 * t97;
t14 = -pkin(9) * t95 + t326;
t354 = t298 * t12 - t293 * t14;
t309 = -qJD(5) * t27 + t354;
t2 = t228 * pkin(5) - t33 * pkin(10) + t309;
t347 = -t293 * t12 - t298 * t14 - t54 * t372 + t60 * t373;
t3 = pkin(10) * t306 - t347;
t362 = t297 * t2 - t292 * t3;
t443 = -g(1) * t189 + g(2) * t187 + t264 * t411 + t62 * t335 + t362;
t307 = qJD(6) * t335 - t292 * t33 + t297 * t306;
t422 = t307 + t445;
t434 = -pkin(9) * t385 + t444;
t145 = t299 * t153;
t384 = -t294 * t254 + t299 * t255;
t441 = -pkin(4) * t382 + pkin(9) * t447 - t243 * qJD(3) - qJD(4) * t384 + t294 * t172 - t145;
t440 = t262 * t420;
t276 = pkin(7) * t365;
t217 = -qJDD(2) * pkin(2) + pkin(7) * t357 + qJDD(3) + t276;
t339 = g(1) * t301 + g(2) * t296;
t410 = g(3) * t300;
t316 = t295 * t339 - t410;
t418 = t217 - t316;
t436 = pkin(10) * t420;
t433 = t165 * t269;
t432 = t166 * t269;
t431 = t420 * t109;
t170 = t298 * t242 + t243 * t293;
t391 = -qJD(5) * t170 - t293 * t385 - t298 * t447;
t171 = -t242 * t293 + t243 * t298;
t390 = qJD(5) * t171 - t293 * t447 + t298 * t385;
t429 = -t109 ^ 2 + t420 ^ 2;
t426 = t33 + t446;
t56 = t293 * t60;
t26 = t298 * t54 - t56;
t18 = t26 + t436;
t16 = -pkin(5) * t262 + t18;
t402 = t297 * t19;
t5 = t292 * t16 + t402;
t425 = -qJD(6) * t5 + t443;
t271 = sin(t285);
t272 = cos(t285);
t196 = t271 * t301 - t272 * t393;
t198 = t271 * t296 + t272 * t392;
t424 = g(1) * t198 - g(2) * t196 - t109 * t122 + t272 * t411 + t347;
t195 = t271 * t393 + t272 * t301;
t197 = -t271 * t392 + t272 * t296;
t423 = -g(1) * t197 + g(2) * t195 + t122 * t420 + t271 * t411 + t309;
t421 = t306 + t440;
t419 = t441 * t298;
t342 = -t299 * t254 - t255 * t294;
t146 = -pkin(9) * t243 + t342;
t147 = -pkin(9) * t242 + t384;
t389 = t293 * t146 + t298 * t147;
t229 = t290 * t364 + t278;
t356 = pkin(4) * t385 - t229;
t233 = t291 * t330;
t173 = -pkin(8) * t396 - t233 + (-pkin(7) * t290 - pkin(3)) * t300;
t200 = pkin(7) * t395 - t290 * t330;
t398 = t290 * t295;
t180 = -pkin(8) * t398 + t200;
t387 = t294 * t173 + t299 * t180;
t417 = t146 * t372 - t147 * t373 + t441 * t293 + t298 * t434;
t338 = g(1) * t296 - g(2) * t301;
t415 = pkin(7) * t234;
t414 = pkin(7) * t236;
t113 = t297 * t170 + t171 * t292;
t408 = -qJD(6) * t113 - t292 * t390 + t297 * t391;
t114 = -t170 * t292 + t171 * t297;
t407 = qJD(6) * t114 + t292 * t391 + t297 * t390;
t406 = t298 * t59 - t56;
t404 = pkin(5) * t390 + t356;
t216 = t242 * t295;
t344 = t299 * t173 - t180 * t294;
t87 = -pkin(4) * t300 + pkin(9) * t216 + t344;
t215 = t243 * t295;
t93 = -pkin(9) * t215 + t387;
t403 = t293 * t87 + t298 * t93;
t222 = qJDD(6) + t228;
t401 = t222 * t293;
t275 = pkin(4) * t298 + pkin(5);
t400 = t275 * t222;
t394 = t293 * t297;
t379 = qJD(2) * t295;
t363 = pkin(7) * t379;
t178 = t291 * t219 + t290 * t363;
t378 = qJD(2) * t300;
t230 = (pkin(3) * t290 + pkin(7)) * t378;
t246 = pkin(3) * t398 + t295 * pkin(7);
t288 = t295 ^ 2;
t383 = -t300 ^ 2 + t288;
t375 = qJD(4) * t295;
t368 = qJD(3) - t249;
t273 = -pkin(3) * t291 - pkin(2);
t359 = qJ(3) * t284;
t355 = qJD(6) * t16 + t3;
t148 = -qJD(2) * t323 - t243 * t375;
t142 = qJD(2) * t331 + t178;
t205 = t290 * t219;
t155 = qJD(2) * t325 + t205;
t346 = t299 * t142 - t294 * t155;
t41 = pkin(4) * t379 - t148 * pkin(9) - qJD(4) * t387 + t346;
t149 = qJD(2) * t324 + t374 * t396 - t375 * t399;
t320 = t294 * t142 + t299 * t155 + t173 * t374 - t180 * t376;
t43 = -pkin(9) * t149 + t320;
t352 = -t293 * t43 + t298 * t41;
t351 = -t293 * t59 - t58;
t350 = -t293 * t93 + t298 * t87;
t345 = t298 * t146 - t147 * t293;
t123 = pkin(4) * t149 + t230;
t177 = pkin(4) * t215 + t246;
t66 = -pkin(10) * t171 + t345;
t341 = -pkin(10) * t390 + qJD(6) * t66 + t417;
t67 = -pkin(10) * t170 + t389;
t340 = pkin(5) * t382 + pkin(10) * t391 + t389 * qJD(5) + qJD(6) * t67 + t293 * t434 - t419;
t140 = t298 * t215 - t216 * t293;
t141 = -t215 * t293 - t216 * t298;
t84 = t297 * t140 + t141 * t292;
t85 = -t140 * t292 + t141 * t297;
t204 = pkin(4) * t242 + t273;
t329 = -0.2e1 * pkin(1) * t367 - pkin(7) * qJDD(2);
t328 = t293 * t41 + t298 * t43 + t87 * t372 - t373 * t93;
t303 = qJD(1) ^ 2;
t318 = pkin(1) * t303 + t339;
t143 = pkin(3) * t193 + t217;
t314 = -t300 * t339 - t411;
t302 = qJD(2) ^ 2;
t312 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t302 + t338;
t65 = pkin(4) * t95 + t143;
t283 = cos(t287);
t282 = sin(t287);
t211 = t282 * t296 + t283 * t392;
t210 = -t282 * t392 + t283 * t296;
t209 = t282 * t301 - t283 * t393;
t208 = t282 * t393 + t283 * t301;
t199 = -pkin(7) * t397 - t233;
t186 = -pkin(7) * t360 + t225;
t179 = -t291 * t363 + t205;
t126 = pkin(5) * t170 + t204;
t103 = pkin(5) * t140 + t177;
t68 = -pkin(4) * t165 - pkin(5) * t420;
t64 = qJD(5) * t141 + t293 * t148 + t298 * t149;
t63 = -qJD(5) * t140 + t298 * t148 - t293 * t149;
t44 = pkin(5) * t64 + t123;
t36 = -pkin(10) * t140 + t403;
t35 = -pkin(5) * t300 - pkin(10) * t141 + t350;
t23 = t406 + t436;
t22 = t351 - t437;
t21 = qJD(6) * t85 + t292 * t63 + t297 * t64;
t20 = -qJD(6) * t84 - t292 * t64 + t297 * t63;
t15 = -pkin(5) * t306 + t65;
t7 = -pkin(10) * t64 + t328;
t6 = pkin(5) * t379 - t63 * pkin(10) - qJD(5) * t403 + t352;
t4 = t297 * t16 - t19 * t292;
t1 = [qJDD(1), t338, t339, qJDD(1) * t288 + 0.2e1 * t295 * t357, 0.2e1 * t284 * t295 - 0.2e1 * t367 * t383, qJDD(2) * t295 + t300 * t302, qJDD(2) * t300 - t295 * t302, 0, t295 * t329 + t300 * t312, -t295 * t312 + t300 * t329, -t339 * t290 + (pkin(7) * t193 + t217 * t290 + (qJD(1) * t199 + t174) * qJD(2)) * t295 + (-t178 * qJD(1) - t199 * qJDD(1) - t118 + t338 * t291 + (t249 * t290 + t415) * qJD(2)) * t300, -t339 * t291 + (pkin(7) * t194 + t217 * t291 + (-qJD(1) * t200 - t175) * qJD(2)) * t295 + (t179 * qJD(1) + t200 * qJDD(1) + t119 - t338 * t290 + (t249 * t291 + t414) * qJD(2)) * t300, -t178 * t236 - t179 * t234 - t200 * t193 - t199 * t194 + (-t174 * t291 - t175 * t290) * t378 + (-t118 * t291 - t119 * t290 + t338) * t295, t118 * t199 + t119 * t200 + t174 * t178 + t175 * t179 + (t217 * t295 + t249 * t378 - t339) * pkin(7) + t338 * t330, -t148 * t165 - t216 * t94, -t148 * t166 + t149 * t165 - t215 * t94 + t216 * t95, -t148 * t269 - t165 * t379 - t216 * t241 - t300 * t94, t149 * t269 - t166 * t379 - t215 * t241 + t300 * t95, -t241 * t300 - t269 * t379, -t346 * t269 + t344 * t241 - t348 * t300 + t72 * t379 + t230 * t166 + t246 * t95 + t143 * t215 + t184 * t149 - g(1) * t209 - g(2) * t211 + (t269 * t387 + t300 * t73) * qJD(4), -g(1) * t208 - g(2) * t210 - t143 * t216 + t184 * t148 - t165 * t230 - t241 * t387 + t246 * t94 + t269 * t320 + t300 * t326 - t379 * t73, t141 * t33 - t420 * t63, t109 * t63 - t140 * t33 + t141 * t306 + t420 * t64, t141 * t228 - t262 * t63 - t300 * t33 - t379 * t420, t109 * t379 - t140 * t228 + t262 * t64 - t300 * t306, -t228 * t300 - t262 * t379, -t352 * t262 + t350 * t228 - t354 * t300 + t26 * t379 - t123 * t109 - t177 * t306 + t65 * t140 + t122 * t64 - g(1) * t196 - g(2) * t198 + (t262 * t403 + t27 * t300) * qJD(5), -g(1) * t195 - g(2) * t197 + t122 * t63 - t123 * t420 + t65 * t141 + t177 * t33 - t228 * t403 + t262 * t328 - t27 * t379 - t300 * t347, -t20 * t335 + t8 * t85, t20 * t442 + t21 * t335 + t307 * t85 - t8 * t84, -t20 * t256 + t222 * t85 - t300 * t8 - t335 * t379, t21 * t256 - t222 * t84 - t300 * t307 + t379 * t442, -t222 * t300 - t256 * t379 -(-t292 * t7 + t297 * t6) * t256 + (-t292 * t36 + t297 * t35) * t222 - t362 * t300 + t4 * t379 - t44 * t442 - t103 * t307 + t15 * t84 + t62 * t21 - g(1) * t188 - g(2) * t190 + (-(-t292 * t35 - t297 * t36) * t256 + t5 * t300) * qJD(6), -t5 * t379 - g(1) * t187 - g(2) * t189 + t103 * t8 + t15 * t85 - t17 * t300 + t62 * t20 - t44 * t335 + ((-qJD(6) * t36 + t6) * t256 - t35 * t222 + t2 * t300) * t292 + ((qJD(6) * t35 + t7) * t256 - t36 * t222 + t355 * t300) * t297; 0, 0, 0, -t295 * t303 * t300, t383 * t303, t365, t284, qJDD(2), t295 * t318 - t276 - t410, t411 + (-pkin(7) * qJDD(1) + t318) * t300, t290 * t359 - pkin(2) * t193 - t418 * t291 + ((-qJ(3) * t380 - t174) * t295 + (t290 * t368 + t185 - t415) * t300) * qJD(1), t291 * t359 - pkin(2) * t194 + t418 * t290 + ((-qJ(3) * t369 + t175) * t295 + (t291 * t368 - t186 - t414) * t300) * qJD(1), t185 * t236 + t186 * t234 + (-qJ(3) * t193 - qJD(3) * t234 + t174 * t381 + t119) * t291 + (qJ(3) * t194 + qJD(3) * t236 + t175 * t381 - t118) * t290 + t314, -t249 * t278 - t174 * t185 - t175 * t186 + (-t174 * t290 + t175 * t291) * qJD(3) - t418 * pkin(2) + (-t118 * t290 + t119 * t291 + t314) * qJ(3), t165 * t447 + t94 * t243, t165 * t385 + t166 * t447 - t94 * t242 - t243 * t95, t165 * t382 + t243 * t241 + t269 * t447, t166 * t382 - t242 * t241 + t269 * t385, t269 * t382, t342 * t241 + t273 * t95 + t143 * t242 - t72 * t382 - t229 * t166 + (t145 + t332 * t299 + (-qJD(4) * t254 + t448) * t294) * t269 + t385 * t184 + t316 * t283, t143 * t243 + t229 * t165 - t184 * t447 - t384 * t241 + t269 * t444 + t273 * t94 - t316 * t282 + t73 * t382, t171 * t33 - t391 * t420, t109 * t391 - t170 * t33 + t171 * t306 + t390 * t420, t171 * t228 - t262 * t391 + t382 * t420, -t109 * t382 - t170 * t228 + t262 * t390, t262 * t382, t345 * t228 - t204 * t306 + t65 * t170 - t26 * t382 + (t147 * t372 + (qJD(5) * t146 + t434) * t293 - t419) * t262 + t390 * t122 - t356 * t109 + t316 * t272, t391 * t122 + t65 * t171 + t204 * t33 - t389 * t228 + t262 * t417 + t27 * t382 - t271 * t316 - t356 * t420, t114 * t8 - t335 * t408, -t113 * t8 + t114 * t307 + t335 * t407 + t408 * t442, t114 * t222 - t256 * t408 + t335 * t382, -t113 * t222 + t256 * t407 - t382 * t442, t256 * t382 (-t292 * t67 + t297 * t66) * t222 - t126 * t307 + t15 * t113 - t4 * t382 + t407 * t62 - t404 * t442 + (t292 * t341 + t297 * t340) * t256 + t316 * t265 -(t292 * t66 + t297 * t67) * t222 + t126 * t8 + t15 * t114 + t5 * t382 + t408 * t62 - t404 * t335 + (-t292 * t340 + t297 * t341) * t256 - t316 * t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290 * t365 - t280 + (-t236 + t380) * t381, t291 * t365 + t366 + (t234 + t369) * t381, -t234 ^ 2 - t236 ^ 2, t174 * t236 + t175 * t234 + t418, 0, 0, 0, 0, 0, t95 + t433, t94 + t432, 0, 0, 0, 0, 0, -t306 + t440, t33 - t446, 0, 0, 0, 0, 0, -t307 + t445, t8 - t449; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165 * t166, t165 ^ 2 - t166 ^ 2, t94 - t432, -t95 + t433, t241, -g(1) * t210 + g(2) * t208 + t165 * t184 - t73 * t269 + t282 * t411 + t308, g(1) * t211 - g(2) * t209 + t166 * t184 - t269 * t72 + t283 * t411 - t326, t431, t429, t426, t421, t228, t351 * t262 + (-t109 * t165 + t228 * t298 + t262 * t373) * pkin(4) + t423, -t406 * t262 + (-t165 * t420 - t228 * t293 + t262 * t372) * pkin(4) + t424, t435, t430, t428, t422, t222, t297 * t400 + (t22 * t297 - t23 * t292) * t256 + t68 * t442 + (-t292 * t401 - (-t292 * t298 - t394) * t256 * qJD(5)) * pkin(4) + (-(-pkin(4) * t394 - t275 * t292) * t256 - t5) * qJD(6) + t443, t68 * t335 + (-t400 - t2 + (-t22 + (-qJD(5) - qJD(6)) * t293 * pkin(4)) * t256) * t292 + (-pkin(4) * t401 + (pkin(4) * t372 + qJD(6) * t275 - t23) * t256 - t355) * t297 + t427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t431, t429, t426, t421, t228, -t27 * t262 + t423, -t26 * t262 + t424, t435, t430, t428, t422, t222 (-t18 * t292 - t402) * t256 + (t222 * t297 + t256 * t371 - t420 * t442) * pkin(5) + t425 (t19 * t256 - t2) * t292 + (-t18 * t256 - t355) * t297 + (-t222 * t292 + t256 * t370 - t335 * t420) * pkin(5) + t427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t435, t430, t428, t422, t222, -t5 * t256 + t425, -t292 * t2 - t4 * t256 - t297 * t355 + t427;];
tau_reg  = t1;
