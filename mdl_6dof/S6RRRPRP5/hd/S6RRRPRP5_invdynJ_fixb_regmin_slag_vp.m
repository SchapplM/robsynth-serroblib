% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:38
% EndTime: 2019-03-09 16:51:58
% DurationCPUTime: 10.04s
% Computational Cost: add. (12915->615), mult. (28929->773), div. (0->0), fcn. (20974->14), ass. (0->313)
t264 = sin(qJ(3));
t265 = sin(qJ(2));
t378 = qJD(1) * t265;
t351 = t264 * t378;
t267 = cos(qJ(3));
t364 = t267 * qJD(2);
t204 = t351 - t364;
t375 = qJD(2) * t264;
t206 = t267 * t378 + t375;
t261 = sin(pkin(10));
t405 = cos(pkin(10));
t135 = -t261 * t204 + t206 * t405;
t263 = sin(qJ(5));
t268 = cos(qJ(2));
t377 = qJD(1) * t268;
t337 = qJD(3) + t377;
t360 = t265 * qJDD(1);
t291 = qJD(2) * t337 + t360;
t284 = t291 * t264;
t362 = qJD(1) * qJD(3);
t329 = t265 * t362 - qJDD(2);
t127 = t267 * t329 + t284;
t310 = t329 * t264;
t279 = t291 * t267 - t310;
t275 = -t261 * t127 + t279 * t405;
t276 = -t127 * t405 - t261 * t279;
t312 = t204 * t405 + t261 * t206;
t433 = cos(qJ(5));
t288 = t433 * t312;
t367 = qJD(5) * t263;
t22 = qJD(5) * t288 + t135 * t367 - t263 * t276 - t433 * t275;
t234 = -qJD(3) + t377;
t225 = -qJD(5) + t234;
t86 = t135 * t263 + t288;
t412 = t225 * t86;
t13 = -t22 - t412;
t392 = t267 * t268;
t431 = pkin(3) * t265;
t324 = -qJ(4) * t392 + t431;
t262 = -qJ(4) - pkin(8);
t342 = qJD(3) * t262;
t335 = pkin(2) * t265 - pkin(8) * t268;
t208 = t335 * qJD(1);
t385 = pkin(7) * t351 + t267 * t208;
t470 = -qJD(1) * t324 - qJD(4) * t264 + t267 * t342 - t385;
t188 = t264 * t208;
t368 = qJD(4) * t267;
t395 = t265 * t267;
t398 = t264 * t268;
t469 = t188 + (-pkin(7) * t395 - qJ(4) * t398) * qJD(1) - t264 * t342 - t368;
t420 = t86 ^ 2;
t450 = t135 * t433 - t263 * t312;
t463 = t450 ^ 2;
t468 = -t420 + t463;
t414 = qJD(2) * pkin(2);
t223 = pkin(7) * t378 - t414;
t148 = t204 * pkin(3) + qJD(4) + t223;
t95 = pkin(4) * t312 + t148;
t29 = t86 * pkin(5) - qJ(6) * t450 + t95;
t467 = t29 * t86;
t466 = t86 * t95;
t419 = t450 * t86;
t196 = t261 * t267 + t264 * t405;
t300 = t268 * t196;
t160 = qJD(1) * t300;
t298 = qJD(3) * t196;
t457 = t160 - t298;
t311 = t261 * t264 - t267 * t405;
t301 = t268 * t311;
t161 = qJD(1) * t301;
t438 = qJD(3) * t311;
t465 = t161 - t438;
t44 = pkin(5) * t450 + qJ(6) * t86;
t449 = t261 * t469 + t405 * t470;
t448 = t261 * t470 - t405 * t469;
t248 = t268 * qJDD(1);
t363 = qJD(1) * qJD(2);
t441 = -t265 * t363 + t248;
t195 = qJDD(3) - t441;
t192 = qJDD(5) + t195;
t258 = qJ(3) + pkin(10);
t249 = qJ(5) + t258;
t239 = sin(t249);
t256 = g(3) * t268;
t266 = sin(qJ(1));
t269 = cos(qJ(1));
t332 = g(1) * t269 + g(2) * t266;
t293 = t332 * t265 - t256;
t218 = t262 * t264;
t219 = t262 * t267;
t146 = t405 * t218 + t219 * t261;
t120 = -pkin(9) * t196 + t146;
t147 = t261 * t218 - t405 * t219;
t121 = -pkin(9) * t311 + t147;
t72 = t263 * t120 + t121 * t433;
t464 = t72 * t192 + t239 * t293;
t461 = pkin(4) * t378 + pkin(9) * t465 - t449;
t460 = pkin(9) * t457 + t448;
t459 = t450 * t95;
t23 = qJD(5) * t450 + t263 * t275 - t433 * t276;
t410 = t450 * t225;
t458 = -t23 - t410;
t183 = t192 * pkin(5);
t437 = qJDD(6) - t183;
t456 = t29 * t450 + t437;
t455 = pkin(9) * t135;
t454 = pkin(9) * t312;
t287 = t433 * t311;
t407 = -t196 * t367 + (-qJD(3) - qJD(5)) * t287 + t161 * t433 + t457 * t263;
t303 = t263 * t311;
t347 = qJD(5) * t433;
t353 = t433 * t196;
t406 = qJD(3) * t353 - qJD(5) * t303 - t160 * t433 + t196 * t347 + t263 * t465;
t245 = pkin(7) * t377;
t359 = pkin(3) * t398;
t371 = qJD(3) * t264;
t452 = pkin(3) * t371 - qJD(1) * t359 - t245;
t394 = t265 * t269;
t396 = t265 * t266;
t451 = g(1) * t394 + g(2) * t396 - t256;
t352 = t405 * pkin(3);
t241 = t352 + pkin(4);
t432 = pkin(3) * t261;
t356 = t263 * t432;
t217 = -pkin(2) * t268 - pkin(8) * t265 - pkin(1);
t193 = t217 * qJD(1);
t224 = qJD(2) * pkin(8) + t245;
t142 = t267 * t193 - t224 * t264;
t109 = -qJ(4) * t206 + t142;
t143 = t264 * t193 + t267 * t224;
t110 = -qJ(4) * t204 + t143;
t341 = t405 * t110;
t64 = -t109 * t261 - t341;
t55 = t64 + t454;
t105 = t261 * t110;
t65 = t405 * t109 - t105;
t56 = t65 - t455;
t447 = -qJD(5) * t356 + t241 * t347 - t263 * t55 - t433 * t56;
t315 = t120 * t433 - t263 * t121;
t446 = qJD(5) * t315 - t263 * t461 + t433 * t460;
t445 = qJD(5) * t72 + t263 * t460 + t433 * t461;
t388 = -pkin(4) * t457 + t452;
t179 = t192 * qJ(6);
t216 = t225 * qJD(6);
t443 = t179 - t216;
t384 = t263 * t241 + t433 * t432;
t243 = pkin(7) * t360;
t345 = t268 * t363;
t403 = qJDD(2) * pkin(2);
t181 = pkin(7) * t345 + t243 - t403;
t440 = qJD(3) * pkin(8) * t234 - t181;
t393 = t266 * t268;
t184 = t264 * t393 + t267 * t269;
t391 = t268 * t269;
t186 = -t264 * t391 + t266 * t267;
t439 = -g(1) * t186 + g(2) * t184;
t286 = t265 * t298;
t236 = pkin(7) * t392;
t209 = t335 * qJD(2);
t374 = qJD(2) * t265;
t430 = pkin(7) * t264;
t386 = t267 * t209 + t374 * t430;
t81 = -t265 * t368 + t324 * qJD(2) + (-t236 + (qJ(4) * t265 - t217) * t264) * qJD(3) + t386;
t369 = qJD(3) * t267;
t387 = t264 * t209 + t217 * t369;
t91 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t395 + (-qJD(4) * t265 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t268) * t264 + t387;
t41 = -t261 * t91 + t405 * t81;
t31 = pkin(9) * t286 + (t265 * pkin(4) + pkin(9) * t301) * qJD(2) + t41;
t173 = t311 * t265;
t198 = t267 * t217;
t139 = -qJ(4) * t395 + t198 + (-pkin(3) - t430) * t268;
t383 = t264 * t217 + t236;
t399 = t264 * t265;
t145 = -qJ(4) * t399 + t383;
t92 = t405 * t139 - t145 * t261;
t68 = -pkin(4) * t268 + pkin(9) * t173 + t92;
t302 = t265 * t196;
t93 = t261 * t139 + t405 * t145;
t73 = -pkin(9) * t302 + t93;
t319 = t263 * t68 + t433 * t73;
t277 = -qJD(2) * t300 + t265 * t438;
t42 = t261 * t81 + t405 * t91;
t34 = pkin(9) * t277 + t42;
t436 = -qJD(5) * t319 - t263 * t34 + t31 * t433;
t435 = -0.2e1 * pkin(1);
t428 = g(1) * t266;
t424 = g(2) * t269;
t423 = g(3) * t265;
t422 = t264 * pkin(3);
t251 = t265 * pkin(7);
t253 = t267 * pkin(3);
t242 = t253 + pkin(2);
t418 = -qJ(6) * t378 + t446;
t417 = pkin(5) * t378 + t445;
t132 = t353 - t303;
t416 = pkin(5) * t406 - qJ(6) * t407 - t132 * qJD(6) + t388;
t144 = qJD(1) * t209 + qJDD(1) * t217;
t138 = t267 * t144;
t180 = pkin(7) * t441 + qJDD(2) * pkin(8);
t361 = t264 * qJDD(2);
t43 = -t264 * t180 + t138 - (t361 + (t345 + t360) * t267) * qJ(4) - t206 * qJD(4) + t195 * pkin(3) - t110 * qJD(3);
t308 = t264 * t144 + t267 * t180 + t193 * t369 - t224 * t371;
t51 = -qJ(4) * t127 - qJD(4) * t204 + t308;
t17 = t261 * t43 + t405 * t51;
t101 = -pkin(3) * t234 + t109;
t59 = t405 * t101 - t105;
t45 = -pkin(4) * t234 - t455 + t59;
t60 = t261 * t101 + t341;
t49 = t60 - t454;
t19 = t263 * t45 + t433 * t49;
t413 = t19 * t225;
t409 = qJD(6) + t447;
t408 = qJD(5) * t384 - t263 * t56 + t433 * t55;
t402 = t206 * t234;
t401 = t239 * t265;
t240 = cos(t249);
t400 = t240 * t265;
t397 = t264 * t269;
t390 = t269 * t239;
t18 = -t263 * t49 + t433 * t45;
t389 = qJD(6) - t18;
t211 = pkin(4) * cos(t258) + t253;
t381 = pkin(3) * t399 + t251;
t380 = t269 * pkin(1) + t266 * pkin(7);
t259 = t265 ^ 2;
t379 = -t268 ^ 2 + t259;
t376 = qJD(2) * t204;
t373 = qJD(2) * t268;
t372 = qJD(3) * t204;
t370 = qJD(3) * t265;
t366 = t206 * qJD(2);
t365 = t223 * qJD(3);
t354 = pkin(7) * t373 + qJD(2) * t359 + t369 * t431;
t350 = t234 * t369;
t349 = t234 * t371;
t348 = t264 * t370;
t12 = pkin(9) * t276 + t17;
t16 = -t261 * t51 + t405 * t43;
t8 = t195 * pkin(4) - pkin(9) * t275 + t16;
t344 = t433 * t12 + t263 * t8 + t45 * t347 - t49 * t367;
t343 = t263 * t12 + t49 * t347 + t45 * t367 - t433 * t8;
t154 = t239 * t393 + t240 * t269;
t155 = t240 * t393 - t390;
t340 = -t154 * pkin(5) + qJ(6) * t155;
t156 = -t266 * t240 + t268 * t390;
t157 = t239 * t266 + t240 * t391;
t339 = -t156 * pkin(5) + qJ(6) * t157;
t338 = -qJD(3) * t193 - t180;
t237 = g(1) * t396;
t336 = -g(2) * t394 + t237;
t102 = pkin(3) * t206 + pkin(4) * t135;
t334 = -g(1) * t154 + g(2) * t156;
t333 = g(1) * t155 - g(2) * t157;
t331 = t224 * t369 - t138;
t330 = -pkin(8) * t195 + t365;
t328 = t242 * t268 - t262 * t265;
t325 = qJD(2) * qJD(3) + t360;
t207 = pkin(2) + t211;
t323 = pkin(5) * t240 + qJ(6) * t239 + t207;
t321 = -t263 * t73 + t433 * t68;
t317 = -pkin(7) * qJDD(2) + t363 * t435;
t316 = t263 * t31 + t433 * t34 + t68 * t347 - t367 * t73;
t314 = t195 * t264 - t350;
t313 = t267 * t195 + t349;
t309 = qJDD(1) * t267 - t264 * t362;
t271 = qJD(1) ^ 2;
t307 = pkin(1) * t271 + t332;
t306 = t241 * t433 - t356;
t270 = qJD(2) ^ 2;
t305 = pkin(7) * t270 + qJDD(1) * t435 + t424;
t304 = t315 * t192 + t240 * t451;
t100 = t127 * pkin(3) + qJDD(4) + t181;
t296 = t22 - t412;
t295 = g(1) * t156 + g(2) * t154 + g(3) * t401 - t343;
t294 = g(1) * t157 + g(2) * t155 + g(3) * t400 - t344;
t290 = t263 * t302;
t162 = pkin(4) * t311 - t242;
t285 = t265 * t353;
t282 = -t18 * t225 + t294;
t140 = pkin(4) * t302 + t381;
t281 = t225 * t408 + t295;
t280 = -t295 + t456;
t278 = -qJD(2) * t301 - t286;
t94 = -pkin(4) * t277 + t354;
t52 = -pkin(4) * t276 + t100;
t3 = t23 * pkin(5) + t22 * qJ(6) - qJD(6) * t450 + t52;
t272 = t23 - t410;
t257 = -pkin(9) + t262;
t254 = t269 * pkin(7);
t213 = qJ(6) * t400;
t210 = pkin(4) * sin(t258) + t422;
t187 = t264 * t266 + t267 * t391;
t185 = -t266 * t392 + t397;
t174 = -pkin(5) - t306;
t172 = qJ(6) + t384;
t131 = t196 * t263 + t287;
t115 = -t173 * t433 - t290;
t114 = -t173 * t263 + t285;
t69 = t131 * pkin(5) - t132 * qJ(6) + t162;
t57 = t114 * pkin(5) - t115 * qJ(6) + t140;
t54 = -qJD(5) * t290 - t173 * t347 + t263 * t278 - t277 * t433;
t53 = qJD(5) * t285 - t173 * t367 - t263 * t277 - t278 * t433;
t36 = t268 * pkin(5) - t321;
t35 = -qJ(6) * t268 + t319;
t33 = t102 + t44;
t15 = -t225 * qJ(6) + t19;
t14 = t225 * pkin(5) + t389;
t9 = t54 * pkin(5) + t53 * qJ(6) - t115 * qJD(6) + t94;
t5 = -pkin(5) * t374 - t436;
t4 = qJ(6) * t374 - qJD(6) * t268 + t316;
t2 = t343 + t437;
t1 = t344 + t443;
t6 = [qJDD(1), -t424 + t428, t332, qJDD(1) * t259 + 0.2e1 * t265 * t345, 0.2e1 * t248 * t265 - 0.2e1 * t363 * t379, qJDD(2) * t265 + t268 * t270, qJDD(2) * t268 - t265 * t270, 0, t317 * t265 + (-t305 + t428) * t268, t265 * t305 + t268 * t317 - t237, -t206 * t348 + (t268 * t366 + (t325 + t345) * t395 - t265 * t310) * t267 (-t204 * t267 - t206 * t264) * t373 + ((-t206 * qJD(3) - t127) * t267 + (-t279 + t372) * t264) * t265 (-t361 + (-t234 - t337) * t364) * t268 + (-t268 * t309 + t313 + t366) * t265 (t234 * t375 + t127) * t268 + (-t314 - t376) * t265, -t195 * t268 - t234 * t374 -(-t217 * t371 + t386) * t234 + t198 * t195 - g(1) * t185 - g(2) * t187 + ((t350 + t376) * pkin(7) + (-pkin(7) * t195 + qJD(2) * t223 - t338) * t264 + t331) * t268 + (pkin(7) * t127 + qJD(2) * t142 + t181 * t264 + t267 * t365) * t265, t387 * t234 - t383 * t195 - g(1) * t184 - g(2) * t186 + (t223 * t364 + (-t349 + t366) * pkin(7) + t308) * t268 + (-t264 * t365 - t143 * qJD(2) + t181 * t267 + (t361 + t309 * t265 + (-t234 + t337) * t364) * pkin(7)) * t265, -t42 * t312 + t93 * t276 - t17 * t302 + t60 * t277 - t41 * t135 - t92 * t275 + t16 * t173 + t59 * (t196 * t370 + t311 * t373) + t336, t17 * t93 + t60 * t42 + t16 * t92 + t59 * t41 + t100 * t381 + t148 * t354 - g(1) * (pkin(3) * t397 + t254) - g(2) * (t242 * t391 - t262 * t394 + t380) + (-g(1) * (-pkin(1) - t328) - g(2) * t422) * t266, -t115 * t22 - t450 * t53, t114 * t22 - t115 * t23 - t450 * t54 + t53 * t86, t115 * t192 + t22 * t268 + t225 * t53 + t374 * t450, -t114 * t192 + t225 * t54 + t23 * t268 - t374 * t86, -t192 * t268 - t225 * t374, t52 * t114 + t140 * t23 + t18 * t374 + t321 * t192 - t225 * t436 + t343 * t268 + t95 * t54 + t94 * t86 + t333, t52 * t115 - t140 * t22 - t19 * t374 - t192 * t319 + t225 * t316 + t268 * t344 + t450 * t94 - t95 * t53 + t334, t114 * t3 - t14 * t374 - t192 * t36 + t2 * t268 + t225 * t5 + t23 * t57 + t29 * t54 + t86 * t9 + t333, -t1 * t114 + t115 * t2 - t14 * t53 - t15 * t54 - t22 * t36 - t23 * t35 - t4 * t86 + t450 * t5 + t336, -t1 * t268 - t115 * t3 + t15 * t374 + t192 * t35 + t22 * t57 - t225 * t4 + t29 * t53 - t450 * t9 - t334, t1 * t35 + t15 * t4 + t3 * t57 + t29 * t9 + t2 * t36 + t14 * t5 - g(1) * (-pkin(5) * t155 - qJ(6) * t154 + t210 * t269 + t254) - g(2) * (pkin(5) * t157 + qJ(6) * t156 + t207 * t391 - t257 * t394 + t380) + (-g(1) * (-t207 * t268 + t257 * t265 - pkin(1)) - g(2) * t210) * t266; 0, 0, 0, -t265 * t271 * t268, t379 * t271, t360, t248, qJDD(2), t265 * t307 - t243 - t256, t423 + (-pkin(7) * qJDD(1) + t307) * t268, -t329 * t264 ^ 2 + (t284 - t402) * t267 (-t127 + t402) * t264 + (-t372 + t361 + t325 * t267 + (-t348 + (t204 + t364) * t268) * qJD(1)) * t267 (-t206 * t265 + t234 * t392) * qJD(1) + t314 (t204 * t265 - t234 * t398) * qJD(1) + t313, t234 * t378, -pkin(2) * t127 + t385 * t234 + t330 * t264 + (-t142 * t265 + (-pkin(7) * t204 - t223 * t264) * t268) * qJD(1) + (t293 + t440) * t267, -t188 * t234 + (-t268 * pkin(7) * t206 + t143 * t265) * qJD(1) + (-pkin(2) * t325 + (t234 * t251 + (-t223 - t414) * t268) * qJD(1) + t330) * t267 + (-t403 + t256 + (pkin(2) * t362 - t332) * t265 - t440) * t264, -g(1) * t391 - g(2) * t393 - t449 * t135 - t146 * t275 + t147 * t276 - t16 * t196 - t17 * t311 - t448 * t312 + t457 * t60 - t465 * t59 - t423, t17 * t147 + t16 * t146 - t100 * t242 - g(3) * t328 + t448 * t60 + t449 * t59 + t452 * t148 + t332 * (t242 * t265 + t262 * t268) -t132 * t22 + t407 * t450, t131 * t22 - t132 * t23 - t406 * t450 - t407 * t86, t132 * t192 - t225 * t407 - t378 * t450, -t131 * t192 + t225 * t406 + t378 * t86, t225 * t378, t52 * t131 + t162 * t23 - t18 * t378 + t225 * t445 + t388 * t86 + t406 * t95 + t304, t52 * t132 - t162 * t22 + t19 * t378 + t225 * t446 + t388 * t450 + t407 * t95 - t464, t131 * t3 + t14 * t378 + t225 * t417 + t23 * t69 + t29 * t406 + t416 * t86 + t304, -t1 * t131 + t132 * t2 + t14 * t407 - t15 * t406 + t22 * t315 - t23 * t72 - t268 * t332 + t417 * t450 - t418 * t86 - t423, -t132 * t3 - t15 * t378 + t22 * t69 - t225 * t418 - t29 * t407 - t416 * t450 + t464, t1 * t72 - t2 * t315 + t3 * t69 + t416 * t29 + t418 * t15 + t417 * t14 + (-g(3) * t323 + t257 * t332) * t268 + (g(3) * t257 + t323 * t332) * t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206 * t204, -t204 ^ 2 + t206 ^ 2, -t204 * t234 + t279, -t402 - t127, t195, -t143 * t234 - t206 * t223 + (t338 + t423) * t264 - t331 + t439, g(1) * t187 - g(2) * t185 + g(3) * t395 - t142 * t234 + t204 * t223 - t308, -t275 * t352 + t276 * t432 + (-t59 + t65) * t312 + (t60 + t64) * t135, -t59 * t64 - t60 * t65 + (g(3) * t399 - t148 * t206 + t16 * t405 + t17 * t261 + t439) * pkin(3), t419, t468, t13, t458, t192, -t102 * t86 + t192 * t306 + t281 - t459, -t102 * t450 - t384 * t192 + t225 * t447 + t294 + t466, -t174 * t192 - t33 * t86 + t281 - t456, -t172 * t23 - t174 * t22 + (t15 + t408) * t450 + (t14 - t409) * t86, t172 * t192 - t225 * t409 + t33 * t450 - t294 + t443 - t467, t1 * t172 + t2 * t174 - t29 * t33 - g(1) * (-t210 * t391 + t211 * t266 + t339) - g(2) * (-t210 * t393 - t211 * t269 + t340) - g(3) * (t213 + (-pkin(5) * t239 - t210) * t265) + t409 * t15 + t408 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135 ^ 2 - t312 ^ 2, t135 * t59 + t312 * t60 + t100 - t293, 0, 0, 0, 0, 0, t272, -t296, t272, -t420 - t463, t296, -t14 * t450 + t15 * t86 + t3 - t451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, t468, t13, t458, t192, t295 - t413 - t459, t282 + t466, -t44 * t86 + t183 - t280 - t413, pkin(5) * t22 - qJ(6) * t23 + (t15 - t19) * t450 + (t14 - t389) * t86, t44 * t450 + 0.2e1 * t179 - 0.2e1 * t216 - t282 - t467, t1 * qJ(6) - t2 * pkin(5) - t29 * t44 - t14 * t19 - g(1) * t339 - g(2) * t340 - g(3) * (-pkin(5) * t401 + t213) + t389 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192 + t419, t13, -t225 ^ 2 - t463, t15 * t225 + t280;];
tau_reg  = t6;
