% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRRP6
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 10:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:17:03
% EndTime: 2019-05-05 10:17:29
% DurationCPUTime: 12.29s
% Computational Cost: add. (38012->507), mult. (78664->736), div. (0->0), fcn. (63843->14), ass. (0->352)
t283 = cos(pkin(7));
t277 = qJD(2) * t283 + qJD(3);
t285 = sin(qJ(4));
t289 = cos(qJ(4));
t286 = sin(qJ(3));
t281 = sin(pkin(7));
t373 = qJD(2) * t281;
t357 = t286 * t373;
t255 = -t289 * t277 + t285 * t357;
t251 = qJD(5) + t255;
t421 = t251 ^ 2;
t257 = t277 * t285 + t289 * t357;
t290 = cos(qJ(3));
t372 = qJD(2) * t290;
t356 = t281 * t372;
t273 = -qJD(4) + t356;
t284 = sin(qJ(5));
t288 = cos(qJ(5));
t236 = t257 * t284 + t288 * t273;
t422 = t236 ^ 2;
t205 = t422 - t421;
t367 = qJDD(2) * t286;
t368 = qJD(2) * qJD(3);
t262 = (t290 * t368 + t367) * t281;
t276 = qJDD(2) * t283 + qJDD(3);
t352 = t285 * t262 - t289 * t276;
t219 = -qJD(4) * t257 - t352;
t218 = qJDD(5) - t219;
t238 = t257 * t288 - t273 * t284;
t387 = t238 * t236;
t159 = -t387 - t218;
t396 = t159 * t288;
t121 = -t205 * t284 + t396;
t397 = t159 * t284;
t125 = -t205 * t288 - t397;
t213 = t238 * t251;
t322 = -t289 * t262 - t285 * t276;
t220 = -qJD(4) * t255 - t322;
t366 = qJDD(2) * t290;
t316 = t286 * t368 - t366;
t310 = t316 * t281;
t306 = qJDD(4) + t310;
t353 = -t284 * t220 + t288 * t306;
t318 = qJD(5) * t238 - t353;
t139 = -t213 + t318;
t82 = t125 * t285 - t139 * t289;
t505 = t283 * t82 + (t286 * (t125 * t289 + t139 * t285) - t290 * t121) * t281;
t235 = t238 ^ 2;
t438 = t235 - t422;
t300 = -t288 * t220 - t284 * t306;
t296 = -t236 * qJD(5) - t300;
t388 = t236 * t251;
t434 = -t388 + t296;
t400 = t434 * t284;
t440 = t213 + t318;
t90 = t440 * t288 + t400;
t66 = t285 * t90 + t289 * t438;
t89 = -t440 * t284 + t288 * t434;
t503 = t281 * (t286 * (-t285 * t438 + t289 * t90) + t290 * t89) + t283 * t66;
t435 = -t387 + t218;
t395 = t435 * t284;
t432 = -t421 - t422;
t447 = t288 * t432 - t395;
t463 = t285 * t447 - t289 * t440;
t394 = t435 * t288;
t448 = t284 * t432 + t394;
t462 = t285 * t440 + t289 * t447;
t476 = t286 * t462 - t290 * t448;
t485 = -t281 * t463 + t283 * t476;
t502 = pkin(2) * t485;
t439 = -t235 - t421;
t109 = -t288 * t439 - t397;
t501 = pkin(3) * t109;
t500 = pkin(4) * t109;
t499 = pkin(11) * t109;
t111 = t284 * t439 - t396;
t498 = pkin(11) * t111;
t497 = t109 * t286;
t496 = t109 * t290;
t495 = t111 * t285;
t494 = t111 * t289;
t282 = sin(pkin(6));
t287 = sin(qJ(2));
t410 = cos(pkin(6));
t417 = cos(qJ(2));
t475 = t286 * t448 + t290 * t462;
t490 = t410 * (t281 * t476 + t283 * t463) + (t287 * t475 + t417 * t485) * t282;
t489 = pkin(9) * t475;
t433 = t388 + t296;
t206 = -t235 + t421;
t464 = -t206 * t284 + t394;
t465 = t206 * t288 + t395;
t477 = t285 * t464 - t289 * t433;
t486 = t283 * t477 + (t286 * (t285 * t433 + t289 * t464) - t290 * t465) * t281;
t483 = pkin(3) * t463;
t482 = pkin(10) * t463;
t478 = -pkin(3) * t448 + pkin(10) * t462;
t472 = pkin(4) * t448;
t471 = pkin(11) * t447;
t470 = pkin(11) * t448;
t437 = t235 + t422;
t461 = pkin(4) * t437;
t460 = qJ(6) * t434;
t457 = t285 * t437;
t452 = t289 * t437;
t408 = sin(pkin(12));
t409 = cos(pkin(12));
t313 = -g(1) * t409 - g(2) * t408;
t312 = g(1) * t408 - g(2) * t409;
t305 = t410 * t312;
t376 = -g(3) + qJDD(1);
t436 = t282 * t376 + t305;
t231 = -t287 * t313 + t417 * t436;
t291 = qJD(2) ^ 2;
t413 = pkin(9) * t281;
t295 = qJDD(2) * pkin(2) + t291 * t413 + t231;
t299 = -t282 * t312 + t376 * t410;
t449 = t281 * t299 + t283 * t295;
t386 = t251 * t284;
t133 = t236 * t386 - t288 * t318;
t385 = t251 * t288;
t363 = t236 * t385;
t319 = t284 * t318 + t363;
t365 = t285 * t387;
t364 = t289 * t387;
t425 = t285 * t319 + t364;
t446 = t283 * t425 + (-t290 * t133 + t286 * (t289 * t319 - t365)) * t281;
t203 = t238 * t386;
t345 = t203 - t363;
t427 = -t289 * t218 + t285 * t345;
t441 = (t236 * t284 + t238 * t288) * t251;
t445 = t283 * t427 + (t290 * t441 + t286 * (t218 * t285 + t289 * t345)) * t281;
t278 = t281 ^ 2;
t444 = t278 * t290;
t384 = t257 * t255;
t302 = t306 - t384;
t443 = t285 * t302;
t442 = t289 * t302;
t245 = t255 * t273;
t187 = t220 + t245;
t193 = pkin(5) * t236 - qJ(6) * t238;
t293 = t449 * t286;
t232 = t287 * t436 + t417 * t313;
t224 = -t291 * pkin(2) + qJDD(2) * t413 + t232;
t416 = pkin(3) * t290;
t348 = -pkin(10) * t286 - t416;
t341 = t348 * t373 ^ 2 + t224;
t419 = t277 ^ 2;
t152 = -pkin(3) * t419 + t276 * pkin(10) + t290 * t341 + t293;
t297 = t283 * t299;
t350 = qJD(2) * (qJD(3) + t277);
t292 = t297 - t262 * pkin(10) + (-t277 * pkin(10) * t372 + (t286 * t350 - t366) * pkin(3) - t295) * t281;
t101 = t289 * t152 + t285 * t292;
t225 = pkin(4) * t255 - pkin(11) * t257;
t420 = t273 ^ 2;
t75 = -pkin(4) * t420 + pkin(11) * t306 - t255 * t225 + t101;
t374 = t449 * t290;
t151 = -t276 * pkin(3) - t419 * pkin(10) + t286 * t341 - t374;
t97 = -t187 * pkin(11) + (-t257 * t273 - t219) * pkin(4) + t151;
t52 = t284 * t97 + t288 * t75;
t351 = t218 * qJ(6) - t236 * t193 + t52;
t431 = -pkin(5) * (t439 + t421) - qJ(6) * t159 + t351;
t371 = qJD(6) * t251;
t247 = 0.2e1 * t371;
t344 = t247 + t351;
t42 = -pkin(5) * t421 + t344;
t51 = t284 * t75 - t288 * t97;
t43 = -t218 * pkin(5) - qJ(6) * t421 + t193 * t238 + qJDD(6) + t51;
t21 = t284 * t43 + t288 * t42;
t355 = qJ(6) * t284 + pkin(4);
t414 = pkin(5) * t288;
t100 = t285 * t152 - t289 * t292;
t74 = -t306 * pkin(4) - t420 * pkin(11) + t257 * t225 + t100;
t304 = t318 * pkin(5) - t460 + t74;
t49 = (pkin(5) * t251 - 0.2e1 * qJD(6)) * t238 + t304;
t430 = -(t355 + t414) * t49 + pkin(11) * t21;
t303 = 0.2e1 * qJD(6) * t238 - t304;
t40 = (-t440 - t213) * pkin(5) + t303;
t429 = t288 * t40 - t355 * t440 + t471;
t39 = -pkin(5) * t213 + t303 + t460;
t428 = t498 + t434 * (pkin(4) + t414) + t284 * t39;
t184 = (qJD(4) + t273) * t257 + t352;
t136 = t238 * t385 + t284 * t296;
t137 = t288 * t296 - t203;
t346 = t285 * t137 - t364;
t379 = t281 * t290;
t380 = t281 * t286;
t426 = t283 * t346 + (t289 * t137 + t365) * t380 - t136 * t379;
t253 = t255 ^ 2;
t254 = t257 ^ 2;
t415 = pkin(4) * t285;
t412 = t284 * t74;
t411 = t288 * t74;
t407 = qJ(6) * t288;
t402 = t433 * t284;
t401 = t433 * t288;
t399 = t151 * t285;
t398 = t151 * t289;
t208 = -t306 - t384;
t391 = t208 * t285;
t390 = t208 * t289;
t383 = t273 * t285;
t382 = t273 * t289;
t381 = t278 * t291;
t272 = t286 * t290 * t381;
t259 = t272 + t276;
t378 = t286 * t259;
t260 = -t272 + t276;
t377 = t290 * t260;
t370 = qJD(4) - t273;
t279 = t286 ^ 2;
t362 = t279 * t381;
t280 = t290 ^ 2;
t361 = t280 * t381;
t360 = t290 * t384;
t359 = -pkin(4) * t289 - pkin(3);
t24 = t284 * t51 + t288 * t52;
t58 = t100 * t285 + t289 * t101;
t349 = -pkin(4) * t74 + pkin(11) * t24;
t347 = -pkin(5) * t43 + qJ(6) * t42;
t340 = -pkin(5) * t433 - qJ(6) * t139;
t14 = t21 * t289 + t285 * t49;
t20 = t284 * t42 - t288 * t43;
t339 = t14 * t286 - t20 * t290;
t19 = t24 * t289 + t285 * t74;
t23 = t284 * t52 - t288 * t51;
t338 = t19 * t286 - t23 * t290;
t141 = (-qJD(5) + t251) * t238 + t353;
t91 = t141 * t288 + t402;
t64 = t289 * t91 - t457;
t87 = t141 * t284 - t401;
t337 = t286 * t64 - t290 * t87;
t92 = -t139 * t288 + t402;
t65 = t289 * t92 - t457;
t88 = -t139 * t284 - t401;
t336 = t286 * t65 - t290 * t88;
t71 = -t285 * t434 + t494;
t335 = t286 * t71 - t496;
t146 = (qJD(5) + t251) * t236 + t300;
t78 = -t146 * t285 - t494;
t333 = t286 * t78 + t496;
t331 = -t151 * t290 + t286 * t58;
t57 = -t100 * t289 + t101 * t285;
t188 = t220 - t245;
t149 = -t184 * t289 + t188 * t285;
t202 = t253 + t254;
t330 = t149 * t286 + t202 * t290;
t160 = t286 * t224 - t374;
t161 = t290 * t224 + t293;
t329 = -t160 * t290 + t161 * t286;
t107 = t160 * t286 + t161 * t290;
t221 = -t420 - t253;
t168 = t221 * t289 - t443;
t185 = -t257 * t370 - t352;
t328 = t168 * t286 + t185 * t290;
t228 = -t254 - t420;
t175 = -t228 * t285 + t390;
t189 = t255 * t370 + t322;
t327 = t175 * t286 + t189 * t290;
t267 = t277 * t356;
t242 = -t267 + t262;
t266 = t277 * t357;
t243 = t266 - t310;
t325 = -t242 * t290 + t243 * t286;
t250 = -t362 - t419;
t324 = t250 * t290 - t260 * t286;
t263 = -t361 - t419;
t323 = t259 * t290 + t263 * t286;
t321 = qJD(2) * t277 - t283 * t291;
t315 = -pkin(4) * t440 - t411 + t471;
t314 = pkin(4) * t146 + t412 - t498;
t311 = pkin(11) * t91 + t24 + t461;
t32 = (t437 - t421) * pkin(5) + t344;
t35 = qJ(6) * t437 + t43;
t308 = pkin(11) * t92 + t284 * t35 + t288 * t32 + t461;
t301 = pkin(5) * t435 + qJ(6) * t432 - t43;
t265 = (-t279 - t280) * t381;
t264 = (t279 - t280) * t381;
t244 = -t266 - t310;
t241 = (t290 * t350 + t367) * t281;
t240 = -t254 + t420;
t239 = t253 - t420;
t229 = t263 * t290 - t378;
t227 = t254 - t253;
t223 = -t250 * t286 - t377;
t201 = (t255 * t285 + t257 * t289) * t273;
t200 = t242 * t286 + t243 * t290;
t199 = t281 * t295 - t297;
t192 = t281 * t244 + t283 * t323;
t183 = -t281 * t241 + t283 * t324;
t182 = t220 * t285 - t257 * t382;
t181 = t219 * t289 - t255 * t383;
t180 = -t281 * t265 + t283 * t325;
t177 = t239 * t285 - t390;
t176 = t240 * t289 + t443;
t174 = t228 * t289 + t391;
t167 = t221 * t285 + t442;
t148 = -t184 * t285 - t188 * t289;
t147 = t185 * t285 + t187 * t289;
t118 = t175 * t290 - t189 * t286;
t117 = t168 * t290 - t185 * t286;
t108 = t149 * t290 - t202 * t286;
t102 = -t281 * t174 + t283 * t327;
t98 = -t281 * t167 + t283 * t328;
t94 = t281 * t199 + t283 * t329;
t85 = pkin(3) * t189 + pkin(10) * t175 + t399;
t84 = pkin(3) * t185 + pkin(10) * t168 - t398;
t76 = t146 * t289 - t495;
t69 = t289 * t434 + t495;
t68 = -t281 * t148 + t283 * t330;
t63 = t285 * t92 + t452;
t62 = t285 * t91 + t452;
t60 = t290 * t78 - t497;
t59 = t411 + t499;
t56 = t412 - t470;
t54 = t290 * t71 + t497;
t53 = -pkin(4) * t88 - t340;
t48 = -pkin(3) * t151 + pkin(10) * t58;
t47 = t151 * t286 + t290 * t58;
t46 = t286 * t88 + t290 * t65;
t45 = t286 * t87 + t290 * t64;
t44 = pkin(3) * t202 + pkin(10) * t149 + t58;
t41 = t52 + t500;
t38 = t51 - t472;
t36 = -t281 * t76 + t283 * t333;
t33 = -t281 * t69 + t283 * t335;
t31 = -t281 * t63 + t283 * t336;
t30 = -t281 * t62 + t283 * t337;
t29 = -t301 - t472;
t28 = -t284 * t40 - t407 * t440 - t470;
t27 = -pkin(5) * t400 + t288 * t39 - t499;
t26 = -0.2e1 * t371 - t431 - t500;
t25 = -t281 * t57 + t283 * t331;
t22 = -pkin(11) * t87 - t23;
t18 = t24 * t285 - t289 * t74;
t17 = pkin(10) * t78 + t285 * t59 + t289 * t41 + t501;
t16 = t285 * t56 + t289 * t38 + t478;
t15 = -pkin(11) * t88 - t284 * t32 + t288 * t35;
t13 = t21 * t285 - t289 * t49;
t12 = pkin(10) * t64 + t285 * t22 + t359 * t87;
t11 = t28 * t285 + t289 * t29 + t478;
t10 = t19 * t290 + t23 * t286;
t9 = pkin(10) * t71 + t26 * t289 + t27 * t285 - t501;
t8 = -pkin(11) * t20 + (pkin(5) * t284 - t407) * t49;
t7 = -pkin(4) * t20 - t347;
t6 = -pkin(3) * t88 + pkin(10) * t65 + t15 * t285 + t289 * t53;
t5 = t14 * t290 + t20 * t286;
t4 = -t281 * t18 + t283 * t338;
t3 = pkin(10) * t19 + (-pkin(11) * t285 + t359) * t23;
t2 = -t281 * t13 + t283 * t339;
t1 = -pkin(3) * t20 + pkin(10) * t14 + t285 * t8 + t289 * t7;
t34 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t376, 0, 0, 0, 0, 0, 0, (qJDD(2) * t417 - t287 * t291) * t282, (-qJDD(2) * t287 - t291 * t417) * t282, 0, t410 ^ 2 * t376 + (t231 * t417 + t287 * t232 - t305) * t282, 0, 0, 0, 0, 0, 0, t410 * (-t283 * t244 + t281 * t323) + (t192 * t417 + t287 * t229) * t282, t410 * (t283 * t241 + t281 * t324) + (t183 * t417 + t287 * t223) * t282, t410 * (t283 * t265 + t281 * t325) + (t180 * t417 + t287 * t200) * t282, t410 * (-t283 * t199 + t281 * t329) + (t287 * t107 + t417 * t94) * t282, 0, 0, 0, 0, 0, 0, t410 * (t283 * t167 + t281 * t328) + (t287 * t117 + t417 * t98) * t282, t410 * (t283 * t174 + t281 * t327) + (t102 * t417 + t287 * t118) * t282, t410 * (t283 * t148 + t281 * t330) + (t287 * t108 + t417 * t68) * t282, t410 * (t281 * t331 + t283 * t57) + (t25 * t417 + t287 * t47) * t282, 0, 0, 0, 0, 0, 0, t490, t410 * (t281 * t333 + t283 * t76) + (t287 * t60 + t36 * t417) * t282, t410 * (t281 * t337 + t283 * t62) + (t287 * t45 + t30 * t417) * t282, t410 * (t283 * t18 + t281 * t338) + (t287 * t10 + t4 * t417) * t282, 0, 0, 0, 0, 0, 0, t490, t410 * (t281 * t336 + t283 * t63) + (t287 * t46 + t31 * t417) * t282, t410 * (t281 * t335 + t283 * t69) + (t287 * t54 + t33 * t417) * t282, t410 * (t283 * t13 + t281 * t339) + (t2 * t417 + t287 * t5) * t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t231, -t232, 0, 0, (t262 * t281 + t321 * t444) * t286, t283 * t264 + (t286 * t244 + (t262 + t267) * t290) * t281, t283 * t242 + (t378 + t290 * (-t362 + t419)) * t281, (-t286 * t321 - t316) * t444, t283 * t243 + (t286 * (t361 - t419) + t377) * t281, t283 * t276, pkin(2) * t192 - t283 * t160 + (pkin(9) * t229 + t199 * t290) * t281, pkin(2) * t183 - t283 * t161 + (pkin(9) * t223 - t199 * t286) * t281, pkin(2) * t180 + (pkin(9) * t200 + t107) * t281, pkin(2) * t94 + t107 * t413, t283 * t182 + (t286 * (t220 * t289 + t257 * t383) - t360) * t281, t283 * t147 + (t286 * (t185 * t289 - t187 * t285) - t290 * t227) * t281, t283 * t176 + (t286 * (-t240 * t285 + t442) - t290 * t188) * t281, t283 * t181 + (t286 * (-t219 * t285 - t255 * t382) + t360) * t281, t283 * t177 + (t286 * (t239 * t289 + t391) + t290 * t184) * t281, -t306 * t379 + t283 * t201 + (t255 * t289 - t257 * t285) * t273 * t380, pkin(2) * t98 + t283 * t84 + (t286 * (-pkin(10) * t167 + t399) + t290 * (-pkin(3) * t167 + t100) + pkin(9) * t117) * t281, pkin(2) * t102 + t283 * t85 + (t286 * (-pkin(10) * t174 + t398) + t290 * (-pkin(3) * t174 + t101) + pkin(9) * t118) * t281, pkin(2) * t68 + t283 * t44 + (t286 * (-pkin(10) * t148 - t57) - t148 * t416 + pkin(9) * t108) * t281, pkin(2) * t25 + t283 * t48 + (pkin(9) * t47 + t348 * t57) * t281, t426, -t503, t486, t446, -t505, t445, t502 + t283 * t16 + (t286 * (-t285 * t38 + t289 * t56 - t482) + t290 * (-t315 - t483) + t489) * t281, pkin(2) * t36 + t283 * t17 + (t286 * (-pkin(10) * t76 - t285 * t41 + t289 * t59) + t290 * (-pkin(3) * t76 - t314) + pkin(9) * t60) * t281, pkin(2) * t30 + t283 * t12 + (t286 * (-pkin(10) * t62 + t22 * t289 + t415 * t87) + t290 * (-pkin(3) * t62 - t311) + pkin(9) * t45) * t281, pkin(2) * t4 + t283 * t3 + (t286 * (-pkin(10) * t18 + (-pkin(11) * t289 + t415) * t23) + t290 * (-pkin(3) * t18 - t349) + pkin(9) * t10) * t281, t426, t486, t503, t445, t505, t446, t502 + t283 * t11 + (t286 * (t28 * t289 - t285 * t29 - t482) + t290 * (-t429 - t483) + t489) * t281, pkin(2) * t31 + t283 * t6 + (t286 * (-pkin(10) * t63 + t15 * t289 - t285 * t53) + t290 * (-pkin(3) * t63 - t308) + pkin(9) * t46) * t281, pkin(2) * t33 + t283 * t9 + (t286 * (-pkin(10) * t69 - t26 * t285 + t27 * t289) + t290 * (-pkin(3) * t69 - t428) + pkin(9) * t54) * t281, pkin(2) * t2 + t283 * t1 + (t286 * (-pkin(10) * t13 - t285 * t7 + t289 * t8) + pkin(9) * t5 + (-pkin(3) * t13 - t430) * t290) * t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t272, t264, t242, t272, t243, t276, -t160, -t161, 0, 0, t182, t147, t176, t181, t177, t201, t84, t85, t44, t48, t346, -t66, t477, t425, -t82, t427, t16, t17, t12, t3, t346, t477, t66, t427, t82, t425, t11, t6, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t384, t227, t188, -t384, -t184, t306, -t100, -t101, 0, 0, t136, t89, t465, t133, -t121, -t441, t315, t314, t311, t349, t136, t465, -t89, -t441, t121, t133, t429, t308, t428, t430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t387, t438, t433, -t387, -t139, t218, -t51, -t52, 0, 0, t387, t433, -t438, t218, t139, -t387, t301, t340, t247 + t431, t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t435, t433, t439, t43;];
tauJ_reg  = t34;
