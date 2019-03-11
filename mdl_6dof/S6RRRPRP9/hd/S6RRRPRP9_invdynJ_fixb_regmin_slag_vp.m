% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:19
% EndTime: 2019-03-09 17:27:38
% DurationCPUTime: 8.52s
% Computational Cost: add. (7414->680), mult. (15612->809), div. (0->0), fcn. (10606->8), ass. (0->295)
t254 = sin(qJ(2));
t239 = t254 * pkin(8);
t258 = cos(qJ(2));
t244 = t258 * pkin(2);
t349 = -pkin(1) - t244;
t300 = t349 - t239;
t151 = t300 * qJD(1);
t370 = qJD(1) * t258;
t230 = pkin(7) * t370;
t186 = qJD(2) * pkin(8) + t230;
t253 = sin(qJ(3));
t257 = cos(qJ(3));
t96 = t257 * t151 - t253 * t186;
t385 = qJD(4) - t96;
t357 = t257 * qJD(2);
t371 = qJD(1) * t254;
t163 = t253 * t371 - t357;
t346 = t257 * t371;
t368 = qJD(2) * t253;
t165 = t346 + t368;
t252 = sin(qJ(5));
t256 = cos(qJ(5));
t303 = -t256 * t163 + t165 * t252;
t251 = qJD(2) * pkin(2);
t185 = pkin(7) * t371 - t251;
t72 = t163 * pkin(3) - t165 * qJ(4) + t185;
t54 = -pkin(4) * t163 - t72;
t91 = t163 * t252 + t165 * t256;
t19 = pkin(5) * t303 - qJ(6) * t91 + t54;
t233 = t258 * qJDD(1);
t356 = qJD(1) * qJD(2);
t441 = -t254 * t356 + t233;
t161 = qJDD(3) - t441;
t319 = pkin(2) * t254 - pkin(8) * t258;
t176 = t319 * qJD(2);
t104 = qJD(1) * t176 + qJDD(1) * t300;
t140 = pkin(7) * t441 + qJDD(2) * pkin(8);
t362 = qJD(3) * t257;
t364 = qJD(3) * t253;
t326 = -t257 * t104 + t253 * t140 + t151 * t364 + t186 * t362;
t306 = qJDD(4) + t326;
t431 = pkin(3) + pkin(4);
t354 = t254 * qJDD(1);
t198 = t257 * t354;
t335 = t258 * t356;
t363 = qJD(3) * t254;
t340 = t253 * t363;
t355 = qJD(2) * qJD(3);
t84 = qJD(1) * t340 - t253 * qJDD(2) - t198 + (-t335 - t355) * t257;
t10 = pkin(9) * t84 - t161 * t431 + t306;
t294 = t354 + t355;
t275 = t294 + t335;
t232 = t257 * qJDD(2);
t334 = qJD(1) * t362;
t291 = t254 * t334 - t232;
t353 = t253 * t104 + t257 * t140 + t151 * t362;
t149 = t161 * qJ(4);
t208 = -qJD(3) + t370;
t189 = t208 * qJD(4);
t443 = t149 - t189;
t13 = t291 * pkin(9) + (pkin(9) * t275 - qJD(3) * t186) * t253 + t353 + t443;
t360 = qJD(5) * t256;
t361 = qJD(5) * t252;
t462 = -pkin(9) * t165 + t385;
t42 = t208 * t431 + t462;
t191 = t208 * qJ(4);
t97 = t253 * t151 + t257 * t186;
t58 = pkin(9) * t163 + t97;
t52 = -t191 + t58;
t333 = -t252 * t10 - t256 * t13 - t42 * t360 + t52 * t361;
t167 = t252 * t253 + t256 * t257;
t446 = t167 * t254;
t259 = cos(qJ(1));
t389 = t257 * t259;
t255 = sin(qJ(1));
t391 = t255 * t258;
t143 = t253 * t391 + t389;
t386 = t259 * t253;
t390 = t257 * t258;
t144 = t255 * t390 - t386;
t77 = t143 * t252 + t144 * t256;
t392 = t255 * t257;
t145 = t258 * t386 - t392;
t387 = t258 * t259;
t146 = t253 * t255 + t257 * t387;
t83 = t145 * t252 + t146 * t256;
t273 = g(1) * t83 + g(2) * t77 + g(3) * t446 + t333;
t470 = -t19 * t303 - t273;
t469 = t54 * t303 + t273;
t468 = t257 * t370 - t362;
t347 = t253 * t370;
t467 = t347 - t364;
t173 = t319 * qJD(1);
t147 = t253 * t173;
t381 = qJ(4) * t371 + t147;
t394 = t254 * t257;
t395 = t253 * t258;
t430 = pkin(8) - pkin(9);
t466 = t430 * t364 + (-pkin(7) * t394 + pkin(9) * t395) * qJD(1) + t381;
t188 = t430 * t257;
t348 = -pkin(7) * t253 - pkin(3);
t279 = -pkin(9) * t390 + (-pkin(4) + t348) * t254;
t398 = t173 * t257;
t465 = -qJD(1) * t279 + qJD(3) * t188 + t398;
t434 = t91 ^ 2;
t464 = t303 ^ 2 - t434;
t463 = t54 * t91;
t133 = t167 * t258;
t99 = t167 * qJD(5) - t252 * t364 - t256 * t362;
t404 = -qJD(1) * t133 - t99;
t383 = -t252 * t468 + t253 * t360 + t256 * t467 - t257 * t361;
t461 = qJ(4) * t468 - t253 * qJD(4) - t230;
t38 = pkin(5) * t91 + qJ(6) * t303;
t396 = t253 * t256;
t130 = t252 * t394 - t254 * t396;
t304 = -t145 * t256 + t146 * t252;
t332 = -t256 * t10 + t252 * t13 + t52 * t360 + t42 * t361;
t444 = -t143 * t256 + t144 * t252;
t274 = g(1) * t304 + g(2) * t444 + g(3) * t130 - t332;
t196 = qJD(5) + t208;
t325 = qJD(3) + t370;
t299 = t325 * qJD(2);
t264 = (-t299 - t354) * t253 - t291;
t21 = qJD(5) * t91 - t252 * t84 + t256 * t264;
t451 = t196 * t91 - t21;
t460 = -pkin(5) * t304 + qJ(6) * t83;
t459 = t19 * t91 + qJDD(6);
t150 = -qJDD(5) + t161;
t432 = t196 ^ 2;
t458 = -t252 * t150 + t165 * t91 + t256 * t432;
t457 = -pkin(5) * t444 + qJ(6) * t77;
t418 = t91 * t303;
t352 = t431 * t253;
t407 = -qJD(3) * t352 + t347 * t431 - t461;
t20 = -t163 * t360 + t165 * t361 + t252 * t264 + t256 * t84;
t452 = -t196 * t303 + t20;
t450 = qJ(4) * t361 + t252 * t58 - t462 * t256 + t360 * t431;
t187 = t430 * t253;
t107 = t187 * t252 + t188 * t256;
t449 = -qJD(5) * t107 + t466 * t252 + t465 * t256;
t302 = t187 * t256 - t188 * t252;
t448 = -qJD(5) * t302 - t465 * t252 + t466 * t256;
t301 = t252 * t257 - t396;
t285 = t254 * t301;
t316 = g(1) * t259 + g(2) * t255;
t445 = t254 * t316;
t374 = t256 * qJ(4) - t252 * t431;
t237 = t253 * qJ(4);
t440 = t257 * pkin(3) + pkin(2) + t237;
t425 = pkin(8) * t161;
t438 = t208 * t72 + t425;
t419 = g(3) * t258;
t424 = pkin(8) * t208;
t437 = -qJD(3) * t424 + t419 - t445;
t373 = -t239 - t244;
t178 = -pkin(1) + t373;
t210 = pkin(7) * t395;
t243 = t258 * pkin(3);
t85 = pkin(4) * t258 + t210 + t243 + (-pkin(9) * t254 - t178) * t257;
t212 = pkin(7) * t390;
t376 = t253 * t178 + t212;
t108 = -qJ(4) * t258 + t376;
t397 = t253 * t254;
t95 = pkin(9) * t397 + t108;
t309 = t252 * t85 + t256 * t95;
t314 = -qJD(3) * t212 + t176 * t257 - t178 * t364;
t36 = pkin(9) * t340 + qJD(2) * t279 - t314;
t367 = qJD(2) * t254;
t380 = t253 * t176 + t178 * t362;
t351 = qJ(4) * t367 + t380;
t37 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t394 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t253) * t258 + t351;
t436 = -qJD(5) * t309 - t252 * t37 + t256 * t36;
t435 = -0.2e1 * pkin(1);
t433 = t165 ^ 2;
t262 = qJD(1) ^ 2;
t429 = pkin(1) * t262;
t428 = pkin(3) * t161;
t427 = pkin(5) * t150;
t426 = pkin(7) * t163;
t423 = g(1) * t255;
t420 = g(2) * t259;
t417 = pkin(5) * t383 - qJ(6) * t404 + qJD(6) * t301 + t407;
t416 = qJ(6) * t371 - t448;
t415 = -pkin(5) * t371 - t449;
t18 = t252 * t42 + t256 * t52;
t16 = qJ(6) * t196 + t18;
t413 = t16 * t196;
t412 = t18 * t196;
t65 = -t191 + t97;
t409 = t208 * t65;
t408 = t208 * t97;
t406 = -qJD(6) - t450;
t405 = t374 * qJD(5) + t252 * t462 + t256 * t58;
t403 = qJ(6) * t150;
t401 = t163 * t208;
t400 = t165 * t163;
t399 = t165 * t208;
t393 = t254 * t259;
t388 = t258 * t208;
t17 = -t252 * t52 + t256 * t42;
t384 = qJD(6) - t17;
t382 = -pkin(3) * t467 + t461;
t103 = t165 * pkin(3) + t163 * qJ(4);
t379 = (g(1) * t389 + g(2) * t392) * t254;
t342 = t258 * t357;
t377 = qJ(4) * t342 + qJD(4) * t394;
t248 = t254 ^ 2;
t372 = -t258 ^ 2 + t248;
t369 = qJD(2) * t165;
t366 = qJD(2) * t258;
t365 = qJD(3) * t163;
t359 = t163 * qJD(2);
t358 = t185 * qJD(3);
t227 = pkin(7) * t354;
t141 = -qJDD(2) * pkin(2) + pkin(7) * t335 + t227;
t350 = g(1) * t387 + g(2) * t391 + g(3) * t254;
t344 = t208 * t357;
t343 = t253 * t366;
t341 = t208 * t364;
t339 = t254 * t362;
t330 = -t143 * pkin(3) + qJ(4) * t144;
t329 = -t145 * pkin(3) + qJ(4) * t146;
t328 = t178 * t257 - t210;
t157 = t257 * pkin(4) + t440;
t324 = pkin(3) * t390 + qJ(4) * t395 - t373;
t323 = -pkin(7) - t352;
t322 = -g(1) * t444 + g(2) * t304;
t321 = g(1) * t77 - g(2) * t83;
t59 = -pkin(4) * t165 - t103;
t320 = t348 * t254;
t318 = -g(1) * t143 + g(2) * t145;
t317 = g(1) * t144 - g(2) * t146;
t315 = -t144 * pkin(3) + t259 * pkin(7) - qJ(4) * t143;
t313 = -pkin(5) * t130 + qJ(6) * t446;
t310 = -t252 * t95 + t256 * t85;
t308 = t358 - t425;
t305 = -qJ(4) * t252 - t256 * t431;
t297 = -g(1) * t145 - g(2) * t143 - g(3) * t397;
t24 = -t264 * pkin(3) + t84 * qJ(4) - t165 * qJD(4) + t141;
t290 = -t257 * t431 - t237;
t289 = -pkin(7) * qJDD(2) + t356 * t435;
t288 = t252 * t36 + t256 * t37 + t85 * t360 - t361 * t95;
t287 = t257 * t161 + t341;
t286 = -t186 * t364 + t353;
t261 = qJD(2) ^ 2;
t281 = pkin(7) * t261 + qJDD(1) * t435 - t423;
t204 = qJ(4) * t394;
t105 = t254 * t323 + t204;
t280 = t259 * pkin(1) + pkin(2) * t387 + t146 * pkin(3) + t255 * pkin(7) + pkin(8) * t393 + qJ(4) * t145;
t112 = t255 * t446;
t114 = t259 * t446;
t278 = g(1) * t114 + g(2) * t112 - g(3) * t133 - t150 * t302;
t111 = t255 * t285;
t113 = t259 * t285;
t132 = t252 * t390 - t256 * t395;
t277 = g(1) * t113 + g(2) * t111 - g(3) * t132 - t107 * t150;
t276 = -t297 - t326;
t271 = -t256 * t150 - t165 * t303 - t252 * t432;
t268 = t165 * t72 + qJDD(4) - t276;
t267 = -t274 + t459;
t266 = -t196 * t405 - t274;
t265 = g(1) * t146 + g(2) * t144 + g(3) * t394 - t208 * t96 - t286;
t43 = t290 * t363 + t323 * t366 + t377;
t14 = pkin(4) * t264 - t24;
t219 = g(2) * t393;
t215 = pkin(8) * t387;
t211 = pkin(8) * t391;
t172 = pkin(5) - t305;
t171 = -qJ(6) + t374;
t123 = -t204 + (pkin(3) * t253 + pkin(7)) * t254;
t109 = t243 - t328;
t102 = qJD(1) * t320 - t398;
t101 = -pkin(7) * t346 + t381;
t63 = pkin(3) * t208 + t385;
t61 = pkin(5) * t167 + qJ(6) * t301 + t157;
t55 = qJ(4) * t340 + pkin(7) * t366 + (t339 + t343) * pkin(3) - t377;
t53 = qJD(2) * t320 - t314;
t51 = -t84 - t401;
t50 = -qJD(4) * t258 + (-t254 * t357 - t258 * t364) * pkin(7) + t351;
t45 = qJD(2) * t133 + (qJD(3) - qJD(5)) * t285;
t44 = t252 * t342 + t254 * t99 - t256 * t343;
t41 = t105 - t313;
t32 = -pkin(5) * t258 - t310;
t31 = qJ(6) * t258 + t309;
t25 = t306 - t428;
t23 = t286 + t443;
t22 = -t38 + t59;
t15 = -pkin(5) * t196 + t384;
t6 = pkin(5) * t44 - qJ(6) * t45 - qJD(6) * t446 + t43;
t5 = pkin(5) * t367 - t436;
t4 = -qJ(6) * t367 + qJD(6) * t258 + t288;
t3 = t21 * pkin(5) + t20 * qJ(6) - t91 * qJD(6) + t14;
t2 = qJDD(6) + t332 + t427;
t1 = qJD(6) * t196 - t333 - t403;
t7 = [qJDD(1), -t420 + t423, t316, qJDD(1) * t248 + 0.2e1 * t254 * t335, 0.2e1 * t233 * t254 - 0.2e1 * t356 * t372, qJDD(2) * t254 + t258 * t261, qJDD(2) * t258 - t254 * t261, 0, t289 * t254 + (-t281 - t420) * t258, t254 * t281 + t258 * t289 + t219, -t84 * t394 + (-t340 + t342) * t165 (-t258 * t359 + (t232 + (-t165 - t346) * qJD(3)) * t254) * t257 + (-t165 * t366 + (-t257 * t299 - t198 + t365 + t84) * t254) * t253 (t84 - t344) * t258 + (t287 + t369) * t254, -t232 * t258 + (-t359 + (t208 + t370) * t362) * t254 + ((-t161 + t233) * t254 + (t208 + t325) * t366) * t253, -t161 * t258 - t208 * t367, -t314 * t208 + t328 * t161 + ((t185 * t253 + t426) * qJD(2) + t326) * t258 + (t257 * t358 + t96 * qJD(2) + t141 * t253 + (-t232 + (qJDD(1) * t253 + t334) * t254 + (-t208 + t325) * t368) * pkin(7)) * t254 + t317, t380 * t208 - t376 * t161 + (t185 * t357 + (-t341 + t369) * pkin(7) + t286) * t258 + (-t253 * t358 - t97 * qJD(2) + t141 * t257 + (-t84 - t344) * pkin(7)) * t254 + t318, -t109 * t161 - t123 * t232 + t55 * t163 + t53 * t208 + t25 * t258 + (-t63 * qJD(2) + (qJD(1) * t123 + t72) * t362) * t254 + ((qJDD(1) * t123 + t24) * t254 + (t123 * t325 + t72 * t258) * qJD(2)) * t253 + t317, t108 * t232 - t109 * t84 - t50 * t163 + t53 * t165 - t219 + (t63 * t390 + (-t108 * t325 - t65 * t258) * t253) * qJD(2) + (t423 + t25 * t257 + (-qJDD(1) * t108 - t23) * t253 + (-t63 * t253 + (-qJD(1) * t108 - t65) * t257) * qJD(3)) * t254, t108 * t161 + t123 * t84 - t165 * t55 - t208 * t50 + (-t357 * t72 - t23) * t258 + (qJD(2) * t65 - t24 * t257 + t364 * t72) * t254 - t318, -g(1) * t315 - g(2) * t280 + t23 * t108 + t25 * t109 + t24 * t123 - t300 * t423 + t65 * t50 + t63 * t53 + t72 * t55, -t20 * t446 + t45 * t91, t130 * t20 - t21 * t446 - t303 * t45 - t44 * t91, -t150 * t446 + t196 * t45 - t20 * t258 - t367 * t91, t130 * t150 - t196 * t44 - t21 * t258 + t303 * t367, -t150 * t258 - t196 * t367, t105 * t21 + t14 * t130 - t310 * t150 - t17 * t367 + t196 * t436 - t332 * t258 + t43 * t303 + t54 * t44 + t321, -t105 * t20 + t14 * t446 + t150 * t309 + t18 * t367 - t196 * t288 + t258 * t333 + t43 * t91 + t54 * t45 + t322, t130 * t3 + t15 * t367 + t150 * t32 + t19 * t44 - t196 * t5 - t2 * t258 + t21 * t41 + t303 * t6 + t321, -t1 * t130 + t15 * t45 - t16 * t44 + t2 * t446 - t20 * t32 - t21 * t31 - t254 * t423 - t303 * t4 + t5 * t91 + t219, t1 * t258 - t150 * t31 - t16 * t367 - t19 * t45 + t196 * t4 + t20 * t41 - t3 * t446 - t6 * t91 - t322, t1 * t31 + t16 * t4 + t3 * t41 + t19 * t6 + t2 * t32 + t15 * t5 - g(1) * (-pkin(4) * t144 - pkin(5) * t77 - qJ(6) * t444 + t315) - g(2) * (pkin(4) * t146 + pkin(5) * t83 - pkin(9) * t393 + qJ(6) * t304 + t280) - (-t254 * t430 + t349) * t423; 0, 0, 0, -t254 * t262 * t258, t372 * t262, t354, t233, qJDD(2), -t419 - t227 + (t316 + t429) * t254 (-pkin(7) * qJDD(1) + t429) * t258 + t350, -t84 * t253 - t257 * t399 (-t84 + t401) * t257 + (-t165 * qJD(3) + t232 - t294 * t253 + (-t339 + (t165 - t368) * t258) * qJD(1)) * t253, -t208 * t362 + t253 * t161 + (-t165 * t254 + t257 * t388) * qJD(1) (t163 * t254 - t253 * t388) * qJD(1) + t287, t208 * t371, pkin(2) * t232 + (-t96 * t254 - t258 * t426) * qJD(1) + (-t419 + t173 * t208 - t141 + (-pkin(2) * t371 + t424) * qJD(3)) * t257 + (-pkin(2) * t294 + (pkin(7) * t208 * t254 + (-t185 - t251) * t258) * qJD(1) + t308) * t253 + t379, pkin(2) * t84 - t147 * t208 + t308 * t257 + (-t185 * t390 + t97 * t254 + (-t165 * t258 + t208 * t394) * pkin(7)) * qJD(1) + (t141 + t437) * t253, t63 * t371 - t102 * t208 + t440 * t232 + t382 * t163 + (-t419 - t24 + (-t371 * t440 + t424) * qJD(3)) * t257 + (-t275 * t440 - t438) * t253 + t379, t101 * t163 - t102 * t165 + (t23 - t208 * t63 + (t232 + (t165 - t346) * qJD(3)) * pkin(8)) * t257 + (t25 + t409 + (-t257 * t275 + t365 - t84) * pkin(8)) * t253 - t350, -t65 * t371 + t101 * t208 - t440 * t84 - t382 * t165 + t438 * t257 + (-t24 - t437) * t253, -t65 * t101 - t63 * t102 - g(1) * t215 - g(2) * t211 - g(3) * t324 + t382 * t72 + (t23 * t257 + t25 * t253 + (-t253 * t65 + t257 * t63) * qJD(3)) * pkin(8) + (-t24 + t445) * t440, t20 * t301 + t404 * t91, t167 * t20 + t21 * t301 - t303 * t404 - t383 * t91, t150 * t301 + t196 * t404 + t371 * t91, t150 * t167 - t196 * t383 - t303 * t371, t196 * t371, t14 * t167 + t157 * t21 + t17 * t371 + t196 * t449 + t407 * t303 + t383 * t54 + t278, -t14 * t301 - t157 * t20 - t18 * t371 + t196 * t448 + t404 * t54 + t407 * t91 - t277, -t15 * t371 + t167 * t3 + t19 * t383 - t196 * t415 + t21 * t61 + t303 * t417 + t278, -t1 * t167 - t107 * t21 + t15 * t404 - t16 * t383 - t2 * t301 + t20 * t302 - t303 * t416 + t415 * t91 + t350, t16 * t371 - t19 * t404 + t196 * t416 + t20 * t61 + t3 * t301 - t417 * t91 + t277, t1 * t107 + t3 * t61 - t2 * t302 - g(1) * (-pkin(5) * t114 - pkin(9) * t387 - qJ(6) * t113 + t215) - g(2) * (-pkin(5) * t112 - pkin(9) * t391 - qJ(6) * t111 + t211) - g(3) * (pkin(4) * t390 + pkin(5) * t133 + qJ(6) * t132 + t324) + t417 * t19 + t416 * t16 + t415 * t15 + (g(3) * pkin(9) + t316 * (pkin(2) - t290)) * t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t400, -t163 ^ 2 + t433, t51, t264 - t399, t161, -t165 * t185 + t276 - t408, t163 * t185 + t265, -t103 * t163 - t268 - t408 + 0.2e1 * t428, pkin(3) * t84 + (t65 - t97) * t165 + (t63 - t385) * t163 + t264 * qJ(4), t103 * t165 - t163 * t72 + 0.2e1 * t149 - 0.2e1 * t189 - t265, t23 * qJ(4) - t25 * pkin(3) - t72 * t103 - t63 * t97 - g(1) * t329 - g(2) * t330 - g(3) * (-pkin(3) * t397 + t204) + t385 * t65, -t418, t464, t452, -t451, t150, -t150 * t305 - t303 * t59 + t266 + t463, t374 * t150 + t196 * t450 - t59 * t91 - t469, -t22 * t303 + (pkin(5) + t172) * t150 + t266 + t459, -t171 * t21 - t172 * t20 + (-t15 - t406) * t303 + (-t16 + t405) * t91, t22 * t91 + (qJ(6) - t171) * t150 + (-qJD(6) + t406) * t196 - t470, t1 * t171 + t2 * t172 - t19 * t22 - g(1) * (-pkin(4) * t145 + t329 - t460) - g(2) * (-pkin(4) * t143 + t330 - t457) - g(3) * (-t254 * t352 + t204 - t313) + t406 * t16 + t405 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161 + t400, t51, -t208 ^ 2 - t433, t268 + t409 - t428, 0, 0, 0, 0, 0, t271, -t458, t271, t252 * t451 + t256 * t452, t458, -t165 * t19 + (-t2 + t413) * t256 + (t15 * t196 + t1) * t252 + t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, -t464, -t452, t451, -t150, t274 + t412 - t463, t17 * t196 + t469, -t303 * t38 - t267 + t412 - 0.2e1 * t427, pkin(5) * t20 - qJ(6) * t21 + (t16 - t18) * t91 + (t15 - t384) * t303, -0.2e1 * t403 + t38 * t91 + (0.2e1 * qJD(6) - t17) * t196 + t470, -t2 * pkin(5) - g(1) * t460 - g(2) * t457 - g(3) * t313 + t1 * qJ(6) - t15 * t18 + t384 * t16 - t19 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150 + t418, -t452, -t432 - t434, t267 - t413 + t427;];
tau_reg  = t7;
