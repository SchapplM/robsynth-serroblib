% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRP5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:05
% EndTime: 2019-12-31 16:45:09
% DurationCPUTime: 4.15s
% Computational Cost: add. (6236->315), mult. (15550->426), div. (0->0), fcn. (10322->6), ass. (0->242)
t405 = sin(pkin(6));
t406 = cos(pkin(6));
t407 = sin(qJ(3));
t409 = cos(qJ(3));
t422 = t405 * t409 + t406 * t407;
t378 = t422 * qJD(1);
t373 = t378 ^ 2;
t411 = qJD(3) ^ 2;
t327 = t411 + t373;
t447 = t405 * t407;
t376 = (-t406 * t409 + t447) * qJD(1);
t342 = t378 * t376;
t478 = qJDD(3) + t342;
t461 = t478 * t407;
t274 = t327 * t409 + t461;
t460 = t478 * t409;
t295 = t327 * t407 - t460;
t234 = t274 * t406 - t295 * t405;
t541 = qJ(2) * t234;
t253 = t274 * t405 + t295 * t406;
t540 = qJ(2) * t253;
t408 = sin(qJ(1));
t539 = t253 * t408;
t410 = cos(qJ(1));
t538 = t253 * t410;
t537 = -pkin(1) * t234 - pkin(2) * t274;
t472 = t376 ^ 2;
t358 = t472 - t411;
t284 = t358 * t407 + t460;
t290 = t358 * t409 - t461;
t249 = t284 * t405 - t290 * t406;
t442 = qJDD(1) * t406;
t374 = qJDD(1) * t447 - t409 * t442;
t536 = t249 * t408 - t374 * t410;
t535 = t249 * t410 + t374 * t408;
t307 = -t472 - t373;
t480 = t422 * qJDD(1);
t491 = -t374 * t409 + t480 * t407;
t492 = -t374 * t407 - t409 * t480;
t498 = -t405 * t492 + t406 * t491;
t516 = t307 * t408 + t410 * t498;
t533 = pkin(4) * t516;
t445 = t378 * qJD(3);
t334 = t374 + 0.2e1 * t445;
t479 = qJDD(3) - t342;
t459 = t479 * t407;
t481 = -t472 - t411;
t487 = t409 * t481 - t459;
t321 = t409 * t479;
t493 = t407 * t481 + t321;
t500 = -t405 * t493 + t406 * t487;
t517 = t334 * t408 + t410 * t500;
t532 = pkin(4) * t517;
t518 = -t307 * t410 + t408 * t498;
t531 = pkin(4) * t518;
t519 = -t410 * t334 + t408 * t500;
t530 = pkin(4) * t519;
t529 = pkin(5) * t274;
t528 = pkin(5) * t295;
t359 = -t373 + t411;
t502 = t409 * t359 + t459;
t503 = -t359 * t407 + t321;
t514 = -t405 * t502 + t406 * t503;
t527 = t408 * t514 - t480 * t410;
t526 = t480 * t408 + t410 * t514;
t525 = t284 * t406 + t290 * t405;
t497 = t405 * t491 + t406 * t492;
t524 = qJ(2) * t497;
t499 = t405 * t487 + t406 * t493;
t523 = qJ(2) * t499;
t212 = -pkin(1) * t497 - pkin(2) * t492;
t522 = -pkin(1) * t499 - pkin(2) * t493;
t521 = -pkin(1) * t334 + qJ(2) * t500;
t520 = -pkin(1) * t307 + qJ(2) * t498;
t515 = t405 * t503 + t406 * t502;
t510 = pkin(5) * t487;
t509 = pkin(5) * t492;
t508 = pkin(5) * t493;
t339 = t373 - t472;
t505 = t339 * t408;
t504 = t339 * t410;
t390 = g(1) * t410 + g(2) * t408;
t412 = qJD(1) ^ 2;
t494 = -pkin(1) * t412 + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t390;
t501 = -pkin(2) * t307 + pkin(5) * t491;
t496 = 2 * qJD(4);
t368 = qJD(3) * t376;
t337 = t480 - t368;
t483 = -t368 + t337;
t495 = t483 * qJ(4);
t419 = (-t376 * t407 - t378 * t409) * qJD(3);
t356 = t407 * t445;
t435 = t409 * t368;
t424 = t356 - t435;
t475 = -t405 * t419 + t406 * t424;
t490 = t408 * qJDD(3) + t410 * t475;
t436 = t410 * t342;
t335 = -t374 - t445;
t421 = -t335 * t407 + t435;
t425 = t409 * t335 + t407 * t368;
t474 = -t405 * t425 + t406 * t421;
t489 = t408 * t474 + t436;
t488 = -qJDD(3) * t410 + t408 * t475;
t437 = t408 * t342;
t486 = t410 * t474 - t437;
t402 = t405 ^ 2;
t403 = t406 ^ 2;
t482 = t402 + t403;
t468 = g(3) * t406;
t471 = pkin(2) * t406;
t313 = -t468 + (-pkin(5) * qJDD(1) + t412 * t471 - t494) * t405;
t344 = -g(3) * t405 + t494 * t406;
t397 = t403 * t412;
t314 = -pkin(2) * t397 + pkin(5) * t442 + t344;
t264 = t407 * t313 + t409 * t314;
t324 = pkin(3) * t376 - qJ(4) * t378;
t423 = qJDD(3) * qJ(4) + qJD(3) * t496 - t376 * t324 + t264;
t477 = t412 * t482;
t476 = t405 * t424 + t406 * t419;
t473 = t405 * t421 + t406 * t425;
t470 = pkin(3) * t335;
t469 = pkin(3) * t409;
t467 = qJDD(1) * pkin(1);
t263 = -t409 * t313 + t407 * t314;
t224 = -t263 * t409 + t264 * t407;
t466 = t224 * t405;
t465 = t224 * t406;
t464 = t483 * t407;
t389 = g(1) * t408 - t410 * g(2);
t426 = -qJDD(2) + t389;
t328 = (pkin(1) + t471) * qJDD(1) + (t482 * pkin(5) + qJ(2)) * t412 + t426;
t463 = t328 * t407;
t462 = t328 * t409;
t458 = t334 * t407;
t456 = t334 * t409;
t369 = qJ(2) * t412 + t426 + t467;
t454 = t369 * t408;
t453 = t369 * t410;
t449 = t402 * t412;
t448 = t405 * t406;
t446 = -t307 - t411;
t441 = qJDD(1) * t408;
t440 = qJDD(1) * t410;
t434 = -qJ(4) * t407 - pkin(2);
t433 = t369 + t467;
t225 = t263 * t407 + t409 * t264;
t343 = t494 * t405 + t468;
t283 = t343 * t405 + t406 * t344;
t352 = -t389 * t408 - t410 * t390;
t430 = t378 * t324 + qJDD(4) + t263;
t388 = -t408 * t412 + t440;
t429 = -pkin(4) * t388 - g(3) * t408;
t300 = t337 * t407 + t409 * t445;
t301 = t337 * t409 - t356;
t259 = -t300 * t405 + t301 * t406;
t428 = t408 * t259 - t436;
t427 = t410 * t259 + t437;
t282 = t343 * t406 - t344 * t405;
t351 = t389 * t410 - t390 * t408;
t387 = t410 * t412 + t441;
t381 = t406 * t477;
t348 = -t381 * t408 + t406 * t440;
t420 = t381 * t410 + t406 * t441;
t418 = -qJDD(3) * pkin(3) + t430;
t416 = -pkin(3) * t445 + t378 * t496 + t328;
t415 = t416 + t495;
t396 = t403 * qJDD(1);
t395 = t402 * qJDD(1);
t386 = t397 - t449;
t385 = t397 + t449;
t384 = t396 - t395;
t383 = t396 + t395;
t380 = t405 * t477;
t370 = -pkin(4) * t387 + g(3) * t410;
t355 = t388 * t448;
t354 = t387 * t448;
t349 = t380 * t410 + t405 * t441;
t347 = t380 * t408 - t405 * t440;
t346 = t383 * t410 - t385 * t408;
t345 = t383 * t408 + t385 * t410;
t336 = t480 - 0.2e1 * t368;
t280 = -t336 * t407 - t456;
t278 = t336 * t409 - t458;
t268 = t283 * t410 - t454;
t267 = t283 * t408 + t453;
t266 = t456 + t464;
t265 = -t409 * t483 + t458;
t261 = -t462 + t529;
t260 = -t463 - t508;
t256 = t300 * t406 + t301 * t405;
t245 = -pkin(2) * t336 - t463 + t528;
t242 = -pkin(2) * t334 + t462 + t510;
t241 = qJ(4) * t411 - t418;
t240 = -pkin(3) * t411 + t423;
t239 = -t278 * t405 + t280 * t406;
t232 = t415 + t470;
t231 = pkin(3) * t480 + qJ(4) * t374 + t212;
t230 = t336 * t408 + t538;
t228 = -t336 * t410 + t539;
t226 = -t265 * t405 + t266 * t406;
t223 = (-t334 + t335) * pkin(3) + t415;
t220 = t446 * qJ(4) + t418;
t219 = t446 * pkin(3) + t423;
t216 = t416 + t470 + 0.2e1 * t495;
t215 = -t408 * t483 - t538;
t214 = t410 * t483 - t539;
t213 = pkin(2) * t328 + pkin(5) * t225;
t211 = -t224 - t509;
t210 = t264 - t537;
t209 = -qJ(4) * t456 - t223 * t407 - t508;
t208 = t240 * t409 - t241 * t407;
t207 = t240 * t407 + t241 * t409;
t206 = t225 + t501;
t205 = t263 + t522;
t204 = t223 * t409 + t434 * t334 + t510;
t203 = -pkin(3) * t464 + t216 * t409 - t529;
t202 = -t245 * t405 + t261 * t406 + t541;
t201 = t225 * t406 - t466;
t200 = t225 * t405 + t465;
t199 = (-t481 - t411) * qJ(4) + (-qJDD(3) - t479) * pkin(3) + t430 + t522;
t198 = -t528 + t216 * t407 + (pkin(2) + t469) * t483;
t197 = -t242 * t405 + t260 * t406 - t523;
t196 = t201 * t410 - t328 * t408;
t195 = t201 * t408 + t328 * t410;
t194 = -t219 * t407 + t220 * t409 - t509;
t193 = -qJ(4) * t478 + (-t327 + t411) * pkin(3) - t423 + t537;
t192 = t219 * t409 + t220 * t407 + t501;
t191 = -pkin(1) * t200 - pkin(2) * t224;
t190 = -t207 * t405 + t208 * t406;
t189 = t207 * t406 + t208 * t405;
t188 = -pkin(5) * t207 + (-pkin(3) * t407 + qJ(4) * t409) * t232;
t187 = -t206 * t405 + t211 * t406 - t524;
t186 = -t204 * t405 + t209 * t406 - t523;
t185 = t190 * t410 - t232 * t408;
t184 = t190 * t408 + t232 * t410;
t183 = pkin(5) * t208 + (-t434 + t469) * t232;
t182 = -pkin(5) * t465 - qJ(2) * t200 - t213 * t405;
t181 = -t198 * t405 + t203 * t406 - t541;
t180 = -t192 * t405 + t194 * t406 - t524;
t179 = -pkin(1) * t189 - pkin(2) * t207 - pkin(3) * t241 - qJ(4) * t240;
t178 = -qJ(2) * t189 - t183 * t405 + t188 * t406;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t387, -t388, 0, t352, 0, 0, 0, 0, 0, 0, -t420, t349, t346, t268, 0, 0, 0, 0, 0, 0, t517, t230, t516, t196, 0, 0, 0, 0, 0, 0, t517, t516, t215, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t388, -t387, 0, t351, 0, 0, 0, 0, 0, 0, t348, t347, t345, t267, 0, 0, 0, 0, 0, 0, t519, t228, t518, t195, 0, 0, 0, 0, 0, 0, t519, t518, t214, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t282, 0, 0, 0, 0, 0, 0, t499, -t234, t497, t200, 0, 0, 0, 0, 0, 0, t499, t497, t234, t189; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t388, 0, -t387, 0, t429, -t370, -t351, -pkin(4) * t351, t355, t384 * t410 - t386 * t408, t349, -t355, t420, 0, -pkin(4) * t348 - t343 * t408 - t405 * t453, -pkin(4) * t347 - t344 * t408 - t406 * t453, -pkin(4) * t345 + t282 * t410, -pkin(4) * t267 - (pkin(1) * t408 - qJ(2) * t410) * t282, t427, t239 * t410 + t505, t526, t486, -t535, t490, t197 * t410 - t205 * t408 - t530, -pkin(4) * t228 + t202 * t410 - t210 * t408, t187 * t410 - t212 * t408 - t531, -pkin(4) * t195 + t182 * t410 - t191 * t408, t427, t526, t226 * t410 - t505, t490, t535, t486, t186 * t410 - t199 * t408 - t530, t180 * t410 - t231 * t408 - t531, -pkin(4) * t214 + t181 * t410 - t193 * t408, -pkin(4) * t184 + t178 * t410 - t179 * t408; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t387, 0, t388, 0, t370, t429, t352, pkin(4) * t352, t354, t384 * t408 + t386 * t410, t347, -t354, -t348, 0, -pkin(4) * t420 + t343 * t410 - t405 * t454, pkin(4) * t349 + t344 * t410 - t406 * t454, pkin(4) * t346 + t282 * t408, pkin(4) * t268 - (-pkin(1) * t410 - qJ(2) * t408) * t282, t428, t239 * t408 - t504, t527, t489, -t536, t488, t197 * t408 + t205 * t410 + t532, pkin(4) * t230 + t202 * t408 + t210 * t410, t187 * t408 + t212 * t410 + t533, pkin(4) * t196 + t182 * t408 + t191 * t410, t428, t527, t226 * t408 + t504, t488, t536, t489, t186 * t408 + t199 * t410 + t532, t180 * t408 + t231 * t410 + t533, pkin(4) * t215 + t181 * t408 + t193 * t410, pkin(4) * t185 + t178 * t408 + t179 * t410; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t389, t390, 0, 0, t395, 0.2e1 * t405 * t442, 0, t396, 0, 0, -qJ(2) * t381 + t433 * t406, qJ(2) * t380 - t433 * t405, pkin(1) * t385 + qJ(2) * t383 + t283, pkin(1) * t369 + qJ(2) * t283, t256, t278 * t406 + t280 * t405, t515, t473, t525, t476, t242 * t406 + t260 * t405 + t521, -pkin(1) * t336 + t245 * t406 + t261 * t405 + t540, t206 * t406 + t211 * t405 + t520, pkin(1) * t328 - pkin(5) * t466 + qJ(2) * t201 + t213 * t406, t256, t515, t265 * t406 + t266 * t405, t476, -t525, t473, t204 * t406 + t209 * t405 + t521, t192 * t406 + t194 * t405 + t520, pkin(1) * t483 + t198 * t406 + t203 * t405 - t540, pkin(1) * t232 + qJ(2) * t190 + t183 * t406 + t188 * t405;];
tauB_reg = t1;
