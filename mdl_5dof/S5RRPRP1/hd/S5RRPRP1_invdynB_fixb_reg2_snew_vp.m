% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRP1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:14
% EndTime: 2020-01-03 11:59:21
% DurationCPUTime: 5.24s
% Computational Cost: add. (18067->381), mult. (25982->523), div. (0->0), fcn. (14576->8), ass. (0->279)
t424 = qJD(1) + qJD(2);
t422 = t424 ^ 2;
t423 = qJDD(1) + qJDD(2);
t428 = sin(pkin(8));
t429 = cos(pkin(8));
t383 = t429 * t422 + t428 * t423;
t386 = t428 * t422 - t429 * t423;
t431 = sin(qJ(2));
t434 = cos(qJ(2));
t324 = t434 * t383 - t431 * t386;
t427 = g(1) - qJDD(3);
t362 = qJ(3) * t383 - t429 * t427;
t508 = qJ(3) * t386 - t428 * t427;
t272 = pkin(6) * t324 + t434 * t362 - t431 * t508;
t328 = t431 * t383 + t434 * t386;
t432 = sin(qJ(1));
t435 = cos(qJ(1));
t287 = t435 * t324 - t432 * t328;
t516 = pkin(6) * t328 + t431 * t362 + t434 * t508;
t525 = -pkin(5) * t287 - t435 * t272 + t432 * t516;
t507 = t432 * t324 + t435 * t328;
t524 = pkin(5) * t507 + t432 * t272 + t435 * t516;
t408 = t435 * g(2) + t432 * g(3);
t439 = qJDD(1) * pkin(1) - t408;
t407 = t432 * g(2) - t435 * g(3);
t498 = qJD(1) ^ 2;
t441 = -t498 * pkin(1) - t407;
t340 = t431 * t439 + t434 * t441;
t332 = -t422 * pkin(2) + t340;
t438 = -t431 * t441 + t434 * t439;
t437 = t423 * pkin(2) + t438;
t290 = t428 * t332 - t429 * t437;
t291 = t429 * t332 + t428 * t437;
t451 = t428 * t290 + t429 * t291;
t235 = t429 * t290 - t428 * t291;
t468 = t434 * t235;
t194 = -t431 * t451 + t468;
t474 = t431 * t235;
t510 = t434 * t451 + t474;
t521 = t432 * t194 + t435 * t510;
t175 = -t435 * t194 + t432 * t510;
t390 = t434 * t422 + t431 * t423;
t393 = t431 * t422 - t434 * t423;
t336 = t435 * t390 - t432 * t393;
t366 = pkin(6) * t390 - t434 * g(1);
t509 = pkin(6) * t393 - t431 * g(1);
t518 = -pkin(5) * t336 - t435 * t366 + t432 * t509;
t442 = t432 * t390 + t435 * t393;
t517 = pkin(5) * t442 + t432 * t366 + t435 * t509;
t450 = t434 * t340 - t431 * t438;
t295 = -t431 * t340 - t434 * t438;
t467 = t435 * t295;
t243 = t432 * t450 - t467;
t473 = t432 * t295;
t511 = t435 * t450 + t473;
t430 = sin(qJ(4));
t433 = cos(qJ(4));
t406 = t433 * t422 * t430;
t463 = qJDD(4) + t406;
t506 = t463 * pkin(4);
t280 = -t422 * pkin(3) + t423 * pkin(7) + t291;
t263 = t430 * t280 + t433 * t427;
t464 = qJD(4) * t433;
t457 = t424 * t464;
t475 = t430 * t423;
t378 = t457 + t475;
t369 = t378 * qJ(5);
t499 = -t369 - t263 + t506;
t436 = qJD(4) ^ 2;
t425 = t430 ^ 2;
t482 = t425 * t422;
t403 = -t436 - t482;
t399 = qJDD(4) - t406;
t476 = t430 * t399;
t346 = t433 * t403 - t476;
t497 = pkin(3) * t346;
t426 = t433 ^ 2;
t413 = t426 * t422;
t405 = -t413 - t436;
t477 = t430 * t463;
t349 = t433 * t405 - t477;
t465 = qJD(4) * t424;
t458 = t430 * t465;
t469 = t433 * t423;
t380 = -0.2e1 * t458 + t469;
t305 = t428 * t349 + t429 * t380;
t307 = t429 * t349 - t428 * t380;
t258 = t434 * t305 + t431 * t307;
t260 = -t431 * t305 + t434 * t307;
t212 = t432 * t258 - t435 * t260;
t496 = pkin(5) * t212;
t470 = t433 * t399;
t352 = -t430 * t403 - t470;
t377 = 0.2e1 * t457 + t475;
t306 = t428 * t352 - t429 * t377;
t308 = t429 * t352 + t428 * t377;
t259 = t434 * t306 + t431 * t308;
t261 = -t431 * t306 + t434 * t308;
t213 = t432 * t259 - t435 * t261;
t495 = pkin(5) * t213;
t466 = t425 + t426;
t388 = t466 * t423;
t394 = t413 + t482;
t330 = t428 * t388 + t429 * t394;
t331 = t429 * t388 - t428 * t394;
t288 = t434 * t330 + t431 * t331;
t289 = -t431 * t330 + t434 * t331;
t237 = t432 * t288 - t435 * t289;
t494 = pkin(5) * t237;
t493 = pkin(6) * t258;
t492 = pkin(6) * t259;
t491 = pkin(6) * t288;
t389 = t433 * t463;
t345 = t430 * t405 + t389;
t490 = pkin(7) * t345;
t489 = pkin(7) * t346;
t486 = qJ(3) * t305;
t485 = qJ(3) * t306;
t484 = qJ(3) * t330;
t483 = t424 * t430;
t240 = (qJ(5) * t464 - 0.2e1 * qJD(5) * t430) * t424 + t499;
t479 = t430 * t240;
t279 = -t423 * pkin(3) - t422 * pkin(7) + t290;
t478 = t430 * t279;
t472 = t433 * t240;
t471 = t433 * t279;
t264 = t433 * t280 - t430 * t427;
t462 = 0.2e1 * qJD(5) * t424;
t461 = t428 * t475;
t460 = t429 * t475;
t456 = -pkin(1) * t345 + pkin(6) * t260;
t455 = -pkin(1) * t346 + pkin(6) * t261;
t454 = -pkin(2) * t345 + qJ(3) * t307;
t453 = -pkin(2) * t346 + qJ(3) * t308;
t400 = -t432 * qJDD(1) - t435 * t498;
t452 = pkin(5) * t400 + t435 * g(1);
t218 = t430 * t263 + t433 * t264;
t357 = -t432 * t407 - t435 * t408;
t448 = t428 * t406;
t447 = t429 * t406;
t446 = pkin(1) * t258 + pkin(2) * t305 + pkin(3) * t380 + pkin(7) * t349;
t445 = pkin(1) * t259 + pkin(2) * t306 - pkin(3) * t377 + pkin(7) * t352;
t444 = pkin(1) * t288 + pkin(2) * t330 + pkin(3) * t394 + pkin(7) * t388;
t246 = -pkin(3) * t345 + t263;
t217 = t433 * t263 - t430 * t264;
t358 = t435 * t407 - t432 * t408;
t379 = -t458 + t469;
t397 = qJD(4) * pkin(4) - qJ(5) * t483;
t440 = t379 * qJ(5) - qJD(4) * t397 + t433 * t462 + t264;
t245 = -t379 * pkin(4) - qJ(5) * t413 + t397 * t483 + qJDD(5) + t279;
t410 = t430 * t462;
t404 = t413 - t436;
t402 = t436 - t482;
t401 = t435 * qJDD(1) - t432 * t498;
t395 = t413 - t482;
t375 = pkin(5) * t401 + t432 * g(1);
t374 = t466 * t465;
t356 = t428 * qJDD(4) + t429 * t374;
t355 = -t429 * qJDD(4) + t428 * t374;
t354 = t433 * t378 - t425 * t465;
t353 = -t430 * t379 - t426 * t465;
t351 = t430 * t404 + t470;
t350 = -t430 * t402 + t389;
t348 = t433 * t404 - t476;
t347 = t433 * t402 + t477;
t342 = (t379 - t458) * t433;
t341 = (t378 + t457) * t430;
t338 = -pkin(4) * t377 - qJ(5) * t399;
t322 = qJ(3) * t331;
t321 = -t430 * t377 + t433 * t380;
t320 = t433 * t377 + t430 * t380;
t316 = t429 * t350 + t461;
t315 = t429 * t348 + t428 * t469;
t314 = t428 * t350 - t460;
t313 = t428 * t348 - t429 * t469;
t312 = t429 * t354 - t448;
t311 = t429 * t353 + t448;
t310 = t428 * t354 + t447;
t309 = t428 * t353 - t447;
t300 = -t431 * t355 + t434 * t356;
t299 = t434 * t355 + t431 * t356;
t298 = t429 * t321 - t428 * t395;
t297 = t428 * t321 + t429 * t395;
t292 = pkin(1) * g(1) + pkin(6) * t450;
t282 = pkin(6) * t289;
t278 = -t431 * t314 + t434 * t316;
t277 = -t431 * t313 + t434 * t315;
t276 = t434 * t314 + t431 * t316;
t275 = t434 * t313 + t431 * t315;
t268 = -t431 * t310 + t434 * t312;
t267 = -t431 * t309 + t434 * t311;
t266 = t434 * t310 + t431 * t312;
t265 = t434 * t309 + t431 * t311;
t253 = t471 - t489;
t252 = t478 - t490;
t251 = t432 * t299 - t435 * t300;
t250 = t435 * t299 + t432 * t300;
t249 = -t431 * t297 + t434 * t298;
t248 = t434 * t297 + t431 * t298;
t247 = t264 - t497;
t242 = -qJ(5) * t403 + t245;
t241 = -pkin(4) * t413 + t440;
t239 = t410 + (-t457 + t475) * qJ(5) - t499;
t238 = pkin(4) * t380 + qJ(5) * t405 - t245;
t234 = t435 * t288 + t432 * t289;
t231 = pkin(5) * t234;
t230 = qJ(5) * t469 + (t394 - t413) * pkin(4) + t440;
t229 = pkin(2) * t427 + qJ(3) * t451;
t228 = -t497 + (-t403 - t413) * pkin(4) + t440;
t227 = t432 * t276 - t435 * t278;
t226 = t432 * t275 - t435 * t277;
t225 = t435 * t276 + t432 * t278;
t224 = t435 * t275 + t432 * t277;
t223 = -qJ(5) * t457 + t246 + t369 + t410 - 0.2e1 * t506;
t222 = t432 * t266 - t435 * t268;
t221 = t432 * t265 - t435 * t267;
t220 = t435 * t266 + t432 * t268;
t219 = t435 * t265 + t432 * t267;
t215 = -qJ(5) * t389 - t430 * t238 - t490;
t214 = t433 * t242 - t430 * t338 - t489;
t211 = t435 * t259 + t432 * t261;
t210 = t435 * t258 + t432 * t260;
t209 = pkin(5) * t211;
t208 = pkin(5) * t210;
t207 = t429 * t217 - t484;
t206 = t428 * t217 + t322;
t205 = t432 * t248 - t435 * t249;
t204 = t435 * t248 + t432 * t249;
t203 = -pkin(4) * t245 + qJ(5) * t241;
t202 = t433 * t241 - t479;
t201 = t430 * t241 + t472;
t200 = t429 * t218 + t428 * t279;
t199 = t428 * t218 - t429 * t279;
t198 = -t428 * t247 + t429 * t253 - t485;
t197 = -t428 * t246 + t429 * t252 - t486;
t196 = -t430 * t230 + t433 * t239;
t191 = t429 * t247 + t428 * t253 + t453;
t190 = t429 * t246 + t428 * t252 + t454;
t189 = -pkin(4) * t461 + t429 * t196 - t484;
t188 = pkin(4) * t460 + t428 * t196 + t322;
t187 = t429 * t202 + t428 * t245;
t186 = t428 * t202 - t429 * t245;
t185 = -pkin(3) * t201 - pkin(4) * t240;
t184 = t429 * t214 - t428 * t228 - t485;
t183 = t429 * t215 - t428 * t223 - t486;
t182 = t428 * t214 + t429 * t228 + t453;
t181 = t428 * t215 + t429 * t223 + t454;
t180 = -t431 * t206 + t434 * t207 - t491;
t179 = t434 * t206 + t431 * t207 + t282;
t178 = -t431 * t199 + t434 * t200;
t177 = t434 * t199 + t431 * t200;
t174 = pkin(6) * t194 + qJ(3) * t468 - t431 * t229;
t173 = pkin(1) * t427 + pkin(6) * t510 + qJ(3) * t474 + t434 * t229;
t172 = -qJ(3) * t199 - (pkin(3) * t428 - pkin(7) * t429) * t217;
t171 = -pkin(7) * t201 - qJ(5) * t472 - t430 * t203;
t170 = -t431 * t191 + t434 * t198 - t492;
t169 = -t431 * t190 + t434 * t197 - t493;
t168 = t434 * t191 + t431 * t198 + t455;
t167 = t434 * t190 + t431 * t197 + t456;
t166 = -t431 * t188 + t434 * t189 - t491;
t165 = t434 * t188 + t431 * t189 + t282;
t164 = -t431 * t186 + t434 * t187;
t163 = t434 * t186 + t431 * t187;
t162 = qJ(3) * t200 - (-pkin(3) * t429 - pkin(7) * t428 - pkin(2)) * t217;
t161 = -t431 * t182 + t434 * t184 - t492;
t160 = -t431 * t181 + t434 * t183 - t493;
t159 = t434 * t182 + t431 * t184 + t455;
t158 = t434 * t181 + t431 * t183 + t456;
t157 = t432 * t177 - t435 * t178;
t156 = t435 * t177 + t432 * t178;
t155 = -qJ(3) * t186 + t429 * t171 - t428 * t185;
t154 = t432 * t163 - t435 * t164;
t153 = t435 * t163 + t432 * t164;
t152 = -pkin(2) * t201 + qJ(3) * t187 + t428 * t171 + t429 * t185;
t151 = -pkin(6) * t177 - t431 * t162 + t434 * t172;
t150 = pkin(1) * t217 + pkin(6) * t178 + t434 * t162 + t431 * t172;
t149 = -pkin(6) * t163 - t431 * t152 + t434 * t155;
t148 = -pkin(1) * t201 + pkin(6) * t164 + t434 * t152 + t431 * t155;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t427, 0, 0, 0, 0, 0, 0, t345, t346, 0, -t217, 0, 0, 0, 0, 0, 0, t345, t346, 0, t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t401, t400, 0, t357, 0, 0, 0, 0, 0, 0, -t442, -t336, 0, t243, 0, 0, 0, 0, 0, 0, -t507, -t287, 0, t175, 0, 0, 0, 0, 0, 0, t210, t211, t234, t156, 0, 0, 0, 0, 0, 0, t210, t211, t234, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t400, t401, 0, t358, 0, 0, 0, 0, 0, 0, t336, -t442, 0, -t511, 0, 0, 0, 0, 0, 0, t287, -t507, 0, -t521, 0, 0, 0, 0, 0, 0, t212, t213, t237, t157, 0, 0, 0, 0, 0, 0, t212, t213, t237, t154; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t408, t407, 0, 0, 0, 0, 0, 0, 0, t423, -pkin(1) * t393 + t438, -pkin(1) * t390 - t340, 0, -pkin(1) * t295, 0, 0, 0, 0, 0, t423, -pkin(1) * t328 - pkin(2) * t386 - t290, -pkin(1) * t324 - pkin(2) * t383 - t291, 0, -pkin(1) * t194 - pkin(2) * t235, t341, t320, t347, t342, t351, 0, t446 - t471, t445 + t478, t218 + t444, pkin(1) * t177 + pkin(2) * t199 - pkin(3) * t279 + pkin(7) * t218, t341, t320, t347, t342, t351, 0, -qJ(5) * t477 + t433 * t238 + t446, t430 * t242 + t433 * t338 + t445, t433 * t230 + t430 * t239 + t444, pkin(1) * t163 + pkin(2) * t186 - pkin(3) * t245 + pkin(7) * t202 - qJ(5) * t479 + t433 * t203; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t400, 0, t401, 0, t452, -t375, -t358, -pkin(5) * t358, 0, 0, t336, 0, -t442, 0, t518, t517, t511, pkin(5) * t511 + pkin(6) * t473 + t435 * t292, 0, 0, t287, 0, -t507, 0, t525, t524, t521, pkin(5) * t521 + t435 * t173 + t432 * t174, t220, t204, t225, t219, t224, t250, t435 * t167 + t432 * t169 - t496, t435 * t168 + t432 * t170 - t495, t435 * t179 + t432 * t180 - t494, -pkin(5) * t157 + t435 * t150 + t432 * t151, t220, t204, t225, t219, t224, t250, t435 * t158 + t432 * t160 - t496, t435 * t159 + t432 * t161 - t495, t435 * t165 + t432 * t166 - t494, -pkin(5) * t154 + t435 * t148 + t432 * t149; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t401, 0, -t400, 0, t375, t452, t357, pkin(5) * t357, 0, 0, t442, 0, t336, 0, -t517, t518, t243, pkin(5) * t243 - pkin(6) * t467 + t432 * t292, 0, 0, t507, 0, t287, 0, -t524, t525, t175, pkin(5) * t175 + t432 * t173 - t435 * t174, t222, t205, t227, t221, t226, t251, t432 * t167 - t435 * t169 + t208, t432 * t168 - t435 * t170 + t209, t432 * t179 - t435 * t180 + t231, pkin(5) * t156 + t432 * t150 - t435 * t151, t222, t205, t227, t221, t226, t251, t432 * t158 - t435 * t160 + t208, t432 * t159 - t435 * t161 + t209, t432 * t165 - t435 * t166 + t231, pkin(5) * t153 + t432 * t148 - t435 * t149;];
tauB_reg = t1;
