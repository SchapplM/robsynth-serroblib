% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:24
% EndTime: 2022-01-20 10:34:33
% DurationCPUTime: 7.98s
% Computational Cost: add. (36284->388), mult. (48058->557), div. (0->0), fcn. (28562->10), ass. (0->257)
t425 = qJD(1) + qJD(2);
t420 = qJD(4) + t425;
t418 = t420 ^ 2;
t436 = cos(qJ(4));
t424 = qJDD(1) + qJDD(2);
t419 = qJDD(4) + t424;
t432 = sin(qJ(4));
t465 = t432 * t419;
t378 = t436 * t418 + t465;
t459 = t436 * t419;
t381 = t432 * t418 - t459;
t429 = sin(pkin(9));
t430 = cos(pkin(9));
t320 = t430 * t378 - t429 * t381;
t428 = g(3) - qJDD(3);
t358 = pkin(7) * t378 - t436 * t428;
t492 = pkin(7) * t381 - t432 * t428;
t270 = qJ(3) * t320 + t430 * t358 - t429 * t492;
t324 = t429 * t378 + t430 * t381;
t433 = sin(qJ(2));
t437 = cos(qJ(2));
t279 = t437 * t320 - t433 * t324;
t505 = qJ(3) * t324 + t429 * t358 + t430 * t492;
t214 = pkin(6) * t279 + t437 * t270 - t433 * t505;
t283 = t433 * t320 + t437 * t324;
t434 = sin(qJ(1));
t438 = cos(qJ(1));
t239 = t434 * t279 + t438 * t283;
t520 = pkin(6) * t283 + t433 * t270 + t437 * t505;
t529 = pkin(5) * t239 + t434 * t214 + t438 * t520;
t509 = t438 * t279 - t434 * t283;
t528 = pkin(5) * t509 + t438 * t214 - t434 * t520;
t410 = t434 * g(1) - t438 * g(2);
t400 = qJDD(1) * pkin(1) + t410;
t411 = t438 * g(1) + t434 * g(2);
t440 = qJD(1) ^ 2;
t401 = -t440 * pkin(1) - t411;
t345 = -t437 * t400 + t433 * t401;
t336 = t424 * pkin(2) - t345;
t346 = t433 * t400 + t437 * t401;
t423 = t425 ^ 2;
t337 = -t423 * pkin(2) + t346;
t297 = t429 * t336 + t430 * t337;
t290 = -t423 * pkin(3) + t297;
t443 = t430 * t336 - t429 * t337;
t441 = t424 * pkin(3) + t443;
t246 = t432 * t290 - t436 * t441;
t247 = t436 * t290 + t432 * t441;
t451 = t432 * t246 + t436 * t247;
t205 = t436 * t246 - t432 * t247;
t471 = t430 * t205;
t179 = -t429 * t451 + t471;
t473 = t429 * t205;
t494 = t430 * t451 + t473;
t166 = t437 * t179 - t433 * t494;
t515 = t433 * t179 + t437 * t494;
t152 = t434 * t166 + t438 * t515;
t525 = t438 * t166 - t434 * t515;
t388 = t430 * t423 + t429 * t424;
t391 = t429 * t423 - t430 * t424;
t330 = t437 * t388 - t433 * t391;
t364 = qJ(3) * t388 - t430 * t428;
t491 = qJ(3) * t391 - t429 * t428;
t278 = pkin(6) * t330 + t437 * t364 - t433 * t491;
t334 = t433 * t388 + t437 * t391;
t294 = t434 * t330 + t438 * t334;
t506 = pkin(6) * t334 + t433 * t364 + t437 * t491;
t522 = pkin(5) * t294 + t434 * t278 + t438 * t506;
t489 = t438 * t330 - t434 * t334;
t521 = pkin(5) * t489 + t438 * t278 - t434 * t506;
t450 = t430 * t297 - t429 * t443;
t251 = -t429 * t297 - t430 * t443;
t457 = t437 * t251;
t209 = -t433 * t450 + t457;
t463 = t433 * t251;
t495 = t437 * t450 + t463;
t184 = t434 * t209 + t438 * t495;
t514 = t438 * t209 - t434 * t495;
t393 = t437 * t423 + t433 * t424;
t396 = t433 * t423 - t437 * t424;
t340 = t434 * t393 + t438 * t396;
t369 = pkin(6) * t393 - t437 * g(3);
t493 = pkin(6) * t396 - t433 * g(3);
t508 = pkin(5) * t340 + t434 * t369 + t438 * t493;
t442 = t438 * t393 - t434 * t396;
t507 = pkin(5) * t442 + t438 * t369 - t434 * t493;
t449 = t433 * t345 + t437 * t346;
t301 = t437 * t345 - t433 * t346;
t456 = t438 * t301;
t496 = -t434 * t449 + t456;
t462 = t434 * t301;
t254 = t438 * t449 + t462;
t478 = pkin(1) * t428;
t477 = pkin(2) * t428;
t431 = sin(qJ(5));
t426 = t431 ^ 2;
t474 = t426 * t418;
t241 = -t419 * pkin(4) - t418 * pkin(8) + t246;
t469 = t431 * t241;
t435 = cos(qJ(5));
t406 = t435 * t418 * t431;
t397 = qJDD(5) + t406;
t468 = t431 * t397;
t398 = qJDD(5) - t406;
t467 = t431 * t398;
t466 = t431 * t419;
t461 = t435 * t241;
t460 = t435 * t398;
t412 = t435 * t419;
t242 = -t418 * pkin(4) + t419 * pkin(8) + t247;
t235 = t435 * t242 - t431 * t428;
t427 = t435 ^ 2;
t455 = t426 + t427;
t454 = qJD(5) * t420;
t453 = t431 * t454;
t452 = t435 * t454;
t234 = t431 * t242 + t435 * t428;
t199 = t431 * t234 + t435 * t235;
t360 = -t434 * t410 - t438 * t411;
t446 = t432 * t406;
t445 = t436 * t406;
t408 = t438 * qJDD(1) - t434 * t440;
t444 = -pkin(5) * t408 - t434 * g(3);
t198 = t435 * t234 - t431 * t235;
t359 = t438 * t410 - t434 * t411;
t439 = qJD(5) ^ 2;
t413 = t427 * t418;
t407 = t434 * qJDD(1) + t438 * t440;
t405 = -t413 - t439;
t404 = t413 - t439;
t403 = -t439 - t474;
t402 = t439 - t474;
t385 = t435 * t397;
t384 = -pkin(5) * t407 + t438 * g(3);
t383 = t413 - t474;
t382 = t413 + t474;
t376 = t455 * t419;
t374 = t412 - 0.2e1 * t453;
t373 = t412 - t453;
t372 = t452 + t466;
t371 = 0.2e1 * t452 + t466;
t370 = t455 * t454;
t354 = t432 * qJDD(5) + t436 * t370;
t353 = -t436 * qJDD(5) + t432 * t370;
t352 = -t431 * t403 - t460;
t351 = -t431 * t402 + t385;
t350 = t435 * t405 - t468;
t349 = t435 * t404 - t467;
t348 = t435 * t403 - t467;
t347 = t431 * t405 + t385;
t344 = t435 * t372 - t426 * t454;
t343 = -t431 * t373 - t427 * t454;
t327 = t436 * t376 - t432 * t382;
t326 = t432 * t376 + t436 * t382;
t319 = -t431 * t371 + t435 * t374;
t318 = t436 * t351 + t431 * t465;
t317 = t436 * t349 + t432 * t412;
t316 = t432 * t351 - t431 * t459;
t315 = t432 * t349 - t435 * t459;
t314 = t436 * t344 - t446;
t313 = t436 * t343 + t446;
t312 = t432 * t344 + t445;
t311 = t432 * t343 - t445;
t310 = t436 * t352 + t432 * t371;
t309 = t436 * t350 - t432 * t374;
t308 = t432 * t352 - t436 * t371;
t307 = t432 * t350 + t436 * t374;
t306 = -t429 * t353 + t430 * t354;
t305 = t430 * t353 + t429 * t354;
t304 = t436 * t319 - t432 * t383;
t303 = t432 * t319 + t436 * t383;
t298 = pkin(1) * g(3) + pkin(6) * t449;
t286 = -t429 * t326 + t430 * t327;
t285 = t430 * t326 + t429 * t327;
t274 = -t429 * t316 + t430 * t318;
t273 = -t429 * t315 + t430 * t317;
t272 = t430 * t316 + t429 * t318;
t271 = t430 * t315 + t429 * t317;
t266 = -t429 * t312 + t430 * t314;
t265 = -t429 * t311 + t430 * t313;
t264 = t430 * t312 + t429 * t314;
t263 = t430 * t311 + t429 * t313;
t262 = -t429 * t308 + t430 * t310;
t261 = -t429 * t307 + t430 * t309;
t260 = t430 * t308 + t429 * t310;
t259 = t430 * t307 + t429 * t309;
t258 = -t433 * t305 + t437 * t306;
t257 = t437 * t305 + t433 * t306;
t256 = -t429 * t303 + t430 * t304;
t255 = t430 * t303 + t429 * t304;
t248 = qJ(3) * t450 + t477;
t244 = -t433 * t285 + t437 * t286;
t243 = t437 * t285 + t433 * t286;
t232 = -t433 * t272 + t437 * t274;
t231 = -t433 * t271 + t437 * t273;
t230 = t437 * t272 + t433 * t274;
t229 = t437 * t271 + t433 * t273;
t228 = -t433 * t264 + t437 * t266;
t227 = -t433 * t263 + t437 * t265;
t226 = t437 * t264 + t433 * t266;
t225 = t437 * t263 + t433 * t265;
t224 = -pkin(8) * t348 + t461;
t223 = -pkin(8) * t347 + t469;
t222 = -pkin(4) * t348 + t235;
t221 = -pkin(4) * t347 + t234;
t220 = -t433 * t260 + t437 * t262;
t219 = -t433 * t259 + t437 * t261;
t218 = t437 * t260 + t433 * t262;
t217 = t437 * t259 + t433 * t261;
t216 = -t433 * t255 + t437 * t256;
t215 = t437 * t255 + t433 * t256;
t202 = pkin(3) * t428 + pkin(7) * t451;
t201 = -t434 * t243 + t438 * t244;
t200 = t438 * t243 + t434 * t244;
t196 = -pkin(7) * t326 + t436 * t198;
t195 = pkin(7) * t327 + t432 * t198;
t194 = -t434 * t218 + t438 * t220;
t193 = -t434 * t217 + t438 * t219;
t192 = t438 * t218 + t434 * t220;
t191 = t438 * t217 + t434 * t219;
t190 = -pkin(7) * t308 - t432 * t222 + t436 * t224;
t189 = -pkin(7) * t307 - t432 * t221 + t436 * t223;
t188 = -pkin(3) * t348 + pkin(7) * t310 + t436 * t222 + t432 * t224;
t187 = -pkin(3) * t347 + pkin(7) * t309 + t436 * t221 + t432 * t223;
t186 = t436 * t199 + t432 * t241;
t185 = t432 * t199 - t436 * t241;
t182 = pkin(6) * t209 + qJ(3) * t457 - t433 * t248;
t181 = pkin(6) * t495 + qJ(3) * t463 + t437 * t248 + t478;
t176 = -qJ(3) * t285 - t429 * t195 + t430 * t196;
t175 = qJ(3) * t286 + t430 * t195 + t429 * t196;
t174 = -qJ(3) * t260 - t429 * t188 + t430 * t190;
t173 = -qJ(3) * t259 - t429 * t187 + t430 * t189;
t172 = -t429 * t185 + t430 * t186;
t171 = t430 * t185 + t429 * t186;
t170 = -pkin(2) * t348 + qJ(3) * t262 + t430 * t188 + t429 * t190;
t169 = -pkin(2) * t347 + qJ(3) * t261 + t430 * t187 + t429 * t189;
t168 = -pkin(7) * t185 - (pkin(4) * t432 - pkin(8) * t436) * t198;
t163 = pkin(7) * t471 + qJ(3) * t179 - t429 * t202;
t162 = pkin(7) * t473 + qJ(3) * t494 + t430 * t202 + t477;
t161 = pkin(7) * t186 - (-pkin(4) * t436 - pkin(8) * t432 - pkin(3)) * t198;
t160 = -pkin(6) * t243 - t433 * t175 + t437 * t176;
t159 = pkin(6) * t244 + t437 * t175 + t433 * t176;
t158 = -t433 * t171 + t437 * t172;
t157 = t437 * t171 + t433 * t172;
t156 = -pkin(6) * t218 - t433 * t170 + t437 * t174;
t155 = -pkin(6) * t217 - t433 * t169 + t437 * t173;
t154 = -pkin(1) * t348 + pkin(6) * t220 + t437 * t170 + t433 * t174;
t153 = -pkin(1) * t347 + pkin(6) * t219 + t437 * t169 + t433 * t173;
t150 = pkin(6) * t166 - t433 * t162 + t437 * t163;
t149 = pkin(6) * t515 + t437 * t162 + t433 * t163 + t478;
t148 = -qJ(3) * t171 - t429 * t161 + t430 * t168;
t147 = -t434 * t157 + t438 * t158;
t146 = t438 * t157 + t434 * t158;
t145 = pkin(2) * t198 + qJ(3) * t172 + t430 * t161 + t429 * t168;
t144 = -pkin(6) * t157 - t433 * t145 + t437 * t148;
t143 = pkin(1) * t198 + pkin(6) * t158 + t437 * t145 + t433 * t148;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t407, -t408, 0, t360, 0, 0, 0, 0, 0, 0, -t442, t340, 0, t254, 0, 0, 0, 0, 0, 0, -t489, t294, 0, t184, 0, 0, 0, 0, 0, 0, -t509, t239, 0, t152, 0, 0, 0, 0, 0, 0, t193, t194, t201, t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t408, -t407, 0, t359, 0, 0, 0, 0, 0, 0, -t340, -t442, 0, -t496, 0, 0, 0, 0, 0, 0, -t294, -t489, 0, -t514, 0, 0, 0, 0, 0, 0, -t239, -t509, 0, -t525, 0, 0, 0, 0, 0, 0, t191, t192, t200, t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t428, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t428, 0, 0, 0, 0, 0, 0, t347, t348, 0, -t198; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t408, 0, -t407, 0, t444, -t384, -t359, -pkin(5) * t359, 0, 0, -t340, 0, -t442, 0, t508, t507, t496, pkin(5) * t496 + pkin(6) * t456 - t434 * t298, 0, 0, -t294, 0, -t489, 0, t522, t521, t514, pkin(5) * t514 - t434 * t181 + t438 * t182, 0, 0, -t239, 0, -t509, 0, t529, t528, t525, pkin(5) * t525 - t434 * t149 + t438 * t150, -t434 * t226 + t438 * t228, -t434 * t215 + t438 * t216, -t434 * t230 + t438 * t232, -t434 * t225 + t438 * t227, -t434 * t229 + t438 * t231, -t434 * t257 + t438 * t258, -pkin(5) * t191 - t434 * t153 + t438 * t155, -pkin(5) * t192 - t434 * t154 + t438 * t156, -pkin(5) * t200 - t434 * t159 + t438 * t160, -pkin(5) * t146 - t434 * t143 + t438 * t144; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t407, 0, t408, 0, t384, t444, t360, pkin(5) * t360, 0, 0, t442, 0, -t340, 0, -t507, t508, t254, pkin(5) * t254 + pkin(6) * t462 + t438 * t298, 0, 0, t489, 0, -t294, 0, -t521, t522, t184, pkin(5) * t184 + t438 * t181 + t434 * t182, 0, 0, t509, 0, -t239, 0, -t528, t529, t152, pkin(5) * t152 + t438 * t149 + t434 * t150, t438 * t226 + t434 * t228, t438 * t215 + t434 * t216, t438 * t230 + t434 * t232, t438 * t225 + t434 * t227, t438 * t229 + t434 * t231, t438 * t257 + t434 * t258, pkin(5) * t193 + t438 * t153 + t434 * t155, pkin(5) * t194 + t438 * t154 + t434 * t156, pkin(5) * t201 + t438 * t159 + t434 * t160, pkin(5) * t147 + t438 * t143 + t434 * t144; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t410, t411, 0, 0, 0, 0, 0, 0, 0, t424, -pkin(1) * t396 - t345, -pkin(1) * t393 - t346, 0, -pkin(1) * t301, 0, 0, 0, 0, 0, t424, -pkin(1) * t334 - pkin(2) * t391 + t443, -pkin(1) * t330 - pkin(2) * t388 - t297, 0, -pkin(1) * t209 - pkin(2) * t251, 0, 0, 0, 0, 0, t419, -pkin(1) * t283 - pkin(2) * t324 - pkin(3) * t381 - t246, -pkin(1) * t279 - pkin(2) * t320 - pkin(3) * t378 - t247, 0, -pkin(1) * t166 - pkin(2) * t179 - pkin(3) * t205, (t372 + t452) * t431, t435 * t371 + t431 * t374, t435 * t402 + t468, (t373 - t453) * t435, t431 * t404 + t460, 0, pkin(1) * t217 + pkin(2) * t259 + pkin(3) * t307 + pkin(4) * t374 + pkin(8) * t350 - t461, pkin(1) * t218 + pkin(2) * t260 + pkin(3) * t308 - pkin(4) * t371 + pkin(8) * t352 + t469, pkin(1) * t243 + pkin(2) * t285 + pkin(3) * t326 + pkin(4) * t382 + pkin(8) * t376 + t199, pkin(1) * t157 + pkin(2) * t171 + pkin(3) * t185 - pkin(4) * t241 + pkin(8) * t199;];
tauB_reg = t1;
