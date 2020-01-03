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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:00:38
% EndTime: 2020-01-03 12:00:47
% DurationCPUTime: 7.25s
% Computational Cost: add. (36284->383), mult. (48058->557), div. (0->0), fcn. (28562->10), ass. (0->257)
t426 = qJD(1) + qJD(2);
t419 = qJD(4) + t426;
t417 = t419 ^ 2;
t437 = cos(qJ(4));
t425 = qJDD(1) + qJDD(2);
t418 = qJDD(4) + t425;
t433 = sin(qJ(4));
t465 = t433 * t418;
t377 = t437 * t417 + t465;
t459 = t437 * t418;
t380 = t433 * t417 - t459;
t430 = sin(pkin(9));
t431 = cos(pkin(9));
t319 = t431 * t377 - t430 * t380;
t429 = g(1) - qJDD(3);
t357 = pkin(7) * t377 - t437 * t429;
t493 = pkin(7) * t380 - t433 * t429;
t269 = qJ(3) * t319 + t431 * t357 - t430 * t493;
t323 = t430 * t377 + t431 * t380;
t434 = sin(qJ(2));
t438 = cos(qJ(2));
t278 = t438 * t319 - t434 * t323;
t506 = qJ(3) * t323 + t430 * t357 + t431 * t493;
t213 = pkin(6) * t278 + t438 * t269 - t434 * t506;
t282 = t434 * t319 + t438 * t323;
t435 = sin(qJ(1));
t439 = cos(qJ(1));
t239 = t439 * t278 - t435 * t282;
t521 = pkin(6) * t282 + t434 * t269 + t438 * t506;
t530 = -pkin(5) * t239 - t439 * t213 + t435 * t521;
t510 = t435 * t278 + t439 * t282;
t529 = pkin(5) * t510 + t435 * t213 + t439 * t521;
t410 = t439 * g(2) + t435 * g(3);
t399 = qJDD(1) * pkin(1) - t410;
t409 = t435 * g(2) - t439 * g(3);
t479 = qJD(1) ^ 2;
t400 = -t479 * pkin(1) - t409;
t344 = -t438 * t399 + t434 * t400;
t335 = t425 * pkin(2) - t344;
t345 = t434 * t399 + t438 * t400;
t424 = t426 ^ 2;
t336 = -t424 * pkin(2) + t345;
t296 = t430 * t335 + t431 * t336;
t289 = -t424 * pkin(3) + t296;
t443 = t431 * t335 - t430 * t336;
t441 = t425 * pkin(3) + t443;
t245 = t433 * t289 - t437 * t441;
t246 = t437 * t289 + t433 * t441;
t450 = t433 * t245 + t437 * t246;
t204 = t437 * t245 - t433 * t246;
t471 = t431 * t204;
t178 = -t430 * t450 + t471;
t473 = t430 * t204;
t496 = t431 * t450 + t473;
t165 = t438 * t178 - t434 * t496;
t516 = t434 * t178 + t438 * t496;
t526 = t435 * t165 + t439 * t516;
t150 = -t439 * t165 + t435 * t516;
t387 = t431 * t424 + t430 * t425;
t390 = t430 * t424 - t431 * t425;
t329 = t438 * t387 - t434 * t390;
t363 = qJ(3) * t387 - t431 * t429;
t492 = qJ(3) * t390 - t430 * t429;
t277 = pkin(6) * t329 + t438 * t363 - t434 * t492;
t333 = t434 * t387 + t438 * t390;
t294 = t439 * t329 - t435 * t333;
t507 = pkin(6) * t333 + t434 * t363 + t438 * t492;
t523 = -pkin(5) * t294 - t439 * t277 + t435 * t507;
t490 = t435 * t329 + t439 * t333;
t522 = pkin(5) * t490 + t435 * t277 + t439 * t507;
t449 = t431 * t296 - t430 * t443;
t250 = -t430 * t296 - t431 * t443;
t457 = t438 * t250;
t208 = -t434 * t449 + t457;
t463 = t434 * t250;
t497 = t438 * t449 + t463;
t515 = t435 * t208 + t439 * t497;
t182 = -t439 * t208 + t435 * t497;
t392 = t438 * t424 + t434 * t425;
t395 = t434 * t424 - t438 * t425;
t340 = t439 * t392 - t435 * t395;
t368 = pkin(6) * t392 - t438 * g(1);
t494 = pkin(6) * t395 - t434 * g(1);
t509 = -pkin(5) * t340 - t439 * t368 + t435 * t494;
t442 = t435 * t392 + t439 * t395;
t508 = pkin(5) * t442 + t435 * t368 + t439 * t494;
t448 = t434 * t344 + t438 * t345;
t300 = t438 * t344 - t434 * t345;
t456 = t439 * t300;
t252 = t435 * t448 - t456;
t462 = t435 * t300;
t495 = t439 * t448 + t462;
t478 = pkin(1) * t429;
t477 = pkin(2) * t429;
t432 = sin(qJ(5));
t427 = t432 ^ 2;
t474 = t427 * t417;
t240 = -t418 * pkin(4) - t417 * pkin(8) + t245;
t469 = t432 * t240;
t436 = cos(qJ(5));
t405 = t436 * t417 * t432;
t396 = qJDD(5) + t405;
t468 = t432 * t396;
t397 = qJDD(5) - t405;
t467 = t432 * t397;
t466 = t432 * t418;
t461 = t436 * t240;
t460 = t436 * t397;
t411 = t436 * t418;
t241 = -t417 * pkin(4) + t418 * pkin(8) + t246;
t234 = t436 * t241 - t432 * t429;
t428 = t436 ^ 2;
t455 = t427 + t428;
t454 = qJD(5) * t419;
t453 = t432 * t454;
t452 = t436 * t454;
t406 = -t435 * qJDD(1) - t439 * t479;
t451 = pkin(5) * t406 + t439 * g(1);
t233 = t432 * t241 + t436 * t429;
t198 = t432 * t233 + t436 * t234;
t358 = -t435 * t409 - t439 * t410;
t445 = t433 * t405;
t444 = t437 * t405;
t197 = t436 * t233 - t432 * t234;
t359 = t439 * t409 - t435 * t410;
t440 = qJD(5) ^ 2;
t412 = t428 * t417;
t407 = t439 * qJDD(1) - t435 * t479;
t404 = -t412 - t440;
t403 = t412 - t440;
t402 = -t440 - t474;
t401 = t440 - t474;
t384 = t436 * t396;
t383 = pkin(5) * t407 + t435 * g(1);
t382 = t412 - t474;
t381 = t412 + t474;
t375 = t455 * t418;
t373 = t411 - 0.2e1 * t453;
t372 = t411 - t453;
t371 = t452 + t466;
t370 = 0.2e1 * t452 + t466;
t369 = t455 * t454;
t353 = t433 * qJDD(5) + t437 * t369;
t352 = -t437 * qJDD(5) + t433 * t369;
t351 = -t432 * t402 - t460;
t350 = -t432 * t401 + t384;
t349 = t436 * t404 - t468;
t348 = t436 * t403 - t467;
t347 = t436 * t402 - t467;
t346 = t432 * t404 + t384;
t343 = t436 * t371 - t427 * t454;
t342 = -t432 * t372 - t428 * t454;
t326 = t437 * t375 - t433 * t381;
t325 = t433 * t375 + t437 * t381;
t318 = -t432 * t370 + t436 * t373;
t317 = t437 * t350 + t432 * t465;
t316 = t437 * t348 + t433 * t411;
t315 = t433 * t350 - t432 * t459;
t314 = t433 * t348 - t436 * t459;
t313 = t437 * t343 - t445;
t312 = t437 * t342 + t445;
t311 = t433 * t343 + t444;
t310 = t433 * t342 - t444;
t309 = t437 * t351 + t433 * t370;
t308 = t437 * t349 - t433 * t373;
t307 = t433 * t351 - t437 * t370;
t306 = t433 * t349 + t437 * t373;
t305 = -t430 * t352 + t431 * t353;
t304 = t431 * t352 + t430 * t353;
t303 = t437 * t318 - t433 * t382;
t302 = t433 * t318 + t437 * t382;
t297 = pkin(1) * g(1) + pkin(6) * t448;
t285 = -t430 * t325 + t431 * t326;
t284 = t431 * t325 + t430 * t326;
t273 = -t430 * t315 + t431 * t317;
t272 = -t430 * t314 + t431 * t316;
t271 = t431 * t315 + t430 * t317;
t270 = t431 * t314 + t430 * t316;
t265 = -t430 * t311 + t431 * t313;
t264 = -t430 * t310 + t431 * t312;
t263 = t431 * t311 + t430 * t313;
t262 = t431 * t310 + t430 * t312;
t261 = -t430 * t307 + t431 * t309;
t260 = -t430 * t306 + t431 * t308;
t259 = t431 * t307 + t430 * t309;
t258 = t431 * t306 + t430 * t308;
t257 = -t434 * t304 + t438 * t305;
t256 = t438 * t304 + t434 * t305;
t255 = -t430 * t302 + t431 * t303;
t254 = t431 * t302 + t430 * t303;
t247 = qJ(3) * t449 + t477;
t243 = -t434 * t284 + t438 * t285;
t242 = t438 * t284 + t434 * t285;
t231 = -t434 * t271 + t438 * t273;
t230 = -t434 * t270 + t438 * t272;
t229 = t438 * t271 + t434 * t273;
t228 = t438 * t270 + t434 * t272;
t227 = -t434 * t263 + t438 * t265;
t226 = -t434 * t262 + t438 * t264;
t225 = t438 * t263 + t434 * t265;
t224 = t438 * t262 + t434 * t264;
t223 = -pkin(8) * t347 + t461;
t222 = -pkin(8) * t346 + t469;
t221 = -pkin(4) * t347 + t234;
t220 = -pkin(4) * t346 + t233;
t219 = -t434 * t259 + t438 * t261;
t218 = -t434 * t258 + t438 * t260;
t217 = t438 * t259 + t434 * t261;
t216 = t438 * t258 + t434 * t260;
t215 = -t434 * t254 + t438 * t255;
t214 = t438 * t254 + t434 * t255;
t201 = pkin(3) * t429 + pkin(7) * t450;
t200 = t435 * t242 - t439 * t243;
t199 = t439 * t242 + t435 * t243;
t195 = -pkin(7) * t325 + t437 * t197;
t194 = pkin(7) * t326 + t433 * t197;
t193 = t435 * t217 - t439 * t219;
t192 = t435 * t216 - t439 * t218;
t191 = t439 * t217 + t435 * t219;
t190 = t439 * t216 + t435 * t218;
t189 = -pkin(7) * t307 - t433 * t221 + t437 * t223;
t188 = -pkin(7) * t306 - t433 * t220 + t437 * t222;
t187 = -pkin(3) * t347 + pkin(7) * t309 + t437 * t221 + t433 * t223;
t186 = -pkin(3) * t346 + pkin(7) * t308 + t437 * t220 + t433 * t222;
t185 = t437 * t198 + t433 * t240;
t184 = t433 * t198 - t437 * t240;
t181 = pkin(6) * t208 + qJ(3) * t457 - t434 * t247;
t180 = pkin(6) * t497 + qJ(3) * t463 + t438 * t247 + t478;
t175 = -qJ(3) * t284 - t430 * t194 + t431 * t195;
t174 = qJ(3) * t285 + t431 * t194 + t430 * t195;
t173 = -qJ(3) * t259 - t430 * t187 + t431 * t189;
t172 = -qJ(3) * t258 - t430 * t186 + t431 * t188;
t171 = -t430 * t184 + t431 * t185;
t170 = t431 * t184 + t430 * t185;
t169 = -pkin(2) * t347 + qJ(3) * t261 + t431 * t187 + t430 * t189;
t168 = -pkin(2) * t346 + qJ(3) * t260 + t431 * t186 + t430 * t188;
t167 = -pkin(7) * t184 - (pkin(4) * t433 - pkin(8) * t437) * t197;
t162 = pkin(7) * t471 + qJ(3) * t178 - t430 * t201;
t161 = pkin(7) * t473 + qJ(3) * t496 + t431 * t201 + t477;
t160 = pkin(7) * t185 - (-pkin(4) * t437 - pkin(8) * t433 - pkin(3)) * t197;
t159 = -pkin(6) * t242 - t434 * t174 + t438 * t175;
t158 = pkin(6) * t243 + t438 * t174 + t434 * t175;
t157 = -t434 * t170 + t438 * t171;
t156 = t438 * t170 + t434 * t171;
t155 = -pkin(6) * t217 - t434 * t169 + t438 * t173;
t154 = -pkin(6) * t216 - t434 * t168 + t438 * t172;
t153 = -pkin(1) * t347 + pkin(6) * t219 + t438 * t169 + t434 * t173;
t152 = -pkin(1) * t346 + pkin(6) * t218 + t438 * t168 + t434 * t172;
t149 = pkin(6) * t165 - t434 * t161 + t438 * t162;
t148 = pkin(6) * t516 + t438 * t161 + t434 * t162 + t478;
t147 = -qJ(3) * t170 - t430 * t160 + t431 * t167;
t146 = t435 * t156 - t439 * t157;
t145 = t439 * t156 + t435 * t157;
t144 = pkin(2) * t197 + qJ(3) * t171 + t431 * t160 + t430 * t167;
t143 = -pkin(6) * t156 - t434 * t144 + t438 * t147;
t142 = pkin(1) * t197 + pkin(6) * t157 + t438 * t144 + t434 * t147;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t429, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t429, 0, 0, 0, 0, 0, 0, t346, t347, 0, -t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t407, t406, 0, t358, 0, 0, 0, 0, 0, 0, -t442, -t340, 0, t252, 0, 0, 0, 0, 0, 0, -t490, -t294, 0, t182, 0, 0, 0, 0, 0, 0, -t510, -t239, 0, t150, 0, 0, 0, 0, 0, 0, t190, t191, t199, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t406, t407, 0, t359, 0, 0, 0, 0, 0, 0, t340, -t442, 0, -t495, 0, 0, 0, 0, 0, 0, t294, -t490, 0, -t515, 0, 0, 0, 0, 0, 0, t239, -t510, 0, -t526, 0, 0, 0, 0, 0, 0, t192, t193, t200, t146; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t410, t409, 0, 0, 0, 0, 0, 0, 0, t425, -pkin(1) * t395 - t344, -pkin(1) * t392 - t345, 0, -pkin(1) * t300, 0, 0, 0, 0, 0, t425, -pkin(1) * t333 - pkin(2) * t390 + t443, -pkin(1) * t329 - pkin(2) * t387 - t296, 0, -pkin(1) * t208 - pkin(2) * t250, 0, 0, 0, 0, 0, t418, -pkin(1) * t282 - pkin(2) * t323 - pkin(3) * t380 - t245, -pkin(1) * t278 - pkin(2) * t319 - pkin(3) * t377 - t246, 0, -pkin(1) * t165 - pkin(2) * t178 - pkin(3) * t204, (t371 + t452) * t432, t436 * t370 + t432 * t373, t436 * t401 + t468, (t372 - t453) * t436, t432 * t403 + t460, 0, pkin(1) * t216 + pkin(2) * t258 + pkin(3) * t306 + pkin(4) * t373 + pkin(8) * t349 - t461, pkin(1) * t217 + pkin(2) * t259 + pkin(3) * t307 - pkin(4) * t370 + pkin(8) * t351 + t469, pkin(1) * t242 + pkin(2) * t284 + pkin(3) * t325 + pkin(4) * t381 + pkin(8) * t375 + t198, pkin(1) * t156 + pkin(2) * t170 + pkin(3) * t184 - pkin(4) * t240 + pkin(8) * t198; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t406, 0, t407, 0, t451, -t383, -t359, -pkin(5) * t359, 0, 0, t340, 0, -t442, 0, t509, t508, t495, pkin(5) * t495 + pkin(6) * t462 + t439 * t297, 0, 0, t294, 0, -t490, 0, t523, t522, t515, pkin(5) * t515 + t439 * t180 + t435 * t181, 0, 0, t239, 0, -t510, 0, t530, t529, t526, pkin(5) * t526 + t439 * t148 + t435 * t149, t439 * t225 + t435 * t227, t439 * t214 + t435 * t215, t439 * t229 + t435 * t231, t439 * t224 + t435 * t226, t439 * t228 + t435 * t230, t439 * t256 + t435 * t257, -pkin(5) * t192 + t439 * t152 + t435 * t154, -pkin(5) * t193 + t439 * t153 + t435 * t155, -pkin(5) * t200 + t439 * t158 + t435 * t159, -pkin(5) * t146 + t439 * t142 + t435 * t143; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t407, 0, -t406, 0, t383, t451, t358, pkin(5) * t358, 0, 0, t442, 0, t340, 0, -t508, t509, t252, pkin(5) * t252 - pkin(6) * t456 + t435 * t297, 0, 0, t490, 0, t294, 0, -t522, t523, t182, pkin(5) * t182 + t435 * t180 - t439 * t181, 0, 0, t510, 0, t239, 0, -t529, t530, t150, pkin(5) * t150 + t435 * t148 - t439 * t149, t435 * t225 - t439 * t227, t435 * t214 - t439 * t215, t435 * t229 - t439 * t231, t435 * t224 - t439 * t226, t435 * t228 - t439 * t230, t435 * t256 - t439 * t257, pkin(5) * t190 + t435 * t152 - t439 * t154, pkin(5) * t191 + t435 * t153 - t439 * t155, pkin(5) * t199 + t435 * t158 - t439 * t159, pkin(5) * t145 + t435 * t142 - t439 * t143;];
tauB_reg = t1;
