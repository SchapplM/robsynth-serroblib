% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRP2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:14
% EndTime: 2022-01-23 09:28:20
% DurationCPUTime: 5.29s
% Computational Cost: add. (16938->382), mult. (25982->522), div. (0->0), fcn. (14576->8), ass. (0->282)
t418 = qJD(1) + qJD(3);
t416 = t418 ^ 2;
t428 = cos(qJ(3));
t417 = qJDD(1) + qJDD(3);
t425 = sin(qJ(3));
t474 = t425 * t417;
t381 = t428 * t416 + t474;
t465 = t428 * t417;
t384 = t425 * t416 - t465;
t422 = sin(pkin(8));
t423 = cos(pkin(8));
t320 = t423 * t381 - t422 * t384;
t421 = g(3) - qJDD(2);
t358 = pkin(6) * t381 - t428 * t421;
t511 = pkin(6) * t384 - t425 * t421;
t268 = qJ(2) * t320 + t423 * t358 - t422 * t511;
t324 = t422 * t381 + t423 * t384;
t426 = sin(qJ(1));
t429 = cos(qJ(1));
t281 = t426 * t320 + t429 * t324;
t517 = qJ(2) * t324 + t422 * t358 + t423 * t511;
t524 = pkin(5) * t281 + t426 * t268 + t429 * t517;
t510 = t429 * t320 - t426 * t324;
t523 = pkin(5) * t510 + t429 * t268 - t426 * t517;
t403 = t429 * g(1) + t426 * g(2);
t431 = qJD(1) ^ 2;
t434 = -t431 * pkin(1) - t403;
t402 = t426 * g(1) - t429 * g(2);
t436 = qJDD(1) * pkin(1) + t402;
t331 = t422 * t436 + t423 * t434;
t329 = -t431 * pkin(2) + t331;
t433 = -t422 * t434 + t423 * t436;
t432 = qJDD(1) * pkin(2) + t433;
t287 = t425 * t329 - t428 * t432;
t288 = t428 * t329 + t425 * t432;
t448 = t425 * t287 + t428 * t288;
t233 = t428 * t287 - t425 * t288;
t480 = t423 * t233;
t190 = -t422 * t448 + t480;
t481 = t422 * t233;
t513 = t423 * t448 + t481;
t172 = t426 * t190 + t429 * t513;
t520 = t429 * t190 - t426 * t513;
t447 = t423 * t331 - t422 * t433;
t291 = -t422 * t331 - t423 * t433;
t463 = t429 * t291;
t514 = -t426 * t447 + t463;
t472 = t426 * t291;
t237 = t429 * t447 + t472;
t393 = t422 * qJDD(1) + t423 * t431;
t394 = t423 * qJDD(1) - t422 * t431;
t337 = -t426 * t393 + t429 * t394;
t362 = qJ(2) * t393 - t423 * t421;
t437 = -qJ(2) * t394 - t422 * t421;
t512 = -pkin(5) * t337 + t426 * t362 + t429 * t437;
t499 = t429 * t393 + t426 * t394;
t508 = pkin(5) * t499 + t429 * t362 - t426 * t437;
t424 = sin(qJ(4));
t427 = cos(qJ(4));
t401 = t427 * t416 * t424;
t459 = qJDD(4) + t401;
t507 = t459 * pkin(4);
t276 = -t416 * pkin(3) + t417 * pkin(7) + t288;
t259 = t424 * t276 + t427 * t421;
t460 = qJD(4) * t427;
t453 = t418 * t460;
t475 = t424 * t417;
t373 = t453 + t475;
t364 = t373 * qJ(5);
t498 = -t364 - t259 + t507;
t496 = pkin(1) * t421;
t430 = qJD(4) ^ 2;
t419 = t424 ^ 2;
t482 = t419 * t416;
t398 = -t430 - t482;
t392 = qJDD(4) - t401;
t476 = t424 * t392;
t344 = t427 * t398 - t476;
t495 = pkin(3) * t344;
t420 = t427 ^ 2;
t408 = t420 * t416;
t400 = -t408 - t430;
t477 = t424 * t459;
t346 = t427 * t400 - t477;
t461 = qJD(4) * t418;
t454 = t424 * t461;
t466 = t427 * t417;
t375 = -0.2e1 * t454 + t466;
t301 = t425 * t346 + t428 * t375;
t303 = t428 * t346 - t425 * t375;
t254 = t423 * t301 + t422 * t303;
t256 = -t422 * t301 + t423 * t303;
t206 = t429 * t254 + t426 * t256;
t494 = pkin(5) * t206;
t467 = t427 * t392;
t348 = -t424 * t398 - t467;
t372 = 0.2e1 * t453 + t475;
t302 = t425 * t348 - t428 * t372;
t304 = t428 * t348 + t425 * t372;
t255 = t423 * t302 + t422 * t304;
t257 = -t422 * t302 + t423 * t304;
t207 = t429 * t255 + t426 * t257;
t493 = pkin(5) * t207;
t462 = t419 + t420;
t379 = t462 * t417;
t385 = t408 + t482;
t326 = t425 * t379 + t428 * t385;
t327 = t428 * t379 - t425 * t385;
t283 = t423 * t326 + t422 * t327;
t284 = -t422 * t326 + t423 * t327;
t228 = t429 * t283 + t426 * t284;
t492 = pkin(5) * t228;
t491 = pkin(6) * t301;
t490 = pkin(6) * t302;
t489 = pkin(6) * t326;
t380 = t427 * t459;
t342 = t424 * t400 + t380;
t488 = pkin(7) * t342;
t487 = pkin(7) * t344;
t486 = qJ(2) * t254;
t485 = qJ(2) * t255;
t484 = qJ(2) * t283;
t483 = t418 * t424;
t238 = (qJ(5) * t460 - 0.2e1 * qJD(5) * t424) * t418 + t498;
t479 = t424 * t238;
t275 = -t417 * pkin(3) - t416 * pkin(7) + t287;
t478 = t424 * t275;
t469 = t427 * t238;
t468 = t427 * t275;
t260 = t427 * t276 - t424 * t421;
t458 = 0.2e1 * qJD(5) * t418;
t457 = t424 * t474;
t456 = t424 * t465;
t452 = -pkin(1) * t342 + qJ(2) * t256;
t451 = -pkin(1) * t344 + qJ(2) * t257;
t450 = -pkin(2) * t342 + pkin(6) * t303;
t449 = -pkin(2) * t344 + pkin(6) * t304;
t214 = t424 * t259 + t427 * t260;
t354 = -t426 * t402 - t429 * t403;
t444 = t425 * t401;
t443 = t428 * t401;
t442 = pkin(1) * t254 + pkin(2) * t301 + pkin(3) * t375 + pkin(7) * t346;
t441 = pkin(1) * t255 + pkin(2) * t302 - pkin(3) * t372 + pkin(7) * t348;
t440 = pkin(1) * t283 + pkin(2) * t326 + pkin(3) * t385 + pkin(7) * t379;
t244 = -pkin(3) * t342 + t259;
t396 = t429 * qJDD(1) - t426 * t431;
t439 = -pkin(5) * t396 - t426 * g(3);
t213 = t427 * t259 - t424 * t260;
t353 = t429 * t402 - t426 * t403;
t374 = -t454 + t466;
t390 = qJD(4) * pkin(4) - qJ(5) * t483;
t435 = t374 * qJ(5) - qJD(4) * t390 + t427 * t458 + t260;
t241 = -t374 * pkin(4) - qJ(5) * t408 + t390 * t483 + qJDD(5) + t275;
t405 = t424 * t458;
t399 = t408 - t430;
t397 = t430 - t482;
t395 = t426 * qJDD(1) + t429 * t431;
t386 = t408 - t482;
t370 = -pkin(5) * t395 + t429 * g(3);
t369 = t462 * t461;
t352 = t425 * qJDD(4) + t428 * t369;
t351 = -t428 * qJDD(4) + t425 * t369;
t350 = t427 * t373 - t419 * t461;
t349 = -t424 * t374 - t420 * t461;
t347 = -t424 * t397 + t380;
t345 = t427 * t399 - t476;
t343 = t427 * t397 + t477;
t341 = t424 * t399 + t467;
t334 = (t373 + t453) * t424;
t333 = (t374 - t454) * t427;
t332 = -pkin(4) * t372 - qJ(5) * t392;
t318 = pkin(6) * t327;
t314 = -t424 * t372 + t427 * t375;
t313 = t427 * t372 + t424 * t375;
t312 = t428 * t347 + t457;
t311 = t428 * t345 + t425 * t466;
t310 = t425 * t347 - t456;
t309 = t425 * t345 - t427 * t465;
t308 = t428 * t350 - t444;
t307 = t428 * t349 + t444;
t306 = t425 * t350 + t443;
t305 = t425 * t349 - t443;
t296 = -t422 * t351 + t423 * t352;
t295 = t423 * t351 + t422 * t352;
t294 = t428 * t314 - t425 * t386;
t293 = t425 * t314 + t428 * t386;
t286 = qJ(2) * t447 + t496;
t277 = qJ(2) * t284;
t272 = -t422 * t310 + t423 * t312;
t271 = -t422 * t309 + t423 * t311;
t270 = t423 * t310 + t422 * t312;
t269 = t423 * t309 + t422 * t311;
t264 = -t422 * t306 + t423 * t308;
t263 = -t422 * t305 + t423 * t307;
t262 = t423 * t306 + t422 * t308;
t261 = t423 * t305 + t422 * t307;
t249 = t468 - t487;
t248 = t478 - t488;
t247 = -t426 * t295 + t429 * t296;
t246 = t429 * t295 + t426 * t296;
t245 = t260 - t495;
t243 = -t422 * t293 + t423 * t294;
t242 = t423 * t293 + t422 * t294;
t240 = -qJ(5) * t398 + t241;
t239 = -pkin(4) * t408 + t435;
t235 = t405 + (-t453 + t475) * qJ(5) - t498;
t230 = pkin(4) * t375 + qJ(5) * t400 - t241;
t229 = -t426 * t283 + t429 * t284;
t227 = pkin(2) * t421 + pkin(6) * t448;
t226 = qJ(5) * t466 + (t385 - t408) * pkin(4) + t435;
t225 = pkin(5) * t229;
t224 = -t495 + (-t398 - t408) * pkin(4) + t435;
t223 = -qJ(5) * t453 + t244 + t364 + t405 - 0.2e1 * t507;
t222 = -t426 * t270 + t429 * t272;
t221 = -t426 * t269 + t429 * t271;
t220 = t429 * t270 + t426 * t272;
t219 = t429 * t269 + t426 * t271;
t218 = -t426 * t262 + t429 * t264;
t217 = -t426 * t261 + t429 * t263;
t216 = t429 * t262 + t426 * t264;
t215 = t429 * t261 + t426 * t263;
t211 = -qJ(5) * t380 - t424 * t230 - t488;
t210 = t427 * t240 - t424 * t332 - t487;
t209 = -t426 * t255 + t429 * t257;
t208 = -t426 * t254 + t429 * t256;
t205 = pkin(5) * t209;
t204 = pkin(5) * t208;
t203 = t428 * t213 - t489;
t202 = t425 * t213 + t318;
t201 = -t426 * t242 + t429 * t243;
t200 = t429 * t242 + t426 * t243;
t199 = -pkin(4) * t241 + qJ(5) * t239;
t198 = t428 * t214 + t425 * t275;
t197 = t425 * t214 - t428 * t275;
t196 = t427 * t239 - t479;
t195 = t424 * t239 + t469;
t194 = -t425 * t245 + t428 * t249 - t490;
t193 = -t425 * t244 + t428 * t248 - t491;
t192 = -t424 * t226 + t427 * t235;
t187 = t428 * t245 + t425 * t249 + t449;
t186 = t428 * t244 + t425 * t248 + t450;
t185 = -pkin(4) * t457 + t428 * t192 - t489;
t184 = pkin(4) * t456 + t425 * t192 + t318;
t183 = t428 * t196 + t425 * t241;
t182 = t425 * t196 - t428 * t241;
t181 = -pkin(3) * t195 - pkin(4) * t238;
t180 = t428 * t210 - t425 * t224 - t490;
t179 = t428 * t211 - t425 * t223 - t491;
t178 = t425 * t210 + t428 * t224 + t449;
t177 = t425 * t211 + t428 * t223 + t450;
t176 = -t422 * t202 + t423 * t203 - t484;
t175 = t423 * t202 + t422 * t203 + t277;
t174 = -t422 * t197 + t423 * t198;
t173 = t423 * t197 + t422 * t198;
t170 = pkin(6) * t480 + qJ(2) * t190 - t422 * t227;
t169 = pkin(6) * t481 + qJ(2) * t513 + t423 * t227 + t496;
t168 = -pkin(6) * t197 - (pkin(3) * t425 - pkin(7) * t428) * t213;
t167 = -pkin(7) * t195 - qJ(5) * t469 - t424 * t199;
t166 = -t422 * t187 + t423 * t194 - t485;
t165 = -t422 * t186 + t423 * t193 - t486;
t164 = t423 * t187 + t422 * t194 + t451;
t163 = t423 * t186 + t422 * t193 + t452;
t162 = -t422 * t184 + t423 * t185 - t484;
t161 = t423 * t184 + t422 * t185 + t277;
t160 = -t422 * t182 + t423 * t183;
t159 = t423 * t182 + t422 * t183;
t158 = pkin(6) * t198 - (-pkin(3) * t428 - pkin(7) * t425 - pkin(2)) * t213;
t157 = -t422 * t178 + t423 * t180 - t485;
t156 = -t422 * t177 + t423 * t179 - t486;
t155 = t423 * t178 + t422 * t180 + t451;
t154 = t423 * t177 + t422 * t179 + t452;
t153 = -t426 * t173 + t429 * t174;
t152 = t429 * t173 + t426 * t174;
t151 = -pkin(6) * t182 + t428 * t167 - t425 * t181;
t150 = -t426 * t159 + t429 * t160;
t149 = t429 * t159 + t426 * t160;
t148 = -pkin(2) * t195 + pkin(6) * t183 + t425 * t167 + t428 * t181;
t147 = -qJ(2) * t173 - t422 * t158 + t423 * t168;
t146 = pkin(1) * t213 + qJ(2) * t174 + t423 * t158 + t422 * t168;
t145 = -qJ(2) * t159 - t422 * t148 + t423 * t151;
t144 = -pkin(1) * t195 + qJ(2) * t160 + t423 * t148 + t422 * t151;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t395, -t396, 0, t354, 0, 0, 0, 0, 0, 0, -t499, -t337, 0, t237, 0, 0, 0, 0, 0, 0, -t510, t281, 0, t172, 0, 0, 0, 0, 0, 0, t208, t209, t229, t153, 0, 0, 0, 0, 0, 0, t208, t209, t229, t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t396, -t395, 0, t353, 0, 0, 0, 0, 0, 0, t337, -t499, 0, -t514, 0, 0, 0, 0, 0, 0, -t281, -t510, 0, -t520, 0, 0, 0, 0, 0, 0, t206, t207, t228, t152, 0, 0, 0, 0, 0, 0, t206, t207, t228, t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421, 0, 0, 0, 0, 0, 0, t342, t344, 0, -t213, 0, 0, 0, 0, 0, 0, t342, t344, 0, t195; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t396, 0, -t395, 0, t439, -t370, -t353, -pkin(5) * t353, 0, 0, t337, 0, -t499, 0, t512, t508, t514, pkin(5) * t514 + qJ(2) * t463 - t426 * t286, 0, 0, -t281, 0, -t510, 0, t524, t523, t520, pkin(5) * t520 - t426 * t169 + t429 * t170, t218, t201, t222, t217, t221, t247, -t426 * t163 + t429 * t165 - t494, -t426 * t164 + t429 * t166 - t493, -t426 * t175 + t429 * t176 - t492, -pkin(5) * t152 - t426 * t146 + t429 * t147, t218, t201, t222, t217, t221, t247, -t426 * t154 + t429 * t156 - t494, -t426 * t155 + t429 * t157 - t493, -t426 * t161 + t429 * t162 - t492, -pkin(5) * t149 - t426 * t144 + t429 * t145; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t395, 0, t396, 0, t370, t439, t354, pkin(5) * t354, 0, 0, t499, 0, t337, 0, -t508, t512, t237, pkin(5) * t237 + qJ(2) * t472 + t429 * t286, 0, 0, t510, 0, -t281, 0, -t523, t524, t172, pkin(5) * t172 + t429 * t169 + t426 * t170, t216, t200, t220, t215, t219, t246, t429 * t163 + t426 * t165 + t204, t429 * t164 + t426 * t166 + t205, t429 * t175 + t426 * t176 + t225, pkin(5) * t153 + t429 * t146 + t426 * t147, t216, t200, t220, t215, t219, t246, t429 * t154 + t426 * t156 + t204, t429 * t155 + t426 * t157 + t205, t429 * t161 + t426 * t162 + t225, pkin(5) * t150 + t429 * t144 + t426 * t145; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t402, t403, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t394 + t433, -pkin(1) * t393 - t331, 0, -pkin(1) * t291, 0, 0, 0, 0, 0, t417, -pkin(1) * t324 - pkin(2) * t384 - t287, -pkin(1) * t320 - pkin(2) * t381 - t288, 0, -pkin(1) * t190 - pkin(2) * t233, t334, t313, t343, t333, t341, 0, t442 - t468, t441 + t478, t214 + t440, pkin(1) * t173 + pkin(2) * t197 - pkin(3) * t275 + pkin(7) * t214, t334, t313, t343, t333, t341, 0, -qJ(5) * t477 + t427 * t230 + t442, t424 * t240 + t427 * t332 + t441, t427 * t226 + t424 * t235 + t440, pkin(1) * t159 + pkin(2) * t182 - pkin(3) * t241 + pkin(7) * t196 - qJ(5) * t479 + t427 * t199;];
tauB_reg = t1;
