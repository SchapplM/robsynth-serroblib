% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:26
% EndTime: 2020-01-03 11:52:36
% DurationCPUTime: 8.18s
% Computational Cost: add. (33385->385), mult. (48058->556), div. (0->0), fcn. (28562->10), ass. (0->257)
t414 = qJD(1) + qJD(3);
t407 = qJD(4) + t414;
t404 = t407 ^ 2;
t425 = cos(qJ(4));
t413 = qJDD(1) + qJDD(3);
t405 = qJDD(4) + t413;
t421 = sin(qJ(4));
t458 = t421 * t405;
t364 = t425 * t404 + t458;
t451 = t425 * t405;
t367 = t421 * t404 - t451;
t422 = sin(qJ(3));
t426 = cos(qJ(3));
t308 = t426 * t364 - t422 * t367;
t417 = g(1) - qJDD(2);
t345 = pkin(7) * t364 - t425 * t417;
t486 = pkin(7) * t367 - t421 * t417;
t257 = pkin(6) * t308 + t426 * t345 - t422 * t486;
t313 = t422 * t364 + t426 * t367;
t418 = sin(pkin(9));
t419 = cos(pkin(9));
t266 = t419 * t308 - t418 * t313;
t499 = pkin(6) * t313 + t422 * t345 + t426 * t486;
t201 = qJ(2) * t266 + t419 * t257 - t418 * t499;
t270 = t418 * t308 + t419 * t313;
t423 = sin(qJ(1));
t427 = cos(qJ(1));
t226 = t427 * t266 - t423 * t270;
t511 = qJ(2) * t270 + t418 * t257 + t419 * t499;
t520 = -pkin(5) * t226 - t427 * t201 + t423 * t511;
t500 = t423 * t266 + t427 * t270;
t519 = pkin(5) * t500 + t423 * t201 + t427 * t511;
t397 = t427 * g(2) + t423 * g(3);
t384 = qJDD(1) * pkin(1) - t397;
t396 = t423 * g(2) - t427 * g(3);
t468 = qJD(1) ^ 2;
t385 = -t468 * pkin(1) - t396;
t328 = -t419 * t384 + t418 * t385;
t324 = qJDD(1) * pkin(2) - t328;
t329 = t418 * t384 + t419 * t385;
t325 = -t468 * pkin(2) + t329;
t285 = t422 * t324 + t426 * t325;
t412 = t414 ^ 2;
t277 = -t412 * pkin(3) + t285;
t430 = t426 * t324 - t422 * t325;
t429 = t413 * pkin(3) + t430;
t233 = t421 * t277 - t425 * t429;
t234 = t425 * t277 + t421 * t429;
t439 = t421 * t233 + t425 * t234;
t192 = t425 * t233 - t421 * t234;
t449 = t426 * t192;
t166 = -t422 * t439 + t449;
t456 = t422 * t192;
t491 = t426 * t439 + t456;
t153 = t419 * t166 - t418 * t491;
t506 = t418 * t166 + t419 * t491;
t516 = t423 * t153 + t427 * t506;
t138 = -t427 * t153 + t423 * t506;
t375 = t426 * t412 + t422 * t413;
t378 = t422 * t412 - t426 * t413;
t317 = t419 * t375 - t418 * t378;
t351 = pkin(6) * t375 - t426 * t417;
t487 = pkin(6) * t378 - t422 * t417;
t265 = qJ(2) * t317 + t419 * t351 - t418 * t487;
t321 = t418 * t375 + t419 * t378;
t281 = t427 * t317 - t423 * t321;
t498 = qJ(2) * t321 + t418 * t351 + t419 * t487;
t513 = -pkin(5) * t281 - t427 * t265 + t423 * t498;
t484 = t423 * t317 + t427 * t321;
t512 = pkin(5) * t484 + t423 * t265 + t427 * t498;
t438 = t426 * t285 - t422 * t430;
t238 = -t422 * t285 - t426 * t430;
t463 = t419 * t238;
t196 = -t418 * t438 + t463;
t464 = t418 * t238;
t490 = t419 * t438 + t464;
t505 = t423 * t196 + t427 * t490;
t170 = -t427 * t196 + t423 * t490;
t437 = t418 * t328 + t419 * t329;
t288 = t419 * t328 - t418 * t329;
t447 = t427 * t288;
t240 = t423 * t437 - t447;
t454 = t423 * t288;
t489 = t427 * t437 + t454;
t386 = t418 * qJDD(1) + t419 * t468;
t387 = t419 * qJDD(1) - t418 * t468;
t337 = -t427 * t386 - t423 * t387;
t356 = qJ(2) * t386 - t419 * t417;
t431 = -qJ(2) * t387 - t418 * t417;
t488 = pkin(5) * t337 - t427 * t356 + t423 * t431;
t470 = t423 * t386 - t427 * t387;
t483 = pkin(5) * t470 + t423 * t356 + t427 * t431;
t467 = pkin(1) * t417;
t466 = pkin(2) * t417;
t420 = sin(qJ(5));
t415 = t420 ^ 2;
t465 = t415 * t404;
t230 = -t405 * pkin(4) - t404 * pkin(8) + t233;
t462 = t420 * t230;
t424 = cos(qJ(5));
t392 = t424 * t404 * t420;
t381 = qJDD(5) + t392;
t461 = t420 * t381;
t382 = qJDD(5) - t392;
t460 = t420 * t382;
t459 = t420 * t405;
t453 = t424 * t230;
t452 = t424 * t382;
t398 = t424 * t405;
t231 = -t404 * pkin(4) + t405 * pkin(8) + t234;
t222 = t424 * t231 - t420 * t417;
t416 = t424 ^ 2;
t444 = t415 + t416;
t443 = qJD(5) * t407;
t442 = t420 * t443;
t441 = t424 * t443;
t393 = -t423 * qJDD(1) - t427 * t468;
t440 = pkin(5) * t393 + t427 * g(1);
t221 = t420 * t231 + t424 * t417;
t186 = t420 * t221 + t424 * t222;
t346 = -t423 * t396 - t427 * t397;
t433 = t421 * t392;
t432 = t425 * t392;
t185 = t424 * t221 - t420 * t222;
t347 = t427 * t396 - t423 * t397;
t428 = qJD(5) ^ 2;
t399 = t416 * t404;
t394 = t427 * qJDD(1) - t423 * t468;
t391 = -t399 - t428;
t390 = t399 - t428;
t389 = -t428 - t465;
t388 = t428 - t465;
t371 = t424 * t381;
t370 = pkin(5) * t394 + t423 * g(1);
t369 = t399 - t465;
t368 = t399 + t465;
t363 = t444 * t405;
t361 = t398 - 0.2e1 * t442;
t360 = t398 - t442;
t359 = t441 + t459;
t358 = 0.2e1 * t441 + t459;
t357 = t444 * t443;
t341 = t421 * qJDD(5) + t425 * t357;
t340 = -t425 * qJDD(5) + t421 * t357;
t335 = -t420 * t389 - t452;
t334 = -t420 * t388 + t371;
t333 = t424 * t391 - t461;
t332 = t424 * t390 - t460;
t331 = t424 * t389 - t460;
t330 = t420 * t391 + t371;
t327 = t424 * t359 - t415 * t443;
t326 = -t420 * t360 - t416 * t443;
t311 = t425 * t363 - t421 * t368;
t307 = t421 * t363 + t425 * t368;
t306 = -t420 * t358 + t424 * t361;
t305 = t425 * t334 + t420 * t458;
t304 = t425 * t332 + t421 * t398;
t303 = t421 * t334 - t420 * t451;
t302 = t421 * t332 - t424 * t451;
t301 = t425 * t327 - t433;
t300 = t425 * t326 + t433;
t299 = t421 * t327 + t432;
t298 = t421 * t326 - t432;
t297 = t425 * t335 + t421 * t358;
t296 = t425 * t333 - t421 * t361;
t295 = t421 * t335 - t425 * t358;
t294 = t421 * t333 + t425 * t361;
t293 = -t422 * t340 + t426 * t341;
t292 = t426 * t340 + t422 * t341;
t291 = t425 * t306 - t421 * t369;
t290 = t421 * t306 + t425 * t369;
t283 = qJ(2) * t437 + t467;
t273 = -t422 * t307 + t426 * t311;
t272 = t426 * t307 + t422 * t311;
t261 = -t422 * t303 + t426 * t305;
t260 = -t422 * t302 + t426 * t304;
t259 = t426 * t303 + t422 * t305;
t258 = t426 * t302 + t422 * t304;
t253 = -t422 * t299 + t426 * t301;
t252 = -t422 * t298 + t426 * t300;
t251 = t426 * t299 + t422 * t301;
t250 = t426 * t298 + t422 * t300;
t249 = -t422 * t295 + t426 * t297;
t248 = -t422 * t294 + t426 * t296;
t247 = t426 * t295 + t422 * t297;
t246 = t426 * t294 + t422 * t296;
t245 = -t418 * t292 + t419 * t293;
t244 = t419 * t292 + t418 * t293;
t243 = -t422 * t290 + t426 * t291;
t242 = t426 * t290 + t422 * t291;
t235 = pkin(6) * t438 + t466;
t228 = -t418 * t272 + t419 * t273;
t227 = t419 * t272 + t418 * t273;
t219 = -t418 * t259 + t419 * t261;
t218 = -t418 * t258 + t419 * t260;
t217 = t419 * t259 + t418 * t261;
t216 = t419 * t258 + t418 * t260;
t215 = -t418 * t251 + t419 * t253;
t214 = -t418 * t250 + t419 * t252;
t213 = t419 * t251 + t418 * t253;
t212 = t419 * t250 + t418 * t252;
t211 = -pkin(8) * t331 + t453;
t210 = -pkin(8) * t330 + t462;
t209 = -pkin(4) * t331 + t222;
t208 = -pkin(4) * t330 + t221;
t207 = -t418 * t247 + t419 * t249;
t206 = -t418 * t246 + t419 * t248;
t205 = t419 * t247 + t418 * t249;
t204 = t419 * t246 + t418 * t248;
t203 = -t418 * t242 + t419 * t243;
t202 = t419 * t242 + t418 * t243;
t189 = pkin(3) * t417 + pkin(7) * t439;
t188 = t423 * t227 - t427 * t228;
t187 = t427 * t227 + t423 * t228;
t183 = -pkin(7) * t307 + t425 * t185;
t182 = pkin(7) * t311 + t421 * t185;
t181 = t423 * t205 - t427 * t207;
t180 = t423 * t204 - t427 * t206;
t179 = t427 * t205 + t423 * t207;
t178 = t427 * t204 + t423 * t206;
t177 = -pkin(7) * t295 - t421 * t209 + t425 * t211;
t176 = -pkin(7) * t294 - t421 * t208 + t425 * t210;
t175 = -pkin(3) * t331 + pkin(7) * t297 + t425 * t209 + t421 * t211;
t174 = -pkin(3) * t330 + pkin(7) * t296 + t425 * t208 + t421 * t210;
t173 = t425 * t186 + t421 * t230;
t172 = t421 * t186 - t425 * t230;
t169 = pkin(6) * t463 + qJ(2) * t196 - t418 * t235;
t168 = pkin(6) * t464 + qJ(2) * t490 + t419 * t235 + t467;
t163 = -pkin(6) * t272 - t422 * t182 + t426 * t183;
t162 = pkin(6) * t273 + t426 * t182 + t422 * t183;
t161 = -t422 * t172 + t426 * t173;
t160 = t426 * t172 + t422 * t173;
t159 = -pkin(6) * t247 - t422 * t175 + t426 * t177;
t158 = -pkin(6) * t246 - t422 * t174 + t426 * t176;
t157 = -pkin(2) * t331 + pkin(6) * t249 + t426 * t175 + t422 * t177;
t156 = -pkin(2) * t330 + pkin(6) * t248 + t426 * t174 + t422 * t176;
t155 = -pkin(7) * t172 - (pkin(4) * t421 - pkin(8) * t425) * t185;
t150 = pkin(6) * t166 + pkin(7) * t449 - t422 * t189;
t149 = pkin(6) * t491 + pkin(7) * t456 + t426 * t189 + t466;
t148 = pkin(7) * t173 - (-pkin(4) * t425 - pkin(8) * t421 - pkin(3)) * t185;
t147 = -qJ(2) * t227 - t418 * t162 + t419 * t163;
t146 = qJ(2) * t228 + t419 * t162 + t418 * t163;
t145 = -t418 * t160 + t419 * t161;
t144 = t419 * t160 + t418 * t161;
t143 = -qJ(2) * t205 - t418 * t157 + t419 * t159;
t142 = -qJ(2) * t204 - t418 * t156 + t419 * t158;
t141 = -pkin(1) * t331 + qJ(2) * t207 + t419 * t157 + t418 * t159;
t140 = -pkin(1) * t330 + qJ(2) * t206 + t419 * t156 + t418 * t158;
t137 = qJ(2) * t153 - t418 * t149 + t419 * t150;
t136 = qJ(2) * t506 + t419 * t149 + t418 * t150 + t467;
t135 = -pkin(6) * t160 - t422 * t148 + t426 * t155;
t134 = t423 * t144 - t427 * t145;
t133 = t427 * t144 + t423 * t145;
t132 = pkin(2) * t185 + pkin(6) * t161 + t426 * t148 + t422 * t155;
t131 = -qJ(2) * t144 - t418 * t132 + t419 * t135;
t130 = pkin(1) * t185 + qJ(2) * t145 + t419 * t132 + t418 * t135;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, 0, 0, 0, 0, 0, 0, t330, t331, 0, -t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t394, t393, 0, t346, 0, 0, 0, 0, 0, 0, -t470, t337, 0, t240, 0, 0, 0, 0, 0, 0, -t484, -t281, 0, t170, 0, 0, 0, 0, 0, 0, -t500, -t226, 0, t138, 0, 0, 0, 0, 0, 0, t178, t179, t187, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t393, t394, 0, t347, 0, 0, 0, 0, 0, 0, -t337, -t470, 0, -t489, 0, 0, 0, 0, 0, 0, t281, -t484, 0, -t505, 0, 0, 0, 0, 0, 0, t226, -t500, 0, -t516, 0, 0, 0, 0, 0, 0, t180, t181, t188, t134; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t397, t396, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t387 - t328, -pkin(1) * t386 - t329, 0, -pkin(1) * t288, 0, 0, 0, 0, 0, t413, -pkin(1) * t321 - pkin(2) * t378 + t430, -pkin(1) * t317 - pkin(2) * t375 - t285, 0, -pkin(1) * t196 - pkin(2) * t238, 0, 0, 0, 0, 0, t405, -pkin(1) * t270 - pkin(2) * t313 - pkin(3) * t367 - t233, -pkin(1) * t266 - pkin(2) * t308 - pkin(3) * t364 - t234, 0, -pkin(1) * t153 - pkin(2) * t166 - pkin(3) * t192, (t359 + t441) * t420, t424 * t358 + t420 * t361, t424 * t388 + t461, (t360 - t442) * t424, t420 * t390 + t452, 0, pkin(1) * t204 + pkin(2) * t246 + pkin(3) * t294 + pkin(4) * t361 + pkin(8) * t333 - t453, pkin(1) * t205 + pkin(2) * t247 + pkin(3) * t295 - pkin(4) * t358 + pkin(8) * t335 + t462, pkin(1) * t227 + pkin(2) * t272 + pkin(3) * t307 + pkin(4) * t368 + pkin(8) * t363 + t186, pkin(1) * t144 + pkin(2) * t160 + pkin(3) * t172 - pkin(4) * t230 + pkin(8) * t186; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t393, 0, t394, 0, t440, -t370, -t347, -pkin(5) * t347, 0, 0, -t337, 0, -t470, 0, t488, t483, t489, pkin(5) * t489 + qJ(2) * t454 + t427 * t283, 0, 0, t281, 0, -t484, 0, t513, t512, t505, pkin(5) * t505 + t427 * t168 + t423 * t169, 0, 0, t226, 0, -t500, 0, t520, t519, t516, pkin(5) * t516 + t427 * t136 + t423 * t137, t427 * t213 + t423 * t215, t427 * t202 + t423 * t203, t427 * t217 + t423 * t219, t427 * t212 + t423 * t214, t427 * t216 + t423 * t218, t427 * t244 + t423 * t245, -pkin(5) * t180 + t427 * t140 + t423 * t142, -pkin(5) * t181 + t427 * t141 + t423 * t143, -pkin(5) * t188 + t427 * t146 + t423 * t147, -pkin(5) * t134 + t427 * t130 + t423 * t131; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t394, 0, -t393, 0, t370, t440, t346, pkin(5) * t346, 0, 0, t470, 0, -t337, 0, -t483, t488, t240, pkin(5) * t240 - qJ(2) * t447 + t423 * t283, 0, 0, t484, 0, t281, 0, -t512, t513, t170, pkin(5) * t170 + t423 * t168 - t427 * t169, 0, 0, t500, 0, t226, 0, -t519, t520, t138, pkin(5) * t138 + t423 * t136 - t427 * t137, t423 * t213 - t427 * t215, t423 * t202 - t427 * t203, t423 * t217 - t427 * t219, t423 * t212 - t427 * t214, t423 * t216 - t427 * t218, t423 * t244 - t427 * t245, pkin(5) * t178 + t423 * t140 - t427 * t142, pkin(5) * t179 + t423 * t141 - t427 * t143, pkin(5) * t187 + t423 * t146 - t427 * t147, pkin(5) * t133 + t423 * t130 - t427 * t131;];
tauB_reg = t1;
