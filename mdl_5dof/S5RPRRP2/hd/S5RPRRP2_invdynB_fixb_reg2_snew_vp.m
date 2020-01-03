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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:45:32
% EndTime: 2020-01-03 11:45:42
% DurationCPUTime: 6.10s
% Computational Cost: add. (16938->383), mult. (25982->522), div. (0->0), fcn. (14576->8), ass. (0->282)
t419 = qJD(1) + qJD(3);
t417 = t419 ^ 2;
t429 = cos(qJ(3));
t418 = qJDD(1) + qJDD(3);
t426 = sin(qJ(3));
t474 = t426 * t418;
t380 = t429 * t417 + t474;
t467 = t429 * t418;
t383 = t426 * t417 - t467;
t423 = sin(pkin(8));
t424 = cos(pkin(8));
t319 = t424 * t380 - t423 * t383;
t422 = g(1) - qJDD(2);
t357 = pkin(6) * t380 - t429 * t422;
t512 = pkin(6) * t383 - t426 * t422;
t267 = qJ(2) * t319 + t424 * t357 - t423 * t512;
t323 = t423 * t380 + t424 * t383;
t427 = sin(qJ(1));
t430 = cos(qJ(1));
t281 = t430 * t319 - t427 * t323;
t518 = qJ(2) * t323 + t423 * t357 + t424 * t512;
t525 = -pkin(5) * t281 - t430 * t267 + t427 * t518;
t511 = t427 * t319 + t430 * t323;
t524 = pkin(5) * t511 + t427 * t267 + t430 * t518;
t402 = t430 * g(2) + t427 * g(3);
t434 = qJDD(1) * pkin(1) - t402;
t401 = t427 * g(2) - t430 * g(3);
t497 = qJD(1) ^ 2;
t436 = -t497 * pkin(1) - t401;
t330 = t423 * t434 + t424 * t436;
t328 = -t497 * pkin(2) + t330;
t433 = -t423 * t436 + t424 * t434;
t432 = qJDD(1) * pkin(2) + t433;
t286 = t426 * t328 - t429 * t432;
t287 = t429 * t328 + t426 * t432;
t447 = t426 * t286 + t429 * t287;
t232 = t429 * t286 - t426 * t287;
t480 = t424 * t232;
t189 = -t423 * t447 + t480;
t481 = t423 * t232;
t514 = t424 * t447 + t481;
t521 = t427 * t189 + t430 * t514;
t170 = -t430 * t189 + t427 * t514;
t446 = t424 * t330 - t423 * t433;
t290 = -t423 * t330 - t424 * t433;
t465 = t430 * t290;
t235 = t427 * t446 - t465;
t472 = t427 * t290;
t515 = t430 * t446 + t472;
t392 = t423 * qJDD(1) + t424 * t497;
t393 = t424 * qJDD(1) - t423 * t497;
t337 = -t430 * t392 - t427 * t393;
t361 = qJ(2) * t392 - t424 * t422;
t437 = -qJ(2) * t393 - t423 * t422;
t513 = pkin(5) * t337 - t430 * t361 + t427 * t437;
t500 = t427 * t392 - t430 * t393;
t510 = pkin(5) * t500 + t427 * t361 + t430 * t437;
t425 = sin(qJ(4));
t428 = cos(qJ(4));
t400 = t428 * t417 * t425;
t459 = qJDD(4) + t400;
t508 = t459 * pkin(4);
t275 = -t417 * pkin(3) + t418 * pkin(7) + t287;
t258 = t425 * t275 + t428 * t422;
t460 = qJD(4) * t428;
t453 = t419 * t460;
t475 = t425 * t418;
t372 = t453 + t475;
t363 = t372 * qJ(5);
t499 = -t363 - t258 + t508;
t496 = pkin(1) * t422;
t431 = qJD(4) ^ 2;
t420 = t425 ^ 2;
t482 = t420 * t417;
t397 = -t431 - t482;
t391 = qJDD(4) - t400;
t476 = t425 * t391;
t341 = t428 * t397 - t476;
t495 = pkin(3) * t341;
t421 = t428 ^ 2;
t407 = t421 * t417;
t399 = -t407 - t431;
t477 = t425 * t459;
t344 = t428 * t399 - t477;
t461 = qJD(4) * t419;
t454 = t425 * t461;
t468 = t428 * t418;
t374 = -0.2e1 * t454 + t468;
t300 = t426 * t344 + t429 * t374;
t302 = t429 * t344 - t426 * t374;
t253 = t424 * t300 + t423 * t302;
t255 = -t423 * t300 + t424 * t302;
t207 = t427 * t253 - t430 * t255;
t494 = pkin(5) * t207;
t469 = t428 * t391;
t347 = -t425 * t397 - t469;
t371 = 0.2e1 * t453 + t475;
t301 = t426 * t347 - t429 * t371;
t303 = t429 * t347 + t426 * t371;
t254 = t424 * t301 + t423 * t303;
t256 = -t423 * t301 + t424 * t303;
t208 = t427 * t254 - t430 * t256;
t493 = pkin(5) * t208;
t462 = t420 + t421;
t378 = t462 * t418;
t384 = t407 + t482;
t325 = t426 * t378 + t429 * t384;
t326 = t429 * t378 - t426 * t384;
t282 = t424 * t325 + t423 * t326;
t283 = -t423 * t325 + t424 * t326;
t228 = t427 * t282 - t430 * t283;
t492 = pkin(5) * t228;
t491 = pkin(6) * t300;
t490 = pkin(6) * t301;
t489 = pkin(6) * t325;
t379 = t428 * t459;
t340 = t425 * t399 + t379;
t488 = pkin(7) * t340;
t487 = pkin(7) * t341;
t486 = qJ(2) * t253;
t485 = qJ(2) * t254;
t484 = qJ(2) * t282;
t483 = t419 * t425;
t237 = (qJ(5) * t460 - 0.2e1 * qJD(5) * t425) * t419 + t499;
t479 = t425 * t237;
t274 = -t418 * pkin(3) - t417 * pkin(7) + t286;
t478 = t425 * t274;
t471 = t428 * t237;
t470 = t428 * t274;
t259 = t428 * t275 - t425 * t422;
t458 = 0.2e1 * qJD(5) * t419;
t457 = t425 * t474;
t456 = t425 * t467;
t452 = -pkin(1) * t340 + qJ(2) * t255;
t451 = -pkin(1) * t341 + qJ(2) * t256;
t450 = -pkin(2) * t340 + pkin(6) * t302;
t449 = -pkin(2) * t341 + pkin(6) * t303;
t394 = -t427 * qJDD(1) - t430 * t497;
t448 = pkin(5) * t394 + t430 * g(1);
t213 = t425 * t258 + t428 * t259;
t352 = -t427 * t401 - t430 * t402;
t443 = t426 * t400;
t442 = t429 * t400;
t441 = pkin(1) * t253 + pkin(2) * t300 + pkin(3) * t374 + pkin(7) * t344;
t440 = pkin(1) * t254 + pkin(2) * t301 - pkin(3) * t371 + pkin(7) * t347;
t439 = pkin(1) * t282 + pkin(2) * t325 + pkin(3) * t384 + pkin(7) * t378;
t243 = -pkin(3) * t340 + t258;
t212 = t428 * t258 - t425 * t259;
t353 = t430 * t401 - t427 * t402;
t373 = -t454 + t468;
t389 = qJD(4) * pkin(4) - qJ(5) * t483;
t435 = t373 * qJ(5) - qJD(4) * t389 + t428 * t458 + t259;
t240 = -t373 * pkin(4) - qJ(5) * t407 + t389 * t483 + qJDD(5) + t274;
t404 = t425 * t458;
t398 = t407 - t431;
t396 = t431 - t482;
t395 = t430 * qJDD(1) - t427 * t497;
t385 = t407 - t482;
t369 = pkin(5) * t395 + t427 * g(1);
t368 = t462 * t461;
t351 = t426 * qJDD(4) + t429 * t368;
t350 = -t429 * qJDD(4) + t426 * t368;
t349 = t428 * t372 - t420 * t461;
t348 = -t425 * t373 - t421 * t461;
t346 = t425 * t398 + t469;
t345 = -t425 * t396 + t379;
t343 = t428 * t398 - t476;
t342 = t428 * t396 + t477;
t333 = (t373 - t454) * t428;
t332 = (t372 + t453) * t425;
t331 = -pkin(4) * t371 - qJ(5) * t391;
t317 = pkin(6) * t326;
t313 = -t425 * t371 + t428 * t374;
t312 = t428 * t371 + t425 * t374;
t311 = t429 * t345 + t457;
t310 = t429 * t343 + t426 * t468;
t309 = t426 * t345 - t456;
t308 = t426 * t343 - t428 * t467;
t307 = t429 * t349 - t443;
t306 = t429 * t348 + t443;
t305 = t426 * t349 + t442;
t304 = t426 * t348 - t442;
t295 = -t423 * t350 + t424 * t351;
t294 = t424 * t350 + t423 * t351;
t293 = t429 * t313 - t426 * t385;
t292 = t426 * t313 + t429 * t385;
t285 = qJ(2) * t446 + t496;
t276 = qJ(2) * t283;
t271 = -t423 * t309 + t424 * t311;
t270 = -t423 * t308 + t424 * t310;
t269 = t424 * t309 + t423 * t311;
t268 = t424 * t308 + t423 * t310;
t263 = -t423 * t305 + t424 * t307;
t262 = -t423 * t304 + t424 * t306;
t261 = t424 * t305 + t423 * t307;
t260 = t424 * t304 + t423 * t306;
t248 = t470 - t487;
t247 = t478 - t488;
t246 = t427 * t294 - t430 * t295;
t245 = t430 * t294 + t427 * t295;
t244 = t259 - t495;
t242 = -t423 * t292 + t424 * t293;
t241 = t424 * t292 + t423 * t293;
t239 = -qJ(5) * t397 + t240;
t238 = -pkin(4) * t407 + t435;
t234 = t404 + (-t453 + t475) * qJ(5) - t499;
t229 = pkin(4) * t374 + qJ(5) * t399 - t240;
t227 = t430 * t282 + t427 * t283;
t226 = pkin(2) * t422 + pkin(6) * t447;
t225 = qJ(5) * t468 + (t384 - t407) * pkin(4) + t435;
t224 = pkin(5) * t227;
t223 = -t495 + (-t397 - t407) * pkin(4) + t435;
t222 = -qJ(5) * t453 + t243 + t363 + t404 - 0.2e1 * t508;
t221 = t427 * t269 - t430 * t271;
t220 = t427 * t268 - t430 * t270;
t219 = t430 * t269 + t427 * t271;
t218 = t430 * t268 + t427 * t270;
t217 = t427 * t261 - t430 * t263;
t216 = t427 * t260 - t430 * t262;
t215 = t430 * t261 + t427 * t263;
t214 = t430 * t260 + t427 * t262;
t210 = -qJ(5) * t379 - t425 * t229 - t488;
t209 = t428 * t239 - t425 * t331 - t487;
t206 = t430 * t254 + t427 * t256;
t205 = t430 * t253 + t427 * t255;
t204 = pkin(5) * t206;
t203 = pkin(5) * t205;
t202 = t429 * t212 - t489;
t201 = t426 * t212 + t317;
t200 = t427 * t241 - t430 * t242;
t199 = t430 * t241 + t427 * t242;
t198 = -pkin(4) * t240 + qJ(5) * t238;
t197 = t429 * t213 + t426 * t274;
t196 = t426 * t213 - t429 * t274;
t195 = t428 * t238 - t479;
t194 = t425 * t238 + t471;
t193 = -t426 * t244 + t429 * t248 - t490;
t192 = -t426 * t243 + t429 * t247 - t491;
t191 = -t425 * t225 + t428 * t234;
t186 = t429 * t244 + t426 * t248 + t449;
t185 = t429 * t243 + t426 * t247 + t450;
t184 = -pkin(4) * t457 + t429 * t191 - t489;
t183 = pkin(4) * t456 + t426 * t191 + t317;
t182 = t429 * t195 + t426 * t240;
t181 = t426 * t195 - t429 * t240;
t180 = -pkin(3) * t194 - pkin(4) * t237;
t179 = t429 * t209 - t426 * t223 - t490;
t178 = t429 * t210 - t426 * t222 - t491;
t177 = t426 * t209 + t429 * t223 + t449;
t176 = t426 * t210 + t429 * t222 + t450;
t175 = -t423 * t201 + t424 * t202 - t484;
t174 = t424 * t201 + t423 * t202 + t276;
t173 = -t423 * t196 + t424 * t197;
t172 = t424 * t196 + t423 * t197;
t169 = pkin(6) * t480 + qJ(2) * t189 - t423 * t226;
t168 = pkin(6) * t481 + qJ(2) * t514 + t424 * t226 + t496;
t167 = -pkin(6) * t196 - (pkin(3) * t426 - pkin(7) * t429) * t212;
t166 = -pkin(7) * t194 - qJ(5) * t471 - t425 * t198;
t165 = -t423 * t186 + t424 * t193 - t485;
t164 = -t423 * t185 + t424 * t192 - t486;
t163 = t424 * t186 + t423 * t193 + t451;
t162 = t424 * t185 + t423 * t192 + t452;
t161 = -t423 * t183 + t424 * t184 - t484;
t160 = t424 * t183 + t423 * t184 + t276;
t159 = -t423 * t181 + t424 * t182;
t158 = t424 * t181 + t423 * t182;
t157 = pkin(6) * t197 - (-pkin(3) * t429 - pkin(7) * t426 - pkin(2)) * t212;
t156 = -t423 * t177 + t424 * t179 - t485;
t155 = -t423 * t176 + t424 * t178 - t486;
t154 = t424 * t177 + t423 * t179 + t451;
t153 = t424 * t176 + t423 * t178 + t452;
t152 = t427 * t172 - t430 * t173;
t151 = t430 * t172 + t427 * t173;
t150 = -pkin(6) * t181 + t429 * t166 - t426 * t180;
t149 = t427 * t158 - t430 * t159;
t148 = t430 * t158 + t427 * t159;
t147 = -pkin(2) * t194 + pkin(6) * t182 + t426 * t166 + t429 * t180;
t146 = -qJ(2) * t172 - t423 * t157 + t424 * t167;
t145 = pkin(1) * t212 + qJ(2) * t173 + t424 * t157 + t423 * t167;
t144 = -qJ(2) * t158 - t423 * t147 + t424 * t150;
t143 = -pkin(1) * t194 + qJ(2) * t159 + t424 * t147 + t423 * t150;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t422, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t422, 0, 0, 0, 0, 0, 0, t340, t341, 0, -t212, 0, 0, 0, 0, 0, 0, t340, t341, 0, t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t395, t394, 0, t352, 0, 0, 0, 0, 0, 0, -t500, t337, 0, t235, 0, 0, 0, 0, 0, 0, -t511, -t281, 0, t170, 0, 0, 0, 0, 0, 0, t205, t206, t227, t151, 0, 0, 0, 0, 0, 0, t205, t206, t227, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t394, t395, 0, t353, 0, 0, 0, 0, 0, 0, -t337, -t500, 0, -t515, 0, 0, 0, 0, 0, 0, t281, -t511, 0, -t521, 0, 0, 0, 0, 0, 0, t207, t208, t228, t152, 0, 0, 0, 0, 0, 0, t207, t208, t228, t149; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t402, t401, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t393 + t433, -pkin(1) * t392 - t330, 0, -pkin(1) * t290, 0, 0, 0, 0, 0, t418, -pkin(1) * t323 - pkin(2) * t383 - t286, -pkin(1) * t319 - pkin(2) * t380 - t287, 0, -pkin(1) * t189 - pkin(2) * t232, t332, t312, t342, t333, t346, 0, t441 - t470, t440 + t478, t213 + t439, pkin(1) * t172 + pkin(2) * t196 - pkin(3) * t274 + pkin(7) * t213, t332, t312, t342, t333, t346, 0, -qJ(5) * t477 + t428 * t229 + t441, t425 * t239 + t428 * t331 + t440, t428 * t225 + t425 * t234 + t439, pkin(1) * t158 + pkin(2) * t181 - pkin(3) * t240 + pkin(7) * t195 - qJ(5) * t479 + t428 * t198; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t394, 0, t395, 0, t448, -t369, -t353, -pkin(5) * t353, 0, 0, -t337, 0, -t500, 0, t513, t510, t515, pkin(5) * t515 + qJ(2) * t472 + t430 * t285, 0, 0, t281, 0, -t511, 0, t525, t524, t521, pkin(5) * t521 + t430 * t168 + t427 * t169, t215, t199, t219, t214, t218, t245, t430 * t162 + t427 * t164 - t494, t430 * t163 + t427 * t165 - t493, t430 * t174 + t427 * t175 - t492, -pkin(5) * t152 + t430 * t145 + t427 * t146, t215, t199, t219, t214, t218, t245, t430 * t153 + t427 * t155 - t494, t430 * t154 + t427 * t156 - t493, t430 * t160 + t427 * t161 - t492, -pkin(5) * t149 + t430 * t143 + t427 * t144; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t395, 0, -t394, 0, t369, t448, t352, pkin(5) * t352, 0, 0, t500, 0, -t337, 0, -t510, t513, t235, pkin(5) * t235 - qJ(2) * t465 + t427 * t285, 0, 0, t511, 0, t281, 0, -t524, t525, t170, pkin(5) * t170 + t427 * t168 - t430 * t169, t217, t200, t221, t216, t220, t246, t427 * t162 - t430 * t164 + t203, t427 * t163 - t430 * t165 + t204, t427 * t174 - t430 * t175 + t224, pkin(5) * t151 + t427 * t145 - t430 * t146, t217, t200, t221, t216, t220, t246, t427 * t153 - t430 * t155 + t203, t427 * t154 - t430 * t156 + t204, t427 * t160 - t430 * t161 + t224, pkin(5) * t148 + t427 * t143 - t430 * t144;];
tauB_reg = t1;
