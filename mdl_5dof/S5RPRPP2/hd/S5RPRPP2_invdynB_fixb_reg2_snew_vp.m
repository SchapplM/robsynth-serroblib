% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRPP2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRPP2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:15
% EndTime: 2019-12-31 18:11:21
% DurationCPUTime: 4.88s
% Computational Cost: add. (7053->387), mult. (14817->472), div. (0->0), fcn. (7453->6), ass. (0->255)
t417 = qJD(3) ^ 2;
t413 = sin(qJ(3));
t403 = t413 ^ 2;
t418 = qJD(1) ^ 2;
t481 = t403 * t418;
t382 = t417 + t481;
t415 = cos(qJ(3));
t450 = t415 * t418 * t413;
t378 = qJDD(3) - t450;
t461 = t415 * t378;
t331 = -t413 * t382 + t461;
t455 = qJD(1) * qJD(3);
t441 = t415 * t455;
t453 = t413 * qJDD(1);
t365 = 0.2e1 * t441 + t453;
t409 = sin(pkin(7));
t410 = cos(pkin(7));
t278 = t409 * t331 + t410 * t365;
t281 = t410 * t331 - t409 * t365;
t414 = sin(qJ(1));
t416 = cos(qJ(1));
t232 = t416 * t278 + t414 * t281;
t493 = pkin(5) * t232;
t236 = t414 * t278 - t416 * t281;
t230 = pkin(5) * t236;
t487 = qJ(2) * t278;
t448 = pkin(1) * t278 + pkin(2) * t365 + pkin(6) * t331;
t471 = t413 * t378;
t325 = t415 * t382 + t471;
t439 = -pkin(1) * t325 + qJ(2) * t281;
t377 = qJDD(3) + t450;
t357 = t415 * t377;
t383 = -t417 + t481;
t330 = t413 * t383 + t357;
t289 = t409 * t330 - t410 * t453;
t293 = t410 * t330 + t409 * t453;
t241 = t416 * t289 + t414 * t293;
t517 = t414 * t289 - t416 * t293;
t404 = t415 ^ 2;
t458 = t403 + t404;
t371 = t458 * qJDD(1);
t374 = t458 * t418;
t311 = t409 * t371 + t410 * t374;
t313 = t410 * t371 - t409 * t374;
t260 = t416 * t311 + t414 * t313;
t492 = pkin(5) * t260;
t263 = t414 * t311 - t416 * t313;
t259 = pkin(5) * t263;
t321 = t413 * t365;
t394 = t413 * t455;
t452 = t415 * qJDD(1);
t367 = -0.2e1 * t394 + t452;
t462 = t415 * t367;
t304 = -t462 + t321;
t375 = (t403 - t404) * t418;
t268 = t409 * t304 + t410 * t375;
t271 = t410 * t304 - t409 * t375;
t223 = t416 * t268 + t414 * t271;
t516 = t414 * t268 - t416 * t271;
t480 = t404 * t418;
t384 = -t417 + t480;
t329 = -t415 * t384 + t471;
t288 = t409 * t329 + t410 * t452;
t292 = t410 * t329 - t409 * t452;
t515 = t416 * t288 + t414 * t292;
t242 = t414 * t288 - t416 * t292;
t380 = t416 * g(1) + t414 * g(2);
t363 = -t418 * pkin(1) - t380;
t379 = t414 * g(1) - t416 * g(2);
t426 = qJDD(1) * pkin(1) + t379;
t300 = t409 * t363 - t410 * t426;
t301 = t410 * t363 + t409 * t426;
t438 = t409 * t300 + t410 * t301;
t255 = t410 * t300 - t409 * t301;
t460 = t416 * t255;
t514 = -t414 * t438 + t460;
t467 = t414 * t255;
t207 = t416 * t438 + t467;
t369 = t409 * qJDD(1) + t410 * t418;
t370 = t410 * qJDD(1) - t409 * t418;
t307 = -t414 * t369 + t416 * t370;
t405 = g(3) - qJDD(2);
t341 = qJ(2) * t369 - t410 * t405;
t428 = -qJ(2) * t370 - t409 * t405;
t513 = -pkin(5) * t307 + t414 * t341 + t416 * t428;
t510 = pkin(2) * t325;
t490 = pkin(6) * t325;
t486 = qJ(2) * t311;
t305 = qJ(2) * t313;
t447 = pkin(1) * t311 + pkin(2) * t374 + pkin(6) * t371;
t504 = t416 * t369 + t414 * t370;
t508 = pkin(5) * t504 + t416 * t341 - t414 * t428;
t322 = t413 * t384 + t461;
t385 = -t417 - t480;
t323 = t413 * t385 + t357;
t496 = pkin(2) * t323;
t502 = -qJ(4) * t385 - t496;
t457 = qJD(1) * t413;
t376 = -qJD(3) * pkin(4) - qJ(5) * t457;
t501 = t376 * t457 + qJDD(5);
t454 = qJD(4) * qJD(3);
t399 = -0.2e1 * t454;
t500 = -qJ(4) * t378 + t399 - t510;
t366 = -t394 + t452;
t499 = t366 * pkin(4) + t501;
t393 = t415 * t405;
t433 = qJDD(3) * pkin(3) + t417 * qJ(4) - qJDD(4) - t393;
t286 = -t418 * pkin(2) + qJDD(1) * pkin(6) + t301;
t431 = -pkin(3) * t415 - qJ(4) * t413;
t362 = t431 * qJD(1);
t436 = qJD(1) * t362 + t286;
t420 = t436 * t413 - t433;
t497 = pkin(3) + pkin(4);
t472 = t413 * t377;
t328 = t415 * t385 - t472;
t277 = t409 * t328 + t410 * t367;
t280 = t410 * t328 - t409 * t367;
t231 = t416 * t277 + t414 * t280;
t494 = pkin(5) * t231;
t491 = pkin(6) * t323;
t488 = qJ(2) * t277;
t485 = qJ(4) * t374;
t482 = qJ(4) * t415;
t285 = -qJDD(1) * pkin(2) - t418 * pkin(6) + t300;
t475 = t413 * t285;
t474 = t413 * t286;
t473 = t413 * t367;
t464 = t415 * t285;
t463 = t415 * t365;
t459 = -t415 * t286 + t413 * t405;
t456 = qJD(1) * t415;
t451 = 0.2e1 * t457;
t449 = pkin(1) * t277 + pkin(2) * t367 + pkin(6) * t328;
t446 = qJD(5) * t456;
t440 = -pkin(1) * t323 + qJ(2) * t280;
t269 = t393 + t474;
t225 = t413 * t269 - t415 * t459;
t319 = -t414 * t379 - t416 * t380;
t435 = t409 * t450;
t434 = t410 * t450;
t373 = t416 * qJDD(1) - t414 * t418;
t432 = -pkin(5) * t373 - t414 * g(3);
t430 = pkin(3) * t413 - t482;
t429 = t417 * pkin(3) - qJDD(3) * qJ(4) - t362 * t456 + t459;
t224 = t415 * t269 + t413 * t459;
t427 = t463 + t473;
t324 = -t415 * t383 + t472;
t318 = t416 * t379 - t414 * t380;
t398 = 0.2e1 * t454;
t249 = t398 - t429;
t425 = t441 + t453;
t424 = pkin(4) * t377 + t425 * qJ(5) + t433;
t423 = -t366 * pkin(3) + t285 + (-t425 - t441) * qJ(4);
t422 = pkin(4) * t480 - qJD(3) * t376 + t429;
t421 = qJD(4) * t451 - t423;
t419 = t366 * qJ(5) + t422;
t239 = (pkin(3) * qJD(3) - 0.2e1 * qJD(4)) * t457 + t423;
t228 = (t367 - t394) * pkin(3) + t421;
t227 = -pkin(3) * t394 + qJ(4) * t365 + t421;
t221 = t474 + (qJ(5) * qJD(3) * t415 + (-0.2e1 * qJD(5) + t362) * t413) * qJD(1) - t424;
t390 = 0.2e1 * t446;
t372 = t414 * qJDD(1) + t416 * t418;
t360 = t430 * qJDD(1);
t356 = t458 * t455;
t344 = -pkin(5) * t372 + t416 * g(3);
t342 = (-t497 * t413 + t482) * qJDD(1);
t337 = -t403 * t455 + t415 * t425;
t336 = -t413 * t366 - t404 * t455;
t335 = t409 * qJDD(3) + t410 * t356;
t334 = -t410 * qJDD(3) + t409 * t356;
t320 = (t366 - t394) * t415;
t317 = qJ(4) * t367 + qJ(5) * t377;
t298 = t410 * t337 - t435;
t297 = t410 * t336 + t435;
t296 = t409 * t337 + t434;
t295 = t409 * t336 - t434;
t284 = -qJ(5) * t378 + t497 * t365;
t265 = -t414 * t334 + t416 * t335;
t264 = t416 * t334 + t414 * t335;
t258 = t464 + t490;
t257 = t475 - t491;
t251 = -t459 + t510;
t250 = t269 - t496;
t248 = pkin(1) * t405 + qJ(2) * t438;
t247 = -t414 * t296 + t416 * t298;
t246 = -t414 * t295 + t416 * t297;
t245 = t416 * t296 + t414 * t298;
t244 = t416 * t295 + t414 * t297;
t238 = t420 + t485;
t237 = pkin(3) * t374 + t249;
t234 = -t414 * t277 + t416 * t280;
t229 = pkin(5) * t234;
t220 = t398 - t419 - 0.2e1 * t446;
t219 = -pkin(3) * t377 + t420 + t502;
t218 = qJ(5) * t480 + t239 - t499;
t217 = -pkin(3) * t382 + t429 + t500;
t216 = -qJ(5) * t441 - t485 + qJD(5) * t451 + (qJ(5) * qJDD(1) - t436) * t413 + t424;
t215 = t410 * t224 - t486;
t214 = t409 * t224 + t305;
t213 = t390 + t399 - t497 * t374 + (t366 + t452) * qJ(5) + t422;
t212 = (t382 - t480) * qJ(5) + t227 + t499;
t211 = -pkin(3) * t321 + t415 * t227 - t490;
t210 = qJ(4) * t462 - t413 * t228 - t491;
t209 = t410 * t225 + t409 * t285;
t208 = t409 * t225 - t410 * t285;
t205 = (-t385 - t480) * qJ(5) + (t366 + t367) * pkin(4) + t228 + t501;
t204 = t415 * t249 + t413 * t420;
t203 = t413 * t249 - t415 * t420;
t202 = -t497 * t382 + t390 + t419 + t500;
t201 = -t497 * t377 + t221 + t502;
t200 = -t413 * t237 + t415 * t238;
t199 = -t409 * t251 + t410 * t258 + t487;
t198 = -t409 * t250 + t410 * t257 - t488;
t197 = -t413 * t205 + t415 * t317 - t491;
t196 = t415 * t212 - t413 * t284 - t490;
t195 = t410 * t251 + t409 * t258 - t439;
t194 = t410 * t250 + t409 * t257 + t440;
t193 = t410 * t200 - t409 * t360 - t486;
t192 = t409 * t200 + t410 * t360 + t305;
t191 = t415 * t220 + t413 * t221;
t190 = t413 * t220 - t415 * t221;
t189 = -qJ(4) * t218 - qJ(5) * t221;
t188 = t410 * t204 + t409 * t239;
t187 = t409 * t204 - t410 * t239;
t186 = -t413 * t213 + t415 * t216;
t185 = -t414 * t208 + t416 * t209;
t184 = t416 * t208 + t414 * t209;
t183 = -pkin(2) * t203 + pkin(3) * t420 - qJ(4) * t249;
t182 = t410 * t210 - t409 * t219 - t488;
t181 = t410 * t211 - t409 * t217 - t487;
t180 = t410 * t186 - t409 * t342 + t486;
t179 = t409 * t186 + t410 * t342 - t305;
t178 = -pkin(6) * t203 + t430 * t239;
t177 = t409 * t210 + t410 * t219 + t440;
t176 = t409 * t211 + t410 * t217 + t439;
t175 = -qJ(5) * t220 - t497 * t218;
t174 = t410 * t191 + t409 * t218;
t173 = t409 * t191 - t410 * t218;
t172 = -qJ(2) * t208 - (pkin(2) * t409 - pkin(6) * t410) * t224;
t171 = t410 * t196 - t409 * t202 - t487;
t170 = t410 * t197 - t409 * t201 - t488;
t169 = t409 * t196 + t410 * t202 + t439;
t168 = t409 * t197 + t410 * t201 + t440;
t167 = qJ(2) * t209 - (-pkin(2) * t410 - pkin(6) * t409 - pkin(1)) * t224;
t166 = -t414 * t187 + t416 * t188;
t165 = t416 * t187 + t414 * t188;
t164 = -pkin(2) * t190 - qJ(4) * t220 + t497 * t221;
t163 = -t414 * t173 + t416 * t174;
t162 = t416 * t173 + t414 * t174;
t161 = -pkin(6) * t190 - t413 * t175 + t415 * t189;
t160 = -qJ(2) * t187 + t410 * t178 - t409 * t183;
t159 = -pkin(1) * t203 + qJ(2) * t188 + t409 * t178 + t410 * t183;
t158 = -qJ(2) * t173 + t410 * t161 - t409 * t164;
t157 = -pkin(1) * t190 + qJ(2) * t174 + t409 * t161 + t410 * t164;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t372, -t373, 0, t319, 0, 0, 0, 0, 0, 0, -t504, -t307, 0, t207, 0, 0, 0, 0, 0, 0, t234, t236, -t263, t185, 0, 0, 0, 0, 0, 0, t234, -t263, -t236, t166, 0, 0, 0, 0, 0, 0, t234, -t236, t263, t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t373, -t372, 0, t318, 0, 0, 0, 0, 0, 0, t307, -t504, 0, -t514, 0, 0, 0, 0, 0, 0, t231, -t232, t260, t184, 0, 0, 0, 0, 0, 0, t231, t260, t232, t165, 0, 0, 0, 0, 0, 0, t231, t232, -t260, t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t405, 0, 0, 0, 0, 0, 0, t323, -t325, 0, -t224, 0, 0, 0, 0, 0, 0, t323, 0, t325, t203, 0, 0, 0, 0, 0, 0, t323, t325, 0, t190; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t373, 0, -t372, 0, t432, -t344, -t318, -pkin(5) * t318, 0, 0, t307, 0, -t504, 0, t513, t508, t514, pkin(5) * t514 + qJ(2) * t460 - t414 * t248, t247, t516, -t517, t246, t242, t265, -t414 * t194 + t416 * t198 - t494, -t414 * t195 + t416 * t199 + t493, -t414 * t214 + t416 * t215 - t492, -pkin(5) * t184 - t414 * t167 + t416 * t172, t247, -t517, -t516, t265, -t242, t246, -t414 * t177 + t416 * t182 - t494, -t414 * t192 + t416 * t193 - t492, -t414 * t176 + t416 * t181 - t493, -pkin(5) * t165 - t414 * t159 + t416 * t160, t247, -t516, t517, t246, t242, t265, -t414 * t168 + t416 * t170 - t494, -t414 * t169 + t416 * t171 - t493, -t414 * t179 + t416 * t180 + t492, -pkin(5) * t162 - t414 * t157 + t416 * t158; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t372, 0, t373, 0, t344, t432, t319, pkin(5) * t319, 0, 0, t504, 0, t307, 0, -t508, t513, t207, pkin(5) * t207 + qJ(2) * t467 + t416 * t248, t245, -t223, t241, t244, -t515, t264, t416 * t194 + t414 * t198 + t229, t416 * t195 + t414 * t199 + t230, t416 * t214 + t414 * t215 - t259, pkin(5) * t185 + t416 * t167 + t414 * t172, t245, t241, t223, t264, t515, t244, t416 * t177 + t414 * t182 + t229, t416 * t192 + t414 * t193 - t259, t416 * t176 + t414 * t181 - t230, pkin(5) * t166 + t416 * t159 + t414 * t160, t245, t223, -t241, t244, -t515, t264, t416 * t168 + t414 * t170 + t229, t416 * t169 + t414 * t171 - t230, t416 * t179 + t414 * t180 + t259, pkin(5) * t163 + t416 * t157 + t414 * t158; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t379, t380, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t370 - t300, -pkin(1) * t369 - t301, 0, -pkin(1) * t255, t321, t427, t324, t320, t322, 0, t449 - t464, -t448 + t475, t225 + t447, pkin(1) * t208 - pkin(2) * t285 + pkin(6) * t225, t321, t324, -t427, 0, -t322, t320, qJ(4) * t473 + t415 * t228 + t449, t415 * t237 + t413 * t238 + t447, pkin(3) * t463 + t413 * t227 + t448, pkin(1) * t187 + pkin(6) * t204 + (-pkin(2) + t431) * t239, t321, -t427, -t324, t320, t322, 0, t415 * t205 + t413 * t317 + t449, t413 * t212 + t415 * t284 + t448, t415 * t213 + t413 * t216 - t447, pkin(1) * t173 - pkin(2) * t218 + pkin(6) * t191 + t415 * t175 + t413 * t189;];
tauB_reg = t1;
