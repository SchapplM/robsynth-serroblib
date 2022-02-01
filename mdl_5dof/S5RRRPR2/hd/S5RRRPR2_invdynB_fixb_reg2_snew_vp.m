% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:56
% EndTime: 2022-01-20 11:31:04
% DurationCPUTime: 7.71s
% Computational Cost: add. (37767->388), mult. (48058->558), div. (0->0), fcn. (28562->10), ass. (0->254)
t426 = qJD(1) + qJD(2);
t421 = qJD(3) + t426;
t419 = t421 ^ 2;
t425 = qJDD(1) + qJDD(2);
t420 = qJDD(3) + t425;
t430 = sin(pkin(9));
t431 = cos(pkin(9));
t377 = t419 * t431 + t420 * t430;
t380 = t419 * t430 - t420 * t431;
t433 = sin(qJ(3));
t437 = cos(qJ(3));
t321 = t377 * t437 - t380 * t433;
t429 = g(3) - qJDD(4);
t359 = qJ(4) * t377 - t429 * t431;
t489 = qJ(4) * t380 - t429 * t430;
t271 = pkin(7) * t321 + t359 * t437 - t433 * t489;
t325 = t377 * t433 + t380 * t437;
t434 = sin(qJ(2));
t438 = cos(qJ(2));
t280 = t321 * t438 - t325 * t434;
t503 = pkin(7) * t325 + t359 * t433 + t437 * t489;
t215 = pkin(6) * t280 + t271 * t438 - t434 * t503;
t284 = t321 * t434 + t325 * t438;
t435 = sin(qJ(1));
t439 = cos(qJ(1));
t240 = t280 * t435 + t284 * t439;
t518 = pkin(6) * t284 + t271 * t434 + t438 * t503;
t527 = pkin(5) * t240 + t215 * t435 + t439 * t518;
t507 = t280 * t439 - t284 * t435;
t526 = pkin(5) * t507 + t215 * t439 - t435 * t518;
t411 = g(1) * t435 - g(2) * t439;
t401 = qJDD(1) * pkin(1) + t411;
t412 = g(1) * t439 + g(2) * t435;
t442 = qJD(1) ^ 2;
t402 = -pkin(1) * t442 - t412;
t346 = -t401 * t438 + t402 * t434;
t337 = pkin(2) * t425 - t346;
t347 = t401 * t434 + t402 * t438;
t424 = t426 ^ 2;
t338 = -pkin(2) * t424 + t347;
t298 = t337 * t433 + t338 * t437;
t295 = -pkin(3) * t419 + t298;
t445 = t337 * t437 - t338 * t433;
t443 = pkin(3) * t420 + t445;
t247 = t295 * t430 - t431 * t443;
t248 = t295 * t431 + t430 * t443;
t453 = t247 * t430 + t248 * t431;
t206 = t247 * t431 - t248 * t430;
t471 = t206 * t437;
t180 = -t433 * t453 + t471;
t472 = t206 * t433;
t492 = t437 * t453 + t472;
t167 = t180 * t438 - t434 * t492;
t513 = t180 * t434 + t438 * t492;
t153 = t167 * t435 + t439 * t513;
t523 = t167 * t439 - t435 * t513;
t383 = t419 * t437 + t420 * t433;
t386 = t419 * t433 - t420 * t437;
t329 = t383 * t438 - t386 * t434;
t365 = pkin(7) * t383 - g(3) * t437;
t490 = pkin(7) * t386 - g(3) * t433;
t279 = pkin(6) * t329 + t365 * t438 - t434 * t490;
t333 = t383 * t434 + t386 * t438;
t290 = t329 * t435 + t333 * t439;
t504 = pkin(6) * t333 + t365 * t434 + t438 * t490;
t520 = pkin(5) * t290 + t279 * t435 + t439 * t504;
t487 = t329 * t439 - t333 * t435;
t519 = pkin(5) * t487 + t279 * t439 - t435 * t504;
t452 = t298 * t437 - t433 * t445;
t252 = -t298 * t433 - t437 * t445;
t467 = t252 * t438;
t210 = -t434 * t452 + t467;
t468 = t252 * t434;
t493 = t438 * t452 + t468;
t184 = t210 * t435 + t439 * t493;
t512 = t210 * t439 - t435 * t493;
t394 = t424 * t438 + t425 * t434;
t397 = t424 * t434 - t425 * t438;
t341 = t394 * t435 + t397 * t439;
t370 = pkin(6) * t394 - g(3) * t438;
t491 = pkin(6) * t397 - g(3) * t434;
t506 = pkin(5) * t341 + t370 * t435 + t439 * t491;
t444 = t394 * t439 - t397 * t435;
t505 = pkin(5) * t444 + t370 * t439 - t435 * t491;
t451 = t346 * t434 + t347 * t438;
t304 = t346 * t438 - t347 * t434;
t465 = t304 * t439;
t494 = -t435 * t451 + t465;
t466 = t304 * t435;
t257 = t439 * t451 + t466;
t242 = -pkin(4) * t420 - pkin(8) * t419 + t247;
t432 = sin(qJ(5));
t470 = t242 * t432;
t436 = cos(qJ(5));
t469 = t242 * t436;
t407 = t436 * t419 * t432;
t398 = qJDD(5) + t407;
t464 = t398 * t432;
t399 = qJDD(5) - t407;
t463 = t399 * t432;
t462 = t399 * t436;
t427 = t432 ^ 2;
t461 = t419 * t427;
t460 = t420 * t432;
t413 = t436 * t420;
t243 = -pkin(4) * t419 + pkin(8) * t420 + t248;
t236 = t243 * t436 - t429 * t432;
t428 = t436 ^ 2;
t457 = t427 + t428;
t456 = qJD(5) * t421;
t455 = t432 * t456;
t454 = t436 * t456;
t235 = t243 * t432 + t429 * t436;
t200 = t235 * t432 + t236 * t436;
t361 = -t411 * t435 - t412 * t439;
t448 = t430 * t407;
t447 = t431 * t407;
t409 = qJDD(1) * t439 - t435 * t442;
t446 = -pkin(5) * t409 - g(3) * t435;
t199 = t235 * t436 - t236 * t432;
t360 = t411 * t439 - t412 * t435;
t441 = qJD(5) ^ 2;
t440 = pkin(1) * g(3);
t414 = t428 * t419;
t408 = qJDD(1) * t435 + t439 * t442;
t406 = -t414 - t441;
t405 = t414 - t441;
t404 = -t441 - t461;
t403 = t441 - t461;
t390 = t436 * t398;
t389 = -pkin(5) * t408 + g(3) * t439;
t388 = t414 - t461;
t387 = t414 + t461;
t381 = t457 * t420;
t376 = t413 - 0.2e1 * t455;
t375 = t413 - t455;
t374 = t454 + t460;
t373 = 0.2e1 * t454 + t460;
t372 = t457 * t456;
t355 = qJDD(5) * t430 + t372 * t431;
t354 = -qJDD(5) * t431 + t372 * t430;
t353 = -t404 * t432 - t462;
t352 = -t403 * t432 + t390;
t351 = t406 * t436 - t464;
t350 = t405 * t436 - t463;
t349 = t404 * t436 - t463;
t348 = t406 * t432 + t390;
t345 = t374 * t436 - t427 * t456;
t344 = -t375 * t432 - t428 * t456;
t328 = t381 * t431 - t387 * t430;
t327 = t381 * t430 + t387 * t431;
t320 = -t373 * t432 + t376 * t436;
t319 = t352 * t431 + t430 * t460;
t318 = t350 * t431 + t413 * t430;
t317 = t352 * t430 - t431 * t460;
t316 = t350 * t430 - t413 * t431;
t315 = t345 * t431 - t448;
t314 = t344 * t431 + t448;
t313 = t345 * t430 + t447;
t312 = t344 * t430 - t447;
t311 = t353 * t431 + t373 * t430;
t310 = t351 * t431 - t376 * t430;
t309 = t353 * t430 - t373 * t431;
t308 = t351 * t430 + t376 * t431;
t307 = -t354 * t433 + t355 * t437;
t306 = t354 * t437 + t355 * t433;
t303 = t320 * t431 - t388 * t430;
t300 = t320 * t430 + t388 * t431;
t299 = pkin(6) * t451 + t440;
t287 = -t327 * t433 + t328 * t437;
t286 = t327 * t437 + t328 * t433;
t275 = -t317 * t433 + t319 * t437;
t274 = -t316 * t433 + t318 * t437;
t273 = t317 * t437 + t319 * t433;
t272 = t316 * t437 + t318 * t433;
t267 = -t313 * t433 + t315 * t437;
t266 = -t312 * t433 + t314 * t437;
t265 = t313 * t437 + t315 * t433;
t264 = t312 * t437 + t314 * t433;
t263 = -t309 * t433 + t311 * t437;
t262 = -t308 * t433 + t310 * t437;
t261 = t309 * t437 + t311 * t433;
t260 = t308 * t437 + t310 * t433;
t259 = -t306 * t434 + t307 * t438;
t258 = t306 * t438 + t307 * t434;
t256 = -t300 * t433 + t303 * t437;
t254 = t300 * t437 + t303 * t433;
t249 = pkin(2) * g(3) + pkin(7) * t452;
t245 = -t286 * t434 + t287 * t438;
t244 = t286 * t438 + t287 * t434;
t233 = -t273 * t434 + t275 * t438;
t232 = -t272 * t434 + t274 * t438;
t231 = t273 * t438 + t275 * t434;
t230 = t272 * t438 + t274 * t434;
t229 = -t265 * t434 + t267 * t438;
t228 = -t264 * t434 + t266 * t438;
t227 = t265 * t438 + t267 * t434;
t226 = t264 * t438 + t266 * t434;
t225 = -pkin(8) * t349 + t469;
t224 = -pkin(8) * t348 + t470;
t223 = -pkin(4) * t349 + t236;
t222 = -pkin(4) * t348 + t235;
t221 = -t261 * t434 + t263 * t438;
t220 = -t260 * t434 + t262 * t438;
t219 = t261 * t438 + t263 * t434;
t218 = t260 * t438 + t262 * t434;
t217 = -t254 * t434 + t256 * t438;
t216 = t254 * t438 + t256 * t434;
t203 = pkin(3) * t429 + qJ(4) * t453;
t202 = -t244 * t435 + t245 * t439;
t201 = t244 * t439 + t245 * t435;
t197 = -qJ(4) * t327 + t199 * t431;
t196 = qJ(4) * t328 + t199 * t430;
t195 = -t219 * t435 + t221 * t439;
t194 = -t218 * t435 + t220 * t439;
t193 = t219 * t439 + t221 * t435;
t192 = t218 * t439 + t220 * t435;
t191 = -qJ(4) * t309 - t223 * t430 + t225 * t431;
t190 = -qJ(4) * t308 - t222 * t430 + t224 * t431;
t189 = -pkin(3) * t349 + qJ(4) * t311 + t223 * t431 + t225 * t430;
t188 = -pkin(3) * t348 + qJ(4) * t310 + t222 * t431 + t224 * t430;
t187 = t200 * t431 + t242 * t430;
t186 = t200 * t430 - t242 * t431;
t185 = pkin(6) * t210 + pkin(7) * t467 - t249 * t434;
t182 = pkin(6) * t493 + pkin(7) * t468 + t249 * t438 + t440;
t177 = -pkin(7) * t286 - t196 * t433 + t197 * t437;
t176 = pkin(7) * t287 + t196 * t437 + t197 * t433;
t175 = -pkin(7) * t261 - t189 * t433 + t191 * t437;
t174 = -pkin(7) * t260 - t188 * t433 + t190 * t437;
t173 = -t186 * t433 + t187 * t437;
t172 = t186 * t437 + t187 * t433;
t171 = -pkin(2) * t349 + pkin(7) * t263 + t189 * t437 + t191 * t433;
t170 = -pkin(2) * t348 + pkin(7) * t262 + t188 * t437 + t190 * t433;
t169 = -qJ(4) * t186 - (pkin(4) * t430 - pkin(8) * t431) * t199;
t164 = pkin(7) * t180 + qJ(4) * t471 - t203 * t433;
t163 = pkin(2) * t429 + pkin(7) * t492 + qJ(4) * t472 + t203 * t437;
t162 = qJ(4) * t187 - (-pkin(4) * t431 - pkin(8) * t430 - pkin(3)) * t199;
t161 = -pkin(6) * t244 - t176 * t434 + t177 * t438;
t160 = pkin(6) * t245 + t176 * t438 + t177 * t434;
t159 = -t172 * t434 + t173 * t438;
t158 = t172 * t438 + t173 * t434;
t157 = -pkin(6) * t219 - t171 * t434 + t175 * t438;
t156 = -pkin(6) * t218 - t170 * t434 + t174 * t438;
t155 = -pkin(1) * t349 + pkin(6) * t221 + t171 * t438 + t175 * t434;
t154 = -pkin(1) * t348 + pkin(6) * t220 + t170 * t438 + t174 * t434;
t151 = pkin(6) * t167 - t163 * t434 + t164 * t438;
t150 = pkin(1) * t429 + pkin(6) * t513 + t163 * t438 + t164 * t434;
t149 = -pkin(7) * t172 - t162 * t433 + t169 * t437;
t148 = -t158 * t435 + t159 * t439;
t147 = t158 * t439 + t159 * t435;
t146 = pkin(2) * t199 + pkin(7) * t173 + t162 * t437 + t169 * t433;
t145 = -pkin(6) * t158 - t146 * t434 + t149 * t438;
t144 = pkin(1) * t199 + pkin(6) * t159 + t146 * t438 + t149 * t434;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t408, -t409, 0, t361, 0, 0, 0, 0, 0, 0, -t444, t341, 0, t257, 0, 0, 0, 0, 0, 0, -t487, t290, 0, t184, 0, 0, 0, 0, 0, 0, -t507, t240, 0, t153, 0, 0, 0, 0, 0, 0, t194, t195, t202, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t409, -t408, 0, t360, 0, 0, 0, 0, 0, 0, -t341, -t444, 0, -t494, 0, 0, 0, 0, 0, 0, -t290, -t487, 0, -t512, 0, 0, 0, 0, 0, 0, -t240, -t507, 0, -t523, 0, 0, 0, 0, 0, 0, t192, t193, t201, t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t429, 0, 0, 0, 0, 0, 0, t348, t349, 0, -t199; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t409, 0, -t408, 0, t446, -t389, -t360, -pkin(5) * t360, 0, 0, -t341, 0, -t444, 0, t506, t505, t494, pkin(5) * t494 + pkin(6) * t465 - t299 * t435, 0, 0, -t290, 0, -t487, 0, t520, t519, t512, pkin(5) * t512 - t182 * t435 + t185 * t439, 0, 0, -t240, 0, -t507, 0, t527, t526, t523, pkin(5) * t523 - t150 * t435 + t151 * t439, -t227 * t435 + t229 * t439, -t216 * t435 + t217 * t439, -t231 * t435 + t233 * t439, -t226 * t435 + t228 * t439, -t230 * t435 + t232 * t439, -t258 * t435 + t259 * t439, -pkin(5) * t192 - t154 * t435 + t156 * t439, -pkin(5) * t193 - t155 * t435 + t157 * t439, -pkin(5) * t201 - t160 * t435 + t161 * t439, -pkin(5) * t147 - t144 * t435 + t145 * t439; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t408, 0, t409, 0, t389, t446, t361, pkin(5) * t361, 0, 0, t444, 0, -t341, 0, -t505, t506, t257, pkin(5) * t257 + pkin(6) * t466 + t299 * t439, 0, 0, t487, 0, -t290, 0, -t519, t520, t184, pkin(5) * t184 + t182 * t439 + t185 * t435, 0, 0, t507, 0, -t240, 0, -t526, t527, t153, pkin(5) * t153 + t150 * t439 + t151 * t435, t227 * t439 + t229 * t435, t216 * t439 + t217 * t435, t231 * t439 + t233 * t435, t226 * t439 + t228 * t435, t230 * t439 + t232 * t435, t258 * t439 + t259 * t435, pkin(5) * t194 + t154 * t439 + t156 * t435, pkin(5) * t195 + t155 * t439 + t157 * t435, pkin(5) * t202 + t160 * t439 + t161 * t435, pkin(5) * t148 + t144 * t439 + t145 * t435; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t411, t412, 0, 0, 0, 0, 0, 0, 0, t425, -pkin(1) * t397 - t346, -pkin(1) * t394 - t347, 0, -pkin(1) * t304, 0, 0, 0, 0, 0, t420, -pkin(1) * t333 - pkin(2) * t386 + t445, -pkin(1) * t329 - pkin(2) * t383 - t298, 0, -pkin(1) * t210 - pkin(2) * t252, 0, 0, 0, 0, 0, t420, -pkin(1) * t284 - pkin(2) * t325 - pkin(3) * t380 - t247, -pkin(1) * t280 - pkin(2) * t321 - pkin(3) * t377 - t248, 0, -pkin(1) * t167 - pkin(2) * t180 - pkin(3) * t206, (t374 + t454) * t432, t373 * t436 + t376 * t432, t403 * t436 + t464, (t375 - t455) * t436, t405 * t432 + t462, 0, pkin(1) * t218 + pkin(2) * t260 + pkin(3) * t308 + pkin(4) * t376 + pkin(8) * t351 - t469, pkin(1) * t219 + pkin(2) * t261 + pkin(3) * t309 - pkin(4) * t373 + pkin(8) * t353 + t470, pkin(1) * t244 + pkin(2) * t286 + pkin(3) * t327 + pkin(4) * t387 + pkin(8) * t381 + t200, pkin(1) * t158 + pkin(2) * t172 + pkin(3) * t186 - pkin(4) * t242 + pkin(8) * t200;];
tauB_reg = t1;
