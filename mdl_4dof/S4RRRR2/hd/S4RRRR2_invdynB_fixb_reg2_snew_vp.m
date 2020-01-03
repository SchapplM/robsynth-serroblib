% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRR2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:21
% EndTime: 2019-12-31 17:23:27
% DurationCPUTime: 4.24s
% Computational Cost: add. (15584->375), mult. (20830->559), div. (0->0), fcn. (13277->8), ass. (0->261)
t413 = qJD(1) + qJD(2);
t409 = t413 ^ 2;
t422 = cos(qJ(2));
t411 = qJDD(1) + qJDD(2);
t418 = sin(qJ(2));
t450 = t418 * t411;
t379 = t422 * t409 + t450;
t443 = t422 * t411;
t382 = t418 * t409 - t443;
t419 = sin(qJ(1));
t423 = cos(qJ(1));
t334 = t419 * t379 + t423 * t382;
t359 = pkin(5) * t379 - t422 * g(3);
t470 = pkin(5) * t382 - t418 * g(3);
t477 = pkin(4) * t334 + t419 * t359 + t423 * t470;
t428 = t423 * t379 - t419 * t382;
t476 = pkin(4) * t428 + t423 * t359 - t419 * t470;
t399 = t423 * g(1) + t419 * g(2);
t425 = qJD(1) ^ 2;
t386 = -t425 * pkin(1) - t399;
t398 = t419 * g(1) - t423 * g(2);
t427 = qJDD(1) * pkin(1) + t398;
t337 = t418 * t386 - t422 * t427;
t338 = t422 * t386 + t418 * t427;
t433 = t418 * t337 + t422 * t338;
t296 = t422 * t337 - t418 * t338;
t442 = t423 * t296;
t473 = -t419 * t433 + t442;
t449 = t419 * t296;
t252 = t423 * t433 + t449;
t416 = sin(qJ(4));
t420 = cos(qJ(4));
t421 = cos(qJ(3));
t417 = sin(qJ(3));
t460 = t413 * t417;
t362 = -t420 * t421 * t413 + t416 * t460;
t364 = (t416 * t421 + t417 * t420) * t413;
t328 = t364 * t362;
t438 = qJDD(3) + qJDD(4);
t467 = -t328 + t438;
t472 = t416 * t467;
t471 = t420 * t467;
t412 = qJD(3) + qJD(4);
t355 = t412 * t362;
t440 = qJD(3) * t413;
t434 = t421 * t440;
t451 = t417 * t411;
t371 = t434 + t451;
t402 = t421 * t411;
t435 = t417 * t440;
t372 = t402 - t435;
t426 = t362 * qJD(4) - t420 * t371 - t416 * t372;
t466 = -t355 - t426;
t432 = t416 * t371 - t420 * t372;
t276 = (qJD(4) - t412) * t364 + t432;
t360 = t362 ^ 2;
t361 = t364 ^ 2;
t408 = t412 ^ 2;
t463 = t409 * t417;
t462 = t412 * t416;
t461 = t412 * t420;
t414 = t417 ^ 2;
t459 = t414 * t409;
t415 = t421 ^ 2;
t403 = t415 * t409;
t325 = -t411 * pkin(2) - t409 * pkin(6) + t337;
t389 = qJD(3) * pkin(3) - pkin(7) * t460;
t274 = -t372 * pkin(3) - pkin(7) * t403 + t389 * t460 + t325;
t458 = t416 * t274;
t319 = t328 + t438;
t457 = t416 * t319;
t326 = -t409 * pkin(2) + t411 * pkin(6) + t338;
t454 = t417 * t326;
t272 = qJDD(3) * pkin(3) - t371 * pkin(7) - t454 + (pkin(3) * t463 + pkin(7) * t440 - g(3)) * t421;
t311 = -t417 * g(3) + t421 * t326;
t273 = -pkin(3) * t403 + t372 * pkin(7) - qJD(3) * t389 + t311;
t236 = -t420 * t272 + t416 * t273;
t237 = t416 * t272 + t420 * t273;
t205 = -t420 * t236 + t416 * t237;
t456 = t417 * t205;
t455 = t417 * t325;
t397 = t421 * t463;
t387 = qJDD(3) + t397;
t453 = t417 * t387;
t388 = qJDD(3) - t397;
t452 = t417 * t388;
t448 = t420 * t274;
t447 = t420 * t319;
t446 = t421 * t205;
t445 = t421 * t325;
t444 = t421 * t388;
t441 = t414 + t415;
t437 = t418 * t328;
t436 = t422 * t328;
t206 = t416 * t236 + t420 * t237;
t310 = t421 * g(3) + t454;
t262 = t417 * t310 + t421 * t311;
t351 = -t419 * t398 - t423 * t399;
t431 = t418 * t397;
t430 = t422 * t397;
t391 = t423 * qJDD(1) - t419 * t425;
t429 = -pkin(4) * t391 - t419 * g(3);
t261 = t421 * t310 - t417 * t311;
t350 = t423 * t398 - t419 * t399;
t424 = qJD(3) ^ 2;
t395 = -t403 - t424;
t394 = t403 - t424;
t393 = -t424 - t459;
t392 = t424 - t459;
t390 = t419 * qJDD(1) + t423 * t425;
t384 = t403 - t459;
t383 = t403 + t459;
t378 = t421 * t387;
t377 = t441 * t411;
t373 = t402 - 0.2e1 * t435;
t370 = 0.2e1 * t434 + t451;
t369 = -pkin(4) * t390 + t423 * g(3);
t368 = t441 * t440;
t353 = -t361 + t408;
t352 = t360 - t408;
t349 = t418 * qJDD(3) + t422 * t368;
t348 = -t422 * qJDD(3) + t418 * t368;
t347 = t421 * t371 - t414 * t440;
t346 = -t417 * t372 - t415 * t440;
t345 = -t361 - t408;
t344 = -t417 * t393 - t444;
t343 = -t417 * t392 + t378;
t342 = t421 * t395 - t453;
t341 = t421 * t394 - t452;
t340 = t421 * t393 - t452;
t339 = t417 * t395 + t378;
t333 = t422 * t377 - t418 * t383;
t330 = t418 * t377 + t422 * t383;
t329 = -t417 * t370 + t421 * t373;
t327 = -t361 + t360;
t324 = t422 * t343 + t417 * t450;
t323 = t422 * t341 + t418 * t402;
t322 = t418 * t343 - t417 * t443;
t321 = t418 * t341 - t421 * t443;
t317 = t422 * t347 - t431;
t316 = t422 * t346 + t431;
t315 = t418 * t347 + t430;
t314 = t418 * t346 - t430;
t313 = -t408 - t360;
t309 = t422 * t344 + t418 * t370;
t308 = t422 * t342 - t418 * t373;
t307 = t418 * t344 - t422 * t370;
t306 = t418 * t342 + t422 * t373;
t305 = (-t362 * t420 + t364 * t416) * t412;
t304 = (-t362 * t416 - t364 * t420) * t412;
t302 = -t360 - t361;
t301 = t422 * t329 - t418 * t384;
t300 = t418 * t329 + t422 * t384;
t298 = -t364 * qJD(4) - t432;
t293 = pkin(1) * g(3) + pkin(5) * t433;
t292 = t420 * t352 - t457;
t291 = -t416 * t353 + t471;
t290 = t416 * t352 + t447;
t289 = t420 * t353 + t472;
t288 = -pkin(6) * t340 + t445;
t287 = -pkin(6) * t339 + t455;
t286 = -t416 * t345 - t447;
t285 = t420 * t345 - t457;
t284 = -pkin(2) * t340 + t311;
t283 = -pkin(2) * t339 + t310;
t282 = -t419 * t330 + t423 * t333;
t281 = t423 * t330 + t419 * t333;
t280 = -t355 + t426;
t275 = (qJD(4) + t412) * t364 + t432;
t271 = -t364 * t462 - t420 * t426;
t270 = t364 * t461 - t416 * t426;
t269 = -t416 * t298 + t362 * t461;
t268 = t420 * t298 + t362 * t462;
t264 = t420 * t313 - t472;
t263 = t416 * t313 + t471;
t259 = -t419 * t307 + t423 * t309;
t258 = -t419 * t306 + t423 * t308;
t257 = t423 * t307 + t419 * t309;
t256 = t423 * t306 + t419 * t308;
t255 = -t417 * t304 + t421 * t305;
t254 = t422 * t255 + t418 * t438;
t253 = t418 * t255 - t422 * t438;
t250 = -pkin(5) * t330 + t422 * t261;
t249 = pkin(5) * t333 + t418 * t261;
t248 = t422 * t262 + t418 * t325;
t247 = t418 * t262 - t422 * t325;
t246 = -t417 * t290 + t421 * t292;
t245 = -t417 * t289 + t421 * t291;
t244 = -t417 * t285 + t421 * t286;
t243 = t421 * t285 + t417 * t286;
t242 = -pkin(7) * t285 + t448;
t241 = -t276 * t420 - t416 * t280;
t240 = -t420 * t275 - t416 * t466;
t239 = -t276 * t416 + t420 * t280;
t238 = -t416 * t275 + t420 * t466;
t235 = -pkin(7) * t263 + t458;
t233 = -t417 * t270 + t421 * t271;
t232 = -t417 * t268 + t421 * t269;
t231 = -t417 * t263 + t421 * t264;
t230 = t421 * t263 + t417 * t264;
t229 = -pkin(5) * t307 - t418 * t284 + t422 * t288;
t228 = -pkin(5) * t306 - t418 * t283 + t422 * t287;
t227 = t422 * t233 + t437;
t226 = t422 * t232 - t437;
t225 = t418 * t233 - t436;
t224 = t418 * t232 + t436;
t223 = -pkin(1) * t340 + pkin(5) * t309 + t422 * t284 + t418 * t288;
t222 = -pkin(1) * t339 + pkin(5) * t308 + t422 * t283 + t418 * t287;
t221 = t422 * t246 - t418 * t276;
t220 = t422 * t245 - t418 * t280;
t219 = t418 * t246 + t422 * t276;
t218 = t418 * t245 + t422 * t280;
t217 = -pkin(3) * t466 + pkin(7) * t286 + t458;
t216 = t422 * t244 + t418 * t466;
t215 = t418 * t244 - t422 * t466;
t214 = -pkin(3) * t275 + pkin(7) * t264 - t448;
t213 = t422 * t231 + t418 * t275;
t212 = t418 * t231 - t422 * t275;
t211 = -t419 * t247 + t423 * t248;
t210 = t423 * t247 + t419 * t248;
t209 = -t417 * t239 + t421 * t241;
t208 = -t417 * t238 + t421 * t240;
t207 = t421 * t239 + t417 * t241;
t204 = -pkin(5) * t247 - (pkin(2) * t418 - pkin(6) * t422) * t261;
t203 = t422 * t208 - t418 * t327;
t202 = t418 * t208 + t422 * t327;
t201 = t422 * t209 + t418 * t302;
t200 = t418 * t209 - t422 * t302;
t199 = -pkin(2) * t243 - pkin(3) * t285 + t237;
t198 = -pkin(3) * t274 + pkin(7) * t206;
t197 = -pkin(2) * t230 - pkin(3) * t263 + t236;
t196 = pkin(5) * t248 - (-pkin(2) * t422 - pkin(6) * t418 - pkin(1)) * t261;
t195 = -t419 * t215 + t423 * t216;
t194 = t423 * t215 + t419 * t216;
t193 = -pkin(2) * t207 - pkin(3) * t239;
t192 = -pkin(7) * t239 - t205;
t191 = -t419 * t212 + t423 * t213;
t190 = t423 * t212 + t419 * t213;
t189 = -pkin(6) * t243 - t417 * t217 + t421 * t242;
t188 = -pkin(3) * t302 + pkin(7) * t241 + t206;
t187 = -pkin(6) * t230 - t417 * t214 + t421 * t235;
t186 = t421 * t206 - t456;
t185 = t417 * t206 + t446;
t184 = t422 * t186 + t418 * t274;
t183 = t418 * t186 - t422 * t274;
t182 = -t419 * t200 + t423 * t201;
t181 = t423 * t200 + t419 * t201;
t180 = -pkin(2) * t185 - pkin(3) * t205;
t179 = -pkin(5) * t215 + t422 * t189 - t418 * t199;
t178 = -pkin(5) * t212 + t422 * t187 - t418 * t197;
t177 = -pkin(1) * t243 + pkin(5) * t216 + t418 * t189 + t422 * t199;
t176 = -pkin(6) * t207 - t417 * t188 + t421 * t192;
t175 = -pkin(6) * t185 - pkin(7) * t446 - t417 * t198;
t174 = -pkin(1) * t230 + pkin(5) * t213 + t418 * t187 + t422 * t197;
t173 = -t419 * t183 + t423 * t184;
t172 = t423 * t183 + t419 * t184;
t171 = -pkin(5) * t200 + t422 * t176 - t418 * t193;
t170 = -pkin(1) * t207 + pkin(5) * t201 + t418 * t176 + t422 * t193;
t169 = -pkin(5) * t183 + t422 * t175 - t418 * t180;
t168 = -pkin(1) * t185 + pkin(5) * t184 + t418 * t175 + t422 * t180;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t390, -t391, 0, t351, 0, 0, 0, 0, 0, 0, -t428, t334, 0, t252, 0, 0, 0, 0, 0, 0, t258, t259, t282, t211, 0, 0, 0, 0, 0, 0, t191, t195, t182, t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t391, -t390, 0, t350, 0, 0, 0, 0, 0, 0, -t334, -t428, 0, -t473, 0, 0, 0, 0, 0, 0, t256, t257, t281, t210, 0, 0, 0, 0, 0, 0, t190, t194, t181, t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t339, t340, 0, -t261, 0, 0, 0, 0, 0, 0, t230, t243, t207, t185; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t391, 0, -t390, 0, t429, -t369, -t350, -pkin(4) * t350, 0, 0, -t334, 0, -t428, 0, t477, t476, t473, pkin(4) * t473 + pkin(5) * t442 - t419 * t293, -t419 * t315 + t423 * t317, -t419 * t300 + t423 * t301, -t419 * t322 + t423 * t324, -t419 * t314 + t423 * t316, -t419 * t321 + t423 * t323, -t419 * t348 + t423 * t349, -pkin(4) * t256 - t419 * t222 + t423 * t228, -pkin(4) * t257 - t419 * t223 + t423 * t229, -pkin(4) * t281 - t419 * t249 + t423 * t250, -pkin(4) * t210 - t419 * t196 + t423 * t204, -t419 * t225 + t423 * t227, -t419 * t202 + t423 * t203, -t419 * t218 + t423 * t220, -t419 * t224 + t423 * t226, -t419 * t219 + t423 * t221, -t419 * t253 + t423 * t254, -pkin(4) * t190 - t419 * t174 + t423 * t178, -pkin(4) * t194 - t419 * t177 + t423 * t179, -pkin(4) * t181 - t419 * t170 + t423 * t171, -pkin(4) * t172 - t419 * t168 + t423 * t169; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t390, 0, t391, 0, t369, t429, t351, pkin(4) * t351, 0, 0, t428, 0, -t334, 0, -t476, t477, t252, pkin(4) * t252 + pkin(5) * t449 + t423 * t293, t423 * t315 + t419 * t317, t423 * t300 + t419 * t301, t423 * t322 + t419 * t324, t423 * t314 + t419 * t316, t423 * t321 + t419 * t323, t423 * t348 + t419 * t349, pkin(4) * t258 + t423 * t222 + t419 * t228, pkin(4) * t259 + t423 * t223 + t419 * t229, pkin(4) * t282 + t423 * t249 + t419 * t250, pkin(4) * t211 + t423 * t196 + t419 * t204, t423 * t225 + t419 * t227, t423 * t202 + t419 * t203, t423 * t218 + t419 * t220, t423 * t224 + t419 * t226, t423 * t219 + t419 * t221, t423 * t253 + t419 * t254, pkin(4) * t191 + t423 * t174 + t419 * t178, pkin(4) * t195 + t423 * t177 + t419 * t179, pkin(4) * t182 + t423 * t170 + t419 * t171, pkin(4) * t173 + t423 * t168 + t419 * t169; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t398, t399, 0, 0, 0, 0, 0, 0, 0, t411, -pkin(1) * t382 - t337, -pkin(1) * t379 - t338, 0, -pkin(1) * t296, (t371 + t434) * t417, t421 * t370 + t417 * t373, t421 * t392 + t453, (t372 - t435) * t421, t417 * t394 + t444, 0, pkin(1) * t306 + pkin(2) * t373 + pkin(6) * t342 - t445, pkin(1) * t307 - pkin(2) * t370 + pkin(6) * t344 + t455, pkin(1) * t330 + pkin(2) * t383 + pkin(6) * t377 + t262, pkin(1) * t247 - pkin(2) * t325 + pkin(6) * t262, t421 * t270 + t417 * t271, t421 * t238 + t417 * t240, t421 * t289 + t417 * t291, t421 * t268 + t417 * t269, t421 * t290 + t417 * t292, t421 * t304 + t417 * t305, pkin(1) * t212 - pkin(2) * t275 + pkin(6) * t231 + t421 * t214 + t417 * t235, pkin(1) * t215 - pkin(2) * t466 + pkin(6) * t244 + t421 * t217 + t417 * t242, pkin(1) * t200 - pkin(2) * t302 + pkin(6) * t209 + t421 * t188 + t417 * t192, pkin(1) * t183 - pkin(2) * t274 + pkin(6) * t186 - pkin(7) * t456 + t421 * t198;];
tauB_reg = t1;
