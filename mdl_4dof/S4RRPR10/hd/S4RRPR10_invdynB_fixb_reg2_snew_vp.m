% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPR10_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:58
% EndTime: 2019-12-31 17:12:03
% DurationCPUTime: 3.48s
% Computational Cost: add. (5746->340), mult. (12299->450), div. (0->0), fcn. (6859->6), ass. (0->247)
t380 = sin(qJ(2));
t375 = t380 ^ 2;
t386 = qJD(1) ^ 2;
t367 = t375 * t386;
t385 = qJD(2) ^ 2;
t355 = -t367 - t385;
t383 = cos(qJ(2));
t411 = t380 * t383 * t386;
t350 = qJDD(2) - t411;
t426 = t383 * t350;
t303 = t355 * t380 + t426;
t420 = qJD(1) * qJD(2);
t404 = t383 * t420;
t419 = qJDD(1) * t380;
t339 = 0.2e1 * t404 + t419;
t381 = sin(qJ(1));
t384 = cos(qJ(1));
t264 = t303 * t381 + t339 * t384;
t474 = pkin(4) * t264;
t268 = t303 * t384 - t339 * t381;
t473 = pkin(4) * t268;
t376 = t383 ^ 2;
t369 = t376 * t386;
t357 = -t369 - t385;
t349 = qJDD(2) + t411;
t436 = t349 * t380;
t300 = -t357 * t383 + t436;
t405 = t380 * t420;
t417 = qJDD(1) * t383;
t342 = -0.2e1 * t405 + t417;
t263 = t300 * t381 - t342 * t384;
t472 = pkin(4) * t263;
t267 = t300 * t384 + t342 * t381;
t471 = pkin(4) * t267;
t427 = t383 * t349;
t294 = t357 * t380 + t427;
t470 = pkin(1) * t294;
t469 = pkin(5) * t294;
t468 = pkin(5) * t303;
t354 = -t367 + t385;
t302 = -t354 * t380 + t427;
t416 = qJDD(1) * t384;
t467 = t302 * t381 - t380 * t416;
t418 = qJDD(1) * t381;
t466 = t302 * t384 + t380 * t418;
t465 = 2 * qJD(3);
t435 = t350 * t380;
t296 = -t355 * t383 + t435;
t464 = pkin(1) * t296;
t463 = pkin(5) * t296;
t462 = pkin(5) * t300;
t379 = sin(qJ(4));
t382 = cos(qJ(4));
t421 = qJD(1) * t383;
t332 = qJD(2) * t379 + t382 * t421;
t334 = qJD(2) * t382 - t379 * t421;
t289 = t334 * t332;
t340 = t404 + t419;
t325 = qJDD(4) + t340;
t456 = -t289 + t325;
t461 = t379 * t456;
t460 = t382 * t456;
t356 = t369 - t385;
t301 = -t356 * t383 + t435;
t459 = t301 * t381 + t383 * t416;
t458 = t301 * t384 - t381 * t417;
t397 = t340 + t404;
t457 = t397 * qJ(3);
t422 = qJD(1) * t380;
t455 = -pkin(2) * t405 + t422 * t465;
t454 = t356 * t380 + t426;
t323 = t332 ^ 2;
t324 = t334 ^ 2;
t361 = qJD(4) + t422;
t358 = t361 ^ 2;
t453 = pkin(2) + pkin(6);
t452 = pkin(2) * t383;
t423 = t375 + t376;
t344 = t423 * qJDD(1);
t347 = t367 + t369;
t285 = t344 * t381 + t347 * t384;
t451 = pkin(4) * t285;
t450 = pkin(5) * t386;
t449 = qJ(3) * t380;
t341 = -t405 + t417;
t351 = pkin(3) * t422 - qJD(2) * pkin(6);
t399 = -t449 - t452;
t337 = t399 * qJD(1);
t353 = g(1) * t384 + g(2) * t381;
t319 = -pkin(1) * t386 + qJDD(1) * pkin(5) - t353;
t424 = g(3) * t380 - t319 * t383;
t392 = pkin(2) * t385 - t337 * t421 + t424;
t415 = qJDD(2) * qJ(3);
t222 = t415 + t341 * pkin(3) - pkin(6) * t369 + (t465 + t351) * qJD(2) - t392;
t448 = t222 * t379;
t447 = t222 * t382;
t273 = t289 + t325;
t446 = t273 * t382;
t352 = g(1) * t381 - g(2) * t384;
t394 = -qJDD(1) * pkin(1) - t352;
t318 = -t394 + t450;
t445 = t318 * t380;
t444 = t318 * t383;
t443 = t332 * t361;
t442 = t339 * t380;
t438 = t342 * t383;
t430 = t361 * t379;
t429 = t361 * t382;
t428 = t379 * t273;
t306 = g(3) * t383 + t319 * t380;
t425 = pkin(1) * t347 + pkin(5) * t344;
t414 = -t324 - t358;
t413 = t380 * t289;
t412 = t383 * t289;
t279 = -qJD(4) * t332 + qJDD(2) * t382 - t341 * t379;
t251 = t306 * t380 - t383 * t424;
t291 = -t352 * t381 - t353 * t384;
t403 = qJDD(2) * t379 + t341 * t382;
t402 = t381 * t411;
t401 = t384 * t411;
t346 = -t381 * t386 + t416;
t400 = -pkin(4) * t346 - g(3) * t381;
t398 = pkin(2) * t380 - qJ(3) * t383;
t388 = t394 - t455 - t457;
t216 = -t351 * t422 + (-pkin(3) * t376 - pkin(5)) * t386 - t453 * t341 + t388;
t393 = -qJDD(2) * pkin(2) - qJ(3) * t385 + t337 * t422 + qJDD(3) + t306;
t223 = -t349 * pkin(6) + (t340 - t404) * pkin(3) + t393;
t193 = t216 * t379 - t223 * t382;
t194 = t216 * t382 + t223 * t379;
t177 = -t193 * t382 + t194 * t379;
t178 = t379 * t193 + t194 * t382;
t250 = t306 * t383 + t380 * t424;
t396 = t354 * t383 + t436;
t290 = t352 * t384 - t353 * t381;
t395 = t279 - t443;
t391 = (-qJD(4) + t361) * t334 - t403;
t390 = -0.2e1 * qJD(3) * qJD(2) + t392;
t252 = -t390 + t415;
t387 = pkin(2) * t341 + t318 + t455;
t348 = -t367 + t369;
t345 = t384 * t386 + t418;
t335 = t398 * qJDD(1);
t327 = t423 * t420;
t317 = -pkin(4) * t345 + g(3) * t384;
t313 = -t324 + t358;
t312 = t323 - t358;
t311 = qJDD(2) * t381 + t327 * t384;
t310 = t340 * t383 - t375 * t420;
t309 = -qJDD(2) * t384 + t327 * t381;
t308 = -t341 * t380 - t376 * t420;
t293 = t397 * t380;
t292 = (t341 - t405) * t383;
t287 = t324 - t323;
t286 = t344 * t384 - t347 * t381;
t284 = pkin(4) * t286;
t283 = t438 - t442;
t282 = t339 * t383 + t342 * t380;
t280 = -t358 - t323;
t278 = -qJD(4) * t334 - t403;
t277 = t310 * t384 - t402;
t276 = t308 * t384 + t402;
t275 = t310 * t381 + t401;
t274 = t308 * t381 - t401;
t271 = -t323 - t324;
t261 = -t444 + t463;
t260 = -t445 - t469;
t259 = (t332 * t382 - t334 * t379) * t361;
t258 = (t332 * t379 + t334 * t382) * t361;
t257 = t283 * t384 - t348 * t381;
t256 = t283 * t381 + t348 * t384;
t254 = -t424 + t464;
t253 = t306 - t470;
t248 = t279 + t443;
t243 = (qJD(4) + t361) * t334 + t403;
t241 = qJ(3) * t347 + t393;
t240 = pkin(2) * t347 + t252;
t239 = -t279 * t379 - t334 * t429;
t238 = t278 * t379 - t332 * t429;
t237 = -t279 * t382 + t334 * t430;
t236 = -t278 * t382 - t332 * t430;
t235 = t387 + t457;
t234 = -t258 * t380 + t325 * t383;
t233 = t313 * t379 - t460;
t232 = -t312 * t379 - t446;
t231 = -t312 * t382 + t428;
t230 = -t313 * t382 - t461;
t229 = -t379 * t414 - t446;
t228 = t382 * t414 - t428;
t227 = -t450 + (-t341 - t342) * pkin(2) + t388;
t226 = (t339 + t397) * qJ(3) + t387;
t225 = t251 * t384 - t318 * t381;
t224 = t251 * t381 + t318 * t384;
t221 = t280 * t382 - t461;
t220 = t280 * t379 + t460;
t219 = pkin(2) * t349 + qJ(3) * t357 - t393 + t470;
t218 = -t239 * t380 + t412;
t217 = -t236 * t380 - t412;
t215 = -t464 + pkin(2) * t355 + (-qJDD(2) - t350) * qJ(3) + t390;
t214 = t252 * t383 + t380 * t393;
t213 = t252 * t380 - t383 * t393;
t212 = t248 * t379 + t382 * t391;
t211 = t243 * t382 + t379 * t395;
t210 = -t248 * t382 + t379 * t391;
t209 = t243 * t379 - t382 * t395;
t208 = -pkin(2) * t442 + t226 * t383 - t463;
t207 = -qJ(3) * t438 - t227 * t380 + t469;
t206 = -t240 * t380 + t241 * t383;
t205 = -t230 * t380 + t248 * t383;
t204 = -t232 * t380 + t383 * t391;
t203 = t228 * t380 + t383 * t395;
t202 = -t228 * t383 + t380 * t395;
t201 = t220 * t380 + t243 * t383;
t200 = -t220 * t383 + t243 * t380;
t199 = -t209 * t380 + t287 * t383;
t198 = t210 * t380 + t271 * t383;
t197 = -t210 * t383 + t271 * t380;
t196 = t214 * t384 - t235 * t381;
t195 = t214 * t381 + t235 * t384;
t191 = -pkin(1) * t213 + pkin(2) * t393 - qJ(3) * t252;
t190 = t203 * t384 + t229 * t381;
t189 = t203 * t381 - t229 * t384;
t188 = t201 * t384 + t221 * t381;
t187 = t201 * t381 - t221 * t384;
t186 = pkin(3) * t210 - qJ(3) * t212;
t185 = -pkin(5) * t213 - t235 * t398;
t184 = t198 * t384 + t212 * t381;
t183 = t198 * t381 - t212 * t384;
t182 = pkin(3) * t395 - t229 * t453 - t448;
t181 = pkin(3) * t243 - t221 * t453 + t447;
t180 = pkin(3) * t228 - qJ(3) * t229 - t194;
t179 = pkin(3) * t220 - qJ(3) * t221 - t193;
t176 = -pkin(1) * t202 - qJ(3) * t395 + t228 * t453 - t447;
t175 = -pkin(1) * t200 - qJ(3) * t243 + t220 * t453 - t448;
t174 = t177 * t380 + t222 * t383;
t173 = -t177 * t383 + t222 * t380;
t172 = pkin(3) * t271 - t212 * t453 - t178;
t171 = pkin(3) * t177 - qJ(3) * t178;
t170 = -pkin(5) * t202 + t180 * t383 - t182 * t380;
t169 = -pkin(5) * t200 + t179 * t383 - t181 * t380;
t168 = -pkin(1) * t197 - qJ(3) * t271 + t210 * t453 + t177;
t167 = pkin(3) * t222 - t178 * t453;
t166 = t174 * t384 + t178 * t381;
t165 = t174 * t381 - t178 * t384;
t164 = -pkin(5) * t197 - t172 * t380 + t186 * t383;
t163 = -pkin(1) * t173 - qJ(3) * t222 + t177 * t453;
t162 = -pkin(5) * t173 - t167 * t380 + t171 * t383;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t345, -t346, 0, t291, 0, 0, 0, 0, 0, 0, -t267, -t268, t286, t225, 0, 0, 0, 0, 0, 0, t286, t267, t268, t196, 0, 0, 0, 0, 0, 0, t188, t190, t184, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t346, -t345, 0, t290, 0, 0, 0, 0, 0, 0, -t263, -t264, t285, t224, 0, 0, 0, 0, 0, 0, t285, t263, t264, t195, 0, 0, 0, 0, 0, 0, t187, t189, t183, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t294, -t296, 0, -t250, 0, 0, 0, 0, 0, 0, 0, -t294, t296, t213, 0, 0, 0, 0, 0, 0, t200, t202, t197, t173; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t346, 0, -t345, 0, t400, -t317, -t290, -pkin(4) * t290, t277, t257, t466, t276, -t458, t311, -t253 * t381 + t260 * t384 + t472, -t254 * t381 + t261 * t384 + t474, t250 * t384 - t451, -pkin(4) * t224 - (pkin(1) * t381 - pkin(5) * t384) * t250, t311, -t466, t458, t277, t257, t276, t206 * t384 - t335 * t381 - t451, t207 * t384 - t219 * t381 - t472, t208 * t384 - t215 * t381 - t474, -pkin(4) * t195 + t185 * t384 - t191 * t381, t218 * t384 - t237 * t381, t199 * t384 - t211 * t381, t205 * t384 - t233 * t381, t217 * t384 - t238 * t381, t204 * t384 - t231 * t381, t234 * t384 - t259 * t381, -pkin(4) * t187 + t169 * t384 - t175 * t381, -pkin(4) * t189 + t170 * t384 - t176 * t381, -pkin(4) * t183 + t164 * t384 - t168 * t381, -pkin(4) * t165 + t162 * t384 - t163 * t381; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t345, 0, t346, 0, t317, t400, t291, pkin(4) * t291, t275, t256, t467, t274, -t459, t309, t253 * t384 + t260 * t381 - t471, t254 * t384 + t261 * t381 - t473, t250 * t381 + t284, pkin(4) * t225 - (-pkin(1) * t384 - pkin(5) * t381) * t250, t309, -t467, t459, t275, t256, t274, t206 * t381 + t335 * t384 + t284, t207 * t381 + t219 * t384 + t471, t208 * t381 + t215 * t384 + t473, pkin(4) * t196 + t185 * t381 + t191 * t384, t218 * t381 + t237 * t384, t199 * t381 + t211 * t384, t205 * t381 + t233 * t384, t217 * t381 + t238 * t384, t204 * t381 + t231 * t384, t234 * t381 + t259 * t384, pkin(4) * t188 + t169 * t381 + t175 * t384, pkin(4) * t190 + t170 * t381 + t176 * t384, pkin(4) * t184 + t164 * t381 + t168 * t384, pkin(4) * t166 + t162 * t381 + t163 * t384; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t352, t353, 0, 0, t293, t282, t396, t292, t454, 0, pkin(1) * t342 + t444 - t462, -pkin(1) * t339 - t445 - t468, t251 + t425, pkin(1) * t318 + pkin(5) * t251, 0, -t396, -t454, t293, t282, t292, t240 * t383 + t241 * t380 + t425, t462 + t227 * t383 + (-pkin(1) - t449) * t342, t468 + t226 * t380 + (pkin(1) + t452) * t339, pkin(5) * t214 + (pkin(1) - t399) * t235, t239 * t383 + t413, t209 * t383 + t287 * t380, t230 * t383 + t248 * t380, t236 * t383 - t413, t232 * t383 + t380 * t391, t258 * t383 + t325 * t380, -pkin(1) * t221 + pkin(5) * t201 + t179 * t380 + t181 * t383, -pkin(1) * t229 + pkin(5) * t203 + t180 * t380 + t182 * t383, -pkin(1) * t212 + pkin(5) * t198 + t172 * t383 + t186 * t380, -pkin(1) * t178 + pkin(5) * t174 + t167 * t383 + t171 * t380;];
tauB_reg = t1;
