% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRPRR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:16
% EndTime: 2019-12-05 15:45:25
% DurationCPUTime: 5.79s
% Computational Cost: add. (19103->377), mult. (27918->562), div. (0->0), fcn. (19096->10), ass. (0->257)
t369 = qJD(2) + qJD(4);
t367 = t369 ^ 2;
t381 = cos(qJ(4));
t368 = qJDD(2) + qJDD(4);
t378 = sin(qJ(4));
t417 = t378 * t368;
t326 = t367 * t381 + t417;
t409 = t381 * t368;
t329 = t367 * t378 - t409;
t373 = sin(pkin(9));
t375 = cos(pkin(9));
t270 = t326 * t375 - t329 * t373;
t274 = t326 * t373 + t329 * t375;
t379 = sin(qJ(2));
t382 = cos(qJ(2));
t236 = t270 * t379 + t274 * t382;
t374 = sin(pkin(8));
t459 = t374 * t236;
t376 = cos(pkin(8));
t458 = t376 * t236;
t349 = g(1) * t374 - g(2) * t376;
t344 = -qJDD(3) + t349;
t285 = pkin(6) * t326 - t344 * t381;
t447 = pkin(6) * t329 - t344 * t378;
t218 = qJ(3) * t270 + t285 * t375 - t373 * t447;
t453 = qJ(3) * t274 + t285 * t373 + t375 * t447;
t159 = pkin(5) * t236 + t218 * t379 + t382 * t453;
t445 = t270 * t382 - t274 * t379;
t160 = pkin(5) * t445 + t218 * t382 - t379 * t453;
t350 = g(1) * t376 + g(2) * t374;
t372 = g(3) - qJDD(1);
t316 = -t350 * t379 + t372 * t382;
t386 = qJDD(2) * pkin(2) - t316;
t317 = -t350 * t382 - t372 * t379;
t434 = qJD(2) ^ 2;
t388 = -pkin(2) * t434 + t317;
t252 = t373 * t386 + t375 * t388;
t250 = -pkin(3) * t434 + t252;
t385 = -t373 * t388 + t375 * t386;
t384 = qJDD(2) * pkin(3) + t385;
t205 = t250 * t378 - t381 * t384;
t206 = t250 * t381 + t378 * t384;
t395 = t205 * t378 + t206 * t381;
t176 = t205 * t381 - t206 * t378;
t425 = t375 * t176;
t147 = -t373 * t395 + t425;
t431 = t373 * t176;
t448 = t375 * t395 + t431;
t133 = t147 * t379 + t382 * t448;
t132 = t147 * t382 - t379 * t448;
t345 = qJDD(2) * t373 + t375 * t434;
t346 = qJDD(2) * t375 - t373 * t434;
t392 = -t345 * t379 + t346 * t382;
t451 = t374 * t392;
t449 = t376 * t392;
t394 = t252 * t375 - t373 * t385;
t209 = -t252 * t373 - t375 * t385;
t408 = t382 * t209;
t179 = -t379 * t394 + t408;
t416 = t379 * t209;
t180 = t382 * t394 + t416;
t300 = qJ(3) * t345 - t344 * t375;
t387 = -qJ(3) * t346 - t344 * t373;
t225 = -pkin(5) * t392 + t300 * t379 + t382 * t387;
t436 = t345 * t382 + t346 * t379;
t226 = pkin(5) * t436 + t300 * t382 - t379 * t387;
t338 = t376 * t349;
t307 = -t350 * t374 + t338;
t433 = pkin(2) * t344;
t377 = sin(qJ(5));
t370 = t377 ^ 2;
t432 = t370 * t367;
t430 = t374 * t344;
t348 = qJDD(2) * t382 - t379 * t434;
t429 = t374 * t348;
t428 = t374 * t349;
t426 = t374 * t372;
t424 = t376 * t348;
t423 = t376 * t372;
t198 = -pkin(4) * t368 - pkin(7) * t367 + t205;
t422 = t377 * t198;
t380 = cos(qJ(5));
t355 = t377 * t367 * t380;
t342 = qJDD(5) + t355;
t421 = t377 * t342;
t343 = qJDD(5) - t355;
t420 = t377 * t343;
t419 = t377 * t368;
t413 = t380 * t198;
t412 = t380 * t342;
t411 = t380 * t343;
t357 = t380 * t368;
t199 = -pkin(4) * t367 + pkin(7) * t368 + t206;
t196 = t199 * t380 - t344 * t377;
t371 = t380 ^ 2;
t407 = t370 + t371;
t406 = qJD(5) * t369;
t405 = t376 * qJDD(2);
t404 = t377 * t406;
t403 = t380 * t406;
t195 = t199 * t377 + t344 * t380;
t163 = t195 * t377 + t196 * t380;
t324 = t407 * t368;
t358 = t371 * t367;
t332 = t358 + t432;
t276 = t324 * t378 + t332 * t381;
t277 = t324 * t381 - t332 * t378;
t238 = t276 * t375 + t277 * t373;
t239 = -t276 * t373 + t277 * t375;
t191 = t238 * t382 + t239 * t379;
t192 = -t238 * t379 + t239 * t382;
t402 = -pkin(1) * t191 - pkin(2) * t238 - pkin(3) * t276 - pkin(4) * t332 - pkin(7) * t324 + qJ(1) * t192 - t163;
t401 = pkin(1) * t445 + pkin(2) * t270 + pkin(3) * t326 + qJ(1) * t236 + t206;
t400 = pkin(1) * t236 + pkin(2) * t274 + pkin(3) * t329 - qJ(1) * t445 + t205;
t399 = pkin(1) * t436 + pkin(2) * t345 - qJ(1) * t392 + t252;
t398 = -pkin(1) * t392 - pkin(2) * t346 - qJ(1) * t436 - t385;
t347 = qJDD(2) * t379 + t382 * t434;
t397 = -pkin(1) * t347 + qJ(1) * t348 - t317;
t396 = pkin(1) * t348 + qJ(1) * t347 - t316;
t259 = t316 * t379 + t317 * t382;
t308 = -t350 * t376 - t428;
t390 = t378 * t355;
t389 = t381 * t355;
t314 = pkin(5) * t347 - t349 * t382;
t313 = -pkin(5) * t348 - t349 * t379;
t162 = t195 * t380 - t196 * t377;
t258 = t316 * t382 - t317 * t379;
t383 = qJD(5) ^ 2;
t362 = t374 * qJDD(2);
t354 = -t358 - t383;
t353 = t358 - t383;
t352 = -t383 - t432;
t351 = t383 - t432;
t341 = pkin(1) * t344;
t337 = t376 * t347;
t336 = t374 * t347;
t333 = t358 - t432;
t325 = t376 * t344;
t323 = t357 - 0.2e1 * t404;
t322 = t357 - t404;
t321 = t403 + t419;
t320 = 0.2e1 * t403 + t419;
t319 = t407 * t406;
t312 = qJDD(5) * t378 + t319 * t381;
t311 = -qJDD(5) * t381 + t319 * t378;
t306 = t321 * t380 - t370 * t406;
t305 = -t322 * t377 - t371 * t406;
t304 = -t352 * t377 - t411;
t303 = -t351 * t377 + t412;
t302 = t354 * t380 - t421;
t301 = t353 * t380 - t420;
t298 = t352 * t380 - t420;
t297 = -t351 * t380 - t421;
t296 = t354 * t377 + t412;
t295 = -t353 * t377 - t411;
t287 = (-t321 - t403) * t377;
t286 = (-t322 + t404) * t380;
t279 = t376 * t436;
t278 = t374 * t436;
t269 = -t320 * t377 + t323 * t380;
t268 = -t320 * t380 - t323 * t377;
t267 = t303 * t381 + t377 * t417;
t266 = t301 * t381 + t357 * t378;
t265 = t303 * t378 - t377 * t409;
t264 = t301 * t378 - t380 * t409;
t263 = t306 * t381 - t390;
t262 = t305 * t381 + t390;
t261 = t306 * t378 + t389;
t260 = t305 * t378 - t389;
t256 = t304 * t381 + t320 * t378;
t255 = t302 * t381 - t323 * t378;
t254 = t304 * t378 - t320 * t381;
t253 = t302 * t378 + t323 * t381;
t248 = -t311 * t373 + t312 * t375;
t247 = t311 * t375 + t312 * t373;
t243 = t269 * t381 - t333 * t378;
t242 = t269 * t378 + t333 * t381;
t241 = t259 * t376 - t428;
t240 = t259 * t374 + t338;
t232 = t376 * t445;
t231 = t374 * t445;
t230 = -t265 * t373 + t267 * t375;
t229 = -t264 * t373 + t266 * t375;
t228 = t265 * t375 + t267 * t373;
t227 = t264 * t375 + t266 * t373;
t224 = -t261 * t373 + t263 * t375;
t223 = -t260 * t373 + t262 * t375;
t222 = t261 * t375 + t263 * t373;
t221 = t260 * t375 + t262 * t373;
t214 = -t254 * t373 + t256 * t375;
t213 = -t253 * t373 + t255 * t375;
t212 = t254 * t375 + t256 * t373;
t211 = t253 * t375 + t255 * t373;
t203 = -t247 * t379 + t248 * t382;
t202 = -t242 * t373 + t243 * t375;
t201 = t242 * t375 + t243 * t373;
t200 = qJ(3) * t394 + t433;
t194 = -pkin(7) * t298 + t413;
t193 = -pkin(7) * t296 + t422;
t190 = -pkin(4) * t298 + t196;
t189 = -pkin(4) * t296 + t195;
t188 = -t228 * t379 + t230 * t382;
t187 = -t227 * t379 + t229 * t382;
t186 = -t222 * t379 + t224 * t382;
t185 = -t221 * t379 + t223 * t382;
t184 = -t212 * t379 + t214 * t382;
t183 = -t211 * t379 + t213 * t382;
t182 = t212 * t382 + t214 * t379;
t181 = t211 * t382 + t213 * t379;
t173 = t180 * t376 - t430;
t172 = t180 * t374 + t325;
t169 = -t201 * t379 + t202 * t382;
t168 = pkin(3) * t344 + pkin(6) * t395;
t167 = t184 * t376 + t298 * t374;
t166 = t183 * t376 + t296 * t374;
t165 = t184 * t374 - t298 * t376;
t164 = t183 * t374 - t296 * t376;
t158 = -pkin(6) * t276 + t162 * t381;
t157 = pkin(6) * t277 + t162 * t378;
t156 = pkin(1) * t179 + pkin(2) * t209;
t155 = -pkin(6) * t254 - t190 * t378 + t194 * t381;
t154 = -pkin(6) * t253 - t189 * t378 + t193 * t381;
t153 = -pkin(3) * t298 + pkin(6) * t256 + t190 * t381 + t194 * t378;
t152 = -pkin(3) * t296 + pkin(6) * t255 + t189 * t381 + t193 * t378;
t151 = t163 * t381 + t198 * t378;
t150 = t163 * t378 - t198 * t381;
t149 = pkin(5) * t179 + qJ(3) * t408 - t200 * t379;
t144 = -pkin(1) * t182 - pkin(2) * t212 - pkin(3) * t254 + pkin(4) * t320 - pkin(7) * t304 - t422;
t143 = -pkin(1) * t181 - pkin(2) * t211 - pkin(3) * t253 - pkin(4) * t323 - pkin(7) * t302 + t413;
t141 = -qJ(3) * t238 - t157 * t373 + t158 * t375;
t140 = qJ(3) * t239 + t157 * t375 + t158 * t373;
t139 = -t150 * t373 + t151 * t375;
t138 = t150 * t375 + t151 * t373;
t137 = -qJ(3) * t212 - t153 * t373 + t155 * t375;
t136 = -qJ(3) * t211 - t152 * t373 + t154 * t375;
t135 = -pkin(2) * t298 + qJ(3) * t214 + t153 * t375 + t155 * t373;
t134 = -pkin(2) * t296 + qJ(3) * t213 + t152 * t375 + t154 * t373;
t130 = pkin(6) * t425 + qJ(3) * t147 - t168 * t373;
t129 = t133 * t376 - t430;
t128 = t133 * t374 + t325;
t127 = pkin(6) * t431 + qJ(3) * t448 + t168 * t375 + t433;
t126 = -pkin(6) * t150 - (pkin(4) * t378 - pkin(7) * t381) * t162;
t125 = pkin(6) * t151 - (-pkin(4) * t381 - pkin(7) * t378 - pkin(3)) * t162;
t124 = -pkin(5) * t191 - t140 * t379 + t141 * t382;
t123 = pkin(1) * t132 + pkin(2) * t147 + pkin(3) * t176;
t122 = -t138 * t379 + t139 * t382;
t121 = t138 * t382 + t139 * t379;
t120 = t122 * t376 - t162 * t374;
t119 = t122 * t374 + t162 * t376;
t118 = -pkin(5) * t182 - t135 * t379 + t137 * t382;
t117 = -pkin(5) * t181 - t134 * t379 + t136 * t382;
t116 = pkin(5) * t132 - t127 * t379 + t130 * t382;
t115 = -qJ(3) * t138 - t125 * t373 + t126 * t375;
t114 = -pkin(1) * t121 - pkin(2) * t138 - pkin(3) * t150 + pkin(4) * t198 - pkin(7) * t163;
t113 = pkin(2) * t162 + qJ(3) * t139 + t125 * t375 + t126 * t373;
t112 = -pkin(5) * t121 - t113 * t379 + t115 * t382;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t308, 0, 0, 0, 0, 0, 0, -t337, -t424, 0, t241, 0, 0, 0, 0, 0, 0, -t279, -t449, 0, t173, 0, 0, 0, 0, 0, 0, -t232, t458, 0, t129, 0, 0, 0, 0, 0, 0, t166, t167, t376 * t192, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, 0, 0, 0, 0, 0, 0, -t336, -t429, 0, t240, 0, 0, 0, 0, 0, 0, -t278, -t451, 0, t172, 0, 0, 0, 0, 0, 0, -t231, t459, 0, t128, 0, 0, 0, 0, 0, 0, t164, t165, t374 * t192, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t372, 0, 0, 0, 0, 0, 0, t348, -t347, 0, -t258, 0, 0, 0, 0, 0, 0, t392, -t436, 0, -t179, 0, 0, 0, 0, 0, 0, -t236, -t445, 0, -t132, 0, 0, 0, 0, 0, 0, t181, t182, t191, t121; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t426, -t423, -t307, -qJ(1) * t307, 0, 0, t424, 0, -t337, t362, t376 * t313 + t374 * t396, t376 * t314 + t374 * t397, t376 * t258, -qJ(1) * t240 - (pkin(1) * t374 - pkin(5) * t376) * t258, 0, 0, t449, 0, -t279, t362, t376 * t225 - t374 * t398, t376 * t226 - t374 * t399, t376 * t179, -qJ(1) * t172 + t149 * t376 - t156 * t374, 0, 0, -t458, 0, -t232, t374 * t368, t376 * t159 - t374 * t400, t376 * t160 - t374 * t401, t376 * t132, -qJ(1) * t128 + t116 * t376 - t123 * t374, t186 * t376 - t287 * t374, t169 * t376 - t268 * t374, t188 * t376 - t297 * t374, t185 * t376 - t286 * t374, t187 * t376 - t295 * t374, t376 * t203, -qJ(1) * t164 + t117 * t376 - t143 * t374, -qJ(1) * t165 + t118 * t376 - t144 * t374, t376 * t124 - t374 * t402, -qJ(1) * t119 + t112 * t376 - t114 * t374; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t423, -t426, t308, qJ(1) * t308, 0, 0, t429, 0, -t336, -t405, t374 * t313 - t376 * t396, t374 * t314 - t376 * t397, t374 * t258, qJ(1) * t241 - (-pkin(1) * t376 - pkin(5) * t374) * t258, 0, 0, t451, 0, -t278, -t405, t374 * t225 + t376 * t398, t374 * t226 + t376 * t399, t374 * t179, qJ(1) * t173 + t149 * t374 + t156 * t376, 0, 0, -t459, 0, -t231, -t376 * t368, t374 * t159 + t376 * t400, t374 * t160 + t376 * t401, t374 * t132, qJ(1) * t129 + t116 * t374 + t123 * t376, t186 * t374 + t287 * t376, t169 * t374 + t268 * t376, t188 * t374 + t297 * t376, t185 * t374 + t286 * t376, t187 * t374 + t295 * t376, t374 * t203, qJ(1) * t166 + t117 * t374 + t143 * t376, qJ(1) * t167 + t118 * t374 + t144 * t376, t374 * t124 + t376 * t402, qJ(1) * t120 + t112 * t374 + t114 * t376; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t349, t350, 0, 0, 0, 0, t347, 0, t348, 0, -t314, t313, t259, pkin(1) * t349 + pkin(5) * t259, 0, 0, t436, 0, t392, 0, -t226, t225, t180, pkin(5) * t180 + qJ(3) * t416 + t200 * t382 + t341, 0, 0, t445, 0, -t236, 0, -t160, t159, t133, pkin(5) * t133 + t127 * t382 + t130 * t379 + t341, t222 * t382 + t224 * t379, t201 * t382 + t202 * t379, t228 * t382 + t230 * t379, t221 * t382 + t223 * t379, t227 * t382 + t229 * t379, t247 * t382 + t248 * t379, -pkin(1) * t296 + pkin(5) * t183 + t134 * t382 + t136 * t379, -pkin(1) * t298 + pkin(5) * t184 + t135 * t382 + t137 * t379, pkin(5) * t192 + t140 * t382 + t141 * t379, pkin(1) * t162 + pkin(5) * t122 + t113 * t382 + t115 * t379;];
tauB_reg = t1;
