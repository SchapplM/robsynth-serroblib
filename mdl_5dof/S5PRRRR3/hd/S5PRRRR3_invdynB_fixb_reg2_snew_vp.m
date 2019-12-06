% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRRRR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRRRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:21
% EndTime: 2019-12-05 17:06:32
% DurationCPUTime: 7.44s
% Computational Cost: add. (27622->378), mult. (39390->548), div. (0->0), fcn. (28538->10), ass. (0->253)
t381 = qJD(2) + qJD(3);
t376 = qJD(4) + t381;
t374 = t376 ^ 2;
t392 = cos(qJ(4));
t380 = qJDD(2) + qJDD(3);
t375 = qJDD(4) + t380;
t388 = sin(qJ(4));
t424 = t388 * t375;
t339 = t392 * t374 + t424;
t417 = t392 * t375;
t342 = t388 * t374 - t417;
t389 = sin(qJ(3));
t393 = cos(qJ(3));
t283 = t393 * t339 - t389 * t342;
t384 = g(3) - qJDD(1);
t322 = pkin(7) * t339 - t392 * t384;
t456 = pkin(7) * t342 - t388 * t384;
t232 = pkin(6) * t283 + t393 * t322 - t389 * t456;
t288 = t389 * t339 + t393 * t342;
t390 = sin(qJ(2));
t394 = cos(qJ(2));
t238 = t394 * t283 - t390 * t288;
t468 = pkin(6) * t288 + t389 * t322 + t393 * t456;
t172 = pkin(5) * t238 + t394 * t232 - t390 * t468;
t243 = t390 * t283 + t394 * t288;
t385 = sin(pkin(9));
t386 = cos(pkin(9));
t197 = t385 * t238 + t386 * t243;
t483 = pkin(5) * t243 + t390 * t232 + t394 * t468;
t490 = qJ(1) * t197 + t385 * t172 + t386 * t483;
t470 = t386 * t238 - t385 * t243;
t489 = qJ(1) * t470 + t386 * t172 - t385 * t483;
t366 = t385 * g(1) - t386 * g(2);
t367 = t386 * g(1) + t385 * g(2);
t401 = t394 * t366 + t390 * t367;
t307 = qJDD(2) * pkin(2) + t401;
t412 = -t390 * t366 + t394 * t367;
t438 = qJD(2) ^ 2;
t308 = -pkin(2) * t438 - t412;
t262 = t389 * t307 + t393 * t308;
t379 = t381 ^ 2;
t256 = -t379 * pkin(3) + t262;
t397 = t393 * t307 - t389 * t308;
t396 = t380 * pkin(3) + t397;
t208 = t388 * t256 - t392 * t396;
t209 = t392 * t256 + t388 * t396;
t407 = t388 * t208 + t392 * t209;
t167 = t392 * t208 - t388 * t209;
t415 = t393 * t167;
t141 = -t389 * t407 + t415;
t422 = t389 * t167;
t460 = t393 * t407 + t422;
t128 = t394 * t141 - t390 * t460;
t476 = t390 * t141 + t394 * t460;
t116 = t385 * t128 + t386 * t476;
t486 = t386 * t128 - t385 * t476;
t346 = t393 * t379 + t389 * t380;
t349 = t389 * t379 - t393 * t380;
t290 = t394 * t346 - t390 * t349;
t326 = pkin(6) * t346 - t393 * t384;
t457 = pkin(6) * t349 - t389 * t384;
t248 = pkin(5) * t290 + t394 * t326 - t390 * t457;
t294 = t390 * t346 + t394 * t349;
t251 = t385 * t290 + t386 * t294;
t469 = pkin(5) * t294 + t390 * t326 + t394 * t457;
t482 = qJ(1) * t251 + t385 * t248 + t386 * t469;
t454 = t386 * t290 - t385 * t294;
t481 = qJ(1) * t454 + t386 * t248 - t385 * t469;
t406 = t393 * t262 - t389 * t397;
t215 = -t389 * t262 - t393 * t397;
t413 = t394 * t215;
t177 = -t390 * t406 + t413;
t420 = t390 * t215;
t461 = t394 * t406 + t420;
t146 = t385 * t177 + t386 * t461;
t475 = t386 * t177 - t385 * t461;
t405 = -t390 * t401 - t394 * t412;
t267 = t390 * t412 - t394 * t401;
t434 = t385 * t267;
t220 = t386 * t405 + t434;
t430 = t386 * t267;
t459 = -t385 * t405 + t430;
t363 = t390 * qJDD(2) + t394 * t438;
t364 = t394 * qJDD(2) - t390 * t438;
t312 = -t385 * t363 + t386 * t364;
t336 = pkin(5) * t363 - t394 * t384;
t398 = -pkin(5) * t364 - t390 * t384;
t458 = -qJ(1) * t312 + t385 * t336 + t386 * t398;
t440 = t386 * t363 + t385 * t364;
t452 = qJ(1) * t440 + t386 * t336 - t385 * t398;
t437 = pkin(1) * t384;
t436 = pkin(2) * t384;
t387 = sin(qJ(5));
t382 = t387 ^ 2;
t435 = t382 * t374;
t431 = t385 * t384;
t429 = t386 * t384;
t205 = -t375 * pkin(4) - t374 * pkin(8) + t208;
t428 = t387 * t205;
t391 = cos(qJ(5));
t362 = t387 * t374 * t391;
t353 = qJDD(5) + t362;
t427 = t387 * t353;
t354 = qJDD(5) - t362;
t426 = t387 * t354;
t425 = t387 * t375;
t419 = t391 * t205;
t418 = t391 * t354;
t368 = t391 * t375;
t206 = -t374 * pkin(4) + t375 * pkin(8) + t209;
t203 = t391 * t206 - t387 * t384;
t383 = t391 ^ 2;
t411 = t382 + t383;
t410 = qJD(5) * t376;
t409 = t387 * t410;
t408 = t391 * t410;
t202 = t387 * t206 + t391 * t384;
t163 = t387 * t202 + t391 * t203;
t316 = -t385 * t366 - t386 * t367;
t400 = t388 * t362;
t399 = t392 * t362;
t162 = t391 * t202 - t387 * t203;
t315 = t386 * t366 - t385 * t367;
t395 = qJD(5) ^ 2;
t369 = t383 * t374;
t361 = -t369 - t395;
t360 = t369 - t395;
t359 = -t395 - t435;
t358 = t395 - t435;
t345 = t391 * t353;
t344 = t369 - t435;
t343 = t369 + t435;
t337 = t411 * t375;
t334 = t368 - 0.2e1 * t409;
t333 = t368 - t409;
t330 = t408 + t425;
t329 = 0.2e1 * t408 + t425;
t328 = t411 * t410;
t310 = t388 * qJDD(5) + t392 * t328;
t309 = -t392 * qJDD(5) + t388 * t328;
t305 = -t387 * t359 - t418;
t304 = -t387 * t358 + t345;
t303 = t391 * t361 - t427;
t302 = t391 * t360 - t426;
t301 = t391 * t359 - t426;
t300 = t387 * t361 + t345;
t299 = t391 * t330 - t382 * t410;
t298 = -t387 * t333 - t383 * t410;
t286 = t392 * t337 - t388 * t343;
t282 = t388 * t337 + t392 * t343;
t281 = -t387 * t329 + t391 * t334;
t280 = t392 * t304 + t387 * t424;
t279 = t392 * t302 + t368 * t388;
t278 = t388 * t304 - t387 * t417;
t277 = t388 * t302 - t391 * t417;
t276 = t392 * t299 - t400;
t275 = t392 * t298 + t400;
t274 = t388 * t299 + t399;
t273 = t388 * t298 - t399;
t272 = t392 * t305 + t388 * t329;
t271 = t392 * t303 - t388 * t334;
t270 = t388 * t305 - t392 * t329;
t269 = t388 * t303 + t392 * t334;
t264 = -t389 * t309 + t393 * t310;
t263 = t393 * t309 + t389 * t310;
t260 = pkin(5) * t405 + t437;
t259 = t392 * t281 - t388 * t344;
t258 = t388 * t281 + t392 * t344;
t241 = -t389 * t282 + t393 * t286;
t237 = t393 * t282 + t389 * t286;
t236 = -t389 * t278 + t393 * t280;
t235 = -t389 * t277 + t393 * t279;
t234 = t393 * t278 + t389 * t280;
t233 = t393 * t277 + t389 * t279;
t228 = -t389 * t274 + t393 * t276;
t227 = -t389 * t273 + t393 * t275;
t226 = t393 * t274 + t389 * t276;
t225 = t393 * t273 + t389 * t275;
t224 = -t389 * t270 + t393 * t272;
t223 = -t389 * t269 + t393 * t271;
t222 = t393 * t270 + t389 * t272;
t221 = t393 * t269 + t389 * t271;
t218 = -t390 * t263 + t394 * t264;
t217 = t394 * t263 + t390 * t264;
t212 = -t389 * t258 + t393 * t259;
t211 = t393 * t258 + t389 * t259;
t210 = pkin(6) * t406 + t436;
t200 = -t390 * t237 + t394 * t241;
t199 = t394 * t237 + t390 * t241;
t194 = -pkin(8) * t301 + t419;
t193 = -pkin(8) * t300 + t428;
t192 = -t390 * t234 + t394 * t236;
t191 = -t390 * t233 + t394 * t235;
t190 = t394 * t234 + t390 * t236;
t189 = t394 * t233 + t390 * t235;
t188 = -pkin(4) * t301 + t203;
t187 = -pkin(4) * t300 + t202;
t186 = -t390 * t226 + t394 * t228;
t185 = -t390 * t225 + t394 * t227;
t184 = t394 * t226 + t390 * t228;
t183 = t394 * t225 + t390 * t227;
t182 = -t390 * t222 + t394 * t224;
t181 = -t390 * t221 + t394 * t223;
t180 = t394 * t222 + t390 * t224;
t179 = t394 * t221 + t390 * t223;
t174 = -t390 * t211 + t394 * t212;
t173 = t394 * t211 + t390 * t212;
t164 = pkin(3) * t384 + pkin(7) * t407;
t160 = -t385 * t199 + t386 * t200;
t159 = t386 * t199 + t385 * t200;
t158 = -pkin(7) * t282 + t392 * t162;
t157 = pkin(7) * t286 + t388 * t162;
t156 = -pkin(7) * t270 - t388 * t188 + t392 * t194;
t155 = -pkin(7) * t269 - t388 * t187 + t392 * t193;
t154 = -pkin(3) * t301 + pkin(7) * t272 + t392 * t188 + t388 * t194;
t153 = -pkin(3) * t300 + pkin(7) * t271 + t392 * t187 + t388 * t193;
t152 = -t385 * t180 + t386 * t182;
t151 = -t385 * t179 + t386 * t181;
t150 = t386 * t180 + t385 * t182;
t149 = t386 * t179 + t385 * t181;
t148 = t392 * t163 + t388 * t205;
t147 = t388 * t163 - t392 * t205;
t144 = pkin(5) * t177 + pkin(6) * t413 - t390 * t210;
t143 = pkin(5) * t461 + pkin(6) * t420 + t394 * t210 + t437;
t138 = -pkin(6) * t237 - t389 * t157 + t393 * t158;
t137 = pkin(6) * t241 + t393 * t157 + t389 * t158;
t136 = -t389 * t147 + t393 * t148;
t135 = t393 * t147 + t389 * t148;
t134 = -pkin(6) * t222 - t389 * t154 + t393 * t156;
t133 = -pkin(6) * t221 - t389 * t153 + t393 * t155;
t132 = -pkin(2) * t301 + pkin(6) * t224 + t393 * t154 + t389 * t156;
t131 = -pkin(2) * t300 + pkin(6) * t223 + t393 * t153 + t389 * t155;
t130 = -pkin(7) * t147 - (pkin(4) * t388 - pkin(8) * t392) * t162;
t125 = pkin(6) * t141 + pkin(7) * t415 - t389 * t164;
t124 = pkin(6) * t460 + pkin(7) * t422 + t393 * t164 + t436;
t123 = pkin(7) * t148 - (-pkin(4) * t392 - pkin(8) * t388 - pkin(3)) * t162;
t122 = -pkin(5) * t199 - t390 * t137 + t394 * t138;
t121 = pkin(5) * t200 + t394 * t137 + t390 * t138;
t120 = -t390 * t135 + t394 * t136;
t119 = t394 * t135 + t390 * t136;
t118 = -pkin(5) * t180 - t390 * t132 + t394 * t134;
t117 = -pkin(5) * t179 - t390 * t131 + t394 * t133;
t114 = -pkin(1) * t301 + pkin(5) * t182 + t394 * t132 + t390 * t134;
t113 = -pkin(1) * t300 + pkin(5) * t181 + t394 * t131 + t390 * t133;
t112 = pkin(5) * t128 - t390 * t124 + t394 * t125;
t111 = pkin(5) * t476 + t394 * t124 + t390 * t125 + t437;
t110 = -pkin(6) * t135 - t389 * t123 + t393 * t130;
t109 = -t385 * t119 + t386 * t120;
t108 = t386 * t119 + t385 * t120;
t107 = pkin(2) * t162 + pkin(6) * t136 + t393 * t123 + t389 * t130;
t106 = -pkin(5) * t119 - t390 * t107 + t394 * t110;
t105 = pkin(1) * t162 + pkin(5) * t120 + t394 * t107 + t390 * t110;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, 0, 0, 0, 0, 0, 0, -t440, -t312, 0, t220, 0, 0, 0, 0, 0, 0, -t454, t251, 0, t146, 0, 0, 0, 0, 0, 0, -t470, t197, 0, t116, 0, 0, 0, 0, 0, 0, t151, t152, t160, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, 0, 0, 0, 0, 0, 0, t312, -t440, 0, -t459, 0, 0, 0, 0, 0, 0, -t251, -t454, 0, -t475, 0, 0, 0, 0, 0, 0, -t197, -t470, 0, -t486, 0, 0, 0, 0, 0, 0, t149, t150, t159, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t384, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t384, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t384, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t384, 0, 0, 0, 0, 0, 0, t300, t301, 0, -t162; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t431, -t429, -t315, -qJ(1) * t315, 0, 0, t312, 0, -t440, 0, t458, t452, t459, pkin(5) * t430 + qJ(1) * t459 - t385 * t260, 0, 0, -t251, 0, -t454, 0, t482, t481, t475, qJ(1) * t475 - t385 * t143 + t386 * t144, 0, 0, -t197, 0, -t470, 0, t490, t489, t486, qJ(1) * t486 - t385 * t111 + t386 * t112, -t385 * t184 + t386 * t186, -t385 * t173 + t386 * t174, -t385 * t190 + t386 * t192, -t385 * t183 + t386 * t185, -t385 * t189 + t386 * t191, -t385 * t217 + t386 * t218, -qJ(1) * t149 - t385 * t113 + t386 * t117, -qJ(1) * t150 - t385 * t114 + t386 * t118, -qJ(1) * t159 - t385 * t121 + t386 * t122, -qJ(1) * t108 - t385 * t105 + t386 * t106; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t429, -t431, t316, qJ(1) * t316, 0, 0, t440, 0, t312, 0, -t452, t458, t220, pkin(5) * t434 + qJ(1) * t220 + t386 * t260, 0, 0, t454, 0, -t251, 0, -t481, t482, t146, qJ(1) * t146 + t386 * t143 + t385 * t144, 0, 0, t470, 0, -t197, 0, -t489, t490, t116, qJ(1) * t116 + t386 * t111 + t385 * t112, t386 * t184 + t385 * t186, t386 * t173 + t385 * t174, t386 * t190 + t385 * t192, t386 * t183 + t385 * t185, t386 * t189 + t385 * t191, t386 * t217 + t385 * t218, qJ(1) * t151 + t386 * t113 + t385 * t117, qJ(1) * t152 + t386 * t114 + t385 * t118, qJ(1) * t160 + t386 * t121 + t385 * t122, qJ(1) * t109 + t386 * t105 + t385 * t106; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t366, t367, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(1) * t364 + t401, -pkin(1) * t363 + t412, 0, -pkin(1) * t267, 0, 0, 0, 0, 0, t380, -pkin(1) * t294 - pkin(2) * t349 + t397, -pkin(1) * t290 - pkin(2) * t346 - t262, 0, -pkin(1) * t177 - pkin(2) * t215, 0, 0, 0, 0, 0, t375, -pkin(1) * t243 - pkin(2) * t288 - pkin(3) * t342 - t208, -pkin(1) * t238 - pkin(2) * t283 - pkin(3) * t339 - t209, 0, -pkin(1) * t128 - pkin(2) * t141 - pkin(3) * t167, (t330 + t408) * t387, t391 * t329 + t387 * t334, t391 * t358 + t427, (t333 - t409) * t391, t387 * t360 + t418, 0, pkin(1) * t179 + pkin(2) * t221 + pkin(3) * t269 + pkin(4) * t334 + pkin(8) * t303 - t419, pkin(1) * t180 + pkin(2) * t222 + pkin(3) * t270 - pkin(4) * t329 + pkin(8) * t305 + t428, pkin(1) * t199 + pkin(2) * t237 + pkin(3) * t282 + pkin(4) * t343 + pkin(8) * t337 + t163, pkin(1) * t119 + pkin(2) * t135 + pkin(3) * t147 - pkin(4) * t205 + pkin(8) * t163;];
tauB_reg = t1;
