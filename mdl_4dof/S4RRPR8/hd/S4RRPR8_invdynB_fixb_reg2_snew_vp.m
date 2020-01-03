% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPR8
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:26
% EndTime: 2019-12-31 17:08:31
% DurationCPUTime: 3.44s
% Computational Cost: add. (5857->338), mult. (12735->454), div. (0->0), fcn. (7251->6), ass. (0->246)
t377 = qJD(2) ^ 2;
t372 = sin(qJ(2));
t369 = t372 ^ 2;
t378 = qJD(1) ^ 2;
t416 = t369 * t378;
t349 = t377 + t416;
t375 = cos(qJ(2));
t354 = t375 * t378 * t372;
t345 = qJDD(2) - t354;
t414 = t375 * t345;
t298 = -t349 * t372 + t414;
t407 = qJD(1) * qJD(2);
t396 = t375 * t407;
t405 = qJDD(1) * t372;
t334 = 0.2e1 * t396 + t405;
t373 = sin(qJ(1));
t376 = cos(qJ(1));
t263 = t298 * t373 + t334 * t376;
t459 = pkin(4) * t263;
t266 = t298 * t376 - t334 * t373;
t458 = pkin(4) * t266;
t457 = pkin(5) * t298;
t445 = 2 * qJD(3);
t421 = t345 * t372;
t292 = t349 * t375 + t421;
t456 = pkin(1) * t292;
t455 = pkin(5) * t292;
t371 = sin(qJ(4));
t374 = cos(qJ(4));
t318 = (-t371 * t372 - t374 * t375) * qJD(1);
t335 = t396 + t405;
t397 = t372 * t407;
t403 = qJDD(1) * t375;
t336 = -t397 + t403;
t255 = qJD(4) * t318 + t335 * t374 - t336 * t371;
t366 = qJD(2) - qJD(4);
t432 = t318 * t366;
t454 = t255 - t432;
t337 = -0.2e1 * t397 + t403;
t425 = t337 * t375;
t429 = t334 * t372;
t280 = -t425 + t429;
t370 = t375 ^ 2;
t343 = (t369 - t370) * t378;
t453 = t280 * t373 + t343 * t376;
t452 = t280 * t376 - t343 * t373;
t415 = t370 * t378;
t351 = -t377 + t415;
t296 = -t351 * t375 + t421;
t402 = qJDD(1) * t376;
t451 = t296 * t373 + t375 * t402;
t450 = t296 * t376 - t373 * t403;
t365 = -qJDD(2) + qJDD(4);
t408 = qJD(1) * t375;
t409 = qJD(1) * t372;
t320 = -t371 * t408 + t374 * t409;
t433 = t318 * t320;
t386 = t365 + t433;
t449 = t371 * t386;
t448 = t374 * t386;
t442 = pkin(2) * t375;
t390 = -qJ(3) * t372 - t442;
t332 = t390 * qJD(1);
t348 = g(1) * t376 + g(2) * t373;
t322 = -pkin(1) * t378 + qJDD(1) * pkin(5) - t348;
t395 = g(3) * t372 - t375 * t322;
t384 = qJDD(2) * qJ(3) + qJD(2) * t445 + t332 * t408 - t395;
t447 = t351 * t372 + t414;
t446 = t332 * t409 + qJDD(3);
t316 = t318 ^ 2;
t317 = t320 ^ 2;
t364 = t366 ^ 2;
t444 = pkin(2) + pkin(3);
t443 = pkin(2) * t336;
t352 = -t377 - t415;
t344 = qJDD(2) + t354;
t422 = t344 * t372;
t295 = t352 * t375 - t422;
t262 = t295 * t373 + t337 * t376;
t441 = pkin(4) * t262;
t410 = t369 + t370;
t339 = t410 * qJDD(1);
t342 = t410 * t378;
t282 = t339 * t373 + t342 * t376;
t440 = pkin(4) * t282;
t327 = t375 * t344;
t290 = t352 * t372 + t327;
t439 = pkin(5) * t290;
t438 = qJ(3) * t375;
t346 = -qJD(2) * pkin(3) - pkin(6) * t409;
t347 = g(1) * t373 - t376 * g(2);
t321 = qJDD(1) * pkin(1) + pkin(5) * t378 + t347;
t381 = -pkin(2) * t397 + t321;
t217 = t443 + qJ(3) * t335 + t336 * pkin(3) - pkin(6) * t415 + (qJD(2) * t438 + (t445 + t346) * t372) * qJD(1) + t381;
t437 = t217 * t371;
t436 = t217 * t374;
t275 = -t365 + t433;
t435 = t275 * t371;
t434 = t275 * t374;
t431 = t321 * t372;
t430 = t321 * t375;
t426 = t337 * t372;
t418 = t366 * t371;
t417 = t366 * t374;
t413 = pkin(1) * t337 + pkin(5) * t295;
t301 = t375 * g(3) + t372 * t322;
t412 = pkin(1) * t342 + pkin(5) * t339;
t411 = t342 - t377;
t404 = qJDD(1) * t373;
t401 = t373 * t433;
t400 = t376 * t433;
t250 = -pkin(2) * t377 + t384;
t220 = -pkin(3) * t415 - pkin(6) * t336 + qJD(2) * t346 + t250;
t383 = -qJDD(2) * pkin(2) + t301 + t446;
t253 = qJ(3) * t377 - t383;
t221 = (-t335 + t396) * pkin(6) - t344 * pkin(3) - t253;
t196 = t220 * t371 - t374 * t221;
t249 = t301 * t372 - t375 * t395;
t394 = t335 * t371 + t374 * t336;
t287 = -t347 * t373 - t376 * t348;
t393 = t373 * t354;
t392 = t376 * t354;
t251 = -pkin(1) * t290 + t301;
t341 = -t373 * t378 + t402;
t391 = -pkin(4) * t341 - g(3) * t373;
t389 = pkin(2) * t372 - t438;
t388 = t335 + t396;
t197 = t220 * t374 + t221 * t371;
t179 = -t196 * t374 + t197 * t371;
t180 = t196 * t371 + t197 * t374;
t248 = t301 * t375 + t372 * t395;
t387 = t334 * t375 + t426;
t286 = t347 * t376 - t348 * t373;
t382 = (-qJD(4) - t366) * t320 - t394;
t380 = t409 * t445 + t381;
t379 = t388 * qJ(3) + t380;
t350 = t377 - t416;
t340 = t376 * t378 + t404;
t330 = t389 * qJDD(1);
t326 = t410 * t407;
t315 = -pkin(4) * t340 + g(3) * t376;
t308 = -t317 + t364;
t307 = t316 - t364;
t306 = qJDD(2) * t373 + t326 * t376;
t305 = t335 * t375 - t369 * t407;
t304 = -qJDD(2) * t376 + t326 * t373;
t303 = -t336 * t372 - t370 * t407;
t300 = -t317 - t364;
t297 = -t350 * t372 + t327;
t291 = t350 * t375 + t422;
t289 = t388 * t372;
t288 = (t336 - t397) * t375;
t283 = t339 * t376 - t342 * t373;
t281 = pkin(4) * t283;
t278 = t317 - t316;
t274 = t305 * t376 - t393;
t273 = t303 * t376 + t393;
t272 = t305 * t373 + t392;
t271 = t303 * t373 - t392;
t270 = t297 * t376 + t372 * t404;
t269 = t297 * t373 - t372 * t402;
t268 = -t364 - t316;
t265 = t295 * t376 - t337 * t373;
t261 = pkin(4) * t265;
t260 = (-t318 * t374 - t320 * t371) * t366;
t259 = (t318 * t371 - t320 * t374) * t366;
t258 = -t430 + t455;
t257 = -t431 - t439;
t256 = -t316 - t317;
t254 = -qJD(4) * t320 - t394;
t252 = -t395 + t456;
t246 = t411 * qJ(3) + t383;
t245 = t411 * pkin(2) + t384;
t244 = t379 + t443;
t243 = t307 * t374 + t435;
t242 = -t308 * t371 + t448;
t241 = -t307 * t371 + t434;
t240 = -t308 * t374 - t449;
t239 = -t300 * t371 + t434;
t238 = t300 * t374 + t435;
t237 = t255 + t432;
t232 = (qJD(4) - t366) * t320 + t394;
t231 = t255 * t374 + t320 * t418;
t230 = -t255 * t371 + t320 * t417;
t229 = -t254 * t371 + t318 * t417;
t228 = -t254 * t374 - t318 * t418;
t227 = (t336 + t337) * pkin(2) + t379;
t226 = t443 + (t334 + t388) * qJ(3) + t380;
t225 = t249 * t376 - t321 * t373;
t224 = t249 * t373 + t321 * t376;
t223 = t268 * t374 - t449;
t222 = t268 * t371 + t448;
t218 = (-t352 - t377) * qJ(3) + (-qJDD(2) - t344) * pkin(2) + t251 + t446;
t216 = -t456 - qJ(3) * t345 + (-t349 + t377) * pkin(2) - t384;
t215 = -t259 * t372 + t260 * t375;
t214 = t250 * t375 - t253 * t372;
t213 = t250 * t372 + t253 * t375;
t212 = -pkin(2) * t429 + t226 * t375 - t455;
t211 = qJ(3) * t425 - t227 * t372 - t439;
t210 = -t245 * t372 + t246 * t375;
t209 = -t241 * t372 + t243 * t375;
t208 = -t240 * t372 + t242 * t375;
t207 = t238 * t372 + t239 * t375;
t206 = -t238 * t375 + t239 * t372;
t205 = t237 * t371 + t374 * t382;
t204 = -t232 * t374 - t371 * t454;
t203 = -t237 * t374 + t371 * t382;
t202 = t232 * t371 - t374 * t454;
t201 = -t230 * t372 + t231 * t375;
t200 = -t228 * t372 + t229 * t375;
t199 = t222 * t372 + t223 * t375;
t198 = -t222 * t375 + t223 * t372;
t195 = t214 * t376 - t244 * t373;
t194 = t214 * t373 + t244 * t376;
t193 = -pkin(1) * t213 - pkin(2) * t253 - qJ(3) * t250;
t192 = t207 * t376 - t373 * t454;
t191 = t207 * t373 + t376 * t454;
t190 = -pkin(6) * t238 + qJ(3) * t454 + t436;
t189 = t199 * t376 - t232 * t373;
t188 = t199 * t373 + t232 * t376;
t187 = -pkin(5) * t213 - t389 * t244;
t186 = -pkin(6) * t222 + qJ(3) * t232 + t437;
t185 = t203 * t372 + t205 * t375;
t184 = -t202 * t372 + t204 * t375;
t183 = -t203 * t375 + t205 * t372;
t182 = -pkin(6) * t239 + t444 * t454 - t437;
t181 = -pkin(6) * t223 + t444 * t232 + t436;
t178 = t185 * t376 - t256 * t373;
t177 = t185 * t373 + t256 * t376;
t176 = -pkin(6) * t179 + qJ(3) * t217;
t175 = -pkin(1) * t206 - qJ(3) * t239 + t444 * t238 - t197;
t174 = -pkin(6) * t203 + qJ(3) * t256 - t179;
t173 = -pkin(6) * t180 + t444 * t217;
t172 = -pkin(1) * t198 - qJ(3) * t223 + t444 * t222 - t196;
t171 = -pkin(6) * t205 + t444 * t256 - t180;
t170 = t179 * t372 + t180 * t375;
t169 = -t179 * t375 + t180 * t372;
t168 = -pkin(5) * t206 - t182 * t372 + t190 * t375;
t167 = -pkin(1) * t183 - qJ(3) * t205 + t444 * t203;
t166 = -pkin(5) * t198 - t181 * t372 + t186 * t375;
t165 = t170 * t376 - t217 * t373;
t164 = t170 * t373 + t217 * t376;
t163 = -pkin(5) * t183 - t171 * t372 + t174 * t375;
t162 = -pkin(5) * t169 - t173 * t372 + t176 * t375;
t161 = -pkin(1) * t169 - qJ(3) * t180 + t444 * t179;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t340, -t341, 0, t287, 0, 0, 0, 0, 0, 0, t265, -t266, t283, t225, 0, 0, 0, 0, 0, 0, t265, t283, t266, t195, 0, 0, 0, 0, 0, 0, t189, t192, t178, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t341, -t340, 0, t286, 0, 0, 0, 0, 0, 0, t262, -t263, t282, t224, 0, 0, 0, 0, 0, 0, t262, t282, t263, t194, 0, 0, 0, 0, 0, 0, t188, t191, t177, t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t290, -t292, 0, -t248, 0, 0, 0, 0, 0, 0, t290, 0, t292, t213, 0, 0, 0, 0, 0, 0, t198, t206, t183, t169; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t341, 0, -t340, 0, t391, -t315, -t286, -pkin(4) * t286, t274, -t452, t270, t273, -t450, t306, -t251 * t373 + t257 * t376 - t441, -t252 * t373 + t258 * t376 + t459, t248 * t376 - t440, -pkin(4) * t224 - (pkin(1) * t373 - pkin(5) * t376) * t248, t274, t270, t452, t306, t450, t273, t211 * t376 - t218 * t373 - t441, t210 * t376 - t330 * t373 - t440, t212 * t376 - t216 * t373 - t459, -pkin(4) * t194 + t187 * t376 - t193 * t373, t201 * t376 + t401, t184 * t376 - t278 * t373, t208 * t376 - t237 * t373, t200 * t376 - t401, t209 * t376 - t373 * t382, t215 * t376 - t365 * t373, -pkin(4) * t188 + t166 * t376 - t172 * t373, -pkin(4) * t191 + t168 * t376 - t175 * t373, -pkin(4) * t177 + t163 * t376 - t167 * t373, -pkin(4) * t164 - t161 * t373 + t162 * t376; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t340, 0, t341, 0, t315, t391, t287, pkin(4) * t287, t272, -t453, t269, t271, -t451, t304, t251 * t376 + t257 * t373 + t261, t252 * t376 + t258 * t373 - t458, t248 * t373 + t281, pkin(4) * t225 - (-pkin(1) * t376 - pkin(5) * t373) * t248, t272, t269, t453, t304, t451, t271, t211 * t373 + t218 * t376 + t261, t210 * t373 + t330 * t376 + t281, t212 * t373 + t216 * t376 + t458, pkin(4) * t195 + t187 * t373 + t193 * t376, t201 * t373 - t400, t184 * t373 + t278 * t376, t208 * t373 + t237 * t376, t200 * t373 + t400, t209 * t373 + t376 * t382, t215 * t373 + t365 * t376, pkin(4) * t189 + t166 * t373 + t172 * t376, pkin(4) * t192 + t168 * t373 + t175 * t376, pkin(4) * t178 + t163 * t373 + t167 * t376, pkin(4) * t165 + t161 * t376 + t162 * t373; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t347, t348, 0, 0, t289, t387, t291, t288, t447, 0, t413 + t430, -pkin(1) * t334 - t431 - t457, t249 + t412, pkin(1) * t321 + pkin(5) * t249, t289, t291, -t387, 0, -t447, t288, qJ(3) * t426 + t227 * t375 + t413, t245 * t375 + t246 * t372 + t412, t457 + t226 * t372 + (pkin(1) + t442) * t334, pkin(5) * t214 + (pkin(1) - t390) * t244, t230 * t375 + t231 * t372, t202 * t375 + t204 * t372, t240 * t375 + t242 * t372, t228 * t375 + t229 * t372, t241 * t375 + t243 * t372, t259 * t375 + t260 * t372, pkin(1) * t232 + pkin(5) * t199 + t181 * t375 + t186 * t372, pkin(1) * t454 + pkin(5) * t207 + t182 * t375 + t190 * t372, pkin(1) * t256 + pkin(5) * t185 + t171 * t375 + t174 * t372, pkin(1) * t217 + pkin(5) * t170 + t173 * t375 + t176 * t372;];
tauB_reg = t1;
