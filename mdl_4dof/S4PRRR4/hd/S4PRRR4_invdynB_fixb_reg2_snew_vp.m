% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRR4
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:41
% EndTime: 2019-12-31 16:32:46
% DurationCPUTime: 3.83s
% Computational Cost: add. (9363->362), mult. (19353->546), div. (0->0), fcn. (13253->8), ass. (0->253)
t354 = sin(pkin(7));
t359 = sin(qJ(2));
t362 = cos(qJ(2));
t355 = cos(pkin(7));
t329 = t355 * g(1) + t354 * g(2);
t373 = t354 * g(1) - t355 * g(2);
t369 = t359 * t329 + t362 * t373;
t386 = t362 * t329 - t359 * t373;
t372 = -t359 * t369 - t362 * t386;
t241 = t359 * t386 - t362 * t369;
t399 = t355 * t241;
t420 = -t354 * t372 + t399;
t404 = t354 * t241;
t195 = t355 * t372 + t404;
t364 = qJD(2) ^ 2;
t380 = t359 * qJDD(2);
t325 = t362 * t364 + t380;
t379 = t362 * qJDD(2);
t326 = -t359 * t364 + t379;
t277 = -t354 * t325 + t355 * t326;
t352 = g(3) - qJDD(1);
t304 = pkin(4) * t325 - t362 * t352;
t366 = -pkin(4) * t326 - t359 * t352;
t419 = -qJ(1) * t277 + t354 * t304 + t355 * t366;
t357 = sin(qJ(4));
t360 = cos(qJ(4));
t361 = cos(qJ(3));
t358 = sin(qJ(3));
t384 = qJD(2) * t358;
t307 = -t360 * t361 * qJD(2) + t357 * t384;
t309 = (t357 * t361 + t358 * t360) * qJD(2);
t273 = t309 * t307;
t378 = qJDD(3) + qJDD(4);
t411 = -t273 + t378;
t418 = t357 * t411;
t417 = t360 * t411;
t409 = t355 * t325 + t354 * t326;
t415 = qJ(1) * t409 + t355 * t304 - t354 * t366;
t349 = qJD(3) + qJD(4);
t300 = t349 * t307;
t382 = qJD(2) * qJD(3);
t374 = t361 * t382;
t381 = t358 * qJDD(2);
t321 = t374 + t381;
t343 = t361 * qJDD(2);
t375 = t358 * t382;
t322 = t343 - t375;
t365 = t307 * qJD(4) - t360 * t321 - t357 * t322;
t410 = -t300 - t365;
t371 = t357 * t321 - t360 * t322;
t223 = (qJD(4) - t349) * t309 + t371;
t283 = -t355 * t329 - t354 * t373;
t282 = -t354 * t329 + t355 * t373;
t305 = t307 ^ 2;
t306 = t309 ^ 2;
t347 = t349 ^ 2;
t407 = t349 * t357;
t406 = t349 * t360;
t350 = t358 ^ 2;
t405 = t350 * t364;
t351 = t361 ^ 2;
t345 = t351 * t364;
t400 = t354 * t352;
t398 = t355 * t352;
t270 = -qJDD(2) * pkin(2) - t364 * pkin(5) - t369;
t332 = qJD(3) * pkin(3) - pkin(6) * t384;
t230 = -t322 * pkin(3) - pkin(6) * t345 + t332 * t384 + t270;
t397 = t357 * t230;
t268 = t273 + t378;
t396 = t357 * t268;
t271 = -t364 * pkin(2) + qJDD(2) * pkin(5) - t386;
t251 = t358 * t271 + t361 * t352;
t337 = t358 * t364 * t361;
t330 = qJDD(3) + t337;
t213 = (-t321 + t374) * pkin(6) + t330 * pkin(3) - t251;
t254 = t361 * t271 - t358 * t352;
t218 = -pkin(3) * t345 + t322 * pkin(6) - qJD(3) * t332 + t254;
t178 = -t360 * t213 + t357 * t218;
t179 = t357 * t213 + t360 * t218;
t150 = -t360 * t178 + t357 * t179;
t395 = t358 * t150;
t394 = t358 * t270;
t393 = t358 * t330;
t331 = qJDD(3) - t337;
t392 = t358 * t331;
t391 = t360 * t230;
t390 = t360 * t268;
t389 = t361 * t150;
t388 = t361 * t270;
t387 = t361 * t331;
t385 = t350 + t351;
t377 = t359 * t273;
t376 = t362 * t273;
t151 = t357 * t178 + t360 * t179;
t207 = t358 * t251 + t361 * t254;
t368 = t359 * t337;
t367 = t362 * t337;
t206 = t361 * t251 - t358 * t254;
t363 = qJD(3) ^ 2;
t336 = -t345 - t363;
t335 = t345 - t363;
t334 = -t363 - t405;
t333 = t363 - t405;
t328 = t345 - t405;
t327 = t345 + t405;
t324 = t385 * qJDD(2);
t323 = t343 - 0.2e1 * t375;
t320 = 0.2e1 * t374 + t381;
t319 = t361 * t330;
t318 = t385 * t382;
t298 = -t306 + t347;
t297 = t305 - t347;
t296 = t359 * qJDD(3) + t362 * t318;
t295 = t361 * t321 - t350 * t382;
t294 = -t362 * qJDD(3) + t359 * t318;
t293 = -t358 * t322 - t351 * t382;
t292 = -t306 - t347;
t291 = -t358 * t334 - t387;
t290 = -t358 * t333 + t319;
t289 = t361 * t336 - t393;
t288 = t361 * t335 - t392;
t287 = t361 * t334 - t392;
t286 = t358 * t336 + t319;
t281 = t362 * t324 - t359 * t327;
t280 = t359 * t324 + t362 * t327;
t274 = -t358 * t320 + t361 * t323;
t272 = -t306 + t305;
t266 = t362 * t295 - t368;
t265 = t362 * t293 + t368;
t264 = t359 * t295 + t367;
t263 = t359 * t293 - t367;
t262 = t362 * t290 + t358 * t380;
t261 = t362 * t288 + t359 * t343;
t260 = t359 * t290 - t358 * t379;
t259 = t359 * t288 - t361 * t379;
t257 = -t347 - t305;
t256 = t362 * t291 + t359 * t320;
t255 = t362 * t289 - t359 * t323;
t253 = t359 * t291 - t362 * t320;
t252 = t359 * t289 + t362 * t323;
t250 = (-t307 * t360 + t309 * t357) * t349;
t249 = (-t307 * t357 - t309 * t360) * t349;
t247 = -t305 - t306;
t246 = t362 * t274 - t359 * t328;
t245 = t359 * t274 + t362 * t328;
t243 = -t309 * qJD(4) - t371;
t238 = -pkin(5) * t287 + t388;
t237 = -pkin(5) * t286 + t394;
t236 = t360 * t297 - t396;
t235 = -t357 * t298 + t417;
t234 = t357 * t297 + t390;
t233 = t360 * t298 + t418;
t232 = -t354 * t280 + t355 * t281;
t231 = t355 * t280 + t354 * t281;
t229 = -t357 * t292 - t390;
t228 = t360 * t292 - t396;
t227 = -t300 + t365;
t222 = (qJD(4) + t349) * t309 + t371;
t221 = pkin(1) * t352 + pkin(4) * t372;
t220 = -pkin(2) * t287 + t254;
t219 = -pkin(2) * t286 + t251;
t217 = -t309 * t407 - t360 * t365;
t216 = t309 * t406 - t357 * t365;
t215 = -t357 * t243 + t307 * t406;
t214 = t360 * t243 + t307 * t407;
t209 = t360 * t257 - t418;
t208 = t357 * t257 + t417;
t204 = -t354 * t253 + t355 * t256;
t203 = -t354 * t252 + t355 * t255;
t202 = t355 * t253 + t354 * t256;
t201 = t355 * t252 + t354 * t255;
t200 = -t358 * t249 + t361 * t250;
t199 = t362 * t200 + t359 * t378;
t198 = t359 * t200 - t362 * t378;
t197 = -pkin(4) * t280 + t362 * t206;
t196 = pkin(4) * t281 + t359 * t206;
t193 = -t358 * t234 + t361 * t236;
t192 = -t358 * t233 + t361 * t235;
t191 = t362 * t207 + t359 * t270;
t190 = t359 * t207 - t362 * t270;
t189 = -pkin(6) * t228 + t391;
t188 = -t358 * t228 + t361 * t229;
t187 = t361 * t228 + t358 * t229;
t186 = -t223 * t360 - t357 * t227;
t185 = -t360 * t222 - t357 * t410;
t184 = -t223 * t357 + t360 * t227;
t183 = -t357 * t222 + t360 * t410;
t182 = -pkin(6) * t208 + t397;
t181 = -t358 * t216 + t361 * t217;
t180 = -t358 * t214 + t361 * t215;
t176 = -t358 * t208 + t361 * t209;
t175 = t361 * t208 + t358 * t209;
t174 = -pkin(4) * t253 - t359 * t220 + t362 * t238;
t173 = -pkin(4) * t252 - t359 * t219 + t362 * t237;
t172 = t362 * t181 + t377;
t171 = t362 * t180 - t377;
t170 = t359 * t181 - t376;
t169 = t359 * t180 + t376;
t168 = -pkin(3) * t410 + pkin(6) * t229 + t397;
t167 = t362 * t193 - t359 * t223;
t166 = t362 * t192 - t359 * t227;
t165 = t359 * t193 + t362 * t223;
t164 = t359 * t192 + t362 * t227;
t163 = -pkin(1) * t287 + pkin(4) * t256 + t362 * t220 + t359 * t238;
t162 = -pkin(1) * t286 + pkin(4) * t255 + t362 * t219 + t359 * t237;
t161 = t362 * t188 + t359 * t410;
t160 = t359 * t188 - t362 * t410;
t159 = -pkin(3) * t222 + pkin(6) * t209 - t391;
t158 = t362 * t176 + t359 * t222;
t157 = t359 * t176 - t362 * t222;
t156 = -t354 * t190 + t355 * t191;
t155 = t355 * t190 + t354 * t191;
t154 = -t358 * t184 + t361 * t186;
t153 = -t358 * t183 + t361 * t185;
t152 = t361 * t184 + t358 * t186;
t149 = t362 * t153 - t359 * t272;
t148 = t359 * t153 + t362 * t272;
t147 = -pkin(4) * t190 - (pkin(2) * t359 - pkin(5) * t362) * t206;
t146 = t362 * t154 + t359 * t247;
t145 = t359 * t154 - t362 * t247;
t144 = -pkin(2) * t187 - pkin(3) * t228 + t179;
t143 = -pkin(3) * t230 + pkin(6) * t151;
t142 = -pkin(2) * t175 - pkin(3) * t208 + t178;
t141 = pkin(4) * t191 - (-pkin(2) * t362 - pkin(5) * t359 - pkin(1)) * t206;
t140 = -pkin(2) * t152 - pkin(3) * t184;
t139 = -t354 * t160 + t355 * t161;
t138 = t355 * t160 + t354 * t161;
t137 = -pkin(6) * t184 - t150;
t136 = -t354 * t157 + t355 * t158;
t135 = t355 * t157 + t354 * t158;
t134 = -pkin(5) * t187 - t358 * t168 + t361 * t189;
t133 = -pkin(3) * t247 + pkin(6) * t186 + t151;
t132 = -pkin(5) * t175 - t358 * t159 + t361 * t182;
t131 = t361 * t151 - t395;
t130 = t358 * t151 + t389;
t129 = -t354 * t145 + t355 * t146;
t128 = t355 * t145 + t354 * t146;
t127 = t362 * t131 + t359 * t230;
t126 = t359 * t131 - t362 * t230;
t125 = -pkin(2) * t130 - pkin(3) * t150;
t124 = -pkin(4) * t160 + t362 * t134 - t359 * t144;
t123 = -pkin(4) * t157 + t362 * t132 - t359 * t142;
t122 = -pkin(1) * t187 + pkin(4) * t161 + t359 * t134 + t362 * t144;
t121 = -pkin(5) * t152 - t358 * t133 + t361 * t137;
t120 = -pkin(1) * t175 + pkin(4) * t158 + t359 * t132 + t362 * t142;
t119 = -pkin(5) * t130 - pkin(6) * t389 - t358 * t143;
t118 = -t354 * t126 + t355 * t127;
t117 = t355 * t126 + t354 * t127;
t116 = -pkin(4) * t145 + t362 * t121 - t359 * t140;
t115 = -pkin(1) * t152 + pkin(4) * t146 + t359 * t121 + t362 * t140;
t114 = -pkin(4) * t126 + t362 * t119 - t359 * t125;
t113 = -pkin(1) * t130 + pkin(4) * t127 + t359 * t119 + t362 * t125;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, 0, 0, 0, 0, 0, 0, -t409, -t277, 0, t195, 0, 0, 0, 0, 0, 0, t203, t204, t232, t156, 0, 0, 0, 0, 0, 0, t136, t139, t129, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, 0, 0, 0, 0, 0, 0, t277, -t409, 0, -t420, 0, 0, 0, 0, 0, 0, t201, t202, t231, t155, 0, 0, 0, 0, 0, 0, t135, t138, t128, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t352, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t352, 0, 0, 0, 0, 0, 0, t286, t287, 0, -t206, 0, 0, 0, 0, 0, 0, t175, t187, t152, t130; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t400, -t398, -t282, -qJ(1) * t282, 0, 0, t277, 0, -t409, 0, t419, t415, t420, pkin(4) * t399 + qJ(1) * t420 - t354 * t221, -t354 * t264 + t355 * t266, -t354 * t245 + t355 * t246, -t354 * t260 + t355 * t262, -t354 * t263 + t355 * t265, -t354 * t259 + t355 * t261, -t354 * t294 + t355 * t296, -qJ(1) * t201 - t354 * t162 + t355 * t173, -qJ(1) * t202 - t354 * t163 + t355 * t174, -qJ(1) * t231 - t354 * t196 + t355 * t197, -qJ(1) * t155 - t354 * t141 + t355 * t147, -t354 * t170 + t355 * t172, -t354 * t148 + t355 * t149, -t354 * t164 + t355 * t166, -t354 * t169 + t355 * t171, -t354 * t165 + t355 * t167, -t354 * t198 + t355 * t199, -qJ(1) * t135 - t354 * t120 + t355 * t123, -qJ(1) * t138 - t354 * t122 + t355 * t124, -qJ(1) * t128 - t354 * t115 + t355 * t116, -qJ(1) * t117 - t354 * t113 + t355 * t114; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t398, -t400, t283, qJ(1) * t283, 0, 0, t409, 0, t277, 0, -t415, t419, t195, pkin(4) * t404 + qJ(1) * t195 + t355 * t221, t355 * t264 + t354 * t266, t355 * t245 + t354 * t246, t355 * t260 + t354 * t262, t355 * t263 + t354 * t265, t355 * t259 + t354 * t261, t355 * t294 + t354 * t296, qJ(1) * t203 + t355 * t162 + t354 * t173, qJ(1) * t204 + t355 * t163 + t354 * t174, qJ(1) * t232 + t355 * t196 + t354 * t197, qJ(1) * t156 + t355 * t141 + t354 * t147, t355 * t170 + t354 * t172, t355 * t148 + t354 * t149, t355 * t164 + t354 * t166, t355 * t169 + t354 * t171, t355 * t165 + t354 * t167, t355 * t198 + t354 * t199, qJ(1) * t136 + t355 * t120 + t354 * t123, qJ(1) * t139 + t355 * t122 + t354 * t124, qJ(1) * t129 + t355 * t115 + t354 * t116, qJ(1) * t118 + t355 * t113 + t354 * t114; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t373, t329, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(1) * t326 + t369, -pkin(1) * t325 + t386, 0, -pkin(1) * t241, (t321 + t374) * t358, t361 * t320 + t358 * t323, t361 * t333 + t393, (t322 - t375) * t361, t358 * t335 + t387, 0, pkin(1) * t252 + pkin(2) * t323 + pkin(5) * t289 - t388, pkin(1) * t253 - pkin(2) * t320 + pkin(5) * t291 + t394, pkin(1) * t280 + pkin(2) * t327 + pkin(5) * t324 + t207, pkin(1) * t190 - pkin(2) * t270 + pkin(5) * t207, t361 * t216 + t358 * t217, t361 * t183 + t358 * t185, t361 * t233 + t358 * t235, t361 * t214 + t358 * t215, t361 * t234 + t358 * t236, t361 * t249 + t358 * t250, pkin(1) * t157 - pkin(2) * t222 + pkin(5) * t176 + t361 * t159 + t358 * t182, pkin(1) * t160 - pkin(2) * t410 + pkin(5) * t188 + t361 * t168 + t358 * t189, pkin(1) * t145 - pkin(2) * t247 + pkin(5) * t154 + t361 * t133 + t358 * t137, pkin(1) * t126 - pkin(2) * t230 + pkin(5) * t131 - pkin(6) * t395 + t361 * t143;];
tauB_reg = t1;
