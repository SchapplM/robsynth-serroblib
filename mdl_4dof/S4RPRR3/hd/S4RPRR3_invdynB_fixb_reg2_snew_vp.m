% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRR3
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:21
% EndTime: 2019-12-31 16:49:25
% DurationCPUTime: 4.01s
% Computational Cost: add. (10332->370), mult. (20830->554), div. (0->0), fcn. (13277->8), ass. (0->255)
t377 = sin(qJ(1));
t380 = cos(qJ(1));
t349 = t380 * g(1) + t377 * g(2);
t382 = qJD(1) ^ 2;
t333 = -t382 * pkin(1) - t349;
t372 = sin(pkin(7));
t373 = cos(pkin(7));
t348 = t377 * g(1) - t380 * g(2);
t384 = qJDD(1) * pkin(1) + t348;
t286 = t372 * t333 - t373 * t384;
t287 = t373 * t333 + t372 * t384;
t391 = t372 * t286 + t373 * t287;
t236 = t373 * t286 - t372 * t287;
t402 = t380 * t236;
t432 = -t377 * t391 + t402;
t410 = t377 * t236;
t196 = t380 * t391 + t410;
t338 = t372 * qJDD(1) + t373 * t382;
t339 = t373 * qJDD(1) - t372 * t382;
t290 = -t377 * t338 + t380 * t339;
t370 = g(3) - qJDD(2);
t314 = qJ(2) * t338 - t373 * t370;
t385 = -qJ(2) * t339 - t372 * t370;
t431 = -pkin(4) * t290 + t377 * t314 + t380 * t385;
t375 = sin(qJ(4));
t378 = cos(qJ(4));
t379 = cos(qJ(3));
t376 = sin(qJ(3));
t400 = qJD(1) * t376;
t323 = -t378 * t379 * qJD(1) + t375 * t400;
t325 = (t375 * t379 + t376 * t378) * qJD(1);
t288 = t325 * t323;
t396 = qJDD(3) + qJDD(4);
t423 = -t288 + t396;
t430 = t375 * t423;
t429 = t378 * t423;
t421 = t380 * t338 + t377 * t339;
t427 = pkin(4) * t421 + t380 * t314 - t377 * t385;
t367 = qJD(3) + qJD(4);
t316 = t367 * t323;
t398 = qJD(1) * qJD(3);
t392 = t379 * t398;
t397 = t376 * qJDD(1);
t335 = t392 + t397;
t361 = t379 * qJDD(1);
t393 = t376 * t398;
t336 = t361 - t393;
t383 = t323 * qJD(4) - t378 * t335 - t375 * t336;
t422 = -t316 - t383;
t390 = t375 * t335 - t378 * t336;
t242 = (qJD(4) - t367) * t325 + t390;
t318 = t323 ^ 2;
t319 = t325 ^ 2;
t365 = t367 ^ 2;
t419 = t367 * t375;
t418 = t367 * t378;
t368 = t376 ^ 2;
t417 = t368 * t382;
t369 = t379 ^ 2;
t363 = t369 * t382;
t270 = -qJDD(1) * pkin(2) - t382 * pkin(5) + t286;
t347 = qJD(3) * pkin(3) - pkin(6) * t400;
t238 = -t336 * pkin(3) - pkin(6) * t363 + t347 * t400 + t270;
t416 = t375 * t238;
t283 = t288 + t396;
t415 = t375 * t283;
t271 = -t382 * pkin(2) + qJDD(1) * pkin(5) + t287;
t257 = t376 * t271 + t379 * t370;
t355 = t379 * t382 * t376;
t345 = qJDD(3) + t355;
t223 = (-t335 + t392) * pkin(6) + t345 * pkin(3) - t257;
t259 = t379 * t271 - t376 * t370;
t225 = -pkin(3) * t363 + t336 * pkin(6) - qJD(3) * t347 + t259;
t188 = -t378 * t223 + t375 * t225;
t189 = t375 * t223 + t378 * t225;
t160 = -t378 * t188 + t375 * t189;
t414 = t376 * t160;
t413 = t376 * t270;
t412 = t376 * t345;
t346 = qJDD(3) - t355;
t411 = t376 * t346;
t407 = t378 * t238;
t406 = t378 * t283;
t405 = t379 * t160;
t404 = t379 * t270;
t403 = t379 * t346;
t401 = t368 + t369;
t395 = t372 * t288;
t394 = t373 * t288;
t161 = t375 * t188 + t378 * t189;
t214 = t376 * t257 + t379 * t259;
t297 = -t377 * t348 - t380 * t349;
t388 = t372 * t355;
t387 = t373 * t355;
t342 = t380 * qJDD(1) - t377 * t382;
t386 = -pkin(4) * t342 - t377 * g(3);
t213 = t379 * t257 - t376 * t259;
t296 = t380 * t348 - t377 * t349;
t381 = qJD(3) ^ 2;
t353 = -t363 - t381;
t352 = t363 - t381;
t351 = -t381 - t417;
t350 = t381 - t417;
t344 = t363 - t417;
t343 = t363 + t417;
t341 = t377 * qJDD(1) + t380 * t382;
t340 = t401 * qJDD(1);
t337 = t361 - 0.2e1 * t393;
t334 = 0.2e1 * t392 + t397;
t331 = t379 * t345;
t330 = t401 * t398;
t317 = -pkin(4) * t341 + t380 * g(3);
t310 = -t319 + t365;
t309 = t318 - t365;
t308 = t379 * t335 - t368 * t398;
t307 = -t376 * t336 - t369 * t398;
t306 = t372 * qJDD(3) + t373 * t330;
t305 = -t373 * qJDD(3) + t372 * t330;
t304 = -t319 - t365;
t303 = -t376 * t351 - t403;
t302 = -t376 * t350 + t331;
t301 = t379 * t353 - t412;
t300 = t379 * t352 - t411;
t299 = t379 * t351 - t411;
t298 = t376 * t353 + t331;
t295 = t373 * t340 - t372 * t343;
t294 = t372 * t340 + t373 * t343;
t289 = -t376 * t334 + t379 * t337;
t285 = -t319 + t318;
t280 = -t365 - t318;
t279 = t373 * t308 - t388;
t278 = t373 * t307 + t388;
t277 = t372 * t308 + t387;
t276 = t372 * t307 - t387;
t275 = t373 * t302 + t372 * t397;
t274 = t373 * t300 + t372 * t361;
t273 = t372 * t302 - t373 * t397;
t272 = t372 * t300 - t373 * t361;
t268 = t373 * t303 + t372 * t334;
t267 = t373 * t301 - t372 * t337;
t266 = t372 * t303 - t373 * t334;
t265 = t372 * t301 + t373 * t337;
t264 = (-t323 * t378 + t325 * t375) * t367;
t263 = (-t323 * t375 - t325 * t378) * t367;
t262 = -t318 - t319;
t260 = -t325 * qJD(4) - t390;
t258 = t373 * t289 - t372 * t344;
t256 = t372 * t289 + t373 * t344;
t254 = t378 * t309 - t415;
t253 = -t375 * t310 + t429;
t252 = t375 * t309 + t406;
t251 = t378 * t310 + t430;
t250 = -t377 * t294 + t380 * t295;
t249 = t380 * t294 + t377 * t295;
t248 = -t375 * t304 - t406;
t247 = t378 * t304 - t415;
t246 = -t316 + t383;
t241 = (qJD(4) + t367) * t325 + t390;
t240 = -pkin(5) * t299 + t404;
t239 = -pkin(5) * t298 + t413;
t233 = -t325 * t419 - t378 * t383;
t232 = t325 * t418 - t375 * t383;
t231 = -t375 * t260 + t323 * t418;
t230 = t378 * t260 + t323 * t419;
t229 = -pkin(2) * t299 + t259;
t228 = -pkin(2) * t298 + t257;
t227 = t378 * t280 - t430;
t226 = t375 * t280 + t429;
t224 = pkin(1) * t370 + qJ(2) * t391;
t219 = -t377 * t266 + t380 * t268;
t218 = -t377 * t265 + t380 * t267;
t217 = t380 * t266 + t377 * t268;
t216 = t380 * t265 + t377 * t267;
t215 = -t376 * t263 + t379 * t264;
t211 = t373 * t215 + t372 * t396;
t210 = t372 * t215 - t373 * t396;
t209 = -t376 * t252 + t379 * t254;
t208 = -t376 * t251 + t379 * t253;
t207 = -qJ(2) * t294 + t373 * t213;
t206 = qJ(2) * t295 + t372 * t213;
t205 = -t376 * t247 + t379 * t248;
t204 = t379 * t247 + t376 * t248;
t203 = -t242 * t378 - t375 * t246;
t202 = -t378 * t241 - t375 * t422;
t201 = -t242 * t375 + t378 * t246;
t200 = -t375 * t241 + t378 * t422;
t199 = -pkin(6) * t247 + t407;
t198 = t373 * t214 + t372 * t270;
t197 = t372 * t214 - t373 * t270;
t194 = -pkin(6) * t226 + t416;
t193 = -t376 * t232 + t379 * t233;
t192 = -t376 * t230 + t379 * t231;
t191 = -t376 * t226 + t379 * t227;
t190 = t379 * t226 + t376 * t227;
t186 = t373 * t193 + t395;
t185 = t373 * t192 - t395;
t184 = t372 * t193 - t394;
t183 = t372 * t192 + t394;
t182 = -qJ(2) * t266 - t372 * t229 + t373 * t240;
t181 = -qJ(2) * t265 - t372 * t228 + t373 * t239;
t180 = t373 * t209 - t372 * t242;
t179 = t373 * t208 - t372 * t246;
t178 = t372 * t209 + t373 * t242;
t177 = t372 * t208 + t373 * t246;
t176 = -pkin(3) * t422 + pkin(6) * t248 + t416;
t175 = t373 * t205 + t372 * t422;
t174 = t372 * t205 - t373 * t422;
t173 = -pkin(1) * t299 + qJ(2) * t268 + t373 * t229 + t372 * t240;
t172 = -pkin(1) * t298 + qJ(2) * t267 + t373 * t228 + t372 * t239;
t171 = -pkin(3) * t241 + pkin(6) * t227 - t407;
t170 = t373 * t191 + t372 * t241;
t169 = t372 * t191 - t373 * t241;
t168 = -t376 * t201 + t379 * t203;
t167 = -t376 * t200 + t379 * t202;
t166 = t379 * t201 + t376 * t203;
t165 = -t377 * t197 + t380 * t198;
t164 = t380 * t197 + t377 * t198;
t163 = t373 * t167 - t372 * t285;
t162 = t372 * t167 + t373 * t285;
t159 = t373 * t168 + t372 * t262;
t158 = t372 * t168 - t373 * t262;
t157 = -qJ(2) * t197 - (pkin(2) * t372 - pkin(5) * t373) * t213;
t156 = -pkin(2) * t204 - pkin(3) * t247 + t189;
t155 = -pkin(2) * t190 - pkin(3) * t226 + t188;
t154 = -pkin(3) * t238 + pkin(6) * t161;
t153 = -pkin(2) * t166 - pkin(3) * t201;
t152 = -t377 * t174 + t380 * t175;
t151 = t380 * t174 + t377 * t175;
t150 = qJ(2) * t198 - (-pkin(2) * t373 - pkin(5) * t372 - pkin(1)) * t213;
t149 = -pkin(6) * t201 - t160;
t148 = -t377 * t169 + t380 * t170;
t147 = t380 * t169 + t377 * t170;
t146 = -pkin(5) * t204 - t376 * t176 + t379 * t199;
t145 = -pkin(3) * t262 + pkin(6) * t203 + t161;
t144 = -pkin(5) * t190 - t376 * t171 + t379 * t194;
t143 = t379 * t161 - t414;
t142 = t376 * t161 + t405;
t141 = -t377 * t158 + t380 * t159;
t140 = t380 * t158 + t377 * t159;
t139 = t373 * t143 + t372 * t238;
t138 = t372 * t143 - t373 * t238;
t137 = -pkin(2) * t142 - pkin(3) * t160;
t136 = -qJ(2) * t174 + t373 * t146 - t372 * t156;
t135 = -qJ(2) * t169 + t373 * t144 - t372 * t155;
t134 = -pkin(1) * t204 + qJ(2) * t175 + t372 * t146 + t373 * t156;
t133 = -pkin(5) * t166 - t376 * t145 + t379 * t149;
t132 = -pkin(1) * t190 + qJ(2) * t170 + t372 * t144 + t373 * t155;
t131 = -pkin(5) * t142 - pkin(6) * t405 - t376 * t154;
t130 = -t377 * t138 + t380 * t139;
t129 = t380 * t138 + t377 * t139;
t128 = -qJ(2) * t158 + t373 * t133 - t372 * t153;
t127 = -pkin(1) * t166 + qJ(2) * t159 + t372 * t133 + t373 * t153;
t126 = -qJ(2) * t138 + t373 * t131 - t372 * t137;
t125 = -pkin(1) * t142 + qJ(2) * t139 + t372 * t131 + t373 * t137;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t341, -t342, 0, t297, 0, 0, 0, 0, 0, 0, -t421, -t290, 0, t196, 0, 0, 0, 0, 0, 0, t218, t219, t250, t165, 0, 0, 0, 0, 0, 0, t148, t152, t141, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t342, -t341, 0, t296, 0, 0, 0, 0, 0, 0, t290, -t421, 0, -t432, 0, 0, 0, 0, 0, 0, t216, t217, t249, t164, 0, 0, 0, 0, 0, 0, t147, t151, t140, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t370, 0, 0, 0, 0, 0, 0, t298, t299, 0, -t213, 0, 0, 0, 0, 0, 0, t190, t204, t166, t142; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t342, 0, -t341, 0, t386, -t317, -t296, -pkin(4) * t296, 0, 0, t290, 0, -t421, 0, t431, t427, t432, pkin(4) * t432 + qJ(2) * t402 - t377 * t224, -t377 * t277 + t380 * t279, -t377 * t256 + t380 * t258, -t377 * t273 + t380 * t275, -t377 * t276 + t380 * t278, -t377 * t272 + t380 * t274, -t377 * t305 + t380 * t306, -pkin(4) * t216 - t377 * t172 + t380 * t181, -pkin(4) * t217 - t377 * t173 + t380 * t182, -pkin(4) * t249 - t377 * t206 + t380 * t207, -pkin(4) * t164 - t377 * t150 + t380 * t157, -t377 * t184 + t380 * t186, -t377 * t162 + t380 * t163, -t377 * t177 + t380 * t179, -t377 * t183 + t380 * t185, -t377 * t178 + t380 * t180, -t377 * t210 + t380 * t211, -pkin(4) * t147 - t377 * t132 + t380 * t135, -pkin(4) * t151 - t377 * t134 + t380 * t136, -pkin(4) * t140 - t377 * t127 + t380 * t128, -pkin(4) * t129 - t377 * t125 + t380 * t126; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t341, 0, t342, 0, t317, t386, t297, pkin(4) * t297, 0, 0, t421, 0, t290, 0, -t427, t431, t196, pkin(4) * t196 + qJ(2) * t410 + t380 * t224, t380 * t277 + t377 * t279, t380 * t256 + t377 * t258, t380 * t273 + t377 * t275, t380 * t276 + t377 * t278, t380 * t272 + t377 * t274, t380 * t305 + t377 * t306, pkin(4) * t218 + t380 * t172 + t377 * t181, pkin(4) * t219 + t380 * t173 + t377 * t182, pkin(4) * t250 + t380 * t206 + t377 * t207, pkin(4) * t165 + t380 * t150 + t377 * t157, t380 * t184 + t377 * t186, t380 * t162 + t377 * t163, t380 * t177 + t377 * t179, t380 * t183 + t377 * t185, t380 * t178 + t377 * t180, t380 * t210 + t377 * t211, pkin(4) * t148 + t380 * t132 + t377 * t135, pkin(4) * t152 + t380 * t134 + t377 * t136, pkin(4) * t141 + t380 * t127 + t377 * t128, pkin(4) * t130 + t380 * t125 + t377 * t126; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t348, t349, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t339 - t286, -pkin(1) * t338 - t287, 0, -pkin(1) * t236, (t335 + t392) * t376, t379 * t334 + t376 * t337, t379 * t350 + t412, (t336 - t393) * t379, t376 * t352 + t403, 0, pkin(1) * t265 + pkin(2) * t337 + pkin(5) * t301 - t404, pkin(1) * t266 - pkin(2) * t334 + pkin(5) * t303 + t413, pkin(1) * t294 + pkin(2) * t343 + pkin(5) * t340 + t214, pkin(1) * t197 - pkin(2) * t270 + pkin(5) * t214, t379 * t232 + t376 * t233, t379 * t200 + t376 * t202, t379 * t251 + t376 * t253, t379 * t230 + t376 * t231, t379 * t252 + t376 * t254, t379 * t263 + t376 * t264, pkin(1) * t169 - pkin(2) * t241 + pkin(5) * t191 + t379 * t171 + t376 * t194, pkin(1) * t174 - pkin(2) * t422 + pkin(5) * t205 + t379 * t176 + t376 * t199, pkin(1) * t158 - pkin(2) * t262 + pkin(5) * t168 + t379 * t145 + t376 * t149, pkin(1) * t138 - pkin(2) * t238 + pkin(5) * t143 - pkin(6) * t414 + t379 * t154;];
tauB_reg = t1;
