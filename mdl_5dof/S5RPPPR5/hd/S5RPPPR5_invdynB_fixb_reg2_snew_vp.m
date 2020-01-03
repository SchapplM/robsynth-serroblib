% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPPR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:34
% EndTime: 2019-12-31 17:46:39
% DurationCPUTime: 4.26s
% Computational Cost: add. (12630->364), mult. (24696->526), div. (0->0), fcn. (13081->8), ass. (0->237)
t363 = sin(qJ(5));
t361 = cos(pkin(8));
t365 = cos(qJ(5));
t359 = sin(pkin(8));
t405 = t359 * t363;
t320 = (-t361 * t365 + t405) * qJD(1);
t376 = t359 * t365 + t361 * t363;
t321 = t376 * qJD(1);
t408 = t320 * t321;
t426 = qJDD(5) - t408;
t428 = t363 * t426;
t427 = t365 * t426;
t360 = sin(pkin(7));
t388 = t360 * qJDD(1);
t362 = cos(pkin(7));
t413 = qJD(1) ^ 2;
t392 = t362 * t413;
t330 = -t388 + t392;
t387 = t362 * qJDD(1);
t331 = t360 * t413 + t387;
t364 = sin(qJ(1));
t366 = cos(qJ(1));
t285 = t366 * t330 + t364 * t331;
t358 = g(3) + qJDD(3);
t305 = qJ(3) * t331 + t360 * t358;
t379 = qJ(3) * t330 + t362 * t358;
t425 = -pkin(5) * t285 + t364 * t305 + t366 * t379;
t356 = qJDD(1) * qJ(2);
t338 = t366 * g(1) + t364 * g(2);
t375 = 0.2e1 * qJD(2) * qJD(1) - t338;
t372 = t356 + t375;
t412 = pkin(1) + pkin(2);
t303 = -t412 * t413 + t372;
t337 = t364 * g(1) - t366 * g(2);
t374 = -qJDD(2) + t337;
t371 = -t413 * qJ(2) - t374;
t370 = -t412 * qJDD(1) + t371;
t265 = t362 * t303 + t360 * t370;
t262 = -t413 * pkin(3) - qJDD(1) * qJ(4) + t265;
t391 = qJD(1) * qJD(4);
t424 = t262 - 0.2e1 * t391;
t354 = t359 ^ 2;
t355 = t361 ^ 2;
t416 = t413 * (t354 + t355);
t326 = t361 * t416;
t383 = t361 * t388;
t294 = -t362 * t326 + t383;
t380 = -t360 * t326 - t361 * t387;
t254 = t364 * t294 - t366 * t380;
t256 = t366 * t294 + t364 * t380;
t381 = -t364 * t330 + t366 * t331;
t422 = -pkin(5) * t381 + t366 * t305 - t364 * t379;
t264 = t360 * t303 - t362 * t370;
t225 = t362 * t264 - t360 * t265;
t226 = t360 * t264 + t362 * t265;
t174 = t366 * t225 + t364 * t226;
t421 = t364 * t225 - t366 * t226;
t351 = t355 * t413;
t393 = t354 * t413;
t417 = -t393 - t351;
t319 = t376 * qJDD(1);
t316 = t320 ^ 2;
t317 = t321 ^ 2;
t414 = 0.2e1 * t359;
t411 = pkin(4) * t361;
t410 = qJDD(1) * pkin(1);
t409 = qJDD(1) * pkin(3);
t346 = t361 * t358;
t233 = t391 * t414 + t346 + (pkin(6) * qJDD(1) + t413 * t411 - t262) * t359;
t243 = t359 * t358 + t424 * t361;
t390 = qJDD(1) * t361;
t234 = -pkin(4) * t351 - pkin(6) * t390 + t243;
t188 = -t365 * t233 + t363 * t234;
t189 = t363 * t233 + t365 * t234;
t169 = -t365 * t188 + t363 * t189;
t407 = t359 * t169;
t406 = t359 * t361;
t373 = -t413 * qJ(4) + qJDD(4) + t264;
t261 = t373 + t409;
t404 = t360 * t261;
t403 = t361 * t169;
t402 = t362 * t261;
t241 = (pkin(3) + t411) * qJDD(1) + t373 + t417 * pkin(6);
t401 = t363 * t241;
t273 = qJDD(5) + t408;
t400 = t363 * t273;
t397 = t365 * t241;
t396 = t365 * t273;
t263 = qJDD(1) * t405 - t365 * t390;
t395 = qJD(5) * t321;
t394 = t320 * qJD(5);
t389 = t355 * qJDD(1);
t385 = t360 * t408;
t384 = t362 * t408;
t382 = t261 + t409;
t170 = t363 * t188 + t365 * t189;
t311 = -t413 * pkin(1) + t372;
t313 = -t371 + t410;
t268 = t366 * t311 - t364 * t313;
t297 = -t364 * t337 - t366 * t338;
t335 = t364 * qJDD(1) + t366 * t413;
t315 = -pkin(5) * t335 + t366 * g(3);
t336 = t366 * qJDD(1) - t364 * t413;
t314 = pkin(5) * t336 + t364 * g(3);
t242 = t424 * t359 - t346;
t208 = t361 * t242 - t359 * t243;
t209 = t359 * t242 + t361 * t243;
t300 = t359 * t383 - t392 * t406;
t301 = t331 * t406;
t378 = t366 * t300 - t364 * t301;
t377 = t364 * t300 + t366 * t301;
t267 = t364 * t311 + t366 * t313;
t296 = t366 * t337 - t364 * t338;
t367 = qJD(5) ^ 2;
t349 = t354 * qJDD(1);
t334 = t351 - t393;
t329 = -t349 - t389;
t328 = t349 - t389;
t325 = t359 * t416;
t310 = -t317 - t367;
t309 = -t317 + t367;
t308 = t316 - t367;
t292 = t362 * t325 - t359 * t388;
t289 = t360 * t325 + t359 * t387;
t284 = t362 * t329 + t360 * t417;
t283 = t362 * t328 - t360 * t334;
t282 = t360 * t329 - t362 * t417;
t281 = -t360 * t328 - t362 * t334;
t279 = -t317 + t316;
t278 = t394 - t319;
t277 = t319 - 0.2e1 * t394;
t276 = 0.2e1 * t395 + t263;
t275 = t395 + t263;
t271 = -t367 - t316;
t270 = (t320 * t365 - t321 * t363) * qJD(5);
t269 = (t320 * t363 + t321 * t365) * qJD(5);
t266 = -t316 - t317;
t260 = t365 * t278 + t363 * t395;
t259 = t363 * t278 - t365 * t395;
t258 = -t363 * t275 - t365 * t394;
t257 = t365 * t275 - t363 * t394;
t255 = t364 * t289 + t366 * t292;
t253 = -t366 * t289 + t364 * t292;
t251 = -t363 * t310 - t396;
t250 = -t363 * t309 + t427;
t249 = t365 * t308 - t400;
t248 = t365 * t310 - t400;
t247 = t365 * t309 + t428;
t246 = t363 * t308 + t396;
t245 = t364 * t282 + t366 * t284;
t244 = -t366 * t282 + t364 * t284;
t240 = t365 * t276 + t363 * t277;
t239 = t365 * t263 - t363 * t319;
t238 = t363 * t276 - t365 * t277;
t237 = t363 * t263 + t365 * t319;
t236 = t365 * t271 - t428;
t235 = t363 * t271 + t427;
t231 = -t359 * t269 + t361 * t270;
t228 = t360 * qJDD(5) + t362 * t231;
t227 = t362 * qJDD(5) - t360 * t231;
t222 = qJ(2) * t358 + qJ(3) * t225;
t221 = -qJ(3) * t226 + t412 * t358;
t220 = -t359 * t259 + t361 * t260;
t219 = -t359 * t257 + t361 * t258;
t218 = -t359 * t248 + t361 * t251;
t217 = -t359 * t247 + t361 * t250;
t216 = -t359 * t246 + t361 * t249;
t215 = t361 * t248 + t359 * t251;
t214 = -pkin(6) * t248 + t397;
t213 = t362 * t217 - t360 * t319;
t212 = t362 * t216 + t263 * t360;
t211 = -t360 * t217 - t362 * t319;
t210 = -t360 * t216 + t263 * t362;
t206 = -pkin(6) * t235 + t401;
t205 = -t359 * t238 + t361 * t240;
t204 = -t359 * t237 + t361 * t239;
t203 = t361 * t237 + t359 * t239;
t202 = -t359 * t235 + t361 * t236;
t201 = t361 * t235 + t359 * t236;
t200 = t362 * t220 + t385;
t199 = t362 * t219 - t385;
t198 = -t360 * t220 + t384;
t197 = -t360 * t219 - t384;
t196 = -qJ(3) * t289 - t360 * t243 + t361 * t402;
t195 = -qJ(3) * t380 - t360 * t242 + t359 * t402;
t194 = -qJ(3) * t292 - t362 * t243 - t361 * t404;
t193 = -qJ(3) * t294 - t362 * t242 - t359 * t404;
t192 = t362 * t218 - t360 * t277;
t191 = t360 * t218 + t362 * t277;
t190 = pkin(4) * t277 + pkin(6) * t251 + t401;
t186 = -qJ(3) * t282 + t362 * t208;
t185 = -qJ(3) * t284 - t360 * t208;
t184 = pkin(4) * t276 + pkin(6) * t236 - t397;
t183 = t362 * t205 - t360 * t279;
t182 = -t360 * t205 - t362 * t279;
t181 = t362 * t202 - t360 * t276;
t180 = t360 * t202 + t362 * t276;
t179 = t362 * t204 + t360 * t266;
t178 = t360 * t204 - t362 * t266;
t177 = t362 * t209 + t404;
t176 = t360 * t209 - t402;
t173 = -pkin(3) * t203 - pkin(4) * t237;
t172 = t364 * t191 + t366 * t192;
t171 = -t366 * t191 + t364 * t192;
t168 = -pkin(3) * t215 - pkin(4) * t248 + t189;
t167 = t364 * t180 + t366 * t181;
t166 = -t366 * t180 + t364 * t181;
t165 = -pkin(3) * t201 - pkin(4) * t235 + t188;
t164 = t364 * t178 + t366 * t179;
t163 = -t366 * t178 + t364 * t179;
t162 = -pkin(6) * t237 - t169;
t161 = t364 * t176 + t366 * t177;
t160 = -t366 * t176 + t364 * t177;
t159 = -pkin(4) * t241 + pkin(6) * t170;
t158 = -qJ(4) * t215 - t359 * t190 + t361 * t214;
t157 = -pkin(4) * t266 + pkin(6) * t239 + t170;
t156 = -qJ(4) * t201 - t359 * t184 + t361 * t206;
t155 = -qJ(3) * t176 - (pkin(3) * t360 - qJ(4) * t362 + qJ(2)) * t208;
t154 = t361 * t170 - t407;
t153 = t359 * t170 + t403;
t152 = t362 * t154 + t360 * t241;
t151 = t360 * t154 - t362 * t241;
t150 = -qJ(3) * t177 - (pkin(3) * t362 + qJ(4) * t360 + t412) * t208;
t149 = -qJ(4) * t203 - t359 * t157 + t361 * t162;
t148 = -pkin(3) * t153 - pkin(4) * t169;
t147 = qJ(2) * t215 - qJ(3) * t191 + t362 * t158 - t360 * t168;
t146 = -qJ(3) * t192 - t360 * t158 - t362 * t168 + t412 * t215;
t145 = qJ(2) * t201 - qJ(3) * t180 + t362 * t156 - t360 * t165;
t144 = -qJ(3) * t181 - t360 * t156 - t362 * t165 + t412 * t201;
t143 = -pkin(6) * t403 - qJ(4) * t153 - t359 * t159;
t142 = t364 * t151 + t366 * t152;
t141 = -t366 * t151 + t364 * t152;
t140 = qJ(2) * t203 - qJ(3) * t178 + t362 * t149 - t360 * t173;
t139 = -qJ(3) * t179 - t360 * t149 - t362 * t173 + t412 * t203;
t138 = qJ(2) * t153 - qJ(3) * t151 + t362 * t143 - t360 * t148;
t137 = -qJ(3) * t152 - t360 * t143 - t362 * t148 + t412 * t153;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t335, -t336, 0, t297, 0, 0, 0, 0, 0, 0, -t335, 0, t336, t268, 0, 0, 0, 0, 0, 0, -t285, t381, 0, -t421, 0, 0, 0, 0, 0, 0, t256, t255, t245, t161, 0, 0, 0, 0, 0, 0, t167, t172, t164, t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t336, -t335, 0, t296, 0, 0, 0, 0, 0, 0, t336, 0, t335, t267, 0, 0, 0, 0, 0, 0, t381, t285, 0, t174, 0, 0, 0, 0, 0, 0, t254, t253, t244, t160, 0, 0, 0, 0, 0, 0, t166, t171, t163, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t358, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, 0, 0, 0, 0, 0, 0, -t201, -t215, -t203, -t153; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t336, 0, -t335, 0, -t314, -t315, -t296, -pkin(5) * t296, 0, t336, 0, 0, t335, 0, -t314, -t267, t315, -pkin(5) * t267 + (-pkin(1) * t364 + qJ(2) * t366) * g(3), 0, 0, -t381, 0, -t285, 0, t422, t425, t174, -pkin(5) * t174 - t364 * t221 + t366 * t222, -t377, -t364 * t281 + t366 * t283, t255, t377, -t256, 0, -pkin(5) * t254 - t364 * t193 + t366 * t195, -pkin(5) * t253 - t364 * t194 + t366 * t196, -pkin(5) * t244 - t364 * t185 + t366 * t186, -pkin(5) * t160 - t364 * t150 + t366 * t155, -t364 * t198 + t366 * t200, -t364 * t182 + t366 * t183, -t364 * t211 + t366 * t213, -t364 * t197 + t366 * t199, -t364 * t210 + t366 * t212, -t364 * t227 + t366 * t228, -pkin(5) * t166 - t364 * t144 + t366 * t145, -pkin(5) * t171 - t364 * t146 + t366 * t147, -pkin(5) * t163 - t364 * t139 + t366 * t140, -pkin(5) * t141 - t364 * t137 + t366 * t138; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t335, 0, t336, 0, t315, -t314, t297, pkin(5) * t297, 0, t335, 0, 0, -t336, 0, t315, t268, t314, pkin(5) * t268 + (pkin(1) * t366 + qJ(2) * t364) * g(3), 0, 0, -t285, 0, t381, 0, t425, -t422, t421, -pkin(5) * t421 + t366 * t221 + t364 * t222, t378, t366 * t281 + t364 * t283, t253, -t378, -t254, 0, pkin(5) * t256 + t366 * t193 + t364 * t195, pkin(5) * t255 + t366 * t194 + t364 * t196, pkin(5) * t245 + t366 * t185 + t364 * t186, pkin(5) * t161 + t366 * t150 + t364 * t155, t366 * t198 + t364 * t200, t366 * t182 + t364 * t183, t366 * t211 + t364 * t213, t366 * t197 + t364 * t199, t366 * t210 + t364 * t212, t366 * t227 + t364 * t228, pkin(5) * t167 + t366 * t144 + t364 * t145, pkin(5) * t172 + t366 * t146 + t364 * t147, pkin(5) * t164 + t366 * t139 + t364 * t140, pkin(5) * t142 + t366 * t137 + t364 * t138; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t337, t338, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t374 + 0.2e1 * t410, 0, 0.2e1 * t356 + t375, pkin(1) * t313 + qJ(2) * t311, 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t330 + t412 * t331 + t264, qJ(2) * t331 + t412 * t330 + t265, 0, qJ(2) * t226 + t225 * t412, t349, t390 * t414, 0, t389, 0, 0, qJ(2) * t294 + qJ(4) * t326 + t382 * t361 - t380 * t412, qJ(2) * t292 - qJ(4) * t325 - t412 * t289 - t382 * t359, pkin(3) * t417 + qJ(2) * t284 - qJ(4) * t329 - t412 * t282 - t209, pkin(3) * t261 + qJ(2) * t177 - qJ(4) * t209 - t412 * t176, -t361 * t259 - t359 * t260, -t361 * t238 - t359 * t240, -t361 * t247 - t359 * t250, -t361 * t257 - t359 * t258, -t361 * t246 - t359 * t249, -t361 * t269 - t359 * t270, -pkin(3) * t276 + qJ(2) * t181 - qJ(4) * t202 - t412 * t180 - t361 * t184 - t359 * t206, -pkin(3) * t277 + qJ(2) * t192 - qJ(4) * t218 - t361 * t190 - t412 * t191 - t359 * t214, pkin(3) * t266 + qJ(2) * t179 - qJ(4) * t204 - t361 * t157 - t359 * t162 - t412 * t178, pkin(3) * t241 + pkin(6) * t407 + qJ(2) * t152 - qJ(4) * t154 - t412 * t151 - t361 * t159;];
tauB_reg = t1;
