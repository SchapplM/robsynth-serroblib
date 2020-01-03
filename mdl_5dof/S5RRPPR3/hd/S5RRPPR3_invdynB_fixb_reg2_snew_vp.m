% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPPR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:37
% EndTime: 2019-12-31 19:26:44
% DurationCPUTime: 4.49s
% Computational Cost: add. (12640->322), mult. (17357->449), div. (0->0), fcn. (9986->8), ass. (0->210)
t346 = (qJD(1) + qJD(2));
t344 = t346 ^ 2;
t345 = qJDD(1) + qJDD(2);
t350 = sin(pkin(8));
t351 = cos(pkin(8));
t312 = t350 * t344 - t351 * t345;
t349 = g(3) - qJDD(3);
t289 = qJ(3) * t312 - t350 * t349;
t353 = sin(qJ(2));
t356 = cos(qJ(2));
t309 = t351 * t344 + t350 * t345;
t366 = t353 * t309 + t356 * t312;
t371 = qJ(3) * t309 - t351 * t349;
t211 = pkin(6) * t366 + t356 * t289 + t353 * t371;
t257 = t356 * t309 - t353 * t312;
t354 = sin(qJ(1));
t357 = cos(qJ(1));
t222 = t357 * t257 - t354 * t366;
t421 = pkin(6) * t257 - t353 * t289 + t356 * t371;
t437 = pkin(5) * t222 - t354 * t211 + t357 * t421;
t417 = t354 * t257 + t357 * t366;
t425 = pkin(5) * t417 + t357 * t211 + t354 * t421;
t332 = t354 * g(1) - t357 * g(2);
t322 = qJDD(1) * pkin(1) + t332;
t333 = t357 * g(1) + t354 * g(2);
t359 = qJD(1) ^ 2;
t323 = -t359 * pkin(1) - t333;
t274 = -t356 * t322 + t353 * t323;
t267 = t345 * pkin(2) - t274;
t275 = t353 * t322 + t356 * t323;
t268 = -t344 * pkin(2) + t275;
t232 = -t351 * t267 + t350 * t268;
t233 = t350 * t267 + t351 * t268;
t369 = t350 * t232 + t351 * t233;
t184 = t351 * t232 - t350 * t233;
t381 = t356 * t184;
t159 = -t353 * t369 + t381;
t389 = t353 * t184;
t429 = t356 * t369 + t389;
t142 = t354 * t159 + t357 * t429;
t443 = t357 * t159 - t354 * t429;
t440 = -pkin(1) * t257 - pkin(2) * t309;
t439 = -pkin(1) * t366 - pkin(2) * t312;
t360 = (2 * qJD(4) * t346) + t233;
t400 = t345 * qJ(4);
t203 = -t344 * pkin(3) + t360 + t400;
t338 = t345 * pkin(3);
t361 = qJDD(4) + t232 - t338;
t216 = -t344 * qJ(4) + t361;
t175 = t350 * t203 - t351 * t216;
t370 = t351 * t203 + t350 * t216;
t154 = t356 * t175 + t353 * t370;
t419 = -t353 * t175 + t356 * t370;
t138 = -t354 * t154 + t357 * t419;
t137 = t357 * t154 + t354 * t419;
t315 = t356 * t344 + t353 * t345;
t318 = t353 * t344 - t356 * t345;
t271 = t354 * t315 + t357 * t318;
t297 = pkin(6) * t315 - t356 * g(3);
t420 = pkin(6) * t318 - t353 * g(3);
t436 = pkin(5) * t271 + t354 * t297 + t357 * t420;
t362 = t357 * t315 - t354 * t318;
t435 = pkin(5) * t362 + t357 * t297 - t354 * t420;
t368 = t353 * t274 + t356 * t275;
t237 = t356 * t274 - t353 * t275;
t378 = t357 * t237;
t430 = -t354 * t368 + t378;
t386 = t354 * t237;
t190 = t357 * t368 + t386;
t404 = pkin(3) + pkin(7);
t403 = pkin(1) * t349;
t352 = sin(qJ(5));
t347 = t352 ^ 2;
t399 = t347 * t344;
t355 = cos(qJ(5));
t348 = t355 ^ 2;
t398 = t348 * t344;
t377 = t347 + t348;
t314 = t377 * t345;
t397 = t350 * t314;
t395 = t351 * t314;
t199 = -t344 * pkin(7) + t203;
t393 = t352 * t199;
t375 = t352 * t344 * t355;
t324 = qJDD(5) + t375;
t392 = t352 * t324;
t325 = qJDD(5) - t375;
t391 = t352 * t325;
t390 = t352 * t345;
t385 = t355 * t199;
t384 = t355 * t324;
t383 = t355 * t325;
t382 = t355 * t345;
t376 = qJD(5) * t346;
t374 = t352 * t376;
t373 = t355 * t376;
t202 = -t345 * pkin(7) + t216;
t197 = -t355 * t202 - t352 * t349;
t287 = -t354 * t332 - t357 * t333;
t365 = t350 * t375;
t364 = t351 * t375;
t327 = t357 * qJDD(1) - t354 * t359;
t363 = -pkin(5) * t327 - t354 * g(3);
t198 = t352 * t202 - t355 * t349;
t170 = -t355 * t197 + t352 * t198;
t171 = t352 * t197 + t355 * t198;
t286 = t357 * t332 - t354 * t333;
t358 = qJD(5) ^ 2;
t331 = -t358 - t398;
t330 = t358 - t398;
t329 = -t358 - t399;
t328 = -t358 + t399;
t326 = t354 * qJDD(1) + t357 * t359;
t320 = (-t347 + t348) * t344;
t319 = t377 * t344;
t306 = -0.2e1 * t374 + t382;
t305 = -t374 + t382;
t304 = -t373 - t390;
t303 = 0.2e1 * t373 + t390;
t302 = -pkin(5) * t326 + t357 * g(3);
t301 = t377 * t376;
t285 = t351 * qJDD(5) - t350 * t301;
t284 = t350 * qJDD(5) + t351 * t301;
t283 = -t352 * t305 - t348 * t376;
t282 = -t355 * t304 - t347 * t376;
t281 = -t352 * t331 - t384;
t280 = t355 * t329 - t391;
t279 = t355 * t331 - t392;
t278 = -t355 * t330 - t391;
t277 = t352 * t329 + t383;
t276 = -t352 * t328 - t384;
t266 = -t351 * t319 - t397;
t265 = -t350 * t319 + t395;
t256 = t352 * t303 - t355 * t306;
t254 = -t350 * t278 + t351 * t382;
t253 = -t350 * t276 - t351 * t390;
t252 = t351 * t278 + t350 * t382;
t251 = t351 * t276 - t350 * t390;
t250 = -t350 * t282 - t364;
t249 = -t350 * t283 + t364;
t248 = t351 * t282 - t365;
t247 = t351 * t283 + t365;
t246 = t350 * t279 + t351 * t306;
t245 = t350 * t277 + t351 * t303;
t244 = -t351 * t279 + t350 * t306;
t243 = -t351 * t277 + t350 * t303;
t242 = -t353 * t284 + t356 * t285;
t241 = t356 * t284 + t353 * t285;
t240 = -t350 * t256 + t351 * t320;
t239 = t351 * t256 + t350 * t320;
t234 = pkin(1) * g(3) + pkin(6) * t368;
t231 = -t353 * t265 + t356 * t266;
t230 = t356 * t265 + t353 * t266;
t220 = -t353 * t252 + t356 * t254;
t219 = -t353 * t251 + t356 * t253;
t218 = t356 * t252 + t353 * t254;
t217 = t356 * t251 + t353 * t253;
t207 = -t353 * t248 + t356 * t250;
t206 = -t353 * t247 + t356 * t249;
t205 = t356 * t248 + t353 * t250;
t204 = t356 * t247 + t353 * t249;
t196 = -t353 * t244 + t356 * t246;
t195 = -t353 * t243 + t356 * t245;
t194 = t356 * t244 + t353 * t246;
t193 = t356 * t243 + t353 * t245;
t192 = -t353 * t239 + t356 * t240;
t191 = t356 * t239 + t353 * t240;
t188 = pkin(4) * t279 - qJ(4) * t281 - t198;
t187 = pkin(4) * t277 - qJ(4) * t280 - t197;
t186 = -t354 * t230 + t357 * t231;
t183 = t357 * t230 + t354 * t231;
t180 = pkin(2) * t349 + qJ(3) * t369;
t179 = pkin(4) * t303 - t404 * t280 + t385;
t178 = pkin(4) * t306 - t404 * t281 - t393;
t173 = -qJ(3) * t175 + (-pkin(3) * t350 + qJ(4) * t351) * t349;
t172 = qJ(3) * t370 + (pkin(3) * t351 + qJ(4) * t350 + pkin(2)) * t349;
t169 = -t354 * t194 + t357 * t196;
t168 = -t354 * t193 + t357 * t195;
t167 = t357 * t194 + t354 * t196;
t166 = t357 * t193 + t354 * t195;
t165 = -pkin(4) * t319 - t171;
t164 = -pkin(4) * t395 - qJ(3) * t265 - t350 * t165;
t163 = -pkin(4) * t397 + qJ(3) * t266 + t351 * t165;
t162 = t350 * t170 + t351 * t199;
t161 = -t351 * t170 + t350 * t199;
t152 = -qJ(3) * t244 - t350 * t178 + t351 * t188;
t151 = -qJ(3) * t243 - t350 * t179 + t351 * t187;
t150 = -pkin(2) * t281 + qJ(3) * t246 + t351 * t178 + t350 * t188;
t149 = -pkin(2) * t280 + qJ(3) * t245 + t351 * t179 + t350 * t187;
t148 = pkin(4) * t170 - qJ(4) * t171;
t147 = pkin(4) * t199 - t404 * t171;
t146 = -t353 * t161 + t356 * t162;
t145 = t356 * t161 + t353 * t162;
t144 = -pkin(6) * t230 - t353 * t163 + t356 * t164;
t143 = pkin(6) * t231 + t356 * t163 + t353 * t164;
t140 = pkin(6) * t159 + qJ(3) * t381 - t353 * t180;
t139 = pkin(6) * t429 + qJ(3) * t389 + t356 * t180 + t403;
t136 = -pkin(6) * t154 - t353 * t172 + t356 * t173;
t135 = pkin(6) * t419 + t356 * t172 + t353 * t173 + t403;
t134 = -pkin(6) * t194 - t353 * t150 + t356 * t152;
t133 = -pkin(6) * t193 - t353 * t149 + t356 * t151;
t132 = -pkin(1) * t281 + pkin(6) * t196 + t356 * t150 + t353 * t152;
t131 = -pkin(1) * t280 + pkin(6) * t195 + t356 * t149 + t353 * t151;
t130 = -qJ(3) * t161 - t350 * t147 + t351 * t148;
t129 = -t354 * t145 + t357 * t146;
t128 = t357 * t145 + t354 * t146;
t127 = -pkin(2) * t171 + qJ(3) * t162 + t351 * t147 + t350 * t148;
t126 = -pkin(6) * t145 - t353 * t127 + t356 * t130;
t125 = -pkin(1) * t171 + pkin(6) * t146 + t356 * t127 + t353 * t130;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t326, -t327, 0, t287, 0, 0, 0, 0, 0, 0, -t362, t271, 0, t190, 0, 0, 0, 0, 0, 0, -t222, t417, 0, t142, 0, 0, 0, 0, 0, 0, 0, t222, -t417, t138, 0, 0, 0, 0, 0, 0, t168, t169, t186, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t327, -t326, 0, t286, 0, 0, 0, 0, 0, 0, -t271, -t362, 0, -t430, 0, 0, 0, 0, 0, 0, -t417, -t222, 0, -t443, 0, 0, 0, 0, 0, 0, 0, t417, t222, t137, 0, 0, 0, 0, 0, 0, t166, t167, t183, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t349, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t349, 0, 0, 0, 0, 0, 0, t280, t281, 0, t171; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t327, 0, -t326, 0, t363, -t302, -t286, -pkin(5) * t286, 0, 0, -t271, 0, -t362, 0, t436, t435, t430, pkin(5) * t430 + pkin(6) * t378 - t354 * t234, 0, 0, -t417, 0, -t222, 0, t425, t437, t443, pkin(5) * t443 - t354 * t139 + t357 * t140, 0, t417, t222, 0, 0, 0, -t137, -t425, -t437, -pkin(5) * t137 - t354 * t135 + t357 * t136, -t354 * t204 + t357 * t206, -t354 * t191 + t357 * t192, -t354 * t218 + t357 * t220, -t354 * t205 + t357 * t207, -t354 * t217 + t357 * t219, -t354 * t241 + t357 * t242, -pkin(5) * t166 - t354 * t131 + t357 * t133, -pkin(5) * t167 - t354 * t132 + t357 * t134, -pkin(5) * t183 - t354 * t143 + t357 * t144, -pkin(5) * t128 - t354 * t125 + t357 * t126; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t326, 0, t327, 0, t302, t363, t287, pkin(5) * t287, 0, 0, t362, 0, -t271, 0, -t435, t436, t190, pkin(5) * t190 + pkin(6) * t386 + t357 * t234, 0, 0, t222, 0, -t417, 0, -t437, t425, t142, pkin(5) * t142 + t357 * t139 + t354 * t140, 0, -t222, t417, 0, 0, 0, t138, t437, -t425, pkin(5) * t138 + t357 * t135 + t354 * t136, t357 * t204 + t354 * t206, t357 * t191 + t354 * t192, t357 * t218 + t354 * t220, t357 * t205 + t354 * t207, t357 * t217 + t354 * t219, t357 * t241 + t354 * t242, pkin(5) * t168 + t357 * t131 + t354 * t133, pkin(5) * t169 + t357 * t132 + t354 * t134, pkin(5) * t186 + t357 * t143 + t354 * t144, pkin(5) * t129 + t357 * t125 + t354 * t126; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t332, t333, 0, 0, 0, 0, 0, 0, 0, t345, -pkin(1) * t318 - t274, -pkin(1) * t315 - t275, 0, -pkin(1) * t237, 0, 0, 0, 0, 0, t345, -t232 + t439, -t233 + t440, 0, -pkin(1) * t159 - pkin(2) * t184, t345, 0, 0, 0, 0, 0, 0, -t338 + t361 - t439, t360 + 0.2e1 * t400 - t440, pkin(1) * t154 + pkin(2) * t175 - pkin(3) * t216 + qJ(4) * t203, (t305 - t374) * t355, -t355 * t303 - t352 * t306, -t352 * t330 + t383, (-t304 + t373) * t352, t355 * t328 - t392, 0, pkin(1) * t193 + pkin(2) * t243 + qJ(4) * t303 - t404 * t277 + t393, pkin(1) * t194 + pkin(2) * t244 + qJ(4) * t306 - t404 * t279 + t385, pkin(1) * t230 + pkin(2) * t265 - qJ(4) * t319 + t404 * t314 - t170, pkin(1) * t145 + pkin(2) * t161 + qJ(4) * t199 - t404 * t170;];
tauB_reg = t1;
