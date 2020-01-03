% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRR6
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRR6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:08
% EndTime: 2019-12-31 16:35:13
% DurationCPUTime: 3.22s
% Computational Cost: add. (8616->358), mult. (17642->548), div. (0->0), fcn. (12008->8), ass. (0->251)
t353 = sin(qJ(4));
t356 = cos(qJ(4));
t357 = cos(qJ(3));
t354 = sin(qJ(3));
t378 = qJD(2) * t354;
t303 = -t356 * t357 * qJD(2) + t353 * t378;
t305 = (t353 * t357 + t354 * t356) * qJD(2);
t268 = t305 * t303;
t371 = qJDD(3) + qJDD(4);
t401 = -t268 + t371;
t405 = t353 * t401;
t404 = t356 * t401;
t350 = sin(pkin(7));
t380 = g(3) - qJDD(1);
t403 = t350 * t380;
t351 = cos(pkin(7));
t402 = t351 * t380;
t346 = qJD(3) + qJD(4);
t300 = t346 * t303;
t375 = qJD(2) * qJD(3);
t367 = t357 * t375;
t374 = t354 * qJDD(2);
t317 = t367 + t374;
t340 = t357 * qJDD(2);
t368 = t354 * t375;
t318 = t340 - t368;
t360 = t303 * qJD(4) - t356 * t317 - t353 * t318;
t400 = -t300 - t360;
t325 = t350 * g(1) - t351 * g(2);
t309 = t351 * t325;
t326 = t351 * g(1) + t350 * g(2);
t274 = -t350 * t326 + t309;
t363 = t353 * t317 - t356 * t318;
t226 = (qJD(4) - t346) * t305 + t363;
t301 = t303 ^ 2;
t302 = t305 ^ 2;
t344 = t346 ^ 2;
t399 = qJD(2) ^ 2;
t398 = t346 * t353;
t397 = t346 * t356;
t355 = sin(qJ(2));
t358 = cos(qJ(2));
t372 = t358 * qJDD(2);
t322 = -t355 * t399 + t372;
t396 = t350 * t322;
t395 = t350 * t325;
t393 = t351 * t322;
t297 = -t355 * t326 + t358 * t380;
t276 = -qJDD(2) * pkin(2) - t399 * pkin(5) + t297;
t329 = qJD(3) * pkin(3) - pkin(6) * t378;
t348 = t357 ^ 2;
t342 = t348 * t399;
t239 = -t318 * pkin(3) - pkin(6) * t342 + t329 * t378 + t276;
t392 = t353 * t239;
t263 = t268 + t371;
t391 = t353 * t263;
t298 = -t358 * t326 - t355 * t380;
t277 = -t399 * pkin(2) + qJDD(2) * pkin(5) + t298;
t251 = t354 * t277 + t357 * t325;
t334 = t354 * t399 * t357;
t327 = qJDD(3) + t334;
t217 = (-t317 + t367) * pkin(6) + t327 * pkin(3) - t251;
t254 = t357 * t277 - t354 * t325;
t222 = -pkin(3) * t342 + t318 * pkin(6) - qJD(3) * t329 + t254;
t178 = -t356 * t217 + t353 * t222;
t181 = t353 * t217 + t356 * t222;
t154 = -t356 * t178 + t353 * t181;
t390 = t354 * t154;
t389 = t354 * t276;
t388 = t354 * t327;
t328 = qJDD(3) - t334;
t387 = t354 * t328;
t386 = t356 * t239;
t385 = t356 * t263;
t384 = t357 * t154;
t383 = t357 * t276;
t382 = t357 * t327;
t381 = t357 * t328;
t347 = t354 ^ 2;
t379 = t347 + t348;
t377 = t347 * t399;
t373 = t355 * qJDD(2);
t370 = t355 * t268;
t369 = t358 * t268;
t207 = t354 * t251 + t357 * t254;
t320 = t379 * qJDD(2);
t323 = t342 + t377;
t272 = t355 * t320 + t358 * t323;
t273 = t358 * t320 - t355 * t323;
t366 = -pkin(1) * t272 - pkin(2) * t323 - pkin(5) * t320 + qJ(1) * t273 - t207;
t321 = t358 * t399 + t373;
t365 = pkin(1) * t322 + qJ(1) * t321 - t297;
t364 = -pkin(1) * t321 + qJ(1) * t322 - t298;
t155 = t353 * t178 + t356 * t181;
t248 = t355 * t297 + t358 * t298;
t275 = -t351 * t326 - t395;
t362 = t355 * t334;
t361 = t358 * t334;
t279 = pkin(4) * t321 - t358 * t325;
t278 = -pkin(4) * t322 - t355 * t325;
t206 = t357 * t251 - t354 * t254;
t247 = t358 * t297 - t355 * t298;
t359 = qJD(3) ^ 2;
t333 = -t342 - t359;
t332 = t342 - t359;
t331 = -t359 - t377;
t330 = t359 - t377;
t324 = t342 - t377;
t319 = t340 - 0.2e1 * t368;
t316 = 0.2e1 * t367 + t374;
t315 = t379 * t375;
t308 = t351 * t321;
t307 = t350 * t321;
t296 = -t302 + t344;
t295 = t301 - t344;
t294 = t355 * qJDD(3) + t358 * t315;
t293 = t357 * t317 - t347 * t375;
t292 = -t354 * t318 - t348 * t375;
t290 = -t302 - t344;
t289 = -t354 * t331 - t381;
t288 = -t354 * t330 + t382;
t287 = t357 * t333 - t388;
t286 = t357 * t332 - t387;
t285 = t357 * t331 - t387;
t284 = -t357 * t330 - t388;
t283 = t354 * t333 + t382;
t282 = -t354 * t332 - t381;
t281 = (-t317 - t367) * t354;
t280 = (-t318 + t368) * t357;
t270 = -t354 * t316 + t357 * t319;
t269 = -t357 * t316 - t354 * t319;
t265 = -t302 + t301;
t261 = t358 * t293 - t362;
t260 = t358 * t292 + t362;
t259 = t358 * t288 + t354 * t373;
t258 = t358 * t286 + t355 * t340;
t257 = -t344 - t301;
t256 = t358 * t289 + t355 * t316;
t255 = t358 * t287 - t355 * t319;
t253 = t355 * t289 - t358 * t316;
t252 = t355 * t287 + t358 * t319;
t250 = (-t303 * t356 + t305 * t353) * t346;
t249 = (-t303 * t353 - t305 * t356) * t346;
t245 = -t301 - t302;
t244 = t358 * t270 - t355 * t324;
t242 = -t305 * qJD(4) - t363;
t241 = -pkin(5) * t285 + t383;
t240 = -pkin(5) * t283 + t389;
t238 = t356 * t295 - t391;
t237 = -t353 * t296 + t404;
t236 = t353 * t295 + t385;
t235 = t356 * t296 + t405;
t234 = t351 * t248 - t395;
t233 = t350 * t248 + t309;
t232 = -t353 * t290 - t385;
t231 = t356 * t290 - t391;
t230 = -t300 + t360;
t225 = (qJD(4) + t346) * t305 + t363;
t224 = -pkin(2) * t285 + t254;
t223 = -pkin(2) * t283 + t251;
t221 = -t305 * t398 - t356 * t360;
t220 = t305 * t397 - t353 * t360;
t219 = -t353 * t242 + t303 * t397;
t218 = t356 * t242 + t303 * t398;
t216 = t351 * t256 + t350 * t285;
t215 = t351 * t255 + t350 * t283;
t214 = t350 * t256 - t351 * t285;
t213 = t350 * t255 - t351 * t283;
t209 = t356 * t257 - t405;
t208 = t353 * t257 + t404;
t204 = -t354 * t249 + t357 * t250;
t203 = -t357 * t249 - t354 * t250;
t202 = t358 * t204 + t355 * t371;
t201 = -pkin(1) * t252 - pkin(2) * t319 - pkin(5) * t287 + t383;
t200 = -pkin(1) * t253 + pkin(2) * t316 - pkin(5) * t289 - t389;
t199 = -pkin(4) * t272 + t358 * t206;
t198 = t358 * t207 + t355 * t276;
t197 = t355 * t207 - t358 * t276;
t196 = -pkin(6) * t231 + t386;
t195 = -t354 * t236 + t357 * t238;
t194 = -t354 * t235 + t357 * t237;
t193 = -t357 * t236 - t354 * t238;
t192 = -t357 * t235 - t354 * t237;
t191 = -t354 * t231 + t357 * t232;
t190 = t357 * t231 + t354 * t232;
t189 = -t226 * t356 - t353 * t230;
t188 = -t356 * t225 - t353 * t400;
t187 = -t226 * t353 + t356 * t230;
t186 = -t353 * t225 + t356 * t400;
t185 = -pkin(6) * t208 + t392;
t183 = -t354 * t220 + t357 * t221;
t182 = -t354 * t218 + t357 * t219;
t180 = -t357 * t220 - t354 * t221;
t179 = -t357 * t218 - t354 * t219;
t176 = -t354 * t208 + t357 * t209;
t175 = t357 * t208 + t354 * t209;
t174 = -pkin(4) * t253 - t355 * t224 + t358 * t241;
t173 = -pkin(4) * t252 - t355 * t223 + t358 * t240;
t172 = t358 * t183 + t370;
t171 = t358 * t182 - t370;
t170 = -pkin(3) * t400 + pkin(6) * t232 + t392;
t169 = t358 * t195 - t355 * t226;
t168 = t358 * t194 - t355 * t230;
t167 = -pkin(3) * t225 + pkin(6) * t209 - t386;
t166 = t358 * t191 + t355 * t400;
t165 = t355 * t191 - t358 * t400;
t164 = t351 * t198 - t206 * t350;
t163 = t350 * t198 + t206 * t351;
t162 = t358 * t176 + t355 * t225;
t161 = t355 * t176 - t358 * t225;
t160 = -pkin(1) * t197 + pkin(2) * t276 - pkin(5) * t207;
t159 = -t354 * t187 + t357 * t189;
t158 = -t354 * t186 + t357 * t188;
t157 = t357 * t187 + t354 * t189;
t156 = -t357 * t186 - t354 * t188;
t153 = -pkin(4) * t197 - (pkin(2) * t355 - pkin(5) * t358) * t206;
t152 = t358 * t158 - t355 * t265;
t151 = t358 * t159 + t355 * t245;
t150 = t355 * t159 - t358 * t245;
t149 = t351 * t166 + t350 * t190;
t148 = t350 * t166 - t351 * t190;
t147 = -pkin(2) * t190 - pkin(3) * t231 + t181;
t146 = -pkin(3) * t239 + pkin(6) * t155;
t145 = -pkin(2) * t175 - pkin(3) * t208 + t178;
t144 = t351 * t162 + t350 * t175;
t143 = t350 * t162 - t351 * t175;
t142 = -pkin(2) * t157 - pkin(3) * t187;
t141 = -pkin(6) * t187 - t154;
t140 = -pkin(5) * t190 - t354 * t170 + t357 * t196;
t139 = -pkin(3) * t245 + pkin(6) * t189 + t155;
t138 = -pkin(5) * t175 - t354 * t167 + t357 * t185;
t137 = t357 * t155 - t390;
t136 = t354 * t155 + t384;
t135 = t351 * t151 + t350 * t157;
t134 = t350 * t151 - t351 * t157;
t133 = t358 * t137 + t355 * t239;
t132 = t355 * t137 - t358 * t239;
t131 = -pkin(1) * t165 + pkin(2) * t400 - pkin(5) * t191 - t357 * t170 - t354 * t196;
t130 = -pkin(1) * t161 + pkin(2) * t225 - pkin(5) * t176 - t357 * t167 - t354 * t185;
t129 = -pkin(2) * t136 - pkin(3) * t154;
t128 = -pkin(4) * t165 + t358 * t140 - t355 * t147;
t127 = -pkin(4) * t161 + t358 * t138 - t355 * t145;
t126 = -pkin(5) * t157 - t354 * t139 + t357 * t141;
t125 = -pkin(5) * t136 - pkin(6) * t384 - t354 * t146;
t124 = t351 * t133 + t350 * t136;
t123 = t350 * t133 - t351 * t136;
t122 = -pkin(1) * t150 + pkin(2) * t245 - pkin(5) * t159 - t357 * t139 - t354 * t141;
t121 = -pkin(4) * t150 + t358 * t126 - t355 * t142;
t120 = -pkin(1) * t132 + pkin(2) * t239 - pkin(5) * t137 + pkin(6) * t390 - t357 * t146;
t119 = -pkin(4) * t132 + t358 * t125 - t355 * t129;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, 0, 0, 0, 0, 0, 0, -t308, -t393, 0, t234, 0, 0, 0, 0, 0, 0, t215, t216, t351 * t273, t164, 0, 0, 0, 0, 0, 0, t144, t149, t135, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t274, 0, 0, 0, 0, 0, 0, -t307, -t396, 0, t233, 0, 0, 0, 0, 0, 0, t213, t214, t350 * t273, t163, 0, 0, 0, 0, 0, 0, t143, t148, t134, t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t380, 0, 0, 0, 0, 0, 0, t322, -t321, 0, -t247, 0, 0, 0, 0, 0, 0, t252, t253, t272, t197, 0, 0, 0, 0, 0, 0, t161, t165, t150, t132; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t403, -t402, -t274, -qJ(1) * t274, 0, 0, t393, 0, -t308, t350 * qJDD(2), t351 * t278 + t365 * t350, t351 * t279 + t364 * t350, t351 * t247, -qJ(1) * t233 - (pkin(1) * t350 - pkin(4) * t351) * t247, t351 * t261 - t350 * t281, t351 * t244 - t350 * t269, t351 * t259 - t350 * t284, t351 * t260 - t350 * t280, t351 * t258 - t350 * t282, t351 * t294, -qJ(1) * t213 + t351 * t173 - t350 * t201, -qJ(1) * t214 + t351 * t174 - t350 * t200, t351 * t199 - t366 * t350, -qJ(1) * t163 + t351 * t153 - t350 * t160, t351 * t172 - t350 * t180, t351 * t152 - t350 * t156, t351 * t168 - t350 * t192, t351 * t171 - t350 * t179, t351 * t169 - t350 * t193, t351 * t202 - t350 * t203, -qJ(1) * t143 + t351 * t127 - t350 * t130, -qJ(1) * t148 + t351 * t128 - t350 * t131, -qJ(1) * t134 + t351 * t121 - t350 * t122, -qJ(1) * t123 + t351 * t119 - t350 * t120; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t402, -t403, t275, qJ(1) * t275, 0, 0, t396, 0, -t307, -t351 * qJDD(2), t350 * t278 - t365 * t351, t350 * t279 - t364 * t351, t350 * t247, qJ(1) * t234 - (-pkin(1) * t351 - pkin(4) * t350) * t247, t350 * t261 + t351 * t281, t350 * t244 + t351 * t269, t350 * t259 + t351 * t284, t350 * t260 + t351 * t280, t350 * t258 + t351 * t282, t350 * t294, qJ(1) * t215 + t350 * t173 + t351 * t201, qJ(1) * t216 + t350 * t174 + t351 * t200, t350 * t199 + t366 * t351, qJ(1) * t164 + t350 * t153 + t351 * t160, t350 * t172 + t351 * t180, t350 * t152 + t351 * t156, t350 * t168 + t351 * t192, t350 * t171 + t351 * t179, t350 * t169 + t351 * t193, t350 * t202 + t351 * t203, qJ(1) * t144 + t350 * t127 + t351 * t130, qJ(1) * t149 + t350 * t128 + t351 * t131, qJ(1) * t135 + t350 * t121 + t351 * t122, qJ(1) * t124 + t350 * t119 + t351 * t120; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t325, t326, 0, 0, 0, 0, t321, 0, t322, 0, -t279, t278, t248, pkin(1) * t325 + pkin(4) * t248, t355 * t293 + t361, t355 * t270 + t358 * t324, t355 * t288 - t354 * t372, t355 * t292 - t361, t355 * t286 - t357 * t372, -t358 * qJDD(3) + t355 * t315, -pkin(1) * t283 + pkin(4) * t255 + t358 * t223 + t355 * t240, -pkin(1) * t285 + pkin(4) * t256 + t358 * t224 + t355 * t241, pkin(4) * t273 + t355 * t206, pkin(4) * t198 - (-pkin(2) * t358 - pkin(5) * t355 - pkin(1)) * t206, t355 * t183 - t369, t355 * t158 + t358 * t265, t355 * t194 + t358 * t230, t355 * t182 + t369, t355 * t195 + t358 * t226, t355 * t204 - t358 * t371, -pkin(1) * t175 + pkin(4) * t162 + t355 * t138 + t358 * t145, -pkin(1) * t190 + pkin(4) * t166 + t355 * t140 + t358 * t147, -pkin(1) * t157 + pkin(4) * t151 + t355 * t126 + t358 * t142, -pkin(1) * t136 + pkin(4) * t133 + t355 * t125 + t358 * t129;];
tauB_reg = t1;
