% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRP8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:32
% EndTime: 2019-12-31 18:47:39
% DurationCPUTime: 3.60s
% Computational Cost: add. (7633->316), mult. (11041->393), div. (0->0), fcn. (4802->6), ass. (0->216)
t350 = qJD(4) ^ 2;
t339 = -qJD(1) + qJD(3);
t337 = t339 ^ 2;
t344 = sin(qJ(4));
t342 = t344 ^ 2;
t401 = t342 * t337;
t316 = t350 + t401;
t347 = cos(qJ(4));
t321 = t347 * t337 * t344;
t313 = qJDD(4) - t321;
t388 = t347 * t313;
t273 = -t344 * t316 + t388;
t381 = qJD(4) * t339;
t376 = t347 * t381;
t338 = qJDD(1) - qJDD(3);
t395 = t344 * t338;
t302 = 0.2e1 * t376 - t395;
t345 = sin(qJ(3));
t348 = cos(qJ(3));
t234 = t345 * t273 + t348 * t302;
t237 = t348 * t273 - t345 * t302;
t346 = sin(qJ(1));
t349 = cos(qJ(1));
t196 = t349 * t234 - t346 * t237;
t437 = pkin(5) * t196;
t198 = t346 * t234 + t349 * t237;
t436 = pkin(5) * t198;
t396 = t344 * t313;
t267 = t347 * t316 + t396;
t435 = -pkin(6) * t234 + qJ(2) * t267;
t412 = pkin(1) + pkin(2);
t434 = -pkin(6) * t237 + t412 * t267;
t398 = t344 * t302;
t377 = t344 * t381;
t387 = t347 * t338;
t359 = 0.2e1 * t377 + t387;
t420 = t359 * t347;
t254 = t420 + t398;
t343 = t347 ^ 2;
t310 = (t342 - t343) * t337;
t223 = t345 * t254 + t348 * t310;
t226 = t348 * t254 - t345 * t310;
t433 = t349 * t223 - t346 * t226;
t432 = t346 * t223 + t349 * t226;
t400 = t343 * t337;
t318 = -t350 + t400;
t271 = -t347 * t318 + t396;
t384 = t348 * t338;
t247 = t345 * t271 - t347 * t384;
t251 = t348 * t271 + t345 * t387;
t431 = t349 * t247 - t346 * t251;
t430 = t346 * t247 + t349 * t251;
t429 = -pkin(7) * t273 + qJ(2) * t237 - t412 * t234;
t365 = t345 * t337 + t384;
t288 = pkin(6) * t365 + t345 * g(3);
t392 = t345 * t338;
t375 = -t348 * t337 + t392;
t414 = t346 * t375 + t349 * t365;
t415 = -pkin(6) * t375 + t348 * g(3);
t428 = -pkin(5) * t414 + t349 * t288 - t346 * t415;
t258 = t346 * t365 - t349 * t375;
t427 = -pkin(5) * t258 + t346 * t288 + t349 * t415;
t426 = 2 * qJD(5);
t425 = pkin(3) * t267;
t424 = pkin(7) * t267;
t341 = qJDD(1) * qJ(2);
t323 = t349 * g(1) + t346 * g(2);
t364 = (2 * qJD(2) * qJD(1)) - t323;
t357 = t341 + t364;
t413 = qJD(1) ^ 2;
t287 = -t412 * t413 + t357;
t322 = t346 * g(1) - t349 * g(2);
t362 = -qJDD(2) + t322;
t355 = -t413 * qJ(2) - t362;
t351 = -t412 * qJDD(1) + t355;
t239 = t345 * t287 - t348 * t351;
t240 = t348 * t287 + t345 * t351;
t202 = t348 * t239 - t345 * t240;
t203 = t345 * t239 + t348 * t240;
t171 = t349 * t202 + t346 * t203;
t417 = t346 * t202 - t349 * t203;
t228 = -t337 * pkin(3) - t338 * pkin(7) + t240;
t220 = t344 * g(3) + t347 * t228;
t404 = qJ(5) * t344;
t411 = pkin(4) * t347;
t371 = -t404 - t411;
t301 = t371 * t339;
t368 = t347 * t339 * t301 + qJDD(4) * qJ(5) + (qJD(4) * t426) + t220;
t402 = t339 * t344;
t416 = t301 * t402 + qJDD(5);
t319 = -t350 - t400;
t312 = qJDD(4) + t321;
t397 = t344 * t312;
t270 = t347 * t319 - t397;
t233 = t345 * t270 - t348 * t359;
t236 = t348 * t270 + t345 * t359;
t194 = -t349 * t233 + t346 * t236;
t410 = pkin(5) * t194;
t382 = t342 + t343;
t304 = t382 * t338;
t309 = t382 * t337;
t257 = -t345 * t304 + t348 * t309;
t260 = -t348 * t304 - t345 * t309;
t217 = -t349 * t257 + t346 * t260;
t409 = pkin(5) * t217;
t408 = pkin(6) * t257;
t407 = pkin(6) * t260;
t389 = t347 * t312;
t265 = t344 * t319 + t389;
t406 = pkin(7) * t265;
t403 = qJDD(1) * pkin(1);
t227 = t338 * pkin(3) - t337 * pkin(7) + t239;
t399 = t344 * t227;
t391 = t347 * t227;
t219 = -t347 * g(3) + t344 * t228;
t383 = t309 - t350;
t292 = -t413 * pkin(1) + t357;
t296 = -t355 + t403;
t246 = t349 * t292 - t346 * t296;
t280 = -t346 * t322 - t349 * t323;
t374 = t345 * t321;
t373 = t348 * t321;
t206 = -pkin(3) * t265 + t219;
t314 = t346 * qJDD(1) + t349 * t413;
t300 = -pkin(5) * t314 + t349 * g(3);
t315 = t349 * qJDD(1) - t346 * t413;
t299 = pkin(5) * t315 + t346 * g(3);
t370 = pkin(4) * t344 - qJ(5) * t347;
t369 = -pkin(6) * t233 + qJ(2) * t265;
t184 = t347 * t219 - t344 * t220;
t185 = t344 * t219 + t347 * t220;
t245 = t346 * t292 + t349 * t296;
t367 = t347 * t302 - t344 * t359;
t366 = t344 * t318 + t388;
t279 = t349 * t322 - t346 * t323;
t363 = -pkin(6) * t236 + t412 * t265;
t361 = t376 - t395;
t360 = -t377 - t387;
t358 = -pkin(7) * t270 + qJ(2) * t236 - t412 * t233;
t356 = -qJDD(4) * pkin(4) + t219 + t416;
t354 = -t360 * pkin(4) + t227 + (-t361 - t376) * qJ(5);
t353 = -pkin(3) * t309 + pkin(7) * t304 + qJ(2) * t260 - t412 * t257;
t352 = t402 * t426 - t354;
t317 = t350 - t401;
t298 = t370 * t338;
t297 = t382 * t381;
t278 = t345 * qJDD(4) + t348 * t297;
t277 = t348 * qJDD(4) - t345 * t297;
t276 = -t342 * t381 + t347 * t361;
t275 = -t343 * t381 - t344 * t360;
t272 = -t344 * t317 + t389;
t266 = -t347 * t317 - t397;
t252 = t348 * t272 - t344 * t392;
t249 = -t345 * t272 - t344 * t384;
t244 = t348 * t276 - t374;
t243 = t348 * t275 + t374;
t242 = -t345 * t276 - t373;
t241 = -t345 * t275 + t373;
t230 = -t346 * t277 + t349 * t278;
t229 = t349 * t277 + t346 * t278;
t218 = t346 * t257 + t349 * t260;
t216 = pkin(5) * t218;
t215 = -t346 * t249 + t349 * t252;
t214 = t349 * t249 + t346 * t252;
t213 = t391 + t424;
t212 = t399 - t406;
t211 = -t346 * t242 + t349 * t244;
t210 = -t346 * t241 + t349 * t243;
t209 = t349 * t242 + t346 * t244;
t208 = t349 * t241 + t346 * t243;
t207 = t220 + t425;
t205 = t350 * qJ(5) - t356;
t204 = -t350 * pkin(4) + t368;
t197 = t346 * t233 + t349 * t236;
t193 = pkin(6) * t202 + qJ(2) * g(3);
t192 = pkin(5) * t197;
t191 = -pkin(6) * t203 + t412 * g(3);
t190 = t383 * qJ(5) + t356;
t189 = t383 * pkin(4) + t368;
t188 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t402 + t354;
t187 = (-t359 - t377) * pkin(4) + t352;
t186 = -pkin(4) * t377 + qJ(5) * t302 + t352;
t182 = (-t319 - t350) * qJ(5) + (-qJDD(4) - t312) * pkin(4) + t206 + t416;
t181 = -t425 - qJ(5) * t313 + (-t316 + t350) * pkin(4) - t368;
t180 = t348 * t184 - t408;
t179 = -t345 * t184 - t407;
t178 = -pkin(4) * t398 + t347 * t186 - t424;
t177 = -qJ(5) * t420 - t344 * t187 - t406;
t176 = t347 * t204 - t344 * t205;
t175 = t344 * t204 + t347 * t205;
t174 = t348 * t185 + t345 * t227;
t173 = t345 * t185 - t348 * t227;
t170 = -t344 * t189 + t347 * t190;
t169 = -t345 * t207 + t348 * t213 - t435;
t168 = -t345 * t206 + t348 * t212 + t369;
t167 = t348 * t170 + t345 * t298 - t408;
t166 = -t345 * t170 + t348 * t298 - t407;
t165 = -t348 * t207 - t345 * t213 - t434;
t164 = -t348 * t206 - t345 * t212 + t363;
t163 = t348 * t176 + t345 * t188;
t162 = t345 * t176 - t348 * t188;
t161 = -pkin(3) * t175 - pkin(4) * t205 - qJ(5) * t204;
t160 = t346 * t173 + t349 * t174;
t159 = -t349 * t173 + t346 * t174;
t158 = t348 * t177 - t345 * t182 + t369;
t157 = t348 * t178 - t345 * t181 + t435;
t156 = -pkin(7) * t175 + t370 * t188;
t155 = -t345 * t177 - t348 * t182 + t363;
t154 = -t345 * t178 - t348 * t181 + t434;
t153 = -pkin(6) * t173 - (pkin(3) * t345 - pkin(7) * t348 + qJ(2)) * t184;
t152 = t346 * t162 + t349 * t163;
t151 = -t349 * t162 + t346 * t163;
t150 = -pkin(6) * t174 - (pkin(3) * t348 + pkin(7) * t345 + t412) * t184;
t149 = -pkin(6) * t162 + qJ(2) * t175 + t348 * t156 - t345 * t161;
t148 = -pkin(6) * t163 - t345 * t156 - t348 * t161 + t412 * t175;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t314, -t315, 0, t280, 0, 0, 0, 0, 0, 0, -t314, 0, t315, t246, 0, 0, 0, 0, 0, 0, -t258, t414, 0, -t417, 0, 0, 0, 0, 0, 0, t197, -t198, t218, t160, 0, 0, 0, 0, 0, 0, t197, t218, t198, t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t315, -t314, 0, t279, 0, 0, 0, 0, 0, 0, t315, 0, t314, t245, 0, 0, 0, 0, 0, 0, t414, t258, 0, t171, 0, 0, 0, 0, 0, 0, t194, t196, t217, t159, 0, 0, 0, 0, 0, 0, t194, t217, -t196, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t265, t267, 0, t184, 0, 0, 0, 0, 0, 0, -t265, 0, -t267, -t175; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t315, 0, -t314, 0, -t299, -t300, -t279, -pkin(5) * t279, 0, t315, 0, 0, t314, 0, -t299, -t245, t300, -pkin(5) * t245 + (-pkin(1) * t346 + qJ(2) * t349) * g(3), 0, 0, -t414, 0, -t258, 0, t428, t427, t171, -pkin(5) * t171 - t346 * t191 + t349 * t193, t211, -t432, t215, t210, -t430, t230, -t346 * t164 + t349 * t168 - t410, -t346 * t165 + t349 * t169 - t437, -t346 * t179 + t349 * t180 - t409, -pkin(5) * t159 - t346 * t150 + t349 * t153, t211, t215, t432, t230, t430, t210, -t346 * t155 + t349 * t158 - t410, -t346 * t166 + t349 * t167 - t409, -t346 * t154 + t349 * t157 + t437, -pkin(5) * t151 - t346 * t148 + t349 * t149; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t314, 0, t315, 0, t300, -t299, t280, pkin(5) * t280, 0, t314, 0, 0, -t315, 0, t300, t246, t299, pkin(5) * t246 + (pkin(1) * t349 + qJ(2) * t346) * g(3), 0, 0, -t258, 0, t414, 0, t427, -t428, t417, -pkin(5) * t417 + t349 * t191 + t346 * t193, t209, t433, t214, t208, t431, t229, t349 * t164 + t346 * t168 + t192, t349 * t165 + t346 * t169 - t436, t349 * t179 + t346 * t180 + t216, pkin(5) * t160 + t349 * t150 + t346 * t153, t209, t214, -t433, t229, -t431, t208, t349 * t155 + t346 * t158 + t192, t349 * t166 + t346 * t167 + t216, t349 * t154 + t346 * t157 + t436, pkin(5) * t152 + t349 * t148 + t346 * t149; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t322, t323, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t362 + 0.2e1 * t403, 0, 0.2e1 * t341 + t364, pkin(1) * t296 + qJ(2) * t292, 0, 0, 0, 0, 0, t338, qJ(2) * t375 + t365 * t412 + t239, qJ(2) * t365 - t375 * t412 + t240, 0, qJ(2) * t203 + t202 * t412, -t398, -t367, t266, t420, -t366, 0, pkin(3) * t359 + t358 + t391, pkin(3) * t302 - t399 - t429, -t185 + t353, pkin(3) * t227 - pkin(7) * t185 + qJ(2) * t174 - t412 * t173, -t398, t266, t367, 0, t366, t420, -t347 * t187 - (-pkin(3) - t404) * t359 + t358, -t347 * t189 - t344 * t190 + t353, -t344 * t186 + (-pkin(3) - t411) * t302 + t429, -pkin(7) * t176 + qJ(2) * t163 - t412 * t162 + (pkin(3) - t371) * t188;];
tauB_reg = t1;
