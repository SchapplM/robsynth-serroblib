% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPR4
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:36
% EndTime: 2019-12-31 17:02:40
% DurationCPUTime: 3.32s
% Computational Cost: add. (12081->322), mult. (17700->486), div. (0->0), fcn. (11435->8), ass. (0->224)
t374 = cos(qJ(2));
t365 = qJDD(1) + qJDD(2);
t371 = sin(qJ(2));
t402 = t371 * t365;
t367 = (qJD(1) + qJD(2));
t418 = t367 ^ 2;
t340 = t418 * t374 + t402;
t398 = t374 * t365;
t401 = t371 * t418;
t343 = -t398 + t401;
t372 = sin(qJ(1));
t375 = cos(qJ(1));
t298 = t340 * t372 + t375 * t343;
t321 = pkin(5) * t340 - t374 * g(3);
t425 = pkin(5) * t343 - g(3) * t371;
t434 = pkin(4) * t298 + t372 * t321 + t375 * t425;
t384 = t375 * t340 - t372 * t343;
t433 = pkin(4) * t384 + t375 * t321 - t372 * t425;
t354 = t375 * g(1) + t372 * g(2);
t377 = qJD(1) ^ 2;
t346 = -t377 * pkin(1) - t354;
t353 = t372 * g(1) - t375 * g(2);
t382 = qJDD(1) * pkin(1) + t353;
t308 = t374 * t346 + t371 * t382;
t427 = -pkin(2) * t418 + qJ(3) * t365 + (2 * qJD(3) * t367) + t308;
t307 = t346 * t371 - t374 * t382;
t388 = t307 * t371 + t374 * t308;
t269 = t374 * t307 - t371 * t308;
t397 = t375 * t269;
t430 = -t372 * t388 + t397;
t400 = t372 * t269;
t228 = t375 * t388 + t400;
t370 = sin(qJ(4));
t368 = sin(pkin(7));
t369 = cos(pkin(7));
t373 = cos(qJ(4));
t426 = t368 * t370 - t369 * t373;
t326 = t426 * t367;
t383 = t368 * t373 + t369 * t370;
t328 = t383 * t367;
t291 = t328 * t326;
t420 = qJDD(4) - t291;
t429 = t370 * t420;
t428 = t373 * t420;
t366 = t368 ^ 2;
t379 = t369 ^ 2;
t419 = t418 * (t366 + t379);
t330 = t369 * t419;
t390 = t369 * t398;
t303 = -t371 * t330 + t390;
t305 = t374 * t330 + t369 * t402;
t264 = t375 * t303 - t372 * t305;
t424 = t372 * t303 + t375 * t305;
t421 = t383 * t365;
t360 = t379 * t418;
t408 = t366 * t418;
t338 = t360 + t408;
t324 = t326 ^ 2;
t325 = t328 ^ 2;
t417 = g(3) * t369;
t361 = t365 * pkin(2);
t250 = -t417 + (pkin(3) * t418 * t369 - pkin(6) * t365 - t427) * t368;
t276 = -g(3) * t368 + t427 * t369;
t409 = t365 * t369;
t251 = -pkin(3) * t360 + pkin(6) * t409 + t276;
t218 = -t373 * t250 + t251 * t370;
t219 = t370 * t250 + t373 * t251;
t186 = -t218 * t373 + t370 * t219;
t414 = t186 * t368;
t413 = t186 * t369;
t280 = -t418 * qJ(3) + qJDD(3) + t307 - t361;
t271 = -pkin(3) * t409 - t338 * pkin(6) + t280;
t412 = t271 * t370;
t411 = t271 * t373;
t283 = qJDD(4) + t291;
t410 = t283 * t373;
t407 = t368 * t369;
t404 = t370 * t283;
t403 = t371 * t280;
t399 = t374 * t280;
t395 = t326 * qJD(4);
t394 = t328 * qJD(4);
t392 = t371 * t291;
t391 = t374 * t291;
t389 = -t280 + t361;
t187 = t218 * t370 + t373 * t219;
t275 = t427 * t368 + t417;
t235 = t275 * t368 + t369 * t276;
t312 = -t353 * t372 - t375 * t354;
t348 = t375 * qJDD(1) - t372 * t377;
t387 = -pkin(4) * t348 - g(3) * t372;
t322 = t426 * t365;
t234 = t275 * t369 - t276 * t368;
t309 = t340 * t407;
t310 = t368 * t390 - t401 * t407;
t386 = t375 * t309 + t372 * t310;
t385 = t372 * t309 - t375 * t310;
t311 = t375 * t353 - t372 * t354;
t376 = qJD(4) ^ 2;
t359 = t379 * t365;
t358 = t366 * t365;
t347 = t372 * qJDD(1) + t375 * t377;
t339 = t360 - t408;
t334 = t359 - t358;
t333 = t359 + t358;
t332 = -pkin(4) * t347 + t375 * g(3);
t329 = t368 * t419;
t315 = -t325 - t376;
t314 = -t325 + t376;
t313 = t324 - t376;
t304 = t374 * t329 + t368 * t402;
t301 = t371 * t329 - t368 * t398;
t295 = t374 * t334 - t371 * t339;
t294 = t374 * t333 - t371 * t338;
t293 = t371 * t334 + t374 * t339;
t292 = t371 * t333 + t374 * t338;
t289 = -t325 + t324;
t288 = t421 - t395;
t287 = t421 - 0.2e1 * t395;
t286 = -t322 - t394;
t285 = t322 + 0.2e1 * t394;
t281 = -t376 - t324;
t278 = (-t326 * t373 + t328 * t370) * qJD(4);
t277 = (-t326 * t370 - t328 * t373) * qJD(4);
t274 = -t324 - t325;
t265 = -t372 * t301 + t375 * t304;
t263 = t375 * t301 + t372 * t304;
t262 = t288 * t373 - t370 * t394;
t261 = t370 * t288 + t373 * t394;
t260 = -t370 * t286 + t373 * t395;
t259 = t286 * t373 + t370 * t395;
t258 = -t370 * t315 - t410;
t257 = -t370 * t314 + t428;
t256 = t313 * t373 - t404;
t255 = t315 * t373 - t404;
t254 = t314 * t373 + t429;
t253 = t370 * t313 + t410;
t252 = pkin(1) * g(3) + pkin(5) * t388;
t246 = -t372 * t292 + t375 * t294;
t245 = t375 * t292 + t372 * t294;
t244 = -t285 * t373 - t370 * t287;
t243 = -t322 * t373 + t370 * t421;
t242 = -t370 * t285 + t287 * t373;
t241 = -t322 * t370 - t373 * t421;
t240 = t281 * t373 - t429;
t239 = t370 * t281 + t428;
t238 = -t277 * t368 + t278 * t369;
t237 = t371 * qJDD(4) + t374 * t238;
t236 = -t374 * qJDD(4) + t371 * t238;
t232 = -pkin(5) * t301 - t371 * t276 + t369 * t399;
t231 = -pkin(5) * t303 - t371 * t275 + t368 * t399;
t230 = pkin(5) * t304 + t374 * t276 + t369 * t403;
t229 = -pkin(5) * t305 + t374 * t275 + t368 * t403;
t226 = -pkin(6) * t255 + t411;
t225 = -t261 * t368 + t262 * t369;
t224 = -t259 * t368 + t260 * t369;
t223 = -t255 * t368 + t258 * t369;
t222 = -t254 * t368 + t257 * t369;
t221 = -t253 * t368 + t256 * t369;
t220 = t255 * t369 + t258 * t368;
t216 = -pkin(6) * t239 + t412;
t215 = -pkin(5) * t292 + t374 * t234;
t214 = pkin(5) * t294 + t234 * t371;
t213 = t374 * t235 + t403;
t212 = t371 * t235 - t399;
t211 = t374 * t222 + t371 * t421;
t210 = t374 * t221 - t371 * t322;
t209 = t371 * t222 - t374 * t421;
t208 = t371 * t221 + t374 * t322;
t207 = -pkin(3) * t287 + pkin(6) * t258 + t412;
t206 = -t242 * t368 + t244 * t369;
t205 = -t241 * t368 + t243 * t369;
t204 = t241 * t369 + t243 * t368;
t203 = -t239 * t368 + t240 * t369;
t202 = t239 * t369 + t240 * t368;
t201 = t374 * t225 + t392;
t200 = t374 * t224 - t392;
t199 = t371 * t225 - t391;
t198 = t371 * t224 + t391;
t197 = t374 * t223 + t371 * t287;
t196 = t371 * t223 - t374 * t287;
t195 = -pkin(3) * t285 + pkin(6) * t240 - t411;
t194 = t374 * t206 - t371 * t289;
t193 = t371 * t206 + t374 * t289;
t192 = t374 * t203 + t371 * t285;
t191 = t371 * t203 - t374 * t285;
t190 = t374 * t205 + t371 * t274;
t189 = t371 * t205 - t374 * t274;
t188 = -pkin(2) * t204 - pkin(3) * t241;
t185 = -t372 * t212 + t375 * t213;
t184 = t375 * t212 + t372 * t213;
t183 = -pkin(2) * t220 - pkin(3) * t255 + t219;
t182 = -pkin(3) * t271 + pkin(6) * t187;
t181 = -t372 * t196 + t375 * t197;
t180 = t375 * t196 + t372 * t197;
t179 = -pkin(6) * t241 - t186;
t178 = -pkin(5) * t212 - (pkin(2) * t371 - qJ(3) * t374) * t234;
t177 = -pkin(2) * t202 - pkin(3) * t239 + t218;
t176 = -pkin(3) * t274 + pkin(6) * t243 + t187;
t175 = -t372 * t191 + t375 * t192;
t174 = t375 * t191 + t372 * t192;
t173 = -qJ(3) * t220 - t207 * t368 + t226 * t369;
t172 = -t372 * t189 + t375 * t190;
t171 = t375 * t189 + t372 * t190;
t170 = pkin(5) * t213 - (-pkin(2) * t374 - qJ(3) * t371 - pkin(1)) * t234;
t169 = -qJ(3) * t202 - t195 * t368 + t216 * t369;
t168 = t187 * t369 - t414;
t167 = t187 * t368 + t413;
t166 = t374 * t168 + t371 * t271;
t165 = t371 * t168 - t374 * t271;
t164 = -pkin(2) * t167 - pkin(3) * t186;
t163 = -pkin(5) * t196 + t374 * t173 - t371 * t183;
t162 = -qJ(3) * t204 - t176 * t368 + t179 * t369;
t161 = -pkin(1) * t220 + pkin(5) * t197 + t371 * t173 + t374 * t183;
t160 = -pkin(5) * t191 + t374 * t169 - t371 * t177;
t159 = -pkin(6) * t413 - qJ(3) * t167 - t182 * t368;
t158 = -pkin(1) * t202 + pkin(5) * t192 + t371 * t169 + t374 * t177;
t157 = -t372 * t165 + t375 * t166;
t156 = t375 * t165 + t372 * t166;
t155 = -pkin(5) * t189 + t374 * t162 - t371 * t188;
t154 = -pkin(1) * t204 + pkin(5) * t190 + t371 * t162 + t374 * t188;
t153 = -pkin(5) * t165 + t374 * t159 - t371 * t164;
t152 = -pkin(1) * t167 + pkin(5) * t166 + t371 * t159 + t374 * t164;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t347, -t348, 0, t312, 0, 0, 0, 0, 0, 0, -t384, t298, 0, t228, 0, 0, 0, 0, 0, 0, -t424, t265, t246, t185, 0, 0, 0, 0, 0, 0, t175, t181, t172, t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t348, -t347, 0, t311, 0, 0, 0, 0, 0, 0, -t298, -t384, 0, -t430, 0, 0, 0, 0, 0, 0, t264, t263, t245, t184, 0, 0, 0, 0, 0, 0, t174, t180, t171, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, 0, 0, 0, 0, 0, 0, t202, t220, t204, t167; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t348, 0, -t347, 0, t387, -t332, -t311, -pkin(4) * t311, 0, 0, -t298, 0, -t384, 0, t434, t433, t430, pkin(4) * t430 + pkin(5) * t397 - t372 * t252, -t385, -t372 * t293 + t375 * t295, t265, t385, t424, 0, -pkin(4) * t264 - t372 * t229 + t375 * t231, -pkin(4) * t263 - t372 * t230 + t375 * t232, -pkin(4) * t245 - t372 * t214 + t375 * t215, -pkin(4) * t184 - t372 * t170 + t375 * t178, -t372 * t199 + t375 * t201, -t372 * t193 + t375 * t194, -t372 * t209 + t375 * t211, -t372 * t198 + t375 * t200, -t372 * t208 + t375 * t210, -t372 * t236 + t375 * t237, -pkin(4) * t174 - t372 * t158 + t375 * t160, -pkin(4) * t180 - t372 * t161 + t375 * t163, -pkin(4) * t171 - t372 * t154 + t375 * t155, -pkin(4) * t156 - t372 * t152 + t375 * t153; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t347, 0, t348, 0, t332, t387, t312, pkin(4) * t312, 0, 0, t384, 0, -t298, 0, -t433, t434, t228, pkin(4) * t228 + pkin(5) * t400 + t375 * t252, t386, t375 * t293 + t372 * t295, t263, -t386, -t264, 0, -pkin(4) * t424 + t375 * t229 + t372 * t231, pkin(4) * t265 + t375 * t230 + t372 * t232, pkin(4) * t246 + t375 * t214 + t372 * t215, pkin(4) * t185 + t375 * t170 + t372 * t178, t375 * t199 + t372 * t201, t375 * t193 + t372 * t194, t375 * t209 + t372 * t211, t375 * t198 + t372 * t200, t375 * t208 + t372 * t210, t375 * t236 + t372 * t237, pkin(4) * t175 + t375 * t158 + t372 * t160, pkin(4) * t181 + t375 * t161 + t372 * t163, pkin(4) * t172 + t375 * t154 + t372 * t155, pkin(4) * t157 + t375 * t152 + t372 * t153; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t353, t354, 0, 0, 0, 0, 0, 0, 0, t365, -pkin(1) * t343 - t307, -pkin(1) * t340 - t308, 0, -pkin(1) * t269, t358, 0.2e1 * t365 * t407, 0, t359, 0, 0, pkin(1) * t303 - qJ(3) * t330 + t389 * t369, pkin(1) * t301 + qJ(3) * t329 - t389 * t368, pkin(1) * t292 + pkin(2) * t338 + qJ(3) * t333 + t235, pkin(1) * t212 - pkin(2) * t280 + qJ(3) * t235, t261 * t369 + t262 * t368, t242 * t369 + t244 * t368, t254 * t369 + t257 * t368, t259 * t369 + t260 * t368, t253 * t369 + t256 * t368, t277 * t369 + t278 * t368, pkin(1) * t191 - pkin(2) * t285 + qJ(3) * t203 + t195 * t369 + t216 * t368, pkin(1) * t196 - pkin(2) * t287 + qJ(3) * t223 + t207 * t369 + t226 * t368, pkin(1) * t189 - pkin(2) * t274 + qJ(3) * t205 + t176 * t369 + t179 * t368, pkin(1) * t165 - pkin(2) * t271 - pkin(6) * t414 + qJ(3) * t168 + t369 * t182;];
tauB_reg = t1;
