% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRPP5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:44
% EndTime: 2019-12-31 18:16:49
% DurationCPUTime: 2.56s
% Computational Cost: add. (3421->315), mult. (7231->342), div. (0->0), fcn. (3194->4), ass. (0->203)
t317 = qJD(3) ^ 2;
t315 = cos(qJ(3));
t311 = t315 ^ 2;
t318 = qJD(1) ^ 2;
t382 = t311 * t318;
t291 = t317 + t382;
t313 = sin(qJ(3));
t353 = t313 * t318 * t315;
t282 = qJDD(3) + t353;
t379 = t313 * t282;
t234 = t315 * t291 + t379;
t363 = qJD(1) * qJD(3);
t350 = t313 * t363;
t358 = t315 * qJDD(1);
t272 = -0.2e1 * t350 + t358;
t314 = sin(qJ(1));
t316 = cos(qJ(1));
t196 = t316 * t234 + t314 * t272;
t392 = pkin(5) * t196;
t201 = t314 * t234 - t316 * t272;
t193 = pkin(5) * t201;
t372 = t315 * t282;
t241 = -t313 * t291 + t372;
t398 = -pkin(6) - pkin(1);
t354 = t398 * t241;
t412 = pkin(2) * t272 - t354;
t340 = -qJ(2) * t272 + t234 * t398;
t310 = t313 ^ 2;
t365 = t310 + t311;
t276 = t365 * qJDD(1);
t279 = t365 * t318;
t339 = -qJ(2) * t279 - t398 * t276;
t369 = t316 * t276;
t219 = -t314 * t279 + t369;
t391 = pkin(5) * t219;
t376 = t314 * t276;
t220 = t316 * t279 + t376;
t215 = pkin(5) * t220;
t347 = t315 * t363;
t360 = t313 * qJDD(1);
t269 = 0.2e1 * t347 + t360;
t212 = t313 * t269 - t315 * t272;
t280 = (-t310 + t311) * t318;
t191 = t314 * t212 - t316 * t280;
t411 = t316 * t212 + t314 * t280;
t383 = t310 * t318;
t289 = -t317 + t383;
t233 = t313 * t289 + t372;
t357 = t316 * qJDD(1);
t205 = t314 * t233 - t313 * t357;
t359 = t314 * qJDD(1);
t410 = t316 * t233 + t313 * t359;
t292 = -t317 + t382;
t283 = qJDD(3) - t353;
t378 = t313 * t283;
t236 = -t315 * t292 + t378;
t206 = t314 * t236 + t315 * t357;
t409 = t316 * t236 - t314 * t358;
t408 = -pkin(2) * t234 + qJ(2) * t241;
t407 = 2 * qJD(1);
t406 = (t279 - t317) * pkin(3);
t290 = -t317 - t383;
t371 = t315 * t283;
t232 = t313 * t290 + t371;
t239 = t315 * t290 - t378;
t405 = pkin(2) * t232 - qJ(2) * t239;
t270 = t347 + t360;
t364 = qJD(1) * t315;
t281 = -qJD(3) * pkin(4) - qJ(5) * t364;
t404 = t270 * qJ(5) + qJD(3) * t281;
t271 = -t350 + t358;
t390 = t270 * pkin(3);
t403 = (t271 + t272) * qJ(4) - t390;
t285 = t314 * g(1) - t316 * g(2);
t336 = qJDD(2) - t285;
t328 = -t318 * qJ(2) + t336;
t252 = t398 * qJDD(1) + t328;
t366 = t315 * g(3) - t313 * t252;
t402 = -qJDD(3) * qJ(4) - 0.2e1 * qJD(4) * qJD(3) + t366;
t356 = qJD(5) * t407;
t401 = t313 * t356 + t404;
t400 = pkin(3) * t283 + qJ(4) * t290 + t405;
t399 = pkin(3) + pkin(4);
t396 = pkin(2) * t276;
t395 = pkin(2) * t279;
t394 = pkin(3) * t315;
t195 = -t316 * t232 + t314 * t269;
t393 = pkin(5) * t195;
t388 = qJ(4) * t279;
t387 = qJ(4) * t313;
t386 = qJDD(1) * pkin(1);
t385 = qJDD(3) * pkin(4);
t384 = t271 * qJ(4);
t286 = t316 * g(1) + t314 * g(2);
t308 = qJDD(1) * qJ(2);
t331 = t286 - t308;
t326 = -t398 * t318 + t331;
t362 = qJD(2) * qJD(1);
t247 = t326 - (2 * t362);
t381 = t313 * t247;
t380 = t313 * t272;
t374 = t315 * t247;
t373 = t315 * t269;
t216 = t313 * g(3) + t315 * t252;
t355 = t398 * t239;
t351 = (-t269 - t270) * pkin(3);
t303 = 2 * t362;
t253 = -t318 * pkin(1) + t303 - t331;
t254 = -t328 + t386;
t198 = t316 * t253 - t314 * t254;
t228 = -t314 * t285 - t316 * t286;
t343 = t314 * t353;
t342 = t316 * t353;
t341 = qJ(2) * t269 + t398 * t232;
t277 = -t314 * t318 + t357;
t338 = pkin(5) * t277 + t314 * g(3);
t278 = t316 * t318 + t359;
t337 = -pkin(5) * t278 + t316 * g(3);
t335 = -t387 - t394;
t334 = pkin(3) * t313 - qJ(4) * t315;
t333 = t271 + t350;
t186 = t315 * t216 - t313 * t366;
t187 = -t313 * t216 - t315 * t366;
t194 = t314 * t253 + t316 * t254;
t213 = t373 + t380;
t332 = -t315 * t289 + t379;
t242 = t313 * t292 + t371;
t227 = t316 * t285 - t314 * t286;
t330 = pkin(2) * t269 + t355;
t266 = t334 * qJD(1);
t329 = -qJDD(3) * pkin(3) - t317 * qJ(4) + t266 * t364 + qJDD(4) - t216;
t327 = -t313 * qJD(1) * t266 - t402;
t188 = -t317 * pkin(3) + t327;
t325 = -pkin(3) * t347 - qJ(4) * t350 + t326;
t324 = -pkin(4) * t353 + t315 * t356 - t329;
t164 = pkin(3) * t291 + qJ(4) * t282 + t188 - t408;
t323 = t333 * qJ(5) + t324;
t322 = (qJD(4) * t315 - qJD(2)) * t407 + t325;
t321 = t322 + t384;
t320 = qJDD(5) + (-(2 * qJD(2)) + (0.2e1 * qJD(4) + t281) * t315) * qJD(1) - t270 * pkin(4) - qJ(5) * t383 + t325;
t319 = t320 + t384;
t264 = t365 * t363;
t251 = t316 * qJDD(3) - t314 * t264;
t250 = t314 * qJDD(3) + t316 * t264;
t249 = -t313 * t271 - t311 * t363;
t248 = t315 * t270 - t310 * t363;
t240 = (t271 - t350) * t315;
t230 = (t270 + t347) * t313;
t229 = t335 * qJDD(1) - t396;
t226 = -qJ(4) * t269 + qJ(5) * t283;
t214 = t396 + (t399 * t315 + t387) * qJDD(1);
t210 = -t314 * t248 - t342;
t209 = -t314 * t249 + t342;
t208 = t316 * t248 - t343;
t207 = t316 * t249 + t343;
t202 = -qJ(5) * t282 + t399 * t272;
t199 = t314 * t232 + t316 * t269;
t192 = pkin(5) * t199;
t185 = t329 + t388;
t184 = t327 + t406;
t183 = t321 - t390;
t182 = -t187 - t395;
t181 = t366 + t408;
t180 = t216 + t405;
t179 = t351 + t321;
t178 = t322 + t403;
t177 = t330 - t374;
t176 = t381 + t412;
t175 = t323 + t385;
t174 = -pkin(4) * t383 + t188 + t401;
t173 = t314 * t186 - t316 * t247;
t172 = -t316 * t186 - t314 * t247;
t171 = t319 - t390;
t170 = t385 - t388 + (t333 + t358) * qJ(5) + t324;
t169 = t315 * t188 + t313 * t329;
t168 = t313 * t188 - t315 * t329;
t167 = (-t279 + t383) * pkin(4) - t406 + (-qJ(5) * qJDD(1) + (-(2 * qJD(5)) + t266) * qJD(1)) * t313 + t402 - t404;
t166 = -t329 + t400;
t165 = qJ(5) * t291 + t320 + t403;
t163 = pkin(2) * t186 - qJ(2) * t187;
t162 = -pkin(4) * t269 - qJ(5) * t290 + t319 + t351;
t161 = -t315 * t184 - t313 * t185 - t395;
t160 = -pkin(2) * t247 + t398 * t187;
t159 = -t313 * t178 + (-pkin(2) - t394) * t272 + t354;
t158 = -t315 * t179 + (pkin(2) + t387) * t269 + t355;
t157 = (t291 - t383) * pkin(4) + t164 + t401;
t156 = (qJDD(3) + t283) * pkin(4) + t323 + t400;
t155 = t315 * t174 - t313 * t175;
t154 = t313 * t174 + t315 * t175;
t153 = t314 * t168 - t316 * t183;
t152 = -t316 * t168 - t314 * t183;
t151 = qJ(4) * t171 + qJ(5) * t175;
t150 = -t315 * t162 - t313 * t226 + t330;
t149 = -t313 * t165 - t315 * t202 - t412;
t148 = -t315 * t167 - t313 * t170 + t395;
t147 = t314 * t154 - t316 * t171;
t146 = -t316 * t154 - t314 * t171;
t145 = -qJ(5) * t174 + t399 * t171;
t144 = pkin(2) * t168 - pkin(3) * t329 - qJ(2) * t169 + qJ(4) * t188;
t143 = t398 * t169 + (-pkin(2) + t335) * t183;
t142 = pkin(2) * t154 - qJ(2) * t155 + qJ(4) * t174 + t399 * t175;
t141 = -pkin(2) * t171 - t315 * t145 - t313 * t151 + t398 * t155;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t278, -t277, 0, t228, 0, 0, 0, 0, 0, 0, 0, t278, t277, t198, 0, 0, 0, 0, 0, 0, t199, -t201, -t220, t173, 0, 0, 0, 0, 0, 0, t199, -t220, t201, t153, 0, 0, 0, 0, 0, 0, t199, t201, t220, t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t277, -t278, 0, t227, 0, 0, 0, 0, 0, 0, 0, -t277, t278, t194, 0, 0, 0, 0, 0, 0, t195, t196, t219, t172, 0, 0, 0, 0, 0, 0, t195, t219, -t196, t152, 0, 0, 0, 0, 0, 0, t195, -t196, -t219, t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t239, -t241, 0, t187, 0, 0, 0, 0, 0, 0, t239, 0, t241, t169, 0, 0, 0, 0, 0, 0, t239, t241, 0, t155; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t277, 0, -t278, 0, -t338, -t337, -t227, -pkin(5) * t227, 0, -t277, t278, 0, 0, 0, -t194, t338, t337, -pkin(5) * t194 + (-pkin(1) * t314 + qJ(2) * t316) * g(3), t209, -t191, t206, t210, t205, t251, -t314 * t177 + t316 * t180 - t393, -t314 * t176 + t316 * t181 - t392, -pkin(2) * t369 - t314 * t182 - t391, -pkin(5) * t172 - t314 * t160 + t316 * t163, t209, t206, t191, t251, -t205, t210, -t314 * t158 + t316 * t166 - t393, -t314 * t161 + t316 * t229 - t391, -t314 * t159 + t316 * t164 + t392, -pkin(5) * t152 - t314 * t143 + t316 * t144, t209, t191, -t206, t210, t205, t251, -t314 * t150 + t316 * t156 - t393, -t314 * t149 + t316 * t157 + t392, -t314 * t148 + t316 * t214 + t391, -pkin(5) * t146 - t314 * t141 + t316 * t142; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t278, 0, t277, 0, t337, -t338, t228, pkin(5) * t228, 0, -t278, -t277, 0, 0, 0, t198, -t337, t338, pkin(5) * t198 + (pkin(1) * t316 + qJ(2) * t314) * g(3), t207, t411, -t409, t208, -t410, t250, t316 * t177 + t314 * t180 + t192, t316 * t176 + t314 * t181 - t193, -pkin(2) * t376 + t316 * t182 - t215, pkin(5) * t173 + t316 * t160 + t314 * t163, t207, -t409, -t411, t250, t410, t208, t316 * t158 + t314 * t166 + t192, t316 * t161 + t314 * t229 - t215, t316 * t159 + t314 * t164 + t193, pkin(5) * t153 + t316 * t143 + t314 * t144, t207, -t411, t409, t208, -t410, t250, t316 * t150 + t314 * t156 + t192, t316 * t149 + t314 * t157 + t193, t316 * t148 + t314 * t214 + t215, pkin(5) * t147 + t316 * t141 + t314 * t142; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t285, t286, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t336 - 0.2e1 * t386, -t286 + t303 + 0.2e1 * t308, pkin(1) * t254 + qJ(2) * t253, t240, -t213, t242, t230, -t332, 0, t341 - t381, -t340 - t374, -t186 + t339, -qJ(2) * t247 + t398 * t186, t240, t242, t213, 0, t332, t230, -qJ(4) * t373 - t313 * t179 + t341, -t313 * t184 + t315 * t185 + t339, -pkin(3) * t380 + t315 * t178 + t340, t398 * t168 + (-qJ(2) - t334) * t183, t240, t213, -t242, t230, -t332, 0, -t313 * t162 + t315 * t226 + t341, t315 * t165 - t313 * t202 + t340, -t313 * t167 + t315 * t170 - t339, -qJ(2) * t171 - t313 * t145 + t315 * t151 + t398 * t154;];
tauB_reg = t1;
