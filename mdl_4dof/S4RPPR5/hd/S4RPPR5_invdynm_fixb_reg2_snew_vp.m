% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RPPR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:52
% EndTime: 2019-12-31 16:39:54
% DurationCPUTime: 2.13s
% Computational Cost: add. (4950->275), mult. (8805->307), div. (0->0), fcn. (3880->6), ass. (0->158)
t357 = sin(pkin(6));
t358 = cos(pkin(6));
t365 = qJD(1) ^ 2;
t329 = -qJDD(1) * t357 + t358 * t365;
t356 = g(3) + qJDD(3);
t313 = qJ(3) * t329 + t356 * t358;
t361 = sin(qJ(1));
t363 = cos(qJ(1));
t330 = qJDD(1) * t358 + t357 * t365;
t381 = qJ(3) * t330 + t356 * t357;
t412 = t363 * t329 + t361 * t330;
t420 = -pkin(4) * t412 + t363 * t313 + t361 * t381;
t391 = (qJD(2) * qJD(1));
t352 = 2 * t391;
t339 = t363 * g(1) + t361 * g(2);
t353 = qJDD(1) * qJ(2);
t376 = -t339 + t353;
t410 = pkin(1) + pkin(2);
t308 = -t365 * t410 + t352 + t376;
t338 = t361 * g(1) - t363 * g(2);
t375 = -qJDD(2) + t338;
t367 = -t365 * qJ(2) - t375;
t366 = -qJDD(1) * t410 + t367;
t269 = t357 * t308 - t358 * t366;
t270 = t308 * t358 + t357 * t366;
t254 = t269 * t358 - t270 * t357;
t380 = t269 * t357 + t270 * t358;
t419 = t254 * t361 - t363 * t380;
t418 = t254 * t363 + t361 * t380;
t288 = -t361 * t329 + t363 * t330;
t417 = pkin(4) * t288 + t361 * t313 - t363 * t381;
t268 = -pkin(3) * t365 - qJDD(1) * pkin(5) + t270;
t360 = sin(qJ(4));
t362 = cos(qJ(4));
t265 = t268 * t360 - t356 * t362;
t266 = t268 * t362 + t356 * t360;
t248 = t265 * t362 - t266 * t360;
t409 = pkin(3) * t248;
t407 = pkin(5) * t358;
t249 = t360 * t265 + t266 * t362;
t267 = qJDD(1) * pkin(3) - pkin(5) * t365 + t269;
t243 = t249 * t358 + t267 * t357;
t406 = qJ(3) * t243;
t405 = qJ(3) * t254;
t404 = qJ(3) * t380;
t403 = qJDD(1) * pkin(1);
t354 = t360 ^ 2;
t402 = t354 * t365;
t401 = t360 * t267;
t344 = t362 * t365 * t360;
t336 = qJDD(4) + t344;
t400 = t360 * t336;
t337 = qJDD(4) - t344;
t399 = t360 * t337;
t396 = t362 * t267;
t395 = t362 * t336;
t394 = t362 * t337;
t355 = t362 ^ 2;
t393 = -t354 - t355;
t392 = qJD(1) * qJD(4);
t390 = t360 * qJDD(1);
t389 = t362 * qJDD(1);
t347 = t360 * t392;
t388 = t362 * t392;
t242 = t249 * t357 - t267 * t358;
t387 = -qJ(3) * t242 - t357 * t409;
t315 = pkin(1) * t365 - t376 - (2 * t391);
t316 = -t367 + t403;
t386 = -t315 * t363 - t316 * t361;
t385 = -t338 * t361 - t339 * t363;
t384 = t357 * t344;
t383 = t358 * t344;
t382 = pkin(3) * t267 - pkin(5) * t249;
t332 = qJDD(1) * t361 + t363 * t365;
t318 = -pkin(4) * t332 + g(3) * t363;
t333 = qJDD(1) * t363 - t361 * t365;
t317 = pkin(4) * t333 + g(3) * t361;
t379 = t315 * t361 - t316 * t363;
t378 = t338 * t363 - t339 * t361;
t377 = -pkin(3) * t358 - pkin(5) * t357 - pkin(2);
t351 = t355 * t365;
t364 = qJD(4) ^ 2;
t343 = -t351 - t364;
t296 = t343 * t360 + t395;
t256 = -pkin(3) * t296 + t265;
t262 = -pkin(5) * t296 + t401;
t300 = t343 * t362 - t400;
t327 = 0.2e1 * t347 - t389;
t273 = t300 * t357 + t327 * t358;
t374 = -qJ(3) * t273 - t256 * t357 + t262 * t358;
t341 = -t364 - t402;
t298 = t341 * t362 - t399;
t257 = -pkin(3) * t298 + t266;
t263 = -pkin(5) * t298 + t396;
t302 = -t341 * t360 - t394;
t325 = 0.2e1 * t388 + t390;
t274 = t302 * t357 + t325 * t358;
t373 = -qJ(3) * t274 - t257 * t357 + t263 * t358;
t372 = pkin(3) * t327 + pkin(5) * t300 - t396;
t371 = pkin(3) * t325 + pkin(5) * t302 + t401;
t275 = t300 * t358 - t327 * t357;
t370 = -qJ(3) * t275 - t358 * t256 - t357 * t262;
t276 = t302 * t358 - t325 * t357;
t369 = -qJ(3) * t276 - t358 * t257 - t357 * t263;
t331 = t393 * qJDD(1);
t334 = t351 + t402;
t368 = pkin(3) * t334 + pkin(5) * t331 + t249;
t342 = t351 - t364;
t340 = t364 - t402;
t335 = -t351 + t402;
t328 = t347 - t389;
t326 = -t388 - t390;
t323 = t393 * t392;
t322 = t375 + 0.2e1 * t403;
t319 = -t339 + t352 + 0.2e1 * t353;
t310 = t326 * t362 + t354 * t392;
t309 = -t328 * t360 + t355 * t392;
t306 = qJDD(4) * t357 + t323 * t358;
t305 = -qJDD(4) * t358 + t323 * t357;
t301 = -t340 * t360 + t395;
t299 = t342 * t362 - t399;
t297 = t340 * t362 + t400;
t295 = t342 * t360 + t394;
t294 = (t326 - t388) * t360;
t293 = (t328 + t347) * t362;
t291 = t331 * t358 - t334 * t357;
t290 = t331 * t357 + t334 * t358;
t287 = t325 * t360 + t327 * t362;
t286 = -t325 * t362 + t327 * t360;
t285 = t310 * t358 - t384;
t284 = t309 * t358 + t384;
t283 = t310 * t357 + t383;
t282 = t309 * t357 - t383;
t281 = t301 * t358 - t357 * t390;
t280 = t299 * t358 - t357 * t389;
t279 = t301 * t357 + t358 * t390;
t278 = t299 * t357 + t358 * t389;
t277 = pkin(1) * t316 - qJ(2) * t315;
t272 = t287 * t358 + t335 * t357;
t271 = t287 * t357 - t335 * t358;
t261 = -qJ(2) * t329 + t330 * t410 + t269;
t260 = qJ(2) * t330 + t329 * t410 + t270;
t251 = qJ(2) * t356 + t405;
t250 = t356 * t410 - t404;
t245 = -qJ(3) * t290 + t248 * t358;
t244 = qJ(3) * t291 + t248 * t357;
t241 = qJ(2) * t276 - t274 * t410 - t371;
t240 = qJ(2) * t275 - t273 * t410 - t372;
t239 = qJ(2) * t291 - t290 * t410 - t368;
t238 = qJ(2) * t298 + t373;
t237 = qJ(2) * t296 + t374;
t236 = t298 * t410 + t369;
t235 = t296 * t410 + t370;
t234 = qJ(2) * t380 + t254 * t410;
t233 = -(qJ(2) - t407) * t248 + t387;
t232 = -t406 - (pkin(1) - t377) * t248;
t231 = qJ(2) * t243 - t242 * t410 + t382;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t333, 0, -t332, 0, -t317, -t318, -t378, -pkin(4) * t378, 0, t333, 0, 0, t332, 0, -t317, t379, t318, pkin(4) * t379 + (-pkin(1) * t361 + qJ(2) * t363) * g(3), 0, 0, -t288, 0, -t412, 0, -t417, t420, t418, -pkin(4) * t418 - t361 * t250 + t363 * t251, t283 * t361 + t285 * t363, t271 * t361 + t272 * t363, t279 * t361 + t281 * t363, t282 * t361 + t284 * t363, t278 * t361 + t280 * t363, t305 * t361 + t306 * t363, t363 * t237 - t361 * t235 - pkin(4) * (-t273 * t363 + t275 * t361), t363 * t238 - t361 * t236 - pkin(4) * (-t274 * t363 + t276 * t361), t363 * t245 + t361 * t244 - pkin(4) * (-t290 * t363 + t291 * t361), t363 * t233 - t361 * t232 - pkin(4) * (-t242 * t363 + t243 * t361); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t332, 0, t333, 0, t318, -t317, t385, pkin(4) * t385, 0, t332, 0, 0, -t333, 0, t318, t386, t317, pkin(4) * t386 + (pkin(1) * t363 + qJ(2) * t361) * g(3), 0, 0, -t412, 0, t288, 0, t420, t417, t419, -pkin(4) * t419 + t363 * t250 + t361 * t251, -t283 * t363 + t285 * t361, -t271 * t363 + t272 * t361, -t279 * t363 + t281 * t361, -t282 * t363 + t284 * t361, -t278 * t363 + t280 * t361, -t305 * t363 + t306 * t361, t361 * t237 + t363 * t235 + pkin(4) * (t273 * t361 + t275 * t363), t361 * t238 + t363 * t236 + pkin(4) * (t274 * t361 + t276 * t363), t361 * t245 - t363 * t244 + pkin(4) * (t290 * t361 + t291 * t363), t361 * t233 + t363 * t232 + pkin(4) * (t242 * t361 + t243 * t363); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t338, t339, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t322, 0, t319, t277, 0, 0, 0, 0, 0, qJDD(1), t261, t260, 0, t234, -t294, -t286, -t297, -t293, -t295, 0, t240, t241, t239, t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t365, 0, 0, -g(3), -t338, 0, 0, qJDD(1), 0, 0, t365, 0, 0, -t316, g(3), qJ(2) * g(3), 0, 0, -t330, 0, -t329, 0, t381, t313, t254, t251, t285, t272, t281, t284, t280, t306, t237, t238, t245, t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t365, 0, qJDD(1), 0, g(3), 0, -t339, 0, 0, t365, 0, 0, -qJDD(1), 0, g(3), -t315, 0, pkin(1) * g(3), 0, 0, -t329, 0, t330, 0, t313, -t381, -t380, t250, -t283, -t271, -t279, -t282, -t278, -t305, t235, t236, -t244, t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t338, t339, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t322, 0, t319, t277, 0, 0, 0, 0, 0, qJDD(1), t261, t260, 0, t234, -t294, -t286, -t297, -t293, -t295, 0, t240, t241, t239, t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t365, 0, 0, -t316, g(3), 0, 0, 0, -t330, 0, -t329, 0, t381, t313, t254, t405, t285, t272, t281, t284, t280, t306, t374, t373, t245, t248 * t407 + t387; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t316, 0, -t315, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(2) * t330 + t269, pkin(2) * t329 + t270, 0, pkin(2) * t254, -t294, -t286, -t297, -t293, -t295, 0, -pkin(2) * t273 - t372, -pkin(2) * t274 - t371, -pkin(2) * t290 - t368, -pkin(2) * t242 + t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t365, 0, 0, qJDD(1), 0, -g(3), t315, 0, 0, 0, 0, t329, 0, -t330, 0, -t313, t381, t380, -pkin(2) * t356 + t404, t283, t271, t279, t282, t278, t305, -pkin(2) * t296 - t370, -pkin(2) * t298 - t369, t244, -t248 * t377 + t406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t365, 0, 0, t356, t269, 0, t310, t287, t301, t309, t299, t323, t262, t263, t248, pkin(5) * t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t365, 0, -qJDD(1), 0, -t356, 0, t270, 0, t344, -t335, t390, -t344, t389, -qJDD(4), t256, t257, 0, t409; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t269, -t270, 0, 0, t294, t286, t297, t293, t295, 0, t372, t371, t368, -t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t326, t327, t336, t388, t342, -t388, 0, t267, t265, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t347, -t325, t340, t328, t337, t347, -t267, 0, t266, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t344, t335, -t390, t344, -t389, qJDD(4), -t265, -t266, 0, 0;];
m_new_reg = t1;
