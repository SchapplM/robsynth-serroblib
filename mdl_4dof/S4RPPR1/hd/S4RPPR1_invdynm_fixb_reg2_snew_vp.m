% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RPPR1
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RPPR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:08:34
% EndTime: 2019-05-04 19:08:36
% DurationCPUTime: 1.65s
% Computational Cost: add. (4421->171), mult. (6991->172), div. (0->0), fcn. (3592->6), ass. (0->103)
t330 = -qJD(1) + qJD(4);
t327 = t330 ^ 2;
t328 = qJDD(1) - qJDD(4);
t336 = sin(qJ(4));
t338 = cos(qJ(4));
t305 = t327 * t338 - t328 * t336;
t334 = sin(pkin(6));
t335 = cos(pkin(6));
t347 = t327 * t336 + t338 * t328;
t276 = t305 * t334 - t335 * t347;
t332 = g(3) - qJDD(2);
t288 = pkin(5) * t305 + t332 * t338;
t384 = pkin(5) * t347 + t332 * t336;
t243 = qJ(2) * t276 - t288 * t334 + t335 * t384;
t337 = sin(qJ(1));
t339 = cos(qJ(1));
t274 = t335 * t305 + t334 * t347;
t395 = t274 * t337 + t276 * t339;
t407 = -qJ(2) * t274 + t288 * t335 + t334 * t384;
t415 = pkin(4) * t395 + t243 * t339 - t337 * t407;
t408 = t274 * t339 - t276 * t337;
t414 = -pkin(4) * t408 + t243 * t337 + t339 * t407;
t340 = qJD(1) ^ 2;
t315 = qJDD(1) * t334 + t335 * t340;
t316 = qJDD(1) * t335 - t334 * t340;
t283 = t339 * t315 + t316 * t337;
t290 = qJ(2) * t316 + t332 * t334;
t360 = -qJ(2) * t315 + t335 * t332;
t252 = -pkin(4) * t283 - t290 * t337 + t339 * t360;
t320 = g(1) * t339 + g(2) * t337;
t314 = -pkin(1) * t340 - t320;
t319 = g(1) * t337 - t339 * g(2);
t346 = qJDD(1) * pkin(1) + t319;
t281 = t335 * t314 + t334 * t346;
t325 = 2 * qJD(3) * qJD(1);
t361 = t325 + t281;
t362 = qJDD(1) * qJ(3);
t344 = t361 + t362;
t381 = pkin(2) + pkin(3);
t265 = -t381 * t340 + t344;
t333 = qJDD(1) * pkin(2);
t280 = t314 * t334 - t335 * t346;
t350 = -qJDD(3) - t280;
t272 = -t340 * qJ(3) - t333 - t350;
t269 = -qJDD(1) * pkin(3) + t272;
t239 = t265 * t336 - t338 * t269;
t240 = t338 * t265 + t336 * t269;
t232 = t239 * t338 - t240 * t336;
t349 = t239 * t336 + t240 * t338;
t222 = t232 * t334 - t335 * t349;
t398 = t232 * t335 + t334 * t349;
t410 = -t222 * t337 + t339 * t398;
t409 = t222 * t339 + t337 * t398;
t358 = t280 * t334 + t335 * t281;
t257 = t280 * t335 - t281 * t334;
t372 = t257 * t339;
t404 = -t337 * t358 + t372;
t373 = t257 * t337;
t403 = t339 * t358 + t373;
t356 = -t315 * t337 + t339 * t316;
t382 = pkin(4) * t356 + t290 * t339 + t337 * t360;
t270 = -pkin(2) * t340 + t344;
t246 = t270 * t334 - t272 * t335;
t359 = t335 * t270 + t272 * t334;
t397 = -t246 * t337 + t339 * t359;
t396 = t246 * t339 + t337 * t359;
t380 = pkin(1) * t332;
t379 = pkin(3) * t232;
t377 = pkin(5) * t232;
t376 = pkin(5) * t349;
t374 = qJ(3) * t332;
t363 = -pkin(2) * t272 + qJ(3) * t270;
t354 = -t319 * t337 - t339 * t320;
t266 = -pkin(1) * t315 - t281;
t353 = pkin(2) * t232 + qJ(3) * t349 + t379;
t352 = pkin(3) * t305 + t240;
t318 = qJDD(1) * t339 - t337 * t340;
t351 = -pkin(4) * t318 - g(3) * t337;
t348 = t319 * t339 - t320 * t337;
t345 = pkin(3) * t347 + t239;
t343 = 0.2e1 * t333 + t350;
t342 = pkin(2) * t305 + qJ(3) * t347 + t352;
t341 = pkin(2) * t347 - qJ(3) * t305 + t345;
t324 = 0.2e1 * t362;
t317 = qJDD(1) * t337 + t339 * t340;
t312 = pkin(1) * t316;
t296 = -pkin(4) * t317 + g(3) * t339;
t267 = -t280 + t312;
t260 = t312 + t343;
t259 = -t266 + t324 + t325;
t254 = pkin(1) * t257;
t251 = qJ(2) * t358 + t380;
t238 = -qJ(2) * t246 + (-pkin(2) * t334 + qJ(3) * t335) * t332;
t237 = qJ(2) * t359 + (pkin(2) * t335 + qJ(3) * t334 + pkin(1)) * t332;
t236 = pkin(1) * t246 + t363;
t235 = -pkin(1) * t276 + t341;
t234 = pkin(1) * t274 + t342;
t227 = t374 + t377;
t226 = t381 * t332 - t376;
t221 = -qJ(2) * t398 - t226 * t334 + t227 * t335;
t220 = -qJ(2) * t222 + t226 * t335 + t227 * t334 + t380;
t219 = pkin(1) * t398 + t353;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t318, 0, -t317, 0, t351, -t296, -t348, -pkin(4) * t348, 0, 0, t356, 0, -t283, 0, -t382, -t252, t404, pkin(4) * t404 + qJ(2) * t372 - t337 * t251, 0, t356, 0, 0, t283, 0, -t382, -t396, t252, -pkin(4) * t396 - t337 * t237 + t339 * t238, 0, 0, t395, 0, -t408, 0, t415, t414, t410, -pkin(4) * t410 - t337 * t220 + t339 * t221; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t317, 0, t318, 0, t296, t351, t354, pkin(4) * t354, 0, 0, t283, 0, t356, 0, t252, -t382, t403, pkin(4) * t403 + qJ(2) * t373 + t339 * t251, 0, t283, 0, 0, -t356, 0, t252, t397, t382, pkin(4) * t397 + t339 * t237 + t337 * t238, 0, 0, -t408, 0, -t395, 0, t414, -t415, t409, -pkin(4) * t409 + t339 * t220 + t337 * t221; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t319, t320, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t267, t266, 0, -t254, 0, 0, 0, qJDD(1), 0, 0, t260, 0, t259, t236, 0, 0, 0, 0, 0, t328, t235, t234, 0, t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t340, 0, 0, -g(3), -t319, 0, 0, 0, t316, 0, -t315, 0, -t290, -t360, t257, qJ(2) * t257, 0, t316, 0, 0, t315, 0, -t290, -t246, t360, t238, 0, 0, t276, 0, -t274, 0, t243, t407, t398, t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, 0, qJDD(1), 0, g(3), 0, -t320, 0, 0, 0, t315, 0, t316, 0, t360, -t290, t358, t251, 0, t315, 0, 0, -t316, 0, t360, t359, t290, t237, 0, 0, -t274, 0, -t276, 0, t407, -t243, t222, t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t319, t320, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t267, t266, 0, -t254, 0, 0, 0, qJDD(1), 0, 0, t260, 0, t259, t236, 0, 0, 0, 0, 0, t328, t235, t234, 0, t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t340, 0, 0, -t332, t280, 0, 0, qJDD(1), 0, 0, t340, 0, 0, t272, t332, t374, 0, 0, -t347, 0, -t305, 0, t384, t288, t232, t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, 0, qJDD(1), 0, t332, 0, t281, 0, 0, t340, 0, 0, -qJDD(1), 0, t332, t270, 0, pkin(2) * t332, 0, 0, -t305, 0, t347, 0, t288, -t384, -t349, t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t280, -t281, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t343, 0, t324 + t361, t363, 0, 0, 0, 0, 0, t328, t341, t342, 0, t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t340, 0, 0, t272, t332, 0, 0, 0, -t347, 0, -t305, 0, t384, t288, t232, t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t272, 0, t270, 0, 0, 0, 0, 0, 0, t328, t345, t352, 0, t379; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, 0, 0, qJDD(1), 0, -t332, -t270, 0, 0, 0, 0, t305, 0, -t347, 0, -t288, t384, t349, -pkin(3) * t332 + t376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t328, 0, -t327, 0, 0, t332, t239, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t327, 0, -t328, 0, -t332, 0, t240, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t328, -t239, -t240, 0, 0;];
m_new_reg  = t1;
