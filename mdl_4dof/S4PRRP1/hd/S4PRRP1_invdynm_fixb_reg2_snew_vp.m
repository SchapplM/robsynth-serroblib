% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRRP1
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRRP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:02:30
% EndTime: 2019-05-04 19:02:32
% DurationCPUTime: 1.70s
% Computational Cost: add. (4526->164), mult. (6887->169), div. (0->0), fcn. (5092->6), ass. (0->101)
t330 = (qJD(2) + qJD(3));
t328 = t330 ^ 2;
t329 = qJDD(2) + qJDD(3);
t334 = sin(qJ(3));
t336 = cos(qJ(3));
t311 = t334 * t328 - t336 * t329;
t331 = g(3) - qJDD(1);
t295 = pkin(5) * t311 - t334 * t331;
t308 = t336 * t328 + t334 * t329;
t298 = pkin(5) * t308 - t336 * t331;
t335 = sin(qJ(2));
t337 = cos(qJ(2));
t350 = t335 * t308 + t337 * t311;
t252 = pkin(4) * t350 + t337 * t295 + t335 * t298;
t277 = t337 * t308 - t335 * t311;
t255 = pkin(4) * t277 - t335 * t295 + t337 * t298;
t332 = sin(pkin(6));
t333 = cos(pkin(6));
t259 = t333 * t277 - t332 * t350;
t412 = qJ(1) * t259 - t332 * t252 + t333 * t255;
t392 = t332 * t277 + t333 * t350;
t411 = qJ(1) * t392 + t333 * t252 + t332 * t255;
t318 = t332 * g(1) - t333 * g(2);
t319 = t333 * g(1) + t332 * g(2);
t293 = t335 * t318 - t337 * t319;
t338 = qJD(2) ^ 2;
t289 = -t338 * pkin(2) + t293;
t348 = t337 * t318 + t335 * t319;
t340 = qJDD(2) * pkin(2) + t348;
t267 = t334 * t289 - t336 * t340;
t268 = t336 * t289 + t334 * t340;
t353 = t334 * t267 + t336 * t268;
t244 = t336 * t267 - t334 * t268;
t361 = t337 * t244;
t229 = -t335 * t353 + t361;
t363 = t335 * t244;
t401 = t337 * t353 + t363;
t410 = t332 * t229 + t333 * t401;
t409 = t333 * t229 - t332 * t401;
t408 = pkin(1) * t277;
t357 = (2 * qJD(4) * t330) + t268;
t371 = t329 * qJ(4);
t258 = -t328 * pkin(3) + t357 + t371;
t324 = t329 * pkin(3);
t344 = qJDD(4) + t267;
t264 = -t328 * qJ(4) - t324 + t344;
t237 = t334 * t258 - t336 * t264;
t354 = t336 * t258 + t334 * t264;
t224 = t337 * t237 + t335 * t354;
t394 = -t335 * t237 + t337 * t354;
t405 = -t332 * t224 + t333 * t394;
t404 = t333 * t224 + t332 * t394;
t352 = t337 * t293 - t335 * t348;
t272 = -t335 * t293 - t337 * t348;
t366 = t333 * t272;
t400 = -t332 * t352 + t366;
t370 = t332 * t272;
t399 = t333 * t352 + t370;
t396 = pkin(2) * t308;
t316 = t335 * qJDD(2) + t337 * t338;
t302 = pkin(4) * t316 - t337 * t331;
t317 = t337 * qJDD(2) - t335 * t338;
t345 = -pkin(4) * t317 - t335 * t331;
t376 = t333 * t316 + t332 * t317;
t390 = qJ(1) * t376 + t333 * t302 - t332 * t345;
t290 = -t332 * t316 + t333 * t317;
t389 = -qJ(1) * t290 + t332 * t302 + t333 * t345;
t373 = pkin(1) * t331;
t367 = t332 * t331;
t365 = t333 * t331;
t359 = -pkin(3) * t264 + qJ(4) * t258;
t358 = pkin(2) * t237 + t359;
t356 = -t268 - t396;
t349 = -t332 * t318 - t333 * t319;
t347 = 0.2e1 * t371 + t357;
t307 = pkin(2) * t311;
t346 = -t267 - t307;
t343 = t347 + t396;
t342 = t333 * t318 - t332 * t319;
t341 = 0.2e1 * t324 - t344;
t339 = -t307 + t341;
t276 = pkin(1) * t350;
t275 = -pkin(1) * t316 - t293;
t274 = pkin(1) * t317 + t348;
t269 = pkin(1) * t272;
t266 = pkin(4) * t352 + t373;
t249 = -t276 + t346;
t248 = t356 - t408;
t247 = -t276 + t339;
t246 = t343 + t408;
t241 = pkin(2) * t244;
t240 = pkin(2) * t331 + pkin(5) * t353;
t234 = -pkin(5) * t237 + (-pkin(3) * t334 + qJ(4) * t336) * t331;
t233 = pkin(5) * t354 + (pkin(3) * t336 + qJ(4) * t334 + pkin(2)) * t331;
t222 = -pkin(1) * t229 - t241;
t221 = pkin(4) * t229 + pkin(5) * t361 - t335 * t240;
t220 = pkin(4) * t401 + pkin(5) * t363 + t337 * t240 + t373;
t219 = pkin(1) * t224 + t358;
t218 = -pkin(4) * t224 - t335 * t233 + t337 * t234;
t217 = pkin(4) * t394 + t337 * t233 + t335 * t234 + t373;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t367, -t365, -t342, -qJ(1) * t342, 0, 0, t290, 0, -t376, 0, t389, t390, t400, pkin(4) * t366 + qJ(1) * t400 - t332 * t266, 0, 0, -t392, 0, -t259, 0, t411, t412, t409, qJ(1) * t409 - t332 * t220 + t333 * t221, 0, -t392, 0, 0, t259, 0, t411, -t404, -t412, -qJ(1) * t404 - t332 * t217 + t333 * t218; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t365, -t367, t349, qJ(1) * t349, 0, 0, t376, 0, t290, 0, -t390, t389, t399, pkin(4) * t370 + qJ(1) * t399 + t333 * t266, 0, 0, t259, 0, -t392, 0, -t412, t411, t410, qJ(1) * t410 + t333 * t220 + t332 * t221, 0, t259, 0, 0, t392, 0, -t412, t405, -t411, qJ(1) * t405 + t333 * t217 + t332 * t218; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t318, t319, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t274, t275, 0, -t269, 0, 0, 0, 0, 0, t329, t249, t248, 0, t222, 0, 0, 0, t329, 0, 0, t247, 0, t246, t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t331, -t318, 0, 0, 0, t317, 0, -t316, 0, t345, t302, t272, pkin(4) * t272, 0, 0, -t350, 0, -t277, 0, t252, t255, t229, t221, 0, -t350, 0, 0, t277, 0, t252, -t224, -t255, t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331, 0, -t319, 0, 0, 0, t316, 0, t317, 0, -t302, t345, t352, t266, 0, 0, t277, 0, -t350, 0, -t255, t252, t401, t220, 0, t277, 0, 0, t350, 0, -t255, t394, -t252, t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, t319, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t274, t275, 0, -t269, 0, 0, 0, 0, 0, t329, t249, t248, 0, t222, 0, 0, 0, t329, 0, 0, t247, 0, t246, t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t338, 0, 0, -t331, -t348, 0, 0, 0, -t311, 0, -t308, 0, t295, t298, t244, pkin(5) * t244, 0, -t311, 0, 0, t308, 0, t295, -t237, -t298, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, 0, qJDD(2), 0, t331, 0, t293, 0, 0, 0, t308, 0, -t311, 0, -t298, t295, t353, t240, 0, t308, 0, 0, t311, 0, -t298, t354, -t295, t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t348, -t293, 0, 0, 0, 0, 0, 0, 0, t329, t346, t356, 0, -t241, 0, 0, 0, t329, 0, 0, t339, 0, t343, t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329, 0, -t328, 0, 0, -t331, t267, 0, 0, t329, 0, 0, t328, 0, 0, t264, t331, qJ(4) * t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, 0, t329, 0, t331, 0, t268, 0, 0, t328, 0, 0, -t329, 0, t331, t258, 0, pkin(3) * t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329, -t267, -t268, 0, 0, 0, 0, 0, t329, 0, 0, t341, 0, t347, t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329, 0, 0, t328, 0, 0, t264, t331, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329, 0, 0, -t264, 0, t258, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t328, 0, 0, t329, 0, -t331, -t258, 0, 0;];
m_new_reg  = t1;
