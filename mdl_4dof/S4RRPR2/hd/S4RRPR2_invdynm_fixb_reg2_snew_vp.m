% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RRPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRPR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:38
% EndTime: 2019-07-18 18:16:40
% DurationCPUTime: 1.63s
% Computational Cost: add. (4832->172), mult. (5823->166), div. (0->0), fcn. (3048->6), ass. (0->107)
t349 = (qJD(1) + qJD(2));
t342 = qJD(4) - t349;
t336 = t342 ^ 2;
t348 = qJDD(1) + qJDD(2);
t337 = -qJDD(4) + t348;
t351 = sin(qJ(4));
t354 = cos(qJ(4));
t311 = t336 * t354 - t337 * t351;
t352 = sin(qJ(2));
t355 = cos(qJ(2));
t367 = t336 * t351 + t354 * t337;
t288 = t311 * t352 - t355 * t367;
t270 = pkin(5) * t288 + (t351 * t355 - t352 * t354) * g(3);
t353 = sin(qJ(1));
t356 = cos(qJ(1));
t286 = t355 * t311 + t352 * t367;
t408 = t286 * t353 + t288 * t356;
t343 = t351 * g(3);
t345 = t355 * g(3);
t422 = -pkin(5) * t286 + t352 * t343 + t354 * t345;
t427 = pkin(4) * t408 + t270 * t356 - t353 * t422;
t419 = t286 * t356 - t288 * t353;
t426 = -pkin(4) * t419 + t270 * t353 + t356 * t422;
t347 = t349 ^ 2;
t328 = g(1) * t356 + g(2) * t353;
t359 = qJD(1) ^ 2;
t324 = -pkin(1) * t359 - t328;
t327 = g(1) * t353 - t356 * g(2);
t364 = qJDD(1) * pkin(1) + t327;
t297 = t355 * t324 + t352 * t364;
t383 = (2 * qJD(3) * t349) + t297;
t387 = t348 * qJ(3);
t365 = t383 + t387;
t274 = (-pkin(2) - pkin(3)) * t347 + t365;
t340 = t348 * pkin(2);
t296 = t324 * t352 - t355 * t364;
t371 = -qJDD(3) - t296;
t284 = -t347 * qJ(3) - t340 - t371;
t279 = -pkin(3) * t348 + t284;
t251 = t274 * t351 - t354 * t279;
t252 = t354 * t274 + t351 * t279;
t244 = t251 * t354 - t252 * t351;
t370 = t251 * t351 + t252 * t354;
t236 = t244 * t352 - t355 * t370;
t411 = t244 * t355 + t352 * t370;
t421 = -t236 * t353 + t356 * t411;
t420 = t236 * t356 + t353 * t411;
t319 = t347 * t355 + t348 * t352;
t322 = t347 * t352 - t348 * t355;
t291 = t319 * t356 - t322 * t353;
t392 = g(3) * t352;
t299 = pkin(5) * t322 - t392;
t302 = pkin(5) * t319 - t345;
t418 = pkin(4) * t291 - t353 * t299 + t356 * t302;
t377 = t319 * t353 + t356 * t322;
t417 = pkin(4) * t377 + t356 * t299 + t353 * t302;
t380 = t296 * t352 + t355 * t297;
t265 = t296 * t355 - t297 * t352;
t390 = t265 * t356;
t414 = -t353 * t380 + t390;
t391 = t265 * t353;
t413 = t356 * t380 + t391;
t412 = pkin(1) * t319;
t280 = -pkin(2) * t347 + t365;
t254 = t280 * t352 - t284 * t355;
t381 = t355 * t280 + t284 * t352;
t410 = -t254 * t353 + t356 * t381;
t409 = t254 * t356 + t353 * t381;
t394 = pkin(3) * g(3);
t393 = pkin(3) * t244;
t350 = qJ(3) * g(3);
t386 = -pkin(2) * t284 + qJ(3) * t280;
t358 = pkin(1) * g(3);
t384 = t352 * t350 + t358;
t376 = -t327 * t353 - t356 * t328;
t375 = 0.2e1 * t387 + t383;
t374 = pkin(2) * t244 + qJ(3) * t370 + t393;
t373 = pkin(3) * t311 + t252;
t326 = qJDD(1) * t356 - t353 * t359;
t372 = -pkin(4) * t326 - g(3) * t353;
t368 = t327 * t356 - t328 * t353;
t363 = pkin(3) * t367 + t251;
t362 = 0.2e1 * t340 + t371;
t361 = pkin(2) * t311 + qJ(3) * t367 + t373;
t360 = pkin(2) * t367 - qJ(3) * t311 + t363;
t357 = pkin(2) * g(3);
t344 = t354 * g(3);
t339 = qJ(3) * t345;
t335 = t357 + t394;
t325 = qJDD(1) * t353 + t356 * t359;
t315 = pkin(1) * t322;
t314 = -pkin(4) * t325 + g(3) * t356;
t282 = -t296 - t315;
t281 = -t297 - t412;
t273 = -t315 + t362;
t267 = t375 + t412;
t262 = pkin(1) * t265;
t261 = pkin(5) * t380 + t358;
t250 = -pkin(2) * t392 - pkin(5) * t254 + t339;
t249 = pkin(2) * t345 + pkin(5) * t381 + t384;
t248 = pkin(1) * t254 + t386;
t247 = -pkin(1) * t288 + t360;
t246 = pkin(1) * t286 + t361;
t235 = -pkin(5) * t411 - t335 * t352 + t339;
t234 = -pkin(5) * t236 + t335 * t355 + t384;
t233 = pkin(1) * t411 + t374;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t326, 0, -t325, 0, t372, -t314, -t368, -pkin(4) * t368, 0, 0, -t377, 0, -t291, 0, t417, t418, t414, pkin(4) * t414 + pkin(5) * t390 - t353 * t261, 0, -t377, 0, 0, t291, 0, t417, -t409, -t418, -pkin(4) * t409 - t353 * t249 + t356 * t250, 0, 0, t408, 0, -t419, 0, t427, t426, t421, -pkin(4) * t421 - t353 * t234 + t356 * t235; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t325, 0, t326, 0, t314, t372, t376, pkin(4) * t376, 0, 0, t291, 0, -t377, 0, -t418, t417, t413, pkin(4) * t413 + pkin(5) * t391 + t356 * t261, 0, t291, 0, 0, t377, 0, -t418, t410, -t417, pkin(4) * t410 + t356 * t249 + t353 * t250, 0, 0, -t419, 0, -t408, 0, t426, -t427, t420, -pkin(4) * t420 + t356 * t234 + t353 * t235; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t327, t328, 0, 0, 0, 0, 0, 0, 0, t348, t282, t281, 0, -t262, 0, 0, 0, t348, 0, 0, t273, 0, t267, t248, 0, 0, 0, 0, 0, t337, t247, t246, 0, t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t359, 0, 0, -g(3), -t327, 0, 0, 0, -t322, 0, -t319, 0, t299, t302, t265, pkin(5) * t265, 0, -t322, 0, 0, t319, 0, t299, -t254, -t302, t250, 0, 0, t288, 0, -t286, 0, t270, t422, t411, t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t359, 0, qJDD(1), 0, g(3), 0, -t328, 0, 0, 0, t319, 0, -t322, 0, -t302, t299, t380, t261, 0, t319, 0, 0, t322, 0, -t302, t381, -t299, t249, 0, 0, -t286, 0, -t288, 0, t422, -t270, t236, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t327, t328, 0, 0, 0, 0, 0, 0, 0, t348, t282, t281, 0, -t262, 0, 0, 0, t348, 0, 0, t273, 0, t267, t248, 0, 0, 0, 0, 0, t337, t247, t246, 0, t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, 0, -t347, 0, 0, -g(3), t296, 0, 0, t348, 0, 0, t347, 0, 0, t284, g(3), t350, 0, 0, -t367, 0, -t311, 0, t343, t344, t244, t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t347, 0, t348, 0, g(3), 0, t297, 0, 0, t347, 0, 0, -t348, 0, g(3), t280, 0, t357, 0, 0, -t311, 0, t367, 0, t344, -t343, -t370, t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, -t296, -t297, 0, 0, 0, 0, 0, t348, 0, 0, t362, 0, t375, t386, 0, 0, 0, 0, 0, t337, t360, t361, 0, t374; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, 0, 0, t347, 0, 0, t284, g(3), 0, 0, 0, -t367, 0, -t311, 0, t343, t344, t244, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, 0, 0, -t284, 0, t280, 0, 0, 0, 0, 0, 0, t337, t363, t373, 0, t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t347, 0, 0, t348, 0, -g(3), -t280, 0, 0, 0, 0, t311, 0, -t367, 0, -t344, t343, t370, -t394; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t337, 0, -t336, 0, 0, g(3), t251, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t336, 0, -t337, 0, -g(3), 0, t252, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t337, -t251, -t252, 0, 0;];
m_new_reg  = t1;
