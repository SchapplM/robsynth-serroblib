% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RPRP1
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RPRP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:14:09
% EndTime: 2019-05-04 19:14:10
% DurationCPUTime: 1.82s
% Computational Cost: add. (5629->170), mult. (8561->177), div. (0->0), fcn. (5108->6), ass. (0->105)
t355 = (qJD(1) + qJD(3));
t353 = t355 ^ 2;
t354 = qJDD(1) + qJDD(3);
t359 = sin(qJ(3));
t361 = cos(qJ(3));
t333 = t359 * t353 - t361 * t354;
t356 = g(3) - qJDD(2);
t315 = pkin(5) * t333 - t359 * t356;
t330 = t361 * t353 + t359 * t354;
t318 = pkin(5) * t330 - t361 * t356;
t357 = sin(pkin(6));
t358 = cos(pkin(6));
t375 = t357 * t330 + t358 * t333;
t275 = qJ(2) * t375 + t358 * t315 + t357 * t318;
t300 = t358 * t330 - t357 * t333;
t278 = qJ(2) * t300 - t357 * t315 + t358 * t318;
t360 = sin(qJ(1));
t362 = cos(qJ(1));
t281 = t362 * t300 - t360 * t375;
t435 = pkin(4) * t281 - t360 * t275 + t362 * t278;
t415 = t360 * t300 + t362 * t375;
t434 = pkin(4) * t415 + t362 * t275 + t360 * t278;
t343 = t360 * g(1) - t362 * g(2);
t337 = qJDD(1) * pkin(1) + t343;
t344 = t362 * g(1) + t360 * g(2);
t363 = qJD(1) ^ 2;
t338 = -t363 * pkin(1) - t344;
t311 = t357 * t337 + t358 * t338;
t309 = -t363 * pkin(2) + t311;
t310 = -t358 * t337 + t357 * t338;
t365 = qJDD(1) * pkin(2) - t310;
t287 = t359 * t309 - t361 * t365;
t288 = t361 * t309 + t359 * t365;
t378 = t359 * t287 + t361 * t288;
t264 = t361 * t287 - t359 * t288;
t391 = t358 * t264;
t249 = -t357 * t378 + t391;
t393 = t357 * t264;
t422 = t358 * t378 + t393;
t433 = t360 * t249 + t362 * t422;
t432 = t362 * t249 - t360 * t422;
t431 = pkin(1) * t300;
t382 = (2 * qJD(4) * t355) + t288;
t394 = t354 * qJ(4);
t272 = -t353 * pkin(3) + t382 + t394;
t349 = t354 * pkin(3);
t370 = qJDD(4) + t287;
t280 = -t353 * qJ(4) - t349 + t370;
t257 = t359 * t272 - t361 * t280;
t379 = t361 * t272 + t359 * t280;
t244 = t358 * t257 + t357 * t379;
t417 = -t357 * t257 + t358 * t379;
t428 = -t360 * t244 + t362 * t417;
t427 = t362 * t244 + t360 * t417;
t377 = t357 * t310 + t358 * t311;
t292 = t358 * t310 - t357 * t311;
t385 = t362 * t292;
t424 = -t360 * t377 + t385;
t388 = t360 * t292;
t423 = t362 * t377 + t388;
t419 = pkin(2) * t330;
t339 = t357 * qJDD(1) + t358 * t363;
t322 = qJ(2) * t339 - t358 * t356;
t340 = t358 * qJDD(1) - t357 * t363;
t368 = -qJ(2) * t340 - t357 * t356;
t399 = t362 * t339 + t360 * t340;
t413 = pkin(4) * t399 + t362 * t322 - t360 * t368;
t312 = -t360 * t339 + t362 * t340;
t412 = -pkin(4) * t312 + t360 * t322 + t362 * t368;
t396 = pkin(1) * t356;
t384 = -pkin(3) * t280 + qJ(4) * t272;
t383 = pkin(2) * t257 + t384;
t381 = -t288 - t419;
t374 = -t360 * t343 - t362 * t344;
t373 = 0.2e1 * t394 + t382;
t342 = t362 * qJDD(1) - t360 * t363;
t372 = -pkin(4) * t342 - t360 * g(3);
t328 = pkin(2) * t333;
t371 = -t287 - t328;
t369 = t373 + t419;
t367 = t362 * t343 - t360 * t344;
t366 = 0.2e1 * t349 - t370;
t364 = -t328 + t366;
t341 = t360 * qJDD(1) + t362 * t363;
t325 = -pkin(4) * t341 + t362 * g(3);
t296 = pkin(1) * t375;
t295 = pkin(1) * t340 - t310;
t294 = -pkin(1) * t339 - t311;
t289 = pkin(1) * t292;
t286 = qJ(2) * t377 + t396;
t269 = -t296 + t371;
t268 = t381 - t431;
t267 = -t296 + t364;
t266 = t369 + t431;
t261 = pkin(2) * t264;
t260 = pkin(2) * t356 + pkin(5) * t378;
t254 = -pkin(5) * t257 + (-pkin(3) * t359 + qJ(4) * t361) * t356;
t253 = pkin(5) * t379 + (pkin(3) * t361 + qJ(4) * t359 + pkin(2)) * t356;
t242 = -pkin(1) * t249 - t261;
t241 = pkin(5) * t391 + qJ(2) * t249 - t357 * t260;
t240 = pkin(5) * t393 + qJ(2) * t422 + t358 * t260 + t396;
t239 = pkin(1) * t244 + t383;
t238 = -qJ(2) * t244 - t357 * t253 + t358 * t254;
t237 = qJ(2) * t417 + t358 * t253 + t357 * t254 + t396;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t342, 0, -t341, 0, t372, -t325, -t367, -pkin(4) * t367, 0, 0, t312, 0, -t399, 0, t412, t413, t424, pkin(4) * t424 + qJ(2) * t385 - t360 * t286, 0, 0, -t415, 0, -t281, 0, t434, t435, t432, pkin(4) * t432 - t360 * t240 + t362 * t241, 0, -t415, 0, 0, t281, 0, t434, -t427, -t435, -pkin(4) * t427 - t360 * t237 + t362 * t238; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t341, 0, t342, 0, t325, t372, t374, pkin(4) * t374, 0, 0, t399, 0, t312, 0, -t413, t412, t423, pkin(4) * t423 + qJ(2) * t388 + t362 * t286, 0, 0, t281, 0, -t415, 0, -t435, t434, t433, pkin(4) * t433 + t362 * t240 + t360 * t241, 0, t281, 0, 0, t415, 0, -t435, t428, -t434, pkin(4) * t428 + t362 * t237 + t360 * t238; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t343, t344, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t295, t294, 0, -t289, 0, 0, 0, 0, 0, t354, t269, t268, 0, t242, 0, 0, 0, t354, 0, 0, t267, 0, t266, t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t363, 0, 0, -g(3), -t343, 0, 0, 0, t340, 0, -t339, 0, t368, t322, t292, qJ(2) * t292, 0, 0, -t375, 0, -t300, 0, t275, t278, t249, t241, 0, -t375, 0, 0, t300, 0, t275, -t244, -t278, t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t363, 0, qJDD(1), 0, g(3), 0, -t344, 0, 0, 0, t339, 0, t340, 0, -t322, t368, t377, t286, 0, 0, t300, 0, -t375, 0, -t278, t275, t422, t240, 0, t300, 0, 0, t375, 0, -t278, t417, -t275, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t343, t344, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t295, t294, 0, -t289, 0, 0, 0, 0, 0, t354, t269, t268, 0, t242, 0, 0, 0, t354, 0, 0, t267, 0, t266, t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t363, 0, 0, -t356, t310, 0, 0, 0, -t333, 0, -t330, 0, t315, t318, t264, pkin(5) * t264, 0, -t333, 0, 0, t330, 0, t315, -t257, -t318, t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t363, 0, qJDD(1), 0, t356, 0, t311, 0, 0, 0, t330, 0, -t333, 0, -t318, t315, t378, t260, 0, t330, 0, 0, t333, 0, -t318, t379, -t315, t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t310, -t311, 0, 0, 0, 0, 0, 0, 0, t354, t371, t381, 0, -t261, 0, 0, 0, t354, 0, 0, t364, 0, t369, t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, 0, -t353, 0, 0, -t356, t287, 0, 0, t354, 0, 0, t353, 0, 0, t280, t356, qJ(4) * t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, t354, 0, t356, 0, t288, 0, 0, t353, 0, 0, -t354, 0, t356, t272, 0, pkin(3) * t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, -t287, -t288, 0, 0, 0, 0, 0, t354, 0, 0, t366, 0, t373, t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, 0, 0, t353, 0, 0, t280, t356, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, 0, 0, -t280, 0, t272, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, 0, 0, t354, 0, -t356, -t272, 0, 0;];
m_new_reg  = t1;
