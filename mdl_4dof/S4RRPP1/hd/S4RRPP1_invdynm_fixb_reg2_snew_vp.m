% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRPP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:19:16
% EndTime: 2019-05-04 19:19:18
% DurationCPUTime: 1.88s
% Computational Cost: add. (6196->174), mult. (8561->178), div. (0->0), fcn. (5108->6), ass. (0->105)
t355 = (qJD(1) + qJD(2));
t353 = t355 ^ 2;
t354 = qJDD(1) + qJDD(2);
t357 = sin(pkin(6));
t358 = cos(pkin(6));
t334 = t353 * t357 - t354 * t358;
t356 = g(3) - qJDD(3);
t315 = qJ(3) * t334 - t356 * t357;
t331 = t353 * t358 + t354 * t357;
t318 = qJ(3) * t331 - t358 * t356;
t359 = sin(qJ(2));
t361 = cos(qJ(2));
t375 = t331 * t359 + t361 * t334;
t275 = pkin(5) * t375 + t315 * t361 + t318 * t359;
t300 = t331 * t361 - t334 * t359;
t278 = pkin(5) * t300 - t315 * t359 + t318 * t361;
t360 = sin(qJ(1));
t362 = cos(qJ(1));
t282 = t300 * t362 - t360 * t375;
t435 = pkin(4) * t282 - t360 * t275 + t362 * t278;
t409 = t300 * t360 + t362 * t375;
t434 = pkin(4) * t409 + t362 * t275 + t360 * t278;
t344 = g(1) * t360 - t362 * g(2);
t340 = qJDD(1) * pkin(1) + t344;
t345 = g(1) * t362 + g(2) * t360;
t363 = qJD(1) ^ 2;
t341 = -pkin(1) * t363 - t345;
t313 = t340 * t359 + t341 * t361;
t308 = -pkin(2) * t353 + t313;
t312 = -t361 * t340 + t341 * t359;
t366 = pkin(2) * t354 - t312;
t286 = t308 * t357 - t358 * t366;
t287 = t358 * t308 + t357 * t366;
t378 = t286 * t357 + t358 * t287;
t264 = t286 * t358 - t287 * t357;
t391 = t264 * t361;
t249 = -t359 * t378 + t391;
t392 = t264 * t359;
t419 = t361 * t378 + t392;
t433 = t249 * t360 + t362 * t419;
t432 = t249 * t362 - t360 * t419;
t431 = pkin(1) * t300;
t382 = (2 * qJD(4) * t355) + t287;
t386 = t354 * qJ(4);
t272 = -pkin(3) * t353 + t382 + t386;
t350 = t354 * pkin(3);
t370 = qJDD(4) + t286;
t280 = -t353 * qJ(4) - t350 + t370;
t257 = t272 * t357 - t280 * t358;
t379 = t358 * t272 + t280 * t357;
t244 = t257 * t361 + t359 * t379;
t411 = -t257 * t359 + t361 * t379;
t428 = -t244 * t360 + t362 * t411;
t427 = t244 * t362 + t360 * t411;
t335 = t353 * t361 + t354 * t359;
t322 = pkin(5) * t335 - g(3) * t361;
t338 = t353 * t359 - t354 * t361;
t368 = t335 * t362 - t338 * t360;
t412 = pkin(5) * t338 - g(3) * t359;
t426 = pkin(4) * t368 + t362 * t322 - t360 * t412;
t407 = t335 * t360 + t338 * t362;
t425 = pkin(4) * t407 + t360 * t322 + t362 * t412;
t377 = t312 * t359 + t361 * t313;
t292 = t312 * t361 - t313 * t359;
t389 = t292 * t362;
t421 = -t360 * t377 + t389;
t390 = t292 * t360;
t420 = t362 * t377 + t390;
t415 = pkin(2) * t331;
t395 = pkin(1) * t356;
t384 = -pkin(3) * t280 + qJ(4) * t272;
t383 = pkin(2) * t257 + t384;
t381 = -t287 - t415;
t374 = -t344 * t360 - t362 * t345;
t373 = 0.2e1 * t386 + t382;
t343 = qJDD(1) * t362 - t360 * t363;
t372 = -pkin(4) * t343 - g(3) * t360;
t328 = pkin(2) * t334;
t371 = -t286 - t328;
t369 = t373 + t415;
t367 = t344 * t362 - t345 * t360;
t365 = 0.2e1 * t350 - t370;
t364 = -t328 + t365;
t342 = qJDD(1) * t360 + t362 * t363;
t329 = -pkin(4) * t342 + g(3) * t362;
t299 = pkin(1) * t375;
t295 = -pkin(1) * t338 - t312;
t294 = -pkin(1) * t335 - t313;
t289 = pkin(1) * t292;
t288 = pkin(1) * g(3) + pkin(5) * t377;
t269 = -t299 + t371;
t268 = t381 - t431;
t267 = -t299 + t364;
t266 = t369 + t431;
t261 = pkin(2) * t264;
t260 = pkin(2) * t356 + qJ(3) * t378;
t254 = -qJ(3) * t257 + (-pkin(3) * t357 + qJ(4) * t358) * t356;
t253 = qJ(3) * t379 + (pkin(3) * t358 + qJ(4) * t357 + pkin(2)) * t356;
t242 = -pkin(1) * t249 - t261;
t241 = pkin(5) * t249 + qJ(3) * t391 - t260 * t359;
t240 = pkin(5) * t419 + qJ(3) * t392 + t260 * t361 + t395;
t239 = pkin(1) * t244 + t383;
t238 = -pkin(5) * t244 - t253 * t359 + t254 * t361;
t237 = pkin(5) * t411 + t253 * t361 + t254 * t359 + t395;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t343, 0, -t342, 0, t372, -t329, -t367, -pkin(4) * t367, 0, 0, -t407, 0, -t368, 0, t425, t426, t421, pkin(4) * t421 + pkin(5) * t389 - t360 * t288, 0, 0, -t409, 0, -t282, 0, t434, t435, t432, pkin(4) * t432 - t360 * t240 + t362 * t241, 0, -t409, 0, 0, t282, 0, t434, -t427, -t435, -pkin(4) * t427 - t360 * t237 + t362 * t238; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t342, 0, t343, 0, t329, t372, t374, pkin(4) * t374, 0, 0, t368, 0, -t407, 0, -t426, t425, t420, pkin(4) * t420 + pkin(5) * t390 + t362 * t288, 0, 0, t282, 0, -t409, 0, -t435, t434, t433, pkin(4) * t433 + t362 * t240 + t360 * t241, 0, t282, 0, 0, t409, 0, -t435, t428, -t434, pkin(4) * t428 + t362 * t237 + t360 * t238; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t344, t345, 0, 0, 0, 0, 0, 0, 0, t354, t295, t294, 0, -t289, 0, 0, 0, 0, 0, t354, t269, t268, 0, t242, 0, 0, 0, t354, 0, 0, t267, 0, t266, t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t363, 0, 0, -g(3), -t344, 0, 0, 0, -t338, 0, -t335, 0, t412, t322, t292, pkin(5) * t292, 0, 0, -t375, 0, -t300, 0, t275, t278, t249, t241, 0, -t375, 0, 0, t300, 0, t275, -t244, -t278, t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t363, 0, qJDD(1), 0, g(3), 0, -t345, 0, 0, 0, t335, 0, -t338, 0, -t322, t412, t377, t288, 0, 0, t300, 0, -t375, 0, -t278, t275, t419, t240, 0, t300, 0, 0, t375, 0, -t278, t411, -t275, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t344, t345, 0, 0, 0, 0, 0, 0, 0, t354, t295, t294, 0, -t289, 0, 0, 0, 0, 0, t354, t269, t268, 0, t242, 0, 0, 0, t354, 0, 0, t267, 0, t266, t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, 0, -t353, 0, 0, -g(3), t312, 0, 0, 0, -t334, 0, -t331, 0, t315, t318, t264, qJ(3) * t264, 0, -t334, 0, 0, t331, 0, t315, -t257, -t318, t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, t354, 0, g(3), 0, t313, 0, 0, 0, t331, 0, -t334, 0, -t318, t315, t378, t260, 0, t331, 0, 0, t334, 0, -t318, t379, -t315, t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, -t312, -t313, 0, 0, 0, 0, 0, 0, 0, t354, t371, t381, 0, -t261, 0, 0, 0, t354, 0, 0, t364, 0, t369, t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, 0, -t353, 0, 0, -t356, t286, 0, 0, t354, 0, 0, t353, 0, 0, t280, t356, qJ(4) * t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, t354, 0, t356, 0, t287, 0, 0, t353, 0, 0, -t354, 0, t356, t272, 0, pkin(3) * t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, -t286, -t287, 0, 0, 0, 0, 0, t354, 0, 0, t365, 0, t373, t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, 0, 0, t353, 0, 0, t280, t356, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, 0, 0, -t280, 0, t272, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, 0, 0, t354, 0, -t356, -t272, 0, 0;];
m_new_reg  = t1;
