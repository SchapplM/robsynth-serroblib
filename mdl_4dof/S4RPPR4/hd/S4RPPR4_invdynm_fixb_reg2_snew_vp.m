% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RPPR4
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RPPR4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:58
% EndTime: 2019-12-31 16:39:00
% DurationCPUTime: 2.77s
% Computational Cost: add. (4012->263), mult. (7379->311), div. (0->0), fcn. (3984->6), ass. (0->170)
t423 = sin(pkin(6));
t424 = cos(pkin(6));
t432 = qJD(1) ^ 2;
t393 = -t424 * qJDD(1) + t423 * t432;
t421 = g(3) - qJDD(2);
t369 = qJ(2) * t393 - t423 * t421;
t428 = sin(qJ(1));
t430 = cos(qJ(1));
t392 = t423 * qJDD(1) + t424 * t432;
t455 = t430 * t392 - t428 * t393;
t459 = -qJ(2) * t392 + t424 * t421;
t515 = -pkin(4) * t455 + t428 * t369 + t430 * t459;
t401 = t428 * g(1) - t430 * g(2);
t384 = qJDD(1) * pkin(1) + t401;
t402 = t430 * g(1) + t428 * g(2);
t385 = -t432 * pkin(1) - t402;
t338 = -t424 * t384 + t423 * t385;
t339 = t423 * t384 + t424 * t385;
t456 = t423 * t338 + t424 * t339;
t304 = t424 * t338 - t423 * t339;
t472 = t430 * t304;
t512 = -t428 * t456 + t472;
t479 = t428 * t304;
t511 = t430 * t456 + t479;
t453 = t428 * t392 + t430 * t393;
t497 = pkin(4) * t453 + t430 * t369 - t428 * t459;
t434 = (2 * qJD(3) * qJD(1)) + t339;
t464 = qJDD(1) * qJ(3);
t326 = -t432 * pkin(2) + t434 + t464;
t422 = qJDD(1) * pkin(2);
t443 = qJDD(3) + t338;
t328 = -t432 * qJ(3) - t422 + t443;
t297 = t423 * t326 - t424 * t328;
t457 = t424 * t326 + t423 * t328;
t508 = -t428 * t297 + t430 * t457;
t507 = t430 * t297 + t428 * t457;
t325 = -qJDD(1) * pkin(5) + t328;
t427 = sin(qJ(4));
t429 = cos(qJ(4));
t307 = -t429 * t325 - t427 * t421;
t308 = t427 * t325 - t429 * t421;
t288 = -t429 * t307 + t427 * t308;
t496 = -pkin(2) - pkin(5);
t495 = pkin(1) * t392;
t494 = pkin(1) * t393;
t493 = pkin(3) * t288;
t317 = -t432 * pkin(5) + t326;
t492 = pkin(3) * t317;
t419 = t427 ^ 2;
t420 = t429 ^ 2;
t466 = t419 + t420;
t394 = t466 * qJDD(1);
t491 = pkin(3) * t394;
t490 = pkin(5) * t288;
t488 = t419 * t432;
t487 = t420 * t432;
t486 = t423 * t394;
t485 = t424 * t394;
t482 = t427 * t317;
t407 = t427 * t432 * t429;
t399 = qJDD(4) + t407;
t481 = t427 * t399;
t400 = qJDD(4) - t407;
t480 = t427 * t400;
t313 = t429 * t317;
t474 = t429 * t399;
t473 = t429 * t400;
t467 = -pkin(2) * t328 + qJ(3) * t326;
t465 = qJD(1) * qJD(4);
t463 = t427 * qJDD(1);
t413 = t429 * qJDD(1);
t462 = t427 * t465;
t461 = t429 * t465;
t431 = qJD(4) ^ 2;
t406 = -t431 - t487;
t358 = t429 * t406 - t481;
t460 = -pkin(5) * t358 + t313;
t451 = -t428 * t401 - t430 * t402;
t450 = t423 * t407;
t449 = t424 * t407;
t448 = -pkin(2) * t288 + qJ(3) * t317 - t490;
t396 = t430 * qJDD(1) - t428 * t432;
t447 = -pkin(4) * t396 - t428 * g(3);
t388 = 0.2e1 * t461 + t463;
t446 = pkin(3) * t388 + t313;
t391 = t413 - 0.2e1 * t462;
t445 = pkin(3) * t391 - t482;
t404 = -t431 - t488;
t356 = t427 * t404 + t473;
t444 = -pkin(5) * t356 + t482;
t289 = t427 * t307 + t429 * t308;
t442 = t430 * t401 - t428 * t402;
t441 = -pkin(2) * t358 + qJ(3) * t391 + t460;
t440 = -pkin(3) * t356 + t307;
t439 = pkin(5) * t394 - t288;
t438 = -0.2e1 * t422 + t443;
t437 = -pkin(2) * t356 + qJ(3) * t388 + t444;
t436 = -pkin(3) * t358 + t308;
t397 = t466 * t432;
t435 = pkin(2) * t394 - qJ(3) * t397 + t439;
t433 = 0.2e1 * t464 + t434;
t405 = t431 - t487;
t403 = -t431 + t488;
t398 = (-t419 + t420) * t432;
t395 = t428 * qJDD(1) + t430 * t432;
t390 = t413 - t462;
t389 = -t461 - t463;
t382 = t466 * t465;
t374 = -pkin(4) * t395 + t430 * g(3);
t367 = t427 * t390 + t420 * t465;
t366 = t429 * t389 + t419 * t465;
t365 = t424 * qJDD(4) - t423 * t382;
t364 = t423 * qJDD(4) + t424 * t382;
t363 = -t427 * t406 - t474;
t362 = -t427 * t405 + t473;
t361 = (t390 - t462) * t429;
t360 = t429 * t404 - t480;
t359 = t429 * t403 - t481;
t357 = t429 * t405 + t480;
t355 = t427 * t403 + t474;
t354 = (-t389 + t461) * t427;
t351 = -t424 * t397 - t486;
t350 = -t423 * t397 + t485;
t341 = -t429 * t388 - t427 * t391;
t340 = -t427 * t388 + t429 * t391;
t336 = t423 * t366 - t449;
t335 = t423 * t367 + t449;
t334 = -t424 * t366 - t450;
t333 = -t424 * t367 + t450;
t332 = t423 * t357 + t424 * t413;
t331 = t423 * t355 - t424 * t463;
t330 = -t424 * t357 + t423 * t413;
t329 = -t424 * t355 - t423 * t463;
t323 = t423 * t358 + t424 * t391;
t322 = t423 * t356 + t424 * t388;
t321 = -t424 * t358 + t423 * t391;
t320 = -t424 * t356 + t423 * t388;
t319 = -t338 - t494;
t318 = -t339 - t495;
t312 = t423 * t340 + t424 * t398;
t311 = -t424 * t340 + t423 * t398;
t310 = t438 + t494;
t309 = t433 + t495;
t301 = pkin(1) * t304;
t300 = pkin(1) * t421 + qJ(2) * t456;
t295 = -qJ(3) * t363 - t436;
t294 = -qJ(3) * t360 - t440;
t293 = -qJ(2) * t297 + (-pkin(2) * t423 + qJ(3) * t424) * t421;
t292 = qJ(2) * t457 + (pkin(2) * t424 + qJ(3) * t423 + pkin(1)) * t421;
t291 = t496 * t360 + t446;
t290 = t496 * t363 + t445;
t286 = pkin(3) * t397 + t289;
t285 = pkin(1) * t321 + t441;
t284 = pkin(1) * t320 + t437;
t283 = t423 * t288 + t424 * t317;
t282 = -t424 * t288 + t423 * t317;
t281 = -pkin(3) * t485 - qJ(2) * t350 + t423 * t286;
t280 = -pkin(3) * t486 + qJ(2) * t351 - t424 * t286;
t279 = pkin(1) * t297 + t467;
t278 = pkin(1) * t350 + t435;
t277 = -qJ(3) * t289 + t493;
t276 = -qJ(2) * t321 - t423 * t290 + t424 * t295;
t275 = -qJ(2) * t320 - t423 * t291 + t424 * t294;
t274 = t496 * t289 + t492;
t273 = -pkin(1) * t363 + qJ(2) * t323 + t424 * t290 + t423 * t295;
t272 = -pkin(1) * t360 + qJ(2) * t322 + t424 * t291 + t423 * t294;
t271 = pkin(1) * t282 + t448;
t270 = -qJ(2) * t282 - t423 * t274 + t424 * t277;
t269 = -pkin(1) * t289 + qJ(2) * t283 + t424 * t274 + t423 * t277;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t396, 0, -t395, 0, t447, -t374, -t442, -pkin(4) * t442, 0, 0, -t453, 0, -t455, 0, t497, -t515, t512, pkin(4) * t512 + qJ(2) * t472 - t428 * t300, 0, t453, t455, 0, 0, 0, -t507, -t497, t515, -pkin(4) * t507 - t428 * t292 + t430 * t293, -t428 * t333 + t430 * t335, -t428 * t311 + t430 * t312, -t428 * t330 + t430 * t332, -t428 * t334 + t430 * t336, -t428 * t329 + t430 * t331, -t428 * t364 + t430 * t365, t430 * t275 - t428 * t272 - pkin(4) * (t430 * t320 + t428 * t322), t430 * t276 - t428 * t273 - pkin(4) * (t430 * t321 + t428 * t323), t430 * t281 - t428 * t280 - pkin(4) * (t430 * t350 + t428 * t351), t430 * t270 - t428 * t269 - pkin(4) * (t430 * t282 + t428 * t283); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t395, 0, t396, 0, t374, t447, t451, pkin(4) * t451, 0, 0, t455, 0, -t453, 0, t515, t497, t511, pkin(4) * t511 + qJ(2) * t479 + t430 * t300, 0, -t455, t453, 0, 0, 0, t508, -t515, -t497, pkin(4) * t508 + t430 * t292 + t428 * t293, t430 * t333 + t428 * t335, t430 * t311 + t428 * t312, t430 * t330 + t428 * t332, t430 * t334 + t428 * t336, t430 * t329 + t428 * t331, t430 * t364 + t428 * t365, t428 * t275 + t430 * t272 + pkin(4) * (-t428 * t320 + t430 * t322), t428 * t276 + t430 * t273 + pkin(4) * (-t428 * t321 + t430 * t323), t428 * t281 + t430 * t280 + pkin(4) * (-t428 * t350 + t430 * t351), t428 * t270 + t430 * t269 + pkin(4) * (-t428 * t282 + t430 * t283); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t401, t402, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t319, t318, 0, -t301, qJDD(1), 0, 0, 0, 0, 0, 0, t310, t309, t279, t361, t341, t362, t354, t359, 0, t284, t285, t278, t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t432, 0, 0, -g(3), -t401, 0, 0, 0, -t393, 0, -t392, 0, t369, -t459, t304, qJ(2) * t304, 0, t393, t392, 0, 0, 0, -t297, -t369, t459, t293, t335, t312, t332, t336, t331, t365, t275, t276, t281, t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t432, 0, qJDD(1), 0, g(3), 0, -t402, 0, 0, 0, t392, 0, -t393, 0, t459, t369, t456, t300, 0, -t392, t393, 0, 0, 0, t457, -t459, -t369, t292, t333, t311, t330, t334, t329, t364, t272, t273, t280, t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t401, t402, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t319, t318, 0, -t301, qJDD(1), 0, 0, 0, 0, 0, 0, t310, t309, t279, t361, t341, t362, t354, t359, 0, t284, t285, t278, t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t432, 0, 0, -t421, t338, 0, 0, -qJDD(1), t432, 0, 0, 0, t328, 0, t421, qJ(3) * t421, t407, t398, t413, -t407, -t463, qJDD(4), t294, t295, -t491, t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t432, 0, qJDD(1), 0, t421, 0, t339, 0, 0, -t432, -qJDD(1), 0, 0, 0, t326, -t421, 0, pkin(2) * t421, -t367, -t340, -t357, -t366, -t355, t382, t291, t290, -t286, t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t338, -t339, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t438, t433, t467, t361, t341, t362, t354, t359, 0, t437, t441, t435, t448; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t328, t326, 0, t361, t341, t362, t354, t359, 0, t444, t460, t439, -t490; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t432, 0, 0, 0, -t328, 0, -t421, 0, -t407, -t398, -t413, t407, t463, -qJDD(4), t440, t436, t491, -t493; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t432, qJDD(1), 0, 0, 0, -t326, t421, 0, 0, t367, t340, t357, t366, t355, -t382, pkin(5) * t360 - t446, pkin(5) * t363 - t445, t286, pkin(5) * t289 - t492; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t390, -t388, t400, t462, t403, -t462, 0, t317, t307, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t461, t391, t405, t389, t399, -t461, -t317, 0, t308, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t407, t398, t413, -t407, -t463, qJDD(4), -t307, -t308, 0, 0;];
m_new_reg = t1;
