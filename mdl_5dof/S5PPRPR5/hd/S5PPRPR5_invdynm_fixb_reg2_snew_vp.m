% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PPRPR5
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PPRPR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:35
% EndTime: 2019-12-31 17:33:39
% DurationCPUTime: 3.92s
% Computational Cost: add. (5244->350), mult. (9029->344), div. (0->0), fcn. (5532->6), ass. (0->191)
t425 = sin(qJ(3));
t429 = qJD(3) ^ 2;
t427 = cos(qJ(3));
t476 = t427 * qJDD(3);
t393 = t425 * t429 - t476;
t418 = g(3) - qJDD(1);
t489 = t425 * t418;
t366 = pkin(5) * t393 + t489;
t420 = sin(pkin(7));
t421 = cos(pkin(7));
t477 = t425 * qJDD(3);
t392 = t427 * t429 + t477;
t485 = t427 * t418;
t452 = pkin(5) * t392 + t485;
t467 = -t420 * t392 + t421 * t393;
t531 = -qJ(1) * t467 + t421 * t366 - t420 * t452;
t469 = t421 * t392 + t420 * t393;
t510 = qJ(1) * t469 - t420 * t366 - t421 * t452;
t396 = t420 * g(1) - t421 * g(2);
t388 = -qJDD(2) + t396;
t397 = t421 * g(1) + t420 * g(2);
t448 = t425 * t388 + t427 * t397;
t480 = (qJD(4) * qJD(3));
t431 = -t448 + (2 * t480);
t479 = qJDD(3) * qJ(4);
t323 = -t429 * pkin(3) + t431 + t479;
t419 = qJDD(3) * pkin(3);
t470 = -t427 * t388 + t425 * t397;
t451 = -qJDD(4) + t470;
t324 = -t429 * qJ(4) - t419 - t451;
t296 = t425 * t323 - t427 * t324;
t450 = t427 * t323 + t425 * t324;
t526 = t420 * t296 + t421 * t450;
t525 = -t421 * t296 + t420 * t450;
t307 = t425 * t448 - t427 * t470;
t449 = -t425 * t470 - t427 * t448;
t524 = t420 * t307 - t421 * t449;
t523 = t421 * t307 + t420 * t449;
t322 = -qJDD(3) * pkin(6) + t324;
t424 = sin(qJ(5));
t426 = cos(qJ(5));
t311 = -t426 * t322 + t424 * t418;
t312 = t424 * t322 + t426 * t418;
t288 = -t426 * t311 + t424 * t312;
t514 = -0.2e1 * t479 - (2 * t480);
t513 = 0.2e1 * t419 - qJDD(4);
t509 = pkin(1) + pkin(2);
t303 = qJ(2) * t393 + t509 * t392 - t448;
t304 = -qJ(2) * t392 + t509 * t393 - t470;
t508 = -pkin(3) - pkin(6);
t507 = pkin(1) * t418;
t506 = pkin(4) * t288;
t321 = -t429 * pkin(6) + t323;
t505 = pkin(4) * t321;
t416 = t424 ^ 2;
t417 = t426 ^ 2;
t482 = t416 + t417;
t391 = t482 * qJDD(3);
t504 = pkin(4) * t391;
t503 = pkin(5) * t307;
t502 = qJ(2) * t418;
t501 = t416 * t429;
t500 = t417 * t429;
t495 = t420 * t418;
t406 = t421 * t418;
t315 = t424 * t321;
t403 = t426 * t429 * t424;
t398 = qJDD(5) + t403;
t491 = t424 * t398;
t490 = t425 * t391;
t316 = t426 * t321;
t488 = t426 * t398;
t399 = qJDD(5) - t403;
t487 = t426 * t399;
t486 = t427 * t391;
t428 = qJD(5) ^ 2;
t401 = -t428 - t501;
t352 = t424 * t401 + t487;
t484 = -pkin(6) * t352 + t315;
t461 = -t428 - t500;
t354 = t426 * t461 - t491;
t483 = -pkin(6) * t354 + t316;
t481 = qJD(3) * qJD(5);
t478 = t424 * qJDD(3);
t408 = t426 * qJDD(3);
t475 = t424 * t481;
t474 = t426 * t481;
t473 = pkin(3) * t489 - pkin(5) * t296;
t410 = pkin(2) * t418;
t472 = -pkin(5) * t450 + t410;
t471 = pkin(5) * t449 - t410;
t379 = t420 * t397;
t465 = t421 * t388 - t379;
t464 = t421 * t396 - t379;
t380 = t421 * t397;
t463 = -t420 * t388 - t380;
t462 = -t420 * t396 - t380;
t460 = t428 - t500;
t459 = t425 * t403;
t458 = t427 * t403;
t457 = -pkin(6) * t391 + t288;
t456 = pkin(3) * t324 - qJ(4) * t323;
t455 = -pkin(3) * t427 - qJ(4) * t425;
t385 = 0.2e1 * t474 + t478;
t454 = pkin(4) * t385 + t316;
t387 = t408 - 0.2e1 * t475;
t453 = pkin(4) * t387 - t315;
t289 = t424 * t311 + t426 * t312;
t447 = pkin(2) * t393 - t470;
t286 = pkin(6) * t288;
t446 = pkin(3) * t288 - qJ(4) * t321 + t286;
t445 = -pkin(4) * t352 + t311;
t272 = t508 * t289 + t505;
t277 = -qJ(4) * t289 + t506;
t281 = -t427 * t288 + t425 * t321;
t444 = -pkin(5) * t281 - t425 * t272 + t427 * t277;
t359 = -t424 * t461 - t488;
t284 = t508 * t359 + t453;
t436 = -pkin(4) * t354 + t312;
t294 = -qJ(4) * t359 - t436;
t318 = -t427 * t354 + t425 * t387;
t443 = -pkin(5) * t318 - t425 * t284 + t427 * t294;
t383 = t424 * t399;
t356 = t426 * t401 - t383;
t285 = t508 * t356 + t454;
t293 = -qJ(4) * t356 - t445;
t317 = -t427 * t352 + t425 * t385;
t442 = -pkin(5) * t317 - t425 * t285 + t427 * t293;
t440 = -t474 - t478;
t439 = pkin(3) * t352 - qJ(4) * t385 - t484;
t438 = pkin(3) * t354 - qJ(4) * t387 - t483;
t437 = pkin(2) * t392 - t448;
t282 = t425 * t288 + t427 * t321;
t435 = -pkin(5) * t282 - t427 * t272 - t425 * t277;
t320 = t425 * t354 + t427 * t387;
t434 = -pkin(5) * t320 - t427 * t284 - t425 * t294;
t319 = t425 * t352 + t427 * t385;
t433 = -pkin(5) * t319 - t427 * t285 - t425 * t293;
t394 = t482 * t429;
t430 = pkin(3) * t391 - qJ(4) * t394 - t457;
t400 = -t428 + t501;
t395 = (-t416 + t417) * t429;
t386 = t408 - t475;
t381 = t482 * t481;
t365 = t427 * qJDD(5) - t425 * t381;
t364 = t425 * qJDD(5) + t427 * t381;
t363 = t424 * t386 + t417 * t481;
t362 = t416 * t481 + t426 * t440;
t358 = -t424 * t460 + t487;
t357 = (t386 - t475) * t426;
t355 = t426 * t400 - t491;
t353 = t426 * t460 + t383;
t351 = t424 * t400 + t488;
t350 = (-t440 + t474) * t424;
t349 = pkin(1) * t388 - qJ(2) * t397;
t346 = -t427 * t394 - t490;
t345 = -t425 * t394 + t486;
t334 = -t426 * t385 - t424 * t387;
t333 = -t424 * t385 + t426 * t387;
t332 = t425 * t362 - t458;
t331 = t425 * t363 + t458;
t330 = -t427 * t363 + t459;
t329 = t427 * t362 + t459;
t328 = t425 * t353 + t426 * t476;
t327 = t425 * t351 - t424 * t476;
t326 = -t427 * t353 + t425 * t408;
t325 = t427 * t351 + t424 * t477;
t314 = t425 * t333 + t427 * t395;
t313 = -t427 * t333 + t425 * t395;
t302 = -t304 + t513;
t301 = -t303 + t514;
t300 = t502 + t503;
t299 = -t471 + t507;
t292 = (-qJ(4) * t427 + qJ(2)) * t418 + t473;
t287 = (pkin(1) - t455) * t418 + t472;
t283 = pkin(4) * t394 + t289;
t280 = -pkin(4) * t486 - pkin(5) * t345 + t425 * t283;
t279 = -pkin(4) * t490 + pkin(5) * t346 - t427 * t283;
t278 = qJ(2) * t449 + t307 * t509;
t275 = qJ(2) * t346 - t509 * t345 - t430;
t274 = qJ(2) * t320 - t509 * t318 + t438;
t273 = qJ(2) * t319 - t509 * t317 + t439;
t271 = qJ(2) * t359 + t443;
t270 = qJ(2) * t356 + t442;
t269 = t509 * t359 + t434;
t268 = t509 * t356 + t433;
t267 = qJ(2) * t450 - t509 * t296 + t456;
t266 = qJ(2) * t282 - t509 * t281 + t446;
t265 = qJ(2) * t289 + t444;
t264 = t509 * t289 + t435;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t495, -t406, -t464, -qJ(1) * t464, 0, 0, 0, 0, 0, 0, -t495, -t465, t406, -qJ(1) * t465 + (-t420 * pkin(1) + t421 * qJ(2)) * t418, 0, 0, -t467, 0, -t469, 0, t531, -t510, t523, -qJ(1) * t523 - t420 * t299 + t421 * t300, 0, t467, t469, 0, 0, 0, t525, -t531, t510, -qJ(1) * t525 - t420 * t287 + t421 * t292, t420 * t330 + t421 * t331, t420 * t313 + t421 * t314, t420 * t326 + t421 * t328, -t420 * t329 + t421 * t332, -t420 * t325 + t421 * t327, t420 * t364 + t421 * t365, t421 * t270 - t420 * t268 - qJ(1) * (-t421 * t317 + t420 * t319), t421 * t271 - t420 * t269 - qJ(1) * (-t421 * t318 + t420 * t320), t421 * t280 + t420 * t279 - qJ(1) * (-t421 * t345 + t420 * t346), t421 * t265 - t420 * t264 - qJ(1) * (-t421 * t281 + t420 * t282); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t406, -t495, t462, qJ(1) * t462, 0, 0, 0, 0, 0, 0, t406, t463, t495, qJ(1) * t463 + (t421 * pkin(1) + t420 * qJ(2)) * t418, 0, 0, -t469, 0, t467, 0, -t510, -t531, t524, -qJ(1) * t524 + t421 * t299 + t420 * t300, 0, t469, -t467, 0, 0, 0, -t526, t510, t531, qJ(1) * t526 + t421 * t287 + t420 * t292, -t421 * t330 + t420 * t331, -t421 * t313 + t420 * t314, -t421 * t326 + t420 * t328, t421 * t329 + t420 * t332, t421 * t325 + t420 * t327, -t421 * t364 + t420 * t365, t420 * t270 + t421 * t268 + qJ(1) * (t420 * t317 + t421 * t319), t420 * t271 + t421 * t269 + qJ(1) * (t420 * t318 + t421 * t320), t420 * t280 - t421 * t279 + qJ(1) * (t420 * t345 + t421 * t346), t420 * t265 + t421 * t264 + qJ(1) * (t420 * t281 + t421 * t282); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t396, t397, 0, 0, 0, 0, 0, 0, 0, 0, t388, 0, -t397, t349, 0, 0, 0, 0, 0, -qJDD(3), t304, t303, 0, t278, -qJDD(3), 0, 0, 0, 0, 0, 0, t302, t301, t267, -t357, -t334, -t358, -t350, -t355, 0, t273, t274, t275, t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t418, -t396, 0, 0, 0, 0, 0, 0, 0, 0, -t388, t418, t502, 0, 0, -t393, 0, -t392, 0, t366, t452, t307, t300, 0, t393, t392, 0, 0, 0, -t296, -t366, -t452, t292, t331, t314, t328, t332, t327, t365, t270, t271, t280, t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, 0, -t397, 0, 0, 0, 0, 0, 0, 0, t418, -t397, 0, t507, 0, 0, -t392, 0, t393, 0, t452, -t366, -t449, t299, 0, t392, -t393, 0, 0, 0, -t450, -t452, t366, t287, -t330, -t313, -t326, t329, t325, -t364, t268, t269, -t279, t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, t397, 0, 0, 0, 0, 0, 0, 0, 0, t388, 0, -t397, t349, 0, 0, 0, 0, 0, -qJDD(3), t304, t303, 0, t278, -qJDD(3), 0, 0, 0, 0, 0, 0, t302, t301, t267, -t357, -t334, -t358, -t350, -t355, 0, t273, t274, t275, t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t388, t418, 0, 0, 0, -t393, 0, -t392, 0, t366, t452, t307, t503, 0, t393, t392, 0, 0, 0, -t296, -t366, -t452, -qJ(4) * t485 + t473, t331, t314, t328, t332, t327, t365, t442, t443, t280, t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t388, 0, -t397, 0, 0, 0, 0, 0, 0, -qJDD(3), t447, t437, 0, pkin(2) * t307, -qJDD(3), 0, 0, 0, 0, 0, 0, -t447 + t513, -t437 + t514, -pkin(2) * t296 + t456, -t357, -t334, -t358, -t350, -t355, 0, -pkin(2) * t317 + t439, -pkin(2) * t318 + t438, -pkin(2) * t345 - t430, -pkin(2) * t281 + t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t418, t397, 0, 0, 0, 0, t392, 0, -t393, 0, -t452, t366, t449, t471, 0, -t392, t393, 0, 0, 0, t450, t452, -t366, t418 * t455 - t472, t330, t313, t326, -t329, -t325, t364, -pkin(2) * t356 - t433, -pkin(2) * t359 - t434, t279, -pkin(2) * t289 - t435; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t429, 0, 0, t418, -t470, 0, 0, -qJDD(3), t429, 0, 0, 0, t324, 0, -t418, -qJ(4) * t418, t403, t395, t408, -t403, -t478, qJDD(5), t293, t294, -t504, t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t429, 0, qJDD(3), 0, -t418, 0, -t448, 0, 0, -t429, -qJDD(3), 0, 0, 0, t323, t418, 0, -pkin(3) * t418, -t363, -t333, -t353, -t362, -t351, t381, t285, t284, -t283, t272; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t470, t448, 0, 0, qJDD(3), 0, 0, 0, 0, 0, 0, -0.2e1 * t419 - t451, t431 + 0.2e1 * t479, -t456, t357, t334, t358, t350, t355, 0, -t439, -t438, t430, -t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, 0, 0, 0, 0, 0, t324, t323, 0, t357, t334, t358, t350, t355, 0, t484, t483, -t457, -t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t429, 0, 0, 0, -t324, 0, t418, 0, -t403, -t395, -t408, t403, t478, -qJDD(5), t445, t436, t504, -t506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t429, qJDD(3), 0, 0, 0, -t323, -t418, 0, 0, t363, t333, t353, t362, t351, -t381, pkin(6) * t356 - t454, pkin(6) * t359 - t453, t283, pkin(6) * t289 - t505; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t386, -t385, t399, t475, t400, -t475, 0, t321, t311, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t474, t387, t460, t440, t398, -t474, -t321, 0, t312, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t403, t395, t408, -t403, -t478, qJDD(5), -t311, -t312, 0, 0;];
m_new_reg = t1;
