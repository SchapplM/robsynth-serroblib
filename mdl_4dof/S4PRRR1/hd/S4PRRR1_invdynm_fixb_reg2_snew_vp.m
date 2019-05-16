% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRRR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:04:41
% EndTime: 2019-05-04 19:04:44
% DurationCPUTime: 2.68s
% Computational Cost: add. (11132->189), mult. (17202->214), div. (0->0), fcn. (13328->8), ass. (0->122)
t405 = qJD(2) + qJD(3);
t400 = qJD(4) + t405;
t398 = t400 ^ 2;
t404 = qJDD(2) + qJDD(3);
t399 = qJDD(4) + t404;
t409 = sin(qJ(4));
t412 = cos(qJ(4));
t380 = t398 * t412 + t399 * t409;
t383 = t398 * t409 - t399 * t412;
t410 = sin(qJ(3));
t413 = cos(qJ(3));
t342 = t380 * t413 - t383 * t410;
t406 = g(3) - qJDD(1);
t367 = pkin(6) * t380 - t406 * t412;
t467 = pkin(6) * t383 - t406 * t409;
t315 = pkin(5) * t342 + t367 * t413 - t410 * t467;
t346 = t380 * t410 + t383 * t413;
t411 = sin(qJ(2));
t414 = cos(qJ(2));
t316 = t342 * t414 - t346 * t411;
t479 = pkin(5) * t346 + t367 * t410 + t413 * t467;
t292 = pkin(4) * t316 + t315 * t414 - t411 * t479;
t407 = sin(pkin(7));
t408 = cos(pkin(7));
t320 = t342 * t411 + t346 * t414;
t482 = t316 * t408 - t320 * t407;
t496 = pkin(4) * t320 + t315 * t411 + t414 * t479;
t507 = qJ(1) * t482 + t408 * t292 - t407 * t496;
t497 = t316 * t407 + t320 * t408;
t506 = qJ(1) * t497 + t407 * t292 + t408 * t496;
t396 = g(1) * t407 - t408 * g(2);
t397 = g(1) * t408 + g(2) * t407;
t424 = t414 * t396 + t397 * t411;
t358 = qJDD(2) * pkin(2) + t424;
t363 = t396 * t411 - t397 * t414;
t415 = qJD(2) ^ 2;
t359 = -pkin(2) * t415 + t363;
t333 = -t413 * t358 + t359 * t410;
t329 = pkin(3) * t404 - t333;
t334 = t358 * t410 + t359 * t413;
t403 = t405 ^ 2;
t330 = -pkin(3) * t403 + t334;
t302 = -t412 * t329 + t330 * t409;
t303 = t329 * t409 + t330 * t412;
t430 = t302 * t409 + t412 * t303;
t287 = t302 * t412 - t303 * t409;
t444 = t287 * t413;
t278 = -t410 * t430 + t444;
t445 = t287 * t410;
t469 = t413 * t430 + t445;
t273 = t278 * t414 - t411 * t469;
t491 = t278 * t411 + t414 * t469;
t504 = t273 * t407 + t408 * t491;
t503 = t273 * t408 - t407 * t491;
t386 = t403 * t413 + t404 * t410;
t389 = t403 * t410 - t404 * t413;
t350 = t386 * t414 - t389 * t411;
t371 = pkin(5) * t386 - t406 * t413;
t468 = pkin(5) * t389 - t406 * t410;
t325 = pkin(4) * t350 + t371 * t414 - t411 * t468;
t354 = t386 * t411 + t389 * t414;
t465 = t350 * t408 - t354 * t407;
t480 = pkin(4) * t354 + t371 * t411 + t414 * t468;
t498 = qJ(1) * t465 + t408 * t325 - t407 * t480;
t481 = t350 * t407 + t354 * t408;
t495 = qJ(1) * t481 + t407 * t325 + t408 * t480;
t429 = t333 * t410 + t413 * t334;
t308 = t333 * t413 - t334 * t410;
t442 = t308 * t414;
t297 = -t411 * t429 + t442;
t443 = t308 * t411;
t470 = t414 * t429 + t443;
t490 = t297 * t407 + t408 * t470;
t489 = t297 * t408 - t407 * t470;
t428 = t414 * t363 - t411 * t424;
t338 = -t363 * t411 - t414 * t424;
t440 = t338 * t408;
t472 = -t407 * t428 + t440;
t441 = t338 * t407;
t471 = t408 * t428 + t441;
t394 = qJDD(2) * t411 + t414 * t415;
t378 = pkin(4) * t394 - t406 * t414;
t395 = qJDD(2) * t414 - t411 * t415;
t421 = -pkin(4) * t395 - t406 * t411;
t450 = t408 * t394 + t395 * t407;
t464 = qJ(1) * t450 + t408 * t378 - t407 * t421;
t360 = -t394 * t407 + t408 * t395;
t463 = -qJ(1) * t360 + t407 * t378 + t408 * t421;
t448 = pkin(1) * t406;
t447 = pkin(2) * t406;
t433 = t407 * t406;
t432 = t408 * t406;
t284 = pkin(3) * t287;
t431 = -pkin(2) * t278 - t284;
t425 = -t396 * t407 - t408 * t397;
t423 = -pkin(3) * t383 - t302;
t422 = -pkin(2) * t389 - t333;
t420 = t396 * t408 - t397 * t407;
t419 = -pkin(2) * t346 + t423;
t418 = -pkin(3) * t380 - t303;
t417 = -pkin(2) * t386 - t334;
t416 = -pkin(2) * t342 + t418;
t349 = -pkin(1) * t394 - t363;
t348 = pkin(1) * t395 + t424;
t335 = pkin(1) * t338;
t332 = pkin(4) * t428 + t448;
t311 = -pkin(1) * t354 + t422;
t310 = -pkin(1) * t350 + t417;
t305 = pkin(2) * t308;
t304 = pkin(5) * t429 + t447;
t294 = -pkin(1) * t320 + t419;
t293 = -pkin(1) * t316 + t416;
t283 = pkin(3) * t406 + pkin(6) * t430;
t282 = -pkin(1) * t297 - t305;
t281 = pkin(4) * t297 + pkin(5) * t442 - t304 * t411;
t280 = pkin(4) * t470 + pkin(5) * t443 + t304 * t414 + t448;
t270 = pkin(5) * t278 + pkin(6) * t444 - t283 * t410;
t269 = pkin(5) * t469 + pkin(6) * t445 + t283 * t413 + t447;
t268 = -pkin(1) * t273 + t431;
t267 = pkin(4) * t273 - t269 * t411 + t270 * t414;
t266 = pkin(4) * t491 + t269 * t414 + t270 * t411 + t448;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t433, -t432, -t420, -qJ(1) * t420, 0, 0, t360, 0, -t450, 0, t463, t464, t472, pkin(4) * t440 + qJ(1) * t472 - t407 * t332, 0, 0, -t481, 0, -t465, 0, t495, t498, t489, qJ(1) * t489 - t407 * t280 + t408 * t281, 0, 0, -t497, 0, -t482, 0, t506, t507, t503, qJ(1) * t503 - t407 * t266 + t408 * t267; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t432, -t433, t425, qJ(1) * t425, 0, 0, t450, 0, t360, 0, -t464, t463, t471, pkin(4) * t441 + qJ(1) * t471 + t408 * t332, 0, 0, t465, 0, -t481, 0, -t498, t495, t490, qJ(1) * t490 + t408 * t280 + t407 * t281, 0, 0, t482, 0, -t497, 0, -t507, t506, t504, qJ(1) * t504 + t408 * t266 + t407 * t267; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t396, t397, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t348, t349, 0, -t335, 0, 0, 0, 0, 0, t404, t311, t310, 0, t282, 0, 0, 0, 0, 0, t399, t294, t293, 0, t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t406, -t396, 0, 0, 0, t395, 0, -t394, 0, t421, t378, t338, pkin(4) * t338, 0, 0, -t354, 0, -t350, 0, t480, t325, t297, t281, 0, 0, -t320, 0, -t316, 0, t496, t292, t273, t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t406, 0, -t397, 0, 0, 0, t394, 0, t395, 0, -t378, t421, t428, t332, 0, 0, t350, 0, -t354, 0, -t325, t480, t470, t280, 0, 0, t316, 0, -t320, 0, -t292, t496, t491, t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, t397, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t348, t349, 0, -t335, 0, 0, 0, 0, 0, t404, t311, t310, 0, t282, 0, 0, 0, 0, 0, t399, t294, t293, 0, t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t415, 0, 0, -t406, -t424, 0, 0, 0, -t389, 0, -t386, 0, t468, t371, t308, pkin(5) * t308, 0, 0, -t346, 0, -t342, 0, t479, t315, t278, t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t415, 0, qJDD(2), 0, t406, 0, t363, 0, 0, 0, t386, 0, -t389, 0, -t371, t468, t429, t304, 0, 0, t342, 0, -t346, 0, -t315, t479, t469, t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t424, -t363, 0, 0, 0, 0, 0, 0, 0, t404, t422, t417, 0, -t305, 0, 0, 0, 0, 0, t399, t419, t416, 0, t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t404, 0, -t403, 0, 0, -t406, t333, 0, 0, 0, -t383, 0, -t380, 0, t467, t367, t287, pkin(6) * t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t403, 0, t404, 0, t406, 0, t334, 0, 0, 0, t380, 0, -t383, 0, -t367, t467, t430, t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t404, -t333, -t334, 0, 0, 0, 0, 0, 0, 0, t399, t423, t418, 0, -t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t399, 0, -t398, 0, 0, -t406, t302, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t398, 0, t399, 0, t406, 0, t303, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t399, -t302, -t303, 0, 0;];
m_new_reg  = t1;
