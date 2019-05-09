% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RPRR1
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RPRR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:16:49
% EndTime: 2019-05-04 19:16:52
% DurationCPUTime: 2.87s
% Computational Cost: add. (14139->195), mult. (21732->222), div. (0->0), fcn. (13344->8), ass. (0->126)
t432 = qJD(1) + qJD(3);
t428 = qJD(4) + t432;
t425 = t428 ^ 2;
t431 = qJDD(1) + qJDD(3);
t426 = qJDD(4) + t431;
t436 = sin(qJ(4));
t439 = cos(qJ(4));
t402 = t439 * t425 + t436 * t426;
t405 = t436 * t425 - t439 * t426;
t437 = sin(qJ(3));
t440 = cos(qJ(3));
t366 = t440 * t402 - t437 * t405;
t433 = g(3) - qJDD(2);
t389 = pkin(6) * t402 - t439 * t433;
t492 = pkin(6) * t405 - t436 * t433;
t337 = pkin(5) * t366 + t440 * t389 - t437 * t492;
t370 = t437 * t402 + t440 * t405;
t434 = sin(pkin(7));
t435 = cos(pkin(7));
t342 = t435 * t366 - t434 * t370;
t505 = pkin(5) * t370 + t437 * t389 + t440 * t492;
t320 = qJ(2) * t342 + t435 * t337 - t434 * t505;
t438 = sin(qJ(1));
t441 = cos(qJ(1));
t346 = t434 * t366 + t435 * t370;
t507 = t441 * t342 - t438 * t346;
t520 = qJ(2) * t346 + t434 * t337 + t435 * t505;
t532 = pkin(4) * t507 + t441 * t320 - t438 * t520;
t522 = t438 * t342 + t441 * t346;
t531 = pkin(4) * t522 + t438 * t320 + t441 * t520;
t423 = t438 * g(1) - t441 * g(2);
t417 = qJDD(1) * pkin(1) + t423;
t424 = t441 * g(1) + t438 * g(2);
t442 = qJD(1) ^ 2;
t418 = -t442 * pkin(1) - t424;
t382 = -t435 * t417 + t434 * t418;
t380 = qJDD(1) * pkin(2) - t382;
t383 = t434 * t417 + t435 * t418;
t381 = -t442 * pkin(2) + t383;
t355 = -t440 * t380 + t437 * t381;
t349 = t431 * pkin(3) - t355;
t356 = t437 * t380 + t440 * t381;
t430 = t432 ^ 2;
t350 = -t430 * pkin(3) + t356;
t324 = -t439 * t349 + t436 * t350;
t325 = t436 * t349 + t439 * t350;
t457 = t436 * t324 + t439 * t325;
t309 = t439 * t324 - t436 * t325;
t461 = t440 * t309;
t300 = -t437 * t457 + t461;
t467 = t437 * t309;
t495 = t440 * t457 + t467;
t295 = t435 * t300 - t434 * t495;
t516 = t434 * t300 + t435 * t495;
t529 = t438 * t295 + t441 * t516;
t528 = t441 * t295 - t438 * t516;
t410 = t440 * t430 + t437 * t431;
t413 = t437 * t430 - t440 * t431;
t373 = t435 * t410 - t434 * t413;
t393 = pkin(5) * t410 - t440 * t433;
t493 = pkin(5) * t413 - t437 * t433;
t341 = qJ(2) * t373 + t435 * t393 - t434 * t493;
t377 = t434 * t410 + t435 * t413;
t490 = t441 * t373 - t438 * t377;
t504 = qJ(2) * t377 + t434 * t393 + t435 * t493;
t523 = pkin(4) * t490 + t441 * t341 - t438 * t504;
t506 = t438 * t373 + t441 * t377;
t521 = pkin(4) * t506 + t438 * t341 + t441 * t504;
t456 = t437 * t355 + t440 * t356;
t330 = t440 * t355 - t437 * t356;
t469 = t435 * t330;
t313 = -t434 * t456 + t469;
t470 = t434 * t330;
t494 = t435 * t456 + t470;
t515 = t438 * t313 + t441 * t494;
t514 = t441 * t313 - t438 * t494;
t455 = t434 * t382 + t435 * t383;
t360 = t435 * t382 - t434 * t383;
t459 = t441 * t360;
t497 = -t438 * t455 + t459;
t465 = t438 * t360;
t496 = t441 * t455 + t465;
t419 = t434 * qJDD(1) + t435 * t442;
t398 = qJ(2) * t419 - t435 * t433;
t420 = t435 * qJDD(1) - t434 * t442;
t448 = -qJ(2) * t420 - t434 * t433;
t475 = t441 * t419 + t438 * t420;
t489 = pkin(4) * t475 + t441 * t398 - t438 * t448;
t384 = -t438 * t419 + t441 * t420;
t488 = -pkin(4) * t384 + t438 * t398 + t441 * t448;
t473 = pkin(1) * t433;
t472 = pkin(2) * t433;
t306 = pkin(3) * t309;
t458 = -pkin(2) * t300 - t306;
t452 = -t438 * t423 - t441 * t424;
t422 = t441 * qJDD(1) - t438 * t442;
t451 = -pkin(4) * t422 - t438 * g(3);
t450 = -pkin(3) * t405 - t324;
t449 = -pkin(2) * t413 - t355;
t447 = t441 * t423 - t438 * t424;
t446 = -pkin(2) * t370 + t450;
t445 = -pkin(3) * t402 - t325;
t444 = -pkin(2) * t410 - t356;
t443 = -pkin(2) * t366 + t445;
t421 = t438 * qJDD(1) + t441 * t442;
t406 = -pkin(4) * t421 + t441 * g(3);
t363 = pkin(1) * t420 - t382;
t362 = -pkin(1) * t419 - t383;
t357 = pkin(1) * t360;
t354 = qJ(2) * t455 + t473;
t333 = -pkin(1) * t377 + t449;
t332 = -pkin(1) * t373 + t444;
t327 = pkin(2) * t330;
t326 = pkin(5) * t456 + t472;
t316 = -pkin(1) * t346 + t446;
t315 = -pkin(1) * t342 + t443;
t305 = pkin(3) * t433 + pkin(6) * t457;
t304 = -pkin(1) * t313 - t327;
t303 = pkin(5) * t469 + qJ(2) * t313 - t434 * t326;
t302 = pkin(5) * t470 + qJ(2) * t494 + t435 * t326 + t473;
t292 = pkin(5) * t300 + pkin(6) * t461 - t437 * t305;
t291 = pkin(5) * t495 + pkin(6) * t467 + t440 * t305 + t472;
t290 = -pkin(1) * t295 + t458;
t289 = qJ(2) * t295 - t434 * t291 + t435 * t292;
t288 = qJ(2) * t516 + t435 * t291 + t434 * t292 + t473;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t422, 0, -t421, 0, t451, -t406, -t447, -pkin(4) * t447, 0, 0, t384, 0, -t475, 0, t488, t489, t497, pkin(4) * t497 + qJ(2) * t459 - t438 * t354, 0, 0, -t506, 0, -t490, 0, t521, t523, t514, pkin(4) * t514 - t438 * t302 + t441 * t303, 0, 0, -t522, 0, -t507, 0, t531, t532, t528, pkin(4) * t528 - t438 * t288 + t441 * t289; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t421, 0, t422, 0, t406, t451, t452, pkin(4) * t452, 0, 0, t475, 0, t384, 0, -t489, t488, t496, pkin(4) * t496 + qJ(2) * t465 + t441 * t354, 0, 0, t490, 0, -t506, 0, -t523, t521, t515, pkin(4) * t515 + t441 * t302 + t438 * t303, 0, 0, t507, 0, -t522, 0, -t532, t531, t529, pkin(4) * t529 + t441 * t288 + t438 * t289; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t423, t424, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t363, t362, 0, -t357, 0, 0, 0, 0, 0, t431, t333, t332, 0, t304, 0, 0, 0, 0, 0, t426, t316, t315, 0, t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t442, 0, 0, -g(3), -t423, 0, 0, 0, t420, 0, -t419, 0, t448, t398, t360, qJ(2) * t360, 0, 0, -t377, 0, -t373, 0, t504, t341, t313, t303, 0, 0, -t346, 0, -t342, 0, t520, t320, t295, t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t442, 0, qJDD(1), 0, g(3), 0, -t424, 0, 0, 0, t419, 0, t420, 0, -t398, t448, t455, t354, 0, 0, t373, 0, -t377, 0, -t341, t504, t494, t302, 0, 0, t342, 0, -t346, 0, -t320, t520, t516, t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t423, t424, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t363, t362, 0, -t357, 0, 0, 0, 0, 0, t431, t333, t332, 0, t304, 0, 0, 0, 0, 0, t426, t316, t315, 0, t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t442, 0, 0, -t433, t382, 0, 0, 0, -t413, 0, -t410, 0, t493, t393, t330, pkin(5) * t330, 0, 0, -t370, 0, -t366, 0, t505, t337, t300, t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t442, 0, qJDD(1), 0, t433, 0, t383, 0, 0, 0, t410, 0, -t413, 0, -t393, t493, t456, t326, 0, 0, t366, 0, -t370, 0, -t337, t505, t495, t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t382, -t383, 0, 0, 0, 0, 0, 0, 0, t431, t449, t444, 0, -t327, 0, 0, 0, 0, 0, t426, t446, t443, 0, t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t431, 0, -t430, 0, 0, -t433, t355, 0, 0, 0, -t405, 0, -t402, 0, t492, t389, t309, pkin(6) * t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t430, 0, t431, 0, t433, 0, t356, 0, 0, 0, t402, 0, -t405, 0, -t389, t492, t457, t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t431, -t355, -t356, 0, 0, 0, 0, 0, 0, 0, t426, t450, t445, 0, -t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t426, 0, -t425, 0, 0, -t433, t324, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t425, 0, t426, 0, t433, 0, t325, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t426, -t324, -t325, 0, 0;];
m_new_reg  = t1;
