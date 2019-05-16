% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RRPR1
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRPR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:21:45
% EndTime: 2019-05-04 19:21:48
% DurationCPUTime: 3.01s
% Computational Cost: add. (15658->199), mult. (21732->223), div. (0->0), fcn. (13344->8), ass. (0->126)
t439 = qJD(1) + qJD(2);
t435 = qJD(4) + t439;
t433 = t435 ^ 2;
t438 = qJDD(1) + qJDD(2);
t434 = qJDD(4) + t438;
t443 = sin(qJ(4));
t446 = cos(qJ(4));
t410 = t433 * t446 + t434 * t443;
t413 = t433 * t443 - t434 * t446;
t441 = sin(pkin(7));
t442 = cos(pkin(7));
t373 = t410 * t442 - t413 * t441;
t440 = g(3) - qJDD(3);
t396 = pkin(6) * t410 - t440 * t446;
t494 = pkin(6) * t413 - t440 * t443;
t344 = qJ(3) * t373 + t396 * t442 - t441 * t494;
t377 = t410 * t441 + t413 * t442;
t444 = sin(qJ(2));
t447 = cos(qJ(2));
t349 = t373 * t447 - t377 * t444;
t509 = qJ(3) * t377 + t396 * t441 + t442 * t494;
t327 = pkin(5) * t349 + t344 * t447 - t444 * t509;
t445 = sin(qJ(1));
t448 = cos(qJ(1));
t353 = t373 * t444 + t377 * t447;
t514 = t349 * t448 - t353 * t445;
t527 = pkin(5) * t353 + t344 * t444 + t447 * t509;
t539 = pkin(4) * t514 + t448 * t327 - t445 * t527;
t529 = t349 * t445 + t353 * t448;
t538 = pkin(4) * t529 + t445 * t327 + t448 * t527;
t431 = g(1) * t445 - g(2) * t448;
t427 = qJDD(1) * pkin(1) + t431;
t432 = g(1) * t448 + g(2) * t445;
t449 = qJD(1) ^ 2;
t428 = -pkin(1) * t449 - t432;
t391 = -t427 * t447 + t428 * t444;
t386 = pkin(2) * t438 - t391;
t392 = t427 * t444 + t428 * t447;
t437 = t439 ^ 2;
t387 = -pkin(2) * t437 + t392;
t361 = -t386 * t442 + t387 * t441;
t356 = pkin(3) * t438 - t361;
t362 = t386 * t441 + t387 * t442;
t357 = -pkin(3) * t437 + t362;
t331 = -t356 * t446 + t357 * t443;
t332 = t356 * t443 + t357 * t446;
t464 = t331 * t443 + t332 * t446;
t316 = t331 * t446 - t332 * t443;
t474 = t316 * t442;
t307 = -t441 * t464 + t474;
t475 = t316 * t441;
t498 = t442 * t464 + t475;
t302 = t307 * t447 - t444 * t498;
t523 = t307 * t444 + t447 * t498;
t536 = t302 * t445 + t448 * t523;
t535 = t302 * t448 - t445 * t523;
t418 = t437 * t442 + t438 * t441;
t421 = t437 * t441 - t438 * t442;
t380 = t418 * t447 - t421 * t444;
t400 = qJ(3) * t418 - t440 * t442;
t493 = qJ(3) * t421 - t440 * t441;
t348 = pkin(5) * t380 + t400 * t447 - t444 * t493;
t384 = t418 * t444 + t421 * t447;
t491 = t380 * t448 - t384 * t445;
t510 = pkin(5) * t384 + t400 * t444 + t447 * t493;
t530 = pkin(4) * t491 + t448 * t348 - t445 * t510;
t512 = t380 * t445 + t384 * t448;
t528 = pkin(4) * t512 + t445 * t348 + t448 * t510;
t463 = t361 * t441 + t362 * t442;
t337 = t361 * t442 - t362 * t441;
t472 = t337 * t447;
t320 = -t444 * t463 + t472;
t473 = t337 * t444;
t499 = t447 * t463 + t473;
t522 = t320 * t445 + t448 * t499;
t521 = t320 * t448 - t445 * t499;
t422 = t437 * t447 + t438 * t444;
t405 = pkin(5) * t422 - g(3) * t447;
t425 = t437 * t444 - t438 * t447;
t455 = t422 * t448 - t425 * t445;
t495 = pkin(5) * t425 - g(3) * t444;
t513 = pkin(4) * t455 + t448 * t405 - t445 * t495;
t490 = t422 * t445 + t425 * t448;
t511 = pkin(4) * t490 + t445 * t405 + t448 * t495;
t462 = t391 * t444 + t392 * t447;
t367 = t391 * t447 - t392 * t444;
t470 = t367 * t448;
t501 = -t445 * t462 + t470;
t471 = t367 * t445;
t500 = t448 * t462 + t471;
t479 = pkin(1) * t440;
t478 = pkin(2) * t440;
t313 = pkin(3) * t316;
t465 = -pkin(2) * t307 - t313;
t459 = -t431 * t445 - t432 * t448;
t430 = qJDD(1) * t448 - t445 * t449;
t458 = -pkin(4) * t430 - g(3) * t445;
t457 = -pkin(3) * t413 - t331;
t456 = -pkin(2) * t421 - t361;
t454 = t431 * t448 - t432 * t445;
t453 = -pkin(2) * t377 + t457;
t452 = -pkin(3) * t410 - t332;
t451 = -pkin(2) * t418 - t362;
t450 = -pkin(2) * t373 + t452;
t429 = qJDD(1) * t445 + t448 * t449;
t416 = -pkin(4) * t429 + g(3) * t448;
t370 = -pkin(1) * t425 - t391;
t369 = -pkin(1) * t422 - t392;
t364 = pkin(1) * t367;
t363 = pkin(1) * g(3) + pkin(5) * t462;
t340 = -pkin(1) * t384 + t456;
t339 = -pkin(1) * t380 + t451;
t334 = pkin(2) * t337;
t333 = qJ(3) * t463 + t478;
t323 = -pkin(1) * t353 + t453;
t322 = -pkin(1) * t349 + t450;
t312 = pkin(3) * t440 + pkin(6) * t464;
t311 = -pkin(1) * t320 - t334;
t310 = pkin(5) * t320 + qJ(3) * t472 - t333 * t444;
t309 = pkin(5) * t499 + qJ(3) * t473 + t333 * t447 + t479;
t299 = pkin(6) * t474 + qJ(3) * t307 - t312 * t441;
t298 = pkin(6) * t475 + qJ(3) * t498 + t312 * t442 + t478;
t297 = -pkin(1) * t302 + t465;
t296 = pkin(5) * t302 - t298 * t444 + t299 * t447;
t295 = pkin(5) * t523 + t298 * t447 + t299 * t444 + t479;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t430, 0, -t429, 0, t458, -t416, -t454, -pkin(4) * t454, 0, 0, -t490, 0, -t455, 0, t511, t513, t501, pkin(4) * t501 + pkin(5) * t470 - t445 * t363, 0, 0, -t512, 0, -t491, 0, t528, t530, t521, pkin(4) * t521 - t445 * t309 + t448 * t310, 0, 0, -t529, 0, -t514, 0, t538, t539, t535, pkin(4) * t535 - t445 * t295 + t448 * t296; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t429, 0, t430, 0, t416, t458, t459, pkin(4) * t459, 0, 0, t455, 0, -t490, 0, -t513, t511, t500, pkin(4) * t500 + pkin(5) * t471 + t448 * t363, 0, 0, t491, 0, -t512, 0, -t530, t528, t522, pkin(4) * t522 + t448 * t309 + t445 * t310, 0, 0, t514, 0, -t529, 0, -t539, t538, t536, pkin(4) * t536 + t448 * t295 + t445 * t296; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t431, t432, 0, 0, 0, 0, 0, 0, 0, t438, t370, t369, 0, -t364, 0, 0, 0, 0, 0, t438, t340, t339, 0, t311, 0, 0, 0, 0, 0, t434, t323, t322, 0, t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t449, 0, 0, -g(3), -t431, 0, 0, 0, -t425, 0, -t422, 0, t495, t405, t367, pkin(5) * t367, 0, 0, -t384, 0, -t380, 0, t510, t348, t320, t310, 0, 0, -t353, 0, -t349, 0, t527, t327, t302, t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t449, 0, qJDD(1), 0, g(3), 0, -t432, 0, 0, 0, t422, 0, -t425, 0, -t405, t495, t462, t363, 0, 0, t380, 0, -t384, 0, -t348, t510, t499, t309, 0, 0, t349, 0, -t353, 0, -t327, t527, t523, t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t431, t432, 0, 0, 0, 0, 0, 0, 0, t438, t370, t369, 0, -t364, 0, 0, 0, 0, 0, t438, t340, t339, 0, t311, 0, 0, 0, 0, 0, t434, t323, t322, 0, t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t438, 0, -t437, 0, 0, -g(3), t391, 0, 0, 0, -t421, 0, -t418, 0, t493, t400, t337, qJ(3) * t337, 0, 0, -t377, 0, -t373, 0, t509, t344, t307, t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t437, 0, t438, 0, g(3), 0, t392, 0, 0, 0, t418, 0, -t421, 0, -t400, t493, t463, t333, 0, 0, t373, 0, -t377, 0, -t344, t509, t498, t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t438, -t391, -t392, 0, 0, 0, 0, 0, 0, 0, t438, t456, t451, 0, -t334, 0, 0, 0, 0, 0, t434, t453, t450, 0, t465; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t438, 0, -t437, 0, 0, -t440, t361, 0, 0, 0, -t413, 0, -t410, 0, t494, t396, t316, pkin(6) * t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t437, 0, t438, 0, t440, 0, t362, 0, 0, 0, t410, 0, -t413, 0, -t396, t494, t464, t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t438, -t361, -t362, 0, 0, 0, 0, 0, 0, 0, t434, t457, t452, 0, -t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t434, 0, -t433, 0, 0, -t440, t331, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t433, 0, t434, 0, t440, 0, t332, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t434, -t331, -t332, 0, 0;];
m_new_reg  = t1;
