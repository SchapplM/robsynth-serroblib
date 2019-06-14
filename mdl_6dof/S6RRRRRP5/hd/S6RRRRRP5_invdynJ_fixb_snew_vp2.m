% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 05:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:58:10
% EndTime: 2019-05-08 04:58:31
% DurationCPUTime: 8.55s
% Computational Cost: add. (87951->326), mult. (180469->397), div. (0->0), fcn. (132171->10), ass. (0->129)
t534 = sin(qJ(3));
t539 = cos(qJ(3));
t535 = sin(qJ(2));
t561 = qJD(1) * t535;
t519 = qJD(2) * t534 + t539 * t561;
t540 = cos(qJ(2));
t559 = qJD(1) * qJD(2);
t556 = t540 * t559;
t522 = qJDD(1) * t535 + t556;
t491 = -qJD(3) * t519 + qJDD(2) * t539 - t522 * t534;
t518 = qJD(2) * t539 - t534 * t561;
t492 = qJD(3) * t518 + qJDD(2) * t534 + t522 * t539;
t533 = sin(qJ(4));
t538 = cos(qJ(4));
t495 = t518 * t533 + t519 * t538;
t457 = -qJD(4) * t495 + t491 * t538 - t492 * t533;
t494 = t518 * t538 - t519 * t533;
t458 = qJD(4) * t494 + t491 * t533 + t492 * t538;
t532 = sin(qJ(5));
t537 = cos(qJ(5));
t473 = t494 * t537 - t495 * t532;
t431 = qJD(5) * t473 + t457 * t532 + t458 * t537;
t474 = t494 * t532 + t495 * t537;
t450 = -mrSges(7,1) * t473 + mrSges(7,2) * t474;
t543 = qJD(1) ^ 2;
t536 = sin(qJ(1));
t541 = cos(qJ(1));
t555 = t536 * g(1) - t541 * g(2);
t513 = -qJDD(1) * pkin(1) - t543 * pkin(7) - t555;
t529 = t535 * t559;
t523 = qJDD(1) * t540 - t529;
t478 = (-t522 - t556) * pkin(8) + (-t523 + t529) * pkin(2) + t513;
t550 = -g(1) * t541 - g(2) * t536;
t514 = -pkin(1) * t543 + qJDD(1) * pkin(7) + t550;
t500 = -g(3) * t535 + t540 * t514;
t521 = (-pkin(2) * t540 - pkin(8) * t535) * qJD(1);
t542 = qJD(2) ^ 2;
t560 = qJD(1) * t540;
t481 = -pkin(2) * t542 + qJDD(2) * pkin(8) + t521 * t560 + t500;
t459 = t539 * t478 - t534 * t481;
t517 = qJDD(3) - t523;
t528 = qJD(3) - t560;
t437 = (t518 * t528 - t492) * pkin(9) + (t518 * t519 + t517) * pkin(3) + t459;
t460 = t534 * t478 + t539 * t481;
t501 = pkin(3) * t528 - pkin(9) * t519;
t516 = t518 ^ 2;
t439 = -pkin(3) * t516 + pkin(9) * t491 - t501 * t528 + t460;
t417 = t538 * t437 - t533 * t439;
t515 = qJDD(4) + t517;
t527 = qJD(4) + t528;
t413 = (t494 * t527 - t458) * pkin(10) + (t494 * t495 + t515) * pkin(4) + t417;
t418 = t533 * t437 + t538 * t439;
t484 = pkin(4) * t527 - pkin(10) * t495;
t493 = t494 ^ 2;
t415 = -pkin(4) * t493 + pkin(10) * t457 - t484 * t527 + t418;
t407 = t537 * t413 - t532 * t415;
t509 = qJDD(5) + t515;
t524 = qJD(5) + t527;
t402 = -0.2e1 * qJD(6) * t474 + (t473 * t524 - t431) * qJ(6) + (t473 * t474 + t509) * pkin(5) + t407;
t462 = -mrSges(7,2) * t524 + mrSges(7,3) * t473;
t558 = m(7) * t402 + t509 * mrSges(7,1) + t524 * t462;
t399 = -t431 * mrSges(7,3) - t474 * t450 + t558;
t408 = t532 * t413 + t537 * t415;
t430 = -qJD(5) * t474 + t457 * t537 - t458 * t532;
t464 = pkin(5) * t524 - qJ(6) * t474;
t472 = t473 ^ 2;
t405 = -pkin(5) * t472 + qJ(6) * t430 + 0.2e1 * qJD(6) * t473 - t464 * t524 + t408;
t564 = Ifges(6,4) + Ifges(7,4);
t569 = Ifges(6,5) + Ifges(7,5);
t570 = Ifges(6,1) + Ifges(7,1);
t565 = t564 * t473 + t474 * t570 + t569 * t524;
t568 = Ifges(6,6) + Ifges(7,6);
t572 = Ifges(6,2) + Ifges(7,2);
t566 = t572 * t473 + t474 * t564 + t568 * t524;
t567 = Ifges(6,3) + Ifges(7,3);
t574 = mrSges(6,1) * t407 + mrSges(7,1) * t402 - mrSges(6,2) * t408 - mrSges(7,2) * t405 + pkin(5) * t399 + t568 * t430 + t431 * t569 - t565 * t473 + t566 * t474 + t567 * t509;
t451 = -mrSges(6,1) * t473 + mrSges(6,2) * t474;
t463 = -mrSges(6,2) * t524 + mrSges(6,3) * t473;
t394 = m(6) * t407 + t509 * mrSges(6,1) + t524 * t463 + (-t450 - t451) * t474 + (-mrSges(6,3) - mrSges(7,3)) * t431 + t558;
t465 = mrSges(7,1) * t524 - mrSges(7,3) * t474;
t466 = mrSges(6,1) * t524 - mrSges(6,3) * t474;
t557 = m(7) * t405 + t430 * mrSges(7,3) + t473 * t450;
t397 = m(6) * t408 + t430 * mrSges(6,3) + t473 * t451 + (-t465 - t466) * t524 + (-mrSges(6,2) - mrSges(7,2)) * t509 + t557;
t392 = t537 * t394 + t532 * t397;
t468 = Ifges(5,4) * t495 + Ifges(5,2) * t494 + Ifges(5,6) * t527;
t469 = Ifges(5,1) * t495 + Ifges(5,4) * t494 + Ifges(5,5) * t527;
t573 = mrSges(5,1) * t417 - mrSges(5,2) * t418 + Ifges(5,5) * t458 + Ifges(5,6) * t457 + Ifges(5,3) * t515 + pkin(4) * t392 + t495 * t468 - t494 * t469 + t574;
t475 = -mrSges(5,1) * t494 + mrSges(5,2) * t495;
t482 = -mrSges(5,2) * t527 + mrSges(5,3) * t494;
t388 = m(5) * t417 + mrSges(5,1) * t515 - mrSges(5,3) * t458 - t475 * t495 + t482 * t527 + t392;
t483 = mrSges(5,1) * t527 - mrSges(5,3) * t495;
t551 = -t394 * t532 + t537 * t397;
t389 = m(5) * t418 - mrSges(5,2) * t515 + mrSges(5,3) * t457 + t475 * t494 - t483 * t527 + t551;
t383 = t538 * t388 + t533 * t389;
t486 = Ifges(4,4) * t519 + Ifges(4,2) * t518 + Ifges(4,6) * t528;
t487 = Ifges(4,1) * t519 + Ifges(4,4) * t518 + Ifges(4,5) * t528;
t571 = mrSges(4,1) * t459 - mrSges(4,2) * t460 + Ifges(4,5) * t492 + Ifges(4,6) * t491 + Ifges(4,3) * t517 + pkin(3) * t383 + t519 * t486 - t518 * t487 + t573;
t563 = -t568 * t473 - t474 * t569 - t567 * t524;
t499 = -t540 * g(3) - t535 * t514;
t496 = -mrSges(4,1) * t518 + mrSges(4,2) * t519;
t497 = -mrSges(4,2) * t528 + mrSges(4,3) * t518;
t381 = m(4) * t459 + mrSges(4,1) * t517 - mrSges(4,3) * t492 - t496 * t519 + t497 * t528 + t383;
t498 = mrSges(4,1) * t528 - mrSges(4,3) * t519;
t552 = -t388 * t533 + t538 * t389;
t382 = m(4) * t460 - mrSges(4,2) * t517 + mrSges(4,3) * t491 + t496 * t518 - t498 * t528 + t552;
t553 = -t381 * t534 + t539 * t382;
t480 = -qJDD(2) * pkin(2) - pkin(8) * t542 + t521 * t561 - t499;
t452 = -pkin(3) * t491 - pkin(9) * t516 + t519 * t501 + t480;
t420 = -pkin(4) * t457 - pkin(10) * t493 + t495 * t484 + t452;
t410 = -pkin(5) * t430 - qJ(6) * t472 + t464 * t474 + qJDD(6) + t420;
t403 = m(7) * t410 - t430 * mrSges(7,1) + t431 * mrSges(7,2) - t473 * t462 + t474 * t465;
t377 = t539 * t381 + t534 * t382;
t549 = m(6) * t420 - t430 * mrSges(6,1) + t431 * mrSges(6,2) - t473 * t463 + t474 * t466 + t403;
t547 = m(5) * t452 - t457 * mrSges(5,1) + t458 * mrSges(5,2) - t494 * t482 + t495 * t483 + t549;
t545 = -m(4) * t480 + t491 * mrSges(4,1) - t492 * mrSges(4,2) + t518 * t497 - t519 * t498 - t547;
t526 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t560;
t525 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t561;
t520 = (-mrSges(3,1) * t540 + mrSges(3,2) * t535) * qJD(1);
t512 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t535 + Ifges(3,4) * t540) * qJD(1);
t511 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t535 + Ifges(3,2) * t540) * qJD(1);
t485 = Ifges(4,5) * t519 + Ifges(4,6) * t518 + Ifges(4,3) * t528;
t467 = Ifges(5,5) * t495 + Ifges(5,6) * t494 + Ifges(5,3) * t527;
t390 = mrSges(6,2) * t420 + mrSges(7,2) * t410 - mrSges(6,3) * t407 - mrSges(7,3) * t402 - qJ(6) * t399 + t564 * t430 + t431 * t570 - t563 * t473 + t569 * t509 - t566 * t524;
t384 = -mrSges(6,1) * t420 + mrSges(6,3) * t408 - mrSges(7,1) * t410 + mrSges(7,3) * t405 - pkin(5) * t403 + qJ(6) * t557 + (-qJ(6) * t465 + t565) * t524 + (-mrSges(7,2) * qJ(6) + t568) * t509 + t563 * t474 + t564 * t431 + t572 * t430;
t379 = mrSges(5,2) * t452 - mrSges(5,3) * t417 + Ifges(5,1) * t458 + Ifges(5,4) * t457 + Ifges(5,5) * t515 - pkin(10) * t392 - t384 * t532 + t390 * t537 + t467 * t494 - t468 * t527;
t378 = -mrSges(5,1) * t452 + mrSges(5,3) * t418 + Ifges(5,4) * t458 + Ifges(5,2) * t457 + Ifges(5,6) * t515 - pkin(4) * t549 + pkin(10) * t551 + t537 * t384 + t532 * t390 - t495 * t467 + t527 * t469;
t376 = mrSges(4,2) * t480 - mrSges(4,3) * t459 + Ifges(4,1) * t492 + Ifges(4,4) * t491 + Ifges(4,5) * t517 - pkin(9) * t383 - t378 * t533 + t379 * t538 + t485 * t518 - t486 * t528;
t375 = -mrSges(4,1) * t480 + mrSges(4,3) * t460 + Ifges(4,4) * t492 + Ifges(4,2) * t491 + Ifges(4,6) * t517 - pkin(3) * t547 + pkin(9) * t552 + t538 * t378 + t533 * t379 - t519 * t485 + t528 * t487;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t555 - mrSges(2,2) * t550 + t535 * (mrSges(3,2) * t513 - mrSges(3,3) * t499 + Ifges(3,1) * t522 + Ifges(3,4) * t523 + Ifges(3,5) * qJDD(2) - pkin(8) * t377 - qJD(2) * t511 - t534 * t375 + t539 * t376) + t540 * (-mrSges(3,1) * t513 + mrSges(3,3) * t500 + Ifges(3,4) * t522 + Ifges(3,2) * t523 + Ifges(3,6) * qJDD(2) - pkin(2) * t377 + qJD(2) * t512 - t571) + pkin(1) * (-m(3) * t513 + t523 * mrSges(3,1) - t522 * mrSges(3,2) + (-t525 * t535 + t526 * t540) * qJD(1) - t377) + pkin(7) * (t540 * (m(3) * t500 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t523 - qJD(2) * t525 + t520 * t560 + t553) - t535 * (m(3) * t499 + qJDD(2) * mrSges(3,1) - t522 * mrSges(3,3) + qJD(2) * t526 - t520 * t561 + t545)); Ifges(3,5) * t522 + Ifges(3,6) * t523 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t499 - mrSges(3,2) * t500 + t534 * t376 + t539 * t375 + pkin(2) * t545 + pkin(8) * t553 + (t535 * t511 - t540 * t512) * qJD(1); t571; t573; t574; t403;];
tauJ  = t1;
