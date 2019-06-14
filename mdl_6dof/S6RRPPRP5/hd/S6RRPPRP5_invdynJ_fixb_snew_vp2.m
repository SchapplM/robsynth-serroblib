% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-05-06 09:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:25:27
% EndTime: 2019-05-06 09:25:32
% DurationCPUTime: 3.15s
% Computational Cost: add. (19268->310), mult. (41976->366), div. (0->0), fcn. (25188->8), ass. (0->123)
t573 = Ifges(6,1) + Ifges(7,1);
t557 = Ifges(6,4) - Ifges(7,5);
t569 = Ifges(7,4) + Ifges(6,5);
t572 = Ifges(6,2) + Ifges(7,3);
t567 = Ifges(6,6) - Ifges(7,6);
t571 = -2 * qJD(3);
t570 = Ifges(3,1) + Ifges(4,2);
t558 = Ifges(3,4) + Ifges(4,6);
t556 = Ifges(3,5) - Ifges(4,4);
t568 = Ifges(3,2) + Ifges(4,3);
t555 = Ifges(3,6) - Ifges(4,5);
t566 = Ifges(3,3) + Ifges(4,1);
t565 = Ifges(6,3) + Ifges(7,2);
t520 = sin(pkin(9));
t521 = cos(pkin(9));
t525 = cos(qJ(2));
t547 = qJD(1) * t525;
t493 = -qJD(2) * t520 - t521 * t547;
t494 = qJD(2) * t521 - t520 * t547;
t522 = sin(qJ(5));
t562 = cos(qJ(5));
t460 = -t493 * t562 + t494 * t522;
t461 = t522 * t493 + t494 * t562;
t523 = sin(qJ(2));
t548 = qJD(1) * t523;
t511 = qJD(5) + t548;
t564 = t572 * t460 - t557 * t461 - t567 * t511;
t563 = -t557 * t460 + t573 * t461 + t569 * t511;
t528 = qJD(1) ^ 2;
t524 = sin(qJ(1));
t526 = cos(qJ(1));
t537 = -g(1) * t526 - g(2) * t524;
t487 = -pkin(1) * t528 + qJDD(1) * pkin(7) + t537;
t465 = -t523 * g(3) + t525 * t487;
t497 = (-pkin(2) * t525 - qJ(3) * t523) * qJD(1);
t527 = qJD(2) ^ 2;
t449 = t527 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t571 - t497 * t547 - t465;
t561 = t528 * pkin(7);
t560 = mrSges(3,1) - mrSges(4,2);
t559 = -mrSges(6,3) - mrSges(7,2);
t546 = qJD(1) * qJD(2);
t542 = t523 * t546;
t501 = qJDD(1) * t525 - t542;
t506 = pkin(3) * t548 - qJD(2) * qJ(4);
t519 = t525 ^ 2;
t543 = t525 * t546;
t500 = qJDD(1) * t523 + t543;
t541 = t524 * g(1) - t526 * g(2);
t535 = -qJDD(1) * pkin(1) - t541;
t532 = pkin(2) * t542 + t548 * t571 + (-t500 - t543) * qJ(3) + t535;
t421 = -t506 * t548 + (-pkin(3) * t519 - pkin(7)) * t528 + (-pkin(2) - qJ(4)) * t501 + t532;
t464 = -t525 * g(3) - t523 * t487;
t450 = -qJDD(2) * pkin(2) - t527 * qJ(3) + t497 * t548 + qJDD(3) - t464;
t440 = (-t523 * t525 * t528 - qJDD(2)) * qJ(4) + (t500 - t543) * pkin(3) + t450;
t414 = -0.2e1 * qJD(4) * t494 - t520 * t421 + t521 * t440;
t470 = qJDD(2) * t521 - t501 * t520;
t411 = (t493 * t548 - t470) * pkin(8) + (t493 * t494 + t500) * pkin(4) + t414;
t415 = 0.2e1 * qJD(4) * t493 + t521 * t421 + t520 * t440;
t469 = -qJDD(2) * t520 - t501 * t521;
t471 = pkin(4) * t548 - pkin(8) * t494;
t492 = t493 ^ 2;
t413 = -pkin(4) * t492 + pkin(8) * t469 - t471 * t548 + t415;
t407 = t522 * t411 + t413 * t562;
t428 = qJD(5) * t461 - t469 * t562 + t470 * t522;
t453 = mrSges(6,1) * t511 - mrSges(6,3) * t461;
t496 = qJDD(5) + t500;
t443 = pkin(5) * t460 - qJ(6) * t461;
t509 = t511 ^ 2;
t403 = -pkin(5) * t509 + qJ(6) * t496 + 0.2e1 * qJD(6) * t511 - t443 * t460 + t407;
t454 = -mrSges(7,1) * t511 + mrSges(7,2) * t461;
t544 = m(7) * t403 + t496 * mrSges(7,3) + t511 * t454;
t444 = mrSges(7,1) * t460 - mrSges(7,3) * t461;
t553 = -mrSges(6,1) * t460 - mrSges(6,2) * t461 - t444;
t393 = m(6) * t407 - t496 * mrSges(6,2) + t428 * t559 - t511 * t453 + t460 * t553 + t544;
t406 = t411 * t562 - t522 * t413;
t429 = -t460 * qJD(5) + t522 * t469 + t470 * t562;
t451 = -mrSges(6,2) * t511 - mrSges(6,3) * t460;
t404 = -t496 * pkin(5) - t509 * qJ(6) + t461 * t443 + qJDD(6) - t406;
t452 = -mrSges(7,2) * t460 + mrSges(7,3) * t511;
t538 = -m(7) * t404 + t496 * mrSges(7,1) + t511 * t452;
t395 = m(6) * t406 + t496 * mrSges(6,1) + t429 * t559 + t511 * t451 + t461 * t553 + t538;
t388 = t522 * t393 + t395 * t562;
t554 = t460 * t567 - t461 * t569 - t511 * t565;
t552 = t566 * qJD(2) + (t523 * t556 + t525 * t555) * qJD(1);
t551 = t555 * qJD(2) + (t523 * t558 + t525 * t568) * qJD(1);
t550 = t556 * qJD(2) + (t523 * t570 + t525 * t558) * qJD(1);
t507 = -mrSges(4,1) * t547 - qJD(2) * mrSges(4,3);
t549 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t547 - t507;
t462 = -mrSges(5,1) * t493 + mrSges(5,2) * t494;
t467 = -mrSges(5,2) * t548 + mrSges(5,3) * t493;
t386 = m(5) * t414 + mrSges(5,1) * t500 - mrSges(5,3) * t470 - t462 * t494 + t467 * t548 + t388;
t468 = mrSges(5,1) * t548 - mrSges(5,3) * t494;
t539 = t393 * t562 - t395 * t522;
t387 = m(5) * t415 - mrSges(5,2) * t500 + mrSges(5,3) * t469 + t462 * t493 - t468 * t548 + t539;
t540 = -t520 * t386 + t521 * t387;
t384 = t521 * t386 + t520 * t387;
t446 = -t501 * pkin(2) + t532 - t561;
t536 = m(4) * t446 + t540;
t433 = -t519 * t528 * qJ(4) + t501 * pkin(3) + qJD(2) * t506 + qJDD(4) - t449;
t417 = -t469 * pkin(4) - t492 * pkin(8) + t494 * t471 + t433;
t409 = (t461 * t511 + t428) * pkin(5) + (t460 * t511 - t429) * qJ(6) - 0.2e1 * qJD(6) * t461 + t417;
t400 = m(7) * t409 + t428 * mrSges(7,1) - t429 * mrSges(7,3) + t460 * t452 - t461 * t454;
t533 = m(4) * t450 + t500 * mrSges(4,1) + t384;
t531 = m(6) * t417 + t428 * mrSges(6,1) + t429 * mrSges(6,2) + t460 * t451 + t461 * t453 + t400;
t399 = t429 * mrSges(7,2) + t461 * t444 - t538;
t530 = mrSges(6,1) * t406 - mrSges(7,1) * t404 - mrSges(6,2) * t407 + mrSges(7,3) * t403 - pkin(5) * t399 + qJ(6) * t544 + t565 * t496 - t564 * t461 + (-qJ(6) * t444 + t563) * t460 + t569 * t429 + (-qJ(6) * mrSges(7,2) - t567) * t428;
t396 = m(5) * t433 - t469 * mrSges(5,1) + t470 * mrSges(5,2) - t493 * t467 + t494 * t468 + t531;
t498 = (mrSges(4,2) * t525 - mrSges(4,3) * t523) * qJD(1);
t508 = mrSges(4,1) * t548 + qJD(2) * mrSges(4,2);
t529 = -m(4) * t449 + qJDD(2) * mrSges(4,3) + qJD(2) * t508 + t498 * t547 + t396;
t504 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t548;
t499 = (-mrSges(3,1) * t525 + mrSges(3,2) * t523) * qJD(1);
t486 = t535 - t561;
t457 = Ifges(5,1) * t494 + Ifges(5,4) * t493 + Ifges(5,5) * t548;
t456 = Ifges(5,4) * t494 + Ifges(5,2) * t493 + Ifges(5,6) * t548;
t455 = Ifges(5,5) * t494 + Ifges(5,6) * t493 + Ifges(5,3) * t548;
t390 = mrSges(6,2) * t417 + mrSges(7,2) * t404 - mrSges(6,3) * t406 - mrSges(7,3) * t409 - qJ(6) * t400 - t557 * t428 + t573 * t429 + t554 * t460 + t569 * t496 + t564 * t511;
t389 = -mrSges(6,1) * t417 - mrSges(7,1) * t409 + mrSges(7,2) * t403 + mrSges(6,3) * t407 - pkin(5) * t400 - t572 * t428 + t557 * t429 + t554 * t461 + t567 * t496 + t563 * t511;
t383 = qJDD(2) * mrSges(4,2) + qJD(2) * t507 + t498 * t548 + t533;
t382 = t501 * mrSges(4,2) - t500 * mrSges(4,3) + (t507 * t525 - t508 * t523) * qJD(1) + t536;
t381 = mrSges(5,2) * t433 - mrSges(5,3) * t414 + Ifges(5,1) * t470 + Ifges(5,4) * t469 + Ifges(5,5) * t500 - pkin(8) * t388 - t522 * t389 + t390 * t562 + t493 * t455 - t456 * t548;
t380 = -mrSges(5,1) * t433 + mrSges(5,3) * t415 + Ifges(5,4) * t470 + Ifges(5,2) * t469 + Ifges(5,6) * t500 - pkin(4) * t531 + pkin(8) * t539 + t389 * t562 + t522 * t390 - t494 * t455 + t457 * t548;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t541 - mrSges(2,2) * t537 + t523 * (t530 + t552 * t547 + t556 * qJDD(2) + (Ifges(5,3) + t570) * t500 - t493 * t457 + t494 * t456 + Ifges(5,5) * t470 + mrSges(3,2) * t486 - mrSges(3,3) * t464 + Ifges(5,6) * t469 + mrSges(4,1) * t450 - mrSges(4,3) * t446 + mrSges(5,1) * t414 - mrSges(5,2) * t415 + pkin(4) * t388 + pkin(3) * t384 - qJ(3) * t382 - t551 * qJD(2) + t558 * t501) + t525 * (-mrSges(3,1) * t486 - mrSges(4,1) * t449 + mrSges(4,2) * t446 + mrSges(3,3) * t465 - pkin(2) * t382 + pkin(3) * t396 - qJ(4) * t540 + t550 * qJD(2) + t555 * qJDD(2) - t521 * t380 - t520 * t381 + t558 * t500 + t568 * t501 - t552 * t548) + pkin(1) * (-m(3) * t486 + t560 * t501 + (-mrSges(3,2) + mrSges(4,3)) * t500 + (t549 * t525 + (-t504 + t508) * t523) * qJD(1) - t536) + pkin(7) * (t525 * (t499 * t547 + t529 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t501 - qJD(2) * t504 + m(3) * t465) + (-m(3) * t464 + t500 * mrSges(3,3) - t560 * qJDD(2) - t549 * qJD(2) + (t498 + t499) * t548 + t533) * t523); mrSges(3,1) * t464 - mrSges(3,2) * t465 + mrSges(4,2) * t450 - mrSges(4,3) * t449 + t521 * t381 - t520 * t380 - qJ(4) * t384 - pkin(2) * t383 + qJ(3) * t529 + (mrSges(4,1) * qJ(3) + t555) * t501 + t556 * t500 + t566 * qJDD(2) + (t551 * t523 - t550 * t525) * qJD(1); t383; t396; t530; t399;];
tauJ  = t1;
