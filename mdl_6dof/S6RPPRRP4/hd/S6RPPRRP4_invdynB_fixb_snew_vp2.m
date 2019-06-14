% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:54:10
% EndTime: 2019-05-05 14:54:15
% DurationCPUTime: 2.82s
% Computational Cost: add. (28698->273), mult. (50787->316), div. (0->0), fcn. (24767->8), ass. (0->110)
t566 = Ifges(6,1) + Ifges(7,1);
t556 = Ifges(6,4) - Ifges(7,5);
t565 = -Ifges(6,5) - Ifges(7,4);
t564 = -Ifges(6,2) - Ifges(7,3);
t553 = Ifges(6,6) - Ifges(7,6);
t563 = -Ifges(6,3) - Ifges(7,2);
t562 = 2 * qJD(6);
t561 = -pkin(1) - pkin(2);
t560 = sin(qJ(5));
t559 = mrSges(2,1) + mrSges(3,1);
t558 = -mrSges(6,3) - mrSges(7,2);
t557 = Ifges(3,4) + Ifges(2,5);
t554 = Ifges(2,6) - Ifges(3,6);
t524 = sin(qJ(1));
t527 = cos(qJ(1));
t507 = -t527 * g(1) - t524 * g(2);
t529 = qJD(1) ^ 2;
t535 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t507;
t488 = -t529 * pkin(1) + t535;
t485 = t561 * t529 + t535;
t506 = t524 * g(1) - t527 * g(2);
t534 = -t529 * qJ(2) + qJDD(2) - t506;
t487 = t561 * qJDD(1) + t534;
t521 = sin(pkin(9));
t522 = cos(pkin(9));
t460 = t522 * t485 + t521 * t487;
t458 = -t529 * pkin(3) - qJDD(1) * pkin(7) + t460;
t518 = g(3) + qJDD(3);
t523 = sin(qJ(4));
t526 = cos(qJ(4));
t454 = t526 * t458 + t523 * t518;
t500 = (mrSges(5,1) * t526 - mrSges(5,2) * t523) * qJD(1);
t545 = qJD(1) * qJD(4);
t543 = t523 * t545;
t503 = -t526 * qJDD(1) + t543;
t547 = qJD(1) * t523;
t504 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t547;
t459 = -t521 * t485 + t522 * t487;
t457 = qJDD(1) * pkin(3) - t529 * pkin(7) - t459;
t542 = t526 * t545;
t502 = -t523 * qJDD(1) - t542;
t449 = (-t502 + t542) * pkin(8) + (-t503 - t543) * pkin(4) + t457;
t501 = (pkin(4) * t526 + pkin(8) * t523) * qJD(1);
t528 = qJD(4) ^ 2;
t546 = t526 * qJD(1);
t452 = -t528 * pkin(4) + qJDD(4) * pkin(8) - t501 * t546 + t454;
t525 = cos(qJ(5));
t447 = t560 * t449 + t525 * t452;
t499 = -t560 * qJD(4) + t525 * t547;
t471 = -t499 * qJD(5) - t525 * qJDD(4) + t560 * t502;
t509 = qJD(5) + t546;
t482 = t509 * mrSges(6,1) + t499 * mrSges(6,3);
t497 = qJDD(5) - t503;
t498 = t525 * qJD(4) + t560 * t547;
t475 = -t498 * pkin(5) + t499 * qJ(6);
t508 = t509 ^ 2;
t443 = -t508 * pkin(5) + t497 * qJ(6) + t498 * t475 + t509 * t562 + t447;
t483 = -t509 * mrSges(7,1) - t499 * mrSges(7,2);
t544 = m(7) * t443 + t497 * mrSges(7,3) + t509 * t483;
t476 = -t498 * mrSges(7,1) + t499 * mrSges(7,3);
t548 = -t498 * mrSges(6,1) - t499 * mrSges(6,2) + t476;
t439 = m(6) * t447 - t497 * mrSges(6,2) + t558 * t471 - t509 * t482 + t548 * t498 + t544;
t446 = t525 * t449 - t560 * t452;
t472 = t498 * qJD(5) + t560 * qJDD(4) + t525 * t502;
t481 = -t509 * mrSges(6,2) + t498 * mrSges(6,3);
t444 = -t497 * pkin(5) - t508 * qJ(6) - t499 * t475 + qJDD(6) - t446;
t484 = t498 * mrSges(7,2) + t509 * mrSges(7,3);
t538 = -m(7) * t444 + t497 * mrSges(7,1) + t509 * t484;
t440 = m(6) * t446 + t497 * mrSges(6,1) + t558 * t472 + t509 * t481 + t548 * t499 + t538;
t537 = t525 * t439 - t560 * t440;
t433 = m(5) * t454 - qJDD(4) * mrSges(5,2) + t503 * mrSges(5,3) - qJD(4) * t504 - t500 * t546 + t537;
t453 = -t523 * t458 + t526 * t518;
t505 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t546;
t451 = -qJDD(4) * pkin(4) - t528 * pkin(8) - t501 * t547 - t453;
t445 = t499 * t562 + (-t498 * t509 - t472) * qJ(6) + (-t499 * t509 + t471) * pkin(5) + t451;
t441 = m(7) * t445 + t471 * mrSges(7,1) - t472 * mrSges(7,3) + t499 * t483 - t498 * t484;
t530 = -m(6) * t451 - t471 * mrSges(6,1) - t472 * mrSges(6,2) + t498 * t481 + t499 * t482 - t441;
t437 = m(5) * t453 + qJDD(4) * mrSges(5,1) - t502 * mrSges(5,3) + qJD(4) * t505 + t500 * t547 + t530;
t539 = t526 * t433 - t523 * t437;
t428 = m(4) * t460 - t529 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t539;
t436 = t560 * t439 + t525 * t440;
t531 = -m(5) * t457 + t503 * mrSges(5,1) - t502 * mrSges(5,2) + t504 * t547 - t505 * t546 - t436;
t431 = m(4) * t459 - qJDD(1) * mrSges(4,1) - t529 * mrSges(4,2) + t531;
t540 = t522 * t428 - t521 * t431;
t536 = m(3) * t488 + qJDD(1) * mrSges(3,3) + t540;
t423 = m(2) * t507 - qJDD(1) * mrSges(2,2) - t559 * t529 + t536;
t425 = t521 * t428 + t522 * t431;
t489 = -qJDD(1) * pkin(1) + t534;
t532 = -m(3) * t489 + qJDD(1) * mrSges(3,1) + t529 * mrSges(3,3) - t425;
t424 = m(2) * t506 + qJDD(1) * mrSges(2,1) - t529 * mrSges(2,2) + t532;
t552 = t524 * t423 + t527 * t424;
t551 = t498 * t564 + t499 * t556 - t509 * t553;
t550 = t498 * t553 + t499 * t565 - t509 * t563;
t549 = t556 * t498 - t499 * t566 - t565 * t509;
t541 = t527 * t423 - t524 * t424;
t430 = t523 * t433 + t526 * t437;
t533 = -m(4) * t518 - t430;
t492 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t523 - Ifges(5,4) * t526) * qJD(1);
t491 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t523 - Ifges(5,2) * t526) * qJD(1);
t490 = (Ifges(5,3) * qJD(4)) + (-Ifges(5,5) * t523 - Ifges(5,6) * t526) * qJD(1);
t435 = mrSges(6,2) * t451 + mrSges(7,2) * t444 - mrSges(6,3) * t446 - mrSges(7,3) * t445 - qJ(6) * t441 - t556 * t471 + t472 * t566 - t497 * t565 + t550 * t498 + t551 * t509;
t434 = -mrSges(6,1) * t451 - mrSges(7,1) * t445 + mrSges(7,2) * t443 + mrSges(6,3) * t447 - pkin(5) * t441 + t471 * t564 + t556 * t472 + t553 * t497 + t550 * t499 + t549 * t509;
t429 = -m(3) * g(3) + t533;
t426 = Ifges(5,4) * t502 + Ifges(5,2) * t503 + Ifges(5,6) * qJDD(4) + t490 * t547 + qJD(4) * t492 - mrSges(5,1) * t457 + mrSges(5,3) * t454 - mrSges(6,1) * t446 + mrSges(6,2) * t447 + mrSges(7,1) * t444 - mrSges(7,3) * t443 - pkin(5) * t538 - qJ(6) * t544 - pkin(4) * t436 + (-pkin(5) * t476 - t551) * t499 + (-qJ(6) * t476 + t549) * t498 + t563 * t497 + (pkin(5) * mrSges(7,2) + t565) * t472 + (qJ(6) * mrSges(7,2) + t553) * t471;
t419 = mrSges(5,2) * t457 - mrSges(5,3) * t453 + Ifges(5,1) * t502 + Ifges(5,4) * t503 + Ifges(5,5) * qJDD(4) - pkin(8) * t436 - qJD(4) * t491 - t560 * t434 + t525 * t435 - t490 * t546;
t418 = -Ifges(4,6) * qJDD(1) + t529 * Ifges(4,5) - mrSges(4,1) * t518 + mrSges(4,3) * t460 - Ifges(5,5) * t502 - Ifges(5,6) * t503 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t453 + mrSges(5,2) * t454 - t560 * t435 - t525 * t434 - pkin(4) * t530 - pkin(8) * t537 - pkin(3) * t430 + (t523 * t491 - t526 * t492) * qJD(1);
t417 = mrSges(4,2) * t518 - mrSges(4,3) * t459 - Ifges(4,5) * qJDD(1) - t529 * Ifges(4,6) - pkin(7) * t430 + t526 * t419 - t523 * t426;
t416 = mrSges(3,2) * t489 - mrSges(2,3) * t506 - qJ(2) * t429 - qJ(3) * t425 + t522 * t417 - t521 * t418 - t554 * t529 + t557 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t415 = mrSges(3,2) * t488 + mrSges(2,3) * t507 - pkin(1) * t429 - pkin(2) * t533 + t559 * g(3) - qJ(3) * t540 + t554 * qJDD(1) - t521 * t417 - t522 * t418 + t557 * t529;
t1 = [-m(1) * g(1) + t541; -m(1) * g(2) + t552; (-m(1) - m(2) - m(3)) * g(3) + t533; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t552 - t524 * t415 + t527 * t416; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t541 + t527 * t415 + t524 * t416; qJ(2) * (-t529 * mrSges(3,1) + t536) + mrSges(2,1) * t506 + pkin(1) * t532 - mrSges(2,2) * t507 - pkin(2) * t425 + mrSges(3,3) * t488 - mrSges(3,1) * t489 - pkin(7) * t539 - t523 * t419 - t526 * t426 - pkin(3) * t531 - mrSges(4,1) * t459 + mrSges(4,2) * t460 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB  = t1;
