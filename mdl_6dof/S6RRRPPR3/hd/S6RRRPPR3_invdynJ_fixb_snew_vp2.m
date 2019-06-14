% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-05-07 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:31:10
% EndTime: 2019-05-07 04:31:15
% DurationCPUTime: 3.20s
% Computational Cost: add. (17264->316), mult. (36173->362), div. (0->0), fcn. (23672->8), ass. (0->126)
t641 = -Ifges(6,4) + Ifges(5,5) - Ifges(4,4);
t642 = Ifges(5,1) + Ifges(4,1) + Ifges(6,2);
t635 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t640 = Ifges(6,1) + Ifges(5,3) + Ifges(4,2);
t634 = Ifges(6,5) + Ifges(4,6) - Ifges(5,6);
t584 = sin(qJ(3));
t588 = cos(qJ(2));
t610 = qJD(1) * t588;
t585 = sin(qJ(2));
t611 = qJD(1) * t585;
t624 = cos(qJ(3));
t551 = t584 * t611 - t610 * t624;
t552 = (t584 * t588 + t585 * t624) * qJD(1);
t578 = qJD(2) + qJD(3);
t636 = t551 * t641 + t552 * t642 + t578 * t635;
t633 = -Ifges(5,2) - Ifges(4,3) - Ifges(6,3);
t632 = t551 * t640 + t552 * t641 - t578 * t634;
t609 = qJD(1) * qJD(2);
t559 = qJDD(1) * t585 + t588 * t609;
t560 = qJDD(1) * t588 - t585 * t609;
t505 = qJD(3) * t552 + t559 * t584 - t560 * t624;
t506 = -qJD(3) * t551 + t559 * t624 + t560 * t584;
t590 = qJD(1) ^ 2;
t586 = sin(qJ(1));
t589 = cos(qJ(1));
t612 = g(1) * t586 - g(2) * t589;
t553 = -qJDD(1) * pkin(1) - pkin(7) * t590 - t612;
t563 = qJD(2) * pkin(2) - pkin(8) * t611;
t618 = t588 ^ 2 * t590;
t507 = -pkin(2) * t560 - pkin(8) * t618 + t563 * t611 + t553;
t619 = t551 * t578;
t600 = t505 * pkin(3) + t507 + (-t506 + t619) * qJ(4);
t623 = pkin(3) * t578;
t628 = -2 * qJD(4);
t464 = (t628 + t623) * t552 + t600;
t541 = -mrSges(5,2) * t551 + mrSges(5,3) * t578;
t631 = m(5) * t464 + mrSges(5,1) * t505 + t541 * t551;
t535 = -mrSges(6,1) * t578 - mrSges(6,3) * t551;
t536 = -mrSges(4,2) * t578 - mrSges(4,3) * t551;
t538 = mrSges(4,1) * t578 - mrSges(4,3) * t552;
t539 = -mrSges(5,1) * t578 + mrSges(5,2) * t552;
t537 = -pkin(4) * t578 - qJ(5) * t552;
t547 = t551 ^ 2;
t627 = 2 * qJD(4);
t593 = -qJ(5) * t547 + qJDD(5) - t600 + (t537 + t627) * t552;
t625 = -pkin(4) - pkin(9);
t455 = t593 + (-pkin(5) * t551 + (-pkin(3) - pkin(9)) * t552) * t578 + t625 * t505 + pkin(5) * t506;
t526 = pkin(5) * t552 - pkin(9) * t551;
t577 = qJDD(2) + qJDD(3);
t605 = -g(1) * t589 - g(2) * t586;
t554 = -pkin(1) * t590 + qJDD(1) * pkin(7) + t605;
t617 = t585 * t554;
t489 = qJDD(2) * pkin(2) - t559 * pkin(8) - t617 + (pkin(2) * t585 * t590 + pkin(8) * t609 - g(3)) * t588;
t534 = -g(3) * t585 + t554 * t588;
t490 = -pkin(2) * t618 + pkin(8) * t560 - qJD(2) * t563 + t534;
t470 = t489 * t624 - t490 * t584;
t523 = pkin(3) * t551 - qJ(4) * t552;
t629 = t578 ^ 2;
t601 = -qJ(4) * t629 + t523 * t552 + qJDD(4) - t470;
t626 = -2 * qJD(5);
t592 = (-t506 - t619) * qJ(5) + t601 + (pkin(4) * t551 + t626) * t552;
t456 = -t629 * pkin(5) - t552 * t526 + (-pkin(3) + t625) * t577 + t592;
t583 = sin(qJ(6));
t587 = cos(qJ(6));
t452 = t455 * t587 - t456 * t583;
t531 = -t551 * t583 - t578 * t587;
t477 = qJD(6) * t531 + t505 * t587 - t577 * t583;
t532 = t551 * t587 - t578 * t583;
t485 = -mrSges(7,1) * t531 + mrSges(7,2) * t532;
t503 = qJDD(6) + t506;
t546 = qJD(6) + t552;
t508 = -mrSges(7,2) * t546 + mrSges(7,3) * t531;
t449 = m(7) * t452 + mrSges(7,1) * t503 - mrSges(7,3) * t477 - t485 * t532 + t508 * t546;
t453 = t455 * t583 + t456 * t587;
t476 = -qJD(6) * t532 - t505 * t583 - t577 * t587;
t509 = mrSges(7,1) * t546 - mrSges(7,3) * t532;
t450 = m(7) * t453 - mrSges(7,2) * t503 + mrSges(7,3) * t476 + t485 * t531 - t509 * t546;
t437 = t449 * t587 + t450 * t583;
t459 = -pkin(4) * t505 - t552 * t623 + t593;
t540 = mrSges(6,2) * t578 - mrSges(6,3) * t552;
t603 = m(6) * t459 + mrSges(6,1) * t506 + t540 * t552 + t437;
t630 = m(4) * t507 + (mrSges(4,1) - mrSges(6,2)) * t505 + (mrSges(4,2) - mrSges(5,3)) * t506 - (-t536 + t535) * t551 + (t538 - t539) * t552 - t603 + t631;
t620 = -mrSges(4,3) - mrSges(5,2);
t522 = mrSges(6,1) * t552 + mrSges(6,2) * t551;
t468 = -t577 * pkin(3) + t601;
t438 = -t449 * t583 + t450 * t587;
t461 = (-pkin(3) - pkin(4)) * t577 + t592;
t604 = m(6) * t461 + mrSges(6,2) * t577 + t535 * t578 + t438;
t598 = m(5) * t468 - mrSges(5,1) * t577 - t541 * t578 + t604;
t524 = mrSges(5,1) * t551 - mrSges(5,3) * t552;
t615 = -mrSges(4,1) * t551 - mrSges(4,2) * t552 - t524;
t431 = m(4) * t470 + mrSges(4,1) * t577 + t536 * t578 + (t522 + t615) * t552 + (mrSges(6,3) + t620) * t506 - t598;
t471 = t489 * t584 + t490 * t624;
t565 = t578 * t627;
t602 = pkin(3) * t629 - qJ(4) * t577 + t523 * t551 - t471;
t467 = t565 - t602;
t597 = pkin(4) * t547 - qJ(5) * t505 + t602;
t458 = pkin(5) * t577 - pkin(9) * t629 + t537 * t578 + t565 + ((2 * qJD(5)) + t526) * t551 - t597;
t454 = -m(7) * t458 + mrSges(7,1) * t476 - mrSges(7,2) * t477 + t508 * t531 - t509 * t532;
t462 = t551 * t626 + (t628 - t537) * t578 + t597;
t595 = -m(6) * t462 + mrSges(6,3) * t505 + t522 * t551 - t454;
t594 = m(5) * t467 + mrSges(5,3) * t577 + t539 * t578 + t595;
t441 = t594 + (-t538 + t540) * t578 + (-mrSges(4,2) + mrSges(6,1)) * t577 + t615 * t551 + t620 * t505 + m(4) * t471;
t429 = t431 * t624 + t441 * t584;
t607 = t551 * t634 - t552 * t635 + t578 * t633;
t606 = -t431 * t584 + t441 * t624;
t479 = Ifges(7,4) * t532 + Ifges(7,2) * t531 + Ifges(7,6) * t546;
t480 = Ifges(7,1) * t532 + Ifges(7,4) * t531 + Ifges(7,5) * t546;
t599 = mrSges(7,1) * t452 - mrSges(7,2) * t453 + Ifges(7,5) * t477 + Ifges(7,6) * t476 + Ifges(7,3) * t503 + t479 * t532 - t531 * t480;
t435 = t505 * mrSges(6,2) + t551 * t535 + t603;
t434 = (-t522 + t524) * t552 + (mrSges(5,2) - mrSges(6,3)) * t506 + t598;
t436 = -mrSges(6,3) * t506 - t522 * t552 + t604;
t478 = Ifges(7,5) * t532 + Ifges(7,6) * t531 + Ifges(7,3) * t546;
t442 = -mrSges(7,1) * t458 + mrSges(7,3) * t453 + Ifges(7,4) * t477 + Ifges(7,2) * t476 + Ifges(7,6) * t503 - t478 * t532 + t480 * t546;
t443 = mrSges(7,2) * t458 - mrSges(7,3) * t452 + Ifges(7,1) * t477 + Ifges(7,4) * t476 + Ifges(7,5) * t503 + t478 * t531 - t479 * t546;
t591 = qJ(4) * (t540 * t578 + t594) - t583 * t443 - t587 * t442 - mrSges(5,1) * t468 + mrSges(4,1) * t470 - mrSges(4,2) * t471 - mrSges(6,1) * t462 + mrSges(5,3) * t467 + mrSges(6,2) * t461 - pkin(5) * t454 - pkin(9) * t438 - pkin(4) * t436 - pkin(3) * t434 + (-qJ(4) * t524 + t636) * t551 + (mrSges(6,1) * qJ(4) - t633) * t577 - t632 * t552 + t635 * t506 + (-mrSges(5,2) * qJ(4) - t634) * t505;
t562 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t610;
t561 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t611;
t558 = (-mrSges(3,1) * t588 + mrSges(3,2) * t585) * qJD(1);
t550 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t585 + Ifges(3,4) * t588) * qJD(1);
t549 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t585 + Ifges(3,2) * t588) * qJD(1);
t533 = -t588 * g(3) - t617;
t432 = -t506 * mrSges(5,3) - t552 * t539 - t435 + t631;
t428 = mrSges(6,1) * t459 + mrSges(4,2) * t507 + mrSges(5,2) * t468 - mrSges(4,3) * t470 - mrSges(5,3) * t464 - mrSges(6,3) * t461 + pkin(5) * t437 - qJ(4) * t432 - qJ(5) * t436 + t505 * t641 + t506 * t642 + t551 * t607 + t577 * t635 + t578 * t632 + t599;
t427 = t583 * t442 - t587 * t443 - qJ(5) * t595 - mrSges(4,1) * t507 + mrSges(4,3) * t471 - mrSges(5,1) * t464 + mrSges(5,2) * t467 - mrSges(6,2) * t459 + mrSges(6,3) * t462 + pkin(9) * t437 + pkin(4) * t435 - pkin(3) * t432 + (-qJ(5) * t540 + t636) * t578 + (-mrSges(6,1) * qJ(5) + t634) * t577 + t607 * t552 - t641 * t506 - t640 * t505;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t612 - mrSges(2,2) * t605 + t585 * (mrSges(3,2) * t553 - mrSges(3,3) * t533 + Ifges(3,1) * t559 + Ifges(3,4) * t560 + Ifges(3,5) * qJDD(2) - pkin(8) * t429 - qJD(2) * t549 - t584 * t427 + t428 * t624) + t588 * (-mrSges(3,1) * t553 + mrSges(3,3) * t534 + Ifges(3,4) * t559 + Ifges(3,2) * t560 + Ifges(3,6) * qJDD(2) - pkin(2) * t630 + pkin(8) * t606 + qJD(2) * t550 + t624 * t427 + t584 * t428) + pkin(1) * ((-t561 * t585 + t562 * t588) * qJD(1) - t559 * mrSges(3,2) + t560 * mrSges(3,1) - m(3) * t553 - t630) + pkin(7) * (t588 * (m(3) * t534 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t560 - qJD(2) * t561 + t558 * t610 + t606) - t585 * (m(3) * t533 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t559 + qJD(2) * t562 - t558 * t611 + t429)); t591 + (t585 * t549 - t588 * t550) * qJD(1) + Ifges(3,5) * t559 + Ifges(3,6) * t560 + mrSges(3,1) * t533 - mrSges(3,2) * t534 + Ifges(3,3) * qJDD(2) + pkin(2) * t429; t591; t434; t435; t599;];
tauJ  = t1;
