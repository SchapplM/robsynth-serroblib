% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 17:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR15_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR15_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:08:01
% EndTime: 2019-05-07 17:08:14
% DurationCPUTime: 13.23s
% Computational Cost: add. (200322->364), mult. (501575->474), div. (0->0), fcn. (416700->14), ass. (0->166)
t646 = Ifges(4,5) - Ifges(5,4);
t654 = -Ifges(4,2) - Ifges(5,3);
t653 = Ifges(5,2) + Ifges(4,1);
t645 = Ifges(4,6) - Ifges(5,5);
t644 = -Ifges(5,6) - Ifges(4,4);
t652 = Ifges(4,3) + Ifges(5,1);
t597 = cos(qJ(2));
t588 = sin(pkin(7));
t643 = cos(pkin(6));
t612 = qJD(1) * t643 + qJD(2);
t608 = t612 * t588;
t589 = sin(pkin(6));
t632 = qJD(1) * t589;
t642 = cos(pkin(7));
t614 = t642 * t632;
t651 = t597 * t614 + t608;
t568 = t651 * pkin(10);
t593 = sin(qJ(2));
t648 = pkin(10) * t593;
t574 = (-pkin(2) * t597 - t588 * t648) * t632;
t628 = qJD(1) * qJD(2);
t580 = (qJDD(1) * t593 + t597 * t628) * t589;
t586 = qJDD(1) * t643 + qJDD(2);
t599 = qJD(1) ^ 2;
t594 = sin(qJ(1));
t598 = cos(qJ(1));
t622 = t594 * g(1) - g(2) * t598;
t649 = pkin(9) * t589;
t577 = qJDD(1) * pkin(1) + t599 * t649 + t622;
t613 = -g(1) * t598 - g(2) * t594;
t578 = -pkin(1) * t599 + qJDD(1) * t649 + t613;
t619 = t597 * t643;
t616 = t577 * t619 - t593 * t578;
t625 = t642 * pkin(10);
t631 = qJD(1) * t593;
t508 = -t580 * t625 + t586 * pkin(2) + t612 * t568 + (-g(3) * t597 - t574 * t631) * t589 + t616;
t573 = pkin(2) * t612 - t614 * t648;
t581 = (qJDD(1) * t597 - t593 * t628) * t589;
t610 = t581 * t642 + t586 * t588;
t630 = qJD(1) * t597;
t620 = t593 * t643;
t633 = t577 * t620 + t597 * t578;
t509 = -t612 * t573 + (-g(3) * t593 + t574 * t630) * t589 + t610 * pkin(10) + t633;
t626 = t643 * g(3);
t514 = -t580 * t588 * pkin(10) - t626 - t581 * pkin(2) + (-t577 + (-t568 * t597 + t573 * t593) * qJD(1)) * t589;
t592 = sin(qJ(3));
t621 = t592 * t642;
t640 = t588 * t592;
t650 = cos(qJ(3));
t485 = t508 * t621 + t650 * t509 + t514 * t640;
t624 = t589 * t631;
t555 = t592 * t624 - t650 * t651;
t556 = t592 * t608 + (t593 * t650 + t597 * t621) * t632;
t537 = pkin(3) * t555 - qJ(4) * t556;
t559 = -t588 * t581 + t586 * t642 + qJDD(3);
t623 = t589 * t630;
t570 = t588 * t623 - t612 * t642 - qJD(3);
t567 = t570 ^ 2;
t479 = pkin(3) * t567 - t559 * qJ(4) + 0.2e1 * qJD(4) * t570 + t555 * t537 - t485;
t647 = mrSges(4,1) - mrSges(5,2);
t641 = t555 * t570;
t639 = t589 * t593;
t638 = t589 * t597;
t615 = t642 * t650;
t627 = t588 * t650;
t484 = t508 * t615 - t592 * t509 + t514 * t627;
t480 = -t559 * pkin(3) - t567 * qJ(4) + t556 * t537 + qJDD(4) - t484;
t534 = -t555 * qJD(3) + t580 * t650 + t610 * t592;
t473 = (t555 * t556 - t559) * pkin(11) + (t534 - t641) * pkin(4) + t480;
t533 = qJD(3) * t556 + t580 * t592 - t581 * t615 - t586 * t627;
t547 = pkin(4) * t556 + pkin(11) * t570;
t554 = t555 ^ 2;
t490 = -t508 * t588 + t642 * t514;
t603 = (-t534 - t641) * qJ(4) + t490 + (-pkin(3) * t570 - 0.2e1 * qJD(4)) * t556;
t474 = -pkin(4) * t554 - t547 * t556 + (pkin(3) + pkin(11)) * t533 + t603;
t591 = sin(qJ(5));
t596 = cos(qJ(5));
t469 = t591 * t473 + t596 * t474;
t541 = t555 * t596 + t570 * t591;
t542 = t555 * t591 - t570 * t596;
t511 = -pkin(5) * t541 - pkin(12) * t542;
t532 = qJDD(5) + t534;
t553 = qJD(5) + t556;
t552 = t553 ^ 2;
t466 = -pkin(5) * t552 + pkin(12) * t532 + t511 * t541 + t469;
t476 = -pkin(4) * t533 - pkin(11) * t554 - t570 * t547 - t479;
t497 = -qJD(5) * t542 + t533 * t596 - t559 * t591;
t498 = qJD(5) * t541 + t533 * t591 + t559 * t596;
t470 = (-t541 * t553 - t498) * pkin(12) + (t542 * t553 - t497) * pkin(5) + t476;
t590 = sin(qJ(6));
t595 = cos(qJ(6));
t462 = -t466 * t590 + t470 * t595;
t518 = -t542 * t590 + t553 * t595;
t483 = qJD(6) * t518 + t498 * t595 + t532 * t590;
t519 = t542 * t595 + t553 * t590;
t492 = -mrSges(7,1) * t518 + mrSges(7,2) * t519;
t496 = qJDD(6) - t497;
t540 = qJD(6) - t541;
t499 = -mrSges(7,2) * t540 + mrSges(7,3) * t518;
t460 = m(7) * t462 + mrSges(7,1) * t496 - mrSges(7,3) * t483 - t492 * t519 + t499 * t540;
t463 = t466 * t595 + t470 * t590;
t482 = -qJD(6) * t519 - t498 * t590 + t532 * t595;
t500 = mrSges(7,1) * t540 - mrSges(7,3) * t519;
t461 = m(7) * t463 - mrSges(7,2) * t496 + mrSges(7,3) * t482 + t492 * t518 - t500 * t540;
t451 = t595 * t460 + t590 * t461;
t637 = t645 * t555 - t646 * t556 + t652 * t570;
t636 = t654 * t555 - t644 * t556 - t645 * t570;
t635 = -t644 * t555 - t653 * t556 + t646 * t570;
t544 = mrSges(5,1) * t555 + mrSges(5,3) * t570;
t634 = mrSges(4,2) * t570 - mrSges(4,3) * t555 - t544;
t546 = -mrSges(4,1) * t570 - mrSges(4,3) * t556;
t478 = pkin(3) * t533 + t603;
t545 = mrSges(5,1) * t556 - mrSges(5,2) * t570;
t510 = -mrSges(6,1) * t541 + mrSges(6,2) * t542;
t521 = mrSges(6,1) * t553 - mrSges(6,3) * t542;
t617 = -t460 * t590 + t595 * t461;
t449 = m(6) * t469 - mrSges(6,2) * t532 + mrSges(6,3) * t497 + t510 * t541 - t521 * t553 + t617;
t468 = t473 * t596 - t474 * t591;
t520 = -mrSges(6,2) * t553 + mrSges(6,3) * t541;
t465 = -pkin(5) * t532 - pkin(12) * t552 + t511 * t542 - t468;
t606 = -m(7) * t465 + t482 * mrSges(7,1) - mrSges(7,2) * t483 + t518 * t499 - t500 * t519;
t456 = m(6) * t468 + mrSges(6,1) * t532 - mrSges(6,3) * t498 - t510 * t542 + t520 * t553 + t606;
t618 = t596 * t449 - t456 * t591;
t609 = m(5) * t478 - t534 * mrSges(5,3) - t556 * t545 + t618;
t437 = m(4) * t490 + mrSges(4,2) * t534 + t533 * t647 + t546 * t556 + t555 * t634 + t609;
t538 = mrSges(4,1) * t555 + mrSges(4,2) * t556;
t443 = t449 * t591 + t456 * t596;
t539 = -mrSges(5,2) * t555 - mrSges(5,3) * t556;
t605 = -m(5) * t480 - t534 * mrSges(5,1) - t556 * t539 - t443;
t440 = m(4) * t484 - mrSges(4,3) * t534 - t538 * t556 + t559 * t647 - t570 * t634 + t605;
t604 = -m(6) * t476 + mrSges(6,1) * t497 - t498 * mrSges(6,2) + t520 * t541 - t542 * t521 - t451;
t601 = -m(5) * t479 + t559 * mrSges(5,3) - t570 * t545 - t604;
t447 = t601 + (-t538 - t539) * t555 + (-mrSges(4,3) - mrSges(5,1)) * t533 + m(4) * t485 - mrSges(4,2) * t559 + t546 * t570;
t430 = t642 * t437 + t440 * t627 + t447 * t640;
t433 = -t592 * t440 + t650 * t447;
t431 = -t588 * t437 + t440 * t615 + t447 * t621;
t486 = Ifges(7,5) * t519 + Ifges(7,6) * t518 + Ifges(7,3) * t540;
t488 = Ifges(7,1) * t519 + Ifges(7,4) * t518 + Ifges(7,5) * t540;
t454 = -mrSges(7,1) * t465 + mrSges(7,3) * t463 + Ifges(7,4) * t483 + Ifges(7,2) * t482 + Ifges(7,6) * t496 - t486 * t519 + t488 * t540;
t487 = Ifges(7,4) * t519 + Ifges(7,2) * t518 + Ifges(7,6) * t540;
t455 = mrSges(7,2) * t465 - mrSges(7,3) * t462 + Ifges(7,1) * t483 + Ifges(7,4) * t482 + Ifges(7,5) * t496 + t486 * t518 - t487 * t540;
t502 = Ifges(6,4) * t542 + Ifges(6,2) * t541 + Ifges(6,6) * t553;
t503 = Ifges(6,1) * t542 + Ifges(6,4) * t541 + Ifges(6,5) * t553;
t602 = mrSges(6,1) * t468 - mrSges(6,2) * t469 + Ifges(6,5) * t498 + Ifges(6,6) * t497 + Ifges(6,3) * t532 + pkin(5) * t606 + pkin(12) * t617 + t595 * t454 + t590 * t455 + t542 * t502 - t541 * t503;
t600 = mrSges(7,1) * t462 - mrSges(7,2) * t463 + Ifges(7,5) * t483 + Ifges(7,6) * t482 + Ifges(7,3) * t496 + t487 * t519 - t488 * t518;
t579 = (-mrSges(3,1) * t597 + mrSges(3,2) * t593) * t632;
t576 = -mrSges(3,2) * t612 + mrSges(3,3) * t623;
t575 = mrSges(3,1) * t612 - mrSges(3,3) * t624;
t563 = -t589 * t577 - t626;
t562 = Ifges(3,5) * qJD(2) + (Ifges(3,5) * t643 + (Ifges(3,1) * t593 + Ifges(3,4) * t597) * t589) * qJD(1);
t561 = Ifges(3,6) * qJD(2) + (Ifges(3,6) * t643 + (Ifges(3,4) * t593 + Ifges(3,2) * t597) * t589) * qJD(1);
t560 = Ifges(3,3) * qJD(2) + (Ifges(3,3) * t643 + (Ifges(3,5) * t593 + Ifges(3,6) * t597) * t589) * qJD(1);
t550 = -g(3) * t639 + t633;
t549 = -g(3) * t638 + t616;
t501 = Ifges(6,5) * t542 + Ifges(6,6) * t541 + Ifges(6,3) * t553;
t442 = mrSges(5,2) * t559 - t544 * t570 - t605;
t441 = -mrSges(5,2) * t533 - t544 * t555 + t609;
t435 = -mrSges(6,1) * t476 + mrSges(6,3) * t469 + Ifges(6,4) * t498 + Ifges(6,2) * t497 + Ifges(6,6) * t532 - pkin(5) * t451 - t501 * t542 + t503 * t553 - t600;
t434 = mrSges(6,2) * t476 - mrSges(6,3) * t468 + Ifges(6,1) * t498 + Ifges(6,4) * t497 + Ifges(6,5) * t532 - pkin(12) * t451 - t454 * t590 + t455 * t595 + t501 * t541 - t502 * t553;
t432 = m(3) * t550 - t586 * mrSges(3,2) + t581 * mrSges(3,3) - t575 * t612 + t579 * t623 + t433;
t429 = m(3) * t549 + t586 * mrSges(3,1) - t580 * mrSges(3,3) + t576 * t612 - t579 * t624 + t431;
t428 = mrSges(5,1) * t480 + mrSges(4,2) * t490 - mrSges(4,3) * t484 - mrSges(5,3) * t478 + pkin(4) * t443 - qJ(4) * t441 + t644 * t533 + t653 * t534 + t637 * t555 + t646 * t559 + t636 * t570 + t602;
t427 = -mrSges(4,1) * t490 - mrSges(5,1) * t479 + mrSges(5,2) * t478 + mrSges(4,3) * t485 - pkin(3) * t441 - pkin(4) * t604 - pkin(11) * t618 - t591 * t434 - t596 * t435 + t654 * t533 - t644 * t534 + t637 * t556 + t645 * t559 + t635 * t570;
t426 = mrSges(4,1) * t484 - mrSges(4,2) * t485 + mrSges(5,2) * t480 - mrSges(5,3) * t479 + t596 * t434 - t591 * t435 - pkin(11) * t443 - pkin(3) * t442 + qJ(4) * t601 + t652 * t559 + t636 * t556 + (-qJ(4) * t539 - t635) * t555 + t646 * t534 + (-mrSges(5,1) * qJ(4) - t645) * t533;
t425 = mrSges(3,1) * t549 - mrSges(3,2) * t550 + t642 * t426 + Ifges(3,5) * t580 + Ifges(3,6) * t581 + Ifges(3,3) * t586 + pkin(2) * t431 + (t561 * t593 - t562 * t597) * t632 + (pkin(10) * t433 + t427 * t650 + t428 * t592) * t588;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t622 - mrSges(2,2) * t613 + (Ifges(3,1) * t580 + Ifges(3,4) * t581 + Ifges(3,5) * t586 + t560 * t623 - t612 * t561 + mrSges(3,2) * t563 - mrSges(3,3) * t549 + t650 * t428 - t592 * t427 + (-t588 * t430 - t431 * t642) * pkin(10)) * t639 + (-mrSges(3,1) * t563 + mrSges(3,3) * t550 + Ifges(3,4) * t580 + Ifges(3,2) * t581 + Ifges(3,6) * t586 - pkin(2) * t430 - t588 * t426 + t427 * t615 + t428 * t621 + t433 * t625 - t560 * t624 + t562 * t612) * t638 + t643 * t425 + pkin(1) * (t429 * t619 + t432 * t620 + (-m(3) * t563 + t581 * mrSges(3,1) - t580 * mrSges(3,2) + (-t575 * t593 + t576 * t597) * t632 - t430) * t589) + (-t429 * t593 + t432 * t597) * t649; t425; t426; t442; t602; t600;];
tauJ  = t1;
