% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:36:40
% EndTime: 2019-12-31 20:36:52
% DurationCPUTime: 12.26s
% Computational Cost: add. (187084->324), mult. (427592->426), div. (0->0), fcn. (334609->12), ass. (0->134)
t599 = cos(pkin(5));
t628 = t599 * g(3);
t597 = sin(pkin(5));
t602 = sin(qJ(2));
t627 = t597 * t602;
t606 = cos(qJ(2));
t626 = t597 * t606;
t625 = t599 * t602;
t624 = t599 * t606;
t603 = sin(qJ(1));
t607 = cos(qJ(1));
t588 = t603 * g(1) - t607 * g(2);
t608 = qJD(1) ^ 2;
t580 = t608 * t597 * pkin(7) + qJDD(1) * pkin(1) + t588;
t589 = -t607 * g(1) - t603 * g(2);
t619 = qJDD(1) * t597;
t581 = -t608 * pkin(1) + pkin(7) * t619 + t589;
t622 = t580 * t625 + t606 * t581;
t553 = -g(3) * t627 + t622;
t593 = t599 * qJD(1) + qJD(2);
t621 = qJD(1) * t597;
t618 = t602 * t621;
t578 = t593 * mrSges(3,1) - mrSges(3,3) * t618;
t583 = (-mrSges(3,1) * t606 + mrSges(3,2) * t602) * t621;
t585 = -qJD(2) * t618 + t606 * t619;
t592 = t599 * qJDD(1) + qJDD(2);
t582 = (-pkin(2) * t606 - qJ(3) * t602) * t621;
t591 = t593 ^ 2;
t620 = qJD(1) * t606;
t541 = -t591 * pkin(2) + t592 * qJ(3) + (-g(3) * t602 + t582 * t620) * t597 + t622;
t584 = (qJD(2) * t620 + qJDD(1) * t602) * t597;
t542 = -t585 * pkin(2) - t628 - t584 * qJ(3) + (-t580 + (pkin(2) * t602 - qJ(3) * t606) * t593 * qJD(1)) * t597;
t596 = sin(pkin(10));
t598 = cos(pkin(10));
t574 = t596 * t593 + t598 * t618;
t514 = -0.2e1 * qJD(3) * t574 - t596 * t541 + t598 * t542;
t562 = t598 * t584 + t596 * t592;
t573 = t598 * t593 - t596 * t618;
t617 = t597 * t620;
t511 = (-t573 * t617 - t562) * pkin(8) + (t573 * t574 - t585) * pkin(3) + t514;
t515 = 0.2e1 * qJD(3) * t573 + t598 * t541 + t596 * t542;
t561 = -t596 * t584 + t598 * t592;
t563 = -pkin(3) * t617 - t574 * pkin(8);
t572 = t573 ^ 2;
t513 = -t572 * pkin(3) + t561 * pkin(8) + t563 * t617 + t515;
t601 = sin(qJ(4));
t605 = cos(qJ(4));
t508 = t601 * t511 + t605 * t513;
t556 = t601 * t573 + t605 * t574;
t527 = -t556 * qJD(4) + t605 * t561 - t601 * t562;
t555 = t605 * t573 - t601 * t574;
t535 = -t555 * mrSges(5,1) + t556 * mrSges(5,2);
t587 = qJD(4) - t617;
t546 = t587 * mrSges(5,1) - t556 * mrSges(5,3);
t577 = qJDD(4) - t585;
t536 = -t555 * pkin(4) - t556 * pkin(9);
t586 = t587 ^ 2;
t506 = -t586 * pkin(4) + t577 * pkin(9) + t555 * t536 + t508;
t552 = -g(3) * t626 + t580 * t624 - t602 * t581;
t540 = -t592 * pkin(2) - t591 * qJ(3) + t582 * t618 + qJDD(3) - t552;
t519 = -t561 * pkin(3) - t572 * pkin(8) + t574 * t563 + t540;
t528 = t555 * qJD(4) + t601 * t561 + t605 * t562;
t509 = (-t555 * t587 - t528) * pkin(9) + (t556 * t587 - t527) * pkin(4) + t519;
t600 = sin(qJ(5));
t604 = cos(qJ(5));
t503 = -t600 * t506 + t604 * t509;
t543 = -t600 * t556 + t604 * t587;
t518 = t543 * qJD(5) + t604 * t528 + t600 * t577;
t544 = t604 * t556 + t600 * t587;
t524 = -t543 * mrSges(6,1) + t544 * mrSges(6,2);
t526 = qJDD(5) - t527;
t554 = qJD(5) - t555;
t529 = -t554 * mrSges(6,2) + t543 * mrSges(6,3);
t501 = m(6) * t503 + t526 * mrSges(6,1) - t518 * mrSges(6,3) - t544 * t524 + t554 * t529;
t504 = t604 * t506 + t600 * t509;
t517 = -t544 * qJD(5) - t600 * t528 + t604 * t577;
t530 = t554 * mrSges(6,1) - t544 * mrSges(6,3);
t502 = m(6) * t504 - t526 * mrSges(6,2) + t517 * mrSges(6,3) + t543 * t524 - t554 * t530;
t613 = -t600 * t501 + t604 * t502;
t492 = m(5) * t508 - t577 * mrSges(5,2) + t527 * mrSges(5,3) + t555 * t535 - t587 * t546 + t613;
t507 = t605 * t511 - t601 * t513;
t545 = -t587 * mrSges(5,2) + t555 * mrSges(5,3);
t505 = -t577 * pkin(4) - t586 * pkin(9) + t556 * t536 - t507;
t611 = -m(6) * t505 + t517 * mrSges(6,1) - t518 * mrSges(6,2) + t543 * t529 - t544 * t530;
t497 = m(5) * t507 + t577 * mrSges(5,1) - t528 * mrSges(5,3) - t556 * t535 + t587 * t545 + t611;
t486 = t601 * t492 + t605 * t497;
t557 = -t573 * mrSges(4,1) + t574 * mrSges(4,2);
t559 = mrSges(4,2) * t617 + t573 * mrSges(4,3);
t484 = m(4) * t514 - t585 * mrSges(4,1) - t562 * mrSges(4,3) - t574 * t557 - t559 * t617 + t486;
t560 = -mrSges(4,1) * t617 - t574 * mrSges(4,3);
t614 = t605 * t492 - t601 * t497;
t485 = m(4) * t515 + t585 * mrSges(4,2) + t561 * mrSges(4,3) + t573 * t557 + t560 * t617 + t614;
t615 = -t596 * t484 + t598 * t485;
t475 = m(3) * t553 - t592 * mrSges(3,2) + t585 * mrSges(3,3) - t593 * t578 + t583 * t617 + t615;
t478 = t598 * t484 + t596 * t485;
t567 = -t597 * t580 - t628;
t579 = -t593 * mrSges(3,2) + mrSges(3,3) * t617;
t477 = m(3) * t567 - t585 * mrSges(3,1) + t584 * mrSges(3,2) + (t578 * t602 - t579 * t606) * t621 + t478;
t493 = t604 * t501 + t600 * t502;
t610 = m(5) * t519 - t527 * mrSges(5,1) + t528 * mrSges(5,2) - t555 * t545 + t556 * t546 + t493;
t609 = -m(4) * t540 + t561 * mrSges(4,1) - t562 * mrSges(4,2) + t573 * t559 - t574 * t560 - t610;
t489 = m(3) * t552 + t592 * mrSges(3,1) - t584 * mrSges(3,3) + t593 * t579 - t583 * t618 + t609;
t465 = t475 * t625 - t597 * t477 + t489 * t624;
t463 = m(2) * t588 + qJDD(1) * mrSges(2,1) - t608 * mrSges(2,2) + t465;
t471 = t606 * t475 - t602 * t489;
t470 = m(2) * t589 - t608 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t471;
t623 = t607 * t463 + t603 * t470;
t464 = t475 * t627 + t599 * t477 + t489 * t626;
t616 = -t603 * t463 + t607 * t470;
t520 = Ifges(6,5) * t544 + Ifges(6,6) * t543 + Ifges(6,3) * t554;
t522 = Ifges(6,1) * t544 + Ifges(6,4) * t543 + Ifges(6,5) * t554;
t494 = -mrSges(6,1) * t505 + mrSges(6,3) * t504 + Ifges(6,4) * t518 + Ifges(6,2) * t517 + Ifges(6,6) * t526 - t544 * t520 + t554 * t522;
t521 = Ifges(6,4) * t544 + Ifges(6,2) * t543 + Ifges(6,6) * t554;
t495 = mrSges(6,2) * t505 - mrSges(6,3) * t503 + Ifges(6,1) * t518 + Ifges(6,4) * t517 + Ifges(6,5) * t526 + t543 * t520 - t554 * t521;
t531 = Ifges(5,5) * t556 + Ifges(5,6) * t555 + Ifges(5,3) * t587;
t532 = Ifges(5,4) * t556 + Ifges(5,2) * t555 + Ifges(5,6) * t587;
t479 = mrSges(5,2) * t519 - mrSges(5,3) * t507 + Ifges(5,1) * t528 + Ifges(5,4) * t527 + Ifges(5,5) * t577 - pkin(9) * t493 - t600 * t494 + t604 * t495 + t555 * t531 - t587 * t532;
t533 = Ifges(5,1) * t556 + Ifges(5,4) * t555 + Ifges(5,5) * t587;
t480 = -mrSges(5,1) * t519 - mrSges(6,1) * t503 + mrSges(6,2) * t504 + mrSges(5,3) * t508 + Ifges(5,4) * t528 - Ifges(6,5) * t518 + Ifges(5,2) * t527 + Ifges(5,6) * t577 - Ifges(6,6) * t517 - Ifges(6,3) * t526 - pkin(4) * t493 - t544 * t521 + t543 * t522 - t556 * t531 + t587 * t533;
t547 = Ifges(4,5) * t574 + Ifges(4,6) * t573 - Ifges(4,3) * t617;
t549 = Ifges(4,1) * t574 + Ifges(4,4) * t573 - Ifges(4,5) * t617;
t466 = -mrSges(4,1) * t540 + mrSges(4,3) * t515 + Ifges(4,4) * t562 + Ifges(4,2) * t561 - Ifges(4,6) * t585 - pkin(3) * t610 + pkin(8) * t614 + t601 * t479 + t605 * t480 - t574 * t547 - t549 * t617;
t548 = Ifges(4,4) * t574 + Ifges(4,2) * t573 - Ifges(4,6) * t617;
t467 = mrSges(4,2) * t540 - mrSges(4,3) * t514 + Ifges(4,1) * t562 + Ifges(4,4) * t561 - Ifges(4,5) * t585 - pkin(8) * t486 + t605 * t479 - t601 * t480 + t573 * t547 + t548 * t617;
t564 = Ifges(3,3) * t593 + (Ifges(3,5) * t602 + Ifges(3,6) * t606) * t621;
t565 = Ifges(3,6) * t593 + (Ifges(3,4) * t602 + Ifges(3,2) * t606) * t621;
t460 = mrSges(3,2) * t567 - mrSges(3,3) * t552 + Ifges(3,1) * t584 + Ifges(3,4) * t585 + Ifges(3,5) * t592 - qJ(3) * t478 - t596 * t466 + t598 * t467 + t564 * t617 - t593 * t565;
t566 = Ifges(3,5) * t593 + (Ifges(3,1) * t602 + Ifges(3,4) * t606) * t621;
t461 = (Ifges(3,2) + Ifges(4,3)) * t585 - pkin(4) * t611 - t604 * t494 - t600 * t495 + Ifges(3,6) * t592 + t593 * t566 + Ifges(3,4) * t584 + t573 * t549 - t574 * t548 - Ifges(5,3) * t577 - Ifges(4,6) * t561 - Ifges(4,5) * t562 - mrSges(3,1) * t567 + mrSges(3,3) * t553 + t555 * t533 - t556 * t532 - Ifges(5,6) * t527 - Ifges(5,5) * t528 - mrSges(4,1) * t514 + mrSges(4,2) * t515 + mrSges(5,2) * t508 - mrSges(5,1) * t507 - pkin(3) * t486 - pkin(2) * t478 - pkin(9) * t613 - t564 * t618;
t612 = pkin(7) * t471 + t460 * t602 + t461 * t606;
t459 = Ifges(3,5) * t584 + Ifges(3,6) * t585 + Ifges(3,3) * t592 + mrSges(3,1) * t552 - mrSges(3,2) * t553 + t596 * t467 + t598 * t466 + pkin(2) * t609 + qJ(3) * t615 + (t565 * t602 - t566 * t606) * t621;
t458 = -mrSges(2,2) * g(3) - mrSges(2,3) * t588 + Ifges(2,5) * qJDD(1) - t608 * Ifges(2,6) + t606 * t460 - t602 * t461 + (-t464 * t597 - t465 * t599) * pkin(7);
t457 = mrSges(2,1) * g(3) + mrSges(2,3) * t589 + t608 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t464 - t597 * t459 + t612 * t599;
t1 = [-m(1) * g(1) + t616; -m(1) * g(2) + t623; (-m(1) - m(2)) * g(3) + t464; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t623 - t603 * t457 + t607 * t458; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t616 + t607 * t457 + t603 * t458; -mrSges(1,1) * g(2) + mrSges(2,1) * t588 + mrSges(1,2) * g(1) - mrSges(2,2) * t589 + Ifges(2,3) * qJDD(1) + pkin(1) * t465 + t599 * t459 + t612 * t597;];
tauB = t1;
