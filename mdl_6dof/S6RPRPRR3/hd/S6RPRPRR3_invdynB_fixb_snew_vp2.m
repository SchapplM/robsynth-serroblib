% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:33:43
% EndTime: 2019-05-05 18:33:57
% DurationCPUTime: 14.53s
% Computational Cost: add. (243130->342), mult. (508215->429), div. (0->0), fcn. (345771->12), ass. (0->134)
t639 = sin(qJ(1));
t643 = cos(qJ(1));
t623 = t639 * g(1) - t643 * g(2);
t615 = qJDD(1) * pkin(1) + t623;
t624 = -t643 * g(1) - t639 * g(2);
t645 = qJD(1) ^ 2;
t618 = -t645 * pkin(1) + t624;
t633 = sin(pkin(10));
t635 = cos(pkin(10));
t591 = t633 * t615 + t635 * t618;
t581 = -t645 * pkin(2) + qJDD(1) * pkin(7) + t591;
t631 = -g(3) + qJDD(2);
t638 = sin(qJ(3));
t642 = cos(qJ(3));
t574 = t642 * t581 + t638 * t631;
t617 = (-mrSges(4,1) * t642 + mrSges(4,2) * t638) * qJD(1);
t658 = qJD(1) * qJD(3);
t628 = t638 * t658;
t620 = t642 * qJDD(1) - t628;
t660 = qJD(1) * t638;
t621 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t660;
t590 = t635 * t615 - t633 * t618;
t580 = -qJDD(1) * pkin(2) - t645 * pkin(7) - t590;
t656 = t642 * t658;
t619 = t638 * qJDD(1) + t656;
t566 = (-t619 - t656) * qJ(4) + (-t620 + t628) * pkin(3) + t580;
t616 = (-pkin(3) * t642 - qJ(4) * t638) * qJD(1);
t644 = qJD(3) ^ 2;
t659 = t642 * qJD(1);
t572 = -t644 * pkin(3) + qJDD(3) * qJ(4) + t616 * t659 + t574;
t632 = sin(pkin(11));
t634 = cos(pkin(11));
t612 = t632 * qJD(3) + t634 * t660;
t546 = -0.2e1 * qJD(4) * t612 + t634 * t566 - t632 * t572;
t597 = t632 * qJDD(3) + t634 * t619;
t611 = t634 * qJD(3) - t632 * t660;
t541 = (-t611 * t659 - t597) * pkin(8) + (t611 * t612 - t620) * pkin(4) + t546;
t547 = 0.2e1 * qJD(4) * t611 + t632 * t566 + t634 * t572;
t596 = t634 * qJDD(3) - t632 * t619;
t598 = -pkin(4) * t659 - t612 * pkin(8);
t610 = t611 ^ 2;
t545 = -t610 * pkin(4) + t596 * pkin(8) + t598 * t659 + t547;
t637 = sin(qJ(5));
t641 = cos(qJ(5));
t533 = t641 * t541 - t637 * t545;
t588 = t641 * t611 - t637 * t612;
t559 = t588 * qJD(5) + t637 * t596 + t641 * t597;
t589 = t637 * t611 + t641 * t612;
t614 = qJDD(5) - t620;
t626 = qJD(5) - t659;
t531 = (t588 * t626 - t559) * pkin(9) + (t588 * t589 + t614) * pkin(5) + t533;
t534 = t637 * t541 + t641 * t545;
t558 = -t589 * qJD(5) + t641 * t596 - t637 * t597;
t577 = t626 * pkin(5) - t589 * pkin(9);
t587 = t588 ^ 2;
t532 = -pkin(5) * t587 + pkin(9) * t558 - t577 * t626 + t534;
t636 = sin(qJ(6));
t640 = cos(qJ(6));
t529 = t531 * t640 - t532 * t636;
t568 = t640 * t588 - t636 * t589;
t543 = t568 * qJD(6) + t636 * t558 + t640 * t559;
t569 = t636 * t588 + t640 * t589;
t553 = -t568 * mrSges(7,1) + t569 * mrSges(7,2);
t625 = qJD(6) + t626;
t555 = -t625 * mrSges(7,2) + t568 * mrSges(7,3);
t608 = qJDD(6) + t614;
t525 = m(7) * t529 + mrSges(7,1) * t608 - mrSges(7,3) * t543 - t553 * t569 + t555 * t625;
t530 = t531 * t636 + t532 * t640;
t542 = -t569 * qJD(6) + t640 * t558 - t636 * t559;
t556 = t625 * mrSges(7,1) - t569 * mrSges(7,3);
t526 = m(7) * t530 - mrSges(7,2) * t608 + mrSges(7,3) * t542 + t553 * t568 - t556 * t625;
t519 = t640 * t525 + t636 * t526;
t571 = -t588 * mrSges(6,1) + t589 * mrSges(6,2);
t575 = -t626 * mrSges(6,2) + t588 * mrSges(6,3);
t517 = m(6) * t533 + mrSges(6,1) * t614 - mrSges(6,3) * t559 - t571 * t589 + t575 * t626 + t519;
t576 = t626 * mrSges(6,1) - t589 * mrSges(6,3);
t650 = -t525 * t636 + t640 * t526;
t518 = m(6) * t534 - mrSges(6,2) * t614 + mrSges(6,3) * t558 + t571 * t588 - t576 * t626 + t650;
t513 = t641 * t517 + t637 * t518;
t592 = -t611 * mrSges(5,1) + t612 * mrSges(5,2);
t594 = mrSges(5,2) * t659 + t611 * mrSges(5,3);
t511 = m(5) * t546 - mrSges(5,1) * t620 - mrSges(5,3) * t597 - t592 * t612 - t594 * t659 + t513;
t595 = -mrSges(5,1) * t659 - t612 * mrSges(5,3);
t651 = -t517 * t637 + t641 * t518;
t512 = m(5) * t547 + mrSges(5,2) * t620 + mrSges(5,3) * t596 + t592 * t611 + t595 * t659 + t651;
t652 = -t511 * t632 + t634 * t512;
t506 = m(4) * t574 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t620 - qJD(3) * t621 + t617 * t659 + t652;
t573 = -t638 * t581 + t642 * t631;
t622 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t659;
t570 = -qJDD(3) * pkin(3) - t644 * qJ(4) + t616 * t660 + qJDD(4) - t573;
t554 = -t596 * pkin(4) - t610 * pkin(8) + t612 * t598 + t570;
t536 = -t558 * pkin(5) - t587 * pkin(9) + t589 * t577 + t554;
t649 = m(7) * t536 - t542 * mrSges(7,1) + t543 * mrSges(7,2) - t568 * t555 + t569 * t556;
t648 = m(6) * t554 - t558 * mrSges(6,1) + t559 * mrSges(6,2) - t588 * t575 + t589 * t576 + t649;
t646 = -m(5) * t570 + t596 * mrSges(5,1) - t597 * mrSges(5,2) + t611 * t594 - t612 * t595 - t648;
t528 = m(4) * t573 + qJDD(3) * mrSges(4,1) - t619 * mrSges(4,3) + qJD(3) * t622 - t617 * t660 + t646;
t653 = t642 * t506 - t528 * t638;
t500 = m(3) * t591 - mrSges(3,1) * t645 - qJDD(1) * mrSges(3,2) + t653;
t507 = t511 * t634 + t512 * t632;
t647 = -m(4) * t580 + t620 * mrSges(4,1) - mrSges(4,2) * t619 - t621 * t660 + t622 * t659 - t507;
t503 = m(3) * t590 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t645 + t647;
t494 = t633 * t500 + t635 * t503;
t492 = m(2) * t623 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t645 + t494;
t654 = t635 * t500 - t503 * t633;
t493 = m(2) * t624 - mrSges(2,1) * t645 - qJDD(1) * mrSges(2,2) + t654;
t661 = t643 * t492 + t639 * t493;
t501 = t638 * t506 + t642 * t528;
t657 = m(3) * t631 + t501;
t655 = -t492 * t639 + t643 * t493;
t607 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t638 + Ifges(4,4) * t642) * qJD(1);
t606 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t638 + Ifges(4,2) * t642) * qJD(1);
t605 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t638 + Ifges(4,6) * t642) * qJD(1);
t584 = Ifges(5,1) * t612 + Ifges(5,4) * t611 - Ifges(5,5) * t659;
t583 = Ifges(5,4) * t612 + Ifges(5,2) * t611 - Ifges(5,6) * t659;
t582 = Ifges(5,5) * t612 + Ifges(5,6) * t611 - Ifges(5,3) * t659;
t564 = Ifges(6,1) * t589 + Ifges(6,4) * t588 + Ifges(6,5) * t626;
t563 = Ifges(6,4) * t589 + Ifges(6,2) * t588 + Ifges(6,6) * t626;
t562 = Ifges(6,5) * t589 + Ifges(6,6) * t588 + Ifges(6,3) * t626;
t550 = Ifges(7,1) * t569 + Ifges(7,4) * t568 + Ifges(7,5) * t625;
t549 = Ifges(7,4) * t569 + Ifges(7,2) * t568 + Ifges(7,6) * t625;
t548 = Ifges(7,5) * t569 + Ifges(7,6) * t568 + Ifges(7,3) * t625;
t521 = mrSges(7,2) * t536 - mrSges(7,3) * t529 + Ifges(7,1) * t543 + Ifges(7,4) * t542 + Ifges(7,5) * t608 + t548 * t568 - t549 * t625;
t520 = -mrSges(7,1) * t536 + mrSges(7,3) * t530 + Ifges(7,4) * t543 + Ifges(7,2) * t542 + Ifges(7,6) * t608 - t548 * t569 + t550 * t625;
t509 = mrSges(6,2) * t554 - mrSges(6,3) * t533 + Ifges(6,1) * t559 + Ifges(6,4) * t558 + Ifges(6,5) * t614 - pkin(9) * t519 - t520 * t636 + t521 * t640 + t562 * t588 - t563 * t626;
t508 = -mrSges(6,1) * t554 + mrSges(6,3) * t534 + Ifges(6,4) * t559 + Ifges(6,2) * t558 + Ifges(6,6) * t614 - pkin(5) * t649 + pkin(9) * t650 + t640 * t520 + t636 * t521 - t589 * t562 + t626 * t564;
t497 = mrSges(5,2) * t570 - mrSges(5,3) * t546 + Ifges(5,1) * t597 + Ifges(5,4) * t596 - Ifges(5,5) * t620 - pkin(8) * t513 - t508 * t637 + t509 * t641 + t582 * t611 + t583 * t659;
t496 = -mrSges(5,1) * t570 + mrSges(5,3) * t547 + Ifges(5,4) * t597 + Ifges(5,2) * t596 - Ifges(5,6) * t620 - pkin(4) * t648 + pkin(8) * t651 + t641 * t508 + t637 * t509 - t612 * t582 - t584 * t659;
t495 = Ifges(4,6) * qJDD(3) - t605 * t660 + (Ifges(4,2) + Ifges(5,3)) * t620 + t611 * t584 - t612 * t583 - Ifges(6,3) * t614 + Ifges(4,4) * t619 - Ifges(5,6) * t596 - Ifges(5,5) * t597 + qJD(3) * t607 - Ifges(7,3) * t608 + t588 * t564 - t589 * t563 + mrSges(4,3) * t574 - mrSges(4,1) * t580 + t568 * t550 - t569 * t549 - Ifges(6,6) * t558 - Ifges(6,5) * t559 - Ifges(7,6) * t542 - Ifges(7,5) * t543 - mrSges(5,1) * t546 + mrSges(5,2) * t547 - mrSges(6,1) * t533 + mrSges(6,2) * t534 - mrSges(7,1) * t529 + mrSges(7,2) * t530 - pkin(5) * t519 - pkin(4) * t513 - pkin(3) * t507;
t488 = mrSges(4,2) * t580 - mrSges(4,3) * t573 + Ifges(4,1) * t619 + Ifges(4,4) * t620 + Ifges(4,5) * qJDD(3) - qJ(4) * t507 - qJD(3) * t606 - t496 * t632 + t497 * t634 + t605 * t659;
t487 = Ifges(3,6) * qJDD(1) + t645 * Ifges(3,5) - mrSges(3,1) * t631 + mrSges(3,3) * t591 - Ifges(4,5) * t619 - Ifges(4,6) * t620 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t573 + mrSges(4,2) * t574 - t632 * t497 - t634 * t496 - pkin(3) * t646 - qJ(4) * t652 - pkin(2) * t501 + (-t606 * t638 + t607 * t642) * qJD(1);
t486 = mrSges(3,2) * t631 - mrSges(3,3) * t590 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t645 - pkin(7) * t501 + t488 * t642 - t495 * t638;
t485 = -mrSges(2,2) * g(3) - mrSges(2,3) * t623 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t645 - qJ(2) * t494 + t486 * t635 - t487 * t633;
t484 = mrSges(2,1) * g(3) + mrSges(2,3) * t624 + t645 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t657 + qJ(2) * t654 + t633 * t486 + t635 * t487;
t1 = [-m(1) * g(1) + t655; -m(1) * g(2) + t661; (-m(1) - m(2)) * g(3) + t657; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t661 - t639 * t484 + t643 * t485; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t655 + t643 * t484 + t639 * t485; pkin(1) * t494 + mrSges(2,1) * t623 - mrSges(2,2) * t624 + t638 * t488 + t642 * t495 + pkin(2) * t647 + pkin(7) * t653 + mrSges(3,1) * t590 - mrSges(3,2) * t591 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
