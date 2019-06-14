% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 01:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:58:03
% EndTime: 2019-05-05 00:58:21
% DurationCPUTime: 17.79s
% Computational Cost: add. (295883->321), mult. (659209->412), div. (0->0), fcn. (504669->14), ass. (0->146)
t688 = qJD(2) ^ 2;
t674 = sin(pkin(11));
t677 = cos(pkin(11));
t660 = g(1) * t674 - g(2) * t677;
t661 = -g(1) * t677 - g(2) * t674;
t672 = -g(3) + qJDD(1);
t675 = sin(pkin(6));
t678 = cos(pkin(6));
t682 = sin(qJ(2));
t686 = cos(qJ(2));
t632 = -t682 * t661 + (t660 * t678 + t672 * t675) * t686;
t676 = cos(pkin(12));
t719 = pkin(3) * t676;
t673 = sin(pkin(12));
t718 = mrSges(4,2) * t673;
t692 = qJDD(3) - t632;
t624 = -qJDD(2) * pkin(2) - t688 * qJ(3) + t692;
t670 = t673 ^ 2;
t714 = t678 * t682;
t715 = t675 * t682;
t633 = t660 * t714 + t686 * t661 + t672 * t715;
t628 = -pkin(2) * t688 + qJDD(2) * qJ(3) + t633;
t649 = -t660 * t675 + t672 * t678;
t708 = qJD(2) * qJD(3);
t712 = t676 * t649 - 0.2e1 * t673 * t708;
t605 = (-pkin(8) * qJDD(2) + t688 * t719 - t628) * t673 + t712;
t608 = t673 * t649 + (t628 + 0.2e1 * t708) * t676;
t706 = qJDD(2) * t676;
t671 = t676 ^ 2;
t716 = t671 * t688;
t606 = -pkin(3) * t716 + pkin(8) * t706 + t608;
t681 = sin(qJ(4));
t685 = cos(qJ(4));
t593 = t681 * t605 + t685 * t606;
t710 = qJD(2) * t676;
t711 = qJD(2) * t673;
t653 = -t681 * t711 + t685 * t710;
t696 = t673 * t685 + t676 * t681;
t654 = t696 * qJD(2);
t639 = -pkin(4) * t653 - pkin(9) * t654;
t687 = qJD(4) ^ 2;
t588 = -pkin(4) * t687 + qJDD(4) * pkin(9) + t639 * t653 + t593;
t618 = (-pkin(2) - t719) * qJDD(2) + (-qJ(3) + (-t670 - t671) * pkin(8)) * t688 + t692;
t651 = t654 * qJD(4);
t707 = qJDD(2) * t673;
t640 = -t681 * t707 + t685 * t706 - t651;
t709 = t653 * qJD(4);
t641 = qJDD(2) * t696 + t709;
t596 = (-t641 - t709) * pkin(9) + (-t640 + t651) * pkin(4) + t618;
t680 = sin(qJ(5));
t684 = cos(qJ(5));
t583 = -t680 * t588 + t684 * t596;
t643 = qJD(4) * t684 - t654 * t680;
t617 = qJD(5) * t643 + qJDD(4) * t680 + t641 * t684;
t638 = qJDD(5) - t640;
t644 = qJD(4) * t680 + t654 * t684;
t652 = qJD(5) - t653;
t581 = (t643 * t652 - t617) * pkin(10) + (t643 * t644 + t638) * pkin(5) + t583;
t584 = t684 * t588 + t680 * t596;
t616 = -qJD(5) * t644 + qJDD(4) * t684 - t641 * t680;
t627 = pkin(5) * t652 - pkin(10) * t644;
t642 = t643 ^ 2;
t582 = -pkin(5) * t642 + pkin(10) * t616 - t627 * t652 + t584;
t679 = sin(qJ(6));
t683 = cos(qJ(6));
t579 = t581 * t683 - t582 * t679;
t619 = t643 * t683 - t644 * t679;
t591 = qJD(6) * t619 + t616 * t679 + t617 * t683;
t620 = t643 * t679 + t644 * t683;
t601 = -mrSges(7,1) * t619 + mrSges(7,2) * t620;
t650 = qJD(6) + t652;
t609 = -mrSges(7,2) * t650 + mrSges(7,3) * t619;
t635 = qJDD(6) + t638;
t577 = m(7) * t579 + mrSges(7,1) * t635 - mrSges(7,3) * t591 - t601 * t620 + t609 * t650;
t580 = t581 * t679 + t582 * t683;
t590 = -qJD(6) * t620 + t616 * t683 - t617 * t679;
t610 = mrSges(7,1) * t650 - mrSges(7,3) * t620;
t578 = m(7) * t580 - mrSges(7,2) * t635 + mrSges(7,3) * t590 + t601 * t619 - t610 * t650;
t569 = t683 * t577 + t679 * t578;
t621 = -mrSges(6,1) * t643 + mrSges(6,2) * t644;
t625 = -mrSges(6,2) * t652 + mrSges(6,3) * t643;
t567 = m(6) * t583 + mrSges(6,1) * t638 - mrSges(6,3) * t617 - t621 * t644 + t625 * t652 + t569;
t626 = mrSges(6,1) * t652 - mrSges(6,3) * t644;
t701 = -t577 * t679 + t683 * t578;
t568 = m(6) * t584 - mrSges(6,2) * t638 + mrSges(6,3) * t616 + t621 * t643 - t626 * t652 + t701;
t563 = t684 * t567 + t680 * t568;
t647 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t653;
t648 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t654;
t691 = m(5) * t618 - t640 * mrSges(5,1) + mrSges(5,2) * t641 - t653 * t647 + t648 * t654 + t563;
t690 = -m(4) * t624 + mrSges(4,1) * t706 - t691 + (t670 * t688 + t716) * mrSges(4,3);
t559 = (mrSges(3,1) - t718) * qJDD(2) + m(3) * t632 - mrSges(3,2) * t688 + t690;
t717 = t559 * t686;
t636 = -mrSges(5,1) * t653 + mrSges(5,2) * t654;
t702 = -t567 * t680 + t684 * t568;
t562 = m(5) * t593 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t640 - qJD(4) * t648 + t636 * t653 + t702;
t592 = t605 * t685 - t681 * t606;
t587 = -qJDD(4) * pkin(4) - pkin(9) * t687 + t654 * t639 - t592;
t585 = -pkin(5) * t616 - pkin(10) * t642 + t627 * t644 + t587;
t693 = m(7) * t585 - t590 * mrSges(7,1) + mrSges(7,2) * t591 - t619 * t609 + t610 * t620;
t689 = -m(6) * t587 + t616 * mrSges(6,1) - mrSges(6,2) * t617 + t643 * t625 - t626 * t644 - t693;
t573 = m(5) * t592 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t641 + qJD(4) * t647 - t636 * t654 + t689;
t554 = t681 * t562 + t685 * t573;
t607 = -t628 * t673 + t712;
t695 = mrSges(4,3) * qJDD(2) + t688 * (-mrSges(4,1) * t676 + t718);
t552 = m(4) * t607 - t673 * t695 + t554;
t703 = t685 * t562 - t681 * t573;
t553 = m(4) * t608 + t676 * t695 + t703;
t704 = -t552 * t673 + t676 * t553;
t544 = m(3) * t633 - mrSges(3,1) * t688 - qJDD(2) * mrSges(3,2) + t704;
t547 = t676 * t552 + t673 * t553;
t546 = m(3) * t649 + t547;
t535 = t544 * t714 - t546 * t675 + t678 * t717;
t533 = m(2) * t660 + t535;
t539 = t686 * t544 - t559 * t682;
t538 = m(2) * t661 + t539;
t713 = t677 * t533 + t674 * t538;
t534 = t544 * t715 + t678 * t546 + t675 * t717;
t705 = -t533 * t674 + t677 * t538;
t700 = Ifges(4,1) * t673 + Ifges(4,4) * t676;
t699 = Ifges(4,4) * t673 + Ifges(4,2) * t676;
t698 = Ifges(4,5) * t673 + Ifges(4,6) * t676;
t597 = Ifges(7,5) * t620 + Ifges(7,6) * t619 + Ifges(7,3) * t650;
t599 = Ifges(7,1) * t620 + Ifges(7,4) * t619 + Ifges(7,5) * t650;
t570 = -mrSges(7,1) * t585 + mrSges(7,3) * t580 + Ifges(7,4) * t591 + Ifges(7,2) * t590 + Ifges(7,6) * t635 - t597 * t620 + t599 * t650;
t598 = Ifges(7,4) * t620 + Ifges(7,2) * t619 + Ifges(7,6) * t650;
t571 = mrSges(7,2) * t585 - mrSges(7,3) * t579 + Ifges(7,1) * t591 + Ifges(7,4) * t590 + Ifges(7,5) * t635 + t597 * t619 - t598 * t650;
t611 = Ifges(6,5) * t644 + Ifges(6,6) * t643 + Ifges(6,3) * t652;
t613 = Ifges(6,1) * t644 + Ifges(6,4) * t643 + Ifges(6,5) * t652;
t555 = -mrSges(6,1) * t587 + mrSges(6,3) * t584 + Ifges(6,4) * t617 + Ifges(6,2) * t616 + Ifges(6,6) * t638 - pkin(5) * t693 + pkin(10) * t701 + t683 * t570 + t679 * t571 - t644 * t611 + t652 * t613;
t612 = Ifges(6,4) * t644 + Ifges(6,2) * t643 + Ifges(6,6) * t652;
t556 = mrSges(6,2) * t587 - mrSges(6,3) * t583 + Ifges(6,1) * t617 + Ifges(6,4) * t616 + Ifges(6,5) * t638 - pkin(10) * t569 - t570 * t679 + t571 * t683 + t611 * t643 - t612 * t652;
t629 = Ifges(5,5) * t654 + Ifges(5,6) * t653 + Ifges(5,3) * qJD(4);
t630 = Ifges(5,4) * t654 + Ifges(5,2) * t653 + Ifges(5,6) * qJD(4);
t540 = mrSges(5,2) * t618 - mrSges(5,3) * t592 + Ifges(5,1) * t641 + Ifges(5,4) * t640 + Ifges(5,5) * qJDD(4) - pkin(9) * t563 - qJD(4) * t630 - t555 * t680 + t556 * t684 + t629 * t653;
t631 = Ifges(5,1) * t654 + Ifges(5,4) * t653 + Ifges(5,5) * qJD(4);
t548 = Ifges(5,4) * t641 + Ifges(5,2) * t640 + Ifges(5,6) * qJDD(4) - t654 * t629 + qJD(4) * t631 - mrSges(5,1) * t618 + mrSges(5,3) * t593 - Ifges(6,5) * t617 - Ifges(6,6) * t616 - Ifges(6,3) * t638 - t644 * t612 + t643 * t613 - mrSges(6,1) * t583 + mrSges(6,2) * t584 - Ifges(7,5) * t591 - Ifges(7,6) * t590 - Ifges(7,3) * t635 - t620 * t598 + t619 * t599 - mrSges(7,1) * t579 + mrSges(7,2) * t580 - pkin(5) * t569 - pkin(4) * t563;
t659 = t698 * qJD(2);
t529 = -mrSges(4,1) * t624 + mrSges(4,3) * t608 - pkin(3) * t691 + pkin(8) * t703 + qJDD(2) * t699 + t681 * t540 + t685 * t548 - t659 * t711;
t531 = mrSges(4,2) * t624 - mrSges(4,3) * t607 - pkin(8) * t554 + qJDD(2) * t700 + t685 * t540 - t681 * t548 + t659 * t710;
t528 = mrSges(3,2) * t649 - mrSges(3,3) * t632 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t688 - qJ(3) * t547 - t529 * t673 + t531 * t676;
t530 = -pkin(2) * t547 + mrSges(3,3) * t633 - mrSges(3,1) * t649 - pkin(3) * t554 - mrSges(4,1) * t607 + mrSges(4,2) * t608 - t684 * t555 - pkin(4) * t689 - pkin(9) * t702 - t680 * t556 - Ifges(5,5) * t641 - Ifges(5,6) * t640 - Ifges(5,3) * qJDD(4) - t654 * t630 + t653 * t631 - mrSges(5,1) * t592 + mrSges(5,2) * t593 + (Ifges(3,6) - t698) * qJDD(2) + (-t673 * t699 + t676 * t700 + Ifges(3,5)) * t688;
t694 = pkin(7) * t539 + t528 * t682 + t530 * t686;
t527 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t632 - mrSges(3,2) * t633 + t673 * t531 + t676 * t529 + pkin(2) * (-mrSges(4,2) * t707 + t690) + qJ(3) * t704;
t526 = mrSges(2,2) * t672 - mrSges(2,3) * t660 + t686 * t528 - t682 * t530 + (-t534 * t675 - t535 * t678) * pkin(7);
t525 = -mrSges(2,1) * t672 + mrSges(2,3) * t661 - pkin(1) * t534 - t675 * t527 + t678 * t694;
t1 = [-m(1) * g(1) + t705; -m(1) * g(2) + t713; -m(1) * g(3) + m(2) * t672 + t534; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t713 - t674 * t525 + t677 * t526; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t705 + t677 * t525 + t674 * t526; -mrSges(1,1) * g(2) + mrSges(2,1) * t660 + mrSges(1,2) * g(1) - mrSges(2,2) * t661 + pkin(1) * t535 + t527 * t678 + t675 * t694;];
tauB  = t1;
