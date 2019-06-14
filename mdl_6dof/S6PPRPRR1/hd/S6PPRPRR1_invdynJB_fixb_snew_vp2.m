% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPRPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:53:15
% EndTime: 2019-05-04 19:53:37
% DurationCPUTime: 15.56s
% Computational Cost: add. (280674->252), mult. (483896->332), div. (0->0), fcn. (397437->16), ass. (0->125)
t681 = sin(pkin(11));
t686 = cos(pkin(11));
t673 = -g(1) * t686 - g(2) * t681;
t680 = sin(pkin(12));
t685 = cos(pkin(12));
t672 = g(1) * t681 - g(2) * t686;
t678 = -g(3) + qJDD(1);
t683 = sin(pkin(6));
t688 = cos(pkin(6));
t702 = t672 * t688 + t678 * t683;
t644 = -t680 * t673 + t685 * t702;
t645 = t685 * t673 + t680 * t702;
t659 = -t672 * t683 + t678 * t688 + qJDD(2);
t691 = sin(qJ(3));
t687 = cos(pkin(7));
t694 = cos(qJ(3));
t715 = t687 * t694;
t682 = sin(pkin(7));
t717 = t682 * t694;
t635 = t644 * t715 - t645 * t691 + t659 * t717;
t633 = qJDD(3) * pkin(3) + t635;
t716 = t687 * t691;
t718 = t682 * t691;
t636 = t644 * t716 + t694 * t645 + t659 * t718;
t696 = qJD(3) ^ 2;
t634 = -pkin(3) * t696 + t636;
t679 = sin(pkin(13));
t684 = cos(pkin(13));
t629 = t679 * t633 + t684 * t634;
t627 = -pkin(4) * t696 + qJDD(3) * pkin(9) + t629;
t640 = -t644 * t682 + t687 * t659;
t639 = qJDD(4) + t640;
t690 = sin(qJ(5));
t693 = cos(qJ(5));
t623 = t693 * t627 + t690 * t639;
t669 = (-pkin(5) * t693 - pkin(10) * t690) * qJD(3);
t695 = qJD(5) ^ 2;
t711 = qJD(3) * t693;
t621 = -pkin(5) * t695 + qJDD(5) * pkin(10) + t669 * t711 + t623;
t628 = t684 * t633 - t679 * t634;
t626 = -qJDD(3) * pkin(4) - t696 * pkin(9) - t628;
t710 = qJD(3) * qJD(5);
t708 = t693 * t710;
t670 = qJDD(3) * t690 + t708;
t709 = t690 * t710;
t671 = qJDD(3) * t693 - t709;
t624 = (-t670 - t708) * pkin(10) + (-t671 + t709) * pkin(5) + t626;
t689 = sin(qJ(6));
t692 = cos(qJ(6));
t617 = -t621 * t689 + t624 * t692;
t712 = qJD(3) * t690;
t666 = qJD(5) * t692 - t689 * t712;
t652 = qJD(6) * t666 + qJDD(5) * t689 + t670 * t692;
t667 = qJD(5) * t689 + t692 * t712;
t653 = -mrSges(7,1) * t666 + mrSges(7,2) * t667;
t676 = qJD(6) - t711;
t657 = -mrSges(7,2) * t676 + mrSges(7,3) * t666;
t665 = qJDD(6) - t671;
t615 = m(7) * t617 + mrSges(7,1) * t665 - mrSges(7,3) * t652 - t653 * t667 + t657 * t676;
t618 = t621 * t692 + t624 * t689;
t651 = -qJD(6) * t667 + qJDD(5) * t692 - t670 * t689;
t658 = mrSges(7,1) * t676 - mrSges(7,3) * t667;
t616 = m(7) * t618 - mrSges(7,2) * t665 + mrSges(7,3) * t651 + t653 * t666 - t658 * t676;
t609 = -t615 * t689 + t692 * t616;
t668 = (-mrSges(6,1) * t693 + mrSges(6,2) * t690) * qJD(3);
t674 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t712;
t607 = m(6) * t623 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t671 - qJD(5) * t674 + t668 * t711 + t609;
t714 = t693 * t639;
t620 = -qJDD(5) * pkin(5) - t695 * pkin(10) - t714 + (qJD(3) * t669 + t627) * t690;
t619 = -m(7) * t620 + t651 * mrSges(7,1) - mrSges(7,2) * t652 + t666 * t657 - t658 * t667;
t622 = -t627 * t690 + t714;
t675 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t711;
t613 = m(6) * t622 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t670 + qJD(5) * t675 - t668 * t712 + t619;
t705 = t693 * t607 - t613 * t690;
t598 = m(5) * t629 - mrSges(5,1) * t696 - qJDD(3) * mrSges(5,2) + t705;
t608 = t615 * t692 + t616 * t689;
t699 = -m(6) * t626 + t671 * mrSges(6,1) - mrSges(6,2) * t670 - t674 * t712 + t675 * t711 - t608;
t604 = m(5) * t628 + qJDD(3) * mrSges(5,1) - mrSges(5,2) * t696 + t699;
t593 = t679 * t598 + t684 * t604;
t591 = m(4) * t635 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t696 + t593;
t706 = t684 * t598 - t604 * t679;
t592 = m(4) * t636 - mrSges(4,1) * t696 - qJDD(3) * mrSges(4,2) + t706;
t602 = t690 * t607 + t693 * t613;
t601 = m(5) * t639 + t602;
t600 = m(4) * t640 + t601;
t578 = t591 * t715 + t592 * t716 - t600 * t682;
t574 = m(3) * t644 + t578;
t584 = -t591 * t691 + t694 * t592;
t583 = m(3) * t645 + t584;
t722 = t574 * t685 + t583 * t680;
t646 = Ifges(7,5) * t667 + Ifges(7,6) * t666 + Ifges(7,3) * t676;
t648 = Ifges(7,1) * t667 + Ifges(7,4) * t666 + Ifges(7,5) * t676;
t610 = -mrSges(7,1) * t620 + mrSges(7,3) * t618 + Ifges(7,4) * t652 + Ifges(7,2) * t651 + Ifges(7,6) * t665 - t646 * t667 + t648 * t676;
t647 = Ifges(7,4) * t667 + Ifges(7,2) * t666 + Ifges(7,6) * t676;
t611 = mrSges(7,2) * t620 - mrSges(7,3) * t617 + Ifges(7,1) * t652 + Ifges(7,4) * t651 + Ifges(7,5) * t665 + t646 * t666 - t647 * t676;
t661 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t690 + Ifges(6,2) * t693) * qJD(3);
t662 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t690 + Ifges(6,4) * t693) * qJD(3);
t721 = mrSges(6,1) * t622 - mrSges(6,2) * t623 + Ifges(6,5) * t670 + Ifges(6,6) * t671 + Ifges(6,3) * qJDD(5) + pkin(5) * t619 + pkin(10) * t609 + t692 * t610 + t689 * t611 + (t661 * t690 - t662 * t693) * qJD(3);
t577 = t591 * t717 + t592 * t718 + t687 * t600;
t576 = m(3) * t659 + t577;
t564 = -t576 * t683 + t722 * t688;
t562 = m(2) * t672 + t564;
t570 = -t574 * t680 + t685 * t583;
t569 = m(2) * t673 + t570;
t713 = t686 * t562 + t681 * t569;
t563 = t688 * t576 + t722 * t683;
t707 = -t562 * t681 + t686 * t569;
t704 = m(2) * t678 + t563;
t660 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t690 + Ifges(6,6) * t693) * qJD(3);
t594 = mrSges(6,2) * t626 - mrSges(6,3) * t622 + Ifges(6,1) * t670 + Ifges(6,4) * t671 + Ifges(6,5) * qJDD(5) - pkin(10) * t608 - qJD(5) * t661 - t610 * t689 + t611 * t692 + t660 * t711;
t698 = mrSges(7,1) * t617 - mrSges(7,2) * t618 + Ifges(7,5) * t652 + Ifges(7,6) * t651 + Ifges(7,3) * t665 + t647 * t667 - t648 * t666;
t595 = -mrSges(6,1) * t626 + mrSges(6,3) * t623 + Ifges(6,4) * t670 + Ifges(6,2) * t671 + Ifges(6,6) * qJDD(5) - pkin(5) * t608 + qJD(5) * t662 - t660 * t712 - t698;
t579 = mrSges(5,2) * t639 - mrSges(5,3) * t628 + Ifges(5,5) * qJDD(3) - Ifges(5,6) * t696 - pkin(9) * t602 + t594 * t693 - t595 * t690;
t585 = -mrSges(5,1) * t639 + mrSges(5,3) * t629 + t696 * Ifges(5,5) + Ifges(5,6) * qJDD(3) - pkin(4) * t602 - t721;
t565 = -mrSges(4,1) * t640 + mrSges(4,3) * t636 + t696 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t601 + qJ(4) * t706 + t679 * t579 + t684 * t585;
t566 = mrSges(4,2) * t640 - mrSges(4,3) * t635 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t696 - qJ(4) * t593 + t579 * t684 - t585 * t679;
t701 = pkin(8) * t584 + t565 * t694 + t566 * t691;
t571 = mrSges(4,1) * t635 - mrSges(4,2) * t636 + mrSges(5,1) * t628 - mrSges(5,2) * t629 + t690 * t594 + t693 * t595 + pkin(4) * t699 + pkin(9) * t705 + pkin(3) * t593 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t559 = -mrSges(3,1) * t659 + mrSges(3,3) * t645 - pkin(2) * t577 - t682 * t571 + t687 * t701;
t560 = mrSges(3,2) * t659 - mrSges(3,3) * t644 - t691 * t565 + t694 * t566 + (-t577 * t682 - t578 * t687) * pkin(8);
t700 = qJ(2) * t570 + t559 * t685 + t560 * t680;
t558 = mrSges(3,1) * t644 - mrSges(3,2) * t645 + pkin(2) * t578 + t687 * t571 + t682 * t701;
t557 = mrSges(2,2) * t678 - mrSges(2,3) * t672 - t680 * t559 + t685 * t560 + (-t563 * t683 - t564 * t688) * qJ(2);
t556 = -mrSges(2,1) * t678 + mrSges(2,3) * t673 - pkin(1) * t563 - t683 * t558 + t688 * t700;
t1 = [-m(1) * g(1) + t707; -m(1) * g(2) + t713; -m(1) * g(3) + t704; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t713 - t681 * t556 + t686 * t557; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t707 + t686 * t556 + t681 * t557; -mrSges(1,1) * g(2) + mrSges(2,1) * t672 + mrSges(1,2) * g(1) - mrSges(2,2) * t673 + pkin(1) * t564 + t688 * t558 + t683 * t700; t704; t576; t571; t601; t721; t698;];
tauJB  = t1;
