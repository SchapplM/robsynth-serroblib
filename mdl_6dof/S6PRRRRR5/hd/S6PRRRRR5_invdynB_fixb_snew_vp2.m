% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:57:38
% EndTime: 2019-05-05 11:58:25
% DurationCPUTime: 46.28s
% Computational Cost: add. (833328->353), mult. (1713773->467), div. (0->0), fcn. (1388311->16), ass. (0->156)
t698 = sin(pkin(7));
t706 = sin(qJ(3));
t711 = cos(qJ(3));
t727 = qJD(2) * qJD(3);
t682 = (-qJDD(2) * t711 + t706 * t727) * t698;
t697 = sin(pkin(13));
t700 = cos(pkin(13));
t689 = g(1) * t697 - g(2) * t700;
t690 = -g(1) * t700 - g(2) * t697;
t696 = -g(3) + qJDD(1);
t707 = sin(qJ(2));
t702 = cos(pkin(6));
t712 = cos(qJ(2));
t730 = t702 * t712;
t699 = sin(pkin(6));
t733 = t699 * t712;
t656 = t689 * t730 - t690 * t707 + t696 * t733;
t713 = qJD(2) ^ 2;
t737 = pkin(9) * t698;
t652 = qJDD(2) * pkin(2) + t713 * t737 + t656;
t731 = t702 * t707;
t734 = t699 * t707;
t657 = t689 * t731 + t712 * t690 + t696 * t734;
t653 = -pkin(2) * t713 + qJDD(2) * t737 + t657;
t674 = -t689 * t699 + t696 * t702;
t701 = cos(pkin(7));
t621 = -t706 * t653 + (t652 * t701 + t674 * t698) * t711;
t695 = qJD(2) * t701 + qJD(3);
t728 = qJD(2) * t698;
t725 = t711 * t728;
t678 = -mrSges(4,2) * t695 + mrSges(4,3) * t725;
t679 = (-mrSges(4,1) * t711 + mrSges(4,2) * t706) * t728;
t681 = (qJDD(2) * t706 + t711 * t727) * t698;
t694 = qJDD(2) * t701 + qJDD(3);
t732 = t701 * t706;
t735 = t698 * t706;
t622 = t652 * t732 + t711 * t653 + t674 * t735;
t680 = (-pkin(3) * t711 - pkin(10) * t706) * t728;
t693 = t695 ^ 2;
t617 = -pkin(3) * t693 + pkin(10) * t694 + t680 * t725 + t622;
t668 = t701 * t674;
t620 = t682 * pkin(3) - t681 * pkin(10) + t668 + (-t652 + (pkin(3) * t706 - pkin(10) * t711) * t695 * qJD(2)) * t698;
t705 = sin(qJ(4));
t710 = cos(qJ(4));
t609 = t710 * t617 + t705 * t620;
t726 = t706 * t728;
t672 = t695 * t710 - t705 * t726;
t673 = t695 * t705 + t710 * t726;
t655 = -pkin(4) * t672 - pkin(11) * t673;
t676 = qJDD(4) + t682;
t688 = qJD(4) - t725;
t686 = t688 ^ 2;
t601 = -pkin(4) * t686 + pkin(11) * t676 + t655 * t672 + t609;
t616 = -t694 * pkin(3) - t693 * pkin(10) + t680 * t726 - t621;
t647 = -t673 * qJD(4) - t705 * t681 + t694 * t710;
t648 = qJD(4) * t672 + t681 * t710 + t694 * t705;
t604 = (-t672 * t688 - t648) * pkin(11) + (t673 * t688 - t647) * pkin(4) + t616;
t704 = sin(qJ(5));
t709 = cos(qJ(5));
t596 = -t704 * t601 + t709 * t604;
t659 = -t673 * t704 + t688 * t709;
t625 = qJD(5) * t659 + t648 * t709 + t676 * t704;
t645 = qJDD(5) - t647;
t660 = t673 * t709 + t688 * t704;
t671 = qJD(5) - t672;
t594 = (t659 * t671 - t625) * pkin(12) + (t659 * t660 + t645) * pkin(5) + t596;
t597 = t709 * t601 + t704 * t604;
t624 = -qJD(5) * t660 - t648 * t704 + t676 * t709;
t639 = pkin(5) * t671 - pkin(12) * t660;
t658 = t659 ^ 2;
t595 = -pkin(5) * t658 + pkin(12) * t624 - t639 * t671 + t597;
t703 = sin(qJ(6));
t708 = cos(qJ(6));
t592 = t594 * t708 - t595 * t703;
t632 = t659 * t708 - t660 * t703;
t607 = qJD(6) * t632 + t624 * t703 + t625 * t708;
t633 = t659 * t703 + t660 * t708;
t618 = -mrSges(7,1) * t632 + mrSges(7,2) * t633;
t669 = qJD(6) + t671;
t626 = -mrSges(7,2) * t669 + mrSges(7,3) * t632;
t640 = qJDD(6) + t645;
t590 = m(7) * t592 + mrSges(7,1) * t640 - mrSges(7,3) * t607 - t618 * t633 + t626 * t669;
t593 = t594 * t703 + t595 * t708;
t606 = -qJD(6) * t633 + t624 * t708 - t625 * t703;
t627 = mrSges(7,1) * t669 - mrSges(7,3) * t633;
t591 = m(7) * t593 - mrSges(7,2) * t640 + mrSges(7,3) * t606 + t618 * t632 - t627 * t669;
t582 = t708 * t590 + t703 * t591;
t634 = -mrSges(6,1) * t659 + mrSges(6,2) * t660;
t637 = -mrSges(6,2) * t671 + mrSges(6,3) * t659;
t580 = m(6) * t596 + mrSges(6,1) * t645 - mrSges(6,3) * t625 - t634 * t660 + t637 * t671 + t582;
t638 = mrSges(6,1) * t671 - mrSges(6,3) * t660;
t721 = -t590 * t703 + t708 * t591;
t581 = m(6) * t597 - mrSges(6,2) * t645 + mrSges(6,3) * t624 + t634 * t659 - t638 * t671 + t721;
t578 = t580 * t709 + t581 * t704;
t661 = -mrSges(5,2) * t688 + mrSges(5,3) * t672;
t662 = mrSges(5,1) * t688 - mrSges(5,3) * t673;
t715 = -m(5) * t616 + t647 * mrSges(5,1) - mrSges(5,2) * t648 + t672 * t661 - t662 * t673 - t578;
t574 = m(4) * t621 + mrSges(4,1) * t694 - mrSges(4,3) * t681 + t678 * t695 - t679 * t726 + t715;
t736 = t574 * t711;
t677 = mrSges(4,1) * t695 - mrSges(4,3) * t726;
t654 = -mrSges(5,1) * t672 + mrSges(5,2) * t673;
t722 = -t580 * t704 + t709 * t581;
t577 = m(5) * t609 - mrSges(5,2) * t676 + mrSges(5,3) * t647 + t654 * t672 - t662 * t688 + t722;
t608 = -t705 * t617 + t620 * t710;
t600 = -pkin(4) * t676 - pkin(11) * t686 + t673 * t655 - t608;
t598 = -pkin(5) * t624 - pkin(12) * t658 + t639 * t660 + t600;
t716 = m(7) * t598 - t606 * mrSges(7,1) + mrSges(7,2) * t607 - t632 * t626 + t627 * t633;
t714 = -m(6) * t600 + t624 * mrSges(6,1) - mrSges(6,2) * t625 + t659 * t637 - t638 * t660 - t716;
t586 = m(5) * t608 + mrSges(5,1) * t676 - mrSges(5,3) * t648 - t654 * t673 + t661 * t688 + t714;
t723 = t710 * t577 - t586 * t705;
t566 = m(4) * t622 - mrSges(4,2) * t694 - mrSges(4,3) * t682 - t677 * t695 + t679 * t725 + t723;
t569 = t705 * t577 + t710 * t586;
t635 = -t698 * t652 + t668;
t568 = m(4) * t635 + t682 * mrSges(4,1) + t681 * mrSges(4,2) + (t677 * t706 - t678 * t711) * t728 + t569;
t555 = t566 * t732 - t568 * t698 + t701 * t736;
t551 = m(3) * t656 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t713 + t555;
t554 = t566 * t735 + t701 * t568 + t698 * t736;
t553 = m(3) * t674 + t554;
t561 = t711 * t566 - t574 * t706;
t560 = m(3) * t657 - mrSges(3,1) * t713 - qJDD(2) * mrSges(3,2) + t561;
t541 = t551 * t730 - t553 * t699 + t560 * t731;
t539 = m(2) * t689 + t541;
t547 = -t551 * t707 + t712 * t560;
t546 = m(2) * t690 + t547;
t729 = t700 * t539 + t697 * t546;
t540 = t551 * t733 + t702 * t553 + t560 * t734;
t724 = -t539 * t697 + t700 * t546;
t610 = Ifges(7,5) * t633 + Ifges(7,6) * t632 + Ifges(7,3) * t669;
t612 = Ifges(7,1) * t633 + Ifges(7,4) * t632 + Ifges(7,5) * t669;
t583 = -mrSges(7,1) * t598 + mrSges(7,3) * t593 + Ifges(7,4) * t607 + Ifges(7,2) * t606 + Ifges(7,6) * t640 - t610 * t633 + t612 * t669;
t611 = Ifges(7,4) * t633 + Ifges(7,2) * t632 + Ifges(7,6) * t669;
t584 = mrSges(7,2) * t598 - mrSges(7,3) * t592 + Ifges(7,1) * t607 + Ifges(7,4) * t606 + Ifges(7,5) * t640 + t610 * t632 - t611 * t669;
t628 = Ifges(6,5) * t660 + Ifges(6,6) * t659 + Ifges(6,3) * t671;
t630 = Ifges(6,1) * t660 + Ifges(6,4) * t659 + Ifges(6,5) * t671;
t570 = -mrSges(6,1) * t600 + mrSges(6,3) * t597 + Ifges(6,4) * t625 + Ifges(6,2) * t624 + Ifges(6,6) * t645 - pkin(5) * t716 + pkin(12) * t721 + t708 * t583 + t703 * t584 - t660 * t628 + t671 * t630;
t629 = Ifges(6,4) * t660 + Ifges(6,2) * t659 + Ifges(6,6) * t671;
t571 = mrSges(6,2) * t600 - mrSges(6,3) * t596 + Ifges(6,1) * t625 + Ifges(6,4) * t624 + Ifges(6,5) * t645 - pkin(12) * t582 - t583 * t703 + t584 * t708 + t628 * t659 - t629 * t671;
t641 = Ifges(5,5) * t673 + Ifges(5,6) * t672 + Ifges(5,3) * t688;
t642 = Ifges(5,4) * t673 + Ifges(5,2) * t672 + Ifges(5,6) * t688;
t556 = mrSges(5,2) * t616 - mrSges(5,3) * t608 + Ifges(5,1) * t648 + Ifges(5,4) * t647 + Ifges(5,5) * t676 - pkin(11) * t578 - t570 * t704 + t571 * t709 + t641 * t672 - t642 * t688;
t643 = Ifges(5,1) * t673 + Ifges(5,4) * t672 + Ifges(5,5) * t688;
t562 = Ifges(5,4) * t648 + Ifges(5,2) * t647 + Ifges(5,6) * t676 - t673 * t641 + t688 * t643 - mrSges(5,1) * t616 + mrSges(5,3) * t609 - Ifges(6,5) * t625 - Ifges(6,6) * t624 - Ifges(6,3) * t645 - t660 * t629 + t659 * t630 - mrSges(6,1) * t596 + mrSges(6,2) * t597 - Ifges(7,5) * t607 - Ifges(7,6) * t606 - Ifges(7,3) * t640 - t633 * t611 + t632 * t612 - mrSges(7,1) * t592 + mrSges(7,2) * t593 - pkin(5) * t582 - pkin(4) * t578;
t665 = Ifges(4,6) * t695 + (Ifges(4,4) * t706 + Ifges(4,2) * t711) * t728;
t666 = Ifges(4,5) * t695 + (Ifges(4,1) * t706 + Ifges(4,4) * t711) * t728;
t542 = Ifges(4,5) * t681 - Ifges(4,6) * t682 + Ifges(4,3) * t694 + mrSges(4,1) * t621 - mrSges(4,2) * t622 + t705 * t556 + t710 * t562 + pkin(3) * t715 + pkin(10) * t723 + (t665 * t706 - t666 * t711) * t728;
t664 = Ifges(4,3) * t695 + (Ifges(4,5) * t706 + Ifges(4,6) * t711) * t728;
t543 = mrSges(4,2) * t635 - mrSges(4,3) * t621 + Ifges(4,1) * t681 - Ifges(4,4) * t682 + Ifges(4,5) * t694 - pkin(10) * t569 + t556 * t710 - t562 * t705 + t664 * t725 - t665 * t695;
t548 = Ifges(4,4) * t681 - Ifges(4,2) * t682 + Ifges(4,6) * t694 - t664 * t726 + t695 * t666 - mrSges(4,1) * t635 + mrSges(4,3) * t622 - Ifges(5,5) * t648 - Ifges(5,6) * t647 - Ifges(5,3) * t676 - t673 * t642 + t672 * t643 - mrSges(5,1) * t608 + mrSges(5,2) * t609 - t704 * t571 - t709 * t570 - pkin(4) * t714 - pkin(11) * t722 - pkin(3) * t569;
t717 = pkin(9) * t561 + t543 * t706 + t548 * t711;
t536 = -mrSges(3,1) * t674 + mrSges(3,3) * t657 + t713 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t554 - t698 * t542 + t701 * t717;
t537 = mrSges(3,2) * t674 - mrSges(3,3) * t656 + Ifges(3,5) * qJDD(2) - t713 * Ifges(3,6) + t711 * t543 - t706 * t548 + (-t554 * t698 - t555 * t701) * pkin(9);
t718 = pkin(8) * t547 + t536 * t712 + t537 * t707;
t535 = mrSges(3,1) * t656 - mrSges(3,2) * t657 + Ifges(3,3) * qJDD(2) + pkin(2) * t555 + t701 * t542 + t698 * t717;
t534 = mrSges(2,2) * t696 - mrSges(2,3) * t689 - t707 * t536 + t712 * t537 + (-t540 * t699 - t541 * t702) * pkin(8);
t533 = -mrSges(2,1) * t696 + mrSges(2,3) * t690 - pkin(1) * t540 - t699 * t535 + t702 * t718;
t1 = [-m(1) * g(1) + t724; -m(1) * g(2) + t729; -m(1) * g(3) + m(2) * t696 + t540; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t729 - t697 * t533 + t700 * t534; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t724 + t700 * t533 + t697 * t534; -mrSges(1,1) * g(2) + mrSges(2,1) * t689 + mrSges(1,2) * g(1) - mrSges(2,2) * t690 + pkin(1) * t541 + t702 * t535 + t699 * t718;];
tauB  = t1;
