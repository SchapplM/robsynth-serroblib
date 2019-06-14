% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 08:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:44:06
% EndTime: 2019-05-05 08:44:53
% DurationCPUTime: 46.20s
% Computational Cost: add. (795050->352), mult. (1668115->468), div. (0->0), fcn. (1343797->16), ass. (0->154)
t689 = sin(pkin(7));
t697 = sin(qJ(3));
t700 = cos(qJ(3));
t716 = qJD(2) * qJD(3);
t671 = (-qJDD(2) * t700 + t697 * t716) * t689;
t688 = sin(pkin(12));
t692 = cos(pkin(12));
t679 = g(1) * t688 - g(2) * t692;
t680 = -g(1) * t692 - g(2) * t688;
t686 = -g(3) + qJDD(1);
t698 = sin(qJ(2));
t694 = cos(pkin(6));
t701 = cos(qJ(2));
t719 = t694 * t701;
t690 = sin(pkin(6));
t722 = t690 * t701;
t649 = t679 * t719 - t680 * t698 + t686 * t722;
t702 = qJD(2) ^ 2;
t726 = pkin(9) * t689;
t642 = qJDD(2) * pkin(2) + t702 * t726 + t649;
t720 = t694 * t698;
t723 = t690 * t698;
t650 = t679 * t720 + t701 * t680 + t686 * t723;
t643 = -pkin(2) * t702 + qJDD(2) * t726 + t650;
t664 = -t679 * t690 + t686 * t694;
t693 = cos(pkin(7));
t612 = -t697 * t643 + (t642 * t693 + t664 * t689) * t700;
t727 = cos(qJ(4));
t685 = qJD(2) * t693 + qJD(3);
t717 = qJD(2) * t689;
t714 = t700 * t717;
t667 = -mrSges(4,2) * t685 + mrSges(4,3) * t714;
t668 = (-mrSges(4,1) * t700 + mrSges(4,2) * t697) * t717;
t670 = (qJDD(2) * t697 + t700 * t716) * t689;
t684 = qJDD(2) * t693 + qJDD(3);
t721 = t693 * t697;
t724 = t689 * t697;
t613 = t642 * t721 + t700 * t643 + t664 * t724;
t669 = (-pkin(3) * t700 - pkin(10) * t697) * t717;
t683 = t685 ^ 2;
t608 = -pkin(3) * t683 + pkin(10) * t684 + t669 * t714 + t613;
t660 = t693 * t664;
t611 = t671 * pkin(3) - t670 * pkin(10) + t660 + (-t642 + (pkin(3) * t697 - pkin(10) * t700) * t685 * qJD(2)) * t689;
t696 = sin(qJ(4));
t597 = t608 * t727 + t696 * t611;
t715 = t697 * t717;
t662 = -t685 * t727 + t696 * t715;
t663 = t696 * t685 + t715 * t727;
t644 = pkin(4) * t662 - qJ(5) * t663;
t665 = qJDD(4) + t671;
t678 = qJD(4) - t714;
t677 = t678 ^ 2;
t592 = -pkin(4) * t677 + qJ(5) * t665 - t644 * t662 + t597;
t607 = -t684 * pkin(3) - t683 * pkin(10) + t669 * t715 - t612;
t637 = qJD(4) * t663 + t670 * t696 - t684 * t727;
t638 = -t662 * qJD(4) + t670 * t727 + t696 * t684;
t595 = (t662 * t678 - t638) * qJ(5) + (t663 * t678 + t637) * pkin(4) + t607;
t687 = sin(pkin(13));
t691 = cos(pkin(13));
t652 = t663 * t691 + t678 * t687;
t587 = -0.2e1 * qJD(5) * t652 - t687 * t592 + t691 * t595;
t624 = t638 * t691 + t665 * t687;
t651 = -t663 * t687 + t678 * t691;
t585 = (t651 * t662 - t624) * pkin(11) + (t651 * t652 + t637) * pkin(5) + t587;
t588 = 0.2e1 * qJD(5) * t651 + t691 * t592 + t687 * t595;
t623 = -t638 * t687 + t665 * t691;
t630 = pkin(5) * t662 - pkin(11) * t652;
t648 = t651 ^ 2;
t586 = -pkin(5) * t648 + pkin(11) * t623 - t630 * t662 + t588;
t695 = sin(qJ(6));
t699 = cos(qJ(6));
t583 = t585 * t699 - t586 * t695;
t620 = t651 * t699 - t652 * t695;
t600 = qJD(6) * t620 + t623 * t695 + t624 * t699;
t621 = t651 * t695 + t652 * t699;
t609 = -mrSges(7,1) * t620 + mrSges(7,2) * t621;
t661 = qJD(6) + t662;
t614 = -mrSges(7,2) * t661 + mrSges(7,3) * t620;
t635 = qJDD(6) + t637;
t581 = m(7) * t583 + mrSges(7,1) * t635 - mrSges(7,3) * t600 - t609 * t621 + t614 * t661;
t584 = t585 * t695 + t586 * t699;
t599 = -qJD(6) * t621 + t623 * t699 - t624 * t695;
t615 = mrSges(7,1) * t661 - mrSges(7,3) * t621;
t582 = m(7) * t584 - mrSges(7,2) * t635 + mrSges(7,3) * t599 + t609 * t620 - t615 * t661;
t573 = t699 * t581 + t695 * t582;
t625 = -mrSges(6,1) * t651 + mrSges(6,2) * t652;
t628 = -mrSges(6,2) * t662 + mrSges(6,3) * t651;
t571 = m(6) * t587 + mrSges(6,1) * t637 - mrSges(6,3) * t624 - t625 * t652 + t628 * t662 + t573;
t629 = mrSges(6,1) * t662 - mrSges(6,3) * t652;
t710 = -t581 * t695 + t699 * t582;
t572 = m(6) * t588 - mrSges(6,2) * t637 + mrSges(6,3) * t623 + t625 * t651 - t629 * t662 + t710;
t569 = t571 * t691 + t572 * t687;
t653 = -mrSges(5,2) * t678 - mrSges(5,3) * t662;
t654 = mrSges(5,1) * t678 - mrSges(5,3) * t663;
t704 = -m(5) * t607 - t637 * mrSges(5,1) - mrSges(5,2) * t638 - t662 * t653 - t654 * t663 - t569;
t565 = m(4) * t612 + mrSges(4,1) * t684 - mrSges(4,3) * t670 + t667 * t685 - t668 * t715 + t704;
t725 = t565 * t700;
t666 = mrSges(4,1) * t685 - mrSges(4,3) * t715;
t645 = mrSges(5,1) * t662 + mrSges(5,2) * t663;
t711 = -t571 * t687 + t691 * t572;
t568 = m(5) * t597 - mrSges(5,2) * t665 - mrSges(5,3) * t637 - t645 * t662 - t654 * t678 + t711;
t596 = -t696 * t608 + t611 * t727;
t591 = -t665 * pkin(4) - t677 * qJ(5) + t663 * t644 + qJDD(5) - t596;
t589 = -t623 * pkin(5) - t648 * pkin(11) + t652 * t630 + t591;
t705 = m(7) * t589 - t599 * mrSges(7,1) + mrSges(7,2) * t600 - t620 * t614 + t615 * t621;
t703 = -m(6) * t591 + t623 * mrSges(6,1) - mrSges(6,2) * t624 + t651 * t628 - t629 * t652 - t705;
t577 = m(5) * t596 + mrSges(5,1) * t665 - mrSges(5,3) * t638 - t645 * t663 + t653 * t678 + t703;
t712 = t568 * t727 - t577 * t696;
t557 = m(4) * t613 - mrSges(4,2) * t684 - mrSges(4,3) * t671 - t666 * t685 + t668 * t714 + t712;
t560 = t696 * t568 + t577 * t727;
t626 = -t689 * t642 + t660;
t559 = m(4) * t626 + t671 * mrSges(4,1) + t670 * mrSges(4,2) + (t666 * t697 - t667 * t700) * t717 + t560;
t546 = t557 * t721 - t559 * t689 + t693 * t725;
t542 = m(3) * t649 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t702 + t546;
t545 = t557 * t724 + t693 * t559 + t689 * t725;
t544 = m(3) * t664 + t545;
t552 = t700 * t557 - t565 * t697;
t551 = m(3) * t650 - mrSges(3,1) * t702 - qJDD(2) * mrSges(3,2) + t552;
t532 = t542 * t719 - t544 * t690 + t551 * t720;
t530 = m(2) * t679 + t532;
t538 = -t542 * t698 + t701 * t551;
t537 = m(2) * t680 + t538;
t718 = t692 * t530 + t688 * t537;
t531 = t542 * t722 + t694 * t544 + t551 * t723;
t713 = -t530 * t688 + t692 * t537;
t601 = Ifges(7,5) * t621 + Ifges(7,6) * t620 + Ifges(7,3) * t661;
t603 = Ifges(7,1) * t621 + Ifges(7,4) * t620 + Ifges(7,5) * t661;
t574 = -mrSges(7,1) * t589 + mrSges(7,3) * t584 + Ifges(7,4) * t600 + Ifges(7,2) * t599 + Ifges(7,6) * t635 - t601 * t621 + t603 * t661;
t602 = Ifges(7,4) * t621 + Ifges(7,2) * t620 + Ifges(7,6) * t661;
t575 = mrSges(7,2) * t589 - mrSges(7,3) * t583 + Ifges(7,1) * t600 + Ifges(7,4) * t599 + Ifges(7,5) * t635 + t601 * t620 - t602 * t661;
t616 = Ifges(6,5) * t652 + Ifges(6,6) * t651 + Ifges(6,3) * t662;
t618 = Ifges(6,1) * t652 + Ifges(6,4) * t651 + Ifges(6,5) * t662;
t561 = -mrSges(6,1) * t591 + mrSges(6,3) * t588 + Ifges(6,4) * t624 + Ifges(6,2) * t623 + Ifges(6,6) * t637 - pkin(5) * t705 + pkin(11) * t710 + t699 * t574 + t695 * t575 - t652 * t616 + t662 * t618;
t617 = Ifges(6,4) * t652 + Ifges(6,2) * t651 + Ifges(6,6) * t662;
t562 = mrSges(6,2) * t591 - mrSges(6,3) * t587 + Ifges(6,1) * t624 + Ifges(6,4) * t623 + Ifges(6,5) * t637 - pkin(11) * t573 - t574 * t695 + t575 * t699 + t616 * t651 - t617 * t662;
t631 = Ifges(5,5) * t663 - Ifges(5,6) * t662 + Ifges(5,3) * t678;
t632 = Ifges(5,4) * t663 - Ifges(5,2) * t662 + Ifges(5,6) * t678;
t547 = mrSges(5,2) * t607 - mrSges(5,3) * t596 + Ifges(5,1) * t638 - Ifges(5,4) * t637 + Ifges(5,5) * t665 - qJ(5) * t569 - t561 * t687 + t562 * t691 - t631 * t662 - t632 * t678;
t633 = Ifges(5,1) * t663 - Ifges(5,4) * t662 + Ifges(5,5) * t678;
t553 = Ifges(5,4) * t638 + Ifges(5,6) * t665 - t663 * t631 + t678 * t633 - mrSges(5,1) * t607 + mrSges(5,3) * t597 - Ifges(6,5) * t624 - Ifges(6,6) * t623 - t652 * t617 + t651 * t618 - mrSges(6,1) * t587 + mrSges(6,2) * t588 - Ifges(7,5) * t600 - Ifges(7,6) * t599 - Ifges(7,3) * t635 - t621 * t602 + t620 * t603 - mrSges(7,1) * t583 + mrSges(7,2) * t584 - pkin(5) * t573 - pkin(4) * t569 + (-Ifges(5,2) - Ifges(6,3)) * t637;
t657 = Ifges(4,6) * t685 + (Ifges(4,4) * t697 + Ifges(4,2) * t700) * t717;
t658 = Ifges(4,5) * t685 + (Ifges(4,1) * t697 + Ifges(4,4) * t700) * t717;
t533 = Ifges(4,5) * t670 - Ifges(4,6) * t671 + Ifges(4,3) * t684 + mrSges(4,1) * t612 - mrSges(4,2) * t613 + t696 * t547 + t727 * t553 + pkin(3) * t704 + pkin(10) * t712 + (t657 * t697 - t658 * t700) * t717;
t656 = Ifges(4,3) * t685 + (Ifges(4,5) * t697 + Ifges(4,6) * t700) * t717;
t534 = mrSges(4,2) * t626 - mrSges(4,3) * t612 + Ifges(4,1) * t670 - Ifges(4,4) * t671 + Ifges(4,5) * t684 - pkin(10) * t560 + t547 * t727 - t696 * t553 + t656 * t714 - t685 * t657;
t539 = Ifges(4,4) * t670 - Ifges(4,2) * t671 + Ifges(4,6) * t684 - t656 * t715 + t685 * t658 - mrSges(4,1) * t626 + mrSges(4,3) * t613 - Ifges(5,5) * t638 + Ifges(5,6) * t637 - Ifges(5,3) * t665 - t663 * t632 - t662 * t633 - mrSges(5,1) * t596 + mrSges(5,2) * t597 - t687 * t562 - t691 * t561 - pkin(4) * t703 - qJ(5) * t711 - pkin(3) * t560;
t706 = pkin(9) * t552 + t534 * t697 + t539 * t700;
t527 = -mrSges(3,1) * t664 + mrSges(3,3) * t650 + t702 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t545 - t689 * t533 + t693 * t706;
t528 = mrSges(3,2) * t664 - mrSges(3,3) * t649 + Ifges(3,5) * qJDD(2) - t702 * Ifges(3,6) + t700 * t534 - t697 * t539 + (-t545 * t689 - t546 * t693) * pkin(9);
t707 = pkin(8) * t538 + t527 * t701 + t528 * t698;
t526 = mrSges(3,1) * t649 - mrSges(3,2) * t650 + Ifges(3,3) * qJDD(2) + pkin(2) * t546 + t693 * t533 + t689 * t706;
t525 = mrSges(2,2) * t686 - mrSges(2,3) * t679 - t698 * t527 + t701 * t528 + (-t531 * t690 - t532 * t694) * pkin(8);
t524 = -mrSges(2,1) * t686 + mrSges(2,3) * t680 - pkin(1) * t531 - t690 * t526 + t694 * t707;
t1 = [-m(1) * g(1) + t713; -m(1) * g(2) + t718; -m(1) * g(3) + m(2) * t686 + t531; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t718 - t688 * t524 + t692 * t525; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t713 + t692 * t524 + t688 * t525; -mrSges(1,1) * g(2) + mrSges(2,1) * t679 + mrSges(1,2) * g(1) - mrSges(2,2) * t680 + pkin(1) * t532 + t694 * t526 + t690 * t707;];
tauB  = t1;
