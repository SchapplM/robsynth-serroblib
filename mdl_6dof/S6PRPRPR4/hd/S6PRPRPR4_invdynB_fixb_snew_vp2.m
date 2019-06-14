% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:46:23
% EndTime: 2019-05-04 22:46:42
% DurationCPUTime: 18.00s
% Computational Cost: add. (278474->320), mult. (640024->414), div. (0->0), fcn. (486255->14), ass. (0->144)
t684 = qJD(2) ^ 2;
t672 = sin(pkin(10));
t676 = cos(pkin(10));
t657 = g(1) * t672 - g(2) * t676;
t658 = -g(1) * t676 - g(2) * t672;
t669 = -g(3) + qJDD(1);
t673 = sin(pkin(6));
t677 = cos(pkin(6));
t680 = sin(qJ(2));
t682 = cos(qJ(2));
t630 = -t680 * t658 + (t657 * t677 + t669 * t673) * t682;
t717 = cos(qJ(4));
t675 = cos(pkin(11));
t716 = pkin(3) * t675;
t671 = sin(pkin(11));
t715 = mrSges(4,2) * t671;
t688 = qJDD(3) - t630;
t619 = -qJDD(2) * pkin(2) - t684 * qJ(3) + t688;
t667 = t671 ^ 2;
t711 = t677 * t680;
t712 = t673 * t680;
t631 = t657 * t711 + t682 * t658 + t669 * t712;
t621 = -pkin(2) * t684 + qJDD(2) * qJ(3) + t631;
t648 = -t657 * t673 + t669 * t677;
t705 = qJD(2) * qJD(3);
t709 = t675 * t648 - 0.2e1 * t671 * t705;
t603 = (-pkin(8) * qJDD(2) + t684 * t716 - t621) * t671 + t709;
t606 = t671 * t648 + (t621 + 0.2e1 * t705) * t675;
t703 = qJDD(2) * t675;
t668 = t675 ^ 2;
t713 = t668 * t684;
t604 = -pkin(3) * t713 + pkin(8) * t703 + t606;
t679 = sin(qJ(4));
t588 = t679 * t603 + t717 * t604;
t702 = t675 * t717;
t708 = qJD(2) * t671;
t650 = -qJD(2) * t702 + t679 * t708;
t691 = t671 * t717 + t675 * t679;
t651 = t691 * qJD(2);
t633 = pkin(4) * t650 - qJ(5) * t651;
t683 = qJD(4) ^ 2;
t586 = -pkin(4) * t683 + qJDD(4) * qJ(5) - t633 * t650 + t588;
t613 = (-pkin(2) - t716) * qJDD(2) + (-qJ(3) + (-t667 - t668) * pkin(8)) * t684 + t688;
t704 = qJDD(2) * t671;
t707 = qJD(4) * t651;
t637 = -qJDD(2) * t702 + t679 * t704 + t707;
t706 = t650 * qJD(4);
t638 = qJDD(2) * t691 - t706;
t591 = (-t638 + t706) * qJ(5) + (t637 + t707) * pkin(4) + t613;
t670 = sin(pkin(12));
t674 = cos(pkin(12));
t643 = qJD(4) * t670 + t651 * t674;
t581 = -0.2e1 * qJD(5) * t643 - t670 * t586 + t674 * t591;
t626 = qJDD(4) * t670 + t638 * t674;
t642 = qJD(4) * t674 - t651 * t670;
t579 = (t642 * t650 - t626) * pkin(9) + (t642 * t643 + t637) * pkin(5) + t581;
t582 = 0.2e1 * qJD(5) * t642 + t674 * t586 + t670 * t591;
t624 = pkin(5) * t650 - pkin(9) * t643;
t625 = qJDD(4) * t674 - t638 * t670;
t641 = t642 ^ 2;
t580 = -pkin(5) * t641 + pkin(9) * t625 - t624 * t650 + t582;
t678 = sin(qJ(6));
t681 = cos(qJ(6));
t577 = t579 * t681 - t580 * t678;
t614 = t642 * t681 - t643 * t678;
t594 = qJD(6) * t614 + t625 * t678 + t626 * t681;
t615 = t642 * t678 + t643 * t681;
t599 = -mrSges(7,1) * t614 + mrSges(7,2) * t615;
t649 = qJD(6) + t650;
t607 = -mrSges(7,2) * t649 + mrSges(7,3) * t614;
t636 = qJDD(6) + t637;
t575 = m(7) * t577 + mrSges(7,1) * t636 - mrSges(7,3) * t594 - t599 * t615 + t607 * t649;
t578 = t579 * t678 + t580 * t681;
t593 = -qJD(6) * t615 + t625 * t681 - t626 * t678;
t608 = mrSges(7,1) * t649 - mrSges(7,3) * t615;
t576 = m(7) * t578 - mrSges(7,2) * t636 + mrSges(7,3) * t593 + t599 * t614 - t608 * t649;
t567 = t681 * t575 + t678 * t576;
t616 = -mrSges(6,1) * t642 + mrSges(6,2) * t643;
t622 = -mrSges(6,2) * t650 + mrSges(6,3) * t642;
t565 = m(6) * t581 + mrSges(6,1) * t637 - mrSges(6,3) * t626 - t616 * t643 + t622 * t650 + t567;
t623 = mrSges(6,1) * t650 - mrSges(6,3) * t643;
t697 = -t575 * t678 + t681 * t576;
t566 = m(6) * t582 - mrSges(6,2) * t637 + mrSges(6,3) * t625 + t616 * t642 - t623 * t650 + t697;
t561 = t674 * t565 + t670 * t566;
t646 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t650;
t647 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t651;
t687 = m(5) * t613 + t637 * mrSges(5,1) + t638 * mrSges(5,2) + t650 * t646 + t651 * t647 + t561;
t686 = -m(4) * t619 + mrSges(4,1) * t703 - t687 + (t667 * t684 + t713) * mrSges(4,3);
t557 = -t684 * mrSges(3,2) + (mrSges(3,1) - t715) * qJDD(2) + t686 + m(3) * t630;
t714 = t557 * t682;
t634 = mrSges(5,1) * t650 + mrSges(5,2) * t651;
t698 = -t565 * t670 + t674 * t566;
t560 = m(5) * t588 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t637 - qJD(4) * t647 - t634 * t650 + t698;
t587 = t603 * t717 - t679 * t604;
t585 = -qJDD(4) * pkin(4) - t683 * qJ(5) + t651 * t633 + qJDD(5) - t587;
t583 = -t625 * pkin(5) - t641 * pkin(9) + t643 * t624 + t585;
t689 = m(7) * t583 - t593 * mrSges(7,1) + mrSges(7,2) * t594 - t614 * t607 + t608 * t615;
t685 = -m(6) * t585 + t625 * mrSges(6,1) - mrSges(6,2) * t626 + t642 * t622 - t623 * t643 - t689;
t571 = m(5) * t587 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t638 + qJD(4) * t646 - t634 * t651 + t685;
t552 = t679 * t560 + t717 * t571;
t605 = -t621 * t671 + t709;
t692 = mrSges(4,3) * qJDD(2) + t684 * (-mrSges(4,1) * t675 + t715);
t550 = m(4) * t605 - t671 * t692 + t552;
t699 = t717 * t560 - t679 * t571;
t551 = m(4) * t606 + t675 * t692 + t699;
t700 = -t550 * t671 + t675 * t551;
t542 = m(3) * t631 - mrSges(3,1) * t684 - qJDD(2) * mrSges(3,2) + t700;
t545 = t675 * t550 + t671 * t551;
t544 = m(3) * t648 + t545;
t533 = t542 * t711 - t544 * t673 + t677 * t714;
t531 = m(2) * t657 + t533;
t537 = t682 * t542 - t557 * t680;
t536 = m(2) * t658 + t537;
t710 = t676 * t531 + t672 * t536;
t532 = t542 * t712 + t677 * t544 + t673 * t714;
t701 = -t531 * t672 + t676 * t536;
t696 = Ifges(4,1) * t671 + Ifges(4,4) * t675;
t695 = Ifges(4,4) * t671 + Ifges(4,2) * t675;
t694 = Ifges(4,5) * t671 + Ifges(4,6) * t675;
t595 = Ifges(7,5) * t615 + Ifges(7,6) * t614 + Ifges(7,3) * t649;
t597 = Ifges(7,1) * t615 + Ifges(7,4) * t614 + Ifges(7,5) * t649;
t568 = -mrSges(7,1) * t583 + mrSges(7,3) * t578 + Ifges(7,4) * t594 + Ifges(7,2) * t593 + Ifges(7,6) * t636 - t595 * t615 + t597 * t649;
t596 = Ifges(7,4) * t615 + Ifges(7,2) * t614 + Ifges(7,6) * t649;
t569 = mrSges(7,2) * t583 - mrSges(7,3) * t577 + Ifges(7,1) * t594 + Ifges(7,4) * t593 + Ifges(7,5) * t636 + t595 * t614 - t596 * t649;
t609 = Ifges(6,5) * t643 + Ifges(6,6) * t642 + Ifges(6,3) * t650;
t611 = Ifges(6,1) * t643 + Ifges(6,4) * t642 + Ifges(6,5) * t650;
t553 = -mrSges(6,1) * t585 + mrSges(6,3) * t582 + Ifges(6,4) * t626 + Ifges(6,2) * t625 + Ifges(6,6) * t637 - pkin(5) * t689 + pkin(9) * t697 + t681 * t568 + t678 * t569 - t643 * t609 + t650 * t611;
t610 = Ifges(6,4) * t643 + Ifges(6,2) * t642 + Ifges(6,6) * t650;
t554 = mrSges(6,2) * t585 - mrSges(6,3) * t581 + Ifges(6,1) * t626 + Ifges(6,4) * t625 + Ifges(6,5) * t637 - pkin(9) * t567 - t568 * t678 + t569 * t681 + t609 * t642 - t610 * t650;
t627 = Ifges(5,5) * t651 - Ifges(5,6) * t650 + Ifges(5,3) * qJD(4);
t628 = Ifges(5,4) * t651 - Ifges(5,2) * t650 + Ifges(5,6) * qJD(4);
t538 = mrSges(5,2) * t613 - mrSges(5,3) * t587 + Ifges(5,1) * t638 - Ifges(5,4) * t637 + Ifges(5,5) * qJDD(4) - qJ(5) * t561 - qJD(4) * t628 - t553 * t670 + t554 * t674 - t627 * t650;
t629 = Ifges(5,1) * t651 - Ifges(5,4) * t650 + Ifges(5,5) * qJD(4);
t546 = Ifges(5,4) * t638 + Ifges(5,6) * qJDD(4) - t651 * t627 + qJD(4) * t629 - mrSges(5,1) * t613 + mrSges(5,3) * t588 - Ifges(6,5) * t626 - Ifges(6,6) * t625 - t643 * t610 + t642 * t611 - mrSges(6,1) * t581 + mrSges(6,2) * t582 - Ifges(7,5) * t594 - Ifges(7,6) * t593 - Ifges(7,3) * t636 - t615 * t596 + t614 * t597 - mrSges(7,1) * t577 + mrSges(7,2) * t578 - pkin(5) * t567 - pkin(4) * t561 + (-Ifges(5,2) - Ifges(6,3)) * t637;
t656 = t694 * qJD(2);
t527 = -mrSges(4,1) * t619 + mrSges(4,3) * t606 - pkin(3) * t687 + pkin(8) * t699 + qJDD(2) * t695 + t679 * t538 + t546 * t717 - t656 * t708;
t529 = t675 * qJD(2) * t656 + mrSges(4,2) * t619 - mrSges(4,3) * t605 - pkin(8) * t552 + qJDD(2) * t696 + t538 * t717 - t679 * t546;
t526 = mrSges(3,2) * t648 - mrSges(3,3) * t630 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t684 - qJ(3) * t545 - t527 * t671 + t529 * t675;
t528 = -pkin(2) * t545 + mrSges(3,3) * t631 - mrSges(3,1) * t648 - pkin(3) * t552 - mrSges(4,1) * t605 + mrSges(4,2) * t606 - t670 * t554 - t674 * t553 - pkin(4) * t685 - qJ(5) * t698 - Ifges(5,5) * t638 + Ifges(5,6) * t637 - Ifges(5,3) * qJDD(4) - t651 * t628 - t650 * t629 - mrSges(5,1) * t587 + mrSges(5,2) * t588 + (Ifges(3,6) - t694) * qJDD(2) + (-t671 * t695 + t675 * t696 + Ifges(3,5)) * t684;
t690 = pkin(7) * t537 + t526 * t680 + t528 * t682;
t525 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t630 - mrSges(3,2) * t631 + t671 * t529 + t675 * t527 + pkin(2) * (-mrSges(4,2) * t704 + t686) + qJ(3) * t700;
t524 = mrSges(2,2) * t669 - mrSges(2,3) * t657 + t682 * t526 - t680 * t528 + (-t532 * t673 - t533 * t677) * pkin(7);
t523 = -mrSges(2,1) * t669 + mrSges(2,3) * t658 - pkin(1) * t532 - t673 * t525 + t677 * t690;
t1 = [-m(1) * g(1) + t701; -m(1) * g(2) + t710; -m(1) * g(3) + m(2) * t669 + t532; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t710 - t672 * t523 + t676 * t524; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t701 + t676 * t523 + t672 * t524; -mrSges(1,1) * g(2) + mrSges(2,1) * t657 + mrSges(1,2) * g(1) - mrSges(2,2) * t658 + pkin(1) * t533 + t677 * t525 + t673 * t690;];
tauB  = t1;
