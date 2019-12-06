% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynJB_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:24:15
% EndTime: 2019-12-05 18:24:20
% DurationCPUTime: 2.59s
% Computational Cost: add. (22117->281), mult. (47676->340), div. (0->0), fcn. (28769->8), ass. (0->115)
t729 = Ifges(3,1) + Ifges(4,1);
t720 = Ifges(3,4) + Ifges(4,4);
t719 = Ifges(3,5) + Ifges(4,5);
t728 = Ifges(3,2) + Ifges(4,2);
t718 = Ifges(3,6) + Ifges(4,6);
t727 = Ifges(3,3) + Ifges(4,3);
t682 = sin(qJ(2));
t686 = cos(qJ(2));
t654 = (-mrSges(4,1) * t686 + mrSges(4,2) * t682) * qJD(1);
t705 = qJD(1) * qJD(2);
t703 = t686 * t705;
t656 = t682 * qJDD(1) + t703;
t704 = qJD(1) * qJD(3);
t688 = qJD(1) ^ 2;
t713 = t682 * t688;
t683 = sin(qJ(1));
t687 = cos(qJ(1));
t667 = -t687 * g(1) - t683 * g(2);
t714 = t682 * t667;
t694 = qJ(3) * t703 - 0.2e1 * t682 * t704 - t714 + (t713 * t686 + qJDD(2)) * pkin(1);
t716 = -pkin(3) - qJ(3);
t606 = qJDD(2) * pkin(2) + t716 * t656 + (pkin(2) * t713 + pkin(3) * t705 - g(3)) * t686 + t694;
t657 = t686 * qJDD(1) - t682 * t705;
t707 = qJD(1) * t682;
t660 = qJD(2) * pkin(1) - qJ(3) * t707;
t665 = qJD(2) * pkin(2) - pkin(3) * t707;
t638 = -t682 * g(3) + t686 * t667;
t698 = t657 * qJ(3) + 0.2e1 * t686 * t704 + t638;
t715 = t686 ^ 2 * t688;
t723 = -pkin(1) - pkin(2);
t607 = t657 * pkin(3) + t723 * t715 + (-t660 - t665) * qJD(2) + t698;
t681 = sin(qJ(4));
t685 = cos(qJ(4));
t600 = t681 * t606 + t685 * t607;
t648 = (t681 * t686 + t682 * t685) * qJD(1);
t623 = -t648 * qJD(4) - t681 * t656 + t685 * t657;
t647 = (t681 * t682 - t685 * t686) * qJD(1);
t632 = t647 * mrSges(5,1) + t648 * mrSges(5,2);
t675 = qJD(2) + qJD(4);
t636 = t675 * mrSges(5,1) - t648 * mrSges(5,3);
t674 = qJDD(2) + qJDD(4);
t596 = (t647 * t648 + t674) * pkin(4) + t600;
t666 = t683 * g(1) - t687 * g(2);
t697 = t660 * t707 + qJDD(3) - t666;
t609 = t723 * t657 + t665 * t707 + t716 * t715 + t697;
t624 = -t647 * qJD(4) + t685 * t656 + t681 * t657;
t601 = (t647 * t675 - t624) * pkin(4) + t609;
t680 = sin(qJ(5));
t684 = cos(qJ(5));
t594 = -t680 * t596 + t684 * t601;
t633 = -t680 * t648 + t684 * t675;
t605 = t633 * qJD(5) + t684 * t624 + t680 * t674;
t634 = t684 * t648 + t680 * t675;
t614 = -t633 * mrSges(6,1) + t634 * mrSges(6,2);
t622 = qJDD(5) - t623;
t640 = qJD(5) + t647;
t625 = -t640 * mrSges(6,2) + t633 * mrSges(6,3);
t592 = m(6) * t594 + t622 * mrSges(6,1) - t605 * mrSges(6,3) - t634 * t614 + t640 * t625;
t595 = t684 * t596 + t680 * t601;
t604 = -t634 * qJD(5) - t680 * t624 + t684 * t674;
t626 = t640 * mrSges(6,1) - t634 * mrSges(6,3);
t593 = m(6) * t595 - t622 * mrSges(6,2) + t604 * mrSges(6,3) + t633 * t614 - t640 * t626;
t701 = -t680 * t592 + t684 * t593;
t579 = m(5) * t600 - t674 * mrSges(5,2) + t623 * mrSges(5,3) - t647 * t632 - t675 * t636 + t701;
t599 = t685 * t606 - t681 * t607;
t597 = (-t648 ^ 2 - t675 ^ 2) * pkin(4) - t599;
t635 = -t675 * mrSges(5,2) - t647 * mrSges(5,3);
t588 = m(5) * t599 - m(6) * t597 + t674 * mrSges(5,1) + t604 * mrSges(6,1) - t605 * mrSges(6,2) - t624 * mrSges(5,3) + t633 * t625 - t634 * t626 - t648 * t632 + t675 * t635;
t574 = t681 * t579 + t685 * t588;
t722 = t686 * g(3);
t617 = -t656 * qJ(3) + t694 - t722;
t706 = qJD(1) * t686;
t663 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t706;
t700 = m(4) * t617 + qJDD(2) * mrSges(4,1) + qJD(2) * t663 + t574;
t573 = -t656 * mrSges(4,3) - t654 * t707 + t700;
t618 = -pkin(1) * t715 - qJD(2) * t660 + t698;
t637 = -t714 - t722;
t610 = Ifges(6,5) * t634 + Ifges(6,6) * t633 + Ifges(6,3) * t640;
t612 = Ifges(6,1) * t634 + Ifges(6,4) * t633 + Ifges(6,5) * t640;
t585 = -mrSges(6,1) * t597 + mrSges(6,3) * t595 + Ifges(6,4) * t605 + Ifges(6,2) * t604 + Ifges(6,6) * t622 - t634 * t610 + t640 * t612;
t611 = Ifges(6,4) * t634 + Ifges(6,2) * t633 + Ifges(6,6) * t640;
t586 = mrSges(6,2) * t597 - mrSges(6,3) * t594 + Ifges(6,1) * t605 + Ifges(6,4) * t604 + Ifges(6,5) * t622 + t633 * t610 - t640 * t611;
t629 = Ifges(5,4) * t648 - Ifges(5,2) * t647 + Ifges(5,6) * t675;
t630 = Ifges(5,1) * t648 - Ifges(5,4) * t647 + Ifges(5,5) * t675;
t691 = -mrSges(5,1) * t599 + mrSges(5,2) * t600 - Ifges(5,5) * t624 - Ifges(5,6) * t623 - Ifges(5,3) * t674 - pkin(4) * t701 - t684 * t585 - t680 * t586 - t648 * t629 - t647 * t630;
t709 = t719 * qJD(2) + (t729 * t682 + t720 * t686) * qJD(1);
t710 = -t718 * qJD(2) + (-t720 * t682 - t728 * t686) * qJD(1);
t726 = mrSges(3,1) * t637 + mrSges(4,1) * t617 - mrSges(3,2) * t638 - mrSges(4,2) * t618 + pkin(1) * t573 + pkin(2) * t574 - (t710 * t682 + t709 * t686) * qJD(1) + t727 * qJDD(2) + t719 * t656 + t718 * t657 - t691;
t721 = -mrSges(3,2) - mrSges(4,2);
t712 = t684 * t592 + t680 * t593;
t711 = t727 * qJD(2) + (t719 * t682 + t718 * t686) * qJD(1);
t661 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t707;
t708 = -qJD(2) * mrSges(3,1) + mrSges(3,3) * t707 - t661;
t702 = t685 * t579 - t681 * t588;
t628 = Ifges(5,5) * t648 - Ifges(5,6) * t647 + Ifges(5,3) * t675;
t570 = mrSges(5,2) * t609 - mrSges(5,3) * t599 + Ifges(5,1) * t624 + Ifges(5,4) * t623 + Ifges(5,5) * t674 - pkin(4) * t712 - t680 * t585 + t684 * t586 - t647 * t628 - t675 * t629;
t627 = -t657 * pkin(1) - qJ(3) * t715 + t697;
t693 = -m(5) * t609 + t623 * mrSges(5,1) - t624 * mrSges(5,2) - t647 * t635 - t648 * t636 - t712;
t692 = m(4) * t627 - t693;
t576 = -t657 * mrSges(4,1) + t656 * mrSges(4,2) + (t661 * t682 - t663 * t686) * qJD(1) + t692;
t690 = mrSges(6,1) * t594 - mrSges(6,2) * t595 + Ifges(6,5) * t605 + Ifges(6,6) * t604 + Ifges(6,3) * t622 + t634 * t611 - t633 * t612;
t580 = -mrSges(5,1) * t609 + mrSges(5,3) * t600 + Ifges(5,4) * t624 + Ifges(5,2) * t623 + Ifges(5,6) * t674 - t648 * t628 + t675 * t630 - t690;
t695 = m(4) * t618 + t657 * mrSges(4,3) + t654 * t706 + t702;
t565 = mrSges(3,1) * t666 + mrSges(3,3) * t638 - mrSges(4,1) * t627 + mrSges(4,3) * t618 + t681 * t570 + t685 * t580 + pkin(2) * t693 + pkin(3) * t702 - pkin(1) * t576 + qJ(3) * t695 + t728 * t657 + t720 * t656 + (-qJ(3) * mrSges(4,2) + t718) * qJDD(2) + (-qJ(3) * t661 + t709) * qJD(2) - t711 * t707;
t567 = -mrSges(3,2) * t666 + mrSges(4,2) * t627 - mrSges(3,3) * t637 - mrSges(4,3) * t617 - pkin(3) * t574 - qJ(3) * t573 + t710 * qJD(2) + t719 * qJDD(2) + t685 * t570 - t681 * t580 + t729 * t656 + t720 * t657 + t711 * t706;
t696 = mrSges(2,1) * t666 - mrSges(2,2) * t667 + Ifges(2,3) * qJDD(1) + t686 * t565 + t682 * t567;
t664 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t706;
t655 = (-mrSges(3,1) * t686 + mrSges(3,2) * t682) * qJD(1);
t575 = qJDD(1) * mrSges(2,1) - t688 * mrSges(2,2) + (m(2) + m(3)) * t666 + (mrSges(3,1) + mrSges(4,1)) * t657 + t721 * t656 + ((t663 + t664) * t686 + t708 * t682) * qJD(1) - t692;
t572 = m(3) * t638 + t657 * mrSges(3,3) + t708 * qJD(2) + t721 * qJDD(2) + t655 * t706 + t695;
t571 = m(3) * t637 + qJDD(2) * mrSges(3,1) + qJD(2) * t664 + (-mrSges(3,3) - mrSges(4,3)) * t656 + (-t654 - t655) * t707 + t700;
t569 = m(2) * t667 - t688 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t682 * t571 + t686 * t572;
t568 = mrSges(2,1) * g(3) + mrSges(2,3) * t667 + t688 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - t726;
t563 = -mrSges(2,2) * g(3) - mrSges(2,3) * t666 + Ifges(2,5) * qJDD(1) - t688 * Ifges(2,6) - t682 * t565 + t686 * t567;
t1 = [-m(1) * g(1) + t687 * t569 - t683 * t575; -m(1) * g(2) + t683 * t569 + t687 * t575; t686 * t571 + t682 * t572 + (-m(1) - m(2)) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t687 * t563 - t683 * t568; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t683 * t563 + t687 * t568; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t696; t696; t726; t576; -t691; t690;];
tauJB = t1;
