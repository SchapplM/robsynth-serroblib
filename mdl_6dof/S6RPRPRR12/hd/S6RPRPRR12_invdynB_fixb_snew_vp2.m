% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-05-05 20:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:45:03
% EndTime: 2019-05-05 20:45:09
% DurationCPUTime: 3.66s
% Computational Cost: add. (38158->321), mult. (76169->377), div. (0->0), fcn. (42329->8), ass. (0->132)
t717 = -2 * qJD(4);
t716 = Ifges(4,1) + Ifges(5,2);
t704 = (Ifges(4,5) - Ifges(5,4));
t715 = Ifges(4,2) + Ifges(5,3);
t702 = (Ifges(4,6) - Ifges(5,5));
t701 = -Ifges(5,6) - Ifges(4,4);
t714 = (-Ifges(4,3) - Ifges(5,1));
t663 = sin(qJ(1));
t667 = cos(qJ(1));
t642 = -t667 * g(1) - t663 * g(2);
t681 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t642;
t669 = qJD(1) ^ 2;
t711 = (-pkin(1) - pkin(7));
t692 = t711 * t669;
t608 = t692 + t681;
t662 = sin(qJ(3));
t666 = cos(qJ(3));
t695 = qJD(1) * qJD(3);
t691 = t666 * t695;
t632 = qJDD(1) * t662 + t691;
t648 = t662 * t695;
t633 = qJDD(1) * t666 - t648;
t696 = qJD(1) * t662;
t636 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t696;
t649 = t666 * qJD(1);
t637 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t649;
t674 = pkin(3) * t691 + t649 * t717 + t681 + (-t633 + t648) * qJ(4);
t586 = t632 * pkin(3) + t674 + t692;
t638 = mrSges(5,1) * t696 - (qJD(3) * mrSges(5,3));
t639 = mrSges(5,1) * t649 + (qJD(3) * mrSges(5,2));
t640 = pkin(4) * t649 - (qJD(3) * pkin(8));
t657 = t662 ^ 2;
t710 = pkin(3) + pkin(8);
t574 = -t640 * t649 + t710 * t632 + (-pkin(4) * t657 + t711) * t669 + t674;
t629 = (pkin(3) * t662 - qJ(4) * t666) * qJD(1);
t668 = qJD(3) ^ 2;
t641 = t663 * g(1) - t667 * g(2);
t680 = -t669 * qJ(2) + qJDD(2) - t641;
t609 = t711 * qJDD(1) + t680;
t700 = t666 * t609;
t679 = -t668 * qJ(4) + t629 * t649 + qJDD(4) - t700;
t709 = pkin(8) * t669;
t579 = t633 * pkin(4) - t710 * qJDD(3) + (pkin(4) * t695 + t666 * t709 - g(3)) * t662 + t679;
t661 = sin(qJ(5));
t665 = cos(qJ(5));
t566 = -t661 * t574 + t665 * t579;
t627 = -qJD(3) * t661 + t665 * t696;
t597 = qJD(5) * t627 + qJDD(3) * t665 + t632 * t661;
t626 = qJDD(5) + t633;
t628 = qJD(3) * t665 + t661 * t696;
t646 = t649 + qJD(5);
t564 = (t627 * t646 - t597) * pkin(9) + (t627 * t628 + t626) * pkin(5) + t566;
t567 = t665 * t574 + t661 * t579;
t596 = -qJD(5) * t628 - qJDD(3) * t661 + t632 * t665;
t607 = pkin(5) * t646 - pkin(9) * t628;
t625 = t627 ^ 2;
t565 = -pkin(5) * t625 + pkin(9) * t596 - t607 * t646 + t567;
t660 = sin(qJ(6));
t664 = cos(qJ(6));
t562 = t564 * t664 - t565 * t660;
t598 = t627 * t664 - t628 * t660;
t572 = qJD(6) * t598 + t596 * t660 + t597 * t664;
t599 = t627 * t660 + t628 * t664;
t585 = -mrSges(7,1) * t598 + mrSges(7,2) * t599;
t643 = qJD(6) + t646;
t589 = -mrSges(7,2) * t643 + mrSges(7,3) * t598;
t619 = qJDD(6) + t626;
t560 = m(7) * t562 + mrSges(7,1) * t619 - mrSges(7,3) * t572 - t585 * t599 + t589 * t643;
t563 = t564 * t660 + t565 * t664;
t571 = -qJD(6) * t599 + t596 * t664 - t597 * t660;
t590 = mrSges(7,1) * t643 - mrSges(7,3) * t599;
t561 = m(7) * t563 - mrSges(7,2) * t619 + mrSges(7,3) * t571 + t585 * t598 - t590 * t643;
t552 = t664 * t560 + t660 * t561;
t600 = -mrSges(6,1) * t627 + mrSges(6,2) * t628;
t603 = -mrSges(6,2) * t646 + mrSges(6,3) * t627;
t550 = m(6) * t566 + mrSges(6,1) * t626 - mrSges(6,3) * t597 - t600 * t628 + t603 * t646 + t552;
t604 = mrSges(6,1) * t646 - mrSges(6,3) * t628;
t687 = -t560 * t660 + t664 * t561;
t551 = m(6) * t567 - mrSges(6,2) * t626 + mrSges(6,3) * t596 + t600 * t627 - t604 * t646 + t687;
t688 = -t661 * t550 + t665 * t551;
t673 = m(5) * t586 - t633 * mrSges(5,3) - (t638 * t662 + t639 * t666) * qJD(1) + t688;
t706 = mrSges(4,1) - mrSges(5,2);
t713 = -m(4) * t608 - t633 * mrSges(4,2) - t706 * t632 - t636 * t696 - t637 * t649 - t673;
t602 = -g(3) * t666 + t662 * t609;
t587 = pkin(3) * t668 - qJDD(3) * qJ(4) + (qJD(3) * t717) + t629 * t696 - t602;
t708 = t662 * g(3);
t707 = mrSges(2,1) - mrSges(3,2);
t705 = Ifges(2,5) - Ifges(3,4);
t703 = (-Ifges(2,6) + Ifges(3,5));
t601 = t700 + t708;
t548 = t665 * t550 + t661 * t551;
t588 = -qJDD(3) * pkin(3) + t679 - t708;
t676 = -m(5) * t588 - t633 * mrSges(5,1) - t548;
t630 = (-mrSges(5,2) * t662 - mrSges(5,3) * t666) * qJD(1);
t685 = qJD(1) * (-t630 - (mrSges(4,1) * t662 + mrSges(4,2) * t666) * qJD(1));
t546 = m(4) * t601 - t633 * mrSges(4,3) + t706 * qJDD(3) + (t636 - t638) * qJD(3) + t666 * t685 + t676;
t578 = -pkin(4) * t632 + qJD(3) * t640 - t657 * t709 - t587;
t569 = -pkin(5) * t596 - pkin(9) * t625 + t607 * t628 + t578;
t678 = m(7) * t569 - mrSges(7,1) * t571 + t572 * mrSges(7,2) - t589 * t598 + t599 * t590;
t672 = -m(6) * t578 + mrSges(6,1) * t596 - t597 * mrSges(6,2) + t603 * t627 - t628 * t604 - t678;
t671 = -m(5) * t587 + qJDD(3) * mrSges(5,3) + qJD(3) * t639 - t672;
t556 = m(4) * t602 - qJD(3) * t637 - qJDD(3) * mrSges(4,2) + (-mrSges(4,3) - mrSges(5,1)) * t632 + t662 * t685 + t671;
t540 = t666 * t546 + t662 * t556;
t611 = -qJDD(1) * pkin(1) + t680;
t677 = -m(3) * t611 + (t669 * mrSges(3,3)) - t540;
t538 = m(2) * t641 - (t669 * mrSges(2,2)) + t707 * qJDD(1) + t677;
t610 = t669 * pkin(1) - t681;
t670 = -m(3) * t610 + (t669 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t713;
t545 = m(2) * t642 - (mrSges(2,1) * t669) - qJDD(1) * mrSges(2,2) + t670;
t699 = t667 * t538 + t663 * t545;
t698 = -(t702 * qJD(3)) + (t715 * t662 + t701 * t666) * qJD(1);
t697 = (t704 * qJD(3)) + (t701 * t662 + t716 * t666) * qJD(1);
t690 = -t538 * t663 + t667 * t545;
t689 = -t662 * t546 + t666 * t556;
t686 = qJD(1) * ((t714 * qJD(3)) + (t702 * t662 - t704 * t666) * qJD(1));
t593 = Ifges(6,1) * t628 + Ifges(6,4) * t627 + Ifges(6,5) * t646;
t592 = Ifges(6,4) * t628 + Ifges(6,2) * t627 + Ifges(6,6) * t646;
t591 = Ifges(6,5) * t628 + Ifges(6,6) * t627 + Ifges(6,3) * t646;
t582 = Ifges(7,1) * t599 + Ifges(7,4) * t598 + Ifges(7,5) * t643;
t581 = Ifges(7,4) * t599 + Ifges(7,2) * t598 + Ifges(7,6) * t643;
t580 = Ifges(7,5) * t599 + Ifges(7,6) * t598 + Ifges(7,3) * t643;
t554 = mrSges(7,2) * t569 - mrSges(7,3) * t562 + Ifges(7,1) * t572 + Ifges(7,4) * t571 + Ifges(7,5) * t619 + t580 * t598 - t581 * t643;
t553 = -mrSges(7,1) * t569 + mrSges(7,3) * t563 + Ifges(7,4) * t572 + Ifges(7,2) * t571 + Ifges(7,6) * t619 - t580 * t599 + t582 * t643;
t547 = -t632 * mrSges(5,2) + t673;
t542 = mrSges(6,2) * t578 - mrSges(6,3) * t566 + Ifges(6,1) * t597 + Ifges(6,4) * t596 + Ifges(6,5) * t626 - pkin(9) * t552 - t553 * t660 + t554 * t664 + t591 * t627 - t592 * t646;
t541 = -mrSges(6,1) * t578 + mrSges(6,3) * t567 + Ifges(6,4) * t597 + Ifges(6,2) * t596 + Ifges(6,6) * t626 - pkin(5) * t678 + pkin(9) * t687 + t664 * t553 + t660 * t554 - t628 * t591 + t646 * t593;
t539 = -m(3) * g(3) + t689;
t536 = t698 * qJD(3) + t716 * t633 + t701 * t632 + t704 * qJDD(3) + t628 * t592 + Ifges(7,3) * t619 + Ifges(6,3) * t626 - t627 * t593 + t599 * t581 - mrSges(4,3) * t601 + mrSges(4,2) * t608 + Ifges(6,6) * t596 + Ifges(6,5) * t597 - t598 * t582 - mrSges(5,3) * t586 + mrSges(5,1) * t588 + mrSges(6,1) * t566 - mrSges(6,2) * t567 + Ifges(7,6) * t571 + Ifges(7,5) * t572 + mrSges(7,1) * t562 - mrSges(7,2) * t563 + pkin(5) * t552 + pkin(4) * t548 - qJ(4) * t547 + t662 * t686;
t535 = -mrSges(4,1) * t608 - mrSges(5,1) * t587 + mrSges(5,2) * t586 + mrSges(4,3) * t602 - pkin(3) * t547 - pkin(4) * t672 - pkin(8) * t688 + t697 * qJD(3) + t702 * qJDD(3) - t665 * t541 - t661 * t542 - t715 * t632 - t701 * t633 + t666 * t686;
t534 = qJ(4) * t671 + pkin(3) * (-qJD(3) * t638 + t676) - t661 * t541 + t665 * t542 - mrSges(2,3) * t641 + mrSges(4,1) * t601 - mrSges(4,2) * t602 + mrSges(3,1) * t611 - mrSges(5,3) * t587 + mrSges(5,2) * t588 - pkin(8) * t548 + pkin(2) * t540 - qJ(2) * t539 + (t703 * t669) + t704 * t633 + (-mrSges(5,1) * qJ(4) - t702) * t632 + (-mrSges(5,2) * pkin(3) - t714) * qJDD(3) + t705 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + ((-pkin(3) * t630 - t698) * t666 + (-qJ(4) * t630 + t697) * t662) * qJD(1);
t533 = -mrSges(3,1) * t610 + mrSges(2,3) * t642 - pkin(1) * t539 - pkin(2) * t713 - pkin(7) * t689 + t707 * g(3) - t703 * qJDD(1) - t666 * t535 - t662 * t536 + t705 * t669;
t1 = [-m(1) * g(1) + t690; -m(1) * g(2) + t699; (-m(1) - m(2) - m(3)) * g(3) + t689; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t699 - t663 * t533 + t667 * t534; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t690 + t667 * t533 + t663 * t534; pkin(1) * t677 + qJ(2) * t670 - t662 * t535 - pkin(7) * t540 + mrSges(2,1) * t641 - mrSges(2,2) * t642 + t666 * t536 + mrSges(3,2) * t611 - mrSges(3,3) * t610 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
