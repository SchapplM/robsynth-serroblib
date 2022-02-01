% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:46
% EndTime: 2022-01-23 09:22:51
% DurationCPUTime: 5.37s
% Computational Cost: add. (62783->270), mult. (135182->344), div. (0->0), fcn. (85870->10), ass. (0->109)
t669 = sin(qJ(1));
t672 = cos(qJ(1));
t650 = t669 * g(1) - t672 * g(2);
t641 = qJDD(1) * pkin(1) + t650;
t651 = -t672 * g(1) - t669 * g(2);
t673 = qJD(1) ^ 2;
t643 = -t673 * pkin(1) + t651;
t664 = sin(pkin(8));
t666 = cos(pkin(8));
t621 = t664 * t641 + t666 * t643;
t615 = -t673 * pkin(2) + qJDD(1) * pkin(6) + t621;
t662 = -g(3) + qJDD(2);
t668 = sin(qJ(3));
t671 = cos(qJ(3));
t604 = -t668 * t615 + t671 * t662;
t687 = qJD(1) * qJD(3);
t686 = t671 * t687;
t644 = t668 * qJDD(1) + t686;
t601 = (-t644 + t686) * qJ(4) + (t668 * t671 * t673 + qJDD(3)) * pkin(3) + t604;
t605 = t671 * t615 + t668 * t662;
t645 = t671 * qJDD(1) - t668 * t687;
t689 = qJD(1) * t668;
t647 = qJD(3) * pkin(3) - qJ(4) * t689;
t661 = t671 ^ 2;
t602 = -t661 * t673 * pkin(3) + t645 * qJ(4) - qJD(3) * t647 + t605;
t663 = sin(pkin(9));
t665 = cos(pkin(9));
t631 = (t663 * t671 + t665 * t668) * qJD(1);
t581 = -0.2e1 * qJD(4) * t631 + t665 * t601 - t663 * t602;
t623 = t665 * t644 + t663 * t645;
t630 = (-t663 * t668 + t665 * t671) * qJD(1);
t579 = (qJD(3) * t630 - t623) * pkin(7) + (t630 * t631 + qJDD(3)) * pkin(4) + t581;
t582 = 0.2e1 * qJD(4) * t630 + t663 * t601 + t665 * t602;
t622 = -t663 * t644 + t665 * t645;
t626 = qJD(3) * pkin(4) - t631 * pkin(7);
t629 = t630 ^ 2;
t580 = -t629 * pkin(4) + t622 * pkin(7) - qJD(3) * t626 + t582;
t667 = sin(qJ(5));
t670 = cos(qJ(5));
t577 = t670 * t579 - t667 * t580;
t612 = t670 * t630 - t667 * t631;
t591 = t612 * qJD(5) + t667 * t622 + t670 * t623;
t613 = t667 * t630 + t670 * t631;
t600 = -t612 * mrSges(6,1) + t613 * mrSges(6,2);
t658 = qJD(3) + qJD(5);
t606 = -t658 * mrSges(6,2) + t612 * mrSges(6,3);
t657 = qJDD(3) + qJDD(5);
t573 = m(6) * t577 + t657 * mrSges(6,1) - t591 * mrSges(6,3) - t613 * t600 + t658 * t606;
t578 = t667 * t579 + t670 * t580;
t590 = -t613 * qJD(5) + t670 * t622 - t667 * t623;
t607 = t658 * mrSges(6,1) - t613 * mrSges(6,3);
t574 = m(6) * t578 - t657 * mrSges(6,2) + t590 * mrSges(6,3) + t612 * t600 - t658 * t607;
t564 = t670 * t573 + t667 * t574;
t617 = -t630 * mrSges(5,1) + t631 * mrSges(5,2);
t624 = -qJD(3) * mrSges(5,2) + t630 * mrSges(5,3);
t562 = m(5) * t581 + qJDD(3) * mrSges(5,1) - t623 * mrSges(5,3) + qJD(3) * t624 - t631 * t617 + t564;
t625 = qJD(3) * mrSges(5,1) - t631 * mrSges(5,3);
t681 = -t667 * t573 + t670 * t574;
t563 = m(5) * t582 - qJDD(3) * mrSges(5,2) + t622 * mrSges(5,3) - qJD(3) * t625 + t630 * t617 + t681;
t558 = t665 * t562 + t663 * t563;
t610 = Ifges(5,4) * t631 + Ifges(5,2) * t630 + Ifges(5,6) * qJD(3);
t611 = Ifges(5,1) * t631 + Ifges(5,4) * t630 + Ifges(5,5) * qJD(3);
t636 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t668 + Ifges(4,2) * t671) * qJD(1);
t637 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t668 + Ifges(4,4) * t671) * qJD(1);
t593 = Ifges(6,4) * t613 + Ifges(6,2) * t612 + Ifges(6,6) * t658;
t594 = Ifges(6,1) * t613 + Ifges(6,4) * t612 + Ifges(6,5) * t658;
t677 = -mrSges(6,1) * t577 + mrSges(6,2) * t578 - Ifges(6,5) * t591 - Ifges(6,6) * t590 - Ifges(6,3) * t657 - t613 * t593 + t612 * t594;
t692 = mrSges(4,1) * t604 + mrSges(5,1) * t581 - mrSges(4,2) * t605 - mrSges(5,2) * t582 + Ifges(4,5) * t644 + Ifges(5,5) * t623 + Ifges(4,6) * t645 + Ifges(5,6) * t622 + pkin(3) * t558 + pkin(4) * t564 + (t668 * t636 - t671 * t637) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t631 * t610 - t630 * t611 - t677;
t642 = (-mrSges(4,1) * t671 + mrSges(4,2) * t668) * qJD(1);
t688 = qJD(1) * t671;
t649 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t688;
t556 = m(4) * t604 + qJDD(3) * mrSges(4,1) - t644 * mrSges(4,3) + qJD(3) * t649 - t642 * t689 + t558;
t648 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t689;
t682 = -t663 * t562 + t665 * t563;
t557 = m(4) * t605 - qJDD(3) * mrSges(4,2) + t645 * mrSges(4,3) - qJD(3) * t648 + t642 * t688 + t682;
t683 = -t668 * t556 + t671 * t557;
t547 = m(3) * t621 - t673 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t683;
t620 = t666 * t641 - t664 * t643;
t678 = -qJDD(1) * pkin(2) - t620;
t603 = -t645 * pkin(3) + qJDD(4) + t647 * t689 + (-qJ(4) * t661 - pkin(6)) * t673 + t678;
t584 = -t622 * pkin(4) - t629 * pkin(7) + t631 * t626 + t603;
t680 = m(6) * t584 - t590 * mrSges(6,1) + t591 * mrSges(6,2) - t612 * t606 + t613 * t607;
t575 = m(5) * t603 - t622 * mrSges(5,1) + t623 * mrSges(5,2) - t630 * t624 + t631 * t625 + t680;
t614 = -t673 * pkin(6) + t678;
t675 = -m(4) * t614 + t645 * mrSges(4,1) - t644 * mrSges(4,2) - t648 * t689 + t649 * t688 - t575;
t568 = m(3) * t620 + qJDD(1) * mrSges(3,1) - t673 * mrSges(3,2) + t675;
t544 = t664 * t547 + t666 * t568;
t541 = m(2) * t650 + qJDD(1) * mrSges(2,1) - t673 * mrSges(2,2) + t544;
t684 = t666 * t547 - t664 * t568;
t542 = m(2) * t651 - t673 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t684;
t690 = t672 * t541 + t669 * t542;
t550 = t671 * t556 + t668 * t557;
t548 = m(3) * t662 + t550;
t685 = -t669 * t541 + t672 * t542;
t592 = Ifges(6,5) * t613 + Ifges(6,6) * t612 + Ifges(6,3) * t658;
t565 = -mrSges(6,1) * t584 + mrSges(6,3) * t578 + Ifges(6,4) * t591 + Ifges(6,2) * t590 + Ifges(6,6) * t657 - t613 * t592 + t658 * t594;
t566 = mrSges(6,2) * t584 - mrSges(6,3) * t577 + Ifges(6,1) * t591 + Ifges(6,4) * t590 + Ifges(6,5) * t657 + t612 * t592 - t658 * t593;
t609 = Ifges(5,5) * t631 + Ifges(5,6) * t630 + Ifges(5,3) * qJD(3);
t551 = -mrSges(5,1) * t603 + mrSges(5,3) * t582 + Ifges(5,4) * t623 + Ifges(5,2) * t622 + Ifges(5,6) * qJDD(3) - pkin(4) * t680 + pkin(7) * t681 + qJD(3) * t611 + t670 * t565 + t667 * t566 - t631 * t609;
t552 = mrSges(5,2) * t603 - mrSges(5,3) * t581 + Ifges(5,1) * t623 + Ifges(5,4) * t622 + Ifges(5,5) * qJDD(3) - pkin(7) * t564 - qJD(3) * t610 - t667 * t565 + t670 * t566 + t630 * t609;
t635 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t668 + Ifges(4,6) * t671) * qJD(1);
t535 = -mrSges(4,1) * t614 + mrSges(4,3) * t605 + Ifges(4,4) * t644 + Ifges(4,2) * t645 + Ifges(4,6) * qJDD(3) - pkin(3) * t575 + qJ(4) * t682 + qJD(3) * t637 + t665 * t551 + t663 * t552 - t635 * t689;
t537 = mrSges(4,2) * t614 - mrSges(4,3) * t604 + Ifges(4,1) * t644 + Ifges(4,4) * t645 + Ifges(4,5) * qJDD(3) - qJ(4) * t558 - qJD(3) * t636 - t663 * t551 + t665 * t552 + t635 * t688;
t676 = mrSges(2,1) * t650 + mrSges(3,1) * t620 - mrSges(2,2) * t651 - mrSges(3,2) * t621 + pkin(1) * t544 + pkin(2) * t675 + pkin(6) * t683 + t671 * t535 + t668 * t537 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t533 = -mrSges(3,1) * t662 + mrSges(3,3) * t621 + t673 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t550 - t692;
t532 = mrSges(3,2) * t662 - mrSges(3,3) * t620 + Ifges(3,5) * qJDD(1) - t673 * Ifges(3,6) - pkin(6) * t550 - t668 * t535 + t671 * t537;
t531 = -mrSges(2,2) * g(3) - mrSges(2,3) * t650 + Ifges(2,5) * qJDD(1) - t673 * Ifges(2,6) - qJ(2) * t544 + t666 * t532 - t664 * t533;
t530 = mrSges(2,1) * g(3) + mrSges(2,3) * t651 + t673 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t548 + qJ(2) * t684 + t664 * t532 + t666 * t533;
t1 = [-m(1) * g(1) + t685; -m(1) * g(2) + t690; (-m(1) - m(2)) * g(3) + t548; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t690 - t669 * t530 + t672 * t531; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t685 + t672 * t530 + t669 * t531; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t676; t676; t548; t692; t575; -t677;];
tauJB = t1;
