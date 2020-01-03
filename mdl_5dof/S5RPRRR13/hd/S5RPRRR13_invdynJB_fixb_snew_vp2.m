% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR13
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:52
% EndTime: 2019-12-31 19:14:56
% DurationCPUTime: 3.03s
% Computational Cost: add. (33982->265), mult. (65535->321), div. (0->0), fcn. (39727->8), ass. (0->110)
t646 = sin(qJ(1));
t650 = cos(qJ(1));
t628 = -t650 * g(1) - t646 * g(2);
t677 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t628;
t676 = -pkin(1) - pkin(6);
t675 = mrSges(2,1) - mrSges(3,2);
t674 = -Ifges(3,4) + Ifges(2,5);
t673 = (Ifges(3,5) - Ifges(2,6));
t627 = t646 * g(1) - t650 * g(2);
t652 = qJD(1) ^ 2;
t662 = -t652 * qJ(2) + qJDD(2) - t627;
t599 = t676 * qJDD(1) + t662;
t645 = sin(qJ(3));
t649 = cos(qJ(3));
t592 = -t649 * g(3) + t645 * t599;
t620 = (mrSges(4,1) * t645 + mrSges(4,2) * t649) * qJD(1);
t670 = qJD(1) * qJD(3);
t631 = t649 * t670;
t622 = -t645 * qJDD(1) - t631;
t671 = qJD(1) * t649;
t626 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t671;
t633 = t645 * qJD(1);
t598 = t676 * t652 - t677;
t668 = t645 * t670;
t623 = t649 * qJDD(1) - t668;
t574 = (-t623 + t668) * pkin(7) + (-t622 + t631) * pkin(3) + t598;
t621 = (pkin(3) * t645 - pkin(7) * t649) * qJD(1);
t651 = qJD(3) ^ 2;
t577 = -t651 * pkin(3) + qJDD(3) * pkin(7) - t621 * t633 + t592;
t644 = sin(qJ(4));
t648 = cos(qJ(4));
t559 = t648 * t574 - t644 * t577;
t618 = t648 * qJD(3) - t644 * t671;
t586 = t618 * qJD(4) + t644 * qJDD(3) + t648 * t623;
t617 = qJDD(4) - t622;
t619 = t644 * qJD(3) + t648 * t671;
t630 = t633 + qJD(4);
t556 = (t618 * t630 - t586) * pkin(8) + (t618 * t619 + t617) * pkin(4) + t559;
t560 = t644 * t574 + t648 * t577;
t585 = -t619 * qJD(4) + t648 * qJDD(3) - t644 * t623;
t597 = t630 * pkin(4) - t619 * pkin(8);
t616 = t618 ^ 2;
t557 = -t616 * pkin(4) + t585 * pkin(8) - t630 * t597 + t560;
t643 = sin(qJ(5));
t647 = cos(qJ(5));
t554 = t647 * t556 - t643 * t557;
t587 = t647 * t618 - t643 * t619;
t566 = t587 * qJD(5) + t643 * t585 + t647 * t586;
t588 = t643 * t618 + t647 * t619;
t571 = -t587 * mrSges(6,1) + t588 * mrSges(6,2);
t629 = qJD(5) + t630;
t578 = -t629 * mrSges(6,2) + t587 * mrSges(6,3);
t610 = qJDD(5) + t617;
t551 = m(6) * t554 + t610 * mrSges(6,1) - t566 * mrSges(6,3) - t588 * t571 + t629 * t578;
t555 = t643 * t556 + t647 * t557;
t565 = -t588 * qJD(5) + t647 * t585 - t643 * t586;
t579 = t629 * mrSges(6,1) - t588 * mrSges(6,3);
t552 = m(6) * t555 - t610 * mrSges(6,2) + t565 * mrSges(6,3) + t587 * t571 - t629 * t579;
t543 = t647 * t551 + t643 * t552;
t590 = -t618 * mrSges(5,1) + t619 * mrSges(5,2);
t593 = -t630 * mrSges(5,2) + t618 * mrSges(5,3);
t541 = m(5) * t559 + t617 * mrSges(5,1) - t586 * mrSges(5,3) - t619 * t590 + t630 * t593 + t543;
t594 = t630 * mrSges(5,1) - t619 * mrSges(5,3);
t664 = -t643 * t551 + t647 * t552;
t542 = m(5) * t560 - t617 * mrSges(5,2) + t585 * mrSges(5,3) + t618 * t590 - t630 * t594 + t664;
t665 = -t644 * t541 + t648 * t542;
t535 = m(4) * t592 - qJDD(3) * mrSges(4,2) + t622 * mrSges(4,3) - qJD(3) * t626 - t620 * t633 + t665;
t591 = t645 * g(3) + t649 * t599;
t625 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t633;
t576 = -qJDD(3) * pkin(3) - t651 * pkin(7) + t621 * t671 - t591;
t558 = -t585 * pkin(4) - t616 * pkin(8) + t619 * t597 + t576;
t659 = m(6) * t558 - t565 * mrSges(6,1) + t566 * mrSges(6,2) - t587 * t578 + t588 * t579;
t654 = -m(5) * t576 + t585 * mrSges(5,1) - t586 * mrSges(5,2) + t618 * t593 - t619 * t594 - t659;
t546 = m(4) * t591 + qJDD(3) * mrSges(4,1) - t623 * mrSges(4,3) + qJD(3) * t625 - t620 * t671 + t654;
t527 = t645 * t535 + t649 * t546;
t604 = -qJDD(1) * pkin(1) + t662;
t661 = -m(3) * t604 + (t652 * mrSges(3,3)) - t527;
t523 = m(2) * t627 - (t652 * mrSges(2,2)) + t675 * qJDD(1) + t661;
t602 = t652 * pkin(1) + t677;
t537 = t648 * t541 + t644 * t542;
t660 = -m(4) * t598 + t622 * mrSges(4,1) - t623 * mrSges(4,2) - t625 * t633 - t626 * t671 - t537;
t656 = -m(3) * t602 + (t652 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t660;
t532 = m(2) * t628 - (t652 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t656;
t672 = t650 * t523 + t646 * t532;
t667 = -t646 * t523 + t650 * t532;
t666 = t649 * t535 - t645 * t546;
t568 = Ifges(6,4) * t588 + Ifges(6,2) * t587 + Ifges(6,6) * t629;
t569 = Ifges(6,1) * t588 + Ifges(6,4) * t587 + Ifges(6,5) * t629;
t658 = -mrSges(6,1) * t554 + mrSges(6,2) * t555 - Ifges(6,5) * t566 - Ifges(6,6) * t565 - Ifges(6,3) * t610 - t588 * t568 + t587 * t569;
t567 = Ifges(6,5) * t588 + Ifges(6,6) * t587 + Ifges(6,3) * t629;
t544 = -mrSges(6,1) * t558 + mrSges(6,3) * t555 + Ifges(6,4) * t566 + Ifges(6,2) * t565 + Ifges(6,6) * t610 - t588 * t567 + t629 * t569;
t545 = mrSges(6,2) * t558 - mrSges(6,3) * t554 + Ifges(6,1) * t566 + Ifges(6,4) * t565 + Ifges(6,5) * t610 + t587 * t567 - t629 * t568;
t580 = Ifges(5,5) * t619 + Ifges(5,6) * t618 + Ifges(5,3) * t630;
t582 = Ifges(5,1) * t619 + Ifges(5,4) * t618 + Ifges(5,5) * t630;
t521 = -mrSges(5,1) * t576 + mrSges(5,3) * t560 + Ifges(5,4) * t586 + Ifges(5,2) * t585 + Ifges(5,6) * t617 - pkin(4) * t659 + pkin(8) * t664 + t647 * t544 + t643 * t545 - t619 * t580 + t630 * t582;
t581 = Ifges(5,4) * t619 + Ifges(5,2) * t618 + Ifges(5,6) * t630;
t529 = mrSges(5,2) * t576 - mrSges(5,3) * t559 + Ifges(5,1) * t586 + Ifges(5,4) * t585 + Ifges(5,5) * t617 - pkin(8) * t543 - t643 * t544 + t647 * t545 + t618 * t580 - t630 * t581;
t608 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t649 - Ifges(4,2) * t645) * qJD(1);
t609 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t649 - Ifges(4,4) * t645) * qJD(1);
t657 = mrSges(4,1) * t591 - mrSges(4,2) * t592 + Ifges(4,5) * t623 + Ifges(4,6) * t622 + Ifges(4,3) * qJDD(3) + pkin(3) * t654 + pkin(7) * t665 + t648 * t521 + t644 * t529 + t608 * t671 + t609 * t633;
t607 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t649 - Ifges(4,6) * t645) * qJD(1);
t518 = mrSges(4,2) * t598 - mrSges(4,3) * t591 + Ifges(4,1) * t623 + Ifges(4,4) * t622 + Ifges(4,5) * qJDD(3) - pkin(7) * t537 - qJD(3) * t608 - t644 * t521 + t648 * t529 - t607 * t633;
t653 = mrSges(5,1) * t559 - mrSges(5,2) * t560 + Ifges(5,5) * t586 + Ifges(5,6) * t585 + Ifges(5,3) * t617 + pkin(4) * t543 + t619 * t581 - t618 * t582 - t658;
t519 = -mrSges(4,1) * t598 + mrSges(4,3) * t592 + Ifges(4,4) * t623 + Ifges(4,2) * t622 + Ifges(4,6) * qJDD(3) - pkin(3) * t537 + qJD(3) * t609 - t607 * t671 - t653;
t525 = qJDD(1) * mrSges(3,2) - t661;
t655 = mrSges(2,1) * t627 - mrSges(2,2) * t628 + mrSges(3,2) * t604 - mrSges(3,3) * t602 - pkin(1) * t525 - pkin(6) * t527 + qJ(2) * t656 + t649 * t518 - t645 * t519 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t526 = -m(3) * g(3) + t666;
t516 = t657 + mrSges(3,1) * t604 + pkin(2) * t527 - qJ(2) * t526 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t674 * qJDD(1) - mrSges(2,3) * t627 + (t673 * t652);
t515 = -mrSges(3,1) * t602 + mrSges(2,3) * t628 - pkin(1) * t526 - pkin(2) * t660 - pkin(6) * t666 + t675 * g(3) - t673 * qJDD(1) - t645 * t518 - t649 * t519 + t674 * t652;
t1 = [-m(1) * g(1) + t667; -m(1) * g(2) + t672; (-m(1) - m(2) - m(3)) * g(3) + t666; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t672 - t646 * t515 + t650 * t516; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t667 + t650 * t515 + t646 * t516; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t655; t655; t525; t657; t653; -t658;];
tauJB = t1;
