% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:25:40
% EndTime: 2019-05-05 21:25:46
% DurationCPUTime: 3.74s
% Computational Cost: add. (37462->299), mult. (70438->344), div. (0->0), fcn. (41025->8), ass. (0->117)
t665 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t651 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t650 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t664 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t663 = -Ifges(5,6) + Ifges(6,5) + Ifges(7,4);
t662 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t661 = -2 * qJD(5);
t660 = 2 * qJD(6);
t659 = cos(qJ(4));
t658 = -mrSges(7,1) - mrSges(5,3);
t622 = sin(qJ(4));
t623 = sin(qJ(3));
t654 = qJD(1) * t623;
t599 = -qJD(3) * t659 + t622 * t654;
t625 = cos(qJ(3));
t653 = t625 * qJD(1);
t612 = -qJD(4) + t653;
t657 = t599 * t612;
t624 = sin(qJ(1));
t626 = cos(qJ(1));
t609 = t624 * g(1) - t626 * g(2);
t601 = qJDD(1) * pkin(1) + t609;
t610 = -t626 * g(1) - t624 * g(2);
t628 = qJD(1) ^ 2;
t603 = -t628 * pkin(1) + t610;
t620 = sin(pkin(9));
t621 = cos(pkin(9));
t568 = t620 * t601 + t621 * t603;
t546 = -t628 * pkin(2) + qJDD(1) * pkin(7) + t568;
t619 = -g(3) + qJDD(2);
t540 = t625 * t546 + t623 * t619;
t602 = (-mrSges(4,1) * t625 + mrSges(4,2) * t623) * qJD(1);
t652 = qJD(1) * qJD(3);
t643 = t623 * t652;
t606 = t625 * qJDD(1) - t643;
t607 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t654;
t567 = t621 * t601 - t620 * t603;
t545 = -qJDD(1) * pkin(2) - t628 * pkin(7) - t567;
t642 = t625 * t652;
t605 = t623 * qJDD(1) + t642;
t534 = (-t605 - t642) * pkin(8) + (-t606 + t643) * pkin(3) + t545;
t604 = (-pkin(3) * t625 - pkin(8) * t623) * qJD(1);
t627 = qJD(3) ^ 2;
t538 = -t627 * pkin(3) + qJDD(3) * pkin(8) + t604 * t653 + t540;
t531 = t534 * t659 - t622 * t538;
t566 = -t599 * qJD(4) + t622 * qJDD(3) + t605 * t659;
t600 = t622 * qJD(3) + t654 * t659;
t572 = -t600 * mrSges(7,2) + t599 * mrSges(7,3);
t574 = t599 * mrSges(5,1) + t600 * mrSges(5,2);
t576 = t612 * mrSges(5,2) - t599 * mrSges(5,3);
t580 = t599 * mrSges(6,1) + t612 * mrSges(6,3);
t598 = qJDD(4) - t606;
t573 = t599 * pkin(4) - t600 * qJ(5);
t611 = t612 ^ 2;
t529 = -t598 * pkin(4) - t611 * qJ(5) + t600 * t573 + qJDD(5) - t531;
t575 = -t599 * mrSges(6,2) - t600 * mrSges(6,3);
t523 = t612 * t660 + (t599 * t600 - t598) * qJ(6) + (t566 - t657) * pkin(5) + t529;
t581 = -t599 * mrSges(7,1) - t612 * mrSges(7,2);
t637 = m(7) * t523 - t598 * mrSges(7,3) + t612 * t581;
t634 = -m(6) * t529 - t566 * mrSges(6,1) - t600 * t575 - t637;
t518 = m(5) * t531 + (-t576 + t580) * t612 + (-t572 - t574) * t600 + (mrSges(5,1) - mrSges(6,2)) * t598 + t658 * t566 + t634;
t532 = t622 * t534 + t659 * t538;
t565 = t600 * qJD(4) - qJDD(3) * t659 + t622 * t605;
t577 = -t612 * mrSges(5,1) - t600 * mrSges(5,3);
t632 = -t611 * pkin(4) + t598 * qJ(5) - t599 * t573 + t532;
t528 = 0.2e1 * qJD(5) * t612 - t632;
t582 = t600 * mrSges(6,1) - t612 * mrSges(6,2);
t578 = t600 * pkin(5) + t612 * qJ(6);
t597 = t599 ^ 2;
t525 = -t565 * pkin(5) - t597 * qJ(6) + qJDD(6) + (t661 - t578) * t612 + t632;
t579 = t600 * mrSges(7,1) + t612 * mrSges(7,3);
t645 = -m(7) * t525 - t598 * mrSges(7,2) + t612 * t579;
t635 = -m(6) * t528 + t598 * mrSges(6,3) - t612 * t582 - t645;
t655 = -t572 - t575;
t520 = m(5) * t532 - t598 * mrSges(5,2) + t612 * t577 + (-t574 + t655) * t599 + (-mrSges(6,1) + t658) * t565 + t635;
t638 = -t622 * t518 + t659 * t520;
t514 = m(4) * t540 - qJDD(3) * mrSges(4,2) + t606 * mrSges(4,3) - qJD(3) * t607 + t602 * t653 + t638;
t539 = -t623 * t546 + t625 * t619;
t608 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t653;
t537 = -qJDD(3) * pkin(3) - t627 * pkin(8) + t604 * t654 - t539;
t630 = (-t566 - t657) * qJ(5) + t537 + (-t612 * pkin(4) + t661) * t600;
t530 = t565 * pkin(4) + t630;
t527 = -t597 * pkin(5) + t599 * t660 - t600 * t578 + (pkin(4) + qJ(6)) * t565 + t630;
t636 = m(7) * t527 - t566 * mrSges(7,2) + t565 * mrSges(7,3) - t600 * t579 + t599 * t581;
t633 = -m(6) * t530 + t565 * mrSges(6,2) + t599 * t580 - t636;
t629 = -m(5) * t537 - t565 * mrSges(5,1) - t599 * t576 + (-t577 + t582) * t600 + (-mrSges(5,2) + mrSges(6,3)) * t566 + t633;
t517 = m(4) * t539 + qJDD(3) * mrSges(4,1) - t605 * mrSges(4,3) + qJD(3) * t608 - t602 * t654 + t629;
t639 = t625 * t514 - t623 * t517;
t507 = m(3) * t568 - t628 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t639;
t515 = t518 * t659 + t622 * t520;
t631 = -m(4) * t545 + t606 * mrSges(4,1) - t605 * mrSges(4,2) - t607 * t654 + t608 * t653 - t515;
t511 = m(3) * t567 + qJDD(1) * mrSges(3,1) - t628 * mrSges(3,2) + t631;
t503 = t620 * t507 + t621 * t511;
t501 = m(2) * t609 + qJDD(1) * mrSges(2,1) - t628 * mrSges(2,2) + t503;
t640 = t621 * t507 - t620 * t511;
t502 = m(2) * t610 - t628 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t640;
t656 = t626 * t501 + t624 * t502;
t508 = t623 * t514 + t625 * t517;
t648 = -t599 * t663 - t600 * t650 + t612 * t662;
t647 = t599 * t664 + t600 * t651 + t612 * t663;
t646 = t651 * t599 - t600 * t665 + t650 * t612;
t644 = m(3) * t619 + t508;
t641 = -t624 * t501 + t626 * t502;
t591 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t623 + Ifges(4,4) * t625) * qJD(1);
t590 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t623 + Ifges(4,2) * t625) * qJD(1);
t589 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t623 + Ifges(4,6) * t625) * qJD(1);
t522 = t566 * mrSges(7,1) + t600 * t572 + t637;
t521 = -t566 * mrSges(6,3) - t600 * t582 - t633;
t509 = mrSges(6,1) * t529 + mrSges(7,1) * t523 + mrSges(5,2) * t537 - mrSges(7,2) * t527 - mrSges(5,3) * t531 - mrSges(6,3) * t530 + pkin(5) * t522 - qJ(5) * t521 - t651 * t565 + t566 * t665 + t650 * t598 + t648 * t599 + t647 * t612;
t504 = -mrSges(5,1) * t537 + mrSges(5,3) * t532 - mrSges(6,1) * t528 + mrSges(6,2) * t530 + mrSges(7,1) * t525 - mrSges(7,3) * t527 - pkin(5) * (t599 * t572 + t645) - qJ(6) * t636 - pkin(4) * t521 + t646 * t612 + t648 * t600 - t663 * t598 + t651 * t566 + (-pkin(5) * mrSges(7,1) + t664) * t565;
t497 = -t589 * t654 - qJ(5) * t635 + Ifges(4,6) * qJDD(3) - pkin(4) * (t612 * t580 + t634) + Ifges(4,4) * t605 + Ifges(4,2) * t606 + qJD(3) * t591 + mrSges(4,3) * t540 - mrSges(4,1) * t545 + mrSges(6,3) * t528 - mrSges(6,2) * t529 - mrSges(5,1) * t531 + mrSges(5,2) * t532 + qJ(6) * t522 + mrSges(7,3) * t523 - mrSges(7,2) * t525 - pkin(3) * t515 + (pkin(4) * t572 - t647) * t600 + (pkin(4) * mrSges(6,2) - t662) * t598 + (pkin(4) * mrSges(7,1) - t650) * t566 + (-qJ(5) * t655 + t646) * t599 + (-qJ(5) * (-mrSges(6,1) - mrSges(7,1)) - t663) * t565;
t496 = mrSges(4,2) * t545 - mrSges(4,3) * t539 + Ifges(4,1) * t605 + Ifges(4,4) * t606 + Ifges(4,5) * qJDD(3) - pkin(8) * t515 - qJD(3) * t590 - t622 * t504 + t509 * t659 + t589 * t653;
t495 = Ifges(3,6) * qJDD(1) + t628 * Ifges(3,5) - mrSges(3,1) * t619 + mrSges(3,3) * t568 - Ifges(4,5) * t605 - Ifges(4,6) * t606 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t539 + mrSges(4,2) * t540 - t622 * t509 - t659 * t504 - pkin(3) * t629 - pkin(8) * t638 - pkin(2) * t508 + (-t623 * t590 + t625 * t591) * qJD(1);
t494 = mrSges(3,2) * t619 - mrSges(3,3) * t567 + Ifges(3,5) * qJDD(1) - t628 * Ifges(3,6) - pkin(7) * t508 + t625 * t496 - t623 * t497;
t493 = -mrSges(2,2) * g(3) - mrSges(2,3) * t609 + Ifges(2,5) * qJDD(1) - t628 * Ifges(2,6) - qJ(2) * t503 + t621 * t494 - t620 * t495;
t492 = mrSges(2,1) * g(3) + mrSges(2,3) * t610 + t628 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t644 + qJ(2) * t640 + t620 * t494 + t621 * t495;
t1 = [-m(1) * g(1) + t641; -m(1) * g(2) + t656; (-m(1) - m(2)) * g(3) + t644; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t656 - t624 * t492 + t626 * t493; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t641 + t626 * t492 + t624 * t493; pkin(1) * t503 + mrSges(2,1) * t609 - mrSges(2,2) * t610 + t623 * t496 + t625 * t497 + pkin(2) * t631 + pkin(7) * t639 + mrSges(3,1) * t567 - mrSges(3,2) * t568 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
