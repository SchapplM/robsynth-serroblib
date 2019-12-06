% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:20
% EndTime: 2019-12-05 16:12:25
% DurationCPUTime: 2.33s
% Computational Cost: add. (18275->238), mult. (37249->289), div. (0->0), fcn. (21851->8), ass. (0->103)
t645 = Ifges(5,1) + Ifges(6,1);
t636 = Ifges(5,4) - Ifges(6,5);
t644 = Ifges(5,5) + Ifges(6,4);
t643 = Ifges(5,2) + Ifges(6,3);
t642 = Ifges(6,2) + Ifges(5,3);
t641 = Ifges(5,6) - Ifges(6,6);
t603 = sin(qJ(3));
t605 = cos(qJ(3));
t621 = qJD(2) * qJD(3);
t618 = t605 * t621;
t587 = t603 * qJDD(2) + t618;
t601 = sin(pkin(8));
t632 = cos(pkin(8));
t567 = -t632 * qJDD(3) + t601 * t587;
t568 = t601 * qJDD(3) + t632 * t587;
t623 = qJD(2) * t603;
t577 = -t632 * qJD(3) + t601 * t623;
t578 = t601 * qJD(3) + t632 * t623;
t602 = sin(pkin(7));
t633 = cos(pkin(7));
t589 = t602 * g(1) - t633 * g(2);
t590 = -t633 * g(1) - t602 * g(2);
t600 = -g(3) + qJDD(1);
t604 = sin(qJ(2));
t606 = cos(qJ(2));
t570 = t606 * t590 + t604 * t600;
t608 = qJD(2) ^ 2;
t561 = -t608 * pkin(2) + qJDD(2) * pkin(6) + t570;
t558 = t603 * t561;
t585 = (-pkin(3) * t605 - qJ(4) * t603) * qJD(2);
t607 = qJD(3) ^ 2;
t611 = -qJDD(3) * pkin(3) - t607 * qJ(4) + t585 * t623 + qJDD(4) + t558;
t638 = -2 * qJD(5);
t537 = t567 * pkin(4) - t568 * qJ(5) + t578 * t638 + (t589 + (-pkin(4) * t578 - qJ(5) * t577) * qJD(2)) * t605 + t611;
t622 = t605 * qJD(2);
t564 = -t577 * mrSges(6,2) - mrSges(6,3) * t622;
t566 = mrSges(6,1) * t622 + t578 * mrSges(6,2);
t533 = m(6) * t537 + t567 * mrSges(6,1) - t568 * mrSges(6,3) + t577 * t564 - t578 * t566;
t569 = -t604 * t590 + t606 * t600;
t560 = -qJDD(2) * pkin(2) - t608 * pkin(6) - t569;
t619 = t603 * t621;
t588 = t605 * qJDD(2) - t619;
t543 = (-t587 - t618) * qJ(4) + (-t588 + t619) * pkin(3) + t560;
t546 = t605 * t561 - t603 * t589;
t544 = -t607 * pkin(3) + qJDD(3) * qJ(4) + t585 * t622 + t546;
t639 = -2 * qJD(4);
t539 = t601 * t543 + t632 * t544 + t577 * t639;
t555 = t577 * pkin(4) - t578 * qJ(5);
t631 = t605 ^ 2 * t608;
t535 = -pkin(4) * t631 - t588 * qJ(5) - t577 * t555 + t622 * t638 + t539;
t630 = t605 * t589;
t542 = t611 + t630;
t625 = t636 * t577 - t645 * t578 + t644 * t622;
t627 = t641 * t577 - t644 * t578 + t642 * t622;
t523 = -mrSges(5,1) * t542 - mrSges(6,1) * t537 + mrSges(6,2) * t535 + mrSges(5,3) * t539 - pkin(4) * t533 - t643 * t567 + t636 * t568 + t627 * t578 - t588 * t641 + t625 * t622;
t613 = t632 * t543 - t601 * t544;
t536 = -qJ(5) * t631 + t588 * pkin(4) + qJDD(5) + ((2 * qJD(4)) + t555) * t578 - t613;
t538 = t578 * t639 + t613;
t626 = t643 * t577 - t636 * t578 + t641 * t622;
t524 = mrSges(5,2) * t542 + mrSges(6,2) * t536 - mrSges(5,3) * t538 - mrSges(6,3) * t537 - qJ(5) * t533 - t636 * t567 + t645 * t568 + t627 * t577 - t588 * t644 - t626 * t622;
t565 = -mrSges(5,1) * t622 - t578 * mrSges(5,3);
t556 = t577 * mrSges(6,1) - t578 * mrSges(6,3);
t624 = -t577 * mrSges(5,1) - t578 * mrSges(5,2) - t556;
t628 = m(6) * t535 - t588 * mrSges(6,3);
t637 = -mrSges(5,3) - mrSges(6,2);
t528 = m(5) * t539 + t588 * mrSges(5,2) + t624 * t577 + t637 * t567 + (t565 - t566) * t622 + t628;
t612 = -m(6) * t536 - t588 * mrSges(6,1) - t564 * t622;
t614 = mrSges(5,2) * t622 - t577 * mrSges(5,3);
t529 = m(5) * t538 - t588 * mrSges(5,1) + t637 * t568 + t624 * t578 - t614 * t622 + t612;
t526 = t632 * t528 - t601 * t529;
t531 = m(5) * t542 + t567 * mrSges(5,1) + t568 * mrSges(5,2) + t578 * t565 + t577 * t614 + t533;
t545 = -t558 - t630;
t574 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t603 + Ifges(4,2) * t605) * qJD(2);
t575 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t603 + Ifges(4,4) * t605) * qJD(2);
t640 = mrSges(4,1) * t545 - mrSges(4,2) * t546 + Ifges(4,5) * t587 + Ifges(4,6) * t588 + Ifges(4,3) * qJDD(3) - pkin(3) * t531 + qJ(4) * t526 + (t574 * t603 - t575 * t605) * qJD(2) + t632 * t523 + t601 * t524;
t586 = (-mrSges(4,1) * t605 + mrSges(4,2) * t603) * qJD(2);
t591 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t623;
t521 = m(4) * t546 - qJDD(3) * mrSges(4,2) + t588 * mrSges(4,3) - qJD(3) * t591 + t586 * t622 + t526;
t592 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t622;
t530 = m(4) * t545 + qJDD(3) * mrSges(4,1) - t587 * mrSges(4,3) + qJD(3) * t592 - t586 * t623 - t531;
t517 = t605 * t521 - t603 * t530;
t513 = m(3) * t570 - t608 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t517;
t525 = t601 * t528 + t632 * t529;
t522 = -m(4) * t560 + t588 * mrSges(4,1) - t587 * mrSges(4,2) - t591 * t623 + t592 * t622 - t525;
t519 = m(3) * t569 + qJDD(2) * mrSges(3,1) - t608 * mrSges(3,2) + t522;
t616 = t606 * t513 - t604 * t519;
t508 = m(2) * t590 + t616;
t516 = t603 * t521 + t605 * t530;
t515 = (m(2) + m(3)) * t589 - t516;
t629 = t602 * t508 + t633 * t515;
t509 = t604 * t513 + t606 * t519;
t620 = m(2) * t600 + t509;
t617 = t633 * t508 - t602 * t515;
t573 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t603 + Ifges(4,6) * t605) * qJD(2);
t505 = mrSges(4,2) * t560 - mrSges(4,3) * t545 + Ifges(4,1) * t587 + Ifges(4,4) * t588 + Ifges(4,5) * qJDD(3) - qJ(4) * t525 - qJD(3) * t574 - t601 * t523 + t632 * t524 + t573 * t622;
t532 = t568 * mrSges(6,2) + t578 * t556 - t612;
t510 = Ifges(4,4) * t587 + Ifges(4,6) * qJDD(3) - t573 * t623 + qJD(3) * t575 - mrSges(4,1) * t560 + mrSges(4,3) * t546 - mrSges(5,1) * t538 + mrSges(5,2) * t539 + mrSges(6,1) * t536 - mrSges(6,3) * t535 + pkin(4) * t532 - qJ(5) * (-t566 * t622 + t628) - pkin(3) * t525 + t626 * t578 + (qJ(5) * t556 + t625) * t577 - t644 * t568 + (qJ(5) * mrSges(6,2) + t641) * t567 + (Ifges(4,2) + t642) * t588;
t610 = mrSges(3,1) * t569 - mrSges(3,2) * t570 + Ifges(3,3) * qJDD(2) + pkin(2) * t522 + pkin(6) * t517 + t603 * t505 + t605 * t510;
t504 = mrSges(3,1) * t589 + mrSges(3,3) * t570 + t608 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t516 - t640;
t503 = -mrSges(3,2) * t589 - mrSges(3,3) * t569 + Ifges(3,5) * qJDD(2) - t608 * Ifges(3,6) - pkin(6) * t516 + t605 * t505 - t603 * t510;
t502 = -mrSges(2,1) * t600 + mrSges(2,3) * t590 - pkin(1) * t509 - t610;
t501 = mrSges(2,2) * t600 - mrSges(2,3) * t589 - pkin(5) * t509 + t606 * t503 - t604 * t504;
t1 = [-m(1) * g(1) + t617; -m(1) * g(2) + t629; -m(1) * g(3) + t620; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t629 + t633 * t501 - t602 * t502; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t617 + t602 * t501 + t633 * t502; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t589 - mrSges(2,2) * t590 + t604 * t503 + t606 * t504 + pkin(1) * (m(3) * t589 - t516) + pkin(5) * t616; t620; t610; t640; t531; t532;];
tauJB = t1;
