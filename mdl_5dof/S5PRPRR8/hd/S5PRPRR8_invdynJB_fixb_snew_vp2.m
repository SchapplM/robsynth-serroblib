% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:44
% EndTime: 2019-12-05 16:02:47
% DurationCPUTime: 2.09s
% Computational Cost: add. (23293->219), mult. (42064->276), div. (0->0), fcn. (25965->10), ass. (0->103)
t584 = sin(pkin(9));
t586 = cos(pkin(9));
t570 = g(1) * t584 - g(2) * t586;
t571 = -g(1) * t586 - g(2) * t584;
t581 = -g(3) + qJDD(1);
t593 = cos(qJ(2));
t587 = cos(pkin(5));
t590 = sin(qJ(2));
t617 = t587 * t590;
t585 = sin(pkin(5));
t618 = t585 * t590;
t536 = t570 * t617 + t593 * t571 + t581 * t618;
t625 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t536;
t535 = -t590 * t571 + (t570 * t587 + t581 * t585) * t593;
t624 = -pkin(2) - pkin(7);
t623 = mrSges(3,1) - mrSges(4,2);
t622 = (-Ifges(4,4) + Ifges(3,5));
t621 = Ifges(4,5) - Ifges(3,6);
t595 = qJD(2) ^ 2;
t599 = -qJ(3) * t595 + qJDD(3) - t535;
t532 = t624 * qJDD(2) + t599;
t549 = -t570 * t585 + t581 * t587;
t589 = sin(qJ(4));
t592 = cos(qJ(4));
t528 = t589 * t532 + t592 * t549;
t566 = (mrSges(5,1) * t589 + mrSges(5,2) * t592) * qJD(2);
t613 = qJD(2) * qJD(4);
t610 = t592 * t613;
t568 = -qJDD(2) * t589 - t610;
t614 = qJD(2) * t592;
t573 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t614;
t567 = (pkin(4) * t589 - pkin(8) * t592) * qJD(2);
t594 = qJD(4) ^ 2;
t615 = qJD(2) * t589;
t524 = -pkin(4) * t594 + qJDD(4) * pkin(8) - t567 * t615 + t528;
t531 = t624 * t595 - t625;
t611 = t589 * t613;
t569 = qJDD(2) * t592 - t611;
t525 = (-t569 + t611) * pkin(8) + (-t568 + t610) * pkin(4) + t531;
t588 = sin(qJ(5));
t591 = cos(qJ(5));
t521 = -t524 * t588 + t525 * t591;
t564 = qJD(4) * t591 - t588 * t614;
t543 = qJD(5) * t564 + qJDD(4) * t588 + t569 * t591;
t565 = qJD(4) * t588 + t591 * t614;
t544 = -mrSges(6,1) * t564 + mrSges(6,2) * t565;
t575 = qJD(5) + t615;
t546 = -mrSges(6,2) * t575 + mrSges(6,3) * t564;
t561 = qJDD(5) - t568;
t518 = m(6) * t521 + t561 * mrSges(6,1) - t543 * mrSges(6,3) - t544 * t565 + t546 * t575;
t522 = t524 * t591 + t525 * t588;
t542 = -qJD(5) * t565 + qJDD(4) * t591 - t569 * t588;
t547 = mrSges(6,1) * t575 - mrSges(6,3) * t565;
t519 = m(6) * t522 - t561 * mrSges(6,2) + t542 * mrSges(6,3) + t544 * t564 - t547 * t575;
t607 = -t518 * t588 + t591 * t519;
t507 = m(5) * t528 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t568 - qJD(4) * t573 - t566 * t615 + t607;
t619 = t549 * t589;
t527 = t532 * t592 - t619;
t572 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t615;
t523 = -qJDD(4) * pkin(4) - pkin(8) * t594 + t619 + (qJD(2) * t567 - t532) * t592;
t600 = -m(6) * t523 + t542 * mrSges(6,1) - t543 * mrSges(6,2) + t564 * t546 - t547 * t565;
t514 = m(5) * t527 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t569 + qJD(4) * t572 - t566 * t614 + t600;
t501 = t507 * t589 + t514 * t592;
t534 = -qJDD(2) * pkin(2) + t599;
t602 = -m(4) * t534 + (t595 * mrSges(4,3)) - t501;
t496 = m(3) * t535 - (mrSges(3,2) * t595) + t623 * qJDD(2) + t602;
t620 = t496 * t593;
t608 = t592 * t507 - t514 * t589;
t500 = m(4) * t549 + t608;
t499 = m(3) * t549 + t500;
t533 = pkin(2) * t595 + t625;
t509 = t591 * t518 + t588 * t519;
t601 = -m(5) * t531 + t568 * mrSges(5,1) - t569 * mrSges(5,2) - t572 * t615 - t573 * t614 - t509;
t597 = -m(4) * t533 + (t595 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t601;
t505 = m(3) * t536 - (mrSges(3,1) * t595) - qJDD(2) * mrSges(3,2) + t597;
t487 = -t499 * t585 + t505 * t617 + t587 * t620;
t485 = m(2) * t570 + t487;
t491 = -t496 * t590 + t593 * t505;
t490 = m(2) * t571 + t491;
t616 = t586 * t485 + t584 * t490;
t486 = t587 * t499 + t505 * t618 + t585 * t620;
t609 = -t485 * t584 + t586 * t490;
t605 = m(2) * t581 + t486;
t537 = Ifges(6,5) * t565 + Ifges(6,6) * t564 + Ifges(6,3) * t575;
t539 = Ifges(6,1) * t565 + Ifges(6,4) * t564 + Ifges(6,5) * t575;
t512 = -mrSges(6,1) * t523 + mrSges(6,3) * t522 + Ifges(6,4) * t543 + Ifges(6,2) * t542 + Ifges(6,6) * t561 - t537 * t565 + t539 * t575;
t538 = Ifges(6,4) * t565 + Ifges(6,2) * t564 + Ifges(6,6) * t575;
t513 = mrSges(6,2) * t523 - mrSges(6,3) * t521 + Ifges(6,1) * t543 + Ifges(6,4) * t542 + Ifges(6,5) * t561 + t537 * t564 - t538 * t575;
t553 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t592 - Ifges(5,6) * t589) * qJD(2);
t554 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t592 - Ifges(5,2) * t589) * qJD(2);
t492 = mrSges(5,2) * t531 - mrSges(5,3) * t527 + Ifges(5,1) * t569 + Ifges(5,4) * t568 + Ifges(5,5) * qJDD(4) - pkin(8) * t509 - qJD(4) * t554 - t512 * t588 + t513 * t591 - t553 * t615;
t555 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t592 - Ifges(5,4) * t589) * qJD(2);
t596 = mrSges(6,1) * t521 - mrSges(6,2) * t522 + Ifges(6,5) * t543 + Ifges(6,6) * t542 + Ifges(6,3) * t561 + t538 * t565 - t539 * t564;
t493 = -mrSges(5,1) * t531 + mrSges(5,3) * t528 + Ifges(5,4) * t569 + Ifges(5,2) * t568 + Ifges(5,6) * qJDD(4) - pkin(4) * t509 + qJD(4) * t555 - t553 * t614 - t596;
t482 = -mrSges(4,1) * t533 + mrSges(3,3) * t536 - pkin(2) * t500 - pkin(3) * t601 - pkin(7) * t608 - t621 * qJDD(2) - t589 * t492 - t592 * t493 - t623 * t549 + (t622 * t595);
t598 = mrSges(5,1) * t527 - mrSges(5,2) * t528 + Ifges(5,5) * t569 + Ifges(5,6) * t568 + Ifges(5,3) * qJDD(4) + pkin(4) * t600 + pkin(8) * t607 + t591 * t512 + t588 * t513 + t554 * t614 + t555 * t615;
t483 = t621 * t595 + (mrSges(3,2) - mrSges(4,3)) * t549 + mrSges(4,1) * t534 - mrSges(3,3) * t535 + pkin(3) * t501 - qJ(3) * t500 + t598 + t622 * qJDD(2);
t603 = pkin(6) * t491 + t482 * t593 + t483 * t590;
t497 = qJDD(2) * mrSges(4,2) - t602;
t481 = mrSges(3,1) * t535 - mrSges(3,2) * t536 + mrSges(4,2) * t534 - mrSges(4,3) * t533 + t592 * t492 - t589 * t493 - pkin(7) * t501 - pkin(2) * t497 + qJ(3) * t597 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t480 = mrSges(2,2) * t581 - mrSges(2,3) * t570 - t482 * t590 + t483 * t593 + (-t486 * t585 - t487 * t587) * pkin(6);
t479 = -mrSges(2,1) * t581 + mrSges(2,3) * t571 - pkin(1) * t486 - t481 * t585 + t603 * t587;
t1 = [-m(1) * g(1) + t609; -m(1) * g(2) + t616; -m(1) * g(3) + t605; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t616 - t584 * t479 + t586 * t480; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t609 + t586 * t479 + t584 * t480; -mrSges(1,1) * g(2) + mrSges(2,1) * t570 + mrSges(1,2) * g(1) - mrSges(2,2) * t571 + pkin(1) * t487 + t481 * t587 + t603 * t585; t605; t481; t497; t598; t596;];
tauJB = t1;
