% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:38
% EndTime: 2019-12-05 17:44:45
% DurationCPUTime: 6.56s
% Computational Cost: add. (59829->279), mult. (163866->381), div. (0->0), fcn. (112502->10), ass. (0->128)
t598 = sin(qJ(1));
t601 = cos(qJ(1));
t577 = t598 * g(2) - t601 * g(3);
t602 = qJD(1) ^ 2;
t642 = -t602 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t577;
t578 = t601 * g(2) + t598 * g(3);
t608 = -t602 * qJ(2) + qJDD(2) - t578;
t593 = sin(pkin(8));
t595 = cos(pkin(8));
t616 = -pkin(2) * t595 - qJ(3) * t593;
t633 = t593 * qJD(1);
t641 = (-pkin(1) + t616) * qJDD(1) + t608 - 0.2e1 * qJD(3) * t633;
t552 = -t595 * g(1) - t642 * t593;
t640 = mrSges(3,2) * t593;
t639 = Ifges(3,6) * t595;
t591 = t593 ^ 2;
t638 = t591 * t602;
t592 = sin(pkin(9));
t637 = t592 * t593;
t594 = cos(pkin(9));
t636 = t593 * t594;
t553 = -t593 * g(1) + t642 * t595;
t572 = (-mrSges(3,1) * t595 + t640) * qJD(1);
t571 = t616 * qJD(1);
t632 = t595 * qJD(1);
t542 = t571 * t632 + t553;
t613 = -pkin(3) * t595 - pkin(6) * t636;
t635 = t641 * t594;
t521 = t613 * qJDD(1) + (-t542 + (-pkin(3) * t591 * t594 + pkin(6) * t593 * t595) * t602) * t592 + t635;
t530 = t594 * t542 + t641 * t592;
t570 = t613 * qJD(1);
t630 = qJDD(1) * t593;
t626 = t592 * t630;
t628 = t592 ^ 2 * t638;
t522 = -pkin(3) * t628 - pkin(6) * t626 + t570 * t632 + t530;
t597 = sin(qJ(4));
t600 = cos(qJ(4));
t510 = t600 * t521 - t597 * t522;
t610 = (-t592 * t600 - t594 * t597) * t593;
t560 = qJD(1) * t610;
t609 = (-t592 * t597 + t594 * t600) * t593;
t546 = t560 * qJD(4) + qJDD(1) * t609;
t561 = qJD(1) * t609;
t629 = t595 * qJDD(1);
t581 = qJDD(4) - t629;
t582 = qJD(4) - t632;
t508 = (t560 * t582 - t546) * pkin(7) + (t560 * t561 + t581) * pkin(4) + t510;
t511 = t597 * t521 + t600 * t522;
t545 = -t561 * qJD(4) + qJDD(1) * t610;
t551 = t582 * pkin(4) - t561 * pkin(7);
t559 = t560 ^ 2;
t509 = -t559 * pkin(4) + t545 * pkin(7) - t582 * t551 + t511;
t596 = sin(qJ(5));
t599 = cos(qJ(5));
t506 = t599 * t508 - t596 * t509;
t539 = t599 * t560 - t596 * t561;
t517 = t539 * qJD(5) + t596 * t545 + t599 * t546;
t540 = t596 * t560 + t599 * t561;
t528 = -t539 * mrSges(6,1) + t540 * mrSges(6,2);
t580 = qJD(5) + t582;
t532 = -t580 * mrSges(6,2) + t539 * mrSges(6,3);
t576 = qJDD(5) + t581;
t504 = m(6) * t506 + t576 * mrSges(6,1) - t517 * mrSges(6,3) - t540 * t528 + t580 * t532;
t507 = t596 * t508 + t599 * t509;
t516 = -t540 * qJD(5) + t599 * t545 - t596 * t546;
t533 = t580 * mrSges(6,1) - t540 * mrSges(6,3);
t505 = m(6) * t507 - t576 * mrSges(6,2) + t516 * mrSges(6,3) + t539 * t528 - t580 * t533;
t496 = t599 * t504 + t596 * t505;
t543 = -t560 * mrSges(5,1) + t561 * mrSges(5,2);
t547 = -t582 * mrSges(5,2) + t560 * mrSges(5,3);
t494 = m(5) * t510 + t581 * mrSges(5,1) - t546 * mrSges(5,3) - t561 * t543 + t582 * t547 + t496;
t548 = t582 * mrSges(5,1) - t561 * mrSges(5,3);
t621 = -t596 * t504 + t599 * t505;
t495 = m(5) * t511 - t581 * mrSges(5,2) + t545 * mrSges(5,3) + t560 * t543 - t582 * t548 + t621;
t490 = t600 * t494 + t597 * t495;
t529 = -t592 * t542 + t635;
t619 = mrSges(4,1) * t592 + mrSges(4,2) * t594;
t564 = t619 * t633;
t611 = mrSges(4,2) * t595 - mrSges(4,3) * t637;
t567 = t611 * qJD(1);
t612 = -mrSges(4,1) * t595 - mrSges(4,3) * t636;
t488 = m(4) * t529 + t612 * qJDD(1) + (-t564 * t636 - t567 * t595) * qJD(1) + t490;
t568 = t612 * qJD(1);
t622 = -t597 * t494 + t600 * t495;
t489 = m(4) * t530 + t611 * qJDD(1) + (-t564 * t637 + t568 * t595) * qJD(1) + t622;
t623 = -t592 * t488 + t594 * t489;
t483 = m(3) * t553 + (qJDD(1) * mrSges(3,3) + qJD(1) * t572) * t595 + t623;
t541 = t571 * t633 + qJDD(3) - t552;
t531 = t594 * t570 * t633 + pkin(3) * t626 - pkin(6) * t628 + t541;
t513 = -t545 * pkin(4) - t559 * pkin(7) + t561 * t551 + t531;
t615 = m(6) * t513 - t516 * mrSges(6,1) + t517 * mrSges(6,2) - t539 * t532 + t540 * t533;
t604 = m(5) * t531 - t545 * mrSges(5,1) + t546 * mrSges(5,2) - t560 * t547 + t561 * t548 + t615;
t603 = -m(4) * t541 - t604;
t500 = t603 + ((-mrSges(3,3) - t619) * qJDD(1) + (-t567 * t592 - t568 * t594 - t572) * qJD(1)) * t593 + m(3) * t552;
t479 = t593 * t483 + t595 * t500;
t624 = t595 * t483 - t593 * t500;
t478 = m(2) * t577 - t602 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t624;
t484 = t594 * t488 + t592 * t489;
t566 = -qJDD(1) * pkin(1) + t608;
t605 = -m(3) * t566 + mrSges(3,1) * t629 - t484 + (t595 ^ 2 * t602 + t638) * mrSges(3,3);
t480 = m(2) * t578 - t602 * mrSges(2,2) + (mrSges(2,1) - t640) * qJDD(1) + t605;
t625 = t601 * t478 - t598 * t480;
t618 = Ifges(3,1) * t593 + Ifges(3,4) * t595;
t617 = Ifges(4,5) * t594 - Ifges(4,6) * t592;
t614 = -t598 * t478 - t601 * t480;
t607 = -Ifges(4,5) * t595 + (Ifges(4,1) * t594 - Ifges(4,4) * t592) * t593;
t606 = -Ifges(4,6) * t595 + (Ifges(4,4) * t594 - Ifges(4,2) * t592) * t593;
t573 = (Ifges(3,5) * t593 + t639) * qJD(1);
t557 = t607 * qJD(1);
t556 = t606 * qJD(1);
t555 = (-Ifges(4,3) * t595 + t617 * t593) * qJD(1);
t536 = Ifges(5,1) * t561 + Ifges(5,4) * t560 + Ifges(5,5) * t582;
t535 = Ifges(5,4) * t561 + Ifges(5,2) * t560 + Ifges(5,6) * t582;
t534 = Ifges(5,5) * t561 + Ifges(5,6) * t560 + Ifges(5,3) * t582;
t525 = Ifges(6,1) * t540 + Ifges(6,4) * t539 + Ifges(6,5) * t580;
t524 = Ifges(6,4) * t540 + Ifges(6,2) * t539 + Ifges(6,6) * t580;
t523 = Ifges(6,5) * t540 + Ifges(6,6) * t539 + Ifges(6,3) * t580;
t498 = mrSges(6,2) * t513 - mrSges(6,3) * t506 + Ifges(6,1) * t517 + Ifges(6,4) * t516 + Ifges(6,5) * t576 + t539 * t523 - t580 * t524;
t497 = -mrSges(6,1) * t513 + mrSges(6,3) * t507 + Ifges(6,4) * t517 + Ifges(6,2) * t516 + Ifges(6,6) * t576 - t540 * t523 + t580 * t525;
t486 = mrSges(5,2) * t531 - mrSges(5,3) * t510 + Ifges(5,1) * t546 + Ifges(5,4) * t545 + Ifges(5,5) * t581 - pkin(7) * t496 - t596 * t497 + t599 * t498 + t560 * t534 - t582 * t535;
t485 = -mrSges(5,1) * t531 + mrSges(5,3) * t511 + Ifges(5,4) * t546 + Ifges(5,2) * t545 + Ifges(5,6) * t581 - pkin(4) * t615 + pkin(7) * t621 + t599 * t497 + t596 * t498 - t561 * t534 + t582 * t536;
t476 = mrSges(4,2) * t541 - mrSges(4,3) * t529 - pkin(6) * t490 - t597 * t485 + t600 * t486 + (-t555 * t637 + t556 * t595) * qJD(1) + t607 * qJDD(1);
t475 = -mrSges(4,1) * t541 + mrSges(4,3) * t530 + t597 * t486 + t600 * t485 - pkin(3) * t604 + pkin(6) * t622 + (-t555 * t636 - t595 * t557) * qJD(1) + t606 * qJDD(1);
t474 = -Ifges(5,3) * t581 - Ifges(6,3) * t576 - t561 * t535 - mrSges(3,1) * t566 + mrSges(3,3) * t553 + t560 * t536 - t540 * t524 - Ifges(5,6) * t545 - Ifges(5,5) * t546 + t539 * t525 - mrSges(4,1) * t529 + mrSges(4,2) * t530 - Ifges(6,6) * t516 - Ifges(6,5) * t517 - mrSges(5,1) * t510 + mrSges(5,2) * t511 - mrSges(6,1) * t506 + mrSges(6,2) * t507 - pkin(4) * t496 - pkin(3) * t490 - pkin(2) * t484 + (Ifges(3,2) + Ifges(4,3)) * t629 + ((Ifges(3,4) - t617) * qJDD(1) + (-t556 * t594 - t557 * t592 - t573) * qJD(1)) * t593;
t473 = mrSges(3,2) * t566 - mrSges(3,3) * t552 - qJ(3) * t484 + t618 * qJDD(1) - t592 * t475 + t594 * t476 + t573 * t632;
t472 = t602 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t577 - mrSges(3,1) * t552 + mrSges(3,2) * t553 - t592 * t476 - t594 * t475 - pkin(2) * t603 - qJ(3) * t623 - pkin(1) * t479 + (-t639 + Ifges(2,6) + (pkin(2) * t619 - Ifges(3,5)) * t593) * qJDD(1) + (-pkin(2) * (-t567 * t637 - t568 * t636) + (-t593 * (Ifges(3,4) * t593 + Ifges(3,2) * t595) + t595 * t618) * qJD(1)) * qJD(1);
t471 = -mrSges(2,2) * g(1) - mrSges(2,3) * t578 + Ifges(2,5) * qJDD(1) - t602 * Ifges(2,6) - qJ(2) * t479 + t595 * t473 - t593 * t474;
t1 = [(-m(1) - m(2)) * g(1) + t479; -m(1) * g(2) + t614; -m(1) * g(3) + t625; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t578 - mrSges(2,2) * t577 + t593 * t473 + t595 * t474 + pkin(1) * (-mrSges(3,2) * t630 + t605) + qJ(2) * t624; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t625 - t598 * t471 - t601 * t472; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t614 + t601 * t471 - t598 * t472;];
tauB = t1;
