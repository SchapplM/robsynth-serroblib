% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 01:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:16:01
% EndTime: 2019-05-05 01:16:10
% DurationCPUTime: 7.20s
% Computational Cost: add. (97773->294), mult. (184968->372), div. (0->0), fcn. (125839->12), ass. (0->126)
t593 = sin(pkin(11));
t595 = cos(pkin(11));
t577 = g(1) * t593 - g(2) * t595;
t578 = -g(1) * t595 - g(2) * t593;
t590 = -g(3) + qJDD(1);
t604 = cos(qJ(2));
t596 = cos(pkin(6));
t600 = sin(qJ(2));
t626 = t596 * t600;
t594 = sin(pkin(6));
t627 = t594 * t600;
t550 = t577 * t626 + t604 * t578 + t590 * t627;
t615 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t550;
t598 = sin(qJ(5));
t599 = sin(qJ(4));
t602 = cos(qJ(5));
t603 = cos(qJ(4));
t566 = (t598 * t603 + t599 * t602) * qJD(2);
t549 = -t600 * t578 + (t577 * t596 + t590 * t594) * t604;
t633 = -pkin(2) - pkin(8);
t605 = qJD(2) ^ 2;
t632 = pkin(4) * t605;
t631 = mrSges(3,1) - mrSges(4,2);
t630 = (Ifges(3,5) - Ifges(4,4));
t629 = -Ifges(3,6) + Ifges(4,5);
t608 = -t605 * qJ(3) + qJDD(3) - t549;
t537 = t633 * qJDD(2) + t608;
t534 = t603 * t537;
t560 = -t577 * t594 + t590 * t596;
t622 = qJD(2) * qJD(4);
t576 = qJDD(2) * t603 - t599 * t622;
t519 = (qJDD(4) * pkin(4)) - t576 * pkin(9) + t534 + (-pkin(9) * t622 - t603 * t632 - t560) * t599;
t527 = t599 * t537 + t603 * t560;
t575 = -qJDD(2) * t599 - t603 * t622;
t623 = qJD(2) * t603;
t582 = (qJD(4) * pkin(4)) - pkin(9) * t623;
t589 = t599 ^ 2;
t520 = pkin(9) * t575 - qJD(4) * t582 - t589 * t632 + t527;
t515 = t598 * t519 + t602 * t520;
t567 = (-t598 * t599 + t602 * t603) * qJD(2);
t541 = -qJD(5) * t567 + t575 * t602 - t576 * t598;
t552 = mrSges(6,1) * t566 + mrSges(6,2) * t567;
t586 = qJD(4) + qJD(5);
t558 = mrSges(6,1) * t586 - mrSges(6,3) * t567;
t585 = qJDD(4) + qJDD(5);
t553 = pkin(5) * t566 - pkin(10) * t567;
t584 = t586 ^ 2;
t513 = -pkin(5) * t584 + pkin(10) * t585 - t553 * t566 + t515;
t525 = -t575 * pkin(4) + t582 * t623 + (-pkin(9) * t589 + t633) * t605 + t615;
t542 = -qJD(5) * t566 + t575 * t598 + t576 * t602;
t516 = (t566 * t586 - t542) * pkin(10) + (t567 * t586 - t541) * pkin(5) + t525;
t597 = sin(qJ(6));
t601 = cos(qJ(6));
t510 = -t513 * t597 + t516 * t601;
t554 = -t567 * t597 + t586 * t601;
t523 = qJD(6) * t554 + t542 * t601 + t585 * t597;
t555 = t567 * t601 + t586 * t597;
t532 = -mrSges(7,1) * t554 + mrSges(7,2) * t555;
t539 = qJDD(6) - t541;
t561 = qJD(6) + t566;
t544 = -mrSges(7,2) * t561 + mrSges(7,3) * t554;
t508 = m(7) * t510 + mrSges(7,1) * t539 - mrSges(7,3) * t523 - t532 * t555 + t544 * t561;
t511 = t513 * t601 + t516 * t597;
t522 = -qJD(6) * t555 - t542 * t597 + t585 * t601;
t545 = mrSges(7,1) * t561 - mrSges(7,3) * t555;
t509 = m(7) * t511 - mrSges(7,2) * t539 + mrSges(7,3) * t522 + t532 * t554 - t545 * t561;
t617 = -t508 * t597 + t601 * t509;
t499 = m(6) * t515 - mrSges(6,2) * t585 + mrSges(6,3) * t541 - t552 * t566 - t558 * t586 + t617;
t514 = t519 * t602 - t520 * t598;
t557 = -mrSges(6,2) * t586 - mrSges(6,3) * t566;
t512 = -pkin(5) * t585 - pkin(10) * t584 + t553 * t567 - t514;
t610 = -m(7) * t512 + t522 * mrSges(7,1) - mrSges(7,2) * t523 + t554 * t544 - t545 * t555;
t504 = m(6) * t514 + mrSges(6,1) * t585 - mrSges(6,3) * t542 - t552 * t567 + t557 * t586 + t610;
t492 = t598 * t499 + t602 * t504;
t526 = -t599 * t560 + t534;
t574 = (mrSges(5,1) * t599 + mrSges(5,2) * t603) * qJD(2);
t624 = qJD(2) * t599;
t579 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t624;
t490 = m(5) * t526 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t576 + qJD(4) * t579 - t574 * t623 + t492;
t580 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t623;
t618 = t602 * t499 - t504 * t598;
t491 = m(5) * t527 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t575 - qJD(4) * t580 - t574 * t624 + t618;
t486 = t603 * t490 + t599 * t491;
t543 = -qJDD(2) * pkin(2) + t608;
t611 = -m(4) * t543 + (t605 * mrSges(4,3)) - t486;
t482 = m(3) * t549 - (t605 * mrSges(3,2)) + t631 * qJDD(2) + t611;
t628 = t482 * t604;
t619 = -t490 * t599 + t603 * t491;
t485 = m(4) * t560 + t619;
t484 = m(3) * t560 + t485;
t540 = t605 * pkin(2) - t615;
t536 = t633 * t605 + t615;
t500 = t601 * t508 + t597 * t509;
t609 = m(6) * t525 - mrSges(6,1) * t541 + t542 * mrSges(6,2) + t557 * t566 + t567 * t558 + t500;
t607 = -m(5) * t536 + mrSges(5,1) * t575 - t576 * mrSges(5,2) - t579 * t624 - t580 * t623 - t609;
t606 = -m(4) * t540 + (t605 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t607;
t496 = m(3) * t550 - (mrSges(3,1) * t605) - qJDD(2) * mrSges(3,2) + t606;
t473 = -t484 * t594 + t496 * t626 + t596 * t628;
t471 = m(2) * t577 + t473;
t479 = -t482 * t600 + t604 * t496;
t478 = m(2) * t578 + t479;
t625 = t595 * t471 + t593 * t478;
t472 = t596 * t484 + t496 * t627 + t594 * t628;
t620 = -t471 * t593 + t595 * t478;
t528 = Ifges(7,5) * t555 + Ifges(7,6) * t554 + Ifges(7,3) * t561;
t530 = Ifges(7,1) * t555 + Ifges(7,4) * t554 + Ifges(7,5) * t561;
t501 = -mrSges(7,1) * t512 + mrSges(7,3) * t511 + Ifges(7,4) * t523 + Ifges(7,2) * t522 + Ifges(7,6) * t539 - t528 * t555 + t530 * t561;
t529 = Ifges(7,4) * t555 + Ifges(7,2) * t554 + Ifges(7,6) * t561;
t502 = mrSges(7,2) * t512 - mrSges(7,3) * t510 + Ifges(7,1) * t523 + Ifges(7,4) * t522 + Ifges(7,5) * t539 + t528 * t554 - t529 * t561;
t546 = Ifges(6,5) * t567 - Ifges(6,6) * t566 + Ifges(6,3) * t586;
t547 = Ifges(6,4) * t567 - Ifges(6,2) * t566 + Ifges(6,6) * t586;
t487 = mrSges(6,2) * t525 - mrSges(6,3) * t514 + Ifges(6,1) * t542 + Ifges(6,4) * t541 + Ifges(6,5) * t585 - pkin(10) * t500 - t501 * t597 + t502 * t601 - t546 * t566 - t547 * t586;
t548 = Ifges(6,1) * t567 - Ifges(6,4) * t566 + Ifges(6,5) * t586;
t488 = -mrSges(6,1) * t525 - mrSges(7,1) * t510 + mrSges(7,2) * t511 + mrSges(6,3) * t515 + Ifges(6,4) * t542 - Ifges(7,5) * t523 + Ifges(6,2) * t541 + Ifges(6,6) * t585 - Ifges(7,6) * t522 - Ifges(7,3) * t539 - pkin(5) * t500 - t529 * t555 + t530 * t554 - t546 * t567 + t548 * t586;
t563 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t603 - Ifges(5,6) * t599) * qJD(2);
t565 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t603 - Ifges(5,4) * t599) * qJD(2);
t474 = -mrSges(5,1) * t536 + mrSges(5,3) * t527 + Ifges(5,4) * t576 + Ifges(5,2) * t575 + Ifges(5,6) * qJDD(4) - pkin(4) * t609 + pkin(9) * t618 + qJD(4) * t565 + t598 * t487 + t602 * t488 - t563 * t623;
t564 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t603 - Ifges(5,2) * t599) * qJD(2);
t475 = mrSges(5,2) * t536 - mrSges(5,3) * t526 + Ifges(5,1) * t576 + Ifges(5,4) * t575 + Ifges(5,5) * qJDD(4) - pkin(9) * t492 - qJD(4) * t564 + t487 * t602 - t488 * t598 - t563 * t624;
t468 = -mrSges(4,1) * t540 + mrSges(3,3) * t550 - pkin(2) * t485 - pkin(3) * t607 - pkin(8) * t619 - t629 * qJDD(2) - t603 * t474 - t599 * t475 - t631 * t560 + (t630 * t605);
t469 = t629 * t605 + t630 * qJDD(2) + (Ifges(5,3) * qJDD(4)) + (mrSges(3,2) - mrSges(4,3)) * t560 + (t564 * t603 + t565 * t599) * qJD(2) + pkin(5) * t610 + t601 * t501 + t597 * t502 + Ifges(6,3) * t585 + t567 * t547 + Ifges(5,6) * t575 + Ifges(5,5) * t576 + t566 * t548 - mrSges(3,3) * t549 + Ifges(6,6) * t541 + Ifges(6,5) * t542 + mrSges(4,1) * t543 + mrSges(5,1) * t526 - mrSges(5,2) * t527 + mrSges(6,1) * t514 - mrSges(6,2) * t515 + pkin(4) * t492 - qJ(3) * t485 + pkin(3) * t486 + pkin(10) * t617;
t612 = pkin(7) * t479 + t468 * t604 + t469 * t600;
t467 = mrSges(3,1) * t549 - mrSges(3,2) * t550 + mrSges(4,2) * t543 - mrSges(4,3) * t540 + t603 * t475 - t599 * t474 - pkin(8) * t486 + pkin(2) * t611 + qJ(3) * t606 + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(3,3)) * qJDD(2);
t466 = mrSges(2,2) * t590 - mrSges(2,3) * t577 - t600 * t468 + t604 * t469 + (-t472 * t594 - t473 * t596) * pkin(7);
t465 = -mrSges(2,1) * t590 + mrSges(2,3) * t578 - pkin(1) * t472 - t594 * t467 + t612 * t596;
t1 = [-m(1) * g(1) + t620; -m(1) * g(2) + t625; -m(1) * g(3) + m(2) * t590 + t472; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t625 - t593 * t465 + t595 * t466; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t620 + t595 * t465 + t593 * t466; -mrSges(1,1) * g(2) + mrSges(2,1) * t577 + mrSges(1,2) * g(1) - mrSges(2,2) * t578 + pkin(1) * t473 + t596 * t467 + t612 * t594;];
tauB  = t1;
