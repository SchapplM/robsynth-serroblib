% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:08
% EndTime: 2019-12-31 19:10:16
% DurationCPUTime: 6.83s
% Computational Cost: add. (74430->292), mult. (175630->368), div. (0->0), fcn. (128560->10), ass. (0->125)
t577 = qJD(1) ^ 2;
t567 = cos(pkin(9));
t603 = pkin(2) * t567;
t566 = sin(pkin(9));
t602 = mrSges(3,2) * t566;
t565 = t567 ^ 2;
t601 = t565 * t577;
t571 = sin(qJ(1));
t575 = cos(qJ(1));
t555 = -t575 * g(1) - t571 * g(2);
t551 = -t577 * pkin(1) + qJDD(1) * qJ(2) + t555;
t596 = qJD(1) * qJD(2);
t593 = -t567 * g(3) - 0.2e1 * t566 * t596;
t523 = (-pkin(6) * qJDD(1) + t577 * t603 - t551) * t566 + t593;
t538 = -t566 * g(3) + (t551 + 0.2e1 * t596) * t567;
t594 = qJDD(1) * t567;
t524 = -pkin(2) * t601 + pkin(6) * t594 + t538;
t570 = sin(qJ(3));
t574 = cos(qJ(3));
t503 = t570 * t523 + t574 * t524;
t597 = t567 * qJD(1);
t598 = t566 * qJD(1);
t549 = -t570 * t598 + t574 * t597;
t583 = t566 * t574 + t567 * t570;
t550 = t583 * qJD(1);
t530 = -t549 * mrSges(4,1) + t550 * mrSges(4,2);
t546 = t550 * qJD(3);
t595 = qJDD(1) * t566;
t535 = -t570 * t595 + t574 * t594 - t546;
t543 = qJD(3) * mrSges(4,1) - t550 * mrSges(4,3);
t533 = -t549 * pkin(3) - t550 * pkin(7);
t576 = qJD(3) ^ 2;
t493 = -t576 * pkin(3) + qJDD(3) * pkin(7) + t549 * t533 + t503;
t564 = t566 ^ 2;
t554 = t571 * g(1) - t575 * g(2);
t587 = qJDD(2) - t554;
t534 = (-pkin(1) - t603) * qJDD(1) + (-qJ(2) + (-t564 - t565) * pkin(6)) * t577 + t587;
t599 = t549 * qJD(3);
t536 = t583 * qJDD(1) + t599;
t496 = (-t536 - t599) * pkin(7) + (-t535 + t546) * pkin(3) + t534;
t569 = sin(qJ(4));
t573 = cos(qJ(4));
t486 = -t569 * t493 + t573 * t496;
t540 = t573 * qJD(3) - t569 * t550;
t512 = t540 * qJD(4) + t569 * qJDD(3) + t573 * t536;
t532 = qJDD(4) - t535;
t541 = t569 * qJD(3) + t573 * t550;
t547 = qJD(4) - t549;
t483 = (t540 * t547 - t512) * pkin(8) + (t540 * t541 + t532) * pkin(4) + t486;
t487 = t573 * t493 + t569 * t496;
t511 = -t541 * qJD(4) + t573 * qJDD(3) - t569 * t536;
t522 = t547 * pkin(4) - t541 * pkin(8);
t539 = t540 ^ 2;
t484 = -t539 * pkin(4) + t511 * pkin(8) - t547 * t522 + t487;
t568 = sin(qJ(5));
t572 = cos(qJ(5));
t481 = t572 * t483 - t568 * t484;
t513 = t572 * t540 - t568 * t541;
t490 = t513 * qJD(5) + t568 * t511 + t572 * t512;
t514 = t568 * t540 + t572 * t541;
t501 = -t513 * mrSges(6,1) + t514 * mrSges(6,2);
t545 = qJD(5) + t547;
t504 = -t545 * mrSges(6,2) + t513 * mrSges(6,3);
t529 = qJDD(5) + t532;
t479 = m(6) * t481 + t529 * mrSges(6,1) - t490 * mrSges(6,3) - t514 * t501 + t545 * t504;
t482 = t568 * t483 + t572 * t484;
t489 = -t514 * qJD(5) + t572 * t511 - t568 * t512;
t505 = t545 * mrSges(6,1) - t514 * mrSges(6,3);
t480 = m(6) * t482 - t529 * mrSges(6,2) + t489 * mrSges(6,3) + t513 * t501 - t545 * t505;
t471 = t572 * t479 + t568 * t480;
t515 = -t540 * mrSges(5,1) + t541 * mrSges(5,2);
t518 = -t547 * mrSges(5,2) + t540 * mrSges(5,3);
t469 = m(5) * t486 + t532 * mrSges(5,1) - t512 * mrSges(5,3) - t541 * t515 + t547 * t518 + t471;
t519 = t547 * mrSges(5,1) - t541 * mrSges(5,3);
t588 = -t568 * t479 + t572 * t480;
t470 = m(5) * t487 - t532 * mrSges(5,2) + t511 * mrSges(5,3) + t540 * t515 - t547 * t519 + t588;
t589 = -t569 * t469 + t573 * t470;
t464 = m(4) * t503 - qJDD(3) * mrSges(4,2) + t535 * mrSges(4,3) - qJD(3) * t543 + t549 * t530 + t589;
t502 = t574 * t523 - t570 * t524;
t542 = -qJD(3) * mrSges(4,2) + t549 * mrSges(4,3);
t492 = -qJDD(3) * pkin(3) - t576 * pkin(7) + t550 * t533 - t502;
t485 = -t511 * pkin(4) - t539 * pkin(8) + t541 * t522 + t492;
t581 = m(6) * t485 - t489 * mrSges(6,1) + t490 * mrSges(6,2) - t513 * t504 + t514 * t505;
t578 = -m(5) * t492 + t511 * mrSges(5,1) - t512 * mrSges(5,2) + t540 * t518 - t541 * t519 - t581;
t475 = m(4) * t502 + qJDD(3) * mrSges(4,1) - t536 * mrSges(4,3) + qJD(3) * t542 - t550 * t530 + t578;
t458 = t570 * t464 + t574 * t475;
t537 = -t566 * t551 + t593;
t582 = mrSges(3,3) * qJDD(1) + t577 * (-mrSges(3,1) * t567 + t602);
t456 = m(3) * t537 - t582 * t566 + t458;
t590 = t574 * t464 - t570 * t475;
t457 = m(3) * t538 + t582 * t567 + t590;
t591 = -t566 * t456 + t567 * t457;
t449 = m(2) * t555 - t577 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t591;
t548 = -qJDD(1) * pkin(1) - t577 * qJ(2) + t587;
t465 = t573 * t469 + t569 * t470;
t580 = m(4) * t534 - t535 * mrSges(4,1) + t536 * mrSges(4,2) - t549 * t542 + t550 * t543 + t465;
t579 = -m(3) * t548 + mrSges(3,1) * t594 - t580 + (t564 * t577 + t601) * mrSges(3,3);
t461 = -t577 * mrSges(2,2) + t579 + m(2) * t554 + (mrSges(2,1) - t602) * qJDD(1);
t600 = t571 * t449 + t575 * t461;
t450 = t567 * t456 + t566 * t457;
t592 = t575 * t449 - t571 * t461;
t586 = Ifges(3,1) * t566 + Ifges(3,4) * t567;
t585 = Ifges(3,4) * t566 + Ifges(3,2) * t567;
t584 = Ifges(3,5) * t566 + Ifges(3,6) * t567;
t553 = t584 * qJD(1);
t527 = Ifges(4,1) * t550 + Ifges(4,4) * t549 + Ifges(4,5) * qJD(3);
t526 = Ifges(4,4) * t550 + Ifges(4,2) * t549 + Ifges(4,6) * qJD(3);
t525 = Ifges(4,5) * t550 + Ifges(4,6) * t549 + Ifges(4,3) * qJD(3);
t508 = Ifges(5,1) * t541 + Ifges(5,4) * t540 + Ifges(5,5) * t547;
t507 = Ifges(5,4) * t541 + Ifges(5,2) * t540 + Ifges(5,6) * t547;
t506 = Ifges(5,5) * t541 + Ifges(5,6) * t540 + Ifges(5,3) * t547;
t499 = Ifges(6,1) * t514 + Ifges(6,4) * t513 + Ifges(6,5) * t545;
t498 = Ifges(6,4) * t514 + Ifges(6,2) * t513 + Ifges(6,6) * t545;
t497 = Ifges(6,5) * t514 + Ifges(6,6) * t513 + Ifges(6,3) * t545;
t473 = mrSges(6,2) * t485 - mrSges(6,3) * t481 + Ifges(6,1) * t490 + Ifges(6,4) * t489 + Ifges(6,5) * t529 + t513 * t497 - t545 * t498;
t472 = -mrSges(6,1) * t485 + mrSges(6,3) * t482 + Ifges(6,4) * t490 + Ifges(6,2) * t489 + Ifges(6,6) * t529 - t514 * t497 + t545 * t499;
t459 = mrSges(5,2) * t492 - mrSges(5,3) * t486 + Ifges(5,1) * t512 + Ifges(5,4) * t511 + Ifges(5,5) * t532 - pkin(8) * t471 - t568 * t472 + t572 * t473 + t540 * t506 - t547 * t507;
t455 = -mrSges(5,1) * t492 + mrSges(5,3) * t487 + Ifges(5,4) * t512 + Ifges(5,2) * t511 + Ifges(5,6) * t532 - pkin(4) * t581 + pkin(8) * t588 + t572 * t472 + t568 * t473 - t541 * t506 + t547 * t508;
t451 = Ifges(4,4) * t536 + Ifges(4,2) * t535 + Ifges(4,6) * qJDD(3) - t550 * t525 + qJD(3) * t527 - mrSges(4,1) * t534 + mrSges(4,3) * t503 - Ifges(5,5) * t512 - Ifges(5,6) * t511 - Ifges(5,3) * t532 - t541 * t507 + t540 * t508 - mrSges(5,1) * t486 + mrSges(5,2) * t487 - Ifges(6,5) * t490 - Ifges(6,6) * t489 - Ifges(6,3) * t529 - t514 * t498 + t513 * t499 - mrSges(6,1) * t481 + mrSges(6,2) * t482 - pkin(4) * t471 - pkin(3) * t465;
t446 = mrSges(4,2) * t534 - mrSges(4,3) * t502 + Ifges(4,1) * t536 + Ifges(4,4) * t535 + Ifges(4,5) * qJDD(3) - pkin(7) * t465 - qJD(3) * t526 - t569 * t455 + t573 * t459 + t549 * t525;
t445 = mrSges(3,2) * t548 - mrSges(3,3) * t537 - pkin(6) * t458 + t586 * qJDD(1) + t574 * t446 - t570 * t451 + t553 * t597;
t444 = mrSges(2,1) * g(3) - pkin(1) * t450 + mrSges(2,3) * t555 - pkin(2) * t458 - mrSges(3,1) * t537 + mrSges(3,2) * t538 - t569 * t459 - t573 * t455 - pkin(3) * t578 - pkin(7) * t589 - mrSges(4,1) * t502 + mrSges(4,2) * t503 - Ifges(4,5) * t536 - Ifges(4,6) * t535 - Ifges(4,3) * qJDD(3) - t550 * t526 + t549 * t527 + (Ifges(2,6) - t584) * qJDD(1) + (-t566 * t585 + t567 * t586 + Ifges(2,5)) * t577;
t443 = -mrSges(3,1) * t548 + mrSges(3,3) * t538 - pkin(2) * t580 + pkin(6) * t590 + t585 * qJDD(1) + t570 * t446 + t574 * t451 - t553 * t598;
t442 = -mrSges(2,2) * g(3) - mrSges(2,3) * t554 + Ifges(2,5) * qJDD(1) - t577 * Ifges(2,6) - qJ(2) * t450 - t566 * t443 + t567 * t445;
t1 = [-m(1) * g(1) + t592; -m(1) * g(2) + t600; (-m(1) - m(2)) * g(3) + t450; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t600 + t575 * t442 - t571 * t444; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t592 + t571 * t442 + t575 * t444; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t554 - mrSges(2,2) * t555 + t566 * t445 + t567 * t443 + pkin(1) * (-mrSges(3,2) * t595 + t579) + qJ(2) * t591;];
tauB = t1;
