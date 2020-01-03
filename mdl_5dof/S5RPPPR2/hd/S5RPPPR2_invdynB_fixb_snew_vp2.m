% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:43
% EndTime: 2020-01-03 11:22:48
% DurationCPUTime: 4.82s
% Computational Cost: add. (36210->282), mult. (101411->386), div. (0->0), fcn. (68373->10), ass. (0->128)
t571 = sin(pkin(7));
t574 = cos(pkin(7));
t576 = sin(qJ(1));
t578 = cos(qJ(1));
t553 = -t576 * g(2) + t578 * g(3);
t579 = qJD(1) ^ 2;
t624 = -t579 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t553;
t525 = -t571 * g(1) + t624 * t574;
t590 = -pkin(2) * t574 - qJ(3) * t571;
t548 = t590 * qJD(1);
t609 = t574 * qJD(1);
t516 = t548 * t609 + t525;
t570 = sin(pkin(8));
t573 = cos(pkin(8));
t554 = -t578 * g(2) - t576 * g(3);
t586 = -t579 * qJ(2) + qJDD(2) - t554;
t610 = t571 * qJD(1);
t623 = (-pkin(1) + t590) * qJDD(1) + t586 - 0.2e1 * qJD(3) * t610;
t496 = -t570 * t516 + t623 * t573;
t524 = -t574 * g(1) - t624 * t571;
t569 = sin(pkin(9));
t572 = cos(pkin(9));
t614 = t571 * t573;
t585 = t569 * t614 + t572 * t574;
t537 = t585 * qJD(1);
t535 = t585 * qJDD(1);
t622 = 2 * qJD(4);
t621 = mrSges(3,2) * t571;
t620 = Ifges(4,4) * t573;
t619 = Ifges(3,6) * t574;
t618 = Ifges(4,6) * t574;
t617 = t571 ^ 2 * t579;
t616 = t574 ^ 2 * t579;
t615 = t570 * t571;
t613 = t574 * t579;
t549 = (-mrSges(3,1) * t574 + t621) * qJD(1);
t497 = t573 * t516 + t623 * t570;
t593 = mrSges(4,1) * t570 + mrSges(4,2) * t573;
t542 = t593 * t610;
t588 = -mrSges(4,1) * t574 - mrSges(4,3) * t614;
t546 = t588 * qJD(1);
t587 = mrSges(4,2) * t574 - mrSges(4,3) * t615;
t541 = (pkin(3) * t570 - qJ(4) * t573) * t610;
t604 = t570 * t610;
t606 = qJDD(1) * t574;
t494 = -pkin(3) * t616 - qJ(4) * t606 - t541 * t604 + t497;
t515 = t548 * t610 + qJDD(3) - t524;
t502 = ((-qJDD(1) * t573 - t570 * t613) * qJ(4) + (qJDD(1) * t570 - t573 * t613) * pkin(3)) * t571 + t515;
t490 = t572 * t494 + t569 * t502 - t537 * t622;
t603 = t573 * t610;
t538 = -t569 * t609 + t572 * t603;
t517 = t537 * mrSges(5,1) + t538 * mrSges(5,2);
t523 = mrSges(5,1) * t604 - t538 * mrSges(5,3);
t518 = t537 * pkin(4) - t538 * pkin(6);
t607 = qJDD(1) * t571;
t600 = t570 * t607;
t605 = t570 ^ 2 * t617;
t488 = -pkin(4) * t605 + pkin(6) * t600 - t537 * t518 + t490;
t493 = pkin(3) * t606 - qJ(4) * t616 + t541 * t603 + qJDD(4) - t496;
t536 = (-t569 * t574 + t572 * t614) * qJDD(1);
t491 = (t537 * t604 - t536) * pkin(6) + (t538 * t604 + t535) * pkin(4) + t493;
t575 = sin(qJ(5));
t577 = cos(qJ(5));
t485 = -t575 * t488 + t577 * t491;
t519 = -t575 * t538 + t577 * t604;
t520 = t577 * t538 + t575 * t604;
t504 = -t519 * mrSges(6,1) + t520 * mrSges(6,2);
t506 = t519 * qJD(5) + t577 * t536 + t575 * t600;
t534 = qJD(5) + t537;
t507 = -t534 * mrSges(6,2) + t519 * mrSges(6,3);
t533 = qJDD(5) + t535;
t483 = m(6) * t485 + t533 * mrSges(6,1) - t506 * mrSges(6,3) - t520 * t504 + t534 * t507;
t486 = t577 * t488 + t575 * t491;
t505 = -t520 * qJD(5) - t575 * t536 + t577 * t600;
t508 = t534 * mrSges(6,1) - t520 * mrSges(6,3);
t484 = m(6) * t486 - t533 * mrSges(6,2) + t505 * mrSges(6,3) + t519 * t504 - t534 * t508;
t595 = -t575 * t483 + t577 * t484;
t477 = m(5) * t490 - t535 * mrSges(5,3) - t537 * t517 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t523) * t615 + t595;
t589 = t569 * t494 - t572 * t502;
t489 = -0.2e1 * qJD(4) * t538 - t589;
t522 = -mrSges(5,2) * t604 - t537 * mrSges(5,3);
t487 = -pkin(4) * t600 - pkin(6) * t605 + (t622 + t518) * t538 + t589;
t583 = -m(6) * t487 + t505 * mrSges(6,1) - t506 * mrSges(6,2) + t519 * t507 - t520 * t508;
t481 = m(5) * t489 - t536 * mrSges(5,3) - t538 * t517 + (mrSges(5,1) * qJDD(1) + qJD(1) * t522) * t615 + t583;
t596 = t572 * t477 - t569 * t481;
t473 = m(4) * t497 + t587 * qJDD(1) + (-t542 * t615 + t546 * t574) * qJD(1) + t596;
t545 = t587 * qJD(1);
t478 = t577 * t483 + t575 * t484;
t580 = -m(5) * t493 - t535 * mrSges(5,1) - t536 * mrSges(5,2) - t537 * t522 - t538 * t523 - t478;
t475 = m(4) * t496 + t588 * qJDD(1) + (-t542 * t614 - t545 * t574) * qJD(1) + t580;
t597 = t573 * t473 - t570 * t475;
t466 = m(3) * t525 + (qJDD(1) * mrSges(3,3) + qJD(1) * t549) * t574 + t597;
t474 = t569 * t477 + t572 * t481;
t584 = -m(4) * t515 - t474;
t471 = m(3) * t524 + ((-mrSges(3,3) - t593) * qJDD(1) + (-t545 * t570 - t546 * t573 - t549) * qJD(1)) * t571 + t584;
t598 = t574 * t466 - t571 * t471;
t459 = m(2) * t553 - t579 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t598;
t467 = t570 * t473 + t573 * t475;
t544 = -qJDD(1) * pkin(1) + t586;
t581 = -m(3) * t544 + mrSges(3,1) * t606 - t467 + (t616 + t617) * mrSges(3,3);
t463 = m(2) * t554 - t579 * mrSges(2,2) + (mrSges(2,1) - t621) * qJDD(1) + t581;
t612 = t576 * t459 + t578 * t463;
t460 = t571 * t466 + t574 * t471;
t599 = -t578 * t459 + t576 * t463;
t592 = Ifges(3,1) * t571 + Ifges(3,4) * t574;
t591 = Ifges(4,5) * t573 - Ifges(4,6) * t570;
t582 = -Ifges(4,5) * t574 + (Ifges(4,1) * t573 - Ifges(4,4) * t570) * t571;
t550 = (Ifges(3,5) * t571 + t619) * qJD(1);
t530 = t582 * qJD(1);
t529 = (-t618 + (-Ifges(4,2) * t570 + t620) * t571) * qJD(1);
t528 = (-Ifges(4,3) * t574 + t591 * t571) * qJD(1);
t511 = Ifges(5,1) * t538 - Ifges(5,4) * t537 + Ifges(5,5) * t604;
t510 = Ifges(5,4) * t538 - Ifges(5,2) * t537 + Ifges(5,6) * t604;
t509 = Ifges(5,5) * t538 - Ifges(5,6) * t537 + Ifges(5,3) * t604;
t500 = Ifges(6,1) * t520 + Ifges(6,4) * t519 + Ifges(6,5) * t534;
t499 = Ifges(6,4) * t520 + Ifges(6,2) * t519 + Ifges(6,6) * t534;
t498 = Ifges(6,5) * t520 + Ifges(6,6) * t519 + Ifges(6,3) * t534;
t480 = mrSges(6,2) * t487 - mrSges(6,3) * t485 + Ifges(6,1) * t506 + Ifges(6,4) * t505 + Ifges(6,5) * t533 + t519 * t498 - t534 * t499;
t479 = -mrSges(6,1) * t487 + mrSges(6,3) * t486 + Ifges(6,4) * t506 + Ifges(6,2) * t505 + Ifges(6,6) * t533 - t520 * t498 + t534 * t500;
t469 = -mrSges(5,1) * t493 - mrSges(6,1) * t485 + mrSges(6,2) * t486 + mrSges(5,3) * t490 + Ifges(5,4) * t536 - Ifges(6,5) * t506 - Ifges(5,2) * t535 - Ifges(6,6) * t505 - Ifges(6,3) * t533 - pkin(4) * t478 - t520 * t499 + t519 * t500 - t538 * t509 + (Ifges(5,6) * qJDD(1) + qJD(1) * t511) * t615;
t468 = mrSges(5,2) * t493 - mrSges(5,3) * t489 + Ifges(5,1) * t536 - Ifges(5,4) * t535 - pkin(6) * t478 - t575 * t479 + t577 * t480 - t537 * t509 + (Ifges(5,5) * qJDD(1) - qJD(1) * t510) * t615;
t457 = -mrSges(4,1) * t515 + mrSges(4,3) * t497 - Ifges(5,5) * t536 + Ifges(5,6) * t535 - t538 * t510 - t537 * t511 - mrSges(5,1) * t489 + mrSges(5,2) * t490 - t575 * t480 - t577 * t479 - pkin(4) * t583 - pkin(6) * t595 - pkin(3) * t474 + (-t528 * t614 - t574 * t530) * qJD(1) + (-t618 + (t620 + (-Ifges(4,2) - Ifges(5,3)) * t570) * t571) * qJDD(1);
t456 = mrSges(4,2) * t515 - mrSges(4,3) * t496 - qJ(4) * t474 + t572 * t468 - t569 * t469 + (-t528 * t615 + t529 * t574) * qJD(1) + t582 * qJDD(1);
t455 = -mrSges(3,1) * t544 + mrSges(3,3) * t525 - mrSges(4,1) * t496 + mrSges(4,2) * t497 - t569 * t468 - t572 * t469 - pkin(3) * t580 - qJ(4) * t596 - pkin(2) * t467 + (Ifges(3,2) + Ifges(4,3)) * t606 + ((Ifges(3,4) - t591) * qJDD(1) + (-t529 * t573 - t530 * t570 - t550) * qJD(1)) * t571;
t454 = mrSges(3,2) * t544 - mrSges(3,3) * t524 - qJ(3) * t467 + t592 * qJDD(1) + t573 * t456 - t570 * t457 + t550 * t609;
t453 = t579 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t553 - mrSges(3,1) * t524 + mrSges(3,2) * t525 - t570 * t456 - t573 * t457 - pkin(2) * t584 - qJ(3) * t597 - pkin(1) * t460 + (-t619 + Ifges(2,6) + (pkin(2) * t593 - Ifges(3,5)) * t571) * qJDD(1) + (-pkin(2) * (-t545 * t615 - t546 * t614) + (-t571 * (Ifges(3,4) * t571 + Ifges(3,2) * t574) + t574 * t592) * qJD(1)) * qJD(1);
t452 = -mrSges(2,2) * g(1) - mrSges(2,3) * t554 + Ifges(2,5) * qJDD(1) - t579 * Ifges(2,6) - qJ(2) * t460 + t574 * t454 - t571 * t455;
t1 = [(-m(1) - m(2)) * g(1) + t460; -m(1) * g(2) + t612; -m(1) * g(3) + t599; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t554 - mrSges(2,2) * t553 + t571 * t454 + t574 * t455 + pkin(1) * (-mrSges(3,2) * t607 + t581) + qJ(2) * t598; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t599 + t576 * t452 + t578 * t453; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t612 - t578 * t452 + t576 * t453;];
tauB = t1;
