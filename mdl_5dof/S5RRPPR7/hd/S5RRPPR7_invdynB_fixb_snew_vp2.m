% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:35:14
% EndTime: 2019-12-31 19:35:20
% DurationCPUTime: 3.12s
% Computational Cost: add. (27257->296), mult. (62678->357), div. (0->0), fcn. (39622->8), ass. (0->117)
t618 = -2 * qJD(3);
t617 = Ifges(4,1) + Ifges(5,2);
t616 = -Ifges(5,1) - Ifges(4,3);
t612 = Ifges(4,4) + Ifges(5,6);
t611 = Ifges(4,5) - Ifges(5,4);
t615 = Ifges(4,2) + Ifges(5,3);
t610 = Ifges(4,6) - Ifges(5,5);
t578 = sin(qJ(2));
t581 = cos(qJ(2));
t597 = qJD(1) * qJD(2);
t565 = t578 * qJDD(1) + t581 * t597;
t579 = sin(qJ(1));
t582 = cos(qJ(1));
t571 = -t582 * g(1) - t579 * g(2);
t584 = qJD(1) ^ 2;
t560 = -t584 * pkin(1) + qJDD(1) * pkin(6) + t571;
t608 = t578 * t560;
t613 = pkin(2) * t584;
t509 = qJDD(2) * pkin(2) - t565 * qJ(3) - t608 + (qJ(3) * t597 + t578 * t613 - g(3)) * t581;
t540 = -t578 * g(3) + t581 * t560;
t566 = t581 * qJDD(1) - t578 * t597;
t601 = qJD(1) * t578;
t567 = qJD(2) * pkin(2) - qJ(3) * t601;
t575 = t581 ^ 2;
t510 = t566 * qJ(3) - qJD(2) * t567 - t575 * t613 + t540;
t576 = sin(pkin(8));
t609 = cos(pkin(8));
t554 = (t576 * t581 + t609 * t578) * qJD(1);
t497 = t609 * t509 - t576 * t510 + t554 * t618;
t614 = -2 * qJD(4);
t600 = qJD(1) * t581;
t553 = t576 * t601 - t609 * t600;
t527 = t553 * mrSges(4,1) + t554 * mrSges(4,2);
t535 = t609 * t565 + t576 * t566;
t541 = -qJD(2) * mrSges(4,2) - t553 * mrSges(4,3);
t543 = t553 * mrSges(5,1) - qJD(2) * mrSges(5,3);
t526 = t553 * pkin(3) - t554 * qJ(4);
t583 = qJD(2) ^ 2;
t494 = -qJDD(2) * pkin(3) - t583 * qJ(4) + t554 * t526 + qJDD(4) - t497;
t599 = qJD(2) * t553;
t489 = (t553 * t554 - qJDD(2)) * pkin(7) + (t535 + t599) * pkin(4) + t494;
t534 = t576 * t565 - t609 * t566;
t545 = t554 * pkin(4) - qJD(2) * pkin(7);
t552 = t553 ^ 2;
t570 = t579 * g(1) - t582 * g(2);
t592 = -qJDD(1) * pkin(1) - t570;
t514 = -t566 * pkin(2) + qJDD(3) + t567 * t601 + (-qJ(3) * t575 - pkin(6)) * t584 + t592;
t586 = (-t535 + t599) * qJ(4) + t514 + (qJD(2) * pkin(3) + t614) * t554;
t492 = -t552 * pkin(4) - t554 * t545 + (pkin(3) + pkin(7)) * t534 + t586;
t577 = sin(qJ(5));
t580 = cos(qJ(5));
t487 = t580 * t489 - t577 * t492;
t536 = -t577 * qJD(2) + t580 * t553;
t505 = t536 * qJD(5) + t580 * qJDD(2) + t577 * t534;
t537 = t580 * qJD(2) + t577 * t553;
t511 = -t536 * mrSges(6,1) + t537 * mrSges(6,2);
t551 = qJD(5) + t554;
t515 = -t551 * mrSges(6,2) + t536 * mrSges(6,3);
t533 = qJDD(5) + t535;
t485 = m(6) * t487 + t533 * mrSges(6,1) - t505 * mrSges(6,3) - t537 * t511 + t551 * t515;
t488 = t577 * t489 + t580 * t492;
t504 = -t537 * qJD(5) - t577 * qJDD(2) + t580 * t534;
t516 = t551 * mrSges(6,1) - t537 * mrSges(6,3);
t486 = m(6) * t488 - t533 * mrSges(6,2) + t504 * mrSges(6,3) + t536 * t511 - t551 * t516;
t477 = t580 * t485 + t577 * t486;
t528 = -t553 * mrSges(5,2) - t554 * mrSges(5,3);
t589 = -m(5) * t494 - t535 * mrSges(5,1) - t554 * t528 - t477;
t475 = m(4) * t497 - t535 * mrSges(4,3) - t554 * t527 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t541 - t543) * qJD(2) + t589;
t548 = t553 * t618;
t605 = t576 * t509 + t609 * t510;
t498 = t548 + t605;
t542 = qJD(2) * mrSges(4,1) - t554 * mrSges(4,3);
t591 = t583 * pkin(3) - qJDD(2) * qJ(4) - t605;
t493 = qJD(2) * t614 + ((2 * qJD(3)) + t526) * t553 + t591;
t544 = t554 * mrSges(5,1) + qJD(2) * mrSges(5,2);
t491 = -t534 * pkin(4) - t552 * pkin(7) - t553 * t526 + t548 + ((2 * qJD(4)) + t545) * qJD(2) - t591;
t590 = -m(6) * t491 + t504 * mrSges(6,1) - t505 * mrSges(6,2) + t536 * t515 - t537 * t516;
t588 = -m(5) * t493 + qJDD(2) * mrSges(5,3) + qJD(2) * t544 - t590;
t482 = m(4) * t498 - qJDD(2) * mrSges(4,2) - qJD(2) * t542 + (-t527 - t528) * t553 + (-mrSges(4,3) - mrSges(5,1)) * t534 + t588;
t471 = t609 * t475 + t576 * t482;
t539 = -t581 * g(3) - t608;
t564 = (-mrSges(3,1) * t581 + mrSges(3,2) * t578) * qJD(1);
t569 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t600;
t469 = m(3) * t539 + qJDD(2) * mrSges(3,1) - t565 * mrSges(3,3) + qJD(2) * t569 - t564 * t601 + t471;
t568 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t601;
t594 = -t576 * t475 + t609 * t482;
t470 = m(3) * t540 - qJDD(2) * mrSges(3,2) + t566 * mrSges(3,3) - qJD(2) * t568 + t564 * t600 + t594;
t595 = -t578 * t469 + t581 * t470;
t463 = m(2) * t571 - t584 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t595;
t559 = -t584 * pkin(6) + t592;
t496 = t534 * pkin(3) + t586;
t606 = -t577 * t485 + t580 * t486;
t476 = m(5) * t496 - t534 * mrSges(5,2) - t535 * mrSges(5,3) - t553 * t543 - t554 * t544 + t606;
t587 = m(4) * t514 + t534 * mrSges(4,1) + t535 * mrSges(4,2) + t553 * t541 + t554 * t542 + t476;
t585 = -m(3) * t559 + t566 * mrSges(3,1) - t565 * mrSges(3,2) - t568 * t601 + t569 * t600 - t587;
t473 = m(2) * t570 + qJDD(1) * mrSges(2,1) - t584 * mrSges(2,2) + t585;
t607 = t579 * t463 + t582 * t473;
t464 = t581 * t469 + t578 * t470;
t604 = t616 * qJD(2) + t610 * t553 - t611 * t554;
t603 = -t610 * qJD(2) + t615 * t553 - t612 * t554;
t602 = t611 * qJD(2) - t612 * t553 + t617 * t554;
t596 = t582 * t463 - t579 * t473;
t557 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t578 + Ifges(3,4) * t581) * qJD(1);
t556 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t578 + Ifges(3,2) * t581) * qJD(1);
t555 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t578 + Ifges(3,6) * t581) * qJD(1);
t501 = Ifges(6,1) * t537 + Ifges(6,4) * t536 + Ifges(6,5) * t551;
t500 = Ifges(6,4) * t537 + Ifges(6,2) * t536 + Ifges(6,6) * t551;
t499 = Ifges(6,5) * t537 + Ifges(6,6) * t536 + Ifges(6,3) * t551;
t479 = mrSges(6,2) * t491 - mrSges(6,3) * t487 + Ifges(6,1) * t505 + Ifges(6,4) * t504 + Ifges(6,5) * t533 + t536 * t499 - t551 * t500;
t478 = -mrSges(6,1) * t491 + mrSges(6,3) * t488 + Ifges(6,4) * t505 + Ifges(6,2) * t504 + Ifges(6,6) * t533 - t537 * t499 + t551 * t501;
t465 = mrSges(5,1) * t494 + mrSges(6,1) * t487 + mrSges(4,2) * t514 - mrSges(6,2) * t488 - mrSges(4,3) * t497 - mrSges(5,3) * t496 + Ifges(6,5) * t505 + Ifges(6,6) * t504 + Ifges(6,3) * t533 + pkin(4) * t477 - qJ(4) * t476 + t537 * t500 - t536 * t501 + t604 * t553 + t617 * t535 - t612 * t534 + t611 * qJDD(2) + t603 * qJD(2);
t460 = -mrSges(4,1) * t514 - mrSges(5,1) * t493 + mrSges(5,2) * t496 + mrSges(4,3) * t498 - pkin(3) * t476 - pkin(4) * t590 - pkin(7) * t606 + t602 * qJD(2) + t610 * qJDD(2) - t580 * t478 - t577 * t479 - t615 * t534 + t612 * t535 + t604 * t554;
t459 = mrSges(3,2) * t559 - mrSges(3,3) * t539 + Ifges(3,1) * t565 + Ifges(3,4) * t566 + Ifges(3,5) * qJDD(2) - qJ(3) * t471 - qJD(2) * t556 - t576 * t460 + t609 * t465 + t555 * t600;
t458 = -qJ(4) * t588 - pkin(3) * (-qJD(2) * t543 + t589) + (qJ(4) * t528 - t602) * t553 + t603 * t554 + (pkin(3) * mrSges(5,2) - Ifges(3,3) + t616) * qJDD(2) + (qJ(4) * mrSges(5,1) + t610) * t534 - t611 * t535 + t584 * Ifges(2,5) + t577 * t478 - t580 * t479 + pkin(7) * t477 - mrSges(3,1) * t539 + mrSges(3,2) * t540 - mrSges(4,1) * t497 + mrSges(4,2) * t498 + mrSges(5,3) * t493 - mrSges(5,2) * t494 - Ifges(3,5) * t565 - Ifges(3,6) * t566 + mrSges(2,3) * t571 - pkin(2) * t471 + Ifges(2,6) * qJDD(1) - pkin(1) * t464 + mrSges(2,1) * g(3) + (-t578 * t556 + t581 * t557) * qJD(1);
t457 = -mrSges(3,1) * t559 + mrSges(3,3) * t540 + Ifges(3,4) * t565 + Ifges(3,2) * t566 + Ifges(3,6) * qJDD(2) - pkin(2) * t587 + qJ(3) * t594 + qJD(2) * t557 + t609 * t460 + t576 * t465 - t555 * t601;
t456 = -mrSges(2,2) * g(3) - mrSges(2,3) * t570 + Ifges(2,5) * qJDD(1) - t584 * Ifges(2,6) - pkin(6) * t464 - t578 * t457 + t581 * t459;
t1 = [-m(1) * g(1) + t596; -m(1) * g(2) + t607; (-m(1) - m(2)) * g(3) + t464; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t607 + t582 * t456 - t579 * t458; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t596 + t579 * t456 + t582 * t458; -mrSges(1,1) * g(2) + mrSges(2,1) * t570 + mrSges(1,2) * g(1) - mrSges(2,2) * t571 + Ifges(2,3) * qJDD(1) + pkin(1) * t585 + pkin(6) * t595 + t581 * t457 + t578 * t459;];
tauB = t1;
