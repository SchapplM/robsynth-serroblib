% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPPRR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 13:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:36:02
% EndTime: 2019-05-05 13:36:07
% DurationCPUTime: 3.98s
% Computational Cost: add. (47297->271), mult. (98363->328), div. (0->0), fcn. (62511->10), ass. (0->119)
t570 = sin(qJ(1));
t573 = cos(qJ(1));
t545 = t570 * g(1) - t573 * g(2);
t543 = qJDD(1) * pkin(1) + t545;
t546 = -t573 * g(1) - t570 * g(2);
t575 = qJD(1) ^ 2;
t544 = -t575 * pkin(1) + t546;
t565 = sin(pkin(9));
t567 = cos(pkin(9));
t529 = t567 * t543 - t565 * t544;
t581 = -t575 * qJ(3) + qJDD(3) - t529;
t607 = -pkin(2) - qJ(4);
t613 = -(2 * qJD(1) * qJD(4)) + t607 * qJDD(1) + t581;
t564 = sin(pkin(10));
t557 = t564 ^ 2;
t566 = cos(pkin(10));
t604 = t566 ^ 2 + t557;
t595 = t604 * mrSges(5,3);
t530 = t565 * t543 + t567 * t544;
t612 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t530;
t569 = sin(qJ(5));
t572 = cos(qJ(5));
t586 = t564 * t572 + t566 * t569;
t536 = t586 * qJD(1);
t585 = -t564 * t569 + t566 * t572;
t537 = t585 * qJD(1);
t601 = t537 * qJD(5);
t527 = -t586 * qJDD(1) - t601;
t611 = pkin(4) * t575;
t610 = mrSges(3,1) - mrSges(4,2);
t609 = -Ifges(4,4) + Ifges(3,5);
t608 = Ifges(3,6) - Ifges(4,5);
t561 = -g(3) + qJDD(2);
t598 = qJDD(1) * t566;
t605 = t613 * t566;
t496 = -pkin(7) * t598 + (-t566 * t611 - t561) * t564 + t605;
t508 = t566 * t561 + t613 * t564;
t599 = qJDD(1) * t564;
t497 = -pkin(7) * t599 - t557 * t611 + t508;
t493 = t569 * t496 + t572 * t497;
t523 = t536 * mrSges(6,1) + t537 * mrSges(6,2);
t534 = qJD(5) * mrSges(6,1) - t537 * mrSges(6,3);
t526 = t536 * pkin(5) - t537 * pkin(8);
t574 = qJD(5) ^ 2;
t490 = -t574 * pkin(5) + qJDD(5) * pkin(8) - t536 * t526 + t493;
t584 = qJDD(4) + t612;
t506 = pkin(4) * t599 + (-t604 * pkin(7) + t607) * t575 + t584;
t602 = t536 * qJD(5);
t528 = t585 * qJDD(1) - t602;
t491 = (-t528 + t602) * pkin(8) + (-t527 + t601) * pkin(5) + t506;
t568 = sin(qJ(6));
t571 = cos(qJ(6));
t487 = -t568 * t490 + t571 * t491;
t531 = t571 * qJD(5) - t568 * t537;
t504 = t531 * qJD(6) + t568 * qJDD(5) + t571 * t528;
t532 = t568 * qJD(5) + t571 * t537;
t509 = -t531 * mrSges(7,1) + t532 * mrSges(7,2);
t535 = qJD(6) + t536;
t514 = -t535 * mrSges(7,2) + t531 * mrSges(7,3);
t525 = qJDD(6) - t527;
t485 = m(7) * t487 + t525 * mrSges(7,1) - t504 * mrSges(7,3) - t532 * t509 + t535 * t514;
t488 = t571 * t490 + t568 * t491;
t503 = -t532 * qJD(6) + t571 * qJDD(5) - t568 * t528;
t515 = t535 * mrSges(7,1) - t532 * mrSges(7,3);
t486 = m(7) * t488 - t525 * mrSges(7,2) + t503 * mrSges(7,3) + t531 * t509 - t535 * t515;
t590 = -t568 * t485 + t571 * t486;
t476 = m(6) * t493 - qJDD(5) * mrSges(6,2) + t527 * mrSges(6,3) - qJD(5) * t534 - t536 * t523 + t590;
t492 = t572 * t496 - t569 * t497;
t533 = -qJD(5) * mrSges(6,2) - t536 * mrSges(6,3);
t489 = -qJDD(5) * pkin(5) - t574 * pkin(8) + t537 * t526 - t492;
t579 = -m(7) * t489 + t503 * mrSges(7,1) - t504 * mrSges(7,2) + t531 * t514 - t532 * t515;
t481 = m(6) * t492 + qJDD(5) * mrSges(6,1) - t528 * mrSges(6,3) + qJD(5) * t533 - t537 * t523 + t579;
t470 = t569 * t476 + t572 * t481;
t507 = -t564 * t561 + t605;
t583 = -qJDD(1) * mrSges(5,3) - t575 * (mrSges(5,1) * t564 + mrSges(5,2) * t566);
t468 = m(5) * t507 + t583 * t566 + t470;
t591 = t572 * t476 - t569 * t481;
t469 = m(5) * t508 + t583 * t564 + t591;
t464 = t566 * t468 + t564 * t469;
t518 = -qJDD(1) * pkin(2) + t581;
t580 = -m(4) * t518 + t575 * mrSges(4,3) - t464;
t462 = m(3) * t529 - t575 * mrSges(3,2) + t610 * qJDD(1) + t580;
t517 = t575 * pkin(2) - t612;
t513 = t607 * t575 + t584;
t477 = t571 * t485 + t568 * t486;
t578 = m(6) * t506 - t527 * mrSges(6,1) + t528 * mrSges(6,2) + t536 * t533 + t537 * t534 + t477;
t577 = -m(5) * t513 - mrSges(5,1) * t599 - mrSges(5,2) * t598 - t578;
t576 = -m(4) * t517 + t575 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t577;
t473 = t576 + m(3) * t530 - qJDD(1) * mrSges(3,2) + (-mrSges(3,1) - t595) * t575;
t460 = t567 * t462 + t565 * t473;
t458 = m(2) * t545 + qJDD(1) * mrSges(2,1) - t575 * mrSges(2,2) + t460;
t593 = -t565 * t462 + t567 * t473;
t459 = m(2) * t546 - t575 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t593;
t606 = t573 * t458 + t570 * t459;
t587 = Ifges(5,5) * t566 - Ifges(5,6) * t564;
t603 = t575 * t587;
t594 = -t570 * t458 + t573 * t459;
t592 = -t564 * t468 + t566 * t469;
t463 = m(4) * t561 + t592;
t589 = Ifges(5,1) * t566 - Ifges(5,4) * t564;
t588 = Ifges(5,4) * t566 - Ifges(5,2) * t564;
t582 = m(3) * t561 + t463;
t521 = Ifges(6,1) * t537 - Ifges(6,4) * t536 + Ifges(6,5) * qJD(5);
t520 = Ifges(6,4) * t537 - Ifges(6,2) * t536 + Ifges(6,6) * qJD(5);
t519 = Ifges(6,5) * t537 - Ifges(6,6) * t536 + Ifges(6,3) * qJD(5);
t500 = Ifges(7,1) * t532 + Ifges(7,4) * t531 + Ifges(7,5) * t535;
t499 = Ifges(7,4) * t532 + Ifges(7,2) * t531 + Ifges(7,6) * t535;
t498 = Ifges(7,5) * t532 + Ifges(7,6) * t531 + Ifges(7,3) * t535;
t479 = mrSges(7,2) * t489 - mrSges(7,3) * t487 + Ifges(7,1) * t504 + Ifges(7,4) * t503 + Ifges(7,5) * t525 + t531 * t498 - t535 * t499;
t478 = -mrSges(7,1) * t489 + mrSges(7,3) * t488 + Ifges(7,4) * t504 + Ifges(7,2) * t503 + Ifges(7,6) * t525 - t532 * t498 + t535 * t500;
t466 = -mrSges(6,1) * t506 - mrSges(7,1) * t487 + mrSges(7,2) * t488 + mrSges(6,3) * t493 + Ifges(6,4) * t528 - Ifges(7,5) * t504 + Ifges(6,2) * t527 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t503 - Ifges(7,3) * t525 - pkin(5) * t477 + qJD(5) * t521 - t532 * t499 + t531 * t500 - t537 * t519;
t465 = mrSges(6,2) * t506 - mrSges(6,3) * t492 + Ifges(6,1) * t528 + Ifges(6,4) * t527 + Ifges(6,5) * qJDD(5) - pkin(8) * t477 - qJD(5) * t520 - t568 * t478 + t571 * t479 - t536 * t519;
t454 = mrSges(5,2) * t513 - mrSges(5,3) * t507 - pkin(7) * t470 + t589 * qJDD(1) + t572 * t465 - t569 * t466 - t564 * t603;
t453 = -mrSges(5,1) * t513 + mrSges(5,3) * t508 - pkin(4) * t578 + pkin(7) * t591 + t588 * qJDD(1) + t569 * t465 + t572 * t466 - t566 * t603;
t452 = -qJ(3) * t463 + pkin(4) * t470 + pkin(3) * t464 + mrSges(6,1) * t492 - mrSges(6,2) * t493 + t571 * t478 + pkin(8) * t590 + t568 * t479 + t536 * t521 + t537 * t520 - mrSges(3,3) * t529 + pkin(5) * t579 + Ifges(6,6) * t527 + Ifges(6,5) * t528 + Ifges(6,3) * qJDD(5) + mrSges(5,1) * t507 - mrSges(5,2) * t508 + mrSges(4,1) * t518 + (mrSges(3,2) - mrSges(4,3)) * t561 + (t587 + t609) * qJDD(1) + (t564 * t589 + t566 * t588 - t608) * t575;
t451 = mrSges(3,3) * t530 - mrSges(4,1) * t517 - t564 * t454 - t566 * t453 - pkin(3) * t577 - qJ(4) * t592 - pkin(2) * t463 - t610 * t561 + t608 * qJDD(1) + (-pkin(3) * t595 + t609) * t575;
t450 = -mrSges(2,2) * g(3) - mrSges(2,3) * t545 + Ifges(2,5) * qJDD(1) - t575 * Ifges(2,6) - qJ(2) * t460 - t565 * t451 + t567 * t452;
t449 = mrSges(2,1) * g(3) + mrSges(2,3) * t546 + t575 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t582 + qJ(2) * t593 + t567 * t451 + t565 * t452;
t1 = [-m(1) * g(1) + t594; -m(1) * g(2) + t606; (-m(1) - m(2)) * g(3) + t582; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t606 - t570 * t449 + t573 * t450; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t594 + t573 * t449 + t570 * t450; pkin(1) * t460 - mrSges(2,2) * t546 + mrSges(2,1) * t545 + pkin(2) * t580 + qJ(3) * (-t575 * t595 + t576) - qJ(4) * t464 + mrSges(3,1) * t529 - mrSges(3,2) * t530 + t566 * t454 - t564 * t453 + mrSges(4,2) * t518 - mrSges(4,3) * t517 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
