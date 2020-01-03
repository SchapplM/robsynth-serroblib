% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:13
% EndTime: 2019-12-31 17:45:14
% DurationCPUTime: 1.41s
% Computational Cost: add. (13790->197), mult. (26765->239), div. (0->0), fcn. (14909->8), ass. (0->93)
t570 = sin(qJ(1));
t572 = cos(qJ(1));
t543 = t570 * g(1) - t572 * g(2);
t540 = qJDD(1) * pkin(1) + t543;
t544 = -t572 * g(1) - t570 * g(2);
t573 = qJD(1) ^ 2;
t541 = -t573 * pkin(1) + t544;
t566 = sin(pkin(7));
t568 = cos(pkin(7));
t529 = t568 * t540 - t566 * t541;
t579 = -t573 * qJ(3) + qJDD(3) - t529;
t601 = -pkin(2) - qJ(4);
t609 = -(2 * qJD(1) * qJD(4)) + t601 * qJDD(1) + t579;
t530 = t566 * t540 + t568 * t541;
t608 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t530;
t515 = t573 * pkin(2) - t608;
t607 = -m(4) * t515 + t573 * mrSges(4,2) + qJDD(1) * mrSges(4,3);
t605 = pkin(4) * t573;
t604 = mrSges(3,1) - mrSges(4,2);
t603 = -Ifges(4,4) + Ifges(3,5);
t602 = Ifges(3,6) - Ifges(4,5);
t562 = -g(3) + qJDD(2);
t565 = sin(pkin(8));
t567 = cos(pkin(8));
t594 = qJDD(1) * t567;
t599 = t609 * t567;
t503 = -pkin(6) * t594 + (-t567 * t605 - t562) * t565 + t599;
t508 = t567 * t562 + t609 * t565;
t554 = t565 ^ 2;
t595 = qJDD(1) * t565;
t504 = -pkin(6) * t595 - t554 * t605 + t508;
t569 = sin(qJ(5));
t571 = cos(qJ(5));
t501 = t571 * t503 - t569 * t504;
t583 = -t565 * t571 - t567 * t569;
t533 = t583 * qJD(1);
t582 = -t565 * t569 + t567 * t571;
t534 = t582 * qJD(1);
t522 = -t533 * mrSges(6,1) + t534 * mrSges(6,2);
t527 = t533 * qJD(5) + t582 * qJDD(1);
t531 = -qJD(5) * mrSges(6,2) + t533 * mrSges(6,3);
t498 = m(6) * t501 + qJDD(5) * mrSges(6,1) - t527 * mrSges(6,3) + qJD(5) * t531 - t534 * t522;
t502 = t569 * t503 + t571 * t504;
t526 = -t534 * qJD(5) + t583 * qJDD(1);
t532 = qJD(5) * mrSges(6,1) - t534 * mrSges(6,3);
t499 = m(6) * t502 - qJDD(5) * mrSges(6,2) + t526 * mrSges(6,3) - qJD(5) * t532 + t533 * t522;
t487 = t571 * t498 + t569 * t499;
t507 = -t565 * t562 + t599;
t580 = -qJDD(1) * mrSges(5,3) - t573 * (mrSges(5,1) * t565 + mrSges(5,2) * t567);
t485 = m(5) * t507 + t580 * t567 + t487;
t587 = -t569 * t498 + t571 * t499;
t486 = m(5) * t508 + t580 * t565 + t587;
t483 = t567 * t485 + t565 * t486;
t517 = -qJDD(1) * pkin(2) + t579;
t577 = -m(4) * t517 + t573 * mrSges(4,3) - t483;
t478 = m(3) * t529 - t573 * mrSges(3,2) + t604 * qJDD(1) + t577;
t581 = qJDD(4) + t608;
t513 = t601 * t573 + t581;
t598 = -t567 ^ 2 - t554;
t506 = pkin(4) * t595 + (t598 * pkin(6) + t601) * t573 + t581;
t578 = m(6) * t506 - t526 * mrSges(6,1) + t527 * mrSges(6,2) - t533 * t531 + t534 * t532;
t575 = m(5) * t513 + mrSges(5,1) * t595 + mrSges(5,2) * t594 + t578;
t591 = t598 * mrSges(5,3);
t492 = m(3) * t530 - qJDD(1) * mrSges(3,2) + (-mrSges(3,1) + t591) * t573 + t575 + t607;
t473 = t568 * t478 + t566 * t492;
t470 = m(2) * t543 + qJDD(1) * mrSges(2,1) - t573 * mrSges(2,2) + t473;
t589 = -t566 * t478 + t568 * t492;
t471 = m(2) * t544 - t573 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t589;
t600 = t572 * t470 + t570 * t471;
t584 = Ifges(5,5) * t567 - Ifges(5,6) * t565;
t597 = t573 * t584;
t590 = -t570 * t470 + t572 * t471;
t588 = -t565 * t485 + t567 * t486;
t482 = m(4) * t562 + t588;
t586 = Ifges(5,1) * t567 - Ifges(5,4) * t565;
t585 = Ifges(5,4) * t567 - Ifges(5,2) * t565;
t481 = m(3) * t562 + t482;
t519 = Ifges(6,4) * t534 + Ifges(6,2) * t533 + Ifges(6,6) * qJD(5);
t520 = Ifges(6,1) * t534 + Ifges(6,4) * t533 + Ifges(6,5) * qJD(5);
t576 = mrSges(6,1) * t501 - mrSges(6,2) * t502 + Ifges(6,5) * t527 + Ifges(6,6) * t526 + Ifges(6,3) * qJDD(5) + t534 * t519 - t533 * t520;
t494 = t573 * t591 + t575;
t518 = Ifges(6,5) * t534 + Ifges(6,6) * t533 + Ifges(6,3) * qJD(5);
t488 = -mrSges(6,1) * t506 + mrSges(6,3) * t502 + Ifges(6,4) * t527 + Ifges(6,2) * t526 + Ifges(6,6) * qJDD(5) + qJD(5) * t520 - t534 * t518;
t489 = mrSges(6,2) * t506 - mrSges(6,3) * t501 + Ifges(6,1) * t527 + Ifges(6,4) * t526 + Ifges(6,5) * qJDD(5) - qJD(5) * t519 + t533 * t518;
t474 = -mrSges(5,1) * t513 + mrSges(5,3) * t508 - pkin(4) * t578 + pkin(6) * t587 + t585 * qJDD(1) + t571 * t488 + t569 * t489 - t567 * t597;
t476 = mrSges(5,2) * t513 - mrSges(5,3) * t507 - pkin(6) * t487 + t586 * qJDD(1) - t569 * t488 + t571 * t489 - t565 * t597;
t480 = qJDD(1) * mrSges(4,2) - t577;
t574 = -mrSges(2,2) * t544 - mrSges(3,2) * t530 - mrSges(4,3) * t515 + pkin(1) * t473 - pkin(2) * t480 - qJ(4) * t483 - t565 * t474 + t567 * t476 + qJ(3) * (t494 + t607) + mrSges(4,2) * t517 + mrSges(3,1) * t529 + mrSges(2,1) * t543 + (Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t466 = t576 + (mrSges(3,2) - mrSges(4,3)) * t562 + (t584 + t603) * qJDD(1) + mrSges(4,1) * t517 - mrSges(3,3) * t529 + mrSges(5,1) * t507 - mrSges(5,2) * t508 + pkin(4) * t487 + pkin(3) * t483 - qJ(3) * t482 + (t565 * t586 + t567 * t585 - t602) * t573;
t465 = -mrSges(4,1) * t515 + mrSges(3,3) * t530 - pkin(2) * t482 + pkin(3) * t494 - qJ(4) * t588 + t602 * qJDD(1) - t567 * t474 - t565 * t476 - t604 * t562 + t603 * t573;
t464 = -mrSges(2,2) * g(3) - mrSges(2,3) * t543 + Ifges(2,5) * qJDD(1) - t573 * Ifges(2,6) - qJ(2) * t473 - t566 * t465 + t568 * t466;
t463 = mrSges(2,1) * g(3) + mrSges(2,3) * t544 + t573 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t481 + qJ(2) * t589 + t568 * t465 + t566 * t466;
t1 = [-m(1) * g(1) + t590; -m(1) * g(2) + t600; (-m(1) - m(2)) * g(3) + t481; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t600 - t570 * t463 + t572 * t464; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t590 + t572 * t463 + t570 * t464; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t574; t574; t481; t480; t494; t576;];
tauJB = t1;
