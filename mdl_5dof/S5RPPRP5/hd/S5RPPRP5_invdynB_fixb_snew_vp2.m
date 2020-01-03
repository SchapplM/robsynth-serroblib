% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRP5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:28
% EndTime: 2019-12-31 17:53:30
% DurationCPUTime: 1.60s
% Computational Cost: add. (9261->228), mult. (21819->274), div. (0->0), fcn. (12831->6), ass. (0->95)
t549 = cos(pkin(7));
t601 = (Ifges(3,6) - Ifges(4,6)) * t549;
t551 = sin(qJ(1));
t553 = cos(qJ(1));
t524 = -t553 * g(1) - t551 * g(2);
t555 = qJD(1) ^ 2;
t600 = -t555 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t524;
t599 = Ifges(5,1) + Ifges(6,1);
t598 = Ifges(4,4) + Ifges(3,5);
t590 = Ifges(5,4) - Ifges(6,5);
t589 = Ifges(5,5) + Ifges(6,4);
t597 = Ifges(5,2) + Ifges(6,3);
t587 = Ifges(5,6) - Ifges(6,6);
t595 = Ifges(5,3) + Ifges(6,2);
t548 = sin(pkin(7));
t542 = t548 ^ 2;
t543 = t549 ^ 2;
t584 = t543 * t555;
t594 = t542 * t555 + t584;
t523 = t551 * g(1) - t553 * g(2);
t514 = -qJDD(1) * pkin(1) - t555 * qJ(2) + qJDD(2) - t523;
t570 = qJDD(1) * t549;
t571 = qJDD(1) * t548;
t574 = t548 * qJD(1);
t498 = -pkin(2) * t570 - qJ(3) * t571 - 0.2e1 * qJD(3) * t574 + t514;
t501 = -t549 * g(3) - t600 * t548;
t593 = -mrSges(5,3) - mrSges(6,2);
t592 = Ifges(3,1) + Ifges(4,1);
t591 = Ifges(3,4) - Ifges(4,5);
t588 = Ifges(3,2) + Ifges(4,3);
t586 = mrSges(3,2) * t548;
t519 = (-mrSges(4,1) * t549 - mrSges(4,3) * t548) * qJD(1);
t520 = (-mrSges(3,1) * t549 + t586) * qJD(1);
t518 = (-pkin(2) * t549 - qJ(3) * t548) * qJD(1);
t478 = t518 * t574 + qJDD(3) - t501;
t473 = (-pkin(3) * t549 * t555 - pkin(6) * qJDD(1)) * t548 + t478;
t502 = -t548 * g(3) + t600 * t549;
t573 = t549 * qJD(1);
t480 = t518 * t573 + t502;
t475 = -pkin(3) * t584 - pkin(6) * t570 + t480;
t550 = sin(qJ(4));
t552 = cos(qJ(4));
t471 = t550 * t473 + t552 * t475;
t561 = t548 * t550 + t549 * t552;
t562 = t548 * t552 - t549 * t550;
t516 = t562 * qJD(1);
t575 = t516 * qJD(4);
t499 = t561 * qJDD(1) + t575;
t506 = qJD(4) * mrSges(5,1) - t516 * mrSges(5,3);
t515 = t561 * qJD(1);
t491 = t515 * pkin(4) - t516 * qJ(5);
t554 = qJD(4) ^ 2;
t466 = -t554 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t515 * t491 + t471;
t507 = -qJD(4) * mrSges(6,1) + t516 * mrSges(6,2);
t569 = m(6) * t466 + qJDD(4) * mrSges(6,3) + qJD(4) * t507;
t492 = t515 * mrSges(6,1) - t516 * mrSges(6,3);
t579 = -t515 * mrSges(5,1) - t516 * mrSges(5,2) - t492;
t462 = m(5) * t471 - qJDD(4) * mrSges(5,2) - qJD(4) * t506 + t593 * t499 + t579 * t515 + t569;
t470 = t552 * t473 - t550 * t475;
t576 = t515 * qJD(4);
t500 = t562 * qJDD(1) - t576;
t505 = -qJD(4) * mrSges(5,2) - t515 * mrSges(5,3);
t467 = -qJDD(4) * pkin(4) - t554 * qJ(5) + t516 * t491 + qJDD(5) - t470;
t508 = -t515 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t563 = -m(6) * t467 + qJDD(4) * mrSges(6,1) + qJD(4) * t508;
t463 = m(5) * t470 + qJDD(4) * mrSges(5,1) + qJD(4) * t505 + t593 * t500 + t579 * t516 + t563;
t457 = t550 * t462 + t552 * t463;
t558 = -m(4) * t478 - t457;
t453 = m(3) * t501 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t519 - t520) * qJD(1)) * t548 + t558;
t564 = t552 * t462 - t550 * t463;
t560 = m(4) * t480 + mrSges(4,2) * t570 + t519 * t573 + t564;
t454 = m(3) * t502 + (qJDD(1) * mrSges(3,3) + qJD(1) * t520) * t549 + t560;
t565 = -t548 * t453 + t549 * t454;
t448 = m(2) * t524 - t555 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t565;
t477 = pkin(3) * t570 + (-t542 - t543) * t555 * pkin(6) - t498;
t469 = -0.2e1 * qJD(5) * t516 + (-t500 + t576) * qJ(5) + (t499 + t575) * pkin(4) + t477;
t464 = m(6) * t469 + t499 * mrSges(6,1) - t500 * mrSges(6,3) - t516 * t507 + t515 * t508;
t557 = -m(5) * t477 - t499 * mrSges(5,1) - t500 * mrSges(5,2) - t515 * t505 - t516 * t506 - t464;
t460 = m(4) * t498 - mrSges(4,1) * t570 - t594 * mrSges(4,2) - mrSges(4,3) * t571 + t557;
t556 = -m(3) * t514 + mrSges(3,1) * t570 + t594 * mrSges(3,3) - t460;
t459 = (mrSges(2,1) - t586) * qJDD(1) + t556 - t555 * mrSges(2,2) + m(2) * t523;
t583 = t551 * t448 + t553 * t459;
t449 = t549 * t453 + t548 * t454;
t582 = -t587 * qJD(4) + t597 * t515 - t590 * t516;
t581 = -t595 * qJD(4) + t587 * t515 - t589 * t516;
t580 = t589 * qJD(4) - t590 * t515 + t599 * t516;
t577 = (t598 * t548 + t601) * qJD(1);
t566 = t553 * t448 - t551 * t459;
t456 = mrSges(5,2) * t477 + mrSges(6,2) * t467 - mrSges(5,3) * t470 - mrSges(6,3) * t469 - qJ(5) * t464 + t582 * qJD(4) + t589 * qJDD(4) - t590 * t499 + t599 * t500 + t581 * t515;
t455 = -mrSges(5,1) * t477 - mrSges(6,1) * t469 + mrSges(6,2) * t466 + mrSges(5,3) * t471 - pkin(4) * t464 + t580 * qJD(4) + t587 * qJDD(4) - t597 * t499 + t590 * t500 + t581 * t516;
t445 = mrSges(3,2) * t514 + mrSges(4,2) * t478 - mrSges(3,3) * t501 - mrSges(4,3) * t498 - pkin(6) * t457 - qJ(3) * t460 - t550 * t455 + t552 * t456 + t577 * t573 + (t592 * t548 + t591 * t549) * qJDD(1);
t444 = -mrSges(3,1) * t514 + mrSges(3,3) * t502 - mrSges(4,1) * t498 + mrSges(4,2) * t480 - t550 * t456 - t552 * t455 - pkin(3) * t557 - pkin(6) * t564 - pkin(2) * t460 - t577 * t574 + (t591 * t548 + t588 * t549) * qJDD(1);
t443 = mrSges(2,1) * g(3) - pkin(2) * t558 + t555 * Ifges(2,5) - qJ(3) * t560 + mrSges(2,3) * t524 + qJ(5) * t569 + pkin(4) * t563 - mrSges(3,1) * t501 + mrSges(3,2) * t502 - mrSges(4,3) * t480 - mrSges(5,2) * t471 + mrSges(4,1) * t478 + mrSges(5,1) * t470 - mrSges(6,1) * t467 + mrSges(6,3) * t466 + pkin(3) * t457 - pkin(1) * t449 + (-pkin(4) * t492 - t582) * t516 + (-qJ(5) * t492 + t580) * t515 + (-pkin(4) * mrSges(6,2) + t589) * t500 + (-qJ(5) * mrSges(6,2) - t587) * t499 + t595 * qJDD(4) + (Ifges(2,6) - t601 + (mrSges(4,2) * pkin(2) - t598) * t548) * qJDD(1) + (t591 * t543 * qJD(1) + (pkin(2) * t519 - t591 * t574 + (-t588 + t592) * t573) * t548) * qJD(1);
t442 = -mrSges(2,2) * g(3) - mrSges(2,3) * t523 + Ifges(2,5) * qJDD(1) - t555 * Ifges(2,6) - qJ(2) * t449 - t548 * t444 + t549 * t445;
t1 = [-m(1) * g(1) + t566; -m(1) * g(2) + t583; (-m(1) - m(2)) * g(3) + t449; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t583 + t553 * t442 - t551 * t443; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t566 + t551 * t442 + t553 * t443; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t523 - mrSges(2,2) * t524 + t548 * t445 + t549 * t444 + pkin(1) * (-mrSges(3,2) * t571 + t556) + qJ(2) * t565;];
tauB = t1;
