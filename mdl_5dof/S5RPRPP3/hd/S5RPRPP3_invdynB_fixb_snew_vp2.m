% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:29
% EndTime: 2019-12-31 18:12:32
% DurationCPUTime: 1.86s
% Computational Cost: add. (12830->246), mult. (30764->287), div. (0->0), fcn. (19260->6), ass. (0->103)
t602 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t585 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t584 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t601 = Ifges(4,2) + Ifges(6,2) + Ifges(5,3);
t583 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t600 = -Ifges(4,3) - Ifges(5,1) - Ifges(6,1);
t559 = qJD(1) ^ 2;
t599 = -2 * qJD(4);
t598 = cos(qJ(3));
t554 = cos(pkin(7));
t597 = pkin(2) * t554;
t596 = -mrSges(6,1) - mrSges(4,3);
t553 = sin(pkin(7));
t595 = mrSges(3,2) * t553;
t549 = t554 ^ 2;
t594 = t549 * t559;
t556 = sin(qJ(1));
t557 = cos(qJ(1));
t539 = -t557 * g(1) - t556 * g(2);
t535 = -t559 * pkin(1) + qJDD(1) * qJ(2) + t539;
t588 = qJD(1) * qJD(2);
t577 = -t554 * g(3) - 0.2e1 * t553 * t588;
t487 = (-pkin(6) * qJDD(1) + t559 * t597 - t535) * t553 + t577;
t517 = -t553 * g(3) + (t535 + 0.2e1 * t588) * t554;
t586 = qJDD(1) * t554;
t488 = -pkin(2) * t594 + pkin(6) * t586 + t517;
t555 = sin(qJ(3));
t480 = t598 * t487 - t555 * t488;
t578 = t554 * t598;
t589 = t553 * qJD(1);
t533 = -qJD(1) * t578 + t555 * t589;
t566 = t598 * t553 + t554 * t555;
t534 = t566 * qJD(1);
t503 = -t534 * mrSges(6,2) + t533 * mrSges(6,3);
t505 = t533 * mrSges(4,1) + t534 * mrSges(4,2);
t591 = t533 * qJD(3);
t515 = t566 * qJDD(1) - t591;
t521 = -qJD(3) * mrSges(4,2) - t533 * mrSges(4,3);
t525 = t533 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t504 = t533 * pkin(3) - t534 * qJ(4);
t558 = qJD(3) ^ 2;
t479 = -qJDD(3) * pkin(3) - t558 * qJ(4) + t534 * t504 + qJDD(4) - t480;
t506 = -t533 * mrSges(5,2) - t534 * mrSges(5,3);
t473 = -0.2e1 * qJD(5) * qJD(3) + (t533 * t534 - qJDD(3)) * qJ(5) + (t515 + t591) * pkin(4) + t479;
t526 = -t533 * mrSges(6,1) + qJD(3) * mrSges(6,2);
t573 = m(6) * t473 - qJDD(3) * mrSges(6,3) - qJD(3) * t526;
t564 = -m(5) * t479 - t515 * mrSges(5,1) - t534 * t506 - t573;
t465 = m(4) * t480 + (-t503 - t505) * t534 + t596 * t515 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t521 - t525) * qJD(3) + t564;
t481 = t555 * t487 + t598 * t488;
t587 = qJDD(1) * t553;
t590 = t534 * qJD(3);
t514 = -qJDD(1) * t578 + t555 * t587 + t590;
t522 = qJD(3) * mrSges(4,1) - t534 * mrSges(4,3);
t563 = -t558 * pkin(3) + qJDD(3) * qJ(4) - t533 * t504 + t481;
t478 = qJD(3) * t599 - t563;
t527 = t534 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t523 = t534 * pkin(4) - qJD(3) * qJ(5);
t532 = t533 ^ 2;
t475 = -t514 * pkin(4) - t532 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t523) * qJD(3) + t563;
t524 = t534 * mrSges(6,1) - qJD(3) * mrSges(6,3);
t582 = m(6) * t475 + qJDD(3) * mrSges(6,2) + qJD(3) * t524;
t565 = -m(5) * t478 + qJDD(3) * mrSges(5,3) + qJD(3) * t527 + t582;
t592 = -t503 - t506;
t468 = m(4) * t481 - qJDD(3) * mrSges(4,2) - qJD(3) * t522 + (-t505 + t592) * t533 + (-mrSges(5,1) + t596) * t514 + t565;
t461 = t598 * t465 + t555 * t468;
t516 = -t553 * t535 + t577;
t567 = mrSges(3,3) * qJDD(1) + t559 * (-mrSges(3,1) * t554 + t595);
t459 = m(3) * t516 - t567 * t553 + t461;
t574 = -t555 * t465 + t598 * t468;
t460 = m(3) * t517 + t567 * t554 + t574;
t575 = -t553 * t459 + t554 * t460;
t452 = m(2) * t539 - t559 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t575;
t538 = t556 * g(1) - t557 * g(2);
t572 = qJDD(2) - t538;
t531 = -qJDD(1) * pkin(1) - t559 * qJ(2) + t572;
t548 = t553 ^ 2;
t513 = (-pkin(1) - t597) * qJDD(1) + (-qJ(2) + (-t548 - t549) * pkin(6)) * t559 + t572;
t560 = pkin(3) * t590 + t534 * t599 + (-t515 + t591) * qJ(4) + t513;
t477 = t514 * pkin(3) + t560;
t472 = -t532 * pkin(4) + 0.2e1 * qJD(5) * t533 - t534 * t523 + (pkin(3) + qJ(5)) * t514 + t560;
t568 = m(6) * t472 - t515 * mrSges(6,2) + t514 * mrSges(6,3) - t534 * t524 + t533 * t526;
t469 = m(5) * t477 - t514 * mrSges(5,2) - t515 * mrSges(5,3) - t533 * t525 - t534 * t527 + t568;
t562 = m(4) * t513 + t514 * mrSges(4,1) + t515 * mrSges(4,2) + t533 * t521 + t534 * t522 + t469;
t561 = -m(3) * t531 + mrSges(3,1) * t586 - t562 + (t548 * t559 + t594) * mrSges(3,3);
t463 = (mrSges(2,1) - t595) * qJDD(1) + t561 - t559 * mrSges(2,2) + m(2) * t538;
t593 = t556 * t452 + t557 * t463;
t453 = t554 * t459 + t553 * t460;
t581 = t600 * qJD(3) + t583 * t533 - t584 * t534;
t580 = -t583 * qJD(3) + t601 * t533 - t585 * t534;
t579 = t584 * qJD(3) - t585 * t533 + t602 * t534;
t576 = t557 * t452 - t556 * t463;
t571 = Ifges(3,1) * t553 + Ifges(3,4) * t554;
t570 = Ifges(3,4) * t553 + Ifges(3,2) * t554;
t569 = Ifges(3,5) * t553 + Ifges(3,6) * t554;
t537 = t569 * qJD(1);
t470 = t515 * mrSges(6,1) + t534 * t503 + t573;
t455 = mrSges(5,1) * t479 + mrSges(6,1) * t473 + mrSges(4,2) * t513 - mrSges(6,2) * t472 - mrSges(4,3) * t480 - mrSges(5,3) * t477 + pkin(4) * t470 - qJ(4) * t469 + t580 * qJD(3) + t584 * qJDD(3) - t585 * t514 + t602 * t515 + t581 * t533;
t454 = -mrSges(4,1) * t513 + mrSges(4,3) * t481 - mrSges(5,1) * t478 + mrSges(5,2) * t477 + mrSges(6,1) * t475 - mrSges(6,3) * t472 - pkin(4) * (t533 * t503 - t582) - qJ(5) * t568 - pkin(3) * t469 + t581 * t534 + t585 * t515 + (-pkin(4) * mrSges(6,1) - t601) * t514 + t583 * qJDD(3) + t579 * qJD(3);
t449 = t554 * qJD(1) * t537 + mrSges(3,2) * t531 - mrSges(3,3) * t516 - pkin(6) * t461 + t571 * qJDD(1) - t555 * t454 + t598 * t455;
t448 = -mrSges(3,1) * t531 + mrSges(3,3) * t517 - pkin(2) * t562 + pkin(6) * t574 + t570 * qJDD(1) + t598 * t454 + t555 * t455 - t537 * t589;
t447 = mrSges(2,1) * g(3) - qJ(4) * t565 - pkin(3) * (-qJD(3) * t525 + t564) + mrSges(2,3) * t539 - mrSges(3,1) * t516 + mrSges(3,2) * t517 + mrSges(5,3) * t478 - mrSges(5,2) * t479 - mrSges(4,1) * t480 + mrSges(4,2) * t481 + mrSges(6,3) * t473 - mrSges(6,2) * t475 + qJ(5) * t470 - pkin(1) * t453 - pkin(2) * t461 + (pkin(3) * t503 + t580) * t534 + (pkin(3) * mrSges(6,1) - t584) * t515 + (pkin(3) * mrSges(5,2) + t600) * qJDD(3) + (Ifges(2,6) - t569) * qJDD(1) + (-qJ(4) * t592 - t579) * t533 + (-qJ(4) * (-mrSges(5,1) - mrSges(6,1)) + t583) * t514 + (-t553 * t570 + t554 * t571 + Ifges(2,5)) * t559;
t446 = -mrSges(2,2) * g(3) - mrSges(2,3) * t538 + Ifges(2,5) * qJDD(1) - t559 * Ifges(2,6) - qJ(2) * t453 - t553 * t448 + t554 * t449;
t1 = [-m(1) * g(1) + t576; -m(1) * g(2) + t593; (-m(1) - m(2)) * g(3) + t453; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t593 + t557 * t446 - t556 * t447; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t576 + t556 * t446 + t557 * t447; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t538 - mrSges(2,2) * t539 + t553 * t449 + t554 * t448 + pkin(1) * (-mrSges(3,2) * t587 + t561) + qJ(2) * t575;];
tauB = t1;
