% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:19:03
% EndTime: 2019-05-06 09:19:09
% DurationCPUTime: 2.85s
% Computational Cost: add. (17225->305), mult. (37593->359), div. (0->0), fcn. (24288->8), ass. (0->121)
t597 = Ifges(6,1) + Ifges(7,1);
t579 = Ifges(6,4) - Ifges(7,5);
t592 = Ifges(7,4) + Ifges(6,5);
t596 = Ifges(6,2) + Ifges(7,3);
t587 = Ifges(6,6) - Ifges(7,6);
t595 = -2 * qJD(3);
t594 = -2 * qJD(4);
t593 = Ifges(4,1) + Ifges(5,1);
t580 = Ifges(4,4) - Ifges(5,5);
t591 = Ifges(4,5) + Ifges(5,4);
t590 = Ifges(4,2) + Ifges(5,3);
t589 = Ifges(5,2) + Ifges(4,3);
t588 = Ifges(4,6) - Ifges(5,6);
t586 = Ifges(6,3) + Ifges(7,2);
t542 = sin(pkin(9));
t544 = sin(qJ(2));
t568 = qJD(1) * t544;
t576 = cos(pkin(9));
t520 = -t576 * qJD(2) + t542 * t568;
t521 = t542 * qJD(2) + t576 * t568;
t543 = sin(qJ(5));
t583 = cos(qJ(5));
t482 = -t583 * t520 + t521 * t543;
t483 = t543 * t520 + t583 * t521;
t546 = cos(qJ(2));
t567 = qJD(1) * t546;
t532 = qJD(5) + t567;
t585 = t482 * t596 - t483 * t579 - t532 * t587;
t584 = -t579 * t482 + t483 * t597 + t592 * t532;
t549 = qJD(1) ^ 2;
t545 = sin(qJ(1));
t547 = cos(qJ(1));
t556 = -g(1) * t547 - g(2) * t545;
t513 = -pkin(1) * t549 + qJDD(1) * pkin(7) + t556;
t488 = -t546 * g(3) - t544 * t513;
t525 = (-pkin(2) * t546 - qJ(3) * t544) * qJD(1);
t548 = qJD(2) ^ 2;
t465 = -qJDD(2) * pkin(2) - t548 * qJ(3) + t525 * t568 + qJDD(3) - t488;
t564 = qJD(1) * qJD(2);
t561 = t546 * t564;
t527 = qJDD(1) * t544 + t561;
t498 = -t576 * qJDD(2) + t527 * t542;
t499 = t542 * qJDD(2) + t576 * t527;
t562 = t520 * t567;
t431 = t465 - (t499 + t562) * qJ(4) - (t521 * t567 - t498) * pkin(3) + t521 * t594;
t560 = t545 * g(1) - t547 * g(2);
t512 = -qJDD(1) * pkin(1) - t549 * pkin(7) - t560;
t534 = t544 * t564;
t528 = qJDD(1) * t546 - t534;
t460 = (-t527 - t561) * qJ(3) + (-t528 + t534) * pkin(2) + t512;
t489 = -g(3) * t544 + t546 * t513;
t466 = -pkin(2) * t548 + qJDD(2) * qJ(3) + t525 * t567 + t489;
t432 = t576 * t460 - t542 * t466 + t521 * t595;
t582 = -mrSges(4,3) - mrSges(5,2);
t581 = -mrSges(6,3) - mrSges(7,2);
t575 = t546 ^ 2 * t549;
t484 = pkin(3) * t520 - qJ(4) * t521;
t430 = t528 * pkin(3) - qJ(4) * t575 + t521 * t484 + qJDD(4) - t432;
t423 = (-t499 + t562) * pkin(8) + (t520 * t521 + t528) * pkin(4) + t430;
t433 = t542 * t460 + t576 * t466 + t520 * t595;
t429 = -pkin(3) * t575 - t528 * qJ(4) - t520 * t484 + t567 * t594 + t433;
t500 = pkin(4) * t567 - pkin(8) * t521;
t518 = t520 ^ 2;
t425 = -pkin(4) * t518 + pkin(8) * t498 - t500 * t567 + t429;
t419 = t543 * t423 + t583 * t425;
t574 = t587 * t482 - t592 * t483 - t586 * t532;
t457 = mrSges(7,1) * t482 - mrSges(7,3) * t483;
t573 = -mrSges(6,1) * t482 - mrSges(6,2) * t483 - t457;
t572 = t588 * t520 - t591 * t521 + t589 * t567;
t571 = t590 * t520 - t580 * t521 + t588 * t567;
t570 = t580 * t520 - t593 * t521 + t591 * t567;
t485 = mrSges(5,1) * t520 - mrSges(5,3) * t521;
t569 = -mrSges(4,1) * t520 - mrSges(4,2) * t521 - t485;
t456 = pkin(5) * t482 - qJ(6) * t483;
t524 = qJDD(5) + t528;
t531 = t532 ^ 2;
t415 = -pkin(5) * t531 + qJ(6) * t524 + 0.2e1 * qJD(6) * t532 - t456 * t482 + t419;
t469 = -mrSges(7,1) * t532 + mrSges(7,2) * t483;
t563 = m(7) * t415 + t524 * mrSges(7,3) + t532 * t469;
t496 = -mrSges(4,1) * t567 - mrSges(4,3) * t521;
t497 = mrSges(5,1) * t567 + mrSges(5,2) * t521;
t444 = qJD(5) * t483 - t583 * t498 + t499 * t543;
t468 = mrSges(6,1) * t532 - mrSges(6,3) * t483;
t406 = m(6) * t419 - t524 * mrSges(6,2) + t581 * t444 - t532 * t468 + t573 * t482 + t563;
t418 = t583 * t423 - t543 * t425;
t445 = -t482 * qJD(5) + t543 * t498 + t583 * t499;
t467 = -mrSges(6,2) * t532 - mrSges(6,3) * t482;
t416 = -t524 * pkin(5) - t531 * qJ(6) + t483 * t456 + qJDD(6) - t418;
t470 = -mrSges(7,2) * t482 + mrSges(7,3) * t532;
t557 = -m(7) * t416 + t524 * mrSges(7,1) + t532 * t470;
t407 = m(6) * t418 + t524 * mrSges(6,1) + t581 * t445 + t532 * t467 + t573 * t483 + t557;
t558 = t583 * t406 - t543 * t407;
t554 = m(5) * t429 - t528 * mrSges(5,3) + t558;
t398 = m(4) * t433 + t528 * mrSges(4,2) + t569 * t520 + t582 * t498 + (t496 - t497) * t567 + t554;
t494 = -mrSges(5,2) * t520 - mrSges(5,3) * t567;
t495 = mrSges(4,2) * t567 - mrSges(4,3) * t520;
t403 = t543 * t406 + t583 * t407;
t553 = m(5) * t430 + t528 * mrSges(5,1) + t403;
t399 = m(4) * t432 - t528 * mrSges(4,1) + t569 * t521 + t582 * t499 + (-t494 - t495) * t567 - t553;
t559 = t576 * t398 - t399 * t542;
t427 = -t498 * pkin(4) - t518 * pkin(8) + t521 * t500 - t431;
t421 = t427 + (t482 * t532 - t445) * qJ(6) - 0.2e1 * qJD(6) * t483 + (t483 * t532 + t444) * pkin(5);
t412 = m(7) * t421 + t444 * mrSges(7,1) - t445 * mrSges(7,3) - t483 * t469 + t482 * t470;
t396 = t542 * t398 + t576 * t399;
t552 = -m(6) * t427 - t444 * mrSges(6,1) - t445 * mrSges(6,2) - t482 * t467 - t483 * t468 - t412;
t408 = m(5) * t431 + t498 * mrSges(5,1) - t499 * mrSges(5,3) + t520 * t494 - t521 * t497 + t552;
t411 = t445 * mrSges(7,2) + t483 * t457 - t557;
t550 = mrSges(6,1) * t418 - mrSges(7,1) * t416 - mrSges(6,2) * t419 + mrSges(7,3) * t415 - pkin(5) * t411 + qJ(6) * t563 + t586 * t524 - t585 * t483 + (-qJ(6) * t457 + t584) * t482 + t592 * t445 + (-qJ(6) * mrSges(7,2) - t587) * t444;
t404 = m(4) * t465 + t498 * mrSges(4,1) + t499 * mrSges(4,2) + t520 * t495 + t521 * t496 + t408;
t530 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t567;
t529 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t568;
t526 = (-mrSges(3,1) * t546 + mrSges(3,2) * t544) * qJD(1);
t511 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t544 + Ifges(3,4) * t546) * qJD(1);
t510 = Ifges(3,6) * qJD(2) + (t544 * Ifges(3,4) + Ifges(3,2) * t546) * qJD(1);
t509 = Ifges(3,3) * qJD(2) + (t544 * Ifges(3,5) + t546 * Ifges(3,6)) * qJD(1);
t402 = mrSges(6,2) * t427 + mrSges(7,2) * t416 - mrSges(6,3) * t418 - mrSges(7,3) * t421 - qJ(6) * t412 - t579 * t444 + t445 * t597 + t574 * t482 + t592 * t524 + t585 * t532;
t401 = -mrSges(6,1) * t427 - mrSges(7,1) * t421 + mrSges(7,2) * t415 + mrSges(6,3) * t419 - pkin(5) * t412 - t444 * t596 + t579 * t445 + t574 * t483 + t587 * t524 + t584 * t532;
t400 = t499 * mrSges(5,2) + t521 * t485 + t494 * t567 + t553;
t395 = mrSges(4,2) * t465 + mrSges(5,2) * t430 - mrSges(4,3) * t432 - mrSges(5,3) * t431 - pkin(8) * t403 - qJ(4) * t408 - t543 * t401 + t583 * t402 - t580 * t498 + t593 * t499 + t572 * t520 - t528 * t591 - t571 * t567;
t394 = -mrSges(4,1) * t465 - mrSges(5,1) * t431 + mrSges(5,2) * t429 + mrSges(4,3) * t433 - pkin(3) * t408 - pkin(4) * t552 - pkin(8) * t558 - t583 * t401 - t543 * t402 - t590 * t498 + t580 * t499 + t572 * t521 - t528 * t588 + t570 * t567;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t560 - mrSges(2,2) * t556 + t544 * (mrSges(3,2) * t512 - mrSges(3,3) * t488 + Ifges(3,1) * t527 + Ifges(3,4) * t528 + Ifges(3,5) * qJDD(2) - qJ(3) * t396 - qJD(2) * t510 - t542 * t394 + t576 * t395 + t509 * t567) + t546 * (t550 + Ifges(3,6) * qJDD(2) - qJ(4) * (-t497 * t567 + t554) - t591 * t499 + t571 * t521 + (Ifges(3,2) + t589) * t528 - t509 * t568 + (mrSges(5,2) * qJ(4) + t588) * t498 + Ifges(3,4) * t527 + qJD(2) * t511 - mrSges(3,1) * t512 + mrSges(3,3) * t489 + mrSges(4,2) * t433 + (qJ(4) * t485 + t570) * t520 - mrSges(5,3) * t429 + mrSges(5,1) * t430 - mrSges(4,1) * t432 + pkin(4) * t403 + pkin(3) * t400 - pkin(2) * t396) + pkin(1) * (-m(3) * t512 + t528 * mrSges(3,1) - t527 * mrSges(3,2) + (-t529 * t544 + t530 * t546) * qJD(1) - t396) + pkin(7) * (t546 * (m(3) * t489 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t528 - qJD(2) * t529 + t526 * t567 + t559) - t544 * (m(3) * t488 + qJDD(2) * mrSges(3,1) - t527 * mrSges(3,3) + qJD(2) * t530 - t526 * t568 - t404)); Ifges(3,5) * t527 + Ifges(3,6) * t528 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t488 - mrSges(3,2) * t489 + t542 * t395 + t576 * t394 - pkin(2) * t404 + qJ(3) * t559 + (t544 * t510 - t546 * t511) * qJD(1); t404; t400; t550; t411;];
tauJ  = t1;
