% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 08:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:29:39
% EndTime: 2019-05-07 08:29:46
% DurationCPUTime: 3.13s
% Computational Cost: add. (23289->308), mult. (45984->360), div. (0->0), fcn. (30413->8), ass. (0->124)
t579 = Ifges(6,4) + Ifges(7,4);
t596 = Ifges(6,2) + Ifges(7,2);
t590 = Ifges(6,6) + Ifges(7,6);
t542 = sin(qJ(3));
t543 = sin(qJ(2));
t568 = qJD(1) * t543;
t583 = cos(qJ(3));
t523 = t542 * qJD(2) + t583 * t568;
t546 = cos(qJ(2));
t566 = qJD(1) * qJD(2);
t562 = t546 * t566;
t526 = qJDD(1) * t543 + t562;
t486 = qJD(3) * t523 - t583 * qJDD(2) + t526 * t542;
t522 = -t583 * qJD(2) + t542 * t568;
t487 = -t522 * qJD(3) + t542 * qJDD(2) + t583 * t526;
t541 = sin(qJ(5));
t545 = cos(qJ(5));
t489 = t522 * t545 - t523 * t541;
t443 = qJD(5) * t489 + t486 * t541 + t487 * t545;
t490 = t522 * t541 + t523 * t545;
t457 = -mrSges(7,1) * t489 + mrSges(7,2) * t490;
t549 = qJD(1) ^ 2;
t544 = sin(qJ(1));
t547 = cos(qJ(1));
t561 = t544 * g(1) - t547 * g(2);
t514 = -qJDD(1) * pkin(1) - t549 * pkin(7) - t561;
t563 = t543 * t566;
t527 = t546 * qJDD(1) - t563;
t460 = (-t526 - t562) * pkin(8) + (-t527 + t563) * pkin(2) + t514;
t557 = -g(1) * t547 - g(2) * t544;
t515 = -pkin(1) * t549 + qJDD(1) * pkin(7) + t557;
t501 = -g(3) * t543 + t546 * t515;
t525 = (-pkin(2) * t546 - pkin(8) * t543) * qJD(1);
t548 = qJD(2) ^ 2;
t567 = qJD(1) * t546;
t465 = -pkin(2) * t548 + qJDD(2) * pkin(8) + t525 * t567 + t501;
t444 = t583 * t460 - t542 * t465;
t493 = pkin(3) * t522 - qJ(4) * t523;
t521 = qJDD(3) - t527;
t533 = qJD(3) - t567;
t532 = t533 ^ 2;
t430 = -t521 * pkin(3) - t532 * qJ(4) + t523 * t493 + qJDD(4) - t444;
t575 = t522 * t533;
t422 = (-t487 - t575) * pkin(9) + (t522 * t523 - t521) * pkin(4) + t430;
t445 = t542 * t460 + t583 * t465;
t584 = 2 * qJD(4);
t428 = -pkin(3) * t532 + t521 * qJ(4) - t522 * t493 + t533 * t584 + t445;
t502 = -pkin(4) * t533 - pkin(9) * t523;
t520 = t522 ^ 2;
t425 = -pkin(4) * t520 + pkin(9) * t486 + t502 * t533 + t428;
t416 = t545 * t422 - t541 * t425;
t519 = qJDD(5) - t521;
t531 = qJD(5) - t533;
t412 = -0.2e1 * qJD(6) * t490 + (t489 * t531 - t443) * qJ(6) + (t489 * t490 + t519) * pkin(5) + t416;
t466 = -mrSges(7,2) * t531 + mrSges(7,3) * t489;
t565 = m(7) * t412 + t519 * mrSges(7,1) + t531 * t466;
t408 = -t443 * mrSges(7,3) - t490 * t457 + t565;
t417 = t541 * t422 + t545 * t425;
t442 = -qJD(5) * t490 + t486 * t545 - t487 * t541;
t468 = pkin(5) * t531 - qJ(6) * t490;
t488 = t489 ^ 2;
t414 = -pkin(5) * t488 + qJ(6) * t442 + 0.2e1 * qJD(6) * t489 - t468 * t531 + t417;
t592 = Ifges(6,5) + Ifges(7,5);
t593 = Ifges(6,1) + Ifges(7,1);
t573 = t579 * t489 + t593 * t490 + t592 * t531;
t587 = t596 * t489 + t579 * t490 + t590 * t531;
t588 = Ifges(6,3) + Ifges(7,3);
t595 = mrSges(6,1) * t416 + mrSges(7,1) * t412 - mrSges(6,2) * t417 - mrSges(7,2) * t414 + pkin(5) * t408 + t590 * t442 + t592 * t443 - t573 * t489 + t587 * t490 + t588 * t519;
t594 = Ifges(4,1) + Ifges(5,1);
t580 = Ifges(4,4) - Ifges(5,5);
t578 = Ifges(4,5) + Ifges(5,4);
t591 = Ifges(4,2) + Ifges(5,3);
t577 = Ifges(4,6) - Ifges(5,6);
t589 = Ifges(4,3) + Ifges(5,2);
t494 = mrSges(5,1) * t522 - mrSges(5,3) * t523;
t458 = -mrSges(6,1) * t489 + mrSges(6,2) * t490;
t467 = -mrSges(6,2) * t531 + mrSges(6,3) * t489;
t403 = m(6) * t416 + t519 * mrSges(6,1) + t531 * t467 + (-t457 - t458) * t490 + (-mrSges(6,3) - mrSges(7,3)) * t443 + t565;
t469 = mrSges(7,1) * t531 - mrSges(7,3) * t490;
t470 = mrSges(6,1) * t531 - mrSges(6,3) * t490;
t564 = m(7) * t414 + t442 * mrSges(7,3) + t489 * t457;
t405 = m(6) * t417 + t442 * mrSges(6,3) + t489 * t458 + (-t469 - t470) * t531 + (-mrSges(6,2) - mrSges(7,2)) * t519 + t564;
t402 = t545 * t403 + t541 * t405;
t499 = -mrSges(5,2) * t522 + mrSges(5,3) * t533;
t554 = -m(5) * t430 + t521 * mrSges(5,1) + t533 * t499 - t402;
t400 = t487 * mrSges(5,2) + t523 * t494 - t554;
t498 = -mrSges(5,1) * t533 + mrSges(5,2) * t523;
t558 = -t541 * t403 + t545 * t405;
t556 = m(5) * t428 + t521 * mrSges(5,3) + t533 * t498 + t558;
t570 = -t580 * t522 + t594 * t523 + t578 * t533;
t572 = t591 * t522 - t580 * t523 - t577 * t533;
t586 = -t577 * t486 + t578 * t487 + t589 * t521 + t570 * t522 - t572 * t523 + mrSges(4,1) * t444 - mrSges(5,1) * t430 - mrSges(4,2) * t445 + mrSges(5,3) * t428 - pkin(3) * t400 - pkin(4) * t402 + qJ(4) * (-t486 * mrSges(5,2) - t522 * t494 + t556) - t595;
t582 = pkin(3) * t533;
t581 = -mrSges(4,3) - mrSges(5,2);
t574 = -t590 * t489 - t592 * t490 - t588 * t531;
t571 = t577 * t522 - t578 * t523 - t589 * t533;
t569 = -mrSges(4,1) * t522 - mrSges(4,2) * t523 - t494;
t500 = -t546 * g(3) - t543 * t515;
t497 = mrSges(4,1) * t533 - mrSges(4,3) * t523;
t397 = m(4) * t445 - t521 * mrSges(4,2) + t581 * t486 - t533 * t497 + t569 * t522 + t556;
t496 = -mrSges(4,2) * t533 - mrSges(4,3) * t522;
t398 = m(4) * t444 + t521 * mrSges(4,1) + t581 * t487 + t533 * t496 + t569 * t523 + t554;
t559 = t583 * t397 - t398 * t542;
t464 = -qJDD(2) * pkin(2) - t548 * pkin(8) + t525 * t568 - t500;
t555 = t486 * pkin(3) + t464 + (-t487 + t575) * qJ(4);
t426 = -pkin(4) * t486 - pkin(9) * t520 - t555 + (t502 - t582 + t584) * t523;
t419 = -pkin(5) * t442 - qJ(6) * t488 + t468 * t490 + qJDD(6) + t426;
t409 = m(7) * t419 - t442 * mrSges(7,1) + t443 * mrSges(7,2) - t489 * t466 + t490 * t469;
t394 = t542 * t397 + t583 * t398;
t553 = -m(6) * t426 + t442 * mrSges(6,1) - t443 * mrSges(6,2) + t489 * t467 - t490 * t470 - t409;
t429 = (-(2 * qJD(4)) + t582) * t523 + t555;
t406 = m(5) * t429 + t486 * mrSges(5,1) - t487 * mrSges(5,3) - t523 * t498 + t522 * t499 + t553;
t551 = -m(4) * t464 - t486 * mrSges(4,1) - t487 * mrSges(4,2) - t522 * t496 - t523 * t497 - t406;
t529 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t567;
t528 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t568;
t524 = (-mrSges(3,1) * t546 + mrSges(3,2) * t543) * qJD(1);
t513 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t543 + Ifges(3,4) * t546) * qJD(1);
t512 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t543 + Ifges(3,2) * t546) * qJD(1);
t511 = Ifges(3,3) * qJD(2) + (t543 * Ifges(3,5) + t546 * Ifges(3,6)) * qJD(1);
t401 = mrSges(6,2) * t426 + mrSges(7,2) * t419 - mrSges(6,3) * t416 - mrSges(7,3) * t412 - qJ(6) * t408 + t579 * t442 + t593 * t443 - t574 * t489 + t592 * t519 - t587 * t531;
t395 = -mrSges(6,1) * t426 + mrSges(6,3) * t417 - mrSges(7,1) * t419 + mrSges(7,3) * t414 - pkin(5) * t409 + qJ(6) * t564 + (-qJ(6) * t469 + t573) * t531 + (-mrSges(7,2) * qJ(6) + t590) * t519 + t574 * t490 + t579 * t443 + t596 * t442;
t393 = mrSges(4,2) * t464 + mrSges(5,2) * t430 - mrSges(4,3) * t444 - mrSges(5,3) * t429 - pkin(9) * t402 - qJ(4) * t406 - t395 * t541 + t401 * t545 - t580 * t486 + t594 * t487 + t578 * t521 + t571 * t522 + t572 * t533;
t392 = -mrSges(4,1) * t464 - mrSges(5,1) * t429 + mrSges(5,2) * t428 + mrSges(4,3) * t445 - pkin(3) * t406 - pkin(4) * t553 - pkin(9) * t558 - t545 * t395 - t541 * t401 - t591 * t486 + t580 * t487 + t577 * t521 + t571 * t523 + t570 * t533;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t561 - mrSges(2,2) * t557 + t543 * (mrSges(3,2) * t514 - mrSges(3,3) * t500 + Ifges(3,1) * t526 + Ifges(3,4) * t527 + Ifges(3,5) * qJDD(2) - pkin(8) * t394 - qJD(2) * t512 - t542 * t392 + t583 * t393 + t511 * t567) + t546 * (-mrSges(3,1) * t514 + mrSges(3,3) * t501 + Ifges(3,4) * t526 + Ifges(3,2) * t527 + Ifges(3,6) * qJDD(2) - pkin(2) * t394 + qJD(2) * t513 - t511 * t568 - t586) + pkin(1) * (-m(3) * t514 + t527 * mrSges(3,1) - t526 * mrSges(3,2) + (-t528 * t543 + t529 * t546) * qJD(1) - t394) + pkin(7) * (t546 * (m(3) * t501 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t527 - qJD(2) * t528 + t524 * t567 + t559) - t543 * (m(3) * t500 + qJDD(2) * mrSges(3,1) - t526 * mrSges(3,3) + qJD(2) * t529 - t524 * t568 + t551)); Ifges(3,5) * t526 + Ifges(3,6) * t527 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t500 - mrSges(3,2) * t501 + t542 * t393 + t583 * t392 + pkin(2) * t551 + pkin(8) * t559 + (t543 * t512 - t546 * t513) * qJD(1); t586; t400; t595; t409;];
tauJ  = t1;
