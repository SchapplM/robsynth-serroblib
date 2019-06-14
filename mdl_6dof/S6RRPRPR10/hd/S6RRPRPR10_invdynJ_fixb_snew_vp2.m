% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:38:43
% EndTime: 2019-05-06 15:38:56
% DurationCPUTime: 6.67s
% Computational Cost: add. (82691->342), mult. (188380->428), div. (0->0), fcn. (150239->12), ass. (0->145)
t596 = Ifges(5,1) + Ifges(6,2);
t588 = Ifges(5,4) + Ifges(6,6);
t587 = Ifges(5,5) - Ifges(6,4);
t595 = -Ifges(5,2) - Ifges(6,3);
t586 = Ifges(5,6) - Ifges(6,5);
t594 = Ifges(5,3) + Ifges(6,1);
t550 = cos(pkin(6));
t544 = qJD(1) * t550 + qJD(2);
t547 = sin(pkin(11));
t549 = cos(pkin(11));
t553 = sin(qJ(2));
t548 = sin(pkin(6));
t575 = qJD(1) * t548;
t572 = t553 * t575;
t524 = t544 * t549 - t547 * t572;
t525 = t544 * t547 + t549 * t572;
t552 = sin(qJ(4));
t591 = cos(qJ(4));
t504 = -t524 * t591 + t552 * t525;
t556 = cos(qJ(2));
t574 = qJD(1) * t556;
t571 = t548 * t574;
t539 = -qJD(4) + t571;
t488 = mrSges(6,1) * t504 + mrSges(6,3) * t539;
t573 = qJDD(1) * t548;
t537 = -qJD(2) * t572 + t556 * t573;
t529 = qJDD(4) - t537;
t534 = (-pkin(2) * t556 - qJ(3) * t553) * t575;
t542 = t544 ^ 2;
t543 = qJDD(1) * t550 + qJDD(2);
t558 = qJD(1) ^ 2;
t554 = sin(qJ(1));
t557 = cos(qJ(1));
t567 = -g(1) * t557 - g(2) * t554;
t533 = -pkin(1) * t558 + pkin(8) * t573 + t567;
t570 = t554 * g(1) - g(2) * t557;
t590 = pkin(8) * t548;
t532 = qJDD(1) * pkin(1) + t558 * t590 + t570;
t583 = t532 * t550;
t576 = t556 * t533 + t553 * t583;
t483 = -t542 * pkin(2) + t543 * qJ(3) + (-g(3) * t553 + t534 * t574) * t548 + t576;
t536 = (qJD(2) * t574 + qJDD(1) * t553) * t548;
t589 = t550 * g(3);
t484 = -t537 * pkin(2) - t589 - t536 * qJ(3) + (-t532 + (pkin(2) * t553 - qJ(3) * t556) * t544 * qJD(1)) * t548;
t442 = -0.2e1 * qJD(3) * t525 - t547 * t483 + t549 * t484;
t512 = t536 * t549 + t543 * t547;
t437 = (-t524 * t571 - t512) * pkin(9) + (t524 * t525 - t537) * pkin(3) + t442;
t443 = 0.2e1 * qJD(3) * t524 + t549 * t483 + t547 * t484;
t511 = -t536 * t547 + t543 * t549;
t513 = -pkin(3) * t571 - pkin(9) * t525;
t523 = t524 ^ 2;
t440 = -pkin(3) * t523 + pkin(9) * t511 + t513 * t571 + t443;
t432 = t437 * t591 - t552 * t440;
t505 = t552 * t524 + t525 * t591;
t476 = pkin(4) * t504 - qJ(5) * t505;
t538 = t539 ^ 2;
t431 = -t529 * pkin(4) - t538 * qJ(5) + t505 * t476 + qJDD(5) - t432;
t464 = -t504 * qJD(4) + t552 * t511 + t512 * t591;
t584 = t504 * t539;
t426 = (t504 * t505 - t529) * pkin(10) + (t464 - t584) * pkin(5) + t431;
t463 = t505 * qJD(4) - t511 * t591 + t552 * t512;
t492 = pkin(5) * t505 + pkin(10) * t539;
t503 = t504 ^ 2;
t581 = t548 * t556;
t500 = -g(3) * t581 - t553 * t533 + t556 * t583;
t482 = -t543 * pkin(2) - t542 * qJ(3) + t534 * t572 + qJDD(3) - t500;
t449 = -t511 * pkin(3) - t523 * pkin(9) + t525 * t513 + t482;
t592 = -2 * qJD(5);
t559 = (-t464 - t584) * qJ(5) + t449 + (-t539 * pkin(4) + t592) * t505;
t429 = t559 + (pkin(4) + pkin(10)) * t463 - t505 * t492 - t503 * pkin(5);
t551 = sin(qJ(6));
t555 = cos(qJ(6));
t424 = t426 * t555 - t429 * t551;
t486 = t504 * t555 + t539 * t551;
t448 = qJD(6) * t486 + t463 * t551 + t529 * t555;
t487 = t504 * t551 - t539 * t555;
t456 = -mrSges(7,1) * t486 + mrSges(7,2) * t487;
t462 = qJDD(6) + t464;
t502 = qJD(6) + t505;
t465 = -mrSges(7,2) * t502 + mrSges(7,3) * t486;
t421 = m(7) * t424 + mrSges(7,1) * t462 - mrSges(7,3) * t448 - t456 * t487 + t465 * t502;
t425 = t426 * t551 + t429 * t555;
t447 = -qJD(6) * t487 + t463 * t555 - t529 * t551;
t466 = mrSges(7,1) * t502 - mrSges(7,3) * t487;
t422 = m(7) * t425 - mrSges(7,2) * t462 + mrSges(7,3) * t447 + t456 * t486 - t466 * t502;
t412 = t555 * t421 + t551 * t422;
t478 = -mrSges(6,2) * t504 - mrSges(6,3) * t505;
t565 = -m(6) * t431 - t464 * mrSges(6,1) - t505 * t478 - t412;
t410 = t529 * mrSges(6,2) - t539 * t488 - t565;
t433 = t552 * t437 + t591 * t440;
t564 = -t538 * pkin(4) + t529 * qJ(5) - t504 * t476 + t433;
t428 = -t463 * pkin(5) - t503 * pkin(10) + (t592 - t492) * t539 + t564;
t450 = Ifges(7,5) * t487 + Ifges(7,6) * t486 + Ifges(7,3) * t502;
t452 = Ifges(7,1) * t487 + Ifges(7,4) * t486 + Ifges(7,5) * t502;
t413 = -mrSges(7,1) * t428 + mrSges(7,3) * t425 + Ifges(7,4) * t448 + Ifges(7,2) * t447 + Ifges(7,6) * t462 - t450 * t487 + t452 * t502;
t451 = Ifges(7,4) * t487 + Ifges(7,2) * t486 + Ifges(7,6) * t502;
t414 = mrSges(7,2) * t428 - mrSges(7,3) * t424 + Ifges(7,1) * t448 + Ifges(7,4) * t447 + Ifges(7,5) * t462 + t450 * t486 - t451 * t502;
t430 = 0.2e1 * qJD(5) * t539 - t564;
t489 = mrSges(6,1) * t505 - mrSges(6,2) * t539;
t566 = -m(7) * t428 + t447 * mrSges(7,1) - t448 * mrSges(7,2) + t486 * t465 - t487 * t466;
t562 = -m(6) * t430 + t529 * mrSges(6,3) - t539 * t489 - t566;
t577 = -t588 * t504 + t596 * t505 - t587 * t539;
t578 = t595 * t504 + t588 * t505 - t586 * t539;
t593 = t586 * t463 - t587 * t464 - t577 * t504 - t578 * t505 - t594 * t529 - mrSges(5,1) * t432 + mrSges(5,2) * t433 - mrSges(6,2) * t431 + mrSges(6,3) * t430 + pkin(4) * t410 + pkin(10) * t412 - qJ(5) * (-t463 * mrSges(6,1) - t504 * t478 + t562) + t551 * t413 - t555 * t414;
t582 = t548 * t553;
t477 = mrSges(5,1) * t504 + mrSges(5,2) * t505;
t490 = mrSges(5,2) * t539 - mrSges(5,3) * t504;
t407 = m(5) * t432 - t464 * mrSges(5,3) - t505 * t477 + (t488 - t490) * t539 + (mrSges(5,1) - mrSges(6,2)) * t529 + t565;
t491 = -mrSges(5,1) * t539 - mrSges(5,3) * t505;
t417 = m(5) * t433 - t529 * mrSges(5,2) + t539 * t491 + (-t477 - t478) * t504 + (-mrSges(5,3) - mrSges(6,1)) * t463 + t562;
t405 = t591 * t407 + t552 * t417;
t506 = -mrSges(4,1) * t524 + mrSges(4,2) * t525;
t509 = mrSges(4,2) * t571 + mrSges(4,3) * t524;
t403 = m(4) * t442 - mrSges(4,1) * t537 - mrSges(4,3) * t512 - t506 * t525 - t509 * t571 + t405;
t510 = -mrSges(4,1) * t571 - mrSges(4,3) * t525;
t568 = -t407 * t552 + t591 * t417;
t404 = m(4) * t443 + mrSges(4,2) * t537 + mrSges(4,3) * t511 + t506 * t524 + t510 * t571 + t568;
t398 = t549 * t403 + t547 * t404;
t580 = -t551 * t421 + t555 * t422;
t579 = t586 * t504 - t587 * t505 + t594 * t539;
t569 = -t403 * t547 + t549 * t404;
t435 = t463 * pkin(4) + t559;
t411 = m(6) * t435 - t463 * mrSges(6,2) - t464 * mrSges(6,3) - t504 * t488 - t505 * t489 + t580;
t563 = mrSges(7,1) * t424 - mrSges(7,2) * t425 + Ifges(7,5) * t448 + Ifges(7,6) * t447 + Ifges(7,3) * t462 + t487 * t451 - t486 * t452;
t561 = m(5) * t449 + t463 * mrSges(5,1) + t464 * mrSges(5,2) + t504 * t490 + t505 * t491 + t411;
t409 = m(4) * t482 - t511 * mrSges(4,1) + t512 * mrSges(4,2) - t524 * t509 + t525 * t510 + t561;
t535 = (-mrSges(3,1) * t556 + mrSges(3,2) * t553) * t575;
t531 = -mrSges(3,2) * t544 + mrSges(3,3) * t571;
t530 = mrSges(3,1) * t544 - mrSges(3,3) * t572;
t517 = -t548 * t532 - t589;
t516 = Ifges(3,5) * t544 + (Ifges(3,1) * t553 + Ifges(3,4) * t556) * t575;
t515 = Ifges(3,6) * t544 + (t553 * Ifges(3,4) + Ifges(3,2) * t556) * t575;
t514 = Ifges(3,3) * t544 + (Ifges(3,5) * t553 + Ifges(3,6) * t556) * t575;
t501 = -g(3) * t582 + t576;
t496 = Ifges(4,1) * t525 + Ifges(4,4) * t524 - Ifges(4,5) * t571;
t495 = Ifges(4,4) * t525 + Ifges(4,2) * t524 - Ifges(4,6) * t571;
t494 = Ifges(4,5) * t525 + Ifges(4,6) * t524 - Ifges(4,3) * t571;
t408 = m(3) * t500 + t543 * mrSges(3,1) - t536 * mrSges(3,3) + t544 * t531 - t535 * t572 - t409;
t399 = mrSges(6,1) * t431 + mrSges(5,2) * t449 - mrSges(5,3) * t432 - mrSges(6,3) * t435 + pkin(5) * t412 - qJ(5) * t411 - t588 * t463 + t596 * t464 + t579 * t504 + t587 * t529 + t578 * t539 + t563;
t397 = m(3) * t501 - mrSges(3,2) * t543 + mrSges(3,3) * t537 - t530 * t544 + t535 * t571 + t569;
t396 = -mrSges(5,1) * t449 - mrSges(6,1) * t430 + mrSges(6,2) * t435 + mrSges(5,3) * t433 - pkin(4) * t411 - pkin(5) * t566 - pkin(10) * t580 - t555 * t413 - t551 * t414 + t595 * t463 + t588 * t464 + t579 * t505 + t586 * t529 - t577 * t539;
t395 = mrSges(4,2) * t482 - mrSges(4,3) * t442 + Ifges(4,1) * t512 + Ifges(4,4) * t511 - Ifges(4,5) * t537 - pkin(9) * t405 - t552 * t396 + t399 * t591 + t524 * t494 + t495 * t571;
t394 = -mrSges(4,1) * t482 + mrSges(4,3) * t443 + Ifges(4,4) * t512 + Ifges(4,2) * t511 - Ifges(4,6) * t537 - pkin(3) * t561 + pkin(9) * t568 + t396 * t591 + t552 * t399 - t525 * t494 - t496 * t571;
t393 = Ifges(3,5) * t536 + Ifges(3,6) * t537 + Ifges(3,3) * t543 + mrSges(3,1) * t500 - mrSges(3,2) * t501 + t547 * t395 + t549 * t394 - pkin(2) * t409 + qJ(3) * t569 + (t515 * t553 - t516 * t556) * t575;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t570 - mrSges(2,2) * t567 + (mrSges(3,2) * t517 - mrSges(3,3) * t500 + Ifges(3,1) * t536 + Ifges(3,4) * t537 + Ifges(3,5) * t543 - qJ(3) * t398 - t394 * t547 + t395 * t549 + t514 * t571 - t515 * t544) * t582 + (-t514 * t572 - pkin(2) * t398 + t593 + (Ifges(3,2) + Ifges(4,3)) * t537 + Ifges(3,6) * t543 + t544 * t516 + Ifges(3,4) * t536 - t525 * t495 - Ifges(4,5) * t512 - mrSges(3,1) * t517 + t524 * t496 - Ifges(4,6) * t511 + mrSges(3,3) * t501 - mrSges(4,1) * t442 + mrSges(4,2) * t443 - pkin(3) * t405) * t581 + t550 * t393 + pkin(1) * ((t397 * t553 + t408 * t556) * t550 + (-m(3) * t517 + t537 * mrSges(3,1) - t536 * mrSges(3,2) + (-t530 * t553 + t531 * t556) * t575 - t398) * t548) + (t397 * t556 - t408 * t553) * t590; t393; t409; -t593; t410; t563;];
tauJ  = t1;
