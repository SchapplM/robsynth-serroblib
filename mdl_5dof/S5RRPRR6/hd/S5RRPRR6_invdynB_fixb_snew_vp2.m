% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:45
% EndTime: 2020-01-03 12:05:49
% DurationCPUTime: 3.88s
% Computational Cost: add. (54801->249), mult. (75094->329), div. (0->0), fcn. (46096->10), ass. (0->111)
t506 = sin(qJ(1));
t510 = cos(qJ(1));
t490 = -g(2) * t510 - g(3) * t506;
t481 = qJDD(1) * pkin(1) + t490;
t489 = -g(2) * t506 + t510 * g(3);
t511 = qJD(1) ^ 2;
t482 = -pkin(1) * t511 + t489;
t505 = sin(qJ(2));
t509 = cos(qJ(2));
t464 = t505 * t481 + t509 * t482;
t499 = (qJD(1) + qJD(2));
t496 = t499 ^ 2;
t497 = qJDD(1) + qJDD(2);
t540 = -pkin(2) * t496 + qJ(3) * t497 + (2 * qJD(3) * t499) + t464;
t501 = sin(pkin(9));
t502 = cos(pkin(9));
t451 = -t502 * g(1) - t501 * t540;
t539 = mrSges(4,2) * t501;
t538 = mrSges(4,3) * t497;
t537 = t496 * t501 ^ 2;
t536 = t499 * t501;
t504 = sin(qJ(4));
t535 = t501 * t504;
t508 = cos(qJ(4));
t534 = t501 * t508;
t533 = t502 * t497;
t532 = t502 * t499;
t452 = -t501 * g(1) + t502 * t540;
t475 = (-mrSges(4,1) * t502 + t539) * t499;
t518 = -pkin(3) * t502 - pkin(7) * t501;
t477 = t518 * t499;
t440 = t477 * t532 + t452;
t463 = t509 * t481 - t505 * t482;
t515 = -t496 * qJ(3) + qJDD(3) - t463;
t450 = (-pkin(2) + t518) * t497 + t515;
t449 = t508 * t450;
t528 = qJD(4) * t499;
t474 = (t497 * t508 - t504 * t528) * t501;
t485 = qJDD(4) - t533;
t486 = qJD(4) - t532;
t433 = t485 * pkin(4) - t474 * pkin(8) + t449 + (-pkin(4) * t508 * t537 - pkin(8) * t486 * t536 - t440) * t504;
t436 = t440 * t508 + t504 * t450;
t525 = t499 * t534;
t472 = pkin(4) * t486 - pkin(8) * t525;
t473 = (-t497 * t504 - t508 * t528) * t501;
t527 = t504 ^ 2 * t537;
t434 = -pkin(4) * t527 + pkin(8) * t473 - t472 * t486 + t436;
t503 = sin(qJ(5));
t507 = cos(qJ(5));
t431 = t433 * t507 - t434 * t503;
t465 = (-t503 * t508 - t504 * t507) * t536;
t443 = qJD(5) * t465 + t473 * t503 + t474 * t507;
t466 = (-t503 * t504 + t507 * t508) * t536;
t453 = -mrSges(6,1) * t465 + mrSges(6,2) * t466;
t484 = qJD(5) + t486;
t458 = -mrSges(6,2) * t484 + mrSges(6,3) * t465;
t483 = qJDD(5) + t485;
t429 = m(6) * t431 + mrSges(6,1) * t483 - mrSges(6,3) * t443 - t453 * t466 + t458 * t484;
t432 = t433 * t503 + t434 * t507;
t442 = -qJD(5) * t466 + t473 * t507 - t474 * t503;
t459 = mrSges(6,1) * t484 - mrSges(6,3) * t466;
t430 = m(6) * t432 - mrSges(6,2) * t483 + mrSges(6,3) * t442 + t453 * t465 - t459 * t484;
t421 = t429 * t507 + t430 * t503;
t435 = -t504 * t440 + t449;
t526 = t499 * t535;
t469 = -mrSges(5,2) * t486 - mrSges(5,3) * t526;
t471 = (mrSges(5,1) * t504 + mrSges(5,2) * t508) * t536;
t419 = m(5) * t435 + mrSges(5,1) * t485 - mrSges(5,3) * t474 + t469 * t486 - t471 * t525 + t421;
t470 = mrSges(5,1) * t486 - mrSges(5,3) * t525;
t519 = -t429 * t503 + t430 * t507;
t420 = m(5) * t436 - mrSges(5,2) * t485 + mrSges(5,3) * t473 - t470 * t486 - t471 * t526 + t519;
t520 = -t504 * t419 + t420 * t508;
t416 = m(4) * t452 + (t475 * t499 + t538) * t502 + t520;
t439 = t477 * t536 - t451;
t437 = -pkin(4) * t473 - pkin(8) * t527 + t472 * t525 + t439;
t514 = m(6) * t437 - mrSges(6,1) * t442 + t443 * mrSges(6,2) - t465 * t458 + t466 * t459;
t512 = -m(5) * t439 + t473 * mrSges(5,1) - t474 * mrSges(5,2) - t514;
t425 = m(4) * t451 + (-t538 + (-t469 * t504 - t470 * t508 - t475) * t499) * t501 + t512;
t521 = t416 * t502 - t425 * t501;
t408 = m(3) * t464 - mrSges(3,1) * t496 - mrSges(3,2) * t497 + t521;
t417 = t508 * t419 + t504 * t420;
t456 = -pkin(2) * t497 + t515;
t513 = -m(4) * t456 + mrSges(4,1) * t533 - t417 + (t496 * t502 ^ 2 + t537) * mrSges(4,3);
t413 = m(3) * t463 - t496 * mrSges(3,2) + (mrSges(3,1) - t539) * t497 + t513;
t522 = t408 * t509 - t413 * t505;
t402 = m(2) * t489 - mrSges(2,1) * t511 - qJDD(1) * mrSges(2,2) + t522;
t404 = t408 * t505 + t413 * t509;
t403 = m(2) * t490 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t511 + t404;
t531 = t402 * t506 + t403 * t510;
t410 = t416 * t501 + t425 * t502;
t523 = -t402 * t510 + t403 * t506;
t517 = Ifges(4,1) * t501 + Ifges(4,4) * t502;
t516 = Ifges(4,5) * t501 + Ifges(4,6) * t502;
t476 = t516 * t499;
t462 = Ifges(5,5) * t486 + (Ifges(5,1) * t508 - Ifges(5,4) * t504) * t536;
t461 = Ifges(5,6) * t486 + (Ifges(5,4) * t508 - Ifges(5,2) * t504) * t536;
t460 = Ifges(5,3) * t486 + (Ifges(5,5) * t508 - Ifges(5,6) * t504) * t536;
t446 = Ifges(6,1) * t466 + Ifges(6,4) * t465 + Ifges(6,5) * t484;
t445 = Ifges(6,4) * t466 + Ifges(6,2) * t465 + Ifges(6,6) * t484;
t444 = Ifges(6,5) * t466 + Ifges(6,6) * t465 + Ifges(6,3) * t484;
t423 = mrSges(6,2) * t437 - mrSges(6,3) * t431 + Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t483 + t444 * t465 - t445 * t484;
t422 = -mrSges(6,1) * t437 + mrSges(6,3) * t432 + Ifges(6,4) * t443 + Ifges(6,2) * t442 + Ifges(6,6) * t483 - t444 * t466 + t446 * t484;
t411 = mrSges(5,2) * t439 - mrSges(5,3) * t435 + Ifges(5,1) * t474 + Ifges(5,4) * t473 + Ifges(5,5) * t485 - pkin(8) * t421 - t422 * t503 + t423 * t507 - t460 * t526 - t461 * t486;
t409 = -mrSges(5,1) * t439 + mrSges(5,3) * t436 + Ifges(5,4) * t474 + Ifges(5,2) * t473 + Ifges(5,6) * t485 - pkin(4) * t514 + pkin(8) * t519 + t507 * t422 + t503 * t423 - t460 * t525 + t486 * t462;
t405 = Ifges(4,2) * t533 - mrSges(4,1) * t456 + mrSges(4,3) * t452 - Ifges(5,5) * t474 - Ifges(5,6) * t473 - Ifges(5,3) * t485 - mrSges(5,1) * t435 + mrSges(5,2) * t436 - Ifges(6,5) * t443 - Ifges(6,6) * t442 - Ifges(6,3) * t483 - t466 * t445 + t465 * t446 - mrSges(6,1) * t431 + mrSges(6,2) * t432 - pkin(4) * t421 - pkin(3) * t417 + (Ifges(4,4) * t497 + (-t461 * t508 - t462 * t504 - t476) * t499) * t501;
t398 = mrSges(4,2) * t456 - mrSges(4,3) * t451 - pkin(7) * t417 - t504 * t409 + t508 * t411 + t476 * t532 + t497 * t517;
t397 = t496 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t464 - mrSges(4,1) * t451 + mrSges(4,2) * t452 - t504 * t411 - t508 * t409 - pkin(3) * t512 - pkin(7) * t520 - pkin(2) * t410 + (Ifges(3,6) - t516) * t497 + (-pkin(3) * (-t469 * t535 - t470 * t534) + (-t501 * (Ifges(4,4) * t501 + Ifges(4,2) * t502) + t502 * t517) * t499) * t499;
t396 = -mrSges(3,2) * g(1) - mrSges(3,3) * t463 + Ifges(3,5) * t497 - Ifges(3,6) * t496 - qJ(3) * t410 + t398 * t502 - t405 * t501;
t395 = -mrSges(2,2) * g(1) - mrSges(2,3) * t490 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t511 - pkin(6) * t404 + t396 * t509 - t397 * t505;
t394 = Ifges(2,6) * qJDD(1) + t511 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t489 + t505 * t396 + t509 * t397 - pkin(1) * (-m(3) * g(1) + t410) + pkin(6) * t522;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t410; -m(1) * g(2) + t531; -m(1) * g(3) + t523; pkin(1) * t404 + t501 * t398 + t502 * t405 + pkin(2) * (-t497 * t539 + t513) + qJ(3) * t521 - mrSges(3,2) * t464 + mrSges(3,1) * t463 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t497 + mrSges(2,1) * t490 - mrSges(2,2) * t489 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t523 + t510 * t394 + t506 * t395; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t531 + t394 * t506 - t395 * t510;];
tauB = t1;
