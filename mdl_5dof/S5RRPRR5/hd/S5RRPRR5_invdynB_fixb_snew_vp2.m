% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR5
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:52
% EndTime: 2019-12-05 18:33:57
% DurationCPUTime: 4.65s
% Computational Cost: add. (77172->250), mult. (106816->313), div. (0->0), fcn. (71884->10), ass. (0->109)
t514 = qJD(1) + qJD(2);
t508 = t514 ^ 2;
t516 = cos(pkin(9));
t549 = pkin(3) * t516;
t515 = sin(pkin(9));
t548 = mrSges(4,2) * t515;
t512 = t516 ^ 2;
t547 = t508 * t512;
t510 = qJDD(1) + qJDD(2);
t546 = t510 * t516;
t534 = Ifges(4,5) * t515 + Ifges(4,6) * t516;
t545 = t508 * t534;
t520 = sin(qJ(1));
t524 = cos(qJ(1));
t501 = t524 * g(2) + t520 * g(3);
t496 = qJDD(1) * pkin(1) + t501;
t500 = t520 * g(2) - t524 * g(3);
t525 = qJD(1) ^ 2;
t497 = -t525 * pkin(1) + t500;
t519 = sin(qJ(2));
t523 = cos(qJ(2));
t484 = t519 * t496 + t523 * t497;
t482 = -t508 * pkin(2) + t510 * qJ(3) + t484;
t544 = qJD(3) * t514;
t542 = -t516 * g(1) - 0.2e1 * t515 * t544;
t463 = (-pkin(7) * t510 + t508 * t549 - t482) * t515 + t542;
t467 = -t515 * g(1) + (t482 + 0.2e1 * t544) * t516;
t464 = -pkin(3) * t547 + pkin(7) * t546 + t467;
t518 = sin(qJ(4));
t522 = cos(qJ(4));
t448 = t522 * t463 - t518 * t464;
t529 = t515 * t522 + t516 * t518;
t528 = -t515 * t518 + t516 * t522;
t489 = t528 * t514;
t543 = t489 * qJD(4);
t481 = t529 * t510 + t543;
t490 = t529 * t514;
t444 = (-t481 + t543) * pkin(8) + (t489 * t490 + qJDD(4)) * pkin(4) + t448;
t449 = t518 * t463 + t522 * t464;
t480 = -t490 * qJD(4) + t528 * t510;
t487 = qJD(4) * pkin(4) - t490 * pkin(8);
t488 = t489 ^ 2;
t445 = -t488 * pkin(4) + t480 * pkin(8) - qJD(4) * t487 + t449;
t517 = sin(qJ(5));
t521 = cos(qJ(5));
t442 = t521 * t444 - t517 * t445;
t473 = t521 * t489 - t517 * t490;
t453 = t473 * qJD(5) + t517 * t480 + t521 * t481;
t474 = t517 * t489 + t521 * t490;
t459 = -t473 * mrSges(6,1) + t474 * mrSges(6,2);
t513 = qJD(4) + qJD(5);
t468 = -t513 * mrSges(6,2) + t473 * mrSges(6,3);
t509 = qJDD(4) + qJDD(5);
t440 = m(6) * t442 + t509 * mrSges(6,1) - t453 * mrSges(6,3) - t474 * t459 + t513 * t468;
t443 = t517 * t444 + t521 * t445;
t452 = -t474 * qJD(5) + t521 * t480 - t517 * t481;
t469 = t513 * mrSges(6,1) - t474 * mrSges(6,3);
t441 = m(6) * t443 - t509 * mrSges(6,2) + t452 * mrSges(6,3) + t473 * t459 - t513 * t469;
t432 = t521 * t440 + t517 * t441;
t478 = -t489 * mrSges(5,1) + t490 * mrSges(5,2);
t485 = -qJD(4) * mrSges(5,2) + t489 * mrSges(5,3);
t430 = m(5) * t448 + qJDD(4) * mrSges(5,1) - t481 * mrSges(5,3) + qJD(4) * t485 - t490 * t478 + t432;
t486 = qJD(4) * mrSges(5,1) - t490 * mrSges(5,3);
t537 = -t517 * t440 + t521 * t441;
t431 = m(5) * t449 - qJDD(4) * mrSges(5,2) + t480 * mrSges(5,3) - qJD(4) * t486 + t489 * t478 + t537;
t426 = t522 * t430 + t518 * t431;
t466 = -t515 * t482 + t542;
t531 = mrSges(4,3) * t510 + (-mrSges(4,1) * t516 + t548) * t508;
t424 = m(4) * t466 - t531 * t515 + t426;
t538 = -t518 * t430 + t522 * t431;
t425 = m(4) * t467 + t531 * t516 + t538;
t539 = -t515 * t424 + t516 * t425;
t417 = m(3) * t484 - t508 * mrSges(3,1) - t510 * mrSges(3,2) + t539;
t483 = t523 * t496 - t519 * t497;
t533 = qJDD(3) - t483;
t479 = -t510 * pkin(2) - t508 * qJ(3) + t533;
t511 = t515 ^ 2;
t465 = (-pkin(2) - t549) * t510 + (-qJ(3) + (-t511 - t512) * pkin(7)) * t508 + t533;
t447 = -t480 * pkin(4) - t488 * pkin(8) + t490 * t487 + t465;
t532 = m(6) * t447 - t452 * mrSges(6,1) + t453 * mrSges(6,2) - t473 * t468 + t474 * t469;
t527 = m(5) * t465 - t480 * mrSges(5,1) + t481 * mrSges(5,2) - t489 * t485 + t490 * t486 + t532;
t526 = -m(4) * t479 + mrSges(4,1) * t546 - t527 + (t508 * t511 + t547) * mrSges(4,3);
t436 = t526 + (mrSges(3,1) - t548) * t510 - t508 * mrSges(3,2) + m(3) * t483;
t414 = t519 * t417 + t523 * t436;
t418 = t516 * t424 + t515 * t425;
t540 = t523 * t417 - t519 * t436;
t412 = m(2) * t500 - t525 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t540;
t413 = m(2) * t501 + qJDD(1) * mrSges(2,1) - t525 * mrSges(2,2) + t414;
t541 = t524 * t412 - t520 * t413;
t536 = Ifges(4,1) * t515 + Ifges(4,4) * t516;
t535 = Ifges(4,4) * t515 + Ifges(4,2) * t516;
t530 = -t520 * t412 - t524 * t413;
t472 = Ifges(5,1) * t490 + Ifges(5,4) * t489 + Ifges(5,5) * qJD(4);
t471 = Ifges(5,4) * t490 + Ifges(5,2) * t489 + Ifges(5,6) * qJD(4);
t470 = Ifges(5,5) * t490 + Ifges(5,6) * t489 + Ifges(5,3) * qJD(4);
t456 = Ifges(6,1) * t474 + Ifges(6,4) * t473 + Ifges(6,5) * t513;
t455 = Ifges(6,4) * t474 + Ifges(6,2) * t473 + Ifges(6,6) * t513;
t454 = Ifges(6,5) * t474 + Ifges(6,6) * t473 + Ifges(6,3) * t513;
t434 = mrSges(6,2) * t447 - mrSges(6,3) * t442 + Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t509 + t473 * t454 - t513 * t455;
t433 = -mrSges(6,1) * t447 + mrSges(6,3) * t443 + Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t509 - t474 * t454 + t513 * t456;
t420 = mrSges(5,2) * t465 - mrSges(5,3) * t448 + Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * qJDD(4) - pkin(8) * t432 - qJD(4) * t471 - t517 * t433 + t521 * t434 + t489 * t470;
t419 = -mrSges(5,1) * t465 + mrSges(5,3) * t449 + Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * qJDD(4) - pkin(4) * t532 + pkin(8) * t537 + qJD(4) * t472 + t521 * t433 + t517 * t434 - t490 * t470;
t410 = mrSges(4,2) * t479 - mrSges(4,3) * t466 - pkin(7) * t426 - t518 * t419 + t522 * t420 + t536 * t510 + t516 * t545;
t409 = -mrSges(4,1) * t479 + mrSges(4,3) * t467 - pkin(3) * t527 + pkin(7) * t538 + t522 * t419 + t518 * t420 + t535 * t510 - t515 * t545;
t408 = mrSges(3,1) * g(1) - Ifges(5,3) * qJDD(4) - Ifges(6,3) * t509 + t489 * t472 - t490 * t471 - Ifges(5,6) * t480 - Ifges(5,5) * t481 + mrSges(3,3) * t484 + t473 * t456 - t474 * t455 - mrSges(4,1) * t466 + mrSges(4,2) * t467 - Ifges(6,6) * t452 - Ifges(6,5) * t453 - mrSges(5,1) * t448 + mrSges(5,2) * t449 + mrSges(6,2) * t443 - mrSges(6,1) * t442 - pkin(4) * t432 - pkin(2) * t418 - pkin(3) * t426 + (Ifges(3,6) - t534) * t510 + (-t515 * t535 + t516 * t536 + Ifges(3,5)) * t508;
t407 = -mrSges(3,2) * g(1) - mrSges(3,3) * t483 + Ifges(3,5) * t510 - t508 * Ifges(3,6) - qJ(3) * t418 - t515 * t409 + t516 * t410;
t406 = -mrSges(2,2) * g(1) - mrSges(2,3) * t501 + Ifges(2,5) * qJDD(1) - t525 * Ifges(2,6) - pkin(6) * t414 + t523 * t407 - t519 * t408;
t405 = Ifges(2,6) * qJDD(1) + t525 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t500 + t519 * t407 + t523 * t408 - pkin(1) * (-m(3) * g(1) + t418) + pkin(6) * t540;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t418; -m(1) * g(2) + t530; -m(1) * g(3) + t541; pkin(1) * t414 + t515 * t410 + t516 * t409 + pkin(2) * (-t510 * t548 + t526) + qJ(3) * t539 + mrSges(3,1) * t483 - mrSges(3,2) * t484 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t510 + mrSges(2,1) * t501 - mrSges(2,2) * t500 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t541 - t524 * t405 - t520 * t406; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t530 - t520 * t405 + t524 * t406;];
tauB = t1;
