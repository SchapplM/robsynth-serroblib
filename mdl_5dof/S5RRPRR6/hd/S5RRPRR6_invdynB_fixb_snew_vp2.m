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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:35:56
% EndTime: 2019-12-05 18:36:01
% DurationCPUTime: 3.92s
% Computational Cost: add. (54801->249), mult. (75094->329), div. (0->0), fcn. (46096->10), ass. (0->111)
t512 = sin(qJ(1));
t516 = cos(qJ(1));
t494 = t516 * g(2) + t512 * g(3);
t485 = qJDD(1) * pkin(1) + t494;
t493 = t512 * g(2) - g(3) * t516;
t517 = qJD(1) ^ 2;
t486 = -pkin(1) * t517 + t493;
t511 = sin(qJ(2));
t515 = cos(qJ(2));
t468 = t511 * t485 + t515 * t486;
t505 = (qJD(1) + qJD(2));
t502 = t505 ^ 2;
t503 = qJDD(1) + qJDD(2);
t546 = -pkin(2) * t502 + qJ(3) * t503 + (2 * qJD(3) * t505) + t468;
t507 = sin(pkin(9));
t508 = cos(pkin(9));
t455 = -t508 * g(1) - t546 * t507;
t545 = mrSges(4,2) * t507;
t544 = mrSges(4,3) * t503;
t543 = t502 * t507 ^ 2;
t542 = t503 * t508;
t541 = t507 * t505;
t510 = sin(qJ(4));
t540 = t507 * t510;
t514 = cos(qJ(4));
t539 = t507 * t514;
t538 = t508 * t505;
t456 = -g(1) * t507 + t546 * t508;
t479 = (-mrSges(4,1) * t508 + t545) * t505;
t525 = -pkin(3) * t508 - pkin(7) * t507;
t481 = t525 * t505;
t444 = t481 * t538 + t456;
t467 = t485 * t515 - t511 * t486;
t521 = -qJ(3) * t502 + qJDD(3) - t467;
t454 = (-pkin(2) + t525) * t503 + t521;
t453 = t514 * t454;
t535 = qJD(4) * t505;
t478 = (t503 * t514 - t510 * t535) * t507;
t489 = qJDD(4) - t542;
t490 = qJD(4) - t538;
t437 = pkin(4) * t489 - pkin(8) * t478 + t453 + (-pkin(4) * t514 * t543 - pkin(8) * t490 * t541 - t444) * t510;
t440 = t514 * t444 + t510 * t454;
t532 = t505 * t539;
t476 = pkin(4) * t490 - pkin(8) * t532;
t477 = (-t503 * t510 - t514 * t535) * t507;
t534 = t510 ^ 2 * t543;
t438 = -pkin(4) * t534 + pkin(8) * t477 - t476 * t490 + t440;
t509 = sin(qJ(5));
t513 = cos(qJ(5));
t435 = t437 * t513 - t438 * t509;
t469 = (-t509 * t514 - t510 * t513) * t541;
t447 = qJD(5) * t469 + t477 * t509 + t478 * t513;
t470 = (-t509 * t510 + t513 * t514) * t541;
t457 = -mrSges(6,1) * t469 + mrSges(6,2) * t470;
t488 = qJD(5) + t490;
t462 = -mrSges(6,2) * t488 + mrSges(6,3) * t469;
t487 = qJDD(5) + t489;
t433 = m(6) * t435 + mrSges(6,1) * t487 - mrSges(6,3) * t447 - t457 * t470 + t462 * t488;
t436 = t437 * t509 + t438 * t513;
t446 = -qJD(5) * t470 + t477 * t513 - t478 * t509;
t463 = mrSges(6,1) * t488 - mrSges(6,3) * t470;
t434 = m(6) * t436 - mrSges(6,2) * t487 + mrSges(6,3) * t446 + t457 * t469 - t463 * t488;
t425 = t513 * t433 + t509 * t434;
t439 = -t444 * t510 + t453;
t533 = t505 * t540;
t473 = -mrSges(5,2) * t490 - mrSges(5,3) * t533;
t475 = (mrSges(5,1) * t510 + mrSges(5,2) * t514) * t541;
t423 = m(5) * t439 + mrSges(5,1) * t489 - mrSges(5,3) * t478 + t473 * t490 - t475 * t532 + t425;
t474 = mrSges(5,1) * t490 - mrSges(5,3) * t532;
t526 = -t433 * t509 + t513 * t434;
t424 = m(5) * t440 - mrSges(5,2) * t489 + mrSges(5,3) * t477 - t474 * t490 - t475 * t533 + t526;
t527 = -t423 * t510 + t514 * t424;
t420 = m(4) * t456 + (t479 * t505 + t544) * t508 + t527;
t443 = t481 * t541 - t455;
t441 = -pkin(4) * t477 - pkin(8) * t534 + t476 * t532 + t443;
t520 = m(6) * t441 - t446 * mrSges(6,1) + mrSges(6,2) * t447 - t469 * t462 + t463 * t470;
t518 = -m(5) * t443 + t477 * mrSges(5,1) - mrSges(5,2) * t478 - t520;
t429 = m(4) * t455 + (-t544 + (-t473 * t510 - t474 * t514 - t479) * t505) * t507 + t518;
t528 = t508 * t420 - t429 * t507;
t412 = m(3) * t468 - mrSges(3,1) * t502 - mrSges(3,2) * t503 + t528;
t421 = t423 * t514 + t424 * t510;
t460 = -pkin(2) * t503 + t521;
t519 = -m(4) * t460 + mrSges(4,1) * t542 - t421 + (t502 * t508 ^ 2 + t543) * mrSges(4,3);
t417 = m(3) * t467 - mrSges(3,2) * t502 + (mrSges(3,1) - t545) * t503 + t519;
t408 = t511 * t412 + t515 * t417;
t414 = t507 * t420 + t508 * t429;
t529 = t515 * t412 - t511 * t417;
t406 = m(2) * t493 - mrSges(2,1) * t517 - qJDD(1) * mrSges(2,2) + t529;
t407 = m(2) * t494 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t517 + t408;
t530 = t516 * t406 - t407 * t512;
t524 = Ifges(4,1) * t507 + Ifges(4,4) * t508;
t523 = Ifges(4,5) * t507 + Ifges(4,6) * t508;
t522 = -t406 * t512 - t407 * t516;
t480 = t523 * t505;
t466 = Ifges(5,5) * t490 + (Ifges(5,1) * t514 - Ifges(5,4) * t510) * t541;
t465 = Ifges(5,6) * t490 + (Ifges(5,4) * t514 - Ifges(5,2) * t510) * t541;
t464 = Ifges(5,3) * t490 + (Ifges(5,5) * t514 - Ifges(5,6) * t510) * t541;
t450 = Ifges(6,1) * t470 + Ifges(6,4) * t469 + Ifges(6,5) * t488;
t449 = Ifges(6,4) * t470 + Ifges(6,2) * t469 + Ifges(6,6) * t488;
t448 = Ifges(6,5) * t470 + Ifges(6,6) * t469 + Ifges(6,3) * t488;
t427 = mrSges(6,2) * t441 - mrSges(6,3) * t435 + Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t487 + t448 * t469 - t449 * t488;
t426 = -mrSges(6,1) * t441 + mrSges(6,3) * t436 + Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t487 - t448 * t470 + t450 * t488;
t415 = mrSges(5,2) * t443 - mrSges(5,3) * t439 + Ifges(5,1) * t478 + Ifges(5,4) * t477 + Ifges(5,5) * t489 - pkin(8) * t425 - t426 * t509 + t427 * t513 - t464 * t533 - t465 * t490;
t413 = -mrSges(5,1) * t443 + mrSges(5,3) * t440 + Ifges(5,4) * t478 + Ifges(5,2) * t477 + Ifges(5,6) * t489 - pkin(4) * t520 + pkin(8) * t526 + t513 * t426 + t509 * t427 - t464 * t532 + t490 * t466;
t409 = Ifges(4,2) * t542 - mrSges(4,1) * t460 + mrSges(4,3) * t456 - Ifges(5,5) * t478 - Ifges(5,6) * t477 - Ifges(5,3) * t489 - mrSges(5,1) * t439 + mrSges(5,2) * t440 - Ifges(6,5) * t447 - Ifges(6,6) * t446 - Ifges(6,3) * t487 - t470 * t449 + t469 * t450 - mrSges(6,1) * t435 + mrSges(6,2) * t436 - pkin(4) * t425 - pkin(3) * t421 + (Ifges(4,4) * t503 + (-t465 * t514 - t466 * t510 - t480) * t505) * t507;
t404 = mrSges(4,2) * t460 - mrSges(4,3) * t455 - pkin(7) * t421 - t413 * t510 + t415 * t514 + t480 * t538 + t503 * t524;
t403 = t502 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t468 - mrSges(4,1) * t455 + mrSges(4,2) * t456 - t510 * t415 - t514 * t413 - pkin(3) * t518 - pkin(7) * t527 - pkin(2) * t414 + (Ifges(3,6) - t523) * t503 + (-pkin(3) * (-t473 * t540 - t474 * t539) + (-t507 * (Ifges(4,4) * t507 + Ifges(4,2) * t508) + t508 * t524) * t505) * t505;
t402 = -mrSges(3,2) * g(1) - mrSges(3,3) * t467 + Ifges(3,5) * t503 - Ifges(3,6) * t502 - qJ(3) * t414 + t404 * t508 - t409 * t507;
t401 = -mrSges(2,2) * g(1) - mrSges(2,3) * t494 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t517 - pkin(6) * t408 + t402 * t515 - t403 * t511;
t400 = Ifges(2,6) * qJDD(1) + t517 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t493 + t511 * t402 + t515 * t403 - pkin(1) * (-m(3) * g(1) + t414) + pkin(6) * t529;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t414; -m(1) * g(2) + t522; -m(1) * g(3) + t530; pkin(1) * t408 + t508 * t409 + pkin(2) * (-t503 * t545 + t519) + qJ(3) * t528 + t507 * t404 + mrSges(3,1) * t467 - mrSges(3,2) * t468 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t503 + mrSges(2,1) * t494 - mrSges(2,2) * t493 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t530 - t516 * t400 - t512 * t401; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t522 - t512 * t400 + t516 * t401;];
tauB = t1;
