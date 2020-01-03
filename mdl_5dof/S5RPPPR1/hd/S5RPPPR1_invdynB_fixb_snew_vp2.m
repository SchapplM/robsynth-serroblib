% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:24
% EndTime: 2020-01-03 11:20:27
% DurationCPUTime: 3.40s
% Computational Cost: add. (31387->236), mult. (72890->326), div. (0->0), fcn. (44132->10), ass. (0->114)
t514 = sin(qJ(1));
t516 = cos(qJ(1));
t493 = -t516 * g(2) - t514 * g(3);
t488 = qJDD(1) * pkin(1) + t493;
t492 = -t514 * g(2) + t516 * g(3);
t517 = qJD(1) ^ 2;
t489 = -t517 * pkin(1) + t492;
t509 = sin(pkin(7));
t512 = cos(pkin(7));
t471 = t509 * t488 + t512 * t489;
t557 = -t517 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t471;
t470 = t512 * t488 - t509 * t489;
t523 = -t517 * qJ(3) + qJDD(3) - t470;
t508 = sin(pkin(8));
t511 = cos(pkin(8));
t529 = -pkin(3) * t511 - qJ(4) * t508;
t547 = qJD(1) * t508;
t556 = (-pkin(2) + t529) * qJDD(1) + t523 - 0.2e1 * qJD(4) * t547;
t506 = -g(1) + qJDD(2);
t456 = t511 * t506 - t557 * t508;
t555 = pkin(6) * t508;
t554 = mrSges(4,2) * t508;
t553 = Ifges(4,6) * t511;
t505 = t508 ^ 2;
t552 = t505 * t517;
t507 = sin(pkin(9));
t551 = t507 * t508;
t510 = cos(pkin(9));
t550 = t508 * t510;
t457 = t508 * t506 + t557 * t511;
t486 = (-mrSges(4,1) * t511 + t554) * qJD(1);
t485 = t529 * qJD(1);
t546 = t511 * qJD(1);
t450 = t485 * t546 + t457;
t528 = -pkin(4) * t511 - pkin(6) * t550;
t548 = t556 * t510;
t443 = t528 * qJDD(1) + (-t450 + (-pkin(4) * t505 * t510 + t511 * t555) * t517) * t507 + t548;
t446 = t510 * t450 + t556 * t507;
t484 = t528 * qJD(1);
t542 = t507 ^ 2 * t552;
t544 = qJDD(1) * t507;
t444 = -pkin(4) * t542 + t484 * t546 - t544 * t555 + t446;
t513 = sin(qJ(5));
t515 = cos(qJ(5));
t441 = t515 * t443 - t513 * t444;
t525 = (-t507 * t515 - t510 * t513) * t508;
t475 = qJD(1) * t525;
t524 = (-t507 * t513 + t510 * t515) * t508;
t476 = qJD(1) * t524;
t460 = -t475 * mrSges(6,1) + t476 * mrSges(6,2);
t463 = t475 * qJD(5) + qJDD(1) * t524;
t495 = qJD(5) - t546;
t468 = -t495 * mrSges(6,2) + t475 * mrSges(6,3);
t543 = t511 * qJDD(1);
t494 = qJDD(5) - t543;
t439 = m(6) * t441 + t494 * mrSges(6,1) - t463 * mrSges(6,3) - t476 * t460 + t495 * t468;
t442 = t513 * t443 + t515 * t444;
t462 = -t476 * qJD(5) + qJDD(1) * t525;
t469 = t495 * mrSges(6,1) - t476 * mrSges(6,3);
t440 = m(6) * t442 - t494 * mrSges(6,2) + t462 * mrSges(6,3) + t475 * t460 - t495 * t469;
t431 = t515 * t439 + t513 * t440;
t445 = -t507 * t450 + t548;
t532 = mrSges(5,1) * t507 + mrSges(5,2) * t510;
t477 = t532 * t547;
t526 = mrSges(5,2) * t511 - mrSges(5,3) * t551;
t479 = t526 * qJD(1);
t527 = -mrSges(5,1) * t511 - mrSges(5,3) * t550;
t429 = m(5) * t445 + t527 * qJDD(1) + (-t477 * t550 - t479 * t511) * qJD(1) + t431;
t480 = t527 * qJD(1);
t535 = -t513 * t439 + t515 * t440;
t430 = m(5) * t446 + t526 * qJDD(1) + (-t477 * t551 + t480 * t511) * qJD(1) + t535;
t536 = -t507 * t429 + t510 * t430;
t426 = m(4) * t457 + (qJDD(1) * mrSges(4,3) + qJD(1) * t486) * t511 + t536;
t449 = t485 * t547 + qJDD(4) - t456;
t447 = -pkin(6) * t542 + (qJD(1) * t484 * t510 + pkin(4) * t544) * t508 + t449;
t522 = m(6) * t447 - t462 * mrSges(6,1) + t463 * mrSges(6,2) - t475 * t468 + t476 * t469;
t518 = -m(5) * t449 - t522;
t435 = m(4) * t456 + ((-mrSges(4,3) - t532) * qJDD(1) + (-t479 * t507 - t480 * t510 - t486) * qJD(1)) * t508 + t518;
t537 = t511 * t426 - t508 * t435;
t418 = m(3) * t471 - t517 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t537;
t427 = t510 * t429 + t507 * t430;
t466 = -qJDD(1) * pkin(2) + t523;
t519 = -m(4) * t466 + mrSges(4,1) * t543 - t427 + (t511 ^ 2 * t517 + t552) * mrSges(4,3);
t423 = m(3) * t470 - t517 * mrSges(3,2) + (mrSges(3,1) - t554) * qJDD(1) + t519;
t538 = t512 * t418 - t509 * t423;
t412 = m(2) * t492 - t517 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t538;
t414 = t509 * t418 + t512 * t423;
t413 = m(2) * t493 + qJDD(1) * mrSges(2,1) - t517 * mrSges(2,2) + t414;
t549 = t514 * t412 + t516 * t413;
t419 = t508 * t426 + t511 * t435;
t541 = m(3) * t506 + t419;
t539 = -t516 * t412 + t514 * t413;
t531 = Ifges(4,1) * t508 + Ifges(4,4) * t511;
t530 = Ifges(5,5) * t510 - Ifges(5,6) * t507;
t521 = -Ifges(5,5) * t511 + (Ifges(5,1) * t510 - Ifges(5,4) * t507) * t508;
t520 = -Ifges(5,6) * t511 + (Ifges(5,4) * t510 - Ifges(5,2) * t507) * t508;
t487 = (Ifges(4,5) * t508 + t553) * qJD(1);
t474 = t521 * qJD(1);
t473 = t520 * qJD(1);
t472 = (-Ifges(5,3) * t511 + t530 * t508) * qJD(1);
t453 = Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t495;
t452 = Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t495;
t451 = Ifges(6,5) * t476 + Ifges(6,6) * t475 + Ifges(6,3) * t495;
t433 = mrSges(6,2) * t447 - mrSges(6,3) * t441 + Ifges(6,1) * t463 + Ifges(6,4) * t462 + Ifges(6,5) * t494 + t475 * t451 - t495 * t452;
t432 = -mrSges(6,1) * t447 + mrSges(6,3) * t442 + Ifges(6,4) * t463 + Ifges(6,2) * t462 + Ifges(6,6) * t494 - t476 * t451 + t495 * t453;
t421 = mrSges(5,2) * t449 - mrSges(5,3) * t445 - pkin(6) * t431 - t513 * t432 + t515 * t433 + (-t472 * t551 + t473 * t511) * qJD(1) + t521 * qJDD(1);
t420 = -mrSges(5,1) * t449 + mrSges(5,3) * t446 + t513 * t433 + t515 * t432 - pkin(4) * t522 + pkin(6) * t535 + (-t472 * t550 - t511 * t474) * qJD(1) + t520 * qJDD(1);
t415 = -mrSges(4,1) * t466 - mrSges(5,1) * t445 - mrSges(6,1) * t441 + mrSges(5,2) * t446 + mrSges(6,2) * t442 + mrSges(4,3) * t457 - Ifges(6,5) * t463 - Ifges(6,6) * t462 - Ifges(6,3) * t494 - pkin(3) * t427 - pkin(4) * t431 - t476 * t452 + t475 * t453 + (Ifges(4,2) + Ifges(5,3)) * t543 + ((Ifges(4,4) - t530) * qJDD(1) + (-t473 * t510 - t474 * t507 - t487) * qJD(1)) * t508;
t408 = mrSges(4,2) * t466 - mrSges(4,3) * t456 - qJ(4) * t427 + t531 * qJDD(1) - t507 * t420 + t510 * t421 + t487 * t546;
t407 = t517 * Ifges(3,5) - mrSges(3,1) * t506 + mrSges(3,3) * t471 - mrSges(4,1) * t456 + mrSges(4,2) * t457 - t507 * t421 - t510 * t420 - pkin(3) * t518 - qJ(4) * t536 - pkin(2) * t419 + (-t553 + Ifges(3,6) + (pkin(3) * t532 - Ifges(4,5)) * t508) * qJDD(1) + (-pkin(3) * (-t479 * t551 - t480 * t550) + (-t508 * (Ifges(4,4) * t508 + Ifges(4,2) * t511) + t511 * t531) * qJD(1)) * qJD(1);
t406 = mrSges(3,2) * t506 - mrSges(3,3) * t470 + Ifges(3,5) * qJDD(1) - t517 * Ifges(3,6) - qJ(3) * t419 + t511 * t408 - t508 * t415;
t405 = -mrSges(2,2) * g(1) - mrSges(2,3) * t493 + Ifges(2,5) * qJDD(1) - t517 * Ifges(2,6) - qJ(2) * t414 + t512 * t406 - t509 * t407;
t404 = mrSges(2,1) * g(1) + mrSges(2,3) * t492 + t517 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t541 + qJ(2) * t538 + t509 * t406 + t512 * t407;
t1 = [(-m(1) - m(2)) * g(1) + t541; -m(1) * g(2) + t549; -m(1) * g(3) + t539; pkin(1) * t414 + t511 * t415 + pkin(2) * t519 + qJ(3) * t537 + t508 * t408 + mrSges(3,1) * t470 - mrSges(3,2) * t471 + mrSges(2,1) * t493 - mrSges(2,2) * t492 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(2) * t554 + Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t539 + t516 * t404 + t514 * t405; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t549 + t514 * t404 - t516 * t405;];
tauB = t1;
