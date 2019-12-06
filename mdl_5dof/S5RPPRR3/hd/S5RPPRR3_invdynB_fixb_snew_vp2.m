% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:55
% EndTime: 2019-12-05 17:41:59
% DurationCPUTime: 4.10s
% Computational Cost: add. (48242->248), mult. (106816->311), div. (0->0), fcn. (71884->10), ass. (0->108)
t521 = qJD(1) ^ 2;
t513 = cos(pkin(9));
t546 = pkin(3) * t513;
t511 = sin(pkin(9));
t545 = mrSges(4,2) * t511;
t508 = t513 ^ 2;
t544 = t508 * t521;
t517 = sin(qJ(1));
t520 = cos(qJ(1));
t494 = t520 * g(2) + t517 * g(3);
t491 = qJDD(1) * pkin(1) + t494;
t493 = t517 * g(2) - t520 * g(3);
t492 = -t521 * pkin(1) + t493;
t512 = sin(pkin(8));
t514 = cos(pkin(8));
t479 = t512 * t491 + t514 * t492;
t472 = -t521 * pkin(2) + qJDD(1) * qJ(3) + t479;
t510 = -g(1) + qJDD(2);
t540 = qJD(1) * qJD(3);
t543 = t513 * t510 - 0.2e1 * t511 * t540;
t458 = (-pkin(6) * qJDD(1) + t521 * t546 - t472) * t511 + t543;
t462 = t511 * t510 + (t472 + 0.2e1 * t540) * t513;
t539 = qJDD(1) * t513;
t459 = -pkin(3) * t544 + pkin(6) * t539 + t462;
t516 = sin(qJ(4));
t519 = cos(qJ(4));
t443 = t519 * t458 - t516 * t459;
t526 = t511 * t519 + t513 * t516;
t525 = -t511 * t516 + t513 * t519;
t484 = t525 * qJD(1);
t541 = t484 * qJD(4);
t477 = t526 * qJDD(1) + t541;
t485 = t526 * qJD(1);
t439 = (-t477 + t541) * pkin(7) + (t484 * t485 + qJDD(4)) * pkin(4) + t443;
t444 = t516 * t458 + t519 * t459;
t476 = -t485 * qJD(4) + t525 * qJDD(1);
t482 = qJD(4) * pkin(4) - t485 * pkin(7);
t483 = t484 ^ 2;
t440 = -t483 * pkin(4) + t476 * pkin(7) - qJD(4) * t482 + t444;
t515 = sin(qJ(5));
t518 = cos(qJ(5));
t437 = t518 * t439 - t515 * t440;
t470 = t518 * t484 - t515 * t485;
t448 = t470 * qJD(5) + t515 * t476 + t518 * t477;
t471 = t515 * t484 + t518 * t485;
t454 = -t470 * mrSges(6,1) + t471 * mrSges(6,2);
t509 = qJD(4) + qJD(5);
t463 = -t509 * mrSges(6,2) + t470 * mrSges(6,3);
t506 = qJDD(4) + qJDD(5);
t435 = m(6) * t437 + t506 * mrSges(6,1) - t448 * mrSges(6,3) - t471 * t454 + t509 * t463;
t438 = t515 * t439 + t518 * t440;
t447 = -t471 * qJD(5) + t518 * t476 - t515 * t477;
t464 = t509 * mrSges(6,1) - t471 * mrSges(6,3);
t436 = m(6) * t438 - t506 * mrSges(6,2) + t447 * mrSges(6,3) + t470 * t454 - t509 * t464;
t427 = t518 * t435 + t515 * t436;
t474 = -t484 * mrSges(5,1) + t485 * mrSges(5,2);
t480 = -qJD(4) * mrSges(5,2) + t484 * mrSges(5,3);
t425 = m(5) * t443 + qJDD(4) * mrSges(5,1) - t477 * mrSges(5,3) + qJD(4) * t480 - t485 * t474 + t427;
t481 = qJD(4) * mrSges(5,1) - t485 * mrSges(5,3);
t533 = -t515 * t435 + t518 * t436;
t426 = m(5) * t444 - qJDD(4) * mrSges(5,2) + t476 * mrSges(5,3) - qJD(4) * t481 + t484 * t474 + t533;
t421 = t519 * t425 + t516 * t426;
t461 = -t511 * t472 + t543;
t524 = mrSges(4,3) * qJDD(1) + t521 * (-mrSges(4,1) * t513 + t545);
t419 = m(4) * t461 - t524 * t511 + t421;
t534 = -t516 * t425 + t519 * t426;
t420 = m(4) * t462 + t524 * t513 + t534;
t535 = -t511 * t419 + t513 * t420;
t412 = m(3) * t479 - t521 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t535;
t478 = t514 * t491 - t512 * t492;
t529 = qJDD(3) - t478;
t466 = -qJDD(1) * pkin(2) - t521 * qJ(3) + t529;
t507 = t511 ^ 2;
t460 = (-pkin(2) - t546) * qJDD(1) + (-qJ(3) + (-t507 - t508) * pkin(6)) * t521 + t529;
t442 = -t476 * pkin(4) - t483 * pkin(7) + t485 * t482 + t460;
t528 = m(6) * t442 - t447 * mrSges(6,1) + t448 * mrSges(6,2) - t470 * t463 + t471 * t464;
t523 = m(5) * t460 - t476 * mrSges(5,1) + t477 * mrSges(5,2) - t484 * t480 + t485 * t481 + t528;
t522 = -m(4) * t466 + mrSges(4,1) * t539 - t523 + (t507 * t521 + t544) * mrSges(4,3);
t431 = t522 + (mrSges(3,1) - t545) * qJDD(1) - t521 * mrSges(3,2) + m(3) * t478;
t409 = t512 * t412 + t514 * t431;
t413 = t513 * t419 + t511 * t420;
t530 = Ifges(4,5) * t511 + Ifges(4,6) * t513;
t542 = t521 * t530;
t538 = m(3) * t510 + t413;
t536 = t514 * t412 - t512 * t431;
t407 = m(2) * t493 - t521 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t536;
t408 = m(2) * t494 + qJDD(1) * mrSges(2,1) - t521 * mrSges(2,2) + t409;
t537 = t520 * t407 - t517 * t408;
t532 = Ifges(4,1) * t511 + Ifges(4,4) * t513;
t531 = Ifges(4,4) * t511 + Ifges(4,2) * t513;
t527 = -t517 * t407 - t520 * t408;
t469 = Ifges(5,1) * t485 + Ifges(5,4) * t484 + Ifges(5,5) * qJD(4);
t468 = Ifges(5,4) * t485 + Ifges(5,2) * t484 + Ifges(5,6) * qJD(4);
t467 = Ifges(5,5) * t485 + Ifges(5,6) * t484 + Ifges(5,3) * qJD(4);
t451 = Ifges(6,1) * t471 + Ifges(6,4) * t470 + Ifges(6,5) * t509;
t450 = Ifges(6,4) * t471 + Ifges(6,2) * t470 + Ifges(6,6) * t509;
t449 = Ifges(6,5) * t471 + Ifges(6,6) * t470 + Ifges(6,3) * t509;
t429 = mrSges(6,2) * t442 - mrSges(6,3) * t437 + Ifges(6,1) * t448 + Ifges(6,4) * t447 + Ifges(6,5) * t506 + t470 * t449 - t509 * t450;
t428 = -mrSges(6,1) * t442 + mrSges(6,3) * t438 + Ifges(6,4) * t448 + Ifges(6,2) * t447 + Ifges(6,6) * t506 - t471 * t449 + t509 * t451;
t415 = mrSges(5,2) * t460 - mrSges(5,3) * t443 + Ifges(5,1) * t477 + Ifges(5,4) * t476 + Ifges(5,5) * qJDD(4) - pkin(7) * t427 - qJD(4) * t468 - t515 * t428 + t518 * t429 + t484 * t467;
t414 = -mrSges(5,1) * t460 + mrSges(5,3) * t444 + Ifges(5,4) * t477 + Ifges(5,2) * t476 + Ifges(5,6) * qJDD(4) - pkin(4) * t528 + pkin(7) * t533 + qJD(4) * t469 + t518 * t428 + t515 * t429 - t485 * t467;
t405 = mrSges(4,2) * t466 - mrSges(4,3) * t461 - pkin(6) * t421 + t532 * qJDD(1) - t516 * t414 + t519 * t415 + t513 * t542;
t404 = -mrSges(4,1) * t466 + mrSges(4,3) * t462 - pkin(3) * t523 + pkin(6) * t534 + t531 * qJDD(1) + t519 * t414 + t516 * t415 - t511 * t542;
t403 = -Ifges(5,3) * qJDD(4) - mrSges(3,1) * t510 - Ifges(6,3) * t506 - t485 * t468 - Ifges(5,6) * t476 - Ifges(5,5) * t477 + mrSges(3,3) * t479 + t484 * t469 + t470 * t451 - t471 * t450 - mrSges(4,1) * t461 + mrSges(4,2) * t462 - Ifges(6,6) * t447 - Ifges(6,5) * t448 - mrSges(5,1) * t443 + mrSges(5,2) * t444 - mrSges(6,1) * t437 + mrSges(6,2) * t438 - pkin(2) * t413 - pkin(3) * t421 - pkin(4) * t427 + (Ifges(3,6) - t530) * qJDD(1) + (-t511 * t531 + t513 * t532 + Ifges(3,5)) * t521;
t402 = mrSges(3,2) * t510 - mrSges(3,3) * t478 + Ifges(3,5) * qJDD(1) - t521 * Ifges(3,6) - qJ(3) * t413 - t511 * t404 + t513 * t405;
t401 = -mrSges(2,2) * g(1) - mrSges(2,3) * t494 + Ifges(2,5) * qJDD(1) - t521 * Ifges(2,6) - qJ(2) * t409 + t514 * t402 - t512 * t403;
t400 = mrSges(2,1) * g(1) + mrSges(2,3) * t493 + t521 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t538 + qJ(2) * t536 + t512 * t402 + t514 * t403;
t1 = [(-m(1) - m(2)) * g(1) + t538; -m(1) * g(2) + t527; -m(1) * g(3) + t537; pkin(1) * t409 + t511 * t405 + t513 * t404 + pkin(2) * t522 + qJ(3) * t535 + mrSges(3,1) * t478 - mrSges(3,2) * t479 + mrSges(2,1) * t494 - mrSges(2,2) * t493 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(2) * t545 + Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t537 - t520 * t400 - t517 * t401; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t527 - t517 * t400 + t520 * t401;];
tauB = t1;
