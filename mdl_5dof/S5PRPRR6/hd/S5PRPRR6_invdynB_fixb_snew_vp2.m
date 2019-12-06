% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:50
% EndTime: 2019-12-05 15:56:58
% DurationCPUTime: 5.70s
% Computational Cost: add. (63533->246), mult. (138025->319), div. (0->0), fcn. (100282->12), ass. (0->116)
t525 = qJD(2) ^ 2;
t512 = sin(pkin(10));
t515 = cos(pkin(10));
t519 = sin(qJ(4));
t522 = cos(qJ(4));
t532 = t512 * t519 - t515 * t522;
t494 = t532 * qJD(2);
t513 = sin(pkin(9));
t516 = cos(pkin(9));
t501 = t513 * g(1) - t516 * g(2);
t502 = -t516 * g(1) - t513 * g(2);
t511 = -g(3) + qJDD(1);
t514 = sin(pkin(5));
t517 = cos(pkin(5));
t520 = sin(qJ(2));
t523 = cos(qJ(2));
t477 = -t520 * t502 + (t501 * t517 + t511 * t514) * t523;
t533 = t512 * t522 + t515 * t519;
t495 = t533 * qJD(2);
t544 = t495 * qJD(4);
t484 = -t532 * qJDD(2) - t544;
t554 = pkin(3) * t515;
t553 = mrSges(4,2) * t512;
t528 = qJDD(3) - t477;
t470 = -qJDD(2) * pkin(2) - t525 * qJ(3) + t528;
t509 = t512 ^ 2;
t549 = t517 * t520;
t550 = t514 * t520;
t478 = t501 * t549 + t523 * t502 + t511 * t550;
t473 = -t525 * pkin(2) + qJDD(2) * qJ(3) + t478;
t492 = -t514 * t501 + t517 * t511;
t543 = qJD(2) * qJD(3);
t547 = t515 * t492 - 0.2e1 * t512 * t543;
t456 = (-pkin(7) * qJDD(2) + t525 * t554 - t473) * t512 + t547;
t459 = t512 * t492 + (t473 + 0.2e1 * t543) * t515;
t542 = qJDD(2) * t515;
t510 = t515 ^ 2;
t551 = t510 * t525;
t457 = -pkin(3) * t551 + pkin(7) * t542 + t459;
t452 = t519 * t456 + t522 * t457;
t483 = t494 * pkin(4) - t495 * pkin(8);
t524 = qJD(4) ^ 2;
t450 = -t524 * pkin(4) + qJDD(4) * pkin(8) - t494 * t483 + t452;
t467 = (-pkin(2) - t554) * qJDD(2) + (-qJ(3) + (-t509 - t510) * pkin(7)) * t525 + t528;
t545 = t494 * qJD(4);
t485 = t533 * qJDD(2) - t545;
t453 = (-t485 + t545) * pkin(8) + (-t484 + t544) * pkin(4) + t467;
t518 = sin(qJ(5));
t521 = cos(qJ(5));
t447 = -t518 * t450 + t521 * t453;
t486 = t521 * qJD(4) - t518 * t495;
t466 = t486 * qJD(5) + t518 * qJDD(4) + t521 * t485;
t487 = t518 * qJD(4) + t521 * t495;
t468 = -t486 * mrSges(6,1) + t487 * mrSges(6,2);
t493 = qJD(5) + t494;
t471 = -t493 * mrSges(6,2) + t486 * mrSges(6,3);
t482 = qJDD(5) - t484;
t445 = m(6) * t447 + t482 * mrSges(6,1) - t466 * mrSges(6,3) - t487 * t468 + t493 * t471;
t448 = t521 * t450 + t518 * t453;
t465 = -t487 * qJD(5) + t521 * qJDD(4) - t518 * t485;
t472 = t493 * mrSges(6,1) - t487 * mrSges(6,3);
t446 = m(6) * t448 - t482 * mrSges(6,2) + t465 * mrSges(6,3) + t486 * t468 - t493 * t472;
t437 = t521 * t445 + t518 * t446;
t490 = -qJD(4) * mrSges(5,2) - t494 * mrSges(5,3);
t491 = qJD(4) * mrSges(5,1) - t495 * mrSges(5,3);
t527 = m(5) * t467 - t484 * mrSges(5,1) + t485 * mrSges(5,2) + t494 * t490 + t495 * t491 + t437;
t526 = -m(4) * t470 + mrSges(4,1) * t542 - t527 + (t509 * t525 + t551) * mrSges(4,3);
t433 = (mrSges(3,1) - t553) * qJDD(2) + t526 - t525 * mrSges(3,2) + m(3) * t477;
t552 = t433 * t523;
t480 = t494 * mrSges(5,1) + t495 * mrSges(5,2);
t538 = -t518 * t445 + t521 * t446;
t436 = m(5) * t452 - qJDD(4) * mrSges(5,2) + t484 * mrSges(5,3) - qJD(4) * t491 - t494 * t480 + t538;
t451 = t522 * t456 - t519 * t457;
t449 = -qJDD(4) * pkin(4) - t524 * pkin(8) + t495 * t483 - t451;
t529 = -m(6) * t449 + t465 * mrSges(6,1) - t466 * mrSges(6,2) + t486 * t471 - t487 * t472;
t441 = m(5) * t451 + qJDD(4) * mrSges(5,1) - t485 * mrSges(5,3) + qJD(4) * t490 - t495 * t480 + t529;
t430 = t519 * t436 + t522 * t441;
t458 = -t512 * t473 + t547;
t531 = mrSges(4,3) * qJDD(2) + t525 * (-mrSges(4,1) * t515 + t553);
t428 = m(4) * t458 - t531 * t512 + t430;
t539 = t522 * t436 - t519 * t441;
t429 = m(4) * t459 + t531 * t515 + t539;
t540 = -t512 * t428 + t515 * t429;
t419 = m(3) * t478 - t525 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t540;
t422 = t515 * t428 + t512 * t429;
t421 = m(3) * t492 + t422;
t409 = t419 * t549 - t514 * t421 + t517 * t552;
t407 = m(2) * t501 + t409;
t415 = t523 * t419 - t520 * t433;
t414 = m(2) * t502 + t415;
t548 = t516 * t407 + t513 * t414;
t535 = Ifges(4,5) * t512 + Ifges(4,6) * t515;
t546 = t525 * t535;
t408 = t419 * t550 + t517 * t421 + t514 * t552;
t541 = -t513 * t407 + t516 * t414;
t537 = Ifges(4,1) * t512 + Ifges(4,4) * t515;
t536 = Ifges(4,4) * t512 + Ifges(4,2) * t515;
t460 = Ifges(6,5) * t487 + Ifges(6,6) * t486 + Ifges(6,3) * t493;
t462 = Ifges(6,1) * t487 + Ifges(6,4) * t486 + Ifges(6,5) * t493;
t438 = -mrSges(6,1) * t449 + mrSges(6,3) * t448 + Ifges(6,4) * t466 + Ifges(6,2) * t465 + Ifges(6,6) * t482 - t487 * t460 + t493 * t462;
t461 = Ifges(6,4) * t487 + Ifges(6,2) * t486 + Ifges(6,6) * t493;
t439 = mrSges(6,2) * t449 - mrSges(6,3) * t447 + Ifges(6,1) * t466 + Ifges(6,4) * t465 + Ifges(6,5) * t482 + t486 * t460 - t493 * t461;
t474 = Ifges(5,5) * t495 - Ifges(5,6) * t494 + Ifges(5,3) * qJD(4);
t475 = Ifges(5,4) * t495 - Ifges(5,2) * t494 + Ifges(5,6) * qJD(4);
t423 = mrSges(5,2) * t467 - mrSges(5,3) * t451 + Ifges(5,1) * t485 + Ifges(5,4) * t484 + Ifges(5,5) * qJDD(4) - pkin(8) * t437 - qJD(4) * t475 - t518 * t438 + t521 * t439 - t494 * t474;
t476 = Ifges(5,1) * t495 - Ifges(5,4) * t494 + Ifges(5,5) * qJD(4);
t424 = -mrSges(5,1) * t467 - mrSges(6,1) * t447 + mrSges(6,2) * t448 + mrSges(5,3) * t452 + Ifges(5,4) * t485 - Ifges(6,5) * t466 + Ifges(5,2) * t484 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t465 - Ifges(6,3) * t482 - pkin(4) * t437 + qJD(4) * t476 - t487 * t461 + t486 * t462 - t495 * t474;
t410 = -mrSges(4,1) * t470 + mrSges(4,3) * t459 - pkin(3) * t527 + pkin(7) * t539 + t536 * qJDD(2) + t519 * t423 + t522 * t424 - t512 * t546;
t411 = mrSges(4,2) * t470 - mrSges(4,3) * t458 - pkin(7) * t430 + t537 * qJDD(2) + t522 * t423 - t519 * t424 + t515 * t546;
t404 = mrSges(3,2) * t492 - mrSges(3,3) * t477 + Ifges(3,5) * qJDD(2) - t525 * Ifges(3,6) - qJ(3) * t422 - t512 * t410 + t515 * t411;
t405 = -pkin(2) * t422 + mrSges(3,3) * t478 - mrSges(3,1) * t492 - pkin(3) * t430 - mrSges(4,1) * t458 + mrSges(4,2) * t459 - t518 * t439 - t521 * t438 - pkin(4) * t529 - pkin(8) * t538 - Ifges(5,5) * t485 - Ifges(5,6) * t484 - Ifges(5,3) * qJDD(4) - t495 * t475 - t494 * t476 - mrSges(5,1) * t451 + mrSges(5,2) * t452 + (Ifges(3,6) - t535) * qJDD(2) + (-t512 * t536 + t515 * t537 + Ifges(3,5)) * t525;
t530 = pkin(6) * t415 + t404 * t520 + t405 * t523;
t403 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t477 - mrSges(3,2) * t478 + t512 * t411 + t515 * t410 + pkin(2) * (-qJDD(2) * t553 + t526) + qJ(3) * t540;
t402 = mrSges(2,2) * t511 - mrSges(2,3) * t501 + t523 * t404 - t520 * t405 + (-t408 * t514 - t409 * t517) * pkin(6);
t401 = -mrSges(2,1) * t511 + mrSges(2,3) * t502 - pkin(1) * t408 - t514 * t403 + t530 * t517;
t1 = [-m(1) * g(1) + t541; -m(1) * g(2) + t548; -m(1) * g(3) + m(2) * t511 + t408; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t548 - t513 * t401 + t516 * t402; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t541 + t516 * t401 + t513 * t402; -mrSges(1,1) * g(2) + mrSges(2,1) * t501 + mrSges(1,2) * g(1) - mrSges(2,2) * t502 + pkin(1) * t409 + t517 * t403 + t530 * t514;];
tauB = t1;
