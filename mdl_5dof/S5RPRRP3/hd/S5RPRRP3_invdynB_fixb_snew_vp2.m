% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:26
% EndTime: 2019-12-05 18:03:30
% DurationCPUTime: 2.43s
% Computational Cost: add. (23998->248), mult. (47662->301), div. (0->0), fcn. (28134->8), ass. (0->98)
t547 = Ifges(5,1) + Ifges(6,1);
t544 = Ifges(5,4) + Ifges(6,4);
t543 = Ifges(5,5) + Ifges(6,5);
t546 = Ifges(5,2) + Ifges(6,2);
t542 = -Ifges(5,6) - Ifges(6,6);
t545 = -Ifges(5,3) - Ifges(6,3);
t518 = sin(qJ(1));
t521 = cos(qJ(1));
t502 = t521 * g(2) + t518 * g(3);
t493 = qJDD(1) * pkin(1) + t502;
t501 = t518 * g(2) - t521 * g(3);
t522 = qJD(1) ^ 2;
t495 = -t522 * pkin(1) + t501;
t514 = sin(pkin(8));
t515 = cos(pkin(8));
t474 = t514 * t493 + t515 * t495;
t470 = -t522 * pkin(2) + qJDD(1) * pkin(6) + t474;
t513 = -g(1) + qJDD(2);
t517 = sin(qJ(3));
t520 = cos(qJ(3));
t455 = -t517 * t470 + t520 * t513;
t536 = qJD(1) * qJD(3);
t532 = t520 * t536;
t496 = t517 * qJDD(1) + t532;
t447 = (-t496 + t532) * pkin(7) + (t517 * t520 * t522 + qJDD(3)) * pkin(3) + t455;
t456 = t520 * t470 + t517 * t513;
t497 = t520 * qJDD(1) - t517 * t536;
t538 = qJD(1) * t517;
t500 = qJD(3) * pkin(3) - pkin(7) * t538;
t512 = t520 ^ 2;
t448 = -t512 * t522 * pkin(3) + t497 * pkin(7) - qJD(3) * t500 + t456;
t516 = sin(qJ(4));
t519 = cos(qJ(4));
t442 = t519 * t447 - t516 * t448;
t488 = (-t516 * t517 + t519 * t520) * qJD(1);
t458 = t488 * qJD(4) + t519 * t496 + t516 * t497;
t489 = (t516 * t520 + t517 * t519) * qJD(1);
t471 = -t488 * mrSges(6,1) + t489 * mrSges(6,2);
t472 = -t488 * mrSges(5,1) + t489 * mrSges(5,2);
t511 = qJD(3) + qJD(4);
t477 = -t511 * mrSges(5,2) + t488 * mrSges(5,3);
t510 = qJDD(3) + qJDD(4);
t437 = -0.2e1 * qJD(5) * t489 + (t488 * t511 - t458) * qJ(5) + (t488 * t489 + t510) * pkin(4) + t442;
t476 = -t511 * mrSges(6,2) + t488 * mrSges(6,3);
t535 = m(6) * t437 + t510 * mrSges(6,1) + t511 * t476;
t431 = m(5) * t442 + t510 * mrSges(5,1) + t511 * t477 + (-t471 - t472) * t489 + (-mrSges(5,3) - mrSges(6,3)) * t458 + t535;
t443 = t516 * t447 + t519 * t448;
t457 = -t489 * qJD(4) - t516 * t496 + t519 * t497;
t479 = t511 * mrSges(6,1) - t489 * mrSges(6,3);
t480 = t511 * mrSges(5,1) - t489 * mrSges(5,3);
t478 = t511 * pkin(4) - t489 * qJ(5);
t481 = t488 ^ 2;
t439 = -t481 * pkin(4) + t457 * qJ(5) + 0.2e1 * qJD(5) * t488 - t511 * t478 + t443;
t534 = m(6) * t439 + t457 * mrSges(6,3) + t488 * t471;
t434 = m(5) * t443 + t457 * mrSges(5,3) + t488 * t472 + (-t479 - t480) * t511 + (-mrSges(5,2) - mrSges(6,2)) * t510 + t534;
t427 = t519 * t431 + t516 * t434;
t494 = (-mrSges(4,1) * t520 + mrSges(4,2) * t517) * qJD(1);
t537 = qJD(1) * t520;
t499 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t537;
t425 = m(4) * t455 + qJDD(3) * mrSges(4,1) - t496 * mrSges(4,3) + qJD(3) * t499 - t494 * t538 + t427;
t498 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t538;
t528 = -t516 * t431 + t519 * t434;
t426 = m(4) * t456 - qJDD(3) * mrSges(4,2) + t497 * mrSges(4,3) - qJD(3) * t498 + t494 * t537 + t528;
t529 = -t517 * t425 + t520 * t426;
t418 = m(3) * t474 - t522 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t529;
t473 = t515 * t493 - t514 * t495;
t525 = -qJDD(1) * pkin(2) - t473;
t469 = -t522 * pkin(6) + t525;
t449 = -t497 * pkin(3) + t500 * t538 + (-pkin(7) * t512 - pkin(6)) * t522 + t525;
t441 = -t457 * pkin(4) - t481 * qJ(5) + t489 * t478 + qJDD(5) + t449;
t527 = m(6) * t441 - t457 * mrSges(6,1) + t458 * mrSges(6,2) - t488 * t476 + t489 * t479;
t524 = m(5) * t449 - t457 * mrSges(5,1) + t458 * mrSges(5,2) - t488 * t477 + t489 * t480 + t527;
t523 = -m(4) * t469 + t497 * mrSges(4,1) - t496 * mrSges(4,2) - t498 * t538 + t499 * t537 - t524;
t429 = m(3) * t473 + qJDD(1) * mrSges(3,1) - t522 * mrSges(3,2) + t523;
t415 = t514 * t418 + t515 * t429;
t419 = t520 * t425 + t517 * t426;
t541 = t542 * t488 - t543 * t489 + t545 * t511;
t540 = -t546 * t488 - t544 * t489 + t542 * t511;
t539 = t544 * t488 + t547 * t489 + t543 * t511;
t533 = m(3) * t513 + t419;
t530 = t515 * t418 - t514 * t429;
t413 = m(2) * t501 - t522 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t530;
t414 = m(2) * t502 + qJDD(1) * mrSges(2,1) - t522 * mrSges(2,2) + t415;
t531 = t521 * t413 - t518 * t414;
t526 = -t518 * t413 - t521 * t414;
t487 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t517 + Ifges(4,4) * t520) * qJD(1);
t486 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t517 + Ifges(4,2) * t520) * qJD(1);
t485 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t517 + Ifges(4,6) * t520) * qJD(1);
t435 = -t458 * mrSges(6,3) - t489 * t471 + t535;
t421 = mrSges(5,2) * t449 + mrSges(6,2) * t441 - mrSges(5,3) * t442 - mrSges(6,3) * t437 - qJ(5) * t435 + t544 * t457 + t547 * t458 - t541 * t488 + t543 * t510 + t540 * t511;
t420 = -mrSges(5,1) * t449 + mrSges(5,3) * t443 - mrSges(6,1) * t441 + mrSges(6,3) * t439 - pkin(4) * t527 + qJ(5) * t534 + (-qJ(5) * t479 + t539) * t511 + (-qJ(5) * mrSges(6,2) - t542) * t510 + t541 * t489 + t544 * t458 + t546 * t457;
t411 = mrSges(4,2) * t469 - mrSges(4,3) * t455 + Ifges(4,1) * t496 + Ifges(4,4) * t497 + Ifges(4,5) * qJDD(3) - pkin(7) * t427 - qJD(3) * t486 - t516 * t420 + t519 * t421 + t485 * t537;
t410 = -mrSges(4,1) * t469 + mrSges(4,3) * t456 + Ifges(4,4) * t496 + Ifges(4,2) * t497 + Ifges(4,6) * qJDD(3) - pkin(3) * t524 + pkin(7) * t528 + qJD(3) * t487 + t519 * t420 + t516 * t421 - t485 * t538;
t409 = -Ifges(4,3) * qJDD(3) + Ifges(3,6) * qJDD(1) + t522 * Ifges(3,5) - mrSges(3,1) * t513 - Ifges(4,5) * t496 - Ifges(4,6) * t497 + mrSges(3,3) * t474 - mrSges(4,1) * t455 + mrSges(4,2) * t456 - mrSges(5,1) * t442 + mrSges(5,2) * t443 - mrSges(6,1) * t437 + mrSges(6,2) * t439 - pkin(4) * t435 - pkin(3) * t427 - pkin(2) * t419 + t545 * t510 + t540 * t489 + t539 * t488 - t543 * t458 + t542 * t457 + (-t517 * t486 + t520 * t487) * qJD(1);
t408 = mrSges(3,2) * t513 - mrSges(3,3) * t473 + Ifges(3,5) * qJDD(1) - t522 * Ifges(3,6) - pkin(6) * t419 - t517 * t410 + t520 * t411;
t407 = -mrSges(2,2) * g(1) - mrSges(2,3) * t502 + Ifges(2,5) * qJDD(1) - t522 * Ifges(2,6) - qJ(2) * t415 + t515 * t408 - t514 * t409;
t406 = mrSges(2,1) * g(1) + mrSges(2,3) * t501 + t522 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t533 + qJ(2) * t530 + t514 * t408 + t515 * t409;
t1 = [(-m(1) - m(2)) * g(1) + t533; -m(1) * g(2) + t526; -m(1) * g(3) + t531; pkin(1) * t415 + t517 * t411 + t520 * t410 + pkin(2) * t523 + pkin(6) * t529 + mrSges(3,1) * t473 - mrSges(3,2) * t474 + mrSges(2,1) * t502 - mrSges(2,2) * t501 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t531 - t521 * t406 - t518 * t407; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t526 - t518 * t406 + t521 * t407;];
tauB = t1;
