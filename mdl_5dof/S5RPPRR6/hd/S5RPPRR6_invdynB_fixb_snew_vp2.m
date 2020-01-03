% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR6
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:54
% EndTime: 2019-12-31 17:57:57
% DurationCPUTime: 3.36s
% Computational Cost: add. (36696->248), mult. (80060->310), div. (0->0), fcn. (52533->10), ass. (0->109)
t499 = qJD(1) ^ 2;
t488 = sin(pkin(9));
t490 = cos(pkin(9));
t493 = sin(qJ(4));
t496 = cos(qJ(4));
t504 = t488 * t493 - t490 * t496;
t465 = t504 * qJD(1);
t505 = t488 * t496 + t490 * t493;
t466 = t505 * qJD(1);
t518 = t466 * qJD(4);
t456 = -t504 * qJDD(1) - t518;
t525 = pkin(3) * t490;
t524 = mrSges(4,2) * t488;
t486 = t490 ^ 2;
t523 = t486 * t499;
t494 = sin(qJ(1));
t497 = cos(qJ(1));
t474 = t494 * g(1) - t497 * g(2);
t472 = qJDD(1) * pkin(1) + t474;
t475 = -t497 * g(1) - t494 * g(2);
t473 = -t499 * pkin(1) + t475;
t489 = sin(pkin(8));
t491 = cos(pkin(8));
t459 = t489 * t472 + t491 * t473;
t450 = -t499 * pkin(2) + qJDD(1) * qJ(3) + t459;
t487 = -g(3) + qJDD(2);
t517 = qJD(1) * qJD(3);
t521 = t490 * t487 - 0.2e1 * t488 * t517;
t435 = (-pkin(6) * qJDD(1) + t499 * t525 - t450) * t488 + t521;
t441 = t488 * t487 + (t450 + 0.2e1 * t517) * t490;
t516 = qJDD(1) * t490;
t438 = -pkin(3) * t523 + pkin(6) * t516 + t441;
t427 = t493 * t435 + t496 * t438;
t452 = t465 * mrSges(5,1) + t466 * mrSges(5,2);
t463 = qJD(4) * mrSges(5,1) - t466 * mrSges(5,3);
t455 = t465 * pkin(4) - t466 * pkin(7);
t498 = qJD(4) ^ 2;
t424 = -t498 * pkin(4) + qJDD(4) * pkin(7) - t465 * t455 + t427;
t485 = t488 ^ 2;
t458 = t491 * t472 - t489 * t473;
t506 = qJDD(3) - t458;
t439 = (-pkin(2) - t525) * qJDD(1) + (-qJ(3) + (-t485 - t486) * pkin(6)) * t499 + t506;
t519 = t465 * qJD(4);
t457 = t505 * qJDD(1) - t519;
t425 = (-t457 + t519) * pkin(7) + (-t456 + t518) * pkin(4) + t439;
t492 = sin(qJ(5));
t495 = cos(qJ(5));
t421 = -t492 * t424 + t495 * t425;
t460 = t495 * qJD(4) - t492 * t466;
t437 = t460 * qJD(5) + t492 * qJDD(4) + t495 * t457;
t461 = t492 * qJD(4) + t495 * t466;
t442 = -t460 * mrSges(6,1) + t461 * mrSges(6,2);
t464 = qJD(5) + t465;
t443 = -t464 * mrSges(6,2) + t460 * mrSges(6,3);
t454 = qJDD(5) - t456;
t419 = m(6) * t421 + t454 * mrSges(6,1) - t437 * mrSges(6,3) - t461 * t442 + t464 * t443;
t422 = t495 * t424 + t492 * t425;
t436 = -t461 * qJD(5) + t495 * qJDD(4) - t492 * t457;
t444 = t464 * mrSges(6,1) - t461 * mrSges(6,3);
t420 = m(6) * t422 - t454 * mrSges(6,2) + t436 * mrSges(6,3) + t460 * t442 - t464 * t444;
t510 = -t492 * t419 + t495 * t420;
t410 = m(5) * t427 - qJDD(4) * mrSges(5,2) + t456 * mrSges(5,3) - qJD(4) * t463 - t465 * t452 + t510;
t426 = t496 * t435 - t493 * t438;
t462 = -qJD(4) * mrSges(5,2) - t465 * mrSges(5,3);
t423 = -qJDD(4) * pkin(4) - t498 * pkin(7) + t466 * t455 - t426;
t502 = -m(6) * t423 + t436 * mrSges(6,1) - t437 * mrSges(6,2) + t460 * t443 - t461 * t444;
t415 = m(5) * t426 + qJDD(4) * mrSges(5,1) - t457 * mrSges(5,3) + qJD(4) * t462 - t466 * t452 + t502;
t405 = t493 * t410 + t496 * t415;
t440 = -t488 * t450 + t521;
t503 = mrSges(4,3) * qJDD(1) + t499 * (-mrSges(4,1) * t490 + t524);
t403 = m(4) * t440 - t503 * t488 + t405;
t511 = t496 * t410 - t493 * t415;
t404 = m(4) * t441 + t503 * t490 + t511;
t512 = -t488 * t403 + t490 * t404;
t396 = m(3) * t459 - t499 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t512;
t446 = -qJDD(1) * pkin(2) - t499 * qJ(3) + t506;
t411 = t495 * t419 + t492 * t420;
t501 = m(5) * t439 - t456 * mrSges(5,1) + t457 * mrSges(5,2) + t465 * t462 + t466 * t463 + t411;
t500 = -m(4) * t446 + mrSges(4,1) * t516 - t501 + (t485 * t499 + t523) * mrSges(4,3);
t407 = -t499 * mrSges(3,2) + t500 + m(3) * t458 + (mrSges(3,1) - t524) * qJDD(1);
t393 = t489 * t396 + t491 * t407;
t391 = m(2) * t474 + qJDD(1) * mrSges(2,1) - t499 * mrSges(2,2) + t393;
t513 = t491 * t396 - t489 * t407;
t392 = m(2) * t475 - t499 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t513;
t522 = t497 * t391 + t494 * t392;
t397 = t490 * t403 + t488 * t404;
t507 = Ifges(4,5) * t488 + Ifges(4,6) * t490;
t520 = t499 * t507;
t515 = m(3) * t487 + t397;
t514 = -t494 * t391 + t497 * t392;
t509 = Ifges(4,1) * t488 + Ifges(4,4) * t490;
t508 = Ifges(4,4) * t488 + Ifges(4,2) * t490;
t449 = Ifges(5,1) * t466 - Ifges(5,4) * t465 + Ifges(5,5) * qJD(4);
t448 = Ifges(5,4) * t466 - Ifges(5,2) * t465 + Ifges(5,6) * qJD(4);
t447 = Ifges(5,5) * t466 - Ifges(5,6) * t465 + Ifges(5,3) * qJD(4);
t431 = Ifges(6,1) * t461 + Ifges(6,4) * t460 + Ifges(6,5) * t464;
t430 = Ifges(6,4) * t461 + Ifges(6,2) * t460 + Ifges(6,6) * t464;
t429 = Ifges(6,5) * t461 + Ifges(6,6) * t460 + Ifges(6,3) * t464;
t413 = mrSges(6,2) * t423 - mrSges(6,3) * t421 + Ifges(6,1) * t437 + Ifges(6,4) * t436 + Ifges(6,5) * t454 + t460 * t429 - t464 * t430;
t412 = -mrSges(6,1) * t423 + mrSges(6,3) * t422 + Ifges(6,4) * t437 + Ifges(6,2) * t436 + Ifges(6,6) * t454 - t461 * t429 + t464 * t431;
t399 = -mrSges(5,1) * t439 - mrSges(6,1) * t421 + mrSges(6,2) * t422 + mrSges(5,3) * t427 + Ifges(5,4) * t457 - Ifges(6,5) * t437 + Ifges(5,2) * t456 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t436 - Ifges(6,3) * t454 - pkin(4) * t411 + qJD(4) * t449 - t461 * t430 + t460 * t431 - t466 * t447;
t398 = mrSges(5,2) * t439 - mrSges(5,3) * t426 + Ifges(5,1) * t457 + Ifges(5,4) * t456 + Ifges(5,5) * qJDD(4) - pkin(7) * t411 - qJD(4) * t448 - t492 * t412 + t495 * t413 - t465 * t447;
t387 = mrSges(4,2) * t446 - mrSges(4,3) * t440 - pkin(6) * t405 + t509 * qJDD(1) + t496 * t398 - t493 * t399 + t490 * t520;
t386 = -mrSges(4,1) * t446 + mrSges(4,3) * t441 - pkin(3) * t501 + pkin(6) * t511 + t508 * qJDD(1) + t493 * t398 + t496 * t399 - t488 * t520;
t385 = -pkin(2) * t397 + mrSges(3,3) * t459 - mrSges(3,1) * t487 - pkin(3) * t405 - mrSges(4,1) * t440 + mrSges(4,2) * t441 - t492 * t413 - t495 * t412 - pkin(4) * t502 - pkin(7) * t510 - Ifges(5,5) * t457 - Ifges(5,6) * t456 - Ifges(5,3) * qJDD(4) - t466 * t448 - t465 * t449 - mrSges(5,1) * t426 + mrSges(5,2) * t427 + (Ifges(3,6) - t507) * qJDD(1) + (-t488 * t508 + t490 * t509 + Ifges(3,5)) * t499;
t384 = mrSges(3,2) * t487 - mrSges(3,3) * t458 + Ifges(3,5) * qJDD(1) - t499 * Ifges(3,6) - qJ(3) * t397 - t488 * t386 + t490 * t387;
t383 = -mrSges(2,2) * g(3) - mrSges(2,3) * t474 + Ifges(2,5) * qJDD(1) - t499 * Ifges(2,6) - qJ(2) * t393 + t491 * t384 - t489 * t385;
t382 = mrSges(2,1) * g(3) + mrSges(2,3) * t475 + t499 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t515 + qJ(2) * t513 + t489 * t384 + t491 * t385;
t1 = [-m(1) * g(1) + t514; -m(1) * g(2) + t522; (-m(1) - m(2)) * g(3) + t515; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t522 - t494 * t382 + t497 * t383; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t514 + t497 * t382 + t494 * t383; pkin(1) * t393 + mrSges(2,1) * t474 - mrSges(2,2) * t475 + t488 * t387 + t490 * t386 + pkin(2) * t500 + qJ(3) * t512 + mrSges(3,1) * t458 - mrSges(3,2) * t459 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t524 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
