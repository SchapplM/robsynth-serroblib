% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:44
% EndTime: 2019-12-31 18:44:47
% DurationCPUTime: 2.27s
% Computational Cost: add. (20792->245), mult. (38917->296), div. (0->0), fcn. (22129->8), ass. (0->98)
t510 = Ifges(5,1) + Ifges(6,1);
t504 = Ifges(5,4) - Ifges(6,5);
t509 = -Ifges(5,5) - Ifges(6,4);
t508 = Ifges(5,2) + Ifges(6,3);
t502 = Ifges(5,6) - Ifges(6,6);
t507 = -Ifges(5,3) - Ifges(6,2);
t506 = cos(qJ(4));
t505 = -mrSges(5,3) - mrSges(6,2);
t478 = sin(qJ(1));
t480 = cos(qJ(1));
t464 = t478 * g(1) - t480 * g(2);
t456 = qJDD(1) * pkin(1) + t464;
t465 = -t480 * g(1) - t478 * g(2);
t482 = qJD(1) ^ 2;
t458 = -t482 * pkin(1) + t465;
t474 = sin(pkin(8));
t475 = cos(pkin(8));
t432 = t474 * t456 + t475 * t458;
t418 = -t482 * pkin(2) + qJDD(1) * pkin(6) + t432;
t473 = -g(3) + qJDD(2);
t477 = sin(qJ(3));
t479 = cos(qJ(3));
t414 = t479 * t418 + t477 * t473;
t457 = (-mrSges(4,1) * t479 + mrSges(4,2) * t477) * qJD(1);
t494 = qJD(1) * qJD(3);
t491 = t477 * t494;
t461 = t479 * qJDD(1) - t491;
t496 = qJD(1) * t477;
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t496;
t431 = t475 * t456 - t474 * t458;
t417 = -qJDD(1) * pkin(2) - t482 * pkin(6) - t431;
t490 = t479 * t494;
t460 = t477 * qJDD(1) + t490;
t409 = (-t460 - t490) * pkin(7) + (-t461 + t491) * pkin(3) + t417;
t459 = (-pkin(3) * t479 - pkin(7) * t477) * qJD(1);
t481 = qJD(3) ^ 2;
t495 = t479 * qJD(1);
t412 = -pkin(3) * t481 + qJDD(3) * pkin(7) + t459 * t495 + t414;
t476 = sin(qJ(4));
t407 = t476 * t409 + t506 * t412;
t455 = t476 * qJD(3) + t506 * t496;
t429 = t455 * qJD(4) - t506 * qJDD(3) + t476 * t460;
t467 = qJD(4) - t495;
t439 = t467 * mrSges(5,1) - t455 * mrSges(5,3);
t453 = qJDD(4) - t461;
t454 = -t506 * qJD(3) + t476 * t496;
t435 = t454 * pkin(4) - t455 * qJ(5);
t466 = t467 ^ 2;
t403 = -pkin(4) * t466 + qJ(5) * t453 + 0.2e1 * qJD(5) * t467 - t435 * t454 + t407;
t440 = -t467 * mrSges(6,1) + t455 * mrSges(6,2);
t493 = m(6) * t403 + t453 * mrSges(6,3) + t467 * t440;
t436 = t454 * mrSges(6,1) - t455 * mrSges(6,3);
t497 = -t454 * mrSges(5,1) - t455 * mrSges(5,2) - t436;
t399 = m(5) * t407 - t453 * mrSges(5,2) + t505 * t429 - t467 * t439 + t497 * t454 + t493;
t406 = t506 * t409 - t476 * t412;
t430 = -t454 * qJD(4) + t476 * qJDD(3) + t506 * t460;
t438 = -t467 * mrSges(5,2) - t454 * mrSges(5,3);
t404 = -t453 * pkin(4) - t466 * qJ(5) + t455 * t435 + qJDD(5) - t406;
t441 = -t454 * mrSges(6,2) + t467 * mrSges(6,3);
t485 = -m(6) * t404 + t453 * mrSges(6,1) + t467 * t441;
t400 = m(5) * t406 + t453 * mrSges(5,1) + t505 * t430 + t467 * t438 + t497 * t455 + t485;
t486 = t506 * t399 - t400 * t476;
t394 = m(4) * t414 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t461 - qJD(3) * t462 + t457 * t495 + t486;
t413 = -t477 * t418 + t479 * t473;
t463 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t495;
t411 = -qJDD(3) * pkin(3) - t481 * pkin(7) + t459 * t496 - t413;
t405 = -0.2e1 * qJD(5) * t455 + (t454 * t467 - t430) * qJ(5) + (t455 * t467 + t429) * pkin(4) + t411;
t401 = m(6) * t405 + mrSges(6,1) * t429 - t430 * mrSges(6,3) - t455 * t440 + t441 * t454;
t483 = -m(5) * t411 - t429 * mrSges(5,1) - mrSges(5,2) * t430 - t454 * t438 - t439 * t455 - t401;
t397 = m(4) * t413 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t460 + qJD(3) * t463 - t457 * t496 + t483;
t487 = t479 * t394 - t397 * t477;
t386 = m(3) * t432 - mrSges(3,1) * t482 - qJDD(1) * mrSges(3,2) + t487;
t395 = t476 * t399 + t506 * t400;
t484 = -m(4) * t417 + t461 * mrSges(4,1) - t460 * mrSges(4,2) - t462 * t496 + t463 * t495 - t395;
t389 = m(3) * t431 + qJDD(1) * mrSges(3,1) - t482 * mrSges(3,2) + t484;
t382 = t474 * t386 + t475 * t389;
t380 = m(2) * t464 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t482 + t382;
t488 = t475 * t386 - t389 * t474;
t381 = m(2) * t465 - mrSges(2,1) * t482 - qJDD(1) * mrSges(2,2) + t488;
t501 = t480 * t380 + t478 * t381;
t387 = t477 * t394 + t479 * t397;
t500 = t508 * t454 - t504 * t455 - t502 * t467;
t499 = t502 * t454 + t509 * t455 + t507 * t467;
t498 = -t504 * t454 + t510 * t455 - t509 * t467;
t492 = m(3) * t473 + t387;
t489 = -t380 * t478 + t480 * t381;
t448 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t477 + Ifges(4,4) * t479) * qJD(1);
t447 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t477 + Ifges(4,2) * t479) * qJD(1);
t446 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t477 + Ifges(4,6) * t479) * qJD(1);
t391 = mrSges(5,2) * t411 + mrSges(6,2) * t404 - mrSges(5,3) * t406 - mrSges(6,3) * t405 - qJ(5) * t401 - t504 * t429 + t510 * t430 - t453 * t509 + t499 * t454 + t500 * t467;
t390 = -mrSges(5,1) * t411 - mrSges(6,1) * t405 + mrSges(6,2) * t403 + mrSges(5,3) * t407 - pkin(4) * t401 - t508 * t429 + t504 * t430 + t502 * t453 + t499 * t455 + t498 * t467;
t383 = Ifges(4,4) * t460 + Ifges(4,2) * t461 + Ifges(4,6) * qJDD(3) - t446 * t496 + qJD(3) * t448 - mrSges(4,1) * t417 + mrSges(4,3) * t414 - mrSges(5,1) * t406 + mrSges(5,2) * t407 + mrSges(6,1) * t404 - mrSges(6,3) * t403 - pkin(4) * t485 - qJ(5) * t493 - pkin(3) * t395 + (pkin(4) * t436 + t500) * t455 + (qJ(5) * t436 - t498) * t454 + t507 * t453 + (mrSges(6,2) * pkin(4) + t509) * t430 + (mrSges(6,2) * qJ(5) + t502) * t429;
t376 = mrSges(4,2) * t417 - mrSges(4,3) * t413 + Ifges(4,1) * t460 + Ifges(4,4) * t461 + Ifges(4,5) * qJDD(3) - pkin(7) * t395 - qJD(3) * t447 - t476 * t390 + t506 * t391 + t446 * t495;
t375 = Ifges(3,6) * qJDD(1) + t482 * Ifges(3,5) - mrSges(3,1) * t473 + mrSges(3,3) * t432 - Ifges(4,5) * t460 - Ifges(4,6) * t461 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t413 + mrSges(4,2) * t414 - t476 * t391 - t506 * t390 - pkin(3) * t483 - pkin(7) * t486 - pkin(2) * t387 + (-t447 * t477 + t448 * t479) * qJD(1);
t374 = mrSges(3,2) * t473 - mrSges(3,3) * t431 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t482 - pkin(6) * t387 + t376 * t479 - t383 * t477;
t373 = -mrSges(2,2) * g(3) - mrSges(2,3) * t464 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t482 - qJ(2) * t382 + t374 * t475 - t375 * t474;
t372 = mrSges(2,1) * g(3) + mrSges(2,3) * t465 + t482 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t492 + qJ(2) * t488 + t474 * t374 + t475 * t375;
t1 = [-m(1) * g(1) + t489; -m(1) * g(2) + t501; (-m(1) - m(2)) * g(3) + t492; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t501 - t478 * t372 + t480 * t373; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t489 + t480 * t372 + t478 * t373; pkin(1) * t382 + mrSges(2,1) * t464 - mrSges(2,2) * t465 + t477 * t376 + t479 * t383 + pkin(2) * t484 + pkin(6) * t487 + mrSges(3,1) * t431 - mrSges(3,2) * t432 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
