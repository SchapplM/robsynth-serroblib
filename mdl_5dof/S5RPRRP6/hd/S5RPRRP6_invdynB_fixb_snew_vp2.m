% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:20
% EndTime: 2019-12-31 18:42:23
% DurationCPUTime: 2.29s
% Computational Cost: add. (21346->247), mult. (40307->296), div. (0->0), fcn. (23079->8), ass. (0->99)
t512 = Ifges(5,1) + Ifges(6,1);
t507 = Ifges(5,4) + Ifges(6,4);
t506 = Ifges(5,5) + Ifges(6,5);
t511 = Ifges(5,2) + Ifges(6,2);
t510 = Ifges(5,6) + Ifges(6,6);
t509 = Ifges(5,3) + Ifges(6,3);
t508 = -mrSges(5,2) - mrSges(6,2);
t479 = sin(qJ(1));
t482 = cos(qJ(1));
t468 = t479 * g(1) - t482 * g(2);
t460 = qJDD(1) * pkin(1) + t468;
t469 = -t482 * g(1) - t479 * g(2);
t484 = qJD(1) ^ 2;
t462 = -t484 * pkin(1) + t469;
t475 = sin(pkin(8));
t476 = cos(pkin(8));
t437 = t475 * t460 + t476 * t462;
t422 = -t484 * pkin(2) + qJDD(1) * pkin(6) + t437;
t474 = -g(3) + qJDD(2);
t478 = sin(qJ(3));
t481 = cos(qJ(3));
t417 = t481 * t422 + t478 * t474;
t461 = (-mrSges(4,1) * t481 + mrSges(4,2) * t478) * qJD(1);
t497 = qJD(1) * qJD(3);
t493 = t478 * t497;
t465 = t481 * qJDD(1) - t493;
t499 = qJD(1) * t478;
t466 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t499;
t436 = t476 * t460 - t475 * t462;
t421 = -qJDD(1) * pkin(2) - t484 * pkin(6) - t436;
t492 = t481 * t497;
t464 = t478 * qJDD(1) + t492;
t412 = (-t464 - t492) * pkin(7) + (-t465 + t493) * pkin(3) + t421;
t463 = (-pkin(3) * t481 - pkin(7) * t478) * qJD(1);
t483 = qJD(3) ^ 2;
t498 = t481 * qJD(1);
t415 = -t483 * pkin(3) + qJDD(3) * pkin(7) + t463 * t498 + t417;
t477 = sin(qJ(4));
t480 = cos(qJ(4));
t408 = t480 * t412 - t477 * t415;
t458 = t480 * qJD(3) - t477 * t499;
t435 = t458 * qJD(4) + t477 * qJDD(3) + t480 * t464;
t459 = t477 * qJD(3) + t480 * t499;
t439 = -t458 * mrSges(6,1) + t459 * mrSges(6,2);
t440 = -t458 * mrSges(5,1) + t459 * mrSges(5,2);
t470 = qJD(4) - t498;
t442 = -t470 * mrSges(5,2) + t458 * mrSges(5,3);
t457 = qJDD(4) - t465;
t404 = -0.2e1 * qJD(5) * t459 + (t458 * t470 - t435) * qJ(5) + (t458 * t459 + t457) * pkin(4) + t408;
t441 = -t470 * mrSges(6,2) + t458 * mrSges(6,3);
t496 = m(6) * t404 + t457 * mrSges(6,1) + t470 * t441;
t397 = m(5) * t408 + t457 * mrSges(5,1) + t470 * t442 + (-t439 - t440) * t459 + (-mrSges(5,3) - mrSges(6,3)) * t435 + t496;
t409 = t477 * t412 + t480 * t415;
t434 = -t459 * qJD(4) + t480 * qJDD(3) - t477 * t464;
t443 = t470 * pkin(4) - t459 * qJ(5);
t456 = t458 ^ 2;
t406 = -t456 * pkin(4) + t434 * qJ(5) + 0.2e1 * qJD(5) * t458 - t470 * t443 + t409;
t495 = m(6) * t406 + t434 * mrSges(6,3) + t458 * t439;
t444 = t470 * mrSges(6,1) - t459 * mrSges(6,3);
t500 = -t470 * mrSges(5,1) + t459 * mrSges(5,3) - t444;
t399 = m(5) * t409 + t434 * mrSges(5,3) + t458 * t440 + t508 * t457 + t500 * t470 + t495;
t488 = -t477 * t397 + t480 * t399;
t394 = m(4) * t417 - qJDD(3) * mrSges(4,2) + t465 * mrSges(4,3) - qJD(3) * t466 + t461 * t498 + t488;
t416 = -t478 * t422 + t481 * t474;
t467 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t498;
t414 = -qJDD(3) * pkin(3) - t483 * pkin(7) + t463 * t499 - t416;
t407 = -t434 * pkin(4) - t456 * qJ(5) + t459 * t443 + qJDD(5) + t414;
t487 = m(6) * t407 - t434 * mrSges(6,1) - t458 * t441;
t485 = -m(5) * t414 + t434 * mrSges(5,1) + t508 * t435 + t458 * t442 + t500 * t459 - t487;
t401 = m(4) * t416 + qJDD(3) * mrSges(4,1) - t464 * mrSges(4,3) + qJD(3) * t467 - t461 * t499 + t485;
t489 = t481 * t394 - t478 * t401;
t387 = m(3) * t437 - t484 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t489;
t396 = t480 * t397 + t477 * t399;
t486 = -m(4) * t421 + t465 * mrSges(4,1) - t464 * mrSges(4,2) - t466 * t499 + t467 * t498 - t396;
t391 = m(3) * t436 + qJDD(1) * mrSges(3,1) - t484 * mrSges(3,2) + t486;
t383 = t475 * t387 + t476 * t391;
t381 = m(2) * t468 + qJDD(1) * mrSges(2,1) - t484 * mrSges(2,2) + t383;
t490 = t476 * t387 - t475 * t391;
t382 = m(2) * t469 - t484 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t490;
t504 = t482 * t381 + t479 * t382;
t388 = t478 * t394 + t481 * t401;
t503 = t510 * t458 + t506 * t459 + t509 * t470;
t502 = -t511 * t458 - t507 * t459 - t510 * t470;
t501 = t507 * t458 + t512 * t459 + t506 * t470;
t494 = m(3) * t474 + t388;
t491 = -t479 * t381 + t482 * t382;
t452 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t478 + Ifges(4,4) * t481) * qJD(1);
t451 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t478 + Ifges(4,2) * t481) * qJD(1);
t450 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t478 + Ifges(4,6) * t481) * qJD(1);
t402 = -t435 * mrSges(6,3) - t459 * t439 + t496;
t395 = mrSges(5,2) * t414 + mrSges(6,2) * t407 - mrSges(5,3) * t408 - mrSges(6,3) * t404 - qJ(5) * t402 + t507 * t434 + t512 * t435 + t506 * t457 + t503 * t458 + t502 * t470;
t389 = -mrSges(5,1) * t414 + mrSges(5,3) * t409 - mrSges(6,1) * t407 + mrSges(6,3) * t406 - pkin(4) * t487 + qJ(5) * t495 + (-qJ(5) * t444 + t501) * t470 + (-pkin(4) * t444 - t503) * t459 + (-qJ(5) * mrSges(6,2) + t510) * t457 + (-pkin(4) * mrSges(6,2) + t507) * t435 + t511 * t434;
t384 = -t450 * t499 - mrSges(4,1) * t421 - mrSges(5,1) * t408 - mrSges(6,1) * t404 + mrSges(5,2) * t409 + mrSges(6,2) * t406 + mrSges(4,3) * t417 + Ifges(4,4) * t464 + Ifges(4,2) * t465 + Ifges(4,6) * qJDD(3) - pkin(3) * t396 - pkin(4) * t402 + qJD(3) * t452 + t502 * t459 + t501 * t458 - t509 * t457 - t506 * t435 - t510 * t434;
t377 = mrSges(4,2) * t421 - mrSges(4,3) * t416 + Ifges(4,1) * t464 + Ifges(4,4) * t465 + Ifges(4,5) * qJDD(3) - pkin(7) * t396 - qJD(3) * t451 - t477 * t389 + t480 * t395 + t450 * t498;
t376 = Ifges(3,6) * qJDD(1) + t484 * Ifges(3,5) - mrSges(3,1) * t474 + mrSges(3,3) * t437 - Ifges(4,5) * t464 - Ifges(4,6) * t465 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t416 + mrSges(4,2) * t417 - t477 * t395 - t480 * t389 - pkin(3) * t485 - pkin(7) * t488 - pkin(2) * t388 + (-t478 * t451 + t481 * t452) * qJD(1);
t375 = mrSges(3,2) * t474 - mrSges(3,3) * t436 + Ifges(3,5) * qJDD(1) - t484 * Ifges(3,6) - pkin(6) * t388 + t481 * t377 - t478 * t384;
t374 = -mrSges(2,2) * g(3) - mrSges(2,3) * t468 + Ifges(2,5) * qJDD(1) - t484 * Ifges(2,6) - qJ(2) * t383 + t476 * t375 - t475 * t376;
t373 = mrSges(2,1) * g(3) + mrSges(2,3) * t469 + t484 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t494 + qJ(2) * t490 + t475 * t375 + t476 * t376;
t1 = [-m(1) * g(1) + t491; -m(1) * g(2) + t504; (-m(1) - m(2)) * g(3) + t494; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t504 - t479 * t373 + t482 * t374; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t491 + t482 * t373 + t479 * t374; pkin(1) * t383 + mrSges(2,1) * t468 - mrSges(2,2) * t469 + t478 * t377 + t481 * t384 + pkin(2) * t486 + pkin(6) * t489 + mrSges(3,1) * t436 - mrSges(3,2) * t437 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
