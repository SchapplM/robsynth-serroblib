% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:14
% EndTime: 2019-12-05 16:12:20
% DurationCPUTime: 2.32s
% Computational Cost: add. (15600->239), mult. (31786->289), div. (0->0), fcn. (18640->8), ass. (0->99)
t500 = Ifges(5,1) + Ifges(6,1);
t492 = Ifges(5,4) - Ifges(6,5);
t499 = Ifges(5,5) + Ifges(6,4);
t498 = -Ifges(5,2) - Ifges(6,3);
t497 = Ifges(6,2) + Ifges(5,3);
t496 = Ifges(5,6) - Ifges(6,6);
t495 = -2 * qJD(4);
t494 = -2 * qJD(5);
t493 = -mrSges(5,3) - mrSges(6,2);
t489 = cos(pkin(7));
t488 = cos(pkin(8));
t462 = cos(qJ(3));
t465 = qJD(2) ^ 2;
t487 = t462 ^ 2 * t465;
t459 = sin(pkin(7));
t447 = t459 * g(1) - t489 * g(2);
t486 = t462 * t447;
t448 = -t489 * g(1) - t459 * g(2);
t457 = -g(3) + qJDD(1);
t461 = sin(qJ(2));
t463 = cos(qJ(2));
t428 = t463 * t448 + t461 * t457;
t418 = -t465 * pkin(2) + qJDD(2) * pkin(6) + t428;
t460 = sin(qJ(3));
t403 = t462 * t418 - t460 * t447;
t444 = (-mrSges(4,1) * t462 + mrSges(4,2) * t460) * qJD(2);
t477 = qJD(2) * qJD(3);
t476 = t460 * t477;
t446 = t462 * qJDD(2) - t476;
t479 = qJD(2) * t460;
t449 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t479;
t427 = -t461 * t448 + t463 * t457;
t417 = -qJDD(2) * pkin(2) - t465 * pkin(6) - t427;
t475 = t462 * t477;
t445 = t460 * qJDD(2) + t475;
t400 = (-t445 - t475) * qJ(4) + (-t446 + t476) * pkin(3) + t417;
t443 = (-pkin(3) * t462 - qJ(4) * t460) * qJD(2);
t464 = qJD(3) ^ 2;
t478 = qJD(2) * t462;
t401 = -t464 * pkin(3) + qJDD(3) * qJ(4) + t443 * t478 + t403;
t458 = sin(pkin(8));
t435 = -t488 * qJD(3) + t458 * t479;
t396 = t458 * t400 + t488 * t401 + t435 * t495;
t436 = t458 * qJD(3) + t488 * t479;
t423 = -mrSges(5,1) * t478 - t436 * mrSges(5,3);
t424 = mrSges(6,1) * t478 + t436 * mrSges(6,2);
t425 = -t488 * qJDD(3) + t458 * t445;
t413 = t435 * mrSges(6,1) - t436 * mrSges(6,3);
t480 = -t435 * mrSges(5,1) - t436 * mrSges(5,2) - t413;
t412 = t435 * pkin(4) - t436 * qJ(5);
t392 = -pkin(4) * t487 - t446 * qJ(5) - t435 * t412 + t478 * t494 + t396;
t484 = m(6) * t392 - t446 * mrSges(6,3);
t387 = m(5) * t396 + t446 * mrSges(5,2) + t480 * t435 + t493 * t425 + (t423 - t424) * t478 + t484;
t469 = t488 * t400 - t458 * t401;
t395 = t436 * t495 + t469;
t421 = -t435 * mrSges(6,2) - mrSges(6,3) * t478;
t422 = mrSges(5,2) * t478 - t435 * mrSges(5,3);
t426 = t458 * qJDD(3) + t488 * t445;
t393 = -qJ(5) * t487 + t446 * pkin(4) + qJDD(5) + ((2 * qJD(4)) + t412) * t436 - t469;
t474 = -m(6) * t393 - t446 * mrSges(6,1);
t388 = m(5) * t395 - t446 * mrSges(5,1) + t480 * t436 + t493 * t426 + (-t421 - t422) * t478 + t474;
t470 = t488 * t387 - t458 * t388;
t382 = m(4) * t403 - qJDD(3) * mrSges(4,2) + t446 * mrSges(4,3) - qJD(3) * t449 + t444 * t478 + t470;
t415 = t460 * t418;
t402 = -t415 - t486;
t450 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t478;
t468 = -qJDD(3) * pkin(3) - t464 * qJ(4) + t443 * t479 + qJDD(4) + t415;
t394 = t425 * pkin(4) - t426 * qJ(5) + t436 * t494 + (t447 + (-pkin(4) * t436 - qJ(5) * t435) * qJD(2)) * t462 + t468;
t390 = m(6) * t394 + t425 * mrSges(6,1) - t426 * mrSges(6,3) + t435 * t421 - t436 * t424;
t399 = t468 + t486;
t466 = -m(5) * t399 - t425 * mrSges(5,1) - t426 * mrSges(5,2) - t435 * t422 - t436 * t423 - t390;
t389 = m(4) * t402 + qJDD(3) * mrSges(4,1) - t445 * mrSges(4,3) + qJD(3) * t450 - t444 * t479 + t466;
t471 = t462 * t382 - t389 * t460;
t375 = m(3) * t428 - mrSges(3,1) * t465 - qJDD(2) * mrSges(3,2) + t471;
t385 = t458 * t387 + t488 * t388;
t467 = -m(4) * t417 + t446 * mrSges(4,1) - t445 * mrSges(4,2) - t449 * t479 + t450 * t478 - t385;
t380 = m(3) * t427 + qJDD(2) * mrSges(3,1) - t465 * mrSges(3,2) + t467;
t472 = t463 * t375 - t380 * t461;
t370 = m(2) * t448 + t472;
t378 = t460 * t382 + t462 * t389;
t377 = (m(2) + m(3)) * t447 - t378;
t485 = t459 * t370 + t489 * t377;
t371 = t461 * t375 + t463 * t380;
t483 = t498 * t435 + t492 * t436 - t496 * t478;
t482 = t496 * t435 - t499 * t436 + t497 * t478;
t481 = t492 * t435 - t500 * t436 + t499 * t478;
t473 = t489 * t370 - t377 * t459;
t433 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t460 + Ifges(4,4) * t462) * qJD(2);
t432 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t460 + Ifges(4,2) * t462) * qJD(2);
t431 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t460 + Ifges(4,6) * t462) * qJD(2);
t384 = mrSges(5,2) * t399 + mrSges(6,2) * t393 - mrSges(5,3) * t395 - mrSges(6,3) * t394 - qJ(5) * t390 - t492 * t425 + t500 * t426 + t482 * t435 - t446 * t499 + t483 * t478;
t383 = -mrSges(5,1) * t399 - mrSges(6,1) * t394 + mrSges(6,2) * t392 + mrSges(5,3) * t396 - pkin(4) * t390 + t498 * t425 + t492 * t426 + t482 * t436 - t446 * t496 + t481 * t478;
t372 = Ifges(4,4) * t445 + Ifges(4,6) * qJDD(3) - t431 * t479 + qJD(3) * t433 - mrSges(4,1) * t417 + mrSges(4,3) * t403 - mrSges(5,1) * t395 + mrSges(5,2) * t396 + mrSges(6,1) * t393 - mrSges(6,3) * t392 - pkin(4) * (-t421 * t478 + t474) - qJ(5) * (-t424 * t478 + t484) - pkin(3) * t385 + (pkin(4) * t413 - t483) * t436 + (qJ(5) * t413 + t481) * t435 + (mrSges(6,2) * pkin(4) - t499) * t426 + (mrSges(6,2) * qJ(5) + t496) * t425 + (Ifges(4,2) + t497) * t446;
t367 = mrSges(4,2) * t417 - mrSges(4,3) * t402 + Ifges(4,1) * t445 + Ifges(4,4) * t446 + Ifges(4,5) * qJDD(3) - qJ(4) * t385 - qJD(3) * t432 - t458 * t383 + t488 * t384 + t431 * t478;
t366 = Ifges(3,6) * qJDD(2) + t465 * Ifges(3,5) + mrSges(3,1) * t447 + mrSges(3,3) * t428 - Ifges(4,5) * t445 - Ifges(4,6) * t446 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t402 + mrSges(4,2) * t403 - t458 * t384 - t488 * t383 - pkin(3) * t466 - qJ(4) * t470 - pkin(2) * t378 + (-t432 * t460 + t433 * t462) * qJD(2);
t365 = -mrSges(3,2) * t447 - mrSges(3,3) * t427 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t465 - pkin(6) * t378 + t367 * t462 - t372 * t460;
t364 = -mrSges(2,1) * t457 - mrSges(3,1) * t427 + mrSges(3,2) * t428 + mrSges(2,3) * t448 - Ifges(3,3) * qJDD(2) - pkin(1) * t371 - pkin(2) * t467 - pkin(6) * t471 - t460 * t367 - t462 * t372;
t363 = mrSges(2,2) * t457 - mrSges(2,3) * t447 - pkin(5) * t371 + t365 * t463 - t366 * t461;
t1 = [-m(1) * g(1) + t473; -m(1) * g(2) + t485; -m(1) * g(3) + m(2) * t457 + t371; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t485 + t489 * t363 - t459 * t364; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t473 + t459 * t363 + t489 * t364; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t447 - mrSges(2,2) * t448 + t461 * t365 + t463 * t366 + pkin(1) * (m(3) * t447 - t378) + pkin(5) * t472;];
tauB = t1;
