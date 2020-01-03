% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR9
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:49
% EndTime: 2019-12-31 19:40:51
% DurationCPUTime: 1.43s
% Computational Cost: add. (3106->244), mult. (6451->287), div. (0->0), fcn. (3001->6), ass. (0->103)
t511 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t489 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t488 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t510 = Ifges(3,2) + Ifges(5,1) + Ifges(4,3);
t487 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t509 = Ifges(3,3) + Ifges(4,2) + Ifges(5,3);
t462 = sin(qJ(2));
t493 = qJD(1) * t462;
t468 = qJD(1) ^ 2;
t463 = sin(qJ(1));
t466 = cos(qJ(1));
t480 = -g(1) * t466 - g(2) * t463;
t416 = -pkin(1) * t468 + qJDD(1) * pkin(6) + t480;
t465 = cos(qJ(2));
t399 = -g(3) * t462 + t465 * t416;
t427 = (-pkin(2) * t465 - qJ(3) * t462) * qJD(1);
t492 = qJD(1) * t465;
t508 = qJDD(2) * qJ(3) + t427 * t492 + t399;
t491 = qJD(1) * qJD(2);
t483 = t465 * t491;
t432 = qJDD(1) * t462 + t483;
t490 = qJD(1) * qJD(4);
t507 = -0.2e1 * t462 * t490 + (-t432 + t483) * qJ(4);
t458 = t465 ^ 2;
t506 = -2 * qJD(3);
t505 = 2 * qJD(3);
t504 = -pkin(2) - pkin(7);
t503 = pkin(3) + pkin(7);
t467 = qJD(2) ^ 2;
t502 = pkin(2) * t467;
t501 = -mrSges(4,1) + mrSges(5,2);
t500 = mrSges(4,2) - mrSges(5,3);
t498 = t458 * t468;
t497 = t465 * t468;
t482 = t462 * t491;
t433 = qJDD(1) * t465 - t482;
t437 = -qJD(2) * pkin(3) - qJ(4) * t493;
t494 = t463 * g(1) - t466 * g(2);
t415 = -qJDD(1) * pkin(1) - t468 * pkin(6) - t494;
t478 = -t433 * pkin(2) + t415 + (-t432 - t483) * qJ(3);
t471 = -qJ(4) * t498 + qJDD(4) - t478 + (t437 + t505) * t493;
t374 = t471 + pkin(4) * t432 + t503 * t433 + (pkin(4) * t465 + t504 * t462) * t491;
t431 = (pkin(4) * t462 + pkin(7) * t465) * qJD(1);
t398 = -t465 * g(3) - t462 * t416;
t481 = t427 * t493 + qJDD(3) - t398;
t377 = (-pkin(4) - qJ(3)) * t467 + (-pkin(3) * t497 - qJD(1) * t431) * t462 + (-pkin(2) - t503) * qJDD(2) + t481 + t507;
t461 = sin(qJ(5));
t464 = cos(qJ(5));
t372 = t374 * t464 - t377 * t461;
t425 = -qJD(2) * t464 + t461 * t492;
t394 = qJD(5) * t425 - qJDD(2) * t461 - t433 * t464;
t426 = -qJD(2) * t461 - t464 * t492;
t395 = -mrSges(6,1) * t425 + mrSges(6,2) * t426;
t446 = qJD(5) + t493;
t396 = -mrSges(6,2) * t446 + mrSges(6,3) * t425;
t423 = qJDD(5) + t432;
t369 = m(6) * t372 + mrSges(6,1) * t423 - mrSges(6,3) * t394 - t395 * t426 + t396 * t446;
t373 = t374 * t461 + t377 * t464;
t393 = -qJD(5) * t426 - qJDD(2) * t464 + t433 * t461;
t397 = mrSges(6,1) * t446 - mrSges(6,3) * t426;
t370 = m(6) * t373 - mrSges(6,2) * t423 + mrSges(6,3) * t393 + t395 * t425 - t397 * t446;
t362 = t464 * t369 + t461 * t370;
t496 = -t461 * t369 + t464 * t370;
t438 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t493;
t440 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t493;
t495 = -t438 - t440;
t486 = -t509 * qJD(2) + (-t462 * t488 - t465 * t487) * qJD(1);
t485 = -t487 * qJD(2) + (-t462 * t489 - t510 * t465) * qJD(1);
t484 = t488 * qJD(2) + (t462 * t511 + t465 * t489) * qJD(1);
t378 = -pkin(2) * t482 + pkin(3) * t433 + t471;
t441 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t492;
t479 = m(5) * t378 + t432 * mrSges(5,1) - t441 * t492 + t362;
t450 = qJD(2) * t505;
t473 = pkin(3) * t498 + qJ(4) * t433 - t508;
t376 = qJDD(2) * pkin(4) + qJD(2) * t437 + t450 + t504 * t467 + (-0.2e1 * qJD(4) - t431) * t492 - t473;
t477 = -m(6) * t376 + t393 * mrSges(6,1) - t394 * mrSges(6,2) + t425 * t396 - t426 * t397;
t381 = (pkin(2) * qJD(2) + t506) * t493 + t478;
t476 = m(4) * t381 - t479;
t385 = -qJDD(2) * pkin(2) - qJ(3) * t467 + t481;
t387 = Ifges(6,4) * t426 + Ifges(6,2) * t425 + Ifges(6,6) * t446;
t388 = Ifges(6,1) * t426 + Ifges(6,4) * t425 + Ifges(6,5) * t446;
t475 = mrSges(6,1) * t372 - mrSges(6,2) * t373 + Ifges(6,5) * t394 + Ifges(6,6) * t393 + Ifges(6,3) * t423 + t426 * t387 - t425 * t388;
t380 = (-t462 * t497 - qJDD(2)) * pkin(3) + t385 + t507;
t430 = (mrSges(5,1) * t462 - mrSges(5,2) * t465) * qJD(1);
t474 = m(5) * t380 + qJDD(2) * mrSges(5,2) + qJD(2) * t441 - t430 * t493 + t496;
t379 = 0.2e1 * t465 * t490 + t502 + (t506 - t437) * qJD(2) + t473;
t472 = -m(5) * t379 + qJDD(2) * mrSges(5,1) - t433 * mrSges(5,3) + qJD(2) * t438 - t477;
t443 = mrSges(4,2) * t492 + qJD(2) * mrSges(4,3);
t470 = m(4) * t385 - qJDD(2) * mrSges(4,1) - qJD(2) * t443 + t474;
t384 = t450 - t502 + t508;
t428 = (-mrSges(4,1) * t465 - mrSges(4,3) * t462) * qJD(1);
t469 = m(4) * t384 + qJDD(2) * mrSges(4,3) + qJD(2) * t440 + t428 * t492 + t472;
t442 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t492;
t439 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t493;
t429 = (-mrSges(3,1) * t465 + mrSges(3,2) * t462) * qJD(1);
t386 = Ifges(6,5) * t426 + Ifges(6,6) * t425 + Ifges(6,3) * t446;
t364 = mrSges(6,2) * t376 - mrSges(6,3) * t372 + Ifges(6,1) * t394 + Ifges(6,4) * t393 + Ifges(6,5) * t423 + t386 * t425 - t387 * t446;
t363 = -mrSges(6,1) * t376 + mrSges(6,3) * t373 + Ifges(6,4) * t394 + Ifges(6,2) * t393 + Ifges(6,6) * t423 - t386 * t426 + t388 * t446;
t361 = -t432 * mrSges(5,3) + t474;
t360 = -mrSges(5,2) * t433 + t438 * t493 + t479;
t359 = t428 * t493 + t500 * t432 + t470;
t358 = -mrSges(4,3) * t432 + t501 * t433 + (-t443 * t465 + t495 * t462) * qJD(1) + t476;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t494 - mrSges(2,2) * t480 + t462 * (mrSges(5,1) * t378 + mrSges(3,2) * t415 + mrSges(4,2) * t385 - mrSges(3,3) * t398 - mrSges(4,3) * t381 - mrSges(5,3) * t380 + pkin(4) * t362 - qJ(3) * t358 - qJ(4) * t361 + t485 * qJD(2) + t488 * qJDD(2) + t511 * t432 + t489 * t433 - t486 * t492 + t475) + t465 * (-pkin(2) * t358 + pkin(3) * t360 + pkin(7) * t362 - mrSges(5,2) * t378 + mrSges(5,3) * t379 - mrSges(4,1) * t381 + mrSges(4,2) * t384 + mrSges(3,3) * t399 - mrSges(3,1) * t415 + t461 * t363 - t464 * t364 - qJ(4) * t472 + t510 * t433 + t489 * t432 + t487 * qJDD(2) + t484 * qJD(2) + (qJ(4) * t430 * t465 + t486 * t462) * qJD(1)) + pkin(1) * (-m(3) * t415 + (-mrSges(3,2) + mrSges(4,3)) * t432 + (mrSges(3,1) - t501) * t433 + ((t442 + t443) * t465 + (-t439 - t495) * t462) * qJD(1) - t476) + pkin(6) * (t465 * (m(3) * t399 - qJDD(2) * mrSges(3,2) - qJD(2) * t439 + t469 + (mrSges(3,3) + mrSges(4,2)) * t433) + (t429 - t430) * t458 * qJD(1) + (-m(3) * t398 - qJDD(2) * mrSges(3,1) - qJD(2) * t442 + t470 - (-mrSges(3,3) - t500) * t432 + (t428 + t429) * t493) * t462); -pkin(2) * t359 - pkin(3) * t361 - pkin(7) * t496 - mrSges(5,1) * t379 + mrSges(5,2) * t380 + mrSges(4,3) * t384 - mrSges(4,1) * t385 + mrSges(3,1) * t398 - mrSges(3,2) * t399 - pkin(4) * t477 - t461 * t364 - t464 * t363 + qJ(3) * t469 + (mrSges(4,2) * qJ(3) + t487) * t433 + t488 * t432 + t509 * qJDD(2) + (-t485 * t462 + (-qJ(3) * t430 - t484) * t465) * qJD(1); t359; t360; t475;];
tauJ = t1;
