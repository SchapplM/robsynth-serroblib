% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:49
% EndTime: 2019-12-31 18:08:51
% DurationCPUTime: 2.18s
% Computational Cost: add. (20934->246), mult. (43897->301), div. (0->0), fcn. (25110->8), ass. (0->98)
t526 = Ifges(5,1) + Ifges(6,1);
t521 = Ifges(5,4) - Ifges(6,5);
t520 = Ifges(5,5) + Ifges(6,4);
t525 = Ifges(5,2) + Ifges(6,3);
t524 = -Ifges(6,2) - Ifges(5,3);
t519 = Ifges(5,6) - Ifges(6,6);
t523 = -2 * qJD(4);
t522 = -mrSges(5,3) - mrSges(6,2);
t518 = cos(pkin(8));
t493 = sin(qJ(1));
t495 = cos(qJ(1));
t478 = t493 * g(1) - t495 * g(2);
t470 = qJDD(1) * pkin(1) + t478;
t479 = -t495 * g(1) - t493 * g(2);
t497 = qJD(1) ^ 2;
t472 = -t497 * pkin(1) + t479;
t490 = sin(pkin(7));
t491 = cos(pkin(7));
t448 = t490 * t470 + t491 * t472;
t437 = -t497 * pkin(2) + qJDD(1) * pkin(6) + t448;
t488 = -g(3) + qJDD(2);
t492 = sin(qJ(3));
t494 = cos(qJ(3));
t427 = -t492 * t437 + t494 * t488;
t510 = qJD(1) * qJD(3);
t507 = t494 * t510;
t473 = t492 * qJDD(1) + t507;
t424 = (-t473 + t507) * qJ(4) + (t492 * t494 * t497 + qJDD(3)) * pkin(3) + t427;
t428 = t494 * t437 + t492 * t488;
t474 = t494 * qJDD(1) - t492 * t510;
t512 = qJD(1) * t492;
t475 = qJD(3) * pkin(3) - qJ(4) * t512;
t487 = t494 ^ 2;
t425 = -t487 * t497 * pkin(3) + t474 * qJ(4) - qJD(3) * t475 + t428;
t489 = sin(pkin(8));
t511 = qJD(1) * t494;
t458 = t489 * t512 - t518 * t511;
t421 = t489 * t424 + t518 * t425 + t458 * t523;
t449 = t489 * t473 - t518 * t474;
t459 = (t489 * t494 + t518 * t492) * qJD(1);
t454 = qJD(3) * mrSges(5,1) - t459 * mrSges(5,3);
t441 = t458 * pkin(4) - t459 * qJ(5);
t496 = qJD(3) ^ 2;
t416 = -t496 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t458 * t441 + t421;
t455 = -qJD(3) * mrSges(6,1) + t459 * mrSges(6,2);
t509 = m(6) * t416 + qJDD(3) * mrSges(6,3) + qJD(3) * t455;
t442 = t458 * mrSges(6,1) - t459 * mrSges(6,3);
t513 = -t458 * mrSges(5,1) - t459 * mrSges(5,2) - t442;
t412 = m(5) * t421 - qJDD(3) * mrSges(5,2) - qJD(3) * t454 + t522 * t449 + t513 * t458 + t509;
t500 = t518 * t424 - t489 * t425;
t420 = t459 * t523 + t500;
t450 = t518 * t473 + t489 * t474;
t453 = -qJD(3) * mrSges(5,2) - t458 * mrSges(5,3);
t417 = -qJDD(3) * pkin(4) - t496 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t441) * t459 - t500;
t456 = -t458 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t502 = -m(6) * t417 + qJDD(3) * mrSges(6,1) + qJD(3) * t456;
t413 = m(5) * t420 + qJDD(3) * mrSges(5,1) + qJD(3) * t453 + t522 * t450 + t513 * t459 + t502;
t406 = t489 * t412 + t518 * t413;
t471 = (-mrSges(4,1) * t494 + mrSges(4,2) * t492) * qJD(1);
t477 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t511;
t402 = m(4) * t427 + qJDD(3) * mrSges(4,1) - t473 * mrSges(4,3) + qJD(3) * t477 - t471 * t512 + t406;
t476 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t512;
t503 = t518 * t412 - t489 * t413;
t403 = m(4) * t428 - qJDD(3) * mrSges(4,2) + t474 * mrSges(4,3) - qJD(3) * t476 + t471 * t511 + t503;
t504 = -t492 * t402 + t494 * t403;
t397 = m(3) * t448 - t497 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t504;
t447 = t491 * t470 - t490 * t472;
t501 = -qJDD(1) * pkin(2) - t447;
t436 = -t497 * pkin(6) + t501;
t426 = -t474 * pkin(3) + qJDD(4) + t475 * t512 + (-qJ(4) * t487 - pkin(6)) * t497 + t501;
t419 = -0.2e1 * qJD(5) * t459 + (qJD(3) * t458 - t450) * qJ(5) + (qJD(3) * t459 + t449) * pkin(4) + t426;
t414 = m(6) * t419 + t449 * mrSges(6,1) - t450 * mrSges(6,3) - t459 * t455 + t458 * t456;
t499 = m(5) * t426 + t449 * mrSges(5,1) + t450 * mrSges(5,2) + t458 * t453 + t459 * t454 + t414;
t498 = -m(4) * t436 + t474 * mrSges(4,1) - t473 * mrSges(4,2) - t476 * t512 + t477 * t511 - t499;
t408 = m(3) * t447 + qJDD(1) * mrSges(3,1) - t497 * mrSges(3,2) + t498;
t394 = t490 * t397 + t491 * t408;
t392 = m(2) * t478 + qJDD(1) * mrSges(2,1) - t497 * mrSges(2,2) + t394;
t505 = t491 * t397 - t490 * t408;
t393 = m(2) * t479 - t497 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t505;
t517 = t495 * t392 + t493 * t393;
t398 = t494 * t402 + t492 * t403;
t516 = -t519 * qJD(3) + t525 * t458 - t521 * t459;
t515 = t524 * qJD(3) + t519 * t458 - t520 * t459;
t514 = t520 * qJD(3) - t521 * t458 + t526 * t459;
t508 = m(3) * t488 + t398;
t506 = -t493 * t392 + t495 * t393;
t465 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t492 + Ifges(4,4) * t494) * qJD(1);
t464 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t492 + Ifges(4,2) * t494) * qJD(1);
t463 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t492 + Ifges(4,6) * t494) * qJD(1);
t405 = mrSges(5,2) * t426 + mrSges(6,2) * t417 - mrSges(5,3) * t420 - mrSges(6,3) * t419 - qJ(5) * t414 + t516 * qJD(3) + t520 * qJDD(3) - t521 * t449 + t526 * t450 + t515 * t458;
t404 = -mrSges(5,1) * t426 - mrSges(6,1) * t419 + mrSges(6,2) * t416 + mrSges(5,3) * t421 - pkin(4) * t414 + t514 * qJD(3) + t519 * qJDD(3) - t525 * t449 + t521 * t450 + t515 * t459;
t388 = mrSges(4,2) * t436 - mrSges(4,3) * t427 + Ifges(4,1) * t473 + Ifges(4,4) * t474 + Ifges(4,5) * qJDD(3) - qJ(4) * t406 - qJD(3) * t464 - t489 * t404 + t518 * t405 + t463 * t511;
t387 = -mrSges(4,1) * t436 + mrSges(4,3) * t428 + Ifges(4,4) * t473 + Ifges(4,2) * t474 + Ifges(4,6) * qJDD(3) - pkin(3) * t499 + qJ(4) * t503 + qJD(3) * t465 + t518 * t404 + t489 * t405 - t463 * t512;
t386 = Ifges(3,6) * qJDD(1) + t497 * Ifges(3,5) - qJ(5) * t509 - pkin(4) * t502 - mrSges(3,1) * t488 - Ifges(4,5) * t473 - Ifges(4,6) * t474 + mrSges(3,3) * t448 - mrSges(5,1) * t420 + mrSges(5,2) * t421 - mrSges(4,1) * t427 + mrSges(4,2) * t428 - mrSges(6,3) * t416 + mrSges(6,1) * t417 - pkin(3) * t406 - pkin(2) * t398 + (pkin(4) * t442 + t516) * t459 + (qJ(5) * t442 - t514) * t458 + (pkin(4) * mrSges(6,2) - t520) * t450 + (qJ(5) * mrSges(6,2) + t519) * t449 + (-t492 * t464 + t494 * t465) * qJD(1) + (-Ifges(4,3) + t524) * qJDD(3);
t385 = mrSges(3,2) * t488 - mrSges(3,3) * t447 + Ifges(3,5) * qJDD(1) - t497 * Ifges(3,6) - pkin(6) * t398 - t492 * t387 + t494 * t388;
t384 = -mrSges(2,2) * g(3) - mrSges(2,3) * t478 + Ifges(2,5) * qJDD(1) - t497 * Ifges(2,6) - qJ(2) * t394 + t491 * t385 - t490 * t386;
t383 = mrSges(2,1) * g(3) + mrSges(2,3) * t479 + t497 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t508 + qJ(2) * t505 + t490 * t385 + t491 * t386;
t1 = [-m(1) * g(1) + t506; -m(1) * g(2) + t517; (-m(1) - m(2)) * g(3) + t508; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t517 - t493 * t383 + t495 * t384; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t506 + t495 * t383 + t493 * t384; pkin(1) * t394 + mrSges(2,1) * t478 - mrSges(2,2) * t479 + t492 * t388 + t494 * t387 + pkin(2) * t498 + pkin(6) * t504 + mrSges(3,1) * t447 - mrSges(3,2) * t448 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
