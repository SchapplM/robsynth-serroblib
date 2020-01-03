% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:46
% EndTime: 2019-12-31 19:27:47
% DurationCPUTime: 1.25s
% Computational Cost: add. (17200->182), mult. (19848->217), div. (0->0), fcn. (8392->8), ass. (0->80)
t457 = sin(qJ(1));
t460 = cos(qJ(1));
t436 = t457 * g(1) - g(2) * t460;
t431 = qJDD(1) * pkin(1) + t436;
t437 = -g(1) * t460 - g(2) * t457;
t461 = qJD(1) ^ 2;
t432 = -pkin(1) * t461 + t437;
t456 = sin(qJ(2));
t459 = cos(qJ(2));
t419 = t456 * t431 + t459 * t432;
t448 = qJDD(1) + qJDD(2);
t449 = (qJD(1) + qJD(2));
t468 = t448 * qJ(3) + (2 * qJD(3) * t449) + t419;
t479 = (-pkin(2) - pkin(3));
t481 = t449 ^ 2;
t410 = (t479 * t481) + t468;
t418 = t459 * t431 - t456 * t432;
t465 = -qJ(3) * t481 + qJDD(3) - t418;
t414 = t448 * t479 + t465;
t453 = sin(pkin(8));
t454 = cos(pkin(8));
t408 = t454 * t410 + t453 * t414;
t405 = -(pkin(4) * t481) - pkin(7) * t448 + t408;
t452 = g(3) + qJDD(4);
t455 = sin(qJ(5));
t458 = cos(qJ(5));
t402 = -t405 * t455 + t452 * t458;
t403 = t405 * t458 + t452 * t455;
t421 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t455 - Ifges(6,2) * t458) * t449;
t422 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t455 - Ifges(6,4) * t458) * t449;
t472 = qJD(5) * t449;
t426 = -t448 * t455 - t458 * t472;
t427 = -t448 * t458 + t455 * t472;
t482 = mrSges(6,1) * t402 - mrSges(6,2) * t403 + Ifges(6,5) * t426 + Ifges(6,6) * t427 + Ifges(6,3) * qJDD(5) - (t421 * t455 - t422 * t458) * t449;
t480 = -m(3) - m(4);
t478 = (mrSges(3,1) + mrSges(4,1));
t477 = Ifges(4,4) + Ifges(3,5);
t476 = (Ifges(3,6) - Ifges(4,6));
t475 = t449 * t455;
t474 = t449 * t458;
t415 = -(pkin(2) * t481) + t468;
t425 = (mrSges(6,1) * t458 - mrSges(6,2) * t455) * t449;
t434 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t474;
t400 = m(6) * t402 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t426 + qJD(5) * t434 + t425 * t475;
t433 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t475;
t401 = m(6) * t403 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t427 - qJD(5) * t433 - t425 * t474;
t394 = -t400 * t455 + t458 * t401;
t390 = m(5) * t408 - (mrSges(5,1) * t481) + mrSges(5,2) * t448 + t394;
t407 = -t410 * t453 + t414 * t454;
t404 = pkin(4) * t448 - pkin(7) * t481 - t407;
t398 = -m(6) * t404 + t427 * mrSges(6,1) - mrSges(6,2) * t426 + t433 * t475 - t434 * t474;
t397 = m(5) * t407 - mrSges(5,1) * t448 - mrSges(5,2) * t481 + t398;
t469 = t454 * t390 - t453 * t397;
t466 = m(4) * t415 + t448 * mrSges(4,3) + t469;
t381 = m(3) * t419 - t448 * mrSges(3,2) - (t478 * t481) + t466;
t388 = t390 * t453 + t397 * t454;
t416 = -pkin(2) * t448 + t465;
t386 = m(4) * t416 - t448 * mrSges(4,1) - mrSges(4,3) * t481 + t388;
t383 = m(3) * t418 + mrSges(3,1) * t448 - mrSges(3,2) * t481 - t386;
t377 = t381 * t456 + t383 * t459;
t374 = m(2) * t436 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t461 + t377;
t470 = t381 * t459 - t383 * t456;
t375 = m(2) * t437 - mrSges(2,1) * t461 - qJDD(1) * mrSges(2,2) + t470;
t473 = t374 * t460 + t375 * t457;
t471 = -t374 * t457 + t375 * t460;
t393 = t458 * t400 + t455 * t401;
t392 = m(5) * t452 + t393;
t420 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t455 - Ifges(6,6) * t458) * t449;
t395 = -mrSges(6,1) * t404 + mrSges(6,3) * t403 + Ifges(6,4) * t426 + Ifges(6,2) * t427 + Ifges(6,6) * qJDD(5) + qJD(5) * t422 + t420 * t475;
t396 = mrSges(6,2) * t404 - mrSges(6,3) * t402 + Ifges(6,1) * t426 + Ifges(6,4) * t427 + Ifges(6,5) * qJDD(5) - qJD(5) * t421 - t420 * t474;
t463 = -mrSges(4,1) * t416 - mrSges(5,1) * t407 - mrSges(3,2) * t419 - pkin(3) * t388 - pkin(4) * t398 - pkin(7) * t394 - t458 * t395 - t455 * t396 + qJ(3) * (-mrSges(4,1) * t481 + t466) - pkin(2) * t386 + mrSges(5,2) * t408 + mrSges(4,3) * t415 + mrSges(3,1) * t418 + (Ifges(5,3) + Ifges(3,3) + Ifges(4,2)) * t448;
t462 = mrSges(2,1) * t436 - mrSges(2,2) * t437 + Ifges(2,3) * qJDD(1) + pkin(1) * t377 + t463;
t391 = -m(4) * g(3) - t392;
t387 = -mrSges(5,1) * t452 + mrSges(5,3) * t408 + (Ifges(5,5) * t481) - Ifges(5,6) * t448 - pkin(4) * t393 - t482;
t378 = mrSges(5,2) * t452 - mrSges(5,3) * t407 - Ifges(5,5) * t448 - Ifges(5,6) * t481 - pkin(7) * t393 - t395 * t455 + t396 * t458;
t370 = mrSges(4,2) * t416 - mrSges(3,3) * t418 - qJ(3) * t391 - qJ(4) * t388 + t454 * t378 - t453 * t387 + t477 * t448 - (t476 * t481) + (-mrSges(3,2) + mrSges(4,3)) * g(3);
t369 = mrSges(4,2) * t415 + mrSges(3,3) * t419 - pkin(2) * t391 + pkin(3) * t392 + g(3) * t478 - qJ(4) * t469 - t453 * t378 - t454 * t387 + t448 * t476 + t477 * t481;
t368 = -mrSges(2,2) * g(3) - mrSges(2,3) * t436 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t461 - pkin(6) * t377 - t369 * t456 + t370 * t459;
t367 = Ifges(2,6) * qJDD(1) + t461 * Ifges(2,5) + mrSges(2,3) * t437 + t456 * t370 + t459 * t369 + pkin(1) * t392 + pkin(6) * t470 + (-pkin(1) * t480 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t471; -m(1) * g(2) + t473; (-m(1) - m(2) + t480) * g(3) - t392; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t473 - t367 * t457 + t368 * t460; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t471 + t460 * t367 + t457 * t368; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t462; t462; t463; t386; t392; t482;];
tauJB = t1;
