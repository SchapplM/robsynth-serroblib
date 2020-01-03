% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:06
% EndTime: 2019-12-31 17:14:07
% DurationCPUTime: 0.86s
% Computational Cost: add. (7611->174), mult. (9590->214), div. (0->0), fcn. (4288->6), ass. (0->74)
t490 = Ifges(4,1) + Ifges(5,1);
t482 = Ifges(4,4) - Ifges(5,5);
t481 = Ifges(4,5) + Ifges(5,4);
t489 = Ifges(4,2) + Ifges(5,3);
t480 = Ifges(4,6) - Ifges(5,6);
t488 = Ifges(4,3) + Ifges(5,2);
t451 = qJD(1) + qJD(2);
t455 = sin(qJ(3));
t458 = cos(qJ(3));
t428 = (-mrSges(5,1) * t458 - mrSges(5,3) * t455) * t451;
t450 = qJDD(1) + qJDD(2);
t472 = qJD(3) * t451;
t430 = t455 * t450 + t458 * t472;
t457 = sin(qJ(1));
t460 = cos(qJ(1));
t444 = t457 * g(1) - t460 * g(2);
t437 = qJDD(1) * pkin(1) + t444;
t445 = -t460 * g(1) - t457 * g(2);
t462 = qJD(1) ^ 2;
t438 = -t462 * pkin(1) + t445;
t456 = sin(qJ(2));
t459 = cos(qJ(2));
t414 = t456 * t437 + t459 * t438;
t449 = t451 ^ 2;
t411 = -t449 * pkin(2) + t450 * pkin(6) + t414;
t427 = (-pkin(3) * t458 - qJ(4) * t455) * t451;
t461 = qJD(3) ^ 2;
t484 = t458 * g(3);
t406 = -qJDD(3) * pkin(3) + t484 - t461 * qJ(4) + qJDD(4) + (t427 * t451 + t411) * t455;
t477 = t451 * t458;
t442 = mrSges(5,2) * t477 + qJD(3) * mrSges(5,3);
t467 = -m(5) * t406 + qJDD(3) * mrSges(5,1) + qJD(3) * t442;
t478 = t451 * t455;
t402 = t430 * mrSges(5,2) + t428 * t478 - t467;
t408 = -t455 * g(3) + t458 * t411;
t405 = -t461 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t427 * t477 + t408;
t407 = -t455 * t411 - t484;
t431 = t458 * t450 - t455 * t472;
t440 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t478;
t468 = m(5) * t405 + qJDD(3) * mrSges(5,3) + qJD(3) * t440 + t428 * t477;
t473 = (t490 * t455 + t482 * t458) * t451 + t481 * qJD(3);
t475 = (-t482 * t455 - t489 * t458) * t451 - t480 * qJD(3);
t487 = -(t475 * t455 + t473 * t458) * t451 + t488 * qJDD(3) + t481 * t430 + t480 * t431 + mrSges(4,1) * t407 - mrSges(5,1) * t406 - mrSges(4,2) * t408 + mrSges(5,3) * t405 - pkin(3) * t402 + qJ(4) * (t431 * mrSges(5,2) + t468);
t483 = mrSges(4,3) + mrSges(5,2);
t429 = (-mrSges(4,1) * t458 + mrSges(4,2) * t455) * t451;
t439 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t478;
t398 = m(4) * t408 - qJDD(3) * mrSges(4,2) - qJD(3) * t439 + t429 * t477 + t483 * t431 + t468;
t441 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t477;
t399 = m(4) * t407 + qJDD(3) * mrSges(4,1) + qJD(3) * t441 + (-t428 - t429) * t478 - t483 * t430 + t467;
t469 = t458 * t398 - t455 * t399;
t389 = m(3) * t414 - t449 * mrSges(3,1) - t450 * mrSges(3,2) + t469;
t413 = t459 * t437 - t456 * t438;
t410 = -t450 * pkin(2) - t449 * pkin(6) - t413;
t403 = -t431 * pkin(3) - t430 * qJ(4) + (-0.2e1 * qJD(4) * t455 + (pkin(3) * t455 - qJ(4) * t458) * qJD(3)) * t451 + t410;
t400 = m(5) * t403 - t431 * mrSges(5,1) - t430 * mrSges(5,3) - t440 * t478 - t442 * t477;
t463 = -m(4) * t410 + t431 * mrSges(4,1) - t430 * mrSges(4,2) - t439 * t478 + t441 * t477 - t400;
t393 = m(3) * t413 + t450 * mrSges(3,1) - t449 * mrSges(3,2) + t463;
t382 = t456 * t389 + t459 * t393;
t379 = m(2) * t444 + qJDD(1) * mrSges(2,1) - t462 * mrSges(2,2) + t382;
t470 = t459 * t389 - t456 * t393;
t380 = m(2) * t445 - t462 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t470;
t476 = t460 * t379 + t457 * t380;
t391 = t455 * t398 + t458 * t399;
t474 = (t481 * t455 + t480 * t458) * t451 + t488 * qJD(3);
t471 = -t457 * t379 + t460 * t380;
t385 = -mrSges(4,1) * t410 - mrSges(5,1) * t403 + mrSges(5,2) * t405 + mrSges(4,3) * t408 - pkin(3) * t400 + t473 * qJD(3) + t480 * qJDD(3) + t482 * t430 + t489 * t431 - t474 * t478;
t386 = mrSges(4,2) * t410 + mrSges(5,2) * t406 - mrSges(4,3) * t407 - mrSges(5,3) * t403 - qJ(4) * t400 + t475 * qJD(3) + t481 * qJDD(3) + t490 * t430 + t482 * t431 + t474 * t477;
t466 = mrSges(3,1) * t413 - mrSges(3,2) * t414 + Ifges(3,3) * t450 + pkin(2) * t463 + pkin(6) * t469 + t458 * t385 + t455 * t386;
t465 = mrSges(2,1) * t444 - mrSges(2,2) * t445 + Ifges(2,3) * qJDD(1) + pkin(1) * t382 + t466;
t375 = mrSges(3,1) * g(3) + mrSges(3,3) * t414 + t449 * Ifges(3,5) + Ifges(3,6) * t450 - pkin(2) * t391 - t487;
t374 = -mrSges(3,2) * g(3) - mrSges(3,3) * t413 + Ifges(3,5) * t450 - t449 * Ifges(3,6) - pkin(6) * t391 - t455 * t385 + t458 * t386;
t373 = -mrSges(2,2) * g(3) - mrSges(2,3) * t444 + Ifges(2,5) * qJDD(1) - t462 * Ifges(2,6) - pkin(5) * t382 + t459 * t374 - t456 * t375;
t372 = Ifges(2,6) * qJDD(1) + t462 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t445 + t456 * t374 + t459 * t375 - pkin(1) * (-m(3) * g(3) + t391) + pkin(5) * t470;
t1 = [-m(1) * g(1) + t471; -m(1) * g(2) + t476; (-m(1) - m(2) - m(3)) * g(3) + t391; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t476 - t457 * t372 + t460 * t373; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t471 + t460 * t372 + t457 * t373; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t465; t465; t466; t487; t402;];
tauJB = t1;
