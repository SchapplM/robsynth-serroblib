% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:01
% EndTime: 2019-12-31 18:26:02
% DurationCPUTime: 1.28s
% Computational Cost: add. (17627->182), mult. (23191->217), div. (0->0), fcn. (8320->8), ass. (0->80)
t456 = qJD(1) ^ 2;
t452 = sin(qJ(1));
t455 = cos(qJ(1));
t435 = -t455 * g(1) - t452 * g(2);
t461 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t435;
t475 = -pkin(1) - pkin(2);
t417 = t456 * t475 + t461;
t434 = t452 * g(1) - t455 * g(2);
t460 = -t456 * qJ(2) + qJDD(2) - t434;
t420 = qJDD(1) * t475 + t460;
t451 = sin(qJ(3));
t454 = cos(qJ(3));
t412 = -t417 * t451 + t454 * t420;
t439 = -qJDD(1) + qJDD(3);
t409 = pkin(3) * t439 + t412;
t413 = t454 * t417 + t451 * t420;
t440 = -qJD(1) + qJD(3);
t438 = t440 ^ 2;
t410 = -pkin(3) * t438 + t413;
t448 = sin(pkin(8));
t449 = cos(pkin(8));
t406 = t448 * t409 + t449 * t410;
t403 = -(pkin(4) * t438) + pkin(7) * t439 + t406;
t446 = g(3) + qJDD(4);
t450 = sin(qJ(5));
t453 = cos(qJ(5));
t400 = -t403 * t450 + t446 * t453;
t428 = (-mrSges(6,1) * t453 + mrSges(6,2) * t450) * t440;
t467 = qJD(5) * t440;
t429 = t439 * t450 + t453 * t467;
t469 = t440 * t453;
t432 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t469;
t470 = t440 * t450;
t398 = m(6) * t400 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t429 + qJD(5) * t432 - t428 * t470;
t401 = t403 * t453 + t446 * t450;
t430 = t439 * t453 - t450 * t467;
t431 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t470;
t399 = m(6) * t401 - qJDD(5) * mrSges(6,2) + t430 * mrSges(6,3) - qJD(5) * t431 + t428 * t469;
t389 = -t398 * t450 + t453 * t399;
t385 = m(5) * t406 - (mrSges(5,1) * t438) - mrSges(5,2) * t439 + t389;
t405 = t409 * t449 - t410 * t448;
t402 = -pkin(4) * t439 - pkin(7) * t438 - t405;
t394 = -m(6) * t402 + t430 * mrSges(6,1) - mrSges(6,2) * t429 - t431 * t470 + t432 * t469;
t393 = m(5) * t405 + mrSges(5,1) * t439 - mrSges(5,2) * t438 + t394;
t382 = t385 * t448 + t393 * t449;
t421 = (Ifges(6,3) * qJD(5)) + (Ifges(6,5) * t450 + Ifges(6,6) * t453) * t440;
t423 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t450 + Ifges(6,4) * t453) * t440;
t390 = -mrSges(6,1) * t402 + mrSges(6,3) * t401 + Ifges(6,4) * t429 + Ifges(6,2) * t430 + Ifges(6,6) * qJDD(5) + qJD(5) * t423 - t421 * t470;
t422 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t450 + Ifges(6,2) * t453) * t440;
t391 = mrSges(6,2) * t402 - mrSges(6,3) * t400 + Ifges(6,1) * t429 + Ifges(6,4) * t430 + Ifges(6,5) * qJDD(5) - qJD(5) * t422 + t421 * t469;
t478 = mrSges(4,1) * t412 + mrSges(5,1) * t405 - mrSges(4,2) * t413 - mrSges(5,2) * t406 + pkin(3) * t382 + pkin(4) * t394 + pkin(7) * t389 + t453 * t390 + t450 * t391 + (Ifges(4,3) + Ifges(5,3)) * t439;
t477 = mrSges(6,1) * t400 - mrSges(6,2) * t401 + Ifges(6,5) * t429 + Ifges(6,6) * t430 + Ifges(6,3) * qJDD(5) + (t422 * t450 - t423 * t453) * t440;
t476 = -m(3) - m(4);
t474 = -mrSges(2,1) - mrSges(3,1);
t473 = Ifges(3,4) + Ifges(2,5);
t472 = Ifges(2,6) - Ifges(3,6);
t424 = -pkin(1) * t456 + t461;
t379 = m(4) * t412 + mrSges(4,1) * t439 - mrSges(4,2) * t438 + t382;
t464 = t385 * t449 - t393 * t448;
t380 = m(4) * t413 - mrSges(4,1) * t438 - mrSges(4,2) * t439 + t464;
t465 = -t451 * t379 + t380 * t454;
t462 = m(3) * t424 + qJDD(1) * mrSges(3,3) + t465;
t371 = m(2) * t435 - qJDD(1) * mrSges(2,2) + t456 * t474 + t462;
t376 = t379 * t454 + t380 * t451;
t427 = -qJDD(1) * pkin(1) + t460;
t375 = m(3) * t427 - qJDD(1) * mrSges(3,1) - t456 * mrSges(3,3) + t376;
t372 = m(2) * t434 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t456 - t375;
t468 = t371 * t452 + t372 * t455;
t388 = t453 * t398 + t450 * t399;
t466 = t371 * t455 - t372 * t452;
t387 = m(5) * t446 + t388;
t457 = -mrSges(3,1) * t427 - mrSges(2,2) * t435 - pkin(2) * t376 + qJ(2) * (-mrSges(3,1) * t456 + t462) - pkin(1) * t375 + mrSges(3,3) * t424 + mrSges(2,1) * t434 - t478 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);
t386 = g(3) * t476 - t387;
t381 = -mrSges(5,1) * t446 + mrSges(5,3) * t406 + (t438 * Ifges(5,5)) + Ifges(5,6) * t439 - pkin(4) * t388 - t477;
t377 = mrSges(5,2) * t446 - mrSges(5,3) * t405 + Ifges(5,5) * t439 - Ifges(5,6) * t438 - pkin(7) * t388 - t390 * t450 + t391 * t453;
t367 = mrSges(4,2) * g(3) - mrSges(4,3) * t412 + Ifges(4,5) * t439 - (Ifges(4,6) * t438) - qJ(4) * t382 + t377 * t449 - t381 * t448;
t366 = -mrSges(4,1) * g(3) + mrSges(4,3) * t413 + t438 * Ifges(4,5) + Ifges(4,6) * t439 - pkin(3) * t387 + qJ(4) * t464 + t448 * t377 + t449 * t381;
t365 = mrSges(3,2) * t427 - mrSges(2,3) * t434 - pkin(6) * t376 - qJ(2) * t386 - t451 * t366 + t454 * t367 - t472 * t456 + t473 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t364 = mrSges(2,3) * t435 + mrSges(3,2) * t424 - t451 * t367 - t454 * t366 + pkin(2) * t387 - pkin(6) * t465 - pkin(1) * t386 + t473 * t456 + t472 * qJDD(1) + (m(4) * pkin(2) - t474) * g(3);
t1 = [-m(1) * g(1) + t466; -m(1) * g(2) + t468; (-m(1) - m(2) + t476) * g(3) - t387; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t468 - t364 * t452 + t365 * t455; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t466 + t455 * t364 + t452 * t365; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t457; t457; t375; t478; t387; t477;];
tauJB = t1;
