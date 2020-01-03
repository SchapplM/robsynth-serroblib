% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRPR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:52
% EndTime: 2019-12-31 16:20:53
% DurationCPUTime: 0.98s
% Computational Cost: add. (8490->163), mult. (18343->213), div. (0->0), fcn. (11630->8), ass. (0->81)
t449 = qJD(2) ^ 2;
t443 = cos(pkin(7));
t474 = pkin(3) * t443;
t441 = sin(pkin(7));
t473 = mrSges(4,2) * t441;
t438 = t443 ^ 2;
t472 = t438 * t449;
t442 = sin(pkin(6));
t444 = cos(pkin(6));
t426 = t442 * g(1) - t444 * g(2);
t427 = -t444 * g(1) - t442 * g(2);
t446 = sin(qJ(2));
t448 = cos(qJ(2));
t416 = t446 * t426 + t448 * t427;
t413 = -t449 * pkin(2) + qJDD(2) * qJ(3) + t416;
t440 = -g(3) + qJDD(1);
t468 = qJD(2) * qJD(3);
t470 = t443 * t440 - 0.2e1 * t441 * t468;
t398 = (-pkin(5) * qJDD(2) + t449 * t474 - t413) * t441 + t470;
t402 = t441 * t440 + (t413 + 0.2e1 * t468) * t443;
t467 = qJDD(2) * t443;
t399 = -pkin(3) * t472 + pkin(5) * t467 + t402;
t445 = sin(qJ(4));
t447 = cos(qJ(4));
t396 = t447 * t398 - t445 * t399;
t455 = -t441 * t445 + t443 * t447;
t419 = t455 * qJD(2);
t456 = t441 * t447 + t443 * t445;
t420 = t456 * qJD(2);
t408 = -t419 * mrSges(5,1) + t420 * mrSges(5,2);
t412 = t419 * qJD(4) + t456 * qJDD(2);
t417 = -qJD(4) * mrSges(5,2) + t419 * mrSges(5,3);
t394 = m(5) * t396 + qJDD(4) * mrSges(5,1) - t412 * mrSges(5,3) + qJD(4) * t417 - t420 * t408;
t397 = t445 * t398 + t447 * t399;
t411 = -t420 * qJD(4) + t455 * qJDD(2);
t418 = qJD(4) * mrSges(5,1) - t420 * mrSges(5,3);
t395 = m(5) * t397 - qJDD(4) * mrSges(5,2) + t411 * mrSges(5,3) - qJD(4) * t418 + t419 * t408;
t384 = t447 * t394 + t445 * t395;
t401 = -t441 * t413 + t470;
t454 = mrSges(4,3) * qJDD(2) + t449 * (-mrSges(4,1) * t443 + t473);
t382 = m(4) * t401 - t454 * t441 + t384;
t462 = -t445 * t394 + t447 * t395;
t383 = m(4) * t402 + t454 * t443 + t462;
t463 = -t441 * t382 + t443 * t383;
t376 = m(3) * t416 - t449 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t463;
t415 = t448 * t426 - t446 * t427;
t457 = qJDD(3) - t415;
t410 = -qJDD(2) * pkin(2) - t449 * qJ(3) + t457;
t437 = t441 ^ 2;
t400 = (-pkin(2) - t474) * qJDD(2) + (-qJ(3) + (-t437 - t438) * pkin(5)) * t449 + t457;
t452 = m(5) * t400 - t411 * mrSges(5,1) + t412 * mrSges(5,2) - t419 * t417 + t420 * t418;
t451 = -m(4) * t410 + mrSges(4,1) * t467 - t452 + (t437 * t449 + t472) * mrSges(4,3);
t388 = m(3) * t415 - t449 * mrSges(3,2) + (mrSges(3,1) - t473) * qJDD(2) + t451;
t371 = t446 * t376 + t448 * t388;
t369 = m(2) * t426 + t371;
t464 = t448 * t376 - t446 * t388;
t370 = m(2) * t427 + t464;
t471 = t444 * t369 + t442 * t370;
t378 = t443 * t382 + t441 * t383;
t458 = Ifges(4,5) * t441 + Ifges(4,6) * t443;
t469 = t449 * t458;
t466 = m(3) * t440 + t378;
t465 = -t442 * t369 + t444 * t370;
t461 = m(2) * t440 + t466;
t460 = Ifges(4,1) * t441 + Ifges(4,4) * t443;
t459 = Ifges(4,4) * t441 + Ifges(4,2) * t443;
t403 = Ifges(5,5) * t420 + Ifges(5,6) * t419 + Ifges(5,3) * qJD(4);
t405 = Ifges(5,1) * t420 + Ifges(5,4) * t419 + Ifges(5,5) * qJD(4);
t385 = -mrSges(5,1) * t400 + mrSges(5,3) * t397 + Ifges(5,4) * t412 + Ifges(5,2) * t411 + Ifges(5,6) * qJDD(4) + qJD(4) * t405 - t420 * t403;
t404 = Ifges(5,4) * t420 + Ifges(5,2) * t419 + Ifges(5,6) * qJD(4);
t386 = mrSges(5,2) * t400 - mrSges(5,3) * t396 + Ifges(5,1) * t412 + Ifges(5,4) * t411 + Ifges(5,5) * qJDD(4) - qJD(4) * t404 + t419 * t403;
t365 = -mrSges(4,1) * t410 + mrSges(4,3) * t402 - pkin(3) * t452 + pkin(5) * t462 + t459 * qJDD(2) + t447 * t385 + t445 * t386 - t441 * t469;
t373 = mrSges(4,2) * t410 - mrSges(4,3) * t401 - pkin(5) * t384 + t460 * qJDD(2) - t445 * t385 + t447 * t386 + t443 * t469;
t390 = qJDD(2) * t473 - t451;
t453 = mrSges(3,1) * t415 - mrSges(3,2) * t416 + Ifges(3,3) * qJDD(2) - pkin(2) * t390 + qJ(3) * t463 + t443 * t365 + t441 * t373;
t450 = mrSges(5,1) * t396 - mrSges(5,2) * t397 + Ifges(5,5) * t412 + Ifges(5,6) * t411 + Ifges(5,3) * qJDD(4) + t420 * t404 - t419 * t405;
t363 = -mrSges(3,1) * t440 - mrSges(4,1) * t401 + mrSges(4,2) * t402 + mrSges(3,3) * t416 - pkin(2) * t378 - pkin(3) * t384 + (Ifges(3,6) - t458) * qJDD(2) - t450 + (-t441 * t459 + t443 * t460 + Ifges(3,5)) * t449;
t362 = mrSges(3,2) * t440 - mrSges(3,3) * t415 + Ifges(3,5) * qJDD(2) - t449 * Ifges(3,6) - qJ(3) * t378 - t441 * t365 + t443 * t373;
t361 = mrSges(2,2) * t440 - mrSges(2,3) * t426 - pkin(4) * t371 + t448 * t362 - t446 * t363;
t360 = -mrSges(2,1) * t440 + mrSges(2,3) * t427 - pkin(1) * t466 + pkin(4) * t464 + t446 * t362 + t448 * t363;
t1 = [-m(1) * g(1) + t465; -m(1) * g(2) + t471; -m(1) * g(3) + t461; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t471 - t442 * t360 + t444 * t361; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t465 + t444 * t360 + t442 * t361; -mrSges(1,1) * g(2) + mrSges(2,1) * t426 + mrSges(1,2) * g(1) - mrSges(2,2) * t427 + pkin(1) * t371 + t453; t461; t453; t390; t450;];
tauJB = t1;
