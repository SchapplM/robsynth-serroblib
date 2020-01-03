% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRR9
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:13
% EndTime: 2019-12-31 16:56:15
% DurationCPUTime: 0.81s
% Computational Cost: add. (5999->191), mult. (11183->232), div. (0->0), fcn. (5689->6), ass. (0->82)
t452 = sin(qJ(1));
t455 = cos(qJ(1));
t439 = -t455 * g(1) - t452 * g(2);
t482 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t439;
t481 = -pkin(1) - pkin(5);
t451 = sin(qJ(3));
t480 = t451 * g(3);
t479 = mrSges(2,1) - mrSges(3,2);
t478 = -Ifges(3,4) + Ifges(2,5);
t477 = (Ifges(3,5) - Ifges(2,6));
t438 = t452 * g(1) - t455 * g(2);
t457 = qJD(1) ^ 2;
t465 = -t457 * qJ(2) + qJDD(2) - t438;
t414 = t481 * qJDD(1) + t465;
t454 = cos(qJ(3));
t408 = -t454 * g(3) + t451 * t414;
t431 = (mrSges(4,1) * t451 + mrSges(4,2) * t454) * qJD(1);
t473 = qJD(1) * qJD(3);
t470 = t454 * t473;
t433 = -t451 * qJDD(1) - t470;
t475 = qJD(1) * t454;
t437 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t475;
t413 = t481 * t457 - t482;
t471 = t451 * t473;
t434 = t454 * qJDD(1) - t471;
t395 = (-t434 + t471) * pkin(6) + (-t433 + t470) * pkin(3) + t413;
t432 = (pkin(3) * t451 - pkin(6) * t454) * qJD(1);
t456 = qJD(3) ^ 2;
t474 = t451 * qJD(1);
t397 = -t456 * pkin(3) + qJDD(3) * pkin(6) - t432 * t474 + t408;
t450 = sin(qJ(4));
t453 = cos(qJ(4));
t393 = t453 * t395 - t450 * t397;
t429 = t453 * qJD(3) - t450 * t475;
t404 = t429 * qJD(4) + t450 * qJDD(3) + t453 * t434;
t430 = t450 * qJD(3) + t453 * t475;
t406 = -t429 * mrSges(5,1) + t430 * mrSges(5,2);
t440 = qJD(4) + t474;
t409 = -t440 * mrSges(5,2) + t429 * mrSges(5,3);
t428 = qJDD(4) - t433;
t390 = m(5) * t393 + t428 * mrSges(5,1) - t404 * mrSges(5,3) - t430 * t406 + t440 * t409;
t394 = t450 * t395 + t453 * t397;
t403 = -t430 * qJD(4) + t453 * qJDD(3) - t450 * t434;
t410 = t440 * mrSges(5,1) - t430 * mrSges(5,3);
t391 = m(5) * t394 - t428 * mrSges(5,2) + t403 * mrSges(5,3) + t429 * t406 - t440 * t410;
t467 = -t450 * t390 + t453 * t391;
t379 = m(4) * t408 - qJDD(3) * mrSges(4,2) + t433 * mrSges(4,3) - qJD(3) * t437 - t431 * t474 + t467;
t407 = t454 * t414 + t480;
t436 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t474;
t396 = -qJDD(3) * pkin(3) - t456 * pkin(6) - t480 + (qJD(1) * t432 - t414) * t454;
t462 = -m(5) * t396 + t403 * mrSges(5,1) - t404 * mrSges(5,2) + t429 * t409 - t430 * t410;
t386 = m(4) * t407 + qJDD(3) * mrSges(4,1) - t434 * mrSges(4,3) + qJD(3) * t436 - t431 * t475 + t462;
t373 = t451 * t379 + t454 * t386;
t419 = -qJDD(1) * pkin(1) + t465;
t464 = -m(3) * t419 + (t457 * mrSges(3,3)) - t373;
t369 = m(2) * t438 - (t457 * mrSges(2,2)) + t479 * qJDD(1) + t464;
t417 = t457 * pkin(1) + t482;
t381 = t453 * t390 + t450 * t391;
t463 = -m(4) * t413 + t433 * mrSges(4,1) - t434 * mrSges(4,2) - t436 * t474 - t437 * t475 - t381;
t460 = -m(3) * t417 + (t457 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t463;
t376 = m(2) * t439 - (t457 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t460;
t476 = t455 * t369 + t452 * t376;
t469 = -t452 * t369 + t455 * t376;
t468 = t454 * t379 - t451 * t386;
t398 = Ifges(5,5) * t430 + Ifges(5,6) * t429 + Ifges(5,3) * t440;
t400 = Ifges(5,1) * t430 + Ifges(5,4) * t429 + Ifges(5,5) * t440;
t384 = -mrSges(5,1) * t396 + mrSges(5,3) * t394 + Ifges(5,4) * t404 + Ifges(5,2) * t403 + Ifges(5,6) * t428 - t430 * t398 + t440 * t400;
t399 = Ifges(5,4) * t430 + Ifges(5,2) * t429 + Ifges(5,6) * t440;
t385 = mrSges(5,2) * t396 - mrSges(5,3) * t393 + Ifges(5,1) * t404 + Ifges(5,4) * t403 + Ifges(5,5) * t428 + t429 * t398 - t440 * t399;
t421 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t454 - Ifges(4,2) * t451) * qJD(1);
t422 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t454 - Ifges(4,4) * t451) * qJD(1);
t461 = mrSges(4,1) * t407 - mrSges(4,2) * t408 + Ifges(4,5) * t434 + Ifges(4,6) * t433 + Ifges(4,3) * qJDD(3) + pkin(3) * t462 + pkin(6) * t467 + t453 * t384 + t450 * t385 + t421 * t475 + t422 * t474;
t420 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t454 - Ifges(4,6) * t451) * qJD(1);
t366 = mrSges(4,2) * t413 - mrSges(4,3) * t407 + Ifges(4,1) * t434 + Ifges(4,4) * t433 + Ifges(4,5) * qJDD(3) - pkin(6) * t381 - qJD(3) * t421 - t450 * t384 + t453 * t385 - t420 * t474;
t458 = mrSges(5,1) * t393 - mrSges(5,2) * t394 + Ifges(5,5) * t404 + Ifges(5,6) * t403 + Ifges(5,3) * t428 + t430 * t399 - t429 * t400;
t367 = -mrSges(4,1) * t413 + mrSges(4,3) * t408 + Ifges(4,4) * t434 + Ifges(4,2) * t433 + Ifges(4,6) * qJDD(3) - pkin(3) * t381 + qJD(3) * t422 - t420 * t475 - t458;
t371 = qJDD(1) * mrSges(3,2) - t464;
t459 = mrSges(2,1) * t438 - mrSges(2,2) * t439 + mrSges(3,2) * t419 - mrSges(3,3) * t417 - pkin(1) * t371 - pkin(5) * t373 + qJ(2) * t460 + t454 * t366 - t451 * t367 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t372 = -m(3) * g(3) + t468;
t364 = t461 + (t477 * t457) - mrSges(2,3) * t438 + mrSges(3,1) * t419 + pkin(2) * t373 - qJ(2) * t372 + t478 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t363 = -mrSges(3,1) * t417 + mrSges(2,3) * t439 - pkin(1) * t372 - pkin(2) * t463 - pkin(5) * t468 + t479 * g(3) - t477 * qJDD(1) - t451 * t366 - t454 * t367 + t478 * t457;
t1 = [-m(1) * g(1) + t469; -m(1) * g(2) + t476; (-m(1) - m(2) - m(3)) * g(3) + t468; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t476 - t452 * t363 + t455 * t364; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t469 + t455 * t363 + t452 * t364; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t459; t459; t371; t461; t458;];
tauJB = t1;
