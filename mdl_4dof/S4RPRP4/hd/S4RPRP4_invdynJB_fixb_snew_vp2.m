% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRP4
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:42
% EndTime: 2019-12-31 16:43:43
% DurationCPUTime: 0.84s
% Computational Cost: add. (4926->172), mult. (9135->212), div. (0->0), fcn. (4085->6), ass. (0->72)
t484 = Ifges(4,1) + Ifges(5,1);
t477 = Ifges(4,4) - Ifges(5,5);
t476 = -Ifges(4,5) - Ifges(5,4);
t483 = Ifges(4,2) + Ifges(5,3);
t475 = Ifges(4,6) - Ifges(5,6);
t482 = Ifges(4,3) + Ifges(5,2);
t452 = sin(qJ(3));
t454 = cos(qJ(3));
t427 = (-mrSges(5,1) * t454 - mrSges(5,3) * t452) * qJD(1);
t466 = qJD(1) * qJD(3);
t430 = qJDD(1) * t452 + t454 * t466;
t453 = sin(qJ(1));
t455 = cos(qJ(1));
t439 = t453 * g(1) - g(2) * t455;
t425 = qJDD(1) * pkin(1) + t439;
t440 = -g(1) * t455 - g(2) * t453;
t457 = qJD(1) ^ 2;
t429 = -pkin(1) * t457 + t440;
t450 = sin(pkin(6));
t451 = cos(pkin(6));
t409 = t450 * t425 + t451 * t429;
t406 = -pkin(2) * t457 + qJDD(1) * pkin(5) + t409;
t426 = (-pkin(3) * t454 - qJ(4) * t452) * qJD(1);
t456 = qJD(3) ^ 2;
t449 = -g(3) + qJDD(2);
t473 = t454 * t449;
t401 = -qJDD(3) * pkin(3) - t456 * qJ(4) - t473 + qJDD(4) + (qJD(1) * t426 + t406) * t452;
t467 = qJD(1) * t454;
t438 = mrSges(5,2) * t467 + qJD(3) * mrSges(5,3);
t461 = -m(5) * t401 + qJDD(3) * mrSges(5,1) + qJD(3) * t438;
t468 = qJD(1) * t452;
t397 = t430 * mrSges(5,2) + t427 * t468 - t461;
t403 = t454 * t406 + t452 * t449;
t400 = -pkin(3) * t456 + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t426 * t467 + t403;
t402 = -t406 * t452 + t473;
t431 = qJDD(1) * t454 - t452 * t466;
t436 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t468;
t462 = m(5) * t400 + qJDD(3) * mrSges(5,3) + qJD(3) * t436 + t427 * t467;
t469 = -t476 * qJD(3) + (t452 * t484 + t454 * t477) * qJD(1);
t471 = -t475 * qJD(3) + (-t452 * t477 - t454 * t483) * qJD(1);
t481 = -(t452 * t471 + t454 * t469) * qJD(1) + t482 * qJDD(3) - t476 * t430 + t475 * t431 + mrSges(4,1) * t402 - mrSges(5,1) * t401 - mrSges(4,2) * t403 + mrSges(5,3) * t400 - pkin(3) * t397 + qJ(4) * (mrSges(5,2) * t431 + t462);
t478 = mrSges(4,3) + mrSges(5,2);
t428 = (-mrSges(4,1) * t454 + mrSges(4,2) * t452) * qJD(1);
t435 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t468;
t393 = m(4) * t403 - qJDD(3) * mrSges(4,2) - qJD(3) * t435 + t428 * t467 + t431 * t478 + t462;
t437 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t467;
t394 = m(4) * t402 + qJDD(3) * mrSges(4,1) + qJD(3) * t437 - t478 * t430 + (-t427 - t428) * t468 + t461;
t463 = t454 * t393 - t394 * t452;
t383 = m(3) * t409 - mrSges(3,1) * t457 - qJDD(1) * mrSges(3,2) + t463;
t408 = t451 * t425 - t450 * t429;
t405 = -qJDD(1) * pkin(2) - t457 * pkin(5) - t408;
t398 = -t431 * pkin(3) - t430 * qJ(4) + (-0.2e1 * qJD(4) * t452 + (pkin(3) * t452 - qJ(4) * t454) * qJD(3)) * qJD(1) + t405;
t395 = m(5) * t398 - mrSges(5,1) * t431 - t430 * mrSges(5,3) - t436 * t468 - t438 * t467;
t458 = -m(4) * t405 + t431 * mrSges(4,1) - mrSges(4,2) * t430 - t435 * t468 + t437 * t467 - t395;
t388 = m(3) * t408 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t457 + t458;
t376 = t383 * t450 + t451 * t388;
t373 = m(2) * t439 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t457 + t376;
t464 = t383 * t451 - t388 * t450;
t374 = m(2) * t440 - mrSges(2,1) * t457 - qJDD(1) * mrSges(2,2) + t464;
t472 = t373 * t455 + t374 * t453;
t386 = t452 * t393 + t454 * t394;
t470 = t482 * qJD(3) + (-t452 * t476 + t454 * t475) * qJD(1);
t384 = m(3) * t449 + t386;
t465 = -t373 * t453 + t374 * t455;
t379 = -mrSges(4,1) * t405 - mrSges(5,1) * t398 + mrSges(5,2) * t400 + mrSges(4,3) * t403 - pkin(3) * t395 + t469 * qJD(3) + t475 * qJDD(3) + t477 * t430 + t431 * t483 - t470 * t468;
t380 = mrSges(4,2) * t405 + mrSges(5,2) * t401 - mrSges(4,3) * t402 - mrSges(5,3) * t398 - qJ(4) * t395 + t471 * qJD(3) - t476 * qJDD(3) + t430 * t484 + t477 * t431 + t470 * t467;
t460 = mrSges(2,1) * t439 + mrSges(3,1) * t408 - mrSges(2,2) * t440 - mrSges(3,2) * t409 + pkin(1) * t376 + pkin(2) * t458 + pkin(5) * t463 + t454 * t379 + t452 * t380 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t369 = -mrSges(3,1) * t449 + mrSges(3,3) * t409 + t457 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t386 - t481;
t368 = mrSges(3,2) * t449 - mrSges(3,3) * t408 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t457 - pkin(5) * t386 - t379 * t452 + t380 * t454;
t367 = -mrSges(2,2) * g(3) - mrSges(2,3) * t439 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t457 - qJ(2) * t376 + t368 * t451 - t369 * t450;
t366 = mrSges(2,1) * g(3) + mrSges(2,3) * t440 + t457 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t384 + qJ(2) * t464 + t450 * t368 + t451 * t369;
t1 = [-m(1) * g(1) + t465; -m(1) * g(2) + t472; (-m(1) - m(2)) * g(3) + t384; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t472 - t366 * t453 + t367 * t455; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t465 + t455 * t366 + t453 * t367; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t460; t460; t384; t481; t397;];
tauJB = t1;
