% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:10
% EndTime: 2019-12-31 18:01:12
% DurationCPUTime: 1.33s
% Computational Cost: add. (16540->182), mult. (22883->217), div. (0->0), fcn. (8210->8), ass. (0->79)
t454 = qJD(1) ^ 2;
t450 = sin(qJ(1));
t453 = cos(qJ(1));
t432 = -t453 * g(1) - t450 * g(2);
t459 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t432;
t473 = -pkin(1) - pkin(2);
t414 = t473 * t454 + t459;
t431 = t450 * g(1) - t453 * g(2);
t458 = -t454 * qJ(2) + qJDD(2) - t431;
t417 = t473 * qJDD(1) + t458;
t446 = sin(pkin(8));
t447 = cos(pkin(8));
t409 = -t446 * t414 + t447 * t417;
t406 = -qJDD(1) * pkin(3) + t409;
t410 = t447 * t414 + t446 * t417;
t407 = -t454 * pkin(3) + t410;
t449 = sin(qJ(4));
t452 = cos(qJ(4));
t403 = t449 * t406 + t452 * t407;
t437 = -qJD(1) + qJD(4);
t435 = t437 ^ 2;
t436 = -qJDD(1) + qJDD(4);
t400 = -(t435 * pkin(4)) + t436 * pkin(7) + t403;
t444 = g(3) + qJDD(3);
t448 = sin(qJ(5));
t451 = cos(qJ(5));
t397 = -t448 * t400 + t451 * t444;
t425 = (-mrSges(6,1) * t451 + mrSges(6,2) * t448) * t437;
t466 = qJD(5) * t437;
t426 = t448 * t436 + t451 * t466;
t468 = t437 * t451;
t429 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t468;
t469 = t437 * t448;
t395 = m(6) * t397 + qJDD(5) * mrSges(6,1) - t426 * mrSges(6,3) + qJD(5) * t429 - t425 * t469;
t398 = t451 * t400 + t448 * t444;
t427 = t451 * t436 - t448 * t466;
t428 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t469;
t396 = m(6) * t398 - qJDD(5) * mrSges(6,2) + t427 * mrSges(6,3) - qJD(5) * t428 + t425 * t468;
t385 = t451 * t395 + t448 * t396;
t384 = (m(4) + m(5)) * t444 + t385;
t419 = (Ifges(6,6) * qJD(5)) + (Ifges(6,4) * t448 + Ifges(6,2) * t451) * t437;
t420 = (Ifges(6,5) * qJD(5)) + (Ifges(6,1) * t448 + Ifges(6,4) * t451) * t437;
t475 = mrSges(6,1) * t397 - mrSges(6,2) * t398 + Ifges(6,5) * t426 + Ifges(6,6) * t427 + Ifges(6,3) * qJDD(5) + (t419 * t448 - t420 * t451) * t437;
t472 = mrSges(2,1) + mrSges(3,1);
t471 = Ifges(3,4) + Ifges(2,5);
t470 = Ifges(2,6) - Ifges(3,6);
t421 = -t454 * pkin(1) + t459;
t386 = -t448 * t395 + t451 * t396;
t382 = m(5) * t403 - (t435 * mrSges(5,1)) - t436 * mrSges(5,2) + t386;
t402 = t452 * t406 - t449 * t407;
t399 = -t436 * pkin(4) - t435 * pkin(7) - t402;
t391 = -m(6) * t399 + t427 * mrSges(6,1) - t426 * mrSges(6,2) - t428 * t469 + t429 * t468;
t390 = m(5) * t402 + t436 * mrSges(5,1) - t435 * mrSges(5,2) + t391;
t379 = t449 * t382 + t452 * t390;
t376 = m(4) * t409 - qJDD(1) * mrSges(4,1) - t454 * mrSges(4,2) + t379;
t463 = t452 * t382 - t449 * t390;
t377 = m(4) * t410 - t454 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t463;
t464 = -t446 * t376 + t447 * t377;
t460 = m(3) * t421 + qJDD(1) * mrSges(3,3) + t464;
t368 = m(2) * t432 - qJDD(1) * mrSges(2,2) - t472 * t454 + t460;
t373 = t447 * t376 + t446 * t377;
t424 = -qJDD(1) * pkin(1) + t458;
t372 = m(3) * t424 - qJDD(1) * mrSges(3,1) - t454 * mrSges(3,3) + t373;
t369 = m(2) * t431 + qJDD(1) * mrSges(2,1) - t454 * mrSges(2,2) - t372;
t467 = t450 * t368 + t453 * t369;
t465 = t453 * t368 - t450 * t369;
t418 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t448 + Ifges(6,6) * t451) * t437;
t387 = -mrSges(6,1) * t399 + mrSges(6,3) * t398 + Ifges(6,4) * t426 + Ifges(6,2) * t427 + Ifges(6,6) * qJDD(5) + qJD(5) * t420 - t418 * t469;
t388 = mrSges(6,2) * t399 - mrSges(6,3) * t397 + Ifges(6,1) * t426 + Ifges(6,4) * t427 + Ifges(6,5) * qJDD(5) - qJD(5) * t419 + t418 * t468;
t456 = mrSges(5,1) * t402 - mrSges(5,2) * t403 + Ifges(5,3) * t436 + pkin(4) * t391 + pkin(7) * t386 + t451 * t387 + t448 * t388;
t455 = -mrSges(3,1) * t424 - mrSges(4,1) * t409 - mrSges(2,2) * t432 - pkin(2) * t373 - pkin(3) * t379 + qJ(2) * (-t454 * mrSges(3,1) + t460) - pkin(1) * t372 + mrSges(4,2) * t410 + mrSges(3,3) * t421 + mrSges(2,1) * t431 - t456 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t383 = -m(3) * g(3) - t384;
t378 = -mrSges(5,1) * t444 + mrSges(5,3) * t403 + (t435 * Ifges(5,5)) + Ifges(5,6) * t436 - pkin(4) * t385 - t475;
t374 = mrSges(5,2) * t444 - mrSges(5,3) * t402 + Ifges(5,5) * t436 - t435 * Ifges(5,6) - pkin(7) * t385 - t448 * t387 + t451 * t388;
t364 = mrSges(4,2) * t444 - mrSges(4,3) * t409 - Ifges(4,5) * qJDD(1) - t454 * Ifges(4,6) - pkin(6) * t379 + t452 * t374 - t449 * t378;
t363 = -Ifges(4,6) * qJDD(1) + t454 * Ifges(4,5) - mrSges(4,1) * t444 + mrSges(4,3) * t410 + t449 * t374 + t452 * t378 - pkin(3) * (m(5) * t444 + t385) + pkin(6) * t463;
t362 = mrSges(3,2) * t424 - mrSges(2,3) * t431 - qJ(2) * t383 - qJ(3) * t373 - t446 * t363 + t447 * t364 - t470 * t454 + t471 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t361 = mrSges(3,2) * t421 + mrSges(2,3) * t432 - pkin(1) * t383 + pkin(2) * t384 + t472 * g(3) - qJ(3) * t464 + t470 * qJDD(1) - t447 * t363 - t446 * t364 + t471 * t454;
t1 = [-m(1) * g(1) + t465; -m(1) * g(2) + t467; (-m(1) - m(2) - m(3)) * g(3) - t384; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t467 - t450 * t361 + t453 * t362; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t465 + t453 * t361 + t450 * t362; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t455; t455; t372; t384; t456; t475;];
tauJB = t1;
