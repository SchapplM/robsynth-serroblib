% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:27
% EndTime: 2019-12-05 15:26:29
% DurationCPUTime: 0.92s
% Computational Cost: add. (8023->167), mult. (12451->203), div. (0->0), fcn. (6658->8), ass. (0->76)
t425 = sin(pkin(7));
t427 = cos(pkin(7));
t412 = -t427 * g(1) - t425 * g(2);
t421 = -g(3) + qJDD(1);
t429 = sin(qJ(2));
t431 = cos(qJ(2));
t393 = -t429 * t412 + t431 * t421;
t391 = qJDD(2) * pkin(2) + t393;
t394 = t431 * t412 + t429 * t421;
t432 = qJD(2) ^ 2;
t392 = -t432 * pkin(2) + t394;
t424 = sin(pkin(8));
t426 = cos(pkin(8));
t386 = t426 * t391 - t424 * t392;
t437 = -t432 * qJ(4) + qJDD(4) - t386;
t452 = -pkin(3) - pkin(6);
t383 = t452 * qJDD(2) + t437;
t411 = t425 * g(1) - t427 * g(2);
t410 = qJDD(3) - t411;
t428 = sin(qJ(5));
t430 = cos(qJ(5));
t379 = t430 * t383 - t428 * t410;
t407 = (mrSges(6,1) * t428 + mrSges(6,2) * t430) * qJD(2);
t444 = qJD(2) * qJD(5);
t409 = t430 * qJDD(2) - t428 * t444;
t446 = qJD(2) * t428;
t413 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t446;
t445 = qJD(2) * t430;
t376 = m(6) * t379 + qJDD(5) * mrSges(6,1) - t409 * mrSges(6,3) + qJD(5) * t413 - t407 * t445;
t380 = t428 * t383 + t430 * t410;
t408 = -t428 * qJDD(2) - t430 * t444;
t414 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t445;
t377 = m(6) * t380 - qJDD(5) * mrSges(6,2) + t408 * mrSges(6,3) - qJD(5) * t414 - t407 * t446;
t367 = t430 * t376 + t428 * t377;
t385 = -qJDD(2) * pkin(3) + t437;
t435 = -m(5) * t385 + t432 * mrSges(5,3) - t367;
t451 = mrSges(4,1) - mrSges(5,2);
t361 = m(4) * t386 - t432 * mrSges(4,2) + t451 * qJDD(2) + t435;
t387 = t424 * t391 + t426 * t392;
t436 = qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) + t387;
t384 = t432 * pkin(3) - t436;
t382 = t452 * t432 + t436;
t438 = -m(6) * t382 + t408 * mrSges(6,1) - t409 * mrSges(6,2) - t413 * t446 - t414 * t445;
t373 = -m(5) * t384 + t432 * mrSges(5,2) + qJDD(2) * mrSges(5,3) - t438;
t370 = m(4) * t387 - t432 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t373;
t359 = t426 * t361 + t424 * t370;
t362 = qJDD(2) * mrSges(5,2) - t435;
t397 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t430 - Ifges(6,6) * t428) * qJD(2);
t399 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t430 - Ifges(6,4) * t428) * qJD(2);
t371 = -mrSges(6,1) * t382 + mrSges(6,3) * t380 + Ifges(6,4) * t409 + Ifges(6,2) * t408 + Ifges(6,6) * qJDD(5) + qJD(5) * t399 - t397 * t445;
t398 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t430 - Ifges(6,2) * t428) * qJD(2);
t372 = mrSges(6,2) * t382 - mrSges(6,3) * t379 + Ifges(6,1) * t409 + Ifges(6,4) * t408 + Ifges(6,5) * qJDD(5) - qJD(5) * t398 - t397 * t446;
t453 = mrSges(3,1) * t393 + mrSges(4,1) * t386 - mrSges(3,2) * t394 - mrSges(4,2) * t387 + mrSges(5,2) * t385 - mrSges(5,3) * t384 + pkin(2) * t359 - pkin(3) * t362 - pkin(6) * t367 + qJ(4) * t373 - t428 * t371 + t430 * t372 + (Ifges(3,3) + Ifges(4,3) + Ifges(5,1)) * qJDD(2);
t450 = -Ifges(5,4) + Ifges(4,5);
t449 = Ifges(5,5) - Ifges(4,6);
t357 = m(3) * t393 + qJDD(2) * mrSges(3,1) - t432 * mrSges(3,2) + t359;
t439 = -t424 * t361 + t426 * t370;
t358 = m(3) * t394 - t432 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t439;
t440 = -t429 * t357 + t431 * t358;
t350 = m(2) * t412 + t440;
t447 = -t428 * t376 + t430 * t377;
t366 = m(5) * t410 + t447;
t365 = m(4) * t410 + t366;
t364 = (m(2) + m(3)) * t411 - t365;
t448 = t425 * t350 + t427 * t364;
t351 = t431 * t357 + t429 * t358;
t442 = m(2) * t421 + t351;
t441 = t427 * t350 - t425 * t364;
t434 = mrSges(6,1) * t379 - mrSges(6,2) * t380 + Ifges(6,5) * t409 + Ifges(6,6) * t408 + Ifges(6,3) * qJDD(5) + t398 * t445 + t399 * t446;
t353 = mrSges(5,1) * t385 - mrSges(4,3) * t386 + pkin(4) * t367 - qJ(4) * t366 + t449 * t432 + (mrSges(4,2) - mrSges(5,3)) * t410 + t450 * qJDD(2) + t434;
t352 = -mrSges(5,1) * t384 + mrSges(4,3) * t387 - pkin(3) * t366 - pkin(4) * t438 - pkin(6) * t447 - t449 * qJDD(2) - t430 * t371 - t428 * t372 - t451 * t410 + t450 * t432;
t347 = -mrSges(3,2) * t411 - mrSges(3,3) * t393 + Ifges(3,5) * qJDD(2) - t432 * Ifges(3,6) - qJ(3) * t359 - t424 * t352 + t426 * t353;
t346 = mrSges(3,1) * t411 + mrSges(3,3) * t394 + t432 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t365 + qJ(3) * t439 + t426 * t352 + t424 * t353;
t345 = -mrSges(2,1) * t421 + mrSges(2,3) * t412 - pkin(1) * t351 - t453;
t344 = mrSges(2,2) * t421 - mrSges(2,3) * t411 - pkin(5) * t351 - t429 * t346 + t431 * t347;
t1 = [-m(1) * g(1) + t441; -m(1) * g(2) + t448; -m(1) * g(3) + t442; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t448 + t427 * t344 - t425 * t345; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t441 + t425 * t344 + t427 * t345; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t411 - mrSges(2,2) * t412 + t429 * t347 + t431 * t346 + pkin(1) * (m(3) * t411 - t365) + pkin(5) * t440; t442; t453; t365; t362; t434;];
tauJB = t1;
