% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRP4
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:50
% EndTime: 2019-12-31 16:27:50
% DurationCPUTime: 0.71s
% Computational Cost: add. (4236->161), mult. (8041->203), div. (0->0), fcn. (3968->6), ass. (0->71)
t464 = Ifges(4,1) + Ifges(5,1);
t457 = Ifges(4,4) - Ifges(5,5);
t456 = Ifges(5,4) + Ifges(4,5);
t463 = Ifges(4,2) + Ifges(5,3);
t455 = Ifges(5,6) - Ifges(4,6);
t462 = Ifges(4,3) + Ifges(5,2);
t430 = sin(qJ(3));
t432 = cos(qJ(3));
t408 = (-mrSges(5,1) * t432 - mrSges(5,3) * t430) * qJD(2);
t446 = qJD(2) * qJD(3);
t410 = t430 * qJDD(2) + t432 * t446;
t428 = sin(pkin(6));
t429 = cos(pkin(6));
t414 = t428 * g(1) - t429 * g(2);
t415 = -t429 * g(1) - t428 * g(2);
t431 = sin(qJ(2));
t433 = cos(qJ(2));
t391 = t431 * t414 + t433 * t415;
t435 = qJD(2) ^ 2;
t388 = -t435 * pkin(2) + qJDD(2) * pkin(5) + t391;
t407 = (-pkin(3) * t432 - qJ(4) * t430) * qJD(2);
t434 = qJD(3) ^ 2;
t427 = -g(3) + qJDD(1);
t453 = t432 * t427;
t383 = -qJDD(3) * pkin(3) - t434 * qJ(4) - t453 + qJDD(4) + (qJD(2) * t407 + t388) * t430;
t447 = qJD(2) * t432;
t419 = mrSges(5,2) * t447 + qJD(3) * mrSges(5,3);
t439 = -m(5) * t383 + qJDD(3) * mrSges(5,1) + qJD(3) * t419;
t448 = qJD(2) * t430;
t379 = t410 * mrSges(5,2) + t408 * t448 - t439;
t385 = t432 * t388 + t430 * t427;
t382 = -t434 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t407 * t447 + t385;
t384 = -t430 * t388 + t453;
t411 = t432 * qJDD(2) - t430 * t446;
t417 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t448;
t441 = m(5) * t382 + qJDD(3) * mrSges(5,3) + qJD(3) * t417 + t408 * t447;
t449 = t456 * qJD(3) + (t464 * t430 + t457 * t432) * qJD(2);
t451 = t455 * qJD(3) + (-t457 * t430 - t463 * t432) * qJD(2);
t461 = -(t430 * t451 + t432 * t449) * qJD(2) + t462 * qJDD(3) + t456 * t410 - t455 * t411 + mrSges(4,1) * t384 - mrSges(5,1) * t383 - mrSges(4,2) * t385 + mrSges(5,3) * t382 - pkin(3) * t379 + qJ(4) * (t411 * mrSges(5,2) + t441);
t458 = mrSges(4,3) + mrSges(5,2);
t409 = (-mrSges(4,1) * t432 + mrSges(4,2) * t430) * qJD(2);
t416 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t448;
t375 = m(4) * t385 - qJDD(3) * mrSges(4,2) - qJD(3) * t416 + t409 * t447 + t458 * t411 + t441;
t418 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t447;
t376 = m(4) * t384 + qJDD(3) * mrSges(4,1) + qJD(3) * t418 - t458 * t410 + (-t408 - t409) * t448 + t439;
t442 = t432 * t375 - t430 * t376;
t366 = m(3) * t391 - t435 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t442;
t390 = t433 * t414 - t431 * t415;
t387 = -qJDD(2) * pkin(2) - t435 * pkin(5) - t390;
t380 = -t411 * pkin(3) - t410 * qJ(4) + (-0.2e1 * qJD(4) * t430 + (pkin(3) * t430 - qJ(4) * t432) * qJD(3)) * qJD(2) + t387;
t377 = m(5) * t380 - t411 * mrSges(5,1) - t410 * mrSges(5,3) - t417 * t448 - t419 * t447;
t436 = -m(4) * t387 + t411 * mrSges(4,1) - t410 * mrSges(4,2) - t416 * t448 + t418 * t447 - t377;
t370 = m(3) * t390 + qJDD(2) * mrSges(3,1) - t435 * mrSges(3,2) + t436;
t359 = t431 * t366 + t433 * t370;
t357 = m(2) * t414 + t359;
t443 = t433 * t366 - t431 * t370;
t358 = m(2) * t415 + t443;
t452 = t429 * t357 + t428 * t358;
t368 = t430 * t375 + t432 * t376;
t450 = t462 * qJD(3) + (t456 * t430 - t455 * t432) * qJD(2);
t445 = m(3) * t427 + t368;
t444 = -t428 * t357 + t429 * t358;
t440 = m(2) * t427 + t445;
t362 = -mrSges(4,1) * t387 - mrSges(5,1) * t380 + mrSges(5,2) * t382 + mrSges(4,3) * t385 - pkin(3) * t377 + t449 * qJD(3) - t455 * qJDD(3) + t457 * t410 + t463 * t411 - t450 * t448;
t363 = mrSges(4,2) * t387 + mrSges(5,2) * t383 - mrSges(4,3) * t384 - mrSges(5,3) * t380 - qJ(4) * t377 + t451 * qJD(3) + t456 * qJDD(3) + t464 * t410 + t457 * t411 + t450 * t447;
t438 = mrSges(3,1) * t390 - mrSges(3,2) * t391 + Ifges(3,3) * qJDD(2) + pkin(2) * t436 + pkin(5) * t442 + t432 * t362 + t430 * t363;
t353 = -mrSges(3,1) * t427 + mrSges(3,3) * t391 + t435 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t368 - t461;
t352 = mrSges(3,2) * t427 - mrSges(3,3) * t390 + Ifges(3,5) * qJDD(2) - t435 * Ifges(3,6) - pkin(5) * t368 - t430 * t362 + t432 * t363;
t351 = mrSges(2,2) * t427 - mrSges(2,3) * t414 - pkin(4) * t359 + t433 * t352 - t431 * t353;
t350 = -mrSges(2,1) * t427 + mrSges(2,3) * t415 - pkin(1) * t445 + pkin(4) * t443 + t431 * t352 + t433 * t353;
t1 = [-m(1) * g(1) + t444; -m(1) * g(2) + t452; -m(1) * g(3) + t440; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t452 - t428 * t350 + t429 * t351; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t444 + t429 * t350 + t428 * t351; -mrSges(1,1) * g(2) + mrSges(2,1) * t414 + mrSges(1,2) * g(1) - mrSges(2,2) * t415 + pkin(1) * t359 + t438; t440; t438; t461; t379;];
tauJB = t1;
