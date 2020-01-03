% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:34
% EndTime: 2019-12-31 17:56:35
% DurationCPUTime: 1.16s
% Computational Cost: add. (13687->181), mult. (19070->217), div. (0->0), fcn. (8062->8), ass. (0->79)
t446 = qJD(1) ^ 2;
t442 = sin(qJ(1));
t445 = cos(qJ(1));
t421 = t442 * g(1) - t445 * g(2);
t416 = qJDD(1) * pkin(1) + t421;
t422 = -t445 * g(1) - t442 * g(2);
t417 = -t446 * pkin(1) + t422;
t438 = sin(pkin(8));
t439 = cos(pkin(8));
t404 = t438 * t416 + t439 * t417;
t453 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t404;
t464 = -pkin(2) - pkin(3);
t396 = t464 * t446 + t453;
t403 = t439 * t416 - t438 * t417;
t450 = -t446 * qJ(3) + qJDD(3) - t403;
t399 = t464 * qJDD(1) + t450;
t441 = sin(qJ(4));
t444 = cos(qJ(4));
t393 = t444 * t396 + t441 * t399;
t429 = -qJD(1) + qJD(4);
t427 = t429 ^ 2;
t428 = -qJDD(1) + qJDD(4);
t390 = -(t427 * pkin(4)) + t428 * pkin(7) + t393;
t436 = g(3) - qJDD(2);
t440 = sin(qJ(5));
t443 = cos(qJ(5));
t387 = -t440 * t390 + t443 * t436;
t388 = t443 * t390 + t440 * t436;
t406 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t440 + Ifges(6,2) * t443) * t429;
t407 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t440 + Ifges(6,4) * t443) * t429;
t457 = qJD(5) * t429;
t411 = t440 * t428 + t443 * t457;
t412 = t443 * t428 - t440 * t457;
t465 = mrSges(6,1) * t387 - mrSges(6,2) * t388 + Ifges(6,5) * t411 + Ifges(6,6) * t412 + Ifges(6,3) * qJDD(5) + (t406 * t440 - t407 * t443) * t429;
t463 = -mrSges(3,1) - mrSges(4,1);
t462 = Ifges(4,4) + Ifges(3,5);
t461 = Ifges(3,6) - Ifges(4,6);
t460 = t429 * t440;
t459 = t429 * t443;
t400 = -t446 * pkin(2) + t453;
t410 = (-mrSges(6,1) * t443 + mrSges(6,2) * t440) * t429;
t419 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t459;
t385 = m(6) * t387 + qJDD(5) * mrSges(6,1) - t411 * mrSges(6,3) + qJD(5) * t419 - t410 * t460;
t418 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t460;
t386 = m(6) * t388 - qJDD(5) * mrSges(6,2) + t412 * mrSges(6,3) - qJD(5) * t418 + t410 * t459;
t379 = -t440 * t385 + t443 * t386;
t375 = m(5) * t393 - (t427 * mrSges(5,1)) - t428 * mrSges(5,2) + t379;
t392 = -t441 * t396 + t444 * t399;
t389 = -t428 * pkin(4) - t427 * pkin(7) - t392;
t383 = -m(6) * t389 + t412 * mrSges(6,1) - t411 * mrSges(6,2) - t418 * t460 + t419 * t459;
t382 = m(5) * t392 + t428 * mrSges(5,1) - t427 * mrSges(5,2) + t383;
t454 = t444 * t375 - t441 * t382;
t451 = m(4) * t400 + qJDD(1) * mrSges(4,3) + t454;
t367 = m(3) * t404 - qJDD(1) * mrSges(3,2) + t463 * t446 + t451;
t373 = t441 * t375 + t444 * t382;
t401 = -qJDD(1) * pkin(2) + t450;
t371 = m(4) * t401 - qJDD(1) * mrSges(4,1) - t446 * mrSges(4,3) + t373;
t368 = m(3) * t403 + qJDD(1) * mrSges(3,1) - t446 * mrSges(3,2) - t371;
t362 = t438 * t367 + t439 * t368;
t359 = m(2) * t421 + qJDD(1) * mrSges(2,1) - t446 * mrSges(2,2) + t362;
t455 = t439 * t367 - t438 * t368;
t360 = m(2) * t422 - t446 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t455;
t458 = t445 * t359 + t442 * t360;
t456 = -t442 * t359 + t445 * t360;
t378 = t443 * t385 + t440 * t386;
t377 = -t378 + (-m(4) - m(5)) * t436;
t376 = -m(3) * t436 + t377;
t405 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t440 + Ifges(6,6) * t443) * t429;
t380 = -mrSges(6,1) * t389 + mrSges(6,3) * t388 + Ifges(6,4) * t411 + Ifges(6,2) * t412 + Ifges(6,6) * qJDD(5) + qJD(5) * t407 - t405 * t460;
t381 = mrSges(6,2) * t389 - mrSges(6,3) * t387 + Ifges(6,1) * t411 + Ifges(6,4) * t412 + Ifges(6,5) * qJDD(5) - qJD(5) * t406 + t405 * t459;
t448 = mrSges(5,1) * t392 - mrSges(5,2) * t393 + Ifges(5,3) * t428 + pkin(4) * t383 + pkin(7) * t379 + t443 * t380 + t440 * t381;
t447 = -mrSges(4,1) * t401 - mrSges(2,2) * t422 - mrSges(3,2) * t404 + pkin(1) * t362 - pkin(3) * t373 + qJ(3) * (-t446 * mrSges(4,1) + t451) - pkin(2) * t371 + mrSges(4,3) * t400 + mrSges(3,1) * t403 + mrSges(2,1) * t421 - t448 + (Ifges(2,3) + Ifges(3,3) + Ifges(4,2)) * qJDD(1);
t372 = -mrSges(5,1) * t436 + mrSges(5,3) * t393 + (t427 * Ifges(5,5)) + Ifges(5,6) * t428 - pkin(4) * t378 - t465;
t363 = mrSges(5,2) * t436 - mrSges(5,3) * t392 + Ifges(5,5) * t428 - t427 * Ifges(5,6) - pkin(7) * t378 - t440 * t380 + t443 * t381;
t355 = mrSges(4,2) * t401 - mrSges(3,3) * t403 - pkin(6) * t373 - qJ(3) * t377 + t444 * t363 - t441 * t372 - t461 * t446 + (-mrSges(3,2) + mrSges(4,3)) * t436 + t462 * qJDD(1);
t354 = mrSges(3,3) * t404 + mrSges(4,2) * t400 - t441 * t363 - t444 * t372 + pkin(3) * t378 - pkin(6) * t454 - pkin(2) * t377 + t462 * t446 + (pkin(3) * m(5) - t463) * t436 + t461 * qJDD(1);
t353 = -mrSges(2,2) * g(3) - mrSges(2,3) * t421 + Ifges(2,5) * qJDD(1) - t446 * Ifges(2,6) - qJ(2) * t362 - t438 * t354 + t439 * t355;
t352 = mrSges(2,1) * g(3) + mrSges(2,3) * t422 + t446 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t376 + qJ(2) * t455 + t439 * t354 + t438 * t355;
t1 = [-m(1) * g(1) + t456; -m(1) * g(2) + t458; (-m(1) - m(2)) * g(3) + t376; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t458 - t442 * t352 + t445 * t353; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t456 + t445 * t352 + t442 * t353; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t447; t447; t376; t371; t448; t465;];
tauJB = t1;
