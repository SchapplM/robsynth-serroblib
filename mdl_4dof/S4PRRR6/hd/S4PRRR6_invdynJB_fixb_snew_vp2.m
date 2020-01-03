% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRR6
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:46
% EndTime: 2019-12-31 16:34:47
% DurationCPUTime: 1.12s
% Computational Cost: add. (10038->187), mult. (19668->244), div. (0->0), fcn. (11878->8), ass. (0->81)
t437 = sin(pkin(7));
t460 = cos(pkin(7));
t425 = -g(1) * t460 - g(2) * t437;
t436 = -g(3) + qJDD(1);
t440 = sin(qJ(2));
t443 = cos(qJ(2));
t409 = t443 * t425 + t440 * t436;
t444 = qJD(2) ^ 2;
t405 = -pkin(2) * t444 + qJDD(2) * pkin(5) + t409;
t424 = g(1) * t437 - g(2) * t460;
t439 = sin(qJ(3));
t442 = cos(qJ(3));
t396 = -t439 * t405 - t442 * t424;
t456 = qJD(2) * qJD(3);
t454 = t442 * t456;
t422 = qJDD(2) * t439 + t454;
t387 = (-t422 + t454) * pkin(6) + (t439 * t442 * t444 + qJDD(3)) * pkin(3) + t396;
t397 = t442 * t405 - t439 * t424;
t423 = qJDD(2) * t442 - t439 * t456;
t458 = qJD(2) * t439;
t428 = qJD(3) * pkin(3) - pkin(6) * t458;
t435 = t442 ^ 2;
t388 = -pkin(3) * t435 * t444 + pkin(6) * t423 - qJD(3) * t428 + t397;
t438 = sin(qJ(4));
t441 = cos(qJ(4));
t385 = t387 * t441 - t388 * t438;
t413 = (-t438 * t439 + t441 * t442) * qJD(2);
t395 = qJD(4) * t413 + t422 * t441 + t423 * t438;
t414 = (t438 * t442 + t439 * t441) * qJD(2);
t402 = -mrSges(5,1) * t413 + mrSges(5,2) * t414;
t434 = qJD(3) + qJD(4);
t406 = -mrSges(5,2) * t434 + mrSges(5,3) * t413;
t433 = qJDD(3) + qJDD(4);
t382 = m(5) * t385 + mrSges(5,1) * t433 - t395 * mrSges(5,3) - t402 * t414 + t406 * t434;
t386 = t387 * t438 + t388 * t441;
t394 = -qJD(4) * t414 - t422 * t438 + t423 * t441;
t407 = mrSges(5,1) * t434 - mrSges(5,3) * t414;
t383 = m(5) * t386 - mrSges(5,2) * t433 + t394 * mrSges(5,3) + t402 * t413 - t407 * t434;
t373 = t441 * t382 + t438 * t383;
t411 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t439 + Ifges(4,2) * t442) * qJD(2);
t412 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t439 + Ifges(4,4) * t442) * qJD(2);
t399 = Ifges(5,4) * t414 + Ifges(5,2) * t413 + Ifges(5,6) * t434;
t400 = Ifges(5,1) * t414 + Ifges(5,4) * t413 + Ifges(5,5) * t434;
t447 = -mrSges(5,1) * t385 + mrSges(5,2) * t386 - Ifges(5,5) * t395 - Ifges(5,6) * t394 - Ifges(5,3) * t433 - t414 * t399 + t413 * t400;
t461 = mrSges(4,1) * t396 - mrSges(4,2) * t397 + Ifges(4,5) * t422 + Ifges(4,6) * t423 + Ifges(4,3) * qJDD(3) + pkin(3) * t373 + (t411 * t439 - t412 * t442) * qJD(2) - t447;
t421 = (-mrSges(4,1) * t442 + mrSges(4,2) * t439) * qJD(2);
t457 = qJD(2) * t442;
t427 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t457;
t371 = m(4) * t396 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t422 + qJD(3) * t427 - t421 * t458 + t373;
t426 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t458;
t451 = -t382 * t438 + t441 * t383;
t372 = m(4) * t397 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t423 - qJD(3) * t426 + t421 * t457 + t451;
t369 = -t371 * t439 + t442 * t372;
t365 = m(3) * t409 - mrSges(3,1) * t444 - qJDD(2) * mrSges(3,2) + t369;
t408 = -t440 * t425 + t443 * t436;
t449 = -qJDD(2) * pkin(2) - t408;
t404 = -pkin(5) * t444 + t449;
t389 = t428 * t458 - t423 * pkin(3) + (-pkin(6) * t435 - pkin(5)) * t444 + t449;
t448 = m(5) * t389 - t394 * mrSges(5,1) + t395 * mrSges(5,2) - t413 * t406 + t407 * t414;
t378 = -m(4) * t404 + t423 * mrSges(4,1) - mrSges(4,2) * t422 - t426 * t458 + t427 * t457 - t448;
t377 = m(3) * t408 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t444 + t378;
t452 = t443 * t365 - t377 * t440;
t360 = m(2) * t425 + t452;
t368 = t442 * t371 + t439 * t372;
t367 = (m(2) + m(3)) * t424 - t368;
t459 = t437 * t360 + t460 * t367;
t361 = t440 * t365 + t443 * t377;
t455 = m(2) * t436 + t361;
t453 = t460 * t360 - t367 * t437;
t398 = Ifges(5,5) * t414 + Ifges(5,6) * t413 + Ifges(5,3) * t434;
t374 = -mrSges(5,1) * t389 + mrSges(5,3) * t386 + Ifges(5,4) * t395 + Ifges(5,2) * t394 + Ifges(5,6) * t433 - t398 * t414 + t400 * t434;
t375 = mrSges(5,2) * t389 - mrSges(5,3) * t385 + Ifges(5,1) * t395 + Ifges(5,4) * t394 + Ifges(5,5) * t433 + t398 * t413 - t399 * t434;
t410 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t439 + Ifges(4,6) * t442) * qJD(2);
t357 = -mrSges(4,1) * t404 + mrSges(4,3) * t397 + Ifges(4,4) * t422 + Ifges(4,2) * t423 + Ifges(4,6) * qJDD(3) - pkin(3) * t448 + pkin(6) * t451 + qJD(3) * t412 + t441 * t374 + t438 * t375 - t410 * t458;
t362 = mrSges(4,2) * t404 - mrSges(4,3) * t396 + Ifges(4,1) * t422 + Ifges(4,4) * t423 + Ifges(4,5) * qJDD(3) - pkin(6) * t373 - qJD(3) * t411 - t374 * t438 + t375 * t441 + t410 * t457;
t446 = mrSges(3,1) * t408 - mrSges(3,2) * t409 + Ifges(3,3) * qJDD(2) + pkin(2) * t378 + pkin(5) * t369 + t357 * t442 + t362 * t439;
t356 = mrSges(3,1) * t424 + mrSges(3,3) * t409 + t444 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t368 - t461;
t355 = -mrSges(3,2) * t424 - mrSges(3,3) * t408 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t444 - pkin(5) * t368 - t357 * t439 + t362 * t442;
t354 = -mrSges(2,1) * t436 + mrSges(2,3) * t425 - pkin(1) * t361 - t446;
t353 = mrSges(2,2) * t436 - mrSges(2,3) * t424 - pkin(4) * t361 + t355 * t443 - t356 * t440;
t1 = [-m(1) * g(1) + t453; -m(1) * g(2) + t459; -m(1) * g(3) + t455; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t459 + t353 * t460 - t354 * t437; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t453 + t437 * t353 + t354 * t460; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t424 - mrSges(2,2) * t425 + t440 * t355 + t443 * t356 + pkin(1) * (m(3) * t424 - t368) + pkin(4) * t452; t455; t446; t461; -t447;];
tauJB = t1;
