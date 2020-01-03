% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:17
% EndTime: 2019-12-31 16:30:18
% DurationCPUTime: 0.69s
% Computational Cost: add. (3928->162), mult. (7291->203), div. (0->0), fcn. (3509->6), ass. (0->70)
t454 = Ifges(4,1) + Ifges(5,1);
t447 = Ifges(4,4) - Ifges(5,5);
t446 = Ifges(5,4) + Ifges(4,5);
t453 = Ifges(4,2) + Ifges(5,3);
t445 = Ifges(5,6) - Ifges(4,6);
t452 = Ifges(4,3) + Ifges(5,2);
t422 = sin(qJ(3));
t424 = cos(qJ(3));
t402 = (-mrSges(5,1) * t424 - mrSges(5,3) * t422) * qJD(2);
t435 = qJD(2) * qJD(3);
t404 = t422 * qJDD(2) + t424 * t435;
t421 = sin(pkin(6));
t443 = cos(pkin(6));
t409 = -t443 * g(1) - t421 * g(2);
t420 = -g(3) + qJDD(1);
t423 = sin(qJ(2));
t425 = cos(qJ(2));
t385 = t425 * t409 + t423 * t420;
t427 = qJD(2) ^ 2;
t383 = -t427 * pkin(2) + qJDD(2) * pkin(5) + t385;
t401 = (-pkin(3) * t424 - qJ(4) * t422) * qJD(2);
t426 = qJD(3) ^ 2;
t408 = t421 * g(1) - t443 * g(2);
t442 = t424 * t408;
t378 = -qJDD(3) * pkin(3) - t426 * qJ(4) + t442 + qJDD(4) + (qJD(2) * t401 + t383) * t422;
t436 = qJD(2) * t424;
t413 = mrSges(5,2) * t436 + qJD(3) * mrSges(5,3);
t430 = -m(5) * t378 + qJDD(3) * mrSges(5,1) + qJD(3) * t413;
t437 = qJD(2) * t422;
t374 = t404 * mrSges(5,2) + t402 * t437 - t430;
t380 = t424 * t383 - t422 * t408;
t377 = -t426 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t401 * t436 + t380;
t379 = -t422 * t383 - t442;
t405 = t424 * qJDD(2) - t422 * t435;
t411 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t437;
t431 = m(5) * t377 + qJDD(3) * mrSges(5,3) + qJD(3) * t411 + t402 * t436;
t438 = t446 * qJD(3) + (t454 * t422 + t447 * t424) * qJD(2);
t440 = t445 * qJD(3) + (-t447 * t422 - t453 * t424) * qJD(2);
t451 = -(t440 * t422 + t438 * t424) * qJD(2) + t452 * qJDD(3) + t446 * t404 - t445 * t405 + mrSges(4,1) * t379 - mrSges(5,1) * t378 - mrSges(4,2) * t380 + mrSges(5,3) * t377 - pkin(3) * t374 + qJ(4) * (t405 * mrSges(5,2) + t431);
t448 = mrSges(4,3) + mrSges(5,2);
t403 = (-mrSges(4,1) * t424 + mrSges(4,2) * t422) * qJD(2);
t410 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t437;
t370 = m(4) * t380 - qJDD(3) * mrSges(4,2) - qJD(3) * t410 + t403 * t436 + t448 * t405 + t431;
t412 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t436;
t371 = m(4) * t379 + qJDD(3) * mrSges(4,1) + qJD(3) * t412 - t448 * t404 + (-t402 - t403) * t437 + t430;
t365 = t424 * t370 - t422 * t371;
t361 = m(3) * t385 - t427 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t365;
t384 = -t423 * t409 + t425 * t420;
t382 = -qJDD(2) * pkin(2) - t427 * pkin(5) - t384;
t375 = -t405 * pkin(3) - t404 * qJ(4) + (-0.2e1 * qJD(4) * t422 + (pkin(3) * t422 - qJ(4) * t424) * qJD(3)) * qJD(2) + t382;
t372 = m(5) * t375 - t405 * mrSges(5,1) - t404 * mrSges(5,3) - t411 * t437 - t413 * t436;
t368 = -m(4) * t382 + t405 * mrSges(4,1) - t404 * mrSges(4,2) - t410 * t437 + t412 * t436 - t372;
t367 = m(3) * t384 + qJDD(2) * mrSges(3,1) - t427 * mrSges(3,2) + t368;
t432 = t425 * t361 - t423 * t367;
t355 = m(2) * t409 + t432;
t364 = t422 * t370 + t424 * t371;
t363 = (m(2) + m(3)) * t408 - t364;
t441 = t421 * t355 + t443 * t363;
t356 = t423 * t361 + t425 * t367;
t439 = t452 * qJD(3) + (t446 * t422 - t445 * t424) * qJD(2);
t434 = m(2) * t420 + t356;
t433 = t443 * t355 - t421 * t363;
t357 = -mrSges(4,1) * t382 - mrSges(5,1) * t375 + mrSges(5,2) * t377 + mrSges(4,3) * t380 - pkin(3) * t372 + t438 * qJD(3) - t445 * qJDD(3) + t447 * t404 + t453 * t405 - t439 * t437;
t358 = mrSges(4,2) * t382 + mrSges(5,2) * t378 - mrSges(4,3) * t379 - mrSges(5,3) * t375 - qJ(4) * t372 + t440 * qJD(3) + t446 * qJDD(3) + t454 * t404 + t447 * t405 + t439 * t436;
t428 = mrSges(3,1) * t384 - mrSges(3,2) * t385 + Ifges(3,3) * qJDD(2) + pkin(2) * t368 + pkin(5) * t365 + t424 * t357 + t422 * t358;
t352 = mrSges(3,1) * t408 + mrSges(3,3) * t385 + t427 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t364 - t451;
t351 = -mrSges(3,2) * t408 - mrSges(3,3) * t384 + Ifges(3,5) * qJDD(2) - t427 * Ifges(3,6) - pkin(5) * t364 - t422 * t357 + t424 * t358;
t350 = -mrSges(2,1) * t420 + mrSges(2,3) * t409 - pkin(1) * t356 - t428;
t349 = mrSges(2,2) * t420 - mrSges(2,3) * t408 - pkin(4) * t356 + t425 * t351 - t423 * t352;
t1 = [-m(1) * g(1) + t433; -m(1) * g(2) + t441; -m(1) * g(3) + t434; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t441 + t443 * t349 - t421 * t350; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t433 + t421 * t349 + t443 * t350; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t408 - mrSges(2,2) * t409 + t423 * t351 + t425 * t352 + pkin(1) * (m(3) * t408 - t364) + pkin(4) * t432; t434; t428; t451; t374;];
tauJB = t1;
