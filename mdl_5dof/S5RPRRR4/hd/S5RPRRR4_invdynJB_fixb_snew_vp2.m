% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:16
% EndTime: 2020-01-03 11:52:19
% DurationCPUTime: 2.84s
% Computational Cost: add. (40644->185), mult. (52826->233), div. (0->0), fcn. (28272->10), ass. (0->86)
t499 = sin(qJ(1));
t503 = cos(qJ(1));
t478 = -t503 * g(2) - t499 * g(3);
t474 = qJDD(1) * pkin(1) + t478;
t477 = -t499 * g(2) + t503 * g(3);
t504 = qJD(1) ^ 2;
t475 = -t504 * pkin(1) + t477;
t494 = sin(pkin(9));
t495 = cos(pkin(9));
t459 = t495 * t474 - t494 * t475;
t456 = qJDD(1) * pkin(2) + t459;
t460 = t494 * t474 + t495 * t475;
t457 = -t504 * pkin(2) + t460;
t498 = sin(qJ(3));
t502 = cos(qJ(3));
t451 = t502 * t456 - t498 * t457;
t489 = qJDD(1) + qJDD(3);
t448 = t489 * pkin(3) + t451;
t452 = t498 * t456 + t502 * t457;
t490 = qJD(1) + qJD(3);
t488 = t490 ^ 2;
t449 = -t488 * pkin(3) + t452;
t497 = sin(qJ(4));
t501 = cos(qJ(4));
t445 = t497 * t448 + t501 * t449;
t483 = qJD(4) + t490;
t481 = t483 ^ 2;
t482 = qJDD(4) + t489;
t442 = -t481 * pkin(4) + t482 * pkin(8) + t445;
t493 = -g(1) + qJDD(2);
t496 = sin(qJ(5));
t500 = cos(qJ(5));
t439 = -t496 * t442 + t500 * t493;
t440 = t500 * t442 + t496 * t493;
t462 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t496 + Ifges(6,2) * t500) * t483;
t463 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t496 + Ifges(6,4) * t500) * t483;
t518 = qJD(5) * t483;
t467 = t496 * t482 + t500 * t518;
t468 = t500 * t482 - t496 * t518;
t522 = mrSges(6,1) * t439 - mrSges(6,2) * t440 + Ifges(6,5) * t467 + Ifges(6,6) * t468 + Ifges(6,3) * qJDD(5) + (t462 * t496 - t463 * t500) * t483;
t521 = t483 * t496;
t520 = t483 * t500;
t466 = (-mrSges(6,1) * t500 + mrSges(6,2) * t496) * t483;
t473 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t520;
t437 = m(6) * t439 + qJDD(5) * mrSges(6,1) - t467 * mrSges(6,3) + qJD(5) * t473 - t466 * t521;
t472 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t521;
t438 = m(6) * t440 - qJDD(5) * mrSges(6,2) + t468 * mrSges(6,3) - qJD(5) * t472 + t466 * t520;
t512 = -t496 * t437 + t500 * t438;
t423 = m(5) * t445 - t481 * mrSges(5,1) - t482 * mrSges(5,2) + t512;
t444 = t501 * t448 - t497 * t449;
t441 = -t482 * pkin(4) - t481 * pkin(8) - t444;
t508 = -m(6) * t441 + t468 * mrSges(6,1) - t467 * mrSges(6,2) - t472 * t521 + t473 * t520;
t432 = m(5) * t444 + t482 * mrSges(5,1) - t481 * mrSges(5,2) + t508;
t420 = t497 * t423 + t501 * t432;
t416 = m(4) * t451 + t489 * mrSges(4,1) - t488 * mrSges(4,2) + t420;
t513 = t501 * t423 - t497 * t432;
t417 = m(4) * t452 - t488 * mrSges(4,1) - t489 * mrSges(4,2) + t513;
t411 = t502 * t416 + t498 * t417;
t408 = m(3) * t459 + qJDD(1) * mrSges(3,1) - t504 * mrSges(3,2) + t411;
t514 = -t498 * t416 + t502 * t417;
t409 = m(3) * t460 - t504 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t514;
t515 = -t494 * t408 + t495 * t409;
t399 = m(2) * t477 - t504 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t515;
t402 = t495 * t408 + t494 * t409;
t400 = m(2) * t478 + qJDD(1) * mrSges(2,1) - t504 * mrSges(2,2) + t402;
t519 = t499 * t399 + t503 * t400;
t426 = t500 * t437 + t496 * t438;
t517 = m(5) * t493 + t426;
t516 = -t503 * t399 + t499 * t400;
t511 = m(4) * t493 + t517;
t424 = m(3) * t493 + t511;
t461 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t496 + Ifges(6,6) * t500) * t483;
t429 = -mrSges(6,1) * t441 + mrSges(6,3) * t440 + Ifges(6,4) * t467 + Ifges(6,2) * t468 + Ifges(6,6) * qJDD(5) + qJD(5) * t463 - t461 * t521;
t430 = mrSges(6,2) * t441 - mrSges(6,3) * t439 + Ifges(6,1) * t467 + Ifges(6,4) * t468 + Ifges(6,5) * qJDD(5) - qJD(5) * t462 + t461 * t520;
t509 = mrSges(5,1) * t444 - mrSges(5,2) * t445 + Ifges(5,3) * t482 + pkin(4) * t508 + pkin(8) * t512 + t500 * t429 + t496 * t430;
t506 = mrSges(4,1) * t451 - mrSges(4,2) * t452 + Ifges(4,3) * t489 + pkin(3) * t420 + t509;
t505 = mrSges(2,1) * t478 + mrSges(3,1) * t459 - mrSges(2,2) * t477 - mrSges(3,2) * t460 + pkin(1) * t402 + pkin(2) * t411 + t506 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t418 = -mrSges(5,1) * t493 + mrSges(5,3) * t445 + t481 * Ifges(5,5) + Ifges(5,6) * t482 - pkin(4) * t426 - t522;
t412 = mrSges(5,2) * t493 - mrSges(5,3) * t444 + Ifges(5,5) * t482 - t481 * Ifges(5,6) - pkin(8) * t426 - t496 * t429 + t500 * t430;
t404 = mrSges(4,2) * t493 - mrSges(4,3) * t451 + Ifges(4,5) * t489 - t488 * Ifges(4,6) - pkin(7) * t420 + t501 * t412 - t497 * t418;
t403 = -mrSges(4,1) * t493 + mrSges(4,3) * t452 + t488 * Ifges(4,5) + Ifges(4,6) * t489 - pkin(3) * t517 + pkin(7) * t513 + t497 * t412 + t501 * t418;
t395 = mrSges(3,2) * t493 - mrSges(3,3) * t459 + Ifges(3,5) * qJDD(1) - t504 * Ifges(3,6) - pkin(6) * t411 - t498 * t403 + t502 * t404;
t394 = -mrSges(3,1) * t493 + mrSges(3,3) * t460 + t504 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t511 + pkin(6) * t514 + t502 * t403 + t498 * t404;
t393 = -mrSges(2,2) * g(1) - mrSges(2,3) * t478 + Ifges(2,5) * qJDD(1) - t504 * Ifges(2,6) - qJ(2) * t402 - t494 * t394 + t495 * t395;
t392 = mrSges(2,1) * g(1) + mrSges(2,3) * t477 + t504 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t424 + qJ(2) * t515 + t495 * t394 + t494 * t395;
t1 = [(-m(1) - m(2)) * g(1) + t424; -m(1) * g(2) + t519; -m(1) * g(3) + t516; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t505; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t516 + t503 * t392 + t499 * t393; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t519 + t499 * t392 - t503 * t393; t505; t424; t506; t509; t522;];
tauJB = t1;
