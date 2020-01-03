% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR9
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:28
% EndTime: 2019-12-31 18:02:29
% DurationCPUTime: 1.50s
% Computational Cost: add. (16851->224), mult. (29440->270), div. (0->0), fcn. (13406->8), ass. (0->95)
t515 = qJD(1) ^ 2;
t510 = sin(qJ(1));
t513 = cos(qJ(1));
t492 = -t513 * g(1) - t510 * g(2);
t520 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t492;
t535 = -pkin(1) - pkin(2);
t470 = t515 * t535 + t520;
t491 = t510 * g(1) - t513 * g(2);
t519 = -t515 * qJ(2) + qJDD(2) - t491;
t473 = qJDD(1) * t535 + t519;
t506 = sin(pkin(8));
t507 = cos(pkin(8));
t455 = -t506 * t470 + t507 * t473;
t452 = qJDD(1) * pkin(3) - t515 * pkin(6) - t455;
t509 = sin(qJ(4));
t512 = cos(qJ(4));
t527 = qJD(1) * qJD(4);
t525 = t512 * t527;
t486 = -t509 * qJDD(1) - t525;
t526 = t509 * t527;
t487 = -t512 * qJDD(1) + t526;
t446 = (-t486 + t525) * pkin(7) + (-t487 - t526) * pkin(4) + t452;
t456 = t507 * t470 + t506 * t473;
t453 = -t515 * pkin(3) - qJDD(1) * pkin(6) + t456;
t503 = g(3) + qJDD(3);
t450 = t512 * t453 + t509 * t503;
t485 = (pkin(4) * t512 + pkin(7) * t509) * qJD(1);
t514 = qJD(4) ^ 2;
t528 = t512 * qJD(1);
t448 = -t514 * pkin(4) + qJDD(4) * pkin(7) - t485 * t528 + t450;
t508 = sin(qJ(5));
t511 = cos(qJ(5));
t444 = t511 * t446 - t508 * t448;
t529 = qJD(1) * t509;
t482 = t511 * qJD(4) + t508 * t529;
t463 = t482 * qJD(5) + t508 * qJDD(4) + t511 * t486;
t483 = t508 * qJD(4) - t511 * t529;
t464 = -t482 * mrSges(6,1) + t483 * mrSges(6,2);
t493 = qJD(5) + t528;
t468 = -t493 * mrSges(6,2) + t482 * mrSges(6,3);
t481 = qJDD(5) - t487;
t441 = m(6) * t444 + t481 * mrSges(6,1) - t463 * mrSges(6,3) - t483 * t464 + t493 * t468;
t445 = t508 * t446 + t511 * t448;
t462 = -t483 * qJD(5) + t511 * qJDD(4) - t508 * t486;
t469 = t493 * mrSges(6,1) - t483 * mrSges(6,3);
t442 = m(6) * t445 - t481 * mrSges(6,2) + t462 * mrSges(6,3) + t482 * t464 - t493 * t469;
t436 = -t508 * t441 + t511 * t442;
t531 = t512 * t503;
t447 = -qJDD(4) * pkin(4) - t514 * pkin(7) - t531 + (-qJD(1) * t485 + t453) * t509;
t457 = Ifges(6,5) * t483 + Ifges(6,6) * t482 + Ifges(6,3) * t493;
t459 = Ifges(6,1) * t483 + Ifges(6,4) * t482 + Ifges(6,5) * t493;
t437 = -mrSges(6,1) * t447 + mrSges(6,3) * t445 + Ifges(6,4) * t463 + Ifges(6,2) * t462 + Ifges(6,6) * t481 - t483 * t457 + t493 * t459;
t458 = Ifges(6,4) * t483 + Ifges(6,2) * t482 + Ifges(6,6) * t493;
t438 = mrSges(6,2) * t447 - mrSges(6,3) * t444 + Ifges(6,1) * t463 + Ifges(6,4) * t462 + Ifges(6,5) * t481 + t482 * t457 - t493 * t458;
t443 = -m(6) * t447 + t462 * mrSges(6,1) - t463 * mrSges(6,2) + t482 * t468 - t483 * t469;
t449 = -t509 * t453 + t531;
t477 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t509 - Ifges(5,2) * t512) * qJD(1);
t478 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t509 - Ifges(5,4) * t512) * qJD(1);
t536 = mrSges(5,1) * t449 - mrSges(5,2) * t450 + Ifges(5,5) * t486 + Ifges(5,6) * t487 + Ifges(5,3) * qJDD(4) + pkin(4) * t443 + pkin(7) * t436 - (t509 * t477 - t512 * t478) * qJD(1) + t511 * t437 + t508 * t438;
t534 = mrSges(2,1) + mrSges(3,1);
t533 = Ifges(3,4) + Ifges(2,5);
t532 = Ifges(2,6) - Ifges(3,6);
t474 = -t515 * pkin(1) + t520;
t484 = (mrSges(5,1) * t512 - mrSges(5,2) * t509) * qJD(1);
t489 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t529;
t433 = m(5) * t450 - qJDD(4) * mrSges(5,2) + t487 * mrSges(5,3) - qJD(4) * t489 - t484 * t528 + t436;
t490 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t528;
t439 = m(5) * t449 + qJDD(4) * mrSges(5,1) - t486 * mrSges(5,3) + qJD(4) * t490 + t484 * t529 + t443;
t430 = t512 * t433 - t509 * t439;
t426 = m(4) * t456 - t515 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t430;
t435 = t511 * t441 + t508 * t442;
t434 = -m(5) * t452 + t487 * mrSges(5,1) - t486 * mrSges(5,2) + t489 * t529 - t490 * t528 - t435;
t431 = m(4) * t455 - qJDD(1) * mrSges(4,1) - t515 * mrSges(4,2) + t434;
t523 = t507 * t426 - t506 * t431;
t521 = m(3) * t474 + qJDD(1) * mrSges(3,3) + t523;
t417 = m(2) * t492 - qJDD(1) * mrSges(2,2) - t515 * t534 + t521;
t422 = t506 * t426 + t507 * t431;
t475 = -qJDD(1) * pkin(1) + t519;
t421 = m(3) * t475 - qJDD(1) * mrSges(3,1) - t515 * mrSges(3,3) + t422;
t418 = m(2) * t491 + qJDD(1) * mrSges(2,1) - t515 * mrSges(2,2) - t421;
t530 = t510 * t417 + t513 * t418;
t524 = t513 * t417 - t510 * t418;
t429 = t509 * t433 + t512 * t439;
t428 = m(4) * t503 + t429;
t518 = mrSges(6,1) * t444 - mrSges(6,2) * t445 + Ifges(6,5) * t463 + Ifges(6,6) * t462 + Ifges(6,3) * t481 + t483 * t458 - t482 * t459;
t476 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t509 - Ifges(5,6) * t512) * qJD(1);
t423 = mrSges(5,2) * t452 - mrSges(5,3) * t449 + Ifges(5,1) * t486 + Ifges(5,4) * t487 + Ifges(5,5) * qJDD(4) - pkin(7) * t435 - qJD(4) * t477 - t508 * t437 + t511 * t438 - t476 * t528;
t424 = -mrSges(5,1) * t452 + mrSges(5,3) * t450 + Ifges(5,4) * t486 + Ifges(5,2) * t487 + Ifges(5,6) * qJDD(4) - pkin(4) * t435 + qJD(4) * t478 + t476 * t529 - t518;
t516 = -mrSges(3,1) * t475 - mrSges(4,1) * t455 - mrSges(2,2) * t492 - pkin(2) * t422 - pkin(3) * t434 - pkin(6) * t430 - t509 * t423 - t512 * t424 + qJ(2) * (-t515 * mrSges(3,1) + t521) - pkin(1) * t421 + mrSges(4,2) * t456 + mrSges(3,3) * t474 + mrSges(2,1) * t491 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t427 = -m(3) * g(3) - t428;
t413 = -mrSges(4,1) * t503 + mrSges(4,3) * t456 + t515 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t429 - t536;
t412 = mrSges(4,2) * t503 - mrSges(4,3) * t455 - Ifges(4,5) * qJDD(1) - t515 * Ifges(4,6) - pkin(6) * t429 + t512 * t423 - t509 * t424;
t411 = mrSges(3,2) * t475 - mrSges(2,3) * t491 - qJ(2) * t427 - qJ(3) * t422 + t507 * t412 - t506 * t413 - t532 * t515 + t533 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t410 = mrSges(3,2) * t474 + mrSges(2,3) * t492 - pkin(1) * t427 + pkin(2) * t428 + g(3) * t534 - qJ(3) * t523 + qJDD(1) * t532 - t506 * t412 - t507 * t413 + t515 * t533;
t1 = [-m(1) * g(1) + t524; -m(1) * g(2) + t530; (-m(1) - m(2) - m(3)) * g(3) - t428; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t530 - t510 * t410 + t513 * t411; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t524 + t513 * t410 + t510 * t411; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t516; t516; t421; t428; t536; t518;];
tauJB = t1;
