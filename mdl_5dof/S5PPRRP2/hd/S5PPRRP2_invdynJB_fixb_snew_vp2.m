% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:42
% EndTime: 2019-12-05 15:08:44
% DurationCPUTime: 1.20s
% Computational Cost: add. (10168->182), mult. (17503->228), div. (0->0), fcn. (10040->8), ass. (0->82)
t527 = Ifges(5,1) + Ifges(6,1);
t519 = Ifges(5,4) - Ifges(6,5);
t518 = -Ifges(5,5) - Ifges(6,4);
t526 = Ifges(5,2) + Ifges(6,3);
t517 = Ifges(5,6) - Ifges(6,6);
t525 = Ifges(5,3) + Ifges(6,2);
t494 = sin(qJ(4));
t496 = cos(qJ(4));
t472 = (-mrSges(6,1) * t496 - mrSges(6,3) * t494) * qJD(3);
t508 = qJD(3) * qJD(4);
t474 = t494 * qJDD(3) + t496 * t508;
t491 = sin(pkin(7));
t493 = cos(pkin(7));
t480 = -t493 * g(1) - t491 * g(2);
t489 = -g(3) + qJDD(1);
t490 = sin(pkin(8));
t492 = cos(pkin(8));
t456 = -t490 * t480 + t492 * t489;
t457 = t492 * t480 + t490 * t489;
t495 = sin(qJ(3));
t497 = cos(qJ(3));
t452 = t495 * t456 + t497 * t457;
t499 = qJD(3) ^ 2;
t450 = -t499 * pkin(3) + qJDD(3) * pkin(6) + t452;
t471 = (-pkin(4) * t496 - qJ(5) * t494) * qJD(3);
t498 = qJD(4) ^ 2;
t479 = t491 * g(1) - t493 * g(2);
t478 = qJDD(2) - t479;
t515 = t496 * t478;
t445 = -qJDD(4) * pkin(4) - t498 * qJ(5) - t515 + qJDD(5) + (qJD(3) * t471 + t450) * t494;
t509 = qJD(3) * t496;
t484 = mrSges(6,2) * t509 + qJD(4) * mrSges(6,3);
t502 = -m(6) * t445 + qJDD(4) * mrSges(6,1) + qJD(4) * t484;
t510 = qJD(3) * t494;
t441 = t474 * mrSges(6,2) + t472 * t510 - t502;
t447 = t496 * t450 + t494 * t478;
t444 = -t498 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t471 * t509 + t447;
t446 = -t494 * t450 + t515;
t475 = t496 * qJDD(3) - t494 * t508;
t482 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t510;
t503 = m(6) * t444 + qJDD(4) * mrSges(6,3) + qJD(4) * t482 + t472 * t509;
t511 = -t518 * qJD(4) + (t527 * t494 + t519 * t496) * qJD(3);
t513 = -t517 * qJD(4) + (-t519 * t494 - t526 * t496) * qJD(3);
t524 = -(t513 * t494 + t511 * t496) * qJD(3) + t525 * qJDD(4) - t518 * t474 + t517 * t475 + mrSges(5,1) * t446 - mrSges(6,1) * t445 - mrSges(5,2) * t447 + mrSges(6,3) * t444 - pkin(4) * t441 + qJ(5) * (t475 * mrSges(6,2) + t503);
t473 = (-mrSges(5,1) * t496 + mrSges(5,2) * t494) * qJD(3);
t481 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t510;
t520 = mrSges(5,3) + mrSges(6,2);
t437 = m(5) * t447 - qJDD(4) * mrSges(5,2) - qJD(4) * t481 + t473 * t509 + t520 * t475 + t503;
t483 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t509;
t438 = m(5) * t446 + qJDD(4) * mrSges(5,1) + qJD(4) * t483 - t520 * t474 + (-t472 - t473) * t510 + t502;
t429 = t494 * t437 + t496 * t438;
t428 = (m(3) + m(4)) * t478 + t429;
t430 = t496 * t437 - t494 * t438;
t425 = m(4) * t452 - t499 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t430;
t451 = t497 * t456 - t495 * t457;
t449 = -qJDD(3) * pkin(3) - t499 * pkin(6) - t451;
t442 = -t475 * pkin(4) - t474 * qJ(5) + (-0.2e1 * qJD(5) * t494 + (pkin(4) * t494 - qJ(5) * t496) * qJD(4)) * qJD(3) + t449;
t439 = m(6) * t442 - t475 * mrSges(6,1) - t474 * mrSges(6,3) - t482 * t510 - t484 * t509;
t433 = -m(5) * t449 + t475 * mrSges(5,1) - t474 * mrSges(5,2) - t481 * t510 + t483 * t509 - t439;
t432 = m(4) * t451 + qJDD(3) * mrSges(4,1) - t499 * mrSges(4,2) + t433;
t420 = t495 * t425 + t497 * t432;
t418 = m(3) * t456 + t420;
t504 = t497 * t425 - t495 * t432;
t419 = m(3) * t457 + t504;
t505 = -t490 * t418 + t492 * t419;
t411 = m(2) * t480 + t505;
t427 = m(2) * t479 - t428;
t514 = t491 * t411 + t493 * t427;
t412 = t492 * t418 + t490 * t419;
t512 = t525 * qJD(4) + (-t518 * t494 + t517 * t496) * qJD(3);
t507 = m(2) * t489 + t412;
t506 = t493 * t411 - t491 * t427;
t421 = -mrSges(5,1) * t449 - mrSges(6,1) * t442 + mrSges(6,2) * t444 + mrSges(5,3) * t447 - pkin(4) * t439 + t511 * qJD(4) + t517 * qJDD(4) + t519 * t474 + t526 * t475 - t512 * t510;
t422 = mrSges(5,2) * t449 + mrSges(6,2) * t445 - mrSges(5,3) * t446 - mrSges(6,3) * t442 - qJ(5) * t439 + t513 * qJD(4) - t518 * qJDD(4) + t527 * t474 + t519 * t475 + t512 * t509;
t500 = mrSges(4,1) * t451 - mrSges(4,2) * t452 + Ifges(4,3) * qJDD(3) + pkin(3) * t433 + pkin(6) * t430 + t496 * t421 + t494 * t422;
t414 = -mrSges(4,1) * t478 + mrSges(4,3) * t452 + t499 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t429 - t524;
t413 = mrSges(4,2) * t478 - mrSges(4,3) * t451 + Ifges(4,5) * qJDD(3) - t499 * Ifges(4,6) - pkin(6) * t429 - t494 * t421 + t496 * t422;
t408 = mrSges(3,2) * t478 - mrSges(3,3) * t456 - pkin(5) * t420 + t497 * t413 - t495 * t414;
t407 = -mrSges(3,1) * t478 + mrSges(3,3) * t457 + t495 * t413 + t497 * t414 - pkin(2) * (m(4) * t478 + t429) + pkin(5) * t504;
t406 = -mrSges(2,1) * t489 - mrSges(3,1) * t456 + mrSges(3,2) * t457 + mrSges(2,3) * t480 - pkin(1) * t412 - pkin(2) * t420 - t500;
t405 = mrSges(2,2) * t489 - mrSges(2,3) * t479 - qJ(2) * t412 - t490 * t407 + t492 * t408;
t1 = [-m(1) * g(1) + t506; -m(1) * g(2) + t514; -m(1) * g(3) + t507; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t514 + t493 * t405 - t491 * t406; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t506 + t491 * t405 + t493 * t406; -mrSges(1,1) * g(2) + mrSges(2,1) * t479 + mrSges(1,2) * g(1) - mrSges(2,2) * t480 - pkin(1) * t428 + qJ(2) * t505 + t492 * t407 + t490 * t408; t507; t428; t500; t524; t441;];
tauJB = t1;
