% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:47
% EndTime: 2019-12-05 15:06:49
% DurationCPUTime: 1.24s
% Computational Cost: add. (10322->185), mult. (18125->226), div. (0->0), fcn. (10364->8), ass. (0->87)
t547 = Ifges(5,1) + Ifges(6,1);
t538 = Ifges(5,4) + Ifges(6,4);
t537 = Ifges(5,5) + Ifges(6,5);
t546 = Ifges(5,2) + Ifges(6,2);
t536 = Ifges(5,6) + Ifges(6,6);
t545 = Ifges(5,3) + Ifges(6,3);
t510 = sin(qJ(4));
t512 = cos(qJ(4));
t489 = (-mrSges(6,1) * t512 + mrSges(6,2) * t510) * qJD(3);
t527 = qJD(3) * qJD(4);
t522 = t512 * t527;
t491 = t510 * qJDD(3) + t522;
t507 = sin(pkin(7));
t509 = cos(pkin(7));
t496 = -t509 * g(1) - t507 * g(2);
t505 = -g(3) + qJDD(1);
t506 = sin(pkin(8));
t508 = cos(pkin(8));
t473 = -t506 * t496 + t508 * t505;
t474 = t508 * t496 + t506 * t505;
t511 = sin(qJ(3));
t513 = cos(qJ(3));
t469 = t511 * t473 + t513 * t474;
t514 = qJD(3) ^ 2;
t467 = -t514 * pkin(3) + qJDD(3) * pkin(6) + t469;
t495 = t507 * g(1) - t509 * g(2);
t494 = qJDD(2) - t495;
t483 = t512 * t494;
t526 = qJD(3) * qJD(5);
t540 = pkin(4) * t514;
t460 = qJDD(4) * pkin(4) + t483 + (-t491 + t522) * qJ(5) + (t512 * t540 - t467 - 0.2e1 * t526) * t510;
t528 = qJD(3) * t512;
t500 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t528;
t525 = m(6) * t460 + qJDD(4) * mrSges(6,1) + qJD(4) * t500;
t529 = qJD(3) * t510;
t456 = -t491 * mrSges(6,3) - t489 * t529 + t525;
t464 = t512 * t467 + t510 * t494;
t492 = t512 * qJDD(3) - t510 * t527;
t497 = qJD(4) * pkin(4) - qJ(5) * t529;
t504 = t512 ^ 2;
t461 = t492 * qJ(5) - qJD(4) * t497 - t504 * t540 + 0.2e1 * t512 * t526 + t464;
t463 = -t510 * t467 + t483;
t531 = t537 * qJD(4) + (t547 * t510 + t538 * t512) * qJD(3);
t532 = t536 * qJD(4) + (t538 * t510 + t546 * t512) * qJD(3);
t544 = mrSges(5,1) * t463 + mrSges(6,1) * t460 - mrSges(5,2) * t464 - mrSges(6,2) * t461 + pkin(4) * t456 + (t532 * t510 - t531 * t512) * qJD(3) + t545 * qJDD(4) + t537 * t491 + t536 * t492;
t490 = (-mrSges(5,1) * t512 + mrSges(5,2) * t510) * qJD(3);
t501 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t528;
t453 = m(5) * t463 + qJDD(4) * mrSges(5,1) + qJD(4) * t501 + (-mrSges(5,3) - mrSges(6,3)) * t491 + (-t489 - t490) * t529 + t525;
t524 = m(6) * t461 + t492 * mrSges(6,3) + t489 * t528;
t498 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t529;
t530 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t529 - t498;
t539 = -mrSges(5,2) - mrSges(6,2);
t454 = m(5) * t464 + t492 * mrSges(5,3) + t530 * qJD(4) + t539 * qJDD(4) + t490 * t528 + t524;
t446 = t512 * t453 + t510 * t454;
t445 = (m(3) + m(4)) * t494 + t446;
t447 = -t510 * t453 + t512 * t454;
t442 = m(4) * t469 - t514 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t447;
t468 = t513 * t473 - t511 * t474;
t517 = -qJDD(3) * pkin(3) - t468;
t466 = -t514 * pkin(6) + t517;
t462 = t497 * t529 - t492 * pkin(4) + qJDD(5) + (-qJ(5) * t504 - pkin(6)) * t514 + t517;
t518 = -m(6) * t462 + t492 * mrSges(6,1) + t500 * t528;
t455 = -m(5) * t466 + t492 * mrSges(5,1) + t539 * t491 + t501 * t528 + t530 * t529 + t518;
t449 = m(4) * t468 + qJDD(3) * mrSges(4,1) - t514 * mrSges(4,2) + t455;
t437 = t511 * t442 + t513 * t449;
t435 = m(3) * t473 + t437;
t519 = t513 * t442 - t511 * t449;
t436 = m(3) * t474 + t519;
t520 = -t506 * t435 + t508 * t436;
t429 = m(2) * t496 + t520;
t444 = m(2) * t495 - t445;
t534 = t507 * t429 + t509 * t444;
t430 = t508 * t435 + t506 * t436;
t533 = t545 * qJD(4) + (t537 * t510 + t536 * t512) * qJD(3);
t523 = m(2) * t505 + t430;
t521 = t509 * t429 - t507 * t444;
t457 = t491 * mrSges(6,2) + t498 * t529 - t518;
t438 = -mrSges(5,1) * t466 + mrSges(5,3) * t464 - mrSges(6,1) * t462 + mrSges(6,3) * t461 - pkin(4) * t457 + qJ(5) * t524 + t546 * t492 + t538 * t491 + (-qJ(5) * mrSges(6,2) + t536) * qJDD(4) + (-qJ(5) * t498 + t531) * qJD(4) - t533 * t529;
t439 = mrSges(5,2) * t466 + mrSges(6,2) * t462 - mrSges(5,3) * t463 - mrSges(6,3) * t460 - qJ(5) * t456 - t532 * qJD(4) + t537 * qJDD(4) + t547 * t491 + t538 * t492 + t533 * t528;
t515 = mrSges(4,1) * t468 - mrSges(4,2) * t469 + Ifges(4,3) * qJDD(3) + pkin(3) * t455 + pkin(6) * t447 + t512 * t438 + t510 * t439;
t431 = -mrSges(4,1) * t494 + mrSges(4,3) * t469 + t514 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t446 - t544;
t426 = mrSges(4,2) * t494 - mrSges(4,3) * t468 + Ifges(4,5) * qJDD(3) - t514 * Ifges(4,6) - pkin(6) * t446 - t510 * t438 + t512 * t439;
t425 = mrSges(3,2) * t494 - mrSges(3,3) * t473 - pkin(5) * t437 + t513 * t426 - t511 * t431;
t424 = -mrSges(3,1) * t494 + mrSges(3,3) * t474 + t511 * t426 + t513 * t431 - pkin(2) * (m(4) * t494 + t446) + pkin(5) * t519;
t423 = -mrSges(2,1) * t505 - mrSges(3,1) * t473 + mrSges(3,2) * t474 + mrSges(2,3) * t496 - pkin(1) * t430 - pkin(2) * t437 - t515;
t422 = mrSges(2,2) * t505 - mrSges(2,3) * t495 - qJ(2) * t430 - t506 * t424 + t508 * t425;
t1 = [-m(1) * g(1) + t521; -m(1) * g(2) + t534; -m(1) * g(3) + t523; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t534 + t509 * t422 - t507 * t423; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t521 + t507 * t422 + t509 * t423; -mrSges(1,1) * g(2) + mrSges(2,1) * t495 + mrSges(1,2) * g(1) - mrSges(2,2) * t496 - pkin(1) * t445 + qJ(2) * t520 + t508 * t424 + t506 * t425; t523; t445; t515; t544; t457;];
tauJB = t1;
