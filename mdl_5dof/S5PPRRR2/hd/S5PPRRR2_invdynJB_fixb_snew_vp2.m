% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRR2
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:25
% EndTime: 2019-12-05 15:14:27
% DurationCPUTime: 2.20s
% Computational Cost: add. (24223->207), mult. (43995->269), div. (0->0), fcn. (28948->10), ass. (0->93)
t509 = sin(pkin(8));
t511 = cos(pkin(8));
t498 = -t511 * g(1) - t509 * g(2);
t507 = -g(3) + qJDD(1);
t508 = sin(pkin(9));
t510 = cos(pkin(9));
t482 = -t508 * t498 + t510 * t507;
t483 = t510 * t498 + t508 * t507;
t514 = sin(qJ(3));
t517 = cos(qJ(3));
t469 = t514 * t482 + t517 * t483;
t518 = qJD(3) ^ 2;
t464 = -t518 * pkin(3) + qJDD(3) * pkin(6) + t469;
t497 = t509 * g(1) - t511 * g(2);
t496 = qJDD(2) - t497;
t513 = sin(qJ(4));
t516 = cos(qJ(4));
t459 = -t513 * t464 + t516 * t496;
t531 = qJD(3) * qJD(4);
t529 = t516 * t531;
t494 = t513 * qJDD(3) + t529;
t456 = (-t494 + t529) * pkin(7) + (t513 * t516 * t518 + qJDD(4)) * pkin(4) + t459;
t460 = t516 * t464 + t513 * t496;
t495 = t516 * qJDD(3) - t513 * t531;
t533 = qJD(3) * t513;
t501 = qJD(4) * pkin(4) - pkin(7) * t533;
t506 = t516 ^ 2;
t457 = -t506 * t518 * pkin(4) + t495 * pkin(7) - qJD(4) * t501 + t460;
t512 = sin(qJ(5));
t515 = cos(qJ(5));
t454 = t515 * t456 - t512 * t457;
t487 = (-t512 * t513 + t515 * t516) * qJD(3);
t471 = t487 * qJD(5) + t515 * t494 + t512 * t495;
t488 = (t512 * t516 + t513 * t515) * qJD(3);
t476 = -t487 * mrSges(6,1) + t488 * mrSges(6,2);
t505 = qJD(4) + qJD(5);
t480 = -t505 * mrSges(6,2) + t487 * mrSges(6,3);
t504 = qJDD(4) + qJDD(5);
t451 = m(6) * t454 + t504 * mrSges(6,1) - t471 * mrSges(6,3) - t488 * t476 + t505 * t480;
t455 = t512 * t456 + t515 * t457;
t470 = -t488 * qJD(5) - t512 * t494 + t515 * t495;
t481 = t505 * mrSges(6,1) - t488 * mrSges(6,3);
t452 = m(6) * t455 - t504 * mrSges(6,2) + t470 * mrSges(6,3) + t487 * t476 - t505 * t481;
t442 = t515 * t451 + t512 * t452;
t493 = (-mrSges(5,1) * t516 + mrSges(5,2) * t513) * qJD(3);
t532 = qJD(3) * t516;
t500 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t532;
t440 = m(5) * t459 + qJDD(4) * mrSges(5,1) - t494 * mrSges(5,3) + qJD(4) * t500 - t493 * t533 + t442;
t499 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t533;
t525 = -t512 * t451 + t515 * t452;
t441 = m(5) * t460 - qJDD(4) * mrSges(5,2) + t495 * mrSges(5,3) - qJD(4) * t499 + t493 * t532 + t525;
t435 = t516 * t440 + t513 * t441;
t434 = (m(3) + m(4)) * t496 + t435;
t485 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t513 + Ifges(5,2) * t516) * qJD(3);
t486 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t513 + Ifges(5,4) * t516) * qJD(3);
t473 = Ifges(6,4) * t488 + Ifges(6,2) * t487 + Ifges(6,6) * t505;
t474 = Ifges(6,1) * t488 + Ifges(6,4) * t487 + Ifges(6,5) * t505;
t521 = -mrSges(6,1) * t454 + mrSges(6,2) * t455 - Ifges(6,5) * t471 - Ifges(6,6) * t470 - Ifges(6,3) * t504 - t488 * t473 + t487 * t474;
t536 = mrSges(5,1) * t459 - mrSges(5,2) * t460 + Ifges(5,5) * t494 + Ifges(5,6) * t495 + Ifges(5,3) * qJDD(4) + pkin(4) * t442 + (t513 * t485 - t516 * t486) * qJD(3) - t521;
t436 = -t513 * t440 + t516 * t441;
t431 = m(4) * t469 - t518 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t436;
t468 = t517 * t482 - t514 * t483;
t523 = -qJDD(3) * pkin(3) - t468;
t463 = -t518 * pkin(6) + t523;
t458 = t501 * t533 - t495 * pkin(4) + (-pkin(7) * t506 - pkin(6)) * t518 + t523;
t522 = m(6) * t458 - t470 * mrSges(6,1) + t471 * mrSges(6,2) - t487 * t480 + t488 * t481;
t447 = -m(5) * t463 + t495 * mrSges(5,1) - t494 * mrSges(5,2) - t499 * t533 + t500 * t532 - t522;
t446 = m(4) * t468 + qJDD(3) * mrSges(4,1) - t518 * mrSges(4,2) + t447;
t427 = t514 * t431 + t517 * t446;
t425 = m(3) * t482 + t427;
t526 = t517 * t431 - t514 * t446;
t426 = m(3) * t483 + t526;
t527 = -t508 * t425 + t510 * t426;
t418 = m(2) * t498 + t527;
t433 = m(2) * t497 - t434;
t534 = t509 * t418 + t511 * t433;
t419 = t510 * t425 + t508 * t426;
t530 = m(2) * t507 + t419;
t528 = t511 * t418 - t509 * t433;
t472 = Ifges(6,5) * t488 + Ifges(6,6) * t487 + Ifges(6,3) * t505;
t443 = -mrSges(6,1) * t458 + mrSges(6,3) * t455 + Ifges(6,4) * t471 + Ifges(6,2) * t470 + Ifges(6,6) * t504 - t488 * t472 + t505 * t474;
t444 = mrSges(6,2) * t458 - mrSges(6,3) * t454 + Ifges(6,1) * t471 + Ifges(6,4) * t470 + Ifges(6,5) * t504 + t487 * t472 - t505 * t473;
t484 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t513 + Ifges(5,6) * t516) * qJD(3);
t421 = -mrSges(5,1) * t463 + mrSges(5,3) * t460 + Ifges(5,4) * t494 + Ifges(5,2) * t495 + Ifges(5,6) * qJDD(4) - pkin(4) * t522 + pkin(7) * t525 + qJD(4) * t486 + t515 * t443 + t512 * t444 - t484 * t533;
t428 = mrSges(5,2) * t463 - mrSges(5,3) * t459 + Ifges(5,1) * t494 + Ifges(5,4) * t495 + Ifges(5,5) * qJDD(4) - pkin(7) * t442 - qJD(4) * t485 - t512 * t443 + t515 * t444 + t484 * t532;
t520 = mrSges(4,1) * t468 - mrSges(4,2) * t469 + Ifges(4,3) * qJDD(3) + pkin(3) * t447 + pkin(6) * t436 + t516 * t421 + t513 * t428;
t420 = -mrSges(4,1) * t496 + mrSges(4,3) * t469 + t518 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t435 - t536;
t415 = mrSges(4,2) * t496 - mrSges(4,3) * t468 + Ifges(4,5) * qJDD(3) - t518 * Ifges(4,6) - pkin(6) * t435 - t513 * t421 + t516 * t428;
t414 = mrSges(3,2) * t496 - mrSges(3,3) * t482 - pkin(5) * t427 + t517 * t415 - t514 * t420;
t413 = -mrSges(2,1) * t507 - mrSges(3,1) * t482 + mrSges(3,2) * t483 + mrSges(2,3) * t498 - pkin(1) * t419 - pkin(2) * t427 - t520;
t412 = -mrSges(3,1) * t496 + mrSges(3,3) * t483 + t514 * t415 + t517 * t420 - pkin(2) * (m(4) * t496 + t435) + pkin(5) * t526;
t411 = mrSges(2,2) * t507 - mrSges(2,3) * t497 - qJ(2) * t419 - t508 * t412 + t510 * t414;
t1 = [-m(1) * g(1) + t528; -m(1) * g(2) + t534; -m(1) * g(3) + t530; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t534 + t511 * t411 - t509 * t413; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t528 + t509 * t411 + t511 * t413; -mrSges(1,1) * g(2) + mrSges(2,1) * t497 + mrSges(1,2) * g(1) - mrSges(2,2) * t498 - pkin(1) * t434 + qJ(2) * t527 + t510 * t412 + t508 * t414; t530; t434; t520; t536; -t521;];
tauJB = t1;
