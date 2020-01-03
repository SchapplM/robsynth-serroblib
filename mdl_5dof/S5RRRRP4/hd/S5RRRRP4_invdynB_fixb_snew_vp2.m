% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:52
% EndTime: 2019-12-31 21:50:55
% DurationCPUTime: 2.51s
% Computational Cost: add. (35580->249), mult. (44947->303), div. (0->0), fcn. (26206->8), ass. (0->100)
t542 = Ifges(5,1) + Ifges(6,1);
t536 = Ifges(5,4) - Ifges(6,5);
t535 = Ifges(6,4) + Ifges(5,5);
t541 = Ifges(5,2) + Ifges(6,3);
t540 = -Ifges(6,2) - Ifges(5,3);
t534 = Ifges(5,6) - Ifges(6,6);
t539 = cos(qJ(4));
t506 = qJD(1) + qJD(2);
t502 = t506 ^ 2;
t538 = pkin(3) * t502;
t537 = -mrSges(5,3) - mrSges(6,2);
t509 = sin(qJ(3));
t533 = t506 * t509;
t512 = cos(qJ(3));
t532 = t506 * t512;
t511 = sin(qJ(1));
t514 = cos(qJ(1));
t496 = t511 * g(1) - t514 * g(2);
t490 = qJDD(1) * pkin(1) + t496;
t497 = -t514 * g(1) - t511 * g(2);
t515 = qJD(1) ^ 2;
t491 = -t515 * pkin(1) + t497;
t510 = sin(qJ(2));
t513 = cos(qJ(2));
t468 = t510 * t490 + t513 * t491;
t504 = qJDD(1) + qJDD(2);
t466 = -t502 * pkin(2) + t504 * pkin(7) + t468;
t531 = t509 * t466;
t525 = qJD(3) * t506;
t485 = t509 * t504 + t512 * t525;
t442 = qJDD(3) * pkin(3) - t485 * pkin(8) - t531 + (pkin(8) * t525 + t509 * t538 - g(3)) * t512;
t451 = -t509 * g(3) + t512 * t466;
t486 = t512 * t504 - t509 * t525;
t494 = qJD(3) * pkin(3) - pkin(8) * t533;
t507 = t512 ^ 2;
t443 = t486 * pkin(8) - qJD(3) * t494 - t507 * t538 + t451;
t508 = sin(qJ(4));
t439 = t508 * t442 + t539 * t443;
t479 = (t508 * t512 + t539 * t509) * t506;
t448 = t479 * qJD(4) + t508 * t485 - t539 * t486;
t505 = qJD(3) + qJD(4);
t472 = t505 * mrSges(5,1) - t479 * mrSges(5,3);
t478 = t508 * t533 - t539 * t532;
t503 = qJDD(3) + qJDD(4);
t462 = t478 * pkin(4) - t479 * qJ(5);
t501 = t505 ^ 2;
t436 = -t501 * pkin(4) + t503 * qJ(5) + 0.2e1 * qJD(5) * t505 - t478 * t462 + t439;
t473 = -t505 * mrSges(6,1) + t479 * mrSges(6,2);
t524 = m(6) * t436 + t503 * mrSges(6,3) + t505 * t473;
t463 = t478 * mrSges(6,1) - t479 * mrSges(6,3);
t526 = -t478 * mrSges(5,1) - t479 * mrSges(5,2) - t463;
t429 = m(5) * t439 - t503 * mrSges(5,2) + t537 * t448 - t505 * t472 + t526 * t478 + t524;
t438 = t539 * t442 - t508 * t443;
t449 = -t478 * qJD(4) + t539 * t485 + t508 * t486;
t471 = -t505 * mrSges(5,2) - t478 * mrSges(5,3);
t437 = -t503 * pkin(4) - t501 * qJ(5) + t479 * t462 + qJDD(5) - t438;
t474 = -t478 * mrSges(6,2) + t505 * mrSges(6,3);
t519 = -m(6) * t437 + t503 * mrSges(6,1) + t505 * t474;
t431 = m(5) * t438 + t503 * mrSges(5,1) + t537 * t449 + t505 * t471 + t526 * t479 + t519;
t424 = t508 * t429 + t539 * t431;
t450 = -t512 * g(3) - t531;
t484 = (-mrSges(4,1) * t512 + mrSges(4,2) * t509) * t506;
t493 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t532;
t422 = m(4) * t450 + qJDD(3) * mrSges(4,1) - t485 * mrSges(4,3) + qJD(3) * t493 - t484 * t533 + t424;
t492 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t533;
t520 = t539 * t429 - t508 * t431;
t423 = m(4) * t451 - qJDD(3) * mrSges(4,2) + t486 * mrSges(4,3) - qJD(3) * t492 + t484 * t532 + t520;
t521 = -t509 * t422 + t512 * t423;
t415 = m(3) * t468 - t502 * mrSges(3,1) - t504 * mrSges(3,2) + t521;
t467 = t513 * t490 - t510 * t491;
t518 = -t504 * pkin(2) - t467;
t465 = -t502 * pkin(7) + t518;
t444 = -t486 * pkin(3) + t494 * t533 + (-pkin(8) * t507 - pkin(7)) * t502 + t518;
t434 = -0.2e1 * qJD(5) * t479 + (t478 * t505 - t449) * qJ(5) + (t479 * t505 + t448) * pkin(4) + t444;
t432 = m(6) * t434 + t448 * mrSges(6,1) - t449 * mrSges(6,3) - t479 * t473 + t478 * t474;
t517 = m(5) * t444 + t448 * mrSges(5,1) + t449 * mrSges(5,2) + t478 * t471 + t479 * t472 + t432;
t516 = -m(4) * t465 + t486 * mrSges(4,1) - t485 * mrSges(4,2) - t492 * t533 + t493 * t532 - t517;
t426 = m(3) * t467 + t504 * mrSges(3,1) - t502 * mrSges(3,2) + t516;
t412 = t510 * t415 + t513 * t426;
t410 = m(2) * t496 + qJDD(1) * mrSges(2,1) - t515 * mrSges(2,2) + t412;
t522 = t513 * t415 - t510 * t426;
t411 = m(2) * t497 - t515 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t522;
t530 = t514 * t410 + t511 * t411;
t416 = t512 * t422 + t509 * t423;
t529 = t541 * t478 - t536 * t479 - t534 * t505;
t528 = t534 * t478 - t535 * t479 + t540 * t505;
t527 = -t536 * t478 + t542 * t479 + t535 * t505;
t523 = -t511 * t410 + t514 * t411;
t477 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t509 + Ifges(4,4) * t512) * t506;
t476 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t509 + Ifges(4,2) * t512) * t506;
t475 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t509 + Ifges(4,6) * t512) * t506;
t418 = mrSges(5,2) * t444 + mrSges(6,2) * t437 - mrSges(5,3) * t438 - mrSges(6,3) * t434 - qJ(5) * t432 - t536 * t448 + t542 * t449 + t528 * t478 + t535 * t503 + t529 * t505;
t417 = -mrSges(5,1) * t444 - mrSges(6,1) * t434 + mrSges(6,2) * t436 + mrSges(5,3) * t439 - pkin(4) * t432 - t541 * t448 + t536 * t449 + t528 * t479 + t534 * t503 + t527 * t505;
t406 = mrSges(4,2) * t465 - mrSges(4,3) * t450 + Ifges(4,1) * t485 + Ifges(4,4) * t486 + Ifges(4,5) * qJDD(3) - pkin(8) * t424 - qJD(3) * t476 - t508 * t417 + t539 * t418 + t475 * t532;
t405 = -mrSges(4,1) * t465 + mrSges(4,3) * t451 + Ifges(4,4) * t485 + Ifges(4,2) * t486 + Ifges(4,6) * qJDD(3) - pkin(3) * t517 + pkin(8) * t520 + qJD(3) * t477 + t539 * t417 + t508 * t418 - t475 * t533;
t404 = mrSges(3,1) * g(3) - Ifges(4,3) * qJDD(3) - qJ(5) * t524 - pkin(4) * t519 + t502 * Ifges(3,5) + Ifges(3,6) * t504 - Ifges(4,5) * t485 - Ifges(4,6) * t486 + mrSges(3,3) * t468 - mrSges(4,1) * t450 + mrSges(4,2) * t451 - mrSges(6,3) * t436 + mrSges(6,1) * t437 - mrSges(5,1) * t438 + mrSges(5,2) * t439 - pkin(3) * t424 - pkin(2) * t416 + (-t509 * t476 + t512 * t477) * t506 + t540 * t503 + (pkin(4) * t463 + t529) * t479 + (qJ(5) * t463 - t527) * t478 + (pkin(4) * mrSges(6,2) - t535) * t449 + (qJ(5) * mrSges(6,2) + t534) * t448;
t403 = -mrSges(3,2) * g(3) - mrSges(3,3) * t467 + Ifges(3,5) * t504 - t502 * Ifges(3,6) - pkin(7) * t416 - t509 * t405 + t512 * t406;
t402 = -mrSges(2,2) * g(3) - mrSges(2,3) * t496 + Ifges(2,5) * qJDD(1) - t515 * Ifges(2,6) - pkin(6) * t412 + t513 * t403 - t510 * t404;
t401 = Ifges(2,6) * qJDD(1) + t515 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t497 + t510 * t403 + t513 * t404 - pkin(1) * (-m(3) * g(3) + t416) + pkin(6) * t522;
t1 = [-m(1) * g(1) + t523; -m(1) * g(2) + t530; (-m(1) - m(2) - m(3)) * g(3) + t416; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t530 - t511 * t401 + t514 * t402; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t523 + t514 * t401 + t511 * t402; -mrSges(1,1) * g(2) + mrSges(2,1) * t496 + mrSges(3,1) * t467 + mrSges(1,2) * g(1) - mrSges(2,2) * t497 - mrSges(3,2) * t468 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t504 + pkin(1) * t412 + pkin(2) * t516 + pkin(7) * t521 + t512 * t405 + t509 * t406;];
tauB = t1;
