% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP2
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:18
% EndTime: 2020-01-03 12:11:21
% DurationCPUTime: 2.61s
% Computational Cost: add. (37574->251), mult. (47662->304), div. (0->0), fcn. (28134->8), ass. (0->100)
t547 = Ifges(5,1) + Ifges(6,1);
t543 = Ifges(5,4) + Ifges(6,4);
t542 = Ifges(5,5) + Ifges(6,5);
t546 = Ifges(5,2) + Ifges(6,2);
t541 = -Ifges(5,6) - Ifges(6,6);
t545 = -Ifges(5,3) - Ifges(6,3);
t512 = qJD(1) + qJD(2);
t508 = t512 ^ 2;
t544 = pkin(3) * t508;
t515 = sin(qJ(3));
t540 = t512 * t515;
t519 = cos(qJ(3));
t539 = t512 * t519;
t517 = sin(qJ(1));
t521 = cos(qJ(1));
t505 = -t521 * g(2) - t517 * g(3);
t499 = qJDD(1) * pkin(1) + t505;
t504 = -t517 * g(2) + t521 * g(3);
t522 = qJD(1) ^ 2;
t500 = -t522 * pkin(1) + t504;
t516 = sin(qJ(2));
t520 = cos(qJ(2));
t477 = t516 * t499 + t520 * t500;
t510 = qJDD(1) + qJDD(2);
t475 = -t508 * pkin(2) + t510 * pkin(7) + t477;
t538 = t515 * t475;
t533 = qJD(3) * t512;
t494 = t515 * t510 + t519 * t533;
t450 = qJDD(3) * pkin(3) - t494 * pkin(8) - t538 + (pkin(8) * t533 + t515 * t544 - g(1)) * t519;
t461 = -t515 * g(1) + t519 * t475;
t495 = t519 * t510 - t515 * t533;
t503 = qJD(3) * pkin(3) - pkin(8) * t540;
t513 = t519 ^ 2;
t451 = t495 * pkin(8) - qJD(3) * t503 - t513 * t544 + t461;
t514 = sin(qJ(4));
t518 = cos(qJ(4));
t445 = t518 * t450 - t514 * t451;
t488 = (-t514 * t515 + t518 * t519) * t512;
t459 = t488 * qJD(4) + t518 * t494 + t514 * t495;
t489 = (t514 * t519 + t515 * t518) * t512;
t472 = -t488 * mrSges(6,1) + t489 * mrSges(6,2);
t473 = -t488 * mrSges(5,1) + t489 * mrSges(5,2);
t511 = qJD(3) + qJD(4);
t480 = -t511 * mrSges(5,2) + t488 * mrSges(5,3);
t509 = qJDD(3) + qJDD(4);
t440 = -0.2e1 * qJD(5) * t489 + (t488 * t511 - t459) * qJ(5) + (t488 * t489 + t509) * pkin(4) + t445;
t479 = -t511 * mrSges(6,2) + t488 * mrSges(6,3);
t532 = m(6) * t440 + t509 * mrSges(6,1) + t511 * t479;
t434 = m(5) * t445 + t509 * mrSges(5,1) + t511 * t480 + (-t472 - t473) * t489 + (-mrSges(5,3) - mrSges(6,3)) * t459 + t532;
t446 = t514 * t450 + t518 * t451;
t458 = -t489 * qJD(4) - t514 * t494 + t518 * t495;
t482 = t511 * mrSges(6,1) - t489 * mrSges(6,3);
t483 = t511 * mrSges(5,1) - t489 * mrSges(5,3);
t481 = t511 * pkin(4) - t489 * qJ(5);
t484 = t488 ^ 2;
t442 = -t484 * pkin(4) + t458 * qJ(5) + 0.2e1 * qJD(5) * t488 - t511 * t481 + t446;
t531 = m(6) * t442 + t458 * mrSges(6,3) + t488 * t472;
t437 = m(5) * t446 + t458 * mrSges(5,3) + t488 * t473 + (-t482 - t483) * t511 + (-mrSges(5,2) - mrSges(6,2)) * t509 + t531;
t430 = t518 * t434 + t514 * t437;
t460 = -t519 * g(1) - t538;
t493 = (-mrSges(4,1) * t519 + mrSges(4,2) * t515) * t512;
t502 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t539;
t428 = m(4) * t460 + qJDD(3) * mrSges(4,1) - t494 * mrSges(4,3) + qJD(3) * t502 - t493 * t540 + t430;
t501 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t540;
t527 = -t514 * t434 + t518 * t437;
t429 = m(4) * t461 - qJDD(3) * mrSges(4,2) + t495 * mrSges(4,3) - qJD(3) * t501 + t493 * t539 + t527;
t528 = -t515 * t428 + t519 * t429;
t421 = m(3) * t477 - t508 * mrSges(3,1) - t510 * mrSges(3,2) + t528;
t476 = t520 * t499 - t516 * t500;
t525 = -t510 * pkin(2) - t476;
t474 = -t508 * pkin(7) + t525;
t452 = -t495 * pkin(3) + t503 * t540 + (-pkin(8) * t513 - pkin(7)) * t508 + t525;
t444 = -t458 * pkin(4) - t484 * qJ(5) + t489 * t481 + qJDD(5) + t452;
t526 = m(6) * t444 - t458 * mrSges(6,1) + t459 * mrSges(6,2) - t488 * t479 + t489 * t482;
t524 = m(5) * t452 - t458 * mrSges(5,1) + t459 * mrSges(5,2) - t488 * t480 + t489 * t483 + t526;
t523 = -m(4) * t474 + t495 * mrSges(4,1) - t494 * mrSges(4,2) - t501 * t540 + t502 * t539 - t524;
t432 = m(3) * t476 + t510 * mrSges(3,1) - t508 * mrSges(3,2) + t523;
t529 = t520 * t421 - t516 * t432;
t416 = m(2) * t504 - t522 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t529;
t418 = t516 * t421 + t520 * t432;
t417 = m(2) * t505 + qJDD(1) * mrSges(2,1) - t522 * mrSges(2,2) + t418;
t537 = t517 * t416 + t521 * t417;
t422 = t519 * t428 + t515 * t429;
t536 = t541 * t488 - t542 * t489 + t545 * t511;
t535 = -t546 * t488 - t543 * t489 + t541 * t511;
t534 = t543 * t488 + t547 * t489 + t542 * t511;
t530 = -t521 * t416 + t517 * t417;
t487 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t515 + Ifges(4,4) * t519) * t512;
t486 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t515 + Ifges(4,2) * t519) * t512;
t485 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t515 + Ifges(4,6) * t519) * t512;
t438 = -t459 * mrSges(6,3) - t489 * t472 + t532;
t424 = mrSges(5,2) * t452 + mrSges(6,2) * t444 - mrSges(5,3) * t445 - mrSges(6,3) * t440 - qJ(5) * t438 + t543 * t458 + t547 * t459 - t536 * t488 + t542 * t509 + t535 * t511;
t423 = -mrSges(5,1) * t452 + mrSges(5,3) * t446 - mrSges(6,1) * t444 + mrSges(6,3) * t442 - pkin(4) * t526 + qJ(5) * t531 + (-qJ(5) * t482 + t534) * t511 + (-qJ(5) * mrSges(6,2) - t541) * t509 + t536 * t489 + t543 * t459 + t546 * t458;
t412 = mrSges(4,2) * t474 - mrSges(4,3) * t460 + Ifges(4,1) * t494 + Ifges(4,4) * t495 + Ifges(4,5) * qJDD(3) - pkin(8) * t430 - qJD(3) * t486 - t514 * t423 + t518 * t424 + t485 * t539;
t411 = -mrSges(4,1) * t474 + mrSges(4,3) * t461 + Ifges(4,4) * t494 + Ifges(4,2) * t495 + Ifges(4,6) * qJDD(3) - pkin(3) * t524 + pkin(8) * t527 + qJD(3) * t487 + t518 * t423 + t514 * t424 - t485 * t540;
t410 = mrSges(3,1) * g(1) - Ifges(4,3) * qJDD(3) + Ifges(3,6) * t510 + t508 * Ifges(3,5) - Ifges(4,5) * t494 - Ifges(4,6) * t495 - mrSges(4,1) * t460 + mrSges(4,2) * t461 + mrSges(3,3) * t477 + mrSges(6,2) * t442 - mrSges(5,1) * t445 + mrSges(5,2) * t446 - pkin(4) * t438 - mrSges(6,1) * t440 - pkin(3) * t430 - pkin(2) * t422 + (-t515 * t486 + t519 * t487) * t512 + t545 * t509 + t535 * t489 + t534 * t488 - t542 * t459 + t541 * t458;
t409 = -mrSges(3,2) * g(1) - mrSges(3,3) * t476 + Ifges(3,5) * t510 - t508 * Ifges(3,6) - pkin(7) * t422 - t515 * t411 + t519 * t412;
t408 = -mrSges(2,2) * g(1) - mrSges(2,3) * t505 + Ifges(2,5) * qJDD(1) - t522 * Ifges(2,6) - pkin(6) * t418 + t520 * t409 - t516 * t410;
t407 = Ifges(2,6) * qJDD(1) + t522 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t504 + t516 * t409 + t520 * t410 - pkin(1) * (-m(3) * g(1) + t422) + pkin(6) * t529;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t422; -m(1) * g(2) + t537; -m(1) * g(3) + t530; mrSges(2,1) * t505 + mrSges(3,1) * t476 - mrSges(1,2) * g(3) - mrSges(2,2) * t504 - mrSges(3,2) * t477 + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t510 + pkin(1) * t418 + pkin(2) * t523 + pkin(7) * t528 + t519 * t411 + t515 * t412; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t530 + t521 * t407 + t517 * t408; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t537 + t517 * t407 - t521 * t408;];
tauB = t1;
