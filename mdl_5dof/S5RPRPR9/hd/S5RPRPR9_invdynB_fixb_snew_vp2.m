% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:52
% EndTime: 2019-12-31 18:23:55
% DurationCPUTime: 1.91s
% Computational Cost: add. (15491->255), mult. (30073->303), div. (0->0), fcn. (15118->8), ass. (0->108)
t547 = Ifges(4,1) + Ifges(5,2);
t540 = Ifges(4,4) + Ifges(5,6);
t539 = Ifges(4,5) - Ifges(5,4);
t546 = Ifges(4,2) + Ifges(5,3);
t538 = Ifges(4,6) - Ifges(5,5);
t545 = Ifges(4,3) + Ifges(5,1);
t544 = -2 * qJD(4);
t543 = -pkin(3) - pkin(7);
t513 = qJD(1) ^ 2;
t542 = pkin(7) * t513;
t541 = t513 * pkin(6);
t503 = -g(3) + qJDD(2);
t510 = cos(qJ(3));
t537 = t510 * t503;
t508 = sin(qJ(1));
t511 = cos(qJ(1));
t492 = t508 * g(1) - t511 * g(2);
t479 = qJDD(1) * pkin(1) + t492;
t493 = -t511 * g(1) - t508 * g(2);
t483 = -t513 * pkin(1) + t493;
t504 = sin(pkin(8));
t505 = cos(pkin(8));
t456 = t504 * t479 + t505 * t483;
t447 = -t513 * pkin(2) + qJDD(1) * pkin(6) + t456;
t507 = sin(qJ(3));
t444 = t507 * t447;
t442 = -t444 + t537;
t481 = (mrSges(5,2) * t510 - mrSges(5,3) * t507) * qJD(1);
t482 = (-mrSges(4,1) * t510 + mrSges(4,2) * t507) * qJD(1);
t529 = qJD(1) * qJD(3);
t526 = t510 * t529;
t484 = t507 * qJDD(1) + t526;
t531 = qJD(1) * t510;
t488 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t531;
t489 = -mrSges(5,1) * t531 - qJD(3) * mrSges(5,3);
t527 = t507 * t529;
t485 = t510 * qJDD(1) - t527;
t530 = t507 * qJD(1);
t491 = pkin(4) * t530 - qJD(3) * pkin(7);
t502 = t510 ^ 2;
t455 = t505 * t479 - t504 * t483;
t521 = -qJDD(1) * pkin(2) - t455;
t515 = pkin(3) * t527 + t530 * t544 + (-t484 - t526) * qJ(4) + t521;
t435 = -t491 * t530 + (-pkin(4) * t502 - pkin(6)) * t513 + t543 * t485 + t515;
t480 = (-pkin(3) * t510 - qJ(4) * t507) * qJD(1);
t512 = qJD(3) ^ 2;
t522 = -t512 * qJ(4) + t480 * t530 + qJDD(4) + t444;
t438 = t484 * pkin(4) + t543 * qJDD(3) + (-pkin(4) * t529 - t507 * t542 - t503) * t510 + t522;
t506 = sin(qJ(5));
t509 = cos(qJ(5));
t433 = -t506 * t435 + t509 * t438;
t477 = -t506 * qJD(3) - t509 * t531;
t454 = t477 * qJD(5) + t509 * qJDD(3) - t506 * t485;
t478 = t509 * qJD(3) - t506 * t531;
t457 = -t477 * mrSges(6,1) + t478 * mrSges(6,2);
t495 = qJD(5) + t530;
t458 = -t495 * mrSges(6,2) + t477 * mrSges(6,3);
t476 = qJDD(5) + t484;
t431 = m(6) * t433 + t476 * mrSges(6,1) - t454 * mrSges(6,3) - t478 * t457 + t495 * t458;
t434 = t509 * t435 + t506 * t438;
t453 = -t478 * qJD(5) - t506 * qJDD(3) - t509 * t485;
t459 = t495 * mrSges(6,1) - t478 * mrSges(6,3);
t432 = m(6) * t434 - t476 * mrSges(6,2) + t453 * mrSges(6,3) + t477 * t457 - t495 * t459;
t423 = t509 * t431 + t506 * t432;
t441 = -qJDD(3) * pkin(3) + t522 - t537;
t518 = -m(5) * t441 - t484 * mrSges(5,1) - t423;
t421 = m(4) * t442 - t484 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t488 - t489) * qJD(3) + (-t481 - t482) * t530 + t518;
t443 = t510 * t447 + t507 * t503;
t487 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t530;
t517 = -t512 * pkin(3) + qJDD(3) * qJ(4) + t480 * t531 + t443;
t440 = qJD(3) * t544 - t517;
t490 = mrSges(5,1) * t530 + qJD(3) * mrSges(5,2);
t437 = -t502 * t542 + t485 * pkin(4) + ((2 * qJD(4)) + t491) * qJD(3) + t517;
t519 = -m(6) * t437 + t453 * mrSges(6,1) - t454 * mrSges(6,2) + t477 * t458 - t478 * t459;
t516 = -m(5) * t440 + qJDD(3) * mrSges(5,3) + qJD(3) * t490 + t481 * t531 - t519;
t428 = t482 * t531 + m(4) * t443 - qJDD(3) * mrSges(4,2) - qJD(3) * t487 + (mrSges(4,3) + mrSges(5,1)) * t485 + t516;
t523 = -t507 * t421 + t510 * t428;
t416 = m(3) * t456 - t513 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t523;
t446 = t521 - t541;
t439 = -t485 * pkin(3) + t515 - t541;
t535 = -t506 * t431 + t509 * t432;
t520 = -m(5) * t439 - t485 * mrSges(5,2) + t490 * t530 - t535;
t514 = -m(4) * t446 + t488 * t531 + t485 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t484 + (-t487 * t507 - t489 * t510) * qJD(1) + t520;
t419 = m(3) * t455 + qJDD(1) * mrSges(3,1) - t513 * mrSges(3,2) + t514;
t413 = t504 * t416 + t505 * t419;
t411 = m(2) * t492 + qJDD(1) * mrSges(2,1) - t513 * mrSges(2,2) + t413;
t524 = t505 * t416 - t504 * t419;
t412 = m(2) * t493 - t513 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t524;
t536 = t511 * t411 + t508 * t412;
t417 = t510 * t421 + t507 * t428;
t534 = t545 * qJD(3) + (t539 * t507 + t538 * t510) * qJD(1);
t533 = -t538 * qJD(3) + (-t540 * t507 - t546 * t510) * qJD(1);
t532 = t539 * qJD(3) + (t547 * t507 + t540 * t510) * qJD(1);
t528 = m(3) * t503 + t417;
t525 = -t508 * t411 + t511 * t412;
t450 = Ifges(6,1) * t478 + Ifges(6,4) * t477 + Ifges(6,5) * t495;
t449 = Ifges(6,4) * t478 + Ifges(6,2) * t477 + Ifges(6,6) * t495;
t448 = Ifges(6,5) * t478 + Ifges(6,6) * t477 + Ifges(6,3) * t495;
t425 = mrSges(6,2) * t437 - mrSges(6,3) * t433 + Ifges(6,1) * t454 + Ifges(6,4) * t453 + Ifges(6,5) * t476 + t477 * t448 - t495 * t449;
t424 = -mrSges(6,1) * t437 + mrSges(6,3) * t434 + Ifges(6,4) * t454 + Ifges(6,2) * t453 + Ifges(6,6) * t476 - t478 * t448 + t495 * t450;
t422 = -t484 * mrSges(5,3) + t489 * t531 - t520;
t407 = mrSges(5,1) * t441 + mrSges(6,1) * t433 + mrSges(4,2) * t446 - mrSges(6,2) * t434 - mrSges(4,3) * t442 - mrSges(5,3) * t439 + Ifges(6,5) * t454 + Ifges(6,6) * t453 + Ifges(6,3) * t476 + pkin(4) * t423 - qJ(4) * t422 + t478 * t449 - t477 * t450 + t540 * t485 + t547 * t484 + t539 * qJDD(3) + t533 * qJD(3) + t534 * t531;
t406 = -mrSges(4,1) * t446 - mrSges(5,1) * t440 + mrSges(5,2) * t439 + mrSges(4,3) * t443 - pkin(3) * t422 - pkin(4) * t519 - pkin(7) * t535 + t532 * qJD(3) + t538 * qJDD(3) - t509 * t424 - t506 * t425 + t540 * t484 + t546 * t485 - t534 * t530;
t405 = -pkin(2) * t417 - mrSges(3,1) * t503 + mrSges(3,3) * t456 + t506 * t424 + pkin(7) * t423 - mrSges(4,1) * t442 + mrSges(4,2) * t443 - pkin(3) * (-qJD(3) * t489 + t518) - qJ(4) * t516 - mrSges(5,2) * t441 + mrSges(5,3) * t440 - t509 * t425 + t513 * Ifges(3,5) + Ifges(3,6) * qJDD(1) + (-qJ(4) * mrSges(5,1) - t538) * t485 - t539 * t484 + (pkin(3) * mrSges(5,2) - t545) * qJDD(3) + (t532 * t510 + (pkin(3) * t481 + t533) * t507) * qJD(1);
t404 = mrSges(3,2) * t503 - mrSges(3,3) * t455 + Ifges(3,5) * qJDD(1) - t513 * Ifges(3,6) - pkin(6) * t417 - t507 * t406 + t510 * t407;
t403 = -mrSges(2,2) * g(3) - mrSges(2,3) * t492 + Ifges(2,5) * qJDD(1) - t513 * Ifges(2,6) - qJ(2) * t413 + t505 * t404 - t504 * t405;
t402 = mrSges(2,1) * g(3) + mrSges(2,3) * t493 + t513 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t528 + qJ(2) * t524 + t504 * t404 + t505 * t405;
t1 = [-m(1) * g(1) + t525; -m(1) * g(2) + t536; (-m(1) - m(2)) * g(3) + t528; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t536 - t508 * t402 + t511 * t403; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t525 + t511 * t402 + t508 * t403; pkin(1) * t413 + mrSges(2,1) * t492 - mrSges(2,2) * t493 + pkin(6) * t523 + t507 * t407 + t510 * t406 + pkin(2) * t514 + mrSges(3,1) * t455 - mrSges(3,2) * t456 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
