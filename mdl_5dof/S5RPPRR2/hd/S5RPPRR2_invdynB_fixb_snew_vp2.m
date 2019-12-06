% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR2
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:48
% EndTime: 2019-12-05 17:39:51
% DurationCPUTime: 2.67s
% Computational Cost: add. (25042->240), mult. (56213->293), div. (0->0), fcn. (37678->8), ass. (0->103)
t506 = sin(qJ(1));
t509 = cos(qJ(1));
t484 = t506 * g(1) - t509 * g(2);
t510 = qJD(1) ^ 2;
t517 = -t510 * qJ(2) + qJDD(2) - t484;
t538 = -pkin(1) - qJ(3);
t544 = -(2 * qJD(1) * qJD(3)) + t538 * qJDD(1) + t517;
t502 = sin(pkin(8));
t495 = t502 ^ 2;
t503 = cos(pkin(8));
t535 = t503 ^ 2 + t495;
t529 = t535 * mrSges(4,3);
t485 = -t509 * g(1) - t506 * g(2);
t543 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t485;
t542 = pkin(3) * t510;
t541 = mrSges(2,1) - mrSges(3,2);
t540 = -Ifges(3,4) + Ifges(2,5);
t539 = -Ifges(2,6) + Ifges(3,5);
t537 = mrSges(4,2) * t503;
t465 = t502 * g(3) + t544 * t503;
t451 = (-pkin(6) * qJDD(1) - t502 * t542) * t503 + t465;
t466 = -t503 * g(3) + t544 * t502;
t531 = qJDD(1) * t502;
t452 = -pkin(6) * t531 - t495 * t542 + t466;
t505 = sin(qJ(4));
t508 = cos(qJ(4));
t441 = t508 * t451 - t505 * t452;
t520 = -t502 * t505 + t503 * t508;
t521 = -t502 * t508 - t503 * t505;
t480 = t521 * qJD(1);
t533 = t480 * qJD(4);
t468 = t520 * qJDD(1) + t533;
t481 = t520 * qJD(1);
t434 = (-t468 + t533) * pkin(7) + (t480 * t481 + qJDD(4)) * pkin(4) + t441;
t442 = t505 * t451 + t508 * t452;
t467 = -t481 * qJD(4) + t521 * qJDD(1);
t475 = qJD(4) * pkin(4) - t481 * pkin(7);
t479 = t480 ^ 2;
t435 = -t479 * pkin(4) + t467 * pkin(7) - qJD(4) * t475 + t442;
t504 = sin(qJ(5));
t507 = cos(qJ(5));
t432 = t507 * t434 - t504 * t435;
t460 = t507 * t480 - t504 * t481;
t440 = t460 * qJD(5) + t504 * t467 + t507 * t468;
t461 = t504 * t480 + t507 * t481;
t447 = -t460 * mrSges(6,1) + t461 * mrSges(6,2);
t497 = qJD(4) + qJD(5);
t453 = -t497 * mrSges(6,2) + t460 * mrSges(6,3);
t494 = qJDD(4) + qJDD(5);
t430 = m(6) * t432 + t494 * mrSges(6,1) - t440 * mrSges(6,3) - t461 * t447 + t497 * t453;
t433 = t504 * t434 + t507 * t435;
t439 = -t461 * qJD(5) + t507 * t467 - t504 * t468;
t454 = t497 * mrSges(6,1) - t461 * mrSges(6,3);
t431 = m(6) * t433 - t494 * mrSges(6,2) + t439 * mrSges(6,3) + t460 * t447 - t497 * t454;
t421 = t507 * t430 + t504 * t431;
t463 = -t480 * mrSges(5,1) + t481 * mrSges(5,2);
t473 = -qJD(4) * mrSges(5,2) + t480 * mrSges(5,3);
t419 = m(5) * t441 + qJDD(4) * mrSges(5,1) - t468 * mrSges(5,3) + qJD(4) * t473 - t481 * t463 + t421;
t474 = qJD(4) * mrSges(5,1) - t481 * mrSges(5,3);
t525 = -t504 * t430 + t507 * t431;
t420 = m(5) * t442 - qJDD(4) * mrSges(5,2) + t467 * mrSges(5,3) - qJD(4) * t474 + t480 * t463 + t525;
t415 = t508 * t419 + t505 * t420;
t519 = -qJDD(1) * mrSges(4,3) - t510 * (mrSges(4,1) * t502 + t537);
t413 = m(4) * t465 + t519 * t503 + t415;
t526 = -t505 * t419 + t508 * t420;
t414 = m(4) * t466 + t519 * t502 + t526;
t409 = t503 * t413 + t502 * t414;
t478 = -qJDD(1) * pkin(1) + t517;
t514 = -m(3) * t478 + t510 * mrSges(3,3) - t409;
t407 = m(2) * t484 - t510 * mrSges(2,2) + t541 * qJDD(1) + t514;
t477 = t510 * pkin(1) - t543;
t516 = qJDD(3) + t543;
t472 = t538 * t510 + t516;
t456 = pkin(3) * t531 + (-t535 * pkin(6) + t538) * t510 + t516;
t437 = -t467 * pkin(4) - t479 * pkin(7) + t481 * t475 + t456;
t515 = m(6) * t437 - t439 * mrSges(6,1) + t440 * mrSges(6,2) - t460 * t453 + t461 * t454;
t513 = m(5) * t456 - t467 * mrSges(5,1) + t468 * mrSges(5,2) - t480 * t473 + t481 * t474 + t515;
t512 = -m(4) * t472 - mrSges(4,1) * t531 - qJDD(1) * t537 - t513;
t511 = -m(3) * t477 + t510 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t512;
t426 = m(2) * t485 + t511 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t529) * t510;
t536 = t509 * t407 + t506 * t426;
t522 = Ifges(4,5) * t503 - Ifges(4,6) * t502;
t534 = t510 * t522;
t528 = -t506 * t407 + t509 * t426;
t527 = -t502 * t413 + t503 * t414;
t524 = Ifges(4,1) * t503 - Ifges(4,4) * t502;
t523 = Ifges(4,4) * t503 - Ifges(4,2) * t502;
t459 = Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * qJD(4);
t458 = Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * qJD(4);
t457 = Ifges(5,5) * t481 + Ifges(5,6) * t480 + Ifges(5,3) * qJD(4);
t445 = Ifges(6,1) * t461 + Ifges(6,4) * t460 + Ifges(6,5) * t497;
t444 = Ifges(6,4) * t461 + Ifges(6,2) * t460 + Ifges(6,6) * t497;
t443 = Ifges(6,5) * t461 + Ifges(6,6) * t460 + Ifges(6,3) * t497;
t423 = mrSges(6,2) * t437 - mrSges(6,3) * t432 + Ifges(6,1) * t440 + Ifges(6,4) * t439 + Ifges(6,5) * t494 + t460 * t443 - t497 * t444;
t422 = -mrSges(6,1) * t437 + mrSges(6,3) * t433 + Ifges(6,4) * t440 + Ifges(6,2) * t439 + Ifges(6,6) * t494 - t461 * t443 + t497 * t445;
t411 = mrSges(5,2) * t456 - mrSges(5,3) * t441 + Ifges(5,1) * t468 + Ifges(5,4) * t467 + Ifges(5,5) * qJDD(4) - pkin(7) * t421 - qJD(4) * t458 - t504 * t422 + t507 * t423 + t480 * t457;
t410 = -mrSges(5,1) * t456 + mrSges(5,3) * t442 + Ifges(5,4) * t468 + Ifges(5,2) * t467 + Ifges(5,6) * qJDD(4) - pkin(4) * t515 + pkin(7) * t525 + qJD(4) * t459 + t507 * t422 + t504 * t423 - t481 * t457;
t408 = -m(3) * g(3) + t527;
t405 = mrSges(4,2) * t472 - mrSges(4,3) * t465 - pkin(6) * t415 + t524 * qJDD(1) - t505 * t410 + t508 * t411 - t502 * t534;
t404 = -mrSges(4,1) * t472 + mrSges(4,3) * t466 - pkin(3) * t513 + pkin(6) * t526 + t523 * qJDD(1) + t508 * t410 + t505 * t411 - t503 * t534;
t403 = Ifges(6,3) * t494 + mrSges(3,1) * t478 - t480 * t459 + t481 * t458 - mrSges(2,3) * t484 + mrSges(4,1) * t465 - mrSges(4,2) * t466 + Ifges(5,6) * t467 + Ifges(5,5) * t468 - t460 * t445 + t461 * t444 + Ifges(6,6) * t439 + Ifges(6,5) * t440 + mrSges(5,1) * t441 - mrSges(5,2) * t442 - mrSges(6,2) * t433 + mrSges(6,1) * t432 + pkin(4) * t421 + pkin(3) * t415 + pkin(2) * t409 - qJ(2) * t408 + (t502 * t524 + t503 * t523 + t539) * t510 + (t522 + t540) * qJDD(1) + Ifges(5,3) * qJDD(4) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t402 = mrSges(2,3) * t485 - mrSges(3,1) * t477 - t502 * t405 - t503 * t404 - pkin(2) * t512 - qJ(3) * t527 - pkin(1) * t408 - t539 * qJDD(1) + t541 * g(3) + (-pkin(2) * t529 + t540) * t510;
t1 = [-m(1) * g(1) + t528; -m(1) * g(2) + t536; (-m(1) - m(2) - m(3)) * g(3) + t527; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t536 - t506 * t402 + t509 * t403; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t528 + t509 * t402 + t506 * t403; pkin(1) * t514 + qJ(2) * (-t510 * t529 + t511) + t503 * t405 - t502 * t404 - qJ(3) * t409 + mrSges(2,1) * t484 - mrSges(2,2) * t485 + mrSges(3,2) * t478 - mrSges(3,3) * t477 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
