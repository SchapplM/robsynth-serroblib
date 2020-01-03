% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:59
% EndTime: 2020-01-03 11:28:07
% DurationCPUTime: 4.70s
% Computational Cost: add. (48242->248), mult. (106816->311), div. (0->0), fcn. (71884->10), ass. (0->108)
t515 = qJD(1) ^ 2;
t507 = cos(pkin(9));
t540 = pkin(3) * t507;
t505 = sin(pkin(9));
t539 = mrSges(4,2) * t505;
t502 = t507 ^ 2;
t538 = t502 * t515;
t511 = sin(qJ(1));
t514 = cos(qJ(1));
t489 = -t511 * g(2) + t514 * g(3);
t490 = -t514 * g(2) - t511 * g(3);
t487 = qJDD(1) * pkin(1) + t490;
t488 = -t515 * pkin(1) + t489;
t506 = sin(pkin(8));
t508 = cos(pkin(8));
t475 = t506 * t487 + t508 * t488;
t468 = -t515 * pkin(2) + qJDD(1) * qJ(3) + t475;
t504 = -g(1) + qJDD(2);
t533 = qJD(1) * qJD(3);
t536 = t507 * t504 - 0.2e1 * t505 * t533;
t454 = (-pkin(6) * qJDD(1) + t515 * t540 - t468) * t505 + t536;
t458 = t505 * t504 + (t468 + 0.2e1 * t533) * t507;
t532 = qJDD(1) * t507;
t455 = -pkin(3) * t538 + pkin(6) * t532 + t458;
t510 = sin(qJ(4));
t513 = cos(qJ(4));
t439 = t513 * t454 - t510 * t455;
t520 = t505 * t513 + t507 * t510;
t519 = -t505 * t510 + t507 * t513;
t480 = t519 * qJD(1);
t534 = t480 * qJD(4);
t473 = t520 * qJDD(1) + t534;
t481 = t520 * qJD(1);
t435 = (-t473 + t534) * pkin(7) + (t480 * t481 + qJDD(4)) * pkin(4) + t439;
t440 = t510 * t454 + t513 * t455;
t472 = -t481 * qJD(4) + t519 * qJDD(1);
t478 = qJD(4) * pkin(4) - t481 * pkin(7);
t479 = t480 ^ 2;
t436 = -t479 * pkin(4) + t472 * pkin(7) - qJD(4) * t478 + t440;
t509 = sin(qJ(5));
t512 = cos(qJ(5));
t433 = t512 * t435 - t509 * t436;
t466 = t512 * t480 - t509 * t481;
t444 = t466 * qJD(5) + t509 * t472 + t512 * t473;
t467 = t509 * t480 + t512 * t481;
t450 = -t466 * mrSges(6,1) + t467 * mrSges(6,2);
t503 = qJD(4) + qJD(5);
t459 = -t503 * mrSges(6,2) + t466 * mrSges(6,3);
t500 = qJDD(4) + qJDD(5);
t431 = m(6) * t433 + t500 * mrSges(6,1) - t444 * mrSges(6,3) - t467 * t450 + t503 * t459;
t434 = t509 * t435 + t512 * t436;
t443 = -t467 * qJD(5) + t512 * t472 - t509 * t473;
t460 = t503 * mrSges(6,1) - t467 * mrSges(6,3);
t432 = m(6) * t434 - t500 * mrSges(6,2) + t443 * mrSges(6,3) + t466 * t450 - t503 * t460;
t423 = t512 * t431 + t509 * t432;
t470 = -t480 * mrSges(5,1) + t481 * mrSges(5,2);
t476 = -qJD(4) * mrSges(5,2) + t480 * mrSges(5,3);
t421 = m(5) * t439 + qJDD(4) * mrSges(5,1) - t473 * mrSges(5,3) + qJD(4) * t476 - t481 * t470 + t423;
t477 = qJD(4) * mrSges(5,1) - t481 * mrSges(5,3);
t526 = -t509 * t431 + t512 * t432;
t422 = m(5) * t440 - qJDD(4) * mrSges(5,2) + t472 * mrSges(5,3) - qJD(4) * t477 + t480 * t470 + t526;
t417 = t513 * t421 + t510 * t422;
t457 = -t505 * t468 + t536;
t518 = mrSges(4,3) * qJDD(1) + t515 * (-mrSges(4,1) * t507 + t539);
t415 = m(4) * t457 - t518 * t505 + t417;
t527 = -t510 * t421 + t513 * t422;
t416 = m(4) * t458 + t518 * t507 + t527;
t528 = -t505 * t415 + t507 * t416;
t408 = m(3) * t475 - t515 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t528;
t474 = t508 * t487 - t506 * t488;
t522 = qJDD(3) - t474;
t462 = -qJDD(1) * pkin(2) - t515 * qJ(3) + t522;
t501 = t505 ^ 2;
t456 = (-pkin(2) - t540) * qJDD(1) + (-qJ(3) + (-t501 - t502) * pkin(6)) * t515 + t522;
t438 = -t472 * pkin(4) - t479 * pkin(7) + t481 * t478 + t456;
t521 = m(6) * t438 - t443 * mrSges(6,1) + t444 * mrSges(6,2) - t466 * t459 + t467 * t460;
t517 = m(5) * t456 - t472 * mrSges(5,1) + t473 * mrSges(5,2) - t480 * t476 + t481 * t477 + t521;
t516 = -m(4) * t462 + mrSges(4,1) * t532 - t517 + (t501 * t515 + t538) * mrSges(4,3);
t427 = t516 + (mrSges(3,1) - t539) * qJDD(1) - t515 * mrSges(3,2) + m(3) * t474;
t529 = t508 * t408 - t506 * t427;
t403 = m(2) * t489 - t515 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t529;
t405 = t506 * t408 + t508 * t427;
t404 = m(2) * t490 + qJDD(1) * mrSges(2,1) - t515 * mrSges(2,2) + t405;
t537 = t511 * t403 + t514 * t404;
t409 = t507 * t415 + t505 * t416;
t523 = Ifges(4,5) * t505 + Ifges(4,6) * t507;
t535 = t515 * t523;
t531 = m(3) * t504 + t409;
t530 = -t514 * t403 + t511 * t404;
t525 = Ifges(4,1) * t505 + Ifges(4,4) * t507;
t524 = Ifges(4,4) * t505 + Ifges(4,2) * t507;
t465 = Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * qJD(4);
t464 = Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * qJD(4);
t463 = Ifges(5,5) * t481 + Ifges(5,6) * t480 + Ifges(5,3) * qJD(4);
t447 = Ifges(6,1) * t467 + Ifges(6,4) * t466 + Ifges(6,5) * t503;
t446 = Ifges(6,4) * t467 + Ifges(6,2) * t466 + Ifges(6,6) * t503;
t445 = Ifges(6,5) * t467 + Ifges(6,6) * t466 + Ifges(6,3) * t503;
t425 = mrSges(6,2) * t438 - mrSges(6,3) * t433 + Ifges(6,1) * t444 + Ifges(6,4) * t443 + Ifges(6,5) * t500 + t466 * t445 - t503 * t446;
t424 = -mrSges(6,1) * t438 + mrSges(6,3) * t434 + Ifges(6,4) * t444 + Ifges(6,2) * t443 + Ifges(6,6) * t500 - t467 * t445 + t503 * t447;
t411 = mrSges(5,2) * t456 - mrSges(5,3) * t439 + Ifges(5,1) * t473 + Ifges(5,4) * t472 + Ifges(5,5) * qJDD(4) - pkin(7) * t423 - qJD(4) * t464 - t509 * t424 + t512 * t425 + t480 * t463;
t410 = -mrSges(5,1) * t456 + mrSges(5,3) * t440 + Ifges(5,4) * t473 + Ifges(5,2) * t472 + Ifges(5,6) * qJDD(4) - pkin(4) * t521 + pkin(7) * t526 + qJD(4) * t465 + t512 * t424 + t509 * t425 - t481 * t463;
t399 = mrSges(4,2) * t462 - mrSges(4,3) * t457 - pkin(6) * t417 + t525 * qJDD(1) - t510 * t410 + t513 * t411 + t507 * t535;
t398 = -mrSges(4,1) * t462 + mrSges(4,3) * t458 - pkin(3) * t517 + pkin(6) * t527 + t524 * qJDD(1) + t513 * t410 + t510 * t411 - t505 * t535;
t397 = -pkin(2) * t409 - Ifges(5,3) * qJDD(4) - mrSges(3,1) * t504 - Ifges(6,3) * t500 - t481 * t464 - pkin(3) * t417 - pkin(4) * t423 - mrSges(6,1) * t433 + mrSges(6,2) * t434 - mrSges(5,1) * t439 + mrSges(5,2) * t440 - Ifges(6,6) * t443 - Ifges(6,5) * t444 - mrSges(4,1) * t457 + mrSges(4,2) * t458 + t466 * t447 - t467 * t446 - Ifges(5,6) * t472 - Ifges(5,5) * t473 + mrSges(3,3) * t475 + t480 * t465 + (Ifges(3,6) - t523) * qJDD(1) + (-t505 * t524 + t507 * t525 + Ifges(3,5)) * t515;
t396 = mrSges(3,2) * t504 - mrSges(3,3) * t474 + Ifges(3,5) * qJDD(1) - t515 * Ifges(3,6) - qJ(3) * t409 - t505 * t398 + t507 * t399;
t395 = -mrSges(2,2) * g(1) - mrSges(2,3) * t490 + Ifges(2,5) * qJDD(1) - t515 * Ifges(2,6) - qJ(2) * t405 + t508 * t396 - t506 * t397;
t394 = mrSges(2,1) * g(1) + mrSges(2,3) * t489 + t515 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t531 + qJ(2) * t529 + t506 * t396 + t508 * t397;
t1 = [(-m(1) - m(2)) * g(1) + t531; -m(1) * g(2) + t537; -m(1) * g(3) + t530; pkin(1) * t405 + t505 * t399 + t507 * t398 + pkin(2) * t516 + qJ(3) * t528 + mrSges(3,1) * t474 - mrSges(3,2) * t475 + mrSges(2,1) * t490 - mrSges(2,2) * t489 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(2) * t539 + Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t530 + t514 * t394 + t511 * t395; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t537 + t511 * t394 - t514 * t395;];
tauB = t1;
