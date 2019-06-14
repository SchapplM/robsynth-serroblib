% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:15:12
% EndTime: 2019-05-05 14:15:16
% DurationCPUTime: 4.36s
% Computational Cost: add. (55452->298), mult. (109701->365), div. (0->0), fcn. (60013->10), ass. (0->119)
t532 = sin(pkin(10));
t534 = cos(pkin(10));
t537 = sin(qJ(4));
t540 = cos(qJ(4));
t506 = (t532 * t537 - t534 * t540) * qJD(1);
t569 = 2 * qJD(5);
t568 = -pkin(1) - pkin(2);
t567 = mrSges(2,1) + mrSges(3,1);
t566 = Ifges(3,4) + Ifges(2,5);
t565 = Ifges(2,6) - Ifges(3,6);
t538 = sin(qJ(1));
t541 = cos(qJ(1));
t520 = -t541 * g(1) - t538 * g(2);
t543 = qJD(1) ^ 2;
t550 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t520;
t502 = -pkin(1) * t543 + t550;
t497 = t568 * t543 + t550;
t519 = t538 * g(1) - t541 * g(2);
t549 = -t543 * qJ(2) + qJDD(2) - t519;
t501 = t568 * qJDD(1) + t549;
t533 = sin(pkin(9));
t535 = cos(pkin(9));
t479 = t535 * t497 + t533 * t501;
t476 = -pkin(3) * t543 - qJDD(1) * pkin(7) + t479;
t529 = g(3) + qJDD(3);
t472 = -t537 * t476 + t540 * t529;
t561 = qJD(1) * qJD(4);
t560 = t540 * t561;
t514 = -qJDD(1) * t537 - t560;
t462 = (-t514 - t560) * qJ(5) + (t537 * t540 * t543 + qJDD(4)) * pkin(4) + t472;
t473 = t540 * t476 + t537 * t529;
t515 = -qJDD(1) * t540 + t537 * t561;
t563 = qJD(1) * t537;
t516 = qJD(4) * pkin(4) + qJ(5) * t563;
t528 = t540 ^ 2;
t463 = -pkin(4) * t528 * t543 + qJ(5) * t515 - qJD(4) * t516 + t473;
t458 = t532 * t462 + t534 * t463 + t506 * t569;
t507 = (t532 * t540 + t534 * t537) * qJD(1);
t486 = -mrSges(6,1) * t506 - mrSges(6,2) * t507;
t490 = -t514 * t532 + t515 * t534;
t499 = qJD(4) * mrSges(6,1) + mrSges(6,3) * t507;
t487 = -pkin(5) * t506 + pkin(8) * t507;
t542 = qJD(4) ^ 2;
t456 = -pkin(5) * t542 + qJDD(4) * pkin(8) + t487 * t506 + t458;
t478 = -t533 * t497 + t535 * t501;
t554 = qJDD(1) * pkin(3) - t478;
t464 = -t516 * t563 - t515 * pkin(4) + qJDD(5) + (-qJ(5) * t528 - pkin(7)) * t543 + t554;
t491 = t514 * t534 + t515 * t532;
t459 = (-qJD(4) * t506 - t491) * pkin(8) + (-qJD(4) * t507 - t490) * pkin(5) + t464;
t536 = sin(qJ(6));
t539 = cos(qJ(6));
t453 = -t456 * t536 + t459 * t539;
t494 = qJD(4) * t539 + t507 * t536;
t471 = qJD(6) * t494 + qJDD(4) * t536 + t491 * t539;
t495 = qJD(4) * t536 - t507 * t539;
t477 = -mrSges(7,1) * t494 + mrSges(7,2) * t495;
t504 = qJD(6) - t506;
t480 = -mrSges(7,2) * t504 + mrSges(7,3) * t494;
t489 = qJDD(6) - t490;
t451 = m(7) * t453 + mrSges(7,1) * t489 - mrSges(7,3) * t471 - t477 * t495 + t480 * t504;
t454 = t456 * t539 + t459 * t536;
t470 = -qJD(6) * t495 + qJDD(4) * t539 - t491 * t536;
t481 = mrSges(7,1) * t504 - mrSges(7,3) * t495;
t452 = m(7) * t454 - mrSges(7,2) * t489 + mrSges(7,3) * t470 + t477 * t494 - t481 * t504;
t555 = -t451 * t536 + t539 * t452;
t442 = m(6) * t458 - qJDD(4) * mrSges(6,2) + mrSges(6,3) * t490 - qJD(4) * t499 + t486 * t506 + t555;
t553 = -t534 * t462 + t532 * t463;
t457 = t507 * t569 - t553;
t498 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t506;
t455 = -qJDD(4) * pkin(5) - t542 * pkin(8) + (-(2 * qJD(5)) - t487) * t507 + t553;
t547 = -m(7) * t455 + t470 * mrSges(7,1) - mrSges(7,2) * t471 + t494 * t480 - t481 * t495;
t447 = m(6) * t457 + qJDD(4) * mrSges(6,1) - mrSges(6,3) * t491 + qJD(4) * t498 + t486 * t507 + t547;
t438 = t532 * t442 + t534 * t447;
t513 = (mrSges(5,1) * t540 - mrSges(5,2) * t537) * qJD(1);
t562 = qJD(1) * t540;
t518 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t562;
t436 = m(5) * t472 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t514 + qJD(4) * t518 + t513 * t563 + t438;
t517 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t563;
t556 = t534 * t442 - t447 * t532;
t437 = m(5) * t473 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t515 - qJD(4) * t517 - t513 * t562 + t556;
t557 = -t436 * t537 + t540 * t437;
t430 = m(4) * t479 - mrSges(4,1) * t543 + qJDD(1) * mrSges(4,2) + t557;
t475 = -t543 * pkin(7) + t554;
t443 = t539 * t451 + t536 * t452;
t545 = m(6) * t464 - t490 * mrSges(6,1) + mrSges(6,2) * t491 - t506 * t498 - t499 * t507 + t443;
t544 = -m(5) * t475 + t515 * mrSges(5,1) - mrSges(5,2) * t514 + t517 * t563 - t518 * t562 - t545;
t439 = m(4) * t478 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t543 + t544;
t558 = t535 * t430 - t533 * t439;
t551 = m(3) * t502 + qJDD(1) * mrSges(3,3) + t558;
t426 = m(2) * t520 - qJDD(1) * mrSges(2,2) - t567 * t543 + t551;
t428 = t430 * t533 + t439 * t535;
t505 = -qJDD(1) * pkin(1) + t549;
t546 = -m(3) * t505 + qJDD(1) * mrSges(3,1) + t543 * mrSges(3,3) - t428;
t427 = m(2) * t519 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t543 + t546;
t564 = t538 * t426 + t541 * t427;
t559 = t541 * t426 - t427 * t538;
t432 = t540 * t436 + t537 * t437;
t548 = -m(4) * t529 - t432;
t510 = Ifges(5,5) * qJD(4) + (-Ifges(5,1) * t537 - Ifges(5,4) * t540) * qJD(1);
t509 = Ifges(5,6) * qJD(4) + (-Ifges(5,4) * t537 - Ifges(5,2) * t540) * qJD(1);
t508 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t537 - Ifges(5,6) * t540) * qJD(1);
t484 = -Ifges(6,1) * t507 + Ifges(6,4) * t506 + Ifges(6,5) * qJD(4);
t483 = -Ifges(6,4) * t507 + Ifges(6,2) * t506 + Ifges(6,6) * qJD(4);
t482 = -Ifges(6,5) * t507 + Ifges(6,6) * t506 + Ifges(6,3) * qJD(4);
t467 = Ifges(7,1) * t495 + Ifges(7,4) * t494 + Ifges(7,5) * t504;
t466 = Ifges(7,4) * t495 + Ifges(7,2) * t494 + Ifges(7,6) * t504;
t465 = Ifges(7,5) * t495 + Ifges(7,6) * t494 + Ifges(7,3) * t504;
t445 = mrSges(7,2) * t455 - mrSges(7,3) * t453 + Ifges(7,1) * t471 + Ifges(7,4) * t470 + Ifges(7,5) * t489 + t465 * t494 - t466 * t504;
t444 = -mrSges(7,1) * t455 + mrSges(7,3) * t454 + Ifges(7,4) * t471 + Ifges(7,2) * t470 + Ifges(7,6) * t489 - t465 * t495 + t467 * t504;
t434 = -mrSges(6,1) * t464 - mrSges(7,1) * t453 + mrSges(7,2) * t454 + mrSges(6,3) * t458 + Ifges(6,4) * t491 - Ifges(7,5) * t471 + Ifges(6,2) * t490 + Ifges(6,6) * qJDD(4) - Ifges(7,6) * t470 - Ifges(7,3) * t489 - pkin(5) * t443 + qJD(4) * t484 - t466 * t495 + t467 * t494 + t482 * t507;
t433 = mrSges(6,2) * t464 - mrSges(6,3) * t457 + Ifges(6,1) * t491 + Ifges(6,4) * t490 + Ifges(6,5) * qJDD(4) - pkin(8) * t443 - qJD(4) * t483 - t444 * t536 + t445 * t539 + t482 * t506;
t431 = -m(3) * g(3) + t548;
t422 = mrSges(5,2) * t475 - mrSges(5,3) * t472 + Ifges(5,1) * t514 + Ifges(5,4) * t515 + Ifges(5,5) * qJDD(4) - qJ(5) * t438 - qJD(4) * t509 + t433 * t534 - t434 * t532 - t508 * t562;
t421 = -mrSges(5,1) * t475 + mrSges(5,3) * t473 + Ifges(5,4) * t514 + Ifges(5,2) * t515 + Ifges(5,6) * qJDD(4) - pkin(4) * t545 + qJ(5) * t556 + qJD(4) * t510 + t532 * t433 + t534 * t434 + t508 * t563;
t420 = -pkin(3) * t432 - mrSges(4,1) * t529 + mrSges(4,3) * t479 - pkin(4) * t438 - Ifges(5,5) * t514 - Ifges(5,6) * t515 - mrSges(5,1) * t472 + mrSges(5,2) * t473 - t536 * t445 - t539 * t444 - pkin(5) * t547 - pkin(8) * t555 - Ifges(6,5) * t491 - Ifges(6,6) * t490 - mrSges(6,1) * t457 + mrSges(6,2) * t458 + t543 * Ifges(4,5) + t507 * t483 + t506 * t484 - Ifges(4,6) * qJDD(1) + (-Ifges(5,3) - Ifges(6,3)) * qJDD(4) + (t509 * t537 - t510 * t540) * qJD(1);
t419 = mrSges(4,2) * t529 - mrSges(4,3) * t478 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t543 - pkin(7) * t432 - t421 * t537 + t422 * t540;
t418 = mrSges(3,2) * t505 - mrSges(2,3) * t519 - qJ(2) * t431 - qJ(3) * t428 + t535 * t419 - t533 * t420 - t565 * t543 + t566 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t417 = mrSges(3,2) * t502 + mrSges(2,3) * t520 - pkin(1) * t431 - pkin(2) * t548 + t567 * g(3) - qJ(3) * t558 + t565 * qJDD(1) - t533 * t419 - t535 * t420 + t566 * t543;
t1 = [-m(1) * g(1) + t559; -m(1) * g(2) + t564; (-m(1) - m(2) - m(3)) * g(3) + t548; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t564 - t538 * t417 + t541 * t418; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t559 + t541 * t417 + t538 * t418; pkin(1) * t546 + qJ(2) * (-mrSges(3,1) * t543 + t551) - mrSges(2,2) * t520 + mrSges(2,1) * t519 - pkin(2) * t428 - mrSges(3,1) * t505 + mrSges(3,3) * t502 - t537 * t422 - t540 * t421 - pkin(3) * t544 - pkin(7) * t557 - mrSges(4,1) * t478 + mrSges(4,2) * t479 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB  = t1;
