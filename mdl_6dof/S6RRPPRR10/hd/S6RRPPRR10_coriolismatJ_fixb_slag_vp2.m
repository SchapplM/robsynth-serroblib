% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:16
% EndTime: 2019-03-09 09:34:36
% DurationCPUTime: 11.04s
% Computational Cost: add. (26958->640), mult. (51181->854), div. (0->0), fcn. (54944->8), ass. (0->320)
t405 = sin(qJ(6));
t408 = cos(qJ(6));
t402 = sin(pkin(10));
t403 = cos(pkin(10));
t406 = sin(qJ(5));
t559 = cos(qJ(5));
t431 = t402 * t559 + t406 * t403;
t619 = -t406 * t402 + t559 * t403;
t296 = t405 * t431 - t408 * t619;
t450 = -t405 * t619 - t408 * t431;
t475 = Ifges(7,5) * t450 + Ifges(7,6) * t296;
t404 = -pkin(2) - qJ(4);
t550 = -pkin(8) + t404;
t375 = t550 * t402;
t376 = t550 * t403;
t307 = t375 * t559 + t406 * t376;
t237 = -pkin(9) * t431 + t307;
t306 = -t375 * t406 + t559 * t376;
t436 = -pkin(9) * t619 + t306;
t128 = t237 * t408 + t405 * t436;
t612 = -t237 * t405 + t408 * t436;
t669 = -t128 * mrSges(7,1) - t612 * mrSges(7,2);
t28 = t475 + t669;
t675 = t28 * qJD(6);
t407 = sin(qJ(2));
t409 = cos(qJ(2));
t339 = t619 * t409;
t341 = t431 * t409;
t452 = -t408 * t339 + t341 * t405;
t248 = -t339 * t405 - t341 * t408;
t544 = Ifges(7,4) * t248;
t131 = Ifges(7,2) * t452 + t407 * Ifges(7,6) + t544;
t232 = Ifges(7,4) * t452;
t133 = Ifges(7,1) * t248 + t407 * Ifges(7,5) + t232;
t288 = Ifges(7,4) * t450;
t179 = -Ifges(7,1) * t296 + t288;
t205 = mrSges(7,1) * t407 - mrSges(7,3) * t248;
t601 = pkin(3) + pkin(7);
t382 = t601 * t409;
t482 = t403 * t409;
t351 = pkin(4) * t482 + t382;
t281 = pkin(5) * t339 + t351;
t389 = t402 * pkin(4) + qJ(3);
t330 = pkin(5) * t431 + t389;
t203 = -mrSges(7,2) * t407 + mrSges(7,3) * t452;
t596 = t203 / 0.2e1;
t599 = -t128 / 0.2e1;
t600 = -t612 / 0.2e1;
t636 = Ifges(7,2) * t296 + t288;
t637 = -Ifges(7,2) * t248 + t232;
t645 = t248 * mrSges(7,1);
t653 = t452 * mrSges(7,2) + t645;
t625 = t450 * mrSges(7,2);
t654 = -t296 * mrSges(7,1) + t625;
t655 = Ifges(7,1) * t452 - t544;
t674 = (t131 / 0.4e1 - t655 / 0.4e1) * t296 + t205 * t599 + t653 * t330 / 0.2e1 + t654 * t281 / 0.2e1 + t612 * t596 + (t636 / 0.4e1 + mrSges(7,3) * t600 + t179 / 0.4e1) * t452 + (t637 / 0.4e1 + t133 / 0.4e1) * t450;
t658 = t450 * mrSges(7,1);
t464 = t658 / 0.2e1;
t502 = qJ(3) * t407;
t359 = t404 * t409 - pkin(1) - t502;
t381 = t601 * t407;
t369 = t403 * t381;
t250 = pkin(4) * t407 + t369 + (pkin(8) * t409 - t359) * t402;
t303 = t403 * t359 + t402 * t381;
t266 = -pkin(8) * t482 + t303;
t145 = t559 * t250 - t266 * t406;
t113 = pkin(9) * t341 + t145;
t109 = pkin(5) * t407 + t113;
t146 = t406 * t250 + t266 * t559;
t114 = -t339 * pkin(9) + t146;
t501 = t114 * t405;
t54 = t109 * t408 - t501;
t65 = t113 * t408 - t501;
t673 = t54 - t65;
t500 = t114 * t408;
t55 = t109 * t405 + t500;
t64 = -t113 * t405 - t500;
t660 = t55 + t64;
t662 = t450 / 0.2e1;
t672 = t205 * t662;
t671 = t281 * t653;
t670 = t330 * t654;
t643 = t296 * mrSges(7,2);
t668 = t658 + t643;
t543 = Ifges(7,4) * t296;
t176 = Ifges(7,2) * t450 - t543;
t656 = Ifges(7,1) * t450 + t543;
t667 = -t176 / 0.4e1 + t656 / 0.4e1;
t664 = -Ifges(7,3) / 0.2e1;
t663 = -t450 / 0.2e1;
t466 = t645 / 0.2e1;
t661 = m(7) * t330;
t338 = t619 * t407;
t340 = t431 * t407;
t243 = t338 * t408 - t340 * t405;
t451 = t338 * t405 + t408 * t340;
t602 = m(7) * pkin(5);
t469 = t602 / 0.2e1;
t590 = t451 / 0.2e1;
t593 = t243 / 0.2e1;
t478 = mrSges(7,1) * t590 + mrSges(7,2) * t593;
t571 = t340 / 0.2e1;
t573 = t338 / 0.2e1;
t657 = -mrSges(6,2) * t573 - mrSges(6,1) * t571 - (-t243 * t405 + t408 * t451) * t469 - t478;
t434 = Ifges(7,5) * t590 + Ifges(7,6) * t593;
t453 = t407 * pkin(2) - qJ(3) * t409;
t360 = qJ(4) * t407 + t453;
t370 = t403 * t382;
t251 = pkin(4) * t409 + t370 + (-pkin(8) * t407 - t360) * t402;
t305 = t403 * t360 + t402 * t382;
t483 = t403 * t407;
t267 = pkin(8) * t483 + t305;
t147 = t559 * t251 - t267 * t406;
t110 = pkin(5) * t409 - pkin(9) * t340 + t147;
t148 = t406 * t251 + t559 * t267;
t116 = pkin(9) * t338 + t148;
t60 = t110 * t408 - t116 * t405;
t61 = t110 * t405 + t116 * t408;
t618 = t61 * mrSges(7,2) / 0.2e1 - t60 * mrSges(7,1) / 0.2e1;
t448 = -t409 * t664 + t434 - t618;
t626 = Ifges(7,5) * t452;
t647 = Ifges(7,6) * t248;
t477 = t626 - t647;
t652 = t296 ^ 2 + t450 ^ 2;
t460 = t647 / 0.2e1 - t626 / 0.2e1;
t651 = t296 * t60 + t450 * t61;
t578 = -t296 / 0.2e1;
t648 = t296 / 0.2e1;
t465 = t625 / 0.2e1;
t635 = t296 * t405 - t408 * t450;
t631 = t452 / 0.2e1;
t630 = -Ifges(3,4) - Ifges(4,6);
t621 = t408 * t452;
t620 = t452 * t648;
t617 = -Ifges(6,5) * t431 - Ifges(6,6) * t619 + t475;
t616 = -Ifges(6,5) * t339 + Ifges(6,6) * t341 + t477;
t615 = -t147 * t619 - t148 * t431;
t545 = Ifges(6,4) * t619;
t310 = -Ifges(6,2) * t431 + t545;
t311 = -Ifges(6,1) * t431 - t545;
t614 = t310 / 0.4e1 - t311 / 0.4e1;
t613 = -t148 * mrSges(6,2) / 0.2e1 + t147 * mrSges(6,1) / 0.2e1;
t527 = t248 * mrSges(7,2);
t533 = t452 * mrSges(7,1);
t138 = t527 - t533;
t569 = t619 / 0.2e1;
t570 = -t341 / 0.2e1;
t604 = m(7) / 0.2e1;
t611 = (t281 * t619 - t330 * t341) * t604 + t138 * t569 - t668 * t570;
t546 = Ifges(6,4) * t341;
t239 = -Ifges(6,2) * t339 + t407 * Ifges(6,6) - t546;
t335 = Ifges(6,4) * t339;
t241 = -Ifges(6,1) * t341 + t407 * Ifges(6,5) - t335;
t253 = -t341 * mrSges(6,1) - t339 * mrSges(6,2);
t255 = Ifges(6,2) * t341 - t335;
t256 = -Ifges(6,1) * t339 + t546;
t308 = mrSges(6,1) * t619 - mrSges(6,2) * t431;
t358 = Ifges(6,4) * t431;
t309 = -Ifges(6,2) * t619 - t358;
t312 = Ifges(6,1) * t619 - t358;
t316 = -mrSges(6,2) * t407 - t339 * mrSges(6,3);
t560 = t407 / 0.4e1;
t594 = -t248 / 0.2e1;
t610 = -(t241 / 0.4e1 + t255 / 0.4e1) * t431 - (-t306 * mrSges(6,3) / 0.2e1 + t312 / 0.4e1 + t309 / 0.4e1) * t339 + (t128 * t594 + t54 * t663 + t648 * t660 + t662 * t65) * mrSges(7,3) + t306 * t316 / 0.2e1 + t351 * t308 / 0.2e1 + t389 * t253 / 0.2e1 + t617 * t560 + (-t128 * t673 + t612 * t660) * t604 - (t239 / 0.4e1 - t256 / 0.4e1) * t619 + t667 * t248 + t674;
t401 = t403 ^ 2;
t609 = -m(5) / 0.2e1;
t608 = m(5) / 0.2e1;
t607 = -m(6) / 0.2e1;
t606 = m(6) / 0.2e1;
t605 = -m(7) / 0.2e1;
t589 = t248 / 0.2e1;
t514 = t341 * mrSges(6,3);
t318 = mrSges(6,1) * t407 + t514;
t575 = -t318 / 0.2e1;
t572 = -t339 / 0.2e1;
t568 = -t619 / 0.2e1;
t567 = -t431 / 0.2e1;
t566 = -t402 / 0.2e1;
t565 = t402 / 0.2e1;
t564 = -t403 / 0.2e1;
t563 = t403 / 0.2e1;
t562 = -t405 / 0.2e1;
t561 = t407 / 0.2e1;
t558 = t341 * pkin(5);
t556 = t54 * mrSges(7,2);
t555 = t55 * mrSges(7,1);
t552 = t64 * mrSges(7,1);
t551 = t65 * mrSges(7,2);
t548 = Ifges(5,4) * t402;
t547 = Ifges(5,4) * t403;
t541 = pkin(5) * qJD(5);
t534 = t243 * mrSges(7,1);
t529 = t451 * mrSges(7,2);
t130 = Ifges(7,4) * t451 + Ifges(7,2) * t243 + Ifges(7,6) * t409;
t132 = Ifges(7,1) * t451 + Ifges(7,4) * t243 + Ifges(7,5) * t409;
t137 = t529 - t534;
t202 = -mrSges(7,2) * t409 + mrSges(7,3) * t243;
t204 = mrSges(7,1) * t409 - mrSges(7,3) * t451;
t238 = Ifges(6,4) * t340 + Ifges(6,2) * t338 + Ifges(6,6) * t409;
t240 = Ifges(6,1) * t340 + Ifges(6,4) * t338 + Ifges(6,5) * t409;
t516 = t340 * mrSges(6,2);
t518 = t338 * mrSges(6,1);
t254 = t516 - t518;
t350 = (-pkin(4) * t403 - t601) * t407;
t280 = -pkin(5) * t338 + t350;
t302 = -t359 * t402 + t369;
t304 = -t360 * t402 + t370;
t315 = -mrSges(6,2) * t409 + mrSges(6,3) * t338;
t317 = mrSges(6,1) * t409 - t340 * mrSges(6,3);
t336 = Ifges(5,6) * t409 + (Ifges(5,2) * t403 + t548) * t407;
t512 = t402 * Ifges(5,1);
t337 = Ifges(5,5) * t409 + (t512 + t547) * t407;
t510 = t403 * mrSges(5,1);
t513 = t402 * mrSges(5,2);
t447 = t510 - t513;
t352 = t447 * t407;
t484 = t402 * t407;
t506 = t409 * mrSges(5,1);
t371 = -mrSges(5,3) * t484 + t506;
t508 = t407 * mrSges(5,1);
t372 = mrSges(5,3) * t402 * t409 + t508;
t505 = t409 * mrSges(5,2);
t373 = mrSges(5,3) * t483 - t505;
t507 = t407 * mrSges(5,2);
t374 = -mrSges(5,3) * t482 - t507;
t446 = -pkin(2) * t409 - t502;
t378 = -pkin(1) + t446;
t435 = Ifges(6,5) * t571 + Ifges(6,6) * t573;
t425 = t434 + t435;
t379 = t409 * mrSges(4,2) - t407 * mrSges(4,3);
t455 = m(4) * t378 + t379;
t509 = t403 * Ifges(5,6);
t511 = t402 * Ifges(5,5);
t515 = t341 * mrSges(6,2);
t517 = t339 * mrSges(6,1);
t3 = (-pkin(1) * mrSges(3,1) - t378 * mrSges(4,2) + (t509 + t511 + t630) * t407 + (-t401 * Ifges(5,2) / 0.2e1 + Ifges(6,3) + Ifges(7,3) - Ifges(3,2) + Ifges(3,1) + Ifges(4,2) - Ifges(4,3) + Ifges(5,3) + (-t547 - t512 / 0.2e1) * t402) * t409 + t425) * t407 + t350 * (-t515 + t517) + t455 * t453 + (-pkin(1) * mrSges(3,2) - t378 * mrSges(4,3) + Ifges(6,5) * t570 + Ifges(7,5) * t589 + Ifges(6,6) * t572 + Ifges(7,6) * t631 + t336 * t564 + t337 * t566 - t381 * t447 + (-t511 / 0.2e1 - t509 / 0.2e1 - t630) * t409) * t409 + t130 * t631 + t305 * t374 - t382 * t352 + t302 * t371 + t304 * t372 + t303 * t373 + m(5) * (t302 * t304 + t303 * t305 - t381 * t382) + m(6) * (t145 * t147 + t146 * t148 + t350 * t351) + m(7) * (t280 * t281 + t54 * t60 + t55 * t61) + t351 * t254 + t145 * t317 + t147 * t318 + t146 * t315 + t148 * t316 + t280 * t138 + t281 * t137 + t60 * t205 + t55 * t202 + t61 * t203 + t54 * t204 + t240 * t570 + t241 * t571 + t238 * t572 + t239 * t573 + t132 * t589 + t133 * t590 + t131 * t593;
t522 = t3 * qJD(1);
t6 = m(7) * (-t281 * t558 + t54 * t64 + t55 * t65) + t65 * t203 + t64 * t205 + t131 * t594 - t138 * t558 + t671 + t655 * t589 + t145 * t316 + t256 * t570 + t341 * t239 / 0.2e1 + t351 * t253 + (-t318 + t514) * t146 + (-t248 * t55 - t452 * t54) * mrSges(7,3) - (-t145 * mrSges(6,3) + t241 / 0.2e1 + t255 / 0.2e1) * t339 + t616 * t561 + (t637 + t133) * t631;
t504 = t6 * qJD(1);
t7 = t54 * t203 - t55 * t205 + t671 + t477 * t561 + (-t55 * mrSges(7,3) - t131 / 0.2e1 + t655 / 0.2e1) * t248 + (-t54 * mrSges(7,3) + t637 / 0.2e1 + t133 / 0.2e1) * t452;
t503 = t7 * qJD(1);
t422 = (t248 * t662 + t620) * mrSges(7,3) - t296 * t596 + t672;
t17 = t422 - t478;
t497 = t17 * qJD(1);
t424 = m(5) * (t302 * t402 - t303 * t403) + t402 * t372 - t403 * t374;
t21 = -t243 * t203 + t451 * t205 - t338 * t316 + t340 * t318 + m(7) * (-t243 * t55 + t451 * t54) + m(6) * (t145 * t340 - t146 * t338) + (t424 - t455) * t407;
t496 = t21 * qJD(1);
t493 = t248 * t408;
t492 = t452 * t405;
t491 = t296 * t408;
t490 = t450 * t405;
t489 = t338 * t431;
t488 = t339 * t431;
t487 = t340 * t619;
t486 = t341 * t619;
t481 = t405 * t202;
t479 = t408 * t204;
t471 = t402 ^ 2 + t401;
t470 = -t602 / 0.2e1;
t457 = t431 ^ 2 + t619 ^ 2;
t456 = m(5) * t471;
t454 = mrSges(6,1) * t431 + mrSges(6,2) * t619;
t449 = t609 + t607 + t605;
t22 = t670 - (-t176 / 0.2e1 + t656 / 0.2e1) * t296 + (t636 / 0.2e1 + t179 / 0.2e1) * t450;
t416 = (mrSges(7,3) * t599 + t667) * t248 + t475 * t560 + t674;
t4 = t416 - t448;
t445 = t4 * qJD(1) + t22 * qJD(2);
t414 = (-t248 * t648 + t452 * t662) * mrSges(7,3) + (t488 / 0.2e1 - t486 / 0.2e1) * mrSges(6,3) + (-t302 * t403 - t303 * t402) * t608 + (-t145 * t619 - t146 * t431 + t306 * t341 - t307 * t339) * t606 + (t128 * t452 - t248 * t612 + t296 * t54 + t450 * t55) * t604 + t205 * t648 + t203 * t662 + t318 * t568 + t316 * t567;
t419 = t381 * t608 + t350 * t607 + t280 * t605 + t534 / 0.2e1 - t529 / 0.2e1 + t518 / 0.2e1 - t516 / 0.2e1;
t12 = t414 + (-t372 / 0.2e1 + t508 / 0.2e1) * t403 + (-t374 / 0.2e1 - t507 / 0.2e1) * t402 + t419;
t33 = t652 * mrSges(7,3) + t457 * mrSges(6,3) + t471 * mrSges(5,3) + m(7) * (t128 * t450 + t296 * t612) + m(6) * (-t306 * t619 - t307 * t431) - t404 * t456;
t444 = qJD(1) * t12 + qJD(2) * t33;
t415 = (-t243 * t662 + t451 * t648) * mrSges(7,3) + (t489 / 0.2e1 - t487 / 0.2e1) * mrSges(6,3) + t382 * t608 + (t306 * t340 - t307 * t338 + t351) * t606 + (-t128 * t243 + t451 * t612 + t281) * t604 - t533 / 0.2e1 + t527 / 0.2e1 + t517 / 0.2e1 - t515 / 0.2e1;
t437 = t304 * t403 + t305 * t402;
t418 = t202 * t662 + t204 * t648 + t315 * t567 + t317 * t568 + t437 * t609 - t605 * t651 - t615 * t607;
t14 = (-t371 / 0.2e1 + t506 / 0.2e1) * t403 + (-t373 / 0.2e1 - t505 / 0.2e1) * t402 + t415 + t418;
t430 = t668 - t454;
t76 = t402 * mrSges(5,1) + t403 * mrSges(5,2) + mrSges(4,3) + (m(5) + m(4)) * qJ(3) + m(6) * t389 + t661 - t430;
t443 = qJD(1) * t14 + qJD(2) * t76;
t49 = (t570 + t493 / 0.2e1 - t492 / 0.2e1) * t602 + t653 + t253;
t71 = (t491 / 0.2e1 + t490 / 0.2e1 + t568) * t602 - t654 - t308;
t442 = qJD(1) * t49 - qJD(2) * t71;
t53 = (t248 * t296 - t450 * t452) * t604 + (t486 - t488) * t606;
t421 = -t652 * t604 + t457 * t607 - t456 / 0.2e1;
t67 = t421 + t449;
t441 = qJD(1) * t53 + qJD(2) * t67;
t72 = 0.2e1 * t631 * mrSges(7,2) + 0.2e1 * t466;
t81 = 0.2e1 * mrSges(7,1) * t578 + 0.2e1 * t465;
t440 = qJD(1) * t72 + qJD(2) * t81;
t23 = t452 * t203 - t248 * t205 - t339 * t316 + t341 * t318 + m(7) * (-t248 * t54 + t452 * t55) + m(6) * (t145 * t341 - t146 * t339) + t424 * t409;
t439 = qJD(1) * t23 + qJD(3) * t53;
t432 = m(7) * (t405 * t61 + t408 * t60);
t1 = (t307 * mrSges(6,3) / 0.2e1 + t614) * t341 + (-t479 / 0.2e1 - t481 / 0.2e1 - t432 / 0.2e1 + t611) * pkin(5) + (-Ifges(6,3) / 0.2e1 + t664) * t409 + t307 * t575 - t425 + t610 - t613 + t618;
t10 = t670 + t176 * t648 + t656 * t578 + t389 * t308 + t311 * t569 + t310 * t568 - (t312 / 0.2e1 + t309 / 0.2e1) * t431 + (t636 + t179) * t662 + (-t668 + t661) * t619 * pkin(5);
t429 = t1 * qJD(1) + t10 * qJD(2);
t417 = (-t248 * t663 + t620) * mrSges(7,3) + (-t339 * t568 - t431 * t570) * mrSges(6,3) + (-t296 * t660 + t450 * t673) * t604 + t672 + t203 * t578 + t316 * t569 + t318 * t567;
t9 = t417 + t657;
t428 = t9 * qJD(1);
t420 = (t408 * t596 + t205 * t562 + (t248 * t562 - t621 / 0.2e1) * mrSges(7,3)) * pkin(5) - t460;
t15 = (-t54 / 0.2e1 + t65 / 0.2e1) * mrSges(7,2) + (-t55 / 0.2e1 - t64 / 0.2e1) * mrSges(7,1) + t420 + t460;
t29 = (t600 + t612 / 0.2e1) * mrSges(7,2) + (t599 + t128 / 0.2e1) * mrSges(7,1);
t377 = (mrSges(7,1) * t405 + mrSges(7,2) * t408) * pkin(5);
t62 = t464 - t658 / 0.2e1 + (t648 + t578) * mrSges(7,2);
t426 = -qJD(1) * t15 - qJD(2) * t29 - qJD(3) * t62 + qJD(5) * t377;
t365 = t377 * qJD(6);
t129 = -t619 * t470 + (t490 + t491) * t469;
t97 = t341 * t470 + (t492 - t493) * t469;
t82 = t465 - t625 / 0.2e1;
t73 = t466 - t645 / 0.2e1;
t66 = t421 - t449;
t63 = t643 + 0.2e1 * t464;
t50 = t53 * qJD(4);
t18 = t422 + t478;
t16 = -t556 / 0.2e1 - t555 / 0.2e1 - t551 / 0.2e1 + t552 / 0.2e1 + t420 - t460;
t13 = (m(4) * pkin(7) + mrSges(4,1) - t513 / 0.2e1 + t510 / 0.2e1) * t409 + t415 + t373 * t565 + t371 * t563 - t418;
t11 = t414 - mrSges(5,1) * t483 / 0.2e1 + mrSges(5,2) * t484 / 0.2e1 + t374 * t566 + t372 * t564 - t419;
t8 = t417 - t657;
t5 = t416 + t448;
t2 = t448 + t435 + (t514 / 0.2e1 + t575) * t307 + Ifges(6,3) * t409 / 0.2e1 + t614 * t341 + t611 * pkin(5) + (t479 + t481 + t432) * pkin(5) / 0.2e1 + t610 + t613;
t19 = [qJD(2) * t3 + qJD(3) * t21 + qJD(4) * t23 + qJD(5) * t6 + qJD(6) * t7, t13 * qJD(3) + t11 * qJD(4) + t2 * qJD(5) + t5 * qJD(6) + t522 + (t615 * mrSges(6,3) + t350 * t454 + t651 * mrSges(7,3) + t612 * t204 + 0.2e1 * (t128 * t61 + t280 * t330 + t60 * t612) * t604 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t563 + Ifges(6,5) * t569 + Ifges(7,5) * t578 + Ifges(5,6) * t566 + Ifges(6,6) * t567 + Ifges(7,6) * t662 - Ifges(4,4) + Ifges(3,5)) * t409 + t130 * t662 + t389 * t254 - qJ(3) * t352 + t330 * t137 + t306 * t317 + t307 * t315 + t128 * t202 + (-t304 * mrSges(5,3) + t404 * t371 - t381 * mrSges(5,2) + t337 / 0.2e1) * t403 - t280 * t668 + (t404 * t373 - t305 * mrSges(5,3) - t381 * mrSges(5,1) - t336 / 0.2e1) * t402 + (-qJ(3) * mrSges(4,1) + (Ifges(5,1) * t403 - t548) * t565 + (-Ifges(5,2) * t402 + t547) * t563 + Ifges(4,5) - Ifges(3,6)) * t407 + (m(4) * t446 - t409 * mrSges(3,1) + t407 * mrSges(3,2) + t379) * pkin(7) + t238 * t567 + t240 * t569 + t312 * t571 + t310 * t573 + t132 * t578 + t179 * t590 + t176 * t593 + 0.2e1 * (t147 * t306 + t148 * t307 + t350 * t389) * t606 + 0.2e1 * (-qJ(3) * t381 + t404 * t437) * t608) * qJD(2), t496 + t13 * qJD(2) + 0.2e1 * ((t243 * t450 - t296 * t451) * t604 + (t487 - t489) * t606) * qJD(3) + t50 + t8 * qJD(5) + t18 * qJD(6), qJD(2) * t11 + qJD(5) * t97 + qJD(6) * t73 + t439, t504 + t2 * qJD(2) + t8 * qJD(3) + t97 * qJD(4) + (-t146 * mrSges(6,1) - t145 * mrSges(6,2) - t551 + t552 + t616) * qJD(5) + t16 * qJD(6) + (m(7) * (t405 * t65 + t408 * t64) + (-t248 * t405 - t621) * mrSges(7,3)) * t541, t503 + t5 * qJD(2) + t18 * qJD(3) + t73 * qJD(4) + t16 * qJD(5) + (t477 - t555 - t556) * qJD(6); qJD(3) * t14 + qJD(4) * t12 + qJD(5) * t1 + qJD(6) * t4 - t522, qJD(3) * t76 + qJD(4) * t33 + qJD(5) * t10 + qJD(6) * t22, qJD(4) * t66 + t443, qJD(3) * t66 + qJD(5) * t129 + qJD(6) * t82 + t444, t129 * qJD(4) + (-t307 * mrSges(6,1) - t306 * mrSges(6,2) + t617 + t669) * qJD(5) + t675 + (m(7) * (-t128 * t408 + t405 * t612) + t635 * mrSges(7,3)) * t541 + t429, t82 * qJD(4) + t28 * qJD(5) + t445 + t675; -qJD(2) * t14 + qJD(5) * t9 + qJD(6) * t17 - t496 + t50, qJD(4) * t67 - t443, 0, t441 (-t602 * t635 + t430) * qJD(5) + t63 * qJD(6) + t428, t63 * qJD(5) + qJD(6) * t668 + t497; -qJD(2) * t12 + qJD(5) * t49 + qJD(6) * t72 - t439, -qJD(3) * t67 - qJD(5) * t71 + qJD(6) * t81 - t444, -t441, 0, t442, t440; -qJD(2) * t1 - qJD(3) * t9 - qJD(4) * t49 + qJD(6) * t15 - t504, qJD(4) * t71 + qJD(6) * t29 - t429, qJD(6) * t62 - t428, -t442, -t365, -t365 - t426; -qJD(2) * t4 - qJD(3) * t17 - qJD(4) * t72 - qJD(5) * t15 - t503, -qJD(4) * t81 - qJD(5) * t29 - t445, -t62 * qJD(5) - t497, -t440, t426, 0;];
Cq  = t19;
