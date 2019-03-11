% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:08:01
% EndTime: 2019-03-09 10:08:29
% DurationCPUTime: 15.79s
% Computational Cost: add. (43839->631), mult. (85132->846), div. (0->0), fcn. (103381->10), ass. (0->342)
t393 = sin(qJ(4));
t390 = sin(pkin(10));
t392 = cos(pkin(10));
t394 = sin(qJ(2));
t555 = -qJ(3) - pkin(7);
t462 = t555 * t394;
t395 = cos(qJ(2));
t616 = t555 * t395;
t338 = -t390 * t462 + t392 * t616;
t362 = -t390 * t394 + t392 * t395;
t422 = -t362 * pkin(8) + t338;
t564 = cos(qJ(4));
t363 = -t390 * t395 - t392 * t394;
t634 = t390 * t616 + t392 * t462;
t658 = t363 * pkin(8) + t634;
t206 = t393 * t422 + t564 * t658;
t389 = sin(pkin(11));
t391 = cos(pkin(11));
t562 = sin(qJ(6));
t563 = cos(qJ(6));
t364 = -t562 * t389 + t563 * t391;
t426 = t389 * t563 + t391 * t562;
t328 = -mrSges(7,1) * t364 + mrSges(7,2) * t426;
t549 = Ifges(7,4) * t426;
t331 = Ifges(7,2) * t364 + t549;
t359 = Ifges(7,4) * t364;
t333 = Ifges(7,1) * t426 + t359;
t570 = t426 / 0.2e1;
t323 = t564 * t362 + t363 * t393;
t645 = t364 * t323;
t591 = t645 / 0.2e1;
t642 = t426 * t323;
t593 = -t642 / 0.2e1;
t647 = t364 / 0.2e1;
t433 = t393 * t362 - t363 * t564;
t651 = Ifges(7,4) * t645 - Ifges(7,2) * t642 + Ifges(7,6) * t433;
t652 = Ifges(7,1) * t645 - Ifges(7,4) * t642 + Ifges(7,5) * t433;
t644 = t389 * t323;
t688 = t393 * t658 - t564 * t422;
t690 = pkin(5) * t644 + t688;
t702 = -t206 * mrSges(5,2) + t690 * t328 + t331 * t593 + t333 * t591 + t652 * t570 + t651 * t647;
t501 = t433 * t389;
t143 = pkin(5) * t501 - t206;
t700 = t143 * t690;
t560 = pkin(2) * t392;
t379 = pkin(3) + t560;
t561 = pkin(2) * t390;
t351 = t393 * t379 + t561 * t564;
t699 = t206 * t351;
t698 = t206 * t389;
t509 = t206 * t688;
t350 = t379 * t564 - t393 * t561;
t349 = -pkin(4) - t350;
t556 = t391 * pkin(5);
t343 = t349 - t556;
t697 = t343 * t690;
t380 = -pkin(4) - t556;
t696 = t380 * t690;
t695 = t391 * t206;
t420 = t426 * t433;
t618 = t364 * t433;
t120 = mrSges(7,1) * t420 + mrSges(7,2) * t618;
t525 = t391 * mrSges(6,2);
t526 = t389 * mrSges(6,1);
t450 = t525 + t526;
t247 = t450 * t433;
t500 = t433 * t391;
t381 = -pkin(2) * t395 - pkin(1);
t346 = -t362 * pkin(3) + t381;
t202 = -pkin(4) * t323 - qJ(5) * t433 + t346;
t95 = t391 * t202 - t389 * t688;
t69 = -pkin(5) * t323 - pkin(9) * t500 + t95;
t96 = t389 * t202 + t391 * t688;
t74 = -pkin(9) * t501 + t96;
t42 = -t562 * t74 + t563 * t69;
t43 = t562 * t69 + t563 * t74;
t590 = t618 / 0.2e1;
t592 = -t420 / 0.2e1;
t640 = t450 * t323;
t653 = mrSges(7,1) * t433 - mrSges(7,3) * t645;
t654 = -mrSges(7,2) * t433 - mrSges(7,3) * t642;
t643 = t391 * t323;
t655 = mrSges(6,1) * t433 - mrSges(6,3) * t643;
t656 = -mrSges(6,2) * t433 - mrSges(6,3) * t644;
t664 = t645 * mrSges(7,2);
t665 = t642 * mrSges(7,1);
t675 = t664 + t665;
t693 = t690 * t120 + t143 * t675 - t206 * t640 + t688 * t247 + t42 * t653 + t43 * t654 + t652 * t590 + t651 * t592 + t95 * t655 + t96 * t656;
t692 = pkin(4) * t688;
t691 = t349 * t688;
t687 = t380 * t675;
t486 = -t665 / 0.2e1 - t664 / 0.2e1;
t686 = (t525 / 0.2e1 + t526 / 0.2e1) * t323 - t486;
t685 = -t363 * mrSges(4,1) + t362 * mrSges(4,2);
t518 = qJ(5) * t655;
t682 = qJ(5) * t656;
t348 = qJ(5) + t351;
t681 = t348 * t656;
t546 = Ifges(6,6) * t389;
t548 = Ifges(6,5) * t391;
t447 = -t546 + t548;
t565 = t391 / 0.2e1;
t566 = -t389 / 0.2e1;
t583 = t433 / 0.2e1;
t630 = -t433 / 0.2e1;
t552 = Ifges(6,4) * t389;
t449 = Ifges(6,1) * t391 - t552;
t635 = Ifges(6,5) * t433 + t323 * t449;
t551 = Ifges(6,4) * t391;
t448 = -Ifges(6,2) * t389 + t551;
t636 = Ifges(6,6) * t433 + t323 * t448;
t649 = -t323 / 0.2e1;
t676 = 0.2e1 * Ifges(5,4) * t630 + t635 * t565 + t636 * t566 + Ifges(7,5) * t590 + Ifges(7,6) * t592 + t447 * t583 + (Ifges(5,2) + Ifges(6,3) + Ifges(7,3)) * t649;
t672 = t635 / 0.2e1;
t671 = t635 / 0.4e1;
t670 = t636 / 0.2e1;
t669 = t636 / 0.4e1;
t668 = pkin(4) * t640;
t666 = pkin(9) * t644;
t662 = t349 * t640;
t445 = -t389 * t95 + t391 * t96;
t148 = mrSges(7,2) * t323 - mrSges(7,3) * t420;
t151 = -mrSges(7,1) * t323 - mrSges(7,3) * t618;
t487 = t148 * t647 - t426 * t151 / 0.2e1;
t254 = mrSges(6,2) * t323 - mrSges(6,3) * t501;
t493 = t391 * t254;
t257 = -mrSges(6,1) * t323 - mrSges(6,3) * t500;
t497 = t389 * t257;
t606 = -m(6) / 0.2e1;
t659 = t445 * t606 - t493 / 0.2e1 + t497 / 0.2e1 - t487;
t657 = pkin(5) * t433 - pkin(9) * t643;
t580 = -t323 / 0.4e1;
t648 = t323 / 0.2e1;
t577 = t323 / 0.4e1;
t387 = t389 ^ 2;
t388 = t391 ^ 2;
t480 = t387 + t388;
t637 = mrSges(6,3) * t480;
t646 = Ifges(7,3) * t630;
t529 = t426 * mrSges(7,3);
t470 = -t529 / 0.2e1;
t639 = t470 * t618;
t638 = t618 * t570;
t559 = pkin(4) * t433;
t258 = -qJ(5) * t323 + t559;
t633 = -mrSges(7,1) / 0.2e1;
t550 = Ifges(7,4) * t618;
t86 = -Ifges(7,2) * t420 - Ifges(7,6) * t323 + t550;
t632 = -t86 / 0.2e1;
t370 = -mrSges(6,1) * t391 + mrSges(6,2) * t389;
t619 = t370 - mrSges(5,1);
t329 = Ifges(7,5) * t426 + Ifges(7,6) * t364;
t372 = Ifges(6,5) * t389 + Ifges(6,6) * t391;
t615 = t372 / 0.4e1 + t329 / 0.4e1;
t617 = (-Ifges(5,6) / 0.2e1 + t615) * t433;
t336 = (-pkin(9) - t348) * t389;
t384 = t391 * pkin(9);
t337 = t348 * t391 + t384;
t265 = t336 * t563 - t337 * t562;
t266 = t336 * t562 + t337 * t563;
t484 = -t265 * t426 + t266 * t364;
t369 = (-pkin(9) - qJ(5)) * t389;
t517 = qJ(5) * t391;
t371 = t384 + t517;
t334 = t369 * t563 - t371 * t562;
t335 = t369 * t562 + t371 * t563;
t482 = -t334 * t426 + t335 * t364;
t614 = t649 + t648;
t613 = Ifges(5,4) * t323 + t346 * mrSges(5,2) + Ifges(5,1) * t583 + (-Ifges(6,6) * t323 + t433 * t448) * t566 + (-Ifges(6,5) * t323 + t433 * t449) * t565;
t281 = t426 * t350;
t282 = t364 * t350;
t612 = (t281 * t426 + t282 * t364) * mrSges(7,3) + (t328 + t619) * t351 + (-mrSges(5,2) + t637) * t350;
t504 = t420 * t364;
t17 = (-t504 / 0.2e1 + t638) * mrSges(7,3) - t487 + t486;
t611 = t17 * qJD(1);
t327 = mrSges(7,1) * t426 + mrSges(7,2) * t364;
t64 = 0.2e1 * mrSges(7,1) * t590 + 0.2e1 * t592 * mrSges(7,2);
t610 = t327 * (qJD(2) + qJD(4)) + qJD(1) * t64;
t608 = 0.2e1 * qJD(4);
t607 = m(5) / 0.2e1;
t605 = m(6) / 0.2e1;
t604 = -m(7) / 0.2e1;
t603 = m(7) / 0.2e1;
t602 = -mrSges(7,2) / 0.2e1;
t601 = -mrSges(7,3) / 0.2e1;
t386 = t394 * pkin(2);
t347 = -pkin(3) * t363 + t386;
t203 = t258 + t347;
t97 = t391 * t203 - t698;
t70 = t657 + t97;
t98 = t389 * t203 + t695;
t75 = t98 - t666;
t44 = -t562 * t75 + t563 * t70;
t600 = -t44 / 0.2e1;
t111 = t391 * t258 - t698;
t72 = t111 + t657;
t112 = t389 * t258 + t695;
t82 = t112 - t666;
t52 = t562 * t72 + t563 * t82;
t599 = t52 / 0.2e1;
t224 = Ifges(7,4) * t420;
t89 = Ifges(7,1) * t618 - Ifges(7,5) * t323 - t224;
t598 = t89 / 0.2e1;
t597 = t690 / 0.2e1;
t596 = t654 / 0.2e1;
t595 = t148 / 0.2e1;
t594 = t688 / 0.2e1;
t589 = t656 / 0.2e1;
t588 = -t655 / 0.2e1;
t587 = t265 / 0.2e1;
t586 = -t266 / 0.2e1;
t585 = -t281 / 0.2e1;
t574 = -t334 / 0.2e1;
t573 = -t335 / 0.2e1;
t572 = t343 / 0.2e1;
t571 = -t364 / 0.2e1;
t569 = t370 / 0.2e1;
t567 = t380 / 0.2e1;
t554 = m(7) * qJD(3);
t547 = Ifges(7,5) * t645;
t545 = Ifges(7,6) * t642;
t544 = t111 * mrSges(6,3);
t543 = t112 * mrSges(6,3);
t542 = t688 * mrSges(5,1);
t537 = t420 * mrSges(7,2);
t534 = t618 * mrSges(7,1);
t319 = t433 * mrSges(5,1);
t45 = t562 * t70 + t563 * t75;
t455 = Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,3) / 0.2e1;
t3 = m(7) * (t42 * t44 + t43 * t45 + t700) + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t395 + (Ifges(3,1) - Ifges(3,2)) * t394) * t395 + t98 * t254 + t97 * t257 + t44 * t151 + t45 * t148 + t645 * t598 + (-t545 + t547) * t649 + (-t347 * mrSges(5,1) - t455 * t433 + t447 * t649 + t613) * t323 + t642 * t632 + m(6) * (t95 * t97 + t96 * t98 - t509) + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t394) * t394 + (t347 * mrSges(5,2) + Ifges(5,1) * t648 + t676) * t433 + (m(4) * t386 + t685) * t381 + (-mrSges(4,2) * t386 - Ifges(4,4) * t363) * t363 + (-mrSges(4,1) * t386 + Ifges(4,4) * t362 + (-Ifges(4,1) + Ifges(4,2)) * t363) * t362 + t693 + (m(5) * t347 + t319) * t346;
t533 = t3 * qJD(1);
t532 = t364 * mrSges(7,3);
t439 = Ifges(7,5) * t591 + Ifges(7,6) * t593;
t51 = -t562 * t82 + t563 * t72;
t4 = -(-(-t548 / 0.2e1 + t546 / 0.2e1) * t323 + t439 - t613) * t323 + t112 * t254 + t111 * t257 + t51 * t151 + t52 * t148 + (t346 * mrSges(5,1) - (-Ifges(5,1) / 0.2e1 + t455) * t323 + t676) * t433 + t89 * t591 + t86 * t593 + m(6) * (t111 * t95 + t112 * t96 - t509) + m(7) * (t42 * t51 + t43 * t52 + t700) + t693;
t524 = t4 * qJD(1);
t523 = t51 * t426;
t522 = t52 * t364;
t117 = t534 - t537;
t121 = -Ifges(7,2) * t618 - t224;
t122 = -Ifges(7,1) * t420 - t550;
t485 = -Ifges(7,5) * t420 - Ifges(7,6) * t618;
t9 = t485 * t649 - t43 * t151 + t143 * t117 + t42 * t148 + (t632 - t43 * mrSges(7,3) + t122 / 0.2e1) * t618 - (t121 / 0.2e1 - t42 * mrSges(7,3) + t598) * t420;
t521 = t9 * qJD(1);
t520 = t97 * t389;
t519 = t98 * t391;
t459 = t480 * t348;
t478 = m(4) * pkin(2) / 0.2e1;
t399 = (t323 * t351 - t350 * t433) * t607 + (t323 * t459 + t349 * t433) * t605 + (-t265 * t642 + t266 * t645 + t343 * t433) * t603 + (t362 * t390 + t363 * t392) * t478;
t401 = t347 * t607 + (t389 * t98 + t391 * t97) * t605 + (t364 * t44 + t426 * t45) * t603 + t653 * t647 + t654 * t570 + t389 * t589 + t655 * t565 + t394 * t478;
t11 = t319 + (-t328 / 0.2e1 - t370 / 0.2e1) * t433 + (-t570 * t642 + t571 * t645) * mrSges(7,3) + (mrSges(5,2) + (-t388 / 0.2e1 - t387 / 0.2e1) * mrSges(6,3)) * t323 - t399 + t401 + t685;
t516 = qJD(1) * t11;
t23 = -t618 * t151 + m(7) * (-t42 * t618 - t420 * t43) - t420 * t148 + (-t389 * t254 - t391 * t257 + m(6) * (-t389 * t96 - t391 * t95)) * t433;
t515 = qJD(1) * t23;
t507 = t206 * t433;
t10 = t645 * t148 - t642 * t151 + (t362 ^ 2 + t363 ^ 2) * mrSges(4,3) + (t323 * mrSges(5,3) + t493 - t497) * t323 + (mrSges(5,3) * t433 + t120 + t247) * t433 + m(7) * (t143 * t433 - t42 * t642 + t43 * t645) + m(6) * (t323 * t445 - t507) + m(5) * (t323 * t688 - t507) + m(4) * (-t338 * t362 + t363 * t634);
t513 = t10 * qJD(1);
t508 = t688 * t370;
t373 = Ifges(6,2) * t391 + t552;
t494 = t389 * t373;
t374 = Ifges(6,1) * t389 + t551;
t490 = t391 * t374;
t412 = (-t364 * t618 - t420 * t426) * t603 + m(6) * t480 * t630;
t477 = m(7) / 0.4e1 + m(6) / 0.4e1;
t62 = -0.2e1 * t433 * t477 + t412;
t489 = t62 * qJD(1);
t481 = Ifges(7,5) * t364 - Ifges(7,6) * t426;
t325 = t327 * qJD(6);
t145 = -t281 * t364 + t282 * t426;
t476 = t145 * t603;
t475 = qJD(4) * t476;
t471 = t532 / 0.2e1;
t466 = t577 + t580;
t330 = -Ifges(7,2) * t426 + t359;
t465 = t330 / 0.4e1 + t333 / 0.4e1;
t332 = Ifges(7,1) * t364 - t549;
t464 = -t331 / 0.4e1 + t332 / 0.4e1;
t460 = t480 * qJ(5);
t458 = t480 * t350;
t457 = t372 / 0.2e1 + t329 / 0.2e1 - Ifges(5,6);
t451 = 0.2e1 * t477 * t351;
t446 = t364 * t43 - t42 * t426;
t444 = t519 - t520;
t403 = (-t420 * t571 - t638) * mrSges(7,3) + t659;
t408 = t603 * t690 + t605 * t688 + t686;
t416 = m(7) * (-t265 * t618 - t266 * t420 + t446);
t13 = -t416 / 0.2e1 + t403 + t408;
t440 = (t364 ^ 2 + t426 ^ 2) * mrSges(7,3) + t637;
t71 = m(6) * t459 + m(7) * t484 + t440;
t443 = -qJD(1) * t13 + qJD(2) * t71;
t442 = -t111 * t389 + t112 * t391;
t441 = t691 - t699;
t438 = mrSges(7,1) * t585 + t282 * t602;
t436 = t328 * t583 + t433 * t569 - t470 * t642 + t471 * t645 + t637 * t648;
t402 = (t111 * t391 + t112 * t389) * t606 + (t364 * t51 + t426 * t52) * t604 + mrSges(5,1) * t630 + t653 * t571 - t426 * t596 + t656 * t566 + t391 * t588;
t404 = -t319 / 0.2e1 + (t323 * t460 - t559) * t605 + (-t334 * t642 + t335 * t645 + t380 * t433) * t603 + t436;
t19 = 0.2e1 * t649 * mrSges(5,2) + t402 + t404;
t429 = -t19 * qJD(1) + qJD(2) * t476;
t428 = -t420 * t471 - t639 - t659;
t427 = -Ifges(5,5) - t490 / 0.2e1 + t494 / 0.2e1;
t425 = t669 + t681 / 0.2e1 + t350 * t254 / 0.2e1;
t424 = t671 + t348 * t588 - t350 * t257 / 0.2e1;
t126 = m(6) * t460 + m(7) * t482 + t440;
t406 = m(6) * t594 + m(7) * t597 + t686;
t415 = m(7) * (-t334 * t618 - t335 * t420 + t446);
t15 = -t415 / 0.2e1 + t403 + t406;
t405 = (t459 + t460) * t606 + (t482 + t484) * t604 - t440;
t55 = t451 + t405;
t423 = -qJD(1) * t15 - qJD(2) * t55 + qJD(4) * t126;
t417 = (t330 / 0.2e1 + t333 / 0.2e1) * t364 - (t331 / 0.2e1 - t332 / 0.2e1) * t426;
t35 = t343 * t327 + t417;
t407 = (t121 / 0.4e1 + t89 / 0.4e1) * t364 - (t86 / 0.4e1 - t122 / 0.4e1) * t426 + t143 * t327 / 0.2e1 + t481 * t580;
t397 = -(t265 * t601 + t465) * t420 + (mrSges(7,3) * t586 + t464) * t618 + t148 * t587 + t151 * t586 + t117 * t572 + t407;
t410 = -t547 / 0.2e1 + t545 / 0.2e1 + t646 + mrSges(7,1) * t600 + t45 * mrSges(7,2) / 0.2e1;
t6 = t397 + t410;
t421 = -t6 * qJD(1) - t35 * qJD(2);
t396 = (t247 / 0.2e1 + t120 / 0.2e1) * t351 + t617 + (t143 * t351 + t265 * t51 + t266 * t52 - t281 * t42 + t282 * t43 + t697) * t603 + t653 * t587 + t266 * t596 + t151 * t585 + t282 * t595 + t675 * t572 + t662 / 0.2e1;
t400 = (qJ(5) * t444 - t692) * t606 + (t334 * t44 + t335 * t45 + t696) * t604 + t668 / 0.2e1 + t653 * t574 + t654 * t573 - t687 / 0.2e1;
t2 = (t597 - t690 / 0.2e1) * t328 - t617 + t614 * Ifges(5,5) + (t348 * t442 + t350 * t445 + t441) * t605 + (-t635 / 0.4e1 + t518 / 0.2e1 + t466 * t373 + (-t111 / 0.2e1 + t97 / 0.2e1) * mrSges(6,3) + t424) * t389 + (-t636 / 0.4e1 - t682 / 0.2e1 - t466 * t374 + (t112 / 0.2e1 - t98 / 0.2e1) * mrSges(6,3) + t425) * t391 + (-(t51 / 0.2e1 + t600) * t426 + (t599 - t45 / 0.2e1) * t364) * mrSges(7,3) + t396 + t400 + t619 * (t594 - t688 / 0.2e1);
t34 = m(7) * (-t265 * t281 + t266 * t282 + t343 * t351) + m(6) * (t348 * t458 + t349 * t351) + t612;
t419 = -t2 * qJD(1) - t34 * qJD(2) - t145 * t554 / 0.2e1;
t409 = (t572 + t567) * t327 + t417;
t31 = t409 - t438;
t54 = t380 * t327 + t417;
t398 = (mrSges(7,3) * t573 + t464) * t618 - (mrSges(7,3) * t574 + t465) * t420 + t334 * t595 + t151 * t573 + t117 * t567 + t407;
t411 = mrSges(7,2) * t599 + t51 * t633 - t439 + t646;
t8 = t398 + t411;
t414 = -t8 * qJD(1) - t31 * qJD(2) - t54 * qJD(4);
t326 = t327 * qJD(5);
t65 = -t537 / 0.2e1 + t534 / 0.2e1 - t420 * t602 + t618 * t633;
t61 = t412 + (m(6) + m(7)) * t583;
t56 = t451 - t405;
t32 = t409 + t438;
t20 = mrSges(5,2) * t614 - t402 + t404;
t18 = -t504 * t601 + t486 + t487 + t639;
t16 = t415 / 0.2e1 + t406 + t428;
t14 = t416 / 0.2e1 + t408 + t428;
t12 = t399 + t401 + t436;
t7 = t398 - t411;
t5 = t397 - t410;
t1 = t615 * t433 + t702 + (-t523 / 0.2e1 + t522 / 0.2e1) * mrSges(7,3) + Ifges(5,6) * t630 + (t569 - mrSges(5,1) / 0.2e1) * t688 + 0.2e1 * t648 * Ifges(5,5) + t508 / 0.2e1 + t494 * t580 + t396 - t542 / 0.2e1 + (t669 + t374 * t577 + t543 / 0.2e1 + (t112 * t348 + t350 * t96) * t605 + t425) * t391 + (t671 + (-t111 * t348 - t350 * t95) * t605 - t544 / 0.2e1 + t373 * t580 + t424) * t389 + t44 * t470 + t45 * t471 + (t519 / 0.2e1 - t520 / 0.2e1) * mrSges(6,3) - t400 + t518 * t566 + t490 * t577 + t517 * t589 + t441 * t605;
t21 = [qJD(2) * t3 + qJD(3) * t10 + qJD(4) * t4 + qJD(5) * t23 + qJD(6) * t9, t12 * qJD(3) + t1 * qJD(4) + t14 * qJD(5) + t5 * qJD(6) + t533 + (0.2e1 * (-t350 * t688 + t699) * t607 + 0.2e1 * (t265 * t44 + t266 * t45 + t697) * t603 + t343 * t675 + t702 + (-t351 * mrSges(5,3) + t457) * t433 - t44 * t529 + t45 * t532 + t662 + (-mrSges(3,1) * t395 + mrSges(3,2) * t394) * pkin(7) + (-t362 * t560 + t363 * t561) * mrSges(4,3) - Ifges(3,6) * t394 + Ifges(3,5) * t395 + Ifges(4,6) * t363 + Ifges(4,5) * t362 + t338 * mrSges(4,1) + (-t350 * mrSges(5,3) - t427) * t323 + 0.2e1 * (t348 * t444 + t691) * t605 - t634 * mrSges(4,2) + 0.2e1 * (t338 * t392 + t390 * t634) * t478 + t508 - t542 + t266 * t654 + t265 * t653 + (-t97 * mrSges(6,3) - t348 * t655 + t672) * t389 + (t98 * mrSges(6,3) + t670 + t681) * t391) * qJD(2), t513 + t12 * qJD(2) + (-t364 * t642 + t426 * t645) * t554 + t20 * qJD(4) + t61 * qJD(5) + t18 * qJD(6), t524 + t1 * qJD(2) + t20 * qJD(3) + t16 * qJD(5) + t7 * qJD(6) + ((qJ(5) * t442 - t692) * t605 + (t334 * t51 + t335 * t52 + t696) * t603) * t608 + (t687 + t334 * t653 + t335 * t654 - t668 + (t543 + t670 + t682) * t391 + (-t544 + t672 - t518) * t389 + t457 * t433 - t427 * t323 + t619 * t688 + (-t523 + t522) * mrSges(7,3) + t702) * qJD(4), qJD(2) * t14 + qJD(3) * t61 + qJD(4) * t16 + qJD(6) * t65 + t515, t521 + t5 * qJD(2) + t18 * qJD(3) + t7 * qJD(4) + t65 * qJD(5) + (-mrSges(7,1) * t43 - mrSges(7,2) * t42 + t485) * qJD(6); -qJD(3) * t11 + qJD(4) * t2 - qJD(5) * t13 + qJD(6) * t6 - t533, qJD(4) * t34 + qJD(5) * t71 + qJD(6) * t35, t475 - t516, t56 * qJD(5) + t32 * qJD(6) + ((-t281 * t334 + t282 * t335 + t351 * t380) * t603 + (-pkin(4) * t351 + qJ(5) * t458) * t605) * t608 - t419 + t612 * qJD(4), qJD(4) * t56 + t443, t32 * qJD(4) + (-mrSges(7,1) * t266 - mrSges(7,2) * t265 + t481) * qJD(6) - t421; qJD(2) * t11 - qJD(4) * t19 + qJD(5) * t62 - qJD(6) * t17 - t513, t475 + t516, 0, t429, t489, -t325 - t611; -qJD(2) * t2 + qJD(3) * t19 - qJD(5) * t15 + qJD(6) * t8 - t524, -t55 * qJD(5) + t31 * qJD(6) + t419, -t429, qJD(5) * t126 + qJD(6) * t54, t423 (-mrSges(7,1) * t335 - mrSges(7,2) * t334 + t481) * qJD(6) - t414; qJD(2) * t13 - qJD(3) * t62 + qJD(4) * t15 + qJD(6) * t64 - t515, qJD(4) * t55 + t325 - t443, -t489, -t423 + t325, 0, t610; -qJD(2) * t6 + qJD(3) * t17 - qJD(4) * t8 - qJD(5) * t64 - t521, -t31 * qJD(4) - t326 + t421, t611, -t326 + t414, -t610, 0;];
Cq  = t21;
