% Calculate vector of inverse dynamics joint torques for
% S5RPRRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:28
% DurationCPUTime: 16.73s
% Computational Cost: add. (10355->575), mult. (22694->679), div. (0->0), fcn. (23688->6), ass. (0->295)
t331 = sin(qJ(4));
t544 = sin(qJ(3));
t545 = sin(qJ(1));
t546 = cos(qJ(3));
t547 = cos(qJ(1));
t272 = -t544 * t545 - t546 * t547;
t273 = t544 * t547 - t545 * t546;
t332 = cos(qJ(4));
t519 = Icges(6,5) * t331;
t393 = Icges(6,1) * t332 + t519;
t135 = Icges(6,4) * t272 + t273 * t393;
t525 = Icges(5,4) * t331;
t395 = Icges(5,1) * t332 - t525;
t139 = Icges(5,5) * t272 + t273 * t395;
t668 = t135 + t139;
t676 = t331 * t668;
t330 = qJD(1) - qJD(3);
t389 = Icges(6,4) * t331 - Icges(6,6) * t332;
t616 = t389 * t273;
t387 = Icges(5,5) * t331 + Icges(5,6) * t332;
t618 = t387 * t273;
t518 = Icges(6,5) * t332;
t289 = -Icges(6,1) * t331 + t518;
t385 = -Icges(6,3) * t332 + t519;
t377 = t289 * t332 - t331 * t385;
t391 = Icges(5,2) * t332 + t525;
t524 = Icges(5,4) * t332;
t394 = Icges(5,1) * t331 + t524;
t665 = -t331 * t391 + t332 * t394 - t377;
t595 = t272 * t665 - t616 - t618;
t622 = t595 * t330;
t617 = t389 * t272;
t619 = t387 * t272;
t596 = t273 * t665 + t617 + t619;
t623 = t596 * t330;
t388 = Icges(5,5) * t332 - Icges(5,6) * t331;
t122 = -Icges(5,3) * t273 + t272 * t388;
t390 = Icges(6,4) * t332 + Icges(6,6) * t331;
t126 = -Icges(6,2) * t273 + t272 * t390;
t637 = t122 + t126;
t123 = Icges(5,3) * t272 + t273 * t388;
t127 = Icges(6,2) * t272 + t273 * t390;
t636 = t123 + t127;
t134 = -Icges(6,4) * t273 + t272 * t393;
t138 = -Icges(5,5) * t273 + t272 * t395;
t584 = t134 + t138;
t497 = t273 * t332;
t225 = Icges(6,5) * t497;
t498 = t273 * t331;
t513 = Icges(6,6) * t272;
t120 = -Icges(6,3) * t498 - t225 - t513;
t675 = -t331 * t120 + t139 * t332;
t674 = t331 * t584;
t673 = Icges(5,1) + Icges(6,1);
t672 = Icges(5,4) - Icges(6,5);
t657 = Icges(6,4) + Icges(5,5);
t671 = -Icges(5,2) - Icges(6,3);
t655 = Icges(5,6) - Icges(6,6);
t392 = -Icges(5,2) * t331 + t524;
t131 = Icges(5,6) * t272 + t273 * t392;
t670 = t120 + t131;
t501 = t272 * t332;
t226 = Icges(6,5) * t501;
t502 = t272 * t331;
t512 = Icges(6,6) * t273;
t121 = -Icges(6,3) * t502 - t226 + t512;
t130 = -Icges(5,6) * t273 + t272 * t392;
t669 = t130 + t121;
t667 = -t121 * t498 + t637 * t272 + t497 * t584;
t666 = t135 * t501 + t272 * t675 - t636 * t273;
t573 = t135 * t497 + t636 * t272 + t273 * t675;
t198 = t330 * t272;
t199 = t330 * t273;
t103 = -pkin(3) * t198 + pkin(7) * t199;
t151 = qJD(4) * t199 - qJDD(4) * t272;
t329 = qJDD(1) - qJDD(3);
t333 = qJD(1) ^ 2;
t323 = t547 * qJ(2);
t446 = t545 * pkin(1);
t296 = t446 - t323;
t321 = qJD(2) * t547;
t459 = pkin(1) * t547 + qJ(2) * t545;
t360 = -qJDD(1) * t296 + qJDD(2) * t545 + (-qJD(1) * t459 + 0.2e1 * t321) * qJD(1);
t335 = (-qJDD(1) * t545 - t333 * t547) * pkin(2) + t360;
t403 = rSges(5,1) * t331 + rSges(5,2) * t332;
t301 = -rSges(5,1) * t332 + rSges(5,2) * t331;
t276 = t301 * qJD(4);
t456 = qJD(4) * t276;
t147 = -rSges(5,1) * t497 + rSges(5,2) * t498 - rSges(5,3) * t272;
t419 = pkin(3) * t273 + pkin(7) * t272;
t478 = -t147 + t419;
t454 = qJD(4) * t331;
t427 = t273 * t454;
t507 = t198 * t332;
t374 = t427 - t507;
t369 = rSges(5,1) * t374 + rSges(5,3) * t199;
t453 = qJD(4) * t332;
t375 = t198 * t331 + t273 * t453;
t593 = rSges(5,2) * t375;
t56 = t369 + t593;
t16 = -t272 * t456 - t151 * t403 + (-t103 - t56) * t330 + t478 * t329 + t335;
t664 = -g(1) + t16;
t319 = qJD(5) * t332;
t548 = rSges(6,1) + pkin(4);
t621 = rSges(6,3) + qJ(5);
t463 = t331 * t621 + t332 * t548;
t473 = qJD(4) * t463 - t319;
t349 = -qJDD(5) * t331 + (-t319 + t473) * qJD(4);
t452 = qJD(5) * t331;
t464 = -t331 * t548 + t332 * t621;
t217 = t273 * t452;
t575 = t199 * rSges(6,2) - t375 * t621 - t507 * t548 - t217;
t540 = t427 * t548 + t575;
t480 = rSges(6,2) * t272 + t497 * t548 + t498 * t621;
t581 = -t419 - t480;
t10 = t199 * t452 + t464 * t151 + (-t103 - t540) * t330 - t581 * t329 + t349 * t272 + t335;
t663 = t10 - g(1);
t601 = -t130 * t498 + t667;
t493 = t131 * t331;
t600 = -t272 * t493 + t666;
t107 = t273 * t126;
t38 = -t121 * t502 + t134 * t501 - t107;
t491 = -t122 * t273 + t138 * t501;
t40 = -t130 * t502 + t491;
t599 = t38 + t40;
t662 = -t199 * t655 - t374 * t672 + t375 * t671;
t429 = t272 * t454;
t505 = t199 * t332;
t372 = t429 + t505;
t373 = -t199 * t331 + t272 * t453;
t661 = -t198 * t655 - t372 * t672 + t373 * t671;
t660 = t199 * t657 + t374 * t673 + t375 * t672;
t659 = -t198 * t657 - t372 * t673 - t373 * t672;
t602 = t273 * t493 - t573;
t598 = t332 * t670 + t676;
t597 = t332 * t669 + t674;
t658 = t332 * t584;
t656 = Icges(6,2) + Icges(5,3);
t654 = (t385 - t391) * t332 + (t289 - t394) * t331;
t386 = Icges(6,3) * t331 + t518;
t653 = (t386 - t392) * qJD(4);
t652 = (t393 + t395) * qJD(4);
t474 = pkin(3) * t199 + pkin(7) * t198;
t58 = rSges(5,1) * t372 + rSges(5,2) * t373 + rSges(5,3) * t198;
t590 = t330 * t147;
t651 = -t474 - t58 - t590;
t101 = -rSges(4,1) * t198 - rSges(4,2) * t199;
t470 = rSges(4,1) * t273 - rSges(4,2) * t272;
t41 = -t101 * t330 + t329 * t470 + t335;
t648 = -g(1) + t41;
t150 = qJD(4) * t198 + qJDD(4) * t273;
t327 = t547 * pkin(2);
t320 = qJD(2) * t545;
t424 = qJD(1) * t545;
t425 = qJD(1) * t547;
t462 = qJ(2) * t425 + t320;
t365 = qJDD(1) * t459 - qJDD(2) * t547 + (-pkin(1) * t424 + t320 + t462) * qJD(1);
t445 = t545 * pkin(2);
t339 = qJDD(1) * t327 - t333 * t445 + t365;
t570 = pkin(3) * t272 - pkin(7) * t273;
t336 = -t329 * t570 + t330 * t474 + t339;
t218 = t272 * t452;
t574 = -t198 * rSges(6,2) + t373 * t621 - t505 * t548 + t218;
t539 = t429 * t548 - t574;
t563 = t273 * rSges(6,2) - t501 * t548 - t502 * t621;
t11 = -t150 * t464 - t198 * t452 + t273 * t349 + t329 * t563 + t330 * t539 + t336;
t647 = t11 - g(2);
t169 = t403 * t273;
t189 = t330 * t419;
t359 = t320 + (-t445 - t296) * qJD(1);
t352 = t189 + t359;
t455 = qJD(4) * t403;
t61 = t272 * t455 + t352 - t590;
t646 = t169 * t61;
t471 = rSges(5,1) * t501 - rSges(5,3) * t273;
t149 = rSges(5,2) * t502 - t471;
t614 = t149 - t570;
t416 = qJD(4) * t464;
t567 = t330 * t480 - t218;
t644 = -t272 * t416 + t567;
t612 = t570 - t563;
t568 = -t330 * t612 - t217;
t643 = -t273 * t416 + t568;
t102 = rSges(4,1) * t199 - rSges(4,2) * t198;
t200 = rSges(4,1) * t272 + rSges(4,2) * t273;
t42 = t102 * t330 - t200 * t329 + t339;
t641 = -g(2) + t42;
t640 = t199 * t656 + t374 * t657 + t375 * t655;
t639 = t198 * t656 + t372 * t657 + t373 * t655;
t638 = t330 * t470;
t635 = (t388 + t390) * qJD(4);
t634 = qJD(4) * t654 + t331 * t653 + t332 * t652;
t633 = -qJD(4) * t597 - t331 * t661 + t332 * t659;
t632 = qJD(4) * t598 + t331 * t662 + t332 * t660;
t367 = -t446 - t445;
t357 = t323 + t367;
t631 = t357 + t419;
t509 = t130 * t331;
t629 = -t121 * t331 - t509 + t658;
t628 = t135 * t332 - t493 + t675;
t627 = t387 + t389;
t626 = -t272 * t600 + t599 * t273;
t625 = t602 * t272 + t273 * t601;
t624 = qJD(1) * t296 - t320 + t462;
t457 = qJD(4) * t273;
t432 = t327 + t459;
t358 = qJD(1) * t432 - t321;
t620 = t200 * t330;
t615 = t272 * t386 - t130 - t512;
t587 = t273 * t386 - t131 + t513;
t613 = pkin(2) * t424 + qJD(1) * t367 + t624;
t458 = qJD(4) * t272;
t609 = qJD(4) * t625 + t623;
t608 = qJD(4) * t626 + t622;
t607 = -qJD(4) * t628 + t331 * t660 - t332 * t662;
t606 = qJD(4) * t629 + t331 * t659 + t332 * t661;
t605 = (-t130 * t273 - t131 * t272) * t331 + t666 + t667;
t604 = t198 * t665 - t199 * t627 + t272 * t635 + t273 * t634;
t603 = -t198 * t627 - t199 * t665 + t272 * t634 - t273 * t635;
t589 = t454 * t548;
t580 = t621 * t497;
t579 = t621 * t501;
t577 = (-t637 * t198 - t629 * t199 + t639 * t273) * t273 + (t632 * t272 + t628 * t199 + t636 * t198 + (t633 - t640) * t273) * t272;
t576 = (t198 * t629 - t199 * t637 + t273 * t633) * t273 + (t640 * t272 + t636 * t199 - t628 * t198 + (t632 - t639) * t273) * t272;
t460 = rSges(3,1) * t547 + rSges(3,3) * t545;
t572 = t459 + t460;
t571 = t457 * t598 - t458 * t597;
t564 = t613 - t189;
t484 = t273 * t391 - t139;
t488 = -t273 * t394 - t131;
t562 = t331 * t484 + t332 * t488;
t486 = -t273 * t385 - t135;
t490 = Icges(6,1) * t498 + t120 - t225;
t561 = t331 * t486 - t332 * t490;
t560 = t272 ^ 2;
t550 = -t330 / 0.2e1;
t543 = g(3) * t332;
t32 = t358 + t643;
t534 = t272 * t32;
t31 = t352 + t644;
t531 = t273 * t31;
t28 = t319 + (t272 * t563 - t273 * t480) * qJD(4);
t530 = t28 * t331;
t62 = t273 * t455 + t330 * t614 + t358;
t529 = t62 * t403;
t496 = t388 * t330;
t495 = t390 * t330;
t489 = -Icges(6,1) * t502 - t121 + t226;
t487 = -t272 * t394 - t130;
t485 = -t272 * t385 - t134;
t483 = t272 * t391 - t138;
t476 = t498 * t548 - t580;
t475 = t502 * t548 - t579;
t468 = t385 + t393;
t467 = -t386 - t289;
t466 = -t391 + t395;
t465 = -t392 - t394;
t449 = -t547 / 0.2e1;
t448 = t545 / 0.2e1;
t436 = t545 * rSges(3,1);
t431 = t272 * t548;
t430 = t273 * t548;
t422 = t458 / 0.2e1;
t421 = -t457 / 0.2e1;
t420 = t457 / 0.2e1;
t418 = t31 * t464;
t417 = t32 * t464;
t402 = -t272 * t31 - t273 * t32;
t397 = t272 * t58 + t273 * t56;
t396 = -t272 * t61 - t273 * t62;
t378 = t147 * t273 + t149 * t272;
t371 = t402 * t332;
t366 = -t436 - t446;
t304 = rSges(2,1) * t547 - rSges(2,2) * t545;
t298 = rSges(2,1) * t545 + rSges(2,2) * t547;
t356 = t331 * t485 + t332 * t489;
t355 = t331 * t483 + t332 * t487;
t354 = t369 + t103;
t351 = (-t331 * t467 + t332 * t468) * t330;
t350 = (t331 * t465 + t332 * t466) * t330;
t347 = t103 + t575;
t346 = -t474 + t574;
t334 = t609 * t420 + (qJD(4) * t377 + t391 * t454 + (-qJD(4) * t394 + t653) * t332 - t652 * t331) * t330 + (-Icges(4,3) + t654) * t329 - (t595 + t597) * t150 / 0.2e1 - (t596 + t598) * t151 / 0.2e1 + (t603 + t606) * t421 + (t604 - t607 + t608) * t422;
t325 = t547 * rSges(3,3);
t317 = rSges(3,3) * t425;
t297 = t436 - t325;
t173 = t403 * t272;
t117 = t358 - t620;
t116 = t359 + t638;
t87 = qJDD(1) * t460 + qJD(1) * (-rSges(3,1) * t424 + t317) + t365;
t86 = -qJDD(1) * t297 - t333 * t460 + t360;
t63 = t378 * qJD(4);
t17 = t149 * t329 + t150 * t403 - t273 * t456 + t330 * t58 + t336;
t1 = qJDD(5) * t332 - t563 * t151 - t480 * t150 + (t539 * t272 + t540 * t273 - t452) * qJD(4);
t2 = [-t334 - m(2) * (-g(1) * t298 + g(2) * t304) + (((t332 * t587 - t676) * t273 + t597 * t272) * qJD(4) + t571) * t550 + (((t601 - t605) * t272 + (t491 + t38 + (-t331 * t587 - t332 * t668) * t273 + (-t509 - t636) * t272 - t602) * t273) * qJD(4) + t622) * t422 + (-t623 + ((t40 + (t123 + t509) * t272 - t491) * t272 + t560 * t127 + ((t636 - t629) * t273 + ((-t587 - t670) * t331 - t637) * t272 + t600) * t273) * qJD(4)) * t421 + (m(2) * (t298 ^ 2 + t304 ^ 2) + Icges(2,3) + Icges(3,2)) * qJDD(1) + (-t418 * t457 - t531 * t589 + t647 * (-t612 + t432) + (t431 * t454 - t346 + t564 - t644) * t32 + (pkin(2) * t425 - t347 + (t459 - t432) * qJD(1) + t568) * t31 + t663 * (t631 + t480)) * m(6) + (-t529 * t458 + (-g(2) + t17) * (t432 + t614) + (-t354 - t358 - t593) * t61 + (t564 + t61 - t651) * t62 + t664 * (-t147 + t631)) * m(5) + (t641 * (-t200 + t432) + (-t101 - t358) * t116 + (t102 + t116 + t613 - t638) * t117 + t648 * (t357 + t470)) * m(4) + ((-g(2) + t87) * t572 + (-g(1) + t86) * (t323 + t325 + t366) + (t317 + (t297 + t366) * qJD(1) + t624) * (qJD(1) * t572 - t321)) * m(3); (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t545 - g(2) * t547) + 0.2e1 * (t10 * t448 + t11 * t449) * m(6) + 0.2e1 * (t16 * t448 + t17 * t449) * m(5) + 0.2e1 * (t41 * t448 + t42 * t449) * m(4) + 0.2e1 * (t448 * t86 + t449 * t87) * m(3); t334 + ((-t598 * t273 + (-t332 * t615 + t674) * t272) * qJD(4) + t571) * t550 + (-t622 + ((t107 + (t122 - t493) * t273 + t573 + t602) * t273 + ((t637 + t628) * t272 + ((t615 + t669) * t331 - t636) * t273 - t601) * t272) * qJD(4)) * t422 + (((-t600 + t605) * t273 + ((t331 * t615 + t658) * t272 + (t493 - t637) * t273 - t573 - t599) * t272) * qJD(4) + t623) * t421 + (-t534 * t589 - t417 * t458 + t663 * t581 + (t430 * t454 + t347 - t643) * t31 + (t189 + t346 + t567) * t32 + t647 * t612) * m(6) + (-(-t272 * t529 + t646) * qJD(4) + t61 * t354 + t17 * (t570 + t471) + (-t17 * t502 + t375 * t61) * rSges(5,2) - t664 * t478 + (-t330 * t61 + g(2)) * t614 + (t189 + t651) * t62) * m(5) + (-t102 * t117 + (t101 + t620) * t116 - (-t117 * t330 + t648) * t470 + t641 * t200) * m(4); t626 * t150 / 0.2e1 + t625 * t151 / 0.2e1 + t608 * t198 / 0.2e1 + t609 * t199 / 0.2e1 - (qJD(4) * t576 + t150 * t601 - t151 * t602 + t329 * t596 + t330 * t604) * t272 / 0.2e1 + (qJD(4) * t577 + t150 * t599 + t151 * t600 + t329 * t595 + t330 * t603) * t273 / 0.2e1 + (-t272 * t598 + t273 * t597) * t329 / 0.2e1 + (((-t465 + t467) * t332 + (t466 + t468) * t331) * t330 + (((-t483 - t485) * t273 + (t484 + t486) * t272) * t332 + ((t487 + t489) * t273 + (-t488 + t490) * t272) * t331) * qJD(4)) * t550 + (t198 * t597 + t199 * t598 + t272 * t607 + t273 * t606) * t330 / 0.2e1 - (t198 * t601 - t199 * t602 + t576) * t458 / 0.2e1 + ((t458 * t618 + t496) * t272 + (t350 + (t355 * t273 + (-t619 - t562) * t272) * qJD(4)) * t273 + (t458 * t616 + t495) * t272 + (t351 + (t356 * t273 + (-t617 - t561) * t272) * qJD(4)) * t273) * t422 + ((t457 * t619 - t496) * t273 + (t350 + (-t562 * t272 + (-t618 + t355) * t273) * qJD(4)) * t272 + (t457 * t617 - t495) * t273 + (t351 + (-t561 * t272 + (-t616 + t356) * t273) * qJD(4)) * t272) * t421 + (t198 * t599 + t199 * t600 + t577) * t420 + (g(1) * t579 + g(2) * t580 + t548 * t543 - (g(1) * t431 + g(2) * t430 - g(3) * t621) * t331 + (-t28 * t563 + t418) * t199 - (t28 * t480 + t417) * t198 + (-t1 * t480 - t11 * t464 + t28 * t540 + t32 * t473) * t273 + (t1 * t563 - t10 * t464 + t28 * t539 + t31 * t473) * t272 - (t371 - t530) * qJD(5) - (-t31 * t476 + t32 * t475) * t330 - ((t28 * t476 + t32 * t463) * t273 + (t28 * t475 + t31 * t463) * t272) * qJD(4)) * m(6) + (-g(1) * t173 - g(2) * t169 - g(3) * t301 - (t173 * t62 - t646) * t330 - (t63 * (t169 * t273 + t173 * t272) + t396 * t301) * qJD(4) + (qJD(4) * t397 + t147 * t150 - t149 * t151) * t378 + t63 * (t147 * t198 - t149 * t199 + t397) + t396 * t276 - (-t16 * t272 - t17 * t273 - t198 * t62 + t199 * t61) * t403) * m(5); ((qJD(4) * t402 + t1) * t332 - (t371 + (-t273 ^ 2 - t560) * t530) * qJD(4) - t543 + (-qJD(4) * t28 - t198 * t32 + t199 * t31 - (t531 - t534) * t330 - t647 * t273 - t663 * t272) * t331) * m(6);];
tau = t2;
