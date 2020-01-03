% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:51
% EndTime: 2020-01-03 12:09:07
% DurationCPUTime: 9.36s
% Computational Cost: add. (49536->530), mult. (32605->653), div. (0->0), fcn. (29520->10), ass. (0->337)
t499 = qJD(1) + qJD(2);
t451 = qJ(1) + qJ(2);
t444 = sin(t451);
t450 = qJ(3) + pkin(9);
t443 = qJ(5) + t450;
t433 = sin(t443);
t434 = cos(t443);
t380 = rSges(6,1) * t433 + rSges(6,2) * t434;
t441 = sin(t450);
t453 = sin(qJ(3));
t612 = t453 * pkin(3);
t405 = -pkin(4) * t441 - t612;
t470 = t380 - t405;
t275 = t470 * t444;
t664 = m(6) / 0.2e1;
t445 = cos(t451);
t277 = t470 * t445;
t709 = t277 * t445;
t665 = m(5) / 0.2e1;
t442 = cos(t450);
t496 = rSges(5,1) * t441 + rSges(5,2) * t442 + t612;
t694 = t496 * t445;
t695 = t496 * t444;
t710 = (t444 * t695 + t445 * t694) * t665;
t530 = (-t275 * t444 - t709) * t664 - t710;
t332 = t380 * t444;
t379 = t444 * t405;
t279 = t379 - t332;
t531 = (-t279 * t444 + t709) * t664 + t710;
t51 = t531 - t530;
t712 = t499 * t51;
t455 = cos(qJ(3));
t447 = t455 * pkin(3);
t438 = t447 + pkin(2);
t613 = pkin(4) * t442;
t404 = t438 + t613;
t373 = t445 * t404;
t452 = -qJ(4) - pkin(7);
t449 = -pkin(8) + t452;
t558 = t434 * t445;
t562 = t433 * t445;
t488 = rSges(6,1) * t558 - rSges(6,2) * t562;
t256 = t373 + (rSges(6,3) - t449) * t444 + t488;
t448 = cos(qJ(1)) * pkin(1);
t247 = t448 + t256;
t238 = t247 * t444;
t245 = t256 * t444;
t563 = t433 * t444;
t506 = -rSges(6,2) * t563 - t445 * rSges(6,3);
t608 = rSges(6,1) * t434;
t255 = t445 * t449 + (t404 + t608) * t444 + t506;
t614 = sin(qJ(1)) * pkin(1);
t246 = t255 + t614;
t411 = t445 * t438;
t550 = t442 * t445;
t554 = t441 * t445;
t489 = rSges(5,1) * t550 - rSges(5,2) * t554;
t269 = t411 + (rSges(5,3) - t452) * t444 + t489;
t265 = t448 + t269;
t257 = t265 * t444;
t263 = t269 * t444;
t555 = t441 * t444;
t505 = -rSges(5,2) * t555 - t445 * rSges(5,3);
t542 = t445 * t452;
t609 = rSges(5,1) * t442;
t268 = t542 + (t438 + t609) * t444 + t505;
t264 = t268 + t614;
t602 = (t257 + t263 + (-t264 - t268) * t445) * t665 + (t238 + t245 + (-t246 - t255) * t445) * t664;
t525 = -t246 + t255;
t706 = -t264 + t268;
t603 = (t706 * t445 + t257 - t263) * t665 + (t445 * t525 + t238 - t245) * t664;
t16 = t603 - t602;
t711 = t16 * qJD(1);
t593 = Icges(6,4) * t433;
t378 = Icges(6,1) * t434 - t593;
t304 = -Icges(6,5) * t445 + t378 * t444;
t559 = t434 * t444;
t272 = t304 * t559;
t301 = Icges(6,5) * t558 - Icges(6,6) * t562 + Icges(6,3) * t444;
t401 = Icges(6,4) * t562;
t305 = Icges(6,1) * t558 + Icges(6,5) * t444 - t401;
t303 = Icges(6,4) * t558 - Icges(6,2) * t562 + Icges(6,6) * t444;
t581 = t303 * t433;
t477 = t305 * t434 - t581;
t708 = -t301 * t444 - t477 * t445 - t272;
t707 = -t247 + t256;
t374 = Icges(6,5) * t434 - Icges(6,6) * t433;
t572 = t374 * t444;
t300 = -Icges(6,3) * t445 + t572;
t705 = t300 * t444 + t304 * t558;
t427 = Icges(6,4) * t434;
t376 = -Icges(6,2) * t433 + t427;
t684 = Icges(6,1) * t433 + t427;
t704 = t376 + t684;
t428 = Icges(5,4) * t442;
t392 = -Icges(5,2) * t441 + t428;
t683 = Icges(5,1) * t441 + t428;
t703 = t392 + t683;
t446 = Icges(4,4) * t455;
t414 = -Icges(4,2) * t453 + t446;
t682 = Icges(4,1) * t453 + t446;
t702 = t414 + t682;
t701 = Icges(4,5) * t453 + Icges(5,5) * t441 + Icges(4,6) * t455 + Icges(5,6) * t442;
t666 = m(4) / 0.2e1;
t647 = -t444 / 0.2e1;
t646 = t444 / 0.2e1;
t645 = -t445 / 0.2e1;
t642 = m(3) * (-t448 * (rSges(3,1) * t444 + rSges(3,2) * t445) + t614 * (t445 * rSges(3,1) - rSges(3,2) * t444));
t610 = rSges(4,1) * t455;
t419 = -rSges(4,2) * t453 + t610;
t696 = t419 * t666;
t439 = t444 ^ 2;
t440 = t445 ^ 2;
t500 = t439 + t440;
t644 = t445 / 0.2e1;
t693 = t644 + t645;
t594 = Icges(5,4) * t441;
t391 = Icges(5,2) * t442 + t594;
t394 = Icges(5,1) * t442 - t594;
t595 = Icges(4,4) * t453;
t413 = Icges(4,2) * t455 + t595;
t416 = Icges(4,1) * t455 - t595;
t692 = t702 * t453 + (t413 - t416) * t455 + t703 * t441 + (t391 - t394) * t442;
t690 = m(6) * t380;
t312 = t444 * t332;
t333 = t380 * t445;
t492 = t445 * t333 + t312;
t381 = -rSges(6,2) * t433 + t608;
t569 = t381 * t445;
t570 = t381 * t444;
t479 = Icges(6,5) * t433 + Icges(6,6) * t434;
t326 = t444 * t479;
t327 = t479 * t445;
t519 = Icges(6,2) * t558 - t305 + t401;
t375 = Icges(6,2) * t434 + t593;
t520 = -t375 * t444 + t304;
t521 = t445 * t684 + t303;
t302 = -Icges(6,6) * t445 + t376 * t444;
t522 = -t444 * t684 - t302;
t680 = t433 * (-t444 * t519 - t445 * t520) + t434 * (t444 * t521 + t445 * t522);
t611 = (-t440 * t326 + (t445 * t327 - t680) * t444) * t645 + (t439 * t327 + (-t444 * t326 + t680) * t445) * t647;
t436 = t445 * pkin(7);
t288 = t444 * (t542 + t436 + (-pkin(2) + t438) * t444);
t293 = t444 * (rSges(6,1) * t559 + t506);
t501 = t445 * pkin(2) + t444 * pkin(7);
t306 = t444 * t452 - t411 + t501;
t309 = rSges(6,3) * t444 + t488;
t98 = t288 + t293 + (t404 - t438) * t439 + (-t306 + t373 - t411 + t309) * t445;
t15 = t611 + m(6) * (t275 * t570 + t277 * t569 - t492 * t98);
t688 = t15 * qJD(5);
t417 = rSges(4,1) * t453 + rSges(4,2) * t455;
t369 = t417 * t444;
t370 = t417 * t445;
t472 = -t500 * t690 / 0.2e1;
t486 = m(6) * t492;
t164 = -t486 / 0.2e1 + t472;
t687 = t499 * t164;
t102 = t246 * t256 - t255 * t247;
t120 = t264 * t269 - t268 * t265;
t547 = t444 * t453;
t502 = -rSges(4,2) * t547 - t445 * rSges(4,3);
t291 = -t436 + (pkin(2) + t610) * t444 + t502;
t281 = t291 + t614;
t540 = t445 * t455;
t541 = t445 * t453;
t469 = rSges(4,1) * t540 - rSges(4,2) * t541 + rSges(4,3) * t444;
t292 = t469 + t501;
t282 = t292 + t448;
t129 = t281 * t292 - t291 * t282;
t686 = t701 * t444;
t685 = t701 * t445;
t498 = ((-t282 + t292) * t445 + (-t281 + t291) * t444) * t417 * t666 + (t525 * t275 + t707 * t277) * t664 + ((-t265 + t269) * t694 + t706 * t695) * t665;
t117 = t246 * t279 - t277 * t247;
t121 = t255 * t279 - t277 * t256;
t124 = -t264 * t695 - t265 * t694;
t127 = -t268 * t695 - t269 * t694;
t150 = -t281 * t369 - t282 * t370;
t166 = -t291 * t369 - t292 * t370;
t681 = (t166 + t150) * t666 + (t121 + t117) * t664 + (t127 + t124) * t665;
t679 = t704 * t433 + (t375 - t378) * t434;
t408 = Icges(5,4) * t554;
t318 = Icges(5,1) * t550 + Icges(5,5) * t444 - t408;
t515 = Icges(5,2) * t550 - t318 + t408;
t317 = -Icges(5,5) * t445 + t394 * t444;
t516 = -t391 * t444 + t317;
t316 = Icges(5,4) * t550 - Icges(5,2) * t554 + Icges(5,6) * t444;
t517 = t445 * t683 + t316;
t315 = -Icges(5,6) * t445 + t392 * t444;
t518 = -t444 * t683 - t315;
t678 = -t441 * (-t444 * t515 - t445 * t516) - t442 * (t444 * t517 + t445 * t518);
t423 = Icges(4,4) * t541;
t339 = Icges(4,1) * t540 + Icges(4,5) * t444 - t423;
t511 = -Icges(4,2) * t540 + t339 - t423;
t337 = Icges(4,4) * t540 - Icges(4,2) * t541 + Icges(4,6) * t444;
t513 = t445 * t682 + t337;
t675 = -t453 * t511 - t455 * t513;
t338 = -Icges(4,5) * t445 + t416 * t444;
t512 = -t413 * t444 + t338;
t336 = -Icges(4,6) * t445 + t414 * t444;
t514 = t444 * t682 + t336;
t674 = -t453 * t512 - t455 * t514;
t491 = t704 * t434 / 0.2e1 + (-t375 / 0.2e1 + t378 / 0.2e1) * t433;
t144 = -t300 * t445 - t302 * t563 + t272;
t273 = t305 * t559;
t145 = t301 * t445 + t303 * t563 - t273;
t580 = t304 * t434;
t582 = t302 * t433;
t494 = ((t273 + (t300 - t581) * t444 - t705) * t444 + ((t300 + t477) * t445 + (t580 + t582) * t444 + t708) * t445) * t647 + (-t144 * t445 - t145 * t444) * t646 + ((-t145 + (t301 - t580) * t445 + t705) * t445 + (t144 + (t301 + t582) * t444 + t708) * t444) * t645;
t673 = t703 * t442 / 0.2e1 + t702 * t455 / 0.2e1 + (t416 / 0.2e1 - t413 / 0.2e1) * t453 + (t394 / 0.2e1 - t391 / 0.2e1) * t441;
t671 = 0.4e1 * qJD(1);
t669 = 4 * qJD(2);
t668 = 2 * qJD(3);
t123 = -t246 * t332 - t247 * t333;
t125 = -t255 * t332 - t256 * t333;
t655 = m(6) * (t125 + t123);
t654 = (t525 * t444 + t707 * t445) * t690;
t210 = t246 * t569;
t526 = t275 * t333 - t277 * t332;
t653 = m(6) * (-t247 * t570 + t210 + t526);
t215 = t255 * t569;
t652 = m(6) * (-t256 * t570 + t215 + t526);
t242 = t279 * t333;
t584 = t277 * t380;
t588 = t247 * t381;
t651 = m(6) * (t210 + t242 + (t584 - t588) * t444);
t586 = t256 * t381;
t650 = m(6) * (t215 + t242 + (t584 - t586) * t444);
t551 = t442 * t444;
t283 = t317 * t551;
t390 = Icges(5,5) * t442 - Icges(5,6) * t441;
t568 = t390 * t444;
t313 = -Icges(5,3) * t445 + t568;
t160 = -t313 * t445 - t315 * t555 + t283;
t284 = t318 * t551;
t314 = Icges(5,5) * t550 - Icges(5,6) * t554 + Icges(5,3) * t444;
t161 = t314 * t445 + t316 * t555 - t284;
t103 = -t160 * t445 - t161 * t444;
t285 = t315 * t554;
t162 = -t313 * t444 - t317 * t550 + t285;
t578 = t316 * t441;
t475 = t318 * t442 - t578;
t163 = t314 * t444 + t475 * t445;
t104 = -t162 * t445 - t163 * t444;
t546 = t444 * t455;
t296 = t338 * t546;
t412 = Icges(4,5) * t455 - Icges(4,6) * t453;
t566 = t412 * t444;
t334 = -Icges(4,3) * t445 + t566;
t187 = -t334 * t445 - t336 * t547 + t296;
t297 = t339 * t546;
t335 = Icges(4,5) * t540 - Icges(4,6) * t541 + Icges(4,3) * t444;
t188 = t335 * t445 + t337 * t547 - t297;
t118 = -t187 * t445 - t188 * t444;
t298 = t336 * t541;
t189 = -t334 * t444 - t338 * t540 + t298;
t574 = t337 * t453;
t473 = t339 * t455 - t574;
t190 = t335 * t444 + t473 * t445;
t119 = -t189 * t445 - t190 * t444;
t577 = t317 * t442;
t579 = t315 * t441;
t36 = (t162 + t284 - t285 + (t313 - t578) * t444) * t444 + (-t283 - t163 + (t313 + t475) * t445 + (t577 + t579) * t444) * t445;
t37 = (-t161 + t285 + (t314 - t577) * t445) * t445 + (t160 - t283 + (t314 + t579) * t444) * t444;
t573 = t338 * t455;
t575 = t336 * t453;
t45 = (t189 + t297 - t298 + (t334 - t574) * t444) * t444 + (-t296 - t190 + (t334 + t473) * t445 + (t573 + t575) * t444) * t445;
t46 = (-t188 + t298 + (t335 - t573) * t445) * t445 + (t187 - t296 + (t335 + t575) * t444) * t444;
t2 = (-t46 / 0.2e1 - t119 / 0.2e1 - t104 / 0.2e1 - t37 / 0.2e1) * t445 + (t103 / 0.2e1 + t118 / 0.2e1 - t45 / 0.2e1 - t36 / 0.2e1) * t444 + t494;
t643 = t2 * qJD(3) - qJD(4) * t51;
t638 = m(4) * t129;
t636 = m(4) * t150;
t635 = m(4) * t166;
t632 = m(5) * t120;
t630 = m(5) * t124;
t629 = m(5) * t127;
t628 = m(5) * (-t264 * t445 + t257);
t627 = m(5) * (-t268 * t445 + t263);
t624 = m(6) * t102;
t622 = m(6) * t117;
t621 = m(6) * t121;
t620 = m(6) * t123;
t619 = m(6) * t125;
t618 = m(6) * (-t246 * t445 + t238);
t617 = m(6) * (-t255 * t445 + t245);
t601 = t164 * qJD(4) + qJD(5) * t494;
t165 = t486 / 0.2e1 + t472;
t53 = t530 + t531;
t600 = t53 * qJD(3) + t165 * qJD(5);
t599 = t51 * qJD(3) - t164 * qJD(5);
t576 = t333 * t380;
t529 = (t275 + t279) * t277;
t495 = t381 + t613;
t493 = t332 * t333;
t485 = t655 / 0.2e1 + t491;
t467 = -t494 + (t433 * t521 + t434 * t519 + t445 * t679 - t572) * t647 + (-t374 * t445 + t433 * t522 + t434 * t520 - t444 * t679) * t645;
t466 = -t491 + t693 * (t303 * t434 + t305 * t433);
t460 = t491 + t673;
t459 = t460 + t681;
t458 = t466 - t673 + (t316 * t442 + t318 * t441 + t337 * t455 + t339 * t453) * t693;
t457 = t53 * qJD(4) + (t467 + (t441 * t518 + t442 * t516 - t453 * t514 + t455 * t512 + (-t390 - t412) * t445 - t692 * t444) * t645 + (t104 + t119 + t37 + t46) * t644 + (t441 * t517 + t442 * t515 + t692 * t445 + t453 * t513 - t455 * t511 + t103 + t118 - t566 - t568) * t647 + (t36 + t45) * t646) * qJD(3);
t426 = pkin(3) * t540;
t396 = -rSges(5,2) * t441 + t609;
t323 = t396 * t445 + t426;
t321 = (-t396 - t447) * t444;
t278 = t445 * t495 + t426;
t276 = (-t495 - t447) * t444;
t228 = t445 * t309 + t293;
t155 = t165 * qJD(4);
t134 = t444 * t379 - t312 + (-t500 + t439) * t612 + (-t333 + (t405 + t612) * t445) * t445;
t101 = t617 + t627;
t96 = t618 + t628;
t85 = t491 + t619;
t84 = t491 + t620;
t82 = t650 / 0.2e1;
t79 = t651 / 0.2e1;
t78 = t652 / 0.2e1;
t76 = t653 / 0.2e1;
t68 = t654 / 0.2e1;
t38 = t624 + t632 + t638 + t642;
t33 = t460 + t621 + t629 + t635;
t32 = t460 + t622 + t630 + t636;
t22 = -t654 / 0.2e1 + t485;
t21 = t68 + t485;
t20 = m(6) * (t500 * t380 * t381 - t228 * t492) + t611;
t19 = t20 * qJD(5);
t18 = t602 + t603;
t14 = t68 - t655 / 0.2e1 + t466;
t11 = t459 + t498;
t10 = t459 - t498;
t9 = t78 - t650 / 0.2e1 + t494;
t8 = t82 - t652 / 0.2e1 + t494;
t7 = t76 - t651 / 0.2e1 + t494;
t6 = t79 - t653 / 0.2e1 + t494;
t5 = t458 + t498 - t681;
t4 = t78 + t82 + t467;
t3 = t76 + t79 + t467;
t1 = [qJD(2) * t38 + qJD(3) * t32 + qJD(4) * t96 + qJD(5) * t84, t38 * qJD(1) + t11 * qJD(3) + t18 * qJD(4) + t21 * qJD(5) + 0.2e1 * (t642 / 0.2e1 + t102 * t664 + t120 * t665 + t129 * t666) * qJD(2), t32 * qJD(1) + t11 * qJD(2) + t3 * qJD(5) + ((t281 * t445 - t282 * t444) * t696 + (t264 * t323 + t265 * t321) * t665 + (t246 * t278 + t247 * t276 + t529) * t664) * t668 + t457, qJD(1) * t96 + qJD(2) * t18 + t600, t84 * qJD(1) + t21 * qJD(2) + t3 * qJD(3) + t155 + ((t210 + (t576 - t588) * t444 - t493) * m(6) + t467) * qJD(5); t10 * qJD(3) - t16 * qJD(4) + t22 * qJD(5) + (-t642 / 0.4e1 - t638 / 0.4e1 - t632 / 0.4e1 - t624 / 0.4e1) * t671, qJD(3) * t33 + qJD(4) * t101 + qJD(5) * t85, t10 * qJD(1) + t33 * qJD(2) + t4 * qJD(5) + ((t291 * t445 - t292 * t444) * t696 + (t268 * t323 + t269 * t321) * t665 + (t255 * t278 + t256 * t276 + t529) * t664) * t668 + t457, qJD(2) * t101 + t600 - t711, t22 * qJD(1) + t85 * qJD(2) + t4 * qJD(3) + t155 + ((t215 + (t576 - t586) * t444 - t493) * m(6) + t467) * qJD(5); t458 * qJD(1) + t5 * qJD(2) + t7 * qJD(5) + (-t636 / 0.4e1 - t630 / 0.4e1 - t622 / 0.4e1) * t671 + t643, t5 * qJD(1) + t458 * qJD(2) + t9 * qJD(5) + (-t635 / 0.4e1 - t629 / 0.4e1 - t621 / 0.4e1) * t669 + t643, (m(6) * (t134 * t98 - t275 * t276 + t277 * t278) + m(5) * (-t695 * t321 + t694 * t323 - (t288 + t444 * (rSges(5,1) * t551 + t505) + (rSges(5,3) * t444 - t306 + t489) * t445) * t496 * t500) + m(4) * ((t445 * t469 + t444 * (rSges(4,1) * t546 + t502)) * (-t369 * t444 - t370 * t445) + t500 * t419 * t417) + t611 + ((t674 * t445 + (-t675 - t686) * t444 - t678) * t445 + t685 * t439) * t647 + ((t675 * t444 + (-t674 + t685) * t445 + t678) * t444 - t686 * t440) * t645) * qJD(3) + t688 + t499 * t2, -t712, t7 * qJD(1) + t9 * qJD(2) + t15 * qJD(3) + t688; (-t628 / 0.4e1 - t618 / 0.4e1) * t671 + t16 * qJD(2) + t599, t711 + (-t617 / 0.4e1 - t627 / 0.4e1) * t669 + t599, ((-t276 * t445 - t278 * t444) * t664 + (-t321 * t445 - t323 * t444) * t665) * t668 + t712, 0, -t687; (t466 - t620) * qJD(1) + t14 * qJD(2) + t6 * qJD(3) + t601, t14 * qJD(1) + (t466 - t619) * qJD(2) + t8 * qJD(3) + t601, t6 * qJD(1) + t8 * qJD(2) + ((t134 * t228 + (-t276 * t444 + t278 * t445) * t380) * m(6) + t611) * qJD(3) + t19, t687, qJD(3) * t20 + t494 * t499 + t19;];
Cq = t1;
