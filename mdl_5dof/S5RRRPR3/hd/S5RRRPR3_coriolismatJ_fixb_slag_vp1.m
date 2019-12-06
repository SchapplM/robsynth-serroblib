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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:42:19
% EndTime: 2019-12-05 18:42:35
% DurationCPUTime: 9.38s
% Computational Cost: add. (49536->542), mult. (32605->655), div. (0->0), fcn. (29520->10), ass. (0->337)
t471 = qJ(1) + qJ(2);
t466 = cos(t471);
t470 = qJ(3) + pkin(9);
t464 = qJ(5) + t470;
t456 = sin(t464);
t457 = cos(t464);
t394 = rSges(6,1) * t456 + rSges(6,2) * t457;
t473 = sin(qJ(3));
t629 = t473 * pkin(3);
t462 = sin(t470);
t631 = pkin(4) * t462;
t422 = -t629 - t631;
t702 = (t394 - t422) * t466;
t276 = t466 * t702;
t465 = sin(t471);
t569 = t465 * t473;
t447 = pkin(3) * t569;
t282 = t447 + (t394 + t631) * t465;
t463 = cos(t470);
t409 = rSges(5,1) * t462 + rSges(5,2) * t463;
t335 = t409 * t465 + t447;
t684 = m(6) / 0.2e1;
t685 = m(5) / 0.2e1;
t502 = (t409 + t629) * t466;
t713 = t466 * t502;
t551 = (-t282 * t465 - t276) * t684 + (-t335 * t465 - t713) * t685;
t580 = t457 * t465;
t584 = t456 * t465;
t347 = rSges(6,1) * t584 + rSges(6,2) * t580;
t286 = -t422 * t465 + t347;
t572 = t463 * t465;
t576 = t462 * t465;
t520 = rSges(5,1) * t576 + rSges(5,2) * t572;
t318 = t447 + t520;
t552 = (t286 * t465 + t276) * t684 + (t318 * t465 + t713) * t685;
t51 = t552 - t551;
t514 = qJD(1) + qJD(2);
t732 = t514 * t51;
t475 = cos(qJ(3));
t468 = t475 * pkin(3);
t459 = t468 + pkin(2);
t630 = pkin(4) * t463;
t421 = t459 + t630;
t472 = -qJ(4) - pkin(7);
t469 = -pkin(8) + t472;
t522 = -rSges(6,2) * t584 - rSges(6,3) * t466;
t625 = rSges(6,1) * t457;
t258 = t466 * t469 + (t421 + t625) * t465 + t522;
t633 = sin(qJ(1)) * pkin(1);
t249 = t258 + t633;
t240 = t466 * t249;
t388 = t466 * t421;
t579 = t457 * t466;
t583 = t456 * t466;
t504 = -rSges(6,1) * t579 + rSges(6,2) * t583;
t259 = -t388 + (-rSges(6,3) + t469) * t465 + t504;
t632 = cos(qJ(1)) * pkin(1);
t250 = t259 - t632;
t521 = -rSges(5,2) * t576 - rSges(5,3) * t466;
t562 = t466 * t472;
t626 = rSges(5,1) * t463;
t274 = t562 + (t459 + t626) * t465 + t521;
t269 = t274 + t633;
t260 = t466 * t269;
t571 = t463 * t466;
t575 = t462 * t466;
t485 = rSges(5,1) * t571 - rSges(5,2) * t575 + rSges(5,3) * t465;
t516 = -t459 * t466 + t465 * t472;
t275 = -t485 + t516;
t270 = t275 - t632;
t563 = t466 * t274;
t565 = t466 * t258;
t619 = (-t260 - t563 + (-t270 - t275) * t465) * t685 + (-t240 - t565 + (-t250 - t259) * t465) * t684;
t546 = t250 - t259;
t620 = (-t260 + t563 + (-t270 + t275) * t465) * t685 + (-t465 * t546 - t240 + t565) * t684;
t16 = t620 - t619;
t731 = t16 * qJD(1);
t460 = t465 ^ 2;
t461 = t466 ^ 2;
t515 = t460 + t461;
t730 = t394 * t515;
t450 = Icges(6,4) * t457;
t705 = Icges(6,2) * t456 - t450;
t312 = Icges(6,6) * t466 + t465 * t705;
t729 = t312 * t583;
t311 = Icges(6,5) * t579 - Icges(6,6) * t583 + Icges(6,3) * t465;
t389 = Icges(6,5) * t457 - Icges(6,6) * t456;
t594 = t389 * t465;
t544 = t466 * (Icges(6,3) * t466 - t594) + t312 * t584;
t313 = Icges(6,4) * t579 - Icges(6,2) * t583 + Icges(6,6) * t465;
t416 = Icges(6,4) * t583;
t315 = Icges(6,1) * t579 + Icges(6,5) * t465 - t416;
t714 = t466 * (t313 * t456 - t315 * t457);
t728 = t311 * t465 + t544 - t714;
t727 = -t249 + t258;
t726 = Icges(4,5) * t473 + Icges(5,5) * t462 + Icges(4,6) * t475 + Icges(5,6) * t463;
t686 = m(4) / 0.2e1;
t668 = -t465 / 0.2e1;
t667 = t465 / 0.2e1;
t665 = t466 / 0.2e1;
t663 = m(3) * (-t632 * (rSges(3,1) * t465 + rSges(3,2) * t466) - (-rSges(3,1) * t466 + rSges(3,2) * t465) * t633);
t725 = t250 * t702;
t724 = t259 * t702;
t666 = -t466 / 0.2e1;
t723 = t665 + t666;
t611 = Icges(5,4) * t462;
t405 = Icges(5,2) * t463 + t611;
t408 = Icges(5,1) * t463 - t611;
t612 = Icges(4,4) * t473;
t433 = Icges(4,2) * t475 + t612;
t436 = Icges(4,1) * t475 - t612;
t467 = Icges(4,4) * t475;
t703 = Icges(4,2) * t473 - t467;
t451 = Icges(5,4) * t463;
t704 = Icges(5,2) * t462 - t451;
t706 = Icges(4,1) * t473 + t467;
t707 = Icges(5,1) * t462 + t451;
t722 = (-t703 + t706) * t473 + (t433 - t436) * t475 + (-t704 + t707) * t462 + (t405 - t408) * t463;
t535 = -Icges(6,2) * t579 + t315 - t416;
t415 = Icges(6,4) * t584;
t314 = -Icges(6,1) * t580 + Icges(6,5) * t466 + t415;
t536 = Icges(6,2) * t580 + t314 + t415;
t708 = Icges(6,1) * t456 + t450;
t537 = t466 * t708 + t313;
t538 = -t465 * t708 + t312;
t721 = -(t465 * t535 + t466 * t536) * t456 - (t465 * t537 + t466 * t538) * t457;
t348 = t394 * t466;
t125 = -t258 * t347 + t259 * t348;
t720 = m(6) * t125;
t395 = -rSges(6,2) * t456 + t625;
t591 = t395 * t465;
t495 = Icges(6,5) * t456 + Icges(6,6) * t457;
t341 = t465 * t495;
t342 = t495 * t466;
t628 = (t461 * t341 + (-t466 * t342 - t721) * t465) * t665 + (-t460 * t342 + (t465 * t341 + t721) * t466) * t667;
t324 = t466 * t348;
t711 = t465 * t347 + t324;
t298 = t466 * (t466 * pkin(2) + pkin(7) * t465 + t516);
t301 = t466 * (rSges(6,3) * t465 - t504);
t458 = t466 * pkin(7);
t316 = -t562 - t458 + (pkin(2) - t459) * t465;
t319 = -rSges(6,1) * t580 - t522;
t98 = -t466 * (-t388 - t516) - t298 + t301 + (-t316 - t319 + (t421 - t459) * t465 - t562) * t465;
t15 = t628 + m(6) * (t276 * t395 + t282 * t591 - t711 * t98);
t718 = t15 * qJD(5);
t627 = rSges(4,1) * t475;
t439 = -rSges(4,2) * t473 + t627;
t717 = t439 * t686;
t560 = t466 * t475;
t561 = t466 * t473;
t352 = Icges(4,4) * t560 - Icges(4,2) * t561 + Icges(4,6) * t465;
t444 = Icges(4,4) * t561;
t354 = Icges(4,1) * t560 + Icges(4,5) * t465 - t444;
t716 = (t352 * t473 - t354 * t475) * t466;
t331 = Icges(5,4) * t571 - Icges(5,2) * t575 + Icges(5,6) * t465;
t426 = Icges(5,4) * t575;
t333 = Icges(5,1) * t571 + Icges(5,5) * t465 - t426;
t715 = (t331 * t462 - t333 * t463) * t466;
t213 = t502 * t270;
t225 = t502 * t275;
t486 = m(6) * t711;
t488 = -m(6) * t730 / 0.2e1;
t164 = -t486 / 0.2e1 + t488;
t712 = t514 * t164;
t102 = -t249 * t259 + t250 * t258;
t120 = -t269 * t275 + t270 * t274;
t511 = pkin(2) + t627;
t517 = -rSges(4,2) * t569 - rSges(4,3) * t466;
t299 = t465 * t511 - t458 + t517;
t288 = t299 + t633;
t446 = rSges(4,2) * t561;
t300 = t446 - t511 * t466 + (-rSges(4,3) - pkin(7)) * t465;
t289 = t300 - t632;
t129 = -t288 * t300 + t289 * t299;
t710 = t726 * t465;
t709 = t726 * t466;
t437 = rSges(4,1) * t473 + rSges(4,2) * t475;
t385 = t437 * t466;
t384 = t437 * t465;
t600 = t299 * t384;
t601 = t288 * t384;
t513 = (-t601 + t600 + (t289 - t300) * t385) * t686 + (t282 * t727 - t724 + t725) * t684 + (t213 - t225 + (-t269 + t274) * t335) * t685;
t117 = -t249 * t286 + t725;
t121 = -t258 * t286 + t724;
t124 = -t269 * t318 + t213;
t127 = -t274 * t318 + t225;
t150 = t289 * t385 - t601;
t166 = t300 * t385 - t600;
t701 = (t166 + t150) * t686 + (t121 + t117) * t684 + (t127 + t124) * t685;
t610 = Icges(6,4) * t456;
t390 = Icges(6,2) * t457 + t610;
t393 = Icges(6,1) * t457 - t610;
t700 = (-t705 + t708) * t456 + (t390 - t393) * t457;
t506 = (-t705 / 0.2e1 + t708 / 0.2e1) * t457 + (-t390 / 0.2e1 + t393 / 0.2e1) * t456;
t144 = -t314 * t580 + t544;
t145 = t466 * t311 + t313 * t584 - t315 * t580;
t509 = -t314 * t457 - t311;
t598 = t312 * t456;
t13 = ((t714 + t728) * t466 + ((t509 - t598) * t466 + t145 + t729) * t465) * t667 + (t144 * t466 + t145 * t465) * t668 + ((t145 + (-t311 + t598) * t466 - t729) * t466 + (t465 * t509 - t144 + t728) * t465) * t665;
t527 = -Icges(4,2) * t560 + t354 - t444;
t529 = t466 * t706 + t352;
t531 = -Icges(5,2) * t571 + t333 - t426;
t533 = t466 * t707 + t331;
t695 = t462 * t531 + t463 * t533 + t473 * t527 + t475 * t529;
t443 = Icges(4,4) * t569;
t568 = t465 * t475;
t353 = -Icges(4,1) * t568 + Icges(4,5) * t466 + t443;
t528 = Icges(4,2) * t568 + t353 + t443;
t351 = Icges(4,6) * t466 + t465 * t703;
t530 = -t465 * t706 + t351;
t425 = Icges(5,4) * t576;
t332 = -Icges(5,1) * t572 + Icges(5,5) * t466 + t425;
t532 = Icges(5,2) * t572 + t332 + t425;
t330 = Icges(5,6) * t466 + t465 * t704;
t534 = -t465 * t707 + t330;
t694 = -t462 * t532 - t463 * t534 - t473 * t528 - t475 * t530;
t693 = (t706 / 0.2e1 - t703 / 0.2e1) * t475 + (t436 / 0.2e1 - t433 / 0.2e1) * t473 + (t707 / 0.2e1 - t704 / 0.2e1) * t463 + (t408 / 0.2e1 - t405 / 0.2e1) * t462;
t691 = 0.4e1 * qJD(1);
t689 = 4 * qJD(2);
t688 = 2 * qJD(3);
t123 = -t249 * t347 + t250 * t348;
t676 = m(6) * (t125 + t123);
t592 = t394 * t465;
t675 = m(6) * (t546 * t348 + t592 * t727);
t211 = t250 * t591;
t547 = t282 * t348 - t347 * t702;
t674 = m(6) * (t240 * t395 + t211 + t547);
t217 = t259 * t591;
t673 = m(6) * (t395 * t565 + t217 + t547);
t244 = t702 * t592;
t602 = t286 * t394;
t606 = t249 * t395;
t672 = m(6) * (t211 + t244 + (-t602 + t606) * t466);
t605 = t258 * t395;
t671 = m(6) * (t217 + t244 + (-t602 + t605) * t466);
t404 = Icges(5,5) * t463 - Icges(5,6) * t462;
t590 = t404 * t465;
t328 = Icges(5,3) * t466 - t590;
t542 = t328 * t466 + t330 * t576;
t160 = -t332 * t572 + t542;
t329 = Icges(5,5) * t571 - Icges(5,6) * t575 + Icges(5,3) * t465;
t161 = t466 * t329 + t331 * t576 - t333 * t572;
t103 = t160 * t466 + t161 * t465;
t541 = t328 * t465 + t332 * t571;
t162 = -t330 * t575 + t541;
t163 = t329 * t465 - t715;
t104 = t162 * t466 + t163 * t465;
t432 = Icges(4,5) * t475 - Icges(4,6) * t473;
t587 = t432 * t465;
t349 = Icges(4,3) * t466 - t587;
t540 = t349 * t466 + t351 * t569;
t187 = -t353 * t568 + t540;
t350 = Icges(4,5) * t560 - Icges(4,6) * t561 + Icges(4,3) * t465;
t188 = t466 * t350 + t352 * t569 - t354 * t568;
t118 = t187 * t466 + t188 * t465;
t539 = t349 * t465 + t353 * t560;
t189 = -t351 * t561 + t539;
t190 = t350 * t465 - t716;
t119 = t189 * t466 + t190 * t465;
t508 = -t332 * t463 - t329;
t597 = t330 * t462;
t36 = (t163 + t542 + t715) * t466 + (-t162 + (t508 - t597) * t466 + t161 + t541) * t465;
t37 = (t161 + (-t329 + t597) * t466 - t541) * t466 + (t465 * t508 - t160 + t542) * t465;
t507 = -t353 * t475 - t350;
t595 = t351 * t473;
t45 = (t190 + t540 + t716) * t466 + (-t189 + (t507 - t595) * t466 + t188 + t539) * t465;
t46 = (t188 + (-t350 + t595) * t466 - t539) * t466 + (t465 * t507 - t187 + t540) * t465;
t2 = (t104 / 0.2e1 + t119 / 0.2e1 + t37 / 0.2e1 + t46 / 0.2e1) * t466 + (t36 / 0.2e1 + t45 / 0.2e1 - t103 / 0.2e1 - t118 / 0.2e1) * t465 + t13;
t664 = qJD(3) * t2 - qJD(4) * t51;
t659 = m(4) * t129;
t657 = m(4) * t150;
t656 = m(4) * t166;
t653 = m(5) * t120;
t651 = m(5) * t124;
t650 = m(5) * t127;
t649 = m(5) * (-t270 * t465 - t260);
t648 = m(5) * (-t275 * t465 - t563);
t644 = m(6) * t102;
t642 = m(6) * t117;
t641 = m(6) * t121;
t640 = m(6) * t123;
t638 = m(6) * (-t250 * t465 - t240);
t637 = m(6) * (-t259 * t465 - t565);
t618 = qJD(4) * t164 + qJD(5) * t13;
t165 = t486 / 0.2e1 + t488;
t53 = t551 + t552;
t617 = qJD(3) * t53 + qJD(5) * t165;
t616 = qJD(3) * t51 - qJD(5) * t164;
t596 = t347 * t394;
t550 = (t282 - t286) * t702;
t545 = (-t318 + t335) * t502;
t510 = t395 + t630;
t501 = t676 / 0.2e1 + t506;
t484 = -t13 + (-t456 * t537 + t457 * t535 - t466 * t700 + t594) * t667 + (t389 * t466 - t456 * t538 + t457 * t536 + t465 * t700) * t665;
t483 = -t506 + t723 * (t313 * t457 + t315 * t456);
t480 = t506 + t693;
t479 = t480 + t701;
t478 = t483 - t693 + (t331 * t463 + t333 * t462 + t352 * t475 + t354 * t473) * t723;
t477 = t53 * qJD(4) + (t484 + (t104 + t119 + t37 + t46) * t666 + (-t462 * t534 + t463 * t532 - t473 * t530 + t475 * t528 + (t404 + t432) * t466 + t722 * t465) * t665 + (-t462 * t533 + t463 * t531 - t466 * t722 - t473 * t529 + t475 * t527 + t103 + t118 + t587 + t590) * t667 + (t36 + t45) * t668) * qJD(3);
t448 = pkin(3) * t568;
t410 = -rSges(5,2) * t462 + t626;
t338 = (-t410 - t468) * t466;
t336 = t410 * t465 + t448;
t285 = (-t510 - t468) * t466;
t283 = t465 * t510 + t448;
t277 = t348 * t592;
t230 = -t319 * t465 + t301;
t155 = t165 * qJD(4);
t134 = t422 * t461 - t324 + (t447 - t286) * t465 + (-t515 + t461) * t629;
t101 = t637 + t648;
t96 = t638 + t649;
t85 = t506 + t720;
t84 = t506 + t640;
t82 = t671 / 0.2e1;
t79 = t672 / 0.2e1;
t78 = t673 / 0.2e1;
t76 = t674 / 0.2e1;
t68 = t675 / 0.2e1;
t38 = t644 + t653 + t659 + t663;
t33 = t480 + t641 + t650 + t656;
t32 = t480 + t642 + t651 + t657;
t22 = -t675 / 0.2e1 + t501;
t21 = t68 + t501;
t20 = m(6) * (-t230 * t711 + t395 * t730) + t628;
t19 = t20 * qJD(5);
t18 = t619 + t620;
t14 = t68 - t676 / 0.2e1 + t483;
t11 = t479 + t513;
t10 = t479 - t513;
t9 = t78 - t671 / 0.2e1 + t13;
t8 = t82 - t673 / 0.2e1 + t13;
t7 = t76 - t672 / 0.2e1 + t13;
t6 = t79 - t674 / 0.2e1 + t13;
t5 = t478 + t513 - t701;
t4 = t78 + t82 + t484;
t3 = t76 + t79 + t484;
t1 = [qJD(2) * t38 + qJD(3) * t32 + qJD(4) * t96 + qJD(5) * t84, t38 * qJD(1) + t11 * qJD(3) + t18 * qJD(4) + t21 * qJD(5) + 0.2e1 * (t663 / 0.2e1 + t102 * t684 + t120 * t685 + t129 * t686) * qJD(2), t32 * qJD(1) + t11 * qJD(2) + t3 * qJD(5) + ((t288 * t466 + t289 * t465) * t717 + (-t269 * t338 + t270 * t336 + t545) * t685 + (-t249 * t285 + t250 * t283 + t550) * t684) * t688 + t477, qJD(1) * t96 + qJD(2) * t18 + t617, t84 * qJD(1) + t21 * qJD(2) + t3 * qJD(3) + t155 + ((t211 + t277 + (-t596 + t606) * t466) * m(6) + t484) * qJD(5); t10 * qJD(3) - t16 * qJD(4) + t22 * qJD(5) + (-t663 / 0.4e1 - t659 / 0.4e1 - t653 / 0.4e1 - t644 / 0.4e1) * t691, qJD(3) * t33 + qJD(4) * t101 + qJD(5) * t85, t10 * qJD(1) + t33 * qJD(2) + t4 * qJD(5) + ((t299 * t466 + t300 * t465) * t717 + (-t274 * t338 + t275 * t336 + t545) * t685 + (-t258 * t285 + t259 * t283 + t550) * t684) * t688 + t477, qJD(2) * t101 + t617 - t731, t22 * qJD(1) + t85 * qJD(2) + t4 * qJD(3) + t155 + ((t217 + t277 + (-t596 + t605) * t466) * m(6) + t484) * qJD(5); t478 * qJD(1) + t5 * qJD(2) + t7 * qJD(5) + (-t657 / 0.4e1 - t651 / 0.4e1 - t642 / 0.4e1) * t691 + t664, t5 * qJD(1) + t478 * qJD(2) + t9 * qJD(5) + (-t656 / 0.4e1 - t650 / 0.4e1 - t641 / 0.4e1) * t689 + t664, (m(4) * ((t466 * (rSges(4,1) * t560 + rSges(4,3) * t465 - t446) - t465 * (-rSges(4,1) * t568 - t517)) * (-t384 * t465 - t385 * t466) + t515 * t439 * t437) + m(6) * (t134 * t98 + t282 * t283 - t285 * t702) + m(5) * ((t466 * t485 - t298 + (rSges(5,1) * t572 - t316 + t521) * t465) * (-t409 * t461 - t465 * t520 - t515 * t629) + t335 * t336 - t502 * t338) + t628 + ((t694 * t466 + (-t695 + t710) * t465) * t466 - t709 * t460) * t667 + ((t695 * t465 + (-t694 - t709) * t466) * t465 + t710 * t461) * t665) * qJD(3) + t718 + t514 * t2, -t732, qJD(1) * t7 + qJD(2) * t9 + qJD(3) * t15 + t718; t16 * qJD(2) + (-t649 / 0.4e1 - t638 / 0.4e1) * t691 + t616, t731 + (-t637 / 0.4e1 - t648 / 0.4e1) * t689 + t616, ((t283 * t466 + t285 * t465) * t684 + (t336 * t466 + t338 * t465) * t685) * t688 + t732, 0, -t712; (t483 - t640) * qJD(1) + t14 * qJD(2) + t6 * qJD(3) + t618, t14 * qJD(1) + (t483 - t720) * qJD(2) + t8 * qJD(3) + t618, t6 * qJD(1) + t8 * qJD(2) + ((t134 * t230 + (t283 * t465 - t285 * t466) * t394) * m(6) + t628) * qJD(3) + t19, t712, qJD(3) * t20 + t13 * t514 + t19;];
Cq = t1;
