% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:36
% DurationCPUTime: 24.69s
% Computational Cost: add. (4794->593), mult. (12419->701), div. (0->0), fcn. (9450->4), ass. (0->323)
t719 = Icges(6,4) + Icges(5,5);
t718 = Icges(5,1) + Icges(6,1);
t712 = -Icges(6,2) - Icges(5,3);
t717 = Icges(5,6) - Icges(6,6);
t336 = sin(qJ(1));
t338 = cos(qJ(1));
t337 = cos(qJ(3));
t335 = sin(qJ(3));
t555 = Icges(4,4) * t335;
t409 = Icges(4,2) * t337 + t555;
t716 = -t338 * t409 + (Icges(4,6) - Icges(5,6)) * t336;
t715 = t719 * t337;
t554 = Icges(4,4) * t337;
t412 = Icges(4,1) * t335 + t554;
t714 = -t338 * t412 + (Icges(5,4) + Icges(4,5)) * t336;
t713 = Icges(5,4) - Icges(6,5);
t527 = t337 * t338;
t530 = t335 * t338;
t709 = Icges(5,5) * t530 - Icges(5,3) * t527 + t716;
t291 = Icges(5,5) * t527;
t711 = -Icges(5,1) * t530 + t291 + t714;
t248 = Icges(4,1) * t337 - t555;
t710 = t719 * t335;
t687 = t337 * t718 + t710;
t676 = t248 + t687;
t684 = t335 * t718 - t715;
t293 = Icges(6,4) * t527;
t549 = Icges(6,5) * t336;
t150 = Icges(6,1) * t530 - t293 + t549;
t680 = -t150 + t711;
t531 = t335 * t336;
t708 = t719 * t531;
t685 = t337 * t712 + t710;
t707 = -t409 + t685;
t706 = -t412 - t684;
t705 = t717 * t338;
t683 = t335 * t712 - t715;
t528 = t336 * t337;
t698 = t528 * t712 + t705 + t708;
t405 = Icges(4,5) * t335 + Icges(4,6) * t337;
t141 = Icges(4,3) * t338 + t336 * t405;
t408 = Icges(5,4) * t335 - Icges(5,6) * t337;
t145 = Icges(5,2) * t338 + t336 * t408;
t704 = -t141 - t145;
t147 = Icges(4,6) * t338 + t336 * t409;
t703 = t147 - t698;
t544 = Icges(6,6) * t336;
t144 = Icges(6,4) * t530 - Icges(6,2) * t527 + t544;
t702 = -t144 - t709;
t690 = t336 * t684 + t338 * t713;
t701 = t676 * t338;
t700 = t335 * t711 + t337 * t709;
t699 = (-Icges(4,5) - t713) * t337 + (Icges(4,6) - t717) * t335;
t142 = -Icges(4,3) * t336 + t338 * t405;
t548 = Icges(5,2) * t336;
t146 = Icges(5,4) * t530 - Icges(5,6) * t527 - t548;
t697 = -t142 - t146;
t402 = Icges(6,5) * t335 - Icges(6,6) * t337;
t669 = -t402 + t408;
t695 = -t405 - t669;
t694 = t697 * t338 + t528 * t709 + t531 * t680;
t242 = -Icges(4,2) * t335 + t554;
t388 = t242 * t337 + t248 * t335;
t675 = t335 * t687 + t337 * t683 + t388;
t541 = Icges(6,3) * t336;
t138 = Icges(6,5) * t530 - Icges(6,6) * t527 + t541;
t631 = t138 * t338 + t144 * t528 + t694;
t294 = Icges(4,4) * t528;
t551 = Icges(4,5) * t338;
t153 = Icges(4,1) * t531 + t294 + t551;
t54 = t141 * t338 + t147 * t528 + t153 * t531;
t137 = -Icges(6,3) * t338 + t336 * t402;
t644 = t338 * t145 + t531 * t690;
t647 = t137 * t338 + t528 * t698 - t644;
t693 = t54 - t647;
t689 = t683 * t336;
t192 = t248 * t336;
t686 = (t528 * t718 + t192 - t703 + t708) * t338 + (-t701 + t702) * t336;
t682 = t700 * t338;
t626 = t336 * t704 - t527 * t698;
t681 = -t153 - t690;
t625 = t699 * t336;
t679 = t707 * qJD(3);
t678 = t706 * qJD(3);
t677 = -t242 - t683;
t620 = t699 * t338;
t674 = t698 * t337;
t395 = t147 * t337 + t153 * t335;
t673 = -t137 * t336 - t395 * t338 - t530 * t690 - t626;
t537 = t144 * t337;
t397 = -t150 * t335 + t537;
t629 = -t397 * t338 - t682 + (t138 + t697) * t336;
t672 = t336 * t675 - t620;
t671 = t338 * t388 + t527 * t683 + t530 * t687 + t625;
t622 = t335 * t690 + t395;
t645 = rSges(6,3) + qJ(5);
t670 = -t138 + t146;
t317 = t338 * qJ(2);
t253 = pkin(1) * t336 - t317;
t230 = qJD(1) * t253;
t314 = qJD(2) * t336;
t476 = qJD(1) * t338;
t490 = qJ(2) * t476 + t314;
t668 = t490 - t314 + t230;
t187 = t242 * t338;
t661 = (Icges(4,2) * t531 - t294 + t681 - t689) * t338 + (t530 * t712 + t187 - t291 - t293 - t680) * t336;
t667 = -t335 * t686 + t337 * t661;
t472 = qJD(3) * t338;
t666 = qJD(3) * t187 + t683 * t472 + (t336 * t685 - t147 + t705) * qJD(1);
t474 = qJD(3) * t336;
t665 = -t242 * t474 - t689 * qJD(3) + (t338 * t685 + t544 + t716) * qJD(1);
t628 = t335 * t703 + t337 * t681;
t664 = -t701 * qJD(3) + (t336 * t412 + t551 + t690) * qJD(1);
t627 = t335 * t702 + t337 * t680;
t663 = qJD(3) * t192 + t687 * t474 + (t338 * t684 + t549 - t714) * qJD(1);
t662 = t675 * qJD(1) + qJD(3) * t695;
t660 = -t397 - t700;
t659 = -t622 + t674;
t658 = t677 + t706;
t657 = t676 + t707;
t656 = t679 * t337 + t678 * t335 + (t335 * t677 + t337 * t676) * qJD(3) + t699 * qJD(1);
t654 = (t335 * t658 + t337 * t657) * qJD(1);
t653 = t672 * qJD(1);
t652 = (t336 * t629 + t338 * t673) * qJD(3);
t651 = (t336 * t631 + t338 * t693) * qJD(3);
t650 = (-t137 - t704) * qJD(1);
t649 = t671 * qJD(1);
t648 = rSges(6,1) + pkin(4);
t329 = t338 * rSges(4,3);
t159 = rSges(4,1) * t531 + rSges(4,2) * t528 + t329;
t260 = pkin(1) * t338 + qJ(2) * t336;
t331 = t338 * pkin(6);
t602 = t331 + t260;
t643 = t159 + t602;
t301 = pkin(4) * t531;
t642 = -t338 * t645 + t301;
t640 = t146 - t674;
t473 = qJD(3) * t337;
t448 = t336 * t473;
t454 = t335 * t476;
t639 = t448 + t454;
t477 = qJD(1) * t336;
t638 = t335 * t472 + t337 * t477;
t637 = t651 + t653;
t636 = -t649 + t652;
t635 = qJD(3) * t659 + t335 * t665 + t337 * t663;
t634 = qJD(3) * t660 + t335 * t666 + t337 * t664;
t633 = t336 * t662 - t338 * t656;
t632 = t336 * t656 + t338 * t662;
t482 = qJD(1) * t142;
t624 = t482 - t625 * qJD(3) + (t338 * t669 - t541 - t548 - t659) * qJD(1);
t623 = -qJD(1) * t660 + qJD(3) * t620 + t650;
t621 = t695 * qJD(1);
t618 = qJD(3) * t628 - t335 * t663 + t337 * t665 + t650;
t617 = qJD(1) * t670 + qJD(3) * t627 + t335 * t664 - t337 * t666 + t482;
t615 = 0.2e1 * qJD(3);
t257 = rSges(6,1) * t337 + rSges(6,2) * t335;
t201 = t257 * t338;
t258 = rSges(5,1) * t337 + rSges(5,3) * t335;
t313 = qJD(4) * t337;
t609 = (qJD(3) * t258 - t313) * t338;
t295 = rSges(6,1) * t531;
t607 = rSges(6,2) * t528 - t295 - t642;
t285 = qJ(4) * t527;
t199 = pkin(3) * t530 - t285;
t606 = qJD(1) * t199 - t230;
t577 = pkin(6) * t336;
t436 = -t253 - t577;
t433 = -rSges(3,2) * t338 + rSges(3,3) * t336;
t605 = t260 + t433;
t491 = -rSges(5,2) * t336 - rSges(5,3) * t527;
t450 = t335 * t474;
t601 = rSges(6,1) * t639 + rSges(6,2) * t450 + pkin(4) * t454 + t477 * t645;
t333 = t336 ^ 2;
t586 = m(5) / 0.2e1;
t585 = m(6) / 0.2e1;
t584 = -pkin(1) - pkin(6);
t583 = t336 / 0.2e1;
t581 = rSges(3,2) - pkin(1);
t580 = pkin(3) * t335;
t579 = pkin(4) * t335;
t578 = pkin(4) * t337;
t573 = rSges(6,2) * t337;
t572 = rSges(3,3) * t338;
t415 = rSges(6,1) * t335 - t573;
t224 = t415 * qJD(3);
t437 = t257 + t578;
t358 = qJD(3) * t437 - t313;
t340 = qJD(1) ^ 2;
t468 = qJD(1) * qJD(2);
t502 = qJD(1) * (-pkin(1) * t477 + t490) + t336 * t468;
t386 = -t340 * t577 + t502;
t256 = pkin(3) * t337 + qJ(4) * t335;
t467 = qJD(1) * qJD(3);
t445 = t256 * t467;
t459 = pkin(3) * t639 + qJ(4) * t450;
t471 = qJD(4) * t336;
t86 = (-qJ(4) * t476 - t471) * t337 + t459;
t375 = qJD(1) * t86 + t336 * t445 + t386;
t249 = qJ(4) * t337 - t580;
t312 = qJD(4) * t335;
t165 = qJD(3) * t249 + t312;
t431 = -t165 - t312;
t463 = qJD(3) ^ 2 * t579;
t469 = qJD(5) * t338;
t275 = pkin(4) * t448;
t430 = t275 - t469;
t452 = t337 * t476;
t522 = rSges(6,2) * t452 - t430 - t601;
t15 = (t463 + (t224 + t431) * qJD(3)) * t338 + (t336 * t358 - t469 - t522) * qJD(1) + t375;
t570 = t15 * t338;
t311 = qJD(5) * t336;
t284 = t335 * t471;
t305 = t338 * t468;
t426 = -t331 * t340 + t305;
t374 = qJD(3) * t284 + t165 * t474 + t338 * t445 + t426;
t447 = t337 * t472;
t493 = -pkin(4) * t447 - t311;
t523 = -qJD(3) * t201 + t493 + (t336 * t415 + t642) * qJD(1);
t315 = qJD(2) * t338;
t166 = qJD(1) * t260 - t315;
t302 = pkin(3) * t531;
t458 = pkin(3) * t447 + qJ(4) * t638;
t470 = qJD(4) * t338;
t85 = qJD(1) * t302 + t337 * t470 - t458;
t559 = -t166 - t85;
t16 = (-qJD(3) * t224 - t463) * t336 + (t338 * t358 + t311 - t523 + t559) * qJD(1) + t374;
t569 = t16 * t336;
t457 = rSges(5,1) * t639 + rSges(5,3) * t450;
t111 = qJD(1) * t491 + t457;
t416 = rSges(5,1) * t335 - rSges(5,3) * t337;
t225 = t416 * qJD(3);
t446 = t336 * t313;
t23 = (t111 - t446) * qJD(1) + (t258 * t477 + (t225 + t431) * t338) * qJD(3) + t375;
t568 = t23 * t338;
t202 = t258 * t338;
t330 = t338 * rSges(5,2);
t108 = -qJD(3) * t202 + (t336 * t416 + t330) * qJD(1);
t24 = -t225 * t474 + (-t108 + t559 + t609) * qJD(1) + t374;
t567 = t24 * t336;
t456 = rSges(4,1) * t639 + rSges(4,2) * t452;
t475 = qJD(3) * t335;
t112 = (-rSges(4,2) * t475 - rSges(4,3) * qJD(1)) * t336 + t456;
t259 = rSges(4,1) * t337 - rSges(4,2) * t335;
t208 = t259 * t474;
t417 = rSges(4,1) * t335 + rSges(4,2) * t337;
t226 = t417 * qJD(3);
t43 = t226 * t472 + (t112 + t208) * qJD(1) + t386;
t566 = t338 * t43;
t326 = t336 * rSges(4,3);
t162 = t338 * t417 - t326;
t64 = t208 + t314 + (t162 + t436) * qJD(1);
t565 = t338 * t64;
t203 = t259 * t338;
t109 = -qJD(3) * t203 + (t336 * t417 + t329) * qJD(1);
t451 = t259 * t472;
t44 = -t226 * t474 + (-t109 - t166 + t451) * qJD(1) + t426;
t564 = t44 * t336;
t563 = -rSges(6,2) - qJ(4);
t562 = -rSges(5,3) - qJ(4);
t560 = -t111 - t86;
t558 = qJD(3) * t313 + t472 * t85;
t529 = t336 * t256;
t519 = t165 * t336 + t256 * t476;
t200 = t256 * t338;
t506 = -t200 * t472 + t313;
t492 = rSges(5,1) * t531 + t330;
t158 = -rSges(5,3) * t528 + t492;
t194 = -qJ(4) * t528 + t302;
t505 = -t158 - t194;
t300 = rSges(6,2) * t527;
t504 = t336 * t645 + t530 * t648 - t300;
t161 = rSges(5,1) * t530 + t491;
t503 = t161 + t199;
t501 = -t256 * t472 - t315;
t489 = rSges(3,2) * t477 + rSges(3,3) * t476;
t488 = -t199 * t472 + t312;
t486 = t317 - t285;
t485 = t338 ^ 2 + t333;
t466 = pkin(3) + t648;
t465 = -rSges(5,2) + t584;
t464 = -rSges(4,3) + t584;
t462 = -t86 + t522;
t461 = -t194 + t607;
t460 = t199 + t504;
t442 = -t474 / 0.2e1;
t440 = -t472 / 0.2e1;
t439 = t472 / 0.2e1;
t438 = t645 + t584;
t429 = t315 + t458;
t427 = t302 + t602;
t425 = t194 + t602;
t424 = t199 + t436;
t420 = qJD(1) * t529 - t335 * t470;
t419 = qJD(1) * t200 + t249 * t474 + t284;
t65 = qJD(1) * t643 - t315 - t451;
t414 = t336 * t64 - t338 * t65;
t413 = t459 + t490;
t384 = t256 * t474 + t314 - t446;
t376 = t258 * t474 + t384;
t66 = (-t159 * t336 - t162 * t338) * qJD(3);
t348 = t257 * t474 + t384 + t430;
t254 = rSges(3,2) * t336 + t572;
t210 = t256 * t477;
t198 = t259 * t336;
t197 = t258 * t336;
t174 = t485 * t473;
t173 = t450 - t452;
t164 = t338 * t199;
t116 = qJD(1) * t605 - t315;
t115 = t314 + (-t253 + t254) * qJD(1);
t106 = t305 + (-qJD(1) * t433 - t166) * qJD(1);
t105 = qJD(1) * t489 + t502;
t83 = t338 * t85;
t49 = (-t161 * t338 + t336 * t505) * qJD(3) + t488;
t48 = -t609 + (t158 + t425) * qJD(1) + t501;
t47 = (t161 + t424) * qJD(1) + t376;
t42 = (-qJD(3) * t257 + t313) * t338 + (t425 - t607) * qJD(1) + t493 + t501;
t41 = (t424 + t504) * qJD(1) + t348;
t40 = (t336 * t461 - t338 * t504) * qJD(3) + t488;
t2 = (t108 * t338 + t560 * t336 + (t336 * t503 + t338 * t505) * qJD(1)) * qJD(3) + t558;
t1 = (t523 * t338 + t462 * t336 + (t336 * t460 + t338 * t461) * qJD(1)) * qJD(3) + t558;
t3 = [(-(-t41 + t348 + t606) * t42 + t16 * (-t300 + t486) + t41 * (t429 - t493) + t15 * (t295 + t301 + t427) + t42 * (t275 + t413 + t601) + (t16 * t438 + (-qJD(4) * t42 + t15 * t563) * t337) * t336 + (-t15 * t645 - t42 * qJD(5) + t41 * (rSges(6,1) * qJD(3) - qJD(4)) * t337 + (rSges(6,2) * qJD(3) * t41 + t16 * t466) * t335) * t338) * m(6) + (-(-t47 + t376 + t606) * t48 + t24 * (t486 + t491) + t47 * t429 + t23 * (t427 + t492) + t48 * (t413 + t457) + (t24 * t584 + (-qJD(4) * t48 + t23 * t562) * t337) * t336 + (t47 * (rSges(5,1) * qJD(3) - qJD(4)) * t337 + (t24 * (rSges(5,1) + pkin(3)) + t47 * rSges(5,3) * qJD(3)) * t335) * t338) * m(5) + (t44 * (-t326 + t436) + t64 * t315 + t43 * t643 + (qJD(3) * t259 * t64 + t417 * t44) * t338 + (-rSges(4,2) * t450 - t208 + t456 + t64 + t668) * t65) * m(4) + (t106 * (t336 * t581 + t317 + t572) + t115 * t315 + t105 * t605 + (t489 + t115 + t668) * t116) * m(3) + (((t54 + (-t137 + t397) * t338 + t629 + t644 + t682) * t338 + ((-t137 + t537) * t336 + (t142 - t622 + t640) * t338 - t626 - t673 + t694) * t336) * qJD(3) + t653) * t442 + ((t142 * t333 + ((-t138 + t640) * t336 + t644 + t647) * t336 + ((t142 + t622 + t670) * t338 + t626 + t631) * t338) * qJD(3) + t636 + t649) * t440 + (t632 + t635) * t439 + (t633 + t634 + t637) * t474 / 0.2e1 - (t671 * t338 + (-t628 + t672) * t336) * t467 / 0.2e1 + (t678 * t337 - t679 * t335 - t675 * qJD(3) + (-(t504 - t577) * t42 + (t337 * t42 * t563 + t41 * t438) * t338 + (t42 * t584 + (-t335 * t466 - qJ(2) + t573) * t41) * t336) * m(6) + (-(t161 - t577) * t48 + (t337 * t48 * t562 + t465 * t47) * t338 + (t47 * (-qJ(2) - t416 - t580) + t48 * t465) * t336) * m(5) + (t464 * t565 + (t64 * (-qJ(2) - t417) + t65 * t464) * t336 - (t162 - t577) * t65) * m(4) + (t115 * t581 * t338 + (t115 * (-rSges(3,3) - qJ(2)) - t116 * pkin(1)) * t336 - t254 * t116) * m(3) + t627 * t439) * qJD(1); 0.2e1 * (-t570 / 0.2e1 + t569 / 0.2e1) * m(6) + 0.2e1 * (-t568 / 0.2e1 + t567 / 0.2e1) * m(5) + 0.2e1 * (t564 / 0.2e1 - t566 / 0.2e1) * m(4) + 0.2e1 * (-t105 * t338 / 0.2e1 + t106 * t583) * m(3); (t16 * t529 - t1 * t164 + (t15 * (-t256 - t437) - t1 * t504) * t338 + (t1 * t461 + t16 * t437) * t336 + (-t419 + t519 + (-(-t415 - t579) * qJD(3) - pkin(4) * t475 - t224) * t336) * t41 + (-t420 + t210 + (-(-t249 + t415) * qJD(3) - t165 + t224) * t338) * t42 + (-t506 - (-t485 * t578 - t201 * t338 + (-t257 * t336 - t529) * t336) * qJD(3) + t83 + (qJD(1) * t461 + t523) * t338 + (qJD(1) * t460 + t462) * t336) * t40) * m(6) + (t24 * t529 + t47 * t519 + t48 * t210 - t2 * t164 + t49 * t83 + (t24 * t258 - t47 * t225 + t2 * t505 + t49 * t560 + (t258 * t48 + t49 * t503) * qJD(1)) * t336 + (t23 * (-t256 - t258) + t48 * (-t165 + t225) - t2 * t161 + t49 * t108 + (t258 * t47 + t49 * t505) * qJD(1)) * t338 - t47 * (qJD(1) * t202 + t419) - t48 * (qJD(1) * t197 + t420) - t49 * t506 - ((t48 * (-t249 + t416) - t49 * t202) * t338 + (-t47 * t416 + t49 * (-t529 - t197)) * t336) * qJD(3)) * m(5) + (-(t198 * t65 + t203 * t64) * qJD(1) - (t66 * (-t198 * t336 - t203 * t338) - t414 * t417) * qJD(3) + 0.2e1 * t66 * (t109 * t338 - t112 * t336 + (-t159 * t338 + t162 * t336) * qJD(1)) - t414 * t226 + (t564 - t566 + (t336 * t65 + t565) * qJD(1)) * t259) * m(4) - ((t661 * t335 + t337 * t686) * qJD(3) + (-t657 * t335 + t658 * t337) * qJD(1)) * qJD(1) / 0.2e1 + (t635 * t338 + t634 * t336 + (t628 * t336 + t627 * t338) * qJD(1)) * qJD(1) / 0.2e1 + ((t474 * t620 + t621) * t336 + ((-t625 * t336 + t667) * qJD(3) - t654) * t338) * t442 + ((-t472 * t625 + t621) * t338 + ((t338 * t620 - t667) * qJD(3) + t654) * t336) * t440 + (t633 * qJD(1) + ((t629 * qJD(1) + t618 * t338) * t338 + (t623 * t336 - t673 * qJD(1) + (-t617 + t624) * t338) * t336) * t615) * t583 + (t632 * qJD(1) + ((t631 * qJD(1) + t624 * t338) * t338 + (t617 * t336 - t693 * qJD(1) + (-t618 + t623) * t338) * t336) * t615) * t338 / 0.2e1 - (t637 + t651) * t477 / 0.2e1 + (t636 + t652) * t476 / 0.2e1; -m(5) * (t173 * t47 + t174 * t49 - t48 * t638) - m(6) * (t173 * t41 + t174 * t40 - t42 * t638) + 0.2e1 * ((t47 * t474 - t472 * t48 + t2) * t586 + (t41 * t474 - t42 * t472 + t1) * t585) * t335 + 0.2e1 * ((qJD(3) * t49 - t47 * t476 - t477 * t48 - t567 + t568) * t586 + (qJD(3) * t40 - t41 * t476 - t42 * t477 - t569 + t570) * t585) * t337; m(6) * (-t15 * t336 - t16 * t338);];
tauc = t3(:);
