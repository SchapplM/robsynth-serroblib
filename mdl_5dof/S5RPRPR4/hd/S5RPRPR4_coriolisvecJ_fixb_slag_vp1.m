% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:09
% EndTime: 2019-12-05 17:53:50
% DurationCPUTime: 28.06s
% Computational Cost: add. (16646->711), mult. (13380->887), div. (0->0), fcn. (10340->10), ass. (0->420)
t689 = Icges(4,3) + Icges(5,3);
t323 = qJ(3) + pkin(9);
t314 = sin(t323);
t316 = cos(t323);
t234 = Icges(5,5) * t316 - Icges(5,6) * t314;
t326 = sin(qJ(3));
t328 = cos(qJ(3));
t281 = Icges(4,5) * t328 - Icges(4,6) * t326;
t684 = t234 + t281;
t324 = qJ(1) + pkin(8);
t315 = sin(t324);
t317 = cos(t324);
t688 = -t315 * t684 + t317 * t689;
t522 = t316 * t317;
t527 = t314 * t317;
t547 = Icges(5,6) * t315;
t162 = Icges(5,4) * t522 - Icges(5,2) * t527 + t547;
t266 = Icges(5,4) * t527;
t555 = Icges(5,5) * t315;
t164 = Icges(5,1) * t522 - t266 + t555;
t520 = t317 * t328;
t521 = t317 * t326;
t548 = Icges(4,6) * t315;
t179 = Icges(4,4) * t520 - Icges(4,2) * t521 + t548;
t295 = Icges(4,4) * t521;
t557 = Icges(4,5) * t315;
t181 = Icges(4,1) * t520 - t295 + t557;
t674 = t162 * t314 - t164 * t316 + t179 * t326 - t181 * t328;
t687 = t689 * t315;
t686 = Icges(4,5) * t520 + Icges(5,5) * t522 - Icges(4,6) * t521 - Icges(5,6) * t527 + t687;
t559 = Icges(5,4) * t314;
t235 = Icges(5,2) * t316 + t559;
t560 = Icges(4,4) * t326;
t282 = Icges(4,2) * t328 + t560;
t319 = Icges(4,4) * t328;
t614 = Icges(4,1) * t326 + t319;
t305 = Icges(5,4) * t316;
t616 = Icges(5,1) * t314 + t305;
t683 = t235 * t314 + t282 * t326 - t316 * t616 - t328 * t614;
t233 = Icges(5,5) * t314 + Icges(5,6) * t316;
t185 = t233 * t315;
t280 = Icges(4,5) * t326 + Icges(4,6) * t328;
t200 = t280 * t315;
t685 = t200 + t185;
t528 = t314 * t315;
t265 = Icges(5,4) * t528;
t525 = t315 * t316;
t554 = Icges(5,5) * t317;
t163 = -Icges(5,1) * t525 + t265 + t554;
t524 = t315 * t326;
t294 = Icges(4,4) * t524;
t523 = t315 * t328;
t556 = Icges(4,5) * t317;
t180 = -Icges(4,1) * t523 + t294 + t556;
t636 = -t163 * t522 - t180 * t520 - t315 * t688;
t402 = -Icges(5,2) * t314 + t305;
t161 = Icges(5,6) * t317 - t315 * t402;
t403 = -Icges(4,2) * t326 + t319;
t178 = Icges(4,6) * t317 - t315 * t403;
t635 = t161 * t528 + t178 * t524 + t317 * t688;
t534 = t280 * t317;
t535 = t233 * t317;
t682 = -t534 - t535;
t681 = -t163 * t316 - t180 * t328;
t633 = t161 * t314 + t178 * t326;
t680 = t674 * t317;
t679 = -t163 * t525 - t180 * t523 + t635;
t647 = t162 * t528 - t164 * t525 + t179 * t524 - t181 * t523 + t317 * t686;
t678 = -t161 * t527 - t178 * t521 - t636;
t668 = t315 * t686 - t680;
t646 = -t161 * t316 - t163 * t314 - t178 * t328 - t180 * t326;
t645 = t162 * t316 + t164 * t314 + t179 * t328 + t181 * t326;
t677 = t315 * t683 - t682;
t676 = -t317 * t683 + t685;
t675 = -t633 - t681;
t673 = t683 * qJD(1) + qJD(3) * t684;
t672 = t688 * qJD(1);
t318 = qJ(5) + t323;
t309 = sin(t318);
t310 = cos(t318);
t223 = rSges(6,1) * t309 + rSges(6,2) * t310;
t213 = t402 * qJD(3);
t238 = Icges(5,1) * t316 - t559;
t214 = t238 * qJD(3);
t257 = t403 * qJD(3);
t285 = Icges(4,1) * t328 - t560;
t258 = t285 * qJD(3);
t671 = t213 * t314 - t214 * t316 + t257 * t326 - t258 * t328 + (t235 * t316 + t282 * t328 + t314 * t616 + t326 * t614) * qJD(3) + (-t233 - t280) * qJD(1);
t667 = rSges(5,2) * t316;
t666 = t677 * qJD(1);
t664 = t685 * qJD(3) + (-t317 * t684 - t675 - t687) * qJD(1);
t663 = qJD(1) * t674 + qJD(3) * t682 + t672;
t662 = (t315 * t668 + t317 * t678) * qJD(3);
t661 = (t315 * t647 + t317 * t679) * qJD(3);
t660 = t676 * qJD(1);
t475 = qJD(3) * t317;
t100 = qJD(1) * t161 - t235 * t475;
t190 = t616 * t317;
t102 = -qJD(3) * t190 + (-t238 * t315 + t554) * qJD(1);
t114 = qJD(1) * t178 - t282 * t475;
t205 = t614 * t317;
t116 = -qJD(3) * t205 + (-t285 * t315 + t556) * qJD(1);
t659 = -qJD(1) * t686 + qJD(3) * t645 + t100 * t314 - t102 * t316 + t114 * t326 - t116 * t328;
t476 = qJD(3) * t315;
t101 = t235 * t476 + (-t317 * t402 - t547) * qJD(1);
t189 = t616 * t315;
t103 = qJD(3) * t189 + (-t238 * t317 - t555) * qJD(1);
t115 = t282 * t476 + (-t317 * t403 - t548) * qJD(1);
t204 = t614 * t315;
t117 = qJD(3) * t204 + (-t285 * t317 - t557) * qJD(1);
t658 = qJD(3) * t646 - t101 * t314 + t103 * t316 - t115 * t326 + t117 * t328 + t672;
t531 = t309 * t317;
t254 = rSges(6,2) * t531;
t529 = t310 * t317;
t383 = -rSges(6,1) * t529 - rSges(6,3) * t315;
t153 = -t254 - t383;
t569 = rSges(6,2) * t309;
t570 = rSges(6,1) * t310;
t224 = -t569 + t570;
t608 = -t317 * rSges(6,3) + t224 * t315;
t657 = t153 * t317 + t315 * t608;
t532 = t309 * t315;
t250 = Icges(6,4) * t532;
t530 = t310 * t315;
t552 = Icges(6,5) * t317;
t151 = -Icges(6,1) * t530 + t250 + t552;
t251 = Icges(6,4) * t531;
t553 = Icges(6,5) * t315;
t152 = Icges(6,1) * t529 - t251 + t553;
t322 = qJD(3) + qJD(5);
t231 = t315 * t322;
t232 = t317 * t322;
t301 = Icges(6,4) * t310;
t401 = -Icges(6,2) * t309 + t301;
t617 = Icges(6,1) * t309 + t301;
t631 = t617 + t401;
t340 = qJD(1) * t631 + t231 * (-Icges(6,2) * t529 + t152 - t251) + t232 * (Icges(6,2) * t530 + t151 + t250);
t558 = Icges(6,4) * t309;
t218 = Icges(6,2) * t310 + t558;
t221 = Icges(6,1) * t310 - t558;
t640 = t218 - t221;
t546 = Icges(6,6) * t315;
t150 = Icges(6,4) * t529 - Icges(6,2) * t531 + t546;
t641 = t317 * t617 + t150;
t149 = Icges(6,6) * t317 - t315 * t401;
t642 = -t315 * t617 + t149;
t595 = qJD(1) * t640 + t231 * t641 + t232 * t642;
t656 = t309 * t340 + t310 * t595;
t358 = t315 * (-Icges(5,2) * t522 + t164 - t266) + t317 * (Icges(5,2) * t525 + t163 + t265);
t602 = t315 * (t162 + t190) + t317 * (t161 - t189);
t655 = -t358 * t314 - t316 * t602;
t654 = 0.2e1 * qJD(3);
t653 = t661 + t666;
t652 = t660 + t662;
t651 = qJD(3) * t675 + t101 * t316 + t103 * t314 + t115 * t328 + t117 * t326;
t650 = -qJD(3) * t674 + t100 * t316 + t102 * t314 + t114 * t328 + t116 * t326;
t649 = t315 * t673 - t317 * t671;
t648 = t315 * t671 + t317 * t673;
t581 = pkin(4) * t316;
t583 = pkin(3) * t328;
t416 = -t581 - t583;
t390 = pkin(2) - t416;
t209 = t317 * t390;
t452 = pkin(2) + t583;
t270 = t317 * t452;
t325 = -qJ(4) - pkin(6);
t321 = -pkin(7) + t325;
t488 = t325 - t321;
t123 = -t315 * t488 - t209 + t270;
t644 = t123 - t153;
t543 = Icges(6,3) * t315;
t148 = Icges(6,5) * t529 - Icges(6,6) * t531 + t543;
t643 = t148 * t317 - t152 * t530;
t410 = rSges(3,1) * t315 + rSges(3,2) * t317;
t327 = sin(qJ(1));
t579 = t327 * pkin(1);
t639 = t410 + t579;
t638 = -t686 + t681;
t637 = t223 * t322;
t196 = t224 * t322;
t634 = t196 * t315 - t224 * t231;
t307 = t317 * rSges(5,3);
t494 = rSges(5,2) * t528 + t307;
t632 = t579 - t494;
t490 = t614 + t403;
t491 = t282 - t285;
t630 = (t326 * t490 + t328 * t491) * qJD(1);
t496 = t616 + t402;
t497 = t235 - t238;
t629 = (t314 * t496 + t316 * t497) * qJD(1);
t182 = t223 * t315;
t628 = qJD(1) * t182 + t196 * t317 - t224 * t232;
t572 = rSges(5,1) * t316;
t242 = -rSges(5,2) * t314 + t572;
t627 = t242 * t315;
t575 = pkin(6) + t325;
t624 = t317 * t575;
t241 = rSges(5,1) * t314 + t667;
t584 = pkin(3) * t326;
t435 = t241 + t584;
t621 = t435 * t475;
t585 = pkin(2) * t317;
t146 = t315 * t575 - t270 + t585;
t580 = pkin(6) * t315;
t244 = t580 + t585;
t230 = qJD(1) * t244;
t619 = -qJD(1) * t146 + t230;
t329 = cos(qJ(1));
t578 = t329 * pkin(1);
t618 = -t270 - t578;
t415 = -rSges(3,1) * t317 - t578;
t615 = rSges(3,2) * t315 + t415;
t217 = Icges(6,5) * t310 - Icges(6,6) * t309;
t147 = Icges(6,3) * t317 - t217 * t315;
t516 = -t147 * t315 - t151 * t529;
t613 = t516 - t643;
t308 = t317 * rSges(4,3);
t489 = rSges(4,2) * t524 + t308;
t389 = rSges(4,1) * t523 - t489;
t465 = qJD(1) * t579;
t418 = -qJD(1) * (-pkin(2) * t315 + pkin(6) * t317) + t465;
t288 = rSges(4,1) * t326 + rSges(4,2) * t328;
t450 = t288 * t475;
t94 = qJD(1) * t389 + t418 + t450;
t297 = rSges(4,2) * t521;
t385 = -rSges(4,1) * t520 - rSges(4,3) * t315;
t184 = -t297 - t385;
t225 = t288 * t476;
t436 = -t244 - t578;
t95 = t225 + (-t184 + t436) * qJD(1);
t612 = t315 * t95 + t317 * t94;
t313 = t317 ^ 2;
t611 = (-t315 ^ 2 - t313) * t584;
t500 = -Icges(4,2) * t520 + t181 - t295;
t502 = t179 + t205;
t600 = t326 * t500 + t328 * t502;
t501 = Icges(4,2) * t523 + t180 + t294;
t503 = t178 - t204;
t599 = -t326 * t501 - t328 * t503;
t216 = Icges(6,5) * t309 + Icges(6,6) * t310;
t424 = t640 * t322;
t425 = t631 * t322;
t598 = -qJD(1) * t216 + t309 * t425 + t310 * t424;
t429 = qJD(1) * t149 + t152 * t322 - t218 * t232;
t431 = -(-t221 * t315 + t552) * qJD(1) + t641 * t322;
t597 = -qJD(1) * t148 + t309 * t429 + t310 * t431;
t430 = t151 * t322 + t218 * t231 + (-t317 * t401 - t546) * qJD(1);
t432 = -(-t221 * t317 - t553) * qJD(1) + t642 * t322;
t487 = qJD(1) * t147;
t596 = t309 * t430 + t310 * t432 - t487;
t331 = qJD(1) ^ 2;
t480 = qJD(1) * t315;
t210 = t322 * t480;
t594 = -t210 / 0.2e1;
t211 = qJD(1) * t232;
t593 = t211 / 0.2e1;
t592 = -t231 / 0.2e1;
t591 = t231 / 0.2e1;
t590 = -t232 / 0.2e1;
t589 = t232 / 0.2e1;
t588 = t315 / 0.2e1;
t587 = t317 / 0.2e1;
t586 = -rSges(4,3) - pkin(6);
t582 = pkin(4) * t314;
t577 = -qJD(1) / 0.2e1;
t576 = qJD(1) / 0.2e1;
t573 = rSges(4,1) * t328;
t564 = -t321 + rSges(6,3);
t170 = t216 * t315;
t537 = t216 * t317;
t536 = t223 * t317;
t191 = t241 * t315;
t206 = t288 * t315;
t533 = t288 * t317;
t526 = t315 * t148;
t395 = t218 * t309 - t310 * t617;
t78 = -t317 * t395 + t170;
t519 = t78 * qJD(1);
t517 = t147 * t317 + t149 * t532;
t477 = qJD(1) * t325;
t495 = t317 * t477 + t452 * t480;
t474 = qJD(3) * t326;
t447 = t315 * t474;
t275 = pkin(3) * t447;
t304 = qJD(4) * t317;
t493 = t275 + t304;
t482 = qJD(1) * t234;
t481 = qJD(1) * t281;
t479 = qJD(1) * t317;
t478 = qJD(1) * t321;
t473 = qJD(3) * t328;
t303 = qJD(4) * t315;
t471 = pkin(4) * t525;
t298 = pkin(3) * t524;
t470 = pkin(3) * t523;
t469 = t331 * t578;
t466 = pkin(3) * t474;
t464 = rSges(5,1) * t525;
t268 = rSges(5,2) * t527;
t463 = rSges(5,1) * t522;
t462 = t149 * t531;
t199 = t241 * t476;
t460 = t199 + t493;
t459 = qJD(1) * t254 + t315 * t637;
t449 = t314 * t476;
t458 = rSges(5,1) * t449 + qJD(1) * t268 + t476 * t667;
t456 = -t303 + t495;
t311 = t331 * t579;
t330 = qJD(3) ^ 2;
t423 = qJD(1) * t466;
t455 = t317 * t423 + t330 * t470 + t311;
t454 = rSges(4,1) * t447 + (t315 * t473 + t326 * t479) * rSges(4,2);
t286 = t315 * t477;
t453 = t286 + t493;
t445 = -pkin(2) - t573;
t444 = qJD(1) * t476;
t443 = -t480 / 0.2e1;
t442 = t479 / 0.2e1;
t441 = -t476 / 0.2e1;
t439 = -t475 / 0.2e1;
t434 = t223 + t582;
t343 = qJD(1) * t608;
t91 = -t322 * t536 - t343;
t92 = qJD(1) * t383 + t459;
t433 = -t315 * t92 + t317 * t91;
t253 = pkin(4) * t449;
t422 = -t244 * t331 - t469;
t421 = t146 + t436;
t417 = t582 + t584;
t414 = t572 + t583;
t302 = pkin(2) * t480;
t109 = t302 + (-pkin(6) * qJD(1) - t466) * t317 - t456;
t222 = pkin(6) * t479 - t302;
t413 = -t109 - t222 - t303;
t409 = -rSges(4,2) * t326 + t573;
t408 = qJD(1) * t434;
t376 = t470 + t624;
t352 = t315 * t376;
t342 = qJD(3) * t352 - t146 * t475 + qJD(2);
t373 = t317 * t488 - t471;
t351 = t315 * t373;
t33 = -qJD(3) * t351 - t123 * t475 + t153 * t232 + t231 * t608 + t342;
t407 = t33 * (-t182 * t231 - t232 * t536);
t76 = t150 * t310 + t152 * t309;
t400 = t150 * t309 - t152 * t310;
t388 = t464 - t494;
t386 = t223 * t231 + t253 + t493;
t384 = -rSges(5,3) * t315 - t463;
t192 = t241 * t317;
t140 = qJD(1) * t376;
t382 = t140 - t303 + t418;
t377 = qJD(1) * t388;
t110 = (-pkin(3) * t520 + t580) * qJD(1) + t453;
t374 = t315 * t423 + t422 + (t110 + t304) * qJD(1);
t372 = -t390 - t570;
t367 = qJD(1) * t217 + t170 * t232 - t231 * t537;
t364 = t317 * t377;
t363 = qJD(1) * t400 - t322 * t537 + t487;
t362 = t322 * t170 + (t149 * t309 - t151 * t310 - t217 * t317 - t543) * qJD(1);
t355 = qJD(1) * t395 + t217 * t322;
t350 = qJD(1) * t373;
t10 = t315 * t363 - t317 * t597;
t11 = t315 * t596 + t317 * t362;
t12 = t315 * t597 + t317 * t363;
t56 = -t151 * t530 + t517;
t125 = t150 * t532;
t57 = t125 + t643;
t77 = t315 * t395 + t537;
t74 = t77 * qJD(1);
t21 = t231 * t57 + t232 * t56 + t74;
t58 = -t462 - t516;
t59 = -t317 * t400 + t526;
t22 = t231 * t59 + t232 * t58 + t519;
t36 = t315 * t355 - t317 * t598;
t37 = t315 * t598 + t317 * t355;
t38 = -t309 * t432 + t310 * t430;
t39 = -t309 * t431 + t310 * t429;
t75 = t149 * t310 + t151 * t309;
t9 = t315 * t362 - t317 * t596;
t349 = (qJD(1) * t36 + t10 * t231 - t210 * t58 + t211 * t59 + t232 * t9) * t588 + (-t309 * t595 + t310 * t340) * t577 + t21 * t443 + (qJD(1) * t37 + t11 * t232 + t12 * t231 - t210 * t56 + t211 * t57) * t587 + t22 * t442 + (t10 * t315 + t317 * t9 + (-t315 * t58 + t317 * t59) * qJD(1)) * t591 + (t315 * t57 + t317 * t56) * t594 + (t315 * t59 + t317 * t58) * t593 + (t11 * t317 + t12 * t315 + (-t315 * t56 + t317 * t57) * qJD(1)) * t589 + (t315 * t39 + t317 * t38 + (-t315 * t75 + t317 * t76) * qJD(1)) * t576 + (t367 * t315 - t317 * t656) * t592 + (t315 * t656 + t367 * t317) * t590;
t347 = t317 * t140;
t346 = t317 * t350;
t345 = t184 * t317 + t315 * t389;
t344 = -t372 - t569;
t335 = t150 * t531 - t526 + (-t151 * t315 - t152 * t317) * t310 + t517;
t120 = -qJD(3) * t533 + (-t315 * t409 + t308) * qJD(1);
t121 = qJD(1) * t385 + t454;
t334 = t317 * t120 - t315 * t121 + (-t184 * t315 + t317 * t389) * qJD(1);
t333 = t109 * t317 - t110 * t315 + t146 * t480 + t347;
t332 = qJD(3) * t347 + t109 * t475 - t110 * t476 + t146 * t444;
t300 = rSges(3,2) * t480;
t274 = t315 * t478;
t259 = t409 * qJD(3);
t239 = t417 * qJD(3);
t215 = t242 * qJD(3);
t197 = qJD(1) * t410 + t465;
t169 = (-t417 + t584) * t317;
t168 = t315 * t417 - t298;
t165 = -t268 - t384;
t133 = t317 * t146;
t105 = qJD(1) * t384 + t458;
t104 = -qJD(3) * t192 + (t307 - t627) * qJD(1);
t93 = qJD(3) * t345 + qJD(2);
t73 = t315 * t239 - t479 * t581 + t274 - t275 - t286;
t72 = -t390 * t480 + (-t239 + t466 - t478) * t317 + t495;
t67 = -t259 * t475 + (t121 + t225) * qJD(1) + t422;
t66 = t259 * t476 + t311 + (-t120 - t222 + t450) * qJD(1);
t61 = (-t165 + t421) * qJD(1) + t460;
t60 = t377 + t382 + t621;
t52 = t165 * t475 + t388 * t476 + t342;
t49 = t334 * qJD(3);
t48 = (t421 + t644) * qJD(1) + t386;
t47 = t223 * t232 + t417 * t475 + t343 - t350 + t382;
t41 = (-qJD(3) * t215 - t330 * t583) * t317 + (t105 + t199) * qJD(1) + t374;
t40 = t215 * t476 + (t241 * t475 - t104 + t413) * qJD(1) + t455;
t24 = -t196 * t232 + t210 * t223 + t416 * t330 * t317 + (t73 + t92 + t253) * qJD(1) + t374;
t23 = t330 * t471 + t196 * t231 + t211 * t223 + (t475 * t582 + t413 - t72 - t91) * qJD(1) + t455;
t14 = qJD(3) * t364 + t104 * t475 - t105 * t476 - t165 * t444 + t332;
t5 = -qJD(3) * t346 + t123 * t444 - t153 * t210 + t211 * t608 - t231 * t92 + t232 * t91 + t475 * t72 - t476 * t73 + t332;
t1 = [(t74 + (t59 + t335) * t232 + (t125 - t58 - t462 - t613) * t231) * t592 + m(3) * ((t331 * t410 + t311) * t615 - t197 * t300 + t469 * t639) + (m(3) * (-t197 * t415 + (rSges(3,1) * t479 + qJD(1) * t615 - t300) * t639) + t257 * t328 + t258 * t326 + t213 * t316 + t214 * t314 - t309 * t424 + t310 * t425 - t683 * qJD(3)) * qJD(1) + (t75 + t77) * t594 + (t76 + t78) * t593 + (-t519 + (t57 + (t149 * t317 - t150 * t315) * t309 + t613) * t232 + (-t56 + t335) * t231 + t22) * t590 + (t38 + t37) * t589 + (t23 * (-t209 + t254 - t578) - t48 * t303 - t24 * t579 - t47 * (t274 + t304 + t459) + (-t23 * t570 + t48 * (t239 + t637) + t24 * t564) * t317 + (-t23 * t564 - t47 * t239 - t24 * t344) * t315 + ((t327 * t48 + t329 * t47) * pkin(1) + (-t372 * t47 - t48 * t564) * t317 + (rSges(6,3) * t47 + t344 * t48) * t315) * qJD(1) - (t48 + (t578 - t644) * qJD(1) - t386 + t619) * t47) * m(6) + (t40 * (t268 - t463 + t618) + t41 * (-t317 * t325 - t632) + (t40 * (-rSges(5,3) + t325) + t41 * (-pkin(2) - t414)) * t315 + (t456 + t621 + (t464 + t632) * qJD(1)) * t61 + (-t453 - t458 - t61 + t460 - t619 + (-t384 - t618 - t165 - t578) * qJD(1)) * t60) * m(5) + (t66 * (t297 - t578) + t95 * t302 + t67 * (t489 - t579) - t94 * t454 + (qJD(3) * t288 * t95 + t67 * pkin(6) + t445 * t66) * t317 + (t67 * t445 + t66 * t586) * t315 + ((t327 * t95 + t329 * t94) * pkin(1) + (-t445 * t94 + t586 * t95) * t317 + (t409 * t95 - t586 * t94) * t315) * qJD(1) - (-t225 + t230 + t95 + (t184 + t578) * qJD(1)) * t94) * m(4) + (t39 + t36 + t21) * t591 + ((((t633 - t686) * t317 + t636 + t647) * t317 + (t315 * t638 + t635 - t679) * t315) * qJD(3) + t652 - t660) * t439 + (t649 + t650 + t653) * t476 / 0.2e1 + (((t635 + t668 + t680) * t317 + ((-t633 + t638) * t317 - t636 + t647 - t678) * t315) * qJD(3) + (-t646 + t677) * qJD(1) + t666) * t441 + (t648 + t651 + (t645 + t676) * qJD(1)) * t475 / 0.2e1; m(4) * t49 + m(5) * t14 + m(6) * t5; t349 + (((t315 * t500 + t317 * t501) * t328 + (-t315 * t502 - t317 * t503) * t326 - t314 * t602 + t316 * t358) * qJD(3) + (-t314 * t497 + t316 * t496 - t326 * t491 + t328 * t490) * qJD(1)) * t577 + (t651 * t317 + t650 * t315 + (t315 * t646 + t317 * t645) * qJD(1)) * t576 + ((-t476 * t535 + t482) * t315 + (-t629 + (t315 * t185 + t655) * qJD(3)) * t317 + (-t476 * t534 + t481) * t315 + (-t630 + (t599 * t317 + (t200 - t600) * t315) * qJD(3)) * t317) * t441 + ((t185 * t475 + t482) * t317 + (t629 + (-t317 * t535 - t655) * qJD(3)) * t315 + (t200 * t475 + t481) * t317 + (t630 + (t600 * t315 + (-t534 - t599) * t317) * qJD(3)) * t315) * t439 + (t5 * (-t133 - t351 + t352 + t657) + t23 * (t315 * t434 + t298) + (-t5 * t123 + t24 * (-t223 - t417)) * t317 - t407 + (qJD(1) * t168 - t315 * t408 + t628) * t47 + (t408 * t317 - (-t169 + t536) * qJD(1) + t634) * t48 + (-t315 * t73 + t333 - t346 + t433 + (t343 + t72) * t317 - (-t168 * t315 + t169 * t317 + t611) * qJD(3) + t644 * t480) * t33) * m(6) + (t14 * (t317 * t165 - t133 + (t315 * t414 - t494 + t624) * t315) + t52 * (t317 * t104 - t315 * t105 - t165 * t480 + t333 + t364) + t40 * (t298 + t191) - t41 * t435 * t317 - t60 * (t241 * t480 + (-pkin(3) * t473 - t215) * t317) + t60 * qJD(1) * t191 - (t52 * (-t191 * t315 + t611) + (-t52 * t192 - t60 * (-t242 - t583)) * t317) * qJD(3) + (-qJD(1) * t192 - qJD(3) * t627 + t215 * t315 + t241 * t479) * t61) * m(5) + (-(-t206 * t94 + t533 * t95) * qJD(1) - (t93 * (-t206 * t315 - t317 * t533) + t612 * t409) * qJD(3) + t206 * t66 + t334 * t93 + t345 * t49 - t533 * t67 + (t479 * t95 - t480 * t94) * t288 + t612 * t259) * m(4) + (t649 * qJD(1) + ((t668 * qJD(1) + t658 * t317) * t317 + (t663 * t315 - t678 * qJD(1) + (-t659 + t664) * t317) * t315) * t654) * t588 + (t648 * qJD(1) + ((t647 * qJD(1) + t664 * t317) * t317 + (t659 * t315 - t679 * qJD(1) + (-t658 + t663) * t317) * t315) * t654) * t587 + (t653 + t661) * t443 + (t652 + t662) * t442; m(5) * (t315 * t41 + t317 * t40) + m(6) * (t23 * t317 + t24 * t315); t349 + (t5 * t657 + t33 * ((-t313 * rSges(6,3) + (t224 * t317 - t153) * t315) * qJD(1) + t433) + t23 * t182 - t24 * t536 - t407 + (-qJD(1) * t536 + t223 * t479 + t634) * t48 + (-t223 * t480 + t628) * t47) * m(6);];
tauc = t1(:);
