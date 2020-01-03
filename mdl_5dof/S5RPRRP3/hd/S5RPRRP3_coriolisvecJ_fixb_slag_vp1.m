% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:55
% EndTime: 2020-01-03 11:47:35
% DurationCPUTime: 30.76s
% Computational Cost: add. (17917->661), mult. (15561->823), div. (0->0), fcn. (12116->8), ass. (0->365)
t720 = Icges(5,4) + Icges(6,4);
t723 = Icges(5,1) + Icges(6,1);
t719 = -Icges(5,2) - Icges(6,2);
t716 = Icges(5,6) + Icges(6,6);
t344 = qJ(3) + qJ(4);
t335 = cos(t344);
t722 = t720 * t335;
t334 = sin(t344);
t721 = t720 * t334;
t717 = Icges(5,5) + Icges(6,5);
t718 = t335 * t723 - t721;
t711 = t719 * t334 + t722;
t707 = t719 * t335 - t721;
t343 = qJ(1) + pkin(8);
t332 = sin(t343);
t715 = t716 * t332;
t709 = t723 * t334 + t722;
t714 = Icges(5,3) + Icges(6,3);
t333 = cos(t343);
t695 = t711 * t332 - t716 * t333;
t531 = t333 * t335;
t532 = t333 * t334;
t694 = t720 * t531 + t719 * t532 + t715;
t693 = t718 * t332 - t717 * t333;
t690 = -t716 * t334 + t717 * t335;
t712 = t720 * t532;
t710 = t717 * t332;
t700 = t723 * t531 + t710 - t712;
t708 = t717 * t334 + t716 * t335;
t706 = t709 + t711;
t705 = -t718 - t707;
t704 = t707 * t332 + t693;
t703 = t709 * t333 + t694;
t702 = t709 * t332 + t695;
t701 = t690 * t332 - t333 * t714;
t699 = t334 * t707 + t335 * t709;
t698 = t694 * t334;
t697 = t714 * t332;
t696 = t717 * t531 - t716 * t532 + t697;
t691 = t708 * t333;
t689 = t700 * t335 - t698;
t342 = qJD(3) + qJD(4);
t688 = t705 * t342;
t687 = t706 * t342;
t686 = t702 * t342 + (-t333 * t718 - t710) * qJD(1);
t685 = t693 * qJD(1) + t703 * t342;
t684 = t704 * t342 + (t333 * t711 + t715) * qJD(1);
t530 = t333 * t342;
t683 = -t695 * qJD(1) + t700 * t342 + t707 * t530;
t682 = t708 * t332;
t681 = t693 * t335;
t680 = t695 * t334;
t679 = t332 * t699 - t691;
t678 = t531 * t709 + t707 * t532 + t682;
t536 = t332 * t335;
t677 = t693 * t536;
t676 = t700 * t536;
t675 = t695 * t532;
t674 = t701 * qJD(1);
t537 = t332 * t334;
t649 = -t333 * t701 - t695 * t537 + t677;
t648 = t333 * t696 + t537 * t694 - t676;
t647 = -t332 * t701 - t693 * t531 + t675;
t646 = t332 * t696 + t333 * t689;
t248 = t332 * t342;
t673 = -t704 * t530 + (t719 * t531 + t700 - t712) * t248 + t706 * qJD(1);
t672 = qJD(1) * t699 - t690 * t342;
t671 = -t682 * t342 + (t333 * t690 + t680 - t681 + t697) * qJD(1);
t670 = t689 * qJD(1) + t342 * t691 + t674;
t669 = -t696 * qJD(1) + t334 * t683 + t335 * t685;
t668 = t334 * t684 + t335 * t686 - t674;
t667 = -t708 * qJD(1) + t687 * t334 + t688 * t335;
t631 = -t705 * qJD(1) - t703 * t248 + t702 * t530;
t666 = -t673 * t334 + t335 * t631;
t665 = t679 * qJD(1);
t347 = cos(qJ(3));
t339 = t347 * pkin(3);
t327 = t339 + pkin(2);
t585 = pkin(4) * t335;
t278 = t327 + t585;
t663 = -rSges(6,2) * t532 + t333 * t278;
t349 = -pkin(7) - pkin(6);
t341 = -qJ(5) + t349;
t662 = rSges(6,1) * t536 + t332 * t278 + t333 * t341;
t345 = sin(qJ(3));
t473 = qJD(3) * t345;
t465 = pkin(3) * t473;
t527 = t334 * t342;
t246 = -pkin(4) * t527 - t465;
t476 = qJD(1) * t333;
t459 = t335 * t476;
t477 = qJD(1) * t332;
t661 = rSges(6,1) * t459 + rSges(6,3) * t477 + t332 * t246 + t278 * t476;
t660 = t678 * qJD(1);
t659 = t332 * t671 - t333 * t668;
t658 = t332 * t670 + t333 * t669;
t657 = t332 * t668 + t333 * t671;
t656 = -t332 * t669 + t333 * t670;
t655 = -t648 * t248 - t530 * t649 + t665;
t654 = -t248 * t646 - t530 * t647 - t660;
t653 = t334 * t686 - t335 * t684;
t652 = -t334 * t685 + t335 * t683;
t651 = t332 * t672 + t333 * t667;
t650 = -t332 * t667 + t333 * t672;
t645 = t334 * t693 + t335 * t695;
t644 = -t334 * t700 - t335 * t694;
t270 = t327 * t476;
t526 = t335 * t342;
t381 = -t332 * t526 - t334 * t476;
t462 = t332 * t527;
t471 = qJD(5) * t333;
t486 = t341 - t349;
t570 = -rSges(6,1) * t462 + rSges(6,2) * t381 - t471 - t270 + (-qJD(1) * t486 + t465) * t332 + t661;
t317 = t333 * t349;
t493 = t332 * t327 + t317;
t574 = rSges(6,3) * t333;
t520 = -rSges(6,2) * t537 - t493 - t574 + t662;
t286 = t333 * t327;
t519 = -rSges(6,1) * t531 + t286 - t663 + (-rSges(6,3) + t486) * t332;
t577 = rSges(6,2) * t334;
t580 = rSges(6,1) * t335;
t268 = -t577 + t580;
t227 = t268 * t342;
t635 = pkin(4) * t526 + t227;
t632 = -t268 - t278;
t630 = -t690 * qJD(1) + t248 * t691 - t530 * t682;
t338 = sin(qJ(1)) * pkin(1);
t415 = -t332 * rSges(3,1) - t333 * rSges(3,2);
t629 = t415 - t338;
t336 = Icges(4,4) * t347;
t407 = -Icges(4,2) * t345 + t336;
t616 = Icges(4,1) * t345 + t336;
t489 = t616 + t407;
t566 = Icges(4,4) * t345;
t298 = Icges(4,2) * t347 + t566;
t301 = Icges(4,1) * t347 - t566;
t490 = t298 - t301;
t625 = (t345 * t489 + t347 * t490) * qJD(1);
t614 = -t701 - t689;
t624 = -t332 * (-t696 - t680) - t333 * t614 - t677;
t623 = 0.2e1 * qJD(3);
t588 = pkin(2) * t333;
t167 = t588 - t286 + (pkin(6) + t349) * t332;
t162 = qJD(1) * t167;
t253 = pkin(6) * t332 + t588;
t247 = qJD(1) * t253;
t619 = t247 - t162;
t395 = -t345 * t298 + t347 * t616;
t267 = rSges(5,1) * t334 + rSges(5,2) * t335;
t214 = t267 * t332;
t216 = t267 * t333;
t578 = rSges(5,2) * t334;
t581 = rSges(5,1) * t335;
t269 = -t578 + t581;
t288 = rSges(5,1) * t536;
t575 = rSges(5,3) * t333;
t181 = -rSges(5,2) * t537 + t288 - t575;
t290 = rSges(5,2) * t532;
t183 = rSges(5,1) * t531 + rSges(5,3) * t332 - t290;
t324 = t332 * pkin(2);
t252 = -pkin(6) * t333 + t324;
t166 = -t252 + t493;
t474 = qJD(3) * t333;
t475 = qJD(3) * t332;
t392 = t166 * t475 - t167 * t474 + qJD(2);
t64 = t181 * t248 + t183 * t530 + t392;
t447 = t252 + t338;
t420 = t166 + t447;
t424 = t333 * t465;
t67 = t424 + t530 * t267 + (t181 + t420) * qJD(1);
t340 = cos(qJ(1)) * pkin(1);
t329 = qJD(1) * t340;
t455 = t332 * t473;
t423 = pkin(3) * t455;
t374 = -t248 * t267 + t329 - t423;
t517 = -t167 + t183;
t68 = (t253 + t517) * qJD(1) + t374;
t612 = -(-qJD(1) * t214 + t269 * t530) * t67 - t64 * (-t214 * t248 - t216 * t530) - t68 * (-qJD(1) * t216 - t248 * t269);
t188 = -Icges(4,6) * t333 + t332 * t407;
t190 = -Icges(4,5) * t333 + t301 * t332;
t117 = t188 * t347 + t190 * t345;
t231 = t298 * t332;
t557 = Icges(4,6) * t332;
t128 = -qJD(3) * t231 + (t333 * t407 + t557) * qJD(1);
t233 = t616 * t332;
t563 = Icges(4,5) * t332;
t130 = -qJD(3) * t233 + (t301 * t333 + t563) * qJD(1);
t297 = Icges(4,5) * t347 - Icges(4,6) * t345;
t186 = -Icges(4,3) * t333 + t297 * t332;
t481 = qJD(1) * t186;
t611 = qJD(3) * t117 + t128 * t345 - t130 * t347 - t481;
t272 = t407 * qJD(3);
t273 = t301 * qJD(3);
t296 = Icges(4,5) * t345 + Icges(4,6) * t347;
t610 = -qJD(1) * t296 + qJD(3) * (t298 * t347 + t345 * t616) + t272 * t345 - t273 * t347;
t127 = qJD(1) * t188 + t298 * t474;
t234 = t616 * t333;
t129 = qJD(1) * t190 + qJD(3) * t234;
t528 = t333 * t347;
t529 = t333 * t345;
t554 = Icges(4,3) * t332;
t187 = Icges(4,5) * t528 - Icges(4,6) * t529 + t554;
t189 = Icges(4,4) * t528 - Icges(4,2) * t529 + t557;
t308 = Icges(4,4) * t529;
t191 = Icges(4,1) * t528 - t308 + t563;
t401 = t189 * t347 + t191 * t345;
t609 = -qJD(1) * t187 + qJD(3) * t401 - t127 * t345 + t129 * t347;
t607 = (-t696 + t681) * t333 - t675;
t351 = qJD(1) ^ 2;
t241 = qJD(1) * t248;
t598 = t241 / 0.2e1;
t242 = qJD(1) * t530;
t596 = -t242 / 0.2e1;
t595 = -t248 / 0.2e1;
t594 = t248 / 0.2e1;
t593 = t530 / 0.2e1;
t592 = -t530 / 0.2e1;
t591 = -t332 / 0.2e1;
t590 = -t333 / 0.2e1;
t589 = rSges(4,3) + pkin(6);
t587 = pkin(3) * t345;
t586 = pkin(4) * t334;
t584 = -qJD(1) / 0.2e1;
t583 = qJD(1) / 0.2e1;
t582 = rSges(4,1) * t347;
t579 = rSges(4,2) * t345;
t576 = rSges(4,3) * t333;
t573 = rSges(5,3) - t349;
t572 = rSges(6,3) - t341;
t321 = qJD(5) * t332;
t266 = rSges(6,1) * t334 + rSges(6,2) * t335;
t533 = t333 * t266;
t571 = t342 * t533 - t321 + (-t246 - t465) * t333 + (t486 * t333 - t574 + (-t327 - t632) * t332) * qJD(1);
t310 = rSges(4,2) * t529;
t466 = rSges(4,1) * t528;
t195 = rSges(4,3) * t332 - t310 + t466;
t302 = rSges(4,1) * t345 + rSges(4,2) * t347;
t458 = t302 * t475;
t416 = t329 - t458;
t100 = (t195 + t253) * qJD(1) + t416;
t551 = t100 * t332;
t544 = t188 * t345;
t543 = t189 * t345;
t542 = t190 * t347;
t541 = t530 * t335;
t538 = t296 * t332;
t230 = t296 * t333;
t535 = t332 * t345;
t534 = t332 * t347;
t521 = t520 * t332;
t508 = (-t268 - t585) * t248;
t507 = t188 + t233;
t506 = t189 + t234;
t505 = t190 - t231;
t504 = -Icges(4,2) * t528 + t191 - t308;
t503 = t635 * t333;
t328 = t351 * t340;
t488 = pkin(2) * t476 + pkin(6) * t477;
t501 = qJD(1) * t488 + t328;
t494 = rSges(5,1) * t459 + rSges(5,3) * t477;
t492 = pkin(4) * t532 + t533;
t491 = rSges(4,3) * t477 + qJD(1) * t466;
t478 = qJD(1) * t297;
t472 = qJD(3) * t347;
t132 = t333 * t395 + t538;
t470 = t132 * qJD(1);
t469 = qJD(3) ^ 2 * t339;
t312 = pkin(3) * t529;
t467 = t351 * t338;
t464 = pkin(3) * t472;
t116 = -rSges(5,1) * t462 + rSges(5,2) * t381 + t494;
t463 = t332 * t116 + t181 * t476 - t183 * t477;
t137 = t270 + (-qJD(1) * t349 - t465) * t332 - t488;
t152 = t166 * t476;
t461 = t332 * t137 + t167 * t477 + t152;
t460 = -t167 - t519;
t457 = t302 * t474;
t454 = t332 * t472;
t453 = t477 / 0.2e1;
t452 = -t476 / 0.2e1;
t450 = t475 / 0.2e1;
t448 = t474 / 0.2e1;
t446 = -t267 - t587;
t445 = -t266 - t586;
t213 = t266 * t332;
t431 = -qJD(1) * t213 + t268 * t530;
t425 = qJD(3) * (-t332 ^ 2 - t333 ^ 2);
t422 = qJD(1) * t137 + t333 * t469 + t501;
t417 = rSges(3,1) * t476 - rSges(3,2) * t477;
t291 = -t586 - t587;
t251 = rSges(3,1) * t333 - rSges(3,2) * t332;
t414 = -t579 + t582;
t413 = qJD(1) * t446;
t412 = qJD(1) * t445;
t309 = rSges(4,1) * t534;
t194 = -rSges(4,2) * t535 + t309 - t576;
t99 = t457 + (t194 + t447) * qJD(1);
t411 = t333 * t99 - t551;
t402 = t542 - t544;
t400 = -t191 * t347 + t543;
t399 = t194 * t332 + t195 * t333;
t394 = -t266 + t291;
t393 = t332 * t570 + t476 * t520 + t477 * t519;
t236 = t302 * t333;
t319 = pkin(6) * t476;
t136 = t424 + t319 + (t317 + (-pkin(2) + t327) * t332) * qJD(1);
t245 = pkin(2) * t477 - t319;
t389 = -t136 - t245 - t424;
t388 = qJD(3) * t302;
t387 = (-t332 * t67 - t333 * t68) * qJD(1);
t163 = t190 * t534;
t79 = -t186 * t333 - t188 * t535 + t163;
t164 = t191 * t534;
t80 = t187 * t333 + t189 * t535 - t164;
t386 = (-t332 * t80 - t333 * t79) * qJD(3);
t165 = t188 * t529;
t81 = -t186 * t332 - t190 * t528 + t165;
t82 = t187 * t332 - t400 * t333;
t385 = (-t332 * t82 - t333 * t81) * qJD(3);
t384 = qJD(3) * t152 - t136 * t474 + (t137 + t162) * t475;
t383 = -t423 - t471;
t382 = qJD(1) * t394;
t380 = -t332 * t469 - t467;
t373 = t345 * t505 + t347 * t507;
t372 = t345 * t504 + t347 * t506;
t367 = qJD(1) * t400 - qJD(3) * t230 - t481;
t366 = qJD(3) * t538 + (-t297 * t333 + t402 - t554) * qJD(1);
t363 = t395 * qJD(1) - t297 * qJD(3);
t362 = -rSges(5,1) * t527 - rSges(5,2) * t526 - t465;
t133 = qJD(3) * t236 + (t332 * t414 - t576) * qJD(1);
t134 = -rSges(4,1) * t455 + (-t345 * t476 - t454) * rSges(4,2) + t491;
t360 = -t133 * t333 + t134 * t332 + (t194 * t333 - t195 * t332) * qJD(1);
t355 = t248 * t445 + t329 + t383;
t352 = (-t332 * t648 - t333 * t649) * t598 + (-t332 * t646 - t333 * t647) * t596 + (t659 * t333 + t658 * t332 + (t332 * t647 - t333 * t646) * qJD(1)) * t595 + (t630 * t332 - t333 * t666) * t594 + (t332 * t666 + t630 * t333) * t593 + (t657 * t333 + t656 * t332 + (t332 * t649 - t333 * t648) * qJD(1)) * t592 + (t651 * qJD(1) + t647 * t241 - t646 * t242 + t658 * t248 + t530 * t659) * t591 + (qJD(1) * t650 + t241 * t649 - t242 * t648 + t248 * t656 + t530 * t657) * t590 + (t631 * t334 + t335 * t673) * t584 + (t653 * t333 + t652 * t332 + (t332 * t645 - t333 * t644) * qJD(1)) * t583 + t655 * t453 + t654 * t452;
t276 = t414 * qJD(3);
t235 = t302 * t332;
t228 = t269 * t342;
t193 = -t291 * t333 - t312;
t192 = (t291 + t587) * t332;
t161 = t332 * t181;
t154 = t530 * t533;
t153 = t332 * t166;
t131 = t332 * t395 - t230;
t121 = t131 * qJD(1);
t114 = t342 * t216 + (t269 * t332 - t575) * qJD(1);
t94 = qJD(3) * t399 + qJD(2);
t78 = t276 * t474 + (t134 - t458) * qJD(1) + t501;
t77 = -t467 - t276 * t475 + (-t133 - t245 - t457) * qJD(1);
t66 = -t332 * t610 + t363 * t333;
t65 = t363 * t332 + t333 * t610;
t63 = qJD(3) * t400 + t127 * t347 + t129 * t345;
t62 = qJD(3) * t402 + t128 * t347 + t130 * t345;
t61 = t360 * qJD(3);
t60 = (t253 + t460) * qJD(1) + t355;
t59 = t424 - t321 - t445 * t530 + (t420 + t520) * qJD(1);
t49 = t228 * t530 - t241 * t267 + (t116 - t423) * qJD(1) + t422;
t48 = -t228 * t248 - t242 * t267 + (-t114 + t389) * qJD(1) + t380;
t43 = t385 - t470;
t42 = t121 + t386;
t41 = t248 * t520 - t519 * t530 + t392;
t26 = t227 * t530 - t241 * t266 + (-t241 * t334 + t526 * t530) * pkin(4) + (t383 + t570) * qJD(1) + t422;
t25 = -t227 * t248 - t242 * t266 + (-t242 * t334 - t248 * t526) * pkin(4) + (t389 + t321 - t571) * qJD(1) + t380;
t16 = -t114 * t530 + t116 * t248 + t181 * t242 - t183 * t241 + t384;
t9 = t241 * t519 + t242 * t520 + t248 * t570 - t530 * t571 + t384;
t1 = [(t121 + ((t81 + t164 - t165 + (t186 - t543) * t332) * t332 + (-t163 - t82 + (t186 - t400) * t333 + (t542 + t544) * t332) * t333) * qJD(3)) * t450 + m(3) * ((t351 * t415 - t467) * (t251 + t340) - (t328 + qJD(1) * (t329 + t417)) * t629) + t678 * t242 / 0.2e1 + t644 * t596 + (-(-t624 + t646) * t530 + ((t701 - t698) * t332 + t607 + t647 + t676) * t248 + t665) * t594 + (t43 + t470 + ((t165 - t80 + (t187 - t542) * t333) * t333 + (-t163 + t79 + (t187 + t544) * t332) * t332) * qJD(3)) * t448 + (qJD(3) * t395 + t272 * t347 + t273 * t345 + m(3) * (qJD(1) * t251 + t329 - t417) * t629 + t687 * t335 - t688 * t334) * qJD(1) + (-(-qJD(1) * t519 + t355 - t60 + t619) * t59 + t25 * (t340 + t663) + t60 * t321 + t26 * (t338 + t662) + t59 * (t329 + t661) + (t25 * t580 + t60 * (-rSges(6,1) * t527 - rSges(6,2) * t526 + t246) - t26 * rSges(6,3) - t59 * qJD(5)) * t333 + (-t266 * t342 * t59 + t25 * t572 - t26 * t577) * t332 + (-t60 * t338 + (t572 * t60 - t577 * t59) * t333 + (-t59 * t341 + t60 * t632) * t332) * qJD(1)) * m(6) + (-(qJD(1) * t183 + t374 + t619 - t68) * t67 + t48 * (t286 - t290 + t340) + t49 * (t288 + t338 + t493) + t67 * (t270 + t329 + t494) + (-t49 * rSges(5,3) + t362 * t68 + t48 * t581) * t333 + (t362 * t67 + t48 * t573 - t49 * t578) * t332 + (-t68 * t338 + (t573 * t68 - t578 * t67) * t333 + (t68 * (-t269 - t327) - t67 * t349) * t332) * qJD(1)) * m(5) + (t77 * (-t310 + t340) + t100 * t319 + t78 * (t309 + t324 + t338) + (-t579 * t78 + t589 * t77) * t332 + (t77 * (pkin(2) + t582) - t78 * t589 - t100 * t388) * t333 + (t100 * (t576 - t338) + (-pkin(2) - t414) * t551) * qJD(1) + (-t388 * t332 + t100 - t247 + t329 - t416 + t488 + t491 + (-t195 - t310) * qJD(1)) * t99) * m(4) - (t63 + t65 + t42) * t475 / 0.2e1 + (t645 + t679) * t598 + (-(-t332 * t614 + t607 + t648) * t530 + (t624 + t649) * t248 + t654 + t660) * t593 + (t650 - t653) * t592 + (t651 - t652 + t655) * t595 - (-qJD(1) * t401 + t62 + t66) * t474 / 0.2e1 + (t333 * t132 + (t117 + t131) * t332) * qJD(3) * t583; m(4) * t61 + m(5) * t16 + m(6) * t9; ((-t345 * t490 + t347 * t489) * qJD(1) + ((t332 * t504 - t333 * t505) * t347 + (-t332 * t506 + t333 * t507) * t345) * qJD(3)) * t584 + ((t230 * t475 - t478) * t332 + (t625 + (-t373 * t333 + (-t538 + t372) * t332) * qJD(3)) * t333) * t450 + t352 + (-t332 * t63 - t333 * t62 + (t117 * t332 + t333 * t401) * qJD(1)) * t583 + ((-t474 * t538 - t478) * t333 + (-t625 + (-t372 * t332 + (t230 + t373) * t333) * qJD(3)) * t332) * t448 + (qJD(1) * t65 + (-(t366 * t332 + t333 * t611) * t333 - (t367 * t332 - t333 * t609) * t332 + (t81 * t332 - t82 * t333) * qJD(1)) * t623) * t591 + (qJD(1) * t66 + (-(-t332 * t611 + t366 * t333) * t333 - (t332 * t609 + t367 * t333) * t332 + (t79 * t332 - t80 * t333) * qJD(1)) * t623) * t590 + (t386 + t42) * t453 + (t385 + t43) * t452 + (t26 * (t312 + t492) + t25 * t394 * t332 + (t333 * t460 + t153 + t521) * t9 + (pkin(3) * t454 - t508 - (-t193 - t533 - t312) * qJD(1) + t382 * t333 + (-t635 - t464) * t332) * t60 + (-pkin(4) * t541 - t431 - (-pkin(3) * t535 + t192) * qJD(1) + t503 + t382 * t332) * t59 + (t193 * t530 + t154 - (t192 - t213) * t248 - t425 * t587 + t393 + t461 + (-t136 - t571) * t333) * t41) * m(6) + (-(-t68 * t454 + (t425 * t64 + t387) * t345) * pkin(3) + t16 * (t153 + t161) + t64 * (t461 + t463) + t49 * t312 + (t16 * t517 + t64 * (-t114 - t136) + t49 * t267 + t67 * t228 + t68 * t413) * t333 + (t48 * t446 + t68 * (-t228 - t464) + t67 * t413) * t332 + t612) * m(5) + (-(-t100 * t236 - t235 * t99) * qJD(1) - (t94 * (-t235 * t332 - t236 * t333) + t411 * t414) * qJD(3) + t61 * t399 + t94 * t360 + t411 * t276 + (-t77 * t332 + t78 * t333 + (-t100 * t333 - t332 * t99) * qJD(1)) * t302) * m(4); t352 + (-t41 * (-t213 * t248 - t154) - t60 * (-qJD(1) * t533 + t508) - t59 * t431 - (t59 * t541 + (t41 * (-t248 * t332 - t333 * t530) + (-t332 * t59 - t333 * t60) * qJD(1)) * t334) * pkin(4) + t9 * t521 + t41 * t393 + t26 * t492 + t59 * t503 + (-t41 * t571 + t412 * t60 - t519 * t9) * t333 + (t25 * t445 + t412 * t59 - t60 * t635) * t332) * m(6) + (t16 * (t183 * t333 + t161) + t64 * (-t114 * t333 + t463) + (-t332 * t68 + t333 * t67) * t228 + (-t48 * t332 + t49 * t333 + t387) * t267 + t612) * m(5); 0.2e1 * (t25 * t590 + t26 * t591 - t41 * (-t248 * t333 + t332 * t530) / 0.2e1) * m(6);];
tauc = t1(:);
