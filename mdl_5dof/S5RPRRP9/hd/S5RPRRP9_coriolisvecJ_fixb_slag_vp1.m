% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP9_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP9_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:36
% EndTime: 2019-12-31 18:49:08
% DurationCPUTime: 24.15s
% Computational Cost: add. (17681->692), mult. (17149->887), div. (0->0), fcn. (13362->8), ass. (0->401)
t384 = pkin(8) + qJ(3);
t359 = qJ(4) + t384;
t348 = cos(t359);
t347 = sin(t359);
t617 = Icges(5,4) * t347;
t278 = Icges(5,2) * t348 + t617;
t340 = Icges(6,5) * t347;
t446 = Icges(6,3) * t348 - t340;
t740 = t278 + t446;
t613 = Icges(6,5) * t348;
t280 = Icges(6,1) * t347 - t613;
t341 = Icges(5,4) * t348;
t282 = Icges(5,1) * t347 + t341;
t744 = t280 + t282;
t273 = Icges(6,3) * t347 + t613;
t391 = sin(qJ(1));
t392 = cos(qJ(1));
t195 = -Icges(6,6) * t392 + t273 * t391;
t585 = t348 * t391;
t587 = t347 * t391;
t608 = Icges(5,6) * t392;
t201 = Icges(5,4) * t585 - Icges(5,2) * t587 - t608;
t743 = -t195 + t201;
t584 = t348 * t392;
t316 = Icges(6,5) * t584;
t586 = t347 * t392;
t607 = Icges(6,6) * t391;
t196 = Icges(6,3) * t586 + t316 + t607;
t447 = -Icges(5,2) * t347 + t341;
t202 = Icges(5,6) * t391 + t392 * t447;
t742 = -t196 + t202;
t449 = Icges(6,1) * t348 + t340;
t203 = -Icges(6,4) * t392 + t391 * t449;
t317 = Icges(5,4) * t587;
t614 = Icges(5,5) * t392;
t205 = Icges(5,1) * t585 - t317 - t614;
t741 = t203 + t205;
t204 = Icges(6,4) * t391 + t392 * t449;
t283 = Icges(5,1) * t348 - t617;
t206 = Icges(5,5) * t391 + t283 * t392;
t729 = t204 + t206;
t738 = -t446 * t391 + t741;
t737 = (Icges(5,6) - Icges(6,6)) * t348 + (Icges(6,4) + Icges(5,5)) * t347;
t736 = t283 + t449 - t740;
t735 = t273 - t447 - t744;
t734 = t740 * t392 - t729;
t733 = -t282 * t392 - t742;
t732 = t744 * t391 + t743;
t275 = Icges(5,5) * t348 - Icges(5,6) * t347;
t198 = Icges(5,3) * t391 + t275 * t392;
t277 = Icges(6,4) * t348 + Icges(6,6) * t347;
t200 = Icges(6,2) * t391 + t277 * t392;
t731 = -t198 - t200;
t730 = t740 * t347 - t348 * t744;
t728 = t275 + t277;
t443 = t201 * t347 - t205 * t348;
t445 = t195 * t347 + t203 * t348;
t727 = -t443 + t445;
t385 = qJD(3) + qJD(4);
t726 = t736 * t385;
t725 = t735 * t385;
t724 = -t729 * qJD(1) + t732 * t385;
t326 = t385 * t392;
t723 = -t280 * t326 + t733 * t385 + (-t283 * t391 - t203 + t614) * qJD(1);
t325 = t385 * t391;
t722 = -t278 * t325 + t738 * t385 + (-t273 * t392 + t202 - t607) * qJD(1);
t721 = -t734 * t385 + (-t391 * t447 + t195 + t608) * qJD(1);
t720 = t737 * t392;
t719 = t737 * t391;
t718 = t730 * t391 + t720;
t717 = -t730 * t392 + t719;
t714 = rSges(6,3) + qJ(5);
t716 = t731 * qJD(1);
t715 = rSges(6,1) + pkin(4);
t570 = -t196 * t587 - t204 * t585;
t163 = t206 * t585;
t474 = t198 * t392 - t163;
t76 = -t202 * t587 - t474;
t690 = -t200 * t392 - t570 + t76;
t78 = t196 * t586 + t391 * t200 + t204 * t584;
t568 = t391 * t198 + t206 * t584;
t80 = -t202 * t586 + t568;
t688 = t78 + t80;
t284 = pkin(4) * t347 - qJ(5) * t348;
t285 = rSges(6,1) * t347 - rSges(6,3) * t348;
t548 = t284 + t285;
t288 = rSges(6,1) * t348 + rSges(6,3) * t347;
t713 = pkin(4) * t348 + qJ(5) * t347 + t288;
t712 = t737 * qJD(1) + t725 * t347 + t726 * t348;
t711 = -t721 * t347 + t723 * t348 - t716;
t604 = Icges(5,3) * t392;
t197 = Icges(5,5) * t585 - Icges(5,6) * t587 - t604;
t199 = -Icges(6,2) * t392 + t277 * t391;
t662 = qJD(1) * t199;
t710 = -qJD(1) * t197 + t722 * t347 + t724 * t348 - t662;
t709 = t732 * t326 + (-Icges(6,1) * t586 + t316 + t733) * t325 + t736 * qJD(1);
t708 = (-Icges(5,2) * t585 - t317 + t738) * t326 + t734 * t325 + t735 * qJD(1);
t707 = -t730 * qJD(1) - t728 * t385;
t706 = t727 * qJD(1) + t719 * t385 + t716;
t598 = t202 * t347;
t705 = t662 + t720 * t385 + (t196 * t347 + t275 * t391 + t729 * t348 - t598 - t604) * qJD(1);
t183 = t391 * t199;
t77 = t195 * t586 + t203 * t584 + t183;
t704 = t717 * qJD(1) - t326 * t77;
t521 = qJD(5) * t392;
t310 = t347 * t521;
t512 = t348 * t326;
t526 = qJD(1) * t392;
t703 = rSges(6,2) * t526 + t714 * t512 + t310;
t73 = -t199 * t392 + t391 * t445;
t702 = t718 * qJD(1) + t326 * t73;
t701 = t706 * t391 + t392 * t710;
t700 = -t705 * t391 + t392 * t711;
t699 = t710 * t391 - t392 * t706;
t698 = t711 * t391 + t392 * t705;
t428 = t443 * t391;
t600 = t197 * t392;
t75 = -t428 - t600;
t697 = t690 * t325 - t326 * t75 - t702;
t569 = -t391 * t197 - t205 * t584;
t79 = -t201 * t586 - t569;
t696 = t688 * t325 - t326 * t79 + t704;
t695 = -t707 * t391 + t392 * t712;
t694 = t712 * t391 + t392 * t707;
t693 = t724 * t347 - t722 * t348;
t692 = t723 * t347 + t721 * t348;
t691 = t73 + t75;
t689 = t77 + t79;
t687 = t741 * t347 + t743 * t348;
t686 = t729 * t347 + t742 * t348;
t527 = qJD(1) * t391;
t423 = -t326 * t347 - t348 * t527;
t503 = t347 * t527;
t575 = t715 * t423 - t714 * t503 + t703;
t248 = t285 * t391;
t377 = t391 * rSges(6,2);
t484 = pkin(4) * t385 - qJD(5);
t574 = (pkin(4) * t526 + qJ(5) * t325) * t348 + (qJ(5) * t526 - t391 * t484) * t347 - t385 * t248 + (t288 * t392 + t377) * qJD(1);
t381 = t392 * rSges(6,2);
t558 = t713 * t391 - t381;
t557 = t715 * t584 + t714 * t586 + t377;
t679 = t284 * t391 + t248;
t678 = t548 * t392;
t673 = t347 * t708 + t348 * t709;
t672 = t728 * qJD(1) - t720 * t325 + t719 * t326;
t671 = t679 * qJD(1) - t326 * t713 + t348 * t521;
t670 = 2 * qJD(3);
t357 = sin(t384);
t669 = rSges(4,2) * t357;
t463 = qJD(1) * t385;
t302 = t391 * t463;
t358 = cos(t384);
t579 = t358 * (qJD(3) ^ 2);
t517 = pkin(3) * t579;
t523 = qJD(5) * t385;
t424 = t348 * t523 - t517;
t522 = qJD(5) * t391;
t496 = t347 * t522;
t629 = pkin(3) * qJD(3);
t514 = t357 * t629;
t318 = t391 * t514;
t390 = -pkin(6) - qJ(2);
t383 = -pkin(7) + t390;
t338 = t383 * t527;
t342 = t390 * t527;
t389 = cos(pkin(8));
t349 = t389 * pkin(2) + pkin(1);
t304 = pkin(3) * t358 + t349;
t541 = t304 - t349;
t135 = t526 * t541 - t318 - t338 + t342;
t362 = t391 * qJ(2);
t329 = t392 * pkin(1) + t362;
t361 = qJD(2) * t392;
t270 = qJD(1) * t329 - t361;
t635 = pkin(1) - t349;
t564 = t342 - (-t392 * t635 - t362) * qJD(1) - t270;
t509 = -t135 + t564;
t519 = qJD(1) * qJD(2);
t351 = t392 * t519;
t540 = qJD(1) * t318 + t351;
t588 = t347 * t385;
t565 = -qJ(5) * t588 - t288 * t385 - t348 * t484;
t18 = t424 * t392 + t565 * t326 + t548 * t302 + (-t496 + t509 - t574) * qJD(1) + t540;
t291 = t392 * t304;
t334 = t392 * t349;
t532 = -t383 + t390;
t158 = t391 * t532 + t291 - t334;
t466 = -t390 * t391 + t334;
t507 = t158 + t466;
t537 = t318 + t361;
t57 = t496 - t548 * t325 + (t507 + t557) * qJD(1) - t537;
t668 = qJD(1) * t57 + t18;
t364 = t392 * qJ(2);
t578 = t390 * t392;
t432 = t391 * t635 - t578;
t212 = -t364 + t432;
t327 = pkin(1) * t391 - t364;
t308 = qJD(1) * t327;
t666 = qJD(1) * t212 - t308;
t375 = t391 * rSges(4,3);
t580 = t358 * t392;
t582 = t357 * t392;
t232 = rSges(4,1) * t580 - rSges(4,2) * t582 + t375;
t665 = t232 + t466;
t631 = rSges(3,2) * sin(pkin(8));
t634 = rSges(3,1) * t389;
t435 = t391 * rSges(3,3) + (-t631 + t634) * t392;
t664 = t329 + t435;
t346 = Icges(4,4) * t358;
t448 = -Icges(4,2) * t357 + t346;
t296 = Icges(4,1) * t357 + t346;
t546 = -t391 * t304 - t392 * t383;
t157 = t349 * t391 + t546 + t578;
t661 = qJD(1) * t157 + t666;
t286 = rSges(5,1) * t347 + rSges(5,2) * t348;
t249 = t286 * t391;
t253 = t286 * t392;
t632 = rSges(5,1) * t348;
t289 = -rSges(5,2) * t347 + t632;
t209 = rSges(5,1) * t585 - rSges(5,2) * t587 - t392 * rSges(5,3);
t374 = t391 * rSges(5,3);
t211 = rSges(5,1) * t584 - rSges(5,2) * t586 + t374;
t524 = qJD(3) * t392;
t525 = qJD(3) * t391;
t572 = -t157 * t525 + t158 * t524;
t66 = t209 * t325 + t211 * t326 + t572;
t360 = qJD(2) * t391;
t498 = t357 * t524;
t464 = pkin(3) * t498;
t436 = t360 - t464;
t419 = -t286 * t326 + t436;
t556 = t212 - t327;
t508 = t157 + t556;
t71 = (-t209 + t508) * qJD(1) + t419;
t72 = -t286 * t325 + (t211 + t507) * qJD(1) - t537;
t660 = -(qJD(1) * t249 - t326 * t289) * t71 - t66 * (-t325 * t249 - t253 * t326) - t72 * (-qJD(1) * t253 - t289 * t325);
t581 = t358 * t391;
t583 = t357 * t391;
t605 = Icges(4,3) * t392;
t216 = Icges(4,5) * t581 - Icges(4,6) * t583 - t605;
t333 = Icges(4,4) * t583;
t615 = Icges(4,5) * t392;
t220 = Icges(4,1) * t581 - t333 - t615;
t609 = Icges(4,6) * t392;
t218 = Icges(4,4) * t581 - Icges(4,2) * t583 - t609;
t596 = t218 * t357;
t441 = -t220 * t358 + t596;
t83 = -t216 * t392 - t391 * t441;
t293 = Icges(4,5) * t358 - Icges(4,6) * t357;
t292 = Icges(4,5) * t357 + Icges(4,6) * t358;
t425 = qJD(3) * t292;
t618 = Icges(4,4) * t357;
t297 = Icges(4,1) * t358 - t618;
t221 = Icges(4,5) * t391 + t297 * t392;
t219 = Icges(4,6) * t391 + t392 * t448;
t595 = t219 * t357;
t440 = -t221 * t358 + t595;
t655 = -t392 * t425 + (-t293 * t391 + t440 + t605) * qJD(1);
t217 = Icges(4,3) * t391 + t293 * t392;
t529 = qJD(1) * t217;
t654 = qJD(1) * t441 - t391 * t425 + t529;
t294 = Icges(4,2) * t358 + t618;
t437 = t294 * t357 - t296 * t358;
t651 = qJD(1) * t437 + t293 * qJD(3);
t650 = t391 * (-t294 * t392 + t221) - t392 * (-Icges(4,2) * t581 + t220 - t333);
t647 = t302 / 0.2e1;
t303 = t392 * t463;
t646 = t303 / 0.2e1;
t645 = -t325 / 0.2e1;
t644 = t325 / 0.2e1;
t643 = -t326 / 0.2e1;
t642 = t326 / 0.2e1;
t641 = t391 / 0.2e1;
t640 = -t392 / 0.2e1;
t638 = pkin(3) * t357;
t637 = -qJD(1) / 0.2e1;
t636 = qJD(1) / 0.2e1;
t633 = rSges(4,1) * t358;
t630 = rSges(4,2) * t358;
t299 = rSges(4,1) * t357 + t630;
t265 = t299 * t392;
t500 = t299 * t525;
t88 = qJD(1) * t665 - t361 - t500;
t628 = t265 * t88;
t536 = rSges(4,2) * t583 + t392 * rSges(4,3);
t231 = rSges(4,1) * t581 - t536;
t499 = t299 * t524;
t456 = t360 - t499;
t87 = (-t231 + t556) * qJD(1) + t456;
t624 = t391 * t87;
t623 = t71 * t286;
t45 = -qJD(5) * t348 + t325 * t558 + t326 * t557 + t572;
t602 = qJD(1) * t45;
t590 = t292 * t391;
t589 = t292 * t392;
t571 = -t391 * t157 + t392 * t158;
t567 = -t391 * t216 - t220 * t580;
t566 = t391 * t217 + t221 * t580;
t563 = t391 * t209 + t392 * t211;
t552 = t548 * t527;
t352 = qJ(2) * t526;
t533 = t352 + t360;
t551 = qJD(1) * (-pkin(1) * t527 + t533) + t391 * t519;
t545 = -t294 + t297;
t544 = t296 + t448;
t539 = rSges(5,2) * t503 + rSges(5,3) * t526;
t538 = rSges(4,3) * t526 + t527 * t669;
t343 = t391 * t631;
t535 = rSges(3,3) * t526 + qJD(1) * t343;
t534 = t392 * rSges(3,3) + t343;
t528 = qJD(1) * t293;
t339 = qJD(5) * t347;
t102 = -t391 * t437 - t589;
t520 = t102 * qJD(1);
t518 = qJD(1) * qJD(3);
t516 = t391 * t634;
t513 = t358 * t629;
t129 = rSges(5,1) * t423 - rSges(5,2) * t512 + t539;
t131 = -t385 * t249 + (t289 * t392 + t374) * qJD(1);
t511 = t392 * t129 + t391 * t131 + t209 * t526;
t134 = -t464 + (-t391 * t541 + t392 * t532) * qJD(1);
t510 = t392 * t134 + t391 * t135 - t157 * t526;
t506 = qJD(1) * (qJD(1) * t432 - t352) + t551;
t504 = t338 + t537;
t501 = t357 * t526;
t495 = -pkin(1) - t634;
t494 = t392 * t518;
t493 = t527 / 0.2e1;
t492 = t526 / 0.2e1;
t491 = -t525 / 0.2e1;
t488 = t524 / 0.2e1;
t487 = -t286 - t638;
t485 = -t349 - t633;
t483 = t357 * (-t391 ^ 2 - t392 ^ 2);
t174 = t221 * t581;
t473 = t217 * t392 - t174;
t472 = -t197 + t598;
t471 = -t216 + t595;
t462 = qJD(1) * t134 + t506;
t461 = t558 * t391 + t557 * t392;
t460 = -t548 - t638;
t234 = t289 * t385;
t457 = -t234 - t513;
t454 = t633 - t669;
t453 = -t391 * t72 - t392 * t71;
t452 = -t391 * t88 - t392 * t87;
t110 = t218 * t358 + t220 * t357;
t111 = t219 * t358 + t221 * t357;
t434 = -t513 + t565;
t433 = t574 * t391 + t575 * t392 + t558 * t526;
t264 = t299 * t391;
t84 = -t219 * t583 - t473;
t430 = (t391 * t84 - t392 * t83) * qJD(3);
t85 = -t218 * t582 - t567;
t86 = -t219 * t582 + t566;
t429 = (t391 * t86 - t392 * t85) * qJD(3);
t427 = qJD(3) * t296;
t426 = qJD(3) * t294;
t132 = (t231 * t391 + t232 * t392) * qJD(3);
t420 = -t158 * t391 * t518 + t134 * t524 + t135 * t525 - t157 * t494;
t418 = t218 * t392 - t219 * t391;
t417 = (-t357 * t544 + t358 * t545) * qJD(1);
t416 = -t326 * t548 + t310 + t436;
t415 = -t304 - t713;
t141 = qJD(1) * t219 - t391 * t426;
t143 = qJD(1) * t221 - t391 * t427;
t402 = qJD(1) * t216 - qJD(3) * t110 - t141 * t357 + t143 * t358;
t140 = -t392 * t426 + (-t391 * t448 + t609) * qJD(1);
t142 = -t392 * t427 + (-t297 * t391 + t615) * qJD(1);
t401 = -qJD(3) * t111 - t140 * t357 + t142 * t358 + t529;
t268 = t448 * qJD(3);
t269 = t297 * qJD(3);
t400 = qJD(1) * t292 - t268 * t357 + t269 * t358 + (-t294 * t358 - t296 * t357) * qJD(3);
t399 = t45 * (-t679 * t325 - t678 * t326 + t339) + t57 * (-t678 * qJD(1) - t325 * t713 + t348 * t522);
t398 = (t690 * t391 - t691 * t392) * t647 + (t688 * t391 - t689 * t392) * t646 + (t672 * t391 + t673 * t392) * t645 + (t701 * t392 + t700 * t391 + (t689 * t391 + t688 * t392) * qJD(1)) * t644 + (t699 * t392 + t698 * t391 + (t691 * t391 + t690 * t392) * qJD(1)) * t643 + (t673 * t391 - t672 * t392) * t642 + (t695 * qJD(1) + t689 * t302 + t688 * t303 + t700 * t325 + t701 * t326) * t641 + (t694 * qJD(1) + t691 * t302 + t690 * t303 + t698 * t325 + t699 * t326) * t640 + (t347 * t709 - t348 * t708) * t637 + (t693 * t392 + t692 * t391 + (t687 * t391 + t686 * t392) * qJD(1)) * t636 + t697 * t493 + t696 * t492;
t397 = -t357 * t650 + t418 * t358;
t271 = t454 * qJD(3);
t257 = t516 - t534;
t160 = qJD(1) * t664 - t361;
t159 = t360 + (-t257 - t327) * qJD(1);
t145 = -qJD(3) * t264 + (t392 * t454 + t375) * qJD(1);
t144 = -t524 * t630 + (-t358 * t527 - t498) * rSges(4,1) + t538;
t137 = t351 + (-qJD(1) * t435 - t270) * qJD(1);
t136 = qJD(1) * (-qJD(1) * t516 + t535) + t551;
t103 = -t392 * t437 + t590;
t101 = t103 * qJD(1);
t70 = -t271 * t524 + t351 + (-t145 + t500 + t564) * qJD(1);
t69 = -t271 * t525 + (t144 - t499) * qJD(1) + t506;
t65 = -qJD(3) * t440 + t140 * t358 + t142 * t357;
t64 = -t441 * qJD(3) + t141 * t358 + t143 * t357;
t63 = t400 * t391 - t392 * t651;
t62 = t391 * t651 + t400 * t392;
t56 = (t508 - t558) * qJD(1) + t416;
t44 = t101 + t429;
t43 = t430 + t520;
t42 = -t392 * t517 - t234 * t326 + t286 * t302 + (-t131 + t509) * qJD(1) + t540;
t41 = qJD(1) * t129 - t234 * t325 - t286 * t303 + (-t357 * t494 - t391 * t579) * pkin(3) + t462;
t17 = t424 * t391 + t565 * t325 - t548 * t303 + ((-t514 + t339) * t392 + t575) * qJD(1) + t462;
t16 = t129 * t326 + t131 * t325 + t209 * t303 - t211 * t302 + t420;
t9 = -t302 * t557 + t303 * t558 + t325 * t574 + t326 * t575 + t347 * t523 + t420;
t1 = [(t101 + ((t84 - t174 + (t217 + t596) * t392 + t567) * t392 + t566 * t391) * qJD(3)) * t488 + ((t76 + (t201 * t392 + t202 * t391) * t347 + t474 + t569) * t326 + (-t205 * t585 + t600 + t75 + (t201 * t391 - t202 * t392) * t347 + t568 + t78) * t325 + t704) * t642 + (t43 - t520 + ((t392 * t471 - t566 + t86) * t392 + (t391 * t471 + t473 + t85) * t391) * qJD(3)) * t491 + (t65 + t62) * t525 / 0.2e1 + (-qJD(3) * t437 + t268 * t358 + t269 * t357 + t726 * t347 - t725 * t348) * qJD(1) + (t18 * (t381 + t546) + t56 * t504 + t17 * (t291 + t557) + t57 * (t360 + t703) + (t57 * (-t588 * t715 - t514) + (-t57 * t383 + t415 * t56) * qJD(1)) * t392 + (-t17 * t383 + (-t385 * t56 * t714 - t18 * t715) * t348 + (-t18 * t714 + t56 * (rSges(6,1) * t385 + t484)) * t347 + (-t56 * rSges(6,2) + t415 * t57) * qJD(1)) * t391 - (-qJD(1) * t558 + t416 - t56 + t661) * t57) * m(6) + (t42 * (-t209 + t546) + t71 * t504 + t41 * (-t383 * t391 + t211 + t291) + t72 * (t436 + t539) + (-t253 * t72 + t391 * t623) * t385 + ((-t71 * rSges(5,3) + t72 * (-t304 - t632)) * t391 + (t71 * (-t289 - t304) - t72 * t383) * t392) * qJD(1) - (-qJD(1) * t209 + t419 + t661 - t71) * t72) * m(5) + (t70 * (t391 * t485 + t536 - t578) + t87 * (t342 + t361) + t69 * t665 + t88 * (t360 + t538) + (t299 * t624 - t628) * qJD(3) + ((-t87 * rSges(4,3) + t485 * t88) * t391 + (t87 * (-t349 - t454) - t88 * t390) * t392) * qJD(1) - (-qJD(1) * t231 + t456 + t666 - t87) * t88) * m(4) + (t137 * (t391 * t495 + t364 + t534) + t159 * t361 + t136 * t664 + t160 * (t533 + t535) + (t159 * (t495 + t631) * t392 + (t159 * (-rSges(3,3) - qJ(2)) + t160 * t495) * t391) * qJD(1) - (-qJD(1) * t257 - t159 - t308 + t360) * t160) * m(3) - (t64 + t63 + t44) * t524 / 0.2e1 + (t687 - t718) * t647 + (t686 + t717) * t646 + ((t472 * t392 - t428 - t568 + t80) * t326 + (t391 * t472 - t163 - t183 + t570 + (-t727 - t731) * t392 + t689) * t325 + t697 + t702) * t645 + (t692 + t695) * t644 + (-t693 + t694 + t696) * t643 + ((t110 + t102) * t391 + (t111 + t103) * t392) * t518 / 0.2e1; 0.2e1 * (t17 * t640 + t18 * t641) * m(6) + 0.2e1 * (t41 * t640 + t42 * t641) * m(5) + 0.2e1 * (t640 * t69 + t641 * t70) * m(4) + 0.2e1 * (t136 * t640 + t137 * t641) * m(3); ((-t524 * t590 - t528) * t392 + (t417 + (t392 * t589 + t397) * qJD(3)) * t391) * t488 + ((-t525 * t589 + t528) * t391 + (t417 + (t391 * t590 + t397) * qJD(3)) * t392) * t491 + (t391 * t65 - t392 * t64 + (t110 * t391 + t111 * t392) * qJD(1)) * t636 + t398 + ((t357 * t545 + t358 * t544) * qJD(1) + (t418 * t357 + t358 * t650) * qJD(3)) * t637 + (qJD(1) * t62 + (-(t391 * t654 + t402 * t392) * t392 + (t391 * t655 + t401 * t392) * t391 + (t85 * t391 + t86 * t392) * qJD(1)) * t670) * t641 + (qJD(1) * t63 + (-(t402 * t391 - t392 * t654) * t392 + (t401 * t391 - t392 * t655) * t391 + (t83 * t391 + t84 * t392) * qJD(1)) * t670) * t640 + (t43 + t430) * t493 + (t44 + t429) * t492 + (t56 * t552 + t9 * (t461 + t571) + t45 * (t433 + t510) + (t434 * t56 + t460 * t668) * t392 + (t17 * t460 + t57 * t434 + (-t158 - t557) * t602) * t391 - t56 * t671 - (-t57 * t501 + ((-t391 * t57 - t392 * t56) * t358 + t45 * t483) * qJD(3)) * pkin(3) - t399) * m(6) + (t16 * (t563 + t571) + t66 * (t510 + t511) + (t457 * t71 + (qJD(1) * t72 + t42) * t487) * t392 + (t41 * t487 + t72 * t457 + (t623 + t66 * (-t158 - t211)) * qJD(1)) * t391 - (-t72 * t501 + (t358 * t453 + t483 * t66) * qJD(3)) * pkin(3) + t660) * m(5) + (-(t264 * t87 - t628) * qJD(1) - (t132 * (-t264 * t391 - t265 * t392) + t452 * t454) * qJD(3) + 0.2e1 * t132 * (t144 * t392 + t145 * t391 + (t231 * t392 - t232 * t391) * qJD(1)) + t452 * t271 + (-t69 * t391 - t70 * t392 + (-t392 * t88 + t624) * qJD(1)) * t299) * m(4); t398 + (-t399 + t9 * t461 + t45 * t433 + (-t557 * t602 + t565 * t57) * t391 + (-t17 * t391 - t392 * t668) * t548 + (t565 * t392 + t552 - t671) * t56) * m(6) + (t16 * t563 + t66 * (-t211 * t527 + t511) + t453 * t234 + (-t41 * t391 - t42 * t392 + (t391 * t71 - t392 * t72) * qJD(1)) * t286 + t660) * m(5); (t17 * t587 + t18 * t586 + t56 * t512 + (t588 - (t325 * t391 + t326 * t392) * t347) * t45 + (-t326 * t56 - t9) * t348) * m(6);];
tauc = t1(:);
