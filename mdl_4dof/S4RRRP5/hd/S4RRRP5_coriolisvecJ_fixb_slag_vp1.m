% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:17:02
% DurationCPUTime: 20.02s
% Computational Cost: add. (10708->625), mult. (15213->821), div. (0->0), fcn. (11978->6), ass. (0->366)
t351 = qJ(2) + qJ(3);
t328 = cos(t351);
t327 = sin(t351);
t563 = Icges(4,4) * t327;
t255 = Icges(4,2) * t328 + t563;
t319 = Icges(5,5) * t327;
t408 = Icges(5,3) * t328 - t319;
t681 = t255 + t408;
t559 = Icges(5,5) * t328;
t257 = Icges(5,1) * t327 - t559;
t320 = Icges(4,4) * t328;
t259 = Icges(4,1) * t327 + t320;
t685 = t257 + t259;
t250 = Icges(5,3) * t327 + t559;
t353 = sin(qJ(1));
t355 = cos(qJ(1));
t178 = -Icges(5,6) * t355 + t250 * t353;
t531 = t328 * t353;
t533 = t327 * t353;
t554 = Icges(4,6) * t355;
t184 = Icges(4,4) * t531 - Icges(4,2) * t533 - t554;
t684 = -t178 + t184;
t530 = t328 * t355;
t300 = Icges(5,5) * t530;
t532 = t327 * t355;
t553 = Icges(5,6) * t353;
t179 = Icges(5,3) * t532 + t300 + t553;
t409 = -Icges(4,2) * t327 + t320;
t185 = Icges(4,6) * t353 + t355 * t409;
t683 = -t179 + t185;
t411 = Icges(5,1) * t328 + t319;
t186 = -Icges(5,4) * t355 + t353 * t411;
t301 = Icges(4,4) * t533;
t560 = Icges(4,5) * t355;
t188 = Icges(4,1) * t531 - t301 - t560;
t682 = t186 + t188;
t187 = Icges(5,4) * t353 + t355 * t411;
t260 = Icges(4,1) * t328 - t563;
t189 = Icges(4,5) * t353 + t260 * t355;
t670 = t187 + t189;
t679 = -t408 * t353 + t682;
t678 = (Icges(4,6) - Icges(5,6)) * t328 + (Icges(5,4) + Icges(4,5)) * t327;
t677 = t260 + t411 - t681;
t676 = t250 - t409 - t685;
t675 = -t681 * t355 + t670;
t674 = -t259 * t355 - t683;
t673 = t685 * t353 + t684;
t252 = Icges(4,5) * t328 - Icges(4,6) * t327;
t181 = Icges(4,3) * t353 + t252 * t355;
t254 = Icges(5,4) * t328 + Icges(5,6) * t327;
t183 = Icges(5,2) * t353 + t254 * t355;
t672 = -t181 - t183;
t671 = t681 * t327 - t328 * t685;
t669 = t252 + t254;
t404 = t184 * t327 - t188 * t328;
t406 = t178 * t327 + t186 * t328;
t668 = -t404 + t406;
t348 = qJD(2) + qJD(3);
t667 = t677 * t348;
t666 = t676 * t348;
t665 = -qJD(1) * t670 + t673 * t348;
t280 = t348 * t355;
t664 = -t257 * t280 + t674 * t348 + (-t260 * t353 - t186 + t560) * qJD(1);
t279 = t348 * t353;
t663 = -t255 * t279 + t679 * t348 + (-t250 * t355 + t185 - t553) * qJD(1);
t662 = t675 * t348 + (-t353 * t409 + t178 + t554) * qJD(1);
t661 = t678 * t355;
t660 = t678 * t353;
t659 = t353 * t671 + t661;
t658 = -t355 * t671 + t660;
t655 = rSges(5,3) + qJ(4);
t657 = t672 * qJD(1);
t656 = rSges(5,1) + pkin(3);
t517 = -t179 * t533 - t187 * t531;
t140 = t189 * t531;
t433 = t355 * t181 - t140;
t70 = -t185 * t533 - t433;
t631 = -t183 * t355 - t517 + t70;
t72 = t179 * t532 + t353 * t183 + t187 * t530;
t515 = t353 * t181 + t189 * t530;
t74 = -t185 * t532 + t515;
t629 = t72 + t74;
t261 = pkin(3) * t327 - qJ(4) * t328;
t262 = rSges(5,1) * t327 - rSges(5,3) * t328;
t493 = t261 + t262;
t265 = rSges(5,1) * t328 + rSges(5,3) * t327;
t654 = pkin(3) * t328 + qJ(4) * t327 + t265;
t653 = qJD(1) * t678 + t666 * t327 + t667 * t328;
t550 = Icges(4,3) * t355;
t180 = Icges(4,5) * t531 - Icges(4,6) * t533 - t550;
t182 = -Icges(5,2) * t355 + t254 * t353;
t605 = qJD(1) * t182;
t652 = -qJD(1) * t180 + t663 * t327 + t665 * t328 - t605;
t651 = -t662 * t327 + t664 * t328 - t657;
t650 = t673 * t280 + (-Icges(5,1) * t532 + t300 + t674) * t279 + t677 * qJD(1);
t649 = (-Icges(4,2) * t531 - t301 + t679) * t280 - t675 * t279 + t676 * qJD(1);
t648 = -qJD(1) * t671 - t669 * t348;
t647 = t668 * qJD(1) + t660 * t348 + t657;
t543 = t185 * t327;
t646 = t605 + t661 * t348 + (t179 * t327 + t252 * t353 + t670 * t328 - t543 - t550) * qJD(1);
t163 = t353 * t182;
t71 = t178 * t532 + t186 * t530 + t163;
t645 = t658 * qJD(1) - t280 * t71;
t472 = qJD(4) * t355;
t293 = t327 * t472;
t465 = t328 * t280;
t477 = qJD(1) * t355;
t644 = rSges(5,2) * t477 + t655 * t465 + t293;
t67 = -t355 * t182 + t353 * t406;
t643 = t659 * qJD(1) + t280 * t67;
t642 = t647 * t353 + t652 * t355;
t641 = -t646 * t353 + t651 * t355;
t640 = t652 * t353 - t647 * t355;
t639 = t651 * t353 + t646 * t355;
t392 = t404 * t353;
t524 = t355 * t180;
t69 = -t392 - t524;
t638 = t631 * t279 - t280 * t69 - t643;
t516 = -t353 * t180 - t188 * t530;
t73 = -t184 * t532 - t516;
t637 = t629 * t279 - t280 * t73 + t645;
t636 = -t648 * t353 + t653 * t355;
t635 = t653 * t353 + t648 * t355;
t634 = t665 * t327 - t663 * t328;
t633 = t664 * t327 + t662 * t328;
t632 = t67 + t69;
t630 = t71 + t73;
t628 = t682 * t327 + t684 * t328;
t627 = t670 * t327 + t683 * t328;
t478 = qJD(1) * t353;
t386 = -t280 * t327 - t328 * t478;
t461 = t327 * t478;
t569 = t656 * t386 - t655 * t461 + t644;
t231 = t262 * t353;
t341 = t353 * rSges(5,2);
t443 = pkin(3) * t348 - qJD(4);
t568 = -t348 * t231 + (t265 * t355 + t341) * qJD(1) + (pkin(3) * t477 + qJ(4) * t279) * t328 + (qJ(4) * t477 - t353 * t443) * t327;
t344 = t355 * rSges(5,2);
t501 = t654 * t353 - t344;
t500 = t656 * t530 + t655 * t532 + t341;
t620 = t261 * t353 + t231;
t619 = t493 * t355;
t614 = t649 * t327 + t650 * t328;
t613 = t669 * qJD(1) - t661 * t279 + t660 * t280;
t612 = qJD(1) * t620 - t280 * t654 + t328 * t472;
t611 = 2 * qJD(2);
t352 = sin(qJ(2));
t610 = rSges(3,2) * t352;
t422 = qJD(1) * t348;
t267 = t353 * t422;
t575 = pkin(2) * qJD(2);
t467 = t352 * t575;
t310 = t353 * t467;
t287 = qJD(1) * t310;
t354 = cos(qJ(2));
t525 = t354 * (qJD(2) ^ 2);
t469 = pkin(2) * t525;
t474 = qJD(4) * t348;
t387 = t328 * t474 - t469;
t473 = qJD(4) * t353;
t454 = t327 * t473;
t534 = t327 * t348;
t511 = -qJ(4) * t534 - t265 * t348 - t328 * t443;
t345 = t353 * pkin(5);
t356 = -pkin(6) - pkin(5);
t484 = t356 * t478 + t310;
t321 = pkin(2) * t354 + pkin(1);
t579 = pkin(1) - t321;
t137 = (-t355 * t579 - t345) * qJD(1) - t484;
t296 = t355 * pkin(1) + t345;
t276 = t296 * qJD(1);
t518 = -t137 - t276;
t26 = t287 + t387 * t355 + t511 * t280 + t493 * t267 + (-t454 + t518 - t568) * qJD(1);
t303 = t355 * t321;
t425 = -t353 * t356 + t303;
t177 = t425 - t296;
t463 = -t177 - t500;
t59 = t454 - t310 - t493 * t279 + (t296 - t463) * qJD(1);
t609 = qJD(1) * t59 + t26;
t346 = t355 * pkin(5);
t295 = pkin(1) * t353 - t346;
t325 = t355 * t356;
t486 = -t353 * t321 - t325;
t176 = t295 + t486;
t277 = qJD(1) * t295;
t607 = qJD(1) * t176 - t277;
t338 = Icges(3,4) * t354;
t410 = -Icges(3,2) * t352 + t338;
t285 = Icges(3,1) * t352 + t338;
t263 = rSges(4,1) * t327 + rSges(4,2) * t328;
t232 = t263 * t353;
t236 = t263 * t355;
t577 = rSges(4,1) * t328;
t266 = -rSges(4,2) * t327 + t577;
t192 = rSges(4,1) * t531 - rSges(4,2) * t533 - t355 * rSges(4,3);
t339 = t353 * rSges(4,3);
t194 = rSges(4,1) * t530 - rSges(4,2) * t532 + t339;
t475 = qJD(2) * t355;
t476 = qJD(2) * t353;
t513 = -t176 * t476 + t177 * t475;
t64 = t192 * t279 + t194 * t280 + t513;
t456 = t352 * t475;
t423 = pkin(2) * t456;
t388 = -t263 * t280 - t423;
t507 = t176 - t295;
t65 = (-t192 + t507) * qJD(1) + t388;
t506 = -t177 - t194;
t66 = -t263 * t279 - t310 + (t296 - t506) * qJD(1);
t604 = -t65 * (qJD(1) * t232 - t280 * t266) - t64 * (-t279 * t232 - t236 * t280) - t66 * (-qJD(1) * t236 - t266 * t279);
t527 = t353 * t354;
t529 = t352 * t353;
t551 = Icges(3,3) * t355;
t210 = Icges(3,5) * t527 - Icges(3,6) * t529 - t551;
t315 = Icges(3,4) * t529;
t561 = Icges(3,5) * t355;
t214 = Icges(3,1) * t527 - t315 - t561;
t555 = Icges(3,6) * t355;
t212 = Icges(3,4) * t527 - Icges(3,2) * t529 - t555;
t542 = t212 * t352;
t402 = -t214 * t354 + t542;
t77 = -t355 * t210 - t353 * t402;
t282 = Icges(3,5) * t354 - Icges(3,6) * t352;
t281 = Icges(3,5) * t352 + Icges(3,6) * t354;
t389 = qJD(2) * t281;
t564 = Icges(3,4) * t352;
t286 = Icges(3,1) * t354 - t564;
t215 = Icges(3,5) * t353 + t286 * t355;
t213 = Icges(3,6) * t353 + t355 * t410;
t541 = t213 * t352;
t401 = -t215 * t354 + t541;
t599 = -t355 * t389 + (-t282 * t353 + t401 + t551) * qJD(1);
t211 = Icges(3,3) * t353 + t282 * t355;
t480 = qJD(1) * t211;
t598 = qJD(1) * t402 - t353 * t389 + t480;
t283 = Icges(3,2) * t354 + t564;
t398 = t352 * t283 - t285 * t354;
t595 = t398 * qJD(1) + t282 * qJD(2);
t594 = t353 * (-t283 * t355 + t215) - t355 * (-Icges(3,2) * t527 + t214 - t315);
t591 = t267 / 0.2e1;
t268 = t355 * t422;
t590 = t268 / 0.2e1;
t589 = -t279 / 0.2e1;
t588 = t279 / 0.2e1;
t587 = -t280 / 0.2e1;
t586 = t280 / 0.2e1;
t585 = t353 / 0.2e1;
t584 = -t355 / 0.2e1;
t582 = pkin(2) * t352;
t581 = -qJD(1) / 0.2e1;
t580 = qJD(1) / 0.2e1;
t578 = rSges(3,1) * t354;
t576 = rSges(3,2) * t354;
t340 = t353 * rSges(3,3);
t571 = t65 * t263;
t51 = -qJD(4) * t328 + t279 * t501 + t280 * t500 + t513;
t548 = qJD(1) * t51;
t483 = rSges(3,2) * t529 + t355 * rSges(3,3);
t216 = rSges(3,1) * t527 - t483;
t288 = rSges(3,1) * t352 + t576;
t457 = t288 * t475;
t117 = -t457 + (-t216 - t295) * qJD(1);
t546 = t117 * t353;
t545 = t117 * t355;
t458 = t288 * t476;
t526 = t354 * t355;
t528 = t352 * t355;
t217 = rSges(3,1) * t526 - rSges(3,2) * t528 + t340;
t497 = t217 + t296;
t118 = qJD(1) * t497 - t458;
t247 = t288 * t355;
t544 = t118 * t247;
t536 = t281 * t353;
t535 = t281 * t355;
t326 = pkin(5) * t477;
t136 = -t423 - t326 + (t353 * t579 - t325) * qJD(1);
t248 = qJD(1) * (-pkin(1) * t478 + t326);
t519 = qJD(1) * t136 + t248;
t512 = -t353 * t176 + t355 * t177;
t510 = t353 * t192 + t355 * t194;
t509 = -t353 * t210 - t214 * t526;
t508 = t353 * t211 + t215 * t526;
t496 = t493 * t478;
t489 = rSges(4,2) * t461 + rSges(4,3) * t477;
t488 = -t283 + t286;
t487 = t285 + t410;
t485 = rSges(3,3) * t477 + t478 * t610;
t479 = qJD(1) * t282;
t318 = qJD(4) * t327;
t121 = -t353 * t398 - t535;
t471 = t121 * qJD(1);
t470 = qJD(1) * qJD(2);
t112 = rSges(4,1) * t386 - rSges(4,2) * t465 + t489;
t114 = -t348 * t232 + (t266 * t355 + t339) * qJD(1);
t468 = t355 * t112 + t353 * t114 + t192 * t477;
t466 = t354 * t575;
t464 = t355 * t136 + t353 * t137 - t176 * t477;
t459 = t352 * t477;
t453 = -pkin(1) - t578;
t452 = t355 * t470;
t451 = t478 / 0.2e1;
t450 = t477 / 0.2e1;
t449 = -t476 / 0.2e1;
t446 = t475 / 0.2e1;
t445 = -t263 - t582;
t441 = t352 * (-t353 ^ 2 - t355 ^ 2);
t170 = t215 * t527;
t432 = t355 * t211 - t170;
t431 = -t180 + t543;
t426 = -t210 + t541;
t421 = t353 * t501 + t355 * t500;
t420 = -t493 - t582;
t209 = t266 * t348;
t417 = -t209 - t466;
t415 = t578 - t610;
t414 = -t353 * t66 - t355 * t65;
t407 = -t118 * t353 - t545;
t119 = t212 * t354 + t214 * t352;
t120 = t213 * t354 + t215 * t352;
t397 = t568 * t353 + t569 * t355 + t501 * t477;
t396 = -t466 + t511;
t246 = t288 * t353;
t78 = -t213 * t529 - t432;
t394 = (t353 * t78 - t355 * t77) * qJD(2);
t79 = -t212 * t528 - t509;
t80 = -t213 * t528 + t508;
t393 = (t353 * t80 - t355 * t79) * qJD(2);
t391 = qJD(2) * t285;
t390 = qJD(2) * t283;
t115 = (t216 * t353 + t217 * t355) * qJD(2);
t383 = -t177 * t353 * t470 + t136 * t475 + t137 * t476 - t176 * t452;
t382 = t212 * t355 - t213 * t353;
t381 = -t280 * t493 + t293 - t423;
t380 = (-t352 * t487 + t354 * t488) * qJD(1);
t379 = -t321 - t654;
t130 = qJD(1) * t213 - t353 * t390;
t132 = qJD(1) * t215 - t353 * t391;
t366 = qJD(1) * t210 - qJD(2) * t119 - t130 * t352 + t132 * t354;
t129 = -t355 * t390 + (-t353 * t410 + t555) * qJD(1);
t131 = -t355 * t391 + (-t286 * t353 + t561) * qJD(1);
t365 = -qJD(2) * t120 - t129 * t352 + t131 * t354 + t480;
t270 = t410 * qJD(2);
t271 = t286 * qJD(2);
t364 = qJD(1) * t281 - t270 * t352 + t271 * t354 + (-t283 * t354 - t285 * t352) * qJD(2);
t363 = t51 * (-t279 * t620 - t280 * t619 + t318) + t59 * (-qJD(1) * t619 - t279 * t654 + t328 * t473);
t362 = (t631 * t353 - t632 * t355) * t591 + (t629 * t353 - t630 * t355) * t590 + (t353 * t613 + t355 * t614) * t589 + (t642 * t355 + t641 * t353 + (t630 * t353 + t629 * t355) * qJD(1)) * t588 + (t640 * t355 + t639 * t353 + (t632 * t353 + t631 * t355) * qJD(1)) * t587 + (t353 * t614 - t355 * t613) * t586 + (t636 * qJD(1) + t630 * t267 + t629 * t268 + t641 * t279 + t642 * t280) * t585 + (t635 * qJD(1) + t632 * t267 + t631 * t268 + t639 * t279 + t640 * t280) * t584 + (t650 * t327 - t649 * t328) * t581 + (t634 * t355 + t633 * t353 + (t628 * t353 + t627 * t355) * qJD(1)) * t580 + t638 * t451 + t637 * t450;
t361 = -t352 * t594 + t382 * t354;
t274 = t415 * qJD(2);
t134 = -qJD(2) * t246 + (t355 * t415 + t340) * qJD(1);
t133 = -t475 * t576 + (-t354 * t478 - t456) * rSges(3,1) + t485;
t122 = -t355 * t398 + t536;
t116 = t122 * qJD(1);
t76 = -t274 * t475 + (-t134 - t276 + t458) * qJD(1);
t75 = -t274 * t476 + t248 + (t133 - t457) * qJD(1);
t63 = t364 * t353 - t355 * t595;
t62 = t353 * t595 + t364 * t355;
t61 = -qJD(2) * t401 + t129 * t354 + t131 * t352;
t60 = -t402 * qJD(2) + t130 * t354 + t132 * t352;
t58 = (-t501 + t507) * qJD(1) + t381;
t53 = -t355 * t469 - t209 * t280 + t263 * t267 + t287 + (-t114 + t518) * qJD(1);
t52 = qJD(1) * t112 - t209 * t279 - t263 * t268 + (-t352 * t452 - t353 * t525) * pkin(2) + t519;
t42 = t116 + t393;
t41 = t394 + t471;
t25 = t387 * t353 + t511 * t279 - t493 * t268 + ((-t467 + t318) * t355 + t569) * qJD(1) + t519;
t16 = t112 * t280 + t114 * t279 + t192 * t268 - t194 * t267 + t383;
t9 = -t267 * t500 + t268 * t501 + t279 * t568 + t280 * t569 + t327 * t474 + t383;
t1 = [(t116 + ((t78 - t170 + (t211 + t542) * t355 + t509) * t355 + t508 * t353) * qJD(2)) * t446 + ((t70 + (t184 * t355 + t185 * t353) * t327 + t433 + t516) * t280 + (-t188 * t531 + t524 + t69 + (t184 * t353 - t185 * t355) * t327 + t515 + t72) * t279 + t645) * t586 + (-t471 + ((t355 * t426 - t508 + t80) * t355 + (t353 * t426 + t432 + t79) * t353) * qJD(2) + t41) * t449 + (t61 + t62) * t476 / 0.2e1 + (-qJD(2) * t398 + t270 * t354 + t271 * t352 + t667 * t327 - t666 * t328) * qJD(1) + (t26 * (t344 + t486) + t58 * t484 + t25 * (t303 + t500) + t59 * t644 + (t59 * (-t534 * t656 - t467) + (-t59 * t356 + t379 * t58) * qJD(1)) * t355 + (-t25 * t356 + (-t348 * t58 * t655 - t26 * t656) * t328 + (-t26 * t655 + t58 * (rSges(5,1) * t348 + t443)) * t327 + (-t58 * rSges(5,2) + t379 * t59) * qJD(1)) * t353 - (-qJD(1) * t501 + t381 - t58 + t607) * t59) * m(5) + (t53 * (-t192 + t486) + t65 * t484 + t52 * (t194 + t425) + t66 * (-t423 + t489) + (-t236 * t66 + t353 * t571) * t348 + ((-t65 * rSges(4,3) + t66 * (-t321 - t577)) * t353 + (t65 * (-t266 - t321) - t66 * t356) * t355) * qJD(1) - (-qJD(1) * t192 + t388 + t607 - t65) * t66) * m(4) + (-(-qJD(1) * t216 - t117 - t277 - t457) * t118 + t76 * (t353 * t453 + t346 + t483) + t75 * t497 + t118 * (t326 + t485) + (t288 * t546 - t544) * qJD(2) + ((-pkin(1) - t415) * t545 + (t117 * (-rSges(3,3) - pkin(5)) + t118 * t453) * t353) * qJD(1)) * m(3) - (t60 + t63 + t42) * t475 / 0.2e1 + (t628 - t659) * t591 + (t627 + t658) * t590 + ((t431 * t355 - t392 - t515 + t74) * t280 + (t431 * t353 - t140 - t163 + t517 + (-t668 - t672) * t355 + t630) * t279 + t638 + t643) * t589 + (t633 + t636) * t588 + (-t634 + t635 + t637) * t587 + ((t119 + t121) * t353 + (t120 + t122) * t355) * t470 / 0.2e1; t362 + ((-t475 * t536 - t479) * t355 + (t380 + (t355 * t535 + t361) * qJD(2)) * t353) * t446 + ((-t476 * t535 + t479) * t353 + (t380 + (t353 * t536 + t361) * qJD(2)) * t355) * t449 + ((t352 * t488 + t354 * t487) * qJD(1) + (t382 * t352 + t354 * t594) * qJD(2)) * t581 + (t353 * t61 - t355 * t60 + (t119 * t353 + t120 * t355) * qJD(1)) * t580 + (qJD(1) * t62 + (-(t353 * t598 + t366 * t355) * t355 + (t353 * t599 + t365 * t355) * t353 + (t79 * t353 + t80 * t355) * qJD(1)) * t611) * t585 + (qJD(1) * t63 + (-(t366 * t353 - t355 * t598) * t355 + t353 * (t365 * t353 - t355 * t599) + (t77 * t353 + t78 * t355) * qJD(1)) * t611) * t584 + (t394 + t41) * t451 + (t393 + t42) * t450 + (t58 * t496 + t9 * (t421 + t512) + t51 * (t397 + t464) + (t396 * t58 + t420 * t609) * t355 + (t25 * t420 + t396 * t59 + t463 * t548) * t353 - t58 * t612 - (-t59 * t459 + ((-t353 * t59 - t355 * t58) * t354 + t51 * t441) * qJD(2)) * pkin(2) - t363) * m(5) + (t16 * (t510 + t512) + t64 * (t464 + t468) + (t417 * t65 + (qJD(1) * t66 + t53) * t445) * t355 + (t52 * t445 + t66 * t417 + (t506 * t64 + t571) * qJD(1)) * t353 - (-t66 * t459 + (t354 * t414 + t441 * t64) * qJD(2)) * pkin(2) + t604) * m(4) + (0.2e1 * t115 * (t355 * t133 + t353 * t134 + (t216 * t355 - t217 * t353) * qJD(1)) + t407 * t274 + (-t75 * t353 - t76 * t355 + (-t118 * t355 + t546) * qJD(1)) * t288 - (t117 * t246 - t544) * qJD(1) - (t115 * (-t246 * t353 - t247 * t355) + t407 * t415) * qJD(2)) * m(3); t362 + (t9 * t421 + t51 * t397 + (-t500 * t548 + t511 * t59) * t353 - t363 + (-t25 * t353 - t355 * t609) * t493 + (t511 * t355 + t496 - t612) * t58) * m(5) + (t16 * t510 + t64 * (-t194 * t478 + t468) + t414 * t209 + (-t52 * t353 - t53 * t355 + (t353 * t65 - t355 * t66) * qJD(1)) * t263 + t604) * m(4); (t25 * t533 + t26 * t532 + t58 * t465 + (t534 - (t279 * t353 + t280 * t355) * t327) * t51 + (-t280 * t58 - t9) * t328) * m(5);];
tauc = t1(:);
