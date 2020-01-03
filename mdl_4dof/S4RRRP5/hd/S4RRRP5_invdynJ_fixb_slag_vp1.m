% Calculate vector of inverse dynamics joint torques for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:17:04
% DurationCPUTime: 22.34s
% Computational Cost: add. (11556->699), mult. (16171->888), div. (0->0), fcn. (12814->6), ass. (0->388)
t377 = qJ(2) + qJ(3);
t352 = sin(t377);
t353 = cos(t377);
t593 = Icges(5,5) * t353;
t269 = Icges(5,1) * t352 - t593;
t338 = Icges(4,4) * t353;
t271 = Icges(4,1) * t352 + t338;
t726 = t269 + t271;
t597 = Icges(4,4) * t352;
t267 = Icges(4,2) * t353 + t597;
t337 = Icges(5,5) * t352;
t435 = Icges(5,3) * t353 - t337;
t725 = t435 + t267;
t262 = Icges(5,3) * t352 + t593;
t379 = sin(qJ(1));
t381 = cos(qJ(1));
t187 = -Icges(5,6) * t381 + t262 * t379;
t565 = t353 * t379;
t567 = t352 * t379;
t588 = Icges(4,6) * t381;
t193 = Icges(4,4) * t565 - Icges(4,2) * t567 - t588;
t729 = -t187 + t193;
t564 = t353 * t381;
t315 = Icges(5,5) * t564;
t566 = t352 * t381;
t587 = Icges(5,6) * t379;
t188 = Icges(5,3) * t566 + t315 + t587;
t436 = -Icges(4,2) * t352 + t338;
t194 = Icges(4,6) * t379 + t381 * t436;
t728 = -t188 + t194;
t438 = Icges(5,1) * t353 + t337;
t195 = -Icges(5,4) * t381 + t379 * t438;
t316 = Icges(4,4) * t567;
t594 = Icges(4,5) * t381;
t197 = Icges(4,1) * t565 - t316 - t594;
t727 = t195 + t197;
t196 = Icges(5,4) * t379 + t381 * t438;
t272 = Icges(4,1) * t353 - t597;
t198 = Icges(4,5) * t379 + t272 * t381;
t715 = t196 + t198;
t724 = -t435 * t379 + t727;
t723 = (Icges(4,6) - Icges(5,6)) * t353 + (Icges(5,4) + Icges(4,5)) * t352;
t722 = t272 + t438 - t725;
t721 = t262 - t436 - t726;
t720 = t725 * t381 - t715;
t719 = -t271 * t381 - t728;
t718 = t726 * t379 + t729;
t264 = Icges(4,5) * t353 - Icges(4,6) * t352;
t190 = Icges(4,3) * t379 + t264 * t381;
t266 = Icges(5,4) * t353 + Icges(5,6) * t352;
t192 = Icges(5,2) * t379 + t266 * t381;
t717 = -t190 - t192;
t716 = -t725 * t352 + t726 * t353;
t714 = t264 + t266;
t430 = t193 * t352 - t197 * t353;
t432 = t187 * t352 + t195 * t353;
t713 = -t430 + t432;
t374 = qJD(2) + qJD(3);
t712 = t722 * t374;
t711 = t721 * t374;
t710 = -t715 * qJD(1) + t718 * t374;
t294 = t374 * t381;
t709 = -t269 * t294 + t719 * t374 + (-t272 * t379 - t195 + t594) * qJD(1);
t293 = t374 * t379;
t708 = -t267 * t293 + t724 * t374 + (-t262 * t381 + t194 - t587) * qJD(1);
t707 = -t720 * t374 + (-t379 * t436 + t187 + t588) * qJD(1);
t706 = t723 * t381;
t705 = t723 * t379;
t675 = t716 * t379 - t706;
t674 = t716 * t381 + t705;
t702 = rSges(5,3) + qJ(4);
t648 = t353 * rSges(5,1) + t352 * rSges(5,3);
t704 = -t352 * qJ(4) - t648;
t703 = t717 * qJD(1);
t619 = rSges(5,1) + pkin(3);
t552 = -t188 * t567 - t196 * t565;
t147 = t198 * t565;
t467 = t381 * t190 - t147;
t74 = -t194 * t567 - t467;
t678 = -t192 * t381 - t552 + t74;
t76 = t188 * t566 + t379 * t192 + t196 * t564;
t550 = t379 * t190 + t198 * t564;
t78 = -t194 * t566 + t550;
t676 = t76 + t78;
t701 = t723 * qJD(1) + t711 * t352 + t712 * t353;
t700 = -t707 * t352 + t709 * t353 - t703;
t584 = Icges(4,3) * t381;
t189 = Icges(4,5) * t565 - Icges(4,6) * t567 - t584;
t191 = -Icges(5,2) * t381 + t266 * t379;
t643 = qJD(1) * t191;
t699 = -qJD(1) * t189 + t708 * t352 + t710 * t353 - t643;
t698 = t718 * t294 + (-Icges(5,1) * t566 + t315 + t719) * t293 + t722 * qJD(1);
t697 = (-Icges(4,2) * t565 - t316 + t724) * t294 + t720 * t293 + t721 * qJD(1);
t696 = t716 * qJD(1) - t714 * t374;
t695 = t713 * qJD(1) + t705 * t374 + t703;
t577 = t194 * t352;
t694 = t643 + t706 * t374 + (t188 * t352 + t264 * t379 + t715 * t353 - t577 - t584) * qJD(1);
t693 = t353 * pkin(3) - t704;
t171 = t379 * t191;
t75 = t187 * t566 + t195 * t564 + t171;
t692 = qJD(1) * t674 - t294 * t75;
t312 = qJ(4) * t564;
t498 = t353 * t294;
t510 = qJD(1) * t381;
t691 = rSges(5,2) * t510 + rSges(5,3) * t498 + t374 * t312;
t521 = t702 * t565;
t520 = rSges(5,3) * t564 + t312;
t71 = -t381 * t191 + t379 * t432;
t690 = -qJD(1) * t675 + t294 * t71;
t689 = t695 * t379 + t699 * t381;
t688 = -t694 * t379 + t700 * t381;
t687 = t699 * t379 - t695 * t381;
t686 = t700 * t379 + t694 * t381;
t416 = t430 * t379;
t558 = t381 * t189;
t73 = -t416 - t558;
t685 = t678 * t293 - t294 * t73 - t690;
t551 = -t379 * t189 - t197 * t564;
t77 = -t193 * t566 - t551;
t684 = t676 * t293 - t294 * t77 + t692;
t683 = -t696 * t379 + t701 * t381;
t682 = t701 * t379 + t696 * t381;
t681 = t710 * t352 - t708 * t353;
t680 = t709 * t352 + t707 * t353;
t679 = t71 + t73;
t677 = t75 + t77;
t673 = t727 * t352 + t729 * t353;
t672 = t715 * t352 + t728 * t353;
t378 = sin(qJ(2));
t671 = rSges(3,2) * t378;
t505 = qJD(4) * t381;
t306 = t352 * t505;
t511 = qJD(1) * t379;
t410 = -t294 * t352 - t353 * t511;
t491 = t352 * t511;
t603 = t619 * t410 - t702 * t491 + t306 + t691;
t274 = rSges(5,1) * t352 - rSges(5,3) * t353;
t366 = t379 * rSges(5,2);
t477 = pkin(3) * t374 - qJD(4);
t602 = -t274 * t293 + (t381 * t648 + t366) * qJD(1) + (pkin(3) * t510 + qJ(4) * t293) * t353 + (qJ(4) * t510 - t379 * t477) * t352;
t503 = qJD(1) * qJD(2);
t342 = t379 * t503;
t502 = qJD(1) * qJD(3);
t207 = t379 * t502 + t342 + (-qJDD(2) - qJDD(3)) * t381;
t289 = -qJDD(2) * t381 + t342;
t618 = pkin(2) * t378;
t251 = t289 * t618;
t380 = cos(qJ(2));
t559 = t380 * qJD(2) ^ 2;
t500 = pkin(2) * t559;
t507 = qJD(4) * t374;
t405 = qJDD(4) * t352 + t353 * t507 - t500;
t369 = t381 * rSges(5,2);
t536 = t693 * t379 - t369;
t372 = t381 * pkin(5);
t308 = pkin(1) * t379 - t372;
t371 = t380 * pkin(2);
t345 = t371 + pkin(1);
t382 = -pkin(6) - pkin(5);
t349 = t381 * t382;
t519 = -t379 * t345 - t349;
t185 = t308 + t519;
t542 = t185 - t308;
t455 = -t536 + t542;
t506 = qJD(4) * t379;
t485 = t352 * t506;
t527 = pkin(3) * t352 - qJ(4) * t353 + t274;
t546 = -t353 * t477 + t374 * t704;
t370 = t379 * pkin(5);
t563 = t378 * t379;
t501 = pkin(2) * t563;
t327 = qJD(2) * t501;
t329 = t382 * t511;
t517 = t327 + t329;
t613 = pkin(1) - t345;
t141 = (-t381 * t613 - t370) * qJD(1) - t517;
t309 = t381 * pkin(1) + t370;
t287 = t309 * qJD(1);
t553 = -t141 - t287;
t12 = t251 + t546 * t294 + t527 * t207 + t405 * t381 + t455 * qJDD(1) + (-t485 + t553 - t602) * qJD(1);
t450 = -t327 + t485;
t459 = t381 * t345 - t379 * t382;
t186 = t459 - t309;
t535 = t619 * t564 + t702 * t566 + t366;
t494 = -t186 - t535;
t59 = -t527 * t293 + (t309 - t494) * qJD(1) + t450;
t670 = (qJD(1) * t59 + t12) * t381;
t663 = -t619 * t567 + t521;
t662 = -t619 * t566 + t520;
t657 = t697 * t352 + t698 * t353;
t492 = t381 * t619;
t508 = qJD(2) * t381;
t509 = qJD(2) * t379;
t548 = -t185 * t509 + t186 * t508;
t53 = -qJD(4) * t353 + t293 * t536 + t294 * t535 + t548;
t616 = g(2) * t379;
t656 = -t53 * (qJD(4) * t352 + t293 * t663 + t294 * t662) - t59 * (qJD(1) * t662 - t293 * t693 + t353 * t506) - (-g(1) * t492 - t616 * t619) * t352;
t655 = t714 * qJD(1) - t706 * t293 + t705 * t294;
t654 = qJD(1) * t663 + t294 * t693 - t353 * t505 + t527 * t511;
t275 = rSges(4,1) * t352 + rSges(4,2) * t353;
t243 = t275 * t379;
t247 = t275 * t381;
t341 = t353 * rSges(4,1);
t647 = -rSges(4,2) * t352 + t341;
t201 = rSges(4,1) * t565 - rSges(4,2) * t567 - t381 * rSges(4,3);
t364 = t379 * rSges(4,3);
t203 = rSges(4,1) * t564 - rSges(4,2) * t566 + t364;
t68 = t201 * t293 + t203 * t294 + t548;
t487 = t378 * t508;
t457 = pkin(2) * t487;
t412 = -t275 * t294 - t457;
t495 = -t201 + t542;
t69 = qJD(1) * t495 + t412;
t541 = -t186 - t203;
t70 = -t275 * t293 - t327 + (t309 - t541) * qJD(1);
t653 = -t69 * (qJD(1) * t243 - t294 * t647) - t68 * (-t293 * t243 - t247 * t294) - t70 * (-qJD(1) * t247 - t293 * t647);
t288 = qJDD(2) * t379 + t381 * t503;
t206 = qJDD(3) * t379 + t381 * t502 + t288;
t350 = pkin(5) * t510;
t140 = -t457 - t350 + (t379 * t613 - t349) * qJD(1);
t530 = qJD(1) * (-pkin(1) * t511 + t350) + qJDD(1) * t309;
t456 = qJD(1) * t140 + qJDD(1) * t186 + t530;
t570 = t288 * t378;
t13 = -pkin(2) * t570 + t546 * t293 - t527 * t206 + t535 * qJDD(1) + t405 * t379 + (t306 + t603) * qJD(1) + t456;
t652 = t13 * t379;
t291 = qJD(1) * t308;
t649 = qJD(1) * t185 - t291;
t363 = Icges(3,4) * t380;
t437 = -Icges(3,2) * t378 + t363;
t299 = Icges(3,1) * t378 + t363;
t644 = g(1) * t381 + t616;
t561 = t379 * t380;
t585 = Icges(3,3) * t381;
t221 = Icges(3,5) * t561 - Icges(3,6) * t563 - t585;
t332 = Icges(3,4) * t563;
t595 = Icges(3,5) * t381;
t225 = Icges(3,1) * t561 - t332 - t595;
t589 = Icges(3,6) * t381;
t223 = Icges(3,4) * t561 - Icges(3,2) * t563 - t589;
t576 = t223 * t378;
t428 = -t225 * t380 + t576;
t79 = -t381 * t221 - t379 * t428;
t296 = Icges(3,5) * t380 - Icges(3,6) * t378;
t295 = Icges(3,5) * t378 + Icges(3,6) * t380;
t413 = qJD(2) * t295;
t598 = Icges(3,4) * t378;
t300 = Icges(3,1) * t380 - t598;
t226 = Icges(3,5) * t379 + t300 * t381;
t224 = Icges(3,6) * t379 + t381 * t437;
t575 = t224 * t378;
t427 = -t226 * t380 + t575;
t638 = -t381 * t413 + (-t296 * t379 + t427 + t585) * qJD(1);
t222 = Icges(3,3) * t379 + t296 * t381;
t513 = qJD(1) * t222;
t637 = qJD(1) * t428 - t379 * t413 + t513;
t297 = Icges(3,2) * t380 + t598;
t422 = t378 * t297 - t299 * t380;
t634 = t422 * qJD(1) + t296 * qJD(2);
t633 = t379 * (-t297 * t381 + t226) - t381 * (-Icges(3,2) * t561 + t225 - t332);
t629 = t206 / 0.2e1;
t628 = t207 / 0.2e1;
t627 = t288 / 0.2e1;
t626 = t289 / 0.2e1;
t625 = -t293 / 0.2e1;
t624 = t293 / 0.2e1;
t623 = -t294 / 0.2e1;
t622 = t294 / 0.2e1;
t621 = t379 / 0.2e1;
t620 = -t381 / 0.2e1;
t615 = -qJD(1) / 0.2e1;
t614 = qJD(1) / 0.2e1;
t612 = rSges(3,1) * t380;
t611 = rSges(3,2) * t380;
t365 = t379 * rSges(3,3);
t606 = t69 * t275;
t605 = qJDD(1) / 0.2e1;
t582 = qJD(1) * t53;
t301 = rSges(3,1) * t378 + t611;
t488 = t301 * t508;
t516 = rSges(3,2) * t563 + t381 * rSges(3,3);
t227 = rSges(3,1) * t561 - t516;
t532 = -t227 - t308;
t121 = qJD(1) * t532 - t488;
t580 = t121 * t379;
t579 = t121 * t381;
t560 = t380 * t381;
t562 = t378 * t381;
t228 = rSges(3,1) * t560 - rSges(3,2) * t562 + t365;
t165 = t228 + t309;
t122 = qJD(1) * t165 - t301 * t509;
t259 = t301 * t381;
t578 = t122 * t259;
t569 = t295 * t379;
t568 = t295 * t381;
t547 = -t379 * t185 + t381 * t186;
t545 = t379 * t201 + t381 * t203;
t544 = -t379 * t221 - t225 * t560;
t543 = t379 * t222 + t226 * t560;
t524 = rSges(4,2) * t491 + rSges(4,3) * t510;
t523 = -t297 + t300;
t522 = t299 + t437;
t518 = rSges(3,3) * t510 + t511 * t671;
t512 = qJD(1) * t296;
t125 = -t379 * t422 - t568;
t504 = t125 * qJD(1);
t499 = qJD(2) * t371;
t116 = rSges(4,1) * t410 - rSges(4,2) * t498 + t524;
t118 = -t374 * t243 + (t381 * t647 + t364) * qJD(1);
t497 = t381 * t116 + t379 * t118 + t201 * t510;
t496 = t381 * t140 + t379 * t141 - t185 * t510;
t489 = t378 * t510;
t484 = -pkin(1) - t612;
t483 = t511 / 0.2e1;
t482 = t510 / 0.2e1;
t481 = -t509 / 0.2e1;
t480 = t509 / 0.2e1;
t479 = -t508 / 0.2e1;
t478 = t508 / 0.2e1;
t411 = -t275 - t618;
t476 = t378 * (-t379 ^ 2 - t381 ^ 2);
t179 = t226 * t561;
t466 = t381 * t222 - t179;
t465 = -t189 + t577;
t460 = -t221 + t575;
t454 = t379 * t536 + t381 * t535;
t452 = -t527 - t618;
t220 = t647 * t374;
t451 = -t220 - t499;
t304 = rSges(2,1) * t381 - rSges(2,2) * t379;
t302 = rSges(2,1) * t379 + rSges(2,2) * t381;
t303 = t612 - t671;
t124 = t224 * t380 + t226 * t378;
t414 = qJD(2) * t297;
t133 = -t381 * t414 + (-t379 * t437 + t589) * qJD(1);
t415 = qJD(2) * t299;
t135 = -t381 * t415 + (-t300 * t379 + t595) * qJD(1);
t391 = -qJD(2) * t124 - t133 * t378 + t135 * t380 + t513;
t123 = t223 * t380 + t225 * t378;
t134 = qJD(1) * t224 - t379 * t414;
t136 = qJD(1) * t226 - t379 * t415;
t392 = qJD(1) * t221 - qJD(2) * t123 - t134 * t378 + t136 * t380;
t446 = -(t379 * t637 + t392 * t381) * t381 + (t379 * t638 + t391 * t381) * t379;
t445 = -(t392 * t379 - t381 * t637) * t381 + (t391 * t379 - t381 * t638) * t379;
t421 = t306 - t457;
t406 = -t294 * t527 + t421;
t58 = qJD(1) * t455 + t406;
t444 = t379 * t59 + t381 * t58;
t443 = -t379 * t70 - t381 * t69;
t80 = -t224 * t563 - t466;
t442 = t379 * t80 - t381 * t79;
t81 = -t223 * t562 - t544;
t82 = -t224 * t562 + t543;
t441 = t379 * t82 - t381 * t81;
t434 = -t122 * t379 - t579;
t137 = -t508 * t611 + (-t380 * t511 - t487) * rSges(3,1) + t518;
t258 = t301 * t379;
t138 = -qJD(2) * t258 + (t303 * t381 + t365) * qJD(1);
t433 = t137 * t381 + t138 * t379;
t426 = t227 * t379 + t228 * t381;
t423 = t297 * t380 + t299 * t378;
t420 = -t499 + t546;
t419 = t379 * t602 + t381 * t603 + t510 * t536;
t418 = t140 * t508 + t141 * t509 - t288 * t185 - t186 * t289;
t407 = t223 * t381 - t224 * t379;
t404 = (-t378 * t522 + t380 * t523) * qJD(1);
t402 = -t345 - t693;
t280 = t437 * qJD(2);
t281 = t300 * qJD(2);
t390 = qJD(1) * t295 - qJD(2) * t423 - t280 * t378 + t281 * t380;
t388 = -t378 * t633 + t407 * t380;
t387 = (t379 * t676 - t381 * t677) * t629 + (t379 * t678 - t381 * t679) * t628 + (t379 * t655 + t381 * t657) * t625 + (t689 * t381 + t688 * t379 + (t379 * t677 + t381 * t676) * qJD(1)) * t624 + (t687 * t381 + t686 * t379 + (t379 * t679 + t381 * t678) * qJD(1)) * t623 + (t379 * t657 - t381 * t655) * t622 + (t683 * qJD(1) + t674 * qJDD(1) + t676 * t206 + t677 * t207 + t688 * t293 + t294 * t689) * t621 + (qJD(1) * t682 + qJDD(1) * t675 + t206 * t678 + t207 * t679 + t293 * t686 + t294 * t687) * t620 + (t352 * t698 - t353 * t697) * t615 + (t681 * t381 + t680 * t379 + (t379 * t673 + t381 * t672) * qJD(1)) * t614 + (t379 * t672 - t381 * t673) * t605 + t685 * t483 + t684 * t482;
t284 = t303 * qJD(2);
t126 = -t381 * t422 + t569;
t120 = t126 * qJD(1);
t119 = t426 * qJD(2);
t65 = qJD(1) * t137 + qJDD(1) * t228 - t284 * t509 - t288 * t301 + t530;
t64 = -t284 * t508 + t289 * t301 + t532 * qJDD(1) + (-t138 - t287) * qJD(1);
t63 = t390 * t379 - t381 * t634;
t62 = t379 * t634 + t390 * t381;
t61 = -qJD(2) * t427 + t133 * t380 + t135 * t378;
t60 = -t428 * qJD(2) + t134 * t380 + t136 * t378;
t48 = qJD(2) * t441 + t120;
t47 = qJD(2) * t442 + t504;
t42 = qJD(1) * t116 + qJDD(1) * t203 - t206 * t275 - t220 * t293 + (-t379 * t559 - t570) * pkin(2) + t456;
t41 = -t381 * t500 + t207 * t275 - t220 * t294 + t251 + t495 * qJDD(1) + (-t118 + t553) * qJD(1);
t18 = t116 * t294 + t118 * t293 + t201 * t206 - t203 * t207 + t418;
t9 = -qJDD(4) * t353 + t206 * t536 - t207 * t535 + t293 * t602 + t294 * t603 + t352 * t507 + t418;
t1 = [(t120 + ((t80 - t179 + (t222 + t576) * t381 + t544) * t381 + t543 * t379) * qJD(2)) * t478 - m(2) * (-g(1) * t302 + g(2) * t304) + (t126 + t124) * t627 + (t125 + t123) * t626 + ((t74 + (t193 * t381 + t194 * t379) * t352 + t467 + t551) * t294 + (-t197 * t565 + t558 + t73 + (t193 * t379 - t194 * t381) * t352 + t550 + t76) * t293 + t692) * t622 + (t47 - t504 + ((t381 * t460 - t543 + t82) * t381 + (t379 * t460 + t466 + t81) * t379) * qJD(2)) * t481 + (t61 + t62) * t480 + (t60 + t63 + t48) * t479 + (-qJD(2) * t422 + t280 * t380 + t281 * t378 + t712 * t352 - t711 * t353) * qJD(1) + (-(-qJD(1) * t536 + t406 - t58 + t649) * t59 + t58 * (t329 - t450) + t59 * (t421 + t691) + (-t59 * t352 * t492 + (t352 * t619 - t353 * t702) * t379 * t58) * t374 + ((-t59 * t382 + t402 * t58) * t381 + (-t58 * rSges(5,2) + t402 * t59) * t379) * qJD(1) + (-g(2) + t13) * (t459 + t535) + (-g(1) + t12) * (t369 + (-t352 * t702 - t353 * t619) * t379 + t519)) * m(5) + (-(-qJD(1) * t201 + t412 + t649 - t69) * t70 + t69 * t517 + t70 * (-t457 + t524) + (-t247 * t70 + t379 * t606) * t374 + ((-t69 * rSges(4,3) + t70 * (-t345 - t341)) * t379 + (t69 * (-t345 - t647) - t70 * t382) * t381) * qJD(1) + (-g(2) + t42) * (t203 + t459) + (-g(1) + t41) * (-t201 + t519)) * m(4) + (-(-qJD(1) * t227 - t121 - t291 - t488) * t122 + t122 * (t350 + t518) + (t301 * t580 - t578) * qJD(2) + ((-pkin(1) - t303) * t579 + (t121 * (-rSges(3,3) - pkin(5)) + t122 * t484) * t379) * qJD(1) + (-g(2) + t65) * t165 + (-g(1) + t64) * (t484 * t379 + t372 + t516)) * m(3) + (t672 + t674) * t629 + (t673 + t675) * t628 + ((t465 * t381 - t416 - t550 + t78) * t294 + (t465 * t379 - t147 - t171 + t552 + (-t713 - t717) * t381 + t677) * t293 + t685 + t690) * t625 + (t680 + t683) * t624 + (m(2) * (t302 ^ 2 + t304 ^ 2) + t423 + Icges(2,3) + t725 * t353 + t726 * t352) * qJDD(1) + (-t681 + t682 + t684) * t623; (qJD(1) * t63 + qJD(2) * t445 + qJDD(1) * t125 + t288 * t80 + t289 * t79) * t620 + ((-t509 * t568 + t512) * t379 + (t404 + (t379 * t569 + t388) * qJD(2)) * t381) * t481 + t387 + (qJD(1) * t62 + qJD(2) * t446 + qJDD(1) * t126 + t288 * t82 + t289 * t81) * t621 + (t379 * t61 - t381 * t60 + (t123 * t379 + t124 * t381) * qJD(1)) * t614 + ((-t508 * t569 - t512) * t381 + (t404 + (t381 * t568 + t388) * qJD(2)) * t379) * t478 + ((t378 * t523 + t380 * t522) * qJD(1) + (t407 * t378 + t380 * t633) * qJD(2)) * t615 + t48 * t482 + t47 * t483 + ((t81 * t379 + t82 * t381) * qJD(1) + t446) * t480 + ((t79 * t379 + t80 * t381) * qJD(1) + t445) * t479 + (-t123 * t381 + t124 * t379) * t605 + t441 * t627 + t442 * t626 + (-g(1) * (-pkin(2) * t562 + t520) - g(2) * (-t501 + t521) - g(3) * (t371 + t693) - (-t59 * t489 + (-t380 * t444 + t476 * t53) * qJD(2)) * pkin(2) + t9 * (t454 + t547) + t53 * (t419 + t496) + t452 * t670 + (t13 * t452 + t420 * t59 + t494 * t582) * t379 + (t420 * t381 + t654) * t58 + t656) * m(5) + (-g(3) * (t647 + t371) - t644 * t411 - (-t70 * t489 + (t380 * t443 + t476 * t68) * qJD(2)) * pkin(2) + t18 * (t545 + t547) + t68 * (t496 + t497) + (t451 * t69 + (qJD(1) * t70 + t41) * t411) * t381 + (t42 * t411 + t70 * t451 + (t541 * t68 + t606) * qJD(1)) * t379 + t653) * m(4) + (-(t121 * t258 - t578) * qJD(1) - (t119 * (-t258 * t379 - t259 * t381) + t434 * t303) * qJD(2) + (qJD(2) * t433 + t227 * t288 - t228 * t289) * t426 + t119 * ((t227 * t381 - t228 * t379) * qJD(1) + t433) + t434 * t284 + (-t65 * t379 - t64 * t381 + (-t122 * t381 + t580) * qJD(1)) * t301 + g(1) * t259 + g(2) * t258 - g(3) * t303) * m(3); t387 + (-g(1) * t520 - g(2) * t521 - g(3) * t693 + t9 * t454 + t53 * t419 + (-t535 * t582 + t546 * t59) * t379 + (-t652 - t670) * t527 + (t546 * t381 + t654) * t58 + t656) * m(5) + (t18 * t545 + t68 * (-t203 * t511 + t497) + t443 * t220 + (-t42 * t379 - t41 * t381 + (t379 * t69 - t381 * t70) * qJD(1)) * t275 + g(1) * t247 + g(2) * t243 - g(3) * t647 + t653) * m(4); ((-t293 * t59 - t294 * t58 + t374 * t444 + g(3) - t9) * t353 + (t12 * t381 + t652 + (-t293 * t379 - t294 * t381 + t374) * t53 - t644) * t352) * m(5);];
tau = t1;
