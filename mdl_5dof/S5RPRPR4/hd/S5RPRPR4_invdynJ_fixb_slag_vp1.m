% Calculate vector of inverse dynamics joint torques for
% S5RPRPR4
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:09
% EndTime: 2019-12-05 17:53:56
% DurationCPUTime: 33.86s
% Computational Cost: add. (18025->793), mult. (14294->958), div. (0->0), fcn. (11105->10), ass. (0->447)
t734 = Icges(4,3) + Icges(5,3);
t356 = qJ(3) + pkin(9);
t346 = sin(t356);
t348 = cos(t356);
t252 = Icges(5,5) * t348 - Icges(5,6) * t346;
t359 = sin(qJ(3));
t361 = cos(qJ(3));
t304 = Icges(4,5) * t361 - Icges(4,6) * t359;
t727 = t252 + t304;
t357 = qJ(1) + pkin(8);
t347 = sin(t357);
t349 = cos(t357);
t733 = -t727 * t347 + t734 * t349;
t558 = t348 * t349;
t563 = t346 * t349;
t584 = Icges(5,6) * t347;
t173 = Icges(5,4) * t558 - Icges(5,2) * t563 + t584;
t287 = Icges(5,4) * t563;
t592 = Icges(5,5) * t347;
t175 = Icges(5,1) * t558 - t287 + t592;
t555 = t349 * t361;
t556 = t349 * t359;
t585 = Icges(4,6) * t347;
t191 = Icges(4,4) * t555 - Icges(4,2) * t556 + t585;
t318 = Icges(4,4) * t556;
t594 = Icges(4,5) * t347;
t193 = Icges(4,1) * t555 - t318 + t594;
t724 = t173 * t346 - t175 * t348 + t191 * t359 - t193 * t361;
t732 = t734 * t347;
t731 = Icges(4,5) * t555 + Icges(5,5) * t558 - Icges(4,6) * t556 - Icges(5,6) * t563 + t732;
t596 = Icges(5,4) * t346;
t253 = Icges(5,2) * t348 + t596;
t597 = Icges(4,4) * t359;
t305 = Icges(4,2) * t361 + t597;
t351 = Icges(4,4) * t361;
t651 = Icges(4,1) * t359 + t351;
t334 = Icges(5,4) * t348;
t653 = Icges(5,1) * t346 + t334;
t726 = t253 * t346 + t305 * t359 - t348 * t653 - t361 * t651;
t564 = t346 * t347;
t286 = Icges(5,4) * t564;
t561 = t347 * t348;
t591 = Icges(5,5) * t349;
t174 = -Icges(5,1) * t561 + t286 + t591;
t560 = t347 * t359;
t317 = Icges(4,4) * t560;
t559 = t347 * t361;
t593 = Icges(4,5) * t349;
t192 = -Icges(4,1) * t559 + t317 + t593;
t674 = -t174 * t558 - t192 * t555 - t733 * t347;
t423 = -Icges(5,2) * t346 + t334;
t172 = Icges(5,6) * t349 - t347 * t423;
t424 = -Icges(4,2) * t359 + t351;
t190 = Icges(4,6) * t349 - t347 * t424;
t673 = t172 * t564 + t190 * t560 + t733 * t349;
t730 = -t174 * t348 - t192 * t361;
t672 = t172 * t346 + t190 * t359;
t729 = t724 * t349;
t700 = -t174 * t561 - t192 * t559 + t673;
t684 = t173 * t564 - t175 * t561 + t191 * t560 - t193 * t559 + t731 * t349;
t699 = -t172 * t563 - t190 * t556 - t674;
t698 = t731 * t347 - t729;
t697 = t172 * t348 + t174 * t346 + t190 * t361 + t192 * t359;
t696 = t173 * t348 + t175 * t346 + t191 * t361 + t193 * t359;
t251 = Icges(5,5) * t346 + Icges(5,6) * t348;
t198 = t251 * t347;
t303 = Icges(4,5) * t359 + Icges(4,6) * t361;
t214 = t303 * t347;
t728 = t214 + t198;
t725 = t253 * t348 + t305 * t361 + t346 * t653 + t359 * t651;
t723 = -t672 - t730;
t569 = t303 * t349;
t571 = t251 * t349;
t722 = -t569 - t571;
t695 = t726 * t347 - t722;
t694 = -t726 * t349 + t728;
t721 = t733 * qJD(1);
t720 = t349 ^ 2;
t719 = t726 * qJD(1) + t727 * qJD(3);
t718 = t728 * qJD(3) + (-t727 * t349 - t723 - t732) * qJD(1);
t717 = t724 * qJD(1) + t722 * qJD(3) + t721;
t716 = t698 * t347 + t699 * t349;
t715 = t684 * t347 + t700 * t349;
t227 = t423 * qJD(3);
t256 = Icges(5,1) * t348 - t596;
t228 = t256 * qJD(3);
t279 = t424 * qJD(3);
t308 = Icges(4,1) * t361 - t597;
t280 = t308 * qJD(3);
t714 = t227 * t346 - t228 * t348 + t279 * t359 - t280 * t361 + t725 * qJD(3) + (-t251 - t303) * qJD(1);
t506 = qJD(3) * t349;
t103 = qJD(1) * t172 - t253 * t506;
t203 = t653 * t349;
t105 = -qJD(3) * t203 + (-t256 * t347 + t591) * qJD(1);
t119 = qJD(1) * t190 - t305 * t506;
t219 = t651 * t349;
t121 = -qJD(3) * t219 + (-t308 * t347 + t593) * qJD(1);
t713 = -qJD(1) * t731 + t696 * qJD(3) + t103 * t346 - t105 * t348 + t119 * t359 - t121 * t361;
t507 = qJD(3) * t347;
t104 = t253 * t507 + (-t349 * t423 - t584) * qJD(1);
t202 = t653 * t347;
t106 = qJD(3) * t202 + (-t256 * t349 - t592) * qJD(1);
t120 = t305 * t507 + (-t349 * t424 - t585) * qJD(1);
t218 = t651 * t347;
t122 = qJD(3) * t218 + (-t308 * t349 - t594) * qJD(1);
t712 = -t697 * qJD(3) - t104 * t346 + t106 * t348 - t120 * t359 + t122 * t361 + t721;
t711 = rSges(5,2) * t348;
t710 = t695 * qJD(1);
t709 = t694 * qJD(1);
t331 = qJDD(3) * t349;
t355 = qJD(3) + qJD(5);
t511 = qJD(1) * t347;
t168 = qJDD(5) * t349 - t355 * t511 + t331;
t350 = qJ(5) + t356;
t341 = cos(t350);
t327 = t341 * rSges(6,1);
t340 = sin(t350);
t603 = rSges(6,2) * t340;
t654 = t327 - t603;
t209 = t654 * t355;
t238 = rSges(6,1) * t340 + rSges(6,2) * t341;
t501 = qJD(1) * qJD(3);
t242 = -t347 * t501 + t331;
t250 = t349 * t355;
t362 = cos(qJ(1));
t364 = qJD(1) ^ 2;
t553 = t362 * t364;
t499 = pkin(1) * t553;
t360 = sin(qJ(1));
t500 = qJDD(1) * t360;
t614 = pkin(6) * t347;
t618 = pkin(2) * t349;
t262 = t614 + t618;
t613 = t349 * pkin(6);
t449 = -pkin(2) * t347 + t613;
t532 = qJDD(1) * t449 - t364 * t262;
t363 = qJD(3) ^ 2;
t557 = t348 * t363;
t616 = pkin(4) * t346;
t353 = t361 * pkin(3);
t342 = t353 + pkin(2);
t608 = pkin(2) - t342;
t476 = t349 * t608;
t358 = -qJ(4) - pkin(6);
t508 = qJD(1) * t358;
t309 = t347 * t508;
t505 = qJD(3) * t359;
t480 = t347 * t505;
t298 = pkin(3) * t480;
t333 = qJD(4) * t349;
t524 = t298 + t333;
t484 = t309 + t524;
t113 = (t476 + t614) * qJD(1) + t484;
t607 = pkin(6) + t358;
t663 = t349 * t607;
t667 = t347 * t608;
t395 = -t663 + t667;
t554 = t361 * t363;
t632 = qJDD(1) * t395 + qJDD(4) * t347 + (-t242 * t359 - t349 * t554) * pkin(3) + (t113 + t333) * qJD(1);
t466 = pkin(4) * t348 + t353;
t453 = pkin(2) + t466;
t428 = t342 - t453;
t354 = -pkin(7) + t358;
t519 = t358 - t354;
t371 = t347 * t428 + t349 * t519;
t337 = t349 * rSges(6,3);
t568 = t340 * t347;
t527 = rSges(6,2) * t568 + t337;
t566 = t341 * t347;
t408 = rSges(6,1) * t566 - t527;
t690 = t371 - t408;
t617 = pkin(3) * t359;
t448 = t616 + t617;
t257 = t448 * qJD(3);
t677 = qJD(1) * t354 + t257;
t455 = t677 * t347;
t510 = qJD(1) * t349;
t74 = t428 * t510 - t298 - t309 + t455;
t565 = t341 * t349;
t494 = rSges(6,1) * t565;
t403 = -t347 * rSges(6,3) - t494;
t194 = rSges(6,1) * t568 + rSges(6,2) * t566;
t567 = t340 * t349;
t276 = rSges(6,2) * t567;
t489 = qJD(1) * t276 + t194 * t355;
t93 = qJD(1) * t403 + t489;
t9 = -pkin(4) * t349 * t557 - pkin(1) * t500 - t168 * t238 - t250 * t209 - t242 * t616 - t499 + t532 + t690 * qJDD(1) + (t74 + t93) * qJD(1) + t632;
t708 = t9 - g(3);
t497 = rSges(5,1) * t558;
t404 = rSges(5,3) * t347 + t497;
t291 = rSges(5,2) * t563;
t482 = t346 * t507;
t488 = rSges(5,1) * t482 + qJD(1) * t291 + t507 * t711;
t108 = -qJD(1) * t404 + t488;
t498 = rSges(5,1) * t561;
t338 = t349 * rSges(5,3);
t525 = rSges(5,2) * t564 + t338;
t176 = -t498 + t525;
t336 = t348 * rSges(5,1);
t652 = -rSges(5,2) * t346 + t336;
t230 = t652 * qJD(3);
t259 = rSges(5,1) * t346 + t711;
t373 = (-t500 - t553) * pkin(1) + t532;
t29 = qJD(1) * t108 + qJDD(1) * t176 - t230 * t506 - t242 * t259 + t373 + t632;
t707 = t29 - g(3);
t706 = qJD(3) * t716 + t709;
t705 = qJD(3) * t715 + t710;
t704 = t723 * qJD(3) + t104 * t348 + t106 * t346 + t120 * t361 + t122 * t359;
t703 = -t724 * qJD(3) + t103 * t348 + t105 * t346 + t119 * t361 + t121 * t359;
t702 = t347 * t719 - t349 * t714;
t701 = t347 * t714 + t349 * t719;
t249 = t347 * t355;
t693 = t209 * t347 - t249 * t654;
t691 = -t371 - t395;
t689 = t718 * t720 + (t713 * t347 + (-t712 + t717) * t349) * t347;
t688 = t712 * t720 + (t717 * t347 + (-t713 + t718) * t349) * t347;
t687 = qJD(1) * t194 + t209 * t349 - t250 * t654;
t269 = Icges(6,4) * t568;
t589 = Icges(6,5) * t349;
t160 = -Icges(6,1) * t566 + t269 + t589;
t270 = Icges(6,4) * t567;
t590 = Icges(6,5) * t347;
t161 = Icges(6,1) * t565 - t270 + t590;
t325 = Icges(6,4) * t341;
t422 = -Icges(6,2) * t340 + t325;
t655 = Icges(6,1) * t340 + t325;
t670 = t655 + t422;
t369 = qJD(1) * t670 + t249 * (-Icges(6,2) * t565 + t161 - t270) + t250 * (Icges(6,2) * t566 + t160 + t269);
t595 = Icges(6,4) * t340;
t233 = Icges(6,2) * t341 + t595;
t236 = Icges(6,1) * t341 - t595;
t678 = t233 - t236;
t583 = Icges(6,6) * t347;
t159 = Icges(6,4) * t565 - Icges(6,2) * t567 + t583;
t679 = t655 * t349 + t159;
t158 = Icges(6,6) * t349 - t347 * t422;
t680 = -t655 * t347 + t158;
t631 = qJD(1) * t678 + t249 * t679 + t250 * t680;
t686 = t369 * t340 + t341 * t631;
t384 = t347 * (-Icges(5,2) * t558 + t175 - t287) + t349 * (Icges(5,2) * t561 + t174 + t286);
t639 = t347 * (t173 + t203) + t349 * (t172 - t202);
t685 = -t384 * t346 - t348 * t639;
t683 = -t113 - t74;
t682 = t347 * t652;
t580 = Icges(6,3) * t347;
t157 = Icges(6,5) * t565 - Icges(6,6) * t567 + t580;
t681 = t349 * t157 - t161 * t566;
t676 = -t731 + t730;
t612 = t360 * pkin(1);
t671 = t612 - t525;
t521 = t651 + t424;
t522 = t305 - t308;
t669 = (t359 * t521 + t361 * t522) * qJD(1);
t528 = t653 + t423;
t529 = t253 - t256;
t668 = (t346 * t528 + t348 * t529) * qJD(1);
t396 = t259 + t617;
t660 = t396 * t506;
t225 = t349 * t453;
t293 = t349 * t342;
t129 = -t347 * t519 - t225 + t293;
t162 = -t276 - t403;
t659 = t129 - t162;
t155 = t347 * t607 - t293 + t618;
t248 = qJD(1) * t262;
t658 = -qJD(1) * t155 + t248;
t657 = t652 + t353;
t611 = t362 * pkin(1);
t656 = -t293 - t611;
t606 = rSges(3,1) * t349;
t446 = -t606 - t611;
t224 = t347 * rSges(3,2) + t446;
t650 = t691 * t347;
t649 = t690 * qJD(1);
t232 = Icges(6,5) * t341 - Icges(6,6) * t340;
t156 = Icges(6,3) * t349 - t232 * t347;
t549 = -t347 * t156 - t160 * t565;
t648 = t549 - t681;
t647 = (-t347 ^ 2 - t720) * t617;
t533 = -Icges(4,2) * t555 + t193 - t318;
t535 = t191 + t219;
t637 = t359 * t533 + t361 * t535;
t534 = Icges(4,2) * t559 + t192 + t317;
t536 = t190 - t218;
t636 = -t359 * t534 - t361 * t536;
t231 = Icges(6,5) * t340 + Icges(6,6) * t341;
t456 = t678 * t355;
t457 = t670 * t355;
t635 = -qJD(1) * t231 + t340 * t457 + t341 * t456;
t461 = qJD(1) * t158 + t161 * t355 - t233 * t250;
t463 = -(-t236 * t347 + t589) * qJD(1) + t679 * t355;
t634 = -qJD(1) * t157 + t340 * t461 + t341 * t463;
t462 = t160 * t355 + t233 * t249 + (-t349 * t422 - t583) * qJD(1);
t464 = -(-t236 * t349 - t590) * qJD(1) + t680 * t355;
t518 = qJD(1) * t156;
t633 = t340 * t462 + t341 * t464 - t518;
t630 = -m(5) - m(6);
t241 = qJDD(3) * t347 + t349 * t501;
t167 = qJD(5) * t510 + qJDD(5) * t347 + t241;
t629 = t167 / 0.2e1;
t628 = t168 / 0.2e1;
t627 = t241 / 0.2e1;
t626 = t242 / 0.2e1;
t625 = -t249 / 0.2e1;
t624 = t249 / 0.2e1;
t623 = -t250 / 0.2e1;
t622 = t250 / 0.2e1;
t621 = t347 / 0.2e1;
t620 = t349 / 0.2e1;
t619 = -rSges(4,3) - pkin(6);
t610 = -qJD(1) / 0.2e1;
t609 = qJD(1) / 0.2e1;
t605 = rSges(4,1) * t361;
t339 = t349 * rSges(4,3);
t602 = qJDD(1) / 0.2e1;
t601 = t354 - rSges(6,3);
t182 = t231 * t347;
t574 = t231 * t349;
t573 = t238 * t349;
t570 = t448 * t347;
t562 = t347 * t157;
t414 = t233 * t340 - t341 * t655;
t79 = -t349 * t414 + t182;
t552 = t79 * qJD(1);
t550 = t349 * t156 + t158 * t568;
t526 = t342 * t511 + t349 * t508;
t204 = rSges(5,1) * t564 + rSges(5,2) * t561;
t520 = rSges(4,2) * t560 + t339;
t513 = qJD(1) * t252;
t512 = qJD(1) * t304;
t504 = qJD(3) * t361;
t332 = qJD(4) * t347;
t502 = -t155 * t506 + qJD(2);
t321 = pkin(3) * t560;
t496 = pkin(3) * t505;
t495 = qJD(1) * t612;
t493 = t158 * t567;
t490 = t259 * t507 + t524;
t486 = -t332 + t526;
t485 = rSges(4,1) * t480 + rSges(4,2) * (t347 * t504 + t359 * t510);
t326 = pkin(2) * t511;
t112 = t326 + (-pkin(6) * qJD(1) - t496) * t349 - t486;
t478 = t112 * t506 - t242 * t155 + qJDD(2);
t477 = -pkin(2) - t605;
t475 = -t511 / 0.2e1;
t474 = t510 / 0.2e1;
t473 = -t507 / 0.2e1;
t472 = t507 / 0.2e1;
t471 = -t506 / 0.2e1;
t470 = t506 / 0.2e1;
t468 = -t262 - t611;
t467 = t238 + t616;
t92 = -t355 * t573 + (-t347 * t654 + t337) * qJD(1);
t465 = -t347 * t93 + t349 * t92;
t343 = t364 * t612;
t454 = t347 * pkin(3) * t554 + qJDD(4) * t349 + t241 * t617 + t343;
t452 = t155 + t468;
t320 = rSges(4,2) * t556;
t405 = -rSges(4,1) * t555 - rSges(4,3) * t347;
t197 = -t320 - t405;
t451 = -t197 + t468;
t450 = -qJD(1) * t449 + t495;
t237 = pkin(6) * t510 - t326;
t445 = -t112 - t237 - t332;
t444 = -t176 + t663;
t314 = rSges(2,1) * t362 - t360 * rSges(2,2);
t443 = rSges(2,1) * t360 + rSges(2,2) * t362;
t442 = rSges(3,1) * t347 + rSges(3,2) * t349;
t313 = -rSges(4,2) * t359 + t605;
t311 = rSges(4,1) * t359 + rSges(4,2) * t361;
t196 = -rSges(4,1) * t559 + t520;
t95 = -qJD(1) * t196 + t311 * t506 + t450;
t240 = t311 * t507;
t96 = qJD(1) * t451 + t240;
t431 = t347 * t96 + t349 * t95;
t430 = qJD(1) * t467;
t35 = qJD(3) * t650 - t129 * t506 + t250 * t162 + t249 * t408 + t502;
t429 = t35 * (-t194 * t249 - t250 * t573);
t221 = t311 * t349;
t125 = -qJD(3) * t221 + (-t313 * t347 + t339) * qJD(1);
t126 = qJD(1) * t405 + t485;
t421 = t125 * t349 - t126 * t347;
t77 = t159 * t341 + t161 * t340;
t420 = t159 * t340 - t161 * t341;
t415 = -t196 * t347 + t197 * t349;
t177 = -t291 + t404;
t409 = -t177 + t452;
t406 = pkin(4) * t482 + t249 * t238 + t524;
t402 = t452 + t659;
t205 = t259 * t349;
t401 = -t453 - t327;
t149 = qJD(1) * t395;
t400 = -t149 - t332 + t450;
t107 = -qJD(3) * t205 + (t338 - t682) * qJD(1);
t399 = t349 * t107 + (-t108 - t113) * t347;
t397 = t442 + t612;
t390 = qJD(1) * t232 + t182 * t250 - t249 * t574;
t387 = qJD(1) * t420 - t355 * t574 + t518;
t386 = t355 * t182 + (t158 * t340 - t160 * t341 - t232 * t349 - t580) * qJD(1);
t379 = qJD(1) * t414 + t232 * t355;
t376 = t349 * t162 + t347 * t408;
t375 = t444 - t667;
t11 = t386 * t347 - t349 * t633;
t12 = t387 * t347 - t349 * t634;
t13 = t347 * t633 + t386 * t349;
t14 = t347 * t634 + t387 * t349;
t59 = -t160 * t566 + t550;
t131 = t159 * t568;
t60 = t131 + t681;
t78 = t347 * t414 + t574;
t75 = t78 * qJD(1);
t23 = t249 * t60 + t250 * t59 + t75;
t61 = -t493 - t549;
t62 = -t349 * t420 + t562;
t24 = t249 * t62 + t250 * t61 + t552;
t38 = t379 * t347 - t349 * t635;
t39 = t347 * t635 + t379 * t349;
t40 = -t340 * t464 + t341 * t462;
t41 = -t340 * t463 + t341 * t461;
t76 = t158 * t341 + t160 * t340;
t372 = (qJD(1) * t38 + qJDD(1) * t79 + t11 * t250 + t12 * t249 + t167 * t62 + t168 * t61) * t621 + (-t340 * t631 + t341 * t369) * t610 + (qJD(1) * t39 + qJDD(1) * t78 + t13 * t250 + t14 * t249 + t167 * t60 + t168 * t59) * t620 + t23 * t475 + t24 * t474 + (t11 * t349 + t12 * t347 + (-t347 * t61 + t349 * t62) * qJD(1)) * t624 + (t347 * t62 + t349 * t61) * t629 + (t347 * t60 + t349 * t59) * t628 + (t13 * t349 + t14 * t347 + (-t347 * t59 + t349 * t60) * qJD(1)) * t622 + (t347 * t77 + t349 * t76) * t602 + (t347 * t41 + t349 * t40 + (-t347 * t76 + t349 * t77) * qJD(1)) * t609 + (t390 * t347 - t686 * t349) * t625 + (t347 * t686 + t390 * t349) * t623;
t366 = t349 * t177 + t347 * t375;
t365 = t159 * t567 - t562 + (-t160 * t347 - t161 * t349) * t341 + t550;
t324 = rSges(3,2) * t511;
t281 = t313 * qJD(3);
t229 = t349 * t448;
t220 = t311 * t347;
t210 = t442 * qJD(1) + t495;
t181 = pkin(3) * t556 - t229;
t180 = -t321 + t570;
t141 = t349 * t155;
t100 = t349 * t112;
t94 = qJD(3) * t415 + qJD(2);
t73 = -t453 * t511 + (t496 - t677) * t349 + t526;
t64 = qJD(1) * t409 + t490;
t63 = -qJD(1) * t176 + t400 + t660;
t55 = qJD(3) * t366 + t502;
t54 = t281 * t507 + t241 * t311 + t343 + (-t125 - t237) * qJD(1) + t451 * qJDD(1);
t53 = qJD(1) * t126 + qJDD(1) * t196 - t242 * t311 - t281 * t506 + t373;
t49 = qJD(1) * t402 + t406;
t48 = t238 * t250 + t448 * t506 + t400 - t649;
t47 = qJD(3) * t421 - t196 * t241 + t197 * t242 + qJDD(2);
t30 = t230 * t507 + t241 * t259 + t409 * qJDD(1) + (-t107 + t445) * qJD(1) + t454;
t15 = qJD(3) * t399 + t242 * t177 + t241 * t375 + t478;
t10 = t167 * t238 + t209 * t249 + (t241 * t346 + t347 * t557) * pkin(4) + t402 * qJDD(1) + (t445 - t73 - t92) * qJD(1) + t454;
t5 = -t242 * t129 + t168 * t162 + t167 * t408 + t241 * t691 - t249 * t93 + t250 * t92 + t506 * t73 + t683 * t507 + t478;
t1 = [(t75 + (t62 + t365) * t250 + (t131 - t61 - t493 - t648) * t249) * t625 + (t77 + t79) * t629 + (t76 + t78) * t628 + (t24 - t552 + (t60 + (t158 * t349 - t159 * t347) * t340 + t648) * t250 + (-t59 + t365) * t249) * t623 + (t40 + t39) * t622 + (((t673 + t698 + t729) * t349 + ((-t672 + t676) * t349 - t674 + t684 - t699) * t347) * qJD(3) + t710) * t473 + (-t726 * qJD(3) + t227 * t348 + t228 * t346 + t279 * t361 + t280 * t359 - t340 * t456 + t341 * t457) * qJD(1) + (-t210 * t324 + (qJD(1) * t224 * t397 - t210 * t446) * qJD(1) + (t364 * t442 - g(2) + t343) * t224 + (t499 + (t606 * qJD(1) - t324) * qJD(1) + g(3)) * t397) * m(3) + (g(2) * t314 + g(3) * t443) * m(2) + (t41 + t38 + t23) * t624 + (-t49 * t332 - t48 * (t333 + t455 + t489) + t49 * (t238 * t355 + t257) * t349 + ((t360 * t49 + t362 * t48) * pkin(1) + (-t401 * t48 + t49 * t601) * t349 + (t49 * (-t401 - t603) + t48 * rSges(6,3)) * t347) * qJD(1) - (t49 + (t611 - t659) * qJD(1) - t406 + t658) * t48 + (t10 - g(2)) * (t347 * t601 - t225 + t276 - t494 - t611) + t708 * (t347 * t401 - t349 * t354 + t527 - t612)) * m(6) + ((t486 + t660 + (t498 + t671) * qJD(1)) * t64 + (-t484 - t488 - t64 + t490 - t658 + (t404 - t656 - t177 - t611) * qJD(1)) * t63 + (t30 - g(2)) * (-t497 + t291 + (-rSges(5,3) + t358) * t347 + t656) + t707 * (-t349 * t358 + (-t342 - t336) * t347 - t671)) * m(5) + (t96 * (t326 + (rSges(4,1) * t505 + rSges(4,2) * t504) * t349) - t95 * t485 + ((t360 * t96 + t362 * t95) * pkin(1) + (-t477 * t95 + t619 * t96) * t349 + (t313 * t96 - t619 * t95) * t347) * qJD(1) - (-t240 + t248 + t96 + (t197 + t611) * qJD(1)) * t95 + (t54 - g(2)) * (t619 * t347 + t477 * t349 + t320 - t611) + (t53 - g(3)) * (t347 * t477 + t520 - t612 + t613)) * m(4) + (t694 + t696) * t627 + (t695 + t697) * t626 + ((((t672 - t731) * t349 + t674 + t684) * t349 + (t347 * t676 + t673 - t700) * t347) * qJD(3) + t706 - t709) * t471 + (t701 + t704) * t470 + (m(2) * (t314 ^ 2 + t443 ^ 2) + m(3) * (t224 ^ 2 + t397 ^ 2) + t233 * t341 + t655 * t340 + Icges(2,3) + Icges(3,3) + t725) * qJDD(1) + (t702 + t703 + t705) * t472; m(3) * qJDD(2) + m(4) * t47 + m(5) * t15 + m(6) * t5 + (-m(3) - m(4) + t630) * g(1); t372 + t716 * t627 + t715 * t626 + (t702 * qJD(1) + t688 * qJD(3) + t694 * qJDD(1) + t698 * t241 + t699 * t242) * t621 + (t701 * qJD(1) + t689 * qJD(3) + t695 * qJDD(1) + t684 * t241 + t700 * t242) * t620 + ((-t346 * t639 + t348 * t384 + (t347 * t533 + t349 * t534) * t361 + (-t347 * t535 - t349 * t536) * t359) * qJD(3) + (-t346 * t529 + t348 * t528 - t359 * t522 + t361 * t521) * qJD(1)) * t610 + (t704 * t349 + t703 * t347 + (-t347 * t697 + t696 * t349) * qJD(1)) * t609 + (t696 * t347 + t697 * t349) * t602 + t705 * t475 + t706 * t474 + ((-t507 * t571 + t513) * t347 + (-t668 + (t347 * t198 + t685) * qJD(3)) * t349 + (-t507 * t569 + t512) * t347 + (-t669 + (t636 * t349 + (t214 - t637) * t347) * qJD(3)) * t349) * t473 + ((-t699 * t347 + t698 * t349) * qJD(1) + t688) * t472 + ((t198 * t506 + t513) * t349 + (t668 + (-t349 * t571 - t685) * qJD(3)) * t347 + (t214 * t506 + t512) * t349 + (t669 + (t637 * t347 + (-t569 - t636) * t349) * qJD(3)) * t347) * t471 + ((-t700 * t347 + t684 * t349) * qJD(1) + t689) * t470 + (-g(1) * (t654 + t466) - g(2) * (t194 + t570) - g(3) * (-t229 - t573) - t429 + t5 * (-t141 + t376 + t650) + t10 * (t347 * t467 + t321) + (-t5 * t129 + t9 * (-t238 - t448)) * t349 + (-(-t181 + t573) * qJD(1) + t430 * t349 + t693) * t49 + (qJD(1) * t180 - t347 * t430 + t687) * t48 + (-(-t180 * t347 + t181 * t349 + t647) * qJD(3) + t100 + t465 + (t155 + t659) * t511 + (-t149 + t73 - t649) * t349 + t683 * t347) * t35) * m(6) + (t15 * (-t141 + t366) + t55 * (t100 + (t444 * t349 + (t155 - t476 - t177) * t347) * qJD(1) + t399) + t30 * (t259 * t347 + t321) - t63 * (t259 * t511 + (-pkin(3) * t504 - t230) * t349) - g(1) * t657 - g(2) * (t321 + t204) + t63 * qJD(1) * t204 - (t55 * (-t204 * t347 + t647) + (-t55 * t205 + t63 * t657) * t349) * qJD(3) - t707 * t349 * t396 + (-qJD(1) * t205 - t682 * qJD(3) + t230 * t347 + t259 * t510) * t64) * m(5) + (-(-t220 * t95 + t221 * t96) * qJD(1) - (t94 * (-t220 * t347 - t221 * t349) + t431 * t313) * qJD(3) - g(1) * t313 - g(2) * t220 + g(3) * t221 + t47 * t415 + t94 * ((-t196 * t349 - t197 * t347) * qJD(1) + t421) + t431 * t281 + (t54 * t347 - t53 * t349 + (-t347 * t95 + t349 * t96) * qJD(1)) * t311) * m(4); t630 * (g(2) * t349 + g(3) * t347) + m(5) * (t29 * t347 + t30 * t349) + m(6) * (t10 * t349 + t347 * t9); t372 + (t5 * t376 + t35 * ((-t347 * t162 + t349 * t408) * qJD(1) + t465) - t429 - g(1) * t654 - g(2) * t194 + t693 * t49 + t687 * t48 + (t10 * t347 - t48 * t511 + t49 * t510) * t238 + (-t49 * qJD(1) - t708) * t573) * m(6);];
tau = t1;
