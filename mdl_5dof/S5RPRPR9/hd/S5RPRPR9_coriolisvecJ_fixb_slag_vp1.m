% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR9_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:24:31
% DurationCPUTime: 41.81s
% Computational Cost: add. (17131->863), mult. (22683->1133), div. (0->0), fcn. (20500->8), ass. (0->426)
t716 = Icges(5,1) + Icges(4,3);
t351 = sin(qJ(3));
t354 = cos(qJ(3));
t290 = Icges(4,5) * t354 - Icges(4,6) * t351;
t432 = Icges(5,4) * t354 - Icges(5,5) * t351;
t715 = t290 - t432;
t349 = qJ(1) + pkin(8);
t345 = cos(t349);
t714 = t716 * t345;
t344 = sin(t349);
t572 = t344 * t354;
t573 = t344 * t351;
t709 = t714 + (-Icges(5,5) + Icges(4,6)) * t573 + (Icges(5,4) - Icges(4,5)) * t572;
t696 = t716 * t344 + t715 * t345;
t598 = Icges(4,4) * t351;
t296 = Icges(4,1) * t354 - t598;
t587 = Icges(5,6) * t351;
t428 = Icges(5,2) * t354 - t587;
t713 = t296 + t428;
t586 = Icges(5,6) * t354;
t426 = -Icges(5,3) * t351 + t586;
t347 = Icges(4,4) * t354;
t433 = -Icges(4,2) * t351 + t347;
t712 = t426 + t433;
t588 = Icges(4,6) * t345;
t178 = Icges(4,4) * t572 - Icges(4,2) * t573 - t588;
t311 = Icges(5,6) * t573;
t596 = Icges(5,4) * t345;
t185 = Icges(5,2) * t572 - t311 + t596;
t711 = t178 * t351 - t185 * t354;
t181 = Icges(4,5) * t344 + t296 * t345;
t182 = Icges(5,5) * t344 - t345 * t426;
t710 = -t181 * t572 - t182 * t573;
t293 = Icges(4,2) * t354 + t598;
t295 = Icges(4,1) * t351 + t347;
t413 = t293 * t351 - t295 * t354;
t425 = Icges(5,3) * t354 + t587;
t427 = Icges(5,2) * t351 + t586;
t703 = -t351 * t425 + t354 * t427 - t413;
t316 = Icges(4,4) * t573;
t592 = Icges(4,5) * t345;
t180 = Icges(4,1) * t572 - t316 - t592;
t591 = Icges(5,5) * t345;
t183 = Icges(5,6) * t572 - Icges(5,3) * t573 + t591;
t691 = -t180 * t354 + t183 * t351 + t711;
t708 = t712 * qJD(3);
t707 = t713 * qJD(3);
t705 = -t293 - t425;
t704 = t295 + t427;
t289 = Icges(4,5) * t351 + Icges(4,6) * t354;
t431 = Icges(5,4) * t351 + Icges(5,5) * t354;
t702 = -t431 + t289;
t701 = t696 * t345 + t710;
t570 = t345 * t354;
t571 = t345 * t351;
t671 = t181 * t570 + t182 * t571 + t696 * t344;
t700 = -t180 * t570 + t183 * t571 + t709 * t344;
t699 = -t344 * t691 + t345 * t709;
t179 = Icges(4,6) * t344 + t345 * t433;
t312 = Icges(5,6) * t571;
t597 = Icges(5,4) * t344;
t184 = -Icges(5,2) * t570 + t312 + t597;
t678 = -t179 * t573 - t184 * t572 - t701;
t677 = -t178 * t571 + t185 * t570 - t700;
t676 = -t179 * t571 - t184 * t570 + t671;
t229 = t431 * t344;
t575 = t289 * t344;
t698 = t345 * t703 - t229 + t575;
t697 = t179 * t351 + t184 * t354;
t695 = t179 - t182;
t694 = t707 * t354 - t708 * t351 + (-t351 * t704 + t354 * t705) * qJD(3) + t702 * qJD(1);
t693 = t705 * qJD(3);
t692 = t704 * qJD(3);
t690 = t181 * t354 + t182 * t351 - t697;
t689 = t703 * qJD(1) - qJD(3) * t715;
t688 = t698 * qJD(1);
t687 = t702 * qJD(3);
t686 = (t678 * t344 - t345 * t699) * qJD(3);
t685 = (t344 * t676 - t345 * t677) * qJD(3);
t352 = sin(qJ(1));
t627 = pkin(1) * t352;
t355 = cos(qJ(1));
t348 = t355 * pkin(1);
t230 = t431 * t345;
t513 = (t425 * t573 - t427 * t572 - t230) * qJD(1);
t574 = t289 * t345;
t133 = -t344 * t413 - t574;
t514 = t133 * qJD(1);
t684 = t514 - t513 + t686;
t683 = t685 + t688;
t682 = (-t695 * qJD(1) - t693 * t344) * t354 + (t692 * t344 + (-t345 * t428 - t181 + t597) * qJD(1)) * t351 + t691 * qJD(3);
t681 = (t693 * t345 + (-t712 * t344 + t588 - t591) * qJD(1)) * t354 + (-t692 * t345 + (-t713 * t344 + t592 - t596) * qJD(1)) * t351 + t690 * qJD(3);
t680 = -t689 * t344 + t694 * t345;
t679 = t694 * t344 + t689 * t345;
t675 = (t178 + t183) * t354 + (t180 + t185) * t351;
t674 = t695 * t354 + (t181 - t184) * t351;
t673 = (-t687 * t344 + (t691 + t696) * qJD(1)) * t345;
t672 = t697 + t709;
t524 = qJD(1) * t345;
t334 = t344 * rSges(5,1);
t197 = -rSges(5,2) * t570 + rSges(5,3) * t571 + t334;
t568 = t351 * qJ(4);
t242 = pkin(3) * t570 + t345 * t568;
t263 = t345 * pkin(2) + t344 * pkin(6);
t653 = t348 + t263;
t660 = t242 + t653;
t670 = t197 + t660;
t667 = -t687 * t345 + (-t344 * t715 - t690 + t714) * qJD(1);
t515 = qJD(5) * t354;
t521 = qJD(3) * t344;
t248 = t345 * t515 + t521;
t520 = qJD(3) * t345;
t249 = -t344 * t515 + t520;
t516 = qJD(5) * t351;
t329 = qJD(1) + t516;
t350 = sin(qJ(5));
t353 = cos(qJ(5));
t567 = t351 * t353;
t213 = -t344 * t350 + t345 * t567;
t569 = t350 * t351;
t214 = t344 * t353 + t345 * t569;
t107 = Icges(6,5) * t214 + Icges(6,6) * t213 + Icges(6,3) * t570;
t595 = Icges(6,4) * t214;
t110 = Icges(6,2) * t213 + Icges(6,6) * t570 + t595;
t203 = Icges(6,4) * t213;
t113 = Icges(6,1) * t214 + Icges(6,5) * t570 + t203;
t35 = t107 * t570 + t213 * t110 + t214 * t113;
t215 = t344 * t567 + t345 * t350;
t216 = -t344 * t569 + t345 * t353;
t109 = -Icges(6,5) * t216 + Icges(6,6) * t215 + Icges(6,3) * t572;
t205 = Icges(6,4) * t216;
t112 = Icges(6,2) * t215 + Icges(6,6) * t572 - t205;
t204 = Icges(6,4) * t215;
t114 = Icges(6,1) * t216 - Icges(6,5) * t572 - t204;
t36 = t109 * t570 + t213 * t112 - t114 * t214;
t429 = Icges(6,5) * t350 + Icges(6,6) * t353;
t377 = -Icges(6,3) * t351 + t354 * t429;
t594 = Icges(6,4) * t350;
t430 = Icges(6,2) * t353 + t594;
t378 = -Icges(6,6) * t351 + t354 * t430;
t593 = Icges(6,4) * t353;
t434 = Icges(6,1) * t350 + t593;
t379 = -Icges(6,5) * t351 + t354 * t434;
t66 = -t213 * t378 - t214 * t379 - t377 * t570;
t10 = t248 * t35 - t249 * t36 + t66 * t329;
t37 = t107 * t572 + t215 * t110 - t216 * t113;
t38 = t109 * t572 + t112 * t215 + t114 * t216;
t67 = -t215 * t378 + t216 * t379 - t377 * t572;
t11 = t248 * t37 - t249 * t38 + t329 * t67;
t423 = t112 * t353 - t114 * t350;
t43 = t109 * t351 - t354 * t423;
t357 = qJD(1) ^ 2;
t665 = 0.2e1 * qJD(3);
t333 = t344 * rSges(4,3);
t196 = rSges(4,1) * t570 - rSges(4,2) * t571 + t333;
t661 = t196 + t653;
t519 = qJD(3) * t351;
t495 = t345 * t519;
t522 = qJD(1) * t354;
t501 = t344 * t522;
t659 = t495 + t501;
t658 = t10 * t345 + t11 * t344;
t525 = qJD(1) * t344;
t446 = rSges(6,1) * t216 - rSges(6,2) * t215;
t120 = rSges(6,3) * t572 - t446;
t346 = qJD(4) * t351;
t299 = t345 * t346;
t445 = rSges(6,1) * t350 + rSges(6,2) * t353;
t382 = -rSges(6,3) * t351 + t354 * t445;
t301 = pkin(3) * t351 - qJ(4) * t354;
t473 = -pkin(7) * t351 - t301;
t451 = t473 * t345;
t657 = qJD(3) * t451 - t120 * t329 + t249 * t382 + t299;
t604 = t351 * rSges(5,2);
t444 = rSges(5,3) * t354 + t604;
t655 = (-qJD(3) * t444 - t346) * t344;
t305 = pkin(3) * t354 + t568;
t238 = t305 * t344;
t208 = qJD(1) * t238;
t339 = t345 * pkin(6);
t262 = pkin(2) * t344 - t339;
t256 = qJD(1) * t262;
t654 = -t208 - t256;
t469 = t345 * rSges(3,1) - rSges(3,2) * t344;
t323 = pkin(7) * t570;
t254 = t344 * pkin(4) + t323;
t652 = t348 + t469;
t510 = -rSges(6,3) - pkin(3) - pkin(7);
t651 = t510 * t354 - pkin(2) - t568;
t118 = t214 * rSges(6,1) + t213 * rSges(6,2) + rSges(6,3) * t570;
t252 = t301 * t521;
t498 = t344 * t519;
t282 = pkin(7) * t498;
t493 = t344 * t346;
t40 = t493 + t118 * t329 + t382 * t248 - t252 - t282 + (t254 + t660) * qJD(1);
t322 = pkin(7) * t572;
t255 = pkin(4) * t345 - t322;
t477 = -t262 - t627;
t457 = -t238 + t477;
t39 = (t255 + t457) * qJD(1) + t657;
t606 = t345 * t39;
t650 = t344 * t40 + t606;
t548 = -Icges(5,3) * t572 + t185 - t311;
t550 = t427 * t344 + t183;
t641 = -t351 * t548 - t354 * t550;
t553 = -Icges(4,2) * t572 + t180 - t316;
t555 = t295 * t344 + t178;
t640 = -t351 * t553 - t354 * t555;
t258 = (Icges(6,2) * t350 - t593) * t354;
t372 = t248 * (-Icges(6,2) * t214 + t113 + t203) - t249 * (Icges(6,2) * t216 - t114 + t204) + t329 * (-t379 + t258);
t259 = (-Icges(6,1) * t353 + t594) * t354;
t373 = t248 * (-Icges(6,1) * t213 + t110 + t595) - t249 * (-Icges(6,1) * t215 + t112 - t205) + t329 * (-t378 - t259);
t639 = m(5) / 0.2e1;
t638 = m(6) / 0.2e1;
t511 = qJD(3) * qJD(5);
t488 = t351 * t511;
t171 = qJD(1) * t248 - t344 * t488;
t637 = t171 / 0.2e1;
t172 = qJD(1) * t249 - t345 * t488;
t636 = t172 / 0.2e1;
t635 = -t248 / 0.2e1;
t634 = t248 / 0.2e1;
t633 = -t249 / 0.2e1;
t632 = t249 / 0.2e1;
t631 = -t329 / 0.2e1;
t630 = t329 / 0.2e1;
t424 = t110 * t353 + t113 * t350;
t518 = qJD(3) * t354;
t380 = -t329 * t350 + t353 * t518;
t523 = qJD(1) * t351;
t460 = qJD(5) + t523;
t411 = t344 * t460;
t102 = t345 * t380 - t353 * t411;
t381 = t329 * t353 + t350 * t518;
t103 = t345 * t381 - t350 * t411;
t53 = Icges(6,5) * t103 + Icges(6,6) * t102 - Icges(6,3) * t659;
t55 = Icges(6,4) * t103 + Icges(6,2) * t102 - Icges(6,6) * t659;
t57 = Icges(6,1) * t103 + Icges(6,4) * t102 - Icges(6,5) * t659;
t8 = (qJD(3) * t424 + t53) * t351 + (qJD(3) * t107 - t350 * t57 - t353 * t55 + (t110 * t350 - t113 * t353) * qJD(5)) * t354;
t626 = t8 * t248;
t410 = t345 * t460;
t100 = t344 * t380 + t353 * t410;
t101 = t344 * t381 + t350 * t410;
t500 = t345 * t522;
t389 = -t498 + t500;
t52 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t389;
t54 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t389;
t56 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t389;
t9 = (qJD(3) * t423 + t52) * t351 + (qJD(3) * t109 - t350 * t56 - t353 * t54 + (t112 * t350 + t114 * t353) * qJD(5)) * t354;
t625 = t9 * t249;
t623 = qJD(1) / 0.2e1;
t217 = Icges(6,3) * t354 + t351 * t429;
t257 = (-Icges(6,5) * t353 + Icges(6,6) * t350) * t354;
t149 = qJD(3) * t217 + qJD(5) * t257;
t219 = Icges(6,6) * t354 + t351 * t430;
t150 = qJD(3) * t219 + qJD(5) * t258;
t221 = Icges(6,5) * t354 + t351 * t434;
t151 = qJD(3) * t221 + qJD(5) * t259;
t415 = -t350 * t379 - t353 * t378;
t31 = (qJD(3) * t415 + t149) * t351 + (-qJD(3) * t377 - t150 * t353 - t151 * t350 + (-t350 * t378 + t353 * t379) * qJD(5)) * t354;
t487 = t354 * t511;
t77 = -t351 * t377 - t354 * t415;
t622 = t31 * t329 + t77 * t487;
t621 = t103 * rSges(6,1) + t102 * rSges(6,2);
t620 = rSges(4,1) * t354;
t618 = rSges(5,2) * t354;
t617 = rSges(5,3) * t351;
t615 = rSges(6,3) * t354;
t303 = rSges(4,1) * t351 + rSges(4,2) * t354;
t241 = t303 * t345;
t499 = t303 * t521;
t86 = qJD(1) * t661 - t499;
t612 = t241 * t86;
t531 = rSges(4,2) * t573 + t345 * rSges(4,3);
t195 = rSges(4,1) * t572 - t531;
t496 = t303 * t520;
t85 = -t496 + (-t195 + t477) * qJD(1);
t607 = t344 * t85;
t605 = t345 * t85;
t42 = t107 * t351 - t354 * t424;
t603 = t42 * t172;
t602 = t43 * t171;
t601 = -rSges(5,3) - qJ(4);
t517 = qJD(4) * t354;
t406 = t238 * t521 + t242 * t520 + qJD(2) - t517;
t32 = t118 * t249 + t120 * t248 + (t254 * t345 - t255 * t344) * qJD(3) + t406;
t582 = qJD(3) * t32;
t211 = t344 * t518 + t345 * t523;
t283 = pkin(3) * t498;
t117 = pkin(3) * t500 + qJ(4) * t211 - t283 + t493;
t251 = t263 * qJD(1);
t562 = -t117 - t251;
t561 = t118 + t254;
t560 = t120 - t255;
t554 = -t295 * t345 - t179;
t552 = -t293 * t345 + t181;
t551 = -t427 * t345 + t182;
t549 = Icges(5,3) * t570 + t184 + t312;
t546 = -t197 - t242;
t545 = t344 * t238 + t345 * t242;
t239 = t301 * t345;
t544 = -qJD(1) * t239 + t344 * t517;
t541 = -t242 - t254;
t250 = qJD(3) * t305 - t517;
t306 = t617 - t618;
t540 = -t306 * qJD(3) - t250;
t494 = t345 * t518;
t539 = qJ(4) * t494 + t299;
t502 = t344 * t523;
t538 = rSges(4,2) * t502 + rSges(4,3) * t524;
t537 = -t425 + t428;
t536 = -t426 - t427;
t535 = -t293 + t296;
t534 = t295 + t433;
t533 = -t301 + t444;
t532 = -t305 - t306;
t530 = t344 ^ 2 + t345 ^ 2;
t527 = qJD(1) * t290;
t526 = qJD(1) * t432;
t512 = qJD(3) * qJD(4);
t509 = t357 * t627;
t508 = t357 * t348;
t116 = -pkin(3) * t659 - qJ(4) * t502 + t539;
t507 = t345 * t116 + t344 * t117 + t238 * t524;
t506 = t40 * t524;
t235 = t301 * t344;
t505 = -t235 * t521 - t239 * t520 + t346;
t327 = pkin(6) * t524;
t504 = t327 + t539;
t492 = t118 * t515;
t491 = t120 * t515;
t490 = -pkin(2) - t620;
t489 = t354 * t512;
t485 = t524 / 0.2e1;
t484 = -t521 / 0.2e1;
t481 = t520 / 0.2e1;
t479 = t518 / 0.2e1;
t475 = t339 - t627;
t470 = rSges(5,1) * t345 - rSges(5,3) * t573;
t468 = t533 * t345;
t463 = t117 * t521 + t351 * t512 + (t116 + t208) * t520;
t462 = qJD(1) * t235 + t345 * t517;
t459 = rSges(5,1) * t524 + rSges(5,2) * t659 + rSges(5,3) * t494;
t458 = qJD(1) * (-pkin(2) * t525 + t327) - t509;
t455 = t382 + t473;
t452 = qJD(5) * t479;
t261 = rSges(3,1) * t344 + rSges(3,2) * t345;
t448 = -rSges(4,2) * t351 + t620;
t447 = rSges(6,1) * t101 + rSges(6,2) * t100;
t443 = t344 * t36 + t345 * t35;
t442 = t344 * t35 - t345 * t36;
t441 = t344 * t38 + t345 * t37;
t440 = t344 * t37 - t345 * t38;
t439 = t344 * t43 + t345 * t42;
t438 = t344 * t42 - t345 * t43;
t437 = -t344 * t86 - t605;
t422 = t118 * t344 - t120 * t345;
t416 = t195 * t344 + t196 * t345;
t412 = qJD(1) * t252 + t345 * t489 - t508;
t243 = t351 * t445 + t615;
t260 = (-rSges(6,1) * t353 + rSges(6,2) * t350) * t354;
t152 = qJD(3) * t243 + qJD(5) * t260;
t408 = -pkin(7) * t518 - t152 - t250;
t407 = -pkin(2) - t305;
t403 = qJD(3) * t468 + t299;
t237 = t303 * t344;
t236 = t444 * t344;
t390 = t344 * t489 + t458 + (t116 + t299) * qJD(1);
t387 = -t107 * t248 + t109 * t249 + t329 * t377;
t386 = (Icges(6,5) * t213 - Icges(6,6) * t214) * t248 - (Icges(6,5) * t215 + Icges(6,6) * t216) * t249 + t257 * t329;
t385 = -t351 * t552 + t354 * t554;
t384 = t351 * t549 + t354 * t551;
t383 = t354 * t386;
t376 = (t351 * t536 + t354 * t537) * qJD(1);
t375 = (-t351 * t534 + t354 * t535) * qJD(1);
t137 = -rSges(4,1) * t659 - rSges(4,2) * t494 + t538;
t138 = -qJD(3) * t237 + (t345 * t448 + t333) * qJD(1);
t371 = t137 * t345 + t138 * t344 + (t195 * t345 - t196 * t344) * qJD(1);
t366 = t32 * t422 - (-t344 * t39 + t345 * t40) * t382;
t359 = (t377 * t345 + t424) * t248 - (t377 * t344 + t423) * t249 + (t217 + t415) * t329;
t358 = t359 * t354;
t356 = qJD(3) ^ 2;
t328 = pkin(4) * t524;
t271 = t448 * qJD(3);
t253 = t301 * t525;
t240 = t444 * t345;
t212 = t494 - t502;
t210 = t530 * t519;
t198 = rSges(5,2) * t572 + t470;
t174 = -pkin(7) * t659 + t328;
t173 = qJD(1) * t254 - t282;
t166 = t382 * t345;
t165 = t382 * t344;
t164 = t379 * t345;
t163 = t379 * t344;
t162 = t378 * t345;
t161 = t378 * t344;
t148 = rSges(6,1) * t215 + rSges(6,2) * t216;
t147 = rSges(6,1) * t213 - rSges(6,2) * t214;
t140 = -rSges(5,3) * t502 + t459;
t139 = qJD(3) * t236 + (t306 * t345 + t334) * qJD(1);
t84 = qJD(3) * t416 + qJD(2);
t65 = qJD(1) * t670 - t252 - t655;
t64 = (t198 + t457) * qJD(1) + t403;
t63 = -t508 - t271 * t520 + (-t138 - t251 + t499) * qJD(1);
t62 = -t271 * t521 + (t137 - t496) * qJD(1) + t458;
t60 = (t197 * t345 - t198 * t344) * qJD(3) + t406;
t59 = -rSges(6,3) * t659 + t621;
t58 = rSges(6,3) * t389 + t447;
t41 = t371 * qJD(3);
t34 = t540 * t520 + (-t139 + t562 + t655) * qJD(1) + t412;
t33 = qJD(1) * t140 + (qJD(1) * t468 + t344 * t540) * qJD(3) + t390;
t17 = (t139 * t344 + t140 * t345 + (-t198 * t345 + t344 * t546) * qJD(1)) * qJD(3) + t463;
t16 = -t102 * t378 - t103 * t379 + t149 * t570 + t150 * t213 + t151 * t214 + t377 * t659;
t15 = -t100 * t378 - t101 * t379 + t149 * t572 + t150 * t215 - t151 * t216 - t377 * t389;
t14 = t248 * t42 - t249 * t43 + t329 * t77;
t13 = -t356 * t323 - t152 * t249 - t171 * t382 - t329 * t58 + (-t250 * t345 - t491) * qJD(3) + (-t173 + (pkin(7) * qJD(3) - qJD(4)) * t573 + t562) * qJD(1) + t412;
t12 = -t356 * t322 + qJD(1) * t174 - t152 * t248 + t172 * t382 + t329 * t59 + (qJD(1) * t451 - t250 * t344 + t492) * qJD(3) + t390;
t7 = -t118 * t171 + t120 * t172 + t248 * t58 + t249 * t59 + (t173 * t344 + t174 * t345 + (-t255 * t345 + t344 * t541) * qJD(1)) * qJD(3) + t463;
t6 = t102 * t112 - t103 * t114 - t109 * t659 + t213 * t54 + t214 * t56 + t52 * t570;
t5 = t102 * t110 + t103 * t113 - t107 * t659 + t213 * t55 + t214 * t57 + t53 * t570;
t4 = t100 * t112 - t101 * t114 + t109 * t389 + t215 * t54 - t216 * t56 + t52 * t572;
t3 = t100 * t110 + t101 * t113 + t107 * t389 + t215 * t55 - t216 * t57 + t53 * t572;
t2 = t16 * t329 + t171 * t36 + t172 * t35 + t248 * t5 - t249 * t6 + t487 * t66;
t1 = t15 * t329 + t171 * t38 + t172 * t37 + t248 * t3 - t249 * t4 + t487 * t67;
t18 = [m(3) * ((-t261 * t357 - t509) * t652 + (-t508 + (-0.2e1 * t469 - t348 + t652) * t357) * (-t261 - t627)) + t626 / 0.2e1 - t625 / 0.2e1 + t622 + t602 / 0.2e1 + t603 / 0.2e1 + t15 * t633 + t16 * t634 + t66 * t636 + t67 * t637 + ((t671 * t344 + ((t696 + t711) * t345 + t678 + t700 + t710) * t345) * qJD(3) + t688) * t481 + (t632 + t633) * t10 + (t703 * qJD(3) + t707 * t351 + t708 * t354) * qJD(1) + (t12 * (t660 + t561) + t651 * t606 * qJD(1) + (t255 + t446 + t475 + (t407 - t615) * t344) * t13 + (t282 + t283 - t447 + (rSges(6,3) * t519 - qJ(4) * t518 - t346) * t344 + (-t348 + (-pkin(4) - pkin(6)) * t344) * qJD(1)) * t39 + (t510 * t495 + t328 + t39 + t504 + t621 - t654 - t657 + (t651 * t344 - t255) * qJD(1)) * t40) * m(6) + (-(-t64 + (t198 - t627) * qJD(1) + t403 + t654) * t65 + t34 * (t470 + t475) + t64 * t283 + t33 * t670 + t65 * (-pkin(3) * t495 + t459 + t504) + (t34 * (t407 + t618) + (-t346 + (t354 * t601 - t604) * qJD(3)) * t64) * t344 + ((-t352 * t65 - t355 * t64) * pkin(1) + t64 * (-pkin(2) + (rSges(5,2) - pkin(3)) * t354 + t601 * t351) * t345 + (t64 * (-rSges(5,1) - pkin(6)) + t65 * (t407 - t617)) * t344) * qJD(1)) * m(5) + (-(-t496 - t256 - t85 + (-t195 - t627) * qJD(1)) * t86 + t63 * (t344 * t490 + t475 + t531) + t62 * t661 + t86 * (t327 + t538) + (t303 * t607 - t612) * qJD(3) + ((-t352 * t86 - t355 * t85) * pkin(1) + (-pkin(2) - t448) * t605 + (t85 * (-rSges(4,3) - pkin(6)) + t86 * t490) * t344) * qJD(1)) * m(4) + (t680 + t681) * t521 / 0.2e1 + (-t514 + 0.2e1 * t513 + ((t345 * t672 - t671 + t676) * t345 + (t344 * t672 + t677 + t701) * t344) * qJD(3) + t684) * t484 - (t679 - t682 + t683) * t520 / 0.2e1 + ((t133 + t675) * t344 + (t674 + t698) * t345) * qJD(3) * t623; m(4) * t41 + m(5) * t17 + m(6) * t7; -t14 * t515 / 0.2e1 + t438 * t452 + (qJD(1) * t439 + t344 * t8 - t345 * t9) * t630 + (((-t162 * t353 - t164 * t350 + t107) * t248 - (-t161 * t353 - t163 * t350 + t109) * t249 + (-t219 * t353 - t221 * t350 - t377) * t329 + t77 * qJD(5)) * t354 + (-qJD(5) * t439 + t359) * t351) * t631 + ((t162 * t215 - t164 * t216) * t248 - (t161 * t215 - t163 * t216) * t249 + (t215 * t219 - t216 * t221) * t329 + (t354 * t67 - t37 * t571) * qJD(5) + ((-qJD(5) * t38 + t387) * t351 + t358) * t344) * t632 + (qJD(1) * t441 + t3 * t344 - t345 * t4) * t633 + (qJD(1) * t443 + t344 * t5 - t345 * t6) * t634 + ((t162 * t213 + t164 * t214) * t248 - (t161 * t213 + t163 * t214) * t249 + (t213 * t219 + t214 * t221) * t329 + (t354 * t66 - t36 * t573) * qJD(5) + ((-qJD(5) * t35 + t387) * t351 + t358) * t345) * t635 + t442 * t636 + t440 * t637 - ((((-t548 - t553) * t345 + (-t549 + t552) * t344) * t354 + ((t550 + t555) * t345 + (t551 + t554) * t344) * t351) * qJD(3) + ((t534 - t536) * t354 + (t535 + t537) * t351) * qJD(1)) * qJD(1) / 0.2e1 + (t682 * t345 + t681 * t344 + (t344 * t675 + t345 * t674) * qJD(1)) * t623 + ((-t521 * t574 + t527) * t344 + (t375 + (-t640 * t345 + (t575 + t385) * t344) * qJD(3)) * t345 + (t230 * t521 - t526) * t344 + (t376 + (-t641 * t345 + (-t229 + t384) * t344) * qJD(3)) * t345) * t484 + ((-t520 * t575 - t527) * t345 + (t375 + (t385 * t344 + (t574 - t640) * t345) * qJD(3)) * t344 + (t229 * t520 + t526) * t345 + (t376 + (t384 * t344 + (-t230 - t641) * t345) * qJD(3)) * t344) * t481 + t658 * t516 / 0.2e1 + (-t40 * (t166 * t329 - t243 * t248 + t492 + t544) - ((-t530 * t582 - t506) * pkin(7) + t366 * qJD(5)) * t351 - t650 * qJD(3) * (-pkin(7) * t354 - t305) + t7 * t545 + (t12 * t455 + t40 * t408 + t7 * t560) * t344 + (t7 * t561 + (qJD(1) * t40 + t13) * t455) * t345 + (t165 * t329 + t243 * t249 + t345 * t408 - t382 * t525 + t253 - t462 + t491) * t39 + (-t165 * t248 - t166 * t249 - t505 + t507 + (t173 + t58 + (-t118 + t541) * qJD(1)) * t344 + (qJD(1) * t560 + t174 + t59) * t345) * t32) * m(6) + (-t64 * (-qJD(1) * t236 + t462) - t65 * (qJD(1) * t240 + t544) - t60 * t505 - ((t60 * t240 + t532 * t64) * t345 + (t60 * t236 + t532 * t65) * t344) * qJD(3) + t64 * t253 + t17 * t545 + t60 * t507 + (t34 * t533 + t64 * t540 + t17 * t197 + t60 * t140 + (-t60 * t198 + t533 * t65) * qJD(1)) * t345 + (t33 * t533 + t65 * t540 - t17 * t198 + t60 * t139 + (-t444 * t64 + t546 * t60) * qJD(1)) * t344) * m(5) + (t41 * t416 + t84 * t371 + t437 * t271 + (-t62 * t344 - t63 * t345 + (-t345 * t86 + t607) * qJD(1)) * t303 - (t237 * t85 - t612) * qJD(1) - (t84 * (-t237 * t344 - t241 * t345) + t437 * t448) * qJD(3)) * m(4) + (t2 + t680 * qJD(1) + (t676 * t524 + (qJD(1) * t677 + t667 * t344 - t673) * t344) * t665) * t344 / 0.2e1 - (t1 + t679 * qJD(1) + ((qJD(1) * t678 + t673) * t345 + (qJD(1) * t699 - t667 * t345) * t344) * t665) * t345 / 0.2e1 + (t11 + t684 + t686) * t525 / 0.2e1 + (t10 + t683 + t685) * t485; -m(5) * (t210 * t60 + t211 * t65 + t212 * t64) - m(6) * (t210 * t32 + t211 * t40 + t212 * t39) + 0.2e1 * ((t520 * t64 + t521 * t65 - t17) * t639 + (t39 * t520 + t40 * t521 - t7) * t638) * t354 + 0.2e1 * ((qJD(3) * t60 + t33 * t344 + t34 * t345 + t524 * t65 - t525 * t64) * t639 + (t12 * t344 + t13 * t345 - t39 * t525 + t506 + t582) * t638) * t351; -t10 * t501 / 0.2e1 + t2 * t570 / 0.2e1 + (t351 * t66 + t354 * t443) * t636 + ((-qJD(3) * t443 + t16) * t351 + (-qJD(1) * t442 + qJD(3) * t66 + t344 * t6 + t345 * t5) * t354) * t634 + t354 * t11 * t485 + t1 * t572 / 0.2e1 + (t351 * t67 + t354 * t441) * t637 + ((-qJD(3) * t441 + t15) * t351 + (-qJD(1) * t440 + qJD(3) * t67 + t3 * t345 + t344 * t4) * t354) * t633 + t14 * t479 + t351 * (t602 + t603 + t622 - t625 + t626) / 0.2e1 + (t351 * t77 + t354 * t439) * t452 + ((-qJD(3) * t439 + t31) * t351 + (-qJD(1) * t438 + qJD(3) * t77 + t344 * t9 + t345 * t8) * t354) * t630 + (t372 * t213 - t214 * t373 + t345 * t383) * t635 + (t215 * t372 + t216 * t373 + t344 * t383) * t632 + (t386 * t351 + (t373 * t350 - t353 * t372) * t354) * t631 - t658 * t519 / 0.2e1 + ((qJD(3) * t366 + t12 * t118 - t13 * t120 - t39 * t58 + t40 * t59) * t351 + (t39 * (-qJD(3) * t120 + t152 * t344) + t40 * (qJD(3) * t118 - t152 * t345) - t7 * t422 + t32 * (-t118 * t524 - t120 * t525 - t344 * t59 + t345 * t58) - (qJD(1) * t650 - t12 * t345 + t13 * t344) * t382) * t354 - t39 * (-t148 * t329 - t249 * t260) - t40 * (t147 * t329 - t248 * t260) - t32 * (t147 * t249 + t148 * t248)) * m(6);];
tauc = t18(:);
