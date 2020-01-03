% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR11
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR11_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:12
% EndTime: 2019-12-31 18:27:29
% DurationCPUTime: 11.80s
% Computational Cost: add. (28381->507), mult. (36261->677), div. (0->0), fcn. (39573->7), ass. (0->322)
t368 = sin(qJ(1));
t573 = -t368 / 0.4e1;
t364 = pkin(8) + qJ(3);
t352 = sin(t364);
t353 = cos(t364);
t565 = sin(qJ(5));
t566 = cos(qJ(5));
t310 = t352 * t566 - t353 * t565;
t385 = t352 * t565 + t353 * t566;
t219 = -Icges(6,5) * t385 - Icges(6,6) * t310;
t369 = cos(qJ(1));
t286 = t310 * t369;
t287 = t385 * t369;
t308 = Icges(6,4) * t310;
t223 = -Icges(6,2) * t385 + t308;
t227 = Icges(6,1) * t385 + t308;
t486 = t223 + t227;
t307 = Icges(6,4) * t385;
t224 = Icges(6,2) * t310 + t307;
t226 = Icges(6,1) * t310 - t307;
t488 = t224 - t226;
t79 = t219 * t368 + t488 * t286 + t486 * t287;
t259 = Icges(6,4) * t287;
t191 = Icges(6,2) * t286 - Icges(6,6) * t368 + t259;
t494 = -Icges(6,1) * t286 + t191 + t259;
t258 = Icges(6,4) * t286;
t193 = Icges(6,1) * t287 - Icges(6,5) * t368 + t258;
t604 = -Icges(6,2) * t287 + t193 + t258;
t88 = t494 * t310 + t385 * t604;
t679 = (t79 + t88) * t573;
t667 = (t227 / 0.2e1 + t223 / 0.2e1) * t310;
t678 = (-t226 / 0.2e1 + t224 / 0.2e1) * t385 - t667;
t284 = t310 * t368;
t285 = t385 * t368;
t77 = -t219 * t369 + t488 * t284 + t486 * t285;
t677 = -t77 / 0.2e1;
t676 = -t77 / 0.4e1;
t674 = t88 / 0.2e1 + t79 / 0.2e1;
t257 = Icges(6,4) * t285;
t522 = Icges(6,2) * t284;
t409 = t522 + t257;
t519 = Icges(6,6) * t369;
t381 = t519 + t409;
t256 = Icges(6,4) * t284;
t529 = Icges(6,1) * t285;
t413 = t256 + t529;
t524 = Icges(6,5) * t369;
t382 = t524 + t413;
t398 = -t284 * t191 - t285 * t193;
t405 = Icges(6,5) * t285 + Icges(6,6) * t284;
t379 = Icges(6,3) * t369 + t405;
t612 = t379 - t405;
t366 = t369 ^ 2;
t464 = 0.2e1 * t257;
t408 = t464 + 0.2e1 * t519;
t463 = 0.2e1 * t524;
t601 = -(t408 + t522) * t284 - (t463 + t529) * t285 - Icges(6,3) * t366;
t658 = t601 * t369;
t15 = -t658 + (t612 * t368 + (-t382 + t413) * t287 - (t381 - t409) * t286 + t398) * t368;
t378 = Icges(6,5) * t287 + Icges(6,6) * t286 - Icges(6,3) * t368;
t375 = t369 * t378;
t100 = t375 - t398;
t101 = t286 * t191 + t287 * t193 - t368 * t378;
t16 = (-t286 * t409 - t287 * t413 + t100 - 0.2e1 * t375 + t398) * t369 + (-t284 * t381 - t285 * t382 - t612 * t369 + t101 - t601) * t368;
t51 = t100 * t368 + t658;
t52 = t101 * t368 - (t286 * t381 + t287 * t382 - t368 * t379) * t369;
t568 = t369 / 0.4e1;
t569 = -t369 / 0.4e1;
t673 = 0.2e1 * t16 * t568 + 0.2e1 * t52 * t569 + 0.2e1 * (t15 + t51) * t573 + (t79 / 0.4e1 + t88 / 0.4e1) * t368;
t230 = rSges(6,1) * t310 - rSges(6,2) * t385;
t516 = qJ(4) * t353;
t323 = pkin(3) * t352 - t516;
t546 = pkin(4) * t352;
t419 = t230 + t323 + t546;
t174 = t419 * t368;
t176 = t419 * t369;
t541 = rSges(5,1) * t352;
t324 = -rSges(5,3) * t353 + t541;
t473 = t323 + t324;
t235 = t473 * t368;
t238 = t473 * t369;
t362 = t369 * rSges(5,2);
t545 = -pkin(6) - qJ(2);
t445 = t369 * t545;
t350 = cos(pkin(8)) * pkin(2) + pkin(1);
t533 = rSges(5,3) + qJ(4);
t567 = rSges(5,1) + pkin(3);
t600 = t533 * t352 + t567 * t353 + t350;
t185 = -t600 * t368 + t362 - t445;
t351 = t368 * t545;
t540 = t368 * rSges(5,2);
t186 = t600 * t369 - t351 + t540;
t504 = t353 * t369;
t505 = t353 * t368;
t496 = t185 * t504 + t186 * t505;
t510 = t352 * qJ(4);
t594 = pkin(3) + pkin(4);
t603 = t594 * t353 + t350 + t510;
t196 = t287 * rSges(6,1) + t286 * rSges(6,2) - t368 * rSges(6,3);
t611 = -t368 * pkin(7) + t196;
t145 = t603 * t369 - t351 + t611;
t610 = -t285 * rSges(6,1) - t284 * rSges(6,2);
t638 = -t603 * t368 + (-pkin(7) - rSges(6,3) - t545) * t369 + t610;
t499 = t145 * t505 + t504 * t638;
t508 = t352 * t369;
t509 = t352 * t368;
t595 = m(6) / 0.2e1;
t596 = m(5) / 0.2e1;
t543 = (-t174 * t508 + t176 * t509 + t499) * t595 + (-t235 * t508 + t238 * t509 + t496) * t596;
t343 = pkin(3) * t509;
t232 = t343 + (-t533 * t353 + t541) * t368;
t335 = qJ(4) * t504;
t342 = rSges(5,3) * t504;
t233 = -t567 * t508 + t335 + t342;
t205 = -rSges(6,1) * t284 + rSges(6,2) * t285;
t171 = t343 + (-t516 + t546) * t368 - t205;
t208 = t286 * rSges(6,1) - t287 * rSges(6,2);
t172 = -t594 * t508 - t208 + t335;
t401 = t171 * t369 + t172 * t368;
t544 = (t401 * t352 + t499) * t595 + ((t232 * t369 + t233 * t368) * t352 + t496) * t596;
t633 = t543 - t544;
t672 = t633 * qJD(1);
t348 = Icges(5,5) * t352;
t530 = Icges(5,1) * t353;
t414 = t348 + t530;
t269 = Icges(5,4) * t368 + t414 * t369;
t526 = Icges(4,4) * t352;
t321 = Icges(4,1) * t353 - t526;
t271 = Icges(4,5) * t368 + t321 * t369;
t668 = t269 + t271;
t501 = Icges(6,1) - Icges(6,2);
t440 = t501 * t285;
t387 = t440 + t524;
t629 = -0.2e1 * t256;
t656 = t387 - t629;
t666 = t656 * t286 + t368 * (-Icges(6,5) * t284 + Icges(6,6) * t285);
t648 = -m(6) / 0.2e1;
t194 = t369 * rSges(6,3) - t610;
t133 = -t368 * t194 - t369 * t196;
t652 = -t368 * t205 + t208 * t369;
t664 = t133 * t652;
t400 = t174 * t368 + t176 * t369;
t325 = rSges(4,1) * t352 + rSges(4,2) * t353;
t365 = t368 ^ 2;
t469 = t365 + t366;
t616 = t469 * t325;
t452 = t400 * t648 + (-t235 * t368 - t238 * t369) * t596 - m(4) * t616 / 0.2e1;
t304 = t325 * t368;
t306 = t325 * t369;
t211 = t368 * t304 + t306 * t369;
t453 = (t368 * t171 - t172 * t369) * t595 + (t368 * t232 - t233 * t369) * t596 + m(4) * t211 / 0.2e1;
t24 = t453 - t452;
t663 = t24 * qJD(1);
t199 = -Icges(6,5) * t286 + Icges(6,6) * t287;
t662 = (-t199 * t368 - t286 * t604 + t494 * t287) * t368;
t661 = (t199 * t369 - t284 * t604 + t494 * t285) * t368;
t660 = t652 * t353;
t229 = -rSges(6,1) * t385 - rSges(6,2) * t310;
t615 = t469 * t352;
t550 = m(6) * (t229 * t615 + t660);
t338 = Icges(5,5) * t504;
t261 = Icges(5,6) * t368 + Icges(5,3) * t508 + t338;
t314 = Icges(4,5) * t353 - Icges(4,6) * t352;
t263 = Icges(4,3) * t368 + t314 * t369;
t315 = Icges(5,4) * t353 + Icges(5,6) * t352;
t265 = Icges(5,2) * t368 + t315 * t369;
t657 = t261 * t508 + t668 * t504 + (t263 + t265) * t368;
t655 = (-Icges(4,6) + Icges(5,6)) * t353 + (-Icges(5,4) - Icges(4,5)) * t352;
t316 = Icges(4,2) * t353 + t526;
t518 = Icges(5,3) * t353;
t406 = t518 - t348;
t654 = (-t316 - t406) * t369 + t668;
t613 = t145 * t368 + t369 * t638;
t649 = -t488 * t385 / 0.2e1 + t667;
t572 = t368 / 0.2e1;
t623 = -t369 / 0.2e1;
t644 = t229 * t613;
t266 = Icges(4,4) * t505 - Icges(4,2) * t509 - Icges(4,6) * t369;
t349 = Icges(4,4) * t353;
t523 = Icges(4,2) * t352;
t267 = Icges(4,6) * t368 + (t349 - t523) * t369;
t242 = t271 * t505;
t424 = t263 * t369 - t242;
t262 = Icges(4,5) * t505 - Icges(4,6) * t509 - Icges(4,3) * t369;
t339 = Icges(4,4) * t509;
t270 = Icges(4,1) * t505 - Icges(4,5) * t369 - t339;
t484 = -t368 * t262 - t270 * t504;
t640 = -t266 * t508 - t267 * t509 - t424 - t484;
t639 = -t267 * t508 + t657;
t512 = (-Icges(5,2) * t369 + t368 * t315) * t369;
t636 = t512 + t657;
t439 = t501 * t284;
t635 = t439 - t519;
t634 = t469 * t353;
t628 = -0.2e1 * t615;
t624 = -t368 / 0.2e1;
t525 = Icges(5,5) * t353;
t313 = Icges(5,3) * t352 + t525;
t260 = -Icges(5,6) * t369 + t313 * t368;
t268 = -Icges(5,4) * t369 + t414 * t368;
t618 = (t260 * t352 + t268 * t353) * t368;
t462 = m(5) / 0.4e1 + m(6) / 0.4e1;
t471 = t634 * t352;
t617 = t462 * (-t352 * t353 + t471);
t429 = t469 * t230;
t609 = t655 * t368;
t608 = t655 * t369;
t607 = t654 * t368;
t531 = Icges(4,1) * t352;
t415 = -t349 - t531;
t435 = (t415 * t369 - t267) * t368;
t437 = (-Icges(5,1) * t508 + t261 + t338) * t368;
t606 = t435 + t437;
t426 = -0.2e1 * t257;
t372 = t426 + t635;
t393 = (-t662 - (t372 * t287 + t666) * t369) * t572 + (-t661 - ((t426 - 0.2e1 * t519) * t285 - (-0.2e1 * t440 + t629 - 0.2e1 * t524) * t284) * t369) * t623;
t374 = t464 - t635;
t394 = (t662 - (t374 * t287 - t666) * t369) * t572 + (t661 - (-(-t629 + t463) * t284 + (t408 - 0.2e1 * t439) * t285) * t369) * t623;
t318 = Icges(5,1) * t352 - t525;
t602 = (t321 / 0.2e1 - t316 / 0.2e1 + t348 + t530 / 0.2e1 - t518 / 0.2e1) * t352 + (t349 + t531 / 0.2e1 - t523 / 0.2e1 + t318 / 0.2e1 - t313 / 0.2e1) * t353;
t599 = 0.4e1 * qJD(1);
t598 = 2 * qJD(3);
t588 = m(6) * (-t174 * t208 - t176 * t205 - t644);
t587 = m(6) * (t401 * t230 - t644);
t428 = t469 * t229;
t490 = t634 * t230;
t586 = m(6) * (-t660 + (t133 - t428) * t352 + t490);
t326 = t353 * pkin(3) + t510;
t474 = t469 * t326;
t109 = (pkin(4) * t504 + t611) * t369 + (pkin(4) * t505 + t369 * pkin(7) + t194) * t368 + t474;
t497 = -t174 * t505 - t176 * t504;
t583 = m(6) * (t109 * t615 + t497);
t582 = m(6) * (t145 * t172 + t171 * t638);
t580 = m(6) * (t145 * t208 + t205 * t638);
t579 = m(6) * (t133 * t615 + t490);
t578 = m(6) * (t145 * t508 - t509 * t638);
t564 = m(3) * t469 * (rSges(3,3) + qJ(2));
t542 = rSges(4,1) * t353;
t442 = t350 + t542;
t470 = rSges(4,2) * t509 + t369 * rSges(4,3);
t217 = -t442 * t368 - t445 + t470;
t441 = -rSges(4,2) * t508 + t368 * rSges(4,3);
t218 = t442 * t369 - t351 + t441;
t563 = m(4) * (t217 * t304 - t218 * t306);
t562 = m(4) * (t217 * t369 + t218 * t368);
t327 = t353 * rSges(5,1) + t352 * rSges(5,3);
t146 = t368 * (t327 * t368 - t362) + (t327 * t369 + t540) * t369 + t474;
t489 = -t235 * t505 - t238 * t504;
t559 = m(5) * (t146 * t615 + t489);
t557 = m(5) * (t185 * t232 + t186 * t233);
t556 = m(5) * (-t185 * t509 + t186 * t508);
t555 = m(5) * (t185 * t369 + t186 * t368);
t551 = m(6) * t613;
t511 = t266 * t352;
t388 = (t205 * t369 + t208 * t368) * t352;
t89 = t388 * t648;
t502 = t89 * qJD(4);
t482 = -t318 * t368 + t260;
t481 = -t415 * t368 + t266;
t480 = -t406 * t368 + t268;
t478 = -Icges(4,2) * t505 + t270 - t339;
t475 = t368 * (qJ(4) * t505 - t343) + t369 * (-pkin(3) * t508 + t335);
t472 = -t326 - t327;
t384 = t652 * t648;
t422 = m(6) * t429;
t107 = t384 - t422 / 0.2e1;
t467 = t107 * qJD(1);
t461 = t595 + t596;
t160 = t461 * t628;
t466 = t160 * qJD(1);
t459 = t15 / 0.2e1 + t51 / 0.2e1;
t458 = t52 / 0.2e1 - t16 / 0.2e1;
t448 = t315 / 0.2e1 + t314 / 0.2e1;
t443 = -t394 - t393;
t438 = t482 * t369;
t436 = t481 * t369;
t434 = t480 * t369;
t432 = t478 * t369;
t423 = t267 * t352 - t262;
t85 = t372 * t310 - (t387 + 0.2e1 * t256) * t385;
t421 = t588 / 0.2e1 + t679 + (-t77 + t85) * t569;
t87 = t374 * t310 + t385 * t656;
t420 = t587 / 0.2e1 + t679 + (t77 + t87) * t568;
t418 = -pkin(4) * t353 + t229 - t326;
t417 = -t261 * t509 + t265 * t369 - t269 * t505;
t175 = t418 * t368;
t177 = t418 * t369;
t399 = t175 * t368 + t177 * t369;
t148 = -t512 + t618;
t377 = t458 + (t148 - t618 + t636) * t624 + t639 * t572 + (-t242 + (t263 + t511) * t369 + t484 + t640) * t623;
t376 = t417 * t624 + t459 + (-(-t270 * t353 + t511) * t368 + t148 - t262 * t369) * t623 + (t423 * t369 - t636 + t639) * t369 / 0.2e1 + (t260 * t508 + t268 * t504 + t423 * t368 + t417 + t424 + t640) * t572;
t373 = t459 * t368 + t458 * t369;
t328 = -rSges(4,2) * t352 + t542;
t239 = t472 * t369;
t236 = t472 * t368;
t162 = 0.4e1 * t617;
t159 = t462 * t628 + (m(5) + m(6)) * t615 / 0.2e1;
t156 = t369 * (-rSges(5,1) * t508 + t342) - t324 * t365 + t475;
t117 = -pkin(4) * t615 + t475 - t652;
t106 = t384 + t422 / 0.2e1;
t102 = -t550 / 0.2e1;
t91 = t579 / 0.2e1;
t90 = t388 * t595;
t63 = t556 + t578;
t59 = -t230 * t428 + t664;
t55 = t586 / 0.2e1;
t46 = t551 + t555 + t562 + t564;
t40 = t580 + t678;
t35 = t559 + t583;
t23 = t452 + t453;
t20 = t91 + t55 + t550 / 0.2e1;
t19 = t102 + t91 - t586 / 0.2e1;
t18 = t102 + t55 - t579 / 0.2e1;
t17 = t563 + t557 + t582 + t602 - t678;
t9 = t543 + t544;
t7 = m(6) * t59 + t394;
t6 = m(6) * (t109 * t652 + t400 * t229) + t393;
t4 = t376 * t368 + t377 * t369;
t3 = -t588 / 0.2e1 + (t676 + t85 / 0.4e1) * t369 + t420 + t673;
t2 = t373 + t420 + t421;
t1 = -t587 / 0.2e1 + (t676 - t87 / 0.4e1) * t369 + t421 + t673;
t5 = [t46 * qJD(2) + t17 * qJD(3) + t63 * qJD(4) + t40 * qJD(5), qJD(1) * t46 + qJD(3) * t23 + qJD(4) * t159 + qJD(5) * t106, t17 * qJD(1) + t23 * qJD(2) + t9 * qJD(4) + t2 * qJD(5) + ((t185 * t239 + t186 * t236 - t232 * t238 - t233 * t235) * t596 + (t145 * t175 - t171 * t176 - t172 * t174 + t177 * t638) * t595) * t598 + ((m(4) * (-t217 * t328 - t304 * t325) - t87 / 0.2e1 + t677 + t448 * t369 - t377) * t369 + (m(4) * (-t218 * t328 + t306 * t325) + t448 * t368 - t376 + t674) * t368 + (-t438 / 0.2e1 + t436 / 0.2e1 + t437 / 0.2e1 + t435 / 0.2e1) * t352 + (-t434 / 0.2e1 - t432 / 0.2e1 + t654 * t572) * t353) * qJD(3), qJD(1) * t63 + qJD(2) * t159 + qJD(3) * t9 - qJD(5) * t89, t40 * qJD(1) + t106 * qJD(2) + t2 * qJD(3) - t502 + ((m(6) * (t205 * t230 + t229 * t638) + t677 + t85 / 0.2e1 - t458) * t369 + (m(6) * (t145 * t229 + t208 * t230) - t459 + t674) * t368) * qJD(5); t24 * qJD(3) + t160 * qJD(4) + t107 * qJD(5) + (-t551 / 0.4e1 - t555 / 0.4e1 - t562 / 0.4e1 - t564 / 0.4e1) * t599, 0, t663 + ((-t236 * t369 + t239 * t368) * t596 + (-t175 * t369 + t177 * t368) * t595) * t598, t466, t467; -t24 * qJD(2) + t4 * qJD(3) + t633 * qJD(4) + t1 * qJD(5) + (-t557 / 0.4e1 - t563 / 0.4e1 - t582 / 0.4e1) * t599 + (-t602 - t649) * qJD(1), -t663, t4 * qJD(1) + (m(6) * (t109 * t117 - t174 * t175 - t176 * t177) + m(5) * (t146 * t156 - t235 * t236 - t238 * t239) + m(4) * (t328 * t616 - (t368 * (rSges(4,1) * t505 - t470) + t369 * (rSges(4,1) * t504 + t441)) * t211) + t394 + ((-t609 * t368 + (t432 + t434 - t607) * t352 + ((t481 - t482) * t369 + t606) * t353) * t369 + t608 * t365) * t572 + ((-t608 * t369 + (t436 - t438 + t606) * t353 + ((t478 + t480) * t369 - t607) * t352) * t368 + t609 * t366) * t623) * qJD(3) + t35 * qJD(4) + t6 * qJD(5), t672 + t35 * qJD(3) + t19 * qJD(5) + (-0.4e1 * t617 + 0.2e1 * t461 * (-t353 * t615 + t471)) * qJD(4), t1 * qJD(1) + t6 * qJD(3) + t19 * qJD(4) + ((-t59 + (-t109 + t133) * t652 + (-t400 - t429) * t229) * m(6) - t393 + t443) * qJD(5); -t160 * qJD(2) - t633 * qJD(3) + t90 * qJD(5) + (-t578 / 0.4e1 - t556 / 0.4e1) * t599, -t466, -t672 + t162 * qJD(4) + t18 * qJD(5) + 0.4e1 * (-t583 / 0.4e1 - t559 / 0.4e1) * qJD(3) + ((-t353 * t117 + t497) * t595 + (-t353 * t156 + t489) * t596 + ((t109 + t399) * t595 + (t236 * t368 + t239 * t369 + t146) * t596) * t352) * t598, t162 * qJD(3), t90 * qJD(1) + t18 * qJD(3) + qJD(5) * t550; (-t580 + t649) * qJD(1) - t107 * qJD(2) + t3 * qJD(3) + t502 + t373 * qJD(5), -t467, t3 * qJD(1) + ((t133 * t117 + t230 * t399) * m(6) - t394 + t443) * qJD(3) + t20 * qJD(4) + t7 * qJD(5), qJD(1) * t89 + qJD(3) * t20, t373 * qJD(1) + t7 * qJD(3) + (m(6) * (t229 * t429 - t664) + t393) * qJD(5);];
Cq = t5;
