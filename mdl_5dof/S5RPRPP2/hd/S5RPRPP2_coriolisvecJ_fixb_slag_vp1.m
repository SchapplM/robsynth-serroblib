% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:11:08
% DurationCPUTime: 23.62s
% Computational Cost: add. (9311->547), mult. (11747->657), div. (0->0), fcn. (8998->6), ass. (0->308)
t661 = Icges(5,1) + Icges(6,1);
t660 = Icges(6,2) + Icges(5,3);
t659 = Icges(5,6) - Icges(6,6);
t658 = Icges(5,4) - Icges(6,5);
t306 = sin(qJ(3));
t302 = Icges(6,4) * t306;
t308 = cos(qJ(3));
t380 = Icges(6,1) * t308 + t302;
t301 = Icges(5,5) * t306;
t381 = Icges(5,1) * t308 + t301;
t657 = t380 + t381;
t527 = Icges(4,4) * t306;
t250 = Icges(4,1) * t308 - t527;
t305 = qJ(1) + pkin(7);
t298 = sin(t305);
t299 = cos(t305);
t147 = Icges(4,5) * t298 + t250 * t299;
t651 = t658 * t298 + t657 * t299;
t627 = t147 + t651;
t614 = -t660 * t308 + t301 + t302;
t632 = -Icges(4,2) * t308 - t527 + t614;
t493 = t299 * t308;
t656 = (Icges(6,4) + Icges(5,5)) * t493;
t303 = Icges(4,4) * t308;
t379 = -Icges(4,2) * t306 + t303;
t523 = Icges(5,5) * t308;
t236 = Icges(5,3) * t306 + t523;
t525 = Icges(6,4) * t308;
t240 = Icges(6,2) * t306 + t525;
t633 = t236 + t240;
t655 = t379 - t633;
t654 = t659 * t298;
t249 = Icges(4,1) * t306 + t303;
t603 = t661 * t306 + t249 - t523 - t525;
t494 = t299 * t306;
t647 = t660 * t494 + t654 + t656;
t238 = Icges(4,5) * t308 - Icges(4,6) * t306;
t135 = Icges(4,3) * t298 + t238 * t299;
t242 = Icges(5,4) * t308 + Icges(5,6) * t306;
t139 = Icges(5,2) * t298 + t242 * t299;
t653 = t135 + t139;
t641 = (Icges(4,6) - t659) * t308 + (Icges(4,5) + t658) * t306;
t652 = t250 + t657;
t650 = -t306 * t632 - t308 * t603;
t495 = t298 * t308;
t496 = t298 * t306;
t513 = Icges(4,3) * t299;
t134 = Icges(4,5) * t495 - Icges(4,6) * t496 - t513;
t269 = Icges(4,4) * t496;
t524 = Icges(4,5) * t299;
t146 = Icges(4,1) * t495 - t269 - t524;
t518 = Icges(4,6) * t299;
t140 = Icges(4,4) * t495 - Icges(4,2) * t496 - t518;
t504 = t140 * t306;
t368 = -t146 * t308 + t504;
t136 = Icges(6,6) * t299 + t240 * t298;
t142 = Icges(6,5) * t299 + t298 * t380;
t371 = t136 * t306 + t142 * t308;
t138 = -Icges(5,2) * t299 + t242 * t298;
t505 = t138 * t299;
t234 = Icges(6,5) * t308 + Icges(6,6) * t306;
t130 = Icges(6,3) * t299 + t234 * t298;
t507 = t130 * t299;
t132 = -Icges(5,6) * t299 + t236 * t298;
t144 = -Icges(5,4) * t299 + t298 * t381;
t375 = t132 * t306 + t144 * t308;
t582 = t298 * t375;
t616 = t298 * t371 - t505 + t507 + t582;
t649 = -t134 * t299 - t298 * t368 + t616;
t141 = Icges(4,6) * t298 + t299 * t379;
t113 = t147 * t495;
t402 = t135 * t299 - t113;
t51 = -t141 * t496 - t402;
t131 = -Icges(6,3) * t298 + t234 * t299;
t128 = t299 * t131;
t615 = -t139 * t299 + t651 * t495 + t647 * t496 + t128;
t648 = t51 + t615;
t646 = -t298 * t134 - t136 * t494 + (-t142 - t146) * t493;
t645 = t653 * t298 + t493 * t627 + t647 * t494;
t644 = t655 * qJD(3);
t643 = t652 * qJD(3);
t642 = t234 - t238 - t242;
t571 = t641 * t299;
t572 = t641 * t298;
t508 = t130 * t298;
t638 = t140 * t494 + t508 + t646;
t606 = -t131 * t298 - t141 * t494 + t645;
t637 = t650 * t298 + t571;
t636 = -t650 * t299 + t572;
t635 = -t132 - t136;
t634 = t142 + t144;
t631 = t368 - t375;
t630 = t140 + t635;
t629 = t141 - t647;
t628 = t146 + t634;
t626 = t643 * t308 - t644 * t306 + (-t603 * t306 + t308 * t632) * qJD(3) + t641 * qJD(1);
t625 = t632 * qJD(3);
t624 = t603 * qJD(3);
t503 = t141 * t306;
t623 = t647 * t306 + t627 * t308 - t503;
t622 = -qJD(1) * t650 + t642 * qJD(3);
t621 = -t131 - t371;
t620 = t636 * qJD(1);
t126 = t298 * t138;
t54 = t132 * t494 + t144 * t493 + t126;
t593 = t299 * t54;
t619 = (t606 * t298 + t638 * t299 - t593) * qJD(3);
t618 = (t648 * t298 - t649 * t299) * qJD(3);
t617 = t637 * qJD(1);
t532 = -rSges(6,3) - qJ(5);
t613 = -t617 + t618;
t612 = t619 + t620;
t611 = (-t625 * t298 + (t633 * t299 - t141 + t654) * qJD(1)) * t308 + (-t627 * qJD(1) + t624 * t298) * t306 + (-t371 + t631) * qJD(3);
t610 = (t625 * t299 + (-t298 * t379 + t518 - t635) * qJD(1)) * t308 + (-t624 * t299 + (-t250 * t298 + t524 - t634) * qJD(1)) * t306 + t623 * qJD(3);
t609 = -t622 * t298 + t626 * t299;
t608 = t626 * t298 + t622 * t299;
t607 = t54 - t638;
t605 = t627 * t306 + t629 * t308;
t604 = t628 * t306 + t630 * t308;
t602 = t641 * qJD(3);
t601 = t505 + t645;
t432 = -rSges(6,1) - pkin(3) - pkin(4);
t510 = qJ(4) * t306;
t541 = rSges(6,2) * t306;
t600 = qJD(1) * (t432 * t308 - pkin(2) - t510 - t541) - qJD(5);
t599 = (-t602 * t298 + (t621 + t631 + t653) * qJD(1)) * t299;
t288 = t298 * rSges(5,2);
t159 = rSges(5,1) * t493 + rSges(5,3) * t494 + t288;
t277 = pkin(3) * t493;
t197 = qJ(4) * t494 + t277;
t209 = t299 * pkin(2) + t298 * pkin(6);
t309 = cos(qJ(1));
t304 = t309 * pkin(1);
t578 = t304 + t209;
t591 = t197 + t578;
t598 = t591 + t159;
t442 = qJD(1) * t299;
t597 = t603 + t655;
t596 = t632 + t652;
t595 = t632 * t299 + t627;
t594 = -t249 * t299 - t661 * t494 - t629 + t656;
t311 = qJD(1) ^ 2;
t287 = t298 * rSges(4,3);
t160 = rSges(4,1) * t493 - rSges(4,2) * t494 + t287;
t592 = t160 + t578;
t300 = qJD(4) * t306;
t252 = t299 * t300;
t254 = pkin(3) * t306 - qJ(4) * t308;
t255 = rSges(6,1) * t306 - rSges(6,2) * t308;
t548 = pkin(4) * t306;
t406 = t255 + t548;
t394 = -t254 - t406;
t435 = qJD(5) * t298;
t439 = qJD(3) * t299;
t328 = t394 * t439 - t435;
t590 = t252 + t328;
t276 = pkin(4) * t493;
t589 = t532 * t298 + t276;
t585 = -t602 * t299 + (-t238 * t298 + t130 - t138 + t513 - t623) * qJD(1);
t584 = 0.2e1 * qJD(3);
t583 = -rSges(5,1) - pkin(3);
t189 = t255 * t298;
t543 = rSges(5,1) * t306;
t256 = -rSges(5,3) * t308 + t543;
t581 = (qJD(3) * t256 - t300) * t298;
t259 = pkin(3) * t308 + t510;
t192 = t259 * t298;
t165 = qJD(1) * t192;
t293 = t299 * pkin(6);
t208 = pkin(2) * t298 - t293;
t206 = qJD(1) * t208;
t580 = -t165 - t206;
t404 = t299 * rSges(3,1) - rSges(3,2) * t298;
t577 = t304 + t404;
t570 = (-t306 * t597 + t308 * t596) * qJD(1);
t569 = -t306 * t595 + t308 * t594;
t568 = t642 * qJD(1);
t567 = -Icges(4,2) * t495 + t614 * t298 - t269 + t628;
t566 = t603 * t298 + t630;
t260 = rSges(6,1) * t308 + t541;
t559 = -rSges(6,3) * t299 - t260 * t298;
t437 = qJD(3) * t308;
t436 = qJD(4) * t308;
t200 = qJD(3) * t259 - t436;
t465 = -t260 * qJD(3) - t200;
t547 = pkin(4) * t308;
t555 = -pkin(4) * t437 - qJD(3) * (-t259 - t260 - t547) + t465;
t554 = t306 * t567 + t308 * t566;
t553 = m(5) / 0.2e1;
t552 = m(6) / 0.2e1;
t307 = sin(qJ(1));
t549 = pkin(1) * t307;
t545 = qJD(1) / 0.2e1;
t544 = rSges(4,1) * t308;
t538 = pkin(1) * qJD(1);
t257 = rSges(4,1) * t306 + rSges(4,2) * t308;
t196 = t257 * t299;
t440 = qJD(3) * t298;
t423 = t257 * t440;
t60 = qJD(1) * t592 - t423;
t537 = t196 * t60;
t452 = rSges(4,2) * t496 + t299 * rSges(4,3);
t157 = rSges(4,1) * t495 - t452;
t409 = -t208 - t549;
t420 = t257 * t439;
t59 = -t420 + (-t157 + t409) * qJD(1);
t536 = t298 * t59;
t535 = t299 * t59;
t534 = -rSges(6,2) - qJ(4);
t533 = -rSges(5,3) - qJ(4);
t203 = t209 * qJD(1);
t441 = qJD(1) * t306;
t168 = t298 * t437 + t299 * t441;
t438 = qJD(3) * t306;
t422 = t298 * t438;
t231 = pkin(3) * t422;
t76 = qJ(4) * t168 + qJD(1) * t277 + t298 * t300 - t231;
t531 = -t203 - t76;
t418 = t299 * t437;
t228 = rSges(6,2) * t418;
t419 = t299 * t438;
t443 = qJD(1) * t298;
t336 = -t308 * t443 - t419;
t489 = -rSges(6,1) * t419 + pkin(4) * t336 - qJ(5) * t442 + qJD(1) * t559 + t228 - t435;
t230 = pkin(4) * t422;
t434 = qJD(5) * t299;
t399 = -t230 + t434;
t488 = -qJD(3) * t189 + t399 + (t260 * t299 + t589) * qJD(1);
t470 = pkin(4) * t495 + qJ(5) * t299 - t559;
t451 = rSges(6,1) * t493 + rSges(6,2) * t494;
t469 = t451 + t589;
t468 = -t159 - t197;
t467 = t298 * t192 + t299 * t197;
t193 = t254 * t299;
t466 = -qJD(1) * t193 + t298 * t436;
t261 = rSges(5,1) * t308 + rSges(5,3) * t306;
t464 = -t261 * qJD(3) - t200;
t463 = qJ(4) * t418 + t252;
t462 = rSges(5,2) * t442 + rSges(5,3) * t418;
t424 = t298 * t441;
t461 = rSges(4,2) * t424 + rSges(4,3) * t442;
t454 = -t254 - t256;
t453 = -t259 - t261;
t450 = t298 ^ 2 + t299 ^ 2;
t433 = qJD(3) * qJD(4);
t431 = t311 * t549;
t430 = t311 * t304;
t75 = pkin(3) * t336 - qJ(4) * t424 + t463;
t429 = t192 * t442 + t298 * t76 + t299 * t75;
t188 = t254 * t298;
t428 = -t188 * t440 - t193 * t439 + t300;
t427 = -t197 - t469;
t280 = pkin(6) * t442;
t426 = t280 + t463;
t204 = t254 * t440;
t417 = -pkin(2) - t544;
t416 = t308 * t433;
t413 = -t440 / 0.2e1;
t410 = t439 / 0.2e1;
t407 = t293 - t549;
t403 = t299 * t454;
t401 = -t134 + t503;
t400 = t306 * t433 + t76 * t440 + (t165 + t75) * t439;
t397 = qJD(1) * (-pkin(2) * t443 + t280) - t431;
t396 = -t192 + t409;
t389 = t308 * t583 - pkin(2);
t207 = rSges(3,1) * t298 + rSges(3,2) * t299;
t385 = -rSges(4,2) * t306 + t544;
t384 = -t298 * t60 - t535;
t366 = t157 * t298 + t160 * t299;
t362 = qJD(1) * t204 + t299 * t416 - t430;
t357 = t192 * t440 + t197 * t439 + qJD(2) - t436;
t356 = qJD(3) * t403 + t252;
t191 = t257 * t298;
t190 = t256 * t298;
t337 = t298 * t416 + t397 + (t252 + t75) * qJD(1);
t335 = (-t547 * qJD(3) + t465) * qJD(3);
t103 = rSges(4,1) * t336 - rSges(4,2) * t418 + t461;
t106 = -qJD(3) * t191 + (t299 * t385 + t287) * qJD(1);
t327 = t103 * t299 + t106 * t298 + (t157 * t299 - t160 * t298) * qJD(1);
t290 = t299 * rSges(5,2);
t253 = t299 * t436;
t221 = t385 * qJD(3);
t205 = t254 * t443;
t195 = t256 * t299;
t194 = t255 * t299;
t169 = t418 - t424;
t167 = t450 * t438;
t156 = t261 * t298 - t290;
t105 = -qJD(3) * t190 + (t261 * t299 + t288) * qJD(1);
t102 = rSges(5,1) * t336 - rSges(5,3) * t424 + t462;
t58 = qJD(3) * t366 + qJD(2);
t45 = qJD(1) * t598 - t204 - t581;
t44 = (-t156 + t396) * qJD(1) + t356;
t43 = -t430 - t221 * t439 + (-t106 - t203 + t423) * qJD(1);
t42 = -t221 * t440 + (t103 - t420) * qJD(1) + t397;
t41 = (t156 * t298 + t159 * t299) * qJD(3) + t357;
t40 = -t204 + (-qJD(3) * t255 + t300) * t298 + (t591 + t469) * qJD(1) + t399;
t39 = (t396 - t470) * qJD(1) + t590;
t32 = (t298 * t470 + t299 * t469) * qJD(3) + t357;
t25 = t327 * qJD(3);
t24 = t464 * t439 + (-t105 + t531 + t581) * qJD(1) + t362;
t23 = qJD(1) * t102 + (qJD(1) * t403 + t298 * t464) * qJD(3) + t337;
t16 = t335 * t299 + (-t434 + (qJD(3) * t406 - t300) * t298 - t488 + t531) * qJD(1) + t362;
t15 = t335 * t298 + (t328 + t489) * qJD(1) + t337;
t2 = (t102 * t299 + t105 * t298 + (t156 * t299 + t298 * t468) * qJD(1)) * qJD(3) + t400;
t1 = (t489 * t299 + t488 * t298 + (t298 * t427 + t299 * t470) * qJD(1)) * qJD(3) + t400;
t3 = [m(3) * ((-t207 * t311 - t431) * t577 + (-t430 + (-0.2e1 * t404 - t304 + t577) * t311) * (-t207 - t549)) + (-qJD(3) * t650 + t643 * t306 + t644 * t308) * qJD(1) + (t16 * t407 + t39 * (-t309 * t538 + t230 + t231) + t15 * (t276 + t591 + t451) + (t16 * t532 + t600 * t39) * t299 + (-t16 * pkin(2) + t15 * t532 + (t16 * t534 + t39 * (rSges(6,1) * qJD(3) - qJD(4))) * t306 + (qJD(3) * t39 * t534 + t16 * t432) * t308 + t39 * (-pkin(6) - t532) * qJD(1)) * t298 + (-t307 * t538 + t228 + t426 + (qJD(1) * t532 + t432 * t438) * t299 + t600 * t298 + t39 - (-t470 - t549) * qJD(1) - t580 - t590) * t40) * m(6) + (-(-t44 + (-t156 - t549) * qJD(1) + t356 + t580) * t45 + t24 * (t290 + t407) + t44 * t231 + t23 * t598 + t45 * (t419 * t583 + t426 + t462) + ((-t307 * t45 - t309 * t44) * pkin(1) + t44 * (t306 * t533 + t389) * t299) * qJD(1) + (t24 * t389 + (-t44 * qJD(4) + t24 * t533) * t306 + t44 * (t308 * t533 + t543) * qJD(3) + (t44 * (-rSges(5,2) - pkin(6)) + t45 * (-pkin(2) + t453)) * qJD(1)) * t298) * m(5) + (t43 * (t298 * t417 + t407 + t452) + t42 * t592 + t60 * (t280 + t461) + (t257 * t536 - t537) * qJD(3) + ((-t307 * t60 - t309 * t59) * pkin(1) + (-pkin(2) - t385) * t535 + (t59 * (-rSges(4,3) - pkin(6)) + t60 * t417) * t298) * qJD(1) - (-t420 - t206 - t59 + (-t157 - t549) * qJD(1)) * t60) * m(4) + ((-t593 + (t51 - t113 + (t135 + t504) * t299 + t646) * t299 + (t298 * t621 - t582 + t601 + t616) * t298) * qJD(3) + t620) * t410 + (((t299 * t401 + t507 - t601 + t606) * t299 + (t298 * t401 - t126 + t128 + t402 + t508 + t607 - t615) * t298) * qJD(3) + t613 + t617) * t413 + (t609 + t610) * t440 / 0.2e1 - (t608 - t611 + t612) * t439 / 0.2e1 + ((t604 - t637) * t298 + (t605 + t636) * t299) * qJD(3) * t545; m(4) * t25 + m(5) * t2 + m(6) * t1; (t1 * t467 + (t1 * t470 + t15 * t394) * t298 + (t1 * t469 + t16 * t394) * t299 + (-t466 + t555 * t298 + (pkin(4) * t494 + t299 * t394 + t194) * qJD(1)) * t40 + (-qJD(1) * t188 + t299 * t555 + t205 - t253) * t39 + (-t428 - (-t189 * t298 - t194 * t299 - t450 * t548) * qJD(3) + t429 + (qJD(1) * t427 + t488) * t298 + (qJD(1) * t470 + t489) * t299) * t32) * m(6) + (-t44 * (t253 + (t188 + t190) * qJD(1)) - t45 * (-qJD(1) * t195 + t466) - t41 * t428 - ((-t41 * t195 + t44 * t453) * t299 + (-t41 * t190 + t45 * t453) * t298) * qJD(3) + t44 * t205 + t2 * t467 + t41 * t429 + (t24 * t454 + t44 * t464 + t2 * t159 + t41 * t102 + (t41 * t156 + t45 * t454) * qJD(1)) * t299 + (t23 * t454 + t45 * t464 + t2 * t156 + t41 * t105 + (t44 * t256 + t41 * t468) * qJD(1)) * t298) * m(5) + (-(t191 * t59 - t537) * qJD(1) - (t58 * (-t191 * t298 - t196 * t299) + t384 * t385) * qJD(3) + t25 * t366 + t58 * t327 + t384 * t221 + (-t42 * t298 - t43 * t299 + (-t299 * t60 + t536) * qJD(1)) * t257) * m(4) - (((t298 * t595 - t567 * t299) * t308 + (t298 * t594 + t566 * t299) * t306) * qJD(3) + (t596 * t306 + t597 * t308) * qJD(1)) * qJD(1) / 0.2e1 + (t611 * t299 + t610 * t298 + (t604 * t298 + t605 * t299) * qJD(1)) * t545 + ((-t571 * t440 - t568) * t298 + ((t554 * t299 + (t569 + t572) * t298) * qJD(3) + t570) * t299) * t413 + ((-t572 * t439 + t568) * t299 + ((t569 * t298 + (t554 + t571) * t299) * qJD(3) + t570) * t298) * t410 + (t609 * qJD(1) + (t606 * t442 + (t607 * qJD(1) + t585 * t298 - t599) * t298) * t584) * t298 / 0.2e1 - (t608 * qJD(1) + ((t648 * qJD(1) + t599) * t299 + (t649 * qJD(1) - t585 * t299) * t298) * t584) * t299 / 0.2e1 + (t613 + t618) * t443 / 0.2e1 + (t612 + t619) * t442 / 0.2e1; -m(5) * (t167 * t41 + t168 * t45 + t169 * t44) - m(6) * (t167 * t32 + t168 * t40 + t169 * t39) + 0.2e1 * ((t439 * t44 + t440 * t45 - t2) * t553 + (t39 * t439 + t40 * t440 - t1) * t552) * t308 + 0.2e1 * ((qJD(3) * t41 + t23 * t298 + t24 * t299 - t44 * t443 + t442 * t45) * t553 + (qJD(3) * t32 + t15 * t298 + t16 * t299 - t39 * t443 + t40 * t442) * t552) * t306; m(6) * (t15 * t299 - t16 * t298);];
tauc = t3(:);
