% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP9_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:37
% EndTime: 2019-12-31 18:48:53
% DurationCPUTime: 8.81s
% Computational Cost: add. (32360->457), mult. (26831->602), div. (0->0), fcn. (24664->7), ass. (0->292)
t419 = pkin(8) + qJ(3);
t405 = qJ(4) + t419;
t400 = cos(t405);
t399 = sin(t405);
t552 = Icges(5,4) * t399;
t358 = Icges(5,1) * t400 - t552;
t423 = sin(qJ(1));
t424 = cos(qJ(1));
t287 = Icges(5,5) * t423 + t358 * t424;
t393 = Icges(6,5) * t399;
t641 = Icges(6,1) * t400 + t393;
t683 = Icges(6,4) * t423 + t424 * t641 + t287;
t394 = Icges(5,4) * t400;
t357 = Icges(5,1) * t399 + t394;
t551 = Icges(6,5) * t400;
t682 = Icges(6,1) * t399 + t357 - t551;
t441 = Icges(6,3) * t400 - t393;
t676 = Icges(5,2) * t400 + t441 + t552;
t362 = pkin(4) * t400 + qJ(5) * t399;
t363 = rSges(6,1) * t400 + rSges(6,3) * t399;
t666 = t362 + t363;
t229 = t666 * t423;
t418 = t424 * rSges(6,2);
t421 = t424 ^ 2;
t567 = t423 * rSges(6,2);
t132 = t424 * (t363 * t424 + t567) + t421 * t362 + (-t418 + t229) * t423;
t401 = cos(pkin(8)) * pkin(2) + pkin(1);
t404 = cos(t419);
t576 = pkin(3) * t404;
t375 = t401 + t576;
t575 = -pkin(6) - qJ(2);
t478 = -pkin(7) + t575;
t396 = t423 * t478;
t402 = t423 * t575;
t471 = t424 * t575;
t485 = -t423 * t375 - t424 * t478;
t504 = -t423 * (t423 * t401 + t471 + t485) + t424 * (-t396 + t402 + (t375 - t401) * t424);
t102 = t132 + t504;
t534 = t399 * t423;
t599 = rSges(6,1) + pkin(4);
t482 = t599 * t534;
t527 = t400 * t424;
t557 = rSges(6,3) + qJ(5);
t483 = t557 * t527;
t528 = t400 * t423;
t533 = t399 * t424;
t145 = (-t599 * t533 + t483) * t424 + (t557 * t528 - t482) * t423;
t468 = t557 * t400;
t473 = t599 * t399;
t487 = t473 - t468;
t403 = sin(t419);
t577 = pkin(3) * t403;
t452 = t487 + t577;
t209 = t452 * t423;
t211 = t452 * t424;
t232 = t666 * t424;
t46 = t102 * t145 + t209 * t229 + t211 * t232;
t289 = rSges(5,1) * t528 - rSges(5,2) * t534 - t424 * rSges(5,3);
t464 = -rSges(5,2) * t533 + t423 * rSges(5,3);
t188 = t423 * t289 + t424 * (rSges(5,1) * t527 + t464);
t110 = t188 + t504;
t568 = rSges(5,1) * t400;
t364 = -rSges(5,2) * t399 + t568;
t361 = rSges(5,1) * t399 + rSges(5,2) * t400;
t431 = t361 + t577;
t646 = t431 * t424;
t647 = t431 * t423;
t639 = t423 * t647 + t424 * t646;
t337 = t361 * t423;
t339 = t361 * t424;
t642 = t423 * t337 + t424 * t339;
t68 = -t110 * t642 + t639 * t364;
t681 = -m(5) * t68 - m(6) * t46;
t626 = m(6) / 0.2e1;
t373 = rSges(4,1) * t403 + rSges(4,2) * t404;
t420 = t423 ^ 2;
t480 = t420 + t421;
t645 = t480 * t373;
t663 = -m(5) / 0.2e1;
t474 = (-t209 * t423 - t211 * t424) * t626 + t639 * t663 - m(4) * t645 / 0.2e1;
t198 = (-t468 + t577) * t423 + t482;
t199 = (-t473 - t577) * t424 + t483;
t347 = t373 * t423;
t348 = t373 * t424;
t217 = t423 * t347 + t348 * t424;
t627 = m(5) / 0.2e1;
t664 = m(4) / 0.2e1;
t475 = (t423 * t198 - t199 * t424) * t626 + t639 * t627 + t217 * t664;
t21 = t475 - t474;
t680 = t21 * qJD(1);
t228 = t487 * t423;
t231 = t487 * t424;
t669 = t361 * t480;
t509 = (-t228 * t423 - t231 * t424) * t626 + t663 * t669;
t221 = -t423 * t468 + t482;
t222 = -t424 * t473 + t483;
t510 = (t423 * t221 - t222 * t424) * t626 + t642 * t627;
t56 = t510 - t509;
t679 = t56 * qJD(1);
t379 = Icges(6,5) * t527;
t277 = Icges(6,6) * t423 + Icges(6,3) * t533 + t379;
t351 = Icges(5,5) * t400 - Icges(5,6) * t399;
t542 = t351 * t424;
t279 = Icges(5,3) * t423 + t542;
t352 = Icges(6,4) * t400 + Icges(6,6) * t399;
t541 = t352 * t424;
t678 = t277 * t533 + t683 * t527 + (Icges(6,2) * t423 + t279 + t541) * t423;
t354 = -Icges(5,2) * t399 + t394;
t677 = t354 + t682;
t675 = (Icges(5,6) - Icges(6,6)) * t400 + (Icges(6,4) + Icges(5,5)) * t399;
t674 = -t676 * t424 + t683;
t284 = -Icges(6,4) * t424 + t423 * t641;
t380 = Icges(5,4) * t534;
t286 = Icges(5,1) * t528 - Icges(5,5) * t424 - t380;
t673 = -Icges(5,2) * t528 - t441 * t423 + t284 + t286 - t380;
t283 = Icges(5,6) * t423 + t354 * t424;
t672 = -Icges(6,1) * t533 - t357 * t424 + t277 - t283 + t379;
t350 = Icges(6,3) * t399 + t551;
t276 = -Icges(6,6) * t424 + t350 * t423;
t282 = Icges(5,4) * t528 - Icges(5,2) * t534 - Icges(5,6) * t424;
t671 = t682 * t423 - t276 + t282;
t670 = t641 + t358;
t668 = -t283 * t533 + t678;
t648 = (t276 * t399 + t284 * t400) * t423;
t665 = t648 + t678;
t601 = t423 / 0.2e1;
t600 = -t424 / 0.2e1;
t659 = t424 / 0.2e1;
t520 = t423 * t352;
t258 = t423 * (-Icges(6,2) * t424 + t520);
t141 = t276 * t533 + t284 * t527 + t258;
t658 = t141 * t424;
t553 = Icges(4,4) * t403;
t369 = Icges(4,2) * t404 + t553;
t372 = Icges(4,1) * t404 - t553;
t657 = (t372 / 0.2e1 - t369 / 0.2e1) * t403;
t656 = t675 * t423;
t655 = t675 * t424;
t654 = (t670 - t676) * t400 + (t350 - t677) * t399;
t213 = -t289 + t485;
t214 = -t396 + (t375 + t568) * t424 + t464;
t644 = t213 * t424 + t214 * t423;
t433 = t644 * t364;
t635 = -t557 * t399 - t599 * t400;
t180 = t635 * t423 + t418 + t485;
t181 = t567 - t396 + (t375 - t635) * t424;
t511 = -t232 * t180 - t229 * t181;
t573 = (-t198 * t231 - t199 * t228 + t511) * t626 + (-t433 + (t423 * t646 - t424 * t647) * t361) * t627;
t574 = (-t209 * t222 - t211 * t221 + t511) * t626 + (-t337 * t646 + t339 * t647 - t433) * t627;
t653 = t573 - t574;
t651 = (t672 * t423 + t671 * t424) * t400 + (-t674 * t423 + t673 * t424) * t399;
t569 = rSges(4,1) * t404;
t466 = t401 + t569;
t526 = t403 * t423;
t481 = rSges(4,2) * t526 + t424 * rSges(4,3);
t234 = -t423 * t466 - t471 + t481;
t525 = t403 * t424;
t465 = -rSges(4,2) * t525 + t423 * rSges(4,3);
t235 = t424 * t466 - t402 + t465;
t643 = t234 * t424 + t235 * t423;
t398 = Icges(4,4) * t404;
t370 = -Icges(4,2) * t403 + t398;
t371 = Icges(4,1) * t403 + t398;
t106 = -t188 * t642 + t364 * t669;
t59 = t132 * t145 + t228 * t229 + t231 * t232;
t638 = -m(5) * t106 - m(6) * t59;
t323 = Icges(4,5) * t423 + t372 * t424;
t488 = -t369 * t424 + t323;
t390 = Icges(4,4) * t526;
t524 = t404 * t423;
t322 = Icges(4,1) * t524 - Icges(4,5) * t424 - t390;
t489 = -Icges(4,2) * t524 + t322 - t390;
t321 = Icges(4,6) * t423 + t370 * t424;
t490 = -t371 * t424 - t321;
t320 = Icges(4,4) * t524 - Icges(4,2) * t526 - Icges(4,6) * t424;
t491 = t371 * t423 + t320;
t634 = (-t488 * t423 + t424 * t489) * t403 + (t490 * t423 + t424 * t491) * t404;
t430 = (-t350 / 0.2e1 + t677 / 0.2e1) * t400 + (-t676 / 0.2e1 + t670 / 0.2e1) * t399;
t241 = t287 * t528;
t457 = t424 * t279 - t241;
t140 = -t283 * t534 - t457;
t278 = Icges(5,5) * t528 - Icges(5,6) * t534 - Icges(5,3) * t424;
t503 = -t423 * t278 - t286 * t527;
t143 = -t282 * t533 - t503;
t455 = t283 * t399 - t278;
t544 = t282 * t399;
t432 = (-t143 * t424 + t668 * t423 - t658) * t659 + (-t658 + (t140 - t241 + (t279 + t544) * t424 + t503) * t424 + (-t648 + t665) * t423) * t600 + (((t455 + t278) * t424 - t665 + t668) * t424 + (t141 - t258 + t143 + t457 - (t286 * t400 - t544) * t424 + t140 + t423 * t455) * t423) * t601;
t631 = 0.4e1 * qJD(1);
t630 = 2 * qJD(3);
t628 = 2 * qJD(4);
t346 = t480 * t399;
t507 = -t209 * t528 - t211 * t527;
t64 = t102 * t346 + t507;
t506 = -t228 * t528 - t231 * t527;
t82 = t132 * t346 + t506;
t618 = m(6) * (t82 + t64);
t450 = -t145 * t400 - t229 * t534 - t232 * t533;
t476 = t399 * t102 + t507;
t617 = m(6) * (t450 + t476);
t48 = t132 * t399 + t450 + t506;
t615 = m(6) * t48;
t508 = t180 * t527 + t181 * t528;
t610 = m(6) * ((t198 * t424 + t199 * t423) * t399 + t508);
t609 = m(6) * (-t209 * t533 + t211 * t534 + t508);
t608 = m(6) * ((t221 * t424 + t222 * t423) * t399 + t508);
t607 = m(6) * (-t228 * t533 + t231 * t534 + t508);
t606 = m(6) * (t180 * t198 + t181 * t199);
t605 = m(6) * (t180 * t221 + t181 * t222);
t602 = -t423 / 0.2e1;
t598 = m(3) * t480 * (rSges(3,3) + qJ(2));
t597 = m(4) * (t234 * t347 - t235 * t348);
t596 = m(4) * t643;
t592 = m(5) * (t213 * t647 - t214 * t646);
t591 = m(5) * (t213 * t337 - t214 * t339);
t590 = m(5) * t644;
t584 = m(6) * (t180 * t424 + t181 * t423);
t572 = m(6) * qJD(3);
t571 = m(6) * qJD(4);
t570 = m(6) * qJD(5);
t543 = t320 * t403;
t535 = t399 * t400;
t523 = t404 * t424;
t318 = Icges(4,5) * t524 - Icges(4,6) * t526 - Icges(4,3) * t424;
t501 = -t423 * t318 - t322 * t523;
t444 = Icges(4,5) * t404 - Icges(4,6) * t403;
t319 = Icges(4,3) * t423 + t424 * t444;
t500 = t423 * t319 + t323 * t523;
t484 = t480 * t535;
t215 = m(6) * t346;
t479 = t215 * qJD(1);
t108 = -t180 * t534 + t181 * t533;
t477 = m(6) * t108 * qJD(1);
t469 = -t364 - t576;
t467 = ((t656 * t423 + t651) * t424 - t655 * t420) * t601 + ((t655 * t424 + t651) * t423 - t656 * t421) * t600;
t250 = t323 * t524;
t456 = t424 * t319 - t250;
t454 = t321 * t403 - t318;
t453 = t480 * t577;
t451 = -t666 - t576;
t443 = -Icges(4,5) * t403 - Icges(4,6) * t404;
t426 = -t432 + (t423 * t351 + t672 * t399 + t674 * t400 + t654 * t424 + t520) * t601 + (-t671 * t399 + t673 * t400 + t654 * t423 - t541 - t542) * t600;
t425 = -t430 + ((t276 + t282) * t400 + (-t284 + t286) * t399) * (t601 + t602);
t374 = -rSges(4,2) * t403 + t569;
t341 = t443 * t424;
t340 = t443 * t423;
t263 = t469 * t424;
t261 = t469 * t423;
t219 = t484 - t535;
t212 = t451 * t424;
t210 = t451 * t423;
t200 = t219 * t570;
t182 = -t453 - t642;
t159 = t229 * t424 - t232 * t423;
t151 = -t321 * t525 + t500;
t150 = -t320 * t525 - t501;
t149 = -t321 * t526 - t456;
t147 = t159 * t571;
t125 = (-t346 * t400 - t219 + t484) * t570;
t122 = -t453 + t145;
t105 = -t150 * t424 + t151 * t423;
t104 = -(-t423 * (-t322 * t404 + t543) - t424 * t318) * t424 + t149 * t423;
t75 = t607 / 0.2e1;
t73 = t608 / 0.2e1;
t71 = t609 / 0.2e1;
t69 = t610 / 0.2e1;
t57 = t509 + t510;
t49 = t584 + t590 + t596 + t598;
t47 = t615 / 0.2e1;
t42 = t617 / 0.2e1;
t40 = t618 / 0.2e1;
t39 = (t149 - t250 + (t319 + t543) * t424 + t501) * t424 + t500 * t423;
t38 = (t424 * t454 + t151 - t500) * t424 + (t423 * t454 + t150 + t456) * t423;
t37 = t430 + t591 + t605;
t22 = t474 + t475;
t20 = t75 - t608 / 0.2e1;
t19 = t75 + t73;
t18 = t73 - t607 / 0.2e1;
t17 = (t371 / 0.2e1 + t370 / 0.2e1) * t404 + t657 + t597 + t592 + t606 + t430;
t16 = t71 - t610 / 0.2e1;
t15 = t71 + t69;
t14 = t69 - t609 / 0.2e1;
t11 = t40 + t47 - t617 / 0.2e1;
t10 = t40 + t42 - t615 / 0.2e1;
t9 = t42 + t47 - t618 / 0.2e1;
t8 = t467 - t638;
t7 = t8 * qJD(4);
t6 = t467 - t681;
t4 = (-t39 / 0.2e1 + t105 / 0.2e1) * t424 + (t104 / 0.2e1 + t38 / 0.2e1) * t423 + t432;
t3 = t432 + t653;
t2 = t432 - t653;
t1 = t426 + t573 + t574;
t5 = [t49 * qJD(2) + t17 * qJD(3) + t37 * qJD(4) + t108 * t570, qJD(1) * t49 + qJD(3) * t22 + qJD(4) * t57, t17 * qJD(1) + t22 * qJD(2) + t1 * qJD(4) + t15 * qJD(5) + ((-t643 * t374 + (-t347 * t424 + t348 * t423) * t373) * t664 + (t213 * t263 + t214 * t261) * t627 + (t180 * t212 + t181 * t210 - t198 * t211 - t199 * t209) * t626) * t630 + (t426 + (t403 * t490 + t404 * t488) * t601 + t39 * t659 + (t104 + t38) * t602 + (-t403 * t491 + t489 * t404 + t105) * t600 + (t420 / 0.2e1 + t421 / 0.2e1) * t444) * qJD(3), t37 * qJD(1) + t57 * qJD(2) + t1 * qJD(3) + t426 * qJD(4) + t19 * qJD(5) + ((-t433 + (-t337 * t424 + t339 * t423) * t361) * t627 + (-t221 * t231 - t222 * t228 + t511) * t626) * t628, t15 * qJD(3) + t19 * qJD(4) + t477; t21 * qJD(3) + t56 * qJD(4) - t215 * qJD(5) + (-t590 / 0.4e1 - t584 / 0.4e1 - t596 / 0.4e1 - t598 / 0.4e1) * t631, 0, t680 + ((-t261 * t424 + t263 * t423) * t627 + (-t210 * t424 + t212 * t423) * t626) * t630 + t147, t159 * t572 + t147 + t679, -t479; (t425 - (t371 + t370) * t404 / 0.2e1 - t657) * qJD(1) - t21 * qJD(2) + t4 * qJD(3) + t2 * qJD(4) + t16 * qJD(5) + (-t597 / 0.4e1 - t592 / 0.4e1 - t606 / 0.4e1) * t631, -t680, t4 * qJD(1) + (m(6) * (t102 * t122 - t209 * t210 - t211 * t212) + m(5) * (t110 * t182 - t261 * t647 - t263 * t646) + m(4) * (t374 * t645 - (t423 * (rSges(4,1) * t524 - t481) + t424 * (rSges(4,1) * t523 + t465)) * t217) + (t421 * t340 + (-t424 * t341 + t634) * t423) * t600 + (t420 * t341 + (-t423 * t340 + t634) * t424) * t601 + t467) * qJD(3) + t6 * qJD(4) + t64 * t570, t2 * qJD(1) + t6 * qJD(3) + t10 * qJD(5) + ((t59 + t46) * t626 + (t68 + t106) * t627) * t628 + (t467 + t638) * qJD(4), t16 * qJD(1) + t10 * qJD(4) + t64 * t572 + t125; t425 * qJD(1) - t56 * qJD(2) + t3 * qJD(3) + t432 * qJD(4) + t20 * qJD(5) + (-t591 / 0.4e1 - t605 / 0.4e1) * t631, -t679, t3 * qJD(1) + t7 + t11 * qJD(5) + ((t122 * t132 - t210 * t228 - t212 * t231 + t46) * t626 + (t188 * t182 + (-t261 * t423 - t263 * t424) * t361 + t68) * t627) * t630 + (t467 + t681) * qJD(3), qJD(1) * t432 + t8 * qJD(3) + t82 * t570 + t7, t20 * qJD(1) + t11 * qJD(3) + t82 * t571 + t125; t215 * qJD(2) + t14 * qJD(3) + t18 * qJD(4) - t477, t479, t14 * qJD(1) + (-t400 * t122 + (t210 * t423 + t212 * t424) * t399 - t64 + t476) * t572 + t9 * qJD(4) + t200, t18 * qJD(1) + t9 * qJD(3) + (t48 - t82) * t571 + t200, 0.4e1 * (qJD(3) / 0.4e1 + qJD(4) / 0.4e1) * t219 * m(6);];
Cq = t5;
