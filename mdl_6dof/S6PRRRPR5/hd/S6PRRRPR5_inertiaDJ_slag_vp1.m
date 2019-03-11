% Calculate time derivative of joint inertia matrix for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:23:07
% EndTime: 2019-03-08 23:24:02
% DurationCPUTime: 30.54s
% Computational Cost: add. (239225->1464), mult. (627548->2008), div. (0->0), fcn. (749137->16), ass. (0->577)
t590 = m(6) / 0.2e1 + m(7) / 0.2e1;
t674 = 0.2e1 * t590;
t489 = sin(pkin(7));
t496 = sin(qJ(3));
t488 = sin(pkin(12));
t491 = cos(pkin(12));
t497 = sin(qJ(2));
t492 = cos(pkin(6));
t500 = cos(qJ(2));
t626 = t492 * t500;
t520 = -t488 * t497 + t491 * t626;
t643 = cos(pkin(7));
t510 = t520 * t643;
t627 = t492 * t497;
t479 = t488 * t500 + t491 * t627;
t647 = cos(qJ(3));
t576 = t479 * t647;
t490 = sin(pkin(6));
t629 = t490 * t491;
t425 = t576 + (-t489 * t629 + t510) * t496;
t572 = t490 * t643;
t458 = -t489 * t520 - t491 * t572;
t592 = qJ(4) + pkin(13);
t487 = sin(t592);
t563 = cos(t592);
t403 = t425 * t563 + t458 * t487;
t575 = t489 * t647;
t550 = t490 * t575;
t424 = t479 * t496 + t491 * t550 - t510 * t647;
t494 = sin(qJ(6));
t498 = cos(qJ(6));
t353 = -t403 * t494 + t424 * t498;
t354 = t403 * t498 + t424 * t494;
t514 = -t425 * t487 + t458 * t563;
t251 = Icges(7,5) * t354 + Icges(7,6) * t353 - Icges(7,3) * t514;
t253 = Icges(7,4) * t354 + Icges(7,2) * t353 - Icges(7,6) * t514;
t255 = Icges(7,1) * t354 + Icges(7,4) * t353 - Icges(7,5) * t514;
t129 = -t251 * t514 + t253 * t353 + t255 * t354;
t519 = t488 * t626 + t491 * t497;
t509 = t519 * t643;
t518 = t488 * t627 - t491 * t500;
t633 = t489 * t490;
t427 = -t518 * t647 + (t488 * t633 - t509) * t496;
t459 = t488 * t572 + t489 * t519;
t405 = t427 * t563 + t459 * t487;
t426 = -t488 * t550 - t496 * t518 + t509 * t647;
t355 = -t405 * t494 + t426 * t498;
t356 = t405 * t498 + t426 * t494;
t513 = -t427 * t487 + t459 * t563;
t252 = Icges(7,5) * t356 + Icges(7,6) * t355 - Icges(7,3) * t513;
t254 = Icges(7,4) * t356 + Icges(7,2) * t355 - Icges(7,6) * t513;
t256 = Icges(7,1) * t356 + Icges(7,4) * t355 - Icges(7,5) * t513;
t130 = -t252 * t514 + t254 * t353 + t256 * t354;
t571 = t496 * t643;
t511 = t497 * t647 + t500 * t571;
t631 = t489 * t496;
t457 = t490 * t511 + t492 * t631;
t478 = t492 * t643 - t500 * t633;
t418 = t457 * t563 + t478 * t487;
t549 = t492 * t575;
t533 = t643 * t647;
t670 = -t496 * t497 + t500 * t533;
t456 = -t490 * t670 - t549;
t400 = -t418 * t494 + t456 * t498;
t401 = t418 * t498 + t456 * t494;
t512 = -t457 * t487 + t478 * t563;
t301 = Icges(7,5) * t401 + Icges(7,6) * t400 - Icges(7,3) * t512;
t302 = Icges(7,4) * t401 + Icges(7,2) * t400 - Icges(7,6) * t512;
t303 = Icges(7,1) * t401 + Icges(7,4) * t400 - Icges(7,5) * t512;
t156 = -t301 * t514 + t302 * t353 + t303 * t354;
t47 = -t129 * t514 - t130 * t513 - t156 * t512;
t673 = t47 / 0.2e1;
t672 = t490 / 0.2e1;
t648 = t492 / 0.2e1;
t470 = t520 * qJD(2);
t471 = t479 * qJD(2);
t573 = qJD(3) * t631;
t396 = t470 * t496 + t471 * t533 - t573 * t629 + (t496 * t510 + t576) * qJD(3);
t397 = -qJD(3) * t424 + t470 * t647 - t471 * t571;
t495 = sin(qJ(4));
t499 = cos(qJ(4));
t407 = -t425 * t495 + t458 * t499;
t632 = t489 * t495;
t504 = qJD(4) * t407 + t471 * t632;
t645 = pkin(4) * t499;
t213 = pkin(4) * t504 + qJ(5) * t396 + qJD(5) * t424 + t397 * t645;
t640 = t458 * t495;
t298 = pkin(4) * t640 + qJ(5) * t424 + t425 * t645;
t472 = t519 * qJD(2);
t473 = t518 * qJD(2);
t398 = qJD(3) * t427 - t472 * t496 - t473 * t533;
t621 = t426 * t213 + t398 * t298;
t638 = t471 * t489;
t442 = pkin(2) * t470 + pkin(9) * t638;
t637 = t473 * t489;
t443 = -pkin(2) * t472 - pkin(9) * t637;
t634 = t488 * t490;
t596 = t442 * t634 + t443 * t629;
t399 = -qJD(3) * t426 - t472 * t647 + t473 * t571;
t331 = rSges(4,1) * t399 - rSges(4,2) * t398 - rSges(4,3) * t637;
t528 = t489 * t563;
t334 = qJD(4) * t405 + t399 * t487 + t473 * t528;
t635 = t487 * t489;
t335 = qJD(4) * t513 + t399 * t563 - t473 * t635;
t225 = rSges(6,1) * t335 - rSges(6,2) * t334 + rSges(6,3) * t398;
t639 = t459 * t495;
t410 = t427 * t499 + t639;
t630 = t489 * t499;
t339 = -qJD(4) * t410 - t399 * t495 - t473 * t630;
t409 = -t427 * t495 + t459 * t499;
t503 = qJD(4) * t409 - t473 * t632;
t340 = t399 * t499 + t503;
t238 = rSges(5,1) * t340 + rSges(5,2) * t339 + rSges(5,3) * t398;
t243 = -qJD(6) * t356 - t335 * t494 + t398 * t498;
t244 = qJD(6) * t355 + t335 * t498 + t398 * t494;
t154 = rSges(7,1) * t244 + rSges(7,2) * t243 + rSges(7,3) * t334;
t623 = pkin(5) * t335 + pkin(11) * t334 + t154;
t505 = -m(5) * t238 - m(6) * t225 - m(7) * t623;
t671 = -m(4) * t331 + t505;
t312 = rSges(6,1) * t405 + rSges(6,2) * t513 + rSges(6,3) * t426;
t323 = rSges(5,1) * t410 + rSges(5,2) * t409 + rSges(5,3) * t426;
t258 = rSges(7,1) * t356 + rSges(7,2) * t355 - rSges(7,3) * t513;
t615 = pkin(5) * t405 - pkin(11) * t513 + t258;
t669 = -m(5) * t323 - m(6) * t312 - m(7) * t615;
t311 = rSges(6,1) * t403 + rSges(6,2) * t514 + rSges(6,3) * t424;
t408 = t425 * t499 + t640;
t322 = rSges(5,1) * t408 + rSges(5,2) * t407 + rSges(5,3) * t424;
t257 = rSges(7,1) * t354 + rSges(7,2) * t353 - rSges(7,3) * t514;
t616 = pkin(5) * t403 - pkin(11) * t514 + t257;
t668 = m(5) * t322 + m(6) * t311 + m(7) * t616;
t420 = t492 * t573 + (t511 * qJD(3) + (t496 * t500 + t497 * t533) * qJD(2)) * t490;
t667 = -0.2e1 * t458;
t666 = m(5) / 0.2e1;
t333 = qJD(4) * t514 + t397 * t563 + t471 * t635;
t241 = -qJD(6) * t354 - t333 * t494 + t396 * t498;
t242 = qJD(6) * t353 + t333 * t498 + t396 * t494;
t332 = qJD(4) * t403 + t397 * t487 - t471 * t528;
t147 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t332;
t149 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t332;
t151 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t332;
t33 = -t147 * t514 + t149 * t353 + t151 * t354 + t241 * t253 + t242 * t255 + t251 * t332;
t148 = Icges(7,5) * t244 + Icges(7,6) * t243 + Icges(7,3) * t334;
t150 = Icges(7,4) * t244 + Icges(7,2) * t243 + Icges(7,6) * t334;
t152 = Icges(7,1) * t244 + Icges(7,4) * t243 + Icges(7,5) * t334;
t34 = -t148 * t514 + t150 * t353 + t152 * t354 + t241 * t254 + t242 * t256 + t252 * t332;
t421 = qJD(3) * t549 + (t670 * qJD(3) + (-t497 * t571 + t500 * t647) * qJD(2)) * t490;
t595 = qJD(2) * t490;
t574 = t497 * t595;
t376 = qJD(4) * t418 + t421 * t487 - t528 * t574;
t548 = t489 * t574;
t377 = qJD(4) * t512 + t421 * t563 + t487 * t548;
t283 = -qJD(6) * t401 - t377 * t494 + t420 * t498;
t284 = qJD(6) * t400 + t377 * t498 + t420 * t494;
t186 = Icges(7,5) * t284 + Icges(7,6) * t283 + Icges(7,3) * t376;
t187 = Icges(7,4) * t284 + Icges(7,2) * t283 + Icges(7,6) * t376;
t188 = Icges(7,1) * t284 + Icges(7,4) * t283 + Icges(7,5) * t376;
t49 = -t186 * t514 + t187 * t353 + t188 * t354 + t241 * t302 + t242 * t303 + t301 * t332;
t1 = t129 * t332 + t130 * t334 + t156 * t376 - t33 * t514 - t34 * t513 - t49 * t512;
t663 = t1 / 0.2e1;
t131 = -t251 * t513 + t253 * t355 + t255 * t356;
t132 = -t252 * t513 + t254 * t355 + t256 * t356;
t157 = -t301 * t513 + t302 * t355 + t303 * t356;
t35 = -t147 * t513 + t149 * t355 + t151 * t356 + t243 * t253 + t244 * t255 + t251 * t334;
t36 = -t148 * t513 + t150 * t355 + t152 * t356 + t243 * t254 + t244 * t256 + t252 * t334;
t50 = -t186 * t513 + t187 * t355 + t188 * t356 + t243 * t302 + t244 * t303 + t301 * t334;
t2 = t131 * t332 + t132 * t334 + t157 * t376 - t35 * t514 - t36 * t513 - t50 * t512;
t662 = t2 / 0.2e1;
t139 = -t251 * t512 + t253 * t400 + t255 * t401;
t140 = -t252 * t512 + t254 * t400 + t256 * t401;
t163 = -t301 * t512 + t302 * t400 + t303 * t401;
t37 = -t147 * t512 + t149 * t400 + t151 * t401 + t251 * t376 + t253 * t283 + t255 * t284;
t38 = -t148 * t512 + t150 * t400 + t152 * t401 + t252 * t376 + t254 * t283 + t256 * t284;
t67 = -t186 * t512 + t187 * t400 + t188 * t401 + t283 * t302 + t284 * t303 + t301 * t376;
t7 = t139 * t332 + t140 * t334 + t163 * t376 - t37 * t514 - t38 * t513 - t512 * t67;
t661 = t7 / 0.2e1;
t58 = -t139 * t514 - t140 * t513 - t163 * t512;
t660 = t58 / 0.2e1;
t659 = t332 / 0.2e1;
t658 = t334 / 0.2e1;
t657 = t376 / 0.2e1;
t656 = -t514 / 0.2e1;
t655 = -t513 / 0.2e1;
t654 = -t512 / 0.2e1;
t653 = -t488 / 0.2e1;
t652 = t488 / 0.2e1;
t651 = -t491 / 0.2e1;
t650 = t491 / 0.2e1;
t649 = -t492 / 0.2e1;
t642 = Icges(3,4) * t497;
t641 = Icges(3,4) * t500;
t636 = t478 * t495;
t446 = Icges(3,5) * t470 - Icges(3,6) * t471;
t628 = t491 * t446;
t153 = rSges(7,1) * t242 + rSges(7,2) * t241 + rSges(7,3) * t332;
t624 = pkin(5) * t333 + pkin(11) * t332 + t153;
t189 = rSges(7,1) * t284 + rSges(7,2) * t283 + rSges(7,3) * t376;
t622 = pkin(5) * t377 + pkin(11) * t376 + t189;
t214 = pkin(4) * t503 + qJ(5) * t398 + qJD(5) * t426 + t399 * t645;
t299 = pkin(4) * t639 + qJ(5) * t426 + t427 * t645;
t620 = t456 * t214 + t420 * t299;
t224 = rSges(6,1) * t333 - rSges(6,2) * t332 + rSges(6,3) * t396;
t619 = -t213 - t224;
t618 = -t214 - t225;
t337 = -qJD(4) * t408 - t397 * t495 + t471 * t630;
t338 = t397 * t499 + t504;
t237 = rSges(5,1) * t338 + rSges(5,2) * t337 + rSges(5,3) * t396;
t346 = pkin(3) * t397 + pkin(10) * t396;
t617 = -t237 - t346;
t428 = -t457 * t495 + t478 * t499;
t501 = qJD(4) * t428 + t495 * t548;
t269 = pkin(4) * t501 + qJ(5) * t420 + qJD(5) * t456 + t421 * t645;
t357 = pkin(4) * t636 + qJ(5) * t456 + t457 * t645;
t614 = t424 * t269 + t396 * t357;
t390 = pkin(3) * t421 + pkin(10) * t420;
t362 = t458 * t390;
t613 = t458 * t269 + t362;
t273 = rSges(6,1) * t377 - rSges(6,2) * t376 + rSges(6,3) * t420;
t612 = -t269 - t273;
t429 = t457 * t499 + t636;
t385 = -qJD(4) * t429 - t421 * t495 + t499 * t548;
t386 = t421 * t499 + t501;
t278 = rSges(5,1) * t386 + rSges(5,2) * t385 + rSges(5,3) * t420;
t611 = -t278 - t390;
t394 = pkin(3) * t425 + pkin(10) * t424;
t375 = t459 * t394;
t610 = t459 * t298 + t375;
t395 = pkin(3) * t427 + pkin(10) * t426;
t387 = t478 * t395;
t609 = t478 * t299 + t387;
t608 = -t298 - t311;
t607 = -t299 - t312;
t304 = rSges(7,1) * t401 + rSges(7,2) * t400 - rSges(7,3) * t512;
t606 = pkin(5) * t418 - pkin(11) * t512 + t304;
t315 = t459 * t346;
t380 = t395 * t638;
t605 = t315 - t380;
t604 = -t322 - t394;
t347 = pkin(3) * t399 + pkin(10) * t398;
t437 = t492 * t443;
t603 = t492 * t347 + t437;
t416 = pkin(3) * t457 + pkin(10) * t456;
t406 = t458 * t416;
t602 = t458 * t357 + t406;
t361 = rSges(6,1) * t418 + rSges(6,2) * t512 + rSges(6,3) * t456;
t601 = -t357 - t361;
t374 = rSges(5,1) * t429 + rSges(5,2) * t428 + rSges(5,3) * t456;
t600 = -t374 - t416;
t436 = -pkin(2) * t518 + pkin(9) * t459;
t430 = t492 * t436;
t599 = t492 * t395 + t430;
t435 = t479 * pkin(2) + pkin(9) * t458;
t598 = t435 * t634 + t436 * t629;
t597 = 0.2e1 * t596;
t591 = -0.2e1 * t637;
t589 = -t213 - t624;
t588 = -t214 - t623;
t587 = -t269 - t622;
t586 = t492 * t214 + t603;
t585 = -t346 + t619;
t584 = -t298 - t616;
t583 = -t299 - t615;
t582 = -t390 + t612;
t581 = t492 * t299 + t599;
t580 = -t394 + t608;
t579 = -t357 - t606;
t578 = t478 * t347 + t395 * t548 + t416 * t637;
t577 = -t416 + t601;
t570 = 0.2e1 * m(4);
t568 = 0.2e1 * m(5);
t566 = 0.2e1 * m(6);
t564 = 0.2e1 * m(7);
t384 = rSges(4,1) * t421 - rSges(4,2) * t420 + rSges(4,3) * t548;
t466 = (pkin(9) * t489 * t497 + pkin(2) * t500) * t595;
t562 = t490 * (-t384 - t466);
t415 = rSges(4,1) * t457 - rSges(4,2) * t456 + rSges(4,3) * t478;
t460 = t490 * t497 * pkin(2) + pkin(9) * t478;
t561 = (-t415 - t460) * t490;
t560 = -t346 + t589;
t559 = -t390 + t587;
t206 = t459 * t213;
t289 = t299 * t638;
t557 = t206 - t289 + t605;
t556 = -t394 + t584;
t555 = -t416 + t579;
t554 = t347 * t667 + t394 * t591 + 0.2e1 * t315 - 0.2e1 * t380;
t343 = t346 * t634;
t344 = t347 * t629;
t553 = 0.2e1 * t343 + 0.2e1 * t344 + t597;
t552 = t343 + t344 + t596;
t551 = t394 * t634 + t395 * t629 + t598;
t270 = Icges(6,5) * t377 - Icges(6,6) * t376 + Icges(6,3) * t420;
t271 = Icges(6,4) * t377 - Icges(6,2) * t376 + Icges(6,6) * t420;
t272 = Icges(6,1) * t377 - Icges(6,4) * t376 + Icges(6,5) * t420;
t358 = Icges(6,5) * t418 + Icges(6,6) * t512 + Icges(6,3) * t456;
t359 = Icges(6,4) * t418 + Icges(6,2) * t512 + Icges(6,6) * t456;
t360 = Icges(6,1) * t418 + Icges(6,4) * t512 + Icges(6,5) * t456;
t110 = t270 * t424 + t271 * t514 + t272 * t403 - t332 * t359 + t333 * t360 + t358 * t396;
t305 = Icges(6,5) * t403 + Icges(6,6) * t514 + Icges(6,3) * t424;
t307 = Icges(6,4) * t403 + Icges(6,2) * t514 + Icges(6,6) * t424;
t309 = Icges(6,1) * t403 + Icges(6,4) * t514 + Icges(6,5) * t424;
t166 = t305 * t424 + t307 * t514 + t309 * t403;
t306 = Icges(6,5) * t405 + Icges(6,6) * t513 + Icges(6,3) * t426;
t308 = Icges(6,4) * t405 + Icges(6,2) * t513 + Icges(6,6) * t426;
t310 = Icges(6,1) * t405 + Icges(6,4) * t513 + Icges(6,5) * t426;
t167 = t306 * t424 + t308 * t514 + t310 * t403;
t190 = t358 * t424 + t359 * t514 + t360 * t403;
t218 = Icges(6,5) * t333 - Icges(6,6) * t332 + Icges(6,3) * t396;
t220 = Icges(6,4) * t333 - Icges(6,2) * t332 + Icges(6,6) * t396;
t222 = Icges(6,1) * t333 - Icges(6,4) * t332 + Icges(6,5) * t396;
t72 = t218 * t424 + t220 * t514 + t222 * t403 + t305 * t396 - t307 * t332 + t309 * t333;
t219 = Icges(6,5) * t335 - Icges(6,6) * t334 + Icges(6,3) * t398;
t221 = Icges(6,4) * t335 - Icges(6,2) * t334 + Icges(6,6) * t398;
t223 = Icges(6,1) * t335 - Icges(6,4) * t334 + Icges(6,5) * t398;
t73 = t219 * t424 + t221 * t514 + t223 * t403 + t306 * t396 - t308 * t332 + t310 * t333;
t13 = t110 * t456 + t166 * t396 + t167 * t398 + t190 * t420 + t424 * t72 + t426 * t73;
t275 = Icges(5,5) * t386 + Icges(5,6) * t385 + Icges(5,3) * t420;
t276 = Icges(5,4) * t386 + Icges(5,2) * t385 + Icges(5,6) * t420;
t277 = Icges(5,1) * t386 + Icges(5,4) * t385 + Icges(5,5) * t420;
t369 = Icges(5,5) * t429 + Icges(5,6) * t428 + Icges(5,3) * t456;
t370 = Icges(5,4) * t429 + Icges(5,2) * t428 + Icges(5,6) * t456;
t371 = Icges(5,1) * t429 + Icges(5,4) * t428 + Icges(5,5) * t456;
t113 = t275 * t424 + t276 * t407 + t277 * t408 + t337 * t370 + t338 * t371 + t369 * t396;
t316 = Icges(5,5) * t408 + Icges(5,6) * t407 + Icges(5,3) * t424;
t318 = Icges(5,4) * t408 + Icges(5,2) * t407 + Icges(5,6) * t424;
t320 = Icges(5,1) * t408 + Icges(5,4) * t407 + Icges(5,5) * t424;
t172 = t316 * t424 + t318 * t407 + t320 * t408;
t317 = Icges(5,5) * t410 + Icges(5,6) * t409 + Icges(5,3) * t426;
t319 = Icges(5,4) * t410 + Icges(5,2) * t409 + Icges(5,6) * t426;
t321 = Icges(5,1) * t410 + Icges(5,4) * t409 + Icges(5,5) * t426;
t173 = t317 * t424 + t319 * t407 + t321 * t408;
t194 = t369 * t424 + t370 * t407 + t371 * t408;
t231 = Icges(5,5) * t338 + Icges(5,6) * t337 + Icges(5,3) * t396;
t233 = Icges(5,4) * t338 + Icges(5,2) * t337 + Icges(5,6) * t396;
t235 = Icges(5,1) * t338 + Icges(5,4) * t337 + Icges(5,5) * t396;
t76 = t231 * t424 + t233 * t407 + t235 * t408 + t316 * t396 + t318 * t337 + t320 * t338;
t232 = Icges(5,5) * t340 + Icges(5,6) * t339 + Icges(5,3) * t398;
t234 = Icges(5,4) * t340 + Icges(5,2) * t339 + Icges(5,6) * t398;
t236 = Icges(5,1) * t340 + Icges(5,4) * t339 + Icges(5,5) * t398;
t77 = t232 * t424 + t234 * t407 + t236 * t408 + t317 * t396 + t319 * t337 + t321 * t338;
t17 = t113 * t456 + t172 * t396 + t173 * t398 + t194 * t420 + t424 * t76 + t426 * t77;
t3 = t129 * t396 + t130 * t398 + t156 * t420 + t33 * t424 + t34 * t426 + t456 * t49;
t546 = t17 / 0.2e1 + t13 / 0.2e1 + t3 / 0.2e1;
t111 = t270 * t426 + t271 * t513 + t272 * t405 - t334 * t359 + t335 * t360 + t358 * t398;
t168 = t305 * t426 + t307 * t513 + t309 * t405;
t169 = t306 * t426 + t308 * t513 + t310 * t405;
t191 = t358 * t426 + t359 * t513 + t360 * t405;
t74 = t218 * t426 + t220 * t513 + t222 * t405 + t305 * t398 - t307 * t334 + t309 * t335;
t75 = t219 * t426 + t221 * t513 + t223 * t405 + t306 * t398 - t308 * t334 + t310 * t335;
t14 = t111 * t456 + t168 * t396 + t169 * t398 + t191 * t420 + t424 * t74 + t426 * t75;
t114 = t275 * t426 + t276 * t409 + t277 * t410 + t339 * t370 + t340 * t371 + t369 * t398;
t174 = t316 * t426 + t318 * t409 + t320 * t410;
t175 = t317 * t426 + t319 * t409 + t321 * t410;
t195 = t369 * t426 + t370 * t409 + t371 * t410;
t78 = t231 * t426 + t233 * t409 + t235 * t410 + t316 * t398 + t318 * t339 + t320 * t340;
t79 = t232 * t426 + t234 * t409 + t236 * t410 + t317 * t398 + t319 * t339 + t321 * t340;
t18 = t114 * t456 + t174 * t396 + t175 * t398 + t195 * t420 + t424 * t78 + t426 * t79;
t4 = t131 * t396 + t132 * t398 + t157 * t420 + t35 * t424 + t36 * t426 + t456 * t50;
t545 = t18 / 0.2e1 + t14 / 0.2e1 + t4 / 0.2e1;
t15 = t110 * t478 + t458 * t72 + t459 * t73 + (t166 * t471 - t167 * t473 + t190 * t574) * t489;
t19 = t113 * t478 + t458 * t76 + t459 * t77 + (t172 * t471 - t173 * t473 + t194 * t574) * t489;
t5 = t33 * t458 + t34 * t459 + t478 * t49 + (t129 * t471 - t130 * t473 + t156 * t574) * t489;
t544 = t19 / 0.2e1 + t15 / 0.2e1 + t5 / 0.2e1;
t16 = t111 * t478 + t458 * t74 + t459 * t75 + (t168 * t471 - t169 * t473 + t191 * t574) * t489;
t20 = t114 * t478 + t458 * t78 + t459 * t79 + (t174 * t471 - t175 * t473 + t195 * t574) * t489;
t6 = t35 * t458 + t36 * t459 + t478 * t50 + (t131 * t471 - t132 * t473 + t157 * t574) * t489;
t543 = t20 / 0.2e1 + t16 / 0.2e1 + t6 / 0.2e1;
t117 = t270 * t456 + t271 * t512 + t272 * t418 + t358 * t420 - t359 * t376 + t360 * t377;
t176 = t305 * t456 + t307 * t512 + t309 * t418;
t177 = t306 * t456 + t308 * t512 + t310 * t418;
t215 = t358 * t456 + t359 * t512 + t360 * t418;
t84 = t218 * t456 + t220 * t512 + t222 * t418 + t305 * t420 - t307 * t376 + t309 * t377;
t85 = t219 * t456 + t221 * t512 + t223 * t418 + t306 * t420 - t308 * t376 + t310 * t377;
t21 = t117 * t456 + t176 * t396 + t177 * t398 + t215 * t420 + t424 * t84 + t426 * t85;
t121 = t275 * t456 + t276 * t428 + t277 * t429 + t369 * t420 + t370 * t385 + t371 * t386;
t179 = t316 * t456 + t318 * t428 + t320 * t429;
t180 = t317 * t456 + t319 * t428 + t321 * t429;
t230 = t369 * t456 + t370 * t428 + t371 * t429;
t87 = t231 * t456 + t233 * t428 + t235 * t429 + t316 * t420 + t318 * t385 + t320 * t386;
t88 = t232 * t456 + t234 * t428 + t236 * t429 + t317 * t420 + t319 * t385 + t321 * t386;
t23 = t121 * t456 + t179 * t396 + t180 * t398 + t230 * t420 + t424 * t87 + t426 * t88;
t8 = t139 * t396 + t140 * t398 + t163 * t420 + t37 * t424 + t38 * t426 + t456 * t67;
t542 = t23 / 0.2e1 + t21 / 0.2e1 + t8 / 0.2e1;
t22 = t117 * t478 + t458 * t84 + t459 * t85 + (t176 * t471 - t177 * t473 + t215 * t574) * t489;
t24 = t121 * t478 + t458 * t87 + t459 * t88 + (t179 * t471 - t180 * t473 + t230 * t574) * t489;
t9 = t37 * t458 + t38 * t459 + t478 * t67 + (t139 * t471 - t140 * t473 + t163 * t574) * t489;
t541 = t24 / 0.2e1 + t22 / 0.2e1 + t9 / 0.2e1;
t10 = t49 * t492 + (-t33 * t491 + t34 * t488) * t490;
t25 = t110 * t492 + (t488 * t73 - t491 * t72) * t490;
t27 = t113 * t492 + (t488 * t77 - t491 * t76) * t490;
t540 = t27 / 0.2e1 + t25 / 0.2e1 + t10 / 0.2e1;
t11 = t492 * t50 + (-t35 * t491 + t36 * t488) * t490;
t26 = t111 * t492 + (t488 * t75 - t491 * t74) * t490;
t28 = t114 * t492 + (t488 * t79 - t491 * t78) * t490;
t539 = t28 / 0.2e1 + t26 / 0.2e1 + t11 / 0.2e1;
t12 = t492 * t67 + (-t37 * t491 + t38 * t488) * t490;
t29 = t117 * t492 + (t488 * t85 - t491 * t84) * t490;
t30 = t121 * t492 + (t488 * t88 - t491 * t87) * t490;
t538 = t30 / 0.2e1 + t29 / 0.2e1 + t12 / 0.2e1;
t55 = t156 * t492 + (-t129 * t491 + t130 * t488) * t490;
t537 = t55 / 0.2e1 + ((-t166 - t172) * t491 + (t167 + t173) * t488) * t672 + (t194 + t190) * t648;
t56 = t157 * t492 + (-t131 * t491 + t132 * t488) * t490;
t536 = t56 / 0.2e1 + ((-t168 - t174) * t491 + (t169 + t175) * t488) * t672 + (t195 + t191) * t648;
t535 = (-t466 + t611) * t490;
t534 = (-t460 + t600) * t490;
t66 = t163 * t492 + (-t139 * t491 + t140 * t488) * t490;
t532 = t66 / 0.2e1 + ((-t176 - t179) * t491 + (t177 + t180) * t488) * t672 + (t230 + t215) * t648;
t527 = (-t466 + t582) * t490;
t526 = (-t460 + t577) * t490;
t524 = t478 * t214 + t299 * t548 + t357 * t637 + t578;
t210 = t213 * t634;
t211 = t214 * t629;
t522 = t210 + t211 + t552;
t521 = t298 * t634 + t299 * t629 + t551;
t517 = (-t466 + t559) * t490;
t516 = (-t460 + t555) * t490;
t68 = -t153 * t513 + t154 * t514 + t257 * t334 - t258 * t332;
t506 = m(5) * t237 + m(6) * t224 + m(7) * t624;
t330 = rSges(4,1) * t397 - rSges(4,2) * t396 + rSges(4,3) * t638;
t502 = m(4) * t330 + t506;
t477 = (rSges(3,1) * t500 - rSges(3,2) * t497) * t595;
t476 = (Icges(3,1) * t500 - t642) * t595;
t475 = (-Icges(3,2) * t497 + t641) * t595;
t474 = (Icges(3,5) * t500 - Icges(3,6) * t497) * t595;
t469 = t492 * rSges(3,3) + (rSges(3,1) * t497 + rSges(3,2) * t500) * t490;
t468 = Icges(3,5) * t492 + (Icges(3,1) * t497 + t641) * t490;
t467 = Icges(3,6) * t492 + (Icges(3,2) * t500 + t642) * t490;
t453 = -rSges(3,1) * t472 + rSges(3,2) * t473;
t452 = rSges(3,1) * t470 - rSges(3,2) * t471;
t451 = -Icges(3,1) * t472 + Icges(3,4) * t473;
t450 = Icges(3,1) * t470 - Icges(3,4) * t471;
t449 = -Icges(3,4) * t472 + Icges(3,2) * t473;
t448 = Icges(3,4) * t470 - Icges(3,2) * t471;
t447 = -Icges(3,5) * t472 + Icges(3,6) * t473;
t445 = -rSges(3,1) * t518 - rSges(3,2) * t519 + rSges(3,3) * t634;
t444 = t479 * rSges(3,1) + rSges(3,2) * t520 - rSges(3,3) * t629;
t441 = -Icges(3,1) * t518 - Icges(3,4) * t519 + Icges(3,5) * t634;
t440 = Icges(3,1) * t479 + Icges(3,4) * t520 - Icges(3,5) * t629;
t439 = -Icges(3,4) * t518 - Icges(3,2) * t519 + Icges(3,6) * t634;
t438 = Icges(3,4) * t479 + Icges(3,2) * t520 - Icges(3,6) * t629;
t414 = Icges(4,1) * t457 - Icges(4,4) * t456 + Icges(4,5) * t478;
t413 = Icges(4,4) * t457 - Icges(4,2) * t456 + Icges(4,6) * t478;
t412 = Icges(4,5) * t457 - Icges(4,6) * t456 + Icges(4,3) * t478;
t383 = Icges(4,1) * t421 - Icges(4,4) * t420 + Icges(4,5) * t548;
t382 = Icges(4,4) * t421 - Icges(4,2) * t420 + Icges(4,6) * t548;
t381 = Icges(4,5) * t421 - Icges(4,6) * t420 + Icges(4,3) * t548;
t373 = rSges(4,1) * t427 - rSges(4,2) * t426 + rSges(4,3) * t459;
t372 = rSges(4,1) * t425 - rSges(4,2) * t424 + rSges(4,3) * t458;
t368 = Icges(4,1) * t427 - Icges(4,4) * t426 + Icges(4,5) * t459;
t367 = Icges(4,1) * t425 - Icges(4,4) * t424 + Icges(4,5) * t458;
t366 = Icges(4,4) * t427 - Icges(4,2) * t426 + Icges(4,6) * t459;
t365 = Icges(4,4) * t425 - Icges(4,2) * t424 + Icges(4,6) * t458;
t364 = Icges(4,5) * t427 - Icges(4,6) * t426 + Icges(4,3) * t459;
t363 = Icges(4,5) * t425 - Icges(4,6) * t424 + Icges(4,3) * t458;
t348 = t424 * t357;
t329 = Icges(4,1) * t399 - Icges(4,4) * t398 - Icges(4,5) * t637;
t328 = Icges(4,1) * t397 - Icges(4,4) * t396 + Icges(4,5) * t638;
t327 = Icges(4,4) * t399 - Icges(4,2) * t398 - Icges(4,6) * t637;
t326 = Icges(4,4) * t397 - Icges(4,2) * t396 + Icges(4,6) * t638;
t325 = Icges(4,5) * t399 - Icges(4,6) * t398 - Icges(4,3) * t637;
t324 = Icges(4,5) * t397 - Icges(4,6) * t396 + Icges(4,3) * t638;
t293 = t373 * t478 - t415 * t459;
t292 = -t372 * t478 + t415 * t458;
t285 = t456 * t299;
t282 = t426 * t298;
t281 = (-t372 - t435) * t492 + t491 * t561;
t280 = t373 * t492 + t488 * t561 + t430;
t274 = t412 * t478 - t413 * t456 + t414 * t457;
t268 = t372 * t459 - t373 * t458;
t263 = t412 * t459 - t413 * t426 + t414 * t427;
t262 = t412 * t458 - t413 * t424 + t414 * t425;
t259 = (t372 * t488 + t373 * t491) * t490 + t598;
t250 = (-t330 - t442) * t492 + t491 * t562;
t249 = t331 * t492 + t488 * t562 + t437;
t248 = t323 * t456 - t374 * t426;
t247 = -t322 * t456 + t374 * t424;
t240 = t364 * t478 - t366 * t456 + t368 * t457;
t239 = t363 * t478 - t365 * t456 + t367 * t457;
t229 = t364 * t459 - t366 * t426 + t368 * t427;
t228 = t363 * t459 - t365 * t426 + t367 * t427;
t227 = t364 * t458 - t366 * t424 + t368 * t425;
t226 = t363 * t458 - t365 * t424 + t367 * t425;
t217 = (t330 * t488 + t331 * t491) * t490 + t596;
t216 = t322 * t426 - t323 * t424;
t202 = t323 * t478 + t459 * t600 + t387;
t201 = t374 * t458 + t478 * t604 + t406;
t197 = (-t435 + t604) * t492 + t491 * t534;
t196 = t323 * t492 + t488 * t534 + t599;
t193 = t331 * t478 - t384 * t459 + (t373 * t574 + t415 * t473) * t489;
t192 = -t330 * t478 + t384 * t458 + (-t372 * t574 + t415 * t471) * t489;
t185 = -t258 * t512 + t304 * t513;
t184 = t257 * t512 - t304 * t514;
t183 = t322 * t459 + t375 + (-t323 - t395) * t458;
t182 = (t322 * t488 + t323 * t491) * t490 + t551;
t181 = t381 * t478 - t382 * t456 + t383 * t457 + t412 * t548 - t413 * t420 + t414 * t421;
t178 = t330 * t459 - t331 * t458 + (-t372 * t473 - t373 * t471) * t489;
t171 = t312 * t456 + t426 * t601 + t285;
t170 = t361 * t424 + t456 * t608 + t348;
t165 = t381 * t459 - t382 * t426 + t383 * t427 - t398 * t413 + t399 * t414 - t412 * t637;
t164 = t381 * t458 - t382 * t424 + t383 * t425 - t396 * t413 + t397 * t414 + t412 * t638;
t162 = -t257 * t513 + t258 * t514;
t161 = (-t435 + t580) * t492 + t491 * t526;
t160 = t312 * t492 + t488 * t526 + t581;
t159 = t312 * t478 + t459 * t577 + t609;
t158 = t361 * t458 + t478 * t580 + t602;
t155 = t311 * t426 + t424 * t607 + t282;
t146 = (-t442 + t617) * t492 + t491 * t535;
t145 = t238 * t492 + t488 * t535 + t603;
t144 = t325 * t478 - t327 * t456 + t329 * t457 + t364 * t548 - t366 * t420 + t368 * t421;
t143 = t324 * t478 - t326 * t456 + t328 * t457 + t363 * t548 - t365 * t420 + t367 * t421;
t142 = (t311 * t488 + t312 * t491) * t490 + t521;
t141 = t311 * t459 + (-t395 + t607) * t458 + t610;
t138 = t325 * t459 - t327 * t426 + t329 * t427 - t364 * t637 - t366 * t398 + t368 * t399;
t137 = t324 * t459 - t326 * t426 + t328 * t427 - t363 * t637 - t365 * t398 + t367 * t399;
t136 = t325 * t458 - t327 * t424 + t329 * t425 + t364 * t638 - t366 * t396 + t368 * t397;
t135 = t324 * t458 - t326 * t424 + t328 * t425 + t363 * t638 - t365 * t396 + t367 * t397;
t134 = t238 * t456 - t278 * t426 + t323 * t420 - t374 * t398;
t133 = -t237 * t456 + t278 * t424 - t322 * t420 + t374 * t396;
t128 = t426 * t579 + t456 * t615 + t285;
t127 = t424 * t606 + t456 * t584 + t348;
t126 = (t237 * t488 + t238 * t491) * t490 + t552;
t125 = (-t435 + t556) * t492 + t491 * t516;
t124 = t488 * t516 + t492 * t615 + t581;
t123 = t459 * t555 + t478 * t615 + t609;
t122 = t458 * t606 + t478 * t556 + t602;
t120 = t237 * t426 - t238 * t424 + t322 * t398 - t323 * t396;
t119 = t238 * t478 + (t323 * t574 + t374 * t473) * t489 + t611 * t459 + t578;
t118 = t278 * t458 + t362 + t617 * t478 + (-t471 * t600 + t574 * t604) * t489;
t116 = t424 * t583 + t426 * t616 + t282;
t115 = (t488 * t616 + t491 * t615) * t490 + t521;
t112 = t616 * t459 + (-t395 + t583) * t458 + t610;
t109 = (-t442 + t585) * t492 + t491 * t527;
t108 = t225 * t492 + t488 * t527 + t586;
t106 = t179 * t458 + t180 * t459 + t230 * t478;
t104 = t179 * t424 + t180 * t426 + t230 * t456;
t103 = t237 * t459 + (-t238 - t347) * t458 + (-t323 * t471 + t473 * t604) * t489 + t605;
t102 = t176 * t458 + t177 * t459 + t215 * t478;
t101 = t176 * t424 + t177 * t426 + t215 * t456;
t98 = t174 * t458 + t175 * t459 + t195 * t478;
t97 = t172 * t458 + t173 * t459 + t194 * t478;
t94 = t174 * t424 + t175 * t426 + t195 * t456;
t93 = t172 * t424 + t173 * t426 + t194 * t456;
t92 = t168 * t458 + t169 * t459 + t191 * t478;
t91 = t166 * t458 + t167 * t459 + t190 * t478;
t90 = t168 * t424 + t169 * t426 + t191 * t456;
t89 = t166 * t424 + t167 * t426 + t190 * t456;
t86 = (t224 * t488 + t225 * t491) * t490 + t522;
t83 = t225 * t456 + t312 * t420 + t398 * t601 + t426 * t612 + t620;
t82 = t273 * t424 + t361 * t396 + t420 * t608 + t456 * t619 + t614;
t81 = -t154 * t512 + t189 * t513 + t258 * t376 - t304 * t334;
t80 = t153 * t512 - t189 * t514 - t257 * t376 + t304 * t332;
t71 = t225 * t478 + (t312 * t574 + t361 * t473) * t489 + t582 * t459 + t524;
t70 = t273 * t458 + t585 * t478 + (-t471 * t577 + t574 * t580) * t489 + t613;
t69 = t181 * t492 + (-t143 * t491 + t144 * t488) * t490;
t65 = t224 * t426 + t311 * t398 + t396 * t607 + t424 * t618 + t621;
t64 = t165 * t492 + (-t137 * t491 + t138 * t488) * t490;
t63 = t164 * t492 + (-t135 * t491 + t136 * t488) * t490;
t62 = t139 * t458 + t140 * t459 + t163 * t478;
t61 = t139 * t424 + t140 * t426 + t163 * t456;
t60 = (-t442 + t560) * t492 + t491 * t517;
t59 = t488 * t517 + t492 * t623 + t586;
t57 = t224 * t459 + (-t347 + t618) * t458 + (-t312 * t471 + t473 * t580) * t489 + t557;
t54 = t131 * t458 + t132 * t459 + t157 * t478;
t53 = t129 * t458 + t130 * t459 + t156 * t478;
t52 = t131 * t424 + t132 * t426 + t157 * t456;
t51 = t129 * t424 + t130 * t426 + t156 * t456;
t48 = -t131 * t514 - t132 * t513 - t157 * t512;
t46 = t143 * t458 + t144 * t459 + t181 * t478 + (t239 * t471 - t240 * t473 + t274 * t574) * t489;
t45 = (t488 * t624 + t491 * t623) * t490 + t522;
t44 = t137 * t458 + t138 * t459 + t165 * t478 + (t228 * t471 - t229 * t473 + t263 * t574) * t489;
t43 = t135 * t458 + t136 * t459 + t164 * t478 + (t226 * t471 - t227 * t473 + t262 * t574) * t489;
t42 = t398 * t579 + t420 * t615 + t426 * t587 + t456 * t623 + t620;
t41 = t396 * t606 + t420 * t584 + t424 * t622 + t456 * t589 + t614;
t40 = t623 * t478 + (t473 * t606 + t574 * t615) * t489 + t559 * t459 + t524;
t39 = t622 * t458 + t560 * t478 + (-t471 * t555 + t556 * t574) * t489 + t613;
t32 = t396 * t583 + t398 * t616 + t424 * t588 + t426 * t624 + t621;
t31 = t624 * t459 + (-t347 + t588) * t458 + (-t471 * t615 + t473 * t556) * t489 + t557;
t95 = [0; m(4) * t597 / 0.2e1 + t553 * t666 + (m(3) * t453 - t671) * t629 + (m(3) * t452 + t502) * t634 + t590 * (0.2e1 * t210 + 0.2e1 * t211 + t553); t11 * t634 - t10 * t629 + t28 * t634 - t27 * t629 + t26 * t634 - t25 * t629 + t64 * t634 - t63 * t629 + ((t473 * t439 - t472 * t441 + t447 * t634 - t449 * t519 - t451 * t518) * t634 - (t473 * t438 - t472 * t440 + t446 * t634 - t448 * t519 - t450 * t518) * t629 + (t473 * t467 - t472 * t468 + t474 * t634 - t475 * t519 - t476 * t518) * t492) * t634 - ((-t471 * t439 + t470 * t441 - t447 * t629 + t449 * t520 + t479 * t451) * t634 - (-t471 * t438 + t470 * t440 + t448 * t520 + t479 * t450 - t490 * t628) * t629 + (-t471 * t467 + t470 * t468 - t474 * t629 + t475 * t520 + t479 * t476) * t492) * t629 + (t115 * t45 + t124 * t59 + t125 * t60) * t564 + t492 * t12 + (t108 * t160 + t109 * t161 + t142 * t86) * t566 + t492 * t29 + t492 * t30 + (t126 * t182 + t145 * t196 + t146 * t197) * t568 + (t217 * t259 + t249 * t280 + t250 * t281) * t570 + t492 * t69 + t492 * (t492 ^ 2 * t474 + (((t449 * t500 + t451 * t497) * t488 - (t448 * t500 + t450 * t497) * t491 + ((-t439 * t497 + t441 * t500) * t488 - (-t438 * t497 + t440 * t500) * t491) * qJD(2)) * t490 + (-t628 + t488 * t447 + t475 * t500 + t476 * t497 + (-t467 * t497 + t468 * t500) * qJD(2)) * t492) * t490) + 0.2e1 * m(3) * ((-t444 * t492 - t469 * t629) * (-t452 * t492 - t477 * t629) + (t445 * t492 - t469 * t634) * (t453 * t492 - t477 * t634) + (t444 * t488 + t445 * t491) * t490 ^ 2 * (t452 * t488 + t453 * t491)); t554 * t666 + t502 * t459 + ((-m(4) * t372 - t668) * t473 - (m(4) * t373 - t669) * t471) * t489 + t671 * t458 + t590 * (t214 * t667 + t298 * t591 + 0.2e1 * t206 - 0.2e1 * t289 + t554); m(4) * (t178 * t259 + t192 * t281 + t193 * t280 + t217 * t268 + t249 * t293 + t250 * t292) + (t103 * t182 + t118 * t197 + t119 * t196 + t126 * t183 + t145 * t202 + t146 * t201) * m(5) + (t108 * t159 + t109 * t158 + t141 * t86 + t142 * t57 + t160 * t71 + t161 * t70) * m(6) + (t112 * t45 + t115 * t31 + t122 * t60 + t123 * t59 + t124 * t40 + t125 * t39) * m(7) + (t46 / 0.2e1 + t541) * t492 + (t69 / 0.2e1 + t538) * t478 + (t64 / 0.2e1 + t539) * t459 + (t63 / 0.2e1 + t540) * t458 + ((-t43 / 0.2e1 - t544) * t491 + (t44 / 0.2e1 + t543) * t488) * t490 + ((t263 * t649 + (t228 * t650 + t229 * t653) * t490 - t536) * t473 - (t262 * t649 + (t226 * t650 + t227 * t653) * t490 - t537) * t471 + (t274 * t648 + (t239 * t651 + t240 * t652) * t490 + t532) * t574) * t489; (t178 * t268 + t192 * t292 + t193 * t293) * t570 + (t103 * t183 + t118 * t201 + t119 * t202) * t568 + (t141 * t57 + t158 * t70 + t159 * t71) * t566 + (t112 * t31 + t122 * t39 + t123 * t40) * t564 + (t24 + t22 + t9 + t46) * t478 + (t20 + t16 + t6 + t44) * t459 + (t19 + t15 + t5 + t43) * t458 + ((-t228 * t458 - t229 * t459 - t263 * t478 - t54 - t92 - t98) * t473 - (-t226 * t458 - t227 * t459 - t262 * t478 - t53 - t91 - t97) * t471 + (t239 * t458 + t240 * t459 + t274 * t478 + t102 + t106 + t62) * t574) * t489; t506 * t426 + t505 * t424 + t668 * t398 + t669 * t396 + (-t424 * t214 - t396 * t299 + t621) * t674; (t120 * t182 + t126 * t216 + t133 * t197 + t134 * t196 + t145 * t248 + t146 * t247) * m(5) + (t108 * t171 + t109 * t170 + t142 * t65 + t155 * t86 + t160 * t83 + t161 * t82) * m(6) + (t115 * t32 + t116 * t45 + t124 * t42 + t125 * t41 + t127 * t60 + t128 * t59) * m(7) + t542 * t492 + t538 * t456 + t539 * t426 + t540 * t424 + t532 * t420 + t536 * t398 + t537 * t396 + (t488 * t545 - t491 * t546) * t490; (t103 * t216 + t118 * t247 + t119 * t248 + t120 * t183 + t133 * t201 + t134 * t202) * m(5) + (t141 * t65 + t155 * t57 + t158 * t82 + t159 * t83 + t170 * t70 + t171 * t71) * m(6) + (t112 * t32 + t116 * t31 + t122 * t41 + t123 * t42 + t127 * t39 + t128 * t40) * m(7) + t542 * t478 + t545 * t459 + t546 * t458 + t541 * t456 + t543 * t426 + t544 * t424 + (t106 / 0.2e1 + t102 / 0.2e1 + t62 / 0.2e1) * t420 + (t98 / 0.2e1 + t92 / 0.2e1 + t54 / 0.2e1) * t398 + (t97 / 0.2e1 + t91 / 0.2e1 + t53 / 0.2e1) * t396 + ((-t90 / 0.2e1 - t94 / 0.2e1 - t52 / 0.2e1) * t473 - (-t89 / 0.2e1 - t93 / 0.2e1 - t51 / 0.2e1) * t471 + (t61 / 0.2e1 + t101 / 0.2e1 + t104 / 0.2e1) * t574) * t489; (t116 * t32 + t127 * t41 + t128 * t42) * t564 + (t155 * t65 + t170 * t82 + t171 * t83) * t566 + (t120 * t216 + t133 * t247 + t134 * t248) * t568 + (t8 + t23 + t21) * t456 + (t4 + t18 + t14) * t426 + (t3 + t17 + t13) * t424 + (t61 + t104 + t101) * t420 + (t52 + t94 + t90) * t398 + (t51 + t93 + t89) * t396; t420 * t674; (m(6) * t86 + m(7) * t45) * t456 + (m(6) * t109 + m(7) * t60) * t426 + (m(6) * t108 + m(7) * t59) * t424 + (m(6) * t142 + m(7) * t115) * t420 + (m(6) * t161 + m(7) * t125) * t398 + (m(6) * t160 + m(7) * t124) * t396; (m(6) * t57 + m(7) * t31) * t456 + (m(6) * t70 + m(7) * t39) * t426 + (m(6) * t71 + m(7) * t40) * t424 + (m(6) * t141 + m(7) * t112) * t420 + (m(6) * t158 + m(7) * t122) * t398 + (m(6) * t159 + m(7) * t123) * t396; (m(6) * t65 + m(7) * t32) * t456 + (m(6) * t82 + m(7) * t41) * t426 + (m(6) * t83 + m(7) * t42) * t424 + (m(6) * t155 + m(7) * t116) * t420 + (m(6) * t170 + m(7) * t127) * t398 + (m(6) * t171 + m(7) * t128) * t396; 0.4e1 * t590 * (t396 * t424 + t398 * t426 + t420 * t456); t68 * m(7); (t115 * t68 + t124 * t81 + t125 * t80 + t162 * t45 + t184 * t60 + t185 * t59) * m(7) + t55 * t659 + t10 * t656 + t56 * t658 + t11 * t655 + t66 * t657 + t12 * t654 + t7 * t648 + (t1 * t651 + t2 * t652) * t490; t459 * t662 + t478 * t661 + t62 * t657 + t9 * t654 + t53 * t659 + t5 * t656 + t458 * t663 + t54 * t658 + t6 * t655 + (t112 * t68 + t122 * t80 + t123 * t81 + t162 * t31 + t184 * t39 + t185 * t40) * m(7) + (-t473 * t48 / 0.2e1 + t574 * t660 + t471 * t673) * t489; (t116 * t68 + t127 * t80 + t128 * t81 + t162 * t32 + t184 * t41 + t185 * t42) * m(7) + t420 * t660 + t456 * t661 + t61 * t657 + t8 * t654 + t51 * t659 + t3 * t656 + t52 * t658 + t4 * t655 + t398 * t48 / 0.2e1 + t426 * t662 + t396 * t673 + t424 * t663; (t162 * t420 + t184 * t398 + t185 * t396 + t424 * t81 + t426 * t80 + t456 * t68) * m(7); t332 * t47 - t514 * t1 + t376 * t58 - t512 * t7 + t334 * t48 - t513 * t2 + (t162 * t68 + t184 * t80 + t185 * t81) * t564;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t95(1) t95(2) t95(4) t95(7) t95(11) t95(16); t95(2) t95(3) t95(5) t95(8) t95(12) t95(17); t95(4) t95(5) t95(6) t95(9) t95(13) t95(18); t95(7) t95(8) t95(9) t95(10) t95(14) t95(19); t95(11) t95(12) t95(13) t95(14) t95(15) t95(20); t95(16) t95(17) t95(18) t95(19) t95(20) t95(21);];
Mq  = res;
