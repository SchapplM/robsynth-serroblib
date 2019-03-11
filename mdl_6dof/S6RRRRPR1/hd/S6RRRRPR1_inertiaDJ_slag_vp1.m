% Calculate time derivative of joint inertia matrix for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:52:57
% EndTime: 2019-03-09 21:53:43
% DurationCPUTime: 26.60s
% Computational Cost: add. (61923->1046), mult. (48097->1402), div. (0->0), fcn. (43323->12), ass. (0->556)
t428 = sin(qJ(1));
t431 = cos(qJ(1));
t425 = qJ(2) + qJ(3);
t414 = qJ(4) + t425;
t400 = pkin(11) + t414;
t395 = sin(t400);
t396 = cos(t400);
t508 = Icges(6,5) * t396 - Icges(6,6) * t395;
t274 = Icges(6,3) * t428 + t431 * t508;
t401 = sin(t414);
t402 = cos(t414);
t509 = Icges(5,5) * t402 - Icges(5,6) * t401;
t303 = Icges(5,3) * t428 + t431 * t509;
t741 = t274 + t303;
t349 = Icges(6,5) * t395 + Icges(6,6) * t396;
t359 = Icges(5,5) * t401 + Icges(5,6) * t402;
t735 = t349 + t359;
t669 = Icges(6,4) * t395;
t350 = Icges(6,2) * t396 + t669;
t668 = Icges(6,4) * t396;
t351 = Icges(6,1) * t395 + t668;
t671 = Icges(5,4) * t401;
t360 = Icges(5,2) * t402 + t671;
t670 = Icges(5,4) * t402;
t361 = Icges(5,1) * t401 + t670;
t740 = t350 * t395 - t351 * t396 + t360 * t401 - t361 * t402;
t273 = -Icges(6,3) * t431 + t428 * t508;
t302 = -Icges(5,3) * t431 + t428 * t509;
t513 = -Icges(6,2) * t395 + t668;
t275 = -Icges(6,6) * t431 + t428 * t513;
t519 = Icges(6,1) * t396 - t669;
t277 = -Icges(6,5) * t431 + t428 * t519;
t501 = t275 * t395 - t277 * t396;
t720 = t431 * t501;
t514 = -Icges(5,2) * t401 + t670;
t304 = -Icges(5,6) * t431 + t428 * t514;
t520 = Icges(5,1) * t402 - t671;
t306 = -Icges(5,5) * t431 + t428 * t520;
t499 = t304 * t401 - t306 * t402;
t721 = t431 * t499;
t739 = t720 + t721 + (-t273 - t302) * t428;
t416 = t428 * rSges(5,3);
t685 = rSges(5,1) * t402;
t532 = -rSges(5,2) * t401 + t685;
t313 = t431 * t532 + t416;
t417 = t428 * rSges(4,3);
t411 = sin(t425);
t412 = cos(t425);
t686 = rSges(4,1) * t412;
t533 = -rSges(4,2) * t411 + t686;
t329 = t431 * t533 + t417;
t305 = Icges(5,6) * t428 + t431 * t514;
t307 = Icges(5,5) * t428 + t431 * t520;
t498 = t305 * t401 - t307 * t402;
t276 = Icges(6,6) * t428 + t431 * t513;
t278 = Icges(6,5) * t428 + t431 * t519;
t500 = t276 * t395 - t278 * t396;
t736 = (-t498 - t500) * t431 + t741 * t428;
t421 = qJD(2) + qJD(3);
t408 = qJD(4) + t421;
t280 = t513 * t408;
t281 = t519 * t408;
t309 = t514 * t408;
t310 = t520 * t408;
t639 = t361 * t408;
t640 = t360 * t408;
t641 = t351 * t408;
t642 = t350 * t408;
t734 = (t310 - t640) * t402 + (-t309 - t639) * t401 + (t281 - t642) * t396 + (-t280 - t641) * t395 + t735 * qJD(1);
t432 = -pkin(8) - pkin(7);
t427 = sin(qJ(2));
t591 = qJD(2) * t427;
t584 = pkin(2) * t591;
t732 = qJD(1) * t432 + t584;
t429 = cos(qJ(6));
t426 = sin(qJ(6));
t666 = Icges(7,4) * t429;
t512 = -Icges(7,2) * t426 + t666;
t633 = t396 * t408;
t667 = Icges(7,4) * t426;
t161 = t512 * t633 + (Icges(7,6) * t408 + (-Icges(7,2) * t429 - t667) * qJD(6)) * t395;
t263 = -Icges(7,6) * t396 + t395 * t512;
t518 = Icges(7,1) * t429 - t667;
t264 = -Icges(7,5) * t396 + t395 * t518;
t731 = -t161 * t426 + (-t263 * t429 - t264 * t426) * qJD(6);
t730 = (t508 + t509) * t408 + t740 * qJD(1);
t355 = rSges(6,1) * t395 + rSges(6,2) * t396;
t478 = t355 * t408;
t679 = rSges(5,2) * t402;
t363 = rSges(5,1) * t401 + t679;
t479 = t363 * t408;
t430 = cos(qJ(2));
t387 = rSges(3,1) * t427 + rSges(3,2) * t430;
t476 = qJD(2) * t387;
t729 = t428 * t476;
t674 = Icges(3,4) * t430;
t517 = -Icges(3,2) * t427 + t674;
t345 = Icges(3,6) * t428 + t431 * t517;
t675 = Icges(3,4) * t427;
t523 = Icges(3,1) * t430 - t675;
t347 = Icges(3,5) * t428 + t431 * t523;
t494 = t345 * t427 - t347 * t430;
t728 = t428 * t494;
t672 = Icges(4,4) * t412;
t515 = -Icges(4,2) * t411 + t672;
t324 = Icges(4,6) * t428 + t431 * t515;
t673 = Icges(4,4) * t411;
t521 = Icges(4,1) * t412 - t673;
t326 = Icges(4,5) * t428 + t431 * t521;
t496 = t324 * t411 - t326 * t412;
t727 = t428 * t496;
t726 = t428 * t498;
t725 = t428 * t500;
t403 = t430 * pkin(2) + pkin(1);
t689 = pkin(1) - t403;
t724 = t428 * t689;
t344 = -Icges(3,6) * t431 + t428 * t517;
t346 = -Icges(3,5) * t431 + t428 * t523;
t495 = t344 * t427 - t346 * t430;
t723 = t431 * t495;
t323 = -Icges(4,6) * t431 + t428 * t515;
t325 = -Icges(4,5) * t431 + t428 * t521;
t497 = t323 * t411 - t325 * t412;
t722 = t431 * t497;
t507 = Icges(7,5) * t429 - Icges(7,6) * t426;
t160 = t507 * t633 + (Icges(7,3) * t408 + (-Icges(7,5) * t426 - Icges(7,6) * t429) * qJD(6)) * t395;
t652 = t263 * t426;
t719 = -t408 * t652 - t160;
t528 = rSges(7,1) * t429 - rSges(7,2) * t426;
t167 = t528 * t633 + (rSges(7,3) * t408 + (-rSges(7,1) * t426 - rSges(7,2) * t429) * qJD(6)) * t395;
t691 = pkin(5) * t396;
t535 = pkin(10) * t395 + t691;
t718 = -t535 * t408 - t167;
t618 = t429 * t431;
t621 = t426 * t428;
t338 = -t396 * t621 - t618;
t619 = t428 * t429;
t620 = t426 * t431;
t339 = t396 * t619 - t620;
t529 = -rSges(7,1) * t339 - rSges(7,2) * t338;
t635 = t395 * t428;
t221 = rSges(7,3) * t635 - t529;
t340 = -t396 * t620 + t619;
t341 = t396 * t618 + t621;
t634 = t395 * t431;
t222 = t341 * rSges(7,1) + t340 * rSges(7,2) + rSges(7,3) * t634;
t717 = -t428 * t221 - t431 * t222;
t632 = t396 * t431;
t319 = pkin(5) * t632 + pkin(10) * t634;
t716 = -t222 - t319;
t265 = -rSges(7,3) * t396 + t395 * t528;
t692 = pkin(5) * t395;
t356 = -pkin(10) * t396 + t692;
t715 = -t265 - t356;
t714 = qJD(1) * t273;
t713 = qJD(1) * t302;
t510 = Icges(4,5) * t412 - Icges(4,6) * t411;
t321 = -Icges(4,3) * t431 + t428 * t510;
t712 = qJD(1) * t321;
t511 = Icges(3,5) * t430 - Icges(3,6) * t427;
t342 = -Icges(3,3) * t431 + t428 * t511;
t540 = qJD(1) * t396 - qJD(6);
t627 = t408 * t431;
t580 = t395 * t627;
t711 = t428 * t540 + t580;
t366 = Icges(4,2) * t412 + t673;
t367 = Icges(4,1) * t411 + t672;
t491 = t366 * t411 - t367 * t412;
t707 = qJD(1) * t491 + t510 * t421;
t422 = -pkin(9) + t432;
t600 = t422 - t432;
t374 = pkin(3) * t412 + t403;
t607 = t374 - t403;
t258 = t428 * t607 + t431 * t600;
t413 = -qJ(5) + t422;
t601 = t413 - t422;
t348 = pkin(4) * t402 + t374;
t609 = t348 - t374;
t226 = t428 * t609 + t431 * t601;
t706 = 2 * m(3);
t705 = 2 * m(4);
t704 = 2 * m(5);
t703 = 2 * m(6);
t702 = 2 * m(7);
t423 = t428 ^ 2;
t424 = t431 ^ 2;
t701 = t428 / 0.2e1;
t700 = -t431 / 0.2e1;
t699 = -rSges(7,3) - pkin(10);
t698 = m(3) * t387;
t681 = rSges(4,2) * t412;
t368 = rSges(4,1) * t411 + t681;
t697 = m(4) * t368;
t696 = m(5) * t363;
t695 = pkin(2) * t427;
t694 = pkin(3) * t411;
t693 = pkin(4) * t401;
t690 = t428 * pkin(7);
t420 = t431 * pkin(7);
t688 = -pkin(7) - t432;
t687 = rSges(3,1) * t430;
t684 = rSges(6,1) * t396;
t683 = rSges(3,2) * t427;
t678 = rSges(3,3) * t431;
t541 = -qJD(6) * t396 + qJD(1);
t488 = t541 * t429;
t629 = t408 * t428;
t581 = t395 * t629;
t192 = t428 * t488 + (-t431 * t540 + t581) * t426;
t489 = t541 * t426;
t628 = t408 * t429;
t193 = t540 * t618 + (-t395 * t628 + t489) * t428;
t593 = qJD(1) * t431;
t457 = t395 * t593 + t396 * t629;
t101 = Icges(7,5) * t193 + Icges(7,6) * t192 + Icges(7,3) * t457;
t103 = Icges(7,4) * t193 + Icges(7,2) * t192 + Icges(7,6) * t457;
t105 = Icges(7,1) * t193 + Icges(7,4) * t192 + Icges(7,5) * t457;
t213 = Icges(7,5) * t339 + Icges(7,6) * t338 + Icges(7,3) * t635;
t215 = Icges(7,4) * t339 + Icges(7,2) * t338 + Icges(7,6) * t635;
t217 = Icges(7,1) * t339 + Icges(7,4) * t338 + Icges(7,5) * t635;
t506 = -t215 * t426 + t217 * t429;
t27 = (t408 * t506 - t101) * t396 + (-t103 * t426 + t105 * t429 + t213 * t408 + (-t215 * t429 - t217 * t426) * qJD(6)) * t395;
t677 = t27 * t431;
t190 = t426 * t711 + t431 * t488;
t191 = -t429 * t711 + t431 * t489;
t594 = qJD(1) * t428;
t565 = t395 * t594;
t579 = t396 * t627;
t456 = -t565 + t579;
t100 = Icges(7,5) * t191 + Icges(7,6) * t190 + Icges(7,3) * t456;
t102 = Icges(7,4) * t191 + Icges(7,2) * t190 + Icges(7,6) * t456;
t104 = Icges(7,1) * t191 + Icges(7,4) * t190 + Icges(7,5) * t456;
t214 = Icges(7,5) * t341 + Icges(7,6) * t340 + Icges(7,3) * t634;
t216 = Icges(7,4) * t341 + Icges(7,2) * t340 + Icges(7,6) * t634;
t218 = Icges(7,1) * t341 + Icges(7,4) * t340 + Icges(7,5) * t634;
t505 = -t216 * t426 + t218 * t429;
t28 = (t408 * t505 - t100) * t396 + (-t102 * t426 + t104 * t429 + t214 * t408 + (-t216 * t429 - t218 * t426) * qJD(6)) * t395;
t676 = t28 * t428;
t418 = t428 * rSges(3,3);
t415 = t428 * rSges(6,3);
t311 = t532 * t408;
t649 = t311 * t428;
t648 = t311 * t431;
t337 = t533 * t421;
t647 = t337 * t428;
t646 = t344 * t430;
t645 = t345 * t430;
t644 = t346 * t427;
t643 = t347 * t427;
t638 = t366 * t421;
t637 = t367 * t421;
t636 = t395 * t408;
t631 = t401 * t408;
t630 = t402 * t408;
t626 = t411 * t421;
t625 = t412 * t421;
t624 = t413 * t431;
t623 = t421 * t428;
t622 = t421 * t431;
t330 = t431 * t348;
t364 = t431 * t374;
t227 = -t428 * t601 + t330 - t364;
t617 = t428 * t226 + t431 * t227;
t287 = rSges(6,1) * t632 - rSges(6,2) * t634 + t415;
t616 = -t227 - t287;
t389 = t431 * t403;
t259 = -t428 * t600 + t364 - t389;
t615 = t428 * t258 + t431 * t259;
t614 = -t259 - t313;
t362 = -pkin(3) * t626 - t584;
t292 = -pkin(4) * t631 + t362;
t613 = qJD(5) * t428 + t431 * t292;
t312 = -rSges(5,3) * t431 + t428 * t532;
t211 = t428 * t312 + t431 * t313;
t316 = t431 * t432 + t420 - t724;
t317 = -t431 * pkin(1) + t428 * t688 + t389;
t612 = t428 * t316 + t431 * t317;
t328 = -t431 * rSges(4,3) + t428 * t533;
t225 = t428 * t328 + t431 * t329;
t564 = t401 * t594;
t380 = pkin(4) * t564;
t611 = t355 * t594 + t380;
t563 = t411 * t594;
t382 = pkin(3) * t563;
t610 = t363 * t594 + t382;
t608 = rSges(6,2) * t565 + rSges(6,3) * t593;
t606 = rSges(5,2) * t564 + rSges(5,3) * t593;
t605 = rSges(4,2) * t563 + rSges(4,3) * t593;
t604 = qJD(5) * t431 + t413 * t594;
t603 = t732 * t428;
t602 = t431 * t687 + t418;
t599 = t423 + t424;
t598 = qJD(1) * t274;
t597 = qJD(1) * t303;
t322 = Icges(4,3) * t428 + t431 * t510;
t596 = qJD(1) * t322;
t343 = Icges(3,3) * t428 + t431 * t511;
t595 = qJD(1) * t343;
t590 = qJD(2) * t430;
t122 = -t273 * t431 - t428 * t501;
t123 = -t274 * t431 - t725;
t468 = t408 * t349;
t178 = -t431 * t468 - t714;
t179 = -t428 * t468 + t598;
t180 = -qJD(1) * t275 - t431 * t642;
t182 = -qJD(1) * t277 - t431 * t641;
t183 = qJD(1) * t278 - t428 * t641;
t552 = t275 * t408 - t183;
t181 = qJD(1) * t276 - t428 * t642;
t554 = t277 * t408 + t181;
t14 = (t431 * t179 + (t123 + t720) * qJD(1)) * t431 + (t122 * qJD(1) + (-t180 * t395 + t182 * t396 - t276 * t633 - t278 * t636 + t598) * t428 + (-t178 + t552 * t396 + t554 * t395 + (-t273 - t500) * qJD(1)) * t431) * t428;
t135 = -t302 * t431 - t428 * t499;
t136 = -t303 * t431 - t726;
t469 = t408 * t359;
t198 = -t431 * t469 - t713;
t199 = -t428 * t469 + t597;
t200 = -qJD(1) * t304 - t431 * t640;
t202 = -qJD(1) * t306 - t431 * t639;
t203 = qJD(1) * t307 - t428 * t639;
t548 = t304 * t408 - t203;
t201 = qJD(1) * t305 - t428 * t640;
t550 = t306 * t408 + t201;
t16 = (t431 * t199 + (t136 + t721) * qJD(1)) * t431 + (t135 * qJD(1) + (-t200 * t401 + t202 * t402 - t305 * t630 - t307 * t631 + t597) * t428 + (-t198 + t548 * t402 + t550 * t401 + (-t302 - t498) * qJD(1)) * t431) * t428;
t21 = t101 * t635 + t103 * t338 + t105 * t339 + t192 * t215 + t193 * t217 + t213 * t457;
t22 = t100 * t635 + t102 * t338 + t104 * t339 + t192 * t216 + t193 * t218 + t214 * t457;
t72 = t213 * t635 + t215 * t338 + t217 * t339;
t73 = t214 * t635 + t216 * t338 + t218 * t339;
t527 = t428 * t72 + t431 * t73;
t9 = qJD(1) * t527 - t21 * t431 + t22 * t428;
t588 = -t14 - t16 - t9;
t587 = pkin(3) * t625;
t586 = pkin(4) * t630;
t585 = t431 * t683;
t583 = pkin(2) * t590;
t262 = -Icges(7,3) * t396 + t395 * t507;
t110 = t262 * t635 + t263 * t338 + t264 * t339;
t86 = -t213 * t396 + t395 * t506;
t578 = t86 / 0.2e1 + t110 / 0.2e1;
t111 = t262 * t634 + t263 * t340 + t264 * t341;
t87 = -t214 * t396 + t395 * t505;
t577 = t87 / 0.2e1 + t111 / 0.2e1;
t352 = t431 * t362;
t394 = t422 * t594;
t576 = t428 * (t394 + (t292 - t362) * t428 + t609 * t593 - t604) + t431 * (-qJD(1) * t226 - t352 + t613) + t226 * t593;
t162 = t518 * t633 + (Icges(7,5) * t408 + (-Icges(7,1) * t426 - t666) * qJD(6)) * t395;
t575 = t395 * t429 * t162 + t396 * t264 * t628 + t262 * t636;
t542 = t431 * t584;
t574 = t428 * (t362 * t428 + t593 * t607 - t394 + t603) + t431 * (-qJD(1) * t258 + t352 + t542) + t258 * t593;
t573 = t191 * rSges(7,1) + t190 * rSges(7,2) + rSges(7,3) * t579;
t572 = t428 * (qJD(1) * t313 - t428 * t479) + t431 * (-t627 * t679 + (-t401 * t627 - t402 * t594) * rSges(5,1) + t606) + t312 * t593;
t571 = -t227 + t716;
t477 = t368 * t421;
t570 = t428 * (qJD(1) * t329 - t428 * t477) + t431 * (-t622 * t681 + (-t411 * t622 - t412 * t594) * rSges(4,1) + t605) + t328 * t593;
t569 = -t259 + t616;
t568 = t428 * ((-t431 * t689 - t690) * qJD(1) - t603) + t431 * (-t542 + (t431 * t688 + t724) * qJD(1)) + t316 * t593;
t255 = t265 * t594;
t567 = t356 * t594 + t255 + t380;
t566 = t382 + t611;
t562 = t427 * t594;
t561 = t633 / 0.2e1;
t560 = t594 / 0.2e1;
t559 = t593 / 0.2e1;
t558 = -t368 - t695;
t557 = -t363 - t694;
t556 = -t355 - t693;
t555 = t278 * t408 + t180;
t553 = -t276 * t408 + t182;
t551 = t307 * t408 + t200;
t549 = -t305 * t408 + t202;
t230 = -qJD(1) * t323 - t431 * t638;
t547 = t326 * t421 + t230;
t231 = qJD(1) * t324 - t428 * t638;
t546 = t325 * t421 + t231;
t232 = -qJD(1) * t325 - t431 * t637;
t545 = -t324 * t421 + t232;
t233 = qJD(1) * t326 - t428 * t637;
t544 = t323 * t421 - t233;
t543 = -t413 * t428 + t330;
t531 = -rSges(6,2) * t395 + t684;
t286 = -rSges(6,3) * t431 + t428 * t531;
t96 = t428 * t286 + t431 * t287 + t617;
t539 = -t259 + t571;
t113 = t211 + t615;
t538 = t382 + t567;
t537 = -t693 + t715;
t536 = -t693 - t694;
t534 = -t683 + t687;
t530 = rSges(7,1) * t193 + rSges(7,2) * t192;
t48 = t428 * t73 - t431 * t72;
t74 = t213 * t634 + t215 * t340 + t217 * t341;
t75 = t214 * t634 + t216 * t340 + t218 * t341;
t49 = t428 * t75 - t431 * t74;
t526 = t428 * t74 + t431 * t75;
t525 = t87 * t428 - t86 * t431;
t524 = t428 * t86 + t431 * t87;
t522 = Icges(3,1) * t427 + t674;
t516 = Icges(3,2) * t430 + t675;
t365 = Icges(4,5) * t411 + Icges(4,6) * t412;
t504 = t221 * t431 - t222 * t428;
t490 = -t586 + t718;
t487 = -pkin(1) - t534;
t486 = t557 - t695;
t485 = -t355 + t536;
t177 = t537 * t431;
t484 = -t403 - t533;
t483 = -t374 - t532;
t482 = -t348 - t531;
t458 = -t396 * t594 - t580;
t481 = t428 * (-t428 * t478 + (t431 * t531 + t415) * qJD(1)) + t431 * (rSges(6,1) * t458 - rSges(6,2) * t579 + t608) + t286 * t593 + t576;
t480 = t572 + t574;
t318 = t535 * t428;
t60 = t428 * t318 + t431 * t319 + t617 - t717;
t63 = t96 + t615;
t475 = -t586 - t587;
t474 = t536 + t715;
t465 = t421 * t365;
t464 = -t583 - t587;
t463 = qJD(2) * t522;
t462 = qJD(2) * t516;
t461 = qJD(2) * (-Icges(3,5) * t427 - Icges(3,6) * t430);
t460 = t536 - t695;
t459 = t395 * t699 - t348 - t691;
t261 = t486 * t431;
t244 = t485 * t431;
t19 = t101 * t634 + t103 * t340 + t105 * t341 + t190 * t215 + t191 * t217 + t213 * t456;
t20 = t100 * t634 + t102 * t340 + t104 * t341 + t190 * t216 + t191 * t218 + t214 * t456;
t8 = qJD(1) * t526 - t19 * t431 + t20 * t428;
t455 = t48 * t594 + t49 * t593 + ((-t122 - t135) * t594 + t739 * t593) * t431 + (t8 + ((t178 + t198) * t428 + (t725 + t726 - t739) * qJD(1)) * t428 + (t123 + t136) * t594 + t736 * t593 + ((-t555 * t395 + t553 * t396 - t551 * t401 + t549 * t402 - t179 - t199) * t428 + (t181 * t395 - t183 * t396 + t201 * t401 - t203 * t402 + t275 * t633 + t277 * t636 + t304 * t630 + t306 * t631 - t713 - t714) * t431 + ((-t501 - t499 + t741) * t428 + t736) * qJD(1)) * t431) * t428;
t282 = t531 * t408;
t454 = -t282 + t475;
t453 = -t311 + t464;
t452 = -t355 + t460;
t55 = t60 + t615;
t157 = t474 * t431;
t451 = t475 + t718;
t108 = -rSges(7,3) * t565 + t573;
t109 = rSges(7,3) * t457 + t530;
t358 = pkin(10) * t579;
t450 = t576 + (t221 + t318) * t593 + (pkin(5) * t458 - pkin(10) * t565 + t108 + t358) * t431 + (t109 + t457 * pkin(10) + (t396 * t593 - t581) * pkin(5)) * t428;
t449 = t460 + t715;
t448 = t481 + t574;
t238 = t452 * t431;
t37 = t160 * t634 + t161 * t340 + t162 * t341 + t190 * t263 + t191 * t264 + t262 * t456;
t3 = (t408 * t526 - t37) * t396 + (-qJD(1) * t49 + t111 * t408 + t19 * t428 + t20 * t431) * t395;
t32 = -t110 * t396 + t395 * t527;
t33 = -t111 * t396 + t395 * t526;
t38 = t160 * t635 + t161 * t338 + t162 * t339 + t192 * t263 + t193 * t264 + t262 * t457;
t4 = (t408 * t527 - t38) * t396 + (-qJD(1) * t48 + t110 * t408 + t21 * t428 + t22 * t431) * t395;
t447 = t3 * t701 + t4 * t700 - t396 * (qJD(1) * t524 + t676 - t677) / 0.2e1 + t32 * t560 + t33 * t559 + t9 * t635 / 0.2e1 + t525 * t636 / 0.2e1 + t8 * t634 / 0.2e1 + (t431 * t561 - t565 / 0.2e1) * t49 + (t395 * t559 + t428 * t561) * t48;
t446 = t464 - t586;
t445 = rSges(3,2) * t562 + rSges(3,3) * t593 - t431 * t476;
t139 = -t321 * t431 - t428 * t497;
t140 = -t322 * t431 - t727;
t141 = t321 * t428 - t722;
t142 = t322 * t428 - t431 * t496;
t228 = -t431 * t465 - t712;
t229 = -t428 * t465 + t596;
t444 = t428 * ((t428 * t228 + (t141 + t727) * qJD(1)) * t428 + (t142 * qJD(1) + (t231 * t411 - t233 * t412 + t323 * t625 + t325 * t626 - t712) * t431 + (-t229 + t545 * t412 - t547 * t411 + (t322 - t497) * qJD(1)) * t428) * t431) + (-t139 * t431 + t140 * t428) * t594 + (-t141 * t431 + t142 * t428) * t593 + t455;
t149 = t449 * t431;
t443 = -t282 + t446;
t442 = t428 * t459 - t624;
t441 = t450 + t574;
t440 = t446 + t718;
t439 = t431 * t588 + t455;
t334 = t515 * t421;
t335 = t521 * t421;
t436 = qJD(1) * t365 + (t335 - t638) * t412 + (-t334 - t637) * t411;
t24 = (t431 * t229 + (t140 + t722) * qJD(1)) * t431 + (t139 * qJD(1) + (-t230 * t411 + t232 * t412 - t324 * t625 - t326 * t626 + t596) * t428 + (-t228 + t544 * t412 + t546 * t411 + (-t321 - t496) * qJD(1)) * t431) * t428;
t435 = t444 + (-t24 + t588) * t431;
t434 = -t677 / 0.2e1 + t676 / 0.2e1 + (t395 * t553 + t396 * t555 + t401 * t549 + t402 * t551 + t428 * t730 + t431 * t734 + t37) * t701 + (-t395 * t552 + t396 * t554 - t401 * t548 + t402 * t550 + t428 * t734 - t431 * t730 + t38) * t700 + (t275 * t396 + t277 * t395 + t304 * t402 + t306 * t401 - t428 * t740 - t431 * t735 + t110 + t86) * t560 + (t276 * t396 + t278 * t395 + t305 * t402 + t307 * t401 + t428 * t735 - t431 * t740 + t111 + t87) * t559;
t433 = t434 + (t411 * t545 + t412 * t547 + t428 * t707 + t436 * t431) * t701 + (-t411 * t544 + t412 * t546 + t436 * t428 - t431 * t707) * t700 + (t323 * t412 + t325 * t411 - t365 * t431 - t428 * t491) * t560 + (t324 * t412 + t326 * t411 + t365 * t428 - t431 * t491) * t559;
t393 = pkin(2) * t562;
t373 = t534 * qJD(2);
t354 = -t585 + t602;
t353 = t428 * t534 - t678;
t315 = t558 * t431;
t314 = t558 * t428;
t295 = t690 + (pkin(1) - t683) * t431 + t602;
t294 = t428 * t487 + t420 + t678;
t289 = t557 * t431;
t288 = t557 * t428;
t268 = t556 * t431;
t267 = t556 * t428;
t260 = t486 * t428;
t257 = -t428 * t432 + t329 + t389;
t256 = (rSges(4,3) - t432) * t431 + t484 * t428;
t248 = t428 * t461 + t595;
t247 = -qJD(1) * t342 + t431 * t461;
t243 = t485 * t428;
t240 = -t428 * t422 + t313 + t364;
t239 = (rSges(5,3) - t422) * t431 + t483 * t428;
t237 = t452 * t428;
t235 = t729 + ((-rSges(3,3) - pkin(7)) * t428 + t487 * t431) * qJD(1);
t234 = (t420 + (-pkin(1) - t687) * t428) * qJD(1) + t445;
t210 = -t368 * t593 - t647 + (-t427 * t593 - t428 * t590) * pkin(2);
t209 = t368 * t594 + t393 + (-t337 - t583) * t431;
t208 = t287 + t543;
t207 = (rSges(6,3) - t413) * t431 + t482 * t428;
t176 = t537 * t428;
t171 = -t363 * t593 - t649 + (-t411 * t593 - t412 * t623) * pkin(3);
t170 = (-t311 - t587) * t431 + t610;
t166 = t343 * t428 - t431 * t494;
t165 = t342 * t428 - t723;
t164 = -t343 * t431 - t728;
t163 = -t342 * t431 - t428 * t495;
t156 = t474 * t428;
t155 = t368 * t623 + (t431 * t484 - t417) * qJD(1) + t603;
t154 = (-t403 - t686) * t594 + (-t477 - t732) * t431 + t605;
t151 = -t355 * t593 - t282 * t428 + (-t401 * t593 - t402 * t629) * pkin(4);
t150 = (-t282 - t586) * t431 + t611;
t148 = t449 * t428;
t146 = qJD(1) * t261 + t428 * t453;
t145 = t431 * t453 + t393 + t610;
t134 = t394 + (-t362 + t479) * t428 + (t431 * t483 - t416) * qJD(1);
t133 = t352 - t431 * t479 + (-t431 * t422 + (-t374 - t685) * t428) * qJD(1) + t606;
t132 = qJD(1) * t244 + t428 * t454;
t131 = t431 * t454 + t566;
t130 = t225 + t612;
t129 = -t222 * t396 - t265 * t634;
t128 = t221 * t396 + t265 * t635;
t127 = t543 - t716;
t126 = t442 + t529;
t119 = qJD(1) * t238 + t428 * t443;
t118 = t431 * t443 + t393 + t566;
t117 = -t262 * t396 + (t264 * t429 - t652) * t395;
t116 = t504 * t395;
t115 = (-t292 + t478) * t428 + (t431 * t482 - t415) * qJD(1) + t604;
t114 = -t431 * t478 + (-t624 + (-t348 - t684) * t428) * qJD(1) + t608 + t613;
t112 = t117 * t636;
t97 = -t329 * t594 + t570;
t93 = -t313 * t594 + t572;
t90 = t113 + t612;
t89 = qJD(1) * t177 + t428 * t490;
t88 = t431 * t490 + t567;
t77 = qJD(1) * t157 + t428 * t451;
t76 = t431 * t451 + t538;
t71 = qJD(1) * t149 + t428 * t440;
t70 = t431 * t440 + t393 + t538;
t59 = (-t292 + (t396 * t699 + t692) * t408) * t428 + t459 * t593 - t530 + t604;
t58 = -pkin(5) * t580 + qJD(1) * t442 + t358 + t573 + t613;
t57 = (-t317 - t329) * t594 + t568 + t570;
t56 = t63 + t612;
t54 = (t265 * t629 + t109) * t396 + (t167 * t428 - t221 * t408 + t265 * t593) * t395;
t53 = (-t265 * t627 - t108) * t396 + (-t167 * t431 + t222 * t408 + t255) * t395;
t52 = t594 * t614 + t480;
t50 = t55 + t612;
t41 = t395 * t731 + t719 * t396 + t575;
t40 = t594 * t616 + t481;
t39 = (-t317 + t614) * t594 + t480 + t568;
t34 = t504 * t633 + (qJD(1) * t717 - t108 * t428 + t109 * t431) * t395;
t29 = t569 * t594 + t448;
t18 = (-t317 + t569) * t594 + t448 + t568;
t17 = t571 * t594 + t450;
t12 = t539 * t594 + t441;
t11 = (-t317 + t539) * t594 + t441 + t568;
t1 = [t575 + (t126 * t59 + t127 * t58) * t702 + (t114 * t208 + t115 * t207) * t703 + (t133 * t240 + t134 * t239) * t704 + (t154 * t257 + t155 * t256) * t705 + (t234 * t295 + t235 * t294) * t706 + t351 * t633 - t350 * t636 - t360 * t631 + t412 * t334 - t366 * t626 + t402 * t309 + t411 * t335 + t401 * t310 + t367 * t625 + t361 * t630 + (-t516 + t523) * t591 + (t517 + t522) * t590 + (t280 + t719) * t396 + (t281 + t731) * t395; (t423 / 0.2e1 + t424 / 0.2e1) * t511 * qJD(2) + (-qJD(2) * t495 + (qJD(1) * t345 - t428 * t462) * t430 + (qJD(1) * t347 - t428 * t463) * t427) * t700 + (-qJD(2) * t494 + (-qJD(1) * t344 - t431 * t462) * t430 + (-qJD(1) * t346 - t431 * t463) * t427) * t701 + t433 + ((t645 / 0.2e1 + t643 / 0.2e1 - t295 * t698) * t431 + (t294 * t698 + t646 / 0.2e1 + t644 / 0.2e1) * t428) * qJD(1) + m(3) * ((-t234 * t428 - t235 * t431) * t387 + (-t294 * t431 - t295 * t428) * t373) + m(4) * (t154 * t314 + t155 * t315 + t209 * t256 + t210 * t257) + m(5) * (t133 * t260 + t134 * t261 + t145 * t239 + t146 * t240) + m(6) * (t114 * t237 + t115 * t238 + t118 * t207 + t119 * t208) + m(7) * (t126 * t70 + t127 * t71 + t148 * t58 + t149 * t59); t444 + t428 * ((t428 * t247 + (t165 + t728) * qJD(1)) * t428 + (t166 * qJD(1) + (t344 * t590 + t346 * t591) * t431 + (-t248 + (-t643 - t645) * qJD(2) + (t343 - t495) * qJD(1)) * t428) * t431) - t431 * ((t431 * t248 + (t164 + t723) * qJD(1)) * t431 + (t163 * qJD(1) + (-t345 * t590 - t347 * t591 + t595) * t428 + (-t247 + (t644 + t646) * qJD(2) - t494 * qJD(1)) * t431) * t428) + ((t353 * t428 + t354 * t431) * ((qJD(1) * t353 + t445) * t431 + (-t729 + (-t354 - t585 + t418) * qJD(1)) * t428) + t599 * t387 * t373) * t706 + (-t163 * t431 + t164 * t428) * t594 + (-t165 * t431 + t166 * t428) * t593 - t431 * t9 - t431 * t16 - t431 * t14 - t431 * t24 + (t118 * t238 + t119 * t237 + t18 * t56) * t703 + (t11 * t50 + t148 * t71 + t149 * t70) * t702 + (t130 * t57 + t209 * t315 + t210 * t314) * t705 + (t145 * t261 + t146 * t260 + t39 * t90) * t704; t433 + m(5) * (t133 * t288 + t134 * t289 + t170 * t239 + t171 * t240) + m(6) * (t114 * t243 + t115 * t244 + t131 * t207 + t132 * t208) + m(7) * (t126 * t76 + t127 * t77 + t156 * t58 + t157 * t59) + m(4) * (-t256 * t431 - t257 * t428) * t337 + (-t154 * t428 - t155 * t431 + (t256 * t428 - t257 * t431) * qJD(1)) * t697; t435 + (-t209 * t431 - t210 * t428 + (-t314 * t431 + t315 * t428) * qJD(1)) * t697 + m(4) * (-t315 * t337 * t431 + t130 * t97 + t225 * t57 - t314 * t647) + m(7) * (t11 * t55 + t12 * t50 + t148 * t77 + t149 * t76 + t156 * t71 + t157 * t70) + m(6) * (t118 * t244 + t119 * t243 + t131 * t238 + t132 * t237 + t18 * t63 + t29 * t56) + m(5) * (t113 * t39 + t145 * t289 + t146 * t288 + t170 * t261 + t171 * t260 + t52 * t90); (t12 * t55 + t156 * t77 + t157 * t76) * t702 + (t131 * t244 + t132 * t243 + t29 * t63) * t703 + (t113 * t52 + t170 * t289 + t171 * t288) * t704 + t435 + (t337 * t368 * t599 + t225 * t97) * t705; (-t133 * t428 - t134 * t431 + (t239 * t428 - t240 * t431) * qJD(1)) * t696 + m(5) * (-t239 * t431 - t240 * t428) * t311 + m(6) * (t114 * t267 + t115 * t268 + t150 * t207 + t151 * t208) + m(7) * (t126 * t88 + t127 * t89 + t176 * t58 + t177 * t59) + t434; (-t145 * t431 - t146 * t428 + (-t260 * t431 + t261 * t428) * qJD(1)) * t696 + m(5) * (t211 * t39 - t260 * t649 - t261 * t648 + t90 * t93) + m(7) * (t11 * t60 + t148 * t89 + t149 * t88 + t17 * t50 + t176 * t71 + t177 * t70) + m(6) * (t118 * t268 + t119 * t267 + t150 * t238 + t151 * t237 + t18 * t96 + t40 * t56) + t439; m(5) * (t113 * t93 + t211 * t52 - t288 * t649 - t289 * t648) + m(7) * (t12 * t60 + t156 * t89 + t157 * t88 + t17 * t55 + t176 * t77 + t177 * t76) + m(6) * (t131 * t268 + t132 * t267 + t150 * t244 + t151 * t243 + t29 * t96 + t40 * t63) + t439 + (-t170 * t431 - t171 * t428 + (-t288 * t431 + t289 * t428) * qJD(1)) * t696; (t17 * t60 + t176 * t89 + t177 * t88) * t702 + (t150 * t268 + t151 * t267 + t40 * t96) * t703 + (t311 * t363 * t599 + t211 * t93) * t704 + t439; m(7) * (t428 * t59 - t431 * t58 + (t126 * t431 + t127 * t428) * qJD(1)) + m(6) * (-t114 * t431 + t115 * t428 + (t207 * t431 + t208 * t428) * qJD(1)); m(7) * (t428 * t70 - t431 * t71 + (t148 * t428 + t149 * t431) * qJD(1)) + m(6) * (t118 * t428 - t119 * t431 + (t237 * t428 + t238 * t431) * qJD(1)); m(7) * (t428 * t76 - t431 * t77 + (t156 * t428 + t157 * t431) * qJD(1)) + m(6) * (t131 * t428 - t132 * t431 + (t243 * t428 + t244 * t431) * qJD(1)); m(7) * (t428 * t88 - t431 * t89 + (t176 * t428 + t177 * t431) * qJD(1)) + m(6) * (t150 * t428 - t151 * t431 + (t267 * t428 + t268 * t431) * qJD(1)); 0; m(7) * (t126 * t54 + t127 * t53 + t128 * t59 + t129 * t58) + t112 + (-t41 + (t428 * t578 + t431 * t577) * t408) * t396 + ((t28 / 0.2e1 + t37 / 0.2e1) * t431 + (t27 / 0.2e1 + t38 / 0.2e1) * t428 + (-t428 * t577 + t431 * t578) * qJD(1)) * t395; m(7) * (t11 * t116 + t128 * t70 + t129 * t71 + t148 * t53 + t149 * t54 + t34 * t50) + t447; m(7) * (t116 * t12 + t128 * t76 + t129 * t77 + t156 * t53 + t157 * t54 + t34 * t55) + t447; m(7) * (t116 * t17 + t128 * t88 + t129 * t89 + t176 * t53 + t177 * t54 + t34 * t60) + t447; m(7) * (t428 * t54 - t431 * t53 + (t128 * t431 + t129 * t428) * qJD(1)); (t116 * t34 + t128 * t54 + t129 * t53) * t702 + (t41 * t396 - t112 + (t428 * t32 + t431 * t33 - t396 * t524) * t408) * t396 + (t431 * t3 + t428 * t4 + t524 * t636 + (-t117 * t408 - t27 * t428 - t28 * t431) * t396 + (t431 * t32 - t428 * t33 + t396 * t525) * qJD(1)) * t395;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
