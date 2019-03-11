% Calculate time derivative of joint inertia matrix for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:12
% EndTime: 2019-03-09 09:17:56
% DurationCPUTime: 27.62s
% Computational Cost: add. (46259->1238), mult. (128972->1682), div. (0->0), fcn. (142129->10), ass. (0->500)
t459 = cos(pkin(6));
t465 = cos(qJ(2));
t466 = cos(qJ(1));
t566 = t465 * t466;
t533 = t459 * t566;
t462 = sin(qJ(2));
t463 = sin(qJ(1));
t569 = t462 * t463;
t430 = -t533 + t569;
t567 = t463 * t465;
t568 = t462 * t466;
t431 = t459 * t568 + t567;
t458 = sin(pkin(6));
t572 = t458 * t466;
t306 = Icges(5,5) * t430 - Icges(5,6) * t431 + Icges(5,3) * t572;
t310 = Icges(3,5) * t431 - Icges(3,6) * t430 - Icges(3,3) * t572;
t314 = Icges(4,4) * t431 - Icges(4,2) * t572 + Icges(4,6) * t430;
t623 = t306 - t310 - t314;
t432 = t459 * t567 + t568;
t531 = t459 * t569;
t433 = -t531 + t566;
t574 = t458 * t463;
t307 = Icges(5,5) * t432 - Icges(5,6) * t433 - Icges(5,3) * t574;
t311 = Icges(3,5) * t433 - Icges(3,6) * t432 + Icges(3,3) * t574;
t315 = Icges(4,4) * t433 + Icges(4,2) * t574 + Icges(4,6) * t432;
t622 = t307 - t311 - t315;
t308 = Icges(4,5) * t431 - Icges(4,6) * t572 + Icges(4,3) * t430;
t316 = Icges(3,4) * t431 - Icges(3,2) * t430 - Icges(3,6) * t572;
t318 = Icges(5,1) * t430 - Icges(5,4) * t431 + Icges(5,5) * t572;
t621 = t308 - t316 + t318;
t309 = Icges(4,5) * t433 + Icges(4,6) * t574 + Icges(4,3) * t432;
t317 = Icges(3,4) * t433 - Icges(3,2) * t432 + Icges(3,6) * t574;
t319 = Icges(5,1) * t432 - Icges(5,4) * t433 - Icges(5,5) * t574;
t620 = t317 - t309 - t319;
t312 = Icges(5,4) * t430 - Icges(5,2) * t431 + Icges(5,6) * t572;
t320 = Icges(4,1) * t431 - Icges(4,4) * t572 + Icges(4,5) * t430;
t322 = Icges(3,1) * t431 - Icges(3,4) * t430 - Icges(3,5) * t572;
t619 = t322 + t320 - t312;
t313 = Icges(5,4) * t432 - Icges(5,2) * t433 - Icges(5,6) * t574;
t321 = Icges(4,1) * t433 + Icges(4,4) * t574 + Icges(4,5) * t432;
t323 = Icges(3,1) * t433 - Icges(3,4) * t432 + Icges(3,5) * t574;
t618 = -t323 - t321 + t313;
t577 = Icges(4,5) * t462;
t390 = Icges(4,6) * t459 + (-Icges(4,3) * t465 + t577) * t458;
t581 = Icges(3,4) * t462;
t394 = Icges(3,6) * t459 + (Icges(3,2) * t465 + t581) * t458;
t579 = Icges(5,4) * t462;
t395 = -Icges(5,5) * t459 + (-Icges(5,1) * t465 - t579) * t458;
t576 = Icges(4,5) * t465;
t396 = Icges(4,4) * t459 + (Icges(4,1) * t462 - t576) * t458;
t580 = Icges(3,4) * t465;
t397 = Icges(3,5) * t459 + (Icges(3,1) * t462 + t580) * t458;
t544 = qJD(2) * t458;
t402 = (Icges(5,5) * t462 - Icges(5,6) * t465) * t544;
t403 = (Icges(4,3) * t462 + t576) * t544;
t404 = (Icges(3,5) * t465 - Icges(3,6) * t462) * t544;
t406 = (Icges(4,4) * t465 + Icges(4,6) * t462) * t544;
t407 = (-Icges(3,2) * t462 + t580) * t544;
t409 = (Icges(4,1) * t465 + t577) * t544;
t410 = (Icges(3,1) * t465 - t581) * t544;
t542 = qJD(2) * t465;
t512 = t458 * t542;
t543 = qJD(2) * t462;
t513 = t458 * t543;
t573 = t458 * t465;
t575 = t458 * t462;
t634 = (t409 + t410) * t575 + (-t403 + t407) * t573 + (t396 + t397) * t512 + (t390 - t394 + t395) * t513 + (-t402 + t404 + t406) * t459;
t545 = qJD(1) * t466;
t514 = t458 * t545;
t546 = qJD(1) * t463;
t515 = t458 * t546;
t457 = t466 * pkin(1);
t547 = pkin(8) * t574 + t457;
t416 = t430 * qJ(3);
t455 = pkin(8) * t572;
t633 = t455 - t416;
t461 = sin(qJ(5));
t593 = cos(qJ(5));
t516 = t458 * t593;
t429 = -t459 * t461 - t465 * t516;
t374 = t432 * t593 - t461 * t574;
t632 = -t430 * t621 - t431 * t619 - t572 * t623;
t631 = -t430 * t620 - t431 * t618 + t572 * t622;
t630 = -t432 * t620 - t433 * t618 - t574 * t622;
t629 = t432 * t621 + t433 * t619 - t574 * t623;
t578 = Icges(5,4) * t465;
t392 = -Icges(5,6) * t459 + (-Icges(5,2) * t462 - t578) * t458;
t408 = (Icges(5,1) * t462 - t578) * t544;
t405 = (-Icges(5,2) * t465 + t579) * t544;
t570 = t462 * t405;
t628 = ((-t570 + (-qJD(2) * t392 - t408) * t465) * t458 + t634) * t459;
t506 = qJD(2) * t459 + qJD(1);
t350 = -qJD(1) * t533 - t466 * t542 + t506 * t569;
t351 = qJD(1) * t431 + qJD(2) * t432;
t227 = -Icges(4,5) * t351 + Icges(4,6) * t514 - Icges(4,3) * t350;
t235 = -Icges(3,4) * t351 + Icges(3,2) * t350 + Icges(3,6) * t514;
t237 = -Icges(5,1) * t350 + Icges(5,4) * t351 - Icges(5,5) * t514;
t627 = -t237 + t235 - t227;
t352 = qJD(1) * t432 + qJD(2) * t431;
t353 = -qJD(1) * t531 - t463 * t543 + t506 * t566;
t228 = Icges(4,5) * t353 + Icges(4,6) * t515 + Icges(4,3) * t352;
t236 = Icges(3,4) * t353 - Icges(3,2) * t352 + Icges(3,6) * t515;
t238 = Icges(5,1) * t352 - Icges(5,4) * t353 - Icges(5,5) * t515;
t626 = t238 - t236 + t228;
t231 = -Icges(5,4) * t350 + Icges(5,2) * t351 - Icges(5,6) * t514;
t239 = -Icges(4,1) * t351 + Icges(4,4) * t514 - Icges(4,5) * t350;
t241 = -Icges(3,1) * t351 + Icges(3,4) * t350 + Icges(3,5) * t514;
t625 = t241 + t239 - t231;
t232 = Icges(5,4) * t352 - Icges(5,2) * t353 - Icges(5,6) * t515;
t240 = Icges(4,1) * t353 + Icges(4,4) * t515 + Icges(4,5) * t352;
t242 = Icges(3,1) * t353 - Icges(3,4) * t352 + Icges(3,5) * t515;
t624 = -t242 - t240 + t232;
t225 = -Icges(5,5) * t350 + Icges(5,6) * t351 - Icges(5,3) * t514;
t226 = Icges(5,5) * t352 - Icges(5,6) * t353 - Icges(5,3) * t515;
t229 = -Icges(3,5) * t351 + Icges(3,6) * t350 + Icges(3,3) * t514;
t230 = Icges(3,5) * t353 - Icges(3,6) * t352 + Icges(3,3) * t515;
t233 = -Icges(4,4) * t351 + Icges(4,2) * t514 - Icges(4,6) * t350;
t234 = Icges(4,4) * t353 + Icges(4,2) * t515 + Icges(4,6) * t352;
t617 = (-t225 + t229 + t233) * t463 + (-t234 - t230 + t226) * t466;
t616 = 2 * m(3);
t615 = 2 * m(4);
t614 = 2 * m(5);
t613 = 2 * m(6);
t612 = 2 * m(7);
t611 = t458 ^ 2;
t610 = 0.2e1 * t458;
t609 = m(5) / 0.2e1;
t608 = -m(6) / 0.2e1;
t607 = m(6) / 0.2e1;
t606 = -m(7) / 0.2e1;
t605 = m(7) / 0.2e1;
t604 = -pkin(2) - pkin(3);
t501 = qJD(1) * t516;
t255 = qJD(5) * t374 - t350 * t461 + t466 * t501;
t603 = t255 / 0.2e1;
t372 = t430 * t593 + t461 * t572;
t257 = qJD(5) * t372 + t352 * t461 + t463 * t501;
t602 = t257 / 0.2e1;
t366 = qJD(5) * t429 + t461 * t513;
t601 = t366 / 0.2e1;
t476 = -t430 * t461 + t466 * t516;
t600 = -t476 / 0.2e1;
t477 = -t432 * t461 - t463 * t516;
t599 = -t477 / 0.2e1;
t428 = -t459 * t593 + t461 * t573;
t598 = -t428 / 0.2e1;
t597 = t459 / 0.2e1;
t596 = t463 / 0.2e1;
t595 = -rSges(4,1) - pkin(2);
t594 = -rSges(7,3) - pkin(10);
t258 = qJD(5) * t476 + t352 * t593 - t461 * t515;
t592 = t258 * pkin(5);
t591 = t372 * pkin(5);
t590 = t463 * pkin(1);
t589 = t352 * rSges(5,1);
t588 = t352 * rSges(4,3);
t587 = t430 * rSges(5,1);
t460 = sin(qJ(6));
t464 = cos(qJ(6));
t369 = -t429 * t460 + t464 * t575;
t370 = t429 * t464 + t460 * t575;
t259 = Icges(7,5) * t370 + Icges(7,6) * t369 - Icges(7,3) * t428;
t260 = Icges(7,4) * t370 + Icges(7,2) * t369 - Icges(7,6) * t428;
t261 = Icges(7,1) * t370 + Icges(7,4) * t369 - Icges(7,5) * t428;
t107 = -t259 * t428 + t260 * t369 + t261 * t370;
t367 = qJD(5) * t428 + t513 * t593;
t275 = -qJD(6) * t370 - t367 * t460 + t464 * t512;
t276 = qJD(6) * t369 + t367 * t464 + t460 * t512;
t147 = Icges(7,5) * t276 + Icges(7,6) * t275 + Icges(7,3) * t366;
t148 = Icges(7,4) * t276 + Icges(7,2) * t275 + Icges(7,6) * t366;
t149 = Icges(7,1) * t276 + Icges(7,4) * t275 + Icges(7,5) * t366;
t48 = -t147 * t428 + t148 * t369 + t149 * t370 + t259 * t366 + t260 * t275 + t261 * t276;
t586 = t107 * t366 - t428 * t48;
t585 = t107 * t512 + t48 * t575;
t302 = Icges(6,5) * t429 + Icges(6,6) * t428 + Icges(6,3) * t575;
t303 = Icges(6,4) * t429 + Icges(6,2) * t428 + Icges(6,6) * t575;
t304 = Icges(6,1) * t429 + Icges(6,4) * t428 + Icges(6,5) * t575;
t163 = t302 * t575 + t303 * t428 + t304 * t429;
t281 = Icges(6,5) * t367 - Icges(6,6) * t366 + Icges(6,3) * t512;
t282 = Icges(6,4) * t367 - Icges(6,2) * t366 + Icges(6,6) * t512;
t283 = Icges(6,1) * t367 - Icges(6,4) * t366 + Icges(6,5) * t512;
t80 = t281 * t575 + t282 * t428 + t283 * t429 + t302 * t512 - t303 * t366 + t304 * t367;
t584 = t163 * t512 + t575 * t80;
t256 = qJD(5) * t477 - t350 * t593 - t461 * t514;
t301 = t374 * t464 + t433 * t460;
t159 = -qJD(6) * t301 - t256 * t460 - t351 * t464;
t300 = -t374 * t460 + t433 * t464;
t160 = qJD(6) * t300 + t256 * t464 - t351 * t460;
t87 = rSges(7,1) * t160 + rSges(7,2) * t159 + rSges(7,3) * t255;
t583 = pkin(5) * t256 + pkin(10) * t255 + t87;
t299 = t372 * t464 + t431 * t460;
t161 = -qJD(6) * t299 - t258 * t460 + t353 * t464;
t298 = -t372 * t460 + t431 * t464;
t162 = qJD(6) * t298 + t258 * t464 + t353 * t460;
t495 = -t162 * rSges(7,1) - t161 * rSges(7,2);
t88 = rSges(7,3) * t257 - t495;
t582 = pkin(10) * t257 + t592 + t88;
t150 = rSges(7,1) * t276 + rSges(7,2) * t275 + rSges(7,3) * t366;
t565 = pkin(5) * t367 + pkin(10) * t366 + t150;
t494 = -rSges(7,1) * t299 - rSges(7,2) * t298;
t197 = -rSges(7,3) * t476 - t494;
t564 = -pkin(10) * t476 + t197 + t591;
t198 = rSges(7,1) * t301 + rSges(7,2) * t300 - rSges(7,3) * t477;
t563 = pkin(5) * t374 - pkin(10) * t477 + t198;
t219 = -pkin(2) * t351 - t350 * qJ(3) + t432 * qJD(3);
t218 = t459 * t219;
t347 = t351 * pkin(3);
t541 = qJD(4) * t463;
t291 = -t347 + (-qJ(4) * t545 - t541) * t458;
t562 = t291 * t459 + t218;
t557 = -qJ(3) * t352 - qJD(3) * t430;
t220 = t353 * pkin(2) - t557;
t499 = qJ(4) * t515 - qJD(4) * t572;
t292 = t353 * pkin(3) - t499;
t561 = -t220 - t292;
t262 = rSges(7,1) * t370 + rSges(7,2) * t369 - rSges(7,3) * t428;
t560 = pkin(5) * t429 - pkin(10) * t428 + t262;
t355 = t431 * pkin(2) + t416;
t356 = pkin(2) * t433 + qJ(3) * t432;
t559 = t355 * t574 + t356 * t572;
t336 = t459 * t356;
t426 = t433 * pkin(3);
t535 = qJ(4) * t574;
t381 = t426 - t535;
t558 = t381 * t459 + t336;
t556 = -rSges(5,1) * t350 + rSges(5,2) * t351;
t534 = qJ(4) * t572;
t380 = t431 * pkin(3) + t534;
t555 = -t355 - t380;
t554 = -t356 - t381;
t379 = (-qJD(3) * t465 + (pkin(2) * t465 + qJ(3) * t462) * qJD(2)) * t458;
t553 = -t379 - (rSges(4,1) * t465 + rSges(4,3) * t462) * t544;
t552 = -pkin(3) * t512 + qJD(4) * t459 - t379;
t434 = (pkin(2) * t462 - qJ(3) * t465) * t458;
t385 = t434 * t515;
t437 = pkin(3) * t575 - qJ(4) * t459;
t551 = t437 * t515 + t385;
t399 = rSges(4,2) * t459 + (rSges(4,1) * t462 - rSges(4,3) * t465) * t458;
t550 = -t399 - t434;
t549 = rSges(5,1) * t432 - rSges(5,2) * t433;
t548 = -t434 - t437;
t540 = rSges(5,2) + t604;
t539 = -rSges(6,3) + t604;
t538 = pkin(1) * t546;
t191 = Icges(7,5) * t299 + Icges(7,6) * t298 - Icges(7,3) * t476;
t193 = Icges(7,4) * t299 + Icges(7,2) * t298 - Icges(7,6) * t476;
t195 = Icges(7,1) * t299 + Icges(7,4) * t298 - Icges(7,5) * t476;
t82 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t257;
t84 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t257;
t86 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t257;
t21 = t191 * t366 + t193 * t275 + t195 * t276 + t369 * t84 + t370 * t86 - t428 * t82;
t35 = -t147 * t476 + t148 * t298 + t149 * t299 + t161 * t260 + t162 * t261 + t257 * t259;
t537 = t21 / 0.2e1 + t35 / 0.2e1;
t192 = Icges(7,5) * t301 + Icges(7,6) * t300 - Icges(7,3) * t477;
t194 = Icges(7,4) * t301 + Icges(7,2) * t300 - Icges(7,6) * t477;
t196 = Icges(7,1) * t301 + Icges(7,4) * t300 - Icges(7,5) * t477;
t81 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t255;
t83 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t255;
t85 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t255;
t22 = t192 * t366 + t194 * t275 + t196 * t276 + t369 * t83 + t370 * t85 - t428 * t81;
t34 = -t147 * t477 + t148 * t300 + t149 * t301 + t159 * t260 + t160 * t261 + t255 * t259;
t536 = t34 / 0.2e1 + t22 / 0.2e1;
t101 = -t259 * t477 + t260 * t300 + t261 * t301;
t74 = -t192 * t428 + t194 * t369 + t196 * t370;
t530 = t101 / 0.2e1 + t74 / 0.2e1;
t100 = -t259 * t476 + t260 * t298 + t261 * t299;
t73 = -t191 * t428 + t193 * t369 + t195 * t370;
t529 = t73 / 0.2e1 + t100 / 0.2e1;
t528 = t219 * t572 + t220 * t574 + t355 * t514;
t273 = -pkin(4) * t350 - t351 * pkin(9);
t527 = t273 * t459 + t562;
t274 = t352 * pkin(4) + pkin(9) * t353;
t526 = -t274 + t561;
t141 = rSges(6,1) * t256 - rSges(6,2) * t255 - rSges(6,3) * t351;
t358 = pkin(4) * t432 + pkin(9) * t433;
t525 = t358 * t459 + t558;
t244 = -rSges(4,1) * t351 + rSges(4,2) * t514 - rSges(4,3) * t350;
t245 = -rSges(3,1) * t351 + rSges(3,2) * t350 + rSges(3,3) * t514;
t357 = pkin(4) * t430 + pkin(9) * t431;
t524 = -t357 + t555;
t523 = -t358 + t554;
t271 = rSges(6,1) * t374 + rSges(6,2) * t477 + rSges(6,3) * t433;
t522 = -(rSges(5,1) * t462 - rSges(5,2) * t465) * t544 + t552;
t521 = -(pkin(4) * t462 + pkin(9) * t465) * t544 + t552;
t435 = (-pkin(4) * t465 + pkin(9) * t462) * t458;
t520 = t435 * t515 + t551;
t401 = -rSges(5,3) * t459 + (-rSges(5,1) * t465 - rSges(5,2) * t462) * t458;
t519 = -t401 + t548;
t329 = rSges(4,1) * t433 + rSges(4,2) * t574 + rSges(4,3) * t432;
t330 = rSges(3,1) * t433 - rSges(3,2) * t432 + rSges(3,3) * t574;
t518 = -t435 + t548;
t511 = t455 - t590;
t510 = t458 * (-rSges(5,3) - qJ(4));
t509 = t466 * t550;
t507 = t609 + t607 + t605;
t284 = rSges(6,1) * t367 - rSges(6,2) * t366 + rSges(6,3) * t512;
t505 = -t284 + t521;
t305 = rSges(6,1) * t429 + rSges(6,2) * t428 + rSges(6,3) * t575;
t504 = -t305 + t518;
t503 = t380 * t574 + t381 * t572 + t559;
t500 = t519 * t466;
t498 = -rSges(3,1) * t353 + rSges(3,2) * t352;
t497 = -rSges(6,1) * t258 + rSges(6,2) * t257;
t496 = -rSges(6,1) * t372 - rSges(6,2) * t476;
t493 = t521 - t565;
t492 = t518 - t560;
t491 = t504 * t466;
t490 = t356 + t547;
t136 = Icges(6,5) * t258 - Icges(6,6) * t257 + Icges(6,3) * t353;
t138 = Icges(6,4) * t258 - Icges(6,2) * t257 + Icges(6,6) * t353;
t140 = Icges(6,1) * t258 - Icges(6,4) * t257 + Icges(6,5) * t353;
t264 = Icges(6,5) * t372 + Icges(6,6) * t476 + Icges(6,3) * t431;
t266 = Icges(6,4) * t372 + Icges(6,2) * t476 + Icges(6,6) * t431;
t268 = Icges(6,1) * t372 + Icges(6,4) * t476 + Icges(6,5) * t431;
t51 = t138 * t428 + t140 * t429 - t266 * t366 + t268 * t367 + (t136 * t462 + t264 * t542) * t458;
t65 = -t257 * t303 + t258 * t304 + t281 * t431 + t282 * t476 + t283 * t372 + t302 * t353;
t489 = t65 / 0.2e1 + t51 / 0.2e1 + t537;
t135 = Icges(6,5) * t256 - Icges(6,6) * t255 - Icges(6,3) * t351;
t137 = Icges(6,4) * t256 - Icges(6,2) * t255 - Icges(6,6) * t351;
t139 = Icges(6,1) * t256 - Icges(6,4) * t255 - Icges(6,5) * t351;
t265 = Icges(6,5) * t374 + Icges(6,6) * t477 + Icges(6,3) * t433;
t267 = Icges(6,4) * t374 + Icges(6,2) * t477 + Icges(6,6) * t433;
t269 = Icges(6,1) * t374 + Icges(6,4) * t477 + Icges(6,5) * t433;
t52 = t137 * t428 + t139 * t429 - t267 * t366 + t269 * t367 + (t135 * t462 + t265 * t542) * t458;
t64 = -t255 * t303 + t256 * t304 + t281 * t433 + t282 * t477 + t283 * t374 - t302 * t351;
t488 = t64 / 0.2e1 + t52 / 0.2e1 + t536;
t487 = rSges(4,2) * t572 - rSges(4,3) * t430;
t486 = -t534 - t590;
t485 = t291 * t572 + t292 * t574 + t380 * t514 + t528;
t484 = t357 * t574 + t358 * t572 + t503;
t119 = t265 * t575 + t267 * t428 + t269 * t429;
t134 = t302 * t433 + t303 * t477 + t304 * t374;
t483 = -t134 / 0.2e1 - t119 / 0.2e1 - t530;
t118 = t264 * t575 + t266 * t428 + t268 * t429;
t133 = t302 * t431 + t303 * t476 + t304 * t372;
t482 = t133 / 0.2e1 + t118 / 0.2e1 + t529;
t481 = t492 * t466;
t480 = t426 + t490;
t479 = t499 + t557;
t478 = t466 * t510 - t590;
t448 = pkin(8) * t514;
t475 = t219 + t448;
t327 = rSges(3,1) * t431 - rSges(3,2) * t430 - rSges(3,3) * t572;
t472 = t273 * t572 + t274 * t574 + t357 * t514 + t485;
t471 = -t357 + t486 + t633;
t470 = -t458 * t541 - t347 + t475;
t469 = t358 + t480 - t535;
t468 = -qJD(1) * t547 - t274 + t479;
t467 = qJD(1) * t486 + t273 + t470;
t413 = (rSges(3,1) * t465 - rSges(3,2) * t462) * t544;
t400 = rSges(3,3) * t459 + (rSges(3,1) * t462 + rSges(3,2) * t465) * t458;
t393 = Icges(4,2) * t459 + (Icges(4,4) * t462 - Icges(4,6) * t465) * t458;
t391 = Icges(3,3) * t459 + (Icges(3,5) * t462 + Icges(3,6) * t465) * t458;
t389 = -Icges(5,3) * t459 + (-Icges(5,5) * t465 - Icges(5,6) * t462) * t458;
t328 = -rSges(5,3) * t574 + t549;
t326 = rSges(4,1) * t431 - t487;
t325 = -t431 * rSges(5,2) + rSges(5,3) * t572 + t587;
t294 = t330 + t547;
t293 = -t327 + t511;
t280 = -t327 * t459 - t400 * t572;
t279 = t330 * t459 - t400 * t574;
t270 = rSges(6,3) * t431 - t496;
t248 = rSges(3,3) * t515 - t498;
t247 = rSges(4,1) * t353 + rSges(4,2) * t515 + t588;
t246 = -rSges(5,2) * t353 - rSges(5,3) * t515 + t589;
t243 = -rSges(5,3) * t514 + t556;
t222 = t490 + t329;
t221 = t431 * t595 - t416 + t487 + t511;
t217 = (-t457 + (-rSges(3,3) - pkin(8)) * t574) * qJD(1) + t498;
t216 = t245 + t448 - t538;
t213 = t391 * t574 - t394 * t432 + t397 * t433;
t212 = t390 * t432 + t393 * t574 + t396 * t433;
t211 = -t389 * t574 - t392 * t433 + t395 * t432;
t210 = -t391 * t572 - t394 * t430 + t397 * t431;
t209 = t390 * t430 - t393 * t572 + t396 * t431;
t208 = t389 * t572 - t392 * t431 + t395 * t430;
t206 = t463 * t510 + t480 + t549;
t205 = t431 * t540 + t478 - t587 + t633;
t200 = (-t326 - t355) * t459 + t458 * t509;
t199 = t329 * t459 + t550 * t574 + t336;
t190 = t459 * t245 + (-t400 * t545 - t413 * t463) * t458;
t189 = -t459 * t248 + (t400 * t546 - t413 * t466) * t458;
t187 = t271 * t575 - t305 * t433;
t186 = -t270 * t575 + t305 * t431;
t185 = -t307 * t459 + (-t313 * t462 - t319 * t465) * t458;
t184 = -t306 * t459 + (-t312 * t462 - t318 * t465) * t458;
t183 = t311 * t459 + (t317 * t465 + t323 * t462) * t458;
t182 = t310 * t459 + (t316 * t465 + t322 * t462) * t458;
t181 = t315 * t459 + (-t309 * t465 + t321 * t462) * t458;
t180 = t314 * t459 + (-t308 * t465 + t320 * t462) * t458;
t155 = (t326 * t463 + t329 * t466) * t458 + t559;
t154 = (-t325 + t555) * t459 + t458 * t500;
t153 = t328 * t459 + t519 * t574 + t558;
t152 = t469 + t271;
t151 = t431 * t539 + t471 + t496;
t146 = t270 * t433 - t271 * t431;
t145 = -t588 + t595 * t353 + (-t457 + (-rSges(4,2) - pkin(8)) * t574) * qJD(1) + t557;
t144 = t244 + t475 - t538;
t142 = rSges(6,3) * t353 - t497;
t130 = -t352 * t394 + t353 * t397 - t430 * t407 + t431 * t410 + (t391 * t546 - t404 * t466) * t458;
t129 = t352 * t390 + t353 * t396 + t430 * t403 + t431 * t409 + (t393 * t546 - t406 * t466) * t458;
t128 = t352 * t395 - t353 * t392 - t431 * t405 + t430 * t408 + (-t389 * t546 + t402 * t466) * t458;
t127 = t350 * t394 - t351 * t397 - t432 * t407 + t433 * t410 + (t391 * t545 + t404 * t463) * t458;
t126 = -t350 * t390 - t351 * t396 + t432 * t403 + t433 * t409 + (t393 * t545 + t406 * t463) * t458;
t125 = -t350 * t395 + t351 * t392 - t433 * t405 + t432 * t408 + (-t389 * t545 - t402 * t463) * t458;
t124 = (t325 * t463 + t328 * t466) * t458 + t503;
t123 = -t589 + t540 * t353 + (-t457 + (rSges(5,3) - pkin(8)) * t574) * qJD(1) + t479;
t122 = qJD(1) * t478 + t470 + t556;
t121 = -t198 * t428 + t262 * t477;
t120 = t197 * t428 - t262 * t476;
t117 = t459 * t244 + t218 + (qJD(1) * t509 + t463 * t553) * t458;
t116 = t385 + (-t220 - t247) * t459 + (t399 * t546 + t466 * t553) * t458;
t115 = (-t270 + t524) * t459 + t458 * t491;
t114 = t271 * t459 + t504 * t574 + t525;
t113 = t469 + t563;
t112 = t431 * t604 - t476 * t594 + t471 + t494 - t591;
t111 = t265 * t433 + t267 * t477 + t269 * t374;
t110 = t264 * t433 + t266 * t477 + t268 * t374;
t109 = t265 * t431 + t267 * t476 + t269 * t372;
t108 = t264 * t431 + t266 * t476 + t268 * t372;
t105 = -t197 * t477 + t198 * t476;
t103 = -t433 * t560 + t563 * t575;
t102 = t431 * t560 - t564 * t575;
t99 = t229 * t459 + (t235 * t465 + t241 * t462 + (-t317 * t462 + t323 * t465) * qJD(2)) * t458;
t98 = t230 * t459 + (t236 * t465 + t242 * t462 + (-t316 * t462 + t322 * t465) * qJD(2)) * t458;
t97 = t233 * t459 + (-t227 * t465 + t239 * t462 + (t309 * t462 + t321 * t465) * qJD(2)) * t458;
t96 = t234 * t459 + (-t228 * t465 + t240 * t462 + (t308 * t462 + t320 * t465) * qJD(2)) * t458;
t95 = -t225 * t459 + (-t231 * t462 - t237 * t465 + (-t313 * t465 + t319 * t462) * qJD(2)) * t458;
t94 = -t226 * t459 + (-t232 * t462 - t238 * t465 + (-t312 * t465 + t318 * t462) * qJD(2)) * t458;
t93 = t459 * t243 + (qJD(1) * t500 + t463 * t522) * t458 + t562;
t92 = (-t246 + t561) * t459 + (t401 * t546 + t466 * t522) * t458 + t551;
t91 = (t270 * t463 + t271 * t466) * t458 + t484;
t90 = t353 * t539 + t468 + t497;
t89 = t467 + t141;
t79 = t80 * t459;
t77 = -t431 * t563 + t433 * t564;
t76 = t284 * t431 + t305 * t353 + (-t142 * t462 - t270 * t542) * t458;
t75 = -t284 * t433 + t305 * t351 + (t141 * t462 + t271 * t542) * t458;
t72 = (t524 - t564) * t459 + t458 * t481;
t71 = t459 * t563 + t492 * t574 + t525;
t70 = -t192 * t477 + t194 * t300 + t196 * t301;
t69 = -t191 * t477 + t193 * t300 + t195 * t301;
t68 = -t192 * t476 + t194 * t298 + t196 * t299;
t67 = -t191 * t476 + t193 * t298 + t195 * t299;
t66 = (t244 * t466 + t247 * t463 + (t326 * t466 + (-t329 - t356) * t463) * qJD(1)) * t458 + t528;
t63 = (t463 * t564 + t466 * t563) * t458 + t484;
t62 = -t141 * t431 + t142 * t433 - t270 * t351 - t271 * t353;
t61 = t459 * t141 + (qJD(1) * t491 + t463 * t505) * t458 + t527;
t60 = (-t142 + t526) * t459 + (t305 * t546 + t466 * t505) * t458 + t520;
t59 = (t243 * t466 + t246 * t463 + (t325 * t466 + (-t328 + t554) * t463) * qJD(1)) * t458 + t485;
t58 = t257 * t594 + t353 * t604 + t468 + t495 - t592;
t57 = t467 + t583;
t56 = t134 * t459 + (-t110 * t466 + t111 * t463) * t458;
t55 = t133 * t459 + (-t108 * t466 + t109 * t463) * t458;
t54 = t110 * t431 + t111 * t433 + t134 * t575;
t53 = t108 * t431 + t109 * t433 + t133 * t575;
t50 = -t150 * t476 - t197 * t366 + t257 * t262 + t428 * t88;
t49 = t150 * t477 + t198 * t366 - t255 * t262 - t428 * t87;
t47 = t48 * t459;
t44 = t135 * t431 + t137 * t476 + t139 * t372 - t257 * t267 + t258 * t269 + t265 * t353;
t43 = t136 * t431 + t138 * t476 + t140 * t372 - t257 * t266 + t258 * t268 + t264 * t353;
t42 = t135 * t433 + t137 * t477 + t139 * t374 - t255 * t267 + t256 * t269 - t265 * t351;
t41 = t136 * t433 + t138 * t477 + t140 * t374 - t255 * t266 + t256 * t268 - t264 * t351;
t40 = t565 * t431 + t560 * t353 + (-t462 * t582 - t542 * t564) * t458;
t39 = -t565 * t433 + t560 * t351 + (t462 * t583 + t542 * t563) * t458;
t38 = t107 * t459 + (t463 * t74 - t466 * t73) * t458;
t37 = (t141 * t466 + t142 * t463 + (t270 * t466 + (-t271 + t523) * t463) * qJD(1)) * t458 + t472;
t36 = t107 * t575 + t431 * t73 + t433 * t74;
t33 = -t107 * t428 - t476 * t73 - t477 * t74;
t32 = t197 * t255 - t198 * t257 + t476 * t87 - t477 * t88;
t31 = t583 * t459 + (qJD(1) * t481 + t463 * t493) * t458 + t527;
t30 = (t526 - t582) * t459 + (t466 * t493 + t546 * t560) * t458 + t520;
t29 = t101 * t459 + (t463 * t70 - t466 * t69) * t458;
t28 = t100 * t459 + (t463 * t68 - t466 * t67) * t458;
t27 = t101 * t575 + t431 * t69 + t433 * t70;
t26 = t100 * t575 + t431 * t67 + t433 * t68;
t25 = -t101 * t428 - t476 * t69 - t477 * t70;
t24 = -t100 * t428 - t476 * t67 - t477 * t68;
t23 = -t351 * t564 - t353 * t563 - t431 * t583 + t433 * t582;
t20 = t161 * t194 + t162 * t196 + t192 * t257 + t298 * t83 + t299 * t85 - t476 * t81;
t19 = t161 * t193 + t162 * t195 + t191 * t257 + t298 * t84 + t299 * t86 - t476 * t82;
t18 = t159 * t194 + t160 * t196 + t192 * t255 + t300 * t83 + t301 * t85 - t477 * t81;
t17 = t159 * t193 + t160 * t195 + t191 * t255 + t300 * t84 + t301 * t86 - t477 * t82;
t16 = (t583 * t466 + t582 * t463 + (t564 * t466 + (t523 - t563) * t463) * qJD(1)) * t458 + t472;
t15 = t79 + (t52 * t463 - t51 * t466 + (t118 * t463 + t119 * t466) * qJD(1)) * t458;
t14 = t118 * t353 - t119 * t351 + t431 * t51 + t433 * t52 + t584;
t13 = t65 * t459 + (-t43 * t466 + t44 * t463 + (t108 * t463 + t109 * t466) * qJD(1)) * t458;
t12 = t64 * t459 + (-t41 * t466 + t42 * t463 + (t110 * t463 + t111 * t466) * qJD(1)) * t458;
t11 = t108 * t353 - t109 * t351 + t43 * t431 + t433 * t44 + (t133 * t542 + t462 * t65) * t458;
t10 = t110 * t353 - t111 * t351 + t41 * t431 + t42 * t433 + (t134 * t542 + t462 * t64) * t458;
t9 = t47 + (-t21 * t466 + t22 * t463 + (t463 * t73 + t466 * t74) * qJD(1)) * t458;
t8 = t21 * t431 + t22 * t433 - t351 * t74 + t353 * t73 + t585;
t7 = -t21 * t476 - t22 * t477 + t255 * t74 + t257 * t73 + t586;
t6 = t35 * t459 + (-t19 * t466 + t20 * t463 + (t463 * t67 + t466 * t68) * qJD(1)) * t458;
t5 = t34 * t459 + (-t17 * t466 + t18 * t463 + (t463 * t69 + t466 * t70) * qJD(1)) * t458;
t4 = t19 * t431 + t20 * t433 - t351 * t68 + t353 * t67 + (t100 * t542 + t35 * t462) * t458;
t3 = t17 * t431 + t18 * t433 - t351 * t70 + t353 * t69 + (t101 * t542 + t34 * t462) * t458;
t2 = t100 * t366 - t19 * t476 - t20 * t477 + t255 * t68 + t257 * t67 - t35 * t428;
t1 = t101 * t366 - t17 * t476 - t18 * t477 + t255 * t70 + t257 * t69 - t34 * t428;
t45 = [t48 + t80 - t408 * t573 - t458 * t570 - t392 * t512 + (t216 * t294 + t217 * t293) * t616 + (t144 * t222 + t145 * t221) * t615 + (t122 * t206 + t123 * t205) * t614 + (t151 * t90 + t152 * t89) * t613 + (t112 * t58 + t113 * t57) * t612 + t634; t47 + t79 + m(3) * (t189 * t293 + t190 * t294 + t216 * t279 + t217 * t280) + m(4) * (t116 * t221 + t117 * t222 + t144 * t199 + t145 * t200) + m(5) * (t122 * t153 + t123 * t154 + t205 * t92 + t206 * t93) + m(6) * (t114 * t89 + t115 * t90 + t151 * t60 + t152 * t61) + m(7) * (t112 * t30 + t113 * t31 + t57 * t71 + t58 * t72) + ((-t96 / 0.2e1 - t98 / 0.2e1 - t94 / 0.2e1 - t129 / 0.2e1 - t130 / 0.2e1 - t128 / 0.2e1 - t489) * t466 + (t97 / 0.2e1 + t99 / 0.2e1 + t95 / 0.2e1 + t125 / 0.2e1 + t127 / 0.2e1 + t126 / 0.2e1 + t488) * t463 + ((t181 / 0.2e1 + t183 / 0.2e1 + t185 / 0.2e1 + t211 / 0.2e1 + t212 / 0.2e1 + t213 / 0.2e1 - t483) * t466 + (t180 / 0.2e1 + t182 / 0.2e1 + t184 / 0.2e1 + t208 / 0.2e1 + t209 / 0.2e1 + t210 / 0.2e1 + t482) * t463) * qJD(1)) * t458 + t628; (t16 * t63 + t30 * t72 + t31 * t71) * t612 + (t114 * t61 + t115 * t60 + t37 * t91) * t613 + (t124 * t59 + t153 * t93 + t154 * t92) * t614 + (t116 * t200 + t117 * t199 + t155 * t66) * t615 + (t280 * t189 + t279 * t190 + (t327 * t463 + t330 * t466) * (t245 * t466 + t248 * t463 + (t327 * t466 - t330 * t463) * qJD(1)) * t611) * t616 + (t12 + t5) * t574 + (-t13 - t6) * t572 + (t55 + t28) * t515 + (t56 + t29) * t514 + ((t463 * t631 + t466 * t632) * t515 + (t463 * t630 - t466 * t629) * t514 + (t629 * t546 + t630 * t545 + (t350 * t621 + t351 * t619 - t432 * t626 + t433 * t624 + t514 * t623) * t466 + ((-t545 * t622 + t617) * t458 + t625 * t433 - t627 * t432 + t618 * t351 + t620 * t350) * t463) * t574 + (t632 * t546 - t631 * t545 + (t352 * t620 + t353 * t618 + t430 * t627 - t431 * t625 + t515 * t622) * t463 + ((-t546 * t623 + t617) * t458 - t624 * t431 + t626 * t430 + t619 * t353 + t621 * t352) * t466) * t572) * t458 + (t15 + t9 + (t127 + t126 + t125) * t574 + (-t130 - t129 - t128) * t572 + (t209 + t208 + t210) * t515 + (t213 + t212 + t211) * t514 + ((-t94 - t96 - t98) * t466 + (t95 + t97 + t99) * t463 + ((t181 + t183 + t185) * t466 + (t180 + t182 + t184) * t463) * qJD(1)) * t458 + t628) * t459; m(7) * (-t112 * t350 + t113 * t352 + t430 * t57 + t432 * t58) + m(6) * (-t151 * t350 + t152 * t352 + t430 * t89 + t432 * t90) + m(4) * (t144 * t430 + t145 * t432 - t221 * t350 + t222 * t352) + m(5) * (t122 * t430 + t123 * t432 - t205 * t350 + t206 * t352); m(7) * (t30 * t432 + t31 * t430 - t350 * t72 + t352 * t71 + (-t16 * t465 + t543 * t63) * t458) + m(6) * (t114 * t352 - t115 * t350 + t430 * t61 + t432 * t60 + (-t37 * t465 + t543 * t91) * t458) + m(5) * (t153 * t352 - t154 * t350 + t430 * t93 + t432 * t92 + (t124 * t543 - t465 * t59) * t458) + m(4) * (t116 * t432 + t117 * t430 + t199 * t352 - t200 * t350 + (t155 * t543 - t465 * t66) * t458); 0.4e1 * (m(4) / 0.2e1 + t507) * (-t462 * t542 * t611 - t350 * t432 + t352 * t430); ((-t112 * t545 - t113 * t546 - t463 * t58 + t466 * t57) * t605 + (-t151 * t545 - t152 * t546 - t463 * t90 + t466 * t89) * t607 + (t122 * t466 - t123 * t463 - t205 * t545 - t206 * t546) * t609) * t610; 0.2e1 * (t16 * t606 + t37 * t608 - m(5) * t59 / 0.2e1) * t459 + 0.2e1 * ((-t30 * t463 + t31 * t466 - t545 * t72 - t546 * t71) * t605 + (-t114 * t546 - t115 * t545 - t463 * t60 + t466 * t61) * t607 + (-t153 * t546 - t154 * t545 - t463 * t92 + t466 * t93) * t609) * t458; t507 * (-t459 * t543 + t350 * t463 + t352 * t466 + (-t430 * t463 - t432 * t466) * qJD(1)) * t610; 0; m(7) * (t102 * t58 + t103 * t57 + t112 * t40 + t113 * t39) + m(6) * (t151 * t76 + t152 * t75 + t186 * t90 + t187 * t89) + t488 * t433 + t489 * t431 + t482 * t353 + t483 * t351 + t584 + t585; (t8 / 0.2e1 + t14 / 0.2e1) * t459 + (t5 / 0.2e1 + t12 / 0.2e1) * t433 + (t6 / 0.2e1 + t13 / 0.2e1) * t431 + (t28 / 0.2e1 + t55 / 0.2e1) * t353 + (-t29 / 0.2e1 - t56 / 0.2e1) * t351 + m(7) * (t102 * t30 + t103 * t31 + t16 * t77 + t23 * t63 + t39 * t71 + t40 * t72) + m(6) * (t114 * t75 + t115 * t76 + t146 * t37 + t186 * t60 + t187 * t61 + t62 * t91) + ((-t4 / 0.2e1 - t11 / 0.2e1) * t466 + (t3 / 0.2e1 + t10 / 0.2e1) * t463 + (t9 / 0.2e1 + t15 / 0.2e1) * t462 + (t38 / 0.2e1 + t163 * t597 + (-t118 * t466 + t119 * t463) * t458 / 0.2e1) * t542 + ((t27 / 0.2e1 + t54 / 0.2e1) * t466 + (t26 / 0.2e1 + t53 / 0.2e1) * t463) * qJD(1)) * t458; m(6) * (-t186 * t350 + t187 * t352 + t430 * t75 + t432 * t76 + (t146 * t543 - t465 * t62) * t458) + m(7) * (-t102 * t350 + t103 * t352 + t39 * t430 + t40 * t432 + (-t23 * t465 + t543 * t77) * t458); 0.2e1 * (t23 * t606 + t608 * t62) * t459 + 0.2e1 * ((-t186 * t545 - t187 * t546 - t463 * t76 + t466 * t75) * t607 + (-t102 * t545 - t103 * t546 + t39 * t466 - t40 * t463) * t605) * t458; (t14 + t8) * t575 + (t10 + t3) * t433 + (t11 + t4) * t431 + (t26 + t53) * t353 + (-t27 - t54) * t351 + (t118 * t431 + t119 * t433 + t163 * t575 + t36) * t512 + (t102 * t40 + t103 * t39 + t23 * t77) * t612 + (t146 * t62 + t186 * t76 + t187 * t75) * t613; m(7) * (t112 * t50 + t113 * t49 + t120 * t58 + t121 * t57) - t536 * t477 - t537 * t476 + t529 * t257 + t530 * t255 + t586; m(7) * (t105 * t16 + t120 * t30 + t121 * t31 + t32 * t63 + t49 * t71 + t50 * t72) + t28 * t602 + t6 * t600 + t29 * t603 + t5 * t599 + t7 * t597 + t38 * t601 + t9 * t598 + (-t466 * t2 / 0.2e1 + t1 * t596 + (t24 * t596 + t466 * t25 / 0.2e1) * qJD(1)) * t458; m(7) * (-t120 * t350 + t121 * t352 + t430 * t49 + t432 * t50 + (t105 * t543 - t32 * t465) * t458); m(7) * (-t32 * t459 + (-t463 * t50 + t466 * t49 + (-t120 * t466 - t121 * t463) * qJD(1)) * t458); m(7) * (t102 * t50 + t103 * t49 + t105 * t23 + t120 * t40 + t121 * t39 + t32 * t77) + t36 * t601 + t8 * t598 + t353 * t24 / 0.2e1 + t431 * t2 / 0.2e1 + t27 * t603 + t3 * t599 + t26 * t602 + t4 * t600 - t351 * t25 / 0.2e1 + t433 * t1 / 0.2e1 + (t33 * t542 / 0.2e1 + t462 * t7 / 0.2e1) * t458; t255 * t25 - t477 * t1 + t257 * t24 - t476 * t2 + t366 * t33 - t428 * t7 + (t105 * t32 + t120 * t50 + t121 * t49) * t612;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t45(1) t45(2) t45(4) t45(7) t45(11) t45(16); t45(2) t45(3) t45(5) t45(8) t45(12) t45(17); t45(4) t45(5) t45(6) t45(9) t45(13) t45(18); t45(7) t45(8) t45(9) t45(10) t45(14) t45(19); t45(11) t45(12) t45(13) t45(14) t45(15) t45(20); t45(16) t45(17) t45(18) t45(19) t45(20) t45(21);];
Mq  = res;
