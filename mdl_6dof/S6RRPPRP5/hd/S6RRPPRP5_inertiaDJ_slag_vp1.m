% Calculate time derivative of joint inertia matrix for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:32
% EndTime: 2019-03-09 08:42:12
% DurationCPUTime: 27.27s
% Computational Cost: add. (23270->1050), mult. (37700->1457), div. (0->0), fcn. (35400->8), ass. (0->503)
t358 = pkin(9) + qJ(5);
t349 = sin(t358);
t366 = sin(qJ(2));
t368 = cos(qJ(2));
t350 = cos(t358);
t581 = Icges(7,5) * t350;
t454 = Icges(7,1) * t349 - t581;
t238 = Icges(7,4) * t366 - t368 * t454;
t585 = Icges(6,4) * t350;
t455 = Icges(6,1) * t349 + t585;
t239 = Icges(6,5) * t366 - t368 * t455;
t656 = -t239 - t238;
t676 = t656 * t349;
t579 = Icges(4,6) * t366;
t589 = Icges(3,4) * t366;
t675 = t579 + t589 + (Icges(3,2) + Icges(4,3)) * t368;
t578 = Icges(4,6) * t368;
t588 = Icges(3,4) * t368;
t674 = t578 + t588 + (Icges(3,1) + Icges(4,2)) * t366;
t445 = Icges(6,5) * t349 + Icges(6,6) * t350;
t235 = Icges(6,3) * t366 - t368 * t445;
t448 = Icges(7,4) * t349 - Icges(7,6) * t350;
t236 = Icges(7,2) * t366 - t368 * t448;
t673 = t235 + t236;
t367 = sin(qJ(1));
t533 = qJD(5) * t366;
t490 = qJD(1) + t533;
t422 = t367 * t490;
t369 = cos(qJ(1));
t570 = t350 * t369;
t537 = qJD(2) * t368;
t503 = t367 * t537;
t540 = qJD(1) * t369;
t638 = t366 * t540 + t503;
t145 = -qJD(5) * t570 + t349 * t422 - t350 * t638;
t489 = qJD(1) * t366 + qJD(5);
t146 = t350 * t422 + (t369 * t489 + t503) * t349;
t567 = t367 * t350;
t256 = t349 * t369 + t366 * t567;
t591 = rSges(7,3) + qJ(6);
t661 = rSges(7,1) + pkin(5);
t672 = t256 * qJD(6) - t591 * t145 - t146 * t661;
t573 = t349 * t367;
t257 = t366 * t573 - t570;
t671 = -t256 * t591 + t257 * t661;
t669 = t675 * qJD(2);
t668 = t674 * qJD(2);
t666 = -Icges(5,3) / 0.2e1;
t364 = cos(pkin(9));
t363 = sin(pkin(9));
t566 = t367 * t363;
t568 = t366 * t369;
t283 = t364 * t568 - t566;
t220 = qJD(1) * t283 + t364 * t503;
t665 = -t220 / 0.2e1;
t518 = t363 * t568;
t565 = t367 * t364;
t284 = t518 + t565;
t221 = qJD(1) * t284 + t363 * t503;
t664 = -t221 / 0.2e1;
t663 = t367 / 0.2e1;
t662 = -t369 / 0.2e1;
t660 = -qJD(1) / 0.2e1;
t659 = qJD(1) / 0.2e1;
t582 = Icges(7,5) * t349;
t444 = -Icges(7,3) * t350 + t582;
t234 = Icges(7,6) * t366 - t368 * t444;
t586 = Icges(6,4) * t349;
t449 = Icges(6,2) * t350 + t586;
t237 = Icges(6,6) * t366 - t368 * t449;
t658 = ((t234 - t237) * t350 + t676) * t368 + t673 * t366;
t532 = qJD(5) * t368;
t174 = (-Icges(7,1) * t350 - t582) * t532 + (Icges(7,4) * t368 + t366 * t454) * qJD(2);
t175 = (-Icges(6,1) * t350 + t586) * t532 + (Icges(6,5) * t368 + t366 * t455) * qJD(2);
t657 = -t175 - t174;
t170 = (-Icges(7,3) * t349 - t581) * t532 + (Icges(7,6) * t368 + t366 * t444) * qJD(2);
t171 = (-Icges(6,5) * t350 + Icges(6,6) * t349) * t532 + (Icges(6,3) * t368 + t366 * t445) * qJD(2);
t172 = (-Icges(7,4) * t350 - Icges(7,6) * t349) * t532 + (Icges(7,2) * t368 + t366 * t448) * qJD(2);
t501 = t349 * t532;
t539 = qJD(2) * t366;
t506 = t350 * t539;
t571 = t350 * t368;
t655 = t170 * t571 - t234 * t506 + (t501 + t506) * t237 + t673 * t537 - t539 * t676 + (t172 + t171) * t366;
t564 = t367 * t368;
t557 = rSges(7,2) * t564 + t671;
t495 = t557 * t369;
t254 = -t350 * t568 + t573;
t255 = t349 * t568 + t567;
t563 = t368 * t369;
t558 = rSges(7,2) * t563 + t254 * t591 + t255 * t661;
t654 = -t367 * t558 + t495;
t536 = qJD(2) * t369;
t504 = t366 * t536;
t541 = qJD(1) * t367;
t399 = t368 * t541 + t504;
t159 = Icges(7,4) * t257 + Icges(7,2) * t564 - Icges(7,6) * t256;
t155 = Icges(7,5) * t257 + Icges(7,6) * t564 - Icges(7,3) * t256;
t163 = Icges(7,1) * t257 + Icges(7,4) * t564 - Icges(7,5) * t256;
t436 = t155 * t350 - t163 * t349;
t65 = t159 * t366 + t368 * t436;
t157 = Icges(6,5) * t257 + Icges(6,6) * t256 + Icges(6,3) * t564;
t161 = Icges(6,4) * t257 + Icges(6,2) * t256 + Icges(6,6) * t564;
t165 = Icges(6,1) * t257 + Icges(6,4) * t256 + Icges(6,5) * t564;
t434 = t161 * t350 + t165 * t349;
t67 = t157 * t366 - t368 * t434;
t599 = t65 + t67;
t158 = Icges(7,4) * t255 + Icges(7,2) * t563 + Icges(7,6) * t254;
t154 = Icges(7,5) * t255 + Icges(7,6) * t563 + Icges(7,3) * t254;
t162 = Icges(7,1) * t255 + Icges(7,4) * t563 + Icges(7,5) * t254;
t437 = t154 * t350 - t162 * t349;
t64 = t158 * t366 + t368 * t437;
t156 = Icges(6,5) * t255 - Icges(6,6) * t254 + Icges(6,3) * t563;
t160 = Icges(6,4) * t255 - Icges(6,2) * t254 + Icges(6,6) * t563;
t164 = Icges(6,1) * t255 - Icges(6,4) * t254 + Icges(6,5) * t563;
t435 = t160 * t350 + t164 * t349;
t66 = t156 * t366 - t368 * t435;
t600 = t64 + t66;
t653 = t367 * t600 - t369 * t599;
t652 = t367 * t599 + t369 * t600;
t544 = t367 ^ 2 + t369 ^ 2;
t651 = -0.1e1 + t544;
t650 = qJD(2) / 0.2e1;
t502 = t368 * t536;
t147 = qJD(1) * t256 + qJD(5) * t255 - t350 * t502;
t148 = t490 * t570 + (-t367 * t489 + t502) * t349;
t82 = Icges(7,5) * t148 - Icges(7,6) * t399 + Icges(7,3) * t147;
t86 = Icges(7,4) * t148 - Icges(7,2) * t399 + Icges(7,6) * t147;
t90 = Icges(7,1) * t148 - Icges(7,4) * t399 + Icges(7,5) * t147;
t20 = (-qJD(2) * t437 + t86) * t366 + (qJD(2) * t158 - t349 * t90 + t350 * t82 + (-t154 * t349 - t162 * t350) * qJD(5)) * t368;
t84 = Icges(6,5) * t148 - Icges(6,6) * t147 - Icges(6,3) * t399;
t88 = Icges(6,4) * t148 - Icges(6,2) * t147 - Icges(6,6) * t399;
t92 = Icges(6,1) * t148 - Icges(6,4) * t147 - Icges(6,5) * t399;
t22 = (qJD(2) * t435 + t84) * t366 + (qJD(2) * t156 - t349 * t92 - t350 * t88 + (t160 * t349 - t164 * t350) * qJD(5)) * t368;
t649 = t20 + t22;
t538 = qJD(2) * t367;
t505 = t366 * t538;
t508 = t368 * t540;
t400 = -t505 + t508;
t81 = Icges(7,5) * t146 + Icges(7,6) * t400 + Icges(7,3) * t145;
t85 = Icges(7,4) * t146 + Icges(7,2) * t400 + Icges(7,6) * t145;
t89 = Icges(7,1) * t146 + Icges(7,4) * t400 + Icges(7,5) * t145;
t21 = (-qJD(2) * t436 + t85) * t366 + (qJD(2) * t159 - t349 * t89 + t350 * t81 + (-t155 * t349 - t163 * t350) * qJD(5)) * t368;
t83 = Icges(6,5) * t146 - Icges(6,6) * t145 + Icges(6,3) * t400;
t87 = Icges(6,4) * t146 - Icges(6,2) * t145 + Icges(6,6) * t400;
t91 = Icges(6,1) * t146 - Icges(6,4) * t145 + Icges(6,5) * t400;
t23 = (qJD(2) * t434 + t83) * t366 + (qJD(2) * t157 - t349 * t91 - t350 * t87 + (t161 * t349 - t165 * t350) * qJD(5)) * t368;
t648 = t21 + t23;
t56 = t156 * t563 - t254 * t160 + t255 * t164;
t57 = t157 * t563 - t254 * t161 + t255 * t165;
t463 = t367 * t57 + t369 * t56;
t54 = t254 * t154 + t158 * t563 + t255 * t162;
t55 = t254 * t155 + t159 * t563 + t255 * t163;
t464 = t367 * t55 + t369 * t54;
t97 = t254 * t234 + t236 * t563 + t255 * t238;
t98 = t235 * t563 - t254 * t237 + t255 * t239;
t647 = (t463 + t464) * t368 + (t97 + t98) * t366;
t100 = t235 * t564 + t237 * t256 + t239 * t257;
t60 = t156 * t564 + t160 * t256 + t164 * t257;
t61 = t157 * t564 + t161 * t256 + t165 * t257;
t461 = t367 * t61 + t369 * t60;
t58 = -t154 * t256 + t158 * t564 + t162 * t257;
t59 = -t155 * t256 + t159 * t564 + t163 * t257;
t462 = t367 * t59 + t369 * t58;
t99 = -t234 * t256 + t236 * t564 + t238 * t257;
t601 = (t461 + t462) * t368 + (t100 + t99) * t366;
t646 = t254 * qJD(6) + t591 * t147 + t148 * t661;
t441 = -Icges(4,3) * t366 + t578;
t264 = Icges(4,5) * t367 - t369 * t441;
t443 = Icges(4,2) * t368 - t579;
t266 = Icges(4,4) * t367 - t369 * t443;
t424 = t264 * t366 - t266 * t368;
t645 = t367 * t424;
t453 = -Icges(3,2) * t366 + t588;
t261 = Icges(3,6) * t367 + t369 * t453;
t458 = Icges(3,1) * t368 - t589;
t263 = Icges(3,5) * t367 + t369 * t458;
t425 = t261 * t366 - t263 * t368;
t644 = t367 * t425;
t634 = Icges(4,4) * t369 + t367 * t443;
t635 = Icges(4,5) * t369 + t367 * t441;
t423 = -t366 * t635 + t368 * t634;
t642 = t369 * t423;
t260 = -Icges(3,6) * t369 + t367 * t453;
t262 = -Icges(3,5) * t369 + t367 * t458;
t426 = t260 * t366 - t262 * t368;
t641 = t369 * t426;
t640 = t367 * rSges(4,1) - rSges(4,2) * t563;
t303 = t367 * pkin(3) + qJ(4) * t563;
t361 = t368 ^ 2;
t639 = -t366 ^ 2 + t361;
t637 = -rSges(3,2) * t568 + t367 * rSges(3,3);
t348 = pkin(3) * t540;
t636 = -qJ(4) * t399 + t348;
t447 = Icges(3,5) * t368 - Icges(3,6) * t366;
t258 = -Icges(3,3) * t369 + t367 * t447;
t451 = Icges(4,4) * t368 - Icges(4,5) * t366;
t633 = Icges(4,1) * t369 + t367 * t451;
t632 = -t367 * t557 - t369 * t558;
t631 = 2 * m(3);
t630 = 2 * m(4);
t629 = 2 * m(5);
t628 = 2 * m(6);
t627 = 2 * m(7);
t626 = m(4) / 0.2e1;
t625 = -m(5) / 0.2e1;
t624 = m(5) / 0.2e1;
t623 = -m(6) / 0.2e1;
t622 = m(6) / 0.2e1;
t621 = -m(7) / 0.2e1;
t620 = m(7) / 0.2e1;
t450 = Icges(5,4) * t363 + Icges(5,2) * t364;
t251 = Icges(5,6) * t366 - t368 * t450;
t619 = t251 / 0.2e1;
t456 = Icges(5,1) * t363 + Icges(5,4) * t364;
t252 = Icges(5,5) * t366 - t368 * t456;
t618 = t252 / 0.2e1;
t617 = t283 / 0.2e1;
t616 = t284 / 0.2e1;
t615 = -t363 / 0.2e1;
t614 = -t364 / 0.2e1;
t609 = -rSges(7,2) - pkin(2);
t608 = -rSges(6,3) - pkin(2);
t321 = rSges(3,1) * t366 + rSges(3,2) * t368;
t607 = m(3) * t321;
t606 = pkin(2) * t368;
t605 = pkin(4) * t363;
t598 = rSges(7,2) * t400 - t672;
t597 = -rSges(7,2) * t399 + t646;
t596 = rSges(4,1) * t369;
t595 = rSges(4,2) * t366;
t594 = rSges(3,3) * t369;
t593 = rSges(5,3) * t366;
t592 = -rSges(4,3) - qJ(3);
t576 = qJ(3) * t366;
t575 = qJ(3) * t368;
t470 = -rSges(6,1) * t257 - rSges(6,2) * t256;
t169 = rSges(6,3) * t564 - t470;
t574 = t169 * t369;
t569 = t363 * t368;
t365 = -pkin(8) - qJ(4);
t562 = -qJ(4) - t365;
t560 = t148 * rSges(6,1) - t147 * rSges(6,2);
t468 = rSges(7,1) * t349 - rSges(7,3) * t350;
t559 = (pkin(5) * t539 - qJ(6) * t532) * t349 + (-qJ(6) * t539 + (-pkin(5) * qJD(5) + qJD(6)) * t368) * t350 + (-rSges(7,1) * t350 - rSges(7,3) * t349) * t532 + (rSges(7,2) * t368 + t366 * t468) * qJD(2);
t285 = t363 * t369 + t366 * t565;
t222 = -qJD(1) * t285 + t364 * t502;
t286 = -t364 * t369 + t366 * t566;
t483 = t363 * t502;
t223 = -qJD(1) * t286 + t483;
t556 = t223 * rSges(5,1) + t222 * rSges(5,2);
t555 = rSges(7,2) * t366 + (-pkin(5) * t349 + qJ(6) * t350 - t468) * t368;
t465 = t576 + t606;
t289 = t465 * t367;
t290 = pkin(2) * t563 + qJ(3) * t568;
t554 = t367 * t289 + t369 * t290;
t281 = qJD(2) * t465 - qJD(3) * t368;
t467 = -rSges(4,2) * t368 + rSges(4,3) * t366;
t553 = -t467 * qJD(2) - t281;
t552 = -t290 - t303;
t319 = pkin(2) * t366 - t575;
t291 = t319 * t541;
t512 = t366 * t541;
t551 = qJ(4) * t512 + t291;
t310 = t365 * t508;
t332 = pkin(2) * t505;
t550 = t310 + t332;
t466 = rSges(4,3) * t368 + t595;
t549 = -t319 + t466;
t344 = pkin(4) * t364 + pkin(3);
t548 = -t369 * t344 - t365 * t564;
t535 = qJD(3) * t366;
t547 = qJ(3) * t502 + t369 * t535;
t546 = rSges(3,2) * t512 + rSges(3,3) * t540;
t545 = t369 * pkin(1) + t367 * pkin(7);
t259 = Icges(3,3) * t367 + t369 * t447;
t543 = qJD(1) * t259;
t268 = Icges(4,1) * t367 - t369 * t451;
t542 = qJD(1) * t268;
t534 = qJD(4) * t368;
t530 = -rSges(5,3) - pkin(2) - qJ(4);
t528 = t366 * t605;
t527 = pkin(4) * t569;
t36 = t54 * t367 - t369 * t55;
t37 = t56 * t367 - t369 * t57;
t524 = -t36 / 0.2e1 - t37 / 0.2e1;
t38 = t58 * t367 - t369 * t59;
t39 = t60 * t367 - t369 * t61;
t523 = t38 / 0.2e1 + t39 / 0.2e1;
t185 = Icges(5,5) * t284 + Icges(5,6) * t283 + Icges(5,3) * t563;
t522 = t185 * t564;
t521 = t185 * t563;
t186 = Icges(5,5) * t286 + Icges(5,6) * t285 + Icges(5,3) * t564;
t520 = t186 * t564;
t519 = t186 * t563;
t517 = t367 * (pkin(2) * t508 + qJ(3) * t638 + t367 * t535 - t332) + t369 * (-pkin(2) * t399 - qJ(3) * t512 + t547) + t289 * t540;
t419 = pkin(4) * t518 + t367 * t344 - t365 * t563;
t209 = t419 - t303;
t516 = -t209 + t552;
t276 = t366 * t562 - t527;
t515 = t276 * t541 + t551;
t167 = t255 * rSges(6,1) - t254 * rSges(6,2) + rSges(6,3) * t563;
t193 = t284 * rSges(5,1) + t283 * rSges(5,2) + rSges(5,3) * t563;
t355 = t369 * pkin(7);
t514 = t355 - t548;
t347 = pkin(7) * t540;
t513 = t347 + t547;
t469 = rSges(6,1) * t349 + rSges(6,2) * t350;
t242 = rSges(6,3) * t366 - t368 * t469;
t510 = t242 * t541;
t499 = t361 * qJD(5) * t349;
t498 = t366 * t537;
t497 = -t451 * qJD(2) / 0.2e1 + t447 * t650;
t496 = -qJ(3) - t605;
t230 = t549 * t369;
t494 = -qJ(4) * t366 - t319;
t493 = qJD(1) * t555;
t173 = (Icges(6,2) * t349 - t585) * t532 + (Icges(6,6) * t368 + t366 * t449) * qJD(2);
t492 = t658 * t537 + (((qJD(5) * t656 - t173) * t350 + (-qJD(5) * t234 + t657) * t349) * t368 + t655) * t366;
t491 = t624 + t622 + t620;
t356 = t369 * pkin(3);
t304 = qJ(4) * t564 - t356;
t488 = t369 * t303 + t367 * t304 + t554;
t487 = pkin(4) * t483 + t344 * t540 + t365 * t399;
t334 = t369 * t534;
t486 = t334 + t513;
t485 = rSges(4,1) * t540 + rSges(4,2) * t399 + rSges(4,3) * t502;
t484 = t545 + t290;
t472 = rSges(5,1) * t363 + rSges(5,2) * t364;
t253 = -t368 * t472 + t593;
t482 = -t253 + t494;
t481 = -t276 + t494;
t480 = t367 * t493;
t479 = t366 * t592 - pkin(1);
t478 = t263 / 0.2e1 - t266 / 0.2e1 + t185 / 0.2e1;
t477 = -t634 / 0.2e1 - t186 / 0.2e1 - t262 / 0.2e1;
t476 = t496 * t368;
t475 = rSges(3,1) * t368 - rSges(3,2) * t366;
t474 = -t221 * rSges(5,1) - t220 * rSges(5,2);
t473 = -rSges(5,1) * t286 - rSges(5,2) * t285;
t471 = t146 * rSges(6,1) - t145 * rSges(6,2);
t410 = t366 * t496 - pkin(1);
t378 = t368 * t609 + t410;
t77 = t367 * t378 + t514 - t671;
t379 = t419 + t484;
t78 = t379 + t558;
t460 = t367 * t78 + t369 * t77;
t446 = Icges(5,5) * t363 + Icges(5,6) * t364;
t377 = t368 * t608 + t410;
t104 = t367 * t377 + t470 + t514;
t105 = t379 + t167;
t439 = t104 * t369 + t105 * t367;
t396 = t368 * t530 - pkin(1) - t576;
t380 = t396 * t367;
t125 = t355 + t356 + t380 + t473;
t126 = t484 + t193 + t303;
t438 = t125 * t369 + t126 * t367;
t433 = -t167 * t369 - t367 * t169;
t432 = t167 * t367 - t574;
t427 = t254 * t369 - t256 * t367;
t421 = -t534 - t535;
t273 = rSges(3,1) * t563 + t637;
t274 = rSges(4,3) * t568 + t640;
t420 = -pkin(1) - t475;
t418 = -t242 + t481;
t34 = t147 * t234 + t148 * t238 + t254 * t170 + t172 * t563 + t255 * t174 - t236 * t399;
t35 = -t147 * t237 + t148 * t239 + t171 * t563 - t254 * t173 + t255 * t175 - t235 * t399;
t417 = t22 / 0.2e1 + t20 / 0.2e1 + t34 / 0.2e1 + t35 / 0.2e1;
t32 = t145 * t234 + t146 * t238 - t256 * t170 + t172 * t564 + t257 * t174 + t236 * t400;
t33 = -t145 * t237 + t146 * t239 + t171 * t564 + t256 * t173 + t257 * t175 + t235 * t400;
t416 = t23 / 0.2e1 + t21 / 0.2e1 + t32 / 0.2e1 + t33 / 0.2e1;
t415 = -t66 / 0.2e1 - t64 / 0.2e1 - t97 / 0.2e1 - t98 / 0.2e1;
t414 = -t99 / 0.2e1 - t100 / 0.2e1 - t67 / 0.2e1 - t65 / 0.2e1;
t413 = -qJ(4) * t368 + t528;
t325 = qJ(4) * t505;
t412 = t367 * (qJD(1) * t303 + t367 * t534 - t325) + t369 * (t334 + t636) + t304 * t540 + t517;
t210 = t367 * t413 + t356 + t548;
t411 = t369 * t209 + t367 * t210 + t488;
t182 = t482 * t369;
t408 = qJD(2) * t321;
t407 = t481 - t555;
t404 = qJD(2) * (Icges(4,4) * t366 + Icges(4,5) * t368);
t403 = qJD(2) * (-Icges(3,5) * t366 - Icges(3,6) * t368);
t136 = t418 * t369;
t398 = -qJ(4) * t537 - qJD(4) * t366 - t281;
t395 = t486 + t487;
t10 = -t654 * t539 + (qJD(1) * t632 - t597 * t367 + t598 * t369) * t368;
t75 = -t366 * t557 + t555 * t564;
t76 = t366 * t558 - t555 * t563;
t394 = t536 * t75 + t538 * t76 - t10;
t113 = t407 * t369;
t112 = t407 * t367;
t385 = t367 * (-t310 + t325 + (t365 * t366 + t527) * t538 + (t413 * t369 + (-pkin(3) + t344) * t367) * qJD(1)) + t369 * (-t512 * t605 + t487 - t636) + t210 * t540 + t412;
t9 = t597 * t369 + t598 * t367 + (t495 + (t516 - t558) * t367) * qJD(1) + t385;
t393 = t112 * t538 + t113 * t536 - t9;
t94 = rSges(6,3) * t400 + t471;
t96 = -rSges(6,3) * t399 + t560;
t11 = t367 * t94 + t369 * t96 + (t574 + (-t167 + t516) * t367) * qJD(1) + t385;
t135 = t418 * t367;
t392 = t135 * t538 + t136 * t536 - t11;
t114 = -t169 * t366 + t242 * t564;
t115 = t366 * t167 - t242 * t563;
t30 = t432 * t539 + (qJD(1) * t433 - t367 * t96 + t369 * t94) * t368;
t391 = t114 * t536 + t115 * t538 - t30;
t181 = t482 * t367;
t194 = rSges(5,3) * t564 - t473;
t31 = t367 * (-rSges(5,3) * t505 - t474) + t369 * (-rSges(5,3) * t504 + t556) + (t369 * t194 + (-t193 + t552) * t367) * qJD(1) + t412;
t390 = t181 * t538 + t182 * t536 - t31;
t389 = (rSges(4,2) - pkin(2)) * t368 + t479;
t388 = -(rSges(5,3) * t368 + t366 * t472) * qJD(2) + t398;
t387 = -(t368 * t562 + t528) * qJD(2) + t398;
t384 = (-pkin(7) - t344) * qJD(1) + t421;
t187 = Icges(5,4) * t284 + Icges(5,2) * t283 + Icges(5,6) * t563;
t189 = Icges(5,1) * t284 + Icges(5,4) * t283 + Icges(5,5) * t563;
t383 = t261 / 0.2e1 - t264 / 0.2e1 + t187 * t614 + t189 * t615;
t188 = Icges(5,4) * t286 + Icges(5,2) * t285 + Icges(5,6) * t564;
t190 = Icges(5,1) * t286 + Icges(5,4) * t285 + Icges(5,5) * t564;
t382 = t635 / 0.2e1 + t188 * t614 + t190 * t615 + t260 / 0.2e1;
t177 = (-rSges(6,1) * t350 + rSges(6,2) * t349) * t532 + (rSges(6,3) * t368 + t366 * t469) * qJD(2);
t381 = -t177 + t387;
t376 = t387 - t559;
t375 = qJD(1) * t378;
t374 = qJD(1) * t377;
t373 = t145 * t367 + t147 * t369 + (-t254 * t367 - t256 * t369) * qJD(1);
t101 = t432 * t368;
t24 = (t536 * t555 + t597) * t366 + (qJD(2) * t558 - t369 * t559 + t480) * t368;
t25 = (-t538 * t555 - t598) * t366 + (-qJD(2) * t557 + t367 * t559 + t369 * t493) * t368;
t45 = (t242 * t536 + t96) * t366 + (qJD(2) * t167 - t177 * t369 + t510) * t368;
t46 = (-t242 * t538 - t94) * t366 + (-qJD(2) * t169 + t367 * t177 + t242 * t540) * t368;
t53 = t654 * t368;
t372 = (-qJD(2) * t101 - t114 * t541 + t115 * t540 + t367 * t45 + t369 * t46) * t622 + (qJD(2) * t53 + t24 * t367 + t25 * t369 + t540 * t76 - t541 * t75) * t620;
t40 = t367 * t375 + t504 * t609 + t395 + t646;
t41 = t369 * t375 + ((t476 + (rSges(7,2) - t365) * t366) * qJD(2) + t384) * t367 + t550 + t672;
t51 = t367 * t374 + t504 * t608 + t395 + t560;
t52 = t369 * t374 + ((t476 + (rSges(6,3) - t365) * t366) * qJD(2) + t384) * t367 - t471 + t550;
t72 = qJD(1) * t380 + t504 * t530 + t348 + t486 + t556;
t73 = t325 + t332 + ((-t575 + t593) * qJD(2) + t421) * t367 + ((-pkin(3) - pkin(7)) * t367 + t396 * t369) * qJD(1) + t474;
t371 = (t367 * t40 + t369 * t41 + t540 * t78 - t541 * t77) * t620 + (-t104 * t541 + t105 * t540 + t367 * t51 + t369 * t52) * t622 + (-t125 * t541 + t126 * t540 + t367 * t72 + t369 * t73) * t624;
t102 = t253 * t541 + t369 * t388 + t551;
t103 = qJD(1) * t182 + t367 * t388;
t42 = t411 - t632;
t47 = t369 * t376 + t480 + t515;
t48 = qJD(1) * t113 + t367 * t376;
t50 = t411 - t433;
t62 = t369 * t381 + t510 + t515;
t63 = qJD(1) * t136 + t367 * t381;
t74 = t193 * t369 + t367 * t194 + t488;
t370 = (qJD(2) * t42 + t112 * t540 - t113 * t541 + t367 * t48 + t369 * t47) * t620 + (qJD(2) * t50 + t135 * t540 - t136 * t541 + t367 * t63 + t369 * t62) * t622 + (qJD(2) * t74 + t102 * t369 + t103 * t367 + t181 * t540 - t182 * t541) * t624;
t299 = t475 * qJD(2);
t275 = t367 * t467 - t596;
t272 = t367 * t475 - t594;
t233 = (Icges(5,5) * t368 + t366 * t456) * qJD(2);
t232 = (Icges(5,6) * t368 + t366 * t450) * qJD(2);
t229 = t549 * t367;
t228 = t273 + t545;
t227 = t367 * t420 + t355 + t594;
t226 = t651 * t498;
t208 = qJD(1) * t633 + t369 * t404;
t207 = t367 * t404 + t542;
t198 = t367 * t403 + t543;
t197 = -qJD(1) * t258 + t369 * t403;
t196 = t274 + t484;
t195 = t367 * t389 + t355 + t596;
t153 = t321 * t538 + ((-rSges(3,3) - pkin(7)) * t367 + t420 * t369) * qJD(1);
t152 = -rSges(3,1) * t399 - rSges(3,2) * t502 - pkin(1) * t541 + t347 + t546;
t138 = qJD(1) * t230 + t367 * t553;
t137 = t369 * t553 - t466 * t541 + t291;
t134 = t367 * t423 + t369 * t633;
t133 = -t268 * t369 + t645;
t132 = -t367 * t633 + t642;
t131 = t367 * t268 + t369 * t424;
t130 = t367 * t259 - t369 * t425;
t129 = t367 * t258 - t641;
t128 = -t259 * t369 - t644;
t127 = -t258 * t369 - t367 * t426;
t124 = Icges(5,1) * t223 + Icges(5,4) * t222 - Icges(5,5) * t399;
t123 = Icges(5,1) * t221 + Icges(5,4) * t220 + Icges(5,5) * t400;
t122 = Icges(5,4) * t223 + Icges(5,2) * t222 - Icges(5,6) * t399;
t121 = Icges(5,4) * t221 + Icges(5,2) * t220 + Icges(5,6) * t400;
t116 = t274 * t369 + t367 * t275 + t554;
t111 = t332 + (-t535 + (t368 * t592 - t595) * qJD(2)) * t367 + ((-rSges(4,1) - pkin(7)) * t367 + t389 * t369) * qJD(1);
t110 = -pkin(2) * t504 + (t479 - t606) * t541 + t485 + t513;
t71 = t188 * t285 + t190 * t286 + t520;
t70 = t187 * t285 + t189 * t286 + t522;
t69 = t283 * t188 + t284 * t190 + t519;
t68 = t283 * t187 + t284 * t189 + t521;
t49 = (qJD(1) * t275 + t485) * t369 + (t466 * t538 + (-t274 - t290 + t640) * qJD(1)) * t367 + t517;
t19 = -t147 * t161 + t148 * t165 - t157 * t399 - t254 * t87 + t255 * t91 + t563 * t83;
t18 = -t147 * t160 + t148 * t164 - t156 * t399 - t254 * t88 + t255 * t92 + t563 * t84;
t17 = t147 * t155 + t148 * t163 - t159 * t399 + t254 * t81 + t255 * t89 + t563 * t85;
t16 = t147 * t154 + t148 * t162 - t158 * t399 + t254 * t82 + t255 * t90 + t563 * t86;
t15 = -t145 * t161 + t146 * t165 + t157 * t400 + t256 * t87 + t257 * t91 + t564 * t83;
t14 = -t145 * t160 + t146 * t164 + t156 * t400 + t256 * t88 + t257 * t92 + t564 * t84;
t13 = t145 * t155 + t146 * t163 + t159 * t400 - t256 * t81 + t257 * t89 + t564 * t85;
t12 = t145 * t154 + t146 * t162 + t158 * t400 - t256 * t82 + t257 * t90 + t564 * t86;
t8 = qJD(1) * t463 + t18 * t367 - t19 * t369;
t7 = qJD(1) * t464 + t16 * t367 - t17 * t369;
t6 = qJD(1) * t461 + t14 * t367 - t15 * t369;
t5 = qJD(1) * t462 + t12 * t367 - t13 * t369;
t4 = (-qJD(2) * t463 + t35) * t366 + (-qJD(1) * t37 + qJD(2) * t98 + t18 * t369 + t19 * t367) * t368;
t3 = (-qJD(2) * t464 + t34) * t366 + (-qJD(1) * t36 + qJD(2) * t97 + t16 * t369 + t17 * t367) * t368;
t2 = (-qJD(2) * t461 + t33) * t366 + (-qJD(1) * t39 + qJD(2) * t100 + t14 * t369 + t15 * t367) * t368;
t1 = (-qJD(2) * t462 + t32) * t366 + (-qJD(1) * t38 + qJD(2) * t99 + t12 * t369 + t13 * t367) * t368;
t26 = [(t110 * t196 + t111 * t195) * t630 + (t152 * t228 + t153 * t227) * t631 + (t40 * t78 + t41 * t77) * t627 + (t104 * t52 + t105 * t51) * t628 + (t125 * t73 + t126 * t72) * t629 - t173 * t571 - t233 * t569 - t234 * t501 + t656 * t350 * t532 + (Icges(5,3) * t366 + t441 + t453 + t674) * t537 + (t251 * t364 + t252 * t363 + t366 * t446 + t443 + t458 - t675) * t539 + (Icges(5,3) * t539 - t364 * t232 + t349 * t657 - t446 * t537) * t368 + t655; m(7) * (t112 * t40 + t113 * t41 + t47 * t77 + t48 * t78) + m(6) * (t104 * t62 + t105 * t63 + t135 * t51 + t136 * t52) + m(5) * (t102 * t125 + t103 * t126 + t181 * t72 + t182 * t73) + m(4) * (t110 * t229 + t111 * t230 + t137 * t195 + t138 * t196) + (t251 * t665 + t252 * t664 - t285 * t232 / 0.2e1 - t286 * t233 / 0.2e1 + m(3) * (-t153 * t321 - t227 * t299) + t497 * t369 + (Icges(5,5) * t664 + Icges(5,6) * t665 + t263 * t660 + t266 * t659 + t400 * t666 + t663 * t668) * t366 - t416) * t369 + (t222 * t619 + t223 * t618 + t232 * t617 + t233 * t616 + m(3) * (-t152 * t321 - t228 * t299) + t497 * t367 + (Icges(5,5) * t223 / 0.2e1 + Icges(5,6) * t222 / 0.2e1 + t399 * t666 + t668 * t662 + (t262 + t634) * t660) * t366 + t417) * t367 + ((t264 * t659 + t121 * t364 / 0.2e1 + t123 * t363 / 0.2e1 + t261 * t660 + t669 * t663) * t369 + (t122 * t614 + t124 * t615 + t669 * t662 + (t260 + t635) * t660) * t367) * t368 + ((t367 * t478 + t369 * t477) * t368 + (-t367 * t383 + t369 * t382) * t366) * qJD(2) + ((-t228 * t607 + t251 * t617 + t252 * t616 + t366 * t478 + t368 * t383 - t415) * t369 + (t227 * t607 + t285 * t619 + t286 * t618 - t366 * t477 + t368 * t382 - t414) * t367) * qJD(1); ((t367 * t272 + t273 * t369) * ((qJD(1) * t272 - t369 * t408 + t546) * t369 + (-t367 * t408 + (-t273 + t637) * qJD(1)) * t367) + t544 * t321 * t299) * t631 + (t112 * t48 + t113 * t47 + t42 * t9) * t627 + (t50 * t11 + t135 * t63 + t136 * t62) * t628 + (t102 * t182 + t103 * t181 + t31 * t74) * t629 + (t116 * t49 + t137 * t230 + t138 * t229) * t630 - t369 * t5 - t369 * t6 + t367 * t8 + t367 * t7 - t369 * ((t207 * t369 + (t133 - t642) * qJD(1)) * t369 + (t134 * qJD(1) + (t264 * t537 + t266 * t539 + t542) * t367 + (-t208 + (t366 * t634 + t368 * t635) * qJD(2) + t424 * qJD(1)) * t369) * t367) - t369 * ((t369 * t198 + (t128 + t641) * qJD(1)) * t369 + (t127 * qJD(1) + (-t261 * t537 - t263 * t539 + t543) * t367 + (-t197 + (t260 * t368 + t262 * t366) * qJD(2) - t425 * qJD(1)) * t369) * t367) + t367 * ((t367 * t208 + (t132 - t645) * qJD(1)) * t367 + (t131 * qJD(1) + (t537 * t635 + t539 * t634) * t369 + (-t207 + (t264 * t368 + t266 * t366) * qJD(2) + (t268 + t423) * qJD(1)) * t367) * t369) + t367 * ((t367 * t197 + (t129 + t644) * qJD(1)) * t367 + (t130 * qJD(1) + (t260 * t537 + t262 * t539) * t369 + (-t198 + (-t261 * t368 - t263 * t366) * qJD(2) + (t259 - t426) * qJD(1)) * t367) * t369) - t369 * ((-t285 * t121 - t286 * t123 - t220 * t188 - t221 * t190 + (t70 - t519) * qJD(1)) * t369 + (t285 * t122 + t286 * t124 + t220 * t187 + t221 * t189 + (t71 + t521) * qJD(1)) * t367) + t367 * ((t283 * t122 + t284 * t124 + t222 * t187 + t223 * t189 + (t69 - t522) * qJD(1)) * t367 + (-t283 * t121 - t284 * t123 - t222 * t188 - t223 * t190 + (t68 + t520) * qJD(1)) * t369) + (t38 + t39 + (-t127 - t134 - t71) * t369 + (t128 + t133 + t70) * t367) * t541 + (t36 + t37 + (-t129 - t132 - t69) * t369 + (t130 + t131 + t68) * t367) * t540; 0.2e1 * (t460 * t620 + t439 * t622 + t438 * t624 + (t195 * t369 + t196 * t367) * t626) * t537 + 0.2e1 * ((t110 * t367 + t111 * t369 - t195 * t541 + t196 * t540) * t626 + t371) * t366; 0.2e1 * (t393 * t620 + t392 * t622 + t390 * t624 + (t229 * t538 + t230 * t536 - t49) * t626) * t368 + 0.2e1 * ((qJD(2) * t116 + t137 * t369 + t138 * t367 + t229 * t540 - t230 * t541) * t626 + t370) * t366; 0.4e1 * (t626 + t491) * t226; 0.2e1 * (t438 * t625 + t439 * t623 + t460 * t621) * t539 + 0.2e1 * t371 * t368; 0.2e1 * (t390 * t625 + t392 * t623 + t393 * t621) * t366 + 0.2e1 * t370 * t368; 0.2e1 * t491 * t639 * t651 * qJD(2); -0.4e1 * t491 * t226; m(7) * (t24 * t78 + t25 * t77 + t40 * t76 + t41 * t75) + m(6) * (t104 * t46 + t105 * t45 + t114 * t52 + t115 * t51) + (t367 * t414 + t369 * t415) * t539 + (t417 * t369 + t416 * t367 + (t367 * t415 - t369 * t414) * qJD(1)) * t368 + t492; m(7) * (t10 * t42 + t112 * t24 + t113 * t25 + t47 * t75 + t48 * t76 + t53 * t9) + m(6) * (-t101 * t11 + t114 * t62 + t115 * t63 + t135 * t45 + t136 * t46 + t30 * t50) + (-t2 / 0.2e1 - t1 / 0.2e1 + t524 * t539) * t369 + (t4 / 0.2e1 + t3 / 0.2e1 - t523 * t539) * t367 + ((t367 * t524 + t369 * t523) * qJD(1) + (t5 + t6) * t663 + (t7 + t8) * t369 / 0.2e1 + t653 * t650) * t368 + (qJD(1) * t652 + t649 * t367 - t648 * t369) * t366 / 0.2e1 + (t601 * t367 + t647 * t369) * t659; 0.2e1 * (t391 * t622 + t394 * t620) * t368 + 0.2e1 * t372 * t366; 0.2e1 * (t391 * t623 + t394 * t621) * t366 + 0.2e1 * t372 * t368; (t53 * t10 + t24 * t76 + t25 * t75) * t627 + (-t101 * t30 + t114 * t46 + t115 * t45) * t628 + (((-t366 * t600 - t647) * t369 + (-t366 * t599 - t601) * t367) * qJD(2) + t492) * t366 + ((t3 + t4) * t369 + (t1 + t2) * t367 + (t648 * t367 + t649 * t369) * t366 + (t366 * t658 + t368 * t652) * qJD(2) + (-t366 * t653 - t367 * t647 + t601 * t369) * qJD(1)) * t368; m(7) * (t145 * t78 + t147 * t77 + t254 * t41 - t256 * t40); m(7) * (-t42 * t501 + t112 * t145 + t113 * t147 + t254 * t47 - t256 * t48 + (t368 * t9 - t42 * t539) * t350); m(7) * (t499 + t427 * t537 + (0.2e1 * t350 * t537 + t373) * t366); m(7) * ((t350 * t639 - t427 * t366) * qJD(2) + (-t349 * t533 + t373) * t368); m(7) * (-t53 * t501 + t145 * t76 + t147 * t75 - t24 * t256 + t25 * t254 + (t10 * t368 - t53 * t539) * t350); (-t145 * t256 + t147 * t254 + (-t350 * t498 - t499) * t350) * t627;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t26(1) t26(2) t26(4) t26(7) t26(11) t26(16); t26(2) t26(3) t26(5) t26(8) t26(12) t26(17); t26(4) t26(5) t26(6) t26(9) t26(13) t26(18); t26(7) t26(8) t26(9) t26(10) t26(14) t26(19); t26(11) t26(12) t26(13) t26(14) t26(15) t26(20); t26(16) t26(17) t26(18) t26(19) t26(20) t26(21);];
Mq  = res;
