% Calculate time derivative of joint inertia matrix for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP10_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:29
% EndTime: 2019-03-09 06:31:11
% DurationCPUTime: 29.11s
% Computational Cost: add. (45314->1103), mult. (68484->1487), div. (0->0), fcn. (66231->8), ass. (0->519)
t466 = sin(qJ(3));
t469 = cos(qJ(3));
t464 = qJ(4) + qJ(5);
t454 = sin(t464);
t455 = cos(t464);
t525 = Icges(6,5) * t455 - Icges(6,6) * t454;
t366 = Icges(6,3) * t466 + t469 * t525;
t529 = Icges(7,4) * t455 + Icges(7,6) * t454;
t367 = Icges(7,2) * t466 + t469 * t529;
t761 = t367 + t366;
t745 = rSges(7,1) + pkin(5);
t694 = rSges(7,3) + qJ(6);
t470 = cos(qJ(1));
t660 = t466 * t470;
t467 = sin(qJ(1));
t667 = t455 * t467;
t392 = t454 * t660 + t667;
t654 = t470 * t455;
t607 = t466 * t654;
t670 = t454 * t467;
t393 = -t607 + t670;
t656 = t469 * t470;
t641 = rSges(7,2) * t656 - t392 * t694 + t393 * t745;
t461 = qJD(4) + qJD(5);
t621 = qJD(1) * t470;
t624 = qJD(1) * t466;
t559 = t461 + t624;
t617 = qJD(3) * t470;
t578 = t469 * t617;
t725 = t559 * t467 - t578;
t245 = t454 * t725 - t455 * t621 - t461 * t607;
t560 = t461 * t466 + qJD(1);
t668 = t454 * t470;
t246 = t455 * t725 + t560 * t668;
t685 = Icges(7,5) * t455;
t524 = Icges(7,3) * t454 + t685;
t665 = t461 * t469;
t686 = Icges(7,5) * t454;
t264 = (Icges(7,3) * t455 - t686) * t665 + (Icges(7,6) * t469 - t466 * t524) * qJD(3);
t266 = (-Icges(7,4) * t454 + Icges(7,6) * t455) * t665 + (Icges(7,2) * t469 - t466 * t529) * qJD(3);
t534 = Icges(7,1) * t455 + t686;
t268 = (-Icges(7,1) * t454 + t685) * t665 + (Icges(7,4) * t469 - t466 * t534) * qJD(3);
t365 = Icges(7,6) * t466 + t469 * t524;
t369 = Icges(7,4) * t466 + t469 * t534;
t581 = t466 * t617;
t622 = qJD(1) * t469;
t586 = t467 * t622;
t732 = t581 + t586;
t80 = t245 * t365 + t246 * t369 - t264 * t392 + t266 * t656 + t268 * t393 - t367 * t732;
t265 = (-Icges(6,5) * t454 - Icges(6,6) * t455) * t665 + (Icges(6,3) * t469 - t466 * t525) * qJD(3);
t688 = Icges(6,4) * t455;
t530 = -Icges(6,2) * t454 + t688;
t689 = Icges(6,4) * t454;
t267 = (-Icges(6,2) * t455 - t689) * t665 + (Icges(6,6) * t469 - t466 * t530) * qJD(3);
t535 = Icges(6,1) * t455 - t689;
t269 = (-Icges(6,1) * t454 - t688) * t665 + (Icges(6,5) * t469 - t466 * t535) * qJD(3);
t368 = Icges(6,6) * t466 + t469 * t530;
t370 = Icges(6,5) * t466 + t469 * t535;
t81 = -t245 * t368 + t246 * t370 + t265 * t656 + t267 * t392 + t269 * t393 - t366 * t732;
t760 = t80 + t81;
t618 = qJD(3) * t469;
t579 = t467 * t618;
t477 = t470 * t559 + t579;
t247 = t454 * t477 + t560 * t667;
t248 = t455 * t477 - t560 * t670;
t661 = t466 * t467;
t390 = t454 * t661 - t654;
t391 = t455 * t661 + t668;
t619 = qJD(3) * t467;
t582 = t466 * t619;
t585 = t469 * t621;
t485 = t582 - t585;
t658 = t467 * t469;
t82 = t247 * t365 + t248 * t369 + t264 * t390 - t266 * t658 + t268 * t391 + t367 * t485;
t83 = -t247 * t368 + t248 * t370 - t265 * t658 - t267 * t390 + t269 * t391 + t366 * t485;
t759 = t82 + t83;
t274 = Icges(7,5) * t391 - Icges(7,6) * t658 + Icges(7,3) * t390;
t278 = Icges(7,4) * t391 - Icges(7,2) * t658 + Icges(7,6) * t390;
t282 = Icges(7,1) * t391 - Icges(7,4) * t658 + Icges(7,5) * t390;
t123 = t274 * t390 - t278 * t658 + t282 * t391;
t275 = Icges(7,5) * t393 + Icges(7,6) * t656 - Icges(7,3) * t392;
t279 = Icges(7,4) * t393 + Icges(7,2) * t656 - Icges(7,6) * t392;
t283 = Icges(7,1) * t393 + Icges(7,4) * t656 - Icges(7,5) * t392;
t124 = t275 * t390 - t279 * t658 + t283 * t391;
t276 = Icges(6,5) * t391 - Icges(6,6) * t390 - Icges(6,3) * t658;
t280 = Icges(6,4) * t391 - Icges(6,2) * t390 - Icges(6,6) * t658;
t284 = Icges(6,1) * t391 - Icges(6,4) * t390 - Icges(6,5) * t658;
t125 = -t276 * t658 - t280 * t390 + t284 * t391;
t277 = Icges(6,5) * t393 + Icges(6,6) * t392 + Icges(6,3) * t656;
t281 = Icges(6,4) * t393 + Icges(6,2) * t392 + Icges(6,6) * t656;
t285 = Icges(6,1) * t393 + Icges(6,4) * t392 + Icges(6,5) * t656;
t126 = -t277 * t658 - t281 * t390 + t285 * t391;
t739 = (t123 + t125) * t470 + (t124 + t126) * t467;
t127 = -t274 * t392 + t278 * t656 + t282 * t393;
t128 = -t275 * t392 + t279 * t656 + t283 * t393;
t129 = t276 * t656 + t280 * t392 + t284 * t393;
t130 = t277 * t656 + t281 * t392 + t285 * t393;
t738 = (t127 + t129) * t470 + (t128 + t130) * t467;
t513 = t274 * t454 + t282 * t455;
t141 = t278 * t466 + t469 * t513;
t511 = t280 * t454 - t284 * t455;
t143 = t276 * t466 - t469 * t511;
t758 = t141 + t143;
t512 = t275 * t454 + t283 * t455;
t142 = t279 * t466 + t469 * t512;
t510 = t281 * t454 - t285 * t455;
t144 = t277 * t466 - t469 * t510;
t757 = t142 + t144;
t189 = t365 * t390 - t367 * t658 + t369 * t391;
t190 = -t366 * t658 - t368 * t390 + t370 * t391;
t756 = t189 + t190;
t191 = -t365 * t392 + t367 * t656 + t369 * t393;
t192 = t366 * t656 + t368 * t392 + t370 * t393;
t755 = t191 + t192;
t620 = qJD(3) * t466;
t583 = t455 * t620;
t584 = t454 * t620;
t608 = t455 * t665;
t609 = t454 * t665;
t669 = t454 * t469;
t754 = t264 * t669 + t365 * t608 + t368 * t584 - t369 * t609 - t370 * t583 + (t268 + t269) * t455 * t469 + t761 * t618 + (t265 + t266) * t466;
t520 = t129 * t467 - t130 * t470;
t521 = t127 * t467 - t128 * t470;
t753 = t520 + t521;
t522 = t125 * t467 - t126 * t470;
t523 = t123 * t467 - t124 * t470;
t752 = t522 + t523;
t504 = t365 * t454 + t369 * t455;
t751 = ((-t368 * t454 + t370 * t455 + t504) * t469 + t761 * t466) * t618;
t750 = t392 * qJD(6) - t694 * t245 - t745 * t246;
t736 = -rSges(7,2) * t658 + t694 * t390 + t391 * t745;
t569 = -t620 / 0.2e1;
t571 = -t622 / 0.2e1;
t749 = t467 * t571 + t470 * t569;
t568 = t620 / 0.2e1;
t748 = t467 * t568 + t470 * t571;
t163 = Icges(7,5) * t248 + Icges(7,6) * t485 + Icges(7,3) * t247;
t167 = Icges(7,4) * t248 + Icges(7,2) * t485 + Icges(7,6) * t247;
t171 = Icges(7,1) * t248 + Icges(7,4) * t485 + Icges(7,5) * t247;
t34 = -t163 * t392 + t167 * t656 + t171 * t393 + t245 * t274 + t246 * t282 - t278 * t732;
t162 = Icges(7,5) * t246 - Icges(7,6) * t732 + Icges(7,3) * t245;
t166 = Icges(7,4) * t246 - Icges(7,2) * t732 + Icges(7,6) * t245;
t170 = Icges(7,1) * t246 - Icges(7,4) * t732 + Icges(7,5) * t245;
t35 = -t162 * t392 + t166 * t656 + t170 * t393 + t245 * t275 + t246 * t283 - t279 * t732;
t165 = Icges(6,5) * t248 - Icges(6,6) * t247 + Icges(6,3) * t485;
t169 = Icges(6,4) * t248 - Icges(6,2) * t247 + Icges(6,6) * t485;
t173 = Icges(6,1) * t248 - Icges(6,4) * t247 + Icges(6,5) * t485;
t36 = t165 * t656 + t169 * t392 + t173 * t393 - t245 * t280 + t246 * t284 - t276 * t732;
t164 = Icges(6,5) * t246 - Icges(6,6) * t245 - Icges(6,3) * t732;
t168 = Icges(6,4) * t246 - Icges(6,2) * t245 - Icges(6,6) * t732;
t172 = Icges(6,1) * t246 - Icges(6,4) * t245 - Icges(6,5) * t732;
t37 = t164 * t656 + t168 * t392 + t172 * t393 - t245 * t281 + t246 * t285 - t277 * t732;
t746 = ((t35 + t37) * t470 + (-t34 - t36) * t467 + t755 * qJD(3) - t738 * qJD(1)) * t469 + (qJD(3) * t753 + t760) * t466;
t38 = t163 * t390 - t167 * t658 + t171 * t391 + t247 * t274 + t248 * t282 + t278 * t485;
t39 = t162 * t390 - t166 * t658 + t170 * t391 + t247 * t275 + t248 * t283 + t279 * t485;
t40 = -t165 * t658 - t169 * t390 + t173 * t391 - t247 * t280 + t248 * t284 + t276 * t485;
t41 = -t164 * t658 - t168 * t390 + t172 * t391 - t247 * t281 + t248 * t285 + t277 * t485;
t744 = ((t39 + t41) * t470 + (-t38 - t40) * t467 + t756 * qJD(3) - t739 * qJD(1)) * t469 + (qJD(3) * t752 + t759) * t466;
t46 = (-qJD(3) * t513 + t167) * t466 + (qJD(3) * t278 + (t274 * t461 + t171) * t455 + (-t282 * t461 + t163) * t454) * t469;
t48 = (qJD(3) * t511 + t165) * t466 + (qJD(3) * t276 + (-t280 * t461 + t173) * t455 + (-t284 * t461 - t169) * t454) * t469;
t743 = t46 + t48;
t47 = (-qJD(3) * t512 + t166) * t466 + (qJD(3) * t279 + (t275 * t461 + t170) * t455 + (-t283 * t461 + t162) * t454) * t469;
t49 = (qJD(3) * t510 + t164) * t466 + (qJD(3) * t277 + (-t281 * t461 + t172) * t455 + (-t285 * t461 - t168) * t454) * t469;
t742 = t47 + t49;
t741 = t466 * t756 - t469 * t752;
t740 = t466 * t755 - t469 * t753;
t737 = t485 * rSges(7,2) + t390 * qJD(6) + t694 * t247 + t248 * t745;
t539 = rSges(7,1) * t455 + rSges(7,3) * t454;
t644 = (-pkin(5) * t620 + qJ(6) * t665) * t455 + (-qJ(6) * t620 + (-pkin(5) * t461 + qJD(6)) * t469) * t454 + (-rSges(7,1) * t454 + rSges(7,3) * t455) * t665 + (rSges(7,2) * t469 - t466 * t539) * qJD(3);
t631 = rSges(7,2) * t466 + (pkin(5) * t455 + qJ(6) * t454 + t539) * t469;
t735 = t467 * t757 + t470 * t758;
t734 = t467 * t758 - t470 * t757;
t733 = t466 * t621 + t579;
t606 = t461 * t455 * t368;
t731 = t751 + ((-t606 + (-t370 * t461 - t267) * t454) * t469 - t504 * t620 + t754) * t466;
t465 = sin(qJ(4));
t468 = cos(qJ(4));
t690 = Icges(5,4) * t468;
t531 = -Icges(5,2) * t465 + t690;
t381 = Icges(5,6) * t466 + t469 * t531;
t691 = Icges(5,4) * t465;
t536 = Icges(5,1) * t468 - t691;
t384 = Icges(5,5) * t466 + t469 * t536;
t730 = -t381 * t465 + t384 * t468;
t451 = pkin(4) * t468 + pkin(3);
t704 = pkin(3) - t451;
t727 = t466 * t704;
t692 = Icges(4,4) * t469;
t537 = Icges(4,1) * t466 + t692;
t385 = Icges(4,5) * t470 + t467 * t537;
t674 = t385 * t466;
t693 = Icges(4,4) * t466;
t532 = Icges(4,2) * t469 + t693;
t382 = Icges(4,6) * t470 + t467 * t532;
t678 = t382 * t469;
t503 = t674 + t678;
t726 = t470 * t503;
t545 = rSges(4,1) * t466 + rSges(4,2) * t469;
t492 = t470 * t545;
t563 = t631 * t467;
t556 = qJD(4) + t624;
t724 = t467 * t556 - t578;
t527 = Icges(4,5) * t466 + Icges(4,6) * t469;
t723 = -Icges(4,3) * t467 + t470 * t527;
t722 = -Icges(4,6) * t467 + t470 * t532;
t721 = -Icges(4,5) * t467 + t470 * t537;
t720 = 2 * m(4);
t719 = 2 * m(5);
t718 = 2 * m(6);
t717 = 2 * m(7);
t462 = t467 ^ 2;
t463 = t470 ^ 2;
t716 = -pkin(1) - pkin(7);
t715 = -t466 / 0.2e1;
t714 = t466 / 0.2e1;
t713 = -t467 / 0.2e1;
t712 = t467 / 0.2e1;
t711 = t469 / 0.2e1;
t710 = t470 / 0.2e1;
t708 = rSges(3,2) - pkin(1);
t700 = rSges(4,2) * t466;
t426 = rSges(4,1) * t469 - t700;
t706 = m(4) * t426;
t705 = pkin(3) * t466;
t471 = -pkin(9) - pkin(8);
t703 = -pkin(8) - t471;
t699 = rSges(5,3) * t469;
t698 = pkin(4) * qJD(4);
t697 = t467 * rSges(4,3);
t458 = t470 * rSges(4,3);
t696 = rSges(7,2) - t471;
t695 = rSges(6,3) - t471;
t659 = t467 * t468;
t402 = t465 * t660 + t659;
t657 = t468 * t470;
t663 = t465 * t467;
t403 = -t466 * t657 + t663;
t543 = -rSges(5,1) * t403 - rSges(5,2) * t402;
t320 = rSges(5,3) * t656 - t543;
t682 = t320 * t470;
t680 = t381 * t468;
t679 = t382 * t466;
t677 = t722 * t466;
t676 = t722 * t469;
t673 = t385 * t469;
t672 = t721 * t466;
t671 = t721 * t469;
t615 = qJD(4) * t469;
t327 = (-Icges(5,2) * t468 - t691) * t615 + (Icges(5,6) * t469 - t466 * t531) * qJD(3);
t664 = t465 * t327;
t662 = t465 * t470;
t655 = t469 * t471;
t526 = Icges(5,5) * t468 - Icges(5,6) * t465;
t378 = Icges(5,3) * t466 + t469 * t526;
t224 = t378 * t466 + t730 * t469;
t324 = (-Icges(5,5) * t465 - Icges(5,6) * t468) * t615 + (Icges(5,3) * t469 - t466 * t526) * qJD(3);
t330 = (-Icges(5,1) * t465 - t690) * t615 + (Icges(5,5) * t469 - t466 * t536) * qJD(3);
t479 = t469 * t468 * t330 + t466 * t324 + t378 * t618 - t730 * t620;
t651 = ((-t664 + (-t384 * t465 - t680) * qJD(4)) * t469 + t479) * t466 + t224 * t618;
t650 = -rSges(7,2) * t732 - t750;
t541 = t246 * rSges(6,1) - t245 * rSges(6,2);
t175 = -rSges(6,3) * t732 + t541;
t411 = t451 * t578;
t448 = pkin(4) * t662;
t616 = qJD(3) * t471;
t580 = t466 * t616;
t591 = pkin(3) * t578 + t732 * pkin(8);
t204 = t470 * t580 - t411 + t402 * t698 + (t448 + (t655 - t727) * t467) * qJD(1) + t591;
t648 = -t175 - t204;
t177 = t248 * rSges(6,1) - t247 * rSges(6,2) + t485 * rSges(6,3);
t592 = t733 * pkin(3) + pkin(8) * t582;
t483 = pkin(8) * t585 - t592;
t557 = qJD(4) * t466 + qJD(1);
t497 = t557 * t465;
t610 = t468 * t698;
t551 = t733 * t451 + t470 * t610 + t471 * t585;
t205 = (-pkin(4) * t497 - t580) * t467 + t483 + t551;
t647 = -t177 - t205;
t447 = pkin(3) * t661;
t335 = t470 * (qJD(1) * t447 - t591);
t646 = t470 * t204 + t335;
t645 = t644 * t656;
t290 = t391 * rSges(6,1) - t390 * rSges(6,2) - rSges(6,3) * t658;
t405 = -pkin(8) * t658 + t447;
t594 = t451 * t661 + t467 * t655 + t448;
t321 = -t405 + t594;
t642 = -t290 - t321;
t292 = rSges(6,1) * t393 + rSges(6,2) * t392 + rSges(6,3) * t656;
t449 = pkin(3) * t660;
t406 = pkin(8) * t656 - t449;
t628 = t451 * t660 + t470 * t655;
t322 = pkin(4) * t663 - t406 - t628;
t640 = -t292 - t322;
t364 = t466 * t703 - t469 * t704;
t639 = t466 * t321 + t364 * t658;
t397 = t470 * t406;
t638 = t470 * t322 + t397;
t400 = -t465 * t661 + t657;
t401 = t466 * t659 + t662;
t630 = t401 * rSges(5,1) + t400 * rSges(5,2);
t319 = -rSges(5,3) * t658 + t630;
t637 = -t319 - t405;
t577 = t465 * t615;
t334 = -pkin(4) * t577 + (t469 * t703 + t727) * qJD(3);
t418 = (pkin(8) * t469 - t705) * qJD(3);
t636 = -t334 - t418;
t540 = rSges(6,1) * t455 - rSges(6,2) * t454;
t372 = rSges(6,3) * t466 + t469 * t540;
t219 = t466 * t290 + t372 * t658;
t635 = t631 * t656;
t428 = t469 * pkin(3) + t466 * pkin(8);
t623 = qJD(1) * t467;
t408 = t428 * t623;
t634 = t364 * t623 + t408;
t416 = t467 * t428;
t633 = t467 * t364 + t416;
t632 = -t364 - t372;
t629 = t467 * t418 + t428 * t621;
t627 = qJ(2) * t621 + qJD(2) * t467;
t626 = t470 * pkin(1) + t467 * qJ(2);
t379 = Icges(4,3) * t470 + t467 * t527;
t625 = qJD(1) * t379;
t613 = -rSges(4,3) + t716;
t611 = t465 * t698;
t304 = -t557 * t659 + (-t470 * t556 - t579) * t465;
t305 = t556 * t657 + (t468 * t618 - t497) * t467;
t102 = t304 * t381 + t305 * t384 - t324 * t658 + t327 * t400 + t330 * t401 + t378 * t485;
t194 = Icges(5,5) * t305 + Icges(5,6) * t304 + Icges(5,3) * t485;
t196 = Icges(5,4) * t305 + Icges(5,2) * t304 + Icges(5,6) * t485;
t198 = Icges(5,1) * t305 + Icges(5,4) * t304 + Icges(5,5) * t485;
t313 = Icges(5,5) * t401 + Icges(5,6) * t400 - Icges(5,3) * t658;
t315 = Icges(5,4) * t401 + Icges(5,2) * t400 - Icges(5,6) * t658;
t317 = Icges(5,1) * t401 + Icges(5,4) * t400 - Icges(5,5) * t658;
t509 = t315 * t465 - t317 * t468;
t56 = (qJD(3) * t509 + t194) * t466 + (qJD(3) * t313 - t196 * t465 + t198 * t468 + (-t315 * t468 - t317 * t465) * qJD(4)) * t469;
t605 = -t102 / 0.2e1 - t56 / 0.2e1;
t498 = t470 * t557;
t302 = -t465 * t724 + t468 * t498;
t303 = t465 * t498 + t468 * t724;
t101 = t302 * t381 + t303 * t384 + t324 * t656 + t327 * t402 + t330 * t403 - t378 * t732;
t193 = Icges(5,5) * t303 + Icges(5,6) * t302 - Icges(5,3) * t732;
t195 = Icges(5,4) * t303 + Icges(5,2) * t302 - Icges(5,6) * t732;
t197 = Icges(5,1) * t303 + Icges(5,4) * t302 - Icges(5,5) * t732;
t314 = Icges(5,5) * t403 + Icges(5,6) * t402 + Icges(5,3) * t656;
t316 = Icges(5,4) * t403 + Icges(5,2) * t402 + Icges(5,6) * t656;
t318 = Icges(5,1) * t403 + Icges(5,4) * t402 + Icges(5,5) * t656;
t508 = t316 * t465 - t318 * t468;
t57 = (qJD(3) * t508 + t193) * t466 + (qJD(3) * t314 - t195 * t465 + t197 * t468 + (-t316 * t468 - t318 * t465) * qJD(4)) * t469;
t604 = t57 / 0.2e1 + t101 / 0.2e1;
t603 = -t204 - t650;
t602 = -t205 - t737;
t601 = t732 * t290 + t292 * t582;
t600 = -t321 - t736;
t599 = -t405 + t642;
t598 = -t322 - t641;
t597 = t732 * t321 + t322 * t582;
t596 = t305 * rSges(5,1) + t304 * rSges(5,2) + rSges(5,3) * t582;
t595 = -t364 - t631;
t593 = t733 * rSges(4,1) + rSges(4,2) * t585;
t387 = rSges(4,1) * t661 + rSges(4,2) * t658 + t458;
t590 = t470 * pkin(7) + t626;
t589 = t469 * (-rSges(5,3) - pkin(8));
t542 = rSges(5,1) * t468 - rSges(5,2) * t465;
t388 = rSges(5,3) * t466 + t469 * t542;
t588 = t388 * t623;
t353 = t364 * t621;
t576 = -t658 / 0.2e1;
t575 = t656 / 0.2e1;
t160 = t313 * t466 - t469 * t509;
t206 = -t378 * t658 + t381 * t400 + t384 * t401;
t574 = t160 / 0.2e1 + t206 / 0.2e1;
t161 = t314 * t466 - t469 * t508;
t207 = t378 * t656 + t381 * t402 + t384 * t403;
t573 = -t161 / 0.2e1 - t207 / 0.2e1;
t147 = -t313 * t658 + t315 * t400 + t317 * t401;
t148 = -t314 * t658 + t316 * t400 + t318 * t401;
t515 = t147 * t467 - t148 * t470;
t73 = t466 * t206 - t469 * t515;
t566 = t160 * t466 + t73;
t417 = t545 * qJD(3);
t565 = t417 * (t462 + t463);
t564 = t632 * t467;
t562 = -t451 * t466 - qJ(2);
t561 = qJD(1) * t631;
t558 = -pkin(4) * t465 + t716;
t273 = (-rSges(6,1) * t454 - rSges(6,2) * t455) * t665 + (rSges(6,3) * t469 - t466 * t540) * qJD(3);
t555 = t466 * t177 + t273 * t658 + t290 * t618 + t372 * t585;
t554 = t466 * t205 + t321 * t618 + t334 * t658 + t469 * t353;
t553 = -t405 + t600;
t552 = t467 * t334 + t353 + t629;
t153 = t736 * t466 + t631 * t658;
t546 = t595 * t467;
t544 = rSges(5,1) * t303 + rSges(5,2) * t302;
t538 = Icges(4,1) * t469 - t693;
t533 = -Icges(4,2) * t466 + t692;
t528 = Icges(4,5) * t469 - Icges(4,6) * t466;
t105 = t147 * t470 + t148 * t467;
t149 = t313 * t656 + t315 * t402 + t317 * t403;
t150 = t314 * t656 + t316 * t402 + t318 * t403;
t106 = t149 * t470 + t150 * t467;
t514 = t149 * t467 - t150 * t470;
t507 = t319 * t470 + t320 * t467;
t502 = -t672 - t676;
t74 = t466 * t207 - t469 * t514;
t501 = -t161 * t466 - t74 - t740;
t499 = t558 * t470;
t496 = t641 * t582 + t736 * t732;
t495 = t590 + t594;
t494 = rSges(3,3) * t470 + t467 * t708;
t333 = (-rSges(5,1) * t465 - rSges(5,2) * t468) * t615 + (-t466 * t542 + t699) * qJD(3);
t493 = t333 * t467 + t388 * t621;
t453 = qJD(2) * t470;
t491 = -t467 * t610 + t411 + t453;
t490 = t746 * t656 + t741 * t582 - t734 * t618 * t469 + (t734 * t620 + t751 + (-t735 * qJD(1) - t743 * t467 + t742 * t470) * t469 + t731) * t466;
t489 = t502 * t467;
t488 = qJD(3) * t538;
t487 = qJD(3) * t533;
t486 = -t584 + t608;
t482 = t737 * t466 + t631 * t585 + t736 * t618 + t644 * t658;
t481 = t467 * t716 + t470 * t589;
t457 = t470 * qJ(2);
t478 = t467 * t558 + t457 + t628;
t19 = -qJD(1) * t521 + t34 * t470 + t35 * t467;
t20 = -qJD(1) * t520 + t36 * t470 + t37 * t467;
t21 = -qJD(1) * t523 + t38 * t470 + t39 * t467;
t22 = -qJD(1) * t522 + t40 * t470 + t41 * t467;
t475 = (-t734 * qJD(1) + t742 * t467 + t743 * t470) * t714 + t746 * t712 + t744 * t710 + (t21 + t22) * t576 + (t19 + t20) * t575 - t741 * t623 / 0.2e1 + t740 * t621 / 0.2e1 + t735 * t618 / 0.2e1 + t748 * t739 + t749 * t738;
t474 = ((-t611 - t616) * t466 + t558 * qJD(1)) * t467 + t551 + t627;
t473 = (t743 + t759) * t576 + (t742 + t760) * t575 + t731 + t748 * (t756 + t758) + t749 * (t755 + t757);
t472 = (-t744 * t467 + (-t467 * t740 - t470 * t741) * qJD(1)) * t469 - t740 * t581 + t490;
t413 = t527 * qJD(3);
t396 = -rSges(3,2) * t470 + rSges(3,3) * t467 + t626;
t395 = t457 + t494;
t389 = t697 - t492;
t355 = t453 + (t708 * t470 + (-rSges(3,3) - qJ(2)) * t467) * qJD(1);
t354 = qJD(1) * t494 + t627;
t351 = t372 * t656;
t347 = t364 * t656;
t342 = t590 + t387;
t341 = t467 * t613 + t457 + t492;
t337 = (-t388 - t428) * t470;
t336 = t388 * t467 + t416;
t326 = qJD(1) * t723 + t528 * t619;
t325 = -t528 * t617 + t625;
t310 = t334 * t656;
t252 = t273 * t656;
t234 = t453 + t426 * t617 + (t613 * t470 + (-qJ(2) - t545) * t467) * qJD(1);
t233 = (-rSges(4,2) * t620 + qJD(1) * t613) * t467 + t593 + t627;
t232 = (-t428 + t632) * t470;
t231 = t372 * t467 + t633;
t230 = t467 * t589 + t447 + t590 + t630;
t229 = t449 + t457 + t481 + t543;
t228 = -t320 * t466 + t388 * t656;
t227 = t319 * t466 + t388 * t658;
t226 = -t467 * t723 - t470 * t502;
t225 = t379 * t467 - t726;
t223 = -t470 * t723 + t489;
t222 = t379 * t470 + t467 * t503;
t220 = -t292 * t466 + t351;
t218 = (-t428 + t595) * t470;
t217 = t633 + t563;
t216 = t290 + t495;
t215 = -t292 + t478;
t210 = t507 * t469;
t209 = t493 + t629;
t208 = t588 + t408 + (-t333 - t418) * t470;
t201 = (-t290 * t470 - t292 * t467) * t469;
t200 = -rSges(5,3) * t585 + t596;
t199 = -rSges(5,3) * t732 + t544;
t180 = t467 * t637 + t397 + t682;
t179 = t495 + t736;
t178 = t478 - t641;
t154 = -t466 * t641 + t635;
t152 = t466 * t640 + t347 + t351;
t151 = t219 + t639;
t146 = rSges(5,3) * t581 + t453 + (t716 * t470 + (-qJ(2) + t699 - t705) * t467) * qJD(1) - t544 + t591;
t145 = qJD(1) * t481 + t592 + t596 + t627;
t132 = t273 * t467 + t372 * t621 + t552;
t131 = t372 * t623 + (-t273 + t636) * t470 + t634;
t122 = (t467 * t640 + t470 * t642) * t469;
t121 = (-t467 * t641 - t470 * t736) * t469;
t120 = t466 * t598 + t347 + t635;
t119 = t153 + t639;
t118 = t292 * t470 + t467 * t599 + t638;
t117 = (qJD(3) * t695 - t611) * t660 + (t499 + (t469 * t695 + t562) * t467) * qJD(1) + t491 - t541;
t116 = t177 + t474;
t115 = (-t388 * t619 + t200) * t466 + (qJD(3) * t319 + t493) * t469;
t114 = (-t388 * t617 - t199) * t466 + (-qJD(3) * t320 + t333 * t470 - t588) * t469;
t113 = t467 * t644 + t470 * t561 + t552;
t112 = t467 * t561 + (t636 - t644) * t470 + t634;
t110 = (t467 * t598 + t470 * t600) * t469;
t109 = -t372 * t582 + t555;
t108 = -t372 * t586 - t175 * t466 + t252 + (-t292 * t469 - t372 * t660) * qJD(3);
t107 = t467 * t553 + t470 * t641 + t638;
t98 = (qJD(3) * t696 - t611) * t660 + (t499 + (t469 * t696 + t562) * t467) * qJD(1) + t491 + t750;
t97 = t474 + t737;
t84 = t507 * t620 + (-t199 * t467 - t200 * t470 + (t319 * t467 - t682) * qJD(1)) * t469;
t79 = t199 * t470 + t335 + (-t200 + t483) * t467 + (t637 * t470 + (-t320 - t406) * t467) * qJD(1);
t70 = (-t175 * t467 + (-qJD(1) * t292 - t177) * t470) * t469 + t601;
t59 = t564 * t620 + t554 + t555;
t58 = t252 + t310 + t648 * t466 + t564 * t622 + (t469 * t640 + t632 * t660) * qJD(3);
t55 = -t563 * t620 + t482;
t54 = -t650 * t466 - t563 * t622 + (-t469 * t641 - t631 * t660) * qJD(3) + t645;
t53 = -t193 * t658 + t195 * t400 + t197 * t401 + t304 * t316 + t305 * t318 + t314 * t485;
t52 = -t194 * t658 + t196 * t400 + t198 * t401 + t304 * t315 + t305 * t317 + t313 * t485;
t51 = t193 * t656 + t195 * t402 + t197 * t403 + t302 * t316 + t303 * t318 - t314 * t732;
t50 = t194 * t656 + t196 * t402 + t198 * t403 + t302 * t315 + t303 * t317 - t313 * t732;
t33 = t175 * t470 + (t483 + t647) * t467 + (t599 * t470 + (-t406 + t640) * t467) * qJD(1) + t646;
t32 = t546 * t620 + t482 + t554;
t31 = t310 + t603 * t466 + t546 * t622 + (t469 * t598 + t595 * t660) * qJD(3) + t645;
t30 = (t648 * t467 + (qJD(1) * t640 + t647) * t470) * t469 + t597 + t601;
t29 = (-t650 * t467 + (-qJD(1) * t641 - t737) * t470) * t469 + t496;
t28 = t650 * t470 + (t483 + t602) * t467 + (t553 * t470 + (-t406 + t598) * t467) * qJD(1) + t646;
t27 = (t603 * t467 + (qJD(1) * t598 + t602) * t470) * t469 + t496 + t597;
t26 = -qJD(1) * t515 + t467 * t53 + t470 * t52;
t25 = -qJD(1) * t514 + t467 * t51 + t470 * t50;
t14 = (qJD(3) * t515 + t102) * t466 + (-qJD(1) * t105 + qJD(3) * t206 - t467 * t52 + t470 * t53) * t469;
t13 = (qJD(3) * t514 + t101) * t466 + (-qJD(1) * t106 + qJD(3) * t207 - t467 * t50 + t470 * t51) * t469;
t1 = [-t466 * t488 + t532 * t620 - t384 * t577 - t615 * t680 - t370 * t609 - t369 * t583 - t365 * t584 + t479 - t537 * t618 - t267 * t669 + 0.2e1 * m(3) * (t354 * t396 + t355 * t395) + (t233 * t342 + t234 * t341) * t720 + (t145 * t230 + t146 * t229) * t719 + (t116 * t216 + t117 * t215) * t718 + (t178 * t98 + t179 * t97) * t717 + (-t606 - t487 - t664) * t469 + t754; m(7) * (t467 * t98 - t470 * t97 + (t178 * t470 + t179 * t467) * qJD(1)) + m(6) * (-t116 * t470 + t117 * t467 + (t215 * t470 + t216 * t467) * qJD(1)) + m(5) * (-t145 * t470 + t146 * t467 + (t229 * t470 + t230 * t467) * qJD(1)) + m(4) * (-t233 * t470 + t234 * t467 + (t341 * t470 + t342 * t467) * qJD(1)) + m(3) * (-t354 * t470 + t355 * t467 + (t395 * t470 + t396 * t467) * qJD(1)); 0; m(5) * (t145 * t337 + t146 * t336 + t208 * t230 + t209 * t229) + m(6) * (t116 * t232 + t117 * t231 + t131 * t216 + t132 * t215) + m(7) * (t112 * t179 + t113 * t178 + t217 * t98 + t218 * t97) + ((qJD(1) * t722 + t467 * t487) * t715 + (qJD(1) * t721 + t467 * t488) * t711 + t82 / 0.2e1 + t83 / 0.2e1 + t48 / 0.2e1 + t46 / 0.2e1 + m(4) * (-t233 * t426 + t342 * t417) - t413 * t710 + (-t678 / 0.2e1 - t674 / 0.2e1) * qJD(3) - t605) * t470 + ((qJD(1) * t382 - t533 * t617) * t715 + (qJD(1) * t385 - t538 * t617) * t711 + t80 / 0.2e1 + t81 / 0.2e1 + t49 / 0.2e1 + t47 / 0.2e1 + m(4) * (t234 * t426 - t341 * t417) - t413 * t712 + (t676 / 0.2e1 + t672 / 0.2e1) * qJD(3) + t604) * t467 + ((t342 * t706 + t679 / 0.2e1 - t673 / 0.2e1 - t189 / 0.2e1 - t190 / 0.2e1 - t143 / 0.2e1 - t141 / 0.2e1 - t574) * t467 + (t341 * t706 + t191 / 0.2e1 + t192 / 0.2e1 + t677 / 0.2e1 - t671 / 0.2e1 + t144 / 0.2e1 + t142 / 0.2e1 - t573) * t470) * qJD(1); m(5) * (-t208 * t470 + t209 * t467 + (t336 * t470 + t337 * t467) * qJD(1)) + m(6) * (-t131 * t470 + t132 * t467 + (t231 * t470 + t232 * t467) * qJD(1)) + m(7) * (-t112 * t470 + t113 * t467 + (t217 * t470 + t218 * t467) * qJD(1)) - m(4) * t565; (t107 * t28 + t112 * t218 + t113 * t217) * t717 + (t118 * t33 + t131 * t232 + t132 * t231) * t718 + t470 * t21 + t467 * t20 + t467 * t19 + t470 * t22 + (t180 * t79 + t208 * t337 + t209 * t336) * t719 + t470 * t26 + t467 * t25 + ((-t387 * t467 + t389 * t470) * (-t467 * t593 + (-t426 * t463 + t462 * t700) * qJD(3) + ((-t387 + t458) * t470 + (-t389 + t492 + t697) * t467) * qJD(1)) - t426 * t565) * t720 + t470 * ((t470 * t326 + (t223 + t726) * qJD(1)) * t470 + (-t222 * qJD(1) + (-t618 * t721 + t620 * t722) * t467 + (t325 + (t673 - t679) * qJD(3) + (-t379 + t502) * qJD(1)) * t470) * t467) + t467 * ((t467 * t325 + (-t225 + t489) * qJD(1)) * t467 + (t226 * qJD(1) + (t382 * t620 - t385 * t618 + t625) * t470 + (t326 + (t671 - t677) * qJD(3) + t503 * qJD(1)) * t467) * t470) + (-t222 * t470 - t223 * t467 - t105 - t739) * t623 + (t225 * t470 + t226 * t467 + t106 + t738) * t621; m(5) * (t114 * t229 + t115 * t230 + t145 * t227 + t146 * t228) + m(6) * (t116 * t151 + t117 * t152 + t215 * t58 + t216 * t59) + m(7) * (t119 * t97 + t120 * t98 + t178 * t31 + t179 * t32) + (t604 * t470 + t605 * t467 + (t467 * t573 - t470 * t574) * qJD(1)) * t469 + t473 + (t467 * t574 + t470 * t573) * t620 + t651; m(5) * (t114 * t467 - t115 * t470 + (t227 * t467 + t228 * t470) * qJD(1)) + m(6) * (t467 * t58 - t470 * t59 + (t151 * t467 + t152 * t470) * qJD(1)) + m(7) * (t31 * t467 - t32 * t470 + (t119 * t467 + t120 * t470) * qJD(1)); m(5) * (t114 * t336 + t115 * t337 + t180 * t84 + t208 * t227 + t209 * t228 - t210 * t79) + m(6) * (t118 * t30 + t122 * t33 + t131 * t151 + t132 * t152 + t231 * t58 + t232 * t59) + m(7) * (t107 * t27 + t110 * t28 + t112 * t119 + t113 * t120 + t217 * t31 + t218 * t32) + (t105 * t568 + (-qJD(1) * t160 + t57) * t714 + t13 / 0.2e1 - qJD(1) * t73 / 0.2e1) * t467 + t475 + (t106 * t569 + (qJD(1) * t161 + t56) * t714 + t14 / 0.2e1 + qJD(1) * t74 / 0.2e1) * t470 + (qJD(3) * (t160 * t470 + t161 * t467) / 0.2e1 + t25 * t710 + t26 * t713 + (-t470 * t105 / 0.2e1 + t106 * t713) * qJD(1)) * t469; (t114 * t228 + t115 * t227 - t210 * t84) * t719 + (t122 * t30 + t151 * t59 + t152 * t58) * t718 + (t110 * t27 + t119 * t32 + t120 * t31) * t717 + ((t467 * t566 + t470 * t501) * qJD(3) + t651) * t466 + ((t466 * t57 + t13) * t470 + (-t466 * t56 - t14 - t744) * t467 + (t224 * t466 + (-t160 * t467 + t161 * t470) * t469) * qJD(3) + ((-t566 - t741) * t470 + t501 * t467) * qJD(1)) * t469 + t490; t473 + m(6) * (t108 * t215 + t109 * t216 + t116 * t219 + t117 * t220) + m(7) * (t153 * t97 + t154 * t98 + t178 * t54 + t179 * t55); m(6) * (t108 * t467 - t109 * t470 + (t219 * t467 + t220 * t470) * qJD(1)) + m(7) * (t467 * t54 - t470 * t55 + (t153 * t467 + t154 * t470) * qJD(1)); t475 + m(6) * (t108 * t231 + t109 * t232 + t118 * t70 + t131 * t219 + t132 * t220 + t201 * t33) + m(7) * (t107 * t29 + t112 * t153 + t113 * t154 + t121 * t28 + t217 * t54 + t218 * t55); m(7) * (t110 * t29 + t119 * t55 + t120 * t54 + t121 * t27 + t153 * t32 + t154 * t31) + m(6) * (t108 * t152 + t109 * t151 + t122 * t70 + t201 * t30 + t219 * t59 + t220 * t58) + t472; (t121 * t29 + t153 * t55 + t154 * t54) * t717 + (t108 * t220 + t109 * t219 + t201 * t70) * t718 + t472; m(7) * (t178 * t247 + t179 * t245 + t390 * t98 - t392 * t97); m(7) * (-t245 * t470 + t247 * t467 + (t390 * t470 - t392 * t467) * qJD(1)); m(7) * (t107 * t486 - t112 * t392 + t113 * t390 + t217 * t247 + t218 * t245 + t28 * t669); m(7) * (t110 * t486 + t119 * t245 + t120 * t247 + t27 * t669 + t31 * t390 - t32 * t392); m(7) * (t121 * t486 + t153 * t245 + t154 * t247 + t29 * t669 + t390 * t54 - t392 * t55); (-t245 * t392 + t247 * t390 + t486 * t669) * t717;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
