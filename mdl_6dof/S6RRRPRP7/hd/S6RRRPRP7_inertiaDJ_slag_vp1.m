% Calculate time derivative of joint inertia matrix for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:04:09
% EndTime: 2019-03-09 17:05:19
% DurationCPUTime: 43.08s
% Computational Cost: add. (130441->1548), mult. (252068->2042), div. (0->0), fcn. (281125->12), ass. (0->592)
t697 = qJ(3) + pkin(11);
t564 = sin(t697);
t567 = cos(pkin(6));
t571 = sin(qJ(2));
t566 = sin(pkin(6));
t636 = cos(t697);
t613 = t566 * t636;
t519 = t567 * t564 + t571 * t613;
t754 = cos(qJ(2));
t647 = qJD(2) * t754;
t620 = t566 * t647;
t484 = qJD(3) * t519 + t564 * t620;
t741 = t566 * t571;
t518 = t564 * t741 - t567 * t636;
t485 = -qJD(3) * t518 + t613 * t647;
t701 = qJD(2) * t571;
t648 = t566 * t701;
t385 = Icges(5,5) * t485 - Icges(5,6) * t484 + Icges(5,3) * t648;
t386 = Icges(5,4) * t485 - Icges(5,2) * t484 + Icges(5,6) * t648;
t387 = Icges(5,1) * t485 - Icges(5,4) * t484 + Icges(5,5) * t648;
t656 = t566 * t754;
t445 = Icges(5,5) * t519 - Icges(5,6) * t518 - Icges(5,3) * t656;
t446 = Icges(5,4) * t519 - Icges(5,2) * t518 - Icges(5,6) * t656;
t447 = Icges(5,1) * t519 - Icges(5,4) * t518 - Icges(5,5) * t656;
t147 = -t385 * t656 - t518 * t386 + t519 * t387 + t445 * t648 - t484 * t446 + t485 * t447;
t570 = sin(qJ(3));
t737 = t567 * t570;
t573 = cos(qJ(3));
t739 = t566 * t573;
t537 = t571 * t739 + t737;
t502 = -qJD(3) * t537 - t570 * t620;
t536 = t567 * t573 - t570 * t741;
t503 = qJD(3) * t536 + t573 * t620;
t423 = Icges(4,5) * t503 + Icges(4,6) * t502 + Icges(4,3) * t648;
t424 = Icges(4,4) * t503 + Icges(4,2) * t502 + Icges(4,6) * t648;
t425 = Icges(4,1) * t503 + Icges(4,4) * t502 + Icges(4,5) * t648;
t453 = Icges(4,5) * t537 + Icges(4,6) * t536 - Icges(4,3) * t656;
t454 = Icges(4,4) * t537 + Icges(4,2) * t536 - Icges(4,6) * t656;
t455 = Icges(4,1) * t537 + Icges(4,4) * t536 - Icges(4,5) * t656;
t176 = -t423 * t656 + t536 * t424 + t537 * t425 + t453 * t648 + t502 * t454 + t503 * t455;
t771 = t147 + t176;
t574 = cos(qJ(1));
t654 = t574 * t754;
t572 = sin(qJ(1));
t736 = t572 * t571;
t541 = -t567 * t736 + t654;
t623 = t567 * t654;
t593 = t623 - t736;
t483 = qJD(1) * t541 + qJD(2) * t593;
t655 = t572 * t754;
t735 = t574 * t571;
t594 = -t567 * t735 - t655;
t585 = t564 * t594 - t574 * t613;
t704 = qJD(1) * t572;
t650 = t566 * t704;
t359 = qJD(3) * t585 + t483 * t636 + t564 * t650;
t738 = t566 * t574;
t497 = -t564 * t738 - t594 * t636;
t569 = sin(qJ(5));
t753 = cos(qJ(5));
t440 = t497 * t753 - t569 * t593;
t595 = -t567 * t655 - t735;
t482 = -qJD(1) * t595 - qJD(2) * t594;
t276 = qJD(5) * t440 + t359 * t569 - t482 * t753;
t601 = -t497 * t569 - t593 * t753;
t277 = qJD(5) * t601 + t359 * t753 + t482 * t569;
t603 = qJD(1) * t613;
t358 = qJD(3) * t497 + t483 * t564 - t572 * t603;
t769 = rSges(7,3) + qJ(6);
t770 = rSges(7,1) + pkin(5);
t733 = t358 * rSges(7,2) - t601 * qJD(6) + t769 * t276 + t277 * t770;
t726 = -rSges(7,2) * t585 + t440 * t770 - t769 * t601;
t363 = Icges(5,5) * t497 + Icges(5,6) * t585 - Icges(5,3) * t593;
t365 = Icges(5,4) * t497 + Icges(5,2) * t585 - Icges(5,6) * t593;
t367 = Icges(5,1) * t497 + Icges(5,4) * t585 - Icges(5,5) * t593;
t205 = -t363 * t656 - t518 * t365 + t519 * t367;
t504 = t570 * t594 - t573 * t738;
t673 = t570 * t738;
t602 = t573 * t594 + t673;
t402 = -Icges(4,5) * t602 + Icges(4,6) * t504 - Icges(4,3) * t593;
t404 = -Icges(4,4) * t602 + Icges(4,2) * t504 - Icges(4,6) * t593;
t406 = -Icges(4,1) * t602 + Icges(4,4) * t504 - Icges(4,5) * t593;
t217 = -t402 * t656 + t536 * t404 + t537 * t406;
t763 = t205 + t217;
t740 = t566 * t572;
t499 = t541 * t636 + t564 * t740;
t586 = -t541 * t564 + t572 * t613;
t364 = Icges(5,5) * t499 + Icges(5,6) * t586 - Icges(5,3) * t595;
t366 = Icges(5,4) * t499 + Icges(5,2) * t586 - Icges(5,6) * t595;
t368 = Icges(5,1) * t499 + Icges(5,4) * t586 - Icges(5,5) * t595;
t206 = -t364 * t656 - t518 * t366 + t519 * t368;
t506 = -t541 * t570 + t572 * t739;
t674 = t570 * t740;
t507 = t541 * t573 + t674;
t403 = Icges(4,5) * t507 + Icges(4,6) * t506 - Icges(4,3) * t595;
t405 = Icges(4,4) * t507 + Icges(4,2) * t506 - Icges(4,6) * t595;
t407 = Icges(4,1) * t507 + Icges(4,4) * t506 - Icges(4,5) * t595;
t218 = -t403 * t656 + t536 * t405 + t537 * t407;
t762 = t206 + t218;
t768 = t771 * t567;
t295 = Icges(7,5) * t440 - Icges(7,6) * t585 - Icges(7,3) * t601;
t299 = Icges(7,4) * t440 - Icges(7,2) * t585 - Icges(7,6) * t601;
t303 = Icges(7,1) * t440 - Icges(7,4) * t585 - Icges(7,5) * t601;
t495 = t519 * t753 - t569 * t656;
t589 = -t519 * t569 - t656 * t753;
t148 = -t295 * t589 + t299 * t518 + t303 * t495;
t442 = t499 * t753 - t569 * t595;
t600 = -t499 * t569 - t595 * t753;
t296 = Icges(7,5) * t442 - Icges(7,6) * t586 - Icges(7,3) * t600;
t300 = Icges(7,4) * t442 - Icges(7,2) * t586 - Icges(7,6) * t600;
t304 = Icges(7,1) * t442 - Icges(7,4) * t586 - Icges(7,5) * t600;
t149 = -t296 * t589 + t300 * t518 + t304 * t495;
t481 = qJD(1) * t594 + qJD(2) * t595;
t356 = qJD(3) * t499 + t481 * t564 - t574 * t603;
t155 = Icges(7,5) * t277 + Icges(7,6) * t358 + Icges(7,3) * t276;
t159 = Icges(7,4) * t277 + Icges(7,2) * t358 + Icges(7,6) * t276;
t163 = Icges(7,1) * t277 + Icges(7,4) * t358 + Icges(7,5) * t276;
t376 = qJD(5) * t495 + t485 * t569 - t753 * t648;
t377 = qJD(5) * t589 + t485 * t753 + t569 * t648;
t48 = -t155 * t589 + t159 * t518 + t163 * t495 + t295 * t376 + t299 * t484 + t303 * t377;
t703 = qJD(1) * t574;
t649 = t566 * t703;
t357 = qJD(3) * t586 + t481 * t636 + t564 * t649;
t480 = -qJD(1) * t623 - t574 * t647 + (qJD(2) * t567 + qJD(1)) * t736;
t274 = qJD(5) * t442 + t357 * t569 + t480 * t753;
t275 = qJD(5) * t600 + t357 * t753 - t480 * t569;
t154 = Icges(7,5) * t275 + Icges(7,6) * t356 + Icges(7,3) * t274;
t158 = Icges(7,4) * t275 + Icges(7,2) * t356 + Icges(7,6) * t274;
t162 = Icges(7,1) * t275 + Icges(7,4) * t356 + Icges(7,5) * t274;
t49 = -t154 * t589 + t158 * t518 + t162 * t495 + t296 * t376 + t300 * t484 + t304 * t377;
t336 = Icges(7,5) * t495 + Icges(7,6) * t518 - Icges(7,3) * t589;
t338 = Icges(7,4) * t495 + Icges(7,2) * t518 - Icges(7,6) * t589;
t340 = Icges(7,1) * t495 + Icges(7,4) * t518 - Icges(7,5) * t589;
t194 = -t336 * t589 + t338 * t518 + t340 * t495;
t246 = Icges(7,5) * t377 + Icges(7,6) * t484 + Icges(7,3) * t376;
t248 = Icges(7,4) * t377 + Icges(7,2) * t484 + Icges(7,6) * t376;
t250 = Icges(7,1) * t377 + Icges(7,4) * t484 + Icges(7,5) * t376;
t93 = -t246 * t589 + t518 * t248 + t495 * t250 + t376 * t336 + t484 * t338 + t377 * t340;
t746 = t194 * t484 + t93 * t518;
t13 = t148 * t358 + t149 * t356 - t48 * t585 - t49 * t586 + t746;
t297 = Icges(6,5) * t440 + Icges(6,6) * t601 - Icges(6,3) * t585;
t301 = Icges(6,4) * t440 + Icges(6,2) * t601 - Icges(6,6) * t585;
t305 = Icges(6,1) * t440 + Icges(6,4) * t601 - Icges(6,5) * t585;
t150 = t297 * t518 + t301 * t589 + t305 * t495;
t298 = Icges(6,5) * t442 + Icges(6,6) * t600 - Icges(6,3) * t586;
t302 = Icges(6,4) * t442 + Icges(6,2) * t600 - Icges(6,6) * t586;
t306 = Icges(6,1) * t442 + Icges(6,4) * t600 - Icges(6,5) * t586;
t151 = t298 * t518 + t302 * t589 + t306 * t495;
t157 = Icges(6,5) * t277 - Icges(6,6) * t276 + Icges(6,3) * t358;
t161 = Icges(6,4) * t277 - Icges(6,2) * t276 + Icges(6,6) * t358;
t165 = Icges(6,1) * t277 - Icges(6,4) * t276 + Icges(6,5) * t358;
t50 = t157 * t518 + t161 * t589 + t165 * t495 + t297 * t484 - t301 * t376 + t305 * t377;
t156 = Icges(6,5) * t275 - Icges(6,6) * t274 + Icges(6,3) * t356;
t160 = Icges(6,4) * t275 - Icges(6,2) * t274 + Icges(6,6) * t356;
t164 = Icges(6,1) * t275 - Icges(6,4) * t274 + Icges(6,5) * t356;
t51 = t156 * t518 + t160 * t589 + t164 * t495 + t298 * t484 - t302 * t376 + t306 * t377;
t337 = Icges(6,5) * t495 + Icges(6,6) * t589 + Icges(6,3) * t518;
t339 = Icges(6,4) * t495 + Icges(6,2) * t589 + Icges(6,6) * t518;
t341 = Icges(6,1) * t495 + Icges(6,4) * t589 + Icges(6,5) * t518;
t195 = t337 * t518 + t339 * t589 + t341 * t495;
t247 = Icges(6,5) * t377 - Icges(6,6) * t376 + Icges(6,3) * t484;
t249 = Icges(6,4) * t377 - Icges(6,2) * t376 + Icges(6,6) * t484;
t251 = Icges(6,1) * t377 - Icges(6,4) * t376 + Icges(6,5) * t484;
t94 = t518 * t247 + t249 * t589 + t495 * t251 + t484 * t337 - t376 * t339 + t377 * t341;
t745 = t195 * t484 + t94 * t518;
t14 = t150 * t358 + t151 * t356 - t50 * t585 - t51 * t586 + t745;
t764 = t14 + t13;
t734 = t356 * rSges(7,2) - qJD(6) * t600 + t769 * t274 + t275 * t770;
t725 = -rSges(7,2) * t586 + t442 * t770 - t769 * t600;
t565 = t574 * pkin(1);
t705 = pkin(8) * t740 + t565;
t236 = Icges(5,5) * t359 - Icges(5,6) * t358 + Icges(5,3) * t482;
t238 = Icges(5,4) * t359 - Icges(5,2) * t358 + Icges(5,6) * t482;
t240 = Icges(5,1) * t359 - Icges(5,4) * t358 + Icges(5,5) * t482;
t101 = -t518 * t238 + t519 * t240 - t484 * t365 + t485 * t367 + (-t236 * t754 + t363 * t701) * t566;
t235 = Icges(5,5) * t357 - Icges(5,6) * t356 - Icges(5,3) * t480;
t237 = Icges(5,4) * t357 - Icges(5,2) * t356 - Icges(5,6) * t480;
t239 = Icges(5,1) * t357 - Icges(5,4) * t356 - Icges(5,5) * t480;
t102 = -t518 * t237 + t519 * t239 - t484 * t366 + t485 * t368 + (-t235 * t754 + t364 * t701) * t566;
t398 = qJD(3) * t602 - t483 * t570 + t573 * t650;
t621 = t570 * t650;
t399 = qJD(3) * t504 + t483 * t573 + t621;
t259 = Icges(4,5) * t399 + Icges(4,6) * t398 + Icges(4,3) * t482;
t261 = Icges(4,4) * t399 + Icges(4,2) * t398 + Icges(4,6) * t482;
t263 = Icges(4,1) * t399 + Icges(4,4) * t398 + Icges(4,5) * t482;
t111 = t536 * t261 + t537 * t263 + t502 * t404 + t503 * t406 + (-t259 * t754 + t402 * t701) * t566;
t396 = -qJD(3) * t507 - t481 * t570 + t573 * t649;
t622 = t570 * t649;
t397 = qJD(3) * t506 + t481 * t573 + t622;
t258 = Icges(4,5) * t397 + Icges(4,6) * t396 - Icges(4,3) * t480;
t260 = Icges(4,4) * t397 + Icges(4,2) * t396 - Icges(4,6) * t480;
t262 = Icges(4,1) * t397 + Icges(4,4) * t396 - Icges(4,5) * t480;
t112 = t536 * t260 + t537 * t262 + t502 * t405 + t503 * t407 + (-t258 * t754 + t403 * t701) * t566;
t91 = t93 * t567;
t17 = t91 + (-t48 * t574 + t49 * t572 + (t148 * t572 + t149 * t574) * qJD(1)) * t566;
t92 = t94 * t567;
t18 = t92 + (-t50 * t574 + t51 * t572 + (t150 * t572 + t151 * t574) * qJD(1)) * t566;
t761 = t17 + t18 + t768 + ((-t101 - t111) * t574 + (t102 + t112) * t572 + (t572 * t763 + t574 * t762) * qJD(1)) * t566;
t760 = -t93 - t94 - t771;
t759 = t566 ^ 2;
t758 = m(7) / 0.2e1;
t563 = pkin(3) * t573 + pkin(2);
t752 = -pkin(2) + t563;
t750 = pkin(3) * qJD(3);
t568 = -qJ(4) - pkin(9);
t748 = -rSges(5,3) + t568;
t744 = Icges(3,4) * t571;
t743 = t482 * t568;
t742 = t593 * t568;
t479 = t482 * pkin(9);
t684 = t573 * t750;
t634 = t566 * t684;
t685 = t570 * t750;
t658 = -qJD(4) * t593 - t574 * t634 + t594 * t685;
t244 = pkin(3) * t621 + t483 * t752 - t479 + t658 - t743;
t534 = t593 * pkin(9);
t554 = pkin(3) * t673;
t374 = -t594 * t752 + t534 - t554 + t742;
t732 = -t244 * t595 - t480 * t374;
t415 = t481 * pkin(2) - pkin(9) * t480;
t584 = pkin(3) * t622 - qJD(4) * t595 + t480 * t568 + t481 * t563 - t541 * t685 + t572 * t634;
t243 = -t415 + t584;
t401 = t567 * t415;
t731 = t567 * t243 + t401;
t283 = t357 * pkin(4) + pkin(10) * t356;
t730 = -t243 - t283;
t416 = t483 * pkin(2) + t479;
t729 = -t244 - t416;
t728 = rSges(7,2) * t484 - qJD(6) * t589 + t769 * t376 + t377 * t770;
t265 = t399 * rSges(4,1) + t398 * rSges(4,2) + t482 * rSges(4,3);
t727 = -t265 - t416;
t331 = t595 * t374;
t428 = pkin(4) * t497 - pkin(10) * t585;
t724 = -t428 * t595 - t331;
t487 = t541 * pkin(2) - pkin(9) * t595;
t659 = pkin(3) * t674 + t541 * t563 + t568 * t595;
t375 = -t487 + t659;
t342 = t375 * t648;
t429 = t499 * pkin(4) - pkin(10) * t586;
t723 = t429 * t648 + t342;
t722 = rSges(7,2) * t518 + t495 * t770 - t769 * t589;
t686 = t754 * pkin(9);
t456 = pkin(3) * t737 + (t568 * t754 + t571 * t752 + t686) * t566;
t721 = t374 * t656 - t456 * t593;
t473 = t567 * t487;
t720 = t567 * t375 + t473;
t615 = -rSges(5,1) * t497 - rSges(5,2) * t585;
t371 = -rSges(5,3) * t593 - t615;
t719 = -t371 - t374;
t372 = t499 * rSges(5,1) + rSges(5,2) * t586 - rSges(5,3) * t595;
t718 = -t372 - t375;
t717 = -t374 - t428;
t716 = -t375 - t429;
t390 = rSges(5,1) * t485 - rSges(5,2) * t484 + rSges(5,3) * t648;
t687 = t754 * pkin(2);
t430 = t567 * t684 + (-t571 * t685 - t754 * qJD(4) + (-t687 + t754 * t563 + (-pkin(9) - t568) * t571) * qJD(2)) * t566;
t715 = -t390 - t430;
t408 = -rSges(4,1) * t602 + rSges(4,2) * t504 - rSges(4,3) * t593;
t486 = -pkin(2) * t594 - t534;
t714 = -t408 - t486;
t409 = t507 * rSges(4,1) + t506 * rSges(4,2) - rSges(4,3) * t595;
t713 = -t409 - t487;
t417 = pkin(4) * t485 + pkin(10) * t484;
t712 = -t417 - t430;
t426 = rSges(4,1) * t503 + rSges(4,2) * t502 + rSges(4,3) * t648;
t702 = qJD(2) * t566;
t528 = (pkin(9) * t571 + t687) * t702;
t711 = -t426 - t528;
t448 = t519 * rSges(5,1) - t518 * rSges(5,2) - rSges(5,3) * t656;
t710 = t448 + t456;
t542 = (pkin(2) * t571 - t686) * t566;
t511 = t542 * t650;
t709 = t456 * t650 + t511;
t452 = pkin(4) * t519 + pkin(10) * t518;
t708 = t452 + t456;
t457 = t537 * rSges(4,1) + t536 * rSges(4,2) - rSges(4,3) * t656;
t707 = -t457 - t542;
t706 = t486 * t740 + t487 * t738;
t138 = -t295 * t600 - t299 * t586 + t303 * t442;
t139 = -t296 * t600 - t300 * t586 + t304 * t442;
t181 = -t336 * t600 - t338 * t586 + t340 * t442;
t38 = -t155 * t600 - t159 * t586 + t163 * t442 + t274 * t295 + t275 * t303 + t299 * t356;
t39 = -t154 * t600 - t158 * t586 + t162 * t442 + t274 * t296 + t275 * t304 + t300 * t356;
t72 = -t246 * t600 - t248 * t586 + t250 * t442 + t274 * t336 + t275 * t340 + t338 * t356;
t1 = t138 * t358 + t139 * t356 + t181 * t484 - t38 * t585 - t39 * t586 + t518 * t72;
t140 = -t297 * t586 + t301 * t600 + t305 * t442;
t141 = -t298 * t586 + t302 * t600 + t306 * t442;
t182 = -t337 * t586 + t339 * t600 + t341 * t442;
t40 = -t157 * t586 + t161 * t600 + t165 * t442 - t274 * t301 + t275 * t305 + t297 * t356;
t41 = -t156 * t586 + t160 * t600 + t164 * t442 - t274 * t302 + t275 * t306 + t298 * t356;
t73 = -t247 * t586 + t249 * t600 + t251 * t442 - t274 * t339 + t275 * t341 + t337 * t356;
t2 = t140 * t358 + t141 * t356 + t182 * t484 - t40 * t585 - t41 * t586 + t518 * t73;
t696 = t1 / 0.2e1 + t2 / 0.2e1;
t134 = -t295 * t601 - t299 * t585 + t303 * t440;
t135 = -t296 * t601 - t300 * t585 + t304 * t440;
t179 = -t336 * t601 - t338 * t585 + t340 * t440;
t42 = -t155 * t601 - t159 * t585 + t163 * t440 + t276 * t295 + t277 * t303 + t299 * t358;
t43 = -t154 * t601 - t158 * t585 + t162 * t440 + t276 * t296 + t277 * t304 + t300 * t358;
t74 = -t246 * t601 - t248 * t585 + t250 * t440 + t276 * t336 + t277 * t340 + t338 * t358;
t3 = t134 * t358 + t135 * t356 + t179 * t484 - t42 * t585 - t43 * t586 + t518 * t74;
t136 = -t297 * t585 + t301 * t601 + t305 * t440;
t137 = -t298 * t585 + t302 * t601 + t306 * t440;
t180 = -t337 * t585 + t339 * t601 + t341 * t440;
t44 = -t157 * t585 + t161 * t601 + t165 * t440 - t276 * t301 + t277 * t305 + t297 * t358;
t45 = -t156 * t585 + t160 * t601 + t164 * t440 - t276 * t302 + t277 * t306 + t298 * t358;
t75 = -t247 * t585 + t249 * t601 + t251 * t440 - t276 * t339 + t277 * t341 + t337 * t358;
t4 = t136 * t358 + t137 * t356 + t180 * t484 - t44 * t585 - t45 * t586 + t518 * t75;
t695 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t138 * t482 - t139 * t480 - t38 * t593 - t39 * t595 + (t181 * t701 - t72 * t754) * t566;
t6 = t140 * t482 - t141 * t480 - t40 * t593 - t41 * t595 + (t182 * t701 - t73 * t754) * t566;
t694 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t134 * t482 - t135 * t480 - t42 * t593 - t43 * t595 + (t179 * t701 - t74 * t754) * t566;
t8 = t136 * t482 - t137 * t480 - t44 * t593 - t45 * t595 + (t180 * t701 - t75 * t754) * t566;
t693 = t8 / 0.2e1 + t7 / 0.2e1;
t692 = -t754 / 0.2e1;
t10 = t73 * t567 + (-t40 * t574 + t41 * t572 + (t140 * t572 + t141 * t574) * qJD(1)) * t566;
t9 = t72 * t567 + (-t38 * t574 + t39 * t572 + (t138 * t572 + t139 * t574) * qJD(1)) * t566;
t691 = t9 / 0.2e1 + t10 / 0.2e1;
t690 = m(5) * t754;
t689 = m(6) * t754;
t688 = m(7) * t754;
t11 = t74 * t567 + (-t42 * t574 + t43 * t572 + (t134 * t572 + t135 * t574) * qJD(1)) * t566;
t12 = t75 * t567 + (-t44 * t574 + t45 * t572 + (t136 * t572 + t137 * t574) * qJD(1)) * t566;
t683 = t11 / 0.2e1 + t12 / 0.2e1;
t190 = t194 * t648;
t15 = t148 * t482 - t149 * t480 - t48 * t593 - t49 * t595 - t656 * t93 + t190;
t191 = t195 * t648;
t16 = t150 * t482 - t151 * t480 - t50 * t593 - t51 * t595 - t656 * t94 + t191;
t682 = t16 / 0.2e1 + t15 / 0.2e1;
t56 = -t134 * t585 - t135 * t586 + t179 * t518;
t57 = -t136 * t585 - t137 * t586 + t180 * t518;
t681 = t56 / 0.2e1 + t57 / 0.2e1;
t58 = -t138 * t585 - t139 * t586 + t181 * t518;
t59 = -t140 * t585 - t141 * t586 + t182 * t518;
t680 = t58 / 0.2e1 + t59 / 0.2e1;
t60 = -t134 * t593 - t135 * t595 - t179 * t656;
t61 = -t136 * t593 - t137 * t595 - t180 * t656;
t679 = t60 / 0.2e1 + t61 / 0.2e1;
t62 = -t138 * t593 - t139 * t595 - t181 * t656;
t63 = -t140 * t593 - t141 * t595 - t182 * t656;
t678 = t62 / 0.2e1 + t63 / 0.2e1;
t64 = t179 * t567 + (-t134 * t574 + t135 * t572) * t566;
t65 = t180 * t567 + (-t136 * t574 + t137 * t572) * t566;
t677 = t64 / 0.2e1 + t65 / 0.2e1;
t66 = t181 * t567 + (-t138 * t574 + t139 * t572) * t566;
t67 = t182 * t567 + (-t140 * t574 + t141 * t572) * t566;
t676 = t66 / 0.2e1 + t67 / 0.2e1;
t675 = ((-t148 - t150) * t574 + (t149 + t151) * t572) * t566 / 0.2e1 + (t194 + t195) * t567 / 0.2e1;
t671 = t244 * t656 - t430 * t593 + t482 * t456;
t670 = t567 * t283 + t731;
t284 = t359 * pkin(4) + t358 * pkin(10);
t669 = -t284 + t729;
t253 = rSges(6,1) * t377 - rSges(6,2) * t376 + rSges(6,3) * t484;
t668 = -t253 + t712;
t167 = t275 * rSges(6,1) - t274 * rSges(6,2) + t356 * rSges(6,3);
t308 = rSges(6,1) * t440 + rSges(6,2) * t601 - rSges(6,3) * t585;
t667 = -t308 + t717;
t310 = t442 * rSges(6,1) + rSges(6,2) * t600 - rSges(6,3) * t586;
t666 = -t310 + t716;
t345 = rSges(6,1) * t495 + rSges(6,2) * t589 + rSges(6,3) * t518;
t665 = t345 + t708;
t241 = t357 * rSges(5,1) - t356 * rSges(5,2) - t480 * rSges(5,3);
t664 = t567 * t429 + t720;
t663 = -t528 + t715;
t264 = t397 * rSges(4,1) + t396 * rSges(4,2) - t480 * rSges(4,3);
t662 = t415 * t738 + t416 * t740 + t486 * t649;
t661 = -t542 - t710;
t660 = t452 * t650 + t709;
t388 = t481 * rSges(3,1) + t480 * rSges(3,2) + rSges(3,3) * t649;
t465 = t541 * rSges(3,1) + rSges(3,2) * t595 + rSges(3,3) * t740;
t653 = t754 * Icges(3,4);
t652 = t754 * t243;
t651 = t754 * t375;
t646 = -t572 * pkin(1) + pkin(8) * t738;
t645 = 2 * m(3);
t643 = 2 * m(4);
t641 = 0.2e1 * m(5);
t639 = 0.2e1 * m(6);
t637 = 0.2e1 * m(7);
t635 = t707 * t574;
t633 = -t284 * t595 - t480 * t428 + t732;
t632 = t712 - t728;
t631 = -t528 + t668;
t630 = t717 - t726;
t629 = t716 - t725;
t628 = t708 + t722;
t627 = -t542 - t665;
t626 = t374 * t740 + t375 * t738 + t706;
t625 = t428 * t656 - t452 * t593 + t721;
t619 = -pkin(1) * t704 + pkin(8) * t649;
t618 = t661 * t574;
t617 = -t483 * rSges(3,1) + t482 * rSges(3,2);
t616 = -t359 * rSges(5,1) + t358 * rSges(5,2);
t612 = -t528 + t632;
t611 = -t542 - t628;
t610 = t659 + t705;
t609 = t627 * t574;
t608 = t50 / 0.2e1 + t48 / 0.2e1 + t75 / 0.2e1 + t74 / 0.2e1;
t607 = t72 / 0.2e1 + t73 / 0.2e1 + t51 / 0.2e1 + t49 / 0.2e1;
t606 = t243 * t738 + t244 * t740 + t374 * t649 + t662;
t605 = t284 * t656 - t417 * t593 + t482 * t452 + t671;
t604 = t428 * t740 + t429 * t738 + t626;
t599 = t150 / 0.2e1 + t148 / 0.2e1 + t180 / 0.2e1 + t179 / 0.2e1;
t598 = t181 / 0.2e1 + t182 / 0.2e1 + t151 / 0.2e1 + t149 / 0.2e1;
t597 = t611 * t574;
t596 = t563 * t594 + t554 + t646;
t169 = t277 * rSges(6,1) - t276 * rSges(6,2) + t358 * rSges(6,3);
t592 = -t283 * t754 - t652;
t591 = -t429 * t754 - t651;
t464 = -rSges(3,1) * t594 + rSges(3,2) * t593 - rSges(3,3) * t738;
t590 = t429 + t610;
t516 = Icges(3,6) * t567 + (Icges(3,2) * t754 + t744) * t566;
t517 = Icges(3,5) * t567 + (Icges(3,1) * t571 + t653) * t566;
t524 = (Icges(3,5) * t754 - Icges(3,6) * t571) * t702;
t525 = (-Icges(3,2) * t571 + t653) * t702;
t526 = (Icges(3,1) * t754 - t744) * t702;
t588 = -t516 * t648 + t517 * t620 + t567 * t524 + t525 * t656 + t526 * t741;
t587 = t283 * t738 + t284 * t740 + t428 * t649 + t606;
t127 = -t356 * t446 + t357 * t447 - t385 * t595 + t386 * t586 + t387 * t499 - t445 * t480;
t132 = t396 * t454 + t397 * t455 - t423 * t595 + t424 * t506 + t425 * t507 - t453 * t480;
t583 = t127 / 0.2e1 + t132 / 0.2e1 + t112 / 0.2e1 + t102 / 0.2e1 + t607;
t128 = -t358 * t446 + t359 * t447 - t385 * t593 + t386 * t585 + t387 * t497 + t445 * t482;
t133 = t398 * t454 + t399 * t455 - t423 * t593 + t424 * t504 - t425 * t602 + t453 * t482;
t582 = t128 / 0.2e1 + t133 / 0.2e1 + t111 / 0.2e1 + t101 / 0.2e1 + t608;
t224 = -t445 * t593 + t446 * t585 + t447 * t497;
t256 = -t453 * t593 + t454 * t504 - t455 * t602;
t581 = t256 / 0.2e1 + t224 / 0.2e1 + t217 / 0.2e1 + t205 / 0.2e1 + t599;
t225 = -t445 * t595 + t446 * t586 + t447 * t499;
t257 = -t453 * t595 + t454 * t506 + t455 * t507;
t580 = -t257 / 0.2e1 - t225 / 0.2e1 - t218 / 0.2e1 - t206 / 0.2e1 - t598;
t579 = -t428 + t596 - t742;
t578 = t584 + t619;
t577 = (-t565 + (-pkin(3) * t570 - pkin(8)) * t740) * qJD(1) - t483 * t563 - t658;
t576 = t283 + t578;
t575 = -t284 + t577 + t743;
t527 = (rSges(3,1) * t754 - rSges(3,2) * t571) * t702;
t520 = t567 * rSges(3,3) + (rSges(3,1) * t571 + rSges(3,2) * t754) * t566;
t515 = Icges(3,3) * t567 + (Icges(3,5) * t571 + Icges(3,6) * t754) * t566;
t463 = Icges(3,1) * t541 + Icges(3,4) * t595 + Icges(3,5) * t740;
t462 = -Icges(3,1) * t594 + Icges(3,4) * t593 - Icges(3,5) * t738;
t461 = Icges(3,4) * t541 + Icges(3,2) * t595 + Icges(3,6) * t740;
t460 = -Icges(3,4) * t594 + Icges(3,2) * t593 - Icges(3,6) * t738;
t459 = Icges(3,5) * t541 + Icges(3,6) * t595 + Icges(3,3) * t740;
t458 = -Icges(3,5) * t594 + Icges(3,6) * t593 - Icges(3,3) * t738;
t444 = t465 + t705;
t443 = -t464 + t646;
t422 = -t567 * t464 - t520 * t738;
t421 = t465 * t567 - t520 * t740;
t389 = rSges(3,3) * t650 - t617;
t383 = Icges(3,1) * t483 - Icges(3,4) * t482 + Icges(3,5) * t650;
t382 = Icges(3,1) * t481 + Icges(3,4) * t480 + Icges(3,5) * t649;
t381 = Icges(3,4) * t483 - Icges(3,2) * t482 + Icges(3,6) * t650;
t380 = Icges(3,4) * t481 + Icges(3,2) * t480 + Icges(3,6) * t649;
t379 = Icges(3,5) * t483 - Icges(3,6) * t482 + Icges(3,3) * t650;
t378 = Icges(3,5) * t481 + Icges(3,6) * t480 + Icges(3,3) * t649;
t335 = (-t565 + (-rSges(3,3) - pkin(8)) * t740) * qJD(1) + t617;
t334 = t619 + t388;
t333 = t515 * t740 + t516 * t595 + t517 * t541;
t332 = -t515 * t738 + t516 * t593 - t517 * t594;
t323 = t705 - t713;
t322 = t646 + t714;
t319 = t588 * t567;
t318 = t567 * t388 + (-t520 * t703 - t527 * t572) * t566;
t317 = -t567 * t389 + (t520 * t704 - t527 * t574) * t566;
t316 = -t409 * t656 + t457 * t595;
t315 = t408 * t656 - t457 * t593;
t314 = t567 * t459 + (t461 * t754 + t463 * t571) * t566;
t313 = t567 * t458 + (t460 * t754 + t462 * t571) * t566;
t312 = t610 + t372;
t311 = -t593 * t748 + t596 + t615;
t294 = t459 * t740 + t461 * t595 + t463 * t541;
t293 = t458 * t740 + t460 * t595 + t462 * t541;
t292 = -t459 * t738 + t461 * t593 - t463 * t594;
t291 = -t458 * t738 + t460 * t593 - t462 * t594;
t290 = -t453 * t656 + t536 * t454 + t537 * t455;
t289 = t290 * t648;
t278 = -t408 * t595 + t409 * t593;
t273 = t566 * t635 + t567 * t714;
t272 = t409 * t567 + t707 * t740 + t473;
t255 = -t445 * t656 - t518 * t446 + t519 * t447;
t245 = t255 * t648;
t242 = t482 * rSges(5,3) - t616;
t227 = -t482 * t516 + t483 * t517 + t593 * t525 - t594 * t526 + (t515 * t704 - t524 * t574) * t566;
t226 = t480 * t516 + t481 * t517 + t595 * t525 + t541 * t526 + (t515 * t703 + t524 * t572) * t566;
t219 = (t408 * t572 + t409 * t574) * t566 + t706;
t216 = t590 + t310;
t215 = -t308 + t579;
t214 = -qJD(1) * t705 + t727;
t213 = t415 + t619 + t264;
t212 = t310 * t518 + t345 * t586;
t211 = -t308 * t518 - t345 * t585;
t210 = -t403 * t595 + t405 * t506 + t407 * t507;
t209 = -t402 * t595 + t404 * t506 + t406 * t507;
t208 = -t403 * t593 + t405 * t504 - t407 * t602;
t207 = -t402 * t593 + t404 * t504 - t406 * t602;
t204 = (-t372 * t754 - t651) * t566 + t710 * t595;
t203 = t371 * t656 - t448 * t593 + t721;
t202 = -t364 * t595 + t366 * t586 + t368 * t499;
t201 = -t363 * t595 + t365 * t586 + t367 * t499;
t200 = -t364 * t593 + t366 * t585 + t368 * t497;
t199 = -t363 * t593 + t365 * t585 + t367 * t497;
t198 = (-t486 + t719) * t567 + t566 * t618;
t197 = t372 * t567 + t661 * t740 + t720;
t196 = -t308 * t586 + t310 * t585;
t193 = t567 * t378 + (t754 * t380 + t382 * t571 + (-t461 * t571 + t463 * t754) * qJD(2)) * t566;
t192 = t567 * t379 + (t754 * t381 + t383 * t571 + (-t460 * t571 + t462 * t754) * qJD(2)) * t566;
t189 = t482 * t748 + t577 + t616;
t188 = t578 + t241;
t187 = t590 + t725;
t186 = t579 - t726;
t183 = -t371 * t595 - t593 * t718 - t331;
t175 = (t371 * t572 + t372 * t574) * t566 + t626;
t173 = t567 * t264 + t401 + (qJD(1) * t635 + t572 * t711) * t566;
t172 = t511 + t727 * t567 + (t457 * t704 + t574 * t711) * t566;
t171 = -t593 * t426 + t482 * t457 + (t265 * t754 - t408 * t701) * t566;
t170 = t595 * t426 + t480 * t457 + (-t264 * t754 + t409 * t701) * t566;
t153 = t518 * t725 + t586 * t722;
t152 = -t518 * t726 - t585 * t722;
t145 = (-t310 * t754 + t591) * t566 + t665 * t595;
t144 = t308 * t656 - t345 * t593 + t625;
t143 = (-t486 + t667) * t567 + t566 * t609;
t142 = t310 * t567 + t627 * t740 + t664;
t131 = t585 * t725 - t586 * t726;
t130 = t264 * t593 - t265 * t595 - t408 * t480 - t409 * t482;
t129 = -t308 * t595 - t593 * t666 + t724;
t126 = (t308 * t572 + t310 * t574) * t566 + t604;
t125 = (-t725 * t754 + t591) * t566 + t628 * t595;
t124 = -t593 * t722 + t656 * t726 + t625;
t123 = (-t486 + t630) * t567 + t566 * t597;
t122 = t567 * t725 + t611 * t740 + t664;
t121 = -t169 + t575;
t120 = t576 + t167;
t119 = (t264 * t574 + t265 * t572 + (t408 * t574 + t572 * t713) * qJD(1)) * t566 + t662;
t118 = t257 * t567 + (-t209 * t574 + t210 * t572) * t566;
t117 = t256 * t567 + (-t207 * t574 + t208 * t572) * t566;
t116 = t567 * t241 + (qJD(1) * t618 + t572 * t663) * t566 + t731;
t115 = (-t242 + t729) * t567 + (t448 * t704 + t574 * t663) * t566 + t709;
t114 = -t209 * t593 - t210 * t595 - t257 * t656;
t113 = -t207 * t593 - t208 * t595 - t256 * t656;
t110 = t225 * t567 + (-t201 * t574 + t202 * t572) * t566;
t109 = t224 * t567 + (-t199 * t574 + t200 * t572) * t566;
t108 = -t593 * t629 - t595 * t726 + t724;
t107 = (t572 * t726 + t574 * t725) * t566 + t604;
t106 = -t201 * t593 - t202 * t595 - t225 * t656;
t105 = -t199 * t593 - t200 * t595 - t224 * t656;
t104 = -t593 * t390 + t482 * t448 + (t242 * t754 + t701 * t719) * t566 + t671;
t103 = t342 - t715 * t595 + t710 * t480 + (-t241 * t754 + t372 * t701 - t652) * t566;
t100 = -t169 * t518 - t253 * t585 - t308 * t484 + t345 * t358;
t99 = t167 * t518 + t253 * t586 + t310 * t484 - t345 * t356;
t98 = -t258 * t593 + t260 * t504 - t262 * t602 + t398 * t405 + t399 * t407 + t403 * t482;
t97 = -t259 * t593 + t261 * t504 - t263 * t602 + t398 * t404 + t399 * t406 + t402 * t482;
t96 = -t258 * t595 + t260 * t506 + t262 * t507 + t396 * t405 + t397 * t407 - t403 * t480;
t95 = -t259 * t595 + t261 * t506 + t263 * t507 + t396 * t404 + t397 * t406 - t402 * t480;
t90 = -t235 * t593 + t237 * t585 + t239 * t497 - t358 * t366 + t359 * t368 + t364 * t482;
t89 = -t236 * t593 + t238 * t585 + t240 * t497 - t358 * t365 + t359 * t367 + t363 * t482;
t88 = -t235 * t595 + t237 * t586 + t239 * t499 - t356 * t366 + t357 * t368 - t364 * t480;
t87 = -t236 * t595 + t238 * t586 + t240 * t499 - t356 * t365 + t357 * t367 - t363 * t480;
t84 = t575 - t733;
t83 = t576 + t734;
t82 = -t242 * t595 - t371 * t480 - (-t241 - t243) * t593 + t718 * t482 + t732;
t79 = -t150 * t593 - t151 * t595 - t195 * t656;
t78 = -t148 * t593 - t149 * t595 - t194 * t656;
t77 = t167 * t585 - t169 * t586 + t308 * t356 - t310 * t358;
t76 = (t241 * t574 + t242 * t572 + (t371 * t574 + (-t487 + t718) * t572) * qJD(1)) * t566 + t606;
t71 = -t150 * t585 - t151 * t586 + t195 * t518;
t70 = -t148 * t585 - t149 * t586 + t194 * t518;
t69 = t567 * t167 + (qJD(1) * t609 + t572 * t631) * t566 + t670;
t68 = (-t169 + t669) * t567 + (t345 * t704 + t574 * t631) * t566 + t660;
t55 = -t593 * t253 + t482 * t345 + (t169 * t754 + t667 * t701) * t566 + t605;
t54 = -t668 * t595 + t665 * t480 + (-t167 * t754 + t310 * t701 + t592) * t566 + t723;
t53 = t358 * t722 - t484 * t726 - t518 * t733 - t585 * t728;
t52 = -t356 * t722 + t484 * t725 + t518 * t734 + t586 * t728;
t47 = t734 * t567 + (qJD(1) * t597 + t572 * t612) * t566 + t670;
t46 = (t669 - t733) * t567 + (t574 * t612 + t704 * t722) * t566 + t660;
t37 = -t169 * t595 - t308 * t480 - (-t167 + t730) * t593 + t666 * t482 + t633;
t36 = (t167 * t574 + t169 * t572 + (t308 * t574 + (-t487 + t666) * t572) * qJD(1)) * t566 + t587;
t35 = -t728 * t593 + t722 * t482 + (t630 * t701 + t733 * t754) * t566 + t605;
t34 = -t632 * t595 + t628 * t480 + (t725 * t701 - t734 * t754 + t592) * t566 + t723;
t33 = t356 * t726 - t358 * t725 + t585 * t734 - t586 * t733;
t31 = -t111 * t593 - t112 * t595 - t176 * t656 + t217 * t482 - t218 * t480 + t289;
t29 = t133 * t567 + (t572 * t98 - t574 * t97 + (t207 * t572 + t208 * t574) * qJD(1)) * t566;
t28 = t132 * t567 + (t572 * t96 - t574 * t95 + (t209 * t572 + t210 * t574) * qJD(1)) * t566;
t27 = -t101 * t593 - t102 * t595 - t147 * t656 + t205 * t482 - t206 * t480 + t245;
t26 = t207 * t482 - t208 * t480 - t97 * t593 - t98 * t595 + (-t133 * t754 + t256 * t701) * t566;
t25 = t209 * t482 - t210 * t480 - t95 * t593 - t96 * t595 + (-t132 * t754 + t257 * t701) * t566;
t24 = (t734 * t574 + t733 * t572 + (t726 * t574 + (-t487 + t629) * t572) * qJD(1)) * t566 + t587;
t23 = t128 * t567 + (t572 * t90 - t574 * t89 + (t199 * t572 + t200 * t574) * qJD(1)) * t566;
t22 = t127 * t567 + (t572 * t88 - t574 * t87 + (t201 * t572 + t202 * t574) * qJD(1)) * t566;
t21 = -t733 * t595 - t726 * t480 - (t730 - t734) * t593 + t629 * t482 + t633;
t20 = t199 * t482 - t200 * t480 - t89 * t593 - t90 * t595 + (-t128 * t754 + t224 * t701) * t566;
t19 = t201 * t482 - t202 * t480 - t87 * t593 - t88 * t595 + (-t127 * t754 + t225 * t701) * t566;
t30 = [(t120 * t216 + t121 * t215) * t639 + (t186 * t84 + t187 * t83) * t637 + (t188 * t312 + t189 * t311) * t641 + (t213 * t323 + t214 * t322) * t643 + (t334 * t444 + t335 * t443) * t645 + t588 - t760; t91 + t92 + t319 + m(3) * (t317 * t443 + t318 * t444 + t334 * t421 + t335 * t422) + (t172 * t322 + t173 * t323 + t213 * t272 + t214 * t273) * m(4) + (t115 * t311 + t116 * t312 + t188 * t197 + t189 * t198) * m(5) + (t120 * t142 + t121 * t143 + t215 * t68 + t216 * t69) * m(6) + (t122 * t83 + t123 * t84 + t186 * t46 + t187 * t47) * m(7) + ((-t192 / 0.2e1 - t227 / 0.2e1 - t582) * t574 + (t193 / 0.2e1 + t226 / 0.2e1 + t583) * t572 + ((t314 / 0.2e1 + t333 / 0.2e1 - t580) * t574 + (t313 / 0.2e1 + t332 / 0.2e1 + t581) * t572) * qJD(1)) * t566 + t768; (t422 * t317 + t421 * t318 + (t464 * t572 + t465 * t574) * (t388 * t574 + t389 * t572 + (t464 * t574 - t465 * t572) * qJD(1)) * t759) * t645 + (t107 * t24 + t122 * t47 + t123 * t46) * t637 + (t115 * t198 + t116 * t197 + t175 * t76) * t641 + (t119 * t219 + t172 * t273 + t173 * t272) * t643 + (t126 * t36 + t142 * t69 + t143 * t68) * t639 + (t28 + t22 + t10 + t9 + ((t380 * t595 + t541 * t382 + t480 * t461 + t481 * t463) * t572 + t294 * t703 - (t381 * t595 + t541 * t383 + t480 * t460 + t481 * t462) * t574 + t293 * t704 + ((t378 * t572 + t459 * t703) * t572 - (t379 * t572 + t458 * t703) * t574) * t566) * t566) * t740 + (-t29 - t23 - t12 - t11 - ((t380 * t593 - t382 * t594 - t482 * t461 + t483 * t463) * t572 + t292 * t703 - (t381 * t593 - t383 * t594 - t482 * t460 + t483 * t462) * t574 + t291 * t704 + ((-t378 * t574 + t459 * t704) * t572 - (-t379 * t574 + t458 * t704) * t574) * t566) * t566) * t738 + (t117 + t109 + t65 + t64 + (-t291 * t574 + t292 * t572) * t566) * t650 + ((-t293 * t574 + t294 * t572) * t566 + t118 + t110 + t67 + t66) * t649 + (t226 * t740 + t319 + (-t192 * t574 + t193 * t572 + (t313 * t572 + t314 * t574) * qJD(1)) * t566 - t227 * t738 + t332 * t650 + t333 * t649 + t761) * t567; t289 + t190 + t191 + t245 + (t170 * t323 + t171 * t322 + t213 * t316 + t214 * t315) * m(4) + (t103 * t312 + t104 * t311 + t188 * t204 + t189 * t203) * m(5) + (t120 * t145 + t121 * t144 + t215 * t55 + t216 * t54) * m(6) + (t124 * t84 + t125 * t83 + t186 * t35 + t187 * t34) * m(7) - t583 * t595 - t582 * t593 + t581 * t482 + t580 * t480 + t760 * t656; (t119 * t278 + t130 * t219 + t170 * t272 + t171 * t273 + t172 * t315 + t173 * t316) * m(4) + (t103 * t197 + t104 * t198 + t115 * t203 + t116 * t204 + t175 * t82 + t183 * t76) * m(5) + (t107 * t21 + t108 * t24 + t122 * t34 + t123 * t35 + t124 * t46 + t125 * t47) * m(7) + (t126 * t37 + t129 * t36 + t142 * t54 + t143 * t55 + t144 * t68 + t145 * t69) * m(6) + (t31 / 0.2e1 + t27 / 0.2e1 + t682) * t567 - (t28 / 0.2e1 + t22 / 0.2e1 + t691) * t595 - (t29 / 0.2e1 + t23 / 0.2e1 + t683) * t593 + (t109 / 0.2e1 + t117 / 0.2e1 + t677) * t482 + (-t110 / 0.2e1 - t118 / 0.2e1 - t676) * t480 + ((-t26 / 0.2e1 - t20 / 0.2e1 - t693) * t574 + (t25 / 0.2e1 + t19 / 0.2e1 + t694) * t572 + ((t762 * t572 / 0.2e1 - t763 * t574 / 0.2e1) * t566 + (t290 / 0.2e1 + t255 / 0.2e1) * t567 + t675) * t701 + ((t106 / 0.2e1 + t114 / 0.2e1 + t678) * t574 + (t105 / 0.2e1 + t113 / 0.2e1 + t679) * t572) * qJD(1) + t761 * t692) * t566; (t130 * t278 + t170 * t316 + t171 * t315) * t643 + (t103 * t204 + t104 * t203 + t183 * t82) * t641 + (t108 * t21 + t124 * t35 + t125 * t34) * t637 + (t129 * t37 + t144 * t55 + t145 * t54) * t639 - (t25 + t19 + t6 + t5) * t595 - (t26 + t20 + t8 + t7) * t593 + (t60 + t61 + t105 + t113) * t482 + (-t62 - t63 - t106 - t114) * t480 + ((t78 + t79 - t762 * t595 - t763 * t593 + (-t255 - t290) * t656) * t701 + (-t15 - t16 - t27 - t31) * t754) * t566; -(m(5) * t189 + m(6) * t121 + m(7) * t84) * t595 - (m(5) * t188 + m(6) * t120 + m(7) * t83) * t593 + (m(5) * t312 + m(6) * t216 + m(7) * t187) * t482 + (-m(5) * t311 - m(6) * t215 - m(7) * t186) * t480; (-t24 * t688 - t36 * t689 - t690 * t76) * t566 - (m(5) * t115 + m(6) * t68 + m(7) * t46) * t595 - (m(5) * t116 + m(6) * t69 + m(7) * t47) * t593 + (m(5) * t197 + m(6) * t142 + m(7) * t122) * t482 + (-m(5) * t198 - m(6) * t143 - m(7) * t123) * t480 + (m(5) * t175 + m(6) * t126 + m(7) * t107) * t648; (-t21 * t688 - t37 * t689 - t690 * t82) * t566 - (m(5) * t104 + m(6) * t55 + m(7) * t35) * t595 - (m(5) * t103 + m(6) * t54 + m(7) * t34) * t593 + (m(5) * t204 + m(6) * t145 + m(7) * t125) * t482 + (-m(5) * t203 - m(6) * t144 - m(7) * t124) * t480 + (m(5) * t183 + m(6) * t129 + m(7) * t108) * t648; 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + t758) * (-t571 * t647 * t759 + t480 * t595 - t482 * t593); (t152 * t84 + t153 * t83 + t186 * t53 + t187 * t52) * m(7) + (t100 * t215 + t120 * t212 + t121 * t211 + t216 * t99) * m(6) - t607 * t586 - t608 * t585 + t599 * t358 + t598 * t356 + t745 + t746; (t107 * t33 + t122 * t52 + t123 * t53 + t131 * t24 + t152 * t46 + t153 * t47) * m(7) + (t100 * t143 + t126 * t77 + t142 * t99 + t196 * t36 + t211 * t68 + t212 * t69) * m(6) + (t13 / 0.2e1 + t14 / 0.2e1) * t567 + (t17 / 0.2e1 + t18 / 0.2e1) * t518 - t691 * t586 - t683 * t585 + t675 * t484 + t677 * t358 + t676 * t356 + (-t695 * t574 + t696 * t572 + (t572 * t681 + t574 * t680) * qJD(1)) * t566; (t100 * t144 + t129 * t77 + t145 * t99 + t196 * t37 + t211 * t55 + t212 * t54) * m(6) + (t108 * t33 + t124 * t53 + t125 * t52 + t131 * t21 + t152 * t35 + t153 * t34) * m(7) - t696 * t595 - t695 * t593 + t682 * t518 - t694 * t586 - t693 * t585 + (t78 / 0.2e1 + t79 / 0.2e1) * t484 + t681 * t482 - t680 * t480 + t679 * t358 + t678 * t356 + ((t70 / 0.2e1 + t71 / 0.2e1) * t701 + t764 * t692) * t566; (-t33 * t688 - t689 * t77) * t566 - (m(6) * t100 + m(7) * t53) * t595 - (m(6) * t99 + m(7) * t52) * t593 + (m(6) * t212 + m(7) * t153) * t482 + (-m(6) * t211 - m(7) * t152) * t480 + (m(6) * t196 + m(7) * t131) * t648; (t131 * t33 + t152 * t53 + t153 * t52) * t637 + (t100 * t211 + t196 * t77 + t212 * t99) * t639 + t764 * t518 - (t1 + t2) * t586 - (t3 + t4) * t585 + (t71 + t70) * t484 + (t56 + t57) * t358 + (t58 + t59) * t356; (t186 * t274 + t187 * t276 - t600 * t84 - t601 * t83) * m(7); (t107 * t376 + t122 * t276 + t123 * t274 - t24 * t589 - t46 * t600 - t47 * t601) * m(7); (t108 * t376 + t124 * t274 + t125 * t276 - t21 * t589 - t34 * t601 - t35 * t600) * m(7); 0.2e1 * (-t274 * t595 - t276 * t593 - t601 * t482 + t600 * t480 + (-t376 * t754 - t589 * t701) * t566) * t758; (t131 * t376 + t152 * t274 + t153 * t276 - t33 * t589 - t52 * t601 - t53 * t600) * m(7); (-t274 * t600 - t276 * t601 - t376 * t589) * t637;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t30(1) t30(2) t30(4) t30(7) t30(11) t30(16); t30(2) t30(3) t30(5) t30(8) t30(12) t30(17); t30(4) t30(5) t30(6) t30(9) t30(13) t30(18); t30(7) t30(8) t30(9) t30(10) t30(14) t30(19); t30(11) t30(12) t30(13) t30(14) t30(15) t30(20); t30(16) t30(17) t30(18) t30(19) t30(20) t30(21);];
Mq  = res;
