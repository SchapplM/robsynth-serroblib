% Calculate time derivative of joint inertia matrix for
% S6RRRPRP6
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP6_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:56
% EndTime: 2019-03-09 16:55:13
% DurationCPUTime: 44.10s
% Computational Cost: add. (132668->1567), mult. (253643->2059), div. (0->0), fcn. (281790->12), ass. (0->603)
t754 = -qJ(6) - pkin(10) - rSges(7,3);
t700 = qJ(3) + pkin(11);
t564 = sin(t700);
t567 = cos(pkin(6));
t572 = sin(qJ(2));
t566 = sin(pkin(6));
t639 = cos(t700);
t614 = t566 * t639;
t521 = t564 * t567 + t572 * t614;
t763 = cos(qJ(2));
t650 = qJD(2) * t763;
t623 = t566 * t650;
t485 = qJD(3) * t521 + t564 * t623;
t744 = t566 * t572;
t520 = t564 * t744 - t567 * t639;
t486 = -qJD(3) * t520 + t614 * t650;
t704 = qJD(2) * t572;
t651 = t566 * t704;
t387 = Icges(5,5) * t486 - Icges(5,6) * t485 + Icges(5,3) * t651;
t388 = Icges(5,4) * t486 - Icges(5,2) * t485 + Icges(5,6) * t651;
t389 = Icges(5,1) * t486 - Icges(5,4) * t485 + Icges(5,5) * t651;
t659 = t566 * t763;
t445 = Icges(5,5) * t521 - Icges(5,6) * t520 - Icges(5,3) * t659;
t446 = Icges(5,4) * t521 - Icges(5,2) * t520 - Icges(5,6) * t659;
t447 = Icges(5,1) * t521 - Icges(5,4) * t520 - Icges(5,5) * t659;
t149 = -t387 * t659 - t388 * t520 + t389 * t521 + t445 * t651 - t446 * t485 + t447 * t486;
t571 = sin(qJ(3));
t740 = t567 * t571;
t575 = cos(qJ(3));
t742 = t566 * t575;
t538 = t572 * t742 + t740;
t503 = -qJD(3) * t538 - t571 * t623;
t537 = t567 * t575 - t571 * t744;
t504 = qJD(3) * t537 + t575 * t623;
t425 = Icges(4,5) * t504 + Icges(4,6) * t503 + Icges(4,3) * t651;
t426 = Icges(4,4) * t504 + Icges(4,2) * t503 + Icges(4,6) * t651;
t427 = Icges(4,1) * t504 + Icges(4,4) * t503 + Icges(4,5) * t651;
t453 = Icges(4,5) * t538 + Icges(4,6) * t537 - Icges(4,3) * t659;
t454 = Icges(4,4) * t538 + Icges(4,2) * t537 - Icges(4,6) * t659;
t455 = Icges(4,1) * t538 + Icges(4,4) * t537 - Icges(4,5) * t659;
t178 = -t425 * t659 + t426 * t537 + t427 * t538 + t453 * t651 + t454 * t503 + t455 * t504;
t785 = t149 + t178;
t576 = cos(qJ(1));
t657 = t576 * t763;
t573 = sin(qJ(1));
t739 = t573 * t572;
t542 = -t567 * t739 + t657;
t626 = t567 * t657;
t594 = t626 - t739;
t484 = qJD(1) * t542 + qJD(2) * t594;
t658 = t573 * t763;
t738 = t576 * t572;
t595 = -t567 * t738 - t658;
t588 = t564 * t595 - t576 * t614;
t707 = qJD(1) * t573;
t653 = t566 * t707;
t361 = qJD(3) * t588 + t484 * t639 + t564 * t653;
t741 = t566 * t576;
t498 = -t564 * t741 - t595 * t639;
t574 = cos(qJ(5));
t570 = sin(qJ(5));
t747 = t594 * t570;
t440 = t498 * t574 - t747;
t596 = -t567 * t658 - t738;
t483 = -qJD(1) * t596 - qJD(2) * t595;
t275 = -qJD(5) * t440 - t361 * t570 + t483 * t574;
t610 = t498 * t570 + t574 * t594;
t768 = qJD(5) * t610 - t483 * t570;
t276 = t361 * t574 - t768;
t603 = qJD(1) * t614;
t360 = qJD(3) * t498 + t484 * t564 - t573 * t603;
t784 = t276 * rSges(7,1) + t275 * rSges(7,2) - pkin(5) * t768 - qJD(6) * t588 - t360 * t754;
t783 = -rSges(7,1) * t440 + rSges(7,2) * t610 - t588 * t754;
t365 = Icges(5,5) * t498 + Icges(5,6) * t588 - Icges(5,3) * t594;
t367 = Icges(5,4) * t498 + Icges(5,2) * t588 - Icges(5,6) * t594;
t369 = Icges(5,1) * t498 + Icges(5,4) * t588 - Icges(5,5) * t594;
t203 = -t365 * t659 - t367 * t520 + t369 * t521;
t505 = t571 * t595 - t575 * t741;
t676 = t571 * t741;
t602 = t575 * t595 + t676;
t404 = -Icges(4,5) * t602 + Icges(4,6) * t505 - Icges(4,3) * t594;
t406 = -Icges(4,4) * t602 + Icges(4,2) * t505 - Icges(4,6) * t594;
t408 = -Icges(4,1) * t602 + Icges(4,4) * t505 - Icges(4,5) * t594;
t217 = -t404 * t659 + t406 * t537 + t408 * t538;
t772 = t203 + t217;
t743 = t566 * t573;
t500 = t542 * t639 + t564 * t743;
t589 = -t542 * t564 + t573 * t614;
t366 = Icges(5,5) * t500 + Icges(5,6) * t589 - Icges(5,3) * t596;
t368 = Icges(5,4) * t500 + Icges(5,2) * t589 - Icges(5,6) * t596;
t370 = Icges(5,1) * t500 + Icges(5,4) * t589 - Icges(5,5) * t596;
t204 = -t366 * t659 - t368 * t520 + t370 * t521;
t507 = -t542 * t571 + t573 * t742;
t677 = t571 * t743;
t508 = t542 * t575 + t677;
t405 = Icges(4,5) * t508 + Icges(4,6) * t507 - Icges(4,3) * t596;
t407 = Icges(4,4) * t508 + Icges(4,2) * t507 - Icges(4,6) * t596;
t409 = Icges(4,1) * t508 + Icges(4,4) * t507 - Icges(4,5) * t596;
t218 = -t405 * t659 + t407 * t537 + t409 * t538;
t771 = t204 + t218;
t441 = -t500 * t570 - t574 * t596;
t481 = -qJD(1) * t626 - t576 * t650 + (qJD(2) * t567 + qJD(1)) * t739;
t781 = qJD(5) * t441 - t481 * t570;
t780 = -pkin(10) - t754;
t779 = t785 * t567;
t482 = qJD(1) * t595 + qJD(2) * t596;
t706 = qJD(1) * t576;
t652 = t566 * t706;
t359 = qJD(3) * t589 + t482 * t639 + t564 * t652;
t746 = t596 * t570;
t442 = t500 * t574 - t746;
t273 = -qJD(5) * t442 - t359 * t570 - t481 * t574;
t274 = t359 * t574 + t781;
t358 = qJD(3) * t500 + t482 * t564 - t576 * t603;
t562 = pkin(5) * t574 + pkin(4);
t778 = t274 * rSges(7,1) + t273 * rSges(7,2) + pkin(5) * t781 - qJD(6) * t589 - t358 * t754 + t359 * t562;
t777 = t442 * rSges(7,1) + t441 * rSges(7,2) - pkin(5) * t746 + t500 * t562 + t589 * t754;
t296 = Icges(7,5) * t440 - Icges(7,6) * t610 - Icges(7,3) * t588;
t300 = Icges(7,4) * t440 - Icges(7,2) * t610 - Icges(7,6) * t588;
t304 = Icges(7,1) * t440 - Icges(7,4) * t610 - Icges(7,5) * t588;
t495 = -t521 * t570 - t574 * t659;
t627 = t570 * t659;
t597 = -t521 * t574 + t627;
t150 = t296 * t520 + t300 * t495 - t304 * t597;
t297 = Icges(7,5) * t442 + Icges(7,6) * t441 - Icges(7,3) * t589;
t301 = Icges(7,4) * t442 + Icges(7,2) * t441 - Icges(7,6) * t589;
t305 = Icges(7,1) * t442 + Icges(7,4) * t441 - Icges(7,5) * t589;
t151 = t297 * t520 + t301 * t495 - t305 * t597;
t157 = Icges(7,5) * t276 + Icges(7,6) * t275 + Icges(7,3) * t360;
t161 = Icges(7,4) * t276 + Icges(7,2) * t275 + Icges(7,6) * t360;
t165 = Icges(7,1) * t276 + Icges(7,4) * t275 + Icges(7,5) * t360;
t378 = qJD(5) * t597 - t486 * t570 + t574 * t651;
t580 = qJD(5) * t495 + t570 * t651;
t379 = t486 * t574 + t580;
t50 = t157 * t520 + t161 * t495 - t165 * t597 + t296 * t485 + t300 * t378 + t304 * t379;
t156 = Icges(7,5) * t274 + Icges(7,6) * t273 + Icges(7,3) * t358;
t160 = Icges(7,4) * t274 + Icges(7,2) * t273 + Icges(7,6) * t358;
t164 = Icges(7,1) * t274 + Icges(7,4) * t273 + Icges(7,5) * t358;
t51 = t156 * t520 + t160 * t495 - t164 * t597 + t297 * t485 + t301 * t378 + t305 * t379;
t337 = -Icges(7,5) * t597 + Icges(7,6) * t495 + Icges(7,3) * t520;
t339 = -Icges(7,4) * t597 + Icges(7,2) * t495 + Icges(7,6) * t520;
t341 = -Icges(7,1) * t597 + Icges(7,4) * t495 + Icges(7,5) * t520;
t192 = t337 * t520 + t339 * t495 - t341 * t597;
t247 = Icges(7,5) * t379 + Icges(7,6) * t378 + Icges(7,3) * t485;
t249 = Icges(7,4) * t379 + Icges(7,2) * t378 + Icges(7,6) * t485;
t251 = Icges(7,1) * t379 + Icges(7,4) * t378 + Icges(7,5) * t485;
t91 = t247 * t520 + t249 * t495 - t251 * t597 + t337 * t485 + t339 * t378 + t341 * t379;
t756 = t192 * t485 + t520 * t91;
t13 = t150 * t360 + t151 * t358 - t50 * t588 - t51 * t589 + t756;
t298 = Icges(6,5) * t440 - Icges(6,6) * t610 - Icges(6,3) * t588;
t302 = Icges(6,4) * t440 - Icges(6,2) * t610 - Icges(6,6) * t588;
t306 = Icges(6,1) * t440 - Icges(6,4) * t610 - Icges(6,5) * t588;
t152 = t298 * t520 + t302 * t495 - t306 * t597;
t299 = Icges(6,5) * t442 + Icges(6,6) * t441 - Icges(6,3) * t589;
t303 = Icges(6,4) * t442 + Icges(6,2) * t441 - Icges(6,6) * t589;
t307 = Icges(6,1) * t442 + Icges(6,4) * t441 - Icges(6,5) * t589;
t153 = t299 * t520 + t303 * t495 - t307 * t597;
t159 = Icges(6,5) * t276 + Icges(6,6) * t275 + Icges(6,3) * t360;
t163 = Icges(6,4) * t276 + Icges(6,2) * t275 + Icges(6,6) * t360;
t167 = Icges(6,1) * t276 + Icges(6,4) * t275 + Icges(6,5) * t360;
t52 = t159 * t520 + t163 * t495 - t167 * t597 + t298 * t485 + t302 * t378 + t306 * t379;
t158 = Icges(6,5) * t274 + Icges(6,6) * t273 + Icges(6,3) * t358;
t162 = Icges(6,4) * t274 + Icges(6,2) * t273 + Icges(6,6) * t358;
t166 = Icges(6,1) * t274 + Icges(6,4) * t273 + Icges(6,5) * t358;
t53 = t158 * t520 + t162 * t495 - t166 * t597 + t299 * t485 + t303 * t378 + t307 * t379;
t338 = -Icges(6,5) * t597 + Icges(6,6) * t495 + Icges(6,3) * t520;
t340 = -Icges(6,4) * t597 + Icges(6,2) * t495 + Icges(6,6) * t520;
t342 = -Icges(6,1) * t597 + Icges(6,4) * t495 + Icges(6,5) * t520;
t193 = t338 * t520 + t340 * t495 - t342 * t597;
t248 = Icges(6,5) * t379 + Icges(6,6) * t378 + Icges(6,3) * t485;
t250 = Icges(6,4) * t379 + Icges(6,2) * t378 + Icges(6,6) * t485;
t252 = Icges(6,1) * t379 + Icges(6,4) * t378 + Icges(6,5) * t485;
t92 = t248 * t520 + t250 * t495 - t252 * t597 + t338 * t485 + t340 * t378 + t342 * t379;
t755 = t193 * t485 + t520 * t92;
t14 = t152 * t360 + t153 * t358 - t52 * t588 - t53 * t589 + t755;
t773 = t13 + t14;
t282 = pkin(4) * t359 + pkin(10) * t358;
t737 = -t282 + t778;
t357 = t360 * pkin(10);
t761 = -pkin(4) + t562;
t736 = t361 * t761 - t357 + t784;
t493 = t588 * pkin(10);
t729 = -pkin(5) * t747 + t498 * t761 + t493 - t783;
t430 = pkin(4) * t500 - pkin(10) * t589;
t728 = -t430 + t777;
t565 = t576 * pkin(1);
t708 = pkin(8) * t743 + t565;
t237 = Icges(5,5) * t361 - Icges(5,6) * t360 + Icges(5,3) * t483;
t239 = Icges(5,4) * t361 - Icges(5,2) * t360 + Icges(5,6) * t483;
t241 = Icges(5,1) * t361 - Icges(5,4) * t360 + Icges(5,5) * t483;
t101 = -t520 * t239 + t521 * t241 - t485 * t367 + t486 * t369 + (-t237 * t763 + t365 * t704) * t566;
t236 = Icges(5,5) * t359 - Icges(5,6) * t358 - Icges(5,3) * t481;
t238 = Icges(5,4) * t359 - Icges(5,2) * t358 - Icges(5,6) * t481;
t240 = Icges(5,1) * t359 - Icges(5,4) * t358 - Icges(5,5) * t481;
t102 = -t520 * t238 + t521 * t240 - t485 * t368 + t486 * t370 + (-t236 * t763 + t366 * t704) * t566;
t400 = qJD(3) * t602 - t484 * t571 + t575 * t653;
t624 = t571 * t653;
t401 = qJD(3) * t505 + t484 * t575 + t624;
t259 = Icges(4,5) * t401 + Icges(4,6) * t400 + Icges(4,3) * t483;
t261 = Icges(4,4) * t401 + Icges(4,2) * t400 + Icges(4,6) * t483;
t263 = Icges(4,1) * t401 + Icges(4,4) * t400 + Icges(4,5) * t483;
t111 = t537 * t261 + t538 * t263 + t503 * t406 + t504 * t408 + (-t259 * t763 + t404 * t704) * t566;
t398 = -qJD(3) * t508 - t482 * t571 + t575 * t652;
t625 = t571 * t652;
t399 = qJD(3) * t507 + t482 * t575 + t625;
t258 = Icges(4,5) * t399 + Icges(4,6) * t398 - Icges(4,3) * t481;
t260 = Icges(4,4) * t399 + Icges(4,2) * t398 - Icges(4,6) * t481;
t262 = Icges(4,1) * t399 + Icges(4,4) * t398 - Icges(4,5) * t481;
t112 = t537 * t260 + t538 * t262 + t503 * t407 + t504 * t409 + (-t258 * t763 + t405 * t704) * t566;
t89 = t91 * t567;
t17 = t89 + (-t50 * t576 + t51 * t573 + (t150 * t573 + t151 * t576) * qJD(1)) * t566;
t90 = t92 * t567;
t18 = t90 + (-t52 * t576 + t53 * t573 + (t152 * t573 + t153 * t576) * qJD(1)) * t566;
t770 = t17 + t18 + t779 + ((-t101 - t111) * t576 + (t102 + t112) * t573 + (t573 * t772 + t576 * t771) * qJD(1)) * t566;
t769 = -t91 - t92 - t785;
t767 = t566 ^ 2;
t766 = m(7) / 0.2e1;
t563 = pkin(3) * t575 + pkin(2);
t762 = -pkin(2) + t563;
t759 = pkin(3) * qJD(3);
t569 = -qJ(4) - pkin(9);
t757 = -rSges(5,3) + t569;
t753 = Icges(3,4) * t572;
t751 = t483 * t569;
t748 = t594 * t569;
t479 = t483 * pkin(9);
t687 = t575 * t759;
t637 = t566 * t687;
t688 = t571 * t759;
t660 = -qJD(4) * t594 - t576 * t637 + t595 * t688;
t245 = pkin(3) * t624 + t484 * t762 - t479 + t660 - t751;
t535 = t594 * pkin(9);
t553 = pkin(3) * t676;
t376 = -t595 * t762 + t535 - t553 + t748;
t735 = -t245 * t596 - t376 * t481;
t734 = rSges(7,1) * t379 + rSges(7,2) * t378 + pkin(5) * t580 + t520 * qJD(6) + t485 * t780 + t486 * t761;
t417 = pkin(2) * t482 - pkin(9) * t481;
t586 = pkin(3) * t625 - qJD(4) * t596 + t481 * t569 + t482 * t563 - t542 * t688 + t573 * t637;
t244 = -t417 + t586;
t403 = t567 * t417;
t733 = t244 * t567 + t403;
t732 = -t244 - t282;
t418 = pkin(2) * t484 + t479;
t731 = -t245 - t418;
t265 = rSges(4,1) * t401 + rSges(4,2) * t400 + rSges(4,3) * t483;
t730 = -t265 - t418;
t727 = -rSges(7,1) * t597 + rSges(7,2) * t495 - pkin(5) * t627 + t520 * t780 + t521 * t761;
t331 = t596 * t376;
t429 = pkin(4) * t498 - t493;
t726 = -t429 * t596 - t331;
t488 = pkin(2) * t542 - pkin(9) * t596;
t661 = pkin(3) * t677 + t542 * t563 + t569 * t596;
t377 = -t488 + t661;
t343 = t377 * t651;
t725 = t430 * t651 + t343;
t689 = t763 * pkin(9);
t456 = pkin(3) * t740 + (t569 * t763 + t572 * t762 + t689) * t566;
t724 = t376 * t659 - t456 * t594;
t473 = t567 * t488;
t723 = t377 * t567 + t473;
t618 = -rSges(5,1) * t498 - rSges(5,2) * t588;
t373 = -rSges(5,3) * t594 - t618;
t722 = -t373 - t376;
t374 = rSges(5,1) * t500 + rSges(5,2) * t589 - rSges(5,3) * t596;
t721 = -t374 - t377;
t720 = -t376 - t429;
t719 = -t377 - t430;
t392 = rSges(5,1) * t486 - rSges(5,2) * t485 + rSges(5,3) * t651;
t690 = t763 * pkin(2);
t431 = t567 * t687 + (-t572 * t688 - t763 * qJD(4) + (-t690 + t763 * t563 + (-pkin(9) - t569) * t572) * qJD(2)) * t566;
t718 = -t392 - t431;
t410 = -rSges(4,1) * t602 + rSges(4,2) * t505 - rSges(4,3) * t594;
t487 = -pkin(2) * t595 - t535;
t717 = -t410 - t487;
t411 = rSges(4,1) * t508 + rSges(4,2) * t507 - rSges(4,3) * t596;
t716 = -t411 - t488;
t419 = pkin(4) * t486 + t485 * pkin(10);
t715 = -t419 - t431;
t428 = rSges(4,1) * t504 + rSges(4,2) * t503 + rSges(4,3) * t651;
t705 = qJD(2) * t566;
t529 = (pkin(9) * t572 + t690) * t705;
t714 = -t428 - t529;
t448 = rSges(5,1) * t521 - rSges(5,2) * t520 - rSges(5,3) * t659;
t713 = t448 + t456;
t543 = (pkin(2) * t572 - t689) * t566;
t512 = t543 * t653;
t712 = t456 * t653 + t512;
t452 = pkin(4) * t521 + pkin(10) * t520;
t711 = t452 + t456;
t457 = rSges(4,1) * t538 + rSges(4,2) * t537 - rSges(4,3) * t659;
t710 = -t457 - t543;
t709 = t487 * t743 + t488 * t741;
t140 = -t296 * t589 + t300 * t441 + t304 * t442;
t141 = -t297 * t589 + t301 * t441 + t305 * t442;
t181 = -t337 * t589 + t339 * t441 + t341 * t442;
t38 = -t157 * t589 + t161 * t441 + t165 * t442 + t273 * t300 + t274 * t304 + t296 * t358;
t39 = -t156 * t589 + t160 * t441 + t164 * t442 + t273 * t301 + t274 * t305 + t297 * t358;
t72 = -t247 * t589 + t249 * t441 + t251 * t442 + t273 * t339 + t274 * t341 + t337 * t358;
t1 = t140 * t360 + t141 * t358 + t181 * t485 - t38 * t588 - t39 * t589 + t520 * t72;
t142 = -t298 * t589 + t302 * t441 + t306 * t442;
t143 = -t299 * t589 + t303 * t441 + t307 * t442;
t182 = -t338 * t589 + t340 * t441 + t342 * t442;
t40 = -t159 * t589 + t163 * t441 + t167 * t442 + t273 * t302 + t274 * t306 + t298 * t358;
t41 = -t158 * t589 + t162 * t441 + t166 * t442 + t273 * t303 + t274 * t307 + t299 * t358;
t73 = -t248 * t589 + t250 * t441 + t252 * t442 + t273 * t340 + t274 * t342 + t338 * t358;
t2 = t142 * t360 + t143 * t358 + t182 * t485 - t40 * t588 - t41 * t589 + t520 * t73;
t699 = t2 / 0.2e1 + t1 / 0.2e1;
t136 = -t296 * t588 - t300 * t610 + t304 * t440;
t137 = -t297 * t588 - t301 * t610 + t305 * t440;
t179 = -t337 * t588 - t339 * t610 + t341 * t440;
t42 = -t157 * t588 - t161 * t610 + t165 * t440 + t275 * t300 + t276 * t304 + t296 * t360;
t43 = -t156 * t588 - t160 * t610 + t164 * t440 + t275 * t301 + t276 * t305 + t297 * t360;
t74 = -t247 * t588 - t249 * t610 + t251 * t440 + t275 * t339 + t276 * t341 + t337 * t360;
t3 = t136 * t360 + t137 * t358 + t179 * t485 - t42 * t588 - t43 * t589 + t520 * t74;
t138 = -t298 * t588 - t302 * t610 + t306 * t440;
t139 = -t299 * t588 - t303 * t610 + t307 * t440;
t180 = -t338 * t588 - t340 * t610 + t342 * t440;
t44 = -t159 * t588 - t163 * t610 + t167 * t440 + t275 * t302 + t276 * t306 + t298 * t360;
t45 = -t158 * t588 - t162 * t610 + t166 * t440 + t275 * t303 + t276 * t307 + t299 * t360;
t75 = -t248 * t588 - t250 * t610 + t252 * t440 + t275 * t340 + t276 * t342 + t338 * t360;
t4 = t138 * t360 + t139 * t358 + t180 * t485 - t44 * t588 - t45 * t589 + t520 * t75;
t698 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t140 * t483 - t141 * t481 - t38 * t594 - t39 * t596 + (t181 * t704 - t72 * t763) * t566;
t6 = t142 * t483 - t143 * t481 - t40 * t594 - t41 * t596 + (t182 * t704 - t73 * t763) * t566;
t697 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t136 * t483 - t137 * t481 - t42 * t594 - t43 * t596 + (t179 * t704 - t74 * t763) * t566;
t8 = t138 * t483 - t139 * t481 - t44 * t594 - t45 * t596 + (t180 * t704 - t75 * t763) * t566;
t696 = t8 / 0.2e1 + t7 / 0.2e1;
t695 = -t763 / 0.2e1;
t10 = t73 * t567 + (-t40 * t576 + t41 * t573 + (t142 * t573 + t143 * t576) * qJD(1)) * t566;
t9 = t72 * t567 + (-t38 * t576 + t39 * t573 + (t140 * t573 + t141 * t576) * qJD(1)) * t566;
t694 = t10 / 0.2e1 + t9 / 0.2e1;
t693 = m(5) * t763;
t692 = m(6) * t763;
t691 = m(7) * t763;
t11 = t74 * t567 + (-t42 * t576 + t43 * t573 + (t136 * t573 + t137 * t576) * qJD(1)) * t566;
t12 = t75 * t567 + (-t44 * t576 + t45 * t573 + (t138 * t573 + t139 * t576) * qJD(1)) * t566;
t686 = t12 / 0.2e1 + t11 / 0.2e1;
t188 = t192 * t651;
t15 = t150 * t483 - t151 * t481 - t50 * t594 - t51 * t596 - t659 * t91 + t188;
t189 = t193 * t651;
t16 = t152 * t483 - t153 * t481 - t52 * t594 - t53 * t596 - t659 * t92 + t189;
t685 = t15 / 0.2e1 + t16 / 0.2e1;
t56 = -t136 * t588 - t137 * t589 + t179 * t520;
t57 = -t138 * t588 - t139 * t589 + t180 * t520;
t684 = t56 / 0.2e1 + t57 / 0.2e1;
t58 = -t140 * t588 - t141 * t589 + t181 * t520;
t59 = -t142 * t588 - t143 * t589 + t182 * t520;
t683 = t59 / 0.2e1 + t58 / 0.2e1;
t60 = -t136 * t594 - t137 * t596 - t179 * t659;
t61 = -t138 * t594 - t139 * t596 - t180 * t659;
t682 = t61 / 0.2e1 + t60 / 0.2e1;
t62 = -t140 * t594 - t141 * t596 - t181 * t659;
t63 = -t142 * t594 - t143 * t596 - t182 * t659;
t681 = t63 / 0.2e1 + t62 / 0.2e1;
t64 = t179 * t567 + (-t136 * t576 + t137 * t573) * t566;
t65 = t180 * t567 + (-t138 * t576 + t139 * t573) * t566;
t680 = t65 / 0.2e1 + t64 / 0.2e1;
t66 = t181 * t567 + (-t140 * t576 + t141 * t573) * t566;
t67 = t182 * t567 + (-t142 * t576 + t143 * t573) * t566;
t679 = t67 / 0.2e1 + t66 / 0.2e1;
t678 = ((-t150 - t152) * t576 + (t151 + t153) * t573) * t566 / 0.2e1 + (t192 + t193) * t567 / 0.2e1;
t674 = t245 * t659 - t431 * t594 + t456 * t483;
t673 = t282 * t567 + t733;
t283 = t361 * pkin(4) + t357;
t672 = -t283 + t731;
t254 = rSges(6,1) * t379 + rSges(6,2) * t378 + rSges(6,3) * t485;
t671 = -t254 + t715;
t169 = rSges(6,1) * t274 + rSges(6,2) * t273 + rSges(6,3) * t358;
t309 = rSges(6,1) * t440 - rSges(6,2) * t610 - rSges(6,3) * t588;
t670 = -t309 + t720;
t311 = rSges(6,1) * t442 + rSges(6,2) * t441 - rSges(6,3) * t589;
t669 = -t311 + t719;
t347 = -rSges(6,1) * t597 + rSges(6,2) * t495 + rSges(6,3) * t520;
t668 = t347 + t711;
t242 = rSges(5,1) * t359 - rSges(5,2) * t358 - rSges(5,3) * t481;
t667 = t430 * t567 + t723;
t666 = -t529 + t718;
t264 = rSges(4,1) * t399 + rSges(4,2) * t398 - rSges(4,3) * t481;
t665 = t417 * t741 + t418 * t743 + t487 * t652;
t664 = -t543 - t713;
t663 = t452 * t653 + t712;
t390 = rSges(3,1) * t482 + rSges(3,2) * t481 + rSges(3,3) * t652;
t465 = rSges(3,1) * t542 + rSges(3,2) * t596 + rSges(3,3) * t743;
t656 = t763 * Icges(3,4);
t655 = t763 * t244;
t654 = t763 * t377;
t649 = -t573 * pkin(1) + pkin(8) * t741;
t648 = 2 * m(3);
t646 = 2 * m(4);
t644 = 0.2e1 * m(5);
t642 = 0.2e1 * m(6);
t640 = 0.2e1 * m(7);
t638 = t576 * t710;
t636 = -t283 * t596 - t429 * t481 + t735;
t635 = t715 - t734;
t634 = -t529 + t671;
t633 = t720 - t729;
t632 = t719 - t728;
t631 = t711 + t727;
t630 = -t543 - t668;
t629 = t376 * t743 + t377 * t741 + t709;
t628 = t429 * t659 - t452 * t594 + t724;
t622 = -pkin(1) * t707 + pkin(8) * t652;
t621 = t664 * t576;
t620 = -t484 * rSges(3,1) + t483 * rSges(3,2);
t619 = -rSges(5,1) * t361 + rSges(5,2) * t360;
t613 = -t529 + t635;
t612 = -t543 - t631;
t611 = t661 + t708;
t609 = t630 * t576;
t608 = t50 / 0.2e1 + t52 / 0.2e1 + t74 / 0.2e1 + t75 / 0.2e1;
t607 = t51 / 0.2e1 + t53 / 0.2e1 + t73 / 0.2e1 + t72 / 0.2e1;
t606 = t244 * t741 + t245 * t743 + t376 * t652 + t665;
t605 = t283 * t659 - t419 * t594 + t452 * t483 + t674;
t604 = t429 * t743 + t430 * t741 + t629;
t601 = t150 / 0.2e1 + t152 / 0.2e1 + t179 / 0.2e1 + t180 / 0.2e1;
t600 = t151 / 0.2e1 + t153 / 0.2e1 + t182 / 0.2e1 + t181 / 0.2e1;
t599 = t612 * t576;
t598 = t563 * t595 + t553 + t649;
t171 = t276 * rSges(6,1) + t275 * rSges(6,2) + t360 * rSges(6,3);
t593 = -t282 * t763 - t655;
t592 = -t430 * t763 - t654;
t464 = -rSges(3,1) * t595 + rSges(3,2) * t594 - rSges(3,3) * t741;
t518 = Icges(3,6) * t567 + (Icges(3,2) * t763 + t753) * t566;
t519 = Icges(3,5) * t567 + (Icges(3,1) * t572 + t656) * t566;
t525 = (Icges(3,5) * t763 - Icges(3,6) * t572) * t705;
t526 = (-Icges(3,2) * t572 + t656) * t705;
t527 = (Icges(3,1) * t763 - t753) * t705;
t591 = -t518 * t651 + t519 * t623 + t525 * t567 + t526 * t659 + t527 * t744;
t590 = t282 * t741 + t283 * t743 + t429 * t652 + t606;
t129 = -t360 * t446 + t361 * t447 - t387 * t594 + t388 * t588 + t389 * t498 + t445 * t483;
t135 = t400 * t454 + t401 * t455 - t425 * t594 + t426 * t505 - t427 * t602 + t453 * t483;
t584 = t111 / 0.2e1 + t101 / 0.2e1 + t135 / 0.2e1 + t129 / 0.2e1 + t608;
t128 = -t358 * t446 + t359 * t447 - t387 * t596 + t388 * t589 + t389 * t500 - t445 * t481;
t134 = t398 * t454 + t399 * t455 - t425 * t596 + t426 * t507 + t427 * t508 - t453 * t481;
t583 = t112 / 0.2e1 + t102 / 0.2e1 + t134 / 0.2e1 + t128 / 0.2e1 + t607;
t224 = -t445 * t594 + t446 * t588 + t447 * t498;
t256 = -t453 * t594 + t454 * t505 - t455 * t602;
t582 = t217 / 0.2e1 + t203 / 0.2e1 + t256 / 0.2e1 + t224 / 0.2e1 + t601;
t225 = -t445 * t596 + t446 * t589 + t447 * t500;
t257 = -t453 * t596 + t454 * t507 + t455 * t508;
t581 = -t257 / 0.2e1 - t225 / 0.2e1 - t218 / 0.2e1 - t204 / 0.2e1 - t600;
t579 = t586 + t622;
t578 = (-t565 + (-pkin(3) * t571 - pkin(8)) * t743) * qJD(1) - t484 * t563 - t660;
t577 = t578 + t751;
t528 = (rSges(3,1) * t763 - rSges(3,2) * t572) * t705;
t522 = t567 * rSges(3,3) + (rSges(3,1) * t572 + rSges(3,2) * t763) * t566;
t517 = Icges(3,3) * t567 + (Icges(3,5) * t572 + Icges(3,6) * t763) * t566;
t463 = Icges(3,1) * t542 + Icges(3,4) * t596 + Icges(3,5) * t743;
t462 = -Icges(3,1) * t595 + Icges(3,4) * t594 - Icges(3,5) * t741;
t461 = Icges(3,4) * t542 + Icges(3,2) * t596 + Icges(3,6) * t743;
t460 = -Icges(3,4) * t595 + Icges(3,2) * t594 - Icges(3,6) * t741;
t459 = Icges(3,5) * t542 + Icges(3,6) * t596 + Icges(3,3) * t743;
t458 = -Icges(3,5) * t595 + Icges(3,6) * t594 - Icges(3,3) * t741;
t444 = t465 + t708;
t443 = -t464 + t649;
t424 = -t464 * t567 - t522 * t741;
t423 = t465 * t567 - t522 * t743;
t391 = rSges(3,3) * t653 - t620;
t385 = Icges(3,1) * t484 - Icges(3,4) * t483 + Icges(3,5) * t653;
t384 = Icges(3,1) * t482 + Icges(3,4) * t481 + Icges(3,5) * t652;
t383 = Icges(3,4) * t484 - Icges(3,2) * t483 + Icges(3,6) * t653;
t382 = Icges(3,4) * t482 + Icges(3,2) * t481 + Icges(3,6) * t652;
t381 = Icges(3,5) * t484 - Icges(3,6) * t483 + Icges(3,3) * t653;
t380 = Icges(3,5) * t482 + Icges(3,6) * t481 + Icges(3,3) * t652;
t335 = (-t565 + (-rSges(3,3) - pkin(8)) * t743) * qJD(1) + t620;
t334 = t622 + t390;
t333 = t517 * t743 + t518 * t596 + t519 * t542;
t332 = -t517 * t741 + t518 * t594 - t519 * t595;
t324 = t708 - t716;
t323 = t649 + t717;
t320 = t591 * t567;
t319 = t567 * t390 + (-t522 * t706 - t528 * t573) * t566;
t318 = -t567 * t391 + (t522 * t707 - t528 * t576) * t566;
t317 = -t411 * t659 + t457 * t596;
t316 = t410 * t659 - t457 * t594;
t315 = t567 * t459 + (t461 * t763 + t463 * t572) * t566;
t314 = t567 * t458 + (t460 * t763 + t462 * t572) * t566;
t313 = t611 + t374;
t312 = -t594 * t757 + t598 + t618;
t293 = t459 * t743 + t461 * t596 + t463 * t542;
t292 = t458 * t743 + t460 * t596 + t462 * t542;
t291 = -t459 * t741 + t461 * t594 - t463 * t595;
t290 = -t458 * t741 + t460 * t594 - t462 * t595;
t289 = -t453 * t659 + t454 * t537 + t455 * t538;
t288 = t289 * t651;
t277 = -t410 * t596 + t411 * t594;
t272 = t566 * t638 + t567 * t717;
t271 = t411 * t567 + t710 * t743 + t473;
t255 = -t445 * t659 - t446 * t520 + t447 * t521;
t246 = t255 * t651;
t243 = rSges(5,3) * t483 - t619;
t227 = -t483 * t518 + t484 * t519 + t594 * t526 - t595 * t527 + (t517 * t707 - t525 * t576) * t566;
t226 = t481 * t518 + t482 * t519 + t596 * t526 + t542 * t527 + (t517 * t706 + t525 * t573) * t566;
t219 = (t410 * t573 + t411 * t576) * t566 + t709;
t216 = t430 + t611 + t311;
t215 = -t309 - t429 + t598 - t748;
t214 = -qJD(1) * t708 + t730;
t213 = t417 + t622 + t264;
t212 = t311 * t520 + t347 * t589;
t211 = -t309 * t520 - t347 * t588;
t210 = t611 + t777;
t209 = -t498 * t562 - (-pkin(5) * t570 + t569) * t594 + t598 + t783;
t208 = -t405 * t596 + t407 * t507 + t409 * t508;
t207 = -t404 * t596 + t406 * t507 + t408 * t508;
t206 = -t405 * t594 + t407 * t505 - t409 * t602;
t205 = -t404 * t594 + t406 * t505 - t408 * t602;
t202 = (-t374 * t763 - t654) * t566 + t713 * t596;
t201 = t373 * t659 - t448 * t594 + t724;
t200 = -t366 * t596 + t368 * t589 + t370 * t500;
t199 = -t365 * t596 + t367 * t589 + t369 * t500;
t198 = -t366 * t594 + t368 * t588 + t370 * t498;
t197 = -t365 * t594 + t367 * t588 + t369 * t498;
t196 = (-t487 + t722) * t567 + t566 * t621;
t195 = t374 * t567 + t664 * t743 + t723;
t194 = -t309 * t589 + t311 * t588;
t191 = t567 * t380 + (t763 * t382 + t384 * t572 + (-t461 * t572 + t463 * t763) * qJD(2)) * t566;
t190 = t567 * t381 + (t763 * t383 + t385 * t572 + (-t460 * t572 + t462 * t763) * qJD(2)) * t566;
t187 = t483 * t757 + t578 + t619;
t186 = t579 + t242;
t183 = -t373 * t596 - t594 * t721 - t331;
t177 = (t373 * t573 + t374 * t576) * t566 + t629;
t175 = t567 * t264 + t403 + (qJD(1) * t638 + t573 * t714) * t566;
t174 = t512 + t730 * t567 + (t457 * t707 + t576 * t714) * t566;
t173 = -t594 * t428 + t483 * t457 + (t265 * t763 - t410 * t704) * t566;
t172 = t596 * t428 + t481 * t457 + (-t264 * t763 + t411 * t704) * t566;
t147 = (-t311 * t763 + t592) * t566 + t668 * t596;
t146 = t309 * t659 - t347 * t594 + t628;
t145 = (-t487 + t670) * t567 + t566 * t609;
t144 = t311 * t567 + t630 * t743 + t667;
t133 = t520 * t728 + t589 * t727;
t132 = -t520 * t729 - t588 * t727;
t131 = t264 * t594 - t265 * t596 - t410 * t481 - t411 * t483;
t130 = -t309 * t596 - t594 * t669 + t726;
t127 = (t309 * t573 + t311 * t576) * t566 + t604;
t126 = t588 * t728 - t589 * t729;
t125 = -t171 - t283 + t577;
t124 = t282 + t579 + t169;
t123 = (t264 * t576 + t265 * t573 + (t410 * t576 + t573 * t716) * qJD(1)) * t566 + t665;
t122 = t257 * t567 + (-t207 * t576 + t208 * t573) * t566;
t121 = t256 * t567 + (-t205 * t576 + t206 * t573) * t566;
t120 = (-t728 * t763 + t592) * t566 + t631 * t596;
t119 = -t594 * t727 + t659 * t729 + t628;
t118 = t567 * t242 + (qJD(1) * t621 + t573 * t666) * t566 + t733;
t117 = (-t243 + t731) * t567 + (t448 * t707 + t576 * t666) * t566 + t712;
t116 = (-t487 + t633) * t567 + t566 * t599;
t115 = t567 * t728 + t612 * t743 + t667;
t114 = -t207 * t594 - t208 * t596 - t257 * t659;
t113 = -t205 * t594 - t206 * t596 - t256 * t659;
t110 = t225 * t567 + (-t199 * t576 + t200 * t573) * t566;
t109 = t224 * t567 + (-t197 * t576 + t198 * t573) * t566;
t108 = -t361 * t562 + t577 - t784;
t107 = t579 + t778;
t106 = -t199 * t594 - t200 * t596 - t225 * t659;
t105 = -t197 * t594 - t198 * t596 - t224 * t659;
t104 = -t594 * t392 + t483 * t448 + (t243 * t763 + t704 * t722) * t566 + t674;
t103 = t343 - t718 * t596 + t713 * t481 + (-t242 * t763 + t374 * t704 - t655) * t566;
t100 = -t594 * t632 - t596 * t729 + t726;
t99 = (t573 * t729 + t576 * t728) * t566 + t604;
t98 = -t171 * t520 - t254 * t588 - t309 * t485 + t347 * t360;
t97 = t169 * t520 + t254 * t589 + t311 * t485 - t347 * t358;
t96 = -t258 * t594 + t260 * t505 - t262 * t602 + t400 * t407 + t401 * t409 + t405 * t483;
t95 = -t259 * t594 + t261 * t505 - t263 * t602 + t400 * t406 + t401 * t408 + t404 * t483;
t94 = -t258 * t596 + t260 * t507 + t262 * t508 + t398 * t407 + t399 * t409 - t405 * t481;
t93 = -t259 * t596 + t261 * t507 + t263 * t508 + t398 * t406 + t399 * t408 - t404 * t481;
t88 = -t236 * t594 + t238 * t588 + t240 * t498 - t360 * t368 + t361 * t370 + t366 * t483;
t87 = -t237 * t594 + t239 * t588 + t241 * t498 - t360 * t367 + t361 * t369 + t365 * t483;
t86 = -t236 * t596 + t238 * t589 + t240 * t500 - t358 * t368 + t359 * t370 - t366 * t481;
t85 = -t237 * t596 + t239 * t589 + t241 * t500 - t358 * t367 + t359 * t369 - t365 * t481;
t82 = -t243 * t596 - t373 * t481 - (-t242 - t244) * t594 + t721 * t483 + t735;
t79 = -t152 * t594 - t153 * t596 - t193 * t659;
t78 = -t150 * t594 - t151 * t596 - t192 * t659;
t77 = t169 * t588 - t171 * t589 + t309 * t358 - t311 * t360;
t76 = (t242 * t576 + t243 * t573 + (t373 * t576 + (-t488 + t721) * t573) * qJD(1)) * t566 + t606;
t71 = -t152 * t588 - t153 * t589 + t193 * t520;
t70 = -t150 * t588 - t151 * t589 + t192 * t520;
t69 = t567 * t169 + (qJD(1) * t609 + t573 * t634) * t566 + t673;
t68 = (-t171 + t672) * t567 + (t347 * t707 + t576 * t634) * t566 + t663;
t55 = -t594 * t254 + t483 * t347 + (t171 * t763 + t670 * t704) * t566 + t605;
t54 = -t671 * t596 + t668 * t481 + (-t169 * t763 + t311 * t704 + t593) * t566 + t725;
t49 = t360 * t727 - t485 * t729 - t520 * t736 - t588 * t734;
t48 = -t358 * t727 + t485 * t728 + t520 * t737 + t589 * t734;
t47 = t737 * t567 + (qJD(1) * t599 + t573 * t613) * t566 + t673;
t46 = (t672 - t736) * t567 + (t576 * t613 + t707 * t727) * t566 + t663;
t37 = -t171 * t596 - t309 * t481 - (-t169 + t732) * t594 + t669 * t483 + t636;
t36 = (t169 * t576 + t171 * t573 + (t309 * t576 + (-t488 + t669) * t573) * qJD(1)) * t566 + t590;
t34 = -t734 * t594 + t727 * t483 + (t633 * t704 + t736 * t763) * t566 + t605;
t33 = -t635 * t596 + t631 * t481 + (t704 * t728 - t737 * t763 + t593) * t566 + t725;
t32 = t358 * t729 - t360 * t728 + t588 * t737 - t589 * t736;
t31 = -t111 * t594 - t112 * t596 - t178 * t659 + t217 * t483 - t218 * t481 + t288;
t29 = t135 * t567 + (t573 * t96 - t576 * t95 + (t205 * t573 + t206 * t576) * qJD(1)) * t566;
t28 = t134 * t567 + (t573 * t94 - t576 * t93 + (t207 * t573 + t208 * t576) * qJD(1)) * t566;
t27 = -t101 * t594 - t102 * t596 - t149 * t659 + t203 * t483 - t204 * t481 + t246;
t26 = t205 * t483 - t206 * t481 - t95 * t594 - t96 * t596 + (-t135 * t763 + t256 * t704) * t566;
t25 = t207 * t483 - t208 * t481 - t93 * t594 - t94 * t596 + (-t134 * t763 + t257 * t704) * t566;
t24 = t129 * t567 + (t573 * t88 - t576 * t87 + (t197 * t573 + t198 * t576) * qJD(1)) * t566;
t23 = t128 * t567 + (t573 * t86 - t576 * t85 + (t199 * t573 + t200 * t576) * qJD(1)) * t566;
t22 = (t737 * t576 + t736 * t573 + (t729 * t576 + (-t488 + t632) * t573) * qJD(1)) * t566 + t590;
t21 = -t736 * t596 - t729 * t481 - (t732 - t737) * t594 + t632 * t483 + t636;
t20 = t197 * t483 - t198 * t481 - t87 * t594 - t88 * t596 + (-t129 * t763 + t224 * t704) * t566;
t19 = t199 * t483 - t200 * t481 - t85 * t594 - t86 * t596 + (-t128 * t763 + t225 * t704) * t566;
t30 = [(t107 * t210 + t108 * t209) * t640 + (t124 * t216 + t125 * t215) * t642 + (t186 * t313 + t187 * t312) * t644 + (t213 * t324 + t214 * t323) * t646 + (t334 * t444 + t335 * t443) * t648 + t591 - t769; t320 + t90 + t89 + m(3) * (t318 * t443 + t319 * t444 + t334 * t423 + t335 * t424) + (t174 * t323 + t175 * t324 + t213 * t271 + t214 * t272) * m(4) + (t117 * t312 + t118 * t313 + t186 * t195 + t187 * t196) * m(5) + (t124 * t144 + t125 * t145 + t215 * t68 + t216 * t69) * m(6) + (t107 * t115 + t108 * t116 + t209 * t46 + t210 * t47) * m(7) + ((-t190 / 0.2e1 - t227 / 0.2e1 - t584) * t576 + (t191 / 0.2e1 + t226 / 0.2e1 + t583) * t573 + ((t315 / 0.2e1 + t333 / 0.2e1 - t581) * t576 + (t314 / 0.2e1 + t332 / 0.2e1 + t582) * t573) * qJD(1)) * t566 + t779; (t123 * t219 + t174 * t272 + t175 * t271) * t646 + (t117 * t196 + t118 * t195 + t177 * t76) * t644 + (t127 * t36 + t144 * t69 + t145 * t68) * t642 + (t115 * t47 + t116 * t46 + t22 * t99) * t640 + (t424 * t318 + t423 * t319 + (t464 * t573 + t465 * t576) * (t390 * t576 + t391 * t573 + (t464 * t576 - t465 * t573) * qJD(1)) * t767) * t648 + (((t382 * t596 + t542 * t384 + t481 * t461 + t482 * t463) * t573 + t293 * t706 - (t383 * t596 + t542 * t385 + t481 * t460 + t482 * t462) * t576 + t292 * t707 + ((t380 * t573 + t459 * t706) * t573 - (t381 * t573 + t458 * t706) * t576) * t566) * t566 + t28 + t23 + t10 + t9) * t743 + (-((t382 * t594 - t384 * t595 - t483 * t461 + t484 * t463) * t573 + t291 * t706 - (t383 * t594 - t385 * t595 - t483 * t460 + t484 * t462) * t576 + t290 * t707 + ((-t380 * t576 + t459 * t707) * t573 - (-t381 * t576 + t458 * t707) * t576) * t566) * t566 - t29 - t24 - t12 - t11) * t741 + ((-t290 * t576 + t291 * t573) * t566 + t121 + t109 + t65 + t64) * t653 + ((-t292 * t576 + t293 * t573) * t566 + t122 + t110 + t67 + t66) * t652 + (-t227 * t741 + t320 + (-t190 * t576 + t191 * t573 + (t314 * t573 + t315 * t576) * qJD(1)) * t566 + t226 * t743 + t333 * t652 + t332 * t653 + t770) * t567; t246 + t188 + (t172 * t324 + t173 * t323 + t213 * t317 + t214 * t316) * m(4) + (t103 * t313 + t104 * t312 + t186 * t202 + t187 * t201) * m(5) + (t124 * t147 + t125 * t146 + t215 * t55 + t216 * t54) * m(6) + (t107 * t120 + t108 * t119 + t209 * t34 + t210 * t33) * m(7) + t288 + t189 - t583 * t596 - t584 * t594 + t582 * t483 + t581 * t481 + t769 * t659; (t123 * t277 + t131 * t219 + t172 * t271 + t173 * t272 + t174 * t316 + t175 * t317) * m(4) + (t103 * t195 + t104 * t196 + t117 * t201 + t118 * t202 + t177 * t82 + t183 * t76) * m(5) + (t127 * t37 + t130 * t36 + t144 * t54 + t145 * t55 + t146 * t68 + t147 * t69) * m(6) + (t100 * t22 + t115 * t33 + t116 * t34 + t119 * t46 + t120 * t47 + t21 * t99) * m(7) + (t31 / 0.2e1 + t27 / 0.2e1 + t685) * t567 - (t28 / 0.2e1 + t23 / 0.2e1 + t694) * t596 - (t29 / 0.2e1 + t24 / 0.2e1 + t686) * t594 + (t109 / 0.2e1 + t121 / 0.2e1 + t680) * t483 + (-t110 / 0.2e1 - t122 / 0.2e1 - t679) * t481 + ((-t26 / 0.2e1 - t20 / 0.2e1 - t696) * t576 + (t19 / 0.2e1 + t25 / 0.2e1 + t697) * t573 + ((t771 * t573 / 0.2e1 - t772 * t576 / 0.2e1) * t566 + (t255 / 0.2e1 + t289 / 0.2e1) * t567 + t678) * t704 + ((t106 / 0.2e1 + t114 / 0.2e1 + t681) * t576 + (t105 / 0.2e1 + t113 / 0.2e1 + t682) * t573) * qJD(1) + t770 * t695) * t566; (t131 * t277 + t172 * t317 + t173 * t316) * t646 + (t103 * t202 + t104 * t201 + t183 * t82) * t644 + (t130 * t37 + t146 * t55 + t147 * t54) * t642 + (t100 * t21 + t119 * t34 + t120 * t33) * t640 - (t25 + t19 + t6 + t5) * t596 - (t26 + t20 + t8 + t7) * t594 + (t61 + t60 + t105 + t113) * t483 + (-t63 - t62 - t106 - t114) * t481 + ((t78 + t79 - t771 * t596 - t772 * t594 + (-t255 - t289) * t659) * t704 + (-t15 - t16 - t27 - t31) * t763) * t566; -(m(5) * t187 + m(6) * t125 + m(7) * t108) * t596 - (m(5) * t186 + m(6) * t124 + m(7) * t107) * t594 + (m(5) * t313 + m(6) * t216 + m(7) * t210) * t483 + (-m(5) * t312 - m(6) * t215 - m(7) * t209) * t481; (-t22 * t691 - t36 * t692 - t693 * t76) * t566 - (m(5) * t117 + m(6) * t68 + m(7) * t46) * t596 - (m(5) * t118 + m(6) * t69 + m(7) * t47) * t594 + (m(5) * t195 + m(6) * t144 + m(7) * t115) * t483 + (-m(5) * t196 - m(6) * t145 - m(7) * t116) * t481 + (m(5) * t177 + m(6) * t127 + m(7) * t99) * t651; (-t21 * t691 - t37 * t692 - t693 * t82) * t566 - (m(5) * t104 + m(6) * t55 + m(7) * t34) * t596 - (m(5) * t103 + m(6) * t54 + m(7) * t33) * t594 + (m(5) * t202 + m(6) * t147 + m(7) * t120) * t483 + (-m(5) * t201 - m(6) * t146 - m(7) * t119) * t481 + (m(5) * t183 + m(6) * t130 + m(7) * t100) * t651; 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + t766) * (-t572 * t650 * t767 + t481 * t596 - t483 * t594); (t107 * t133 + t108 * t132 + t209 * t49 + t210 * t48) * m(7) + (t124 * t212 + t125 * t211 + t215 * t98 + t216 * t97) * m(6) - t607 * t589 - t608 * t588 + t601 * t360 + t600 * t358 + t755 + t756; (t115 * t48 + t116 * t49 + t126 * t22 + t132 * t46 + t133 * t47 + t32 * t99) * m(7) + (t127 * t77 + t144 * t97 + t145 * t98 + t194 * t36 + t211 * t68 + t212 * t69) * m(6) + (t13 / 0.2e1 + t14 / 0.2e1) * t567 + (t17 / 0.2e1 + t18 / 0.2e1) * t520 - t694 * t589 - t686 * t588 + t678 * t485 + t680 * t360 + t679 * t358 + (-t698 * t576 + t699 * t573 + (t573 * t684 + t576 * t683) * qJD(1)) * t566; (t130 * t77 + t146 * t98 + t147 * t97 + t194 * t37 + t211 * t55 + t212 * t54) * m(6) + (t100 * t32 + t119 * t49 + t120 * t48 + t126 * t21 + t132 * t34 + t133 * t33) * m(7) - t699 * t596 - t698 * t594 + t685 * t520 - t697 * t589 - t696 * t588 + (t78 / 0.2e1 + t79 / 0.2e1) * t485 + t684 * t483 - t683 * t481 + t682 * t360 + t681 * t358 + ((t70 / 0.2e1 + t71 / 0.2e1) * t704 + t773 * t695) * t566; (-t32 * t691 - t692 * t77) * t566 - (m(6) * t98 + m(7) * t49) * t596 - (m(6) * t97 + m(7) * t48) * t594 + (m(6) * t212 + m(7) * t133) * t483 + (-m(6) * t211 - m(7) * t132) * t481 + (m(6) * t194 + m(7) * t126) * t651; (t126 * t32 + t132 * t49 + t133 * t48) * t640 + (t194 * t77 + t211 * t98 + t212 * t97) * t642 + t773 * t520 - (t2 + t1) * t589 - (t3 + t4) * t588 + (t70 + t71) * t485 + (t56 + t57) * t360 + (t59 + t58) * t358; (-t107 * t588 - t108 * t589 + t209 * t358 + t210 * t360) * m(7); (t115 * t360 + t116 * t358 + t22 * t520 - t46 * t589 - t47 * t588 + t485 * t99) * m(7); (t100 * t485 + t119 * t358 + t120 * t360 + t21 * t520 - t33 * t588 - t34 * t589) * m(7); 0.2e1 * (-t358 * t596 - t360 * t594 + t589 * t481 - t588 * t483 + (-t485 * t763 + t520 * t704) * t566) * t766; (t126 * t485 + t132 * t358 + t133 * t360 + t32 * t520 - t48 * t588 - t49 * t589) * m(7); (-t358 * t589 - t360 * t588 + t485 * t520) * t640;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t30(1) t30(2) t30(4) t30(7) t30(11) t30(16); t30(2) t30(3) t30(5) t30(8) t30(12) t30(17); t30(4) t30(5) t30(6) t30(9) t30(13) t30(18); t30(7) t30(8) t30(9) t30(10) t30(14) t30(19); t30(11) t30(12) t30(13) t30(14) t30(15) t30(20); t30(16) t30(17) t30(18) t30(19) t30(20) t30(21);];
Mq  = res;