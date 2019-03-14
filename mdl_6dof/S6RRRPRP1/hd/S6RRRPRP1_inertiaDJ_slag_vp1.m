% Calculate time derivative of joint inertia matrix for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:45
% EndTime: 2019-03-09 16:31:37
% DurationCPUTime: 35.17s
% Computational Cost: add. (53493->1094), mult. (52386->1448), div. (0->0), fcn. (48716->10), ass. (0->552)
t459 = qJ(2) + qJ(3);
t445 = pkin(10) + t459;
t439 = cos(t445);
t464 = cos(qJ(5));
t466 = cos(qJ(1));
t672 = t464 * t466;
t461 = sin(qJ(5));
t463 = sin(qJ(1));
t675 = t461 * t463;
t388 = -t439 * t675 - t672;
t673 = t463 * t464;
t674 = t461 * t466;
t389 = t439 * t673 - t674;
t438 = sin(t445);
t687 = t438 * t463;
t273 = Icges(7,5) * t389 + Icges(7,6) * t388 + Icges(7,3) * t687;
t275 = Icges(6,5) * t389 + Icges(6,6) * t388 + Icges(6,3) * t687;
t822 = t273 + t275;
t390 = -t439 * t674 + t673;
t391 = t439 * t672 + t675;
t685 = t438 * t466;
t274 = Icges(7,5) * t391 + Icges(7,6) * t390 + Icges(7,3) * t685;
t276 = Icges(6,5) * t391 + Icges(6,6) * t390 + Icges(6,3) * t685;
t821 = t274 + t276;
t277 = Icges(7,4) * t389 + Icges(7,2) * t388 + Icges(7,6) * t687;
t279 = Icges(6,4) * t389 + Icges(6,2) * t388 + Icges(6,6) * t687;
t820 = t277 + t279;
t278 = Icges(7,4) * t391 + Icges(7,2) * t390 + Icges(7,6) * t685;
t280 = Icges(6,4) * t391 + Icges(6,2) * t390 + Icges(6,6) * t685;
t819 = t278 + t280;
t281 = Icges(7,1) * t389 + Icges(7,4) * t388 + Icges(7,5) * t687;
t283 = Icges(6,1) * t389 + Icges(6,4) * t388 + Icges(6,5) * t687;
t818 = t281 + t283;
t282 = Icges(7,1) * t391 + Icges(7,4) * t390 + Icges(7,5) * t685;
t284 = Icges(6,1) * t391 + Icges(6,4) * t390 + Icges(6,5) * t685;
t817 = t282 + t284;
t816 = t820 * t390 + t818 * t391 + t685 * t822;
t815 = t390 * t819 + t391 * t817 + t685 * t821;
t460 = -qJ(6) - pkin(9);
t814 = rSges(7,3) - t460;
t586 = -qJD(5) * t439 + qJD(1);
t517 = t464 * t586;
t585 = qJD(1) * t439 - qJD(5);
t456 = qJD(2) + qJD(3);
t678 = t456 * t463;
t629 = t438 * t678;
t252 = t463 * t517 + (-t466 * t585 + t629) * t461;
t516 = t586 * t461;
t677 = t456 * t464;
t253 = t585 * t672 + (-t438 * t677 + t516) * t463;
t642 = qJD(1) * t466;
t487 = t438 * t642 + t439 * t678;
t148 = Icges(7,5) * t253 + Icges(7,6) * t252 + Icges(7,3) * t487;
t150 = Icges(6,5) * t253 + Icges(6,6) * t252 + Icges(6,3) * t487;
t152 = Icges(7,4) * t253 + Icges(7,2) * t252 + Icges(7,6) * t487;
t154 = Icges(6,4) * t253 + Icges(6,2) * t252 + Icges(6,6) * t487;
t156 = Icges(7,1) * t253 + Icges(7,4) * t252 + Icges(7,5) * t487;
t158 = Icges(6,1) * t253 + Icges(6,4) * t252 + Icges(6,5) * t487;
t676 = t456 * t466;
t628 = t438 * t676;
t755 = t463 * t585 + t628;
t250 = t461 * t755 + t466 * t517;
t251 = -t464 * t755 + t466 * t516;
t643 = qJD(1) * t463;
t615 = t438 * t643;
t625 = t439 * t676;
t486 = -t615 + t625;
t813 = (t148 + t150) * t685 + t822 * t486 + (t156 + t158) * t391 + (t152 + t154) * t390 + t818 * t251 + t820 * t250;
t147 = Icges(7,5) * t251 + Icges(7,6) * t250 + Icges(7,3) * t486;
t149 = Icges(6,5) * t251 + Icges(6,6) * t250 + Icges(6,3) * t486;
t151 = Icges(7,4) * t251 + Icges(7,2) * t250 + Icges(7,6) * t486;
t153 = Icges(6,4) * t251 + Icges(6,2) * t250 + Icges(6,6) * t486;
t155 = Icges(7,1) * t251 + Icges(7,4) * t250 + Icges(7,5) * t486;
t157 = Icges(6,1) * t251 + Icges(6,4) * t250 + Icges(6,5) * t486;
t812 = (t147 + t149) * t685 + t821 * t486 + (t155 + t157) * t391 + (t151 + t153) * t390 + t817 * t251 + t819 * t250;
t547 = Icges(7,5) * t464 - Icges(7,6) * t461;
t684 = t439 * t456;
t221 = t547 * t684 + (Icges(7,3) * t456 + (-Icges(7,5) * t461 - Icges(7,6) * t464) * qJD(5)) * t438;
t708 = Icges(7,4) * t464;
t552 = -Icges(7,2) * t461 + t708;
t709 = Icges(7,4) * t461;
t223 = t552 * t684 + (Icges(7,6) * t456 + (-Icges(7,2) * t464 - t709) * qJD(5)) * t438;
t558 = Icges(7,1) * t464 - t709;
t225 = t558 * t684 + (Icges(7,5) * t456 + (-Icges(7,1) * t461 - t708) * qJD(5)) * t438;
t319 = -Icges(7,3) * t439 + t438 * t547;
t321 = -Icges(7,6) * t439 + t438 * t552;
t323 = -Icges(7,5) * t439 + t438 * t558;
t62 = t221 * t685 + t223 * t390 + t225 * t391 + t250 * t321 + t251 * t323 + t319 * t486;
t548 = Icges(6,5) * t464 - Icges(6,6) * t461;
t222 = t548 * t684 + (Icges(6,3) * t456 + (-Icges(6,5) * t461 - Icges(6,6) * t464) * qJD(5)) * t438;
t710 = Icges(6,4) * t464;
t553 = -Icges(6,2) * t461 + t710;
t711 = Icges(6,4) * t461;
t224 = t553 * t684 + (Icges(6,6) * t456 + (-Icges(6,2) * t464 - t711) * qJD(5)) * t438;
t559 = Icges(6,1) * t464 - t711;
t226 = t559 * t684 + (Icges(6,5) * t456 + (-Icges(6,1) * t461 - t710) * qJD(5)) * t438;
t320 = -Icges(6,3) * t439 + t438 * t548;
t322 = -Icges(6,6) * t439 + t438 * t553;
t324 = -Icges(6,5) * t439 + t438 * t559;
t63 = t222 * t685 + t224 * t390 + t226 * t391 + t250 * t322 + t251 * t324 + t320 * t486;
t811 = -t62 - t63;
t64 = t221 * t687 + t223 * t388 + t225 * t389 + t252 * t321 + t253 * t323 + t319 * t487;
t65 = t222 * t687 + t224 * t388 + t226 * t389 + t252 * t322 + t253 * t324 + t320 * t487;
t810 = -t64 - t65;
t168 = t319 * t687 + t321 * t388 + t323 * t389;
t169 = t320 * t687 + t322 * t388 + t324 * t389;
t809 = t168 + t169;
t170 = t319 * t685 + t321 * t390 + t323 * t391;
t171 = t320 * t685 + t322 * t390 + t324 * t391;
t808 = t170 + t171;
t807 = t319 + t320;
t779 = t322 + t321;
t800 = t323 + t324;
t549 = Icges(5,5) * t439 - Icges(5,6) * t438;
t338 = Icges(5,3) * t463 + t466 * t549;
t448 = sin(t459);
t449 = cos(t459);
t550 = Icges(4,5) * t449 - Icges(4,6) * t448;
t358 = Icges(4,3) * t463 + t466 * t550;
t806 = t338 + t358;
t392 = Icges(5,5) * t438 + Icges(5,6) * t439;
t404 = Icges(4,5) * t448 + Icges(4,6) * t449;
t777 = t392 + t404;
t713 = Icges(5,4) * t438;
t393 = Icges(5,2) * t439 + t713;
t712 = Icges(5,4) * t439;
t394 = Icges(5,1) * t438 + t712;
t715 = Icges(4,4) * t448;
t405 = Icges(4,2) * t449 + t715;
t714 = Icges(4,4) * t449;
t406 = Icges(4,1) * t448 + t714;
t805 = t393 * t438 - t394 * t439 + t405 * t448 - t406 * t449;
t804 = t463 * t816 + t815 * t466;
t774 = t815 * t463 - t466 * t816;
t113 = t275 * t687 + t279 * t388 + t283 * t389;
t114 = t276 * t687 + t280 * t388 + t284 * t389;
t542 = t113 * t463 + t114 * t466;
t111 = t273 * t687 + t277 * t388 + t281 * t389;
t112 = t274 * t687 + t278 * t388 + t282 * t389;
t544 = t111 * t463 + t112 * t466;
t803 = t542 + t544;
t773 = (-t111 - t113) * t466 + (t112 + t114) * t463;
t337 = -Icges(5,3) * t466 + t463 * t549;
t357 = -Icges(4,3) * t466 + t463 * t550;
t554 = -Icges(5,2) * t438 + t712;
t339 = -Icges(5,6) * t466 + t463 * t554;
t560 = Icges(5,1) * t439 - t713;
t341 = -Icges(5,5) * t466 + t463 * t560;
t527 = t339 * t438 - t341 * t439;
t762 = t466 * t527;
t555 = -Icges(4,2) * t448 + t714;
t359 = -Icges(4,6) * t466 + t463 * t555;
t561 = Icges(4,1) * t449 - t715;
t361 = -Icges(4,5) * t466 + t463 * t561;
t525 = t359 * t448 - t361 * t449;
t763 = t466 * t525;
t802 = t762 + t763 + (-t337 - t357) * t463;
t440 = pkin(5) * t464 + pkin(4);
t683 = t439 * t466;
t801 = t391 * rSges(7,1) + t390 * rSges(7,2) + pkin(5) * t675 + t440 * t683 + t685 * t814;
t451 = t463 * rSges(4,3);
t728 = rSges(4,1) * t449;
t571 = -rSges(4,2) * t448 + t728;
t366 = t466 * t571 + t451;
t637 = qJD(5) * t464;
t630 = pkin(5) * t637;
t738 = pkin(5) * t461;
t632 = qJD(1) * t738;
t799 = t251 * rSges(7,1) + t250 * rSges(7,2) + rSges(7,3) * t625 + qJD(6) * t685 + t460 * t615 + t463 * t630 + t466 * t632;
t533 = -t278 * t461 + t282 * t464;
t133 = -t274 * t439 + t438 * t533;
t531 = -t280 * t461 + t284 * t464;
t135 = -t276 * t439 + t438 * t531;
t669 = t133 + t135;
t534 = -t277 * t461 + t281 * t464;
t132 = -t273 * t439 + t438 * t534;
t532 = -t279 * t461 + t283 * t464;
t134 = -t275 * t439 + t438 * t532;
t670 = t132 + t134;
t798 = t463 * t669 - t466 * t670;
t797 = t463 * t670 + t466 * t669;
t421 = pkin(4) * t683;
t368 = pkin(9) * t685 + t421;
t663 = -t368 + t801;
t733 = pkin(9) + t460;
t735 = pkin(4) - t440;
t483 = -t438 * t733 - t439 * t735;
t565 = -rSges(7,1) * t389 - rSges(7,2) * t388;
t664 = rSges(7,3) * t687 - pkin(5) * t674 + t463 * t483 - t565;
t796 = -t463 * t664 - t466 * t663;
t795 = (t456 * t804 + t811) * t439 + (-t774 * qJD(1) + t808 * t456 + t463 * t813 + t812 * t466) * t438;
t33 = t148 * t687 + t152 * t388 + t156 * t389 + t252 * t277 + t253 * t281 + t273 * t487;
t34 = t147 * t687 + t151 * t388 + t155 * t389 + t252 * t278 + t253 * t282 + t274 * t487;
t35 = t150 * t687 + t154 * t388 + t158 * t389 + t252 * t279 + t253 * t283 + t275 * t487;
t36 = t149 * t687 + t153 * t388 + t157 * t389 + t252 * t280 + t253 * t284 + t276 * t487;
t794 = (t456 * t803 + t810) * t439 + ((t34 + t36) * t466 + (t33 + t35) * t463 + t809 * t456 - t773 * qJD(1)) * t438;
t793 = t804 * qJD(1) + t812 * t463 - t466 * t813;
t17 = qJD(1) * t544 - t33 * t466 + t34 * t463;
t18 = qJD(1) * t542 - t35 * t466 + t36 * t463;
t792 = t17 + t18;
t732 = t438 * t803 - t439 * t809;
t791 = t438 * t804 - t439 * t808;
t403 = pkin(9) * t625;
t638 = qJD(5) * t461;
t631 = pkin(5) * t638;
t507 = -t456 * t460 - t631;
t783 = -t403 + (pkin(9) * t643 + t676 * t735) * t438 + (t466 * t507 + t643 * t735) * t439 - rSges(7,3) * t615 + t799;
t402 = pkin(4) * t629;
t566 = rSges(7,1) * t253 + rSges(7,2) * t252;
t682 = t440 * t456;
t756 = -(-qJD(6) + t682) * t438 + t632;
t667 = t402 + (qJD(1) * t483 - t630) * t466 + ((-t456 * t733 - t631) * t439 + t756) * t463 + rSges(7,3) * t487 + t566;
t782 = t807 * t439 + (t461 * t779 - t464 * t800) * t438;
t360 = Icges(4,6) * t463 + t466 * t555;
t362 = Icges(4,5) * t463 + t466 * t561;
t524 = t360 * t448 - t362 * t449;
t340 = Icges(5,6) * t463 + t466 * t554;
t342 = Icges(5,5) * t463 + t466 * t560;
t526 = t340 * t438 - t342 * t439;
t781 = (-t524 - t526) * t466 + t806 * t463;
t780 = -t221 - t222;
t564 = rSges(7,1) * t464 - rSges(7,2) * t461;
t661 = (-rSges(7,3) + t733) * t439 + (t564 - t735) * t438;
t355 = t554 * t456;
t356 = t560 * t456;
t372 = t555 * t456;
t373 = t561 * t456;
t690 = t406 * t456;
t691 = t405 * t456;
t692 = t394 * t456;
t693 = t393 * t456;
t776 = (t373 - t691) * t449 + (-t372 - t690) * t448 + (t356 - t693) * t439 + (-t355 - t692) * t438 + t777 * qJD(1);
t467 = -pkin(8) - pkin(7);
t462 = sin(qJ(2));
t640 = qJD(2) * t462;
t634 = pkin(2) * t640;
t772 = qJD(1) * t467 + t634;
t771 = (t549 + t550) * t456 + t805 * qJD(1);
t770 = (-t224 - t223) * t461;
t395 = rSges(5,1) * t438 + rSges(5,2) * t439;
t503 = t395 * t456;
t465 = cos(qJ(2));
t426 = rSges(3,1) * t462 + rSges(3,2) * t465;
t501 = qJD(2) * t426;
t769 = t463 * t501;
t716 = Icges(3,4) * t465;
t557 = -Icges(3,2) * t462 + t716;
t382 = Icges(3,6) * t463 + t466 * t557;
t717 = Icges(3,4) * t462;
t563 = Icges(3,1) * t465 - t717;
t384 = Icges(3,5) * t463 + t466 * t563;
t522 = t382 * t462 - t384 * t465;
t768 = t463 * t522;
t767 = t463 * t524;
t766 = t463 * t526;
t441 = t465 * pkin(2) + pkin(1);
t736 = pkin(1) - t441;
t765 = t463 * t736;
t381 = -Icges(3,6) * t466 + t463 * t557;
t383 = -Icges(3,5) * t466 + t463 * t563;
t523 = t381 * t462 - t383 * t465;
t764 = t466 * t523;
t568 = -rSges(6,1) * t389 - rSges(6,2) * t388;
t286 = rSges(6,3) * t687 - t568;
t288 = t391 * rSges(6,1) + t390 * rSges(6,2) + rSges(6,3) * t685;
t761 = -t463 * t286 - t466 * t288;
t689 = t438 * t456;
t759 = t807 * t689 + (t225 + t226) * t438 * t464 + t800 * t439 * t677;
t758 = qJD(1) * t337;
t757 = qJD(1) * t357;
t551 = Icges(3,5) * t465 - Icges(3,6) * t462;
t379 = -Icges(3,3) * t466 + t463 * t551;
t455 = -qJ(4) + t467;
t648 = t455 - t467;
t413 = pkin(3) * t449 + t441;
t654 = t413 - t441;
t312 = t463 * t654 + t466 * t648;
t752 = 2 * m(3);
t751 = 2 * m(4);
t750 = 2 * m(5);
t749 = 2 * m(6);
t748 = 2 * m(7);
t457 = t463 ^ 2;
t458 = t466 ^ 2;
t746 = t463 / 0.2e1;
t745 = -t466 / 0.2e1;
t744 = -rSges(6,3) - pkin(9);
t743 = m(3) * t426;
t724 = rSges(4,2) * t449;
t408 = rSges(4,1) * t448 + t724;
t742 = m(4) * t408;
t741 = pkin(2) * t462;
t740 = pkin(3) * t448;
t739 = pkin(4) * t439;
t737 = t463 * pkin(7);
t454 = t466 * pkin(7);
t734 = -pkin(7) - t467;
t679 = t456 * t461;
t730 = t759 + (-t679 * t779 + t780) * t439 + (t770 + (-t461 * t800 - t464 * t779) * qJD(5)) * t438;
t729 = rSges(3,1) * t465;
t727 = rSges(5,1) * t439;
t726 = rSges(3,2) * t462;
t723 = rSges(3,3) * t466;
t42 = (t456 * t534 - t148) * t439 + (-t152 * t461 + t156 * t464 + t273 * t456 + (-t277 * t464 - t281 * t461) * qJD(5)) * t438;
t722 = t42 * t466;
t43 = (t456 * t533 - t147) * t439 + (-t151 * t461 + t155 * t464 + t274 * t456 + (-t278 * t464 - t282 * t461) * qJD(5)) * t438;
t721 = t43 * t463;
t44 = (t456 * t532 - t150) * t439 + (-t154 * t461 + t158 * t464 + t275 * t456 + (-t279 * t464 - t283 * t461) * qJD(5)) * t438;
t720 = t44 * t466;
t45 = (t456 * t531 - t149) * t439 + (-t153 * t461 + t157 * t464 + t276 * t456 + (-t280 * t464 - t284 * t461) * qJD(5)) * t438;
t719 = t45 * t463;
t452 = t463 * rSges(3,3);
t450 = t463 * rSges(5,3);
t378 = t571 * t456;
t698 = t378 * t463;
t697 = t381 * t465;
t696 = t382 * t465;
t695 = t383 * t462;
t694 = t384 * t462;
t681 = t448 * t456;
t680 = t449 * t456;
t671 = t466 * t455;
t666 = t782 * t689;
t612 = t438 * t638;
t665 = -pkin(5) * t612 - qJD(6) * t439 + t456 * t483 + t564 * t684 + (rSges(7,3) * t456 + (-rSges(7,1) * t461 - rSges(7,2) * t464) * qJD(5)) * t438;
t399 = t466 * t413;
t427 = t466 * t441;
t313 = -t463 * t648 + t399 - t427;
t662 = t463 * t312 + t466 * t313;
t346 = rSges(5,1) * t683 - rSges(5,2) * t685 + t450;
t660 = -t313 - t346;
t659 = -t313 - t368;
t352 = t466 * t467 + t454 - t765;
t353 = -t466 * pkin(1) + t463 * t734 + t427;
t658 = t463 * t352 + t466 * t353;
t365 = -t466 * rSges(4,3) + t463 * t571;
t266 = t463 * t365 + t466 * t366;
t614 = t448 * t643;
t422 = pkin(3) * t614;
t657 = t395 * t643 + t422;
t396 = pkin(4) * t438 - pkin(9) * t439;
t656 = t396 * t643 + t422;
t397 = -pkin(3) * t681 - t634;
t655 = qJD(4) * t463 + t466 * t397;
t653 = rSges(5,2) * t615 + rSges(5,3) * t642;
t652 = rSges(4,2) * t614 + rSges(4,3) * t642;
t651 = qJD(4) * t466 + t455 * t643;
t650 = t772 * t463;
t649 = t466 * t729 + t452;
t647 = t457 + t458;
t646 = qJD(1) * t338;
t645 = qJD(1) * t358;
t380 = Icges(3,3) * t463 + t466 * t551;
t644 = qJD(1) * t380;
t639 = qJD(2) * t465;
t636 = pkin(3) * t680;
t635 = t466 * t726;
t633 = pkin(2) * t639;
t587 = t466 * t634;
t624 = t463 * (t397 * t463 + t642 * t654 + t650 - t651) + t466 * (-qJD(1) * t312 + t587 + t655) + t312 * t642;
t620 = t251 * rSges(6,1) + t250 * rSges(6,2) + rSges(6,3) * t625;
t502 = t408 * t456;
t619 = t463 * (qJD(1) * t366 - t463 * t502) + t466 * (-t676 * t724 + (-t448 * t676 - t449 * t643) * rSges(4,1) + t652) + t365 * t642;
t618 = -t288 + t659;
t617 = t463 * ((-t466 * t736 - t737) * qJD(1) - t650) + t466 * (-t587 + (t466 * t734 + t765) * qJD(1)) + t352 * t642;
t567 = rSges(6,1) * t464 - rSges(6,2) * t461;
t326 = -rSges(6,3) * t439 + t438 * t567;
t317 = t326 * t643;
t616 = t317 + t656;
t613 = t462 * t643;
t606 = t643 / 0.2e1;
t605 = t642 / 0.2e1;
t604 = -t408 - t741;
t603 = -t395 - t740;
t602 = -t396 - t740;
t601 = t463 * t661;
t600 = t661 * t466;
t242 = -qJD(1) * t339 - t466 * t693;
t599 = t342 * t456 + t242;
t243 = qJD(1) * t340 - t463 * t693;
t598 = t341 * t456 + t243;
t244 = -qJD(1) * t341 - t466 * t692;
t597 = -t340 * t456 + t244;
t245 = qJD(1) * t342 - t463 * t692;
t596 = t339 * t456 - t245;
t269 = -qJD(1) * t359 - t466 * t691;
t595 = t362 * t456 + t269;
t270 = qJD(1) * t360 - t463 * t691;
t594 = t361 * t456 + t270;
t271 = -qJD(1) * t361 - t466 * t690;
t593 = -t360 * t456 + t271;
t272 = qJD(1) * t362 - t463 * t690;
t592 = t359 * t456 - t272;
t591 = -t463 * t455 + t399;
t590 = -t439 * t440 - t413;
t589 = qJD(1) * t661;
t584 = t659 - t663;
t583 = t643 * t661 + t656;
t570 = -rSges(5,2) * t438 + t727;
t345 = -rSges(5,3) * t466 + t463 * t570;
t167 = t463 * t345 + t466 * t346 + t662;
t573 = pkin(9) * t438 + t739;
t367 = t573 * t463;
t582 = t463 * t367 + t466 * t368 + t662;
t364 = t573 * t456;
t576 = -t364 - t636;
t575 = -t326 + t602;
t574 = -t740 - t741;
t572 = -t726 + t729;
t569 = rSges(6,1) * t253 + rSges(6,2) * t252;
t562 = Icges(3,1) * t462 + t716;
t556 = Icges(3,2) * t465 + t717;
t109 = t438 * t601 + t439 * t664;
t110 = -t438 * t600 - t439 * t663;
t546 = t109 * t466 + t110 * t463;
t513 = -t396 + t574;
t484 = t513 - t661;
t177 = t484 * t463;
t178 = t484 * t466;
t537 = t177 * t463 + t178 * t466;
t485 = -t438 * t814 + t590;
t179 = (-t455 + t738) * t466 + t485 * t463 + t565;
t180 = t591 + t801;
t536 = t179 * t466 + t180 * t463;
t518 = t602 - t661;
t184 = t518 * t463;
t185 = t518 * t466;
t535 = t184 * t463 + t185 * t466;
t530 = t286 * t466 - t288 * t463;
t230 = t567 * t684 + (rSges(6,3) * t456 + (-rSges(6,1) * t461 - rSges(6,2) * t464) * qJD(5)) * t438;
t519 = -t230 + t576;
t515 = -pkin(1) - t572;
t514 = -t395 + t574;
t232 = t575 * t466;
t511 = -t441 - t571;
t510 = -t413 - t570;
t488 = -t439 * t643 - t628;
t509 = t463 * (-t463 * t503 + (t466 * t570 + t450) * qJD(1)) + t466 * (rSges(5,1) * t488 - rSges(5,2) * t625 + t653) + t345 * t642 + t624;
t508 = t463 * (pkin(9) * t487 + qJD(1) * t421 - t402) + t466 * (pkin(4) * t488 - pkin(9) * t615 + t403) + t367 * t642 + t624;
t102 = t582 - t761;
t506 = t134 / 0.2e1 + t132 / 0.2e1 + t168 / 0.2e1 + t169 / 0.2e1;
t505 = t135 / 0.2e1 + t133 / 0.2e1 + t170 / 0.2e1 + t171 / 0.2e1;
t504 = t576 - t665;
t500 = -t326 + t513;
t495 = t456 * t404;
t494 = t456 * t392;
t493 = -t633 - t636;
t492 = qJD(2) * t562;
t491 = qJD(2) * t556;
t490 = qJD(2) * (-Icges(3,5) * t462 - Icges(3,6) * t465);
t489 = t438 * t744 - t413 - t739;
t315 = t514 * t466;
t363 = t570 * t456;
t482 = -t363 + t493;
t481 = -t364 + t493;
t90 = t582 - t796;
t218 = t500 * t466;
t480 = -t463 * t663 + t466 * t664;
t479 = -t230 + t481;
t160 = -rSges(6,3) * t615 + t620;
t162 = rSges(6,3) * t487 + t569;
t478 = t466 * t160 + t463 * t162 + t286 * t642 + t508;
t477 = t481 - t665;
t476 = rSges(3,2) * t613 + rSges(3,3) * t642 - t466 * t501;
t186 = -t337 * t466 - t463 * t527;
t187 = -t338 * t466 - t766;
t194 = -t357 * t466 - t463 * t525;
t195 = -t358 * t466 - t767;
t240 = -t466 * t494 - t758;
t241 = -t463 * t494 + t646;
t267 = -t466 * t495 - t757;
t268 = -t463 * t495 + t645;
t475 = t773 * t643 + t774 * t642 + ((-t186 - t194) * t643 + t802 * t642) * t466 + (((t240 + t267) * t463 + (t766 + t767 - t802) * qJD(1)) * t463 + (t187 + t195) * t643 + t781 * t642 + ((-t599 * t438 + t597 * t439 - t595 * t448 + t593 * t449 - t241 - t268) * t463 + (t243 * t438 - t245 * t439 + t270 * t448 - t272 * t449 + t339 * t684 + t341 * t689 + t359 * t680 + t361 * t681 - t757 - t758) * t466 + ((-t527 - t525 + t806) * t463 + t781) * qJD(1)) * t466 + t793) * t463;
t474 = t463 * t489 - t671;
t473 = t463 * t667 + t466 * t783 + t642 * t664 + t508;
t25 = (t466 * t241 + (t187 + t762) * qJD(1)) * t466 + (t186 * qJD(1) + (-t242 * t438 + t244 * t439 - t340 * t684 - t342 * t689 + t646) * t463 + (-t240 + t596 * t439 + t598 * t438 + (-t337 - t526) * qJD(1)) * t466) * t463;
t28 = (t466 * t268 + (t195 + t763) * qJD(1)) * t466 + (t194 * qJD(1) + (-t269 * t448 + t271 * t449 - t360 * t680 - t362 * t681 + t645) * t463 + (-t267 + t592 * t449 + t594 * t448 + (-t357 - t524) * qJD(1)) * t466) * t463;
t470 = (-t28 - t25 - t792) * t466 + t475;
t469 = -(qJD(1) * t797 + t719 - t720 + t721 - t722) * t439 / 0.2e1 + t795 * t746 + t794 * t745 + t798 * t689 / 0.2e1 + t792 * t687 / 0.2e1 + t793 * t685 / 0.2e1 + t732 * t606 - t774 * t615 / 0.2e1 + (t463 * t773 + t466 * t774) * t684 / 0.2e1 + (t438 * t773 + t791) * t605;
t468 = t719 / 0.2e1 - t720 / 0.2e1 + t721 / 0.2e1 - t722 / 0.2e1 + (t438 * t597 + t439 * t599 + t448 * t593 + t449 * t595 + t463 * t771 + t466 * t776 - t811) * t746 + (-t438 * t596 + t439 * t598 - t448 * t592 + t449 * t594 + t463 * t776 - t466 * t771 - t810) * t745 + (t339 * t439 + t341 * t438 + t359 * t449 + t361 * t448 - t463 * t805 - t466 * t777 + t670 + t809) * t606 + (t340 * t439 + t342 * t438 + t360 * t449 + t362 * t448 + t463 * t777 - t466 * t805 + t669 + t808) * t605;
t433 = pkin(2) * t613;
t412 = t572 * qJD(2);
t387 = -t635 + t649;
t386 = t463 * t572 - t723;
t351 = t604 * t466;
t350 = t604 * t463;
t334 = t737 + (pkin(1) - t726) * t466 + t649;
t333 = t463 * t515 + t454 + t723;
t330 = t603 * t466;
t329 = t603 * t463;
t314 = t514 * t463;
t311 = -t463 * t467 + t366 + t427;
t310 = (rSges(4,3) - t467) * t466 + t511 * t463;
t302 = t463 * t490 + t644;
t301 = -qJD(1) * t379 + t466 * t490;
t292 = t346 + t591;
t291 = (rSges(5,3) - t455) * t466 + t510 * t463;
t290 = t769 + ((-rSges(3,3) - pkin(7)) * t463 + t515 * t466) * qJD(1);
t289 = (t454 + (-pkin(1) - t729) * t463) * qJD(1) + t476;
t249 = -t408 * t642 - t698 + (-t462 * t642 - t463 * t639) * pkin(2);
t248 = t408 * t643 + t433 + (-t378 - t633) * t466;
t231 = t575 * t463;
t217 = t500 * t463;
t216 = -t395 * t642 - t363 * t463 + (-t448 * t642 - t449 * t678) * pkin(3);
t215 = (-t363 - t636) * t466 + t657;
t209 = t380 * t463 - t466 * t522;
t208 = t379 * t463 - t764;
t207 = -t380 * t466 - t768;
t206 = -t379 * t466 - t463 * t523;
t205 = t408 * t678 + (t466 * t511 - t451) * qJD(1) + t650;
t204 = (-t441 - t728) * t643 + (-t502 - t772) * t466 + t652;
t199 = qJD(1) * t315 + t463 * t482;
t198 = t466 * t482 + t433 + t657;
t193 = t591 + t288 + t368;
t192 = t474 + t568;
t191 = -t288 * t439 - t326 * t685;
t190 = t286 * t439 + t326 * t687;
t183 = (-t397 + t503) * t463 + (t466 * t510 - t450) * qJD(1) + t651;
t182 = -t466 * t503 + (-t671 + (-t413 - t727) * t463) * qJD(1) + t653 + t655;
t181 = t266 + t658;
t174 = t530 * t438;
t138 = -t366 * t643 + t619;
t137 = qJD(1) * t232 + t463 * t519;
t136 = t466 * t519 + t616;
t129 = qJD(1) * t218 + t463 * t479;
t128 = t466 * t479 + t433 + t616;
t121 = t167 + t658;
t101 = t402 + (t684 * t744 - t397) * t463 + t489 * t642 - t569 + t651;
t100 = -pkin(4) * t628 + qJD(1) * t474 + t403 + t620 + t655;
t99 = t480 * t438;
t98 = qJD(1) * t185 + t463 * t504;
t97 = t466 * t504 + t583;
t96 = t102 + t658;
t95 = (qJD(1) * t485 + t630) * t466 + (-t397 + (-t456 * t814 + t631) * t439 - t756) * t463 - t566 + t651;
t94 = (-t438 * t682 + t439 * t507) * t466 + (-t671 + (-rSges(7,3) * t438 + t590) * t463) * qJD(1) + t655 + t799;
t93 = (-t353 - t366) * t643 + t617 + t619;
t92 = qJD(1) * t178 + t463 * t477;
t91 = t466 * t477 + t433 + t583;
t89 = (t326 * t678 + t162) * t439 + (t230 * t463 - t286 * t456 + t326 * t642) * t438;
t88 = (-t326 * t676 - t160) * t439 + (-t230 * t466 + t288 * t456 + t317) * t438;
t85 = t90 + t658;
t76 = t643 * t660 + t509;
t57 = t530 * t684 + (qJD(1) * t761 - t160 * t463 + t162 * t466) * t438;
t56 = (-t353 + t660) * t643 + t509 + t617;
t47 = (t456 * t601 + t667) * t439 + (-t456 * t664 + t463 * t665 + t466 * t589) * t438;
t46 = (-t456 * t600 - t783) * t439 + (t456 * t663 + t463 * t589 - t466 * t665) * t438;
t37 = t618 * t643 + t478;
t26 = (-t353 + t618) * t643 + t478 + t617;
t23 = t480 * t684 + (qJD(1) * t796 - t463 * t783 + t667 * t466) * t438;
t22 = t584 * t643 + t473;
t21 = (-t353 + t584) * t643 + t473 + t617;
t1 = [(t179 * t95 + t180 * t94) * t748 + (t100 * t193 + t101 * t192) * t749 + (t182 * t292 + t183 * t291) * t750 + (t204 * t311 + t205 * t310) * t751 + (t289 * t334 + t290 * t333) * t752 + t448 * t373 + t449 * t372 + t394 * t684 + t406 * t680 - t393 * t689 - t405 * t681 + (t563 - t556) * t640 + (t557 + t562) * t639 - t800 * t612 + (t355 + t780) * t439 + t759 + t779 * (-t438 * t637 - t439 * t679) + (t356 + t770) * t438; (t457 / 0.2e1 + t458 / 0.2e1) * t551 * qJD(2) + t468 + (-qJD(2) * t523 + (qJD(1) * t382 - t463 * t491) * t465 + (qJD(1) * t384 - t463 * t492) * t462) * t745 + ((t696 / 0.2e1 + t694 / 0.2e1 - t334 * t743) * t466 + (t333 * t743 + t697 / 0.2e1 + t695 / 0.2e1) * t463) * qJD(1) + (-qJD(2) * t522 + (-qJD(1) * t381 - t466 * t491) * t465 + (-qJD(1) * t383 - t466 * t492) * t462) * t746 + m(3) * ((-t289 * t463 - t290 * t466) * t426 + (-t333 * t466 - t334 * t463) * t412) + m(4) * (t204 * t350 + t205 * t351 + t248 * t310 + t249 * t311) + m(5) * (t182 * t314 + t183 * t315 + t198 * t291 + t199 * t292) + m(6) * (t100 * t217 + t101 * t218 + t128 * t192 + t129 * t193) + m(7) * (t177 * t94 + t178 * t95 + t179 * t91 + t180 * t92); ((t386 * t463 + t387 * t466) * ((qJD(1) * t386 + t476) * t466 + (-t769 + (-t387 - t635 + t452) * qJD(1)) * t463) + t647 * t426 * t412) * t752 + (t181 * t93 + t248 * t351 + t249 * t350) * t751 + (t121 * t56 + t198 * t315 + t199 * t314) * t750 + (t128 * t218 + t129 * t217 + t26 * t96) * t749 + (t177 * t92 + t178 * t91 + t21 * t85) * t748 - t466 * t18 - t466 * t17 - t466 * t28 - t466 * t25 + (-t206 * t466 + t207 * t463) * t643 + (-t208 * t466 + t209 * t463) * t642 + t475 + t463 * ((t463 * t301 + (t208 + t768) * qJD(1)) * t463 + (t209 * qJD(1) + (t381 * t639 + t383 * t640) * t466 + (-t302 + (-t694 - t696) * qJD(2) + (t380 - t523) * qJD(1)) * t463) * t466) - t466 * ((t466 * t302 + (t207 + t764) * qJD(1)) * t466 + (t206 * qJD(1) + (-t382 * t639 - t384 * t640 + t644) * t463 + (-t301 + (t695 + t697) * qJD(2) - t522 * qJD(1)) * t466) * t463); m(5) * (t182 * t329 + t183 * t330 + t215 * t291 + t216 * t292) + m(6) * (t100 * t231 + t101 * t232 + t136 * t192 + t137 * t193) + m(7) * (t179 * t97 + t180 * t98 + t184 * t94 + t185 * t95) + t468 + (-t204 * t463 - t205 * t466 + (t310 * t463 - t311 * t466) * qJD(1)) * t742 + m(4) * (-t310 * t466 - t311 * t463) * t378; m(4) * (-t351 * t378 * t466 + t138 * t181 + t266 * t93 - t350 * t698) + m(7) * (t177 * t98 + t178 * t97 + t184 * t92 + t185 * t91 + t21 * t90 + t22 * t85) + m(6) * (t102 * t26 + t128 * t232 + t129 * t231 + t136 * t218 + t137 * t217 + t37 * t96) + m(5) * (t121 * t76 + t167 * t56 + t198 * t330 + t199 * t329 + t215 * t315 + t216 * t314) + t470 + (-t248 * t466 - t249 * t463 + (-t350 * t466 + t351 * t463) * qJD(1)) * t742; (t184 * t98 + t185 * t97 + t22 * t90) * t748 + (t102 * t37 + t136 * t232 + t137 * t231) * t749 + (t167 * t76 + t215 * t330 + t216 * t329) * t750 + t470 + (t378 * t408 * t647 + t138 * t266) * t751; m(7) * (qJD(1) * t536 + t463 * t95 - t466 * t94) + m(6) * (-t100 * t466 + t101 * t463 + (t192 * t466 + t193 * t463) * qJD(1)) + m(5) * (-t182 * t466 + t183 * t463 + (t291 * t466 + t292 * t463) * qJD(1)); m(7) * (qJD(1) * t537 + t463 * t91 - t466 * t92) + m(6) * (t128 * t463 - t129 * t466 + (t217 * t463 + t218 * t466) * qJD(1)) + m(5) * (t198 * t463 - t199 * t466 + (t314 * t463 + t315 * t466) * qJD(1)); m(7) * (qJD(1) * t535 + t463 * t97 - t466 * t98) + m(6) * (t136 * t463 - t137 * t466 + (t231 * t463 + t232 * t466) * qJD(1)) + m(5) * (t215 * t463 - t216 * t466 + (t329 * t463 + t330 * t466) * qJD(1)); 0; m(6) * (t100 * t191 + t101 * t190 + t192 * t89 + t193 * t88) + m(7) * (t109 * t95 + t110 * t94 + t179 * t47 + t180 * t46) + ((t463 * t506 + t466 * t505) * t456 - t730) * t439 + ((t45 / 0.2e1 + t43 / 0.2e1 + t62 / 0.2e1 + t63 / 0.2e1) * t466 + (t42 / 0.2e1 + t65 / 0.2e1 + t64 / 0.2e1 + t44 / 0.2e1) * t463 + (-t463 * t505 + t466 * t506) * qJD(1)) * t438 - t666; t469 + m(6) * (t128 * t190 + t129 * t191 + t174 * t26 + t217 * t88 + t218 * t89 + t57 * t96) + m(7) * (t109 * t91 + t110 * t92 + t177 * t46 + t178 * t47 + t21 * t99 + t23 * t85); t469 + m(6) * (t102 * t57 + t136 * t190 + t137 * t191 + t174 * t37 + t231 * t88 + t232 * t89) + m(7) * (t109 * t97 + t110 * t98 + t184 * t46 + t185 * t47 + t22 * t99 + t23 * t90); m(6) * (t463 * t89 - t466 * t88 + (t190 * t466 + t191 * t463) * qJD(1)) + m(7) * (qJD(1) * t546 - t46 * t466 + t463 * t47); (t109 * t47 + t110 * t46 + t23 * t99) * t748 + (t174 * t57 + t190 * t89 + t191 * t88) * t749 + (t730 * t439 + ((-t439 * t669 + t791) * t466 + (-t439 * t670 + t732) * t463) * t456 + t666) * t439 + (t795 * t466 + t794 * t463 + t797 * t689 + ((-t43 - t45) * t466 + (-t42 - t44) * t463 + t782 * t456) * t439 + (t439 * t798 - t463 * t791 + t732 * t466) * qJD(1)) * t438; m(7) * (t536 * t684 + (t463 * t94 + t466 * t95 + (-t179 * t463 + t180 * t466) * qJD(1)) * t438); m(7) * ((t456 * t537 - t21) * t439 + (t456 * t85 + t463 * t92 + t466 * t91 + (t177 * t466 - t178 * t463) * qJD(1)) * t438); m(7) * ((t456 * t535 - t22) * t439 + (t456 * t90 + t463 * t98 + t466 * t97 + (t184 * t466 - t185 * t463) * qJD(1)) * t438); 0; m(7) * ((t456 * t546 - t23) * t439 + (t456 * t99 + t46 * t463 + t466 * t47 + (-t109 * t463 + t110 * t466) * qJD(1)) * t438); (-0.1e1 + t647) * t438 * t684 * t748;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;