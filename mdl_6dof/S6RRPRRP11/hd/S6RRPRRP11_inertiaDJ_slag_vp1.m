% Calculate time derivative of joint inertia matrix for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP11_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:28
% EndTime: 2019-03-09 12:46:26
% DurationCPUTime: 41.11s
% Computational Cost: add. (49161->1286), mult. (75886->1696), div. (0->0), fcn. (71779->8), ass. (0->589)
t522 = qJ(4) + qJ(5);
t506 = sin(t522);
t524 = sin(qJ(2));
t527 = cos(qJ(2));
t507 = cos(t522);
t787 = Icges(7,4) * t507;
t619 = Icges(7,1) * t506 + t787;
t391 = Icges(7,5) * t524 - t527 * t619;
t789 = Icges(6,4) * t507;
t620 = Icges(6,1) * t506 + t789;
t392 = Icges(6,5) * t524 - t527 * t620;
t872 = -t392 - t391;
t884 = t872 * t506;
t609 = Icges(7,5) * t506 + Icges(7,6) * t507;
t387 = Icges(7,3) * t524 - t527 * t609;
t610 = Icges(6,5) * t506 + Icges(6,6) * t507;
t388 = Icges(6,3) * t524 - t527 * t610;
t883 = t387 + t388;
t788 = Icges(7,4) * t506;
t613 = Icges(7,2) * t507 + t788;
t389 = Icges(7,6) * t524 - t527 * t613;
t790 = Icges(6,4) * t506;
t614 = Icges(6,2) * t507 + t790;
t390 = Icges(6,6) * t524 - t527 * t614;
t882 = t389 + t390;
t783 = Icges(4,6) * t524;
t795 = Icges(3,4) * t524;
t881 = -t783 - t795 + (-Icges(3,2) - Icges(4,3)) * t527;
t782 = Icges(4,6) * t527;
t794 = Icges(3,4) * t527;
t880 = -t782 - t794 + (-Icges(3,1) - Icges(4,2)) * t524;
t525 = sin(qJ(1));
t528 = cos(qJ(1));
t517 = qJD(4) + qJD(5);
t725 = qJD(1) * t524;
t662 = t517 + t725;
t719 = qJD(2) * t527;
t683 = t525 * t719;
t538 = t528 * t662 + t683;
t663 = t517 * t524 + qJD(1);
t574 = t506 * t663;
t262 = t507 * t538 - t525 * t574;
t573 = t507 * t663;
t263 = t506 * t538 + t525 * t573;
t772 = t517 * t527;
t276 = (-Icges(7,5) * t507 + Icges(7,6) * t506) * t772 + (Icges(7,3) * t527 + t524 * t609) * qJD(2);
t278 = (Icges(7,2) * t506 - t787) * t772 + (Icges(7,6) * t527 + t524 * t613) * qJD(2);
t280 = (-Icges(7,1) * t507 + t788) * t772 + (Icges(7,5) * t527 + t524 * t619) * qJD(2);
t770 = t524 * t525;
t426 = t506 * t528 + t507 * t770;
t427 = t506 * t770 - t507 * t528;
t720 = qJD(2) * t525;
t686 = t524 * t720;
t722 = qJD(1) * t528;
t690 = t527 * t722;
t551 = -t686 + t690;
t767 = t525 * t527;
t82 = t262 * t389 + t263 * t391 + t276 * t767 + t278 * t426 + t280 * t427 + t387 * t551;
t277 = (-Icges(6,5) * t507 + Icges(6,6) * t506) * t772 + (Icges(6,3) * t527 + t524 * t610) * qJD(2);
t279 = (Icges(6,2) * t506 - t789) * t772 + (Icges(6,6) * t527 + t524 * t614) * qJD(2);
t281 = (-Icges(6,1) * t507 + t790) * t772 + (Icges(6,5) * t527 + t524 * t620) * qJD(2);
t83 = t262 * t390 + t263 * t392 + t277 * t767 + t279 * t426 + t281 * t427 + t388 * t551;
t879 = t82 + t83;
t718 = qJD(2) * t528;
t682 = t527 * t718;
t537 = -t525 * t662 + t682;
t264 = t507 * t537 - t528 * t574;
t265 = t506 * t537 + t528 * t573;
t769 = t524 * t528;
t424 = -t506 * t525 + t507 * t769;
t425 = t506 * t769 + t507 * t525;
t765 = t527 * t528;
t685 = t524 * t718;
t723 = qJD(1) * t527;
t691 = t525 * t723;
t844 = t685 + t691;
t84 = t264 * t389 + t265 * t391 + t276 * t765 + t278 * t424 + t280 * t425 - t387 * t844;
t85 = t264 * t390 + t265 * t392 + t277 * t765 + t279 * t424 + t281 * t425 - t388 * t844;
t878 = t84 + t85;
t288 = Icges(7,5) * t425 + Icges(7,6) * t424 + Icges(7,3) * t765;
t292 = Icges(7,4) * t425 + Icges(7,2) * t424 + Icges(7,6) * t765;
t296 = Icges(7,1) * t425 + Icges(7,4) * t424 + Icges(7,5) * t765;
t128 = t288 * t765 + t292 * t424 + t296 * t425;
t289 = Icges(7,5) * t427 + Icges(7,6) * t426 + Icges(7,3) * t767;
t293 = Icges(7,4) * t427 + Icges(7,2) * t426 + Icges(7,6) * t767;
t297 = Icges(7,1) * t427 + Icges(7,4) * t426 + Icges(7,5) * t767;
t129 = t289 * t765 + t293 * t424 + t297 * t425;
t290 = Icges(6,5) * t425 + Icges(6,6) * t424 + Icges(6,3) * t765;
t294 = Icges(6,4) * t425 + Icges(6,2) * t424 + Icges(6,6) * t765;
t298 = Icges(6,1) * t425 + Icges(6,4) * t424 + Icges(6,5) * t765;
t130 = t290 * t765 + t294 * t424 + t298 * t425;
t291 = Icges(6,5) * t427 + Icges(6,6) * t426 + Icges(6,3) * t767;
t295 = Icges(6,4) * t427 + Icges(6,2) * t426 + Icges(6,6) * t767;
t299 = Icges(6,1) * t427 + Icges(6,4) * t426 + Icges(6,5) * t767;
t131 = t291 * t765 + t295 * t424 + t299 * t425;
t849 = (-t129 - t131) * t528 + (t128 + t130) * t525;
t132 = t288 * t767 + t292 * t426 + t296 * t427;
t133 = t289 * t767 + t293 * t426 + t297 * t427;
t134 = t290 * t767 + t294 * t426 + t298 * t427;
t135 = t291 * t767 + t295 * t426 + t299 * t427;
t848 = (-t133 - t135) * t528 + (t132 + t134) * t525;
t592 = t292 * t507 + t296 * t506;
t152 = t288 * t524 - t527 * t592;
t590 = t294 * t507 + t298 * t506;
t154 = t290 * t524 - t527 * t590;
t877 = t152 + t154;
t591 = t293 * t507 + t297 * t506;
t153 = t289 * t524 - t527 * t591;
t589 = t295 * t507 + t299 * t506;
t155 = t291 * t524 - t527 * t589;
t876 = -t153 - t155;
t195 = t387 * t765 + t389 * t424 + t391 * t425;
t196 = t388 * t765 + t390 * t424 + t392 * t425;
t875 = t195 + t196;
t197 = t387 * t767 + t389 * t426 + t391 * t427;
t198 = t388 * t767 + t390 * t426 + t392 * t427;
t874 = t197 + t198;
t873 = -t279 - t278;
t871 = t881 * qJD(2);
t870 = t880 * qJD(2);
t600 = t134 * t528 + t135 * t525;
t601 = t132 * t528 + t133 * t525;
t869 = t600 + t601;
t602 = t130 * t528 + t131 * t525;
t603 = t128 * t528 + t129 * t525;
t868 = t602 + t603;
t721 = qJD(2) * t524;
t867 = t882 * (t506 * t772 + t507 * t721) + t883 * t719 - t721 * t884 + (t276 + t277) * t524;
t866 = ((-t882 * t507 + t884) * t527 + t883 * t524) * t719;
t865 = (-t281 - t280) * t506;
t526 = cos(qJ(4));
t501 = t526 * pkin(4) + pkin(3);
t458 = pkin(5) * t507 + t501;
t523 = sin(qJ(4));
t808 = pkin(4) * t523;
t461 = pkin(5) * t506 + t808;
t864 = t425 * rSges(7,1) + t424 * rSges(7,2) + rSges(7,3) * t765 + t525 * t458 + t461 * t769;
t618 = -Icges(3,2) * t524 + t794;
t406 = -Icges(3,6) * t528 + t525 * t618;
t606 = -Icges(4,3) * t524 + t782;
t830 = Icges(4,5) * t528 + t525 * t606;
t863 = -t406 - t830;
t407 = Icges(3,6) * t525 + t528 * t618;
t411 = Icges(4,5) * t525 - t528 * t606;
t862 = t407 - t411;
t623 = Icges(3,1) * t527 - t795;
t409 = -Icges(3,5) * t528 + t525 * t623;
t608 = Icges(4,2) * t527 - t783;
t829 = Icges(4,4) * t528 + t525 * t608;
t861 = t409 + t829;
t410 = Icges(3,5) * t525 + t528 * t623;
t413 = Icges(4,4) * t525 - t528 * t608;
t860 = t410 - t413;
t708 = qJD(4) * t808;
t807 = pkin(5) * t517;
t438 = -t506 * t807 - t708;
t716 = qJD(4) * t526;
t707 = pkin(4) * t716;
t439 = t507 * t807 + t707;
t529 = -pkin(9) - pkin(8);
t516 = -qJ(6) + t529;
t714 = qJD(6) * t527;
t859 = t265 * rSges(7,1) + t264 * rSges(7,2) + t525 * t438 + t439 * t769 + t458 * t722 + t461 * t682 + t516 * t844 + t528 * t714;
t858 = -rSges(7,1) * t427 - rSges(7,2) * t426 + t458 * t528;
t672 = -t721 / 0.2e1;
t857 = t528 * t672 - t691 / 0.2e1;
t673 = t722 / 0.2e1;
t856 = t525 * t672 + t527 * t673;
t631 = rSges(7,1) * t506 + rSges(7,2) * t507;
t669 = -t461 + t808;
t729 = t516 - t529;
t827 = t524 * t729 - t527 * t669;
t742 = rSges(7,3) * t524 - t527 * t631 - t827;
t170 = Icges(7,5) * t265 + Icges(7,6) * t264 - Icges(7,3) * t844;
t174 = Icges(7,4) * t265 + Icges(7,2) * t264 - Icges(7,6) * t844;
t178 = Icges(7,1) * t265 + Icges(7,4) * t264 - Icges(7,5) * t844;
t36 = t170 * t767 + t174 * t426 + t178 * t427 + t262 * t292 + t263 * t296 + t288 * t551;
t169 = Icges(7,5) * t263 + Icges(7,6) * t262 + Icges(7,3) * t551;
t173 = Icges(7,4) * t263 + Icges(7,2) * t262 + Icges(7,6) * t551;
t177 = Icges(7,1) * t263 + Icges(7,4) * t262 + Icges(7,5) * t551;
t37 = t169 * t767 + t173 * t426 + t177 * t427 + t262 * t293 + t263 * t297 + t289 * t551;
t172 = Icges(6,5) * t265 + Icges(6,6) * t264 - Icges(6,3) * t844;
t176 = Icges(6,4) * t265 + Icges(6,2) * t264 - Icges(6,6) * t844;
t180 = Icges(6,1) * t265 + Icges(6,4) * t264 - Icges(6,5) * t844;
t38 = t172 * t767 + t176 * t426 + t180 * t427 + t262 * t294 + t263 * t298 + t290 * t551;
t171 = Icges(6,5) * t263 + Icges(6,6) * t262 + Icges(6,3) * t551;
t175 = Icges(6,4) * t263 + Icges(6,2) * t262 + Icges(6,6) * t551;
t179 = Icges(6,1) * t263 + Icges(6,4) * t262 + Icges(6,5) * t551;
t39 = t171 * t767 + t175 * t426 + t179 * t427 + t262 * t295 + t263 * t299 + t291 * t551;
t855 = ((t36 + t38) * t528 + (t37 + t39) * t525 + t874 * qJD(2) - t848 * qJD(1)) * t527 + (-qJD(2) * t869 + t879) * t524;
t40 = t170 * t765 + t174 * t424 + t178 * t425 + t264 * t292 + t265 * t296 - t288 * t844;
t41 = t169 * t765 + t173 * t424 + t177 * t425 + t264 * t293 + t265 * t297 - t289 * t844;
t42 = t172 * t765 + t176 * t424 + t180 * t425 + t264 * t294 + t265 * t298 - t290 * t844;
t43 = t171 * t765 + t175 * t424 + t179 * t425 + t264 * t295 + t265 * t299 - t291 * t844;
t854 = ((t40 + t42) * t528 + (t41 + t43) * t525 + t875 * qJD(2) - t849 * qJD(1)) * t527 + (-qJD(2) * t868 + t878) * t524;
t50 = (qJD(2) * t592 + t170) * t524 + (qJD(2) * t288 + (-t296 * t517 - t174) * t507 + (t292 * t517 - t178) * t506) * t527;
t52 = (qJD(2) * t590 + t172) * t524 + (qJD(2) * t290 + (-t298 * t517 - t176) * t507 + (t294 * t517 - t180) * t506) * t527;
t853 = t50 + t52;
t51 = (qJD(2) * t591 + t169) * t524 + (qJD(2) * t289 + (-t297 * t517 - t173) * t507 + (t293 * t517 - t177) * t506) * t527;
t53 = (qJD(2) * t589 + t171) * t524 + (qJD(2) * t291 + (-t299 * t517 - t175) * t507 + (t295 * t517 - t179) * t506) * t527;
t852 = t51 + t53;
t851 = t524 * t875 + t527 * t868;
t850 = t524 * t874 + t527 * t869;
t689 = t529 * t723;
t475 = t528 * t689;
t642 = t669 * t524;
t544 = -t516 * t527 - t642;
t633 = rSges(7,1) * t263 + rSges(7,2) * t262;
t646 = -t439 + t707;
t759 = t475 + (qJD(1) * t544 - t438 - t708) * t528 + (t714 - t646 * t524 + (t458 - t501) * qJD(1) + t827 * qJD(2)) * t525 + rSges(7,3) * t551 + t633;
t660 = qJD(4) + t725;
t562 = t660 * t808;
t681 = t524 * t716;
t684 = t529 * t721;
t624 = t501 * t722 + t525 * t689 + t682 * t808 + (pkin(4) * t681 + t684) * t528;
t758 = (-t461 * t725 + t562) * t525 - t624 - rSges(7,3) * t844 + t859;
t664 = t729 * t527;
t847 = qJD(6) * t524 + t646 * t527 + (-rSges(7,1) * t507 + rSges(7,2) * t506) * t772 + (rSges(7,3) * t527 + t524 * t631 - t642 - t664) * qJD(2);
t703 = t523 * t769;
t734 = -pkin(4) * t703 - t525 * t501;
t752 = -t528 * t664 + t734 + t864;
t764 = t527 * t529;
t733 = t528 * t501 + t525 * t764;
t751 = rSges(7,3) * t767 + t525 * t544 + t733 - t858;
t846 = t525 * t877 + t528 * t876;
t845 = t525 * t876 - t528 * t877;
t842 = pkin(8) * t723 - t562;
t841 = t866 + (((t517 * t872 + t873) * t507 + t865) * t527 + t867) * t524;
t792 = Icges(5,4) * t523;
t615 = Icges(5,2) * t526 + t792;
t405 = Icges(5,6) * t524 - t527 * t615;
t791 = Icges(5,4) * t526;
t621 = Icges(5,1) * t523 + t791;
t408 = Icges(5,5) * t524 - t527 * t621;
t840 = t405 * t526 + t408 * t523;
t837 = qJD(2) / 0.2e1;
t577 = t411 * t524 - t413 * t527;
t836 = t525 * t577;
t578 = t407 * t524 - t410 * t527;
t835 = t525 * t578;
t576 = -t524 * t830 + t527 * t829;
t834 = t528 * t576;
t579 = t406 * t524 - t409 * t527;
t833 = t528 * t579;
t832 = t525 * rSges(4,1) - rSges(4,2) * t765;
t459 = t525 * pkin(3) + pkin(8) * t765;
t831 = -rSges(3,2) * t769 + t525 * rSges(3,3);
t612 = Icges(3,5) * t527 - Icges(3,6) * t524;
t403 = -Icges(3,3) * t528 + t525 * t612;
t616 = Icges(4,4) * t527 - Icges(4,5) * t524;
t828 = Icges(4,1) * t528 + t525 * t616;
t826 = -t525 * t751 - t528 * t752;
t825 = 2 * m(3);
t824 = 2 * m(4);
t823 = 2 * m(5);
t822 = 2 * m(6);
t821 = 2 * m(7);
t820 = m(4) / 0.2e1;
t819 = m(5) / 0.2e1;
t818 = m(6) / 0.2e1;
t817 = m(7) / 0.2e1;
t816 = t524 / 0.2e1;
t815 = t525 / 0.2e1;
t814 = -t528 / 0.2e1;
t813 = t528 / 0.2e1;
t812 = -rSges(6,3) - pkin(2);
t811 = -rSges(7,3) - pkin(2);
t478 = rSges(3,1) * t524 + rSges(3,2) * t527;
t810 = m(3) * t478;
t809 = pkin(2) * t527;
t806 = qJD(1) / 0.2e1;
t805 = -pkin(8) - t529;
t802 = rSges(4,1) * t528;
t801 = rSges(4,2) * t524;
t800 = rSges(3,3) * t528;
t799 = rSges(5,3) * t524;
t798 = -rSges(4,3) - qJ(3);
t797 = rSges(7,3) - t516;
t780 = qJ(3) * t524;
t779 = qJ(3) * t527;
t768 = t525 * t526;
t442 = t523 * t528 + t524 * t768;
t766 = t526 * t528;
t443 = t523 * t770 - t766;
t638 = -rSges(5,1) * t443 - rSges(5,2) * t442;
t332 = rSges(5,3) * t767 - t638;
t778 = t332 * t528;
t771 = t523 * t527;
t763 = -qJ(3) - t461;
t611 = Icges(5,5) * t523 + Icges(5,6) * t526;
t402 = Icges(5,3) * t524 - t527 * t611;
t240 = t402 * t524 - t527 * t840;
t715 = qJD(4) * t527;
t338 = (Icges(5,2) * t523 - t791) * t715 + (Icges(5,6) * t527 + t524 * t615) * qJD(2);
t341 = (-Icges(5,1) * t526 + t792) * t715 + (Icges(5,5) * t527 + t524 * t621) * qJD(2);
t335 = (-Icges(5,5) * t526 + Icges(5,6) * t523) * t715 + (Icges(5,3) * t527 + t524 * t611) * qJD(2);
t625 = t523 * t405 * t715 + t524 * t335 + t402 * t719 + t721 * t840;
t760 = ((-t341 * t523 + (-qJD(4) * t408 - t338) * t526) * t527 + t625) * t524 + t240 * t719;
t636 = t263 * rSges(6,1) + t262 * rSges(6,2);
t182 = rSges(6,3) * t551 + t636;
t301 = t425 * rSges(6,1) + t424 * rSges(6,2) + rSges(6,3) * t765;
t757 = t182 * t765 + t301 * t686;
t488 = pkin(8) * t686;
t223 = -t475 + t488 - t842 * t528 + (t684 + (-pkin(3) + t501) * qJD(1) + (t523 * t719 + t681) * pkin(4)) * t525;
t565 = -t528 * t764 - t734;
t352 = t565 - t459;
t756 = t223 * t765 + t352 * t686;
t755 = t751 * t765;
t754 = t752 * t524;
t749 = t265 * rSges(6,1) + t264 * rSges(6,2);
t634 = rSges(6,1) * t506 + rSges(6,2) * t507;
t285 = (-rSges(6,1) * t507 + rSges(6,2) * t506) * t772 + (rSges(6,3) * t527 + t524 * t634) * qJD(2);
t395 = rSges(6,3) * t524 - t527 * t634;
t748 = t285 * t767 + t395 * t690;
t747 = -t301 - t352;
t635 = -rSges(6,1) * t427 - rSges(6,2) * t426;
t303 = rSges(6,3) * t767 - t635;
t514 = t528 * pkin(3);
t711 = t524 * t808;
t353 = t514 + (-pkin(8) * t527 + t711) * t525 - t733;
t746 = -t303 - t353;
t535 = -t525 * t660 + t682;
t661 = qJD(4) * t524 + qJD(1);
t570 = t661 * t523;
t309 = t526 * t535 - t528 * t570;
t571 = t526 * t661;
t310 = t523 * t535 + t528 * t571;
t745 = t310 * rSges(5,1) + t309 * rSges(5,2);
t744 = t742 * t767;
t680 = t526 * t715;
t368 = -pkin(4) * t680 + (t527 * t805 + t711) * qJD(2);
t431 = -pkin(4) * t771 + t524 * t805;
t743 = t368 * t767 + t431 * t690;
t741 = -t395 - t431;
t628 = t780 + t809;
t447 = t628 * t525;
t448 = pkin(2) * t765 + qJ(3) * t769;
t740 = t525 * t447 + t528 * t448;
t432 = qJD(2) * t628 - qJD(3) * t527;
t630 = -rSges(4,2) * t527 + rSges(4,3) * t524;
t739 = -t630 * qJD(2) - t432;
t737 = -t448 - t459;
t476 = pkin(2) * t524 - t779;
t724 = qJD(1) * t525;
t449 = t476 * t724;
t692 = t524 * t724;
t736 = pkin(8) * t692 + t449;
t629 = rSges(4,3) * t527 + t801;
t735 = -t476 + t629;
t717 = qJD(3) * t524;
t732 = qJ(3) * t682 + t528 * t717;
t731 = rSges(3,2) * t692 + rSges(3,3) * t722;
t730 = t528 * pkin(1) + t525 * pkin(7);
t728 = t525 ^ 2 + t528 ^ 2;
t404 = Icges(3,3) * t525 + t528 * t612;
t727 = qJD(1) * t404;
t415 = Icges(4,1) * t525 - t528 * t616;
t726 = qJD(1) * t415;
t713 = -rSges(5,3) - pkin(2) - pkin(8);
t712 = -0.1e1 + t728;
t440 = -t523 * t525 + t524 * t766;
t441 = t703 + t768;
t104 = t309 * t405 + t310 * t408 + t335 * t765 + t338 * t440 + t341 * t441 - t402 * t844;
t202 = Icges(5,5) * t310 + Icges(5,6) * t309 - Icges(5,3) * t844;
t204 = Icges(5,4) * t310 + Icges(5,2) * t309 - Icges(5,6) * t844;
t206 = Icges(5,1) * t310 + Icges(5,4) * t309 - Icges(5,5) * t844;
t319 = Icges(5,5) * t441 + Icges(5,6) * t440 + Icges(5,3) * t765;
t321 = Icges(5,4) * t441 + Icges(5,2) * t440 + Icges(5,6) * t765;
t323 = Icges(5,1) * t441 + Icges(5,4) * t440 + Icges(5,5) * t765;
t587 = t321 * t526 + t323 * t523;
t58 = (qJD(2) * t587 + t202) * t524 + (qJD(2) * t319 - t204 * t526 - t206 * t523 + (t321 * t523 - t323 * t526) * qJD(4)) * t527;
t702 = t58 / 0.2e1 + t104 / 0.2e1;
t569 = t660 * t528;
t307 = t526 * t569 + (t526 * t719 - t570) * t525;
t308 = t525 * t571 + (t569 + t683) * t523;
t103 = t307 * t405 + t308 * t408 + t335 * t767 + t338 * t442 + t341 * t443 + t402 * t551;
t201 = Icges(5,5) * t308 + Icges(5,6) * t307 + Icges(5,3) * t551;
t203 = Icges(5,4) * t308 + Icges(5,2) * t307 + Icges(5,6) * t551;
t205 = Icges(5,1) * t308 + Icges(5,4) * t307 + Icges(5,5) * t551;
t320 = Icges(5,5) * t443 + Icges(5,6) * t442 + Icges(5,3) * t767;
t322 = Icges(5,4) * t443 + Icges(5,2) * t442 + Icges(5,6) * t767;
t324 = Icges(5,1) * t443 + Icges(5,4) * t442 + Icges(5,5) * t767;
t586 = t322 * t526 + t324 * t523;
t59 = (qJD(2) * t586 + t201) * t524 + (qJD(2) * t320 - t203 * t526 - t205 * t523 + (t322 * t523 - t324 * t526) * qJD(4)) * t527;
t701 = t59 / 0.2e1 + t103 / 0.2e1;
t700 = -t352 - t752;
t699 = -t353 - t751;
t489 = pkin(2) * t686;
t698 = t525 * (pkin(2) * t690 + t525 * t717 - t489 + (t524 * t722 + t683) * qJ(3)) + t528 * (-pkin(2) * t844 - qJ(3) * t692 + t732) + t447 * t722;
t697 = -t352 + t737;
t696 = -t431 - t742;
t695 = t431 * t724 + t736;
t331 = t441 * rSges(5,1) + t440 * rSges(5,2) + rSges(5,3) * t765;
t504 = pkin(7) * t722;
t694 = t504 + t732;
t637 = rSges(5,1) * t523 + rSges(5,2) * t526;
t420 = -t527 * t637 + t799;
t693 = t420 * t724;
t679 = t767 / 0.2e1;
t678 = t765 / 0.2e1;
t167 = t319 * t524 - t527 * t587;
t215 = t402 * t765 + t405 * t440 + t408 * t441;
t677 = -t215 / 0.2e1 - t167 / 0.2e1;
t168 = t320 * t524 - t527 * t586;
t216 = t402 * t767 + t405 * t442 + t408 * t443;
t676 = t216 / 0.2e1 + t168 / 0.2e1;
t675 = t612 * t837 - t616 * qJD(2) / 0.2e1;
t670 = -qJ(3) - t808;
t668 = -t524 * pkin(8) - t476;
t158 = t319 * t767 + t321 * t442 + t323 * t443;
t159 = t320 * t767 + t322 * t442 + t324 * t443;
t594 = t158 * t528 + t159 * t525;
t77 = t524 * t216 + t527 * t594;
t667 = t168 * t524 + t77;
t666 = t751 * t524;
t665 = t746 * t524;
t380 = t735 * t528;
t659 = t686 * t752 + t759 * t765;
t184 = -rSges(6,3) * t844 + t749;
t658 = t524 * t184 + t301 * t719 + t395 * t844;
t505 = pkin(3) * t722;
t224 = pkin(8) * t685 + t525 * t842 - t505 + t624;
t657 = t524 * t224 + t352 * t719 + t431 * t844;
t656 = t690 * t742 + t767 * t847;
t460 = pkin(8) * t767 - t514;
t655 = t528 * t459 + t525 * t460 + t740;
t654 = rSges(4,1) * t722 + rSges(4,2) * t844 + rSges(4,3) * t682;
t653 = t730 + t448;
t648 = -t420 + t668;
t647 = -t431 + t668;
t645 = -pkin(8) * t719 - t432;
t644 = t699 * t524;
t643 = t524 * t798 - pkin(1);
t641 = t524 * t763 - pkin(1);
t640 = rSges(3,1) * t527 - rSges(3,2) * t524;
t639 = rSges(5,1) * t308 + rSges(5,2) * t307;
t111 = t811 * t685 + (t527 * t811 + t641) * t724 + t694 + t859;
t536 = (-pkin(2) - t797) * t527 + t641;
t112 = t489 + (qJD(1) * t536 + t438) * t528 + (-t714 + (-qJD(3) - t439) * t524 + (-pkin(7) - t458) * qJD(1) + (t524 * t797 + t527 * t763) * qJD(2)) * t525 - t633;
t604 = t111 * t525 + t112 * t528;
t156 = t319 * t765 + t321 * t440 + t323 * t441;
t157 = t320 * t765 + t322 * t440 + t324 * t441;
t595 = t156 * t528 + t157 * t525;
t107 = t156 * t525 - t157 * t528;
t108 = t158 * t525 - t159 * t528;
t513 = t528 * pkin(7);
t199 = t525 * t536 + t513 + t858;
t200 = -t516 * t765 + t653 + t864;
t593 = t199 * t528 + t200 * t525;
t588 = -t301 * t528 - t303 * t525;
t585 = -t331 * t528 - t332 * t525;
t584 = t331 * t525 - t778;
t76 = t524 * t215 + t527 * t595;
t575 = -t167 * t524 - t76 - t851;
t572 = -t395 + t647;
t421 = rSges(3,1) * t765 + t831;
t422 = rSges(4,3) * t769 + t832;
t351 = (-rSges(5,1) * t526 + rSges(5,2) * t523) * t715 + (rSges(5,3) * t527 + t524 * t637) * qJD(2);
t568 = -t351 + t645;
t567 = -t368 + t645;
t566 = -pkin(1) - t640;
t318 = t648 * t528;
t564 = t525 * (qJD(1) * t459 - t488) + t528 * (-pkin(8) * t844 + t505) + t460 * t722 + t698;
t563 = t528 * t352 + t525 * t353 + t655;
t561 = t647 - t742;
t560 = -t285 + t567;
t559 = qJD(2) * t478;
t556 = qJD(2) * (Icges(4,4) * t524 + Icges(4,5) * t527);
t555 = qJD(2) * (-Icges(3,5) * t524 - Icges(3,6) * t527);
t255 = t572 * t528;
t549 = t567 - t847;
t313 = t353 * t765;
t102 = t700 * t767 + t313 + t755;
t543 = t758 * t524 + t752 * t719 + t742 * t844;
t32 = (-t368 - t847) * t765 + t543 + t657;
t33 = (-t223 - t759) * t524 + (t527 * t699 + t696 * t770) * qJD(2) + t656 + t743;
t548 = qJD(2) * t102 + t32 * t525 + t33 * t528;
t119 = -t752 * t767 + t755;
t44 = -t765 * t847 + t543;
t45 = -t759 * t524 + (-t527 * t751 - t742 * t770) * qJD(2) + t656;
t547 = qJD(2) * t119 + t44 * t525 + t45 * t528;
t546 = t527 * t713 - pkin(1) - t780;
t113 = t528 * t549 + t724 * t742 + t695;
t211 = t561 * t528;
t114 = qJD(1) * t211 + t525 * t549;
t95 = t563 - t826;
t545 = qJD(2) * t95 + t113 * t528 + t114 * t525;
t542 = t855 * t767 + t854 * t765 + t850 * t690 - t845 * t719 * t527 + (t845 * t721 + t866 + (-qJD(1) * t846 + t525 * t852 + t528 * t853) * t527 + t841) * t524;
t541 = (rSges(4,2) - pkin(2)) * t527 + t643;
t540 = t525 * t223 + t528 * t224 + t353 * t722 + t564;
t539 = t546 * t525;
t534 = t524 * t670 + t527 * t812 - pkin(1);
t533 = qJD(1) * t534 - t708;
t21 = qJD(1) * t601 + t36 * t525 - t37 * t528;
t22 = qJD(1) * t600 + t38 * t525 - t39 * t528;
t23 = qJD(1) * t603 + t40 * t525 - t41 * t528;
t24 = qJD(1) * t602 + t42 * t525 - t43 * t528;
t532 = (-qJD(1) * t845 + t525 * t853 - t528 * t852) * t816 + t854 * t815 + t855 * t814 + (t21 + t22) * t679 + (t23 + t24) * t678 + t850 * t724 / 0.2e1 + t851 * t673 + t846 * t719 / 0.2e1 + t857 * t849 + t856 * t848;
t531 = (-t525 * t850 - t528 * t851) * t721 - t851 * t691 + t542;
t530 = (t852 + t879) * t679 + (t853 + t878) * t678 + t841 + t857 * (t875 + t877) + t856 * (t874 - t876);
t457 = t640 * qJD(2);
t423 = t525 * t630 - t802;
t419 = t525 * t640 - t800;
t393 = t431 * t767;
t379 = t735 * t525;
t375 = t421 + t730;
t374 = t525 * t566 + t513 + t800;
t373 = t712 * t524 * t719;
t370 = t395 * t767;
t349 = qJD(1) * t828 + t528 * t556;
t348 = t525 * t556 + t726;
t337 = t525 * t555 + t727;
t336 = -qJD(1) * t403 + t528 * t555;
t334 = t422 + t653;
t333 = t525 * t541 + t513 + t802;
t328 = t524 * t352;
t317 = t648 * t525;
t287 = t478 * t720 + ((-rSges(3,3) - pkin(7)) * t525 + t566 * t528) * qJD(1);
t286 = -rSges(3,1) * t844 - rSges(3,2) * t682 - pkin(1) * t724 + t504 + t731;
t283 = t524 * t301;
t273 = t303 * t765;
t257 = qJD(1) * t380 + t525 * t739;
t256 = t528 * t739 - t629 * t724 + t449;
t254 = t572 * t525;
t251 = t331 * t524 - t420 * t765;
t250 = -t332 * t524 + t420 * t767;
t246 = t525 * t576 + t528 * t828;
t245 = -t415 * t528 + t836;
t244 = -t525 * t828 + t834;
t243 = t415 * t525 + t528 * t577;
t242 = t404 * t525 - t528 * t578;
t241 = t403 * t525 - t833;
t239 = -t404 * t528 - t835;
t238 = -t403 * t528 - t525 * t579;
t237 = t653 + t331 + t459;
t236 = t513 + t514 + t539 + t638;
t231 = t422 * t528 + t423 * t525 + t740;
t230 = -t395 * t765 + t283;
t229 = -t303 * t524 + t370;
t226 = t489 + (-t717 + (t527 * t798 - t801) * qJD(2)) * t525 + ((-rSges(4,1) - pkin(7)) * t525 + t541 * t528) * qJD(1);
t225 = -pkin(2) * t685 + (t643 - t809) * t724 + t654 + t694;
t220 = t584 * t527;
t214 = t565 + t653 + t301;
t213 = t525 * t534 + t513 + t635 + t733;
t210 = t561 * t525;
t209 = -t301 * t767 + t273;
t208 = -rSges(5,3) * t844 + t745;
t207 = rSges(5,3) * t551 + t639;
t194 = qJD(1) * t318 + t525 * t568;
t193 = t528 * t568 + t693 + t736;
t162 = t741 * t765 + t283 + t328;
t161 = t370 + t393 + t665;
t160 = -t585 + t655;
t143 = t488 + t489 + (-t717 + (-t779 + t799) * qJD(2)) * t525 + ((-pkin(3) - pkin(7)) * t525 + t546 * t528) * qJD(1) - t639;
t142 = qJD(1) * t539 + t685 * t713 + t505 + t694 + t745;
t141 = qJD(1) * t255 + t525 * t560;
t140 = t395 * t724 + t528 * t560 + t695;
t137 = -t742 * t765 + t754;
t136 = -t666 + t744;
t126 = t747 * t767 + t273 + t313;
t124 = t563 - t588;
t123 = (-t420 * t720 - t207) * t524 + (-qJD(2) * t332 + t351 * t525 + t420 * t722) * t527;
t122 = (t420 * t718 + t208) * t524 + (qJD(2) * t331 - t351 * t528 + t693) * t527;
t121 = t696 * t765 + t328 + t754;
t120 = t393 + t644 + t744;
t118 = t475 + t489 + t533 * t528 + ((-qJD(3) - t707) * t524 + (-pkin(7) - t501) * qJD(1) + (t670 * t527 + (rSges(6,3) - t529) * t524) * qJD(2)) * t525 - t636;
t117 = t525 * t533 + t685 * t812 + t624 + t694 + t749;
t116 = (qJD(1) * t423 + t654) * t528 + (t629 * t720 + (-t422 - t448 + t832) * qJD(1)) * t525 + t698;
t110 = -t182 * t524 + (-t303 * t527 - t395 * t770) * qJD(2) + t748;
t109 = -t285 * t765 + t658;
t86 = t584 * t721 + (qJD(1) * t585 + t207 * t528 - t208 * t525) * t527;
t73 = -t303 * t685 + (qJD(1) * t588 - t184 * t525) * t527 + t757;
t62 = t207 * t525 + t208 * t528 + (t778 + (-t331 + t737) * t525) * qJD(1) + t564;
t61 = (-t182 - t223) * t524 + (t527 * t746 + t741 * t770) * qJD(2) + t743 + t748;
t60 = (-t285 - t368) * t765 + t657 + t658;
t57 = t201 * t765 + t203 * t440 + t205 * t441 + t309 * t322 + t310 * t324 - t320 * t844;
t56 = t202 * t765 + t204 * t440 + t206 * t441 + t309 * t321 + t310 * t323 - t319 * t844;
t55 = t201 * t767 + t203 * t442 + t205 * t443 + t307 * t322 + t308 * t324 + t320 * t551;
t54 = t202 * t767 + t204 * t442 + t206 * t443 + t307 * t321 + t308 * t323 + t319 * t551;
t35 = t182 * t525 + t184 * t528 + (t303 * t528 + (-t301 + t697) * t525) * qJD(1) + t540;
t34 = t665 * t718 + ((-t184 - t224) * t525 + (t525 * t746 + t528 * t747) * qJD(1)) * t527 + t756 + t757;
t31 = -t666 * t718 + (qJD(1) * t826 - t758 * t525) * t527 + t659;
t30 = t758 * t528 + t759 * t525 + (t751 * t528 + (t697 - t752) * t525) * qJD(1) + t540;
t29 = t644 * t718 + ((-t224 - t758) * t525 + (t525 * t699 + t528 * t700) * qJD(1)) * t527 + t659 + t756;
t28 = qJD(1) * t595 + t525 * t56 - t528 * t57;
t27 = qJD(1) * t594 + t525 * t54 - t528 * t55;
t16 = (-qJD(2) * t595 + t104) * t524 + (-qJD(1) * t107 + qJD(2) * t215 + t525 * t57 + t528 * t56) * t527;
t15 = (-qJD(2) * t594 + t103) * t524 + (-qJD(1) * t108 + qJD(2) * t216 + t525 * t55 + t528 * t54) * t527;
t1 = [-t341 * t771 - t408 * t680 + (t111 * t200 + t112 * t199) * t821 + (t117 * t214 + t118 * t213) * t822 + (t142 * t237 + t143 * t236) * t823 + (t225 * t334 + t226 * t333) * t824 + (t286 * t375 + t287 * t374) * t825 + t625 + t872 * t507 * t772 + (t623 + t608 + t881) * t721 + (t618 + t606 - t880) * t719 + (-t526 * t338 + t507 * t873 + t865) * t527 + t867; m(4) * (t225 * t379 + t226 * t380 + t256 * t333 + t257 * t334) + m(5) * (t142 * t317 + t143 * t318 + t193 * t236 + t194 * t237) + m(6) * (t117 * t254 + t118 * t255 + t140 * t213 + t141 * t214) + m(7) * (t111 * t210 + t112 * t211 + t113 * t199 + t114 * t200) + (m(3) * (-t287 * t478 - t374 * t457) - t53 / 0.2e1 - t51 / 0.2e1 - t83 / 0.2e1 - t82 / 0.2e1 + t675 * t528 - t701) * t528 + (m(3) * (-t286 * t478 - t375 * t457) + t52 / 0.2e1 + t50 / 0.2e1 + t84 / 0.2e1 + t85 / 0.2e1 + t675 * t525 + t702) * t525 + ((-t862 * qJD(2) + t528 * t870) * t815 + (t863 * qJD(2) + t525 * t870) * t814 + (t814 * t860 - t815 * t861) * qJD(1)) * t524 + ((t860 * qJD(2) + t528 * t871) * t815 + (t861 * qJD(2) + t525 * t871) * t814 + (t814 * t862 + t815 * t863) * qJD(1)) * t527 + ((t195 / 0.2e1 + t196 / 0.2e1 + t154 / 0.2e1 + t152 / 0.2e1 - t375 * t810 + (t407 / 0.2e1 - t411 / 0.2e1) * t527 + (t410 / 0.2e1 - t413 / 0.2e1) * t524 - t677) * t528 + (t197 / 0.2e1 + t198 / 0.2e1 + t155 / 0.2e1 + t153 / 0.2e1 + t374 * t810 + (t406 / 0.2e1 + t830 / 0.2e1) * t527 + (t409 / 0.2e1 + t829 / 0.2e1) * t524 + t676) * t525) * qJD(1); t525 * ((t525 * t349 + (t244 - t836) * qJD(1)) * t525 + (t243 * qJD(1) + (t719 * t830 + t721 * t829) * t528 + (-t348 + (t411 * t527 + t413 * t524) * qJD(2) + (t415 + t576) * qJD(1)) * t525) * t528) + (t113 * t211 + t114 * t210 + t30 * t95) * t821 + (t124 * t35 + t140 * t255 + t141 * t254) * t822 + (t160 * t62 + t193 * t318 + t194 * t317) * t823 + (t116 * t231 + t256 * t380 + t257 * t379) * t824 + ((t419 * t525 + t421 * t528) * ((qJD(1) * t419 - t528 * t559 + t731) * t528 + (-t525 * t559 + (-t421 + t831) * qJD(1)) * t525) + t728 * t478 * t457) * t825 - t528 * ((t528 * t348 + (t245 - t834) * qJD(1)) * t528 + (t246 * qJD(1) + (t411 * t719 + t413 * t721 + t726) * t525 + (-t349 + (t524 * t829 + t527 * t830) * qJD(2) + t577 * qJD(1)) * t528) * t525) - t528 * ((t528 * t337 + (t239 + t833) * qJD(1)) * t528 + (t238 * qJD(1) + (-t407 * t719 - t410 * t721 + t727) * t525 + (-t336 + (t406 * t527 + t409 * t524) * qJD(2) - t578 * qJD(1)) * t528) * t525) + t525 * ((t525 * t336 + (t241 + t835) * qJD(1)) * t525 + (t242 * qJD(1) + (t406 * t719 + t409 * t721) * t528 + (-t337 + (-t407 * t527 - t410 * t524) * qJD(2) + (t404 - t579) * qJD(1)) * t525) * t528) + t525 * t28 + t525 * t24 + t525 * t23 - t528 * t27 - t528 * t21 - t528 * t22 + (t108 + (-t238 - t246) * t528 + (t239 + t245) * t525 + t848) * t724 + (t107 + (-t241 - t244) * t528 + (t242 + t243) * t525 + t849) * t722; 0.2e1 * (t593 * t817 + (t213 * t528 + t214 * t525) * t818 + (t236 * t528 + t237 * t525) * t819 + (t333 * t528 + t334 * t525) * t820) * t719 + 0.2e1 * ((-t199 * t724 + t200 * t722 + t604) * t817 + (t117 * t525 + t118 * t528 - t213 * t724 + t214 * t722) * t818 + (t142 * t525 + t143 * t528 - t236 * t724 + t237 * t722) * t819 + (t225 * t525 + t226 * t528 - t333 * t724 + t334 * t722) * t820) * t524; 0.2e1 * ((t210 * t720 + t211 * t718 - t30) * t817 + (t254 * t720 + t255 * t718 - t35) * t818 + (t317 * t720 + t318 * t718 - t62) * t819 + (t379 * t720 + t380 * t718 - t116) * t820) * t527 + 0.2e1 * ((t210 * t722 - t211 * t724 + t545) * t817 + (qJD(2) * t124 + t140 * t528 + t141 * t525 + t254 * t722 - t255 * t724) * t818 + (qJD(2) * t160 + t193 * t528 + t194 * t525 + t317 * t722 - t318 * t724) * t819 + (qJD(2) * t231 + t256 * t528 + t257 * t525 + t379 * t722 - t380 * t724) * t820) * t524; 0.4e1 * (t820 + t819 + t818 + t817) * t373; t530 + (t702 * t528 + t701 * t525 + (t525 * t677 + t528 * t676) * qJD(1)) * t527 + m(5) * (t122 * t237 + t123 * t236 + t142 * t251 + t143 * t250) + m(6) * (t117 * t162 + t118 * t161 + t213 * t61 + t214 * t60) + m(7) * (t111 * t121 + t112 * t120 + t199 * t33 + t200 * t32) + (-t525 * t676 + t528 * t677) * t721 + t760; t532 + ((t167 * t525 - t168 * t528) * t837 + t28 * t813 + t27 * t815 + (t108 * t813 - t525 * t107 / 0.2e1) * qJD(1)) * t527 + (t76 * t806 + t107 * t672 - t15 / 0.2e1 + (qJD(1) * t167 - t59) * t816) * t528 + (t77 * t806 + t108 * t672 + t16 / 0.2e1 + (qJD(1) * t168 + t58) * t816) * t525 + m(5) * (t122 * t317 + t123 * t318 + t160 * t86 + t193 * t250 + t194 * t251 - t220 * t62) + m(6) * (t124 * t34 + t126 * t35 + t140 * t161 + t141 * t162 + t254 * t60 + t255 * t61) + m(7) * (t102 * t30 + t113 * t120 + t114 * t121 + t210 * t32 + t211 * t33 + t29 * t95); 0.2e1 * ((t250 * t718 + t251 * t720 - t86) * t819 + (t161 * t718 + t162 * t720 - t34) * t818 + (t120 * t718 + t121 * t720 - t29) * t817) * t527 + 0.2e1 * ((-qJD(2) * t220 + t122 * t525 + t123 * t528 - t250 * t724 + t251 * t722) * t819 + (qJD(2) * t126 - t161 * t724 + t162 * t722 + t525 * t60 + t528 * t61) * t818 + (-t120 * t724 + t121 * t722 + t548) * t817) * t524; ((t575 * t528 + (-t667 - t850) * t525) * qJD(2) + t760) * t524 + t542 + (t122 * t251 + t123 * t250 - t220 * t86) * t823 + (t126 * t34 + t161 * t61 + t162 * t60) * t822 + (t102 * t29 + t120 * t33 + t121 * t32) * t821 + (t528 * t16 + t525 * t15 + t524 * (t525 * t59 + t528 * t58) + (t240 * t524 + (t167 * t528 + t168 * t525) * t527) * qJD(2) + (t525 * t575 + t528 * t667) * qJD(1)) * t527; t530 + m(6) * (t109 * t214 + t110 * t213 + t117 * t230 + t118 * t229) + m(7) * (t111 * t137 + t112 * t136 + t199 * t45 + t200 * t44); t532 + m(6) * (t109 * t254 + t110 * t255 + t124 * t73 + t140 * t229 + t141 * t230 + t209 * t35) + m(7) * (t113 * t136 + t114 * t137 + t119 * t30 + t210 * t44 + t211 * t45 + t31 * t95); 0.2e1 * ((t229 * t718 + t230 * t720 - t73) * t818 + (t136 * t718 + t137 * t720 - t31) * t817) * t527 + 0.2e1 * ((qJD(2) * t209 + t109 * t525 + t110 * t528 - t229 * t724 + t230 * t722) * t818 + (-t136 * t724 + t137 * t722 + t547) * t817) * t524; t531 + m(7) * (t102 * t31 + t119 * t29 + t120 * t45 + t121 * t44 + t136 * t33 + t137 * t32) + m(6) * (t109 * t162 + t110 * t161 + t126 * t73 + t209 * t34 + t229 * t61 + t230 * t60); t531 + (t119 * t31 + t136 * t45 + t137 * t44) * t821 + (t109 * t230 + t110 * t229 + t209 * t73) * t822; m(7) * (-t593 * t721 + ((-t199 * t525 + t200 * t528) * qJD(1) + t604) * t527); m(7) * ((t30 + (-t210 * t525 - t211 * t528) * qJD(2)) * t524 + ((t210 * t528 - t211 * t525) * qJD(1) + t545) * t527); m(7) * (-t524 ^ 2 + t527 ^ 2) * t712 * qJD(2); m(7) * ((t29 + (-t120 * t528 - t121 * t525) * qJD(2)) * t524 + ((-t120 * t525 + t121 * t528) * qJD(1) + t548) * t527); m(7) * ((t31 + (-t136 * t528 - t137 * t525) * qJD(2)) * t524 + ((-t136 * t525 + t137 * t528) * qJD(1) + t547) * t527); -0.2e1 * m(7) * t373;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
