% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:23
% EndTime: 2019-12-31 20:03:54
% DurationCPUTime: 23.92s
% Computational Cost: add. (28614->758), mult. (70188->975), div. (0->0), fcn. (77891->6), ass. (0->455)
t524 = cos(qJ(2));
t827 = sin(qJ(4));
t651 = t524 * t827;
t522 = sin(qJ(2));
t828 = cos(qJ(4));
t653 = t522 * t828;
t462 = t653 - t651;
t652 = t522 * t827;
t557 = t524 * t828 + t652;
t343 = -Icges(6,5) * t557 - Icges(6,6) * t462;
t523 = sin(qJ(1));
t617 = t523 * t651;
t618 = t523 * t653;
t431 = t617 - t618;
t432 = t557 * t523;
t525 = cos(qJ(1));
t457 = Icges(6,4) * t462;
t350 = -Icges(6,2) * t557 + t457;
t357 = Icges(6,1) * t557 + t457;
t705 = t350 + t357;
t456 = Icges(6,4) * t557;
t351 = Icges(6,2) * t462 + t456;
t356 = Icges(6,1) * t462 - t456;
t707 = t351 - t356;
t153 = t343 * t525 + t707 * t431 - t705 * t432;
t346 = -Icges(5,5) * t557 - Icges(5,6) * t462;
t459 = Icges(5,4) * t462;
t353 = -Icges(5,2) * t557 + t459;
t360 = Icges(5,1) * t557 + t459;
t701 = t353 + t360;
t458 = Icges(5,4) * t557;
t354 = Icges(5,2) * t462 + t458;
t359 = Icges(5,1) * t462 - t458;
t703 = t354 - t359;
t154 = t346 * t525 + t703 * t431 - t701 * t432;
t977 = -t153 - t154;
t976 = t557 * (t359 / 0.2e1 + t356 / 0.2e1 - t354 / 0.2e1 - t351 / 0.2e1) + t462 * (t360 / 0.2e1 + t357 / 0.2e1 + t350 / 0.2e1 + t353 / 0.2e1);
t615 = t525 * t651;
t616 = t525 * t653;
t433 = -t616 + t615;
t434 = t557 * t525;
t157 = -t343 * t523 + t707 * t433 - t705 * t434;
t158 = -t346 * t523 + t703 * t433 - t701 * t434;
t396 = Icges(6,4) * t434;
t300 = -Icges(6,2) * t433 - Icges(6,6) * t523 + t396;
t722 = Icges(6,1) * t433 + t300 + t396;
t395 = Icges(6,4) * t433;
t304 = Icges(6,1) * t434 - Icges(6,5) * t523 - t395;
t877 = -Icges(6,2) * t434 + t304 - t395;
t164 = t722 * t462 + t557 * t877;
t400 = Icges(5,4) * t434;
t302 = -Icges(5,2) * t433 - Icges(5,6) * t523 + t400;
t720 = Icges(5,1) * t433 + t302 + t400;
t399 = Icges(5,4) * t433;
t306 = Icges(5,1) * t434 - Icges(5,5) * t523 - t399;
t876 = -Icges(5,2) * t434 + t306 - t399;
t166 = t720 * t462 + t557 * t876;
t975 = t157 + t158 - t164 - t166;
t974 = t153 / 0.2e1 + t154 / 0.2e1;
t973 = t153 / 0.4e1 + t154 / 0.4e1;
t972 = t164 / 0.2e1 + t166 / 0.2e1 - t157 / 0.2e1 - t158 / 0.2e1;
t394 = Icges(6,4) * t432;
t761 = Icges(6,2) * t431;
t598 = -t761 + t394;
t757 = Icges(6,6) * t525;
t551 = t757 + t598;
t393 = Icges(6,4) * t431;
t772 = Icges(6,1) * t432;
t604 = -t393 + t772;
t766 = Icges(6,5) * t525;
t554 = t766 + t604;
t581 = t431 * t300 - t432 * t304;
t592 = Icges(6,5) * t432 - Icges(6,6) * t431;
t547 = Icges(6,3) * t525 + t592;
t885 = t547 - t592;
t520 = t525 ^ 2;
t669 = 0.2e1 * t394;
t597 = t669 + 0.2e1 * t757;
t667 = 0.2e1 * t766;
t868 = t431 * (t597 - t761) - t432 * (t667 + t772) - Icges(6,3) * t520;
t948 = t868 * t525;
t25 = -t948 + (t885 * t523 + (-t554 + t604) * t434 + (t551 - t598) * t433 + t581) * t523;
t398 = Icges(5,4) * t432;
t764 = Icges(5,2) * t431;
t600 = -t764 + t398;
t758 = Icges(5,6) * t525;
t553 = t758 + t600;
t397 = Icges(5,4) * t431;
t775 = Icges(5,1) * t432;
t607 = -t397 + t775;
t767 = Icges(5,5) * t525;
t555 = t767 + t607;
t580 = t431 * t302 - t432 * t306;
t594 = Icges(5,5) * t432 - Icges(5,6) * t431;
t549 = Icges(5,3) * t525 + t594;
t884 = t549 - t594;
t671 = 0.2e1 * t398;
t599 = t671 + 0.2e1 * t758;
t668 = 0.2e1 * t767;
t867 = t431 * (t599 - t764) - t432 * (t668 + t775) - Icges(5,3) * t520;
t949 = t867 * t525;
t26 = -t949 + (t884 * t523 + (-t555 + t607) * t434 + (t553 - t600) * t433 + t580) * t523;
t546 = Icges(6,5) * t434 - Icges(6,6) * t433 - Icges(6,3) * t523;
t544 = t525 * t546;
t178 = t544 - t581;
t180 = -t433 * t300 + t434 * t304 - t523 * t546;
t27 = (t433 * t598 - t434 * t604 + t178 - 0.2e1 * t544 + t581) * t525 + (t431 * t551 - t432 * t554 - t525 * t885 + t180 - t868) * t523;
t548 = Icges(5,5) * t434 - Icges(5,6) * t433 - Icges(5,3) * t523;
t545 = t525 * t548;
t179 = t545 - t580;
t181 = -t433 * t302 + t434 * t306 - t523 * t548;
t28 = (t433 * t600 - t434 * t607 + t179 - 0.2e1 * t545 + t580) * t525 + (t431 * t553 - t432 * t555 - t525 * t884 + t181 - t867) * t523;
t830 = t525 / 0.4e1;
t832 = -t525 / 0.4e1;
t836 = -t523 / 0.4e1;
t94 = t178 * t523 + t948;
t95 = t179 * t523 + t949;
t96 = t180 * t523 - (-t433 * t551 + t434 * t554 - t523 * t547) * t525;
t97 = t181 * t523 - (-t433 * t553 + t434 * t555 - t523 * t549) * t525;
t971 = 0.2e1 * (t27 + t28) * t830 + 0.2e1 * (t96 + t97) * t832 + 0.2e1 * (t25 + t26 + t94 + t95) * t836 + (t166 / 0.4e1 + t164 / 0.4e1 - t158 / 0.4e1 - t157 / 0.4e1) * t523;
t509 = Icges(4,5) * t522;
t776 = Icges(4,1) * t524;
t608 = t509 + t776;
t410 = Icges(4,4) * t523 + t525 * t608;
t769 = Icges(3,4) * t522;
t475 = Icges(3,1) * t524 - t769;
t412 = Icges(3,5) * t523 + t475 * t525;
t968 = t410 + t412;
t752 = qJ(3) * t524;
t476 = pkin(2) * t522 - t752;
t793 = t522 * pkin(3);
t642 = t476 + t793;
t504 = t828 * pkin(4) + pkin(3);
t791 = pkin(3) - t504;
t644 = t522 * t791;
t698 = rSges(6,1) * t462 - rSges(6,2) * t557 - pkin(4) * t651 - t644;
t577 = t642 + t698;
t255 = t577 * t523;
t257 = t577 * t525;
t364 = rSges(5,1) * t462 - rSges(5,2) * t557;
t614 = t364 + t642;
t284 = t614 * t523;
t286 = t614 * t525;
t788 = rSges(4,1) * t522;
t477 = -rSges(4,3) * t524 + t788;
t683 = t476 + t477;
t373 = t683 * t523;
t375 = t683 * t525;
t516 = t525 * rSges(4,2);
t518 = t525 * pkin(6);
t779 = rSges(4,3) + qJ(3);
t829 = rSges(4,1) + pkin(2);
t869 = t779 * t522 + t829 * t524 + pkin(1);
t313 = -t523 * t869 + t516 + t518;
t314 = (rSges(4,2) + pkin(6)) * t523 + t869 * t525;
t740 = t524 * t525;
t742 = t523 * t524;
t724 = t313 * t740 + t314 * t742;
t753 = qJ(3) * t522;
t854 = pkin(2) + pkin(3);
t871 = t854 * t524 + pkin(1) + t753;
t881 = -t432 * rSges(5,1) + t431 * rSges(5,2);
t253 = t518 + (-rSges(5,3) - pkin(7)) * t525 - t871 * t523 + t881;
t312 = t434 * rSges(5,1) - t433 * rSges(5,2) - t523 * rSges(5,3);
t517 = t523 * pkin(7);
t254 = t523 * pkin(6) + t525 * t871 + t312 - t517;
t729 = t253 * t740 + t254 * t742;
t480 = pkin(2) * t524 + t753;
t521 = -qJ(5) - pkin(7);
t622 = pkin(4) * t652;
t943 = -t432 * rSges(6,1) + t431 * rSges(6,2) - t504 * t742 - t523 * t622;
t237 = t518 + (-rSges(6,3) + t521) * t525 + (-pkin(1) - t480) * t523 + t943;
t311 = t434 * rSges(6,1) - t433 * rSges(6,2) - t523 * rSges(6,3);
t792 = pkin(2) + t504;
t238 = (pkin(6) + t521) * t523 + (pkin(1) + t792 * t524 + (t827 * pkin(4) + qJ(3)) * t522) * t525 + t311;
t731 = t237 * t740 + t238 * t742;
t743 = t522 * t525;
t745 = t522 * t523;
t856 = m(6) / 0.2e1;
t858 = m(5) / 0.2e1;
t859 = m(4) / 0.2e1;
t663 = (-t284 * t743 + t286 * t745 + t729) * t858 + (-t373 * t743 + t375 * t745 + t724) * t859 + (-t255 * t743 + t257 * t745 + t731) * t856;
t500 = pkin(2) * t745;
t571 = t431 * rSges(6,1) + t432 * rSges(6,2) + pkin(4) * t617;
t259 = t500 + (t522 * t504 - t752) * t523 - t571;
t492 = qJ(3) * t740;
t883 = t433 * rSges(6,1) + t434 * rSges(6,2) + pkin(4) * t615;
t260 = -t792 * t743 + t492 + t883;
t367 = t500 + (-t779 * t524 + t788) * t523;
t499 = rSges(4,3) * t740;
t368 = -t829 * t743 + t492 + t499;
t328 = rSges(5,1) * t431 + rSges(5,2) * t432;
t273 = t500 + (-t752 + t793) * t523 - t328;
t332 = -t433 * rSges(5,1) - t434 * rSges(5,2);
t274 = -t854 * t743 - t332 + t492;
t586 = t273 * t525 + t274 * t523;
t664 = (t522 * t586 + t729) * t858 + ((t367 * t525 + t368 * t523) * t522 + t724) * t859 + ((t259 * t525 + t260 * t523) * t522 + t731) * t856;
t911 = t663 - t664;
t967 = t911 * qJD(1);
t737 = Icges(5,1) - Icges(5,2);
t640 = t737 * t432;
t563 = t640 + t767;
t902 = 0.2e1 * t397;
t945 = t563 - t902;
t964 = -t945 * t433 + t523 * (Icges(5,5) * t431 + Icges(5,6) * t432);
t736 = Icges(6,1) - Icges(6,2);
t638 = t736 * t432;
t562 = t638 + t766;
t903 = 0.2e1 * t393;
t946 = t562 - t903;
t963 = -t946 * t433 + t523 * (Icges(6,5) * t431 + Icges(6,6) * t432);
t961 = t976 * qJD(1);
t960 = m(6) * (-t237 * t525 - t238 * t523);
t233 = -t523 * (rSges(5,3) * t525 - t881) - t312 * t525;
t940 = -t523 * t328 + t332 * t525;
t959 = t233 * t940;
t315 = Icges(6,5) * t433 + Icges(6,6) * t434;
t956 = (-t315 * t523 + t433 * t877 + t722 * t434) * t523;
t317 = Icges(5,5) * t433 + Icges(5,6) * t434;
t955 = (-t317 * t523 + t433 * t876 + t720 * t434) * t523;
t954 = (t315 * t525 + t431 * t877 + t722 * t432) * t523;
t953 = (t317 * t525 + t431 * t876 + t720 * t432) * t523;
t952 = t940 * t524;
t519 = t523 ^ 2;
t679 = t519 + t520;
t685 = t679 * t480;
t941 = pkin(3) * t742 + pkin(7) * t525;
t620 = t523 * t941 + t525 * (pkin(3) * t740 - t517) + t685;
t190 = -t233 + t620;
t362 = -rSges(5,1) * t557 - rSges(5,2) * t462;
t585 = t284 * t523 + t286 * t525;
t951 = t190 * t940 + t362 * t585;
t950 = qJD(1) * t960;
t495 = Icges(4,5) * t740;
t402 = Icges(4,6) * t523 + Icges(4,3) * t743 + t495;
t468 = Icges(3,5) * t524 - Icges(3,6) * t522;
t404 = Icges(3,3) * t523 + t468 * t525;
t469 = Icges(4,4) * t524 + Icges(4,6) * t522;
t406 = Icges(4,2) * t523 + t469 * t525;
t947 = t402 * t743 + t968 * t740 + (t404 + t406) * t523;
t944 = (-Icges(3,6) + Icges(4,6)) * t524 + (-Icges(4,4) - Icges(3,5)) * t522;
t470 = Icges(3,2) * t524 + t769;
t756 = Icges(4,3) * t524;
t595 = t756 - t509;
t942 = (-t470 - t595) * t525 + t968;
t934 = -m(5) / 0.2e1;
t933 = -m(6) / 0.2e1;
t835 = t523 / 0.2e1;
t896 = -t525 / 0.2e1;
t923 = t362 * (t253 * t525 + t254 * t523);
t893 = t364 * t679;
t407 = Icges(3,4) * t742 - Icges(3,2) * t745 - Icges(3,6) * t525;
t512 = Icges(3,4) * t524;
t765 = Icges(3,2) * t522;
t408 = Icges(3,6) * t523 + (t512 - t765) * t525;
t379 = t412 * t742;
t624 = t404 * t525 - t379;
t403 = Icges(3,5) * t742 - Icges(3,6) * t745 - Icges(3,3) * t525;
t496 = Icges(3,4) * t745;
t411 = Icges(3,1) * t742 - Icges(3,5) * t525 - t496;
t696 = -t523 * t403 - t411 * t740;
t918 = -t407 * t743 - t408 * t745 - t624 - t696;
t917 = -t408 * t743 + t947;
t747 = (-Icges(4,2) * t525 + t523 * t469) * t525;
t914 = t747 + t947;
t637 = t736 * t431;
t913 = -t637 - t757;
t639 = t737 * t431;
t912 = -t639 - t758;
t392 = -t791 * t524 + t622;
t171 = -(rSges(6,3) * t525 - t941 - t943) * t523 - (t392 * t525 + t311 + t517) * t525;
t709 = t524 * t893;
t280 = t698 * t523;
t282 = t698 * t525;
t726 = t280 * t742 + t282 * t740;
t888 = t679 * t522;
t735 = (t171 * t888 + t726) * t856 + (t233 * t888 + t709) * t858;
t203 = -(t523 * t644 + t571) * t523 - (t525 * t644 + t883) * t525;
t361 = -rSges(6,1) * t557 - rSges(6,2) * t462;
t697 = -t361 + t392;
t281 = t697 * t523;
t283 = t697 * t525;
t629 = t679 * t362;
t778 = (-t952 + (t233 - t629) * t522 + t709) * t858 + (-t203 * t524 + (t281 * t523 + t283 * t525 + t171) * t522 + t726) * t856;
t909 = t735 - t778;
t908 = t362 * t888 + t952;
t279 = pkin(4) * t616 - t883;
t537 = pkin(4) * t618 - t571;
t208 = t279 * t525 + t523 * t537;
t837 = -t523 / 0.2e1;
t831 = t525 / 0.2e1;
t699 = -t557 * pkin(4) + t361;
t288 = t699 * t523;
t289 = t699 * t525;
t564 = m(6) * (t288 * t525 - t523 * t289);
t768 = Icges(4,5) * t524;
t467 = Icges(4,3) * t522 + t768;
t401 = -Icges(4,6) * t525 + t467 * t523;
t409 = -Icges(4,4) * t525 + t523 * t608;
t892 = t523 * (t401 * t522 + t409 * t524);
t744 = t522 * t524;
t681 = t679 * t744;
t855 = m(6) / 0.4e1;
t857 = m(5) / 0.4e1;
t889 = (m(4) / 0.4e1 + t857 + t855) * (t681 - t744);
t880 = t944 * t523;
t879 = t944 * t525;
t878 = t942 * t523;
t732 = (m(6) * (t279 * t523 - t525 * t537) + m(5) * (t328 * t525 + t332 * t523)) * t522 / 0.2e1;
t872 = t524 * t208 + (t288 * t523 + t289 * t525) * t522;
t733 = t872 * t933 + t908 * t934;
t588 = -t237 * t289 - t238 * t288;
t875 = (-t255 * t279 + t257 * t537 + t588) * t933 + (-t284 * t332 - t286 * t328 - t923) * t934;
t874 = (t237 * t283 + t238 * t281 + t259 * t282 + t260 * t280) * t933 + (t364 * t586 - t923) * t934;
t49 = t954 - ((-t903 + t667) * t431 + (t597 + 0.2e1 * t637) * t432) * t525;
t50 = t953 - ((-t902 + t668) * t431 + (t599 + 0.2e1 * t639) * t432) * t525;
t542 = t669 - t913;
t53 = t956 - (t542 * t434 - t963) * t525;
t543 = t671 - t912;
t54 = t955 - (t543 * t434 - t964) * t525;
t534 = (t53 + t54) * t837 + (t49 + t50) * t831;
t472 = Icges(4,1) * t522 - t768;
t777 = Icges(3,1) * t522;
t870 = -t522 * (t475 / 0.2e1 - t470 / 0.2e1 + t509 + t776 / 0.2e1 - t756 / 0.2e1) - t524 * (t512 + t777 / 0.2e1 - t765 / 0.2e1 + t472 / 0.2e1 - t467 / 0.2e1);
t625 = -0.2e1 * t394;
t51 = -t954 - ((t625 - 0.2e1 * t757) * t432 + (-0.2e1 * t638 + t903 - 0.2e1 * t766) * t431) * t525;
t627 = -0.2e1 * t398;
t52 = -t953 - ((t627 - 0.2e1 * t758) * t432 + (-0.2e1 * t640 + t902 - 0.2e1 * t767) * t431) * t525;
t538 = t625 + t913;
t55 = -t956 - (t538 * t434 + t963) * t525;
t539 = t627 + t912;
t56 = -t955 - (t539 * t434 + t964) * t525;
t533 = (t51 + t52) * t896 + (t55 + t56) * t835;
t445 = t472 * t523;
t446 = -Icges(4,1) * t743 + t495;
t609 = -t512 - t777;
t447 = t609 * t523;
t448 = t609 * t525;
t866 = ((t407 - t447 - t401 + t445) * t525 + (t402 + t446 - t408 + t448) * t523) * t524;
t864 = 0.4e1 * qJD(1);
t863 = 2 * qJD(2);
t862 = 4 * qJD(2);
t861 = 2 * qJD(4);
t860 = 4 * qJD(4);
t851 = m(5) * t951;
t130 = -t171 + t620;
t846 = m(6) * (t130 * t208 + t255 * t288 + t257 * t289);
t842 = m(6) * (t171 * t203 + t280 * t281 + t282 * t283);
t728 = -t255 * t742 - t257 * t740;
t841 = m(6) * (t130 * t888 + t728);
t789 = rSges(3,1) * t524;
t647 = pkin(1) + t789;
t680 = rSges(3,2) * t745 + t525 * rSges(3,3);
t371 = -t523 * t647 + t518 + t680;
t498 = rSges(3,2) * t743;
t372 = -t498 + t647 * t525 + (rSges(3,3) + pkin(6)) * t523;
t478 = rSges(3,1) * t522 + rSges(3,2) * t524;
t453 = t478 * t523;
t455 = t478 * t525;
t826 = m(3) * (t371 * t453 - t372 * t455);
t481 = rSges(4,1) * t524 + rSges(4,3) * t522;
t250 = t523 * (t481 * t523 - t516) + (t523 * rSges(4,2) + t481 * t525) * t525 + t685;
t708 = -t373 * t742 - t375 * t740;
t822 = m(4) * (t250 * t888 + t708);
t820 = m(4) * (t313 * t367 + t314 * t368);
t819 = m(4) * (-t313 * t745 + t314 * t743);
t817 = m(5) * (-t364 * t629 + t959);
t725 = -t284 * t742 - t286 * t740;
t815 = m(5) * (t190 * t888 + t725);
t812 = m(5) * (t253 * t273 + t254 * t274);
t811 = m(5) * (t253 * t328 + t254 * t332);
t809 = m(5) * (-t253 * t745 + t254 * t743);
t805 = m(6) * (t237 * t259 + t238 * t260);
t802 = m(6) * (-t237 * t537 + t238 * t279);
t800 = m(6) * (-t237 * t745 + t238 * t743);
t799 = m(6) * (t523 * t255 + t257 * t525);
t798 = m(6) * (-t523 * t259 + t260 * t525);
t796 = m(6) * (t280 * t743 - t282 * t745);
t795 = m(6) * t208;
t794 = m(6) * (-t280 * t523 - t282 * t525);
t746 = t407 * t522;
t692 = -t595 * t523 + t409;
t690 = -Icges(3,2) * t742 + t411 - t496;
t686 = t523 * (qJ(3) * t742 - t500) + t525 * (-pkin(2) * t743 + t492);
t682 = -t480 - t481;
t675 = qJD(2) * t522;
t556 = (t281 * t525 - t523 * t283) * t856;
t132 = t564 / 0.2e1 + t556;
t674 = t132 * qJD(5);
t337 = m(6) * t888;
t673 = t337 * qJD(1);
t648 = t469 / 0.2e1 + t468 / 0.2e1;
t641 = -t524 * pkin(3) - t480;
t634 = t692 * t525;
t632 = t690 * t525;
t623 = t408 * t522 - t403;
t613 = t362 + t641;
t612 = -t402 * t745 + t406 * t525 - t410 * t742;
t285 = t613 * t523;
t287 = t613 * t525;
t584 = t285 * t523 + t287 * t525;
t576 = t641 - t697;
t573 = t27 / 0.2e1 - t96 / 0.2e1 + t28 / 0.2e1 - t97 / 0.2e1;
t572 = t94 / 0.2e1 + t25 / 0.2e1 + t95 / 0.2e1 + t26 / 0.2e1;
t559 = -pkin(3) * t888 + t686;
t558 = -t533 + t534;
t163 = t542 * t462 + t557 * t946;
t165 = t543 * t462 + t557 * t945;
t541 = -t874 - t975 * t836 + (t163 + t165 + t977) * t830;
t159 = t538 * t462 - (t562 - 0.2e1 * t393) * t557;
t161 = t539 * t462 - (t563 - 0.2e1 * t397) * t557;
t540 = -t875 + t975 * t523 / 0.4e1 + (t159 + t161 - t977) * t832;
t261 = -t747 + t892;
t536 = t612 * t837 + t572 + (-t523 * (-t411 * t524 + t746) + t261 - t403 * t525) * t896 + (t525 * t623 - t914 + t917) * t831 + (t401 * t743 + t409 * t740 + t523 * t623 + t612 + t624 + t918) * t835;
t535 = -t573 + (t261 - t892 + t914) * t837 + t917 * t835 + (-t379 + (t404 + t746) * t525 + t696 + t918) * t896;
t528 = t523 * t572 - t525 * t573;
t482 = -rSges(3,2) * t522 + t789;
t376 = t682 * t525;
t374 = t682 * t523;
t269 = t525 * (-rSges(4,1) * t743 + t499) - t477 * t519 + t686;
t258 = t576 * t525;
t256 = t576 * t523;
t243 = 0.4e1 * t889;
t219 = t794 / 0.2e1;
t215 = t795 / 0.2e1;
t209 = t796 / 0.2e1;
t206 = -t940 + t559;
t189 = t798 / 0.2e1;
t186 = t799 / 0.2e1;
t167 = -t203 + t559;
t133 = t556 - t564 / 0.2e1;
t81 = t219 - t795 / 0.2e1;
t80 = t219 + t215;
t79 = t215 - t794 / 0.2e1;
t66 = t800 + t809 + t819;
t65 = t186 - t798 / 0.2e1;
t64 = t186 + t189;
t63 = t189 - t799 / 0.2e1;
t34 = t209 - t732;
t33 = -t796 / 0.2e1 + t732;
t32 = t209 + t732;
t30 = t815 + t822 + t841;
t29 = t802 + t811 - t976;
t16 = t805 + t812 + t820 + t826 - t870 + t976;
t13 = t735 + t778 - t733;
t12 = t733 - t909;
t11 = t733 + t909;
t9 = t663 + t664;
t7 = (-t49 / 0.2e1 - t50 / 0.2e1) * t525 + (t53 / 0.2e1 + t54 / 0.2e1) * t523 + t817 + t842;
t6 = (-t52 / 0.2e1 - t51 / 0.2e1) * t525 + (t56 / 0.2e1 + t55 / 0.2e1) * t523 + t851 + t846;
t4 = t523 * t536 + t525 * t535;
t3 = t540 + (-t165 / 0.4e1 - t163 / 0.4e1 + t973) * t525 + t874 + t971;
t2 = t541 + (t161 / 0.4e1 + t159 / 0.4e1 + t973) * t525 + t875 + t971;
t1 = t528 + t541 + t540;
t5 = [t16 * qJD(2) + t66 * qJD(3) + t29 * qJD(4) + qJD(5) * t960, t16 * qJD(1) + t9 * qJD(3) + t1 * qJD(4) + t64 * qJD(5) + ((t313 * t376 + t314 * t374 - t367 * t375 - t368 * t373) * t859 + (t253 * t287 + t254 * t285 - t273 * t286 - t274 * t284) * t858 + (t237 * t258 + t238 * t256 - t255 * t260 - t257 * t259) * t856) * t863 + ((m(3) * (-t371 * t482 - t453 * t478) - t165 / 0.2e1 - t163 / 0.2e1 + t648 * t525 - t535 + t974) * qJD(2) + (-t401 / 0.2e1 + t445 / 0.2e1 + t407 / 0.2e1 - t447 / 0.2e1) * t675) * t525 + ((m(3) * (-t372 * t482 + t455 * t478) + t648 * t523 - t536 + t972) * qJD(2) + (t402 / 0.2e1 + t446 / 0.2e1 - t408 / 0.2e1 + t448 / 0.2e1) * t675) * t523 + (-t634 / 0.2e1 - t632 / 0.2e1 + t942 * t835) * qJD(2) * t524, qJD(1) * t66 + qJD(2) * t9 + qJD(4) * t32, t29 * qJD(1) + t1 * qJD(2) + t32 * qJD(3) + t80 * qJD(5) + (m(6) * (t279 * t280 - t282 * t537 - t588) + (m(5) * (t253 * t362 + t328 * t364) + t161 / 0.2e1 + t159 / 0.2e1 + t573 + t974) * t525 + (m(5) * (t254 * t362 + t332 * t364) - t572 + t972) * t523) * qJD(4), t64 * qJD(2) + t80 * qJD(4) + t950; t4 * qJD(2) + t911 * qJD(3) + t3 * qJD(4) + t65 * qJD(5) + (-t826 / 0.4e1 - t820 / 0.4e1 - t812 / 0.4e1 - t805 / 0.4e1) * t864 + t870 * qJD(1) - t961, t4 * qJD(1) + (m(4) * (t250 * t269 - t373 * t374 - t375 * t376) + m(3) * ((t523 * (rSges(3,1) * t742 - t680) + t525 * (rSges(3,1) * t740 + t523 * rSges(3,3) - t498)) * (-t523 * t453 - t455 * t525) + t679 * t482 * t478) + m(6) * (t130 * t167 - t255 * t256 - t257 * t258) + m(5) * (t190 * t206 - t284 * t285 - t286 * t287) - t534 + ((-t880 * t523 + (t634 + t632 - t878) * t522 + t866) * t525 + t879 * t519) * t835 + ((-t879 * t525 + t866 + ((t690 + t692) * t525 - t878) * t522) * t523 + t880 * t520) * t896) * qJD(2) + t30 * qJD(3) + t6 * qJD(4), t967 + t30 * qJD(2) + t11 * qJD(4) + (-0.4e1 * t889 + 0.2e1 * (t856 + t858 + t859) * (-t524 * t888 + t681)) * qJD(3), t3 * qJD(1) + t6 * qJD(2) + t11 * qJD(3) + (-t533 + t558) * qJD(4) - t674 + (-t817 / 0.4e1 - t842 / 0.4e1) * t860 + (((-t257 - t282) * t289 + (-t255 - t280) * t288 + (-t130 + t171) * t208) * t856 + ((-t190 + t233) * t940 + (-t585 - t893) * t362) * t858) * t861, qJD(1) * t65 - qJD(4) * t132; -t911 * qJD(2) + t33 * qJD(4) - t337 * qJD(5) + (-t800 / 0.4e1 - t809 / 0.4e1 - t819 / 0.4e1) * t864, -t967 + t243 * qJD(3) + t12 * qJD(4) + (-t841 / 0.4e1 - t815 / 0.4e1 - t822 / 0.4e1) * t862 + ((-t524 * t167 + t728) * t856 + (-t524 * t206 + t725) * t858 + (-t524 * t269 + t708) * t859 + ((t256 * t523 + t258 * t525 + t130) * t856 + (t190 + t584) * t858 + (t374 * t523 + t376 * t525 + t250) * t859) * t522) * t863, t243 * qJD(2), t33 * qJD(1) + t12 * qJD(2) + (t872 * t856 + t908 * t858) * t861, -t673; t2 * qJD(2) + t34 * qJD(3) + t528 * qJD(4) + t81 * qJD(5) + (-t802 / 0.4e1 - t811 / 0.4e1) * t864 + t961, t2 * qJD(1) + (t534 + t558) * qJD(2) + t13 * qJD(3) + t7 * qJD(4) + t674 + (-t846 / 0.4e1 - t851 / 0.4e1) * t862 + ((t130 * t203 + t167 * t171 - t255 * t281 + t256 * t280 - t257 * t283 + t258 * t282) * t856 + (t233 * t206 + t364 * t584 + t951) * t858) * t863, qJD(1) * t34 + qJD(2) * t13, t528 * qJD(1) + t7 * qJD(2) + t533 * qJD(4) + ((t362 * t893 - t959) * t857 + (-t171 * t208 + t280 * t288 + t282 * t289) * t855) * t860, qJD(1) * t81 + qJD(2) * t132; t63 * qJD(2) + t337 * qJD(3) + t79 * qJD(4) - t950, t63 * qJD(1) + m(6) * (t256 * t525 - t523 * t258) * qJD(2) + t133 * qJD(4), t673, t79 * qJD(1) + t133 * qJD(2) + qJD(4) * t564, 0;];
Cq = t5;
