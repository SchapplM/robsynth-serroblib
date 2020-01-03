% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP10_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP10_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:57
% EndTime: 2019-12-31 18:51:43
% DurationCPUTime: 36.58s
% Computational Cost: add. (65902->816), mult. (83266->1128), div. (0->0), fcn. (90442->7), ass. (0->478)
t546 = pkin(8) + qJ(3);
t541 = cos(t546);
t553 = cos(qJ(4));
t554 = cos(qJ(1));
t740 = t553 * t554;
t551 = sin(qJ(4));
t552 = sin(qJ(1));
t742 = t552 * t551;
t501 = t541 * t742 + t740;
t737 = t554 * t551;
t741 = t552 * t553;
t502 = t541 * t741 - t737;
t540 = sin(t546);
t747 = t540 * t552;
t353 = Icges(6,5) * t502 - Icges(6,6) * t501 + Icges(6,3) * t747;
t356 = Icges(5,5) * t502 - Icges(5,6) * t501 + Icges(5,3) * t747;
t483 = Icges(6,4) * t502;
t359 = -Icges(6,2) * t501 + Icges(6,6) * t747 + t483;
t486 = Icges(5,4) * t502;
t362 = -Icges(5,2) * t501 + Icges(5,6) * t747 + t486;
t482 = Icges(6,4) * t501;
t366 = -Icges(6,1) * t502 - Icges(6,5) * t747 + t482;
t485 = Icges(5,4) * t501;
t369 = -Icges(5,1) * t502 - Icges(5,5) * t747 + t485;
t961 = (t353 + t356) * t747 + (-t366 - t369) * t502 + (-t359 - t362) * t501;
t503 = -t541 * t737 + t741;
t504 = t541 * t740 + t742;
t746 = t540 * t554;
t355 = Icges(6,5) * t504 + Icges(6,6) * t503 + Icges(6,3) * t746;
t358 = Icges(5,5) * t504 + Icges(5,6) * t503 + Icges(5,3) * t746;
t771 = Icges(6,4) * t504;
t361 = Icges(6,2) * t503 + Icges(6,6) * t746 + t771;
t774 = Icges(5,4) * t504;
t364 = Icges(5,2) * t503 + Icges(5,6) * t746 + t774;
t484 = Icges(6,4) * t503;
t367 = Icges(6,1) * t504 + Icges(6,5) * t746 + t484;
t487 = Icges(5,4) * t503;
t370 = Icges(5,1) * t504 + Icges(5,5) * t746 + t487;
t960 = (t355 + t358) * t747 + (t367 + t370) * t502 + (-t361 - t364) * t501;
t612 = Icges(6,5) * t553 - Icges(6,6) * t551;
t439 = -Icges(6,3) * t541 + t540 * t612;
t769 = Icges(6,4) * t553;
t615 = -Icges(6,2) * t551 + t769;
t443 = -Icges(6,6) * t541 + t540 * t615;
t770 = Icges(6,4) * t551;
t617 = Icges(6,1) * t553 - t770;
t447 = -Icges(6,5) * t541 + t540 * t617;
t260 = t439 * t747 - t443 * t501 + t447 * t502;
t613 = Icges(5,5) * t553 - Icges(5,6) * t551;
t441 = -Icges(5,3) * t541 + t540 * t613;
t772 = Icges(5,4) * t553;
t616 = -Icges(5,2) * t551 + t772;
t445 = -Icges(5,6) * t541 + t540 * t616;
t773 = Icges(5,4) * t551;
t618 = Icges(5,1) * t553 - t773;
t449 = -Icges(5,5) * t541 + t540 * t618;
t261 = t441 * t747 - t445 * t501 + t449 * t502;
t902 = t260 + t261;
t959 = t961 * t552 + t960 * t554;
t834 = t540 / 0.2e1;
t833 = -t541 / 0.2e1;
t404 = -rSges(5,1) * t501 - rSges(5,2) * t502;
t406 = rSges(5,1) * t503 - rSges(5,2) * t504;
t306 = t552 * t404 + t406 * t554;
t852 = m(6) / 0.2e1;
t856 = -m(5) / 0.2e1;
t933 = rSges(6,1) + pkin(4);
t337 = t502 * rSges(6,2) + t933 * t501;
t338 = -t504 * rSges(6,2) + t933 * t503;
t924 = t552 * t337 - t338 * t554;
t722 = t306 * t856 + t924 * t852;
t538 = pkin(4) * t553 + pkin(3);
t791 = pkin(3) - t538;
t880 = t791 * t540;
t783 = rSges(6,1) * t553;
t620 = -rSges(6,2) * t551 + t783;
t550 = -qJ(5) - pkin(7);
t789 = pkin(7) + t550;
t925 = t540 * t620 + (-rSges(6,3) + t789) * t541;
t700 = -t880 + t925;
t636 = t791 * t541;
t793 = pkin(7) * t540;
t875 = t504 * rSges(6,1) + t503 * rSges(6,2) + pkin(4) * t742 - t550 * t746;
t715 = (-t636 - t793) * t554 + rSges(6,3) * t746 + t875;
t226 = t541 * t715 + t700 * t746;
t624 = t504 * rSges(5,1) + t503 * rSges(5,2);
t376 = rSges(5,3) * t746 + t624;
t784 = rSges(5,1) * t553;
t623 = -rSges(5,2) * t551 + t784;
t578 = t623 * t540;
t452 = -rSges(5,3) * t541 + t578;
t304 = t541 * t376 + t452 * t746;
t876 = -t502 * rSges(5,1) + t501 * rSges(5,2);
t372 = rSges(5,3) * t747 - t876;
t934 = t372 * t541 + t452 * t747;
t879 = t304 * t552 - t554 * t934;
t792 = t541 * pkin(3);
t882 = t789 * t540;
t744 = t541 * t552;
t944 = -t502 * rSges(6,1) + t501 * rSges(6,2) - t538 * t744;
t716 = -pkin(4) * t737 - (t792 + t882) * t552 + rSges(6,3) * t747 - t944;
t951 = t541 * t716 + t700 * t747;
t736 = (-t226 * t552 + t554 * t951) * t852 + t856 * t879;
t15 = t736 - t722;
t958 = qJD(1) * t15;
t831 = t552 / 0.2e1;
t829 = -t554 / 0.2e1;
t515 = t792 + t793;
t547 = t552 ^ 2;
t548 = t554 ^ 2;
t679 = t547 + t548;
t687 = t679 * t515;
t173 = t552 * t716 + t554 * t715 + t687;
t248 = t552 * t372 + t376 * t554 + t687;
t514 = pkin(3) * t540 - pkin(7) * t541;
t656 = t514 + t700;
t311 = t656 * t552;
t313 = t656 * t554;
t794 = pkin(4) * t551;
t627 = (rSges(6,1) * t551 + rSges(6,2) * t553 + t794) * t540;
t409 = t627 * t552;
t410 = t627 * t554;
t491 = (-rSges(5,1) * t551 - rSges(5,2) * t553) * t540;
t692 = t452 + t514;
t377 = t692 * t552;
t379 = t692 * t554;
t873 = t377 * t552 + t379 * t554;
t957 = -m(5) * (t248 * t306 + t873 * t491) - m(6) * (-t173 * t924 - t311 * t409 - t313 * t410);
t664 = t902 * t833 + t959 * t834;
t264 = t439 * t746 + t503 * t443 + t504 * t447;
t217 = t353 * t746 + t503 * t359 - t504 * t366;
t218 = t355 * t746 + t503 * t361 + t504 * t367;
t609 = t552 * t217 + t218 * t554;
t945 = -t264 * t541 + t540 * t609;
t265 = t441 * t746 + t503 * t445 + t504 * t449;
t219 = t356 * t746 + t503 * t362 - t504 * t369;
t220 = t358 * t746 + t503 * t364 + t504 * t370;
t608 = t552 * t219 + t220 * t554;
t946 = -t265 * t541 + t540 * t608;
t663 = t946 / 0.2e1 + t945 / 0.2e1;
t389 = -Icges(6,5) * t501 - Icges(6,6) * t502;
t391 = -Icges(5,5) * t501 - Icges(5,6) * t502;
t953 = t389 + t391;
t390 = Icges(6,5) * t503 - Icges(6,6) * t504;
t392 = Icges(5,5) * t503 - Icges(5,6) * t504;
t952 = t390 + t392;
t707 = -Icges(5,2) * t504 + t370 + t487;
t709 = -Icges(6,2) * t504 + t367 + t484;
t950 = t707 + t709;
t708 = -Icges(5,2) * t502 - t369 - t485;
t710 = -Icges(6,2) * t502 - t366 - t482;
t949 = t708 + t710;
t711 = Icges(5,1) * t503 - t364 - t774;
t713 = Icges(6,1) * t503 - t361 - t771;
t948 = t711 + t713;
t712 = -Icges(5,1) * t501 - t362 - t486;
t714 = -Icges(6,1) * t501 - t359 - t483;
t947 = t712 + t714;
t755 = t356 * t541;
t921 = t362 * t551 + t369 * t553;
t241 = t921 * t540 + t755;
t757 = t353 * t541;
t922 = t359 * t551 + t366 * t553;
t238 = t922 * t540 + t757;
t778 = -rSges(6,3) + t550;
t335 = (t778 * t541 + (t538 + t620) * t540) * t552;
t661 = t540 * t737;
t743 = t541 * t554;
t682 = rSges(6,2) * t661 + rSges(6,3) * t743;
t745 = t541 * t550;
t336 = (-t745 + (-t538 - t783) * t540) * t554 + t682;
t530 = pkin(3) * t747;
t825 = rSges(5,3) + pkin(7);
t347 = t530 + (-t825 * t541 + t578) * t552;
t531 = pkin(7) * t743;
t681 = rSges(5,2) * t661 + rSges(5,3) * t743;
t348 = t531 + (-pkin(3) - t784) * t746 + t681;
t511 = rSges(4,1) * t540 + rSges(4,2) * t541;
t492 = t511 * t552;
t493 = t511 * t554;
t381 = t552 * t492 + t493 * t554;
t855 = m(5) / 0.2e1;
t657 = (t552 * t335 - t336 * t554) * t852 + (t552 * t347 - t348 * t554) * t855 + m(4) * t381 / 0.2e1;
t883 = t679 * t511;
t658 = (-t311 * t552 - t313 * t554) * t852 + t873 * t856 - m(4) * t883 / 0.2e1;
t59 = t658 - t657;
t942 = qJD(1) * t59;
t247 = (t337 * t554 + t338 * t552) * t540;
t938 = -t238 - t241;
t604 = -t361 * t551 + t367 * t553;
t756 = t355 * t541;
t239 = t540 * t604 - t756;
t602 = -t364 * t551 + t370 * t553;
t754 = t358 * t541;
t242 = t540 * t602 - t754;
t937 = t239 + t242;
t936 = -t217 * t554 + t218 * t552;
t935 = -t219 * t554 + t220 * t552;
t535 = cos(pkin(8)) * pkin(2) + pkin(1);
t790 = -pkin(6) - qJ(2);
t536 = t552 * t790;
t781 = rSges(6,3) * t540;
t290 = -t536 + (t538 * t541 + t535 + t781) * t554 + t875;
t898 = (t778 * t540 - t535) * t552 + (-t790 + t794) * t554 + t944;
t195 = t290 * t746 - t747 * t898;
t932 = t195 * m(6) * qJD(1);
t931 = -t949 * t501 + t947 * t502 + t953 * t747;
t930 = -t950 * t501 + t948 * t502 + t952 * t747;
t929 = t949 * t503 + t947 * t504 + t953 * t746;
t928 = t950 * t503 + t948 * t504 + t952 * t746;
t470 = (-Icges(6,5) * t551 - Icges(6,6) * t553) * t540;
t474 = (-Icges(6,2) * t553 - t770) * t540;
t694 = t447 + t474;
t478 = (-Icges(6,1) * t551 - t769) * t540;
t696 = -t443 + t478;
t207 = t470 * t747 - t501 * t694 + t502 * t696;
t471 = (-Icges(5,5) * t551 - Icges(5,6) * t553) * t540;
t475 = (-Icges(5,2) * t553 - t773) * t540;
t693 = t449 + t475;
t479 = (-Icges(5,1) * t551 - t772) * t540;
t695 = -t445 + t479;
t208 = t471 * t747 - t501 * t693 + t502 * t695;
t927 = -t208 - t207;
t209 = t470 * t746 + t503 * t694 + t504 * t696;
t210 = t471 * t746 + t503 * t693 + t504 * t695;
t926 = -t210 - t209;
t867 = t825 * t540 + t535 + t792;
t317 = t554 * t867 - t536 + t624;
t637 = t554 * t790;
t897 = -t552 * t867 - t637 + t876;
t878 = -t317 * t552 - t554 * t897;
t853 = -m(6) / 0.2e1;
t910 = t452 * t552;
t415 = t443 * t552;
t419 = t447 * t552;
t588 = -t439 * t552 + t922;
t164 = -t588 * t541 + (t415 * t551 - t419 * t553 + t353) * t540;
t417 = t445 * t552;
t421 = t449 * t552;
t586 = -t441 * t552 + t921;
t166 = -t586 * t541 + (t417 * t551 - t421 * t553 + t356) * t540;
t909 = t164 + t166;
t416 = t443 * t554;
t420 = t447 * t554;
t587 = -t439 * t554 - t604;
t165 = -t587 * t541 + (t416 * t551 - t420 * t553 + t355) * t540;
t418 = t445 * t554;
t422 = t449 * t554;
t585 = -t441 * t554 - t602;
t167 = -t585 * t541 + (t418 * t551 - t422 * t553 + t358) * t540;
t908 = t165 + t167;
t175 = -t389 * t541 + (-t551 * t710 + t553 * t714) * t540;
t177 = -t391 * t541 + (-t551 * t708 + t553 * t712) * t540;
t907 = t175 + t177;
t176 = -t390 * t541 + (-t551 * t709 + t553 * t713) * t540;
t178 = -t392 * t541 + (-t551 * t707 + t553 * t711) * t540;
t906 = t176 + t178;
t444 = Icges(6,6) * t540 + t541 * t615;
t448 = Icges(6,5) * t540 + t541 * t617;
t440 = Icges(6,3) * t540 + t541 * t612;
t600 = -t551 * t443 + t553 * t447;
t584 = t440 - t600;
t752 = t439 * t541;
t562 = t540 * t584 + t752;
t184 = -t444 * t501 + t448 * t502 + t552 * t562;
t446 = Icges(5,6) * t540 + t541 * t616;
t450 = Icges(5,5) * t540 + t541 * t618;
t442 = Icges(5,3) * t540 + t541 * t613;
t599 = -t551 * t445 + t553 * t449;
t583 = t442 - t599;
t751 = t441 * t541;
t561 = t540 * t583 + t751;
t185 = -t446 * t501 + t450 * t502 + t552 * t561;
t905 = t184 + t185;
t186 = t503 * t444 + t504 * t448 + t554 * t562;
t187 = t503 * t446 + t504 * t450 + t554 * t561;
t904 = t186 + t187;
t196 = -t584 * t541 + (-t551 * t444 + t553 * t448 + t439) * t540;
t197 = -t583 * t541 + (-t551 * t446 + t553 * t450 + t441) * t540;
t903 = t196 + t197;
t901 = t264 + t265;
t282 = t540 * t600 - t752;
t283 = t540 * t599 - t751;
t900 = t282 + t283;
t457 = Icges(4,4) * t744 - Icges(4,2) * t747 - Icges(4,6) * t554;
t534 = Icges(4,4) * t541;
t766 = Icges(4,2) * t540;
t458 = Icges(4,6) * t552 + (t534 - t766) * t554;
t775 = Icges(4,4) * t540;
t510 = Icges(4,1) * t541 - t775;
t460 = Icges(4,5) * t552 + t510 * t554;
t431 = t460 * t744;
t506 = Icges(4,5) * t541 - Icges(4,6) * t540;
t456 = Icges(4,3) * t552 + t506 * t554;
t629 = t554 * t456 - t431;
t455 = Icges(4,5) * t744 - Icges(4,6) * t747 - Icges(4,3) * t554;
t525 = Icges(4,4) * t747;
t459 = Icges(4,1) * t744 - Icges(4,5) * t554 - t525;
t698 = -t552 * t455 - t459 * t743;
t899 = -t457 * t746 - t458 * t747 - t629 - t698;
t896 = t938 * t552 + t937 * t554;
t894 = -t540 / 0.2e1;
t893 = t541 / 0.2e1;
t892 = -t552 / 0.2e1;
t827 = t554 / 0.2e1;
t889 = t927 * t541 + (t931 * t552 + t930 * t554) * t540;
t888 = t926 * t541 + (t929 * t552 + t928 * t554) * t540;
t872 = t164 / 0.2e1 + t166 / 0.2e1;
t871 = -t165 / 0.2e1 - t167 / 0.2e1;
t870 = (-t377 * t406 + t379 * t404 + t878 * t491) * t856 + (t290 * t409 - t311 * t338 - t313 * t337 + t410 * t898) * t853;
t699 = t541 * t620 - t636 + t781 - t882;
t704 = t530 + (-t538 * t540 - t925) * t552;
t135 = (t552 * t700 + t704) * t541 + (t552 * t699 - t716) * t540;
t660 = t540 * t740;
t703 = -t531 + (-t745 + t880) * t554 - rSges(6,1) * t660 + t682;
t136 = (-t554 * t700 - t703) * t541 + (-t554 * t699 + t715) * t540;
t454 = rSges(5,3) * t540 + t541 * t623;
t243 = (t454 * t552 - t372) * t540;
t428 = -rSges(5,1) * t660 + t681;
t244 = (-t452 * t554 - t428) * t541 + (-t454 * t554 + t376) * t540;
t869 = (t135 * t898 + t136 * t290 - t226 * t336 + t335 * t951) * t853 + (t243 * t897 + t244 * t317 - t304 * t348 + t347 * t934) * t856;
t641 = -t449 / 0.2e1 - t447 / 0.2e1;
t643 = t445 / 0.2e1 + t443 / 0.2e1;
t866 = -t551 * (t475 / 0.2e1 + t474 / 0.2e1 - t641) - t553 * (-t479 / 0.2e1 - t478 / 0.2e1 + t643);
t776 = Icges(4,1) * t540;
t865 = t551 * t643 + t553 * t641 + t440 / 0.2e1 + t442 / 0.2e1 - t534 + t766 / 0.2e1 - t776 / 0.2e1;
t507 = Icges(4,2) * t541 + t775;
t864 = t551 * (t446 / 0.2e1 + t444 / 0.2e1) + t553 * (-t450 / 0.2e1 - t448 / 0.2e1) - t439 / 0.2e1 - t441 / 0.2e1 + t507 / 0.2e1 - t510 / 0.2e1;
t476 = -Icges(4,2) * t744 - t525;
t477 = t507 * t554;
t619 = -t534 - t776;
t480 = t619 * t552;
t481 = t619 * t554;
t863 = (t554 * (t457 - t480) + (-t458 + t481) * t552) * t541 + (t554 * (t459 + t476) + (-t460 + t477) * t552) * t540;
t861 = 0.4e1 * qJD(1);
t860 = 2 * qJD(3);
t858 = 2 * qJD(4);
t857 = 4 * qJD(4);
t601 = t372 * t554 - t376 * t552;
t199 = t601 * t541 + (-t428 * t552 - t554 * t910) * t540;
t271 = t601 * t540;
t849 = m(5) * (t199 * t271 + t243 * t934 - t244 * t304);
t576 = -t552 * t715 + t554 * t716;
t103 = t576 * t541 + (-t552 * t703 + t554 * t704) * t540;
t189 = t576 * t540;
t727 = -t226 * t744 + t743 * t951;
t844 = m(6) * (-t541 * t103 + (t135 * t554 + t136 * t552 + t189) * t540 + t727);
t843 = m(6) * (t103 * t189 + t135 * t951 - t136 * t226);
t505 = t679 * t540;
t837 = m(6) * (t189 * t505 + t727);
t828 = -t554 / 0.4e1;
t824 = m(3) * t679 * (rSges(3,3) + qJ(2));
t785 = rSges(4,1) * t541;
t633 = t535 + t785;
t680 = rSges(4,2) * t747 + t554 * rSges(4,3);
t401 = -t552 * t633 - t637 + t680;
t632 = -rSges(4,2) * t746 + t552 * rSges(4,3);
t402 = t554 * t633 - t536 + t632;
t823 = m(4) * (t401 * t492 - t402 * t493);
t822 = m(4) * (t401 * t554 + t402 * t552);
t816 = m(5) * (t317 * t348 + t347 * t897);
t815 = m(5) * (t317 * t406 - t404 * t897);
t813 = m(5) * t878;
t721 = t290 * t744 + t743 * t898;
t807 = m(6) * (-t311 * t746 + t313 * t747 + t721);
t806 = m(6) * (-t226 * t746 - t747 * t951);
t804 = m(6) * ((t335 * t554 + t336 * t552) * t540 + t721);
t803 = m(6) * (t290 * t336 + t335 * t898);
t802 = m(6) * (t290 * t338 + t337 * t898);
t801 = m(6) * (t541 * t924 + (t409 * t552 + t410 * t554) * t540);
t800 = m(6) * (t290 * t552 + t554 * t898);
t797 = m(6) * t247;
t787 = m(6) * qJD(3);
t786 = m(6) * qJD(5);
t750 = t457 * t540;
t557 = (t243 * t552 - t244 * t554) * t856 + (t135 * t552 - t136 * t554) * t853;
t577 = (-t409 * t554 + t410 * t552) * t852;
t52 = t577 + t557;
t749 = t52 * qJD(2);
t748 = t540 * t541;
t719 = -t311 * t744 - t313 * t743;
t697 = t552 * t456 + t460 * t743;
t691 = -t454 - t515;
t688 = t552 * (pkin(7) * t744 - t530) + t554 * (-pkin(3) * t746 + t531);
t683 = t679 * t748;
t678 = qJD(1) * t540;
t677 = qJD(1) * t541;
t676 = qJD(3) * t552;
t675 = qJD(3) * t554;
t674 = qJD(4) * t540;
t673 = qJD(4) * t541;
t384 = m(6) * t505;
t672 = t384 * qJD(1);
t566 = t540 * t588 + t757;
t142 = t415 * t501 - t419 * t502 + t552 * t566;
t565 = t540 * t587 + t756;
t143 = t416 * t501 - t420 * t502 + t552 * t565;
t564 = t540 * t586 + t755;
t144 = t417 * t501 - t421 * t502 + t552 * t564;
t563 = t540 * t585 + t754;
t145 = t418 * t501 - t422 * t502 + t552 * t563;
t669 = (-t905 + t959) * t893 + ((t143 + t145) * t554 + (t142 + t144) * t552 + t902) * t834;
t146 = -t503 * t415 - t504 * t419 + t554 * t566;
t147 = -t503 * t416 - t504 * t420 + t554 * t565;
t148 = -t503 * t417 - t504 * t421 + t554 * t564;
t149 = -t503 * t418 - t504 * t422 + t554 * t563;
t668 = (t608 + t609 - t904) * t893 + ((t147 + t149) * t554 + (t146 + t148) * t552 + t901) * t834;
t667 = (t909 * t552 + t908 * t554 + t900) * t894 + (t896 - t903) * t833;
t666 = t931 * t827 + t930 * t892;
t665 = t929 * t829 + t928 * t831;
t662 = t893 * t900 + t894 * t896;
t655 = -t515 - t699;
t653 = t747 / 0.4e1;
t646 = t175 / 0.2e1 + t177 / 0.2e1;
t645 = -t176 / 0.2e1 - t178 / 0.2e1;
t639 = t471 / 0.2e1 + t470 / 0.2e1;
t628 = t458 * t540 - t455;
t614 = -Icges(4,5) * t540 - Icges(4,6) * t541;
t597 = -t666 - t669;
t596 = -t665 + t668;
t593 = t627 * t540;
t295 = -t458 * t746 + t697;
t575 = (t554 * t628 + t295 - t697) * t827 + (-t552 * (-t459 * t541 + t750) - t554 * t455) * t829 + (t552 * t628 + t629 + t899) * t831;
t574 = t295 * t831 + t697 * t892 + (-t431 + (t456 + t750) * t554 + t698 + t899) * t829;
t567 = -t870 + (t906 - t926) * t552 / 0.4e1 + (t907 - t927) * t828;
t556 = (-t747 / 0.4e1 + t653) * (t936 + t935) + (t828 + t554 / 0.4e1) * (t945 + t946);
t555 = -t869 + t900 * t834 + t903 * t833 + (t905 + t909) * t653 + (t904 + t908) * t746 / 0.4e1 + (t902 + t938) * t744 / 0.4e1 + (t901 + t937) * t743 / 0.4e1;
t512 = -rSges(4,2) * t540 + t785;
t473 = t614 * t554;
t472 = t614 * t552;
t407 = t683 - t748;
t380 = t691 * t554;
t378 = t691 * t552;
t320 = -t541 * t406 - t491 * t746;
t319 = t404 * t541 + t491 * t747;
t314 = t655 * t554;
t312 = t655 * t552;
t291 = (t404 * t554 - t406 * t552) * t540;
t274 = -t338 * t541 + t554 * t593;
t273 = -t337 * t541 - t552 * t593;
t267 = t428 * t554 - t552 * t910 + t688;
t250 = t797 / 0.2e1;
t246 = -t541 * t471 + (-t551 * t693 + t553 * t695) * t540;
t245 = -t541 * t470 + (-t551 * t694 + t553 * t696) * t540;
t192 = t552 * t704 + t554 * t703 + t688;
t190 = t801 / 0.2e1;
t137 = t804 / 0.2e1;
t109 = t806 / 0.2e1;
t108 = t807 / 0.2e1;
t100 = t173 * t505 + t719;
t88 = t800 - t813 + t822 + t824;
t76 = t837 / 0.2e1;
t75 = -t148 * t554 + t149 * t552;
t74 = -t146 * t554 + t147 * t552;
t73 = -t144 * t554 + t145 * t552;
t72 = -t142 * t554 + t143 * t552;
t60 = t540 * t866 - t639 * t541 + t802 + t815;
t57 = t657 + t658;
t51 = t577 - t557;
t42 = t108 - t804 / 0.2e1;
t41 = t108 + t137;
t40 = t137 - t807 / 0.2e1;
t39 = -t540 * t864 - t541 * t865 + t803 + t816 + t823;
t21 = t109 + t250;
t20 = t109 - t797 / 0.2e1;
t19 = t250 - t806 / 0.2e1;
t16 = t844 / 0.2e1;
t13 = t722 + t736;
t10 = t76 + t190 - t844 / 0.2e1;
t9 = t76 + t16 - t801 / 0.2e1;
t8 = t190 + t16 - t837 / 0.2e1;
t7 = t552 * t665 + t554 * t666 - t957;
t6 = t552 * t575 + t554 * t574;
t4 = t849 + t843 + (t552 * t664 + t554 * t663 + t667) * t541 + (t552 * t669 + t554 * t668 - t662) * t540;
t3 = t555 + t567;
t2 = t556 + t555 + (-t210 / 0.4e1 - t209 / 0.4e1 - t178 / 0.4e1 - t176 / 0.4e1) * t552 + (t208 / 0.4e1 + t207 / 0.4e1 + t177 / 0.4e1 + t175 / 0.4e1) * t554 + t870;
t1 = (t197 / 0.2e1 + t196 / 0.2e1 + (-t265 / 0.4e1 - t264 / 0.4e1 - t242 / 0.4e1 - t239 / 0.4e1) * t554 + (-t261 / 0.4e1 - t260 / 0.4e1 + t241 / 0.4e1 + t238 / 0.4e1) * t552) * t541 + t556 + (-t283 / 0.2e1 - t282 / 0.2e1 + (-t187 / 0.4e1 - t186 / 0.4e1 - t167 / 0.4e1 - t165 / 0.4e1) * t554 + (-t185 / 0.4e1 - t184 / 0.4e1 - t166 / 0.4e1 - t164 / 0.4e1) * t552) * t540 + t567 + t869;
t5 = [t88 * qJD(2) + t39 * qJD(3) + t60 * qJD(4) + t195 * t786, qJD(1) * t88 + qJD(3) * t57 + qJD(4) * t13, t39 * qJD(1) + t57 * qJD(2) + t3 * qJD(4) + t41 * qJD(5) + ((t317 * t378 - t347 * t379 - t348 * t377 + t380 * t897) * t855 + (t290 * t312 - t311 * t336 - t313 * t335 + t314 * t898) * t852) * t860 + (m(4) * (-t401 * t512 - t492 * t511) - t184 / 0.2e1 - t185 / 0.2e1 + t506 * t827 + (-t459 / 0.2e1 - t476 / 0.2e1) * t541 + (t457 / 0.2e1 - t480 / 0.2e1) * t540 - t574 - t872) * t675 + (m(4) * (-t402 * t512 + t493 * t511) + t186 / 0.2e1 + t187 / 0.2e1 + t506 * t831 + (t460 / 0.2e1 - t477 / 0.2e1) * t541 + (-t458 / 0.2e1 + t481 / 0.2e1) * t540 - t575 - t871) * t676, t60 * qJD(1) + t13 * qJD(2) + t3 * qJD(3) + t21 * qJD(5) + ((-t226 * t338 + t273 * t898 + t274 * t290 + t337 * t951) * t852 + (-t304 * t406 + t317 * t320 + t319 * t897 - t404 * t934) * t855) * t858 + (-t245 - t246) * t673 + ((t210 / 0.2e1 + t209 / 0.2e1 - t645) * t554 + (t208 / 0.2e1 + t207 / 0.2e1 + t646) * t552) * t674, t41 * qJD(3) + t21 * qJD(4) + t932; -t59 * qJD(3) - t15 * qJD(4) - t384 * qJD(5) + (-t800 / 0.4e1 + t813 / 0.4e1 - t822 / 0.4e1 - t824 / 0.4e1) * t861, 0, -t942 + ((-t378 * t554 + t380 * t552) * t855 + (-t312 * t554 + t314 * t552) * t852) * t860 + t51 * qJD(4), -t958 + t51 * qJD(3) + ((t319 * t552 - t320 * t554) * t855 + (t273 * t552 - t274 * t554) * t852) * t858, -t672; t59 * qJD(2) + t6 * qJD(3) + t1 * qJD(4) + t42 * qJD(5) + (-t823 / 0.4e1 - t816 / 0.4e1 - t803 / 0.4e1) * t861 + t865 * t677 + t864 * t678, qJD(4) * t52 + t942, t6 * qJD(1) + t7 * qJD(4) + t100 * t786 + (m(5) * (t248 * t267 - t377 * t378 - t379 * t380) + m(4) * (t512 * t883 - (t552 * (rSges(4,1) * t744 - t680) + t554 * (rSges(4,1) * t743 + t632)) * t381) + m(6) * (t173 * t192 - t311 * t312 - t313 * t314) + (t75 + t74 + t547 * t473 + (-t552 * t472 + t863) * t554) * t831 + (t73 + t72 + t548 * t472 + (-t554 * t473 + t863) * t552) * t829) * qJD(3), t1 * qJD(1) + t749 + t7 * qJD(3) + (t889 * t829 + t888 * t831) * qJD(4) + t10 * qJD(5) + (-t843 / 0.4e1 - t849 / 0.4e1) * t857 + ((t291 * t248 + t271 * t306 - t319 * t379 - t320 * t377 + t879 * t491) * t855 + (-t173 * t247 - t189 * t924 - t226 * t409 - t273 * t313 - t274 * t311 + t410 * t951) * t852) * t858 + ((t646 - t663) * t554 + (t645 - t664) * t552 - t667) * t673 + (t552 * t597 - t554 * t596 + t662) * t674, t42 * qJD(1) + t100 * t787 + t10 * qJD(4) + (-t505 * t541 - t407 + t683) * t786; t15 * qJD(2) + t2 * qJD(3) + t20 * qJD(5) + t639 * t677 + (-t802 / 0.4e1 - t815 / 0.4e1) * t861 - t866 * t678, -qJD(3) * t52 + t958, t2 * qJD(1) - t749 + t4 * qJD(4) + t9 * qJD(5) + ((t103 * t173 - t135 * t313 - t136 * t311 + t189 * t192 - t226 * t312 + t314 * t951) * t852 + (t199 * t248 - t243 * t379 - t244 * t377 + t267 * t271 - t304 * t378 + t380 * t934) * t855) * t860 + t597 * t675 + t596 * t676 + (((t936 / 0.2e1 + t935 / 0.2e1 + t872) * t554 + (t961 * t829 + t960 * t831 + t871) * t552) * t541 + ((t75 / 0.2e1 + t74 / 0.2e1 + t238 / 0.2e1 + t241 / 0.2e1) * t554 + (t72 / 0.2e1 + t73 / 0.2e1 + t239 / 0.2e1 + t242 / 0.2e1) * t552) * t540 + t957) * qJD(3), t4 * qJD(3) + ((-t189 * t247 - t226 * t274 + t273 * t951) * m(6) / 0.4e1 + (t271 * t291 - t304 * t320 + t319 * t934) * m(5) / 0.4e1) * t857 + (t246 / 0.2e1 + t245 / 0.2e1) * qJD(4) * t541 ^ 2 + ((t907 * t552 + t906 * t554) * t833 + t889 * t831 + t888 * t827) * t674, qJD(1) * t20 + qJD(3) * t9; t384 * qJD(2) + t40 * qJD(3) + t19 * qJD(4) - t932, t672, t40 * qJD(1) + (-t541 * t192 + (t312 * t552 + t314 * t554 + t173) * t540 - t100 + t719) * t787 + t8 * qJD(4) + t407 * t786, t19 * qJD(1) + t8 * qJD(3) + m(6) * (t541 * t247 + (t273 * t554 + t274 * t552) * t540) * qJD(4), t407 * t787;];
Cq = t5;
