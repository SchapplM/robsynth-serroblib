% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:41
% EndTime: 2022-01-20 12:08:12
% DurationCPUTime: 16.27s
% Computational Cost: add. (109611->646), mult. (69111->802), div. (0->0), fcn. (63994->10), ass. (0->408)
t588 = cos(qJ(3));
t582 = t588 * pkin(3);
t573 = t582 + pkin(2);
t584 = qJ(3) + qJ(4);
t578 = cos(t584);
t791 = pkin(4) * t578;
t527 = t573 + t791;
t585 = qJ(1) + qJ(2);
t579 = cos(t585);
t861 = -pkin(8) - pkin(7);
t583 = pkin(9) - t861;
t551 = t583 * t579;
t577 = sin(t585);
t581 = qJ(5) + t584;
t571 = sin(t581);
t759 = t571 * t577;
t679 = rSges(6,2) * t759 + t579 * rSges(6,3);
t572 = cos(t581);
t787 = rSges(6,1) * t572;
t882 = t551 + (-t527 - t787) * t577 + t679;
t795 = sin(qJ(1)) * pkin(1);
t890 = t882 - t795;
t909 = t890 - t882;
t754 = t572 * t579;
t758 = t571 * t579;
t668 = rSges(6,1) * t754 - rSges(6,2) * t758 + t577 * rSges(6,3);
t905 = t579 * t527 + t583 * t577;
t346 = t668 + t905;
t794 = cos(qJ(1)) * pkin(1);
t331 = t794 + t346;
t910 = -t331 + t346;
t755 = t572 * t577;
t457 = rSges(6,1) * t759 + rSges(6,2) * t755;
t500 = rSges(6,1) * t571 + rSges(6,2) * t572;
t458 = t500 * t579;
t204 = -t346 * t458 + t457 * t882;
t908 = m(6) * t204;
t202 = -t331 * t458 + t457 * t890;
t907 = m(6) * t202;
t576 = sin(t584);
t792 = pkin(4) * t576;
t665 = t500 + t792;
t642 = t665 * t579;
t886 = t331 * t642;
t885 = t346 * t642;
t554 = Icges(6,4) * t572;
t497 = -Icges(6,2) * t571 + t554;
t498 = Icges(6,1) * t571 + t554;
t904 = t497 + t498;
t564 = Icges(5,4) * t578;
t512 = -Icges(5,2) * t576 + t564;
t513 = Icges(5,1) * t576 + t564;
t903 = t512 + t513;
t580 = Icges(4,4) * t588;
t586 = sin(qJ(3));
t537 = -Icges(4,2) * t586 + t580;
t538 = Icges(4,1) * t586 + t580;
t902 = t537 + t538;
t175 = t331 * t882 - t346 * t890;
t864 = m(4) / 0.2e1;
t863 = m(5) / 0.2e1;
t862 = m(6) / 0.2e1;
t845 = t577 / 0.2e1;
t844 = -t579 / 0.2e1;
t899 = t579 / 0.2e1;
t843 = m(3) * (t794 * (-rSges(3,1) * t577 - rSges(3,2) * t579) + (t579 * rSges(3,1) - t577 * rSges(3,2)) * t795);
t574 = t577 ^ 2;
t575 = t579 ^ 2;
t676 = t574 + t575;
t515 = rSges(5,1) * t576 + rSges(5,2) * t578;
t793 = pkin(3) * t586;
t614 = t515 + t793;
t883 = t614 * t579;
t884 = t614 * t577;
t617 = (t577 * t883 - t579 * t884) * t515;
t744 = t577 * t578;
t751 = t576 * t577;
t435 = rSges(5,1) * t744 - rSges(5,2) * t751 - t579 * rSges(5,3);
t678 = t577 * t573 + t579 * t861;
t362 = -t435 - t678;
t552 = t577 * t861;
t750 = t576 * t579;
t663 = -rSges(5,2) * t750 + t577 * rSges(5,3);
t788 = rSges(5,1) * t578;
t363 = -t552 + (t573 + t788) * t579 + t663;
t517 = -rSges(5,2) * t576 + t788;
t621 = (-t362 * t579 - t363 * t577) * t517;
t528 = -t792 - t793;
t375 = -t528 * t577 + t457;
t494 = t579 * t528;
t376 = t494 - t458;
t405 = t665 * t577;
t706 = -t375 * t642 - t405 * t376;
t501 = -rSges(6,2) * t571 + t787;
t664 = -t501 - t791;
t406 = t664 * t577;
t408 = t664 * t579;
t713 = t406 * t346 + t408 * t882;
t718 = (t706 + t713) * t862 + (t621 + t617) * t863;
t473 = t515 * t577;
t474 = t515 * t579;
t700 = -t473 * t883 + t474 * t884;
t628 = t500 - t528;
t371 = t628 * t577;
t373 = t628 * t579;
t394 = pkin(4) * t751 + t457;
t709 = t371 * t642 - t373 * t394;
t719 = (t709 + t713) * t862 + (t621 + t700) * t863;
t893 = t718 - t719;
t352 = t362 - t795;
t353 = t363 + t794;
t622 = (-t352 * t579 - t353 * t577) * t517;
t714 = t406 * t331 + t408 * t890;
t720 = (t706 + t714) * t862 + (t622 + t617) * t863;
t721 = (t709 + t714) * t862 + (t622 + t700) * t863;
t892 = t720 - t721;
t846 = -t577 / 0.2e1;
t891 = t845 + t846;
t203 = t352 * t884 - t353 * t883;
t205 = t362 * t884 - t363 * t883;
t193 = -t363 * t352 + t353 * t362;
t570 = t579 * pkin(7);
t789 = rSges(4,1) * t588;
t667 = pkin(2) + t789;
t743 = t577 * t586;
t677 = rSges(4,2) * t743 + t579 * rSges(4,3);
t380 = -t577 * t667 + t570 + t677;
t369 = t380 - t795;
t730 = t579 * t586;
t548 = rSges(4,2) * t730;
t381 = -t548 + t667 * t579 + (rSges(4,3) + pkin(7)) * t577;
t370 = t381 + t794;
t206 = -t381 * t369 + t370 * t380;
t739 = t578 * t579;
t317 = t577 * t435 + t579 * (rSges(5,1) * t739 + t663);
t694 = -t577 * (pkin(2) * t577 - t570 - t678) + t579 * (-t577 * pkin(7) - t552 + (-pkin(2) + t573) * t579);
t212 = t317 + t694;
t350 = -t577 * t473 - t579 * t474;
t138 = t212 * t350 + (t577 * t884 + t579 * t883) * t517;
t311 = t577 * (rSges(6,1) * t755 - t679) + t579 * t668;
t201 = t577 * (t527 * t577 - t551 - t678) + t579 * (-t579 * t573 + t552 + t905) + t311;
t164 = t201 + t694;
t334 = -t577 * t457 - t579 * t458;
t287 = -t676 * t792 + t334;
t707 = -t371 * t406 - t373 * t408;
t88 = t164 * t287 + t707;
t881 = m(5) * t138 + m(6) * t88;
t675 = qJD(1) + qJD(2);
t723 = (t909 * t405 + t885 - t886) * t862 + ((-t353 + t363) * t579 + (t352 - t362) * t577) * t515 * t863;
t196 = t394 * t890 - t886;
t200 = t394 * t882 - t885;
t209 = t352 * t473 - t353 * t474;
t214 = t362 * t473 - t363 * t474;
t879 = (t200 + t196) * t862 + (t214 + t209) * t863;
t540 = rSges(4,1) * t586 + rSges(4,2) * t588;
t671 = (t203 - t205) * t863 + ((-t370 + t381) * t579 + (t369 - t380) * t577) * t540 * t864 + (t909 * t371 + t910 * t373) * t862;
t192 = t331 * t376 + t375 * t890;
t195 = t346 * t376 + t375 * t882;
t489 = t540 * t577;
t490 = t540 * t579;
t229 = t369 * t489 - t370 * t490;
t236 = t380 * t489 - t381 * t490;
t877 = (t205 + t203) * t863 + (t236 + t229) * t864 + (t195 + t192) * t862;
t779 = Icges(4,4) * t586;
t536 = Icges(4,2) * t588 + t779;
t539 = Icges(4,1) * t588 - t779;
t876 = t902 * t588 / 0.2e1 + (t539 / 0.2e1 - t536 / 0.2e1) * t586;
t432 = Icges(5,6) * t577 + t512 * t579;
t778 = Icges(5,4) * t576;
t514 = Icges(5,1) * t578 - t778;
t434 = Icges(5,5) * t577 + t514 * t579;
t377 = t434 * t744;
t510 = Icges(5,5) * t578 - Icges(5,6) * t576;
t765 = t510 * t579;
t430 = Icges(5,3) * t577 + t765;
t657 = t579 * t430 - t377;
t238 = -t432 * t751 - t657;
t431 = Icges(5,4) * t744 - Icges(5,2) * t751 - Icges(5,6) * t579;
t429 = Icges(5,5) * t744 - Icges(5,6) * t751 - Icges(5,3) * t579;
t531 = Icges(5,4) * t751;
t433 = Icges(5,1) * t744 - Icges(5,5) * t579 - t531;
t697 = -t577 * t429 - t433 * t739;
t239 = -t431 * t750 - t697;
t696 = t577 * t430 + t434 * t739;
t240 = -t432 * t750 + t696;
t654 = t432 * t576 - t429;
t769 = t431 * t576;
t875 = (-t239 * t579 + t240 * t577) * t899 + ((t238 - t377 + (t430 + t769) * t579 + t697) * t579 + t696 * t577) * t844 + ((t577 * t654 + t238 + t239 + t657) * t577 + (t577 * (-t433 * t578 + t769) + t240 - t696 + (t429 + t654) * t579) * t579) * t845;
t511 = Icges(5,2) * t578 + t778;
t874 = t903 * t578 / 0.2e1 + (t514 / 0.2e1 - t511 / 0.2e1) * t576;
t415 = Icges(6,6) * t577 + t497 * t579;
t777 = Icges(6,4) * t571;
t499 = Icges(6,1) * t572 - t777;
t417 = Icges(6,5) * t577 + t499 * t579;
t366 = t417 * t755;
t495 = Icges(6,5) * t572 - Icges(6,6) * t571;
t767 = t495 * t579;
t413 = Icges(6,3) * t577 + t767;
t658 = t579 * t413 - t366;
t225 = -t415 * t759 - t658;
t414 = Icges(6,4) * t755 - Icges(6,2) * t759 - Icges(6,6) * t579;
t412 = Icges(6,5) * t755 - Icges(6,6) * t759 - Icges(6,3) * t579;
t521 = Icges(6,4) * t759;
t416 = Icges(6,1) * t755 - Icges(6,5) * t579 - t521;
t699 = -t577 * t412 - t416 * t754;
t226 = -t414 * t758 - t699;
t698 = t577 * t413 + t417 * t754;
t227 = -t415 * t758 + t698;
t655 = t415 * t571 - t412;
t770 = t414 * t571;
t652 = (-t226 * t579 + t227 * t577) * t899 + ((t225 - t366 + (t413 + t770) * t579 + t699) * t579 + t698 * t577) * t844 + ((t577 * t655 + t225 + t226 + t658) * t577 + (t577 * (-t416 * t572 + t770) + t227 - t698 + (t412 + t655) * t579) * t579) * t845;
t496 = Icges(6,2) * t572 + t777;
t650 = t904 * t572 / 0.2e1 + (-t496 / 0.2e1 + t499 / 0.2e1) * t571;
t684 = -t511 * t579 + t434;
t685 = -Icges(5,2) * t744 + t433 - t531;
t686 = -t513 * t579 - t432;
t687 = t513 * t577 + t431;
t873 = (-t684 * t577 + t579 * t685) * t576 + (t686 * t577 + t579 * t687) * t578;
t688 = -t496 * t579 + t417;
t689 = -Icges(6,2) * t755 + t416 - t521;
t690 = -t498 * t579 - t415;
t691 = t498 * t577 + t414;
t872 = (-t688 * t577 + t579 * t689) * t571 + (t690 * t577 + t579 * t691) * t572;
t870 = 4 * qJD(1);
t868 = 4 * qJD(2);
t867 = 2 * qJD(3);
t866 = 2 * qJD(4);
t131 = t164 * t334;
t179 = t201 * t334;
t856 = m(6) * (t131 + t179 + ((t373 + t642) * t579 + (t371 + t405) * t577) * t501);
t207 = t311 * t287;
t618 = (-t406 * t577 - t408 * t579) * t500;
t89 = t131 + (t371 * t577 + t373 * t579) * t501;
t853 = m(6) * (t207 + t618 + t89);
t228 = t574 * (t528 + t793) + t579 * (pkin(3) * t730 + t494) + t334;
t113 = t179 + (t405 * t577 + t579 * t642) * t501;
t595 = t618 + t113;
t852 = m(6) * (t228 * t311 + t595);
t703 = -t405 * t406 - t408 * t642;
t849 = m(6) * (t201 * t228 + t703);
t839 = m(4) * t206;
t837 = m(4) * t229;
t836 = m(4) * t236;
t826 = m(5) * t193;
t199 = t676 * t515 * t517 + t317 * t350;
t197 = m(5) * t199;
t824 = m(5) * t203;
t823 = m(5) * t205;
t822 = m(5) * t209;
t821 = m(5) * t214;
t816 = m(6) * (t204 + t202);
t812 = m(6) * (t909 * t577 + t910 * t579) * t500;
t624 = (-t331 * t577 - t579 * t890) * t501;
t704 = t371 * t458 - t373 * t457;
t810 = m(6) * (t624 + t704);
t623 = (-t346 * t577 - t579 * t882) * t501;
t809 = m(6) * (t623 + t704);
t620 = (-t375 * t579 - t376 * t577) * t500;
t808 = m(6) * (t624 + t620);
t701 = t405 * t458 - t457 * t642;
t807 = m(6) * (t624 + t701);
t806 = m(6) * (t623 + t620);
t619 = (-t394 * t579 + t577 * t642) * t500;
t805 = m(6) * (t624 + t619);
t804 = m(6) * (t623 + t701);
t803 = m(6) * (t623 + t619);
t802 = m(6) * t175;
t189 = t676 * t500 * t501 + t311 * t334;
t800 = m(6) * t189;
t799 = m(6) * t192;
t798 = m(6) * t195;
t797 = m(6) * t196;
t796 = m(6) * t200;
t635 = Icges(6,5) * t571 + Icges(6,6) * t572;
t450 = t635 * t577;
t451 = t579 * t635;
t790 = (-t574 * t451 + (t577 * t450 + t872) * t579) * t845 + (-t575 * t450 + (t579 * t451 + t872) * t577) * t844;
t742 = t577 * t588;
t446 = Icges(4,4) * t742 - Icges(4,2) * t743 - Icges(4,6) * t579;
t768 = t446 * t586;
t535 = Icges(4,5) * t588 - Icges(4,6) * t586;
t762 = t535 * t579;
t729 = t579 * t588;
t711 = -t371 * t376 - t373 * t375;
t705 = (-t394 + t405) * t642;
t444 = Icges(4,5) * t742 - Icges(4,6) * t743 - Icges(4,3) * t579;
t546 = Icges(4,4) * t743;
t448 = Icges(4,1) * t742 - Icges(4,5) * t579 - t546;
t693 = -t577 * t444 - t448 * t729;
t445 = Icges(4,3) * t577 + t762;
t449 = Icges(4,5) * t577 + t539 * t579;
t692 = t577 * t445 + t449 * t729;
t683 = t538 * t577 + t446;
t447 = Icges(4,6) * t577 + t537 * t579;
t682 = -t538 * t579 - t447;
t681 = -Icges(4,2) * t742 + t448 - t546;
t680 = -t536 * t579 + t449;
t674 = t852 / 0.2e1 + t790;
t672 = -t800 + t790;
t669 = t201 * t287 + t703;
t666 = -t517 - t582;
t402 = t449 * t742;
t656 = t579 * t445 - t402;
t653 = t447 * t586 - t444;
t636 = Icges(5,5) * t576 + Icges(5,6) * t578;
t467 = t636 * t577;
t468 = t579 * t636;
t651 = (-t574 * t468 + (t577 * t467 + t873) * t579) * t845 + (-t575 * t467 + (t579 * t468 + t873) * t577) * t844 + t790;
t648 = t676 * t793;
t643 = t197 + t651;
t641 = t816 / 0.2e1 + t650;
t637 = Icges(4,5) * t586 + Icges(4,6) * t588;
t627 = t664 - t582;
t616 = (-t473 * t579 + t474 * t577) * t515;
t615 = (-t489 * t579 + t490 * t577) * t540;
t613 = t652 + t875;
t606 = (-t496 + t499) * t572 - t904 * t571;
t611 = -t652 + (t495 * t577 + t571 * t690 + t572 * t688 + t579 * t606) * t845 + (-t571 * t691 + t572 * t689 + t577 * t606 - t767) * t844;
t610 = -t650 + t891 * (t414 * t572 + t416 * t571);
t609 = t650 + t874;
t608 = t586 * t681 + t588 * t683;
t607 = -t586 * t680 + t588 * t682;
t605 = (-t511 + t514) * t578 - t903 * t576;
t604 = (-t536 + t539) * t588 - t902 * t586;
t600 = (-t457 * t579 + t458 * t577) * t500;
t599 = t609 + t879;
t598 = t609 + t876;
t597 = t598 + t877;
t596 = t611 - t875 + (t510 * t577 + t576 * t686 + t578 * t684 + t579 * t605) * t845 + (-t576 * t687 + t577 * t605 + t578 * t685 - t765) * t844;
t594 = t610 - t874 + t891 * (t431 * t578 + t433 * t576);
t592 = t596 * qJD(4);
t591 = t594 - t876 + t891 * (t446 * t588 + t586 * t448);
t256 = -t447 * t743 - t656;
t190 = -(-t577 * (-t448 * t588 + t768) - t579 * t444) * t579 + t256 * t577;
t257 = -t446 * t730 - t693;
t258 = -t447 * t730 + t692;
t191 = -t257 * t579 + t258 * t577;
t79 = (t579 * t653 + t258 - t692) * t579 + (t577 * t653 + t257 + t656) * t577;
t80 = (t256 - t402 + (t445 + t768) * t579 + t693) * t579 + t692 * t577;
t590 = (t80 * t899 + t596 + (t79 + t190) * t846 + (t535 * t577 + t579 * t604 + t586 * t682 + t588 * t680) * t845 + (t577 * t604 - t586 * t683 + t588 * t681 + t191 - t762) * t844) * qJD(3);
t542 = -rSges(4,2) * t586 + t789;
t484 = t579 * t637;
t483 = t637 * t577;
t439 = t666 * t579;
t437 = t666 * t577;
t374 = t627 * t579;
t372 = t627 * t577;
t310 = -t648 + t350;
t215 = -t648 + t228;
t160 = t650 + t908;
t154 = t650 + t907;
t151 = t803 / 0.2e1;
t149 = t804 / 0.2e1;
t144 = t805 / 0.2e1;
t143 = t806 / 0.2e1;
t142 = t807 / 0.2e1;
t137 = t808 / 0.2e1;
t136 = t809 / 0.2e1;
t133 = t810 / 0.2e1;
t124 = t812 / 0.2e1;
t87 = t609 + t796 + t821;
t83 = t609 + t797 + t822;
t72 = t802 + t826 + t839 + t843;
t68 = t853 / 0.2e1;
t60 = t598 + t798 + t823 + t836;
t59 = t598 + t799 + t824 + t837;
t49 = t856 / 0.2e1;
t45 = t790 + t800;
t44 = t45 * qJD(5);
t43 = -t812 / 0.2e1 + t641;
t42 = t124 + t641;
t41 = m(6) * t113 + t790;
t40 = m(6) * t89 + t790;
t36 = t124 - t816 / 0.2e1 + t610;
t35 = t599 - t723;
t34 = t599 + t723;
t33 = t643 + t849;
t31 = t652 * qJD(5);
t30 = t651 + t881;
t29 = t594 + t723 - t879;
t28 = t49 - t853 / 0.2e1 + t674;
t27 = t49 + t68 - t852 / 0.2e1 + t790;
t26 = t68 - t856 / 0.2e1 + t674;
t25 = t597 - t671;
t24 = t597 + t671;
t23 = t149 - t803 / 0.2e1 + t652;
t22 = t151 - t804 / 0.2e1 + t652;
t21 = t142 - t805 / 0.2e1 + t652;
t20 = t144 - t807 / 0.2e1 + t652;
t19 = t136 - t806 / 0.2e1 + t652;
t18 = t143 - t809 / 0.2e1 + t652;
t17 = t133 - t808 / 0.2e1 + t652;
t16 = t137 - t810 / 0.2e1 + t652;
t15 = t591 + t671 - t877;
t14 = t149 + t151 + t611;
t13 = t142 + t144 + t611;
t12 = t136 + t143 + t611;
t11 = t133 + t137 + t611;
t9 = t613 * qJD(4);
t8 = (-t80 / 0.2e1 + t191 / 0.2e1) * t579 + (t190 / 0.2e1 + t79 / 0.2e1) * t577 + t613;
t7 = t8 * qJD(3);
t6 = t613 + t893;
t5 = t613 - t893;
t4 = t613 + t892;
t3 = t613 - t892;
t2 = t596 + t718 + t719;
t1 = t596 + t720 + t721;
t10 = [t72 * qJD(2) + t59 * qJD(3) + t83 * qJD(4) + t154 * qJD(5), t72 * qJD(1) + t24 * qJD(3) + t34 * qJD(4) + t42 * qJD(5) + 0.2e1 * (t175 * t862 + t193 * t863 + t206 * t864 + t843 / 0.2e1) * qJD(2), t59 * qJD(1) + t24 * qJD(2) + t590 + t1 * qJD(4) + t11 * qJD(5) + ((t352 * t439 + t353 * t437) * t863 + ((-t369 * t579 - t370 * t577) * t542 + t615) * t864 + (t331 * t372 + t374 * t890 + t711) * t862) * t867, t83 * qJD(1) + t34 * qJD(2) + t1 * qJD(3) + t592 + t13 * qJD(5) + ((t705 + t714) * t862 + (t622 + t616) * t863) * t866, t154 * qJD(1) + t42 * qJD(2) + t11 * qJD(3) + t13 * qJD(4) + ((t624 + t600) * m(6) + t611) * qJD(5); t25 * qJD(3) + t35 * qJD(4) + t43 * qJD(5) + (-t802 / 0.4e1 - t826 / 0.4e1 - t839 / 0.4e1 - t843 / 0.4e1) * t870, qJD(3) * t60 + qJD(4) * t87 + qJD(5) * t160, t25 * qJD(1) + t60 * qJD(2) + t590 + t2 * qJD(4) + t12 * qJD(5) + ((t362 * t439 + t363 * t437) * t863 + ((-t380 * t579 - t381 * t577) * t542 + t615) * t864 + (t346 * t372 + t374 * t882 + t711) * t862) * t867, t35 * qJD(1) + t87 * qJD(2) + t2 * qJD(3) + t592 + t14 * qJD(5) + ((t705 + t713) * t862 + (t621 + t616) * t863) * t866, t43 * qJD(1) + t160 * qJD(2) + t12 * qJD(3) + t14 * qJD(4) + ((t623 + t600) * m(6) + t611) * qJD(5); t591 * qJD(1) + t15 * qJD(2) + t7 + t3 * qJD(4) + t17 * qJD(5) + (-t799 / 0.4e1 - t824 / 0.4e1 - t837 / 0.4e1) * t870, t15 * qJD(1) + t591 * qJD(2) + t7 + t5 * qJD(4) + t19 * qJD(5) + (-t798 / 0.4e1 - t823 / 0.4e1 - t836 / 0.4e1) * t868, (m(6) * (t164 * t215 - t371 * t372 - t373 * t374) + m(5) * (t212 * t310 - t437 * t884 - t439 * t883) + m(4) * ((t577 * (rSges(4,1) * t742 - t677) + t579 * (rSges(4,1) * t729 + t577 * rSges(4,3) - t548)) * (-t489 * t577 - t490 * t579) + t676 * t542 * t540) + (-t575 * t483 + (t607 * t577 + (t484 + t608) * t579) * t577) * t844 + (-t574 * t484 + (t608 * t579 + (t483 + t607) * t577) * t579) * t845 + t651) * qJD(3) + t30 * qJD(4) + t40 * qJD(5) + t675 * t8, t3 * qJD(1) + t5 * qJD(2) + t30 * qJD(3) + t27 * qJD(5) + ((t669 + t88) * t862 + (t138 + t199) * t863) * t866 + (t651 - t197 - t849) * qJD(4), t17 * qJD(1) + t19 * qJD(2) + t40 * qJD(3) + t27 * qJD(4) + (m(6) * (t89 + t189) + t672) * qJD(5); t594 * qJD(1) + t29 * qJD(2) + t4 * qJD(3) + t9 + t21 * qJD(5) + (-t797 / 0.4e1 - t822 / 0.4e1) * t870, t29 * qJD(1) + t6 * qJD(3) + t9 + t23 * qJD(5) + (-t796 / 0.4e1 - t821 / 0.4e1) * t868 + t594 * qJD(2), t4 * qJD(1) + t6 * qJD(2) + t33 * qJD(4) + t28 * qJD(5) + ((t164 * t228 + t201 * t215 - t372 * t405 - t374 * t642 + t707) * t862 + (t310 * t317 + (-t437 * t577 - t439 * t579) * t515 + t138) * t863) * t867 + (t651 - t881) * qJD(3), t33 * qJD(3) + (m(6) * t669 + t643) * qJD(4) + t41 * qJD(5) + t675 * t613, t21 * qJD(1) + t23 * qJD(2) + t28 * qJD(3) + t41 * qJD(4) + (m(6) * (t113 + t189) + t672) * qJD(5); (t610 - t907) * qJD(1) + t36 * qJD(2) + t16 * qJD(3) + t20 * qJD(4) + t31, t36 * qJD(1) + (t610 - t908) * qJD(2) + t18 * qJD(3) + t22 * qJD(4) + t31, t16 * qJD(1) + t18 * qJD(2) + ((t215 * t311 + (-t372 * t577 - t374 * t579) * t500) * m(6) + t790) * qJD(3) + t26 * qJD(4) + t44, t20 * qJD(1) + t22 * qJD(2) + t26 * qJD(3) + ((t207 - t113 + t595) * m(6) + t790) * qJD(4) + t44, t44 + (qJD(3) + qJD(4)) * t45 + t675 * t652;];
Cq = t10;
