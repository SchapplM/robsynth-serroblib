% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR12_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:31
% EndTime: 2019-12-31 19:12:59
% DurationCPUTime: 22.41s
% Computational Cost: add. (68537->740), mult. (92364->992), div. (0->0), fcn. (99843->8), ass. (0->463)
t589 = cos(qJ(1));
t585 = sin(qJ(3));
t798 = pkin(3) * t585;
t566 = t589 * t798;
t586 = sin(qJ(1));
t479 = t589 * (-pkin(7) * t586 + t566);
t853 = -pkin(7) - pkin(6);
t568 = t589 * t853;
t741 = t585 * t586;
t679 = pkin(3) * t741;
t518 = -pkin(6) * t589 - t568 + t679;
t583 = qJ(3) + qJ(4);
t569 = sin(t583);
t751 = t569 * t586;
t587 = cos(qJ(5));
t733 = t587 * t589;
t584 = sin(qJ(5));
t736 = t586 * t584;
t519 = -t569 * t736 + t733;
t735 = t586 * t587;
t742 = t584 * t589;
t520 = t569 * t735 + t742;
t570 = cos(t583);
t745 = t570 * t586;
t382 = t520 * rSges(6,1) + t519 * rSges(6,2) - rSges(6,3) * t745;
t877 = -pkin(8) * t745 + t382;
t702 = -pkin(4) * t751 - t877;
t743 = t570 * t589;
t521 = t569 * t742 + t735;
t522 = t569 * t733 - t736;
t876 = t522 * rSges(6,1) - t521 * rSges(6,2);
t384 = rSges(6,3) * t743 - t876;
t562 = pkin(8) * t743;
t750 = t569 * t589;
t707 = (-pkin(4) * t750 + t384 + t562) * t589;
t236 = -t479 + (-t518 + t702) * t586 + t707;
t746 = t570 * t584;
t677 = rSges(6,2) * t746;
t744 = t570 * t587;
t678 = rSges(6,1) * t744;
t424 = rSges(6,3) * t751 + (-t677 + t678) * t586;
t700 = -pkin(4) * t745 - pkin(8) * t751 - t424;
t685 = rSges(6,3) * t750 + t589 * t678;
t425 = t589 * t677 - t685;
t684 = -pkin(4) * t743 - pkin(8) * t750;
t893 = t425 + t684;
t281 = t586 * t700 + t893 * t589;
t641 = rSges(6,1) * t587 - rSges(6,2) * t584;
t459 = rSges(6,3) * t569 + t570 * t641;
t440 = t586 * t459;
t541 = pkin(4) * t570 + pkin(8) * t569;
t389 = t586 * t541 + t440;
t588 = cos(qJ(3));
t734 = t586 * t588;
t565 = pkin(3) * t734;
t363 = t565 + t389;
t697 = t459 + t541;
t797 = pkin(3) * t588;
t365 = (t697 + t797) * t589;
t794 = t570 * rSges(6,3);
t458 = -t569 * t641 + t794;
t439 = t586 * t458;
t796 = t569 * pkin(4);
t540 = pkin(8) * t570 - t796;
t388 = t586 * t540 + t439;
t390 = (-t458 - t540) * t589;
t109 = t236 * t281 + t363 * t388 - t365 * t390;
t643 = rSges(5,1) * t569 + rSges(5,2) * t570;
t757 = t589 * t643;
t603 = -t586 * rSges(5,3) + t757;
t445 = t589 * t603;
t758 = t586 * t643;
t470 = rSges(5,3) * t589 + t758;
t307 = -t445 - t479 + (-t470 - t518) * t586;
t509 = rSges(5,1) * t745 - rSges(5,2) * t751;
t539 = rSges(5,1) * t570 - rSges(5,2) * t569;
t511 = t539 * t589;
t393 = -t509 * t586 - t589 * t511;
t461 = t539 * t586 + t565;
t645 = (t539 + t797) * t589;
t207 = t307 * t393 - t461 * t758 - t645 * t757;
t918 = -m(5) * t207 - m(6) * t109;
t653 = qJ(2) + t798;
t799 = pkin(1) * t589;
t322 = t799 - t568 + (t653 + t796) * t586 + t877;
t405 = rSges(6,1) * t519 - rSges(6,2) * t520;
t406 = rSges(6,1) * t521 + rSges(6,2) * t522;
t778 = Icges(6,4) * t587;
t635 = -Icges(6,2) * t584 + t778;
t455 = Icges(6,6) * t569 + t570 * t635;
t779 = Icges(6,4) * t584;
t638 = Icges(6,1) * t587 - t779;
t457 = Icges(6,5) * t569 + t570 * t638;
t501 = (-Icges(6,2) * t587 - t779) * t570;
t504 = (-Icges(6,1) * t584 - t778) * t570;
t498 = (-Icges(6,5) * t584 - Icges(6,6) * t587) * t570;
t754 = t569 * t498;
t571 = t589 * qJ(2);
t606 = t566 + t571 + (-pkin(1) + t853) * t586;
t892 = (-t794 + t796) * t589 - t562 + t606 + t876;
t122 = t570 * (-(t457 / 0.2e1 + t501 / 0.2e1) * t584 + (t504 / 0.2e1 - t455 / 0.2e1) * t587) + m(6) * (t322 * t405 - t406 * t892) + t754 / 0.2e1;
t917 = t122 * qJD(1);
t372 = Icges(6,5) * t520 + Icges(6,6) * t519 - Icges(6,3) * t745;
t780 = Icges(6,4) * t520;
t375 = Icges(6,2) * t519 - Icges(6,6) * t745 + t780;
t495 = Icges(6,4) * t519;
t378 = Icges(6,1) * t520 - Icges(6,5) * t745 + t495;
t229 = t372 * t743 + t521 * t375 - t522 * t378;
t916 = t229 * t589;
t915 = t586 * t229;
t374 = -Icges(6,5) * t522 + Icges(6,6) * t521 + Icges(6,3) * t743;
t765 = t374 * t569;
t497 = Icges(6,4) * t522;
t377 = Icges(6,2) * t521 + Icges(6,6) * t743 - t497;
t496 = Icges(6,4) * t521;
t379 = Icges(6,1) * t522 - Icges(6,5) * t743 - t496;
t908 = t377 * t584 + t379 * t587;
t247 = t570 * t908 - t765;
t838 = t586 / 0.2e1;
t837 = -t589 / 0.2e1;
t450 = t509 + t565;
t880 = t586 * t645;
t353 = -t450 * t589 + t880;
t354 = -t700 + t565;
t355 = (-t677 + t797) * t589 - t684 + t685;
t551 = rSges(4,1) * t588 - rSges(4,2) * t585;
t531 = t551 * t586;
t532 = t551 * t589;
t407 = -t531 * t589 + t586 * t532;
t855 = m(6) / 0.2e1;
t856 = m(5) / 0.2e1;
t857 = m(4) / 0.2e1;
t667 = (-t354 * t589 + t586 * t355) * t855 + t353 * t856 + t407 * t857;
t713 = (t363 * t589 - t365 * t586) * t855 + (t461 * t589 - t880) * t856;
t92 = t713 - t667;
t914 = qJD(1) * t92;
t625 = -t521 * t377 - t522 * t379;
t708 = t519 * t375 + t520 * t378;
t913 = t625 + t708 + (-t372 * t586 - t374 * t589) * t570;
t782 = Icges(5,4) * t569;
t636 = Icges(5,2) * t570 + t782;
t466 = Icges(5,6) * t589 + t586 * t636;
t553 = Icges(5,4) * t745;
t468 = Icges(5,1) * t751 + Icges(5,5) * t589 + t553;
t631 = Icges(5,5) * t569 + Icges(5,6) * t570;
t899 = t586 * t631;
t311 = t589 * (Icges(5,3) * t589 + t899) + t466 * t745 + t468 * t751;
t898 = t589 * t631;
t465 = -Icges(5,3) * t586 + t898;
t467 = -Icges(5,6) * t586 + t589 * t636;
t781 = Icges(5,4) * t570;
t639 = Icges(5,1) * t569 + t781;
t469 = -Icges(5,5) * t586 + t589 * t639;
t312 = -t589 * t465 - t467 * t745 - t469 * t751;
t878 = (t467 * t570 + t469 * t569) * t589;
t314 = -t586 * t465 + t878;
t230 = t374 * t743 - t625;
t722 = t230 + t913;
t228 = -t374 * t745 + t519 * t377 - t379 * t520;
t885 = t228 * t586;
t51 = t722 * t589 + t885;
t227 = -t372 * t745 + t708;
t719 = -t227 + t913;
t52 = t719 * t586 - t916;
t581 = t586 ^ 2;
t620 = -t466 * t570 - t468 * t569;
t836 = t589 / 0.2e1;
t839 = -t586 / 0.2e1;
t138 = t227 * t589 + t885;
t139 = t230 * t586 + t916;
t896 = t138 * t838 + t139 * t837;
t910 = t589 * t620;
t604 = (t311 * t589 + t312 * t586) * t839 + ((t312 - t910) * t586 + (t314 - t878 + (t465 + t620) * t586 + t311) * t589 + t51) * t838 + (t314 * t586 + t581 * t465 + t52 + (t312 + (t465 - t620) * t589 + t910) * t589) * t836 - t896;
t912 = -t384 * t569 + t459 * t743;
t582 = t589 ^ 2;
t683 = t581 + t582;
t911 = t643 * t683;
t535 = -Icges(5,2) * t569 + t781;
t909 = t535 + t639;
t907 = -t569 / 0.2e1;
t841 = t569 / 0.2e1;
t906 = -t570 / 0.2e1;
t905 = t570 / 0.2e1;
t904 = t586 / 0.4e1;
t903 = t589 / 0.4e1;
t902 = t139 + t52;
t630 = Icges(6,5) * t587 - Icges(6,6) * t584;
t453 = Icges(6,3) * t569 + t570 * t630;
t278 = t453 * t743 + t521 * t455 - t522 * t457;
t901 = t278 * t569;
t783 = Icges(4,4) * t588;
t546 = -Icges(4,2) * t585 + t783;
t640 = Icges(4,1) * t585 + t783;
t897 = (t640 / 0.2e1 + t546 / 0.2e1) * t588;
t420 = t455 * t586;
t422 = t457 * t586;
t626 = -t375 * t584 + t378 * t587;
t609 = t453 * t586 - t626;
t183 = (-t420 * t584 + t422 * t587 + t372) * t570 + t609 * t569;
t454 = Icges(6,6) * t570 - t569 * t635;
t456 = Icges(6,5) * t570 - t569 * t638;
t452 = Icges(6,3) * t570 - t569 * t630;
t762 = t457 * t587;
t763 = t455 * t584;
t622 = t762 - t763;
t607 = t452 - t622;
t764 = t453 * t569;
t872 = t570 * t607 - t764;
t205 = t454 * t519 + t456 * t520 - t586 * t872;
t895 = t183 + t205;
t421 = t455 * t589;
t423 = t457 * t589;
t608 = -t453 * t589 + t908;
t184 = (t421 * t584 - t423 * t587 + t374) * t570 + t608 * t569;
t206 = t521 * t454 - t522 * t456 + t589 * t872;
t894 = t184 + t206;
t391 = t697 * t589;
t394 = t603 + t606;
t395 = -t568 + (rSges(5,3) + pkin(1)) * t589 + (t643 + t653) * t586;
t646 = -t394 * t758 + t395 * t757;
t714 = t390 * t322 + t388 * t892;
t724 = (-t354 * t391 + t355 * t389 + t714) * t855 + (t353 * t539 + t646) * t856;
t725 = (-t363 * t893 + t365 * t700 + t714) * t855 + (t461 * t511 - t509 * t645 + t646) * t856;
t891 = t724 - t725;
t693 = t535 * t589 + t469;
t694 = -Icges(5,2) * t751 + t468 + t553;
t537 = Icges(5,1) * t570 - t782;
t695 = -t537 * t589 + t467;
t696 = -t537 * t586 + t466;
t890 = -(t586 * t695 - t589 * t696) * t569 + (t586 * t693 - t589 * t694) * t570;
t493 = -Icges(4,5) * t586 + t589 * t640;
t689 = t546 * t589 + t493;
t564 = Icges(4,4) * t734;
t492 = Icges(4,1) * t741 + Icges(4,5) * t589 + t564;
t690 = -Icges(4,2) * t741 + t492 + t564;
t784 = Icges(4,4) * t585;
t637 = Icges(4,2) * t588 + t784;
t491 = -Icges(4,6) * t586 + t589 * t637;
t548 = Icges(4,1) * t588 - t784;
t691 = -t548 * t589 + t491;
t490 = Icges(4,6) * t589 + t586 * t637;
t692 = -t548 * t586 + t490;
t889 = -(t586 * t691 - t589 * t692) * t585 + (t586 * t689 - t589 * t690) * t588;
t217 = (-t454 * t584 + t456 * t587 + t453) * t570 + t607 * t569;
t767 = t372 * t569;
t246 = t570 * t626 + t767;
t276 = -t453 * t745 + t455 * t519 + t457 * t520;
t299 = t570 * t622 + t764;
t887 = t217 * t841 + t299 * t905 + (t246 + t276) * t751 / 0.4e1 - (-t247 + t278) * t750 / 0.4e1;
t884 = t228 * t589;
t731 = t588 * t491;
t879 = (t585 * t493 + t731) * t589;
t256 = t586 * t702 + t707;
t119 = t256 * t281 + t389 * t388 - t391 * t390;
t361 = -t470 * t586 - t445;
t249 = t361 * t393 - t539 * t911;
t875 = m(5) * t249 + m(6) * t119;
t399 = Icges(6,5) * t519 - Icges(6,6) * t520;
t704 = -Icges(6,2) * t520 + t378 + t495;
t706 = -Icges(6,1) * t519 + t375 + t780;
t173 = -t399 * t745 + t519 * t704 - t520 * t706;
t400 = Icges(6,5) * t521 + Icges(6,6) * t522;
t703 = Icges(6,2) * t522 - t379 + t496;
t705 = -Icges(6,1) * t521 + t377 - t497;
t174 = -t400 * t745 + t519 * t703 - t520 * t705;
t86 = t173 * t589 + t174 * t586;
t175 = t399 * t743 + t521 * t704 + t522 * t706;
t176 = t400 * t743 + t521 * t703 + t522 * t705;
t87 = t175 * t589 + t176 * t586;
t795 = t86 * t836 + t87 * t838;
t392 = -t509 * t589 + t586 * t511;
t711 = (-t586 * t893 + t589 * t700) * t855 + t392 * t856;
t871 = t570 * t608 - t765;
t870 = t570 * t609 - t767;
t869 = t569 * t909 - t570 * (-t636 + t537);
t698 = t457 + t501;
t699 = t455 - t504;
t233 = -t498 * t745 + t519 * t698 - t520 * t699;
t234 = t498 * t743 + t521 * t698 + t522 * t699;
t194 = t400 * t569 + (-t584 * t703 - t587 * t705) * t570;
t771 = t194 * t586;
t193 = t399 * t569 + (-t584 * t704 - t587 * t706) * t570;
t772 = t193 * t589;
t649 = t772 / 0.4e1 + t771 / 0.4e1 + t234 * t904 + t233 * t903;
t272 = t276 * t569;
t35 = t272 + (-t722 * t586 + t884) * t570;
t36 = -t901 + (t719 * t589 + t915) * t570;
t628 = t230 * t589 - t915;
t107 = t570 * t628 + t901;
t730 = t589 * t107;
t629 = -t227 * t586 + t884;
t106 = t570 * t629 + t272;
t740 = t586 * t106;
t863 = t730 / 0.4e1 - t740 / 0.4e1 + t36 * t903 + t35 * t904;
t597 = t456 * t744 / 0.2e1 - t454 * t746 / 0.2e1 + t453 * t905 + (t762 + t537) * t907 + t909 * t906 + (t763 + t452 + t636) * t841;
t861 = 0.4e1 * qJD(1);
t860 = 2 * qJD(3);
t858 = 2 * qJD(4);
t854 = pkin(1) + pkin(6);
t623 = t382 * t589 + t384 * t586;
t283 = t623 * t570;
t323 = t382 * t569 + t570 * t440;
t668 = -t283 * t281 + t323 * t390 + t388 * t912;
t218 = (-t424 * t589 - t425 * t586) * t570 + t623 * t569;
t250 = (t382 + t439) * t570 + (t424 - t440) * t569;
t251 = (t458 * t589 - t384) * t570 + (-t459 * t589 - t425) * t569;
t669 = t218 * t236 - t250 * t365 + t251 * t363;
t852 = m(6) * (t668 + t669);
t40 = t218 * t256 - t250 * t391 + t251 * t389 + t668;
t851 = m(6) * t40;
t847 = m(6) * (-t218 * t283 + t250 * t323 + t251 * t912);
t327 = -t586 * t405 + t406 * t589;
t510 = (-rSges(6,1) * t584 - rSges(6,2) * t587) * t570;
t761 = t510 * t586;
t716 = t256 * t327 + t389 * t761;
t718 = t236 * t327 + t363 * t761;
t760 = t510 * t589;
t846 = m(6) * ((t365 + t391) * t760 + t716 + t718);
t717 = t250 * t322 + t251 * t892;
t845 = m(6) * (t323 * t354 + t355 * t912 + t717);
t844 = m(6) * (-t323 * t700 - t893 * t912 + t717);
t835 = m(3) * ((rSges(3,3) * t589 - t586 * pkin(1) + t571) * t589 + (t799 + (rSges(3,3) + qJ(2)) * t586) * t586);
t644 = rSges(4,1) * t585 + rSges(4,2) * t588;
t602 = -t586 * rSges(4,3) + t589 * t644;
t427 = -t854 * t586 + t571 + t602;
t428 = (rSges(4,3) + t854) * t589 + (qJ(2) + t644) * t586;
t834 = m(4) * (t427 * t532 + t428 * t531);
t833 = m(4) * (t427 * t589 + t428 * t586);
t825 = m(5) * (t394 * t645 + t395 * t450);
t824 = m(5) * (t394 * t511 + t395 * t509);
t823 = m(5) * (t394 * t589 + t395 * t586);
t647 = -t322 * t760 + t761 * t892;
t813 = m(6) * (-t363 * t406 - t365 * t405 + t647);
t812 = m(6) * (-t389 * t406 - t391 * t405 + t647);
t810 = m(6) * (t322 * t354 + t355 * t892);
t809 = m(6) * (-t322 * t700 - t892 * t893);
t808 = m(6) * (t322 * t586 + t589 * t892);
t807 = m(6) * (t323 * t586 + t589 * t912);
t802 = m(6) * (t389 * t589 - t391 * t586);
t326 = -t405 * t589 - t586 * t406;
t800 = m(6) * t326;
t769 = t247 * t589;
t627 = -t246 * t586 - t769;
t118 = t299 * t569 + t570 * t627;
t164 = t420 * t519 + t422 * t520 - t586 * t870;
t165 = -t421 * t519 - t423 * t520 - t586 * t871;
t27 = (-t164 * t586 + t165 * t589 + t276) * t570 + (t205 - t629) * t569;
t166 = t521 * t420 - t522 * t422 + t589 * t870;
t167 = -t521 * t421 + t522 * t423 + t589 * t871;
t28 = (-t166 * t586 + t167 * t589 + t278) * t570 + (t206 - t628) * t569;
t39 = (-t183 * t586 + t184 * t589 + t299) * t570 + (t217 - t627) * t569;
t13 = t847 + (t27 * t839 + t28 * t836 + t118 / 0.2e1) * t570 + (t740 / 0.2e1 - t730 / 0.2e1 + t39 / 0.2e1) * t569;
t605 = m(6) * (-t250 * t589 + t251 * t586);
t614 = t683 * t510 * t855;
t151 = -t605 / 0.2e1 + t614;
t682 = t151 * qJD(2);
t785 = t13 * qJD(5) - t682;
t774 = t183 * t589;
t773 = t184 * t586;
t732 = t588 * t490;
t150 = t605 / 0.2e1 + t614;
t259 = m(6) * (t388 * t586 - t390 * t589) - m(5) * t911;
t723 = t259 * qJD(4) + t150 * qJD(5);
t681 = qJD(3) + qJD(4);
t680 = t846 / 0.2e1 + t795;
t671 = -t106 / 0.2e1 + t35 / 0.2e1;
t670 = -t36 / 0.2e1 - t107 / 0.2e1;
t633 = Icges(4,5) * t585 + Icges(4,6) * t588;
t488 = Icges(4,3) * t589 + t586 * t633;
t328 = t589 * t488 + t492 * t741 + t586 * t732;
t489 = -Icges(4,3) * t586 + t589 * t633;
t329 = -t589 * t489 - t493 * t741 - t586 * t731;
t662 = -t745 / 0.2e1;
t661 = -t745 / 0.4e1;
t660 = t745 / 0.4e1;
t659 = -t743 / 0.4e1;
t658 = t743 / 0.2e1;
t657 = t743 / 0.4e1;
t651 = t683 * t644;
t632 = Icges(5,5) * t570 - Icges(5,6) * t569;
t499 = t586 * t632;
t500 = t632 * t589;
t75 = t164 * t589 + t165 * t586;
t76 = t166 * t589 + t167 * t586;
t650 = (t76 - t581 * t500 + (t586 * t499 + t890) * t589) * t838 + (t75 + t582 * t499 + (-t589 * t500 - t890) * t586) * t836;
t648 = t683 * t797;
t634 = Icges(4,5) * t588 - Icges(4,6) * t585;
t460 = (-t643 - t798) * t586;
t462 = t566 + t757;
t621 = t460 * t586 - t462 * t589;
t618 = t467 * t569 - t469 * t570;
t617 = -t585 * t492 - t732;
t612 = -t283 * t327 - t323 * t760 + t761 * t912;
t601 = t138 * t659 + t51 * t657 + t902 * t661 + t863;
t600 = t27 * t836 + t28 * t838 + t76 * t658 + t75 * t662 + (t246 * t589 - t247 * t586) * t905 + (t773 + t774) * t841 - t795 + t896 * t569;
t598 = t894 * t657 + t895 * t661 + t887;
t61 = t233 * t569 + (-t173 * t586 + t174 * t589) * t570;
t62 = t234 * t569 + (-t175 * t586 + t176 * t589) * t570;
t596 = t118 * t906 + t27 * t745 / 0.2e1 - t28 * t743 / 0.2e1 + t62 * t838 + t61 * t836 - t847 + t86 * t662 + t87 * t658 + (t740 + t39) * t907 + (t730 + t771 + t772) * t841;
t594 = t774 / 0.2e1 + t773 / 0.2e1 + (t569 * t693 + t570 * t695 + t589 * t869 + t206 - t899) * t838 + (-t569 * t694 - t570 * t696 - t586 * t869 + t205 - t898) * t836 - t604;
t593 = -t597 - t769 / 0.2e1 + t618 * t837 + (t618 + t247) * t836;
t592 = t138 * t657 + t51 * t659 + t902 * t660 + t598 + t649 - t863;
t591 = t894 * t659 + t895 * t660 + t601 + t649 - t887;
t590 = t598 + t601 - t649;
t526 = t634 * t589;
t525 = t586 * t634;
t472 = t586 * t488;
t364 = t566 + t390;
t362 = -t679 + t388;
t352 = t393 - t648;
t335 = -t569 * t406 + t510 * t743;
t334 = t405 * t569 + t510 * t745;
t331 = -t586 * t489 + t879;
t330 = t589 * t617 + t472;
t319 = t800 / 0.2e1;
t308 = t326 * t570;
t296 = t802 / 0.2e1;
t269 = -t648 + t281;
t254 = (t754 + (-t584 * t698 - t587 * t699) * t570) * t569;
t253 = t330 * t589 + t331 * t586;
t252 = t328 * t589 + t329 * t586;
t237 = t807 / 0.2e1;
t155 = t391 * t760 + t716;
t154 = t296 - t711;
t153 = -t802 / 0.2e1 + t711;
t152 = t296 + t711;
t143 = t151 * qJD(5);
t131 = t812 / 0.2e1;
t123 = t813 / 0.2e1;
t121 = t808 + t823 + t833 + t835;
t120 = t365 * t760 + t718;
t116 = t581 * t489 + (t329 - t472 + (t489 - t617) * t589) * t589;
t115 = (-t330 + t472 + t329) * t586 + (t331 - t879 + (t489 + t617) * t586 + t328) * t589;
t90 = t667 + t713;
t79 = t237 - t800 / 0.2e1;
t78 = t319 + t237;
t77 = t319 - t807 / 0.2e1;
t70 = t597 + t809 + t824;
t68 = t844 / 0.2e1;
t66 = t845 / 0.2e1;
t57 = t834 + t825 + t810 + t597 - t897 + (-t548 / 0.2e1 + t637 / 0.2e1) * t585;
t38 = t851 / 0.2e1;
t29 = t852 / 0.2e1;
t22 = m(6) * t155 + t795;
t21 = m(6) * t120 + t795;
t20 = t650 + t875;
t19 = t20 * qJD(4);
t18 = t650 - t918;
t17 = t38 - t852 / 0.2e1 + t680;
t16 = t29 - t851 / 0.2e1 + t680;
t15 = (t586 * t670 + t589 * t671) * t570;
t11 = (t116 / 0.2e1 + t253 / 0.2e1) * t589 + (-t252 / 0.2e1 + t115 / 0.2e1) * t586 + t604;
t10 = t604 + t891;
t9 = t604 - t891;
t8 = t594 + t724 + t725;
t7 = t29 + t38 - t846 / 0.2e1 + t600;
t6 = t68 + t131 + t592;
t5 = -t844 / 0.2e1 + t591 + t131;
t4 = -t812 / 0.2e1 + t68 + t590;
t3 = t66 + t123 + t592;
t2 = -t845 / 0.2e1 + t591 + t123;
t1 = -t813 / 0.2e1 + t66 + t590;
t12 = [t121 * qJD(2) + t57 * qJD(3) + t70 * qJD(4) + t122 * qJD(5), qJD(1) * t121 + qJD(3) * t90 + qJD(4) * t152 + qJD(5) * t78, t57 * qJD(1) + t90 * qJD(2) + t8 * qJD(4) + t3 * qJD(5) + ((t407 * t551 - (t427 * t586 - t428 * t589) * t644) * t857 + (t394 * t460 + t395 * t462 + (-t450 + t461) * t645) * t856 + (t322 * t364 - t354 * t365 + t355 * t363 + t362 * t892) * t855) * t860 + ((-t585 * t690 - t588 * t692) * t836 + t594 + t115 * t839 + (t585 * t689 + t588 * t691 + t252) * t838 + (t116 + t253) * t837 - (t581 / 0.2e1 + t582 / 0.2e1) * t633) * qJD(3), t70 * qJD(1) + t152 * qJD(2) + t8 * qJD(3) + t594 * qJD(4) + t6 * qJD(5) + ((t392 * t539 + t646) * t856 + (-t389 * t893 + t391 * t700 + t714) * t855) * t858, t917 + t78 * qJD(2) + t3 * qJD(3) + t6 * qJD(4) + (m(6) * (t322 * t334 + t323 * t405 + t335 * t892 - t406 * t912) + t254 + ((t234 / 0.2e1 + t194 / 0.2e1 - t671) * t589 + (-t233 / 0.2e1 - t193 / 0.2e1 - t670) * t586) * t570) * qJD(5); -t92 * qJD(3) + t153 * qJD(4) + t77 * qJD(5) + (-t808 / 0.4e1 - t823 / 0.4e1 - t833 / 0.4e1 - t835 / 0.4e1) * t861, 0, -t914 + (t621 * t856 + (t362 * t586 - t364 * t589) * t855 - t651 * t857) * t860 + t723, qJD(1) * t153 + qJD(3) * t259 + t723, t77 * qJD(1) + m(6) * (-t334 * t589 + t335 * t586) * qJD(5) + t681 * t150; (t593 + (-t637 + t548) * t585 / 0.2e1 + t897) * qJD(1) + t92 * qJD(2) + t11 * qJD(3) + t9 * qJD(4) + t2 * qJD(5) + (-t834 / 0.4e1 - t825 / 0.4e1 - t810 / 0.4e1) * t861, t143 + t914, t11 * qJD(1) + (m(6) * (t236 * t269 + t362 * t363 - t364 * t365) + m(5) * (t307 * t352 + t460 * t461 - t462 * t645) + m(4) * ((-t589 * t602 + (-t589 * rSges(4,3) - t586 * t644) * t586) * (-t586 * t531 - t589 * t532) - t551 * t651) + (t582 * t525 + (-t589 * t526 - t889) * t586) * t836 + (-t581 * t526 + (t586 * t525 + t889) * t589) * t838 + t650) * qJD(3) + t18 * qJD(4) + t21 * qJD(5), t9 * qJD(1) + t18 * qJD(3) + t16 * qJD(5) + ((t119 + t109) * t855 + (t207 + t249) * t856) * t858 + (t650 - t875) * qJD(4), t2 * qJD(1) + t682 + t21 * qJD(3) + t16 * qJD(4) + (m(6) * (t308 * t236 - t334 * t365 + t335 * t363 + t612) + t596) * qJD(5); t593 * qJD(1) + t154 * qJD(2) + t10 * qJD(3) + t604 * qJD(4) + t5 * qJD(5) + (-t824 / 0.4e1 - t809 / 0.4e1) * t861, qJD(1) * t154 + t143, t10 * qJD(1) + t19 + t17 * qJD(5) + ((t256 * t269 + t362 * t389 - t364 * t391 + t109) * t855 + (t361 * t352 + t539 * t621 + t207) * t856) * t860 + (t650 + t918) * qJD(3), qJD(1) * t604 + qJD(3) * t20 + qJD(5) * t22 + t19, t5 * qJD(1) + t682 + t17 * qJD(3) + t22 * qJD(4) + (m(6) * (t308 * t256 - t334 * t391 + t335 * t389 + t612) + t596) * qJD(5); t79 * qJD(2) + t1 * qJD(3) + t4 * qJD(4) + t15 * qJD(5) - t917, qJD(1) * t79 - t151 * t681, t1 * qJD(1) + ((-t269 * t283 + t323 * t364 + t362 * t912 - t120 + t669) * m(6) + t600) * qJD(3) + t7 * qJD(4) + t785, t4 * qJD(1) + t7 * qJD(3) + ((t40 - t155) * m(6) + t600) * qJD(4) + t785, t15 * qJD(1) + (m(6) * (-t283 * t308 + t323 * t334 + t335 * t912) + t254 * t841 + (t61 * t839 + t62 * t836 + (-t193 * t586 + t194 * t589) * t841) * t570) * qJD(5) + t681 * t13;];
Cq = t12;
