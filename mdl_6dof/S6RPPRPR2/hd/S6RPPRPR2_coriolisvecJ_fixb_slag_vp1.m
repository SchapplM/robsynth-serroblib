% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:34
% EndTime: 2019-03-09 01:42:21
% DurationCPUTime: 44.17s
% Computational Cost: add. (26578->910), mult. (24595->1169), div. (0->0), fcn. (21822->10), ass. (0->466)
t788 = Icges(6,4) - Icges(5,5);
t787 = Icges(6,5) - Icges(5,6);
t786 = Icges(6,1) + Icges(5,3);
t384 = pkin(10) + qJ(4);
t379 = sin(t384);
t381 = cos(t384);
t748 = t787 * t379 - t381 * t788;
t385 = qJ(1) + pkin(9);
t382 = cos(t385);
t785 = t786 * t382;
t380 = sin(t385);
t624 = t380 * t381;
t626 = t379 * t380;
t749 = t624 * t788 - t787 * t626 + t785;
t756 = t380 * t786 + t748 * t382;
t291 = Icges(5,5) * t379 + Icges(5,6) * t381;
t476 = Icges(6,4) * t379 + Icges(6,5) * t381;
t777 = t291 - t476;
t652 = Icges(5,4) * t379;
t295 = Icges(5,2) * t381 + t652;
t641 = Icges(6,6) * t379;
t470 = Icges(6,3) * t381 + t641;
t776 = -t295 - t470;
t367 = Icges(5,4) * t381;
t297 = Icges(5,1) * t379 + t367;
t640 = Icges(6,6) * t381;
t472 = Icges(6,2) * t379 + t640;
t784 = t297 + t472;
t298 = Icges(5,1) * t381 - t652;
t473 = Icges(6,2) * t381 - t641;
t783 = t298 + t473;
t471 = -Icges(6,3) * t379 + t640;
t478 = -Icges(5,2) * t379 + t367;
t782 = t471 + t478;
t642 = Icges(5,6) * t382;
t195 = Icges(5,4) * t624 - Icges(5,2) * t626 - t642;
t327 = Icges(6,6) * t626;
t650 = Icges(6,4) * t382;
t202 = Icges(6,2) * t624 - t327 + t650;
t781 = t195 * t379 - t202 * t381;
t198 = Icges(5,5) * t380 + t298 * t382;
t199 = Icges(6,5) * t380 - t382 * t471;
t780 = -t198 * t624 - t199 * t626;
t458 = t295 * t379 - t297 * t381;
t771 = -t379 * t470 + t381 * t472 - t458;
t332 = Icges(5,4) * t626;
t646 = Icges(5,5) * t382;
t197 = Icges(5,1) * t624 - t332 - t646;
t645 = Icges(6,5) * t382;
t200 = Icges(6,6) * t624 - Icges(6,3) * t626 + t645;
t753 = -t197 * t381 + t200 * t379 + t781;
t196 = Icges(5,6) * t380 + t382 * t478;
t764 = t196 - t199;
t779 = t782 * qJD(4);
t778 = t783 * qJD(4);
t773 = t776 * qJD(4);
t772 = t784 * qJD(4);
t769 = t382 * t756 + t780;
t621 = t381 * t382;
t625 = t379 * t382;
t728 = t198 * t621 + t199 * t625 + t380 * t756;
t768 = -t197 * t621 + t200 * t625 + t380 * t749;
t712 = t777 * t380;
t767 = -t380 * t753 + t382 * t749;
t328 = Icges(6,6) * t625;
t651 = Icges(6,4) * t380;
t201 = -Icges(6,2) * t621 + t328 + t651;
t738 = -t196 * t626 - t201 * t624 - t769;
t737 = -t195 * t625 + t202 * t621 - t768;
t736 = -t196 * t625 - t201 * t621 + t728;
t766 = t382 * t771 + t712;
t765 = t195 + t200;
t763 = t197 + t202;
t762 = t198 - t201;
t388 = -pkin(7) - qJ(3);
t353 = t382 * t388;
t387 = cos(pkin(10));
t376 = pkin(3) * t387 + pkin(2);
t585 = -t380 * t376 - t353;
t390 = sin(qJ(1));
t680 = pkin(1) * t390;
t501 = t585 - t680;
t761 = t196 * t379 + t201 * t381;
t760 = t773 * t382 + (-t380 * t782 + t642 - t645) * qJD(1);
t759 = qJD(1) * t764 + t773 * t380;
t758 = -t772 * t382 + (-t380 * t783 + t646 - t650) * qJD(1);
t757 = t772 * t380 + (-t382 * t473 - t198 + t651) * qJD(1);
t755 = t778 * t381 - t779 * t379 + (-t379 * t784 + t381 * t776) * qJD(4) + t777 * qJD(1);
t556 = rSges(5,1) * t624;
t754 = -t556 + t501;
t752 = t198 * t381 + t199 * t379 - t761;
t751 = t771 * qJD(1) - qJD(4) * t748;
t735 = t379 * t763 + t381 * t765;
t734 = t379 * t762 + t381 * t764;
t750 = t766 * qJD(1);
t747 = t777 * qJD(4);
t746 = (t738 * t380 - t382 * t767) * qJD(4);
t745 = (t380 * t736 - t382 * t737) * qJD(4);
t241 = t476 * t382;
t564 = (t470 * t626 - t472 * t624 - t241) * qJD(1);
t627 = t291 * t382;
t106 = -t380 * t458 - t627;
t565 = t106 * qJD(1);
t744 = t565 - t564 + t746;
t743 = t745 + t750;
t742 = qJD(4) * t753 + t379 * t757 - t381 * t759;
t741 = t752 * qJD(4) + t758 * t379 + t381 * t760;
t740 = -t380 * t751 + t382 * t755;
t739 = t380 * t755 + t382 * t751;
t733 = -t380 * t764 + t382 * t765;
t732 = t749 + t761;
t731 = t756 * qJD(1);
t730 = t782 + t784;
t729 = t783 + t776;
t727 = (t327 + t332 + (Icges(5,2) + Icges(6,3)) * t624 - t763) * t382 + (-Icges(6,3) * t621 - t295 * t382 - t328 + t762) * t380;
t726 = t749 * qJD(1) + t735 * qJD(4) + t379 * t759 + t381 * t757;
t725 = -t734 * qJD(4) - t379 * t760 + t758 * t381 + t731;
t724 = -t747 * t382 + (-t380 * t748 - t752 + t785) * qJD(1);
t723 = qJD(1) * t753 - t747 * t380 + t731;
t566 = qJD(6) * t381;
t572 = qJD(4) * t380;
t273 = t382 * t566 + t572;
t570 = qJD(4) * t382;
t274 = -t380 * t566 + t570;
t567 = qJD(6) * t379;
t349 = qJD(1) + t567;
t391 = cos(qJ(6));
t619 = t382 * t391;
t389 = sin(qJ(6));
t623 = t380 * t389;
t264 = t379 * t619 - t623;
t620 = t382 * t389;
t622 = t380 * t391;
t265 = t379 * t620 + t622;
t139 = Icges(7,5) * t265 + Icges(7,6) * t264 + Icges(7,3) * t621;
t649 = Icges(7,4) * t265;
t142 = Icges(7,2) * t264 + Icges(7,6) * t621 + t649;
t246 = Icges(7,4) * t264;
t145 = Icges(7,1) * t265 + Icges(7,5) * t621 + t246;
t39 = t139 * t621 + t264 * t142 + t265 * t145;
t266 = t379 * t622 + t620;
t267 = -t379 * t623 + t619;
t141 = -Icges(7,5) * t267 + Icges(7,6) * t266 + Icges(7,3) * t624;
t248 = Icges(7,4) * t267;
t144 = Icges(7,2) * t266 + Icges(7,6) * t624 - t248;
t247 = Icges(7,4) * t266;
t146 = Icges(7,1) * t267 - Icges(7,5) * t624 - t247;
t40 = t141 * t621 + t264 * t144 - t146 * t265;
t474 = Icges(7,5) * t389 + Icges(7,6) * t391;
t416 = -Icges(7,3) * t379 + t381 * t474;
t648 = Icges(7,4) * t389;
t475 = Icges(7,2) * t391 + t648;
t417 = -Icges(7,6) * t379 + t381 * t475;
t647 = Icges(7,4) * t391;
t479 = Icges(7,1) * t389 + t647;
t418 = -Icges(7,5) * t379 + t381 * t479;
t70 = -t264 * t417 - t265 * t418 - t416 * t621;
t10 = t273 * t39 - t274 * t40 + t70 * t349;
t41 = t139 * t624 + t266 * t142 - t267 * t145;
t42 = t141 * t624 + t144 * t266 + t146 * t267;
t71 = -t266 * t417 + t267 * t418 - t416 * t624;
t11 = t273 * t41 - t274 * t42 + t349 * t71;
t468 = t144 * t391 - t146 * t389;
t53 = t141 * t379 - t381 * t468;
t394 = qJD(1) ^ 2;
t722 = 0.2e1 * qJD(4);
t718 = -t379 * t727 + t381 * t733;
t717 = (-t379 * t730 + t381 * t729) * qJD(1);
t671 = rSges(4,2) * sin(pkin(10));
t673 = rSges(4,1) * t387;
t451 = t380 * rSges(4,3) + (-t671 + t673) * t382;
t362 = t380 * qJ(3);
t307 = t382 * pkin(2) + t362;
t392 = cos(qJ(1));
t383 = t392 * pkin(1);
t523 = t307 + t383;
t716 = t451 + t523;
t543 = t379 * t570;
t575 = qJD(1) * t380;
t547 = t381 * t575;
t715 = t543 + t547;
t714 = t748 * qJD(1);
t713 = t627 - t241;
t490 = rSges(7,1) * t267 - rSges(7,2) * t266;
t150 = rSges(7,3) * t624 - t490;
t489 = rSges(7,1) * t389 + rSges(7,2) * t391;
t421 = -rSges(7,3) * t379 + t381 * t489;
t299 = pkin(4) * t379 - qJ(5) * t381;
t522 = -pkin(8) * t379 - t299;
t495 = t522 * t382;
t568 = qJD(5) * t382;
t323 = t379 * t568;
t360 = qJD(3) * t380;
t586 = t323 + t360;
t711 = qJD(4) * t495 - t150 * t349 + t274 * t421 + t586;
t374 = t380 * pkin(5);
t359 = qJD(5) * t379;
t659 = t379 * rSges(6,2);
t488 = rSges(6,3) * t381 + t659;
t709 = (-qJD(4) * t488 - t359) * t380;
t363 = t382 * qJ(3);
t302 = pkin(2) * t380 - t363;
t188 = t302 + t585;
t286 = qJD(1) * t302;
t708 = qJD(1) * t188 - t286;
t707 = -rSges(5,2) * t626 - t382 * rSges(5,3);
t518 = t382 * rSges(3,1) - rSges(3,2) * t380;
t341 = pkin(8) * t621;
t276 = t374 + t341;
t706 = t383 + t518;
t560 = -rSges(7,3) - pkin(4) - pkin(8);
t636 = qJ(5) * t379;
t705 = t560 * t381 - t376 - t636;
t340 = pkin(8) * t624;
t277 = pkin(5) * t382 - t340;
t304 = pkin(4) * t381 + t636;
t254 = t304 * t380;
t524 = -t302 - t680;
t503 = t188 + t524;
t456 = -t254 + t503;
t37 = (t277 + t456) * qJD(1) + t711;
t148 = t265 * rSges(7,1) + t264 * rSges(7,2) + rSges(7,3) * t621;
t544 = t379 * t572;
t317 = pkin(8) * t544;
t258 = pkin(4) * t621 + qJ(5) * t625;
t334 = t382 * t376;
t511 = -t380 * t388 + t334;
t502 = t511 - t307 + t523;
t455 = t258 + t502;
t539 = t380 * t359;
t268 = t299 * t572;
t361 = qJD(3) * t382;
t594 = -t268 - t361;
t38 = t539 + t148 * t349 + t421 * t273 - t317 + (t276 + t455) * qJD(1) + t594;
t704 = t37 * t382 + t38 * t380;
t228 = qJD(1) * t254;
t701 = -t228 + t708;
t271 = (Icges(7,2) * t389 - t647) * t381;
t411 = t273 * (-Icges(7,2) * t265 + t145 + t246) - t274 * (Icges(7,2) * t267 - t146 + t247) + t349 * (-t418 + t271);
t272 = (-Icges(7,1) * t391 + t648) * t381;
t412 = t273 * (-Icges(7,1) * t264 + t142 + t649) - t274 * (-Icges(7,1) * t266 + t144 - t248) + t349 * (-t417 - t272);
t692 = m(6) / 0.2e1;
t691 = m(7) / 0.2e1;
t561 = qJD(4) * qJD(6);
t534 = t379 * t561;
t186 = qJD(1) * t273 - t380 * t534;
t690 = t186 / 0.2e1;
t187 = qJD(1) * t274 - t382 * t534;
t689 = t187 / 0.2e1;
t688 = -t273 / 0.2e1;
t687 = t273 / 0.2e1;
t686 = -t274 / 0.2e1;
t685 = t274 / 0.2e1;
t684 = -t349 / 0.2e1;
t683 = t349 / 0.2e1;
t682 = t380 / 0.2e1;
t681 = -t382 / 0.2e1;
t469 = t142 * t391 + t145 * t389;
t571 = qJD(4) * t381;
t419 = -t349 * t389 + t391 * t571;
t508 = qJD(1) * t379 + qJD(6);
t454 = t380 * t508;
t121 = t382 * t419 - t391 * t454;
t420 = t349 * t391 + t389 * t571;
t122 = t382 * t420 - t389 * t454;
t61 = Icges(7,5) * t122 + Icges(7,6) * t121 - Icges(7,3) * t715;
t63 = Icges(7,4) * t122 + Icges(7,2) * t121 - Icges(7,6) * t715;
t65 = Icges(7,1) * t122 + Icges(7,4) * t121 - Icges(7,5) * t715;
t8 = (qJD(4) * t469 + t61) * t379 + (qJD(4) * t139 - t389 * t65 - t391 * t63 + (t142 * t389 - t145 * t391) * qJD(6)) * t381;
t679 = t8 * t273;
t453 = t382 * t508;
t119 = t380 * t419 + t391 * t453;
t120 = t380 * t420 + t389 * t453;
t574 = qJD(1) * t382;
t546 = t381 * t574;
t430 = -t544 + t546;
t60 = Icges(7,5) * t120 + Icges(7,6) * t119 + Icges(7,3) * t430;
t62 = Icges(7,4) * t120 + Icges(7,2) * t119 + Icges(7,6) * t430;
t64 = Icges(7,1) * t120 + Icges(7,4) * t119 + Icges(7,5) * t430;
t9 = (qJD(4) * t468 + t60) * t379 + (qJD(4) * t141 - t389 * t64 - t391 * t62 + (t144 * t389 + t146 * t391) * qJD(6)) * t381;
t678 = t9 * t274;
t676 = qJD(1) / 0.2e1;
t675 = pkin(2) - t376;
t219 = Icges(7,3) * t381 + t379 * t474;
t270 = (-Icges(7,5) * t391 + Icges(7,6) * t389) * t381;
t153 = qJD(4) * t219 + qJD(6) * t270;
t221 = Icges(7,6) * t381 + t379 * t475;
t154 = qJD(4) * t221 + qJD(6) * t271;
t223 = Icges(7,5) * t381 + t379 * t479;
t155 = qJD(4) * t223 + qJD(6) * t272;
t460 = -t389 * t418 - t391 * t417;
t27 = (qJD(4) * t460 + t153) * t379 + (-qJD(4) * t416 - t154 * t391 - t155 * t389 + (-t389 * t417 + t391 * t418) * qJD(6)) * t381;
t533 = t381 * t561;
t85 = -t379 * t416 - t381 * t460;
t674 = t27 * t349 + t85 * t533;
t670 = rSges(6,2) * t381;
t669 = rSges(6,3) * t379;
t667 = rSges(7,3) * t381;
t225 = t379 * t489 + t667;
t275 = (-rSges(7,1) * t391 + rSges(7,2) * t389) * t381;
t156 = qJD(4) * t225 + qJD(6) * t275;
t358 = pkin(5) * t574;
t191 = -pkin(8) * t715 + t358;
t569 = qJD(5) * t381;
t249 = qJD(4) * t304 - t569;
t393 = qJD(4) ^ 2;
t542 = t381 * t570;
t311 = qJ(5) * t542;
t548 = t379 * t575;
t117 = -pkin(4) * t715 - qJ(5) * t548 + t311 + t323;
t354 = qJ(3) * t574;
t559 = t394 * t680;
t563 = qJD(1) * qJD(3);
t581 = t354 + t360;
t457 = qJD(1) * (-pkin(2) * t575 + t581) + t380 * t563 - t559;
t445 = qJD(1) * (-t354 + (t380 * t675 - t353) * qJD(1)) + t457;
t562 = qJD(4) * qJD(5);
t535 = t381 * t562;
t422 = t380 * t535 + t445 + (t117 + t323) * qJD(1);
t538 = t148 * t566;
t618 = t122 * rSges(7,1) + t121 * rSges(7,2);
t67 = -rSges(7,3) * t715 + t618;
t12 = -t393 * t340 + qJD(1) * t191 - t156 * t273 + t187 * t421 + t349 * t67 + (qJD(1) * t495 - t249 * t380 + t538) * qJD(4) + t422;
t666 = t12 * t382;
t190 = qJD(1) * t276 - t317;
t558 = t394 * t383;
t504 = t382 * t563 - t558;
t444 = qJD(1) * t268 + t382 * t535 + t504;
t537 = t150 * t566;
t232 = t379 * t574 + t380 * t571;
t318 = pkin(4) * t544;
t118 = pkin(4) * t546 + qJ(5) * t232 - t318 + t539;
t250 = qJD(1) * t307 - t361;
t344 = t388 * t575;
t607 = t344 - (-t382 * t675 - t362) * qJD(1) - t250;
t552 = -t118 + t607;
t491 = rSges(7,1) * t120 + rSges(7,2) * t119;
t66 = rSges(7,3) * t430 + t491;
t13 = -t393 * t341 - t156 * t274 - t186 * t421 - t349 * t66 + (-t249 * t382 - t537) * qJD(4) + (-t190 + (pkin(8) * qJD(4) - qJD(5)) * t626 + t552) * qJD(1) + t444;
t665 = t13 * t380;
t301 = rSges(5,1) * t379 + rSges(5,2) * t381;
t257 = t301 * t382;
t368 = t380 * rSges(5,3);
t213 = rSges(5,1) * t621 - rSges(5,2) * t625 + t368;
t545 = t301 * t572;
t83 = -t545 - t361 + (t213 + t502) * qJD(1);
t664 = t257 * t83;
t370 = t380 * rSges(6,1);
t212 = t556 + t707;
t540 = t301 * t570;
t496 = t360 - t540;
t82 = (-t212 + t503) * qJD(1) + t496;
t657 = t380 * t82;
t52 = t139 * t379 - t381 * t469;
t656 = t52 * t187;
t655 = t53 * t186;
t448 = t254 * t572 + t258 * t570 + qJD(2) - t569;
t34 = t148 * t274 + t150 * t273 + (t276 * t382 - t277 * t380) * qJD(4) + t448;
t635 = qJD(4) * t34;
t613 = t148 + t276;
t612 = t150 - t277;
t214 = -rSges(6,2) * t621 + rSges(6,3) * t625 + t370;
t601 = -t214 - t258;
t600 = t380 * t254 + t382 * t258;
t255 = t299 * t382;
t597 = -qJD(1) * t255 + t380 * t569;
t305 = t669 - t670;
t596 = -t305 * qJD(4) - t249;
t595 = -t258 - t276;
t589 = -t299 + t488;
t588 = -t304 - t305;
t587 = rSges(5,2) * t548 + rSges(5,3) * t574;
t345 = t380 * t671;
t584 = rSges(4,3) * t574 + qJD(1) * t345;
t583 = t344 + t361;
t582 = t382 * rSges(4,3) + t345;
t580 = t380 ^ 2 + t382 ^ 2;
t573 = qJD(4) * t379;
t557 = t380 * t673;
t554 = t38 * t574;
t553 = t382 * t117 + t380 * t118 + t254 * t574;
t251 = t299 * t380;
t551 = -t251 * t572 - t255 * t570 + t359;
t550 = t311 + t586;
t549 = t318 + t583;
t536 = -pkin(2) - t673;
t531 = t574 / 0.2e1;
t530 = -t572 / 0.2e1;
t528 = t571 / 0.2e1;
t527 = -t570 / 0.2e1;
t526 = t570 / 0.2e1;
t519 = rSges(6,1) * t382 - rSges(6,3) * t626;
t517 = t382 * t589;
t510 = qJD(1) * t251 + t381 * t568;
t507 = t118 * t572 + t379 * t562 + (t117 + t228) * t570;
t506 = rSges(6,1) * t574 + rSges(6,2) * t715 + rSges(6,3) * t542;
t505 = t334 + t383 + t258;
t500 = t421 + t522;
t497 = qJD(6) * t528;
t303 = rSges(3,1) * t380 + rSges(3,2) * t382;
t492 = rSges(5,1) * t381 - rSges(5,2) * t379;
t487 = t380 * t40 + t382 * t39;
t486 = t380 * t39 - t382 * t40;
t485 = t380 * t42 + t382 * t41;
t484 = t380 * t41 - t382 * t42;
t483 = t380 * t53 + t382 * t52;
t482 = t380 * t52 - t382 * t53;
t481 = -t380 * t83 - t382 * t82;
t467 = t148 * t380 - t150 * t382;
t461 = t212 * t380 + t213 * t382;
t450 = -pkin(8) * t571 - t156 - t249;
t253 = t301 * t380;
t252 = t488 * t380;
t431 = qJD(4) * t517 + t586;
t428 = -t139 * t273 + t141 * t274 + t349 * t416;
t427 = (Icges(7,5) * t264 - Icges(7,6) * t265) * t273 - (Icges(7,5) * t266 + Icges(7,6) * t267) * t274 + t270 * t349;
t424 = -t304 - t376 - t669;
t423 = t381 * t427;
t135 = -rSges(5,1) * t715 - rSges(5,2) * t542 + t587;
t136 = -qJD(4) * t253 + (t382 * t492 + t368) * qJD(1);
t410 = t135 * t382 + t136 * t380 + (t212 * t382 - t213 * t380) * qJD(1);
t405 = t34 * t467 - (-t37 * t380 + t38 * t382) * t421;
t396 = (t416 * t382 + t469) * t273 - (t416 * t380 + t468) * t274 + (t219 + t460) * t349;
t395 = t396 * t381;
t285 = t492 * qJD(4);
t269 = t299 * t575;
t256 = t488 * t382;
t233 = t542 - t548;
t231 = t580 * t573;
t218 = t557 - t582;
t215 = rSges(6,2) * t624 + t519;
t179 = t421 * t382;
t178 = t421 * t380;
t177 = t418 * t382;
t176 = t418 * t380;
t175 = t417 * t382;
t174 = t417 * t380;
t164 = rSges(7,1) * t266 + rSges(7,2) * t267;
t163 = rSges(7,1) * t264 - rSges(7,2) * t265;
t152 = qJD(1) * t716 - t361;
t151 = t360 + (-t218 + t524) * qJD(1);
t138 = -rSges(6,3) * t548 + t506;
t137 = qJD(4) * t252 + (t305 * t382 + t370) * qJD(1);
t98 = (-qJD(1) * t451 - t250) * qJD(1) + t504;
t97 = qJD(1) * (-qJD(1) * t557 + t584) + t457;
t92 = qJD(4) * t461 + qJD(2);
t69 = (t214 * t382 - t215 * t380) * qJD(4) + t448;
t59 = -t709 + (t214 + t455) * qJD(1) + t594;
t58 = (t215 + t456) * qJD(1) + t431;
t55 = -t285 * t570 + (-t136 + t545 + t607) * qJD(1) + t504;
t54 = -t285 * t572 + (t135 - t540) * qJD(1) + t445;
t43 = t410 * qJD(4);
t29 = t596 * t570 + (-t137 + t552 + t709) * qJD(1) + t444;
t28 = qJD(1) * t138 + (qJD(1) * t517 + t380 * t596) * qJD(4) + t422;
t17 = (t137 * t380 + t138 * t382 + (-t215 * t382 + t380 * t601) * qJD(1)) * qJD(4) + t507;
t16 = -t121 * t417 - t122 * t418 + t153 * t621 + t154 * t264 + t155 * t265 + t416 * t715;
t15 = -t119 * t417 - t120 * t418 + t153 * t624 + t154 * t266 - t155 * t267 - t416 * t430;
t14 = t273 * t52 - t274 * t53 + t349 * t85;
t7 = -t148 * t186 + t150 * t187 + t273 * t66 + t274 * t67 + (t190 * t380 + t191 * t382 + (-t277 * t382 + t380 * t595) * qJD(1)) * qJD(4) + t507;
t6 = t121 * t144 - t122 * t146 - t141 * t715 + t264 * t62 + t265 * t64 + t60 * t621;
t5 = t121 * t142 + t122 * t145 - t139 * t715 + t264 * t63 + t265 * t65 + t61 * t621;
t4 = t119 * t144 - t120 * t146 + t141 * t430 + t266 * t62 - t267 * t64 + t60 * t624;
t3 = t119 * t142 + t120 * t145 + t139 * t430 + t266 * t63 - t267 * t65 + t61 * t624;
t2 = t16 * t349 + t186 * t40 + t187 * t39 + t273 * t5 - t274 * t6 + t533 * t70;
t1 = t15 * t349 + t186 * t42 + t187 * t41 + t273 * t3 - t274 * t4 + t533 * t71;
t18 = [m(3) * ((-t303 * t394 - t559) * t706 + (-t558 + (-0.2e1 * t518 - t383 + t706) * t394) * (-t303 - t680)) + t655 / 0.2e1 + t15 * t686 + t16 * t687 + t70 * t689 + t71 * t690 - t678 / 0.2e1 + t679 / 0.2e1 + t674 + t656 / 0.2e1 + ((t728 * t380 + ((t756 + t781) * t382 + t738 + t768 + t780) * t382) * qJD(4) + t750) * t526 + (t686 + t685) * t10 + (t771 * qJD(4) + t778 * t379 + t779 * t381) * qJD(1) + (t13 * (t277 + t490 + t501) + t12 * (t505 + t613) + (t13 * (-t304 - t667) - t12 * t388) * t380 + (t317 - t491 + t549 + (rSges(7,3) * t573 - qJ(5) * t571 - t359) * t380 + (t382 * t705 - t374 - t383) * qJD(1)) * t37 + (t560 * t543 + t358 + t37 + t550 + t618 - t701 - t711 + (t380 * t705 - t277 - t353) * qJD(1)) * t38) * m(7) + (t29 * (t501 + t519) + t58 * t549 + t28 * (t214 + t505) + t59 * (-pkin(4) * t543 + t506 + t550) + (t29 * (-t304 + t670) - t28 * t388 + (-t359 + (-t659 + (-rSges(6,3) - qJ(5)) * t381) * qJD(4)) * t58) * t380 + ((-t390 * t59 - t392 * t58) * pkin(1) + (-t58 * rSges(6,1) + t424 * t59) * t380 + (t58 * (t424 + t670) - t59 * t388) * t382) * qJD(1) - (-t58 + (t215 - t680) * qJD(1) + t431 + t701) * t59) * m(6) + (t55 * (-t707 + t754) + t54 * (t213 + t383 + t511) + (t301 * t657 - t664) * qJD(4) + (t583 + (-t368 - t383 + (-t376 - t492) * t382) * qJD(1)) * t82 + (t360 + t587 + t82 - t496 - t708 + (t212 + t680 + t754) * qJD(1)) * t83) * m(5) + (t98 * (t380 * t536 + t363 + t582 - t680) + t151 * t361 + t97 * t716 + t152 * (t581 + t584) + ((-t151 * t392 - t152 * t390) * pkin(1) + t151 * (t536 + t671) * t382 + (t151 * (-rSges(4,3) - qJ(3)) + t152 * t536) * t380) * qJD(1) - (-t151 - t286 + t360 + (-t218 - t680) * qJD(1)) * t152) * m(4) + (t740 + t741) * t572 / 0.2e1 + (0.2e1 * t564 - t565 + ((t382 * t732 - t728 + t736) * t382 + (t380 * t732 + t737 + t769) * t380) * qJD(4) + t744) * t530 + (t739 - t742 + t743) * t527 + ((t106 + t735) * t380 + (t734 + t766) * t382) * qJD(4) * t676; m(5) * t43 + m(6) * t17 + m(7) * t7; 0.2e1 * (-t666 / 0.2e1 + t665 / 0.2e1) * m(7) + 0.2e1 * (t28 * t681 + t29 * t682) * m(6) + 0.2e1 * (t54 * t681 + t55 * t682) * m(5) + 0.2e1 * (t681 * t97 + t682 * t98) * m(4); (qJD(1) * t483 + t380 * t8 - t382 * t9) * t683 + (((-t175 * t391 - t177 * t389 + t139) * t273 - (-t174 * t391 - t176 * t389 + t141) * t274 + (-t221 * t391 - t223 * t389 - t416) * t349 + t85 * qJD(6)) * t381 + (-qJD(6) * t483 + t396) * t379) * t684 + ((t175 * t266 - t177 * t267) * t273 - (t174 * t266 - t176 * t267) * t274 + (t221 * t266 - t223 * t267) * t349 + (t381 * t71 - t41 * t625) * qJD(6) + ((-qJD(6) * t42 + t428) * t379 + t395) * t380) * t685 + (qJD(1) * t485 + t3 * t380 - t382 * t4) * t686 + (qJD(1) * t487 + t380 * t5 - t382 * t6) * t687 + ((t175 * t264 + t177 * t265) * t273 - (t174 * t264 + t176 * t265) * t274 + (t221 * t264 + t223 * t265) * t349 + (t381 * t70 - t40 * t626) * qJD(6) + ((-qJD(6) * t39 + t428) * t379 + t395) * t382) * t688 + t486 * t689 + t484 * t690 + t482 * t497 - t14 * t566 / 0.2e1 - ((t379 * t733 + t381 * t727) * qJD(4) + (t729 * t379 + t730 * t381) * qJD(1)) * qJD(1) / 0.2e1 + (t742 * t382 + t741 * t380 + (t380 * t735 + t382 * t734) * qJD(1)) * t676 + ((-t572 * t713 + t714) * t380 + ((t380 * t712 + t718) * qJD(4) + t717) * t382) * t530 + ((-t570 * t712 - t714) * t382 + ((t382 * t713 + t718) * qJD(4) + t717) * t380) * t526 + (t10 * t382 + t11 * t380) * t567 / 0.2e1 + (t7 * t600 + (t12 * t500 + t38 * t450 + t7 * t612) * t380 + (t7 * t613 + (qJD(1) * t38 + t13) * t500) * t382 - t38 * (t179 * t349 - t225 * t273 + t538 + t597) - ((-t580 * t635 - t554) * pkin(8) + t405 * qJD(6)) * t379 - t704 * qJD(4) * (-pkin(8) * t381 - t304) + (t178 * t349 + t225 * t274 + t382 * t450 - t421 * t575 + t269 - t510 + t537) * t37 + (t553 + (t190 + t66 + (-t148 + t595) * qJD(1)) * t380 + (qJD(1) * t612 + t191 + t67) * t382 - t178 * t273 - t179 * t274 - t551) * t34) * m(7) + (t58 * t269 + t17 * t600 + t69 * t553 + (t29 * t589 + t58 * t596 + t17 * t214 + t69 * t138 + (-t69 * t215 + t589 * t59) * qJD(1)) * t382 + (t28 * t589 + t59 * t596 - t17 * t215 + t69 * t137 + (-t488 * t58 + t601 * t69) * qJD(1)) * t380 - t58 * (-qJD(1) * t252 + t510) - t59 * (qJD(1) * t256 + t597) - t69 * t551 - ((t69 * t256 + t58 * t588) * t382 + (t69 * t252 + t588 * t59) * t380) * qJD(4)) * m(6) + (t43 * t461 + t92 * t410 + t481 * t285 + (-t54 * t380 - t55 * t382 + (-t382 * t83 + t657) * qJD(1)) * t301 - (t253 * t82 - t664) * qJD(1) - (t92 * (-t253 * t380 - t257 * t382) + t481 * t492) * qJD(4)) * m(5) + (t2 + t740 * qJD(1) + ((t736 * qJD(1) + t726 * t382) * t382 + (t724 * t380 + t737 * qJD(1) + (-t723 + t725) * t382) * t380) * t722) * t682 + (t1 + t739 * qJD(1) + ((t738 * qJD(1) + t723 * t382) * t382 + (t725 * t380 + t767 * qJD(1) + (-t724 + t726) * t382) * t380) * t722) * t681 + (t11 + t744 + t746) * t575 / 0.2e1 + (t10 + t743 + t745) * t531; -m(6) * (t231 * t69 + t232 * t59 + t233 * t58) - m(7) * (t231 * t34 + t232 * t38 + t233 * t37) + 0.2e1 * ((t570 * t58 + t572 * t59 - t17) * t692 + (t37 * t570 + t38 * t572 - t7) * t691) * t381 + 0.2e1 * ((qJD(4) * t69 + t28 * t380 + t29 * t382 + t574 * t59 - t575 * t58) * t692 + (t12 * t380 + t13 * t382 - t37 * t575 + t554 + t635) * t691) * t379; t2 * t621 / 0.2e1 + (t379 * t70 + t381 * t487) * t689 + ((-qJD(4) * t487 + t16) * t379 + (-qJD(1) * t486 + qJD(4) * t70 + t380 * t6 + t382 * t5) * t381) * t687 + t1 * t624 / 0.2e1 + (t379 * t71 + t381 * t485) * t690 + ((-qJD(4) * t485 + t15) * t379 + (-qJD(1) * t484 + qJD(4) * t71 + t3 * t382 + t380 * t4) * t381) * t686 + t14 * t528 + t379 * (t655 + t656 + t674 - t678 + t679) / 0.2e1 + (t379 * t85 + t381 * t483) * t497 + ((-qJD(4) * t483 + t27) * t379 + (-qJD(1) * t482 + qJD(4) * t85 + t380 * t9 + t382 * t8) * t381) * t683 + (t411 * t264 - t265 * t412 + t382 * t423) * t688 + (t266 * t411 + t267 * t412 + t380 * t423) * t685 + (t427 * t379 + (t412 * t389 - t391 * t411) * t381) * t684 + (t379 * t530 + t381 * t531) * t11 + (-t547 / 0.2e1 + t379 * t527) * t10 + ((qJD(4) * t405 + t12 * t148 - t13 * t150 - t37 * t66 + t38 * t67) * t379 + (t37 * (-qJD(4) * t150 + t156 * t380) + t38 * (qJD(4) * t148 - t156 * t382) - t7 * t467 + t34 * (-t148 * t574 - t150 * t575 - t380 * t67 + t382 * t66) - (qJD(1) * t704 + t665 - t666) * t421) * t381 - t37 * (-t164 * t349 - t274 * t275) - t38 * (t163 * t349 - t273 * t275) - t34 * (t163 * t274 + t164 * t273)) * m(7);];
tauc  = t18(:);