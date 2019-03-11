% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPPRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:39
% EndTime: 2019-03-08 19:22:41
% DurationCPUTime: 62.73s
% Computational Cost: add. (54227->1481), mult. (151853->2074), div. (0->0), fcn. (181466->12), ass. (0->549)
t786 = Icges(3,4) - Icges(4,5);
t785 = Icges(4,4) + Icges(3,5);
t776 = Icges(3,6) - Icges(4,6);
t507 = cos(pkin(10));
t508 = cos(pkin(6));
t511 = sin(qJ(2));
t667 = t508 * t511;
t505 = sin(pkin(10));
t514 = cos(qJ(2));
t675 = t505 * t514;
t472 = t507 * t667 + t675;
t666 = t508 * t514;
t542 = -t505 * t511 + t507 * t666;
t693 = sin(pkin(11));
t694 = cos(pkin(11));
t376 = -t472 * t694 + t542 * t693;
t513 = cos(qJ(5));
t506 = sin(pkin(6));
t510 = sin(qJ(5));
t673 = t506 * t510;
t310 = -t376 * t513 + t507 * t673;
t451 = t542 * qJD(2);
t452 = t472 * qJD(2);
t328 = t451 * t694 + t452 * t693;
t210 = qJD(5) * t310 + t328 * t510;
t327 = t451 * t693 - t452 * t694;
t623 = qJD(5) * t327;
t156 = qJD(6) * t210 + t623;
t721 = t156 / 0.2e1;
t473 = t505 * t666 + t507 * t511;
t453 = t473 * qJD(2);
t612 = t505 * t667;
t627 = qJD(2) * t514;
t454 = -qJD(2) * t612 + t507 * t627;
t330 = -t453 * t694 + t454 * t693;
t669 = t507 * t514;
t474 = -t612 + t669;
t380 = -t473 * t693 - t474 * t694;
t613 = t505 * t673;
t617 = qJD(5) * t513;
t212 = -qJD(5) * t613 + t330 * t510 - t380 * t617;
t329 = -t453 * t693 - t454 * t694;
t622 = qJD(5) * t329;
t157 = qJD(6) * t212 + t622;
t720 = t157 / 0.2e1;
t671 = t506 * t513;
t311 = -t380 * t510 + t505 * t671;
t630 = qJD(2) * t506;
t497 = t505 * t630;
t528 = -t473 * t694 + t474 * t693;
t313 = qJD(5) * t528 + t497;
t207 = qJD(6) * t311 + t313;
t718 = t207 / 0.2e1;
t529 = t472 * t693 + t542 * t694;
t597 = t507 * t630;
t314 = qJD(5) * t529 - t597;
t543 = t376 * t510 + t507 * t671;
t208 = -qJD(6) * t543 + t314;
t716 = t208 / 0.2e1;
t449 = (t511 * t693 + t514 * t694) * t506;
t434 = qJD(2) * t449;
t670 = t506 * t514;
t672 = t506 * t511;
t450 = -t693 * t670 + t694 * t672;
t668 = t508 * t510;
t303 = -qJD(5) * t668 + t434 * t510 + t450 * t617;
t433 = t450 * qJD(2);
t619 = qJD(5) * t433;
t279 = qJD(6) * t303 - t619;
t715 = t279 / 0.2e1;
t420 = t450 * t510 + t508 * t513;
t503 = qJD(2) * t508;
t424 = qJD(5) * t449 + t503;
t301 = qJD(6) * t420 + t424;
t713 = t301 / 0.2e1;
t787 = Icges(3,1) + Icges(4,1);
t784 = Icges(3,2) + Icges(4,3);
t736 = (-t511 * t776 + t514 * t785) * t506;
t783 = qJD(2) * t736;
t509 = sin(qJ(6));
t512 = cos(qJ(6));
t220 = -t310 * t509 + t512 * t529;
t221 = t310 * t512 + t509 * t529;
t104 = Icges(7,5) * t221 + Icges(7,6) * t220 - Icges(7,3) * t543;
t685 = Icges(7,4) * t221;
t106 = Icges(7,2) * t220 - Icges(7,6) * t543 + t685;
t216 = Icges(7,4) * t220;
t108 = Icges(7,1) * t221 - Icges(7,5) * t543 + t216;
t213 = -qJD(5) * t311 + t330 * t513;
t312 = -t380 * t513 - t613;
t223 = t312 * t512 + t509 * t528;
t112 = -qJD(6) * t223 - t213 * t509 + t329 * t512;
t222 = -t312 * t509 + t512 * t528;
t113 = qJD(6) * t222 + t213 * t512 + t329 * t509;
t211 = qJD(5) * t543 + t328 * t513;
t110 = -qJD(6) * t221 - t211 * t509 + t327 * t512;
t111 = qJD(6) * t220 + t211 * t512 + t327 * t509;
t53 = Icges(7,5) * t111 + Icges(7,6) * t110 + Icges(7,3) * t210;
t55 = Icges(7,4) * t111 + Icges(7,2) * t110 + Icges(7,6) * t210;
t57 = Icges(7,1) * t111 + Icges(7,4) * t110 + Icges(7,5) * t210;
t10 = t104 * t212 + t106 * t112 + t108 * t113 + t222 * t55 + t223 * t57 + t311 * t53;
t105 = Icges(7,5) * t223 + Icges(7,6) * t222 + Icges(7,3) * t311;
t684 = Icges(7,4) * t223;
t107 = Icges(7,2) * t222 + Icges(7,6) * t311 + t684;
t217 = Icges(7,4) * t222;
t109 = Icges(7,1) * t223 + Icges(7,5) * t311 + t217;
t54 = Icges(7,5) * t113 + Icges(7,6) * t112 + Icges(7,3) * t212;
t56 = Icges(7,4) * t113 + Icges(7,2) * t112 + Icges(7,6) * t212;
t58 = Icges(7,1) * t113 + Icges(7,4) * t112 + Icges(7,5) * t212;
t11 = t105 * t212 + t107 * t112 + t109 * t113 + t222 * t56 + t223 * t58 + t311 * t54;
t421 = t450 * t513 - t668;
t305 = -t421 * t509 + t449 * t512;
t306 = t421 * t512 + t449 * t509;
t166 = Icges(7,5) * t306 + Icges(7,6) * t305 + Icges(7,3) * t420;
t683 = Icges(7,4) * t306;
t167 = Icges(7,2) * t305 + Icges(7,6) * t420 + t683;
t302 = Icges(7,4) * t305;
t168 = Icges(7,1) * t306 + Icges(7,5) * t420 + t302;
t304 = -qJD(5) * t420 + t434 * t513;
t170 = -qJD(6) * t306 - t304 * t509 - t433 * t512;
t171 = qJD(6) * t305 + t304 * t512 - t433 * t509;
t78 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t303;
t79 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t303;
t80 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t303;
t20 = t112 * t167 + t113 * t168 + t166 * t212 + t222 * t79 + t223 * t80 + t311 * t78;
t41 = t104 * t311 + t106 * t222 + t108 * t223;
t42 = t105 * t311 + t107 * t222 + t109 * t223;
t48 = t166 * t311 + t167 * t222 + t168 * t223;
t782 = t10 * t716 + t11 * t718 + t20 * t713 + t41 * t721 + t42 * t720 + t48 * t715;
t12 = t104 * t303 + t106 * t170 + t108 * t171 + t305 * t55 + t306 * t57 + t420 * t53;
t13 = t105 * t303 + t107 * t170 + t109 * t171 + t305 * t56 + t306 * t58 + t420 * t54;
t21 = t166 * t303 + t167 * t170 + t168 * t171 + t305 * t79 + t306 * t80 + t420 * t78;
t43 = t104 * t420 + t106 * t305 + t108 * t306;
t44 = t105 * t420 + t107 * t305 + t109 * t306;
t67 = t166 * t420 + t167 * t305 + t168 * t306;
t781 = t12 * t716 + t13 * t718 + t21 * t713 + t43 * t721 + t44 * t720 + t67 * t715;
t711 = t313 / 0.2e1;
t709 = t314 / 0.2e1;
t703 = t424 / 0.2e1;
t780 = qJD(5) / 0.2e1;
t674 = t506 * t507;
t765 = Icges(5,4) * t376;
t241 = -Icges(5,2) * t529 + Icges(5,6) * t674 - t765;
t775 = Icges(5,1) * t529 + t241 - t765;
t676 = t505 * t506;
t764 = Icges(5,4) * t380;
t242 = -Icges(5,2) * t528 - Icges(5,6) * t676 - t764;
t774 = Icges(5,1) * t528 + t242 - t764;
t438 = Icges(5,4) * t450;
t322 = -Icges(5,2) * t449 - Icges(5,6) * t508 + t438;
t773 = Icges(5,1) * t449 + t322 + t438;
t772 = t786 * t474;
t771 = t786 * t542;
t770 = t786 * t473;
t769 = t786 * t472;
t768 = -t773 * t508 + t674 * t775 - t676 * t774;
t100 = Icges(6,1) * t211 - Icges(6,4) * t210 + Icges(6,5) * t327;
t158 = Icges(6,5) * t310 + Icges(6,6) * t543 + Icges(6,3) * t529;
t688 = Icges(6,4) * t310;
t160 = Icges(6,2) * t543 + Icges(6,6) * t529 + t688;
t307 = Icges(6,4) * t543;
t162 = Icges(6,1) * t310 + Icges(6,5) * t529 + t307;
t96 = Icges(6,5) * t211 - Icges(6,6) * t210 + Icges(6,3) * t327;
t98 = Icges(6,4) * t211 - Icges(6,2) * t210 + Icges(6,6) * t327;
t29 = t100 * t421 - t158 * t433 - t160 * t303 + t162 * t304 - t420 * t98 + t449 * t96;
t101 = Icges(6,1) * t213 - Icges(6,4) * t212 + Icges(6,5) * t329;
t159 = Icges(6,5) * t312 - Icges(6,6) * t311 + Icges(6,3) * t528;
t687 = Icges(6,4) * t312;
t161 = -Icges(6,2) * t311 + Icges(6,6) * t528 + t687;
t308 = Icges(6,4) * t311;
t163 = Icges(6,1) * t312 + Icges(6,5) * t528 - t308;
t97 = Icges(6,5) * t213 - Icges(6,6) * t212 + Icges(6,3) * t329;
t99 = Icges(6,4) * t213 - Icges(6,2) * t212 + Icges(6,6) * t329;
t30 = t101 * t421 - t159 * t433 - t161 * t303 + t163 * t304 - t420 * t99 + t449 * t97;
t172 = Icges(6,5) * t304 - Icges(6,6) * t303 - Icges(6,3) * t433;
t173 = Icges(6,4) * t304 - Icges(6,2) * t303 - Icges(6,6) * t433;
t174 = Icges(6,1) * t304 - Icges(6,4) * t303 - Icges(6,5) * t433;
t276 = Icges(6,5) * t421 - Icges(6,6) * t420 + Icges(6,3) * t449;
t686 = Icges(6,4) * t421;
t277 = -Icges(6,2) * t420 + Icges(6,6) * t449 + t686;
t418 = Icges(6,4) * t420;
t278 = Icges(6,1) * t421 + Icges(6,5) * t449 - t418;
t45 = t172 * t449 - t173 * t420 + t174 * t421 - t276 * t433 - t277 * t303 + t278 * t304;
t70 = t158 * t449 - t160 * t420 + t162 * t421;
t71 = t159 * t449 - t161 * t420 + t163 * t421;
t84 = t276 * t449 - t277 * t420 + t278 * t421;
t536 = t327 * t70 + t329 * t71 - t433 * t84;
t767 = t29 * t709 + t30 * t711 + t45 * t703 + t536 * t780 + t781;
t26 = t100 * t312 + t158 * t329 - t160 * t212 + t162 * t213 - t311 * t98 + t528 * t96;
t27 = t101 * t312 + t159 * t329 - t161 * t212 + t163 * t213 - t311 * t99 + t528 * t97;
t38 = t172 * t528 - t173 * t311 + t174 * t312 - t212 * t277 + t213 * t278 + t276 * t329;
t63 = t158 * t528 - t160 * t311 + t162 * t312;
t64 = t159 * t528 - t161 * t311 + t163 * t312;
t76 = t276 * t528 - t277 * t311 + t278 * t312;
t537 = t327 * t63 + t329 * t64 - t433 * t76;
t766 = t26 * t709 + t27 * t711 + t38 * t703 + t537 * t780 + t782;
t753 = Icges(5,4) * t529;
t243 = -Icges(5,1) * t376 + Icges(5,5) * t674 - t753;
t763 = Icges(5,2) * t376 + t243 - t753;
t754 = Icges(5,4) * t528;
t244 = -Icges(5,1) * t380 - Icges(5,5) * t676 - t754;
t762 = Icges(5,2) * t380 + t244 - t754;
t689 = Icges(5,4) * t449;
t323 = Icges(5,1) * t450 - Icges(5,5) * t508 - t689;
t761 = -Icges(5,2) * t450 + t323 - t689;
t752 = -t542 * t784 + t674 * t776 - t769;
t751 = t473 * t784 - t676 * t776 - t772;
t750 = t472 * t787 - t674 * t785 + t771;
t749 = t474 * t787 + t676 * t785 - t770;
t748 = -t451 * t786 + t452 * t784;
t747 = t453 * t786 + t454 * t784;
t746 = t451 * t787 - t452 * t786;
t745 = -t453 * t787 - t454 * t786;
t500 = Icges(4,5) * t672;
t690 = Icges(3,4) * t511;
t742 = -Icges(4,3) * t670 - (Icges(3,2) * t514 + t690) * t506 + t500 - t776 * t508;
t501 = Icges(3,4) * t670;
t680 = Icges(4,5) * t514;
t741 = (Icges(4,1) * t511 - t680) * t506 + Icges(3,1) * t672 + t501 + t785 * t508;
t264 = Icges(5,5) * t529 - Icges(5,6) * t376;
t265 = Icges(5,5) * t528 - Icges(5,6) * t380;
t350 = Icges(5,5) * t449 + Icges(5,6) * t450;
t743 = -t473 * t785 - t474 * t776;
t744 = t472 * t776 - t542 * t785;
t760 = t508 * (t350 - t736) - t674 * (t264 + t744) - t676 * (-t265 + t743);
t759 = t761 * t508 - t674 * t763 + t676 * t762;
t19 = t110 * t167 + t111 * t168 + t166 * t210 + t220 * t79 + t221 * t80 - t543 * t78;
t39 = -t104 * t543 + t106 * t220 + t108 * t221;
t40 = -t105 * t543 + t107 * t220 + t109 * t221;
t47 = -t166 * t543 + t167 * t220 + t168 * t221;
t8 = t104 * t210 + t106 * t110 + t108 * t111 + t220 * t55 + t221 * t57 - t53 * t543;
t9 = t105 * t210 + t107 * t110 + t109 * t111 + t220 * t56 + t221 * t58 - t54 * t543;
t1 = t156 * t39 + t157 * t40 + t19 * t301 + t207 * t9 + t208 * t8 + t279 * t47;
t24 = t100 * t310 + t158 * t327 - t160 * t210 + t162 * t211 + t529 * t96 + t543 * t98;
t25 = t101 * t310 + t159 * t327 - t161 * t210 + t163 * t211 + t529 * t97 + t543 * t99;
t37 = t172 * t529 + t173 * t543 + t174 * t310 - t210 * t277 + t211 * t278 + t276 * t327;
t61 = t158 * t529 + t160 * t543 + t162 * t310;
t62 = t159 * t529 + t161 * t543 + t163 * t310;
t75 = t276 * t529 + t277 * t543 + t278 * t310;
t538 = t327 * t61 + t329 * t62 - t433 * t75;
t758 = qJD(5) * t538 + t314 * t24 + t25 * t313 + t37 * t424 + t1;
t755 = pkin(3) * t473;
t504 = t506 ^ 2;
t740 = t452 - (t473 * t508 + t504 * t675) * qJD(2);
t739 = t454 - (t504 * t669 + t508 * t542) * qJD(2);
t475 = (Icges(4,3) * t511 + t680) * t506;
t738 = qJD(2) * t475 - (Icges(3,4) * t514 - Icges(3,2) * t511) * t630;
t480 = (Icges(3,1) * t514 - t690) * t506;
t737 = (Icges(4,1) * t514 + Icges(4,5) * t511) * t630 + qJD(2) * t480;
t735 = -Icges(3,2) * t672 - t475 + t501 + t741;
t734 = Icges(4,1) * t670 + t480 + t500 + t742;
t733 = -t474 * t784 + t749 - t770;
t732 = -t472 * t784 + t750 + t771;
t731 = -t473 * t787 + t751 - t772;
t730 = -t542 * t787 - t752 + t769;
t729 = Icges(5,5) * t328 - Icges(5,6) * t327 - t451 * t785 + t452 * t776;
t728 = -Icges(5,5) * t330 + Icges(5,6) * t329 - t453 * t785 - t454 * t776;
t726 = -Icges(5,5) * t434 - Icges(5,6) * t433 + t783;
t114 = rSges(7,1) * t221 + rSges(7,2) * t220 - rSges(7,3) * t543;
t115 = rSges(7,1) * t223 + rSges(7,2) * t222 + rSges(7,3) * t311;
t204 = pkin(5) * t310 - pkin(9) * t543;
t206 = pkin(5) * t312 + pkin(9) * t311;
t273 = -pkin(4) * t376 + pkin(8) * t529;
t275 = -pkin(4) * t380 + pkin(8) * t528;
t425 = pkin(3) * t472 + qJ(4) * t674;
t426 = pkin(3) * t474 - qJ(4) * t676;
t402 = pkin(2) * t472 - qJ(3) * t542;
t406 = pkin(2) * t474 + qJ(3) * t473;
t593 = t402 * t497 + t406 * t597 + qJD(1);
t624 = qJD(4) * t508;
t530 = t425 * t497 + t426 * t597 + t593 - t624;
t626 = qJD(3) * t514;
t520 = t273 * t497 + t275 * t597 - t506 * t626 + t530;
t32 = t207 * t114 - t208 * t115 + t313 * t204 - t314 * t206 + t520;
t662 = t115 + t206;
t130 = pkin(5) * t213 + pkin(9) * t212;
t60 = rSges(7,1) * t113 + rSges(7,2) * t112 + rSges(7,3) * t212;
t696 = t130 + t60;
t129 = pkin(5) * t211 + pkin(9) * t210;
t239 = pkin(4) * t328 + pkin(8) * t327;
t240 = pkin(4) * t330 + pkin(8) * t329;
t625 = qJD(4) * t506;
t496 = t507 * t625;
t422 = pkin(3) * t451 + t496;
t596 = t505 * t625;
t423 = -pkin(3) * t453 - t596;
t461 = qJD(3) * t542;
t299 = pkin(2) * t451 + qJ(3) * t452 - t461;
t463 = qJD(3) * t473;
t300 = -pkin(2) * t453 + qJ(3) * t454 + t463;
t499 = qJD(3) * t672;
t605 = qJD(2) * t499 + t299 * t497 + t300 * t597;
t564 = t422 * t497 + t423 * t597 + t605;
t540 = t239 * t497 + t240 * t597 + t564;
t59 = rSges(7,1) * t111 + rSges(7,2) * t110 + rSges(7,3) * t210;
t7 = t114 * t157 - t115 * t156 + t129 * t313 - t130 * t314 + t207 * t59 - t208 * t60 + (t204 * t329 - t206 * t327) * qJD(5) + t540;
t723 = t32 * t696 + t662 * t7;
t722 = qJD(2) ^ 2;
t719 = -t207 / 0.2e1;
t717 = -t208 / 0.2e1;
t714 = -t301 / 0.2e1;
t712 = -t313 / 0.2e1;
t710 = -t314 / 0.2e1;
t704 = -t424 / 0.2e1;
t699 = pkin(3) * t542;
t697 = t129 + t59;
t190 = pkin(5) * t304 + pkin(9) * t303;
t81 = rSges(7,1) * t171 + rSges(7,2) * t170 + rSges(7,3) * t303;
t695 = t190 + t81;
t679 = t529 * t510;
t678 = t528 * t510;
t677 = t449 * t510;
t665 = t509 * t513;
t664 = t512 * t513;
t663 = t114 + t204;
t169 = rSges(7,1) * t306 + rSges(7,2) * t305 + rSges(7,3) * t420;
t298 = pkin(5) * t421 + pkin(9) * t420;
t661 = t169 + t298;
t660 = t299 * t676 + t300 * t674;
t659 = qJD(3) * t452 + t300 * t503;
t292 = t508 * t300;
t658 = t508 * t423 + t292;
t368 = rSges(4,1) * t451 + rSges(4,3) * t452;
t657 = -t299 - t368;
t656 = -t299 - t422;
t344 = rSges(4,1) * t472 - rSges(4,2) * t674 - rSges(4,3) * t542;
t647 = -t344 - t402;
t646 = t402 * t676 + t406 * t674;
t566 = pkin(2) * t514 + qJ(3) * t511;
t482 = t566 * t506;
t645 = -rSges(5,1) * t449 - rSges(5,2) * t450 - t482;
t403 = -pkin(2) * t473 + qJ(3) * t474;
t644 = qJD(3) * t472 + t403 * t503;
t643 = t406 * t503 - t461;
t384 = t508 * t406;
t642 = t508 * t426 + t384;
t641 = -t402 - t425;
t427 = (qJD(2) * t566 - t626) * t506;
t483 = (rSges(4,1) * t514 + rSges(4,3) * t511) * t506;
t640 = -qJD(2) * t483 - t427;
t639 = -pkin(3) * t506 * t627 - t427 + t624;
t481 = (pkin(2) * t511 - qJ(3) * t514) * t506;
t634 = -t508 * rSges(4,2) - (rSges(4,1) * t511 - rSges(4,3) * t514) * t506 - t481;
t633 = -pkin(3) * t672 + qJ(4) * t508 - t481;
t632 = -t482 - t483;
t631 = qJD(2) * t505;
t628 = qJD(2) * t511;
t621 = qJD(5) * t376;
t620 = qJD(5) * t380;
t618 = qJD(5) * t450;
t616 = qJD(6) * t510;
t614 = pkin(3) * t504 * t514;
t611 = t508 * t240 + t658;
t237 = rSges(5,1) * t328 - rSges(5,2) * t327;
t610 = -t237 + t656;
t609 = -t239 + t656;
t245 = -rSges(5,1) * t376 - rSges(5,2) * t529 + rSges(5,3) * t674;
t608 = -t245 + t641;
t607 = t508 * t275 + t642;
t606 = -t273 + t641;
t604 = t423 * t503 + t659;
t603 = -rSges(5,1) * t434 - rSges(5,2) * t433 + t639;
t602 = -pkin(4) * t434 + pkin(8) * t433 + t639;
t601 = -rSges(5,1) * t450 + rSges(5,2) * t449 + rSges(5,3) * t508 + t633;
t399 = pkin(2) * t542 + qJ(3) * t472;
t600 = t399 * t497 + t403 * t597 + t499;
t599 = -pkin(4) * t450 - pkin(8) * t449 + t633;
t598 = -t503 * t755 + t644;
t591 = t623 / 0.2e1;
t590 = t622 / 0.2e1;
t589 = -t619 / 0.2e1;
t588 = -t399 - t699;
t587 = t640 * t506;
t586 = t634 * t506;
t175 = rSges(6,1) * t304 - rSges(6,2) * t303 - rSges(6,3) * t433;
t585 = -t175 + t602;
t274 = pkin(4) * t528 + pkin(8) * t380;
t584 = t274 * t503 + t598;
t280 = rSges(6,1) * t421 - rSges(6,2) * t420 + rSges(6,3) * t449;
t583 = -t280 + t599;
t582 = t422 * t676 + t423 * t674 + t660;
t581 = t425 * t676 + t426 * t674 + t646;
t580 = t426 * t503 + t496 + t643;
t575 = t506 * t603;
t574 = t506 * t601;
t573 = t602 * t506;
t572 = t599 * t506;
t571 = pkin(5) * t513 + pkin(9) * t510;
t570 = t463 - t596;
t569 = t602 - t695;
t568 = rSges(6,1) * t513 - rSges(6,2) * t510;
t567 = -rSges(7,1) * t512 + rSges(7,2) * t509;
t565 = t599 - t661;
t563 = t497 * t699 - t597 * t755 + t600;
t562 = Icges(6,1) * t513 - Icges(6,4) * t510;
t561 = -Icges(7,1) * t512 + Icges(7,4) * t509;
t560 = Icges(6,4) * t513 - Icges(6,2) * t510;
t559 = -Icges(7,4) * t512 + Icges(7,2) * t509;
t558 = Icges(6,5) * t513 - Icges(6,6) * t510;
t557 = -Icges(7,5) * t512 + Icges(7,6) * t509;
t556 = -t160 * t510 + t162 * t513;
t555 = -t161 * t510 + t163 * t513;
t164 = rSges(6,1) * t310 + rSges(6,2) * t543 + rSges(6,3) * t529;
t165 = rSges(6,1) * t312 - rSges(6,2) * t311 + rSges(6,3) * t528;
t554 = t164 * t329 - t165 * t327;
t553 = t164 * t433 + t280 * t327;
t552 = -t165 * t433 - t280 * t329;
t551 = -t277 * t510 + t278 * t513;
t347 = rSges(3,1) * t474 - rSges(3,2) * t473 + rSges(3,3) * t676;
t444 = t508 * rSges(3,3) + (rSges(3,1) * t511 + rSges(3,2) * t514) * t506;
t544 = t347 * t508 - t444 * t676;
t281 = t544 * qJD(2);
t345 = rSges(3,1) * t472 + rSges(3,2) * t542 - rSges(3,3) * t674;
t545 = -t345 * t508 - t444 * t674;
t282 = t545 * qJD(2);
t550 = -t281 * t505 - t282 * t507;
t549 = t345 * t505 + t347 * t507;
t369 = rSges(3,1) * t451 - rSges(3,2) * t452;
t373 = -rSges(3,1) * t453 - rSges(3,2) * t454;
t548 = t369 * t505 + t373 * t507;
t547 = t239 * t676 + t240 * t674 + t582;
t546 = t273 * t676 + t275 * t674 + t581;
t541 = (rSges(3,1) * t514 - rSges(3,2) * t511) * t506;
t272 = pkin(4) * t529 + pkin(8) * t376;
t539 = t272 * t497 + t274 * t597 + t563;
t535 = -t614 + (-pkin(4) * t449 + pkin(8) * t450 - t482) * t506;
t534 = t32 * t697 + t663 * t7;
t533 = (Icges(7,5) * t220 - Icges(7,6) * t221) * t208 + (Icges(7,5) * t222 - Icges(7,6) * t223) * t207 + (Icges(7,5) * t305 - Icges(7,6) * t306) * t301;
t532 = (Icges(6,5) * t543 - Icges(6,6) * t310) * t314 + (-Icges(6,5) * t311 - Icges(6,6) * t312) * t313 + (-Icges(6,5) * t420 - Icges(6,6) * t421) * t424;
t526 = t240 * t503 + t573 * t631 + t604;
t16 = t115 * t279 + t130 * t424 - t157 * t169 - t190 * t313 - t207 * t81 + t301 * t60 + (-t206 * t433 - t298 * t329) * qJD(5) + t526;
t525 = t275 * t503 + t572 * t631 + t580;
t35 = t115 * t301 - t169 * t207 + t206 * t424 - t298 * t313 + t525;
t531 = t16 * t662 + t35 * t696;
t527 = t505 * t535;
t524 = (Icges(7,1) * t222 - t107 - t684) * t207 + (Icges(7,1) * t220 - t106 - t685) * t208 + (Icges(7,1) * t305 - t167 - t683) * t301;
t523 = (-Icges(7,2) * t223 + t109 + t217) * t207 + (-Icges(7,2) * t221 + t108 + t216) * t208 + (-Icges(7,2) * t306 + t168 + t302) * t301;
t522 = (-Icges(6,1) * t311 - t161 - t687) * t313 + (Icges(6,1) * t543 - t160 - t688) * t314 + (-Icges(6,1) * t420 - t277 - t686) * t424;
t521 = (Icges(6,2) * t312 - t163 + t308) * t313 + (Icges(6,2) * t310 - t162 - t307) * t314 + (Icges(6,2) * t421 - t278 + t418) * t424;
t437 = qJD(3) * t454;
t519 = t437 + (t507 * t573 + t508 * t609) * qJD(2);
t518 = (-t272 + t588) * t508 + t535 * t507;
t517 = (t507 * t572 + t508 * t606) * qJD(2) + t570;
t516 = t555 * t313 + t314 * t556 + t551 * t424;
t515 = (Icges(7,3) * t312 + t107 * t509 - t109 * t512 + t311 * t557) * t207 + (Icges(7,3) * t310 + t106 * t509 - t108 * t512 - t543 * t557) * t208 + (Icges(7,3) * t421 + t167 * t509 - t168 * t512 + t420 * t557) * t301;
t466 = qJD(2) * t541;
t464 = qJD(3) * t474;
t405 = -rSges(3,1) * t473 - rSges(3,2) * t474;
t404 = -rSges(4,1) * t473 + rSges(4,3) * t474;
t401 = rSges(3,1) * t542 - rSges(3,2) * t472;
t400 = rSges(4,1) * t542 + rSges(4,3) * t472;
t372 = -rSges(4,1) * t453 + rSges(4,3) * t454;
t346 = rSges(4,1) * t474 + rSges(4,2) * t676 + rSges(4,3) * t473;
t331 = t449 * t616 - t618;
t325 = t571 * t449;
t320 = t449 * t664 - t450 * t509;
t319 = -t449 * t665 - t450 * t512;
t318 = Icges(5,1) * t434 + Icges(5,4) * t433;
t317 = Icges(5,4) * t434 + Icges(5,2) * t433;
t315 = (t473 * t507 - t505 * t542) * t630;
t297 = -pkin(5) * t420 + pkin(9) * t421;
t296 = -rSges(6,1) * t420 - rSges(6,2) * t421;
t286 = -rSges(6,3) * t450 + t449 * t568;
t285 = -Icges(6,5) * t450 + t449 * t562;
t284 = -Icges(6,6) * t450 + t449 * t560;
t283 = -Icges(6,3) * t450 + t449 * t558;
t271 = rSges(5,1) * t528 - rSges(5,2) * t380;
t270 = rSges(5,1) * t529 - rSges(5,2) * t376;
t262 = t528 * t616 + t620;
t261 = t529 * t616 + t621;
t258 = t571 * t528;
t257 = t571 * t529;
t256 = t380 * t509 + t528 * t664;
t255 = t380 * t512 - t528 * t665;
t254 = t376 * t509 + t529 * t664;
t253 = t376 * t512 - t529 * t665;
t246 = -rSges(5,1) * t380 - rSges(5,2) * t528 - rSges(5,3) * t676;
t238 = rSges(5,1) * t330 - rSges(5,2) * t329;
t236 = Icges(5,1) * t330 - Icges(5,4) * t329;
t235 = Icges(5,1) * t328 - Icges(5,4) * t327;
t234 = Icges(5,4) * t330 - Icges(5,2) * t329;
t233 = Icges(5,4) * t328 - Icges(5,2) * t327;
t229 = rSges(7,3) * t421 + t420 * t567;
t228 = Icges(7,5) * t421 + t420 * t561;
t227 = Icges(7,6) * t421 + t420 * t559;
t224 = t548 * t630;
t209 = t549 * t630 + qJD(1);
t205 = -pkin(5) * t311 + pkin(9) * t312;
t203 = pkin(5) * t543 + pkin(9) * t310;
t202 = -rSges(6,1) * t311 - rSges(6,2) * t312;
t201 = rSges(6,1) * t543 - rSges(6,2) * t310;
t194 = rSges(7,1) * t305 - rSges(7,2) * t306;
t189 = rSges(7,1) * t320 + rSges(7,2) * t319 + rSges(7,3) * t677;
t188 = Icges(7,1) * t320 + Icges(7,4) * t319 + Icges(7,5) * t677;
t187 = Icges(7,4) * t320 + Icges(7,2) * t319 + Icges(7,6) * t677;
t186 = Icges(7,5) * t320 + Icges(7,6) * t319 + Icges(7,3) * t677;
t183 = rSges(6,3) * t380 + t528 * t568;
t182 = rSges(6,3) * t376 + t529 * t568;
t181 = Icges(6,5) * t380 + t528 * t562;
t180 = Icges(6,5) * t376 + t529 * t562;
t179 = Icges(6,6) * t380 + t528 * t560;
t178 = Icges(6,6) * t376 + t529 * t560;
t177 = Icges(6,3) * t380 + t528 * t558;
t176 = Icges(6,3) * t376 + t529 * t558;
t155 = t463 + (t507 * t586 + t508 * t647) * qJD(2);
t154 = (t346 * t508 + t505 * t586) * qJD(2) + t643;
t153 = rSges(7,3) * t312 + t311 * t567;
t152 = rSges(7,3) * t310 - t543 * t567;
t151 = Icges(7,5) * t312 + t311 * t561;
t150 = Icges(7,5) * t310 - t543 * t561;
t149 = Icges(7,6) * t312 + t311 * t559;
t148 = Icges(7,6) * t310 - t543 * t559;
t141 = t437 + (t507 * t587 + t657 * t508) * qJD(2);
t140 = (t372 * t508 + t505 * t587) * qJD(2) + t659;
t139 = (-t626 + (t344 * t505 + t346 * t507) * qJD(2)) * t506 + t593;
t138 = rSges(7,1) * t222 - rSges(7,2) * t223;
t137 = rSges(7,1) * t220 - rSges(7,2) * t221;
t128 = rSges(7,1) * t256 + rSges(7,2) * t255 + rSges(7,3) * t678;
t127 = rSges(7,1) * t254 + rSges(7,2) * t253 + rSges(7,3) * t679;
t126 = Icges(7,1) * t256 + Icges(7,4) * t255 + Icges(7,5) * t678;
t125 = Icges(7,1) * t254 + Icges(7,4) * t253 + Icges(7,5) * t679;
t124 = Icges(7,4) * t256 + Icges(7,2) * t255 + Icges(7,6) * t678;
t123 = Icges(7,4) * t254 + Icges(7,2) * t253 + Icges(7,6) * t679;
t122 = Icges(7,5) * t256 + Icges(7,6) * t255 + Icges(7,3) * t678;
t121 = Icges(7,5) * t254 + Icges(7,6) * t253 + Icges(7,3) * t679;
t120 = (t368 * t505 + t372 * t507) * t630 + t605;
t103 = rSges(6,1) * t213 - rSges(6,2) * t212 + rSges(6,3) * t329;
t102 = rSges(6,1) * t211 - rSges(6,2) * t210 + rSges(6,3) * t327;
t95 = (t507 * t574 + t508 * t608) * qJD(2) + t570;
t94 = (t246 * t508 + t505 * t574) * qJD(2) + t580;
t83 = t437 + (t507 * t575 + t610 * t508) * qJD(2);
t82 = (t238 * t508 + t505 * t575) * qJD(2) + t604;
t77 = (-t626 + (t245 * t505 + t246 * t507) * qJD(2)) * t506 + t530;
t72 = (t237 * t505 + t238 * t507) * t630 + t564;
t66 = -t164 * t424 + t280 * t314 + t517;
t65 = t165 * t424 - t280 * t313 + t525;
t46 = t313 * t164 - t314 * t165 + t520;
t36 = -t114 * t301 + t169 * t208 - t204 * t424 + t298 * t314 + t517;
t34 = qJD(5) * t553 - t102 * t424 + t175 * t314 + t519;
t33 = qJD(5) * t552 + t103 * t424 - t175 * t313 + t526;
t31 = qJD(5) * t554 + t102 * t313 - t103 * t314 + t540;
t28 = t313 * t71 + t314 * t70 + t424 * t84;
t23 = t313 * t64 + t314 * t63 + t424 * t76;
t22 = t313 * t62 + t314 * t61 + t424 * t75;
t18 = t207 * t44 + t208 * t43 + t301 * t67;
t17 = -t114 * t279 - t129 * t424 + t156 * t169 + t190 * t314 + t208 * t81 - t301 * t59 + (t204 * t433 + t298 * t327) * qJD(5) + t519;
t15 = t207 * t42 + t208 * t41 + t301 * t48;
t14 = t207 * t40 + t208 * t39 + t301 * t47;
t2 = [m(3) * t224 + m(4) * t120 + m(5) * t72 + m(6) * t31 + m(7) * t7; -(t508 * ((-t508 * t350 + t761 * t449 + t773 * t450) * t508 + ((t264 * t507 - t265 * t505) * t508 + (t505 * t774 - t507 * t775) * t450 + (t505 * t762 - t507 * t763) * t449) * t506) + ((t473 * t732 + t474 * t730) * t674 + (-t473 * t735 + t474 * t734) * t508 + (-t473 * t733 + t474 * t731 - t760) * t676 + t759 * t528 + t768 * t380) * t676) * t722 / 0.2e1 + (t16 * t607 + t35 * t611 + t7 * t546 + (t17 * (t606 - t663) + t36 * (t609 - t697) + t531) * t508 + ((t17 * t565 + t36 * t569 + t723) * t507 + (t16 * t565 + t35 * t569 + t534) * t505) * t506 - t36 * (-t331 * t114 - t301 * t127 + t261 * t169 + t208 * t189 - t424 * t257 + t314 * t325 + t464) - t35 * (t331 * t115 + t301 * t128 - t262 * t169 - t207 * t189 + t424 * t258 - t313 * t325 + t584) - (t36 * (t204 * t450 + t298 * t376) + t35 * (-t206 * t450 - t298 * t380)) * qJD(5) - (t35 * t527 + t36 * t518) * qJD(2) + (t547 - t114 * t262 + t115 * t261 - t127 * t207 + t128 * t208 - t257 * t313 + t258 * t314 - t539 - (t204 * t380 - t206 * t376) * qJD(5)) * t32) * m(7) + t676 * t766 + t508 * t767 + (-t95 * t464 - t94 * t598 - t77 * t563 - ((-t505 * t94 - t507 * t95) * t614 + (t95 * (-t270 + t588) + t94 * t271) * t508 + ((t77 * t271 + t645 * t95) * t507 + (t77 * t270 + t645 * t94) * t505) * t506) * qJD(2) + t82 * t642 + t94 * t658 + t72 * t581 + t77 * t582 + (t94 * t238 + t82 * t246 + t608 * t83 + t610 * t95) * t508 + ((t77 * t238 + t72 * t246 + t601 * t83 + t603 * t95) * t507 + (t77 * t237 + t72 * t245 + t601 * t82 + t603 * t94) * t505) * t506) * m(5) + (-t155 * t464 - t154 * t644 - t139 * t600 - ((t155 * (-t399 - t400) + t154 * t404) * t508 + ((t139 * t404 + t155 * t632) * t507 + (t139 * t400 + t154 * t632) * t505) * t506) * qJD(2) + t140 * t384 + t154 * t292 + t120 * t646 + t139 * t660 + (t140 * t346 + t141 * t647 + t154 * t372 + t155 * t657) * t508 + ((t120 * t346 + t139 * t372 + t141 * t634 + t155 * t640) * t507 + (t120 * t344 + t139 * t368 + t140 * t634 + t154 * t640) * t505) * t506) * m(4) + (t508 * t84 + (t505 * t71 - t507 * t70) * t506) * t589 + (t508 * t76 + (t505 * t64 - t507 * t63) * t506) * t590 + (t508 * t75 + (t505 * t62 - t507 * t61) * t506) * t591 + (-t66 * (-t424 * t182 + t314 * t286 + t464) - t65 * (t424 * t183 - t313 * t286 + t584) - t46 * (t182 * t313 - t183 * t314 + t539) - (t66 * (t164 * t450 + t280 * t376) + t65 * (-t165 * t450 - t280 * t380) + t46 * (t164 * t380 - t165 * t376)) * qJD(5) - (t518 * t66 + t527 * t65) * qJD(2) + t33 * t607 + t65 * t611 + t31 * t546 + t46 * t547 + (t34 * (-t164 + t606) + t66 * (-t102 + t609) + t33 * t165 + t65 * t103) * t508 + ((t46 * t103 + t31 * t165 + t34 * t583 + t585 * t66) * t507 + (t46 * t102 + t31 * t164 + t33 * t583 + t585 * t65) * t505) * t506) * m(6) - t261 * t14 / 0.2e1 - t262 * t15 / 0.2e1 + ((t105 * t679 + t107 * t253 + t109 * t254 - t122 * t543 + t124 * t220 + t126 * t221) * t207 + t40 * t262 + (t104 * t679 + t106 * t253 + t108 * t254 - t121 * t543 + t123 * t220 + t125 * t221) * t208 + t39 * t261 + (t166 * t679 + t167 * t253 + t168 * t254 - t186 * t543 + t187 * t220 + t188 * t221) * t301 + t47 * t331) * t717 + ((t159 * t376 + t177 * t529 + t179 * t543 + t181 * t310) * t313 + (t158 * t376 + t176 * t529 + t178 * t543 + t180 * t310) * t314 + (t276 * t376 + t283 * t529 + t284 * t543 + t285 * t310) * t424 + t516 * t529 + (t376 * t61 + t380 * t62 - t450 * t75) * qJD(5)) * t710 + ((-t317 * t528 - t318 * t380 - t322 * t329 + t323 * t330 - t453 * t741 + t454 * t742 + t473 * t738 + t474 * t737 + t676 * t726) * t508 + ((t233 * t528 + t235 * t380 + t241 * t329 - t243 * t330 + t453 * t750 - t454 * t752 - t473 * t748 - t474 * t746 + t676 * t729) * t507 + (-t234 * t528 - t236 * t380 - t242 * t329 + t244 * t330 - t453 * t749 + t454 * t751 + t473 * t747 + t474 * t745 + t676 * t728) * t505) * t506) * t497 + (t45 * t508 + (-t29 * t507 + t30 * t505) * t506) * t703 + ((-t159 * t450 - t179 * t420 + t181 * t421) * t313 + (-t158 * t450 - t178 * t420 + t180 * t421) * t314 + (-t276 * t450 - t284 * t420 + t285 * t421) * t424 + (t376 * t70 + t380 * t71 - t450 * t84) * qJD(5) + ((t177 + t555) * t313 + (t176 + t556) * t314 + (t283 + t551) * t424) * t449) * t704 + (t37 * t508 + (-t24 * t507 + t25 * t505) * t506) * t709 + (t38 * t508 + (-t26 * t507 + t27 * t505) * t506) * t711 + (t21 * t508 + (-t12 * t507 + t13 * t505) * t506) * t713 + ((t105 * t677 + t107 * t319 + t109 * t320 + t122 * t420 + t124 * t305 + t126 * t306) * t207 + t44 * t262 + (t104 * t677 + t106 * t319 + t108 * t320 + t121 * t420 + t123 * t305 + t125 * t306) * t208 + t43 * t261 + (t166 * t677 + t167 * t319 + t168 * t320 + t186 * t420 + t187 * t305 + t188 * t306) * t301 + t67 * t331) * t714 + (t508 * t67 + (-t43 * t507 + t44 * t505) * t506) * t715 + ((t281 * t373 - t282 * t369) * t508 + (t209 * t548 + t224 * t549 + t466 * t550) * t506 + (-(t281 * t405 - t282 * t401) * t508 - (t209 * (t401 * t505 + t405 * t507) + t550 * t541) * t506 + (-t369 * t508 - t466 * t674) * t545 + (t373 * t508 - t466 * t676) * t544) * qJD(2)) * m(3) + ((t317 * t529 + t318 * t376 + t322 * t327 - t323 * t328 - t451 * t741 - t452 * t742 - t472 * t737 + t542 * t738 + t674 * t726) * t508 + ((-t233 * t529 - t235 * t376 - t241 * t327 + t243 * t328 + t451 * t750 + t452 * t752 + t472 * t746 - t542 * t748 + t674 * t729) * t507 + (t234 * t529 + t236 * t376 + t242 * t327 - t244 * t328 - t451 * t749 - t452 * t751 - t472 * t745 + t542 * t747 + t674 * t728) * t505) * t506) * t597 + t28 * t618 / 0.2e1 - t23 * t620 / 0.2e1 + ((t472 * t731 + t542 * t733) * t676 + (t472 * t734 + t542 * t735) * t508 + (t472 * t730 - t542 * t732 + t760) * t674 + t759 * t529 + t768 * t376) * t722 * t674 / 0.2e1 - t758 * t674 / 0.2e1 + (t19 * t508 + (t505 * t9 - t507 * t8) * t506) * t716 + (t20 * t508 + (-t10 * t507 + t11 * t505) * t506) * t718 + ((t105 * t678 + t107 * t255 + t109 * t256 + t122 * t311 + t124 * t222 + t126 * t223) * t207 + t42 * t262 + (t104 * t678 + t106 * t255 + t108 * t256 + t121 * t311 + t123 * t222 + t125 * t223) * t208 + t41 * t261 + (t166 * t678 + t167 * t255 + t168 * t256 + t186 * t311 + t187 * t222 + t188 * t223) * t301 + t48 * t331) * t719 + (t48 * t508 + (-t41 * t507 + t42 * t505) * t506) * t720 + (t47 * t508 + (-t39 * t507 + t40 * t505) * t506) * t721 - t22 * t621 / 0.2e1 - t331 * t18 / 0.2e1 + ((t159 * t380 + t177 * t528 - t179 * t311 + t181 * t312) * t313 + (t158 * t380 + t176 * t528 - t178 * t311 + t180 * t312) * t314 + (t276 * t380 + t283 * t528 - t284 * t311 + t285 * t312) * t424 + t516 * t528 + (t376 * t63 + t380 * t64 - t450 * t76) * qJD(5)) * t712 - (((t743 * t505 + t744 * t507 + t734 * t511 + t735 * t514) * t508 + ((t733 * t505 - t732 * t507) * t514 + (t731 * t505 + t730 * t507) * t511) * t506) * t630 + t508 ^ 2 * t783) * t503 / 0.2e1 + ((-t317 * t449 + t318 * t450 + t322 * t433 + t323 * t434 + t726 * t508) * t508 + ((t233 * t449 - t235 * t450 - t241 * t433 - t243 * t434 + (t748 * t514 - t746 * t511 + (-t511 * t752 - t514 * t750) * qJD(2)) * t506) * t507 + (-t234 * t449 + t236 * t450 + t242 * t433 + t244 * t434 + (-t747 * t514 + t745 * t511 + (t511 * t751 + t514 * t749) * qJD(2)) * t506) * t505 + (-t738 * t514 + t737 * t511 + (t511 * t742 + t514 * t741) * qJD(2) + t729 * t507 + t728 * t505) * t508) * t506) * t503; (-t16 * t542 + t17 * t473 + (t32 * t628 - t514 * t7) * t506 - t315 * t32 + t739 * t36 + t740 * t35) * m(7) + (-t33 * t542 + t34 * t473 + (-t31 * t514 + t46 * t628) * t506 - t315 * t46 + t739 * t66 + t740 * t65) * m(6) + (-t82 * t542 + t83 * t473 + (-t514 * t72 + t628 * t77) * t506 - t315 * t77 + t739 * t95 + t740 * t94) * m(5) + (-t140 * t542 + t141 * t473 + (-t120 * t514 + t139 * t628) * t506 - t139 * t315 + t739 * t155 + t740 * t154) * m(4); m(5) * (-t508 * t72 + (-t505 * t83 + t507 * t82) * t506) + m(6) * (-t31 * t508 + (t33 * t507 - t34 * t505) * t506) + m(7) * (-t508 * t7 + (t16 * t507 - t17 * t505) * t506); -(t14 * t310 + t15 * t312 + t18 * t421) * qJD(6) / 0.2e1 - (t28 + t18) * t433 / 0.2e1 + (t23 + t15) * t329 / 0.2e1 + (t22 + t14) * t327 / 0.2e1 + t528 * t766 + t449 * t767 + (t310 * t522 - t521 * t543 + t529 * t532) * t710 + ((t105 * t310 + t149 * t220 + t151 * t221) * t207 + (t104 * t310 + t148 * t220 + t150 * t221) * t208 + (t166 * t310 + t220 * t227 + t221 * t228) * t301 + (t310 * t39 + t312 * t40 + t421 * t47) * qJD(6) - t515 * t543) * t717 + (t420 * t521 + t421 * t522 + t449 * t532) * t704 + ((t105 * t421 + t149 * t305 + t151 * t306) * t207 + (t104 * t421 + t148 * t305 + t150 * t306) * t208 + (t166 * t421 + t227 * t305 + t228 * t306) * t301 + (t310 * t43 + t312 * t44 + t421 * t67) * qJD(6) + t515 * t420) * t714 + (t12 * t529 + t13 * t528 + t21 * t449 + t327 * t43 + t329 * t44 - t433 * t67) * t713 + (t19 * t449 + t327 * t39 + t329 * t40 - t433 * t47 + t528 * t9 + t529 * t8) * t716 + (t10 * t529 + t11 * t528 + t20 * t449 + t327 * t41 + t329 * t42 - t433 * t48) * t718 + (-(t35 * t662 - t36 * t663) * t433 + (t32 * t663 - t35 * t661) * t329 + (-t32 * t662 + t36 * t661) * t327 + (-t17 * t663 - t36 * t697 + t531) * t449 + (-t16 * t661 - t35 * t695 + t534) * t528 + (t17 * t661 + t36 * t695 - t723) * t529 - t36 * (-t152 * t301 - t203 * t424 + t208 * t229 + t297 * t314) - t35 * (t153 * t301 + t205 * t424 - t207 * t229 - t297 * t313) - t32 * (t152 * t207 - t153 * t208 + t203 * t313 - t205 * t314) - (t36 * (-t114 * t421 + t169 * t310) + t35 * (t115 * t421 - t169 * t312) + t32 * (t114 * t312 - t115 * t310)) * qJD(6)) * m(7) + t758 * t529 / 0.2e1 + ((t105 * t312 + t149 * t222 + t151 * t223) * t207 + (t104 * t312 + t148 * t222 + t150 * t223) * t208 + (t166 * t312 + t222 * t227 + t223 * t228) * t301 + (t310 * t41 + t312 * t42 + t421 * t48) * qJD(6) + t515 * t311) * t719 + (t449 * t84 + t528 * t71 + t529 * t70) * t589 + (t449 * t76 + t528 * t64 + t529 * t63) * t590 + (t449 * t75 + t528 * t62 + t529 * t61) * t591 + (t29 * t529 + t30 * t528 + t449 * t45 + t536) * t703 + (t24 * t529 + t25 * t528 + t37 * t449 + t538) * t709 + (t26 * t529 + t27 * t528 + t38 * t449 + t537) * t711 + (t311 * t521 + t312 * t522 + t528 * t532) * t712 + (t43 * t529 + t44 * t528 + t449 * t67) * t715 + (t34 * (-t164 * t449 + t280 * t529) + t33 * (t165 * t449 - t280 * t528) + t31 * (t164 * t528 - t165 * t529) + (-t102 * t449 + t175 * t529 + t201 * t424 - t296 * t314 + t553) * t66 + (t103 * t449 - t175 * t528 - t202 * t424 + t296 * t313 + t552) * t65 + (t102 * t528 - t103 * t529 - t201 * t313 + t202 * t314 + t554) * t46) * m(6) + (t41 * t529 + t42 * t528 + t449 * t48) * t720 + (t39 * t529 + t40 * t528 + t449 * t47) * t721; t212 * t15 / 0.2e1 + t311 * t782 + (t311 * t42 - t41 * t543 + t420 * t48) * t720 + (-t10 * t543 + t11 * t311 + t20 * t420 + t210 * t41 + t212 * t42 + t303 * t48) * t718 + t210 * t14 / 0.2e1 - t543 * t1 / 0.2e1 + (t311 * t40 - t39 * t543 + t420 * t47) * t721 + (t19 * t420 + t210 * t39 + t212 * t40 + t303 * t47 + t311 * t9 - t543 * t8) * t716 + t303 * t18 / 0.2e1 + t420 * t781 + (t311 * t44 + t420 * t67 - t43 * t543) * t715 + (-t12 * t543 + t13 * t311 + t21 * t420 + t210 * t43 + t212 * t44 + t303 * t67) * t713 + (t222 * t523 + t223 * t524 + t311 * t533) * t719 + (t220 * t523 + t221 * t524 - t533 * t543) * t717 + (t305 * t523 + t306 * t524 + t420 * t533) * t714 + (t17 * (-t114 * t420 - t169 * t543) + t16 * (t115 * t420 - t169 * t311) + t7 * (t114 * t311 + t115 * t543) + (-t114 * t303 + t137 * t301 + t169 * t210 - t194 * t208 - t420 * t59 - t543 * t81) * t36 + (t115 * t303 - t138 * t301 - t169 * t212 + t194 * t207 - t311 * t81 + t420 * t60) * t35 + (t114 * t212 - t115 * t210 - t137 * t207 + t138 * t208 + t311 * t59 + t543 * t60) * t32) * m(7);];
tauc  = t2(:);
