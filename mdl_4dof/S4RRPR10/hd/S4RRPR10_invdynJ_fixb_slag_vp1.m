% Calculate vector of inverse dynamics joint torques for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR10_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR10_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:54
% DurationCPUTime: 40.04s
% Computational Cost: add. (8738->891), mult. (22838->1157), div. (0->0), fcn. (20819->6), ass. (0->422)
t754 = Icges(4,4) - Icges(3,5);
t753 = Icges(4,5) - Icges(3,6);
t752 = Icges(4,1) + Icges(3,3);
t377 = sin(qJ(2));
t380 = cos(qJ(2));
t725 = t753 * t377 - t754 * t380;
t611 = Icges(3,4) * t377;
t289 = Icges(3,2) * t380 + t611;
t600 = Icges(4,6) * t377;
t443 = Icges(4,3) * t380 + t600;
t742 = -t289 - t443;
t363 = Icges(3,4) * t380;
t291 = Icges(3,1) * t377 + t363;
t599 = Icges(4,6) * t380;
t445 = Icges(4,2) * t377 + t599;
t751 = t291 + t445;
t381 = cos(qJ(1));
t750 = t752 * t381;
t378 = sin(qJ(1));
t583 = t378 * t380;
t586 = t377 * t378;
t728 = t754 * t583 - t753 * t586 + t750;
t736 = t752 * t378 + t725 * t381;
t601 = Icges(3,6) * t381;
t187 = Icges(3,4) * t583 - Icges(3,2) * t586 - t601;
t329 = Icges(4,6) * t586;
t609 = Icges(4,4) * t381;
t196 = Icges(4,2) * t583 - t329 + t609;
t749 = t187 * t377 - t196 * t380;
t292 = Icges(3,1) * t380 - t611;
t192 = Icges(3,5) * t378 + t292 * t381;
t444 = -Icges(4,3) * t377 + t599;
t193 = Icges(4,5) * t378 - t381 * t444;
t748 = -t192 * t583 - t193 * t586;
t604 = Icges(4,5) * t381;
t194 = Icges(4,6) * t583 - Icges(4,3) * t586 + t604;
t747 = t187 + t194;
t451 = -Icges(3,2) * t377 + t363;
t188 = Icges(3,6) * t378 + t381 * t451;
t746 = t188 - t193;
t334 = Icges(3,4) * t586;
t605 = Icges(3,5) * t381;
t191 = Icges(3,1) * t583 - t334 - t605;
t745 = t191 + t196;
t585 = t377 * t381;
t330 = Icges(4,6) * t585;
t581 = t380 * t381;
t610 = Icges(4,4) * t378;
t195 = -Icges(4,2) * t581 + t330 + t610;
t744 = t192 - t195;
t285 = Icges(3,5) * t377 + Icges(3,6) * t380;
t449 = Icges(4,4) * t377 + Icges(4,5) * t380;
t743 = t285 - t449;
t446 = Icges(4,2) * t380 - t600;
t740 = t292 + t446;
t739 = t742 * qJD(2);
t738 = t751 * qJD(2);
t721 = -t191 * t380 + t194 * t377 + t749;
t737 = t444 + t451;
t428 = t289 * t377 - t291 * t380;
t723 = -t377 * t443 + t380 * t445 - t428;
t735 = t736 * t381 + t748;
t674 = t192 * t581 + t193 * t585 + t736 * t378;
t734 = -t191 * t581 + t194 * t585 + t728 * t378;
t733 = t188 * t377 + t195 * t380;
t694 = -t721 * t378 + t728 * t381;
t693 = -t188 * t586 - t195 * t583 - t735;
t692 = -t187 * t585 + t196 * t581 - t734;
t691 = -t188 * t585 - t195 * t581 + t674;
t690 = t745 * t377 + t747 * t380;
t689 = t744 * t377 + t746 * t380;
t732 = t739 * t381 + (-t737 * t378 + t601 - t604) * qJD(1);
t731 = t746 * qJD(1) + t739 * t378;
t730 = -t738 * t381 + (-t740 * t378 + t605 - t609) * qJD(1);
t729 = t738 * t378 + (-t381 * t446 - t192 + t610) * qJD(1);
t727 = t737 * qJD(2);
t726 = t740 * qJD(2);
t724 = t743 * qJD(2);
t722 = -t377 * t751 + t742 * t380;
t720 = t192 * t380 + t193 * t377 - t733;
t672 = t743 * t378;
t587 = t285 * t381;
t104 = -t378 * t428 - t587;
t243 = t449 * t381;
t107 = t443 * t586 - t445 * t583 - t243;
t719 = t104 - t107;
t684 = t723 * t381 + t672;
t718 = t736 * qJD(1);
t376 = sin(qJ(4));
t379 = cos(qJ(4));
t584 = t378 * t379;
t232 = t376 * t381 + t377 * t584;
t582 = t379 * t381;
t233 = -t376 * t586 + t582;
t558 = t233 * rSges(5,1) - t232 * rSges(5,2);
t124 = rSges(5,3) * t583 - t558;
t472 = rSges(5,1) * t376 + rSges(5,2) * t379;
t676 = t380 * t472;
t208 = rSges(5,3) * t377 - t676;
t528 = qJD(4) * t380;
t532 = qJD(2) * t381;
t263 = -t378 * t528 + t532;
t530 = qJD(3) * t381;
t323 = t377 * t530;
t529 = qJD(4) * t377;
t346 = qJD(1) + t529;
t297 = pkin(2) * t377 - qJ(3) * t380;
t490 = -pkin(6) * t377 - t297;
t717 = -t124 * t346 - t263 * t208 + t490 * t532 + t323;
t716 = t381 ^ 2;
t715 = t743 * qJD(1) + t722 * qJD(2) - t727 * t377 + t726 * t380;
t714 = -qJD(2) * t689 - t377 * t732 + t380 * t730 + t718;
t713 = qJD(1) * t728 + qJD(2) * t690 + t377 * t731 + t380 * t729;
t712 = -t746 * t378 + t747 * t381;
t711 = t691 * t378 - t381 * t692;
t710 = t693 * t378 - t694 * t381;
t709 = t737 + t751;
t708 = t740 + t742;
t707 = (t329 + t334 + (Icges(3,2) + Icges(4,3)) * t583 - t745) * t381 + (-Icges(4,3) * t581 - t289 * t381 - t330 + t744) * t378;
t706 = t723 * qJD(1) - qJD(2) * t725;
t705 = qJD(1) * t721 - t378 * t724 + t718;
t704 = -t724 * t381 + (-t378 * t725 - t720 + t750) * qJD(1);
t703 = t684 * qJD(1);
t367 = t380 * rSges(5,3);
t206 = t377 * t472 + t367;
t372 = t381 * pkin(5);
t702 = t372 + t558;
t701 = t719 * qJD(1);
t534 = qJD(2) * t378;
t262 = t381 * t528 + t534;
t230 = -t376 * t378 + t377 * t582;
t231 = t376 * t585 + t584;
t111 = Icges(5,5) * t231 + Icges(5,6) * t230 + Icges(5,3) * t581;
t608 = Icges(5,4) * t231;
t114 = Icges(5,2) * t230 + Icges(5,6) * t581 + t608;
t215 = Icges(5,4) * t230;
t117 = Icges(5,1) * t231 + Icges(5,5) * t581 + t215;
t35 = t111 * t581 + t230 * t114 + t231 * t117;
t113 = -Icges(5,5) * t233 + Icges(5,6) * t232 + Icges(5,3) * t583;
t217 = Icges(5,4) * t233;
t116 = Icges(5,2) * t232 + Icges(5,6) * t583 - t217;
t216 = Icges(5,4) * t232;
t118 = Icges(5,1) * t233 - Icges(5,5) * t583 - t216;
t36 = t113 * t581 + t230 * t116 - t118 * t231;
t447 = Icges(5,5) * t376 + Icges(5,6) * t379;
t400 = -Icges(5,3) * t377 + t380 * t447;
t607 = Icges(5,4) * t376;
t448 = Icges(5,2) * t379 + t607;
t401 = -Icges(5,6) * t377 + t380 * t448;
t606 = Icges(5,4) * t379;
t452 = Icges(5,1) * t376 + t606;
t402 = -Icges(5,5) * t377 + t380 * t452;
t62 = -t230 * t401 - t231 * t402 - t400 * t581;
t12 = t262 * t35 - t263 * t36 + t62 * t346;
t37 = t111 * t583 + t232 * t114 - t233 * t117;
t38 = t113 * t583 + t116 * t232 + t118 * t233;
t63 = -t232 * t401 + t233 * t402 - t400 * t583;
t13 = t262 * t37 - t263 * t38 + t346 * t63;
t441 = t116 * t379 - t118 * t376;
t42 = t113 * t377 - t380 * t441;
t700 = t710 * qJD(2) + t701;
t699 = t711 * qJD(2) + t703;
t698 = qJD(2) * t721 + t377 * t729 - t380 * t731;
t697 = qJD(2) * t720 + t377 * t730 + t380 * t732;
t696 = -t706 * t378 + t715 * t381;
t695 = t715 * t378 + t706 * t381;
t683 = -t707 * t377 + t712 * t380;
t682 = (-t709 * t377 + t708 * t380) * qJD(1);
t681 = t705 * t716 + (t714 * t378 + (-t704 + t713) * t381) * t378;
t680 = t713 * t716 + (t704 * t378 + (-t705 + t714) * t381) * t378;
t679 = t728 + t733;
t524 = qJD(2) * qJD(3);
t678 = qJDD(3) * t377 + t380 * t524;
t505 = t377 * t532;
t537 = qJD(1) * t378;
t509 = t380 * t537;
t677 = t505 + t509;
t675 = t725 * qJD(1);
t673 = t587 - t243;
t366 = t378 * rSges(4,1);
t210 = -rSges(4,2) * t581 + rSges(4,3) * t585 + t366;
t257 = pkin(2) * t581 + qJ(3) * t585;
t306 = t381 * pkin(1) + t378 * pkin(5);
t478 = t257 + t306;
t126 = t210 + t478;
t359 = t377 * qJ(3);
t668 = t380 * pkin(2) + t359;
t252 = t668 * t378;
t305 = pkin(1) * t378 - t372;
t280 = qJD(1) * t305;
t670 = -qJD(1) * t252 - t280;
t669 = -pkin(6) * t380 - t668;
t364 = t377 * rSges(4,3);
t471 = -rSges(4,2) * t380 + t364;
t275 = t378 * pkin(3) + pkin(6) * t581;
t491 = -pkin(1) - t359;
t640 = -rSges(5,3) - pkin(2);
t521 = -pkin(6) + t640;
t667 = t521 * t380 + t491;
t347 = pkin(6) * t583;
t276 = pkin(3) * t381 - t347;
t556 = -t252 - t305;
t512 = t276 + t556;
t39 = qJD(1) * t512 + t717;
t122 = t231 * rSges(5,1) + t230 * rSges(5,2) + rSges(5,3) * t581;
t258 = t297 * t534;
t506 = t377 * t534;
t316 = pkin(6) * t506;
t358 = qJD(3) * t377;
t503 = t378 * t358;
t40 = t503 + t122 * t346 - t208 * t262 - t258 - t316 + (t275 + t478) * qJD(1);
t666 = t378 * t40 + t381 * t39;
t664 = (g(1) * t381 + g(2) * t378) * t377;
t241 = (Icges(5,2) * t376 - t606) * t380;
t395 = t262 * (-Icges(5,2) * t231 + t117 + t215) - t263 * (Icges(5,2) * t233 - t118 + t216) + t346 * (-t402 + t241);
t246 = (-Icges(5,1) * t379 + t607) * t380;
t396 = t262 * (-Icges(5,1) * t230 + t114 + t608) - t263 * (-Icges(5,1) * t232 + t116 - t217) + t346 * (-t401 - t246);
t655 = m(4) / 0.2e1;
t654 = m(5) / 0.2e1;
t525 = qJD(1) * qJD(2);
t277 = qJDD(2) * t378 + t381 * t525;
t522 = qJDD(4) * t380;
t147 = -qJD(4) * t677 + t381 * t522 + t277;
t653 = t147 / 0.2e1;
t278 = -qJDD(2) * t381 + t378 * t525;
t536 = qJD(1) * t381;
t508 = t380 * t536;
t410 = -t506 + t508;
t148 = qJD(4) * t410 + t378 * t522 + t278;
t652 = t148 / 0.2e1;
t261 = qJD(2) * t528 + qJDD(4) * t377 + qJDD(1);
t651 = t261 / 0.2e1;
t650 = -t262 / 0.2e1;
t649 = t262 / 0.2e1;
t648 = -t263 / 0.2e1;
t647 = t263 / 0.2e1;
t646 = t277 / 0.2e1;
t645 = t278 / 0.2e1;
t644 = -t346 / 0.2e1;
t643 = t346 / 0.2e1;
t442 = t114 * t379 + t117 * t376;
t480 = qJD(1) * t377 + qJD(4);
t501 = t380 * t532;
t403 = -t378 * t480 + t501;
t426 = t346 * t376;
t102 = t379 * t403 - t381 * t426;
t427 = t379 * t346;
t103 = t376 * t403 + t381 * t427;
t54 = Icges(5,5) * t103 + Icges(5,6) * t102 - Icges(5,3) * t677;
t56 = Icges(5,4) * t103 + Icges(5,2) * t102 - Icges(5,6) * t677;
t58 = Icges(5,1) * t103 + Icges(5,4) * t102 - Icges(5,5) * t677;
t8 = (qJD(2) * t442 + t54) * t377 + (qJD(2) * t111 - t376 * t58 - t379 * t56 + (t114 * t376 - t117 * t379) * qJD(4)) * t380;
t636 = t8 * t262;
t425 = t480 * t381;
t533 = qJD(2) * t380;
t100 = t379 * t425 + (t379 * t533 - t426) * t378;
t502 = t378 * t533;
t101 = t378 * t427 + (t425 + t502) * t376;
t53 = Icges(5,5) * t101 + Icges(5,6) * t100 + Icges(5,3) * t410;
t55 = Icges(5,4) * t101 + Icges(5,2) * t100 + Icges(5,6) * t410;
t57 = Icges(5,1) * t101 + Icges(5,4) * t100 + Icges(5,5) * t410;
t9 = (qJD(2) * t441 + t53) * t377 + (qJD(2) * t113 - t376 * t57 - t379 * t55 + (t116 * t376 + t118 * t379) * qJD(4)) * t380;
t635 = t9 * t263;
t181 = Icges(5,3) * t380 + t377 * t447;
t238 = (-Icges(5,5) * t379 + Icges(5,6) * t376) * t380;
t127 = qJD(2) * t181 + qJD(4) * t238;
t185 = Icges(5,6) * t380 + t377 * t448;
t130 = qJD(2) * t185 + qJD(4) * t241;
t189 = Icges(5,5) * t380 + t377 * t452;
t133 = qJD(2) * t189 + qJD(4) * t246;
t438 = -t376 * t402 - t379 * t401;
t19 = (qJD(2) * t438 + t127) * t377 + (-qJD(2) * t400 - t130 * t379 - t133 * t376 + (-t376 * t401 + t379 * t402) * qJD(4)) * t380;
t71 = -t377 * t400 - t380 * t438;
t632 = t19 * t346 + t71 * t261;
t631 = t103 * rSges(5,1) + t102 * rSges(5,2);
t630 = rSges(3,1) * t380;
t628 = rSges(4,2) * t377;
t299 = rSges(3,1) * t377 + rSges(3,2) * t380;
t256 = t299 * t381;
t365 = t378 * rSges(3,3);
t209 = rSges(3,1) * t581 - rSges(3,2) * t585 + t365;
t158 = t209 + t306;
t91 = qJD(1) * t158 - t299 * t534;
t625 = t256 * t91;
t507 = t299 * t532;
t543 = rSges(3,2) * t586 + t381 * rSges(3,3);
t207 = rSges(3,1) * t583 - t543;
t562 = -t207 - t305;
t90 = qJD(1) * t562 - t507;
t620 = t378 * t90;
t618 = t381 * t90;
t41 = t111 * t377 - t380 * t442;
t617 = t41 * t147;
t616 = t42 * t148;
t614 = -rSges(4,3) - qJ(3);
t531 = qJD(3) * t380;
t475 = t252 * t534 + t257 * t532 - t531;
t34 = t122 * t263 + t124 * t262 + (t275 * t381 - t276 * t378) * qJD(2) + t475;
t595 = qJD(2) * t34;
t580 = t380 * qJD(2) ^ 2;
t575 = t122 + t275;
t574 = t124 - t276;
t561 = -t210 - t257;
t560 = t378 * t252 + t381 * t257;
t218 = qJD(2) * t668 - t531;
t559 = -t471 * qJD(2) - t218;
t327 = qJ(3) * t581;
t254 = -pkin(2) * t585 + t327;
t557 = qJD(1) * t254 + t378 * t531;
t555 = -t257 - t275;
t355 = pkin(5) * t536;
t554 = qJD(1) * (-pkin(1) * t537 + t355) + qJDD(1) * t306;
t549 = t676 * t378;
t548 = t676 * t381;
t470 = rSges(4,3) * t380 + t628;
t547 = -t297 + t470;
t546 = -t668 - t471;
t545 = qJ(3) * t501 + t323;
t510 = t377 * t537;
t544 = rSges(3,2) * t510 + rSges(3,3) * t536;
t250 = rSges(4,2) * t586 + rSges(4,3) * t583;
t255 = rSges(4,2) * t585 + rSges(4,3) * t581;
t542 = t378 ^ 2 + t716;
t535 = qJD(2) * t377;
t517 = t40 * t536;
t120 = -pkin(2) * t677 - qJ(3) * t510 + t545;
t228 = t377 * t536 + t502;
t317 = pkin(2) * t506;
t121 = pkin(2) * t508 + qJ(3) * t228 - t317 + t503;
t516 = t381 * t120 + t378 * t121 + t252 * t536;
t515 = t278 * t297 + t381 * t678;
t325 = qJ(3) * t583;
t249 = -pkin(2) * t586 + t325;
t514 = t249 * t534 + t254 * t532 + t358;
t488 = rSges(4,1) * t381 - rSges(4,3) * t586;
t211 = rSges(4,2) * t583 + t488;
t513 = t211 + t556;
t511 = t355 + t545;
t500 = -pkin(1) - t630;
t497 = t536 / 0.2e1;
t496 = -t534 / 0.2e1;
t495 = t534 / 0.2e1;
t494 = -t532 / 0.2e1;
t493 = t532 / 0.2e1;
t483 = qJD(2) * t559;
t482 = -qJD(1) * t249 + t380 * t530;
t479 = rSges(4,1) * t536 + rSges(4,2) * t677 + rSges(4,3) * t501;
t477 = -t208 + t490;
t476 = -pkin(1) + (rSges(4,2) - pkin(2)) * t380;
t304 = rSges(2,1) * t381 - rSges(2,2) * t378;
t300 = rSges(2,1) * t378 + rSges(2,2) * t381;
t303 = -rSges(3,2) * t377 + t630;
t474 = rSges(5,1) * t101 + rSges(5,2) * t100;
t464 = t35 * t381 + t36 * t378;
t463 = t35 * t378 - t36 * t381;
t462 = t37 * t381 + t378 * t38;
t461 = t37 * t378 - t38 * t381;
t460 = t378 * t42 + t381 * t41;
t459 = t378 * t41 - t381 * t42;
t454 = -t378 * t91 - t618;
t440 = t122 * t378 - t124 * t381;
t143 = -rSges(3,1) * t677 - rSges(3,2) * t501 + t544;
t251 = t299 * t378;
t144 = -qJD(2) * t251 + (t303 * t381 + t365) * qJD(1);
t439 = t143 * t381 + t144 * t378;
t432 = t207 * t378 + t209 * t381;
t253 = (-rSges(5,1) * t379 + rSges(5,2) * t376) * t380;
t142 = qJD(2) * t206 + qJD(4) * t253;
t424 = -pkin(6) * t533 - t142 - t218;
t423 = -pkin(1) - t668;
t274 = t306 * qJD(1);
t422 = -t121 - t274 - t503;
t419 = t532 * t547 + t323;
t418 = -qJDD(3) * t380 + t120 * t532 + t121 * t534 + t277 * t252 + t377 * t524;
t411 = qJDD(1) * t257 + t554 + t678 * t378 + (t120 + t323) * qJD(1);
t408 = -t111 * t262 + t113 * t263 + t346 * t400;
t407 = (Icges(5,5) * t230 - Icges(5,6) * t231) * t262 - (Icges(5,5) * t232 + Icges(5,6) * t233) * t263 + t238 * t346;
t404 = t380 * t407;
t394 = t34 * t440 + (-t378 * t39 + t381 * t40) * t208;
t385 = (t400 * t381 + t442) * t262 - (t400 * t378 + t441) * t263 + (t181 + t438) * t346;
t384 = t385 * t380;
t356 = pkin(3) * t536;
t271 = t303 * qJD(2);
t259 = t297 * t537;
t229 = t501 - t510;
t227 = t542 * t535;
t179 = -pkin(6) * t677 + t356;
t178 = qJD(1) * t275 - t316;
t167 = -rSges(5,3) * t585 + t548;
t166 = -rSges(5,3) * t586 + t549;
t164 = t402 * t381;
t163 = t402 * t378;
t162 = t401 * t381;
t161 = t401 * t378;
t156 = rSges(5,1) * t232 + rSges(5,2) * t233;
t155 = rSges(5,1) * t230 - rSges(5,2) * t231;
t146 = -rSges(4,3) * t510 + t479;
t145 = t470 * t534 + (t381 * t471 + t366) * qJD(1);
t84 = t432 * qJD(2);
t67 = -t258 + (qJD(2) * t470 + t358) * t378 + t126 * qJD(1);
t66 = qJD(1) * t513 + t419;
t64 = (t210 * t381 - t211 * t378) * qJD(2) + t475;
t60 = -rSges(5,3) * t677 + t631;
t59 = rSges(5,3) * t410 + t474;
t52 = qJD(1) * t143 + qJDD(1) * t209 - t271 * t534 - t277 * t299 + t554;
t51 = -t271 * t532 + t278 * t299 + t562 * qJDD(1) + (-t144 - t274) * qJD(1);
t29 = qJD(1) * t146 + qJDD(1) * t210 + t277 * t547 + t378 * t483 + t411;
t28 = -t278 * t470 + t381 * t483 + t513 * qJDD(1) + (-t145 + t422) * qJD(1) + t515;
t17 = -t211 * t277 + t561 * t278 + (t145 * t378 + t146 * t381) * qJD(2) + t418;
t16 = -t102 * t401 - t103 * t402 + t127 * t581 + t130 * t230 + t133 * t231 + t400 * t677;
t15 = -t100 * t401 - t101 * t402 + t127 * t583 + t130 * t232 - t133 * t233 - t400 * t410;
t14 = t262 * t41 - t263 * t42 + t346 * t71;
t11 = -t218 * t534 + qJD(1) * t179 + qJDD(1) * t275 + t122 * t261 - t142 * t262 - t147 * t208 - t277 * t297 + t346 * t60 + (-t277 * t377 - t378 * t580) * pkin(6) + t411;
t10 = -t218 * t532 - t124 * t261 - t142 * t263 + t148 * t208 - t346 * t59 + (t278 * t377 - t381 * t580) * pkin(6) + t512 * qJDD(1) + (-t178 + t422) * qJD(1) + t515;
t7 = t102 * t116 - t103 * t118 - t113 * t677 + t230 * t55 + t231 * t57 + t53 * t581;
t6 = t102 * t114 + t103 * t117 - t111 * t677 + t230 * t56 + t231 * t58 + t54 * t581;
t5 = t100 * t116 - t101 * t118 + t113 * t410 + t232 * t55 - t233 * t57 + t53 * t583;
t4 = t100 * t114 + t101 * t117 + t111 * t410 + t232 * t56 - t233 * t58 + t54 * t583;
t3 = -t122 * t148 + t124 * t147 + t262 * t59 + t263 * t60 - t276 * t277 + t555 * t278 + (t178 * t378 + t179 * t381) * qJD(2) + t418;
t2 = t147 * t35 + t148 * t36 + t16 * t346 + t261 * t62 + t262 * t6 - t263 * t7;
t1 = t147 * t37 + t148 * t38 + t15 * t346 + t261 * t63 + t262 * t4 - t263 * t5;
t18 = [(t15 + t12) * t648 + (-(qJD(1) * t276 - t39 + t670 + t717) * t40 + t10 * (-t347 + t702) + t39 * (t316 + t317 - t474) + t40 * (t356 + t511 + t631) + (qJD(1) * t39 * t667 + t40 * t521 * t535 + t10 * pkin(3)) * t381 + (t10 * (t423 - t367) + t39 * (rSges(5,3) * t535 - qJ(3) * t533 - t358) + (t39 * (-pkin(3) - pkin(5)) + t667 * t40) * qJD(1)) * t378 - (t276 + (t380 * t640 + t491) * t378 + t702) * g(1) + (t11 - g(2)) * (t478 + t575)) * m(5) + t616 / 0.2e1 + t617 / 0.2e1 - m(2) * (-g(1) * t300 + g(2) * t304) + ((t674 * t378 + ((t736 + t749) * t381 + t693 + t734 + t748) * t381) * qJD(2) + t703) * t493 + (t696 + t697) * t495 + (t695 - t698 + t699) * t494 + t636 / 0.2e1 - t635 / 0.2e1 + t632 + (t684 + t689) * t646 + (t104 + t690) * t645 - t278 * t107 / 0.2e1 + (t723 * qJD(2) + t726 * t377 + t727 * t380) * qJD(1) + (((t381 * t679 - t674 + t691) * t381 + (t378 * t679 + t692 + t735) * t378) * qJD(2) + t700 - t701) * t496 + (m(2) * (t300 ^ 2 + t304 ^ 2) + Icges(2,3) - t722) * qJDD(1) + t16 * t649 + t63 * t652 + t62 * t653 + ((t317 + (-t358 + (t380 * t614 - t628) * qJD(2)) * t378 + ((t377 * t614 + t476) * t381 + (-rSges(4,1) - pkin(5)) * t378) * qJD(1)) * t66 + (-qJD(1) * t211 - t419 + t66 - t670 - pkin(2) * t505 + t479 + t511 + (t423 - t364) * t537) * t67 + (t29 - g(2)) * t126 + (t28 - g(1)) * (t372 + (t476 - t359) * t378 + t488)) * m(4) + (t91 * (t355 + t544) + (t299 * t620 - t625) * qJD(2) + ((-pkin(1) - t303) * t618 + (t90 * (-rSges(3,3) - pkin(5)) + t91 * t500) * t378) * qJD(1) - (-qJD(1) * t207 - t280 - t507 - t90) * t91 + (t52 - g(2)) * t158 + (t51 - g(1)) * (t500 * t378 + t372 + t543)) * m(3) + t12 * t647; (((-t162 * t379 - t164 * t376 + t111) * t262 - (-t161 * t379 - t163 * t376 + t113) * t263 + (-t185 * t379 - t189 * t376 - t400) * t346 + t71 * qJD(4)) * t380 + (-qJD(4) * t460 + t385) * t377) * t644 + (-g(1) * (t327 + t548) - g(2) * (t325 + t549) - g(3) * (t206 - t669) - t521 * t664 - t40 * (t122 * t528 + t167 * t346 - t206 * t262 + t557) - ((-t542 * t595 - t517) * pkin(6) + t394 * qJD(4)) * t377 - t666 * qJD(2) * t669 + t3 * t560 + (t11 * t477 + t3 * t574 + t40 * t424) * t378 + (t3 * t575 + (qJD(1) * t40 + t10) * t477) * t381 + (t124 * t528 + t166 * t346 + t206 * t263 + t208 * t537 + t381 * t424 + t259 - t482) * t39 + (-t166 * t262 - t167 * t263 - t514 + t516 + (t178 + t59 + (-t122 + t555) * qJD(1)) * t378 + (qJD(1) * t574 + t179 + t60) * t381) * t34) * m(5) + ((t378 * t692 + t381 * t691) * qJD(1) + t680) * t495 + ((t378 * t694 + t381 * t693) * qJD(1) + t681) * t494 + (qJD(1) * t696 + qJD(2) * t680 + qJDD(1) * t684 + t277 * t691 + t278 * t692 + t2) * t378 / 0.2e1 + (t698 * t381 + t697 * t378 + (t690 * t378 + t689 * t381) * qJD(1)) * qJD(1) / 0.2e1 + (t12 + t699) * t497 + (g(1) * t256 + g(2) * t251 - g(3) * t303 - (t251 * t90 - t625) * qJD(1) - (t84 * (-t251 * t378 - t256 * t381) + t454 * t303) * qJD(2) + (qJD(2) * t439 + t207 * t277 - t209 * t278) * t432 + t84 * ((t207 * t381 - t209 * t378) * qJD(1) + t439) + t454 * t271 + (-t52 * t378 - t51 * t381 + (-t381 * t91 + t620) * qJD(1)) * t299) * m(3) + (t66 * t259 + t17 * t560 + t64 * t516 + (t28 * t547 + t66 * t559 + t17 * t210 + t64 * t146 + (-t64 * t211 + t547 * t67) * qJD(1)) * t381 + (t29 * t547 + t67 * t559 - t17 * t211 + t64 * t145 + (-t470 * t66 + t561 * t64) * qJD(1)) * t378 - t66 * (-qJD(1) * t250 + t482) - t67 * (qJD(1) * t255 + t557) - t64 * t514 - ((t64 * t255 + t546 * t66) * t381 + (t64 * t250 + t546 * t67) * t378) * qJD(2) - g(1) * (t254 + t255) - g(2) * (t249 + t250) + g(3) * t546) * m(4) + ((-t532 * t672 - t675) * t381 + ((t381 * t673 + t683) * qJD(2) + t682) * t378) * t493 + ((-t534 * t673 + t675) * t378 + ((t378 * t672 + t683) * qJD(2) + t682) * t381) * t496 + (t689 * t378 - t690 * t381) * qJDD(1) / 0.2e1 - t14 * t528 / 0.2e1 + (t13 + t700) * t537 / 0.2e1 + (qJD(1) * t460 + t378 * t8 - t381 * t9) * t643 + ((t162 * t232 - t164 * t233) * t262 - (t161 * t232 - t163 * t233) * t263 + (t185 * t232 - t189 * t233) * t346 + (-t37 * t585 + t380 * t63) * qJD(4) + ((-qJD(4) * t38 + t408) * t377 + t384) * t378) * t647 + (qJD(1) * t462 + t378 * t4 - t381 * t5) * t648 + (qJD(1) * t464 + t378 * t6 - t381 * t7) * t649 + ((t162 * t230 + t164 * t231) * t262 - (t161 * t230 + t163 * t231) * t263 + (t185 * t230 + t189 * t231) * t346 + (-t36 * t586 + t380 * t62) * qJD(4) + ((-qJD(4) * t35 + t408) * t377 + t384) * t381) * t650 + t459 * t651 + t461 * t652 + t463 * t653 + t710 * t645 + t711 * t646 - ((t712 * t377 + t707 * t380) * qJD(2) + (t708 * t377 + t709 * t380) * qJD(1)) * qJD(1) / 0.2e1 - (t695 * qJD(1) + t681 * qJD(2) + t719 * qJDD(1) + t693 * t277 + t694 * t278 + t1) * t381 / 0.2e1 + (t12 * t381 + t13 * t378) * t529 / 0.2e1; (-m(4) - m(5)) * (-g(3) * t380 + t664) - m(4) * (t227 * t64 + t228 * t67 + t229 * t66) - m(5) * (t227 * t34 + t228 * t40 + t229 * t39) + 0.2e1 * ((t532 * t66 + t534 * t67 - t17) * t655 + (t39 * t532 + t40 * t534 - t3) * t654) * t380 + 0.2e1 * ((qJD(2) * t64 + t28 * t381 + t29 * t378 + t536 * t67 - t537 * t66) * t655 + (t10 * t381 + t11 * t378 - t39 * t537 + t517 + t595) * t654) * t377; t2 * t581 / 0.2e1 + (t377 * t62 + t380 * t464) * t653 + ((-qJD(2) * t464 + t16) * t377 + (-qJD(1) * t463 + qJD(2) * t62 + t378 * t7 + t381 * t6) * t380) * t649 + t1 * t583 / 0.2e1 + (t377 * t63 + t380 * t462) * t652 + ((-qJD(2) * t462 + t15) * t377 + (-qJD(1) * t461 + qJD(2) * t63 + t378 * t5 + t381 * t4) * t380) * t648 + t14 * t533 / 0.2e1 + t377 * (t616 + t617 + t632 - t635 + t636) / 0.2e1 + (t377 * t71 + t380 * t460) * t651 + ((-qJD(2) * t460 + t19) * t377 + (-qJD(1) * t459 + qJD(2) * t71 + t378 * t9 + t381 * t8) * t380) * t643 + (t395 * t230 - t231 * t396 + t381 * t404) * t650 + (t232 * t395 + t233 * t396 + t378 * t404) * t647 + (t407 * t377 + (t396 * t376 - t379 * t395) * t380) * t644 + (t377 * t496 + t380 * t497) * t13 + (-t509 / 0.2e1 + t377 * t494) * t12 + ((qJD(2) * t394 - t10 * t124 + t11 * t122 - t39 * t59 + t40 * t60) * t377 + (t39 * (-qJD(2) * t124 + t142 * t378) + t40 * (qJD(2) * t122 - t142 * t381) - t3 * t440 + t34 * (-t122 * t536 - t124 * t537 - t378 * t60 + t381 * t59) + (qJD(1) * t666 + t10 * t378 - t11 * t381) * t208) * t380 - t39 * (-t156 * t346 - t253 * t263) - t40 * (t155 * t346 - t253 * t262) - t34 * (t155 * t263 + t156 * t262) - g(1) * t155 - g(2) * t156 - g(3) * t253) * m(5);];
tau = t18;
