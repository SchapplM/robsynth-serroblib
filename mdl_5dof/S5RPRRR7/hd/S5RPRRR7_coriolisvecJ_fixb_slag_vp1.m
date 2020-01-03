% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:12
% EndTime: 2019-12-31 19:04:03
% DurationCPUTime: 39.93s
% Computational Cost: add. (37907->1144), mult. (40994->1571), div. (0->0), fcn. (39117->10), ass. (0->539)
t429 = qJ(1) + pkin(9);
t422 = sin(t429);
t431 = sin(qJ(4));
t435 = cos(qJ(3));
t662 = t431 * t435;
t423 = cos(t429);
t434 = cos(qJ(4));
t666 = t423 * t434;
t314 = t422 * t662 + t666;
t659 = t434 * t435;
t668 = t423 * t431;
t315 = t422 * t659 - t668;
t432 = sin(qJ(3));
t673 = t422 * t432;
t174 = Icges(5,5) * t315 - Icges(5,6) * t314 + Icges(5,3) * t673;
t304 = Icges(5,4) * t315;
t177 = -Icges(5,2) * t314 + Icges(5,6) * t673 + t304;
t303 = Icges(5,4) * t314;
t181 = -Icges(5,1) * t315 - Icges(5,5) * t673 + t303;
t793 = t177 * t431 + t181 * t434;
t75 = -t174 * t435 - t432 * t793;
t430 = qJ(4) + qJ(5);
t424 = sin(t430);
t664 = t424 * t435;
t425 = cos(t430);
t669 = t423 * t425;
t289 = t422 * t664 + t669;
t663 = t425 * t435;
t670 = t423 * t424;
t290 = t422 * t663 - t670;
t152 = Icges(6,5) * t290 - Icges(6,6) * t289 + Icges(6,3) * t673;
t281 = Icges(6,4) * t290;
t155 = -Icges(6,2) * t289 + Icges(6,6) * t673 + t281;
t280 = Icges(6,4) * t289;
t159 = -Icges(6,1) * t290 - Icges(6,5) * t673 + t280;
t794 = t155 * t424 + t159 * t425;
t71 = -t152 * t435 - t432 * t794;
t517 = Icges(6,5) * t425 - Icges(6,6) * t424;
t293 = -Icges(6,3) * t435 + t432 * t517;
t698 = Icges(6,4) * t425;
t519 = -Icges(6,2) * t424 + t698;
t295 = -Icges(6,6) * t435 + t432 * t519;
t699 = Icges(6,4) * t424;
t522 = Icges(6,1) * t425 - t699;
t297 = -Icges(6,5) * t435 + t432 * t522;
t104 = -t289 * t295 + t290 * t297 + t293 * t673;
t410 = qJD(3) * t422;
t618 = qJD(4) * t432;
t344 = t423 * t618 + t410;
t616 = qJD(5) * t432;
t263 = t423 * t616 + t344;
t621 = qJD(3) * t423;
t345 = -t422 * t618 + t621;
t264 = -t422 * t616 + t345;
t63 = t152 * t673 - t155 * t289 - t159 * t290;
t675 = t422 * t425;
t291 = -t423 * t664 + t675;
t676 = t422 * t424;
t292 = t423 * t663 + t676;
t667 = t423 * t432;
t154 = Icges(6,5) * t292 + Icges(6,6) * t291 + Icges(6,3) * t667;
t700 = Icges(6,4) * t292;
t157 = Icges(6,2) * t291 + Icges(6,6) * t667 + t700;
t282 = Icges(6,4) * t291;
t160 = Icges(6,1) * t292 + Icges(6,5) * t667 + t282;
t64 = t154 * t673 - t289 * t157 + t290 * t160;
t428 = qJD(4) + qJD(5);
t558 = t428 * t435;
t769 = qJD(1) - t558;
t35 = t104 * t769 + t263 * t64 - t264 * t63;
t105 = t291 * t295 + t292 * t297 + t293 * t667;
t65 = t152 * t667 + t291 * t155 - t159 * t292;
t66 = t154 * t667 + t291 * t157 + t292 * t160;
t36 = t105 * t769 + t263 * t66 - t264 * t65;
t518 = Icges(5,5) * t434 - Icges(5,6) * t431;
t318 = -Icges(5,3) * t435 + t432 * t518;
t701 = Icges(5,4) * t434;
t520 = -Icges(5,2) * t431 + t701;
t320 = -Icges(5,6) * t435 + t432 * t520;
t702 = Icges(5,4) * t431;
t523 = Icges(5,1) * t434 - t702;
t322 = -Icges(5,5) * t435 + t432 * t523;
t109 = -t314 * t320 + t315 * t322 + t318 * t673;
t617 = qJD(4) * t435;
t411 = qJD(1) - t617;
t67 = t174 * t673 - t177 * t314 - t181 * t315;
t672 = t422 * t434;
t316 = -t423 * t662 + t672;
t674 = t422 * t431;
t317 = t423 * t659 + t674;
t176 = Icges(5,5) * t317 + Icges(5,6) * t316 + Icges(5,3) * t667;
t703 = Icges(5,4) * t317;
t179 = Icges(5,2) * t316 + Icges(5,6) * t667 + t703;
t305 = Icges(5,4) * t316;
t182 = Icges(5,1) * t317 + Icges(5,5) * t667 + t305;
t68 = t176 * t673 - t314 * t179 + t315 * t182;
t38 = t109 * t411 + t344 * t68 - t345 * t67;
t110 = t316 * t320 + t317 * t322 + t318 * t667;
t69 = t174 * t667 + t316 * t177 - t181 * t317;
t70 = t176 * t667 + t316 * t179 + t317 * t182;
t39 = t110 * t411 + t344 * t70 - t345 * t69;
t438 = qJD(1) ^ 2;
t786 = t422 * t423;
t415 = t422 * rSges(4,3);
t665 = t423 * t435;
t272 = rSges(4,1) * t665 - rSges(4,2) * t667 + t415;
t358 = t423 * pkin(2) + t422 * pkin(6);
t436 = cos(qJ(1));
t427 = t436 * pkin(1);
t774 = t427 + t358;
t785 = t272 + t774;
t619 = qJD(3) * t435;
t573 = t619 / 0.2e1;
t623 = qJD(1) * t432;
t598 = t422 * t623;
t784 = t423 * t573 - t598 / 0.2e1;
t624 = qJD(1) * t423;
t579 = t624 / 0.2e1;
t783 = t422 * t573 + t432 * t579;
t437 = -pkin(8) - pkin(7);
t732 = pkin(7) + t437;
t584 = t732 * t435;
t421 = pkin(4) * t434 + pkin(3);
t733 = pkin(3) - t421;
t585 = t733 * t432;
t287 = t584 - t585;
t597 = t423 * t623;
t479 = t422 * t619 + t597;
t538 = rSges(6,1) * t290 - rSges(6,2) * t289;
t161 = rSges(6,3) * t673 + t538;
t658 = t435 * t421;
t355 = t422 * t658;
t612 = pkin(4) * t668;
t734 = t435 * pkin(3);
t777 = t432 * t732;
t191 = t612 - t355 + (t734 + t777) * t422;
t537 = rSges(6,1) * t425 - rSges(6,2) * t424;
t299 = -rSges(6,3) * t435 + t432 * t537;
t400 = pkin(3) * t432 - pkin(7) * t435;
t595 = t400 * t621;
t782 = -t161 * t769 + t191 * t411 - t264 * t299 - t287 * t345 - t595;
t541 = rSges(5,1) * t315 - rSges(5,2) * t314;
t185 = rSges(5,3) * t673 + t541;
t540 = rSges(5,1) * t434 - rSges(5,2) * t431;
t332 = -rSges(5,3) * t435 + t432 * t540;
t781 = -t185 * t411 - t332 * t345 - t595;
t780 = 0.2e1 * qJD(3);
t163 = t292 * rSges(6,1) + t291 * rSges(6,2) + rSges(6,3) * t667;
t404 = pkin(3) * t665;
t337 = pkin(7) * t667 + t404;
t402 = pkin(4) * t674;
t633 = t423 * t658 + t402;
t660 = t432 * t437;
t192 = -t423 * t660 - t337 + t633;
t735 = pkin(7) * t432;
t401 = t734 + t735;
t335 = t401 * t422;
t587 = t335 * t410 + t337 * t621 + qJD(2);
t50 = t161 * t263 + t163 * t264 - t191 * t344 + t192 * t345 + t587;
t653 = t161 - t191;
t779 = t50 * t653;
t776 = t435 * t733;
t566 = t423 * rSges(3,1) - rSges(3,2) * t422;
t426 = Icges(4,4) * t435;
t521 = -Icges(4,2) * t432 + t426;
t384 = Icges(4,1) * t432 + t426;
t773 = t427 + t566;
t419 = t423 * pkin(6);
t357 = pkin(2) * t422 - t419;
t349 = qJD(1) * t357;
t433 = sin(qJ(1));
t737 = pkin(1) * t433;
t768 = -t349 + (-t335 - t737) * qJD(1);
t671 = t422 * t435;
t690 = Icges(4,3) * t423;
t265 = Icges(4,5) * t671 - Icges(4,6) * t673 - t690;
t395 = Icges(4,4) * t673;
t697 = Icges(4,5) * t423;
t269 = Icges(4,1) * t671 - t395 - t697;
t693 = Icges(4,6) * t423;
t267 = Icges(4,4) * t671 - Icges(4,2) * t673 - t693;
t681 = t267 * t432;
t510 = -t269 * t435 + t681;
t112 = -t265 * t423 - t422 * t510;
t381 = Icges(4,5) * t435 - Icges(4,6) * t432;
t380 = Icges(4,5) * t432 + Icges(4,6) * t435;
t483 = qJD(3) * t380;
t704 = Icges(4,4) * t432;
t385 = Icges(4,1) * t435 - t704;
t270 = Icges(4,5) * t422 + t385 * t423;
t268 = Icges(4,6) * t422 + t423 * t521;
t680 = t268 * t432;
t509 = -t270 * t435 + t680;
t767 = -t423 * t483 + (-t381 * t422 + t509 + t690) * qJD(1);
t266 = Icges(4,3) * t422 + t381 * t423;
t766 = -t422 * t483 + (t266 + t510) * qJD(1);
t382 = Icges(4,2) * t435 + t704;
t505 = t382 * t432 - t384 * t435;
t765 = qJD(1) * t505 + t381 * qJD(3);
t590 = t423 * t619;
t376 = pkin(7) * t590;
t725 = pkin(4) * qJD(4);
t607 = t431 * t725;
t557 = t435 * t607;
t606 = t434 * t725;
t600 = qJD(1) * t612 + t422 * t606 + t437 * t598;
t625 = qJD(1) * t422;
t657 = t435 * t437;
t102 = -t376 + (t735 + t776) * t625 + (-t557 + (t585 - t657) * qJD(3)) * t423 + t600;
t288 = -t776 - t777;
t620 = qJD(3) * t432;
t375 = t422 * pkin(3) * t620;
t661 = t432 * t421;
t103 = -t423 * t606 + t375 + (-t557 + (-t584 - t661) * qJD(3)) * t422 + (t288 * t423 + t402) * qJD(1);
t614 = qJD(1) * qJD(3);
t405 = t422 * t614;
t613 = qJD(3) * qJD(4);
t581 = t435 * t613;
t583 = qJD(1) * t618;
t259 = t422 * t581 + t423 * t583 + t405;
t183 = qJD(5) * t479 + t259;
t406 = t423 * t614;
t632 = t423 * t581 + t406;
t184 = qJD(5) * t590 - t428 * t598 + t632;
t260 = -t422 * t583 + t632;
t591 = t423 * t620;
t622 = qJD(1) * t435;
t478 = -t422 * t622 - t591;
t228 = pkin(3) * t478 - pkin(7) * t598 + t376;
t229 = pkin(7) * t479 + qJD(1) * t404 - t375;
t471 = t228 * t621 + t229 * t410 + t335 * t406 - t337 * t405;
t464 = t424 * t620 + t425 * t769;
t559 = -t428 + t622;
t134 = t423 * t464 + t559 * t676;
t463 = t424 * t769 - t425 * t620;
t135 = t423 * t463 - t559 * t675;
t605 = t135 * rSges(6,1) + t134 * rSges(6,2) + rSges(6,3) * t590;
t87 = -rSges(6,3) * t598 + t605;
t136 = t422 * t464 - t559 * t670;
t137 = t422 * t463 + t559 * t669;
t539 = -rSges(6,1) * t137 - rSges(6,2) * t136;
t88 = rSges(6,3) * t479 - t539;
t11 = t102 * t345 + t103 * t344 + t161 * t184 - t163 * t183 - t191 * t260 - t192 * t259 + t263 * t88 + t264 * t87 + t471;
t652 = t163 + t192;
t709 = t102 + t87;
t764 = -t11 * t652 - t50 * t709;
t642 = -Icges(4,2) * t671 + t269 - t395;
t644 = t384 * t422 + t267;
t763 = -t432 * t642 - t435 * t644;
t294 = Icges(6,3) * t432 + t435 * t517;
t507 = -t295 * t424 + t297 * t425;
t514 = -t157 * t424 + t160 * t425;
t762 = t263 * (-t293 * t423 - t514) - t264 * (-t293 * t422 + t794) + t769 * (t294 - t507);
t339 = (-Icges(6,2) * t425 - t699) * t432;
t761 = t263 * (-Icges(6,2) * t292 + t160 + t282) - t264 * (-Icges(6,2) * t290 - t159 - t280) + t769 * (t297 + t339);
t319 = Icges(5,3) * t432 + t435 * t518;
t506 = -t320 * t431 + t322 * t434;
t512 = -t179 * t431 + t182 * t434;
t760 = t344 * (-t318 * t423 - t512) - t345 * (-t318 * t422 + t793) + t411 * (t319 - t506);
t351 = (-Icges(5,2) * t434 - t702) * t432;
t759 = t344 * (-Icges(5,2) * t317 + t182 + t305) - t345 * (-Icges(5,2) * t315 - t181 - t303) + t411 * (t322 + t351);
t758 = t183 / 0.2e1;
t757 = t184 / 0.2e1;
t756 = t259 / 0.2e1;
t755 = t260 / 0.2e1;
t754 = -t263 / 0.2e1;
t753 = t263 / 0.2e1;
t752 = -t264 / 0.2e1;
t751 = t264 / 0.2e1;
t750 = -t344 / 0.2e1;
t749 = t344 / 0.2e1;
t748 = -t345 / 0.2e1;
t747 = t345 / 0.2e1;
t378 = t428 * t432;
t359 = qJD(3) * t378;
t746 = t359 / 0.2e1;
t745 = -t769 / 0.2e1;
t744 = t769 / 0.2e1;
t743 = -t411 / 0.2e1;
t742 = t411 / 0.2e1;
t739 = -t435 / 0.2e1;
t738 = -rSges(5,3) - pkin(7);
t736 = pkin(4) * t431;
t731 = rSges(4,1) * t435;
t729 = rSges(5,3) * t432;
t727 = rSges(6,3) * t432;
t82 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t479;
t84 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t479;
t86 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t479;
t18 = (-qJD(3) * t794 - t82) * t435 + (qJD(3) * t152 + (-t155 * t428 + t86) * t425 + (t159 * t428 - t84) * t424) * t432;
t724 = t18 * t264;
t477 = t590 - t598;
t81 = Icges(6,5) * t135 + Icges(6,6) * t134 + Icges(6,3) * t477;
t83 = Icges(6,4) * t135 + Icges(6,2) * t134 + Icges(6,6) * t477;
t85 = Icges(6,1) * t135 + Icges(6,4) * t134 + Icges(6,5) * t477;
t19 = (qJD(3) * t514 - t81) * t435 + (qJD(3) * t154 + (-t157 * t428 + t85) * t425 + (-t160 * t428 - t83) * t424) * t432;
t723 = t19 * t263;
t462 = t411 * t434 + t431 * t620;
t555 = -qJD(4) + t622;
t171 = t422 * t462 - t555 * t668;
t461 = t411 * t431 - t434 * t620;
t172 = t422 * t461 + t555 * t666;
t95 = Icges(5,5) * t172 + Icges(5,6) * t171 + Icges(5,3) * t479;
t97 = Icges(5,4) * t172 + Icges(5,2) * t171 + Icges(5,6) * t479;
t99 = Icges(5,1) * t172 + Icges(5,4) * t171 + Icges(5,5) * t479;
t27 = (-qJD(3) * t793 - t95) * t435 + (qJD(3) * t174 - t431 * t97 + t434 * t99 + (-t177 * t434 + t181 * t431) * qJD(4)) * t432;
t720 = t27 * t345;
t169 = t423 * t462 + t555 * t674;
t170 = t423 * t461 - t555 * t672;
t94 = Icges(5,5) * t170 + Icges(5,6) * t169 + Icges(5,3) * t477;
t96 = Icges(5,4) * t170 + Icges(5,2) * t169 + Icges(5,6) * t477;
t98 = Icges(5,1) * t170 + Icges(5,4) * t169 + Icges(5,5) * t477;
t28 = (qJD(3) * t512 - t94) * t435 + (qJD(3) * t176 - t431 * t96 + t434 * t98 + (-t179 * t434 - t182 * t431) * qJD(4)) * t432;
t719 = t28 * t344;
t571 = -t357 - t737;
t503 = (-t335 + t571) * qJD(1);
t89 = t503 + t781;
t716 = t423 * t89;
t476 = (t337 + t774) * qJD(1) - t400 * t410;
t58 = t163 * t769 + t192 * t411 - t263 * t299 - t287 * t344 + t476;
t715 = t58 * t192;
t714 = t71 * t183;
t72 = -t154 * t435 + t432 * t514;
t713 = t72 * t184;
t712 = t75 * t259;
t76 = -t176 * t435 + t432 * t512;
t711 = t76 * t260;
t710 = -rSges(6,3) + t437;
t116 = -t293 * t435 + t432 * t507;
t338 = (-Icges(6,5) * t424 - Icges(6,6) * t425) * t432;
t188 = qJD(3) * t294 + t338 * t428;
t296 = Icges(6,6) * t432 + t435 * t519;
t189 = qJD(3) * t296 + t339 * t428;
t298 = Icges(6,5) * t432 + t435 * t522;
t340 = (-Icges(6,1) * t424 - t698) * t432;
t190 = qJD(3) * t298 + t340 * t428;
t56 = (qJD(3) * t507 - t188) * t435 + (qJD(3) * t293 + (-t295 * t428 + t190) * t425 + (-t297 * t428 - t189) * t424) * t432;
t708 = t116 * t359 + t56 * t769;
t125 = -t318 * t435 + t432 * t506;
t582 = t432 * t613;
t350 = (-Icges(5,5) * t431 - Icges(5,6) * t434) * t432;
t230 = qJD(3) * t319 + qJD(4) * t350;
t321 = Icges(5,6) * t432 + t435 * t520;
t231 = qJD(3) * t321 + qJD(4) * t351;
t323 = Icges(5,5) * t432 + t435 * t523;
t352 = (-Icges(5,1) * t431 - t701) * t432;
t232 = qJD(3) * t323 + qJD(4) * t352;
t62 = (qJD(3) * t506 - t230) * t435 + (qJD(3) * t318 - t231 * t431 + t232 * t434 + (-t320 * t434 - t322 * t431) * qJD(4)) * t432;
t707 = t125 * t582 + t62 * t411;
t706 = t161 * t590 + t88 * t667;
t628 = rSges(4,2) * t673 + t423 * rSges(4,3);
t271 = rSges(4,1) * t671 - t628;
t389 = rSges(4,1) * t432 + rSges(4,2) * t435;
t592 = t389 * t621;
t144 = -t592 + (-t271 + t571) * qJD(1);
t685 = t144 * t422;
t684 = t144 * t423;
t596 = t389 * t410;
t145 = qJD(1) * t785 - t596;
t331 = t389 * t423;
t683 = t145 * t331;
t366 = qJD(3) * t401;
t679 = t366 * t423;
t678 = t380 * t422;
t677 = t380 * t423;
t656 = t163 * t620 + t299 * t598;
t187 = t317 * rSges(5,1) + t316 * rSges(5,2) + rSges(5,3) * t667;
t649 = -t187 - t337;
t300 = t435 * t537 + t727;
t341 = (-rSges(6,1) * t424 - rSges(6,2) * t425) * t432;
t201 = qJD(3) * t300 + t341 * t428;
t240 = qJD(3) * t288 - t432 * t607;
t648 = -t201 - t240;
t333 = t435 * t540 + t729;
t353 = (-rSges(5,1) * t431 - rSges(5,2) * t434) * t432;
t233 = qJD(3) * t333 + qJD(4) * t353;
t647 = -t233 - t366;
t646 = -t422 * t265 - t269 * t665;
t645 = t422 * t266 + t270 * t665;
t643 = -t384 * t423 - t268;
t641 = -t382 * t423 + t270;
t334 = t400 * t422;
t336 = t400 * t423;
t640 = -t334 * t410 - t336 * t621;
t638 = t422 * t335 + t423 * t337;
t637 = t287 + t299;
t634 = -t332 - t400;
t631 = rSges(4,2) * t598 + rSges(4,3) * t624;
t630 = -t382 + t385;
t629 = t384 + t521;
t626 = qJD(1) * t381;
t199 = -t422 * t505 - t677;
t615 = t199 * qJD(1);
t611 = t432 * t736;
t610 = t438 * t737;
t609 = t438 * t427;
t604 = t170 * rSges(5,1) + t169 * rSges(5,2) + rSges(5,3) * t590;
t603 = -t366 + t648;
t602 = t423 * t228 + t422 * t229 + t335 * t624;
t601 = -t400 - t637;
t589 = t673 / 0.2e1;
t588 = t667 / 0.2e1;
t586 = -pkin(2) - t731;
t578 = -t410 / 0.2e1;
t575 = t621 / 0.2e1;
t574 = t620 / 0.2e1;
t569 = t419 - t737;
t568 = -pkin(2) - t658;
t210 = -rSges(6,1) * t289 - rSges(6,2) * t290;
t211 = rSges(6,1) * t291 - rSges(6,2) * t292;
t565 = t263 * t210 + t211 * t264;
t564 = t211 * t769 - t263 * t341;
t563 = -t210 * t769 - t264 * t341;
t243 = t270 * t671;
t562 = t266 * t423 - t243;
t561 = -t265 + t680;
t556 = t201 * t673 + t299 * t479 + t435 * t88;
t409 = pkin(6) * t624;
t550 = qJD(1) * (-pkin(2) * t625 + t409) - t610;
t547 = qJD(4) * t574;
t546 = -qJD(1) * t336 - t401 * t410;
t354 = rSges(3,1) * t422 + rSges(3,2) * t423;
t543 = -rSges(4,2) * t432 + t731;
t542 = rSges(5,1) * t172 + rSges(5,2) * t171;
t536 = t422 * t64 - t423 * t63;
t535 = t422 * t63 + t423 * t64;
t534 = t422 * t66 - t423 * t65;
t533 = t422 * t65 + t423 * t66;
t532 = t422 * t68 - t423 * t67;
t531 = t422 * t67 + t423 * t68;
t530 = t422 * t70 - t423 * t69;
t529 = t422 * t69 + t423 * t70;
t528 = t422 * t72 - t423 * t71;
t527 = t422 * t71 + t423 * t72;
t526 = t422 * t76 - t423 * t75;
t525 = t422 * t75 + t423 * t76;
t516 = -t145 * t422 - t684;
t511 = t185 * t423 - t187 * t422;
t164 = t267 * t435 + t269 * t432;
t165 = t268 * t435 + t270 * t432;
t508 = t271 * t422 + t272 * t423;
t504 = qJD(1) * t228 + t550;
t492 = -t366 * t422 - t400 * t624;
t330 = t389 * t422;
t346 = t358 * qJD(1);
t469 = t400 * t405 + (-t229 - t346) * qJD(1) - t609;
t22 = -t103 * t411 - t161 * t359 + t183 * t299 - t201 * t264 - t240 * t345 + t259 * t287 - t769 * t88 + (t191 * t618 - t679) * qJD(3) + t469;
t490 = t11 * t161 * t667 + t22 * (t435 * t161 + t299 * t673);
t57 = t503 + t782;
t489 = t57 * t287 - t50 * t652;
t488 = qJD(1) * t334 - t401 * t621;
t486 = t432 * t738 - pkin(2) - t734;
t485 = qJD(3) * t384;
t484 = qJD(3) * t382;
t113 = -t268 * t673 - t562;
t482 = (-t112 * t423 + t113 * t422) * qJD(3);
t114 = -t267 * t667 - t646;
t115 = -t268 * t667 + t645;
t481 = (-t114 * t423 + t115 * t422) * qJD(3);
t474 = -t174 * t345 + t176 * t344 + t318 * t411;
t473 = (-Icges(6,5) * t289 - Icges(6,6) * t290) * t264 - (Icges(6,5) * t291 - Icges(6,6) * t292) * t263 - t338 * t769;
t472 = (-Icges(5,5) * t314 - Icges(5,6) * t315) * t345 - (Icges(5,5) * t316 - Icges(5,6) * t317) * t344 - t350 * t411;
t470 = -t432 * t641 + t435 * t643;
t468 = t432 * t473;
t467 = t432 * t472;
t454 = (-t432 * t629 + t435 * t630) * qJD(1);
t452 = (Icges(6,1) * t291 - t157 - t700) * t263 - (-Icges(6,1) * t289 - t155 - t281) * t264 + (-t295 + t340) * t769;
t450 = (Icges(5,1) * t316 - t179 - t703) * t344 - (-Icges(5,1) * t314 - t177 - t304) * t345 + (-t320 + t352) * t411;
t208 = rSges(4,1) * t478 - rSges(4,2) * t590 + t631;
t209 = -qJD(3) * t330 + (t423 * t543 + t415) * qJD(1);
t448 = t208 * t423 + t209 * t422 + (t271 * t423 - t272 * t422) * qJD(1);
t14 = t134 * t155 - t135 * t159 + t152 * t477 + t291 * t84 + t292 * t86 + t667 * t82;
t15 = t134 * t157 + t135 * t160 + t154 * t477 + t291 * t83 + t292 * t85 + t667 * t81;
t16 = t136 * t155 - t137 * t159 + t152 * t479 - t289 * t84 + t290 * t86 + t673 * t82;
t17 = t136 * t157 + t137 * t160 + t154 * t479 - t289 * t83 + t290 * t85 + t673 * t81;
t41 = t116 * t769 + t263 * t72 - t264 * t71;
t44 = t134 * t295 + t135 * t297 + t188 * t667 + t189 * t291 + t190 * t292 + t293 * t477;
t45 = t136 * t295 + t137 * t297 + t188 * t673 - t189 * t289 + t190 * t290 + t293 * t479;
t5 = t105 * t359 - t14 * t264 + t15 * t263 + t183 * t65 + t184 * t66 + t44 * t769;
t6 = t104 * t359 - t16 * t264 + t17 * t263 + t183 * t63 + t184 * t64 + t45 * t769;
t447 = ((qJD(3) * t533 - t44) * t435 + (-qJD(1) * t534 + qJD(3) * t105 + t14 * t422 + t15 * t423) * t432) * t753 + (t291 * t761 + t452 * t292 - t423 * t468) * t754 + (-t289 * t761 + t290 * t452 - t422 * t468) * t751 + ((qJD(3) * t535 - t45) * t435 + (-qJD(1) * t536 + qJD(3) * t104 + t16 * t422 + t17 * t423) * t432) * t752 + (t473 * t435 + (-t424 * t761 + t425 * t452) * t432) * t745 + t6 * t589 + (-t104 * t435 + t432 * t535) * t758 + (-t105 * t435 + t432 * t533) * t757 + t5 * t588 + t41 * t574 + (-t116 * t435 + t432 * t527) * t746 + ((qJD(3) * t527 - t56) * t435 + (-qJD(1) * t528 + qJD(3) * t116 + t18 * t422 + t19 * t423) * t432) * t744 + (t708 + t713 + t714 + t723 - t724) * t739 + t784 * t36 + t783 * t35;
t73 = t185 * t344 + t187 * t345 + t587;
t90 = t187 * t411 - t332 * t344 + t476;
t444 = t73 * t511 + (t422 * t89 - t423 * t90) * t332;
t361 = t521 * qJD(3);
t362 = t385 * qJD(3);
t441 = qJD(1) * t380 - t361 * t432 + t362 * t435 + (-t382 * t435 - t384 * t432) * qJD(3);
t440 = t760 * t432;
t439 = (-t152 * t264 + t154 * t263 + t293 * t769) * t435 + t762 * t432;
t363 = t543 * qJD(3);
t348 = t400 * t625;
t313 = t423 * t558;
t312 = t422 * t558;
t307 = t316 * pkin(4);
t306 = t314 * pkin(4);
t253 = t332 * t423;
t252 = t332 * t422;
t251 = t322 * t423;
t250 = t322 * t422;
t249 = t320 * t423;
t248 = t320 * t422;
t242 = t299 * t423;
t241 = t299 * t422;
t239 = t297 * t423;
t238 = t297 * t422;
t237 = t295 * t423;
t236 = t295 * t422;
t227 = t423 * t287;
t226 = t422 * t287;
t225 = rSges(5,1) * t316 - rSges(5,2) * t317;
t224 = -rSges(5,1) * t314 - rSges(5,2) * t315;
t200 = -t423 * t505 + t678;
t173 = t200 * qJD(1);
t143 = qJD(3) * t508 + qJD(2);
t108 = -t609 - t363 * t621 + (-t209 - t346 + t596) * qJD(1);
t107 = -t363 * t410 + (t208 - t592) * qJD(1) + t550;
t101 = rSges(5,3) * t479 + t542;
t100 = -rSges(5,3) * t598 + t604;
t92 = t441 * t422 - t423 * t765;
t91 = t422 * t765 + t441 * t423;
t79 = -qJD(3) * t509 + (-t423 * t484 + (-t422 * t521 + t693) * qJD(1)) * t435 + (-t423 * t485 + (-t385 * t422 + t697) * qJD(1)) * t432;
t78 = -qJD(3) * t510 + (qJD(1) * t268 - t422 * t484) * t435 + (qJD(1) * t270 - t422 * t485) * t432;
t74 = t448 * qJD(3);
t61 = t173 + t481;
t60 = t482 + t615;
t49 = t171 * t320 + t172 * t322 + t230 * t673 - t231 * t314 + t232 * t315 + t318 * t479;
t48 = t169 * t320 + t170 * t322 + t230 * t667 + t231 * t316 + t232 * t317 + t318 * t477;
t47 = -t101 * t411 - t233 * t345 + t259 * t332 + (-t185 * t618 - t679) * qJD(3) + t469;
t46 = t100 * t411 - t233 * t344 - t260 * t332 + (t187 * t618 + t492) * qJD(3) + t504;
t43 = t125 * t411 + t344 * t76 - t345 * t75;
t37 = t100 * t345 + t101 * t344 + t185 * t260 - t187 * t259 + t471;
t26 = t171 * t179 + t172 * t182 + t176 * t479 - t314 * t96 + t315 * t98 + t673 * t94;
t25 = t171 * t177 - t172 * t181 + t174 * t479 - t314 * t97 + t315 * t99 + t673 * t95;
t24 = t169 * t179 + t170 * t182 + t176 * t477 + t316 * t96 + t317 * t98 + t667 * t94;
t23 = t169 * t177 - t170 * t181 + t174 * t477 + t316 * t97 + t317 * t99 + t667 * t95;
t21 = t102 * t411 + t163 * t359 - t184 * t299 - t201 * t263 - t240 * t344 - t260 * t287 + t769 * t87 + (t192 * t618 + t492) * qJD(3) + t504;
t10 = t109 * t582 - t25 * t345 + t259 * t67 + t26 * t344 + t260 * t68 + t411 * t49;
t9 = t110 * t582 - t23 * t345 + t24 * t344 + t259 * t69 + t260 * t70 + t411 * t48;
t1 = [(t748 + t747) * t39 + (t752 + t751) * t36 + (t47 * (-t541 + t569) + t89 * (t375 - t542) + t46 * (t774 - t649) + t90 * (-pkin(3) * t591 + t376 + t409 + t604) + (t89 * t738 * t619 + t47 * t486) * t422 + ((-t433 * t90 - t436 * t89) * pkin(1) + t486 * t716 + (-t89 * pkin(6) + t90 * (-pkin(2) - t401 - t729)) * t422) * qJD(1) - (t768 + t781 - t89) * t90) * m(5) - (t61 + t78 + t92) * t621 / 0.2e1 + (t79 + t91) * t410 / 0.2e1 + m(3) * ((-t354 * t438 - t610) * t773 + (-t609 + (-0.2e1 * t566 - t427 + t773) * t438) * (-t354 - t737)) + (t22 * (-t355 - t538 + t569) + t57 * t539 + t21 * (t774 + t163 + t633) + t58 * (t409 + t600 + t605) + (-t21 * t660 + t58 * (-t657 - t661) * qJD(3) + (t22 * t431 + (t434 * t57 - t58 * t662) * qJD(4)) * pkin(4)) * t423 + (t22 * (-pkin(2) + t660 - t727) + (t557 + (t435 * t710 + t661) * qJD(3)) * t57) * t422 + ((-t433 * t58 - t436 * t57) * pkin(1) + t57 * (t432 * t710 + t568) * t423 + (t57 * (-pkin(6) - t736) + t58 * (t568 - t727)) * t422) * qJD(1) - (-t57 + t768 + t782) * t58) * m(6) + (t108 * (t422 * t586 + t569 + t628) + t107 * t785 + t145 * (t409 + t631) + (t389 * t685 - t683) * qJD(3) + ((-t144 * t436 - t145 * t433) * pkin(1) + (-pkin(2) - t543) * t684 + (t144 * (-rSges(4,3) - pkin(6)) + t145 * t586) * t422) * qJD(1) - (-t592 - t144 - t349 + (-t271 - t737) * qJD(1)) * t145) * m(4) + ((t165 + t200) * t423 + (t164 + t199) * t422) * t614 / 0.2e1 + (t173 + ((t113 - t243 + (t266 + t681) * t423 + t646) * t423 + t645 * t422) * qJD(3)) * t575 + (-t615 + ((t423 * t561 + t115 - t645) * t423 + (t422 * t561 + t114 + t562) * t422) * qJD(3) + t60) * t578 - t724 / 0.2e1 + t711 / 0.2e1 + t712 / 0.2e1 + t707 + t708 + t49 * t748 + t48 * t749 + (-qJD(3) * t505 + t361 * t435 + t362 * t432) * qJD(1) + t713 / 0.2e1 + t714 / 0.2e1 + t719 / 0.2e1 - t720 / 0.2e1 + t45 * t752 + t44 * t753 + t110 * t755 + t109 * t756 + t105 * t757 + t104 * t758 + t723 / 0.2e1; m(4) * t74 + m(5) * t37 + m(6) * t11; (-t57 * (-t161 * t378 + t226 * t411 + t241 * t769 - t264 * t300 - t288 * t345 + t299 * t312 + t488) - t58 * (t163 * t378 - t227 * t411 - t242 * t769 - t263 * t300 - t288 * t344 - t299 * t313 + t546) - t50 * (t161 * t313 - t163 * t312 - t226 * t344 - t227 * t345 - t241 * t263 - t242 * t264 + t640) - ((t191 * t57 + t715) * t432 + (t50 * (-t191 * t423 - t192 * t422) + (t422 * t57 - t423 * t58) * t287) * t435) * qJD(4) + t57 * t348 + t11 * t638 + t50 * t602 + (t22 * t601 + t57 * t603 + (t58 * t601 + t779) * qJD(1) - t764) * t423 + (t21 * t601 + t58 * t603 + t11 * t653 + t50 * (t103 + t88) + (t57 * t637 + t50 * (-t337 - t652)) * qJD(1)) * t422) * m(6) + ((t237 * t289 - t239 * t290) * t263 + t64 * t313 - (t236 * t289 - t238 * t290) * t264 + t63 * t312 + (-t289 * t296 + t290 * t298) * t769 + t104 * t378 + t439 * t422) * t751 + ((-t237 * t291 - t239 * t292) * t263 + t66 * t313 - (-t236 * t291 - t238 * t292) * t264 + t65 * t312 + (t291 * t296 + t292 * t298) * t769 + t105 * t378 + t439 * t423) * t754 + (t116 * t378 + t71 * t312 + t72 * t313 - t762 * t435 + ((t237 * t424 - t239 * t425 + t154) * t263 - (t236 * t424 - t238 * t425 + t152) * t264 + (-t296 * t424 + t298 * t425 + t293) * t769) * t432) * t745 + (((t249 * t431 - t251 * t434 + t176) * t344 - (t248 * t431 - t250 * t434 + t174) * t345 + (-t321 * t431 + t323 * t434 + t318) * t411 + t125 * qJD(4)) * t432 + (qJD(4) * t525 - t760) * t435) * t743 + ((t249 * t314 - t251 * t315) * t344 - (t248 * t314 - t250 * t315) * t345 + (-t314 * t321 + t315 * t323) * t411 + (t109 * t432 + t665 * t68) * qJD(4) + ((qJD(4) * t67 + t474) * t435 + t440) * t422) * t747 + ((-t249 * t316 - t251 * t317) * t344 - (-t248 * t316 - t250 * t317) * t345 + (t316 * t321 + t317 * t323) * t411 + (t110 * t432 + t671 * t69) * qJD(4) + ((qJD(4) * t70 + t474) * t435 + t440) * t423) * t750 + (t89 * t348 + t37 * t638 + (t37 * t187 + t89 * t647 + (qJD(1) * t90 + t47) * t634) * t423 + (qJD(1) * t332 * t89 + t37 * t185 + t46 * t634 + t90 * t647) * t422 - t89 * (t252 * t411 - t333 * t345 + t488) - t90 * (-t253 * t411 - t333 * t344 + t546) - ((-t185 * t89 + t187 * t90) * t432 + t444 * t435) * qJD(4) + (t602 + (qJD(1) * t185 + t100) * t423 + (qJD(1) * t649 + t101) * t422 + t252 * t344 + t253 * t345 - t640) * t73) * m(5) + (-(t144 * t330 - t683) * qJD(1) - (t143 * (-t330 * t422 - t331 * t423) + t516 * t543) * qJD(3) + t74 * t508 + t143 * t448 + t516 * t363 + (-t107 * t422 - t108 * t423 + (-t145 * t423 + t685) * qJD(1)) * t389) * m(4) - qJD(1) * ((t432 * t630 + t435 * t629) * qJD(1) + ((t422 * t641 - t423 * t642) * t435 + (t422 * t643 + t423 * t644) * t432) * qJD(3)) / 0.2e1 + t526 * t547 + qJD(1) * (t422 * t79 - t423 * t78 + (t164 * t422 + t165 * t423) * qJD(1)) / 0.2e1 + ((-t621 * t678 - t626) * t423 + (t454 + (t470 * t422 + (t677 - t763) * t423) * qJD(3)) * t422) * t575 + ((-t410 * t677 + t626) * t422 + (t454 + (-t763 * t423 + (t678 + t470) * t422) * qJD(3)) * t423) * t578 + (qJD(1) * t91 + (t767 * t422 ^ 2 - t766 * t786 + (t114 * t422 + t115 * t423) * qJD(1)) * t780 + t5 + t9) * t422 / 0.2e1 - (qJD(1) * t92 + t10 + (-t767 * t786 + t423 ^ 2 * t766 + (t112 * t422 + t113 * t423) * qJD(1)) * t780 + t6) * t423 / 0.2e1 - t43 * t618 / 0.2e1 + (qJD(1) * t525 - t27 * t423 + t28 * t422) * t742 + (qJD(1) * t527 - t18 * t423 + t19 * t422) * t744 + t528 * t746 + (qJD(1) * t531 - t25 * t423 + t26 * t422) * t748 + (qJD(1) * t529 - t23 * t423 + t24 * t422) * t749 - (t38 * t422 + t39 * t423) * t617 / 0.2e1 + (t482 + t60 + t38 + t35) * t625 / 0.2e1 + (qJD(1) * t535 - t16 * t423 + t17 * t422) * t752 + (qJD(1) * t533 - t14 * t423 + t15 * t422) * t753 + t530 * t755 + t532 * t756 + t534 * t757 + t536 * t758 - t312 * t35 / 0.2e1 - t313 * t36 / 0.2e1 + (t61 + t39 + t36 + t481) * t579 - t378 * t41 / 0.2e1; ((qJD(3) * t531 - t49) * t435 + (-qJD(1) * t532 + qJD(3) * t109 + t25 * t422 + t26 * t423) * t432) * t748 + t43 * t574 + t10 * t589 + t9 * t588 + (t472 * t435 + (-t431 * t759 + t434 * t450) * t432) * t743 + (t707 + t711 + t712 + t719 - t720) * t739 + ((qJD(3) * t529 - t48) * t435 + (-qJD(1) * t530 + qJD(3) * t110 + t23 * t422 + t24 * t423) * t432) * t749 + ((qJD(3) * t525 - t62) * t435 + (-qJD(1) * t526 + qJD(3) * t125 + t27 * t422 + t28 * t423) * t432) * t742 + (-t110 * t435 + t432 * t529) * t755 + (-t109 * t435 + t432 * t531) * t756 + t447 + (-t125 * t435 + t432 * t525) * t547 + (t316 * t759 + t450 * t317 - t423 * t467) * t750 + (-t314 * t759 + t315 * t450 - t422 * t467) * t747 + t784 * t39 + t783 * t38 + (-t57 * (t306 * t411 + t345 * t611 + t563) - t58 * (t307 * t411 + t344 * t611 + t564) - t50 * (-t306 * t344 + t307 * t345 + t565) + t57 * t556 + t58 * t656 + t50 * t706 + (-t22 * t191 + t57 * t103 - t21 * t652 - t58 * t709 + ((-t50 * t191 - t58 * t637) * t423 + t489 * t422) * qJD(3)) * t435 + ((-t57 * t653 + t715) * qJD(3) + (qJD(1) * t489 + t50 * t103 - t11 * t191 - t21 * t637 + t58 * t648) * t423 + (t22 * t287 + t57 * t240 + (t58 * t287 - t779) * qJD(1) + t764) * t422) * t432 + t490) * m(6) + ((qJD(3) * t444 - t90 * t100 + t89 * t101 + t47 * t185 - t46 * t187) * t435 + (t89 * (-qJD(3) * t185 + t233 * t422) + t90 * (qJD(3) * t187 - t233 * t423) + t37 * t511 + t73 * (-t100 * t422 + t101 * t423 - t185 * t625 - t187 * t624) + (t47 * t422 - t46 * t423 + (t422 * t90 + t716) * qJD(1)) * t332) * t432 - t89 * (-t224 * t411 - t345 * t353) - t90 * (t225 * t411 - t344 * t353) - t73 * (t224 * t344 + t225 * t345)) * m(5); t447 + (t21 * (-t163 * t435 - t299 * t667) - t11 * t163 * t673 + t490 + (-t435 * t87 + (-t201 * t432 - t299 * t619) * t423 + t656 - t564) * t58 + (-t161 * t620 + t556 - t563) * t57 + (-t163 * t597 + t706 + (-t163 * t619 + (-qJD(1) * t161 - t87) * t432) * t422 - t565) * t50) * m(6);];
tauc = t1(:);