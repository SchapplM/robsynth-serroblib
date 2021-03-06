% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:08
% EndTime: 2019-03-09 01:39:52
% DurationCPUTime: 37.44s
% Computational Cost: add. (35030->1044), mult. (28425->1355), div. (0->0), fcn. (26072->12), ass. (0->498)
t402 = pkin(10) + qJ(4);
t398 = cos(t402);
t403 = qJ(1) + pkin(9);
t399 = cos(t403);
t406 = cos(pkin(11));
t622 = t399 * t406;
t396 = sin(t403);
t404 = sin(pkin(11));
t628 = t396 * t404;
t303 = t398 * t628 + t622;
t623 = t399 * t404;
t627 = t396 * t406;
t304 = t398 * t627 - t623;
t395 = sin(t402);
t632 = t395 * t396;
t156 = Icges(6,5) * t304 - Icges(6,6) * t303 + Icges(6,3) * t632;
t629 = t396 * t398;
t642 = Icges(5,3) * t399;
t235 = Icges(5,5) * t629 - Icges(5,6) * t632 - t642;
t351 = Icges(5,4) * t632;
t649 = Icges(5,5) * t399;
t243 = Icges(5,1) * t629 - t351 - t649;
t645 = Icges(5,6) * t399;
t239 = Icges(5,4) * t629 - Icges(5,2) * t632 - t645;
t636 = t239 * t395;
t481 = -t243 * t398 + t636;
t159 = Icges(6,4) * t304 - Icges(6,2) * t303 + Icges(6,6) * t632;
t162 = Icges(6,1) * t304 - Icges(6,4) * t303 + Icges(6,5) * t632;
t485 = -t159 * t303 + t162 * t304;
t744 = t156 * t632 - t235 * t399 - t396 * t481 + t485;
t305 = -t398 * t623 + t627;
t306 = t398 * t622 + t628;
t631 = t395 * t399;
t158 = Icges(6,5) * t306 + Icges(6,6) * t305 + Icges(6,3) * t631;
t161 = Icges(6,4) * t306 + Icges(6,2) * t305 + Icges(6,6) * t631;
t164 = Icges(6,1) * t306 + Icges(6,4) * t305 + Icges(6,5) * t631;
t48 = t158 * t632 - t303 * t161 + t304 * t164;
t385 = Icges(5,4) * t398;
t493 = -Icges(5,2) * t395 + t385;
t240 = Icges(5,6) * t396 + t399 * t493;
t653 = Icges(5,4) * t395;
t327 = Icges(5,1) * t398 - t653;
t244 = Icges(5,5) * t396 + t327 * t399;
t209 = t244 * t629;
t323 = Icges(5,5) * t398 - Icges(5,6) * t395;
t236 = Icges(5,3) * t396 + t323 * t399;
t527 = t236 * t399 - t209;
t88 = -t240 * t632 - t527;
t743 = t48 + t88;
t490 = Icges(6,5) * t406 - Icges(6,6) * t404;
t264 = -Icges(6,3) * t398 + t395 * t490;
t492 = Icges(6,4) * t406 - Icges(6,2) * t404;
t266 = -Icges(6,6) * t398 + t395 * t492;
t495 = Icges(6,1) * t406 - Icges(6,4) * t404;
t268 = -Icges(6,5) * t398 + t395 * t495;
t324 = Icges(5,2) * t398 + t653;
t326 = Icges(5,1) * t395 + t385;
t477 = t324 * t395 - t326 * t398;
t322 = Icges(5,5) * t395 + Icges(5,6) * t398;
t633 = t322 * t399;
t742 = t264 * t632 - t266 * t303 + t268 * t304 - t396 * t477 - t633;
t634 = t322 * t396;
t741 = t264 * t631 + t266 * t305 + t268 * t306 - t399 * t477 + t634;
t718 = t305 * t159 + t306 * t162;
t49 = t156 * t631 + t718;
t50 = t158 * t631 + t305 * t161 + t306 * t164;
t452 = (t396 * t50 - t399 * t49) * qJD(4);
t740 = qJD(1) * t741 + t452;
t739 = (t743 * t396 - t399 * t744) * qJD(4);
t738 = t742 * qJD(1);
t401 = pkin(11) + qJ(6);
t394 = sin(t401);
t397 = cos(t401);
t626 = t397 * t399;
t273 = t394 * t629 + t626;
t624 = t399 * t394;
t274 = t397 * t629 - t624;
t137 = Icges(7,5) * t274 - Icges(7,6) * t273 + Icges(7,3) * t632;
t260 = Icges(7,4) * t274;
t140 = -Icges(7,2) * t273 + Icges(7,6) * t632 + t260;
t259 = Icges(7,4) * t273;
t144 = -Icges(7,1) * t274 - Icges(7,5) * t632 + t259;
t727 = t140 * t394 + t144 * t397;
t51 = -t137 * t398 - t395 * t727;
t410 = sin(qJ(1));
t681 = pkin(1) * t410;
t411 = cos(qJ(1));
t400 = t411 * pkin(1);
t737 = t738 + t739;
t625 = t398 * t399;
t615 = -t396 * t235 - t243 * t625;
t89 = -t239 * t631 - t615;
t614 = t396 * t236 + t244 * t625;
t90 = -t240 * t631 + t614;
t450 = (t396 * t90 - t399 * t89) * qJD(4);
t736 = t450 + t740;
t447 = qJD(4) * t324;
t170 = qJD(1) * t240 - t396 * t447;
t448 = qJD(4) * t326;
t173 = qJD(1) * t244 - t396 * t448;
t717 = t159 * t404 - t162 * t406;
t582 = qJD(4) * t395;
t553 = t404 * t582;
t216 = qJD(1) * t305 + t396 * t553;
t552 = t406 * t582;
t217 = qJD(1) * t306 - t396 * t552;
t580 = qJD(4) * t398;
t583 = qJD(1) * t399;
t281 = t395 * t583 + t396 * t580;
t94 = Icges(6,5) * t217 + Icges(6,6) * t216 + Icges(6,3) * t281;
t96 = Icges(6,4) * t217 + Icges(6,2) * t216 + Icges(6,6) * t281;
t98 = Icges(6,1) * t217 + Icges(6,4) * t216 + Icges(6,5) * t281;
t735 = (t94 - t170) * t398 + (t404 * t96 - t406 * t98 - t173) * t395 + (-t156 * t395 + t398 * t717 + t481) * qJD(4);
t169 = -t399 * t447 + (-t396 * t493 + t645) * qJD(1);
t172 = -t399 * t448 + (-t327 * t396 + t649) * qJD(1);
t635 = t240 * t395;
t480 = -t244 * t398 + t635;
t483 = -t161 * t404 + t164 * t406;
t214 = qJD(1) * t303 + t399 * t553;
t215 = -qJD(1) * t304 - t399 * t552;
t579 = qJD(4) * t399;
t551 = t398 * t579;
t585 = qJD(1) * t396;
t557 = t395 * t585;
t282 = t551 - t557;
t93 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t282;
t95 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t282;
t97 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t282;
t734 = (-t93 + t169) * t398 + (-t404 * t95 + t406 * t97 + t172) * t395 + (t158 * t395 + t398 * t483 - t480) * qJD(4);
t265 = Icges(6,3) * t395 + t398 * t490;
t230 = t265 * qJD(4);
t267 = Icges(6,6) * t395 + t398 * t492;
t231 = t267 * qJD(4);
t269 = Icges(6,5) * t395 + t398 * t495;
t232 = t269 * qJD(4);
t315 = t493 * qJD(4);
t316 = t327 * qJD(4);
t416 = qJD(1) * t322 - t315 * t395 + t316 * t398 + (-t324 * t398 - t326 * t395) * qJD(4);
t698 = t477 * qJD(1) + t323 * qJD(4);
t733 = t214 * t266 + t215 * t268 + t230 * t631 + t231 * t305 + t232 * t306 + t264 * t282 + t396 * t698 + t416 * t399;
t732 = t216 * t266 + t217 * t268 + t230 * t632 - t231 * t303 + t232 * t304 + t264 * t281 + t416 * t396 - t399 * t698;
t731 = t49 + t89;
t123 = t239 * t398 + t243 * t395;
t730 = -t156 * t398 - t395 * t717 + t123;
t124 = t240 * t398 + t244 * t395;
t729 = -t158 * t398 + t395 * t483 + t124;
t409 = -pkin(7) - qJ(3);
t373 = t399 * t409;
t728 = t681 - t373;
t576 = qJD(6) * t395;
t581 = qJD(4) * t396;
t312 = t399 * t576 + t581;
t313 = -t396 * t576 + t579;
t575 = qJD(6) * t398;
t368 = qJD(1) - t575;
t43 = t137 * t632 - t140 * t273 - t144 * t274;
t275 = t396 * t397 - t398 * t624;
t276 = t394 * t396 + t397 * t625;
t139 = Icges(7,5) * t276 + Icges(7,6) * t275 + Icges(7,3) * t631;
t652 = Icges(7,4) * t276;
t142 = Icges(7,2) * t275 + Icges(7,6) * t631 + t652;
t261 = Icges(7,4) * t275;
t145 = Icges(7,1) * t276 + Icges(7,5) * t631 + t261;
t44 = t139 * t632 - t273 * t142 + t274 * t145;
t489 = Icges(7,5) * t397 - Icges(7,6) * t394;
t233 = -Icges(7,3) * t398 + t395 * t489;
t650 = Icges(7,4) * t397;
t491 = -Icges(7,2) * t394 + t650;
t237 = -Icges(7,6) * t398 + t395 * t491;
t651 = Icges(7,4) * t394;
t494 = Icges(7,1) * t397 - t651;
t241 = -Icges(7,5) * t398 + t395 * t494;
t76 = t233 * t632 - t237 * t273 + t241 * t274;
t12 = t312 * t44 - t313 * t43 + t368 * t76;
t45 = t137 * t631 + t275 * t140 - t144 * t276;
t46 = t139 * t631 + t275 * t142 + t276 * t145;
t77 = t233 * t631 + t237 * t275 + t241 * t276;
t13 = t312 * t46 - t313 * t45 + t77 * t368;
t412 = qJD(1) ^ 2;
t670 = rSges(4,2) * sin(pkin(10));
t407 = cos(pkin(10));
t672 = rSges(4,1) * t407;
t469 = t396 * rSges(4,3) + (-t670 + t672) * t399;
t380 = t396 * qJ(3);
t334 = t399 * pkin(2) + t380;
t534 = t334 + t400;
t722 = t469 + t534;
t408 = -pkin(8) - qJ(5);
t621 = qJ(5) + t408;
t530 = t621 * t398;
t392 = pkin(5) * t406 + pkin(4);
t674 = pkin(4) - t392;
t547 = t674 * t395;
t225 = t530 - t547;
t505 = rSges(7,1) * t274 - rSges(7,2) * t273;
t146 = rSges(7,3) * t632 + t505;
t504 = rSges(7,1) * t397 - rSges(7,2) * t394;
t253 = -rSges(7,3) * t398 + t395 * t504;
t328 = pkin(4) * t395 - qJ(5) * t398;
t529 = (-t225 - t328) * t399;
t577 = qJD(5) * t399;
t346 = t395 * t577;
t378 = qJD(3) * t396;
t595 = t346 + t378;
t720 = qJD(4) * t529 - t146 * t368 - t253 * t313 + t595;
t478 = -t266 * t404 + t268 * t406;
t719 = qJD(4) * (-t396 * t483 - t399 * t717) + (t265 - t478) * qJD(1);
t716 = 0.2e1 * qJD(4);
t715 = -rSges(6,3) - qJ(5);
t713 = t395 * t621;
t712 = t398 * t674;
t507 = rSges(6,1) * t406 - rSges(6,2) * t404;
t270 = -rSges(6,3) * t398 + t395 * t507;
t377 = qJD(5) * t395;
t711 = (qJD(4) * t270 - t377) * t396;
t710 = (qJD(4) * t225 - t377) * t396;
t381 = t399 * qJ(3);
t330 = pkin(2) * t396 - t381;
t393 = pkin(3) * t407 + pkin(2);
t352 = t396 * t393;
t594 = -t352 - t373;
t226 = t330 + t594;
t222 = qJD(1) * t226;
t320 = qJD(1) * t330;
t709 = t222 - t320;
t708 = -rSges(5,2) * t632 - t399 * rSges(5,3);
t531 = t399 * rSges(3,1) - rSges(3,2) * t396;
t707 = t400 + t531;
t584 = qJD(1) * t398;
t597 = t326 + t493;
t598 = -t324 + t327;
t706 = t719 * t395 + t264 * t584 + (-t395 * t597 + t398 * t598) * qJD(1);
t439 = t239 * t399 - t240 * t396;
t696 = t396 * (-t324 * t399 + t244) - t399 * (-Icges(5,2) * t629 + t243 - t351);
t705 = -t395 * t696 + (-t156 * t399 + t158 * t396 + t439) * t398;
t639 = qJ(5) * t395;
t680 = pkin(4) * t398;
t332 = t639 + t680;
t299 = t332 * t396;
t278 = qJD(1) * t299;
t523 = -qJD(6) + t584;
t554 = t395 * t579;
t701 = t396 * t523 + t554;
t446 = qJD(4) * t322;
t700 = -t399 * t446 + (-t323 * t396 + t480 + t642) * qJD(1);
t587 = qJD(1) * t236;
t699 = qJD(1) * t481 - t396 * t446 + t587;
t234 = Icges(7,3) * t395 + t398 * t489;
t482 = -t237 * t394 + t241 * t397;
t487 = -t142 * t394 + t145 * t397;
t695 = t312 * (-t233 * t399 - t487) - t313 * (-t233 * t396 + t727) + t368 * (t234 - t482);
t288 = (-Icges(7,2) * t397 - t651) * t395;
t694 = t312 * (-Icges(7,2) * t276 + t145 + t261) - t313 * (-Icges(7,2) * t274 - t144 - t259) + t368 * (t241 + t288);
t693 = m(6) / 0.2e1;
t692 = m(7) / 0.2e1;
t571 = qJD(4) * qJD(6);
t544 = t398 * t571;
t223 = qJD(1) * t312 + t396 * t544;
t691 = t223 / 0.2e1;
t224 = qJD(1) * t313 + t399 * t544;
t690 = t224 / 0.2e1;
t689 = -t312 / 0.2e1;
t688 = t312 / 0.2e1;
t687 = -t313 / 0.2e1;
t686 = t313 / 0.2e1;
t685 = -t368 / 0.2e1;
t684 = t368 / 0.2e1;
t683 = t396 / 0.2e1;
t682 = -t399 / 0.2e1;
t472 = t368 * t397;
t555 = t395 * t581;
t121 = t396 * t472 + (-t399 * t523 + t555) * t394;
t473 = t368 * t394;
t122 = t523 * t626 + (-t397 * t582 + t473) * t396;
t68 = Icges(7,5) * t122 + Icges(7,6) * t121 + Icges(7,3) * t281;
t70 = Icges(7,4) * t122 + Icges(7,2) * t121 + Icges(7,6) * t281;
t72 = Icges(7,1) * t122 + Icges(7,4) * t121 + Icges(7,5) * t281;
t8 = (-qJD(4) * t727 - t68) * t398 + (qJD(4) * t137 - t394 * t70 + t397 * t72 + (-t140 * t397 + t144 * t394) * qJD(6)) * t395;
t679 = t8 * t313;
t119 = t394 * t701 + t399 * t472;
t120 = -t397 * t701 + t399 * t473;
t67 = Icges(7,5) * t120 + Icges(7,6) * t119 + Icges(7,3) * t282;
t69 = Icges(7,4) * t120 + Icges(7,2) * t119 + Icges(7,6) * t282;
t71 = Icges(7,1) * t120 + Icges(7,4) * t119 + Icges(7,5) * t282;
t9 = (qJD(4) * t487 - t67) * t398 + (qJD(4) * t139 - t394 * t69 + t397 * t71 + (-t142 * t397 - t145 * t394) * qJD(6)) * t395;
t678 = t9 * t312;
t676 = qJD(1) / 0.2e1;
t675 = pkin(2) - t393;
t285 = (-Icges(7,5) * t394 - Icges(7,6) * t397) * t395;
t165 = qJD(4) * t234 + qJD(6) * t285;
t238 = Icges(7,6) * t395 + t398 * t491;
t168 = qJD(4) * t238 + qJD(6) * t288;
t242 = Icges(7,5) * t395 + t398 * t494;
t291 = (-Icges(7,1) * t394 - t650) * t395;
t171 = qJD(4) * t242 + qJD(6) * t291;
t27 = (qJD(4) * t482 - t165) * t398 + (qJD(4) * t233 - t168 * t394 + t171 * t397 + (-t237 * t397 - t241 * t394) * qJD(6)) * t395;
t545 = t395 * t571;
t86 = -t233 * t398 + t395 * t482;
t673 = t27 * t368 + t86 * t545;
t669 = rSges(6,3) * t395;
t667 = rSges(7,3) * t395;
t148 = t276 * rSges(7,1) + t275 * rSges(7,2) + rSges(7,3) * t631;
t255 = t398 * t504 + t667;
t296 = (-rSges(7,1) * t394 - rSges(7,2) * t397) * t395;
t177 = qJD(4) * t255 + qJD(6) * t296;
t338 = qJ(5) * t551;
t443 = -t396 * t584 - t554;
t152 = pkin(4) * t443 - qJ(5) * t557 + t338 + t346;
t374 = qJ(3) * t583;
t569 = t412 * t681;
t573 = qJD(1) * qJD(3);
t589 = t374 + t378;
t476 = qJD(1) * (-pkin(2) * t585 + t589) + t396 * t573 - t569;
t456 = qJD(1) * (-t374 + (t396 * t675 - t373) * qJD(1)) + t476;
t572 = qJD(4) * qJD(5);
t546 = t398 * t572;
t437 = t396 * t546 + t456 + (t152 + t346) * qJD(1);
t227 = -t712 - t713;
t578 = qJD(5) * t398;
t294 = qJD(4) * t332 - t578;
t613 = -t227 * qJD(4) - t294;
t564 = t120 * rSges(7,1) + t119 * rSges(7,2) + rSges(7,3) * t551;
t73 = -rSges(7,3) * t557 + t564;
t570 = pkin(5) * t623;
t599 = qJD(1) * t570 + t408 * t557;
t81 = -t338 + (-t398 * t408 + t547) * t579 + (t639 + t712) * t585 + t599;
t10 = qJD(1) * t81 - t177 * t312 - t224 * t253 + t368 * t73 + (qJD(1) * t529 + t148 * t576 + t396 * t613) * qJD(4) + t437;
t665 = t10 * t399;
t310 = t328 * t581;
t568 = t412 * t400;
t521 = t399 * t573 - t568;
t455 = qJD(1) * t310 + t399 * t546 + t521;
t342 = pkin(4) * t555;
t358 = pkin(4) * t625;
t153 = qJ(5) * t281 + qJD(1) * t358 + t377 * t396 - t342;
t379 = qJD(3) * t399;
t295 = qJD(1) * t334 - t379;
t362 = t409 * t585;
t612 = t362 - (-t399 * t675 - t380) * qJD(1) - t295;
t562 = -t153 + t612;
t506 = rSges(7,1) * t122 + rSges(7,2) * t121;
t74 = rSges(7,3) * t281 + t506;
t367 = pkin(5) * t628;
t82 = t342 + (-t392 * t395 - t530) * t581 + (t227 * t399 + t367) * qJD(1);
t11 = -t177 * t313 + t223 * t253 - t368 * t74 + (-t146 * t576 + t399 * t613) * qJD(4) + (t562 - t82 + t710) * qJD(1) + t455;
t664 = t11 * t396;
t329 = rSges(5,1) * t395 + rSges(5,2) * t398;
t301 = t329 * t399;
t386 = t396 * rSges(5,3);
t256 = rSges(5,1) * t625 - rSges(5,2) * t631 + t386;
t353 = t399 * t393;
t525 = -t396 * t409 + t353;
t518 = t525 - t334 + t534;
t556 = t329 * t581;
t92 = -t556 - t379 + (t256 + t518) * qJD(1);
t663 = t301 * t92;
t566 = rSges(5,1) * t629;
t254 = t566 + t708;
t549 = t329 * t579;
t513 = t378 - t549;
t535 = -t330 - t681;
t519 = t226 + t535;
t91 = (-t254 + t519) * qJD(1) + t513;
t659 = t396 * t91;
t657 = t51 * t223;
t52 = -t139 * t398 + t395 * t487;
t656 = t52 * t224;
t655 = -rSges(7,3) + t408;
t630 = t395 * t408;
t520 = -t392 * t629 + t570;
t154 = (t680 + t713) * t396 + t520;
t618 = t146 - t154;
t302 = qJ(5) * t631 + t358;
t468 = t392 * t625 - t399 * t630 + t367;
t155 = t468 - t302;
t617 = -t155 - t302;
t176 = t306 * rSges(6,1) + t305 * rSges(6,2) + rSges(6,3) * t631;
t616 = -t176 - t302;
t611 = t225 + t253;
t610 = -t227 - t332;
t272 = t398 * t507 + t669;
t606 = -t272 * qJD(4) - t294;
t604 = t396 * t299 + t399 * t302;
t603 = -t270 - t328;
t602 = -t272 - t332;
t300 = t328 * t399;
t601 = -qJD(1) * t300 + t396 * t578;
t600 = -t310 - t379;
t596 = rSges(5,2) * t557 + rSges(5,3) * t583;
t593 = t353 + t400;
t363 = t396 * t670;
t592 = rSges(4,3) * t583 + qJD(1) * t363;
t591 = t362 + t379;
t590 = t399 * rSges(4,3) + t363;
t588 = t378 - t320;
t586 = qJD(1) * t323;
t567 = t396 * t672;
t563 = t399 * t152 + t396 * t153 + t299 * t583;
t561 = -t177 + t613;
t560 = t215 * rSges(6,1) + t214 * rSges(6,2) + rSges(6,3) * t551;
t559 = -t328 - t611;
t297 = t328 * t396;
t558 = -t297 * t581 - t300 * t579 + t377;
t548 = -pkin(2) - t672;
t542 = t583 / 0.2e1;
t541 = t582 / 0.2e1;
t540 = -t581 / 0.2e1;
t539 = t581 / 0.2e1;
t537 = t579 / 0.2e1;
t474 = t302 + t518;
t56 = -t711 + (t176 + t474) * qJD(1) + t600;
t533 = t56 * t603;
t528 = t603 * t399;
t526 = -t235 + t635;
t524 = qJD(1) * t297 + t398 * t577;
t522 = t153 * t581 + t395 * t572 + (t152 + t278) * t579;
t517 = t594 - t681;
t514 = qJD(6) * t541;
t331 = rSges(3,1) * t396 + rSges(3,2) * t399;
t510 = rSges(5,1) * t398 - rSges(5,2) * t395;
t509 = rSges(6,1) * t217 + rSges(6,2) * t216;
t508 = rSges(6,1) * t304 - rSges(6,2) * t303;
t503 = t396 * t44 - t399 * t43;
t502 = t396 * t43 + t399 * t44;
t501 = t396 * t46 - t399 * t45;
t500 = t396 * t45 + t399 * t46;
t499 = t396 * t52 - t399 * t51;
t498 = t396 * t51 + t399 * t52;
t497 = -t396 * t92 - t399 * t91;
t486 = t146 * t399 - t148 * t396;
t479 = t254 * t396 + t256 * t399;
t475 = -t299 + t519;
t467 = -t566 - t681;
t466 = t299 * t581 + t302 * t579 + qJD(2) - t578;
t465 = -t392 * t398 - t393 - t667;
t298 = t329 * t396;
t444 = -t332 - t669;
t442 = -t137 * t313 + t139 * t312 + t233 * t368;
t441 = (-Icges(7,5) * t273 - Icges(7,6) * t274) * t313 - (Icges(7,5) * t275 - Icges(7,6) * t276) * t312 - t285 * t368;
t174 = rSges(6,3) * t632 + t508;
t438 = t395 * t441;
t426 = (Icges(7,1) * t275 - t142 - t652) * t312 - (-Icges(7,1) * t273 - t140 - t260) * t313 + (-t237 + t291) * t368;
t178 = rSges(5,1) * t443 - rSges(5,2) * t551 + t596;
t179 = -qJD(4) * t298 + (t399 * t510 + t386) * qJD(1);
t424 = t178 * t399 + t179 * t396 + (t254 * t399 - t256 * t396) * qJD(1);
t34 = t146 * t312 + t148 * t313 + (-t154 * t396 + t155 * t399) * qJD(4) + t466;
t39 = (t154 + t475) * qJD(1) + t720;
t40 = t148 * t368 - t253 * t312 - t710 + (t155 + t474) * qJD(1) + t600;
t419 = t34 * t486 + (t39 * t396 - t399 * t40) * t253;
t418 = qJD(1) * t235 - qJD(4) * t123 - t170 * t395 + t173 * t398;
t417 = -qJD(4) * t124 - t169 * t395 + t172 * t398 + t587;
t414 = t695 * t395;
t317 = t510 * qJD(4);
t311 = t328 * t585;
t280 = (t396 ^ 2 + t399 ^ 2) * t582;
t271 = t567 - t590;
t213 = t270 * t399;
t212 = t270 * t396;
t207 = t268 * t399;
t206 = t268 * t396;
t205 = t266 * t399;
t204 = t266 * t396;
t199 = t253 * t399;
t198 = t253 * t396;
t197 = t241 * t399;
t196 = t241 * t396;
t195 = t237 * t399;
t194 = t237 * t396;
t191 = t225 * t399;
t190 = t225 * t396;
t189 = rSges(7,1) * t275 - rSges(7,2) * t276;
t188 = -rSges(7,1) * t273 - rSges(7,2) * t274;
t181 = qJD(1) * t722 - t379;
t180 = t378 + (-t271 + t535) * qJD(1);
t126 = (-qJD(1) * t469 - t295) * qJD(1) + t521;
t125 = qJD(1) * (-qJD(1) * t567 + t592) + t476;
t116 = qJD(4) * t479 + qJD(2);
t100 = rSges(6,3) * t281 + t509;
t99 = -rSges(6,3) * t557 + t560;
t66 = -t317 * t579 + (-t179 + t556 + t612) * qJD(1) + t521;
t65 = -t317 * t581 + (t178 - t549) * qJD(1) + t456;
t62 = (t174 * t396 + t176 * t399) * qJD(4) + t466;
t57 = t424 * qJD(4);
t55 = qJD(4) * t528 + (-t174 + t475) * qJD(1) + t595;
t33 = t606 * t579 + (-t100 + t562 + t711) * qJD(1) + t455;
t32 = qJD(1) * t99 + (qJD(1) * t528 + t396 * t606) * qJD(4) + t437;
t25 = t121 * t237 + t122 * t241 + t165 * t632 - t168 * t273 + t171 * t274 + t233 * t281;
t24 = t119 * t237 + t120 * t241 + t165 * t631 + t168 * t275 + t171 * t276 + t233 * t282;
t21 = (t100 * t396 + t399 * t99 + (t174 * t399 + t396 * t616) * qJD(1)) * qJD(4) + t522;
t20 = t312 * t52 - t313 * t51 + t368 * t86;
t7 = t121 * t142 + t122 * t145 + t139 * t281 - t273 * t69 + t274 * t71 + t632 * t67;
t6 = t121 * t140 - t122 * t144 + t137 * t281 - t273 * t70 + t274 * t72 + t632 * t68;
t5 = t119 * t142 + t120 * t145 + t139 * t282 + t275 * t69 + t276 * t71 + t631 * t67;
t4 = t119 * t140 - t120 * t144 + t137 * t282 + t275 * t70 + t276 * t72 + t631 * t68;
t3 = t146 * t224 - t148 * t223 + t312 * t74 + t313 * t73 + (t396 * t82 + t399 * t81 + (-t154 * t399 + t396 * t617) * qJD(1)) * qJD(4) + t522;
t2 = t223 * t43 + t224 * t44 + t25 * t368 + t312 * t7 - t313 * t6 + t545 * t76;
t1 = t223 * t45 + t224 * t46 + t24 * t368 + t312 * t5 - t313 * t4 + t545 * t77;
t14 = [m(3) * ((-t331 * t412 - t569) * t707 + (-t568 + (-0.2e1 * t531 - t400 + t707) * t412) * (-t331 - t681)) + t656 / 0.2e1 + t657 / 0.2e1 + t678 / 0.2e1 + t25 * t687 + t24 * t688 + t77 * t690 + t76 * t691 - t679 / 0.2e1 + t673 + (((t88 - t209 + (t236 + t636) * t399 + t615) * t399 + t614 * t396) * qJD(4) + t740) * t537 + (t686 + t687) * t13 + ((t315 - t230) * t398 + (-t231 * t404 + t232 * t406 + t316) * t395 + (t264 * t395 + t398 * t478 - t477) * qJD(4)) * qJD(1) + (t11 * (-t505 + t517 + t520) + t39 * (-t506 + t591) + t10 * (t468 + t148 + t593) + t40 * (-t392 * t554 - t408 * t551 + t564 + t595 + t599) + (-t10 * t409 + t39 * t655 * t580 + (t11 * t655 + t39 * (qJD(4) * t392 - qJD(5))) * t395) * t396 + ((-t39 * t411 - t40 * t410) * pkin(1) + (-t39 * pkin(5) * t404 + t40 * t465) * t396 + (t39 * (t465 + t630) - t40 * t409) * t399) * qJD(1) - (-t278 - t39 + (t154 - t681) * qJD(1) + t709 + t720) * t40) * m(7) + (t33 * (-t508 + t517) + t32 * (t593 - t616) + (-t32 * t409 + t33 * t444) * t396 - t533 * t579 + (t342 - t509 + t591 + (t580 * t715 - t377) * t396 + (-t400 + (-t393 + t444) * t399) * qJD(1)) * t55 + (-pkin(4) * t554 - t222 + t278 + t338 - t346 + t55 + t560 - t588 + t595 + (-t681 + (t395 * t715 - t393 - t680) * t396 + t174 + t728) * qJD(1)) * t56) * m(6) + (t66 * (t467 + t594 - t708) + t65 * (t256 + t400 + t525) + (t329 * t659 - t663) * qJD(4) + (t591 + (-t386 - t400 + (-t393 - t510) * t399) * qJD(1)) * t91 + (t91 - t513 - t709 + t378 + t596 + (t254 - t352 + t467 + t728) * qJD(1)) * t92) * m(5) + (t126 * (t396 * t548 + t381 + t590 - t681) + t180 * t379 + t125 * t722 + t181 * (t589 + t592) + ((-t180 * t411 - t181 * t410) * pkin(1) + t180 * (t548 + t670) * t399 + (t180 * (-rSges(4,3) - qJ(3)) + t181 * t548) * t396) * qJD(1) - (-t180 + (-t271 - t681) * qJD(1) + t588) * t181) * m(4) + (t733 + t734) * t539 + (((t399 * t526 + t485 - t614 + t90) * t399 + (t396 * t526 - t48 + t527 - t718 + t731) * t396) * qJD(4) + t737 - t738) * t540 - (t732 - t735 + t736) * t579 / 0.2e1 + ((t730 + t742) * t396 + (t729 + t741) * t399) * qJD(4) * t676; m(5) * t57 + m(6) * t21 + m(7) * t3; 0.2e1 * (-t665 / 0.2e1 + t664 / 0.2e1) * m(7) + 0.2e1 * (t32 * t682 + t33 * t683) * m(6) + 0.2e1 * (t65 * t682 + t66 * t683) * m(5) + 0.2e1 * (t125 * t682 + t126 * t683) * m(4); t499 * t514 - t20 * t576 / 0.2e1 + (qJD(1) * t502 + t396 * t7 - t399 * t6) * t687 + (qJD(1) * t500 + t396 * t5 - t399 * t4) * t688 + ((-t195 * t275 - t197 * t276) * t312 - (-t194 * t275 - t196 * t276) * t313 + (t238 * t275 + t242 * t276) * t368 + (t395 * t77 + t45 * t629) * qJD(6) + ((qJD(6) * t46 + t442) * t398 + t414) * t399) * t689 + t501 * t690 + t503 * t691 + (qJD(1) * t498 + t396 * t9 - t399 * t8) * t684 + (((t195 * t394 - t197 * t397 + t139) * t312 - (t194 * t394 - t196 * t397 + t137) * t313 + (-t238 * t394 + t242 * t397 + t233) * t368 + t86 * qJD(6)) * t395 + (qJD(6) * t498 - t695) * t398) * t685 + ((t195 * t273 - t197 * t274) * t312 - (t194 * t273 - t196 * t274) * t313 + (-t238 * t273 + t242 * t274) * t368 + (t395 * t76 + t44 * t625) * qJD(6) + ((qJD(6) * t43 + t442) * t398 + t414) * t396) * t686 - (-t719 * t398 + ((-t267 * t404 + t269 * t406 + t264) * qJD(1) + ((t205 * t404 - t207 * t406 + t158) * t396 - (t204 * t404 - t206 * t406 + t156) * t399) * qJD(4)) * t395 + (t395 * t598 + t398 * t597) * qJD(1) + (t439 * t395 + t398 * t696) * qJD(4)) * qJD(1) / 0.2e1 + (t735 * t399 + t734 * t396 + (t396 * t730 + t399 * t729) * qJD(1)) * t676 + ((-t205 * t305 - t207 * t306) * t581 + (t267 * t305 + t269 * t306) * qJD(1) + (-t581 * t633 + t586) * t396 + ((t204 * t305 + t206 * t306 + t396 * t634 + t705) * qJD(4) + t706) * t399) * t540 + (-(t204 * t303 - t206 * t304) * t579 + (-t267 * t303 + t269 * t304) * qJD(1) + (-t579 * t634 - t586) * t399 + ((t205 * t303 - t207 * t304 + t399 * t633 + t705) * qJD(4) + t706) * t396) * t537 - (t12 * t396 + t13 * t399) * t575 / 0.2e1 + (t39 * t311 + t3 * t604 + t34 * t563 + (t11 * t559 + t39 * t561 + t3 * (t148 + t155) + t34 * (t73 + t81) + (t34 * t618 + t40 * t559) * qJD(1)) * t399 + (t10 * t559 + t40 * t561 + t3 * t618 + t34 * (t74 + t82) + (t39 * t611 + t34 * (-t148 + t617)) * qJD(1)) * t396 - t39 * (qJD(1) * t190 + t198 * t368 - t255 * t313 + t524) - t40 * (-qJD(1) * t191 - t199 * t368 - t255 * t312 + t601) - t34 * (-t198 * t312 - t199 * t313 + t558) - ((-t146 * t39 + t148 * t40) * t395 + t419 * t398) * qJD(6) - ((-t34 * t191 + t39 * t610) * t399 + (-t34 * t190 + t40 * t610) * t396) * qJD(4)) * m(7) + (t55 * t311 + t21 * t604 + t62 * t563 + (t33 * t603 + t55 * t606 + t21 * t176 + t62 * t99 + (t62 * t174 + t533) * qJD(1)) * t399 + (t32 * t603 + t56 * t606 + t21 * t174 + t62 * t100 + (t55 * t270 + t616 * t62) * qJD(1)) * t396 - t55 * (qJD(1) * t212 + t524) - t56 * (-qJD(1) * t213 + t601) - t62 * t558 - ((-t62 * t213 + t55 * t602) * t399 + (-t62 * t212 + t56 * t602) * t396) * qJD(4)) * m(6) + (-(t298 * t91 - t663) * qJD(1) - (t116 * (-t298 * t396 - t301 * t399) + t497 * t510) * qJD(4) + t57 * t479 + t116 * t424 + t497 * t317 + (-t65 * t396 - t66 * t399 + (-t399 * t92 + t659) * qJD(1)) * t329) * m(5) + (t1 + t733 * qJD(1) + ((-t156 * t282 - t159 * t214 - t162 * t215 - t305 * t96 - t306 * t98 - t418 * t399 - t631 * t94) * t399 + (t158 * t282 + t161 * t214 + t164 * t215 + t305 * t95 + t306 * t97 + t396 * t700 + t631 * t93 + (t417 - t699) * t399) * t396 + ((t50 + t90) * t399 + t731 * t396) * qJD(1)) * t716) * t683 + (t2 + t732 * qJD(1) + ((-t156 * t281 - t159 * t216 - t162 * t217 + t303 * t96 - t304 * t98 + t399 * t699 - t632 * t94) * t399 + (t158 * t281 + t161 * t216 + t164 * t217 - t303 * t95 + t304 * t97 + t417 * t396 + t632 * t93 + (-t418 - t700) * t399) * t396 + (t396 * t744 + t743 * t399) * qJD(1)) * t716) * t682 + (t12 + t737 + t739) * t585 / 0.2e1 + (t452 + t450 + t13 + t736) * t542; -m(6) * (t280 * t62 + t281 * t56 + t282 * t55) - m(7) * (t280 * t34 + t281 * t40 + t282 * t39) + 0.2e1 * ((t55 * t579 + t56 * t581 - t21) * t693 + (t39 * t579 + t40 * t581 - t3) * t692) * t398 + 0.2e1 * ((qJD(4) * t62 + t32 * t396 + t33 * t399 - t55 * t585 + t56 * t583) * t693 + (qJD(4) * t34 + t10 * t396 + t11 * t399 - t39 * t585 + t40 * t583) * t692) * t395; t1 * t631 / 0.2e1 + (t395 * t500 - t398 * t77) * t690 + ((qJD(4) * t500 - t24) * t398 + (-qJD(1) * t501 + qJD(4) * t77 + t396 * t4 + t399 * t5) * t395) * t688 + t2 * t632 / 0.2e1 + (t395 * t502 - t398 * t76) * t691 + ((qJD(4) * t502 - t25) * t398 + (-qJD(1) * t503 + qJD(4) * t76 + t396 * t6 + t399 * t7) * t395) * t687 + t20 * t541 - t398 * (t656 + t657 + t673 + t678 - t679) / 0.2e1 + (t395 * t498 - t398 * t86) * t514 + ((qJD(4) * t498 - t27) * t398 + (-qJD(1) * t499 + qJD(4) * t86 + t396 * t8 + t399 * t9) * t395) * t684 + (t275 * t694 + t426 * t276 - t399 * t438) * t689 + (-t273 * t694 + t274 * t426 - t396 * t438) * t686 + (t441 * t398 + (-t394 * t694 + t397 * t426) * t395) * t685 + (-t557 / 0.2e1 + t398 * t537) * t13 + (t395 * t542 + t398 * t539) * t12 + ((qJD(4) * t419 - t10 * t148 + t11 * t146 + t39 * t74 - t40 * t73) * t398 + (t39 * (-qJD(4) * t146 + t177 * t396) + t40 * (qJD(4) * t148 - t177 * t399) + t3 * t486 + t34 * (-t146 * t585 - t148 * t583 - t396 * t73 + t399 * t74) + (-t665 + t664 + (t39 * t399 + t396 * t40) * qJD(1)) * t253) * t395 - t39 * (-t188 * t368 - t296 * t313) - t40 * (t189 * t368 - t296 * t312) - t34 * (t188 * t312 + t189 * t313)) * m(7);];
tauc  = t14(:);
