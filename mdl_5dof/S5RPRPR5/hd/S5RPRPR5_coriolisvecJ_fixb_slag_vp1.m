% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:26
% EndTime: 2020-01-03 11:42:52
% DurationCPUTime: 45.25s
% Computational Cost: add. (22119->942), mult. (29651->1230), div. (0->0), fcn. (28321->10), ass. (0->419)
t726 = -Icges(5,3) - Icges(4,3);
t428 = qJ(3) + pkin(9);
t418 = sin(t428);
t419 = cos(t428);
t435 = cos(qJ(1));
t430 = cos(pkin(8));
t433 = sin(qJ(1));
t610 = t430 * t433;
t319 = -t418 * t610 - t419 * t435;
t602 = t435 * t418;
t606 = t433 * t419;
t320 = t430 * t606 - t602;
t432 = sin(qJ(3));
t601 = t435 * t432;
t434 = cos(qJ(3));
t605 = t433 * t434;
t353 = t430 * t605 - t601;
t604 = t434 * t435;
t608 = t432 * t433;
t466 = t430 * t608 + t604;
t429 = sin(pkin(8));
t612 = t429 * t433;
t713 = Icges(4,5) * t353 + Icges(5,5) * t320 - Icges(4,6) * t466 + Icges(5,6) * t319 - t726 * t612;
t321 = t430 * t602 - t606;
t609 = t430 * t435;
t322 = t418 * t433 + t419 * t609;
t354 = t430 * t601 - t605;
t355 = t430 * t604 + t608;
t611 = t429 * t435;
t681 = Icges(4,5) * t355 + Icges(5,5) * t322 - Icges(4,6) * t354 - Icges(5,6) * t321 - t726 * t611;
t712 = t726 * t430 + (Icges(4,5) * t434 + Icges(5,5) * t419 - Icges(4,6) * t432 - Icges(5,6) * t418) * t429;
t665 = (Icges(4,5) * t432 + Icges(5,5) * t418 + Icges(4,6) * t434 + Icges(5,6) * t419) * t429;
t308 = Icges(5,4) * t322;
t187 = Icges(5,2) * t321 - Icges(5,6) * t611 - t308;
t307 = Icges(5,4) * t321;
t189 = Icges(5,1) * t322 + Icges(5,5) * t611 - t307;
t339 = Icges(4,4) * t355;
t238 = Icges(4,2) * t354 - Icges(4,6) * t611 - t339;
t338 = Icges(4,4) * t354;
t240 = Icges(4,1) * t355 + Icges(4,5) * t611 - t338;
t724 = t187 * t321 + t189 * t322 + t238 * t354 + t240 * t355;
t626 = Icges(5,4) * t320;
t185 = Icges(5,2) * t319 + Icges(5,6) * t612 + t626;
t306 = Icges(5,4) * t319;
t188 = Icges(5,1) * t320 + Icges(5,5) * t612 + t306;
t629 = Icges(4,4) * t353;
t236 = -Icges(4,2) * t466 + Icges(4,6) * t612 + t629;
t337 = Icges(4,4) * t466;
t239 = Icges(4,1) * t353 + Icges(4,5) * t612 - t337;
t677 = -t321 * t185 + t322 * t188 - t354 * t236 + t355 * t239;
t723 = -t319 * t187 + t189 * t320 + t238 * t466 + t240 * t353;
t721 = t319 * t185 + t320 * t188 - t236 * t466 + t353 * t239 + t713 * t612;
t720 = -t681 * t612 - t723;
t704 = -t713 * t611 - t677;
t703 = t681 * t611 + t724;
t550 = qJD(3) * t430;
t403 = qJD(1) - t550;
t624 = Icges(5,4) * t419;
t280 = -Icges(5,6) * t430 + (-Icges(5,2) * t418 + t624) * t429;
t625 = Icges(5,4) * t418;
t281 = -Icges(5,5) * t430 + (Icges(5,1) * t419 - t625) * t429;
t627 = Icges(4,4) * t434;
t312 = -Icges(4,6) * t430 + (-Icges(4,2) * t432 + t627) * t429;
t628 = Icges(4,4) * t432;
t313 = -Icges(4,5) * t430 + (Icges(4,1) * t434 - t628) * t429;
t680 = -t280 * t321 + t281 * t322 - t312 * t354 + t313 * t355 + t712 * t611;
t722 = t680 * t403;
t687 = (-t187 * t418 - t189 * t419 - t238 * t432 - t240 * t434) * t429 + t681 * t430;
t719 = t280 * t319 + t281 * t320 - t312 * t466 + t313 * t353 + t612 * t712;
t220 = qJD(1) * t319 + qJD(3) * t322;
t221 = qJD(1) * t320 + qJD(3) * t321;
t253 = -qJD(1) * t466 + qJD(3) * t355;
t254 = qJD(1) * t353 + qJD(3) * t354;
t554 = qJD(1) * t433;
t530 = t429 * t554;
t717 = Icges(4,5) * t254 + Icges(5,5) * t221 + Icges(4,6) * t253 + Icges(5,6) * t220 - t530 * t726;
t222 = -qJD(1) * t321 - qJD(3) * t320;
t223 = qJD(1) * t322 + qJD(3) * t319;
t255 = -qJD(1) * t354 - qJD(3) * t353;
t456 = t466 * qJD(3);
t256 = qJD(1) * t355 - t456;
t553 = qJD(1) * t435;
t532 = t429 * t553;
t714 = Icges(4,5) * t256 + Icges(5,5) * t223 + Icges(4,6) * t255 + Icges(5,6) * t222 - t532 * t726;
t711 = t665 * qJD(3);
t715 = t717 * t435;
t710 = (t704 * t433 - t703 * t435) * t429;
t709 = (t721 * t433 - t720 * t435) * t429;
t620 = qJ(2) * t433;
t377 = pkin(1) * t435 + t620;
t559 = pkin(1) * t553 + qJ(2) * t554;
t708 = -qJD(1) * t377 + t559;
t420 = qJ(5) + t428;
t408 = sin(t420);
t603 = t435 * t408;
t409 = cos(t420);
t607 = t433 * t409;
t303 = t430 * t603 - t607;
t304 = t408 * t433 + t409 * t609;
t165 = Icges(6,5) * t304 - Icges(6,6) * t303 + Icges(6,3) * t611;
t290 = Icges(6,4) * t304;
t169 = Icges(6,2) * t303 - Icges(6,6) * t611 - t290;
t289 = Icges(6,4) * t303;
t171 = Icges(6,1) * t304 + Icges(6,5) * t611 - t289;
t68 = t430 * t165 - t429 * (t169 * t408 + t171 * t409);
t110 = Icges(5,4) * t223 + Icges(5,2) * t222 + Icges(5,6) * t532;
t112 = Icges(5,1) * t223 + Icges(5,4) * t222 + Icges(5,5) * t532;
t121 = Icges(4,4) * t256 + Icges(4,2) * t255 + Icges(4,6) * t532;
t123 = Icges(4,1) * t256 + Icges(4,4) * t255 + Icges(4,5) * t532;
t707 = -t714 * t430 + (-t110 * t418 + t112 * t419 - t121 * t432 + t123 * t434 + (-t185 * t419 - t188 * t418 - t236 * t434 - t239 * t432) * qJD(3)) * t429;
t109 = Icges(5,4) * t221 + Icges(5,2) * t220 + Icges(5,6) * t530;
t111 = Icges(5,1) * t221 + Icges(5,4) * t220 + Icges(5,5) * t530;
t120 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t530;
t122 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t530;
t706 = t717 * t430 + (t109 * t418 - t111 * t419 + t120 * t432 - t122 * t434 + (t187 * t419 - t189 * t418 + t238 * t434 - t240 * t432) * qJD(3)) * t429;
t705 = -t713 * t430 + (-t185 * t418 + t188 * t419 - t236 * t432 + t239 * t434) * t429;
t316 = (-Icges(5,2) * t419 - t625) * t429;
t294 = qJD(3) * t316;
t317 = (-Icges(5,1) * t418 - t624) * t429;
t295 = qJD(3) * t317;
t350 = (-Icges(4,2) * t434 - t628) * t429;
t331 = qJD(3) * t350;
t351 = (-Icges(4,1) * t432 - t627) * t429;
t332 = qJD(3) * t351;
t702 = t711 * t430 + (-t294 * t418 + t295 * t419 - t331 * t432 + t332 * t434 + (-t280 * t419 - t281 * t418 - t312 * t434 - t313 * t432) * qJD(3)) * t429;
t701 = t719 * t403;
t699 = t720 * t433 + t721 * t435;
t642 = pkin(6) * t429;
t644 = pkin(2) * t430;
t497 = t642 + t644;
t358 = t497 * t435;
t698 = -qJD(1) * t358 + t708;
t697 = -t721 + t724;
t427 = qJD(3) + qJD(5);
t514 = t427 * t429;
t346 = t433 * t514;
t347 = t435 * t514;
t365 = -t427 * t430 + qJD(1);
t301 = -t408 * t610 - t409 * t435;
t302 = t430 * t607 - t603;
t164 = Icges(6,5) * t302 + Icges(6,6) * t301 + Icges(6,3) * t612;
t623 = Icges(6,4) * t302;
t167 = Icges(6,2) * t301 + Icges(6,6) * t612 + t623;
t288 = Icges(6,4) * t301;
t170 = Icges(6,1) * t302 + Icges(6,5) * t612 + t288;
t55 = t164 * t612 + t301 * t167 + t302 * t170;
t696 = -t301 * t169 + t302 * t171;
t56 = -t165 * t612 - t696;
t273 = -Icges(6,3) * t430 + (Icges(6,5) * t409 - Icges(6,6) * t408) * t429;
t621 = Icges(6,4) * t409;
t274 = -Icges(6,6) * t430 + (-Icges(6,2) * t408 + t621) * t429;
t622 = Icges(6,4) * t408;
t275 = -Icges(6,5) * t430 + (Icges(6,1) * t409 - t622) * t429;
t86 = t273 * t612 + t274 * t301 + t275 * t302;
t19 = t346 * t55 - t347 * t56 + t86 * t365;
t58 = t165 * t611 + t169 * t303 + t171 * t304;
t690 = t709 * qJD(3) + t701;
t689 = t710 * qJD(3) - t722;
t671 = t220 * t280 + t221 * t281 + t253 * t312 + t254 * t313 + t294 * t321 - t295 * t322 + t331 * t354 - t332 * t355 + (t711 * t435 + t712 * t554) * t429;
t670 = t222 * t280 + t223 * t281 + t255 * t312 + t256 * t313 + t294 * t319 + t295 * t320 - t331 * t466 + t332 * t353 + (-t711 * t433 + t712 * t553) * t429;
t245 = rSges(4,1) * t355 - rSges(4,2) * t354 + rSges(4,3) * t611;
t684 = t245 * t403;
t318 = -rSges(4,3) * t430 + (rSges(4,1) * t434 - rSges(4,2) * t432) * t429;
t551 = qJD(3) * t429;
t683 = t318 * t551;
t549 = qJD(4) * t429;
t395 = t435 * t549;
t638 = pkin(3) * qJD(3);
t543 = t432 * t638;
t512 = t435 * t543;
t431 = -qJ(4) - pkin(6);
t542 = t434 * t638;
t643 = pkin(3) * t432;
t544 = qJD(1) * t643;
t534 = t431 * t530 + t433 * t542 + t435 * t544;
t679 = t430 * t512 - t395 - t534;
t424 = t434 * pkin(3);
t415 = t424 + pkin(2);
t562 = pkin(3) * t608 + t415 * t609;
t640 = pkin(6) + t431;
t243 = (t429 * t640 + t644) * t435 - t562;
t678 = t403 * t243 + t698;
t676 = t681 * t433;
t675 = t713 * t435;
t174 = rSges(6,1) * t304 - rSges(6,2) * t303 + rSges(6,3) * t611;
t276 = -rSges(6,3) * t430 + (rSges(6,1) * t409 - rSges(6,2) * t408) * t429;
t425 = -pkin(7) + t431;
t557 = t425 - t431;
t371 = pkin(4) * t419 + t424;
t364 = pkin(2) + t371;
t563 = t364 - t415;
t247 = t429 * t563 + t430 * t557;
t641 = -pkin(2) + t415;
t278 = t429 * t641 + t430 * t640;
t580 = -t247 - t278;
t674 = t174 * t365 + (t551 * t580 - qJD(2)) * t435 - t276 * t347;
t449 = t273 * t611 - t274 * t303 + t275 * t304;
t673 = t347 * t58 + t449 * t365;
t531 = t430 * t553;
t555 = qJD(1) * t429;
t672 = t435 * (-rSges(3,2) * t555 - qJD(2)) + rSges(3,1) * t531 + rSges(3,3) * t554;
t341 = t354 * pkin(3);
t340 = t466 * pkin(3);
t669 = t702 * t403;
t666 = -t303 * t167 + t304 * t170;
t664 = (t706 * t435 + t707 * t433 + (t433 * t687 + t435 * t705) * qJD(1)) * t429;
t663 = (t699 * qJD(1) + (-t109 * t319 - t111 * t320 + t120 * t466 - t122 * t353 - t187 * t222 + t189 * t223 - t238 * t255 + t240 * t256 + t532 * t681) * t435 + (t110 * t319 + t112 * t320 - t121 * t466 + t123 * t353 + t185 * t222 + t188 * t223 + t236 * t255 + t239 * t256 + (t714 * t433 + t713 * t553 - t715) * t429) * t433) * t429;
t662 = ((-t109 * t321 + t111 * t322 - t120 * t354 + t122 * t355 - t187 * t220 + t189 * t221 - t238 * t253 + t240 * t254 + (t554 * t681 + t715) * t429 + t704 * qJD(1)) * t435 + (t121 * t354 - t123 * t355 + t236 * t253 + t239 * t254 + t110 * t321 - t112 * t322 + t185 * t220 + t188 * t221 + (-t714 * t435 + t713 * t554) * t429 + t703 * qJD(1)) * t433) * t429;
t399 = rSges(3,1) * t610;
t660 = rSges(3,2) * t612 + rSges(3,3) * t435 - t399;
t659 = (-Icges(4,5) * t354 - Icges(5,5) * t321 - Icges(4,6) * t355 - Icges(5,6) * t322) * t435 + (-Icges(4,5) * t466 + Icges(5,5) * t319 - Icges(4,6) * t353 - Icges(5,6) * t320) * t433;
t450 = t433 * (-Icges(4,2) * t353 + t239 - t337) - t435 * (Icges(4,2) * t355 - t240 + t338);
t658 = t433 * (Icges(4,1) * t466 + t236 + t629) - t435 * (-Icges(4,1) * t354 + t238 - t339);
t452 = t433 * (-Icges(5,2) * t320 + t188 + t306) - t435 * (Icges(5,2) * t322 - t189 + t307);
t657 = t433 * (-Icges(5,1) * t319 + t185 + t626) - t435 * (-Icges(5,1) * t321 + t187 - t308);
t299 = (-Icges(6,2) * t409 - t622) * t429;
t444 = t346 * (-Icges(6,2) * t302 + t170 + t288) - t347 * (Icges(6,2) * t304 - t171 + t289) + t365 * (t275 + t299);
t300 = (-Icges(6,1) * t408 - t621) * t429;
t656 = t346 * (-Icges(6,1) * t301 + t167 + t623) - t347 * (-Icges(6,1) * t303 + t169 - t290) + t365 * (t274 - t300);
t486 = qJD(1) * t514;
t328 = t433 * t486;
t654 = t328 / 0.2e1;
t329 = t435 * t486;
t653 = t329 / 0.2e1;
t652 = -t346 / 0.2e1;
t651 = t346 / 0.2e1;
t650 = -t347 / 0.2e1;
t649 = t347 / 0.2e1;
t645 = -t430 / 0.2e1;
t639 = rSges(3,2) * t429;
t180 = -qJD(1) * t303 - t302 * t427;
t181 = qJD(1) * t304 + t301 * t427;
t100 = Icges(6,1) * t181 + Icges(6,4) * t180 + Icges(6,5) * t532;
t96 = Icges(6,5) * t181 + Icges(6,6) * t180 + Icges(6,3) * t532;
t98 = Icges(6,4) * t181 + Icges(6,2) * t180 + Icges(6,6) * t532;
t30 = -t430 * t96 + ((-t167 * t427 + t100) * t409 + (-t170 * t427 - t98) * t408) * t429;
t637 = t30 * t346;
t178 = qJD(1) * t301 + t304 * t427;
t179 = qJD(1) * t302 + t303 * t427;
t95 = Icges(6,5) * t179 + Icges(6,6) * t178 + Icges(6,3) * t530;
t97 = Icges(6,4) * t179 + Icges(6,2) * t178 + Icges(6,6) * t530;
t99 = Icges(6,1) * t179 + Icges(6,4) * t178 + Icges(6,5) * t530;
t31 = -t430 * t95 + ((-t169 * t427 + t99) * t409 + (t171 * t427 - t97) * t408) * t429;
t636 = t31 * t347;
t67 = -t164 * t430 + (-t167 * t408 + t170 * t409) * t429;
t633 = t67 * t329;
t632 = t68 * t328;
t631 = rSges(3,3) + qJ(2);
t490 = rSges(6,1) * t179 + rSges(6,2) * t178;
t101 = rSges(6,3) * t530 + t490;
t630 = t430 * t101 + t276 * t530;
t613 = t425 * t429;
t370 = pkin(4) * t418 + t643;
t359 = t433 * t370;
t138 = (t430 * t641 - t642) * t554 + t679;
t600 = t430 * t138 + t278 * t530;
t333 = t364 * t610;
t367 = t415 * t610;
t519 = t557 * t429;
t521 = -t370 + t643;
t146 = -t433 * t519 + t435 * t521 + t333 - t367;
t357 = pkin(2) * t610 + pkin(6) * t612;
t496 = -t431 * t612 + t367;
t242 = -pkin(3) * t601 - t357 + t496;
t598 = t146 + t242;
t173 = t302 * rSges(6,1) + t301 * rSges(6,2) + rSges(6,3) * t612;
t597 = t173 * t611 - t174 * t612;
t194 = t320 * rSges(5,1) + t319 * rSges(5,2) + rSges(5,3) * t612;
t587 = t194 + t242;
t585 = t242 * t611 + t243 * t612;
t548 = qJD(4) * t430;
t360 = -t429 * t543 - t548;
t361 = t370 * qJD(3);
t498 = -t361 + t543;
t577 = -t498 * t429 - t360;
t287 = -rSges(5,3) * t430 + (rSges(5,1) * t419 - rSges(5,2) * t418) * t429;
t576 = -t278 - t287;
t575 = t280 - t317;
t574 = t281 + t316;
t426 = t429 ^ 2;
t513 = t433 * t543;
t573 = -t403 * t340 + t426 * t513;
t528 = t435 * t551;
t529 = t433 * t551;
t572 = -t340 * t528 + t341 * t529;
t323 = (-rSges(5,1) * t418 - rSges(5,2) * t419) * t429;
t571 = -qJD(3) * t323 - t360;
t570 = t312 - t351;
t569 = t313 + t350;
t552 = qJD(2) * t435;
t327 = qJD(1) * (-t552 + t559);
t560 = pkin(2) * t531 + pkin(6) * t532;
t568 = qJD(1) * t560 + t327;
t567 = -t364 * t609 - t359;
t421 = qJD(2) * t433;
t558 = qJ(2) * t553 + t421;
t345 = pkin(1) * t554 - t558;
t566 = -t497 * t554 - t345;
t564 = t415 * t531 + t433 * t544;
t493 = rSges(3,1) * t430 - t639;
t326 = rSges(3,3) * t433 + t435 * t493;
t250 = -t552 + (t326 + t377) * qJD(1);
t556 = qJD(1) * t250;
t546 = qJD(1) * qJD(2);
t102 = t181 * rSges(6,1) + t180 * rSges(6,2) + rSges(6,3) * t532;
t545 = t101 * t612 + t102 * t611 - t174 * t532;
t394 = t433 * t549;
t139 = -pkin(3) * t456 - t431 * t532 + t394 - t560 + t564;
t541 = t138 * t529 + (qJD(1) * t243 + t139) * t528;
t540 = t138 * t612 + t139 * t611 + t243 * t532;
t539 = -t173 - t598;
t114 = t223 * rSges(5,1) + t222 * rSges(5,2) + rSges(5,3) * t532;
t538 = -t276 + t580;
t125 = t256 * rSges(4,1) + t255 * rSges(4,2) + rSges(4,3) * t532;
t267 = t278 * t529;
t410 = t433 * t546;
t525 = qJD(1) * t549;
t537 = qJD(1) * t267 + t435 * t525 + t410;
t305 = (-rSges(6,1) * t408 - rSges(6,2) * t409) * t429;
t272 = t427 * t305;
t536 = -t272 + t577;
t535 = -t361 * t610 + t364 * t531 + t370 * t554;
t244 = t353 * rSges(4,1) - rSges(4,2) * t466 + rSges(4,3) * t612;
t527 = t612 / 0.2e1;
t526 = -t611 / 0.2e1;
t524 = t555 / 0.2e1;
t523 = -t551 / 0.2e1;
t522 = t551 / 0.2e1;
t423 = t433 * pkin(1);
t375 = -qJ(2) * t435 + t423;
t477 = (t357 + t375) * qJD(1) - t421;
t459 = -t267 - t395 + t477;
t508 = t287 * t529;
t60 = t403 * t587 + t459 - t508;
t520 = t60 * t576;
t213 = rSges(6,1) * t301 - rSges(6,2) * t302;
t214 = rSges(6,1) * t303 + rSges(6,2) * t304;
t518 = t213 * t347 + t346 * t214;
t517 = t365 * t213 - t305 * t346;
t516 = -t214 * t365 - t347 * t305;
t515 = (t358 + t377) * qJD(1);
t509 = t247 * t529;
t507 = t318 * t529;
t505 = t435 * t524;
t504 = t433 * t524;
t503 = t433 * t523;
t502 = t433 * t522;
t501 = t435 * t523;
t500 = t435 * t522;
t499 = qJD(1) * t522;
t494 = t242 * t528 + t243 * t529 - t548;
t492 = rSges(4,1) * t254 + rSges(4,2) * t253;
t491 = rSges(5,1) * t221 + rSges(5,2) * t220;
t124 = rSges(4,3) * t530 + t492;
t465 = t429 * (-rSges(4,1) * t432 - rSges(4,2) * t434);
t335 = qJD(3) * t465;
t70 = -t335 * t528 - t124 * t403 + t410 + (t507 + t566) * qJD(1);
t464 = (-qJD(2) - t683) * t435;
t71 = qJD(1) * t464 + t125 * t403 - t335 * t529 + t568;
t488 = -t71 * t433 - t70 * t435;
t104 = t244 * t403 + t477 - t507;
t105 = t464 + t515 + t684;
t485 = -t104 * t433 - t105 * t435;
t480 = t429 * (t403 + t550);
t478 = t394 + t515;
t476 = t433 * t499;
t475 = t435 * t499;
t474 = pkin(1) + t493;
t57 = -t164 * t611 - t666;
t195 = rSges(5,1) * t322 - rSges(5,2) * t321 + rSges(5,3) * t611;
t298 = (-Icges(6,5) * t408 - Icges(6,6) * t409) * t429;
t458 = t430 * t513 - t564;
t457 = (Icges(6,5) * t301 - Icges(6,6) * t302) * t346 - (Icges(6,5) * t303 + Icges(6,6) * t304) * t347 + t298 * t365;
t117 = (t244 * t435 - t245 * t433) * t551;
t455 = t403 * t139 + t433 * t525 - t435 * t546 + t568;
t11 = -t100 * t304 + t167 * t178 + t170 * t179 + t303 * t98 + (t164 * t554 - t435 * t96) * t429;
t12 = t169 * t178 - t171 * t179 + t303 * t97 - t304 * t99 + (-t165 * t554 - t435 * t95) * t429;
t13 = t100 * t302 + t167 * t180 + t170 * t181 + t301 * t98 + (t164 * t553 + t433 * t96) * t429;
t14 = t169 * t180 - t171 * t181 + t301 * t97 + t302 * t99 + (-t165 * t553 + t433 * t95) * t429;
t20 = t346 * t57 - t673;
t269 = t427 * t298;
t270 = t427 * t299;
t271 = t427 * t300;
t48 = t178 * t274 + t179 * t275 + t270 * t303 - t271 * t304 + (-t269 * t435 + t273 * t554) * t429;
t49 = t180 * t274 + t181 * t275 + t270 * t301 + t271 * t302 + (t269 * t433 + t273 * t553) * t429;
t69 = -t269 * t430 + ((-t274 * t427 + t271) * t409 + (-t275 * t427 - t270) * t408) * t429;
t62 = t69 * t365;
t446 = (t13 * t346 - t14 * t347 + t328 * t56 + t329 * t55 + t365 * t49) * t527 - (-t457 * t430 + (-t408 * t444 - t409 * t656) * t429) * t365 / 0.2e1 + t20 * t504 + t19 * t505 + (t11 * t346 - t12 * t347 + t328 * t58 + t329 * t57 + t365 * t48) * t526 + (t430 * t449 + (t433 * t57 - t435 * t58) * t429) * t654 + (-t430 * t86 + (t433 * t55 - t435 * t56) * t429) * t653 + (-t430 * t49 + (t13 * t433 - t14 * t435 + (t433 * t56 + t435 * t55) * qJD(1)) * t429) * t651 + (-t430 * t48 + (t11 * t433 - t12 * t435 + (t433 * t58 + t435 * t57) * qJD(1)) * t429) * t650 + (t62 + t632 + t633 - t636 + t637) * t645 + t365 * (-t430 * t69 + (t30 * t433 - t31 * t435 + (t433 * t68 + t435 * t67) * qJD(1)) * t429) / 0.2e1 + (t301 * t444 - t302 * t656 + t457 * t612) * t652 + (t444 * t303 + t304 * t656 - t457 * t611) * t649;
t373 = t426 * t512;
t362 = t371 * qJD(3);
t296 = t521 * t429;
t286 = t435 * t480;
t285 = t433 * t480;
t264 = rSges(4,1) * t354 + rSges(4,2) * t355;
t263 = -rSges(4,1) * t466 - rSges(4,2) * t353;
t231 = rSges(5,1) * t321 + rSges(5,2) * t322;
t230 = rSges(5,1) * t319 - rSges(5,2) * t320;
t219 = t430 * t243;
t204 = t370 * t609 - t371 * t433 - t341;
t203 = -t359 * t430 - t371 * t435 + t340;
t193 = qJD(1) * t672 + t327;
t192 = t410 + (qJD(1) * t660 - t345) * qJD(1);
t155 = t430 * t174;
t147 = t557 * t611 + t562 + t567;
t113 = rSges(5,3) * t530 + t491;
t84 = (-qJD(1) * t519 - t362 + t542) * t435 + t458 + t535;
t83 = -t362 * t433 - t498 * t609 + (-t370 * t435 + (t430 * t563 - t613) * t433) * qJD(1) + t534;
t61 = (t195 - t243) * t403 + (t551 * t576 - qJD(2)) * t435 + t478;
t54 = (t194 * t435 - t195 * t433) * t551 + t494;
t47 = (-t147 - t243) * t403 + t478 + t674;
t46 = t173 * t365 - t276 * t346 + t403 * t598 + t459 - t509;
t44 = t114 * t403 + (t433 * t571 + t553 * t576) * t551 + t455;
t43 = (-t113 - t138) * t403 + t571 * t528 + (t508 + t566) * qJD(1) + t537;
t40 = t173 * t347 - t174 * t346 + (t146 * t435 + t147 * t433) * t551 + t494;
t23 = (t113 * t433 + t114 * t435 + (-t195 * t435 - t433 * t587) * qJD(1)) * t551 + t541;
t16 = t102 * t365 - t272 * t346 - t276 * t329 + t403 * t84 + (t433 * t577 + t553 * t580) * t551 + t455;
t15 = -t101 * t365 - t272 * t347 + t276 * t328 + (-t138 - t83) * t403 + t577 * t528 + (t509 + t566) * qJD(1) + t537;
t7 = t101 * t346 + t102 * t347 - t173 * t328 - t174 * t329 + (t433 * t83 + t435 * t84 + (t147 * t435 - t433 * t598) * qJD(1)) * t551 + t541;
t1 = [t669 + t19 * t649 + t49 * t651 + t86 * t653 - t449 * t654 + t62 + t637 / 0.2e1 - t636 / 0.2e1 + t632 / 0.2e1 + t633 / 0.2e1 + (t20 + (t56 + (t164 * t435 + t165 * t433) * t429 + t666 + t696) * t346 + t673) * t652 + (t19 + t48) * t650 + ((((t675 - t676) * t429 + t677 + t704 - t720) * t435 + (-t697 + t703) * t433) * t551 + t701) * t500 + (t15 * (t174 - t567) + t47 * (t395 - t490 + t558) + t16 * (t333 + t423 + t173) + (t15 * qJ(2) - t16 * t613 + (t362 + (-t364 * t430 - pkin(1) + (-rSges(6,3) + t425) * t429) * qJD(1)) * t47) * t433 + (t15 * (pkin(1) - t613) + t47 * (qJD(1) * t370 - t361 * t430) + t16 * (-qJ(2) - t370)) * t435 + (t535 + t102 + (-t425 * t555 - qJD(2) - t362) * t435 + t147 * t403 + t47 - t674 + t678) * t46) * m(6) + (t43 * (t195 + t562 + t620) + t44 * (t423 + t496 + t194) + (-t491 + t558 + (-rSges(5,3) * t429 - t415 * t430 - pkin(1)) * t554 - t679) * t61 + (-t195 * t403 + t114 - t458 + t61 + t678) * t60 + (-t520 * t551 + t43 * (-t429 * t431 + pkin(1)) + t44 * (-qJ(2) - t643) + t60 * (-t431 * t555 - t542)) * t435) * m(5) + (t70 * (t245 + t620) + t71 * (t423 + t244 + t357) + (-t492 + t558 + (-t644 - pkin(1) + (-rSges(4,3) - pkin(6)) * t429) * t554) * t105 + (t105 + t125 + t560 - t684 + t698) * t104 + (t104 * t683 + t70 * (pkin(1) + t497) - t71 * qJ(2)) * t435) * m(4) + (t250 * t558 + t193 * (t399 + t423) + (t192 * t631 - t193 * t639 - t474 * t556) * t433 + (rSges(3,3) * t556 + t192 * t474 - t193 * t631) * t435 + (-qJD(1) * t326 + t250 + t552 + t672 + t708) * (-t421 + (t375 - t660) * qJD(1))) * m(3) + (t670 + t707) * t502 + (-t680 + t687) * t476 + (t719 + t705) * t475 + (t671 + t690 - t706) * t501 + ((t697 * t435 + (t677 + t723) * t433 + (t681 * t435 ^ 2 + (t675 + t676) * t433) * t429 + t699) * t551 + t689 + t722) * t503; m(3) * (-t192 * t435 - t193 * t433) + m(4) * t488 + m(5) * (-t43 * t435 - t433 * t44) + m(6) * (-t15 * t435 - t16 * t433); t446 - (((-t418 * t574 - t419 * t575 - t432 * t569 - t434 * t570) * t403 + ((-t418 * t452 - t419 * t657 - t432 * t450 - t434 * t658) * t429 - t659 * t430) * qJD(3)) * t429 + t665 * t403 * t430) * t403 / 0.2e1 + (-t702 * t430 + t664) * t403 / 0.2e1 + (qJD(3) * t664 + t669) * t645 + (qJD(3) * t663 + t403 * t670) * t527 + (qJD(3) * t662 + t403 * t671) * t526 + t690 * t505 + t689 * t504 + ((t319 * t452 - t320 * t657 - t353 * t658 - t450 * t466 + t612 * t659) * t551 + (t319 * t574 - t320 * t575 - t353 * t570 - t466 * t569 - t612 * t665) * t403) * t503 + (-t430 * t670 + t663) * t502 + (-t430 * t671 + t662) * t501 + ((t452 * t321 + t322 * t657 + t450 * t354 + t355 * t658 - t611 * t659) * t551 + (t321 * t574 + t322 * t575 + t354 * t569 + t355 * t570 + t611 * t665) * t403) * t500 + (t430 * t680 + t710) * t476 + (-t430 * t719 + t709) * t475 + (-t40 * ((t203 * t435 + t204 * t433) * t551 + t518 + t572) - t47 * (-t296 * t528 + t373 + (-t204 - t341) * t403 + t516) - t46 * (t203 * t403 - t296 * t529 + t517 + t573) + t7 * (t585 + t597) + t40 * (t540 + t545) + t15 * (-t155 + t219) + t47 * (t600 + t630) + (t15 * t147 + t47 * t83 + t16 * t539 + t46 * (-t102 - t139 - t84)) * t430 + ((t7 * t146 + t40 * t84 + t15 * t538 + t47 * t536 + (t40 * t147 + t46 * t538) * qJD(1)) * t435 + (t7 * t147 + t40 * t83 + t16 * t538 + t46 * t536 + (t47 * t247 + t40 * t539) * qJD(1)) * t433) * t429) * m(6) + (-t54 * ((t230 * t435 + t231 * t433) * t551 + t572) - t61 * (-t323 * t528 + t373 + (-t231 - t341) * t403) - t60 * (t230 * t403 - t323 * t529 + t573) + t23 * t585 + t54 * t540 + t43 * t219 + t61 * t600 + (-t43 * t195 + t61 * t113 - t44 * t587 + t60 * (-t114 - t139)) * t430 + ((t23 * t194 + t54 * t114 + t43 * t576 + t61 * t571 + (-t195 * t54 + t520) * qJD(1)) * t435 + (-t23 * t195 + t54 * t113 + t44 * t576 + t60 * t571 + (t61 * t287 - t54 * t587) * qJD(1)) * t433) * t429) * m(5) + ((-t104 * t125 + t105 * t124 - t244 * t71 - t70 * t245) * t430 + (0.2e1 * t117 * (t124 * t433 + t125 * t435 + (-t244 * t433 - t245 * t435) * qJD(1)) + t485 * t335 + ((-t104 * t435 + t105 * t433) * qJD(1) + t488) * t318) * t429 - (t104 * t263 - t105 * t264) * t403 - (t117 * (t263 * t435 + t264 * t433) + t485 * t465) * t551) * m(4); (-t285 * t46 - t286 * t47 - t430 * t7) * m(6) + (-t23 * t430 + (-t286 + t532) * t61 + (-t285 + t530) * t60) * m(5) + (m(6) * (t15 * t433 - t16 * t435 + t46 * t554 + t47 * t553) + (t43 * t433 - t435 * t44) * m(5)) * t429; t446 + (t7 * t597 + t15 * (-t276 * t611 - t155) + t16 * (-t173 * t430 - t276 * t612) + (-t272 * t611 - t516 + t630) * t47 + (-t102 * t430 + (-t272 * t433 - t276 * t553) * t429 - t517) * t46 + (-t173 * t530 - t518 + t545) * t40) * m(6);];
tauc = t1(:);
