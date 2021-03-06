% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR13
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR13_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:33
% EndTime: 2019-12-31 19:15:20
% DurationCPUTime: 36.91s
% Computational Cost: add. (23518->1189), mult. (41912->1605), div. (0->0), fcn. (39803->8), ass. (0->540)
t440 = sin(qJ(3));
t442 = cos(qJ(4));
t444 = cos(qJ(1));
t654 = t442 * t444;
t439 = sin(qJ(4));
t441 = sin(qJ(1));
t657 = t441 * t439;
t327 = -t440 * t657 + t654;
t656 = t441 * t442;
t662 = t439 * t444;
t328 = t440 * t656 + t662;
t443 = cos(qJ(3));
t655 = t441 * t443;
t190 = Icges(5,5) * t328 + Icges(5,6) * t327 - Icges(5,3) * t655;
t660 = t440 * t444;
t329 = t439 * t660 + t656;
t330 = t440 * t654 - t657;
t653 = t443 * t444;
t192 = -Icges(5,5) * t330 + Icges(5,6) * t329 + Icges(5,3) * t653;
t316 = Icges(5,4) * t330;
t195 = Icges(5,2) * t329 + Icges(5,6) * t653 - t316;
t315 = Icges(5,4) * t329;
t197 = Icges(5,1) * t330 - Icges(5,5) * t653 - t315;
t508 = -t195 * t329 - t197 * t330;
t686 = Icges(5,4) * t328;
t193 = Icges(5,2) * t327 - Icges(5,6) * t655 + t686;
t314 = Icges(5,4) * t327;
t196 = Icges(5,1) * t328 - Icges(5,5) * t655 + t314;
t650 = t193 * t327 + t196 * t328;
t791 = t508 + t650 + (-t190 * t441 - t192 * t444) * t443;
t438 = qJ(4) + qJ(5);
t426 = cos(t438);
t425 = sin(t438);
t658 = t441 * t425;
t302 = t426 * t444 - t440 * t658;
t661 = t440 * t441;
t303 = t425 * t444 + t426 * t661;
t157 = Icges(6,5) * t303 + Icges(6,6) * t302 - Icges(6,3) * t655;
t304 = t425 * t660 + t426 * t441;
t305 = t426 * t660 - t658;
t159 = -Icges(6,5) * t305 + Icges(6,6) * t304 + Icges(6,3) * t653;
t281 = Icges(6,4) * t305;
t162 = Icges(6,2) * t304 + Icges(6,6) * t653 - t281;
t280 = Icges(6,4) * t304;
t164 = Icges(6,1) * t305 - Icges(6,5) * t653 - t280;
t511 = -t162 * t304 - t164 * t305;
t683 = Icges(6,4) * t303;
t160 = Icges(6,2) * t302 - Icges(6,6) * t655 + t683;
t279 = Icges(6,4) * t302;
t163 = Icges(6,1) * t303 - Icges(6,5) * t655 + t279;
t651 = t160 * t302 + t163 * t303;
t790 = t511 + t651 + (-t157 * t441 - t159 * t444) * t443;
t515 = Icges(5,5) * t442 - Icges(5,6) * t439;
t287 = Icges(5,3) * t440 + t443 * t515;
t684 = Icges(5,4) * t442;
t518 = -Icges(5,2) * t439 + t684;
t291 = Icges(5,6) * t440 + t443 * t518;
t685 = Icges(5,4) * t439;
t521 = Icges(5,1) * t442 - t685;
t295 = Icges(5,5) * t440 + t443 * t521;
t113 = t287 * t653 + t291 * t329 - t295 * t330;
t422 = qJD(3) * t444;
t607 = qJD(4) * t443;
t361 = -t441 * t607 + t422;
t73 = t190 * t653 + t193 * t329 - t196 * t330;
t421 = qJD(4) * t440;
t754 = qJD(1) + t421;
t789 = -t113 * t754 - t361 * t73;
t437 = qJD(4) + qJD(5);
t542 = t437 * t655;
t285 = t422 - t542;
t371 = qJD(5) * t440 + t754;
t67 = t157 * t653 + t160 * t304 - t163 * t305;
t514 = Icges(6,5) * t426 - Icges(6,6) * t425;
t269 = Icges(6,3) * t440 + t443 * t514;
t681 = Icges(6,4) * t426;
t517 = -Icges(6,2) * t425 + t681;
t271 = Icges(6,6) * t440 + t443 * t517;
t682 = Icges(6,4) * t425;
t520 = Icges(6,1) * t426 - t682;
t273 = Icges(6,5) * t440 + t443 * t520;
t97 = t269 * t653 + t271 * t304 - t273 * t305;
t788 = -t285 * t67 - t371 * t97;
t536 = rSges(6,1) * t305 - rSges(6,2) * t304;
t168 = rSges(6,3) * t653 - t536;
t348 = pkin(3) * t660 - pkin(7) * t653;
t418 = pkin(4) * t442 + pkin(3);
t445 = -pkin(8) - pkin(7);
t652 = t443 * t445;
t623 = t418 * t660 + t444 * t652;
t205 = pkin(4) * t657 + t348 - t623;
t535 = rSges(6,1) * t426 - rSges(6,2) * t425;
t275 = rSges(6,3) * t440 + t443 * t535;
t373 = t437 * t443;
t610 = qJD(3) * t441;
t284 = t373 * t444 + t610;
t606 = qJD(4) * t444;
t360 = t443 * t606 + t610;
t714 = pkin(7) + t445;
t715 = pkin(3) - t418;
t744 = t440 * t714 + t443 * t715;
t785 = -t168 * t371 - t205 * t754 + t275 * t284 - t360 * t744;
t611 = qJD(1) * t444;
t717 = pkin(4) * t439;
t741 = -pkin(1) - pkin(6);
t559 = -t717 + t741;
t784 = t559 * qJD(1);
t507 = t195 * t439 + t197 * t442;
t78 = t192 * t440 - t443 * t507;
t510 = t162 * t425 + t164 * t426;
t70 = t159 * t440 - t443 * t510;
t66 = -t159 * t655 + t162 * t302 - t164 * t303;
t72 = -t192 * t655 + t195 * t327 - t197 * t328;
t203 = -rSges(5,1) * t330 + rSges(5,2) * t329 + rSges(5,3) * t653;
t779 = t203 * t754;
t434 = t444 * rSges(4,3);
t299 = rSges(4,1) * t661 + rSges(4,2) * t655 + t434;
t386 = pkin(1) * t444 + qJ(2) * t441;
t435 = t444 * pkin(6);
t756 = t435 + t386;
t777 = t299 + t756;
t571 = -t422 / 0.2e1;
t612 = qJD(1) * t443;
t576 = -t612 / 0.2e1;
t776 = t440 * t571 + t441 * t576;
t573 = t610 / 0.2e1;
t775 = t440 * t573 + t444 * t576;
t609 = qJD(3) * t443;
t583 = t441 * t609;
t774 = t440 * t611 + t583;
t584 = t440 * t422;
t773 = t441 * t612 + t584;
t406 = Icges(4,4) * t655;
t680 = Icges(4,5) * t444;
t296 = Icges(4,1) * t661 + t406 + t680;
t687 = Icges(4,4) * t443;
t522 = Icges(4,1) * t440 + t687;
t297 = -Icges(4,5) * t441 + t444 * t522;
t377 = -Icges(4,2) * t440 + t687;
t339 = t377 * t444;
t466 = t441 * (t297 + t339) - t444 * (-Icges(4,2) * t661 + t296 + t406);
t688 = Icges(4,4) * t440;
t519 = Icges(4,2) * t443 + t688;
t292 = Icges(4,6) * t444 + t441 * t519;
t293 = -Icges(4,6) * t441 + t444 * t519;
t379 = Icges(4,1) * t443 - t688;
t341 = t379 * t441;
t342 = t379 * t444;
t467 = t441 * (t293 - t342) - t444 * (t292 - t341);
t772 = -t467 * t440 + t466 * t443;
t621 = t377 + t522;
t622 = -t519 + t379;
t771 = (t440 * t621 - t443 * t622) * qJD(1);
t767 = 0.2e1 * qJD(3);
t389 = pkin(3) * t443 + pkin(7) * t440;
t352 = t389 * t610;
t423 = qJD(2) * t441;
t428 = t444 * qJ(2);
t382 = pkin(1) * t441 - t428;
t716 = pkin(6) * t441;
t569 = -t382 - t716;
t462 = t352 + t423 + (t348 + t569) * qJD(1);
t61 = t462 + t785;
t636 = -t744 + t275;
t766 = t61 * t636;
t633 = rSges(6,1) * t303 + rSges(6,2) * t302;
t166 = -rSges(6,3) * t655 + t633;
t410 = pkin(3) * t661;
t346 = -pkin(7) * t655 + t410;
t411 = pkin(4) * t662;
t595 = t418 * t661 + t441 * t652 + t411;
t204 = -t346 + t595;
t424 = qJD(2) * t444;
t461 = (t346 + t756) * qJD(1) - t389 * t422 - t424;
t62 = t166 * t371 + t204 * t754 - t275 * t285 + t361 * t744 + t461;
t765 = t62 * t636;
t96 = -t269 * t655 + t271 * t302 + t273 * t303;
t764 = t284 * t66 + t371 * t96;
t763 = t440 * t715;
t286 = Icges(5,3) * t443 - t440 * t515;
t504 = t291 * t439 - t295 * t442;
t509 = t193 * t439 - t196 * t442;
t447 = t360 * (-t287 * t444 + t507) + t361 * (t287 * t441 + t509) + t754 * (t286 + t504);
t762 = t447 * t443;
t501 = t293 * t443 + t297 * t440;
t759 = t501 * t444;
t112 = -t287 * t655 + t291 * t327 + t295 * t328;
t758 = t112 * t754 + t360 * t72;
t566 = -rSges(3,2) * t444 + rSges(3,3) * t441;
t757 = t386 + t566;
t614 = qJD(1) * t440;
t554 = qJD(4) + t614;
t582 = t443 * t422;
t753 = t441 * t554 - t582;
t560 = t437 + t614;
t752 = t441 * t560 - t582;
t538 = rSges(5,1) * t442 - rSges(5,2) * t439;
t300 = rSges(5,3) * t440 + t443 * t538;
t751 = t444 * t560 + t583;
t175 = t293 * t440 - t297 * t443;
t210 = qJD(1) * t292 - qJD(3) * t339;
t213 = -qJD(3) * t342 + (t441 * t522 + t680) * qJD(1);
t516 = Icges(4,5) * t440 + Icges(4,6) * t443;
t289 = -Icges(4,3) * t441 + t444 * t516;
t616 = qJD(1) * t289;
t750 = qJD(3) * t175 + t210 * t443 + t213 * t440 + t616;
t363 = t519 * qJD(3);
t364 = t522 * qJD(3);
t375 = Icges(4,5) * t443 - Icges(4,6) * t440;
t749 = qJD(1) * t375 + qJD(3) * (t377 * t440 - t379 * t443) + t363 * t443 + t364 * t440;
t211 = qJD(1) * t293 + t377 * t610;
t214 = qJD(1) * t297 + qJD(3) * t341;
t502 = t292 * t440 - t296 * t443;
t288 = Icges(4,3) * t444 + t441 * t516;
t617 = qJD(1) * t288;
t748 = qJD(3) * t502 - t211 * t443 - t214 * t440 + t617;
t706 = pkin(4) * qJD(4);
t600 = t439 * t706;
t557 = t440 * t600;
t592 = pkin(3) * t582 + pkin(7) * t773;
t599 = t442 * t706;
t663 = t418 * t443;
t106 = t441 * t599 + (t557 + (t440 * t445 - t663) * qJD(3)) * t444 + (t411 + (t652 - t763) * t441) * qJD(1) + t592;
t587 = t443 * t611;
t585 = t440 * t610;
t593 = pkin(3) * t774 + pkin(7) * t585;
t476 = pkin(7) * t587 - t593;
t494 = t754 * t439;
t553 = t418 * t774 + t444 * t599 + t445 * t587;
t608 = qJD(3) * t445;
t107 = (-pkin(4) * t494 - t440 * t608) * t441 + t476 + t553;
t581 = t441 * t421;
t390 = qJD(3) * t581;
t558 = t444 * t437;
t199 = qJD(5) * t585 + t390 + (-t443 * t558 - t610) * qJD(1);
t326 = t440 * t558;
t603 = qJD(1) * qJD(3);
t413 = t444 * t603;
t200 = -qJD(1) * t542 - qJD(3) * t326 + t413;
t264 = -qJD(1) * t360 + t390;
t265 = -qJD(4) * t773 + t413;
t235 = qJD(1) * t410 - t592;
t453 = t235 * t422 + (t476 * t441 + (-t346 * t444 + t348 * t441) * qJD(1)) * qJD(3);
t561 = t437 * t440 + qJD(1);
t497 = t444 * t561;
t132 = -t425 * t752 + t426 * t497;
t133 = t425 * t497 + t426 * t752;
t537 = rSges(6,1) * t133 + rSges(6,2) * t132;
t85 = -rSges(6,3) * t773 + t537;
t498 = t441 * t561;
t134 = -t425 * t751 - t426 * t498;
t135 = -t425 * t498 + t426 * t751;
t598 = rSges(6,1) * t135 + rSges(6,2) * t134 + rSges(6,3) * t585;
t86 = -rSges(6,3) * t587 + t598;
t11 = t106 * t361 - t107 * t360 - t166 * t200 + t168 * t199 - t204 * t265 + t205 * t264 - t284 * t86 + t285 * t85 + t453;
t545 = -t346 * t610 - t348 * t422;
t52 = -t166 * t284 + t168 * t285 - t204 * t360 + t205 * t361 + t545;
t643 = t168 + t205;
t692 = t106 + t85;
t747 = t11 * t643 + t52 * t692;
t268 = Icges(6,3) * t443 - t440 * t514;
t505 = t271 * t425 - t273 * t426;
t512 = t160 * t425 - t163 * t426;
t448 = t284 * (-t269 * t444 + t510) + t285 * (t269 * t441 + t512) + t371 * (t268 + t505);
t746 = -(t157 * t285 + t159 * t284 + t269 * t371) * t440 + t443 * t448;
t310 = (-Icges(6,2) * t426 - t682) * t443;
t457 = t284 * (Icges(6,2) * t305 - t164 + t280) + t285 * (-Icges(6,2) * t303 + t163 + t279) + t371 * (t273 + t310);
t311 = (-Icges(6,1) * t425 - t681) * t443;
t743 = t284 * (-Icges(6,1) * t304 + t162 - t281) + t285 * (-Icges(6,1) * t302 + t160 + t683) + t371 * (t271 - t311);
t337 = (-Icges(5,2) * t442 - t685) * t443;
t455 = t360 * (Icges(5,2) * t330 - t197 + t315) + t361 * (-Icges(5,2) * t328 + t196 + t314) + t754 * (t295 + t337);
t340 = (-Icges(5,1) * t439 - t684) * t443;
t742 = t360 * (-Icges(5,1) * t329 + t195 - t316) + t361 * (-Icges(5,1) * t327 + t193 + t686) + t754 * (t291 - t340);
t740 = t199 / 0.2e1;
t739 = t200 / 0.2e1;
t738 = t264 / 0.2e1;
t737 = t265 / 0.2e1;
t736 = -t284 / 0.2e1;
t735 = t284 / 0.2e1;
t734 = -t285 / 0.2e1;
t733 = t285 / 0.2e1;
t359 = qJD(3) * t373;
t732 = t359 / 0.2e1;
t731 = -t360 / 0.2e1;
t730 = t360 / 0.2e1;
t729 = -t361 / 0.2e1;
t728 = t361 / 0.2e1;
t727 = -t371 / 0.2e1;
t726 = t371 / 0.2e1;
t725 = -t754 / 0.2e1;
t724 = t754 / 0.2e1;
t723 = t440 / 0.2e1;
t722 = t441 / 0.2e1;
t721 = -t444 / 0.2e1;
t719 = rSges(3,2) - pkin(1);
t718 = pkin(3) * t440;
t712 = rSges(4,2) * t440;
t711 = rSges(3,3) * t444;
t709 = rSges(5,3) * t443;
t707 = rSges(6,3) * t443;
t478 = t585 - t587;
t80 = Icges(6,5) * t135 + Icges(6,6) * t134 + Icges(6,3) * t478;
t82 = Icges(6,4) * t135 + Icges(6,2) * t134 + Icges(6,6) * t478;
t84 = Icges(6,1) * t135 + Icges(6,4) * t134 + Icges(6,5) * t478;
t18 = (qJD(3) * t512 + t80) * t440 + (qJD(3) * t157 + (-t160 * t437 + t84) * t426 + (-t163 * t437 - t82) * t425) * t443;
t705 = t18 * t285;
t79 = Icges(6,5) * t133 + Icges(6,6) * t132 - Icges(6,3) * t773;
t81 = Icges(6,4) * t133 + Icges(6,2) * t132 - Icges(6,6) * t773;
t83 = Icges(6,1) * t133 + Icges(6,4) * t132 - Icges(6,5) * t773;
t19 = (qJD(3) * t510 + t79) * t440 + (qJD(3) * t159 + (-t162 * t437 + t83) * t426 + (t164 * t437 - t81) * t425) * t443;
t704 = t19 * t284;
t178 = -t754 * t656 + (-t444 * t554 - t583) * t439;
t179 = t554 * t654 + (t442 * t609 - t494) * t441;
t101 = Icges(5,4) * t179 + Icges(5,2) * t178 + Icges(5,6) * t478;
t103 = Icges(5,1) * t179 + Icges(5,4) * t178 + Icges(5,5) * t478;
t99 = Icges(5,5) * t179 + Icges(5,6) * t178 + Icges(5,3) * t478;
t27 = (qJD(3) * t509 + t99) * t440 + (qJD(3) * t190 - t101 * t439 + t103 * t442 + (-t193 * t442 - t196 * t439) * qJD(4)) * t443;
t703 = t27 * t361;
t495 = t444 * t754;
t176 = -t439 * t753 + t442 * t495;
t177 = t439 * t495 + t442 * t753;
t100 = Icges(5,4) * t177 + Icges(5,2) * t176 - Icges(5,6) * t773;
t102 = Icges(5,1) * t177 + Icges(5,4) * t176 - Icges(5,5) * t773;
t98 = Icges(5,5) * t177 + Icges(5,6) * t176 - Icges(5,3) * t773;
t28 = (qJD(3) * t507 + t98) * t440 + (qJD(3) * t192 - t100 * t439 + t102 * t442 + (-t195 * t442 + t197 * t439) * qJD(4)) * t443;
t702 = t28 * t360;
t698 = t62 * t204;
t69 = t157 * t440 - t443 * t512;
t697 = t69 * t199;
t696 = t70 * t200;
t77 = t190 * t440 - t443 * t509;
t695 = t77 * t264;
t694 = t78 * t265;
t625 = rSges(5,1) * t328 + rSges(5,2) * t327;
t201 = -rSges(5,3) * t655 + t625;
t94 = t201 * t754 - t300 * t361 + t461;
t693 = t94 * t300;
t691 = -t107 - t86;
t114 = t269 * t440 - t443 * t505;
t309 = (-Icges(6,5) * t425 - Icges(6,6) * t426) * t443;
t152 = qJD(3) * t268 + t309 * t437;
t270 = Icges(6,6) * t443 - t440 * t517;
t153 = qJD(3) * t270 + t310 * t437;
t272 = Icges(6,5) * t443 - t440 * t520;
t154 = qJD(3) * t272 + t311 * t437;
t51 = (qJD(3) * t505 + t152) * t440 + (qJD(3) * t269 + (-t271 * t437 + t154) * t426 + (-t273 * t437 - t153) * t425) * t443;
t690 = t114 * t359 + t371 * t51;
t118 = t287 * t440 - t443 * t504;
t334 = (-Icges(5,5) * t439 - Icges(5,6) * t442) * t443;
t206 = qJD(3) * t286 + qJD(4) * t334;
t290 = Icges(5,6) * t443 - t440 * t518;
t209 = qJD(3) * t290 + qJD(4) * t337;
t294 = Icges(5,5) * t443 - t440 * t521;
t212 = qJD(3) * t294 + qJD(4) * t340;
t54 = (qJD(3) * t504 + t206) * t440 + (qJD(3) * t287 - t209 * t439 + t212 * t442 + (-t291 * t442 - t295 * t439) * qJD(4)) * t443;
t578 = qJD(3) * t607;
t689 = t118 * t578 + t54 * t754;
t594 = rSges(4,1) * t774 + rSges(4,2) * t587;
t219 = (-rSges(4,3) * qJD(1) - qJD(3) * t712) * t441 + t594;
t385 = rSges(4,1) * t443 - t712;
t350 = t385 * t610;
t540 = rSges(4,1) * t440 + rSges(4,2) * t443;
t366 = t540 * qJD(3);
t446 = qJD(1) ^ 2;
t604 = qJD(1) * qJD(2);
t613 = qJD(1) * t441;
t620 = qJ(2) * t611 + t423;
t626 = qJD(1) * (-pkin(1) * t613 + t620) + t441 * t604;
t496 = -t446 * t716 + t626;
t108 = t366 * t422 + (t219 + t350) * qJD(1) + t496;
t673 = t108 * t444;
t345 = t385 * t444;
t218 = -qJD(3) * t345 + (t441 * t540 + t434) * qJD(1);
t317 = qJD(1) * t386 - t424;
t415 = t444 * t604;
t548 = -t435 * t446 + t415;
t586 = t385 * t422;
t109 = -t366 * t610 + (-t218 - t317 + t586) * qJD(1) + t548;
t672 = t109 * t441;
t432 = t441 * rSges(4,3);
t301 = t444 * t540 - t432;
t142 = t350 + t423 + (t301 + t569) * qJD(1);
t670 = t142 * t444;
t388 = pkin(7) * t443 - t718;
t368 = qJD(3) * t388;
t665 = t368 * t444;
t664 = t375 * t444;
t335 = t441 * t375;
t659 = t441 * t389;
t274 = -t440 * t535 + t707;
t312 = (-rSges(6,1) * t425 - rSges(6,2) * t426) * t443;
t156 = qJD(3) * t274 + t312 * t437;
t266 = -t443 * t714 + t763;
t222 = qJD(3) * t266 - t443 * t600;
t649 = t156 + t222;
t644 = -t166 - t204;
t638 = -t201 - t346;
t637 = -t203 + t348;
t624 = t368 * t441 + t389 * t611;
t619 = rSges(3,2) * t613 + rSges(3,3) * t611;
t618 = -qJD(1) * t382 + t423;
t615 = qJD(1) * t516;
t500 = t377 * t443 + t379 * t440;
t181 = t444 * t500 - t335;
t605 = t181 * qJD(1);
t602 = -rSges(4,3) + t741;
t601 = t443 * t717;
t597 = -t346 + t644;
t596 = rSges(5,1) * t179 + rSges(5,2) * t178 + rSges(5,3) * t585;
t116 = t288 * t444 + t292 * t655 + t296 * t661;
t117 = -t289 * t444 - t293 * t655 - t297 * t661;
t590 = t443 * (-rSges(5,3) - pkin(7));
t580 = -t655 / 0.2e1;
t579 = t653 / 0.2e1;
t574 = -t610 / 0.2e1;
t572 = t609 / 0.2e1;
t188 = rSges(6,1) * t302 - rSges(6,2) * t303;
t189 = rSges(6,1) * t304 + rSges(6,2) * t305;
t565 = -t188 * t284 + t189 * t285;
t564 = t188 * t371 - t285 * t312;
t563 = -t189 * t371 + t284 * t312;
t349 = t389 * t444;
t562 = qJD(1) * t349 + t388 * t610;
t556 = t156 * t655 + t166 * t609 + t275 * t587 + t440 * t86;
t546 = qJD(4) * t572;
t544 = -t349 * t422 - t610 * t659;
t543 = qJD(1) * t659 - t388 * t422;
t539 = rSges(5,1) * t177 + rSges(5,2) * t176;
t65 = -t157 * t655 + t651;
t534 = t441 * t66 + t444 * t65;
t533 = t441 * t65 - t444 * t66;
t68 = t159 * t653 - t511;
t532 = t441 * t68 + t444 * t67;
t531 = t441 * t67 - t444 * t68;
t530 = t441 * t70 + t444 * t69;
t529 = t441 * t69 - t444 * t70;
t71 = -t190 * t655 + t650;
t528 = t441 * t72 + t444 * t71;
t527 = t441 * t71 - t444 * t72;
t74 = t192 * t653 - t508;
t526 = t441 * t74 + t444 * t73;
t525 = t441 * t73 - t444 * t74;
t524 = t441 * t78 + t444 * t77;
t523 = t441 * t77 - t444 * t78;
t143 = qJD(1) * t777 - t424 - t586;
t513 = t142 * t441 - t143 * t444;
t506 = t201 * t444 + t203 * t441;
t503 = t292 * t443 + t296 * t440;
t344 = (-rSges(5,1) * t439 - rSges(5,2) * t442) * t443;
t483 = t204 * t52 - t766;
t481 = (t116 * t444 + t117 * t441) * qJD(3);
t276 = t441 * t288;
t119 = -t444 * t503 + t276;
t120 = -t289 * t441 + t759;
t480 = (t119 * t444 + t120 * t441) * qJD(3);
t169 = (-t299 * t441 - t301 * t444) * qJD(3);
t479 = t496 + (-t476 + t352) * qJD(1);
t475 = -pkin(6) * t613 + qJD(1) * t348 + t352 + t618;
t473 = (Icges(6,5) * t302 - Icges(6,6) * t303) * t285 + (Icges(6,5) * t304 + Icges(6,6) * t305) * t284 + t309 * t371;
t472 = t190 * t361 + t192 * t360 + t287 * t754;
t471 = (Icges(5,5) * t327 - Icges(5,6) * t328) * t361 + (Icges(5,5) * t329 + Icges(5,6) * t330) * t360 + t334 * t754;
t465 = -qJD(1) * t501 - qJD(3) * t664 + t617;
t464 = qJD(1) * t503 + qJD(3) * t335 + t616;
t463 = t500 * qJD(1) - qJD(3) * t516;
t21 = t107 * t754 - t156 * t285 + t166 * t359 - t199 * t275 - t222 * t361 + t264 * t744 + t371 * t86 + (t204 * t607 - t665) * qJD(3) + t479;
t460 = t21 * (t166 * t440 + t275 * t655) + t52 * (t166 * t773 + t168 * t585);
t298 = -t440 * t538 + t709;
t459 = t368 * t610 + t389 * t413 + (-t235 - t317) * qJD(1) + t548;
t14 = t132 * t160 + t133 * t163 - t157 * t773 + t304 * t82 - t305 * t84 + t653 * t80;
t15 = t132 * t162 - t133 * t164 - t159 * t773 + t304 * t81 - t305 * t83 + t653 * t79;
t16 = t134 * t160 + t135 * t163 + t157 * t478 + t302 * t82 + t303 * t84 - t655 * t80;
t17 = t134 * t162 - t135 * t164 + t159 * t478 + t302 * t81 + t303 * t83 - t655 * t79;
t35 = t285 * t65 + t764;
t36 = t284 * t68 - t788;
t39 = t114 * t371 + t284 * t70 + t285 * t69;
t43 = t132 * t271 + t133 * t273 + t152 * t653 + t153 * t304 - t154 * t305 - t269 * t773;
t44 = t134 * t271 + t135 * t273 - t152 * t655 + t153 * t302 + t154 * t303 + t269 * t478;
t5 = t14 * t285 + t15 * t284 + t199 * t67 + t200 * t68 + t359 * t97 + t371 * t43;
t6 = t16 * t285 + t17 * t284 + t199 * t65 + t200 * t66 + t359 * t96 + t371 * t44;
t452 = ((qJD(3) * t531 + t43) * t440 + (-qJD(1) * t532 + qJD(3) * t97 - t14 * t441 + t15 * t444) * t443) * t735 + (t304 * t457 + t305 * t743 + t473 * t653) * t736 + (t302 * t457 - t303 * t743 - t473 * t655) * t734 + ((qJD(3) * t533 + t44) * t440 + (-qJD(1) * t534 + qJD(3) * t96 - t16 * t441 + t17 * t444) * t443) * t733 + (t473 * t440 + (-t425 * t457 - t426 * t743) * t443) * t727 + t6 * t580 + (t440 * t96 - t443 * t533) * t740 + (t440 * t97 - t443 * t531) * t739 + t39 * t572 + t5 * t579 + (t114 * t440 - t443 * t529) * t732 + ((qJD(3) * t529 + t51) * t440 + (-qJD(1) * t530 + qJD(3) * t114 - t18 * t441 + t19 * t444) * t443) * t726 + (t690 + t696 + t697 + t704 + t705) * t723 + t776 * t36 + t775 * t35;
t75 = -t201 * t360 + t203 * t361 + t545;
t93 = t300 * t360 + t462 - t779;
t449 = t75 * t506 + (-t441 * t94 - t444 * t93) * t300;
t383 = rSges(3,2) * t441 + t711;
t354 = t389 * t613;
t343 = t385 * t441;
t325 = t437 * t661;
t319 = t329 * pkin(4);
t318 = t327 * pkin(4);
t313 = t444 * t348;
t259 = t300 * t444;
t258 = t300 * t441;
t257 = t295 * t444;
t256 = t295 * t441;
t255 = t291 * t444;
t254 = t291 * t441;
t251 = qJD(1) * t757 - t424;
t250 = t423 + (-t382 + t383) * qJD(1);
t249 = t275 * t653;
t246 = t275 * t444;
t245 = t275 * t441;
t244 = t273 * t444;
t243 = t273 * t441;
t242 = t271 * t444;
t241 = t271 * t441;
t234 = t744 * t444;
t233 = t744 * t441;
t232 = rSges(5,1) * t329 + rSges(5,2) * t330;
t231 = rSges(5,1) * t327 - rSges(5,2) * t328;
t223 = t444 * t235;
t217 = qJD(3) * t298 + qJD(4) * t344;
t216 = t415 + (-qJD(1) * t566 - t317) * qJD(1);
t215 = qJD(1) * t619 + t626;
t180 = t441 * t500 + t664;
t171 = t180 * qJD(1);
t145 = t156 * t653;
t105 = -rSges(5,3) * t587 + t596;
t104 = -rSges(5,3) * t773 + t539;
t92 = -t441 * t749 + t444 * t463;
t91 = t441 * t463 + t444 * t749;
t90 = qJD(3) * t501 - t210 * t440 + t213 * t443;
t89 = -qJD(3) * t503 - t211 * t440 + t214 * t443;
t64 = t480 - t605;
t63 = t171 + t481;
t50 = t178 * t291 + t179 * t295 - t206 * t655 + t209 * t327 + t212 * t328 + t287 * t478;
t49 = t176 * t291 + t177 * t295 + t206 * t653 + t209 * t329 - t212 * t330 - t287 * t773;
t47 = -t104 * t754 - t203 * t578 + t217 * t360 + t265 * t300 + t459;
t46 = t105 * t754 - t217 * t361 - t264 * t300 + (t201 * t607 - t665) * qJD(3) + t479;
t45 = t118 * t754 + t360 * t78 + t361 * t77;
t41 = t360 * t74 - t789;
t40 = t361 * t71 + t758;
t37 = t104 * t361 - t105 * t360 - t201 * t265 + t203 * t264 + t453;
t26 = t100 * t327 + t102 * t328 + t178 * t195 - t179 * t197 + t192 * t478 - t655 * t98;
t25 = t101 * t327 + t103 * t328 + t178 * t193 + t179 * t196 + t190 * t478 - t655 * t99;
t24 = t100 * t329 - t102 * t330 + t176 * t195 - t177 * t197 - t192 * t773 + t653 * t98;
t23 = t101 * t329 - t103 * t330 + t176 * t193 + t177 * t196 - t190 * t773 + t653 * t99;
t22 = -t106 * t754 + t156 * t284 - t168 * t359 + t200 * t275 - t205 * t578 + t222 * t360 - t265 * t744 - t371 * t85 + t459;
t10 = t112 * t578 + t25 * t361 + t26 * t360 + t264 * t71 + t265 * t72 + t50 * t754;
t9 = t113 * t578 + t23 * t361 + t24 * t360 + t264 * t73 + t265 * t74 + t49 * t754;
t1 = [(t43 + t35) * t735 + (t49 + t40) * t730 + (t605 + (t289 * t441 ^ 2 + (-t276 + t117 + (t289 + t503) * t444) * t444) * qJD(3) + t64) * t571 - ((-t502 + t180) * t441 + t444 * t181) * t603 / 0.2e1 + (t90 + t91 + t63) * t573 + t697 / 0.2e1 + (t109 * (-t432 + t569) + t142 * t424 + t108 * t777 + t143 * (-rSges(4,2) * t585 + t594 + t620) + (qJD(3) * t142 * t385 + t109 * t540) * t444 + (t602 * t670 + (t142 * (-qJ(2) - t540) + t143 * t602) * t441) * qJD(1) - (-t142 + t350 + (t301 - t716) * qJD(1) + t618) * t143) * m(4) + t695 / 0.2e1 + t696 / 0.2e1 + t44 * t733 + t113 * t737 + t112 * t738 + t97 * t739 + t96 * t740 + (-(qJD(1) * t383 - t250 + t618) * t251 + t216 * (t441 * t719 + t428 + t711) + t250 * t424 + t215 * t757 + t251 * (t619 + t620) + (t250 * t719 * t444 + (t250 * (-rSges(3,3) - qJ(2)) - t251 * pkin(1)) * t441) * qJD(1)) * m(3) + (-qJD(3) * t500 + t363 * t440 - t364 * t443) * qJD(1) + t704 / 0.2e1 + t705 / 0.2e1 + (t171 + ((-t119 + t276 + t117) * t441 + (t120 - t759 + (t289 - t503) * t441 + t116) * t444) * qJD(3)) * t574 + t702 / 0.2e1 + t703 / 0.2e1 + t50 * t728 + t694 / 0.2e1 + (t47 * (t428 + t637) + t93 * (rSges(5,3) * t584 + t424 - t539 + t592) + t46 * (t410 + t756 + t625) + t94 * (t593 + t596 + t620) + (t46 * t590 + t47 * t741) * t441 + ((t590 * t94 + t741 * t93) * t444 + (t93 * (-qJ(2) + t709 - t718) + t94 * t741) * t441) * qJD(1) - (t475 - t93 - t779) * t94 - t693 * t360) * m(5) + t689 + t690 + (qJD(1) * t175 + t89 + t92) * t422 / 0.2e1 + ((-t441 * t707 + t595 + t633 + t756) * t21 + (t441 * t559 - t444 * t707 + t428 + t536 + t623) * t22 + (t424 - t537 + (-t557 + (t663 + (rSges(6,3) - t445) * t440) * qJD(3) + t784) * t444 + (-t599 + (-t418 * t440 - qJ(2) - t652 + t707) * qJD(1)) * t441) * t61 + (-t475 + t61 + t553 + t598 + t620 - t707 * t611 + ((-t600 - t608) * t440 + t784) * t441 - t785) * t62) * m(6) + (t36 + (-t65 + t790) * t284 + t788) * t734 + ((t68 + t790) * t285 + t764) * t736 + (t41 + (-t71 + t791) * t360 + t789) * t729 + ((t74 + t791) * t361 + t758) * t731; 0.2e1 * (t21 * t721 + t22 * t722) * m(6) + 0.2e1 * (t46 * t721 + t47 * t722) * m(5) + 0.2e1 * (-t673 / 0.2e1 + t672 / 0.2e1) * m(4) + 0.2e1 * (t215 * t721 + t216 * t722) * m(3); t524 * t546 - qJD(1) * ((-t440 * t622 - t443 * t621) * qJD(1) + (t440 * t466 + t443 * t467) * qJD(3)) / 0.2e1 + ((t335 * t422 - t615) * t444 + (-t771 + (-t444 * t664 - t772) * qJD(3)) * t441) * t571 + ((-t610 * t664 - t615) * t441 + (t771 + (t441 * t335 + t772) * qJD(3)) * t444) * t574 - t40 * t581 / 0.2e1 + t530 * t732 + (-qJD(1) * t533 + t16 * t444 + t17 * t441) * t733 + (-qJD(1) * t531 + t14 * t444 + t15 * t441) * t735 + t526 * t737 + t528 * t738 + t532 * t739 + t534 * t740 - t45 * t607 / 0.2e1 + (-(t142 * t345 + t143 * t343) * qJD(1) - (t169 * (-t343 * t441 - t345 * t444) - t513 * t540) * qJD(3) + 0.2e1 * t169 * (t218 * t444 - t219 * t441 + (-t299 * t444 + t301 * t441) * qJD(1)) - t513 * t366 + (-t673 + t672 + (t143 * t441 + t670) * qJD(1)) * t385) * m(4) + (qJD(1) * t91 + (t441 * (t441 * t465 - t444 * t750) + t444 * (t441 * t464 + t444 * t748) + (-t119 * t441 + t120 * t444) * qJD(1)) * t767 + t5 + t9) * t722 + (qJD(1) * t92 + t10 + (t441 * (t441 * t750 + t444 * t465) + t444 * (-t441 * t748 + t444 * t464) + (-t116 * t441 + t117 * t444) * qJD(1)) * t767 + t6) * t444 / 0.2e1 + (t22 * t659 + t61 * t624 + t62 * t354 - t11 * t313 + t52 * t223 + (t22 * t636 + t61 * t649 + t11 * t597 + t52 * (t476 + t691) + (t765 + t52 * (t348 - t643)) * qJD(1)) * t441 + (t21 * (-t389 - t636) + t62 * (-t368 - t649) + (t52 * t597 + t766) * qJD(1) + t747) * t444 - t61 * (-t168 * t373 - t234 * t754 + t246 * t371 + t266 * t360 + t274 * t284 - t275 * t326 + t562) - t62 * (t166 * t373 - t233 * t754 + t245 * t371 - t266 * t361 - t274 * t285 - t275 * t325 + t543) - t52 * (t166 * t326 + t168 * t325 + t233 * t360 + t234 * t361 - t245 * t284 - t246 * t285 + t544) - ((-t205 * t61 + t698) * t443 + (t52 * (t204 * t444 + t205 * t441) - (-t441 * t62 - t444 * t61) * t744) * t440) * qJD(4)) * m(6) + (t480 + t64 + t41 + t36) * t611 / 0.2e1 - (t481 + t63 + t40 + t35) * t613 / 0.2e1 + (-qJD(1) * t523 + t27 * t444 + t28 * t441) * t724 + (-qJD(1) * t529 + t18 * t444 + t19 * t441) * t726 + (-qJD(1) * t527 + t25 * t444 + t26 * t441) * t728 + (-qJD(1) * t525 + t23 * t444 + t24 * t441) * t730 + (((-t254 * t439 + t256 * t442 + t190) * t361 + (t255 * t439 - t257 * t442 + t192) * t360 + (-t290 * t439 + t294 * t442 + t287) * t754 + t118 * qJD(4)) * t443 + (qJD(4) * t523 + t447) * t440) * t725 + ((t254 * t327 + t256 * t328) * t361 + (-t255 * t327 - t257 * t328) * t360 + (t290 * t327 + t294 * t328) * t754 + (t112 * t443 - t660 * t72) * qJD(4) + ((qJD(4) * t71 + t472) * t440 - t762) * t441) * t729 + ((t254 * t329 - t256 * t330) * t361 + (-t255 * t329 + t257 * t330) * t360 + (t290 * t329 - t294 * t330) * t754 + (t113 * t443 + t661 * t73) * qJD(4) + ((-qJD(4) * t74 - t472) * t440 + t762) * t444) * t731 + (t47 * t659 + t93 * t624 + t94 * t354 - t37 * t313 + (t693 * qJD(1) + t93 * t217 + t47 * t300 + t37 * t638) * t441 + (t46 * (-t300 - t389) + t94 * (-t217 - t368) + t37 * t203 + t93 * t300 * qJD(1)) * t444 - t93 * (t259 * t754 + t298 * t360 + t562) - t94 * (t258 * t754 - t298 * t361 + t543) - ((t201 * t94 - t203 * t93) * t443 + t449 * t440) * qJD(4) + (t223 + (qJD(1) * t637 - t105 + t476) * t441 + (qJD(1) * t638 + t104) * t444 + t258 * t360 + t259 * t361 - t544) * t75) * m(5) + ((t241 * t304 - t243 * t305) * t285 + t67 * t325 + (-t242 * t304 + t244 * t305) * t284 - t68 * t326 + (t270 * t304 - t272 * t305) * t371 + t97 * t373 + t746 * t444) * t736 + (t114 * t373 + t69 * t325 - t70 * t326 + ((-t241 * t425 + t243 * t426 + t157) * t285 + (t242 * t425 - t244 * t426 + t159) * t284 + (-t270 * t425 + t272 * t426 + t269) * t371) * t443 + t448 * t440) * t727 + ((t241 * t302 + t243 * t303) * t285 + t65 * t325 + (-t242 * t302 - t244 * t303) * t284 - t66 * t326 + (t270 * t302 + t272 * t303) * t371 + t96 * t373 - t746 * t441) * t734 + qJD(1) * (t441 * t90 + t444 * t89 + (t175 * t444 + t441 * t502) * qJD(1)) / 0.2e1 + t41 * t606 * t723 - t325 * t35 / 0.2e1 + t326 * t36 / 0.2e1 - t373 * t39 / 0.2e1; (t471 * t440 + (-t439 * t455 - t442 * t742) * t443) * t725 + (t455 * t329 + t330 * t742 + t471 * t653) * t731 + ((qJD(3) * t523 + t54) * t440 + (-qJD(1) * t524 + qJD(3) * t118 - t27 * t441 + t28 * t444) * t443) * t724 + (t112 * t440 - t443 * t527) * t738 + t9 * t579 + t10 * t580 + t45 * t572 + ((qJD(3) * t525 + t49) * t440 + (-qJD(1) * t526 + qJD(3) * t113 - t23 * t441 + t24 * t444) * t443) * t730 + (t327 * t455 - t328 * t742 - t471 * t655) * t729 + ((qJD(3) * t527 + t50) * t440 + (-qJD(1) * t528 + qJD(3) * t112 - t25 * t441 + t26 * t444) * t443) * t728 + (t113 * t440 - t443 * t525) * t737 + t452 + (t118 * t440 - t443 * t523) * t546 + (t689 + t694 + t695 + t702 + t703) * t723 + t776 * t41 + t775 * t40 + (t22 * t249 + t61 * t145 + t62 * t556 + (-t22 * t643 - t61 * t692 + t21 * t204 + t62 * t107 + (t483 * t444 + (t205 * t52 - t765) * t441) * qJD(3)) * t440 + ((-t61 * t643 + t698) * qJD(3) + (-t22 * t744 + t61 * t222 + t11 * t644 + t52 * t691 + (-t52 * t643 - t62 * t744) * qJD(1)) * t444 + (qJD(1) * t483 - t21 * t744 + t222 * t62 - t747) * t441) * t443 + t460 - t61 * (-t319 * t754 - t360 * t601 + t563) - t62 * (t318 * t754 + t361 * t601 + t564) - t52 * (-t318 * t360 + t319 * t361 + t565)) * m(6) + ((qJD(3) * t449 - t104 * t93 + t105 * t94 + t201 * t46 - t203 * t47) * t440 + (t93 * (-qJD(3) * t203 + t217 * t444) + t94 * (qJD(3) * t201 + t217 * t441) - t37 * t506 + t75 * (-t104 * t441 - t105 * t444 + t201 * t613 - t203 * t611) + (t46 * t441 + t47 * t444 + (-t441 * t93 + t444 * t94) * qJD(1)) * t300) * t443 - t93 * (-t232 * t754 + t344 * t360) - t94 * (t231 * t754 - t344 * t361) - t75 * (-t231 * t360 + t232 * t361)) * m(5); t452 + (-t52 * t565 + t22 * (-t168 * t440 + t249) + (t11 * (-t166 * t444 - t168 * t441) + t52 * (-t168 * t611 - t441 * t85 - t444 * t86)) * t443 + t460 + (-t275 * t585 + t556 - t564) * t62 + (-t563 - t275 * t584 - t440 * t85 + t145 + (-qJD(3) * t168 - t275 * t613) * t443) * t61) * m(6);];
tauc = t1(:);
