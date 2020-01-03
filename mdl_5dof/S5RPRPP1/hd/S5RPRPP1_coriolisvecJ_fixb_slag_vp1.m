% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:39
% EndTime: 2019-12-31 18:09:06
% DurationCPUTime: 23.34s
% Computational Cost: add. (12142->566), mult. (11535->683), div. (0->0), fcn. (8829->8), ass. (0->325)
t667 = Icges(4,3) + Icges(5,3);
t309 = qJ(3) + pkin(8);
t303 = sin(t309);
t305 = cos(t309);
t312 = sin(qJ(3));
t314 = cos(qJ(3));
t666 = Icges(4,5) * t314 + Icges(5,5) * t305 - Icges(4,6) * t312 - Icges(5,6) * t303;
t213 = Icges(6,4) * t305 + Icges(6,6) * t303;
t310 = qJ(1) + pkin(7);
t304 = sin(t310);
t306 = cos(t310);
t147 = Icges(6,2) * t304 + t213 * t306;
t623 = t667 * t304 + t666 * t306;
t665 = t147 + t623;
t285 = Icges(6,5) * t303;
t384 = Icges(6,3) * t305 - t285;
t530 = Icges(5,4) * t303;
t662 = Icges(5,2) * t305 + t384 + t530;
t526 = Icges(6,5) * t305;
t216 = Icges(6,1) * t303 - t526;
t289 = Icges(5,4) * t305;
t661 = Icges(5,1) * t303 + t216 + t289;
t664 = t667 * t306;
t498 = t304 * t314;
t499 = t304 * t312;
t500 = t304 * t305;
t502 = t303 * t304;
t624 = -Icges(4,5) * t498 - Icges(5,5) * t500 + Icges(4,6) * t499 + Icges(5,6) * t502 + t664;
t387 = Icges(6,1) * t305 + t285;
t151 = Icges(6,4) * t304 + t306 * t387;
t219 = Icges(5,1) * t305 - t530;
t153 = Icges(5,5) * t304 + t219 * t306;
t655 = t151 + t153;
t209 = Icges(6,3) * t303 + t526;
t385 = -Icges(5,2) * t303 + t289;
t663 = -t209 + t385;
t660 = t219 + t387;
t521 = Icges(5,6) * t306;
t148 = Icges(5,4) * t500 - Icges(5,2) * t502 - t521;
t522 = Icges(4,6) * t306;
t160 = Icges(4,4) * t498 - Icges(4,2) * t499 - t522;
t659 = t148 * t303 + t160 * t312;
t645 = Icges(4,5) * t312 + Icges(4,6) * t314 + (Icges(5,6) - Icges(6,6)) * t305 + (Icges(6,4) + Icges(5,5)) * t303;
t531 = Icges(4,4) * t312;
t259 = Icges(4,1) * t314 - t531;
t163 = Icges(4,5) * t304 + t259 * t306;
t658 = -t153 * t500 - t163 * t498;
t657 = t213 + t666;
t243 = Icges(5,4) * t502;
t527 = Icges(5,5) * t306;
t152 = Icges(5,1) * t500 - t243 - t527;
t267 = Icges(4,4) * t499;
t528 = Icges(4,5) * t306;
t162 = Icges(4,1) * t498 - t267 - t528;
t656 = t152 * t305 + t162 * t314 - t659;
t256 = Icges(4,2) * t314 + t531;
t307 = Icges(4,4) * t314;
t258 = Icges(4,1) * t312 + t307;
t644 = t256 * t312 - t258 * t314 + t303 * t662 - t305 * t661;
t654 = t663 * qJD(3);
t653 = t660 * qJD(3);
t651 = t662 * qJD(3);
t650 = t661 * qJD(3);
t497 = t305 * t306;
t242 = Icges(6,5) * t497;
t501 = t303 * t306;
t520 = Icges(6,6) * t304;
t143 = Icges(6,3) * t501 + t242 + t520;
t395 = -t143 * t502 + t147 * t306 - t151 * t500;
t149 = Icges(5,6) * t304 + t306 * t385;
t386 = -Icges(4,2) * t312 + t307;
t161 = Icges(4,6) * t304 + t306 * t386;
t648 = -t623 * t306 - t658;
t616 = -t149 * t502 - t161 * t499 + t648;
t649 = -t395 + t616;
t495 = t306 * t314;
t647 = -t152 * t497 - t162 * t495 + t624 * t304;
t646 = t143 * t501 + t163 * t495 + t665 * t304 + t655 * t497;
t569 = t645 * t306;
t570 = t645 * t304;
t496 = t306 * t312;
t643 = t148 * t501 + t160 * t496 + t647;
t311 = -qJ(4) - pkin(6);
t272 = t306 * t311;
t549 = pkin(3) * t314;
t300 = pkin(2) + t549;
t463 = -t304 * t300 - t272;
t313 = sin(qJ(1));
t551 = pkin(1) * t313;
t642 = t463 - t551;
t641 = t149 * t303 + t161 * t312;
t146 = -Icges(6,2) * t306 + t213 * t304;
t514 = t146 * t306;
t142 = -Icges(6,6) * t306 + t209 * t304;
t150 = -Icges(6,4) * t306 + t304 * t387;
t382 = t142 * t303 + t150 * t305;
t582 = t304 * t382;
t48 = -t514 + t582;
t640 = t656 * t304 + t306 * t624 + t48;
t606 = -t149 * t501 - t161 * t496 + t646;
t639 = t644 * t304 + t569;
t638 = -t306 * t644 + t570;
t637 = t651 * t306 + (t304 * t385 - t142 - t521) * qJD(1);
t636 = t651 * t304 + (t209 * t306 - t149 + t520) * qJD(1);
t635 = -t650 * t306 + (-t219 * t304 - t150 + t527) * qJD(1);
t634 = -qJD(1) * t655 + t650 * t304;
t633 = -t142 + t148;
t632 = t143 - t149;
t631 = t150 + t152;
t435 = rSges(5,1) * t500;
t629 = -t435 + t642;
t233 = t386 * qJD(3);
t234 = t259 * qJD(3);
t628 = -t233 * t312 + t234 * t314 + t653 * t305 - t654 * t303 + (-t256 * t314 - t258 * t312 - t303 * t661 - t305 * t662) * qJD(3) + t645 * qJD(1);
t627 = t143 * t303 + t163 * t314 + t305 * t655 - t641;
t626 = -t382 - t656;
t625 = t644 * qJD(1) + qJD(3) * t657;
t315 = cos(qJ(1));
t308 = t315 * pkin(1);
t405 = t306 * t300 - t304 * t311;
t621 = t308 + t405;
t620 = t638 * qJD(1);
t132 = t304 * t146;
t52 = t142 * t501 + t150 * t497 + t132;
t593 = t306 * t52;
t619 = (t606 * t304 + t643 * t306 - t593) * qJD(3);
t618 = (t304 * t649 - t640 * t306) * qJD(3);
t617 = t639 * qJD(1);
t292 = t304 * rSges(6,2);
t614 = t624 + t641;
t445 = qJD(5) * t303;
t613 = -t617 + t618;
t612 = t619 + t620;
t353 = qJD(3) * t256;
t106 = qJD(1) * t161 - t304 * t353;
t356 = qJD(3) * t258;
t108 = qJD(1) * t163 - t304 * t356;
t611 = qJD(3) * t626 - t106 * t314 - t108 * t312 + t303 * t634 + t305 * t636;
t105 = -t306 * t353 + (-t304 * t386 + t522) * qJD(1);
t107 = -t306 * t356 + (-t259 * t304 + t528) * qJD(1);
t610 = t627 * qJD(3) + t105 * t314 + t107 * t312 + t635 * t303 - t305 * t637;
t609 = t304 * t625 + t306 * t628;
t608 = t304 * t628 - t306 * t625;
t607 = t52 - t643;
t605 = -t161 * t314 - t163 * t312 - t303 * t655 + t305 * t632;
t604 = t160 * t314 + t162 * t312 + t303 * t631 + t305 * t633;
t603 = t645 * qJD(3);
t602 = t514 + t646;
t601 = t304 * (-t306 * t662 + t655) - t306 * (-Icges(5,2) * t500 - t384 * t304 - t243 + t631);
t592 = rSges(6,3) + qJ(5);
t600 = t633 * t306 + (-Icges(6,1) * t501 + t216 * t306 + t242 + t632) * t304;
t599 = t661 + t663;
t598 = -t662 + t660;
t595 = t665 * qJD(1);
t317 = qJD(1) ^ 2;
t552 = rSges(6,1) + pkin(4);
t594 = rSges(4,2) * t312;
t220 = pkin(4) * t303 - qJ(5) * t305;
t221 = rSges(6,1) * t303 - rSges(6,3) * t305;
t467 = t220 + t221;
t591 = (qJD(3) * t467 - t445) * t304;
t225 = rSges(6,1) * t305 + rSges(6,3) * t303;
t590 = pkin(4) * t305 + qJ(5) * t303 + t225;
t589 = t605 * qJD(3) - t105 * t312 + t107 * t314 + t303 * t637 + t635 * t305 + t595;
t572 = qJD(1) * t146;
t588 = t624 * qJD(1) + t604 * qJD(3) + t106 * t312 - t108 * t314 - t303 * t636 + t305 * t634 - t572;
t587 = qJD(1) * t626 - t603 * t304 + t595;
t586 = -t572 - t603 * t306 + (-t304 * t666 - t627 + t664) * qJD(1);
t585 = 0.2e1 * qJD(3);
t298 = t306 * pkin(6);
t227 = pkin(2) * t304 - t298;
t138 = t227 + t463;
t297 = t304 * pkin(6);
t228 = t306 * pkin(2) + t297;
t139 = t405 - t228;
t448 = qJD(3) * t306;
t449 = qJD(3) * t304;
t423 = -t138 * t449 + t139 * t448 + qJD(2);
t444 = qJD(5) * t305;
t479 = t497 * t552 + t501 * t592 + t292;
t295 = t306 * rSges(6,2);
t480 = t304 * t590 - t295;
t23 = -t444 + (t304 * t480 + t306 * t479) * qJD(3) + t423;
t584 = qJD(3) * t23;
t583 = t303 * t552;
t134 = qJD(1) * t138;
t207 = qJD(1) * t227;
t581 = t134 - t207;
t291 = t304 * rSges(4,3);
t165 = rSges(4,1) * t495 - rSges(4,2) * t496 + t291;
t414 = t228 + t308;
t580 = t165 + t414;
t426 = t305 * t448;
t450 = qJD(1) * t306;
t578 = rSges(6,2) * t450 + t426 * t592;
t577 = -rSges(5,2) * t502 - t306 * rSges(5,3);
t410 = t306 * rSges(3,1) - rSges(3,2) * t304;
t576 = t304 ^ 2 + t306 ^ 2;
t575 = t308 + t410;
t574 = -t303 * t601 + t600 * t305;
t550 = pkin(3) * t312;
t400 = -t467 - t550;
t367 = t400 * qJD(3);
t238 = t306 * t445;
t277 = qJD(4) * t304;
t464 = t238 + t277;
t573 = t306 * t367 + t464;
t459 = t258 + t386;
t460 = -t256 + t259;
t568 = (-t303 * t599 + t305 * t598 - t312 * t459 + t314 * t460) * qJD(1);
t567 = t657 * qJD(1);
t476 = -Icges(4,2) * t498 + t162 - t267;
t478 = t258 * t304 + t160;
t556 = t312 * t476 + t314 * t478;
t433 = -t139 - t479;
t276 = pkin(6) * t450;
t447 = qJD(3) * t312;
t424 = t306 * t447;
t546 = pkin(2) - t300;
t100 = -pkin(3) * t424 - t276 + t277 + (t304 * t546 - t272) * qJD(1);
t451 = qJD(1) * t304;
t252 = t304 * pkin(3) * t447;
t461 = qJD(4) * t306 + t252;
t432 = t311 * t451 + t461;
t101 = (-t306 * t546 - t297) * qJD(1) - t432;
t436 = t101 * t449 + (t100 - t134) * t448;
t180 = t221 * t304;
t544 = (pkin(4) * t450 + qJ(5) * t449) * t305 + (qJ(5) * t450 + (-pkin(4) * qJD(3) + qJD(5)) * t304) * t303 - qJD(3) * t180 + (t225 * t306 + t292) * qJD(1);
t345 = -t303 * t448 - t305 * t451;
t431 = t303 * t451;
t545 = t345 * t552 - t431 * t592 + t238 + t578;
t1 = (t445 + t545 * t306 + t544 * t304 + (t304 * t433 + t306 * t480) * qJD(1)) * qJD(3) + t436;
t555 = m(6) * t1;
t554 = t304 / 0.2e1;
t553 = -t306 / 0.2e1;
t547 = qJD(1) / 0.2e1;
t543 = rSges(4,1) * t314;
t261 = rSges(4,1) * t312 + rSges(4,2) * t314;
t196 = t261 * t306;
t428 = t261 * t449;
t74 = qJD(1) * t580 - t428;
t541 = t196 * t74;
t290 = t304 * rSges(5,3);
t458 = rSges(4,2) * t499 + t306 * rSges(4,3);
t164 = rSges(4,1) * t498 - t458;
t415 = -t227 - t551;
t427 = t261 * t448;
t73 = -t427 + (-t164 + t415) * qJD(1);
t540 = t304 * t73;
t539 = t306 * t73;
t222 = rSges(5,1) * t303 + rSges(5,2) * t305;
t155 = t435 + t577;
t413 = -t222 - t550;
t366 = t413 * t448;
t344 = t277 + t366;
t46 = (t138 - t155 + t415) * qJD(1) + t344;
t537 = t46 * t222;
t206 = t228 * qJD(1);
t492 = -t101 - t206;
t488 = -t304 * t138 + t306 * t139;
t157 = rSges(5,1) * t497 - rSges(5,2) * t501 + t290;
t485 = -t139 - t157;
t477 = -t258 * t306 - t161;
t475 = -t256 * t306 + t163;
t474 = -qJD(3) * t590 + t444;
t473 = -t220 * t304 - t180;
t472 = t467 * t306;
t465 = rSges(5,2) * t431 + rSges(5,3) * t450;
t462 = rSges(4,3) * t450 + t451 * t594;
t446 = qJD(3) * t314;
t442 = qJD(1) * qJD(4);
t441 = pkin(3) * t496;
t440 = qJD(3) ^ 2 * t549;
t439 = t317 * t551;
t438 = t317 * t308;
t437 = t306 * t100 + t304 * t101 - t138 * t450;
t434 = pkin(3) * t446;
t429 = t222 * t449;
t422 = -pkin(2) - t543;
t419 = -t449 / 0.2e1;
t416 = t448 / 0.2e1;
t226 = rSges(5,1) * t305 - rSges(5,2) * t303;
t412 = -t226 - t549;
t404 = t576 * t550;
t403 = qJD(1) * (-pkin(2) * t451 + t276) - t439;
t402 = t139 + t414;
t401 = -t480 - t551;
t399 = -t549 - t590;
t205 = t226 * qJD(3);
t396 = -t205 - t434;
t47 = -t429 + (t157 + t402) * qJD(1) - t461;
t394 = t47 * t413;
t223 = rSges(3,1) * t304 + rSges(3,2) * t306;
t391 = t543 - t594;
t390 = -t304 * t74 - t539;
t375 = t164 * t304 + t165 * t306;
t371 = qJD(1) * t252 + t306 * t442 - t438;
t370 = -t434 + t474;
t365 = -qJD(3) * t205 - t440;
t364 = qJD(1) * t100 + t304 * t442 + t403;
t195 = t261 * t304;
t181 = t222 * t304;
t341 = -t312 * t475 + t314 * t477;
t337 = -t440 + (t444 + t474) * qJD(3);
t336 = -t300 - t590;
t111 = -rSges(4,2) * t306 * t446 + (-t314 * t451 - t424) * rSges(4,1) + t462;
t112 = -qJD(3) * t195 + (t306 * t391 + t291) * qJD(1);
t335 = t111 * t306 + t112 * t304 + (t164 * t306 - t165 * t304) * qJD(1);
t235 = t391 * qJD(3);
t185 = t222 * t306;
t96 = -qJD(3) * t181 + (t226 * t306 + t290) * qJD(1);
t94 = rSges(5,1) * t345 - rSges(5,2) * t426 + t465;
t72 = qJD(3) * t375 + qJD(2);
t57 = -t438 - t235 * t448 + (-t112 - t206 + t428) * qJD(1);
t56 = -t235 * t449 + (t111 - t427) * qJD(1) + t403;
t41 = (t155 * t304 + t157 * t306) * qJD(3) + t423;
t38 = t335 * qJD(3);
t37 = -t591 + (t402 + t479) * qJD(1) - t461;
t36 = (t138 - t227 + t401) * qJD(1) + t573;
t25 = t365 * t306 + (-t96 + t429 + t492) * qJD(1) + t371;
t24 = t365 * t304 + (t94 + t366) * qJD(1) + t364;
t12 = t337 * t306 + (t492 - t544 + t591) * qJD(1) + t371;
t11 = t337 * t304 + ((t367 + t445) * t306 + t545) * qJD(1) + t364;
t2 = (t304 * t96 + t306 * t94 + (t155 * t306 + t304 * t485) * qJD(1)) * qJD(3) + t436;
t3 = [m(3) * ((-t223 * t317 - t439) * t575 + (-t438 + (-0.2e1 * t410 - t308 + t575) * t317) * (-t223 - t551)) + (-t644 * qJD(3) + t233 * t314 + t234 * t312 + t653 * t303 + t654 * t305) * qJD(1) + ((t479 + t621) * t11 + (t295 + (-t303 * t592 - t552 * t305) * t304 + t642) * t12 + (t432 + (-t445 + (-t305 * t592 + t583) * qJD(3)) * t304 + (t336 * t306 - t292 - t308) * qJD(1)) * t36 + (t36 - t573 - t581 + t464 + t578 + (-t550 - t583) * t448 + (t336 * t304 - t272 - t401 - t551) * qJD(1)) * t37) * m(6) + (t25 * (-t577 + t629) + t24 * (t157 + t621) + (t304 * t537 + t306 * t394) * qJD(3) + (t432 + (-t290 - t308 + (-t226 - t300) * t306) * qJD(1)) * t46 + (t46 - t344 - t581 + t277 + t465 + (t155 + t551 + t629) * qJD(1)) * t47) * m(5) + (-(-t427 - t207 - t73 + (-t164 - t551) * qJD(1)) * t74 + t57 * (t304 * t422 + t298 + t458 - t551) + t56 * t580 + t74 * (t276 + t462) + (t261 * t540 - t541) * qJD(3) + ((-t313 * t74 - t315 * t73) * pkin(1) + (-pkin(2) - t391) * t539 + (t73 * (-rSges(4,3) - pkin(6)) + t74 * t422) * t304) * qJD(1)) * m(4) + ((-t593 + ((t623 + t659) * t306 + t616 + t647 + t658) * t306 + (t48 - t582 + t602) * t304) * qJD(3) + t620) * t416 + (((t614 * t306 - t602 + t606) * t306 + (t304 * t614 - t132 + t395 + t607 - t648) * t304) * qJD(3) + t613 + t617) * t419 + (t609 + t610) * t449 / 0.2e1 - (t608 - t611 + t612) * t448 / 0.2e1 + ((t604 - t639) * t304 + (-t605 + t638) * t306) * qJD(3) * t547; m(4) * t38 + m(5) * t2 + t555; (t38 * t375 + t72 * t335 + t390 * t235 + (-t56 * t304 - t57 * t306 + (-t306 * t74 + t540) * qJD(1)) * t261 - (t195 * t73 - t541) * qJD(1) - (t72 * (-t195 * t304 - t196 * t306) + t390 * t391) * qJD(3)) * m(4) - (((t304 * t475 - t306 * t476) * t314 + (t304 * t477 + t306 * t478) * t312 + t601 * t305 + t600 * t303) * qJD(3) + (t303 * t598 + t305 * t599 + t312 * t460 + t314 * t459) * qJD(1)) * qJD(1) / 0.2e1 + (t611 * t306 + t610 * t304 + (t304 * t604 - t605 * t306) * qJD(1)) * t547 + ((-t569 * t449 + t567) * t304 + ((t556 * t306 + (t341 + t570) * t304 + t574) * qJD(3) + t568) * t306) * t419 + ((-t570 * t448 - t567) * t306 + ((t341 * t304 + (t556 + t569) * t306 + t574) * qJD(3) + t568) * t304) * t416 + (-(t23 * t303 + (t304 * t37 + t306 * t36) * t305) * qJD(5) - (-t36 * t473 + t37 * (-t441 - t472)) * qJD(1) - (-t23 * t404 + (-t23 * t472 + t36 * t399) * t306 + (t23 * t473 + t37 * t399) * t304) * qJD(3) + t1 * t488 + t23 * t437 + (t12 * t400 + t36 * t370 + t1 * t479 + t23 * t545 + (t23 * t480 + t37 * t400) * qJD(1)) * t306 + (t11 * t400 + t37 * t370 + t1 * t480 + t23 * t544 + (t23 * t433 + t36 * t467) * qJD(1)) * t304) * m(6) + (-(t46 * t181 + t47 * (-t185 - t441)) * qJD(1) - (-t41 * t404 + (-t41 * t185 + t412 * t46) * t306 + (-t41 * t181 + t412 * t47) * t304) * qJD(3) + t2 * t488 + t41 * t437 + (t25 * t413 + t46 * t396 + t2 * t157 + t41 * t94 + (t41 * t155 + t394) * qJD(1)) * t306 + (t24 * t413 + t47 * t396 + t2 * t155 + t41 * t96 + (t41 * t485 + t537) * qJD(1)) * t304) * m(5) + (t609 * qJD(1) + ((qJD(1) * t606 + t588 * t306) * t306 + (t586 * t304 + t607 * qJD(1) + (-t587 + t589) * t306) * t304) * t585) * t554 + (t608 * qJD(1) + ((t649 * qJD(1) + t587 * t306) * t306 + (t589 * t304 + t640 * qJD(1) + (-t586 + t588) * t306) * t304) * t585) * t553 + (t613 + t618) * t451 / 0.2e1 + (t612 + t619) * t450 / 0.2e1; 0.2e1 * (t11 * t553 + t12 * t554) * m(6) + 0.2e1 * (t24 * t553 + t25 * t554) * m(5); -t305 * t555 + 0.2e1 * (m(6) * (t11 * t304 + t12 * t306 + t584) / 0.2e1 - m(6) * t576 * t584 / 0.2e1) * t303;];
tauc = t3(:);
