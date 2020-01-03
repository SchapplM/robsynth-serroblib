% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:30
% EndTime: 2019-12-31 20:51:55
% DurationCPUTime: 20.61s
% Computational Cost: add. (13272->571), mult. (14870->667), div. (0->0), fcn. (11494->6), ass. (0->334)
t737 = Icges(5,1) + Icges(6,1);
t736 = Icges(6,2) + Icges(5,3);
t371 = sin(qJ(3));
t366 = Icges(6,4) * t371;
t373 = cos(qJ(3));
t465 = Icges(6,1) * t373 + t366;
t365 = Icges(5,5) * t371;
t466 = Icges(5,1) * t373 + t365;
t735 = -t465 - t466;
t734 = Icges(5,4) - Icges(6,5);
t733 = Icges(5,6) - Icges(6,6);
t732 = Icges(4,5) + t734;
t731 = Icges(4,6) - t733;
t370 = qJ(1) + qJ(2);
t363 = sin(t370);
t606 = Icges(4,4) * t371;
t312 = Icges(4,1) * t373 - t606;
t364 = cos(t370);
t433 = t312 * t364;
t197 = Icges(4,5) * t363 + t433;
t729 = t735 * t364;
t725 = t734 * t363 - t729;
t706 = t197 + t725;
t701 = -t736 * t373 + t365 + t366;
t720 = -Icges(4,2) * t373 - t606 + t701;
t575 = t364 * t373;
t730 = (Icges(6,4) + Icges(5,5)) * t575;
t367 = Icges(4,4) * t373;
t464 = -Icges(4,2) * t371 + t367;
t604 = Icges(5,5) * t373;
t298 = Icges(5,3) * t371 + t604;
t605 = Icges(6,4) * t373;
t302 = Icges(6,2) * t371 + t605;
t712 = t298 + t302;
t705 = -t464 + t712;
t311 = Icges(4,1) * t371 + t367;
t699 = t737 * t371 + t311 - t604 - t605;
t182 = -Icges(5,6) * t364 + t298 * t363;
t186 = Icges(6,6) * t364 + t302 * t363;
t728 = -t182 - t186;
t576 = t364 * t371;
t726 = t733 * t363 + t736 * t576 + t730;
t192 = Icges(6,5) * t364 + t363 * t465;
t194 = -Icges(5,4) * t364 + t363 * t466;
t579 = t363 * t371;
t336 = Icges(4,4) * t579;
t578 = t363 * t373;
t196 = Icges(4,1) * t578 - Icges(4,5) * t364 - t336;
t707 = t192 + t194 + t196;
t715 = t732 * t371 + t731 * t373;
t727 = t312 - t735;
t724 = -t371 * t720 - t373 * t699;
t184 = Icges(4,5) * t578 - Icges(4,6) * t579 - Icges(4,3) * t364;
t190 = Icges(4,4) * t578 - Icges(4,2) * t579 - Icges(4,6) * t364;
t592 = t190 * t371;
t453 = -t196 * t373 + t592;
t456 = t186 * t371 + t192 * t373;
t304 = Icges(5,4) * t373 + Icges(5,6) * t371;
t188 = -Icges(5,2) * t364 + t304 * t363;
t593 = t188 * t364;
t296 = Icges(6,5) * t373 + Icges(6,6) * t371;
t180 = Icges(6,3) * t364 + t296 * t363;
t595 = t180 * t364;
t460 = t182 * t371 + t194 * t373;
t660 = t363 * t460;
t681 = t363 * t456 - t593 + t595 + t660;
t723 = -t184 * t364 - t363 * t453 + t681;
t427 = t296 * t364;
t181 = -Icges(6,3) * t363 + t427;
t174 = t364 * t181;
t429 = t304 * t364;
t189 = Icges(5,2) * t363 + t429;
t680 = -t189 * t364 + t725 * t578 + t726 * t579 + t174;
t430 = t464 * t364;
t191 = Icges(4,6) * t363 + t430;
t157 = t197 * t578;
t300 = Icges(4,5) * t373 - Icges(4,6) * t371;
t428 = t300 * t364;
t185 = Icges(4,3) * t363 + t428;
t487 = t185 * t364 - t157;
t80 = -t191 * t579 - t487;
t722 = t80 + t680;
t719 = t726 * t576 + (t185 + t189) * t363 + t706 * t575;
t172 = t363 * t188;
t718 = -t363 * t184 - t707 * t575 + t728 * t576 - t172;
t717 = t705 * qJD(3);
t716 = t727 * qJD(3);
t700 = t296 - t300 - t304;
t649 = t715 * t364;
t650 = t715 * t363;
t369 = qJD(1) + qJD(2);
t713 = (-t699 * qJD(3) + t732 * t369) * t371 - (-t720 * qJD(3) - t731 * t369) * t373;
t596 = t180 * t363;
t671 = -t190 * t576 - t596 - t718;
t670 = -t181 * t363 - t191 * t576 + t719;
t711 = t724 * t363 + t649;
t710 = -t724 * t364 + t650;
t709 = t190 + t728;
t708 = t191 - t726;
t698 = t716 * t373 + t717 * t371 + t715 * t369 + (-t371 * t699 + t373 * t720) * qJD(3);
t697 = (Icges(5,2) + Icges(4,3) + Icges(6,3)) * t369 - t715 * qJD(3);
t591 = t191 * t371;
t696 = -t371 * t726 - t373 * t706 + t591;
t695 = t453 - t456 - t460;
t694 = qJD(3) * t700 - t369 * t724;
t693 = rSges(6,1) + pkin(4);
t692 = rSges(6,3) + qJ(5);
t577 = t364 * t369;
t691 = ((t428 + t429 - t427 + t695) * t369 + t697 * t363) * t364;
t690 = t710 * t369;
t689 = (t670 * t363 - t671 * t364) * qJD(3);
t688 = (t722 * t363 - t364 * t723) * qJD(3);
t687 = t699 - t705;
t686 = t727 + t720;
t685 = t364 * t720 + t706;
t684 = -t311 * t364 - t737 * t576 - t708 + t730;
t683 = t711 * t369;
t682 = t363 ^ 2;
t526 = qJD(3) * t371;
t502 = t364 * t526;
t573 = t369 * t373;
t422 = -t363 * t573 - t502;
t574 = t369 * t371;
t515 = t363 * t574;
t525 = qJD(3) * t373;
t501 = t364 * t525;
t285 = qJ(4) * t501;
t362 = qJD(4) * t371;
t329 = t364 * t362;
t544 = t285 + t329;
t140 = pkin(3) * t422 - qJ(4) * t515 + t544;
t600 = qJ(4) * t371;
t323 = pkin(3) * t373 + t600;
t247 = t323 * t363;
t215 = t369 * t247;
t679 = t140 + t215;
t218 = t501 - t515;
t504 = t363 * t525;
t217 = t364 * t574 + t504;
t372 = sin(qJ(1));
t615 = pkin(1) * qJD(1);
t516 = t372 * t615;
t262 = rSges(3,1) * t363 + rSges(3,2) * t364;
t590 = t262 * t369;
t211 = -t516 - t590;
t617 = rSges(6,2) * t371;
t324 = rSges(6,1) * t373 + t617;
t204 = rSges(6,3) * t364 + t324 * t363;
t523 = qJD(5) * t363;
t678 = t369 * t204 + t523;
t677 = -t683 + t688;
t676 = t689 + t690;
t675 = t695 * qJD(3) + t712 * t373 * t577 + (-t430 * t373 + (-t433 + t729) * t371) * t369 - t713 * t363;
t580 = t363 * t369;
t674 = -t696 * qJD(3) + (-t371 * t727 + t705 * t373) * t580 + t713 * t364;
t673 = -t694 * t363 + t698 * t364;
t672 = t698 * t363 + t694 * t364;
t667 = t593 + t719;
t666 = t706 * t371 + t708 * t373;
t665 = t707 * t371 + t709 * t373;
t663 = t697 * t364 + t696 * t369 + t700 * t580;
t662 = 0.2e1 * qJD(3);
t505 = t363 * t526;
t287 = rSges(5,1) * t505;
t325 = rSges(5,1) * t373 + rSges(5,3) * t371;
t352 = t363 * rSges(5,2);
t138 = rSges(5,3) * t504 - t287 + (t325 * t364 + t352) * t369;
t320 = rSges(5,1) * t371 - rSges(5,3) * t373;
t318 = pkin(3) * t371 - qJ(4) * t373;
t528 = qJD(3) * t363;
t256 = t318 * t528;
t521 = qJD(3) * qJD(4);
t496 = t373 * t521;
t374 = cos(qJ(1));
t368 = t374 * pkin(1);
t376 = qJD(1) ^ 2;
t518 = t376 * t368;
t448 = t369 * t256 + t364 * t496 - t518;
t527 = qJD(3) * t364;
t524 = qJD(4) * t373;
t255 = qJD(3) * t323 - t524;
t545 = -t325 * qJD(3) - t255;
t293 = pkin(3) * t505;
t500 = t363 * t362;
t513 = t364 * t573;
t141 = pkin(3) * t513 + qJ(4) * t217 - t293 + t500;
t265 = t364 * pkin(2) + t363 * pkin(7);
t222 = t265 * t369;
t570 = -t141 - t222;
t37 = t545 * t527 + (-t138 + (qJD(3) * t320 - t362) * t363 + t570) * t369 + t448;
t532 = -t323 - t325;
t611 = -rSges(5,3) - qJ(4);
t354 = t364 * rSges(5,2);
t205 = t325 * t363 - t354;
t533 = -t318 - t320;
t488 = t364 * t533;
t443 = qJD(3) * t488 + t329;
t411 = t443 - t516;
t358 = t364 * pkin(7);
t264 = pkin(2) * t363 - t358;
t547 = -t247 - t264;
t73 = (-t205 + t547) * t369 + t411;
t613 = t369 * t73;
t622 = -rSges(5,1) - pkin(3);
t473 = -t256 + t500;
t444 = -t320 * t528 + t473;
t517 = t374 * t615;
t208 = rSges(5,1) * t575 + rSges(5,3) * t576 + t352;
t252 = pkin(3) * t575 + qJ(4) * t576;
t482 = t252 + t265;
t653 = t482 + t208;
t74 = t369 * t653 + t444 + t517;
t661 = (t74 * t622 * t526 + (t371 * t611 + t373 * t622 - pkin(2)) * t613) * t364 + (-t37 * pkin(2) + (-t73 * qJD(4) + t37 * t611) * t371 + (qJD(3) * t611 * t73 + t37 * t622) * t373 + (t73 * (-rSges(5,2) - pkin(7)) + t74 * (-pkin(2) + t532)) * t369) * t363;
t254 = t369 * t264;
t658 = t215 + t254;
t260 = pkin(4) * t578 + qJ(5) * t364;
t291 = rSges(6,2) * t501;
t657 = t369 * t260 + t291 - t329 + t658 + t678;
t656 = rSges(4,1) * t575 + t363 * rSges(4,3);
t655 = rSges(6,2) * t576 + t693 * t575;
t319 = rSges(6,1) * t371 - rSges(6,2) * t373;
t620 = pkin(4) * t371;
t489 = t319 + t620;
t477 = -t318 - t489;
t426 = t477 * t527;
t398 = t426 - t523;
t546 = -t324 * qJD(3) - t255;
t619 = pkin(4) * t373;
t419 = (-t619 * qJD(3) + t546) * qJD(3);
t317 = pkin(7) * t577;
t621 = pkin(1) * t372;
t519 = t376 * t621;
t478 = t369 * (-pkin(2) * t580 + t317) - t519;
t421 = t363 * t496 + t478 + (t140 + t329) * t369;
t572 = -rSges(6,1) * t502 + pkin(4) * t422 - qJ(5) * t577 + t291 - t678;
t21 = t419 * t363 + (t398 + t572) * t369 + t421;
t522 = qJD(5) * t364;
t292 = pkin(4) * t505;
t652 = rSges(6,1) * t505 + t692 * t580 + t292;
t571 = rSges(6,1) * t513 + rSges(6,2) * t217 + (pkin(4) * t573 + qJD(5)) * t364 - t652;
t22 = t419 * t364 + (-t522 + (qJD(3) * t489 - t362) * t363 + t570 - t571) * t369 + t448;
t520 = -pkin(3) - t693;
t612 = -rSges(6,2) - qJ(4);
t631 = -qJD(5) + (t520 * t373 - pkin(2) - t600 - t617) * t369;
t475 = t329 - t516;
t552 = t204 + t260;
t64 = (t547 - t552) * t369 + t398 + t475;
t551 = -t363 * t692 + t655;
t510 = -t252 - t551;
t633 = -t319 * t528 + t369 * (t265 - t510) - t292 + t473 + t522;
t65 = t517 + t633;
t654 = (t65 * t520 * t526 - (t65 * t369 + t22) * t692 + t631 * t64) * t364 + (-t22 * pkin(2) - t21 * t692 + (-t64 * qJD(4) + t22 * t612) * t371 + (qJD(3) * t612 * t64 + t22 * t520) * t373 - t64 * pkin(7) * t369 + t631 * t65) * t363 - t426 * t65;
t531 = rSges(4,2) * t579 + t364 * rSges(4,3);
t206 = rSges(4,1) * t578 - t531;
t177 = t369 * t206;
t434 = -rSges(4,2) * t218 + rSges(4,3) * t577;
t651 = -rSges(4,1) * t502 + t177 + t254 + t317 + t434;
t648 = (-t687 * t371 + t686 * t373) * t369;
t647 = -t685 * t371 + t684 * t373;
t646 = -Icges(4,2) * t578 + t701 * t363 - t336 + t707;
t645 = t699 * t363 + t709;
t644 = t700 * t369;
t321 = rSges(4,1) * t371 + rSges(4,2) * t373;
t259 = t321 * t528;
t209 = -rSges(4,2) * t576 + t656;
t420 = t209 + t265;
t637 = t369 * t420 - t259;
t176 = t369 * t205;
t541 = rSges(5,2) * t577 + rSges(5,3) * t501;
t632 = t176 + t541 + t658;
t630 = t371 * t646 + t373 * t645;
t629 = (t552 * t369 + t572) * t364 + (t510 * t369 + t571) * t363;
t628 = m(5) / 0.2e1;
t627 = m(6) / 0.2e1;
t623 = t369 / 0.2e1;
t618 = rSges(4,1) * t373;
t506 = t321 * t527;
t423 = -t506 - t516;
t93 = (-t206 - t264) * t369 + t423;
t614 = t364 * t93;
t550 = -t208 - t252;
t549 = t363 * t247 + t364 * t252;
t248 = t318 * t364;
t548 = -t369 * t248 + t363 * t524;
t542 = t287 + t293;
t529 = t364 ^ 2 + t682;
t512 = t363 * t141 + t364 * t679;
t243 = t318 * t363;
t511 = -t243 * t528 - t248 * t527 + t362;
t508 = rSges(4,1) * t505 + rSges(4,2) * t217;
t507 = t317 + t544;
t497 = -pkin(2) - t618;
t495 = -t528 / 0.2e1;
t492 = t527 / 0.2e1;
t490 = t358 - t621;
t263 = t364 * rSges(3,1) - rSges(3,2) * t363;
t486 = -t184 + t591;
t484 = t141 * t528 + t371 * t521 + t527 * t679;
t220 = rSges(3,1) * t577 - rSges(3,2) * t580;
t472 = t247 * t528 + t252 * t527 - t524;
t470 = -rSges(4,2) * t371 + t618;
t94 = t517 + t637;
t469 = -t363 * t94 - t614;
t468 = t293 + t652;
t447 = -pkin(4) * t525 + t546;
t442 = t285 + t317 + t475;
t98 = (t206 * t363 + t209 * t364) * qJD(3);
t425 = t482 + t655;
t414 = t363 * t497 + t358 + t531;
t380 = (t497 * t614 + (t93 * (-rSges(4,3) - pkin(7)) + t94 * t497) * t363) * t369;
t379 = (((t80 - t157 + (t185 + t592) * t364 + t718) * t364 + (-t660 + (-t181 - t456) * t363 + t667 + t681) * t363) * qJD(3) + t690) * t492 + (-qJD(3) * t724 + t371 * t716 - t373 * t717) * t369 + (((t364 * t486 + t595 - t667 + t670) * t364 + (t363 * t486 - t172 + t174 + t487 + t596 + t671 - t680) * t363) * qJD(3) + t677 + t683) * t495 + (t673 + t674) * t528 / 0.2e1 - (t672 - t675 + t676) * t527 / 0.2e1 + ((t665 - t711) * t363 + (t666 + t710) * t364) * qJD(3) * t623;
t330 = t364 * t524;
t280 = t470 * qJD(3);
t253 = t318 * t580;
t251 = t321 * t364;
t250 = t320 * t364;
t249 = t319 * t364;
t246 = t321 * t363;
t245 = t320 * t363;
t224 = t529 * t526;
t212 = t263 * t369 + t517;
t169 = -t220 * t369 - t518;
t168 = -t369 * t590 - t519;
t139 = t369 * t656 - t508;
t136 = rSges(4,1) * t422 + t434;
t135 = rSges(5,1) * t422 - rSges(5,3) * t515 + t541;
t72 = (t205 * t363 + t208 * t364) * qJD(3) + t472;
t70 = -t518 - t280 * t527 + (-t139 - t222 + t259) * t369;
t69 = t136 * t369 + (-t280 * t363 - t321 * t577) * qJD(3) + t478;
t57 = (t552 * t363 + t551 * t364) * qJD(3) + t472;
t36 = t135 * t369 + (t545 * t363 + t369 * t488) * qJD(3) + t421;
t8 = ((t135 + t176) * t364 + (t550 * t369 + t138) * t363) * qJD(3) + t484;
t7 = qJD(3) * t629 + t484;
t1 = [m(3) * (t169 * (-t262 - t621) + t168 * (t263 + t368) + (-t220 - t517 + t212) * t211) + t379 + (t22 * t490 + t64 * (t468 - t517) + t21 * (t368 + t425) + (t442 + t64 + t516 + t657) * t65 + t654) * m(6) + (t37 * (t354 + t490) + t73 * (-t517 + t542) + t36 * (t368 + t653) + (t442 + t73 - t411 + t632) * t74 + t661) * m(5) + (t70 * (t414 - t621) + t93 * (t508 - t517) + t69 * (t368 + t420) + t380 + (-t516 - t423 + t93 + t651) * t94) * m(4); t379 + (t21 * t425 + t22 * t358 + (t507 + t657) * t65 + (t468 + t633) * t64 + t654) * m(6) + (t37 * (t354 + t358) + (t507 - t443 + t632) * t74 + (t542 + t444) * t73 + t661 + (t36 + t613) * t653) * m(5) + (t70 * t414 + t69 * t420 + t380 + (t506 + t651) * t94 + (t508 + t637) * t93) * m(4) + (t168 * t263 - t169 * t262 - t211 * t220 - t212 * t590 - (-t211 * t263 - t212 * t262) * t369) * m(3); (t7 * t549 + (t21 * t477 + t7 * t552) * t363 + (t22 * t477 + t7 * t551) * t364 + (t447 * t363 - t548 + (pkin(4) * t576 + t477 * t364 + t249) * t369) * t65 + (-t243 * t369 + t447 * t364 + t253 - t330) * t64 - (t363 * t65 + t364 * t64) * qJD(3) * (-t323 - t324 - t619) + (-t511 - (-t249 * t364 - t319 * t682 - t529 * t620) * qJD(3) + t512 + t629) * t57) * m(6) + (-t73 * (t330 + (t243 + t245) * t369) - t74 * (-t250 * t369 + t548) - t72 * t511 - ((-t72 * t250 + t532 * t73) * t364 + (-t72 * t245 + t532 * t74) * t363) * qJD(3) + t73 * t253 + t8 * t549 + t72 * t512 + (t37 * t533 + t73 * t545 + t8 * t208 + t72 * t135 + (t72 * t205 + t533 * t74) * t369) * t364 + (t36 * t533 + t74 * t545 + t8 * t205 + t72 * t138 + (t73 * t320 + t550 * t72) * t369) * t363) * m(5) + (0.2e1 * t98 * ((t136 + t177) * t364 + (-t209 * t369 + t139) * t363) + t469 * t280 + ((-t369 * t94 - t70) * t364 + (t369 * t93 - t69) * t363) * t321 - (t246 * t93 - t251 * t94) * t369 - (t98 * (-t246 * t363 - t251 * t364) + t469 * t470) * qJD(3)) * m(4) - ((t686 * t371 + t687 * t373) * t369 + ((t685 * t363 - t646 * t364) * t373 + (t684 * t363 + t645 * t364) * t371) * qJD(3)) * t369 / 0.2e1 + ((t369 * t666 + t675) * t364 + (t369 * t665 + t674) * t363) * t623 + ((-t649 * t528 - t644) * t363 + ((t630 * t364 + (t647 + t650) * t363) * qJD(3) + t648) * t364) * t495 + ((-t650 * t527 + t644) * t364 + ((t647 * t363 + (t630 + t649) * t364) * qJD(3) + t648) * t363) * t492 + (t673 * t369 + (t670 * t577 + (t663 * t363 + t671 * t369 - t691) * t363) * t662) * t363 / 0.2e1 - (t672 * t369 + ((t369 * t722 + t691) * t364 + (-t663 * t364 + t369 * t723) * t363) * t662) * t364 / 0.2e1 + (t677 + t688) * t580 / 0.2e1 + (t676 + t689) * t577 / 0.2e1; -m(5) * (t217 * t74 + t218 * t73 + t224 * t72) - m(6) * (t217 * t65 + t218 * t64 + t224 * t57) + 0.2e1 * ((t527 * t73 + t528 * t74 - t8) * t628 + (t527 * t64 + t528 * t65 - t7) * t627) * t373 + 0.2e1 * ((qJD(3) * t72 + t36 * t363 + t364 * t37 + t577 * t74 - t580 * t73) * t628 + (qJD(3) * t57 + t21 * t363 + t22 * t364 + t577 * t65 - t580 * t64) * t627) * t371; m(6) * (t21 * t364 - t22 * t363);];
tauc = t1(:);
