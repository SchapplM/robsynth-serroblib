% Calculate matrix of centrifugal and coriolis load on the joints for
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
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:57
% EndTime: 2022-01-23 09:25:25
% DurationCPUTime: 17.79s
% Computational Cost: add. (48789->613), mult. (52403->837), div. (0->0), fcn. (56339->10), ass. (0->367)
t737 = Icges(5,3) + Icges(4,3);
t485 = cos(qJ(1));
t478 = qJ(3) + pkin(9);
t470 = sin(t478);
t482 = sin(qJ(3));
t629 = pkin(3) * t482;
t447 = pkin(4) * t470 + t629;
t481 = qJ(4) + pkin(6);
t476 = -pkin(7) - t481;
t483 = sin(qJ(1));
t479 = sin(pkin(8));
t600 = t479 * t485;
t504 = (qJ(2) + t447) * t483 - t476 * t600;
t471 = cos(t478);
t484 = cos(qJ(3));
t475 = t484 * pkin(3);
t448 = pkin(4) * t471 + t475;
t446 = pkin(2) + t448;
t480 = cos(pkin(8));
t516 = t446 * t480 + pkin(1);
t628 = pkin(2) * t480 + pkin(1);
t421 = t480 * t475 + t479 * t481 + t628;
t468 = qJ(2) + t629;
t707 = t421 * t485 + t468 * t483;
t222 = t516 * t485 + t504 - t707;
t472 = qJ(5) + t478;
t466 = sin(t472);
t587 = t485 * t466;
t467 = cos(t472);
t593 = t483 * t467;
t396 = -t480 * t587 + t593;
t594 = t483 * t466;
t596 = t480 * t485;
t397 = t467 * t596 + t594;
t506 = t397 * rSges(6,1) + t396 * rSges(6,2);
t289 = rSges(6,3) * t600 + t506;
t449 = pkin(6) * t479 + t628;
t595 = t483 * qJ(2);
t685 = -t449 * t485 - t595 + t707;
t292 = t480 * t685;
t344 = (t476 + t481) * t480 + (-pkin(2) + t446 - t475) * t479;
t372 = -rSges(6,3) * t480 + (rSges(6,1) * t467 - rSges(6,2) * t466) * t479;
t422 = t479 * t475 + (pkin(6) - t481) * t480;
t378 = t422 * t600;
t681 = (t344 + t372) * t600 + t378;
t100 = t292 + (t222 + t289) * t480 + t681;
t580 = -t685 - t222;
t533 = -t289 + t580;
t101 = t533 * t480 - t681;
t744 = t100 + t101;
t500 = rSges(6,3) * t479 + t516;
t202 = t485 * t500 + t504 + t506;
t394 = t467 * t485 + t480 * t594;
t395 = t480 * t593 - t587;
t321 = -t394 * rSges(6,1) - t395 * rSges(6,2);
t322 = t396 * rSges(6,1) - t397 * rSges(6,2);
t473 = t485 * qJ(2);
t444 = t485 * t447;
t601 = t479 * t483;
t687 = t476 * t601 + t444;
t530 = t473 + t687;
t689 = -t395 * rSges(6,1) + t394 * rSges(6,2);
t706 = -t483 * t500 + t530 + t689;
t743 = m(6) * (t202 * t322 - t321 * t706);
t616 = Icges(5,4) * t471;
t375 = -Icges(5,6) * t480 + (-Icges(5,2) * t470 + t616) * t479;
t617 = Icges(5,4) * t470;
t376 = -Icges(5,5) * t480 + (Icges(5,1) * t471 - t617) * t479;
t619 = Icges(4,4) * t484;
t410 = -Icges(4,6) * t480 + (-Icges(4,2) * t482 + t619) * t479;
t620 = Icges(4,4) * t482;
t411 = -Icges(4,5) * t480 + (Icges(4,1) * t484 - t620) * t479;
t413 = (-Icges(5,2) * t471 - t617) * t479;
t414 = (-Icges(5,1) * t470 - t616) * t479;
t436 = (-Icges(4,2) * t484 - t620) * t479;
t437 = (-Icges(4,1) * t482 - t619) * t479;
t742 = ((t376 / 0.2e1 + t413 / 0.2e1) * t470 - (t414 / 0.2e1 - t375 / 0.2e1) * t471 + (t411 / 0.2e1 + t436 / 0.2e1) * t482 - (t437 / 0.2e1 - t410 / 0.2e1) * t484) * t479;
t592 = t483 * t470;
t416 = t471 * t485 + t480 * t592;
t586 = t485 * t470;
t591 = t483 * t471;
t417 = t480 * t591 - t586;
t588 = t484 * t485;
t590 = t483 * t482;
t438 = t480 * t590 + t588;
t585 = t485 * t482;
t589 = t483 * t484;
t439 = t480 * t589 - t585;
t741 = Icges(4,5) * t439 + Icges(5,5) * t417 - Icges(4,6) * t438 - Icges(5,6) * t416 + t737 * t601;
t740 = t479 ^ 2;
t287 = rSges(6,3) * t601 - t689;
t686 = -t421 * t483 + t468 * t485;
t510 = -t473 + t686;
t303 = t449 * t483 + t510;
t739 = t303 - t287;
t266 = t480 * t287;
t361 = t372 * t601;
t738 = t361 + t266;
t307 = t417 * rSges(5,1) - t416 * rSges(5,2) + rSges(5,3) * t601;
t377 = t422 * t601;
t379 = -rSges(5,3) * t480 + (rSges(5,1) * t471 - rSges(5,2) * t470) * t479;
t561 = t379 * t601 + t377;
t164 = (-t303 + t307) * t480 + t561;
t291 = t480 * t303;
t165 = t307 * t480 - t291 + t561;
t736 = t164 - t165;
t735 = t741 * t600;
t400 = Icges(5,4) * t417;
t297 = -Icges(5,2) * t416 + Icges(5,6) * t601 + t400;
t399 = Icges(5,4) * t416;
t301 = -Icges(5,1) * t417 - Icges(5,5) * t601 + t399;
t429 = Icges(4,4) * t439;
t335 = -Icges(4,2) * t438 + Icges(4,6) * t601 + t429;
t428 = Icges(4,4) * t438;
t339 = -Icges(4,1) * t439 - Icges(4,5) * t601 + t428;
t418 = -t480 * t586 + t591;
t419 = t471 * t596 + t592;
t441 = t480 * t588 + t590;
t704 = t480 * t585 - t589;
t734 = -t418 * t297 + t419 * t301 + t335 * t704 + t441 * t339;
t688 = -t439 * rSges(4,1) + t438 * rSges(4,2);
t341 = rSges(4,3) * t601 - t688;
t415 = -rSges(4,3) * t480 + (rSges(4,1) * t484 - rSges(4,2) * t482) * t479;
t733 = t341 * t480 + t415 * t601;
t732 = -t479 / 0.2e1;
t231 = -t480 * t444 + t483 * t448 + t322;
t730 = t231 * t483;
t614 = Icges(6,4) * t466;
t371 = -Icges(6,5) * t480 + (Icges(6,1) * t467 - t614) * t479;
t559 = t371 + (-Icges(6,2) * t467 - t614) * t479;
t729 = t559 * t466;
t728 = -t734 + t735;
t727 = -t307 + t686;
t381 = Icges(6,4) * t395;
t281 = -Icges(6,2) * t394 + Icges(6,6) * t601 + t381;
t380 = Icges(6,4) * t394;
t285 = -Icges(6,1) * t395 - Icges(6,5) * t601 + t380;
t725 = t396 * t281 - t397 * t285;
t673 = m(4) / 0.2e1;
t672 = m(5) / 0.2e1;
t670 = m(6) / 0.2e1;
t698 = t479 / 0.2e1;
t278 = Icges(6,5) * t395 - Icges(6,6) * t394 + Icges(6,3) * t601;
t716 = t278 * t600;
t119 = t716 + t725;
t693 = t119 - t716;
t719 = (t693 - t725) * t600;
t567 = -t418 * rSges(5,1) + t419 * rSges(5,2) + t704 * pkin(3);
t711 = t567 * t483;
t618 = Icges(5,4) * t419;
t299 = Icges(5,2) * t418 + Icges(5,6) * t600 + t618;
t401 = Icges(5,4) * t418;
t302 = Icges(5,1) * t419 + Icges(5,5) * t600 + t401;
t621 = Icges(4,4) * t441;
t337 = -Icges(4,2) * t704 + Icges(4,6) * t600 + t621;
t430 = Icges(4,4) * t704;
t340 = Icges(4,1) * t441 + Icges(4,5) * t600 - t430;
t710 = t418 * t299 + t419 * t302 - t704 * t337 + t441 * t340 + (Icges(4,5) * t441 + Icges(5,5) * t419 - Icges(4,6) * t704 + Icges(5,6) * t418 + t737 * t600) * t600;
t515 = rSges(4,3) * t479 + t449;
t705 = -t515 * t483 + t473 + t688;
t508 = t441 * rSges(4,1) - rSges(4,2) * t704;
t343 = rSges(4,3) * t600 + t508;
t243 = t480 * t343 + t415 * t600;
t694 = t516 * t483;
t220 = t510 - t687 + t694;
t531 = t344 * t601 + t361 + t377;
t99 = t220 * t480 + t266 - t291 + t531;
t57 = t101 * t483 + t99 * t485;
t511 = -t379 * t600 - t378;
t309 = t419 * rSges(5,1) + t418 * rSges(5,2) + rSges(5,3) * t600;
t568 = -t685 - t309;
t167 = t568 * t480 + t511;
t90 = t165 * t485 + t167 * t483;
t545 = (-t243 * t483 + t485 * t733) * t673 + t57 * t670 + t90 * t672;
t166 = t480 * t309 + t292 - t511;
t534 = -t530 + t686 + t694 - t739;
t98 = t534 * t480 + t531;
t548 = (t100 * t483 - t485 * t98 + t57) * t670 + (-t164 * t485 + t166 * t483 + t90) * t672;
t703 = t545 - t548;
t155 = t167 * t600;
t91 = t101 * t600;
t626 = (-t99 * t601 + t91) * t670 + (-t165 * t601 + t155) * t672;
t47 = t685 * t600 + ((t222 + t580) * t485 + (-t220 + t534 + t739) * t483) * t479;
t625 = t98 - t99;
t627 = (-t480 * t47 + t91 + (t100 * t485 + t625 * t483) * t479) * t670 + (t155 + (t166 * t485 + t736 * t483) * t479) * t672;
t702 = t626 - t627;
t701 = t728 - t735;
t391 = (-Icges(6,5) * t466 - Icges(6,6) * t467) * t479;
t599 = t480 * t391;
t700 = -t599 / 0.2e1 + t729 * t732;
t649 = -t480 / 0.2e1;
t696 = m(6) * t479;
t290 = t321 * t600;
t195 = -t322 * t601 + t290;
t398 = (-rSges(6,1) * t466 - rSges(6,2) * t467) * t479;
t228 = t480 * t321 + t398 * t601;
t229 = -t480 * t322 - t398 * t600;
t312 = -Icges(6,5) * t394 - Icges(6,6) * t395;
t575 = -Icges(6,2) * t395 - t285 - t380;
t577 = -Icges(6,1) * t394 - t281 - t381;
t103 = -t312 * t480 + (-t575 * t466 + t577 * t467) * t479;
t313 = Icges(6,5) * t396 - Icges(6,6) * t397;
t382 = Icges(6,4) * t396;
t286 = Icges(6,1) * t397 + Icges(6,5) * t600 + t382;
t574 = -Icges(6,2) * t397 + t286 + t382;
t615 = Icges(6,4) * t397;
t283 = Icges(6,2) * t396 + Icges(6,6) * t600 + t615;
t576 = Icges(6,1) * t396 - t283 - t615;
t104 = -t313 * t480 + (-t574 * t466 + t576 * t467) * t479;
t613 = Icges(6,4) * t467;
t370 = -Icges(6,6) * t480 + (-Icges(6,2) * t466 + t613) * t479;
t393 = (-Icges(6,1) * t466 - t613) * t479;
t560 = -t370 + t393;
t122 = t391 * t601 - t559 * t394 + t560 * t395;
t123 = t391 * t600 + t559 * t396 + t560 * t397;
t522 = t600 / 0.2e1;
t524 = t601 / 0.2e1;
t168 = -t599 + (t560 * t467 - t729) * t479;
t609 = t168 * t480;
t547 = ((t313 * t601 - t574 * t394 + t576 * t395) * t600 + (t312 * t601 - t575 * t394 + t577 * t395) * t601 - t122 * t480) * t524 + ((t313 * t600 + t574 * t396 + t576 * t397) * t600 + (t312 * t600 + t575 * t396 + t577 * t397) * t601 - t123 * t480) * t522 + (-t609 + (t103 * t483 + t104 * t485) * t479) * t649;
t258 = t287 * t600;
t276 = t303 * t600;
t75 = t258 - t276 + (t220 * t485 + t533 * t483) * t479;
t9 = t547 + m(6) * (t101 * t229 + t75 * t195 + t99 * t228);
t695 = t9 * qJD(5);
t554 = t411 + t436;
t555 = -t410 + t437;
t557 = t376 + t413;
t558 = -t375 + t414;
t435 = (-Icges(4,5) * t482 - Icges(4,6) * t484) * t479;
t597 = t480 * t435;
t412 = (-Icges(5,5) * t470 - Icges(5,6) * t471) * t479;
t598 = t480 * t412;
t690 = t597 + t598 + (t557 * t470 - t558 * t471 + t554 * t482 - t555 * t484) * t479;
t225 = t309 + t707;
t493 = -t483 * t480 * t447 - t448 * t485;
t230 = -t321 - t493;
t330 = -rSges(5,1) * t416 - rSges(5,2) * t417;
t551 = t438 * pkin(3);
t264 = -t330 + t551;
t583 = (t230 * t485 + t730) * t696 / 0.2e1 + m(5) * (t264 * t485 - t711) * t698;
t355 = -rSges(4,1) * t438 - rSges(4,2) * t439;
t356 = -rSges(4,1) * t704 - rSges(4,2) * t441;
t535 = (t483 * t230 - t231 * t485) * t670 + (t483 * t264 + t485 * t567) * t672 + (-t483 * t355 - t356 * t485) * t673;
t324 = -Icges(5,5) * t416 - Icges(5,6) * t417;
t325 = Icges(5,5) * t418 - Icges(5,6) * t419;
t349 = -Icges(4,5) * t438 - Icges(4,6) * t439;
t350 = -Icges(4,5) * t704 - Icges(4,6) * t441;
t680 = (t325 + t350) * t600 + (t349 + t324) * t601;
t445 = (-t483 ^ 2 - t485 ^ 2) * t479;
t678 = 0.2e1 * t445;
t677 = 2 * qJD(1);
t676 = 4 * qJD(1);
t675 = 2 * qJD(3);
t674 = 4 * qJD(3);
t671 = m(5) / 0.4e1;
t669 = m(6) / 0.4e1;
t667 = m(5) * (t164 * t167 + t165 * t166);
t188 = -t289 * t601 + t258;
t215 = t480 * t289 + t372 * t600;
t662 = m(6) * (t188 * t47 - t625 * t215 + t744 * t738);
t661 = m(6) * (t100 * t99 + t101 * t98 + t47 * t75);
t584 = t229 * t202 + t228 * t706;
t657 = m(6) * (t101 * t322 - t321 * t99 + t584);
t656 = m(6) * (-t215 * t231 + t230 * t738 + t584);
t651 = m(6) * (t202 * t231 + t230 * t706);
t648 = m(3) * ((rSges(3,2) * t601 + rSges(3,3) * t485 + t473) * t485 + (-rSges(3,2) * t600 + (rSges(3,3) + qJ(2)) * t483) * t483);
t261 = t515 * t485 + t508 + t595;
t646 = m(4) * (t261 * t356 - t355 * t705);
t644 = m(4) * (t261 * t483 + t485 * t705);
t642 = m(5) * (-t225 * t567 + t264 * t727);
t216 = t225 * t600;
t641 = m(5) * (-t601 * t727 + t216);
t640 = m(5) * (t225 * t483 + t485 * t727);
t194 = t202 * t600;
t637 = m(6) * (-t601 * t706 + t194);
t636 = m(6) * (t202 * t483 + t485 * t706);
t635 = m(6) * (-t215 * t600 - t601 * t738);
t634 = m(6) * (-t215 * t483 + t485 * t738);
t631 = (-t321 * t485 + t322 * t483) * t696;
t630 = m(6) * (-t483 * t321 - t322 * t485);
t624 = m(6) * qJD(3);
t623 = m(6) * qJD(5);
t608 = ((-Icges(6,3) * t480 + (Icges(6,5) * t467 - Icges(6,6) * t466) * t479) * t600 + t396 * t370 + t397 * t371) * t480;
t603 = t467 * t479;
t573 = -Icges(5,1) * t416 - t297 - t400;
t572 = Icges(5,1) * t418 - t299 - t618;
t571 = -Icges(5,2) * t417 - t301 - t399;
t570 = -Icges(5,2) * t419 + t302 + t401;
t566 = -Icges(4,1) * t438 - t335 - t429;
t565 = -Icges(4,1) * t704 - t337 - t621;
t564 = -Icges(4,2) * t439 - t339 - t428;
t563 = -Icges(4,2) * t441 + t340 - t430;
t251 = (t670 + t672) * t678;
t549 = t251 * qJD(1);
t546 = t740 * t629;
t120 = (Icges(6,5) * t397 + Icges(6,6) * t396 + Icges(6,3) * t600) * t600 + t396 * t283 + t397 * t286;
t527 = -t603 / 0.2e1;
t526 = t603 / 0.2e1;
t525 = -t601 / 0.2e1;
t25 = -t608 + ((t278 * t601 + t120) * t485 + t693 * t483) * t479;
t64 = -t608 + (t119 * t483 + t120 * t485) * t479;
t512 = t25 * t524 + t719 * t522 + t64 * t525;
t505 = t370 * t527 + t393 * t526 + t700;
t503 = t662 / 0.2e1 + t512;
t502 = (t701 + t734) * t698 * t485 + (t480 / 0.2e1 + t649) * (-t375 * t416 + t376 * t417 - t410 * t438 + t411 * t439 + (-t737 * t480 + (Icges(4,5) * t484 + Icges(5,5) * t471 - Icges(4,6) * t482 - Icges(5,6) * t470) * t479) * t601);
t501 = (t728 * t483 + t710 * t485) * t732 + ((t741 * t601 + t710) * t485 + t701 * t483) * t698;
t490 = t25 * t525 - t719 * t600 / 0.2e1 + (t123 + t104) * t522 + (t64 + t122 + t103) * t524;
t489 = t370 * t526 + t393 * t527 - t700;
t488 = t490 - t609;
t453 = t485 * t546;
t442 = (-rSges(4,1) * t482 - rSges(4,2) * t484) * t479;
t420 = (-rSges(5,1) * t470 - rSges(5,2) * t471) * t479;
t408 = t480 * t551;
t390 = (-t447 + t629) * t479;
t387 = t551 * t600;
t310 = t493 + t551;
t260 = -t480 * t356 - t442 * t600;
t259 = t355 * t480 + t442 * t601;
t250 = (t669 + t671) * t678 - (m(6) + m(5)) * t445 / 0.2e1;
t205 = -t420 * t600 + t567 * t480 + t453;
t204 = t330 * t480 - t408 + (t420 * t479 - t546) * t483;
t199 = t630 / 0.2e1;
t193 = t631 / 0.2e1;
t174 = -t387 + (t330 * t485 + t711) * t479;
t170 = t435 * t600 + t555 * t441 - t554 * t704;
t169 = t435 * t601 - t554 * t438 + t555 * t439;
t156 = t228 * t483 - t229 * t485;
t149 = t156 * t623;
t139 = t453 + (-t390 - t398) * t600 - t231 * t480;
t138 = t310 * t480 - t408 + (t390 * t479 - t546) * t483 + t228;
t136 = t634 / 0.2e1;
t135 = t412 * t600 + t557 * t418 + t558 * t419;
t134 = t412 * t601 - t557 * t416 + t558 * t417;
t126 = t635 / 0.2e1;
t125 = -t350 * t480 + (-t563 * t482 + t565 * t484) * t479;
t124 = -t349 * t480 + (-t564 * t482 + t566 * t484) * t479;
t111 = t290 - t387 + (t310 * t485 - t730) * t479;
t110 = -t325 * t480 + (-t570 * t470 + t572 * t471) * t479;
t109 = -t324 * t480 + (-t571 * t470 + t573 * t471) * t479;
t107 = -t480 * t195 + (t228 * t485 + t229 * t483) * t479;
t106 = t107 * t623;
t76 = t637 + t641;
t72 = t505 + t743;
t58 = t636 + t640 + t644 + t648;
t50 = t136 - t630 / 0.2e1;
t49 = t199 + t136;
t48 = t199 - t634 / 0.2e1;
t45 = t656 / 0.2e1;
t43 = t126 - t631 / 0.2e1;
t42 = t193 + t126;
t41 = t193 - t635 / 0.2e1;
t35 = t657 / 0.2e1;
t29 = (-t435 / 0.2e1 - t412 / 0.2e1) * t480 + t646 + t642 + t651 - t742 + t505;
t14 = m(6) * (t188 * t195 - t215 * t229 + t228 * t738) + t547;
t13 = t14 * qJD(5);
t12 = t545 + t548 - t535;
t11 = t535 - t703;
t10 = t535 + t703;
t8 = t583 + t702;
t7 = t626 + t627 - t583;
t6 = t583 - t702;
t4 = t45 - t657 / 0.2e1 + t503;
t3 = t35 - t656 / 0.2e1 + t503;
t2 = t35 + t45 - t662 / 0.2e1 + t488;
t1 = t667 + t661 + (t483 * t501 + t485 * t502) * t479 + t512;
t5 = [t58 * qJD(2) + t29 * qJD(3) + t76 * qJD(4) + t72 * qJD(5), qJD(1) * t58 + qJD(3) * t10 + qJD(4) * t250 + qJD(5) * t49, t29 * qJD(1) + t10 * qJD(2) + t8 * qJD(4) + t2 * qJD(5) + (-t667 / 0.4e1 - t661 / 0.4e1) * t674 + ((t101 * t231 + t138 * t706 + t139 * t202 + t230 * t99) * t670 + (t165 * t264 - t167 * t567 + t204 * t727 + t205 * t225) * t672 + (-t243 * t356 + t259 * t705 + t260 * t261 - t355 * t733) * t673) * t675 + (t490 + (-t168 + t690) * t480 + ((t125 / 0.2e1 + t110 / 0.2e1 + t170 / 0.2e1 + t135 / 0.2e1 - t502) * t485 + (t124 / 0.2e1 + t109 / 0.2e1 + t169 / 0.2e1 + t134 / 0.2e1 - t501) * t483) * t479) * qJD(3), qJD(1) * t76 + qJD(2) * t250 + qJD(3) * t8 + qJD(5) * t42, t72 * qJD(1) + t49 * qJD(2) + t2 * qJD(3) + t42 * qJD(4) + ((-t215 * t322 - t321 * t738 + t584) * m(6) + t488) * qJD(5); t11 * qJD(3) + t251 * qJD(4) + t48 * qJD(5) + (-t636 / 0.4e1 - t640 / 0.4e1 - t644 / 0.4e1 - t648 / 0.4e1) * t676, 0, t11 * qJD(1) + t149 + ((t259 * t483 - t260 * t485) * t673 + (t204 * t483 - t205 * t485) * t672 + (t138 * t483 - t139 * t485) * t670) * t675, t549, t48 * qJD(1) + t156 * t624 + t149; t12 * qJD(2) + t1 * qJD(3) + t7 * qJD(4) + t3 * qJD(5) + (-t646 / 0.4e1 - t642 / 0.4e1 - t651 / 0.4e1) * t676 + ((t202 * t625 + t744 * t706) * t670 + (t736 * t225 + (t166 + t167) * t727) * t672) * t677 + (t597 / 0.2e1 + t598 / 0.2e1 + t489 + t742) * qJD(1), t12 * qJD(1), t1 * qJD(1) + (t547 + (t690 * t480 + ((t110 + t125) * t485 + (t109 + t124) * t483) * t479) * t649 + ((-t416 * t570 + t417 * t572 - t438 * t563 + t439 * t565) * t600 + (-t169 - t134) * t480 + (-t416 * t571 + t417 * t573 - t438 * t564 + t439 * t566 + t680) * t601) * t524 + ((t418 * t571 + t419 * t573 + t441 * t566 - t564 * t704) * t601 + (-t170 - t135) * t480 + (t418 * t570 + t419 * t572 + t441 * t565 - t563 * t704 + t680) * t600) * t522) * qJD(3) + t695 + ((t101 * t139 + t111 * t75 + t138 * t99) * t669 + ((-t276 + (t307 * t485 + t568 * t483) * t479) * t174 + t165 * t204 + t167 * t205) * t671 + m(4) * (t733 * t259 - t243 * t260 + (t341 * t485 - t343 * t483) * t740 * (t355 * t485 - t356 * t483)) / 0.4e1) * t674, t7 * qJD(1), t3 * qJD(1) + t9 * qJD(3) + t695; -t251 * qJD(2) + t6 * qJD(3) + t41 * qJD(5) + (-t637 / 0.4e1 - t641 / 0.4e1) * t676 + (t194 * t670 + t216 * t672 + (-t202 * t670 - t225 * t672) * t600) * t677, -t549, t6 * qJD(1) + ((-t480 * t111 + (t138 * t485 + t139 * t483) * t479) * t670 + (-t480 * t174 + (t204 * t485 + t205 * t483) * t479) * t672) * t675 + t106, 0, t41 * qJD(1) + t107 * t624 + t106; (t489 - t743) * qJD(1) + t50 * qJD(2) + t4 * qJD(3) + t43 * qJD(4) + t512 * qJD(5), t50 * qJD(1), t4 * qJD(1) + ((t111 * t188 + t138 * t738 - t139 * t215) * m(6) + t547) * qJD(3) + t13, t43 * qJD(1), qJD(1) * t512 + qJD(3) * t14 + t13;];
Cq = t5;
