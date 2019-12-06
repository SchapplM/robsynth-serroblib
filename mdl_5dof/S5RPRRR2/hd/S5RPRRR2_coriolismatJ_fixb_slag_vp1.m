% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR2
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:18
% EndTime: 2019-12-05 18:11:42
% DurationCPUTime: 10.14s
% Computational Cost: add. (53906->528), mult. (37229->690), div. (0->0), fcn. (34518->9), ass. (0->338)
t429 = sin(qJ(1));
t426 = t429 ^ 2;
t430 = cos(qJ(1));
t427 = t430 ^ 2;
t499 = t426 + t427;
t425 = pkin(9) + qJ(3);
t410 = qJ(4) + t425;
t403 = sin(t410);
t600 = pkin(4) * t403;
t407 = qJ(5) + t410;
t398 = cos(t407);
t556 = t398 * t429;
t397 = sin(t407);
t560 = t397 * t429;
t329 = rSges(6,1) * t560 + rSges(6,2) * t556;
t355 = rSges(6,1) * t397 + rSges(6,2) * t398;
t330 = t355 * t430;
t669 = t429 * t329 + t430 * t330;
t172 = -t499 * t600 - t669;
t408 = sin(t425);
t601 = pkin(3) * t408;
t365 = -t600 - t601;
t447 = t355 - t365;
t225 = t447 * t429;
t227 = t447 * t430;
t592 = rSges(6,1) * t398;
t356 = -rSges(6,2) * t397 + t592;
t404 = cos(t410);
t599 = pkin(4) * t404;
t487 = -t356 - t599;
t251 = t487 * t429;
t253 = t487 * t430;
t525 = -t225 * t251 - t227 * t253;
t295 = rSges(6,1) * t556 - rSges(6,2) * t560 - t430 * rSges(6,3);
t559 = t397 * t430;
t480 = -rSges(6,2) * t559 + t429 * rSges(6,3);
t555 = t398 * t430;
t191 = t429 * t295 + t430 * (rSges(6,1) * t555 + t480);
t409 = cos(t425);
t402 = pkin(3) * t409;
t405 = cos(pkin(9)) * pkin(2) + pkin(1);
t373 = t402 + t405;
t341 = t373 + t599;
t598 = -pkin(6) - qJ(2);
t424 = -pkin(7) + t598;
t492 = -pkin(8) + t424;
t395 = t429 * t492;
t400 = t429 * t424;
t503 = -t429 * t341 - t430 * t492;
t530 = t430 * t424;
t115 = -t429 * (t429 * t373 + t503 + t530) + t191 + t430 * (-t395 + t400 + (t341 - t373) * t430);
t406 = t429 * t598;
t500 = t405 - t373;
t522 = -t429 * ((-t424 + t598) * t430 + t500 * t429) + t430 * (-t500 * t430 - t400 + t406);
t86 = t115 + t522;
t57 = t86 * t172 + t525;
t551 = t403 * t430;
t481 = -rSges(5,2) * t551 + t429 * rSges(5,3);
t552 = t403 * t429;
t502 = rSges(5,2) * t552 + t430 * rSges(5,3);
t547 = t404 * t430;
t548 = t404 * t429;
t199 = t429 * (rSges(5,1) * t548 - t502) + t430 * (rSges(5,1) * t547 + t481);
t123 = t199 + t522;
t593 = rSges(5,1) * t404;
t364 = -rSges(5,2) * t403 + t593;
t363 = rSges(5,1) * t403 + rSges(5,2) * t404;
t440 = t363 + t601;
t674 = t440 * t430;
t675 = t440 * t429;
t666 = t429 * t675 + t430 * t674;
t339 = t363 * t429;
t340 = t363 * t430;
t668 = t429 * t339 + t430 * t340;
t79 = -t123 * t668 + t666 * t364;
t693 = -m(5) * t79 - m(6) * t57;
t664 = t225 * t429 + t227 * t430;
t371 = rSges(4,1) * t408 + rSges(4,2) * t409;
t673 = t499 * t371;
t684 = -m(6) / 0.2e1;
t685 = -m(5) / 0.2e1;
t493 = t664 * t684 + t666 * t685 - m(4) * t673 / 0.2e1;
t348 = t371 * t429;
t349 = t371 * t430;
t217 = t429 * t348 + t349 * t430;
t229 = -t365 * t429 + t329;
t357 = t430 * t365;
t230 = t357 - t330;
t651 = m(6) / 0.2e1;
t652 = m(5) / 0.2e1;
t686 = m(4) / 0.2e1;
t494 = (t429 * t229 - t230 * t430) * t651 + t666 * t652 + t217 * t686;
t49 = t494 - t493;
t692 = t49 * qJD(1);
t488 = t355 + t600;
t250 = t488 * t429;
t461 = t488 * t430;
t570 = t461 * t430;
t665 = t250 * t429 + t570;
t689 = t363 * t499;
t526 = t665 * t684 + t685 * t689;
t248 = pkin(4) * t552 + t329;
t527 = (t429 * t248 + t570) * t651 + t668 * t652;
t70 = t527 - t526;
t691 = t70 * qJD(1);
t690 = t355 * t499;
t393 = Icges(6,4) * t398;
t352 = -Icges(6,2) * t397 + t393;
t353 = Icges(6,1) * t397 + t393;
t688 = t352 + t353;
t394 = Icges(5,4) * t404;
t360 = -Icges(5,2) * t403 + t394;
t361 = Icges(5,1) * t403 + t394;
t687 = t360 + t361;
t628 = t429 / 0.2e1;
t627 = -t430 / 0.2e1;
t681 = t430 / 0.2e1;
t581 = Icges(4,4) * t408;
t367 = Icges(4,2) * t409 + t581;
t370 = Icges(4,1) * t409 - t581;
t680 = (t370 / 0.2e1 - t367 / 0.2e1) * t408;
t483 = t373 + t593;
t215 = -t429 * t483 + t502 - t530;
t216 = t430 * t483 - t400 + t481;
t671 = t215 * t430 + t216 * t429;
t443 = t671 * t364;
t197 = -t295 + t503;
t198 = -t395 + (t341 + t592) * t430 + t480;
t528 = t253 * t197 + t251 * t198;
t595 = (-t229 * t461 - t230 * t250 + t528) * t651 + (-t443 + (t429 * t674 - t430 * t675) * t363) * t652;
t596 = (t225 * t461 - t227 * t248 + t528) * t651 + (-t339 * t674 + t340 * t675 - t443) * t652;
t679 = t595 - t596;
t629 = -t429 / 0.2e1;
t678 = t628 + t629;
t672 = t197 * t430 + t198 * t429;
t594 = rSges(4,1) * t409;
t484 = t405 + t594;
t545 = t408 * t429;
t501 = rSges(4,2) * t545 + t430 * rSges(4,3);
t238 = -t484 * t429 - t430 * t598 + t501;
t544 = t408 * t430;
t482 = -rSges(4,2) * t544 + t429 * rSges(4,3);
t239 = t430 * t484 - t406 + t482;
t670 = t238 * t430 + t239 * t429;
t401 = Icges(4,4) * t409;
t368 = -Icges(4,2) * t408 + t401;
t369 = Icges(4,1) * t408 + t401;
t303 = Icges(5,6) * t429 + t360 * t430;
t580 = Icges(5,4) * t403;
t362 = Icges(5,1) * t404 - t580;
t305 = Icges(5,5) * t429 + t362 * t430;
t244 = t305 * t548;
t358 = Icges(5,5) * t404 - Icges(5,6) * t403;
t564 = t358 * t430;
t301 = Icges(5,3) * t429 + t564;
t472 = t430 * t301 - t244;
t145 = -t303 * t552 - t472;
t302 = Icges(5,4) * t548 - Icges(5,2) * t552 - Icges(5,6) * t430;
t300 = Icges(5,5) * t548 - Icges(5,6) * t552 - Icges(5,3) * t430;
t383 = Icges(5,4) * t552;
t304 = Icges(5,1) * t548 - Icges(5,5) * t430 - t383;
t519 = -t429 * t300 - t304 * t547;
t146 = -t302 * t551 - t519;
t518 = t429 * t301 + t305 * t547;
t147 = -t303 * t551 + t518;
t469 = t303 * t403 - t300;
t567 = t302 * t403;
t661 = (-t146 * t430 + t147 * t429) * t681 + ((t145 - t244 + (t301 + t567) * t430 + t519) * t430 + t518 * t429) * t627 + ((t429 * t469 + t145 + t146 + t472) * t429 + (t429 * (-t304 * t404 + t567) + t147 - t518 + (t300 + t469) * t430) * t430) * t628;
t359 = Icges(5,2) * t404 + t580;
t660 = t687 * t404 / 0.2e1 + (t362 / 0.2e1 - t359 / 0.2e1) * t403;
t579 = Icges(6,4) * t397;
t351 = Icges(6,2) * t398 + t579;
t354 = Icges(6,1) * t398 - t579;
t467 = t688 * t398 / 0.2e1 + (-t351 / 0.2e1 + t354 / 0.2e1) * t397;
t282 = Icges(6,6) * t429 + t352 * t430;
t284 = Icges(6,5) * t429 + t354 * t430;
t235 = t284 * t556;
t350 = Icges(6,5) * t398 - Icges(6,6) * t397;
t566 = t350 * t430;
t280 = Icges(6,3) * t429 + t566;
t473 = t430 * t280 - t235;
t137 = -t282 * t560 - t473;
t281 = Icges(6,4) * t556 - Icges(6,2) * t560 - Icges(6,6) * t430;
t279 = Icges(6,5) * t556 - Icges(6,6) * t560 - Icges(6,3) * t430;
t376 = Icges(6,4) * t560;
t283 = Icges(6,1) * t556 - Icges(6,5) * t430 - t376;
t521 = -t429 * t279 - t283 * t555;
t138 = -t281 * t559 - t521;
t520 = t429 * t280 + t284 * t555;
t139 = -t282 * t559 + t520;
t470 = t282 * t397 - t279;
t568 = t281 * t397;
t486 = ((t137 - t235 + (t280 + t568) * t430 + t521) * t430 + t520 * t429) * t627 + (-t138 * t430 + t139 * t429) * t681 + ((t429 * t470 + t137 + t138 + t473) * t429 + (t139 - t520 + t429 * (-t283 * t398 + t568) + (t470 + t279) * t430) * t430) * t628;
t328 = Icges(4,5) * t429 + t370 * t430;
t504 = -t367 * t430 + t328;
t390 = Icges(4,4) * t545;
t543 = t409 * t429;
t327 = Icges(4,1) * t543 - Icges(4,5) * t430 - t390;
t505 = -Icges(4,2) * t543 + t327 - t390;
t326 = Icges(4,6) * t429 + t368 * t430;
t506 = -t369 * t430 - t326;
t325 = Icges(4,4) * t543 - Icges(4,2) * t545 - Icges(4,6) * t430;
t507 = t369 * t429 + t325;
t659 = (-t504 * t429 + t430 * t505) * t408 + (t506 * t429 + t430 * t507) * t409;
t508 = -t359 * t430 + t305;
t509 = -Icges(5,2) * t548 + t304 - t383;
t510 = -t361 * t430 - t303;
t511 = t361 * t429 + t302;
t658 = (-t508 * t429 + t509 * t430) * t403 + (t510 * t429 + t511 * t430) * t404;
t512 = -t351 * t430 + t284;
t513 = -Icges(6,2) * t556 + t283 - t376;
t514 = -t353 * t430 - t282;
t515 = t353 * t429 + t281;
t657 = (-t512 * t429 + t513 * t430) * t397 + (t514 * t429 + t515 * t430) * t398;
t656 = 0.4e1 * qJD(1);
t655 = 2 * qJD(3);
t653 = 2 * qJD(4);
t75 = t86 * t669;
t93 = t115 * t669;
t644 = m(6) * (-t75 - t93 + ((t227 + t461) * t430 + (t225 + t250) * t429) * t356);
t126 = t191 * t172;
t442 = (-t251 * t429 - t253 * t430) * t355;
t58 = t664 * t356 - t75;
t642 = m(6) * (t126 + t442 + t58);
t135 = t426 * (t365 + t601) + t430 * (pkin(3) * t544 + t357) - t669;
t64 = t665 * t356 - t93;
t432 = t442 + t64;
t640 = m(6) * (t191 * t135 + t432);
t524 = -t250 * t251 - t253 * t461;
t637 = m(6) * (t115 * t135 + t524);
t444 = t672 * t356;
t634 = m(6) * (t225 * t330 - t227 * t329 - t444);
t633 = m(6) * (-t444 + (-t229 * t430 - t230 * t429) * t355);
t632 = m(6) * (t250 * t330 - t329 * t461 - t444);
t631 = m(6) * (-t444 + (-t248 * t430 + t429 * t461) * t355);
t626 = m(3) * t499 * (rSges(3,3) + qJ(2));
t625 = m(4) * (t238 * t348 - t239 * t349);
t624 = m(4) * t670;
t114 = -t199 * t668 + t364 * t689;
t113 = m(5) * t114;
t621 = m(5) * (t215 * t675 - t216 * t674);
t620 = m(5) * (t215 * t339 - t216 * t340);
t619 = m(5) * t671;
t613 = m(6) * (t197 * t229 + t198 * t230);
t108 = -t191 * t669 + t356 * t690;
t612 = m(6) * t108;
t611 = m(6) * (t197 * t248 - t198 * t461);
t610 = m(6) * (t197 * t329 - t198 * t330);
t609 = m(6) * t672;
t602 = m(6) * (-t251 * t430 + t253 * t429);
t454 = Icges(6,5) * t397 + Icges(6,6) * t398;
t317 = t454 * t429;
t318 = t430 * t454;
t597 = (-t426 * t318 + (t429 * t317 + t657) * t430) * t628 + (-t427 * t317 + (t430 * t318 + t657) * t429) * t627;
t546 = t408 * t325;
t542 = t409 * t430;
t323 = Icges(4,5) * t543 - Icges(4,6) * t545 - Icges(4,3) * t430;
t517 = -t429 * t323 - t327 * t542;
t457 = Icges(4,5) * t409 - Icges(4,6) * t408;
t324 = Icges(4,3) * t429 + t430 * t457;
t516 = t429 * t324 + t328 * t542;
t445 = t669 * t651;
t466 = m(6) * t690;
t142 = t445 + t466 / 0.2e1;
t498 = t142 * qJD(1);
t497 = t640 / 0.2e1 + t597;
t496 = -t612 + t597;
t495 = t115 * t172 + t524;
t489 = -t364 - t402;
t455 = Icges(5,5) * t403 + Icges(5,6) * t404;
t333 = t455 * t429;
t334 = t430 * t455;
t485 = (-t426 * t334 + (t429 * t333 + t658) * t430) * t628 + (-t427 * t333 + (t430 * t334 + t658) * t429) * t627 + t597;
t255 = t328 * t543;
t471 = t430 * t324 - t255;
t468 = t408 * t326 - t323;
t465 = t499 * t601;
t462 = t113 + t485;
t456 = -Icges(4,5) * t408 - Icges(4,6) * t409;
t446 = t487 - t402;
t441 = t486 + t661;
t436 = (-t351 + t354) * t398 - t688 * t397;
t439 = -t486 + (t429 * t350 + t514 * t397 + t512 * t398 + t430 * t436) * t628 + (-t515 * t397 + t513 * t398 + t429 * t436 - t566) * t627;
t438 = -t467 + t678 * (t281 * t398 + t283 * t397);
t437 = t467 + t660;
t435 = (-t359 + t362) * t404 - t687 * t403;
t433 = t439 - t661 + (t429 * t358 + t510 * t403 + t508 * t404 + t430 * t435) * t628 + (-t511 * t403 + t509 * t404 + t429 * t435 - t564) * t627;
t431 = t438 - t660 + t678 * (t302 * t404 + t304 * t403);
t372 = -rSges(4,2) * t408 + t594;
t343 = t456 * t430;
t342 = t456 * t429;
t274 = t489 * t430;
t272 = t489 * t429;
t228 = t446 * t430;
t226 = t446 * t429;
t185 = -t465 - t668;
t165 = qJD(4) * t602;
t156 = -t326 * t544 + t516;
t155 = -t325 * t544 - t517;
t154 = -t326 * t545 - t471;
t141 = t445 - t466 / 0.2e1;
t127 = -t465 + t135;
t111 = -t155 * t430 + t156 * t429;
t110 = -(-t429 * (-t327 * t409 + t546) - t430 * t323) * t430 + t154 * t429;
t83 = t631 / 0.2e1;
t82 = t632 / 0.2e1;
t81 = t467 + t610;
t78 = t633 / 0.2e1;
t76 = t634 / 0.2e1;
t71 = t526 + t527;
t61 = t609 + t619 + t624 + t626;
t53 = t437 + t611 + t620;
t52 = (t154 - t255 + (t324 + t546) * t430 + t517) * t430 + t516 * t429;
t51 = (t430 * t468 + t156 - t516) * t430 + (t429 * t468 + t155 + t471) * t429;
t48 = t493 + t494;
t29 = t642 / 0.2e1;
t27 = (t369 / 0.2e1 + t368 / 0.2e1) * t409 + t680 + t625 + t621 + t613 + t437;
t24 = t644 / 0.2e1;
t23 = t597 + t612;
t22 = t23 * qJD(5);
t21 = m(6) * t64 + t597;
t20 = m(6) * t58 + t597;
t16 = t462 + t637;
t15 = t485 - t693;
t14 = t24 - t642 / 0.2e1 + t497;
t13 = t24 + t29 - t640 / 0.2e1 + t597;
t12 = t29 - t644 / 0.2e1 + t497;
t11 = t82 - t631 / 0.2e1 + t486;
t10 = t83 - t632 / 0.2e1 + t486;
t9 = t76 - t633 / 0.2e1 + t486;
t8 = t78 - t634 / 0.2e1 + t486;
t7 = t82 + t83 + t439;
t6 = t76 + t78 + t439;
t4 = (-t52 / 0.2e1 + t111 / 0.2e1) * t430 + (t110 / 0.2e1 + t51 / 0.2e1) * t429 + t441;
t3 = t441 + t679;
t2 = t441 - t679;
t1 = t433 + t595 + t596;
t5 = [qJD(2) * t61 + qJD(3) * t27 + qJD(4) * t53 + qJD(5) * t81, qJD(1) * t61 + qJD(3) * t48 + qJD(4) * t71 + qJD(5) * t141, t27 * qJD(1) + t48 * qJD(2) + t1 * qJD(4) + t6 * qJD(5) + ((-t670 * t372 + (-t348 * t430 + t349 * t429) * t371) * t686 + (t215 * t274 + t216 * t272) * t652 + (t197 * t228 + t198 * t226 - t225 * t230 - t227 * t229) * t651) * t655 + (t433 + (t408 * t506 + t409 * t504) * t628 + t52 * t681 + (t110 + t51) * t629 + (-t408 * t507 + t409 * t505 + t111) * t627 + (t427 / 0.2e1 + t426 / 0.2e1) * t457) * qJD(3), t53 * qJD(1) + t71 * qJD(2) + t1 * qJD(3) + t433 * qJD(4) + t7 * qJD(5) + ((t528 + (-t248 + t250) * t461) * t651 + (-t443 + (-t339 * t430 + t340 * t429) * t363) * t652) * t653, t81 * qJD(1) + t141 * qJD(2) + t6 * qJD(3) + t7 * qJD(4) + ((-t444 + (-t329 * t430 + t330 * t429) * t355) * m(6) + t439) * qJD(5); t49 * qJD(3) + t70 * qJD(4) + t142 * qJD(5) + (-t609 / 0.4e1 - t619 / 0.4e1 - t624 / 0.4e1 - t626 / 0.4e1) * t656, 0, t692 + ((-t272 * t430 + t274 * t429) * t652 + (-t226 * t430 + t228 * t429) * t651) * t655 + t165, qJD(3) * t602 + t165 + t691, t498; (t431 - (t369 + t368) * t409 / 0.2e1 - t680) * qJD(1) - t49 * qJD(2) + t4 * qJD(3) + t2 * qJD(4) + t9 * qJD(5) + (-t625 / 0.4e1 - t613 / 0.4e1 - t621 / 0.4e1) * t656, -t692, t4 * qJD(1) + (m(6) * (t127 * t86 - t225 * t226 - t227 * t228) + m(5) * (t123 * t185 - t272 * t675 - t274 * t674) + (t427 * t342 + (-t430 * t343 + t659) * t429) * t627 + (t426 * t343 + (-t429 * t342 + t659) * t430) * t628 + m(4) * (t372 * t673 - (t429 * (rSges(4,1) * t543 - t501) + t430 * (rSges(4,1) * t542 + t482)) * t217) + t485) * qJD(3) + t15 * qJD(4) + t20 * qJD(5), t2 * qJD(1) + t15 * qJD(3) + t13 * qJD(5) + ((t495 + t57) * t651 + (t79 + t114) * t652) * t653 + (t485 - t113 - t637) * qJD(4), t9 * qJD(1) + t20 * qJD(3) + t13 * qJD(4) + (m(6) * (t58 + t108) + t496) * qJD(5); t431 * qJD(1) - t70 * qJD(2) + t3 * qJD(3) + t441 * qJD(4) + t11 * qJD(5) + (-t611 / 0.4e1 - t620 / 0.4e1) * t656, -t691, t3 * qJD(1) + t16 * qJD(4) + t14 * qJD(5) + ((t115 * t127 + t135 * t86 - t226 * t250 - t228 * t461 + t525) * t651 + (t199 * t185 + (-t272 * t429 - t274 * t430) * t363 + t79) * t652) * t655 + (t485 + t693) * qJD(3), t441 * qJD(1) + t16 * qJD(3) + (m(6) * t495 + t462) * qJD(4) + t21 * qJD(5), t11 * qJD(1) + t14 * qJD(3) + t21 * qJD(4) + (m(6) * (t64 + t108) + t496) * qJD(5); (t438 - t610) * qJD(1) - t142 * qJD(2) + t8 * qJD(3) + t10 * qJD(4) + t486 * qJD(5), -t498, t8 * qJD(1) + ((t191 * t127 + (-t226 * t429 - t228 * t430) * t355) * m(6) + t597) * qJD(3) + t12 * qJD(4) + t22, t10 * qJD(1) + t12 * qJD(3) + ((t126 - t64 + t432) * m(6) + t597) * qJD(4) + t22, qJD(1) * t486 + t22 + (qJD(3) + qJD(4)) * t23;];
Cq = t5;
