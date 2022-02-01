% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:41
% EndTime: 2022-01-20 12:01:59
% DurationCPUTime: 9.01s
% Computational Cost: add. (76859->480), mult. (43211->599), div. (0->0), fcn. (38666->10), ass. (0->299)
t444 = qJ(1) + qJ(2);
t442 = qJ(3) + t444;
t435 = cos(t442);
t434 = sin(t442);
t443 = qJ(4) + qJ(5);
t439 = cos(t443);
t543 = t434 * t439;
t437 = sin(t443);
t544 = t434 * t437;
t329 = rSges(6,1) * t543 - rSges(6,2) * t544 - t435 * rSges(6,3);
t447 = cos(qJ(4));
t585 = pkin(4) * t447;
t436 = pkin(3) + t585;
t655 = -pkin(9) - pkin(8);
t502 = -t434 * t436 - t435 * t655;
t288 = -t329 + t502;
t438 = sin(t444);
t588 = pkin(2) * t438;
t280 = t288 - t588;
t590 = sin(qJ(1)) * pkin(1);
t274 = t280 - t590;
t418 = t434 * t655;
t536 = t435 * t437;
t493 = -rSges(6,2) * t536 + t434 * rSges(6,3);
t582 = rSges(6,1) * t439;
t289 = -t418 + (t436 + t582) * t435 + t493;
t440 = cos(t444);
t587 = pkin(2) * t440;
t281 = t289 + t587;
t589 = cos(qJ(1)) * pkin(1);
t275 = t281 + t589;
t146 = -t289 * t274 + t275 * t288;
t429 = t435 * pkin(8);
t583 = rSges(5,1) * t447;
t495 = pkin(3) + t583;
t445 = sin(qJ(4));
t542 = t434 * t445;
t501 = rSges(5,2) * t542 + t435 * rSges(5,3);
t299 = -t434 * t495 + t429 + t501;
t290 = t299 - t588;
t286 = t290 - t590;
t534 = t435 * t445;
t408 = rSges(5,2) * t534;
t300 = -t408 + t495 * t435 + (rSges(5,3) + pkin(8)) * t434;
t291 = t300 + t587;
t287 = t291 + t589;
t155 = -t300 * t286 + t287 * t299;
t379 = -rSges(4,1) * t434 - rSges(4,2) * t435;
t359 = t379 - t588;
t349 = t359 - t590;
t380 = t435 * rSges(4,1) - t434 * rSges(4,2);
t360 = t380 + t587;
t350 = t360 + t589;
t242 = -t380 * t349 + t350 * t379;
t513 = t380 * t359 - t360 * t379;
t522 = t300 * t290 - t291 * t299;
t524 = t289 * t280 - t281 * t288;
t656 = m(6) / 0.2e1;
t657 = m(5) / 0.2e1;
t658 = m(4) / 0.2e1;
t497 = (t513 + t242) * t658 + (t524 + t146) * t656 + (t522 + t155) * t657;
t498 = (-t513 + t242) * t658 + (-t524 + t146) * t656 + (-t522 + t155) * t657;
t13 = t498 - t497;
t689 = t13 * qJD(1);
t430 = Icges(6,4) * t439;
t391 = -Icges(6,2) * t437 + t430;
t392 = Icges(6,1) * t437 + t430;
t688 = t391 + t392;
t441 = Icges(5,4) * t447;
t411 = -Icges(5,2) * t445 + t441;
t412 = Icges(5,1) * t445 + t441;
t687 = t412 + t411;
t394 = rSges(6,1) * t437 + rSges(6,2) * t439;
t586 = pkin(4) * t445;
t468 = t394 + t586;
t671 = t468 * t435;
t672 = t468 * t434;
t160 = t288 * t672 - t289 * t671;
t414 = rSges(5,1) * t445 + rSges(5,2) * t447;
t678 = t280 * t672 - t281 * t671;
t577 = ((-t291 + t300) * t435 + (t290 - t299) * t434) * t414 * t657 + (-t160 + t678) * t656;
t367 = t414 * t434;
t368 = t414 * t435;
t180 = t290 * t367 - t291 * t368;
t184 = t299 * t367 - t300 * t368;
t676 = (t184 + t180) * t657 + (t160 + t678) * t656;
t686 = t577 - t676;
t268 = t286 * t367;
t679 = t274 * t672 - t275 * t671;
t578 = (t268 + (-t299 * t434 + (-t287 + t300) * t435) * t414) * t657 + (-t160 + t679) * t656;
t176 = -t287 * t368 + t268;
t677 = (t184 + t176) * t657 + (t160 + t679) * t656;
t685 = t578 - t677;
t637 = t434 / 0.2e1;
t636 = -t435 / 0.2e1;
t684 = t435 / 0.2e1;
t635 = m(3) * (t589 * (-rSges(3,1) * t438 - rSges(3,2) * t440) + (t440 * rSges(3,1) - t438 * rSges(3,2)) * t590);
t431 = t434 ^ 2;
t432 = t435 ^ 2;
t500 = t431 + t432;
t638 = -t434 / 0.2e1;
t680 = t637 + t638;
t571 = Icges(6,4) * t437;
t390 = Icges(6,2) * t439 + t571;
t393 = Icges(6,1) * t439 - t571;
t487 = t688 * t439 / 0.2e1 + (-t390 / 0.2e1 + t393 / 0.2e1) * t437;
t357 = t394 * t434;
t249 = t274 * t357;
t358 = t394 * t435;
t159 = -t275 * t358 + t249;
t171 = t288 * t357 - t289 * t358;
t639 = m(6) * (t171 + t159);
t675 = t487 + t639 / 0.2e1;
t162 = t280 * t357 - t281 * t358;
t611 = m(6) * (t171 + t162);
t674 = t487 + t611 / 0.2e1;
t535 = t435 * t439;
t253 = t434 * t329 + t435 * (rSges(6,1) * t535 + t493);
t165 = -t434 * (pkin(3) * t434 - t429 + t502) + (-t434 * pkin(8) - t418 + (-pkin(3) + t436) * t435) * t435 + t253;
t277 = -t434 * t357 - t435 * t358;
t396 = -rSges(6,2) * t437 + t582;
t478 = Icges(6,5) * t437 + Icges(6,6) * t439;
t351 = t478 * t434;
t352 = t435 * t478;
t328 = Icges(6,5) * t434 + t393 * t435;
t507 = -t390 * t435 + t328;
t326 = Icges(6,6) * t434 + t391 * t435;
t509 = -t392 * t435 - t326;
t463 = -t437 * t507 + t439 * t509;
t400 = Icges(6,4) * t544;
t327 = Icges(6,1) * t543 - Icges(6,5) * t435 - t400;
t508 = -Icges(6,2) * t543 + t327 - t400;
t325 = Icges(6,4) * t543 - Icges(6,2) * t544 - Icges(6,6) * t435;
t510 = t392 * t434 + t325;
t464 = t437 * t508 + t439 * t510;
t584 = (-t431 * t352 + (t464 * t435 + (t351 + t463) * t434) * t435) * t637 + (-t432 * t351 + (t463 * t434 + (t352 + t464) * t435) * t434) * t636;
t29 = t584 + m(6) * (t165 * t277 + (t434 * t672 + t435 * t671) * t396);
t673 = t29 * qJD(5);
t139 = -t281 * t274 + t275 * t280;
t152 = -t291 * t286 + t287 * t290;
t220 = -t360 * t349 + t350 * t359;
t670 = m(5) * t184 + m(6) * t160;
t669 = qJD(1) + qJD(2) + qJD(3);
t579 = (t268 + (-t290 * t434 + (-t287 + t291) * t435) * t414) * t657 + (-t678 + t679) * t656;
t667 = (t180 + t176) * t657 + (t678 + t679) * t656;
t572 = Icges(5,4) * t445;
t410 = Icges(5,2) * t447 + t572;
t413 = Icges(5,1) * t447 - t572;
t666 = t687 * t447 / 0.2e1 + (t413 / 0.2e1 - t410 / 0.2e1) * t445;
t296 = t328 * t543;
t389 = Icges(6,5) * t439 - Icges(6,6) * t437;
t551 = t389 * t435;
t324 = Icges(6,3) * t434 + t551;
t492 = t435 * t324 - t296;
t186 = -t326 * t544 - t492;
t323 = Icges(6,5) * t543 - Icges(6,6) * t544 - Icges(6,3) * t435;
t516 = -t434 * t323 - t327 * t535;
t187 = -t325 * t536 - t516;
t515 = t434 * t324 + t328 * t535;
t188 = -t326 * t536 + t515;
t490 = t326 * t437 - t323;
t559 = t325 * t437;
t488 = (-t187 * t435 + t188 * t434) * t684 + ((t186 - t296 + (t324 + t559) * t435 + t516) * t435 + t515 * t434) * t636 + ((t434 * t490 + t186 + t187 + t492) * t434 + (t434 * (-t327 * t439 + t559) + t188 - t515 + (t323 + t490) * t435) * t435) * t637;
t664 = 0.4e1 * qJD(1);
t662 = 0.4e1 * qJD(2);
t661 = 0.2e1 * qJD(3);
t660 = 2 * qJD(4);
t641 = m(6) * (t162 + t159);
t631 = m(4) * t220;
t629 = m(4) * t242;
t628 = m(4) * t513;
t620 = m(5) * t152;
t618 = m(5) * t155;
t617 = m(5) * t522;
t615 = m(5) * t176;
t614 = m(5) * t180;
t610 = m(6) * (t249 + (-t280 * t434 + (-t275 + t281) * t435) * t394);
t609 = m(6) * (t249 + (-t288 * t434 + (-t275 + t289) * t435) * t394);
t608 = m(6) * ((-t281 + t289) * t435 + (t280 - t288) * t434) * t394;
t473 = (-t274 * t435 - t275 * t434) * t396;
t517 = -t357 * t671 + t358 * t672;
t607 = m(6) * (t473 + t517);
t470 = (t434 * t671 - t435 * t672) * t394;
t606 = m(6) * (t473 + t470);
t472 = (-t280 * t435 - t281 * t434) * t396;
t605 = m(6) * (t472 + t517);
t604 = m(6) * (t472 + t470);
t471 = (-t288 * t435 - t289 * t434) * t396;
t603 = m(6) * (t471 + t517);
t602 = m(6) * (t471 + t470);
t601 = m(6) * t139;
t599 = m(6) * t146;
t598 = m(6) * t524;
t596 = m(6) * t679;
t595 = m(6) * t678;
t594 = m(6) * t159;
t592 = m(6) * t162;
t591 = m(6) * t171;
t541 = t434 * t447;
t342 = Icges(5,4) * t541 - Icges(5,2) * t542 - Icges(5,6) * t435;
t556 = t342 * t445;
t409 = Icges(5,5) * t447 - Icges(5,6) * t445;
t548 = t409 * t435;
t533 = t435 * t447;
t340 = Icges(5,5) * t541 - Icges(5,6) * t542 - Icges(5,3) * t435;
t406 = Icges(5,4) * t542;
t344 = Icges(5,1) * t541 - Icges(5,5) * t435 - t406;
t512 = -t434 * t340 - t344 * t533;
t341 = Icges(5,3) * t434 + t548;
t345 = Icges(5,5) * t434 + t413 * t435;
t511 = t434 * t341 + t345 * t533;
t506 = t412 * t434 + t342;
t343 = Icges(5,6) * t434 + t411 * t435;
t505 = -t412 * t435 - t343;
t504 = -Icges(5,2) * t541 + t344 - t406;
t503 = -t410 * t435 + t345;
t494 = -t396 - t585;
t313 = t345 * t541;
t491 = t435 * t341 - t313;
t489 = t343 * t445 - t340;
t484 = t641 / 0.2e1 + t487;
t479 = Icges(5,5) * t445 + Icges(5,6) * t447;
t469 = (-t367 * t435 + t368 * t434) * t414;
t460 = (-t390 + t393) * t439 - t688 * t437;
t467 = -t488 + (t389 * t434 + t435 * t460 + t437 * t509 + t439 * t507) * t637 + (t434 * t460 - t437 * t510 + t439 * t508 - t551) * t636;
t466 = -t487 + t680 * (t439 * t325 + t437 * t327);
t465 = t487 + t666;
t462 = t445 * t504 + t447 * t506;
t461 = -t445 * t503 + t447 * t505;
t459 = (-t410 + t413) * t447 - t687 * t445;
t455 = (-t357 * t435 + t358 * t434) * t394;
t454 = t465 + t667;
t451 = t466 - t666 + t680 * (t447 * t342 + t445 * t344);
t207 = -t343 * t542 - t491;
t144 = -(-t434 * (-t344 * t447 + t556) - t435 * t340) * t435 + t207 * t434;
t208 = -t342 * t534 - t512;
t209 = -t343 * t534 + t511;
t145 = -t208 * t435 + t209 * t434;
t46 = (t435 * t489 + t209 - t511) * t435 + (t434 * t489 + t208 + t491) * t434;
t47 = (t207 - t313 + (t341 + t556) * t435 + t512) * t435 + t511 * t434;
t449 = (t47 * t684 + t467 + (t46 + t144) * t638 + (t409 * t434 + t435 * t459 + t445 * t505 + t447 * t503) * t637 + (t434 * t459 - t445 * t506 + t447 * t504 + t145 - t548) * t636) * qJD(4);
t416 = -rSges(5,2) * t445 + t583;
t362 = t435 * t479;
t361 = t479 * t434;
t339 = t494 * t435;
t337 = t494 * t434;
t250 = -t500 * t586 + t277;
t142 = t487 + t591;
t136 = t487 + t592;
t132 = t602 / 0.2e1;
t130 = t487 + t594;
t129 = t603 / 0.2e1;
t122 = t604 / 0.2e1;
t120 = t605 / 0.2e1;
t111 = t606 / 0.2e1;
t110 = t607 / 0.2e1;
t105 = t608 / 0.2e1;
t103 = t609 / 0.2e1;
t99 = t610 / 0.2e1;
t68 = t465 + t670;
t57 = t465 + t595 + t614;
t56 = t465 + t596 + t615;
t55 = -t598 - t617 - t628;
t54 = t599 + t618 + t629;
t48 = t601 + t620 + t631 + t635;
t37 = -t608 / 0.2e1 + t674;
t36 = t105 + t674;
t35 = -t609 / 0.2e1 + t675;
t34 = t103 + t675;
t33 = -t610 / 0.2e1 + t484;
t32 = t99 + t484;
t31 = m(6) * (t500 * t394 * t396 + t253 * t277) + t584;
t30 = t31 * qJD(5);
t28 = t105 - t611 / 0.2e1 + t466;
t27 = t103 - t639 / 0.2e1 + t466;
t26 = t99 - t641 / 0.2e1 + t466;
t25 = t465 + t577 + t676;
t24 = t465 - t686;
t23 = t465 - t685;
t22 = t465 + t578 + t677;
t21 = t454 - t579;
t20 = t454 + t579;
t18 = t488 * qJD(5);
t17 = t451 + t686;
t16 = t451 + t685;
t15 = t451 + t579 - t667;
t12 = t497 + t498;
t11 = t129 - t602 / 0.2e1 + t488;
t10 = t132 - t603 / 0.2e1 + t488;
t9 = t120 - t604 / 0.2e1 + t488;
t8 = t122 - t605 / 0.2e1 + t488;
t7 = t110 - t606 / 0.2e1 + t488;
t6 = t111 - t607 / 0.2e1 + t488;
t5 = t129 + t132 + t467;
t4 = t120 + t122 + t467;
t3 = t110 + t111 + t467;
t2 = (t145 / 0.2e1 - t47 / 0.2e1) * t435 + (t46 / 0.2e1 + t144 / 0.2e1) * t434 + t488;
t1 = t2 * qJD(4);
t14 = [qJD(2) * t48 + qJD(3) * t54 + qJD(4) * t56 + qJD(5) * t130, t48 * qJD(1) + t12 * qJD(3) + t20 * qJD(4) + t32 * qJD(5) + 0.2e1 * (t635 / 0.2e1 + t139 * t656 + t152 * t657 + t220 * t658) * qJD(2), t54 * qJD(1) + t12 * qJD(2) + t22 * qJD(4) + t34 * qJD(5) + (t146 * t656 + t155 * t657 + t242 * t658) * t661, t56 * qJD(1) + t20 * qJD(2) + t22 * qJD(3) + t449 + t3 * qJD(5) + ((t274 * t339 + t275 * t337) * t656 + ((-t286 * t435 - t287 * t434) * t416 + t469) * t657) * t660, t130 * qJD(1) + t32 * qJD(2) + t34 * qJD(3) + t3 * qJD(4) + ((t473 + t455) * m(6) + t467) * qJD(5); t13 * qJD(3) + t21 * qJD(4) + t33 * qJD(5) + (-t601 / 0.4e1 - t620 / 0.4e1 - t631 / 0.4e1 - t635 / 0.4e1) * t664, qJD(3) * t55 + qJD(4) * t57 + qJD(5) * t136, t689 + t55 * qJD(2) + t25 * qJD(4) + t36 * qJD(5) + (-t513 * t658 - t522 * t657 - t524 * t656) * t661, t21 * qJD(1) + t57 * qJD(2) + t25 * qJD(3) + t449 + t4 * qJD(5) + ((t280 * t339 + t281 * t337) * t656 + ((-t290 * t435 - t291 * t434) * t416 + t469) * t657) * t660, t33 * qJD(1) + t136 * qJD(2) + t36 * qJD(3) + t4 * qJD(4) + ((t472 + t455) * m(6) + t467) * qJD(5); -t13 * qJD(2) + t23 * qJD(4) + t35 * qJD(5) + (-t599 / 0.4e1 - t618 / 0.4e1 - t629 / 0.4e1) * t664, -t689 + t24 * qJD(4) + t37 * qJD(5) + (t598 / 0.4e1 + t617 / 0.4e1 + t628 / 0.4e1) * t662, qJD(4) * t68 + qJD(5) * t142, t23 * qJD(1) + t24 * qJD(2) + t68 * qJD(3) + t449 + t5 * qJD(5) + ((t288 * t339 + t289 * t337) * t656 + ((-t299 * t435 - t300 * t434) * t416 + t469) * t657) * t660, t35 * qJD(1) + t37 * qJD(2) + t142 * qJD(3) + t5 * qJD(4) + ((t471 + t455) * m(6) + t467) * qJD(5); t451 * qJD(1) + t15 * qJD(2) + t16 * qJD(3) + t1 + t7 * qJD(5) + (-t596 / 0.4e1 - t615 / 0.4e1) * t664, t15 * qJD(1) + t17 * qJD(3) + t1 + t9 * qJD(5) + (-t595 / 0.4e1 - t614 / 0.4e1) * t662 + t451 * qJD(2), t16 * qJD(1) + t17 * qJD(2) + t11 * qJD(5) + t1 + (t451 - t670) * qJD(3), (m(5) * ((t434 * (rSges(5,1) * t541 - t501) + t435 * (rSges(5,1) * t533 + t434 * rSges(5,3) - t408)) * (-t367 * t434 - t435 * t368) + t500 * t416 * t414) + (-t431 * t362 + (t462 * t435 + (t361 + t461) * t434) * t435) * t637 + (-t432 * t361 + (t461 * t434 + (t362 + t462) * t435) * t434) * t636 + m(6) * (t165 * t250 - t337 * t672 - t339 * t671) + t584) * qJD(4) + t673 + t669 * t2, t7 * qJD(1) + t9 * qJD(2) + t11 * qJD(3) + t29 * qJD(4) + t673; (t466 - t594) * qJD(1) + t26 * qJD(2) + t27 * qJD(3) + t6 * qJD(4) + t18, t26 * qJD(1) + (t466 - t592) * qJD(2) + t28 * qJD(3) + t8 * qJD(4) + t18, t27 * qJD(1) + t28 * qJD(2) + (t466 - t591) * qJD(3) + t10 * qJD(4) + t18, t6 * qJD(1) + t8 * qJD(2) + t10 * qJD(3) + ((t250 * t253 + (-t337 * t434 - t339 * t435) * t394) * m(6) + t584) * qJD(4) + t30, qJD(4) * t31 + t488 * t669 + t30;];
Cq = t14;
