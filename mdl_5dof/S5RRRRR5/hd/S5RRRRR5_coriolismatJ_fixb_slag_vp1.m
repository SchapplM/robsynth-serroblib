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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:13:04
% EndTime: 2020-01-03 12:13:16
% DurationCPUTime: 8.58s
% Computational Cost: add. (76859->495), mult. (43211->603), div. (0->0), fcn. (38666->10), ass. (0->311)
t443 = qJ(1) + qJ(2);
t440 = qJ(3) + t443;
t432 = sin(t440);
t446 = cos(qJ(4));
t586 = pkin(4) * t446;
t434 = pkin(3) + t586;
t433 = cos(t440);
t442 = qJ(4) + qJ(5);
t435 = sin(t442);
t532 = t432 * t435;
t492 = -rSges(6,2) * t532 - t433 * rSges(6,3);
t448 = -pkin(9) - pkin(8);
t521 = t433 * t448;
t437 = cos(t442);
t582 = rSges(6,1) * t437;
t292 = t521 + (t434 + t582) * t432 + t492;
t436 = sin(t443);
t588 = pkin(2) * t436;
t282 = t292 + t588;
t585 = sin(qJ(1)) * pkin(1);
t276 = t282 + t585;
t402 = t433 * t434;
t524 = t433 * t437;
t525 = t433 * t435;
t476 = rSges(6,1) * t524 - rSges(6,2) * t525;
t293 = t402 + (rSges(6,3) - t448) * t432 + t476;
t438 = cos(t443);
t431 = pkin(2) * t438;
t283 = t431 + t293;
t441 = cos(qJ(1)) * pkin(1);
t277 = t441 + t283;
t146 = t276 * t293 - t292 * t277;
t425 = t433 * pkin(8);
t444 = sin(qJ(4));
t530 = t432 * t444;
t491 = -rSges(5,2) * t530 - t433 * rSges(5,3);
t583 = rSges(5,1) * t446;
t305 = -t425 + (pkin(3) + t583) * t432 + t491;
t295 = t305 + t588;
t290 = t295 + t585;
t522 = t433 * t446;
t523 = t433 * t444;
t482 = rSges(5,1) * t522 - rSges(5,2) * t523 + t432 * rSges(5,3);
t488 = t433 * pkin(3) + t432 * pkin(8);
t306 = t482 + t488;
t296 = t431 + t306;
t291 = t441 + t296;
t155 = t290 * t306 - t305 * t291;
t380 = rSges(4,1) * t432 + rSges(4,2) * t433;
t360 = t380 + t588;
t350 = t360 + t585;
t381 = t433 * rSges(4,1) - rSges(4,2) * t432;
t361 = t381 + t431;
t351 = t361 + t441;
t242 = t350 * t381 - t380 * t351;
t503 = -t360 * t381 + t380 * t361;
t510 = -t295 * t306 + t305 * t296;
t512 = -t282 * t293 + t292 * t283;
t653 = m(6) / 0.2e1;
t654 = m(5) / 0.2e1;
t655 = m(4) / 0.2e1;
t484 = (t503 + t242) * t655 + (t512 + t146) * t653 + (t510 + t155) * t654;
t485 = (-t503 + t242) * t655 + (-t512 + t146) * t653 + (-t510 + t155) * t654;
t12 = t485 - t484;
t697 = t12 * qJD(1);
t395 = rSges(6,1) * t435 + rSges(6,2) * t437;
t359 = t395 * t433;
t587 = pkin(4) * t444;
t481 = t395 + t587;
t684 = t481 * t432;
t289 = t684 * t359;
t571 = Icges(6,4) * t435;
t394 = Icges(6,1) * t437 - t571;
t329 = -Icges(6,5) * t433 + t394 * t432;
t531 = t432 * t437;
t302 = t329 * t531;
t326 = Icges(6,5) * t524 - Icges(6,6) * t525 + Icges(6,3) * t432;
t401 = Icges(6,4) * t525;
t330 = Icges(6,1) * t524 + Icges(6,5) * t432 - t401;
t328 = Icges(6,4) * t524 - Icges(6,2) * t525 + Icges(6,6) * t432;
t552 = t328 * t435;
t464 = t330 * t437 - t552;
t696 = -t432 * t326 - t433 * t464 - t302;
t390 = Icges(6,5) * t437 - Icges(6,6) * t435;
t540 = t390 * t432;
t325 = -Icges(6,3) * t433 + t540;
t695 = t432 * t325 + t329 * t524;
t427 = Icges(6,4) * t437;
t392 = -Icges(6,2) * t435 + t427;
t675 = Icges(6,1) * t435 + t427;
t694 = t392 + t675;
t439 = Icges(5,4) * t446;
t413 = -Icges(5,2) * t444 + t439;
t674 = Icges(5,1) * t444 + t439;
t693 = t413 + t674;
t416 = rSges(5,1) * t444 + rSges(5,2) * t446;
t683 = t481 * t433;
t507 = t292 * t684 + t293 * t683;
t369 = t416 * t433;
t555 = t306 * t369;
t689 = -t282 * t684 - t283 * t683;
t577 = (t555 + (-t296 * t433 + (-t295 + t305) * t432) * t416) * t654 + (t507 + t689) * t653;
t368 = t416 * t432;
t180 = -t295 * t368 - t296 * t369;
t184 = -t305 * t368 - t555;
t679 = (t184 + t180) * t654 + (-t507 + t689) * t653;
t692 = t577 - t679;
t690 = -t276 * t684 - t277 * t683;
t578 = (t555 + (-t291 * t433 + (-t290 + t305) * t432) * t416) * t654 + (t507 + t690) * t653;
t176 = -t290 * t368 - t291 * t369;
t680 = (t184 + t176) * t654 + (-t507 + t690) * t653;
t691 = t578 - t680;
t636 = -t432 / 0.2e1;
t688 = t432 / 0.2e1;
t635 = -t433 / 0.2e1;
t633 = m(3) * (-t441 * (rSges(3,1) * t436 + rSges(3,2) * t438) + t585 * (t438 * rSges(3,1) - rSges(3,2) * t436));
t418 = -rSges(5,2) * t444 + t583;
t685 = t418 * t654;
t428 = t432 ^ 2;
t429 = t433 ^ 2;
t487 = t428 + t429;
t634 = t433 / 0.2e1;
t681 = t634 + t635;
t391 = Icges(6,2) * t437 + t571;
t477 = t694 * t437 / 0.2e1 + (-t391 / 0.2e1 + t394 / 0.2e1) * t435;
t358 = t395 * t432;
t159 = -t276 * t358 - t277 * t359;
t561 = t293 * t359;
t171 = -t292 * t358 - t561;
t637 = m(6) * (t171 + t159);
t678 = t477 + t637 / 0.2e1;
t162 = -t282 * t358 - t283 * t359;
t609 = m(6) * (t171 + t162);
t677 = t477 + t609 / 0.2e1;
t309 = t432 * (rSges(6,1) * t531 + t492);
t331 = t432 * rSges(6,3) + t476;
t166 = t309 + (-t432 * t448 + t331 + t402 - t488) * t433 + (t521 + t425 + (-pkin(3) + t434) * t432) * t432;
t279 = -t432 * t358 - t359 * t433;
t397 = -rSges(6,2) * t435 + t582;
t536 = t397 * t433;
t537 = t397 * t432;
t466 = Icges(6,5) * t435 + Icges(6,6) * t437;
t352 = t432 * t466;
t353 = t466 * t433;
t500 = -t391 * t432 + t329;
t327 = -Icges(6,6) * t433 + t392 * t432;
t502 = t432 * t675 + t327;
t667 = -t435 * t500 - t437 * t502;
t499 = -Icges(6,2) * t524 + t330 - t401;
t501 = t433 * t675 + t328;
t668 = -t435 * t499 - t437 * t501;
t584 = (-t429 * t352 + (t668 * t432 + (t353 - t667) * t433) * t432) * t635 + (t428 * t353 + (t667 * t433 + (-t352 - t668) * t432) * t433) * t636;
t29 = t584 + m(6) * (t166 * t279 + t536 * t683 + t537 * t684);
t676 = t29 * qJD(5);
t139 = t276 * t283 - t282 * t277;
t152 = t290 * t296 - t295 * t291;
t220 = t350 * t361 - t360 * t351;
t673 = m(5) * t184 - m(6) * t507;
t672 = qJD(1) + qJD(2) + qJD(3);
t579 = ((-t291 + t296) * t433 + (-t290 + t295) * t432) * t416 * t654 + (-t689 + t690) * t653;
t670 = (t180 + t176) * t654 + (t689 + t690) * t653;
t669 = t694 * t435 + (t391 - t394) * t437;
t572 = Icges(5,4) * t444;
t412 = Icges(5,2) * t446 + t572;
t415 = Icges(5,1) * t446 - t572;
t666 = t693 * t444 + (t412 - t415) * t446;
t407 = Icges(5,4) * t523;
t346 = Icges(5,1) * t522 + Icges(5,5) * t432 - t407;
t495 = -Icges(5,2) * t522 + t346 - t407;
t344 = Icges(5,4) * t522 - Icges(5,2) * t523 + Icges(5,6) * t432;
t497 = t433 * t674 + t344;
t665 = -t444 * t495 - t446 * t497;
t345 = -Icges(5,5) * t433 + t415 * t432;
t496 = -t412 * t432 + t345;
t343 = -Icges(5,6) * t433 + t413 * t432;
t498 = t432 * t674 + t343;
t664 = -t444 * t496 - t446 * t498;
t663 = t693 * t446 / 0.2e1 + (t415 / 0.2e1 - t412 / 0.2e1) * t444;
t185 = -t433 * t325 - t327 * t532 + t302;
t303 = t330 * t531;
t186 = t433 * t326 + t328 * t532 - t303;
t551 = t329 * t437;
t553 = t327 * t435;
t478 = (-t185 * t433 - t186 * t432) * t688 + ((t303 + (t325 - t552) * t432 - t695) * t432 + ((t325 + t464) * t433 + (t551 + t553) * t432 + t696) * t433) * t636 + ((-t186 + (t326 - t551) * t433 + t695) * t433 + (t185 + (t326 + t553) * t432 + t696) * t432) * t635;
t661 = 0.4e1 * qJD(1);
t659 = 0.4e1 * qJD(2);
t658 = 0.2e1 * qJD(3);
t657 = 2 * qJD(4);
t639 = m(6) * (t162 + t159);
t629 = m(4) * t220;
t627 = m(4) * t242;
t626 = m(4) * t503;
t618 = m(5) * t152;
t616 = m(5) * t155;
t615 = m(5) * t510;
t613 = m(5) * t176;
t612 = m(5) * t180;
t608 = m(6) * ((-t277 + t283) * t433 + (-t276 + t282) * t432) * t395;
t607 = m(6) * (t561 + (-t277 * t433 + (-t276 + t292) * t432) * t395);
t606 = m(6) * (t561 + (-t283 * t433 + (-t282 + t292) * t432) * t395);
t249 = t276 * t536;
t505 = -t358 * t683 + t289;
t605 = m(6) * (-t277 * t537 + t249 + t505);
t554 = t683 * t395;
t567 = t277 * t397;
t604 = m(6) * (t249 - t289 + (t554 - t567) * t432);
t251 = t282 * t536;
t603 = m(6) * (-t283 * t537 + t251 + t505);
t564 = t283 * t397;
t602 = m(6) * (t251 - t289 + (t554 - t564) * t432);
t264 = t292 * t536;
t601 = m(6) * (-t293 * t537 + t264 + t505);
t560 = t293 * t397;
t600 = m(6) * (t264 - t289 + (t554 - t560) * t432);
t599 = m(6) * t139;
t597 = m(6) * t146;
t596 = m(6) * t512;
t594 = m(6) * t690;
t593 = m(6) * t689;
t592 = m(6) * t159;
t590 = m(6) * t162;
t589 = m(6) * t171;
t548 = t343 * t444;
t547 = t344 * t444;
t546 = t345 * t446;
t545 = t359 * t395;
t411 = Icges(5,5) * t446 - Icges(5,6) * t444;
t535 = t411 * t432;
t529 = t432 * t446;
t480 = t397 + t586;
t479 = t358 * t359;
t472 = t639 / 0.2e1 + t477;
t467 = Icges(5,5) * t444 + Icges(5,6) * t446;
t462 = t346 * t446 - t547;
t460 = -t478 + (t669 * t433 + t435 * t501 - t437 * t499 - t540) * t636 + (-t390 * t433 - t669 * t432 - t435 * t502 + t437 * t500) * t635;
t459 = -t477 + t681 * (t328 * t437 + t330 * t435);
t458 = t477 + t663;
t454 = t458 + t670;
t451 = t459 - t663 + t681 * (t344 * t446 + t444 * t346);
t316 = t345 * t529;
t341 = -Icges(5,3) * t433 + t535;
t206 = -t433 * t341 - t343 * t530 + t316;
t317 = t346 * t529;
t342 = Icges(5,5) * t522 - Icges(5,6) * t523 + Icges(5,3) * t432;
t207 = t433 * t342 + t344 * t530 - t317;
t144 = -t206 * t433 - t207 * t432;
t318 = t343 * t523;
t208 = -t432 * t341 - t345 * t522 + t318;
t209 = t432 * t342 + t433 * t462;
t145 = -t208 * t433 - t209 * t432;
t46 = (t208 + t317 - t318 + (t341 - t547) * t432) * t432 + (-t316 - t209 + (t341 + t462) * t433 + (t546 + t548) * t432) * t433;
t47 = (-t207 + t318 + (t342 - t546) * t433) * t433 + (t206 - t316 + (t342 + t548) * t432) * t432;
t449 = (t460 + t46 * t688 + (-t411 * t433 - t666 * t432 - t444 * t498 + t446 * t496) * t635 + (t145 + t47) * t634 + (t666 * t433 + t444 * t497 - t446 * t495 + t144 - t535) * t636) * qJD(4);
t363 = t467 * t433;
t362 = t432 * t467;
t340 = t480 * t433;
t338 = t480 * t432;
t255 = t433 * t331 + t309;
t250 = -t487 * t587 + t279;
t142 = t477 + t589;
t136 = t477 + t590;
t132 = t600 / 0.2e1;
t130 = t477 + t592;
t129 = t601 / 0.2e1;
t122 = t602 / 0.2e1;
t120 = t603 / 0.2e1;
t111 = t604 / 0.2e1;
t110 = t605 / 0.2e1;
t105 = t606 / 0.2e1;
t103 = t607 / 0.2e1;
t99 = t608 / 0.2e1;
t68 = t458 + t673;
t57 = t458 + t593 + t612;
t56 = t458 + t594 + t613;
t55 = -t596 - t615 - t626;
t54 = t597 + t616 + t627;
t48 = t599 + t618 + t629 + t633;
t37 = -t606 / 0.2e1 + t677;
t36 = t105 + t677;
t35 = -t607 / 0.2e1 + t678;
t34 = t103 + t678;
t33 = -t608 / 0.2e1 + t472;
t32 = t99 + t472;
t31 = m(6) * (t395 * t397 * t487 + t255 * t279) + t584;
t30 = t31 * qJD(5);
t28 = t105 - t609 / 0.2e1 + t459;
t27 = t103 - t637 / 0.2e1 + t459;
t26 = t99 - t639 / 0.2e1 + t459;
t25 = t458 - t692;
t24 = t458 + t577 + t679;
t23 = t458 + t578 + t680;
t22 = t458 - t691;
t21 = t454 + t579;
t20 = t454 - t579;
t18 = t478 * qJD(5);
t17 = t451 + t692;
t16 = t451 + t691;
t15 = t451 + t579 - t670;
t13 = t484 + t485;
t11 = t129 - t600 / 0.2e1 + t478;
t10 = t132 - t601 / 0.2e1 + t478;
t9 = t120 - t602 / 0.2e1 + t478;
t8 = t122 - t603 / 0.2e1 + t478;
t7 = t110 - t604 / 0.2e1 + t478;
t6 = t111 - t605 / 0.2e1 + t478;
t5 = t129 + t132 + t460;
t4 = t120 + t122 + t460;
t3 = t110 + t111 + t460;
t2 = (-t47 / 0.2e1 - t145 / 0.2e1) * t433 + (t144 / 0.2e1 - t46 / 0.2e1) * t432 + t478;
t1 = t2 * qJD(4);
t14 = [qJD(2) * t48 + qJD(3) * t54 + qJD(4) * t56 + qJD(5) * t130, t48 * qJD(1) + t13 * qJD(3) + t21 * qJD(4) + t32 * qJD(5) + 0.2e1 * (t633 / 0.2e1 + t139 * t653 + t152 * t654 + t220 * t655) * qJD(2), t54 * qJD(1) + t13 * qJD(2) + t23 * qJD(4) + t34 * qJD(5) + (t146 * t653 + t155 * t654 + t242 * t655) * t658, t56 * qJD(1) + t21 * qJD(2) + t23 * qJD(3) + t449 + t3 * qJD(5) + ((t290 * t433 - t291 * t432) * t685 + (t276 * t340 - t277 * t338) * t653) * t657, t130 * qJD(1) + t32 * qJD(2) + t34 * qJD(3) + t3 * qJD(4) + ((t249 + (t545 - t567) * t432 - t479) * m(6) + t460) * qJD(5); t12 * qJD(3) + t20 * qJD(4) + t33 * qJD(5) + (-t633 / 0.4e1 - t629 / 0.4e1 - t618 / 0.4e1 - t599 / 0.4e1) * t661, qJD(3) * t55 + qJD(4) * t57 + qJD(5) * t136, t697 + t55 * qJD(2) + t24 * qJD(4) + t36 * qJD(5) + (-t503 * t655 - t510 * t654 - t512 * t653) * t658, t20 * qJD(1) + t57 * qJD(2) + t24 * qJD(3) + t449 + t4 * qJD(5) + ((t282 * t340 - t283 * t338) * t653 + (t295 * t433 - t296 * t432) * t685) * t657, t33 * qJD(1) + t136 * qJD(2) + t36 * qJD(3) + t4 * qJD(4) + ((t251 + (t545 - t564) * t432 - t479) * m(6) + t460) * qJD(5); -t12 * qJD(2) + t22 * qJD(4) + t35 * qJD(5) + (-t627 / 0.4e1 - t616 / 0.4e1 - t597 / 0.4e1) * t661, -t697 + t25 * qJD(4) + t37 * qJD(5) + (t596 / 0.4e1 + t615 / 0.4e1 + t626 / 0.4e1) * t659, qJD(4) * t68 + qJD(5) * t142, t22 * qJD(1) + t25 * qJD(2) + t68 * qJD(3) + t449 + t5 * qJD(5) + ((t292 * t340 - t293 * t338) * t653 + (t305 * t433 - t306 * t432) * t685) * t657, t35 * qJD(1) + t37 * qJD(2) + t142 * qJD(3) + t5 * qJD(4) + ((t264 + (t545 - t560) * t432 - t479) * m(6) + t460) * qJD(5); t451 * qJD(1) + t15 * qJD(2) + t16 * qJD(3) + t1 + t7 * qJD(5) + (-t613 / 0.4e1 - t594 / 0.4e1) * t661, t15 * qJD(1) + t17 * qJD(3) + t1 + t9 * qJD(5) + (-t593 / 0.4e1 - t612 / 0.4e1) * t659 + t451 * qJD(2), t16 * qJD(1) + t17 * qJD(2) + t11 * qJD(5) + t1 + (t451 - t673) * qJD(3), (m(5) * ((t433 * t482 + t432 * (rSges(5,1) * t529 + t491)) * (-t368 * t432 - t369 * t433) + t487 * t418 * t416) + (-t429 * t362 + (t665 * t432 + (t363 - t664) * t433) * t432) * t635 + (t428 * t363 + (t664 * t433 + (-t362 - t665) * t432) * t433) * t636 + m(6) * (t166 * t250 + t338 * t684 + t340 * t683) + t584) * qJD(4) + t676 + t672 * t2, t7 * qJD(1) + t9 * qJD(2) + t11 * qJD(3) + t29 * qJD(4) + t676; (t459 - t592) * qJD(1) + t26 * qJD(2) + t27 * qJD(3) + t6 * qJD(4) + t18, t26 * qJD(1) + (t459 - t590) * qJD(2) + t28 * qJD(3) + t8 * qJD(4) + t18, t27 * qJD(1) + t28 * qJD(2) + (t459 - t589) * qJD(3) + t10 * qJD(4) + t18, t6 * qJD(1) + t8 * qJD(2) + t10 * qJD(3) + ((t250 * t255 + (t338 * t432 + t340 * t433) * t395) * m(6) + t584) * qJD(4) + t30, qJD(4) * t31 + t478 * t672 + t30;];
Cq = t14;
