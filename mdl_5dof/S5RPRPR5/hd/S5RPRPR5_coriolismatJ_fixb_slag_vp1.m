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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:41:29
% EndTime: 2020-01-03 11:42:24
% DurationCPUTime: 20.78s
% Computational Cost: add. (46095->622), mult. (50999->848), div. (0->0), fcn. (54935->10), ass. (0->365)
t721 = Icges(4,3) + Icges(5,3);
t465 = qJ(3) + pkin(9);
t455 = sin(t465);
t466 = sin(pkin(8));
t467 = cos(pkin(8));
t456 = cos(t465);
t600 = Icges(5,4) * t456;
t363 = -Icges(5,6) * t467 + (-Icges(5,2) * t455 + t600) * t466;
t601 = Icges(5,4) * t455;
t364 = -Icges(5,5) * t467 + (Icges(5,1) * t456 - t601) * t466;
t469 = sin(qJ(3));
t471 = cos(qJ(3));
t603 = Icges(4,4) * t471;
t393 = -Icges(4,6) * t467 + (-Icges(4,2) * t469 + t603) * t466;
t604 = Icges(4,4) * t469;
t394 = -Icges(4,5) * t467 + (Icges(4,1) * t471 - t604) * t466;
t396 = (-Icges(5,2) * t456 - t601) * t466;
t397 = (-Icges(5,1) * t455 - t600) * t466;
t418 = (-Icges(4,2) * t471 - t604) * t466;
t419 = (-Icges(4,1) * t469 - t603) * t466;
t720 = ((t364 / 0.2e1 + t396 / 0.2e1) * t455 - (t397 / 0.2e1 - t363 / 0.2e1) * t456 + (t394 / 0.2e1 + t418 / 0.2e1) * t469 - (t419 / 0.2e1 - t393 / 0.2e1) * t471) * t466;
t461 = t471 * pkin(3);
t431 = pkin(4) * t456 + t461;
t428 = pkin(2) + t431;
t472 = cos(qJ(1));
t459 = t472 * qJ(2);
t470 = sin(qJ(1));
t457 = qJ(5) + t465;
t451 = cos(t457);
t450 = sin(t457);
t577 = t470 * t450;
t379 = -t451 * t472 - t467 * t577;
t570 = t472 * t450;
t576 = t470 * t451;
t380 = t467 * t576 - t570;
t494 = rSges(6,1) * t380 + rSges(6,2) * t379;
t615 = pkin(3) * t469;
t430 = pkin(4) * t455 + t615;
t426 = t472 * t430;
t468 = -qJ(4) - pkin(6);
t463 = -pkin(7) + t468;
t584 = t466 * t470;
t530 = -t463 * t584 - t426;
t199 = -t459 + (rSges(6,3) * t466 + t428 * t467 + pkin(1)) * t470 + t494 + t530;
t307 = t379 * rSges(6,1) - t380 * rSges(6,2);
t381 = t467 * t570 - t576;
t579 = t467 * t472;
t382 = t451 * t579 + t577;
t308 = t381 * rSges(6,1) + t382 * rSges(6,2);
t583 = t466 * t472;
t280 = t382 * rSges(6,1) - t381 * rSges(6,2) + rSges(6,3) * t583;
t526 = t472 * pkin(1) + t470 * qJ(2);
t425 = t470 * t430;
t531 = -t428 * t579 - t425;
t671 = -t463 * t583 + t280 + t526 - t531;
t719 = m(6) * (t199 * t307 - t308 * t671);
t525 = t463 - t468;
t573 = t470 * t469;
t446 = pkin(3) * t573;
t454 = t461 + pkin(2);
t528 = t454 * t579 + t446;
t243 = t525 * t583 + t528 + t531;
t261 = t467 * t280;
t612 = pkin(6) + t468;
t616 = pkin(2) * t467;
t332 = (t612 * t466 + t616) * t472 - t528;
t313 = t467 * t332;
t360 = -rSges(6,3) * t467 + (rSges(6,1) * t451 - rSges(6,2) * t450) * t466;
t613 = -pkin(2) + t454;
t361 = t613 * t466 + t612 * t467;
t529 = t428 - t454;
t514 = t466 * t529 + t467 * t525 + t360 + t361;
t488 = t514 * t583;
t104 = t467 * t243 - t261 + t313 - t488;
t105 = (t243 - t280 + t332) * t467 - t488;
t567 = t104 - t105;
t718 = m(6) * t567;
t575 = t470 * t455;
t399 = -t456 * t472 - t467 * t575;
t569 = t472 * t455;
t574 = t470 * t456;
t400 = t467 * t574 - t569;
t568 = t472 * t469;
t572 = t470 * t471;
t421 = t467 * t572 - t568;
t571 = t471 * t472;
t479 = t467 * t573 + t571;
t717 = Icges(4,5) * t421 + Icges(5,5) * t400 - Icges(4,6) * t479 + Icges(5,6) * t399 + t721 * t584;
t401 = t467 * t569 - t574;
t402 = t456 * t579 + t575;
t517 = t467 * t568;
t422 = t517 - t572;
t423 = t467 * t571 + t573;
t716 = Icges(4,5) * t423 + Icges(5,5) * t402 - Icges(4,6) * t422 - Icges(5,6) * t401 + t721 * t583;
t464 = t466 ^ 2;
t270 = Icges(6,5) * t380 + Icges(6,6) * t379 + Icges(6,3) * t584;
t715 = t270 * t583;
t714 = -t360 * t583 - t261;
t271 = Icges(6,5) * t382 - Icges(6,6) * t381 + Icges(6,3) * t583;
t713 = t271 * t583;
t386 = Icges(5,4) * t402;
t290 = Icges(5,2) * t401 - Icges(5,6) * t583 - t386;
t385 = Icges(5,4) * t401;
t292 = Icges(5,1) * t402 + Icges(5,5) * t583 - t385;
t412 = Icges(4,4) * t423;
t327 = Icges(4,2) * t422 - Icges(4,6) * t583 - t412;
t411 = Icges(4,4) * t422;
t329 = Icges(4,1) * t423 + Icges(4,5) * t583 - t411;
t712 = -t399 * t290 + t400 * t292 + t327 * t479 + t421 * t329;
t711 = -t466 / 0.2e1;
t598 = Icges(6,4) * t450;
t359 = -Icges(6,5) * t467 + (Icges(6,1) * t451 - t598) * t466;
t539 = t359 + (-Icges(6,2) * t451 - t598) * t466;
t710 = t450 * t539;
t687 = t716 * t584 + t712;
t709 = t716 * t470;
t708 = t717 * t472;
t334 = t423 * rSges(4,1) - t422 * rSges(4,2) + rSges(4,3) * t583;
t398 = -rSges(4,3) * t467 + (rSges(4,1) * t471 - rSges(4,2) * t469) * t466;
t706 = -t467 * t334 - t398 * t583;
t705 = t717 * t583;
t704 = t716 * t583;
t368 = Icges(6,4) * t382;
t275 = Icges(6,2) * t381 - Icges(6,6) * t583 - t368;
t367 = Icges(6,4) * t381;
t277 = Icges(6,1) * t382 + Icges(6,5) * t583 - t367;
t703 = -t379 * t275 + t380 * t277;
t661 = m(4) / 0.2e1;
t659 = m(5) / 0.2e1;
t657 = m(6) / 0.2e1;
t680 = t466 / 0.2e1;
t115 = -t271 * t584 - t703;
t591 = t271 * t470;
t592 = t270 * t472;
t667 = t472 ^ 2;
t689 = t466 * t667;
t697 = (t271 * t689 - t713 * t472 + (t115 + t703 + (t591 + t592) * t466 - t715) * t470) * t466;
t602 = Icges(5,4) * t400;
t288 = Icges(5,2) * t399 + Icges(5,6) * t584 + t602;
t384 = Icges(5,4) * t399;
t291 = Icges(5,1) * t400 + Icges(5,5) * t584 + t384;
t605 = Icges(4,4) * t421;
t325 = -Icges(4,2) * t479 + Icges(4,6) * t584 + t605;
t410 = Icges(4,4) * t479;
t328 = Icges(4,1) * t421 + Icges(4,5) * t584 - t410;
t688 = t399 * t288 + t400 * t291 - t325 * t479 + t421 * t328 + t717 * t584;
t295 = t402 * rSges(5,1) - t401 * rSges(5,2) + rSges(5,3) * t583;
t102 = t104 * t470;
t448 = pkin(3) * t568;
t527 = t468 * t584 + t448;
t242 = t467 * t470 * t529 + t527 + t530;
t279 = rSges(6,3) * t584 + t494;
t614 = pkin(6) * t466;
t331 = (t613 * t467 - t614) * t470 - t527;
t103 = (t242 + t279 + t331) * t467 + t514 * t584;
t538 = t361 - rSges(5,3) * t467 + (rSges(5,1) * t456 - rSges(5,2) * t455) * t466;
t497 = t538 * t583;
t165 = -t295 * t467 + t313 - t497;
t162 = t165 * t470;
t495 = rSges(5,1) * t400 + rSges(5,2) * t399;
t294 = rSges(5,3) * t584 + t495;
t164 = (t294 + t331) * t467 + t538 * t584;
t496 = rSges(4,1) * t421 - rSges(4,2) * t479;
t333 = rSges(4,3) * t584 + t496;
t233 = t333 * t467 + t398 * t584;
t521 = (t233 * t472 + t706 * t470) * t661 + (t103 * t472 + t102) * t657 + (t164 * t472 + t162) * t659;
t166 = (-t295 + t332) * t467 - t497;
t523 = (-t105 * t470 + t102) * t657 + (-t166 * t470 + t162) * t659;
t685 = t521 - t523;
t61 = -t103 * t584 + t104 * t583;
t87 = -t164 * t584 + t165 * t583;
t610 = t61 * t657 + t659 * t87;
t611 = ((t103 * t470 - t105 * t472) * t466 + t61) * t657 + ((t164 * t470 - t166 * t472) * t466 + t87) * t659;
t684 = t610 - t611;
t314 = Icges(5,5) * t399 - Icges(5,6) * t400;
t315 = Icges(5,5) * t401 + Icges(5,6) * t402;
t341 = -Icges(4,5) * t479 - Icges(4,6) * t421;
t342 = Icges(4,5) * t422 + Icges(4,6) * t423;
t683 = -(t315 + t342) * t583 + (t341 + t314) * t584;
t376 = (-Icges(6,5) * t450 - Icges(6,6) * t451) * t466;
t582 = t467 * t376;
t682 = -t582 / 0.2e1 + t710 * t711;
t634 = -t467 / 0.2e1;
t413 = t479 * pkin(3);
t194 = t307 * t583 + t308 * t584;
t383 = (-rSges(6,1) * t450 - rSges(6,2) * t451) * t466;
t221 = -t307 * t467 - t383 * t584;
t284 = t467 * t308;
t222 = -t383 * t583 + t284;
t597 = Icges(6,4) * t451;
t358 = -Icges(6,6) * t467 + (-Icges(6,2) * t450 + t597) * t466;
t378 = (-Icges(6,1) * t450 - t597) * t466;
t540 = t358 - t378;
t119 = t376 * t584 + t379 * t539 - t380 * t540;
t120 = -t376 * t583 + t381 * t539 + t382 * t540;
t299 = Icges(6,5) * t379 - Icges(6,6) * t380;
t300 = Icges(6,5) * t381 + Icges(6,6) * t382;
t507 = -t583 / 0.2e1;
t508 = t584 / 0.2e1;
t552 = Icges(6,2) * t382 - t277 + t367;
t366 = Icges(6,4) * t379;
t276 = Icges(6,1) * t380 + Icges(6,5) * t584 + t366;
t553 = -Icges(6,2) * t380 + t276 + t366;
t554 = -Icges(6,1) * t381 + t275 - t368;
t599 = Icges(6,4) * t380;
t273 = Icges(6,2) * t379 + Icges(6,6) * t584 + t599;
t555 = -Icges(6,1) * t379 + t273 + t599;
t163 = -t582 + (-t451 * t540 - t710) * t466;
t596 = t163 * t467;
t95 = -t299 * t467 + (-t450 * t553 - t451 * t555) * t466;
t96 = -t300 * t467 + (-t450 * t552 - t451 * t554) * t466;
t522 = (-t120 * t467 + (-t299 * t583 + t381 * t553 + t382 * t555) * t584 - (-t300 * t583 + t381 * t552 + t382 * t554) * t583) * t507 + (-t119 * t467 + (t299 * t584 + t379 * t553 - t380 * t555) * t584 - (t300 * t584 + t379 * t552 - t380 * t554) * t583) * t508 + (-t596 + (t470 * t95 - t472 * t96) * t466) * t634;
t185 = t279 * t583 - t280 * t584;
t547 = t331 * t583 + t332 * t584;
t81 = (t242 * t472 + t243 * t470) * t466 + t547 + t185;
t6 = t522 + m(6) * (-t103 * t221 + t104 * t222 + t81 * t194);
t678 = t6 * qJD(5);
t532 = t394 + t418;
t533 = t393 - t419;
t536 = t364 + t396;
t537 = t363 - t397;
t417 = (-Icges(4,5) * t469 - Icges(4,6) * t471) * t466;
t580 = t467 * t417;
t395 = (-Icges(5,5) * t455 - Icges(5,6) * t456) * t466;
t581 = t467 * t395;
t677 = t580 + t581 + (t455 * t536 + t456 * t537 + t469 * t532 + t471 * t533) * t466;
t670 = -t468 * t583 + t295 + t526 + t528;
t669 = (t614 + t616) * t472 + t526 + t334;
t499 = -t467 * t425 - t431 * t472;
t223 = t499 + t307;
t578 = t470 * t431;
t224 = -t426 * t467 - t308 + t578;
t321 = t401 * rSges(5,1) + t402 * rSges(5,2);
t447 = pkin(3) * t572;
t414 = pkin(3) * t517 - t447;
t260 = -t321 - t414;
t320 = t399 * rSges(5,1) - t400 * rSges(5,2);
t545 = -t320 + t413;
t565 = (m(6) * (-t223 * t472 + t224 * t470) + m(5) * (t260 * t470 + t472 * t545)) * t680;
t347 = -rSges(4,1) * t479 - rSges(4,2) * t421;
t348 = rSges(4,1) * t422 + rSges(4,2) * t423;
t515 = (-t470 * t223 - t224 * t472) * t657 + (-t260 * t472 + t470 * t545) * t659 + (-t470 * t347 + t348 * t472) * t661;
t427 = (t470 ^ 2 + t667) * t466;
t666 = 0.2e1 * t427;
t665 = 2 * qJD(1);
t664 = 4 * qJD(1);
t663 = 2 * qJD(3);
t662 = 4 * qJD(3);
t658 = m(5) / 0.4e1;
t656 = m(6) / 0.4e1;
t564 = t165 - t166;
t653 = m(5) * t164 * t564;
t218 = -t459 + (rSges(5,3) * t466 + t454 * t467 + pkin(1)) * t470 + t495 - t527;
t648 = m(5) * (-t218 * t545 + t260 * t670);
t213 = t279 * t467 + t360 * t584;
t647 = t213 * t718;
t646 = t103 * t718;
t566 = t221 * t199 + t222 * t671;
t642 = m(6) * (-t103 * t307 - t104 * t308 + t566);
t641 = m(6) * (-t213 * t223 + t224 * t714 + t566);
t133 = -t213 * t584 + t583 * t714;
t638 = m(6) * ((t213 * t470 - t472 * t714) * t466 + t133);
t636 = m(6) * (t199 * t223 + t224 * t671);
t633 = m(3) * (-(-rSges(3,2) * t584 - t472 * rSges(3,3) + pkin(1) * t470 - t459) * t472 + (-rSges(3,2) * t583 + t470 * rSges(3,3) + t526) * t470);
t248 = -t459 + (t616 + pkin(1) + (rSges(4,3) + pkin(6)) * t466) * t470 + t496;
t632 = m(4) * (t248 * t347 - t348 * t669);
t630 = m(4) * (-t248 * t472 + t669 * t470);
t144 = t218 * t584 + t583 * t670;
t628 = m(5) * t144;
t627 = m(5) * (-t218 * t472 + t670 * t470);
t113 = t199 * t584 + t583 * t671;
t624 = m(6) * t113;
t623 = m(6) * (-t199 * t472 + t671 * t470);
t622 = m(6) * t133;
t621 = m(6) * (t213 * t472 + t714 * t470);
t618 = m(6) * t194;
t617 = m(6) * (-t470 * t307 + t308 * t472);
t609 = m(6) * qJD(3);
t608 = m(6) * qJD(5);
t595 = ((-Icges(6,3) * t467 + (Icges(6,5) * t451 - Icges(6,6) * t450) * t466) * t584 + t358 * t379 + t359 * t380) * t467;
t585 = t451 * t466;
t551 = -Icges(5,1) * t399 + t288 + t602;
t550 = -Icges(5,1) * t401 + t290 - t386;
t549 = -Icges(5,2) * t400 + t291 + t384;
t548 = Icges(5,2) * t402 - t292 + t385;
t544 = Icges(4,1) * t479 + t325 + t605;
t543 = -Icges(4,1) * t422 + t327 - t412;
t542 = -Icges(4,2) * t421 + t328 - t410;
t541 = Icges(4,2) * t423 - t329 + t411;
t535 = -t413 * t583 + t414 * t584;
t534 = t467 * t414 + t464 * t448;
t244 = (t659 + t657) * t666;
t524 = t244 * qJD(1);
t114 = t270 * t584 + t379 * t273 + t380 * t276;
t511 = -t585 / 0.2e1;
t510 = t585 / 0.2e1;
t506 = t583 / 0.2e1;
t501 = -t430 + t615;
t500 = t466 * (-t501 * t466 - t383);
t24 = -t595 + ((t114 + t713) * t470 + ((-t591 + t592) * t466 - t115 - t715) * t472) * t466;
t59 = -t595 + (t114 * t470 - t115 * t472) * t466;
t498 = t24 * t507 + t59 * t506 + t697 * t508;
t493 = t358 * t511 + t378 * t510 + t682;
t492 = t647 / 0.2e1 + t498;
t487 = (((t708 - t709) * t466 + t687 - t705) * t472 + (t688 + t704) * t470) * t711 + (t470 * t688 + t472 * t687) * t680;
t486 = (t716 * t689 - t704 * t472 + ((t708 + t709) * t466 - t705 - t687 + t712) * t470) * t680 + (-t401 * t363 + t402 * t364 - t422 * t393 + t423 * t394 + (-t721 * t467 + (Icges(4,5) * t471 + Icges(5,5) * t456 - Icges(4,6) * t469 - Icges(5,6) * t455) * t466) * t583) * (t467 / 0.2e1 + t634);
t478 = t24 * t506 - t697 * t584 / 0.2e1 + (t119 + t95) * t508 + (t59 + t120 + t96) * t507;
t475 = t358 * t510 + t378 * t511 - t682;
t473 = t478 - t596;
t433 = t464 * t446;
t424 = (-rSges(4,1) * t469 - rSges(4,2) * t471) * t466;
t403 = (-rSges(5,1) * t455 - rSges(5,2) * t456) * t466;
t298 = -t501 * t579 + t447 - t578;
t297 = t499 + t413;
t258 = t467 * t348 - t424 * t583;
t257 = -t347 * t467 - t424 * t584;
t245 = (t658 + t656) * t666 - (m(5) + m(6)) * t427 / 0.2e1;
t203 = t467 * t321 - t403 * t583 + t534;
t202 = -t403 * t584 + t467 * t545 + t433;
t197 = t617 / 0.2e1;
t191 = -t618 / 0.2e1;
t172 = (t320 * t472 + t321 * t470) * t466 + t535;
t168 = -t417 * t583 + t422 * t532 + t423 * t533;
t167 = t417 * t584 - t421 * t533 - t479 * t532;
t154 = -t221 * t470 - t222 * t472;
t149 = t154 * t608;
t139 = t467 * t298 + t472 * t500 + t284 + t534;
t138 = t433 + t470 * t500 + (-t297 - t307 + t413) * t467;
t136 = t621 / 0.2e1;
t135 = -t395 * t583 + t401 * t536 + t402 * t537;
t134 = t395 * t584 + t399 * t536 - t400 * t537;
t125 = t622 / 0.2e1;
t122 = -t342 * t467 + (-t469 * t541 - t471 * t543) * t466;
t121 = -t341 * t467 + (-t469 * t542 - t471 * t544) * t466;
t109 = (t297 * t472 + t298 * t470) * t466 + t535 + t194;
t108 = -t315 * t467 + (-t455 * t548 - t456 * t550) * t466;
t107 = -t314 * t467 + (-t455 * t549 - t456 * t551) * t466;
t106 = -t467 * t194 + (-t221 * t472 + t222 * t470) * t466;
t101 = t106 * t608;
t74 = t624 + t628;
t71 = t493 + t719;
t67 = t638 / 0.2e1;
t53 = t623 + t627 + t630 + t633;
t49 = t136 - t617 / 0.2e1;
t48 = t197 + t136;
t47 = t197 - t621 / 0.2e1;
t45 = t641 / 0.2e1;
t43 = t125 + t67 + t618 / 0.2e1;
t42 = t191 + t125 - t638 / 0.2e1;
t41 = t191 + t67 - t622 / 0.2e1;
t35 = t642 / 0.2e1;
t28 = (-t417 / 0.2e1 - t395 / 0.2e1) * t467 + t632 + t648 + t636 - t720 + t493;
t14 = m(6) * (t185 * t194 - t213 * t221 + t222 * t714) + t522;
t13 = t14 * qJD(5);
t12 = t521 + t523 - t515;
t11 = t515 + t685;
t10 = t515 - t685;
t9 = t565 + t684;
t8 = t610 + t611 - t565;
t7 = t565 - t684;
t4 = t45 - t642 / 0.2e1 + t492;
t3 = t35 - t641 / 0.2e1 + t492;
t2 = t35 + t45 - t647 / 0.2e1 + t473;
t1 = t653 + t646 + (t470 * t486 + t472 * t487) * t466 + t498;
t5 = [t53 * qJD(2) + t28 * qJD(3) + t74 * qJD(4) + t71 * qJD(5), qJD(1) * t53 + qJD(3) * t11 + qJD(4) * t245 + qJD(5) * t48, t28 * qJD(1) + t11 * qJD(2) + t9 * qJD(4) + t2 * qJD(5) + (-t653 / 0.4e1 - t646 / 0.4e1) * t662 + ((-t233 * t347 + t248 * t257 + t258 * t669 - t348 * t706) * t661 + (t164 * t545 + t165 * t260 + t202 * t218 + t203 * t670) * t659 + (-t103 * t223 + t104 * t224 + t138 * t199 + t139 * t671) * t657) * t663 + (t478 + (-t163 + t677) * t467 + ((-t122 / 0.2e1 - t108 / 0.2e1 - t168 / 0.2e1 - t135 / 0.2e1 - t487) * t472 + (t167 / 0.2e1 + t134 / 0.2e1 + t121 / 0.2e1 + t107 / 0.2e1 - t486) * t470) * t466) * qJD(3), qJD(1) * t74 + qJD(2) * t245 + qJD(3) * t9 + qJD(5) * t42, t71 * qJD(1) + t48 * qJD(2) + t2 * qJD(3) + t42 * qJD(4) + ((-t213 * t307 - t308 * t714 + t566) * m(6) + t473) * qJD(5); t10 * qJD(3) - t244 * qJD(4) + t47 * qJD(5) + (-t633 / 0.4e1 - t630 / 0.4e1 - t627 / 0.4e1 - t623 / 0.4e1) * t664, 0, t10 * qJD(1) + t149 + ((-t257 * t470 - t258 * t472) * t661 + (-t202 * t470 - t203 * t472) * t659 + (-t138 * t470 - t139 * t472) * t657) * t663, -t524, t47 * qJD(1) + t154 * t609 + t149; t12 * qJD(2) + t1 * qJD(3) + t8 * qJD(4) + t3 * qJD(5) + (-t632 / 0.4e1 - t648 / 0.4e1 - t636 / 0.4e1) * t664 + (-t199 * t567 * t657 - t218 * t564 * t659) * t665 + (t580 / 0.2e1 + t581 / 0.2e1 + t475 + t720) * qJD(1), t12 * qJD(1), t1 * qJD(1) + (t522 + (t677 * t467 + ((-t108 - t122) * t472 + (t107 + t121) * t470) * t466) * t634 + ((-t399 * t548 + t400 * t550 + t421 * t543 + t479 * t541) * t583 + (-t167 - t134) * t467 + (t399 * t549 - t400 * t551 - t421 * t544 - t479 * t542 + t683) * t584) * t508 + ((t401 * t549 + t402 * t551 + t422 * t542 + t423 * t544) * t584 + (-t135 - t168) * t467 + (-t401 * t548 - t402 * t550 - t422 * t541 - t423 * t543 - t683) * t583) * t507) * qJD(3) + t678 + ((-t103 * t138 + t104 * t139 + t109 * t81) * t656 + (((t294 * t472 - t295 * t470) * t466 + t547) * t172 - t164 * t202 + t165 * t203) * t658 + (-t233 * t257 + t706 * t258 + (t333 * t472 - t334 * t470) * t464 * (t347 * t472 + t348 * t470)) * m(4) / 0.4e1) * t662, t8 * qJD(1), t3 * qJD(1) + t6 * qJD(3) + t678; t244 * qJD(2) + t7 * qJD(3) + t41 * qJD(5) + (-t628 / 0.4e1 - t624 / 0.4e1) * t664 + (((-t218 * t470 - t472 * t670) * t466 + t144) * t659 + ((-t199 * t470 - t472 * t671) * t466 + t113) * t657) * t665, t524, t7 * qJD(1) + ((-t467 * t109 + (-t138 * t472 + t139 * t470) * t466) * t657 + (-t467 * t172 + (-t202 * t472 + t203 * t470) * t466) * t659) * t663 + t101, 0, t41 * qJD(1) + t106 * t609 + t101; (t475 - t719) * qJD(1) + t49 * qJD(2) + t4 * qJD(3) + t43 * qJD(4) + t498 * qJD(5), t49 * qJD(1), t4 * qJD(1) + ((t109 * t185 - t138 * t213 + t139 * t714) * m(6) + t522) * qJD(3) + t13, t43 * qJD(1), qJD(1) * t498 + qJD(3) * t14 + t13;];
Cq = t5;
