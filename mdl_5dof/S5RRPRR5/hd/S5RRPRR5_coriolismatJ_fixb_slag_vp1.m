% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:12
% EndTime: 2022-01-20 11:02:32
% DurationCPUTime: 8.02s
% Computational Cost: add. (45073->394), mult. (27917->517), div. (0->0), fcn. (25096->10), ass. (0->273)
t393 = qJ(1) + qJ(2);
t388 = sin(t393);
t389 = cos(t393);
t394 = cos(pkin(9));
t610 = rSges(4,2) * sin(pkin(9)) - rSges(4,1) * t394 - pkin(2);
t613 = rSges(4,3) + qJ(3);
t257 = t610 * t388 + t613 * t389;
t526 = sin(qJ(1)) * pkin(1);
t250 = t257 - t526;
t258 = t613 * t388 - t610 * t389;
t529 = cos(qJ(1)) * pkin(1);
t251 = t258 + t529;
t164 = t250 * t389 + t251 * t388;
t174 = t257 * t389 + t258 * t388;
t578 = m(6) / 0.2e1;
t579 = m(5) / 0.2e1;
t382 = t394 * pkin(3) + pkin(2);
t392 = pkin(9) + qJ(4);
t386 = cos(t392);
t521 = t386 * rSges(5,1);
t440 = t382 + t521;
t525 = pkin(7) + qJ(3);
t443 = t389 * t525;
t385 = sin(t392);
t483 = t385 * t388;
t453 = rSges(5,2) * t483 + t389 * rSges(5,3);
t237 = -t388 * t440 + t443 + t453;
t365 = t388 * t525;
t482 = t385 * t389;
t439 = -rSges(5,2) * t482 + t388 * rSges(5,3);
t238 = t389 * t440 + t365 + t439;
t591 = -t237 * t389 - t238 * t388;
t234 = t237 - t526;
t235 = t238 + t529;
t595 = -t234 * t389 - t235 * t388;
t387 = qJ(5) + t392;
t381 = cos(t387);
t486 = t381 * t389;
t380 = sin(t387);
t490 = t380 * t389;
t446 = rSges(6,1) * t486 - rSges(6,2) * t490 + t388 * rSges(6,3);
t527 = pkin(4) * t386;
t346 = t382 + t527;
t391 = pkin(8) + t525;
t621 = t389 * t346 + t391 * t388;
t225 = t446 + t621;
t217 = t529 + t225;
t364 = t391 * t389;
t491 = t380 * t388;
t454 = rSges(6,2) * t491 + t389 * rSges(6,3);
t522 = rSges(6,1) * t381;
t593 = t364 + (-t346 - t522) * t388 + t454;
t607 = t593 - t526;
t597 = -t217 * t388 - t389 * t607;
t617 = m(4) / 0.2e1;
t618 = -t225 * t388 - t389 * t593;
t448 = (-t591 - t595) * t579 + (t174 + t164) * t617 + (-t618 - t597) * t578;
t449 = (-t595 + t591) * t579 + (t164 - t174) * t617 + (t618 - t597) * t578;
t14 = t449 - t448;
t626 = t14 * qJD(1);
t335 = rSges(6,1) * t380 + rSges(6,2) * t381;
t301 = t335 * t388;
t181 = t607 * t301;
t302 = t335 * t389;
t115 = -t217 * t302 + t181;
t625 = m(6) * t115;
t116 = -t225 * t302 + t301 * t593;
t624 = m(6) * t116;
t450 = qJD(1) + qJD(2);
t528 = pkin(4) * t385;
t409 = t335 + t528;
t601 = t409 * t389;
t602 = t409 * t388;
t592 = t388 * t602 + t389 * t601;
t342 = rSges(5,1) * t385 + rSges(5,2) * t386;
t383 = t388 ^ 2;
t384 = t389 ^ 2;
t451 = t383 + t384;
t599 = t451 * t342;
t469 = -m(6) * t592 / 0.2e1 - m(5) * t599 / 0.2e1;
t317 = t342 * t388;
t318 = t342 * t389;
t229 = t317 * t388 + t318 * t389;
t470 = t229 * t579 + t592 * t578;
t56 = t470 - t469;
t623 = t450 * t56;
t622 = t335 * t451;
t366 = Icges(6,4) * t381;
t332 = -Icges(6,2) * t380 + t366;
t333 = Icges(6,1) * t380 + t366;
t620 = t332 + t333;
t367 = Icges(5,4) * t386;
t339 = -Icges(5,2) * t385 + t367;
t340 = Icges(5,1) * t385 + t367;
t619 = t339 + t340;
t471 = t225 * t601 - t593 * t602;
t96 = t217 * t593 - t225 * t607;
t606 = -t217 * t601 + t602 * t607;
t558 = t388 / 0.2e1;
t557 = -t389 / 0.2e1;
t614 = t389 / 0.2e1;
t555 = m(3) * (t529 * (-rSges(3,1) * t388 - rSges(3,2) * t389) + (t389 * rSges(3,1) - t388 * rSges(3,2)) * t526);
t551 = m(4) * (-t258 * t250 + t251 * t257);
t559 = -t388 / 0.2e1;
t608 = t558 + t559;
t598 = (t382 - t346) * t388;
t241 = (-t391 + t525) * t389 - t598;
t605 = m(6) * (t241 + t364 - t443 + t598) * t389;
t487 = t381 * t388;
t199 = t388 * (rSges(6,1) * t487 - t454) + t389 * t446;
t114 = t241 * t388 + (-t389 * t382 - t365 + t621) * t389 + t199;
t336 = -rSges(6,2) * t380 + t522;
t420 = Icges(6,5) * t380 + Icges(6,6) * t381;
t295 = t420 * t388;
t296 = t389 * t420;
t509 = Icges(6,4) * t380;
t334 = Icges(6,1) * t381 - t509;
t280 = Icges(6,5) * t388 + t334 * t389;
t331 = Icges(6,2) * t381 + t509;
t459 = -t331 * t389 + t280;
t349 = Icges(6,4) * t491;
t279 = Icges(6,1) * t487 - Icges(6,5) * t389 - t349;
t460 = -Icges(6,2) * t487 + t279 - t349;
t278 = Icges(6,6) * t388 + t332 * t389;
t461 = -t333 * t389 - t278;
t277 = Icges(6,4) * t487 - Icges(6,2) * t491 - Icges(6,6) * t389;
t462 = t333 * t388 + t277;
t586 = (-t459 * t388 + t460 * t389) * t380 + (t461 * t388 + t462 * t389) * t381;
t524 = (-t383 * t296 + (t388 * t295 + t586) * t389) * t558 + (-t384 * t295 + (t389 * t296 + t586) * t388) * t557;
t594 = t388 * t301 + t389 * t302;
t18 = t524 + m(6) * (-t114 * t594 + t336 * t592);
t603 = t18 * qJD(5);
t415 = t594 * t578;
t427 = m(6) * t622;
t152 = t415 + t427 / 0.2e1;
t600 = t450 * t152;
t106 = -t238 * t234 + t235 * t237;
t523 = (t471 + t606) * t578 + ((-t235 + t238) * t389 + (t234 - t237) * t388) * t342 * t579;
t119 = t234 * t317 - t235 * t318;
t122 = t237 * t317 - t238 * t318;
t589 = (-t471 + t606) * t578 + (t122 + t119) * t579;
t510 = Icges(5,4) * t385;
t338 = Icges(5,2) * t386 + t510;
t341 = Icges(5,1) * t386 - t510;
t588 = t619 * t386 / 0.2e1 + (t341 / 0.2e1 - t338 / 0.2e1) * t385;
t429 = t620 * t381 / 0.2e1 + (-t331 / 0.2e1 + t334 / 0.2e1) * t380;
t244 = t280 * t487;
t330 = Icges(6,5) * t381 - Icges(6,6) * t380;
t497 = t330 * t389;
t276 = Icges(6,3) * t388 + t497;
t433 = t389 * t276 - t244;
t133 = -t278 * t491 - t433;
t275 = Icges(6,5) * t487 - Icges(6,6) * t491 - Icges(6,3) * t389;
t466 = -t388 * t275 - t279 * t486;
t134 = -t277 * t490 - t466;
t465 = t388 * t276 + t280 * t486;
t135 = -t278 * t490 + t465;
t431 = t278 * t380 - t275;
t500 = t277 * t380;
t441 = ((t133 - t244 + (t276 + t500) * t389 + t466) * t389 + t465 * t388) * t557 + (-t134 * t389 + t135 * t388) * t614 + ((t388 * t431 + t133 + t134 + t433) * t388 + (t135 - t465 + t388 * (-t279 * t381 + t500) + (t431 + t275) * t389) * t389) * t558;
t291 = Icges(5,5) * t388 + t341 * t389;
t455 = -t338 * t389 + t291;
t355 = Icges(5,4) * t483;
t479 = t386 * t388;
t290 = Icges(5,1) * t479 - Icges(5,5) * t389 - t355;
t456 = -Icges(5,2) * t479 + t290 - t355;
t289 = Icges(5,6) * t388 + t339 * t389;
t457 = -t340 * t389 - t289;
t288 = Icges(5,4) * t479 - Icges(5,2) * t483 - Icges(5,6) * t389;
t458 = t340 * t388 + t288;
t587 = (-t455 * t388 + t389 * t456) * t385 + (t457 * t388 + t389 * t458) * t386;
t585 = 0.4e1 * qJD(1);
t583 = 0.4e1 * qJD(2);
t582 = 2 * qJD(4);
t575 = t199 * t605;
t574 = t114 * t605;
t569 = m(6) * (t116 + t115);
t568 = m(6) * (t181 + (-t593 * t388 + (-t217 + t225) * t389) * t335);
t413 = t597 * t336;
t467 = -t301 * t601 + t302 * t602;
t566 = m(6) * (t413 + t467);
t411 = (t388 * t601 - t389 * t602) * t335;
t565 = m(6) * (t413 + t411);
t412 = t618 * t336;
t564 = m(6) * (t412 + t467);
t563 = m(6) * (t412 + t411);
t560 = m(6) * t96;
t286 = Icges(5,5) * t479 - Icges(5,6) * t483 - Icges(5,3) * t389;
t478 = t386 * t389;
t464 = -t388 * t286 - t290 * t478;
t149 = -t288 * t482 - t464;
t337 = Icges(5,5) * t386 - Icges(5,6) * t385;
t495 = t337 * t389;
t287 = Icges(5,3) * t388 + t495;
t463 = t388 * t287 + t291 * t478;
t150 = -t289 * t482 + t463;
t430 = t289 * t385 - t286;
t252 = t291 * t479;
t432 = t389 * t287 - t252;
t34 = (t389 * t430 + t150 - t463) * t389 + (t388 * t430 + t149 + t432) * t388;
t148 = -t289 * t483 - t432;
t498 = t288 * t385;
t35 = (t148 - t252 + (t287 + t498) * t389 + t464) * t389 + t463 * t388;
t97 = -(-t388 * (-t290 * t386 + t498) - t286 * t389) * t389 + t148 * t388;
t98 = -t149 * t389 + t150 * t388;
t2 = t574 + (t98 / 0.2e1 - t35 / 0.2e1) * t389 + (t34 / 0.2e1 + t97 / 0.2e1) * t388 + t441;
t556 = -qJD(3) * t56 + t2 * qJD(4);
t549 = m(4) * t164;
t548 = m(4) * t174;
t545 = m(5) * t106;
t543 = m(5) * t119;
t542 = m(5) * t122;
t541 = m(5) * t595;
t540 = m(5) * t591;
t536 = m(6) * t606;
t535 = m(6) * t471;
t534 = m(6) * t597;
t533 = m(6) * t618;
t516 = -t152 * qJD(3) + qJD(5) * t441;
t151 = t415 - t427 / 0.2e1;
t57 = t469 + t470;
t514 = t57 * qJD(4) + t151 * qJD(5);
t513 = t56 * qJD(4) + t152 * qJD(5);
t442 = -t336 - t527;
t426 = t575 / 0.2e1 + t441;
t424 = t569 / 0.2e1 + t429;
t421 = Icges(5,5) * t385 + Icges(5,6) * t386;
t410 = (-t317 * t389 + t318 * t388) * t342;
t404 = (-t331 + t334) * t381 - t620 * t380;
t408 = -t441 + (t330 * t388 + t461 * t380 + t459 * t381 + t389 * t404) * t558 + (-t462 * t380 + t460 * t381 + t388 * t404 - t497) * t557;
t406 = -t429 + t608 * (t277 * t381 + t279 * t380);
t405 = t429 + t588;
t403 = (-t338 + t341) * t386 - t619 * t385;
t402 = -t575 / 0.2e1 + t408;
t400 = t405 + t589;
t399 = (-t301 * t389 + t302 * t388) * t335;
t398 = t406 - t588 + t608 * (t288 * t386 + t290 * t385);
t397 = t57 * qJD(3) + (t35 * t614 - t574 + t408 + (-t458 * t385 + t456 * t386 + t388 * t403 - t495 + t98) * t557 + (t34 + t97) * t559 + (t337 * t388 + t457 * t385 + t455 * t386 + t389 * t403) * t558) * qJD(4);
t343 = -rSges(5,2) * t385 + t521;
t312 = t389 * t421;
t311 = t421 * t388;
t274 = t442 * t389;
t272 = t442 * t388;
t182 = -t451 * t528 - t594;
t143 = t151 * qJD(3);
t81 = t429 + t624;
t78 = t429 + t625;
t76 = t563 / 0.2e1;
t74 = t564 / 0.2e1;
t70 = t565 / 0.2e1;
t69 = t566 / 0.2e1;
t63 = t568 / 0.2e1;
t62 = -t533 - t540 + t548;
t51 = -t534 - t541 + t549;
t42 = t405 - t535 + t542;
t41 = t405 + t536 + t543;
t36 = t545 + t551 + t555 + t560;
t22 = m(6) * (-t199 * t594 + t336 * t622) + t524;
t21 = t22 * qJD(5);
t20 = -t568 / 0.2e1 + t424;
t19 = t63 + t424;
t17 = t63 - t569 / 0.2e1 + t406;
t15 = t448 + t449;
t13 = t400 - t523;
t12 = t400 + t523;
t9 = t398 + t523 - t589;
t8 = t74 - t563 / 0.2e1 + t426;
t7 = t76 - t564 / 0.2e1 + t426;
t6 = t69 - t565 / 0.2e1 + t426;
t5 = t70 - t566 / 0.2e1 + t426;
t4 = t74 + t76 + t402;
t3 = t69 + t70 + t402;
t1 = [t36 * qJD(2) + t51 * qJD(3) + t41 * qJD(4) + t78 * qJD(5), t36 * qJD(1) + t15 * qJD(3) + t12 * qJD(4) + t19 * qJD(5) + 0.2e1 * (t96 * t578 + t106 * t579 + t551 / 0.2e1 + t555 / 0.2e1) * qJD(2), qJD(1) * t51 + qJD(2) * t15 + t514, t41 * qJD(1) + t12 * qJD(2) + t3 * qJD(5) + ((t217 * t272 + t274 * t607) * t578 + (t343 * t595 + t410) * t579) * t582 + t397, t78 * qJD(1) + t19 * qJD(2) + t143 + t3 * qJD(4) + ((t413 + t399) * m(6) + t408) * qJD(5); -t14 * qJD(3) + t13 * qJD(4) + t20 * qJD(5) + (-t560 / 0.4e1 - t545 / 0.4e1 - t551 / 0.4e1 - t555 / 0.4e1) * t585, qJD(3) * t62 + qJD(4) * t42 + qJD(5) * t81, qJD(2) * t62 + t514 - t626, t13 * qJD(1) + t42 * qJD(2) + t4 * qJD(5) + ((t225 * t272 + t274 * t593) * t578 + (t343 * t591 + t410) * t579) * t582 + t397, t20 * qJD(1) + t81 * qJD(2) + t143 + t4 * qJD(4) + ((t412 + t399) * m(6) + t408) * qJD(5); t14 * qJD(2) + (t534 / 0.4e1 + t541 / 0.4e1 - t549 / 0.4e1) * t585 + t513, t626 + (t533 / 0.4e1 + t540 / 0.4e1 - t548 / 0.4e1) * t583 + t513, 0, m(6) * (-t272 * t389 + t274 * t388) * qJD(4) + t623, t600; t398 * qJD(1) + t9 * qJD(2) + t6 * qJD(5) + (-t536 / 0.4e1 - t543 / 0.4e1) * t585 + t556, t9 * qJD(1) + t8 * qJD(5) + (t535 / 0.4e1 - t542 / 0.4e1) * t583 + t556 + t398 * qJD(2), -t623, (m(5) * (t343 * t599 - (t388 * (rSges(5,1) * t479 - t453) + t389 * (rSges(5,1) * t478 + t439)) * t229) + (-t383 * t312 + (t388 * t311 + t587) * t389) * t558 + (-t384 * t311 + (t389 * t312 + t587) * t388) * t557 + m(6) * (t114 * t182 - t272 * t602 - t274 * t601) + t524) * qJD(4) + t603 + t450 * t2, t6 * qJD(1) + t8 * qJD(2) + t18 * qJD(4) + t603; (t406 - t625) * qJD(1) + t17 * qJD(2) + t5 * qJD(4) + t516, t17 * qJD(1) + (t406 - t624) * qJD(2) + t7 * qJD(4) + t516, -t600, t5 * qJD(1) + t7 * qJD(2) + ((t182 * t199 + (-t272 * t388 - t274 * t389) * t335) * m(6) + t524) * qJD(4) + t21, qJD(4) * t22 + t441 * t450 + t21;];
Cq = t1;
