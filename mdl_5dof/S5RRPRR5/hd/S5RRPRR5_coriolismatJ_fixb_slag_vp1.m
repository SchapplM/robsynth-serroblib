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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:33:37
% EndTime: 2019-12-05 18:33:48
% DurationCPUTime: 6.05s
% Computational Cost: add. (43773->419), mult. (27225->526), div. (0->0), fcn. (24452->10), ass. (0->281)
t424 = qJD(1) + qJD(2);
t384 = qJ(1) + qJ(2);
t380 = cos(t384);
t383 = pkin(9) + qJ(4);
t378 = qJ(5) + t383;
t371 = sin(t378);
t372 = cos(t378);
t332 = rSges(6,1) * t371 + rSges(6,2) * t372;
t376 = sin(t383);
t511 = pkin(4) * t376;
t409 = (t332 + t511) * t380;
t254 = t409 * t380;
t379 = sin(t384);
t463 = t376 * t379;
t358 = pkin(4) * t463;
t477 = t332 * t379;
t268 = t358 + t477;
t562 = m(6) / 0.2e1;
t377 = cos(t383);
t339 = rSges(5,1) * t376 + rSges(5,2) * t377;
t374 = t379 ^ 2;
t375 = t380 ^ 2;
t425 = t374 + t375;
t586 = t425 * t339;
t448 = (-t268 * t379 - t254) * t562 - m(5) * t586 / 0.2e1;
t314 = t339 * t379;
t315 = t339 * t380;
t228 = t314 * t379 + t315 * t380;
t467 = t372 * t379;
t471 = t371 * t379;
t298 = rSges(6,1) * t471 + rSges(6,2) * t467;
t263 = t358 + t298;
t563 = m(5) / 0.2e1;
t449 = (t263 * t379 + t254) * t562 + t228 * t563;
t57 = t449 - t448;
t606 = t424 * t57;
t385 = cos(pkin(9));
t373 = t385 * pkin(3) + pkin(2);
t510 = pkin(4) * t377;
t343 = t373 + t510;
t386 = -pkin(7) - qJ(3);
t382 = -pkin(8) + t386;
t427 = -rSges(6,2) * t471 - t380 * rSges(6,3);
t506 = rSges(6,1) * t372;
t224 = t380 * t382 + (t343 + t506) * t379 + t427;
t513 = sin(qJ(1)) * pkin(1);
t218 = t224 + t513;
t213 = t380 * t218;
t324 = t380 * t343;
t466 = t372 * t380;
t470 = t371 * t380;
t411 = -rSges(6,1) * t466 + rSges(6,2) * t470;
t225 = -t324 + (-rSges(6,3) + t382) * t379 + t411;
t512 = cos(qJ(1)) * pkin(1);
t219 = t225 - t512;
t507 = rSges(5,1) * t377;
t417 = t373 + t507;
t426 = -rSges(5,2) * t463 - t380 * rSges(5,3);
t590 = t380 * t386;
t236 = t417 * t379 + t426 + t590;
t233 = t236 + t513;
t227 = t380 * t233;
t363 = t379 * t386;
t462 = t376 * t380;
t416 = rSges(5,2) * t462 - t379 * rSges(5,3);
t237 = -t417 * t380 + t363 + t416;
t234 = t237 - t512;
t597 = rSges(4,2) * sin(pkin(9)) - rSges(4,1) * t385 - pkin(2);
t598 = -rSges(4,3) - qJ(3);
t255 = -t597 * t379 + t598 * t380;
t247 = t255 + t513;
t238 = t380 * t247;
t256 = t598 * t379 + t597 * t380;
t248 = t256 - t512;
t452 = t380 * t255;
t453 = t380 * t236;
t454 = t380 * t224;
t600 = m(4) / 0.2e1;
t422 = (-t238 - t452 + (-t248 - t256) * t379) * t600 + (-t213 - t454 + (-t219 - t225) * t379) * t562 + (-t227 - t453 + (-t234 - t237) * t379) * t563;
t444 = t234 - t237;
t446 = t219 - t225;
t423 = (-t238 + t452 + (-t248 + t256) * t379) * t600 + (-t446 * t379 - t213 + t454) * t562 + (-t444 * t379 - t227 + t453) * t563;
t14 = t423 - t422;
t605 = t14 * qJD(1);
t604 = t332 * t425;
t364 = Icges(6,4) * t372;
t580 = Icges(6,2) * t371 - t364;
t274 = Icges(6,6) * t380 + t580 * t379;
t603 = t274 * t470;
t273 = Icges(6,5) * t466 - Icges(6,6) * t470 + Icges(6,3) * t379;
t327 = Icges(6,5) * t372 - Icges(6,6) * t371;
t479 = t327 * t379;
t443 = t380 * (Icges(6,3) * t380 - t479) + t274 * t471;
t275 = Icges(6,4) * t466 - Icges(6,2) * t470 + Icges(6,6) * t379;
t347 = Icges(6,4) * t470;
t277 = Icges(6,1) * t466 + Icges(6,5) * t379 - t347;
t588 = (t275 * t371 - t277 * t372) * t380;
t602 = t273 * t379 + t443 - t588;
t601 = -t218 + t224;
t599 = -t379 / 0.2e1;
t542 = t379 / 0.2e1;
t540 = t380 / 0.2e1;
t538 = m(3) * (-t512 * (rSges(3,1) * t379 + rSges(3,2) * t380) - (-t380 * rSges(3,1) + t379 * rSges(3,2)) * t513);
t534 = m(4) * (-t256 * t247 + t248 * t255);
t541 = -t380 / 0.2e1;
t596 = t540 + t541;
t436 = -Icges(6,2) * t466 + t277 - t347;
t346 = Icges(6,4) * t471;
t276 = -Icges(6,1) * t467 + Icges(6,5) * t380 + t346;
t437 = Icges(6,2) * t467 + t276 + t346;
t583 = Icges(6,1) * t371 + t364;
t438 = t380 * t583 + t275;
t439 = -t379 * t583 + t274;
t595 = -(t436 * t379 + t437 * t380) * t371 - (t438 * t379 + t439 * t380) * t372;
t355 = Icges(5,4) * t462;
t458 = t377 * t380;
t288 = Icges(5,1) * t458 + Icges(5,5) * t379 - t355;
t432 = -Icges(5,2) * t458 + t288 - t355;
t354 = Icges(5,4) * t463;
t459 = t377 * t379;
t287 = -Icges(5,1) * t459 + Icges(5,5) * t380 + t354;
t433 = Icges(5,2) * t459 + t287 + t354;
t286 = Icges(5,4) * t458 - Icges(5,2) * t462 + Icges(5,6) * t379;
t365 = Icges(5,4) * t377;
t582 = Icges(5,1) * t376 + t365;
t434 = t380 * t582 + t286;
t579 = Icges(5,2) * t376 - t365;
t285 = Icges(5,6) * t380 + t579 * t379;
t435 = -t379 * t582 + t285;
t594 = -(t379 * t432 + t380 * t433) * t376 - (t379 * t434 + t380 * t435) * t377;
t299 = t332 * t380;
t112 = -t224 * t298 + t225 * t299;
t593 = m(6) * t112;
t262 = t380 * (t379 * rSges(6,3) - t411);
t278 = -rSges(6,1) * t467 - t427;
t110 = -t380 * (t380 * t373 - t324 - t363) + t262 + (-t278 + (t343 - t373) * t379 - t590) * t379;
t333 = -rSges(6,2) * t371 + t506;
t476 = t333 * t379;
t404 = Icges(6,5) * t371 + Icges(6,6) * t372;
t292 = t379 * t404;
t293 = t404 * t380;
t509 = (t375 * t292 + (-t380 * t293 - t595) * t379) * t540 + (-t374 * t293 + (t379 * t292 + t595) * t380) * t542;
t584 = t379 * t298 + t380 * t299;
t18 = t509 + m(6) * (-t110 * t584 + t333 * t254 + t268 * t476);
t592 = t18 * qJD(5);
t340 = -rSges(5,2) * t376 + t507;
t591 = t340 * t563;
t589 = (t286 * t376 - t288 * t377) * t380;
t163 = t409 * t219;
t168 = t409 * t225;
t397 = m(6) * t584;
t399 = -m(6) * t604 / 0.2e1;
t151 = -t397 / 0.2e1 + t399;
t587 = t424 * t151;
t91 = -t225 * t218 + t219 * t224;
t101 = -t237 * t233 + t234 * t236;
t585 = t234 * t379 + t227;
t581 = t237 * t379 + t453;
t487 = t236 * t314;
t489 = t233 * t314;
t508 = (t601 * t268 + t163 - t168) * t562 + (t444 * t315 + t487 - t489) * t563;
t102 = -t218 * t263 + t163;
t104 = -t224 * t263 + t168;
t117 = t234 * t315 - t489;
t120 = t237 * t315 - t487;
t577 = (t104 + t102) * t562 + (t120 + t117) * t563;
t495 = Icges(6,4) * t371;
t328 = Icges(6,2) * t372 + t495;
t331 = Icges(6,1) * t372 - t495;
t576 = (-t580 + t583) * t371 + (t328 - t331) * t372;
t496 = Icges(5,4) * t376;
t335 = Icges(5,2) * t377 + t496;
t338 = Icges(5,1) * t377 - t496;
t575 = (-t579 + t582) * t376 + (t335 - t338) * t377;
t570 = (t582 / 0.2e1 - t579 / 0.2e1) * t377 + (t338 / 0.2e1 - t335 / 0.2e1) * t376;
t413 = (-t580 / 0.2e1 + t583 / 0.2e1) * t372 + (-t328 / 0.2e1 + t331 / 0.2e1) * t371;
t132 = -t276 * t467 + t443;
t133 = t380 * t273 + t275 * t471 - t277 * t467;
t415 = -t276 * t372 - t273;
t482 = t274 * t371;
t11 = ((t588 + t602) * t380 + ((t415 - t482) * t380 + t133 + t603) * t379) * t542 + (t132 * t380 + t133 * t379) * t599 + ((t133 + (-t273 + t482) * t380 - t603) * t380 + (t415 * t379 - t132 + t602) * t379) * t540;
t569 = 0.4e1 * qJD(1);
t567 = 4 * qJD(2);
t566 = 2 * qJD(4);
t111 = -t218 * t298 + t219 * t299;
t552 = m(6) * (t112 + t111);
t551 = m(6) * (t446 * t299 + t601 * t477);
t182 = t219 * t476;
t445 = t268 * t299 - t298 * t409;
t550 = m(6) * (t333 * t213 + t182 + t445);
t226 = t409 * t477;
t483 = t263 * t332;
t492 = t218 * t333;
t549 = m(6) * (t182 + t226 + (-t483 + t492) * t380);
t190 = t225 * t476;
t548 = m(6) * (t333 * t454 + t190 + t445);
t491 = t224 * t333;
t547 = m(6) * (t190 + t226 + (-t483 + t491) * t380);
t544 = m(6) * t91;
t284 = Icges(5,5) * t458 - Icges(5,6) * t462 + Icges(5,3) * t379;
t148 = t380 * t284 + t286 * t463 - t288 * t459;
t334 = Icges(5,5) * t377 - Icges(5,6) * t376;
t475 = t334 * t379;
t283 = Icges(5,3) * t380 - t475;
t440 = t379 * t283 + t287 * t458;
t149 = -t285 * t462 + t440;
t150 = t284 * t379 - t589;
t414 = -t287 * t377 - t284;
t441 = t380 * t283 + t285 * t463;
t481 = t285 * t376;
t31 = (t150 + t441 + t589) * t380 + (-t149 + (t414 - t481) * t380 + t148 + t440) * t379;
t147 = -t287 * t459 + t441;
t32 = (t148 + (-t284 + t481) * t380 - t440) * t380 + (t414 * t379 - t147 + t441) * t379;
t92 = t147 * t380 + t148 * t379;
t93 = t149 * t380 + t150 * t379;
t2 = (t32 / 0.2e1 + t93 / 0.2e1) * t380 + (-t92 / 0.2e1 + t31 / 0.2e1) * t379 + t11;
t539 = -qJD(3) * t57 + t2 * qJD(4);
t532 = m(4) * (-t248 * t379 - t238);
t531 = m(4) * (-t256 * t379 - t452);
t530 = m(5) * t101;
t528 = m(5) * t117;
t527 = m(5) * t120;
t526 = m(5) * t585;
t525 = m(5) * t581;
t522 = m(6) * t102;
t521 = m(6) * t104;
t520 = m(6) * t111;
t518 = m(6) * (-t219 * t379 - t213);
t517 = m(6) * (-t225 * t379 - t454);
t501 = t151 * qJD(3) + t11 * qJD(5);
t152 = t397 / 0.2e1 + t399;
t58 = t448 + t449;
t500 = t58 * qJD(4) + t152 * qJD(5);
t499 = t57 * qJD(4) - t151 * qJD(5);
t480 = t298 * t332;
t447 = (-t263 + t268) * t409;
t418 = t333 + t510;
t408 = t552 / 0.2e1 + t413;
t405 = Icges(5,5) * t376 + Icges(5,6) * t377;
t308 = t379 * t405;
t396 = -t11 + (-t438 * t371 + t436 * t372 - t576 * t380 + t479) * t542 + (t327 * t380 - t439 * t371 + t437 * t372 + t576 * t379) * t540;
t395 = -t413 + t596 * (t275 * t372 + t277 * t371);
t394 = t413 + t570;
t391 = t394 + t577;
t390 = t395 - t570 + t596 * (t286 * t377 + t288 * t376);
t389 = t58 * qJD(3) + (t31 * t599 + t396 + (-t434 * t376 + t432 * t377 - t575 * t380 + t475 + t92) * t542 + (t93 + t32) * t541 + (t334 * t380 - t435 * t376 + t433 * t377 + t575 * t379) * t540) * qJD(4);
t309 = t405 * t380;
t271 = t418 * t380;
t269 = t418 * t379;
t239 = t299 * t477;
t203 = -t379 * t278 + t262;
t183 = -t425 * t511 - t584;
t144 = t152 * qJD(3);
t78 = t413 + t593;
t75 = t413 + t520;
t73 = t547 / 0.2e1;
t72 = t548 / 0.2e1;
t68 = t549 / 0.2e1;
t67 = t550 / 0.2e1;
t63 = t551 / 0.2e1;
t62 = t517 - t525 + t531;
t51 = t518 - t526 + t532;
t40 = t394 + t521 + t527;
t39 = t394 + t522 + t528;
t33 = t530 + t534 + t538 + t544;
t22 = -t551 / 0.2e1 + t408;
t21 = t63 + t408;
t20 = m(6) * (-t203 * t584 + t333 * t604) + t509;
t19 = t20 * qJD(5);
t17 = t63 - t552 / 0.2e1 + t395;
t16 = t422 + t423;
t13 = t391 + t508;
t12 = t391 - t508;
t9 = t390 + t508 - t577;
t8 = t72 - t547 / 0.2e1 + t11;
t7 = t73 - t548 / 0.2e1 + t11;
t6 = t67 - t549 / 0.2e1 + t11;
t5 = t68 - t550 / 0.2e1 + t11;
t4 = t72 + t73 + t396;
t3 = t67 + t68 + t396;
t1 = [qJD(2) * t33 + qJD(3) * t51 + qJD(4) * t39 + qJD(5) * t75, t33 * qJD(1) + t16 * qJD(3) + t13 * qJD(4) + t21 * qJD(5) + 0.2e1 * (t534 / 0.2e1 + t101 * t563 + t538 / 0.2e1 + t91 * t562) * qJD(2), qJD(1) * t51 + qJD(2) * t16 + t500, t39 * qJD(1) + t13 * qJD(2) + t3 * qJD(5) + (t585 * t591 + (t218 * t271 + t219 * t269 + t447) * t562) * t566 + t389, t75 * qJD(1) + t21 * qJD(2) + t144 + t3 * qJD(4) + ((t182 + t239 + (-t480 + t492) * t380) * m(6) + t396) * qJD(5); -t14 * qJD(3) + t12 * qJD(4) + t22 * qJD(5) + (-t538 / 0.4e1 - t534 / 0.4e1 - t530 / 0.4e1 - t544 / 0.4e1) * t569, qJD(3) * t62 + qJD(4) * t40 + qJD(5) * t78, qJD(2) * t62 + t500 - t605, t12 * qJD(1) + t40 * qJD(2) + t4 * qJD(5) + ((t224 * t271 + t225 * t269 + t447) * t562 + t581 * t591) * t566 + t389, t22 * qJD(1) + t78 * qJD(2) + t144 + t4 * qJD(4) + ((t190 + t239 + (-t480 + t491) * t380) * m(6) + t396) * qJD(5); t14 * qJD(2) + (-t532 / 0.4e1 + t526 / 0.4e1 - t518 / 0.4e1) * t569 + t499, t605 + (-t517 / 0.4e1 + t525 / 0.4e1 - t531 / 0.4e1) * t567 + t499, 0, m(6) * (t269 * t380 - t271 * t379) * qJD(4) + t606, -t587; t9 * qJD(2) + t6 * qJD(5) + (-t522 / 0.4e1 - t528 / 0.4e1) * t569 + t539 + t390 * qJD(1), t9 * qJD(1) + t390 * qJD(2) + t8 * qJD(5) + (-t521 / 0.4e1 - t527 / 0.4e1) * t567 + t539, -t606, (m(5) * (t340 * t586 - (t380 * (rSges(5,1) * t458 - t416) - t379 * (-rSges(5,1) * t459 - t426)) * t228) + (t375 * t308 + (-t380 * t309 - t594) * t379) * t540 + (-t374 * t309 + (t379 * t308 + t594) * t380) * t542 + m(6) * (t110 * t183 + t268 * t269 + t271 * t409) + t509) * qJD(4) + t592 + t424 * t2, t6 * qJD(1) + t8 * qJD(2) + t18 * qJD(4) + t592; (t395 - t520) * qJD(1) + t17 * qJD(2) + t5 * qJD(4) + t501, t17 * qJD(1) + (t395 - t593) * qJD(2) + t7 * qJD(4) + t501, t587, t5 * qJD(1) + t7 * qJD(2) + ((t183 * t203 + (t269 * t379 + t271 * t380) * t332) * m(6) + t509) * qJD(4) + t19, qJD(4) * t20 + t11 * t424 + t19;];
Cq = t1;
