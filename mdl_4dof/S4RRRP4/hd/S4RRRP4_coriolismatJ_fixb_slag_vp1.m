% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:23
% DurationCPUTime: 8.46s
% Computational Cost: add. (16563->348), mult. (20560->476), div. (0->0), fcn. (19182->6), ass. (0->231)
t578 = Icges(4,4) + Icges(5,4);
t576 = Icges(4,1) + Icges(5,1);
t575 = Icges(4,5) + Icges(5,5);
t570 = Icges(4,6) + Icges(5,6);
t340 = qJ(2) + qJ(3);
t321 = sin(t340);
t579 = t578 * t321;
t574 = Icges(4,2) + Icges(5,2);
t322 = cos(t340);
t577 = -t570 * t321 + t575 * t322;
t556 = t576 * t322 - t579;
t573 = Icges(4,3) + Icges(5,3);
t572 = t578 * t322;
t342 = sin(qJ(1));
t344 = cos(qJ(1));
t565 = t575 * t342 + t556 * t344;
t569 = -t574 * t321 + t572;
t441 = t321 * t342;
t568 = t578 * t441;
t567 = t577 * t344;
t435 = t322 * t342;
t558 = -t573 * t344 + t575 * t435 - t570 * t441;
t566 = t573 * t342 + t567;
t544 = -t575 * t344 + t576 * t435 - t568;
t557 = t574 * t322 + t579;
t564 = t576 * t321 + t572;
t563 = t565 * t435;
t562 = -t570 * t344 + t578 * t435 - t441 * t574;
t561 = t570 * t342 + t569 * t344;
t338 = t342 ^ 2;
t339 = t344 ^ 2;
t394 = t338 + t339;
t472 = rSges(5,2) * t322;
t283 = rSges(5,1) * t321 + t472;
t396 = rSges(5,1) * t441 + rSges(5,2) * t435;
t403 = -t339 * t283 - t342 * t396;
t482 = pkin(3) * t321;
t126 = -t394 * t482 + t403;
t341 = sin(qJ(2));
t483 = pkin(2) * t341;
t288 = -t482 - t483;
t356 = t283 - t288;
t155 = t356 * t342;
t157 = t356 * t344;
t473 = rSges(5,1) * t322;
t481 = pkin(3) * t322;
t388 = rSges(5,2) * t321 - t473 - t481;
t191 = t388 * t342;
t193 = t388 * t344;
t421 = -t155 * t191 - t157 * t193;
t345 = -pkin(6) - pkin(5);
t320 = t342 * t345;
t337 = t344 * pkin(5);
t423 = t344 * t345;
t343 = cos(qJ(2));
t336 = t343 * pkin(2);
t319 = t336 + pkin(1);
t480 = pkin(1) - t319;
t414 = -t342 * (t480 * t342 - t337 - t423) + t344 * (-pkin(5) * t342 - t480 * t344 - t320);
t287 = t319 + t481;
t559 = rSges(5,3) + qJ(4) - t345;
t144 = -rSges(5,1) * t435 + rSges(5,2) * t441 - t342 * t287 + t344 * t559;
t434 = t322 * t344;
t440 = t321 * t344;
t529 = -rSges(5,2) * t440 + t342 * t559;
t75 = (rSges(5,1) * t434 + t320 + (t287 - t319) * t344 + t529) * t344 + (-t342 * t319 - t144 - t423) * t342;
t62 = t75 + t414;
t30 = t62 * t126 + t421;
t284 = rSges(4,1) * t321 + rSges(4,2) * t322;
t261 = t284 * t342;
t262 = t284 * t344;
t146 = -t342 * t261 - t344 * t262;
t474 = rSges(4,1) * t322;
t286 = -rSges(4,2) * t321 + t474;
t352 = t284 + t483;
t521 = t352 * t344;
t522 = t352 * t342;
t385 = -rSges(4,2) * t440 + rSges(4,3) * t342;
t397 = rSges(4,2) * t441 + t344 * rSges(4,3);
t136 = t342 * (rSges(4,1) * t435 - t397) + t344 * (rSges(4,1) * t434 + t385);
t88 = t136 + t414;
t42 = t88 * t146 + (t342 * t522 + t344 * t521) * t286;
t560 = -m(4) * t42 - m(5) * t30;
t555 = -t575 * t321 - t570 * t322;
t554 = -t566 * t344 + t563;
t553 = -t557 * t344 + t565;
t552 = -t435 * t574 + t544 - t568;
t551 = -t564 * t344 - t561;
t550 = t564 * t342 + t562;
t549 = -t558 * t342 - t544 * t434;
t540 = t566 * t342 + t565 * t434;
t548 = -t564 - t569;
t547 = -t561 * t441 + t554;
t546 = t562 * t440 + t549;
t545 = -t561 * t440 + t540;
t541 = t561 * t321 - t558;
t539 = t562 * t321;
t511 = m(4) / 0.2e1;
t510 = m(5) / 0.2e1;
t145 = (t287 + t473) * t344 + t529;
t538 = m(5) * (t144 * t344 + t145 * t342);
t493 = t342 / 0.2e1;
t491 = -t344 / 0.2e1;
t534 = t344 / 0.2e1;
t460 = Icges(3,4) * t341;
t290 = Icges(3,2) * t343 + t460;
t293 = Icges(3,1) * t343 - t460;
t533 = (t293 / 0.2e1 - t290 / 0.2e1) * t341;
t532 = t555 * t342;
t531 = t555 * t344;
t530 = (t556 - t557) * t322 + t548 * t321;
t159 = -t288 * t342 + t396;
t272 = t344 * t288;
t160 = -t283 * t344 + t272;
t389 = t283 + t482;
t190 = t389 * t342;
t192 = t389 * t344;
t386 = t319 + t474;
t151 = -t342 * t386 + t397 - t423;
t152 = t344 * t386 - t320 + t385;
t354 = (-t151 * t344 - t152 * t342) * t286;
t422 = t193 * t144 + t191 * t145;
t478 = (-t159 * t192 - t160 * t190 + t422) * t510 + (t354 + (t342 * t521 - t344 * t522) * t284) * t511;
t172 = pkin(3) * t441 + t396;
t173 = (-t472 + (-rSges(5,1) - pkin(3)) * t321) * t344;
t479 = (-t155 * t173 - t157 * t172 + t422) * t510 + (-t261 * t521 + t262 * t522 + t354) * t511;
t528 = t478 - t479;
t526 = (t342 * t551 + t344 * t550) * t322 + (-t342 * t553 + t344 * t552) * t321;
t523 = qJD(1) * t538;
t332 = Icges(3,4) * t343;
t291 = -Icges(3,2) * t341 + t332;
t292 = Icges(3,1) * t341 + t332;
t247 = Icges(3,5) * t342 + t293 * t344;
t399 = -t290 * t344 + t247;
t433 = t341 * t342;
t312 = Icges(3,4) * t433;
t429 = t342 * t343;
t246 = Icges(3,1) * t429 - Icges(3,5) * t344 - t312;
t400 = -Icges(3,2) * t429 + t246 - t312;
t245 = Icges(3,6) * t342 + t291 * t344;
t401 = -t292 * t344 - t245;
t244 = Icges(3,4) * t429 - Icges(3,2) * t433 - Icges(3,6) * t344;
t402 = t292 * t342 + t244;
t518 = (-t399 * t342 + t344 * t400) * t341 + (t401 * t342 + t344 * t402) * t343;
t351 = -t548 * t322 / 0.2e1 + (-t557 / 0.2e1 + t556 / 0.2e1) * t321;
t353 = (t545 * t342 + t546 * t344) * t534 + (t540 * t342 + ((t539 + t566) * t344 + t547 + t549 - t563) * t344) * t491 + ((t342 * t541 - t546 + t547 - t554) * t342 + ((t541 + t558) * t344 + (-t544 * t322 + t539) * t342 - t540 + t545) * t344) * t493;
t515 = 0.4e1 * qJD(1);
t514 = 2 * qJD(2);
t512 = 2 * qJD(3);
t69 = t394 * t284 * t286 + t136 * t146;
t68 = m(4) * t69;
t503 = m(4) * (t151 * t522 - t152 * t521);
t502 = m(4) * (t151 * t261 - t152 * t262);
t420 = -t190 * t191 - t192 * t193;
t432 = t341 * t344;
t92 = t338 * (t288 + t483) + t344 * (pkin(2) * t432 + t272) + t403;
t499 = m(5) * (t75 * t92 + t420);
t496 = m(5) * (t144 * t159 + t145 * t160);
t495 = m(5) * (t144 * t172 + t145 * t173);
t494 = -t342 / 0.2e1;
t475 = rSges(3,1) * t343;
t391 = pkin(1) + t475;
t395 = rSges(3,2) * t433 + t344 * rSges(3,3);
t170 = -t342 * t391 + t337 + t395;
t314 = rSges(3,2) * t432;
t171 = -t314 + t391 * t344 + (rSges(3,3) + pkin(5)) * t342;
t294 = rSges(3,1) * t341 + rSges(3,2) * t343;
t270 = t294 * t342;
t271 = t294 * t344;
t490 = m(3) * (t170 * t270 - t171 * t271);
t488 = m(5) * (t159 * t342 - t160 * t344);
t487 = m(5) * (-t342 * t155 - t157 * t344);
t486 = m(5) * (t172 * t342 - t173 * t344);
t484 = m(5) * (-t342 * t190 - t192 * t344);
t476 = m(5) * qJD(2);
t449 = t244 * t341;
t428 = t343 * t344;
t242 = Icges(3,5) * t429 - Icges(3,6) * t433 - Icges(3,3) * t344;
t413 = -t342 * t242 - t246 * t428;
t366 = Icges(3,5) * t343 - Icges(3,6) * t341;
t243 = Icges(3,3) * t342 + t344 * t366;
t412 = t342 * t243 + t247 * t428;
t392 = t75 * t126 + t420;
t390 = -t286 - t336;
t387 = ((-t532 * t342 + t526) * t344 + t531 * t338) * t493 + ((-t531 * t344 + t526) * t342 + t532 * t339) * t491;
t194 = t247 * t429;
t375 = t243 * t344 - t194;
t372 = t245 * t341 - t242;
t371 = t394 * t483;
t370 = t68 + t387;
t365 = -Icges(3,5) * t341 - Icges(3,6) * t343;
t355 = t388 - t336;
t347 = -t353 + (t551 * t321 + t553 * t322 + t577 * t342 + t530 * t344) * t493 + (-t321 * t550 + t322 * t552 + t530 * t342 - t567) * t491;
t346 = -t351 + (t544 * t321 + t562 * t322) * (t493 + t494);
t296 = -rSges(3,2) * t341 + t475;
t265 = t365 * t344;
t264 = t365 * t342;
t204 = t390 * t344;
t202 = t390 * t342;
t158 = t355 * t344;
t156 = t355 * t342;
t135 = -t371 + t146;
t133 = -t191 * t344 + t193 * t342;
t121 = t484 / 0.2e1;
t115 = m(5) * t133 * qJD(3);
t114 = t486 / 0.2e1;
t108 = -t245 * t432 + t412;
t107 = -t244 * t432 - t413;
t106 = -t245 * t433 - t375;
t104 = t487 / 0.2e1;
t103 = t488 / 0.2e1;
t83 = -t371 + t92;
t72 = -t107 * t344 + t108 * t342;
t71 = -(-t342 * (-t246 * t343 + t449) - t242 * t344) * t344 + t106 * t342;
t52 = t121 - t486 / 0.2e1;
t51 = t121 + t114;
t50 = t114 - t484 / 0.2e1;
t45 = t104 + t103;
t44 = t104 - t488 / 0.2e1;
t43 = t103 - t487 / 0.2e1;
t26 = t351 + t495 + t502;
t25 = (t106 - t194 + (t243 + t449) * t344 + t413) * t344 + t412 * t342;
t24 = (t344 * t372 + t108 - t412) * t344 + (t342 * t372 + t107 + t375) * t342;
t10 = (t292 / 0.2e1 + t291 / 0.2e1) * t343 + t533 + t490 + t503 + t496 + t351;
t7 = t370 + t499;
t6 = t387 - t560;
t4 = (t72 / 0.2e1 - t25 / 0.2e1) * t344 + (t24 / 0.2e1 + t71 / 0.2e1) * t342 + t353;
t3 = t353 + t528;
t2 = t353 - t528;
t1 = t347 + t478 + t479;
t5 = [t10 * qJD(2) + t26 * qJD(3) + qJD(4) * t538, t10 * qJD(1) + t1 * qJD(3) + t45 * qJD(4) + ((t144 * t158 + t145 * t156 - t155 * t160 - t157 * t159) * t510 + (t151 * t204 + t152 * t202) * t511 + m(3) * ((-t170 * t344 - t171 * t342) * t296 + (-t270 * t344 + t271 * t342) * t294) / 0.2e1) * t514 + (t347 + t25 * t534 + (t341 * t401 + t343 * t399) * t493 + (t24 + t71) * t494 + (-t341 * t402 + t343 * t400 + t72) * t491 + (t338 / 0.2e1 + t339 / 0.2e1) * t366) * qJD(2), t26 * qJD(1) + t1 * qJD(2) + t347 * qJD(3) + t51 * qJD(4) + ((t354 + (-t261 * t344 + t262 * t342) * t284) * t511 + (-t172 * t192 - t173 * t190 + t422) * t510) * t512, t45 * qJD(2) + t51 * qJD(3) + t523; (t346 - (t292 + t291) * t343 / 0.2e1 - t533) * qJD(1) + t4 * qJD(2) + t2 * qJD(3) + t44 * qJD(4) + (-t496 / 0.4e1 - t503 / 0.4e1 - t490 / 0.4e1) * t515, t4 * qJD(1) + (m(5) * (-t155 * t156 - t157 * t158 + t62 * t83) + m(4) * (t135 * t88 - t202 * t522 - t204 * t521) + (t338 * t265 + (-t342 * t264 + t518) * t344) * t493 + m(3) * ((t342 * (rSges(3,1) * t429 - t395) + t344 * (rSges(3,1) * t428 + rSges(3,3) * t342 - t314)) * (-t270 * t342 - t271 * t344) + t394 * t296 * t294) + (t339 * t264 + (-t344 * t265 + t518) * t342) * t491 + t387) * qJD(2) + t6 * qJD(3), t2 * qJD(1) + t6 * qJD(2) + ((t392 + t30) * t510 + (t42 + t69) * t511) * t512 + (t387 - t68 - t499) * qJD(3), t44 * qJD(1); t346 * qJD(1) + t3 * qJD(2) + t353 * qJD(3) + t52 * qJD(4) + (-t502 / 0.4e1 - t495 / 0.4e1) * t515, t3 * qJD(1) + t7 * qJD(3) + ((-t156 * t190 - t158 * t192 + t62 * t92 + t75 * t83 + t421) * t510 + (t135 * t136 + (-t202 * t342 - t204 * t344) * t284 + t42) * t511) * t514 + (t387 + t560) * qJD(2), t353 * qJD(1) + t7 * qJD(2) + (m(5) * t392 + t370) * qJD(3), t52 * qJD(1); t43 * qJD(2) + t50 * qJD(3) - t523, t43 * qJD(1) + (-t156 * t344 + t158 * t342) * t476 + t115, t50 * qJD(1) + t133 * t476 + t115, 0;];
Cq = t5;
