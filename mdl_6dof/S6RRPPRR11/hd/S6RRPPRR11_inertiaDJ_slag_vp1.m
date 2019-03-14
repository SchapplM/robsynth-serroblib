% Calculate time derivative of joint inertia matrix for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:06
% EndTime: 2019-03-09 09:39:52
% DurationCPUTime: 27.86s
% Computational Cost: add. (68991->1308), mult. (133818->1751), div. (0->0), fcn. (147512->12), ass. (0->516)
t487 = cos(pkin(6));
t491 = sin(qJ(1));
t493 = cos(qJ(2));
t597 = t491 * t493;
t490 = sin(qJ(2));
t494 = cos(qJ(1));
t598 = t490 * t494;
t454 = t487 * t597 + t598;
t599 = t490 * t491;
t560 = t487 * t599;
t596 = t493 * t494;
t455 = -t560 + t596;
t485 = sin(pkin(6));
t603 = t485 * t491;
t340 = Icges(4,4) * t603 - Icges(4,2) * t455 + Icges(4,6) * t454;
t349 = Icges(3,1) * t455 - Icges(3,4) * t454 + Icges(3,5) * t603;
t659 = t340 - t349;
t561 = t487 * t596;
t452 = -t561 + t599;
t453 = t487 * t598 + t597;
t601 = t485 * t494;
t341 = -Icges(4,4) * t601 - Icges(4,2) * t453 + Icges(4,6) * t452;
t348 = Icges(3,1) * t453 - Icges(3,4) * t452 - Icges(3,5) * t601;
t658 = t341 - t348;
t342 = Icges(4,1) * t603 - Icges(4,4) * t455 + Icges(4,5) * t454;
t345 = Icges(3,5) * t455 - Icges(3,6) * t454 + Icges(3,3) * t603;
t655 = t342 + t345;
t343 = -Icges(4,1) * t601 - Icges(4,4) * t453 + Icges(4,5) * t452;
t344 = Icges(3,5) * t453 - Icges(3,6) * t452 - Icges(3,3) * t601;
t648 = t343 + t344;
t339 = -Icges(4,5) * t601 - Icges(4,6) * t453 + Icges(4,3) * t452;
t346 = Icges(3,4) * t453 - Icges(3,2) * t452 - Icges(3,6) * t601;
t647 = -t346 + t339;
t338 = Icges(4,5) * t603 - Icges(4,6) * t455 + Icges(4,3) * t454;
t347 = Icges(3,4) * t455 - Icges(3,2) * t454 + Icges(3,6) * t603;
t646 = t347 - t338;
t484 = sin(pkin(11));
t486 = cos(pkin(11));
t602 = t485 * t493;
t450 = -t484 * t487 - t486 * t602;
t451 = -t484 * t602 + t486 * t487;
t604 = t485 * t490;
t333 = Icges(5,5) * t451 + Icges(5,6) * t450 + Icges(5,3) * t604;
t334 = Icges(5,4) * t451 + Icges(5,2) * t450 + Icges(5,6) * t604;
t335 = Icges(5,1) * t451 + Icges(5,4) * t450 + Icges(5,5) * t604;
t573 = qJD(2) * t485;
t385 = (Icges(5,3) * t493 + (Icges(5,5) * t484 + Icges(5,6) * t486) * t490) * t573;
t386 = (Icges(5,6) * t493 + (Icges(5,4) * t484 + Icges(5,2) * t486) * t490) * t573;
t387 = (Icges(5,5) * t493 + (Icges(5,1) * t484 + Icges(5,4) * t486) * t490) * t573;
t612 = Icges(3,4) * t490;
t418 = Icges(3,6) * t487 + (Icges(3,2) * t493 + t612) * t485;
t611 = Icges(3,4) * t493;
t419 = Icges(3,5) * t487 + (Icges(3,1) * t490 + t611) * t485;
t610 = Icges(4,6) * t490;
t420 = Icges(4,5) * t487 + (-Icges(4,3) * t493 - t610) * t485;
t432 = (Icges(3,5) * t493 - Icges(3,6) * t490) * t573;
t433 = (-Icges(4,4) * t493 + Icges(4,5) * t490) * t573;
t434 = (-Icges(3,2) * t490 + t611) * t573;
t435 = (Icges(3,1) * t493 - t612) * t573;
t571 = qJD(2) * t493;
t542 = t485 * t571;
t572 = qJD(2) * t490;
t543 = t485 * t572;
t657 = t450 * t386 + t451 * t387 + t434 * t602 + (t435 + t385) * t604 + (t419 + t333) * t542 + (t432 + t433) * t487 + (t334 * t486 + t335 * t484 - t418 + t420) * t543;
t536 = qJD(2) * t487 + qJD(1);
t373 = -qJD(1) * t561 - t494 * t571 + t536 * t599;
t374 = qJD(1) * t453 + qJD(2) * t454;
t574 = qJD(1) * t494;
t544 = t485 * t574;
t269 = Icges(4,1) * t544 + Icges(4,4) * t374 - Icges(4,5) * t373;
t270 = -Icges(3,5) * t374 + Icges(3,6) * t373 + Icges(3,3) * t544;
t656 = t269 + t270;
t654 = -t646 * t454 - t455 * t659 + t655 * t603;
t653 = -t452 * t647 + t453 * t658 + t601 * t648;
t652 = t646 * t452 + t453 * t659 + t655 * t601;
t651 = t454 * t647 - t455 * t658 + t603 * t648;
t265 = Icges(4,5) * t544 + Icges(4,6) * t374 - Icges(4,3) * t373;
t272 = -Icges(3,4) * t374 + Icges(3,2) * t373 + Icges(3,6) * t544;
t650 = -t272 + t265;
t375 = qJD(1) * t454 + qJD(2) * t453;
t376 = -qJD(1) * t560 - t491 * t572 + t536 * t596;
t575 = qJD(1) * t491;
t545 = t485 * t575;
t264 = Icges(4,5) * t545 - Icges(4,6) * t376 + Icges(4,3) * t375;
t273 = Icges(3,4) * t376 - Icges(3,2) * t375 + Icges(3,6) * t545;
t649 = t273 - t264;
t569 = pkin(11) + qJ(5);
t481 = sin(t569);
t540 = cos(t569);
t424 = -t481 * t602 + t487 * t540;
t609 = Icges(4,6) * t493;
t421 = Icges(4,4) * t487 + (-Icges(4,2) * t490 - t609) * t485;
t430 = (Icges(4,3) * t490 - t609) * t573;
t431 = (-Icges(4,2) * t493 + t610) * t573;
t600 = t490 * t431;
t645 = ((-t600 + (-qJD(2) * t421 - t430) * t493) * t485 + t657) * t487;
t328 = -t373 * t486 - t484 * t544;
t608 = t373 * t484;
t329 = t486 * t544 - t608;
t206 = Icges(5,5) * t329 + Icges(5,6) * t328 - Icges(5,3) * t374;
t267 = Icges(4,4) * t544 + Icges(4,2) * t374 - Icges(4,6) * t373;
t274 = -Icges(3,1) * t374 + Icges(3,4) * t373 + Icges(3,5) * t544;
t644 = t274 - t267 + t206;
t326 = t375 * t486 - t484 * t545;
t607 = t375 * t484;
t327 = t486 * t545 + t607;
t205 = Icges(5,5) * t327 + Icges(5,6) * t326 + Icges(5,3) * t376;
t266 = Icges(4,4) * t545 - Icges(4,2) * t376 + Icges(4,6) * t375;
t275 = Icges(3,1) * t376 - Icges(3,4) * t375 + Icges(3,5) * t545;
t643 = t275 - t266 + t205;
t403 = t452 * t486 + t484 * t601;
t606 = t452 * t484;
t404 = -t486 * t601 + t606;
t289 = Icges(5,5) * t404 + Icges(5,6) * t403 + Icges(5,3) * t453;
t642 = t289 - t658;
t401 = t454 * t486 - t484 * t603;
t605 = t454 * t484;
t402 = t486 * t603 + t605;
t288 = Icges(5,5) * t402 + Icges(5,6) * t401 + Icges(5,3) * t455;
t641 = -t288 + t659;
t268 = Icges(4,1) * t545 - Icges(4,4) * t376 + Icges(4,5) * t375;
t271 = Icges(3,5) * t376 - Icges(3,6) * t375 + Icges(3,3) * t545;
t640 = (-t268 - t271) * t494;
t639 = 2 * m(3);
t638 = 2 * m(4);
t637 = 2 * m(5);
t636 = 2 * m(6);
t635 = 2 * m(7);
t483 = t485 ^ 2;
t520 = t485 * t540;
t511 = t494 * t520;
t246 = -t375 * t540 - qJD(5) * t511 + (qJD(5) * t452 + t545) * t481;
t634 = t246 / 0.2e1;
t391 = t454 * t481 + t491 * t520;
t248 = qJD(5) * t391 + t373 * t540 + t481 * t544;
t633 = t248 / 0.2e1;
t378 = qJD(5) * t424 - t520 * t572;
t632 = t378 / 0.2e1;
t390 = -t454 * t540 + t481 * t603;
t631 = t390 / 0.2e1;
t392 = t452 * t540 + t481 * t601;
t630 = -t392 / 0.2e1;
t499 = -t487 * t481 - t493 * t520;
t629 = -t499 / 0.2e1;
t628 = t487 / 0.2e1;
t627 = t491 / 0.2e1;
t626 = rSges(4,2) - pkin(2);
t625 = -rSges(5,3) - pkin(2);
t624 = rSges(7,3) + pkin(10);
t623 = pkin(4) * t484;
t393 = t452 * t481 - t511;
t622 = pkin(5) * t393;
t505 = qJD(1) * t520;
t247 = qJD(5) * t392 + t375 * t481 + t491 * t505;
t621 = t247 * pkin(5);
t482 = t494 * pkin(1);
t488 = -pkin(9) - qJ(4);
t620 = -pkin(2) + t488;
t480 = pkin(4) * t486 + pkin(3);
t619 = -pkin(3) + t480;
t618 = t375 * rSges(4,3);
t489 = sin(qJ(6));
t492 = cos(qJ(6));
t388 = -t424 * t489 + t492 * t604;
t389 = t424 * t492 + t489 * t604;
t237 = Icges(7,5) * t389 + Icges(7,6) * t388 - Icges(7,3) * t499;
t238 = Icges(7,4) * t389 + Icges(7,2) * t388 - Icges(7,6) * t499;
t239 = Icges(7,1) * t389 + Icges(7,4) * t388 - Icges(7,5) * t499;
t110 = -t237 * t499 + t238 * t388 + t239 * t389;
t377 = qJD(5) * t499 + t481 * t543;
t262 = -qJD(6) * t389 - t377 * t489 + t492 * t542;
t263 = qJD(6) * t388 + t377 * t492 + t489 * t542;
t155 = Icges(7,5) * t263 + Icges(7,6) * t262 + Icges(7,3) * t378;
t156 = Icges(7,4) * t263 + Icges(7,2) * t262 + Icges(7,6) * t378;
t157 = Icges(7,1) * t263 + Icges(7,4) * t262 + Icges(7,5) * t378;
t48 = -t155 * t499 + t388 * t156 + t389 * t157 + t378 * t237 + t262 * t238 + t263 * t239;
t617 = t110 * t378 - t48 * t499;
t616 = t110 * t542 + t48 * t604;
t319 = Icges(6,5) * t424 + Icges(6,6) * t499 + Icges(6,3) * t604;
t320 = Icges(6,4) * t424 + Icges(6,2) * t499 + Icges(6,6) * t604;
t321 = Icges(6,1) * t424 + Icges(6,4) * t499 + Icges(6,5) * t604;
t162 = t319 * t604 + t320 * t499 + t321 * t424;
t276 = Icges(6,5) * t377 - Icges(6,6) * t378 + Icges(6,3) * t542;
t277 = Icges(6,4) * t377 - Icges(6,2) * t378 + Icges(6,6) * t542;
t278 = Icges(6,1) * t377 - Icges(6,4) * t378 + Icges(6,5) * t542;
t77 = t276 * t604 + t277 * t499 + t424 * t278 + t319 * t542 - t378 * t320 + t377 * t321;
t615 = t162 * t542 + t77 * t604;
t314 = t393 * t492 + t453 * t489;
t167 = -qJD(6) * t314 - t247 * t489 + t376 * t492;
t313 = -t393 * t489 + t453 * t492;
t168 = qJD(6) * t313 + t247 * t492 + t376 * t489;
t522 = -t168 * rSges(7,1) - t167 * rSges(7,2);
t92 = t246 * rSges(7,3) - t522;
t614 = t246 * pkin(10) + t621 + t92;
t249 = -qJD(5) * t390 - t373 * t481 + t494 * t505;
t312 = t391 * t492 + t455 * t489;
t169 = -qJD(6) * t312 - t249 * t489 - t374 * t492;
t311 = -t391 * t489 + t455 * t492;
t170 = qJD(6) * t311 + t249 * t492 - t374 * t489;
t93 = t170 * rSges(7,1) + t169 * rSges(7,2) + t248 * rSges(7,3);
t613 = t249 * pkin(5) + pkin(10) * t248 + t93;
t595 = -qJ(4) - t488;
t158 = rSges(7,1) * t263 + rSges(7,2) * t262 + rSges(7,3) * t378;
t594 = pkin(5) * t377 + pkin(10) * t378 + t158;
t194 = t312 * rSges(7,1) + t311 * rSges(7,2) + t390 * rSges(7,3);
t593 = t391 * pkin(5) + pkin(10) * t390 + t194;
t521 = -rSges(7,1) * t314 - rSges(7,2) * t313;
t195 = -rSges(7,3) * t392 - t521;
t592 = -pkin(10) * t392 + t195 + t622;
t235 = -t374 * pkin(2) - qJ(3) * t373 + qJD(3) * t454;
t234 = t487 * t235;
t538 = -pkin(3) * t544 + qJ(4) * t374;
t570 = qJD(4) * t455;
t307 = -t538 + t570;
t591 = t487 * t307 + t234;
t586 = -t375 * qJ(3) - t452 * qJD(3);
t236 = t376 * pkin(2) - t586;
t365 = t376 * qJ(4);
t439 = t453 * qJD(4);
t306 = pkin(3) * t545 + t365 + t439;
t590 = -t236 - t306;
t240 = rSges(7,1) * t389 + rSges(7,2) * t388 - rSges(7,3) * t499;
t589 = pkin(5) * t424 - pkin(10) * t499 + t240;
t440 = t452 * qJ(3);
t379 = t453 * pkin(2) + t440;
t380 = t455 * pkin(2) + qJ(3) * t454;
t588 = t379 * t603 + t380 * t601;
t363 = t487 * t380;
t409 = pkin(3) * t603 + qJ(4) * t455;
t587 = t487 * t409 + t363;
t576 = pkin(3) * t601 - t453 * qJ(4);
t585 = -t379 + t576;
t584 = -t380 - t409;
t408 = (-qJD(3) * t493 + (pkin(2) * t493 + qJ(3) * t490) * qJD(2)) * t485;
t582 = -t408 - (-rSges(4,2) * t493 + rSges(4,3) * t490) * t573;
t581 = -t408 - (qJ(4) * t571 + qJD(4) * t490) * t485;
t456 = (pkin(2) * t490 - qJ(3) * t493) * t485;
t413 = t456 * t545;
t458 = pkin(3) * t487 + qJ(4) * t604;
t580 = t458 * t545 + t413;
t427 = rSges(4,1) * t487 + (-rSges(4,2) * t490 - rSges(4,3) * t493) * t485;
t579 = -t427 - t456;
t578 = -t456 - t458;
t577 = pkin(8) * t603 + t482;
t568 = -rSges(6,3) + t620;
t567 = pkin(4) * t607;
t189 = Icges(7,5) * t314 + Icges(7,6) * t313 - Icges(7,3) * t392;
t191 = Icges(7,4) * t314 + Icges(7,2) * t313 - Icges(7,6) * t392;
t193 = Icges(7,1) * t314 + Icges(7,4) * t313 - Icges(7,5) * t392;
t86 = Icges(7,5) * t168 + Icges(7,6) * t167 + Icges(7,3) * t246;
t88 = Icges(7,4) * t168 + Icges(7,2) * t167 + Icges(7,6) * t246;
t90 = Icges(7,1) * t168 + Icges(7,4) * t167 + Icges(7,5) * t246;
t22 = t189 * t378 + t191 * t262 + t193 * t263 + t388 * t88 + t389 * t90 - t499 * t86;
t34 = -t155 * t392 + t156 * t313 + t157 * t314 + t167 * t238 + t168 * t239 + t237 * t246;
t566 = -t34 / 0.2e1 - t22 / 0.2e1;
t188 = Icges(7,5) * t312 + Icges(7,6) * t311 + Icges(7,3) * t390;
t190 = Icges(7,4) * t312 + Icges(7,2) * t311 + Icges(7,6) * t390;
t192 = Icges(7,1) * t312 + Icges(7,4) * t311 + Icges(7,5) * t390;
t87 = Icges(7,5) * t170 + Icges(7,6) * t169 + Icges(7,3) * t248;
t89 = Icges(7,4) * t170 + Icges(7,2) * t169 + Icges(7,6) * t248;
t91 = Icges(7,1) * t170 + Icges(7,4) * t169 + Icges(7,5) * t248;
t21 = t188 * t378 + t190 * t262 + t192 * t263 + t388 * t89 + t389 * t91 - t499 * t87;
t35 = t155 * t390 + t156 * t311 + t157 * t312 + t169 * t238 + t170 * t239 + t237 * t248;
t565 = t35 / 0.2e1 + t21 / 0.2e1;
t78 = -t188 * t499 + t190 * t388 + t192 * t389;
t98 = t237 * t390 + t238 * t311 + t239 * t312;
t564 = t98 / 0.2e1 + t78 / 0.2e1;
t79 = -t189 * t499 + t191 * t388 + t193 * t389;
t99 = -t237 * t392 + t238 * t313 + t239 * t314;
t563 = t99 / 0.2e1 + t79 / 0.2e1;
t550 = -pkin(4) * t608 + t374 * t488 + t480 * t544;
t216 = t538 + t550;
t559 = t487 * t216 + t591;
t215 = -t376 * t488 + t545 * t619 - t365 + t567;
t558 = -t215 + t590;
t557 = t235 * t601 + t236 * t603 + t379 * t544;
t153 = t249 * rSges(6,1) - t248 * rSges(6,2) - t374 * rSges(6,3);
t546 = pkin(4) * t605 - t455 * t488 + t480 * t603;
t296 = -t409 + t546;
t556 = t487 * t296 + t587;
t555 = -t296 + t584;
t531 = -pkin(4) * t606 + t480 * t601;
t297 = -t453 * t488 - t531 + t576;
t554 = -t297 + t585;
t212 = t329 * rSges(5,1) + t328 * rSges(5,2) - t374 * rSges(5,3);
t355 = t619 * t487 + (t490 * t595 - t493 * t623) * t485;
t553 = t355 * t545 + t580;
t336 = rSges(5,1) * t451 + rSges(5,2) * t450 + rSges(5,3) * t604;
t552 = -t336 + t578;
t551 = -t355 + t578;
t549 = -t439 + t586;
t281 = -t374 * rSges(3,1) + t373 * rSges(3,2) + rSges(3,3) * t544;
t260 = t391 * rSges(6,1) - t390 * rSges(6,2) + t455 * rSges(6,3);
t548 = -(rSges(5,3) * t493 + (rSges(5,1) * t484 + rSges(5,2) * t486) * t490) * t573 + t581;
t547 = -(t490 * t623 + t493 * t595) * t573 + t581;
t294 = t402 * rSges(5,1) + t401 * rSges(5,2) + t455 * rSges(5,3);
t354 = t455 * rSges(3,1) - t454 * rSges(3,2) + rSges(3,3) * t603;
t280 = rSges(4,1) * t544 + t374 * rSges(4,2) - t373 * rSges(4,3);
t351 = rSges(4,1) * t603 - t455 * rSges(4,2) + t454 * rSges(4,3);
t541 = -t491 * pkin(1) + pkin(8) * t601;
t539 = t494 * t579;
t537 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t283 = rSges(6,1) * t377 - rSges(6,2) * t378 + rSges(6,3) * t542;
t535 = -t283 + t547;
t325 = rSges(6,1) * t424 + rSges(6,2) * t499 + rSges(6,3) * t604;
t534 = -t325 + t551;
t533 = t409 * t601 - t576 * t603 + t588;
t532 = t483 * t490 * t571;
t530 = -t440 + t541;
t529 = -pkin(1) * t575 + pkin(8) * t544;
t528 = t552 * t494;
t527 = -t376 * rSges(3,1) + t375 * rSges(3,2);
t526 = -t327 * rSges(5,1) - t326 * rSges(5,2);
t525 = -rSges(5,1) * t404 - rSges(5,2) * t403;
t524 = -t247 * rSges(6,1) + t246 * rSges(6,2);
t523 = -rSges(6,1) * t393 - rSges(6,2) * t392;
t518 = t547 - t594;
t517 = t551 - t589;
t516 = t534 * t494;
t515 = t380 + t577;
t146 = Icges(6,5) * t247 - Icges(6,6) * t246 + Icges(6,3) * t376;
t148 = Icges(6,4) * t247 - Icges(6,2) * t246 + Icges(6,6) * t376;
t150 = Icges(6,1) * t247 - Icges(6,4) * t246 + Icges(6,5) * t376;
t255 = Icges(6,5) * t393 + Icges(6,6) * t392 + Icges(6,3) * t453;
t257 = Icges(6,4) * t393 + Icges(6,2) * t392 + Icges(6,6) * t453;
t259 = Icges(6,1) * t393 + Icges(6,4) * t392 + Icges(6,5) * t453;
t52 = t148 * t499 + t150 * t424 - t257 * t378 + t259 * t377 + (t146 * t490 + t255 * t571) * t485;
t66 = -t246 * t320 + t247 * t321 + t276 * t453 + t277 * t392 + t278 * t393 + t319 * t376;
t514 = t66 / 0.2e1 + t52 / 0.2e1 - t566;
t147 = Icges(6,5) * t249 - Icges(6,6) * t248 - Icges(6,3) * t374;
t149 = Icges(6,4) * t249 - Icges(6,2) * t248 - Icges(6,6) * t374;
t151 = Icges(6,1) * t249 - Icges(6,4) * t248 - Icges(6,5) * t374;
t254 = Icges(6,5) * t391 - Icges(6,6) * t390 + Icges(6,3) * t455;
t256 = Icges(6,4) * t391 - Icges(6,2) * t390 + Icges(6,6) * t455;
t258 = Icges(6,1) * t391 - Icges(6,4) * t390 + Icges(6,5) * t455;
t51 = t149 * t499 + t151 * t424 - t256 * t378 + t258 * t377 + (t147 * t490 + t254 * t571) * t485;
t67 = -t248 * t320 + t249 * t321 + t276 * t455 - t277 * t390 + t278 * t391 - t319 * t374;
t513 = t67 / 0.2e1 + t51 / 0.2e1 + t565;
t512 = rSges(4,1) * t601 - t452 * rSges(4,3);
t125 = t254 * t604 + t256 * t499 + t258 * t424;
t138 = t319 * t455 - t320 * t390 + t321 * t391;
t510 = -t138 / 0.2e1 - t125 / 0.2e1 - t564;
t126 = t255 * t604 + t257 * t499 + t259 * t424;
t139 = t319 * t453 + t320 * t392 + t321 * t393;
t509 = t139 / 0.2e1 + t126 / 0.2e1 + t563;
t508 = t306 * t603 + t307 * t601 - t544 * t576 + t557;
t507 = t296 * t601 + t297 * t603 + t533;
t504 = t517 * t494;
t353 = t453 * rSges(3,1) - t452 * rSges(3,2) - rSges(3,3) * t601;
t503 = t530 + t531;
t502 = t515 + t546;
t500 = t215 * t603 + t216 * t601 + t297 * t544 + t508;
t498 = t235 + t529;
t497 = t498 + t570;
t496 = (-t482 + (-pkin(8) - t480) * t603) * qJD(1) - t567 + t549;
t495 = t497 + t550;
t437 = (rSges(3,1) * t493 - rSges(3,2) * t490) * t573;
t426 = rSges(3,3) * t487 + (rSges(3,1) * t490 + rSges(3,2) * t493) * t485;
t422 = Icges(4,1) * t487 + (-Icges(4,4) * t490 - Icges(4,5) * t493) * t485;
t417 = Icges(3,3) * t487 + (Icges(3,5) * t490 + Icges(3,6) * t493) * t485;
t352 = -t453 * rSges(4,2) - t512;
t318 = t354 + t577;
t317 = -t353 + t541;
t300 = -t487 * t353 - t426 * t601;
t299 = t354 * t487 - t426 * t603;
t295 = rSges(5,3) * t453 - t525;
t293 = Icges(5,1) * t404 + Icges(5,4) * t403 + Icges(5,5) * t453;
t292 = Icges(5,1) * t402 + Icges(5,4) * t401 + Icges(5,5) * t455;
t291 = Icges(5,4) * t404 + Icges(5,2) * t403 + Icges(5,6) * t453;
t290 = Icges(5,4) * t402 + Icges(5,2) * t401 + Icges(5,6) * t455;
t282 = rSges(3,3) * t545 - t527;
t279 = rSges(4,1) * t545 - t376 * rSges(4,2) + t618;
t261 = rSges(6,3) * t453 - t523;
t253 = t515 + t351;
t252 = t453 * t626 + t512 + t530;
t233 = (-t482 + (-rSges(3,3) - pkin(8)) * t603) * qJD(1) + t527;
t232 = t529 + t281;
t229 = t417 * t603 - t418 * t454 + t419 * t455;
t228 = -t417 * t601 - t452 * t418 + t453 * t419;
t227 = t452 * t420 - t453 * t421 - t422 * t601;
t226 = t420 * t454 - t421 * t455 + t422 * t603;
t218 = (-t352 - t379) * t487 + t485 * t539;
t217 = t351 * t487 + t579 * t603 + t363;
t214 = t487 * t281 + (-t426 * t574 - t437 * t491) * t485;
t213 = -t487 * t282 + (t426 * t575 - t437 * t494) * t485;
t211 = t376 * rSges(5,3) - t526;
t210 = Icges(5,1) * t329 + Icges(5,4) * t328 - Icges(5,5) * t374;
t209 = Icges(5,1) * t327 + Icges(5,4) * t326 + Icges(5,5) * t376;
t208 = Icges(5,4) * t329 + Icges(5,2) * t328 - Icges(5,6) * t374;
t207 = Icges(5,4) * t327 + Icges(5,2) * t326 + Icges(5,6) * t376;
t203 = t343 * t487 + (-t339 * t493 - t341 * t490) * t485;
t202 = t342 * t487 + (-t338 * t493 - t340 * t490) * t485;
t201 = t345 * t487 + (t347 * t493 + t349 * t490) * t485;
t200 = t344 * t487 + (t346 * t493 + t348 * t490) * t485;
t197 = t409 + t515 + t294;
t196 = t453 * t625 + t525 + t530 + t576;
t187 = t260 * t604 - t325 * t455;
t186 = -t261 * t604 + t325 * t453;
t175 = (t351 * t494 + t352 * t491) * t485 + t588;
t172 = t502 + t260;
t171 = t453 * t568 + t503 + t523;
t166 = -t618 + t626 * t376 + (-t482 + (-rSges(4,1) - pkin(8)) * t603) * qJD(1) + t586;
t165 = t498 + t280;
t161 = -t260 * t453 + t261 * t455;
t160 = t333 * t453 + t334 * t403 + t335 * t404;
t159 = t333 * t455 + t334 * t401 + t335 * t402;
t152 = t376 * rSges(6,3) - t524;
t143 = -t375 * t418 + t376 * t419 - t452 * t434 + t453 * t435 + (t417 * t575 - t432 * t494) * t485;
t142 = t373 * t418 - t374 * t419 - t454 * t434 + t455 * t435 + (t417 * t574 + t432 * t491) * t485;
t141 = -t373 * t420 + t374 * t421 + t454 * t430 - t455 * t431 + (t422 * t574 + t433 * t491) * t485;
t140 = t375 * t420 - t376 * t421 + t452 * t430 - t453 * t431 + (t422 * t575 - t433 * t494) * t485;
t137 = (-t295 + t585) * t487 + t485 * t528;
t136 = t294 * t487 + t552 * t603 + t587;
t133 = t289 * t604 + t291 * t450 + t293 * t451;
t132 = t288 * t604 + t290 * t450 + t292 * t451;
t131 = t487 * t280 + t234 + (qJD(1) * t539 + t491 * t582) * t485;
t130 = t413 + (-t236 - t279) * t487 + (t427 * t575 + t494 * t582) * t485;
t129 = t195 * t499 - t240 * t392;
t128 = -t194 * t499 - t240 * t390;
t124 = -t365 + t625 * t376 + (-t482 + (-pkin(3) - pkin(8)) * t603) * qJD(1) + t526 + t549;
t123 = t497 - t538 + t212;
t122 = t502 + t593;
t121 = t392 * t624 + t453 * t620 + t503 + t521 - t622;
t120 = t289 * t453 + t291 * t403 + t293 * t404;
t119 = t288 * t453 + t290 * t403 + t292 * t404;
t118 = t289 * t455 + t291 * t401 + t293 * t402;
t117 = t288 * t455 + t290 * t401 + t292 * t402;
t116 = (t294 * t494 + t295 * t491) * t485 + t533;
t115 = t255 * t453 + t257 * t392 + t259 * t393;
t114 = t254 * t453 + t256 * t392 + t258 * t393;
t113 = t255 * t455 - t257 * t390 + t259 * t391;
t112 = t254 * t455 - t256 * t390 + t258 * t391;
t111 = t194 * t392 + t195 * t390;
t109 = t269 * t487 + (-t265 * t493 - t267 * t490 + (t338 * t490 - t340 * t493) * qJD(2)) * t485;
t108 = t268 * t487 + (-t264 * t493 - t266 * t490 + (t339 * t490 - t341 * t493) * qJD(2)) * t485;
t107 = t270 * t487 + (t272 * t493 + t274 * t490 + (-t347 * t490 + t349 * t493) * qJD(2)) * t485;
t106 = t271 * t487 + (t273 * t493 + t275 * t490 + (-t346 * t490 + t348 * t493) * qJD(2)) * t485;
t104 = (-t261 + t554) * t487 + t485 * t516;
t103 = t260 * t487 + t534 * t603 + t556;
t102 = -t455 * t589 + t593 * t604;
t101 = t453 * t589 - t592 * t604;
t97 = t376 * t568 + t496 + t524;
t96 = t495 + t153;
t95 = t328 * t334 + t329 * t335 - t333 * t374 + t385 * t455 + t386 * t401 + t387 * t402;
t94 = t326 * t334 + t327 * t335 + t333 * t376 + t385 * t453 + t386 * t403 + t387 * t404;
t85 = (t260 * t494 + t261 * t491) * t485 + t507;
t84 = t487 * t212 + (qJD(1) * t528 + t491 * t548) * t485 + t591;
t83 = (-t211 + t590) * t487 + (t336 * t575 + t494 * t548) * t485 + t580;
t82 = -t453 * t593 + t455 * t592;
t81 = t283 * t453 + t325 * t376 + (-t152 * t490 - t261 * t571) * t485;
t80 = -t283 * t455 + t325 * t374 + (t153 * t490 + t260 * t571) * t485;
t76 = t77 * t487;
t74 = (t279 * t491 + t280 * t494 + (t352 * t494 + (-t351 - t380) * t491) * qJD(1)) * t485 + t557;
t73 = -t189 * t392 + t191 * t313 + t193 * t314;
t72 = -t188 * t392 + t190 * t313 + t192 * t314;
t71 = t189 * t390 + t191 * t311 + t193 * t312;
t70 = t188 * t390 + t190 * t311 + t192 * t312;
t69 = (t554 - t592) * t487 + t485 * t504;
t68 = t487 * t593 + t517 * t603 + t556;
t65 = t152 * t455 - t153 * t453 - t260 * t376 - t261 * t374;
t64 = t207 * t450 + t209 * t451 + (t205 * t490 + (t289 * t493 + (t291 * t486 + t293 * t484) * t490) * qJD(2)) * t485;
t63 = t208 * t450 + t210 * t451 + (t206 * t490 + (t288 * t493 + (t290 * t486 + t292 * t484) * t490) * qJD(2)) * t485;
t62 = (t491 * t592 + t494 * t593) * t485 + t507;
t61 = t487 * t153 + (qJD(1) * t516 + t491 * t535) * t485 + t559;
t60 = (-t152 + t558) * t487 + (t325 * t575 + t494 * t535) * t485 + t553;
t59 = -t246 * t624 + t376 * t620 + t496 + t522 - t621;
t58 = t495 + t613;
t57 = t139 * t487 + (t114 * t491 - t115 * t494) * t485;
t56 = t138 * t487 + (t112 * t491 - t113 * t494) * t485;
t55 = (t211 * t491 + t212 * t494 + (t295 * t494 + (-t294 + t584) * t491) * qJD(1)) * t485 + t508;
t54 = t114 * t455 + t115 * t453 + t139 * t604;
t53 = t112 * t455 + t113 * t453 + t138 * t604;
t50 = -t158 * t390 + t194 * t378 - t240 * t248 - t499 * t93;
t49 = -t158 * t392 - t195 * t378 + t240 * t246 + t499 * t92;
t47 = t48 * t487;
t45 = t146 * t455 - t148 * t390 + t150 * t391 - t248 * t257 + t249 * t259 - t255 * t374;
t44 = t147 * t455 - t149 * t390 + t151 * t391 - t248 * t256 + t249 * t258 - t254 * t374;
t43 = t146 * t453 + t148 * t392 + t150 * t393 - t246 * t257 + t247 * t259 + t255 * t376;
t42 = t147 * t453 + t149 * t392 + t151 * t393 - t246 * t256 + t247 * t258 + t254 * t376;
t40 = t594 * t453 + t589 * t376 + (-t490 * t614 - t571 * t592) * t485;
t39 = -t594 * t455 + t589 * t374 + (t490 * t613 + t571 * t593) * t485;
t38 = t110 * t487 + (t491 * t78 - t494 * t79) * t485;
t37 = t110 * t604 + t453 * t79 + t455 * t78;
t36 = -t194 * t246 + t195 * t248 + t390 * t92 + t392 * t93;
t33 = -t110 * t499 + t390 * t78 - t392 * t79;
t32 = (t152 * t491 + t153 * t494 + (t261 * t494 + (-t260 + t555) * t491) * qJD(1)) * t485 + t500;
t31 = t613 * t487 + (qJD(1) * t504 + t491 * t518) * t485 + t559;
t30 = (t558 - t614) * t487 + (t494 * t518 + t575 * t589) * t485 + t553;
t29 = t99 * t487 + (t491 * t72 - t494 * t73) * t485;
t28 = t98 * t487 + (t491 * t70 - t494 * t71) * t485;
t27 = t453 * t73 + t455 * t72 + t604 * t99;
t26 = t453 * t71 + t455 * t70 + t604 * t98;
t25 = t390 * t72 - t392 * t73 - t499 * t99;
t24 = t390 * t70 - t392 * t71 - t499 * t98;
t23 = -t374 * t592 - t376 * t593 - t453 * t613 + t455 * t614;
t20 = t169 * t191 + t170 * t193 + t189 * t248 + t311 * t88 + t312 * t90 + t390 * t86;
t19 = t169 * t190 + t170 * t192 + t188 * t248 + t311 * t89 + t312 * t91 + t390 * t87;
t18 = t167 * t191 + t168 * t193 + t189 * t246 + t313 * t88 + t314 * t90 - t392 * t86;
t17 = t167 * t190 + t168 * t192 + t188 * t246 + t313 * t89 + t314 * t91 - t392 * t87;
t16 = (t613 * t494 + t614 * t491 + (t592 * t494 + (t555 - t593) * t491) * qJD(1)) * t485 + t500;
t15 = t76 + (t51 * t491 - t52 * t494 + (t125 * t494 + t126 * t491) * qJD(1)) * t485;
t14 = -t125 * t374 + t126 * t376 + t52 * t453 + t51 * t455 + t615;
t13 = t67 * t487 + (t44 * t491 - t45 * t494 + (t112 * t494 + t113 * t491) * qJD(1)) * t485;
t12 = t66 * t487 + (t42 * t491 - t43 * t494 + (t114 * t494 + t115 * t491) * qJD(1)) * t485;
t11 = -t112 * t374 + t113 * t376 + t44 * t455 + t45 * t453 + (t138 * t571 + t490 * t67) * t485;
t10 = -t114 * t374 + t115 * t376 + t42 * t455 + t43 * t453 + (t139 * t571 + t490 * t66) * t485;
t9 = t47 + (t21 * t491 - t22 * t494 + (t79 * t491 + t78 * t494) * qJD(1)) * t485;
t8 = t21 * t455 + t22 * t453 - t78 * t374 + t79 * t376 + t616;
t7 = t21 * t390 - t22 * t392 + t79 * t246 + t78 * t248 + t617;
t6 = t35 * t487 + (t19 * t491 - t20 * t494 + (t491 * t71 + t494 * t70) * qJD(1)) * t485;
t5 = t34 * t487 + (t17 * t491 - t18 * t494 + (t491 * t73 + t494 * t72) * qJD(1)) * t485;
t4 = t19 * t455 + t20 * t453 - t374 * t70 + t376 * t71 + (t35 * t490 + t571 * t98) * t485;
t3 = t17 * t455 + t18 * t453 - t374 * t72 + t376 * t73 + (t34 * t490 + t571 * t99) * t485;
t2 = t19 * t390 - t20 * t392 + t246 * t71 + t248 * t70 - t35 * t499 + t378 * t98;
t1 = t17 * t390 - t18 * t392 + t246 * t73 + t248 * t72 - t34 * t499 + t378 * t99;
t41 = [t77 - t430 * t602 + t48 - t485 * t600 + (t232 * t318 + t233 * t317) * t639 + (t165 * t253 + t166 * t252) * t638 + (t123 * t197 + t124 * t196) * t637 + (t171 * t97 + t172 * t96) * t636 + (t121 * t59 + t122 * t58) * t635 - t421 * t542 + t657; t76 + t47 + m(3) * (t213 * t317 + t214 * t318 + t232 * t299 + t233 * t300) + m(4) * (t130 * t252 + t131 * t253 + t165 * t217 + t166 * t218) + m(5) * (t123 * t136 + t124 * t137 + t196 * t83 + t197 * t84) + m(6) * (t103 * t96 + t104 * t97 + t171 * t60 + t172 * t61) + m(7) * (t121 * t30 + t122 * t31 + t58 * t68 + t59 * t69) + ((-t106 / 0.2e1 - t108 / 0.2e1 - t64 / 0.2e1 - t140 / 0.2e1 - t143 / 0.2e1 - t94 / 0.2e1 - t514) * t494 + (t107 / 0.2e1 + t109 / 0.2e1 + t63 / 0.2e1 + t141 / 0.2e1 + t142 / 0.2e1 + t95 / 0.2e1 + t513) * t491 + ((t229 / 0.2e1 + t226 / 0.2e1 + t159 / 0.2e1 + t201 / 0.2e1 + t202 / 0.2e1 + t132 / 0.2e1 - t510) * t494 + (t200 / 0.2e1 + t203 / 0.2e1 + t133 / 0.2e1 + t228 / 0.2e1 + t227 / 0.2e1 + t160 / 0.2e1 + t509) * t491) * qJD(1)) * t485 + t645; (t16 * t62 + t30 * t69 + t31 * t68) * t635 + (t103 * t61 + t104 * t60 + t32 * t85) * t636 + (t116 * t55 + t136 * t84 + t137 * t83) * t637 + (t130 * t218 + t131 * t217 + t175 * t74) * t638 + (t300 * t213 + t299 * t214 + (t353 * t491 + t354 * t494) * (t281 * t494 + t282 * t491 + (t353 * t494 - t354 * t491) * qJD(1)) * t483) * t639 + (t13 + t6) * t603 + (-t12 - t5) * t601 + (t57 + t29) * t545 + (t56 + t28) * t544 + (((-t120 + t653) * t494 + (t119 - t652) * t491) * t545 + ((-t118 - t651) * t494 + (t117 + t654) * t491) * t544 + (t651 * t575 + t654 * t574 + (t117 * qJD(1) - t207 * t401 - t209 * t402 - t291 * t328 - t293 * t329 + t373 * t647 + t374 * t642 + t454 * t649 - t455 * t643 - t544 * t648) * t494 + (t118 * qJD(1) + t208 * t401 + t210 * t402 + t290 * t328 + t292 * t329 + (t491 * t656 + t655 * t574 + t640) * t485 + t650 * t454 + t646 * t373 + t644 * t455 + t641 * t374) * t491) * t603 + (t653 * t575 + t652 * t574 + (-qJD(1) * t119 + t207 * t403 + t209 * t404 + t291 * t326 + t293 * t327 + (t575 * t648 + t640) * t485 - t649 * t452 + t647 * t375 + t643 * t453 + t642 * t376) * t494 + (-t208 * t403 - t210 * t404 - t290 * t326 - t292 * t327 - qJD(1) * t120 + (t494 * t656 - t655 * t575) * t485 - t650 * t452 + t646 * t375 - t644 * t453 + t641 * t376) * t491) * t601) * t485 + (t15 + t9 + (t142 + t141 + t95) * t603 + (-t143 - t140 - t94) * t601 + (t228 + t227 + t160) * t545 + (t229 + t226 + t159) * t544 + ((-t106 - t108 - t64) * t494 + (t107 + t109 + t63) * t491 + ((t132 + t201 + t202) * t494 + (t133 + t200 + t203) * t491) * qJD(1)) * t485 + t645) * t487; m(7) * (-t121 * t373 + t122 * t375 + t452 * t58 + t454 * t59) + m(6) * (-t171 * t373 + t172 * t375 + t452 * t96 + t454 * t97) + m(5) * (t123 * t452 + t124 * t454 - t196 * t373 + t197 * t375) + m(4) * (t165 * t452 + t166 * t454 - t252 * t373 + t253 * t375); m(7) * (t30 * t454 + t31 * t452 - t373 * t69 + t375 * t68 + (-t16 * t493 + t572 * t62) * t485) + m(6) * (t103 * t375 - t104 * t373 + t452 * t61 + t454 * t60 + (-t32 * t493 + t572 * t85) * t485) + m(5) * (t136 * t375 - t137 * t373 + t452 * t84 + t454 * t83 + (t116 * t572 - t493 * t55) * t485) + m(4) * (t130 * t454 + t131 * t452 + t217 * t375 - t218 * t373 + (t175 * t572 - t493 * t74) * t485); 0.4e1 * (m(4) / 0.2e1 + t537) * (-t373 * t454 + t375 * t452 - t532); m(7) * (-t121 * t374 + t122 * t376 + t453 * t58 + t455 * t59) + m(6) * (-t171 * t374 + t172 * t376 + t453 * t96 + t455 * t97) + m(5) * (t123 * t453 + t124 * t455 - t196 * t374 + t197 * t376); m(7) * (t30 * t455 + t31 * t453 - t374 * t69 + t376 * t68 + (t16 * t490 + t571 * t62) * t485) + m(6) * (t103 * t376 - t104 * t374 + t453 * t61 + t455 * t60 + (t32 * t490 + t571 * t85) * t485) + m(5) * (t136 * t376 - t137 * t374 + t453 * t84 + t455 * t83 + (t116 * t571 + t490 * t55) * t485); 0.2e1 * t537 * (-t373 * t455 - t374 * t454 + t375 * t453 + t376 * t452 + (t490 ^ 2 - t493 ^ 2) * t483 * qJD(2)); 0.4e1 * t537 * (-t374 * t455 + t376 * t453 + t532); m(7) * (t101 * t59 + t102 * t58 + t121 * t40 + t122 * t39) + m(6) * (t171 * t81 + t172 * t80 + t186 * t97 + t187 * t96) + t513 * t455 + t514 * t453 + t509 * t376 + t510 * t374 + t615 + t616; (t8 / 0.2e1 + t14 / 0.2e1) * t487 + (t6 / 0.2e1 + t13 / 0.2e1) * t455 + (t5 / 0.2e1 + t12 / 0.2e1) * t453 + (t29 / 0.2e1 + t57 / 0.2e1) * t376 + (-t28 / 0.2e1 - t56 / 0.2e1) * t374 + m(7) * (t101 * t30 + t102 * t31 + t16 * t82 + t23 * t62 + t39 * t68 + t40 * t69) + m(6) * (t103 * t80 + t104 * t81 + t161 * t32 + t186 * t60 + t187 * t61 + t65 * t85) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t494 + (t4 / 0.2e1 + t11 / 0.2e1) * t491 + (t9 / 0.2e1 + t15 / 0.2e1) * t490 + (t38 / 0.2e1 + t162 * t628 + (t125 * t491 - t126 * t494) * t485 / 0.2e1) * t571 + ((t26 / 0.2e1 + t53 / 0.2e1) * t494 + (t27 / 0.2e1 + t54 / 0.2e1) * t491) * qJD(1)) * t485; m(6) * (-t186 * t373 + t187 * t375 + t452 * t80 + t454 * t81 + (t161 * t572 - t493 * t65) * t485) + m(7) * (-t101 * t373 + t102 * t375 + t39 * t452 + t40 * t454 + (-t23 * t493 + t572 * t82) * t485); m(6) * (-t186 * t374 + t187 * t376 + t453 * t80 + t455 * t81 + (t161 * t571 + t490 * t65) * t485) + m(7) * (-t101 * t374 + t102 * t376 + t39 * t453 + t40 * t455 + (t23 * t490 + t571 * t82) * t485); (t14 + t8) * t604 + (t11 + t4) * t455 + (t10 + t3) * t453 + (t27 + t54) * t376 + (-t26 - t53) * t374 + (t125 * t455 + t126 * t453 + t162 * t604 + t37) * t542 + (t101 * t40 + t102 * t39 + t23 * t82) * t635 + (t161 * t65 + t186 * t81 + t187 * t80) * t636; m(7) * (t121 * t49 + t122 * t50 + t128 * t58 + t129 * t59) + t566 * t392 + t565 * t390 + t564 * t248 + t563 * t246 + t617; m(7) * (t111 * t16 + t128 * t31 + t129 * t30 + t36 * t62 + t49 * t69 + t50 * t68) + t29 * t634 + t5 * t630 + t38 * t632 + t9 * t629 + t7 * t628 + t28 * t633 + t6 * t631 + (t2 * t627 - t494 * t1 / 0.2e1 + (t494 * t24 / 0.2e1 + t25 * t627) * qJD(1)) * t485; m(7) * (t128 * t375 - t129 * t373 + t452 * t50 + t454 * t49 + (t111 * t572 - t36 * t493) * t485); m(7) * (t128 * t376 - t129 * t374 + t453 * t50 + t455 * t49 + (t111 * t571 + t36 * t490) * t485); t37 * t632 + t8 * t629 + m(7) * (t101 * t49 + t102 * t50 + t111 * t23 + t128 * t39 + t129 * t40 + t36 * t82) - t374 * t24 / 0.2e1 + t455 * t2 / 0.2e1 + t376 * t25 / 0.2e1 + t453 * t1 / 0.2e1 + t27 * t634 + t3 * t630 + t26 * t633 + t4 * t631 + (t33 * t571 / 0.2e1 + t490 * t7 / 0.2e1) * t485; t248 * t24 + t390 * t2 + t246 * t25 - t392 * t1 + t378 * t33 - t499 * t7 + (t111 * t36 + t128 * t50 + t129 * t49) * t635;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t41(1) t41(2) t41(4) t41(7) t41(11) t41(16); t41(2) t41(3) t41(5) t41(8) t41(12) t41(17); t41(4) t41(5) t41(6) t41(9) t41(13) t41(18); t41(7) t41(8) t41(9) t41(10) t41(14) t41(19); t41(11) t41(12) t41(13) t41(14) t41(15) t41(20); t41(16) t41(17) t41(18) t41(19) t41(20) t41(21);];
Mq  = res;