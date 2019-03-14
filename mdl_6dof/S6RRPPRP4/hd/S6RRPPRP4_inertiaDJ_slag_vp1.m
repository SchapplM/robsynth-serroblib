% Calculate time derivative of joint inertia matrix for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:37:03
% EndTime: 2019-03-09 08:37:39
% DurationCPUTime: 23.42s
% Computational Cost: add. (23136->1085), mult. (63986->1445), div. (0->0), fcn. (68829->8), ass. (0->452)
t603 = Icges(5,2) + Icges(4,3);
t611 = t603 / 0.2e1;
t604 = Icges(5,4) + Icges(4,5);
t610 = -Icges(4,6) + Icges(5,6);
t365 = sin(qJ(5));
t366 = sin(qJ(2));
t363 = sin(pkin(9));
t551 = cos(qJ(5));
t462 = t551 * t363;
t364 = cos(pkin(9));
t523 = t364 * t366;
t300 = t365 * t523 - t366 * t462;
t401 = t363 * t365 + t364 * t551;
t301 = t401 * t366;
t368 = cos(qJ(2));
t199 = Icges(7,5) * t301 + Icges(7,6) * t368 + Icges(7,3) * t300;
t202 = Icges(6,4) * t301 - Icges(6,2) * t300 + Icges(6,6) * t368;
t609 = t199 - t202;
t200 = Icges(6,5) * t301 - Icges(6,6) * t300 + Icges(6,3) * t368;
t201 = Icges(7,4) * t301 + Icges(7,2) * t368 + Icges(7,6) * t300;
t608 = -t200 - t201;
t203 = Icges(7,1) * t301 + Icges(7,4) * t368 + Icges(7,5) * t300;
t204 = Icges(6,1) * t301 - Icges(6,4) * t300 + Icges(6,5) * t368;
t607 = t203 + t204;
t369 = cos(qJ(1));
t367 = sin(qJ(1));
t519 = t367 * t368;
t314 = t363 * t519 + t364 * t369;
t315 = -t363 * t369 + t364 * t519;
t242 = -t314 * t551 + t315 * t365;
t243 = t314 * t365 + t315 * t551;
t592 = rSges(7,3) + qJ(6);
t594 = rSges(7,1) + pkin(5);
t606 = t242 * t592 + t243 * t594;
t492 = qJD(2) * t367;
t455 = t366 * t492;
t494 = qJD(1) * t369;
t459 = t368 * t494;
t495 = qJD(1) * t367;
t262 = -t364 * t495 + (-t455 + t459) * t363;
t518 = t368 * t369;
t317 = t367 * t363 + t364 * t518;
t263 = qJD(1) * t317 - t364 * t455;
t132 = qJD(5) * t243 - t262 * t551 + t263 * t365;
t133 = -qJD(5) * t242 + t262 * t365 + t263 * t551;
t605 = -t242 * qJD(6) - t132 * t592 - t133 * t594;
t602 = Icges(5,6) / 0.2e1 - Icges(4,6) / 0.2e1;
t597 = t367 / 0.2e1;
t557 = -t369 / 0.2e1;
t593 = -qJD(1) / 0.2e1;
t521 = t366 * t367;
t149 = Icges(7,5) * t243 - Icges(7,6) * t521 + Icges(7,3) * t242;
t153 = Icges(7,4) * t243 - Icges(7,2) * t521 + Icges(7,6) * t242;
t157 = Icges(7,1) * t243 - Icges(7,4) * t521 + Icges(7,5) * t242;
t491 = qJD(2) * t368;
t456 = t364 * t491;
t457 = t363 * t491;
t225 = qJD(5) * t301 + t365 * t456 - t457 * t551;
t226 = (-t364 * t365 + t462) * t366 * qJD(5) + t401 * t491;
t493 = qJD(2) * t366;
t452 = t367 * t491;
t385 = t366 * t494 + t452;
t75 = Icges(7,5) * t133 - Icges(7,6) * t385 + Icges(7,3) * t132;
t79 = Icges(7,4) * t133 - Icges(7,2) * t385 + Icges(7,6) * t132;
t83 = Icges(7,1) * t133 - Icges(7,4) * t385 + Icges(7,5) * t132;
t20 = t149 * t225 - t153 * t493 + t157 * t226 + t300 * t75 + t301 * t83 + t368 * t79;
t151 = Icges(6,5) * t243 - Icges(6,6) * t242 - Icges(6,3) * t521;
t155 = Icges(6,4) * t243 - Icges(6,2) * t242 - Icges(6,6) * t521;
t159 = Icges(6,1) * t243 - Icges(6,4) * t242 - Icges(6,5) * t521;
t77 = Icges(6,5) * t133 - Icges(6,6) * t132 - Icges(6,3) * t385;
t81 = Icges(6,4) * t133 - Icges(6,2) * t132 - Icges(6,6) * t385;
t85 = Icges(6,1) * t133 - Icges(6,4) * t132 - Icges(6,5) * t385;
t22 = -t151 * t493 - t155 * t225 + t159 * t226 - t300 * t81 + t301 * t85 + t368 * t77;
t585 = -t20 - t22;
t316 = t363 * t518 - t367 * t364;
t244 = -t316 * t551 + t317 * t365;
t245 = t316 * t365 + t317 * t551;
t520 = t366 * t369;
t150 = Icges(7,5) * t245 - Icges(7,6) * t520 + Icges(7,3) * t244;
t154 = Icges(7,4) * t245 - Icges(7,2) * t520 + Icges(7,6) * t244;
t158 = Icges(7,1) * t245 - Icges(7,4) * t520 + Icges(7,5) * t244;
t490 = qJD(2) * t369;
t453 = t366 * t490;
t260 = qJD(1) * t314 + t363 * t453;
t261 = -qJD(1) * t315 - t364 * t453;
t130 = qJD(5) * t245 + t260 * t551 + t261 * t365;
t131 = -qJD(5) * t244 - t260 * t365 + t261 * t551;
t451 = t368 * t490;
t460 = t366 * t495;
t384 = -t451 + t460;
t74 = Icges(7,5) * t131 + Icges(7,6) * t384 + Icges(7,3) * t130;
t78 = Icges(7,4) * t131 + Icges(7,2) * t384 + Icges(7,6) * t130;
t82 = Icges(7,1) * t131 + Icges(7,4) * t384 + Icges(7,5) * t130;
t21 = t150 * t225 - t154 * t493 + t158 * t226 + t300 * t74 + t301 * t82 + t368 * t78;
t152 = Icges(6,5) * t245 - Icges(6,6) * t244 - Icges(6,3) * t520;
t156 = Icges(6,4) * t245 - Icges(6,2) * t244 - Icges(6,6) * t520;
t160 = Icges(6,1) * t245 - Icges(6,4) * t244 - Icges(6,5) * t520;
t76 = Icges(6,5) * t131 - Icges(6,6) * t130 + Icges(6,3) * t384;
t80 = Icges(6,4) * t131 - Icges(6,2) * t130 + Icges(6,6) * t384;
t84 = Icges(6,1) * t131 - Icges(6,4) * t130 + Icges(6,5) * t384;
t23 = -t152 * t493 - t156 * t225 + t160 * t226 - t300 * t80 + t301 * t84 + t368 * t76;
t584 = t21 + t23;
t414 = Icges(5,5) * t364 + Icges(5,3) * t363;
t418 = Icges(4,4) * t364 - Icges(4,2) * t363;
t591 = ((t414 - t418) * t368 + t610 * t366) * qJD(2);
t421 = Icges(5,1) * t364 + Icges(5,5) * t363;
t422 = Icges(4,1) * t364 - Icges(4,4) * t363;
t590 = ((t421 + t422) * t368 + t604 * t366) * qJD(2);
t280 = -Icges(5,6) * t368 + t366 * t414;
t283 = -Icges(4,6) * t368 + t366 * t418;
t589 = t280 - t283;
t284 = -Icges(5,4) * t368 + t366 * t421;
t285 = -Icges(4,5) * t368 + t366 * t422;
t588 = t284 + t285;
t141 = Icges(7,5) * t226 - Icges(7,6) * t493 + Icges(7,3) * t225;
t142 = Icges(6,5) * t226 - Icges(6,6) * t225 - Icges(6,3) * t493;
t143 = Icges(7,4) * t226 - Icges(7,2) * t493 + Icges(7,6) * t225;
t144 = Icges(6,4) * t226 - Icges(6,2) * t225 - Icges(6,6) * t493;
t145 = Icges(7,1) * t226 - Icges(7,4) * t493 + Icges(7,5) * t225;
t146 = Icges(6,1) * t226 - Icges(6,4) * t225 - Icges(6,5) * t493;
t587 = t608 * t493 + (t142 + t143) * t368 + (t145 + t146) * t301 + (t141 - t144) * t300 + t607 * t226 + t609 * t225;
t63 = t150 * t300 + t154 * t368 + t158 * t301;
t65 = t152 * t368 - t156 * t300 + t160 * t301;
t543 = t63 + t65;
t62 = t149 * t300 + t153 * t368 + t157 * t301;
t64 = t151 * t368 - t155 * t300 + t159 * t301;
t544 = t62 + t64;
t586 = t367 * t544 + t369 * t543;
t56 = -t151 * t521 - t155 * t242 + t159 * t243;
t57 = -t152 * t521 - t156 * t242 + t160 * t243;
t429 = t367 * t56 + t369 * t57;
t54 = t149 * t242 - t153 * t521 + t157 * t243;
t55 = t150 * t242 - t154 * t521 + t158 * t243;
t430 = t367 * t54 + t369 * t55;
t90 = t199 * t242 - t201 * t521 + t203 * t243;
t91 = -t200 * t521 - t202 * t242 + t204 * t243;
t583 = (t90 + t91) * t368 + (-t429 - t430) * t366;
t60 = -t151 * t520 - t244 * t155 + t245 * t159;
t61 = -t152 * t520 - t244 * t156 + t245 * t160;
t427 = t367 * t60 + t369 * t61;
t58 = t244 * t149 - t153 * t520 + t245 * t157;
t59 = t244 * t150 - t154 * t520 + t245 * t158;
t428 = t367 * t58 + t369 * t59;
t92 = t244 * t199 - t201 * t520 + t245 * t203;
t93 = -t200 * t520 - t244 * t202 + t245 * t204;
t582 = (t92 + t93) * t368 + (-t427 - t428) * t366;
t581 = rSges(7,2) * t460 + t244 * qJD(6) + t592 * t130 + t131 * t594;
t530 = Icges(3,4) * t368;
t420 = -Icges(3,2) * t366 + t530;
t292 = Icges(3,6) * t367 + t369 * t420;
t531 = Icges(3,4) * t366;
t424 = Icges(3,1) * t368 - t531;
t294 = Icges(3,5) * t367 + t369 * t424;
t408 = t292 * t366 - t294 * t368;
t394 = t408 * t367;
t291 = -Icges(3,6) * t369 + t367 * t420;
t293 = -Icges(3,5) * t369 + t367 * t424;
t409 = t291 * t366 - t293 * t368;
t395 = t409 * t369;
t516 = -rSges(7,2) * t521 + t606;
t580 = t516 * t369;
t578 = t592 * t244 + t245 * t594;
t209 = Icges(5,5) * t315 + Icges(5,6) * t521 + Icges(5,3) * t314;
t215 = Icges(4,4) * t315 - Icges(4,2) * t314 + Icges(4,6) * t521;
t577 = t209 - t215;
t210 = Icges(5,5) * t317 + Icges(5,6) * t520 + Icges(5,3) * t316;
t216 = Icges(4,4) * t317 - Icges(4,2) * t316 + Icges(4,6) * t520;
t576 = t210 - t216;
t217 = Icges(5,1) * t315 + Icges(5,4) * t521 + Icges(5,5) * t314;
t219 = Icges(4,1) * t315 - Icges(4,4) * t314 + Icges(4,5) * t521;
t575 = t217 + t219;
t218 = Icges(5,1) * t317 + Icges(5,4) * t520 + Icges(5,5) * t316;
t220 = Icges(4,1) * t317 - Icges(4,4) * t316 + Icges(4,5) * t520;
t574 = t218 + t220;
t573 = t363 * t610 + t604 * t364;
t572 = -rSges(3,2) * t520 + t367 * rSges(3,3);
t416 = Icges(3,5) * t368 - Icges(3,6) * t366;
t289 = -Icges(3,3) * t369 + t367 * t416;
t571 = t368 * t543 + t582;
t570 = 2 * m(3);
t569 = 2 * m(4);
t568 = 2 * m(5);
t567 = 2 * m(6);
t566 = 2 * m(7);
t361 = t367 ^ 2;
t362 = t369 ^ 2;
t565 = m(4) / 0.2e1;
t564 = m(5) / 0.2e1;
t563 = m(6) / 0.2e1;
t562 = m(7) / 0.2e1;
t561 = -pkin(3) - pkin(4);
t555 = -rSges(5,1) - pkin(3);
t553 = -rSges(7,2) - pkin(8);
t552 = -rSges(6,3) - pkin(8);
t550 = pkin(2) * t366;
t549 = pkin(2) * t368;
t358 = t367 * pkin(7);
t545 = t587 * t368;
t542 = -rSges(7,2) * t451 + t581;
t541 = -rSges(7,2) * t385 - t605;
t540 = -t300 * t609 - t301 * t607 + t368 * t608;
t539 = rSges(3,3) * t369;
t538 = rSges(5,3) * t314;
t537 = t262 * rSges(5,3);
t536 = -rSges(5,2) - qJ(3);
t535 = rSges(7,2) - qJ(3);
t534 = -rSges(4,3) - qJ(3);
t533 = rSges(6,3) - qJ(3);
t433 = -rSges(6,1) * t243 + rSges(6,2) * t242;
t162 = -rSges(6,3) * t521 - t433;
t526 = t162 * t369;
t525 = t363 * t366;
t524 = t363 * t368;
t522 = t364 * t368;
t517 = -rSges(7,2) * t493 + qJD(6) * t300 + t225 * t592 + t226 * t594;
t515 = -rSges(7,2) * t520 + t578;
t514 = rSges(7,2) * t368 + t300 * t592 + t301 * t594;
t512 = t245 * rSges(6,1) - t244 * rSges(6,2);
t249 = t317 * pkin(3) + t316 * qJ(4);
t321 = pkin(2) * t518 + qJ(3) * t520;
t511 = -t249 - t321;
t510 = -t262 * qJ(4) - t314 * qJD(4);
t509 = t261 * pkin(4) + pkin(8) * t460;
t432 = qJ(3) * t366 + t549;
t312 = qJD(2) * t432 - qJD(3) * t368;
t431 = pkin(3) * t364 + qJ(4) * t363;
t508 = -qJD(4) * t525 - t431 * t491 - t312;
t436 = rSges(4,1) * t364 - rSges(4,2) * t363;
t507 = -(rSges(4,3) * t366 + t368 * t436) * qJD(2) - t312;
t318 = t431 * t366;
t333 = -qJ(3) * t368 + t550;
t322 = t333 * t495;
t506 = t318 * t495 + t322;
t288 = -rSges(4,3) * t368 + t366 * t436;
t505 = -t288 - t333;
t320 = t432 * t367;
t504 = t367 * t320 + t369 * t321;
t503 = -t318 - t333;
t489 = qJD(3) * t366;
t502 = qJ(3) * t451 + t369 * t489;
t501 = rSges(3,2) * t460 + rSges(3,3) * t494;
t500 = t385 * pkin(8);
t499 = t369 * pkin(1) + t358;
t305 = t314 * qJ(4);
t359 = t369 * pkin(7);
t498 = t359 - t305;
t497 = t361 + t362;
t290 = Icges(3,3) * t367 + t369 * t416;
t496 = qJD(1) * t290;
t35 = t55 * t367 - t369 * t54;
t36 = t57 * t367 - t369 * t56;
t486 = -t35 / 0.2e1 - t36 / 0.2e1;
t37 = t59 * t367 - t369 * t58;
t38 = t61 * t367 - t369 * t60;
t485 = t37 / 0.2e1 + t38 / 0.2e1;
t211 = Icges(4,5) * t315 - Icges(4,6) * t314 + Icges(4,3) * t521;
t484 = t211 * t521;
t483 = t211 * t520;
t212 = Icges(4,5) * t317 - Icges(4,6) * t316 + Icges(4,3) * t520;
t482 = t212 * t521;
t481 = t212 * t520;
t213 = Icges(5,4) * t315 + Icges(5,2) * t521 + Icges(5,6) * t314;
t480 = t213 * t521;
t479 = t213 * t520;
t214 = Icges(5,4) * t317 + Icges(5,2) * t520 + Icges(5,6) * t316;
t478 = t214 * t521;
t477 = t214 * t520;
t475 = t131 * rSges(6,1) - t130 * rSges(6,2) + rSges(6,3) * t460;
t342 = pkin(2) * t455;
t386 = -t368 * t495 - t453;
t474 = t367 * (pkin(2) * t459 + qJ(3) * t385 + t367 * t489 - t342) + t369 * (pkin(2) * t386 - qJ(3) * t460 + t502) + t320 * t494;
t310 = t317 * pkin(4);
t271 = -pkin(8) * t520 + t310;
t473 = -t271 + t511;
t472 = t261 * rSges(4,1) + t260 * rSges(4,2) + rSges(4,3) * t451;
t471 = t261 * rSges(5,1) + rSges(5,2) * t451 - t260 * rSges(5,3);
t435 = rSges(5,1) * t364 + rSges(5,3) * t363;
t470 = -(rSges(5,2) * t366 + t368 * t435) * qJD(2) + t508;
t469 = -(pkin(4) * t522 - pkin(8) * t366) * qJD(2) + t508;
t327 = pkin(4) * t523 + pkin(8) * t368;
t468 = t327 * t495 + t506;
t287 = -rSges(5,2) * t368 + t366 * t435;
t467 = -t287 + t503;
t223 = t317 * rSges(5,1) + rSges(5,2) * t520 + t316 * rSges(5,3);
t224 = t317 * rSges(4,1) - t316 * rSges(4,2) + rSges(4,3) * t520;
t466 = -t327 + t503;
t356 = pkin(7) * t494;
t465 = t356 + t502;
t464 = t342 + t510;
t463 = -pkin(1) - t549;
t206 = rSges(6,1) * t301 - rSges(6,2) * t300 + rSges(6,3) * t368;
t461 = t206 * t495;
t454 = t366 * t491;
t450 = t280 / 0.2e1 - t283 / 0.2e1;
t449 = t284 / 0.2e1 + t285 / 0.2e1;
t448 = t367 * t514;
t447 = t369 * t514;
t247 = t505 * t369;
t446 = qJD(1) * t514;
t445 = pkin(2) * t453;
t444 = t564 + t563 + t562;
t148 = rSges(6,1) * t226 - rSges(6,2) * t225 - rSges(6,3) * t493;
t443 = -t148 + t469;
t442 = -t206 + t466;
t248 = pkin(3) * t315 + t305;
t441 = t367 * t248 + t369 * t249 + t504;
t440 = t499 + t321;
t194 = t467 * t369;
t439 = rSges(3,1) * t368 - rSges(3,2) * t366;
t334 = rSges(3,1) * t366 + rSges(3,2) * t368;
t438 = -t263 * rSges(4,1) + t262 * rSges(4,2);
t437 = -rSges(4,1) * t315 + rSges(4,2) * t314;
t434 = t133 * rSges(6,1) - t132 * rSges(6,2);
t426 = t469 - t517;
t425 = t466 - t514;
t423 = Icges(3,1) * t366 + t530;
t419 = Icges(3,2) * t368 + t531;
t164 = -rSges(6,3) * t520 + t512;
t413 = t164 * t367 - t526;
t412 = t367 * t162 + t164 * t369;
t299 = rSges(3,1) * t518 + t572;
t140 = t442 * t369;
t407 = -t368 * t544 - t583;
t406 = -pkin(1) - t439;
t405 = -t64 / 0.2e1 - t62 / 0.2e1 - t90 / 0.2e1 - t91 / 0.2e1;
t404 = -t92 / 0.2e1 - t93 / 0.2e1 - t65 / 0.2e1 - t63 / 0.2e1;
t398 = t261 * pkin(3) - t260 * qJ(4) + t316 * qJD(4);
t403 = t367 * (t263 * pkin(3) - t510) + t369 * t398 + t248 * t494 + t474;
t352 = pkin(8) * t521;
t270 = pkin(4) * t315 - t352;
t402 = t367 * t270 + t369 * t271 + t441;
t397 = qJD(2) * t334;
t114 = t425 * t369;
t396 = t315 * t561 + t352 + t498;
t393 = qJD(2) * t423;
t392 = qJD(2) * t419;
t391 = qJD(2) * (-Icges(3,5) * t366 - Icges(3,6) * t368);
t390 = t366 * t536 + t463;
t389 = t366 * t535 + t463;
t388 = t366 * t534 + t463;
t387 = t366 * t533 + t463;
t383 = t249 + t440;
t380 = t367 * t515 - t580;
t379 = t367 * t516 + t369 * t515;
t378 = t390 * t367;
t377 = t388 * t367;
t376 = t263 * t561 + t464 + t500;
t375 = t310 + t383;
t372 = t367 * (t263 * pkin(4) - t500) + t369 * (-pkin(8) * t451 + t509) + t270 * t494 + t403;
t371 = t398 + t465;
t370 = (-pkin(1) - t432) * t495 + t371 + t509;
t326 = t439 * qJD(2);
t298 = t367 * t439 - t539;
t268 = t299 + t499;
t267 = t367 * t406 + t359 + t539;
t246 = t505 * t367;
t229 = t367 * t391 + t496;
t228 = -qJD(1) * t289 + t369 * t391;
t222 = rSges(4,3) * t521 - t437;
t221 = rSges(5,1) * t315 + rSges(5,2) * t521 + t538;
t196 = t334 * t492 + ((-rSges(3,3) - pkin(7)) * t367 + t406 * t369) * qJD(1);
t195 = rSges(3,1) * t386 - rSges(3,2) * t451 - pkin(1) * t495 + t356 + t501;
t193 = t467 * t367;
t192 = t440 + t224;
t191 = t359 + t377 + t437;
t190 = t367 * t290 - t369 * t408;
t189 = t367 * t289 - t395;
t188 = -t290 * t369 - t394;
t187 = -t289 * t369 - t367 * t409;
t186 = Icges(4,1) * t263 - Icges(4,4) * t262 + Icges(4,5) * t385;
t185 = Icges(4,1) * t261 + Icges(4,4) * t260 - Icges(4,5) * t384;
t184 = Icges(5,1) * t263 + Icges(5,4) * t385 + Icges(5,5) * t262;
t183 = Icges(5,1) * t261 - Icges(5,4) * t384 - Icges(5,5) * t260;
t182 = Icges(4,4) * t263 - Icges(4,2) * t262 + Icges(4,6) * t385;
t181 = Icges(4,4) * t261 + Icges(4,2) * t260 - Icges(4,6) * t384;
t176 = Icges(5,5) * t263 + Icges(5,6) * t385 + Icges(5,3) * t262;
t175 = Icges(5,5) * t261 - Icges(5,6) * t384 - Icges(5,3) * t260;
t170 = qJD(1) * t247 + t367 * t507;
t169 = t288 * t495 + t369 * t507 + t322;
t139 = t442 * t367;
t136 = t383 + t223;
t135 = t315 * t555 + t378 + t498 - t538;
t119 = t367 * t222 + t224 * t369 + t504;
t118 = t342 + (t491 * t534 - t489) * t367 + (t369 * t388 - t358) * qJD(1) + t438;
t117 = qJD(1) * t377 - t445 + t465 + t472;
t116 = qJD(1) * t194 + t367 * t470;
t115 = t287 * t495 + t369 * t470 + t506;
t113 = t425 * t367;
t112 = t368 * t164 + t206 * t520;
t111 = -t162 * t368 - t206 * t521;
t110 = -t316 * t216 + t317 * t220 + t481;
t109 = -t316 * t215 + t317 * t219 + t483;
t108 = t316 * t210 + t317 * t218 + t477;
t107 = t316 * t209 + t317 * t217 + t479;
t106 = -t216 * t314 + t220 * t315 + t482;
t105 = -t215 * t314 + t219 * t315 + t484;
t104 = t210 * t314 + t218 * t315 + t478;
t103 = t209 * t314 + t217 * t315 + t480;
t102 = t520 * t552 + t375 + t512;
t101 = t367 * t387 + t396 + t433;
t97 = t413 * t366;
t96 = t367 * t221 + t223 * t369 + t441;
t95 = -t537 + t555 * t263 + (t491 * t536 - t489) * t367 + (t369 * t390 - t358) * qJD(1) + t464;
t94 = qJD(1) * t378 + t371 - t445 + t471;
t89 = -rSges(6,3) * t385 + t434;
t87 = -rSges(6,3) * t451 + t475;
t73 = t520 * t553 + t375 + t578;
t72 = t367 * t389 + t396 - t606;
t69 = qJD(1) * t140 + t367 * t443;
t68 = t369 * t443 + t461 + t468;
t67 = t366 * t447 + t368 * t515;
t66 = -t366 * t448 - t368 * t516;
t53 = t402 + t412;
t52 = t380 * t366;
t51 = t367 * (rSges(4,3) * t452 - t438) + t369 * t472 + (t369 * t222 + (-t224 - t321) * t367) * qJD(1) + t474;
t50 = (t491 * t533 - t489) * t367 + (t369 * t387 - t358) * qJD(1) + t376 - t434;
t49 = (t368 * t552 - t550) * t490 + t370 + t475;
t48 = qJD(1) * t114 + t367 * t426;
t47 = t369 * t426 + t495 * t514 + t468;
t46 = t379 + t402;
t45 = (-t206 * t492 - t89) * t368 + (qJD(2) * t162 - t367 * t148 - t206 * t494) * t366;
t44 = (t206 * t490 + t87) * t368 + (-qJD(2) * t164 + t148 * t369 - t461) * t366;
t41 = t367 * (t263 * rSges(5,1) + rSges(5,2) * t452 + t537) + t369 * t471 + (t369 * t221 + (-t223 + t511) * t367) * qJD(1) + t403;
t40 = (t491 * t535 - t489) * t367 + (t369 * t389 - t358) * qJD(1) + t376 + t605;
t39 = t370 + (t368 * t553 - t550) * t490 + t581;
t34 = -t132 * t202 + t133 * t204 - t142 * t521 - t242 * t144 + t243 * t146 - t200 * t385;
t33 = t132 * t199 + t133 * t203 + t242 * t141 - t143 * t521 + t243 * t145 - t201 * t385;
t32 = -t130 * t202 + t131 * t204 - t142 * t520 - t244 * t144 + t245 * t146 + t200 * t384;
t31 = t130 * t199 + t131 * t203 + t244 * t141 - t143 * t520 + t245 * t145 + t201 * t384;
t30 = t413 * t491 + (qJD(1) * t412 + t367 * t87 - t369 * t89) * t366;
t25 = (-qJD(2) * t448 - t541) * t368 + (qJD(2) * t516 - t367 * t517 - t369 * t446) * t366;
t24 = (qJD(2) * t447 + t542) * t368 + (-qJD(2) * t515 - t367 * t446 + t369 * t517) * t366;
t19 = t367 * t89 + t369 * t87 + (t526 + (-t164 + t473) * t367) * qJD(1) + t372;
t18 = -t132 * t156 + t133 * t160 - t152 * t385 - t242 * t80 + t243 * t84 - t521 * t76;
t17 = -t132 * t155 + t133 * t159 - t151 * t385 - t242 * t81 + t243 * t85 - t521 * t77;
t16 = t132 * t150 + t133 * t158 - t154 * t385 + t242 * t74 + t243 * t82 - t521 * t78;
t15 = t132 * t149 + t133 * t157 - t153 * t385 + t242 * t75 + t243 * t83 - t521 * t79;
t14 = -t130 * t156 + t131 * t160 + t152 * t384 - t244 * t80 + t245 * t84 - t520 * t76;
t13 = -t130 * t155 + t131 * t159 + t151 * t384 - t244 * t81 + t245 * t85 - t520 * t77;
t12 = t130 * t150 + t131 * t158 + t154 * t384 + t244 * t74 + t245 * t82 - t520 * t78;
t11 = t130 * t149 + t131 * t157 + t153 * t384 + t244 * t75 + t245 * t83 - t520 * t79;
t10 = t380 * t491 + (qJD(1) * t379 + t367 * t542 - t369 * t541) * t366;
t9 = t542 * t369 + t541 * t367 + (t580 + (t473 - t515) * t367) * qJD(1) + t372;
t8 = qJD(1) * t429 - t17 * t369 + t18 * t367;
t7 = qJD(1) * t430 - t15 * t369 + t16 * t367;
t6 = qJD(1) * t427 - t13 * t369 + t14 * t367;
t5 = qJD(1) * t428 - t11 * t369 + t12 * t367;
t4 = (-qJD(2) * t429 + t34) * t368 + (qJD(1) * t36 - qJD(2) * t91 - t17 * t367 - t18 * t369) * t366;
t3 = (-qJD(2) * t430 + t33) * t368 + (qJD(1) * t35 - qJD(2) * t90 - t15 * t367 - t16 * t369) * t366;
t2 = (-qJD(2) * t427 + t32) * t368 + (qJD(1) * t38 - qJD(2) * t93 - t13 * t367 - t14 * t369) * t366;
t1 = (-qJD(2) * t428 + t31) * t368 + (qJD(1) * t37 - qJD(2) * t92 - t11 * t367 - t12 * t369) * t366;
t26 = [(t195 * t268 + t196 * t267) * t570 + (t117 * t192 + t118 * t191) * t569 + (t135 * t95 + t136 * t94) * t568 + (t101 * t50 + t102 * t49) * t567 + (t39 * t73 + t40 * t72) * t566 + t591 * t525 + t590 * t523 + t589 * t457 + t588 * t456 + (t366 * t573 - t368 * t603 - t419 + t424) * t493 + (-t366 * t603 - t368 * t573 + t420 + t423) * t491 + t587; m(4) * (t117 * t246 + t118 * t247 + t169 * t191 + t170 * t192) + m(5) * (t115 * t135 + t116 * t136 + t193 * t94 + t194 * t95) + m(6) * (t101 * t68 + t102 * t69 + t139 * t49 + t140 * t50) + m(7) * (t113 * t39 + t114 * t40 + t47 * t72 + t48 * t73) + m(3) * ((-t195 * t367 - t196 * t369) * t334 + (-t267 * t369 - t268 * t367) * t326) + ((t292 * t593 + t392 * t597 + t604 * t263 / 0.2e1 + t385 * t611 + t602 * t262) * t369 + (t291 * t593 + t392 * t557 - t604 * t261 / 0.2e1 + t384 * t611 + t602 * t260) * t367) * t368 + ((t316 * t450 + t317 * t449 - t404) * t369 + (t314 * t450 + t315 * t449 - t405) * t367 + m(3) * (t267 * t367 - t268 * t369) * t334 + ((t292 / 0.2e1 - t214 / 0.2e1 - t212 / 0.2e1) * t369 + (t291 / 0.2e1 - t213 / 0.2e1 - t211 / 0.2e1) * t367) * t368 + (t363 * t576 + t364 * t574 + t294) * t520 / 0.2e1) * qJD(1) + ((t362 / 0.2e1 + t361 / 0.2e1) * t416 + t395 / 0.2e1 - t394 / 0.2e1) * qJD(2) + (t31 + t32 + (t574 * t522 + t576 * t524) * qJD(2) + (-t369 * t393 + (t183 + t185) * t364 + (t175 - t181) * t363 + (t212 + t214) * qJD(2) + (t363 * t577 + t364 * t575) * qJD(1)) * t366 + t590 * t317 + t591 * t316 + t588 * t261 - t589 * t260 + t584) * t597 + (t33 + t34 + (qJD(1) * t294 - t367 * t393 + (t184 + t186) * t364 + (t176 - t182) * t363) * t366 + (t577 * t524 + t575 * t522 + (t211 + t213) * t366) * qJD(2) + t590 * t315 + t591 * t314 + t588 * t263 + t589 * t262 - t585) * t557; -t369 * ((t369 * t229 + (t188 + t395) * qJD(1)) * t369 + (t187 * qJD(1) + (-t292 * t491 - t294 * t493 + t496) * t367 + (-t228 + (t291 * t368 + t293 * t366) * qJD(2) - t408 * qJD(1)) * t369) * t367) + t367 * ((t367 * t228 + (t189 + t394) * qJD(1)) * t367 + (t190 * qJD(1) + (t291 * t491 + t293 * t493) * t369 + (-t229 + (-t292 * t368 - t294 * t366) * qJD(2) + (t290 - t409) * qJD(1)) * t367) * t369) - t369 * t8 - t369 * t7 + t367 * t6 + t367 * t5 - t369 * ((-t314 * t176 - t315 * t184 - t262 * t209 - t263 * t217 + (t104 - t479) * qJD(1)) * t369 + (t314 * t175 + t315 * t183 + t262 * t210 + t263 * t218 + (t103 + t477) * qJD(1)) * t367) + t367 * ((t316 * t175 + t317 * t183 - t260 * t210 + t261 * t218 + (t107 - t478) * qJD(1)) * t367 + (-t316 * t176 - t317 * t184 + t260 * t209 - t261 * t217 + (t108 + t480) * qJD(1)) * t369) - t369 * ((t314 * t182 - t315 * t186 + t262 * t215 - t263 * t219 + (t106 - t483) * qJD(1)) * t369 + (-t314 * t181 + t315 * t185 - t262 * t216 + t263 * t220 + (t105 + t481) * qJD(1)) * t367) + t367 * ((-t316 * t181 + t317 * t185 + t260 * t216 + t261 * t220 + (t109 - t482) * qJD(1)) * t367 + (t316 * t182 - t317 * t186 - t260 * t215 - t261 * t219 + (t110 + t484) * qJD(1)) * t369) + (t113 * t48 + t114 * t47 + t46 * t9) * t566 + (t139 * t69 + t140 * t68 + t19 * t53) * t567 + (t115 * t194 + t116 * t193 + t41 * t96) * t568 + (t119 * t51 + t169 * t247 + t170 * t246) * t569 + ((t367 * t298 + t299 * t369) * ((qJD(1) * t298 - t369 * t397 + t501) * t369 + (-t367 * t397 + (-t299 + t572) * qJD(1)) * t367) + t497 * t334 * t326) * t570 + (t35 + t36 + (-t103 - t105 - t187) * t369 + (t104 + t106 + t188) * t367) * t495 + (t37 + t38 + (-t107 - t109 - t189) * t369 + (t108 + t110 + t190) * t367) * t494; 0.2e1 * ((t367 * t73 + t369 * t72) * t562 + (t101 * t369 + t102 * t367) * t563 + (t191 * t369 + t192 * t367) * t565 + (t135 * t369 + t136 * t367) * t564) * t491 + 0.2e1 * ((t367 * t39 + t369 * t40 + t494 * t73 - t495 * t72) * t562 + (-t101 * t495 + t102 * t494 + t367 * t49 + t369 * t50) * t563 + (t117 * t367 + t118 * t369 - t191 * t495 + t192 * t494) * t565 + (-t135 * t495 + t136 * t494 + t367 * t94 + t369 * t95) * t564) * t366; 0.2e1 * ((t113 * t492 + t114 * t490 - t9) * t562 + (t139 * t492 + t140 * t490 - t19) * t563 + (t193 * t492 + t194 * t490 - t41) * t564 + (t246 * t492 + t247 * t490 - t51) * t565) * t368 + 0.2e1 * ((qJD(2) * t46 + t113 * t494 - t114 * t495 + t367 * t48 + t369 * t47) * t562 + (qJD(2) * t53 + t139 * t494 - t140 * t495 + t367 * t69 + t369 * t68) * t563 + (qJD(2) * t96 + t115 * t369 + t116 * t367 + t193 * t494 - t194 * t495) * t564 + (qJD(2) * t119 + t169 * t369 + t170 * t367 + t246 * t494 - t247 * t495) * t565) * t366; 0.4e1 * (t565 + t444) * (-0.1e1 + t497) * t454; m(7) * (-t260 * t72 + t262 * t73 + t314 * t39 + t316 * t40) + m(6) * (-t101 * t260 + t102 * t262 + t314 * t49 + t316 * t50) + m(5) * (-t135 * t260 + t136 * t262 + t314 * t94 + t316 * t95); m(7) * (t113 * t262 - t114 * t260 + t314 * t48 + t316 * t47 + (t366 * t9 + t46 * t491) * t363) + m(6) * (t139 * t262 - t140 * t260 + t314 * t69 + t316 * t68 + (t19 * t366 + t491 * t53) * t363) + m(5) * (t115 * t316 + t116 * t314 + t193 * t262 - t194 * t260 + (t366 * t41 + t491 * t96) * t363); 0.2e1 * t444 * ((t367 * t314 + t369 * t316 - t524) * t491 + (t363 * t493 - t260 * t369 + t262 * t367 + (t369 * t314 - t367 * t316) * qJD(1)) * t366); 0.4e1 * t444 * (t363 ^ 2 * t454 - t260 * t316 + t262 * t314); m(6) * (t101 * t45 + t102 * t44 + t111 * t50 + t112 * t49) + m(7) * (t24 * t73 + t25 * t72 + t39 * t67 + t40 * t66) + (t367 * t405 + t369 * t404) * t491 + (t540 * qJD(2) + (-t21 / 0.2e1 - t31 / 0.2e1 - t32 / 0.2e1 - t23 / 0.2e1) * t369 + (-t22 / 0.2e1 - t20 / 0.2e1 - t33 / 0.2e1 - t34 / 0.2e1) * t367 + (-t367 * t404 + t369 * t405) * qJD(1)) * t366 + t545; m(6) * (t111 * t68 + t112 * t69 + t139 * t44 + t140 * t45 + t19 * t97 + t30 * t53) + m(7) * (t10 * t46 + t113 * t24 + t114 * t25 + t47 * t66 + t48 * t67 + t52 * t9) + (-t4 / 0.2e1 - t3 / 0.2e1 - t485 * t491) * t369 + (t2 / 0.2e1 + t1 / 0.2e1 + t486 * t491) * t367 + ((t367 * t485 + t369 * t486) * qJD(1) - (t7 + t8) * t367 / 0.2e1 + (t6 + t5) * t557 - (t367 * t543 - t369 * t544) * qJD(2) / 0.2e1) * t366 + (qJD(1) * t586 + t584 * t367 + t585 * t369) * t368 / 0.2e1 + (t583 * t367 + t582 * t369) * qJD(1) / 0.2e1; 0.2e1 * ((t111 * t490 + t112 * t492 - t30) * t563 + (t490 * t66 + t492 * t67 - t10) * t562) * t368 + 0.2e1 * ((qJD(2) * t97 - t111 * t495 + t112 * t494 + t367 * t44 + t369 * t45) * t563 + (qJD(2) * t52 + t24 * t367 + t25 * t369 + t494 * t67 - t495 * t66) * t562) * t366; m(6) * (-t111 * t260 + t112 * t262 + t314 * t44 + t316 * t45 + (t30 * t366 + t491 * t97) * t363) + m(7) * (t24 * t314 + t25 * t316 - t260 * t66 + t262 * t67 + (t10 * t366 + t491 * t52) * t363); (t10 * t52 + t24 * t67 + t25 * t66) * t566 + (t111 * t45 + t112 * t44 + t30 * t97) * t567 + ((t407 * t367 - t369 * t571) * qJD(2) + t545) * t368 + ((-t368 * t584 - t1 - t2) * t369 + (t368 * t585 - t3 - t4) * t367 + (t366 * t586 + 0.2e1 * t540 * t368) * qJD(2) + (t367 * t571 + t407 * t369) * qJD(1)) * t366; m(7) * (t130 * t72 + t132 * t73 + t242 * t39 + t244 * t40); m(7) * (t113 * t132 + t114 * t130 + t225 * t46 + t242 * t48 + t244 * t47 + t300 * t9); m(7) * ((-t225 + (t242 * t367 + t244 * t369) * qJD(2)) * t368 + (qJD(2) * t300 + t130 * t369 + t132 * t367 + (t242 * t369 - t244 * t367) * qJD(1)) * t366); m(7) * (t130 * t316 + t132 * t314 + t242 * t262 - t244 * t260 + (t225 * t366 + t300 * t491) * t363); m(7) * (t10 * t300 + t130 * t66 + t132 * t67 + t225 * t52 + t24 * t242 + t244 * t25); (t130 * t244 + t132 * t242 + t225 * t300) * t566;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t26(1) t26(2) t26(4) t26(7) t26(11) t26(16); t26(2) t26(3) t26(5) t26(8) t26(12) t26(17); t26(4) t26(5) t26(6) t26(9) t26(13) t26(18); t26(7) t26(8) t26(9) t26(10) t26(14) t26(19); t26(11) t26(12) t26(13) t26(14) t26(15) t26(20); t26(16) t26(17) t26(18) t26(19) t26(20) t26(21);];
Mq  = res;