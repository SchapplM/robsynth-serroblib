% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR9_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:44
% EndTime: 2019-12-31 18:24:04
% DurationCPUTime: 15.10s
% Computational Cost: add. (36068->622), mult. (43630->879), div. (0->0), fcn. (46828->8), ass. (0->349)
t390 = qJ(1) + pkin(8);
t387 = sin(t390);
t651 = -t387 / 0.2e1;
t388 = cos(t390);
t581 = t388 / 0.2e1;
t392 = sin(qJ(3));
t471 = m(5) / 0.4e1 + m(6) / 0.4e1;
t385 = t387 ^ 2;
t386 = t388 ^ 2;
t480 = t385 + t386;
t616 = t480 * t392;
t444 = 0.2e1 * t471 * t616;
t603 = m(6) / 0.2e1;
t604 = m(5) / 0.2e1;
t470 = t603 + t604;
t227 = -t392 * t470 + t444;
t395 = cos(qJ(3));
t354 = pkin(3) * t392 - qJ(4) * t395;
t391 = sin(qJ(5));
t394 = cos(qJ(5));
t440 = rSges(6,1) * t391 + rSges(6,2) * t394;
t608 = -rSges(6,3) * t392 + t395 * t440;
t445 = pkin(7) * t392 + t354 - t608;
t233 = t445 * t387;
t235 = t445 * t388;
t562 = rSges(5,2) * t392;
t439 = rSges(5,3) * t395 + t562;
t484 = t354 - t439;
t271 = t484 * t387;
t273 = t484 * t388;
t526 = t387 * t392;
t452 = rSges(5,1) * t388 - rSges(5,3) * t526;
t535 = qJ(4) * t392;
t454 = -pkin(2) - t535;
t455 = -sin(qJ(1)) * pkin(1) + t388 * pkin(6);
t209 = ((rSges(5,2) - pkin(3)) * t395 + t454) * t387 + t452 + t455;
t522 = t388 * t395;
t375 = rSges(5,2) * t522;
t554 = rSges(5,3) + qJ(4);
t567 = pkin(3) * t395;
t568 = cos(qJ(1)) * pkin(1);
t210 = t568 - t375 + (rSges(5,1) + pkin(6)) * t387 + (t392 * t554 + pkin(2) + t567) * t388;
t525 = t387 * t395;
t504 = t209 * t522 + t210 * t525;
t377 = pkin(7) * t525;
t453 = pkin(4) * t388 - t377;
t516 = t392 * t394;
t306 = t387 * t516 + t388 * t391;
t519 = t391 * t392;
t307 = -t387 * t519 + t388 * t394;
t614 = t307 * rSges(6,1) - t306 * rSges(6,2);
t169 = ((-rSges(6,3) - pkin(3)) * t395 + t454) * t387 + t453 + t455 + t614;
t304 = -t387 * t391 + t388 * t516;
t305 = t387 * t394 + t388 * t519;
t442 = t305 * rSges(6,1) + t304 * rSges(6,2);
t579 = rSges(6,3) + pkin(7);
t472 = pkin(3) + t579;
t171 = t568 + (pkin(4) + pkin(6)) * t387 + (t395 * t472 - t454) * t388 + t442;
t508 = t169 * t522 + t171 * t525;
t523 = t388 * t392;
t552 = (-t271 * t523 + t273 * t526 + t504) * t604 + (-t233 * t523 + t235 * t526 + t508) * t603;
t376 = pkin(3) * t526;
t211 = t376 + (t579 * t392 + (-qJ(4) - t440) * t395) * t387;
t362 = qJ(4) * t522;
t486 = t440 * t522;
t212 = -t472 * t523 + t362 + t486;
t251 = t376 + (-t395 * t554 - t562) * t387;
t446 = -pkin(3) * t523 + t362;
t481 = rSges(5,2) * t523 + rSges(5,3) * t522;
t252 = t446 + t481;
t553 = ((t251 * t388 + t252 * t387) * t392 + t504) * t604 + ((t211 * t388 + t212 * t387) * t392 + t508) * t603;
t625 = t552 - t553;
t650 = -t625 * qJD(1) - t227 * qJD(2);
t649 = Icges(5,1) + Icges(4,3);
t648 = Icges(5,4) - Icges(4,5);
t647 = Icges(5,5) - Icges(4,6);
t243 = rSges(6,1) * t304 - rSges(6,2) * t305;
t244 = rSges(6,1) * t306 + rSges(6,2) * t307;
t646 = m(6) * (-t169 * t244 + t171 * t243);
t543 = Icges(6,4) * t391;
t435 = Icges(6,2) * t394 + t543;
t405 = -Icges(6,6) * t392 + t395 * t435;
t542 = Icges(6,4) * t394;
t437 = Icges(6,1) * t391 + t542;
t406 = -Icges(6,5) * t392 + t395 * t437;
t433 = Icges(6,5) * t391 + Icges(6,6) * t394;
t404 = -Icges(6,3) * t392 + t395 * t433;
t512 = t395 * t404;
t166 = -t306 * t405 + t307 * t406 - t387 * t512;
t217 = Icges(6,5) * t305 + Icges(6,6) * t304 + Icges(6,3) * t522;
t544 = Icges(6,4) * t305;
t220 = Icges(6,2) * t304 + Icges(6,6) * t522 + t544;
t299 = Icges(6,4) * t304;
t223 = Icges(6,1) * t305 + Icges(6,5) * t522 + t299;
t119 = t217 * t525 + t220 * t306 - t223 * t307;
t219 = -Icges(6,5) * t307 + Icges(6,6) * t306 + Icges(6,3) * t525;
t301 = Icges(6,4) * t307;
t222 = Icges(6,2) * t306 + Icges(6,6) * t525 - t301;
t300 = Icges(6,4) * t306;
t224 = Icges(6,1) * t307 - Icges(6,5) * t525 - t300;
t120 = t219 * t525 + t222 * t306 + t224 * t307;
t430 = t119 * t388 + t120 * t387;
t52 = t166 * t392 + t395 * t430;
t565 = m(6) * qJD(5);
t164 = -t304 * t405 - t305 * t406 - t388 * t512;
t117 = t217 * t522 + t220 * t304 + t223 * t305;
t118 = t219 * t522 + t222 * t304 - t224 * t305;
t431 = t117 * t388 + t118 * t387;
t644 = t164 * t392 + t395 * t431;
t229 = rSges(6,3) * t525 - t614;
t532 = t229 * t392;
t187 = -t525 * t608 - t532;
t348 = Icges(4,5) * t395 - Icges(4,6) * t392;
t349 = -Icges(5,4) * t395 + Icges(5,5) * t392;
t643 = (t348 + t349) * t388 + t649 * t387;
t639 = t388 * t649 + t525 * t648 - t526 * t647;
t533 = t219 * t392;
t634 = t222 * t394 - t224 * t391;
t137 = t395 * t634 - t533;
t545 = Icges(4,4) * t392;
t353 = Icges(4,1) * t395 - t545;
t282 = Icges(4,5) * t387 + t353 * t388;
t538 = Icges(5,6) * t395;
t345 = Icges(5,3) * t392 - t538;
t283 = Icges(5,5) * t387 + t345 * t388;
t641 = -t282 * t525 - t283 * t526;
t640 = t117 * t387 - t118 * t388;
t638 = t392 * t648 + t395 * t647;
t637 = t388 * t643 + t641;
t626 = t282 * t522 + t283 * t523 + t387 * t643;
t369 = Icges(4,4) * t526;
t281 = Icges(4,1) * t525 - Icges(4,5) * t388 - t369;
t284 = Icges(5,5) * t388 + Icges(5,6) * t525 - Icges(5,3) * t526;
t636 = -t281 * t522 + t284 * t523 + t387 * t639;
t279 = Icges(4,4) * t525 - Icges(4,2) * t526 - Icges(4,6) * t388;
t364 = Icges(5,6) * t526;
t286 = Icges(5,4) * t388 + Icges(5,2) * t525 - t364;
t635 = t279 * t392 - t286 * t395;
t389 = Icges(4,4) * t395;
t540 = Icges(4,2) * t392;
t280 = Icges(4,6) * t387 + (t389 - t540) * t388;
t365 = Icges(5,6) * t523;
t285 = Icges(5,4) * t387 - Icges(5,2) * t522 + t365;
t628 = -t280 * t523 - t285 * t522 + t626;
t627 = t280 * t392 + t285 * t395 + t639;
t624 = -t279 * t523 - t280 * t526 - t285 * t525 + t286 * t522 - t636 - t637;
t585 = t387 / 0.2e1;
t583 = -t388 / 0.2e1;
t622 = t395 / 0.2e1;
t515 = t392 * t395;
t485 = t480 * t515;
t617 = t471 * (t485 - t515);
t613 = t638 * t387;
t612 = t638 * t388;
t315 = Icges(5,3) * t522 + t365;
t432 = Icges(5,2) * t392 + t538;
t317 = t432 * t388;
t350 = Icges(4,2) * t395 + t545;
t323 = t350 * t388;
t546 = Icges(4,1) * t392;
t438 = -t389 - t546;
t325 = t438 * t388;
t611 = (-t280 + t325 + t283 - t317) * t395 + (-t282 + t323 + t285 + t315) * t392;
t314 = Icges(5,3) * t525 + t364;
t316 = t432 * t387;
t322 = -Icges(4,2) * t525 - t369;
t324 = t438 * t387;
t610 = (t279 - t324 + t284 + t316) * t395 + (t281 + t322 + t286 - t314) * t392;
t336 = (Icges(6,2) * t391 - t542) * t395;
t337 = (-Icges(6,1) * t394 + t543) * t395;
t607 = -(t337 / 0.2e1 + t405 / 0.2e1) * t391 - (-t406 / 0.2e1 + t336 / 0.2e1) * t394;
t606 = 0.4e1 * qJD(1);
t605 = 2 * qJD(3);
t602 = t52 / 0.2e1;
t262 = t405 * t388;
t264 = t406 * t388;
t428 = -t220 * t394 - t223 * t391;
t420 = t388 * t404 - t428;
t98 = (-t262 * t394 - t264 * t391 + t217) * t395 + t420 * t392;
t601 = t98 / 0.2e1;
t261 = t405 * t387;
t263 = t406 * t387;
t419 = t387 * t404 + t634;
t99 = (-t261 * t394 - t263 * t391 + t219) * t395 + t419 * t392;
t600 = -t99 / 0.2e1;
t265 = t608 * t387;
t266 = -rSges(6,3) * t523 + t486;
t411 = rSges(6,3) * t522 + t442;
t407 = t392 * t411;
t122 = (t265 * t395 - t532) * t388 + (-t266 * t395 + t407) * t387;
t559 = rSges(6,3) * t395;
t329 = t392 * t440 + t559;
t518 = t392 * t608;
t145 = -t229 * t395 - t265 * t392 + (t329 * t395 + t518) * t387;
t146 = t395 * t442 + t392 * t266 + (-t518 + (-t329 + t559) * t395) * t388;
t152 = (t229 * t388 - t387 * t411) * t395;
t188 = -t522 * t608 - t407;
t506 = t187 * t522 - t188 * t525;
t597 = m(6) * (-t122 * t395 + (t145 * t388 + t146 * t387 + t152) * t392 + t506);
t596 = m(6) * (t145 * t169 + t146 * t171 + t187 * t211 - t188 * t212);
t595 = m(6) * (t122 * t152 + t145 * t187 - t146 * t188);
t338 = (-rSges(6,1) * t394 + rSges(6,2) * t391) * t395;
t592 = m(6) * (-t233 * t243 + t235 * t244 + (-t169 * t388 - t171 * t387) * t338);
t589 = m(6) * (t169 * t211 + t171 * t212);
t587 = m(6) * (t152 * t616 + t506);
t237 = Icges(6,5) * t304 - Icges(6,6) * t305;
t500 = -Icges(6,2) * t305 + t223 + t299;
t502 = -Icges(6,1) * t304 + t220 + t544;
t102 = t237 * t392 + (t391 * t502 - t394 * t500) * t395;
t586 = t102 / 0.2e1;
t584 = t387 / 0.4e1;
t582 = -t388 / 0.4e1;
t580 = t392 / 0.2e1;
t564 = rSges(4,1) * t395;
t456 = pkin(2) + t564;
t482 = rSges(4,2) * t526 + rSges(4,3) * t388;
t245 = -t387 * t456 + t455 + t482;
t373 = rSges(4,2) * t523;
t246 = t568 - t373 + t456 * t388 + (rSges(4,3) + pkin(6)) * t387;
t356 = rSges(4,1) * t392 + rSges(4,2) * t395;
t326 = t356 * t387;
t328 = t356 * t388;
t578 = m(4) * (t245 * t326 - t246 * t328);
t575 = m(5) * (t209 * t251 + t210 * t252);
t573 = m(5) * (-t209 * t526 + t210 * t523);
t572 = m(6) * (-t169 * t526 + t171 * t523);
t571 = m(6) * (-t187 * t526 - t188 * t523);
t182 = t243 * t388 + t244 * t387;
t570 = m(6) * (-t182 * t395 - t338 * t616);
t426 = t243 * t387 - t244 * t388;
t569 = m(6) * t426 * t392;
t558 = t387 * t52;
t557 = t388 * t644;
t534 = t217 * t392;
t527 = t404 * t392;
t312 = Icges(6,5) * t395 + t392 * t437;
t521 = t391 * t312;
t520 = t391 * t406;
t335 = (-Icges(6,5) * t394 + Icges(6,6) * t391) * t395;
t517 = t392 * t335;
t310 = Icges(6,6) * t395 + t392 * t435;
t514 = t394 * t310;
t513 = t394 * t405;
t93 = 0.2e1 * (t122 / 0.4e1 - t182 / 0.4e1) * m(6);
t510 = t93 * qJD(2);
t503 = -t233 * t525 - t235 * t522;
t501 = -Icges(6,1) * t306 + t222 - t301;
t499 = Icges(6,2) * t307 - t224 + t300;
t497 = -t271 * t525 - t273 * t522;
t492 = t387 * (qJ(4) * t525 - t376) + t388 * t446;
t358 = t535 + t567;
t491 = t480 * t358;
t489 = -t405 - t337;
t488 = -t406 + t336;
t487 = -t329 - t358;
t483 = rSges(5,2) * t395 - rSges(5,3) * t392 - t358;
t479 = qJD(1) * t392;
t478 = qJD(1) * t395;
t477 = qJD(3) * t392;
t476 = qJD(3) * t395;
t475 = qJD(5) * t395;
t467 = -t52 / 0.2e1 + t602;
t308 = Icges(6,3) * t395 + t392 * t433;
t423 = t513 + t520;
t418 = t308 - t423;
t400 = t395 * t418 + t527;
t128 = t304 * t310 + t305 * t312 + t388 * t400;
t402 = t395 * t420 - t534;
t81 = t262 * t304 + t264 * t305 + t388 * t402;
t401 = t395 * t419 - t533;
t82 = t261 * t304 + t263 * t305 + t388 * t401;
t12 = (t387 * t82 + t388 * t81 + t164) * t395 + (t128 - t431) * t392;
t88 = t237 * t522 + t304 * t500 - t305 * t502;
t238 = Icges(6,5) * t306 + Icges(6,6) * t307;
t89 = t238 * t522 + t304 * t499 - t305 * t501;
t40 = t387 * t88 - t388 * t89;
t466 = t40 / 0.2e1 - t12 / 0.2e1;
t127 = t306 * t310 - t307 * t312 + t387 * t400;
t79 = t262 * t306 - t264 * t307 + t387 * t402;
t80 = t261 * t306 - t263 * t307 + t387 * t401;
t11 = (t387 * t80 + t388 * t79 + t166) * t395 + (t127 - t430) * t392;
t90 = t237 * t525 + t306 * t500 + t307 * t502;
t91 = t238 * t525 + t306 * t499 + t307 * t501;
t41 = t387 * t90 - t388 * t91;
t465 = t41 / 0.2e1 - t11 / 0.2e1;
t462 = t525 / 0.4e1;
t457 = t349 / 0.2e1 + t348 / 0.2e1;
t115 = (t229 - t453) * t387 + (t387 * pkin(4) + t522 * t579 + t442) * t388 + t491;
t157 = -t387 * (rSges(5,2) * t525 + t452) + t388 * (t387 * rSges(5,1) + rSges(5,3) * t523 - t375) + t491;
t48 = m(5) * (t157 * t616 + t497) + m(6) * (t115 * t616 + t503);
t103 = t238 * t392 + (t391 * t501 - t394 * t499) * t395;
t139 = t304 * t488 - t305 * t489 + t335 * t522;
t140 = t306 * t488 + t307 * t489 + t335 * t525;
t443 = t592 / 0.2e1 + (t102 + t139) * t584 + (t103 + t140) * t582;
t136 = t395 * t428 + t534;
t429 = t136 * t388 - t137 * t387;
t410 = t644 * t582 - t52 * t584 + t558 / 0.4e1 + t557 / 0.4e1 + (t462 - t525 / 0.4e1) * t640;
t409 = (t639 * t388 + (t281 * t395 - t284 * t392 - t635) * t387) * t583 + (t388 * t627 - t626 + t628) * t581 + (t387 * t627 + t624 + t637) * t585;
t408 = t626 * t651 + t628 * t585 + ((t635 + t643) * t388 + t624 + t636 + t641) * t583;
t144 = (-t404 - t514 - t521) * t395 + t418 * t392;
t190 = t395 * t423 - t527;
t403 = t144 * t580 + t190 * t622 + t596 / 0.2e1 - (-t137 + t166) * t526 / 0.4e1 + (t127 + t99) * t462 - (t136 + t164) * t523 / 0.4e1 + (t128 + t98) * t522 / 0.4e1;
t399 = -t521 / 0.2e1 - t514 / 0.2e1 - t404 / 0.2e1 + t389 + t546 / 0.2e1 - t540 / 0.2e1 + t432 / 0.2e1 - t345 / 0.2e1;
t398 = -t520 / 0.2e1 - t513 / 0.2e1 + t308 / 0.2e1 + t353 / 0.2e1 - t350 / 0.2e1 + Icges(5,2) * t622 - Icges(5,6) * t392 - Icges(5,3) * t395 / 0.2e1;
t360 = -rSges(4,2) * t392 + t564;
t274 = t483 * t388;
t272 = t483 * t387;
t236 = (-pkin(7) * t395 + t487) * t388;
t234 = t387 * t487 - t377;
t232 = -t326 * t387 - t328 * t388;
t226 = t444 + (m(5) + m(6)) * t580;
t203 = 0.4e1 * t617;
t198 = t243 * t392 - t338 * t522;
t197 = -t244 * t392 + t338 * t525;
t181 = t385 * t439 + t388 * t481 + t492;
t168 = t426 * t395;
t163 = t569 / 0.2e1;
t156 = (t517 + (t391 * t489 - t394 * t488) * t395) * t392;
t147 = -pkin(7) * t616 + t265 * t387 + t266 * t388 + t492;
t141 = t570 / 0.2e1;
t130 = t571 / 0.2e1;
t94 = (t122 + t182) * t603;
t86 = t587 / 0.2e1;
t70 = t572 + t573;
t67 = t517 / 0.2e1 + t646 + t607 * t395;
t56 = t190 * t392 + t395 * t429;
t55 = t115 * t182 + (t233 * t387 + t235 * t388) * t338;
t39 = t387 * t81 - t388 * t82;
t38 = t387 * t79 - t388 * t80;
t35 = t130 - t569 / 0.2e1;
t34 = t163 + t130;
t33 = t163 - t571 / 0.2e1;
t30 = t597 / 0.2e1;
t29 = t140 * t392 + (t387 * t91 + t388 * t90) * t395;
t28 = t139 * t392 + (t387 * t89 + t388 * t88) * t395;
t27 = t392 * t398 + t395 * t399 + t575 + t578 + t589;
t26 = (t99 * t387 + t98 * t388 + t190) * t395 + (t144 - t429) * t392;
t19 = t86 + t30 - t570 / 0.2e1;
t18 = t141 + t86 - t597 / 0.2e1;
t17 = t141 + t30 - t587 / 0.2e1;
t9 = t552 + t553;
t7 = m(6) * t55 + t40 * t585 + t41 * t583;
t6 = t467 * t522;
t5 = t595 + (t12 * t581 + t11 * t585 + t56 / 0.2e1) * t395 + (-t557 / 0.2e1 - t558 / 0.2e1 + t26 / 0.2e1) * t392;
t4 = t387 * t409 + t388 * t408;
t3 = t403 + (-t139 / 0.4e1 - t102 / 0.4e1) * t387 - t592 / 0.2e1 + t410 + (t103 / 0.4e1 + t140 / 0.4e1) * t388;
t2 = (-t190 / 0.2e1 + (-t98 / 0.4e1 - t128 / 0.4e1) * t388 + (-t99 / 0.4e1 - t127 / 0.4e1) * t387) * t395 - t596 / 0.2e1 + (-t144 / 0.2e1 + (t164 / 0.4e1 + t136 / 0.4e1) * t388 + (t166 / 0.4e1 - t137 / 0.4e1) * t387) * t392 + t410 + t443;
t1 = t403 + t443;
t8 = [qJD(3) * t27 + qJD(4) * t70 + qJD(5) * t67, 0, t27 * qJD(1) + t9 * qJD(4) + t1 * qJD(5) + ((t209 * t274 + t210 * t272 - t251 * t273 - t252 * t271) * t604 + (t169 * t236 + t171 * t234 - t211 * t235 - t212 * t233) * t603) * t605 + ((-t286 / 0.2e1 + t314 / 0.2e1 - t281 / 0.2e1 - t322 / 0.2e1) * t476 + (t284 / 0.2e1 + t316 / 0.2e1 + t279 / 0.2e1 - t324 / 0.2e1) * t477) * t388 + ((t282 / 0.2e1 - t323 / 0.2e1 - t285 / 0.2e1 - t315 / 0.2e1) * t476 + (-t280 / 0.2e1 + t325 / 0.2e1 + t283 / 0.2e1 - t317 / 0.2e1) * t477) * t387 + ((-t127 / 0.2e1 + t600 + m(4) * (-t245 * t360 - t326 * t356) + t457 * t388 - t408) * t388 + (t601 + t128 / 0.2e1 + m(4) * (-t246 * t360 + t328 * t356) + t457 * t387 - t409) * t387) * qJD(3), qJD(1) * t70 + qJD(3) * t9 + qJD(5) * t34, t67 * qJD(1) + t1 * qJD(3) + t34 * qJD(4) + t156 * qJD(5) + (t169 * t197 + t171 * t198 - t187 * t244 - t188 * t243) * t565 + ((t139 / 0.2e1 + t586 - t467) * t388 + (t103 / 0.2e1 + t140 / 0.2e1) * t387) * t475; 0, 0, (m(4) * t232 / 0.2e1 + t181 * t604 + t147 * t603) * t605 + t226 * qJD(4) + t94 * qJD(5), t226 * qJD(3), qJD(3) * t94 - t168 * t565; t4 * qJD(3) + t625 * qJD(4) + t2 * qJD(5) + (-t589 / 0.4e1 - t575 / 0.4e1 - t578 / 0.4e1) * t606 - t399 * t478 - t398 * t479, qJD(4) * t227 - qJD(5) * t93, t4 * qJD(1) + t48 * qJD(4) + t7 * qJD(5) + (m(6) * (t115 * t147 - t233 * t234 - t235 * t236) + m(5) * (t157 * t181 - t271 * t272 - t273 * t274) + m(4) * (t356 * t360 * t480 + (t387 * (rSges(4,1) * t525 - t482) + t388 * (rSges(4,1) * t522 + t387 * rSges(4,3) - t373)) * t232) + (t39 + t612 * t385 + (t610 * t388 + (t611 - t613) * t387) * t388) * t585 + (t38 + t613 * t386 + (t611 * t387 + (t610 - t612) * t388) * t387) * t583) * qJD(3), t48 * qJD(3) + (-0.4e1 * t617 + 0.2e1 * t470 * (-t395 * t616 + t485)) * qJD(4) + t18 * qJD(5) - t650, t2 * qJD(1) - t510 + t7 * qJD(3) + t18 * qJD(4) + (-t56 / 0.2e1 + t466 * t388 + t465 * t387) * t475 + (t29 * t583 + t28 * t585 + m(6) * (-t115 * t168 + t152 * t182 - t197 * t235 - t198 * t233 + (-t187 * t388 + t188 * t387) * t338) - t595 + (-t26 / 0.2e1 + (-t103 / 0.2e1 + t644 / 0.2e1) * t388 + (t586 + t602) * t387) * t392) * qJD(5); -t625 * qJD(3) + t33 * qJD(5) + (-t572 / 0.4e1 - t573 / 0.4e1) * t606, -t227 * qJD(3), (m(6) * (-t147 * t395 + t503) + m(5) * (-t181 * t395 + t497) + 0.2e1 * ((t234 * t387 + t236 * t388 + t115) * t603 + (t272 * t387 + t274 * t388 + t157) * t604) * t392 - t48) * qJD(3) + t203 * qJD(4) + t17 * qJD(5) + t650, t203 * qJD(3), t33 * qJD(1) + t17 * qJD(3) + (t168 * t395 + (t197 * t388 + t198 * t387) * t392) * t565; -t335 * t479 / 0.2e1 + t3 * qJD(3) + t35 * qJD(4) + t6 * qJD(5) - qJD(1) * t646 - t607 * t478, qJD(3) * t93, t3 * qJD(1) + t510 + (((t39 / 0.2e1 + t137 / 0.2e1) * t395 + (-t640 / 0.2e1 + t600) * t392 + t465) * t388 + ((t38 / 0.2e1 + t136 / 0.2e1) * t395 + (t119 * t651 + t120 * t581 + t601) * t392 - t466) * t387 + (t115 * t122 - t145 * t235 - t146 * t233 + t147 * t152 + t187 * t236 - t188 * t234 - t55) * m(6)) * qJD(3) + t19 * qJD(4) + t5 * qJD(5), qJD(1) * t35 + qJD(3) * t19, t6 * qJD(1) + t5 * qJD(3) + (m(6) * (-t152 * t168 + t187 * t197 - t188 * t198) + t156 * t580 + (t28 * t581 + t29 * t585 + (t102 * t388 + t103 * t387) * t580) * t395) * qJD(5);];
Cq = t8;
