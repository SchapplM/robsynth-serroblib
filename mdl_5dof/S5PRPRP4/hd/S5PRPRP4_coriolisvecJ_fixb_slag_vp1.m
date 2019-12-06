% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:35:53
% DurationCPUTime: 34.00s
% Computational Cost: add. (14536->651), mult. (22053->956), div. (0->0), fcn. (21212->8), ass. (0->361)
t605 = Icges(6,4) + Icges(5,5);
t604 = Icges(5,6) - Icges(6,6);
t285 = qJ(2) + pkin(8);
t282 = cos(t285);
t287 = cos(pkin(7));
t291 = cos(qJ(4));
t472 = t287 * t291;
t286 = sin(pkin(7));
t289 = sin(qJ(4));
t475 = t286 * t289;
t247 = t282 * t475 + t472;
t473 = t287 * t289;
t474 = t286 * t291;
t248 = t282 * t474 - t473;
t281 = sin(t285);
t481 = t281 * t286;
t101 = Icges(6,4) * t248 + Icges(6,2) * t481 + Icges(6,6) * t247;
t99 = Icges(5,5) * t248 - Icges(5,6) * t247 + Icges(5,3) * t481;
t596 = t101 + t99;
t497 = Icges(5,4) * t248;
t103 = -Icges(5,2) * t247 + Icges(5,6) * t481 + t497;
t220 = Icges(6,5) * t248;
t97 = Icges(6,6) * t481 + Icges(6,3) * t247 + t220;
t617 = t103 - t97;
t249 = t282 * t473 - t474;
t480 = t281 * t287;
t250 = t282 * t472 + t475;
t496 = Icges(5,4) * t250;
t104 = -Icges(5,2) * t249 + Icges(5,6) * t480 + t496;
t221 = Icges(6,5) * t250;
t98 = Icges(6,6) * t480 + Icges(6,3) * t249 + t221;
t616 = t104 - t98;
t100 = Icges(5,5) * t250 - Icges(5,6) * t249 + Icges(5,3) * t480;
t102 = Icges(6,4) * t250 + Icges(6,2) * t480 + Icges(6,6) * t249;
t591 = t102 + t100;
t491 = Icges(6,5) * t247;
t105 = Icges(6,1) * t248 + Icges(6,4) * t481 + t491;
t222 = Icges(5,4) * t247;
t107 = Icges(5,1) * t248 + Icges(5,5) * t481 - t222;
t615 = t105 + t107;
t490 = Icges(6,5) * t249;
t106 = Icges(6,1) * t250 + Icges(6,4) * t480 + t490;
t223 = Icges(5,4) * t249;
t108 = Icges(5,1) * t250 + Icges(5,5) * t480 - t223;
t614 = t106 + t108;
t580 = Icges(5,2) + Icges(6,3);
t637 = Icges(6,2) + Icges(5,3);
t478 = t281 * t291;
t275 = Icges(6,5) * t478;
t479 = t281 * t289;
t485 = Icges(6,6) * t282;
t164 = Icges(6,3) * t479 + t275 - t485;
t494 = Icges(5,4) * t291;
t375 = -Icges(5,2) * t289 + t494;
t170 = -Icges(5,6) * t282 + t281 * t375;
t613 = t164 - t170;
t371 = Icges(5,5) * t291 - Icges(5,6) * t289;
t166 = -Icges(5,3) * t282 + t281 * t371;
t374 = Icges(6,4) * t291 + Icges(6,6) * t289;
t168 = -Icges(6,2) * t282 + t281 * t374;
t636 = t166 + t168;
t489 = Icges(6,5) * t289;
t380 = Icges(6,1) * t291 + t489;
t172 = -Icges(6,4) * t282 + t281 * t380;
t495 = Icges(5,4) * t289;
t381 = Icges(5,1) * t291 - t495;
t174 = -Icges(5,5) * t282 + t281 * t381;
t578 = t172 + t174;
t635 = (-t605 * t289 - t604 * t291) * t281;
t602 = -t616 * t247 + t614 * t248 + t591 * t481;
t601 = -t617 * t249 + t615 * t250 + t596 * t480;
t603 = -t617 * t247 + t615 * t248 + t596 * t481;
t564 = -t616 * t249 + t614 * t250 + t591 * t480;
t634 = t613 * t247 + t578 * t248 + t636 * t481;
t633 = t613 * t249 + t578 * t250 + t636 * t480;
t443 = qJD(4) * t291;
t430 = t282 * t443;
t445 = qJD(4) * t287;
t454 = qJD(2) * t286;
t150 = -t286 * t430 + (t281 * t454 + t445) * t289;
t449 = qJD(2) * t291;
t434 = t281 * t449;
t151 = -qJD(4) * t247 - t286 * t434;
t433 = t282 * t454;
t632 = t604 * t150 + t605 * t151 + t637 * t433;
t444 = qJD(4) * t289;
t450 = qJD(2) * t289;
t152 = -t286 * t444 - t287 * t430 + t450 * t480;
t153 = -qJD(4) * t249 - t287 * t434;
t452 = qJD(2) * t287;
t432 = t282 * t452;
t631 = t604 * t152 + t605 * t153 + t637 * t432;
t167 = Icges(5,3) * t281 + t282 * t371;
t169 = Icges(6,2) * t281 + t282 * t374;
t630 = t635 * qJD(4) + (t167 + t169) * qJD(2);
t629 = Icges(5,1) + Icges(6,1);
t628 = Icges(5,4) - Icges(6,5);
t488 = Icges(6,5) * t291;
t370 = Icges(6,3) * t289 + t488;
t165 = Icges(6,6) * t281 + t282 * t370;
t171 = Icges(5,6) * t281 + t282 * t375;
t558 = -t165 + t171;
t173 = Icges(6,4) * t281 + t282 * t380;
t175 = Icges(5,5) * t281 + t282 * t381;
t557 = -t173 - t175;
t626 = (-t580 * t291 + t489 - t495) * t281;
t625 = t602 * t287;
t624 = t601 * t286;
t623 = -t580 * t150 - t628 * t151 - t604 * t433;
t622 = -t580 * t152 - t628 * t153 - t604 * t432;
t621 = t628 * t150 + t629 * t151 + t605 * t433;
t620 = t628 * t152 + t629 * t153 + t605 * t432;
t619 = t558 * qJD(2) + t626 * qJD(4);
t231 = (-Icges(5,1) * t289 - t494) * t281;
t447 = qJD(4) * t281;
t618 = -(-Icges(6,1) * t289 + t488) * t447 - qJD(4) * t231 + t557 * qJD(2);
t455 = qJD(2) * t282;
t612 = -t630 * t281 - t455 * t636;
t611 = t631 * t281 + t591 * t455;
t610 = t632 * t281 + t596 * t455;
t290 = sin(qJ(2));
t292 = cos(qJ(2));
t585 = -Icges(3,5) * t290 - Icges(4,5) * t281 - Icges(3,6) * t292 - Icges(4,6) * t282;
t609 = t564 * t287 + t624;
t608 = t603 * t286 + t625;
t607 = t633 * t281;
t606 = t634 * t281;
t261 = t281 * t445 + t454;
t525 = t261 / 0.2e1;
t262 = t286 * t447 - t452;
t523 = t262 / 0.2e1;
t570 = t617 * t150 + t615 * t151 + t623 * t247 + t621 * t248 + t610 * t286;
t569 = t616 * t150 + t614 * t151 + t622 * t247 + t620 * t248 + t611 * t286;
t568 = t617 * t152 + t615 * t153 + t623 * t249 + t621 * t250 + t610 * t287;
t567 = t616 * t152 + t614 * t153 + t622 * t249 + t620 * t250 + t611 * t287;
t388 = t105 * t291 + t289 * t97;
t47 = -t101 * t282 + t281 * t388;
t368 = -t103 * t289 + t107 * t291;
t49 = t281 * t368 - t282 * t99;
t600 = t47 + t49;
t387 = t106 * t291 + t289 * t98;
t48 = -t102 * t282 + t281 * t387;
t367 = -t104 * t289 + t108 * t291;
t50 = -t100 * t282 + t281 * t367;
t599 = t48 + t50;
t365 = t164 * t289 + t172 * t291;
t482 = t168 * t282;
t66 = t281 * t365 - t482;
t364 = -t170 * t289 + t174 * t291;
t483 = t166 * t282;
t67 = t281 * t364 - t483;
t582 = t66 + t67;
t403 = rSges(5,1) * t291 - rSges(5,2) * t289;
t179 = -rSges(5,3) * t282 + t281 * t403;
t593 = t179 * t261;
t592 = t179 * t262;
t588 = (t613 * t152 - t578 * t153 + t619 * t249 + t618 * t250 + t612 * t287) * t282 + (t609 * t282 + t607) * qJD(2);
t587 = (t613 * t150 - t578 * t151 + t619 * t247 + t618 * t248 + t612 * t286) * t282 + (t608 * t282 + t606) * qJD(2);
t586 = t585 * qJD(2);
t501 = Icges(3,4) * t290;
t378 = -Icges(3,2) * t292 - t501;
t500 = Icges(3,4) * t292;
t384 = -Icges(3,1) * t290 - t500;
t498 = Icges(4,4) * t282;
t499 = Icges(4,4) * t281;
t544 = qJD(2) * t281;
t584 = -(-t290 * t378 + t292 * t384) * qJD(2) + (-Icges(4,2) * t282 - t499) * t544 - (-Icges(4,1) * t281 - t498) * t455;
t583 = t568 * t523 + t567 * t525 + t588 * qJD(4) / 0.2e1;
t283 = t286 ^ 2;
t284 = t287 ^ 2;
t456 = t283 + t284;
t446 = qJD(4) * t282;
t577 = t635 * t446 + (t605 * t247 + t604 * t248) * t262 + (t605 * t249 + t604 * t250) * t261;
t393 = t49 * t286 + t50 * t287;
t394 = t47 * t286 + t48 * t287;
t576 = t393 + t394;
t575 = t587 * qJD(4) + t569 * t261 + t570 * t262;
t541 = ((t368 + t388) * qJD(2) - t632) * t282 + (t621 * t291 + t623 * t289 + (-t289 * t615 - t291 * t617) * qJD(4) + t596 * qJD(2)) * t281;
t540 = ((t367 + t387) * qJD(2) - t631) * t282 + (t620 * t291 + t622 * t289 + (-t289 * t614 - t291 * t616) * qJD(4) + t591 * qJD(2)) * t281;
t573 = rSges(6,1) + pkin(4);
t572 = t602 * t261 + t603 * t262 - t446 * t634;
t571 = t564 * t261 + t601 * t262 - t446 * t633;
t566 = t599 * t261 + t600 * t262 - t582 * t446;
t563 = rSges(6,3) + qJ(5);
t315 = -t281 * t370 + t485;
t134 = t315 * t286;
t140 = t170 * t286;
t562 = t134 + t140;
t135 = t315 * t287;
t141 = t170 * t287;
t561 = t135 + t141;
t142 = t172 * t286;
t144 = t174 * t286;
t560 = -t142 - t144;
t143 = t172 * t287;
t145 = t174 * t287;
t559 = -t143 - t145;
t556 = t586 * t286;
t555 = t586 * t287;
t554 = t585 * t286;
t553 = t585 * t287;
t344 = t167 - t364;
t345 = -t169 + t365;
t528 = (t168 * t287 + t387) * t261 + (t168 * t286 + t388) * t262;
t529 = -(-t166 * t287 - t367) * t261 - (-t166 * t286 - t368) * t262;
t552 = (-t528 - t529 + (-t344 + t345) * t446) * t281;
t379 = -Icges(3,2) * t290 + t500;
t209 = Icges(3,6) * t286 + t287 * t379;
t385 = Icges(3,1) * t292 - t501;
t211 = Icges(3,5) * t286 + t287 * t385;
t377 = -Icges(4,2) * t281 + t498;
t383 = Icges(4,1) * t282 - t499;
t542 = -(Icges(4,6) * t286 + t287 * t377) * t282 - (Icges(4,5) * t286 + t287 * t383) * t281;
t551 = t584 * t287 + (t209 * t292 + t211 * t290 - t542) * qJD(2);
t208 = -Icges(3,6) * t287 + t286 * t379;
t210 = -Icges(3,5) * t287 + t286 * t385;
t543 = (-Icges(4,6) * t287 + t286 * t377) * t282 + (-Icges(4,5) * t287 + t286 * t383) * t281;
t550 = t584 * t286 + (t208 * t292 + t210 * t290 + t543) * qJD(2);
t401 = pkin(4) * t291 + qJ(5) * t289;
t402 = rSges(6,1) * t291 + rSges(6,3) * t289;
t549 = rSges(6,2) * t282 + (-t401 - t402) * t281;
t548 = t591 * t261 + t262 * t596;
t547 = -t482 - t483;
t535 = t582 * t544 + (t630 * t282 + (t618 * t291 + t619 * t289 + (t289 * t578 - t291 * t613) * qJD(4)) * t281 + ((-t364 - t365) * t282 - t636 * t281 + t576) * qJD(2)) * t282;
t534 = (t578 + t626) * t446 + (t248 * t580 + t222 - t491 - t615) * t262 + (t250 * t580 + t223 - t490 - t614) * t261;
t533 = (Icges(6,1) * t479 - t231 - t275 - t613) * t446 + (-t247 * t629 + t220 - t497 - t617) * t262 + (-t249 * t629 + t221 - t496 - t616) * t261;
t532 = t577 * t281;
t273 = rSges(3,1) * t290 + rSges(3,2) * t292;
t339 = qJD(2) * t273;
t440 = qJD(5) * t289;
t272 = t281 * t440;
t268 = pkin(3) * t282 + pkin(6) * t281;
t224 = t268 * t286;
t225 = t268 * t287;
t518 = pkin(2) * t292;
t182 = -qJ(3) * t287 + t286 * t518;
t183 = qJ(3) * t286 + t287 * t518;
t427 = t182 * t454 + t183 * t452 + qJD(1);
t386 = t224 * t454 + t225 * t452 + t427;
t470 = rSges(6,2) * t480 + t249 * t563 + t250 * t573;
t471 = rSges(6,2) * t481 + t247 * t563 + t248 * t573;
t27 = t261 * t471 - t262 * t470 + t272 + t386;
t267 = pkin(3) * t281 - pkin(6) * t282;
t342 = qJD(2) * t267;
t206 = t286 * t342;
t207 = t287 * t342;
t510 = pkin(2) * qJD(2);
t437 = t290 * t510;
t448 = qJD(3) * t287;
t263 = -t286 * t437 - t448;
t280 = qJD(3) * t286;
t264 = -t287 * t437 + t280;
t459 = t263 * t454 + t264 * t452;
t414 = -t206 * t454 - t207 * t452 + t459;
t431 = t281 * t443;
t441 = qJD(5) * t249;
t514 = rSges(6,2) * t432 - t152 * t563 + t153 * t573 + t441;
t442 = qJD(5) * t247;
t515 = rSges(6,2) * t433 - t150 * t563 + t151 * t573 + t442;
t5 = qJD(5) * t431 - t514 * t262 + t515 * t261 + (t440 + (-t286 * t470 + t287 * t471) * qJD(4)) * t455 + t414;
t530 = t27 * t514 + t470 * t5;
t527 = -t290 * (t378 * t286 + t210) - t292 * (-t384 * t286 + t208);
t293 = qJD(2) ^ 2;
t526 = -t261 / 0.2e1;
t524 = -t262 / 0.2e1;
t519 = pkin(2) * t290;
t180 = rSges(6,2) * t281 + t282 * t402;
t235 = (-rSges(6,1) * t289 + rSges(6,3) * t291) * t281;
t330 = t282 * t450 + t431;
t513 = t272 + t330 * qJ(5) + (-t281 * t444 + t282 * t449) * pkin(4) + qJD(2) * t180 + qJD(4) * t235;
t477 = t282 * t286;
t476 = t282 * t287;
t469 = -t247 * t573 + t248 * t563;
t468 = t249 * t573 - t250 * t563;
t467 = t549 * t286;
t466 = t549 * t287;
t465 = t286 * t182 + t287 * t183;
t463 = t401 * t282 + t180;
t462 = t456 * t342;
t458 = (-pkin(4) * t289 + qJ(5) * t291) * t281 + t235;
t457 = t286 * t263 + t287 * t264;
t439 = qJD(2) * qJD(4);
t438 = t293 * t518;
t436 = t292 * t510;
t424 = t452 / 0.2e1;
t422 = -t446 / 0.2e1;
t421 = t446 / 0.2e1;
t265 = rSges(4,1) * t281 + rSges(4,2) * t282;
t420 = -t265 - t519;
t266 = rSges(4,1) * t282 - rSges(4,2) * t281;
t419 = -t266 - t518;
t418 = -t267 - t519;
t417 = t439 / 0.2e1;
t416 = t456 * t290;
t415 = t286 * t224 + t287 * t225 + t465;
t413 = -t286 * t206 - t287 * t207 + t457;
t408 = -t179 + t418;
t407 = t281 * t417;
t406 = t282 * t417;
t253 = t266 * qJD(2);
t405 = -t253 - t436;
t254 = t268 * qJD(2);
t404 = -t254 - t436;
t274 = rSges(3,1) * t292 - rSges(3,2) * t290;
t112 = rSges(5,1) * t250 - rSges(5,2) * t249 + rSges(5,3) * t480;
t86 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t432;
t236 = (-rSges(5,1) * t289 - rSges(5,2) * t291) * t281;
t321 = rSges(5,3) * t281 + t282 * t403;
t96 = qJD(2) * t321 + qJD(4) * t236;
t34 = -t286 * t438 - t86 * t446 - t261 * t96 + (-t254 * t286 + (t112 * t281 - t179 * t476) * qJD(4)) * qJD(2);
t110 = rSges(5,1) * t248 - rSges(5,2) * t247 + rSges(5,3) * t481;
t84 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t433;
t35 = -t287 * t438 + t84 * t446 + t262 * t96 + (-t254 * t287 + (-t110 * t281 + t179 * t477) * qJD(4)) * qJD(2);
t400 = t286 * t35 - t287 * t34;
t390 = qJD(2) * t418;
t323 = t286 * t390 - t448;
t37 = t261 * t549 - t446 * t470 + t323 + t442;
t329 = t287 * t390 + t280;
t38 = -t262 * t549 + t446 * t471 + t329 + t441;
t399 = -t286 * t37 - t287 * t38;
t55 = -t112 * t446 + t323 - t593;
t56 = t110 * t446 + t329 + t592;
t392 = -t286 * t55 - t287 * t56;
t391 = qJD(2) * t420;
t366 = t110 * t287 - t112 * t286;
t363 = t456 * t274;
t362 = t456 * t339;
t361 = t404 - t96;
t360 = t418 + t549;
t359 = t286 * t406;
t358 = t287 * t406;
t357 = -qJD(2) * t253 - t438;
t356 = -qJD(2) * t254 - t438;
t343 = t404 - t513;
t190 = -rSges(4,3) * t287 + t266 * t286;
t191 = rSges(4,3) * t286 + t266 * t287;
t65 = (t190 * t286 + t191 * t287) * qJD(2) + t427;
t341 = t65 * t265;
t338 = qJD(2) * t265;
t36 = t110 * t261 - t112 * t262 + t386;
t337 = t36 * t366;
t328 = t27 * t515 + t471 * t5;
t327 = t37 * t470 - t38 * t471;
t326 = (t384 * t287 - t209) * t292 + (-t378 * t287 - t211) * t290;
t297 = (t27 * t471 + t37 * t549) * t287 + (-t27 * t470 - t38 * t549) * t286;
t296 = t286 * t542 + t287 * t543;
t205 = t287 * t338;
t204 = t286 * t338;
t157 = t357 * t287;
t156 = t357 * t286;
t155 = t287 * t391 + t280;
t154 = t286 * t391 - t448;
t132 = -rSges(5,1) * t249 - rSges(5,2) * t250;
t128 = -rSges(5,1) * t247 - rSges(5,2) * t248;
t113 = t362 * qJD(2);
t88 = qJD(2) * t363 + qJD(1);
t68 = (-t204 * t286 - t205 * t287) * qJD(2) + t459;
t26 = t282 * t366 * t439 + t261 * t84 - t262 * t86 + t414;
t11 = -qJD(5) * t152 + t356 * t287 + t513 * t262 + (t515 * t282 + (-t281 * t471 - t477 * t549) * qJD(2)) * qJD(4);
t10 = -qJD(5) * t150 + t356 * t286 - t513 * t261 + (-t514 * t282 + (t281 * t470 + t476 * t549) * qJD(2)) * qJD(4);
t1 = [-m(3) * t113 + m(4) * t68 + m(5) * t26 + m(6) * t5; (((t249 * t558 + t250 * t557 + t624) * t282 + t607) * qJD(4) + (((t547 + t564) * qJD(4) + t548) * t282 + t552) * t287 + (t249 * t562 + t250 * t560) * t262 + (t249 * t561 + t250 * t559) * t261) * t526 + (t286 * t567 - t287 * t568) * t525 + (((t247 * t558 + t248 * t557 + t625) * t282 + t606) * qJD(4) + (((t547 + t603) * qJD(4) + t548) * t282 + t552) * t286 + (t247 * t562 + t248 * t560) * t262 + (t247 * t561 + t248 * t559) * t261) * t524 + (t286 * t569 - t287 * t570) * t523 + t286 * t583 - t575 * t287 / 0.2e1 + (t550 * t284 + (t555 * t286 + (-t551 - t556) * t287) * t286) * t454 + (-t556 * t284 + (t551 * t286 + (-t550 + t555) * t287) * t286) * t452 - (t553 * qJD(2) * t283 + (-t527 * t287 + t296 + (t326 - t554) * t286) * t452) * t454 / 0.2e1 + ((t326 * t286 + t296 + (-t527 - t553) * t287) * t454 + t554 * qJD(2) * t284) * t424 - t566 * t447 / 0.2e1 + (((t135 * t289 - t143 * t291 + t102) * t261 + (t134 * t289 - t142 * t291 + t101) * t262 + t66 * qJD(4)) * t281 + ((-t345 * t282 + (-t165 * t289 - t173 * t291 - t168) * t281 + t394) * qJD(4) + t528) * t282 + ((t141 * t289 - t145 * t291 + t100) * t261 + (t140 * t289 - t144 * t291 + t99) * t262 + t67 * qJD(4)) * t281 + ((t344 * t282 + (t171 * t289 - t175 * t291 - t166) * t281 + t393) * qJD(4) + t529) * t282) * t421 + (t599 * t286 - t600 * t287) * t407 + (t602 * t286 - t287 * t603) * t359 + (t564 * t286 - t601 * t287) * t358 + (t5 * t415 + t27 * t413 + (t11 * t360 + t343 * t38 + t530) * t287 + (t10 * t360 + t343 * t37 + t328) * t286 + t27 * t462 - (t27 * t282 + t281 * t399) * t440 - (-t27 * t466 + t38 * t463) * t262 - (t27 * t467 - t37 * t463) * t261 - (t399 * t268 + (-t27 * t416 + t292 * t399) * pkin(2)) * qJD(2) - (t327 * t281 + (-t37 * t466 + t38 * t467 + t297) * t282) * qJD(4)) * m(6) + (t26 * t415 + t36 * t413 + (t26 * t112 + t35 * t408 + t36 * t86 + t361 * t56) * t287 + (t26 * t110 + t34 * t408 + t36 * t84 + t361 * t55) * t286 - t36 * (-t286 * t593 + t287 * t592 - t462) - (-t261 * t55 + t262 * t56) * t321 - (t392 * t268 + (t292 * t392 - t36 * t416) * pkin(2)) * qJD(2) - ((-t110 * t56 + t112 * t55) * t281 + t337 * t282) * qJD(4)) * m(5) + (t68 * t465 + t65 * t457 + (t155 * t405 + t157 * t420 + t68 * t191 - t65 * t205) * t287 + (t154 * t405 + t156 * t420 + t68 * t190 - t65 * t204) * t286 - (-t65 * pkin(2) * t416 + (t155 * t419 - t287 * t341) * t287 + (t154 * t419 - t286 * t341) * t286) * qJD(2)) * m(4) + (-t113 * t363 - t362 * t88 + (t273 * t274 * t293 + t339 * t88) * t456) * m(3) + ((-t541 + t571) * t287 + (t540 + t572) * t286) * t422; m(4) * (-t156 * t287 + t157 * t286) + m(5) * t400 + m(6) * (-t10 * t287 + t11 * t286); (t249 * t534 + t250 * t533 - t287 * t532) * t526 + ((t286 * t568 + t287 * t567) * t281 + t588) * t525 + (t247 * t534 + t248 * t533 - t286 * t532) * t524 + ((t286 * t570 + t287 * t569) * t281 + t587) * t523 - (qJD(4) * t535 + t261 * t540 + t262 * t541) * t282 / 0.2e1 + t575 * t481 / 0.2e1 + t480 * t583 + t566 * t544 / 0.2e1 + ((t286 * t541 + t287 * t540) * t281 + t535) * t422 + (t577 * t282 + (t289 * t534 + t291 * t533) * t281) * t421 + t572 * t433 / 0.2e1 + t571 * t282 * t424 + (t281 * t576 - t282 * t582) * t407 + (t608 * t281 - t282 * t634) * t359 + (t609 * t281 - t282 * t633) * t358 + ((qJD(2) * t297 - t10 * t470 + t11 * t471 - t37 * t514 + t38 * t515) * t282 + (t327 * qJD(2) + (t10 * t549 - t37 * t513 + t328) * t287 + (-t11 * t549 + t38 * t513 - t530) * t286) * t281 - (t248 * t37 + t250 * t38 + t27 * t478) * qJD(5) - (t27 * t468 + t38 * t458) * t262 - (t27 * t469 - t37 * t458) * t261 - (t37 * t468 + t38 * t469) * t446) * m(6) + ((t35 * t110 - t34 * t112 - t55 * t86 + t56 * t84 + (t337 + (t286 * t56 - t287 * t55) * t179) * qJD(2)) * t282 + (t56 * (-qJD(2) * t110 + t286 * t96) + t55 * (qJD(2) * t112 - t287 * t96) + t26 * t366 + t36 * (-t286 * t86 + t287 * t84) + t400 * t179) * t281 - t56 * (t128 * t446 + t236 * t262) - t55 * (-t132 * t446 - t236 * t261) - t36 * (t128 * t261 - t132 * t262)) * m(5); (t10 * t247 + t11 * t249 + t5 * t479 + (-t247 * t446 - t262 * t479 - t152) * t38 + (t249 * t446 + t261 * t479 - t150) * t37 + (-t247 * t261 + t249 * t262 + t330) * t27) * m(6);];
tauc = t1(:);
