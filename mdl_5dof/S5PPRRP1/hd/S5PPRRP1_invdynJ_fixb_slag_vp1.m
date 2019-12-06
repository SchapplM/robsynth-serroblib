% Calculate vector of inverse dynamics joint torques for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:39
% EndTime: 2019-12-05 15:07:13
% DurationCPUTime: 29.71s
% Computational Cost: add. (15243->658), mult. (21661->906), div. (0->0), fcn. (20781->6), ass. (0->337)
t583 = Icges(5,5) + Icges(6,5);
t581 = Icges(5,6) + Icges(6,6);
t565 = Icges(5,1) + Icges(6,1);
t601 = Icges(5,2) + Icges(6,2);
t607 = Icges(5,3) + Icges(6,3);
t274 = pkin(8) + qJ(3);
t272 = sin(t274);
t278 = sin(qJ(4));
t279 = cos(qJ(4));
t606 = (-t583 * t278 - t581 * t279) * t272;
t273 = cos(t274);
t276 = cos(pkin(7));
t438 = t276 * t278;
t275 = sin(pkin(7));
t439 = t275 * t279;
t235 = t273 * t439 - t438;
t415 = qJD(3) * t272;
t392 = t278 * t415;
t150 = -qJD(4) * t235 + t275 * t392;
t437 = t276 * t279;
t440 = t275 * t278;
t234 = -t273 * t440 - t437;
t391 = t279 * t415;
t151 = qJD(4) * t234 - t275 * t391;
t413 = qJD(3) * t275;
t390 = t273 * t413;
t605 = t581 * t150 + t583 * t151 + t607 * t390;
t237 = t273 * t437 + t440;
t152 = -qJD(4) * t237 + t276 * t392;
t533 = t273 * t438 - t439;
t153 = -qJD(4) * t533 - t276 * t391;
t412 = qJD(3) * t276;
t389 = t273 * t412;
t604 = t581 * t152 + t583 * t153 + t607 * t389;
t347 = Icges(6,5) * t279 - Icges(6,6) * t278;
t161 = Icges(6,3) * t272 + t273 * t347;
t348 = Icges(5,5) * t279 - Icges(5,6) * t278;
t163 = Icges(5,3) * t272 + t273 * t348;
t603 = -t606 * qJD(4) + (-t161 - t163) * qJD(3);
t449 = t272 * t275;
t100 = Icges(5,5) * t235 + Icges(5,6) * t234 + Icges(5,3) * t449;
t98 = Icges(6,5) * t235 + Icges(6,6) * t234 + Icges(6,3) * t449;
t574 = t100 + t98;
t448 = t272 * t276;
t101 = Icges(5,5) * t237 - Icges(5,6) * t533 + Icges(5,3) * t448;
t99 = Icges(6,5) * t237 - Icges(6,6) * t533 + Icges(6,3) * t448;
t573 = t101 + t99;
t602 = Icges(5,4) + Icges(6,4);
t160 = -Icges(6,3) * t273 + t272 * t347;
t162 = -Icges(5,3) * t273 + t272 * t348;
t587 = t160 + t162;
t456 = Icges(6,4) * t279;
t351 = -Icges(6,2) * t278 + t456;
t165 = Icges(6,6) * t272 + t273 * t351;
t460 = Icges(5,4) * t279;
t352 = -Icges(5,2) * t278 + t460;
t167 = Icges(5,6) * t272 + t273 * t352;
t600 = t165 + t167;
t457 = Icges(6,4) * t278;
t355 = Icges(6,1) * t279 - t457;
t169 = Icges(6,5) * t272 + t273 * t355;
t461 = Icges(5,4) * t278;
t356 = Icges(5,1) * t279 - t461;
t171 = Icges(5,5) * t272 + t273 * t356;
t538 = -t169 - t171;
t599 = (-t601 * t279 - t457 - t461) * t272;
t598 = (t565 * t278 + t456 + t460) * t272;
t597 = t601 * t150 + t602 * t151 + t581 * t390;
t596 = t601 * t152 + t602 * t153 + t581 * t389;
t595 = t602 * t150 + t565 * t151 + t583 * t390;
t594 = t602 * t152 + t565 * t153 + t583 * t389;
t593 = t600 * qJD(3) + t599 * qJD(4);
t592 = t538 * qJD(3) + t598 * qJD(4);
t459 = Icges(6,4) * t235;
t102 = Icges(6,2) * t234 + Icges(6,6) * t449 + t459;
t463 = Icges(5,4) * t235;
t104 = Icges(5,2) * t234 + Icges(5,6) * t449 + t463;
t591 = t102 + t104;
t458 = Icges(6,4) * t237;
t103 = -Icges(6,2) * t533 + Icges(6,6) * t448 + t458;
t462 = Icges(5,4) * t237;
t105 = -Icges(5,2) * t533 + Icges(5,6) * t448 + t462;
t590 = t103 + t105;
t212 = Icges(6,4) * t234;
t106 = Icges(6,1) * t235 + Icges(6,5) * t449 + t212;
t214 = Icges(5,4) * t234;
t108 = Icges(5,1) * t235 + Icges(5,5) * t449 + t214;
t589 = t106 + t108;
t213 = Icges(6,4) * t533;
t107 = Icges(6,1) * t237 + Icges(6,5) * t448 - t213;
t215 = Icges(5,4) * t533;
t109 = Icges(5,1) * t237 + Icges(5,5) * t448 - t215;
t588 = t107 + t109;
t164 = -Icges(6,6) * t273 + t272 * t351;
t166 = -Icges(5,6) * t273 + t272 * t352;
t564 = t164 + t166;
t168 = -Icges(6,5) * t273 + t272 * t355;
t170 = -Icges(5,5) * t273 + t272 * t356;
t568 = t168 + t170;
t414 = qJD(3) * t273;
t586 = t603 * t272 - t587 * t414;
t585 = t604 * t272 + t573 * t414;
t584 = t605 * t272 + t574 * t414;
t501 = t276 ^ 2;
t502 = t275 ^ 2;
t526 = t501 + t502;
t555 = rSges(6,1) + pkin(4);
t349 = -Icges(4,5) * t272 - Icges(4,6) * t273;
t552 = t591 * t150 + t589 * t151 + t597 * t234 + t595 * t235 + t584 * t275;
t551 = t590 * t150 + t588 * t151 + t596 * t234 + t594 * t235 + t585 * t275;
t550 = t591 * t152 + t589 * t153 + t595 * t237 + t584 * t276 - t597 * t533;
t549 = t590 * t152 + t588 * t153 + t594 * t237 + t585 * t276 - t596 * t533;
t580 = -t564 * t150 - t568 * t151 - t593 * t234 + t592 * t235 + t586 * t275;
t579 = -t564 * t152 - t568 * t153 + t592 * t237 + t586 * t276 + t593 * t533;
t567 = t591 * t234 + t589 * t235 + t574 * t449;
t578 = t590 * t234 + t588 * t235 + t573 * t449;
t577 = t589 * t237 + t574 * t448 - t591 * t533;
t546 = t588 * t237 + t573 * t448 - t590 * t533;
t345 = -t102 * t278 + t106 * t279;
t47 = t272 * t345 - t273 * t98;
t343 = -t104 * t278 + t108 * t279;
t49 = -t273 * t100 + t272 * t343;
t545 = t47 + t49;
t344 = -t103 * t278 + t107 * t279;
t48 = t272 * t344 - t273 * t99;
t342 = -t105 * t278 + t109 * t279;
t50 = -t273 * t101 + t272 * t342;
t544 = t48 + t50;
t576 = t564 * t234 + t568 * t235 + t587 * t449;
t575 = t568 * t237 + t587 * t448 - t564 * t533;
t340 = -t164 * t278 + t168 * t279;
t447 = t273 * t160;
t62 = t272 * t340 - t447;
t339 = -t166 * t278 + t170 * t279;
t446 = t273 * t162;
t63 = t272 * t339 - t446;
t519 = t62 + t63;
t411 = qJD(4) * t272;
t387 = t276 * t411;
t241 = t387 + t413;
t388 = t275 * t411;
t242 = t388 - t412;
t410 = qJD(4) * t273;
t513 = t241 * (-t237 * t601 - t213 - t215 + t588) + t242 * (-t235 * t601 + t212 + t214 + t589) - t410 * (t568 + t599);
t566 = t275 * t276;
t563 = t606 * t410 + (-t583 * t234 + t581 * t235) * t242 + (t581 * t237 + t583 * t533) * t241;
t559 = 0.2e1 * qJD(3);
t558 = 2 * qJDD(3);
t407 = qJD(3) * qJD(4);
t308 = qJDD(4) * t272 + t273 * t407;
t404 = qJDD(3) * t275;
t156 = t276 * t308 + t404;
t403 = qJDD(3) * t276;
t157 = t275 * t308 - t403;
t240 = -qJDD(4) * t273 + t272 * t407;
t557 = t578 * t156 + t567 * t157 + t576 * t240 + t551 * t241 + t242 * t552 + t580 * t410;
t556 = t156 * t546 + t157 * t577 + t240 * t575 + t241 * t549 + t242 * t550 + t410 * t579;
t522 = ((t343 + t345) * qJD(3) - t605) * t273 + (t595 * t279 - t597 * t278 + (-t278 * t589 - t279 * t591) * qJD(4) + t574 * qJD(3)) * t272;
t521 = ((t342 + t344) * qJD(3) - t604) * t273 + (t594 * t279 - t596 * t278 + (-t278 * t588 - t279 * t590) * qJD(4) + t573 * qJD(3)) * t272;
t554 = t241 * t578 + t242 * t567 - t410 * t576;
t553 = t241 * t546 + t242 * t577 - t410 * t575;
t548 = t241 * t544 + t242 * t545 - t410 * t519;
t138 = t164 * t275;
t140 = t166 * t275;
t543 = -t138 - t140;
t139 = t164 * t276;
t141 = t166 * t276;
t542 = -t139 - t141;
t142 = t168 * t275;
t144 = t170 * t275;
t541 = -t142 - t144;
t143 = t168 * t276;
t145 = t170 * t276;
t540 = -t143 - t145;
t514 = t273 * pkin(3) + t272 * pkin(6);
t239 = t514 * qJD(3);
t408 = qJD(5) * t273;
t537 = -t239 + t408;
t318 = t163 - t339;
t319 = t161 - t340;
t503 = -(-t160 * t276 - t344) * t241 - (-t160 * t275 - t345) * t242;
t504 = -(-t162 * t276 - t342) * t241 - (-t162 * t275 - t343) * t242;
t536 = (-t503 - t504 + (-t318 - t319) * t410) * t272;
t535 = t241 * t573 + t242 * t574;
t269 = pkin(4) * t279 + pkin(3);
t488 = pkin(3) - t269;
t382 = t488 * t272;
t277 = -qJ(5) - pkin(6);
t487 = pkin(6) + t277;
t534 = -t273 * t487 + t382;
t245 = pkin(3) * t272 - pkin(6) * t273;
t317 = qJD(3) * t245;
t204 = t275 * t317;
t205 = t276 * t317;
t445 = t273 * t275;
t259 = pkin(6) * t445;
t444 = t273 * t276;
t260 = pkin(6) * t444;
t532 = -t275 * t204 - t276 * t205 - (-pkin(3) * t449 + t259) * t413 - (-pkin(3) * t448 + t260) * t412;
t531 = -t446 - t447;
t530 = t575 * t272;
t529 = t576 * t272;
t528 = t578 * t276;
t527 = t577 * t275;
t520 = ((-t339 - t340) * qJD(3) - t603) * t273 + (t592 * t279 + t593 * t278 + (t278 * t568 + t564 * t279) * qJD(4) - t587 * qJD(3)) * t272;
t264 = t272 * rSges(6,3);
t441 = t273 * t279;
t442 = t273 * t278;
t515 = rSges(6,1) * t441 - rSges(6,2) * t442 + t273 * t269 + t264;
t512 = (t564 + t598) * t410 + (t234 * t565 - t459 - t463 - t591) * t242 + (-t533 * t565 - t458 - t462 - t590) * t241;
t511 = t563 * t272;
t359 = t49 * t275 + t50 * t276;
t360 = t47 * t275 + t48 * t276;
t510 = t359 + t360;
t509 = t276 * t546 + t527;
t508 = t275 * t567 + t528;
t409 = qJD(5) * t272;
t253 = t276 * t409;
t271 = qJD(2) * t275;
t369 = -t245 * t412 + t271;
t482 = rSges(6,1) * t279;
t367 = -rSges(6,2) * t278 + t482;
t428 = -t273 * rSges(6,3) + t272 * t367 - t534;
t302 = -t272 * t487 - t273 * t488;
t468 = rSges(6,1) * t235 + rSges(6,2) * t234 + rSges(6,3) * t449 - pkin(4) * t438 + t275 * t302;
t36 = t242 * t428 + t410 * t468 + t253 + t369;
t252 = t275 * t409;
t416 = qJD(2) * t276;
t315 = -t245 * t413 - t416;
t467 = rSges(6,1) * t237 - rSges(6,2) * t533 + rSges(6,3) * t448 + pkin(4) * t440 + t276 * t302;
t37 = -t241 * t428 - t410 * t467 + t252 + t315;
t507 = t275 * t37 + t276 * t36;
t506 = g(1) * t276 + g(2) * t275;
t224 = t514 * t275;
t225 = t514 * t276;
t384 = t224 * t413 + t225 * t412 + qJD(1);
t33 = t241 * t468 - t242 * t467 + t384 - t408;
t480 = pkin(4) * qJD(4);
t401 = t278 * t480;
t284 = qJD(3) * t534 - t273 * t401;
t400 = t279 * t480;
t485 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t389 + t275 * t400 + t276 * t284 + t253;
t333 = -t204 * t413 - t205 * t412 + t224 * t404 + t225 * t403 + qJDD(1);
t486 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t390 + t275 * t284 - t276 * t400 + t252;
t5 = qJD(3) * t409 - qJDD(5) * t273 + t156 * t468 - t157 * t467 + t241 * t486 - t242 * t485 + t333;
t505 = t33 * t485 + t467 * t5;
t500 = t156 / 0.2e1;
t499 = t157 / 0.2e1;
t498 = t240 / 0.2e1;
t497 = -t241 / 0.2e1;
t496 = t241 / 0.2e1;
t495 = -t242 / 0.2e1;
t494 = t242 / 0.2e1;
t481 = rSges(6,2) * t279;
t232 = (-rSges(6,1) * t278 - t481) * t272;
t484 = qJD(4) * t232 - t272 * t401 - t408 + (t273 * t367 + t264 + t302) * qJD(3);
t483 = rSges(5,1) * t279;
t265 = t272 * rSges(5,3);
t233 = (-rSges(5,1) * t278 - rSges(5,2) * t279) * t272;
t368 = -rSges(5,2) * t278 + t483;
t97 = qJD(4) * t233 + (t273 * t368 + t265) * qJD(3);
t466 = -t239 - t97;
t443 = t273 * t277;
t316 = t382 - t443;
t397 = t272 * t439;
t398 = t272 * t440;
t420 = rSges(6,2) * t398 + rSges(6,3) * t445;
t432 = -rSges(6,1) * t397 + t275 * t316 - t259 + t420;
t395 = t272 * t437;
t396 = t272 * t438;
t418 = rSges(6,2) * t396 + rSges(6,3) * t444;
t431 = rSges(6,1) * t395 - t276 * t316 + t260 - t418;
t430 = -t235 * rSges(6,2) + t234 * t555;
t429 = t237 * rSges(6,2) + t555 * t533;
t175 = -t273 * rSges(5,3) + t272 * t368;
t423 = -t175 - t245;
t421 = t275 * t224 + t276 * t225;
t419 = rSges(5,2) * t398 + rSges(5,3) * t445;
t417 = rSges(5,2) * t396 + rSges(5,3) * t444;
t406 = qJDD(2) * t276;
t405 = qJDD(3) * t245;
t402 = -t239 - t484;
t399 = -m(3) - m(4) - m(5) - m(6);
t393 = -t245 - t428;
t380 = t412 / 0.2e1;
t378 = -t410 / 0.2e1;
t377 = t410 / 0.2e1;
t371 = t272 * t278 * pkin(4) - t232;
t244 = rSges(4,1) * t273 - rSges(4,2) * t272;
t243 = rSges(4,1) * t272 + rSges(4,2) * t273;
t111 = rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t449;
t270 = qJDD(2) * t275;
t330 = -qJD(3) * t239 - t405;
t83 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t390;
t34 = -t111 * t240 + t157 * t175 + t242 * t97 + t276 * t330 + t410 * t83 + t270;
t113 = rSges(5,1) * t237 - rSges(5,2) * t533 + rSges(5,3) * t448;
t85 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t389;
t35 = t113 * t240 - t156 * t175 - t241 * t97 + t275 * t330 - t410 * t85 - t406;
t365 = t275 * t34 - t276 * t35;
t350 = Icges(4,5) * t273 - Icges(4,6) * t272;
t341 = t111 * t276 - t113 * t275;
t336 = -(-t243 * t412 + t271) * t276 - (-t243 * t413 - t416) * t275;
t335 = t526 * t244;
t334 = t526 * qJD(3) * t243;
t177 = rSges(5,1) * t441 - rSges(5,2) * t442 + t265;
t238 = t244 * qJD(3);
t331 = -qJD(3) * t238 - qJDD(3) * t243;
t46 = t111 * t241 - t113 * t242 + t384;
t313 = t46 * t341;
t305 = qJD(3) * t349;
t304 = t33 * t486 + t468 * t5;
t303 = -t36 * t468 + t37 * t467;
t291 = qJD(3) * t537 + qJDD(5) * t272 - t405;
t283 = (t33 * t468 - t37 * t428) * t276 + (-t33 * t467 + t36 * t428) * t275;
t282 = t526 * t349;
t219 = t243 * t276;
t218 = t243 * t275;
t207 = t349 * t276;
t206 = t349 * t275;
t197 = t276 * t305;
t196 = t275 * t305;
t179 = Icges(4,3) * t275 + t276 * t350;
t178 = -Icges(4,3) * t276 + t275 * t350;
t149 = -rSges(5,1) * t395 + t417;
t147 = -rSges(5,1) * t397 + t419;
t133 = -rSges(5,1) * t533 - rSges(5,2) * t237;
t131 = rSges(5,1) * t234 - rSges(5,2) * t235;
t115 = t275 * t331 - t406;
t114 = t276 * t331 + t270;
t86 = qJD(3) * t335 + qJD(1);
t61 = -qJD(3) * t334 + qJDD(3) * t335 + qJDD(1);
t60 = -t113 * t410 - t175 * t241 + t315;
t59 = t111 * t410 + t175 * t242 + t369;
t26 = t111 * t156 - t113 * t157 + t241 * t83 - t242 * t85 + t333;
t7 = -t156 * t428 + t240 * t467 - t241 * t484 + t275 * t291 - t410 * t485 - t406;
t6 = t157 * t428 - t240 * t468 + t242 * t484 + t276 * t291 + t410 * t486 + t270;
t1 = [(m(2) + m(3)) * qJDD(1) + m(4) * t61 + m(5) * t26 + m(6) * t5 + (-m(2) + t399) * g(3); t399 * (g(1) * t275 - g(2) * t276) + m(4) * (t114 * t275 - t115 * t276) + m(5) * t365 + m(6) * (t275 * t6 - t276 * t7) + m(3) * t526 * qJDD(2); (t206 * qJD(3) * t501 + (-t276 * t207 + t282) * t413) * t380 - (t207 * qJD(3) * t502 + (-t275 * t206 + t282) * t412) * t413 / 0.2e1 + (t546 * t275 - t276 * t577) * t500 + (t275 * t578 - t276 * t567) * t499 + (t275 * t544 - t276 * t545) * t498 + (((t237 * t538 + t533 * t600 + t527) * t273 + t530) * qJD(4) + (((t531 + t546) * qJD(4) + t535) * t273 + t536) * t276 + (t237 * t541 - t533 * t543) * t242 + (t237 * t540 - t533 * t542) * t241) * t497 + (t275 * t549 - t276 * t550) * t496 + (((-t234 * t600 + t235 * t538 + t528) * t273 + t529) * qJD(4) + (((t531 + t567) * qJD(4) + t535) * t273 + t536) * t275 + (t234 * t543 + t235 * t541) * t242 + (t234 * t542 + t235 * t540) * t241) * t495 + (t275 * t551 - t276 * t552) * t494 - t548 * t411 / 0.2e1 + (((t141 * t278 - t145 * t279 + t101) * t241 + (t140 * t278 - t144 * t279 + t100) * t242 + t63 * qJD(4)) * t272 + ((t318 * t273 + (t278 * t167 - t279 * t171 - t162) * t272 + t359) * qJD(4) + t504) * t273 + ((t139 * t278 - t143 * t279 + t99) * t241 + (t138 * t278 - t142 * t279 + t98) * t242 + t62 * qJD(4)) * t272 + ((t319 * t273 + (t278 * t165 - t279 * t169 - t160) * t272 + t360) * qJD(4) + t503) * t273) * t377 + (t5 * t421 + (t36 * t402 + t393 * t6 + t505) * t276 + (t37 * t402 + t393 * t7 + t304) * t275 - g(1) * (-t276 * t443 + t418) - g(2) * (-t275 * t443 + t420) - g(3) * t515 - (-g(3) * t277 + t506 * (-t269 - t482)) * t272 - (t303 * t272 + (t36 * t432 + t37 * t431 + t283) * t273) * qJD(4) - t507 * t537 + (t37 * t241 - t36 * t242) * (-t272 * t277 - t514 + t515) + (-t432 * t241 - t431 * t242 - t409 + t532) * t33) * m(6) + (-t59 * (t177 * t242 - t412 * t514) - t60 * (-t177 * t241 - t413 * t514) - ((-t111 * t59 + t113 * t60) * t272 + (t59 * (t175 * t275 + t147) + t60 * (-t175 * t276 - t149) + t313) * t273) * qJD(4) - g(1) * (t260 + t417) - g(2) * (t259 + t419) - g(3) * (t177 + t514) - t506 * t272 * (-pkin(3) - t483) + t26 * t421 + (t26 * t113 + t34 * t423 + t466 * t59) * t276 + (t26 * t111 + t35 * t423 + t466 * t60) * t275 + (-t147 * t241 + t149 * t242 + t83 * t275 + t85 * t276 + t532) * t46) * m(5) + (-(t86 * (-t218 * t275 - t219 * t276) + t336 * t244) * qJD(3) + t61 * t335 - t86 * t334 + (-t114 * t276 - t115 * t275) * t243 + t336 * t238 + g(1) * t219 + g(2) * t218 - g(3) * t244) * m(4) + ((-t178 * t566 + t179 * t502) * t558 + (-t196 * t566 + t197 * t502) * t559 + t556) * t275 / 0.2e1 - ((t178 * t501 - t179 * t566) * t558 + (t196 * t501 - t197 * t566) * t559 + t557) * t276 / 0.2e1 + ((-t522 + t553) * t276 + (t521 + t554) * t275) * t378; (t509 * t272 - t273 * t575) * t500 + (t508 * t272 - t273 * t576) * t499 + (t272 * t510 - t273 * t519) * t498 + (t237 * t512 - t276 * t511 - t513 * t533) * t497 + (t579 * t273 + (t275 * t550 + t276 * t549) * t272 + (t273 * t509 + t530) * qJD(3)) * t496 + (t234 * t513 + t235 * t512 - t275 * t511) * t495 + (t580 * t273 + (t275 * t552 + t276 * t551) * t272 + (t273 * t508 + t529) * qJD(3)) * t494 - (t156 * t544 + t157 * t545 + t519 * t240 + t521 * t241 + t522 * t242 + t520 * t410) * t273 / 0.2e1 + t557 * t449 / 0.2e1 + t556 * t448 / 0.2e1 + t548 * t415 / 0.2e1 + (t520 * t273 + (t275 * t522 + t276 * t521) * t272 + (t272 * t519 + t273 * t510) * qJD(3)) * t378 + (t563 * t273 + (-t278 * t513 + t512 * t279) * t272) * t377 + t554 * t390 / 0.2e1 + t553 * t273 * t380 + (-(t33 * t429 - t36 * t371) * t242 - (t33 * t430 + t37 * t371) * t241 - (t36 * t430 + t37 * t429) * t410 + g(1) * t429 - g(2) * t430 + (qJD(3) * t283 + t36 * t486 - t37 * t485 - t467 * t7 + t468 * t6) * t273 + (-g(3) * (-t278 * t555 - t481) + t303 * qJD(3) + (-t37 * t484 - t428 * t7 + t304) * t276 + (t36 * t484 + t428 * t6 - t505) * t275) * t272) * m(6) + (-g(1) * t133 - g(2) * t131 - g(3) * t233 - t59 * (t131 * t410 + t233 * t242) - t60 * (-t133 * t410 - t233 * t241) - t46 * (t131 * t241 - t133 * t242) + (t34 * t111 - t35 * t113 + t59 * t83 - t60 * t85 + (t313 + (t275 * t59 - t276 * t60) * t175) * qJD(3)) * t273 + (t59 * (-qJD(3) * t111 + t275 * t97) + t60 * (qJD(3) * t113 - t276 * t97) + t26 * t341 + t46 * (-t275 * t85 + t276 * t83) + t365 * t175) * t272) * m(5); ((-t5 + t507 * qJD(3) - t36 * (-t242 + t388) - t37 * (t241 - t387) + g(3)) * t273 + ((-t241 * t275 + t242 * t276 + qJD(3)) * t33 + t275 * t7 + t276 * t6 - t506) * t272) * m(6);];
tau = t1;
