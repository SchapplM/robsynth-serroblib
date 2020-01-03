% Calculate vector of inverse dynamics joint torques for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR9_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:47
% DurationCPUTime: 23.06s
% Computational Cost: add. (14337->804), mult. (33870->1057), div. (0->0), fcn. (38306->8), ass. (0->366)
t312 = sin(qJ(4));
t314 = cos(qJ(4));
t311 = sin(qJ(5));
t313 = cos(qJ(5));
t491 = sin(pkin(8));
t492 = cos(pkin(8));
t506 = sin(qJ(1));
t507 = cos(qJ(1));
t261 = -t491 * t506 - t492 * t507;
t262 = t507 * t491 - t506 * t492;
t468 = t311 * t314;
t171 = -t261 * t313 + t262 * t468;
t471 = t262 * t312;
t466 = t313 * t314;
t170 = t261 * t311 + t262 * t466;
t488 = Icges(6,4) * t170;
t96 = -Icges(6,2) * t171 + Icges(6,6) * t471 + t488;
t162 = Icges(6,4) * t171;
t99 = Icges(6,1) * t170 + Icges(6,5) * t471 - t162;
t558 = t311 * t96 - t313 * t99;
t93 = Icges(6,5) * t170 - Icges(6,6) * t171 + Icges(6,3) * t471;
t33 = -t558 * t312 - t314 * t93;
t434 = qJD(5) * t312;
t442 = qJD(4) * t261;
t191 = t262 * t434 + t442;
t441 = qJD(4) * t262;
t192 = -t261 * t434 + t441;
t473 = t261 * t312;
t173 = t261 * t468 + t262 * t313;
t174 = -t261 * t466 + t262 * t311;
t559 = t173 * t96 + t174 * t99;
t26 = t473 * t93 - t559;
t163 = Icges(6,4) * t173;
t101 = Icges(6,1) * t174 - Icges(6,5) * t473 + t163;
t95 = Icges(6,5) * t174 + Icges(6,6) * t173 - Icges(6,3) * t473;
t487 = Icges(6,4) * t174;
t98 = Icges(6,2) * t173 - Icges(6,6) * t473 + t487;
t27 = t174 * t101 + t173 * t98 - t95 * t473;
t433 = qJD(5) * t314;
t292 = qJD(1) + t433;
t469 = t311 * t312;
t290 = Icges(6,4) * t469;
t467 = t312 * t313;
t483 = Icges(6,5) * t314;
t214 = -Icges(6,1) * t467 + t290 + t483;
t371 = Icges(6,5) * t313 - Icges(6,6) * t311;
t328 = -Icges(6,3) * t314 + t312 * t371;
t485 = Icges(6,4) * t313;
t374 = -Icges(6,2) * t311 + t485;
t329 = -Icges(6,6) * t314 + t312 * t374;
t60 = -t173 * t329 + t174 * t214 + t328 * t473;
t11 = -t191 * t26 + t192 * t27 + t60 * t292;
t373 = Icges(5,5) * t314 - Icges(5,6) * t312;
t145 = Icges(5,3) * t261 + t262 * t373;
t489 = Icges(5,4) * t314;
t376 = -Icges(5,2) * t312 + t489;
t148 = Icges(5,6) * t261 + t262 * t376;
t490 = Icges(5,4) * t312;
t379 = Icges(5,1) * t314 - t490;
t151 = Icges(5,5) * t261 + t262 * t379;
t546 = t148 * t312 - t151 * t314;
t52 = t261 * t145 - t262 * t546;
t24 = t170 * t99 - t171 * t96 + t471 * t93;
t392 = t170 * rSges(6,1) - t171 * rSges(6,2);
t103 = -rSges(6,3) * t471 - t392;
t501 = rSges(6,1) * t313;
t391 = rSges(6,2) * t311 - t501;
t216 = rSges(6,3) * t314 + t312 * t391;
t302 = qJD(2) * t506;
t505 = t312 * pkin(4);
t286 = pkin(7) * t314 - t505;
t437 = qJD(4) * t286;
t555 = t292 * t103 + t191 * t216 + t261 * t437 - t302;
t537 = t170 * t214 + t171 * t329 - t328 * t471;
t554 = t191 * t24 + t537 * t292;
t147 = Icges(5,3) * t262 - t261 * t373;
t551 = t262 * t147;
t472 = t261 * t314;
t156 = -rSges(5,1) * t472 + rSges(5,2) * t473 + t262 * rSges(5,3);
t201 = -t261 * pkin(3) + pkin(6) * t262;
t309 = t507 * pkin(2);
t443 = t507 * pkin(1) + t506 * qJ(2);
t532 = t309 + t443;
t544 = t201 + t532;
t550 = t156 + t544;
t372 = Icges(5,5) * t312 + Icges(5,6) * t314;
t176 = t372 * t261;
t375 = Icges(5,2) * t314 + t490;
t378 = Icges(5,1) * t312 + t489;
t362 = -t312 * t375 + t314 * t378;
t119 = t262 * t362 + t176;
t548 = qJD(1) * t119;
t369 = -t148 * t314 - t151 * t312;
t516 = -t191 / 0.2e1;
t175 = t372 * t262;
t199 = -t261 * rSges(4,1) - t262 * rSges(4,2);
t545 = t199 + t532;
t243 = t261 * qJD(1);
t244 = t262 * qJD(1);
t543 = t244 * pkin(3) + t243 * pkin(6);
t303 = qJD(2) * t507;
t542 = -t532 * qJD(1) + t303;
t436 = qJD(4) * t312;
t349 = t244 * t314 + t261 * t436;
t435 = qJD(4) * t314;
t476 = t244 * t312;
t350 = t261 * t435 - t476;
t92 = rSges(5,1) * t349 + rSges(5,2) * t350 + t243 * rSges(5,3);
t538 = t170 * t101 - t171 * t98;
t534 = t262 * rSges(4,1) - t261 * rSges(4,2);
t256 = t261 * pkin(6);
t404 = t262 * pkin(3) + t256;
t426 = t506 * pkin(2);
t427 = t506 * pkin(1);
t343 = -t427 - t426;
t531 = t343 + t426;
t444 = t507 * rSges(3,1) + t506 * rSges(3,3);
t530 = t443 + t444;
t305 = t507 * qJ(2);
t279 = t427 - t305;
t271 = qJD(1) * t279;
t446 = t302 - t271;
t412 = qJD(1) * t507;
t447 = qJ(2) * t412 + t302;
t529 = t447 - t446;
t235 = t244 * pkin(6);
t527 = -t235 + t542;
t526 = -qJD(1) * t404 + t447 + t543;
t459 = t375 * t262 - t151;
t461 = -t378 * t262 - t148;
t525 = t312 * t459 + t314 * t461;
t211 = -Icges(6,3) * t312 - t314 * t371;
t364 = -t214 * t313 - t311 * t329;
t380 = -t101 * t313 + t311 * t98;
t524 = t192 * (-t328 * t261 - t380) - t191 * (-t328 * t262 + t558) - t292 * (t211 + t364);
t249 = (Icges(6,1) * t311 + t485) * t312;
t324 = t192 * (Icges(6,1) * t173 - t487 - t98) - t191 * (Icges(6,1) * t171 + t488 + t96) - (-t329 - t249) * t292;
t315 = qJD(1) ^ 2;
t189 = qJD(4) * t243 + qJDD(4) * t262;
t430 = qJDD(5) * t312;
t107 = -qJD(5) * t350 - t261 * t430 + t189;
t523 = t107 / 0.2e1;
t190 = qJD(4) * t244 - qJDD(4) * t261;
t477 = t243 * t312;
t352 = t262 * t435 + t477;
t108 = -qJD(5) * t352 - t262 * t430 + t190;
t522 = t108 / 0.2e1;
t521 = t189 / 0.2e1;
t520 = t190 / 0.2e1;
t519 = -t192 / 0.2e1;
t518 = t192 / 0.2e1;
t517 = t191 / 0.2e1;
t515 = t243 / 0.2e1;
t514 = t244 / 0.2e1;
t260 = -qJD(4) * t434 + qJDD(5) * t314 + qJDD(1);
t513 = t260 / 0.2e1;
t510 = -t292 / 0.2e1;
t509 = t292 / 0.2e1;
t508 = rSges(6,3) + pkin(7);
t504 = t314 * pkin(4);
t351 = -t243 * t314 + t262 * t436;
t333 = -qJD(5) * t261 + t351;
t395 = t262 * t433 + t244;
t81 = -t311 * t333 + t313 * t395;
t82 = t311 * t395 + t313 * t333;
t42 = Icges(6,5) * t82 + Icges(6,6) * t81 - Icges(6,3) * t352;
t44 = Icges(6,4) * t82 + Icges(6,2) * t81 - Icges(6,6) * t352;
t46 = Icges(6,1) * t82 + Icges(6,4) * t81 - Icges(6,5) * t352;
t8 = (-qJD(4) * t558 + t42) * t314 + (qJD(4) * t93 + t311 * t44 - t313 * t46 + (-t311 * t99 - t313 * t96) * qJD(5)) * t312;
t503 = t8 * t191;
t332 = qJD(5) * t262 + t349;
t396 = t261 * t433 + t243;
t83 = -t311 * t332 + t313 * t396;
t84 = t311 * t396 + t313 * t332;
t43 = Icges(6,5) * t84 + Icges(6,6) * t83 - Icges(6,3) * t350;
t45 = Icges(6,4) * t84 + Icges(6,2) * t83 - Icges(6,6) * t350;
t47 = Icges(6,1) * t84 + Icges(6,4) * t83 - Icges(6,5) * t350;
t9 = (qJD(4) * t380 + t43) * t314 + (-qJD(4) * t95 + t311 * t45 - t313 * t47 + (t101 * t311 + t313 * t98) * qJD(5)) * t312;
t502 = t9 * t192;
t500 = rSges(6,3) * t312;
t499 = t244 * rSges(5,3);
t250 = t261 * rSges(5,3);
t498 = t33 * t108;
t34 = t312 * t380 + t314 * t95;
t497 = t34 * t107;
t123 = t312 * t364 - t314 * t328;
t247 = (Icges(6,5) * t311 + Icges(6,6) * t313) * t312;
t158 = qJD(4) * t211 + qJD(5) * t247;
t213 = -Icges(6,6) * t312 - t314 * t374;
t486 = Icges(6,4) * t311;
t159 = (Icges(6,2) * t313 + t486) * t434 + t213 * qJD(4);
t377 = Icges(6,1) * t313 - t486;
t215 = -Icges(6,5) * t312 - t314 * t377;
t160 = qJD(4) * t215 + qJD(5) * t249;
t41 = (qJD(4) * t364 + t158) * t314 + (qJD(4) * t328 + t159 * t311 - t160 * t313 + (t214 * t311 - t313 * t329) * qJD(5)) * t312;
t496 = t123 * t260 + t41 * t292;
t112 = pkin(4) * t349 - pkin(7) * t350;
t49 = t84 * rSges(6,1) + t83 * rSges(6,2) - rSges(6,3) * t350;
t495 = t112 + t49;
t150 = Icges(5,6) * t262 - t261 * t376;
t478 = t150 * t312;
t470 = t262 * t314;
t287 = -pkin(7) * t312 - t504;
t186 = t287 * t262;
t465 = t103 + t186;
t104 = t174 * rSges(6,1) + t173 * rSges(6,2) - rSges(6,3) * t473;
t188 = -pkin(4) * t472 - pkin(7) * t473;
t464 = t104 + t188;
t153 = Icges(5,5) * t262 - t261 * t379;
t463 = -t261 * t147 - t153 * t470;
t462 = -t153 * t472 + t551;
t460 = -t378 * t261 + t150;
t458 = t375 * t261 + t153;
t259 = (rSges(6,1) * t311 + rSges(6,2) * t313) * t312;
t161 = qJD(5) * t259 + (t314 * t391 - t500) * qJD(4);
t269 = qJD(4) * t287;
t457 = -t161 - t269;
t242 = qJD(1) * t443 - t303;
t456 = t243 * pkin(3) - t235 - t242;
t424 = rSges(6,2) * t469;
t455 = -rSges(6,3) * t470 - t262 * t424;
t454 = -rSges(6,3) * t472 - t261 * t424;
t452 = t216 + t286;
t451 = t244 * rSges(4,1) - t243 * rSges(4,2);
t184 = pkin(4) * t470 + pkin(7) * t471;
t450 = -t375 + t379;
t449 = -t376 - t378;
t448 = qJD(1) * t303 + qJDD(2) * t506;
t282 = -rSges(5,1) * t314 + rSges(5,2) * t312;
t266 = t282 * qJD(4);
t440 = qJD(4) * t266;
t439 = qJD(4) * t269;
t393 = rSges(5,1) * t312 + rSges(5,2) * t314;
t438 = qJD(4) * t393;
t432 = t373 * qJD(1);
t431 = m(4) + m(5) + m(6);
t25 = -t95 * t471 - t538;
t429 = -t507 / 0.2e1;
t428 = t506 / 0.2e1;
t425 = rSges(6,1) * t467;
t422 = t506 * rSges(3,1);
t411 = qJD(1) * t506;
t410 = -t442 / 0.2e1;
t409 = t442 / 0.2e1;
t408 = -t441 / 0.2e1;
t407 = t441 / 0.2e1;
t406 = -t435 / 0.2e1;
t405 = t433 / 0.2e1;
t403 = -t279 - t426;
t401 = -rSges(6,1) * t466 + rSges(6,2) * t468;
t398 = t82 * rSges(6,1) + t81 * rSges(6,2);
t394 = t243 * rSges(4,1) + t244 * rSges(4,2);
t367 = t150 * t314 + t153 * t312;
t88 = Icges(5,4) * t349 + Icges(5,2) * t350 + Icges(5,6) * t243;
t90 = Icges(5,1) * t349 + Icges(5,4) * t350 + Icges(5,5) * t243;
t322 = qJD(4) * t367 + t312 * t88 - t314 * t90;
t87 = Icges(5,4) * t351 + Icges(5,2) * t352 + Icges(5,6) * t244;
t89 = Icges(5,1) * t351 + Icges(5,4) * t352 + Icges(5,5) * t244;
t323 = qJD(4) * t369 + t312 * t87 - t314 * t89;
t366 = -t153 * t314 + t478;
t85 = Icges(5,5) * t351 + Icges(5,6) * t352 + Icges(5,3) * t244;
t86 = Icges(5,5) * t349 + Icges(5,6) * t350 + Icges(5,3) * t243;
t390 = -(-t145 * t244 - t243 * t546 - t261 * t85 + t262 * t323) * t261 + (t147 * t244 + t243 * t366 - t261 * t86 + t262 * t322) * t262;
t389 = -(-t145 * t243 + t244 * t546 + t261 * t323 + t262 * t85) * t261 + (t147 * t243 - t244 * t366 + t261 * t322 + t262 * t86) * t262;
t388 = -t24 * t262 - t25 * t261;
t387 = -t26 * t262 - t261 * t27;
t386 = -t261 * t34 - t262 * t33;
t53 = t150 * t471 + t463;
t385 = -t261 * t52 + t262 * t53;
t54 = -t145 * t262 - t148 * t473 + t151 * t472;
t55 = t150 * t473 + t462;
t384 = -t261 * t54 + t262 * t55;
t155 = t262 * t282 - t250;
t361 = t403 + t404;
t346 = -t155 + t361;
t64 = qJD(1) * t346 + t261 * t438 + t302;
t65 = qJD(1) * t550 + t262 * t438 - t303;
t383 = -t261 * t64 - t262 * t65;
t91 = rSges(5,1) * t351 + rSges(5,2) * t352 + t499;
t382 = -t261 * t92 - t262 * t91;
t370 = -t103 * t261 + t104 * t262;
t365 = t155 * t262 + t156 * t261;
t363 = -t312 * t378 - t314 * t375;
t359 = t534 + t403;
t154 = rSges(5,1) * t470 - rSges(5,2) * t471 + t250;
t358 = pkin(3) - t282;
t357 = -t312 * t42 + t435 * t93;
t356 = -t312 * t43 - t435 * t95;
t347 = -t309 * t315 + t448;
t345 = -t186 + t361;
t342 = -t422 - t427;
t341 = qJDD(1) * t443 - qJDD(2) * t507 + (-pkin(1) * t411 + t302 + t447) * qJD(1);
t285 = rSges(2,1) * t507 - rSges(2,2) * t506;
t281 = rSges(2,1) * t506 + rSges(2,2) * t507;
t338 = -t191 * t93 - t192 * t95 + t292 * t328;
t337 = -(Icges(6,5) * t171 + Icges(6,6) * t170) * t191 + (Icges(6,5) * t173 - Icges(6,6) * t174) * t192 + t247 * t292;
t336 = t305 + t343;
t335 = t312 * t458 + t314 * t460;
t334 = t256 + t336;
t331 = t312 * t337;
t330 = t312 * t377 - t483;
t327 = (t312 * t449 + t314 * t450) * qJD(1);
t325 = t336 + t404;
t321 = (-Icges(6,2) * t174 + t101 + t163) * t192 - (Icges(6,2) * t170 + t162 - t99) * t191 + (Icges(6,2) * t467 + t214 + t290) * t292;
t320 = qJDD(1) * t309 - t315 * t426 + t341;
t264 = t376 * qJD(4);
t265 = t379 * qJD(4);
t319 = qJD(4) * t363 - t264 * t312 + t265 * t314;
t28 = t103 * t192 + t104 * t191 - qJD(3) + (t186 * t262 + t188 * t261) * qJD(4);
t39 = qJD(1) * t345 - t555;
t40 = -t262 * t437 + t292 * t104 - t192 * t216 - t303 + (t188 + t544) * qJD(1);
t318 = t28 * t370 + (t261 * t40 - t262 * t39) * t216;
t317 = qJD(1) * t543 + qJDD(1) * t201 + t320;
t316 = t524 * t312;
t307 = t507 * rSges(3,3);
t300 = rSges(3,3) * t412;
t280 = t422 - t307;
t263 = t373 * qJD(4);
t240 = pkin(7) * t472;
t238 = pkin(7) * t470;
t217 = t401 - t500;
t187 = pkin(4) * t473 - t240;
t185 = pkin(4) * t471 - t238;
t182 = t393 * t261;
t181 = t393 * t262;
t143 = qJD(1) * t359 + t302;
t142 = t261 * t425 + t454;
t141 = t262 * t425 + t455;
t140 = t330 * t261;
t139 = t330 * t262;
t138 = t329 * t261;
t137 = t329 * t262;
t126 = qJDD(1) * t444 + qJD(1) * (-rSges(3,1) * t411 + t300) + t341;
t125 = -qJD(1) * t242 - t315 * t444 + t448 + (-t279 - t280) * qJDD(1);
t122 = rSges(6,1) * t173 - rSges(6,2) * t174;
t121 = rSges(6,1) * t171 + rSges(6,2) * t170;
t120 = t261 * t362 - t175;
t111 = pkin(4) * t351 - pkin(7) * t352;
t110 = t120 * qJD(1);
t73 = qJD(1) * t451 + qJDD(1) * t199 + t320;
t72 = (-t242 + t394) * qJD(1) + t359 * qJDD(1) + t347;
t63 = qJD(4) * t365 - qJD(3);
t51 = -t243 * t372 - t244 * t362 + t261 * t319 - t262 * t263;
t50 = t243 * t362 - t244 * t372 + t261 * t263 + t262 * t319;
t48 = -rSges(6,3) * t352 + t398;
t32 = qJD(4) * t366 - t312 * t90 - t314 * t88;
t31 = -qJD(4) * t546 - t312 * t89 - t314 * t87;
t30 = qJD(1) * t92 + qJDD(1) * t156 + t189 * t393 - t262 * t440 + t317;
t29 = -t261 * t440 - t190 * t393 + (-t91 + t456) * qJD(1) + t346 * qJDD(1) + t347;
t23 = qJD(4) * t382 - t155 * t189 + t156 * t190 + qJDD(3);
t22 = qJD(4) * t384 + t110;
t21 = qJD(4) * t385 + t548;
t20 = -t158 * t473 + t159 * t173 + t160 * t174 + t214 * t84 + t328 * t350 - t329 * t83;
t19 = -t158 * t471 + t159 * t171 - t160 * t170 + t214 * t82 + t328 * t352 - t329 * t81;
t14 = t292 * t123 - t191 * t33 + t192 * t34;
t13 = qJD(1) * t112 + qJDD(1) * t188 + t260 * t104 - t107 * t216 - t192 * t161 - t189 * t286 - t262 * t439 + t292 * t49 + t317;
t12 = -t261 * t439 - t260 * t103 + t108 * t216 - t191 * t161 + t190 * t286 - t292 * t48 + (-t111 + t456) * qJD(1) + t345 * qJDD(1) + t347;
t10 = t192 * t25 - t554;
t7 = -t107 * t103 + t104 * t108 - t111 * t441 - t112 * t442 - t189 * t186 + t188 * t190 - t191 * t49 - t192 * t48 + qJDD(3);
t6 = t101 * t84 + t173 * t45 + t174 * t47 + t261 * t356 + t476 * t95 + t83 * t98;
t5 = t173 * t44 + t174 * t46 + t261 * t357 - t476 * t93 - t83 * t96 - t84 * t99;
t4 = t101 * t82 - t170 * t47 + t171 * t45 + t262 * t356 - t477 * t95 + t81 * t98;
t3 = -t170 * t46 + t171 * t44 + t262 * t357 + t477 * t93 - t81 * t96 - t82 * t99;
t2 = t107 * t27 + t108 * t26 - t191 * t5 + t192 * t6 + t20 * t292 + t260 * t60;
t1 = t107 * t25 + t108 * t24 + t19 * t292 - t191 * t3 + t192 * t4 - t260 * t537;
t15 = [((-g(2) + t73) * t545 + (t394 + t542) * t143 + (t143 + t451 + (t531 - t534) * qJD(1) + t529) * (qJD(1) * t545 - t303) + (-g(1) + t72) * (t336 + t534)) * m(4) - t537 * t522 + (m(2) * (t281 ^ 2 + t285 ^ 2) - t363 + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJDD(1) + (t11 + t19) * t516 + (t119 - t369) * t520 + (t120 - t367) * t521 + (t31 + t50 + t22) * t410 + (t32 + t51) * t407 + t11 * t517 + (t110 + ((t53 - t54 - t463) * t261 + t462 * t262) * qJD(4)) * t409 + ((t26 + (-t261 * t93 + t262 * t95) * t312 + t538 + t559) * t192 + t10 + t554) * t519 + (-t63 * (t154 + t155) * t442 - g(1) * (t154 + t325) + (t358 * t262 + t250 + t334) * t29 + (t358 * t243 - t393 * t441 - t499 + t527) * t64 + (-t393 * t442 - t446 + t64 + (-t154 + t531) * qJD(1) + t526 + t92) * t65 + (t30 - g(2)) * t550) * m(5) + (qJD(4) * t362 + t264 * t314 + t265 * t312) * qJD(1) + ((t126 - g(2)) * t530 + (t125 - g(1)) * (t305 + t307 + t342) + (t300 + (t280 + t342) * qJD(1) + t529) * (qJD(1) * t530 - t303)) * m(3) + t497 / 0.2e1 + t498 / 0.2e1 + (t21 + ((t54 + (t145 - t366) * t262) * t262 + (t55 + (t145 - t478) * t261 + t551 - t462) * t261) * qJD(4) - t548) * t408 + t20 * t518 + t60 * t523 + t496 + (-g(1) * (t325 - t103 + t184) - t28 * (t184 + t186) * t442 + (t13 - g(2)) * (t544 + t464) + (t334 + t392 + (t312 * t508 + pkin(3) + t504) * t262) * t12 + (-t398 - (-pkin(3) + t287 - t500) * t243 + (t314 * t508 - t505) * t441 + t527) * t39 + (t271 + t39 + t495 + (-t184 + t531) * qJD(1) + t526 + t555) * t40) * m(6) - t503 / 0.2e1 + t502 / 0.2e1 - m(2) * (-g(1) * t281 + g(2) * t285); (-m(3) - t431) * (g(1) * t506 - g(2) * t507) + 0.2e1 * (t12 * t428 + t13 * t429) * m(6) + 0.2e1 * (t29 * t428 + t30 * t429) * m(5) + 0.2e1 * (t428 * t72 + t429 * t73) * m(4) + 0.2e1 * (t125 * t428 + t126 * t429) * m(3); m(4) * qJDD(3) + m(5) * t23 + m(6) * t7 + g(3) * t431; (t262 * t405 + t514) * t10 + ((t138 * t171 - t140 * t170) * t192 - (t137 * t171 - t139 * t170) * t191 + (-t170 * t215 + t171 * t213) * t292 + (-t25 * t472 + t312 * t537) * qJD(5) + ((-qJD(5) * t24 + t338) * t314 + t316) * t262) * t517 + (((t138 * t311 - t140 * t313 - t95) * t192 - (t137 * t311 - t139 * t313 + t93) * t191 + (t213 * t311 - t215 * t313 + t328) * t292 - t123 * qJD(5)) * t312 + (qJD(5) * t386 - t524) * t314) * t510 + (t261 * t405 + t515) * t11 + (-g(1) * (-t240 + t454) - g(2) * (-t238 + t455) - g(3) * (t401 - t504) - (-g(3) * t508 + (g(1) * t261 + g(2) * t262) * (pkin(4) + t501)) * t312 - t39 * (-qJD(1) * t185 - t141 * t292 - t191 * t217 - t261 * t269) - t40 * (qJD(1) * t187 + t142 * t292 - t192 * t217 - t262 * t269) - t28 * (t141 * t192 + t142 * t191 + t185 * t441 + t187 * t442) - ((t103 * t39 - t104 * t40) * t312 + t318 * t314) * qJD(5) + (-t28 * t464 + t39 * t452) * t244 - (-t28 * t465 + t40 * t452) * t243 + (-t13 * t452 + t40 * t457 - t7 * t465 + t28 * (t111 + t48)) * t262 + (-t12 * t452 + t28 * t495 + t39 * t457 - t464 * t7) * t261) * m(6) + ((t138 * t173 + t140 * t174) * t192 - (t137 * t173 + t139 * t174) * t191 + (t173 * t213 + t174 * t215) * t292 + (-t26 * t470 - t312 * t60) * qJD(5) + ((-qJD(5) * t27 + t338) * t314 + t316) * t261) * t519 + (-(-t181 * t64 + t182 * t65) * qJD(1) - (t63 * (t181 * t262 + t182 * t261) + t383 * t282) * qJD(4) - t23 * t365 + t63 * (t155 * t243 - t156 * t244 - t382) + t383 * t266 - (-t243 * t65 + t244 * t64 - t261 * t29 - t262 * t30) * t393 - g(1) * t182 - g(2) * t181 - g(3) * t282) * m(5) + qJD(1) * (-t243 * t367 - t244 * t369 - t261 * t31 + t262 * t32) / 0.2e1 + qJDD(1) * (t261 * t369 - t262 * t367) / 0.2e1 - qJD(1) * ((t312 * t450 - t314 * t449) * qJD(1) + ((t261 * t459 - t262 * t458) * t314 + (-t261 * t461 + t262 * t460) * t312) * qJD(4)) / 0.2e1 + t14 * t434 / 0.2e1 + ((t175 * t442 + t432) * t261 + (t327 + (t335 * t262 + (-t176 - t525) * t261) * qJD(4)) * t262) * t409 + ((t176 * t441 - t432) * t262 + (t327 + (-t525 * t261 + (-t175 + t335) * t262) * qJD(4)) * t261) * t408 + (t243 * t53 + t244 * t52 + t390) * t410 + (t243 * t34 + t244 * t33 - t261 * t8 + t262 * t9) * t509 + (t243 * t55 + t244 * t54 + t389) * t407 + (t24 * t244 + t243 * t25 - t261 * t3 + t262 * t4) * t516 + (t243 * t27 + t244 * t26 - t261 * t5 + t262 * t6) * t518 + (qJD(1) * t51 + qJD(4) * t389 + qJDD(1) * t120 + t189 * t55 + t190 * t54 + t2) * t262 / 0.2e1 - (qJD(1) * t50 + qJD(4) * t390 + qJDD(1) * t119 + t189 * t53 + t190 * t52 + t1) * t261 / 0.2e1 + (-t261 * t33 + t262 * t34) * t513 + t21 * t514 + t22 * t515 + t385 * t520 + t384 * t521 + (-t24 * t261 + t25 * t262) * t522 + (-t26 * t261 + t262 * t27) * t523; -t2 * t473 / 0.2e1 + (t312 * t387 + t314 * t60) * t523 + ((qJD(4) * t387 + t20) * t314 + (-qJD(4) * t60 - t243 * t26 + t244 * t27 - t261 * t6 - t262 * t5) * t312) * t518 - t1 * t471 / 0.2e1 + (t312 * t388 - t314 * t537) * t522 + ((qJD(4) * t388 + t19) * t314 + (qJD(4) * t537 - t24 * t243 + t244 * t25 - t261 * t4 - t262 * t3) * t312) * t516 - t14 * t436 / 0.2e1 + t314 * (t496 + t497 + t498 + t502 - t503) / 0.2e1 + (t123 * t314 + t312 * t386) * t513 + ((qJD(4) * t386 + t41) * t314 + (-qJD(4) * t123 - t243 * t33 + t244 * t34 - t261 * t9 - t262 * t8) * t312) * t509 + (t173 * t321 + t174 * t324 - t261 * t331) * t519 + (-t170 * t324 + t171 * t321 - t262 * t331) * t517 + (t337 * t314 + (t321 * t311 - t313 * t324) * t312) * t510 + (t476 / 0.2e1 + t261 * t406) * t11 + (-t477 / 0.2e1 + t262 * t406) * t10 + ((qJD(4) * t318 - t12 * t103 + t13 * t104 - t39 * t48 + t40 * t49) * t314 + (t39 * (qJD(4) * t103 - t161 * t262) + t40 * (-qJD(4) * t104 + t161 * t261) - t7 * t370 + t28 * (t103 * t244 + t104 * t243 - t261 * t48 + t262 * t49) + (-t12 * t262 + t13 * t261 - t243 * t39 - t244 * t40) * t216) * t312 - t39 * (-t121 * t292 - t191 * t259) - t40 * (t122 * t292 - t192 * t259) - t28 * (t121 * t192 + t122 * t191) - g(1) * t122 - g(2) * t121 - g(3) * t259) * m(6);];
tau = t15;
