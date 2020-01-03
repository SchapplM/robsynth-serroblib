% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP5
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:39
% EndTime: 2019-12-31 17:16:50
% DurationCPUTime: 6.97s
% Computational Cost: add. (19582->413), mult. (24360->543), div. (0->0), fcn. (22692->6), ass. (0->261)
t363 = qJ(2) + qJ(3);
t346 = cos(t363);
t345 = sin(t363);
t481 = Icges(4,4) * t345;
t306 = Icges(4,1) * t346 - t481;
t365 = sin(qJ(1));
t367 = cos(qJ(1));
t235 = Icges(4,5) * t365 + t306 * t367;
t340 = Icges(5,5) * t345;
t550 = Icges(5,1) * t346 + t340;
t583 = Icges(5,4) * t365 + t367 * t550 + t235;
t341 = Icges(4,4) * t346;
t305 = Icges(4,1) * t345 + t341;
t480 = Icges(5,5) * t346;
t582 = Icges(5,1) * t345 + t305 - t480;
t384 = Icges(5,3) * t346 - t340;
t578 = Icges(4,2) * t346 + t384 + t481;
t310 = pkin(3) * t346 + qJ(4) * t345;
t311 = rSges(5,1) * t346 + rSges(5,3) * t345;
t569 = t310 + t311;
t182 = t569 * t365;
t464 = t345 * t365;
t508 = rSges(5,1) + pkin(3);
t417 = t508 * t464;
t457 = t346 * t367;
t486 = rSges(5,3) + qJ(4);
t419 = t486 * t457;
t458 = t346 * t365;
t463 = t345 * t367;
t129 = (-t508 * t463 + t419) * t367 + (t486 * t458 - t417) * t365;
t409 = t486 * t346;
t412 = t508 * t345;
t422 = t412 - t409;
t364 = sin(qJ(2));
t504 = pkin(2) * t364;
t395 = t422 + t504;
t172 = t395 * t365;
t174 = t395 * t367;
t184 = t569 * t367;
t359 = t367 * rSges(5,2);
t362 = t367 ^ 2;
t495 = t365 * rSges(5,2);
t116 = t367 * (t311 * t367 + t495) + t362 * t310 + (-t359 + t182) * t365;
t366 = cos(qJ(2));
t503 = pkin(2) * t366;
t342 = pkin(1) + t503;
t536 = -pkin(6) - pkin(5);
t343 = t365 * t536;
t360 = t367 * pkin(5);
t418 = -t365 * t342 - t367 * t536;
t437 = -t365 * (pkin(1) * t365 - t360 + t418) + t367 * (-t365 * pkin(5) - t343 + (-pkin(1) + t342) * t367);
t98 = t116 + t437;
t43 = t98 * t129 + t172 * t182 + t174 * t184;
t237 = rSges(4,1) * t458 - rSges(4,2) * t464 - t367 * rSges(4,3);
t407 = -rSges(4,2) * t463 + t365 * rSges(4,3);
t155 = t365 * t237 + t367 * (rSges(4,1) * t457 + t407);
t112 = t155 + t437;
t309 = rSges(4,1) * t345 + rSges(4,2) * t346;
t285 = t309 * t365;
t287 = t309 * t367;
t163 = -t365 * t285 - t367 * t287;
t496 = rSges(4,1) * t346;
t312 = -rSges(4,2) * t345 + t496;
t374 = t309 + t504;
t551 = t374 * t367;
t552 = t374 * t365;
t61 = t112 * t163 + (t365 * t552 + t367 * t551) * t312;
t581 = -m(4) * t61 - m(5) * t43;
t327 = Icges(5,5) * t457;
t225 = Icges(5,6) * t365 + Icges(5,3) * t463 + t327;
t299 = Icges(4,5) * t346 - Icges(4,6) * t345;
t472 = t299 * t367;
t227 = Icges(4,3) * t365 + t472;
t300 = Icges(5,4) * t346 + Icges(5,6) * t345;
t471 = t300 * t367;
t580 = t225 * t463 + t583 * t457 + (Icges(5,2) * t365 + t227 + t471) * t365;
t302 = -Icges(4,2) * t345 + t341;
t579 = t302 + t582;
t577 = (-Icges(4,6) + Icges(5,6)) * t346 + (-Icges(5,4) - Icges(4,5)) * t345;
t576 = -t578 * t367 + t583;
t232 = -Icges(5,4) * t367 + t365 * t550;
t328 = Icges(4,4) * t464;
t234 = Icges(4,1) * t458 - Icges(4,5) * t367 - t328;
t575 = -Icges(4,2) * t458 - t384 * t365 + t232 + t234 - t328;
t231 = Icges(4,6) * t365 + t302 * t367;
t574 = -Icges(5,1) * t463 - t305 * t367 + t225 - t231 + t327;
t298 = Icges(5,3) * t345 + t480;
t224 = -Icges(5,6) * t367 + t298 * t365;
t230 = Icges(4,4) * t458 - Icges(4,2) * t464 - Icges(4,6) * t367;
t573 = t582 * t365 - t224 + t230;
t572 = t550 + t306;
t571 = -t231 * t463 + t580;
t553 = (t224 * t345 + t232 * t346) * t365;
t568 = t553 + t580;
t538 = m(4) / 0.2e1;
t537 = m(5) / 0.2e1;
t511 = t365 / 0.2e1;
t509 = -t367 / 0.2e1;
t564 = t367 / 0.2e1;
t451 = t365 * t300;
t204 = t365 * (-Icges(5,2) * t367 + t451);
t125 = t224 * t463 + t232 * t457 + t204;
t563 = t125 * t367;
t482 = Icges(3,4) * t364;
t316 = Icges(3,2) * t366 + t482;
t319 = Icges(3,1) * t366 - t482;
t562 = (t319 / 0.2e1 - t316 / 0.2e1) * t364;
t561 = t577 * t365;
t560 = t577 * t367;
t361 = t365 ^ 2;
t415 = t361 + t362;
t559 = (t572 - t578) * t346 + (t298 - t579) * t345;
t164 = (-t409 + t504) * t365 + t417;
t165 = (-t412 - t504) * t367 + t419;
t181 = t422 * t365;
t183 = t422 * t367;
t176 = -t237 + t418;
t177 = -t343 + (t342 + t496) * t367 + t407;
t376 = (-t176 * t367 - t177 * t365) * t312;
t546 = -t486 * t345 - t508 * t346;
t153 = t546 * t365 + t359 + t418;
t154 = t495 - t343 + (t342 - t546) * t367;
t444 = -t184 * t153 - t182 * t154;
t501 = (-t164 * t183 - t165 * t181 + t444) * t537 + (t376 + (t365 * t551 - t367 * t552) * t309) * t538;
t178 = -t365 * t409 + t417;
t179 = -t367 * t412 + t419;
t502 = (-t172 * t179 - t174 * t178 + t444) * t537 + (-t285 * t551 + t287 * t552 + t376) * t538;
t558 = t501 - t502;
t556 = (t574 * t365 + t573 * t367) * t346 + (-t576 * t365 + t575 * t367) * t345;
t356 = Icges(3,4) * t366;
t317 = -Icges(3,2) * t364 + t356;
t318 = Icges(3,1) * t364 + t356;
t51 = t116 * t129 + t181 * t182 + t183 * t184;
t93 = t415 * t309 * t312 + t155 * t163;
t549 = m(4) * t93 + m(5) * t51;
t271 = Icges(3,5) * t365 + t319 * t367;
t423 = -t316 * t367 + t271;
t455 = t364 * t365;
t337 = Icges(3,4) * t455;
t450 = t365 * t366;
t270 = Icges(3,1) * t450 - Icges(3,5) * t367 - t337;
t424 = -Icges(3,2) * t450 + t270 - t337;
t269 = Icges(3,6) * t365 + t317 * t367;
t425 = -t318 * t367 - t269;
t268 = Icges(3,4) * t450 - Icges(3,2) * t455 - Icges(3,6) * t367;
t426 = t318 * t365 + t268;
t545 = (-t423 * t365 + t424 * t367) * t364 + (t425 * t365 + t426 * t367) * t366;
t373 = (-t298 / 0.2e1 + t579 / 0.2e1) * t346 + (-t578 / 0.2e1 + t572 / 0.2e1) * t345;
t187 = t235 * t458;
t400 = t227 * t367 - t187;
t124 = -t231 * t464 - t400;
t226 = Icges(4,5) * t458 - Icges(4,6) * t464 - Icges(4,3) * t367;
t440 = -t365 * t226 - t234 * t457;
t127 = -t230 * t463 - t440;
t398 = t231 * t345 - t226;
t474 = t230 * t345;
t375 = (-t127 * t367 + t571 * t365 - t563) * t564 + (-t563 + (t124 - t187 + (t227 + t474) * t367 + t440) * t367 + (-t553 + t568) * t365) * t509 + (((t398 + t226) * t367 - t568 + t571) * t367 + (t125 - t204 + t127 + t400 - (t234 * t346 - t474) * t367 + t124 + t398 * t365) * t365) * t511;
t542 = 4 * qJD(1);
t541 = 2 * qJD(2);
t539 = 2 * qJD(3);
t288 = t415 * t345;
t442 = -t172 * t458 - t174 * t457;
t57 = t98 * t288 + t442;
t441 = -t181 * t458 - t183 * t457;
t72 = t116 * t288 + t441;
t527 = m(5) * (t72 + t57);
t393 = -t129 * t346 - t182 * t464 - t184 * t463;
t413 = t345 * t98 + t442;
t526 = m(5) * (t393 + t413);
t45 = t116 * t345 + t393 + t441;
t524 = m(5) * t45;
t443 = t153 * t457 + t154 * t458;
t519 = m(5) * ((t164 * t367 + t165 * t365) * t345 + t443);
t518 = m(5) * (-t172 * t463 + t174 * t464 + t443);
t517 = m(5) * ((t178 * t367 + t179 * t365) * t345 + t443);
t516 = m(5) * (-t181 * t463 + t183 * t464 + t443);
t515 = m(5) * (t153 * t164 + t154 * t165);
t513 = m(5) * (t153 * t178 + t154 * t179);
t512 = -t365 / 0.2e1;
t497 = rSges(3,1) * t366;
t411 = pkin(1) + t497;
t416 = rSges(3,2) * t455 + t367 * rSges(3,3);
t196 = -t411 * t365 + t360 + t416;
t454 = t364 * t367;
t339 = rSges(3,2) * t454;
t197 = -t339 + t411 * t367 + (rSges(3,3) + pkin(5)) * t365;
t320 = rSges(3,1) * t364 + rSges(3,2) * t366;
t295 = t320 * t365;
t296 = t320 * t367;
t507 = m(3) * (t196 * t295 - t197 * t296);
t506 = m(4) * (t176 * t552 - t177 * t551);
t505 = m(4) * (t176 * t285 - t177 * t287);
t500 = m(5) * qJD(2);
t499 = m(5) * qJD(3);
t498 = m(5) * qJD(4);
t465 = t345 * t346;
t456 = t364 * t268;
t449 = t366 * t367;
t266 = Icges(3,5) * t450 - Icges(3,6) * t455 - Icges(3,3) * t367;
t436 = -t365 * t266 - t270 * t449;
t387 = Icges(3,5) * t366 - Icges(3,6) * t364;
t267 = Icges(3,3) * t365 + t367 * t387;
t435 = t365 * t267 + t271 * t449;
t420 = t415 * t465;
t101 = -t153 * t464 + t154 * t463;
t414 = m(5) * t101 * qJD(1);
t410 = -t312 - t503;
t408 = ((-t561 * t365 + t556) * t367 + t560 * t361) * t511 + ((-t560 * t367 + t556) * t365 + t561 * t362) * t509;
t214 = t271 * t450;
t399 = t267 * t367 - t214;
t397 = t364 * t269 - t266;
t396 = t415 * t504;
t394 = -t569 - t503;
t386 = -Icges(3,5) * t364 - Icges(3,6) * t366;
t369 = -t375 + (t365 * t299 + t574 * t345 + t576 * t346 + t559 * t367 + t451) * t511 + (-t573 * t345 + t575 * t346 + t559 * t365 - t471 - t472) * t509;
t368 = -t373 + ((t224 + t230) * t346 + (-t232 + t234) * t345) * (t511 + t512);
t322 = -rSges(3,2) * t364 + t497;
t290 = t386 * t367;
t289 = t386 * t365;
t222 = t410 * t367;
t220 = t410 * t365;
t180 = t420 - t465;
t175 = t394 * t367;
t173 = t394 * t365;
t171 = t180 * t498;
t152 = -t396 + t163;
t136 = -t269 * t454 + t435;
t135 = -t268 * t454 - t436;
t134 = -t269 * t455 - t399;
t118 = (-t288 * t346 - t180 + t420) * t498;
t114 = -t396 + t129;
t95 = -t135 * t367 + t136 * t365;
t94 = -(-t365 * (-t366 * t270 + t456) - t266 * t367) * t367 + t134 * t365;
t68 = t516 / 0.2e1;
t66 = t517 / 0.2e1;
t64 = t518 / 0.2e1;
t62 = t519 / 0.2e1;
t44 = t524 / 0.2e1;
t39 = t526 / 0.2e1;
t38 = (t134 - t214 + (t267 + t456) * t367 + t436) * t367 + t435 * t365;
t37 = (t397 * t367 + t136 - t435) * t367 + (t397 * t365 + t135 + t399) * t365;
t35 = t373 + t505 + t513;
t33 = t527 / 0.2e1;
t20 = t68 - t517 / 0.2e1;
t19 = t68 + t66;
t18 = t66 - t516 / 0.2e1;
t17 = (t318 / 0.2e1 + t317 / 0.2e1) * t366 + t562 + t507 + t506 + t515 + t373;
t16 = t64 - t519 / 0.2e1;
t15 = t64 + t62;
t14 = t62 - t518 / 0.2e1;
t11 = t33 + t44 - t526 / 0.2e1;
t10 = t33 + t39 - t524 / 0.2e1;
t9 = t39 + t44 - t527 / 0.2e1;
t8 = t408 + t549;
t7 = t8 * qJD(3);
t6 = t408 - t581;
t4 = (t95 / 0.2e1 - t38 / 0.2e1) * t367 + (t37 / 0.2e1 + t94 / 0.2e1) * t365 + t375;
t3 = t375 - t558;
t2 = t375 + t558;
t1 = t369 + t501 + t502;
t5 = [t17 * qJD(2) + t35 * qJD(3) + t101 * t498, t17 * qJD(1) + t1 * qJD(3) + t15 * qJD(4) + (m(3) * ((-t196 * t367 - t197 * t365) * t322 + (-t295 * t367 + t296 * t365) * t320) / 0.2e1 + (t176 * t222 + t177 * t220) * t538 + (t153 * t175 + t154 * t173 - t164 * t174 - t165 * t172) * t537) * t541 + (t369 + (t425 * t364 + t423 * t366) * t511 + t38 * t564 + (t37 + t94) * t512 + (-t426 * t364 + t424 * t366 + t95) * t509 + (t361 / 0.2e1 + t362 / 0.2e1) * t387) * qJD(2), t35 * qJD(1) + t1 * qJD(2) + t369 * qJD(3) + t19 * qJD(4) + ((t376 + (-t285 * t367 + t287 * t365) * t309) * t538 + (-t178 * t183 - t179 * t181 + t444) * t537) * t539, t15 * qJD(2) + t19 * qJD(3) + t414; (t368 - (t318 + t317) * t366 / 0.2e1 - t562) * qJD(1) + t4 * qJD(2) + t3 * qJD(3) + t16 * qJD(4) + (-t515 / 0.4e1 - t506 / 0.4e1 - t507 / 0.4e1) * t542, t4 * qJD(1) + (m(5) * (t114 * t98 - t172 * t173 - t174 * t175) + m(4) * (t112 * t152 - t220 * t552 - t222 * t551) + (t361 * t290 + (-t365 * t289 + t545) * t367) * t511 + (t362 * t289 + (-t367 * t290 + t545) * t365) * t509 + m(3) * ((t365 * (rSges(3,1) * t450 - t416) + t367 * (rSges(3,1) * t449 + t365 * rSges(3,3) - t339)) * (-t365 * t295 - t296 * t367) + t415 * t322 * t320) + t408) * qJD(2) + t6 * qJD(3) + t57 * t498, t3 * qJD(1) + t6 * qJD(2) + t10 * qJD(4) + ((t51 + t43) * t537 + (t61 + t93) * t538) * t539 + (t408 - t549) * qJD(3), t16 * qJD(1) + t10 * qJD(3) + t57 * t500 + t118; t368 * qJD(1) + t2 * qJD(2) + t375 * qJD(3) + t20 * qJD(4) + (-t505 / 0.4e1 - t513 / 0.4e1) * t542, t2 * qJD(1) + t7 + t11 * qJD(4) + ((t114 * t116 - t173 * t181 - t175 * t183 + t43) * t537 + (t155 * t152 + (-t220 * t365 - t222 * t367) * t309 + t61) * t538) * t541 + (t408 + t581) * qJD(2), qJD(1) * t375 + t8 * qJD(2) + t72 * t498 + t7, t20 * qJD(1) + t11 * qJD(2) + t72 * t499 + t118; t14 * qJD(2) + t18 * qJD(3) - t414, t14 * qJD(1) + (-t346 * t114 + (t173 * t365 + t175 * t367) * t345 - t57 + t413) * t500 + t9 * qJD(3) + t171, t18 * qJD(1) + t9 * qJD(2) + (t45 - t72) * t499 + t171, 0.4e1 * (qJD(2) / 0.4e1 + qJD(3) / 0.4e1) * t180 * m(5);];
Cq = t5;
