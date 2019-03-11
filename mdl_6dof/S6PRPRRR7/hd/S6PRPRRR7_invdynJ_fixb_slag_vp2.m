% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_invdynJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:40
% EndTime: 2019-03-08 20:51:41
% DurationCPUTime: 32.31s
% Computational Cost: add. (24514->873), mult. (74732->1271), div. (0->0), fcn. (67064->18), ass. (0->430)
t334 = sin(pkin(7));
t336 = cos(pkin(14));
t492 = t334 * t336;
t332 = sin(pkin(14));
t338 = cos(pkin(7));
t498 = t332 * t338;
t301 = pkin(2) * t498 + qJ(3) * t492;
t337 = cos(pkin(8));
t491 = t334 * t337;
t333 = sin(pkin(8));
t495 = t333 * t338;
t376 = (t336 * t491 + t495) * pkin(10);
t248 = t376 + t301;
t342 = sin(qJ(4));
t346 = cos(qJ(4));
t485 = t336 * t338;
t329 = pkin(2) * t485;
t438 = pkin(10) * t337 + qJ(3);
t499 = t332 * t334;
t534 = pkin(3) * t338;
t254 = -t438 * t499 + t329 + t534;
t533 = pkin(10) * t332;
t394 = -pkin(3) * t336 - t333 * t533;
t274 = (-pkin(2) + t394) * t334;
t397 = t254 * t337 + t274 * t333;
t139 = -t342 * t248 + t397 * t346;
t484 = t337 * t342;
t372 = (-t332 * t484 + t336 * t346) * t334;
t113 = qJD(3) * t372 + qJD(4) * t139;
t335 = sin(pkin(6));
t347 = cos(qJ(2));
t343 = sin(qJ(2));
t482 = t338 * t343;
t276 = (-t332 * t347 - t336 * t482) * t335;
t266 = qJD(1) * t276;
t277 = (-t332 * t482 + t336 * t347) * t335;
t270 = qJD(1) * t277;
t473 = qJD(1) * t335;
t449 = t343 * t473;
t427 = t334 * t449;
t405 = t333 * t427;
t169 = t270 * t346 + (t266 * t337 + t405) * t342;
t606 = -t169 + t113;
t197 = -t254 * t333 + t337 * t274;
t389 = t332 * t346 + t336 * t484;
t494 = t333 * t342;
t251 = t334 * t389 + t338 * t494;
t483 = t337 * t346;
t390 = -t332 * t342 + t336 * t483;
t493 = t333 * t346;
t365 = t334 * t390 + t338 * t493;
t119 = -pkin(4) * t365 - pkin(11) * t251 + t197;
t237 = t346 * t248;
t140 = t254 * t484 + t274 * t494 + t237;
t293 = -t333 * t492 + t337 * t338;
t126 = pkin(11) * t293 + t140;
t244 = t365 * qJD(4);
t245 = t251 * qJD(4);
t496 = t333 * t334;
t424 = qJD(3) * t332 * t496;
t160 = pkin(4) * t245 - pkin(11) * t244 + t424;
t228 = -t266 * t333 + t337 * t427;
t341 = sin(qJ(5));
t345 = cos(qJ(5));
t467 = qJD(5) * t345;
t468 = qJD(5) * t341;
t609 = t119 * t467 - t126 * t468 + t606 * t345 + (t160 - t228) * t341;
t373 = t334 * (t332 * t483 + t336 * t342);
t608 = qJD(3) * t373 + (t342 * t397 + t237) * qJD(4) + t266 * t483 - t270 * t342 + t346 * t405;
t271 = qJD(2) * t372;
t469 = qJD(4) * t346;
t445 = t333 * t469;
t625 = -t271 + t445;
t242 = t365 * qJD(2);
t239 = qJD(5) - t242;
t461 = qJD(2) * qJD(4);
t171 = (qJDD(2) * t342 + t346 * t461) * t495 + (qJDD(2) * t389 + t390 * t461) * t334;
t243 = t251 * qJD(2);
t279 = qJD(2) * t293 + qJD(4);
t195 = -t243 * t341 + t279 * t345;
t278 = qJDD(2) * t293 + qJDD(4);
t100 = qJD(5) * t195 + t171 * t345 + t278 * t341;
t550 = t100 / 0.2e1;
t624 = Ifges(6,4) * t550;
t623 = -pkin(12) * t245 - t609;
t208 = t251 * t341 - t345 * t293;
t149 = -qJD(5) * t208 + t244 * t345;
t209 = t251 * t345 + t293 * t341;
t150 = qJD(5) * t209 + t244 * t341;
t622 = pkin(5) * t150 - pkin(12) * t149 + t608;
t340 = sin(qJ(6));
t344 = cos(qJ(6));
t416 = t340 * mrSges(7,1) + t344 * mrSges(7,2);
t599 = m(6) + m(7);
t570 = -pkin(11) * t599 + mrSges(5,2) - t416;
t196 = t243 * t345 + t279 * t341;
t101 = -qJD(5) * t196 - t171 * t341 + t278 * t345;
t549 = t101 / 0.2e1;
t172 = t334 * (qJDD(2) * t390 - t389 * t461) - (-qJDD(2) * t346 + t342 * t461) * t495;
t167 = qJDD(5) - t172;
t544 = t167 / 0.2e1;
t304 = -t345 * t337 + t341 * t494;
t471 = qJD(2) * t334;
t448 = t332 * t471;
t426 = t333 * t448;
t582 = -qJD(5) * t304 - t341 * t426 + t345 * t625;
t267 = qJD(2) * t373;
t470 = qJD(4) * t342;
t446 = t333 * t470;
t621 = -t267 + t446;
t307 = qJ(3) * t471 + t449;
t312 = qJD(2) * pkin(2) + t347 * t473;
t339 = cos(pkin(6));
t472 = qJD(1) * t339;
t450 = t334 * t472;
t216 = t336 * t307 + t312 * t498 + t332 * t450;
t191 = qJD(2) * t376 + t216;
t462 = t338 * t472 + qJD(3);
t221 = (qJD(2) * t394 - t312) * t334 + t462;
t215 = -t307 * t332 + t312 * t485 + t336 * t450;
t382 = -t491 * t533 + t534;
t192 = qJD(2) * t382 + t215;
t508 = t192 * t337;
t399 = t221 * t333 + t508;
t87 = t191 * t346 + t342 * t399;
t77 = pkin(11) * t279 + t87;
t131 = -t192 * t333 + t337 * t221;
t88 = -pkin(4) * t242 - pkin(11) * t243 + t131;
t40 = t341 * t88 + t345 * t77;
t34 = pkin(12) * t239 + t40;
t86 = -t342 * t191 + t346 * t399;
t76 = -pkin(4) * t279 - t86;
t47 = -pkin(5) * t195 - pkin(12) * t196 + t76;
t15 = -t34 * t340 + t344 * t47;
t16 = t34 * t344 + t340 * t47;
t565 = t15 * mrSges(7,1) - t16 * mrSges(7,2);
t135 = -t196 * t340 + t239 * t344;
t136 = t196 * t344 + t239 * t340;
t194 = qJD(6) - t195;
t56 = t136 * Ifges(7,5) + t135 * Ifges(7,6) + t194 * Ifges(7,3);
t516 = Ifges(6,4) * t196;
t611 = t239 * Ifges(6,6);
t612 = t195 * Ifges(6,2);
t98 = t516 + t611 + t612;
t604 = -t98 / 0.2e1 + t56 / 0.2e1;
t620 = t565 + t604;
t96 = qJDD(6) - t101;
t552 = t96 / 0.2e1;
t46 = -qJD(6) * t136 - t100 * t340 + t167 * t344;
t558 = t46 / 0.2e1;
t45 = qJD(6) * t135 + t100 * t344 + t167 * t340;
t559 = t45 / 0.2e1;
t442 = qJD(2) * t473;
t322 = t347 * t442;
t460 = qJDD(1) * t335;
t292 = t343 * t460 + t322;
t255 = (qJ(3) * qJDD(2) + qJD(2) * qJD(3)) * t334 + t292;
t423 = t343 * t442;
t291 = t347 * t460 - t423;
t280 = qJDD(2) * pkin(2) + t291;
t459 = qJDD(1) * t339;
t441 = t334 * t459;
t180 = -t255 * t332 + t280 * t485 + t336 * t441;
t158 = qJDD(2) * t382 + t180;
t457 = t338 * t459 + qJDD(3);
t205 = (qJDD(2) * t394 - t280) * t334 + t457;
t181 = t336 * t255 + t280 * t498 + t332 * t441;
t607 = qJD(4) * t508 + qJDD(2) * t376 + t181;
t36 = t346 * (t158 * t337 + t205 * t333) - t191 * t469 - t221 * t446 - t607 * t342;
t32 = -pkin(4) * t278 - t36;
t14 = -pkin(5) * t101 - pkin(12) * t100 + t32;
t35 = t158 * t484 - t191 * t470 + t205 * t494 + t221 * t445 + t346 * t607;
t31 = pkin(11) * t278 + t35;
t110 = -t158 * t333 + t337 * t205;
t53 = -pkin(4) * t172 - pkin(11) * t171 + t110;
t7 = t345 * t31 + t341 * t53 + t88 * t467 - t468 * t77;
t5 = pkin(12) * t167 + t7;
t1 = qJD(6) * t15 + t14 * t340 + t344 * t5;
t2 = -qJD(6) * t16 + t14 * t344 - t340 * t5;
t566 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t9 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t96;
t619 = t566 + 0.2e1 * Ifges(6,2) * t549 + 0.2e1 * Ifges(6,6) * t544 + t624 - t9 / 0.2e1 - Ifges(7,5) * t559 - Ifges(7,6) * t558 - Ifges(7,3) * t552;
t173 = pkin(4) * t243 - pkin(11) * t242;
t60 = t341 * t173 + t345 * t86;
t618 = pkin(11) * t468 + pkin(12) * t243 + t60;
t617 = -pkin(11) * qJD(6) * t345 - t87 + t239 * (pkin(5) * t341 - pkin(12) * t345);
t587 = t341 * t119 + t345 * t126;
t55 = -pkin(12) * t365 + t587;
t125 = -pkin(4) * t293 - t139;
t72 = pkin(5) * t208 - pkin(12) * t209 + t125;
t25 = t340 * t72 + t344 * t55;
t616 = -qJD(6) * t25 + t340 * t623 + t344 * t622;
t24 = -t340 * t55 + t344 * t72;
t615 = qJD(6) * t24 + t340 * t622 - t344 * t623;
t614 = mrSges(6,1) * t76;
t613 = t76 * mrSges(6,2);
t489 = t335 * t343;
t452 = t334 * t489;
t233 = -t276 * t333 + t337 * t452;
t479 = t340 * t345;
t161 = -t242 * t479 + t243 * t344;
t605 = t340 * t467 + t161;
t481 = t338 * t347;
t490 = t334 * t339;
t249 = t336 * t489 + (t335 * t481 + t490) * t332;
t388 = -t332 * t343 + t336 * t481;
t364 = t335 * t388 + t336 * t490;
t488 = t335 * t347;
t387 = t334 * t488 - t338 * t339;
t351 = -t387 * t333 + t364 * t337;
t151 = t249 * t342 - t346 * t351;
t417 = -mrSges(7,1) * t344 + mrSges(7,2) * t340;
t379 = m(7) * pkin(5) - t417;
t418 = -mrSges(6,1) * t345 + mrSges(6,2) * t341;
t455 = m(7) * pkin(12) + mrSges(7,3);
t603 = pkin(4) * t599 + t341 * t455 + t345 * t379 + mrSges(5,1) - t418;
t602 = Ifges(6,1) * t550 + Ifges(6,5) * t544;
t560 = Ifges(6,4) * t549 + t602;
t542 = t194 / 0.2e1;
t545 = t136 / 0.2e1;
t547 = t135 / 0.2e1;
t601 = Ifges(7,5) * t545 + Ifges(7,6) * t547 + Ifges(7,3) * t542;
t39 = -t341 * t77 + t345 * t88;
t598 = t39 * mrSges(6,1);
t597 = t40 * mrSges(6,2);
t17 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t70 = mrSges(6,1) * t167 - mrSges(6,3) * t100;
t596 = -t70 + t17;
t595 = mrSges(6,1) + t379;
t431 = mrSges(6,2) - t455;
t518 = mrSges(6,3) * t196;
t138 = mrSges(6,1) * t239 - t518;
t75 = -mrSges(7,1) * t135 + mrSges(7,2) * t136;
t594 = t138 - t75;
t141 = mrSges(5,1) * t278 - mrSges(5,3) * t171;
t48 = -mrSges(6,1) * t101 + mrSges(6,2) * t100;
t593 = t141 - t48;
t318 = -pkin(5) * t345 - pkin(12) * t341 - pkin(4);
t464 = qJD(6) * t344;
t592 = t318 * t464 + t340 * t617 - t344 * t618;
t466 = qJD(6) * t340;
t591 = -t318 * t466 + t340 * t618 + t344 * t617;
t586 = t341 * t464 + t605;
t477 = t344 * t345;
t162 = t242 * t477 + t243 * t340;
t465 = qJD(6) * t341;
t585 = t340 * t465 - t344 * t467 + t162;
t305 = t337 * t341 + t345 * t494;
t391 = -t305 * t344 + t340 * t493;
t584 = qJD(6) * t391 - t340 * t582 + t344 * t621;
t259 = -t305 * t340 - t344 * t493;
t583 = qJD(6) * t259 + t340 * t621 + t344 * t582;
t520 = mrSges(5,3) * t243;
t476 = -mrSges(5,1) * t279 - mrSges(6,1) * t195 + mrSges(6,2) * t196 + t520;
t581 = qJD(5) * t305 + t341 * t625 + t345 * t426;
t509 = sin(pkin(13));
t433 = t509 * t343;
t510 = cos(pkin(13));
t434 = t510 * t347;
t297 = -t339 * t433 + t434;
t432 = t509 * t347;
t435 = t510 * t343;
t368 = t339 * t432 + t435;
t436 = t335 * t509;
t360 = t334 * t436 - t338 * t368;
t353 = t297 * t332 - t336 * t360;
t361 = -t334 * t368 - t338 * t436;
t580 = t361 * t333 + t353 * t337;
t296 = t339 * t435 + t432;
t367 = -t339 * t434 + t433;
t437 = t335 * t510;
t358 = -t334 * t437 - t338 * t367;
t354 = t296 * t332 - t336 * t358;
t359 = -t334 * t367 + t338 * t437;
t579 = t359 * t333 + t354 * t337;
t578 = -mrSges(3,1) * t347 + mrSges(3,2) * t343;
t8 = -qJD(5) * t40 - t31 * t341 + t345 * t53;
t577 = -t341 * t8 + t345 * t7;
t576 = t1 * t344 - t2 * t340;
t575 = -m(5) - t599;
t572 = m(4) - t575;
t571 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t100 + Ifges(6,6) * t101 + Ifges(6,3) * t167;
t29 = -qJD(5) * t587 - t113 * t341 + t160 * t345;
t569 = t36 * mrSges(5,1) - t35 * mrSges(5,2) + Ifges(5,5) * t171 + Ifges(5,6) * t172 + Ifges(5,3) * t278;
t33 = -pkin(5) * t239 - t39;
t568 = -m(7) * t33 + t594;
t567 = -m(6) * t76 - t476;
t543 = -t194 / 0.2e1;
t546 = -t136 / 0.2e1;
t548 = -t135 / 0.2e1;
t563 = Ifges(7,5) * t546 + Ifges(7,6) * t548 + Ifges(7,3) * t543 - t565;
t10 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t96 * Ifges(7,6);
t562 = t10 / 0.2e1;
t561 = Ifges(7,1) * t559 + Ifges(7,4) * t558 + Ifges(7,5) * t552;
t513 = Ifges(7,4) * t136;
t57 = Ifges(7,2) * t135 + Ifges(7,6) * t194 + t513;
t556 = -t57 / 0.2e1;
t555 = t57 / 0.2e1;
t134 = Ifges(7,4) * t135;
t58 = Ifges(7,1) * t136 + Ifges(7,5) * t194 + t134;
t554 = -t58 / 0.2e1;
t553 = t58 / 0.2e1;
t541 = -t195 / 0.2e1;
t540 = -t196 / 0.2e1;
t539 = t196 / 0.2e1;
t538 = -t239 / 0.2e1;
t535 = t243 / 0.2e1;
t6 = -pkin(5) * t167 - t8;
t529 = t341 * t6;
t523 = mrSges(4,1) * t336;
t521 = mrSges(5,3) * t242;
t519 = mrSges(6,3) * t195;
t517 = Ifges(5,4) * t243;
t515 = Ifges(6,4) * t341;
t514 = Ifges(6,4) * t345;
t512 = Ifges(7,4) * t340;
t511 = Ifges(7,4) * t344;
t507 = t195 * t340;
t506 = t195 * t344;
t224 = -t296 * t485 + t332 * t367;
t505 = t224 * t333;
t226 = -t297 * t485 + t332 * t368;
t504 = t226 * t333;
t503 = t242 * t341;
t502 = t242 * t345;
t392 = -mrSges(4,2) * t338 + mrSges(4,3) * t492;
t303 = t392 * qJD(2);
t486 = t336 * t303;
t480 = t340 * t341;
t478 = t341 * t344;
t474 = pkin(2) * t488 + qJ(3) * t452;
t456 = pkin(11) * t467;
t454 = t334 * t493;
t443 = qJDD(2) * t523;
t440 = t467 / 0.2e1;
t439 = -t465 / 0.2e1;
t430 = t333 * t452;
t425 = qJD(2) * t452;
t419 = mrSges(4,2) * t332 - t523;
t415 = Ifges(6,1) * t345 - t515;
t414 = Ifges(7,1) * t344 - t512;
t413 = Ifges(7,1) * t340 + t511;
t412 = -Ifges(6,2) * t341 + t514;
t411 = -Ifges(7,2) * t340 + t511;
t410 = Ifges(7,2) * t344 + t512;
t409 = Ifges(6,5) * t345 - Ifges(6,6) * t341;
t408 = Ifges(7,5) * t344 - Ifges(7,6) * t340;
t407 = Ifges(7,5) * t340 + Ifges(7,6) * t344;
t59 = t173 * t345 - t341 * t86;
t404 = t333 * t425;
t152 = t249 * t346 + t342 * t351;
t352 = -t333 * t364 - t337 * t387;
t103 = t152 * t345 + t341 * t352;
t63 = t103 * t344 + t151 * t340;
t62 = -t103 * t340 + t151 * t344;
t68 = t119 * t345 - t126 * t341;
t155 = t209 * t344 - t340 * t365;
t154 = -t209 * t340 - t344 * t365;
t398 = -t215 * t332 + t216 * t336;
t225 = -t296 * t498 - t336 * t367;
t289 = t367 * pkin(2);
t396 = t225 * pkin(3) - pkin(10) * t505 - t289;
t227 = -t297 * t498 - t336 * t368;
t290 = t368 * pkin(2);
t395 = t227 * pkin(3) - pkin(10) * t504 - t290;
t393 = mrSges(4,1) * t338 - mrSges(4,3) * t499;
t386 = t33 * t416;
t377 = t277 * pkin(3) + pkin(10) * t233 + t474;
t366 = -mrSges(6,3) + t570;
t102 = t152 * t341 - t345 * t352;
t319 = qJDD(2) * mrSges(4,2) * t499;
t302 = t393 * qJD(2);
t300 = -qJ(3) * t499 + t329;
t295 = t392 * qJDD(2);
t294 = t393 * qJDD(2);
t285 = t419 * t471;
t284 = -t334 * t443 + t319;
t282 = pkin(11) * t477 + t318 * t340;
t281 = -pkin(11) * t479 + t318 * t344;
t269 = qJD(2) * t277;
t268 = qJD(2) * t276;
t265 = -t312 * t334 + t462;
t246 = -t280 * t334 + t457;
t238 = Ifges(5,4) * t242;
t229 = -t268 * t333 + t337 * t425;
t207 = t297 * t336 + t332 * t360;
t206 = t296 * t336 + t332 * t358;
t201 = -mrSges(5,2) * t279 + t521;
t193 = Ifges(6,4) * t195;
t190 = t277 * t346 + (t276 * t337 + t430) * t342;
t188 = t297 * t491 - t504;
t187 = t296 * t491 - t505;
t170 = -mrSges(5,1) * t242 + mrSges(5,2) * t243;
t157 = t333 * t353 - t337 * t361;
t156 = t333 * t354 - t337 * t359;
t144 = t243 * Ifges(5,1) + t279 * Ifges(5,5) + t238;
t143 = t242 * Ifges(5,2) + t279 * Ifges(5,6) + t517;
t142 = -mrSges(5,2) * t278 + mrSges(5,3) * t172;
t137 = -mrSges(6,2) * t239 + t519;
t130 = t227 * t346 + (t226 * t337 + t297 * t496) * t342;
t128 = t225 * t346 + (t224 * t337 + t296 * t496) * t342;
t120 = t169 * t341 - t345 * t228;
t117 = pkin(5) * t196 - pkin(12) * t195;
t111 = -mrSges(5,1) * t172 + mrSges(5,2) * t171;
t109 = t207 * t346 - t342 * t580;
t108 = t207 * t342 + t346 * t580;
t107 = t206 * t346 - t342 * t579;
t106 = t206 * t342 + t346 * t579;
t99 = t196 * Ifges(6,1) + t239 * Ifges(6,5) + t193;
t97 = t196 * Ifges(6,5) + t195 * Ifges(6,6) + t239 * Ifges(6,3);
t92 = t269 * t346 + (t268 * t337 + t404) * t342 - t151 * qJD(4);
t91 = qJD(4) * t152 - t268 * t483 + t269 * t342 - t346 * t404;
t90 = mrSges(7,1) * t194 - mrSges(7,3) * t136;
t89 = -mrSges(7,2) * t194 + mrSges(7,3) * t135;
t74 = -qJD(6) * t155 - t149 * t340 + t245 * t344;
t73 = qJD(6) * t154 + t149 * t344 + t245 * t340;
t71 = -mrSges(6,2) * t167 + mrSges(6,3) * t101;
t67 = t109 * t345 + t157 * t341;
t65 = t107 * t345 + t156 * t341;
t54 = pkin(5) * t365 - t68;
t50 = -pkin(5) * t243 - t59;
t42 = -qJD(5) * t102 + t229 * t341 + t92 * t345;
t27 = -pkin(5) * t245 - t29;
t23 = -mrSges(7,2) * t96 + mrSges(7,3) * t46;
t22 = mrSges(7,1) * t96 - mrSges(7,3) * t45;
t21 = t117 * t340 + t344 * t39;
t20 = t117 * t344 - t340 * t39;
t13 = qJD(6) * t62 + t340 * t91 + t344 * t42;
t12 = -qJD(6) * t63 - t340 * t42 + t344 * t91;
t3 = [(m(3) * t339 ^ 2 + m(2)) * qJDD(1) + m(7) * (t1 * t63 + t12 * t15 + t13 * t16 + t2 * t62) + t352 * t111 + m(4) * (t181 * t249 + t215 * t268 + t216 * t269 + (t180 * t492 + t246 * t338) * t339) + (-t578 * qJDD(2) + m(3) * (t291 * t347 + t292 * t343) + (-mrSges(3,1) * t343 - mrSges(3,2) * t347) * qJD(2) ^ 2 + m(4) * (t180 * t388 + (qJD(2) * t265 * t343 - t246 * t347) * t334)) * t335 + m(6) * (t103 * t7 + t40 * t42) + m(5) * (t110 * t352 + t131 * t229 + t35 * t152 + t87 * t92) + (-m(6) * t39 - t568) * (qJD(5) * t103 - t229 * t345 + t92 * t341) + t364 * t294 + (-m(5) * t86 - t567) * t91 + (-m(2) - m(3) - t572) * g(3) + t268 * t302 + t269 * t303 + t249 * t295 + t229 * t170 + t92 * t201 + t152 * t142 + t42 * t137 + t103 * t71 + t13 * t89 + t12 * t90 + (-m(5) * t36 + m(6) * t32 - t593) * t151 + (-m(6) * t8 + m(7) * t6 + t596) * t102 + t62 * t22 + t63 * t23 - t387 * t284 + t285 * t425; (t1 * t154 - t15 * t73 - t155 * t2 + t16 * t74) * mrSges(7,3) - t149 * t39 * mrSges(6,3) + t615 * t89 + (t1 * t25 + t2 * t24 + t54 * t6 + (-t120 + t27) * t33 + t615 * t16 + t616 * t15) * m(7) + t616 * t90 + t587 * t71 + (Ifges(7,4) * t155 + Ifges(7,2) * t154) * t558 + (Ifges(7,4) * t73 + Ifges(7,2) * t74) * t547 + (-t40 * mrSges(6,3) - t611 / 0.2e1 - t612 / 0.2e1 + t614 - Ifges(6,4) * t539 + t601 + t620) * t150 + (-t244 * t86 - t245 * t87) * mrSges(5,3) + t609 * t137 + (t125 * t32 + t68 * t8 + t587 * t7 + t608 * t76 + t609 * t40 + (t120 + t29) * t39) * m(6) + (Ifges(4,3) * t338 ^ 2 + Ifges(3,3) + ((Ifges(4,2) * t336 ^ 2 + (Ifges(4,1) * t332 + 0.2e1 * Ifges(4,4) * t336) * t332) * t334 + 0.2e1 * (Ifges(4,5) * t332 + Ifges(4,6) * t336) * t338) * t334) * qJDD(2) + (Ifges(7,1) * t155 + Ifges(7,4) * t154) * t559 + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t545 + (Ifges(6,1) * t149 + Ifges(6,5) * t245) * t539 - (t110 * mrSges(5,1) - t35 * mrSges(5,3) - Ifges(5,4) * t171 - Ifges(5,2) * t172 - Ifges(5,6) * t278 + t571) * t365 + (mrSges(5,2) * t110 - mrSges(5,3) * t36 + Ifges(5,1) * t171 + Ifges(5,4) * t172 + Ifges(5,5) * t278) * t251 + (mrSges(4,1) * t180 - mrSges(4,2) * t181) * t338 + (-m(4) * t474 - t277 * mrSges(4,1) - t276 * mrSges(4,2) - t190 * mrSges(5,1) - t233 * mrSges(5,3) + t578 * t335 - t595 * (t190 * t345 + t233 * t341) + t366 * (-t276 * t483 + t277 * t342 - t346 * t430) + t431 * (t190 * t341 - t233 * t345) + t599 * (-t190 * pkin(4) - t377)) * g(3) + (m(4) * t289 - t225 * mrSges(4,1) - t224 * mrSges(4,2) - t128 * mrSges(5,1) - t187 * mrSges(5,3) + t367 * mrSges(3,1) + t296 * mrSges(3,2) + t431 * (t128 * t341 - t187 * t345) - t595 * (t128 * t345 + t187 * t341) + t366 * (-t224 * t483 + t225 * t342 - t296 * t454) + t599 * (-t128 * pkin(4) - t396)) * g(2) + (m(4) * t290 - t227 * mrSges(4,1) - t226 * mrSges(4,2) - t130 * mrSges(5,1) - t188 * mrSges(5,3) + t368 * mrSges(3,1) + t297 * mrSges(3,2) + t431 * (t130 * t341 - t188 * t345) - t595 * (t130 * t345 + t188 * t341) + t366 * (-t226 * t483 + t227 * t342 - t297 * t454) + t599 * (-t130 * pkin(4) - t395)) * g(1) + t245 * t598 + (Ifges(7,5) * t155 + Ifges(7,6) * t154) * t552 + (Ifges(7,5) * t73 + Ifges(7,6) * t74) * t542 + t606 * t201 + t608 * t476 + (-g(1) * t395 - g(2) * t396 - g(3) * t377 + t110 * t197 + t139 * t36 + t140 * t35 + t606 * t87 - t608 * t86 + (-t228 + t424) * t131) * m(5) + (mrSges(6,2) * t32 - mrSges(6,3) * t8 + 0.2e1 * t560) * t209 + t569 * t293 + (t246 * t419 - pkin(2) * t284 - t285 * t449 + (t486 + (t170 * t333 - t302) * t332) * qJD(3) + (-g(3) * t489 - t180 * t332 + t181 * t336) * mrSges(4,3) + (g(1) * t297 + g(2) * t296) * (-m(4) * qJ(3) + t438 * t575 - mrSges(4,3))) * t334 + t239 * (Ifges(6,5) * t149 + Ifges(6,3) * t245) / 0.2e1 + t300 * t294 + t301 * t295 - t266 * t302 - t270 * t303 + t279 * (Ifges(5,5) * t244 - Ifges(5,6) * t245) / 0.2e1 + t244 * t144 / 0.2e1 + t245 * t97 / 0.2e1 + t131 * (mrSges(5,1) * t245 + mrSges(5,2) * t244) + t242 * (Ifges(5,4) * t244 - Ifges(5,2) * t245) / 0.2e1 - t245 * t143 / 0.2e1 - t228 * t170 + t197 * t111 + t195 * (Ifges(6,4) * t149 + Ifges(6,6) * t245) / 0.2e1 + t6 * (-mrSges(7,1) * t154 + mrSges(7,2) * t155) + t149 * t99 / 0.2e1 + t29 * t138 + t139 * t141 + t140 * t142 + t125 * t48 - m(4) * (t215 * t266 + t216 * t270 + t265 * t427) + t25 * t23 + t27 * t75 + t24 * t22 + (t291 + t423) * mrSges(3,1) + (-t292 + t322) * mrSges(3,2) + t155 * t561 + t154 * t562 + t73 * t553 + t74 * t555 + (Ifges(5,1) * t244 - Ifges(5,4) * t245) * t535 + t149 * t613 + m(4) * (t180 * t300 + t181 * t301 + (-pkin(2) * t246 + qJD(3) * t398) * t334) + t594 * t120 - t245 * t597 + t54 * t17 + t68 * t70 + t33 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + (t32 * mrSges(6,1) - t7 * mrSges(6,3) - t619 - t624) * t208; (-t170 * t448 + t342 * t142 + t593 * t346 + (t201 * t346 + t342 * t476) * qJD(4)) * t333 - t476 * t267 + t319 + t584 * t90 + (-t443 + (t302 * t332 - t486) * qJD(2)) * t334 + t583 * t89 + t596 * t304 + t337 * t111 + t305 * t71 - t271 * t201 + t259 * t22 - t391 * t23 + t582 * t137 - t594 * t581 + (-t1 * t391 + t15 * t584 + t16 * t583 + t2 * t259 + t304 * t6 + t33 * t581) * m(7) + (-t304 * t8 + t305 * t7 + (-t32 * t346 + t470 * t76) * t333 - t267 * t76 + t582 * t40 - t581 * t39) * m(6) + (t110 * t337 + (t342 * t35 + t346 * t36 + (-t342 * t86 + t346 * t87) * qJD(4)) * t333 - t131 * t426 + t267 * t86 - t271 * t87) * m(5) + (-t398 * t471 + t246) * m(4) + (g(1) * t361 + g(2) * t359 + g(3) * t387) * t572; t591 * t90 + (t439 * t57 + t440 * t58) * t344 + (t151 * t603 + t152 * t570) * g(3) + (t106 * t603 + t107 * t570) * g(2) + (t108 * t603 + t109 * t570) * g(1) + t605 * t556 + (qJD(5) * t601 + t408 * t552 + t411 * t558 + t414 * t559 + t560 + t602) * t341 + t569 + (pkin(11) * t71 + (t408 * t542 + t411 * t547 + t414 * t545) * qJD(5) + t619) * t345 + (-pkin(11) * t137 + t620) * t468 + t239 * t76 * (mrSges(6,1) * t341 + mrSges(6,2) * t345) + (Ifges(7,5) * t162 + Ifges(7,6) * t161) * t543 + (t563 - t604) * t503 + (mrSges(7,1) * t586 - mrSges(7,2) * t585) * t33 + (-t1 * t480 + t15 * t585 - t16 * t586 - t2 * t478) * mrSges(7,3) + (-g(1) * t109 - g(2) * t107 - g(3) * t152 + (-t468 + t503) * t40 + (-t467 + t502) * t39 + t577) * mrSges(6,3) + (-t39 * t59 - t40 * t60 - pkin(4) * t32 + ((-t341 * t40 - t345 * t39) * qJD(5) + t577) * pkin(11)) * m(6) + (-t502 / 0.2e1 + t440) * t99 + (-t201 + t521) * t86 + t515 * t549 + (-t456 - t59) * t138 + (t456 - t50) * t75 + (Ifges(7,4) * t162 + Ifges(7,2) * t161) * t548 + t514 * t550 + t243 * t597 + t596 * pkin(11) * t341 + (-t407 * t542 - t410 * t547 - t413 * t545) * t465 + (Ifges(7,1) * t162 + Ifges(7,4) * t161) * t546 + (t195 * t412 + t196 * t415 + t239 * t409) * qJD(5) / 0.2e1 - (Ifges(5,1) * t242 - t517 + t97) * t243 / 0.2e1 - (-Ifges(5,2) * t243 + t144 + t238) * t242 / 0.2e1 - t10 * t480 / 0.2e1 + (t520 + t567) * t87 + t282 * t23 - t279 * (Ifges(5,5) * t242 - Ifges(5,6) * t243) / 0.2e1 + t281 * t22 - t131 * (mrSges(5,1) * t243 + mrSges(5,2) * t242) - t60 * t137 + t32 * t418 + t478 * t561 + t162 * t554 + (Ifges(6,6) * t243 + t242 * t412) * t541 + t143 * t535 + (Ifges(6,3) * t243 + t242 * t409) * t538 + (Ifges(6,5) * t243 + t242 * t415) * t540 + t416 * t529 + (-t33 * t50 + t1 * t282 + t2 * t281 + (t33 * t467 + t529) * pkin(11) + t592 * t16 + t591 * t15) * m(7) + t592 * t89 + t340 * t58 * t439 - t243 * t598 - pkin(4) * t48; (Ifges(6,1) * t540 + Ifges(6,5) * t538 + t408 * t543 + t411 * t548 + t414 * t546 - t386 - t613) * t195 + (-t137 + t519) * t39 + (-t516 + t56) * t540 + (t193 + t99) * t541 + (-pkin(5) * t6 - t15 * t20 - t16 * t21) * m(7) + (-Ifges(6,2) * t541 - Ifges(6,6) * t538 + t563 - t614) * t196 + ((-t466 + t507) * t16 + (-t464 + t506) * t15 + t576) * mrSges(7,3) + (-t90 * t464 - t89 * t466 + t344 * t23 - t340 * t22 + m(7) * ((-t15 * t344 - t16 * t340) * qJD(6) + t576)) * pkin(12) + qJD(6) * t386 + (t518 + t568) * t40 + (t135 * t411 + t136 * t414 + t194 * t408) * qJD(6) / 0.2e1 + t571 - t21 * t89 - t20 * t90 - pkin(5) * t17 + t6 * t417 + t413 * t559 + t340 * t561 + t344 * t562 + t407 * t552 + t464 * t553 + t506 * t554 + t507 * t555 + t466 * t556 + t410 * t558 + t98 * t539 + (t102 * t595 + t103 * t431) * g(3) + (t431 * t65 - t595 * (-t107 * t341 + t156 * t345)) * g(2) + (t431 * t67 - t595 * (-t109 * t341 + t157 * t345)) * g(1); -t33 * (mrSges(7,1) * t136 + mrSges(7,2) * t135) - t15 * t89 + t16 * t90 + (Ifges(7,1) * t135 - t513) * t546 + t57 * t545 + (Ifges(7,5) * t135 - Ifges(7,6) * t136) * t543 - g(1) * ((t108 * t344 - t340 * t67) * mrSges(7,1) + (-t108 * t340 - t344 * t67) * mrSges(7,2)) - g(2) * ((t106 * t344 - t340 * t65) * mrSges(7,1) + (-t106 * t340 - t344 * t65) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t62 - mrSges(7,2) * t63) + (t135 * t15 + t136 * t16) * mrSges(7,3) + t9 + (-Ifges(7,2) * t136 + t134 + t58) * t548 - t566;];
tau  = t3;
