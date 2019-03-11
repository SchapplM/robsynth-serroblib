% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:59
% EndTime: 2019-03-09 00:48:59
% DurationCPUTime: 34.02s
% Computational Cost: add. (17295->909), mult. (38877->1252), div. (0->0), fcn. (29596->18), ass. (0->403)
t345 = sin(qJ(3));
t350 = cos(qJ(3));
t396 = pkin(3) * t345 - pkin(9) * t350;
t290 = t396 * qJD(3);
t298 = -pkin(3) * t350 - pkin(9) * t345 - pkin(2);
t344 = sin(qJ(4));
t346 = sin(qJ(2));
t349 = cos(qJ(4));
t437 = qJD(4) * t349;
t439 = qJD(4) * t344;
t441 = qJD(3) * t349;
t340 = sin(pkin(6));
t448 = qJD(1) * t340;
t351 = cos(qJ(2));
t457 = t350 * t351;
t566 = -(t344 * t346 + t349 * t457) * t448 + t344 * t290 + t298 * t437 + (-t345 * t441 - t350 * t439) * pkin(8);
t442 = qJD(3) * t345;
t426 = pkin(8) * t442;
t610 = t349 * t290 + t344 * t426 - (-t344 * t457 + t346 * t349) * t448;
t458 = t349 * t350;
t323 = pkin(8) * t458;
t377 = pkin(4) * t345 - pkin(10) * t458;
t609 = t377 * qJD(3) + (-t323 + (pkin(10) * t345 - t298) * t344) * qJD(4) + t610;
t440 = qJD(3) * t350;
t363 = t344 * t440 + t345 * t437;
t608 = -pkin(10) * t363 + t566;
t420 = t346 * t448;
t294 = qJD(2) * pkin(8) + t420;
t341 = cos(pkin(6));
t447 = qJD(1) * t341;
t227 = -t345 * t294 + t350 * t447;
t287 = t396 * qJD(2);
t167 = -t227 * t344 + t349 * t287;
t352 = -pkin(10) - pkin(9);
t421 = qJD(4) * t352;
t607 = -qJD(2) * t377 + t349 * t421 - t167;
t168 = t349 * t227 + t344 * t287;
t443 = qJD(2) * t350;
t416 = t344 * t443;
t606 = -pkin(10) * t416 - t344 * t421 + t168;
t343 = sin(qJ(5));
t348 = cos(qJ(5));
t378 = t343 * t344 - t348 * t349;
t557 = qJD(4) + qJD(5);
t201 = t557 * t378;
t369 = t378 * t350;
t236 = qJD(2) * t369;
t605 = t201 - t236;
t281 = t343 * t349 + t344 * t348;
t202 = t557 * t281;
t370 = t350 * t281;
t235 = qJD(2) * t370;
t597 = t202 - t235;
t395 = mrSges(4,1) * t350 - mrSges(4,2) * t345;
t553 = m(7) * (-pkin(11) + t352) - mrSges(7,3) + m(6) * t352 - mrSges(6,3) - m(5) * pkin(9) - mrSges(5,3);
t338 = qJ(4) + qJ(5);
t334 = cos(t338);
t495 = pkin(4) * t349;
t297 = pkin(5) * t334 + t495;
t335 = qJ(6) + t338;
t327 = sin(t335);
t328 = cos(t335);
t330 = pkin(3) + t495;
t333 = sin(t338);
t394 = -mrSges(5,1) * t349 + mrSges(5,2) * t344;
t595 = m(7) * (pkin(3) + t297) + mrSges(7,1) * t328 - mrSges(7,2) * t327 + m(6) * t330 + mrSges(6,1) * t334 - mrSges(6,2) * t333 + m(5) * pkin(3) - t394;
t604 = t345 * t553 - t350 * t595 - mrSges(3,1) - t395;
t305 = t352 * t344;
t306 = t352 * t349;
t435 = qJD(5) * t348;
t436 = qJD(5) * t343;
t577 = t305 * t435 + t306 * t436 + t343 * t607 - t348 * t606;
t220 = t343 * t305 - t348 * t306;
t576 = -qJD(5) * t220 + t343 * t606 + t348 * t607;
t277 = t349 * t298;
t460 = t345 * t349;
t197 = -pkin(10) * t460 + t277 + (-pkin(8) * t344 - pkin(4)) * t350;
t234 = t344 * t298 + t323;
t462 = t344 * t345;
t212 = -pkin(10) * t462 + t234;
t121 = t343 * t197 + t348 * t212;
t575 = -qJD(5) * t121 - t608 * t343 + t348 * t609;
t574 = t197 * t435 - t212 * t436 + t343 * t609 + t608 * t348;
t445 = qJD(2) * t345;
t602 = -pkin(5) * t445 + pkin(11) * t605 + t576;
t601 = pkin(11) * t597 - t577;
t133 = -qJD(3) * t369 - t202 * t345;
t600 = pkin(5) * t442 - pkin(11) * t133 + t575;
t240 = t378 * t345;
t134 = -qJD(3) * t370 + t240 * t557;
t599 = -pkin(11) * t134 - t574;
t555 = -m(7) - m(5) - m(4) - m(6);
t216 = -qJD(3) * pkin(3) - t227;
t278 = -t344 * t445 + t441;
t169 = -pkin(4) * t278 + t216;
t279 = qJD(3) * t344 + t349 * t445;
t398 = t348 * t278 - t279 * t343;
t101 = -pkin(5) * t398 + t169;
t194 = t278 * t343 + t279 * t348;
t342 = sin(qJ(6));
t347 = cos(qJ(6));
t113 = t194 * t347 + t342 * t398;
t432 = qJD(2) * qJD(3);
t291 = qJDD(2) * t350 - t345 * t432;
t275 = qJDD(4) - t291;
t269 = qJDD(5) + t275;
t250 = qJDD(6) + t269;
t399 = -t194 * t342 + t347 * t398;
t292 = qJDD(2) * t345 + t350 * t432;
t184 = qJD(4) * t278 + qJDD(3) * t344 + t292 * t349;
t185 = -qJD(4) * t279 + qJDD(3) * t349 - t292 * t344;
t78 = qJD(5) * t398 + t184 * t348 + t185 * t343;
t79 = -qJD(5) * t194 - t184 * t343 + t185 * t348;
t26 = qJD(6) * t399 + t342 * t79 + t347 * t78;
t27 = -qJD(6) * t113 - t342 * t78 + t347 * t79;
t429 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t250;
t321 = qJD(4) - t443;
t311 = qJD(5) + t321;
t593 = pkin(11) * t194;
t228 = t350 * t294 + t345 * t447;
t217 = qJD(3) * pkin(9) + t228;
t419 = t351 * t448;
t230 = qJD(2) * t298 - t419;
t147 = t217 * t349 + t230 * t344;
t117 = pkin(10) * t278 + t147;
t105 = t343 * t117;
t146 = -t217 * t344 + t349 * t230;
t116 = -pkin(10) * t279 + t146;
t99 = pkin(4) * t321 + t116;
t61 = t348 * t99 - t105;
t48 = t61 - t593;
t47 = pkin(5) * t311 + t48;
t582 = pkin(11) * t398;
t107 = t348 * t117;
t62 = t343 * t99 + t107;
t49 = t62 + t582;
t479 = t342 * t49;
t16 = t347 * t47 - t479;
t430 = qJDD(1) * t341;
t446 = qJD(2) * t340;
t410 = qJD(1) * t446;
t308 = t351 * t410;
t431 = qJDD(1) * t340;
t244 = t346 * t431 + t308;
t589 = qJDD(2) * pkin(8) + qJD(3) * t447 + t244;
t143 = -t294 * t442 + t345 * t430 + t350 * t589;
t131 = qJDD(3) * pkin(9) + t143;
t307 = t346 * t410;
t243 = t351 * t431 - t307;
t231 = -qJDD(2) * pkin(2) - t243;
t160 = -pkin(3) * t291 - pkin(9) * t292 + t231;
t55 = -qJD(4) * t147 - t131 * t344 + t349 * t160;
t39 = pkin(4) * t275 - pkin(10) * t184 + t55;
t54 = t349 * t131 + t344 * t160 - t217 * t439 + t230 * t437;
t44 = pkin(10) * t185 + t54;
t13 = -qJD(5) * t62 - t343 * t44 + t348 * t39;
t6 = pkin(5) * t269 - pkin(11) * t78 + t13;
t12 = -t117 * t436 + t343 * t39 + t348 * t44 + t99 * t435;
t7 = pkin(11) * t79 + t12;
t2 = qJD(6) * t16 + t342 * t6 + t347 * t7;
t477 = t347 * t49;
t17 = t342 * t47 + t477;
t3 = -qJD(6) * t17 - t342 * t7 + t347 * t6;
t549 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t376 = t429 - t549;
t428 = Ifges(6,5) * t78 + Ifges(6,6) * t79 + Ifges(6,3) * t269;
t493 = t17 * mrSges(7,3);
t499 = mrSges(7,3) * t16;
t299 = qJD(6) + t311;
t506 = -t299 / 0.2e1;
t519 = -t113 / 0.2e1;
t521 = -t399 / 0.2e1;
t102 = Ifges(7,4) * t399;
t60 = Ifges(7,1) * t113 + Ifges(7,5) * t299 + t102;
t530 = -t60 / 0.2e1;
t483 = Ifges(7,4) * t113;
t59 = Ifges(7,2) * t399 + Ifges(7,6) * t299 + t483;
t532 = -t59 / 0.2e1;
t546 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t596 = t376 + t428 - t546 + (-t101 * mrSges(7,2) + Ifges(7,1) * t519 + Ifges(7,4) * t521 + Ifges(7,5) * t506 + t499 + t530) * t399 - (t101 * mrSges(7,1) + Ifges(7,4) * t519 + Ifges(7,2) * t521 + Ifges(7,6) * t506 - t493 + t532) * t113;
t497 = pkin(4) * t344;
t296 = pkin(5) * t333 + t497;
t393 = t344 * mrSges(5,1) + t349 * mrSges(5,2);
t594 = -m(6) * t497 - m(7) * t296 - t333 * mrSges(6,1) - t327 * mrSges(7,1) - t334 * mrSges(6,2) - t328 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - t393;
t494 = pkin(5) * t194;
t118 = -mrSges(6,1) * t398 + mrSges(6,2) * t194;
t63 = -mrSges(7,1) * t399 + mrSges(7,2) * t113;
t592 = t118 + t63;
t568 = -t228 + (-t416 + t439) * pkin(4);
t588 = -pkin(2) * t555 - t604;
t536 = t26 / 0.2e1;
t535 = t27 / 0.2e1;
t528 = t78 / 0.2e1;
t527 = t79 / 0.2e1;
t517 = t184 / 0.2e1;
t516 = t185 / 0.2e1;
t511 = t250 / 0.2e1;
t510 = t269 / 0.2e1;
t509 = t275 / 0.2e1;
t586 = t291 / 0.2e1;
t585 = t292 / 0.2e1;
t120 = t348 * t197 - t212 * t343;
t85 = -pkin(5) * t350 + pkin(11) * t240 + t120;
t239 = t281 * t345;
t92 = -pkin(11) * t239 + t121;
t51 = t342 * t85 + t347 * t92;
t584 = -qJD(6) * t51 + t342 * t599 + t347 * t600;
t50 = -t342 * t92 + t347 * t85;
t583 = qJD(6) * t50 + t342 * t600 - t347 * t599;
t219 = t348 * t305 + t306 * t343;
t170 = -pkin(11) * t281 + t219;
t171 = -pkin(11) * t378 + t220;
t90 = t170 * t347 - t171 * t342;
t580 = qJD(6) * t90 + t342 * t602 - t601 * t347;
t91 = t170 * t342 + t171 * t347;
t579 = -qJD(6) * t91 + t601 * t342 + t347 * t602;
t496 = pkin(4) * t348;
t329 = pkin(5) + t496;
t433 = qJD(6) * t347;
t434 = qJD(6) * t342;
t464 = t342 * t343;
t64 = -t116 * t343 - t107;
t52 = t64 - t582;
t65 = t348 * t116 - t105;
t53 = t65 - t593;
t573 = -t342 * t52 - t347 * t53 + t329 * t433 + (-t343 * t434 + (t347 * t348 - t464) * qJD(5)) * pkin(4);
t463 = t343 * t347;
t572 = t342 * t53 - t347 * t52 - t329 * t434 + (-t343 * t433 + (-t342 * t348 - t463) * qJD(5)) * pkin(4);
t538 = m(6) * pkin(4);
t571 = -t538 - mrSges(5,1);
t569 = pkin(5) * t597 + t568;
t567 = -qJD(4) * t234 + t610;
t100 = -mrSges(5,1) * t185 + mrSges(5,2) * t184;
t565 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t292 + t100;
t424 = mrSges(4,3) * t445;
t564 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t278 + mrSges(5,2) * t279 + t424;
t466 = t340 * t350;
t258 = t341 * t345 + t346 * t466;
t465 = t340 * t351;
t375 = -t258 * t333 - t334 * t465;
t452 = (-t258 * t327 - t328 * t465) * mrSges(7,1) + (-t258 * t328 + t327 * t465) * mrSges(7,2);
t563 = -t375 * mrSges(6,1) - (-t258 * t334 + t333 * t465) * mrSges(6,2) - t452;
t475 = cos(pkin(12));
t401 = t475 * t351;
t339 = sin(pkin(12));
t469 = t339 * t346;
t254 = -t341 * t469 + t401;
t207 = t339 * t340 * t345 + t254 * t350;
t402 = t475 * t346;
t468 = t339 * t351;
t253 = t341 * t468 + t402;
t380 = -t207 * t333 + t253 * t334;
t455 = (-t207 * t327 + t253 * t328) * mrSges(7,1) + (-t207 * t328 - t253 * t327) * mrSges(7,2);
t562 = -t380 * mrSges(6,1) - (-t207 * t334 - t253 * t333) * mrSges(6,2) - t455;
t252 = t341 * t402 + t468;
t403 = t340 * t475;
t205 = t252 * t350 - t345 * t403;
t251 = -t341 * t401 + t469;
t382 = -t205 * t333 + t251 * t334;
t456 = (-t205 * t327 + t251 * t328) * mrSges(7,1) + (-t205 * t328 - t251 * t327) * mrSges(7,2);
t561 = -t382 * mrSges(6,1) - (-t205 * t334 - t251 * t333) * mrSges(6,2) - t456;
t273 = Ifges(5,4) * t278;
t174 = t279 * Ifges(5,1) + t321 * Ifges(5,5) + t273;
t331 = Ifges(4,4) * t443;
t560 = Ifges(4,1) * t445 + Ifges(4,5) * qJD(3) + t349 * t174 + t331;
t144 = -t294 * t440 - t345 * t589 + t350 * t430;
t559 = t143 * t350 - t144 * t345;
t558 = -t344 * t55 + t349 * t54;
t556 = t279 * Ifges(5,5) + t194 * Ifges(6,5) + t113 * Ifges(7,5) + t278 * Ifges(5,6) + Ifges(6,6) * t398 + Ifges(7,6) * t399 + t321 * Ifges(5,3) + t311 * Ifges(6,3) + t299 * Ifges(7,3);
t552 = mrSges(4,1) + t595;
t551 = mrSges(4,2) + t553;
t548 = -m(5) * t216 - t564;
t547 = -t55 * mrSges(5,1) + t54 * mrSges(5,2);
t544 = pkin(8) * t555 + t594;
t488 = Ifges(4,4) * t345;
t388 = t350 * Ifges(4,2) + t488;
t543 = t17 * mrSges(7,2) + t62 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(2) * t388 / 0.2e1 - t16 * mrSges(7,1) - t61 * mrSges(6,1);
t353 = qJD(2) ^ 2;
t540 = Ifges(7,4) * t536 + Ifges(7,2) * t535 + Ifges(7,6) * t511;
t539 = Ifges(7,1) * t536 + Ifges(7,4) * t535 + Ifges(7,5) * t511;
t537 = m(7) * pkin(5);
t534 = Ifges(6,4) * t528 + Ifges(6,2) * t527 + Ifges(6,6) * t510;
t533 = Ifges(6,1) * t528 + Ifges(6,4) * t527 + Ifges(6,5) * t510;
t531 = t59 / 0.2e1;
t529 = t60 / 0.2e1;
t526 = Ifges(5,1) * t517 + Ifges(5,4) * t516 + Ifges(5,5) * t509;
t484 = Ifges(6,4) * t194;
t97 = Ifges(6,2) * t398 + t311 * Ifges(6,6) + t484;
t525 = -t97 / 0.2e1;
t524 = t97 / 0.2e1;
t188 = Ifges(6,4) * t398;
t98 = t194 * Ifges(6,1) + t311 * Ifges(6,5) + t188;
t523 = -t98 / 0.2e1;
t522 = t98 / 0.2e1;
t520 = t399 / 0.2e1;
t518 = t113 / 0.2e1;
t515 = -t398 / 0.2e1;
t514 = t398 / 0.2e1;
t513 = -t194 / 0.2e1;
t512 = t194 / 0.2e1;
t507 = t279 / 0.2e1;
t505 = t299 / 0.2e1;
t504 = -t311 / 0.2e1;
t503 = t311 / 0.2e1;
t501 = mrSges(6,3) * t61;
t500 = mrSges(6,3) * t62;
t498 = pkin(4) * t279;
t336 = t345 * pkin(8);
t490 = mrSges(6,3) * t398;
t489 = mrSges(6,3) * t194;
t487 = Ifges(4,4) * t350;
t486 = Ifges(5,4) * t344;
t485 = Ifges(5,4) * t349;
t482 = t146 * mrSges(5,3);
t481 = t147 * mrSges(5,3);
t480 = t279 * Ifges(5,4);
t132 = -qJDD(3) * pkin(3) - t144;
t474 = t132 * t345;
t467 = t340 * t346;
t461 = t344 * t350;
t293 = pkin(4) * t462 + t336;
t444 = qJD(2) * t346;
t438 = qJD(4) * t345;
t332 = pkin(8) * t440;
t423 = mrSges(4,3) * t443;
t422 = Ifges(5,5) * t184 + Ifges(5,6) * t185 + Ifges(5,3) * t275;
t229 = pkin(4) * t363 + t332;
t418 = t340 * t444;
t417 = t351 * t446;
t173 = t278 * Ifges(5,2) + t321 * Ifges(5,6) + t480;
t413 = -t344 * t173 / 0.2e1;
t400 = t432 / 0.2e1;
t390 = Ifges(5,1) * t349 - t486;
t389 = Ifges(5,1) * t344 + t485;
t387 = -Ifges(5,2) * t344 + t485;
t386 = Ifges(5,2) * t349 + t486;
t385 = Ifges(4,5) * t350 - Ifges(4,6) * t345;
t384 = Ifges(5,5) * t349 - Ifges(5,6) * t344;
t383 = Ifges(5,5) * t344 + Ifges(5,6) * t349;
t210 = -t258 * t344 - t349 * t465;
t374 = -t258 * t349 + t344 * t465;
t123 = t210 * t348 + t343 * t374;
t124 = t210 * t343 - t348 * t374;
t67 = t123 * t347 - t124 * t342;
t68 = t123 * t342 + t124 * t347;
t165 = -t239 * t347 + t240 * t342;
t166 = -t239 * t342 - t240 * t347;
t195 = -t281 * t342 - t347 * t378;
t196 = t281 * t347 - t342 * t378;
t257 = -t341 * t350 + t345 * t467;
t373 = t216 * t393;
t295 = -qJD(2) * pkin(2) - t419;
t372 = t295 * (mrSges(4,1) * t345 + mrSges(4,2) * t350);
t371 = t345 * (Ifges(4,1) * t350 - t488);
t364 = -t344 * t438 + t349 * t440;
t361 = Ifges(5,5) * t345 + t350 * t390;
t360 = Ifges(5,6) * t345 + t350 * t387;
t359 = Ifges(5,3) * t345 + t350 * t384;
t82 = -pkin(4) * t185 + t132;
t301 = -qJD(3) * mrSges(4,2) + t423;
t285 = t395 * qJD(2);
t248 = pkin(4) * t463 + t329 * t342;
t247 = -pkin(4) * t464 + t329 * t347;
t245 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t291;
t237 = pkin(5) * t378 - t330;
t233 = -pkin(8) * t461 + t277;
t226 = mrSges(5,1) * t321 - mrSges(5,3) * t279;
t225 = -mrSges(5,2) * t321 + mrSges(5,3) * t278;
t213 = -mrSges(4,1) * t291 + mrSges(4,2) * t292;
t209 = -qJD(3) * t257 + t350 * t417;
t208 = qJD(3) * t258 + t345 * t417;
t199 = pkin(5) * t239 + t293;
t164 = mrSges(6,1) * t311 - t489;
t163 = -mrSges(6,2) * t311 + t490;
t159 = -t235 * t342 - t236 * t347;
t158 = -t235 * t347 + t236 * t342;
t157 = t498 + t494;
t153 = -mrSges(5,2) * t275 + mrSges(5,3) * t185;
t152 = mrSges(5,1) * t275 - mrSges(5,3) * t184;
t104 = qJD(4) * t210 + t209 * t349 + t344 * t418;
t103 = qJD(4) * t374 - t209 * t344 + t349 * t418;
t93 = -pkin(5) * t134 + t229;
t89 = mrSges(7,1) * t299 - mrSges(7,3) * t113;
t88 = -mrSges(7,2) * t299 + mrSges(7,3) * t399;
t86 = t184 * Ifges(5,4) + t185 * Ifges(5,2) + t275 * Ifges(5,6);
t81 = -qJD(6) * t196 + t201 * t342 - t202 * t347;
t80 = qJD(6) * t195 - t201 * t347 - t202 * t342;
t71 = -mrSges(6,2) * t269 + mrSges(6,3) * t79;
t70 = mrSges(6,1) * t269 - mrSges(6,3) * t78;
t57 = -qJD(6) * t166 - t133 * t342 + t134 * t347;
t56 = qJD(6) * t165 + t133 * t347 + t134 * t342;
t43 = -pkin(5) * t79 + t82;
t41 = -qJD(5) * t124 + t103 * t348 - t104 * t343;
t40 = qJD(5) * t123 + t103 * t343 + t104 * t348;
t36 = -mrSges(6,1) * t79 + mrSges(6,2) * t78;
t23 = -mrSges(7,2) * t250 + mrSges(7,3) * t27;
t22 = mrSges(7,1) * t250 - mrSges(7,3) * t26;
t19 = t347 * t48 - t479;
t18 = -t342 * t48 - t477;
t15 = -qJD(6) * t68 - t342 * t40 + t347 * t41;
t14 = qJD(6) * t67 + t342 * t41 + t347 * t40;
t10 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t1 = [m(2) * qJDD(1) + t103 * t226 + t104 * t225 + t123 * t70 + t124 * t71 + t14 * t88 + t15 * t89 + t210 * t152 - t374 * t153 + t40 * t163 + t41 * t164 + t209 * t301 + t67 * t22 + t68 * t23 + t258 * t245 + (t36 + t10 + t565) * t257 + (t564 + t592) * t208 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t353 - t213) * t351 + (-mrSges(3,1) * t353 - mrSges(3,2) * qJDD(2) - qJD(2) * t285) * t346) * t340 + (-m(2) - m(3) + t555) * g(3) + m(4) * (t143 * t258 - t144 * t257 - t208 * t227 + t209 * t228 + (-t231 * t351 + t295 * t444) * t340) + m(3) * (qJDD(1) * t341 ^ 2 + (t243 * t351 + t244 * t346) * t340) + m(7) * (t101 * t208 + t14 * t17 + t15 * t16 + t2 * t68 + t257 * t43 + t3 * t67) + m(6) * (t12 * t124 + t123 * t13 + t169 * t208 + t257 * t82 + t40 * t62 + t41 * t61) + m(5) * (t103 * t146 + t104 * t147 + t132 * t257 + t208 * t216 + t210 * t55 - t374 * t54); (-t16 * t56 + t165 * t2 - t166 * t3 + t17 * t57) * mrSges(7,3) + (Ifges(7,5) * t166 + Ifges(7,6) * t165) * t511 + (Ifges(7,5) * t56 + Ifges(7,6) * t57) * t505 + (-t146 * t364 - t147 * t363 - t460 * t55 - t462 * t54) * mrSges(5,3) + (Ifges(6,1) * t133 + Ifges(6,4) * t134) * t512 + (t308 - t244) * mrSges(3,2) + (Ifges(6,5) * t133 + Ifges(6,6) * t134) * t503 + (Ifges(6,4) * t133 + Ifges(6,2) * t134) * t514 + (-t227 * t440 + t559) * mrSges(4,3) + (t307 + t243) * mrSges(3,1) + (-m(6) * t169 - m(7) * t101 + t548 - t592) * t345 * t419 + t133 * t522 + t134 * t524 + t460 * t526 + t56 * t529 + t57 * t531 - t240 * t533 + (-t228 * mrSges(4,3) + t556 / 0.2e1 + t146 * mrSges(5,1) - t147 * mrSges(5,2) + Ifges(6,5) * t512 + Ifges(7,5) * t518 + Ifges(6,6) * t514 + Ifges(7,6) * t520 + Ifges(6,3) * t503 + Ifges(7,3) * t505 - t543) * t442 + (t251 * t588 + t252 * t544) * g(2) + t285 * t420 - t86 * t462 / 0.2e1 + t487 * t585 + t388 * t586 + (t253 * t588 + t254 * t544) * g(1) + (-t12 * t239 + t13 * t240 - t133 * t61 + t134 * t62) * mrSges(6,3) + (-Ifges(6,1) * t240 - Ifges(6,4) * t239) * t528 + (-Ifges(6,5) * t240 - Ifges(6,6) * t239) * t510 + (-Ifges(6,4) * t240 - Ifges(6,2) * t239) * t527 + t82 * (mrSges(6,1) * t239 - mrSges(6,2) * t240) + (t555 * (pkin(2) * t465 + pkin(8) * t467) + (t594 * t346 + t351 * t604) * t340) * g(3) - (t349 * t173 + t344 * t174) * t438 / 0.2e1 - (t429 + t428 + t422) * t350 / 0.2e1 + t278 * (qJD(3) * t360 - t386 * t438) / 0.2e1 + t321 * (qJD(3) * t359 - t383 * t438) / 0.2e1 - t301 * t426 - t239 * t534 + t166 * t539 + t165 * t540 + (Ifges(7,1) * t166 + Ifges(7,4) * t165) * t536 + (Ifges(7,1) * t56 + Ifges(7,4) * t57) * t518 - pkin(2) * t213 + t199 * t10 + t583 * t88 + t584 * t89 + (t101 * t93 + t16 * t584 + t17 * t583 + t199 * t43 + t2 * t51 + t3 * t50) * m(7) + (-t301 * t419 + Ifges(4,4) * t585 + Ifges(4,2) * t586 - Ifges(6,6) * t527 - Ifges(6,5) * t528 + (-Ifges(4,2) * t345 + t487) * t400 - Ifges(7,6) * t535 - Ifges(7,5) * t536 + pkin(8) * t245 - Ifges(6,3) * t510 - Ifges(7,3) * t511 - Ifges(5,6) * t516 - Ifges(5,5) * t517 - Ifges(5,3) * t509 + t546 + t547 + t549) * t350 + (Ifges(4,1) * t292 + Ifges(4,4) * t586 + t384 * t509 + t387 * t516 + t390 * t517) * t345 + t43 * (-mrSges(7,1) * t165 + mrSges(7,2) * t166) + t169 * (-mrSges(6,1) * t134 + mrSges(6,2) * t133) + Ifges(3,3) * qJDD(2) + t574 * t163 + t575 * t164 + (t12 * t121 + t120 * t13 + t169 * t229 + t293 * t82 + t574 * t62 + t575 * t61) * m(6) + t120 * t70 + t121 * t71 + t565 * t336 + t566 * t225 + t567 * t226 + (t233 * t55 + t234 * t54 + (t216 * t440 + t474) * pkin(8) + t566 * t147 + t567 * t146) * m(5) + t101 * (-mrSges(7,1) * t57 + mrSges(7,2) * t56) + t560 * t440 / 0.2e1 + t93 * t63 + t564 * t332 + (-pkin(2) * t231 + ((-t227 * t350 - t228 * t345) * qJD(3) + t559) * pkin(8) - (t295 * t346 + (-t227 * t345 + t228 * t350) * t351) * t448) * m(4) + t51 * t23 + t50 * t22 + (Ifges(7,4) * t56 + Ifges(7,2) * t57) * t520 + (Ifges(7,4) * t166 + Ifges(7,2) * t165) * t535 - t231 * t395 + qJD(3) ^ 2 * t385 / 0.2e1 + (qJD(3) * t361 - t389 * t438) * t507 + t393 * t474 + t413 * t440 + qJD(3) * t372 + t229 * t118 + t371 * t400 + t293 * t36 + t233 * t152 + t234 * t153 + qJDD(3) * (Ifges(4,5) * t345 + Ifges(4,6) * t350) + t216 * (mrSges(5,1) * t363 + mrSges(5,2) * t364); (-t482 + t174 / 0.2e1) * t437 + (t423 - t301) * t227 + ((-t159 + t80) * mrSges(7,2) + (t158 - t81) * mrSges(7,1)) * t101 + (t597 * mrSges(6,1) - mrSges(6,2) * t605) * t169 - (Ifges(6,4) * t512 + Ifges(6,2) * t514 + Ifges(6,6) * t503 + t500 + t524) * t202 + (-t147 * (-mrSges(5,2) * t345 - mrSges(5,3) * t461) - t372 - t146 * (mrSges(5,1) * t345 - mrSges(5,3) * t458)) * qJD(2) + (Ifges(7,4) * t518 + Ifges(7,2) * t520 + Ifges(7,6) * t505 + t493 + t531) * t81 + (Ifges(7,1) * t518 + Ifges(7,4) * t520 + Ifges(7,5) * t505 - t499 + t529) * t80 - t353 * t371 / 0.2e1 - t439 * t481 - t236 * t523 - t235 * t525 + t344 * t526 + t159 * t530 + t158 * t532 + t281 * t533 + (-Ifges(6,1) * t236 - Ifges(6,4) * t235) * t513 + (-Ifges(6,4) * t236 - Ifges(6,2) * t235) * t515 + (-Ifges(6,5) * t236 - Ifges(6,6) * t235) * t504 + (-t12 * t378 - t13 * t281 + t235 * t62 - t236 * t61) * mrSges(6,3) + (Ifges(7,4) * t159 + Ifges(7,2) * t158) * t521 + (Ifges(7,5) * t159 + Ifges(7,6) * t158) * t506 - t378 * t534 + (Ifges(6,4) * t281 - Ifges(6,2) * t378) * t527 + (Ifges(6,1) * t281 - Ifges(6,4) * t378) * t528 + (Ifges(6,5) * t281 - Ifges(6,6) * t378) * t510 + t82 * (mrSges(6,1) * t378 + mrSges(6,2) * t281) + t173 * t416 / 0.2e1 + (-t158 * t17 + t159 * t16 + t195 * t2 - t196 * t3) * mrSges(7,3) + (t278 * t387 + t279 * t390 + t321 * t384) * qJD(4) / 0.2e1 - (t278 * t360 + t279 * t361 + t321 * t359) * qJD(2) / 0.2e1 - t373 * t443 - t385 * t432 / 0.2e1 + (Ifges(7,1) * t159 + Ifges(7,4) * t158) * t519 + (-pkin(3) * t132 - t146 * t167 - t147 * t168) * m(5) + t219 * t70 + t220 * t71 - t168 * t225 - t167 * t226 + (Ifges(7,4) * t196 + Ifges(7,2) * t195) * t535 + (Ifges(7,1) * t196 + Ifges(7,4) * t195) * t536 + t196 * t539 + t195 * t540 + t43 * (-mrSges(7,1) * t195 + mrSges(7,2) * t196) + t579 * t89 + t580 * t88 + (t101 * t569 + t16 * t579 + t17 * t580 + t2 * t91 + t237 * t43 + t3 * t90) * m(7) + t576 * t164 + t577 * t163 + (t12 * t220 + t13 * t219 + t169 * t568 - t330 * t82 + t576 * t61 + t577 * t62) * m(6) + Ifges(4,3) * qJDD(3) - t143 * mrSges(4,2) + t144 * mrSges(4,1) - pkin(3) * t100 + t568 * t118 + t569 * t63 + t90 * t22 - (-Ifges(4,2) * t445 + t331 + t560) * t443 / 0.2e1 + t91 * t23 + (m(5) * ((-t146 * t349 - t147 * t344) * qJD(4) + t558) - t226 * t437 - t225 * t439 + t349 * t153 - t344 * t152) * pkin(9) + t558 * mrSges(5,3) - t556 * t445 / 0.2e1 + (t257 * t552 + t258 * t551) * g(3) + (t551 * t205 - t552 * (-t252 * t345 - t350 * t403)) * g(2) + (t551 * t207 - t552 * (-t254 * t345 + t339 * t466)) * g(1) + (t424 + t548) * t228 + (Ifges(6,5) * t513 + Ifges(7,5) * t519 + Ifges(6,6) * t515 + Ifges(7,6) * t521 + Ifges(6,3) * t504 + Ifges(7,3) * t506 + t543) * t445 + t132 * t394 + (t373 + t413) * qJD(4) + t383 * t509 + (Ifges(7,5) * t196 + Ifges(7,6) * t195) * t511 + t386 * t516 + t389 * t517 + Ifges(4,6) * t291 + Ifges(4,5) * t292 - t330 * t36 + t237 * t10 + t349 * t86 / 0.2e1 - (Ifges(6,1) * t512 + Ifges(6,4) * t514 + Ifges(6,5) * t503 - t501 + t522) * t201; t596 + (t163 * t435 - t164 * t436 + t343 * t71) * pkin(4) - t279 * (Ifges(5,1) * t278 - t480) / 0.2e1 - (t169 * mrSges(6,1) + Ifges(6,4) * t513 + Ifges(6,2) * t515 + Ifges(6,6) * t504 - t500 + t525) * t194 + (-m(7) * (-t258 * t296 - t297 * t465) - mrSges(5,2) * t374 + t571 * t210 + t563) * g(3) - (-Ifges(5,2) * t279 + t174 + t273) * t278 / 0.2e1 + (-t169 * mrSges(6,2) + Ifges(6,1) * t513 + Ifges(6,4) * t515 + Ifges(6,5) * t504 + t501 + t523) * t398 - t146 * t225 + t147 * t226 + (t12 * t343 + t13 * t348 + (-t343 * t61 + t348 * t62) * qJD(5)) * t538 - t64 * t164 - t65 * t163 - t157 * t63 + t573 * t88 + (-t101 * t157 + t572 * t16 + t573 * t17 + t2 * t248 + t247 * t3) * m(7) + (-m(7) * (-t205 * t296 + t251 * t297) - (-t205 * t349 - t251 * t344) * mrSges(5,2) + t571 * (-t205 * t344 + t251 * t349) + t561) * g(2) + (-m(7) * (-t207 * t296 + t253 * t297) - (-t207 * t349 - t253 * t344) * mrSges(5,2) + t571 * (-t207 * t344 + t253 * t349) + t562) * g(1) + t572 * t89 - t547 - t118 * t498 - m(6) * (t169 * t498 + t61 * t64 + t62 * t65) + t173 * t507 + t279 * t481 + t278 * t482 + t70 * t496 - t216 * (mrSges(5,1) * t279 + mrSges(5,2) * t278) - t321 * (Ifges(5,5) * t278 - Ifges(5,6) * t279) / 0.2e1 + t247 * t22 + t248 * t23 + t422; -t63 * t494 - m(7) * (t101 * t494 + t16 * t18 + t17 * t19) + (t2 * t342 + t3 * t347 + (-t16 * t342 + t17 * t347) * qJD(6)) * t537 - t169 * (mrSges(6,1) * t194 + mrSges(6,2) * t398) - t19 * t88 - t18 * t89 + t97 * t512 + (Ifges(6,1) * t398 - t484) * t513 + (Ifges(6,5) * t398 - Ifges(6,6) * t194) * t504 + (t164 + t489) * t62 + (-t163 + t490) * t61 + (-Ifges(6,2) * t194 + t188 + t98) * t515 + (-t375 * t537 + t563) * g(3) + (-t382 * t537 + t561) * g(2) + (-t380 * t537 + t562) * g(1) + (t22 * t347 + t23 * t342 + t433 * t88 - t434 * t89) * pkin(5) + t596; -t101 * (mrSges(7,1) * t113 + mrSges(7,2) * t399) + (Ifges(7,1) * t399 - t483) * t519 + t59 * t518 + (Ifges(7,5) * t399 - Ifges(7,6) * t113) * t506 - t16 * t88 + t17 * t89 - g(1) * t455 - g(2) * t456 - g(3) * t452 + (t113 * t17 + t16 * t399) * mrSges(7,3) + t376 + (-Ifges(7,2) * t113 + t102 + t60) * t521;];
tau  = t1;
