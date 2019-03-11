% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:37:04
% EndTime: 2019-03-09 03:37:14
% DurationCPUTime: 7.62s
% Computational Cost: add. (23977->519), mult. (45664->702), div. (0->0), fcn. (50927->10), ass. (0->299)
t483 = sin(pkin(11));
t484 = cos(pkin(11));
t525 = sin(qJ(3));
t527 = cos(qJ(3));
t311 = t483 * t525 - t484 * t527;
t341 = sin(qJ(6));
t342 = sin(qJ(5));
t343 = cos(qJ(5));
t526 = cos(qJ(6));
t371 = t341 * t343 + t342 * t526;
t237 = t371 * t311;
t573 = -t341 * t342 + t526 * t343;
t239 = t573 * t311;
t101 = t237 * t573 - t239 * t371;
t604 = m(7) * t101 * qJD(4);
t331 = sin(pkin(10)) * pkin(1) + pkin(7);
t414 = t525 * t331;
t303 = -qJ(4) * t525 - t414;
t418 = t527 * t331;
t304 = qJ(4) * t527 + t418;
t251 = -t484 * t303 + t304 * t483;
t313 = -t483 * t527 - t484 * t525;
t433 = t525 * pkin(3);
t261 = -t313 * pkin(4) + t311 * pkin(8) + t433;
t146 = t251 * t342 + t343 * t261;
t147 = -t251 * t343 + t342 * t261;
t171 = mrSges(7,2) * t313 + t237 * mrSges(7,3);
t173 = -mrSges(7,1) * t313 + mrSges(7,3) * t239;
t434 = pkin(5) * t526;
t394 = t434 / 0.2e1;
t464 = t311 * t343;
t409 = -t464 / 0.2e1;
t465 = t311 * t342;
t410 = t465 / 0.2e1;
t562 = m(7) * pkin(5);
t438 = t562 / 0.2e1;
t524 = pkin(5) * t341;
t538 = -t313 / 0.2e1;
t554 = -t239 / 0.2e1;
t556 = t237 / 0.2e1;
t374 = Ifges(7,5) * t554 + Ifges(7,6) * t556;
t602 = t313 / 0.2e1;
t100 = -t313 * pkin(5) + pkin(9) * t464 + t146;
t122 = pkin(9) * t465 + t147;
t66 = t100 * t526 - t341 * t122;
t67 = t341 * t100 + t122 * t526;
t570 = Ifges(7,3) * t602 + t67 * mrSges(7,2) / 0.2e1 - t66 * mrSges(7,1) / 0.2e1 - t374;
t603 = Ifges(6,3) * t538 + t146 * mrSges(6,1) / 0.2e1 - t147 * mrSges(6,2) / 0.2e1 + (t341 * t67 + t526 * t66) * t438 + Ifges(6,5) * t409 + Ifges(6,6) * t410 + t171 * t524 / 0.2e1 + t173 * t394 - t570;
t537 = t573 / 0.2e1;
t419 = -cos(pkin(10)) * pkin(1) - pkin(2);
t319 = -pkin(3) * t527 + t419;
t235 = t311 * pkin(4) + t313 * pkin(8) + t319;
t575 = t483 * t303 + t484 * t304;
t137 = t342 * t235 + t343 * t575;
t463 = t313 * t342;
t104 = pkin(9) * t463 + t137;
t459 = t341 * t104;
t136 = t343 * t235 - t342 * t575;
t462 = t313 * t343;
t103 = pkin(9) * t462 + t136;
t96 = t311 * pkin(5) + t103;
t56 = t526 * t96 - t459;
t71 = t103 * t526 - t459;
t588 = t56 - t71;
t563 = m(7) / 0.2e1;
t601 = t101 * t563;
t413 = t483 * pkin(3);
t330 = t413 + pkin(8);
t515 = pkin(9) + t330;
t305 = t515 * t342;
t306 = t515 * t343;
t260 = -t341 * t305 + t306 * t526;
t396 = -t526 * t305 - t306 * t341;
t441 = Ifges(7,5) * t573 - Ifges(7,6) * t371;
t44 = -t260 * mrSges(7,1) - t396 * mrSges(7,2) + t441;
t600 = t44 * qJD(6);
t236 = t573 * t313;
t238 = t371 * t313;
t138 = t236 * mrSges(7,1) - t238 * mrSges(7,2);
t497 = t238 * mrSges(7,3);
t172 = -mrSges(7,2) * t311 + t497;
t499 = t236 * mrSges(7,3);
t174 = mrSges(7,1) * t311 + t499;
t178 = -pkin(5) * t463 + t251;
t267 = mrSges(7,1) * t371 + mrSges(7,2) * t573;
t415 = t484 * pkin(3);
t332 = -t415 - pkin(4);
t522 = t343 * pkin(5);
t318 = t332 - t522;
t534 = -t318 / 0.2e1;
t539 = t311 / 0.4e1;
t551 = -t260 / 0.2e1;
t552 = t396 / 0.2e1;
t599 = t174 * t551 + t172 * t552 + t178 * t267 / 0.2e1 + t441 * t539 + t138 * t534;
t550 = t260 / 0.2e1;
t417 = t526 * t104;
t57 = t341 * t96 + t417;
t70 = -t341 * t103 - t417;
t587 = t57 + t70;
t338 = t342 ^ 2;
t339 = t343 ^ 2;
t440 = t338 + t339;
t596 = t330 * t440;
t594 = qJD(6) * t138;
t593 = t267 * qJD(6);
t592 = -t371 * t174 / 0.2e1 + t172 * t537;
t140 = -mrSges(7,1) * t238 - mrSges(7,2) * t236;
t321 = t342 * mrSges(6,1) + t343 * mrSges(6,2);
t432 = mrSges(6,3) * t463;
t264 = -mrSges(6,2) * t311 + t432;
t266 = t311 * mrSges(6,1) + mrSges(6,3) * t462;
t529 = -t343 / 0.2e1;
t531 = -t342 / 0.2e1;
t372 = t264 * t531 + t266 * t529;
t335 = Ifges(6,5) * t343;
t501 = Ifges(6,6) * t342;
t398 = t335 - t501;
t503 = Ifges(7,4) * t371;
t270 = Ifges(7,2) * t573 + t503;
t271 = Ifges(7,1) * t573 - t503;
t401 = t270 / 0.4e1 - t271 / 0.4e1;
t310 = Ifges(7,4) * t573;
t269 = -Ifges(7,2) * t371 + t310;
t272 = Ifges(7,1) * t371 + t310;
t402 = t269 / 0.4e1 + t272 / 0.4e1;
t505 = Ifges(6,4) * t342;
t325 = t343 * Ifges(6,1) - t505;
t201 = Ifges(6,5) * t311 - t313 * t325;
t453 = t343 * t201;
t336 = Ifges(6,4) * t343;
t383 = t342 * Ifges(6,2) - t336;
t199 = Ifges(6,6) * t311 + t313 * t383;
t456 = t342 * t199;
t523 = pkin(5) * t342;
t590 = t251 * t321 / 0.2e1 + t398 * t539 - t456 / 0.4e1 + t453 / 0.4e1 + t140 * t523 / 0.2e1 + t402 * t238 + t401 * t236 + t372 * t330 + t599;
t489 = t371 * mrSges(7,3);
t572 = -mrSges(6,1) * t343 + t342 * mrSges(6,2);
t585 = mrSges(5,1) - t572;
t584 = Ifges(5,4) - t335 / 0.2e1 + t501 / 0.2e1;
t369 = t572 * t313;
t231 = Ifges(7,4) * t238;
t119 = -Ifges(7,1) * t236 + t311 * Ifges(7,5) + t231;
t141 = Ifges(7,2) * t236 + t231;
t579 = t119 + t141;
t558 = -t236 / 0.2e1;
t392 = t489 * t558;
t467 = t311 * t267;
t578 = t467 / 0.2e1 - t392;
t506 = Ifges(6,1) * t342;
t577 = t336 + t506;
t504 = Ifges(7,4) * t236;
t117 = Ifges(7,2) * t238 + t311 * Ifges(7,6) - t504;
t142 = Ifges(7,1) * t238 + t504;
t388 = (t117 / 0.4e1 - t142 / 0.4e1) * t371;
t576 = (t141 / 0.4e1 + t119 / 0.4e1) * t573 - t388;
t553 = -t396 / 0.2e1;
t574 = t236 * t550 + t238 * t553;
t376 = -t146 * t342 + t147 * t343;
t490 = t573 * mrSges(7,3);
t423 = -t490 / 0.2e1;
t571 = t238 * t423 + t592;
t496 = t239 * mrSges(7,2);
t498 = t237 * mrSges(7,1);
t448 = t498 / 0.2e1 + t496 / 0.2e1;
t567 = t236 ^ 2;
t566 = t238 ^ 2;
t565 = m(5) / 0.2e1;
t564 = m(6) / 0.2e1;
t557 = t236 / 0.2e1;
t555 = t238 / 0.2e1;
t268 = -mrSges(7,1) * t573 + mrSges(7,2) * t371;
t548 = -t268 / 0.2e1;
t546 = t270 / 0.2e1;
t543 = t272 / 0.2e1;
t541 = -t311 / 0.2e1;
t540 = t311 / 0.2e1;
t535 = t371 / 0.2e1;
t533 = -t325 / 0.4e1;
t532 = t332 / 0.2e1;
t530 = t342 / 0.2e1;
t528 = t343 / 0.2e1;
t521 = t56 * mrSges(7,2);
t520 = t57 * mrSges(7,1);
t517 = t70 * mrSges(7,1);
t516 = t71 * mrSges(7,2);
t508 = mrSges(5,3) * t313;
t507 = mrSges(6,3) * t311;
t502 = Ifges(6,2) * t343;
t491 = t311 * mrSges(5,3);
t277 = t311 * t313;
t363 = m(7) * (t236 * t239 + t237 * t238 - t277);
t364 = m(6) * (t440 - 0.1e1) * t277;
t41 = t363 / 0.2e1 + t364 / 0.2e1;
t482 = qJD(1) * t41;
t358 = t138 * t541 + t172 * t555 + t174 * t557;
t400 = t339 / 0.2e1 + t338 / 0.2e1;
t389 = mrSges(6,3) * t400;
t408 = t464 / 0.2e1;
t10 = -m(7) * (t236 * t588 + t238 * t587) / 0.2e1 + (t566 / 0.2e1 + t567 / 0.2e1) * mrSges(7,3) + (t313 * t389 + t408 * t562 - t540 * t572 + t372) * t313 - t358;
t481 = t10 * qJD(1);
t421 = -t489 / 0.2e1;
t180 = t236 * t421;
t351 = (t237 * t526 - t239 * t341) * t438 + mrSges(6,1) * t410 + mrSges(6,2) * t408 + t448;
t349 = t180 + t351;
t451 = t343 * t264;
t454 = t342 * t266;
t373 = t454 / 0.2e1 - t451 / 0.2e1;
t353 = (-t371 * t588 + t573 * t587) * t563 - t373;
t422 = t490 / 0.2e1;
t11 = t238 * t422 + t349 - t353 - t592;
t480 = t11 * qJD(1);
t479 = t137 * t343;
t15 = t358 - (t566 + t567) * mrSges(7,3) / 0.2e1;
t476 = t15 * qJD(1);
t475 = t178 * t342;
t391 = -t392 + t571;
t20 = t391 - t448;
t474 = t20 * qJD(1);
t473 = t237 * t174;
t472 = t237 * t371;
t471 = t239 * t172;
t470 = t239 * t573;
t468 = t251 * t313;
t466 = t311 * t321;
t458 = t341 * t174;
t265 = -t313 * mrSges(6,1) + mrSges(6,3) * t464;
t455 = t342 * t265;
t263 = mrSges(6,2) * t313 + mrSges(6,3) * t465;
t452 = t343 * t263;
t446 = Ifges(7,5) * t238 + Ifges(7,6) * t236;
t439 = mrSges(7,3) * t524;
t437 = qJD(3) * t601;
t436 = pkin(5) * t462;
t435 = pkin(5) * t465;
t429 = t70 / 0.2e1 + t57 / 0.2e1;
t428 = t71 / 0.2e1 - t56 / 0.2e1;
t420 = t180 + t578;
t322 = t502 + t505;
t405 = t322 * t531;
t403 = t548 - t572 / 0.2e1;
t397 = t236 * t439;
t395 = mrSges(7,3) * t434;
t386 = t238 * t395;
t382 = Ifges(6,5) * t342 + Ifges(6,6) * t343;
t116 = -Ifges(7,4) * t239 + Ifges(7,2) * t237 - Ifges(7,6) * t313;
t118 = -Ifges(7,1) * t239 + Ifges(7,4) * t237 - Ifges(7,5) * t313;
t139 = -t496 - t498;
t177 = t575 - t435;
t198 = -Ifges(6,6) * t313 + t311 * t383;
t200 = -Ifges(6,5) * t313 - t311 * t325;
t256 = t321 * t313;
t301 = t311 * mrSges(5,2);
t365 = mrSges(4,1) * t525 + mrSges(4,2) * t527;
t4 = (-Ifges(4,2) + Ifges(4,1)) * t527 * t525 + t137 * t263 + t147 * t264 + t136 * t265 + t146 * t266 + t57 * t171 + t67 * t172 + t56 * t173 + t66 * t174 + t177 * t140 + t178 * t139 + (-mrSges(5,2) * t433 + t200 * t529 + t198 * t530 + Ifges(7,5) * t557 - Ifges(7,6) * t238 / 0.2e1 - t584 * t313 + (-Ifges(5,2) - Ifges(7,3) + Ifges(5,1) - Ifges(6,3)) * t311) * t313 + (-t453 / 0.2e1 + t456 / 0.2e1 + mrSges(5,1) * t433 + t374 + t584 * t311) * t311 + m(7) * (t177 * t178 + t56 * t66 + t57 * t67) + (-t525 ^ 2 + t527 ^ 2) * Ifges(4,4) + t119 * t554 + t116 * t555 + t117 * t556 + t118 * t558 + m(6) * (t136 * t146 + t137 * t147 + t251 * t575) + t419 * t365 - t575 * t256 - t251 * t466 + (m(5) * t433 - t313 * mrSges(5,1) - t301) * t319;
t362 = -t178 * t313 + t237 * t56 - t239 * t57;
t377 = t136 * t342 - t479;
t8 = t473 / 0.2e1 + t173 * t555 - t471 / 0.2e1 + t171 * t558 + (-t452 / 0.2e1 + t455 / 0.2e1 + t256 / 0.2e1 - t140 / 0.2e1) * t313 + (-t466 / 0.2e1 + t139 / 0.2e1 + t373) * t311 + (t177 * t311 - t236 * t67 + t238 * t66 + t362) * t563 + ((-t251 - t376) * t313 + (t575 + t377) * t311) * t564;
t381 = t4 * qJD(1) + t8 * qJD(2);
t368 = -t178 * t138 + t446 * t540 + t57 * t499;
t5 = -t70 * t174 - m(7) * (t56 * t70 + t57 * t71) - t71 * t172 + t140 * t436 - t251 * t369 + t137 * t266 + (-t117 / 0.2e1 + t142 / 0.2e1) * t236 + (t56 * mrSges(7,3) - t141 / 0.2e1 - t119 / 0.2e1) * t238 + (t201 * t531 + t199 * t529 + m(7) * t178 * t522 + t382 * t541 - mrSges(6,3) * t479 + (t528 * t577 + t405) * t313) * t313 - t368 + (-t264 + t432) * t136;
t380 = -t5 * qJD(1) - t10 * qJD(2);
t9 = -t174 * t57 + t117 * t557 + t142 * t558 + (t172 - t497) * t56 + t368 + t579 * t555;
t379 = t9 * qJD(1) + t15 * qJD(2);
t16 = -t471 + t473 + (-t140 + t256 + t508) * t313 + (-t451 + t454 + t491) * t311 + m(7) * t362 + m(6) * (t311 * t377 - t468) + m(5) * (-t311 * t575 - t468);
t378 = -t16 * qJD(1) - t41 * qJD(2);
t347 = (t343 * t146 + t342 * t147) * t564 + (t371 * t67 + t573 * t66) * t563 + t173 * t537 + t171 * t535 + t263 * t530 + t265 * t528 + t433 * t565;
t360 = (-t311 * t596 - t313 * t332) * t564 + (-t311 * t483 + t313 * t484) * pkin(3) * t565 + (t237 * t396 - t239 * t260 - t313 * t318) * t563;
t348 = (-t472 / 0.2e1 - t470 / 0.2e1) * mrSges(7,3) - t400 * t507 + t360;
t17 = t301 + (mrSges(5,1) + t403) * t313 - t347 + t348;
t366 = -t17 * qJD(1) + qJD(2) * t601;
t14 = -t428 * mrSges(7,2) + t429 * mrSges(7,1) + (-t526 * t172 / 0.2e1 + t458 / 0.2e1 + (t341 * t558 + t526 * t555) * mrSges(7,3)) * pkin(5);
t317 = (mrSges(7,1) * t341 + mrSges(7,2) * t526) * pkin(5);
t45 = (t553 + t552) * mrSges(7,2) + (t551 + t550) * mrSges(7,1);
t361 = t14 * qJD(1) - t45 * qJD(3) + t317 * qJD(5);
t40 = t363 + t364;
t359 = -t8 * qJD(1) - t40 * qJD(2) - t604 / 0.2e1;
t357 = -t260 * t588 + t396 * t587;
t350 = -t392 + t351;
t352 = t435 * t563 + t466 / 0.2e1;
t29 = -t467 / 0.2e1 + t180 + t350 - t352;
t3 = (0.2e1 * t577 - t383) * t463 / 0.4e1 + t579 * t573 / 0.4e1 + t587 * t421 + t574 * mrSges(7,3) + t56 * t423 + t71 * t422 + mrSges(6,3) * t596 * t602 + ((-t318 * t462 + t475) * pkin(5) + t357) * t563 + t436 * t548 + t369 * t532 + t462 * t533 - t388 + t322 * t462 / 0.2e1 + t590 - t603;
t36 = t318 * t267 - (-t271 / 0.2e1 + t546) * t371 + (t543 + t269 / 0.2e1) * t573;
t31 = t332 * t321 + t325 * t530 + t405 + (t577 / 0.2e1 - t383 / 0.2e1) * t343 + t36 + (m(7) * t318 + t268) * t523;
t355 = -t3 * qJD(1) + t29 * qJD(2) - t31 * qJD(3);
t33 = t420 - t448;
t344 = (mrSges(7,3) * t553 + t402) * t238 + (mrSges(7,3) * t550 + t401) * t236 + t576 + t599;
t7 = t344 + t570;
t354 = -t7 * qJD(1) - t33 * qJD(2) - t36 * qJD(3);
t314 = t317 * qJD(6);
t32 = t420 + t448;
t30 = t349 + t352 + t578;
t19 = t391 + t448;
t18 = t313 * t403 + t347 + t348;
t13 = -t520 / 0.2e1 + t172 * t394 + t397 / 0.2e1 - pkin(5) * t458 / 0.2e1 - t386 / 0.2e1 - t521 / 0.2e1 + t517 / 0.2e1 - t516 / 0.2e1 + t446;
t12 = t350 + t353 + t571;
t6 = t344 - t570;
t2 = (t330 * t389 + (t577 / 0.4e1 - t383 / 0.4e1 + mrSges(6,2) * t532 + t506 / 0.4e1) * t342 + (t533 + t322 / 0.4e1 + t505 / 0.2e1 - t332 * mrSges(6,1) / 0.2e1 + t502 / 0.4e1 + (m(7) * t534 + t548) * pkin(5)) * t343) * t313 + (pkin(5) * t475 + t357) * t563 + (-t371 * t429 + t428 * t573 + t574) * mrSges(7,3) + t576 + t590 + t603;
t1 = qJD(3) * t8 + qJD(4) * t41 - qJD(5) * t10 + qJD(6) * t15;
t21 = [qJD(3) * t4 + qJD(4) * t16 - qJD(5) * t5 + qJD(6) * t9, t1 (t260 * t171 + t396 * t173 + m(7) * (t177 * t318 + t260 * t67 + t396 * t66) + (Ifges(7,5) * t371 + Ifges(7,6) * t573 + t382) * t538 + t376 * mrSges(6,3) + t322 * t410 + t237 * t546 + t198 * t528 + t200 * t530 + t118 * t535 + t116 * t537 - t239 * t543 + t413 * t508 - t66 * t489 + t67 * t490 + t415 * t491 + Ifges(4,5) * t527 - Ifges(4,6) * t525 - mrSges(4,1) * t418 + t577 * t409 + (-m(5) * t415 + m(6) * t332 - t585) * t575 - (m(5) * t413 - mrSges(5,2)) * t251 - t332 * t466 + t177 * t268 - Ifges(5,5) * t311 + Ifges(5,6) * t313 + t318 * t139 + mrSges(4,2) * t414 + (m(6) * t376 + t452 - t455) * t330) * qJD(3) + t18 * qJD(4) + t2 * qJD(5) + t6 * qJD(6) + t381, t18 * qJD(3) + t12 * qJD(5) + t19 * qJD(6) - t378 + t604, t2 * qJD(3) + t12 * qJD(4) + (Ifges(6,5) * t463 + Ifges(6,6) * t462 + t517 + (t341 * t71 + t526 * t70) * t562 + t397 - t386 - t516 - t136 * mrSges(6,2) - t137 * mrSges(6,1) + t446) * qJD(5) + t13 * qJD(6) + t380, t6 * qJD(3) + t19 * qJD(4) + t13 * qJD(5) + (t446 - t520 - t521) * qJD(6) + t379; t1, t40 * qJD(3), t30 * qJD(5) + t32 * qJD(6) - t359 + (t301 - t365 + (-t268 + t585) * t313 + 0.2e1 * t360 - t440 * t507 + (-t470 - t472) * mrSges(7,3)) * qJD(3), t437 + t482, -t481 + t30 * qJD(3) + (m(7) * (t236 * t434 + t238 * t524) - t369 + t138) * qJD(5) + t594, qJD(3) * t32 + qJD(5) * t138 + t476 + t594; qJD(4) * t17 + qJD(5) * t3 + qJD(6) * t7 - t381, -t29 * qJD(5) + t33 * qJD(6) + t359, qJD(5) * t31 + qJD(6) * t36, -t366 (-t371 * t439 - t573 * t395 + (-t260 * t526 + t341 * t396) * t562 + t398 + t572 * t330 + t44) * qJD(5) + t600 - t355, t44 * qJD(5) - t354 + t600; -qJD(3) * t17 - qJD(5) * t11 + qJD(6) * t20 + t378, t437 - t482, t366, 0, -t480 - t593 + (-t267 - t321 + (-t371 * t434 + t524 * t573) * m(7)) * qJD(5), -qJD(5) * t267 + t474 - t593; -qJD(3) * t3 + qJD(4) * t11 - qJD(6) * t14 - t380, t29 * qJD(3) + t481, t45 * qJD(6) + t355, t480, -t314, -t314 - t361; -qJD(3) * t7 - qJD(4) * t20 + qJD(5) * t14 - t379, -t33 * qJD(3) - t476, -t45 * qJD(5) + t354, -t474, t361, 0;];
Cq  = t21;
