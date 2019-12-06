% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:27
% EndTime: 2019-12-05 17:23:50
% DurationCPUTime: 9.82s
% Computational Cost: add. (14934->717), mult. (40950->1061), div. (0->0), fcn. (44135->12), ass. (0->383)
t625 = m(5) / 0.2e1;
t392 = sin(qJ(3));
t396 = cos(qJ(3));
t393 = sin(qJ(2));
t518 = sin(pkin(5));
t519 = cos(pkin(6));
t447 = t519 * t518;
t436 = t393 * t447;
t562 = cos(qJ(2));
t448 = t562 * t518;
t293 = -t392 * t436 + t396 * t448;
t391 = sin(qJ(4));
t395 = cos(qJ(4));
t389 = sin(pkin(6));
t459 = t393 * t518;
t450 = t389 * t459;
t238 = t293 * t395 + t391 * t450;
t292 = t392 * t448 + t396 * t436;
t390 = sin(qJ(5));
t394 = cos(qJ(5));
t146 = -t238 * t390 + t292 * t394;
t147 = t238 * t394 + t292 * t390;
t354 = -pkin(4) * t395 - pkin(10) * t391 - pkin(3);
t497 = t390 * t395;
t296 = -pkin(9) * t497 + t354 * t394;
t491 = t394 * t395;
t297 = pkin(9) * t491 + t354 * t390;
t528 = t394 * mrSges(6,2);
t533 = t390 * mrSges(6,1);
t357 = t528 + t533;
t338 = t357 * t391;
t526 = t395 * mrSges(5,1);
t356 = mrSges(5,2) * t391 - t526;
t237 = t293 * t391 - t395 * t450;
t510 = t237 * t391;
t496 = t391 * t394;
t347 = -mrSges(6,1) * t395 - mrSges(6,3) * t496;
t583 = -t347 / 0.2e1;
t498 = t390 * t391;
t345 = mrSges(6,2) * t395 - mrSges(6,3) * t498;
t586 = -t345 / 0.2e1;
t607 = -t237 / 0.2e1;
t624 = -m(6) / 0.2e1;
t509 = t238 * t395;
t632 = t509 + t510;
t639 = (t356 / 0.2e1 - mrSges(4,1) / 0.2e1) * t292 + (-pkin(3) * t292 + t632 * pkin(9)) * t625 - (pkin(9) * t510 + t146 * t296 + t147 * t297) * t624 - t146 * t583 - t147 * t586 - t338 * t607 - t293 * mrSges(4,2) / 0.2e1;
t475 = pkin(2) * t519;
t454 = t392 * t475;
t501 = t389 * t396;
t335 = pkin(8) * t501 + t454;
t304 = t519 * pkin(9) + t335;
t305 = (-pkin(3) * t396 - pkin(9) * t392 - pkin(2)) * t389;
t200 = -t391 * t304 + t395 * t305;
t182 = pkin(4) * t501 - t200;
t561 = m(6) * t182;
t502 = t389 * t392;
t331 = t391 * t519 + t395 * t502;
t534 = t331 * mrSges(5,3);
t638 = -mrSges(4,1) + t356;
t530 = t390 * Ifges(6,4);
t365 = Ifges(6,1) * t394 - t530;
t548 = Ifges(6,5) * t395;
t312 = t365 * t391 - t548;
t360 = t394 * Ifges(6,2) + t530;
t340 = t360 * t391;
t637 = t312 - t340;
t381 = Ifges(6,5) * t394;
t544 = Ifges(6,6) * t390;
t636 = t381 - t544;
t383 = Ifges(6,4) * t394;
t635 = -t390 * Ifges(6,2) + t383;
t384 = Ifges(5,4) * t395;
t363 = -Ifges(5,2) * t391 + t384;
t520 = cos(pkin(5));
t634 = t520 * t389 + t562 * t447;
t371 = pkin(4) * t391 - pkin(10) * t395;
t299 = pkin(9) * t498 + t371 * t394;
t300 = -pkin(9) * t496 + t371 * t390;
t633 = -t299 * t390 + t300 * t394;
t631 = -t390 * Ifges(6,1) - t383;
t574 = -t360 / 0.4e1;
t460 = t365 / 0.4e1 + t574;
t630 = m(6) * pkin(10) + mrSges(6,3);
t446 = mrSges(6,1) * t394 - mrSges(6,2) * t390;
t628 = -m(6) * pkin(4) - mrSges(5,1) - t446;
t627 = t389 ^ 2;
t626 = 2 * qJD(3);
t623 = m(6) / 0.2e1;
t622 = -mrSges(5,1) / 0.2e1;
t621 = mrSges(6,1) / 0.2e1;
t620 = mrSges(5,2) / 0.2e1;
t619 = -mrSges(6,2) / 0.2e1;
t333 = -pkin(8) * t502 + t396 * t475;
t334 = (pkin(3) * t392 - pkin(9) * t396) * t389;
t236 = t395 * t333 + t391 * t334;
t209 = pkin(10) * t502 + t236;
t251 = t454 + (pkin(8) + t371) * t501;
t100 = -t209 * t390 + t251 * t394;
t617 = t100 / 0.2e1;
t101 = t209 * t394 + t251 * t390;
t616 = -t101 / 0.2e1;
t267 = -t390 * t331 - t394 * t501;
t266 = Ifges(6,4) * t267;
t268 = t331 * t394 - t390 * t501;
t330 = t391 * t502 - t395 * t519;
t549 = Ifges(6,5) * t330;
t111 = Ifges(6,1) * t268 + t266 + t549;
t615 = -t111 / 0.4e1;
t156 = mrSges(6,1) * t268 + mrSges(6,2) * t267;
t614 = t156 / 0.2e1;
t263 = t634 * t392 + t396 * t459;
t412 = -t389 * t448 + t520 * t519;
t176 = t263 * t391 - t412 * t395;
t613 = t176 / 0.2e1;
t177 = t263 * t395 + t412 * t391;
t612 = t177 / 0.2e1;
t539 = t267 * mrSges(6,3);
t194 = -mrSges(6,2) * t330 + t539;
t611 = t194 / 0.2e1;
t538 = t268 * mrSges(6,3);
t195 = mrSges(6,1) * t330 - t538;
t610 = -t195 / 0.2e1;
t609 = t195 / 0.2e1;
t490 = t395 * t396;
t294 = (-t390 * t490 + t392 * t394) * t389;
t295 = (t390 * t392 + t394 * t490) * t389;
t202 = -mrSges(6,1) * t294 + mrSges(6,2) * t295;
t608 = t202 / 0.2e1;
t495 = t391 * t396;
t483 = t389 * t495;
t253 = -mrSges(6,2) * t483 + t294 * mrSges(6,3);
t606 = t253 / 0.2e1;
t254 = mrSges(6,1) * t483 - t295 * mrSges(6,3);
t605 = t254 / 0.2e1;
t262 = t392 * t459 - t634 * t396;
t604 = t262 / 0.2e1;
t603 = t267 / 0.2e1;
t602 = t267 / 0.4e1;
t601 = -t268 / 0.2e1;
t600 = t268 / 0.2e1;
t599 = t268 / 0.4e1;
t552 = Ifges(5,4) * t391;
t367 = Ifges(5,1) * t395 - t552;
t598 = (Ifges(5,5) * t392 + t367 * t396) * t389 / 0.2e1;
t536 = t330 * mrSges(5,3);
t278 = mrSges(5,2) * t501 - t536;
t597 = -t278 / 0.2e1;
t596 = t294 / 0.2e1;
t595 = t295 / 0.2e1;
t594 = t296 / 0.2e1;
t593 = t297 / 0.2e1;
t592 = t300 / 0.2e1;
t591 = t330 / 0.2e1;
t590 = t330 / 0.4e1;
t337 = t446 * t391;
t589 = -t337 / 0.2e1;
t588 = t338 / 0.2e1;
t339 = t357 * t395;
t587 = t339 / 0.2e1;
t585 = t345 / 0.2e1;
t346 = -mrSges(6,2) * t391 - mrSges(6,3) * t497;
t584 = t346 / 0.2e1;
t582 = t347 / 0.2e1;
t348 = mrSges(6,1) * t391 - mrSges(6,3) * t491;
t581 = t348 / 0.2e1;
t580 = t446 / 0.2e1;
t579 = -t446 / 0.2e1;
t578 = -t357 / 0.2e1;
t577 = t357 / 0.2e1;
t525 = t395 * mrSges(5,2);
t358 = t391 * mrSges(5,1) + t525;
t576 = t358 / 0.2e1;
t359 = Ifges(6,5) * t390 + Ifges(6,6) * t394;
t575 = t359 / 0.2e1;
t573 = t631 / 0.4e1;
t571 = -t381 / 0.4e1;
t570 = -t390 / 0.2e1;
t569 = t390 / 0.2e1;
t568 = t390 / 0.4e1;
t566 = -t394 / 0.2e1;
t565 = t394 / 0.2e1;
t564 = t394 / 0.4e1;
t563 = t395 / 0.2e1;
t560 = pkin(9) * t262;
t559 = pkin(9) * t391;
t558 = pkin(9) * t395;
t557 = pkin(10) * t390;
t556 = pkin(10) * t394;
t555 = Ifges(6,1) * t267;
t554 = Ifges(4,4) * t392;
t553 = Ifges(5,4) * t331;
t551 = Ifges(6,4) * t268;
t382 = Ifges(5,5) * t395;
t550 = Ifges(6,5) * t295;
t546 = Ifges(6,2) * t268;
t545 = Ifges(6,6) * t294;
t543 = Ifges(6,6) * t395;
t542 = Ifges(5,3) * t392;
t541 = Ifges(6,3) * t331;
t540 = Ifges(6,3) * t391;
t537 = t297 * mrSges(6,3);
t535 = t330 * Ifges(6,6);
t532 = t390 * mrSges(6,3);
t527 = t394 * mrSges(6,3);
t303 = -t519 * pkin(3) - t333;
t178 = t330 * pkin(4) - t331 * pkin(10) + t303;
t201 = t395 * t304 + t391 * t305;
t183 = -pkin(10) * t501 + t201;
t67 = t178 * t394 - t183 * t390;
t524 = t67 * t390;
t68 = t178 * t390 + t183 * t394;
t523 = t68 * t394;
t86 = -t177 * t390 + t262 * t394;
t522 = t86 * t390;
t87 = t177 * t394 + t262 * t390;
t521 = t87 * t394;
t124 = t262 * t497 + t263 * t394;
t517 = t124 * t390;
t125 = -t262 * t491 + t263 * t390;
t516 = t125 * t394;
t515 = t146 * t390;
t514 = t147 * t394;
t513 = t176 * t237;
t512 = t176 * t391;
t23 = m(5) * (-t177 * t395 + t263 - t512) * t262 + (t124 * t86 + t125 * t87 - t512 * t262) * m(6);
t511 = t23 * qJD(1);
t506 = t262 * t292;
t24 = m(6) * (t146 * t86 + t147 * t87 + t513) + m(5) * (t177 * t238 + t506 + t513) + m(4) * (t263 * t293 + t412 * t450 + t506);
t508 = t24 * qJD(1);
t25 = m(6) * (t177 - t521 + t522) * t176;
t507 = t25 * qJD(1);
t505 = t262 * t391;
t110 = Ifges(6,2) * t267 + t535 + t551;
t500 = t390 * t110;
t310 = t391 * t635 - t543;
t499 = t390 * t310;
t494 = t394 * t110;
t493 = t394 * t111;
t492 = t394 * t312;
t489 = Ifges(6,5) * t267 - Ifges(6,6) * t268;
t488 = -Ifges(5,5) * t330 - Ifges(5,6) * t331;
t487 = mrSges(4,3) * t502;
t486 = mrSges(4,3) * t501;
t485 = t201 * t623;
t484 = -t557 / 0.2e1;
t482 = -t539 / 0.2e1;
t481 = -t534 / 0.2e1;
t480 = -t532 / 0.2e1;
t479 = t532 / 0.2e1;
t478 = -t527 / 0.2e1;
t477 = t527 / 0.2e1;
t476 = t579 + t622;
t474 = -t502 / 0.2e1;
t473 = t502 / 0.2e1;
t472 = -t501 / 0.2e1;
t471 = t501 / 0.2e1;
t469 = -t497 / 0.2e1;
t468 = t491 / 0.2e1;
t467 = -t396 * t382 / 0.4e1;
t215 = t357 * t330;
t466 = -t215 / 0.2e1 + t597;
t157 = -mrSges(6,1) * t267 + mrSges(6,2) * t268;
t279 = -mrSges(5,1) * t501 - t534;
t465 = t279 / 0.2e1 - t157 / 0.2e1;
t341 = t631 * t391;
t464 = t310 / 0.4e1 - t341 / 0.4e1;
t463 = -t312 / 0.4e1 + t340 / 0.4e1;
t308 = -Ifges(6,3) * t395 + t391 * t636;
t362 = Ifges(5,2) * t395 + t552;
t462 = -t362 / 0.2e1 + t308 / 0.2e1;
t461 = t573 - t635 / 0.4e1;
t457 = t266 - t546;
t453 = t391 * t472;
t452 = t391 * t471;
t451 = t395 * t471;
t449 = t627 * t459;
t445 = -t551 + t555;
t244 = mrSges(5,1) * t330 + mrSges(5,2) * t331;
t317 = (-mrSges(5,2) * t392 - mrSges(5,3) * t495) * t389;
t318 = (mrSges(5,1) * t392 - mrSges(5,3) * t490) * t389;
t332 = (mrSges(4,1) * t392 + mrSges(4,2) * t396) * t389;
t343 = t519 * mrSges(4,1) - t487;
t398 = (-t318 / 0.2e1 + t608) * t176 + (mrSges(4,3) * t474 - t343 / 0.2e1 + t244 / 0.2e1) * t263 + t124 * t609 + t125 * t611 + t317 * t612 + t86 * t605 + t87 * t606 + t412 * t332 / 0.2e1;
t298 = t358 * t501;
t344 = -t519 * mrSges(4,2) + t486;
t408 = mrSges(4,3) * t471 - t344 / 0.2e1 + t298 / 0.2e1 + t395 * t597 + t465 * t391;
t235 = -t333 * t391 + t334 * t395;
t208 = -pkin(4) * t502 - t235;
t410 = t100 * t86 + t101 * t87 + t124 * t67 + t125 * t68 + t176 * t208;
t418 = -t176 * t235 + t177 * t236 + t263 * t303;
t425 = t200 * t391 - t201 * t395 + t335;
t4 = (t425 * t262 + t418) * t625 + t408 * t262 + (-t509 / 0.2e1 - t510 / 0.2e1) * mrSges(5,3) + (-t182 * t505 + t410) * t623 + t398 - t639;
t109 = Ifges(6,5) * t268 + Ifges(6,6) * t267 + Ifges(6,3) * t330;
t171 = Ifges(6,3) * t483 + t545 + t550;
t172 = Ifges(6,4) * t295 + Ifges(6,2) * t294 + Ifges(6,6) * t483;
t173 = Ifges(6,1) * t295 + Ifges(6,4) * t294 + Ifges(6,5) * t483;
t210 = -Ifges(5,2) * t330 - Ifges(5,6) * t501 + t553;
t324 = Ifges(5,4) * t330;
t211 = Ifges(5,1) * t331 - Ifges(5,5) * t501 - t324;
t269 = (Ifges(5,6) * t392 + t363 * t396) * t389;
t376 = Ifges(4,5) * t501;
t5 = (Ifges(4,5) * t519 + 0.2e1 * Ifges(4,4) * t501 + (Ifges(4,1) - Ifges(4,2)) * t502) * t471 + t111 * t595 + t110 * t596 + t331 * t598 + t173 * t600 + t172 * t603 + t171 * t591 + m(6) * (t100 * t67 + t101 * t68 + t182 * t208) + (-t396 * (t542 + t396 * (-Ifges(5,6) * t391 + t382)) / 0.2e1 + t392 * (Ifges(4,1) * t396 - t554) / 0.2e1) * t627 - t330 * t269 / 0.2e1 + t201 * t317 + t200 * t318 + t303 * t298 + t236 * t278 + t235 * t279 + t68 * t253 + t67 * t254 + t208 * t157 + t182 * t202 + t101 * t194 + t100 * t195 + (m(5) * t303 + t244 - t343 - t487) * t335 + (t344 - t486) * t333 + m(5) * (t200 * t235 + t201 * t236) + t210 * t453 + t211 * t451 + t109 * t452 - t389 * pkin(2) * t332 + (Ifges(5,5) * t331 - Ifges(5,6) * t330 - Ifges(5,3) * t501) * t473 + t519 * (-Ifges(4,6) * t502 + t376) / 0.2e1 + (Ifges(4,6) * t519 + (Ifges(4,2) * t396 + t554) * t389) * t474;
t442 = t4 * qJD(1) + t5 * qJD(2);
t233 = -mrSges(6,2) * t331 + t330 * t532;
t234 = mrSges(6,1) * t331 + t330 * t527;
t243 = mrSges(5,1) * t331 - mrSges(5,2) * t330;
t426 = -t536 / 0.2e1 + t466;
t247 = pkin(4) * t331 + pkin(10) * t330;
t94 = -t200 * t390 + t247 * t394;
t95 = t200 * t394 + t247 * t390;
t399 = ((t201 - t523 + t524) * t623 + t194 * t566 + t195 * t569 + t426) * t176 + (t177 * t182 + t86 * t94 + t87 * t95) * t623 + t243 * t604 + t86 * t234 / 0.2e1 + t87 * t233 / 0.2e1;
t409 = (-pkin(4) * t237 + (t514 - t515) * pkin(10)) * t624 + t238 * t620;
t427 = t481 - t465;
t7 = -t476 * t237 + (-t514 / 0.2e1 + t515 / 0.2e1) * mrSges(6,3) + t427 * t177 + t399 + t409;
t168 = -t330 * t636 + t541;
t169 = Ifges(6,6) * t331 - t330 * t635;
t170 = Ifges(6,5) * t331 - t365 * t330;
t245 = -Ifges(5,2) * t331 - t324;
t246 = -Ifges(5,1) * t330 - t553;
t8 = m(6) * (t67 * t94 + t68 * t95) + t94 * t195 + t67 * t234 + t95 * t194 + t68 * t233 + t170 * t600 + t169 * t603 - t182 * t215 + t200 * t278 + t303 * t243 + t488 * t472 + (t246 / 0.2e1 - t210 / 0.2e1 + t109 / 0.2e1) * t331 + (-t245 / 0.2e1 + t168 / 0.2e1 - t211 / 0.2e1 + t500 / 0.2e1 + t200 * mrSges(5,3) - t493 / 0.2e1) * t330 + (t157 - t279 - t534 + t561) * t201;
t441 = t7 * qJD(1) + t8 * qJD(2);
t439 = t390 * t457;
t13 = -t68 * t195 + t67 * t194 + t445 * t600 + t489 * t591 + t110 * t601 + t182 * t156 + (-t267 * t67 - t268 * t68) * mrSges(6,3) + (t111 + t457) * t603;
t407 = (t87 * t601 - t86 * t267 / 0.2e1) * mrSges(6,3) + t156 * t613 + t86 * t611 + t87 * t610;
t432 = t146 * t621 + t147 * t619;
t14 = t407 - t432;
t438 = t14 * qJD(1) + t13 * qJD(2);
t437 = t296 * t390 - t297 * t394;
t435 = pkin(4) * t614 + t182 * t578;
t434 = pkin(4) * t589 + t395 * t571;
t433 = t124 * t621 + t125 * t619;
t431 = -t299 * mrSges(6,1) / 0.2e1 + mrSges(6,2) * t592;
t430 = t528 / 0.2e1 + t533 / 0.2e1;
t429 = pkin(10) * t586 - t464;
t428 = pkin(10) * t583 - t463;
t424 = -t500 / 0.4e1 + t493 / 0.4e1;
t423 = t394 * t445;
t422 = t391 * t359;
t400 = ((t437 + t558) * t623 + t587 + t345 * t566 + t347 * t569) * t176 + (t177 * t559 + t299 * t86 + t300 * t87) * t623 + t177 * t588 + t86 * t581 + t87 * t584;
t413 = m(6) * (pkin(4) * t505 + (t516 - t517) * pkin(10));
t12 = (-t516 / 0.2e1 + t517 / 0.2e1) * mrSges(6,3) + (t576 - t525 / 0.2e1 + t476 * t391) * t262 - t413 / 0.2e1 + t400;
t311 = Ifges(6,6) * t391 + t395 * t635;
t313 = Ifges(6,5) * t391 + t365 * t395;
t397 = (-t362 / 0.4e1 + t367 / 0.4e1 + t308 / 0.4e1) * t331 + (t296 * t94 + t297 * t95 + t299 * t67 + t300 * t68) * t623 - pkin(3) * t243 / 0.2e1 + t182 * t587 + t201 * t588 + t311 * t602 + t313 * t599 + t234 * t594 + t233 * t593 + t299 * t609 + t194 * t592 + t303 * t576 + t67 * t581 + t68 * t584 + t94 * t582 + t95 * t585;
t402 = (-pkin(4) * t208 + (-t100 * t390 + t101 * t394) * pkin(10)) * t624 + pkin(4) * t608 + t208 * t580 + t235 * t622 + t236 * t620 + t294 * t574 + t295 * t573;
t403 = t245 / 0.4e1 + t211 / 0.4e1 - t168 / 0.4e1 + (t561 / 0.2e1 + t427) * pkin(9) + t424;
t309 = t395 * t636 + t540;
t366 = Ifges(5,1) * t391 + t384;
t414 = -t363 / 0.4e1 - t366 / 0.4e1 + t309 / 0.4e1 - t492 / 0.4e1 + t499 / 0.4e1;
t415 = t109 / 0.4e1 - t210 / 0.4e1 + t246 / 0.4e1 - t390 * t169 / 0.4e1 + t170 * t564;
t2 = ((0.3e1 / 0.4e1 * Ifges(5,6) - t359 / 0.4e1) * t501 + (t485 + t426) * pkin(9) + t415) * t391 + (Ifges(5,5) * t472 + t403) * t395 + (-t172 / 0.4e1 - pkin(10) * t253 / 0.2e1 + mrSges(6,3) * t616) * t394 + (-t173 / 0.4e1 + pkin(10) * t605 + mrSges(6,3) * t617) * t390 + (t467 - t542 / 0.2e1) * t389 + t414 * t330 + t397 + t402;
t31 = m(6) * (t296 * t299 + t297 * t300) + t300 * t345 + t297 * t346 + t299 * t347 + t296 * t348 - pkin(3) * t358 + (pkin(9) * t338 + t492 / 0.2e1 - t499 / 0.2e1 - t309 / 0.2e1 + t363 / 0.2e1 + t366 / 0.2e1) * t395 + (t313 * t565 + t311 * t570 + t367 / 0.2e1 + (m(6) * t558 + t339) * pkin(9) + t462) * t391;
t420 = t12 * qJD(1) + t2 * qJD(2) + t31 * qJD(3);
t404 = t182 * t589 - t296 * t194 / 0.2e1 + t195 * t593 + t395 * t489 / 0.4e1 + t67 * t586 + t68 * t582;
t411 = t550 / 0.2e1 + t545 / 0.2e1 + mrSges(6,1) * t617 + mrSges(6,2) * t616;
t10 = (t537 / 0.2e1 + t464) * t268 + (mrSges(6,3) * t594 + t463) * t267 + (Ifges(6,3) * t471 + t359 * t590 - t423 / 0.4e1 + t494 / 0.4e1 - pkin(9) * t156 / 0.2e1 + t439 / 0.4e1 + t111 * t568 + (-t524 / 0.2e1 + t523 / 0.2e1) * mrSges(6,3)) * t391 + t404 + t411;
t406 = (-t521 / 0.2e1 + t522 / 0.2e1) * t391 * mrSges(6,3) + t337 * t613 + t86 * t585 + t87 * t583;
t18 = t406 - t433;
t38 = t296 * t345 - t297 * t347 + (t437 * mrSges(6,3) + pkin(9) * t337 + t310 * t566 + t341 * t565 + t359 * t563 + t637 * t570) * t391;
t419 = t18 * qJD(1) - t10 * qJD(2) + t38 * qJD(3);
t417 = t541 / 0.2e1 + t94 * t621 + t95 * t619;
t151 = -pkin(4) * t357 + (-t631 / 0.2e1 + t635 / 0.2e1) * t394 + (t365 / 0.2e1 - t360 / 0.2e1) * t390;
t16 = t330 * t571 - t460 * t268 + t461 * t267 + (-t549 / 0.2e1 + t546 / 0.4e1 - t266 / 0.4e1 + t615 + (t609 + t538 / 0.2e1) * pkin(10)) * t394 + (0.3e1 / 0.4e1 * t535 + t110 / 0.4e1 - t555 / 0.4e1 + t551 / 0.4e1 + (t611 + t482) * pkin(10)) * t390 + t417 + t435;
t27 = (t578 + t430) * t176;
t385 = t390 ^ 2;
t387 = t394 ^ 2;
t405 = pkin(9) * t577 + t460 * t394 + t461 * t390 + (-t387 / 0.2e1 - t385 / 0.2e1) * pkin(10) * mrSges(6,3);
t35 = (-t548 / 0.2e1 + t428) * t394 + (0.3e1 / 0.4e1 * t543 + t429) * t390 + (-Ifges(6,3) / 0.2e1 + t405) * t391 + t431 + t434;
t416 = t27 * qJD(1) + t16 * qJD(2) - t35 * qJD(3) - t151 * qJD(4);
t388 = t395 ^ 2;
t386 = t391 ^ 2;
t250 = t386 * t560;
t34 = Ifges(6,5) * t468 + Ifges(6,6) * t469 + t540 / 0.2e1 + (t543 / 0.4e1 + t429) * t390 + t428 * t394 + t405 * t391 - t431 + t434;
t28 = (t430 + t577) * t176;
t19 = t406 + t433;
t17 = t556 * t610 + pkin(10) * t267 * t479 + t194 * t484 + t636 * t590 + t457 * t564 + t445 * t568 + (-t381 / 0.2e1 + t544 / 0.2e1) * t330 + t417 + t424 - t435 + (-t631 + t635) * t602 + (pkin(10) * t478 + t460) * t268;
t15 = t407 + t432;
t11 = t262 * t576 + t413 / 0.2e1 + t525 * t604 + t125 * t477 + t124 * t480 + t400 + (t580 + mrSges(5,1) / 0.2e1) * t505;
t9 = Ifges(6,3) * t452 + t411 + t341 * t599 - t268 * t310 / 0.4e1 - t404 + t559 * t614 + t498 * t615 + t296 * t482 + t537 * t601 - t422 * t590 + t637 * t602 - (t494 + t439) * t391 / 0.4e1 + (t67 * t479 + t68 * t478 + t423 / 0.4e1) * t391;
t6 = mrSges(5,1) * t607 + t146 * t480 + t147 * t477 + t157 * t612 + t237 * t579 + t399 - t409 + (t481 - t279 / 0.2e1) * t177;
t3 = (t425 * t625 - t391 * t561 / 0.2e1 + t408) * t262 + t418 * t625 + t410 * t623 + t398 + t632 * mrSges(5,3) / 0.2e1 + t639;
t1 = Ifges(5,5) * t451 + Ifges(5,6) * t453 + (-mrSges(5,3) * t559 / 0.2e1 + t414) * t330 + t389 * t467 + t403 * t395 + t172 * t564 + t173 * t568 + t556 * t606 + t101 * t477 + Ifges(5,3) * t473 + t254 * t484 + t100 * t480 + t397 - t402 + ((t485 + t466) * pkin(9) + t415 + (t359 + Ifges(5,6)) * t501 / 0.4e1) * t391;
t20 = [t24 * qJD(2) + t23 * qJD(3) + t25 * qJD(4), t3 * qJD(3) + t6 * qJD(4) + t15 * qJD(5) + t508 + (-mrSges(3,2) * t448 - mrSges(3,1) * t459 + t237 * t157 + t292 * t244 + t238 * t278 - t237 * t279 + t147 * t194 + t146 * t195 + (-mrSges(4,1) * t396 + mrSges(4,2) * t392) * t449 + t293 * t344 - t292 * t343 + 0.2e1 * (t146 * t67 + t147 * t68 + t182 * t237) * t623 + 0.2e1 * (-t200 * t237 + t201 * t238 + t292 * t303) * t625 + m(4) * (-pkin(2) * t449 - t333 * t292 + t335 * t293)) * qJD(2), t511 + t3 * qJD(2) + t11 * qJD(4) + t19 * qJD(5) + ((t124 * t296 + t125 * t297 - t250) * t623 + (-pkin(3) * t263 - t388 * t560 - t250) * t625) * t626 + (t124 * t347 + t125 * t345 + (-t391 * t338 + mrSges(4,2) + (-t386 - t388) * mrSges(5,3)) * t262 + t638 * t263) * qJD(3), t507 + t6 * qJD(2) + t11 * qJD(3) + (t628 * t177 + (mrSges(5,2) + t630 * (-t385 - t387)) * t176) * qJD(4) + t28 * qJD(5), t15 * qJD(2) + t19 * qJD(3) + t28 * qJD(4) + (-mrSges(6,1) * t87 - mrSges(6,2) * t86) * qJD(5); qJD(3) * t4 + qJD(4) * t7 + qJD(5) * t14 - t508, qJD(3) * t5 + qJD(4) * t8 + qJD(5) * t13, t1 * qJD(4) + t9 * qJD(5) + ((-pkin(3) * t335 + (-t235 * t391 + t236 * t395) * pkin(9)) * t625 + (t100 * t296 + t101 * t297 + t208 * t559) * t623) * t626 + t442 + (-t333 * mrSges(4,2) - pkin(3) * t298 + t100 * t347 + t101 * t345 + t208 * t338 + t297 * t253 + t296 * t254 + t310 * t596 + t312 * t595 + t376 + (t269 / 0.2e1 - t171 / 0.2e1 + pkin(9) * t317 + t236 * mrSges(5,3)) * t395 + (t598 + t173 * t565 + t172 * t570 - t235 * mrSges(5,3) + (-t318 + t202) * pkin(9)) * t391 + ((Ifges(5,5) * t391 / 0.2e1 + Ifges(5,6) * t563 - Ifges(4,6)) * t392 + (t366 * t563 + t391 * t462) * t396) * t389 + t638 * t335) * qJD(3), t1 * qJD(3) + (-t234 * t557 + t233 * t556 + t331 * t575 + t169 * t565 + t170 * t569 + pkin(4) * t215 - t200 * mrSges(5,2) + (t360 * t569 - t566 * t631) * t330 + t488 + t630 * (-t390 * t94 + t394 * t95) + t628 * t201) * qJD(4) + t17 * qJD(5) + t441, t9 * qJD(3) + t17 * qJD(4) + (-mrSges(6,1) * t68 - mrSges(6,2) * t67 + t489) * qJD(5) + t438; -qJD(2) * t4 + qJD(4) * t12 + qJD(5) * t18 - t511, qJD(4) * t2 - qJD(5) * t10 - t442, qJD(4) * t31 + qJD(5) * t38, t34 * qJD(5) + t420 + (-t446 * t558 - t631 * t468 + t360 * t469 - pkin(9) * t526 + m(6) * (-pkin(4) * t558 + t633 * pkin(10)) + t346 * t556 - t348 * t557 + t313 * t569 + t311 * t565 - pkin(4) * t339 + t382 + (pkin(9) * mrSges(5,2) - Ifges(5,6) + t575) * t391 + t633 * mrSges(6,3)) * qJD(4), t34 * qJD(4) + (-mrSges(6,1) * t297 - mrSges(6,2) * t296 - t422) * qJD(5) + t419; -qJD(2) * t7 - qJD(3) * t12 - qJD(5) * t27 - t507, -qJD(3) * t2 - qJD(5) * t16 - t441, qJD(5) * t35 - t420, t151 * qJD(5), (-pkin(10) * t446 + t636) * qJD(5) - t416; -t14 * qJD(2) - t18 * qJD(3) + t27 * qJD(4), qJD(3) * t10 + qJD(4) * t16 - t438, -qJD(4) * t35 - t419, t416, 0;];
Cq = t20;
