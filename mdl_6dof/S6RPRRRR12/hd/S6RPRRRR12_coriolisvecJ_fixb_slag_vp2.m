% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:40:47
% EndTime: 2018-11-23 16:41:21
% DurationCPUTime: 34.00s
% Computational Cost: add. (63891->918), mult. (220232->1335), div. (0->0), fcn. (192757->16), ass. (0->425)
t362 = sin(pkin(14));
t365 = sin(pkin(6));
t366 = cos(pkin(14));
t369 = cos(pkin(6));
t377 = cos(qJ(3));
t368 = cos(pkin(7));
t373 = sin(qJ(3));
t476 = t368 * t373;
t364 = sin(pkin(7));
t481 = t364 * t373;
t301 = t365 * (t362 * t377 + t366 * t476) + t369 * t481;
t294 = t301 * qJD(1);
t372 = sin(qJ(4));
t376 = cos(qJ(4));
t479 = t365 * t366;
t336 = -t364 * t479 + t368 * t369;
t322 = qJD(1) * t336 + qJD(3);
t363 = sin(pkin(8));
t475 = t368 * t377;
t480 = t364 * t377;
t300 = t369 * t480 + (-t362 * t373 + t366 * t475) * t365;
t293 = t300 * qJD(1);
t367 = cos(pkin(8));
t491 = t293 * t367;
t398 = t322 * t363 + t491;
t221 = t294 * t376 + t372 * t398;
t295 = t300 * qJD(3);
t282 = qJD(1) * t295;
t296 = t301 * qJD(3);
t283 = qJD(1) * t296;
t477 = t367 * t376;
t156 = qJD(4) * t221 + t282 * t372 + t283 * t477;
t602 = -t156 / 0.2e1;
t220 = -t294 * t372 + t376 * t398;
t478 = t367 * t372;
t155 = qJD(4) * t220 + t282 * t376 - t283 * t478;
t603 = t155 / 0.2e1;
t621 = Ifges(5,4) * t603;
t625 = Ifges(5,2) * t602 + t621;
t355 = qJ(2) * t479;
t541 = pkin(1) * t369;
t442 = qJD(1) * t541;
t332 = qJD(1) * t355 + t362 * t442;
t389 = (t364 * t369 + t368 * t479) * pkin(10);
t284 = qJD(1) * t389 + t332;
t354 = t366 * t442;
t486 = t362 * t365;
t387 = pkin(2) * t369 + (-pkin(10) * t368 - qJ(2)) * t486;
t291 = qJD(1) * t387 + t354;
t318 = (-pkin(10) * t362 * t364 - pkin(2) * t366 - pkin(1)) * t365;
t311 = qJD(1) * t318 + qJD(2);
t214 = -t284 * t373 + t291 * t475 + t311 * t480;
t539 = pkin(11) * t367;
t190 = -t294 * t539 + t214;
t215 = t373 * (t291 * t368 + t311 * t364) + t284 * t377;
t191 = -pkin(11) * t491 - t215;
t492 = t293 * t363;
t246 = pkin(3) * t294 - pkin(11) * t492;
t484 = t363 * t372;
t358 = pkin(11) * t484;
t343 = pkin(3) * t477 - t358;
t457 = t343 * qJD(4) - t376 * t190 - t191 * t478 - t246 * t484;
t490 = t294 * t363;
t623 = pkin(12) * t490 - t457;
t141 = -t191 * t363 + t367 * t246;
t243 = t293 * t372 + t294 * t477;
t244 = t293 * t376 - t294 * t478;
t622 = -pkin(4) * t243 + pkin(12) * t244 - t141 + (pkin(4) * t372 - pkin(12) * t376) * t363 * qJD(4);
t452 = qJD(2) * t365;
t437 = t366 * t452;
t349 = t377 * t437;
t438 = t362 * t452;
t423 = qJD(1) * t438;
t451 = qJD(3) * t377;
t435 = t368 * t451;
t436 = t364 * t451;
t198 = t291 * t435 + t311 * t436 + qJD(1) * t349 + (-qJD(3) * t284 - t368 * t423) * t373;
t494 = t283 * t367;
t159 = -pkin(11) * t494 + t198;
t388 = (-t362 * t475 - t366 * t373) * t365;
t386 = qJD(2) * t388;
t199 = qJD(1) * t386 - qJD(3) * t215;
t160 = -t282 * t539 + t199;
t165 = pkin(11) * t398 + t215;
t166 = pkin(3) * t322 + t190;
t252 = -t291 * t364 + t368 * t311;
t200 = -pkin(3) * t293 - pkin(11) * t490 + t252;
t394 = t364 * t423;
t540 = pkin(11) * t363;
t228 = pkin(3) * t283 - t282 * t540 + t394;
t450 = qJD(4) * t372;
t432 = t367 * t450;
t434 = t363 * t450;
t449 = qJD(4) * t376;
t44 = t376 * (t160 * t367 + t228 * t363) - t372 * t159 - t165 * t449 - t166 * t432 - t200 * t434;
t495 = t283 * t363;
t39 = -pkin(4) * t495 - t44;
t552 = t156 / 0.2e1;
t253 = t322 * t367 + qJD(4) - t492;
t371 = sin(qJ(5));
t375 = cos(qJ(5));
t170 = t221 * t375 + t253 * t371;
t483 = t363 * t375;
t111 = qJD(5) * t170 + t155 * t371 - t283 * t483;
t561 = -t111 / 0.2e1;
t169 = -t221 * t371 + t253 * t375;
t485 = t363 * t371;
t110 = qJD(5) * t169 + t155 * t375 + t283 * t485;
t562 = t110 / 0.2e1;
t403 = t166 * t367 + t200 * t363;
t87 = t165 * t376 + t372 * t403;
t81 = pkin(12) * t253 + t87;
t121 = -t166 * t363 + t367 * t200;
t84 = -pkin(4) * t220 - pkin(12) * t221 + t121;
t34 = t371 * t84 + t375 * t81;
t431 = t367 * t449;
t433 = t363 * t449;
t43 = t376 * t159 + t160 * t478 - t165 * t450 + t166 * t431 + t200 * t433 + t228 * t484;
t38 = pkin(12) * t495 + t43;
t122 = -t160 * t363 + t367 * t228;
t75 = pkin(4) * t156 - pkin(12) * t155 + t122;
t8 = -qJD(5) * t34 - t371 * t38 + t375 * t75;
t616 = t39 * mrSges(6,2) - t8 * mrSges(6,3) + 0.2e1 * Ifges(6,1) * t562 + 0.2e1 * Ifges(6,4) * t561 + 0.2e1 * Ifges(6,5) * t552;
t560 = t111 / 0.2e1;
t430 = t495 / 0.2e1;
t211 = t244 * t371 - t294 * t483;
t342 = t367 * t371 + t372 * t483;
t310 = qJD(5) * t342 + t371 * t433;
t620 = t211 - t310;
t213 = t244 * t375 + t294 * t485;
t341 = -t375 * t367 + t371 * t484;
t309 = -qJD(5) * t341 + t375 * t433;
t619 = t213 - t309;
t618 = t243 - t434;
t315 = qJD(1) * t388;
t453 = qJD(1) * t365;
t316 = (-t362 * t476 + t366 * t377) * t453;
t439 = t362 * t453;
t426 = t364 * t439;
t405 = t363 * t426;
t469 = t376 * t377;
t473 = t372 * t373;
t582 = t367 * t469 - t473;
t460 = t368 * t433 + (t582 * qJD(4) + (-t367 * t473 + t469) * qJD(3)) * t364 - t316 * t376 - (t315 * t367 + t405) * t372;
t274 = -t315 * t363 + t367 * t426;
t594 = qJD(3) * t363 * t481 - t274;
t482 = t363 * t376;
t344 = pkin(3) * t478 + pkin(11) * t482;
t456 = -t344 * qJD(4) + t372 * t190 - t376 * (t191 * t367 + t246 * t363);
t108 = Ifges(6,6) * t111;
t109 = Ifges(6,5) * t110;
t40 = Ifges(6,3) * t156 - t108 + t109;
t447 = qJD(5) * t375;
t448 = qJD(5) * t371;
t7 = t371 * t75 + t375 * t38 + t84 * t447 - t448 * t81;
t422 = -t8 * mrSges(6,1) + t7 * mrSges(6,2);
t617 = t422 + 0.2e1 * Ifges(5,6) * t430 + t43 * mrSges(5,3) - t40 / 0.2e1 - t122 * mrSges(5,1) + t625;
t219 = qJD(5) - t220;
t32 = pkin(13) * t219 + t34;
t370 = sin(qJ(6));
t374 = cos(qJ(6));
t86 = -t372 * t165 + t376 * t403;
t80 = -pkin(4) * t253 - t86;
t47 = -pkin(5) * t169 - pkin(13) * t170 + t80;
t16 = -t32 * t370 + t374 * t47;
t17 = t32 * t374 + t370 * t47;
t586 = t80 * mrSges(6,1) + t16 * mrSges(7,1) - t17 * mrSges(7,2) - t34 * mrSges(6,3);
t168 = qJD(6) - t169;
t519 = t168 * Ifges(7,3);
t127 = t170 * t374 + t219 * t370;
t525 = t127 * Ifges(7,5);
t126 = -t170 * t370 + t219 * t374;
t526 = t126 * Ifges(7,6);
t67 = t519 + t525 + t526;
t513 = t219 * Ifges(6,6);
t518 = t169 * Ifges(6,2);
t533 = Ifges(6,4) * t170;
t97 = t513 + t518 + t533;
t615 = -t67 / 0.2e1 + t97 / 0.2e1 - t586;
t56 = -qJD(6) * t127 - t110 * t370 + t156 * t374;
t53 = Ifges(7,6) * t56;
t55 = qJD(6) * t126 + t110 * t374 + t156 * t370;
t54 = Ifges(7,5) * t55;
t13 = Ifges(7,3) * t111 + t53 + t54;
t20 = pkin(5) * t111 - pkin(13) * t110 + t39;
t5 = pkin(13) * t156 + t7;
t1 = qJD(6) * t16 + t20 * t370 + t374 * t5;
t2 = -qJD(6) * t17 + t20 * t374 - t370 * t5;
t421 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t612 = Ifges(6,6) * t602 - t110 * Ifges(6,4) / 0.2e1;
t614 = -t39 * mrSges(6,1) + t7 * mrSges(6,3) + t421 - t13 / 0.2e1 - Ifges(6,2) * t560 - t612;
t609 = t87 * mrSges(5,3);
t608 = t121 * mrSges(5,1);
t329 = pkin(12) * t367 + t344;
t330 = (-pkin(4) * t376 - pkin(12) * t372 - pkin(3)) * t363;
t501 = -t329 * t448 + t330 * t447 + t371 * t622 - t375 * t623;
t455 = pkin(4) * t490 - t456;
t33 = -t371 * t81 + t375 * t84;
t31 = -pkin(5) * t219 - t33;
t411 = Ifges(7,5) * t374 - Ifges(7,6) * t370;
t530 = Ifges(7,4) * t374;
t413 = -Ifges(7,2) * t370 + t530;
t531 = Ifges(7,4) * t370;
t415 = Ifges(7,1) * t374 - t531;
t416 = mrSges(7,1) * t370 + mrSges(7,2) * t374;
t542 = t374 / 0.2e1;
t543 = -t370 / 0.2e1;
t549 = t168 / 0.2e1;
t556 = t127 / 0.2e1;
t558 = t126 / 0.2e1;
t532 = Ifges(7,4) * t127;
t68 = Ifges(7,2) * t126 + Ifges(7,6) * t168 + t532;
t125 = Ifges(7,4) * t126;
t69 = Ifges(7,1) * t127 + Ifges(7,5) * t168 + t125;
t607 = -t31 * t416 - t411 * t549 - t413 * t558 - t415 * t556 - t542 * t69 - t543 * t68 + (t16 * t374 + t17 * t370) * mrSges(7,3);
t567 = t56 / 0.2e1;
t568 = t55 / 0.2e1;
t606 = -Ifges(6,4) * t562 + Ifges(7,5) * t568 - Ifges(6,2) * t561 - Ifges(6,6) * t552 + Ifges(7,6) * t567 + Ifges(7,3) * t560 - t614;
t601 = -t220 / 0.2e1;
t600 = -t221 / 0.2e1;
t599 = -t253 / 0.2e1;
t598 = t336 / 0.2e1;
t597 = -pkin(13) * t618 + t501;
t596 = -pkin(5) * t620 + pkin(13) * t619 + t455;
t471 = t373 * t376;
t472 = t372 * t377;
t391 = t367 * t472 + t471;
t306 = t364 * t391 + t368 * t484;
t337 = -t363 * t480 + t367 * t368;
t261 = t306 * t371 - t375 * t337;
t462 = -qJD(5) * t261 + t371 * t594 + t375 * t460;
t595 = t315 * t477 - t316 * t372 + t376 * t405 + t368 * t434 + (t391 * qJD(4) + (t367 * t471 + t472) * qJD(3)) * t364;
t514 = t219 * Ifges(6,5);
t576 = -t80 * mrSges(6,2) + t33 * mrSges(6,3);
t167 = Ifges(6,4) * t169;
t516 = t170 * Ifges(6,1);
t98 = t167 + t514 + t516;
t385 = -t98 / 0.2e1 - t514 / 0.2e1 + t576;
t593 = t385 + t607;
t130 = mrSges(6,1) * t219 - mrSges(6,3) * t170;
t82 = -mrSges(7,1) * t126 + mrSges(7,2) * t127;
t496 = -t82 + t130;
t591 = m(6) * t33 - m(7) * t31 + t496;
t563 = Ifges(5,1) * t603 + Ifges(5,4) * t602 + Ifges(5,5) * t430;
t589 = -pkin(12) * qJD(6) * t375 - t87 + t219 * (pkin(5) * t371 - pkin(13) * t375);
t350 = -pkin(5) * t375 - pkin(13) * t371 - pkin(4);
t146 = pkin(4) * t221 - pkin(12) * t220;
t71 = t371 * t146 + t375 * t86;
t588 = pkin(12) * t448 + pkin(13) * t221 - qJD(6) * t350 + t71;
t587 = t199 * mrSges(4,1) - t198 * mrSges(4,2);
t119 = -mrSges(6,1) * t169 + mrSges(6,2) * t170;
t175 = mrSges(5,1) * t253 - mrSges(5,3) * t221;
t584 = t175 - t119;
t583 = t375 * t329 + t371 * t330;
t218 = Ifges(5,4) * t220;
t507 = t253 * Ifges(5,5);
t509 = t221 * Ifges(5,1);
t135 = t218 + t507 + t509;
t553 = -t135 / 0.2e1;
t581 = -t218 / 0.2e1 + t86 * mrSges(5,3) + t553 - t121 * mrSges(5,2) - t507 / 0.2e1;
t454 = t362 * t541 + t355;
t298 = t389 + t454;
t289 = t377 * t298;
t357 = t366 * t541;
t303 = t357 + t387;
t224 = t303 * t476 + t318 * t481 + t289;
t397 = t300 * t367 + t336 * t363;
t189 = pkin(11) * t397 + t224;
t223 = -t298 * t373 + t303 * t475 + t318 * t480;
t194 = pkin(3) * t336 - t301 * t539 + t223;
t254 = -t303 * t364 + t368 * t318;
t209 = -pkin(3) * t300 - t301 * t540 + t254;
t99 = -t372 * t189 + t376 * (t194 * t367 + t209 * t363);
t205 = t303 * t435 + t318 * t436 + t349 + (-qJD(3) * t298 - t368 * t438) * t373;
t178 = -t296 * t539 + t205;
t206 = t386 + (-t289 + (-t303 * t368 - t318 * t364) * t373) * qJD(3);
t179 = -t295 * t539 + t206;
t425 = t364 * t438;
t239 = pkin(3) * t296 - t295 * t540 + t425;
t58 = t376 * (t179 * t367 + t239 * t363) - t372 * t178 - t189 * t449 - t194 * t432 - t209 * t434;
t488 = t296 * t363;
t57 = t376 * t178 + t179 * t478 - t189 * t450 + t194 * t431 + t209 * t433 + t239 * t484;
t49 = pkin(12) * t488 + t57;
t100 = t376 * t189 + t194 * t478 + t209 * t484;
t260 = -t300 * t363 + t336 * t367;
t90 = pkin(12) * t260 + t100;
t128 = -t194 * t363 + t367 * t209;
t487 = t301 * t372;
t231 = -t300 * t477 - t336 * t482 + t487;
t232 = t301 * t376 + t372 * t397;
t93 = pkin(4) * t231 - pkin(12) * t232 + t128;
t535 = t371 * t93 + t375 * t90;
t136 = -t179 * t363 + t367 * t239;
t171 = -t296 * t478 + t295 * t376 + (t376 * t397 - t487) * qJD(4);
t172 = qJD(4) * t232 + t295 * t372 + t296 * t477;
t79 = pkin(4) * t172 - pkin(12) * t171 + t136;
t12 = -qJD(5) * t535 - t371 * t49 + t375 * t79;
t500 = -qJD(5) * t583 + t371 * t623 + t375 * t622;
t555 = Ifges(5,5) * t600 + Ifges(5,6) * t601 + Ifges(5,3) * t599;
t14 = t55 * Ifges(7,4) + t56 * Ifges(7,2) + t111 * Ifges(7,6);
t572 = t14 / 0.2e1;
t15 = t55 * Ifges(7,1) + t56 * Ifges(7,4) + t111 * Ifges(7,5);
t571 = t15 / 0.2e1;
t566 = -t68 / 0.2e1;
t512 = t219 * Ifges(6,3);
t515 = t170 * Ifges(6,5);
t517 = t169 * Ifges(6,6);
t96 = t512 + t515 + t517;
t565 = t96 / 0.2e1;
t152 = Ifges(5,6) * t156;
t153 = Ifges(5,5) * t155;
t104 = Ifges(5,3) * t495 - t152 + t153;
t564 = t104 / 0.2e1;
t559 = -t126 / 0.2e1;
t557 = -t127 / 0.2e1;
t506 = t253 * Ifges(5,6);
t511 = t220 * Ifges(5,2);
t534 = Ifges(5,4) * t221;
t134 = t506 + t511 + t534;
t554 = -t134 / 0.2e1;
t551 = -t167 / 0.2e1;
t550 = -t168 / 0.2e1;
t548 = t169 / 0.2e1;
t547 = t170 / 0.2e1;
t546 = t219 / 0.2e1;
t544 = t294 / 0.2e1;
t538 = t33 * mrSges(6,1);
t537 = t34 * mrSges(6,2);
t21 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t76 = mrSges(6,1) * t156 - mrSges(6,3) * t110;
t536 = t21 - t76;
t504 = t294 * Ifges(4,4);
t503 = pkin(5) * t618 - t500;
t328 = t358 + (-pkin(3) * t376 - pkin(4)) * t367;
t263 = pkin(5) * t341 - pkin(13) * t342 + t328;
t265 = -pkin(13) * t482 + t583;
t217 = t263 * t370 + t265 * t374;
t502 = -qJD(6) * t217 - t370 * t597 + t374 * t596;
t216 = t263 * t374 - t265 * t370;
t499 = qJD(6) * t216 + t370 * t596 + t374 * t597;
t498 = t370 * t589 - t374 * t588;
t497 = t370 * t588 + t374 * t589;
t474 = t370 * t375;
t470 = t374 * t375;
t262 = t306 * t375 + t337 * t371;
t305 = -t364 * t582 - t368 * t482;
t234 = -t262 * t370 + t305 * t374;
t468 = qJD(6) * t234 + t370 * t595 + t374 * t462;
t235 = t262 * t374 + t305 * t370;
t467 = -qJD(6) * t235 - t370 * t462 + t374 * t595;
t147 = -t213 * t370 + t243 * t374;
t393 = -t342 * t374 + t370 * t482;
t251 = qJD(6) * t393 - t309 * t370 + t374 * t434;
t466 = t147 - t251;
t148 = t213 * t374 + t243 * t370;
t312 = -t342 * t370 - t374 * t482;
t250 = qJD(6) * t312 + t309 * t374 + t370 * t434;
t465 = t148 - t250;
t458 = Ifges(4,5) * t282 - Ifges(4,6) * t283;
t420 = t44 * mrSges(5,1) - t43 * mrSges(5,2);
t418 = t1 * t374 - t2 * t370;
t417 = mrSges(7,1) * t374 - mrSges(7,2) * t370;
t414 = Ifges(7,1) * t370 + t530;
t412 = Ifges(7,2) * t374 + t531;
t410 = Ifges(7,5) * t370 + Ifges(7,6) * t374;
t36 = pkin(13) * t231 + t535;
t192 = t232 * t371 - t375 * t260;
t193 = t232 * t375 + t260 * t371;
t89 = -pkin(4) * t260 - t99;
t64 = pkin(5) * t192 - pkin(13) * t193 + t89;
t19 = t36 * t374 + t370 * t64;
t18 = -t36 * t370 + t374 * t64;
t45 = -t371 * t90 + t375 * t93;
t70 = t146 * t375 - t371 * t86;
t138 = t193 * t374 + t231 * t370;
t137 = -t193 * t370 + t231 * t374;
t269 = -t329 * t371 + t330 * t375;
t11 = t371 * t79 + t375 * t49 + t93 * t447 - t448 * t90;
t339 = (mrSges(3,1) * t369 - mrSges(3,3) * t486) * qJD(1);
t340 = (-mrSges(3,2) * t369 + mrSges(3,3) * t479) * qJD(1);
t50 = -pkin(4) * t488 - t58;
t382 = -t519 / 0.2e1 + t513 / 0.2e1 - t526 / 0.2e1 - t525 / 0.2e1 + t533 / 0.2e1 + t615;
t381 = t608 + t538 + t554 + t565 - t534 / 0.2e1 + t517 / 0.2e1 + t515 / 0.2e1 + t512 / 0.2e1 - t506 / 0.2e1 - t537 - t609;
t380 = t518 / 0.2e1 + t382;
t331 = -qJ(2) * t439 + t354;
t326 = pkin(12) * t470 + t350 * t370;
t325 = -pkin(12) * t474 + t350 * t374;
t290 = Ifges(4,4) * t293;
t264 = pkin(5) * t482 - t269;
t256 = mrSges(4,1) * t322 - mrSges(4,3) * t294;
t255 = -mrSges(4,2) * t322 + mrSges(4,3) * t293;
t249 = -mrSges(4,1) * t293 + mrSges(4,2) * t294;
t241 = mrSges(4,1) * t283 + mrSges(4,2) * t282;
t230 = t294 * Ifges(4,1) + t322 * Ifges(4,5) + t290;
t229 = t293 * Ifges(4,2) + t322 * Ifges(4,6) + t504;
t174 = -mrSges(5,2) * t253 + mrSges(5,3) * t220;
t145 = -mrSges(5,1) * t220 + mrSges(5,2) * t221;
t143 = t220 * t470 + t221 * t370;
t142 = -t220 * t474 + t221 * t374;
t140 = -mrSges(5,2) * t495 - mrSges(5,3) * t156;
t139 = mrSges(5,1) * t495 - mrSges(5,3) * t155;
t129 = -mrSges(6,2) * t219 + mrSges(6,3) * t169;
t120 = pkin(5) * t170 - pkin(13) * t169;
t118 = qJD(5) * t193 + t171 * t371 - t296 * t483;
t117 = -qJD(5) * t192 + t171 * t375 + t296 * t485;
t116 = mrSges(5,1) * t156 + mrSges(5,2) * t155;
t95 = mrSges(7,1) * t168 - mrSges(7,3) * t127;
t94 = -mrSges(7,2) * t168 + mrSges(7,3) * t126;
t77 = -mrSges(6,2) * t156 - mrSges(6,3) * t111;
t66 = -qJD(6) * t138 - t117 * t370 + t172 * t374;
t65 = qJD(6) * t137 + t117 * t374 + t172 * t370;
t63 = mrSges(6,1) * t111 + mrSges(6,2) * t110;
t61 = -pkin(5) * t221 - t70;
t35 = -pkin(5) * t231 - t45;
t30 = -mrSges(7,2) * t111 + mrSges(7,3) * t56;
t29 = mrSges(7,1) * t111 - mrSges(7,3) * t55;
t28 = t120 * t370 + t33 * t374;
t27 = t120 * t374 - t33 * t370;
t22 = pkin(5) * t118 - pkin(13) * t117 + t50;
t10 = -pkin(5) * t172 - t12;
t9 = pkin(13) * t172 + t11;
t6 = -pkin(5) * t156 - t8;
t4 = -qJD(6) * t19 + t22 * t374 - t370 * t9;
t3 = qJD(6) * t18 + t22 * t370 + t374 * t9;
t23 = [t322 * (Ifges(4,5) * t295 - Ifges(4,6) * t296) / 0.2e1 + t293 * (Ifges(4,4) * t295 - Ifges(4,2) * t296) / 0.2e1 + t252 * (mrSges(4,1) * t296 + mrSges(4,2) * t295) + (t198 * t300 - t199 * t301 - t214 * t295 - t215 * t296) * mrSges(4,3) + (Ifges(4,1) * t295 - Ifges(4,4) * t296) * t544 + t253 * (Ifges(5,5) * t171 - Ifges(5,6) * t172 + Ifges(5,3) * t488) / 0.2e1 + t220 * (Ifges(5,4) * t171 - Ifges(5,2) * t172 + Ifges(5,6) * t488) / 0.2e1 + t221 * (Ifges(5,1) * t171 - Ifges(5,4) * t172 + Ifges(5,5) * t488) / 0.2e1 + t86 * (mrSges(5,1) * t488 - mrSges(5,3) * t171) + (-Ifges(6,4) * t547 + Ifges(7,5) * t556 - Ifges(6,2) * t548 - Ifges(6,6) * t546 + Ifges(7,6) * t558 + Ifges(7,3) * t549 - t615) * t118 + (t117 * t80 - t172 * t34) * mrSges(6,2) - t172 * t609 + (t121 * t171 + t122 * t232 - t260 * t43 - t488 * t87) * mrSges(5,2) + (-mrSges(4,1) * t300 + mrSges(4,2) * t301) * t394 + t172 * t608 + (Ifges(6,1) * t117 + Ifges(6,5) * t172) * t547 + t18 * t29 + t19 * t30 + (Ifges(6,5) * t562 + Ifges(6,6) * t561 + Ifges(6,3) * t552 - t617 - t625) * t231 + t535 * t77 + t606 * t192 + t249 * t425 + (Ifges(7,5) * t65 + Ifges(7,6) * t66) * t549 + (Ifges(7,5) * t138 + Ifges(7,6) * t137) * t560 + (Ifges(7,4) * t65 + Ifges(7,2) * t66) * t558 + (Ifges(7,4) * t138 + Ifges(7,2) * t137) * t567 + (Ifges(5,4) * t232 + Ifges(5,6) * t260) * t602 + (t1 * t137 - t138 * t2 - t16 * t65 + t17 * t66) * mrSges(7,3) + (Ifges(6,5) * t117 + Ifges(6,3) * t172) * t546 + m(6) * (t11 * t34 + t12 * t33 + t39 * t89 + t45 * t8 + t50 * t80 + t535 * t7) + m(5) * (t100 * t43 + t121 * t136 + t122 * t128 + t44 * t99 + t57 * t87 + t58 * t86) + m(7) * (t1 * t19 + t10 * t31 + t16 * t4 + t17 * t3 + t18 * t2 + t35 * t6) + (Ifges(5,5) * t232 + Ifges(5,3) * t260) * t430 + m(3) * ((-t331 * t362 + t332 * t366) * t365 + (t454 * t479 - (-qJ(2) * t486 + t357) * t486) * qJD(1)) * qJD(2) + 0.2e1 * t340 * t437 - 0.2e1 * t339 * t438 + (Ifges(6,4) * t117 + Ifges(6,6) * t172) * t548 - t296 * t229 / 0.2e1 + t295 * t230 / 0.2e1 + t44 * (mrSges(5,1) * t260 - mrSges(5,3) * t232) + t254 * t241 + t205 * t255 + t206 * t256 + (Ifges(5,1) * t232 + Ifges(5,5) * t260) * t603 + (t458 / 0.2e1 + t587) * t336 + m(4) * (t198 * t224 + t199 * t223 + t205 * t215 + t206 * t214 + (qJD(1) * t254 + t252) * t425) + t35 * t21 + t616 * t193 + t172 * t554 - t488 * t555 + t232 * t563 + t260 * t564 + t172 * t565 + t138 * t571 + t137 * t572 - (t224 * mrSges(4,3) + Ifges(4,4) * t301 + Ifges(4,2) * t300 + Ifges(4,6) * t598) * t283 + (-t223 * mrSges(4,3) + Ifges(4,1) * t301 + Ifges(4,4) * t300 + Ifges(4,5) * t598) * t282 + t31 * (-mrSges(7,1) * t66 + mrSges(7,2) * t65) + t66 * t68 / 0.2e1 + t65 * t69 / 0.2e1 + t45 * t76 + t10 * t82 + t89 * t63 + t3 * t94 + t4 * t95 + (Ifges(7,1) * t65 + Ifges(7,4) * t66) * t556 + (Ifges(7,1) * t138 + Ifges(7,4) * t137) * t568 + t117 * t98 / 0.2e1 + t50 * t119 + t128 * t116 + t11 * t129 + t12 * t130 + t6 * (-mrSges(7,1) * t137 + mrSges(7,2) * t138) + t99 * t139 + t100 * t140 + t136 * t145 + t171 * t135 / 0.2e1 + t33 * (mrSges(6,1) * t172 - mrSges(6,3) * t117) + t57 * t174 + t58 * t175; t536 * t261 + t467 * t95 + (t1 * t235 + t16 * t467 + t17 * t468 + t2 * t234 + t261 * t6) * m(7) + t468 * t94 + t460 * t174 + (t121 * t594 + t122 * t337 - t305 * t44 + t306 * t43 + t460 * t87) * m(5) + t462 * t129 + (-t261 * t8 + t262 * t7 + t305 * t39 + t34 * t462) * m(6) + (-t214 * t315 - t215 * t316 + (t198 * t373 + t199 * t377 + (-t214 * t373 + t215 * t377) * qJD(3) + (qJD(2) * t368 - t252) * t439) * t364) * m(4) + (-t249 * t439 + (-t282 * t377 - t283 * t373) * mrSges(4,3) + (t255 * t377 + (t145 * t363 - t256) * t373) * qJD(3)) * t364 + ((-m(3) * t332 - t340) * t366 + (m(3) * t331 + t339) * t362) * t453 + t368 * t241 + t337 * t116 - t315 * t256 - t316 * t255 + t306 * t140 - t274 * t145 + t262 * t77 + t234 * t29 + t235 * t30 + (t63 - t139) * t305 + t591 * (-qJD(5) * t262 - t371 * t460 + t375 * t594) + (-m(5) * t86 + m(6) * t80 - t584) * t595; (t269 * t8 + t328 * t39 + t33 * t500 + t34 * t501 + t455 * t80 + t583 * t7) * m(6) + t583 * t77 + (-Ifges(7,4) * t393 + Ifges(7,2) * t312) * t567 + (-Ifges(7,1) * t393 + Ifges(7,4) * t312) * t568 + (-Ifges(7,5) * t393 + Ifges(7,6) * t312) * t560 + (t1 * t312 + t16 * t465 - t17 * t466 + t2 * t393) * mrSges(7,3) + t6 * (-mrSges(7,1) * t312 - mrSges(7,2) * t393) - t393 * t571 + t616 * t342 + t576 * t619 - t586 * t620 + (t214 * t293 + t215 * t294) * mrSges(4,3) + (t67 - t97) * (t310 / 0.2e1 - t211 / 0.2e1) - (-Ifges(4,2) * t294 + t230 + t290) * t293 / 0.2e1 + t587 + (t153 / 0.2e1 - t152 / 0.2e1 + t564 + t420) * t367 + (Ifges(5,5) * t244 - Ifges(5,6) * t243) * t599 + (Ifges(5,1) * t244 - Ifges(5,4) * t243) * t600 + (Ifges(5,4) * t244 - Ifges(5,2) * t243) * t601 + t606 * t341 + (t243 * t87 + t244 * t86) * mrSges(5,3) + t458 + (t251 / 0.2e1 - t147 / 0.2e1) * t68 + (t250 / 0.2e1 - t148 / 0.2e1) * t69 - t243 * t538 - t294 * (Ifges(4,1) * t293 - t504) / 0.2e1 + t499 * t94 + t500 * t130 + t501 * t129 + t502 * t95 + t503 * t82 + (t1 * t217 + t16 * t502 + t17 * t499 + t2 * t216 + t264 * t6 + t31 * t503) * m(7) + (mrSges(7,1) * t466 - mrSges(7,2) * t465) * t31 + t455 * t119 + t456 * t175 + (-pkin(3) * t122 * t363 - t121 * t141 + t343 * t44 + t344 * t43 + t456 * t86 + t457 * t87) * m(5) + t457 * t174 + t343 * t139 + t344 * t140 - t322 * (Ifges(4,5) * t293 - Ifges(4,6) * t294) / 0.2e1 + t328 * t63 - t252 * (mrSges(4,1) * t294 + mrSges(4,2) * t293) + t264 * t21 + t269 * t76 - t214 * t255 + t215 * t256 + t243 * t134 / 0.2e1 - t121 * (mrSges(5,1) * t243 + mrSges(5,2) * t244) - t219 * (Ifges(6,5) * t213 - Ifges(6,6) * t211 + Ifges(6,3) * t243) / 0.2e1 - t170 * (Ifges(6,1) * t213 - Ifges(6,4) * t211 + Ifges(6,5) * t243) / 0.2e1 - t169 * (Ifges(6,4) * t213 - Ifges(6,2) * t211 + Ifges(6,6) * t243) / 0.2e1 - t243 * t96 / 0.2e1 + t216 * t29 + t217 * t30 + (Ifges(5,3) * t494 / 0.2e1 - pkin(3) * t116 + (t122 * mrSges(5,2) - t44 * mrSges(5,3) + 0.2e1 * t563) * t372 + (-t109 / 0.2e1 + t108 / 0.2e1 + t621 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1) * t156 + t617) * t376 + (-t86 * mrSges(5,1) + t87 * mrSges(5,2) + 0.2e1 * t555) * t294 + ((t509 / 0.2e1 - t581) * t376 + (-t511 / 0.2e1 + t381) * t372) * qJD(4)) * t363 + t243 * t537 + t229 * t544 + (Ifges(6,5) * t309 - Ifges(6,6) * t310) * t546 + (Ifges(6,1) * t309 - Ifges(6,4) * t310) * t547 + (Ifges(6,4) * t309 - Ifges(6,2) * t310) * t548 + (Ifges(7,5) * t250 + Ifges(7,6) * t251 + Ifges(7,3) * t310) * t549 + (Ifges(7,5) * t148 + Ifges(7,6) * t147 + Ifges(7,3) * t211) * t550 + t244 * t553 + (Ifges(7,1) * t250 + Ifges(7,4) * t251 + Ifges(7,5) * t310) * t556 + (Ifges(7,1) * t148 + Ifges(7,4) * t147 + Ifges(7,5) * t211) * t557 + (Ifges(7,4) * t250 + Ifges(7,2) * t251 + Ifges(7,6) * t310) * t558 + (Ifges(7,4) * t148 + Ifges(7,2) * t147 + Ifges(7,6) * t211) * t559 + t312 * t572 + (t309 / 0.2e1 - t213 / 0.2e1) * t98 - t141 * t145; (-t54 / 0.2e1 - t53 / 0.2e1 + (-Ifges(6,2) / 0.2e1 - Ifges(7,3) / 0.2e1) * t111 + (m(6) * t7 + t77) * pkin(12) + (t551 - t516 / 0.2e1 + t385) * t220 - t612 + t614) * t375 + (t415 * t568 + t413 * t567 + t411 * t560 + t6 * t416 + t15 * t542 + t14 * t543 + (-t1 * t370 - t2 * t374) * mrSges(7,3) + (-m(6) * t8 + m(7) * t6 + t536) * pkin(12) + t380 * t220 + (t414 * t557 + t410 * t550 + t412 * t559 + t69 * t543 + t374 * t566 + t31 * t417 + (t16 * t370 - t17 * t374) * mrSges(7,3)) * qJD(6) + t616) * t371 + (-t142 * t17 + t143 * t16) * mrSges(7,3) + t584 * t87 + t581 * t220 + t420 + ((-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t220 - t381) * t221 + t104 + (t1 * t326 + t16 * t497 + t17 * t498 + t2 * t325 - t31 * t61) * m(7) + t498 * t94 + t497 * t95 + (-pkin(4) * t39 - t33 * t70 - t34 * t71 - t80 * t87) * m(6) + t325 * t29 + t326 * t30 + (((-m(6) * t34 - t129) * pkin(12) - t380) * t371 + (t516 / 0.2e1 - t591 * pkin(12) + t167 / 0.2e1 - t593) * t375) * qJD(5) + (Ifges(7,5) * t143 + Ifges(7,6) * t142) * t550 + (Ifges(7,1) * t143 + Ifges(7,4) * t142) * t557 + (Ifges(7,4) * t143 + Ifges(7,2) * t142) * t559 + t142 * t566 - pkin(4) * t63 - t61 * t82 - t71 * t129 - t70 * t130 - t143 * t69 / 0.2e1 - t31 * (-mrSges(7,1) * t142 + mrSges(7,2) * t143) - t86 * t174; (-t29 * t370 + t30 * t374) * pkin(13) - t422 + t40 + t382 * t170 - pkin(5) * t21 + t14 * t542 + t370 * t571 + t496 * t34 + (((-m(7) * t16 - t95) * t374 + (-m(7) * t17 - t94) * t370) * pkin(13) - t607) * qJD(6) + (t551 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t170 + t593) * t169 + t410 * t560 + t412 * t567 + t414 * t568 - t28 * t94 - t27 * t95 - t6 * t417 + t418 * mrSges(7,3) + (-pkin(5) * t6 + pkin(13) * t418 - t16 * t27 - t17 * t28 - t31 * t34) * m(7) - t33 * t129; -t31 * (mrSges(7,1) * t127 + mrSges(7,2) * t126) - t16 * t94 + t17 * t95 + (Ifges(7,1) * t126 - t532) * t557 + t68 * t556 + (Ifges(7,5) * t126 - Ifges(7,6) * t127) * t550 + (t126 * t16 + t127 * t17) * mrSges(7,3) - t421 + t13 + (-Ifges(7,2) * t127 + t125 + t69) * t559;];
tauc  = t23(:);
