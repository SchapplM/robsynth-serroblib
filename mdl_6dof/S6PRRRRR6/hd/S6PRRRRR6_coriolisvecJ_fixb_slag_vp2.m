% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:36:13
% EndTime: 2018-11-23 15:36:40
% DurationCPUTime: 26.98s
% Computational Cost: add. (28066->914), mult. (83287->1317), div. (0->0), fcn. (70347->16), ass. (0->448)
t348 = sin(qJ(2));
t340 = sin(pkin(6));
t454 = qJD(1) * t340;
t434 = t348 * t454;
t339 = sin(pkin(7));
t452 = qJD(2) * t339;
t314 = pkin(10) * t452 + t434;
t347 = sin(qJ(3));
t341 = cos(pkin(8));
t528 = pkin(11) * t341;
t440 = t339 * t528;
t353 = cos(qJ(2));
t433 = t353 * t454;
t322 = qJD(2) * pkin(2) + t433;
t352 = cos(qJ(3));
t343 = cos(pkin(6));
t453 = qJD(1) * t343;
t435 = t339 * t453;
t323 = t352 * t435;
t342 = cos(pkin(7));
t471 = t342 * t352;
t455 = t322 * t471 + t323;
t208 = (-qJD(2) * t440 - t314) * t347 + t455;
t472 = t342 * t347;
t480 = t314 * t352;
t380 = -t322 * t472 - t480;
t451 = qJD(2) * t352;
t430 = t341 * t451;
t209 = (-pkin(11) * t430 - t347 * t453) * t339 + t380;
t338 = sin(pkin(8));
t529 = pkin(11) * t338;
t383 = pkin(3) * t347 - t352 * t529;
t285 = t383 * t452;
t346 = sin(qJ(4));
t479 = t338 * t346;
t331 = pkin(11) * t479;
t351 = cos(qJ(4));
t473 = t341 * t351;
t309 = pkin(3) * t473 - t331;
t474 = t341 * t346;
t570 = t309 * qJD(4) - t351 * t208 - t209 * t474 - t285 * t479;
t432 = t347 * t452;
t417 = t338 * t432;
t613 = pkin(12) * t417 - t570;
t162 = -t209 * t338 + t341 * t285;
t466 = t347 * t351;
t468 = t346 * t352;
t378 = t341 * t466 + t468;
t271 = t378 * t452;
t462 = t351 * t352;
t469 = t346 * t347;
t376 = -t341 * t469 + t462;
t273 = t376 * t452;
t448 = qJD(4) * t338;
t612 = -pkin(4) * t271 + pkin(12) * t273 - t162 + (pkin(4) * t346 - pkin(12) * t351) * t448;
t330 = qJD(2) * t342 + qJD(3);
t379 = t341 * t462 - t469;
t478 = t338 * t351;
t235 = t330 * t478 + t379 * t452;
t234 = qJD(5) - t235;
t335 = pkin(2) * t472;
t476 = t339 * t352;
t312 = pkin(10) * t476 + t335;
t246 = (t338 * t342 + t341 * t476) * pkin(11) + t312;
t336 = pkin(2) * t471;
t410 = t339 * (-pkin(10) - t528);
t384 = t347 * t410;
t256 = pkin(3) * t342 + t336 + t384;
t428 = qJD(3) * t471;
t328 = pkin(2) * t428;
t257 = qJD(3) * t384 + t328;
t464 = t348 * t352;
t465 = t347 * t353;
t374 = -t342 * t464 - t465;
t272 = t374 * t454;
t461 = t352 * t353;
t467 = t347 * t348;
t372 = -t342 * t467 + t461;
t274 = t372 * t454;
t382 = -pkin(3) * t352 - t347 * t529;
t281 = (-pkin(2) + t382) * t339;
t419 = t339 * t434;
t393 = t338 * t419;
t447 = qJD(4) * t346;
t425 = t341 * t447;
t427 = t338 * t447;
t446 = qJD(4) * t351;
t611 = -t246 * t446 - t256 * t425 - t272 * t473 - t281 * t427 - t351 * t393 + (-t257 + t274) * t346;
t311 = pkin(3) * t474 + pkin(11) * t478;
t610 = -t311 * qJD(4) + t346 * t208;
t345 = sin(qJ(5));
t350 = cos(qJ(5));
t230 = t273 * t345 - t350 * t417;
t308 = t341 * t345 + t350 * t479;
t426 = t338 * t446;
t260 = qJD(5) * t308 + t345 * t426;
t609 = t230 - t260;
t231 = t273 * t350 + t345 * t417;
t307 = -t350 * t341 + t345 * t479;
t259 = -qJD(5) * t307 + t350 * t426;
t608 = t231 - t259;
t607 = t271 - t427;
t377 = t341 * t468 + t466;
t360 = (qJD(3) * t378 + qJD(4) * t377) * t339;
t194 = qJD(2) * t360 + t330 * t427;
t595 = -t194 / 0.2e1;
t361 = (qJD(3) * t376 + qJD(4) * t379) * t339;
t204 = t342 * t426 + t361;
t205 = t342 * t427 + t360;
t258 = (t352 * t410 - t335) * qJD(3);
t368 = t383 * qJD(3);
t286 = t339 * t368;
t210 = -t258 * t338 + t341 * t286;
t104 = pkin(4) * t205 - pkin(12) * t204 + t210;
t202 = -t256 * t338 + t341 * t281;
t249 = -t339 * t379 - t342 * t478;
t367 = t377 * t339;
t251 = t342 * t479 + t367;
t130 = pkin(4) * t249 - pkin(12) * t251 + t202;
t148 = t351 * t246 + t256 * t474 + t281 * t479;
t303 = -t338 * t476 + t341 * t342;
t135 = pkin(12) * t303 + t148;
t176 = t274 * t351 + (t272 * t341 + t393) * t346;
t228 = -t272 * t338 + t341 * t419;
t444 = qJD(5) * t350;
t445 = qJD(5) * t345;
t449 = qJD(3) * t339;
t429 = t347 * t449;
t415 = t338 * t429;
t424 = t341 * t446;
t92 = -t246 * t447 + t256 * t424 + t351 * t257 + t258 * t474 + t281 * t426 + t286 * t479;
t87 = pkin(12) * t415 + t92;
t579 = t130 * t444 - t135 * t445 + (-t176 + t87) * t350 + (t104 - t228) * t345;
t291 = pkin(12) * t341 + t311;
t292 = (-pkin(4) * t351 - pkin(12) * t346 - pkin(3)) * t338;
t577 = -t291 * t445 + t292 * t444 + t612 * t345 - t350 * t613;
t575 = -t258 * t473 + (-pkin(4) * t429 - t286 * t351) * t338 - t611;
t572 = t209 * t473 - (-pkin(4) * t432 - t285 * t351) * t338 - t610;
t366 = t374 * qJD(2);
t365 = t340 * t366;
t161 = qJD(1) * t365 + qJD(3) * t209;
t244 = (t368 + t434) * t452;
t422 = qJD(2) * t449;
t413 = t347 * t422;
t418 = t342 * t434;
t436 = qJD(3) * t323 + t322 * t428 + t433 * t451;
t450 = qJD(3) * t314;
t160 = (-t450 + (-qJD(3) * t440 - t418) * qJD(2)) * t347 + t436;
t223 = t480 + (t322 * t342 + t435) * t347;
t190 = (t330 * t338 + t339 * t430) * pkin(11) + t223;
t195 = pkin(3) * t330 + t208;
t329 = t342 * t453;
t227 = t329 + (qJD(2) * t382 - t322) * t339;
t421 = -t346 * t160 - t190 * t446 - t195 * t425 - t227 * t427;
t39 = -t161 * t473 + (-pkin(4) * t413 - t244 * t351) * t338 - t421;
t534 = t194 / 0.2e1;
t236 = qJD(2) * t367 + t330 * t479;
t431 = t339 * t451;
t277 = t330 * t341 - t338 * t431 + qJD(4);
t192 = t236 * t350 + t277 * t345;
t193 = qJD(2) * t361 + t330 * t426;
t385 = t338 * t413;
t112 = qJD(5) * t192 + t193 * t345 - t350 * t385;
t549 = -t112 / 0.2e1;
t191 = -t236 * t345 + t277 * t350;
t111 = qJD(5) * t191 + t193 * t350 + t345 * t385;
t550 = t111 / 0.2e1;
t555 = Ifges(6,1) * t550 + Ifges(6,4) * t549 + Ifges(6,5) * t534;
t606 = t39 * mrSges(6,2) + 0.2e1 * t555;
t390 = t195 * t341 + t227 * t338;
t97 = t190 * t351 + t346 * t390;
t605 = -pkin(12) * qJD(6) * t350 - t97 + t234 * (pkin(5) * t345 - pkin(13) * t350);
t325 = -pkin(5) * t350 - pkin(13) * t345 - pkin(4);
t170 = pkin(4) * t236 - pkin(12) * t235;
t96 = -t346 * t190 + t351 * t390;
t69 = t345 * t170 + t350 * t96;
t604 = pkin(12) * t445 + pkin(13) * t236 - qJD(6) * t325 + t69;
t548 = t112 / 0.2e1;
t602 = -pkin(13) * t205 - t579;
t601 = -pkin(13) * t607 + t577;
t600 = -pkin(5) * t609 + pkin(13) * t608 + t572;
t212 = t251 * t345 - t350 * t303;
t126 = -qJD(5) * t212 + t204 * t350 + t345 * t415;
t213 = t251 * t350 + t303 * t345;
t127 = qJD(5) * t213 + t204 * t345 - t350 * t415;
t599 = pkin(5) * t127 - pkin(13) * t126 + t575;
t85 = pkin(12) * t277 + t97;
t136 = -t195 * t338 + t341 * t227;
t95 = -pkin(4) * t235 - pkin(12) * t236 + t136;
t42 = -t345 * t85 + t350 * t95;
t35 = -pkin(5) * t234 - t42;
t344 = sin(qJ(6));
t349 = cos(qJ(6));
t43 = t345 * t95 + t350 * t85;
t36 = pkin(13) * t234 + t43;
t84 = -pkin(4) * t277 - t96;
t47 = -pkin(5) * t191 - pkin(13) * t192 + t84;
t14 = -t344 * t36 + t349 * t47;
t15 = t344 * t47 + t349 * t36;
t397 = t14 * t349 + t15 * t344;
t399 = Ifges(7,5) * t349 - Ifges(7,6) * t344;
t510 = Ifges(7,4) * t349;
t401 = -Ifges(7,2) * t344 + t510;
t511 = Ifges(7,4) * t344;
t403 = Ifges(7,1) * t349 - t511;
t404 = mrSges(7,1) * t344 + mrSges(7,2) * t349;
t530 = t349 / 0.2e1;
t531 = -t344 / 0.2e1;
t186 = qJD(6) - t191;
t537 = t186 / 0.2e1;
t141 = t192 * t349 + t234 * t344;
t542 = t141 / 0.2e1;
t140 = -t192 * t344 + t234 * t349;
t544 = t140 / 0.2e1;
t512 = Ifges(7,4) * t141;
t60 = t140 * Ifges(7,2) + t186 * Ifges(7,6) + t512;
t137 = Ifges(7,4) * t140;
t61 = t141 * Ifges(7,1) + t186 * Ifges(7,5) + t137;
t598 = -t397 * mrSges(7,3) + t35 * t404 + t399 * t537 + t401 * t544 + t403 * t542 + t530 * t61 + t531 * t60;
t597 = Ifges(6,6) * t595 - t111 * Ifges(6,4) / 0.2e1;
t596 = t193 / 0.2e1;
t594 = t277 / 0.2e1;
t567 = t350 * t291 + t345 * t292;
t576 = -qJD(5) * t567 + t345 * t613 + t612 * t350;
t593 = t344 * t604 + t349 * t605;
t592 = t344 * t605 - t349 * t604;
t185 = Ifges(6,4) * t191;
t494 = t234 * Ifges(6,5);
t501 = t192 * Ifges(6,1);
t102 = t185 + t494 + t501;
t582 = t84 * mrSges(6,2);
t369 = -t102 / 0.2e1 - t494 / 0.2e1 - t582;
t583 = t42 * mrSges(6,3);
t591 = t369 + t583 - t598;
t51 = -qJD(6) * t141 - t111 * t344 + t194 * t349;
t48 = Ifges(7,6) * t51;
t50 = qJD(6) * t140 + t111 * t349 + t194 * t344;
t49 = Ifges(7,5) * t50;
t11 = Ifges(7,3) * t112 + t48 + t49;
t16 = pkin(5) * t112 - pkin(13) * t111 + t39;
t40 = t351 * t160 + t161 * t474 - t190 * t447 + t195 * t424 + t227 * t426 + t244 * t479;
t38 = pkin(12) * t385 + t40;
t125 = -t161 * t338 + t341 * t244;
t72 = pkin(4) * t194 - pkin(12) * t193 + t125;
t7 = t345 * t72 + t350 * t38 + t95 * t444 - t445 * t85;
t5 = pkin(13) * t194 + t7;
t1 = qJD(6) * t14 + t16 * t344 + t349 * t5;
t2 = -qJD(6) * t15 + t16 * t349 - t344 * t5;
t411 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t590 = t39 * mrSges(6,1) + t11 / 0.2e1 - t411 + Ifges(6,2) * t548 + t597;
t589 = -t330 * Ifges(4,6) / 0.2e1;
t588 = Ifges(5,3) * t594;
t569 = t345 * t130 + t350 * t135;
t67 = pkin(13) * t249 + t569;
t147 = -t346 * t246 + t351 * (t256 * t341 + t281 * t338);
t134 = -pkin(4) * t303 - t147;
t81 = pkin(5) * t212 - pkin(13) * t213 + t134;
t29 = t344 * t81 + t349 * t67;
t587 = -qJD(6) * t29 + t344 * t602 + t349 * t599;
t28 = -t344 * t67 + t349 * t81;
t586 = qJD(6) * t28 + t344 * t599 - t349 * t602;
t519 = t43 * mrSges(6,3);
t290 = t331 + (-pkin(3) * t351 - pkin(4)) * t341;
t214 = pkin(5) * t307 - pkin(13) * t308 + t290;
t216 = -pkin(13) * t478 + t567;
t138 = t214 * t349 - t216 * t344;
t581 = qJD(6) * t138 + t344 * t600 + t349 * t601;
t139 = t214 * t344 + t216 * t349;
t580 = -qJD(6) * t139 - t344 * t601 + t349 * t600;
t578 = pkin(5) * t607 - t576;
t574 = (t258 * t341 + t286 * t338) * t351 + t611;
t573 = -t176 + t92;
t571 = -(t209 * t341 + t285 * t338) * t351 + t610;
t568 = t210 - t228;
t8 = -qJD(5) * t43 - t345 * t38 + t350 * t72;
t566 = -t345 * t8 + t350 * t7;
t565 = t1 * t349 - t2 * t344;
t232 = Ifges(5,4) * t235;
t487 = t277 * Ifges(5,5);
t489 = t236 * Ifges(5,1);
t146 = t232 + t487 + t489;
t540 = t146 / 0.2e1;
t564 = t136 * mrSges(5,2) + t540 + t487 / 0.2e1 - t96 * mrSges(5,3) + t232 / 0.2e1;
t563 = Ifges(5,1) * t596 + Ifges(5,4) * t595;
t21 = -qJD(5) * t569 + t104 * t350 - t345 * t87;
t562 = -t84 * mrSges(6,1) - t14 * mrSges(7,1) + t15 * mrSges(7,2);
t493 = t234 * Ifges(6,6);
t503 = t191 * Ifges(6,2);
t513 = Ifges(6,4) * t192;
t101 = t493 + t503 + t513;
t504 = t186 * Ifges(7,3);
t505 = t141 * Ifges(7,5);
t506 = t140 * Ifges(7,6);
t59 = t504 + t505 + t506;
t561 = t101 / 0.2e1 - t59 / 0.2e1 + t562;
t553 = t51 / 0.2e1;
t554 = t50 / 0.2e1;
t560 = -mrSges(6,3) * t7 - Ifges(6,4) * t550 + Ifges(7,5) * t554 - Ifges(6,2) * t549 - Ifges(6,6) * t534 + Ifges(7,6) * t553 + Ifges(7,3) * t548 + t590;
t12 = t50 * Ifges(7,4) + t51 * Ifges(7,2) + t112 * Ifges(7,6);
t558 = t12 / 0.2e1;
t13 = t50 * Ifges(7,1) + t51 * Ifges(7,4) + t112 * Ifges(7,5);
t557 = t13 / 0.2e1;
t552 = -t60 / 0.2e1;
t492 = t234 * Ifges(6,3);
t500 = t192 * Ifges(6,5);
t502 = t191 * Ifges(6,6);
t100 = t492 + t500 + t502;
t551 = -t100 / 0.2e1;
t188 = Ifges(5,6) * t194;
t189 = Ifges(5,5) * t193;
t115 = Ifges(5,3) * t385 - t188 + t189;
t547 = t115 / 0.2e1;
t546 = Ifges(5,5) * t385 / 0.2e1 + t563;
t545 = -t140 / 0.2e1;
t543 = -t141 / 0.2e1;
t486 = t277 * Ifges(5,6);
t491 = t235 * Ifges(5,2);
t514 = Ifges(5,4) * t236;
t145 = t486 + t491 + t514;
t541 = t145 / 0.2e1;
t539 = -t185 / 0.2e1;
t538 = -t186 / 0.2e1;
t536 = t191 / 0.2e1;
t535 = t192 / 0.2e1;
t533 = t234 / 0.2e1;
t532 = -t277 / 0.2e1;
t523 = t40 * mrSges(5,2);
t41 = (t161 * t341 + t244 * t338) * t351 + t421;
t522 = t41 * mrSges(5,1);
t521 = t42 * mrSges(6,1);
t520 = t43 * mrSges(6,2);
t518 = qJD(3) / 0.2e1;
t19 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t82 = mrSges(6,1) * t194 - mrSges(6,3) * t111;
t517 = t19 - t82;
t516 = mrSges(5,3) * t236;
t515 = Ifges(4,4) * t347;
t110 = Ifges(6,5) * t111;
t109 = Ifges(6,6) * t112;
t498 = t193 * Ifges(5,4);
t143 = mrSges(6,1) * t234 - mrSges(6,3) * t192;
t80 = -mrSges(7,1) * t140 + mrSges(7,2) * t141;
t482 = t143 - t80;
t373 = t342 * t465 + t464;
t477 = t339 * t347;
t252 = t340 * t373 + t343 * t477;
t481 = t252 * t346;
t475 = t340 * t348;
t470 = t344 * t350;
t463 = t349 * t350;
t261 = -t308 * t344 - t349 * t478;
t180 = qJD(6) * t261 + t259 * t349 + t344 * t427;
t183 = t231 * t349 + t271 * t344;
t460 = t180 - t183;
t381 = -t308 * t349 + t344 * t478;
t181 = qJD(6) * t381 - t259 * t344 + t349 * t427;
t182 = -t231 * t344 + t271 * t349;
t459 = t181 - t182;
t121 = -mrSges(6,1) * t191 + mrSges(6,2) * t192;
t197 = mrSges(5,1) * t277 - t516;
t458 = -t197 + t121;
t438 = t343 * t476;
t44 = Ifges(6,3) * t194 - t109 + t110;
t416 = t452 * t475;
t414 = t343 * t429;
t412 = -t8 * mrSges(6,1) + t7 * mrSges(6,2);
t409 = t522 - t523;
t172 = (-qJD(2) * t418 - t450) * t347 + t436;
t173 = t380 * qJD(3) + (t365 - t414) * qJD(1);
t407 = t173 * mrSges(4,1) - t172 * mrSges(4,2);
t406 = -mrSges(4,1) * t352 + mrSges(4,2) * t347;
t405 = mrSges(7,1) * t349 - mrSges(7,2) * t344;
t402 = Ifges(7,1) * t344 + t510;
t400 = Ifges(7,2) * t349 + t511;
t398 = Ifges(7,5) * t344 + Ifges(7,6) * t349;
t395 = t345 * t43 + t350 * t42;
t68 = t170 * t350 - t345 * t96;
t392 = t338 * t416;
t375 = t342 * t461 - t467;
t250 = t340 * t375 + t438;
t304 = -t339 * t340 * t353 + t342 * t343;
t389 = t250 * t341 + t304 * t338;
t153 = t252 * t351 + t346 * t389;
t211 = -t250 * t338 + t304 * t341;
t114 = t153 * t350 + t211 * t345;
t152 = -t250 * t473 - t304 * t478 + t481;
t74 = t114 * t349 + t152 * t344;
t73 = -t114 * t344 + t152 * t349;
t75 = t130 * t350 - t135 * t345;
t113 = t153 * t345 - t350 * t211;
t157 = t213 * t349 + t249 * t344;
t156 = -t213 * t344 + t249 * t349;
t224 = -t291 * t345 + t292 * t350;
t276 = -t322 * t339 + t329;
t370 = t276 * mrSges(4,1) + t589 - (Ifges(4,2) * t352 + t515) * t452 / 0.2e1 - t223 * mrSges(4,3);
t222 = -t314 * t347 + t455;
t327 = Ifges(4,4) * t431;
t363 = Ifges(4,1) * t432 / 0.2e1 + t327 / 0.2e1 + t330 * Ifges(4,5) - t222 * mrSges(4,3) + t276 * mrSges(4,2);
t362 = t96 * mrSges(5,1) - t97 * mrSges(5,2) + t236 * Ifges(5,5) + t235 * Ifges(5,6) + t588;
t359 = t513 / 0.2e1 - t506 / 0.2e1 - t505 / 0.2e1 - t504 / 0.2e1 + t493 / 0.2e1 + t561;
t358 = t520 - t521 - t136 * mrSges(5,1) + t486 / 0.2e1 - t502 / 0.2e1 - t500 / 0.2e1 - t492 / 0.2e1 + t551 + t541 + t514 / 0.2e1;
t357 = t503 / 0.2e1 + t359;
t324 = Ifges(4,5) * t352 * t422;
t310 = -pkin(10) * t477 + t336;
t301 = t312 * qJD(3);
t300 = -pkin(10) * t429 + t328;
t296 = t406 * t452;
t295 = -mrSges(4,2) * t330 + mrSges(4,3) * t431;
t294 = mrSges(4,1) * t330 - mrSges(4,3) * t432;
t288 = pkin(12) * t463 + t325 * t344;
t287 = -pkin(12) * t470 + t325 * t349;
t284 = (mrSges(4,1) * t347 + mrSges(4,2) * t352) * t422;
t215 = pkin(5) * t478 - t224;
t207 = -t414 + (-qJD(3) * t373 + t366) * t340;
t206 = qJD(3) * t438 + (qJD(2) * t372 + qJD(3) * t375) * t340;
t196 = -mrSges(5,2) * t277 + mrSges(5,3) * t235;
t174 = -t207 * t338 + t341 * t416;
t169 = -mrSges(5,1) * t235 + mrSges(5,2) * t236;
t166 = -mrSges(5,2) * t385 - mrSges(5,3) * t194;
t165 = mrSges(5,1) * t385 - mrSges(5,3) * t193;
t164 = t235 * t463 + t236 * t344;
t163 = -t235 * t470 + t236 * t349;
t142 = -mrSges(6,2) * t234 + mrSges(6,3) * t191;
t131 = t176 * t345 - t350 * t228;
t123 = mrSges(5,1) * t194 + mrSges(5,2) * t193;
t122 = pkin(5) * t192 - pkin(13) * t191;
t116 = -t194 * Ifges(5,2) + Ifges(5,6) * t385 + t498;
t99 = mrSges(7,1) * t186 - mrSges(7,3) * t141;
t98 = -mrSges(7,2) * t186 + mrSges(7,3) * t140;
t83 = -mrSges(6,2) * t194 - mrSges(6,3) * t112;
t79 = qJD(4) * t153 + t206 * t346 - t207 * t473 - t351 * t392;
t78 = t206 * t351 + (t207 * t341 + t392) * t346 + (t351 * t389 - t481) * qJD(4);
t66 = -pkin(5) * t249 - t75;
t65 = -qJD(6) * t157 - t126 * t344 + t205 * t349;
t64 = qJD(6) * t156 + t126 * t349 + t205 * t344;
t56 = mrSges(6,1) * t112 + mrSges(6,2) * t111;
t54 = -pkin(5) * t236 - t68;
t33 = qJD(5) * t114 - t350 * t174 + t345 * t78;
t32 = -qJD(5) * t113 + t174 * t345 + t350 * t78;
t31 = -mrSges(7,2) * t112 + mrSges(7,3) * t51;
t30 = mrSges(7,1) * t112 - mrSges(7,3) * t50;
t27 = t122 * t344 + t349 * t42;
t26 = t122 * t349 - t344 * t42;
t18 = -pkin(5) * t205 - t21;
t10 = -qJD(6) * t74 - t32 * t344 + t349 * t79;
t9 = qJD(6) * t73 + t32 * t349 + t344 * t79;
t6 = -pkin(5) * t194 - t8;
t3 = [t10 * t99 + t114 * t83 + t211 * t123 + t32 * t142 + t153 * t166 + t174 * t169 + t78 * t196 + t206 * t295 + t207 * t294 + t304 * t284 + t73 * t30 + t74 * t31 + t9 * t98 + t458 * t79 + (-mrSges(3,1) * t348 - mrSges(3,2) * t353) * qJD(2) ^ 2 * t340 - t482 * t33 + (-t165 + t56) * t152 + t517 * t113 + (t296 * t475 + (-t250 * t352 - t252 * t347) * qJD(3) * mrSges(4,3)) * t452 + m(4) * (t172 * t252 + t173 * t250 + t206 * t223 + t207 * t222 + (qJD(1) * t304 + t276) * t416) + m(5) * (t125 * t211 + t136 * t174 - t152 * t41 + t153 * t40 + t78 * t97 - t79 * t96) + m(6) * (-t113 * t8 + t114 * t7 + t152 * t39 + t32 * t43 - t33 * t42 + t79 * t84) + m(7) * (t1 * t74 + t10 * t14 + t113 * t6 + t15 * t9 + t2 * t73 + t33 * t35); (t126 * t84 - t205 * t43 + t213 * t39 - t249 * t7) * mrSges(6,2) + (Ifges(7,5) * t64 + Ifges(7,6) * t65) * t537 + (Ifges(7,5) * t157 + Ifges(7,6) * t156) * t548 + (Ifges(7,4) * t64 + Ifges(7,2) * t65) * t544 + (Ifges(7,4) * t157 + Ifges(7,2) * t156) * t553 + t569 * t83 + (t134 * t39 + t7 * t569 + t75 * t8 + t575 * t84 + t579 * t43 + (t131 + t21) * t42) * m(6) + t482 * t131 + t235 * (Ifges(5,4) * t204 - Ifges(5,2) * t205) / 0.2e1 - t303 * t523 + m(4) * (t172 * t312 + t173 * t310 - t222 * t301 + t223 * t300) + (-t272 - t301) * t294 + ((-m(4) * pkin(2) + t406) * t419 + ((Ifges(4,5) * t342 / 0.2e1 - t310 * mrSges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,4) * t476) * t352 + (-Ifges(4,6) * t342 - t312 * mrSges(4,3) + t338 * (Ifges(5,5) * t251 - Ifges(5,6) * t249 + Ifges(5,3) * t303) / 0.2e1 + (-0.3e1 / 0.2e1 * t515 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t352) * t339) * t347) * qJD(3)) * t452 + (Ifges(6,4) * t126 + Ifges(6,6) * t205) * t536 + (Ifges(6,4) * t213 + Ifges(6,6) * t249) * t549 + t202 * t123 + (t1 * t156 - t14 * t64 + t15 * t65 - t157 * t2) * mrSges(7,3) + (Ifges(7,1) * t64 + Ifges(7,4) * t65) * t542 + (Ifges(7,1) * t157 + Ifges(7,4) * t156) * t554 + (-t274 + t300) * t295 + (Ifges(6,1) * t126 + Ifges(6,5) * t205) * t535 + (Ifges(6,1) * t213 + Ifges(6,5) * t249) * t550 + (-t204 * t96 - t205 * t97 - t249 * t40 - t251 * t41) * mrSges(5,3) + (Ifges(6,5) * t126 + Ifges(6,3) * t205) * t533 + (Ifges(6,5) * t213 + Ifges(6,3) * t249) * t534 + t147 * t165 + t148 * t166 + t6 * (-mrSges(7,1) * t156 + mrSges(7,2) * t157) - m(4) * (t222 * t272 + t223 * t274 + t276 * t419) + t21 * t143 + t134 * t56 + t126 * t102 / 0.2e1 + t75 * t82 + t18 * t80 + t35 * (-mrSges(7,1) * t65 + mrSges(7,2) * t64) + t65 * t60 / 0.2e1 + t66 * t19 + t64 * t61 / 0.2e1 + t29 * t31 + t28 * t30 + (t324 / 0.2e1 + t407) * t342 + t586 * t98 + t587 * t99 + (t1 * t29 + t2 * t28 + t6 * t66 + (-t131 + t18) * t35 + t586 * t15 + t587 * t14) * m(7) + (-t296 * t434 - pkin(2) * t284 + (t172 * t352 - t173 * t347) * mrSges(4,3) + (t363 * t352 + (t589 + (t588 + t362) * t338 + t370) * t347) * qJD(3)) * t339 + (-Ifges(6,4) * t535 + Ifges(7,5) * t542 - Ifges(6,2) * t536 - Ifges(6,6) * t533 + Ifges(7,6) * t544 + Ifges(7,3) * t537 - t519 - t561) * t127 + t573 * t196 + t574 * t197 + (t125 * t202 + t136 * t568 + t147 * t41 + t148 * t40 + t573 * t97 + t574 * t96) * m(5) + t575 * t121 + t579 * t142 + t568 * t169 + t236 * (Ifges(5,1) * t204 - Ifges(5,4) * t205) / 0.2e1 + t42 * (mrSges(6,1) * t205 - mrSges(6,3) * t126) + t205 * t100 / 0.2e1 - t205 * t145 / 0.2e1 + t136 * (mrSges(5,1) * t205 + mrSges(5,2) * t204) + t249 * t44 / 0.2e1 + t8 * (mrSges(6,1) * t249 - mrSges(6,3) * t213) - t249 * t116 / 0.2e1 + t125 * (mrSges(5,1) * t249 + mrSges(5,2) * t251) + t303 * t522 + t204 * t540 + t251 * t546 + t303 * t547 + t213 * t555 + t157 * t557 + t156 * t558 + (Ifges(5,5) * t204 - Ifges(5,6) * t205) * t594 + (Ifges(5,4) * t251 - Ifges(5,2) * t249 + Ifges(5,6) * t303) * t595 + (Ifges(5,1) * t251 - Ifges(5,4) * t249 + Ifges(5,5) * t303) * t596 + t560 * t212; (t101 - t59) * (t230 / 0.2e1 - t260 / 0.2e1) + (-mrSges(7,1) * t459 + mrSges(7,2) * t460) * t35 + (t224 * t8 + t290 * t39 + t42 * t576 + t43 * t577 + t567 * t7 + t572 * t84) * m(6) + t567 * t83 + (-Ifges(7,1) * t381 + Ifges(7,4) * t261) * t554 + (-Ifges(7,4) * t381 + Ifges(7,2) * t261) * t553 + (t1 * t261 - t14 * t460 + t15 * t459 + t2 * t381) * mrSges(7,3) + t6 * (-mrSges(7,1) * t261 - mrSges(7,2) * t381) + (-Ifges(7,5) * t381 + Ifges(7,6) * t261) * t548 - t381 * t557 + (t547 + t189 / 0.2e1 - t188 / 0.2e1 + t409) * t341 + (t181 / 0.2e1 - t182 / 0.2e1) * t60 + (-t231 / 0.2e1 + t259 / 0.2e1) * t102 - t271 * t521 + (-t582 + t583) * t608 + (t562 + t519) * t609 + ((-t327 / 0.2e1 - t363) * t352 + ((t515 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t352) * t452 + (-qJD(3) + t330 / 0.2e1) * Ifges(4,6) + ((Ifges(5,5) * t346 + Ifges(5,6) * t351) * t338 * t518 + (t341 * t518 + t532) * Ifges(5,3) - t362) * t338 - t370) * t347) * t452 - t236 * (Ifges(5,1) * t273 - Ifges(5,4) * t271) / 0.2e1 + (-mrSges(6,3) * t8 + t606) * t308 + t324 + t407 - t235 * (Ifges(5,4) * t273 - Ifges(5,2) * t271) / 0.2e1 - t162 * t169 + t139 * t31 + t138 * t30 + t580 * t99 + t581 * t98 + (t1 * t139 + t138 * t2 + t14 * t580 + t15 * t581 + t215 * t6 + t35 * t578) * m(7) + t576 * t143 + t577 * t142 + t578 * t80 + t570 * t196 + t571 * t197 + (-pkin(3) * t125 * t338 - t136 * t162 + t309 * t41 + t311 * t40 + t570 * t97 + t571 * t96) * m(5) + t572 * t121 + ((t489 / 0.2e1 + t564) * t351 + (-t97 * mrSges(5,3) - t491 / 0.2e1 - t358) * t346) * t448 + (-pkin(3) * t123 + (t125 * mrSges(5,2) - t41 * mrSges(5,3) + t546 + t563) * t346 + (t116 / 0.2e1 - t44 / 0.2e1 - t110 / 0.2e1 + t109 / 0.2e1 + t40 * mrSges(5,3) - t125 * mrSges(5,1) + t498 / 0.2e1 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1) * t194 + t412) * t351) * t338 + t215 * t19 + t224 * t82 + (-t183 / 0.2e1 + t180 / 0.2e1) * t61 - t191 * (Ifges(6,4) * t231 - Ifges(6,2) * t230 + Ifges(6,6) * t271) / 0.2e1 - t192 * (Ifges(6,1) * t231 - Ifges(6,4) * t230 + Ifges(6,5) * t271) / 0.2e1 - t234 * (Ifges(6,5) * t231 - Ifges(6,6) * t230 + Ifges(6,3) * t271) / 0.2e1 - t273 * t146 / 0.2e1 - t136 * (mrSges(5,1) * t271 + mrSges(5,2) * t273) + t290 * t56 + t223 * t294 - t222 * t295 + t309 * t165 + t311 * t166 + t271 * t520 + (Ifges(5,5) * t273 - Ifges(5,6) * t271) * t532 + (Ifges(6,5) * t259 - Ifges(6,6) * t260) * t533 + (Ifges(6,1) * t259 - Ifges(6,4) * t260) * t535 + (Ifges(6,4) * t259 - Ifges(6,2) * t260) * t536 + (Ifges(7,5) * t180 + Ifges(7,6) * t181 + Ifges(7,3) * t260) * t537 + (Ifges(7,5) * t183 + Ifges(7,6) * t182 + Ifges(7,3) * t230) * t538 + t271 * t541 + (Ifges(7,1) * t180 + Ifges(7,4) * t181 + Ifges(7,5) * t260) * t542 + (Ifges(7,1) * t183 + Ifges(7,4) * t182 + Ifges(7,5) * t230) * t543 + (Ifges(7,4) * t180 + Ifges(7,2) * t181 + Ifges(7,6) * t260) * t544 + (Ifges(7,4) * t183 + Ifges(7,2) * t182 + Ifges(7,6) * t230) * t545 + t271 * t551 + t261 * t558 + t560 * t307 + (t271 * t97 + t273 * t96) * mrSges(5,3); t592 * t98 + (t1 * t288 + t14 * t593 + t15 * t592 + t2 * t287 - t35 * t54) * m(7) + t593 * t99 + (t83 * pkin(12) + (t185 / 0.2e1 + t501 / 0.2e1 + (m(7) * t35 - t482) * pkin(12) - t591) * qJD(5) - t49 / 0.2e1 - t48 / 0.2e1 + (-Ifges(6,2) / 0.2e1 - Ifges(7,3) / 0.2e1) * t112 + (t539 - t501 / 0.2e1 + t369) * t235 - t590 - t597) * t350 + ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t236 - t564) * t235 + (t14 * t164 - t15 * t163) * mrSges(7,3) + t409 + (-pkin(4) * t39 - t42 * t68 - t43 * t69 - t84 * t97 + (-qJD(5) * t395 + t566) * pkin(12)) * m(6) + ((-t357 - t519) * qJD(5) + (-t1 * t344 - t2 * t349) * mrSges(7,3) + (t349 * t552 + t61 * t531 + t398 * t538 + t400 * t545 + t402 * t543 + t35 * t405 + (t14 * t344 - t15 * t349) * mrSges(7,3)) * qJD(6) + t12 * t531 + t13 * t530 + t399 * t548 + t401 * t553 + t403 * t554 + t6 * t404 + t357 * t235 + (m(7) * t6 - qJD(5) * t142 + t517) * pkin(12) + t606) * t345 - t96 * t196 + t115 + (-t458 + t516) * t97 - t35 * (-mrSges(7,1) * t163 + mrSges(7,2) * t164) - t164 * t61 / 0.2e1 - t69 * t142 - t68 * t143 - t54 * t80 - pkin(4) * t56 + (t235 * t395 + t566) * mrSges(6,3) + t287 * t30 + t288 * t31 + (Ifges(7,5) * t164 + Ifges(7,6) * t163) * t538 + (Ifges(7,1) * t164 + Ifges(7,4) * t163) * t543 + (Ifges(7,4) * t164 + Ifges(7,2) * t163) * t545 + t163 * t552 + t358 * t236; t482 * t43 + t44 - t412 + t565 * mrSges(7,3) - t42 * t142 - t27 * t98 - t26 * t99 - pkin(5) * t19 - t6 * t405 + t402 * t554 + (t539 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t192 + t591) * t191 + t598 * qJD(6) + t398 * t548 + t400 * t553 + t344 * t557 + t12 * t530 + (t359 + t519) * t192 + (-pkin(5) * t6 - t14 * t26 - t15 * t27 - t35 * t43) * m(7) + (-t30 * t344 + t31 * t349 + m(7) * t565 + (-m(7) * t397 - t344 * t98 - t349 * t99) * qJD(6)) * pkin(13); -t35 * (mrSges(7,1) * t141 + mrSges(7,2) * t140) + t60 * t542 + (Ifges(7,5) * t140 - Ifges(7,6) * t141) * t538 + (Ifges(7,1) * t140 - t512) * t543 - t14 * t98 + t15 * t99 + (t14 * t140 + t141 * t15) * mrSges(7,3) - t411 + t11 + (-Ifges(7,2) * t141 + t137 + t61) * t545;];
tauc  = t3(:);
