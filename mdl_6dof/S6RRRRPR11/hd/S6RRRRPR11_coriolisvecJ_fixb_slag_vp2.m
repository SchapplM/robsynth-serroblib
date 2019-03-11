% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:12:23
% EndTime: 2019-03-09 23:13:25
% DurationCPUTime: 33.02s
% Computational Cost: add. (32639->1026), mult. (83676->1435), div. (0->0), fcn. (66958->12), ass. (0->429)
t411 = sin(qJ(2));
t405 = sin(pkin(6));
t459 = qJD(1) * t405;
t448 = t411 * t459;
t407 = cos(pkin(6));
t415 = cos(qJ(2));
t517 = pkin(1) * t415;
t451 = t407 * t517;
t352 = -pkin(8) * t448 + qJD(1) * t451;
t427 = (pkin(2) * t411 - pkin(9) * t415) * t405;
t353 = qJD(1) * t427;
t410 = sin(qJ(3));
t414 = cos(qJ(3));
t263 = t414 * t352 + t410 * t353;
t247 = pkin(10) * t448 + t263;
t396 = t407 * t411 * pkin(1);
t438 = pkin(3) * t410 - pkin(10) * t414;
t478 = t405 * t415;
t272 = (t396 + (pkin(8) + t438) * t478) * qJD(1);
t409 = sin(qJ(4));
t413 = cos(qJ(4));
t176 = -t409 * t247 + t413 * t272;
t473 = t414 * t415;
t314 = (t409 * t411 + t413 * t473) * t459;
t383 = -pkin(3) * t414 - pkin(10) * t410 - pkin(2);
t474 = t413 * t414;
t398 = pkin(9) * t474;
t447 = t415 * t459;
t440 = t410 * t447;
t452 = qJD(5) * t413;
t378 = t438 * qJD(3);
t456 = qJD(3) * t410;
t512 = pkin(9) * t409;
t461 = t413 * t378 + t456 * t512;
t513 = pkin(4) * t410;
t601 = -pkin(4) * t440 + t314 * qJ(5) - t176 - t410 * t452 + (-qJ(5) * t474 + t513) * qJD(3) + (-t398 + (qJ(5) * t410 - t383) * t409) * qJD(4) + t461;
t177 = t413 * t247 + t409 * t272;
t313 = (-t409 * t473 + t411 * t413) * t459;
t453 = qJD(4) * t413;
t462 = t409 * t378 + t383 * t453;
t475 = t410 * t413;
t600 = qJ(5) * t313 + t177 - (-pkin(9) * qJD(3) - qJ(5) * qJD(4)) * t475 - (-qJD(5) * t410 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t414) * t409 - t462;
t404 = sin(pkin(12));
t406 = cos(pkin(12));
t240 = t313 * t406 - t314 * t404;
t428 = t404 * t409 - t406 * t413;
t363 = t428 * qJD(4);
t372 = t404 * t413 + t406 * t409;
t455 = qJD(3) * t414;
t265 = t363 * t410 - t372 * t455;
t466 = t240 - t265;
t241 = t313 * t404 + t314 * t406;
t362 = t372 * qJD(4);
t266 = -t362 * t410 - t428 * t455;
t465 = t241 - t266;
t579 = t600 * t404 + t406 * t601;
t578 = t404 * t601 - t600 * t406;
t392 = qJD(1) * t407 + qJD(2);
t326 = t392 * t414 - t410 * t448;
t238 = t372 * t326;
t599 = t238 - t362;
t239 = t428 * t326;
t467 = -t239 + t363;
t460 = pkin(8) * t478 + t396;
t355 = t460 * qJD(1);
t307 = t392 * pkin(9) + t355;
t347 = (-pkin(2) * t415 - pkin(9) * t411 - pkin(1)) * t405;
t319 = qJD(1) * t347;
t235 = -t410 * t307 + t319 * t414;
t384 = qJD(3) - t447;
t217 = -pkin(3) * t384 - t235;
t327 = t392 * t410 + t414 * t448;
t273 = -t327 * t409 + t384 * t413;
t167 = -pkin(4) * t273 + qJD(5) + t217;
t274 = t327 * t413 + t384 * t409;
t441 = t406 * t273 - t274 * t404;
t101 = -pkin(5) * t441 + t167;
t196 = t273 * t404 + t274 * t406;
t408 = sin(qJ(6));
t412 = cos(qJ(6));
t574 = -t196 * t408 + t412 * t441;
t116 = Ifges(7,4) * t574;
t120 = t196 * t412 + t408 * t441;
t321 = qJD(4) - t326;
t306 = -t392 * pkin(2) - t352;
t215 = -t326 * pkin(3) - t327 * pkin(10) + t306;
t236 = t414 * t307 + t410 * t319;
t218 = pkin(10) * t384 + t236;
t136 = t215 * t409 + t218 * t413;
t112 = qJ(5) * t273 + t136;
t105 = t404 * t112;
t135 = t413 * t215 - t218 * t409;
t111 = -qJ(5) * t274 + t135;
t92 = pkin(4) * t321 + t111;
t55 = t406 * t92 - t105;
t584 = pkin(11) * t196;
t40 = pkin(5) * t321 + t55 - t584;
t477 = t406 * t112;
t56 = t404 * t92 + t477;
t571 = pkin(11) * t441;
t43 = t56 + t571;
t15 = t40 * t412 - t408 * t43;
t16 = t40 * t408 + t412 * t43;
t506 = Ifges(7,4) * t120;
t309 = qJD(6) + t321;
t528 = -t309 / 0.2e1;
t543 = -t120 / 0.2e1;
t545 = -t574 / 0.2e1;
t66 = Ifges(7,1) * t120 + Ifges(7,5) * t309 + t116;
t598 = (Ifges(7,5) * t574 - Ifges(7,6) * t120) * t528 + (t120 * t16 + t15 * t574) * mrSges(7,3) + (-Ifges(7,2) * t120 + t116 + t66) * t545 - t101 * (mrSges(7,1) * t120 + mrSges(7,2) * t574) + (Ifges(7,1) * t574 - t506) * t543;
t597 = t579 + t465 * pkin(11) + (-t440 + t456) * pkin(5);
t596 = pkin(11) * t466 - t578;
t510 = -qJ(5) - pkin(10);
t442 = qJD(4) * t510;
t360 = t409 * t442 + t452;
t361 = -qJD(5) * t409 + t413 * t442;
t267 = -t360 * t404 + t406 * t361;
t255 = pkin(3) * t327 - pkin(10) * t326;
t165 = -t235 * t409 + t413 * t255;
t128 = -qJ(5) * t326 * t413 + pkin(4) * t327 + t165;
t166 = t413 * t235 + t409 * t255;
t480 = t326 * t409;
t144 = -qJ(5) * t480 + t166;
t77 = t406 * t128 - t144 * t404;
t595 = t267 - t77;
t268 = t406 * t360 + t404 * t361;
t78 = t404 * t128 + t406 * t144;
t594 = t268 - t78;
t65 = Ifges(7,2) * t574 + Ifges(7,6) * t309 + t506;
t592 = t65 / 0.2e1;
t587 = -pkin(5) * t327 + pkin(11) * t467 + t595;
t586 = pkin(11) * t599 + t594;
t262 = -t410 * t352 + t353 * t414;
t246 = -pkin(3) * t448 - t262;
t514 = pkin(4) * t409;
t576 = pkin(4) * t313 - t246 + t453 * t513 + (pkin(9) + t514) * t455;
t497 = t274 * Ifges(5,4);
t171 = t273 * Ifges(5,2) + t321 * Ifges(5,6) + t497;
t271 = Ifges(5,4) * t273;
t172 = t274 * Ifges(5,1) + t321 * Ifges(5,5) + t271;
t429 = t135 * t413 + t136 * t409;
t431 = Ifges(5,5) * t413 - Ifges(5,6) * t409;
t507 = Ifges(5,4) * t413;
t433 = -Ifges(5,2) * t409 + t507;
t508 = Ifges(5,4) * t409;
t435 = Ifges(5,1) * t413 - t508;
t436 = mrSges(5,1) * t409 + mrSges(5,2) * t413;
t518 = t413 / 0.2e1;
t519 = -t409 / 0.2e1;
t524 = t321 / 0.2e1;
t531 = t274 / 0.2e1;
t533 = t273 / 0.2e1;
t585 = t429 * mrSges(5,3) - t171 * t519 - t172 * t518 - t217 * t436 - t431 * t524 - t433 * t533 - t435 * t531;
t457 = qJD(2) * t415;
t445 = t414 * t457;
t280 = t392 * t455 + (-t411 * t456 + t445) * t459;
t458 = qJD(2) * t405;
t443 = qJD(1) * t458;
t439 = t411 * t443;
t183 = qJD(4) * t273 + t280 * t413 + t409 * t439;
t184 = -qJD(4) * t274 - t280 * t409 + t413 * t439;
t110 = t183 * t406 + t184 * t404;
t446 = t410 * t457;
t281 = t392 * t456 + (t411 * t455 + t446) * t459;
t354 = qJD(2) * t427;
t337 = qJD(1) * t354;
t479 = t405 * t411;
t393 = pkin(8) * t479;
t368 = -t393 + t451;
t356 = t368 * qJD(2);
t338 = qJD(1) * t356;
t168 = -t307 * t456 + t319 * t455 + t410 * t337 + t414 * t338;
t159 = pkin(10) * t439 + t168;
t357 = t460 * qJD(2);
t339 = qJD(1) * t357;
t180 = t281 * pkin(3) - t280 * pkin(10) + t339;
t69 = -qJD(4) * t136 - t159 * t409 + t413 * t180;
t39 = pkin(4) * t281 - qJ(5) * t183 - qJD(5) * t274 + t69;
t454 = qJD(4) * t409;
t68 = t413 * t159 + t409 * t180 + t215 * t453 - t218 * t454;
t42 = qJ(5) * t184 + qJD(5) * t273 + t68;
t12 = t406 * t39 - t404 * t42;
t6 = pkin(5) * t281 - pkin(11) * t110 + t12;
t109 = -t183 * t404 + t184 * t406;
t13 = t404 * t39 + t406 * t42;
t7 = pkin(11) * t109 + t13;
t2 = qJD(6) * t15 + t408 * t6 + t412 * t7;
t3 = -qJD(6) * t16 - t408 * t7 + t412 * t6;
t583 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t374 = t413 * t383;
t288 = -qJ(5) * t475 + t374 + (-pkin(4) - t512) * t414;
t335 = t409 * t383 + t398;
t476 = t409 * t410;
t300 = -qJ(5) * t476 + t335;
t219 = t406 * t288 - t300 * t404;
t349 = t428 * t410;
t185 = -pkin(5) * t414 + pkin(11) * t349 + t219;
t220 = t404 * t288 + t406 * t300;
t348 = t372 * t410;
t190 = -pkin(11) * t348 + t220;
t114 = t185 * t412 - t190 * t408;
t582 = qJD(6) * t114 + t408 * t597 - t596 * t412;
t115 = t185 * t408 + t190 * t412;
t581 = -qJD(6) * t115 + t596 * t408 + t412 * t597;
t580 = Ifges(6,4) * t196;
t577 = pkin(5) * t466 + t576;
t320 = Ifges(4,4) * t326;
t486 = t384 * Ifges(4,5);
t489 = t327 * Ifges(4,1);
t232 = t320 + t486 + t489;
t422 = t235 * mrSges(4,3) - t232 / 0.2e1 - t306 * mrSges(4,2) - t486 / 0.2e1;
t575 = t422 + t585;
t35 = qJD(6) * t574 + t109 * t408 + t110 * t412;
t555 = t35 / 0.2e1;
t36 = -qJD(6) * t120 + t109 * t412 - t110 * t408;
t554 = t36 / 0.2e1;
t99 = Ifges(6,2) * t441 + Ifges(6,6) * t321 + t580;
t573 = t99 / 0.2e1;
t547 = t109 / 0.2e1;
t546 = t110 / 0.2e1;
t530 = t281 / 0.2e1;
t538 = -t441 / 0.2e1;
t572 = -Ifges(3,6) * t392 / 0.2e1;
t570 = Ifges(5,3) + Ifges(6,3);
t569 = -t69 * mrSges(5,1) - t12 * mrSges(6,1) + t68 * mrSges(5,2) + t13 * mrSges(6,2) - t583;
t385 = t510 * t409;
t386 = t510 * t413;
t304 = t406 * t385 + t386 * t404;
t269 = -pkin(11) * t372 + t304;
t305 = t404 * t385 - t406 * t386;
t270 = -pkin(11) * t428 + t305;
t188 = t269 * t412 - t270 * t408;
t568 = qJD(6) * t188 + t408 * t587 + t412 * t586;
t189 = t269 * t408 + t270 * t412;
t567 = -qJD(6) * t189 - t408 * t586 + t412 * t587;
t566 = Ifges(6,4) * t441;
t400 = pkin(4) * t406 + pkin(5);
t515 = pkin(4) * t404;
t359 = t400 * t408 + t412 * t515;
t62 = -t111 * t404 - t477;
t44 = t62 - t571;
t63 = t406 * t111 - t105;
t45 = t63 - t584;
t565 = -t359 * qJD(6) + t408 * t45 - t412 * t44;
t358 = t400 * t412 - t408 * t515;
t564 = t358 * qJD(6) - t408 * t44 - t412 * t45;
t345 = t393 + (-pkin(2) - t517) * t407;
t364 = -t407 * t414 + t410 * t479;
t365 = t407 * t410 + t414 * t479;
t243 = t364 * pkin(3) - t365 * pkin(10) + t345;
t346 = pkin(9) * t407 + t460;
t257 = t414 * t346 + t410 * t347;
t245 = -pkin(10) * t478 + t257;
t162 = t409 * t243 + t413 * t245;
t198 = pkin(4) * t480 + t236;
t561 = pkin(4) * t454 - pkin(5) * t599 - t198;
t560 = -t409 * t69 + t413 * t68;
t103 = Ifges(6,6) * t109;
t104 = Ifges(6,5) * t110;
t51 = Ifges(6,3) * t281 + t103 + t104;
t33 = Ifges(7,6) * t36;
t34 = Ifges(7,5) * t35;
t8 = Ifges(7,3) * t281 + t33 + t34;
t181 = Ifges(5,6) * t184;
t182 = Ifges(5,5) * t183;
t85 = Ifges(5,3) * t281 + t181 + t182;
t559 = t85 + t51 + t8;
t496 = t274 * Ifges(5,5);
t498 = t273 * Ifges(5,6);
t170 = t321 * Ifges(5,3) + t496 + t498;
t485 = t384 * Ifges(4,6);
t488 = t327 * Ifges(4,4);
t491 = t326 * Ifges(4,2);
t231 = t485 + t488 + t491;
t450 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t492 = t309 * Ifges(7,3);
t499 = t196 * Ifges(6,5);
t500 = t441 * Ifges(6,6);
t503 = t120 * Ifges(7,5);
t504 = t574 * Ifges(7,6);
t64 = t492 + t503 + t504;
t98 = t321 * Ifges(6,3) + t499 + t500;
t416 = -t450 * t321 + t136 * mrSges(5,2) + t16 * mrSges(7,2) + t236 * mrSges(4,3) + t56 * mrSges(6,2) - t170 / 0.2e1 + t231 / 0.2e1 - t64 / 0.2e1 - t98 / 0.2e1 - t504 / 0.2e1 - t503 / 0.2e1 - t135 * mrSges(5,1) - t15 * mrSges(7,1) - t500 / 0.2e1 - t499 / 0.2e1 - t498 / 0.2e1 - t496 / 0.2e1 - t306 * mrSges(4,1) - t492 / 0.2e1 + t488 / 0.2e1 + t485 / 0.2e1 - t55 * mrSges(6,1);
t558 = t416 + t491 / 0.2e1;
t557 = Ifges(7,4) * t555 + Ifges(7,2) * t554 + Ifges(7,6) * t530;
t556 = Ifges(7,1) * t555 + Ifges(7,4) * t554 + Ifges(7,5) * t530;
t553 = Ifges(6,4) * t546 + Ifges(6,2) * t547 + Ifges(6,6) * t530;
t552 = Ifges(6,1) * t546 + Ifges(6,4) * t547 + Ifges(6,5) * t530;
t87 = t183 * Ifges(5,1) + t184 * Ifges(5,4) + t281 * Ifges(5,5);
t551 = t87 / 0.2e1;
t550 = pkin(1) * mrSges(3,1);
t549 = pkin(1) * mrSges(3,2);
t544 = t574 / 0.2e1;
t542 = t120 / 0.2e1;
t541 = -t171 / 0.2e1;
t540 = t183 / 0.2e1;
t539 = t184 / 0.2e1;
t537 = t441 / 0.2e1;
t536 = -t196 / 0.2e1;
t535 = t196 / 0.2e1;
t534 = -t273 / 0.2e1;
t532 = -t274 / 0.2e1;
t527 = t309 / 0.2e1;
t526 = -t320 / 0.2e1;
t525 = -t321 / 0.2e1;
t523 = -t364 / 0.2e1;
t521 = t365 / 0.2e1;
t520 = t407 / 0.2e1;
t516 = pkin(4) * t274;
t511 = pkin(9) * t414;
t295 = -qJD(3) * t364 + t405 * t445;
t296 = -t365 * t409 - t413 * t478;
t444 = t411 * t458;
t214 = qJD(4) * t296 + t295 * t413 + t409 * t444;
t294 = qJD(3) * t365 + t405 * t446;
t426 = -t365 * t413 + t409 * t478;
t186 = -t346 * t456 + t347 * t455 + t410 * t354 + t414 * t356;
t174 = pkin(10) * t444 + t186;
t204 = t294 * pkin(3) - t295 * pkin(10) + t357;
t80 = -qJD(4) * t162 - t174 * t409 + t413 * t204;
t50 = pkin(4) * t294 - qJ(5) * t214 + qJD(5) * t426 + t80;
t213 = qJD(4) * t426 - t295 * t409 + t413 * t444;
t79 = t413 * t174 + t409 * t204 + t243 * t453 - t245 * t454;
t57 = qJ(5) * t213 + qJD(5) * t296 + t79;
t21 = t404 * t50 + t406 * t57;
t509 = Ifges(3,4) * t411;
t502 = t168 * mrSges(4,2);
t169 = -t307 * t455 - t319 * t456 + t337 * t414 - t410 * t338;
t501 = t169 * mrSges(4,1);
t495 = t280 * Ifges(4,1);
t494 = t280 * Ifges(4,4);
t493 = t281 * Ifges(4,4);
t161 = t413 * t243 - t245 * t409;
t127 = pkin(4) * t364 + qJ(5) * t426 + t161;
t137 = qJ(5) * t296 + t162;
t76 = t404 * t127 + t406 * t137;
t258 = -t348 * t412 + t349 * t408;
t145 = qJD(6) * t258 + t265 * t408 + t266 * t412;
t158 = t240 * t408 + t241 * t412;
t472 = t145 - t158;
t259 = -t348 * t408 - t349 * t412;
t146 = -qJD(6) * t259 + t265 * t412 - t266 * t408;
t157 = t240 * t412 - t241 * t408;
t471 = t146 - t157;
t155 = -t238 * t412 + t239 * t408;
t287 = t372 * t412 - t408 * t428;
t211 = -qJD(6) * t287 - t362 * t412 + t363 * t408;
t470 = t155 - t211;
t156 = -t238 * t408 - t239 * t412;
t286 = -t372 * t408 - t412 * t428;
t210 = qJD(6) * t286 - t362 * t408 - t363 * t412;
t469 = t156 - t210;
t464 = -mrSges(3,1) * t392 - mrSges(4,1) * t326 + mrSges(4,2) * t327 + mrSges(3,3) * t448;
t201 = -mrSges(5,1) * t273 + mrSges(5,2) * t274;
t285 = mrSges(4,1) * t384 - mrSges(4,3) * t327;
t463 = t285 - t201;
t379 = pkin(4) * t476 + t410 * pkin(9);
t449 = Ifges(4,5) * t280 - Ifges(4,6) * t281 + Ifges(4,3) * t439;
t401 = -pkin(4) * t413 - pkin(3);
t11 = -t36 * mrSges(7,1) + t35 * mrSges(7,2);
t20 = -t404 * t57 + t406 * t50;
t60 = -t109 * mrSges(6,1) + t110 * mrSges(6,2);
t75 = t406 * t127 - t137 * t404;
t256 = -t410 * t346 + t347 * t414;
t244 = pkin(3) * t478 - t256;
t437 = mrSges(5,1) * t413 - mrSges(5,2) * t409;
t434 = Ifges(5,1) * t409 + t507;
t432 = Ifges(5,2) * t413 + t508;
t430 = Ifges(5,5) * t409 + Ifges(5,6) * t413;
t223 = t296 * t404 - t406 * t426;
t59 = pkin(5) * t364 - pkin(11) * t223 + t75;
t222 = t296 * t406 + t404 * t426;
t67 = pkin(11) * t222 + t76;
t22 = -t408 * t67 + t412 * t59;
t23 = t408 * t59 + t412 * t67;
t142 = t222 * t412 - t223 * t408;
t143 = t222 * t408 + t223 * t412;
t187 = -t346 * t455 - t347 * t456 + t354 * t414 - t410 * t356;
t197 = -pkin(4) * t296 + t244;
t387 = Ifges(3,4) * t447;
t423 = -t352 * mrSges(3,3) + Ifges(3,1) * t448 / 0.2e1 + t387 / 0.2e1 + t392 * Ifges(3,5);
t175 = -pkin(3) * t444 - t187;
t160 = -pkin(3) * t439 - t169;
t121 = -pkin(4) * t213 + t175;
t93 = -pkin(4) * t184 + t160;
t420 = t235 * mrSges(4,1) + t384 * Ifges(4,3) + t327 * Ifges(4,5) + t326 * Ifges(4,6) + t572 - (Ifges(3,2) * t415 + t509) * t459 / 0.2e1 - t236 * mrSges(4,2) - t355 * mrSges(3,3);
t382 = Ifges(3,5) * t415 * t443;
t351 = -t392 * mrSges(3,2) + mrSges(3,3) * t447;
t336 = pkin(5) * t428 + t401;
t334 = -t409 * t511 + t374;
t289 = pkin(5) * t348 + t379;
t284 = -mrSges(4,2) * t384 + mrSges(4,3) * t326;
t254 = -qJD(4) * t335 + t461;
t253 = (-t413 * t456 - t414 * t454) * pkin(9) + t462;
t251 = -mrSges(4,2) * t439 - mrSges(4,3) * t281;
t250 = mrSges(4,1) * t439 - mrSges(4,3) * t280;
t225 = mrSges(5,1) * t321 - mrSges(5,3) * t274;
t224 = -mrSges(5,2) * t321 + mrSges(5,3) * t273;
t208 = mrSges(4,1) * t281 + mrSges(4,2) * t280;
t192 = Ifges(4,5) * t439 - t493 + t495;
t191 = -t281 * Ifges(4,2) + Ifges(4,6) * t439 + t494;
t164 = mrSges(6,1) * t321 - mrSges(6,3) * t196;
t163 = -mrSges(6,2) * t321 + mrSges(6,3) * t441;
t152 = pkin(5) * t196 + t516;
t150 = -mrSges(5,2) * t281 + mrSges(5,3) * t184;
t149 = mrSges(5,1) * t281 - mrSges(5,3) * t183;
t133 = -pkin(5) * t222 + t197;
t132 = t213 * t404 + t214 * t406;
t130 = t213 * t406 - t214 * t404;
t122 = -mrSges(6,1) * t441 + mrSges(6,2) * t196;
t113 = -mrSges(5,1) * t184 + mrSges(5,2) * t183;
t100 = Ifges(6,1) * t196 + Ifges(6,5) * t321 + t566;
t95 = mrSges(7,1) * t309 - mrSges(7,3) * t120;
t94 = -mrSges(7,2) * t309 + mrSges(7,3) * t574;
t86 = t183 * Ifges(5,4) + t184 * Ifges(5,2) + t281 * Ifges(5,6);
t84 = mrSges(6,1) * t281 - mrSges(6,3) * t110;
t83 = -mrSges(6,2) * t281 + mrSges(6,3) * t109;
t73 = -pkin(5) * t130 + t121;
t71 = -mrSges(7,1) * t574 + mrSges(7,2) * t120;
t58 = -pkin(5) * t109 + t93;
t47 = -qJD(6) * t143 + t130 * t412 - t132 * t408;
t46 = qJD(6) * t142 + t130 * t408 + t132 * t412;
t29 = -mrSges(7,2) * t281 + mrSges(7,3) * t36;
t28 = mrSges(7,1) * t281 - mrSges(7,3) * t35;
t17 = pkin(11) * t130 + t21;
t14 = pkin(5) * t294 - pkin(11) * t132 + t20;
t5 = -qJD(6) * t23 + t14 * t412 - t17 * t408;
t4 = qJD(6) * t22 + t14 * t408 + t17 * t412;
t1 = [t160 * (-mrSges(5,1) * t296 - mrSges(5,2) * t426) + t69 * (mrSges(5,1) * t364 + mrSges(5,3) * t426) + (-Ifges(5,4) * t426 + Ifges(5,2) * t296 + Ifges(5,6) * t364) * t539 + (-Ifges(5,1) * t426 + Ifges(5,4) * t296 + Ifges(5,5) * t364) * t540 - t426 * t551 + (-Ifges(5,5) * t426 + Ifges(6,5) * t223 + Ifges(7,5) * t143 + Ifges(5,6) * t296 + Ifges(6,6) * t222 + Ifges(7,6) * t142 + (Ifges(7,3) + t570) * t364) * t530 + t464 * t357 + t306 * (mrSges(4,1) * t294 + mrSges(4,2) * t295) + t559 * t364 / 0.2e1 + t326 * (Ifges(4,4) * t295 - Ifges(4,2) * t294) / 0.2e1 + t55 * (mrSges(6,1) * t294 - mrSges(6,3) * t132) + t56 * (-mrSges(6,2) * t294 + mrSges(6,3) * t130) + t15 * (mrSges(7,1) * t294 - mrSges(7,3) * t46) + t16 * (-mrSges(7,2) * t294 + mrSges(7,3) * t47) - t294 * t231 / 0.2e1 + t295 * t232 / 0.2e1 + t296 * t86 / 0.2e1 + t136 * (-mrSges(5,2) * t294 + mrSges(5,3) * t213) + t135 * (mrSges(5,1) * t294 - mrSges(5,3) * t214) + t186 * t284 + t187 * t285 - t478 * t501 + (t170 + t98 + t64) * t294 / 0.2e1 + ((Ifges(3,5) * t520 - t368 * mrSges(3,3) + (-0.2e1 * t549 + 0.3e1 / 0.2e1 * Ifges(3,4) * t415) * t405) * t415 + (Ifges(4,5) * t521 + Ifges(4,6) * t523 - Ifges(3,6) * t407 - t460 * mrSges(3,3) + (-0.2e1 * t550 - 0.3e1 / 0.2e1 * t509 + (-Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t415) * t405) * t411) * t443 + m(3) * (t338 * t460 - t339 * t368 - t352 * t357 + t355 * t356) + m(4) * (t168 * t257 + t169 * t256 + t186 * t236 + t187 * t235 + t306 * t357 + t339 * t345) + m(5) * (t135 * t80 + t136 * t79 + t160 * t244 + t161 * t69 + t162 * t68 + t175 * t217) + m(6) * (t12 * t75 + t121 * t167 + t13 * t76 + t197 * t93 + t20 * t55 + t21 * t56) + m(7) * (t101 * t73 + t133 * t58 + t15 * t5 + t16 * t4 + t2 * t23 + t22 * t3) + t384 * (Ifges(4,5) * t295 - Ifges(4,6) * t294) / 0.2e1 - t449 * t478 / 0.2e1 - t281 * (Ifges(4,4) * t365 - Ifges(4,2) * t364 - Ifges(4,6) * t478) / 0.2e1 + t280 * (Ifges(4,1) * t365 - Ifges(4,4) * t364 - Ifges(4,5) * t478) / 0.2e1 + t338 * (-t407 * mrSges(3,2) + mrSges(3,3) * t478) + (-mrSges(3,1) * t407 + mrSges(4,1) * t364 + mrSges(4,2) * t365 + mrSges(3,3) * t479) * t339 + t256 * t250 + t257 * t251 + t244 * t113 + t93 * (-mrSges(6,1) * t222 + mrSges(6,2) * t223) + t79 * t224 + t80 * t225 + t217 * (-mrSges(5,1) * t213 + mrSges(5,2) * t214) + t213 * t171 / 0.2e1 + t214 * t172 / 0.2e1 + t175 * t201 + t197 * t60 + t21 * t163 + t20 * t164 + t167 * (-mrSges(6,1) * t130 + mrSges(6,2) * t132) + t161 * t149 + t162 * t150 + t58 * (-mrSges(7,1) * t142 + mrSges(7,2) * t143) + t133 * t11 + t132 * t100 / 0.2e1 + t121 * t122 + t101 * (-mrSges(7,1) * t47 + mrSges(7,2) * t46) + t5 * t95 + t4 * t94 + t76 * t83 + t75 * t84 + t73 * t71 + t46 * t66 / 0.2e1 + t22 * t28 + t23 * t29 + (-t168 * t364 - t169 * t365 - t235 * t295 - t236 * t294) * mrSges(4,3) + t345 * t208 + t356 * t351 + t68 * (-mrSges(5,2) * t364 + mrSges(5,3) * t296) + t13 * (-mrSges(6,2) * t364 + mrSges(6,3) * t222) + t12 * (mrSges(6,1) * t364 - mrSges(6,3) * t223) + t2 * (-mrSges(7,2) * t364 + mrSges(7,3) * t142) + t3 * (mrSges(7,1) * t364 - mrSges(7,3) * t143) + t478 * t502 + t47 * t592 + t382 * t520 + t192 * t521 + t191 * t523 + (Ifges(7,5) * t46 + Ifges(7,6) * t47 + Ifges(7,3) * t294) * t527 + (Ifges(5,1) * t214 + Ifges(5,4) * t213 + Ifges(5,5) * t294) * t531 + (Ifges(5,4) * t214 + Ifges(5,2) * t213 + Ifges(5,6) * t294) * t533 + (Ifges(6,1) * t132 + Ifges(6,4) * t130 + Ifges(6,5) * t294) * t535 + (Ifges(6,4) * t132 + Ifges(6,2) * t130 + Ifges(6,6) * t294) * t537 + (Ifges(7,1) * t46 + Ifges(7,4) * t47 + Ifges(7,5) * t294) * t542 + (Ifges(7,4) * t46 + Ifges(7,2) * t47 + Ifges(7,6) * t294) * t544 + (Ifges(6,1) * t223 + Ifges(6,4) * t222 + Ifges(6,5) * t364) * t546 + (Ifges(6,4) * t223 + Ifges(6,2) * t222 + Ifges(6,6) * t364) * t547 + t223 * t552 + t222 * t553 + (Ifges(7,4) * t143 + Ifges(7,2) * t142 + Ifges(7,6) * t364) * t554 + (Ifges(7,1) * t143 + Ifges(7,4) * t142 + Ifges(7,5) * t364) * t555 + (Ifges(5,5) * t214 + Ifges(6,5) * t132 + Ifges(5,6) * t213 + Ifges(6,6) * t130 + t294 * t570) * t524 + (t423 * t415 + (t572 + t420) * t411) * t458 + t327 * (Ifges(4,1) * t295 - Ifges(4,4) * t294) / 0.2e1 + t143 * t556 + t142 * t557 + t130 * t573; (-Ifges(6,5) * t349 + Ifges(7,5) * t259 - Ifges(6,6) * t348 + Ifges(7,6) * t258 + t410 * t431) * t530 - t464 * t355 + (mrSges(6,1) * t466 - mrSges(6,2) * t465) * t167 + (t266 / 0.2e1 - t241 / 0.2e1) * t100 + (t253 - t177) * t224 + (-t158 / 0.2e1 + t145 / 0.2e1) * t66 - t217 * (-mrSges(5,1) * t313 + mrSges(5,2) * t314) - t314 * t172 / 0.2e1 + ((qJD(2) * (Ifges(4,5) * t410 + Ifges(4,6) * t414) / 0.2e1 + (t550 + t509 / 0.2e1) * t459 + (t392 / 0.2e1 - qJD(2)) * Ifges(3,6) - t420) * t411 + (-t387 / 0.2e1 + (t549 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t411) * t459 + (t526 - t489 / 0.2e1 + t422) * t414 + t558 * t410 - t423) * t415) * t459 - t263 * t284 - t262 * t285 + t289 * t11 + (-t157 / 0.2e1 + t146 / 0.2e1) * t65 + m(4) * (-pkin(2) * t339 + t168 * t511) + (t265 / 0.2e1 - t240 / 0.2e1) * t99 + t58 * (-mrSges(7,1) * t258 + mrSges(7,2) * t259) + (t254 - t176) * t225 + (t135 * t314 - t136 * t313) * mrSges(5,3) + (-mrSges(7,1) * t471 + mrSges(7,2) * t472) * t101 + (-t15 * t472 + t16 * t471 + t2 * t258 - t259 * t3) * mrSges(7,3) - m(4) * (t235 * t262 + t236 * t263 + t306 * t355) - m(5) * (t135 * t176 + t136 * t177 + t217 * t246) - t246 * t201 + t219 * t84 + t220 * t83 - pkin(2) * t208 + t114 * t28 + t115 * t29 + (t192 / 0.2e1 + t87 * t518 - t169 * mrSges(4,3) + t86 * t519 + t160 * t436 + t495 / 0.2e1 - t493 / 0.2e1 + t339 * mrSges(4,2) + t435 * t540 + t433 * t539 + (-t409 * t68 - t413 * t69) * mrSges(5,3) + (-m(4) * t169 + m(5) * t160 + t113 - t250) * pkin(9) + (t430 * t525 + t432 * t534 + t434 * t532 + t413 * t541 + t172 * t519 + t217 * t437 + (t135 * t409 - t136 * t413) * mrSges(5,3)) * qJD(4)) * t410 + m(5) * (t135 * t254 + t136 * t253 + t334 * t69 + t335 * t68) + (t12 * t349 - t13 * t348 + t465 * t55 - t466 * t56) * mrSges(6,3) + t93 * (mrSges(6,1) * t348 - mrSges(6,2) * t349) + (-Ifges(6,1) * t349 - Ifges(6,4) * t348) * t546 + (-Ifges(6,4) * t349 - Ifges(6,2) * t348) * t547 + t334 * t149 + t335 * t150 - t338 * mrSges(3,2) - t339 * mrSges(3,1) - t352 * t351 + t379 * t60 + (Ifges(6,5) * t266 + Ifges(6,6) * t265) * t524 + (Ifges(6,5) * t241 + Ifges(6,6) * t240) * t525 + (Ifges(5,5) * t314 + Ifges(5,6) * t313) * t525 + (Ifges(7,5) * t145 + Ifges(7,6) * t146) * t527 + (Ifges(7,5) * t158 + Ifges(7,6) * t157) * t528 + (Ifges(5,1) * t314 + Ifges(5,4) * t313) * t532 + (Ifges(5,4) * t314 + Ifges(5,2) * t313) * t534 + (Ifges(6,1) * t266 + Ifges(6,4) * t265) * t535 + (Ifges(6,1) * t241 + Ifges(6,4) * t240) * t536 + (Ifges(6,4) * t266 + Ifges(6,2) * t265) * t537 + (Ifges(6,4) * t241 + Ifges(6,2) * t240) * t538 + t313 * t541 + (Ifges(7,1) * t145 + Ifges(7,4) * t146) * t542 + (Ifges(7,1) * t158 + Ifges(7,4) * t157) * t543 + (Ifges(7,4) * t145 + Ifges(7,2) * t146) * t544 + (Ifges(7,4) * t158 + Ifges(7,2) * t157) * t545 - t349 * t552 - t348 * t553 + (Ifges(7,4) * t259 + Ifges(7,2) * t258) * t554 + (Ifges(7,1) * t259 + Ifges(7,4) * t258) * t555 + t259 * t556 + (-t51 / 0.2e1 - t85 / 0.2e1 - t8 / 0.2e1 + t191 / 0.2e1 - t182 / 0.2e1 - t181 / 0.2e1 - t103 / 0.2e1 - t104 / 0.2e1 - t339 * mrSges(4,1) + pkin(9) * t251 + t168 * mrSges(4,3) + t494 / 0.2e1 - t34 / 0.2e1 - t33 / 0.2e1 + (-Ifges(7,3) / 0.2e1 - Ifges(4,2) / 0.2e1 - t450) * t281 + t569) * t414 + (((-m(4) * t235 + m(5) * t217 - t463) * pkin(9) + t320 / 0.2e1 + t489 / 0.2e1 - t575) * t414 + ((-m(4) * t236 - t284) * pkin(9) - t558) * t410) * qJD(3) + t576 * t122 + t577 * t71 + t578 * t163 + t579 * t164 + (t12 * t219 + t13 * t220 + t167 * t576 + t379 * t93 + t55 * t579 + t56 * t578) * m(6) + t382 + t581 * t95 + t582 * t94 + (t101 * t577 + t114 * t3 + t115 * t2 + t15 * t581 + t16 * t582 + t289 * t58) * m(7) + t258 * t557; t594 * t163 + t595 * t164 + (Ifges(6,5) * t372 + Ifges(7,5) * t287 - Ifges(6,6) * t428 + Ifges(7,6) * t286 + t430) * t530 + t93 * (mrSges(6,1) * t428 + mrSges(6,2) * t372) + (Ifges(6,1) * t372 - Ifges(6,4) * t428) * t546 + (Ifges(6,4) * t372 - Ifges(6,2) * t428) * t547 - t428 * t553 + t463 * t236 - m(6) * (t167 * t198 + t55 * t77 + t56 * t78) + (t210 / 0.2e1 - t156 / 0.2e1) * t66 + t304 * t84 + t305 * t83 + (-pkin(3) * t160 - t135 * t165 - t136 * t166 - t217 * t236) * m(5) + (mrSges(7,1) * t470 - mrSges(7,2) * t469) * t101 + (t15 * t469 - t16 * t470 + t2 * t286 - t287 * t3) * mrSges(7,3) - t235 * t284 + t58 * (-mrSges(7,1) * t286 + mrSges(7,2) * t287) + (t211 / 0.2e1 - t155 / 0.2e1) * t65 + t449 + t501 - t502 + (-Ifges(6,5) * t239 - Ifges(6,6) * t238) * t525 + m(6) * (t12 * t304 + t13 * t305 + t267 * t55 + t268 * t56 + t401 * t93) + (-t12 * t372 - t13 * t428 + t467 * t55 + t56 * t599) * mrSges(6,3) + (-mrSges(6,1) * t599 - mrSges(6,2) * t467) * t167 + (-Ifges(6,1) * t239 - Ifges(6,4) * t238) * t536 + (-Ifges(6,4) * t239 - Ifges(6,2) * t238) * t538 - t166 * t224 - t165 * t225 - t198 * t122 + t188 * t28 + t189 * t29 - t160 * t437 - pkin(3) * t113 + t416 * t327 + t336 * t11 + (t238 / 0.2e1 - t362 / 0.2e1) * t99 + (t239 / 0.2e1 - t363 / 0.2e1) * t100 + (-Ifges(6,5) * t363 - Ifges(6,6) * t362) * t524 + (-Ifges(6,1) * t363 - Ifges(6,4) * t362) * t535 + (-Ifges(6,4) * t363 - Ifges(6,2) * t362) * t537 + t401 * t60 + t86 * t518 + (Ifges(7,5) * t210 + Ifges(7,6) * t211) * t527 + (Ifges(7,5) * t156 + Ifges(7,6) * t155) * t528 + t432 * t539 + t434 * t540 + (Ifges(7,1) * t210 + Ifges(7,4) * t211) * t542 + (Ifges(7,1) * t156 + Ifges(7,4) * t155) * t543 + (Ifges(7,4) * t210 + Ifges(7,2) * t211) * t544 + (Ifges(7,4) * t156 + Ifges(7,2) * t155) * t545 + t409 * t551 + t372 * t552 + (Ifges(7,4) * t287 + Ifges(7,2) * t286) * t554 + (Ifges(7,1) * t287 + Ifges(7,4) * t286) * t555 + t287 * t556 + ((m(6) * t167 + t122) * t514 - t585) * qJD(4) + t567 * t95 + t568 * t94 + (t101 * t561 + t15 * t567 + t16 * t568 + t188 * t3 + t189 * t2 + t336 * t58) * m(7) + (t526 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t327 + t575) * t326 + t560 * mrSges(5,3) + (m(5) * t560 + (-m(5) * t429 - t409 * t224 - t413 * t225) * qJD(4) - t149 * t409 + t150 * t413) * pkin(10) + t561 * t71 + t286 * t557; t120 * t592 + t598 + (-Ifges(6,2) * t196 + t100 + t566) * t538 + (Ifges(5,5) * t273 + Ifges(6,5) * t441 - Ifges(5,6) * t274 - Ifges(6,6) * t196) * t525 - t167 * (mrSges(6,1) * t196 + mrSges(6,2) * t441) + (t196 * t56 + t441 * t55) * mrSges(6,3) + t196 * t573 + (-t167 * t516 - t55 * t62 - t56 * t63 + (t12 * t406 + t13 * t404) * pkin(4)) * m(6) + t559 - t217 * (mrSges(5,1) * t274 + mrSges(5,2) * t273) - t569 + (t135 * t273 + t136 * t274) * mrSges(5,3) + (-Ifges(5,2) * t274 + t172 + t271) * t534 - t135 * t224 + t136 * t225 - t62 * t164 - t63 * t163 - t152 * t71 + (-t122 * t274 + t404 * t83 + t406 * t84) * pkin(4) + t358 * t28 + t359 * t29 + t564 * t94 + t565 * t95 + (-t101 * t152 + t15 * t565 + t16 * t564 + t2 * t359 + t3 * t358) * m(7) + t171 * t531 + (Ifges(5,1) * t273 - t497) * t532 + (Ifges(6,1) * t441 - t580) * t536; t120 * t95 - t574 * t94 - t441 * t163 + t196 * t164 + t11 + t60 + (t120 * t15 - t16 * t574 + t58) * m(7) + (t196 * t55 - t441 * t56 + t93) * m(6); -t15 * t94 + t16 * t95 + t65 * t542 + t583 + t598 + t8;];
tauc  = t1(:);
