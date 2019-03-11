% Calculate time derivative of joint inertia matrix for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:25
% EndTime: 2019-03-10 04:49:00
% DurationCPUTime: 16.14s
% Computational Cost: add. (37969->986), mult. (110656->1434), div. (0->0), fcn. (114950->14), ass. (0->402)
t390 = sin(qJ(6));
t395 = cos(qJ(6));
t396 = cos(qJ(5));
t550 = (t390 ^ 2 + t395 ^ 2) * t396;
t387 = sin(pkin(6));
t558 = 0.2e1 * t387;
t399 = cos(qJ(2));
t389 = cos(pkin(6));
t508 = pkin(1) * t389;
t377 = t399 * t508;
t371 = qJD(2) * t377;
t394 = sin(qJ(2));
t388 = cos(pkin(7));
t424 = t387 * (-pkin(10) * t388 - pkin(9));
t411 = t394 * t424;
t284 = pkin(2) * t389 + t377 + t411;
t481 = t284 * t388;
t557 = qJD(2) * t411 + qJD(3) * t481 + t371;
t391 = sin(qJ(5));
t392 = sin(qJ(4));
t397 = cos(qJ(4));
t343 = t391 * t392 - t396 * t397;
t549 = qJD(4) + qJD(5);
t282 = t549 * t343;
t344 = t391 * t397 + t392 * t396;
t453 = qJD(6) * t395;
t408 = -t282 * t390 + t344 * t453;
t454 = qJD(6) * t390;
t482 = t282 * t395;
t407 = t344 * t454 + t482;
t283 = t549 * t344;
t458 = qJD(4) * t392;
t451 = pkin(4) * t458;
t202 = pkin(5) * t283 + pkin(13) * t282 + t451;
t381 = -pkin(4) * t397 - pkin(3);
t266 = pkin(5) * t343 - pkin(13) * t344 + t381;
t537 = -pkin(12) - pkin(11);
t366 = t537 * t397;
t446 = t537 * t392;
t304 = -t396 * t366 + t391 * t446;
t205 = t266 * t395 - t304 * t390;
t303 = -t391 * t366 - t396 * t446;
t443 = qJD(4) * t537;
t354 = t392 * t443;
t427 = t397 * t443;
t221 = -qJD(5) * t303 + t396 * t354 + t391 * t427;
t103 = qJD(6) * t205 + t202 * t390 + t221 * t395;
t206 = t266 * t390 + t304 * t395;
t104 = -qJD(6) * t206 + t202 * t395 - t221 * t390;
t556 = t103 * t395 - t104 * t390;
t386 = sin(pkin(7));
t393 = sin(qJ(3));
t460 = qJD(3) * t393;
t440 = t386 * t460;
t472 = t388 * t393;
t398 = cos(qJ(3));
t474 = t386 * t398;
t338 = pkin(2) * t472 + pkin(10) * t474;
t318 = pkin(11) * t388 + t338;
t319 = (-pkin(3) * t398 - pkin(11) * t393 - pkin(2)) * t386;
t243 = t397 * t318 + t392 * t319;
t461 = qJD(3) * t386;
t325 = (pkin(3) * t393 - pkin(11) * t398) * t461;
t475 = t386 * t393;
t372 = pkin(10) * t475;
t471 = t388 * t398;
t336 = pkin(2) * t471 - t372;
t326 = t336 * qJD(3);
t183 = -qJD(4) * t243 + t397 * t325 - t326 * t392;
t331 = t388 * t397 - t392 * t475;
t459 = qJD(3) * t398;
t439 = t386 * t459;
t293 = qJD(4) * t331 + t397 * t439;
t153 = pkin(4) * t440 - pkin(12) * t293 + t183;
t457 = qJD(4) * t397;
t182 = -t318 * t458 + t319 * t457 + t392 * t325 + t397 * t326;
t332 = t388 * t392 + t397 * t475;
t294 = -qJD(4) * t332 - t392 * t439;
t160 = pkin(12) * t294 + t182;
t242 = -t318 * t392 + t397 * t319;
t198 = -pkin(4) * t474 - pkin(12) * t332 + t242;
t213 = pkin(12) * t331 + t243;
t455 = qJD(5) * t396;
t456 = qJD(5) * t391;
t62 = t391 * t153 + t396 * t160 + t198 * t455 - t213 * t456;
t58 = pkin(13) * t440 + t62;
t553 = t391 * t198 + t396 * t213;
t129 = -pkin(13) * t474 + t553;
t253 = t331 * t391 + t332 * t396;
t317 = t372 + (-pkin(2) * t398 - pkin(3)) * t388;
t258 = -pkin(4) * t331 + t317;
t412 = t396 * t331 - t332 * t391;
t164 = -pkin(5) * t412 - pkin(13) * t253 + t258;
t85 = -t129 * t390 + t164 * t395;
t175 = qJD(5) * t412 + t293 * t396 + t294 * t391;
t176 = qJD(5) * t253 + t293 * t391 - t396 * t294;
t327 = t338 * qJD(3);
t241 = -pkin(4) * t294 + t327;
t91 = pkin(5) * t176 - pkin(13) * t175 + t241;
t22 = qJD(6) * t85 + t390 * t91 + t395 * t58;
t86 = t129 * t395 + t164 * t390;
t23 = -qJD(6) * t86 - t390 * t58 + t395 * t91;
t555 = t22 * t395 - t23 * t390;
t464 = t398 * t399;
t468 = t393 * t394;
t551 = t388 * t464 - t468;
t224 = t389 * t439 + (t551 * qJD(3) + (-t388 * t468 + t464) * qJD(2)) * t387;
t466 = t394 * t398;
t467 = t393 * t399;
t409 = t388 * t467 + t466;
t277 = t387 * t409 + t389 * t475;
t473 = t387 * t399;
t330 = -t386 * t473 + t389 * t388;
t230 = t277 * t397 + t330 * t392;
t462 = qJD(2) * t387;
t442 = t394 * t462;
t426 = t386 * t442;
t150 = -qJD(4) * t230 - t224 * t392 + t397 * t426;
t229 = -t277 * t392 + t330 * t397;
t151 = qJD(4) * t229 + t224 * t397 + t392 * t426;
t414 = t396 * t229 - t230 * t391;
t68 = qJD(5) * t414 + t150 * t391 + t151 * t396;
t169 = t229 * t391 + t230 * t396;
t69 = qJD(5) * t169 - t396 * t150 + t151 * t391;
t376 = t394 * t508;
t288 = (t399 * t424 - t376) * qJD(2);
t507 = pkin(10) * t386;
t314 = (pkin(2) * t394 - t399 * t507) * t462;
t339 = pkin(9) * t473 + t376;
t267 = (t386 * t389 + t388 * t473) * pkin(10) + t339;
t310 = (-pkin(2) * t399 - t394 * t507 - pkin(1)) * t387;
t428 = -t267 * t459 - t310 * t440 - t393 * t557;
t101 = -t288 * t471 + (-pkin(3) * t442 - t314 * t398) * t386 - t428;
t71 = -pkin(4) * t150 + t101;
t16 = pkin(5) * t69 - pkin(13) * t68 + t71;
t276 = -t387 * t551 - t389 * t474;
t216 = -t284 * t386 + t388 * t310;
t158 = pkin(3) * t276 - pkin(11) * t277 + t216;
t179 = t398 * t267 + t284 * t472 + t310 * t475;
t163 = pkin(11) * t330 + t179;
t92 = t397 * t158 - t163 * t392;
t74 = pkin(4) * t276 - pkin(12) * t230 + t92;
t93 = t392 * t158 + t397 * t163;
t80 = pkin(12) * t229 + t93;
t504 = t391 * t74 + t396 * t80;
t34 = pkin(13) * t276 + t504;
t178 = -t393 * t267 + t398 * (t310 * t386 + t481);
t162 = -pkin(3) * t330 - t178;
t116 = -pkin(4) * t229 + t162;
t60 = -pkin(5) * t414 - pkin(13) * t169 + t116;
t17 = -t34 * t390 + t395 * t60;
t223 = t389 * t440 + (t409 * qJD(3) + (t388 * t466 + t467) * qJD(2)) * t387;
t105 = -t267 * t460 + t288 * t472 + t310 * t439 + t314 * t475 + t398 * t557;
t100 = pkin(11) * t426 + t105;
t225 = -t288 * t386 + t388 * t314;
t113 = pkin(3) * t223 - pkin(11) * t224 + t225;
t38 = -qJD(4) * t93 - t100 * t392 + t397 * t113;
t28 = pkin(4) * t223 - pkin(12) * t151 + t38;
t37 = t397 * t100 + t392 * t113 + t158 * t457 - t163 * t458;
t31 = pkin(12) * t150 + t37;
t8 = t391 * t28 + t396 * t31 + t74 * t455 - t456 * t80;
t5 = pkin(13) * t223 + t8;
t2 = qJD(6) * t17 + t16 * t390 + t395 * t5;
t18 = t34 * t395 + t390 * t60;
t3 = -qJD(6) * t18 + t16 * t395 - t390 * t5;
t554 = t2 * t395 - t3 * t390;
t513 = t390 / 0.2e1;
t512 = t392 / 0.2e1;
t510 = t395 / 0.2e1;
t509 = t397 / 0.2e1;
t431 = -t454 / 0.2e1;
t358 = -mrSges(5,1) * t397 + mrSges(5,2) * t392;
t552 = -m(5) * pkin(3) + t358;
t63 = -qJD(5) * t553 + t153 * t396 - t160 * t391;
t9 = -qJD(5) * t504 + t28 * t396 - t31 * t391;
t548 = 2 * m(4);
t547 = 0.2e1 * m(5);
t546 = 2 * m(6);
t545 = 2 * m(7);
t544 = -2 * mrSges(3,3);
t543 = -2 * mrSges(4,3);
t542 = -2 * mrSges(6,3);
t222 = qJD(5) * t304 + t391 * t354 - t396 * t427;
t541 = 0.2e1 * t222;
t540 = 0.2e1 * t303;
t89 = Ifges(6,1) * t169 + Ifges(6,4) * t414 + Ifges(6,5) * t276;
t538 = t89 / 0.2e1;
t124 = -t169 * t390 + t276 * t395;
t536 = t124 / 0.2e1;
t125 = t169 * t395 + t276 * t390;
t535 = t125 / 0.2e1;
t534 = t169 / 0.2e1;
t186 = Ifges(6,1) * t253 + Ifges(6,4) * t412 - Ifges(6,5) * t474;
t533 = t186 / 0.2e1;
t232 = -t253 * t390 - t395 * t474;
t532 = t232 / 0.2e1;
t410 = -t253 * t395 + t390 * t474;
t531 = -t410 / 0.2e1;
t498 = Ifges(7,4) * t395;
t421 = -Ifges(7,2) * t390 + t498;
t237 = Ifges(7,6) * t343 + t344 * t421;
t530 = t237 / 0.2e1;
t499 = Ifges(7,4) * t390;
t422 = Ifges(7,1) * t395 - t499;
t238 = Ifges(7,5) * t343 + t344 * t422;
t529 = t238 / 0.2e1;
t528 = t253 / 0.2e1;
t292 = Ifges(6,1) * t344 - Ifges(6,4) * t343;
t527 = t292 / 0.2e1;
t382 = Ifges(7,5) * t453;
t526 = Ifges(7,6) * t431 + t382 / 0.2e1;
t350 = t421 * qJD(6);
t525 = t350 / 0.2e1;
t500 = Ifges(5,4) * t397;
t351 = (-Ifges(5,2) * t392 + t500) * qJD(4);
t524 = t351 / 0.2e1;
t352 = t422 * qJD(6);
t523 = t352 / 0.2e1;
t501 = Ifges(5,4) * t392;
t353 = (Ifges(5,1) * t397 - t501) * qJD(4);
t522 = t353 / 0.2e1;
t521 = Ifges(7,5) * t513 + Ifges(7,6) * t510;
t361 = Ifges(7,2) * t395 + t499;
t519 = t361 / 0.2e1;
t362 = Ifges(5,2) * t397 + t501;
t518 = t362 / 0.2e1;
t363 = Ifges(7,1) * t390 + t498;
t517 = t363 / 0.2e1;
t364 = Ifges(5,1) * t392 + t500;
t516 = t364 / 0.2e1;
t515 = t388 / 0.2e1;
t514 = -t390 / 0.2e1;
t511 = -t395 / 0.2e1;
t503 = Ifges(4,4) * t393;
t502 = Ifges(4,4) * t398;
t497 = Ifges(7,6) * t390;
t496 = pkin(4) * qJD(5);
t494 = t225 * mrSges(4,1);
t493 = t225 * mrSges(4,2);
t328 = -pkin(9) * t442 + t371;
t491 = t328 * mrSges(3,2);
t329 = t339 * qJD(2);
t490 = t329 * mrSges(3,1);
t489 = t391 * mrSges(6,1);
t488 = t396 * mrSges(6,2);
t127 = mrSges(6,1) * t276 - mrSges(6,3) * t169;
t78 = -mrSges(7,1) * t124 + mrSges(7,2) * t125;
t487 = -t127 + t78;
t484 = t222 * t303;
t480 = t303 * t391;
t479 = t344 * t390;
t478 = t344 * t395;
t470 = t390 * t396;
t357 = -mrSges(7,1) * t395 + mrSges(7,2) * t390;
t469 = t391 * t357;
t465 = t395 * t396;
t171 = -mrSges(7,1) * t232 - mrSges(7,2) * t410;
t240 = -mrSges(6,1) * t474 - mrSges(6,3) * t253;
t463 = t171 - t240;
t210 = -Ifges(6,5) * t282 - Ifges(6,6) * t283;
t43 = qJD(6) * t124 + t223 * t390 + t395 * t68;
t44 = -qJD(6) * t125 + t223 * t395 - t390 * t68;
t12 = Ifges(7,5) * t43 + Ifges(7,6) * t44 + Ifges(7,3) * t69;
t25 = Ifges(6,5) * t68 - Ifges(6,6) * t69 + Ifges(6,3) * t223;
t26 = Ifges(6,4) * t68 - Ifges(6,2) * t69 + Ifges(6,6) * t223;
t450 = -t26 / 0.2e1 + t12 / 0.2e1;
t117 = qJD(6) * t232 + t175 * t395 + t390 * t440;
t118 = qJD(6) * t410 - t175 * t390 + t395 * t440;
t47 = Ifges(7,5) * t117 + Ifges(7,6) * t118 + Ifges(7,3) * t176;
t96 = Ifges(6,4) * t175 - Ifges(6,2) * t176 + Ifges(6,6) * t440;
t449 = t47 / 0.2e1 - t96 / 0.2e1;
t54 = Ifges(7,5) * t125 + Ifges(7,6) * t124 - Ifges(7,3) * t414;
t88 = Ifges(6,4) * t169 + Ifges(6,2) * t414 + Ifges(6,6) * t276;
t448 = t54 / 0.2e1 - t88 / 0.2e1;
t15 = -mrSges(7,1) * t44 + mrSges(7,2) * t43;
t6 = -pkin(5) * t223 - t9;
t447 = m(7) * t6 + t15;
t75 = Ifges(5,5) * t151 + Ifges(5,6) * t150 + Ifges(5,3) * t223;
t95 = Ifges(6,5) * t175 - Ifges(6,6) * t176 + Ifges(6,3) * t440;
t147 = Ifges(4,5) * t224 - Ifges(4,6) * t223 + Ifges(4,3) * t426;
t199 = Ifges(5,5) * t293 + Ifges(5,6) * t294 + Ifges(5,3) * t440;
t59 = -pkin(5) * t440 - t63;
t70 = -mrSges(7,1) * t118 + mrSges(7,2) * t117;
t444 = m(7) * t59 + t70;
t130 = -Ifges(7,5) * t410 + Ifges(7,6) * t232 - Ifges(7,3) * t412;
t185 = Ifges(6,4) * t253 + Ifges(6,2) * t412 - Ifges(6,6) * t474;
t436 = -t185 / 0.2e1 + t130 / 0.2e1;
t140 = -Ifges(7,5) * t407 - Ifges(7,6) * t408 + Ifges(7,3) * t283;
t211 = -Ifges(6,4) * t282 - Ifges(6,2) * t283;
t435 = -t211 / 0.2e1 + t140 / 0.2e1;
t236 = Ifges(7,3) * t343 + (Ifges(7,5) * t395 - t497) * t344;
t291 = Ifges(6,4) * t344 - Ifges(6,2) * t343;
t434 = -t291 / 0.2e1 + t236 / 0.2e1;
t383 = Ifges(5,5) * t457;
t433 = -Ifges(5,6) * t458 / 0.2e1 + t383 / 0.2e1 + t210 / 0.2e1;
t432 = Ifges(5,5) * t512 + Ifges(5,6) * t509 + Ifges(6,5) * t344 / 0.2e1 - Ifges(6,6) * t343 / 0.2e1;
t430 = t453 / 0.2e1;
t177 = mrSges(7,1) * t408 - mrSges(7,2) * t407;
t429 = m(7) * t222 + t177;
t425 = mrSges(7,3) * t550;
t423 = mrSges(7,1) * t390 + mrSges(7,2) * t395;
t419 = t37 * t397 - t38 * t392;
t35 = -t391 * t80 + t396 * t74;
t416 = t182 * t397 - t183 * t392;
t138 = t198 * t396 - t213 * t391;
t406 = t395 * t350 + t390 * t352 - t361 * t454 + t363 * t453;
t19 = mrSges(7,1) * t69 - mrSges(7,3) * t43;
t20 = -mrSges(7,2) * t69 + mrSges(7,3) * t44;
t83 = mrSges(7,2) * t414 + mrSges(7,3) * t124;
t84 = -mrSges(7,1) * t414 - mrSges(7,3) * t125;
t405 = -t84 * t453 - t83 * t454 - t390 * t19 + t395 * t20 + m(7) * (-t17 * t453 - t18 * t454 + t554);
t180 = mrSges(7,2) * t412 + mrSges(7,3) * t232;
t181 = -mrSges(7,1) * t412 + mrSges(7,3) * t410;
t81 = mrSges(7,1) * t176 - mrSges(7,3) * t117;
t82 = -mrSges(7,2) * t176 + mrSges(7,3) * t118;
t404 = -t181 * t453 - t180 * t454 - t390 * t81 + t395 * t82 + m(7) * (-t453 * t85 - t454 * t86 + t555);
t189 = mrSges(7,1) * t283 + mrSges(7,3) * t407;
t190 = -mrSges(7,2) * t283 - mrSges(7,3) * t408;
t271 = -mrSges(7,2) * t343 - mrSges(7,3) * t479;
t272 = mrSges(7,1) * t343 - mrSges(7,3) * t478;
t403 = -t272 * t453 - t271 * t454 - t390 * t189 + t395 * t190 + m(7) * (-t205 * t453 - t206 * t454 + t556);
t13 = Ifges(7,4) * t43 + Ifges(7,2) * t44 + Ifges(7,6) * t69;
t14 = Ifges(7,1) * t43 + Ifges(7,4) * t44 + Ifges(7,5) * t69;
t33 = -pkin(5) * t276 - t35;
t346 = t423 * qJD(6);
t55 = Ifges(7,4) * t125 + Ifges(7,2) * t124 - Ifges(7,6) * t414;
t56 = Ifges(7,1) * t125 + Ifges(7,4) * t124 - Ifges(7,5) * t414;
t402 = t9 * mrSges(6,1) - t8 * mrSges(6,2) + t124 * t525 + t125 * t523 + t13 * t510 + t14 * t513 - t414 * t526 + t33 * t346 + t6 * t357 + t43 * t517 + t56 * t430 + t55 * t431 + t44 * t519 + t69 * t521 + t25 + ((-t17 * t395 - t18 * t390) * qJD(6) + t554) * mrSges(7,3);
t128 = pkin(5) * t474 - t138;
t131 = -Ifges(7,4) * t410 + Ifges(7,2) * t232 - Ifges(7,6) * t412;
t132 = -Ifges(7,1) * t410 + Ifges(7,4) * t232 - Ifges(7,5) * t412;
t48 = Ifges(7,4) * t117 + Ifges(7,2) * t118 + Ifges(7,6) * t176;
t49 = Ifges(7,1) * t117 + Ifges(7,4) * t118 + Ifges(7,5) * t176;
t401 = t63 * mrSges(6,1) - t62 * mrSges(6,2) + t117 * t517 + t118 * t519 + t128 * t346 + t131 * t431 + t132 * t430 + t176 * t521 + t232 * t525 - t410 * t523 - t412 * t526 + t59 * t357 + t48 * t510 + t49 * t513 + t95 + ((-t390 * t86 - t395 * t85) * qJD(6) + t555) * mrSges(7,3);
t141 = -Ifges(7,4) * t407 - Ifges(7,2) * t408 + Ifges(7,6) * t283;
t142 = -Ifges(7,1) * t407 - Ifges(7,4) * t408 + Ifges(7,5) * t283;
t400 = t142 * t513 + t141 * t510 + t238 * t430 - t482 * t517 + t283 * t521 + t303 * t346 - t350 * t479 / 0.2e1 + t478 * t523 + t343 * t526 - t221 * mrSges(6,2) + t210 - t408 * t361 / 0.2e1 + (t344 * t363 + t237) * t431 + (t357 - mrSges(6,1)) * t222 + ((-t205 * t395 - t206 * t390) * qJD(6) + t556) * mrSges(7,3);
t380 = -pkin(4) * t396 - pkin(5);
t379 = pkin(4) * t391 + pkin(13);
t370 = Ifges(3,5) * t399 * t462;
t369 = Ifges(4,5) * t439;
t347 = (mrSges(5,1) * t392 + mrSges(5,2) * t397) * qJD(4);
t342 = -mrSges(4,2) * t388 + mrSges(4,3) * t474;
t341 = mrSges(4,1) * t388 - mrSges(4,3) * t475;
t337 = -pkin(9) * t387 * t394 + t377;
t324 = (Ifges(4,1) * t398 - t503) * t461;
t323 = (-Ifges(4,2) * t393 + t502) * t461;
t322 = -Ifges(4,6) * t440 + t369;
t321 = (mrSges(4,1) * t393 + mrSges(4,2) * t398) * t461;
t316 = Ifges(4,5) * t388 + (Ifges(4,1) * t393 + t502) * t386;
t315 = Ifges(4,6) * t388 + (Ifges(4,2) * t398 + t503) * t386;
t302 = -mrSges(5,1) * t474 - mrSges(5,3) * t332;
t301 = mrSges(5,2) * t474 + mrSges(5,3) * t331;
t289 = mrSges(6,1) * t343 + mrSges(6,2) * t344;
t264 = t423 * t344;
t259 = -mrSges(5,1) * t331 + mrSges(5,2) * t332;
t251 = -mrSges(5,2) * t440 + mrSges(5,3) * t294;
t250 = mrSges(5,1) * t440 - mrSges(5,3) * t293;
t248 = Ifges(5,1) * t332 + Ifges(5,4) * t331 - Ifges(5,5) * t474;
t247 = Ifges(5,4) * t332 + Ifges(5,2) * t331 - Ifges(5,6) * t474;
t246 = Ifges(5,5) * t332 + Ifges(5,6) * t331 - Ifges(5,3) * t474;
t239 = mrSges(6,2) * t474 + mrSges(6,3) * t412;
t235 = mrSges(4,1) * t330 - mrSges(4,3) * t277;
t234 = -mrSges(4,2) * t330 - mrSges(4,3) * t276;
t214 = -mrSges(5,1) * t294 + mrSges(5,2) * t293;
t212 = -Ifges(6,1) * t282 - Ifges(6,4) * t283;
t209 = mrSges(6,1) * t283 - mrSges(6,2) * t282;
t201 = Ifges(5,1) * t293 + Ifges(5,4) * t294 + Ifges(5,5) * t440;
t200 = Ifges(5,4) * t293 + Ifges(5,2) * t294 + Ifges(5,6) * t440;
t197 = mrSges(4,1) * t426 - mrSges(4,3) * t224;
t196 = -mrSges(4,2) * t426 - mrSges(4,3) * t223;
t193 = -mrSges(6,1) * t412 + mrSges(6,2) * t253;
t192 = Ifges(4,1) * t277 - Ifges(4,4) * t276 + Ifges(4,5) * t330;
t191 = Ifges(4,4) * t277 - Ifges(4,2) * t276 + Ifges(4,6) * t330;
t188 = mrSges(5,1) * t276 - mrSges(5,3) * t230;
t187 = -mrSges(5,2) * t276 + mrSges(5,3) * t229;
t184 = Ifges(6,5) * t253 + Ifges(6,6) * t412 - Ifges(6,3) * t474;
t170 = -mrSges(5,1) * t229 + mrSges(5,2) * t230;
t167 = -mrSges(6,2) * t440 - mrSges(6,3) * t176;
t166 = mrSges(6,1) * t440 - mrSges(6,3) * t175;
t161 = mrSges(4,1) * t223 + mrSges(4,2) * t224;
t149 = Ifges(4,1) * t224 - Ifges(4,4) * t223 + Ifges(4,5) * t426;
t148 = Ifges(4,4) * t224 - Ifges(4,2) * t223 + Ifges(4,6) * t426;
t135 = Ifges(5,1) * t230 + Ifges(5,4) * t229 + Ifges(5,5) * t276;
t134 = Ifges(5,4) * t230 + Ifges(5,2) * t229 + Ifges(5,6) * t276;
t133 = Ifges(5,5) * t230 + Ifges(5,6) * t229 + Ifges(5,3) * t276;
t126 = -mrSges(6,2) * t276 + mrSges(6,3) * t414;
t108 = mrSges(5,1) * t223 - mrSges(5,3) * t151;
t107 = -mrSges(5,2) * t223 + mrSges(5,3) * t150;
t106 = (t288 * t388 + t314 * t386) * t398 + t428;
t98 = mrSges(6,1) * t176 + mrSges(6,2) * t175;
t97 = Ifges(6,1) * t175 - Ifges(6,4) * t176 + Ifges(6,5) * t440;
t94 = -mrSges(6,1) * t414 + mrSges(6,2) * t169;
t90 = -mrSges(5,1) * t150 + mrSges(5,2) * t151;
t87 = Ifges(6,5) * t169 + Ifges(6,6) * t414 + Ifges(6,3) * t276;
t77 = Ifges(5,1) * t151 + Ifges(5,4) * t150 + Ifges(5,5) * t223;
t76 = Ifges(5,4) * t151 + Ifges(5,2) * t150 + Ifges(5,6) * t223;
t51 = -mrSges(6,2) * t223 - mrSges(6,3) * t69;
t50 = mrSges(6,1) * t223 - mrSges(6,3) * t68;
t30 = mrSges(6,1) * t69 + mrSges(6,2) * t68;
t27 = Ifges(6,1) * t68 - Ifges(6,4) * t69 + Ifges(6,5) * t223;
t1 = [(0.2e1 * (t328 * t399 + t329 * t394) * mrSges(3,3) + ((t337 * t544 + Ifges(3,5) * t389 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t399) * t558) * t399 + (t386 * (Ifges(4,5) * t277 - Ifges(4,6) * t276 + Ifges(4,3) * t330) + t339 * t544 - 0.2e1 * Ifges(3,6) * t389 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t394 + (Ifges(3,1) - Ifges(3,2)) * t399) * t558) * t394) * qJD(2)) * t387 + 0.2e1 * t504 * t51 + (t116 * t71 + t35 * t9 + t504 * t8) * t546 + 0.2e1 * m(3) * (t328 * t339 - t329 * t337) + (t370 - 0.2e1 * t490 - 0.2e1 * t491) * t389 + t330 * t147 + 0.2e1 * t105 * t234 + 0.2e1 * t106 * t235 + t229 * t76 + t230 * t77 + t224 * t192 + 0.2e1 * t216 * t161 + 0.2e1 * t179 * t196 + 0.2e1 * t178 * t197 + 0.2e1 * t37 * t187 + 0.2e1 * t38 * t188 + t169 * t27 + 0.2e1 * t101 * t170 + 0.2e1 * t162 * t90 + t150 * t134 + t151 * t135 + 0.2e1 * t8 * t126 + 0.2e1 * t9 * t127 + t124 * t13 + t125 * t14 + 0.2e1 * t116 * t30 + 0.2e1 * t93 * t107 + 0.2e1 * t92 * t108 + 0.2e1 * t71 * t94 + t68 * t89 + 0.2e1 * t2 * t83 + 0.2e1 * t3 * t84 + 0.2e1 * t6 * t78 + t44 * t55 + t43 * t56 + 0.2e1 * t35 * t50 + 0.2e1 * t33 * t15 + 0.2e1 * t18 * t20 + 0.2e1 * t17 * t19 + (-t88 + t54) * t69 + (t149 + 0.2e1 * t493) * t277 + (-t148 + t25 + t75 + 0.2e1 * t494) * t276 + (t17 * t3 + t18 * t2 + t33 * t6) * t545 + (t101 * t162 + t37 * t93 + t38 * t92) * t547 + (t105 * t179 + t106 * t178 + t216 * t225) * t548 + (t133 + t87 - t191) * t223 - (-t26 + t12) * t414; t504 * t167 - t490 - t491 + m(5) * (t101 * t317 + t162 * t327 + t182 * t93 + t183 * t92 + t242 * t38 + t243 * t37) + m(7) * (t128 * t6 + t17 * t23 + t18 * t22 + t2 * t86 + t3 * t85 + t33 * t59) + t370 + (-t323 / 0.2e1 + t199 / 0.2e1 + t95 / 0.2e1) * t276 + (-t315 / 0.2e1 + t246 / 0.2e1 + t184 / 0.2e1) * t223 + t336 * t197 + t338 * t196 + t106 * t341 + t105 * t342 + t330 * t322 / 0.2e1 + t331 * t76 / 0.2e1 + t332 * t77 / 0.2e1 + t216 * t321 + t277 * t324 / 0.2e1 + t326 * t234 + t224 * t316 / 0.2e1 + t317 * t90 + t293 * t135 / 0.2e1 + t294 * t134 / 0.2e1 + t37 * t301 + t38 * t302 + t258 * t30 + t101 * t259 + t243 * t107 + t150 * t247 / 0.2e1 + t151 * t248 / 0.2e1 + t92 * t250 + t93 * t251 + t8 * t239 + t9 * t240 + t241 * t94 + t242 * t108 + t229 * t200 / 0.2e1 + t230 * t201 / 0.2e1 + t162 * t214 + t182 * t187 + t183 * t188 + t71 * t193 + t2 * t180 + t3 * t181 + t6 * t171 + t35 * t166 + t138 * t50 + t62 * t126 + t63 * t127 + t128 * t15 + t44 * t131 / 0.2e1 + t43 * t132 / 0.2e1 + t116 * t98 + t117 * t56 / 0.2e1 + t118 * t55 / 0.2e1 + t22 * t83 + t23 * t84 + t85 * t19 + t86 * t20 + t17 * t81 + t18 * t82 + t59 * t78 + t33 * t70 + (t170 - t235) * t327 + m(4) * (t105 * t338 + t106 * t336 - t178 * t327 + t179 * t326) + ((t493 + t149 / 0.2e1) * t393 + (-t494 + t148 / 0.2e1 - t75 / 0.2e1 - t25 / 0.2e1) * t398 + (Ifges(4,3) * t515 + (Ifges(4,5) * t393 + Ifges(4,6) * t398) * t386 / 0.2e1) * t442 + (-m(4) * t225 - t161) * pkin(2) + ((-t178 * mrSges(4,3) + t192 / 0.2e1) * t398 + (-t179 * mrSges(4,3) + t133 / 0.2e1 + t87 / 0.2e1 - t191 / 0.2e1) * t393) * qJD(3)) * t386 + t147 * t515 + t27 * t528 + t14 * t531 + t13 * t532 + t68 * t533 + t97 * t534 + t49 * t535 + t48 * t536 + t175 * t538 + t436 * t69 - Ifges(3,6) * t442 + t448 * t176 + m(6) * (t116 * t241 + t138 * t9 + t258 * t71 + t35 * t63 + t504 * t62 + t553 * t8) + t553 * t51 - t449 * t414 - t450 * t412; t388 * t322 + 0.2e1 * t326 * t342 + t331 * t200 + t332 * t201 + 0.2e1 * t317 * t214 + t293 * t248 + t294 * t247 + 0.2e1 * t182 * t301 + 0.2e1 * t183 * t302 + t253 * t97 + 0.2e1 * t258 * t98 + 0.2e1 * t242 * t250 + 0.2e1 * t243 * t251 + 0.2e1 * t62 * t239 + 0.2e1 * t63 * t240 + 0.2e1 * t241 * t193 + t232 * t48 + 0.2e1 * t22 * t180 + 0.2e1 * t23 * t181 + t175 * t186 + 0.2e1 * t59 * t171 + 0.2e1 * t138 * t166 + 0.2e1 * t128 * t70 + t118 * t131 + t117 * t132 + 0.2e1 * t85 * t81 + 0.2e1 * t86 * t82 + (-t185 + t130) * t176 + (-0.2e1 * pkin(2) * t321 + t393 * t324 + (-t199 + t323 - t95) * t398 + ((t336 * t543 + t316) * t398 + (t338 * t543 + t184 + t246 - t315) * t393) * qJD(3)) * t386 + 0.2e1 * (-t341 + t259) * t327 + (t128 * t59 + t22 * t86 + t23 * t85) * t545 + (t182 * t243 + t183 * t242 + t317 * t327) * t547 + (t326 * t338 - t327 * t336) * t548 + 0.2e1 * t553 * t167 + (t138 * t63 + t241 * t258 + t553 * t62) * t546 - t410 * t49 - (-t96 + t47) * t412; t552 * t101 + (-mrSges(6,3) * t504 + t448) * t283 + m(6) * (t116 * t451 + t221 * t504 - t222 * t35 - t303 * t9 + t304 * t8 + t381 * t71) + ((-t392 * t93 - t397 * t92) * qJD(4) + t419) * mrSges(5,3) + t147 + t381 * t30 + t162 * t347 + t304 * t51 + t71 * t289 + t6 * t264 + t2 * t271 + t3 * t272 + t221 * t126 + t205 * t19 + t206 * t20 + t116 * t209 + t17 * t189 + t18 * t190 + t33 * t177 + t104 * t84 - t105 * mrSges(4,2) + t106 * mrSges(4,1) + t103 * t83 - pkin(3) * t90 + m(7) * (t103 * t18 + t104 * t17 + t2 * t206 + t205 * t3 + t222 * t33 + t303 * t6) + (t15 - t50) * t303 - (-t35 * mrSges(6,3) + t510 * t56 + t514 * t55 + t538) * t282 + (t27 / 0.2e1 - t9 * mrSges(6,3) + t13 * t514 + t14 * t510 + (t511 * t55 + t514 * t56) * qJD(6)) * t344 + (t135 * t509 + (pkin(4) * t94 - t134 / 0.2e1) * t392) * qJD(4) + t487 * t222 + (-t188 * t457 - t187 * t458 - t392 * t108 + t397 * t107 + m(5) * (-t457 * t92 - t458 * t93 + t419)) * pkin(11) + t76 * t509 + t77 * t512 + t151 * t516 + t150 * t518 + t230 * t522 + t229 * t524 + t68 * t527 + t43 * t529 + t44 * t530 + t212 * t534 + t142 * t535 + t141 * t536 + t432 * t223 + t433 * t276 + t434 * t69 + (-t8 * mrSges(6,3) + t450) * t343 - t435 * t414; (-mrSges(4,1) + t552) * t327 + ((-t242 * t397 - t243 * t392) * qJD(4) + t416) * mrSges(5,3) + (-t166 + t70) * t303 + m(7) * (t103 * t86 + t104 * t85 + t128 * t222 + t205 * t23 + t206 * t22 + t303 * t59) + t369 + t381 * t98 + t317 * t347 - t326 * mrSges(4,2) + t304 * t167 + t241 * t289 + t59 * t264 + t22 * t271 + t23 * t272 + t258 * t209 + t221 * t239 + t205 * t81 + t206 * t82 - pkin(3) * t214 + t85 * t189 + t86 * t190 + t103 * t180 + t104 * t181 + t128 * t177 - (-t138 * mrSges(6,3) + t131 * t514 + t132 * t510 + t533) * t282 + (t97 / 0.2e1 - t63 * mrSges(6,3) + t48 * t514 + t49 * t510 + (t131 * t511 + t132 * t514) * qJD(6)) * t344 + (t248 * t509 + (pkin(4) * t193 - t247 / 0.2e1) * t392) * qJD(4) + t463 * t222 + (-t302 * t457 - t301 * t458 - t392 * t250 + t397 * t251 + m(5) * (-t242 * t457 - t243 * t458 + t416)) * pkin(11) + (-t433 * t398 + (-Ifges(4,6) + t432) * t460) * t386 + t200 * t509 + t201 * t512 + t293 * t516 + t294 * t518 + t332 * t522 + t331 * t524 + t175 * t527 + t212 * t528 + t117 * t529 + t118 * t530 + t142 * t531 + t141 * t532 + t434 * t176 + (-t62 * mrSges(6,3) + t449) * t343 + (-mrSges(6,3) * t553 + t436) * t283 + m(6) * (-t138 * t222 + t221 * t553 + t241 * t381 + t258 * t451 - t303 * t63 + t304 * t62) - t435 * t412; -0.2e1 * pkin(3) * t347 + 0.2e1 * t103 * t271 + 0.2e1 * t104 * t272 + t177 * t540 + 0.2e1 * t205 * t189 + 0.2e1 * t206 * t190 + 0.2e1 * t381 * t209 + t264 * t541 + t397 * t351 + t392 * t353 + (t397 * t364 + (0.2e1 * pkin(4) * t289 - t362) * t392) * qJD(4) + (t103 * t206 + t104 * t205 + t484) * t545 + (t221 * t304 + t381 * t451 + t484) * t546 + (t221 * t542 + t140 - t211) * t343 + (t304 * t542 + t236 - t291) * t283 - (mrSges(6,3) * t540 - t237 * t390 + t238 * t395 + t292) * t282 + (mrSges(6,3) * t541 - t390 * t141 + t395 * t142 + t212 + (-t237 * t395 - t238 * t390) * qJD(6)) * t344; (t391 * t51 + t396 * t50 + m(6) * (t391 * t8 + t396 * t9) + (t487 * t391 + (-t390 * t84 + t395 * t83 + t126) * t396 + m(6) * (-t35 * t391 + t396 * t504) + m(7) * (-t17 * t470 + t18 * t465 + t33 * t391)) * qJD(5)) * pkin(4) + t447 * t380 + t405 * t379 + t402 + t38 * mrSges(5,1) - t37 * mrSges(5,2) + t75; t401 + (t391 * t167 + t396 * t166 + m(6) * (t391 * t62 + t396 * t63) + (t463 * t391 + (t180 * t395 - t181 * t390 + t239) * t396 + m(6) * (-t138 * t391 + t396 * t553) + m(7) * (t128 * t391 + t465 * t86 - t470 * t85)) * qJD(5)) * pkin(4) + t444 * t380 + t404 * t379 - t182 * mrSges(5,2) + t183 * mrSges(5,1) + t199; t400 + t429 * t380 + t403 * t379 + (m(6) * (t221 * t391 - t222 * t396) + (t282 * t396 - t283 * t391) * mrSges(6,3) + ((mrSges(6,3) * t344 + t264) * t391 + (-t343 * mrSges(6,3) + t271 * t395 - t272 * t390) * t396 + m(6) * (t304 * t396 + t480) + m(7) * (-t205 * t470 + t206 * t465 + t480)) * qJD(5)) * pkin(4) + t383 + (-Ifges(5,6) * t392 + pkin(11) * t358) * qJD(4); 0.2e1 * t380 * t346 + (-0.2e1 * t489 - 0.2e1 * t488 + 0.2e1 * t469 + (t379 * t550 + t380 * t391) * t545 + 0.2e1 * t425) * t496 + t406; -pkin(5) * t447 + pkin(13) * t405 + t402; -pkin(5) * t444 + pkin(13) * t404 + t401; -pkin(5) * t429 + pkin(13) * t403 + t400; (t380 - pkin(5)) * t346 + (-t489 - t488 + t469 + m(7) * (-pkin(5) * t391 + pkin(13) * t550) + t425) * t496 + t406; -0.2e1 * pkin(5) * t346 + t406; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t23 - mrSges(7,2) * t22 + t47; mrSges(7,1) * t104 - mrSges(7,2) * t103 + t140; t382 - t423 * pkin(4) * t455 + (t357 * t379 - t497) * qJD(6); t382 + (pkin(13) * t357 - t497) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
