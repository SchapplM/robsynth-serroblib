% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:26
% EndTime: 2019-03-09 09:29:18
% DurationCPUTime: 36.05s
% Computational Cost: add. (8955->869), mult. (21527->1142), div. (0->0), fcn. (15520->10), ass. (0->389)
t267 = qJ(3) - pkin(9);
t557 = m(6) + m(7);
t579 = t267 * t557 - mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t552 = Ifges(3,5) + Ifges(5,5);
t570 = -Ifges(4,4) + t552;
t553 = Ifges(5,4) + Ifges(4,5);
t578 = Ifges(3,6) - t553;
t270 = sin(qJ(5));
t274 = cos(qJ(5));
t404 = qJD(5) * t274;
t271 = sin(qJ(2));
t265 = sin(pkin(6));
t414 = qJD(1) * t265;
t383 = t271 * t414;
t224 = qJ(4) * t383;
t235 = pkin(2) * t383;
t275 = cos(qJ(2));
t382 = t275 * t414;
t111 = -t267 * t382 + t224 + t235;
t266 = cos(pkin(6));
t413 = qJD(1) * t266;
t395 = pkin(1) * t413;
t240 = t275 * t395;
t490 = pkin(3) + pkin(8);
t363 = t265 * (-pkin(4) - t490);
t339 = t271 * t363;
t114 = qJD(1) * t339 + t240;
t61 = t274 * t111 + t270 * t114;
t577 = qJD(3) * t270 + t267 * t404 - t61;
t188 = pkin(8) * t382 + t271 * t395;
t147 = pkin(3) * t382 + t188;
t576 = -qJD(4) - t147;
t575 = mrSges(3,1) - mrSges(4,2);
t574 = Ifges(5,2) + Ifges(4,3);
t366 = pkin(10) * t274 - qJ(4);
t207 = pkin(5) * t270 + pkin(2) - t366;
t573 = pkin(10) * t383 + qJD(6) * t207 + t577;
t340 = pkin(5) * t274 + pkin(10) * t270;
t430 = t267 * t270;
t572 = qJD(5) * t340 - qJD(6) * t430 - (-pkin(4) - t340) * t382 - t576;
t399 = qJD(1) * qJD(2);
t192 = (qJDD(1) * t271 + t275 * t399) * t265;
t397 = qJDD(1) * t266;
t243 = qJDD(2) + t397;
t246 = qJD(2) + t413;
t474 = pkin(1) * t266;
t394 = qJD(2) * t474;
t355 = qJD(1) * t394;
t391 = pkin(1) * t397;
t108 = -pkin(8) * t192 - t271 * t355 + t275 * t391;
t89 = -t243 * pkin(2) + qJDD(3) - t108;
t51 = t192 * pkin(3) - t243 * qJ(4) - t246 * qJD(4) + t89;
t42 = pkin(4) * t192 + t51;
t405 = qJD(5) * t270;
t60 = -t111 * t270 + t114 * t274;
t571 = qJD(3) * t274 - t267 * t405 - t60;
t568 = Ifges(3,3) + Ifges(5,1) + Ifges(4,1);
t565 = -pkin(5) * t383 - t571;
t435 = t265 * t271;
t410 = qJD(2) * t265;
t381 = t271 * t410;
t398 = qJDD(1) * t265;
t564 = -qJD(1) * t381 + t275 * t398;
t269 = sin(qJ(6));
t273 = cos(qJ(6));
t333 = -mrSges(7,1) * t273 + mrSges(7,2) * t269;
t300 = m(7) * pkin(5) - t333;
t335 = mrSges(6,1) * t270 + mrSges(6,2) * t274;
t563 = -t270 * t300 - t335;
t211 = qJD(5) + t382;
t268 = pkin(2) + qJ(4);
t341 = -t268 * t275 - pkin(1);
t439 = qJ(3) * t271;
t117 = (t341 - t439) * t414;
t149 = (-pkin(2) * t275 - pkin(1) - t439) * t414;
t561 = (-(-t271 * t578 + t275 * t570) * t246 / 0.2e1 - t149 * (-mrSges(4,2) * t271 - mrSges(4,3) * t275) - t117 * (-mrSges(5,2) * t275 + mrSges(5,3) * t271)) * t265;
t160 = t246 * t274 + t270 * t383;
t104 = -t160 * t269 + t211 * t273;
t105 = t160 * t273 + t211 * t269;
t159 = t246 * t270 - t274 * t383;
t158 = qJD(6) + t159;
t38 = t105 * Ifges(7,5) + t104 * Ifges(7,6) + t158 * Ifges(7,3);
t115 = -pkin(4) * t382 - t147;
t223 = t246 * qJ(3);
t530 = qJD(4) + t223;
t85 = -pkin(9) * t246 - t115 + t530;
t97 = (-t267 * t271 + t341) * t414;
t45 = t270 * t85 + t274 * t97;
t450 = t160 * Ifges(6,4);
t546 = t211 * Ifges(6,6);
t548 = t159 * Ifges(6,2);
t73 = t450 + t546 - t548;
t560 = -t73 / 0.2e1 + t38 / 0.2e1 - t45 * mrSges(6,3);
t173 = qJDD(5) + t564;
t480 = t173 / 0.2e1;
t80 = -qJD(5) * t159 + t192 * t270 + t243 * t274;
t493 = t80 / 0.2e1;
t559 = Ifges(6,1) * t493 + Ifges(6,5) * t480;
t32 = qJD(6) * t104 + t173 * t269 + t273 * t80;
t500 = t32 / 0.2e1;
t33 = -qJD(6) * t105 + t173 * t273 - t269 * t80;
t499 = t33 / 0.2e1;
t81 = -qJD(5) * t160 + t192 * t274 - t243 * t270;
t77 = qJDD(6) - t81;
t494 = t77 / 0.2e1;
t492 = t81 / 0.2e1;
t558 = -m(5) - m(4);
t551 = -t269 * t573 + t273 * t572;
t550 = t269 * t572 + t273 * t573;
t12 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t58 = mrSges(6,1) * t173 - mrSges(6,3) * t80;
t549 = t58 - t12;
t547 = t211 * Ifges(6,5);
t531 = -t246 * pkin(2) + qJD(3);
t522 = t246 * qJ(4) - t531;
t84 = t114 + t522;
t545 = t84 * (mrSges(6,1) * t274 - mrSges(6,2) * t270);
t544 = mrSges(6,1) + t300;
t392 = m(7) * pkin(10) + mrSges(7,3);
t517 = mrSges(6,2) - t392;
t464 = mrSges(6,3) * t160;
t113 = mrSges(6,1) * t211 - t464;
t55 = -mrSges(7,1) * t104 + mrSges(7,2) * t105;
t444 = -t113 + t55;
t121 = t192 * mrSges(5,1) - t243 * mrSges(5,3);
t34 = -t81 * mrSges(6,1) + t80 * mrSges(6,2);
t543 = t34 - t121;
t359 = mrSges(5,1) * t383;
t180 = -mrSges(5,3) * t246 + t359;
t90 = mrSges(6,1) * t159 + mrSges(6,2) * t160;
t542 = t90 - t180;
t432 = t265 * t275;
t233 = Ifges(3,4) * t382;
t451 = Ifges(5,6) * t275;
t304 = (t271 * Ifges(5,3) - t451) * t265;
t540 = Ifges(3,1) * t383 + qJD(1) * t304 + t246 * t552 + t233;
t232 = Ifges(5,6) * t383;
t454 = Ifges(4,6) * t271;
t305 = (-t275 * Ifges(4,3) - t454) * t265;
t539 = -Ifges(5,2) * t382 + qJD(1) * t305 + t246 * t553 + t232;
t427 = t270 * t275;
t152 = (-t269 * t427 - t271 * t273) * t414;
t377 = t269 * t405;
t401 = qJD(6) * t274;
t538 = t273 * t401 + t152 - t377;
t421 = t273 * t275;
t153 = (-t269 * t271 + t270 * t421) * t414;
t537 = t269 * t401 + t273 * t405 + t153;
t360 = mrSges(4,1) * t382;
t181 = -mrSges(4,3) * t246 - t360;
t358 = mrSges(5,1) * t382;
t182 = mrSges(5,2) * t246 + t358;
t536 = -t182 + t181;
t18 = mrSges(7,1) * t77 - mrSges(7,3) * t32;
t19 = -mrSges(7,2) * t77 + mrSges(7,3) * t33;
t532 = -t269 * t18 + t273 * t19;
t400 = -m(5) - t557;
t528 = -mrSges(5,1) - mrSges(4,1) - mrSges(3,3);
t526 = mrSges(5,3) + t575;
t523 = t192 * t570 + t568 * t243 + t564 * t578;
t461 = Ifges(3,4) * t271;
t307 = (t275 * Ifges(3,2) + t461) * t265;
t521 = t271 * (t160 * Ifges(6,5) + t246 * Ifges(3,6) - t159 * Ifges(6,6) + t211 * Ifges(6,3) + qJD(1) * t307);
t107 = pkin(8) * t564 + t271 * t391 + t275 * t355;
t79 = -t243 * qJ(3) - t246 * qJD(3) - t107;
t279 = qJDD(4) - t79;
t491 = -pkin(3) - pkin(4);
t41 = -pkin(9) * t243 - t491 * t564 + t279;
t407 = qJD(3) * t271;
t443 = pkin(1) * qJDD(1);
t473 = pkin(2) * t564;
t278 = -t564 * qJ(4) - t473 + (-t443 + (-qJD(4) * t275 - t407) * qJD(1)) * t265;
t46 = -t192 * t267 + t278;
t8 = -qJD(5) * t45 - t270 * t46 + t274 * t41;
t276 = cos(qJ(1));
t419 = t275 * t276;
t272 = sin(qJ(1));
t425 = t271 * t272;
t199 = -t266 * t419 + t425;
t423 = t272 * t275;
t424 = t271 * t276;
t200 = t266 * t424 + t423;
t431 = t265 * t276;
t312 = -t200 * t270 + t274 * t431;
t519 = -t199 * t273 + t269 * t312;
t518 = t199 * t269 + t273 * t312;
t257 = t271 * t474;
t205 = pkin(8) * t432 + t257;
t162 = -t266 * qJ(3) - t205;
t137 = pkin(3) * t432 - t162;
t103 = pkin(4) * t432 - pkin(9) * t266 + t137;
t245 = qJ(4) * t432;
t416 = pkin(2) * t432 + qJ(3) * t435;
t384 = t245 + t416;
t116 = (pkin(9) * t271 - pkin(1)) * t265 - t384;
t445 = t270 * t103 + t274 * t116;
t237 = pkin(2) * t381;
t418 = qJ(4) * t381 + t237;
t88 = (-t407 + (-qJD(2) * t267 - qJD(4)) * t275) * t265 + t418;
t241 = t275 * t394;
t260 = t266 * qJD(3);
t417 = t241 + t260;
t99 = qJD(2) * t339 + t417;
t21 = -qJD(5) * t445 - t270 * t88 + t274 * t99;
t36 = pkin(10) * t211 + t45;
t48 = pkin(5) * t159 - pkin(10) * t160 + t84;
t13 = -t269 * t36 + t273 * t48;
t15 = -pkin(5) * t81 - pkin(10) * t80 - t42;
t7 = t270 * t41 + t274 * t46 + t85 * t404 - t405 * t97;
t5 = pkin(10) * t173 + t7;
t1 = qJD(6) * t13 + t15 * t269 + t273 * t5;
t14 = t269 * t48 + t273 * t36;
t2 = -qJD(6) * t14 + t15 * t273 - t269 * t5;
t516 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t187 = pkin(8) * t383 - t240;
t135 = t187 + t531;
t357 = mrSges(3,3) * t383;
t361 = mrSges(4,1) * t383;
t515 = -m(4) * t135 + t246 * t575 - t357 - t361;
t332 = t269 * mrSges(7,1) + t273 * mrSges(7,2);
t44 = -t270 * t97 + t274 * t85;
t35 = -pkin(5) * t211 - t44;
t98 = Ifges(7,4) * t104;
t40 = Ifges(7,1) * t105 + Ifges(7,5) * t158 + t98;
t447 = t273 * t40;
t457 = Ifges(7,4) * t105;
t39 = Ifges(7,2) * t104 + Ifges(7,6) * t158 + t457;
t497 = -t39 / 0.2e1;
t514 = t35 * t332 + t447 / 0.2e1 + t269 * t497;
t513 = t332 - t579;
t512 = -t13 * mrSges(7,1) + t14 * mrSges(7,2);
t511 = -m(7) * t366 - t274 * mrSges(7,3) + t526 - t563 + (m(5) + m(6)) * qJ(4);
t507 = t265 ^ 2;
t426 = t271 * t507;
t453 = Ifges(4,6) * t275;
t510 = (Ifges(3,1) * t275 - t461) * t426 / 0.2e1 - (-Ifges(4,2) * t275 + t454) * t426 / 0.2e1 + (-pkin(1) * (mrSges(3,1) * t271 + mrSges(3,2) * t275) - (t574 * t271 + t451 - t453) * t275 / 0.2e1) * t507;
t509 = t546 / 0.2e1 - t84 * mrSges(6,1) + t512;
t9 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t77;
t508 = Ifges(7,5) * t500 + Ifges(7,6) * t499 + Ifges(7,3) * t494 - t80 * Ifges(6,4) / 0.2e1 + t9 / 0.2e1 + (-t480 - t173 / 0.2e1) * Ifges(6,6) + (-t492 - t81 / 0.2e1) * Ifges(6,2);
t10 = t32 * Ifges(7,4) + t33 * Ifges(7,2) + t77 * Ifges(7,6);
t504 = t10 / 0.2e1;
t503 = Ifges(7,1) * t500 + Ifges(7,4) * t499 + Ifges(7,5) * t494;
t501 = Ifges(6,4) * t492 + t559;
t496 = t39 / 0.2e1;
t489 = -t104 / 0.2e1;
t488 = t104 / 0.2e1;
t487 = -t105 / 0.2e1;
t486 = t105 / 0.2e1;
t485 = -t158 / 0.2e1;
t484 = t158 / 0.2e1;
t483 = t159 / 0.2e1;
t482 = -t160 / 0.2e1;
t481 = t160 / 0.2e1;
t475 = pkin(1) * t265;
t471 = t270 * t7;
t6 = -pkin(5) * t173 - t8;
t470 = t274 * t6;
t465 = mrSges(6,3) * t159;
t463 = mrSges(7,3) * t269;
t462 = mrSges(7,3) * t273;
t460 = Ifges(3,4) * t275;
t459 = Ifges(6,4) * t270;
t458 = Ifges(6,4) * t274;
t456 = Ifges(7,4) * t269;
t455 = Ifges(7,4) * t273;
t452 = Ifges(5,6) * t271;
t446 = t274 * mrSges(6,3);
t442 = qJ(3) * t192;
t441 = qJ(3) * t199;
t201 = t266 * t423 + t424;
t440 = qJ(3) * t201;
t434 = t265 * t272;
t433 = t265 * t274;
t428 = t269 * t274;
t422 = t273 * t274;
t204 = -pkin(8) * t435 + t275 * t474;
t415 = t276 * pkin(1) + pkin(8) * t434;
t409 = qJD(2) * t275;
t403 = qJD(6) * t269;
t402 = qJD(6) * t273;
t396 = Ifges(6,5) * t80 + Ifges(6,6) * t81 + Ifges(6,3) * t173;
t202 = -t266 * t425 + t419;
t387 = t202 * pkin(2) + t415;
t164 = -t266 * pkin(2) - t204;
t380 = t265 * t409;
t372 = t414 / 0.2e1;
t367 = -pkin(1) * t272 + pkin(8) * t431;
t124 = t192 * mrSges(4,1) + t243 * mrSges(4,2);
t123 = mrSges(5,1) * t564 + t243 * mrSges(5,2);
t362 = t490 * t435;
t356 = mrSges(3,3) * t382;
t354 = t266 * qJ(4) - t164;
t353 = t274 * t382;
t348 = t275 * t372;
t344 = -t200 * pkin(2) + t367;
t163 = -t416 - t475;
t198 = t266 * t274 + t270 * t435;
t310 = -t266 * t270 + t271 * t433;
t337 = -mrSges(6,1) * t310 + mrSges(6,2) * t198;
t127 = -t198 * t269 + t265 * t421;
t128 = t198 * t273 + t269 * t432;
t334 = mrSges(7,1) * t127 - mrSges(7,2) * t128;
t331 = mrSges(4,2) * t275 - mrSges(4,3) * t271;
t330 = mrSges(5,2) * t271 + mrSges(5,3) * t275;
t329 = Ifges(6,1) * t270 + t458;
t328 = Ifges(7,1) * t273 - t456;
t327 = Ifges(7,1) * t269 + t455;
t326 = Ifges(6,2) * t274 + t459;
t325 = -Ifges(7,2) * t269 + t455;
t324 = Ifges(7,2) * t273 + t456;
t323 = Ifges(6,5) * t270 + Ifges(6,6) * t274;
t322 = Ifges(7,5) * t273 - Ifges(7,6) * t269;
t321 = Ifges(7,5) * t269 + Ifges(7,6) * t273;
t320 = t13 * t269 - t14 * t273;
t50 = pkin(10) * t432 + t445;
t102 = t435 * t491 + t354;
t65 = -pkin(5) * t310 - pkin(10) * t198 + t102;
t23 = t269 * t65 + t273 * t50;
t22 = -t269 * t50 + t273 * t65;
t66 = -mrSges(7,2) * t158 + mrSges(7,3) * t104;
t67 = mrSges(7,1) * t158 - mrSges(7,3) * t105;
t319 = t269 * t67 - t273 * t66;
t317 = t270 * t44 - t274 * t45;
t56 = t103 * t274 - t116 * t270;
t189 = -pkin(8) * t381 + t241;
t185 = -qJ(3) * t382 + t235;
t112 = -mrSges(6,2) * t211 - t465;
t313 = -t112 + t319;
t311 = t200 * t274 + t270 * t431;
t20 = t103 * t404 - t116 * t405 + t270 * t99 + t274 * t88;
t309 = t330 * t265;
t308 = t331 * t265;
t306 = (-t271 * Ifges(4,2) - t453) * t265;
t302 = pkin(3) * t434 + qJ(4) * t202 + t387;
t286 = pkin(4) * t434 + t302;
t285 = pkin(3) * t431 - qJ(4) * t200 + t344;
t283 = pkin(4) * t431 + t285;
t281 = t1 * t273 - t2 * t269 + (-t13 * t273 - t14 * t269) * qJD(6);
t280 = -qJD(5) * t317 + t274 * t8 + t471;
t259 = t266 * qJD(4);
t100 = t259 + (t275 * t363 - t257) * qJD(2);
t203 = (-mrSges(3,1) * t275 + mrSges(3,2) * t271) * t265;
t195 = t201 * pkin(2);
t193 = t199 * pkin(2);
t186 = qJD(1) * t309;
t184 = qJD(1) * t308;
t179 = -mrSges(3,2) * t246 + t356;
t157 = -t189 - t260;
t156 = Ifges(6,4) * t159;
t155 = t207 * t269 + t273 * t430;
t154 = t207 * t273 - t269 * t430;
t151 = t246 * t269 - t273 * t353;
t150 = t246 * t273 + t269 * t353;
t148 = t237 + (-qJ(3) * t409 - t407) * t265;
t146 = -qJD(1) * t362 + t240;
t145 = -t223 - t188;
t144 = t246 * Ifges(4,4) + qJD(1) * t306;
t138 = t185 + t224;
t136 = t163 - t245;
t130 = t202 * t270 + t272 * t433;
t129 = -t202 * t274 + t270 * t434;
t126 = qJD(5) * t198 - t274 * t380;
t125 = qJD(5) * t310 + t270 * t380;
t122 = -mrSges(4,1) * t564 - mrSges(4,3) * t243;
t120 = -t259 + (t432 * t490 + t257) * qJD(2);
t119 = -qJD(2) * t362 + t417;
t118 = pkin(3) * t435 - t354;
t110 = t147 + t530;
t95 = (-t407 + (-qJ(3) * qJD(2) - qJD(4)) * t275) * t265 + t418;
t92 = t383 * t490 - t240 - t522;
t91 = pkin(5) * t160 + pkin(10) * t159;
t87 = t130 * t273 - t201 * t269;
t86 = -t130 * t269 - t201 * t273;
t82 = -t473 - t442 + (-qJD(1) * t407 - t443) * t265;
t74 = t160 * Ifges(6,1) - t156 + t547;
t64 = qJD(6) * t127 + t125 * t273 - t269 * t381;
t63 = -qJD(6) * t128 - t125 * t269 - t273 * t381;
t62 = pkin(3) * t564 + t279;
t59 = -mrSges(6,2) * t173 + mrSges(6,3) * t81;
t54 = t278 - t442;
t49 = -pkin(5) * t432 - t56;
t47 = pkin(5) * t126 - pkin(10) * t125 + t100;
t25 = t269 * t91 + t273 * t44;
t24 = -t269 * t44 + t273 * t91;
t17 = pkin(5) * t381 - t21;
t16 = -pkin(10) * t381 + t20;
t4 = -qJD(6) * t23 - t16 * t269 + t273 * t47;
t3 = qJD(6) * t22 + t16 * t273 + t269 * t47;
t11 = [(t568 * t266 + (t271 * t570 + t275 * t578) * t265) * t243 / 0.2e1 + (-m(7) * (pkin(5) * t130 + t286) - t87 * mrSges(7,1) - t86 * mrSges(7,2) - m(6) * t286 - t130 * mrSges(6,1) - m(5) * (t302 + t440) - m(4) * (t387 + t440) - m(3) * t415 - mrSges(2,1) * t276 + t272 * mrSges(2,2) + t517 * t129 + t528 * t434 - t526 * t202 - t579 * t201) * g(2) + t523 * t266 / 0.2e1 + ((t271 * Ifges(3,1) + t460) * t265 + t304 + t552 * t266) * t192 / 0.2e1 + (t1 * t127 - t128 * t2 - t13 * t64 + t14 * t63) * mrSges(7,3) + (t125 * t84 + t381 * t45 - t432 * t7) * mrSges(6,2) + (t271 * t539 + t275 * t540) * t410 / 0.2e1 - (t275 * t144 + t521) * t410 / 0.2e1 - (-mrSges(6,3) * t7 - Ifges(6,4) * t493 + t508 + t516) * t310 + (Ifges(3,4) * t192 + Ifges(3,2) * t564 + Ifges(3,6) * t243 + t396) * t432 / 0.2e1 + t205 * (-mrSges(3,2) * t243 + mrSges(3,3) * t564) + t136 * (-mrSges(5,2) * t192 - mrSges(5,3) * t564) + t163 * (mrSges(4,2) * t564 - mrSges(4,3) * t192) + t564 * (Ifges(3,6) * t266 + t307) / 0.2e1 - (t553 * t243 + (-Ifges(4,6) + Ifges(5,6)) * t192 - t574 * t564) * t432 / 0.2e1 - ((-t275 * Ifges(5,2) + t452) * t265 + t305 + t553 * t266) * t564 / 0.2e1 + (t552 * t243 + (Ifges(3,1) + Ifges(5,3)) * t192 - (-Ifges(3,4) + Ifges(5,6)) * t564) * t435 / 0.2e1 - (-mrSges(3,1) * t564 + mrSges(3,2) * t192) * t475 - (Ifges(4,4) * t243 - Ifges(4,2) * t192 - Ifges(4,6) * t564) * t435 / 0.2e1 + t82 * t308 - t54 * t309 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t507 + t107 * t205 + t108 * t204 + t188 * t189) + m(4) * (t145 * t157 + t148 * t149 + t162 * t79 + t163 * t82 + t164 * t89) - t159 * (Ifges(6,4) * t125 - Ifges(6,6) * t381) / 0.2e1 + (Ifges(6,4) * t198 + Ifges(6,6) * t432) * t492 - t188 * mrSges(3,3) * t381 - t110 * mrSges(5,1) * t381 + t211 * (Ifges(6,5) * t125 - Ifges(6,3) * t381) / 0.2e1 + (Ifges(6,5) * t198 + Ifges(6,3) * t432) * t480 + t108 * (mrSges(3,1) * t266 - mrSges(3,3) * t435) + t51 * (mrSges(5,1) * t435 - mrSges(5,3) * t266) + t89 * (mrSges(4,1) * t435 + mrSges(4,2) * t266) + t79 * (-mrSges(4,1) * t432 - mrSges(4,3) * t266) + t107 * (-mrSges(3,2) * t266 + mrSges(3,3) * t432) + t8 * (mrSges(6,1) * t432 - mrSges(6,3) * t198) + t62 * (mrSges(5,1) * t432 + mrSges(5,2) * t266) + t22 * t18 + t23 * t19 - t192 * (Ifges(4,4) * t266 + t306) / 0.2e1 + m(6) * (t100 * t84 - t102 * t42 + t20 * t45 + t21 * t44 + t445 * t7 + t56 * t8) + t445 * t59 + (Ifges(7,1) * t64 + Ifges(7,4) * t63) * t486 + (Ifges(7,1) * t128 + Ifges(7,4) * t127) * t500 - pkin(1) * t203 * t398 - t6 * t334 - t42 * t337 + (Ifges(6,1) * t125 - Ifges(6,5) * t381) * t481 + (Ifges(6,1) * t198 + Ifges(6,5) * t432) * t493 + t44 * (-mrSges(6,1) * t381 - mrSges(6,3) * t125) + t145 * mrSges(4,1) * t381 + ((m(3) * t187 - t515) * t205 - t561) * qJD(2) + (Ifges(7,4) * t64 + Ifges(7,2) * t63) * t488 + (Ifges(7,4) * t128 + Ifges(7,2) * t127) * t499 + (t548 / 0.2e1 - Ifges(6,4) * t481 + Ifges(7,3) * t484 + Ifges(7,5) * t486 + Ifges(7,6) * t488 - t509 + t560) * t126 + Ifges(2,3) * qJDD(1) + (-m(6) * t283 - t312 * mrSges(6,1) - m(7) * (pkin(5) * t312 + t283) - t518 * mrSges(7,1) + t519 * mrSges(7,2) - m(5) * (t285 - t441) - m(4) * (t344 - t441) - m(3) * t367 + t272 * mrSges(2,1) + mrSges(2,2) * t276 + t517 * t311 + t528 * t431 + t526 * t200 + t579 * t199) * g(1) + ((t275 * (-Ifges(3,2) * t271 + t460) + t271 * (Ifges(5,3) * t275 + t452)) * t507 / 0.2e1 + t510) * t399 + m(5) * (t110 * t119 + t117 * t95 + t118 * t51 + t120 * t92 + t136 * t54 + t137 * t62) + m(7) * (t1 * t23 + t13 * t4 + t14 * t3 + t17 * t35 + t2 * t22 + t49 * t6) + t204 * (mrSges(3,1) * t243 - mrSges(3,3) * t192) + (mrSges(4,1) * t135 + mrSges(5,1) * t92 + mrSges(3,3) * t187) * t380 + t49 * t12 + t17 * t55 + t56 * t58 + t64 * t40 / 0.2e1 + t35 * (-mrSges(7,1) * t63 + mrSges(7,2) * t64) + t3 * t66 + t4 * t67 + (Ifges(7,5) * t64 + Ifges(7,6) * t63) * t484 + (Ifges(7,5) * t128 + Ifges(7,6) * t127) * t494 + t100 * t90 + t102 * t34 + t20 * t112 + t21 * t113 + t118 * t121 + t125 * t74 / 0.2e1 + t63 * t496 + t198 * t501 + t128 * t503 + t127 * t504 + t137 * t123 + t162 * t122 + t164 * t124 + t120 * t180 + t157 * t181 + t119 * t182 + t148 * t184 - t95 * t186 + t189 * t179; t577 * t112 + (qJ(3) * t62 - t117 * t138 - t268 * t51 + t576 * t92 + (qJD(3) - t146) * t110) * m(5) + (-m(4) * t145 + t179 - t181 - t356) * t187 + t551 * t67 - t536 * qJD(3) + (-t1 * t270 - t14 * t353 - t35 * t537) * mrSges(7,2) + (t13 * t353 + t2 * t270 + t35 * t538) * mrSges(7,1) + (t521 + t159 * (-Ifges(6,6) * t271 + t275 * t326)) * t372 + (-t321 * t484 - t324 * t488 - t327 * t486) * t401 + t549 * t267 * t274 + (-t44 * t60 - t45 * t61 - t268 * t42 + (t270 * t45 + t274 * t44) * qJD(3) + t280 * t267 + (qJD(4) - t115) * t84) * m(6) + (t405 * t44 - t471) * mrSges(6,3) + t550 * t66 + (-t122 + t123) * qJ(3) - ((Ifges(5,3) * t382 + t232 + t539) * t271 + (-Ifges(3,2) * t383 + t270 * t74 + t274 * t73 + t233 + t540) * t275 + t160 * (-Ifges(6,5) * t271 + t275 * t329) + t211 * (-Ifges(6,3) * t271 + t275 * t323)) * t414 / 0.2e1 - (t160 * t329 + t211 * t323) * qJD(5) / 0.2e1 + t382 * t545 + (-t45 * (mrSges(6,2) * t271 + t275 * t446) - t44 * (-mrSges(6,1) * t271 - mrSges(6,3) * t427)) * t414 + t458 * t492 + t59 * t430 + t571 * t113 + (t357 + t515) * t188 - t8 * t446 + t565 * t55 + (t1 * t155 + t551 * t13 + t550 * t14 + t154 * t2 - t470 * t267 + t35 * t565) * m(7) + t542 * qJD(4) + t543 * t268 + t523 - t459 * t493 - (t269 * t40 + t273 * t39) * t401 / 0.2e1 - (t447 + t74) * t405 / 0.2e1 - t10 * t428 / 0.2e1 + t144 * t348 - t42 * t335 + (t545 + t326 * t483 + (Ifges(7,3) * t274 - t270 * t322) * t484 + (Ifges(7,5) * t274 - t270 * t328) * t486 + (Ifges(7,6) * t274 - t270 * t325) * t488) * qJD(5) + t508 * t270 + (-pkin(2) * t89 - qJ(3) * t79 - qJD(3) * t145 - t149 * t185) * m(4) + (t558 * (qJ(3) * t202 - t195) + t557 * t195 + t513 * t202 + t511 * t201) * g(1) + (t558 * (qJ(3) * t200 - t193) + t557 * t193 + t513 * t200 + t511 * t199) * g(2) + (-m(4) * t416 + t203 + t400 * t384 + (-t330 + t331 + (t274 * t392 + t563) * t275 + (pkin(9) * t557 + mrSges(6,3) + t332) * t271) * t265) * g(3) + (-t510 * qJD(1) + t561) * qJD(1) + (-t1 * t428 + t13 * t537 - t14 * t538 - t2 * t422) * mrSges(7,3) + (-t512 + t560) * t404 + (t322 * t494 + t325 * t499 + t328 * t500 + t348 * t38 + t501 + t559) * t274 - t145 * t361 - t135 * t360 - t92 * t358 + t110 * t359 - t51 * mrSges(5,3) + t62 * mrSges(5,2) - t79 * mrSges(4,3) + t89 * mrSges(4,2) - t107 * mrSges(3,2) + t108 * mrSges(3,1) - t115 * t90 - pkin(2) * t124 + t332 * t470 + (Ifges(7,5) * t153 + Ifges(7,6) * t152 - Ifges(7,3) * t353) * t485 + (Ifges(7,1) * t153 + Ifges(7,4) * t152 - Ifges(7,5) * t353) * t487 + (Ifges(7,4) * t153 + Ifges(7,2) * t152 - Ifges(7,6) * t353) * t489 + t377 * t496 + t152 * t497 + t422 * t503 - t153 * t40 / 0.2e1 + t154 * t18 + t155 * t19 - t147 * t180 - t146 * t182 - t185 * t184 + t138 * t186; (t184 - t186) * t383 + t319 * qJD(6) + t313 * t159 + t536 * t246 + t444 * t160 - t273 * t18 - t269 * t19 + t124 + (-g(1) * t201 - g(2) * t199 + g(3) * t432) * (m(4) - t400) + (-t1 * t269 + t158 * t320 + t160 * t35 - t2 * t273) * m(7) + (-t159 * t45 - t160 * t44 + t42) * m(6) + (-t110 * t246 + t117 * t383 + t51) * m(5) + (t145 * t246 + t149 * t383 + t89) * m(4) - t543; -t150 * t67 - t151 * t66 - t542 * t246 - t186 * t382 + (-qJD(5) * t313 + t112 * t382 + t549) * t274 + (t59 + (-t269 * t66 - t273 * t67) * qJD(6) + t211 * t444 + t532) * t270 + t123 + (g(1) * t202 + g(2) * t200 + g(3) * t435) * t400 + ((-qJD(5) * t320 - t6) * t274 - t13 * t150 - t14 * t151 + (t211 * t35 + t281) * t270) * m(7) + (-t246 * t84 - t317 * t382 + t280) * m(6) + (t117 * t382 + t246 * t92 + t62) * m(5); (m(7) * t281 - t402 * t67 - t403 * t66 + t532) * pkin(10) + (-t198 * t392 - t300 * t310 + t337) * g(3) + (-t13 * t402 - t14 * t403) * mrSges(7,3) + t514 * qJD(6) + (-t156 + t74) * t483 + (Ifges(7,5) * t487 - Ifges(6,2) * t483 + Ifges(7,6) * t489 + Ifges(7,3) * t485 + t509) * t160 + (-m(7) * t35 - t444 + t464) * t45 + (-t465 - t112) * t44 - pkin(5) * t12 - t7 * mrSges(6,2) + t8 * mrSges(6,1) + t396 + (t104 * t325 + t105 * t328 + t158 * t322) * qJD(6) / 0.2e1 + t6 * t333 + (t129 * t544 + t130 * t517) * g(1) + (-t311 * t544 - t312 * t517) * g(2) + (t547 / 0.2e1 - Ifges(6,1) * t482 - t322 * t485 - t328 * t487 - t325 * t489 + t84 * mrSges(6,2) - t13 * t462 - t14 * t463 + t514) * t159 + (-pkin(5) * t6 - t13 * t24 - t14 * t25) * m(7) + t1 * t462 - t25 * t66 - t24 * t67 + (-t450 + t38) * t482 + t73 * t481 + t321 * t494 + t324 * t499 + t327 * t500 + t269 * t503 + t273 * t504 - t2 * t463; -t35 * (mrSges(7,1) * t105 + mrSges(7,2) * t104) + (Ifges(7,1) * t104 - t457) * t487 + t39 * t486 + (Ifges(7,5) * t104 - Ifges(7,6) * t105) * t485 - t13 * t66 + t14 * t67 - g(1) * (mrSges(7,1) * t86 - mrSges(7,2) * t87) - g(2) * (mrSges(7,1) * t519 + mrSges(7,2) * t518) - g(3) * t334 + (t104 * t13 + t105 * t14) * mrSges(7,3) + t9 + (-Ifges(7,2) * t105 + t40 + t98) * t489 + t516;];
tau  = t11;
