% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:56
% EndTime: 2019-03-09 18:11:41
% DurationCPUTime: 26.93s
% Computational Cost: add. (16035->788), mult. (34657->1007), div. (0->0), fcn. (24691->12), ass. (0->362)
t598 = mrSges(4,1) + mrSges(5,1);
t345 = qJD(2) + qJD(3);
t339 = qJD(5) - t345;
t349 = sin(qJ(6));
t354 = cos(qJ(6));
t351 = sin(qJ(3));
t356 = cos(qJ(2));
t520 = cos(qJ(3));
t422 = t520 * t356;
t352 = sin(qJ(2));
t449 = qJD(1) * t352;
t252 = -qJD(1) * t422 + t351 * t449;
t272 = t351 * t356 + t352 * t520;
t253 = t272 * qJD(1);
t350 = sin(qJ(5));
t355 = cos(qJ(5));
t550 = t252 * t350 + t355 * t253;
t138 = t339 * t354 - t349 * t550;
t139 = t339 * t349 + t354 * t550;
t498 = mrSges(6,3) * t550;
t552 = mrSges(6,1) * t339 + mrSges(7,1) * t138 - mrSges(7,2) * t139 - t498;
t244 = t253 * pkin(9);
t358 = -pkin(8) - pkin(7);
t304 = t358 * t352;
t282 = qJD(1) * t304;
t263 = qJD(2) * pkin(2) + t282;
t305 = t358 * t356;
t283 = qJD(1) * t305;
t456 = t351 * t283;
t180 = t520 * t263 + t456;
t421 = qJD(4) - t180;
t540 = pkin(3) + pkin(4);
t112 = -t345 * t540 - t244 + t421;
t338 = t345 * qJ(4);
t423 = t520 * t283;
t181 = t351 * t263 - t423;
t513 = pkin(9) * t252;
t407 = t181 + t513;
t122 = t338 + t407;
t74 = t112 * t355 - t122 * t350;
t69 = -pkin(5) * t339 - t74;
t593 = -m(7) * t69 + t552;
t599 = -m(6) * t74 - t593;
t554 = m(6) + m(7);
t581 = Ifges(4,1) + Ifges(5,1);
t580 = Ifges(5,4) + Ifges(4,5);
t579 = Ifges(5,5) - Ifges(4,4);
t346 = qJ(2) + qJ(3);
t340 = sin(t346);
t341 = cos(t346);
t596 = -t598 * t341 + (mrSges(4,2) - mrSges(5,3)) * t340;
t190 = t282 * t351 - t423;
t446 = qJD(3) * t351;
t434 = pkin(2) * t446;
t595 = t434 - t190;
t442 = qJD(1) * qJD(2);
t286 = qJDD(1) * t356 - t352 * t442;
t287 = qJDD(1) * t352 + t356 * t442;
t370 = t272 * qJD(3);
t149 = qJD(1) * t370 - t286 * t520 + t351 * t287;
t378 = -t351 * t352 + t422;
t369 = t378 * qJD(3);
t148 = qJD(1) * t369 + t351 * t286 + t287 * t520;
t347 = qJDD(1) * pkin(1);
t240 = -t286 * pkin(2) - t347;
t68 = t149 * pkin(3) - t148 * qJ(4) - t253 * qJD(4) + t240;
t44 = -pkin(4) * t149 - t68;
t171 = t252 * t355 - t350 * t253;
t64 = qJD(5) * t171 + t148 * t355 + t149 * t350;
t65 = -qJD(5) * t550 - t148 * t350 + t149 * t355;
t11 = -pkin(5) * t65 - pkin(10) * t64 + t44;
t344 = qJDD(2) + qJDD(3);
t269 = t287 * pkin(7);
t196 = qJDD(2) * pkin(2) - pkin(8) * t287 - t269;
t268 = t286 * pkin(7);
t199 = pkin(8) * t286 + t268;
t420 = qJD(3) * t520;
t92 = t196 * t520 - t351 * t199 - t263 * t446 + t283 * t420;
t368 = qJDD(4) - t92;
t51 = -t148 * pkin(9) - t344 * t540 + t368;
t91 = t351 * t196 + t520 * t199 + t263 * t420 + t283 * t446;
t81 = t344 * qJ(4) + t345 * qJD(4) + t91;
t52 = pkin(9) * t149 + t81;
t15 = qJD(5) * t74 + t350 * t51 + t355 * t52;
t337 = qJDD(5) - t344;
t12 = pkin(10) * t337 + t15;
t448 = qJD(1) * t356;
t303 = -qJD(1) * pkin(1) - pkin(2) * t448;
t153 = t252 * pkin(3) - t253 * qJ(4) + t303;
t113 = -pkin(4) * t252 - t153;
t66 = -pkin(5) * t171 - pkin(10) * t550 + t113;
t75 = t112 * t350 + t122 * t355;
t70 = pkin(10) * t339 + t75;
t23 = -t349 * t70 + t354 * t66;
t2 = qJD(6) * t23 + t11 * t349 + t12 * t354;
t24 = t349 * t66 + t354 * t70;
t3 = -qJD(6) * t24 + t11 * t354 - t12 * t349;
t365 = t2 * t354 - t3 * t349 + (-t23 * t354 - t24 * t349) * qJD(6);
t444 = qJD(6) * t354;
t445 = qJD(6) * t349;
t544 = qJD(6) * t138;
t39 = t337 * t349 + t354 * t64 + t544;
t63 = qJDD(6) - t65;
t18 = mrSges(7,1) * t63 - mrSges(7,3) * t39;
t543 = qJD(6) * t139;
t40 = t337 * t354 - t349 * t64 - t543;
t19 = -mrSges(7,2) * t63 + mrSges(7,3) * t40;
t545 = -t349 * t18 + t354 * t19;
t168 = qJD(6) - t171;
t93 = -mrSges(7,2) * t168 + mrSges(7,3) * t138;
t94 = mrSges(7,1) * t168 - mrSges(7,3) * t139;
t594 = m(7) * t365 - t94 * t444 - t93 * t445 + t545;
t526 = t272 / 0.2e1;
t591 = t286 / 0.2e1;
t590 = t349 / 0.2e1;
t521 = t356 / 0.2e1;
t397 = Ifges(7,5) * t354 - Ifges(7,6) * t349;
t589 = t397 / 0.2e1;
t492 = Ifges(7,4) * t354;
t399 = -Ifges(7,2) * t349 + t492;
t588 = t399 / 0.2e1;
t493 = Ifges(7,4) * t349;
t401 = Ifges(7,1) * t354 - t493;
t587 = t401 / 0.2e1;
t586 = Ifges(7,5) * t39;
t585 = Ifges(7,6) * t40;
t584 = Ifges(7,3) * t63;
t583 = t23 * mrSges(7,1);
t582 = t24 * mrSges(7,2);
t578 = Ifges(4,6) - Ifges(5,6);
t17 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t54 = mrSges(6,1) * t337 - mrSges(6,3) * t64;
t577 = t54 - t17;
t576 = t356 * Ifges(3,2);
t242 = Ifges(4,4) * t252;
t486 = t252 * Ifges(5,5);
t574 = t253 * t581 + t345 * t580 - t242 + t486;
t501 = mrSges(4,3) * t252;
t212 = -mrSges(4,2) * t345 - t501;
t503 = mrSges(5,2) * t252;
t215 = mrSges(5,3) * t345 - t503;
t573 = -t212 - t215;
t500 = mrSges(4,3) * t253;
t572 = -mrSges(5,2) * t253 + t345 * t598 - t500;
t514 = pkin(7) * t356;
t515 = pkin(7) * t352;
t571 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t449) * t514 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t448) * t515;
t570 = t268 * t356 + t269 * t352;
t297 = -mrSges(7,1) * t354 + mrSges(7,2) * t349;
t353 = sin(qJ(1));
t357 = cos(qJ(1));
t569 = g(1) * t357 + g(2) * t353;
t432 = m(7) * pkin(10) + mrSges(7,3);
t499 = mrSges(6,3) * t171;
t151 = -mrSges(6,2) * t339 + t499;
t382 = -t349 * t94 + t354 * t93 + t151;
t460 = t341 * t353;
t461 = t340 * t355;
t216 = t350 * t460 - t353 * t461;
t247 = t340 * t350 + t341 * t355;
t217 = t247 * t353;
t454 = t216 * mrSges(6,1) + t217 * mrSges(6,2);
t567 = -mrSges(5,3) * t460 + t217 * mrSges(7,3) + t216 * t297 - t454;
t455 = t355 * t357;
t457 = t350 * t357;
t218 = -t340 * t457 - t341 * t455;
t219 = -t340 * t455 + t341 * t457;
t453 = t219 * mrSges(6,1) - t218 * mrSges(6,2);
t459 = t341 * t357;
t566 = -mrSges(5,3) * t459 - t218 * mrSges(7,3) + t219 * t297 - t453;
t565 = -mrSges(6,2) + t432;
t298 = -mrSges(3,1) * t356 + mrSges(3,2) * t352;
t564 = -m(3) * pkin(1) - mrSges(2,1) + t298 + t596;
t248 = -t341 * t350 + t461;
t452 = t247 * mrSges(6,1) + t248 * mrSges(6,2);
t563 = t248 * mrSges(7,3) + t247 * t297 - t452 + t596;
t562 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t108 = pkin(5) * t550 - pkin(10) * t171;
t560 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3) + mrSges(6,3) - t554 * (-pkin(9) - t358);
t525 = -t339 / 0.2e1;
t534 = -t171 / 0.2e1;
t535 = -t168 / 0.2e1;
t538 = -t139 / 0.2e1;
t539 = -t138 / 0.2e1;
t559 = -t113 * mrSges(6,1) + Ifges(7,5) * t538 - Ifges(6,2) * t534 - Ifges(6,6) * t525 + Ifges(7,6) * t539 + Ifges(7,3) * t535 + t582 - t583;
t491 = t139 * Ifges(7,4);
t60 = t138 * Ifges(7,2) + t168 * Ifges(7,6) + t491;
t429 = t60 * t590;
t133 = Ifges(7,4) * t138;
t61 = Ifges(7,1) * t139 + Ifges(7,5) * t168 + t133;
t479 = t354 * t61;
t481 = t354 * mrSges(7,3);
t497 = mrSges(7,3) * t349;
t533 = -t550 / 0.2e1;
t402 = t349 * mrSges(7,1) + t354 * mrSges(7,2);
t553 = t402 * t69;
t558 = -t113 * mrSges(6,2) + Ifges(6,1) * t533 + Ifges(6,5) * t525 + t23 * t481 + t24 * t497 + t397 * t535 + t399 * t539 + t401 * t538 - t479 / 0.2e1 + t429 - t553;
t59 = Ifges(7,5) * t139 + Ifges(7,6) * t138 + Ifges(7,3) * t168;
t557 = t59 / 0.2e1;
t494 = Ifges(6,4) * t550;
t98 = Ifges(6,2) * t171 + Ifges(6,6) * t339 + t494;
t556 = -t98 / 0.2e1;
t165 = Ifges(6,4) * t171;
t99 = Ifges(6,1) * t550 + t339 * Ifges(6,5) + t165;
t555 = t99 / 0.2e1;
t137 = t244 + t180;
t289 = t355 * qJ(4) - t350 * t540;
t551 = qJD(5) * t289 + t355 * t407 + (qJD(4) - t137) * t350;
t203 = t351 * t304 - t520 * t305;
t437 = t520 * pkin(2);
t327 = -t437 - pkin(3);
t318 = -pkin(4) + t327;
t518 = pkin(2) * t351;
t323 = qJ(4) + t518;
t224 = t350 * t318 + t355 * t323;
t306 = qJ(4) * t460;
t517 = pkin(2) * t352;
t384 = -t340 * t540 - t517;
t549 = t353 * t384 + t306;
t319 = t340 * qJ(4);
t548 = pkin(3) * t459 + t357 * t319;
t308 = qJ(4) * t459;
t547 = t357 * t384 + t308;
t396 = -t23 * t349 + t24 * t354;
t16 = -qJD(5) * t75 - t350 * t52 + t355 * t51;
t13 = -pkin(5) * t337 - t16;
t8 = Ifges(7,4) * t39 + Ifges(7,2) * t40 + Ifges(7,6) * t63;
t9 = t39 * Ifges(7,1) + t40 * Ifges(7,4) + t63 * Ifges(7,5);
t542 = t13 * t297 + t16 * mrSges(6,1) - t15 * mrSges(6,2) + Ifges(6,5) * t64 + Ifges(6,6) * t65 + Ifges(6,3) * t337 + t2 * t481 - t3 * t497 + (-t23 * t444 - t24 * t445) * mrSges(7,3) + t543 * t587 + t544 * t588 + t63 * (Ifges(7,5) * t349 + Ifges(7,6) * t354) / 0.2e1 + t40 * (Ifges(7,2) * t354 + t493) / 0.2e1 + t39 * (Ifges(7,1) * t349 + t492) / 0.2e1 + t354 * t8 / 0.2e1 + t9 * t590 + (t168 * t589 + t553) * qJD(6);
t537 = t139 / 0.2e1;
t532 = t550 / 0.2e1;
t530 = -t252 / 0.2e1;
t529 = t252 / 0.2e1;
t527 = t253 / 0.2e1;
t522 = t345 / 0.2e1;
t519 = mrSges(6,3) * t75;
t516 = pkin(4) * t253;
t325 = t341 * pkin(3);
t343 = t356 * pkin(2);
t508 = t74 * mrSges(6,3);
t328 = t343 + pkin(1);
t504 = mrSges(4,2) * t341;
t496 = Ifges(3,4) * t352;
t495 = Ifges(3,4) * t356;
t169 = t338 + t181;
t490 = t169 * mrSges(5,2);
t489 = t181 * mrSges(4,3);
t485 = t253 * Ifges(4,4);
t192 = qJD(2) * t378 + t369;
t193 = qJD(2) * t272 + t370;
t388 = -t272 * t350 - t355 * t378;
t88 = qJD(5) * t388 + t192 * t355 + t193 * t350;
t483 = t349 * t88;
t185 = t272 * t355 - t350 * t378;
t467 = t185 * t349;
t466 = t185 * t354;
t458 = t345 * t355;
t174 = t253 * pkin(3) + t252 * qJ(4);
t191 = t520 * t282 + t456;
t450 = t325 + t319;
t447 = qJD(2) * t352;
t443 = m(5) + t554;
t440 = t584 + t585 + t586;
t436 = pkin(2) * t449;
t435 = pkin(2) * t447;
t427 = t479 / 0.2e1;
t324 = t341 * pkin(4);
t426 = t324 + t450;
t425 = t343 + t450;
t424 = qJD(2) * t358;
t417 = -t445 / 0.2e1;
t416 = -t444 / 0.2e1;
t415 = t216 * pkin(5) - pkin(10) * t217;
t414 = t219 * pkin(5) + pkin(10) * t218;
t413 = t442 / 0.2e1;
t124 = -t344 * mrSges(5,1) + t148 * mrSges(5,2);
t411 = -t328 - t319;
t202 = -t520 * t304 - t305 * t351;
t309 = t357 * t328;
t410 = -t353 * t358 + t309;
t179 = -pkin(3) * t378 - t272 * qJ(4) - t328;
t409 = pkin(2) * t420;
t408 = pkin(4) * t459 + t309 + t548;
t134 = -t174 - t516;
t405 = mrSges(3,1) * t352 + mrSges(3,2) * t356;
t400 = t496 + t576;
t398 = Ifges(3,5) * t356 - Ifges(3,6) * t352;
t166 = -pkin(9) * t272 + t202;
t167 = -pkin(9) * t378 + t203;
t101 = t166 * t350 + t167 * t355;
t147 = pkin(4) * t378 - t179;
t78 = -pkin(5) * t388 - pkin(10) * t185 + t147;
t38 = t101 * t354 + t349 * t78;
t37 = -t101 * t349 + t354 * t78;
t288 = -t350 * qJ(4) - t355 * t540;
t392 = t355 * t166 - t167 * t350;
t391 = -t217 * t354 - t349 * t357;
t390 = t217 * t349 - t354 * t357;
t223 = t318 * t355 - t323 * t350;
t158 = t436 + t174;
t383 = t190 + t513;
t381 = pkin(1) * t405;
t380 = t185 * t444 + t483;
t379 = t185 * t445 - t354 * t88;
t377 = t247 * pkin(5) - pkin(10) * t248 + t426;
t376 = t352 * (Ifges(3,1) * t356 - t496);
t103 = t193 * pkin(3) - t192 * qJ(4) - t272 * qJD(4) + t435;
t284 = t352 * t424;
t285 = t356 * t424;
t114 = t520 * t284 + t351 * t285 + t304 * t420 + t305 * t446;
t374 = -m(7) * pkin(5) + t297;
t367 = t353 * (-t341 * t540 + t411);
t79 = -pkin(4) * t193 - t103;
t366 = m(5) * (-pkin(3) * t340 - t517) - t340 * mrSges(5,1);
t71 = -t108 + t134;
t364 = (-t349 * t93 - t354 * t94) * qJD(6) + t545;
t115 = qJD(3) * t203 + t351 * t284 - t520 * t285;
t363 = -t192 * pkin(9) + t115;
t362 = t504 + (m(5) * pkin(3) + t540 * t554 + t598) * t340;
t157 = -pkin(3) * t345 + t421;
t241 = Ifges(5,5) * t253;
t161 = t345 * Ifges(5,6) + t252 * Ifges(5,3) + t241;
t162 = -t252 * Ifges(4,2) + t345 * Ifges(4,6) + t485;
t84 = -t344 * pkin(3) + t368;
t360 = -t153 * (mrSges(5,1) * t253 + mrSges(5,3) * t252) - t303 * (mrSges(4,1) * t253 - mrSges(4,2) * t252) + qJD(6) * t429 + t61 * t416 - t542 - t91 * mrSges(4,2) + t92 * mrSges(4,1) - t84 * mrSges(5,1) + t81 * mrSges(5,3) + t253 * t490 - t180 * t501 + t157 * t503 + t162 * t527 + (Ifges(5,3) * t253 - t486) * t530 - (-t252 * t580 - t253 * t578) * t345 / 0.2e1 + (Ifges(5,2) + Ifges(4,3)) * t344 - t578 * t149 + t580 * t148 + (-Ifges(4,2) * t253 - t242 + t574) * t529 - (-t581 * t252 + t161 + t241 - t485) * t253 / 0.2e1 + (Ifges(6,4) * t533 - t519 + t556 + t557 - t559) * t550 + (-Ifges(6,4) * t534 - t508 + t555 - t558) * t171;
t330 = Ifges(3,4) * t448;
t316 = t409 + qJD(4);
t281 = -pkin(10) + t289;
t280 = pkin(5) - t288;
t251 = Ifges(3,1) * t449 + Ifges(3,5) * qJD(2) + t330;
t250 = Ifges(3,6) * qJD(2) + qJD(1) * t400;
t229 = t355 * qJD(4) + qJD(5) * t288;
t220 = pkin(5) - t223;
t198 = t253 * t349 + t354 * t458;
t197 = t253 * t354 - t349 * t458;
t188 = -t218 * t354 - t349 * t353;
t187 = t218 * t349 - t353 * t354;
t178 = mrSges(4,1) * t252 + mrSges(4,2) * t253;
t177 = mrSges(5,1) * t252 - mrSges(5,3) * t253;
t150 = t244 + t191;
t126 = -mrSges(5,2) * t149 + mrSges(5,3) * t344;
t125 = -mrSges(4,2) * t344 - mrSges(4,3) * t149;
t123 = mrSges(4,1) * t344 - mrSges(4,3) * t148;
t120 = -t158 - t516;
t107 = -mrSges(6,1) * t171 + mrSges(6,2) * t550;
t96 = pkin(9) * t193 + t114;
t89 = qJD(5) * t185 + t192 * t350 - t355 * t193;
t87 = t355 * t150 + t350 * t383;
t83 = t355 * t137 + t350 * t407;
t67 = t71 - t436;
t55 = -mrSges(6,2) * t337 + mrSges(6,3) * t65;
t34 = t108 * t349 + t354 * t74;
t33 = t108 * t354 - t349 * t74;
t30 = t349 * t67 + t354 * t87;
t29 = -t349 * t87 + t354 * t67;
t28 = t349 * t71 + t354 * t83;
t27 = -t349 * t83 + t354 * t71;
t25 = qJD(5) * t392 + t350 * t363 + t355 * t96;
t20 = pkin(5) * t89 - pkin(10) * t88 + t79;
t5 = -qJD(6) * t38 + t20 * t354 - t25 * t349;
t4 = qJD(6) * t37 + t20 * t349 + t25 * t354;
t1 = [(m(4) * t91 + m(5) * t81 + t125 + t126) * t203 + (-m(4) * t92 + m(5) * t84 - t123 + t124) * t202 + t599 * (qJD(5) * t101 + t350 * t96 - t355 * t363) + (Ifges(3,4) * t287 + Ifges(3,2) * t286) * t521 + t287 * t495 / 0.2e1 + (-t2 * t467 + t23 * t379 - t24 * t380 - t3 * t466) * mrSges(7,3) + t9 * t466 / 0.2e1 - t8 * t467 / 0.2e1 - t381 * t442 - t250 * t447 / 0.2e1 - t89 * t519 + t88 * t427 + m(5) * (t103 * t153 + t179 * t68) + m(4) * (-t240 * t328 + t303 * t435) + (t356 * t495 + t376) * t413 + t400 * t591 + t178 * t435 - t179 * mrSges(5,3) * t148 + (t574 / 0.2e1 + mrSges(4,2) * t303 + mrSges(5,2) * t157 - mrSges(4,3) * t180 - mrSges(5,3) * t153 + Ifges(4,4) * t530 + Ifges(5,5) * t529 + t522 * t580 + t527 * t581) * t192 + t89 * t583 - t328 * mrSges(4,2) * t148 + m(7) * (t2 * t38 + t23 * t5 + t24 * t4 + t3 * t37) + m(6) * (t101 * t15 + t113 * t79 + t147 * t44 + t25 * t75) + t240 * (-mrSges(4,1) * t378 + mrSges(4,2) * t272) + t68 * (-mrSges(5,1) * t378 - mrSges(5,3) * t272) - t378 * (Ifges(5,5) * t148 + Ifges(5,6) * t344) / 0.2e1 + t378 * (Ifges(4,4) * t148 + Ifges(4,6) * t344) / 0.2e1 + (-t272 * t92 + t378 * t91) * mrSges(4,3) + (t272 * t84 + t378 * t81) * mrSges(5,2) + Ifges(2,3) * qJDD(1) + t339 * (Ifges(6,5) * t88 - Ifges(6,6) * t89) / 0.2e1 + t89 * t556 + t89 * t557 - pkin(1) * (-mrSges(3,1) * t286 + mrSges(3,2) * t287) + t103 * t177 + t171 * (Ifges(6,4) * t88 - Ifges(6,2) * t89) / 0.2e1 + t25 * t151 + t147 * (-mrSges(6,1) * t65 + mrSges(6,2) * t64) + t113 * (mrSges(6,1) * t89 + mrSges(6,2) * t88) + t79 * t107 + t101 * t55 + t4 * t93 + t5 * t94 + t69 * (mrSges(7,1) * t380 - mrSges(7,2) * t379) + t168 * (-Ifges(7,5) * t379 - Ifges(7,6) * t380 + Ifges(7,3) * t89) / 0.2e1 + t138 * (-Ifges(7,4) * t379 - Ifges(7,2) * t380 + Ifges(7,6) * t89) / 0.2e1 + t37 * t18 + t38 * t19 + (t44 * mrSges(6,2) - t16 * mrSges(6,3) + Ifges(6,1) * t64 + Ifges(6,4) * t65 + Ifges(6,5) * t337 + t13 * t402 + t39 * t587 + t40 * t588 + t60 * t416 + t61 * t417 + t589 * t63) * t185 + (-mrSges(3,1) * t515 - mrSges(3,2) * t514 + 0.2e1 * Ifges(3,6) * t521) * qJDD(2) + (Ifges(3,1) * t287 + Ifges(3,4) * t591 + Ifges(3,5) * qJDD(2) - t413 * t576) * t352 - t60 * t483 / 0.2e1 + (-m(6) * t367 + t217 * mrSges(6,1) - t391 * mrSges(7,1) - t390 * mrSges(7,2) - (-pkin(5) * t217 + t367) * m(7) + t565 * t216 + (m(4) * t328 - m(5) * (t411 - t325) - t564) * t353 + ((m(4) + m(5)) * t358 + t560) * t357) * g(1) + (-m(4) * t410 - m(5) * (t410 + t548) - m(6) * t408 + t218 * mrSges(6,1) - m(7) * (-pkin(5) * t218 + t408) - t188 * mrSges(7,1) - t187 * mrSges(7,2) - t565 * t219 + t564 * t357 + t560 * t353) * g(2) + (t286 * t514 + t287 * t515 + t570) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t570) + (t251 * t521 + t398 * qJD(2) / 0.2e1 - t571) * qJD(2) + (-m(4) * t180 + m(5) * t157 - t572) * t115 + (m(4) * t181 + m(5) * t169 - t573) * t114 - t298 * t347 + (Ifges(6,1) * t88 - Ifges(6,4) * t89) * t532 + (-Ifges(7,1) * t379 - Ifges(7,4) * t380 + Ifges(7,5) * t89) * t537 + t88 * t555 + (-t328 * mrSges(4,1) + t179 * mrSges(5,1) + (-Ifges(4,2) - Ifges(5,3)) * t378 + 0.2e1 * t579 * t526) * t149 + (m(6) * t16 - m(7) * t13 + t577) * t392 + (-t489 - t490 + t303 * mrSges(4,1) - t162 / 0.2e1 + t153 * mrSges(5,1) + t161 / 0.2e1 + Ifges(5,3) * t529 - Ifges(4,2) * t530 + t579 * t527 - t578 * t522) * t193 + (t272 * t580 + t378 * t578) * t344 / 0.2e1 + (t272 * t581 - t378 * t579) * t148 / 0.2e1 + (t148 * t581 + t344 * t580) * t526 - t89 * t582 + (Ifges(6,6) * t337 + Ifges(6,4) * t64 + Ifges(6,2) * t65 - t44 * mrSges(6,1) - t585 / 0.2e1 - t586 / 0.2e1 - t584 / 0.2e1 - t440 / 0.2e1 + t15 * mrSges(6,3) + t562) * t388 - t88 * t508; (t13 * t220 - t23 * t29 - t24 * t30) * m(7) + (-t113 * t120 + t15 * t224 + t16 * t223 - t75 * t87) * m(6) + t599 * (qJD(5) * t224 + (t383 - t434) * t355 + (-t150 + t316) * t350) + t360 - t178 * t436 - t398 * t442 / 0.2e1 + (t571 + (t381 - t376 / 0.2e1) * qJD(1)) * qJD(1) - (-Ifges(3,2) * t449 + t251 + t330) * t448 / 0.2e1 + t250 * t449 / 0.2e1 + t123 * t437 + t212 * t409 + (m(4) * t517 + mrSges(4,1) * t340 + t405 + t504) * t569 + (t180 * t190 - t181 * t191 - t303 * t436 + (t520 * t92 + t351 * t91 + (-t180 * t351 + t181 * t520) * qJD(3)) * pkin(2)) * m(4) + t327 * t124 + t323 * t126 + t316 * t215 + Ifges(3,5) * t287 + Ifges(3,6) * t286 - t268 * mrSges(3,2) - t269 * mrSges(3,1) + t220 * t17 + t223 * t54 + t224 * t55 - t158 * t177 + Ifges(3,3) * qJDD(2) - t87 * t151 - t120 * t107 - t30 * t93 - t29 * t94 + (-m(5) * t425 - m(6) * (t324 + t425) + t298 - m(7) * (t343 + t377) - m(4) * t343 + t563) * g(3) + (-(t414 + t547) * m(7) - t547 * m(6) - m(5) * t308 - t357 * t366 + t566) * g(1) + (-(t415 + t549) * m(7) - t549 * m(6) - m(5) * t306 - t353 * t366 + t567) * g(2) + (m(6) * t75 + m(7) * t396 + t382) * (qJD(5) * t223 + t316 * t355 + t350 * t434) + t594 * (-pkin(10) + t224) + (-t153 * t158 + t323 * t81 + t327 * t84 + (t316 - t191) * t169 + t595 * t157) * m(5) - t595 * t572 + t573 * t191 + t253 * t489 + t125 * t518; t360 + t563 * g(3) + (-t308 * t443 + t357 * t362 + t566) * g(1) + (-t306 * t443 + t353 * t362 + t567) * g(2) + t573 * t180 + (t500 + t572) * t181 + t288 * t54 + t289 * t55 + t280 * t17 + qJD(4) * t215 - t174 * t177 - t83 * t151 - t134 * t107 - pkin(3) * t124 + qJ(4) * t126 - t28 * t93 - t27 * t94 + t382 * t229 + t364 * t281 - t552 * t551 + (-t426 * g(3) - t113 * t134 + t15 * t289 + t16 * t288 + (t229 - t83) * t75 - t551 * t74) * m(6) + (-pkin(3) * t84 - t450 * g(3) + qJ(4) * t81 - t153 * t174 - t157 * t181 + t169 * t421) * m(5) + (-t414 * g(1) - t415 * g(2) - t377 * g(3) + t13 * t280 + t229 * t396 - t23 * t27 - t24 * t28 + t365 * t281 + t551 * t69) * m(7); -t197 * t94 - t198 * t93 - t345 * t215 + (-t107 + t177) * t253 + (qJD(5) * t382 - t345 * t151 + t577) * t355 + (-t339 * t552 + t364 + t55) * t350 + t124 + ((qJD(5) * t396 - t13) * t355 - t197 * t23 - t198 * t24 + (t339 * t69 + t365) * t350) * m(7) + (-t113 * t253 + t15 * t350 + t16 * t355 + t339 * (-t350 * t74 + t355 * t75)) * m(6) + (t153 * t253 - t169 * t345 + t84) * m(5) + (t341 * g(3) - t340 * t569) * t443; (t498 + t593) * t75 + (-t151 + t499) * t74 + qJD(6) * t427 + (-pkin(5) * t13 - t23 * t33 - t24 * t34) * m(7) + (t165 + t99) * t534 + t60 * t417 + (-t494 + t59) * t533 + t542 + t559 * t550 + t558 * t171 - t34 * t93 - t33 * t94 + (-t247 * t374 - t248 * t432 + t452) * g(3) + (-t216 * t374 - t217 * t432 + t454) * g(2) + (t218 * t432 - t219 * t374 + t453) * g(1) - pkin(5) * t17 + t594 * pkin(10) + t98 * t532; -t69 * (mrSges(7,1) * t139 + mrSges(7,2) * t138) + (Ifges(7,1) * t138 - t491) * t538 + t60 * t537 + (Ifges(7,5) * t138 - Ifges(7,6) * t139) * t535 - t23 * t93 + t24 * t94 - g(1) * (mrSges(7,1) * t187 - mrSges(7,2) * t188) - g(2) * (-mrSges(7,1) * t390 + mrSges(7,2) * t391) + g(3) * t402 * t248 + (t138 * t23 + t139 * t24) * mrSges(7,3) + t440 + (-Ifges(7,2) * t139 + t133 + t61) * t539 - t562;];
tau  = t1;
