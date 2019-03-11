% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:46
% EndTime: 2019-03-09 07:20:20
% DurationCPUTime: 17.89s
% Computational Cost: add. (15449->772), mult. (29581->1010), div. (0->0), fcn. (20151->14), ass. (0->374)
t609 = -mrSges(7,3) - mrSges(6,3);
t348 = sin(qJ(5));
t519 = cos(qJ(4));
t411 = t519 * qJD(4);
t400 = pkin(3) * t411;
t350 = sin(qJ(3));
t443 = qJD(1) * t350;
t358 = -pkin(1) - pkin(7);
t303 = qJD(1) * t358 + qJD(2);
t464 = t303 * t350;
t240 = -pkin(8) * t443 + t464;
t349 = sin(qJ(4));
t235 = t349 * t240;
t354 = cos(qJ(3));
t282 = t354 * t303;
t442 = qJD(1) * t354;
t241 = -pkin(8) * t442 + t282;
t174 = t241 * t519 - t235;
t413 = qJD(1) * t519;
t255 = -t349 * t442 - t350 * t413;
t257 = -t349 * t443 + t354 * t413;
t201 = pkin(4) * t257 - pkin(9) * t255;
t427 = pkin(3) * t442;
t181 = t201 + t427;
t353 = cos(qJ(5));
t97 = -t174 * t348 + t353 * t181;
t612 = -t348 * t400 - t97;
t98 = t353 * t174 + t348 * t181;
t611 = t353 * t400 - t98;
t437 = qJD(5) * t348;
t470 = t255 * t348;
t610 = t437 - t470;
t517 = pkin(3) * t349;
t323 = pkin(9) + t517;
t506 = -pkin(10) - t323;
t406 = qJD(5) * t506;
t430 = pkin(10) * t470;
t608 = t348 * t406 + t430 + t611;
t469 = t255 * t353;
t398 = t257 * pkin(5) - pkin(10) * t469;
t607 = t353 * t406 - t398 + t612;
t237 = qJD(3) * pkin(3) + t241;
t164 = t237 * t519 - t235;
t104 = t353 * t164 + t348 * t201;
t356 = -pkin(10) - pkin(9);
t419 = qJD(5) * t356;
t606 = t348 * t419 - t104 + t430;
t103 = -t164 * t348 + t353 * t201;
t605 = t353 * t419 - t103 - t398;
t346 = qJ(3) + qJ(4);
t336 = cos(t346);
t345 = qJ(5) + qJ(6);
t333 = sin(t345);
t501 = mrSges(7,2) * t333;
t502 = mrSges(6,2) * t348;
t604 = (t501 + t502) * t336;
t352 = cos(qJ(6));
t347 = sin(qJ(6));
t343 = qJD(3) + qJD(4);
t213 = -t257 * t348 + t343 * t353;
t236 = t519 * t240;
t165 = t349 * t237 + t236;
t152 = t343 * pkin(9) + t165;
t290 = pkin(3) * t443 + qJD(1) * qJ(2);
t168 = -pkin(4) * t255 - pkin(9) * t257 + t290;
t87 = t152 * t353 + t168 * t348;
t73 = pkin(10) * t213 + t87;
t482 = t347 * t73;
t250 = qJD(5) - t255;
t214 = t257 * t353 + t343 * t348;
t86 = -t152 * t348 + t353 * t168;
t72 = -pkin(10) * t214 + t86;
t59 = pkin(5) * t250 + t72;
t21 = t352 * t59 - t482;
t479 = t352 * t73;
t22 = t347 * t59 + t479;
t586 = Ifges(5,6) * t343;
t603 = -t290 * mrSges(5,1) - t86 * mrSges(6,1) - t21 * mrSges(7,1) + t87 * mrSges(6,2) + t22 * mrSges(7,2) + t586 / 0.2e1;
t583 = t343 * Ifges(5,5);
t602 = t290 * mrSges(5,2) + t583 / 0.2e1;
t601 = t354 / 0.2e1;
t138 = t213 * t347 + t214 * t352;
t239 = qJD(6) + t250;
t401 = t352 * t213 - t214 * t347;
t584 = Ifges(6,3) * t250;
t585 = Ifges(6,6) * t213;
t599 = Ifges(6,5) * t214 + Ifges(7,5) * t138 + Ifges(7,6) * t401 + Ifges(7,3) * t239 + t584 + t585;
t487 = t257 * mrSges(5,3);
t574 = mrSges(5,1) * t343 + mrSges(6,1) * t213 - mrSges(6,2) * t214 - t487;
t598 = t610 * pkin(5);
t173 = t241 * t349 + t236;
t438 = qJD(4) * t349;
t597 = pkin(3) * t438 - t173;
t212 = Ifges(6,4) * t213;
t118 = t214 * Ifges(6,1) + t250 * Ifges(6,5) + t212;
t248 = Ifges(5,4) * t255;
t596 = t257 * Ifges(5,1) + t353 * t118 + t248 + t583;
t436 = qJD(5) * t353;
t342 = qJDD(3) + qJDD(4);
t302 = qJDD(1) * t358 + qJDD(2);
t441 = qJD(3) * t350;
t217 = t354 * t302 - t303 * t441;
t433 = qJD(1) * qJD(3);
t278 = qJDD(1) * t354 - t350 * t433;
t177 = qJDD(3) * pkin(3) - pkin(8) * t278 + t217;
t440 = qJD(3) * t354;
t218 = t350 * t302 + t303 * t440;
t279 = -qJDD(1) * t350 - t354 * t433;
t186 = pkin(8) * t279 + t218;
t70 = t349 * t177 + t519 * t186 + t237 * t411 - t240 * t438;
t67 = pkin(9) * t342 + t70;
t371 = -t349 * t354 - t350 * t519;
t157 = qJD(1) * qJD(4) * t371 + t278 * t519 + t349 * t279;
t158 = -qJD(4) * t257 - t349 * t278 + t279 * t519;
t434 = qJD(1) * qJD(2);
t304 = qJDD(1) * qJ(2) + t434;
t230 = -pkin(3) * t279 + t304;
t82 = -pkin(4) * t158 - pkin(9) * t157 + t230;
t15 = -t152 * t437 + t168 * t436 + t348 * t82 + t353 * t67;
t16 = -qJD(5) * t87 - t348 * t67 + t353 * t82;
t595 = t15 * t353 - t16 * t348;
t504 = mrSges(6,1) * t353;
t594 = t502 - t504;
t509 = t353 * pkin(5);
t324 = pkin(4) + t509;
t334 = sin(t346);
t379 = t324 * t334 + t336 * t356;
t407 = -pkin(4) * t334 + t336 * pkin(9);
t593 = -m(6) * t407 + m(7) * t379;
t592 = mrSges(5,1) * t334 + (mrSges(5,2) + t609) * t336;
t549 = m(7) * pkin(5);
t101 = qJD(5) * t213 + t157 * t353 + t342 * t348;
t102 = -qJD(5) * t214 - t157 * t348 + t342 * t353;
t31 = qJD(6) * t401 + t101 * t352 + t102 * t347;
t547 = t31 / 0.2e1;
t32 = -qJD(6) * t138 - t101 * t347 + t102 * t352;
t546 = t32 / 0.2e1;
t541 = -m(3) - m(4);
t155 = qJDD(5) - t158;
t150 = qJDD(6) + t155;
t534 = t150 / 0.2e1;
t259 = t506 * t348;
t340 = t353 * pkin(10);
t463 = t323 * t353;
t260 = t340 + t463;
t198 = t259 * t347 + t260 * t352;
t588 = -qJD(6) * t198 - t347 * t608 + t352 * t607;
t197 = t259 * t352 - t260 * t347;
t587 = qJD(6) * t197 + t347 * t607 + t352 * t608;
t582 = -mrSges(6,1) - t549;
t292 = t356 * t348;
t293 = pkin(9) * t353 + t340;
t215 = t292 * t352 - t293 * t347;
t581 = qJD(6) * t215 + t347 * t605 + t352 * t606;
t216 = t292 * t347 + t293 * t352;
t580 = -qJD(6) * t216 - t347 * t606 + t352 * t605;
t571 = t519 * qJD(3) + t411;
t208 = -t349 * t441 - t350 * t438 + t354 * t571;
t271 = t347 * t353 + t348 * t352;
t378 = t347 * t348 - t352 * t353;
t565 = qJD(5) + qJD(6);
t575 = t378 * t371;
t579 = t378 * qJD(1) - t271 * t208 - t565 * t575;
t206 = t565 * t271;
t578 = -t271 * qJD(1) + t206 * t371 - t208 * t378;
t144 = mrSges(5,1) * t342 - mrSges(5,3) * t157;
t48 = -mrSges(6,1) * t102 + mrSges(6,2) * t101;
t577 = t48 - t144;
t151 = -t343 * pkin(4) - t164;
t392 = mrSges(6,1) * t348 + mrSges(6,2) * t353;
t576 = t151 * t392;
t269 = t349 * t350 - t354 * t519;
t190 = t271 * t269;
t339 = t350 * pkin(3);
t320 = qJ(2) + t339;
t202 = -pkin(4) * t371 + pkin(9) * t269 + t320;
t507 = pkin(8) - t358;
t284 = t507 * t350;
t285 = t507 * t354;
t211 = -t284 * t519 - t349 * t285;
t203 = t353 * t211;
t123 = t348 * t202 + t203;
t573 = t349 * t284 - t519 * t285;
t572 = t597 + t598;
t570 = -t165 + t598;
t380 = t217 * t354 + t218 * t350;
t568 = -t348 * t86 + t353 * t87;
t351 = sin(qJ(1));
t355 = cos(qJ(1));
t567 = -g(1) * t351 + g(2) * t355;
t566 = -m(6) - m(7) - m(5);
t395 = mrSges(4,1) * t354 - mrSges(4,2) * t350;
t496 = Ifges(4,4) * t354;
t564 = (-Ifges(4,1) * t350 - t496) * t601 + qJ(2) * t395;
t460 = t334 * t356;
t461 = t334 * t355;
t563 = -mrSges(5,2) * t461 + (-m(7) * t460 - t604) * t355;
t421 = t336 * t504;
t335 = cos(t345);
t503 = mrSges(7,1) * t335;
t423 = t336 * t503;
t457 = t336 * t351;
t462 = t334 * t351;
t562 = -mrSges(5,1) * t457 + t609 * t462 + (-t423 - t421 + t604) * t351;
t561 = t592 + (-t501 + t503 - t594) * t334;
t560 = mrSges(5,1) * t336 - t334 * t609 + t421;
t559 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t288 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t443;
t289 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t442;
t558 = (t288 * t354 - t289 * t350) * qJD(3);
t207 = -t349 * t440 - t350 * t571 - t354 * t438;
t71 = t177 * t519 - t349 * t186 - t237 * t438 - t240 * t411;
t557 = t164 * t207 + t165 * t208 - t269 * t71 - t371 * t70;
t11 = pkin(5) * t155 - pkin(10) * t101 + t16;
t13 = pkin(10) * t102 + t15;
t3 = qJD(6) * t21 + t11 * t347 + t13 * t352;
t4 = -qJD(6) * t22 + t11 * t352 - t13 * t347;
t556 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t555 = m(4) * t380 + t354 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t278) + t350 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t279);
t554 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t499 = mrSges(6,3) * t213;
t160 = -mrSges(6,2) * t250 + t499;
t498 = mrSges(6,3) * t214;
t161 = mrSges(6,1) * t250 - t498;
t384 = t348 * t87 + t353 * t86;
t361 = -qJD(5) * t384 + t595;
t65 = mrSges(6,1) * t155 - mrSges(6,3) * t101;
t481 = t348 * t65;
t553 = m(6) * t361 - t160 * t437 - t161 * t436 - t481;
t394 = mrSges(4,1) * t350 + mrSges(4,2) * t354;
t552 = -t394 - mrSges(3,3) + mrSges(2,2) - t592 - t593;
t551 = qJD(1) ^ 2;
t550 = Ifges(7,4) * t547 + Ifges(7,2) * t546 + Ifges(7,6) * t534;
t548 = Ifges(7,1) * t547 + Ifges(7,4) * t546 + Ifges(7,5) * t534;
t491 = Ifges(7,4) * t138;
t62 = Ifges(7,2) * t401 + Ifges(7,6) * t239 + t491;
t545 = -t62 / 0.2e1;
t544 = t62 / 0.2e1;
t133 = Ifges(7,4) * t401;
t63 = Ifges(7,1) * t138 + Ifges(7,5) * t239 + t133;
t543 = -t63 / 0.2e1;
t542 = t63 / 0.2e1;
t540 = t101 / 0.2e1;
t539 = t102 / 0.2e1;
t538 = -t401 / 0.2e1;
t537 = t401 / 0.2e1;
t536 = -t138 / 0.2e1;
t535 = t138 / 0.2e1;
t533 = t155 / 0.2e1;
t530 = -t213 / 0.2e1;
t529 = -t214 / 0.2e1;
t528 = t214 / 0.2e1;
t527 = -t239 / 0.2e1;
t526 = t239 / 0.2e1;
t525 = -t250 / 0.2e1;
t523 = t255 / 0.2e1;
t521 = t257 / 0.2e1;
t518 = mrSges(7,3) * t21;
t516 = pkin(3) * t354;
t514 = pkin(5) * t214;
t511 = g(3) * t336;
t510 = t22 * mrSges(7,3);
t500 = mrSges(5,3) * t255;
t497 = Ifges(4,4) * t350;
t495 = Ifges(5,4) * t257;
t494 = Ifges(6,4) * t214;
t493 = Ifges(6,4) * t348;
t492 = Ifges(6,4) * t353;
t490 = pkin(5) * qJD(6);
t68 = -t342 * pkin(4) - t71;
t486 = t269 * t68;
t66 = -mrSges(6,2) * t155 + mrSges(6,3) * t102;
t478 = t353 * t66;
t474 = t207 * t353;
t468 = t269 * t348;
t467 = t269 * t353;
t459 = t335 * t351;
t458 = t335 * t355;
t117 = Ifges(6,2) * t213 + Ifges(6,6) * t250 + t494;
t456 = t348 * t117;
t455 = t348 * t351;
t454 = t348 * t355;
t451 = t351 * t353;
t449 = t353 * t355;
t231 = -t333 * t462 + t458;
t232 = t333 * t355 + t334 * t459;
t447 = t231 * mrSges(7,1) - t232 * mrSges(7,2);
t233 = t333 * t461 + t459;
t234 = -t333 * t351 + t334 * t458;
t446 = t233 * mrSges(7,1) + t234 * mrSges(7,2);
t445 = pkin(4) * t457 + pkin(9) * t462;
t444 = t355 * pkin(1) + t351 * qJ(2);
t435 = qJDD(1) * mrSges(3,2);
t305 = pkin(3) * t440 + qJD(2);
t431 = m(5) * t516;
t429 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t150;
t428 = t519 * pkin(3);
t420 = Ifges(6,5) * t101 + Ifges(6,6) * t102 + Ifges(6,3) * t155;
t414 = -t456 / 0.2e1;
t410 = t436 / 0.2e1;
t338 = t355 * qJ(2);
t408 = -pkin(1) * t351 + t338;
t405 = -t433 / 0.2e1;
t403 = (t304 + t434) * qJ(2);
t121 = pkin(4) * t208 - pkin(9) * t207 + t305;
t263 = t507 * t441;
t264 = qJD(3) * t285;
t126 = qJD(4) * t573 + t349 * t263 - t519 * t264;
t402 = t353 * t121 - t126 * t348;
t122 = t353 * t202 - t211 * t348;
t325 = -t428 - pkin(4);
t397 = -pkin(4) * t336 - pkin(9) * t334;
t396 = t324 * t457 - t351 * t460;
t391 = -mrSges(7,1) * t333 - mrSges(7,2) * t335;
t390 = t354 * Ifges(4,1) - t497;
t389 = Ifges(6,1) * t353 - t493;
t388 = -t350 * Ifges(4,2) + t496;
t387 = -Ifges(6,2) * t348 + t492;
t386 = -Ifges(4,5) * t350 - Ifges(4,6) * t354;
t385 = Ifges(6,5) * t353 - Ifges(6,6) * t348;
t105 = pkin(10) * t468 + t123;
t88 = -pkin(5) * t371 + pkin(10) * t467 + t122;
t45 = t105 * t352 + t347 * t88;
t44 = -t105 * t347 + t352 * t88;
t381 = -t160 * t348 - t161 * t353;
t373 = t429 + t556;
t246 = t334 * t454 + t451;
t244 = -t334 * t455 + t449;
t370 = -t207 * t348 + t269 * t436;
t369 = t269 * t437 + t474;
t367 = t350 * (-Ifges(4,2) * t354 - t497);
t42 = t348 * t121 + t353 * t126 + t202 * t436 - t211 * t437;
t127 = qJD(4) * t211 - t519 * t263 - t349 * t264;
t205 = t565 * t378;
t114 = -t213 * pkin(5) + t151;
t171 = t271 * t255;
t172 = t378 * t255;
t182 = Ifges(5,2) * t255 + t495 + t586;
t39 = t101 * Ifges(6,4) + t102 * Ifges(6,2) + t155 * Ifges(6,6);
t40 = t101 * Ifges(6,1) + t102 * Ifges(6,4) + t155 * Ifges(6,5);
t41 = -t102 * pkin(5) + t68;
t360 = (-Ifges(7,4) * t172 - Ifges(7,2) * t171) * t538 + (t576 + t414) * qJD(5) + (-Ifges(7,5) * t172 - Ifges(7,6) * t171) * t527 + (-t610 * t87 + (-t436 + t469) * t86 + t595) * mrSges(6,3) + (-Ifges(7,1) * t172 - Ifges(7,4) * t171) * t536 - (Ifges(7,4) * t535 + Ifges(7,2) * t537 + Ifges(7,6) * t526 + t510 + t544) * t206 - (Ifges(7,1) * t535 + Ifges(7,4) * t537 + Ifges(7,5) * t526 - t518 + t542) * t205 + (Ifges(6,2) * t353 + t493) * t539 + (Ifges(6,1) * t348 + t492) * t540 - t172 * t543 - t171 * t545 + t271 * t548 + t118 * t410 - (Ifges(5,1) * t255 - t495 + t599) * t257 / 0.2e1 - (-Ifges(5,2) * t257 + t248 + t596) * t255 / 0.2e1 + t68 * t594 + ((t172 - t205) * mrSges(7,2) + (-t171 + t206) * mrSges(7,1)) * t114 + (t171 * t22 - t172 * t21 - t271 * t4 - t3 * t378) * mrSges(7,3) + (Ifges(7,1) * t271 - Ifges(7,4) * t378) * t547 + (Ifges(7,5) * t271 - Ifges(7,6) * t378) * t534 + t41 * (mrSges(7,1) * t378 + mrSges(7,2) * t271) + (Ifges(7,4) * t271 - Ifges(7,2) * t378) * t546 - t378 * t550 + (t213 * t387 + t214 * t389 + t250 * t385) * qJD(5) / 0.2e1 + t353 * t39 / 0.2e1 + t348 * t40 / 0.2e1 + Ifges(5,3) * t342 + (Ifges(6,5) * t348 + Ifges(6,6) * t353) * t533 + t182 * t521 + t456 * t523 + t164 * t500 + t165 * t487 + Ifges(5,5) * t157 + Ifges(5,6) * t158 + (t385 * t525 + t387 * t530 + t389 * t529 - t576 - t602) * t255 + (Ifges(6,5) * t529 + Ifges(7,5) * t536 + Ifges(6,6) * t530 + Ifges(7,6) * t538 + Ifges(6,3) * t525 + Ifges(7,3) * t527 + t603) * t257 - t70 * mrSges(5,2) + t71 * mrSges(5,1);
t357 = -pkin(8) - pkin(7);
t331 = -pkin(1) * qJDD(1) + qJDD(2);
t316 = t351 * t516;
t287 = t325 - t509;
t272 = t394 * qJD(1);
t253 = Ifges(4,5) * qJD(3) + qJD(1) * t390;
t252 = Ifges(4,6) * qJD(3) + qJD(1) * t388;
t247 = t334 * t449 - t455;
t245 = t334 * t451 + t454;
t223 = -mrSges(5,2) * t343 + t500;
t200 = -mrSges(5,1) * t255 + mrSges(5,2) * t257;
t192 = t378 * t269;
t189 = t271 * t371;
t166 = -pkin(5) * t468 - t573;
t145 = -mrSges(5,2) * t342 + mrSges(5,3) * t158;
t111 = mrSges(7,1) * t239 - mrSges(7,3) * t138;
t110 = -mrSges(7,2) * t239 + mrSges(7,3) * t401;
t83 = -pkin(5) * t370 + t127;
t81 = -mrSges(7,1) * t401 + mrSges(7,2) * t138;
t57 = -t205 * t269 - t207 * t271;
t55 = t190 * t565 - t378 * t207;
t43 = -qJD(5) * t123 + t402;
t27 = t352 * t72 - t482;
t26 = -t347 * t72 - t479;
t25 = pkin(10) * t370 + t42;
t20 = -pkin(10) * t474 + pkin(5) * t208 + (-t203 + (-pkin(10) * t269 - t202) * t348) * qJD(5) + t402;
t19 = -mrSges(7,2) * t150 + mrSges(7,3) * t32;
t18 = mrSges(7,1) * t150 - mrSges(7,3) * t31;
t12 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t6 = -qJD(6) * t45 + t20 * t352 - t25 * t347;
t5 = qJD(6) * t44 + t20 * t347 + t25 * t352;
t1 = [(Ifges(7,1) * t192 + Ifges(7,4) * t190) * t547 + (Ifges(7,1) * t55 + Ifges(7,4) * t57) * t535 + (Ifges(6,1) * t369 + Ifges(6,4) * t370) * t528 + (-t230 * mrSges(5,2) - Ifges(5,1) * t157 - Ifges(5,4) * t158 - Ifges(5,5) * t342 + t117 * t410 - t385 * t533 - t387 * t539 - t389 * t540) * t269 + m(6) * (t122 * t16 + t123 * t15 + t42 * t87 + t43 * t86) + m(5) * (t126 * t165 + t211 * t70 + t230 * t320 + t290 * t305) + (t190 * t3 - t192 * t4 - t21 * t55 + t22 * t57) * mrSges(7,3) + (t15 * t468 + t16 * t467 - t369 * t86 + t370 * t87) * mrSges(6,3) + (qJD(5) * t118 + t39) * t468 / 0.2e1 + (Ifges(4,1) * t278 + Ifges(4,4) * t279) * t601 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + (-m(5) * t164 + m(6) * t151 - t574) * t127 - t380 * mrSges(4,3) + t564 * t433 - t557 * mrSges(5,3) + t358 * t558 + (Ifges(7,5) * t192 + Ifges(7,6) * t190) * t534 + (Ifges(7,5) * t55 + Ifges(7,6) * t57) * t526 + t250 * (Ifges(6,5) * t369 + Ifges(6,6) * t370) / 0.2e1 + (Ifges(7,4) * t192 + Ifges(7,2) * t190) * t546 + (Ifges(7,4) * t55 + Ifges(7,2) * t57) * t537 + t213 * (Ifges(6,4) * t369 + Ifges(6,2) * t370) / 0.2e1 - t350 * (Ifges(4,4) * t278 + Ifges(4,2) * t279) / 0.2e1 + m(4) * t403 - (-m(5) * t71 + m(6) * t68 + t577) * t573 + t55 * t542 + t57 * t544 + t192 * t548 + t190 * t550 + t151 * (-mrSges(6,1) * t370 + mrSges(6,2) * t369) - t392 * t486 + m(7) * (t114 * t83 + t166 * t41 + t21 * t6 + t22 * t5 + t3 * t45 + t4 * t44) + t367 * t405 + t555 * t358 - t40 * t467 / 0.2e1 - t252 * t440 / 0.2e1 - t253 * t441 / 0.2e1 - pkin(1) * t435 + m(3) * (-pkin(1) * t331 + t403) - (t429 + t420) * t371 / 0.2e1 - (mrSges(5,1) * t230 - Ifges(5,4) * t157 + Ifges(6,5) * t540 + Ifges(7,5) * t547 - Ifges(5,2) * t158 - Ifges(5,6) * t342 + Ifges(6,6) * t539 + Ifges(7,6) * t546 + Ifges(6,3) * t533 + Ifges(7,3) * t534 + t554 + t556) * t371 + (t394 + 0.2e1 * mrSges(3,3)) * t304 + qJDD(3) * (Ifges(4,5) * t354 - Ifges(4,6) * t350) + t331 * mrSges(3,2) + t320 * (-mrSges(5,1) * t158 + mrSges(5,2) * t157) + t305 * t200 + qJ(2) * (-mrSges(4,1) * t279 + mrSges(4,2) * t278) + qJD(2) * t272 + t126 * t223 + t211 * t145 + t41 * (-mrSges(7,1) * t190 + mrSges(7,2) * t192) + t43 * t161 + t166 * t12 + t42 * t160 + t122 * t65 + t123 * t66 + (t414 + Ifges(5,1) * t521 + Ifges(5,4) * t523 + t596 / 0.2e1 + t602) * t207 + (t599 / 0.2e1 + t584 / 0.2e1 + t585 / 0.2e1 + Ifges(6,5) * t528 + Ifges(7,5) * t535 + Ifges(7,6) * t537 - Ifges(5,4) * t521 - Ifges(5,2) * t523 + Ifges(7,3) * t526 - t182 / 0.2e1 - t603) * t208 + t44 * t18 + t45 * t19 + (-t454 * t549 - t245 * mrSges(6,1) - t232 * mrSges(7,1) - t244 * mrSges(6,2) - t231 * mrSges(7,2) + t541 * t444 + t566 * (t351 * t339 - t355 * t357 + t444) + (-m(4) * pkin(7) + t559) * t355 + t552 * t351) * g(2) + (t455 * t549 - m(3) * t408 - m(4) * t338 - t247 * mrSges(6,1) - t234 * mrSges(7,1) + t246 * mrSges(6,2) + t233 * mrSges(7,2) + t566 * (t355 * t339 + t351 * t357 + t408) + (-m(4) * t358 - t559) * t351 + t552 * t355) * g(1) + qJD(3) ^ 2 * t386 / 0.2e1 + t279 * t388 / 0.2e1 + t278 * t390 / 0.2e1 + t83 * t81 + t5 * t110 + t6 * t111 + t114 * (-mrSges(7,1) * t57 + mrSges(7,2) * t55); t435 + t189 * t18 + t575 * t19 + t579 * t111 + t578 * t110 + t558 + (qJ(2) * t541 - mrSges(3,3)) * t551 + (t12 + t577) * t269 + (t160 * t353 - t161 * t348 + t223) * t208 + (-t81 + t574) * t207 - (qJD(5) * t381 + t145 + t478 - t481) * t371 + (-m(5) * t290 - t200 - t272 + t381) * qJD(1) + m(5) * t557 + m(3) * t331 + t567 * (-t541 - t566) + (-t114 * t207 + t189 * t4 + t21 * t579 + t22 * t578 + t269 * t41 + t3 * t575) * m(7) + (-t384 * qJD(1) - t151 * t207 + t208 * t568 - t361 * t371 + t486) * m(6) + t555; (t164 * t173 - t165 * t174 - t290 * t427 + (t519 * t71 + t349 * t70 + (-t164 * t349 + t165 * t519) * qJD(4)) * pkin(3)) * m(5) + (t400 - t174) * t223 + t572 * t81 + t567 * t395 + (t68 * t325 + (t151 * t349 + t519 * t568) * qJD(4) * pkin(3) - t151 * t173 - t86 * t97 - t87 * t98) * m(6) + (t367 / 0.2e1 - t564) * t551 + (-m(6) * (t316 + t445) - (-mrSges(5,2) * t334 + t431) * t351 - m(7) * (t316 + t396) + t562) * g(1) + ((-m(7) * (-t324 * t336 - t516) + t423 - m(6) * (t397 - t516) + t431 + t560) * t355 + t563) * g(2) + t553 * t323 + t360 + t386 * t405 + t611 * t160 + t612 * t161 - t288 * t282 + t253 * t443 / 0.2e1 + t252 * t442 / 0.2e1 - t200 * t427 - t574 * t597 + (m(5) * t339 - m(6) * (-t339 + t407) + t394 - m(7) * (-t339 - t379) + t561) * g(3) + Ifges(4,3) * qJDD(3) + t145 * t517 + t66 * t463 + t289 * t464 + t325 * t48 + Ifges(4,5) * t278 + Ifges(4,6) * t279 + t287 * t12 + t217 * mrSges(4,1) - t218 * mrSges(4,2) + t197 * t18 + t198 * t19 + t144 * t428 + t587 * t110 + t588 * t111 + (t572 * t114 + t197 * t4 + t198 * t3 + t588 * t21 + t587 * t22 + t287 * t41) * m(7); -pkin(4) * t48 - t103 * t161 - t104 * t160 - t324 * t12 - t164 * t223 + t215 * t18 + t216 * t19 + t360 + t570 * t81 + t574 * t165 + t580 * t111 + t581 * t110 + (t114 * t570 + t21 * t580 + t215 * t4 + t216 * t3 + t22 * t581 - t324 * t41) * m(7) + (-pkin(4) * t68 - t103 * t86 - t104 * t87 - t151 * t165) * m(6) + (t561 + t593) * g(3) + ((-(-m(7) * t324 - t503) * t336 - m(6) * t397 + t560) * t355 + t563) * g(2) + (-m(6) * t445 - m(7) * t396 + mrSges(5,2) * t462 + t562) * g(1) + (t478 + t553) * pkin(9); t554 - (mrSges(7,1) * t114 + Ifges(7,4) * t536 + Ifges(7,2) * t538 + Ifges(7,6) * t527 - t510 + t545) * t138 + (-mrSges(7,2) * t114 + Ifges(7,1) * t536 + Ifges(7,4) * t538 + Ifges(7,5) * t527 + t518 + t543) * t401 + (mrSges(6,2) * t245 + t244 * t582 - t447) * g(1) + (-mrSges(6,2) * t247 + t246 * t582 - t446) * g(2) + t420 + (-t347 * t490 - t26) * t111 + (t3 * t347 + t352 * t4 + (-t21 * t347 + t22 * t352) * qJD(6)) * t549 - t81 * t514 - m(7) * (t114 * t514 + t21 * t26 + t22 * t27) + (t352 * t18 + t347 * t19) * pkin(5) + (t352 * t490 - t27) * t110 + (t498 + t161) * t87 + (t499 - t160) * t86 + (-Ifges(6,2) * t214 + t118 + t212) * t530 + (t348 * t549 - t391 + t392) * t511 + t117 * t528 + (Ifges(6,1) * t213 - t494) * t529 + (Ifges(6,5) * t213 - Ifges(6,6) * t214) * t525 - t151 * (mrSges(6,1) * t214 + mrSges(6,2) * t213) + t373; -t114 * (mrSges(7,1) * t138 + mrSges(7,2) * t401) + (Ifges(7,1) * t401 - t491) * t536 + t62 * t535 + (Ifges(7,5) * t401 - Ifges(7,6) * t138) * t527 - t21 * t110 + t22 * t111 - g(1) * t447 - g(2) * t446 - t391 * t511 + (t138 * t22 + t21 * t401) * mrSges(7,3) + t373 + (-Ifges(7,2) * t138 + t133 + t63) * t538;];
tau  = t1;
