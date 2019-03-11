% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:31
% EndTime: 2019-03-09 10:04:14
% DurationCPUTime: 28.39s
% Computational Cost: add. (4834->721), mult. (9968->882), div. (0->0), fcn. (5587->6), ass. (0->326)
t536 = Ifges(7,4) + Ifges(6,5);
t225 = cos(qJ(2));
t544 = -t225 / 0.2e1;
t537 = Ifges(6,4) + Ifges(5,5);
t527 = -Ifges(7,5) + t537;
t528 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t221 = sin(qJ(4));
t224 = cos(qJ(4));
t376 = qJD(1) * t225;
t152 = qJD(2) * t221 + t224 * t376;
t222 = sin(qJ(2));
t377 = qJD(1) * t222;
t189 = qJD(4) + t377;
t343 = t221 * t376;
t374 = qJD(2) * t224;
t153 = -t343 + t374;
t411 = t153 * Ifges(5,4);
t57 = -t152 * Ifges(5,2) + t189 * Ifges(5,6) + t411;
t543 = t57 / 0.2e1;
t415 = Ifges(4,6) * t222;
t542 = Ifges(4,3) * t544 - t415 / 0.2e1;
t423 = Ifges(3,4) * t222;
t541 = Ifges(3,2) * t544 - t423 / 0.2e1;
t540 = t536 * t224;
t539 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t361 = qJD(1) * qJD(2);
t161 = -t225 * qJDD(1) + t222 * t361;
t371 = qJD(4) * t152;
t72 = qJDD(2) * t224 + t161 * t221 - t371;
t460 = t72 / 0.2e1;
t368 = qJD(4) * t224;
t73 = qJD(2) * t368 - qJD(4) * t343 + qJDD(2) * t221 - t224 * t161;
t458 = t73 / 0.2e1;
t162 = qJDD(1) * t222 + t225 * t361;
t151 = qJDD(4) + t162;
t453 = t151 / 0.2e1;
t538 = t376 / 0.2e1;
t510 = -Ifges(4,4) + Ifges(3,5);
t508 = Ifges(4,5) - Ifges(3,6);
t535 = Ifges(6,2) + Ifges(5,3);
t534 = Ifges(7,2) + Ifges(6,3);
t533 = -Ifges(5,6) + Ifges(6,6);
t532 = Ifges(6,6) - Ifges(7,6);
t417 = Ifges(6,5) * t221;
t279 = -Ifges(6,3) * t224 + t417;
t419 = Ifges(7,4) * t221;
t284 = -Ifges(7,2) * t224 + t419;
t531 = t279 + t284;
t281 = Ifges(5,5) * t221 + Ifges(5,6) * t224;
t286 = Ifges(6,4) * t221 - Ifges(6,6) * t224;
t530 = t281 + t286;
t198 = pkin(7) * t377;
t199 = pkin(3) * t377;
t383 = t198 + t199;
t529 = qJD(3) + t383;
t299 = -t221 * mrSges(7,1) + t224 * mrSges(7,2);
t302 = mrSges(5,1) * t221 + mrSges(5,2) * t224;
t403 = qJ(5) * t224;
t410 = t221 * mrSges(6,1);
t526 = -m(7) * (pkin(5) * t221 - t403) + t299 - t410 - (-m(6) * qJ(5) - mrSges(6,3)) * t224 - t302;
t457 = pkin(2) + pkin(8);
t386 = qJ(6) - t457;
t525 = -m(7) * t386 - t539 + (m(5) + m(6)) * t457;
t420 = Ifges(5,4) * t224;
t471 = t221 * t528 + t420 - t540;
t146 = Ifges(5,4) * t152;
t516 = t536 * t152;
t475 = t528 * t153 + t527 * t189 - t146 + t516;
t519 = t536 * t153;
t504 = t534 * t152 + t532 * t189 + t519;
t524 = t221 * t475 / 0.2e1 + (t543 - t504 / 0.2e1) * t224 + (Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1) * qJD(2) + (t542 + t541) * qJD(1);
t207 = t222 * qJ(3);
t323 = -pkin(1) - t207;
t106 = (-t225 * t457 + t323) * qJD(1);
t362 = qJD(3) + t198;
t108 = -qJD(2) * t457 + t199 + t362;
t369 = qJD(4) * t221;
t372 = qJD(3) * t222;
t402 = qJDD(1) * pkin(1);
t247 = -qJ(3) * t162 - qJD(1) * t372 - t402;
t43 = t161 * t457 + t247;
t150 = t162 * pkin(7);
t324 = qJDD(3) + t150;
t76 = pkin(3) * t162 - qJDD(2) * t457 + t324;
t6 = -t106 * t369 + t108 * t368 + t221 * t76 + t224 * t43;
t3 = t151 * qJ(5) + t189 * qJD(5) + t6;
t2 = qJ(6) * t73 + qJD(6) * t152 + t3;
t523 = t2 * mrSges(7,2);
t522 = t3 * mrSges(6,3);
t521 = t6 * mrSges(5,2);
t298 = t225 * mrSges(4,2) - t222 * mrSges(4,3);
t305 = mrSges(3,1) * t225 - mrSges(3,2) * t222;
t518 = t298 - t305;
t517 = -qJD(5) * t224 + t529;
t454 = -t151 / 0.2e1;
t514 = -t72 * Ifges(5,4) / 0.2e1 + Ifges(5,6) * t454 + t536 * t460 + t532 * t453 + (Ifges(5,2) + t534) * t458;
t180 = t189 * qJ(5);
t45 = t224 * t106 + t221 * t108;
t28 = qJ(6) * t152 + t45;
t25 = t180 + t28;
t32 = t180 + t45;
t513 = -t32 * mrSges(6,2) - t45 * mrSges(5,3) + t25 * mrSges(7,3);
t512 = -m(4) - m(5);
t511 = -m(7) - m(6);
t450 = -t153 / 0.2e1;
t36 = mrSges(5,1) * t151 - mrSges(5,3) * t72;
t37 = -t151 * mrSges(6,1) + t72 * mrSges(6,2);
t506 = t36 - t37;
t39 = -mrSges(5,2) * t151 - mrSges(5,3) * t73;
t40 = -mrSges(6,2) * t73 + mrSges(6,3) * t151;
t505 = t40 + t39;
t404 = qJ(5) * t221;
t455 = pkin(4) + pkin(5);
t260 = t224 * t455 + t404;
t503 = -t189 * t260 - t517;
t425 = mrSges(7,3) * t152;
t94 = mrSges(7,2) * t189 + t425;
t429 = mrSges(6,2) * t152;
t96 = mrSges(6,3) * t189 - t429;
t430 = -t94 - t96;
t427 = mrSges(5,3) * t152;
t95 = -mrSges(5,2) * t189 - t427;
t502 = t96 + t95;
t426 = mrSges(5,3) * t153;
t98 = mrSges(5,1) * t189 - t426;
t428 = mrSges(6,2) * t153;
t99 = -mrSges(6,1) * t189 + t428;
t501 = -t98 + t99;
t414 = Ifges(4,6) * t225;
t276 = -t222 * Ifges(4,2) - t414;
t500 = Ifges(4,4) * qJD(2) + t153 * Ifges(7,5) + t152 * Ifges(7,6) - t189 * Ifges(7,3) + qJD(1) * t276;
t318 = qJD(4) * t386;
t363 = t224 * qJD(6);
t397 = t221 * t222;
t200 = pkin(2) * t377;
t406 = qJ(3) * t225;
t272 = pkin(8) * t222 - t406;
t114 = qJD(1) * t272 + t200;
t201 = pkin(7) * t376;
t158 = pkin(3) * t376 + t201;
t68 = -t221 * t114 + t158 * t224;
t499 = -(-qJ(6) * t397 - t225 * t455) * qJD(1) + t68 + t221 * t318 - t363;
t364 = t221 * qJD(6);
t69 = t224 * t114 + t221 * t158;
t49 = qJ(5) * t376 + t69;
t498 = t364 - t49 + (qJ(6) * t377 + t318) * t224;
t349 = mrSges(4,1) * t376;
t171 = -qJD(2) * mrSges(4,3) - t349;
t84 = mrSges(5,1) * t152 + mrSges(5,2) * t153;
t497 = t84 - t171;
t274 = pkin(4) * t224 + t404;
t496 = t189 * t274 + t517;
t495 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t376 - t171;
t350 = mrSges(4,1) * t377;
t494 = -mrSges(3,3) * t377 - t350 + (mrSges(3,1) - mrSges(4,2)) * qJD(2);
t493 = t530 * t222 + t535 * t225;
t492 = t531 * t222 + t532 * t225;
t491 = t222 * (-Ifges(4,2) * t225 + t415) + t225 * (Ifges(4,3) * t222 - t414);
t490 = t534 * t221 + t540;
t489 = t533 * t221 + t537 * t224;
t488 = t508 * t222 + t510 * t225;
t487 = t535 * t151 + t533 * t73 + t537 * t72;
t149 = t161 * pkin(7);
t486 = -t149 * t225 + t150 * t222;
t107 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t149;
t113 = -qJDD(2) * pkin(2) + t324;
t485 = -t107 * t225 + t113 * t222;
t7 = -t106 * t368 - t108 * t369 - t221 * t43 + t224 * t76;
t483 = t221 * t6 + t224 * t7;
t258 = qJDD(5) - t7;
t5 = -pkin(4) * t151 + t258;
t482 = t221 * t3 - t224 * t5;
t1 = -qJ(6) * t72 - qJD(6) * t153 - t151 * t455 + t258;
t481 = -t1 * t224 + t2 * t221;
t223 = sin(qJ(1));
t226 = cos(qJ(1));
t480 = g(1) * t226 + g(2) * t223;
t479 = -m(5) + t511;
t478 = -mrSges(4,1) - mrSges(3,3) + mrSges(2,2);
t477 = mrSges(7,2) + mrSges(6,3) - mrSges(5,2);
t476 = (-Ifges(5,4) + t536) * t73 + t528 * t72 + t527 * t151;
t197 = Ifges(3,4) * t376;
t474 = Ifges(3,1) * t377 + Ifges(3,5) * qJD(2) + t533 * t152 + t537 * t153 + t535 * t189 + t197;
t473 = -mrSges(2,1) + t518;
t472 = t471 * t222 + t527 * t225;
t421 = Ifges(5,4) * t221;
t470 = t528 * t224 + t417 + t419 - t421;
t468 = t221 * t455 - t403;
t408 = t225 * mrSges(4,3);
t467 = -t408 + t526 * t225 + (m(4) * pkin(2) - mrSges(4,2) + t525) * t222;
t465 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1) + mrSges(7,1);
t464 = m(7) * qJ(6) + t539;
t220 = qJD(2) * qJ(3);
t128 = t220 + t158;
t264 = qJ(5) * t153 - t128;
t29 = -t152 * t455 + qJD(6) + t264;
t300 = t224 * mrSges(7,1) + t221 * mrSges(7,2);
t301 = t224 * mrSges(6,1) + t221 * mrSges(6,3);
t303 = mrSges(5,1) * t224 - mrSges(5,2) * t221;
t47 = pkin(4) * t152 - t264;
t463 = t128 * t303 - t29 * t300 + t301 * t47;
t459 = -t73 / 0.2e1;
t456 = pkin(3) + pkin(7);
t452 = -t152 / 0.2e1;
t451 = t152 / 0.2e1;
t449 = t153 / 0.2e1;
t448 = -t189 / 0.2e1;
t447 = t189 / 0.2e1;
t442 = pkin(7) * t222;
t439 = g(3) * t225;
t213 = t225 * pkin(2);
t211 = t225 * pkin(7);
t82 = mrSges(6,1) * t152 - mrSges(6,3) * t153;
t83 = -mrSges(7,1) * t152 + mrSges(7,2) * t153;
t431 = t82 - t83;
t424 = mrSges(7,3) * t153;
t422 = Ifges(3,4) * t225;
t44 = -t221 * t106 + t224 * t108;
t405 = qJ(5) * t152;
t396 = t221 * t223;
t395 = t221 * t225;
t393 = t222 * t224;
t392 = t222 * t226;
t391 = t223 * t224;
t390 = t224 * t225;
t389 = t224 * t226;
t387 = t225 * t226;
t382 = t213 + t207;
t348 = t225 * pkin(8) + t382;
t143 = -pkin(1) - t348;
t176 = t456 * t222;
t79 = t224 * t143 + t221 * t176;
t177 = t225 * pkin(3) + t211;
t214 = t226 * pkin(7);
t381 = t226 * pkin(3) + t214;
t380 = t226 * pkin(1) + t223 * pkin(7);
t375 = qJD(2) * t222;
t373 = qJD(2) * t225;
t367 = qJD(4) * t225;
t366 = qJD(4) * t457;
t365 = qJD(5) * t221;
t356 = pkin(4) * t395;
t351 = m(4) - t479;
t74 = t222 * qJ(5) + t79;
t347 = t221 * t375;
t346 = t221 * t367;
t21 = -t73 * mrSges(7,1) + t72 * mrSges(7,2);
t35 = -t151 * mrSges(7,1) - t72 * mrSges(7,3);
t322 = -t361 / 0.2e1;
t118 = t162 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t78 = -t221 * t143 + t176 * t224;
t317 = pkin(2) * t375 - t372;
t314 = pkin(2) * t387 + qJ(3) * t392 + t380;
t27 = qJ(6) * t153 + t44;
t304 = mrSges(3,1) * t222 + mrSges(3,2) * t225;
t289 = -Ifges(5,2) * t221 + t420;
t288 = Ifges(5,2) * t224 + t421;
t278 = Ifges(7,5) * t224 + Ifges(7,6) * t221;
t277 = Ifges(7,5) * t221 - Ifges(7,6) * t224;
t273 = -pkin(4) * t221 + t403;
t30 = -pkin(4) * t189 + qJD(5) - t44;
t270 = t221 * t30 + t224 * t32;
t269 = t221 * t44 - t224 * t45;
t165 = -qJD(2) * pkin(2) + t362;
t169 = -t201 - t220;
t267 = t165 * t225 + t169 * t222;
t101 = qJD(2) * t272 + t317;
t160 = t456 * t373;
t24 = -t221 * t101 - t143 * t368 + t160 * t224 - t176 * t369;
t266 = t323 - t213;
t263 = t223 * pkin(3) + pkin(8) * t387 + t314;
t262 = pkin(1) * t304;
t259 = pkin(3) * t161 + t107;
t129 = t266 * qJD(1);
t257 = t129 * (-mrSges(4,2) * t222 - t408);
t256 = t222 * (Ifges(3,1) * t225 - t423);
t23 = t224 * t101 - t143 * t369 + t221 * t160 + t176 * t368;
t253 = Ifges(7,5) * t72 + Ifges(7,6) * t73 - Ifges(7,3) * t151;
t251 = t222 * t374 + t346;
t130 = -t222 * t389 + t396;
t132 = t221 * t226 + t222 * t391;
t248 = -g(1) * t130 + g(2) * t132 - g(3) * t390;
t17 = qJ(5) * t373 + t222 * qJD(5) + t23;
t242 = Ifges(5,6) * t225 + t222 * t288;
t238 = -Ifges(7,3) * t225 + t222 * t277;
t97 = -mrSges(7,1) * t189 - t424;
t237 = (t95 - t430) * t224 + (t97 + t501) * t221;
t235 = qJ(5) * t72 + qJD(5) * t153 + t259;
t233 = qJD(4) * t270 + t482;
t232 = -qJD(4) * t269 + t483;
t185 = qJ(3) * t387;
t184 = t223 * t406;
t167 = -pkin(1) - t382;
t166 = qJ(3) - t273;
t164 = t386 * t224;
t163 = t386 * t221;
t159 = t456 * t375;
t156 = -qJ(3) * t376 + t200;
t155 = t298 * qJD(1);
t142 = t303 * t225;
t141 = -qJ(3) - t468;
t133 = -t222 * t396 + t389;
t131 = t221 * t392 + t391;
t119 = -qJ(3) * t373 + t317;
t117 = mrSges(4,1) * t161 - qJDD(2) * mrSges(4,3);
t100 = t225 * t274 + t177;
t81 = pkin(4) * t153 + t405;
t80 = -t225 * t260 - t177;
t75 = -pkin(4) * t222 - t78;
t71 = pkin(2) * t161 + t247;
t51 = -pkin(4) * t376 - t68;
t50 = qJ(6) * t390 + t74;
t48 = -t153 * t455 - t405;
t46 = qJ(6) * t395 - t222 * t455 - t78;
t41 = (qJD(4) * t273 + t365) * t225 + (-t274 - t456) * t375;
t38 = mrSges(7,2) * t151 + mrSges(7,3) * t73;
t26 = (qJD(4) * t468 - t365) * t225 + (t260 + t456) * t375;
t22 = mrSges(5,1) * t73 + mrSges(5,2) * t72;
t20 = mrSges(6,1) * t73 - mrSges(6,3) * t72;
t19 = -t189 * t455 + qJD(5) - t27;
t18 = -pkin(4) * t373 - t24;
t10 = -qJ(6) * t251 + t225 * t363 + t17;
t9 = -qJ(6) * t347 + (qJ(6) * t368 - qJD(2) * t455 + t364) * t225 - t24;
t8 = pkin(4) * t73 - t235;
t4 = -t455 * t73 + qJDD(6) + t235;
t11 = [m(5) * (-t128 * t159 - t177 * t259 + t23 * t45 + t24 * t44 + t6 * t79 + t7 * t78) - t259 * t142 + t305 * t402 + (mrSges(5,2) * t128 + t30 * mrSges(6,2) + mrSges(7,2) * t29 - t44 * mrSges(5,3) - mrSges(6,3) * t47 - t19 * mrSges(7,3)) * (-t224 * t367 + t347) - (t221 * t504 + t224 * t475) * t367 / 0.2e1 + (-t494 * pkin(7) - t45 * mrSges(5,2) - t19 * mrSges(7,1) - t30 * mrSges(6,1) + t32 * mrSges(6,3) + t25 * mrSges(7,2) + t44 * mrSges(5,1) + t474 / 0.2e1 + t165 * mrSges(4,1) - t500 / 0.2e1) * t373 + (-Ifges(3,4) * t161 + Ifges(3,5) * qJDD(2) + t487) * t222 / 0.2e1 + t162 * t422 / 0.2e1 + t162 * t222 * Ifges(3,1) + (-mrSges(5,1) * t128 - mrSges(6,1) * t47 + mrSges(7,1) * t29 - t513) * t251 + t485 * mrSges(4,1) + (-t3 * mrSges(6,2) - t6 * mrSges(5,3) + t2 * mrSges(7,3) + t514) * t390 - t222 * t521 + t161 * t541 + t161 * t542 + t346 * t543 + m(6) * (t100 * t8 + t17 * t32 + t18 * t30 + t3 * t74 + t41 * t47 + t5 * t75) + m(7) * (t1 * t46 + t10 * t25 + t19 * t9 + t2 * t50 + t26 * t29 + t4 * t80) + (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t162 + Ifges(4,3) * t161) * t544 + (-qJDD(2) * mrSges(3,2) - t117) * t211 + qJD(2) * t257 + (Ifges(5,6) * t222 - t225 * t288) * t459 + (qJD(2) * t238 - t278 * t367) * t448 + (qJD(2) * t242 - t289 * t367) * t452 + (-Ifges(7,3) * t222 - t225 * t277) * t454 + (-qJDD(2) * mrSges(3,1) + t118) * t442 + t8 * t301 * t225 - t162 * t276 / 0.2e1 + t71 * t298 - t262 * t361 - t4 * t300 * t225 + Ifges(2,3) * qJDD(1) + (t256 + t225 * (-Ifges(3,2) * t222 + t422)) * t361 / 0.2e1 - (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t162 + Ifges(4,6) * t161 + t253) * t222 / 0.2e1 + (t169 * mrSges(4,1) - t495 * pkin(7) + t524) * t375 + (-m(5) * t381 + t511 * (t133 * pkin(4) + qJ(5) * t132 + t381) + (-m(3) - m(4)) * t214 + t478 * t226 - t465 * t133 - t477 * t132 + (m(3) * pkin(1) - m(4) * t266 + t525 * t225 + t479 * t323 - t473) * t223) * g(1) + (t222 * t527 - t471 * t225) * t460 + t46 * t35 + t50 * t38 + (t532 * t222 - t531 * t225) * t458 + (t535 * t222 - t530 * t225) * t453 + t74 * t40 + t75 * t37 + t78 * t36 + t79 * t39 + (qJD(2) * t472 - t367 * t470) * t449 + t80 * t21 + t41 * t82 + t26 * t83 + t5 * (-mrSges(6,1) * t222 - mrSges(6,2) * t395) + t1 * (-mrSges(7,1) * t222 + mrSges(7,3) * t395) + t7 * (mrSges(5,1) * t222 + mrSges(5,3) * t395) - t476 * t395 / 0.2e1 + t10 * t94 + t23 * t95 + t17 * t96 + t9 * t97 + t24 * t98 + t18 * t99 + t100 * t20 + m(4) * (t119 * t129 + t167 * t71 + (qJD(2) * t267 + t485) * pkin(7)) + (-t161 * t211 + t162 * t442 + t486) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t486) + t488 * qJD(2) ^ 2 / 0.2e1 + t491 * t322 + (qJD(2) * t492 - t367 * t490) * t451 + (qJD(2) * t493 - t367 * t489) * t447 + t119 * t155 - t159 * t84 - pkin(1) * (mrSges(3,1) * t161 + mrSges(3,2) * t162) + t167 * (-mrSges(4,2) * t161 - mrSges(4,3) * t162) + t177 * t22 + (t222 * t510 - t225 * t508) * qJDD(2) / 0.2e1 + (-m(3) * t380 - m(4) * t314 - m(5) * t263 + t511 * (t131 * pkin(4) + qJ(5) * t130 + t263) + t464 * t387 + t473 * t226 + t478 * t223 - t465 * t131 - t477 * t130) * g(2) + t225 * (Ifges(3,4) * t162 - Ifges(3,2) * t161 + Ifges(3,6) * qJDD(2)) / 0.2e1 + t222 * t522 + t222 * t523; -t259 * t302 + t383 * t84 + (-t502 * t366 - t506 * t457 + t476 / 0.2e1 - t8 * mrSges(6,3)) * t224 + (t8 * t166 - t233 * t457 - t30 * t51 - t32 * t49 + t47 * t496) * m(6) - (t197 + t474) * t376 / 0.2e1 + (t369 * t44 - t483) * mrSges(5,3) + (-t30 * t369 - t482) * mrSges(6,2) + (t19 * t369 + t481) * mrSges(7,3) + ((-t286 / 0.2e1 - t281 / 0.2e1 + t277 / 0.2e1) * t189 + t471 * t450 + t463) * qJD(4) + t500 * t538 + (-t57 / 0.2e1 + t504 / 0.2e1 + t513) * t368 + (t22 - t117) * qJ(3) + (-t501 * t366 - t505 * t457 + t514) * t221 + ((t238 / 0.2e1 - t493 / 0.2e1) * t189 + (t242 / 0.2e1 - t492 / 0.2e1) * t152 + t472 * t450 - m(4) * pkin(7) * t267 - t25 * (mrSges(7,2) * t225 - mrSges(7,3) * t393) - t45 * (-mrSges(5,2) * t225 + mrSges(5,3) * t393) - t32 * (mrSges(6,2) * t393 + mrSges(6,3) * t225) - t19 * (-mrSges(7,1) * t225 - mrSges(7,3) * t397) - t44 * (mrSges(5,1) * t225 - mrSges(5,3) * t397) - t30 * (-mrSges(6,1) * t225 + mrSges(6,2) * t397) - t257 + (t262 - t256 / 0.2e1 + t491 / 0.2e1) * qJD(1)) * qJD(1) + (t288 / 0.2e1 - t284 / 0.2e1 - t279 / 0.2e1) * t371 + t289 * t459 + t278 * t454 + t4 * t299 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (-pkin(2) * t113 - qJ(3) * t107 - qJD(3) * t169 - t129 * t156) * m(4) + (-m(4) * t382 - m(5) * t348 + t511 * (pkin(4) * t397 + t348) + t464 * t225 + t526 * t222 + t518) * g(3) + (-t259 * qJ(3) + t529 * t128 - t457 * t232 - t44 * t68 - t45 * t69) * m(5) + t8 * t410 + t470 * t460 - t475 * t369 / 0.2e1 - t165 * t349 - t169 * t350 + (Ifges(3,2) * t538 + t463 - t524) * t377 - t69 * t95 - t49 * t96 - t68 * t98 - t51 * t99 - t107 * mrSges(4,3) + t113 * mrSges(4,2) - pkin(2) * t118 + t480 * t304 + t488 * t322 + t489 * t453 + t490 * t458 + t494 * t201 + t495 * t198 + t141 * t21 + t149 * mrSges(3,2) - t150 * mrSges(3,1) + t496 * t82 + t497 * qJD(3) + t498 * t94 + t499 * t97 - t156 * t155 + t163 * t38 - t164 * t35 + t166 * t20 + t503 * t83 + (-t1 * t164 + t141 * t4 + t163 * t2 + t19 * t499 + t25 * t498 + t29 * t503) * m(7) + t508 * t161 + t510 * t162 + (t511 * (t226 * t356 + t185) + t512 * t185 + t467 * t226) * g(1) + (t511 * (t223 * t356 + t184) + t512 * t184 + t467 * t223) * g(2); (-t35 + t506) * t224 + (t38 + t505) * t221 + (-t431 - t497) * qJD(2) + t351 * t439 + t237 * qJD(4) + ((t155 + t237) * qJD(1) - t480 * t351) * t222 + t118 + (qJD(2) * t29 + t481 + t189 * (t19 * t221 + t224 * t25)) * m(7) + (-qJD(2) * t47 + t270 * t377 + t233) * m(6) + (-qJD(2) * t128 - t269 * t377 + t232) * m(5) + (qJD(2) * t169 + t129 * t377 + t113) * m(4); (qJ(5) * t2 - t1 * t455 - t19 * t28 - t25 * t27 - t29 * t48) * m(7) - t455 * t35 + t523 - t521 + t522 + t32 * t428 + t30 * t429 - t253 + (t38 + t40) * qJ(5) - t1 * mrSges(7,1) + (-pkin(4) * t5 + qJ(5) * t3 + qJD(5) * t32 + t274 * t439 - t47 * t81) * m(6) + (-Ifges(7,5) * t152 + Ifges(7,6) * t153) * t447 + t57 * t449 + (-t152 * t537 + t533 * t153) * t448 + (t534 * t153 - t516) * t452 + t487 + (t301 - (-m(7) * qJ(5) - mrSges(7,2)) * t221 - (-m(7) * t455 - mrSges(7,1)) * t224) * t439 - t5 * mrSges(6,1) + t7 * mrSges(5,1) + (-t528 * t152 - t411 + t504 + t519) * t450 - pkin(4) * t37 - t81 * t82 - t48 * t83 + (-Ifges(5,2) * t153 - t146 + t475) * t451 - t27 * t94 - t28 * t97 + g(3) * t142 - t29 * (-mrSges(7,1) * t153 - mrSges(7,2) * t152) - t47 * (mrSges(6,1) * t153 + mrSges(6,3) * t152) - t128 * (mrSges(5,1) * t153 - mrSges(5,2) * t152) + (-m(6) * t30 + t426 - t501) * t45 + (-m(6) * t32 - t427 - t502) * t44 + (m(7) * t25 - t430) * qJD(5) - t25 * t424 - t19 * t425 + (t511 * (t132 * pkin(4) - qJ(5) * t133) + t477 * t133 - t465 * t132) * g(2) + (t511 * (-t130 * pkin(4) + qJ(5) * t131) - t477 * t131 + t465 * t130) * g(1); t431 * t153 + t430 * t189 + t35 + t37 + (-t153 * t29 - t189 * t25 + t1 + t248) * m(7) + (t153 * t47 - t189 * t32 + t248 + t5) * m(6); -t152 * t94 + t153 * t97 + (g(3) * t222 - t25 * t152 + t19 * t153 + t225 * t480 + t4) * m(7) + t21;];
tau  = t11;
