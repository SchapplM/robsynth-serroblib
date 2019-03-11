% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:29
% EndTime: 2019-03-09 09:08:18
% DurationCPUTime: 38.51s
% Computational Cost: add. (8957->865), mult. (21337->1144), div. (0->0), fcn. (15389->10), ass. (0->378)
t548 = Ifges(5,4) + Ifges(4,5);
t550 = Ifges(4,1) + Ifges(5,1);
t285 = cos(qJ(2));
t583 = t548 * t285;
t278 = qJ(3) - pkin(9);
t558 = m(6) + m(7);
t582 = t278 * t558 - mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t549 = Ifges(4,4) + Ifges(3,5);
t517 = -Ifges(5,5) + t549;
t284 = cos(qJ(5));
t280 = sin(qJ(5));
t407 = qJD(5) * t280;
t277 = cos(pkin(6));
t475 = pkin(1) * t277;
t264 = t285 * t475;
t281 = sin(qJ(2));
t276 = sin(pkin(6));
t414 = qJD(1) * t276;
t384 = t281 * t414;
t181 = -pkin(8) * t384 + qJD(1) * t264;
t224 = qJ(4) * t384;
t141 = t224 + t181;
t383 = t285 * t414;
t227 = qJ(3) * t383;
t490 = -pkin(3) - pkin(4);
t275 = pkin(2) - t490;
t294 = t276 * (-pkin(9) * t285 - t275 * t281);
t97 = qJD(1) * t294 + t227;
t61 = t284 * t141 + t280 * t97;
t581 = qJD(3) * t284 - t278 * t407 - t61;
t580 = mrSges(3,1) + mrSges(4,1);
t547 = Ifges(5,2) + Ifges(4,3);
t579 = Ifges(3,6) - Ifges(4,6);
t546 = -Ifges(5,6) + Ifges(4,6);
t578 = Ifges(3,3) + Ifges(4,2);
t209 = pkin(5) * t284 + pkin(10) * t280 + t275;
t577 = -pkin(10) * t384 - qJD(6) * t209 - t581;
t226 = qJ(4) * t383;
t263 = t281 * t475;
t345 = -pkin(5) * t280 + pkin(10) * t284;
t431 = t278 * t284;
t434 = t276 * t285;
t576 = t345 * qJD(5) - qJD(6) * t431 - t226 - (-t263 + (-pkin(8) - t345) * t434) * qJD(1);
t399 = qJDD(1) * t277;
t250 = qJDD(2) + t399;
t401 = qJD(1) * qJD(2);
t186 = (qJDD(1) * t281 + t285 * t401) * t276;
t437 = t276 * t281;
t247 = qJD(4) * t437;
t397 = qJD(2) * t475;
t358 = qJD(1) * t397;
t394 = pkin(1) * t399;
t105 = -pkin(8) * t186 - t281 * t358 + t285 * t394;
t89 = -t250 * pkin(2) + qJDD(3) - t105;
t51 = -t250 * pkin(3) - t186 * qJ(4) - qJD(1) * t247 + t89;
t46 = -pkin(4) * t250 + t51;
t406 = qJD(5) * t284;
t60 = -t141 * t280 + t284 * t97;
t575 = -qJD(3) * t280 - t278 * t406 - t60;
t574 = (t281 * t550 - t583) * t276;
t400 = qJDD(1) * t276;
t246 = t285 * t400;
t412 = qJD(2) * t281;
t382 = t276 * t412;
t355 = qJD(1) * t382;
t185 = -t246 + t355;
t248 = qJD(3) * t437;
t266 = pkin(1) * t400;
t82 = t185 * pkin(2) - t186 * qJ(3) - qJD(1) * t248 - t266;
t320 = qJDD(4) - t82;
t38 = -pkin(9) * t186 + t185 * t490 + t320;
t386 = pkin(8) * t246 + t281 * t394 + t285 * t358;
t408 = qJD(4) * t285;
t255 = qJD(1) * t277 + qJD(2);
t531 = t250 * qJ(3) + t255 * qJD(3);
t54 = t185 * qJ(4) + (-pkin(8) * t412 - t408) * t414 + t386 + t531;
t48 = -pkin(9) * t250 + t54;
t147 = -pkin(1) * t414 - pkin(2) * t383 - qJ(3) * t384;
t113 = pkin(3) * t383 + qJD(4) - t147;
t87 = (pkin(4) * t285 - pkin(9) * t281) * t414 + t113;
t416 = pkin(8) * t434 + t263;
t182 = t416 * qJD(1);
t143 = -t226 + t182;
t223 = t255 * qJ(3);
t107 = t223 + t143;
t96 = -pkin(9) * t255 + t107;
t7 = t280 * t38 + t284 * t48 + t87 * t406 - t407 * t96;
t573 = t7 * mrSges(6,2);
t504 = t276 ^ 2;
t571 = t504 / 0.2e1;
t156 = t255 * t284 + t280 * t384;
t150 = Ifges(6,4) * t156;
t570 = t156 * Ifges(6,6);
t213 = qJD(5) + t383;
t569 = t213 * Ifges(6,5);
t568 = t213 * Ifges(6,3);
t566 = -pkin(5) * t384 - t575;
t565 = t548 * t384;
t279 = sin(qJ(6));
t283 = cos(qJ(6));
t337 = -mrSges(7,1) * t283 + mrSges(7,2) * t279;
t308 = m(7) * pkin(5) - t337;
t340 = mrSges(6,1) * t284 - mrSges(6,2) * t280;
t395 = m(7) * pkin(10) + mrSges(7,3);
t564 = -t280 * t395 - t284 * t308 - t340;
t43 = t280 * t87 + t284 * t96;
t37 = pkin(10) * t213 + t43;
t154 = -t255 * t280 + t284 * t384;
t131 = -t255 * pkin(2) + qJD(3) - t181;
t94 = -t255 * pkin(3) + t131 - t224;
t83 = pkin(4) * t255 - t94;
t44 = pkin(5) * t156 - pkin(10) * t154 + t83;
t13 = -t279 * t37 + t283 * t44;
t14 = t279 * t44 + t283 * t37;
t78 = -qJD(5) * t156 + t186 * t284 - t250 * t280;
t79 = -qJD(5) * t154 - t186 * t280 - t250 * t284;
t15 = -pkin(5) * t79 - pkin(10) * t78 - t46;
t170 = qJDD(5) - t185;
t5 = pkin(10) * t170 + t7;
t1 = qJD(6) * t13 + t15 * t279 + t283 * t5;
t2 = -qJD(6) * t14 + t15 * t283 - t279 * t5;
t344 = t1 * t283 - t2 * t279;
t403 = qJD(6) * t283;
t405 = qJD(6) * t279;
t563 = -t13 * t403 - t14 * t405 + t344;
t423 = t284 * t285;
t148 = (t279 * t423 + t281 * t283) * t414;
t562 = t279 * t406 + t148;
t151 = qJD(6) + t156;
t282 = sin(qJ(1));
t425 = t282 * t285;
t286 = cos(qJ(1));
t426 = t281 * t286;
t202 = t277 * t426 + t425;
t433 = t276 * t286;
t126 = t202 * t284 + t280 * t433;
t421 = t285 * t286;
t427 = t281 * t282;
t201 = -t277 * t421 + t427;
t560 = t126 * t279 + t201 * t283;
t559 = -t126 * t283 + t201 * t279;
t42 = -t280 * t96 + t284 * t87;
t36 = -pkin(5) * t213 - t42;
t557 = m(7) * t36;
t556 = t1 * mrSges(7,2);
t555 = t2 * mrSges(7,1);
t552 = t13 * mrSges(7,1);
t551 = t14 * mrSges(7,2);
t102 = -t154 * t279 + t213 * t283;
t32 = qJD(6) * t102 + t170 * t279 + t283 * t78;
t103 = t154 * t283 + t213 * t279;
t33 = -qJD(6) * t103 + t170 * t283 - t279 * t78;
t12 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t58 = mrSges(6,1) * t170 - mrSges(6,3) * t78;
t545 = t12 - t58;
t544 = t279 * t577 + t283 * t576;
t543 = t279 * t576 - t283 * t577;
t542 = t156 * Ifges(6,2);
t541 = t213 * Ifges(6,6);
t540 = mrSges(6,1) + t308;
t512 = mrSges(6,2) - t395;
t158 = t277 * qJ(3) + t416;
t132 = -qJ(4) * t434 + t158;
t114 = -pkin(9) * t277 + t132;
t270 = t276 * pkin(1);
t417 = pkin(2) * t434 + qJ(3) * t437;
t159 = -t270 - t417;
t259 = pkin(3) * t434;
t133 = t259 - t159;
t258 = pkin(4) * t434;
t99 = -pkin(9) * t437 + t133 + t258;
t539 = t284 * t114 + t280 * t99;
t462 = Ifges(3,4) * t281;
t312 = (t285 * Ifges(3,2) + t462) * t276;
t538 = t154 * Ifges(6,5) + t255 * Ifges(3,6) + qJD(1) * t312 + t568 - t570;
t463 = mrSges(6,3) * t154;
t109 = mrSges(6,1) * t213 - t463;
t55 = -mrSges(7,1) * t102 + mrSges(7,2) * t103;
t446 = -t55 + t109;
t360 = mrSges(5,3) * t384;
t537 = mrSges(5,1) * t255 + mrSges(6,1) * t156 + mrSges(6,2) * t154 + t360;
t536 = t255 * t546 - t383 * t547 + t565;
t535 = t280 * t403 + t562;
t424 = t283 * t284;
t149 = -t279 * t384 + t383 * t424;
t404 = qJD(6) * t280;
t534 = -t279 * t404 + t283 * t406 + t149;
t359 = mrSges(5,3) * t383;
t175 = mrSges(5,2) * t255 - t359;
t363 = mrSges(4,2) * t383;
t177 = mrSges(4,3) * t255 + t363;
t533 = t175 + t177;
t184 = t416 * qJD(2);
t411 = qJD(2) * t285;
t381 = t276 * t411;
t532 = -qJ(4) * t381 + t184 - t247;
t438 = t255 * t276;
t530 = (-t281 * t579 + t285 * t549) * t438;
t443 = qJD(5) * t43;
t8 = -t280 * t48 + t284 * t38 - t443;
t528 = -t185 * t579 + t186 * t549 + t578 * t250;
t76 = qJDD(6) - t79;
t18 = mrSges(7,1) * t76 - mrSges(7,3) * t32;
t19 = -mrSges(7,2) * t76 + mrSges(7,3) * t33;
t525 = -t279 * t18 + t283 * t19;
t524 = -t13 * t279 + t14 * t283;
t523 = t548 * t281;
t522 = -t280 * t8 + t284 * t7;
t520 = mrSges(5,1) + t580;
t519 = -mrSges(4,2) + mrSges(5,3) - mrSges(3,3);
t238 = Ifges(3,4) * t383;
t516 = Ifges(3,1) * t384 + t574 * qJD(1) + t255 * t517 + t238;
t514 = (Ifges(5,5) * t285 + Ifges(5,6) * t281) * t438 / 0.2e1 + (-t147 * (mrSges(4,1) * t281 - mrSges(4,3) * t285) - t113 * (-mrSges(5,1) * t281 + mrSges(5,2) * t285)) * t276;
t245 = t285 * t397;
t268 = t277 * qJD(3);
t101 = t245 + t268 + (-t408 + (-pkin(8) + qJ(4)) * t412) * t276;
t419 = qJ(3) * t381 + t248;
t84 = qJD(2) * t294 + t419;
t21 = -qJD(5) * t539 - t101 * t280 + t284 * t84;
t511 = t555 - t556;
t362 = mrSges(3,3) * t384;
t364 = mrSges(4,2) * t384;
t510 = -m(4) * t131 + t255 * t580 - t362 - t364;
t509 = t520 - t564;
t336 = t279 * mrSges(7,1) + t283 * mrSges(7,2);
t508 = t336 - t582;
t507 = -t551 + t552;
t506 = t281 * (Ifges(3,1) * t285 - t462) * t571 + (-pkin(1) * (mrSges(3,1) * t281 + mrSges(3,2) * t285) - (t547 * t281 + t583) * t285 / 0.2e1) * t504;
t39 = t103 * Ifges(7,5) + t102 * Ifges(7,6) + t151 * Ifges(7,3);
t505 = t43 * mrSges(6,3) - t39 / 0.2e1 - t507;
t9 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t76;
t502 = t9 / 0.2e1;
t10 = t32 * Ifges(7,4) + t33 * Ifges(7,2) + t76 * Ifges(7,6);
t501 = t10 / 0.2e1;
t11 = t32 * Ifges(7,1) + t33 * Ifges(7,4) + t76 * Ifges(7,5);
t500 = t11 / 0.2e1;
t499 = -t78 * Ifges(6,4) / 0.2e1 - t79 * Ifges(6,2) / 0.2e1 - t170 * Ifges(6,6) / 0.2e1;
t498 = t32 / 0.2e1;
t497 = t33 / 0.2e1;
t451 = t103 * Ifges(7,4);
t40 = t102 * Ifges(7,2) + t151 * Ifges(7,6) + t451;
t496 = t40 / 0.2e1;
t100 = Ifges(7,4) * t102;
t41 = Ifges(7,1) * t103 + Ifges(7,5) * t151 + t100;
t495 = -t41 / 0.2e1;
t458 = Ifges(6,4) * t154;
t72 = t458 + t541 - t542;
t494 = t72 / 0.2e1;
t493 = t76 / 0.2e1;
t492 = t78 / 0.2e1;
t491 = t79 / 0.2e1;
t287 = -pkin(2) - pkin(3);
t489 = -t102 / 0.2e1;
t488 = t102 / 0.2e1;
t487 = -t103 / 0.2e1;
t486 = t103 / 0.2e1;
t485 = -t151 / 0.2e1;
t484 = t151 / 0.2e1;
t481 = t154 / 0.2e1;
t480 = t170 / 0.2e1;
t6 = -pkin(5) * t170 - t8;
t471 = t280 * t6;
t464 = mrSges(6,3) * t156;
t461 = Ifges(3,4) * t285;
t457 = Ifges(6,4) * t280;
t456 = Ifges(6,4) * t284;
t455 = Ifges(7,4) * t279;
t454 = Ifges(7,4) * t283;
t445 = qJ(3) * t201;
t203 = t277 * t425 + t426;
t444 = qJ(3) * t203;
t442 = t156 * t279;
t441 = t156 * t283;
t436 = t276 * t282;
t435 = t276 * t284;
t430 = t279 * t280;
t429 = t280 * t283;
t428 = t280 * t285;
t418 = mrSges(5,1) * t434 + mrSges(5,2) * t437;
t206 = -pkin(8) * t437 + t264;
t415 = t286 * pkin(1) + pkin(8) * t436;
t402 = m(5) + t558;
t398 = Ifges(6,5) * t78 + Ifges(6,6) * t79 + Ifges(6,3) * t170;
t204 = -t277 * t427 + t421;
t387 = t204 * pkin(2) + t415;
t160 = -t277 * pkin(2) - t206;
t385 = t259 + t417;
t377 = t437 / 0.2e1;
t370 = t403 / 0.2e1;
t369 = -pkin(1) * t282 + pkin(8) * t433;
t119 = -t250 * mrSges(4,1) + t186 * mrSges(4,2);
t368 = -t185 * mrSges(5,1) + t186 * mrSges(5,2);
t118 = -t250 * mrSges(5,1) - t186 * mrSges(5,3);
t189 = t201 * pkin(2);
t367 = qJ(3) * t202 - t189;
t195 = t203 * pkin(2);
t366 = qJ(3) * t204 - t195;
t365 = t287 * t437;
t361 = mrSges(3,3) * t383;
t348 = -t202 * pkin(2) + t369;
t35 = -t79 * mrSges(6,1) + t78 * mrSges(6,2);
t342 = mrSges(4,1) * t285 + mrSges(4,3) * t281;
t199 = t277 * t284 + t280 * t437;
t200 = -t277 * t280 + t281 * t435;
t341 = mrSges(6,1) * t199 + mrSges(6,2) * t200;
t123 = -t200 * t279 + t283 * t434;
t124 = t200 * t283 + t279 * t434;
t338 = mrSges(7,1) * t123 - mrSges(7,2) * t124;
t335 = Ifges(6,1) * t284 - t457;
t334 = Ifges(7,1) * t283 - t455;
t333 = Ifges(7,1) * t279 + t454;
t332 = -Ifges(6,2) * t280 + t456;
t331 = -Ifges(7,2) * t279 + t454;
t330 = Ifges(7,2) * t283 + t455;
t329 = Ifges(6,5) * t284 - Ifges(6,6) * t280;
t328 = Ifges(7,5) * t283 - Ifges(7,6) * t279;
t327 = Ifges(7,5) * t279 + Ifges(7,6) * t283;
t50 = pkin(10) * t434 + t539;
t115 = -t277 * pkin(3) - qJ(4) * t437 + t160;
t106 = t277 * pkin(4) - t115;
t65 = pkin(5) * t199 - pkin(10) * t200 + t106;
t23 = t279 * t65 + t283 * t50;
t22 = -t279 * t50 + t283 * t65;
t66 = -mrSges(7,2) * t151 + mrSges(7,3) * t102;
t67 = mrSges(7,1) * t151 - mrSges(7,3) * t103;
t325 = t279 * t67 - t283 * t66;
t56 = -t114 * t280 + t284 * t99;
t183 = -pkin(8) * t382 + t245;
t108 = -mrSges(6,2) * t213 - t464;
t318 = -t108 + t325;
t317 = -t202 * t280 + t284 * t433;
t20 = t284 * t101 - t114 * t407 + t280 * t84 + t99 * t406;
t315 = t342 * t276;
t306 = Ifges(5,5) * t186 + Ifges(5,6) * t185 - Ifges(5,3) * t250;
t295 = t204 * pkin(3) - qJ(4) * t436 + t387;
t104 = -pkin(8) * t355 + t386;
t292 = t204 * pkin(4) + t295;
t291 = -t202 * pkin(3) - qJ(4) * t433 + t348;
t289 = -t202 * pkin(4) + t291;
t205 = (-mrSges(3,1) * t285 + mrSges(3,2) * t281) * t276;
t194 = t203 * pkin(3);
t188 = t201 * pkin(3);
t180 = pkin(2) * t384 - t227;
t179 = (mrSges(5,1) * t285 + mrSges(5,2) * t281) * t414;
t178 = qJD(1) * t315;
t176 = -mrSges(3,2) * t255 + t361;
t152 = t183 + t268;
t146 = t209 * t279 + t278 * t424;
t145 = t209 * t283 - t279 * t431;
t144 = pkin(2) * t382 - t419;
t142 = qJD(1) * t365 + t227;
t140 = t223 + t182;
t130 = t204 * t284 - t280 * t436;
t129 = t204 * t280 + t282 * t435;
t122 = -qJD(5) * t199 + t284 * t381;
t121 = -t277 * t407 + (t280 * t411 + t281 * t406) * t276;
t120 = -mrSges(4,2) * t185 + mrSges(4,3) * t250;
t117 = mrSges(5,2) * t250 + mrSges(5,3) * t185;
t111 = qJD(2) * t365 + t419;
t92 = pkin(5) * t154 + pkin(10) * t156;
t86 = t130 * t283 - t203 * t279;
t85 = -t130 * t279 - t203 * t283;
t77 = t104 + t531;
t73 = t154 * Ifges(6,1) - t150 + t569;
t64 = qJD(6) * t123 + t122 * t283 - t279 * t382;
t63 = -qJD(6) * t124 - t122 * t279 - t283 * t382;
t62 = -pkin(3) * t185 + t320;
t59 = -mrSges(6,2) * t170 + mrSges(6,3) * t79;
t49 = -pkin(5) * t434 - t56;
t47 = pkin(5) * t121 - pkin(10) * t122 - t532;
t29 = t78 * Ifges(6,1) + t79 * Ifges(6,4) + t170 * Ifges(6,5);
t25 = t279 * t92 + t283 * t42;
t24 = -t279 * t42 + t283 * t92;
t17 = pkin(5) * t382 - t21;
t16 = -pkin(10) * t382 + t20;
t4 = -qJD(6) * t23 - t16 * t279 + t283 * t47;
t3 = qJD(6) * t22 + t16 * t283 + t279 * t47;
t26 = [(-m(6) * t292 - t130 * mrSges(6,1) - m(4) * (t387 + t444) - m(5) * (t295 + t444) - m(3) * t415 - m(7) * (pkin(5) * t130 + t292) - t86 * mrSges(7,1) - t85 * mrSges(7,2) - mrSges(2,1) * t286 + mrSges(2,2) * t282 + t512 * t129 + t519 * t436 - t520 * t204 - t582 * t203) * g(2) + (-m(5) * (t291 - t445) - m(4) * (t348 - t445) - m(3) * t369 - m(6) * t289 + t126 * mrSges(6,1) - m(7) * (-pkin(5) * t126 + t289) - t559 * mrSges(7,1) - t560 * mrSges(7,2) + mrSges(2,1) * t282 + mrSges(2,2) * t286 + t512 * t317 + t519 * t433 + t520 * t202 + t582 * t201) * g(1) + m(6) * (-t106 * t46 + t20 * t43 + t21 * t42 + t539 * t7 + t56 * t8) + t539 * t59 + m(4) * (t140 * t152 + t144 * t147 + t158 * t77 + t159 * t82 + t160 * t89) + ((t281 * t536 + t285 * t516) * t276 + t530) * qJD(2) / 0.2e1 + ((t281 * Ifges(3,1) + t461) * t276 + t517 * t277 + t574) * t186 / 0.2e1 + (t1 * t123 - t124 * t2 - t13 * t64 + t14 * t63) * mrSges(7,3) + (Ifges(7,1) * t64 + Ifges(7,4) * t63) * t486 + (Ifges(7,1) * t124 + Ifges(7,4) * t123) * t498 + (t131 * mrSges(4,2) - t181 * mrSges(3,3) - t94 * mrSges(5,3)) * t381 + (Ifges(6,4) * t200 + Ifges(6,6) * t434) * t491 + (Ifges(3,4) * t186 - Ifges(3,2) * t185 + Ifges(3,6) * t250 + t398) * t434 / 0.2e1 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t504 + t104 * t416 + t105 * t206 + t182 * t183) + t416 * (-mrSges(3,2) * t250 - mrSges(3,3) * t185) + (t546 * t277 + (-t285 * t547 + t523) * t276) * t185 / 0.2e1 - (t185 * t547 + t186 * t548 + t250 * t546) * t434 / 0.2e1 + (t517 * t250 + (Ifges(3,1) + t550) * t186 + (-Ifges(3,4) + t548) * t185) * t377 - t277 * t306 / 0.2e1 + t62 * t418 + (-m(3) * t181 - t510) * t184 + (-mrSges(6,3) * t7 - Ifges(6,4) * t492 + Ifges(7,5) * t498 - Ifges(6,2) * t491 - Ifges(6,6) * t480 + Ifges(7,6) * t497 + Ifges(7,3) * t493 + t499 + t502 + t511) * t199 + (m(5) * t94 - m(6) * t83 - t537) * t532 + t63 * t496 + t124 * t500 + t123 * t501 + (Ifges(7,5) * t64 + Ifges(7,6) * t63) * t484 + (Ifges(7,5) * t124 + Ifges(7,6) * t123) * t493 + t528 * t277 / 0.2e1 + m(5) * (t101 * t107 + t111 * t113 + t115 * t51 + t132 * t54 + t133 * t62) + t17 * t55 + t56 * t58 + t49 * t12 + t22 * t18 + t23 * t19 + (-Ifges(6,4) * t481 + Ifges(7,3) * t484 + Ifges(7,5) * t486 + Ifges(7,6) * t488 - t72 / 0.2e1 + t83 * mrSges(6,1) - t541 / 0.2e1 + t542 / 0.2e1 - t505) * t121 + (Ifges(6,5) * t200 + Ifges(6,3) * t434) * t480 - t514 * qJD(2) - t205 * t266 + t51 * (-mrSges(5,1) * t277 - mrSges(5,3) * t437) + t105 * (mrSges(3,1) * t277 - mrSges(3,3) * t437) + t89 * (-mrSges(4,1) * t277 + mrSges(4,2) * t437) + t104 * (-mrSges(3,2) * t277 + mrSges(3,3) * t434) + t77 * (mrSges(4,2) * t434 + mrSges(4,3) * t277) + t8 * (mrSges(6,1) * t434 - mrSges(6,3) * t200) + t54 * (mrSges(5,2) * t277 - mrSges(5,3) * t434) - t185 * (Ifges(3,6) * t277 + t312) / 0.2e1 - (mrSges(3,1) * t185 + mrSges(3,2) * t186) * t270 - t6 * t338 - t46 * t341 + (Ifges(7,4) * t64 + Ifges(7,2) * t63) * t488 + (Ifges(7,4) * t124 + Ifges(7,2) * t123) * t497 + ((t285 * (-Ifges(3,2) * t281 + t461) + (t285 * t550 + t523) * t281) * t571 + t506) * t401 + Ifges(2,3) * qJDD(1) + m(7) * (t1 * t23 + t13 * t4 + t14 * t3 + t17 * t36 + t2 * t22 + t49 * t6) + t206 * (mrSges(3,1) * t250 - mrSges(3,3) * t186) + t64 * t41 / 0.2e1 + t36 * (-mrSges(7,1) * t63 + mrSges(7,2) * t64) + t3 * t66 + t4 * t67 + (t578 * t277 + (t549 * t281 + t285 * t579) * t276) * t250 / 0.2e1 + (Ifges(6,1) * t200 + Ifges(6,5) * t434) * t492 + t106 * t35 + t20 * t108 + t21 * t109 + t115 * t118 + (t570 / 0.2e1 - t538 / 0.2e1 - t568 / 0.2e1 - t182 * mrSges(3,3) - t140 * mrSges(4,2) + t107 * mrSges(5,3) + t43 * mrSges(6,2) - Ifges(6,5) * t481 - t42 * mrSges(6,1)) * t382 + t132 * t117 - t82 * t315 + (-t150 / 0.2e1 + t569 / 0.2e1 + t83 * mrSges(6,2) + Ifges(6,1) * t481 + t73 / 0.2e1 - t42 * mrSges(6,3)) * t122 + t158 * t120 + t160 * t119 + t101 * t175 + t152 * t177 - t144 * t178 + t111 * t179 + t183 * t176 + t159 * (mrSges(4,1) * t185 - mrSges(4,3) * t186) - t434 * t573 + t200 * t29 / 0.2e1 - t250 * (-Ifges(5,3) * t277 + (Ifges(5,5) * t281 - Ifges(5,6) * t285) * t276) / 0.2e1 + t133 * t368; (t117 + t120) * qJ(3) + (-t42 * (-mrSges(6,1) * t281 - mrSges(6,3) * t423) - t43 * (mrSges(6,2) * t281 - mrSges(6,3) * t428)) * t414 + (t54 * qJ(3) - t113 * t142 - t143 * t94 + t287 * t51 + (qJD(3) - t141) * t107) * m(5) + t528 + t284 * t555 + (-pkin(2) * t89 + qJ(3) * t77 + qJD(3) * t140 - t147 * t180) * m(4) + (-m(4) * t140 - t176 - t177 + t361) * t181 - (t154 * t335 - t156 * t332 + t213 * t329) * qJD(5) / 0.2e1 + t575 * t109 - (t283 * t41 + t73) * t406 / 0.2e1 + (qJD(6) * t41 + t10) * t430 / 0.2e1 + (t494 + t505) * t407 + (t362 + t510) * t182 - t284 * t556 + t537 * t143 + (-Ifges(6,5) * t280 - Ifges(6,6) * t284) * t480 + (t327 * t404 + (-Ifges(7,3) * t280 - t284 * t328) * qJD(5)) * t484 + (t333 * t404 + (-Ifges(7,5) * t280 - t284 * t334) * qJD(5)) * t486 + (-m(5) * (-t194 + t366) - m(4) * t366 - t558 * (-t203 * pkin(4) - t194 - t195) + t508 * t204 + t509 * t203) * g(1) + (-m(5) * (-t188 + t367) - m(4) * t367 - t558 * (-t201 * pkin(4) - t188 - t189) + t508 * t202 + t509 * t201) * g(2) + (t330 * t404 + (-Ifges(7,6) * t280 - t284 * t331) * qJD(5)) * t488 + (-Ifges(6,2) * t284 - t457) * t491 + (-Ifges(6,1) * t280 - t456) * t492 + (Ifges(7,3) * t284 - t280 * t328) * t493 + t149 * t495 + (Ifges(7,6) * t284 - t280 * t331) * t497 + (Ifges(7,5) * t284 - t280 * t334) * t498 + t284 * t499 + t284 * t502 + t533 * qJD(3) + (t1 * t430 + t13 * t534 + t14 * t535 + t2 * t429) * mrSges(7,3) + (-mrSges(7,1) * t535 - mrSges(7,2) * t534) * t36 + t59 * t431 + (-t530 / 0.2e1 + t514 + t538 * t377 - t506 * qJD(1)) * qJD(1) - t107 * t360 + t54 * mrSges(5,2) - t51 * mrSges(5,1) + (Ifges(7,1) * t149 - Ifges(7,4) * t148) * t487 + (Ifges(7,4) * t149 - Ifges(7,2) * t148) * t489 + (t143 * t83 - t42 * t60 - t43 * t61 - t275 * t46 + (-t280 * t42 + t284 * t43) * qJD(3) + ((-t280 * t43 - t284 * t42) * qJD(5) + t522) * t278) * m(6) + (t406 * t42 - t522) * mrSges(6,3) - t306 - t11 * t429 / 0.2e1 - t336 * t471 - t46 * t340 + t543 * t66 + t544 * t67 + t545 * t278 * t280 - t213 * t83 * (mrSges(6,1) * t280 + mrSges(6,2) * t284) - t131 * t363 + t287 * t118 - t280 * t29 / 0.2e1 + t275 * t35 + t562 * t496 + (-m(4) * t417 - m(5) * t385 + t205 - t418 - t558 * (t258 + t385) + (-t342 + t564 * t285 + (pkin(9) * t558 + mrSges(6,3) + t336) * t281) * t276) * g(3) - (t213 * (-Ifges(6,3) * t281 + t285 * t329) + t154 * (-Ifges(6,5) * t281 + t285 * t335) - t156 * (-Ifges(6,6) * t281 + t285 * t332) + (t383 * t550 + t536 + t565) * t281 + (-Ifges(3,2) * t384 + t280 * t39 + t284 * t73 + t238 + t516) * t285) * t414 / 0.2e1 + t77 * mrSges(4,3) - t89 * mrSges(4,1) + t280 * t40 * t370 + (Ifges(7,5) * t487 + Ifges(7,6) * t489 + Ifges(7,3) * t485 + t494 - t507) * t280 * t383 + (Ifges(7,5) * t149 - Ifges(7,6) * t148) * t485 - t104 * mrSges(3,2) + t105 * mrSges(3,1) + t566 * t55 + (t1 * t146 + t544 * t13 + t543 * t14 + t145 * t2 + t471 * t278 + t36 * t566) * m(7) - pkin(2) * t119 + t145 * t18 + t146 * t19 - t141 * t175 - t142 * t179 + t180 * t178 + t94 * t359 + t140 * t364 + t581 * t108; t318 * t156 - t35 + t325 * qJD(6) + (-t178 - t179) * t384 + t118 - t533 * t255 - t446 * t154 - t283 * t18 - t279 * t19 + t119 + (-g(1) * t203 - g(2) * t201 + g(3) * t434) * (m(4) + t402) + (-t1 * t279 - t151 * t524 + t154 * t36 - t2 * t283) * m(7) + (-t154 * t42 - t156 * t43 + t46) * m(6) + (-t107 * t255 - t113 * t384 + t51) * m(5) + (-t140 * t255 + t147 * t384 + t89) * m(4); -t148 * t67 + t149 * t66 + m(5) * t62 - m(7) * (t13 * t148 - t14 * t149) + t402 * t277 * g(3) + (m(6) * (t8 + t443) - m(7) * t6 + (m(7) * t524 - t318) * qJD(5) - t545) * t284 + (t59 + (-t279 * t66 - t283 * t67) * qJD(6) - t446 * qJD(5) + m(6) * (-qJD(5) * t42 + t7) + m(7) * (qJD(5) * t36 + t563) + t525) * t280 + ((-t537 * t281 + (t284 * t108 - t280 * t446 + t175) * t285 + t428 * t557 - m(5) * (-t107 * t285 - t281 * t94) - m(6) * (t281 * t83 + t42 * t428 - t423 * t43)) * qJD(1) + (g(1) * t282 - g(2) * t286) * t402) * t276 + t368; -(-Ifges(6,1) * t156 + t39 - t458) * t154 / 0.2e1 + (Ifges(7,3) * t154 - t156 * t328) * t485 + t151 * t36 * t336 + t154 * t551 + (Ifges(7,5) * t154 - t156 * t334) * t487 + (Ifges(7,6) * t154 - t156 * t331) * t489 - t213 * (-Ifges(6,5) * t156 - Ifges(6,6) * t154) / 0.2e1 - t83 * (mrSges(6,1) * t154 - mrSges(6,2) * t156) + (-Ifges(6,2) * t154 - t150 + t73) * t156 / 0.2e1 + (t102 * t331 + t103 * t334 + t151 * t328) * qJD(6) / 0.2e1 + (-t464 - t108) * t42 - t40 * t405 / 0.2e1 + (t446 + t463 - t557) * t43 + t41 * t370 + t72 * t481 + t327 * t493 - t441 * t495 - t442 * t496 + t330 * t497 + t333 * t498 + t279 * t500 + t283 * t501 + (t129 * t540 + t130 * t512) * g(1) + (t126 * t512 - t317 * t540) * g(2) - t573 + (t199 * t308 - t200 * t395 + t341) * g(3) + (m(7) * ((-t13 * t283 - t14 * t279) * qJD(6) + t344) - t67 * t403 - t66 * t405 + t525) * pkin(10) + t398 - pkin(5) * t12 + t8 * mrSges(6,1) - t154 * t552 + (-pkin(5) * t6 - t13 * t24 - t14 * t25) * m(7) + t6 * t337 - t25 * t66 - t24 * t67 + (-t13 * t441 - t14 * t442 + t563) * mrSges(7,3); -t36 * (mrSges(7,1) * t103 + mrSges(7,2) * t102) + (Ifges(7,1) * t102 - t451) * t487 + t40 * t486 + (Ifges(7,5) * t102 - Ifges(7,6) * t103) * t485 - t13 * t66 + t14 * t67 - g(1) * (mrSges(7,1) * t85 - mrSges(7,2) * t86) - g(2) * (-mrSges(7,1) * t560 + mrSges(7,2) * t559) - g(3) * t338 + (t102 * t13 + t103 * t14) * mrSges(7,3) + t9 + (-Ifges(7,2) * t103 + t100 + t41) * t489 + t511;];
tau  = t26;
