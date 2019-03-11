% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:36
% EndTime: 2019-03-09 08:34:07
% DurationCPUTime: 21.84s
% Computational Cost: add. (4207->663), mult. (8391->804), div. (0->0), fcn. (4415->6), ass. (0->302)
t482 = Ifges(6,1) + Ifges(7,1);
t455 = Ifges(6,5) + Ifges(7,5);
t479 = Ifges(6,2) + Ifges(7,2);
t453 = Ifges(6,6) + Ifges(7,6);
t486 = -mrSges(6,3) - mrSges(7,3);
t480 = Ifges(6,4) + Ifges(7,4);
t199 = sin(qJ(5));
t202 = cos(qJ(5));
t332 = qJD(2) * t202;
t203 = cos(qJ(2));
t334 = qJD(1) * t203;
t125 = t199 * t334 - t332;
t485 = -t125 / 0.2e1;
t200 = sin(qJ(2));
t335 = qJD(1) * t200;
t155 = qJD(5) + t335;
t484 = -t155 / 0.2e1;
t483 = Ifges(5,1) + Ifges(4,3);
t481 = Ifges(4,4) + Ifges(3,5);
t454 = Ifges(3,6) - Ifges(4,6);
t478 = Ifges(4,6) - Ifges(5,5);
t477 = Ifges(6,3) + Ifges(7,3);
t476 = -Ifges(5,5) - t454;
t166 = t202 * pkin(5) + pkin(4);
t260 = mrSges(7,1) * t202 - mrSges(7,2) * t199;
t262 = mrSges(6,1) * t202 - mrSges(6,2) * t199;
t475 = -m(7) * t166 - t260 - t262;
t436 = -t453 * t199 + t455 * t202;
t365 = Ifges(7,4) * t202;
t367 = Ifges(6,4) * t202;
t435 = -t479 * t199 + t365 + t367;
t366 = Ifges(7,4) * t199;
t368 = Ifges(6,4) * t199;
t433 = t482 * t202 - t366 - t368;
t264 = t203 * mrSges(4,1) + t200 * mrSges(4,3);
t266 = mrSges(3,1) * t203 - mrSges(3,2) * t200;
t474 = -t264 - t266;
t409 = -pkin(3) - pkin(8);
t190 = -pkin(2) + t409;
t197 = -qJ(6) - pkin(8);
t205 = -pkin(2) - pkin(3);
t473 = -m(6) * t190 - m(5) * t205 - mrSges(5,2) - m(7) * (t197 + t205) - t486;
t323 = qJD(1) * qJD(2);
t299 = t200 * t323;
t321 = t203 * qJDD(1);
t135 = t299 - t321;
t136 = qJDD(1) * t200 + t203 * t323;
t176 = t200 * qJD(3);
t194 = qJDD(1) * pkin(1);
t56 = t135 * pkin(2) - t136 * qJ(3) - qJD(1) * t176 - t194;
t234 = qJDD(4) - t56;
t17 = pkin(4) * t136 + t135 * t409 + t234;
t328 = qJD(5) * t202;
t329 = qJD(5) * t199;
t119 = t136 * pkin(7);
t286 = qJDD(3) + t119;
t330 = qJD(4) * t200;
t209 = -qJ(4) * t136 - qJD(1) * t330 + t286;
t36 = qJDD(2) * t190 + t209;
t268 = t200 * pkin(4) + t203 * pkin(8);
t105 = -qJD(1) * pkin(1) - pkin(2) * t334 - qJ(3) * t335;
t79 = pkin(3) * t334 + qJD(4) - t105;
t62 = qJD(1) * t268 + t79;
t159 = qJ(4) * t335;
t171 = pkin(7) * t335;
t325 = qJD(3) + t171;
t301 = -t159 + t325;
t77 = qJD(2) * t190 + t301;
t3 = t199 * t17 + t202 * t36 + t62 * t328 - t329 * t77;
t126 = qJD(2) * t199 + t202 * t334;
t426 = qJD(5) * t126;
t59 = -qJDD(2) * t202 - t135 * t199 + t426;
t2 = qJ(6) * t59 + qJD(6) * t125 + t3;
t472 = t2 * mrSges(7,2);
t471 = t3 * mrSges(6,2);
t23 = t199 * t62 + t202 * t77;
t4 = -qJD(5) * t23 + t202 * t17 - t199 * t36;
t470 = t4 * mrSges(6,1);
t461 = t203 / 0.2e1;
t123 = qJDD(5) + t136;
t58 = qJD(5) * t125 - qJDD(2) * t199 + t135 * t202;
t469 = -t480 * t59 / 0.2e1 - t482 * t58 / 0.2e1 - t455 * t123 / 0.2e1;
t201 = sin(qJ(1));
t468 = g(2) * t201;
t467 = Ifges(5,6) + t481;
t466 = t480 * t125;
t465 = t480 * t126;
t169 = Ifges(4,5) * t335;
t370 = Ifges(5,4) * t200;
t257 = -t203 * Ifges(5,1) - t370;
t450 = -t126 * t482 + t155 * t455 + t466;
t464 = -Ifges(4,3) * t334 + qJD(1) * t257 + qJD(2) * t478 + t202 * t450 + t169;
t372 = Ifges(3,4) * t200;
t252 = t203 * Ifges(3,2) + t372;
t451 = t125 * t479 + t155 * t453 - t465;
t463 = Ifges(3,6) * qJD(2) + qJD(1) * t252 + t199 * t451;
t179 = t200 * mrSges(5,1);
t441 = -t203 * mrSges(5,2) + t179;
t462 = t203 * t486 - t441 + t474;
t412 = m(7) * pkin(5);
t460 = t123 * t453 + t479 * t59 + t480 * t58;
t22 = -t199 * t77 + t202 * t62;
t238 = t23 * t199 + t22 * t202;
t429 = -t199 * t4 + t202 * t3;
t459 = m(6) * (-qJD(5) * t238 + t429);
t324 = m(5) + m(6) + m(7);
t312 = m(4) + t324;
t458 = -mrSges(7,1) - mrSges(6,1);
t457 = mrSges(6,2) + mrSges(7,2);
t339 = qJ(6) - t190;
t278 = qJD(5) * t339;
t305 = t199 * t335;
t172 = pkin(7) * t334;
t131 = -qJ(4) * t334 + t172;
t161 = qJ(3) * t334;
t224 = pkin(4) * t203 + t190 * t200;
t68 = qJD(1) * t224 + t161;
t33 = t202 * t131 + t199 * t68;
t449 = qJ(6) * t305 - qJD(6) * t202 + t199 * t278 - t33;
t346 = t200 * t202;
t233 = pkin(5) * t203 - qJ(6) * t346;
t32 = -t131 * t199 + t202 * t68;
t448 = -qJD(1) * t233 + qJD(6) * t199 + t202 * t278 - t32;
t447 = mrSges(7,1) + t412;
t382 = pkin(7) - qJ(4);
t146 = t382 * t200;
t127 = t202 * t146;
t177 = t200 * qJ(3);
t185 = t203 * pkin(2);
t337 = t185 + t177;
t307 = t203 * pkin(3) + t337;
t237 = t307 + t268;
t76 = pkin(1) + t237;
t41 = t199 * t76 + t127;
t309 = mrSges(5,3) * t334;
t445 = qJD(2) * mrSges(5,1) - mrSges(6,1) * t125 - mrSges(6,2) * t126 - t309;
t128 = t171 - t159;
t444 = qJD(3) + t128 + (-t305 - t329) * pkin(5);
t153 = qJ(3) + t166;
t311 = mrSges(4,2) * t335;
t443 = -mrSges(3,3) * t335 - t311 + (mrSges(3,1) + mrSges(4,1)) * qJD(2);
t310 = mrSges(4,2) * t334;
t145 = qJD(2) * mrSges(4,3) + t310;
t442 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t334 + t145;
t440 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t439 = t200 * t436 + t203 * t477;
t438 = t200 * t435 + t203 * t453;
t437 = t200 * t433 + t203 * t455;
t434 = -t200 * t454 + t203 * t481;
t432 = t477 * t123 + t453 * t59 + t455 * t58;
t168 = pkin(7) * t321;
t118 = -pkin(7) * t299 + t168;
t431 = t118 * t203 + t119 * t200;
t80 = t118 + t440;
t87 = -qJDD(2) * pkin(2) + t286;
t430 = t200 * t87 + t203 * t80;
t1 = pkin(5) * t123 - qJ(6) * t58 + qJD(6) * t126 + t4;
t428 = -t1 * t199 + t2 * t202;
t204 = cos(qJ(1));
t427 = g(1) * t204 + t468;
t425 = mrSges(6,1) + t447;
t363 = Ifges(4,5) * t203;
t369 = Ifges(5,4) * t203;
t424 = t200 * t370 + (t363 - t369 + (-Ifges(5,2) + t483) * t200) * t203;
t420 = mrSges(5,3) - mrSges(3,3) - mrSges(4,2) + mrSges(2,2);
t170 = Ifges(3,4) * t334;
t258 = t200 * Ifges(4,1) - t363;
t419 = Ifges(3,1) * t335 + qJD(1) * t258 + qJD(2) * t481 + t125 * t453 - t126 * t455 + t155 * t477 + t170;
t333 = qJD(2) * t200;
t314 = pkin(7) * t333;
t231 = qJD(4) * t203 + t314;
t39 = -t135 * qJ(4) + qJD(1) * t231 - t168 - t440;
t357 = t203 * mrSges(4,3);
t359 = t203 * mrSges(5,1);
t418 = -t359 - t357 + (-m(6) * pkin(4) + t475) * t203 + (m(4) * pkin(2) + mrSges(4,1) + t473) * t200;
t417 = qJ(3) * t312;
t192 = qJD(2) * qJ(3);
t97 = -t131 - t192;
t86 = qJD(2) * pkin(4) - t97;
t415 = m(5) * t97 - m(6) * t86 - t445;
t405 = -t126 / 0.2e1;
t404 = t126 / 0.2e1;
t396 = -t203 / 0.2e1;
t394 = pkin(7) * t200;
t393 = pkin(7) * t203;
t390 = g(3) * t203;
t198 = qJ(3) + pkin(4);
t26 = mrSges(7,1) * t123 - mrSges(7,3) * t58;
t27 = mrSges(6,1) * t123 - mrSges(6,3) * t58;
t381 = t26 + t27;
t28 = -mrSges(7,2) * t123 + mrSges(7,3) * t59;
t29 = -mrSges(6,2) * t123 + mrSges(6,3) * t59;
t380 = t28 + t29;
t374 = mrSges(7,3) * t125;
t72 = -mrSges(7,2) * t155 + t374;
t376 = mrSges(6,3) * t125;
t73 = -mrSges(6,2) * t155 + t376;
t379 = t72 + t73;
t373 = mrSges(7,3) * t126;
t74 = mrSges(7,1) * t155 + t373;
t375 = mrSges(6,3) * t126;
t75 = mrSges(6,1) * t155 + t375;
t378 = -t74 - t75;
t377 = mrSges(7,2) * t202;
t371 = Ifges(3,4) * t203;
t364 = Ifges(4,5) * t200;
t350 = t197 * t203;
t349 = t199 * t200;
t348 = t199 * t203;
t347 = t199 * t204;
t345 = t200 * t204;
t344 = t201 * t199;
t343 = t201 * t202;
t342 = t202 * t203;
t341 = t203 * t204;
t340 = t204 * t202;
t91 = qJDD(2) * mrSges(5,2) - t136 * mrSges(5,3);
t331 = qJD(2) * t203;
t338 = qJ(3) * t331 + t176;
t336 = t204 * pkin(1) + t201 * pkin(7);
t327 = qJD(5) * t203;
t326 = qJD(6) * t203;
t61 = qJD(2) * t224 + t338;
t95 = t331 * t382 - t330;
t318 = t199 * t61 + t202 * t95 + t76 * t328;
t66 = -mrSges(7,1) * t125 - mrSges(7,2) * t126;
t317 = t66 + t445;
t308 = mrSges(5,3) * t335;
t306 = t205 * qJD(2);
t304 = t202 * t327;
t20 = -t59 * mrSges(7,1) + t58 * mrSges(7,2);
t285 = -pkin(1) - t177;
t284 = -t199 * t95 + t202 * t61;
t283 = -t323 / 0.2e1;
t282 = t323 / 0.2e1;
t281 = t136 * mrSges(5,1) + t135 * mrSges(5,2);
t40 = -t146 * t199 + t202 * t76;
t90 = -qJDD(2) * mrSges(4,1) + t136 * mrSges(4,2);
t276 = pkin(2) * t341 + qJ(3) * t345 + t336;
t265 = mrSges(3,1) * t200 + mrSges(3,2) * t203;
t261 = mrSges(6,1) * t199 + mrSges(6,2) * t202;
t259 = mrSges(7,1) * t199 + t377;
t255 = Ifges(6,1) * t199 + t367;
t253 = Ifges(7,1) * t199 + t365;
t250 = -t200 * Ifges(5,2) - t369;
t248 = Ifges(6,2) * t202 + t368;
t246 = Ifges(7,2) * t202 + t366;
t244 = Ifges(5,5) * t200 - Ifges(5,6) * t203;
t242 = Ifges(6,5) * t199 + Ifges(6,6) * t202;
t240 = Ifges(7,5) * t199 + Ifges(7,6) * t202;
t137 = -qJD(2) * pkin(2) + t325;
t142 = t172 + t192;
t235 = t137 * t203 - t142 * t200;
t18 = qJ(6) * t126 + t22;
t232 = pkin(1) * t265;
t108 = -t199 * t345 - t343;
t106 = t200 * t344 - t340;
t230 = t79 * (t200 * mrSges(5,2) + t359);
t229 = t105 * (mrSges(4,1) * t200 - t357);
t228 = t200 * (Ifges(3,1) * t203 - t372);
t93 = -qJ(4) * t333 + t231;
t221 = t199 * t327 + t200 * t332;
t220 = t199 * t333 - t304;
t218 = t199 * t378 + t202 * t379;
t217 = -t199 * t379 + t202 * t378;
t37 = qJDD(2) * pkin(4) - t39;
t186 = t204 * pkin(7);
t178 = t203 * qJ(4);
t139 = qJD(2) * mrSges(5,2) - t308;
t138 = -pkin(1) - t337;
t134 = t339 * t202;
t133 = t339 * t199;
t132 = t441 * qJD(1);
t130 = pkin(2) * t335 - t161;
t129 = t264 * qJD(1);
t113 = pkin(1) + t307;
t112 = t261 * t203;
t109 = t200 * t340 - t344;
t107 = -t200 * t343 - t347;
t100 = -Ifges(5,6) * qJD(2) + qJD(1) * t250;
t96 = -t178 + (-pkin(5) * t199 + pkin(7)) * t203;
t94 = pkin(2) * t333 - t338;
t92 = -mrSges(4,2) * t135 + qJDD(2) * mrSges(4,3);
t89 = -qJDD(2) * mrSges(5,1) - mrSges(5,3) * t135;
t88 = t205 * t335 + t161;
t83 = t306 + t301;
t78 = t200 * t306 + t338;
t63 = pkin(5) * t220 - t93;
t60 = -pkin(5) * t125 + qJD(6) + t86;
t38 = qJDD(2) * t205 + t209;
t35 = qJ(6) * t348 + t41;
t31 = pkin(5) * t200 + qJ(6) * t342 + t40;
t30 = -pkin(3) * t135 + t234;
t21 = -mrSges(6,1) * t59 + mrSges(6,2) * t58;
t19 = qJ(6) * t125 + t23;
t14 = -pkin(5) * t59 + qJDD(6) + t37;
t13 = pkin(5) * t155 + t18;
t8 = -qJD(5) * t41 + t284;
t7 = -t146 * t329 + t318;
t6 = qJ(6) * t304 + (-qJ(6) * t333 - qJD(5) * t146 + t326) * t199 + t318;
t5 = t202 * t326 + t233 * qJD(2) + (-t127 + (-qJ(6) * t203 - t76) * t199) * qJD(5) + t284;
t9 = [(t1 * t342 - t13 * t221 - t19 * t220 + t2 * t348) * mrSges(7,3) + (-t200 * t38 + t203 * t39) * mrSges(5,3) + t430 * mrSges(4,2) + t415 * t93 + (t199 * t450 + t202 * t451) * t327 / 0.2e1 + (-t244 / 0.2e1 + t434 / 0.2e1) * qJD(2) ^ 2 + t266 * t194 + t30 * t441 + (-t463 / 0.2e1 - t97 * mrSges(5,3) - t142 * mrSges(4,2) + t464 / 0.2e1) * t333 + (t344 * t412 - m(3) * t336 - m(4) * t276 - t324 * (pkin(3) * t341 - qJ(4) * t201 + t276) + t458 * t109 - t457 * t108 + t420 * t201 + (-m(7) * (t166 * t200 - t350) - m(6) * t268 - mrSges(2,1) + t462) * t204) * g(2) + (-t443 * pkin(7) - t100 / 0.2e1 - t19 * mrSges(7,2) - t23 * mrSges(6,2) + t13 * mrSges(7,1) + t22 * mrSges(6,1) - t83 * mrSges(5,3) + t419 / 0.2e1 + t137 * mrSges(4,2)) * t331 - t135 * t252 / 0.2e1 + t146 * t91 - pkin(1) * (mrSges(3,1) * t135 + mrSges(3,2) * t136) + t138 * (mrSges(4,1) * t135 - mrSges(4,3) * t136) + t95 * t139 + t78 * t132 + m(6) * (t22 * t8 + t23 * t7 + t3 * t41 + t4 * t40) - t94 * t129 - t37 * t112 + t96 * t20 + t6 * t72 + t7 * t73 + t5 * t74 + t8 * t75 + t63 * t66 + t342 * t469 + (t228 + t200 * (Ifges(4,1) * t203 + t364) + t203 * (-Ifges(3,2) * t200 + t371)) * t282 + t1 * mrSges(7,1) * t200 + t113 * t281 + m(7) * (t1 * t31 + t13 * t5 + t14 * t96 + t19 * t6 + t2 * t35 + t60 * t63) + t40 * t27 + t41 * t29 + t31 * t26 + t35 * t28 - t136 * t250 / 0.2e1 - t200 * t472 - t232 * t323 + (Ifges(3,4) * t136 - Ifges(3,2) * t135) * t461 + (t200 * t477 - t436 * t203) * t123 / 0.2e1 + (-mrSges(3,1) * t394 - mrSges(3,2) * t393 + t467 * t200 + t478 * t396 + (Ifges(3,6) - t476) * t461) * qJDD(2) + t60 * (mrSges(7,1) * t220 + mrSges(7,2) * t221) + t86 * (mrSges(6,1) * t220 + mrSges(6,2) * t221) + ((-Ifges(5,4) + Ifges(4,5)) * t136 + t483 * t135) * t396 + m(5) * (t113 * t30 + t146 * t38 + t78 * t79 + t83 * t95) + t424 * t283 + (t347 * t412 - t324 * (-qJ(4) * t204 + t186) + (-m(3) - m(4)) * t186 + t458 * t107 - t457 * t106 + t420 * t204 + (-m(6) * (-t198 * t200 - pkin(1)) - m(5) * t285 + t179 + m(3) * pkin(1) - m(7) * (-t153 * t200 - pkin(1)) - m(4) * (t285 - t185) + mrSges(2,1) + t473 * t203 - t474) * t201) * g(1) + ((Ifges(3,1) + Ifges(4,1)) * t136 + (-Ifges(3,4) + Ifges(4,5)) * t135 + t432) * t200 / 0.2e1 - t200 * t471 - t56 * t264 + (-t22 * t221 - t220 * t23 + t3 * t348 + t342 * t4) * mrSges(6,3) + Ifges(2,3) * qJDD(1) + t200 * t470 - t14 * t259 * t203 + m(4) * (t105 * t94 + t138 * t56 + (qJD(2) * t235 + t430) * pkin(7)) + (-t135 * t393 + t136 * t394 + t431) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t431) + ((t253 + t255) * t327 + t437 * qJD(2)) * t405 + ((t246 + t248) * t327 + t438 * qJD(2)) * t125 / 0.2e1 + ((t240 + t242) * t327 + t439 * qJD(2)) * t155 / 0.2e1 - t442 * t314 - t200 * (Ifges(5,4) * t135 - Ifges(5,2) * t136) / 0.2e1 + (t200 * t453 - t203 * t435) * t59 / 0.2e1 + (t200 * t455 - t203 * t433) * t58 / 0.2e1 + t92 * t393 + t90 * t394 + qJD(2) * t229 + qJD(2) * t230 + (-m(5) * t39 + m(6) * t37 + t21 - t89) * (-t178 + t393) + (t200 * Ifges(3,1) + t258 + t371) * t136 / 0.2e1 + (-t203 * Ifges(4,3) + t257 + t364) * t135 / 0.2e1 + t460 * t348 / 0.2e1; (-t242 / 0.2e1 - t240 / 0.2e1) * t123 + (Ifges(5,3) + Ifges(4,2) + Ifges(3,3)) * qJDD(2) + (-t39 * qJ(3) - t128 * t97 - t131 * t83 + t205 * t38 - t79 * t88) * m(5) + (-t253 / 0.2e1 - t255 / 0.2e1) * t58 + (m(4) * t142 + t145 - t415) * qJD(3) + (t128 * t86 + t198 * t37 - t22 * t32 - t23 * t33) * m(6) + t155 * (-t259 * t60 - t261 * t86) + t463 * t335 / 0.2e1 - (Ifges(4,1) * t334 + t169 + t464) * t335 / 0.2e1 + t467 * t136 + t14 * t260 + t37 * t262 + t153 * t20 - t131 * t139 + t130 * t129 - t88 * t132 + t133 * t26 - t134 * t28 - t118 * mrSges(3,2) - t119 * mrSges(3,1) - t87 * mrSges(4,1) - pkin(2) * t90 + t80 * mrSges(4,3) - t33 * t73 - t32 * t75 + (-t203 * t417 + t418) * t468 + (t438 * t485 + t439 * t484 - t230 - t229 - m(4) * t235 * pkin(7) - t13 * (mrSges(7,1) * t203 - mrSges(7,3) * t346) - t22 * (mrSges(6,1) * t203 - mrSges(6,3) * t346) - t19 * (-mrSges(7,2) * t203 - mrSges(7,3) * t349) - t23 * (-mrSges(6,2) * t203 - mrSges(6,3) * t349) + t437 * t404 + (-t228 / 0.2e1 + t232 + t424 / 0.2e1) * qJD(1)) * qJD(1) - t137 * t310 + t38 * mrSges(5,2) - t39 * mrSges(5,1) + t100 * t334 / 0.2e1 + (-t246 / 0.2e1 - t248 / 0.2e1) * t59 + (-t89 + t92) * qJ(3) + (t204 * t418 - t341 * t417) * g(1) + t476 * t135 + (-pkin(2) * t87 + qJ(3) * t80 - t105 * t130) * m(4) + (-m(4) * t337 - m(5) * t307 - m(7) * (t307 - t350) - m(6) * t237 + t475 * t200 + t462) * g(3) - (-Ifges(3,2) * t335 + t170 + t419) * t334 / 0.2e1 + t199 * t469 + t198 * t21 + t427 * t265 + (t13 * t328 + t19 * t329 - t428) * mrSges(7,3) + (t22 * t328 + t23 * t329 - t429) * mrSges(6,3) + t433 * t426 / 0.2e1 + t434 * t283 - (t125 * t435 + t155 * t436) * qJD(5) / 0.2e1 + t442 * t171 + t443 * t172 + t97 * t308 + t83 * t309 + t142 * t311 + t444 * t66 + t445 * t128 + t205 * t91 + t448 * t74 + t449 * t72 + (t1 * t133 + t448 * t13 - t134 * t2 + t14 * t153 + t449 * t19 + t444 * t60) * m(7) - t450 * t328 / 0.2e1 + t451 * t329 / 0.2e1 + t244 * t282 + (-t199 * t27 + t202 * t29 - t328 * t75 - t329 * t73 + t459) * t190 - t460 * t202 / 0.2e1; t380 * t202 - t381 * t199 + t217 * qJD(5) + (-t145 - t317) * qJD(2) + t312 * t390 + ((-t129 - t132 + t217) * qJD(1) - t427 * t312) * t200 + t459 - m(6) * (qJD(2) * t86 + t238 * t335) + t90 + t91 + (-qJD(2) * t60 + t428 - t155 * (t13 * t202 + t19 * t199)) * m(7) + (qJD(2) * t97 - t335 * t79 + t38) * m(5) + (-qJD(2) * t142 + t105 * t335 + t87) * m(4); t381 * t202 + t380 * t199 + t218 * qJD(5) + m(5) * t30 + m(7) * (t1 * t202 + t199 * t2 + (-t13 * t199 + t19 * t202) * qJD(5)) + m(6) * (t199 * t3 + t202 * t4 + (-t199 * t22 + t202 * t23) * qJD(5)) + (t317 * t203 + (t139 + t218) * t200 - m(5) * (-t200 * t83 + t203 * t97) - m(6) * (-t203 * t86 + t22 * t349 - t23 * t346) - m(7) * (t13 * t349 - t19 * t346 - t203 * t60)) * qJD(1) + t281 + (g(1) * t201 - g(2) * t204) * t324; (t106 * t425 - t107 * t457) * g(2) + (-t108 * t425 + t109 * t457) * g(1) + (t75 - t375) * t23 + (-t73 + t376) * t22 - t86 * (-mrSges(6,1) * t126 + mrSges(6,2) * t125) - t60 * (-mrSges(7,1) * t126 + mrSges(7,2) * t125) + (t74 - t373 - m(7) * (-t13 + t18)) * t19 - g(3) * t112 - t18 * t72 - t472 - t471 + t470 + t451 * t405 - (t199 * t447 + t377) * t390 + t447 * t1 + t13 * t374 + (t125 * t455 + t126 * t453) * t484 + (t126 * t479 + t450 + t466) * t485 + t432 + (t125 * t482 + t465) * t404 + ((m(7) * t60 + t66) * t126 + t26) * pkin(5); -t125 * t72 - t126 * t74 + (-g(3) * t200 - t125 * t19 - t126 * t13 - t203 * t427 + t14) * m(7) + t20;];
tau  = t9;
