% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:50
% EndTime: 2019-03-09 09:50:33
% DurationCPUTime: 27.93s
% Computational Cost: add. (7799->742), mult. (17565->894), div. (0->0), fcn. (12253->10), ass. (0->329)
t510 = Ifges(7,4) + Ifges(6,5);
t513 = -Ifges(5,4) + t510;
t494 = -m(7) - m(6);
t512 = m(5) - t494;
t491 = Ifges(6,4) + Ifges(5,5);
t452 = Ifges(7,5) - t491;
t239 = sin(qJ(2));
t342 = qJD(1) * qJD(2);
t316 = t239 * t342;
t242 = cos(qJ(2));
t341 = qJDD(1) * t242;
t201 = -t316 + t341;
t202 = qJDD(1) * t239 + t242 * t342;
t236 = sin(pkin(9));
t376 = cos(pkin(9));
t141 = t236 * t201 + t202 * t376;
t238 = sin(qJ(4));
t241 = cos(qJ(4));
t303 = qJD(1) * t376;
t352 = qJD(1) * t242;
t179 = -t236 * t352 - t239 * t303;
t269 = t241 * qJD(2) + t179 * t238;
t497 = qJD(4) * t269;
t81 = qJDD(2) * t238 + t141 * t241 + t497;
t437 = t81 / 0.2e1;
t150 = qJD(2) * t238 - t179 * t241;
t82 = qJD(4) * t150 - t241 * qJDD(2) + t141 * t238;
t435 = t82 / 0.2e1;
t140 = t376 * t201 - t202 * t236;
t136 = qJDD(4) - t140;
t432 = t136 / 0.2e1;
t511 = mrSges(6,2) + mrSges(5,3);
t489 = Ifges(7,2) + Ifges(6,3);
t488 = Ifges(5,6) - Ifges(6,6);
t487 = Ifges(7,6) - Ifges(6,6);
t486 = -Ifges(5,3) - Ifges(6,2);
t353 = qJD(1) * t239;
t178 = -t236 * t353 + t242 * t303;
t170 = qJD(4) - t178;
t453 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t290 = t241 * mrSges(6,1) + t238 * mrSges(6,3);
t292 = mrSges(5,1) * t241 - mrSges(5,2) * t238;
t380 = t238 * mrSges(7,2);
t508 = t290 + t292 + t380;
t507 = m(7) * qJ(6) + mrSges(7,3) - t511;
t377 = qJDD(2) / 0.2e1;
t448 = 0.2e1 * t377;
t506 = -t432 * t452 + t435 * t513 + t453 * t437;
t240 = sin(qJ(1));
t505 = g(2) * t240;
t169 = Ifges(4,4) * t178;
t504 = t178 * Ifges(4,2);
t207 = -mrSges(3,1) * t242 + mrSges(3,2) * t239;
t503 = -m(3) * pkin(1) - mrSges(2,1) + t207;
t237 = -qJ(3) - pkin(7);
t206 = t237 * t239;
t199 = qJD(1) * t206;
t208 = t237 * t242;
t200 = qJD(1) * t208;
t306 = t376 * t200;
t137 = t199 * t236 - t306;
t502 = qJD(5) * t238 + t137;
t501 = t510 * t150;
t359 = t238 * qJ(5);
t408 = t241 * pkin(4);
t276 = t359 + t408;
t500 = t510 * t269;
t499 = t510 * t241;
t498 = t510 * t238;
t182 = t236 * t200;
t189 = qJD(2) * pkin(2) + t199;
t128 = t376 * t189 + t182;
t234 = t242 * pkin(2);
t227 = t234 + pkin(1);
t203 = -qJD(1) * t227 + qJD(3);
t496 = t203 * mrSges(4,2) - t128 * mrSges(4,3);
t129 = t236 * t189 - t306;
t434 = pkin(4) + pkin(5);
t118 = qJD(2) * pkin(8) + t129;
t98 = -pkin(3) * t178 + pkin(8) * t179 + t203;
t46 = -t238 * t118 + t241 * t98;
t31 = qJ(6) * t150 + t46;
t466 = qJD(5) - t31;
t22 = -t170 * t434 + t466;
t158 = t170 * qJ(5);
t47 = t241 * t118 + t238 * t98;
t32 = -qJ(6) * t269 + t47;
t27 = t158 + t32;
t34 = -pkin(4) * t170 + qJD(5) - t46;
t35 = t158 + t47;
t495 = t203 * mrSges(4,1) + t46 * mrSges(5,1) - t34 * mrSges(6,1) - t22 * mrSges(7,1) - t47 * mrSges(5,2) + t27 * mrSges(7,2) - t129 * mrSges(4,3) + t35 * mrSges(6,3);
t493 = t201 / 0.2e1;
t492 = qJD(2) / 0.2e1;
t485 = -t136 * t487 + t489 * t82 + t510 * t81;
t37 = mrSges(5,1) * t136 - mrSges(5,3) * t81;
t38 = -t136 * mrSges(6,1) + t81 * mrSges(6,2);
t484 = -t37 + t38;
t40 = -mrSges(5,2) * t136 - mrSges(5,3) * t82;
t41 = -mrSges(6,2) * t82 + mrSges(6,3) * t136;
t483 = t41 + t40;
t482 = -t170 * t487 - t269 * t489 + t501;
t481 = t150 * t491 - t170 * t486 + t269 * t488;
t402 = mrSges(6,2) * t269;
t92 = mrSges(6,3) * t170 + t402;
t398 = mrSges(7,3) * t269;
t93 = mrSges(7,2) * t170 - t398;
t480 = t92 + t93;
t400 = mrSges(5,3) * t269;
t94 = -mrSges(5,2) * t170 + t400;
t479 = -t94 - t92;
t399 = mrSges(5,3) * t150;
t96 = mrSges(5,1) * t170 - t399;
t401 = mrSges(6,2) * t150;
t97 = -mrSges(6,1) * t170 + t401;
t478 = t97 - t96;
t383 = t179 * Ifges(4,4);
t468 = Ifges(4,6) * qJD(2);
t475 = t150 * Ifges(7,5) - Ifges(7,6) * t269 - t170 * Ifges(7,3) - t383 + t468 + t504;
t374 = qJ(5) * t241;
t265 = -t238 * t434 + t374;
t474 = t170 * t265 + t502;
t346 = qJD(4) * t238;
t416 = pkin(2) * t236;
t223 = pkin(8) + t416;
t354 = qJ(6) - t223;
t372 = t178 * t238;
t335 = pkin(2) * t353;
t108 = -pkin(3) * t179 - pkin(8) * t178 + t335;
t138 = t199 * t376 + t182;
t56 = t238 * t108 + t241 * t138;
t43 = -t179 * qJ(5) + t56;
t473 = -qJ(6) * t372 - qJD(6) * t241 + t346 * t354 - t43;
t125 = t238 * t138;
t193 = t354 * t241;
t472 = -qJD(4) * t193 - qJD(6) * t238 - t125 - (-qJ(6) * t178 - t108) * t241 - t434 * t179;
t275 = pkin(4) * t238 - t374;
t471 = t170 * t275 - t502;
t470 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t269 + mrSges(5,2) * t150 - mrSges(4,3) * t179;
t469 = Ifges(4,5) * qJD(2);
t467 = t241 * t434;
t197 = t236 * t242 + t239 * t376;
t180 = t197 * qJD(2);
t256 = -t236 * t239 + t242 * t376;
t465 = t180 * qJ(5) - qJD(5) * t256;
t464 = -t238 * t488 + t241 * t491;
t463 = t238 * t489 + t499;
t462 = -t136 * t486 - t488 * t82 + t491 * t81;
t345 = qJD(4) * t241;
t371 = t178 * t241;
t461 = t345 - t371;
t460 = t346 - t372;
t228 = pkin(7) * t341;
t194 = -pkin(7) * t316 + t228;
t195 = t202 * pkin(7);
t459 = t194 * t242 + t195 * t239;
t373 = qJDD(1) * pkin(1);
t166 = -pkin(2) * t201 + qJDD(3) - t373;
t51 = -pkin(3) * t140 - pkin(8) * t141 + t166;
t350 = qJD(3) * t239;
t123 = qJDD(2) * pkin(2) - qJ(3) * t202 - qJD(1) * t350 - t195;
t351 = qJD(2) * t239;
t333 = pkin(7) * t351;
t349 = qJD(3) * t242;
t134 = qJ(3) * t201 + t228 + (-t333 + t349) * qJD(1);
t62 = t236 * t123 + t376 * t134;
t58 = qJDD(2) * pkin(8) + t62;
t6 = -t118 * t346 + t238 * t51 + t241 * t58 + t98 * t345;
t7 = -t118 * t345 - t238 * t58 + t241 * t51 - t98 * t346;
t458 = -t238 * t7 + t241 * t6;
t3 = t136 * qJ(5) + t170 * qJD(5) + t6;
t257 = qJDD(5) - t7;
t5 = -pkin(4) * t136 + t257;
t457 = t238 * t5 + t241 * t3;
t243 = cos(qJ(1));
t456 = g(1) * t243 + t505;
t454 = -mrSges(7,2) - mrSges(6,3) + mrSges(5,2);
t148 = Ifges(5,4) * t269;
t450 = t150 * t453 - t170 * t452 + t148 - t500;
t449 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t447 = m(7) * pkin(5) + mrSges(7,1);
t235 = qJ(2) + pkin(9);
t230 = sin(t235);
t231 = cos(t235);
t299 = -m(7) * t434 - mrSges(7,1);
t415 = pkin(2) * t239;
t446 = t512 * t415 + t507 * t231 + (-m(7) * (-pkin(3) - t359) - t299 * t241 - m(6) * (-pkin(3) - t276) + m(5) * pkin(3) + t508) * t230;
t379 = t241 * mrSges(7,2);
t288 = -t238 * mrSges(7,1) + t379;
t378 = t241 * mrSges(6,3);
t289 = mrSges(6,1) * t238 - t378;
t291 = mrSges(5,1) * t238 + mrSges(5,2) * t241;
t296 = qJD(2) * pkin(3) + t128;
t252 = qJ(5) * t150 + t296;
t30 = t269 * t434 + qJD(6) + t252;
t48 = -pkin(4) * t269 - t252;
t445 = t30 * t288 + t48 * t289 - t291 * t296;
t394 = Ifges(5,4) * t238;
t444 = t241 * t453 - t394 + t498;
t293 = t231 * mrSges(4,1) - t230 * mrSges(4,2);
t443 = t230 * t511 + t293;
t442 = pkin(8) * t512;
t441 = mrSges(5,1) + mrSges(6,1) + t447;
t1 = -qJ(6) * t81 - qJD(6) * t150 - t136 * t434 + t257;
t2 = qJ(6) * t82 - qJD(6) * t269 + t3;
t440 = t7 * mrSges(5,1) - t5 * mrSges(6,1) - t1 * mrSges(7,1) - t6 * mrSges(5,2) + t2 * mrSges(7,2) + t3 * mrSges(6,3);
t436 = -t82 / 0.2e1;
t433 = -t136 / 0.2e1;
t431 = t269 / 0.2e1;
t430 = -t269 / 0.2e1;
t429 = -t150 / 0.2e1;
t428 = t150 / 0.2e1;
t427 = -t170 / 0.2e1;
t426 = t170 / 0.2e1;
t425 = -t178 / 0.2e1;
t424 = -t179 / 0.2e1;
t423 = t179 / 0.2e1;
t414 = pkin(7) * t242;
t411 = g(3) * t230;
t220 = t230 * pkin(8);
t221 = t231 * pkin(3);
t85 = -mrSges(6,1) * t269 - mrSges(6,3) * t150;
t86 = mrSges(7,1) * t269 + mrSges(7,2) * t150;
t404 = t85 - t86;
t397 = mrSges(7,3) * t150;
t396 = Ifges(3,4) * t239;
t395 = Ifges(3,4) * t242;
t393 = Ifges(5,4) * t241;
t384 = t150 * Ifges(5,4);
t375 = qJ(5) * t269;
t181 = t256 * qJD(2);
t370 = t181 * t238;
t369 = t181 * t241;
t366 = t197 * t238;
t362 = t230 * t243;
t361 = t231 * t243;
t360 = t237 * t243;
t358 = t238 * t240;
t357 = t240 * t241;
t356 = t241 * t243;
t355 = t243 * t238;
t61 = t376 * t123 - t236 * t134;
t127 = -pkin(3) * t256 - pkin(8) * t197 - t227;
t145 = t236 * t206 - t208 * t376;
t74 = t238 * t127 + t241 * t145;
t343 = qJD(5) * t241;
t338 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t353) * t414;
t334 = pkin(2) * t351;
t307 = qJD(2) * t237;
t176 = t239 * t307 + t349;
t177 = t242 * t307 - t350;
t105 = t176 * t376 + t236 * t177;
t332 = t238 * t105 + t127 * t346 + t145 * t345;
t109 = pkin(3) * t180 - pkin(8) * t181 + t334;
t329 = t241 * t105 + t238 * t109 + t127 * t345;
t50 = -qJ(5) * t256 + t74;
t328 = t221 + t220 + t234;
t327 = t376 * pkin(2);
t326 = t197 * t346;
t325 = t197 * t345;
t25 = -t82 * mrSges(7,1) + t81 * mrSges(7,2);
t311 = -t346 / 0.2e1;
t309 = t345 / 0.2e1;
t36 = -t136 * mrSges(7,1) - t81 * mrSges(7,3);
t308 = -t227 - t221;
t304 = -t140 * mrSges(4,1) + t141 * mrSges(4,2);
t55 = t108 * t241 - t125;
t142 = t238 * t145;
t73 = t127 * t241 - t142;
t104 = t176 * t236 - t376 * t177;
t144 = -t376 * t206 - t208 * t236;
t300 = t243 * t227 - t237 * t240;
t224 = -t327 - pkin(3);
t295 = qJDD(2) * pkin(3) + t61;
t294 = mrSges(3,1) * t239 + mrSges(3,2) * t242;
t284 = t242 * Ifges(3,2) + t396;
t283 = -Ifges(5,2) * t238 + t393;
t280 = Ifges(3,5) * t242 - Ifges(3,6) * t239;
t277 = Ifges(7,5) * t241 + Ifges(7,6) * t238;
t271 = t231 * t276 + t328;
t270 = -qJ(6) * t181 - qJD(6) * t197;
t21 = t109 * t241 - t332;
t267 = pkin(3) * t361 + pkin(8) * t362 + t300;
t266 = pkin(1) * t294;
t258 = t239 * (Ifges(3,1) * t242 - t396);
t20 = -t145 * t346 + t329;
t254 = t224 - t359;
t253 = Ifges(7,5) * t81 + Ifges(7,6) * t82 - Ifges(7,3) * t136;
t171 = t231 * t358 + t356;
t173 = t231 * t355 - t357;
t250 = -g(1) * t173 - g(2) * t171 - t238 * t411;
t248 = qJ(5) * t81 + qJD(5) * t150 + t295;
t229 = Ifges(3,4) * t352;
t209 = t230 * t374;
t205 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t352;
t192 = t354 * t238;
t191 = t254 - t408;
t186 = Ifges(3,1) * t353 + Ifges(3,5) * qJD(2) + t229;
t185 = Ifges(3,6) * qJD(2) + qJD(1) * t284;
t174 = t231 * t356 + t358;
t172 = t231 * t357 - t355;
t159 = -t254 + t467;
t154 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t178;
t124 = -mrSges(4,1) * t178 - mrSges(4,2) * t179;
t120 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t141;
t119 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t140;
t115 = -t179 * Ifges(4,1) + t169 + t469;
t95 = -mrSges(7,1) * t170 - t397;
t84 = pkin(4) * t150 - t375;
t83 = t197 * t275 + t144;
t69 = Ifges(5,2) * t269 + t170 * Ifges(5,6) + t384;
t60 = t197 * t265 - t144;
t59 = -t150 * t434 + t375;
t52 = pkin(4) * t256 - t73;
t44 = pkin(4) * t179 - t55;
t42 = qJ(6) * t366 + t50;
t39 = mrSges(7,2) * t136 + mrSges(7,3) * t82;
t33 = t142 + (-qJ(6) * t197 - t127) * t241 + t434 * t256;
t26 = mrSges(5,1) * t82 + mrSges(5,2) * t81;
t24 = mrSges(6,1) * t82 - mrSges(6,3) * t81;
t23 = t275 * t181 + (qJD(4) * t276 - t343) * t197 + t104;
t16 = t81 * Ifges(5,4) - t82 * Ifges(5,2) + t136 * Ifges(5,6);
t13 = t265 * t181 + (t343 + (-t359 - t467) * qJD(4)) * t197 - t104;
t12 = -pkin(4) * t180 - t21;
t11 = t20 + t465;
t10 = pkin(4) * t82 - t248;
t9 = qJ(6) * t325 + (-qJD(4) * t145 - t270) * t238 + t329 + t465;
t8 = qJ(6) * t326 - t434 * t180 + (-t109 + t270) * t241 + t332;
t4 = -t434 * t82 + qJDD(6) + t248;
t14 = [t450 * t369 / 0.2e1 + (Ifges(4,5) * t492 + Ifges(4,1) * t424 + t169 / 0.2e1 + t115 / 0.2e1 + t496) * t181 - qJDD(2) * mrSges(3,2) * t414 + (t495 - t475 / 0.2e1 - t452 * t428 - t487 * t430 - t486 * t426 - Ifges(4,6) * t492 - Ifges(4,4) * t424 - Ifges(7,3) * t427 + Ifges(5,6) * t431 + t481 / 0.2e1 - t504 / 0.2e1) * t180 + t242 * t186 * t492 + (m(5) * t360 + t494 * (-t172 * pkin(4) - qJ(5) * t171 - t360) + (m(4) * t237 + t449) * t243 + t441 * t172 - t454 * t171 + (-m(7) * t308 - (m(7) * (-pkin(8) + qJ(6)) + mrSges(7,3)) * t230 + m(4) * t227 + (-m(5) - m(6)) * (t308 - t220) + t443 - t503) * t240) * g(1) + (-t16 / 0.2e1 - t6 * mrSges(5,3) - t3 * mrSges(6,2) + t2 * mrSges(7,3) + t485 / 0.2e1) * t366 + t284 * t493 + t202 * t395 / 0.2e1 + (-m(4) * t61 - m(5) * t295 - t120 + t26) * t144 + (-m(4) * t128 - m(5) * t296 + t470) * t104 - t185 * t351 / 0.2e1 + (-t296 * mrSges(5,1) + t48 * mrSges(6,1) + t27 * mrSges(7,3) - t30 * mrSges(7,1) - t47 * mrSges(5,3) - t35 * mrSges(6,2) + Ifges(7,6) * t427 - Ifges(5,2) * t431 - t69 / 0.2e1 + t513 * t428 + t489 * t430 - t488 * t426) * (t325 + t370) - t207 * t373 + ((mrSges(6,2) * t5 - mrSges(5,3) * t7 - mrSges(7,3) * t1 + t506) * t241 + t166 * mrSges(4,2) - t61 * mrSges(4,3) + Ifges(4,1) * t141 + Ifges(4,4) * t140 + Ifges(4,5) * t448 + t10 * t289 + t277 * t433 + t283 * t436 + t4 * t288 - t291 * t295 + t309 * t482 + t311 * t450 + t432 * t464 + t435 * t463 + t437 * t444) * t197 - t227 * t304 - t266 * t342 + (mrSges(5,2) * t296 - mrSges(6,2) * t34 - mrSges(7,2) * t30 + mrSges(5,3) * t46 + mrSges(6,3) * t48 + mrSges(7,3) * t22 - Ifges(5,4) * t431 - Ifges(7,5) * t427 - t426 * t491 - t428 * t453 - t430 * t510) * (t326 - t369) + m(5) * (t20 * t47 + t21 * t46 + t6 * t74 + t7 * t73) + m(4) * (t105 * t129 + t145 * t62 - t166 * t227 + t203 * t334) + (t258 + t242 * (-Ifges(3,2) * t239 + t395)) * t342 / 0.2e1 + m(7) * (t1 * t33 + t13 * t30 + t2 * t42 + t22 * t8 + t27 * t9 + t4 * t60) + m(6) * (t10 * t83 + t11 * t35 + t12 * t34 + t23 * t48 + t3 * t50 + t5 * t52) - t205 * t333 + t124 * t334 + (-m(4) * t300 - m(5) * t267 + t494 * (t174 * pkin(4) + qJ(5) * t173 + t267) + t507 * t362 + (-t293 + t503) * t243 + t449 * t240 - t441 * t174 + t454 * t173) * g(2) + t242 * (Ifges(3,4) * t202 + Ifges(3,2) * t201 + Ifges(3,6) * qJDD(2)) / 0.2e1 + Ifges(2,3) * qJDD(1) - pkin(1) * (-mrSges(3,1) * t201 + mrSges(3,2) * t202) + t105 * t154 + t145 * t119 + t11 * t92 + t9 * t93 + t20 * t94 + t8 * t95 + t21 * t96 + t12 * t97 + t23 * t85 + t13 * t86 + t83 * t24 + t73 * t37 + t74 * t40 + t60 * t25 + t50 * t41 + t52 * t38 + t42 * t39 + t33 * t36 + (t201 * t414 + t459) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t459) + (-t462 / 0.2e1 + t253 / 0.2e1 + Ifges(4,2) * t140 + Ifges(4,4) * t141 + Ifges(7,3) * t433 - Ifges(5,6) * t436 - t166 * mrSges(4,1) + t62 * mrSges(4,3) + t487 * t435 + t486 * t432 + t452 * t437 - t440 + Ifges(4,6) * t448) * t256 + t482 * t370 / 0.2e1 + (t280 * t492 - t338) * qJD(2) + (Ifges(3,1) * t202 + Ifges(3,4) * t493 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t202) + t448 * Ifges(3,5)) * t239 + Ifges(3,6) * t242 * t377; (t346 / 0.2e1 - t372 / 0.2e1) * t482 + (Ifges(4,2) * t425 + Ifges(7,3) * t426 - Ifges(5,6) * t430 - t468 / 0.2e1 + t487 * t431 + t486 * t427 + t452 * t429 + t495) * t179 + (t10 * t191 - t34 * t44 - t35 * t43 + t471 * t48) * m(6) + (t478 * t345 + t479 * t346 + t483 * t241 + t484 * t238 + ((-t238 * t47 - t241 * t46) * qJD(4) + t458) * m(5) + ((-t238 * t35 + t241 * t34) * qJD(4) + t457) * m(6)) * t223 + (t309 - t371 / 0.2e1) * t450 + ((t236 * t62 + t376 * t61) * pkin(2) + t128 * t137 - t129 * t138 - t203 * t335) * m(4) + (t243 * t446 - t361 * t442) * g(1) + (-t231 * t442 + t446) * t505 + (t338 + (t266 - t258 / 0.2e1) * qJD(1)) * qJD(1) + (t428 * t444 + t445) * qJD(4) + (-t241 * t489 + t498) * t435 + (t238 * t453 + t393 - t499) * t437 + t295 * t292 + (t137 * t296 - t224 * t295 - t46 * t55 - t47 * t56) * m(5) + t4 * (t241 * mrSges(7,1) + t380) + (t372 / 0.2e1 + t311) * t69 - t280 * t342 / 0.2e1 + (t115 + t169) * t425 - (-Ifges(3,2) * t353 + t186 + t229) * t352 / 0.2e1 + (m(4) * t415 + mrSges(4,1) * t230 + mrSges(4,2) * t231 + t294) * t456 - t10 * t290 + (Ifges(4,1) * t423 + t277 * t426 + t283 * t430 - t469 / 0.2e1 + t463 * t431 + t464 * t427 + t444 * t429 - t445 - t496) * t178 - t124 * t335 + (-t277 / 0.2e1 + t464 / 0.2e1) * qJD(4) * t170 + (Ifges(7,5) * t238 - Ifges(7,6) * t241) * t433 + (Ifges(5,2) * t241 + t394) * t436 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t241 * t16 / 0.2e1 + t224 * t26 + Ifges(3,6) * t201 + Ifges(3,5) * t202 + t191 * t24 - t192 * t36 - t193 * t39 - t194 * mrSges(3,2) - t195 * mrSges(3,1) + t159 * t25 - t138 * t154 + t119 * t416 + Ifges(4,6) * t140 + Ifges(4,5) * t141 - t43 * t92 - t56 * t94 - t55 * t96 - t44 * t97 + t61 * mrSges(4,1) - t62 * mrSges(4,2) + (-m(7) * (-qJ(6) * t230 + t271) + t230 * mrSges(7,3) - m(4) * t234 - m(5) * t328 - m(6) * t271 + t207 + (-t241 * t447 - t508) * t231 - t443) * g(3) + (t34 * t461 - t35 * t460 + t457) * mrSges(6,2) + (-t1 * t238 - t2 * t241 - t22 * t461 + t27 * t460) * mrSges(7,3) + (-t46 * t461 - t460 * t47 + t458) * mrSges(5,3) - t470 * t137 + t471 * t85 + t472 * t95 + t473 * t93 + t474 * t86 + (-t1 * t192 + t159 * t4 - t193 * t2 + t22 * t472 + t27 * t473 + t30 * t474) * m(7) + t475 * t424 + (t383 + t481) * t423 - t485 * t241 / 0.2e1 + (t238 * t491 + t241 * t488) * t432 + t120 * t327 + (t283 / 0.2e1 - t463 / 0.2e1) * t497 + (t185 / 0.2e1 + pkin(7) * t205) * t353 + t238 * t506; -t178 * t154 + (t404 + t470) * t179 + (-t36 + t170 * (t94 + t480) - t484) * t241 + (t39 + t170 * (t95 + t478) + t483) * t238 + t304 + (-g(1) * t240 + g(2) * t243) * (m(4) + t512) + (-t1 * t241 - t179 * t30 + t2 * t238 + t170 * (t22 * t238 + t241 * t27)) * m(7) + (t179 * t48 + t238 * t3 - t241 * t5 + t170 * (t238 * t34 + t241 * t35)) * m(6) + (-t296 * t179 + t238 * t6 + t241 * t7 + t170 * (-t238 * t46 + t241 * t47)) * m(5) + (-t128 * t179 - t129 * t178 + t166) * m(4); (-pkin(4) * t5 + qJ(5) * t3 + qJD(5) * t35 - t48 * t84) * m(6) + (-m(6) * t209 + (-t378 - t379 + (m(6) * pkin(4) + mrSges(6,1) - t299) * t238) * t230) * g(3) + (t150 * t489 + t500) * t431 + (t269 * t453 - t384 + t482 + t501) * t429 + t296 * (mrSges(5,1) * t150 + mrSges(5,2) * t269) + (Ifges(7,5) * t269 + Ifges(7,6) * t150) * t426 - t48 * (mrSges(6,1) * t150 - mrSges(6,3) * t269) - t30 * (-mrSges(7,1) * t150 + mrSges(7,2) * t269) + (-t150 * t488 + t269 * t491) * t427 - t434 * t36 + (-g(3) * t209 + t2 * qJ(5) - t1 * t434 - t22 * t32 + t27 * t466 - t30 * t59) * m(7) + (-Ifges(5,2) * t150 + t148 + t450) * t430 + t462 + t69 * t428 + (t39 + t41) * qJ(5) - t27 * t397 + t22 * t398 - t253 + t35 * t401 - t34 * t402 + t291 * t411 - t31 * t93 - t32 * t95 - t84 * t85 - t59 * t86 - pkin(4) * t38 + t440 + (-m(6) * t34 + t399 - t478) * t47 + (-m(6) * t35 + t400 + t479) * t46 + t480 * qJD(5) + (t494 * (-t173 * pkin(4) + qJ(5) * t174) + t454 * t174 + t441 * t173) * g(1) + (t494 * (-t171 * pkin(4) + qJ(5) * t172) + t454 * t172 + t441 * t171) * g(2); t404 * t150 - t480 * t170 + t36 + t38 + (-t150 * t30 - t170 * t27 + t1 + t250) * m(7) + (t150 * t48 - t170 * t35 + t250 + t5) * m(6); t269 * t93 + t150 * t95 + (-g(3) * t231 + t22 * t150 + t230 * t456 + t269 * t27 + t4) * m(7) + t25;];
tau  = t14;
