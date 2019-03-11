% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:30
% EndTime: 2019-03-09 04:32:01
% DurationCPUTime: 21.66s
% Computational Cost: add. (4768->670), mult. (9750->799), div. (0->0), fcn. (5738->10), ass. (0->292)
t484 = Ifges(7,4) + Ifges(6,5);
t485 = Ifges(6,4) + Ifges(5,5);
t482 = Ifges(7,2) + Ifges(6,3);
t481 = Ifges(5,6) - Ifges(6,6);
t493 = Ifges(7,5) - t485;
t479 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t205 = cos(qJ(3));
t492 = t205 / 0.2e1;
t202 = sin(qJ(3));
t376 = Ifges(4,4) * t202;
t491 = t376 / 0.2e1;
t204 = cos(qJ(4));
t490 = t484 * t204;
t201 = sin(qJ(4));
t489 = t484 * t201;
t325 = qJD(1) * qJD(3);
t152 = qJDD(1) * t202 + t205 * t325;
t326 = t204 * qJD(3);
t341 = qJD(1) * t202;
t145 = t201 * t341 - t326;
t335 = qJD(4) * t145;
t75 = qJDD(3) * t201 + t152 * t204 - t335;
t414 = t75 / 0.2e1;
t338 = qJD(3) * t201;
t146 = t204 * t341 + t338;
t76 = qJD(4) * t146 - t204 * qJDD(3) + t152 * t201;
t412 = t76 / 0.2e1;
t151 = t205 * qJDD(1) - t202 * t325;
t141 = qJDD(4) - t151;
t409 = t141 / 0.2e1;
t407 = t145 / 0.2e1;
t406 = -t146 / 0.2e1;
t405 = t146 / 0.2e1;
t340 = qJD(1) * t205;
t174 = -qJD(4) + t340;
t403 = t174 / 0.2e1;
t488 = Ifges(4,2) * t492 + t491;
t487 = qJD(3) / 0.2e1;
t486 = -mrSges(6,2) - mrSges(5,3);
t483 = Ifges(6,2) + Ifges(5,3);
t480 = Ifges(6,6) - Ifges(7,6);
t441 = t201 * t482 + t490;
t440 = -t481 * t201 + t204 * t485;
t477 = m(7) * qJ(6) + mrSges(7,3) + t486;
t373 = Ifges(5,4) * t201;
t422 = t204 * t479 - t373 + t489;
t262 = t204 * mrSges(6,1) + t201 * mrSges(6,3);
t264 = mrSges(5,1) * t204 - mrSges(5,2) * t201;
t365 = t201 * mrSges(7,2);
t476 = -t262 - t264 - t365;
t475 = t325 / 0.2e1;
t474 = (-Ifges(5,4) + t484) * t412 + t479 * t414 - t493 * t409;
t198 = qJ(1) + pkin(9);
t187 = cos(t198);
t473 = g(1) * t187;
t186 = sin(t198);
t472 = g(2) * t186;
t471 = t484 * t146;
t199 = sin(pkin(9));
t179 = pkin(1) * t199 + pkin(7);
t470 = qJD(2) * qJD(3) + t179 * qJDD(1);
t156 = t179 * qJD(1);
t339 = qJD(2) * t202;
t118 = t156 * t205 + t339;
t105 = qJD(3) * pkin(8) + t118;
t194 = t202 * pkin(8);
t196 = t205 * pkin(3);
t307 = -pkin(2) - t196;
t237 = t307 - t194;
t200 = cos(pkin(9));
t398 = pkin(1) * t200;
t110 = (t237 - t398) * qJD(1);
t331 = qJD(4) * t204;
t333 = qJD(4) * t201;
t337 = qJD(3) * t202;
t55 = t202 * qJDD(2) - t156 * t337 + t205 * t470;
t51 = qJDD(3) * pkin(8) + t55;
t180 = -pkin(2) - t398;
t155 = t180 * qJDD(1);
t74 = -pkin(3) * t151 - pkin(8) * t152 + t155;
t6 = -t105 * t333 + t110 * t331 + t201 * t74 + t204 * t51;
t7 = -t105 * t331 - t110 * t333 - t201 * t51 + t204 * t74;
t268 = -t201 * t7 + t204 * t6;
t42 = -t201 * t105 + t204 * t110;
t43 = t204 * t105 + t201 * t110;
t469 = -t42 * t331 - t43 * t333 + t268;
t4 = t141 * qJ(5) - t174 * qJD(5) + t6;
t233 = qJDD(5) - t7;
t5 = -pkin(4) * t141 + t233;
t269 = t201 * t5 + t204 * t4;
t30 = pkin(4) * t174 + qJD(5) - t42;
t168 = t174 * qJ(5);
t31 = -t168 + t43;
t468 = t30 * t331 - t31 * t333 + t269;
t467 = -qJD(5) * t201 - t339;
t466 = t484 * t145;
t185 = Ifges(4,4) * t340;
t138 = Ifges(5,4) * t145;
t427 = t479 * t146 + t174 * t493 - t138 + t466;
t455 = t482 * t145 - t480 * t174 + t471;
t465 = Ifges(4,1) * t341 + Ifges(4,5) * qJD(3) + t201 * t455 + t204 * t427 + t185;
t464 = Ifges(7,5) * t405 + Ifges(4,6) * t487 + qJD(1) * t488 + t485 * t406 + (Ifges(7,6) + t481) * t407 + (Ifges(7,3) + t483) * t403;
t411 = pkin(4) + pkin(5);
t1 = -qJ(6) * t75 - qJD(6) * t146 - t141 * t411 + t233;
t2 = qJ(6) * t76 + qJD(6) * t145 + t4;
t463 = -t7 * mrSges(5,1) + t5 * mrSges(6,1) + t1 * mrSges(7,1) + t6 * mrSges(5,2) - t2 * mrSges(7,2) - t4 * mrSges(6,3);
t462 = -m(6) - m(7);
t461 = mrSges(3,2) - mrSges(4,3);
t459 = t480 * t141 + t482 * t76 + t484 * t75;
t34 = mrSges(5,1) * t141 - mrSges(5,3) * t75;
t35 = -t141 * mrSges(6,1) + t75 * mrSges(6,2);
t458 = t35 - t34;
t37 = -mrSges(5,2) * t141 - mrSges(5,3) * t76;
t38 = -mrSges(6,2) * t76 + mrSges(6,3) * t141;
t457 = t38 + t37;
t379 = mrSges(5,3) * t146;
t98 = -mrSges(5,1) * t174 - t379;
t381 = mrSges(6,2) * t146;
t99 = mrSges(6,1) * t174 + t381;
t454 = -t98 + t99;
t266 = mrSges(4,1) * t205 - mrSges(4,2) * t202;
t453 = -mrSges(3,1) - t266;
t382 = mrSges(6,2) * t145;
t100 = -mrSges(6,3) * t174 - t382;
t378 = mrSges(7,3) * t145;
t95 = -mrSges(7,2) * t174 + t378;
t359 = t100 + t95;
t20 = mrSges(5,1) * t76 + mrSges(5,2) * t75;
t451 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t152 + t20;
t384 = pkin(8) - qJ(6);
t162 = t384 * t204;
t347 = t204 * t205;
t117 = t205 * qJD(2) - t202 * t156;
t393 = pkin(8) * t205;
t270 = pkin(3) * t202 - t393;
t149 = t270 * qJD(1);
t66 = -t201 * t117 + t149 * t204;
t450 = -(-qJ(6) * t347 - t202 * t411) * qJD(1) + t66 + qJD(4) * t162 - qJD(6) * t201;
t306 = t201 * t340;
t327 = qJD(6) * t204;
t67 = t204 * t117 + t201 * t149;
t48 = qJ(5) * t341 + t67;
t449 = -qJ(6) * t306 - t333 * t384 - t327 - t48;
t317 = t411 * t201;
t355 = qJ(5) * t204;
t236 = -t317 + t355;
t448 = -(qJD(1) * t236 - t156) * t205 + qJD(4) * t236 - t467;
t396 = pkin(4) * t201;
t239 = -t355 + t396;
t447 = -(qJD(1) * t239 + t156) * t205 + qJD(4) * t239 + t467;
t316 = mrSges(4,3) * t341;
t446 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t145 + mrSges(5,2) * t146 + t316;
t445 = t204 * t411;
t444 = qJ(5) * t337 - qJD(5) * t205;
t443 = t483 * t202 + t440 * t205;
t442 = t480 * t202 + t441 * t205;
t439 = -t482 * t204 + t489;
t438 = t485 * t201 + t481 * t204;
t26 = qJ(6) * t146 + t42;
t437 = -t26 + qJD(5);
t434 = t483 * t141 - t481 * t76 + t485 * t75;
t336 = qJD(3) * t205;
t56 = t205 * qJDD(2) - t156 * t336 - t202 * t470;
t433 = -t202 * t56 + t205 * t55;
t432 = t486 * t202;
t431 = t472 + t473;
t430 = -m(5) + t462;
t429 = mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t426 = m(7) * pkin(5) + mrSges(7,1);
t356 = qJ(5) * t201;
t240 = pkin(4) * t204 + t356;
t153 = -pkin(3) - t240;
t278 = -m(7) * t411 - mrSges(7,1);
t284 = pkin(3) + t356;
t425 = t477 * t205 + (m(5) * pkin(3) - m(6) * t153 + m(7) * t284 - t204 * t278 - t476) * t202;
t424 = -t202 * t493 + t422 * t205;
t271 = qJD(3) * pkin(3) + t117;
t226 = qJ(5) * t146 + t271;
t25 = -t145 * t411 + qJD(6) + t226;
t361 = t204 * mrSges(6,3);
t261 = t201 * mrSges(6,1) - t361;
t263 = mrSges(5,1) * t201 + mrSges(5,2) * t204;
t41 = pkin(4) * t145 - t226;
t423 = t271 * t263 - t41 * t261 - t25 * (-t201 * mrSges(7,1) + mrSges(7,2) * t204);
t372 = Ifges(5,4) * t204;
t421 = t479 * t201 + t372 - t490;
t17 = t174 * t411 + t437;
t27 = qJ(6) * t145 + t43;
t21 = -t168 + t27;
t418 = -t1 * t201 - t17 * t331 - t2 * t204 + t21 * t333;
t416 = mrSges(5,1) + mrSges(6,1) + t426;
t374 = Ifges(5,4) * t146;
t62 = -Ifges(5,2) * t145 - Ifges(5,6) * t174 + t374;
t415 = -t62 / 0.2e1;
t413 = -t76 / 0.2e1;
t410 = -t141 / 0.2e1;
t408 = -t145 / 0.2e1;
t404 = -t174 / 0.2e1;
t203 = sin(qJ(1));
t397 = pkin(1) * t203;
t206 = cos(qJ(1));
t197 = t206 * pkin(1);
t84 = mrSges(6,1) * t145 - mrSges(6,3) * t146;
t85 = -mrSges(7,1) * t145 + mrSges(7,2) * t146;
t383 = t84 - t85;
t380 = mrSges(5,3) * t145;
t377 = mrSges(7,3) * t146;
t375 = Ifges(4,4) * t205;
t357 = qJ(5) * t145;
t354 = qJ(6) * t202;
t352 = t187 * t202;
t351 = t187 * t205;
t350 = t201 * t202;
t349 = t201 * t205;
t348 = t202 * t204;
t344 = t196 + t194;
t134 = t180 - t344;
t148 = t179 * t347;
t346 = qJD(4) * t148 + t134 * t333;
t150 = t270 * qJD(3);
t345 = t134 * t331 + t201 * t150;
t82 = t201 * t134 + t148;
t332 = qJD(4) * t202;
t329 = qJD(5) * t204;
t97 = mrSges(7,1) * t174 - t377;
t323 = t97 + t454;
t96 = mrSges(5,2) * t174 - t380;
t322 = t96 + t359;
t319 = pkin(8) * t333;
t315 = mrSges(4,3) * t340;
t308 = t187 * pkin(2) + t186 * pkin(7) + t197;
t305 = t179 * t337;
t303 = t201 * t336;
t19 = -t76 * mrSges(7,1) + t75 * mrSges(7,2);
t286 = t331 / 0.2e1;
t33 = -t141 * mrSges(7,1) - t75 * mrSges(7,3);
t285 = t187 * pkin(7) - t397;
t283 = -t179 * t201 - pkin(4);
t147 = t179 * t349;
t81 = t134 * t204 - t147;
t279 = pkin(4) * t347 + qJ(5) * t349 + t344;
t78 = -qJ(5) * t205 + t82;
t267 = -t150 * t204 + t346;
t265 = mrSges(4,1) * t202 + mrSges(4,2) * t205;
t253 = -Ifges(5,2) * t201 + t372;
t252 = Ifges(5,2) * t204 + t373;
t247 = Ifges(4,5) * t205 - Ifges(4,6) * t202;
t242 = Ifges(7,5) * t204 + Ifges(7,6) * t201;
t241 = Ifges(7,5) * t201 - Ifges(7,6) * t204;
t238 = pkin(3) * t351 + pkin(8) * t352 + t308;
t231 = t180 * qJD(1) * t265;
t230 = t202 * (Ifges(4,1) * t205 - t376);
t229 = qJDD(3) * pkin(3) + t56;
t227 = Ifges(7,5) * t75 + Ifges(7,6) * t76 - Ifges(7,3) * t141;
t111 = t186 * t349 + t187 * t204;
t113 = -t186 * t204 + t187 * t349;
t223 = -g(1) * t113 - g(2) * t111 - g(3) * t350;
t217 = Ifges(5,6) * t202 + t205 * t253;
t213 = -Ifges(7,3) * t202 + t205 * t242;
t28 = (-t202 * t326 - t205 * t333) * t179 + t345;
t211 = qJ(5) * t75 + qJD(5) * t146 + t229;
t195 = t205 * pkin(4);
t173 = mrSges(7,2) * t348;
t171 = qJ(5) * t348;
t160 = t384 * t201;
t159 = -qJD(3) * mrSges(4,2) + t315;
t135 = t284 + t445;
t133 = t263 * t202;
t120 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t151;
t114 = t186 * t201 + t187 * t347;
t112 = t186 * t347 - t187 * t201;
t92 = -t171 + (t179 + t396) * t202;
t83 = pkin(4) * t146 + t357;
t80 = t171 + (-t179 - t317) * t202;
t79 = t195 - t81;
t53 = qJ(6) * t350 + t78;
t50 = -pkin(4) * t341 - t66;
t47 = -t146 * t411 - t357;
t44 = pkin(5) * t205 + t147 + t195 + (-t134 - t354) * t204;
t39 = (qJD(4) * t240 - t329) * t202 + (t179 + t239) * t336;
t36 = mrSges(7,2) * t141 + mrSges(7,3) * t76;
t29 = t201 * t305 - t267;
t24 = t283 * t337 + t267;
t23 = (t329 + (-t356 - t445) * qJD(4)) * t202 + (-t179 + t236) * t336;
t22 = t28 + t444;
t18 = mrSges(6,1) * t76 - mrSges(6,3) * t75;
t13 = t75 * Ifges(5,4) - t76 * Ifges(5,2) + t141 * Ifges(5,6);
t10 = (qJ(6) * qJD(4) - qJD(3) * t179) * t348 + (qJD(6) * t202 + (qJ(6) * qJD(3) - qJD(4) * t179) * t205) * t201 + t345 + t444;
t9 = (-qJ(6) * t336 - t150) * t204 + (qJ(6) * t333 - t327 + (-pkin(5) + t283) * qJD(3)) * t202 + t346;
t8 = pkin(4) * t76 - t211;
t3 = -t411 * t76 + qJDD(6) + t211;
t11 = [-(t201 * t427 + t204 * t62) * t332 / 0.2e1 + (-t13 / 0.2e1 - t3 * mrSges(7,1) + t459 / 0.2e1 - t6 * mrSges(5,3) - t4 * mrSges(6,2) + t2 * mrSges(7,3)) * t350 + (Ifges(4,1) * t152 + t455 * t286 + t242 * t410 + t253 * t413 + t440 * t409 + t441 * t412 + t422 * t414 + Ifges(4,5) * qJDD(3) + (m(7) * t384 - mrSges(7,3)) * t186 * g(1) + t8 * t261) * t202 + (Ifges(7,3) * t410 - Ifges(5,6) * t413 + (-Ifges(4,2) * t202 + t375) * t475 - t434 / 0.2e1 - t483 * t409 - t480 * t412 + t493 * t414 + t463) * t205 + (-t117 * t336 + t433) * mrSges(4,3) + (-t241 * t403 - t252 * t408 - t404 * t438 - t405 * t421 - t407 * t439) * t332 + (Ifges(4,4) * t152 + Ifges(4,2) * t151 + 0.2e1 * Ifges(4,6) * qJDD(3) + t227) * t492 - t229 * t133 + (-mrSges(5,1) * t271 + mrSges(6,1) * t41 - mrSges(7,1) * t25 - mrSges(6,2) * t31 - mrSges(5,3) * t43 + mrSges(7,3) * t21) * (t202 * t331 + t303) + (-mrSges(5,2) * t271 + mrSges(6,2) * t30 + mrSges(7,2) * t25 - mrSges(5,3) * t42 - mrSges(6,3) * t41 - mrSges(7,3) * t17) * (-t201 * t332 + t205 * t326) + (m(5) * (-t202 * t229 - t271 * t336) + t446 * t336 + t451 * t202 + m(4) * ((-t117 * t205 - t118 * t202) * qJD(3) + t433) + t205 * t120) * t179 + (m(4) * t180 - t266) * t155 + m(6) * (t22 * t31 + t24 * t30 + t39 * t41 + t4 * t78 + t5 * t79 + t8 * t92) + m(7) * (t1 * t44 + t10 * t21 + t17 * t9 + t2 * t53 + t23 * t25 + t3 * t80) + t53 * t36 + t44 * t33 + (-t180 * mrSges(4,1) + t488 + t491) * t151 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t200 - 0.2e1 * mrSges(3,2) * t199 + m(3) * (t199 ^ 2 + t200 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (m(3) * t397 + mrSges(2,1) * t203 + mrSges(2,2) * t206 + (-m(5) - m(4)) * t285 + t462 * (-t112 * pkin(4) - qJ(5) * t111 + t285) + t461 * t187 + t416 * t112 + t429 * t111 + (-m(7) * t307 + m(4) * pkin(2) + (-m(5) - m(6)) * t237 - t432 - t453) * t186) * g(1) + t303 * t415 + t230 * t475 + t3 * t173 - t159 * t305 + (-m(3) * t197 - m(4) * t308 - m(5) * t238 - mrSges(2,1) * t206 + mrSges(2,2) * t203 + t462 * (t114 * pkin(4) + qJ(5) * t113 + t238) + t453 * t187 + t461 * t186 + t477 * t352 - t416 * t114 - t429 * t113) * g(2) + (t5 * mrSges(6,2) - t7 * mrSges(5,3) - t1 * mrSges(7,3) + t474) * t348 + (t42 * mrSges(5,1) - t30 * mrSges(6,1) - t17 * mrSges(7,1) - t43 * mrSges(5,2) + t21 * mrSges(7,2) - t118 * mrSges(4,3) + t31 * mrSges(6,3) - t464) * t337 + t465 * t336 / 0.2e1 + m(5) * (t28 * t43 + t29 * t42 + t6 * t82 + t7 * t81) + t78 * t38 + t79 * t35 + t80 * t19 + t81 * t34 + t82 * t37 + t39 * t84 + t23 * t85 + t92 * t18 + t10 * t95 + t28 * t96 + t9 * t97 + t29 * t98 + t24 * t99 + t22 * t100 + (t213 * t403 + t217 * t408 + t247 * t487 + t443 * t404 + t424 * t405 + t442 * t407 + t231) * qJD(3) + (t375 / 0.2e1 + t180 * mrSges(4,2)) * t152; m(3) * qJDD(2) + (-m(3) - m(4) + t430) * g(3) + (-t18 + t19 + (t201 * t323 + t204 * t322 + t159) * qJD(3) + m(4) * (qJD(3) * t118 + t56) + m(5) * (t326 * t43 - t338 * t42 + t229) + m(6) * (t30 * t338 + t31 * t326 - t8) + m(7) * (t17 * t338 + t21 * t326 + t3) - t451) * t205 + (t120 + (t36 + t457) * t204 + (t33 + t458) * t201 + (t383 + t446) * qJD(3) + (-t201 * t322 + t204 * t323) * qJD(4) + m(4) * (-qJD(3) * t117 + t55) + m(5) * (-qJD(3) * t271 + t469) + m(6) * (qJD(3) * t41 + t468) + m(7) * (-qJD(3) * t25 - t418)) * t202; (t415 + t455 / 0.2e1) * t333 + t427 * t286 + (-t319 - t67) * t96 + t421 * t414 + t423 * t340 + (g(3) * t202 + t418) * mrSges(7,3) + (t242 / 0.2e1 - t440 / 0.2e1) * qJD(4) * t174 + (-m(5) * t344 - m(6) * t279 - m(7) * (t279 - t354) - t266 + (-t204 * t426 + t476) * t205 + t432) * g(3) + (pkin(3) * t229 - t42 * t66 - t43 * t67) * m(5) + t229 * t264 + (m(5) * t271 + t316 - t446) * t118 - t55 * mrSges(4,2) + t56 * mrSges(4,1) - pkin(3) * t20 + (t153 * t8 - t30 * t50 - t31 * t48 + t447 * t41) * m(6) + (m(5) * ((-t201 * t43 - t204 * t42) * qJD(4) + t268) + t454 * t331 + t457 * t204 + t458 * t201 + t351 * t430 * g(1) + m(6) * ((-t201 * t31 + t204 * t30) * qJD(4) + t269)) * pkin(8) + t3 * (t204 * mrSges(7,1) + t365) + (-t319 - t48) * t100 + (t315 - t159) * t117 + (-t253 / 0.2e1 + t441 / 0.2e1) * t335 + t241 * t410 + t252 * t413 + t201 * t474 + (t393 * t430 + t425) * t472 + t425 * t473 + ((t217 / 0.2e1 - t442 / 0.2e1) * t145 + t424 * t406 - t231 - t17 * (-t202 * mrSges(7,1) - mrSges(7,3) * t347) - t42 * (mrSges(5,1) * t202 - mrSges(5,3) * t347) - t30 * (-mrSges(6,1) * t202 + mrSges(6,2) * t347) - t43 * (-mrSges(5,2) * t202 - mrSges(5,3) * t349) - t31 * (-mrSges(6,2) * t349 + mrSges(6,3) * t202) - t21 * (mrSges(7,2) * t202 + mrSges(7,3) * t349) + (t443 / 0.2e1 - t213 / 0.2e1) * t174 - t230 * qJD(1) / 0.2e1) * qJD(1) + t62 * t306 / 0.2e1 + (t405 * t422 - t423) * qJD(4) + t464 * t341 - (-Ifges(4,2) * t341 + t185 + t465) * t340 / 0.2e1 - t8 * t262 - t247 * t325 / 0.2e1 + t468 * mrSges(6,2) + t469 * mrSges(5,3) + Ifges(4,3) * qJDD(3) + t431 * t265 + t438 * t409 + t439 * t412 - t66 * t98 - t50 * t99 + t447 * t84 + t448 * t85 + t449 * t95 + t450 * t97 + (t1 * t160 + t135 * t3 + t162 * t2 + t17 * t450 + t21 * t449 + t25 * t448) * m(7) + t135 * t19 - t459 * t204 / 0.2e1 + Ifges(4,6) * t151 + Ifges(4,5) * t152 + t153 * t18 + t160 * t33 + t162 * t36 + t204 * t13 / 0.2e1; (-Ifges(5,2) * t146 - t138 + t427) * t407 + (-t145 * t485 - t481 * t146) * t403 + (t482 * t146 - t466) * t408 + t271 * (mrSges(5,1) * t146 - mrSges(5,2) * t145) + (qJ(5) * t2 - t1 * t411 - t17 * t27 + t21 * t437 - t25 * t47) * m(7) - t411 * t33 - pkin(4) * t35 - t463 + (-Ifges(7,5) * t145 + Ifges(7,6) * t146) * t404 + t62 * t405 + t31 * t381 + t30 * t382 - t227 - t21 * t377 - t17 * t378 + (-m(6) * t31 - t100 - t380 - t96) * t42 + (-t479 * t145 - t374 + t455 + t471) * t406 + (-pkin(4) * t5 + qJ(5) * t4 + qJD(5) * t31 - t41 * t83) * m(6) + (t36 + t38) * qJ(5) - t83 * t84 - t47 * t85 - t26 * t95 - t27 * t97 + t359 * qJD(5) + (-m(6) * t30 + t379 - t454) * t43 - t41 * (mrSges(6,1) * t146 + mrSges(6,3) * t145) - t25 * (-mrSges(7,1) * t146 - mrSges(7,2) * t145) + (t462 * (-t113 * pkin(4) + qJ(5) * t114) - t429 * t114 + t416 * t113) * g(1) + (t462 * (-t111 * pkin(4) + qJ(5) * t112) - t429 * t112 + t416 * t111) * g(2) + (-(t361 + (-m(6) * pkin(4) - mrSges(6,1)) * t201) * t202 + t133 - t278 * t350 - t173 + t462 * t171) * g(3) + t434; t383 * t146 + t359 * t174 + t33 + t35 + (-t146 * t25 + t174 * t21 + t1 + t223) * m(7) + (t146 * t41 + t174 * t31 + t223 + t5) * m(6); -t145 * t95 + t146 * t97 + (-g(3) * t205 - t21 * t145 + t17 * t146 + t202 * t431 + t3) * m(7) + t19;];
tau  = t11;
