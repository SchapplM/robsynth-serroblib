% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:22
% EndTime: 2019-03-09 04:47:05
% DurationCPUTime: 28.77s
% Computational Cost: add. (7075->667), mult. (14159->852), div. (0->0), fcn. (8909->10), ass. (0->316)
t246 = -pkin(1) - pkin(7);
t215 = qJD(1) * t246 + qJD(2);
t244 = cos(qJ(3));
t354 = t215 * t244;
t168 = -qJD(3) * pkin(3) - t354;
t240 = sin(qJ(4));
t243 = cos(qJ(4));
t331 = qJD(3) * t243;
t333 = qJD(1) * t244;
t184 = -t240 * t333 + t331;
t117 = -pkin(4) * t184 + qJD(5) + t168;
t241 = sin(qJ(3));
t334 = qJD(1) * t241;
t219 = qJD(4) + t334;
t237 = sin(pkin(9));
t289 = pkin(3) * t241 - pkin(8) * t244;
t197 = qJ(2) + t289;
t160 = t197 * qJD(1);
t196 = t241 * t215;
t167 = qJD(3) * pkin(8) + t196;
t100 = t243 * t160 - t167 * t240;
t185 = qJD(3) * t240 + t243 * t333;
t77 = -qJ(5) * t185 + t100;
t70 = pkin(4) * t219 + t77;
t238 = cos(pkin(9));
t101 = t160 * t240 + t167 * t243;
t78 = qJ(5) * t184 + t101;
t72 = t238 * t78;
t24 = t237 * t70 + t72;
t22 = qJ(6) * t219 + t24;
t112 = -t238 * t184 + t185 * t237;
t268 = t184 * t237 + t238 * t185;
t29 = pkin(5) * t112 - qJ(6) * t268 + t117;
t386 = t219 / 0.2e1;
t387 = -t219 / 0.2e1;
t395 = t268 / 0.2e1;
t396 = -t268 / 0.2e1;
t399 = -t112 / 0.2e1;
t502 = Ifges(6,4) * t395 + Ifges(7,5) * t396 + Ifges(6,6) * t386 + Ifges(7,6) * t387 + (Ifges(6,2) + Ifges(7,3)) * t399 - mrSges(6,1) * t117 - t29 * mrSges(7,1) + mrSges(7,2) * t22 + mrSges(6,3) * t24;
t398 = t112 / 0.2e1;
t446 = Ifges(6,6) - Ifges(7,6);
t448 = Ifges(6,4) - Ifges(7,5);
t501 = -Ifges(6,2) * t398 + Ifges(7,3) * t399 - t387 * t446 - t396 * t448 + t502;
t500 = Ifges(6,2) * t399 - Ifges(7,3) * t398 + t395 * t448 + t502;
t470 = Ifges(7,2) + Ifges(6,3);
t488 = Ifges(5,3) + t470;
t499 = -t448 + Ifges(7,5);
t447 = Ifges(7,4) + Ifges(6,5);
t498 = t447 * t396;
t322 = qJD(1) * qJD(3);
t191 = qJDD(1) * t244 - t241 * t322;
t107 = qJD(4) * t184 + qJDD(3) * t240 + t191 * t243;
t108 = -qJD(4) * t185 + qJDD(3) * t243 - t191 * t240;
t57 = t107 * t237 - t238 * t108;
t407 = t57 / 0.2e1;
t58 = t107 * t238 + t108 * t237;
t496 = -t58 / 0.2e1;
t192 = -qJDD(1) * t241 - t244 * t322;
t178 = qJDD(4) - t192;
t495 = -t178 / 0.2e1;
t494 = -t184 / 0.2e1;
t493 = -t185 / 0.2e1;
t492 = mrSges(6,1) + mrSges(7,1);
t491 = -mrSges(6,2) + mrSges(7,3);
t449 = Ifges(6,1) + Ifges(7,1);
t290 = pkin(3) * t244 + pkin(8) * t241;
t189 = t290 * qJD(1);
t341 = t243 * t244;
t120 = t240 * t189 + t215 * t341;
t239 = -qJ(5) - pkin(8);
t294 = qJD(4) * t239;
t310 = t240 * t334;
t325 = qJD(5) * t243;
t490 = -qJ(5) * t310 + t240 * t294 - t120 + t325;
t349 = t240 * t244;
t119 = t243 * t189 - t215 * t349;
t346 = t241 * t243;
t489 = -qJD(5) * t240 + t243 * t294 - (pkin(4) * t244 + qJ(5) * t346) * qJD(1) - t119;
t420 = -t117 * mrSges(6,2) + t29 * mrSges(7,3);
t438 = -t112 * t448 + t219 * t447 + t268 * t449;
t362 = t237 * t78;
t23 = t238 * t70 - t362;
t21 = -pkin(5) * t219 + qJD(6) - t23;
t463 = -mrSges(7,2) * t21 + mrSges(6,3) * t23;
t456 = -t438 / 0.2e1 + t463;
t487 = Ifges(6,4) * t399 + Ifges(7,5) * t398 + t447 * t386 + t449 * t395 - t420 - t456;
t486 = m(6) + m(7);
t485 = t191 / 0.2e1;
t484 = t192 / 0.2e1;
t483 = -pkin(4) * t486 - mrSges(5,1);
t323 = qJD(1) * qJD(2);
t216 = qJDD(1) * qJ(2) + t323;
t109 = -pkin(3) * t192 - pkin(8) * t191 + t216;
t211 = qJDD(1) * t246 + qJDD(2);
t330 = qJD(3) * t244;
t127 = t241 * t211 + t215 * t330;
t123 = qJDD(3) * pkin(8) + t127;
t28 = -qJD(4) * t101 + t243 * t109 - t123 * t240;
t12 = pkin(4) * t178 - qJ(5) * t107 - qJD(5) * t185 + t28;
t327 = qJD(4) * t243;
t328 = qJD(4) * t240;
t27 = t240 * t109 + t243 * t123 + t160 * t327 - t167 * t328;
t19 = qJ(5) * t108 + qJD(5) * t184 + t27;
t4 = t237 * t12 + t238 * t19;
t1 = qJ(6) * t178 + qJD(6) * t219 + t4;
t391 = t178 / 0.2e1;
t406 = t58 / 0.2e1;
t408 = -t57 / 0.2e1;
t332 = qJD(3) * t241;
t126 = t211 * t244 - t215 * t332;
t122 = -qJDD(3) * pkin(3) - t126;
t67 = -pkin(4) * t108 + qJDD(5) + t122;
t5 = pkin(5) * t57 - qJ(6) * t58 - qJD(6) * t268 + t67;
t479 = mrSges(6,1) * t67 + mrSges(7,1) * t5 - mrSges(7,2) * t1 - mrSges(6,3) * t4 + Ifges(6,4) * t496 + Ifges(6,6) * t495 + 0.2e1 * Ifges(7,3) * t407 + (-t408 + t407) * Ifges(6,2) + t499 * t406 + (-t446 + Ifges(7,6)) * t391;
t368 = Ifges(4,4) * t244;
t278 = -t241 * Ifges(4,2) + t368;
t477 = t21 * mrSges(7,1) + t24 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t278 / 0.2e1 + Ifges(5,5) * t493 + Ifges(5,6) * t494 + t498 + t446 * t398 + t488 * t387 - t22 * mrSges(7,3) - t23 * mrSges(6,1);
t3 = t12 * t238 - t19 * t237;
t2 = -pkin(5) * t178 + qJDD(6) - t3;
t476 = t28 * mrSges(5,1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t27 * mrSges(5,2) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t475 = -mrSges(6,2) * t67 - mrSges(7,2) * t2 + mrSges(6,3) * t3 + mrSges(7,3) * t5 - Ifges(6,4) * t408 + (-t406 + t496) * t449 + (-t391 + t495) * t447 - t499 * t407;
t474 = -m(4) - m(5);
t471 = mrSges(6,3) + mrSges(7,2);
t440 = -t237 * t490 + t238 * t489;
t439 = t237 * t489 + t238 * t490;
t62 = -mrSges(5,1) * t108 + mrSges(5,2) * t107;
t468 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t191 - t62;
t431 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t184 - mrSges(5,2) * t185 - mrSges(4,3) * t333;
t236 = qJ(4) + pkin(9);
t229 = sin(t236);
t230 = cos(t236);
t467 = t229 * t491 + t230 * t492;
t326 = qJD(4) * t244;
t304 = t243 * t326;
t309 = t240 * t332;
t466 = t304 - t309;
t242 = sin(qJ(1));
t343 = t242 * t243;
t245 = cos(qJ(1));
t345 = t241 * t245;
t163 = t240 * t345 + t343;
t340 = t243 * t245;
t347 = t241 * t242;
t161 = -t240 * t347 + t340;
t427 = -t196 + (t328 + t310) * pkin(4);
t270 = t126 * t244 + t127 * t241;
t369 = Ifges(4,4) * t241;
t281 = t244 * Ifges(4,1) - t369;
t174 = Ifges(5,4) * t184;
t98 = t185 * Ifges(5,1) + t219 * Ifges(5,5) + t174;
t465 = Ifges(4,5) * qJD(3) + qJD(1) * t281 + t243 * t98;
t287 = mrSges(4,1) * t244 - mrSges(4,2) * t241;
t460 = t244 * (-Ifges(4,1) * t241 - t368) / 0.2e1 + qJ(2) * t287;
t459 = mrSges(3,2) - mrSges(4,3) - mrSges(2,1);
t457 = Ifges(5,5) * t107 + Ifges(5,6) * t108 + t178 * t488 - t446 * t57 + t447 * t58;
t286 = mrSges(4,1) * t241 + mrSges(4,2) * t244;
t455 = -m(5) * t289 + t244 * mrSges(5,3) + mrSges(2,2) - mrSges(3,3) - t286;
t422 = m(7) * pkin(5) + t492;
t421 = m(7) * qJ(6) + t491;
t401 = t107 / 0.2e1;
t400 = t108 / 0.2e1;
t451 = g(2) * t239 * t345;
t375 = g(2) * t245;
t17 = t57 * mrSges(7,1) - t58 * mrSges(7,3);
t18 = t57 * mrSges(6,1) + t58 * mrSges(6,2);
t445 = -t17 - t18;
t444 = -qJ(6) * t333 + t439;
t31 = -mrSges(7,2) * t57 + mrSges(7,3) * t178;
t32 = -mrSges(6,2) * t178 - mrSges(6,3) * t57;
t443 = t31 + t32;
t33 = mrSges(6,1) * t178 - mrSges(6,3) * t58;
t34 = -t178 * mrSges(7,1) + t58 * mrSges(7,2);
t442 = t34 - t33;
t441 = pkin(5) * t333 - t440;
t181 = t237 * t243 + t238 * t240;
t156 = t181 * qJD(1);
t133 = t241 * t156;
t180 = t237 * t240 - t238 * t243;
t433 = t180 * t241;
t134 = qJD(1) * t433;
t157 = t180 * qJD(4);
t423 = t181 * qJD(4);
t437 = -qJD(6) * t181 + t427 + (t157 + t134) * qJ(6) + (t423 + t133) * pkin(5);
t82 = -mrSges(6,2) * t219 - mrSges(6,3) * t112;
t83 = -mrSges(7,2) * t112 + mrSges(7,3) * t219;
t371 = t82 + t83;
t84 = mrSges(6,1) * t219 - mrSges(6,3) * t268;
t85 = -mrSges(7,1) * t219 + mrSges(7,2) * t268;
t370 = t84 - t85;
t306 = t241 * t328;
t435 = -t238 * t241 * t327 + t180 * qJD(1) - t181 * t330 + t237 * t306;
t145 = t180 * t244;
t434 = -qJD(3) * t145 - t241 * t423 - t156;
t376 = g(1) * t242;
t424 = t375 - t376;
t432 = t244 * t424;
t344 = t241 * t246;
t132 = t240 * t197 + t243 * t344;
t227 = pkin(4) * t243 + pkin(3);
t272 = pkin(5) * t230 + qJ(6) * t229;
t430 = -m(7) * (-t227 - t272) + m(6) * t227 + t467;
t428 = t486 * t239 - t471;
t75 = mrSges(5,1) * t178 - mrSges(5,3) * t107;
t76 = -mrSges(5,2) * t178 + mrSges(5,3) * t108;
t426 = -t240 * t75 + t243 * t76;
t425 = -t240 * t28 + t243 * t27;
t285 = -mrSges(5,1) * t243 + mrSges(5,2) * t240;
t259 = m(5) * pkin(3) - t285;
t317 = m(5) * pkin(8) + mrSges(5,3);
t419 = t241 * t317 + t244 * t259;
t417 = Ifges(6,4) * t398 + Ifges(7,5) * t399 + t387 * t447 + t396 * t449 + t420;
t416 = qJD(1) ^ 2;
t413 = Ifges(5,1) * t401 + Ifges(5,4) * t400 + Ifges(5,5) * t391;
t363 = t185 * Ifges(5,4);
t97 = t184 * Ifges(5,2) + t219 * Ifges(5,6) + t363;
t404 = -t97 / 0.2e1;
t403 = -m(3) - m(4);
t388 = t185 / 0.2e1;
t380 = pkin(4) * t185;
t379 = pkin(4) * t237;
t378 = pkin(4) * t238;
t377 = pkin(4) * t240;
t374 = g(3) * t244;
t179 = qJD(3) * t290 + qJD(2);
t253 = -qJD(4) * t132 + t243 * t179;
t295 = -t240 * t246 + pkin(4);
t308 = t241 * t331;
t42 = qJ(5) * t308 + (qJ(5) * t328 + qJD(3) * t295 - t325) * t244 + t253;
t329 = qJD(3) * t246;
t307 = t244 * t329;
t313 = t240 * t179 + t197 * t327 + t243 * t307;
t60 = -qJ(5) * t304 + (-qJD(5) * t244 + (qJ(5) * qJD(3) - qJD(4) * t246) * t241) * t240 + t313;
t15 = t237 * t42 + t238 * t60;
t367 = Ifges(5,4) * t240;
t366 = Ifges(5,4) * t243;
t365 = t100 * mrSges(5,3);
t364 = t101 * mrSges(5,3);
t177 = t243 * t197;
t110 = -qJ(5) * t341 + t241 * t295 + t177;
t118 = -qJ(5) * t349 + t132;
t64 = t237 * t110 + t238 * t118;
t348 = t240 * t245;
t342 = t242 * t244;
t339 = t244 * t245;
t338 = t244 * t246;
t337 = t245 * t230;
t335 = t245 * pkin(1) + t242 * qJ(2);
t324 = qJDD(1) * mrSges(3,2);
t311 = t245 * pkin(7) + t335;
t218 = t241 * t329;
t305 = t240 * t326;
t301 = t239 * t240;
t293 = -t322 / 0.2e1;
t291 = (t216 + t323) * qJ(2);
t220 = pkin(4) * t349;
t182 = t220 - t338;
t284 = mrSges(5,1) * t240 + mrSges(5,2) * t243;
t280 = Ifges(5,1) * t243 - t367;
t279 = Ifges(5,1) * t240 + t366;
t277 = -Ifges(5,2) * t240 + t366;
t276 = Ifges(5,2) * t243 + t367;
t275 = -Ifges(4,5) * t241 - Ifges(4,6) * t244;
t274 = Ifges(5,5) * t243 - Ifges(5,6) * t240;
t273 = Ifges(5,5) * t240 + Ifges(5,6) * t243;
t14 = -t237 * t60 + t238 * t42;
t271 = t100 * t243 + t101 * t240;
t63 = t110 * t238 - t118 * t237;
t129 = -mrSges(5,2) * t219 + mrSges(5,3) * t184;
t130 = mrSges(5,1) * t219 - mrSges(5,3) * t185;
t269 = -t240 * t129 - t243 * t130;
t263 = t241 * (-Ifges(4,2) * t244 - t369);
t128 = pkin(4) * t466 + t218;
t142 = t181 * t241;
t257 = t305 + t308;
t252 = Ifges(5,5) * t244 - t241 * t280;
t251 = Ifges(5,6) * t244 - t241 * t277;
t250 = Ifges(5,3) * t244 - t241 * t274;
t249 = -qJD(4) * t271 + t425;
t233 = t245 * qJ(2);
t228 = -qJDD(1) * pkin(1) + qJDD(2);
t226 = -pkin(5) - t378;
t224 = qJ(6) + t379;
t206 = t239 * t243;
t203 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t334;
t194 = t227 * t342;
t188 = t286 * qJD(1);
t173 = t284 * t244;
t164 = -t240 * t242 + t241 * t340;
t162 = t241 * t343 + t348;
t148 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t192;
t143 = t181 * t244;
t141 = -t229 * t242 + t241 * t337;
t140 = t229 * t345 + t230 * t242;
t139 = t229 * t245 + t230 * t347;
t138 = t229 * t347 - t337;
t131 = -t240 * t344 + t177;
t125 = -t238 * t206 + t237 * t301;
t124 = -t206 * t237 - t238 * t301;
t106 = pkin(5) * t180 - qJ(6) * t181 - t227;
t91 = qJD(3) * t433 - t244 * t423;
t89 = qJD(3) * t142 + t237 * t305 - t238 * t304;
t80 = -t240 * t307 + t253;
t79 = -t246 * t306 + t313;
t74 = pkin(5) * t143 + qJ(6) * t145 + t182;
t66 = mrSges(6,1) * t112 + mrSges(6,2) * t268;
t65 = mrSges(7,1) * t112 - mrSges(7,3) * t268;
t61 = -pkin(5) * t241 - t63;
t59 = qJ(6) * t241 + t64;
t40 = pkin(5) * t268 + qJ(6) * t112 + t380;
t38 = t107 * Ifges(5,4) + t108 * Ifges(5,2) + t178 * Ifges(5,6);
t26 = t238 * t77 - t362;
t25 = t237 * t77 + t72;
t20 = -pkin(5) * t89 - qJ(6) * t91 + qJD(6) * t145 + t128;
t13 = -pkin(5) * t330 - t14;
t11 = qJ(6) * t330 + qJD(6) * t241 + t15;
t6 = [-t431 * t218 - t465 * t332 / 0.2e1 - t270 * mrSges(4,3) + t460 * t322 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t479 * t143 + (t457 / 0.2e1 + Ifges(7,6) * t407 + Ifges(6,6) * t408 + t488 * t391 + t447 * t406 + Ifges(5,5) * t401 + Ifges(5,6) * t400 - Ifges(4,6) * qJDD(3) - Ifges(4,4) * t191 / 0.2e1 - Ifges(4,2) * t192 / 0.2e1 + t476) * t241 + (-m(5) * t122 * t246 + Ifges(4,1) * t485 + Ifges(4,4) * t484 + Ifges(4,5) * qJDD(3) + t274 * t391 + t277 * t400 + t280 * t401) * t244 - (t240 * t98 + t243 * t97) * t326 / 0.2e1 + t475 * t145 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(6,6) * t399 + Ifges(7,6) * t398 + t470 * t386 + t447 * t395 - t477) * t330 + m(5) * (t168 * t246 * t332 + t100 * t80 + t101 * t79 + t131 * t28 + t132 * t27) + t468 * t338 + (t286 + 0.2e1 * mrSges(3,3)) * t216 + t500 * t89 + m(4) * (t246 * t270 + t291) + m(3) * (-pkin(1) * t228 + t291) + t278 * t484 + t281 * t485 + t203 * t307 + (-m(3) * t335 - t162 * mrSges(5,1) - t161 * mrSges(5,2) + t471 * t342 + t474 * t311 - t486 * (pkin(4) * t348 + t227 * t347 + t239 * t342 + t311) - t422 * t139 - t421 * t138 + t459 * t245 + t455 * t242) * g(2) + (-t164 * mrSges(5,1) + t163 * mrSges(5,2) - t486 * ((-t377 + t246) * t242 + t227 * t345 + t239 * t339 + t233) + t471 * t339 - t422 * t141 - t421 * t140 + (-m(3) + t474) * t233 + (m(3) * pkin(1) + t474 * t246 - t459) * t242 + t455 * t245) * g(1) - t38 * t349 / 0.2e1 + (qJD(3) * t252 - t279 * t326) * t388 + t148 * t344 + m(7) * (t1 * t59 + t11 * t22 + t13 * t21 + t2 * t61 + t20 * t29 + t5 * t74) + m(6) * (t117 * t128 + t14 * t23 + t15 * t24 + t182 * t67 + t3 * t63 + t4 * t64) + t263 * t293 + t341 * t413 + t487 * t91 + (qJD(3) * t250 - t273 * t326 + t446 * t89) * t386 + t97 * t309 / 0.2e1 - pkin(1) * t324 + t184 * (qJD(3) * t251 - t276 * t326) / 0.2e1 + t228 * mrSges(3,2) + qJ(2) * (-mrSges(4,1) * t192 + mrSges(4,2) * t191) + qJD(2) * t188 + t182 * t18 + t122 * t173 + t128 * t66 + t79 * t129 + t80 * t130 + t131 * t75 + t132 * t76 + t168 * (mrSges(5,1) * t466 - mrSges(5,2) * t257) + (t100 * t257 - t101 * t466 - t27 * t349 - t28 * t341) * mrSges(5,3) + qJD(3) ^ 2 * t275 / 0.2e1 + t59 * t31 + t61 * t34 + t63 * t33 + t64 * t32 + t20 * t65 + t74 * t17 + t15 * t82 + t11 * t83 + t14 * t84 + t13 * t85; t324 - t443 * t433 + t442 * t142 + (qJ(2) * t403 - mrSges(3,3)) * t416 + (-t188 + t269) * qJD(1) + ((t129 * t243 - t130 * t240 + t203) * qJD(3) + t445 + t468) * t244 + (t148 + t269 * qJD(4) + (t65 + t66 - t431) * qJD(3) + t426) * t241 + m(3) * t228 + m(4) * t270 + t434 * t371 + t435 * t370 + t424 * (m(5) + t486 - t403) + (-t1 * t433 + t142 * t2 - t21 * t435 + t22 * t434 - t244 * t5 + t29 * t332) * m(7) + (t117 * t332 - t142 * t3 + t23 * t435 + t24 * t434 - t244 * t67 - t4 * t433) * m(6) + ((-t122 + (-t100 * t240 + t101 * t243) * qJD(3)) * t244 + (qJD(3) * t168 + t249) * t241 - t271 * qJD(1)) * m(5); (-m(7) * t194 + ((-m(7) * t272 - t467) * t244 + t428 * t241) * t242) * g(1) + t465 * t334 / 0.2e1 + (t263 / 0.2e1 - t460) * t416 + t479 * t180 + (t417 + t456) * t134 + (t286 + (-t317 + t428) * t244 + (t259 + t430) * t241) * g(3) - t475 * t181 + (Ifges(6,6) * t398 + Ifges(7,6) * t399 + t470 * t387 + t477 + t498) * t333 + (t241 * t471 + t430 * t244 + t419) * t375 + t431 * t196 + t424 * t287 + t425 * mrSges(5,3) + (m(5) * t249 - t129 * t328 - t130 * t327 + t426) * pkin(8) + t427 * t66 + (-pkin(3) * t122 - t100 * t119 - t101 * t120 - t168 * t196) * m(5) + (-t386 * t446 - t500) * t423 - t501 * t133 + t122 * t285 + (t404 - t364) * t328 + (t184 * t277 + t185 * t280 + t219 * t274) * qJD(4) / 0.2e1 - (t184 * t251 + t185 * t252 + t219 * t250) * qJD(1) / 0.2e1 + t275 * t293 + t276 * t400 + t279 * t401 + t310 * t404 + t240 * t413 + t273 * t391 + (t98 / 0.2e1 - t365) * t327 + (-t101 * (mrSges(5,3) * t240 * t241 - mrSges(5,2) * t244) - t100 * (mrSges(5,1) * t244 + mrSges(5,3) * t346)) * qJD(1) - t419 * t376 - t487 * t157 + t219 * t284 * t168 - t203 * t354 + t243 * t38 / 0.2e1 - t227 * t18 + Ifges(4,5) * t191 + Ifges(4,6) * t192 + Ifges(4,3) * qJDD(3) + t126 * mrSges(4,1) - t127 * mrSges(4,2) - t120 * t129 - t119 * t130 + t437 * t65 + t439 * t82 + t440 * t84 + t441 * t85 + t442 * t124 + t443 * t125 + t444 * t83 + (-g(1) * t194 + t117 * t427 - t124 * t3 + t125 * t4 - t227 * t67 + t23 * t440 + t24 * t439 - t451) * m(6) + (t1 * t125 + t106 * t5 + t124 * t2 + t21 * t441 + t22 * t444 + t29 * t437 - t451) * m(7) - pkin(3) * t62 + t106 * t17; -(t417 + t463) * t112 + t476 + (Ifges(5,1) * t184 - t363) * t493 + (-Ifges(5,2) * t185 + t174 + t98) * t494 + (-mrSges(5,2) * t164 - t422 * t140 + t421 * t141 + t483 * t163) * g(2) + (mrSges(5,2) * t162 + t422 * t138 - t421 * t139 + t483 * t161) * g(1) + t457 + t438 * t398 + t501 * t268 + (-(-t229 * mrSges(7,1) + t230 * mrSges(7,3)) * t244 + t173) * g(3) + ((t237 * t4 + t238 * t3) * pkin(4) + t377 * t374 - t117 * t380 + t23 * t25 - t24 * t26) * m(6) + t33 * t378 + t32 * t379 + (Ifges(5,5) * t184 - Ifges(5,6) * t185) * t387 + t185 * t364 + t184 * t365 + (-t21 * t25 - t29 * t40 + t1 * t224 + t2 * t226 - g(3) * (-t220 + (-pkin(5) * t229 + qJ(6) * t230) * t244) + (qJD(6) - t26) * t22) * m(7) + t97 * t388 - (-mrSges(6,1) * t229 - mrSges(6,2) * t230) * t374 + t224 * t31 + t226 * t34 - t168 * (mrSges(5,1) * t185 + mrSges(5,2) * t184) - t100 * t129 + t101 * t130 + t370 * t25 - t371 * t26 - t66 * t380 - t40 * t65 + qJD(6) * t83; -t486 * t241 * g(3) + t370 * t268 + t371 * t112 + (t112 * t22 - t21 * t268 - t432 + t5) * m(7) + (t112 * t24 + t23 * t268 - t432 + t67) * m(6) - t445; t268 * t65 - t219 * t83 + (-g(1) * t138 + g(2) * t140 - t22 * t219 - t229 * t374 + t268 * t29 + t2) * m(7) + t34;];
tau  = t6;
