% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:16
% EndTime: 2019-03-09 06:15:06
% DurationCPUTime: 30.36s
% Computational Cost: add. (14960->747), mult. (35779->943), div. (0->0), fcn. (27347->14), ass. (0->322)
t520 = mrSges(6,1) + mrSges(7,1);
t519 = mrSges(6,2) + mrSges(7,2);
t484 = Ifges(6,4) + Ifges(7,4);
t485 = Ifges(6,1) + Ifges(7,1);
t483 = Ifges(6,5) + Ifges(7,5);
t482 = Ifges(6,2) + Ifges(7,2);
t481 = Ifges(6,6) + Ifges(7,6);
t505 = Ifges(7,3) + Ifges(6,3);
t279 = sin(pkin(10));
t280 = cos(pkin(10));
t284 = sin(qJ(3));
t288 = cos(qJ(3));
t498 = -t279 * t284 + t288 * t280;
t228 = t498 * qJD(1);
t238 = t279 * t288 + t280 * t284;
t229 = t238 * qJD(1);
t173 = pkin(3) * t229 - pkin(8) * t228;
t403 = pkin(7) + qJ(2);
t253 = t403 * t279;
t239 = qJD(1) * t253;
t254 = t403 * t280;
t240 = qJD(1) * t254;
t178 = -t239 * t288 - t284 * t240;
t283 = sin(qJ(4));
t287 = cos(qJ(4));
t118 = t283 * t173 + t287 * t178;
t290 = -pkin(9) - pkin(8);
t333 = qJD(4) * t290;
t380 = t228 * t283;
t524 = pkin(9) * t380 + t283 * t333 - t118;
t117 = t287 * t173 - t178 * t283;
t379 = t228 * t287;
t523 = -pkin(4) * t229 + pkin(9) * t379 + t287 * t333 - t117;
t344 = qJD(1) * qJD(2);
t261 = qJ(2) * qJDD(1) + t344;
t522 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t231 = t238 * qJD(3);
t341 = qJDD(1) * t280;
t342 = qJDD(1) * t279;
t177 = -qJD(1) * t231 - t284 * t342 + t288 * t341;
t170 = qJDD(4) - t177;
t167 = qJDD(5) + t170;
t230 = t498 * qJD(3);
t176 = qJD(1) * t230 + qJDD(1) * t238;
t189 = qJD(3) * t287 - t229 * t283;
t126 = qJD(4) * t189 + qJDD(3) * t283 + t176 * t287;
t190 = qJD(3) * t283 + t229 * t287;
t127 = -qJD(4) * t190 + qJDD(3) * t287 - t176 * t283;
t282 = sin(qJ(5));
t286 = cos(qJ(5));
t320 = t286 * t189 - t190 * t282;
t53 = qJD(5) * t320 + t126 * t286 + t127 * t282;
t136 = t189 * t282 + t190 * t286;
t54 = -qJD(5) * t136 - t126 * t282 + t127 * t286;
t521 = -t482 * t54 / 0.2e1 - t484 * t53 / 0.2e1 - t481 * t167 / 0.2e1;
t518 = t484 * t320;
t303 = t282 * t283 - t286 * t287;
t149 = t303 * t228;
t459 = qJD(4) + qJD(5);
t181 = t459 * t303;
t468 = -t181 + t149;
t242 = t282 * t287 + t283 * t286;
t148 = t242 * t228;
t182 = t459 * t242;
t467 = -t182 + t148;
t352 = t279 ^ 2 + t280 ^ 2;
t517 = t484 * t136;
t171 = -qJD(3) * pkin(3) - t178;
t388 = t229 * mrSges(4,3);
t466 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t189 - mrSges(5,2) * t190 - t388;
t516 = -m(5) * t171 + t466;
t265 = pkin(2) * t280 + pkin(1);
t250 = -qJD(1) * t265 + qJD(2);
t471 = Ifges(4,5) * qJD(3);
t515 = -t250 * mrSges(4,2) + t178 * mrSges(4,3) - t471 / 0.2e1;
t144 = -pkin(3) * t228 - pkin(8) * t229 + t250;
t179 = -t284 * t239 + t288 * t240;
t172 = qJD(3) * pkin(8) + t179;
t100 = t144 * t283 + t172 * t287;
t219 = qJD(4) - t228;
t211 = qJD(5) + t219;
t99 = t287 * t144 - t172 * t283;
t84 = -pkin(9) * t190 + t99;
t73 = pkin(4) * t219 + t84;
t85 = pkin(9) * t189 + t100;
t77 = t282 * t85;
t32 = t286 * t73 - t77;
t500 = qJ(6) * t136;
t24 = t32 - t500;
t22 = pkin(5) * t211 + t24;
t79 = t286 * t85;
t33 = t282 * t73 + t79;
t469 = qJ(6) * t320;
t25 = t33 + t469;
t470 = Ifges(4,6) * qJD(3);
t514 = -t250 * mrSges(4,1) - t99 * mrSges(5,1) - t32 * mrSges(6,1) - t22 * mrSges(7,1) + t100 * mrSges(5,2) + t33 * mrSges(6,2) + t25 * mrSges(7,2) + t470 / 0.2e1;
t278 = qJ(4) + qJ(5);
t273 = cos(t278);
t411 = pkin(4) * t287;
t252 = pkin(5) * t273 + t411;
t246 = pkin(3) + t252;
t268 = pkin(3) + t411;
t272 = sin(t278);
t313 = -mrSges(5,1) * t287 + mrSges(5,2) * t283;
t513 = -m(5) * pkin(3) - m(6) * t268 - m(7) * t246 + t519 * t272 - t273 * t520 + t313;
t274 = -qJ(6) + t290;
t512 = -m(5) * pkin(8) + m(6) * t290 + m(7) * t274 - t522;
t455 = m(6) * pkin(4);
t476 = t136 * t485 + t483 * t211 + t518;
t510 = t476 / 0.2e1;
t477 = t211 * t481 + t320 * t482 + t517;
t509 = t477 / 0.2e1;
t503 = t483 * t167 + t484 * t54 + t485 * t53;
t502 = t189 * Ifges(5,6);
t501 = t219 * Ifges(5,3);
t256 = t290 * t283;
t257 = t290 * t287;
t346 = qJD(5) * t286;
t347 = qJD(5) * t282;
t475 = t256 * t346 + t257 * t347 + t282 * t523 + t286 * t524;
t188 = t282 * t256 - t286 * t257;
t474 = -qJD(5) * t188 - t282 * t524 + t286 * t523;
t348 = qJD(4) * t287;
t330 = t238 * t348;
t302 = t230 * t283 + t330;
t349 = qJD(4) * t283;
t464 = -t179 + (t349 - t380) * pkin(4);
t499 = t167 * t505 + t481 * t54 + t483 * t53;
t285 = sin(qJ(1));
t289 = cos(qJ(1));
t497 = g(1) * t289 + g(2) * t285;
t128 = -pkin(4) * t189 + t171;
t76 = -pkin(5) * t320 + qJD(6) + t128;
t496 = -t76 * mrSges(7,1) + mrSges(6,3) * t33 + mrSges(7,3) * t25 + t509;
t328 = m(3) * qJ(2) + mrSges(3,3);
t495 = mrSges(2,2) - t328 - mrSges(4,3);
t494 = t190 * Ifges(5,5) + t136 * t483 + t211 * t505 + t320 * t481 + t501 + t502;
t493 = -m(6) - m(5) - m(7) - m(4);
t413 = pkin(4) * t283;
t251 = pkin(5) * t272 + t413;
t490 = -m(6) * t413 - m(7) * t251;
t249 = -qJDD(1) * t265 + qJDD(2);
t108 = -pkin(3) * t177 - pkin(8) * t176 + t249;
t322 = pkin(7) * qJDD(1) + t261;
t213 = t322 * t279;
t214 = t322 * t280;
t350 = qJD(3) * t288;
t351 = qJD(3) * t284;
t115 = -t284 * t213 + t288 * t214 - t239 * t350 - t240 * t351;
t109 = qJDD(3) * pkin(8) + t115;
t29 = t283 * t108 + t287 * t109 + t144 * t348 - t172 * t349;
t30 = -t100 * qJD(4) + t287 * t108 - t109 * t283;
t489 = t30 * mrSges(5,1) - t29 * mrSges(5,2);
t277 = pkin(10) + qJ(3);
t270 = sin(t277);
t271 = cos(t277);
t315 = mrSges(4,1) * t271 - mrSges(4,2) * t270;
t316 = -mrSges(3,1) * t280 + mrSges(3,2) * t279;
t488 = m(3) * pkin(1) + t270 * t522 + mrSges(2,1) + t315 - t316;
t19 = pkin(4) * t170 - pkin(9) * t126 + t30;
t23 = pkin(9) * t127 + t29;
t6 = -qJD(5) * t33 + t286 * t19 - t23 * t282;
t2 = pkin(5) * t167 - qJ(6) * t53 - qJD(6) * t136 + t6;
t5 = t282 * t19 + t286 * t23 + t73 * t346 - t347 * t85;
t3 = qJ(6) * t54 + qJD(6) * t320 + t5;
t487 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t448 = t126 / 0.2e1;
t447 = t127 / 0.2e1;
t435 = t170 / 0.2e1;
t445 = -t320 / 0.2e1;
t480 = -pkin(5) * t229 - qJ(6) * t468 - qJD(6) * t242 + t474;
t479 = qJ(6) * t467 - qJD(6) * t303 + t475;
t37 = -mrSges(7,2) * t167 + mrSges(7,3) * t54;
t38 = -mrSges(6,2) * t167 + mrSges(6,3) * t54;
t478 = t37 + t38;
t473 = -pkin(5) * t467 + t464;
t472 = t455 + mrSges(5,1);
t162 = t303 * t238;
t175 = -pkin(3) * t498 - pkin(8) * t238 - t265;
t185 = -t253 * t284 + t254 * t288;
t180 = t287 * t185;
t125 = t283 * t175 + t180;
t465 = -t288 * t253 - t254 * t284;
t370 = t273 * t285;
t371 = t272 * t289;
t206 = -t271 * t371 + t370;
t369 = t273 * t289;
t372 = t272 * t285;
t207 = t271 * t369 + t372;
t463 = -t206 * t520 + t207 * t519;
t204 = t271 * t372 + t369;
t205 = -t271 * t370 + t371;
t462 = t204 * t520 - t205 * t519;
t461 = -t283 * t30 + t287 * t29;
t460 = mrSges(6,1) * t272 + t273 * t519;
t456 = t76 * mrSges(7,2) - t32 * mrSges(6,3) - t22 * mrSges(7,3);
t454 = m(7) * pkin(5);
t453 = t53 / 0.2e1;
t452 = t54 / 0.2e1;
t451 = Ifges(5,1) * t448 + Ifges(5,4) * t447 + Ifges(5,5) * t435;
t444 = t320 / 0.2e1;
t442 = -t136 / 0.2e1;
t441 = t136 / 0.2e1;
t436 = t167 / 0.2e1;
t432 = -t189 / 0.2e1;
t431 = -t190 / 0.2e1;
t430 = t190 / 0.2e1;
t429 = -t211 / 0.2e1;
t428 = t211 / 0.2e1;
t427 = -t219 / 0.2e1;
t425 = t228 / 0.2e1;
t423 = t229 / 0.2e1;
t415 = pkin(4) * t190;
t412 = pkin(4) * t286;
t407 = g(3) * t270;
t41 = t286 * t84 - t77;
t376 = t238 * t283;
t101 = -pkin(9) * t376 + t125;
t124 = t287 * t175 - t185 * t283;
t375 = t238 * t287;
t94 = -pkin(4) * t498 - pkin(9) * t375 + t124;
t58 = t286 * t101 + t282 * t94;
t400 = mrSges(5,3) * t189;
t399 = mrSges(6,3) * t320;
t398 = mrSges(6,3) * t136;
t397 = mrSges(7,3) * t320;
t396 = mrSges(7,3) * t136;
t395 = Ifges(5,4) * t283;
t394 = Ifges(5,4) * t287;
t393 = t100 * mrSges(5,3);
t389 = t190 * Ifges(5,4);
t387 = t229 * Ifges(4,4);
t377 = t230 * t287;
t374 = t251 * t285;
t120 = t189 * Ifges(5,2) + t219 * Ifges(5,6) + t389;
t365 = t283 * t120;
t364 = t283 * t285;
t363 = t283 * t289;
t362 = t285 * t287;
t186 = Ifges(5,4) * t189;
t121 = t190 * Ifges(5,1) + t219 * Ifges(5,5) + t186;
t361 = t287 * t121;
t360 = t287 * t289;
t337 = pkin(4) * t347;
t336 = pkin(4) * t346;
t334 = Ifges(5,5) * t126 + Ifges(5,6) * t127 + Ifges(5,3) * t170;
t329 = t361 / 0.2e1;
t15 = -t54 * mrSges(7,1) + t53 * mrSges(7,2);
t325 = -t349 / 0.2e1;
t40 = -t282 * t84 - t79;
t324 = -t177 * mrSges(4,1) + t176 * mrSges(4,2);
t57 = -t101 * t282 + t286 * t94;
t145 = qJD(2) * t498 + qJD(3) * t465;
t174 = pkin(3) * t231 - pkin(8) * t230;
t321 = -t145 * t283 + t287 * t174;
t187 = t286 * t256 + t257 * t282;
t318 = pkin(3) * t271 + pkin(8) * t270;
t147 = pkin(4) * t376 - t465;
t317 = -mrSges(3,1) * t341 + mrSges(3,2) * t342;
t312 = mrSges(5,1) * t283 + mrSges(5,2) * t287;
t310 = Ifges(5,1) * t287 - t395;
t309 = -Ifges(5,2) * t283 + t394;
t308 = Ifges(5,5) * t287 - Ifges(5,6) * t283;
t142 = -mrSges(5,2) * t219 + t400;
t143 = mrSges(5,1) * t219 - mrSges(5,3) * t190;
t306 = t142 * t287 - t143 * t283;
t305 = t246 * t271 - t270 * t274;
t304 = t268 * t271 - t270 * t290;
t116 = -t213 * t288 - t284 * t214 + t239 * t351 - t240 * t350;
t222 = -t271 * t363 + t362;
t220 = t271 * t364 + t360;
t301 = t238 * t349 - t377;
t44 = -pkin(9) * t377 + pkin(4) * t231 + (-t180 + (pkin(9) * t238 - t175) * t283) * qJD(4) + t321;
t61 = t287 * t145 + t283 * t174 + t175 * t348 - t185 * t349;
t56 = -pkin(9) * t302 + t61;
t9 = -t101 * t347 + t282 * t44 + t286 * t56 + t94 * t346;
t300 = t171 * t312;
t110 = -qJDD(3) * pkin(3) - t116;
t63 = -pkin(4) * t127 + t110;
t10 = -qJD(5) * t58 - t282 * t56 + t286 * t44;
t295 = t487 + t499;
t146 = qJD(2) * t238 + qJD(3) * t185;
t107 = pkin(4) * t302 + t146;
t269 = -qJDD(1) * pkin(1) + qJDD(2);
t267 = pkin(5) + t412;
t223 = t271 * t360 + t364;
t221 = -t271 * t362 + t363;
t215 = Ifges(4,4) * t228;
t210 = pkin(5) * t303 - t268;
t199 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t228;
t161 = t242 * t238;
t154 = t229 * Ifges(4,1) + t215 + t471;
t153 = t228 * Ifges(4,2) + t387 + t470;
t152 = -qJ(6) * t303 + t188;
t151 = -qJ(6) * t242 + t187;
t114 = mrSges(6,1) * t211 - t398;
t113 = mrSges(7,1) * t211 - t396;
t112 = -mrSges(6,2) * t211 + t399;
t111 = -mrSges(7,2) * t211 + t397;
t105 = pkin(5) * t136 + t415;
t104 = pkin(5) * t161 + t147;
t88 = -mrSges(5,2) * t170 + mrSges(5,3) * t127;
t87 = mrSges(5,1) * t170 - mrSges(5,3) * t126;
t81 = -mrSges(6,1) * t320 + mrSges(6,2) * t136;
t80 = -mrSges(7,1) * t320 + mrSges(7,2) * t136;
t75 = t162 * t459 - t242 * t230;
t74 = -t182 * t238 - t230 * t303;
t64 = -mrSges(5,1) * t127 + mrSges(5,2) * t126;
t62 = -qJD(4) * t125 + t321;
t59 = t126 * Ifges(5,4) + t127 * Ifges(5,2) + t170 * Ifges(5,6);
t47 = -pkin(5) * t75 + t107;
t39 = -qJ(6) * t161 + t58;
t36 = mrSges(6,1) * t167 - mrSges(6,3) * t53;
t35 = mrSges(7,1) * t167 - mrSges(7,3) * t53;
t34 = -pkin(5) * t498 + qJ(6) * t162 + t57;
t27 = t41 - t500;
t26 = t40 - t469;
t20 = -pkin(5) * t54 + qJDD(6) + t63;
t16 = -mrSges(6,1) * t54 + mrSges(6,2) * t53;
t8 = qJ(6) * t75 - qJD(6) * t161 + t9;
t7 = pkin(5) * t231 - qJ(6) * t74 + qJD(6) * t162 + t10;
t1 = [(-t161 * t5 + t162 * t6 - t32 * t74 + t33 * t75) * mrSges(6,3) + t75 * t509 + t74 * t510 - (-m(4) * t116 + m(5) * t110 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t176 + t64) * t465 + (t249 * mrSges(4,2) - t116 * mrSges(4,3) + Ifges(4,1) * t176 + Ifges(4,4) * t177 + Ifges(4,5) * qJDD(3) + t110 * t312 + t121 * t325 + t308 * t435 + t309 * t447 + t310 * t448) * t238 + t20 * (mrSges(7,1) * t161 - mrSges(7,2) * t162) + t63 * (mrSges(6,1) * t161 - mrSges(6,2) * t162) + (-Ifges(5,1) * t301 - Ifges(5,4) * t302) * t430 + m(7) * (t104 * t20 + t2 * t34 + t22 * t7 + t25 * t8 + t3 * t39 + t47 * t76) + m(6) * (t10 * t32 + t107 * t128 + t147 * t63 + t33 * t9 + t5 * t58 + t57 * t6) + (-t100 * t302 - t29 * t376 - t30 * t375 + t301 * t99) * mrSges(5,3) + t171 * (mrSges(5,1) * t302 - mrSges(5,2) * t301) + (Ifges(3,4) * t279 + Ifges(3,2) * t280) * t341 + (Ifges(3,1) * t279 + Ifges(3,4) * t280) * t342 - (t334 + t499) * t498 / 0.2e1 - (t249 * mrSges(4,1) - t115 * mrSges(4,3) - Ifges(4,4) * t176 + Ifges(5,5) * t448 - Ifges(4,2) * t177 - Ifges(4,6) * qJDD(3) + Ifges(5,6) * t447 + Ifges(5,3) * t435 + t436 * t505 + t452 * t481 + t453 * t483 + t487 + t489) * t498 + (t481 * t75 + t483 * t74) * t428 + (-t161 * t481 - t162 * t483) * t436 + (t482 * t75 + t484 * t74) * t444 + (-t161 * t482 - t162 * t484) * t452 + (t484 * t75 + t485 * t74) * t441 + (-t161 * t484 - t162 * t485) * t453 + m(5) * (t100 * t61 + t124 * t30 + t125 * t29 + t62 * t99) + m(4) * (t115 * t185 + t145 * t179 - t249 * t265) + (-t221 * mrSges(5,1) - t220 * mrSges(5,2) - t520 * t205 - t519 * t204 + (t403 * t493 + t490 + t495) * t289 + (-m(6) * (-t265 - t304) - m(5) * (-t265 - t318) - m(7) * (-t265 - t305) + m(4) * t265 + t488) * t285) * g(1) + (-t364 * t455 - m(7) * t374 - t223 * mrSges(5,1) - t222 * mrSges(5,2) + t493 * (t289 * t265 + t285 * t403) - t520 * t207 - t519 * t206 + t495 * t285 + (-m(5) * t318 - m(6) * t304 - m(7) * t305 - t488) * t289) * g(2) + (-t161 * t3 + t162 * t2 - t22 * t74 + t25 * t75) * mrSges(7,3) + 0.2e1 * t352 * t261 * mrSges(3,3) + (t494 / 0.2e1 + t483 * t441 - Ifges(4,4) * t423 - Ifges(4,2) * t425 + Ifges(5,5) * t430 - t153 / 0.2e1 + t481 * t444 + t501 / 0.2e1 + t502 / 0.2e1 - t179 * mrSges(4,3) + t505 * t428 - t514) * t231 - t503 * t162 / 0.2e1 - (-Ifges(4,4) * t425 - Ifges(4,1) * t423 - t329 + t365 / 0.2e1 - t154 / 0.2e1 + t515) * t230 + (-m(4) * t178 - t516) * t146 + Ifges(2,3) * qJDD(1) + m(3) * (-pkin(1) * t269 + (t261 + t344) * qJ(2) * t352) - t59 * t376 / 0.2e1 + t219 * (-Ifges(5,5) * t301 - Ifges(5,6) * t302) / 0.2e1 + t145 * t199 + t185 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t177) + t147 * t16 + t61 * t142 + t62 * t143 - t120 * t330 / 0.2e1 + t124 * t87 + t125 * t88 + t128 * (-mrSges(6,1) * t75 + mrSges(6,2) * t74) + t189 * (-Ifges(5,4) * t301 - Ifges(5,2) * t302) / 0.2e1 + t269 * t316 - pkin(1) * t317 + t375 * t451 + t161 * t521 + t34 * t35 + t39 * t37 + t57 * t36 + t58 * t38 + t76 * (-mrSges(7,1) * t75 + mrSges(7,2) * t74) + t47 * t80 + t104 * t15 + t107 * t81 + t8 * t111 + t9 * t112 + t7 * t113 + t10 * t114 - t265 * t324; t317 - (t35 + t36) * t303 + (-t80 - t81 + t466) * t229 + t324 + t478 * t242 + m(3) * t269 + t287 * t87 + t283 * t88 + (-t199 - t306) * t228 + t306 * qJD(4) + (t112 + t111) * t468 + (t113 + t114) * t467 + (-g(1) * t285 + g(2) * t289) * (m(3) - t493) - t328 * t352 * qJD(1) ^ 2 + (-t2 * t303 + t22 * t467 - t229 * t76 + t242 * t3 + t25 * t468) * m(7) + (-t128 * t229 + t242 * t5 - t303 * t6 + t32 * t467 + t33 * t468) * m(6) + (-t171 * t229 + t283 * t29 + t287 * t30 + t219 * (t100 * t287 - t283 * t99)) * m(5) + (t178 * t229 - t179 * t228 + t249) * m(4); (t300 + t329) * qJD(4) - (Ifges(4,1) * t228 - t387 + t494) * t229 / 0.2e1 - (t428 * t481 + t441 * t484 + t444 * t482 + t496) * t182 + t148 * t509 + t149 * t510 - t76 * (mrSges(7,1) * t148 - mrSges(7,2) * t149) + (t189 * t309 + t190 * t310 + t219 * t308) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t229 + t154 + t215 + t361) * t228 / 0.2e1 + (t512 * g(3) + t497 * (mrSges(4,1) - t513)) * t270 + (t513 * g(3) + t497 * (mrSges(4,2) + t512)) * t271 + (Ifges(5,5) * t431 + Ifges(5,6) * t432 + Ifges(5,3) * t427 + t505 * t429 + t483 * t442 + t481 * t445 + t514) * t229 + (-pkin(3) * t110 - t100 * t118 - t117 * t99) * m(5) + t503 * t242 / 0.2e1 + t120 * t325 + (t308 * t427 + t309 * t432 + t310 * t431 - t300 + t515) * t228 + (t388 + t516) * t179 + (t242 * t485 - t303 * t484) * t453 + (t148 * t25 - t149 * t22 - t2 * t242 - t3 * t303) * mrSges(7,3) + (t148 * t33 - t149 * t32 - t242 * t6 - t303 * t5) * mrSges(6,3) + t20 * (mrSges(7,1) * t303 + mrSges(7,2) * t242) + t63 * (mrSges(6,1) * t303 + mrSges(6,2) * t242) + (t242 * t483 - t303 * t481) * t436 + (t242 * t484 - t303 * t482) * t452 - t349 * t393 + t287 * t59 / 0.2e1 - t268 * t16 + Ifges(4,3) * qJDD(3) + t210 * t15 - t178 * t199 + t187 * t36 + t188 * t38 + Ifges(4,5) * t176 + Ifges(4,6) * t177 + t151 * t35 + t152 * t37 - t118 * t142 - t117 * t143 + (t287 * t88 + m(5) * ((-t100 * t283 - t99 * t287) * qJD(4) + t461) - t143 * t348 - t142 * t349 - t283 * t87) * pkin(8) + (t100 * t380 + (-t348 + t379) * t99 + t461) * mrSges(5,3) + t464 * t81 + t110 * t313 + (-mrSges(6,1) * t467 + mrSges(6,2) * t468) * t128 - g(3) * t315 + t473 * t80 + t474 * t114 + t475 * t112 + (t128 * t464 + t187 * t6 + t188 * t5 - t268 * t63 + t32 * t474 + t33 * t475) * m(6) + t153 * t423 + t365 * t425 + (Ifges(5,5) * t283 + Ifges(5,6) * t287) * t435 + (Ifges(5,2) * t287 + t395) * t447 + (Ifges(5,1) * t283 + t394) * t448 + t283 * t451 + t303 * t521 + t479 * t111 + t480 * t113 + (t151 * t2 + t152 * t3 + t20 * t210 + t22 * t480 + t25 * t479 + t473 * t76) * m(7) + (-t148 * t481 - t149 * t483) * t429 + (-t148 * t482 - t149 * t484) * t445 + (-t148 * t484 - t149 * t485) * t442 - (t483 * t428 + t485 * t441 + t484 * t444 + t456 + t510) * t181 - pkin(3) * t64 - t115 * mrSges(4,2) + t116 * mrSges(4,1); (-Ifges(5,2) * t190 + t121 + t186) * t432 + (-t128 * mrSges(6,1) - t481 * t429 - t484 * t442 - t482 * t445 + t496) * t136 + (mrSges(7,1) * t272 + t312 + t460 - t490) * t407 + t295 + t489 + t478 * pkin(4) * t282 + t334 + t476 * t445 + (-t142 + t400) * t99 - t81 * t415 - m(6) * (t128 * t415 + t32 * t40 + t33 * t41) + (-t337 - t40) * t114 + t267 * t35 + (-t337 - t26) * t113 + (t336 - t41) * t112 - t171 * (mrSges(5,1) * t190 + mrSges(5,2) * t189) + (t336 - t27) * t111 + t100 * t143 + (-m(7) * (-t252 * t289 - t271 * t374) - mrSges(5,2) * t221 + t472 * t220 + t462) * g(2) + (-m(7) * (-t251 * t271 * t289 + t252 * t285) + mrSges(5,2) * t223 - t472 * t222 + t463) * g(1) + t190 * t393 + t36 * t412 + (Ifges(5,5) * t189 - Ifges(5,6) * t190) * t427 + t120 * t430 + (Ifges(5,1) * t189 - t389) * t431 + (t282 * t5 + t286 * t6 + (-t282 * t32 + t286 * t33) * qJD(5)) * t455 + (t2 * t267 + (t282 * t3 + (-t22 * t282 + t25 * t286) * qJD(5)) * pkin(4) - t105 * t76 - t22 * t26 - t25 * t27) * m(7) + (-t128 * mrSges(6,2) + t429 * t483 + t442 * t485 + t445 * t484 - t456) * t320 - t105 * t80; t295 - t128 * (mrSges(6,1) * t136 + mrSges(6,2) * t320) - t76 * (mrSges(7,1) * t136 + mrSges(7,2) * t320) + t22 * t397 + t2 * t454 - t24 * t111 + (t320 * t485 - t517) * t442 + t477 * t441 + (-t136 * t481 + t320 * t483) * t429 + (-(-mrSges(7,1) - t454) * t272 + t460) * t407 + (t398 + t114) * t33 + (t399 - t112) * t32 + (-m(7) * (-t22 + t24) + t396 + t113) * t25 + (t204 * t454 + t462) * g(2) + (-t206 * t454 + t463) * g(1) + (-t136 * t482 + t476 + t518) * t445 + (t35 + (-m(7) * t76 - t80) * t136) * pkin(5); -t320 * t111 + t136 * t113 + (g(3) * t271 + t136 * t22 - t320 * t25 - t270 * t497 + t20) * m(7) + t15;];
tau  = t1;
