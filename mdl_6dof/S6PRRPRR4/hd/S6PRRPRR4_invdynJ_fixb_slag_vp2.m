% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:20
% EndTime: 2019-03-08 22:10:57
% DurationCPUTime: 21.65s
% Computational Cost: add. (7341->682), mult. (16127->914), div. (0->0), fcn. (11818->12), ass. (0->301)
t479 = mrSges(5,1) + mrSges(4,1);
t481 = -m(7) - m(6);
t477 = Ifges(4,5) + Ifges(5,4);
t476 = Ifges(5,6) - Ifges(4,6);
t263 = sin(qJ(3));
t428 = pkin(8) - pkin(9);
t216 = t428 * t263;
t267 = cos(qJ(3));
t217 = t428 * t267;
t262 = sin(qJ(5));
t266 = cos(qJ(5));
t119 = t216 * t262 + t217 * t266;
t361 = qJD(3) * t263;
t193 = t428 * t361;
t268 = cos(qJ(2));
t257 = sin(pkin(6));
t367 = qJD(1) * t257;
t226 = t268 * t367;
t310 = t267 * t226;
t311 = t263 * t226;
t312 = qJD(3) * t217;
t469 = -qJD(5) * t119 + (-t311 + t312) * t266 + (t193 + t310) * t262;
t181 = t262 * t263 + t266 * t267;
t374 = t257 * t268;
t123 = t181 * t374;
t288 = t266 * t216 - t217 * t262;
t494 = qJD(1) * t123 - qJD(5) * t288 + t266 * t193 - t262 * t312;
t357 = qJD(5) * t266;
t358 = qJD(5) * t262;
t359 = qJD(3) * t267;
t360 = qJD(3) * t266;
t100 = t262 * t359 - t267 * t358 + (t357 - t360) * t263;
t269 = -pkin(3) - pkin(4);
t336 = t269 * qJD(3);
t246 = t263 * qJD(4);
t368 = qJ(4) * t359 + t246;
t129 = t263 * t336 + t368;
t264 = sin(qJ(2));
t335 = t264 * t367;
t277 = t181 * qJD(5);
t99 = qJD(3) * t181 - t277;
t493 = pkin(5) * t100 - pkin(10) * t99 + t129 + t335;
t364 = qJD(2) * t263;
t343 = mrSges(4,3) * t364;
t345 = mrSges(5,2) * t364;
t492 = qJD(3) * t479 - t343 - t345;
t362 = qJD(2) * t267;
t344 = mrSges(5,2) * t362;
t208 = qJD(3) * mrSges(5,3) + t344;
t342 = mrSges(4,3) * t362;
t459 = -qJD(3) * mrSges(4,2) + t208 + t342;
t197 = qJD(2) * pkin(8) + t335;
t259 = cos(pkin(6));
t366 = qJD(1) * t259;
t124 = -t263 * t197 + t267 * t366;
t107 = pkin(9) * t364 + t124;
t491 = qJD(4) - t107;
t383 = sin(pkin(11));
t322 = t383 * t264;
t258 = cos(pkin(11));
t372 = t258 * t268;
t161 = -t259 * t322 + t372;
t323 = t257 * t383;
t103 = t161 * t263 - t267 * t323;
t104 = t161 * t267 + t263 * t323;
t470 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t261 = sin(qJ(6));
t265 = cos(qJ(6));
t209 = -mrSges(7,1) * t265 + mrSges(7,2) * t261;
t471 = m(7) * pkin(5) + mrSges(6,1) - t209;
t52 = t103 * t262 + t104 * t266;
t490 = t470 * t52 - t471 * (t103 * t266 - t104 * t262);
t321 = t383 * t268;
t373 = t258 * t264;
t159 = t259 * t373 + t321;
t375 = t257 * t267;
t101 = t159 * t263 + t258 * t375;
t102 = -t258 * t257 * t263 + t159 * t267;
t48 = t101 * t262 + t102 * t266;
t489 = t470 * t48 - t471 * (t101 * t266 - t102 * t262);
t376 = t257 * t264;
t163 = -t259 * t267 + t263 * t376;
t164 = t259 * t263 + t264 * t375;
t290 = t163 * t266 - t164 * t262;
t88 = t163 * t262 + t164 * t266;
t488 = -t471 * t290 + t470 * t88;
t302 = t261 * mrSges(7,1) + t265 * mrSges(7,2);
t89 = t336 + t491;
t255 = qJD(3) * qJ(4);
t125 = t267 * t197 + t263 * t366;
t286 = -pkin(9) * t362 + t125;
t93 = t255 + t286;
t39 = -t262 * t93 + t266 * t89;
t449 = -qJD(5) + qJD(3);
t33 = pkin(5) * t449 - t39;
t487 = t302 * t33;
t172 = -t262 * t362 + t266 * t364;
t116 = -t172 * t261 - t265 * t449;
t117 = t172 * t265 - t261 * t449;
t409 = mrSges(6,3) * t172;
t385 = -mrSges(6,1) * t449 + mrSges(7,1) * t116 - mrSges(7,2) * t117 - t409;
t365 = qJD(2) * t257;
t331 = qJD(1) * t365;
t219 = t268 * t331;
t352 = qJDD(1) * t257;
t147 = t264 * t352 + t219;
t485 = qJDD(2) * pkin(8) + qJD(3) * t366 + t147;
t305 = mrSges(4,1) * t267 - mrSges(4,2) * t263;
t304 = mrSges(5,1) * t267 + mrSges(5,3) * t263;
t184 = t304 * qJD(2);
t170 = -t262 * t364 - t266 * t362;
t90 = -mrSges(6,1) * t170 + mrSges(6,2) * t172;
t384 = -t184 - t90;
t484 = -t305 * qJD(2) + t384;
t287 = t262 * t267 - t263 * t266;
t447 = -mrSges(3,1) - t304 - t305;
t483 = -t181 * t471 + t287 * t470 + t447;
t198 = -qJD(2) * pkin(2) - t226;
t128 = -pkin(3) * t362 - qJ(4) * t364 + t198;
t108 = pkin(4) * t362 - t128;
t40 = t262 * t89 + t266 * t93;
t34 = -pkin(10) * t449 + t40;
t57 = -pkin(5) * t170 - pkin(10) * t172 + t108;
t15 = -t261 * t34 + t265 * t57;
t154 = Ifges(6,4) * t170;
t155 = Ifges(6,4) * t172;
t16 = t261 * t57 + t265 * t34;
t251 = -qJDD(3) + qJDD(5);
t294 = Ifges(7,5) * t265 - Ifges(7,6) * t261;
t404 = Ifges(7,4) * t265;
t297 = -Ifges(7,2) * t261 + t404;
t405 = Ifges(7,4) * t261;
t300 = Ifges(7,1) * t265 - t405;
t356 = qJD(6) * t261;
t326 = -t356 / 0.2e1;
t416 = t265 / 0.2e1;
t109 = Ifges(7,4) * t116;
t157 = qJD(6) - t170;
t43 = t117 * Ifges(7,1) + t157 * Ifges(7,5) + t109;
t337 = t43 * t416;
t355 = qJD(6) * t265;
t379 = t170 * t265;
t380 = t170 * t261;
t354 = qJD(2) * qJD(3);
t195 = qJDD(2) * t263 + t267 * t354;
t351 = qJDD(1) * t259;
t61 = -t197 * t359 - t263 * t485 + t267 * t351;
t274 = qJDD(4) - t61;
t35 = -pkin(9) * t195 + qJDD(3) * t269 + t274;
t194 = -t267 * qJDD(2) + t263 * t354;
t60 = -t197 * t361 + t263 * t351 + t267 * t485;
t53 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t60;
t36 = pkin(9) * t194 + t53;
t6 = -qJD(5) * t40 - t262 * t36 + t266 * t35;
t4 = -pkin(5) * t251 - t6;
t41 = t117 * Ifges(7,5) + t116 * Ifges(7,6) + t157 * Ifges(7,3);
t410 = mrSges(6,3) * t170;
t419 = t172 / 0.2e1;
t397 = t117 * Ifges(7,4);
t42 = t116 * Ifges(7,2) + t157 * Ifges(7,6) + t397;
t422 = t157 / 0.2e1;
t423 = -t157 / 0.2e1;
t424 = t117 / 0.2e1;
t425 = -t117 / 0.2e1;
t426 = t116 / 0.2e1;
t427 = -t116 / 0.2e1;
t72 = t287 * qJD(2) * qJD(5) + t194 * t266 - t195 * t262;
t70 = qJDD(6) - t72;
t429 = t70 / 0.2e1;
t71 = -qJD(2) * t277 + t194 * t262 + t195 * t266;
t29 = -qJD(6) * t117 + t251 * t265 - t261 * t71;
t431 = t29 / 0.2e1;
t28 = qJD(6) * t116 + t251 * t261 + t265 * t71;
t432 = t28 / 0.2e1;
t433 = Ifges(7,1) * t432 + Ifges(7,4) * t431 + Ifges(7,5) * t429;
t218 = t264 * t331;
t146 = t268 * t352 - t218;
t132 = -qJDD(2) * pkin(2) - t146;
t65 = t194 * pkin(3) - t195 * qJ(4) - qJD(2) * t246 + t132;
t44 = -pkin(4) * t194 - t65;
t14 = -pkin(5) * t72 - pkin(10) * t71 + t44;
t5 = t262 * t35 + t266 * t36 + t89 * t357 - t358 * t93;
t3 = pkin(10) * t251 + t5;
t1 = qJD(6) * t15 + t14 * t261 + t265 * t3;
t2 = -qJD(6) * t16 + t14 * t265 - t261 * t3;
t450 = t1 * t265 - t2 * t261;
t7 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + t70 * Ifges(7,6);
t472 = t449 * Ifges(6,6);
t474 = t170 * Ifges(6,2);
t83 = t155 - t472 + t474;
t473 = t449 * Ifges(6,5);
t84 = t172 * Ifges(6,1) + t154 - t473;
t482 = (Ifges(7,5) * t261 + Ifges(7,6) * t265) * t429 + (Ifges(7,2) * t265 + t405) * t431 + (Ifges(7,1) * t261 + t404) * t432 + (-t15 * t355 - t16 * t356 + t450) * mrSges(7,3) + (Ifges(7,5) * t172 + t170 * t300) * t425 + (Ifges(7,6) * t172 + t170 * t297) * t427 + (Ifges(7,3) * t172 + t170 * t294) * t423 + t449 * (Ifges(6,5) * t170 - Ifges(6,6) * t172) / 0.2e1 - (Ifges(6,1) * t170 - t155 + t41) * t172 / 0.2e1 - (-Ifges(6,2) * t172 + t154 + t84) * t170 / 0.2e1 + t261 * t433 + t7 * t416 + t83 * t419 + t16 * (mrSges(7,2) * t172 + mrSges(7,3) * t380) + t15 * (-mrSges(7,1) * t172 + mrSges(7,3) * t379) - t108 * (mrSges(6,1) * t172 + mrSges(6,2) * t170) + (t380 / 0.2e1 + t326) * t42 - t170 * t487 + (t294 * t422 + t297 * t426 + t300 * t424 + t337 + t487) * qJD(6) - t43 * t379 / 0.2e1 + Ifges(6,3) * t251 + t4 * t209 + Ifges(6,5) * t71 + Ifges(6,6) * t72 + t40 * t409 + t39 * t410 - t5 * mrSges(6,2) + t6 * mrSges(6,1);
t247 = t263 * qJ(4);
t201 = -t267 * pkin(3) - pkin(2) - t247;
t176 = t267 * pkin(4) - t201;
t79 = pkin(5) * t181 + pkin(10) * t287 + t176;
t37 = -t119 * t261 + t265 * t79;
t480 = qJD(6) * t37 + t261 * t493 - t265 * t494;
t478 = -mrSges(5,3) + mrSges(4,2);
t11 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t63 = mrSges(6,1) * t251 - mrSges(6,3) * t71;
t411 = t11 - t63;
t38 = t119 * t265 + t261 * t79;
t475 = -qJD(6) * t38 + t261 * t494 + t265 * t493;
t158 = t259 * t372 - t322;
t382 = t158 * t267;
t466 = pkin(3) * t382 + t158 * t247;
t160 = -t259 * t321 - t373;
t381 = t160 * t267;
t465 = pkin(3) * t381 + t160 * t247;
t151 = -mrSges(5,2) * t194 + qJDD(3) * mrSges(5,3);
t464 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t194 + t151;
t150 = -qJDD(3) * mrSges(5,1) + t195 * mrSges(5,2);
t463 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t195 + t150;
t241 = Ifges(4,4) * t362;
t402 = Ifges(5,5) * t267;
t301 = Ifges(5,1) * t263 - t402;
t461 = Ifges(4,1) * t364 + qJD(2) * t301 + qJD(3) * t477 + t241;
t200 = t266 * qJ(4) + t262 * t269;
t458 = t198 * (mrSges(4,1) * t263 + mrSges(4,2) * t267) + t128 * (mrSges(5,1) * t263 - mrSges(5,3) * t267);
t457 = t263 * t476 + t267 * t477;
t283 = t261 * t99 - t287 * t355;
t454 = -t263 * t61 + t267 * t60;
t56 = -qJDD(3) * pkin(3) + t274;
t453 = t263 * t56 + t267 * t53;
t17 = mrSges(7,1) * t70 - mrSges(7,3) * t28;
t18 = -mrSges(7,2) * t70 + mrSges(7,3) * t29;
t452 = -t261 * t17 + t265 * t18;
t448 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t126 = mrSges(6,2) * t449 + t410;
t73 = -mrSges(7,2) * t157 + mrSges(7,3) * t116;
t74 = mrSges(7,1) * t157 - mrSges(7,3) * t117;
t446 = -t261 * t74 + t265 * t73 + t126;
t439 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t272 = (-t15 * t265 - t16 * t261) * qJD(6) + t450;
t438 = m(7) * t272 - t74 * t355 - t73 * t356 + t452;
t435 = t428 * t481 + mrSges(6,3) + t302 + t448;
t270 = qJD(2) ^ 2;
t417 = t263 / 0.2e1;
t407 = Ifges(4,4) * t263;
t406 = Ifges(4,4) * t267;
t403 = Ifges(5,5) * t263;
t378 = t287 * t261;
t377 = t287 * t265;
t369 = pkin(2) * t374 + pkin(8) * t376;
t363 = qJD(2) * t264;
t349 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t70;
t339 = t263 * t374;
t338 = t267 * t374;
t334 = t257 * t363;
t333 = t268 * t365;
t144 = t158 * pkin(2);
t325 = pkin(8) * t159 + t144;
t145 = t160 * pkin(2);
t324 = pkin(8) * t161 + t145;
t320 = -t354 / 0.2e1;
t318 = -t101 * pkin(3) + qJ(4) * t102;
t317 = -t103 * pkin(3) + qJ(4) * t104;
t316 = -t163 * pkin(3) + qJ(4) * t164;
t313 = pkin(3) * t338 + qJ(4) * t339 + t369;
t91 = pkin(5) * t172 - pkin(10) * t170;
t299 = t267 * Ifges(4,2) + t407;
t293 = -t15 * t261 + t16 * t265;
t199 = -qJ(4) * t262 + t266 * t269;
t238 = qJ(4) * t362;
t143 = t269 * t364 + t238;
t66 = -t261 * t88 + t265 * t374;
t67 = t261 * t374 + t265 * t88;
t282 = -t265 * t99 - t287 * t356;
t279 = t263 * (Ifges(4,1) * t267 - t407);
t278 = t267 * (Ifges(5,3) * t263 + t402);
t275 = -g(1) * t103 - g(2) * t101 - g(3) * t163;
t273 = pkin(4) * t338 - pkin(9) * t376 + t313;
t240 = Ifges(5,5) * t364;
t191 = pkin(5) - t199;
t186 = pkin(3) * t364 - t238;
t171 = t261 * t364 + t265 * t360;
t169 = -t261 * t360 + t265 * t364;
t166 = Ifges(4,6) * qJD(3) + qJD(2) * t299;
t165 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t362 + t240;
t156 = pkin(3) * t361 - t368;
t115 = t255 + t125;
t112 = -qJD(3) * pkin(3) + qJD(4) - t124;
t111 = mrSges(4,1) * t194 + mrSges(4,2) * t195;
t110 = mrSges(5,1) * t194 - mrSges(5,3) * t195;
t106 = -qJD(3) * t163 + t267 * t333;
t105 = qJD(3) * t164 + t263 * t333;
t68 = t143 - t91;
t64 = -mrSges(6,2) * t251 + mrSges(6,3) * t72;
t55 = t266 * t107 + t262 * t286;
t25 = -mrSges(6,1) * t72 + mrSges(6,2) * t71;
t24 = qJD(5) * t290 + t105 * t262 + t106 * t266;
t23 = qJD(5) * t88 - t105 * t266 + t106 * t262;
t22 = t261 * t91 + t265 * t39;
t21 = -t261 * t39 + t265 * t91;
t20 = t261 * t68 + t265 * t55;
t19 = -t261 * t55 + t265 * t68;
t13 = qJD(6) * t66 + t24 * t265 - t261 * t334;
t12 = -qJD(6) * t67 - t24 * t261 - t265 * t334;
t8 = [m(2) * qJDD(1) + t12 * t74 + t24 * t126 + t13 * t73 + t66 * t17 + t67 * t18 + t88 * t64 - t411 * t290 - t385 * t23 + t464 * t164 + t463 * t163 + t459 * t106 - t492 * t105 + (-m(2) - m(3) - m(4) - m(5) + t481) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t270 - t110 - t111 + t25) * t268 + (-mrSges(3,1) * t270 - mrSges(3,2) * qJDD(2) + qJD(2) * t484) * t264) * t257 + m(4) * (-t105 * t124 + t106 * t125 - t163 * t61 + t164 * t60 + (-t132 * t268 + t198 * t363) * t257) + m(5) * (t105 * t112 + t106 * t115 + t163 * t56 + t164 * t53 + (t128 * t363 - t268 * t65) * t257) + m(6) * (-t23 * t39 + t24 * t40 + t5 * t88 + t6 * t290 + (-t108 * t363 + t268 * t44) * t257) + m(3) * (qJDD(1) * t259 ^ 2 + (t146 * t268 + t147 * t264) * t257) + m(7) * (t1 * t67 + t12 * t15 + t13 * t16 + t2 * t66 + t23 * t33 - t290 * t4); (t403 / 0.2e1 - t299 / 0.2e1 + (Ifges(5,5) - Ifges(4,4)) * t417 + (-Ifges(5,3) - Ifges(4,2) / 0.2e1) * t267) * t194 + t267 * (Ifges(4,4) * t195 + Ifges(4,6) * qJDD(3)) / 0.2e1 + (t1 * t378 + t15 * t282 - t16 * t283 + t2 * t377) * mrSges(7,3) + (t279 + t263 * (Ifges(5,1) * t267 + t403) + t267 * (-Ifges(4,2) * t263 + t406)) * t354 / 0.2e1 + (Ifges(4,1) * t263 + t301 + t406) * t195 / 0.2e1 + t475 * t74 + (t477 * t263 - t476 * t267) * qJDD(3) / 0.2e1 + ((Ifges(4,1) + Ifges(5,1)) * t195 + t477 * qJDD(3)) * t417 + (t108 * t129 + t119 * t5 + t176 * t44 + t288 * t6 + t39 * t469 - t40 * t494) * m(6) - t494 * t126 + (((t112 * t267 - t115 * t263) * qJD(3) + t453) * m(5) + ((-t124 * t267 - t125 * t263) * qJD(3) + t454) * m(4) + t463 * t263 + t464 * t267) * pkin(8) + t461 * t359 / 0.2e1 + (t112 * t359 - t115 * t361 + t453) * mrSges(5,2) + (-t124 * t359 - t125 * t361 + t454) * mrSges(4,3) - t283 * t42 / 0.2e1 + (t458 + t457 * qJD(3) / 0.2e1) * qJD(3) + (m(6) * t108 - t484) * t335 + t469 * t385 + (-m(5) * (t325 + t466) - m(4) * t325 + t481 * (pkin(4) * t382 + t144 + t466) + t435 * t159 + t483 * t158) * g(2) + (-m(5) * (t324 + t465) - m(4) * t324 + t481 * (pkin(4) * t381 + t145 + t465) + t435 * t161 + t483 * t160) * g(1) + (t349 / 0.2e1 + Ifges(7,3) * t429 + Ifges(7,6) * t431 + Ifges(7,5) * t432 - Ifges(6,6) * t251 - Ifges(6,4) * t71 - Ifges(6,2) * t72 + t44 * mrSges(6,1) - t5 * mrSges(6,3) + t439) * t181 + t459 * (-pkin(8) * t361 - t310) + t492 * (-pkin(8) * t359 + t311) + (-Ifges(7,4) * t282 - Ifges(7,2) * t283) * t426 + (-(t128 * t264 + (t112 * t263 + t115 * t267) * t268) * t367 + t128 * t156 + t201 * t65) * m(5) + (-(t198 * t264 + (-t124 * t263 + t125 * t267) * t268) * t367 - pkin(2) * t132) * m(4) + (-t166 / 0.2e1 + t165 / 0.2e1) * t361 + t480 * t73 + t278 * t320 - t377 * t433 - t411 * t288 + (t1 * t38 + t475 * t15 + t480 * t16 + t2 * t37 - t288 * t4 - t469 * t33) * m(7) - (t44 * mrSges(6,2) - t6 * mrSges(6,3) + Ifges(6,1) * t71 + Ifges(6,4) * t72 + Ifges(6,5) * t251 + t294 * t429 + t297 * t431 + t300 * t432 + t302 * t4 + t326 * t43) * t287 + (-Ifges(6,4) * t419 + Ifges(7,3) * t422 + Ifges(7,5) * t424 + Ifges(7,6) * t426 + t472 / 0.2e1 - t474 / 0.2e1 + t108 * mrSges(6,1) - t83 / 0.2e1 + t41 / 0.2e1 + t15 * mrSges(7,1) - t16 * mrSges(7,2)) * t100 + (Ifges(6,1) * t419 - t473 / 0.2e1 + t154 / 0.2e1 + t108 * mrSges(6,2) + t84 / 0.2e1 + t337) * t99 + (t218 + t146) * mrSges(3,1) + (-Ifges(7,1) * t282 - Ifges(7,4) * t283) * t424 + (-m(4) * t369 - m(7) * (pkin(5) * t123 + t273) - (t123 * t265 - t261 * t376) * mrSges(7,1) - (-t123 * t261 - t265 * t376) * mrSges(7,2) - m(6) * t273 - t123 * mrSges(6,1) + mrSges(6,3) * t376 - m(5) * t313 + t470 * (t262 * t338 - t266 * t339) + (t264 * t448 + t268 * t447) * t257) * g(3) - t267 * (Ifges(5,5) * t195 + Ifges(5,6) * qJDD(3)) / 0.2e1 + t7 * t378 / 0.2e1 + (-Ifges(7,5) * t282 - Ifges(7,6) * t283) * t422 + Ifges(3,3) * qJDD(2) + t201 * t110 - t156 * t184 + t176 * t25 + t129 * t90 + t119 * t64 - pkin(2) * t111 + (-t100 * t40 - t39 * t99) * mrSges(6,3) + t37 * t17 + t38 * t18 + t33 * (mrSges(7,1) * t283 - mrSges(7,2) * t282) - t65 * t304 - t132 * t305 + (t219 - t147) * mrSges(3,2); (t278 / 0.2e1 - t279 / 0.2e1) * t270 - (Ifges(5,1) * t362 + t165 + t240) * t364 / 0.2e1 + t476 * t194 + t477 * t195 + (Ifges(4,3) + Ifges(5,2)) * qJDD(3) + (-m(5) * t112 + t343 + t492) * t125 - (-Ifges(4,2) * t364 + t241 + t461) * t362 / 0.2e1 + t457 * t320 - t458 * qJD(2) + (-m(5) * t115 + t342 - t459) * t124 + (m(6) * t40 + m(7) * t293 + t446) * (qJD(4) * t266 + qJD(5) * t199) + (-m(5) * t316 + t481 * (-t163 * pkin(4) + t316) + t478 * t164 + t479 * t163 - t488) * g(3) + (-m(5) * t318 + t481 * (-t101 * pkin(4) + t318) + t478 * t102 + t479 * t101 - t489) * g(2) + (-m(5) * t317 + t481 * (-t103 * pkin(4) + t317) + t478 * t104 + t479 * t103 - t490) * g(1) - t482 + (-t39 * m(6) + t33 * m(7) - t385) * (qJD(5) * t200 + t262 * t491 + t266 * t286) - t112 * t344 + t166 * t364 / 0.2e1 + (-t15 * t19 - t16 * t20 + t191 * t4) * m(7) + (-t108 * t143 + t199 * t6 + t200 * t5 - t40 * t55) * m(6) + (-pkin(3) * t56 + qJ(4) * t53 + qJD(4) * t115 - t128 * t186) * m(5) + qJD(4) * t208 + t199 * t63 + t200 * t64 + t191 * t11 + t186 * t184 + t438 * (-pkin(10) + t200) - pkin(3) * t150 + qJ(4) * t151 - t143 * t90 - t55 * t126 - t20 * t73 - t19 * t74 - t60 * mrSges(4,2) + t61 * mrSges(4,1) + t53 * mrSges(5,3) - t56 * mrSges(5,1) + t115 * t345; -qJD(3) * t208 - t169 * t74 - t171 * t73 + t384 * t364 + (-qJD(3) * t126 + t446 * qJD(5) - t411) * t266 + (t64 + (-t261 * t73 - t265 * t74) * qJD(6) + t449 * t385 + t452) * t262 + t150 + (t275 + (qJD(5) * t293 - t4) * t266 - t15 * t169 - t16 * t171 + (-t33 * t449 + t272) * t262) * m(7) + (-t108 * t364 + t262 * t5 + t266 * t6 + t275 - t449 * (-t262 * t39 + t266 * t40)) * m(6) + (-qJD(3) * t115 + t128 * t364 + t275 + t56) * m(5); t488 * g(3) + t489 * g(2) + t490 * g(1) + t385 * t40 + t438 * pkin(10) - t39 * t126 - t22 * t73 - t21 * t74 + (-pkin(5) * t4 - t15 * t21 - t16 * t22 - t33 * t40) * m(7) - pkin(5) * t11 + t482; -t33 * (mrSges(7,1) * t117 + mrSges(7,2) * t116) + (Ifges(7,1) * t116 - t397) * t425 + t42 * t424 + (Ifges(7,5) * t116 - Ifges(7,6) * t117) * t423 - t15 * t73 + t16 * t74 - g(1) * ((t160 * t265 - t261 * t52) * mrSges(7,1) + (-t160 * t261 - t265 * t52) * mrSges(7,2)) - g(2) * ((t158 * t265 - t261 * t48) * mrSges(7,1) + (-t158 * t261 - t265 * t48) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t66 - mrSges(7,2) * t67) + (t116 * t15 + t117 * t16) * mrSges(7,3) + t349 + (-Ifges(7,2) * t117 + t109 + t43) * t427 + t439;];
tau  = t8;
