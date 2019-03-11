% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:20
% EndTime: 2019-03-09 06:27:04
% DurationCPUTime: 30.74s
% Computational Cost: add. (9254->701), mult. (18475->903), div. (0->0), fcn. (11891->10), ass. (0->319)
t514 = -Ifges(4,2) / 0.2e1;
t260 = sin(qJ(3));
t249 = t260 * qJD(1);
t240 = t249 + qJD(4);
t235 = qJD(5) + t240;
t406 = -t235 / 0.2e1;
t259 = sin(qJ(4));
t263 = cos(qJ(4));
t346 = qJD(3) * t263;
t264 = cos(qJ(3));
t348 = qJD(1) * t264;
t206 = -t259 * t348 + t346;
t207 = qJD(3) * t259 + t263 * t348;
t258 = sin(qJ(5));
t262 = cos(qJ(5));
t127 = t206 * t258 + t207 * t262;
t420 = -t127 / 0.2e1;
t467 = Ifges(6,5) + Ifges(7,5);
t486 = Ifges(6,3) + Ifges(7,3);
t513 = t486 * t406 + t467 * t420;
t512 = -t206 / 0.2e1;
t511 = -t207 / 0.2e1;
t510 = -t240 / 0.2e1;
t380 = Ifges(4,4) * t264;
t509 = t260 * t514 + t380 / 0.2e1;
t465 = Ifges(6,6) + Ifges(7,6);
t508 = -t465 / 0.2e1;
t468 = Ifges(6,4) + Ifges(7,4);
t405 = t235 / 0.2e1;
t419 = t127 / 0.2e1;
t307 = t262 * t206 - t207 * t258;
t422 = t307 / 0.2e1;
t267 = -pkin(1) - pkin(7);
t237 = qJD(1) * t267 + qJD(2);
t365 = t237 * t264;
t193 = -qJD(3) * pkin(3) - t365;
t137 = -pkin(4) * t206 + t193;
t305 = pkin(3) * t260 - pkin(8) * t264;
t226 = qJ(2) + t305;
t182 = t226 * qJD(1);
t223 = t260 * t237;
t192 = qJD(3) * pkin(8) + t223;
t110 = t263 * t182 - t192 * t259;
t90 = -pkin(9) * t207 + t110;
t80 = pkin(4) * t240 + t90;
t111 = t182 * t259 + t192 * t263;
t91 = pkin(9) * t206 + t111;
t83 = t258 * t91;
t32 = t262 * t80 - t83;
t483 = qJ(6) * t127;
t22 = t32 - t483;
t20 = pkin(5) * t235 + t22;
t73 = -pkin(5) * t307 + qJD(6) + t137;
t436 = -t137 * mrSges(6,2) - t73 * mrSges(7,2) + t32 * mrSges(6,3) + t20 * mrSges(7,3);
t469 = Ifges(6,1) + Ifges(7,1);
t495 = t468 * t307;
t458 = t127 * t469 + t467 * t235 + t495;
t507 = -t436 + t467 * t405 + t469 * t419 + t468 * t422 + t458 / 0.2e1;
t466 = Ifges(6,2) + Ifges(7,2);
t306 = pkin(3) * t264 + pkin(8) * t260;
t217 = t306 * qJD(1);
t364 = t259 * t264;
t133 = t263 * t217 - t237 * t364;
t266 = -pkin(9) - pkin(8);
t324 = qJD(4) * t266;
t361 = t260 * t263;
t334 = pkin(9) * t361;
t506 = t263 * t324 - (pkin(4) * t264 + t334) * qJD(1) - t133;
t356 = t263 * t264;
t134 = t259 * t217 + t237 * t356;
t323 = t259 * t249;
t505 = pkin(9) * t323 - t259 * t324 + t134;
t257 = qJ(4) + qJ(5);
t251 = cos(t257);
t397 = pkin(4) * t263;
t225 = pkin(5) * t251 + t397;
t245 = pkin(3) + t397;
t504 = -m(6) * t245 - m(7) * (pkin(3) + t225);
t336 = qJD(1) * qJD(3);
t221 = qJDD(1) * t264 - t260 * t336;
t222 = -t260 * qJDD(1) - t264 * t336;
t337 = qJD(1) * qJD(2);
t238 = qJDD(1) * qJ(2) + t337;
t122 = -pkin(3) * t222 - pkin(8) * t221 + t238;
t233 = qJDD(1) * t267 + qJDD(2);
t345 = qJD(3) * t264;
t143 = t260 * t233 + t237 * t345;
t139 = qJDD(3) * pkin(8) + t143;
t342 = qJD(4) * t263;
t343 = qJD(4) * t259;
t45 = t259 * t122 + t263 * t139 + t182 * t342 - t192 * t343;
t46 = -qJD(4) * t111 + t263 * t122 - t139 * t259;
t503 = t46 * mrSges(5,1) - t45 * mrSges(5,2);
t502 = -m(6) * t266 - m(7) * (-qJ(6) + t266) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t203 = qJDD(4) - t222;
t196 = qJDD(5) + t203;
t118 = qJD(4) * t206 + qJDD(3) * t259 + t221 * t263;
t119 = -qJD(4) * t207 + qJDD(3) * t263 - t221 * t259;
t40 = qJD(5) * t307 + t118 * t262 + t119 * t258;
t21 = pkin(4) * t203 - pkin(9) * t118 + t46;
t25 = pkin(9) * t119 + t45;
t85 = t262 * t91;
t33 = t258 * t80 + t85;
t6 = -qJD(5) * t33 + t262 * t21 - t25 * t258;
t2 = pkin(5) * t196 - qJ(6) * t40 - qJD(6) * t127 + t6;
t41 = -qJD(5) * t127 - t118 * t258 + t119 * t262;
t339 = qJD(5) * t262;
t340 = qJD(5) * t258;
t5 = t258 * t21 + t262 * t25 + t80 * t339 - t340 * t91;
t3 = qJ(6) * t41 + qJD(6) * t307 + t5;
t501 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t452 = qJ(6) * t307;
t23 = t33 + t452;
t493 = t468 * t127;
t459 = t235 * t465 + t307 * t466 + t493;
t475 = -t137 * mrSges(6,1) - t73 * mrSges(7,1) + mrSges(6,3) * t33 + mrSges(7,3) * t23 + t459 / 0.2e1;
t500 = t465 * t405 + t468 * t419 + t466 * t422 + t475;
t498 = t23 * mrSges(7,2) + t33 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t509 + Ifges(5,5) * t511 + Ifges(5,6) * t512 + Ifges(5,3) * t510 + t307 * t508 - t20 * mrSges(7,1) - t32 * mrSges(6,1) + t513;
t497 = t264 / 0.2e1;
t496 = -t466 * t41 / 0.2e1 - t468 * t40 / 0.2e1 + t196 * t508;
t488 = -mrSges(7,1) - mrSges(6,1);
t487 = mrSges(7,2) + mrSges(6,2);
t211 = t258 * t263 + t259 * t262;
t441 = qJD(4) + qJD(5);
t132 = t441 * t211;
t181 = t211 * qJD(1);
t157 = t260 * t181;
t494 = t132 + t157;
t261 = sin(qJ(1));
t265 = cos(qJ(1));
t442 = -g(1) * t261 + g(2) * t265;
t434 = m(6) * pkin(4);
t399 = pkin(4) * t259;
t489 = m(6) * t399;
t485 = t467 * t196 + t40 * t469 + t468 * t41;
t231 = t266 * t259;
t232 = t266 * t263;
t141 = t258 * t231 - t262 * t232;
t461 = -qJD(5) * t141 + t258 * t505 + t262 * t506;
t460 = t231 * t339 + t232 * t340 + t258 * t506 - t262 * t505;
t64 = -mrSges(5,1) * t119 + mrSges(5,2) * t118;
t484 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t221 - t64;
t287 = t258 * t259 - t262 * t263;
t482 = t441 * t287;
t481 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t206 + mrSges(5,2) * t207 + mrSges(4,3) * t348;
t347 = qJD(3) * t260;
t322 = t259 * t347;
t341 = qJD(4) * t264;
t480 = t263 * t341 - t322;
t448 = -t223 + (t323 + t343) * pkin(4);
t479 = t196 * t486 + t40 * t467 + t41 * t465;
t198 = Ifges(5,4) * t206;
t107 = t207 * Ifges(5,1) + t240 * Ifges(5,5) + t198;
t381 = Ifges(4,4) * t260;
t299 = Ifges(4,1) * t264 - t381;
t478 = Ifges(4,5) * qJD(3) + qJD(1) * t299 + t263 * t107;
t142 = t233 * t264 - t237 * t347;
t289 = t142 * t264 + t143 * t260;
t477 = m(5) + m(6) + m(7);
t304 = mrSges(4,1) * t264 - mrSges(4,2) * t260;
t476 = qJ(2) * t304 + (-Ifges(4,1) * t260 - t380) * t497;
t473 = -m(4) - t477;
t250 = sin(t257);
t224 = pkin(5) * t250 + t399;
t472 = m(7) * t224 + mrSges(2,1) - mrSges(3,2) + mrSges(4,3);
t303 = mrSges(4,1) * t260 + mrSges(4,2) * t264;
t471 = -m(5) * t305 + t260 * t504 + t264 * t502 + mrSges(2,2) - mrSges(3,3) - t303;
t426 = t118 / 0.2e1;
t425 = t119 / 0.2e1;
t411 = t203 / 0.2e1;
t423 = -t307 / 0.2e1;
t30 = -mrSges(7,2) * t196 + mrSges(7,3) * t41;
t31 = -mrSges(6,2) * t196 + mrSges(6,3) * t41;
t464 = t30 + t31;
t463 = -qJ(6) * t494 - qJD(6) * t287 + t460;
t451 = t287 * t260;
t158 = qJD(1) * t451;
t462 = -pkin(5) * t348 - qJD(6) * t211 + t461 + (t482 + t158) * qJ(6);
t383 = mrSges(7,3) * t307;
t95 = -mrSges(7,2) * t235 + t383;
t385 = mrSges(6,3) * t307;
t96 = -mrSges(6,2) * t235 + t385;
t457 = t95 + t96;
t382 = mrSges(7,3) * t127;
t97 = mrSges(7,1) * t235 - t382;
t384 = mrSges(6,3) * t127;
t98 = mrSges(6,1) * t235 - t384;
t456 = t97 + t98;
t455 = -t434 - mrSges(5,1);
t170 = t287 * t264;
t454 = -qJD(3) * t170 - t132 * t260 - t181;
t168 = t211 * t264;
t453 = t287 * qJD(1) - qJD(3) * t168 + t260 * t482;
t450 = pkin(5) * t494 + t448;
t202 = t263 * t226;
t311 = -t259 * t267 + pkin(4);
t123 = -pkin(9) * t356 + t260 * t311 + t202;
t359 = t260 * t267;
t234 = t263 * t359;
t148 = t259 * t226 + t234;
t136 = -pkin(9) * t364 + t148;
t68 = t258 * t123 + t262 * t136;
t360 = t260 * t265;
t164 = t250 * t360 + t251 * t261;
t165 = -t250 * t261 + t251 * t360;
t447 = t164 * t488 - t165 * t487;
t362 = t260 * t261;
t162 = -t250 * t362 + t251 * t265;
t163 = t250 * t265 + t251 * t362;
t446 = t162 * t488 + t163 * t487;
t86 = mrSges(5,1) * t203 - mrSges(5,3) * t118;
t87 = -mrSges(5,2) * t203 + mrSges(5,3) * t119;
t445 = -t259 * t86 + t263 * t87;
t444 = -t259 * t46 + t263 * t45;
t443 = mrSges(6,1) * t250 + t251 * t487;
t302 = -mrSges(5,1) * t263 + mrSges(5,2) * t259;
t440 = m(5) * pkin(3) - t250 * t487 - t251 * t488 - t302 - t504;
t439 = -m(5) * pkin(8) - t502;
t435 = qJD(1) ^ 2;
t433 = m(7) * pkin(5);
t432 = t40 / 0.2e1;
t431 = t41 / 0.2e1;
t430 = Ifges(5,1) * t426 + Ifges(5,4) * t425 + Ifges(5,5) * t411;
t427 = -m(3) - m(4);
t412 = t196 / 0.2e1;
t409 = t207 / 0.2e1;
t401 = pkin(4) * t207;
t398 = pkin(4) * t262;
t393 = g(3) * t264;
t44 = t262 * t90 - t83;
t379 = Ifges(5,4) * t259;
t378 = Ifges(5,4) * t263;
t377 = t110 * mrSges(5,3);
t376 = t111 * mrSges(5,3);
t373 = t207 * Ifges(5,4);
t363 = t259 * t265;
t358 = t261 * t263;
t355 = t263 * t265;
t354 = t264 * t267;
t349 = t265 * pkin(1) + t261 * qJ(2);
t344 = qJD(3) * t267;
t338 = qJDD(1) * mrSges(3,2);
t327 = t259 * t359;
t326 = Ifges(5,5) * t118 + Ifges(5,6) * t119 + Ifges(5,3) * t203;
t239 = t260 * t344;
t321 = t264 * t344;
t13 = -t41 * mrSges(7,1) + t40 * mrSges(7,2);
t313 = -t341 / 0.2e1;
t43 = -t258 * t90 - t85;
t310 = -t336 / 0.2e1;
t308 = (t238 + t337) * qJ(2);
t67 = t262 * t123 - t136 * t258;
t140 = t262 * t231 + t232 * t258;
t205 = pkin(4) * t364 - t354;
t301 = mrSges(5,1) * t259 + mrSges(5,2) * t263;
t298 = Ifges(5,1) * t263 - t379;
t297 = Ifges(5,1) * t259 + t378;
t295 = -Ifges(5,2) * t259 + t378;
t294 = Ifges(5,2) * t263 + t379;
t293 = -Ifges(4,5) * t260 - Ifges(4,6) * t264;
t292 = Ifges(5,5) * t263 - Ifges(5,6) * t259;
t291 = Ifges(5,5) * t259 + Ifges(5,6) * t263;
t290 = t110 * t263 + t111 * t259;
t145 = -mrSges(5,2) * t240 + mrSges(5,3) * t206;
t146 = mrSges(5,1) * t240 - mrSges(5,3) * t207;
t288 = -t259 * t145 - t263 * t146;
t185 = t259 * t360 + t358;
t183 = -t259 * t362 + t355;
t138 = -qJDD(3) * pkin(3) - t142;
t285 = t260 * (-Ifges(4,2) * t264 - t381);
t144 = pkin(4) * t480 + t239;
t204 = qJD(3) * t306 + qJD(2);
t174 = t263 * t204;
t65 = t174 + (-t234 + (pkin(9) * t264 - t226) * t259) * qJD(4) + (t264 * t311 + t334) * qJD(3);
t88 = -qJD(4) * t327 + t259 * t204 + t226 * t342 + t263 * t321;
t71 = -pkin(9) * t480 + t88;
t15 = t123 * t339 - t136 * t340 + t258 * t65 + t262 * t71;
t279 = t259 * t341 + t260 * t346;
t72 = -pkin(4) * t119 + t138;
t274 = Ifges(5,5) * t264 - t260 * t298;
t273 = Ifges(5,6) * t264 - t260 * t295;
t272 = Ifges(5,3) * t264 - t260 * t292;
t16 = -qJD(5) * t68 - t258 * t71 + t262 * t65;
t271 = t479 + t501;
t270 = -qJD(4) * t290 + t444;
t246 = -pkin(1) * qJDD(1) + qJDD(2);
t244 = pkin(5) + t398;
t228 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t249;
t215 = t303 * qJD(1);
t197 = t301 * t264;
t186 = -t259 * t261 + t260 * t355;
t184 = t260 * t358 + t363;
t172 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t222;
t167 = t211 * t260;
t160 = pkin(5) * t287 - t245;
t147 = t202 - t327;
t129 = pkin(5) * t168 + t205;
t106 = t206 * Ifges(5,2) + t240 * Ifges(5,6) + t373;
t103 = -qJ(6) * t287 + t141;
t102 = -qJ(6) * t211 + t140;
t92 = pkin(5) * t127 + t401;
t89 = -qJD(4) * t148 - t259 * t321 + t174;
t79 = t211 * t347 + t264 * t482;
t77 = qJD(3) * t451 - t132 * t264;
t70 = -mrSges(6,1) * t307 + mrSges(6,2) * t127;
t69 = -mrSges(7,1) * t307 + mrSges(7,2) * t127;
t55 = -pkin(5) * t79 + t144;
t51 = t118 * Ifges(5,4) + t119 * Ifges(5,2) + t203 * Ifges(5,6);
t50 = -qJ(6) * t168 + t68;
t49 = pkin(5) * t260 + qJ(6) * t170 + t67;
t29 = mrSges(6,1) * t196 - mrSges(6,3) * t40;
t28 = mrSges(7,1) * t196 - mrSges(7,3) * t40;
t27 = t44 - t483;
t26 = t43 - t452;
t17 = -pkin(5) * t41 + qJDD(6) + t72;
t14 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t8 = qJ(6) * t79 - qJD(6) * t168 + t15;
t7 = pkin(5) * t345 - qJ(6) * t77 + qJD(6) * t170 + t16;
t1 = [(t322 / 0.2e1 + t263 * t313) * t106 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t481 * t239 - t289 * mrSges(4,3) - t478 * t347 / 0.2e1 + m(5) * (t193 * t347 * t267 + t110 * t89 + t111 * t88 + t147 * t46 + t148 * t45) + t507 * t77 + (t110 * mrSges(5,1) - t111 * mrSges(5,2) + t486 * t405 + t467 * t419 + t465 * t422 - t498) * t345 + (mrSges(6,1) * t72 + mrSges(7,1) * t17 - mrSges(6,3) * t5 - mrSges(7,3) * t3 - t412 * t465 - t431 * t466 - t432 * t468 + t496) * t168 + (-m(5) * t138 * t267 + Ifges(4,5) * qJDD(3) + t292 * t411 + t295 * t425 + t298 * t426) * t264 + t500 * t79 + t484 * t354 + (Ifges(4,1) * t221 + Ifges(4,4) * t222) * t497 + (-t485 / 0.2e1 - t467 * t412 - t468 * t431 - t469 * t432 - t17 * mrSges(7,2) - t72 * mrSges(6,2) + t6 * mrSges(6,3) + t2 * mrSges(7,3)) * t170 + t50 * t30 + t49 * t28 - pkin(1) * t338 + m(7) * (t129 * t17 + t2 * t49 + t20 * t7 + t23 * t8 + t3 * t50 + t55 * t73) + m(6) * (t137 * t144 + t15 * t33 + t16 * t32 + t205 * t72 + t5 * t68 + t6 * t67) + (-t363 * t434 - m(3) * t349 - t184 * mrSges(5,1) - t183 * mrSges(5,2) + t488 * t163 - t487 * t162 + t473 * (t265 * pkin(7) + t349) - t472 * t265 + t471 * t261) * g(2) + m(4) * (t267 * t289 + t308) + m(3) * (-pkin(1) * t246 + t308) + (t110 * t279 - t111 * t480 - t356 * t46 - t364 * t45) * mrSges(5,3) + t193 * (mrSges(5,1) * t480 - mrSges(5,2) * t279) + t206 * (qJD(3) * t273 - t294 * t341) / 0.2e1 + t240 * (qJD(3) * t272 - t291 * t341) / 0.2e1 - t51 * t364 / 0.2e1 + t476 * t336 + t221 * t299 / 0.2e1 + qJD(3) ^ 2 * t293 / 0.2e1 + t356 * t430 + t222 * t509 + (qJD(3) * t274 - t297 * t341) * t409 + t259 * t107 * t313 + t67 * t29 + t68 * t31 + t55 * t69 + (t326 / 0.2e1 + t479 / 0.2e1 - Ifges(4,4) * t221 / 0.2e1 + t222 * t514 + t486 * t412 + t465 * t431 + t467 * t432 + Ifges(5,6) * t425 + Ifges(5,5) * t426 + Ifges(5,3) * t411 - Ifges(4,6) * qJDD(3) + t501 + t503) * t260 + (-t186 * mrSges(5,1) + t185 * mrSges(5,2) + t488 * t165 + t487 * t164 + (m(3) * pkin(1) + t267 * t473 + t472 + t489) * t261 + ((-m(3) + t473) * qJ(2) + t471) * t265) * g(1) + t8 * t95 + t15 * t96 + t7 * t97 + t16 * t98 + t228 * t321 + t129 * t13 + t285 * t310 + t172 * t359 + t144 * t70 + t88 * t145 + t89 * t146 + t147 * t86 + t148 * t87 + t138 * t197 + t205 * t14 + qJD(2) * t215 + qJ(2) * (-mrSges(4,1) * t222 + mrSges(4,2) * t221) + t246 * mrSges(3,2) + (t303 + 0.2e1 * mrSges(3,3)) * t238; t338 - t464 * t451 - (t29 + t28) * t167 + (qJ(2) * t427 - mrSges(3,3)) * t435 + (-t215 + t288) * qJD(1) + (-t13 - t14 + (t145 * t263 - t146 * t259 + t228) * qJD(3) + t484) * t264 + (t172 + t288 * qJD(4) + (t69 + t70 + t481) * qJD(3) + t445) * t260 + m(4) * t289 + m(3) * t246 + t454 * t457 + t453 * t456 + t442 * (-t427 + t477) + (-t167 * t2 - t17 * t264 + t20 * t453 + t23 * t454 - t3 * t451 + t73 * t347) * m(7) + (t137 * t347 - t167 * t6 - t264 * t72 + t32 * t453 + t33 * t454 - t451 * t5) * m(6) + ((-t138 + (-t110 * t259 + t111 * t263) * qJD(3)) * t264 + (qJD(3) * t193 + t270) * t260 - t290 * qJD(1)) * m(5); t485 * t211 / 0.2e1 + t240 * (t193 * t301 - t259 * t106 / 0.2e1) + t478 * t249 / 0.2e1 + t444 * mrSges(5,3) + (m(5) * t270 - t145 * t343 - t146 * t342 + t445) * pkin(8) + t448 * t70 + t450 * t69 + (t465 * t423 + t498 + t513) * t348 + (t157 * t465 + t158 * t467) * t406 - t442 * (t260 * t439 - t264 * t440 - t304) + (t157 * t466 + t158 * t468) * t423 + (t157 * t468 + t158 * t469) * t420 - t458 * t158 / 0.2e1 - t459 * t157 / 0.2e1 + t460 * t96 + (t137 * t448 + t140 * t6 + t141 * t5 - t245 * t72 + t32 * t461 + t33 * t460) * m(6) + t461 * t98 + t462 * t97 + (t102 * t2 + t103 * t3 + t160 * t17 + t20 * t462 + t23 * t463 + t450 * t73) * m(7) + t463 * t95 - t507 * t482 - t500 * t132 + t287 * t496 + (t107 / 0.2e1 - t377) * t342 + (-t111 * (mrSges(5,3) * t259 * t260 - mrSges(5,2) * t264) - t110 * (mrSges(5,1) * t264 + mrSges(5,3) * t361)) * qJD(1) - t481 * t223 + (t211 * t467 - t287 * t465) * t412 + (t211 * t468 - t287 * t466) * t431 + (t211 * t469 - t287 * t468) * t432 + (-t157 * t33 + t158 * t32 - t211 * t6 - t287 * t5) * mrSges(6,3) + (-t157 * t23 + t158 * t20 - t2 * t211 - t287 * t3) * mrSges(7,3) + t17 * (mrSges(7,1) * t287 + mrSges(7,2) * t211) + t72 * (mrSges(6,1) * t287 + mrSges(6,2) * t211) + (t260 * t440 + t264 * t439 + t303) * g(3) - t228 * t365 + (t285 / 0.2e1 - t476) * t435 + t138 * t302 + (-pkin(3) * t138 - t110 * t133 - t111 * t134 - t193 * t223) * m(5) + t259 * t430 + t294 * t425 + t297 * t426 + t291 * t411 + Ifges(4,3) * qJDD(3) + (t206 * t295 + t207 * t298 + t240 * t292) * qJD(4) / 0.2e1 - (t206 * t273 + t207 * t274 + t240 * t272) * qJD(1) / 0.2e1 - pkin(3) * t64 + t102 * t28 + t103 * t30 + t293 * t310 + t140 * t29 + t141 * t31 + t142 * mrSges(4,1) - t143 * mrSges(4,2) - t134 * t145 - t133 * t146 - t73 * (-mrSges(7,1) * t157 + mrSges(7,2) * t158) - t137 * (-mrSges(6,1) * t157 + mrSges(6,2) * t158) + t160 * t13 + Ifges(4,5) * t221 + Ifges(4,6) * t222 - t245 * t14 - t343 * t376 + t263 * t51 / 0.2e1; (-m(7) * (-t224 * t362 + t225 * t265) + mrSges(5,2) * t184 + t455 * t183 + t446) * g(1) + (-m(7) * (t224 * t360 + t225 * t261) - mrSges(5,2) * t186 + t455 * t185 + t447) * g(2) + (t406 * t467 + t420 * t469 + t423 * t468 + t436) * t307 + (mrSges(7,1) * t250 + t443 + t489) * t393 - t70 * t401 - m(6) * (t137 * t401 + t32 * t43 + t33 * t44) + t271 + (t2 * t244 - t20 * t26 + t224 * t393 - t23 * t27 - t73 * t92) * m(7) + ((t258 * t3 + (-t20 * t258 + t23 * t262) * qJD(5)) * m(7) - t456 * t340 + t457 * t339 + t464 * t258) * pkin(4) + t503 + t458 * t423 + t326 + (-t465 * t406 - t468 * t420 - t466 * t423 + t475) * t127 + (Ifges(5,1) * t206 - t373) * t511 + (-Ifges(5,2) * t207 + t107 + t198) * t512 + (t258 * t5 + t262 * t6 + (-t258 * t32 + t262 * t33) * qJD(5)) * t434 + (Ifges(5,5) * t206 - Ifges(5,6) * t207) * t510 + t106 * t409 - t92 * t69 - t27 * t95 - t44 * t96 - t26 * t97 - t43 * t98 + t207 * t376 + t206 * t377 + t29 * t398 - t110 * t145 + t111 * t146 + g(3) * t197 - t193 * (mrSges(5,1) * t207 + mrSges(5,2) * t206) + t244 * t28; t271 + t2 * t433 - t22 * t95 - t73 * (mrSges(7,1) * t127 + mrSges(7,2) * t307) + t20 * t383 - t137 * (mrSges(6,1) * t127 + mrSges(6,2) * t307) + (t307 * t469 - t493) * t420 + t459 * t419 + (-t127 * t465 + t307 * t467) * t406 + (-(-mrSges(7,1) - t433) * t250 + t443) * t393 + (t98 + t384) * t33 + (-t96 + t385) * t32 + (-m(7) * (-t20 + t22) + t97 + t382) * t23 + (-t164 * t433 + t447) * g(2) + (-t162 * t433 + t446) * g(1) + (-t127 * t466 + t458 + t495) * t423 + (t28 + (-m(7) * t73 - t69) * t127) * pkin(5); -t307 * t95 + t127 * t97 + (-g(3) * t260 + t20 * t127 - t23 * t307 - t264 * t442 + t17) * m(7) + t13;];
tau  = t1;
