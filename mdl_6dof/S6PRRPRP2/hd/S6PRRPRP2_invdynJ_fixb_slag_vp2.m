% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:25
% EndTime: 2019-03-08 21:29:05
% DurationCPUTime: 24.24s
% Computational Cost: add. (7274->713), mult. (16780->948), div. (0->0), fcn. (12749->14), ass. (0->322)
t473 = Ifges(6,1) + Ifges(7,1);
t474 = -Ifges(6,4) + Ifges(7,5);
t443 = Ifges(7,4) + Ifges(6,5);
t472 = Ifges(6,6) - Ifges(7,6);
t471 = -Ifges(6,3) - Ifges(7,2);
t241 = sin(pkin(11));
t247 = sin(qJ(3));
t250 = cos(qJ(3));
t378 = cos(pkin(11));
t261 = -t241 * t247 + t250 * t378;
t188 = t261 * qJD(2);
t453 = qJD(5) - t188;
t342 = qJD(2) * qJD(3);
t204 = qJDD(2) * t250 - t247 * t342;
t205 = qJDD(2) * t247 + t250 * t342;
t133 = t241 * t204 + t205 * t378;
t246 = sin(qJ(5));
t249 = cos(qJ(5));
t300 = t378 * t247;
t348 = qJD(2) * t250;
t189 = -qJD(2) * t300 - t241 * t348;
t269 = t249 * qJD(3) + t189 * t246;
t454 = qJD(5) * t269;
t65 = qJDD(3) * t246 + t133 * t249 + t454;
t425 = t65 / 0.2e1;
t147 = qJD(3) * t246 - t189 * t249;
t66 = qJD(5) * t147 - t249 * qJDD(3) + t133 * t246;
t423 = t66 / 0.2e1;
t449 = -m(6) - m(7);
t132 = t204 * t378 - t241 * t205;
t126 = qJDD(5) - t132;
t422 = t126 / 0.2e1;
t470 = t443 * t126 + t473 * t65 + t474 * t66;
t35 = mrSges(6,1) * t126 - mrSges(6,3) * t65;
t36 = -t126 * mrSges(7,1) + t65 * mrSges(7,2);
t401 = t36 - t35;
t37 = -mrSges(6,2) * t126 - mrSges(6,3) * t66;
t38 = -mrSges(7,2) * t66 + mrSges(7,3) * t126;
t400 = t37 + t38;
t469 = t147 * t443 + t269 * t472 - t471 * t453;
t145 = Ifges(6,4) * t269;
t384 = t269 * Ifges(7,5);
t440 = t147 * t473 + t443 * t453 + t145 - t384;
t86 = mrSges(7,2) * t269 + mrSges(7,3) * t453;
t393 = mrSges(6,3) * t269;
t87 = -mrSges(6,2) * t453 + t393;
t399 = t86 + t87;
t392 = mrSges(6,3) * t147;
t88 = mrSges(6,1) * t453 - t392;
t89 = -mrSges(7,1) * t453 + mrSges(7,2) * t147;
t398 = t88 - t89;
t394 = mrSges(5,3) * t189;
t468 = qJD(3) * mrSges(5,1) + mrSges(6,1) * t269 - mrSges(6,2) * t147 + t394;
t275 = pkin(5) * t246 - qJ(6) * t249;
t244 = cos(pkin(6));
t352 = qJD(1) * t244;
t229 = t250 * t352;
t248 = sin(qJ(2));
t243 = sin(pkin(6));
t353 = qJD(1) * t243;
t324 = t248 * t353;
t207 = qJD(2) * pkin(8) + t324;
t298 = qJ(4) * qJD(2) + t207;
t137 = -t247 * t298 + t229;
t322 = t247 * t352;
t138 = t250 * t298 + t322;
t301 = t378 * t138;
t69 = t137 * t241 + t301;
t467 = -qJD(6) * t246 + t275 * t453 - t69;
t466 = Ifges(5,5) * qJD(3);
t465 = Ifges(5,6) * qJD(3);
t251 = cos(qJ(2));
t351 = qJD(2) * t243;
t313 = qJD(1) * t351;
t220 = t251 * t313;
t340 = qJDD(1) * t243;
t177 = t248 * t340 + t220;
t163 = qJDD(2) * pkin(8) + t177;
t464 = qJD(3) * t352 + t163;
t285 = mrSges(7,1) * t246 - mrSges(7,3) * t249;
t287 = mrSges(6,1) * t246 + mrSges(6,2) * t249;
t127 = t241 * t138;
t131 = qJD(3) * pkin(3) + t137;
t58 = t131 * t378 - t127;
t53 = -qJD(3) * pkin(4) - t58;
t30 = -pkin(5) * t269 - t147 * qJ(6) + t53;
t462 = t30 * t285 + t53 * t287;
t461 = -t246 * t472 + t249 * t443;
t386 = Ifges(7,5) * t246;
t388 = Ifges(6,4) * t246;
t460 = t249 * t473 + t386 - t388;
t360 = t243 * t248;
t192 = t244 * t250 - t247 * t360;
t379 = cos(pkin(10));
t302 = t379 * t251;
t242 = sin(pkin(10));
t362 = t242 * t248;
t187 = -t244 * t362 + t302;
t359 = t243 * t250;
t459 = -t187 * t247 + t242 * t359;
t343 = qJD(5) * t249;
t372 = t188 * t249;
t458 = t343 - t372;
t344 = qJD(5) * t246;
t373 = t188 * t246;
t457 = -t344 + t373;
t339 = qJDD(1) * t244;
t226 = t250 * t339;
t341 = qJD(2) * qJD(4);
t346 = qJD(3) * t250;
t45 = -t207 * t346 + qJDD(3) * pkin(3) - qJ(4) * t205 + t226 + (-t341 - t464) * t247;
t347 = qJD(3) * t247;
t74 = -t207 * t347 + t247 * t339 + t250 * t464;
t46 = qJ(4) * t204 + t250 * t341 + t74;
t16 = t241 * t45 + t378 * t46;
t14 = qJDD(3) * pkin(9) + t16;
t219 = t248 * t313;
t176 = t251 * t340 - t219;
t162 = -qJDD(2) * pkin(2) - t176;
t112 = -pkin(3) * t204 + qJDD(4) + t162;
t39 = -pkin(4) * t132 - pkin(9) * t133 + t112;
t59 = t241 * t131 + t301;
t54 = qJD(3) * pkin(9) + t59;
t234 = pkin(3) * t250 + pkin(2);
t323 = t251 * t353;
t174 = -qJD(2) * t234 + qJD(4) - t323;
t83 = -pkin(4) * t188 + pkin(9) * t189 + t174;
t3 = t249 * t14 + t246 * t39 + t83 * t343 - t344 * t54;
t28 = t246 * t83 + t249 * t54;
t4 = -qJD(5) * t28 - t14 * t246 + t249 * t39;
t456 = -t246 * t4 + t249 * t3;
t1 = qJ(6) * t126 + qJD(6) * t453 + t3;
t2 = -pkin(5) * t126 + qJDD(6) - t4;
t455 = t1 * t249 + t2 * t246;
t452 = Ifges(7,5) * t425 + Ifges(7,6) * t422 - t65 * Ifges(6,4) / 0.2e1 - t126 * Ifges(6,6) / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t423;
t276 = t249 * pkin(5) + t246 * qJ(6);
t286 = -t249 * mrSges(7,1) - t246 * mrSges(7,3);
t288 = mrSges(6,1) * t249 - mrSges(6,2) * t246;
t451 = -m(7) * t276 - mrSges(5,1) + t286 - t288;
t303 = t379 * t248;
t361 = t242 * t251;
t185 = t244 * t303 + t361;
t240 = qJ(3) + pkin(11);
t236 = sin(t240);
t237 = cos(t240);
t304 = t243 * t379;
t120 = t185 * t237 - t236 * t304;
t363 = t242 * t243;
t122 = t187 * t237 + t236 * t363;
t166 = t236 * t244 + t237 * t360;
t450 = -g(1) * t122 - g(2) * t120 - g(3) * t166;
t198 = t241 * t250 + t300;
t190 = t198 * qJD(3);
t448 = -t190 / 0.2e1;
t191 = t261 * qJD(3);
t447 = t191 / 0.2e1;
t446 = t204 / 0.2e1;
t258 = -t185 * t247 - t250 * t304;
t254 = t258 * pkin(3);
t445 = mrSges(6,3) + mrSges(7,2);
t380 = qJDD(3) / 0.2e1;
t141 = t261 * t323;
t245 = -qJ(4) - pkin(8);
t305 = qJD(3) * t245;
t182 = qJD(4) * t250 + t247 * t305;
t183 = -qJD(4) * t247 + t250 * t305;
t98 = t182 * t378 + t241 * t183;
t439 = -t141 + t98;
t110 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t133;
t23 = mrSges(6,1) * t66 + mrSges(6,2) * t65;
t437 = t23 - t110;
t27 = -t246 * t54 + t249 * t83;
t436 = qJD(6) - t27;
t116 = -pkin(4) * t261 - pkin(9) * t198 - t234;
t215 = t245 * t247;
t216 = t245 * t250;
t143 = t241 * t215 - t216 * t378;
t435 = t246 * t116 + t249 * t143;
t265 = t191 * t246 + t198 * t343;
t334 = pkin(3) * t347;
t434 = -t324 + t334;
t433 = -t126 * t471 + t443 * t65 - t472 * t66;
t152 = t207 * t250 + t322;
t75 = -t152 * qJD(3) - t163 * t247 + t226;
t432 = -t247 * t75 + t250 * t74;
t333 = m(4) * pkin(8) + mrSges(4,3);
t430 = -mrSges(5,3) - t333 + mrSges(3,2);
t81 = -mrSges(7,1) * t269 - mrSges(7,3) * t147;
t338 = t81 - t468;
t429 = 0.2e1 * t380;
t290 = -mrSges(4,1) * t250 + mrSges(4,2) * t247;
t259 = m(4) * pkin(2) - t290;
t289 = t237 * mrSges(5,1) - mrSges(5,2) * t236;
t428 = t289 + t259 + mrSges(3,1);
t100 = pkin(4) * t190 - pkin(9) * t191 + t334;
t18 = -qJD(5) * t435 + t100 * t249 - t246 * t98;
t297 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t294 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t427 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t252 = qJD(2) ^ 2;
t424 = -t66 / 0.2e1;
t421 = t269 / 0.2e1;
t420 = -t269 / 0.2e1;
t419 = -t147 / 0.2e1;
t418 = t147 / 0.2e1;
t417 = -t453 / 0.2e1;
t415 = -t188 / 0.2e1;
t414 = -t189 / 0.2e1;
t413 = t189 / 0.2e1;
t410 = t246 / 0.2e1;
t409 = pkin(3) * t241;
t408 = pkin(4) * t237;
t407 = g(3) * t243;
t71 = t137 * t378 - t127;
t350 = qJD(2) * t247;
t335 = pkin(3) * t350;
t99 = -pkin(4) * t189 - pkin(9) * t188 + t335;
t32 = t246 * t99 + t249 * t71;
t395 = mrSges(5,3) * t188;
t391 = Ifges(4,4) * t247;
t390 = Ifges(4,4) * t250;
t389 = Ifges(5,4) * t189;
t387 = Ifges(6,4) * t249;
t385 = Ifges(7,5) * t249;
t383 = t147 * Ifges(6,4);
t184 = -t244 * t302 + t362;
t377 = t184 * t236;
t186 = t244 * t361 + t303;
t375 = t186 * t236;
t370 = t191 * t249;
t365 = t237 * t246;
t364 = t237 * t249;
t358 = t243 * t251;
t356 = t249 * t251;
t355 = -t184 * t234 - t185 * t245;
t354 = -t186 * t234 - t187 * t245;
t349 = qJD(2) * t248;
t332 = mrSges(4,3) * t350;
t331 = mrSges(4,3) * t348;
t330 = t236 * t358;
t327 = t243 * t356;
t223 = t246 * t358;
t144 = Ifges(7,5) * t147;
t47 = Ifges(7,6) * t453 - Ifges(7,3) * t269 + t144;
t326 = t47 * t410;
t325 = t378 * pkin(3);
t321 = t243 * t349;
t320 = t251 * t351;
t308 = -t344 / 0.2e1;
t307 = t343 / 0.2e1;
t67 = -t132 * mrSges(5,1) + t133 * mrSges(5,2);
t97 = t182 * t241 - t378 * t183;
t142 = -t378 * t215 - t216 * t241;
t292 = t459 * pkin(3);
t233 = -t325 - pkin(4);
t15 = -t241 * t46 + t378 * t45;
t282 = t250 * Ifges(4,2) + t391;
t281 = -Ifges(6,2) * t246 + t387;
t279 = Ifges(4,5) * t250 - Ifges(4,6) * t247;
t277 = Ifges(7,3) * t246 + t385;
t31 = -t246 * t71 + t249 * t99;
t55 = t116 * t249 - t143 * t246;
t268 = t192 * pkin(3);
t193 = t244 * t247 + t248 * t359;
t105 = t241 * t192 + t193 * t378;
t84 = t105 * t246 + t327;
t264 = t198 * t344 - t370;
t208 = -qJD(2) * pkin(2) - t323;
t263 = t208 * (mrSges(4,1) * t247 + mrSges(4,2) * t250);
t262 = t247 * (Ifges(4,1) * t250 - t391);
t17 = t246 * t100 + t116 * t343 - t143 * t344 + t249 * t98;
t13 = -qJDD(3) * pkin(4) - t15;
t256 = -g(1) * t186 - g(2) * t184 + g(3) * t358;
t235 = Ifges(4,4) * t348;
t214 = -qJD(3) * mrSges(4,2) + t331;
t213 = qJD(3) * mrSges(4,1) - t332;
t203 = t234 * t358;
t200 = t290 * qJD(2);
t196 = -t276 + t233;
t195 = Ifges(4,1) * t350 + Ifges(4,5) * qJD(3) + t235;
t194 = Ifges(4,6) * qJD(3) + qJD(2) * t282;
t180 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t205;
t179 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t204;
t175 = Ifges(5,4) * t188;
t165 = -t236 * t360 + t237 * t244;
t160 = -qJD(3) * mrSges(5,2) + t395;
t151 = -t207 * t247 + t229;
t139 = -mrSges(4,1) * t204 + mrSges(4,2) * t205;
t136 = qJD(3) * t192 + t250 * t320;
t135 = -qJD(3) * t193 - t247 * t320;
t123 = t166 * t246 + t327;
t121 = -t187 * t236 + t237 * t363;
t119 = -t185 * t236 - t237 * t304;
t111 = -mrSges(5,1) * t188 - mrSges(5,2) * t189;
t109 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t132;
t107 = -t189 * Ifges(5,1) + t175 + t466;
t106 = t188 * Ifges(5,2) - t389 + t465;
t104 = -t192 * t378 + t193 * t241;
t102 = t141 * t249 + t246 * t324;
t101 = t141 * t246 - t249 * t324;
t85 = t105 * t249 - t223;
t80 = pkin(5) * t147 - qJ(6) * t269;
t78 = t122 * t246 - t186 * t249;
t76 = t120 * t246 - t184 * t249;
t72 = t198 * t275 + t142;
t70 = t241 * t135 + t136 * t378;
t68 = -t135 * t378 + t136 * t241;
t50 = Ifges(6,2) * t269 + Ifges(6,6) * t453 + t383;
t41 = pkin(5) * t261 - t55;
t40 = -qJ(6) * t261 + t435;
t29 = pkin(5) * t189 - t31;
t26 = -qJ(6) * t189 + t32;
t25 = -qJD(5) * t223 + t105 * t343 + t246 * t70 - t249 * t321;
t24 = -qJD(5) * t84 + t246 * t321 + t249 * t70;
t22 = mrSges(7,1) * t66 - mrSges(7,3) * t65;
t21 = t275 * t191 + (qJD(5) * t276 - qJD(6) * t249) * t198 + t97;
t20 = qJ(6) * t453 + t28;
t19 = -pkin(5) * t453 + t436;
t7 = -pkin(5) * t190 - t18;
t6 = qJ(6) * t190 - qJD(6) * t261 + t17;
t5 = t66 * pkin(5) - t65 * qJ(6) - t147 * qJD(6) + t13;
t8 = [m(2) * qJDD(1) + t105 * t109 + t135 * t213 + t136 * t214 + t70 * t160 + t193 * t179 + t192 * t180 + t400 * t85 + t401 * t84 - t398 * t25 + t399 * t24 + t338 * t68 + (t22 + t437) * t104 + (-m(2) - m(3) - m(4) - m(5) + t449) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t252 - t139 - t67) * t251 + (-mrSges(3,1) * t252 - mrSges(3,2) * qJDD(2) + (t111 + t200) * qJD(2)) * t248) * t243 + m(5) * (-t104 * t15 + t105 * t16 - t58 * t68 + t59 * t70 + (-t112 * t251 + t174 * t349) * t243) + m(4) * (t135 * t151 + t136 * t152 + t192 * t75 + t193 * t74 + (-t162 * t251 + t208 * t349) * t243) + m(3) * (qJDD(1) * t244 ^ 2 + (t176 * t251 + t177 * t248) * t243) + m(7) * (t1 * t85 + t104 * t5 + t19 * t25 + t2 * t84 + t20 * t24 + t30 * t68) + m(6) * (t104 * t13 + t24 * t28 - t25 * t27 + t3 * t85 - t4 * t84 + t53 * t68); (-Ifges(6,4) * t264 - Ifges(6,2) * t265 + Ifges(6,6) * t190) * t421 + (Ifges(5,1) * t191 - Ifges(5,4) * t190) * t414 + (-Ifges(7,5) * t264 + Ifges(7,6) * t190 + Ifges(7,3) * t265) * t420 + (t262 + t250 * (-Ifges(4,2) * t247 + t390)) * t342 / 0.2e1 + (m(4) * ((-t151 * t250 - t152 * t247) * qJD(3) + t432) - t213 * t346 - t214 * t347 + t250 * t179 - t247 * t180) * pkin(8) + (-t151 * t346 - t152 * t347 + t432) * mrSges(4,3) + t434 * t111 - t265 * t50 / 0.2e1 + (-t190 * t59 - t191 * t58) * mrSges(5,3) - t200 * t324 + t469 * t190 / 0.2e1 + (t13 * t142 + t3 * t435 + t4 * t55 + t53 * t97 + (-t102 + t17) * t28 + (t101 + t18) * t27) * m(6) + t435 * t37 + t195 * t346 / 0.2e1 - t194 * t347 / 0.2e1 + t162 * t290 + ((t470 / 0.2e1 + mrSges(7,2) * t2 - mrSges(6,3) * t4) * t249 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) + t452) * t246 + t112 * mrSges(5,2) - t15 * mrSges(5,3) + Ifges(5,1) * t133 + Ifges(5,4) * t132 + t429 * Ifges(5,5) + t13 * t287 + t277 * t423 + t281 * t424 + t5 * t285 + t47 * t307 + t422 * t461 + t425 * t460 + (m(5) * t58 - m(6) * t53 - m(7) * t30 - t338) * t323 + t440 * t308) * t198 - t468 * t97 + (-t190 * t471 - t264 * t443 - t265 * t472) * t453 / 0.2e1 + (-t433 / 0.2e1 - t112 * mrSges(5,1) + t16 * mrSges(5,3) + Ifges(5,4) * t133 + Ifges(5,2) * t132 + Ifges(5,6) * t429 - Ifges(6,6) * t424 - Ifges(7,6) * t423 + t422 * t471 - t425 * t443 - t427) * t261 + t282 * t446 + t107 * t447 + t106 * t448 + t27 * (mrSges(6,1) * t190 + mrSges(6,3) * t264) + t19 * (-mrSges(7,1) * t190 - mrSges(7,2) * t264) + t30 * (mrSges(7,1) * t265 + mrSges(7,3) * t264) + t53 * (mrSges(6,1) * t265 - mrSges(6,2) * t264) + t28 * (-mrSges(6,2) * t190 - mrSges(6,3) * t265) + t20 * (-mrSges(7,2) * t265 + mrSges(7,3) * t190) + (t248 * t407 - t177 + t220) * mrSges(3,2) + (t1 * t40 + t2 * t41 + t21 * t30 + t5 * t72 + (-t102 + t6) * t20 + (-t101 + t7) * t19) * m(7) - t250 * t214 * t323 + t205 * t390 / 0.2e1 + t250 * (Ifges(4,4) * t205 + Ifges(4,2) * t204 + Ifges(4,6) * qJDD(3)) / 0.2e1 + Ifges(3,3) * qJDD(2) - t234 * t67 + (t443 * t190 - t473 * t264 + t265 * t474) * t418 + t188 * (Ifges(5,4) * t191 - Ifges(5,2) * t190) / 0.2e1 + t174 * (mrSges(5,1) * t190 + mrSges(5,2) * t191) + Ifges(4,6) * t250 * t380 + t143 * t109 - pkin(2) * t139 + t17 * t87 + t18 * t88 + t7 * t89 + t6 * t86 + t21 * t81 + t72 * t22 + t55 * t35 + t41 * t36 + t40 * t38 + (-m(5) * t203 - t445 * t330 + t449 * (pkin(9) * t330 - t245 * t360 + t358 * t408 + t203) + t294 * (t223 * t237 - t249 * t360) + (-t289 * t251 - (-m(5) * t245 + mrSges(5,3)) * t248 - t297 * (t237 * t356 + t246 * t248)) * t243) * g(3) + (-t251 * t407 + t176 + t219) * mrSges(3,1) + t191 * t326 - (t248 * t333 + t251 * t259) * t407 + t437 * t142 + (-t112 * t234 - t142 * t15 + t143 * t16 + t174 * t434 + t439 * t59 - t58 * t97) * m(5) + t439 * t160 + t398 * t101 - t399 * t102 + t440 * t370 / 0.2e1 + (-pkin(2) * t162 - (t208 * t248 + (-t151 * t247 + t152 * t250) * t251) * t353) * m(4) + (Ifges(4,1) * t205 + Ifges(4,4) * t446 + Ifges(4,5) * t429 + t213 * t323) * t247 + (t263 + Ifges(5,5) * t447 + Ifges(5,6) * t448 + t279 * qJD(3) / 0.2e1) * qJD(3) + (-m(5) * t354 + t449 * (-pkin(9) * t375 - t186 * t408 + t354) - t297 * (-t186 * t364 + t187 * t246) + t294 * (-t186 * t365 - t187 * t249) + t445 * t375 + t430 * t187 + t428 * t186) * g(1) + (-m(5) * t355 + t449 * (-pkin(9) * t377 - t184 * t408 + t355) - t297 * (-t184 * t364 + t185 * t246) + t294 * (-t184 * t365 - t185 * t249) + t445 * t377 + t430 * t185 + t428 * t184) * g(2); (t13 * t233 - t27 * t31 - t28 * t32 - t53 * t69) * m(6) + (-t19 * t29 + t196 * t5 - t20 * t26 + t467 * t30) * m(7) - (-Ifges(4,2) * t350 + t195 + t235) * t348 / 0.2e1 + (-t372 / 0.2e1 + t307) * t440 + (-t174 * t335 + t58 * t69 - t59 * t71 + (t15 * t378 + t16 * t241) * pkin(3)) * m(5) + (-t385 + t387) * t425 + t109 * t409 + t106 * t414 - t279 * t342 / 0.2e1 + t50 * t308 + (t107 + t175) * t415 + (t332 + t213) * t152 + (t281 / 0.2e1 - t277 / 0.2e1) * t454 + t386 * t423 + t388 * t424 + (t331 - t214) * t151 + t110 * t325 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (t389 + t469) * t413 + t470 * t410 + (Ifges(6,2) * t424 - Ifges(7,3) * t423 + t422 * t472 - t452) * t249 + (-Ifges(7,6) * t421 + Ifges(5,2) * t415 - Ifges(6,6) * t420 + t27 * mrSges(6,1) - t19 * mrSges(7,1) - t28 * mrSges(6,2) + t20 * mrSges(7,3) + t174 * mrSges(5,1) - t465 / 0.2e1 - t443 * t419 + t471 * t417) * t189 + t58 * t395 + (t19 * t458 + t20 * t457 + t450 + t455) * mrSges(7,2) + (-t27 * t458 + t28 * t457 + t450 + t456) * mrSges(6,3) + (t422 * t443 + t425 * t473) * t246 + (t50 / 0.2e1 - t47 / 0.2e1) * t373 + t194 * t350 / 0.2e1 + (t326 + t462) * qJD(5) + (-m(5) * t268 - mrSges(4,1) * t192 + mrSges(4,2) * t193 + t166 * mrSges(5,2) + t449 * (t165 * pkin(4) + pkin(9) * t166 + t268) + t451 * t165) * g(3) + (-t459 * mrSges(4,1) - (-t187 * t250 - t247 * t363) * mrSges(4,2) - m(5) * t292 + t122 * mrSges(5,2) + t449 * (t121 * pkin(4) + pkin(9) * t122 + t292) + t451 * t121) * g(1) + t5 * t286 - t13 * t288 + (((-t246 * t28 - t249 * t27) * qJD(5) + t456) * m(6) - t398 * t343 - t399 * t344 + t400 * t249 + t401 * t246 + ((t19 * t249 - t20 * t246) * qJD(5) + t455) * m(7)) * (pkin(9) + t409) + (t147 * t460 + t453 * t461) * qJD(5) / 0.2e1 - t111 * t335 + (t277 * t421 + Ifges(5,1) * t413 + t281 * t420 - t174 * mrSges(5,2) - t466 / 0.2e1 + t460 * t419 + t461 * t417 - t462) * t188 + t467 * t81 + t468 * t69 - qJD(2) * t263 - t59 * t394 + t233 * t23 + Ifges(4,5) * t205 + Ifges(4,6) * t204 + t196 * t22 - t71 * t160 + Ifges(5,5) * t133 + Ifges(5,6) * t132 - t31 * t88 - t29 * t89 - t26 * t86 - t32 * t87 - t74 * mrSges(4,2) + t75 * mrSges(4,1) + t15 * mrSges(5,1) - t16 * mrSges(5,2) + (-t258 * mrSges(4,1) - (-t185 * t250 + t247 * t304) * mrSges(4,2) - m(5) * t254 + t120 * mrSges(5,2) + t451 * t119 + t449 * (t119 * pkin(4) + t120 * pkin(9) + t254)) * g(2) - t252 * t262 / 0.2e1; -t188 * t160 + t338 * t189 + (t399 * t453 - t401) * t249 + (-t398 * t453 + t400) * t246 + t67 + (t1 * t246 + t189 * t30 - t2 * t249 + t256 + t453 * (t19 * t246 + t20 * t249)) * m(7) + (t189 * t53 + t3 * t246 + t4 * t249 + t256 + t453 * (-t246 * t27 + t249 * t28)) * m(6) + (-t188 * t59 - t189 * t58 + t112 + t256) * m(5); (-Ifges(6,2) * t147 + t145 + t440) * t420 + t50 * t418 + t427 + (t147 * t20 - t19 * t269) * mrSges(7,2) - t53 * (mrSges(6,1) * t147 + mrSges(6,2) * t269) - t30 * (mrSges(7,1) * t147 - mrSges(7,3) * t269) + (-t147 * t472 + t269 * t443) * t417 + (t269 * t473 + t144 - t383 + t47) * t419 + (Ifges(7,3) * t147 + t384) * t421 + (t393 - t399) * t27 + (t294 * (t120 * t249 + t184 * t246) + t297 * t76) * g(2) + (t294 * (t122 * t249 + t186 * t246) + t297 * t78) * g(1) + (t294 * (t166 * t249 - t223) + t297 * t123) * g(3) + (-pkin(5) * t2 + qJ(6) * t1 - t19 * t28 + t20 * t436 - t30 * t80) * m(7) + qJD(6) * t86 - t80 * t81 + qJ(6) * t38 - pkin(5) * t36 + (t392 + t398) * t28 + t433; t147 * t81 - t453 * t86 + (-g(1) * t78 - g(2) * t76 - g(3) * t123 + t147 * t30 - t20 * t453 + t2) * m(7) + t36;];
tau  = t8;
