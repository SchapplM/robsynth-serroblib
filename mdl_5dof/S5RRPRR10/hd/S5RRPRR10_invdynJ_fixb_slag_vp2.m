% Calculate vector of inverse dynamics joint torques for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:49
% EndTime: 2019-12-31 20:24:33
% DurationCPUTime: 22.69s
% Computational Cost: add. (10396->738), mult. (29201->1037), div. (0->0), fcn. (23499->12), ass. (0->346)
t246 = sin(pkin(5));
t245 = sin(pkin(10));
t250 = sin(qJ(2));
t254 = cos(qJ(2));
t368 = cos(pkin(10));
t269 = -t250 * t245 + t254 * t368;
t261 = t246 * t269;
t185 = qJD(1) * t261;
t436 = t185 - qJD(4);
t396 = pkin(2) * t245;
t242 = pkin(8) + t396;
t249 = sin(qJ(4));
t342 = qJD(4) * t249;
t247 = cos(pkin(5));
t358 = t247 * t250;
t238 = pkin(1) * t358;
t360 = t246 * t254;
t389 = pkin(7) + qJ(3);
t177 = (t360 * t389 + t238) * qJD(1);
t167 = t245 * t177;
t397 = pkin(1) * t247;
t239 = t254 * t397;
t232 = qJD(1) * t239;
t317 = t389 * t250;
t301 = t246 * t317;
t176 = -qJD(1) * t301 + t232;
t107 = t176 * t368 - t167;
t313 = t250 * t368;
t346 = qJD(1) * t246;
t324 = t254 * t346;
t186 = -t245 * t324 - t313 * t346;
t325 = t250 * t346;
t306 = pkin(2) * t325;
t118 = -pkin(3) * t186 - pkin(8) * t185 + t306;
t253 = cos(qJ(4));
t59 = t253 * t107 + t249 * t118;
t461 = -t242 * t342 - t59;
t311 = t368 * t177;
t106 = t176 * t245 + t311;
t460 = -t106 - t436 * (pkin(4) * t249 - pkin(9) * t253);
t459 = -pkin(9) * t186 - t461;
t235 = qJD(1) * t247 + qJD(2);
t157 = t186 * t249 + t235 * t253;
t156 = qJD(5) - t157;
t181 = -t269 * t346 + qJD(4);
t248 = sin(qJ(5));
t252 = cos(qJ(5));
t275 = t186 * t253 - t235 * t249;
t94 = t181 * t252 + t248 * t275;
t95 = t181 * t248 - t252 * t275;
t36 = t95 * Ifges(6,5) + t94 * Ifges(6,6) + t156 * Ifges(6,3);
t379 = t275 * Ifges(5,4);
t456 = t181 * Ifges(5,6);
t457 = t157 * Ifges(5,2);
t69 = -t379 + t456 + t457;
t453 = -t69 / 0.2e1 + t36 / 0.2e1;
t159 = pkin(2) * t235 + t176;
t90 = t159 * t368 - t167;
t81 = -t235 * pkin(3) - t90;
t458 = t81 * mrSges(5,2);
t356 = t248 * t253;
t125 = -t185 * t356 - t186 * t252;
t341 = qJD(4) * t253;
t455 = t248 * t341 + t125;
t255 = cos(qJ(1));
t349 = t254 * t255;
t251 = sin(qJ(1));
t352 = t251 * t250;
t454 = t247 * t349 - t352;
t334 = qJDD(1) * t247;
t234 = qJDD(2) + t334;
t347 = pkin(7) * t360 + t238;
t200 = t347 * qJD(2);
t336 = qJD(1) * qJD(2);
t202 = (qJDD(1) * t250 + t254 * t336) * t246;
t330 = pkin(1) * t334;
t230 = t254 * t330;
t362 = t246 * t250;
t322 = qJD(3) * t362;
t335 = qJDD(1) * t246;
t329 = pkin(7) * t335;
t87 = -t250 * t329 + pkin(2) * t234 - qJ(3) * t202 + t230 + (-t200 - t322) * qJD(1);
t201 = (qJDD(1) * t254 - t250 * t336) * t246;
t333 = qJD(2) * t397;
t302 = qJD(1) * t333;
t327 = t250 * t330 + (t302 + t329) * t254;
t343 = qJD(3) * t254;
t344 = qJD(2) * t250;
t96 = qJ(3) * t201 + (-pkin(7) * t344 + t343) * t346 + t327;
t45 = -t245 * t96 + t368 * t87;
t42 = -t234 * pkin(3) - t45;
t140 = t245 * t201 + t202 * t368;
t77 = qJD(4) * t157 + t140 * t253 + t234 * t249;
t78 = qJD(4) * t275 - t140 * t249 + t234 * t253;
t13 = -t78 * pkin(4) - t77 * pkin(9) + t42;
t244 = pkin(2) * t254 + pkin(1);
t204 = -t244 * t346 + qJD(3);
t101 = -pkin(3) * t185 + pkin(8) * t186 + t204;
t91 = t245 * t159 + t311;
t82 = pkin(8) * t235 + t91;
t48 = t101 * t249 + t253 * t82;
t40 = pkin(9) * t181 + t48;
t44 = -t157 * pkin(4) + pkin(9) * t275 + t81;
t14 = -t248 * t40 + t252 * t44;
t46 = t245 * t87 + t368 * t96;
t43 = pkin(8) * t234 + t46;
t139 = t201 * t368 - t245 * t202;
t172 = -pkin(1) * t335 - pkin(2) * t201 + qJDD(3);
t64 = -pkin(3) * t139 - pkin(8) * t140 + t172;
t11 = t101 * t341 + t249 * t64 + t253 * t43 - t342 * t82;
t138 = qJDD(4) - t139;
t8 = pkin(9) * t138 + t11;
t1 = qJD(5) * t14 + t13 * t248 + t252 * t8;
t15 = t248 * t44 + t252 * t40;
t2 = -qJD(5) * t15 + t13 * t252 - t248 * t8;
t452 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t451 = t81 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2);
t409 = t138 / 0.2e1;
t415 = t78 / 0.2e1;
t416 = t77 / 0.2e1;
t426 = Ifges(5,1) * t416 + Ifges(5,4) * t415 + Ifges(5,5) * t409;
t32 = qJD(5) * t94 + t138 * t248 + t252 * t77;
t425 = t32 / 0.2e1;
t33 = -qJD(5) * t95 + t138 * t252 - t248 * t77;
t424 = t33 / 0.2e1;
t76 = qJDD(5) - t78;
t417 = t76 / 0.2e1;
t410 = -m(6) - m(5);
t47 = t101 * t253 - t249 * t82;
t450 = t47 * mrSges(5,1);
t449 = t48 * mrSges(5,2);
t10 = -mrSges(6,1) * t33 + mrSges(6,2) * t32;
t51 = mrSges(5,1) * t138 - mrSges(5,3) * t77;
t448 = t10 - t51;
t291 = -mrSges(6,1) * t252 + mrSges(6,2) * t248;
t266 = m(6) * pkin(4) - t291;
t447 = mrSges(5,1) + t266;
t290 = t248 * mrSges(6,1) + t252 * mrSges(6,2);
t446 = -m(6) * pkin(8) + mrSges(4,2) - t290;
t331 = m(6) * pkin(9) + mrSges(6,3);
t307 = mrSges(5,2) - t331;
t387 = mrSges(5,3) * t275;
t100 = mrSges(5,1) * t181 + t387;
t57 = -mrSges(6,1) * t94 + mrSges(6,2) * t95;
t445 = t100 - t57;
t326 = t368 * pkin(2);
t243 = -t326 - pkin(3);
t212 = -t253 * pkin(4) - t249 * pkin(9) + t243;
t165 = t212 * t252 - t242 * t356;
t444 = qJD(5) * t165 + t248 * t460 - t459 * t252;
t350 = t252 * t253;
t166 = t212 * t248 + t242 * t350;
t443 = -qJD(5) * t166 + t459 * t248 + t252 * t460;
t377 = t186 * mrSges(4,3);
t370 = -mrSges(4,1) * t235 - mrSges(5,1) * t157 - mrSges(5,2) * t275 - t377;
t338 = qJD(5) * t252;
t442 = t249 * t338 + t455;
t126 = t185 * t350 - t186 * t248;
t339 = qJD(5) * t249;
t441 = t248 * t339 - t252 * t341 + t126;
t365 = t185 * t249;
t440 = t342 - t365;
t16 = mrSges(6,1) * t76 - mrSges(6,3) * t32;
t17 = -mrSges(6,2) * t76 + mrSges(6,3) * t33;
t439 = -t248 * t16 + t252 * t17;
t12 = -qJD(4) * t48 - t249 * t43 + t253 * t64;
t438 = t11 * t253 - t12 * t249;
t437 = t1 * t252 - t2 * t248;
t386 = Ifges(3,4) * t250;
t435 = pkin(1) * (mrSges(3,1) * t250 + mrSges(3,2) * t254) - t250 * (Ifges(3,1) * t254 - t386) / 0.2e1;
t293 = mrSges(5,1) * t253 - mrSges(5,2) * t249;
t434 = -t249 * t331 - t253 * t266 - mrSges(4,1) - t293;
t433 = -pkin(8) * t410 - mrSges(4,2) + mrSges(5,3);
t432 = t12 * mrSges(5,1) - t11 * mrSges(5,2) + Ifges(5,5) * t77 + Ifges(5,6) * t78 + Ifges(5,3) * t138;
t214 = -t254 * t245 - t313;
t193 = t214 * t246;
t187 = qJD(2) * t193;
t188 = qJD(2) * t261;
t323 = t246 * t344;
t119 = pkin(2) * t323 - pkin(3) * t187 - pkin(8) * t188;
t173 = pkin(2) * t247 + t239 - t301;
t182 = qJ(3) * t360 + t347;
t114 = t245 * t173 + t368 * t182;
t103 = pkin(8) * t247 + t114;
t236 = pkin(2) * t360;
t348 = pkin(3) * t261 + t236;
t298 = -pkin(8) * t193 + t348;
t398 = pkin(1) * t246;
t124 = -t298 - t398;
t369 = t253 * t103 + t249 * t124;
t233 = t254 * t333;
t160 = t233 + (-qJD(2) * t317 + t343) * t246;
t318 = t389 * t246;
t161 = -t322 + (-t254 * t318 - t238) * qJD(2);
t85 = t160 * t368 + t245 * t161;
t25 = -qJD(4) * t369 + t119 * t253 - t249 * t85;
t431 = m(5) * pkin(8) - t446;
t305 = pkin(7) * t323;
t153 = -qJD(1) * t305 + t327;
t154 = -pkin(7) * t202 - t250 * t302 + t230;
t430 = t154 * mrSges(3,1) + t45 * mrSges(4,1) - t153 * mrSges(3,2) - t46 * mrSges(4,2) + Ifges(3,5) * t202 + Ifges(4,5) * t140 + Ifges(3,6) * t201 + Ifges(4,6) * t139;
t408 = -t156 / 0.2e1;
t412 = -t95 / 0.2e1;
t414 = -t94 / 0.2e1;
t429 = Ifges(6,5) * t412 + Ifges(6,6) * t414 + Ifges(6,3) * t408;
t6 = Ifges(6,4) * t32 + Ifges(6,2) * t33 + Ifges(6,6) * t76;
t428 = t6 / 0.2e1;
t427 = Ifges(6,1) * t425 + Ifges(6,4) * t424 + Ifges(6,5) * t417;
t399 = Ifges(6,4) * t95;
t37 = t94 * Ifges(6,2) + t156 * Ifges(6,6) + t399;
t422 = -t37 / 0.2e1;
t421 = t37 / 0.2e1;
t92 = Ifges(6,4) * t94;
t38 = t95 * Ifges(6,1) + t156 * Ifges(6,5) + t92;
t420 = -t38 / 0.2e1;
t419 = t38 / 0.2e1;
t413 = t94 / 0.2e1;
t411 = t95 / 0.2e1;
t407 = t156 / 0.2e1;
t406 = -t157 / 0.2e1;
t405 = t275 / 0.2e1;
t404 = -t275 / 0.2e1;
t403 = -t181 / 0.2e1;
t401 = -t186 / 0.2e1;
t9 = -pkin(4) * t138 - t12;
t393 = t249 * t9;
t388 = mrSges(5,3) * t157;
t385 = Ifges(5,4) * t249;
t384 = Ifges(5,4) * t253;
t383 = Ifges(6,4) * t248;
t382 = Ifges(6,4) * t252;
t378 = t185 * mrSges(4,3);
t376 = t186 * Ifges(4,4);
t375 = t235 * Ifges(3,5);
t374 = t235 * Ifges(3,6);
t367 = t157 * t248;
t366 = t157 * t252;
t364 = t185 * t253;
t361 = t246 * t251;
t359 = t246 * t255;
t357 = t248 * t249;
t355 = t249 * t252;
t353 = t250 * t255;
t351 = t251 * t254;
t340 = qJD(5) * t248;
t337 = m(4) - t410;
t5 = Ifges(6,5) * t32 + Ifges(6,6) * t33 + Ifges(6,3) * t76;
t332 = -m(3) * pkin(1) - mrSges(2,1);
t320 = t242 * t341;
t315 = t341 / 0.2e1;
t314 = -t339 / 0.2e1;
t310 = -t139 * mrSges(4,1) + t140 * mrSges(4,2);
t194 = t214 * t247;
t147 = -t194 * t255 + t251 * t269;
t128 = t147 * t253 - t249 * t359;
t127 = -t147 * t249 - t253 * t359;
t84 = t160 * t245 - t368 * t161;
t203 = pkin(2) * t358 - t318;
t309 = -t203 * t251 + t255 * t244;
t304 = mrSges(3,3) * t325;
t303 = mrSges(3,3) * t324;
t299 = t454 * pkin(2);
t148 = -t251 * t194 - t255 * t269;
t295 = -pkin(3) * t148 + t309;
t170 = -t193 * t249 - t247 * t253;
t171 = -t193 * t253 + t247 * t249;
t294 = mrSges(5,1) * t170 + mrSges(5,2) * t171;
t120 = -t171 * t248 - t252 * t261;
t121 = t171 * t252 - t248 * t261;
t292 = mrSges(6,1) * t120 - mrSges(6,2) * t121;
t289 = Ifges(5,1) * t253 - t385;
t288 = Ifges(6,1) * t252 - t383;
t287 = Ifges(6,1) * t248 + t382;
t286 = -Ifges(5,2) * t249 + t384;
t285 = -Ifges(6,2) * t248 + t382;
t284 = Ifges(6,2) * t252 + t383;
t283 = Ifges(5,5) * t253 - Ifges(5,6) * t249;
t282 = Ifges(6,5) * t252 - Ifges(6,6) * t248;
t281 = Ifges(6,5) * t248 + Ifges(6,6) * t252;
t56 = -pkin(9) * t261 + t369;
t113 = t173 * t368 - t245 * t182;
t102 = -t247 * pkin(3) - t113;
t60 = t170 * pkin(4) - t171 * pkin(9) + t102;
t21 = t248 * t60 + t252 * t56;
t20 = -t248 * t56 + t252 * t60;
t61 = -t103 * t249 + t124 * t253;
t58 = -t107 * t249 + t118 * t253;
t207 = -t247 * t351 - t353;
t39 = -pkin(4) * t181 - t47;
t272 = t39 * t290;
t24 = -t103 * t342 + t249 * t119 + t124 * t341 + t253 * t85;
t264 = t207 * pkin(2);
t260 = t247 * t269;
t257 = (-t14 * t252 - t15 * t248) * qJD(5) + t437;
t231 = Ifges(3,4) * t324;
t223 = Ifges(3,3) * t234;
t222 = Ifges(4,3) * t234;
t215 = -t236 - t398;
t210 = -pkin(7) * t362 + t239;
t209 = (-mrSges(3,1) * t254 + mrSges(3,2) * t250) * t246;
t208 = -t247 * t352 + t349;
t206 = -t247 * t353 - t351;
t199 = t233 - t305;
t198 = t347 * qJD(1);
t197 = -pkin(7) * t325 + t232;
t196 = -mrSges(3,2) * t235 + t303;
t195 = mrSges(3,1) * t235 - t304;
t180 = Ifges(4,4) * t185;
t175 = Ifges(3,1) * t325 + t231 + t375;
t174 = t374 + (t254 * Ifges(3,2) + t386) * t346;
t162 = -mrSges(4,2) * t235 + t378;
t152 = Ifges(5,4) * t157;
t149 = t214 * t255 - t251 * t260;
t146 = t251 * t214 + t255 * t260;
t133 = -mrSges(4,1) * t185 - mrSges(4,2) * t186;
t132 = -t148 * t253 + t249 * t361;
t131 = -t148 * t249 - t253 * t361;
t123 = mrSges(4,1) * t234 - mrSges(4,3) * t140;
t122 = -mrSges(4,2) * t234 + mrSges(4,3) * t139;
t117 = -t186 * Ifges(4,1) + t235 * Ifges(4,5) + t180;
t116 = t185 * Ifges(4,2) + t235 * Ifges(4,6) - t376;
t112 = -qJD(4) * t170 + t188 * t253;
t111 = qJD(4) * t171 + t188 * t249;
t99 = -mrSges(5,2) * t181 + t388;
t86 = -pkin(4) * t275 - pkin(9) * t157;
t73 = t132 * t252 - t149 * t248;
t72 = -t132 * t248 - t149 * t252;
t70 = -Ifges(5,1) * t275 + t181 * Ifges(5,5) + t152;
t68 = -Ifges(5,5) * t275 + t157 * Ifges(5,6) + t181 * Ifges(5,3);
t67 = mrSges(6,1) * t156 - mrSges(6,3) * t95;
t66 = -mrSges(6,2) * t156 + mrSges(6,3) * t94;
t55 = pkin(4) * t261 - t61;
t54 = qJD(5) * t120 + t112 * t252 - t187 * t248;
t53 = -qJD(5) * t121 - t112 * t248 - t187 * t252;
t52 = -mrSges(5,2) * t138 + mrSges(5,3) * t78;
t49 = pkin(4) * t186 - t58;
t35 = pkin(4) * t111 - pkin(9) * t112 + t84;
t34 = -mrSges(5,1) * t78 + mrSges(5,2) * t77;
t29 = t248 * t86 + t252 * t47;
t28 = -t248 * t47 + t252 * t86;
t26 = t77 * Ifges(5,4) + t78 * Ifges(5,2) + t138 * Ifges(5,6);
t19 = pkin(4) * t187 - t25;
t18 = -pkin(9) * t187 + t24;
t4 = -qJD(5) * t21 - t18 * t248 + t252 * t35;
t3 = qJD(5) * t20 + t18 * t252 + t248 * t35;
t7 = [(-mrSges(5,3) * t12 + 0.2e1 * t426) * t171 + (-m(4) * t309 + t148 * mrSges(4,1) + t251 * mrSges(2,2) - m(6) * (pkin(4) * t132 + t295) - t73 * mrSges(6,1) - t72 * mrSges(6,2) - m(5) * t295 - t132 * mrSges(5,1) - t208 * mrSges(3,1) - t207 * mrSges(3,2) + t332 * t255 + t433 * t149 + t307 * t131) * g(2) + (Ifges(4,1) * t188 + Ifges(4,4) * t187) * t401 + (Ifges(5,1) * t112 - Ifges(5,5) * t187) * t404 + t204 * (-mrSges(4,1) * t187 + mrSges(4,2) * t188) + (t187 * t91 - t188 * t90) * mrSges(4,3) + t235 * (Ifges(4,5) * t188 + Ifges(4,6) * t187) / 0.2e1 + t181 * (Ifges(5,5) * t112 - Ifges(5,3) * t187) / 0.2e1 + t157 * (Ifges(5,4) * t112 - Ifges(5,6) * t187) / 0.2e1 + t185 * (Ifges(4,4) * t188 + Ifges(4,2) * t187) / 0.2e1 - (t172 * mrSges(4,1) - t46 * mrSges(4,3) - Ifges(4,4) * t140 - Ifges(4,2) * t139 - Ifges(4,6) * t234 + t432) * t261 + (-Ifges(5,6) * t409 - t11 * mrSges(5,3) - t26 / 0.2e1 + t5 / 0.2e1 - Ifges(5,2) * t415 - Ifges(5,4) * t416 + Ifges(6,3) * t417 + Ifges(6,6) * t424 + Ifges(6,5) * t425 + t452) * t170 + t112 * t458 + (-t206 * mrSges(3,1) + t454 * mrSges(3,2) + (t203 * t337 + mrSges(2,2)) * t255 + (t244 * t337 - t332) * t251 + t447 * t128 + (-t290 - t433) * t146 + t307 * t127 + (-pkin(3) * t410 + mrSges(4,1)) * t147) * g(1) + ((mrSges(3,1) * t201 - mrSges(3,2) * t202 + (m(3) * t398 - t209) * qJDD(1)) * pkin(1) + (mrSges(3,3) * t153 + Ifges(3,4) * t202 + Ifges(3,2) * t201 + Ifges(3,6) * t234) * t254 + (-mrSges(3,3) * t154 + Ifges(3,1) * t202 + Ifges(3,4) * t201 + Ifges(3,5) * t234) * t250 + ((t375 / 0.2e1 + t175 / 0.2e1 - t197 * mrSges(3,3)) * t254 + (-t374 / 0.2e1 - t174 / 0.2e1 - t198 * mrSges(3,3) + (m(4) * t204 + t133) * pkin(2)) * t250 + (t254 * (Ifges(3,4) * t254 - Ifges(3,2) * t250) / 0.2e1 - t435) * t346) * qJD(2) + (g(1) * t255 + g(2) * t251) * (-m(3) * pkin(7) - mrSges(3,3) - mrSges(4,3))) * t246 - t9 * t292 - (mrSges(4,2) * t172 - mrSges(4,3) * t45 + Ifges(4,1) * t140 + Ifges(4,4) * t139 + Ifges(4,5) * t234) * t193 + t199 * t196 - t200 * t195 + m(3) * (t153 * t347 + t154 * t210 - t197 * t200 + t198 * t199) + m(4) * (t113 * t45 + t114 * t46 + t172 * t215 - t84 * t90 + t85 * t91) + m(5) * (t102 * t42 + t11 * t369 + t12 * t61 + t24 * t48 + t25 * t47 + t81 * t84) + t369 * t52 + t347 * (-mrSges(3,2) * t234 + mrSges(3,3) * t201) + (Ifges(6,1) * t54 + Ifges(6,4) * t53) * t411 + (Ifges(6,1) * t121 + Ifges(6,4) * t120) * t425 + (Ifges(6,5) * t54 + Ifges(6,6) * t53) * t407 + (Ifges(6,5) * t121 + Ifges(6,6) * t120) * t417 + (Ifges(6,4) * t54 + Ifges(6,2) * t53) * t413 + (Ifges(6,4) * t121 + Ifges(6,2) * t120) * t424 + t210 * (mrSges(3,1) * t234 - mrSges(3,3) * t202) + t42 * t294 - t187 * t450 - t112 * t47 * mrSges(5,3) + m(6) * (t1 * t21 + t14 * t4 + t15 * t3 + t19 * t39 + t2 * t20 + t55 * t9) + (t222 / 0.2e1 + t223 / 0.2e1 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t234 + t430) * t247 + t187 * t116 / 0.2e1 + t188 * t117 / 0.2e1 - t187 * t68 / 0.2e1 + t85 * t162 + Ifges(2,3) * qJDD(1) + t114 * t122 + t113 * t123 + t112 * t70 / 0.2e1 + t24 * t99 + t25 * t100 + t102 * t34 + t3 * t66 + t4 * t67 + t61 * t51 + t55 * t10 + t19 * t57 + t39 * (-mrSges(6,1) * t53 + mrSges(6,2) * t54) + t54 * t419 + t53 * t421 + t187 * t449 + (-Ifges(5,4) * t404 + Ifges(6,3) * t407 + Ifges(6,5) * t411 + Ifges(6,6) * t413 - t48 * mrSges(5,3) - t456 / 0.2e1 - t457 / 0.2e1 + t451 + t453) * t111 + t121 * t427 + t120 * t428 + t215 * t310 + (t1 * t120 - t121 * t2 - t14 * t54 + t15 * t53) * mrSges(6,3) + t370 * t84 + t20 * t16 + t21 * t17; (t429 - t453) * t365 - ((-Ifges(3,2) * t325 + t175 + t231) * t254 + t235 * (Ifges(3,5) * t254 - Ifges(3,6) * t250)) * t346 / 0.2e1 + (g(1) * t148 - g(2) * t147 + g(3) * t193 - t440 * t48 + (-t341 + t364) * t47 + t438) * mrSges(5,3) + (t157 * t286 + t181 * t283 - t275 * t289) * qJD(4) / 0.2e1 + (-t106 * t81 + t243 * t42 - t47 * t58 - t48 * t59) * m(5) + (t1 * t166 + t443 * t14 + t444 * t15 + t165 * t2 - t39 * t49) * m(6) + (((-t249 * t48 - t253 * t47) * qJD(4) + t438) * m(5) + t448 * t249 + (t341 * t39 + t393) * m(6) + t253 * t52) * t242 + (-m(4) * t236 - m(5) * t298 - m(6) * t348 - t193 * t446 + t261 * t434 + t209) * g(3) + t116 * t401 + (-Ifges(5,3) * t186 + t185 * t283) * t403 + (-Ifges(5,5) * t186 + t185 * t289) * t405 + (-Ifges(5,6) * t186 + t185 * t286) * t406 + (-t281 * t339 + (Ifges(6,3) * t249 + t253 * t282) * qJD(4)) * t407 + (Ifges(5,5) * t249 + Ifges(5,6) * t253) * t409 + (-t287 * t339 + (Ifges(6,5) * t249 + t253 * t288) * qJD(4)) * t411 + (-t284 * t339 + (Ifges(6,6) * t249 + t253 * t285) * qJD(4)) * t413 + t290 * t393 + t122 * t396 + t222 + t223 + t461 * t99 + (-m(4) * t299 - mrSges(3,1) * t454 - mrSges(3,2) * t206 + t410 * (t146 * pkin(3) + t299) - t431 * t147 + t434 * t146) * g(2) + (t303 - t196) * t197 + (t314 * t37 + t315 * t38) * t252 + t90 * t378 + (Ifges(6,5) * t126 + Ifges(6,6) * t125) * t408 + (t320 - t49) * t57 + (mrSges(6,1) * t440 + mrSges(6,3) * t441) * t14 + (mrSges(6,1) * t442 - mrSges(6,2) * t441) * t39 + (-mrSges(6,2) * t440 - mrSges(6,3) * t442) * t15 + t123 * t326 + t435 * qJD(1) ^ 2 * t246 ^ 2 - t436 * t81 * (mrSges(5,1) * t249 + mrSges(5,2) * t253) - t204 * (-mrSges(4,1) * t186 + mrSges(4,2) * t185) + t430 + (Ifges(4,1) * t185 + t376 + t68) * t186 / 0.2e1 - (Ifges(4,2) * t186 + t117 + t180) * t185 / 0.2e1 - t253 * t5 / 0.2e1 + t253 * t26 / 0.2e1 + (Ifges(6,4) * t126 + Ifges(6,2) * t125) * t414 + ((t245 * t46 + t368 * t45) * pkin(2) + t106 * t90 - t107 * t91 - t204 * t306) * m(4) + t243 * t34 - t235 * (Ifges(4,5) * t185 + Ifges(4,6) * t186) / 0.2e1 - t42 * t293 + t186 * t450 + (-t58 - t320) * t100 + (t315 - t364 / 0.2e1) * t70 + t248 * t38 * t314 + t165 * t16 + t166 * t17 - t107 * t162 + (Ifges(6,1) * t126 + Ifges(6,4) * t125) * t412 + (Ifges(5,2) * t253 + t385) * t415 + (Ifges(5,1) * t249 + t384) * t416 + (-Ifges(6,3) * t253 + t249 * t282) * t417 + t126 * t420 - t370 * t106 + t443 * t67 + t444 * t66 + t453 * t342 + t455 * t422 - t186 * t449 + (-m(4) * t264 - mrSges(3,1) * t207 + mrSges(3,2) * t208 + t410 * (t149 * pkin(3) + t264) + t431 * t148 + t434 * t149) * g(1) + (t304 + t195) * t198 + (-Ifges(6,6) * t253 + t249 * t285) * t424 + (-Ifges(6,5) * t253 + t249 * t288) * t425 + t249 * t426 + t355 * t427 - t133 * t306 + t174 * t325 / 0.2e1 + t2 * (-mrSges(6,1) * t253 - mrSges(6,3) * t355) - t6 * t357 / 0.2e1 + t1 * (mrSges(6,2) * t253 - mrSges(6,3) * t357) - t91 * t377; -t125 * t67 - t126 * t66 - t185 * t162 + t370 * t186 + (-t185 * t99 + (-t248 * t67 + t252 * t66 + t99) * qJD(4) - t448) * t253 + (t52 + (-t248 * t66 - t252 * t67) * qJD(5) + t436 * t445 + t439) * t249 + t310 + (-t247 * g(3) + (-g(1) * t251 + g(2) * t255) * t246) * t337 + (-t125 * t14 - t126 * t15 - t39 * t365 + (-t9 + (-t14 * t248 + t15 * t252) * qJD(4)) * t253 + (qJD(4) * t39 + t257) * t249) * m(6) + (t11 * t249 + t12 * t253 + t186 * t81 - t436 * (-t249 * t47 + t253 * t48)) * m(5) + (-t185 * t91 - t186 * t90 + t172) * m(4); -(-Ifges(5,2) * t406 - Ifges(5,6) * t403 + t429 - t451) * t275 + t69 * t404 + qJD(5) * t272 + t432 + t9 * t291 + (-pkin(4) * t9 - t14 * t28 - t15 * t29) * m(6) + (t388 - t99) * t47 + ((-t340 + t367) * t15 + (-t338 + t366) * t14 + t437) * mrSges(6,3) + (m(6) * t257 - t338 * t67 - t340 * t66 + t439) * pkin(9) + (t36 + t379) * t405 + (t156 * t282 + t285 * t94 + t288 * t95) * qJD(5) / 0.2e1 + (t70 + t152) * t406 - t29 * t66 - t28 * t67 + t281 * t417 + t338 * t419 + t366 * t420 + t367 * t421 + t340 * t422 + (-m(6) * t39 - t387 + t445) * t48 + (t131 * t447 + t132 * t307) * g(1) + (-t127 * t447 + t128 * t307) * g(2) + (Ifges(5,1) * t405 + Ifges(5,5) * t403 + t282 * t408 + t285 * t414 + t288 * t412 - t272 - t458) * t157 + t284 * t424 + t287 * t425 + t248 * t427 + t252 * t428 + (t170 * t266 - t171 * t331 + t294) * g(3) - pkin(4) * t10; -t39 * (mrSges(6,1) * t95 + mrSges(6,2) * t94) + (Ifges(6,1) * t94 - t399) * t412 + t37 * t411 + (Ifges(6,5) * t94 - Ifges(6,6) * t95) * t408 - t14 * t66 + t15 * t67 - g(1) * (mrSges(6,1) * t72 - mrSges(6,2) * t73) - g(2) * ((-t128 * t248 - t146 * t252) * mrSges(6,1) + (-t128 * t252 + t146 * t248) * mrSges(6,2)) - g(3) * t292 + (t14 * t94 + t15 * t95) * mrSges(6,3) + t5 + (-Ifges(6,2) * t95 + t38 + t92) * t414 + t452;];
tau = t7;
