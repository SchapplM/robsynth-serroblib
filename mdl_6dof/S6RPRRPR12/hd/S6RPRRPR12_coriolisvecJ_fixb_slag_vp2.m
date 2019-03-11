% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:10
% EndTime: 2019-03-09 05:48:57
% DurationCPUTime: 23.73s
% Computational Cost: add. (19375->713), mult. (64249->939), div. (0->0), fcn. (54284->12), ass. (0->340)
t244 = sin(pkin(12));
t246 = sin(pkin(6));
t247 = cos(pkin(7));
t250 = sin(qJ(3));
t253 = cos(qJ(3));
t354 = cos(pkin(12));
t311 = t354 * t250;
t245 = sin(pkin(7));
t355 = cos(pkin(6));
t312 = t355 * t245;
t197 = t250 * t312 + (t244 * t253 + t247 * t311) * t246;
t193 = t197 * qJD(1);
t249 = sin(qJ(4));
t252 = cos(qJ(4));
t313 = t246 * t354;
t218 = -t245 * t313 + t247 * t355;
t262 = qJD(1) * t218 + qJD(3);
t259 = t252 * t262;
t162 = t193 * t249 - t259;
t400 = -t162 / 0.2e1;
t341 = t252 * t193;
t163 = t249 * t262 + t341;
t161 = qJD(6) + t163;
t401 = t161 / 0.2e1;
t338 = qJD(1) * t246;
t324 = t244 * t338;
t310 = t354 * t253;
t303 = t247 * t310;
t304 = t253 * t312;
t478 = t246 * t303 + t304;
t190 = qJD(1) * t478 - t250 * t324;
t187 = qJD(4) - t190;
t248 = sin(qJ(6));
t251 = cos(qJ(6));
t129 = t162 * t248 + t187 * t251;
t404 = t129 / 0.2e1;
t128 = t162 * t251 - t187 * t248;
t406 = t128 / 0.2e1;
t491 = Ifges(5,4) * t400 + Ifges(7,5) * t404 + Ifges(7,6) * t406 + Ifges(7,3) * t401;
t494 = Ifges(6,6) * t400 + t491;
t397 = t163 / 0.2e1;
t493 = Ifges(5,1) * t397;
t393 = t187 / 0.2e1;
t394 = -t187 / 0.2e1;
t325 = pkin(1) * t355;
t241 = t354 * t325;
t237 = qJD(1) * t241;
t349 = t244 * t246;
t267 = t355 * pkin(2) + (-pkin(9) * t247 - qJ(2)) * t349;
t188 = qJD(1) * t267 + t237;
t212 = (-pkin(9) * t244 * t245 - pkin(2) * t354 - pkin(1)) * t246;
t204 = qJD(1) * t212 + qJD(2);
t158 = -t188 * t245 + t204 * t247;
t268 = -pkin(3) * t190 - pkin(10) * t193 + t158;
t240 = t244 * t325;
t302 = qJD(1) * t313;
t216 = qJ(2) * t302 + qJD(1) * t240;
t266 = (t247 * t313 + t312) * pkin(9);
t182 = qJD(1) * t266 + t216;
t172 = t253 * t182;
t345 = t247 * t250;
t347 = t245 * t250;
t119 = t188 * t345 + t204 * t347 + t172;
t99 = pkin(10) * t262 + t119;
t45 = t249 * t99 - t252 * t268;
t454 = -qJD(5) - t45;
t42 = -pkin(4) * t187 - t454;
t413 = pkin(4) + pkin(11);
t279 = pkin(5) * t163 + t45;
t476 = qJD(5) + t279;
t31 = -t187 * t413 + t476;
t171 = t250 * t182;
t344 = t247 * t253;
t346 = t245 * t253;
t98 = -pkin(3) * t262 - t188 * t344 - t204 * t346 + t171;
t54 = t162 * pkin(4) - t163 * qJ(5) + t98;
t41 = t162 * pkin(11) + t54;
t5 = -t248 * t41 + t251 * t31;
t6 = t248 * t31 + t251 * t41;
t488 = t42 * mrSges(6,1) + t5 * mrSges(7,1) - t6 * mrSges(7,2) + t45 * mrSges(5,3) + Ifges(6,4) * t394 + Ifges(5,5) * t393 + Ifges(6,2) * t397 + t493 + t494;
t492 = t98 * mrSges(5,2) - t54 * mrSges(6,3) + t488;
t270 = t246 * (-t244 * t345 + t310);
t210 = qJD(1) * t270;
t479 = -qJD(3) * t346 + t210;
t490 = Ifges(6,6) + Ifges(5,4);
t333 = qJD(4) * t249;
t352 = t190 * t249;
t489 = -qJD(5) * t249 - t119 + (t333 - t352) * pkin(4);
t348 = t244 * t250;
t191 = (t304 + (t303 - t348) * t246) * qJD(3);
t258 = qJD(1) * t191;
t126 = -qJD(4) * t259 + t193 * t333 - t252 * t258;
t411 = -t126 / 0.2e1;
t350 = t191 * t249;
t127 = (qJD(3) * t249 + t341) * qJD(4) + (t218 * t333 + t350) * qJD(1);
t409 = -t127 / 0.2e1;
t192 = t197 * qJD(3);
t180 = qJD(1) * t192;
t395 = t180 / 0.2e1;
t486 = Ifges(6,4) - Ifges(5,5);
t470 = Ifges(5,6) - Ifges(6,5);
t469 = Ifges(5,3) + Ifges(6,1);
t485 = t489 + t187 * (pkin(11) * t249 - qJ(5) * t252);
t283 = t188 * t247 + t204 * t245;
t118 = t283 * t253 - t171;
t114 = t249 * t118;
t151 = pkin(3) * t193 - pkin(10) * t190;
t412 = pkin(5) + pkin(10);
t235 = t412 * t252;
t484 = qJD(4) * t235 - t114 - (pkin(5) * t190 - t151) * t252 + t413 * t193;
t222 = t247 * t249 + t252 * t347;
t307 = t245 * t324;
t483 = qJD(4) * t222 - t249 * t479 + t252 * t307;
t269 = t246 * (t244 * t344 + t311);
t209 = qJD(1) * t269;
t336 = qJD(3) * t250;
t322 = t245 * t336;
t482 = t322 - t209;
t46 = t249 * t268 + t252 * t99;
t43 = -t187 * qJ(5) - t46;
t441 = t43 * mrSges(6,1) - t46 * mrSges(5,3);
t481 = -t98 * mrSges(5,1) + t54 * mrSges(6,2) - t441;
t410 = t126 / 0.2e1;
t408 = t127 / 0.2e1;
t265 = qJD(2) * t270;
t92 = qJD(1) * t265 + qJD(3) * t118;
t480 = qJD(4) * t268 + t92;
t339 = qJ(2) * t313 + t240;
t195 = t266 + t339;
t198 = t241 + t267;
t282 = t198 * t247 + t212 * t245;
t134 = -t250 * t195 + t282 * t253;
t57 = qJD(6) * t128 + t127 * t248 + t180 * t251;
t58 = -qJD(6) * t129 + t127 * t251 - t180 * t248;
t13 = Ifges(7,5) * t57 + Ifges(7,6) * t58 - Ifges(7,3) * t126;
t264 = qJD(2) * t269;
t93 = qJD(1) * t264 + (t250 * t283 + t172) * qJD(3);
t255 = t126 * qJ(5) - t163 * qJD(5) + t93;
t32 = t127 * pkin(4) + t255;
t419 = t58 / 0.2e1;
t420 = t57 / 0.2e1;
t25 = t127 * t413 + t255;
t337 = qJD(2) * t246;
t323 = t244 * t337;
t305 = qJD(1) * t323;
t281 = t245 * t305;
t141 = t180 * pkin(3) - pkin(10) * t258 + t281;
t332 = qJD(4) * t252;
t24 = t141 * t252 - t249 * t480 - t99 * t332;
t7 = -pkin(5) * t126 - t180 * t413 - t24;
t1 = qJD(6) * t5 + t248 * t7 + t25 * t251;
t2 = -qJD(6) * t6 - t248 * t25 + t251 * t7;
t437 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t472 = -Ifges(6,4) / 0.2e1;
t475 = t13 / 0.2e1 + mrSges(5,2) * t93 - mrSges(6,3) * t32 + Ifges(7,5) * t420 + Ifges(7,6) * t419 + t180 * t472 + t437 - Ifges(6,2) * t410 + (0.2e1 * Ifges(5,1) + Ifges(7,3) + Ifges(6,2)) * t411 + t490 * t409 + (Ifges(5,5) - t486) * t395;
t473 = 0.2e1 * Ifges(6,3) * t408 + t93 * mrSges(5,1) - t180 * Ifges(5,6) / 0.2e1 - t32 * mrSges(6,2) + t490 * t410 + (t408 - t409) * Ifges(5,2) + (-t470 + Ifges(6,5)) * t395;
t467 = -t162 * t470 - t163 * t486 + t187 * t469;
t88 = mrSges(6,1) * t127 - mrSges(6,3) * t180;
t91 = -mrSges(5,2) * t180 - mrSges(5,3) * t127;
t465 = -t88 + t91;
t89 = -t126 * mrSges(6,1) + t180 * mrSges(6,2);
t90 = mrSges(5,1) * t180 + mrSges(5,3) * t126;
t464 = -t90 + t89;
t186 = Ifges(4,4) * t190;
t463 = Ifges(4,5) * t258;
t462 = Ifges(4,2) * t190;
t461 = Ifges(4,6) * t180;
t460 = t262 * Ifges(4,5);
t459 = t262 * Ifges(4,6);
t314 = -qJ(5) * t249 - pkin(3);
t223 = -t252 * t413 + t314;
t234 = t412 * t249;
t200 = t223 * t251 + t234 * t248;
t458 = -qJD(6) * t200 - t248 * t485 + t251 * t484;
t199 = -t223 * t248 + t234 * t251;
t457 = qJD(6) * t199 + t248 * t484 + t251 * t485;
t351 = t190 * t252;
t456 = (-t332 + t351) * qJ(5) + t489;
t75 = t252 * t118 + t249 * t151;
t65 = -qJ(5) * t193 - t75;
t455 = pkin(5) * t352 - t412 * t333 + t65;
t362 = t193 * mrSges(4,3);
t453 = -mrSges(4,1) * t262 + mrSges(5,1) * t162 + mrSges(5,2) * t163 + t362;
t130 = mrSges(6,1) * t162 - mrSges(6,3) * t187;
t132 = -mrSges(5,2) * t187 - mrSges(5,3) * t162;
t452 = -t130 + t132;
t131 = mrSges(6,1) * t163 + mrSges(6,2) * t187;
t133 = mrSges(5,1) * t187 - mrSges(5,3) * t163;
t451 = -t131 + t133;
t147 = -t193 * t248 + t251 * t352;
t331 = qJD(6) * t252;
t320 = t248 * t331;
t450 = -t251 * t333 + t147 - t320;
t148 = t193 * t251 + t248 * t352;
t449 = -t248 * t333 + t251 * t331 + t148;
t328 = t249 * t347;
t221 = -t252 * t247 + t328;
t278 = -t221 * t248 + t251 * t346;
t448 = qJD(6) * t278 - t248 * t482 + t251 * t483;
t205 = t221 * t251 + t248 * t346;
t447 = qJD(6) * t205 + t248 * t483 + t251 * t482;
t446 = t98 * (mrSges(5,1) * t249 + mrSges(5,2) * t252) + t54 * (-mrSges(6,2) * t249 - mrSges(6,3) * t252);
t445 = -t249 * t470 - t252 * t486;
t444 = t126 * t486 - t127 * t470 + t180 * t469;
t23 = t249 * t141 + t252 * t480 - t333 * t99;
t443 = t23 * t252 - t24 * t249;
t16 = -qJ(5) * t180 - qJD(5) * t187 - t23;
t17 = -pkin(4) * t180 - t24;
t442 = -t16 * t252 + t17 * t249;
t440 = t1 * t248 + t2 * t251;
t80 = t187 * Ifges(6,5) - t163 * Ifges(6,6) + t162 * Ifges(6,3);
t83 = t163 * Ifges(5,4) - t162 * Ifges(5,2) + t187 * Ifges(5,6);
t439 = t83 / 0.2e1 - t80 / 0.2e1;
t438 = -m(5) * t98 - t453;
t435 = -t93 * mrSges(4,1) - t92 * mrSges(4,2);
t432 = t158 * mrSges(4,2) + t460 / 0.2e1;
t430 = -m(4) * t118 - t438;
t429 = t24 * mrSges(5,1) - t23 * mrSges(5,2) + t17 * mrSges(6,2) - t16 * mrSges(6,3);
t289 = Ifges(7,5) * t248 + Ifges(7,6) * t251;
t372 = Ifges(7,4) * t248;
t292 = Ifges(7,2) * t251 + t372;
t371 = Ifges(7,4) * t251;
t296 = Ifges(7,1) * t248 + t371;
t299 = mrSges(7,1) * t251 - mrSges(7,2) * t248;
t300 = t248 * t5 - t251 * t6;
t383 = t162 * pkin(5);
t33 = -t43 - t383;
t373 = Ifges(7,4) * t129;
t51 = t128 * Ifges(7,2) + t161 * Ifges(7,6) + t373;
t358 = t251 * t51;
t125 = Ifges(7,4) * t128;
t52 = t129 * Ifges(7,1) + t161 * Ifges(7,5) + t125;
t359 = t248 * t52;
t402 = -t161 / 0.2e1;
t405 = -t129 / 0.2e1;
t407 = -t128 / 0.2e1;
t426 = t289 * t402 + t292 * t407 + t296 * t405 + t33 * t299 + t300 * mrSges(7,3) - t359 / 0.2e1 - t358 / 0.2e1;
t424 = t43 * mrSges(6,3) + t45 * mrSges(5,1) + t46 * mrSges(5,2) - t158 * mrSges(4,1) + t459 / 0.2e1 - t42 * mrSges(6,2);
t15 = Ifges(7,1) * t57 + Ifges(7,4) * t58 - Ifges(7,5) * t126;
t422 = t15 / 0.2e1;
t421 = t51 / 0.2e1;
t415 = -t83 / 0.2e1;
t399 = t162 / 0.2e1;
t398 = -t163 / 0.2e1;
t392 = -t190 / 0.2e1;
t389 = t193 / 0.2e1;
t377 = mrSges(4,3) * t180;
t376 = Ifges(4,4) * t193;
t375 = Ifges(5,4) * t249;
t374 = Ifges(5,4) * t252;
t370 = Ifges(6,6) * t249;
t369 = Ifges(6,6) * t252;
t363 = t190 * mrSges(4,3);
t357 = t253 * t93;
t73 = -mrSges(7,1) * t128 + mrSges(7,2) * t129;
t356 = t73 - t130;
t353 = qJ(5) * t162;
t343 = t248 * t252;
t342 = t251 * t252;
t185 = t253 * t195;
t164 = -t198 * t245 + t247 * t212;
t196 = t246 * t348 - t478;
t110 = pkin(3) * t196 - pkin(10) * t197 + t164;
t135 = t198 * t345 + t212 * t347 + t185;
t117 = pkin(10) * t218 + t135;
t64 = t249 * t110 + t252 * t117;
t340 = -t461 + t463;
t330 = pkin(10) * t333;
t329 = pkin(10) * t332;
t63 = t110 * t252 - t249 * t117;
t74 = t151 * t252 - t114;
t306 = t245 * t323;
t48 = -qJ(5) * t196 - t64;
t298 = Ifges(5,1) * t252 - t375;
t297 = Ifges(7,1) * t251 - t372;
t295 = -Ifges(5,2) * t249 + t374;
t293 = -Ifges(7,2) * t248 + t371;
t290 = Ifges(7,5) * t251 - Ifges(7,6) * t248;
t288 = -Ifges(6,2) * t252 + t370;
t287 = Ifges(6,3) * t249 - t369;
t168 = t197 * t252 + t218 * t249;
t35 = pkin(5) * t168 - t196 * t413 - t63;
t167 = t197 * t249 - t218 * t252;
t116 = -pkin(3) * t218 - t134;
t260 = -qJ(5) * t168 + t116;
t47 = t167 * t413 + t260;
t9 = -t248 * t47 + t251 * t35;
t10 = t248 * t35 + t251 * t47;
t36 = -mrSges(7,1) * t126 - mrSges(7,3) * t57;
t37 = mrSges(7,2) * t126 + mrSges(7,3) * t58;
t285 = t248 * t37 + t251 * t36;
t78 = -mrSges(7,2) * t161 + mrSges(7,3) * t128;
t79 = mrSges(7,1) * t161 - mrSges(7,3) * t129;
t284 = -t248 * t79 + t251 * t78;
t142 = t167 * t251 - t196 * t248;
t143 = t167 * t248 + t196 * t251;
t102 = qJD(3) * t134 + t265;
t145 = pkin(3) * t192 - pkin(10) * t191 + t306;
t30 = -t249 * t102 - t110 * t333 - t117 * t332 + t145 * t252;
t275 = -(-qJ(2) * t324 + t237) * t244 + t216 * t354;
t29 = t252 * t102 + t110 * t332 - t117 * t333 + t249 * t145;
t219 = (mrSges(3,1) * t355 - mrSges(3,3) * t349) * qJD(1);
t220 = (-mrSges(3,2) * t355 + mrSges(3,3) * t313) * qJD(1);
t22 = -qJ(5) * t192 - qJD(5) * t196 - t29;
t257 = mrSges(4,3) * t258;
t103 = t264 + (t250 * t282 + t185) * qJD(3);
t140 = -t197 * t333 + (qJD(4) * t218 + t191) * t252;
t256 = -t140 * qJ(5) - t168 * qJD(5) + t103;
t230 = -pkin(4) * t252 + t314;
t165 = -mrSges(4,2) * t262 + t363;
t150 = -mrSges(4,1) * t190 + mrSges(4,2) * t193;
t146 = t180 * mrSges(4,1) + mrSges(4,2) * t258;
t139 = qJD(4) * t168 + t350;
t138 = Ifges(4,1) * t193 + t186 + t460;
t137 = t376 + t459 + t462;
t108 = -mrSges(6,2) * t162 - mrSges(6,3) * t163;
t106 = pkin(4) * t163 + t353;
t77 = t163 * t413 + t353;
t72 = mrSges(5,1) * t127 - mrSges(5,2) * t126;
t71 = -mrSges(6,2) * t127 + mrSges(6,3) * t126;
t70 = qJD(6) * t142 + t139 * t248 + t192 * t251;
t69 = -qJD(6) * t143 + t139 * t251 - t192 * t248;
t68 = pkin(4) * t167 + t260;
t67 = -pkin(4) * t193 - t74;
t49 = -pkin(4) * t196 - t63;
t40 = -pkin(5) * t167 - t48;
t39 = t46 - t383;
t34 = t139 * pkin(4) + t256;
t28 = t139 * t413 + t256;
t27 = -pkin(4) * t192 - t30;
t26 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t19 = t248 * t39 + t251 * t77;
t18 = -t248 * t77 + t251 * t39;
t14 = Ifges(7,4) * t57 + Ifges(7,2) * t58 - Ifges(7,6) * t126;
t12 = -pkin(5) * t139 - t22;
t11 = pkin(5) * t140 - t192 * t413 - t30;
t8 = -pkin(5) * t127 - t16;
t4 = -qJD(6) * t10 + t11 * t251 - t248 * t28;
t3 = qJD(6) * t9 + t11 * t248 + t251 * t28;
t20 = [(t16 * mrSges(6,1) - t23 * mrSges(5,3) - Ifges(5,4) * t411 + Ifges(6,6) * t410 + t473) * t167 + (mrSges(4,2) * t281 + mrSges(4,3) * t93 + Ifges(4,1) * t258 - Ifges(4,4) * t180) * t197 + m(3) * ((t354 * t339 + (qJ(2) * t349 - t241) * t244) * qJD(1) + t275) * t337 + (t340 / 0.2e1 + t463 / 0.2e1 - t461 / 0.2e1 + t435) * t218 + (-Ifges(6,2) * t398 - Ifges(6,6) * t399 - t393 * t486 + t491 + t492 + t493) * t140 + t430 * t103 + t150 * t306 + (mrSges(6,1) * t17 - mrSges(5,3) * t24 + Ifges(5,4) * t409 - Ifges(6,6) * t408 + t475) * t168 - t135 * t377 + (Ifges(7,4) * t143 + Ifges(7,2) * t142) * t419 + (Ifges(7,4) * t70 + Ifges(7,2) * t69) * t406 - 0.2e1 * t219 * t323 + (t1 * t142 - t143 * t2 - t5 * t70 + t6 * t69) * mrSges(7,3) + (Ifges(7,1) * t143 + Ifges(7,4) * t142) * t420 + (Ifges(7,1) * t70 + Ifges(7,4) * t69) * t404 + 0.2e1 * qJD(2) * t220 * t313 + m(7) * (t1 * t10 + t12 * t33 + t2 * t9 + t3 * t6 + t4 * t5 + t40 * t8) + m(6) * (t16 * t48 + t17 * t49 + t22 * t43 + t27 * t42 + t32 * t68 + t34 * t54) + m(4) * (t102 * t119 - t134 * t93 + t135 * t92 + (qJD(1) * t164 + t158) * t306) + m(5) * (t116 * t93 + t23 * t64 + t24 * t63 + t29 * t46 - t30 * t45) + t69 * t421 + t143 * t422 + (t467 / 0.2e1 - t137 / 0.2e1 - t462 / 0.2e1 - Ifges(4,4) * t389 + Ifges(5,5) * t397 + Ifges(6,4) * t398 + Ifges(6,5) * t399 + Ifges(5,6) * t400 + t469 * t393 - t424 - t119 * mrSges(4,3)) * t192 + (t415 + t80 / 0.2e1 - Ifges(5,4) * t397 + Ifges(6,6) * t398 + Ifges(6,3) * t399 - Ifges(5,2) * t400 - t470 * t393 - t481) * t139 + (t138 / 0.2e1 + t186 / 0.2e1 + Ifges(4,1) * t389 + t432 - t118 * mrSges(4,3)) * t191 + (mrSges(4,1) * t281 - mrSges(4,3) * t92 - Ifges(4,4) * t258 + Ifges(6,4) * t410 + Ifges(5,5) * t411 + Ifges(6,5) * t408 + Ifges(4,2) * t180 + Ifges(5,6) * t409 + t469 * t395 + t429 + t444 / 0.2e1) * t196 + t9 * t36 + t10 * t37 + t40 * t26 - t134 * t257 + t70 * t52 / 0.2e1 + t33 * (-mrSges(7,1) * t69 + mrSges(7,2) * t70) + t68 * t71 + t12 * t73 + t3 * t78 + t4 * t79 + t48 * t88 + t49 * t89 + t63 * t90 + t64 * t91 + (Ifges(7,5) * t143 + Ifges(7,6) * t142) * t411 + (Ifges(7,5) * t70 + Ifges(7,6) * t69) * t401 + t34 * t108 + t116 * t72 + t22 * t130 + t27 * t131 + t29 * t132 + t30 * t133 + t142 * t14 / 0.2e1 + t8 * (-mrSges(7,1) * t142 + mrSges(7,2) * t143) + t164 * t146 + t102 * t165; -m(3) * t275 * t338 + t247 * t146 - t150 * t307 + t205 * t36 - t278 * t37 + t219 * t324 - t220 * t302 - t347 * t377 + t448 * t79 + t447 * t78 + (m(5) * (t336 * t98 - t357) + m(6) * (-t253 * t32 + t336 * t54)) * t245 + (-m(5) * t24 + m(6) * t17 + t464) * t221 - t479 * t165 + (-t1 * t278 + t2 * t205 + t447 * t6 + t448 * t5) * m(7) + ((t247 * t305 + t250 * t92 - t357 + (-t118 * t250 + t119 * t253) * qJD(3)) * t245 - t119 * t210 - t158 * t307) * m(4) + (-t257 - t72 - t71) * t346 + (t108 + t453) * t322 + (m(5) * t23 - m(6) * t16 + m(7) * t8 + t26 + t465) * t222 + (-m(6) * t54 - t108 - t430) * t209 + t483 * (m(5) * t45 + m(6) * t42 - t451) + (qJD(4) * t328 - t247 * t332 + t249 * t307 + t252 * t479) * (-m(5) * t46 + m(6) * t43 - m(7) * t33 - t452 - t73); (t445 * t393 + (-t295 / 0.2e1 + t287 / 0.2e1) * t162 + (Ifges(7,3) * t252 + t249 * t289) * t401 + (Ifges(7,5) * t252 + t249 * t296) * t404 + (Ifges(7,6) * t252 + t249 * t292) * t406 + t446 + (t298 / 0.2e1 - t288 / 0.2e1) * t163) * qJD(4) + (-t289 * t411 - t292 * t419 - t296 * t420 + t8 * t299 - t473) * t252 + (-t352 * t43 + t442) * mrSges(6,1) + (t352 * t46 + t443) * mrSges(5,3) + (-pkin(3) * t93 + t45 * t74 - t46 * t75) * m(5) + (((t249 * t43 + t252 * t42) * qJD(4) + t442) * m(6) + t464 * t249 + t465 * t252 + ((-t249 * t46 + t252 * t45) * qJD(4) + t443) * m(5)) * pkin(10) + (t230 * t32 - t42 * t67 - t43 * t65 + t456 * t54) * m(6) + (Ifges(7,5) * t405 + Ifges(7,6) * t407 + Ifges(7,3) * t402 - t488) * t351 + t488 * t332 + (-t329 - t74) * t133 + t455 * t73 + t456 * t108 + t457 * t78 + (t1 * t200 + t199 * t2 + t235 * t8 + t33 * t455 + t457 * t6 + t458 * t5) * m(7) + t458 * t79 + (t186 + t138) * t392 + t439 * t352 + (t415 + t441) * t333 + (t80 + t358 + t359) * t333 / 0.2e1 - (qJD(6) * t52 + t14) * t342 / 0.2e1 + (-t290 * t401 - t293 * t406 - t297 * t404) * t331 - t15 * t343 / 0.2e1 - t369 * t410 + t374 * t411 + (t329 - t67) * t131 + t435 + t475 * t249 - (Ifges(4,1) * t190 - t376 + t467) * t193 / 0.2e1 + (Ifges(6,4) * t397 + Ifges(5,5) * t398 + Ifges(6,5) * t400 - Ifges(4,2) * t392 + Ifges(5,6) * t399 + t394 * t469 + t424) * t193 + (t287 * t400 + t288 * t397 + t295 * t399 + t298 * t398 + t394 * t445 - t432 - t446) * t190 + (t330 - t65) * t130 - t370 * t408 + t375 * t409 + (t362 + t438) * t119 + (-t165 + t363) * t118 + (-t330 - t75) * t132 + (mrSges(7,1) * t450 - mrSges(7,2) * t449) * t33 + (-t1 * t342 + t2 * t343 + t449 * t5 - t450 * t6) * mrSges(7,3) + (Ifges(7,5) * t148 + Ifges(7,6) * t147) * t402 + (Ifges(7,4) * t148 + Ifges(7,2) * t147) * t407 + t320 * t421 + (Ifges(7,1) * t148 + Ifges(7,4) * t147) * t405 + t235 * t26 + t230 * t71 + t199 * t36 + t200 * t37 + t340 - pkin(3) * t72 - t147 * t51 / 0.2e1 - t148 * t52 / 0.2e1 + t137 * t389; t356 * qJD(5) - ((-m(7) * t300 + t284) * qJD(6) + t285 + m(7) * t440) * t413 + ((Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t163 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t187 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t162 + t426 + t439 + t481) * t163 + t444 + (t8 * qJ(5) - t18 * t5 - t19 * t6 + t476 * t33) * m(7) - t440 * mrSges(7,3) + t429 + (t26 - t88) * qJ(5) + ((t472 + Ifges(5,5) / 0.2e1) * t187 + t492 + t494) * t162 + t426 * qJD(6) + (-pkin(4) * t17 - qJ(5) * t16 - t106 * t54 - t42 * t46 + t43 * t454) * m(6) + t452 * t45 + t451 * t46 + t251 * t422 + t290 * t411 + t293 * t419 + t297 * t420 + t8 * (mrSges(7,1) * t248 + mrSges(7,2) * t251) - t248 * t14 / 0.2e1 + t279 * t73 - t19 * t78 - t18 * t79 - pkin(4) * t89 - t106 * t108; -t356 * t187 + t284 * qJD(6) + (t108 + t284) * t163 + t285 + t89 + (-t161 * t300 - t187 * t33 + t440) * m(7) + (t163 * t54 + t187 * t43 + t17) * m(6); -t33 * (mrSges(7,1) * t129 + mrSges(7,2) * t128) + (Ifges(7,1) * t128 - t373) * t405 + t51 * t404 + (Ifges(7,5) * t128 - Ifges(7,6) * t129) * t402 - t5 * t78 + t6 * t79 + (t128 * t5 + t129 * t6) * mrSges(7,3) + t13 + (-Ifges(7,2) * t129 + t125 + t52) * t407 + t437;];
tauc  = t20(:);
