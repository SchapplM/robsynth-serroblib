% Calculate vector of inverse dynamics joint torques for
% S5RRRPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:26
% EndTime: 2019-12-31 21:33:06
% DurationCPUTime: 22.37s
% Computational Cost: add. (5331->700), mult. (11591->911), div. (0->0), fcn. (7380->8), ass. (0->331)
t489 = Ifges(4,1) + Ifges(5,1);
t488 = Ifges(5,4) + Ifges(4,5);
t486 = Ifges(4,6) - Ifges(5,6);
t224 = sin(qJ(2));
t216 = t224 * pkin(7);
t228 = cos(qJ(2));
t219 = t228 * pkin(2);
t321 = -pkin(1) - t219;
t256 = t321 - t216;
t153 = t256 * qJD(1);
t341 = qJD(1) * t228;
t212 = pkin(6) * t341;
t185 = qJD(2) * pkin(7) + t212;
t223 = sin(qJ(3));
t227 = cos(qJ(3));
t92 = t227 * t153 - t223 * t185;
t494 = -t92 + qJD(4);
t196 = qJD(3) - t341;
t191 = qJD(5) - t196;
t401 = t191 / 0.2e1;
t342 = qJD(1) * t224;
t319 = t223 * t342;
t332 = t227 * qJD(2);
t165 = t319 - t332;
t318 = t227 * t342;
t166 = qJD(2) * t223 + t318;
t222 = sin(qJ(5));
t226 = cos(qJ(5));
t259 = t165 * t222 + t166 * t226;
t411 = t259 / 0.2e1;
t88 = t165 * t226 - t166 * t222;
t413 = t88 / 0.2e1;
t493 = Ifges(6,5) * t411 + Ifges(6,6) * t413 + Ifges(6,3) * t401;
t405 = t165 / 0.2e1;
t492 = -t166 / 0.2e1;
t491 = -t196 / 0.2e1;
t490 = mrSges(5,2) + mrSges(4,3);
t487 = Ifges(5,2) + Ifges(4,3);
t288 = pkin(2) * t224 - pkin(7) * t228;
t171 = t288 * qJD(1);
t147 = t223 * t171;
t204 = qJ(4) * t342;
t337 = qJD(3) * t223;
t352 = t224 * t227;
t353 = t223 * t228;
t409 = pkin(7) - pkin(8);
t485 = t409 * t337 + t147 + t204 + (-pkin(6) * t352 + pkin(8) * t353) * qJD(1);
t187 = t409 * t227;
t320 = -pkin(6) * t223 - pkin(3);
t349 = t227 * t228;
t241 = -pkin(8) * t349 + (-pkin(4) + t320) * t224;
t356 = t171 * t227;
t484 = -qJD(1) * t241 + qJD(3) * t187 + t356;
t483 = -pkin(8) * t166 + t494;
t449 = -t486 * t223 + t488 * t227;
t376 = Ifges(5,5) * t223;
t378 = Ifges(4,4) * t223;
t446 = t489 * t227 + t376 - t378;
t482 = m(6) * pkin(8) + mrSges(6,3) - t490;
t280 = t227 * mrSges(5,1) + t223 * mrSges(5,3);
t282 = mrSges(4,1) * t227 - mrSges(4,2) * t223;
t257 = t222 * t223 + t226 * t227;
t354 = t223 * t226;
t258 = t222 * t227 - t354;
t475 = t257 * mrSges(6,1) - t258 * mrSges(6,2);
t481 = t280 + t282 + t475;
t229 = cos(qJ(1));
t225 = sin(qJ(1));
t350 = t225 * t228;
t143 = t223 * t350 + t227 * t229;
t347 = t229 * t223;
t144 = t225 * t349 - t347;
t260 = t143 * t222 + t144 * t226;
t443 = t143 * t226 - t144 * t222;
t480 = mrSges(6,1) * t443 - t260 * mrSges(6,2);
t331 = qJD(1) * qJD(2);
t176 = qJDD(1) * t224 + t228 * t331;
t338 = qJD(3) * t165;
t83 = qJDD(2) * t223 + t176 * t227 - t338;
t417 = t83 / 0.2e1;
t84 = qJD(3) * t166 - t227 * qJDD(2) + t176 * t223;
t416 = -t84 / 0.2e1;
t175 = t228 * qJDD(1) - t224 * t331;
t163 = qJDD(3) - t175;
t407 = t163 / 0.2e1;
t472 = t175 / 0.2e1;
t479 = t176 / 0.2e1;
t438 = qJD(3) - qJD(5);
t97 = t438 * t258;
t478 = qJD(4) * t223 + t212;
t210 = Ifges(3,4) * t341;
t159 = Ifges(4,4) * t165;
t373 = t165 * Ifges(5,5);
t463 = t166 * t489 + t196 * t488 - t159 + t373;
t158 = Ifges(5,5) * t166;
t69 = t196 * Ifges(5,6) + t165 * Ifges(5,3) + t158;
t477 = Ifges(3,1) * t342 + Ifges(3,5) * qJD(2) + t223 * t69 + t227 * t463 + t210;
t410 = pkin(3) + pkin(4);
t42 = -t196 * t410 + t483;
t189 = t196 * qJ(4);
t93 = t223 * t153 + t227 * t185;
t60 = pkin(8) * t165 + t93;
t49 = t189 + t60;
t13 = -t222 * t49 + t226 * t42;
t14 = t222 * t42 + t226 * t49;
t380 = Ifges(3,4) * t224;
t462 = Ifges(3,2) * t228;
t274 = t380 + t462;
t476 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t274 / 0.2e1 + t487 * t491 + t488 * t492 + t486 * t405 + t493;
t152 = qJDD(5) - t163;
t17 = qJD(5) * t88 + t222 * t84 + t226 * t83;
t18 = -qJD(5) * t259 - t222 * t83 + t226 * t84;
t474 = -t17 * Ifges(6,4) / 0.2e1 - t18 * Ifges(6,2) / 0.2e1 - t152 * Ifges(6,6) / 0.2e1;
t424 = t17 / 0.2e1;
t423 = t18 / 0.2e1;
t473 = -m(5) - m(6);
t408 = t152 / 0.2e1;
t339 = qJD(2) * t228;
t471 = pkin(6) * t339;
t470 = -mrSges(3,3) + mrSges(2,2);
t469 = -mrSges(5,3) + mrSges(4,2);
t467 = (-Ifges(4,4) + Ifges(5,5)) * t84 + t489 * t83 + t488 * t163;
t186 = t409 * t223;
t105 = t186 * t222 + t187 * t226;
t466 = -qJD(5) * t105 + t222 * t485 + t226 * t484;
t104 = t186 * t226 - t187 * t222;
t465 = qJD(5) * t104 + t222 * t484 - t226 * t485;
t178 = t226 * qJ(4) - t222 * t410;
t459 = -qJD(5) * t178 - t222 * t483 - t226 * t60;
t177 = -t222 * qJ(4) - t226 * t410;
t458 = qJD(5) * t177 - t222 * t60 + t226 * t483;
t284 = mrSges(3,1) * t228 - mrSges(3,2) * t224;
t457 = -t284 - mrSges(2,1);
t322 = t410 * t223;
t361 = qJ(4) * t227;
t253 = -t322 + t361;
t456 = t196 * t253 + t478;
t392 = pkin(3) * t223;
t263 = -t361 + t392;
t455 = t196 * t263 - t478;
t131 = t222 * t352 - t224 * t354;
t132 = t257 * t224;
t299 = t131 * mrSges(6,1) + t132 * mrSges(6,2);
t453 = qJD(2) * mrSges(3,1) - mrSges(4,1) * t165 - mrSges(4,2) * t166 - mrSges(3,3) * t342;
t452 = t224 * t487 + t228 * t449;
t451 = t224 * t488 + t228 * t446;
t211 = pkin(6) * t342;
t184 = -qJD(2) * pkin(2) + t211;
t365 = t227 * mrSges(5,3);
t279 = t223 * mrSges(5,1) - t365;
t281 = mrSges(4,1) * t223 + mrSges(4,2) * t227;
t255 = qJ(4) * t166 - t184;
t68 = pkin(3) * t165 - t255;
t450 = t184 * t281 + t68 * t279;
t448 = t223 * t488 + t227 * t486;
t375 = Ifges(5,5) * t227;
t377 = Ifges(4,4) * t227;
t447 = t223 * t489 - t375 + t377;
t445 = t163 * t487 - t486 * t84 + t488 * t83;
t161 = t175 * pkin(6);
t162 = t176 * pkin(6);
t444 = t161 * t228 + t162 * t224;
t138 = qJDD(2) * pkin(7) + t161;
t335 = qJD(3) * t227;
t360 = qJDD(1) * pkin(1);
t95 = -pkin(2) * t175 - pkin(7) * t176 - t360;
t30 = t227 * t138 + t153 * t335 - t185 * t337 + t223 * t95;
t31 = -t223 * t138 - t153 * t337 - t185 * t335 + t227 * t95;
t442 = -t223 * t31 + t227 * t30;
t441 = t490 * t224;
t19 = t163 * qJ(4) + t196 * qJD(4) + t30;
t247 = qJDD(4) - t31;
t21 = -pkin(3) * t163 + t247;
t440 = t19 * t227 + t21 * t223;
t437 = -t83 * Ifges(5,5) / 0.2e1 - t163 * Ifges(5,6) / 0.2e1 + Ifges(4,4) * t417 + Ifges(4,6) * t407 + (Ifges(5,3) + Ifges(4,2)) * t416;
t434 = pkin(7) * (-m(4) + t473);
t362 = qJ(4) * t223;
t433 = -t227 * t410 - t362;
t432 = m(6) * pkin(4) + mrSges(4,1) + mrSges(5,1);
t264 = pkin(3) * t227 + t362;
t179 = -pkin(2) - t264;
t283 = mrSges(3,1) * t224 + mrSges(3,2) * t228;
t391 = pkin(4) * t227;
t431 = t283 + t482 * t228 + (-m(6) * (t179 - t391) - m(5) * t179 + m(4) * pkin(2) + t481) * t224;
t8 = -pkin(8) * t83 - t163 * t410 + t247;
t9 = pkin(8) * t84 + t19;
t1 = qJD(5) * t13 + t222 * t8 + t226 * t9;
t2 = -qJD(5) * t14 - t222 * t9 + t226 * t8;
t430 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t427 = -t31 * mrSges(4,1) + t21 * mrSges(5,1) + t30 * mrSges(4,2) - t19 * mrSges(5,3);
t425 = Ifges(6,1) * t424 + Ifges(6,4) * t423 + Ifges(6,5) * t408;
t393 = Ifges(6,4) * t259;
t33 = Ifges(6,2) * t88 + Ifges(6,6) * t191 + t393;
t422 = -t33 / 0.2e1;
t421 = t33 / 0.2e1;
t85 = Ifges(6,4) * t88;
t34 = Ifges(6,1) * t259 + Ifges(6,5) * t191 + t85;
t420 = -t34 / 0.2e1;
t419 = t34 / 0.2e1;
t372 = t166 * Ifges(4,4);
t72 = -t165 * Ifges(4,2) + t196 * Ifges(4,6) + t372;
t418 = -t72 / 0.2e1;
t415 = t84 / 0.2e1;
t414 = -t88 / 0.2e1;
t412 = -t259 / 0.2e1;
t406 = -t165 / 0.2e1;
t403 = t166 / 0.2e1;
t402 = -t191 / 0.2e1;
t398 = t223 / 0.2e1;
t395 = mrSges(6,3) * t13;
t394 = mrSges(6,3) * t14;
t390 = pkin(6) * t224;
t387 = pkin(8) * t224;
t386 = -qJD(1) / 0.2e1;
t385 = qJD(3) / 0.2e1;
t384 = mrSges(5,2) * t165;
t383 = mrSges(5,2) * t166;
t382 = mrSges(4,3) * t165;
t381 = mrSges(4,3) * t166;
t379 = Ifges(3,4) * t228;
t363 = qJ(4) * t165;
t355 = t223 * t224;
t351 = t224 * t229;
t348 = t228 * t229;
t174 = t288 * qJD(2);
t345 = t219 + t216;
t180 = -pkin(1) - t345;
t346 = t223 * t174 + t180 * t335;
t198 = pkin(6) * t349;
t117 = t223 * t180 + t198;
t344 = t229 * pkin(1) + t225 * pkin(6);
t340 = qJD(2) * t224;
t336 = qJD(3) * t224;
t333 = qJD(4) * t227;
t329 = Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t152;
t326 = pkin(6) * t340;
t324 = pkin(7) * t337;
t323 = pkin(7) * t335;
t317 = t223 * t339;
t316 = t223 * t336;
t308 = -t341 / 0.2e1;
t301 = t335 / 0.2e1;
t51 = -t163 * mrSges(5,1) + t83 * mrSges(5,2);
t300 = t331 / 0.2e1;
t197 = pkin(6) * t353;
t116 = t180 * t227 - t197;
t295 = pkin(3) * t349 + qJ(4) * t353 + t345;
t294 = pkin(2) * t348 + pkin(7) * t351 + t344;
t290 = t320 * t224;
t145 = -t225 * t227 + t228 * t347;
t146 = t223 * t225 + t227 * t348;
t75 = t145 * t226 - t146 * t222;
t76 = t145 * t222 + t146 * t226;
t289 = mrSges(6,1) * t75 - mrSges(6,2) * t76;
t287 = qJDD(2) * pkin(2) - t162;
t106 = -qJ(4) * t228 + t117;
t285 = qJD(3) * t198 - t174 * t227 + t180 * t337;
t273 = -Ifges(4,2) * t223 + t377;
t272 = Ifges(4,2) * t227 + t378;
t269 = Ifges(3,5) * t228 - Ifges(3,6) * t224;
t266 = Ifges(5,3) * t223 + t375;
t265 = -Ifges(5,3) * t227 + t376;
t61 = -mrSges(6,2) * t191 + mrSges(6,3) * t88;
t62 = mrSges(6,1) * t191 - mrSges(6,3) * t259;
t261 = -t222 * t62 + t226 * t61;
t218 = t228 * pkin(3);
t82 = pkin(4) * t228 + t197 + t218 + (-t180 - t387) * t227;
t91 = pkin(8) * t355 + t106;
t35 = -t222 * t91 + t226 * t82;
t36 = t222 * t82 + t226 * t91;
t109 = -pkin(6) * t318 + t147;
t254 = pkin(1) * t283;
t250 = t224 * (Ifges(3,1) * t228 - t380);
t249 = t258 * t228;
t248 = t257 * t228;
t246 = t329 + t430;
t245 = t228 * t332 - t316;
t244 = t224 * t335 + t317;
t242 = -g(1) * t145 - g(2) * t143 - g(3) * t355;
t237 = Ifges(4,6) * t224 + t228 * t273;
t236 = Ifges(5,6) * t224 + t228 * t266;
t234 = qJ(4) * t83 + qJD(4) * t166 + t287;
t57 = (-t224 * t332 - t228 * t337) * pkin(6) + t346;
t96 = t438 * t257;
t220 = t229 * pkin(6);
t205 = qJ(4) * t340;
t193 = qJ(4) * t352;
t182 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t341;
t157 = pkin(2) - t433;
t154 = t281 * t224;
t125 = -t193 + (pkin(6) + t392) * t224;
t119 = qJD(1) * t248;
t118 = qJD(1) * t249;
t113 = mrSges(5,3) * t196 - t384;
t112 = -mrSges(5,1) * t196 + t383;
t111 = mrSges(4,1) * t196 - t381;
t110 = -mrSges(4,2) * t196 - t382;
t108 = pkin(6) * t319 + t356;
t107 = -t116 + t218;
t103 = t193 + (-pkin(6) - t322) * t224;
t101 = mrSges(5,1) * t165 - mrSges(5,3) * t166;
t100 = pkin(3) * t166 + t363;
t99 = qJD(1) * t290 - t356;
t98 = t109 + t204;
t66 = t189 + t93;
t65 = -pkin(3) * t196 + t494;
t63 = -t166 * t410 - t363;
t58 = t223 * t326 - t285;
t56 = (qJD(3) * t264 - t333) * t224 + (pkin(6) + t263) * t339;
t55 = -t165 * t410 + t255;
t54 = qJD(2) * t290 + t285;
t53 = -mrSges(5,2) * t84 + mrSges(5,3) * t163;
t52 = -mrSges(4,2) * t163 - mrSges(4,3) * t84;
t50 = mrSges(4,1) * t163 - mrSges(4,3) * t83;
t48 = -qJD(4) * t228 + t205 + t57;
t45 = qJD(2) * t248 + t224 * t97;
t44 = -qJD(2) * t249 + t224 * t96;
t43 = (qJD(3) * t433 + t333) * t224 + (-pkin(6) + t253) * t339;
t41 = -mrSges(6,1) * t88 + mrSges(6,2) * t259;
t40 = t205 + (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t352 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t223) * t228 + t346;
t39 = pkin(8) * t316 + qJD(2) * t241 + t285;
t38 = mrSges(4,1) * t84 + mrSges(4,2) * t83;
t37 = mrSges(5,1) * t84 - mrSges(5,3) * t83;
t20 = pkin(3) * t84 - t234;
t12 = -mrSges(6,2) * t152 + mrSges(6,3) * t18;
t11 = mrSges(6,1) * t152 - mrSges(6,3) * t17;
t10 = -t410 * t84 + t234;
t7 = -mrSges(6,1) * t18 + mrSges(6,2) * t17;
t6 = -qJD(5) * t36 - t222 * t40 + t226 * t39;
t5 = qJD(5) * t35 + t222 * t39 + t226 * t40;
t3 = [(t92 * mrSges(4,1) - t65 * mrSges(5,1) - t93 * mrSges(4,2) + t66 * mrSges(5,3) - t476 - t493) * t340 + t284 * t360 - (t223 * t463 + t227 * t72) * t336 / 0.2e1 + (-m(4) * pkin(6) * t287 + Ifges(3,1) * t176 + Ifges(3,4) * t472 + Ifges(3,5) * qJDD(2) + t20 * t279 + t266 * t415 + t273 * t416 - t300 * t462 + t69 * t301 + t449 * t407 + t446 * t417) * t224 - t287 * t154 + (Ifges(3,4) * t479 + Ifges(3,2) * t472 + t329 / 0.2e1 + Ifges(3,6) * qJDD(2) - t445 / 0.2e1 + Ifges(6,6) * t423 + Ifges(6,5) * t424 + Ifges(6,3) * t408 - Ifges(5,6) * t415 - Ifges(4,6) * t416 + pkin(6) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t175) + t379 * t300 - t488 * t417 - t487 * t407 + t427 + t430) * t228 + t379 * t479 + (-m(3) * t344 - m(4) * t294 - t76 * mrSges(6,1) - t75 * mrSges(6,2) + t473 * (t146 * pkin(3) + qJ(4) * t145 + t294) + t457 * t229 + t470 * t225 - t432 * t146 + t469 * t145 + t482 * t351) * g(2) + t44 * t421 + t132 * t425 + (qJD(2) * t236 - t265 * t336) * t405 + (qJD(2) * t237 - t272 * t336) * t406 + t317 * t418 + t45 * t419 + qJD(2) ^ 2 * t269 / 0.2e1 + t274 * t472 + (t260 * mrSges(6,1) + t443 * mrSges(6,2) + t473 * (-t144 * pkin(3) - qJ(4) * t143 + t220) + t470 * t229 + (-m(3) - m(4)) * t220 + t432 * t144 - t469 * t143 + (m(3) * pkin(1) - m(6) * t321 - (-m(6) * t409 + mrSges(6,3)) * t224 + (-m(4) - m(5)) * t256 + t441 - t457) * t225) * g(1) + (-t244 * t93 - t245 * t92 - t31 * t352) * mrSges(4,3) + (Ifges(6,4) * t132 - Ifges(6,2) * t131) * t423 + (Ifges(6,4) * t45 + Ifges(6,2) * t44) * t413 + t10 * t299 + (-t19 * mrSges(5,2) - t30 * mrSges(4,3) - t437) * t355 + (-qJDD(2) * mrSges(3,1) + t38) * t390 + (Ifges(6,1) * t132 - Ifges(6,4) * t131) * t424 + (Ifges(6,1) * t45 + Ifges(6,4) * t44) * t411 - t254 * t331 + (Ifges(6,5) * t132 - Ifges(6,6) * t131) * t408 + (Ifges(6,5) * t45 + Ifges(6,6) * t44) * t401 - t182 * t326 + (-t1 * t131 - t13 * t45 - t132 * t2 + t14 * t44) * mrSges(6,3) + t131 * t474 + (t21 * t352 - t244 * t66 + t245 * t65) * mrSges(5,2) + t68 * (mrSges(5,1) * t244 - mrSges(5,3) * t245) + t184 * (mrSges(4,1) * t244 + mrSges(4,2) * t245) + (t176 * t390 + t444) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t444) + (qJD(2) * t451 - t336 * t447) * t403 + (qJD(2) * t452 - t336 * t448) * t196 / 0.2e1 + t477 * t339 / 0.2e1 + t250 * t300 + t467 * t352 / 0.2e1 - t453 * t471 + m(4) * (t116 * t31 + t117 * t30 + t184 * t471 + t57 * t93 + t58 * t92) + t35 * t11 + t36 * t12 + t43 * t41 + t55 * (-mrSges(6,1) * t44 + mrSges(6,2) * t45) + t5 * t61 + t6 * t62 + t56 * t101 + t103 * t7 + t106 * t53 + t107 * t51 + t57 * t110 + t58 * t111 + t54 * t112 + t48 * t113 + t116 * t50 + t117 * t52 + t125 * t37 - pkin(1) * (-mrSges(3,1) * t175 + mrSges(3,2) * t176) + m(5) * (t106 * t19 + t107 * t21 + t125 * t20 + t48 * t66 + t54 * t65 + t56 * t68) + m(6) * (t1 * t36 + t10 * t103 + t13 * t6 + t14 * t5 + t2 * t35 + t43 * t55) + Ifges(2,3) * qJDD(1); (t385 * t449 + t386 * t452) * t196 - t258 * t425 + (-t1 * t257 + t118 * t14 + t119 * t13 + t2 * t258) * mrSges(6,3) + ((-t119 + t96) * mrSges(6,2) + (-t118 + t97) * mrSges(6,1)) * t55 + (Ifges(6,1) * t119 - Ifges(6,4) * t118) * t412 + (Ifges(6,4) * t119 - Ifges(6,2) * t118) * t414 + (Ifges(6,5) * t119 - Ifges(6,6) * t118) * t402 + t287 * t282 + (pkin(2) * t287 - t108 * t92 - t109 * t93 - t184 * t212) * m(4) + t257 * t474 + (-Ifges(6,4) * t258 - Ifges(6,2) * t257) * t423 + (-Ifges(6,1) * t258 - Ifges(6,4) * t257) * t424 + (-Ifges(6,5) * t258 - Ifges(6,6) * t257) * t408 + (t385 * t446 + t386 * t451) * t166 + t10 * t475 + (-m(4) * t345 - m(5) * t295 - m(6) * (t295 - t387) + t224 * mrSges(6,3) - t284 + (-m(6) * t391 - t481) * t228 - t441) * g(3) - (Ifges(6,4) * t411 + Ifges(6,2) * t413 + Ifges(6,6) * t401 + t394 + t421) * t97 + (-t324 - t109) * t110 + t182 * t211 + ((t237 / 0.2e1 - t236 / 0.2e1) * t165 - t92 * (mrSges(4,1) * t224 - mrSges(4,3) * t349) - t65 * (-mrSges(5,1) * t224 + mrSges(5,2) * t349) - t93 * (-mrSges(4,2) * t224 - mrSges(4,3) * t353) - t66 * (-mrSges(5,2) * t353 + mrSges(5,3) * t224) + (-t250 / 0.2e1 + t254) * qJD(1)) * qJD(1) + (-t324 - t98) * t113 - t20 * t280 - t118 * t422 + t265 * t415 + t272 * t416 + t119 * t420 + (Ifges(6,1) * t411 + Ifges(6,4) * t413 + Ifges(6,5) * t401 - t395 + t419) * t96 + (t179 * t20 + t455 * t68 - t65 * t99 - t66 * t98) * m(5) + ((t53 + t52) * t227 + (-t50 + t51) * t223 + ((-t223 * t93 - t227 * t92) * qJD(3) + t442) * m(4) + ((-t223 * t66 + t227 * t65) * qJD(3) + t440) * m(5)) * pkin(7) + (-t323 - t108) * t111 + t437 * t227 + (t323 - t99) * t112 + (t418 + t69 / 0.2e1) * t337 + (t225 * t431 + t350 * t434) * g(2) + (t229 * t431 + t348 * t434) * g(1) - t269 * t331 / 0.2e1 + (-t273 / 0.2e1 + t266 / 0.2e1) * t338 + (-Ifges(6,5) * t412 - Ifges(3,2) * t308 - Ifges(6,6) * t414 - Ifges(6,3) * t402 + t476) * t342 + (t335 * t65 - t337 * t66 + t440) * mrSges(5,2) + (-t335 * t92 - t337 * t93 + t442) * mrSges(4,3) + t447 * t417 + t448 * t407 + (t398 * t72 - t450) * t341 + t450 * qJD(3) + t453 * t212 + t455 * t101 + (t210 + t477) * t308 + t456 * t41 + t463 * t301 + t465 * t61 + t466 * t62 + (t1 * t105 + t10 * t157 + t104 * t2 + t13 * t466 + t14 * t465 + t456 * t55) * m(6) + t467 * t398 - pkin(2) * t38 + Ifges(3,3) * qJDD(2) + t104 * t11 + t105 * t12 + t157 * t7 - t161 * mrSges(3,2) - t162 * mrSges(3,1) + Ifges(3,6) * t175 + Ifges(3,5) * t176 + t179 * t37; (-t165 * t488 - t166 * t486) * t491 + (t473 * (-t143 * pkin(3) + qJ(4) * t144) + t469 * t144 + t432 * t143 + t480) * g(2) + (Ifges(5,3) * t166 - t373) * t406 + t72 * t403 - t427 + (-m(5) * t65 + t111 - t112 + t381) * t93 + (-m(5) * t66 - t110 - t113 - t382) * t92 + (-t165 * t489 + t158 - t372 + t69) * t492 + (-m(5) * t193 - (t365 + (-m(5) * pkin(3) - mrSges(5,1)) * t223) * t224 - m(6) * (-t224 * t322 + t193) - t299 + t154) * g(3) - t246 + (-pkin(3) * t21 + qJ(4) * t19 + qJD(4) * t66 - t100 * t68) * m(5) + t66 * t383 + t65 * t384 + t445 + (t55 * mrSges(6,1) + Ifges(6,4) * t412 + Ifges(6,2) * t414 + Ifges(6,6) * t402 - t394 + t422) * t259 - (-t55 * mrSges(6,2) + Ifges(6,1) * t412 + Ifges(6,4) * t414 + Ifges(6,5) * t402 + t395 + t420) * t88 + t458 * t61 + (t1 * t178 + t13 * t459 + t14 * t458 + t177 * t2 - t55 * t63) * m(6) + t459 * t62 + (-Ifges(4,2) * t166 - t159 + t463) * t405 + (t289 + t473 * (-t145 * pkin(3) + qJ(4) * t146) + t469 * t146 + t432 * t145) * g(1) - pkin(3) * t51 + qJ(4) * t53 - t63 * t41 - t100 * t101 + qJD(4) * t113 - t68 * (mrSges(5,1) * t166 + mrSges(5,3) * t165) + t177 * t11 + t178 * t12 - t184 * (mrSges(4,1) * t166 - mrSges(4,2) * t165); t226 * t11 + t222 * t12 + (t101 - t41) * t166 + t261 * qJD(5) + (-t113 - t261) * t196 + t51 + (t1 * t222 - t166 * t55 + t2 * t226 + t242 + t191 * (-t13 * t222 + t14 * t226)) * m(6) + (t166 * t68 - t196 * t66 + t21 + t242) * m(5); -t55 * (mrSges(6,1) * t259 + mrSges(6,2) * t88) + (Ifges(6,1) * t88 - t393) * t412 + t33 * t411 + (Ifges(6,5) * t88 - Ifges(6,6) * t259) * t402 - t13 * t61 + t14 * t62 - g(1) * t289 - g(2) * t480 + g(3) * t299 + (t13 * t88 + t14 * t259) * mrSges(6,3) + t246 + (-Ifges(6,2) * t259 + t34 + t85) * t414;];
tau = t3;
