% Calculate vector of inverse dynamics joint torques for
% S5RRRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:40
% EndTime: 2019-12-31 22:00:26
% DurationCPUTime: 25.25s
% Computational Cost: add. (7032->646), mult. (15787->869), div. (0->0), fcn. (10594->10), ass. (0->292)
t447 = -mrSges(6,1) - mrSges(5,1);
t446 = -mrSges(6,2) - mrSges(5,2);
t444 = Ifges(5,4) + Ifges(6,4);
t248 = sin(qJ(2));
t252 = cos(qJ(2));
t318 = qJD(1) * qJD(2);
t205 = qJDD(1) * t252 - t248 * t318;
t188 = t205 * pkin(6);
t470 = t188 * t252;
t289 = pkin(2) * t248 - pkin(7) * t252;
t201 = t289 * qJD(1);
t251 = cos(qJ(3));
t247 = sin(qJ(3));
t328 = qJD(1) * t248;
t306 = t247 * t328;
t134 = pkin(6) * t306 + t251 * t201;
t336 = t251 * t252;
t272 = pkin(3) * t248 - pkin(8) * t336;
t254 = -pkin(8) - pkin(7);
t307 = qJD(3) * t254;
t469 = -qJD(1) * t272 + t251 * t307 - t134;
t174 = t247 * t201;
t339 = t248 * t251;
t341 = t247 * t252;
t468 = t174 + (-pkin(6) * t339 - pkin(8) * t341) * qJD(1) - t247 * t307;
t327 = qJD(1) * t252;
t225 = qJD(3) - t327;
t218 = qJD(4) + t225;
t383 = t218 / 0.2e1;
t325 = qJD(2) * t251;
t193 = -t306 + t325;
t304 = t251 * t328;
t194 = qJD(2) * t247 + t304;
t246 = sin(qJ(4));
t250 = cos(qJ(4));
t120 = t193 * t246 + t194 * t250;
t397 = t120 / 0.2e1;
t291 = t250 * t193 - t194 * t246;
t400 = t291 / 0.2e1;
t235 = pkin(6) * t328;
t214 = -qJD(2) * pkin(2) + t235;
t138 = -pkin(3) * t193 + t214;
t290 = pkin(2) * t252 + pkin(7) * t248;
t210 = -pkin(1) - t290;
t182 = t210 * qJD(1);
t236 = pkin(6) * t327;
t215 = qJD(2) * pkin(7) + t236;
t123 = t251 * t182 - t215 * t247;
t89 = -pkin(8) * t194 + t123;
t79 = pkin(3) * t225 + t89;
t124 = t182 * t247 + t215 * t251;
t90 = pkin(8) * t193 + t124;
t84 = t246 * t90;
t39 = t250 * t79 - t84;
t454 = qJ(5) * t120;
t21 = t39 - t454;
t18 = pkin(4) * t218 + t21;
t74 = -pkin(4) * t291 + qJD(5) + t138;
t413 = -mrSges(5,2) * t138 - mrSges(6,2) * t74 + mrSges(5,3) * t39 + t18 * mrSges(6,3);
t443 = Ifges(5,5) + Ifges(6,5);
t445 = Ifges(5,1) + Ifges(6,1);
t453 = t444 * t291;
t435 = t120 * t445 + t218 * t443 + t453;
t457 = t435 / 0.2e1;
t467 = t443 * t383 + t445 * t397 + t444 * t400 - t413 + t457;
t466 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t245 = qJ(3) + qJ(4);
t239 = cos(t245);
t374 = pkin(3) * t251;
t209 = pkin(4) * t239 + t374;
t200 = pkin(2) + t209;
t231 = pkin(2) + t374;
t238 = sin(t245);
t286 = -mrSges(4,1) * t251 + mrSges(4,2) * t247;
t465 = m(4) * pkin(2) + m(5) * t231 + m(6) * t200 + t446 * t238 - t447 * t239 - t286;
t206 = qJDD(1) * t248 + t252 * t318;
t347 = qJDD(1) * pkin(1);
t127 = -pkin(2) * t205 - pkin(7) * t206 - t347;
t162 = qJDD(2) * pkin(7) + t188;
t321 = qJD(3) * t251;
t323 = qJD(3) * t247;
t50 = t247 * t127 + t251 * t162 + t182 * t321 - t215 * t323;
t464 = t50 * mrSges(4,2);
t51 = -t124 * qJD(3) + t251 * t127 - t162 * t247;
t463 = t51 * mrSges(4,1);
t442 = Ifges(5,2) + Ifges(6,2);
t441 = Ifges(5,6) + Ifges(6,6);
t440 = Ifges(5,3) + Ifges(6,3);
t196 = t246 * t251 + t247 * t250;
t420 = qJD(3) + qJD(4);
t129 = t420 * t196;
t268 = t196 * t252;
t150 = qJD(1) * t268;
t462 = t129 - t150;
t244 = -qJ(5) + t254;
t461 = -m(4) * pkin(7) + m(5) * t254 + m(6) * t244 - t466;
t86 = t250 * t90;
t40 = t246 * t79 + t86;
t431 = qJ(5) * t291;
t22 = t40 + t431;
t452 = t444 * t120;
t436 = t218 * t441 + t291 * t442 + t452;
t456 = t436 / 0.2e1;
t450 = -mrSges(5,1) * t138 - mrSges(6,1) * t74 + mrSges(5,3) * t40 + mrSges(6,3) * t22 + t456;
t460 = t441 * t383 + t444 * t397 + t442 * t400 + t450;
t361 = Ifges(3,4) * t248;
t281 = t252 * Ifges(3,2) + t361;
t458 = t22 * mrSges(6,2) + t40 * mrSges(5,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t281 / 0.2e1 - t18 * mrSges(6,1) - t39 * mrSges(5,1);
t190 = qJDD(3) - t205;
t181 = qJDD(4) + t190;
t109 = qJD(3) * t193 + qJDD(2) * t247 + t206 * t251;
t110 = -qJD(3) * t194 + qJDD(2) * t251 - t206 * t247;
t37 = qJD(4) * t291 + t109 * t250 + t110 * t246;
t38 = -qJD(4) * t120 - t109 * t246 + t110 * t250;
t455 = -t442 * t38 / 0.2e1 - t444 * t37 / 0.2e1 - t441 * t181 / 0.2e1;
t216 = t254 * t247;
t217 = t254 * t251;
t133 = t246 * t216 - t250 * t217;
t434 = -qJD(4) * t133 + t246 * t468 + t250 * t469;
t319 = qJD(4) * t250;
t320 = qJD(4) * t246;
t433 = t216 * t319 + t217 * t320 + t246 * t469 - t250 * t468;
t324 = qJD(2) * t252;
t262 = t247 * t324 + t248 * t321;
t305 = t247 * t327;
t427 = -t236 + (-t305 + t323) * pkin(3);
t249 = sin(qJ(1));
t253 = cos(qJ(1));
t451 = g(1) * t253 + g(2) * t249;
t411 = m(5) * pkin(3);
t404 = t109 / 0.2e1;
t403 = t110 / 0.2e1;
t389 = t190 / 0.2e1;
t401 = -t291 / 0.2e1;
t439 = t181 * t443 + t37 * t445 + t38 * t444;
t273 = t246 * t247 - t250 * t251;
t128 = t420 * t273;
t267 = t273 * t252;
t151 = qJD(1) * t267;
t438 = -pkin(4) * t328 - qJD(5) * t196 + t434 + (t128 - t151) * qJ(5);
t437 = -qJ(5) * t462 - qJD(5) * t273 + t433;
t432 = t411 + mrSges(4,1);
t161 = t273 * t248;
t430 = pkin(4) * t462 + t427;
t192 = t251 * t210;
t122 = -pkin(8) * t339 + t192 + (-pkin(6) * t247 - pkin(3)) * t252;
t227 = pkin(6) * t336;
t141 = t247 * t210 + t227;
t343 = t247 * t248;
t131 = -pkin(8) * t343 + t141;
t69 = t246 * t122 + t250 * t131;
t186 = Ifges(4,4) * t193;
t101 = t194 * Ifges(4,1) + t225 * Ifges(4,5) + t186;
t234 = Ifges(3,4) * t327;
t429 = Ifges(3,1) * t328 + Ifges(3,5) * qJD(2) + t251 * t101 + t234;
t428 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t193 + mrSges(4,2) * t194 + mrSges(3,3) * t328;
t426 = t181 * t440 + t37 * t443 + t38 * t441;
t335 = t252 * t253;
t157 = -t238 * t335 + t239 * t249;
t158 = t238 * t249 + t239 * t335;
t425 = t157 * t447 - t158 * t446;
t338 = t249 * t252;
t155 = t238 * t338 + t239 * t253;
t156 = t238 * t253 - t239 * t338;
t424 = -t155 * t447 + t156 * t446;
t189 = t206 * pkin(6);
t423 = t189 * t248 + t470;
t422 = -t247 * t51 + t251 * t50;
t421 = mrSges(5,1) * t238 - t239 * t446;
t419 = t194 * Ifges(4,5) + t193 * Ifges(4,6) + t225 * Ifges(4,3) + t120 * t443 + t218 * t440 + t291 * t441;
t418 = -m(3) - m(6) - m(5) - m(4);
t376 = pkin(3) * t247;
t208 = pkin(4) * t238 + t376;
t416 = -m(6) * t208 + mrSges(2,2) - mrSges(3,3);
t288 = t252 * mrSges(3,1) - mrSges(3,2) * t248;
t415 = t248 * t466 + mrSges(2,1) + t288;
t410 = m(6) * pkin(4);
t409 = t37 / 0.2e1;
t408 = t38 / 0.2e1;
t407 = Ifges(4,1) * t404 + Ifges(4,4) * t403 + Ifges(4,5) * t389;
t398 = -t120 / 0.2e1;
t390 = t181 / 0.2e1;
t387 = t194 / 0.2e1;
t384 = -t218 / 0.2e1;
t378 = pkin(3) * t194;
t375 = pkin(3) * t250;
t370 = g(3) * t248;
t240 = t248 * pkin(6);
t46 = t250 * t89 - t84;
t365 = mrSges(5,3) * t291;
t364 = mrSges(5,3) * t120;
t363 = mrSges(6,3) * t291;
t362 = mrSges(6,3) * t120;
t360 = Ifges(3,4) * t252;
t359 = Ifges(4,4) * t247;
t358 = Ifges(4,4) * t251;
t355 = t123 * mrSges(4,3);
t354 = t124 * mrSges(4,3);
t353 = t194 * Ifges(4,4);
t342 = t247 * t249;
t340 = t247 * t253;
t204 = t289 * qJD(2);
t326 = qJD(2) * t248;
t310 = pkin(6) * t326;
t330 = t251 * t204 + t247 * t310;
t207 = pkin(3) * t343 + t240;
t322 = qJD(3) * t248;
t237 = pkin(6) * t324;
t308 = Ifges(4,5) * t109 + Ifges(4,6) * t110 + Ifges(4,3) * t190;
t139 = pkin(3) * t262 + t237;
t100 = t193 * Ifges(4,2) + t225 * Ifges(4,6) + t353;
t301 = -t247 * t100 / 0.2e1;
t13 = -t38 * mrSges(6,1) + t37 * mrSges(6,2);
t45 = -t246 * t89 - t86;
t68 = t250 * t122 - t131 * t246;
t132 = t250 * t216 + t217 * t246;
t163 = -qJDD(2) * pkin(2) + t189;
t287 = mrSges(3,1) * t248 + mrSges(3,2) * t252;
t285 = mrSges(4,1) * t247 + mrSges(4,2) * t251;
t283 = Ifges(4,1) * t251 - t359;
t282 = Ifges(4,1) * t247 + t358;
t280 = -Ifges(4,2) * t247 + t358;
t279 = Ifges(4,2) * t251 + t359;
t278 = Ifges(3,5) * t252 - Ifges(3,6) * t248;
t277 = Ifges(4,5) * t251 - Ifges(4,6) * t247;
t276 = Ifges(4,5) * t247 + Ifges(4,6) * t251;
t275 = t200 * t252 - t244 * t248;
t274 = t252 * t231 - t248 * t254;
t271 = pkin(1) * t287;
t171 = -t247 * t335 + t249 * t251;
t169 = t247 * t338 + t251 * t253;
t23 = pkin(3) * t190 - pkin(8) * t109 + t51;
t31 = pkin(8) * t110 + t50;
t5 = t246 * t23 + t250 * t31 + t79 * t319 - t320 * t90;
t270 = t214 * t285;
t269 = t248 * (Ifges(3,1) * t252 - t361);
t67 = t272 * qJD(2) + (-t227 + (pkin(8) * t248 - t210) * t247) * qJD(3) + t330;
t87 = t247 * t204 + t210 * t321 + (-t248 * t325 - t252 * t323) * pkin(6);
t71 = -pkin(8) * t262 + t87;
t15 = t122 * t319 - t131 * t320 + t246 * t67 + t250 * t71;
t80 = -pkin(3) * t110 + t163;
t263 = -t247 * t322 + t251 * t324;
t261 = Ifges(4,5) * t248 + t252 * t283;
t260 = Ifges(4,6) * t248 + t252 * t280;
t259 = Ifges(4,3) * t248 + t252 * t277;
t6 = -qJD(4) * t40 + t250 * t23 - t246 * t31;
t16 = -qJD(4) * t69 - t246 * t71 + t250 * t67;
t2 = pkin(4) * t181 - qJ(5) * t37 - qJD(5) * t120 + t6;
t3 = qJ(5) * t38 + qJD(5) * t291 + t5;
t258 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) + t426;
t230 = pkin(4) + t375;
t212 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t327;
t183 = t285 * t248;
t172 = t251 * t335 + t342;
t170 = -t249 * t336 + t340;
t160 = t196 * t248;
t152 = pkin(4) * t273 - t231;
t140 = -pkin(6) * t341 + t192;
t137 = mrSges(4,1) * t225 - mrSges(4,3) * t194;
t136 = -mrSges(4,2) * t225 + mrSges(4,3) * t193;
t135 = -pkin(6) * t304 + t174;
t125 = pkin(4) * t160 + t207;
t98 = -qJ(5) * t273 + t133;
t97 = -qJ(5) * t196 + t132;
t94 = mrSges(5,1) * t218 - t364;
t93 = mrSges(6,1) * t218 - t362;
t92 = -mrSges(5,2) * t218 + t365;
t91 = -mrSges(6,2) * t218 + t363;
t88 = -qJD(3) * t141 + t330;
t83 = pkin(4) * t120 + t378;
t82 = -mrSges(4,2) * t190 + mrSges(4,3) * t110;
t81 = mrSges(4,1) * t190 - mrSges(4,3) * t109;
t73 = -qJD(2) * t268 + t161 * t420;
t72 = -qJD(2) * t267 - t129 * t248;
t66 = -mrSges(5,1) * t291 + mrSges(5,2) * t120;
t65 = -mrSges(6,1) * t291 + mrSges(6,2) * t120;
t62 = -mrSges(4,1) * t110 + mrSges(4,2) * t109;
t53 = -pkin(4) * t73 + t139;
t52 = -qJ(5) * t160 + t69;
t48 = t109 * Ifges(4,4) + t110 * Ifges(4,2) + t190 * Ifges(4,6);
t47 = -pkin(4) * t252 + qJ(5) * t161 + t68;
t29 = -mrSges(5,2) * t181 + mrSges(5,3) * t38;
t28 = -mrSges(6,2) * t181 + mrSges(6,3) * t38;
t27 = mrSges(5,1) * t181 - mrSges(5,3) * t37;
t26 = mrSges(6,1) * t181 - mrSges(6,3) * t37;
t25 = t46 - t454;
t24 = t45 - t431;
t17 = -pkin(4) * t38 + qJDD(5) + t80;
t14 = -mrSges(5,1) * t38 + mrSges(5,2) * t37;
t8 = qJ(5) * t73 - qJD(5) * t160 + t15;
t7 = pkin(4) * t326 - qJ(5) * t72 + qJD(5) * t161 + t16;
t1 = [t467 * t72 + t288 * t347 + (t361 + t281) * t205 / 0.2e1 + (-t123 * t263 - t124 * t262 - t339 * t51 - t343 * t50) * mrSges(4,3) + (-t160 * t441 - t161 * t443 - t252 * t440) * t390 - t252 * t463 + (-t340 * t411 - t170 * mrSges(4,1) - t169 * mrSges(4,2) + t447 * t156 + t446 * t155 + (m(3) * pkin(1) - m(6) * (-pkin(1) - t275) - m(5) * (-pkin(1) - t274) - m(4) * t210 + t415) * t249 + (t418 * pkin(6) + t416) * t253) * g(1) + t160 * t455 + t252 * (Ifges(3,4) * t206 + Ifges(3,2) * t205) / 0.2e1 + (-mrSges(3,1) * t240 + Ifges(3,5) * t248 + (-mrSges(3,2) * pkin(6) + Ifges(3,6)) * t252) * qJDD(2) + (t206 * t240 + t423 + t470) * mrSges(3,3) + t206 * t360 / 0.2e1 + t252 * t464 + (-t160 * t442 - t161 * t444 - t252 * t441) * t408 + (-t160 * t444 - t161 * t445 - t252 * t443) * t409 + (-t342 * t411 - t172 * mrSges(4,1) - t171 * mrSges(4,2) + t418 * (t253 * pkin(1) + t249 * pkin(6)) + t416 * t249 + t447 * t158 + t446 * t157 + (-m(4) * t290 - m(5) * t274 - m(6) * t275 - t415) * t253) * g(2) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t423) - (t308 + t426) * t252 / 0.2e1 + t428 * t237 + t429 * t324 / 0.2e1 + qJD(2) ^ 2 * t278 / 0.2e1 + t193 * (qJD(2) * t260 - t279 * t322) / 0.2e1 + t225 * (qJD(2) * t259 - t276 * t322) / 0.2e1 + m(4) * (t123 * t88 + t124 * t87 + t140 * t51 + t141 * t50 + (t163 * t248 + t214 * t324) * pkin(6)) + m(6) * (t125 * t17 + t18 * t7 + t2 * t47 + t22 * t8 + t3 * t52 + t53 * t74) + m(5) * (t138 * t139 + t15 * t40 + t16 * t39 + t207 * t80 + t5 * t69 + t6 * t68) - t48 * t343 / 0.2e1 + t301 * t324 + t62 * t240 - t439 * t161 / 0.2e1 + t206 * t248 * Ifges(3,1) + t3 * (t252 * mrSges(6,2) - mrSges(6,3) * t160) + t5 * (t252 * mrSges(5,2) - mrSges(5,3) * t160) + t17 * (mrSges(6,1) * t160 - mrSges(6,2) * t161) + t80 * (mrSges(5,1) * t160 - mrSges(5,2) * t161) + t2 * (-t252 * mrSges(6,1) + mrSges(6,3) * t161) + t6 * (-t252 * mrSges(5,1) + mrSges(5,3) * t161) + (t269 + t252 * (-Ifges(3,2) * t248 + t360)) * t318 / 0.2e1 - (t251 * t100 + t247 * t101) * t322 / 0.2e1 - t212 * t310 - t271 * t318 + t460 * t73 + (qJD(2) * t261 - t282 * t322) * t387 + (-Ifges(4,3) * t252 + t248 * t277) * t389 + (-Ifges(4,6) * t252 + t248 * t280) * t403 + (-Ifges(4,5) * t252 + t248 * t283) * t404 + t339 * t407 + t47 * t26 + t52 * t28 + t53 * t65 + t68 * t27 + t69 * t29 + t8 * t91 + t15 * t92 + t7 * t93 + t16 * t94 + (t440 * t383 + t441 * t400 + t443 * t397 + t419 / 0.2e1 - t124 * mrSges(4,2) + t123 * mrSges(4,1) - t458) * t326 + t125 * t13 + t87 * t136 + t88 * t137 + t139 * t66 + t214 * (mrSges(4,1) * t262 + mrSges(4,2) * t263) + t140 * t81 + t141 * t82 + t163 * t183 - pkin(1) * (-mrSges(3,1) * t205 + mrSges(3,2) * t206) + t207 * t14 + Ifges(2,3) * qJDD(1); -t467 * t128 + (t196 * t443 - t273 * t441) * t390 + (t196 * t444 - t273 * t442) * t408 + (t196 * t445 - t273 * t444) * t409 + t17 * (mrSges(6,1) * t273 + mrSges(6,2) * t196) + t80 * (mrSges(5,1) * t273 + mrSges(5,2) * t196) + (t150 * t22 - t151 * t18 - t196 * t2 - t273 * t3) * mrSges(6,3) + (t150 * t40 - t151 * t39 - t196 * t6 - t273 * t5) * mrSges(5,3) + (-t150 * t441 - t151 * t443) * t384 + (t270 + t301) * qJD(3) + t433 * t92 + (t132 * t6 + t133 * t5 + t138 * t427 - t231 * t80 + t39 * t434 + t40 * t433) * m(5) + t434 * t94 + (-t123 * (mrSges(4,1) * t248 - mrSges(4,3) * t336) - t124 * (-mrSges(4,2) * t248 - mrSges(4,3) * t341)) * qJD(1) + t273 * t455 + t150 * t456 + t151 * t457 + t451 * (t248 * t465 + t287) + (-g(3) * t465 + t451 * t461) * t252 + (t101 / 0.2e1 - t355) * t321 + (-t150 * t442 - t151 * t444) * t401 + (-t150 * t444 - t151 * t445) * t398 - t74 * (mrSges(6,1) * t150 - mrSges(6,2) * t151) - t138 * (mrSges(5,1) * t150 - mrSges(5,2) * t151) + (m(4) * ((-t123 * t251 - t124 * t247) * qJD(3) + t422) - t247 * t81 - t137 * t321 - t136 * t323 + t251 * t82) * pkin(7) + t422 * mrSges(4,3) + t427 * t66 - t428 * t236 - (-Ifges(3,2) * t328 + t234 + t429) * t327 / 0.2e1 + t430 * t65 + (t248 * t461 - t288) * g(3) - t419 * t328 / 0.2e1 + t163 * t286 - t323 * t354 + t437 * t91 + t438 * t93 + (t152 * t17 + t18 * t438 + t2 * t97 + t22 * t437 + t3 * t98 + t430 * t74) * m(6) + t439 * t196 / 0.2e1 + (-pkin(2) * t163 - t123 * t134 - t124 * t135 - t214 * t236) * m(4) + t100 * t305 / 0.2e1 + (t271 - t269 / 0.2e1) * qJD(1) ^ 2 + (t193 * t280 + t194 * t283 + t225 * t277) * qJD(3) / 0.2e1 - (t193 * t260 + t194 * t261 + t225 * t259) * qJD(1) / 0.2e1 - t270 * t327 - t278 * t318 / 0.2e1 - t460 * t129 + t251 * t48 / 0.2e1 + t212 * t235 + t276 * t389 + t279 * t403 + t282 * t404 + t247 * t407 + Ifges(3,3) * qJDD(2) - pkin(2) * t62 + t97 * t26 + t98 * t28 + (t440 * t384 + t443 * t398 + t441 * t401 + t458) * t328 + t132 * t27 + t133 * t29 - t135 * t136 - t134 * t137 + t152 * t13 - t188 * mrSges(3,2) - t189 * mrSges(3,1) + Ifges(3,6) * t205 + Ifges(3,5) * t206 - t231 * t14; (-m(6) * (-t208 * t338 - t209 * t253) - mrSges(4,2) * t170 + t432 * t169 + t424) * g(2) + (-m(6) * (-t208 * t335 + t209 * t249) + mrSges(4,2) * t172 - t432 * t171 + t425) * g(1) + (t443 * t384 + t445 * t398 + t444 * t401 + t413) * t291 - t464 + (m(5) * t376 + mrSges(6,1) * t238 + t421) * t370 + t308 - t66 * t378 - m(5) * (t138 * t378 + t39 * t45 + t40 * t46) + (-t441 * t384 - t444 * t398 - t442 * t401 + t450) * t120 + t463 + t27 * t375 + t194 * t354 + t193 * t355 - (-Ifges(4,2) * t194 + t101 + t186) * t193 / 0.2e1 - t194 * (Ifges(4,1) * t193 - t353) / 0.2e1 + (-t18 * t24 + t2 * t230 + t208 * t370 - t22 * t25 - t74 * t83) * m(6) + ((t246 * t3 + (-t18 * t246 + t22 * t250) * qJD(4)) * m(6) + (t92 + t91) * t319 + (-t94 - t93) * t320 + (t29 + t28) * t246) * pkin(3) + t100 * t387 + (t246 * t5 + t250 * t6 + (-t246 * t39 + t250 * t40) * qJD(4)) * t411 + t435 * t401 + t258 - t83 * t65 - t25 * t91 - t46 * t92 - t24 * t93 - t45 * t94 - t123 * t136 + t124 * t137 + g(3) * t183 - t214 * (mrSges(4,1) * t194 + mrSges(4,2) * t193) - t225 * (Ifges(4,5) * t193 - Ifges(4,6) * t194) / 0.2e1 + t230 * t26; t18 * t363 + t2 * t410 + t258 - t21 * t91 - t74 * (mrSges(6,1) * t120 + mrSges(6,2) * t291) - t138 * (mrSges(5,1) * t120 + mrSges(5,2) * t291) + (t364 + t94) * t40 + (t291 * t445 - t452) * t398 + t436 * t397 + (t365 - t92) * t39 + (-t120 * t441 + t291 * t443) * t384 + (-(-mrSges(6,1) - t410) * t238 + t421) * t370 + (t362 - m(6) * (-t18 + t21) + t93) * t22 + (t155 * t410 + t424) * g(2) + (-t157 * t410 + t425) * g(1) + (-t120 * t442 + t435 + t453) * t401 + (t26 + (-m(6) * t74 - t65) * t120) * pkin(4); -t291 * t91 + t120 * t93 + (g(3) * t252 + t18 * t120 - t22 * t291 - t451 * t248 + t17) * m(6) + t13;];
tau = t1;
