% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR14V3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14V3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:27
% EndTime: 2019-04-12 15:04:22
% DurationCPUTime: 20.88s
% Computational Cost: add. (6216->746), mult. (14750->1016), div. (0->0), fcn. (10864->10), ass. (0->336)
t152 = sin(qJ(4));
t157 = cos(qJ(4));
t269 = qJD(2) * t157;
t153 = sin(qJ(2));
t271 = qJD(1) * t153;
t127 = -t152 * t271 + t269;
t125 = qJD(5) - t127;
t318 = t125 * Ifges(6,6);
t117 = t127 * qJ(3);
t151 = sin(qJ(5));
t156 = cos(qJ(5));
t92 = qJD(3) * t151 + t117 * t156;
t337 = t92 * mrSges(6,3);
t128 = qJD(2) * t152 + t157 * t271;
t111 = t128 * qJ(3);
t150 = sin(qJ(6));
t155 = cos(qJ(6));
t54 = t111 * t155 - t150 * t92;
t55 = t111 * t150 + t155 * t92;
t405 = t54 * mrSges(7,1) - t55 * mrSges(7,2);
t425 = t111 * mrSges(6,1);
t407 = t425 + t405;
t158 = cos(qJ(2));
t270 = qJD(1) * t158;
t141 = qJD(4) - t270;
t87 = -t128 * t151 + t141 * t156;
t78 = qJD(6) - t87;
t342 = t78 * Ifges(7,3);
t88 = t128 * t156 + t141 * t151;
t53 = t125 * t150 + t155 * t88;
t349 = t53 * Ifges(7,5);
t52 = t125 * t155 - t150 * t88;
t350 = t52 * Ifges(7,6);
t17 = t342 + t349 + t350;
t341 = t87 * Ifges(6,2);
t357 = Ifges(6,4) * t88;
t40 = t318 + t341 + t357;
t414 = -t40 / 0.2e1 + t17 / 0.2e1;
t461 = t407 + t414 - t318 / 0.2e1 - t337;
t259 = qJD(1) * qJD(2);
t130 = -t158 * qJDD(1) + t153 * t259;
t126 = qJDD(4) + t130;
t131 = qJDD(1) * t153 + t158 * t259;
t71 = qJD(4) * t127 + qJDD(2) * t152 + t131 * t157;
t176 = -qJD(4) * t128 - t152 * t131;
t72 = qJDD(2) * t157 + t176;
t460 = -t126 * Ifges(5,6) / 0.2e1 - t72 * Ifges(5,2) / 0.2e1 - t71 * Ifges(5,4) / 0.2e1;
t448 = t126 / 0.2e1;
t449 = t72 / 0.2e1;
t450 = t71 / 0.2e1;
t459 = 0.2e1 * Ifges(5,1) * t450 + 0.2e1 * Ifges(5,4) * t449 + 0.2e1 * Ifges(5,5) * t448;
t257 = qJD(2) * qJD(3);
t137 = qJDD(2) * qJ(3) + t257;
t258 = qJD(1) * qJD(3);
t236 = t153 * t258;
t47 = qJ(3) * t176 + t157 * t137 - t152 * t236;
t31 = qJD(5) * t92 - t156 * qJDD(3) + t151 * t47;
t458 = t31 * mrSges(6,3);
t403 = t150 * t55 + t155 * t54;
t457 = mrSges(7,3) * t403;
t456 = t405 + t414;
t35 = qJD(5) * t87 + t126 * t151 + t156 * t71;
t70 = qJDD(5) - t72;
t14 = qJD(6) * t52 + t150 * t70 + t155 * t35;
t15 = -qJD(6) * t53 - t150 * t35 + t155 * t70;
t36 = -qJD(5) * t88 + t126 * t156 - t151 * t71;
t34 = qJDD(6) - t36;
t1 = Ifges(7,5) * t14 + Ifges(7,6) * t15 + Ifges(7,3) * t34;
t202 = t156 * qJD(3) - t117 * t151;
t437 = qJD(5) * t202;
t30 = qJDD(3) * t151 + t156 * t47 + t437;
t382 = t34 / 0.2e1;
t387 = t15 / 0.2e1;
t388 = t14 / 0.2e1;
t298 = qJ(3) * t131;
t438 = qJD(4) * t117;
t48 = t137 * t152 - t157 * (-t236 - t298) + t438;
t11 = -qJD(6) * t55 - t150 * t30 + t155 * t48;
t429 = t11 * mrSges(7,1);
t10 = qJD(6) * t54 + t150 * t48 + t155 * t30;
t430 = t10 * mrSges(7,2);
t406 = t429 - t430;
t8 = Ifges(6,4) * t35 + Ifges(6,2) * t36 + Ifges(6,6) * t70;
t455 = t406 - mrSges(6,3) * t30 + Ifges(7,5) * t388 + Ifges(7,6) * t387 + Ifges(7,3) * t382 - t8 / 0.2e1 + t1 / 0.2e1;
t424 = t111 * mrSges(6,2);
t454 = t202 * mrSges(6,3) - t424;
t229 = qJD(6) * t156 - qJD(4);
t263 = qJD(5) * t152;
t282 = t152 * t155;
t265 = qJD(4) * t156;
t228 = -qJD(6) + t265;
t418 = t228 * t157;
t452 = t150 * (-t151 * t263 + t418) + t229 * t282;
t351 = t30 * mrSges(6,2);
t32 = Ifges(6,6) * t36;
t33 = Ifges(6,5) * t35;
t64 = Ifges(6,3) * t70;
t7 = t33 + t32 + t64;
t451 = -t351 + t7 / 0.2e1 - t31 * mrSges(6,1) - t47 * mrSges(5,3) + t460;
t447 = t141 * Ifges(5,6) / 0.2e1;
t334 = Ifges(4,4) + Ifges(3,5);
t278 = t155 * t156;
t116 = -t150 * t157 + t152 * t278;
t179 = t116 * t158;
t277 = t156 * t157;
t191 = t150 * t152 + t155 * t277;
t201 = t229 * t157;
t272 = qJ(3) * qJD(1);
t261 = qJD(5) * t157;
t411 = t151 * t261 + t152 * t228;
t445 = t191 * qJD(3) + (-t150 * t201 - t155 * t411) * qJ(3) + t179 * t272;
t281 = t152 * t156;
t189 = t150 * t281 + t155 * t157;
t181 = t158 * t189;
t190 = -t150 * t277 + t282;
t444 = t190 * qJD(3) + (t150 * t411 - t155 * t201) * qJ(3) - t181 * t272;
t319 = t125 * Ifges(6,5);
t262 = qJD(5) * t156;
t290 = t150 * t156;
t73 = -t127 * t290 + t128 * t155;
t442 = t150 * t262 + t73;
t145 = Ifges(3,4) * t270;
t317 = t125 * Ifges(6,3);
t338 = t88 * Ifges(6,5);
t340 = t87 * Ifges(6,6);
t39 = t317 + t338 + t340;
t325 = Ifges(5,4) * t128;
t435 = t447 + t127 * Ifges(5,2) / 0.2e1 + t325 / 0.2e1;
t398 = -t202 * mrSges(6,1) + t92 * mrSges(6,2) + t117 * mrSges(5,3) + t435;
t167 = t317 / 0.2e1 + t340 / 0.2e1 + t338 / 0.2e1 + t39 / 0.2e1 - t398 - t435;
t420 = qJD(3) * mrSges(5,1);
t163 = (-t167 - t420) * t152;
t124 = Ifges(5,4) * t127;
t311 = t128 * Ifges(5,1);
t309 = t141 * Ifges(5,5);
t431 = -t309 / 0.2e1;
t67 = t124 + t309 + t311;
t174 = t431 - t124 / 0.2e1 - t311 / 0.2e1 - t111 * mrSges(5,3) - t67 / 0.2e1;
t419 = qJD(3) * mrSges(5,2);
t171 = -t174 + t419;
t428 = qJD(2) / 0.2e1;
t441 = (Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1) * qJD(2) + t171 * t157 + qJD(3) * mrSges(4,2) + (t153 * Ifges(4,1) - Ifges(4,5) * t158) * qJD(1) / 0.2e1 + Ifges(3,1) * t271 / 0.2e1 + t145 / 0.2e1 + t334 * t428 - t163;
t305 = t151 * t31;
t193 = -t202 * t262 + t305;
t264 = qJD(5) * t151;
t440 = -t156 * t30 + t92 * t264 - t193;
t21 = -mrSges(6,2) * t70 + mrSges(6,3) * t36;
t267 = qJD(4) * t111;
t330 = -mrSges(5,1) * t141 - mrSges(6,1) * t87 + mrSges(6,2) * t88 + mrSges(5,3) * t128;
t358 = mrSges(6,3) * t88;
t331 = -mrSges(6,1) * t125 - mrSges(7,1) * t52 + mrSges(7,2) * t53 + t358;
t361 = mrSges(6,1) * t70 + mrSges(7,1) * t15 - mrSges(7,2) * t14 - mrSges(6,3) * t35;
t60 = -mrSges(6,2) * t125 + mrSges(6,3) * t87;
t436 = qJD(5) * (-t151 * t60 + t156 * t331) + m(5) * (t47 + t267) - m(6) * (-t267 + t440) + m(7) * t193 + t330 * qJD(4) - t361 * t151 - mrSges(5,2) * t126 + mrSges(5,3) * t72 + t156 * t21;
t434 = t152 * (t150 * t229 + t155 * t264) - t155 * t418;
t372 = t70 / 0.2e1;
t380 = t36 / 0.2e1;
t381 = t35 / 0.2e1;
t390 = Ifges(6,1) * t381 + Ifges(6,4) * t380 + Ifges(6,5) * t372;
t432 = -Ifges(3,6) / 0.2e1;
t335 = mrSges(6,2) - mrSges(7,3);
t364 = -t125 / 0.2e1;
t427 = Ifges(6,5) * t364;
t426 = Ifges(6,6) * t364;
t260 = qJD(6) * t151;
t423 = t155 * t260 + t442;
t74 = t127 * t278 + t128 * t150;
t422 = -t150 * t260 + t155 * t262 - t74;
t51 = Ifges(7,4) * t52;
t19 = Ifges(7,1) * t53 + Ifges(7,5) * t78 + t51;
t301 = t155 * t19;
t339 = t88 * Ifges(6,1);
t77 = Ifges(6,4) * t87;
t41 = t319 + t77 + t339;
t421 = t41 + t301;
t314 = t127 * mrSges(5,3);
t93 = -mrSges(5,2) * t141 + t314;
t231 = t156 * t60 + t93;
t159 = cos(qJ(1));
t154 = sin(qJ(1));
t355 = g(1) * t154;
t415 = g(2) * t159 - t355;
t129 = (mrSges(4,1) * t158 + mrSges(4,3) * t153) * qJD(1);
t184 = t151 * t331 + t231;
t216 = mrSges(7,1) * t150 + mrSges(7,2) * t155;
t195 = -mrSges(6,3) - t216;
t186 = -mrSges(5,2) - t195;
t233 = mrSges(5,1) * t157 + mrSges(4,1);
t412 = t152 * t186 + mrSges(3,1) + t233;
t289 = t151 * t152;
t253 = t202 * t289;
t294 = t111 * t157;
t168 = m(5) * (-t117 * t152 + t294) + m(6) * (-t281 * t92 + t253 + t294) + t129 + m(7) * t253;
t409 = t157 * t330 + t168;
t408 = t48 * mrSges(5,1) + t47 * mrSges(5,2) - Ifges(5,5) * t71 - Ifges(5,6) * t72 - Ifges(5,3) * t126;
t404 = t10 * t155 - t11 * t150;
t402 = mrSges(5,2) - t216;
t377 = -t41 / 0.2e1;
t401 = t377 + t454;
t399 = mrSges(7,1) * t278 - mrSges(7,2) * t290 + mrSges(7,3) * t151 + mrSges(5,1);
t371 = -t78 / 0.2e1;
t374 = -t53 / 0.2e1;
t376 = -t52 / 0.2e1;
t397 = Ifges(7,5) * t374 + Ifges(7,6) * t376 + Ifges(7,3) * t371;
t395 = t397 - t456;
t2 = t14 * Ifges(7,4) + t15 * Ifges(7,2) + t34 * Ifges(7,6);
t393 = t2 / 0.2e1;
t392 = Ifges(7,1) * t388 + Ifges(7,4) * t387 + Ifges(7,5) * t382;
t389 = Ifges(3,4) / 0.2e1;
t375 = t52 / 0.2e1;
t373 = t53 / 0.2e1;
t370 = t78 / 0.2e1;
t369 = -t87 / 0.2e1;
t368 = t87 / 0.2e1;
t367 = -t88 / 0.2e1;
t366 = t88 / 0.2e1;
t365 = m(4) + m(5);
t363 = -t127 / 0.2e1;
t362 = -t128 / 0.2e1;
t360 = m(4) * qJ(3);
t359 = m(4) * qJ(3) ^ 2;
t356 = Ifges(7,4) * t53;
t354 = g(2) * t154;
t352 = g(3) * t153;
t333 = Ifges(4,5) - Ifges(3,4);
t332 = Ifges(4,6) - Ifges(3,6);
t275 = t157 * t158;
t288 = t151 * t153;
t188 = t156 * t275 + t288;
t103 = t188 * qJD(1);
t241 = t152 * t270;
t76 = t103 * t155 + t150 * t241;
t329 = -t434 - t76;
t75 = -t103 * t150 + t155 * t241;
t328 = -t452 - t75;
t113 = -t151 * t158 + t153 * t277;
t327 = mrSges(6,2) * t113;
t324 = Ifges(6,4) * t151;
t323 = Ifges(6,4) * t156;
t322 = Ifges(7,4) * t150;
t321 = Ifges(7,4) * t155;
t18 = Ifges(7,2) * t52 + Ifges(7,6) * t78 + t356;
t306 = t150 * t18;
t303 = t151 * t202;
t297 = t111 * t128;
t295 = t111 * t152;
t293 = t127 * t151;
t292 = t127 * t156;
t291 = t150 * t151;
t287 = t151 * t155;
t286 = t151 * t157;
t284 = t152 * t153;
t283 = t152 * t154;
t280 = t153 * t156;
t279 = t153 * t159;
t276 = t156 * t158;
t274 = t157 * t159;
t273 = t159 * t152;
t268 = qJD(2) * t158;
t266 = qJD(4) * t153;
t255 = mrSges(2,2) - mrSges(4,2) - mrSges(3,3);
t254 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t252 = t202 * t286;
t250 = -Ifges(4,5) / 0.2e1 + t389;
t249 = Ifges(4,6) / 0.2e1 + t432;
t248 = m(6) + m(7) + t365;
t247 = qJD(4) * t303;
t244 = t151 * t275;
t234 = -t260 / 0.2e1;
t225 = t159 * t248;
t101 = t113 * t159;
t99 = t113 * t154;
t224 = -g(1) * t101 - g(2) * t99;
t79 = -t113 * t150 + t153 * t282;
t80 = t113 * t155 + t150 * t284;
t222 = mrSges(7,1) * t79 - mrSges(7,2) * t80;
t59 = (-qJD(5) + t269) * t276 + (-t152 * t265 + (qJD(2) - t261) * t151) * t153;
t220 = qJD(6) * t284 + t59;
t217 = mrSges(6,1) * t151 + mrSges(6,2) * t156;
t215 = Ifges(6,1) * t156 - t324;
t214 = Ifges(7,1) * t155 - t322;
t213 = Ifges(7,1) * t150 + t321;
t212 = -Ifges(6,2) * t151 + t323;
t211 = -Ifges(7,2) * t150 + t321;
t210 = Ifges(7,2) * t155 + t322;
t209 = Ifges(6,5) * t156 - Ifges(6,6) * t151;
t208 = Ifges(7,5) * t155 - Ifges(7,6) * t150;
t207 = Ifges(7,5) * t150 + Ifges(7,6) * t155;
t205 = -t150 * t54 + t155 * t55;
t204 = t156 * t92 - t303;
t197 = -Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1;
t196 = mrSges(7,1) * t155 - mrSges(7,2) * t150 + mrSges(6,1);
t37 = -mrSges(7,2) * t78 + mrSges(7,3) * t52;
t38 = mrSges(7,1) * t78 - mrSges(7,3) * t53;
t194 = -t150 * t38 + t155 * t37 + t60;
t115 = t154 * t275 - t273;
t192 = -t115 * t151 + t154 * t280;
t82 = t115 * t156 + t154 * t288;
t112 = t153 * t286 + t276;
t187 = t111 * t217;
t180 = t116 * t153;
t178 = t189 * t153;
t177 = -qJD(6) * t113 + t152 * t268 + t157 * t266;
t175 = -g(1) * t225 - t248 * t354;
t172 = -t77 / 0.2e1 - t339 / 0.2e1 - t319 / 0.2e1 + t401;
t144 = Ifges(4,5) * t271;
t170 = Ifges(4,6) * t428 - Ifges(4,3) * t270 / 0.2e1 + t144 / 0.2e1 + qJD(2) * t432 - (Ifges(3,4) * t153 + t158 * Ifges(3,2)) * qJD(1) / 0.2e1 + t141 * Ifges(5,3) + t128 * Ifges(5,5) + t127 * Ifges(5,6) - t111 * mrSges(5,1) - t117 * mrSges(5,2);
t169 = t153 * t412 + t254 * t158;
t166 = t342 / 0.2e1 - t341 / 0.2e1 - t357 / 0.2e1 + t350 / 0.2e1 + t349 / 0.2e1 + t461;
t164 = -mrSges(6,1) * t36 + mrSges(6,2) * t35 - mrSges(5,1) * t126 + mrSges(5,3) * t71 - t184 * qJD(4) + m(7) * t247 + m(6) * (-t265 * t92 + t247 + t48) + m(5) * (t48 - t438);
t161 = qJD(1) ^ 2;
t160 = qJD(2) ^ 2;
t149 = t153 ^ 2;
t134 = mrSges(4,2) * t270 + qJD(2) * mrSges(4,3);
t123 = (-mrSges(5,1) * t152 - mrSges(5,2) * t157) * t153;
t121 = t158 * t274 + t283;
t120 = -t154 * t157 + t158 * t273;
t114 = t158 * t283 + t274;
t105 = t191 * qJ(3);
t104 = t190 * qJ(3);
t98 = t112 * t154;
t97 = qJ(3) * t180;
t96 = qJ(3) * t178;
t86 = t121 * t156 + t151 * t279;
t85 = t121 * t151 - t156 * t279;
t63 = -t111 * t278 + t117 * t150;
t62 = t111 * t290 + t117 * t155;
t50 = t120 * t150 + t155 * t86;
t49 = t120 * t155 - t150 * t86;
t26 = t150 * t177 + t155 * t220;
t25 = -t150 * t220 + t155 * t177;
t23 = qJD(3) * t178 + (qJD(2) * t181 + t153 * t452) * qJ(3);
t22 = -qJD(3) * t180 + (-qJD(2) * t179 + t153 * t434) * qJ(3);
t6 = -mrSges(7,2) * t34 + mrSges(7,3) * t15;
t5 = mrSges(7,1) * t34 - mrSges(7,3) * t14;
t3 = [(t137 * mrSges(4,2) - t333 * t131 + (-Ifges(4,3) - Ifges(3,2)) * t130 - t332 * qJDD(2) - t415 * mrSges(3,1) + (-t415 + t298) * mrSges(4,1) + (t250 * t270 + (-t152 * t184 + t409) * qJ(3) + t441) * qJD(2) + t408) * t158 + (qJDD(3) * mrSges(4,2) + (Ifges(4,1) + Ifges(3,1)) * t131 + t333 * t130 + t334 * qJDD(2) + (t48 * mrSges(5,3) + t459) * t157 + (t33 / 0.2e1 + t32 / 0.2e1 + t64 / 0.2e1 + t451 + t460) * t152 + (t152 * t174 + t157 * t167) * qJD(4) + (t249 * qJD(2) + (-t250 * t153 + (-t197 + t359) * t158) * qJD(1) + t170) * qJD(2) + (t129 + (mrSges(5,1) * qJD(4) + t330) * t157 + (-mrSges(5,2) * qJD(4) - t184) * t152 + t168) * qJD(3) + (-mrSges(4,1) * t130 - t160 * mrSges(4,2) + (0.2e1 * mrSges(4,3) + t360) * t131 + (mrSges(4,3) * t268 + (m(4) * qJD(3) - qJD(2) * mrSges(4,1)) * t153) * qJD(1) - g(2) * t225 + t248 * t355 + t164 * t157 - t436 * t152) * qJ(3) + t415 * t254) * t153 + (-Ifges(6,4) * t366 + Ifges(7,5) * t373 - Ifges(6,2) * t368 + Ifges(7,6) * t375 + Ifges(7,3) * t370 + t461) * (-t266 * t289 - t158 * t264 - qJD(2) * t280 + (t151 * t268 + t153 * t262) * t157) + (Ifges(7,5) * t26 + Ifges(7,6) * t25) * t370 + (Ifges(7,5) * t80 + Ifges(7,6) * t79) * t382 - t202 * (-mrSges(7,1) * t25 + mrSges(7,2) * t26) + m(7) * (-t10 * t97 + t11 * t96 + t22 * t55 + t23 * t54) + (-mrSges(2,1) * t159 - mrSges(5,1) * t121 - mrSges(6,1) * t86 - mrSges(7,1) * t50 - mrSges(7,2) * t49 + t335 * t85 + (mrSges(5,2) - mrSges(6,3)) * t120 + t255 * t154) * g(2) + (mrSges(2,1) * t154 + mrSges(5,1) * t115 + t114 * t186 + t159 * t255 + t192 * t335 + t196 * t82) * g(1) + (t10 * t79 - t11 * t80 + t25 * t55 - t26 * t54) * mrSges(7,3) + (t319 / 0.2e1 + Ifges(6,1) * t366 + Ifges(6,4) * t368 + t41 / 0.2e1 - t454) * t59 + (mrSges(6,1) * t48 - Ifges(6,4) * t381 - Ifges(6,2) * t380 - Ifges(6,6) * t372 + t455) * t112 + t48 * t327 + (Ifges(7,4) * t26 + Ifges(7,2) * t25) * t375 + (Ifges(7,4) * t80 + Ifges(7,2) * t79) * t387 + (0.2e1 * t390 + t458) * t113 + Ifges(2,3) * qJDD(1) - t31 * t222 - qJDD(3) * t123 + t149 * t258 * t360 + t80 * t392 + t79 * t393 + (Ifges(7,1) * t26 + Ifges(7,4) * t25) * t373 + (Ifges(7,1) * t80 + Ifges(7,4) * t79) * t388 + t25 * t18 / 0.2e1 + t26 * t19 / 0.2e1 + t22 * t37 + t23 * t38 + t96 * t5 - t97 * t6; -t233 * qJDD(3) - (-mrSges(7,1) * t328 + mrSges(7,2) * t329) * t202 + (g(2) * t98 - t10 * t189 - t11 * t116 + t328 * t55 - t329 * t54) * mrSges(7,3) + (-t155 * t224 + t189 * t31) * mrSges(7,1) + (Ifges(7,1) * t116 - Ifges(7,4) * t189) * t388 + (Ifges(7,4) * t116 - Ifges(7,2) * t189) * t387 + (Ifges(7,5) * t116 - Ifges(7,6) * t189) * t382 - t189 * t393 + (t335 * (t244 - t280) + t254 * t153 - t196 * t188 - t412 * t158) * g(3) + t332 * t130 + t334 * t131 + ((m(5) * t117 + t184) * qJD(3) + t436 * qJ(3) + Ifges(5,6) * t448 + Ifges(5,2) * t449 + Ifges(5,4) * t450 + (t151 * t166 - t156 * t172 + t171) * qJD(4) - Ifges(6,3) * t372 - Ifges(6,6) * t380 - Ifges(6,5) * t381 - t451) * t157 + (-Ifges(7,5) * t434 - Ifges(7,6) * t452) * t370 + (-Ifges(7,1) * t434 - Ifges(7,4) * t452) * t373 + (-Ifges(7,4) * t434 - Ifges(7,2) * t452) * t375 + (-t452 / 0.2e1 - t75 / 0.2e1) * t18 + (-t434 / 0.2e1 - t76 / 0.2e1) * t19 - t163 * qJD(4) + (t151 * t172 + t156 * t166) * t263 + t455 * t289 + (mrSges(5,2) * qJDD(3) + t330 * qJD(3) + t209 * t372 + t212 * t380 + t215 * t381 + t459 + (t184 * t270 + t164) * qJ(3) + (t217 + mrSges(5,3)) * t48) * t152 + (t390 + t458) * t281 + (t31 * t116 + t150 * t224) * mrSges(7,2) + (mrSges(6,1) * t101 + (-t112 * t335 + t169) * t159) * g(1) + ((-t144 / 0.2e1 + (qJ(3) * mrSges(4,1) + t389) * t271 + (qJ(3) * mrSges(4,2) - t249) * qJD(2) - t170) * t153 + (-t145 / 0.2e1 + Ifges(4,5) * t270 / 0.2e1 + t197 * t271 + (-mrSges(4,3) * t271 - t409) * qJ(3) - t441) * t158) * qJD(1) + (Ifges(7,1) * t76 + Ifges(7,4) * t75) * t374 + t444 * t38 + t445 * t37 + (-qJD(3) * t252 + t10 * t105 + t104 * t11 + t444 * t54 + t445 * t55) * m(7) + (Ifges(7,4) * t76 + Ifges(7,2) * t75) * t376 + t137 * mrSges(4,3) - g(2) * (-mrSges(6,1) * t99 + mrSges(6,2) * t98) + (-Ifges(6,4) * t367 - Ifges(6,2) * t369 + t337 + t395 - t425 - t426) * (qJD(1) * t244 - t156 * t271) + (Ifges(6,1) * t367 + Ifges(6,4) * t369 + t401 + t427) * t103 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) + (Ifges(7,5) * t76 + Ifges(7,6) * t75) * t371 + t116 * t392 + (t134 + m(5) * t295 + m(6) * (t277 * t92 - t252 + t295)) * qJD(3) + (-mrSges(4,2) * t130 + qJDD(2) * mrSges(4,3) + (t137 + t257) * m(4) - t248 * t352 + t175 * t158) * qJ(3) + t169 * t354 - t153 * t161 * t158 * t359 + t104 * t5 + t105 * t6; -qJDD(2) * mrSges(4,1) - t72 * mrSges(5,1) + t131 * mrSges(4,2) + t71 * mrSges(5,2) - qJD(2) * t134 - t127 * t93 - t74 * t37 - t73 * t38 - t330 * t128 + t365 * qJDD(3) + t248 * t158 * g(3) + (qJD(5) * t194 - t127 * t60 + t361) * t156 + (-t150 * t5 + t155 * t6 + t21 + (-t150 * t37 - t155 * t38) * qJD(6) + t125 * t331) * t151 - m(5) * (t117 * t127 + t297) - (t149 * t161 + t160) * t360 + (-qJD(1) * t129 + t175) * t153 + ((qJD(5) * t205 - t31) * t156 + (-qJD(6) * t403 + t404 - t437) * t151 + t293 * t202 - t54 * t73 - t55 * t74) * m(7) + (t125 * t204 + t151 * t30 - t156 * t31 - t297) * m(6); (-t442 / 0.2e1 + t155 * t234) * t18 + (-t314 - m(6) * (t117 - t204) + t231 + (-m(7) * t202 + t331) * t151) * t111 - (mrSges(7,1) * t423 + mrSges(7,2) * t422) * t202 + t421 * t262 / 0.2e1 + (-t10 * t291 - t11 * t287 + t289 * t352 - t422 * t54 - t423 * t55) * mrSges(7,3) + (Ifges(7,5) * t74 + Ifges(7,6) * t73) * t371 + t395 * t293 + t156 * t430 + (t114 * t399 + t115 * t402) * g(2) + (t120 * t399 + t121 * t402) * g(1) + (Ifges(6,5) * t367 - Ifges(5,2) * t363 + Ifges(6,6) * t369 + Ifges(6,3) * t364 + t398 - t420 + t447) * t128 + t456 * t264 + (t150 * t234 - t74 / 0.2e1) * t19 + (-t325 + t39) * t362 - (-mrSges(7,1) * t116 + mrSges(7,2) * t189) * t352 - t2 * t291 / 0.2e1 + (Ifges(7,1) * t74 + Ifges(7,4) * t73) * t374 + qJD(5) * t187 + (-g(1) * t121 - g(2) * t115 - t157 * t352 + t202 * t292 + t293 * t92 - t440) * mrSges(6,3) - m(7) * (t54 * t62 + t55 * t63) + (Ifges(7,4) * t74 + Ifges(7,2) * t73) * t376 - t408 - t156 * t1 / 0.2e1 + t156 * t8 / 0.2e1 + (t124 + t67) * t363 - g(3) * t123 + (t125 * t209 + t212 * t87 + t215 * t88) * qJD(5) / 0.2e1 - t330 * t117 - t156 * t429 + (Ifges(5,1) * t362 + t209 * t364 + t212 * t369 + t215 * t367 - t187 - t419 + t431) * t127 + (-t207 * t260 + (Ifges(7,3) * t151 + t156 * t208) * qJD(5)) * t370 + (Ifges(6,5) * t151 + Ifges(6,6) * t156) * t372 + (-t213 * t260 + (Ifges(7,5) * t151 + t156 * t214) * qJD(5)) * t373 + (-t210 * t260 + (Ifges(7,6) * t151 + t156 * t211) * qJD(5)) * t375 + t292 * t377 + (Ifges(6,2) * t156 + t324) * t380 + (Ifges(6,1) * t151 + t323) * t381 + (-Ifges(7,3) * t156 + t151 * t208) * t382 + (-Ifges(7,6) * t156 + t151 * t211) * t387 + (-Ifges(7,5) * t156 + t151 * t214) * t388 + t151 * t390 + t287 * t392 + t216 * t305 - t62 * t38 - t63 * t37 + (-g(1) * t120 - g(2) * t114 - t152 * t352 + t48) * (-mrSges(6,1) * t156 + mrSges(6,2) * t151); -t351 + t7 + (t112 * t196 + t327) * g(3) + (t301 / 0.2e1 - t306 / 0.2e1 + t211 * t375 + t214 * t373 + t208 * t370 - t202 * t216 - t457) * qJD(6) + (-g(1) * t86 - g(2) * t82 - g(3) * t113 + t404) * mrSges(7,3) + (-t331 + t358) * t92 - (-m(7) * (-t205 + t92) + t194) * t202 + (mrSges(6,2) * t86 + t196 * t85) * g(1) - t196 * t31 + (mrSges(6,2) * t82 - t192 * t196) * g(2) + t40 * t366 + t306 * t368 + t207 * t382 + t210 * t387 + t213 * t388 + t150 * t392 + t155 * t393 + (-t195 * t202 + t208 * t371 + t211 * t376 + t214 * t374 - t424 + t427 + t457) * t87 + (t397 - t407 - t426) * t88 + (Ifges(6,1) * t87 + t17 - t357) * t367 + (-Ifges(6,2) * t88 + t421 + t77) * t369; t202 * (mrSges(7,1) * t53 + mrSges(7,2) * t52) + (Ifges(7,1) * t52 - t356) * t374 + t18 * t373 + (Ifges(7,5) * t52 - Ifges(7,6) * t53) * t371 - t54 * t37 + t55 * t38 - g(1) * (mrSges(7,1) * t49 - mrSges(7,2) * t50) - g(2) * ((t114 * t155 - t150 * t82) * mrSges(7,1) + (-t114 * t150 - t155 * t82) * mrSges(7,2)) - g(3) * t222 + (t52 * t54 + t53 * t55) * mrSges(7,3) + t1 + (-Ifges(7,2) * t53 + t19 + t51) * t376 + t406;];
tau  = t3;
