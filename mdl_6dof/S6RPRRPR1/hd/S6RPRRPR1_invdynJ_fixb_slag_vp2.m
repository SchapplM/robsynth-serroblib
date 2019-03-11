% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:32
% EndTime: 2019-03-09 04:57:54
% DurationCPUTime: 13.00s
% Computational Cost: add. (12411->645), mult. (27056->861), div. (0->0), fcn. (19180->18), ass. (0->302)
t269 = qJ(3) + qJ(4);
t259 = pkin(11) + t269;
t240 = sin(t259);
t241 = cos(t259);
t274 = sin(qJ(6));
t387 = mrSges(7,2) * t274;
t455 = t240 * t387 + t241 * (m(7) * pkin(9) + mrSges(7,3));
t267 = qJD(3) + qJD(4);
t278 = cos(qJ(6));
t275 = sin(qJ(4));
t276 = sin(qJ(3));
t279 = cos(qJ(4));
t280 = cos(qJ(3));
t216 = -t275 * t276 + t279 * t280;
t206 = t216 * qJD(1);
t217 = t275 * t280 + t276 * t279;
t207 = t217 * qJD(1);
t270 = sin(pkin(11));
t272 = cos(pkin(11));
t304 = t206 * t270 + t272 * t207;
t135 = t267 * t278 - t274 * t304;
t136 = t267 * t274 + t278 * t304;
t371 = mrSges(6,1) * t267 + mrSges(7,1) * t135 - mrSges(7,2) * t136 - mrSges(6,3) * t304;
t271 = sin(pkin(10));
t243 = pkin(1) * t271 + pkin(7);
t228 = t243 * qJD(1);
t324 = pkin(8) * qJD(1) + t228;
t349 = qJD(2) * t276;
t170 = t280 * t324 + t349;
t166 = t279 * t170;
t260 = t280 * qJD(2);
t169 = -t276 * t324 + t260;
t167 = qJD(3) * pkin(3) + t169;
t117 = t167 * t275 + t166;
t369 = qJ(5) * t206;
t106 = t117 + t369;
t368 = t106 * t270;
t164 = t275 * t170;
t116 = t279 * t167 - t164;
t199 = t207 * qJ(5);
t105 = t116 - t199;
t99 = pkin(4) * t267 + t105;
t54 = t272 * t99 - t368;
t48 = -pkin(5) * t267 - t54;
t454 = -m(6) * t54 + m(7) * t48 - t371;
t102 = t272 * t106;
t55 = t270 * t99 + t102;
t49 = pkin(9) * t267 + t55;
t273 = cos(pkin(10));
t245 = -pkin(1) * t273 - pkin(2);
t263 = t280 * pkin(3);
t225 = t245 - t263;
t210 = t225 * qJD(1);
t157 = -pkin(4) * t206 + qJD(5) + t210;
t323 = t272 * t206 - t207 * t270;
t72 = -pkin(5) * t323 - pkin(9) * t304 + t157;
t17 = -t274 * t49 + t278 * t72;
t453 = t17 * mrSges(7,1);
t18 = t274 * t72 + t278 * t49;
t452 = t18 * mrSges(7,2);
t261 = sin(t269);
t262 = cos(t269);
t451 = mrSges(5,1) * t261 + mrSges(6,1) * t240 + mrSges(5,2) * t262 + mrSges(6,2) * t241;
t388 = mrSges(7,1) * t278;
t450 = t387 - t388;
t342 = qJD(1) * qJD(3);
t221 = qJDD(1) * t280 - t276 * t342;
t222 = qJDD(1) * t276 + t280 * t342;
t295 = t217 * qJD(4);
t138 = -qJD(1) * t295 + t221 * t279 - t222 * t275;
t227 = t245 * qJDD(1);
t171 = -pkin(3) * t221 + t227;
t108 = -pkin(4) * t138 + qJDD(5) + t171;
t294 = t216 * qJD(4);
t137 = qJD(1) * t294 + t221 * t275 + t222 * t279;
t82 = -t137 * t270 + t138 * t272;
t83 = t137 * t272 + t138 * t270;
t19 = -pkin(5) * t82 - pkin(9) * t83 + t108;
t265 = qJDD(3) + qJDD(4);
t186 = t228 * t280 + t349;
t226 = t243 * qJDD(1);
t146 = -qJD(3) * t186 + t280 * qJDD(2) - t226 * t276;
t123 = qJDD(3) * pkin(3) - pkin(8) * t222 + t146;
t348 = qJD(3) * t276;
t145 = qJD(3) * t260 + t276 * qJDD(2) + t280 * t226 - t228 * t348;
t128 = pkin(8) * t221 + t145;
t46 = -qJD(4) * t117 + t279 * t123 - t128 * t275;
t29 = pkin(4) * t265 - qJ(5) * t137 - qJD(5) * t207 + t46;
t345 = qJD(4) * t279;
t346 = qJD(4) * t275;
t45 = t275 * t123 + t279 * t128 + t167 * t345 - t170 * t346;
t31 = qJ(5) * t138 + qJD(5) * t206 + t45;
t11 = t270 * t29 + t272 * t31;
t8 = pkin(9) * t265 + t11;
t2 = qJD(6) * t17 + t19 * t274 + t278 * t8;
t3 = -qJD(6) * t18 + t19 * t278 - t274 * t8;
t449 = t2 * t278 - t274 * t3;
t268 = qJ(1) + pkin(10);
t255 = sin(t268);
t256 = cos(t268);
t448 = g(1) * t256 + g(2) * t255;
t447 = -t262 * mrSges(5,1) - t241 * mrSges(6,1) + t261 * mrSges(5,2) + (mrSges(6,2) - mrSges(7,3)) * t240;
t52 = qJD(6) * t135 + t265 * t274 + t278 * t83;
t418 = t52 / 0.2e1;
t53 = -qJD(6) * t136 + t265 * t278 - t274 * t83;
t417 = t53 / 0.2e1;
t81 = qJDD(6) - t82;
t415 = t81 / 0.2e1;
t443 = -m(6) - m(7);
t442 = t221 / 0.2e1;
t408 = t265 / 0.2e1;
t441 = -t304 / 0.2e1;
t440 = t304 / 0.2e1;
t439 = -t323 / 0.2e1;
t16 = -mrSges(7,1) * t53 + mrSges(7,2) * t52;
t74 = mrSges(6,1) * t265 - mrSges(6,3) * t83;
t438 = t16 - t74;
t437 = Ifges(6,4) * t304;
t436 = Ifges(6,4) * t323;
t435 = Ifges(5,5) * t217;
t434 = Ifges(5,6) * t216;
t313 = mrSges(7,1) * t274 + mrSges(7,2) * t278;
t433 = t313 * t48;
t139 = -mrSges(6,2) * t267 + mrSges(6,3) * t323;
t142 = qJD(6) - t323;
t88 = -mrSges(7,2) * t142 + mrSges(7,3) * t135;
t89 = mrSges(7,1) * t142 - mrSges(7,3) * t136;
t305 = -t274 * t89 + t278 * t88;
t302 = -t139 - t305;
t389 = pkin(8) + t243;
t211 = t389 * t276;
t212 = t389 * t280;
t153 = -t275 * t211 + t279 * t212;
t320 = t241 * pkin(5) + t240 * pkin(9);
t162 = qJD(3) * t216 + t294;
t163 = -qJD(3) * t217 - t295;
t111 = t162 * t272 + t163 * t270;
t159 = t216 * t270 + t217 * t272;
t343 = qJD(6) * t278;
t299 = t111 * t274 + t159 * t343;
t430 = t145 * t280 - t146 * t276;
t20 = mrSges(7,1) * t81 - mrSges(7,3) * t52;
t21 = -mrSges(7,2) * t81 + mrSges(7,3) * t53;
t429 = -t274 * t20 + t278 * t21;
t427 = -m(5) - m(4) - m(3);
t426 = 0.2e1 * t408;
t425 = t241 * t450 + t447;
t423 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t232 = -mrSges(4,1) * t280 + mrSges(4,2) * t276;
t422 = mrSges(3,1) + m(5) * (t263 + pkin(2)) + m(4) * pkin(2) - t232 - t447;
t286 = (-t17 * t278 - t18 * t274) * qJD(6) + t449;
t344 = qJD(6) * t274;
t421 = m(7) * t286 - t89 * t343 - t88 * t344 + t429;
t282 = -pkin(8) - pkin(7);
t420 = -m(4) * pkin(7) + m(5) * t282 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t419 = Ifges(7,1) * t418 + Ifges(7,4) * t417 + Ifges(7,5) * t415;
t414 = -t135 / 0.2e1;
t413 = -t136 / 0.2e1;
t412 = t136 / 0.2e1;
t411 = -t142 / 0.2e1;
t409 = t207 / 0.2e1;
t407 = -t267 / 0.2e1;
t405 = t278 / 0.2e1;
t277 = sin(qJ(1));
t404 = pkin(1) * t277;
t403 = pkin(3) * t276;
t402 = pkin(3) * t279;
t401 = pkin(4) * t207;
t400 = pkin(4) * t261;
t247 = pkin(4) * t262;
t399 = pkin(4) * t270;
t398 = pkin(4) * t272;
t397 = pkin(5) * t240;
t281 = cos(qJ(1));
t264 = t281 * pkin(1);
t392 = t54 * mrSges(6,3);
t391 = t55 * mrSges(6,3);
t386 = Ifges(4,4) * t276;
t385 = Ifges(4,4) * t280;
t384 = Ifges(7,4) * t274;
t383 = Ifges(7,4) * t278;
t382 = pkin(3) * qJD(4);
t381 = t116 * mrSges(5,3);
t380 = t117 * mrSges(5,3);
t379 = t136 * Ifges(7,4);
t378 = t207 * Ifges(5,4);
t364 = t323 * t274;
t363 = t323 * t278;
t362 = t159 * t274;
t361 = t159 * t278;
t358 = t255 * t274;
t357 = t255 * t278;
t356 = t256 * t274;
t355 = t256 * t278;
t354 = t270 * t275;
t353 = t272 * t275;
t125 = t279 * t169 - t164;
t248 = pkin(4) + t402;
t198 = pkin(3) * t353 + t270 * t248;
t352 = t247 + t263;
t351 = qJD(1) * t276;
t350 = qJD(1) * t280;
t347 = qJD(3) * t280;
t341 = Ifges(7,5) * t52 + Ifges(7,6) * t53 + Ifges(7,3) * t81;
t252 = pkin(3) * t348;
t339 = t240 * t388;
t251 = pkin(3) * t351;
t338 = mrSges(4,3) * t351;
t337 = mrSges(4,3) * t350;
t131 = Ifges(7,4) * t135;
t64 = t136 * Ifges(7,1) + t142 * Ifges(7,5) + t131;
t334 = t64 * t405;
t333 = t247 + t320;
t331 = -t82 * mrSges(6,1) + t83 * mrSges(6,2);
t329 = -t344 / 0.2e1;
t151 = -pkin(4) * t163 + t252;
t328 = qJD(3) * t389;
t124 = -t169 * t275 - t166;
t152 = -t279 * t211 - t212 * t275;
t322 = t455 * t255;
t321 = t455 * t256;
t319 = -g(1) * t255 + g(2) * t256;
t317 = mrSges(4,1) * t276 + mrSges(4,2) * t280;
t312 = Ifges(7,1) * t278 - t384;
t311 = t280 * Ifges(4,2) + t386;
t310 = -Ifges(7,2) * t274 + t383;
t309 = Ifges(4,5) * t280 - Ifges(4,6) * t276;
t308 = Ifges(7,5) * t278 - Ifges(7,6) * t274;
t306 = -t17 * t274 + t18 * t278;
t10 = -t270 * t31 + t272 * t29;
t132 = qJ(5) * t216 + t153;
t300 = -qJ(5) * t217 + t152;
t76 = t272 * t132 + t270 * t300;
t158 = -t272 * t216 + t217 * t270;
t168 = -pkin(4) * t216 + t225;
t87 = pkin(5) * t158 - pkin(9) * t159 + t168;
t33 = t274 * t87 + t278 * t76;
t32 = -t274 * t76 + t278 * t87;
t197 = -pkin(3) * t354 + t248 * t272;
t301 = t124 - t369;
t298 = -t111 * t278 + t159 * t344;
t297 = t245 * qJD(1) * t317;
t296 = t276 * (Ifges(4,1) * t280 - t386);
t200 = t276 * t328;
t201 = t280 * t328;
t100 = -t279 * t200 - t275 * t201 - t211 * t345 - t212 * t346;
t86 = pkin(5) * t304 - pkin(9) * t323 + t401;
t223 = -t400 - t403;
t288 = m(7) * (t223 - t397) - t339;
t287 = m(7) * (-t397 - t400) - t339;
t101 = -qJD(4) * t153 + t200 * t275 - t279 * t201;
t285 = -qJ(5) * t162 - qJD(5) * t217 + t101;
t14 = t52 * Ifges(7,4) + t53 * Ifges(7,2) + t81 * Ifges(7,6);
t143 = t206 * Ifges(5,2) + t267 * Ifges(5,6) + t378;
t202 = Ifges(5,4) * t206;
t144 = t207 * Ifges(5,1) + t267 * Ifges(5,5) + t202;
t62 = t136 * Ifges(7,5) + t135 * Ifges(7,6) + t142 * Ifges(7,3);
t63 = t135 * Ifges(7,2) + t142 * Ifges(7,6) + t379;
t7 = -pkin(5) * t265 - t10;
t91 = Ifges(6,2) * t323 + t267 * Ifges(6,6) + t437;
t92 = Ifges(6,1) * t304 + t267 * Ifges(6,5) + t436;
t284 = (t334 + t433) * qJD(6) + (t364 / 0.2e1 + t329) * t63 + (-t157 * mrSges(6,2) + Ifges(6,1) * t441 + Ifges(6,5) * t407 + t308 * t411 + t310 * t414 + t312 * t413 + t392 - t433) * t323 + (t436 + t92) * t439 + (-t437 + t62) * t441 + (Ifges(5,3) + Ifges(6,3)) * t265 + (Ifges(7,5) * t274 + Ifges(7,6) * t278) * t415 + (Ifges(7,2) * t278 + t384) * t417 + (Ifges(7,1) * t274 + t383) * t418 + t274 * t419 + t14 * t405 + (Ifges(5,5) * t206 - Ifges(5,6) * t207) * t407 + t143 * t409 + t207 * t380 + t206 * t381 + (t135 * t310 + t136 * t312 + t142 * t308) * qJD(6) / 0.2e1 - (-Ifges(5,2) * t207 + t144 + t202) * t206 / 0.2e1 - t210 * (mrSges(5,1) * t207 + mrSges(5,2) * t206) + Ifges(5,5) * t137 + Ifges(5,6) * t138 + Ifges(6,6) * t82 + Ifges(6,5) * t83 + t46 * mrSges(5,1) - t45 * mrSges(5,2) - t11 * mrSges(6,2) + t10 * mrSges(6,1) - t64 * t363 / 0.2e1 + ((-t344 + t364) * t18 + (-t343 + t363) * t17 + t449) * mrSges(7,3) + t7 * t450 + t91 * t440 - t207 * (Ifges(5,1) * t206 - t378) / 0.2e1 + (-t157 * mrSges(6,1) + Ifges(7,5) * t413 - Ifges(6,2) * t439 - Ifges(6,6) * t407 + Ifges(7,6) * t414 + Ifges(7,3) * t411 + t391 + t452 - t453) * t304;
t266 = -qJ(5) + t282;
t250 = Ifges(4,4) * t350;
t244 = -pkin(5) - t398;
t231 = -qJD(3) * mrSges(4,2) + t337;
t229 = qJD(3) * mrSges(4,1) - t338;
t220 = pkin(2) + t352;
t205 = Ifges(4,1) * t351 + Ifges(4,5) * qJD(3) + t250;
t204 = Ifges(4,6) * qJD(3) + qJD(1) * t311;
t195 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t222;
t194 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t221;
t187 = -pkin(5) - t197;
t185 = -t228 * t276 + t260;
t178 = mrSges(5,1) * t267 - mrSges(5,3) * t207;
t177 = -mrSges(5,2) * t267 + mrSges(5,3) * t206;
t176 = t241 * t355 + t358;
t175 = -t241 * t356 + t357;
t174 = t251 + t401;
t173 = -t241 * t357 + t356;
t172 = t241 * t358 + t355;
t155 = -mrSges(5,1) * t206 + mrSges(5,2) * t207;
t121 = -mrSges(5,2) * t265 + mrSges(5,3) * t138;
t120 = mrSges(5,1) * t265 - mrSges(5,3) * t137;
t110 = t162 * t270 - t272 * t163;
t109 = -t199 + t125;
t98 = -mrSges(6,1) * t323 + mrSges(6,2) * t304;
t84 = t251 + t86;
t73 = -mrSges(6,2) * t265 + mrSges(6,3) * t82;
t68 = qJ(5) * t163 + qJD(5) * t216 + t100;
t59 = t272 * t109 + t270 * t301;
t57 = t105 * t272 - t368;
t56 = t105 * t270 + t102;
t44 = pkin(5) * t110 - pkin(9) * t111 + t151;
t27 = t270 * t285 + t272 * t68;
t25 = t274 * t84 + t278 * t59;
t24 = -t274 * t59 + t278 * t84;
t23 = t274 * t86 + t278 * t57;
t22 = -t274 * t57 + t278 * t86;
t5 = -qJD(6) * t33 - t27 * t274 + t278 * t44;
t4 = qJD(6) * t32 + t27 * t278 + t274 * t44;
t1 = [(t297 + t309 * qJD(3) / 0.2e1) * qJD(3) + t168 * t331 + m(7) * (t17 * t5 + t18 * t4 + t2 * t33 + t3 * t32) + m(6) * (t108 * t168 + t11 * t76 + t151 * t157 + t27 * t55) + (m(4) * t245 + t232) * t227 + m(5) * (t100 * t117 + t101 * t116 + t152 * t46 + t153 * t45 + t171 * t225 + t210 * t252) + (mrSges(5,2) * t225 + Ifges(5,1) * t217 + Ifges(5,4) * t216) * t137 + t48 * (mrSges(7,1) * t299 - mrSges(7,2) * t298) + t142 * (-Ifges(7,5) * t298 - Ifges(7,6) * t299 + Ifges(7,3) * t110) / 0.2e1 + t135 * (-Ifges(7,4) * t298 - Ifges(7,2) * t299 + Ifges(7,6) * t110) / 0.2e1 + (t17 * t298 - t18 * t299 - t2 * t362 - t3 * t361) * mrSges(7,3) + (-Ifges(7,1) * t298 - Ifges(7,4) * t299 + Ifges(7,5) * t110) * t412 + t361 * t419 + (Ifges(5,1) * t162 + Ifges(5,4) * t163) * t409 + t163 * t380 + (t435 / 0.2e1 + t434 / 0.2e1) * t265 + (t434 + t435) * t408 + (t216 * t45 - t217 * t46) * mrSges(5,3) + t280 * (Ifges(4,4) * t222 + Ifges(4,2) * t221) / 0.2e1 + (-t229 * t347 - t231 * t348 + t280 * t194 - t276 * t195 + m(4) * ((-t185 * t280 - t186 * t276) * qJD(3) + t430)) * t243 + (-t185 * t347 - t186 * t348 + t430) * mrSges(4,3) - t299 * t63 / 0.2e1 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t273 - 0.2e1 * mrSges(3,2) * t271 + m(3) * (t271 ^ 2 + t273 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (t108 * mrSges(6,2) - t10 * mrSges(6,3) + Ifges(6,1) * t83 + Ifges(6,4) * t82 + Ifges(6,5) * t426 + t308 * t415 + t310 * t417 + t312 * t418 + t7 * t313 + t64 * t329) * t159 + (-Ifges(6,4) * t83 - Ifges(6,2) * t82 + t341 / 0.2e1 + Ifges(7,3) * t415 + Ifges(7,6) * t417 + Ifges(7,5) * t418 + t108 * mrSges(6,1) - t11 * mrSges(6,3) - t426 * Ifges(6,6) + t423) * t158 + t155 * t252 + t245 * (-mrSges(4,1) * t221 + mrSges(4,2) * t222) + (-mrSges(5,1) * t225 + Ifges(5,4) * t217 + Ifges(5,2) * t216) * t138 + (t280 * (-Ifges(4,2) * t276 + t385) + t296) * t342 / 0.2e1 + (Ifges(5,5) * t162 + Ifges(6,5) * t111 + Ifges(5,6) * t163 - Ifges(6,6) * t110) * t267 / 0.2e1 + t171 * (-mrSges(5,1) * t216 + mrSges(5,2) * t217) + t210 * (-mrSges(5,1) * t163 + mrSges(5,2) * t162) + t206 * (Ifges(5,4) * t162 + Ifges(5,2) * t163) / 0.2e1 + t100 * t177 + t101 * t178 + t162 * t144 / 0.2e1 + t163 * t143 / 0.2e1 + t157 * (mrSges(6,1) * t110 + mrSges(6,2) * t111) + t151 * t98 + t152 * t120 + t153 * t121 + t27 * t139 + t111 * t92 / 0.2e1 + t110 * t62 / 0.2e1 - t110 * t91 / 0.2e1 + t4 * t88 + t5 * t89 + t76 * t73 + t32 * t20 + t33 * t21 + t323 * (Ifges(6,4) * t111 - Ifges(6,2) * t110) / 0.2e1 + t205 * t347 / 0.2e1 - t204 * t348 / 0.2e1 + (-m(6) * t10 + m(7) * t7 + t438) * (t132 * t270 - t272 * t300) + t454 * (t270 * t68 - t272 * t285) + qJDD(3) * (Ifges(4,5) * t276 + Ifges(4,6) * t280) + (Ifges(4,1) * t222 + Ifges(4,4) * t442) * t276 + (-mrSges(2,1) * t281 - t176 * mrSges(7,1) + mrSges(2,2) * t277 - t175 * mrSges(7,2) + t443 * (t256 * t220 - t255 * t266 + t264) + t427 * t264 + t420 * t255 + (-m(7) * t320 - t422) * t256) * g(2) + (mrSges(2,1) * t277 - t173 * mrSges(7,1) + mrSges(2,2) * t281 - t172 * mrSges(7,2) + t443 * (-t256 * t266 - t404) - t427 * t404 + t420 * t256 + (m(6) * t220 - m(7) * (-t220 - t320) + t422) * t255) * g(1) - t14 * t362 / 0.2e1 + t111 * t334 + t110 * t453 + (Ifges(6,1) * t111 - Ifges(6,4) * t110) * t440 + t311 * t442 - t162 * t381 + t222 * t385 / 0.2e1 - t110 * t452 - t110 * t391 - t111 * t392; m(3) * qJDD(2) + t216 * t120 + t217 * t121 + t162 * t177 + t163 * t178 + t276 * t194 + t280 * t195 + t438 * t158 - t371 * t110 + (-t229 * t276 + t231 * t280) * qJD(3) - t302 * t111 + (t73 + (-t274 * t88 - t278 * t89) * qJD(6) + t429) * t159 + (t427 + t443) * g(3) + m(4) * (t145 * t276 + t146 * t280 + (-t185 * t276 + t186 * t280) * qJD(3)) + m(7) * (t110 * t48 + t111 * t306 + t158 * t7 + t159 * t286) + m(6) * (-t10 * t158 + t11 * t159 - t110 * t54 + t111 * t55) + m(5) * (t116 * t163 + t117 * t162 + t216 * t46 + t217 * t45); (m(5) * (t275 * t45 + t279 * t46 + (-t116 * t275 + t117 * t279) * qJD(4)) - t178 * t346 + t275 * t121 + t177 * t345) * pkin(3) + (-m(7) * (t263 + t333) + t232 - m(6) * t352 - m(5) * t263 + t425) * g(3) + (t337 - t231) * t185 + (t338 + t229) * t186 - t309 * t342 / 0.2e1 + t120 * t402 + (-t17 * t24 - t18 * t25 + t187 * t7) * m(7) + (t10 * t197 + t11 * t198 - t157 * t174 - t55 * t59) * m(6) + (-t297 - t296 * qJD(1) / 0.2e1) * qJD(1) + t284 + t421 * (pkin(9) + t198) - g(1) * (t256 * t288 + t321) - g(2) * (t255 * t288 + t322) - (-Ifges(4,2) * t351 + t205 + t250) * t350 / 0.2e1 + Ifges(4,5) * t222 + Ifges(4,6) * t221 + t187 * t16 + t197 * t74 + t198 * t73 - t125 * t177 - t124 * t178 - t174 * t98 - t145 * mrSges(4,2) + t146 * mrSges(4,1) - t59 * t139 - t25 * t88 - t24 * t89 - t155 * t251 - m(5) * (t116 * t124 + t117 * t125 + t210 * t251) + Ifges(4,3) * qJDD(3) + t204 * t351 / 0.2e1 + t454 * (-t109 * t270 + t272 * t301 + (t270 * t279 + t353) * t382) + t448 * (m(5) * t403 - m(6) * t223 + t317 + t451) + (m(6) * t55 + m(7) * t306 - t302) * (t272 * t279 - t354) * t382; t74 * t398 + t73 * t399 + t284 - g(1) * (t256 * t287 + t321) - g(2) * (t255 * t287 + t322) + t244 * t16 - t116 * t177 + t117 * t178 - t57 * t139 - t23 * t88 - t22 * t89 - t98 * t401 + t371 * t56 + (-t17 * t22 - t18 * t23 + t244 * t7 - t48 * t56) * m(7) + ((t10 * t272 + t11 * t270) * pkin(4) - t157 * t401 + t54 * t56 - t55 * t57) * m(6) + (-m(6) * t247 - m(7) * t333 + t425) * g(3) + t421 * (pkin(9) + t399) + (m(6) * t400 + t451) * t448; t305 * qJD(6) + t371 * t304 + t302 * t323 + t278 * t20 + t274 * t21 + t331 + (t142 * t306 + t2 * t274 + t278 * t3 - t304 * t48 + t319) * m(7) + (t304 * t54 - t323 * t55 + t108 + t319) * m(6); -t48 * (mrSges(7,1) * t136 + mrSges(7,2) * t135) + (Ifges(7,1) * t135 - t379) * t413 + t63 * t412 + (Ifges(7,5) * t135 - Ifges(7,6) * t136) * t411 - t17 * t88 + t18 * t89 - g(1) * (mrSges(7,1) * t175 - mrSges(7,2) * t176) - g(2) * (-mrSges(7,1) * t172 + mrSges(7,2) * t173) + g(3) * t313 * t240 + (t135 * t17 + t136 * t18) * mrSges(7,3) + t341 + (-Ifges(7,2) * t136 + t131 + t64) * t414 + t423;];
tau  = t1;
