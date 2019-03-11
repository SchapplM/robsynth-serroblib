% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:43
% EndTime: 2019-03-09 04:43:23
% DurationCPUTime: 27.45s
% Computational Cost: add. (7373->682), mult. (17525->827), div. (0->0), fcn. (12664->10), ass. (0->288)
t478 = Ifges(7,4) + Ifges(6,5);
t481 = -Ifges(5,4) + t478;
t443 = Ifges(6,4) + Ifges(5,5);
t480 = -Ifges(7,5) + t443;
t216 = sin(pkin(9));
t217 = cos(pkin(9));
t220 = sin(qJ(3));
t223 = cos(qJ(3));
t457 = -t216 * t220 + t223 * t217;
t173 = t457 * qJD(3);
t180 = t216 * t223 + t217 * t220;
t126 = qJD(1) * t173 + qJDD(1) * t180;
t219 = sin(qJ(4));
t222 = cos(qJ(4));
t172 = t180 * qJD(1);
t247 = t222 * qJD(3) - t172 * t219;
t455 = qJD(4) * t247;
t81 = qJDD(3) * t219 + t126 * t222 + t455;
t403 = t81 / 0.2e1;
t139 = qJD(3) * t219 + t172 * t222;
t82 = qJD(4) * t139 - t222 * qJDD(3) + t126 * t219;
t401 = t82 / 0.2e1;
t174 = t180 * qJD(3);
t311 = qJDD(1) * t217;
t312 = qJDD(1) * t216;
t127 = -qJD(1) * t174 - t220 * t312 + t223 * t311;
t120 = qJDD(4) - t127;
t398 = t120 / 0.2e1;
t479 = -mrSges(6,2) - mrSges(5,3);
t466 = Ifges(6,2) + Ifges(5,3);
t442 = Ifges(7,2) + Ifges(6,3);
t441 = Ifges(5,6) - Ifges(6,6);
t465 = Ifges(6,6) - Ifges(7,6);
t314 = qJD(1) * qJD(2);
t199 = qJDD(1) * qJ(2) + t314;
t412 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t476 = t480 * t398 + t401 * t481 + t412 * t403;
t271 = -mrSges(3,1) * t217 + mrSges(3,2) * t216;
t475 = -m(3) * pkin(1) - mrSges(2,1) + t271;
t474 = t478 * t139;
t323 = t216 ^ 2 + t217 ^ 2;
t473 = t478 * t247;
t472 = t478 * t222;
t471 = t478 * t219;
t171 = t457 * qJD(1);
t161 = qJD(4) - t171;
t470 = m(7) * qJ(6) + mrSges(7,3) + t479;
t266 = t222 * mrSges(6,1) + t219 * mrSges(6,3);
t268 = mrSges(5,1) * t222 - mrSges(5,2) * t219;
t346 = t219 * mrSges(7,2);
t469 = -t266 - t268 - t346;
t444 = -m(7) - m(6);
t221 = sin(qJ(1));
t468 = g(2) * t221;
t464 = t120 * t465 + t442 * t82 + t478 * t81;
t438 = t161 * t465 - t247 * t442 + t474;
t463 = t139 * t443 + t161 * t466 + t247 * t441;
t349 = t172 * Ifges(4,4);
t428 = Ifges(4,6) * qJD(3);
t460 = t139 * Ifges(7,5) + t171 * Ifges(4,2) - Ifges(7,6) * t247 - t161 * Ifges(7,3) + t349 + t428;
t433 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t247 - mrSges(5,2) * t139 - mrSges(4,3) * t172;
t371 = pkin(7) + qJ(2);
t186 = t371 * t216;
t181 = qJD(1) * t186;
t187 = t371 * t217;
t182 = qJD(1) * t187;
t129 = -t181 * t220 + t182 * t223;
t459 = qJD(5) * t219 + t129;
t423 = t219 * t442 + t472;
t422 = -t219 * t441 + t222 * t443;
t317 = qJD(4) * t222;
t299 = t180 * t317;
t338 = t173 * t219;
t240 = t299 + t338;
t458 = t120 * t466 - t441 * t82 + t443 * t81;
t342 = qJ(5) * t219;
t254 = pkin(4) * t222 + t342;
t454 = mrSges(7,2) + mrSges(6,3) - mrSges(5,2);
t137 = Ifges(5,4) * t247;
t411 = t139 * t412 + t161 * t480 + t137 - t473;
t410 = m(7) * pkin(5) + mrSges(7,1);
t451 = mrSges(5,1) + mrSges(6,1) + t410;
t292 = m(3) * qJ(2) + mrSges(3,3);
t450 = mrSges(2,2) - mrSges(4,3) - t292;
t360 = Ifges(5,4) * t219;
t407 = t222 * t412 - t360 + t471;
t215 = pkin(9) + qJ(3);
t209 = sin(t215);
t210 = cos(t215);
t270 = t210 * mrSges(4,1) - t209 * mrSges(4,2);
t449 = t209 * t479 - t270;
t205 = pkin(2) * t217 + pkin(1);
t184 = -qJD(1) * t205 + qJD(2);
t100 = -pkin(3) * t171 - pkin(8) * t172 + t184;
t122 = qJD(3) * pkin(8) + t129;
t318 = qJD(4) * t219;
t183 = -qJDD(1) * t205 + qJDD(2);
t54 = -pkin(3) * t127 - pkin(8) * t126 + t183;
t281 = pkin(7) * qJDD(1) + t199;
t153 = t281 * t216;
t154 = t281 * t217;
t321 = qJD(3) * t223;
t322 = qJD(3) * t220;
t59 = -t220 * t153 + t223 * t154 - t181 * t321 - t182 * t322;
t55 = qJDD(3) * pkin(8) + t59;
t7 = -t100 * t318 - t122 * t317 - t219 * t55 + t222 * t54;
t237 = qJDD(5) - t7;
t400 = pkin(4) + pkin(5);
t1 = -qJ(6) * t81 - qJD(6) * t139 - t120 * t400 + t237;
t6 = t100 * t317 - t122 * t318 + t219 * t54 + t222 * t55;
t4 = t120 * qJ(5) + t161 * qJD(5) + t6;
t2 = qJ(6) * t82 - qJD(6) * t247 + t4;
t5 = -pkin(4) * t120 + t237;
t445 = t7 * mrSges(5,1) - t5 * mrSges(6,1) - t1 * mrSges(7,1) - t6 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3);
t376 = g(3) * t209;
t36 = mrSges(5,1) * t120 - mrSges(5,3) * t81;
t37 = -t120 * mrSges(6,1) + t81 * mrSges(6,2);
t440 = -t36 + t37;
t34 = -mrSges(6,2) * t82 + mrSges(6,3) * t120;
t39 = -mrSges(5,2) * t120 - mrSges(5,3) * t82;
t439 = t39 + t34;
t367 = mrSges(6,2) * t247;
t94 = mrSges(6,3) * t161 + t367;
t362 = mrSges(7,3) * t247;
t95 = mrSges(7,2) * t161 - t362;
t437 = t94 + t95;
t364 = mrSges(5,3) * t247;
t96 = -mrSges(5,2) * t161 + t364;
t436 = -t96 - t94;
t363 = mrSges(5,3) * t139;
t98 = mrSges(5,1) * t161 - t363;
t366 = mrSges(6,2) * t139;
t99 = -mrSges(6,1) * t161 + t366;
t435 = t99 - t98;
t341 = qJ(5) * t222;
t244 = -t219 * t400 + t341;
t434 = t161 * t244 + t459;
t340 = t171 * t219;
t370 = pkin(8) - qJ(6);
t123 = pkin(3) * t172 - pkin(8) * t171;
t128 = -t223 * t181 - t220 * t182;
t62 = t219 * t123 + t222 * t128;
t43 = t172 * qJ(5) + t62;
t432 = -qJ(6) * t340 - qJD(6) * t222 - t318 * t370 - t43;
t253 = pkin(4) * t219 - t341;
t431 = t161 * t253 - t459;
t117 = t219 * t128;
t189 = t370 * t222;
t430 = qJD(4) * t189 - qJD(6) * t219 - t117 - (-qJ(6) * t171 - t123) * t222 + t400 * t172;
t429 = Ifges(4,5) * qJD(3);
t427 = t222 * t400;
t46 = t222 * t100 - t219 * t122;
t31 = qJ(6) * t139 + t46;
t426 = qJD(5) - t31;
t425 = t174 * qJ(5) - qJD(5) * t457;
t424 = -t223 * t186 - t187 * t220;
t339 = t171 * t222;
t421 = t317 - t339;
t420 = t318 - t340;
t419 = -t219 * t7 + t222 * t6;
t418 = t219 * t5 + t222 * t4;
t224 = cos(qJ(1));
t417 = g(1) * t224 + t468;
t134 = -t186 * t220 + t187 * t223;
t415 = t134 * qJD(3);
t414 = -m(5) + t444;
t185 = -pkin(3) - t254;
t276 = -m(7) * t400 - mrSges(7,1);
t285 = pkin(3) + t342;
t409 = t470 * t210 + (m(5) * pkin(3) - m(6) * t185 + m(7) * t285 - t222 * t276 - t469) * t209;
t345 = t222 * mrSges(7,2);
t264 = -t219 * mrSges(7,1) + t345;
t344 = t222 * mrSges(6,3);
t265 = mrSges(6,1) * t219 - t344;
t267 = mrSges(5,1) * t219 + mrSges(5,2) * t222;
t273 = qJD(3) * pkin(3) + t128;
t233 = qJ(5) * t139 + t273;
t30 = t247 * t400 + qJD(6) + t233;
t48 = -pkin(4) * t247 - t233;
t408 = t30 * t264 + t48 * t265 - t267 * t273;
t406 = pkin(8) * t414;
t402 = -t82 / 0.2e1;
t399 = -t120 / 0.2e1;
t397 = t247 / 0.2e1;
t396 = -t247 / 0.2e1;
t395 = -t139 / 0.2e1;
t394 = t139 / 0.2e1;
t393 = -t161 / 0.2e1;
t392 = t161 / 0.2e1;
t391 = -t171 / 0.2e1;
t390 = -t172 / 0.2e1;
t389 = t172 / 0.2e1;
t201 = t209 * pkin(8);
t202 = t210 * pkin(3);
t85 = -mrSges(6,1) * t247 - mrSges(6,3) * t139;
t86 = mrSges(7,1) * t247 + mrSges(7,2) * t139;
t369 = t85 - t86;
t361 = mrSges(7,3) * t139;
t359 = Ifges(5,4) * t222;
t354 = t128 * mrSges(4,3);
t353 = t129 * mrSges(4,3);
t350 = t139 * Ifges(5,4);
t47 = t219 * t100 + t222 * t122;
t343 = qJ(5) * t247;
t337 = t173 * t222;
t336 = t180 * t219;
t333 = t209 * t224;
t332 = t210 * t224;
t330 = t371 * t224;
t329 = t219 * t221;
t328 = t221 * t222;
t327 = t222 * t224;
t325 = t224 * t219;
t125 = -pkin(3) * t457 - pkin(8) * t180 - t205;
t80 = t219 * t125 + t222 * t134;
t324 = t202 + t201;
t315 = qJD(5) * t222;
t101 = qJD(2) * t457 + qJD(3) * t424;
t304 = t219 * t101 + t125 * t318 + t134 * t317;
t124 = pkin(3) * t174 - pkin(8) * t173;
t303 = t222 * t101 + t219 * t124 + t125 * t317;
t50 = -qJ(5) * t457 + t80;
t300 = t180 * t318;
t24 = -t82 * mrSges(7,1) + t81 * mrSges(7,2);
t288 = -t318 / 0.2e1;
t286 = t317 / 0.2e1;
t35 = -t120 * mrSges(7,1) - t81 * mrSges(7,3);
t284 = -t205 - t202;
t283 = -t127 * mrSges(4,1) + t126 * mrSges(4,2);
t61 = t123 * t222 - t117;
t131 = t219 * t134;
t79 = t125 * t222 - t131;
t278 = t224 * t205 + t221 * t371;
t60 = -t223 * t153 - t220 * t154 + t181 * t322 - t182 * t321;
t277 = t210 * t254 + t324;
t32 = -qJ(6) * t247 + t47;
t272 = -mrSges(3,1) * t311 + mrSges(3,2) * t312;
t260 = -Ifges(5,2) * t219 + t359;
t255 = Ifges(7,5) * t222 + Ifges(7,6) * t219;
t248 = -qJ(6) * t173 - qJD(6) * t180;
t21 = t124 * t222 - t304;
t245 = pkin(3) * t332 + pkin(8) * t333 + t278;
t239 = t300 - t337;
t20 = -t134 * t318 + t303;
t236 = qJDD(3) * pkin(3) + t60;
t234 = Ifges(7,5) * t81 + Ifges(7,6) * t82 - Ifges(7,3) * t120;
t162 = t210 * t329 + t327;
t164 = t210 * t325 - t328;
t232 = -g(1) * t164 - g(2) * t162 - t219 * t376;
t229 = qJ(5) * t81 + qJD(5) * t139 + t236;
t102 = t180 * qJD(2) + t415;
t208 = -qJDD(1) * pkin(1) + qJDD(2);
t190 = t209 * t341;
t188 = t370 * t219;
t177 = t285 + t427;
t165 = t210 * t327 + t329;
t163 = t210 * t328 - t325;
t157 = Ifges(4,4) * t171;
t148 = t161 * qJ(5);
t141 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t171;
t107 = t172 * Ifges(4,1) + t157 + t429;
t97 = -mrSges(7,1) * t161 - t361;
t84 = pkin(4) * t139 - t343;
t83 = t180 * t253 - t424;
t69 = Ifges(5,2) * t247 + t161 * Ifges(5,6) + t350;
t58 = -t139 * t400 + t343;
t57 = t180 * t244 + t424;
t51 = pkin(4) * t457 - t79;
t45 = -pkin(4) * t172 - t61;
t42 = qJ(6) * t336 + t50;
t41 = t148 + t47;
t40 = -pkin(4) * t161 + qJD(5) - t46;
t38 = mrSges(7,2) * t120 + mrSges(7,3) * t82;
t29 = t131 + (-qJ(6) * t180 - t125) * t222 + t400 * t457;
t27 = t148 + t32;
t26 = t253 * t173 + (qJD(4) * t254 - t315) * t180 + t102;
t25 = mrSges(5,1) * t82 + mrSges(5,2) * t81;
t23 = mrSges(6,1) * t82 - mrSges(6,3) * t81;
t22 = -t161 * t400 + t426;
t19 = -pkin(4) * t174 - t21;
t18 = -t415 + t244 * t173 + (-qJD(2) + t315 + (-t342 - t427) * qJD(4)) * t180;
t14 = t81 * Ifges(5,4) - t82 * Ifges(5,2) + t120 * Ifges(5,6);
t11 = t20 + t425;
t10 = qJ(6) * t299 + (-qJD(4) * t134 - t248) * t219 + t303 + t425;
t9 = pkin(4) * t82 - t229;
t8 = qJ(6) * t300 - t400 * t174 + (-t124 + t248) * t222 + t304;
t3 = -t400 * t82 + qJDD(6) + t229;
t12 = [-(-t234 / 0.2e1 - Ifges(4,6) * qJDD(3) - Ifges(4,2) * t127 - Ifges(4,4) * t126 - Ifges(7,3) * t399 + Ifges(5,6) * t402 + t183 * mrSges(4,1) - t59 * mrSges(4,3) + t465 * t401 + t466 * t398 + t480 * t403 + t445 + t458 / 0.2e1) * t457 + (t174 * t465 - t239 * t478 + t240 * t442) * t396 + (t174 * t466 - t239 * t443 - t240 * t441) * t392 - (-m(4) * t60 - m(5) * t236 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t126 + t25) * t424 + (t480 * t174 - t412 * t239 + t240 * t481) * t394 - t273 * (mrSges(5,1) * t240 - mrSges(5,2) * t239) + (-m(4) * t128 - m(5) * t273 - t433) * t102 - t460 * t174 / 0.2e1 + m(3) * (-pkin(1) * t208 + (t199 + t314) * qJ(2) * t323) + (t183 * mrSges(4,2) - t60 * mrSges(4,3) + Ifges(4,1) * t126 + Ifges(4,4) * t127 + Ifges(4,5) * qJDD(3) - t267 * t236 + t255 * t399 + t260 * t402 + t264 * t3 + t265 * t9 + t286 * t438 + t288 * t411 + t398 * t422 + t401 * t423 + t403 * t407 + (mrSges(6,2) * t5 - mrSges(5,3) * t7 - mrSges(7,3) * t1 + t476) * t222) * t180 + t463 * t174 / 0.2e1 + t438 * t338 / 0.2e1 - t174 * t353 - t173 * t354 + m(5) * (t20 * t47 + t21 * t46 + t6 * t80 + t7 * t79) + m(4) * (t101 * t129 + t134 * t59 - t183 * t205) + (Ifges(3,4) * t216 + Ifges(3,2) * t217) * t311 + (Ifges(3,1) * t216 + Ifges(3,4) * t217) * t312 + (-Ifges(5,4) * t239 - Ifges(5,2) * t240 + Ifges(5,6) * t174) * t397 + m(7) * (t1 * t29 + t10 * t27 + t18 * t30 + t2 * t42 + t22 * t8 + t3 * t57) + m(6) * (t11 * t41 + t19 * t40 + t26 * t48 + t4 * t50 + t5 * t51 + t83 * t9) + (Ifges(4,1) * t173 - Ifges(4,4) * t174) * t389 + (-Ifges(7,5) * t239 + Ifges(7,6) * t240 - Ifges(7,3) * t174) * t393 - t240 * t69 / 0.2e1 + (t464 / 0.2e1 - t14 / 0.2e1 - t6 * mrSges(5,3) - t4 * mrSges(6,2) + t2 * mrSges(7,3)) * t336 + 0.2e1 * t323 * t199 * mrSges(3,3) + t411 * t337 / 0.2e1 + t184 * (mrSges(4,1) * t174 + mrSges(4,2) * t173) + t171 * (Ifges(4,4) * t173 - Ifges(4,2) * t174) / 0.2e1 + qJD(3) * (Ifges(4,5) * t173 - Ifges(4,6) * t174) / 0.2e1 + t173 * t107 / 0.2e1 + Ifges(2,3) * qJDD(1) + t101 * t141 + t134 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t127) + t19 * t99 + t11 * t94 + t10 * t95 + t20 * t96 + t8 * t97 + t21 * t98 + t83 * t23 + t26 * t85 + t18 * t86 + t79 * t36 + t80 * t39 + t57 * t24 + t50 * t34 + t51 * t37 + t42 * t38 + t29 * t35 + t46 * (mrSges(5,1) * t174 + mrSges(5,3) * t239) + t22 * (-mrSges(7,1) * t174 + mrSges(7,3) * t239) + t40 * (-mrSges(6,1) * t174 - mrSges(6,2) * t239) + t27 * (mrSges(7,2) * t174 + mrSges(7,3) * t240) + t48 * (mrSges(6,1) * t240 + mrSges(6,3) * t239) + t41 * (-mrSges(6,2) * t240 + mrSges(6,3) * t174) + t47 * (-mrSges(5,2) * t174 - mrSges(5,3) * t240) + t30 * (-mrSges(7,1) * t240 - mrSges(7,2) * t239) + (-m(5) * t330 + t444 * (-t163 * pkin(4) - qJ(5) * t162 + t330) + (-m(4) * t371 + t450) * t224 + t451 * t163 + t454 * t162 + (-m(7) * t284 - (-m(7) * t370 + mrSges(7,3)) * t209 + m(4) * t205 + (-m(5) - m(6)) * (t284 - t201) - t449 - t475) * t221) * g(1) + (-m(4) * t278 - m(5) * t245 + t444 * (t165 * pkin(4) + qJ(5) * t164 + t245) + t470 * t333 + (-t270 + t475) * t224 + t450 * t221 - t451 * t165 - t454 * t164) * g(2) - t205 * t283 + t208 * t271 - pkin(1) * t272; -t171 * t141 + (-t369 + t433) * t172 + (-t35 + t161 * (t96 + t437) - t440) * t222 + (t38 + t161 * (t97 + t435) + t439) * t219 + m(3) * t208 + t272 + t283 + (-g(1) * t221 + g(2) * t224) * (m(3) + m(4) - t414) - t292 * t323 * qJD(1) ^ 2 + (-t1 * t222 + t172 * t30 + t2 * t219 + t161 * (t219 * t22 + t222 * t27)) * m(7) + (-t172 * t48 + t219 * t4 - t222 * t5 + t161 * (t219 * t40 + t222 * t41)) * m(6) + (t273 * t172 + t219 * t6 + t222 * t7 + t161 * (-t219 * t46 + t222 * t47)) * m(5) + (t128 * t172 - t129 * t171 + t183) * m(4); (t22 * mrSges(7,1) - t46 * mrSges(5,1) + t40 * mrSges(6,1) + t47 * mrSges(5,2) + Ifges(5,6) * t396 - Ifges(4,2) * t391 - Ifges(7,3) * t392 + t353 - t41 * mrSges(6,3) - t27 * mrSges(7,2) - t184 * mrSges(4,1) + t428 / 0.2e1 + t465 * t397 + t466 * t393 + t480 * t395) * t172 + (-t423 / 0.2e1 + t260 / 0.2e1) * t455 + (t210 * t406 + t409) * t468 + (t260 * t396 + Ifges(4,1) * t390 + t255 * t392 + t354 - t184 * mrSges(4,2) - t429 / 0.2e1 + t423 * t397 + t422 * t393 + t407 * t395 - t408) * t171 + t430 * t97 + t431 * t85 + t432 * t95 + t433 * t129 + (t1 * t188 + t177 * t3 + t189 * t2 + t22 * t430 + t27 * t432 + t30 * t434) * m(7) + t434 * t86 + (t286 - t339 / 0.2e1) * t411 - t464 * t222 / 0.2e1 + (t318 / 0.2e1 - t340 / 0.2e1) * t438 + (t394 * t407 + t408) * qJD(4) + (t422 / 0.2e1 - t255 / 0.2e1) * qJD(4) * t161 + t236 * t268 + (pkin(3) * t236 + t129 * t273 - t46 * t61 - t47 * t62) * m(5) + (t185 * t9 - t40 * t45 - t41 * t43 + t431 * t48) * m(6) + (((-t219 * t41 + t222 * t40) * qJD(4) + t418) * m(6) + t439 * t222 + t440 * t219 + t435 * t317 + t436 * t318 + ((-t219 * t47 - t222 * t46) * qJD(4) + t419) * m(5)) * pkin(8) + t460 * t389 + t3 * (t222 * mrSges(7,1) + t346) + (t107 + t157) * t391 + (-t349 + t463) * t390 + (-m(5) * t324 - m(6) * t277 - m(7) * (-qJ(6) * t209 + t277) + (-t222 * t410 + t469) * t210 + t449) * g(3) + t417 * (mrSges(4,1) * t209 + mrSges(4,2) * t210) + (t40 * t421 - t41 * t420 + t418) * mrSges(6,2) + (-t420 * t47 - t421 * t46 + t419) * mrSges(5,3) + (Ifges(7,5) * t219 - Ifges(7,6) * t222) * t399 + (Ifges(5,2) * t222 + t360) * t402 + (t224 * t409 + t332 * t406) * g(1) + t185 * t23 + t188 * t35 + t189 * t38 + t177 * t24 - t128 * t141 + Ifges(4,5) * t126 + Ifges(4,6) * t127 - t45 * t99 - t43 * t94 - t62 * t96 - t61 * t98 - t59 * mrSges(4,2) + t60 * mrSges(4,1) - pkin(3) * t25 + Ifges(4,3) * qJDD(3) + (-t222 * t442 + t471) * t401 + (t219 * t412 + t359 - t472) * t403 + (t219 * t443 + t222 * t441) * t398 + (-t1 * t219 - t2 * t222 - t22 * t421 + t27 * t420 + t376) * mrSges(7,3) + t219 * t476 + (t340 / 0.2e1 + t288) * t69 + t222 * t14 / 0.2e1 - t9 * t266; t445 + t273 * (mrSges(5,1) * t139 + mrSges(5,2) * t247) + (Ifges(7,5) * t247 + Ifges(7,6) * t139) * t392 - t48 * (mrSges(6,1) * t139 - mrSges(6,3) * t247) + (-m(6) * t190 + (-t344 - t345 + (m(6) * pkin(4) + mrSges(6,1) - t276) * t219) * t209) * g(3) + (-pkin(4) * t5 + qJ(5) * t4 + qJD(5) * t41 - t48 * t84) * m(6) + (t34 + t38) * qJ(5) + t22 * t362 - t27 * t361 + (t444 * (-t164 * pkin(4) + qJ(5) * t165) - t454 * t165 + t451 * t164) * g(1) + (t444 * (-t162 * pkin(4) + qJ(5) * t163) - t454 * t163 + t451 * t162) * g(2) + t69 * t394 + t41 * t366 - t40 * t367 + t267 * t376 + t458 - t31 * t95 - t32 * t97 - t84 * t85 - t58 * t86 - pkin(4) * t37 + (t139 * t442 + t473) * t397 + (-m(6) * t40 + t363 - t435) * t47 + (-m(6) * t41 + t364 + t436) * t46 + (t247 * t412 - t350 + t438 + t474) * t395 + t437 * qJD(5) - t234 - t30 * (-mrSges(7,1) * t139 + mrSges(7,2) * t247) + (-t139 * t441 + t247 * t443) * t393 + (-g(3) * t190 + qJ(5) * t2 - t1 * t400 - t22 * t32 + t27 * t426 - t30 * t58) * m(7) - t400 * t35 + (-Ifges(5,2) * t139 + t137 + t411) * t396; t369 * t139 - t437 * t161 + t35 + t37 + (-t139 * t30 - t161 * t27 + t1 + t232) * m(7) + (t139 * t48 - t161 * t41 + t232 + t5) * m(6); t247 * t95 + t139 * t97 + (-g(3) * t210 + t22 * t139 + t209 * t417 + t247 * t27 + t3) * m(7) + t24;];
tau  = t12;
