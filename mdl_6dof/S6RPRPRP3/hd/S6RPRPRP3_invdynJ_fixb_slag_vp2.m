% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:49
% EndTime: 2019-03-09 03:08:21
% DurationCPUTime: 19.92s
% Computational Cost: add. (6700->629), mult. (14145->798), div. (0->0), fcn. (9437->14), ass. (0->283)
t411 = mrSges(6,1) + mrSges(7,1);
t410 = mrSges(6,2) - mrSges(7,3);
t232 = cos(qJ(3));
t304 = qJD(1) * t232;
t202 = qJD(5) - t304;
t350 = t202 / 0.2e1;
t405 = -Ifges(6,3) - Ifges(7,2);
t435 = t405 * t350;
t229 = sin(qJ(3));
t223 = sin(pkin(10));
t225 = cos(pkin(10));
t271 = -mrSges(5,1) * t225 + mrSges(5,2) * t223;
t247 = m(5) * pkin(3) - t271;
t430 = -mrSges(6,3) - mrSges(7,2);
t434 = t430 * t229 - t232 * t247;
t305 = qJD(1) * t229;
t168 = qJD(3) * t223 + t225 * t305;
t228 = sin(qJ(5));
t231 = cos(qJ(5));
t250 = -t225 * qJD(3) + t223 * t305;
t113 = t228 * t168 + t231 * t250;
t295 = qJD(1) * qJD(3);
t180 = qJDD(1) * t229 + t232 * t295;
t138 = qJDD(3) * t225 - t180 * t223;
t139 = qJDD(3) * t223 + t180 * t225;
t47 = -qJD(5) * t113 + t228 * t138 + t231 * t139;
t433 = -t47 / 0.2e1;
t240 = t231 * t168 - t228 * t250;
t48 = qJD(5) * t240 - t231 * t138 + t228 * t139;
t369 = -t48 / 0.2e1;
t363 = -t113 / 0.2e1;
t407 = Ifges(7,4) + Ifges(6,5);
t432 = -t407 / 0.2e1;
t412 = qJD(3) / 0.2e1;
t208 = pkin(4) * t225 + pkin(3);
t198 = t232 * t208;
t284 = m(5) * qJ(4) + mrSges(5,3);
t340 = pkin(8) + qJ(4);
t306 = t340 * t229;
t417 = -m(6) - m(7);
t431 = t417 * (t198 + t306) - t229 * t284 + t434;
t409 = Ifges(6,1) + Ifges(7,1);
t408 = -Ifges(6,4) + Ifges(7,5);
t406 = -Ifges(6,6) + Ifges(7,6);
t172 = t223 * t231 + t225 * t228;
t249 = t172 * t232;
t136 = qJD(1) * t249;
t157 = t172 * qJD(5);
t395 = t136 - t157;
t311 = t223 * t228;
t171 = -t231 * t225 + t311;
t248 = t171 * t232;
t137 = qJD(1) * t248;
t156 = t171 * qJD(5);
t394 = -t137 + t156;
t221 = pkin(10) + qJ(5);
t213 = sin(t221);
t215 = cos(t221);
t429 = -t410 * t213 + t411 * t215;
t224 = sin(pkin(9));
t207 = pkin(1) * t224 + pkin(7);
t190 = t207 * qJD(1);
t149 = t229 * qJD(2) + t232 * t190;
t133 = qJD(3) * qJ(4) + t149;
t278 = -qJ(4) * t229 - pkin(2);
t226 = cos(pkin(9));
t348 = pkin(1) * t226;
t164 = -pkin(3) * t232 + t278 - t348;
t140 = t164 * qJD(1);
t70 = -t133 * t223 + t225 * t140;
t60 = -pkin(4) * t304 - pkin(8) * t168 + t70;
t71 = t225 * t133 + t223 * t140;
t62 = -pkin(8) * t250 + t71;
t18 = -t228 * t62 + t231 * t60;
t16 = -pkin(5) * t202 + qJD(6) - t18;
t19 = t228 * t60 + t231 * t62;
t17 = qJ(6) * t202 + t19;
t332 = Ifges(4,4) * t229;
t266 = t232 * Ifges(4,2) + t332;
t428 = -t17 * mrSges(7,3) - t18 * mrSges(6,1) + Ifges(4,6) * t412 + qJD(1) * t266 / 0.2e1 + t16 * mrSges(7,1) + t19 * mrSges(6,2) - Ifges(5,5) * t168 / 0.2e1 + Ifges(5,6) * t250 / 0.2e1 + Ifges(5,3) * t304 / 0.2e1 + t406 * t363 + t435 + t240 * t432;
t179 = -t232 * qJDD(1) + t229 * t295;
t170 = qJDD(5) + t179;
t175 = t229 * t190;
t423 = qJD(2) * qJD(3) + t207 * qJDD(1);
t289 = t229 * qJDD(2) + t232 * t423;
t81 = qJDD(3) * qJ(4) + (qJD(4) - t175) * qJD(3) + t289;
t209 = -pkin(2) - t348;
t189 = t209 * qJDD(1);
t301 = qJD(4) * t229;
t91 = pkin(3) * t179 - qJ(4) * t180 - qJD(1) * t301 + t189;
t34 = -t223 * t81 + t225 * t91;
t21 = pkin(4) * t179 - pkin(8) * t139 + t34;
t35 = t223 * t91 + t225 * t81;
t24 = pkin(8) * t138 + t35;
t4 = -qJD(5) * t19 + t21 * t231 - t228 * t24;
t2 = -pkin(5) * t170 + qJDD(6) - t4;
t354 = t170 / 0.2e1;
t368 = t48 / 0.2e1;
t370 = t47 / 0.2e1;
t302 = qJD(3) * t232;
t97 = qJDD(2) * t232 - t190 * t302 - t229 * t423;
t89 = -qJDD(3) * pkin(3) + qJDD(4) - t97;
t61 = -pkin(4) * t138 + t89;
t5 = pkin(5) * t48 - qJ(6) * t47 - qJD(6) * t240 + t61;
t427 = -mrSges(6,2) * t61 - mrSges(7,2) * t2 + mrSges(6,3) * t4 + mrSges(7,3) * t5 - Ifges(7,5) * t368 + t170 * t432 - t354 * t407 + (-t370 + t433) * t409 + (-Ifges(6,4) + t408) * t369;
t360 = t240 / 0.2e1;
t362 = t113 / 0.2e1;
t148 = t232 * qJD(2) - t175;
t131 = -qJD(3) * pkin(3) + qJD(4) - t148;
t101 = pkin(4) * t250 + t131;
t27 = t113 * pkin(5) - qJ(6) * t240 + t101;
t108 = Ifges(7,5) * t240;
t49 = t202 * Ifges(7,6) + t113 * Ifges(7,3) + t108;
t328 = Ifges(6,4) * t240;
t52 = -t113 * Ifges(6,2) + t202 * Ifges(6,6) + t328;
t421 = t27 * mrSges(7,1) + t49 / 0.2e1 - t52 / 0.2e1;
t426 = -Ifges(6,2) * t363 + Ifges(7,3) * t362 + t350 * t406 + t408 * t360 + t421;
t414 = mrSges(7,3) * t27;
t425 = Ifges(6,4) * t363 + Ifges(7,5) * t362 + t407 * t350 + t409 * t360 - t414;
t109 = Ifges(6,4) * t113;
t327 = Ifges(7,5) * t113;
t401 = t202 * t407 + t240 * t409 - t109 + t327;
t424 = -t401 / 0.2e1;
t222 = qJ(1) + pkin(9);
t214 = sin(t222);
t216 = cos(t222);
t422 = g(1) * t216 + g(2) * t214;
t418 = -m(4) - m(5);
t358 = t138 / 0.2e1;
t357 = t139 / 0.2e1;
t416 = -t179 / 0.2e1;
t352 = t179 / 0.2e1;
t415 = t180 / 0.2e1;
t30 = mrSges(6,1) * t170 - mrSges(6,3) * t47;
t31 = -t170 * mrSges(7,1) + t47 * mrSges(7,2);
t403 = t31 - t30;
t32 = -mrSges(6,2) * t170 - mrSges(6,3) * t48;
t33 = -mrSges(7,2) * t48 + mrSges(7,3) * t170;
t402 = t32 + t33;
t287 = t223 * t304;
t121 = pkin(4) * t287 + t149;
t400 = -pkin(5) * t395 + qJ(6) * t394 - qJD(6) * t172 - t121;
t334 = mrSges(6,3) * t113;
t85 = -mrSges(6,2) * t202 - t334;
t336 = mrSges(7,2) * t113;
t88 = mrSges(7,3) * t202 - t336;
t339 = t85 + t88;
t333 = mrSges(6,3) * t240;
t86 = mrSges(6,1) * t202 - t333;
t335 = mrSges(7,2) * t240;
t87 = -mrSges(7,1) * t202 + t335;
t338 = -t87 + t86;
t273 = mrSges(4,1) * t232 - mrSges(4,2) * t229;
t399 = -mrSges(3,1) - t273;
t291 = mrSges(4,3) * t305;
t396 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t250 + t168 * mrSges(5,2) + t291;
t329 = Ifges(5,4) * t225;
t265 = -Ifges(5,2) * t223 + t329;
t330 = Ifges(5,4) * t223;
t267 = Ifges(5,1) * t225 - t330;
t393 = t168 * (Ifges(5,5) * t229 + t232 * t267) - (Ifges(5,6) * t229 + t232 * t265) * t250;
t392 = -t170 * t405 + t406 * t48 + t407 * t47;
t134 = mrSges(5,2) * t304 - mrSges(5,3) * t250;
t135 = -mrSges(5,1) * t304 - mrSges(5,3) * t168;
t391 = t134 * t225 - t135 * t223;
t102 = -mrSges(5,2) * t179 + mrSges(5,3) * t138;
t103 = mrSges(5,1) * t179 - mrSges(5,3) * t139;
t390 = t102 * t225 - t103 * t223;
t303 = qJD(3) * t229;
t96 = -t190 * t303 + t289;
t389 = -t229 * t97 + t232 * t96;
t259 = -t223 * t34 + t225 * t35;
t14 = t48 * mrSges(7,1) - t47 * mrSges(7,3);
t15 = t48 * mrSges(6,1) + t47 * mrSges(6,2);
t78 = -t138 * mrSges(5,1) + t139 * mrSges(5,2);
t387 = -t14 - t15 - t78;
t270 = mrSges(5,1) * t223 + mrSges(5,2) * t225;
t386 = mrSges(3,2) - mrSges(4,3) - t270;
t212 = Ifges(4,4) * t304;
t384 = t225 * (Ifges(5,1) * t168 - Ifges(5,4) * t250 - Ifges(5,5) * t304) + Ifges(4,1) * t305 + Ifges(4,5) * qJD(3) + t212;
t147 = t225 * t164;
t308 = t225 * t229;
t92 = -pkin(8) * t308 + t147 + (-t207 * t223 - pkin(4)) * t232;
t307 = t225 * t232;
t117 = t223 * t164 + t207 * t307;
t310 = t223 * t229;
t98 = -pkin(8) * t310 + t117;
t337 = t228 * t92 + t231 * t98;
t262 = pkin(3) * t229 - qJ(4) * t232;
t153 = qJD(3) * t262 - t301;
t286 = t207 * t303;
t110 = t225 * t153 + t223 * t286;
t253 = pkin(4) * t229 - pkin(8) * t307;
t74 = qJD(3) * t253 + t110;
t141 = t223 * t153;
t197 = t229 * t207;
t309 = t223 * t232;
t90 = t141 + (-pkin(8) * t309 - t197 * t225) * qJD(3);
t13 = -qJD(5) * t337 - t228 * t90 + t231 * t74;
t383 = m(7) * pkin(5) + t411;
t272 = mrSges(4,1) * t229 + mrSges(4,2) * t232;
t382 = -t131 * t232 * t270 - t71 * (-mrSges(5,2) * t229 - mrSges(5,3) * t309) - t70 * (mrSges(5,1) * t229 - mrSges(5,3) * t307) - t209 * qJD(1) * t272;
t380 = m(7) * qJ(6) - t410;
t298 = qJD(5) * t231;
t300 = qJD(5) * t228;
t3 = t228 * t21 + t231 * t24 + t60 * t298 - t300 * t62;
t1 = qJ(6) * t170 + qJD(6) * t202 + t3;
t378 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t374 = mrSges(6,1) * t61 + mrSges(7,1) * t5 + 0.2e1 * Ifges(7,3) * t368 + Ifges(6,4) * t433 - t170 * Ifges(6,6) / 0.2e1 + (t408 + Ifges(7,5)) * t370 + (t406 + Ifges(7,6)) * t354 + (-t369 + t368) * Ifges(6,2);
t365 = Ifges(5,1) * t357 + Ifges(5,4) * t358 + Ifges(5,5) * t352;
t361 = -t240 / 0.2e1;
t351 = -t202 / 0.2e1;
t230 = sin(qJ(1));
t347 = pkin(1) * t230;
t346 = pkin(4) * t223;
t343 = g(3) * t229;
t233 = cos(qJ(1));
t220 = t233 * pkin(1);
t178 = t262 * qJD(1);
t99 = -t148 * t223 + t225 * t178;
t67 = qJD(1) * t253 + t99;
t100 = t225 * t148 + t223 * t178;
t76 = -pkin(8) * t287 + t100;
t29 = t228 * t67 + t231 * t76;
t331 = Ifges(4,4) * t232;
t322 = t229 * t89;
t314 = t214 * t232;
t312 = t216 * t232;
t182 = t207 * t302;
t285 = t223 * t302;
t143 = pkin(4) * t285 + t182;
t152 = pkin(4) * t310 + t197;
t299 = qJD(5) * t229;
t296 = m(5) - t417;
t290 = mrSges(4,3) * t304;
t288 = t216 * pkin(2) + t214 * pkin(7) + t220;
t282 = -t304 / 0.2e1;
t279 = t216 * pkin(7) - t347;
t276 = -t295 / 0.2e1;
t275 = t295 / 0.2e1;
t264 = Ifges(4,5) * t232 - Ifges(4,6) * t229;
t263 = Ifges(5,5) * t225 - Ifges(5,6) * t223;
t261 = pkin(5) * t215 + qJ(6) * t213;
t258 = -t223 * t70 + t225 * t71;
t28 = -t228 * t76 + t231 * t67;
t38 = -t228 * t98 + t231 * t92;
t186 = t340 * t223;
t187 = t340 * t225;
t255 = -t231 * t186 - t187 * t228;
t120 = -t186 * t228 + t187 * t231;
t12 = t228 * t74 + t231 * t90 + t92 * t298 - t300 * t98;
t251 = t229 * (Ifges(4,1) * t232 - t332);
t241 = t232 * (Ifges(5,3) * t229 + t232 * t263);
t193 = -qJD(3) * mrSges(4,2) + t290;
t151 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t180;
t150 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t179;
t145 = t171 * t229;
t144 = t172 * t229;
t129 = t213 * t214 + t215 * t312;
t128 = t213 * t312 - t214 * t215;
t127 = -t216 * t213 + t215 * t314;
t126 = t213 * t314 + t215 * t216;
t116 = -t207 * t309 + t147;
t111 = -t225 * t286 + t141;
t107 = pkin(5) * t171 - qJ(6) * t172 - t208;
t105 = Ifges(5,4) * t168 - Ifges(5,2) * t250 - Ifges(5,6) * t304;
t95 = qJD(3) * t249 + t298 * t308 - t299 * t311;
t94 = -qJD(3) * t248 - t172 * t299;
t84 = qJD(4) * t172 + qJD(5) * t120;
t83 = -qJD(4) * t171 + qJD(5) * t255;
t65 = pkin(5) * t144 + qJ(6) * t145 + t152;
t63 = t139 * Ifges(5,4) + t138 * Ifges(5,2) + t179 * Ifges(5,6);
t59 = mrSges(6,1) * t113 + mrSges(6,2) * t240;
t58 = mrSges(7,1) * t113 - mrSges(7,3) * t240;
t57 = pkin(5) * t240 + qJ(6) * t113;
t37 = pkin(5) * t232 - t38;
t36 = -qJ(6) * t232 + t337;
t26 = -pkin(5) * t305 - t28;
t25 = qJ(6) * t305 + t29;
t22 = pkin(5) * t95 - qJ(6) * t94 + qJD(6) * t145 + t143;
t7 = -pkin(5) * t303 - t13;
t6 = qJ(6) * t303 - qJD(6) * t232 + t12;
t8 = [(-m(3) * t220 - mrSges(2,1) * t233 + mrSges(2,2) * t230 + t418 * t288 + t417 * (t214 * t346 + t288) - t383 * t129 - t380 * t128 + t386 * t214 + (t399 + t431) * t216) * g(2) - t189 * t273 + t308 * t365 + t393 * t412 + t331 * t415 + m(7) * (t1 * t36 + t16 * t7 + t17 * t6 + t2 * t37 + t22 * t27 + t5 * t65) + (-t148 * t302 + t389) * mrSges(4,3) + t396 * t182 + (-t308 * t34 - t310 * t35) * mrSges(5,3) + (t401 / 0.2e1 + mrSges(6,2) * t101 + mrSges(7,2) * t16 - mrSges(6,3) * t18 + t425) * t94 + (t101 * mrSges(6,1) - t17 * mrSges(7,2) - t19 * mrSges(6,3) + t426) * t95 + (-t149 * mrSges(4,3) + Ifges(6,6) * t363 + Ifges(7,6) * t362 + t407 * t360 - t428 - t435) * t303 + (-Ifges(5,5) * t357 - Ifges(5,6) * t358 + (-Ifges(4,2) * t229 + t331) * t275 - Ifges(5,3) * t352 + Ifges(4,4) * t415 + Ifges(4,2) * t416 + t207 * t150 + t35 * mrSges(5,2) - t34 * mrSges(5,1) - Ifges(7,6) * t368 - Ifges(6,6) * t369 - t407 * t370 + t405 * t354 + t378) * t232 + (Ifges(4,1) * t180 + Ifges(4,4) * t416 + t263 * t352 + t265 * t358 + t267 * t357) * t229 + m(4) * (t189 * t209 + t389 * t207) - (Ifges(5,5) * t139 + Ifges(5,6) * t138 + Ifges(5,3) * t179 + t392) * t232 / 0.2e1 + t427 * t145 + (m(3) * t347 + mrSges(2,1) * t230 + mrSges(2,2) * t233 + t418 * t279 + t417 * (-t214 * t306 + t216 * t346 + t279) + t383 * t127 + t380 * t126 + t386 * t216 + (m(4) * pkin(2) - m(5) * t278 + t229 * mrSges(5,3) + t417 * (-pkin(2) - t198) - t399 - t434) * t214) * g(1) + t270 * t322 + (m(4) * (-t148 * t232 - t149 * t229) * t207 + t264 * t412 - t382) * qJD(3) + m(6) * (t101 * t143 + t12 * t19 + t13 * t18 + t152 * t61 + t3 * t337 + t38 * t4) + t337 * t32 + t384 * t302 / 0.2e1 + (-mrSges(7,2) * t1 - mrSges(6,3) * t3 + t374) * t144 + (t78 - t151) * t197 + qJDD(3) * (Ifges(4,5) * t229 + Ifges(4,6) * t232) + t209 * (mrSges(4,1) * t179 + mrSges(4,2) * t180) + t152 * t15 + t143 * t59 + t111 * t134 + t110 * t135 + t116 * t103 + t117 * t102 + (0.2e1 * (mrSges(3,1) * t226 - mrSges(3,2) * t224) * pkin(1) + m(3) * (t224 ^ 2 + t226 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + m(5) * (t110 * t70 + t111 * t71 + t116 * t34 + t117 * t35 + (t131 * t302 + t322) * t207) + t12 * t85 + t13 * t86 + t7 * t87 + t6 * t88 + t65 * t14 + t22 * t58 + t36 * t33 + t37 * t31 + t38 * t30 - t63 * t310 / 0.2e1 - t105 * t285 / 0.2e1 - t193 * t286 + t266 * t416 + t251 * t275 + t241 * t276; m(3) * qJDD(2) - t338 * t95 + t339 * t94 - t402 * t145 + t403 * t144 + (t150 + t390) * t229 + (t151 + t387) * t232 + (-m(3) - m(4) - t296) * g(3) + ((t193 + t391) * t232 + (t58 + t59 + t396) * t229) * qJD(3) + m(7) * (-t1 * t145 + t144 * t2 + t16 * t95 + t17 * t94 - t232 * t5 + t27 * t303) + m(6) * (t101 * t303 - t144 * t4 - t145 * t3 - t18 * t95 + t19 * t94 - t232 * t61) + m(4) * (t229 * t96 + t232 * t97 + (-t148 * t229 + t149 * t232) * qJD(3)) + m(5) * (-t232 * t89 + t259 * t229 + (t131 * t229 + t232 * t258) * qJD(3)); (-t273 + (-m(7) * t261 - t429) * t232 + t431) * g(3) + (t1 * t120 + t107 * t5 - t255 * t2 + t400 * t27 + (t83 - t25) * t17 + (t84 - t26) * t16) * m(7) - t403 * t255 + (-t101 * t121 + t255 * t4 + t120 * t3 - t208 * t61 + (t83 - t29) * t19 + (-t84 - t28) * t18) * m(6) + (-Ifges(6,2) * t362 + Ifges(7,3) * t363 + t406 * t351 + t408 * t361 - t421) * t136 + (-t393 / 0.2e1 + t382 + (-t251 / 0.2e1 + t241 / 0.2e1) * qJD(1)) * qJD(1) + t223 * t365 + (-pkin(3) * t89 + qJ(4) * t259 + qJD(4) * t258 - t100 * t71 - t70 * t99) * m(5) + (Ifges(5,1) * t223 + t329) * t357 + (Ifges(5,2) * t225 + t330) * t358 + (-mrSges(6,1) * t395 - mrSges(6,2) * t394) * t101 + (-m(5) * t131 + t291 - t396) * t149 - t338 * t84 + (t424 - t425) * t156 + t426 * t157 - (Ifges(6,4) * t362 + Ifges(7,5) * t363 + t351 * t407 + t361 * t409 + t414 + t424) * t137 + t259 * mrSges(5,3) + t390 * qJ(4) + t391 * qJD(4) - t427 * t172 + (-Ifges(4,2) * t282 + Ifges(6,6) * t362 + Ifges(7,6) * t363 - t351 * t405 + t361 * t407 + t428) * t305 + t339 * t83 + t400 * t58 + t402 * t120 + (t290 - t193) * t148 + (Ifges(5,5) * t223 + t225 * Ifges(5,6)) * t352 + (t212 + t384) * t282 + t374 * t171 + t422 * (t272 + (t340 * t417 - t284 + t430) * t232 + (-m(7) * (-t208 - t261) + t247 + m(6) * t208 + t429) * t229) + t89 * t271 + (-t171 * t3 + t18 * t394 + t19 * t395) * mrSges(6,3) + (-t1 * t171 - t16 * t394 + t17 * t395) * mrSges(7,2) + t225 * t63 / 0.2e1 - t208 * t15 - Ifges(4,6) * t179 + Ifges(4,5) * t180 - t100 * t134 - t99 * t135 - t121 * t59 + Ifges(4,3) * qJDD(3) + t107 * t14 - t96 * mrSges(4,2) + t97 * mrSges(4,1) - pkin(3) * t78 - t29 * t85 - t28 * t86 - t26 * t87 - t25 * t88 + t105 * t287 / 0.2e1 + t264 * t276; t338 * t240 + t339 * t113 + t168 * t135 + t250 * t134 + (t232 * g(3) - t229 * t422) * t296 + (t113 * t17 - t16 * t240 + t5) * m(7) + (t113 * t19 + t18 * t240 + t61) * m(6) + (t70 * t168 + t250 * t71 + t89) * m(5) - t387; t52 * t360 + (-m(7) * t16 + t333 + t338) * t19 - t378 + t16 * t336 + t17 * t335 + (-t113 * t409 + t108 - t328 + t49) * t361 + (-m(7) * t17 - t334 - t339) * t18 + (Ifges(7,3) * t240 - t327) * t363 + (-Ifges(6,2) * t240 - t109 + t401) * t362 + (-t113 * t407 + t240 * t406) * t351 - t101 * (mrSges(6,1) * t240 - mrSges(6,2) * t113) - t27 * (mrSges(7,1) * t240 + mrSges(7,3) * t113) + (t128 * t411 + t129 * t410) * g(1) + (t126 * t411 + t127 * t410) * g(2) + (t213 * t411 + t215 * t410) * t343 + t392 + (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t17 - t27 * t57 - (-pkin(5) * t213 + qJ(6) * t215) * t343 - g(2) * (-pkin(5) * t126 + qJ(6) * t127) - g(1) * (-pkin(5) * t128 + qJ(6) * t129)) * m(7) + qJD(6) * t88 - t57 * t58 - pkin(5) * t31 + qJ(6) * t33; t240 * t58 - t202 * t88 + (-g(1) * t128 - g(2) * t126 - t17 * t202 - t213 * t343 + t240 * t27 + t2) * m(7) + t31;];
tau  = t8;
