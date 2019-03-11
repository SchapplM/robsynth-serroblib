% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:42
% EndTime: 2019-03-09 02:13:00
% DurationCPUTime: 13.49s
% Computational Cost: add. (6147->521), mult. (12827->639), div. (0->0), fcn. (8713->10), ass. (0->238)
t173 = -qJ(6) - pkin(8);
t389 = -m(6) * pkin(8) + m(7) * t173 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t167 = pkin(9) + qJ(4);
t159 = sin(t167);
t160 = cos(t167);
t388 = mrSges(5,1) * t159 + t389 * t160;
t371 = Ifges(6,4) + Ifges(7,4);
t177 = sin(qJ(4));
t172 = cos(pkin(9));
t313 = cos(qJ(4));
t241 = t313 * t172;
t171 = sin(pkin(9));
t260 = qJD(1) * t171;
t125 = qJD(1) * t241 - t177 * t260;
t176 = sin(qJ(5));
t179 = cos(qJ(5));
t104 = qJD(4) * t179 - t125 * t176;
t190 = -t171 * t313 - t177 * t172;
t124 = t190 * qJD(1);
t189 = -t177 * t171 + t241;
t94 = qJD(4) * t124 + qJDD(1) * t189;
t54 = qJD(5) * t104 + qJDD(4) * t176 + t179 * t94;
t329 = t54 / 0.2e1;
t105 = qJD(4) * t176 + t125 * t179;
t55 = -qJD(5) * t105 + qJDD(4) * t179 - t176 * t94;
t328 = t55 / 0.2e1;
t95 = -t189 * qJD(1) * qJD(4) + qJDD(1) * t190;
t91 = qJDD(5) - t95;
t327 = t91 / 0.2e1;
t372 = Ifges(6,1) + Ifges(7,1);
t370 = Ifges(6,5) + Ifges(7,5);
t369 = Ifges(6,2) + Ifges(7,2);
t368 = Ifges(6,6) + Ifges(7,6);
t367 = Ifges(7,3) + Ifges(6,3);
t261 = t171 ^ 2 + t172 ^ 2;
t220 = t261 * mrSges(4,3);
t178 = sin(qJ(1));
t180 = cos(qJ(1));
t345 = -g(1) * t178 + g(2) * t180;
t175 = -pkin(1) - qJ(3);
t387 = -qJD(1) * qJD(3) + qJDD(1) * t175;
t156 = pkin(5) * t179 + pkin(4);
t386 = -m(6) * pkin(4) - m(7) * t156;
t384 = t370 * t327 + t371 * t328 + t372 * t329;
t383 = t371 * t104;
t257 = qJD(5) * t179;
t258 = qJD(5) * t176;
t133 = qJDD(2) + t387;
t215 = -pkin(7) * qJDD(1) + t133;
t109 = t215 * t171;
t110 = t215 * t172;
t140 = qJD(1) * t175 + qJD(2);
t219 = -pkin(7) * qJD(1) + t140;
t114 = t219 * t171;
t115 = t219 * t172;
t229 = qJD(4) * t313;
t259 = qJD(4) * t177;
t31 = t313 * t109 + t177 * t110 - t114 * t259 + t115 * t229;
t27 = qJDD(4) * pkin(8) + t31;
t169 = qJD(1) * qJD(2);
t141 = -qJDD(1) * qJ(2) - t169;
t139 = qJDD(3) - t141;
t254 = qJDD(1) * t171;
t128 = pkin(3) * t254 + t139;
t39 = -pkin(4) * t95 - pkin(8) * t94 + t128;
t74 = t114 * t313 + t177 * t115;
t67 = qJD(4) * pkin(8) + t74;
t157 = qJD(1) * qJ(2) + qJD(3);
t136 = pkin(3) * t260 + t157;
t68 = -pkin(4) * t124 - pkin(8) * t125 + t136;
t3 = t176 * t39 + t179 * t27 + t68 * t257 - t258 * t67;
t30 = t176 * t68 + t179 * t67;
t4 = -qJD(5) * t30 - t176 * t27 + t179 * t39;
t212 = -t176 * t4 + t179 * t3;
t29 = -t176 * t67 + t179 * t68;
t382 = -t29 * t257 - t30 * t258 + t212;
t381 = t371 * t105;
t380 = t371 * t179;
t379 = t371 * t176;
t378 = qJD(5) - t124;
t330 = m(7) * pkin(5);
t376 = t368 * t91 + t369 * t55 + t371 * t54;
t374 = -mrSges(7,1) - mrSges(6,1);
t373 = mrSges(6,2) + mrSges(7,2);
t16 = -mrSges(6,1) * t55 + mrSges(6,2) * t54;
t366 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t94 + t16;
t365 = t104 * t368 + t105 * t370 + t367 * t378;
t364 = t104 * t369 + t368 * t378 + t381;
t363 = t105 * t372 + t370 * t378 + t383;
t116 = Ifges(5,4) * t124;
t362 = t124 * Ifges(5,2);
t296 = mrSges(5,3) * t125;
t361 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t104 - mrSges(6,2) * t105 - t296;
t221 = qJD(5) * t173;
t277 = t124 * t176;
t73 = -t177 * t114 + t115 * t313;
t93 = pkin(4) * t125 - pkin(8) * t124;
t38 = t176 * t93 + t179 * t73;
t360 = qJ(6) * t277 + qJD(6) * t179 + t176 * t221 - t38;
t276 = t124 * t179;
t37 = -t176 * t73 + t179 * t93;
t359 = -pkin(5) * t125 + qJ(6) * t276 - qJD(6) * t176 + t179 * t221 - t37;
t358 = -t74 + (t258 - t277) * pkin(5);
t357 = t330 + mrSges(7,1);
t356 = Ifges(5,5) * qJD(4);
t355 = Ifges(5,6) * qJD(4);
t303 = -pkin(7) + t175;
t134 = t303 * t171;
t135 = t303 * t172;
t354 = -t177 * t134 + t313 * t135;
t206 = -mrSges(7,1) * t179 + mrSges(7,2) * t176;
t208 = -mrSges(6,1) * t179 + mrSges(6,2) * t176;
t353 = -t206 - t208 - t386;
t298 = mrSges(7,2) * t179;
t205 = mrSges(7,1) * t176 + t298;
t207 = mrSges(6,1) * t176 + mrSges(6,2) * t179;
t66 = -qJD(4) * pkin(4) - t73;
t40 = -t104 * pkin(5) + qJD(6) + t66;
t352 = t40 * t205 + t66 * t207;
t351 = -t176 * t368 + t179 * t370;
t350 = -t176 * t369 + t380;
t349 = t179 * t372 - t379;
t346 = t367 * t91 + t368 * t55 + t370 * t54;
t343 = -m(7) - m(6) - m(5);
t341 = -mrSges(6,1) - t357;
t126 = -t171 * t229 - t172 * t259;
t127 = -t171 * t259 + t172 * t229;
t32 = -t177 * t109 + t110 * t313 - t114 * t229 - t115 * t259;
t339 = -t126 * t73 - t127 * t74 - t189 * t32;
t338 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t1 = pkin(5) * t91 - qJ(6) * t54 - qJD(6) * t105 + t4;
t17 = -qJ(6) * t105 + t29;
t14 = pkin(5) * t378 + t17;
t18 = qJ(6) * t104 + t30;
t2 = qJ(6) * t55 + qJD(6) * t104 + t3;
t337 = -t1 * t176 - t14 * t257 + t179 * t2 - t18 * t258;
t336 = -m(6) * t66 + t361;
t335 = t136 * mrSges(5,2) + t356 / 0.2e1;
t211 = mrSges(4,1) * t171 + mrSges(4,2) * t172;
t334 = t159 * t386 + mrSges(2,2) - mrSges(3,3) - t211 - t388;
t333 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t332 = m(4) * t157 + m(5) * t136 - mrSges(5,1) * t124 + mrSges(5,2) * t125 + t211 * qJD(1);
t331 = t18 * mrSges(7,2) + t30 * mrSges(6,2) + t355 / 0.2e1 - t136 * mrSges(5,1) - t14 * mrSges(7,1) - t29 * mrSges(6,1);
t326 = -t104 / 0.2e1;
t325 = t104 / 0.2e1;
t324 = -t105 / 0.2e1;
t323 = t105 / 0.2e1;
t322 = -t378 / 0.2e1;
t321 = t378 / 0.2e1;
t320 = -t124 / 0.2e1;
t319 = -t125 / 0.2e1;
t318 = t125 / 0.2e1;
t161 = t171 * pkin(3);
t20 = mrSges(7,1) * t91 - mrSges(7,3) * t54;
t21 = mrSges(6,1) * t91 - mrSges(6,3) * t54;
t302 = -t20 - t21;
t22 = -mrSges(7,2) * t91 + mrSges(7,3) * t55;
t23 = -mrSges(6,2) * t91 + mrSges(6,3) * t55;
t301 = t22 + t23;
t293 = mrSges(7,3) * t104;
t69 = -mrSges(7,2) * t378 + t293;
t295 = mrSges(6,3) * t104;
t70 = -mrSges(6,2) * t378 + t295;
t300 = t69 + t70;
t292 = mrSges(7,3) * t105;
t71 = mrSges(7,1) * t378 - t292;
t294 = mrSges(6,3) * t105;
t72 = mrSges(6,1) * t378 - t294;
t299 = t71 + t72;
t149 = qJ(2) + t161;
t92 = -pkin(4) * t190 - pkin(8) * t189 + t149;
t98 = t134 * t313 + t177 * t135;
t96 = t179 * t98;
t42 = t176 * t92 + t96;
t297 = mrSges(5,3) * t124;
t291 = Ifges(5,4) * t125;
t28 = -qJDD(4) * pkin(4) - t32;
t11 = -t55 * pkin(5) + qJDD(6) + t28;
t284 = t11 * t189;
t281 = t189 * t28;
t278 = qJDD(1) * pkin(1);
t275 = t127 * t176;
t274 = t127 * t179;
t273 = t189 * t176;
t272 = t189 * t179;
t271 = t176 * t126;
t270 = t176 * t178;
t269 = t176 * t180;
t266 = t178 * t179;
t265 = t179 * t126;
t264 = t179 * t180;
t253 = qJDD(1) * t172;
t263 = mrSges(4,1) * t254 + mrSges(4,2) * t253;
t262 = t180 * pkin(1) + t178 * qJ(2);
t62 = t190 * qJD(3) + qJD(4) * t354;
t88 = pkin(4) * t127 - pkin(8) * t126 + qJD(2);
t250 = t176 * t88 + t179 * t62 + t92 * t257;
t57 = -mrSges(7,1) * t104 + mrSges(7,2) * t105;
t249 = -t57 + t361;
t247 = -m(4) + t343;
t239 = t189 * t257;
t230 = -t95 * mrSges(5,1) + t94 * mrSges(5,2);
t15 = -t55 * mrSges(7,1) + t54 * mrSges(7,2);
t224 = t257 / 0.2e1;
t163 = t180 * qJ(2);
t223 = -pkin(1) * t178 + t163;
t222 = -t176 * t62 + t179 * t88;
t41 = -t176 * t98 + t179 * t92;
t218 = t261 * t140;
t217 = t261 * t133;
t197 = t30 * t176 + t29 * t179;
t195 = -qJ(6) * t126 - qJD(6) * t189;
t120 = t159 * t269 + t266;
t118 = -t159 * t270 + t264;
t188 = -t239 - t271;
t187 = -t189 * t258 + t265;
t184 = -t176 * t300 - t179 * t299;
t63 = qJD(3) * t189 + qJD(4) * t98;
t181 = qJD(1) ^ 2;
t174 = -pkin(7) - qJ(3);
t158 = qJDD(2) - t278;
t138 = t173 * t179;
t137 = t173 * t176;
t121 = t159 * t264 - t270;
t119 = t159 * t266 + t269;
t112 = -qJD(4) * mrSges(5,2) + t297;
t81 = t125 * Ifges(5,1) + t116 + t356;
t80 = t291 + t355 + t362;
t79 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t95;
t64 = pkin(5) * t273 - t354;
t36 = -pkin(5) * t188 + t63;
t33 = -qJ(6) * t273 + t42;
t25 = -pkin(5) * t190 - qJ(6) * t272 + t41;
t13 = -qJD(5) * t42 + t222;
t12 = -t258 * t98 + t250;
t6 = -qJ(6) * t239 + (-qJD(5) * t98 + t195) * t176 + t250;
t5 = pkin(5) * t127 + t195 * t179 + (-t96 + (qJ(6) * t189 - t92) * t176) * qJD(5) + t222;
t7 = [m(4) * (qJ(2) * t139 - qJD(3) * t218 + t175 * t217) + (-t1 * t272 - t14 * t187 + t18 * t188 - t2 * t273) * mrSges(7,3) + (-t187 * t29 + t188 * t30 - t272 * t4 - t273 * t3) * mrSges(6,3) + (t187 * t370 + t188 * t368) * t321 + m(6) * (t12 * t30 + t13 * t29 + t3 * t42 + t4 * t41) + m(5) * (t128 * t149 + t31 * t98 + t62 * t74) + (-m(5) * t73 - t336) * t63 - (t363 * qJD(5) + t376) * t273 / 0.2e1 + (t128 * mrSges(5,2) + Ifges(5,1) * t94 + Ifges(5,4) * t95 + Ifges(5,5) * qJDD(4) - t224 * t364 + t327 * t351 + t328 * t350 + t329 * t349) * t189 + (-t128 * mrSges(5,1) - t1 * mrSges(7,1) + t31 * mrSges(5,3) + Ifges(5,4) * t94 + Ifges(5,2) * t95 + Ifges(5,6) * qJDD(4) - t327 * t367 - t328 * t368 - t329 * t370 - t333 - t346 / 0.2e1) * t190 + t339 * mrSges(5,3) + (-t269 * t330 + (-m(4) - m(3)) * t262 + t343 * (t178 * t161 - t174 * t180 + t262) + t374 * t119 - t373 * t118 + (-m(4) * qJ(3) + t338) * t180 + t334 * t178) * g(2) + (t270 * t330 - m(3) * t223 - m(4) * t163 + t343 * (t180 * t161 + t178 * t174 + t223) + t374 * t121 + t373 * t120 + (-m(4) * t175 - t338) * t178 + t334 * t180) * g(1) + (Ifges(5,1) * t318 + t116 / 0.2e1 + t81 / 0.2e1 + t335) * t126 + t363 * t265 / 0.2e1 - t364 * t271 / 0.2e1 - (-m(5) * t32 + m(6) * t28 + t366) * t354 + t332 * qJD(2) + t40 * (-mrSges(7,1) * t188 + mrSges(7,2) * t187) + t66 * (-mrSges(6,1) * t188 + mrSges(6,2) * t187) + t149 * t230 + t205 * t284 + t207 * t281 + (t187 * t371 + t188 * t369) * t325 + (t187 * t372 + t188 * t371) * t323 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + (Ifges(4,1) * t172 - Ifges(4,4) * t171) * t253 + (-t133 - t387) * t220 - (Ifges(4,4) * t172 - Ifges(4,2) * t171) * t254 + m(3) * (-pkin(1) * t158 + (-t141 + t169) * qJ(2)) + (-t278 + t158) * mrSges(3,2) + (t365 / 0.2e1 - Ifges(5,4) * t318 - t362 / 0.2e1 - t80 / 0.2e1 + t368 * t325 + t370 * t323 + t367 * t321 - t331) * t127 + qJ(2) * t263 + t139 * t211 - 0.2e1 * t141 * mrSges(3,3) + t62 * t112 + m(7) * (t1 * t25 + t11 * t64 + t14 * t5 + t18 * t6 + t2 * t33 + t36 * t40) + t272 * t384 + t25 * t20 + t33 * t22 + t41 * t21 + t42 * t23 + t36 * t57 + t64 * t15 + t6 * t69 + t12 * t70 + t5 * t71 + t13 * t72 + t98 * t79; (-m(3) * qJ(2) - mrSges(3,3)) * t181 - (t15 + t366) * t189 + t249 * t126 + (mrSges(3,2) - t220) * qJDD(1) + (-t176 * t299 + t179 * t300 + t112) * t127 - m(5) * t339 + m(7) * (-t126 * t40 - t14 * t275 + t18 * t274 - t284) + m(6) * (-t126 * t66 + t274 * t30 - t275 * t29 - t281) + m(3) * t158 + m(4) * t217 - (m(5) * t31 + m(6) * t382 + m(7) * t337 + t184 * qJD(5) + t302 * t176 + t301 * t179 + t79) * t190 + (-m(7) * (t14 * t179 + t18 * t176) - m(6) * t197 + t184 - t332) * qJD(1) + t345 * (m(3) - t247); -t124 * t112 - t181 * t220 + t249 * t125 + (t300 * t378 - t302) * t179 + (-t299 * t378 + t301) * t176 + t230 + t263 + (g(1) * t180 + g(2) * t178) * t247 + (t1 * t179 - t125 * t40 + t2 * t176 + t378 * (-t14 * t176 + t18 * t179)) * m(7) + (-t125 * t66 + t3 * t176 + t4 * t179 + t378 * (-t29 * t176 + t30 * t179)) * m(6) + (-t124 * t74 + t125 * t73 + t128) * m(5) + (qJD(1) * t218 + t139) * m(4); (-t291 + t365) * t319 + (-Ifges(5,2) * t320 + t322 * t367 + t324 * t370 + t326 * t368 + t331) * t125 + (t176 * t370 + t179 * t368) * t327 + (t297 - t112) * t73 + (-t258 / 0.2e1 + t277 / 0.2e1) * t364 + (t296 + t336) * t74 + (t14 * t276 + t18 * t277 + t337) * mrSges(7,3) + (m(6) * (-qJD(5) * t197 + t212) - t72 * t257 - t70 * t258 - t176 * t21 + t179 * t23) * pkin(8) + t376 * t179 / 0.2e1 + (t104 * t350 + t105 * t349 + t351 * t378) * qJD(5) / 0.2e1 + (t159 * t353 + t388) * g(3) + (t179 * t369 + t379) * t328 + (t176 * t372 + t380) * t329 + (t276 * t29 + t277 * t30 + t382) * mrSges(6,3) - t345 * ((-mrSges(5,1) - t353) * t160 + t389 * t159) + t80 * t318 + (Ifges(5,1) * t319 + t322 * t351 + t324 * t349 + t326 * t350 - t335 - t352) * t124 + t352 * qJD(5) + t358 * t57 + t359 * t71 + t360 * t69 + (t1 * t137 - t11 * t156 - t138 * t2 + t14 * t359 + t18 * t360 + t358 * t40) * m(7) + t11 * t206 + t28 * t208 + (-pkin(4) * t28 - t29 * t37 - t30 * t38) * m(6) + (t116 + t81) * t320 - t156 * t15 + Ifges(5,3) * qJDD(4) + t137 * t20 - t138 * t22 + t176 * t384 - pkin(4) * t16 + (-t276 / 0.2e1 + t224) * t363 - t31 * mrSges(5,2) + t32 * mrSges(5,1) - t38 * t70 - t37 * t72 + Ifges(5,5) * t94 + Ifges(5,6) * t95; t364 * t323 + (t294 + t72) * t30 + (t295 - t70) * t29 + (t292 - m(7) * (-t14 + t17) + t71) * t18 + t357 * t1 + (t120 * t341 - t121 * t373) * g(2) + (t104 * t370 - t105 * t368) * t322 + t14 * t293 + (-t105 * t369 + t363 + t383) * t326 + (t176 * t357 + t207 + t298) * g(3) * t160 + (t118 * t341 + t119 * t373) * g(1) + t333 + (t104 * t372 - t381) * t324 - t40 * (mrSges(7,1) * t105 + mrSges(7,2) * t104) - t66 * (mrSges(6,1) * t105 + mrSges(6,2) * t104) - t17 * t69 + t346 + ((-m(7) * t40 - t57) * t105 + t20) * pkin(5); -t104 * t69 + t105 * t71 + (-g(3) * t159 - t18 * t104 + t14 * t105 - t160 * t345 + t11) * m(7) + t15;];
tau  = t7;
