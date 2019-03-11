% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP8
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:17
% EndTime: 2019-03-09 02:15:33
% DurationCPUTime: 12.30s
% Computational Cost: add. (6204->563), mult. (12816->709), div. (0->0), fcn. (8665->10), ass. (0->246)
t367 = Ifges(6,1) + Ifges(7,1);
t379 = Ifges(6,4) - Ifges(7,5);
t366 = Ifges(7,4) + Ifges(6,5);
t365 = Ifges(6,6) - Ifges(7,6);
t364 = Ifges(6,3) + Ifges(7,2);
t171 = sin(pkin(9));
t172 = cos(pkin(9));
t258 = t171 ^ 2 + t172 ^ 2;
t225 = t258 * mrSges(4,3);
t174 = -pkin(1) - qJ(3);
t378 = -qJD(1) * qJD(3) + qJDD(1) * t174;
t176 = sin(qJ(4));
t316 = cos(qJ(4));
t189 = -t171 * t316 - t176 * t172;
t118 = t189 * qJD(1);
t373 = qJD(5) - t118;
t167 = pkin(9) + qJ(4);
t158 = cos(t167);
t309 = g(3) * t158;
t238 = t316 * t172;
t257 = qJD(1) * t171;
t119 = qJD(1) * t238 - t176 * t257;
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t100 = qJD(4) * t175 + t119 * t178;
t188 = -t176 * t171 + t238;
t89 = qJD(4) * t118 + qJDD(1) * t188;
t54 = qJD(5) * t100 - t178 * qJDD(4) + t175 * t89;
t331 = -t54 / 0.2e1;
t377 = Ifges(6,2) * t331;
t211 = t178 * mrSges(7,1) + t175 * mrSges(7,3);
t213 = mrSges(6,1) * t178 - mrSges(6,2) * t175;
t376 = -t211 - t213;
t254 = qJD(5) * t178;
t255 = qJD(5) * t175;
t127 = qJDD(2) + t378;
t220 = -pkin(7) * qJDD(1) + t127;
t104 = t220 * t171;
t105 = t220 * t172;
t137 = qJD(1) * t174 + qJD(2);
t224 = -pkin(7) * qJD(1) + t137;
t109 = t224 * t171;
t110 = t224 * t172;
t231 = qJD(4) * t316;
t256 = qJD(4) * t176;
t31 = t316 * t104 + t176 * t105 - t109 * t256 + t110 * t231;
t26 = qJDD(4) * pkin(8) + t31;
t169 = qJD(1) * qJD(2);
t138 = -qJDD(1) * qJ(2) - t169;
t132 = qJDD(3) - t138;
t251 = qJDD(1) * t171;
t122 = pkin(3) * t251 + t132;
t90 = -t188 * qJD(1) * qJD(4) + qJDD(1) * t189;
t38 = -pkin(4) * t90 - pkin(8) * t89 + t122;
t72 = t109 * t316 + t176 * t110;
t65 = qJD(4) * pkin(8) + t72;
t155 = qJD(1) * qJ(2) + qJD(3);
t130 = pkin(3) * t257 + t155;
t66 = -pkin(4) * t118 - pkin(8) * t119 + t130;
t3 = t175 * t38 + t178 * t26 + t66 * t254 - t255 * t65;
t29 = t175 * t66 + t178 * t65;
t4 = -qJD(5) * t29 - t175 * t26 + t178 * t38;
t217 = -t175 * t4 + t178 * t3;
t28 = -t175 * t65 + t178 * t66;
t375 = -t28 * t254 - t29 * t255 + t217;
t17 = -pkin(5) * t373 + qJD(6) - t28;
t18 = qJ(6) * t373 + t29;
t86 = qJDD(5) - t90;
t1 = qJ(6) * t86 + qJD(6) * t373 + t3;
t2 = -pkin(5) * t86 + qJDD(6) - t4;
t218 = t1 * t178 + t175 * t2;
t374 = t17 * t254 - t18 * t255 + t218;
t194 = t178 * qJD(4) - t119 * t175;
t53 = qJD(5) * t194 + qJDD(4) * t175 + t178 * t89;
t20 = mrSges(6,1) * t86 - mrSges(6,3) * t53;
t21 = -t86 * mrSges(7,1) + t53 * mrSges(7,2);
t299 = -t20 + t21;
t19 = -mrSges(7,2) * t54 + mrSges(7,3) * t86;
t22 = -mrSges(6,2) * t86 - mrSges(6,3) * t54;
t300 = t19 + t22;
t372 = t299 * t175 + t300 * t178;
t332 = t53 / 0.2e1;
t329 = t86 / 0.2e1;
t371 = m(6) + m(7);
t370 = mrSges(6,1) + mrSges(7,1);
t369 = mrSges(6,2) - mrSges(7,3);
t368 = mrSges(6,3) + mrSges(7,2);
t363 = t366 * t86 + t367 * t53 - t379 * t54;
t16 = mrSges(6,1) * t54 + mrSges(6,2) * t53;
t362 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t89 + t16;
t361 = t100 * t366 + t194 * t365 + t364 * t373;
t303 = t194 * Ifges(7,5);
t95 = Ifges(6,4) * t194;
t360 = t100 * t367 + t366 * t373 - t303 + t95;
t315 = mrSges(7,2) * t194;
t67 = mrSges(7,3) * t373 + t315;
t314 = mrSges(6,3) * t194;
t68 = -mrSges(6,2) * t373 + t314;
t298 = t67 + t68;
t292 = mrSges(6,3) * t100;
t69 = mrSges(6,1) * t373 - t292;
t295 = mrSges(7,2) * t100;
t70 = -mrSges(7,1) * t373 + t295;
t297 = t69 - t70;
t293 = mrSges(5,3) * t119;
t359 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t194 - mrSges(6,2) * t100 - t293;
t202 = pkin(5) * t175 - qJ(6) * t178;
t358 = -qJD(6) * t175 + t202 * t373 - t72;
t357 = Ifges(5,5) * qJD(4);
t356 = Ifges(5,6) * qJD(4);
t157 = sin(t167);
t177 = sin(qJ(1));
t272 = t157 * t177;
t271 = t158 * t177;
t301 = -pkin(7) + t174;
t128 = t301 * t171;
t129 = t301 * t172;
t355 = -t176 * t128 + t316 * t129;
t203 = pkin(5) * t178 + qJ(6) * t175;
t131 = -pkin(4) - t203;
t354 = m(6) * pkin(4) - m(7) * t131 - t376;
t210 = t175 * mrSges(7,1) - t178 * mrSges(7,3);
t212 = mrSges(6,1) * t175 + mrSges(6,2) * t178;
t71 = -t176 * t109 + t110 * t316;
t64 = -qJD(4) * pkin(4) - t71;
t30 = -pkin(5) * t194 - t100 * qJ(6) + t64;
t353 = t30 * t210 + t64 * t212;
t352 = -t175 * t365 + t178 * t366;
t287 = Ifges(7,5) * t175;
t289 = Ifges(6,4) * t175;
t351 = t178 * t367 + t287 - t289;
t348 = t364 * t86 - t365 * t54 + t366 * t53;
t347 = Ifges(6,4) * t332 + t377 + Ifges(6,6) * t329 - t53 * Ifges(7,5) / 0.2e1 - t86 * Ifges(7,6) / 0.2e1 + Ifges(7,3) * t331;
t344 = g(1) * t272;
t120 = -t171 * t231 - t172 * t256;
t121 = -t171 * t256 + t172 * t231;
t32 = -t176 * t104 + t105 * t316 - t109 * t231 - t110 * t256;
t342 = t120 * t71 + t121 * t72 + t188 * t32;
t341 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t214 = t157 * mrSges(5,1) + t158 * mrSges(5,2);
t216 = mrSges(4,1) * t171 + mrSges(4,2) * t172;
t340 = -t214 - t216 - mrSges(3,3) + mrSges(2,2);
t161 = t171 * pkin(3);
t146 = qJ(2) + t161;
t87 = -pkin(4) * t189 - pkin(8) * t188 + t146;
t93 = t128 * t316 + t176 * t129;
t296 = t175 * t87 + t178 * t93;
t61 = t189 * qJD(3) + qJD(4) * t355;
t83 = pkin(4) * t121 - pkin(8) * t120 + qJD(2);
t13 = -qJD(5) * t296 - t175 * t61 + t178 * t83;
t339 = -m(7) * pkin(5) - t370;
t338 = -m(6) * t64 + t359;
t337 = -m(7) * qJ(6) + t369;
t336 = m(4) * t155 + m(5) * t130 - mrSges(5,1) * t118 + mrSges(5,2) * t119 + t216 * qJD(1);
t335 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t330 = t54 / 0.2e1;
t328 = t194 / 0.2e1;
t327 = -t194 / 0.2e1;
t326 = -t100 / 0.2e1;
t325 = t100 / 0.2e1;
t324 = -t373 / 0.2e1;
t322 = -t118 / 0.2e1;
t321 = -t119 / 0.2e1;
t320 = t119 / 0.2e1;
t311 = g(1) * t177;
t179 = cos(qJ(1));
t310 = g(2) * t179;
t27 = -qJDD(4) * pkin(4) - t32;
t5 = t54 * pkin(5) - t53 * qJ(6) - t100 * qJD(6) + t27;
t307 = t188 * t5;
t88 = pkin(4) * t119 - pkin(8) * t118;
t37 = t175 * t88 + t178 * t71;
t294 = mrSges(5,3) * t118;
t291 = Ifges(5,4) * t119;
t290 = Ifges(6,4) * t100;
t288 = Ifges(6,4) * t178;
t286 = Ifges(7,5) * t178;
t283 = t188 * t27;
t279 = qJDD(1) * pkin(1);
t278 = t118 * t175;
t277 = t118 * t178;
t276 = t121 * t175;
t275 = t121 * t178;
t273 = t188 * t178;
t270 = t158 * t179;
t269 = t175 * t120;
t268 = t175 * t177;
t267 = t175 * t179;
t264 = t177 * t178;
t263 = t178 * t120;
t262 = t179 * t178;
t261 = pkin(4) * t271 + pkin(8) * t272;
t250 = qJDD(1) * t172;
t260 = mrSges(4,1) * t251 + mrSges(4,2) * t250;
t259 = t179 * pkin(1) + t177 * qJ(2);
t57 = -mrSges(7,1) * t194 - mrSges(7,3) * t100;
t247 = -t57 + t359;
t244 = -m(4) - m(5) - t371;
t45 = Ifges(6,2) * t194 + Ifges(6,6) * t373 + t290;
t239 = -t175 * t45 / 0.2e1;
t232 = -t90 * mrSges(5,1) + t89 * mrSges(5,2);
t228 = t255 / 0.2e1;
t227 = t254 / 0.2e1;
t163 = t179 * qJ(2);
t226 = -pkin(1) * t177 + t163;
t223 = t258 * t137;
t222 = t258 * t127;
t215 = mrSges(5,1) * t158 - mrSges(5,2) * t157;
t207 = -Ifges(6,2) * t175 + t288;
t204 = Ifges(7,3) * t175 + t286;
t201 = t17 * t178 - t18 * t175;
t198 = t29 * t175 + t28 * t178;
t36 = -t175 * t71 + t178 * t88;
t40 = -t175 * t93 + t178 * t87;
t173 = -pkin(7) - qJ(3);
t193 = t179 * t161 + t177 * t173 + t226;
t192 = t177 * t161 - t173 * t179 + t259;
t12 = t175 * t83 + t178 * t61 + t87 * t254 - t255 * t93;
t187 = -t188 * t254 - t269;
t186 = -t188 * t255 + t263;
t184 = -t175 * t298 - t178 * t297;
t62 = qJD(3) * t188 + qJD(4) * t93;
t180 = qJD(1) ^ 2;
t156 = qJDD(2) - t279;
t116 = t157 * t262 - t268;
t115 = t157 * t267 + t264;
t114 = t157 * t264 + t267;
t113 = t157 * t268 - t262;
t111 = Ifges(5,4) * t118;
t107 = -qJD(4) * mrSges(5,2) + t294;
t94 = Ifges(7,5) * t100;
t78 = t119 * Ifges(5,1) + t111 + t357;
t77 = t118 * Ifges(5,2) + t291 + t356;
t76 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t90;
t56 = pkin(5) * t100 - qJ(6) * t194;
t55 = t188 * t202 - t355;
t42 = Ifges(7,6) * t373 - Ifges(7,3) * t194 + t94;
t35 = pkin(5) * t189 - t40;
t34 = -qJ(6) * t189 + t296;
t25 = -pkin(5) * t119 - t36;
t24 = qJ(6) * t119 + t37;
t15 = mrSges(7,1) * t54 - mrSges(7,3) * t53;
t14 = t202 * t120 - (-qJD(5) * t203 + qJD(6) * t178) * t188 + t62;
t7 = -pkin(5) * t121 - t13;
t6 = qJ(6) * t121 - qJD(6) * t189 + t12;
t8 = [(-t279 + t156) * mrSges(3,2) + m(6) * (t12 * t29 + t13 * t28 + t296 * t3 + t4 * t40) + t296 * t22 + t42 * t269 / 0.2e1 + t132 * t216 + t146 * t232 + t18 * (mrSges(7,2) * t187 + mrSges(7,3) * t121) + t29 * (-mrSges(6,2) * t121 + mrSges(6,3) * t187) + t30 * (-mrSges(7,1) * t187 - mrSges(7,3) * t186) + t64 * (-mrSges(6,1) * t187 + mrSges(6,2) * t186) + t17 * (-mrSges(7,1) * t121 + mrSges(7,2) * t186) + t28 * (mrSges(6,1) * t121 - mrSges(6,3) * t186) - (-m(5) * t32 + m(6) * t27 + t362) * t355 + (t121 * t364 + t186 * t366 + t187 * t365) * t373 / 0.2e1 + (-m(5) * t71 - t338) * t62 + m(7) * (t1 * t34 + t14 * t30 + t17 * t7 + t18 * t6 + t2 * t35 + t5 * t55) + (Ifges(4,1) * t172 - Ifges(4,4) * t171) * t250 + m(5) * (t122 * t146 + t31 * t93 + t61 * t72) + t360 * t263 / 0.2e1 + t361 * t121 / 0.2e1 + (qJD(5) * t42 + t363) * t273 / 0.2e1 + (t366 * t121 + t367 * t186 + t187 * t379) * t325 + qJ(2) * t260 + m(3) * (-pkin(1) * t156 + (-t138 + t169) * qJ(2)) + (Ifges(7,5) * t186 + Ifges(7,6) * t121 - Ifges(7,3) * t187) * t327 + (Ifges(6,4) * t186 + Ifges(6,2) * t187 + Ifges(6,6) * t121) * t328 - (Ifges(4,4) * t172 - Ifges(4,2) * t171) * t251 + (Ifges(5,1) * t120 - Ifges(5,4) * t121) * t320 + t336 * qJD(2) + t212 * t283 + m(4) * (qJ(2) * t132 - qJD(3) * t223 + t174 * t222) - t342 * mrSges(5,3) + (-(mrSges(7,2) * t1 + mrSges(6,3) * t3 + t347) * t175 + t122 * mrSges(5,2) + Ifges(5,1) * t89 + Ifges(5,4) * t90 + Ifges(5,5) * qJDD(4) + t204 * t330 + t207 * t331 - t227 * t45 + t329 * t352 + t332 * t351 - t360 * t228) * t188 + (-t348 / 0.2e1 - t122 * mrSges(5,1) + t31 * mrSges(5,3) + Ifges(5,4) * t89 + Ifges(5,2) * t90 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t331 - Ifges(7,6) * t330 - t329 * t364 - t332 * t366 - t335) * t189 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + (-m(3) * t226 - m(4) * t163 - m(5) * t193 + t368 * t270 - t371 * (t179 * t157 * pkin(4) - pkin(8) * t270 + t193) + t339 * t116 + t337 * t115 + t340 * t179 + (-m(4) * t174 - t341) * t177) * g(1) + (-m(5) * t192 + t368 * t271 + (-m(4) - m(3)) * t259 - t371 * (pkin(4) * t272 - pkin(8) * t271 + t192) + t339 * t114 + t337 * t113 + (-m(4) * qJ(3) + t341) * t179 + t340 * t177) * g(2) - (-mrSges(7,2) * t2 + mrSges(6,3) * t4) * t273 - 0.2e1 * t138 * mrSges(3,3) + t130 * (mrSges(5,1) * t121 + mrSges(5,2) * t120) + t210 * t307 + t120 * t239 + (-t127 - t378) * t225 + t34 * t19 + t35 * t21 + t40 * t20 + t55 * t15 + t14 * t57 + t6 * t67 + t12 * t68 + t13 * t69 + t7 * t70 + t93 * t76 + t61 * t107 + t120 * t78 / 0.2e1 - t121 * t77 / 0.2e1 + t118 * (Ifges(5,4) * t120 - Ifges(5,2) * t121) / 0.2e1 + qJD(4) * (Ifges(5,5) * t120 - Ifges(5,6) * t121) / 0.2e1; (-m(3) * qJ(2) - mrSges(3,3)) * t180 - (t15 + t362) * t188 + t247 * t120 + (mrSges(3,2) - t225) * qJDD(1) + (-t175 * t297 + t178 * t298 + t107) * t121 + m(5) * t342 + m(7) * (-t120 * t30 + t17 * t276 + t18 * t275 - t307) + m(6) * (-t120 * t64 + t275 * t29 - t276 * t28 - t283) + m(3) * t156 + m(4) * t222 - (m(5) * t31 + m(6) * t375 + m(7) * t374 + t184 * qJD(5) + t372 + t76) * t189 + (-m(6) * t198 + m(7) * t201 + t184 - t336) * qJD(1) + (t310 - t311) * (m(3) - t244); -t118 * t107 - t180 * t225 + t247 * t119 + (t298 * t373 - t299) * t178 + (-t297 * t373 + t300) * t175 + t232 + t260 + (g(1) * t179 + g(2) * t177) * t244 + (t1 * t175 - t119 * t30 - t2 * t178 + t373 * (t17 * t175 + t18 * t178)) * m(7) + (-t119 * t64 + t3 * t175 + t4 * t178 + t373 * (-t28 * t175 + t29 * t178)) * m(6) + (-t118 * t72 + t119 * t71 + t122) * m(5) + (qJD(1) * t223 + t132) * m(4); (t294 - t107) * t71 - t27 * t213 - t5 * t211 + (t45 / 0.2e1 - t42 / 0.2e1) * t278 + ((qJD(5) * t201 + t218) * m(7) - t297 * t254 - t298 * t255 + (-qJD(5) * t198 + t217) * m(6) + (t157 * t310 - t309) * t371 + t372) * pkin(8) + (-t17 * t277 + t18 * t278 - t344 + t374) * mrSges(7,2) + (t277 * t28 + t278 * t29 - t344 + t375) * mrSges(6,3) + (-m(7) * t261 + (-m(7) * t203 + t376) * t271) * g(1) + (t100 * t351 + t352 * t373) * qJD(5) / 0.2e1 + (t293 + t338) * t72 + (-Ifges(7,3) * t330 + t329 * t365 + t347 + t377) * t178 + (-t277 / 0.2e1 + t227) * t360 + (-t286 + t288) * t332 + (t207 * t327 + t204 * t328 + Ifges(5,1) * t321 - t130 * mrSges(5,2) - t357 / 0.2e1 + t351 * t326 + t352 * t324 - t353) * t118 + t358 * t57 + (-t291 + t361) * t321 + (t363 / 0.2e1 + t366 * t329 + t367 * t332) * t175 + (-t28 * mrSges(6,1) + t17 * mrSges(7,1) + t29 * mrSges(6,2) - t18 * mrSges(7,3) + Ifges(6,6) * t327 + Ifges(7,6) * t328 - Ifges(5,2) * t322 - t130 * mrSges(5,1) + t356 / 0.2e1 + t366 * t326 + t364 * t324) * t119 + t77 * t320 + (t131 * t5 - t17 * t25 - t18 * t24 + t358 * t30) * m(7) + (-pkin(4) * t27 - g(1) * t261 - t28 * t36 - t29 * t37) * m(6) + (t157 * t368 + t354 * t158 + t215) * t310 + (t157 * t354 - t158 * t368 + t214) * g(3) + (t78 + t111) * t322 + (t239 + t353 - (-t207 / 0.2e1 + t204 / 0.2e1) * t194) * qJD(5) + t287 * t330 + t289 * t331 + Ifges(5,3) * qJDD(4) + t131 * t15 - t215 * t311 + t42 * t228 - pkin(4) * t16 - t31 * mrSges(5,2) + t32 * mrSges(5,1) - t24 * t67 - t37 * t68 - t36 * t69 - t25 * t70 + Ifges(5,5) * t89 + Ifges(5,6) * t90; (-Ifges(6,2) * t100 + t360 + t95) * t327 + (-t30 * t56 - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t18 + t202 * t309 - g(1) * (-pkin(5) * t113 + qJ(6) * t114) - g(2) * (pkin(5) * t115 - qJ(6) * t116)) * m(7) + (t212 + t210) * t309 + (-t370 * t115 - t369 * t116) * g(2) + t45 * t325 - t17 * t315 + (Ifges(7,3) * t100 + t303) * t328 + (-m(7) * t18 - t298 + t314) * t28 + t18 * t295 + (-t100 * t365 + t194 * t366) * t324 + (t194 * t367 - t290 + t42 + t94) * t326 - t30 * (mrSges(7,1) * t100 - mrSges(7,3) * t194) - t64 * (mrSges(6,1) * t100 + mrSges(6,2) * t194) + (-m(7) * t17 + t292 + t297) * t29 + (t370 * t113 + t369 * t114) * g(1) + t335 + t348 + qJ(6) * t19 - pkin(5) * t21 - t56 * t57 + qJD(6) * t67; t100 * t57 - t373 * t67 + (-g(1) * t113 + g(2) * t115 + t30 * t100 - t175 * t309 - t18 * t373 + t2) * m(7) + t21;];
tau  = t8;
