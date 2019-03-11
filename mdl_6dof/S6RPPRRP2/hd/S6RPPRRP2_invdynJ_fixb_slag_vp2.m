% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:17
% EndTime: 2019-03-09 02:00:32
% DurationCPUTime: 12.25s
% Computational Cost: add. (6301->552), mult. (13811->697), div. (0->0), fcn. (9701->14), ass. (0->234)
t362 = -mrSges(6,3) - mrSges(7,2);
t349 = Ifges(6,1) + Ifges(7,1);
t361 = -Ifges(6,4) + Ifges(7,5);
t348 = Ifges(7,4) + Ifges(6,5);
t347 = Ifges(6,6) - Ifges(7,6);
t346 = Ifges(6,3) + Ifges(7,2);
t173 = sin(qJ(5));
t176 = cos(qJ(5));
t207 = t176 * mrSges(7,1) + t173 * mrSges(7,3);
t209 = mrSges(6,1) * t176 - mrSges(6,2) * t173;
t360 = -t207 - t209;
t169 = sin(pkin(9));
t146 = pkin(1) * t169 + qJ(3);
t131 = qJD(1) * qJD(3) + qJDD(1) * t146;
t168 = sin(pkin(10));
t170 = cos(pkin(10));
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t133 = t168 * t174 - t177 * t170;
t125 = t133 * qJD(1);
t355 = qJD(5) + t125;
t134 = t168 * t177 + t170 * t174;
t126 = t134 * qJD(1);
t190 = t176 * qJD(4) - t126 * t173;
t127 = t133 * qJD(4);
t91 = -qJD(1) * t127 + qJDD(1) * t134;
t54 = qJD(5) * t190 + qJDD(4) * t173 + t176 * t91;
t312 = t54 / 0.2e1;
t107 = qJD(4) * t173 + t126 * t176;
t55 = qJD(5) * t107 - t176 * qJDD(4) + t173 * t91;
t310 = t55 / 0.2e1;
t128 = t134 * qJD(4);
t92 = -qJD(1) * t128 - qJDD(1) * t133;
t88 = qJDD(5) - t92;
t309 = t88 / 0.2e1;
t212 = -mrSges(4,1) * t170 + mrSges(4,2) * t168;
t359 = -m(4) * pkin(2) - mrSges(3,1) + t212;
t358 = t168 ^ 2 + t170 ^ 2;
t109 = t168 * qJDD(2) + t170 * t131;
t240 = qJDD(1) * t170;
t100 = pkin(7) * t240 + t109;
t143 = t146 * qJD(1);
t156 = t170 * qJD(2);
t267 = pkin(7) * qJD(1);
t104 = t156 + (-t143 - t267) * t168;
t117 = t168 * qJD(2) + t170 * t143;
t105 = t170 * t267 + t117;
t246 = qJD(4) * t177;
t247 = qJD(4) * t174;
t154 = t170 * qJDD(2);
t99 = t154 + (-pkin(7) * qJDD(1) - t131) * t168;
t22 = t177 * t100 + t104 * t246 - t105 * t247 + t174 * t99;
t20 = qJDD(4) * pkin(8) + t22;
t243 = qJD(5) * t176;
t244 = qJD(5) * t173;
t151 = pkin(3) * t170 + pkin(2);
t171 = cos(pkin(9));
t296 = pkin(1) * t171;
t141 = -t151 - t296;
t120 = qJDD(1) * t141 + qJDD(3);
t39 = -pkin(4) * t92 - pkin(8) * t91 + t120;
t61 = t104 * t174 + t105 * t177;
t58 = qJD(4) * pkin(8) + t61;
t123 = qJD(1) * t141 + qJD(3);
t68 = pkin(4) * t125 - pkin(8) * t126 + t123;
t3 = t173 * t39 + t176 * t20 + t68 * t243 - t244 * t58;
t25 = t173 * t68 + t176 * t58;
t4 = -qJD(5) * t25 - t173 * t20 + t176 * t39;
t215 = -t173 * t4 + t176 * t3;
t24 = -t173 * t58 + t176 * t68;
t357 = -t24 * t243 - t25 * t244 + t215;
t17 = -pkin(5) * t355 + qJD(6) - t24;
t18 = qJ(6) * t355 + t25;
t1 = qJ(6) * t88 + qJD(6) * t355 + t3;
t2 = -pkin(5) * t88 + qJDD(6) - t4;
t216 = t1 * t176 + t173 * t2;
t356 = t17 * t243 - t18 * t244 + t216;
t30 = mrSges(6,1) * t88 - mrSges(6,3) * t54;
t31 = -t88 * mrSges(7,1) + t54 * mrSges(7,2);
t284 = -t30 + t31;
t29 = -mrSges(7,2) * t55 + mrSges(7,3) * t88;
t32 = -mrSges(6,2) * t88 - mrSges(6,3) * t55;
t285 = t29 + t32;
t354 = t284 * t173 + t285 * t176;
t353 = -m(4) - m(3);
t352 = -m(6) - m(7);
t351 = mrSges(6,1) + mrSges(7,1);
t350 = mrSges(6,2) - mrSges(7,3);
t345 = t348 * t88 + t349 * t54 + t361 * t55;
t16 = mrSges(6,1) * t55 + mrSges(6,2) * t54;
t344 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t91 + t16;
t343 = t107 * t348 + t190 * t347 + t346 * t355;
t103 = Ifges(6,4) * t190;
t270 = Ifges(7,5) * t190;
t342 = t107 * t349 + t348 * t355 + t103 - t270;
t280 = mrSges(7,2) * t190;
t69 = mrSges(7,3) * t355 + t280;
t276 = mrSges(6,3) * t190;
t70 = -mrSges(6,2) * t355 + t276;
t283 = t69 + t70;
t275 = mrSges(6,3) * t107;
t71 = mrSges(6,1) * t355 - t275;
t279 = mrSges(7,2) * t107;
t72 = -mrSges(7,1) * t355 + t279;
t282 = t71 - t72;
t277 = mrSges(5,3) * t126;
t341 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t190 - mrSges(6,2) * t107 - t277;
t197 = pkin(5) * t173 - qJ(6) * t176;
t340 = -qJD(6) * t173 + t197 * t355 - t61;
t339 = Ifges(5,5) * qJD(4);
t338 = Ifges(5,6) * qJD(4);
t286 = pkin(7) + t146;
t129 = t286 * t168;
t130 = t286 * t170;
t337 = -t177 * t129 - t130 * t174;
t166 = pkin(10) + qJ(4);
t157 = sin(t166);
t159 = cos(t166);
t336 = t159 * pkin(4) + t157 * pkin(8);
t206 = t173 * mrSges(7,1) - t176 * mrSges(7,3);
t208 = mrSges(6,1) * t173 + mrSges(6,2) * t176;
t60 = t104 * t177 - t174 * t105;
t57 = -qJD(4) * pkin(4) - t60;
t28 = -pkin(5) * t190 - qJ(6) * t107 + t57;
t335 = t28 * t206 + t57 * t208;
t334 = -t173 * t347 + t176 * t348;
t269 = Ifges(7,5) * t173;
t272 = Ifges(6,4) * t173;
t333 = t176 * t349 + t269 - t272;
t259 = t127 * t173;
t186 = t134 * t243 - t259;
t329 = t346 * t88 - t347 * t55 + t348 * t54;
t108 = -t131 * t168 + t154;
t328 = -t108 * t168 + t109 * t170;
t325 = -m(5) + t352;
t322 = pkin(8) * t352;
t321 = Ifges(7,5) * t312 + Ifges(7,6) * t309 - t54 * Ifges(6,4) / 0.2e1 - t88 * Ifges(6,6) / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t310;
t198 = pkin(5) * t176 + qJ(6) * t173;
t142 = -pkin(4) - t198;
t320 = (mrSges(5,2) + t362) * t159 + (m(6) * pkin(4) - m(7) * t142 + mrSges(5,1) - t360) * t157;
t211 = mrSges(5,1) * t159 - mrSges(5,2) * t157;
t319 = t157 * t362 - t211;
t76 = pkin(4) * t133 - pkin(8) * t134 + t141;
t84 = -t129 * t174 + t130 * t177;
t281 = t173 * t76 + t176 * t84;
t63 = -t133 * qJD(3) + qJD(4) * t337;
t90 = pkin(4) * t128 + pkin(8) * t127;
t9 = -qJD(5) * t281 - t173 * t63 + t176 * t90;
t318 = m(7) * pkin(5) + t351;
t317 = -m(6) * t57 + t341;
t316 = m(7) * qJ(6) - t350;
t315 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t314 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t311 = -t55 / 0.2e1;
t308 = t190 / 0.2e1;
t307 = -t190 / 0.2e1;
t306 = -t107 / 0.2e1;
t305 = t107 / 0.2e1;
t304 = -t355 / 0.2e1;
t302 = t125 / 0.2e1;
t301 = -t126 / 0.2e1;
t300 = t126 / 0.2e1;
t297 = t173 / 0.2e1;
t175 = sin(qJ(1));
t295 = pkin(1) * t175;
t292 = g(3) * t157;
t178 = cos(qJ(1));
t163 = t178 * pkin(1);
t89 = pkin(4) * t126 + pkin(8) * t125;
t34 = t173 * t89 + t176 * t60;
t278 = mrSges(5,3) * t125;
t274 = Ifges(5,4) * t126;
t273 = Ifges(6,4) * t107;
t271 = Ifges(6,4) * t176;
t268 = Ifges(7,5) * t176;
t261 = t125 * t173;
t260 = t125 * t176;
t167 = qJ(1) + pkin(9);
t160 = cos(t167);
t255 = t157 * t160;
t158 = sin(t167);
t254 = t158 * t173;
t253 = t158 * t176;
t252 = t159 * t160;
t251 = t160 * t173;
t250 = t160 * t176;
t249 = t176 * t127;
t241 = qJDD(1) * t168;
t66 = -mrSges(7,1) * t190 - mrSges(7,3) * t107;
t237 = -t66 + t341;
t234 = m(4) - t325;
t102 = Ifges(7,5) * t107;
t43 = Ifges(7,6) * t355 - Ifges(7,3) * t190 + t102;
t229 = t43 * t297;
t152 = -pkin(2) - t296;
t222 = -t92 * mrSges(5,1) + t91 * mrSges(5,2);
t219 = -t244 / 0.2e1;
t218 = t243 / 0.2e1;
t172 = -pkin(7) - qJ(3);
t214 = t160 * t151 - t158 * t172 + t163;
t213 = -mrSges(4,1) * t240 + mrSges(4,2) * t241;
t203 = -Ifges(6,2) * t173 + t271;
t200 = Ifges(7,3) * t173 + t268;
t33 = -t173 * t60 + t176 * t89;
t40 = -t173 * t84 + t176 * t76;
t23 = -t174 * t100 - t104 * t247 - t105 * t246 + t177 * t99;
t191 = -(-t143 * t168 + t156) * t168 + t117 * t170;
t8 = t173 * t90 + t176 * t63 + t76 * t243 - t244 * t84;
t185 = t134 * t244 + t249;
t21 = -qJDD(4) * pkin(4) - t23;
t64 = qJD(3) * t134 + qJD(4) * t84;
t136 = qJDD(1) * t152 + qJDD(3);
t119 = Ifges(5,4) * t125;
t114 = -qJD(4) * mrSges(5,2) - t278;
t113 = t159 * t250 + t254;
t112 = t159 * t251 - t253;
t111 = t159 * t253 - t251;
t110 = t159 * t254 + t250;
t80 = t126 * Ifges(5,1) - t119 + t339;
t79 = -t125 * Ifges(5,2) + t274 + t338;
t78 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t92;
t65 = pkin(5) * t107 - qJ(6) * t190;
t46 = Ifges(6,2) * t190 + Ifges(6,6) * t355 + t273;
t42 = t134 * t197 - t337;
t36 = -pkin(5) * t133 - t40;
t35 = qJ(6) * t133 + t281;
t27 = -pkin(5) * t126 - t33;
t26 = qJ(6) * t126 + t34;
t15 = mrSges(7,1) * t55 - mrSges(7,3) * t54;
t14 = -t197 * t127 + (qJD(5) * t198 - qJD(6) * t176) * t134 + t64;
t7 = -pkin(5) * t128 - t9;
t6 = qJ(6) * t128 + qJD(6) * t133 + t8;
t5 = pkin(5) * t55 - qJ(6) * t54 - qJD(6) * t107 + t21;
t10 = [(Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t171 - 0.2e1 * mrSges(3,2) * t169 + m(3) * (t169 ^ 2 + t171 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) - t342 * t249 / 0.2e1 + t343 * t128 / 0.2e1 + (t348 * t128 - t349 * t185 + t186 * t361) * t305 + t17 * (-mrSges(7,1) * t128 - mrSges(7,2) * t185) + t57 * (mrSges(6,1) * t186 - mrSges(6,2) * t185) + t28 * (mrSges(7,1) * t186 + mrSges(7,3) * t185) + t25 * (-mrSges(6,2) * t128 - mrSges(6,3) * t186) + t18 * (-mrSges(7,2) * t186 + mrSges(7,3) * t128) + (t120 * mrSges(5,1) - t22 * mrSges(5,3) - Ifges(5,4) * t91 - Ifges(5,2) * t92 - Ifges(5,6) * qJDD(4) + Ifges(6,6) * t311 + Ifges(7,6) * t310 + t309 * t346 + t312 * t348 + t314 + t329 / 0.2e1) * t133 + (-Ifges(5,1) * t127 - Ifges(5,4) * t128) * t300 + (t127 * t60 - t128 * t61) * mrSges(5,3) + qJD(4) * (-Ifges(5,5) * t127 - Ifges(5,6) * t128) / 0.2e1 + t123 * (mrSges(5,1) * t128 - mrSges(5,2) * t127) - t125 * (-Ifges(5,4) * t127 - Ifges(5,2) * t128) / 0.2e1 + m(6) * (t24 * t9 + t25 * t8 + t281 * t3 + t4 * t40) + t281 * t32 - (-m(5) * t23 + m(6) * t21 + t344) * t337 + (t128 * t346 - t185 * t348 - t186 * t347) * t355 / 0.2e1 + (mrSges(2,1) * t175 + mrSges(2,2) * t178 - t353 * t295 + t325 * (-t160 * t172 - t295) + t318 * t111 + t316 * t110 + t315 * t160 + (m(5) * t151 + t352 * (-t151 - t336) - t319 - t359) * t158) * g(1) + t136 * t212 + t152 * t213 + ((-t1 * mrSges(7,2) - t3 * mrSges(6,3) + t321) * t173 + t342 * t219 + t120 * mrSges(5,2) - t23 * mrSges(5,3) + Ifges(5,1) * t91 + Ifges(5,4) * t92 + Ifges(5,5) * qJDD(4) + t200 * t310 + t203 * t311 + t206 * t5 + t208 * t21 + t218 * t43 + t309 * t334 + t312 * t333 + (t345 / 0.2e1 + mrSges(7,2) * t2 - mrSges(6,3) * t4) * t176) * t134 + (-m(5) * t214 - mrSges(2,1) * t178 + mrSges(2,2) * t175 + t362 * t255 + t352 * (pkin(4) * t252 + pkin(8) * t255 + t214) + t353 * t163 - t318 * t113 - t316 * t112 + (-t211 + t359) * t160 + t315 * t158) * g(2) + m(4) * (t191 * qJD(3) + t136 * t152 + t146 * t328) + m(7) * (t1 * t35 + t14 * t28 + t17 * t7 + t18 * t6 + t2 * t36 + t42 * t5) + t24 * (mrSges(6,1) * t128 + mrSges(6,3) * t185) + (Ifges(4,4) * t168 + Ifges(4,2) * t170) * t240 + (Ifges(4,1) * t168 + Ifges(4,4) * t170) * t241 + (t131 * t358 + t328) * mrSges(4,3) + t141 * t222 + m(5) * (t120 * t141 + t22 * t84 + t61 * t63) - t186 * t46 / 0.2e1 + (-m(5) * t60 - t317) * t64 + (-Ifges(7,5) * t185 + Ifges(7,6) * t128 + Ifges(7,3) * t186) * t307 + (-Ifges(6,4) * t185 - Ifges(6,2) * t186 + Ifges(6,6) * t128) * t308 - t127 * t229 + t35 * t29 + t36 * t31 + t40 * t30 + t42 * t15 + t14 * t66 + t6 * t69 + t8 * t70 + t9 * t71 + t7 * t72 + t84 * t78 + t63 * t114 - t127 * t80 / 0.2e1 - t128 * t79 / 0.2e1; m(3) * qJDD(2) + (t15 + t344) * t133 - t237 * t128 - (-t173 * t282 + t176 * t283 + t114) * t127 + (-m(3) - t234) * g(3) + m(5) * (-t127 * t61 - t128 * t60 - t133 * t23) + m(4) * (t108 * t170 + t109 * t168) + m(7) * (t128 * t28 + t133 * t5 - t17 * t259 - t18 * t249) + m(6) * (t128 * t57 + t133 * t21 + t24 * t259 - t249 * t25) + (t78 + (-t173 * t283 - t176 * t282) * qJD(5) + m(5) * t22 + m(7) * t356 + m(6) * t357 + t354) * t134; t125 * t114 - t358 * qJD(1) ^ 2 * mrSges(4,3) + t237 * t126 + (t283 * t355 - t284) * t176 + (-t282 * t355 + t285) * t173 + t213 + t222 + (-g(1) * t158 + g(2) * t160) * t234 + (t1 * t173 - t126 * t28 - t176 * t2 + t355 * (t17 * t173 + t176 * t18)) * m(7) + (-t126 * t57 + t173 * t3 + t176 * t4 + t355 * (-t24 * t173 + t25 * t176)) * m(6) + (t125 * t61 + t126 * t60 + t120) * m(5) + (-qJD(1) * t191 + t136) * m(4); (-t274 + t343) * t301 + t345 * t297 + (Ifges(6,2) * t311 - Ifges(7,3) * t310 + t309 * t347 - t321) * t176 + (-t24 * mrSges(6,1) + t17 * mrSges(7,1) + t25 * mrSges(6,2) - t18 * mrSges(7,3) - Ifges(5,2) * t302 + Ifges(6,6) * t307 + Ifges(7,6) * t308 + t338 / 0.2e1 - t123 * mrSges(5,1) + t348 * t306 + t346 * t304) * t126 + (t309 * t348 + t312 * t349) * t173 + t79 * t300 + (-pkin(4) * t21 - t24 * t33 - t25 * t34) * m(6) + (t142 * t5 - t17 * t27 - t18 * t26 + t340 * t28) * m(7) + (-t278 - t114) * t60 + (t107 * t333 + t334 * t355) * qJD(5) / 0.2e1 - t5 * t207 - t21 * t209 + (-t282 * t243 - t283 * t244 + ((-t173 * t25 - t176 * t24) * qJD(5) + t215) * m(6) + ((t17 * t176 - t173 * t18) * qJD(5) + t216) * m(7) + t354) * pkin(8) + (t160 * t320 + t252 * t322) * g(1) + t269 * t310 + t272 * t311 + (t17 * t260 - t18 * t261 + t356) * mrSges(7,2) + (-t24 * t260 - t25 * t261 + t357) * mrSges(6,3) - (t46 / 0.2e1 - t43 / 0.2e1) * t261 + (t218 + t260 / 0.2e1) * t342 - (Ifges(5,1) * t301 + t203 * t307 + t200 * t308 - t339 / 0.2e1 - t123 * mrSges(5,2) + t333 * t306 + t334 * t304 - t335) * t125 + t340 * t66 + (-(-t203 / 0.2e1 + t200 / 0.2e1) * t190 + t229 + t335) * qJD(5) + (t80 - t119) * t302 + (t352 * t336 + (-m(7) * t198 + t360) * t159 + t319) * g(3) + (t159 * t322 + t320) * g(2) * t158 + (-t268 + t271) * t312 + (t277 + t317) * t61 + t46 * t219 - pkin(4) * t16 - t22 * mrSges(5,2) + t23 * mrSges(5,1) - t26 * t69 - t34 * t70 - t33 * t71 - t27 * t72 + Ifges(5,5) * t91 + Ifges(5,6) * t92 + t142 * t15 + Ifges(5,3) * qJDD(4); (t351 * t110 + t350 * t111) * g(2) + (t208 + t206) * t292 + t18 * t279 - t17 * t280 + (t351 * t112 + t350 * t113) * g(1) + t46 * t305 + (Ifges(7,3) * t107 + t270) * t308 + (-m(7) * t18 + t276 - t283) * t24 + (-t107 * t347 + t190 * t348) * t304 + (t190 * t349 + t102 - t273 + t43) * t306 - t28 * (mrSges(7,1) * t107 - mrSges(7,3) * t190) - t57 * (mrSges(6,1) * t107 + mrSges(6,2) * t190) + (-m(7) * t17 + t275 + t282) * t25 + (-t28 * t65 - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t18 + t197 * t292 - g(2) * (-pkin(5) * t110 + qJ(6) * t111) - g(1) * (-pkin(5) * t112 + qJ(6) * t113)) * m(7) + t314 + t329 + (-Ifges(6,2) * t107 + t103 + t342) * t307 + qJ(6) * t29 - pkin(5) * t31 - t65 * t66 + qJD(6) * t69; t107 * t66 - t355 * t69 + (-g(1) * t112 - g(2) * t110 + t28 * t107 - t173 * t292 - t18 * t355 + t2) * m(7) + t31;];
tau  = t10;
