% Calculate vector of inverse dynamics joint torques for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:23
% EndTime: 2019-12-31 17:29:44
% DurationCPUTime: 11.72s
% Computational Cost: add. (4699->566), mult. (12025->818), div. (0->0), fcn. (9078->10), ass. (0->263)
t184 = cos(qJ(3));
t296 = pkin(7) * t184;
t351 = m(5) + m(4);
t181 = sin(qJ(2));
t177 = sin(pkin(4));
t259 = qJD(1) * t177;
t237 = t181 * t259;
t180 = sin(qJ(3));
t255 = qJD(3) * t180;
t185 = cos(qJ(2));
t178 = cos(pkin(4));
t258 = qJD(1) * t178;
t245 = pkin(1) * t258;
t128 = -pkin(6) * t237 + t185 * t245;
t198 = (pkin(2) * t181 - pkin(7) * t185) * t177;
t129 = qJD(1) * t198;
t79 = t184 * t128 + t180 * t129;
t362 = t255 * pkin(7) + pkin(8) * t237 + t79;
t174 = t178 * t181 * pkin(1);
t218 = pkin(3) * t180 - pkin(8) * t184;
t268 = t177 * t185;
t361 = -qJD(4) * t296 + t218 * qJD(3) - (t174 + (pkin(6) + t218) * t268) * qJD(1);
t179 = sin(qJ(4));
t183 = cos(qJ(4));
t214 = t179 * mrSges(5,1) + t183 * mrSges(5,2);
t360 = -t214 - mrSges(4,3);
t215 = -mrSges(5,1) * t183 + mrSges(5,2) * t179;
t190 = m(5) * pkin(3) - t215;
t216 = mrSges(4,1) * t184 - mrSges(4,2) * t180;
t243 = m(5) * pkin(8) + mrSges(5,3);
t359 = -t180 * t243 - t184 * t190 - t216;
t247 = qJDD(1) * t178;
t249 = qJD(1) * qJD(2);
t358 = -pkin(6) * t177 * t249 + pkin(1) * t247;
t248 = qJDD(1) * t177;
t357 = pkin(6) * t248 + qJD(2) * t245;
t263 = t184 * t185;
t104 = (-t179 * t263 + t181 * t183) * t259;
t254 = qJD(3) * t184;
t356 = t179 * t254 + t104;
t355 = pkin(2) * t351 + mrSges(3,1) - t359;
t167 = qJD(2) + t258;
t101 = -t167 * pkin(2) - t128;
t111 = t167 * t184 - t180 * t237;
t112 = t167 * t180 + t184 * t237;
t44 = -t111 * pkin(3) - t112 * pkin(8) + t101;
t236 = t185 * t259;
t154 = qJD(3) - t236;
t261 = pkin(6) * t268 + t174;
t131 = t261 * qJD(1);
t102 = t167 * pkin(7) + t131;
t106 = (-pkin(2) * t185 - pkin(7) * t181 - pkin(1)) * t259;
t57 = t102 * t184 + t106 * t180;
t46 = pkin(8) * t154 + t57;
t17 = t179 * t44 + t183 * t46;
t347 = t17 * mrSges(5,2);
t16 = -t179 * t46 + t183 * t44;
t348 = t16 * mrSges(5,1);
t288 = Ifges(4,4) * t112;
t341 = Ifges(4,6) * t154;
t342 = Ifges(4,2) * t111;
t51 = t288 + t341 + t342;
t354 = -t347 + t348 - t51 / 0.2e1;
t301 = cos(qJ(1));
t238 = t301 * t185;
t182 = sin(qJ(1));
t265 = t181 * t182;
t142 = -t178 * t238 + t265;
t239 = t301 * t181;
t264 = t182 * t185;
t143 = t178 * t239 + t264;
t240 = t177 * t301;
t93 = t143 * t184 - t180 * t240;
t353 = -t142 * t183 + t179 * t93;
t352 = -t142 * t179 - t183 * t93;
t134 = (-qJDD(1) * t185 + t181 * t249) * t177;
t125 = qJDD(3) + t134;
t302 = t125 / 0.2e1;
t135 = (qJDD(1) * t181 + t185 * t249) * t177;
t166 = qJDD(2) + t247;
t59 = -qJD(3) * t112 - t135 * t180 + t166 * t184;
t312 = t59 / 0.2e1;
t58 = qJD(3) * t111 + t135 * t184 + t166 * t180;
t313 = t58 / 0.2e1;
t323 = Ifges(4,1) * t313 + Ifges(4,4) * t312 + Ifges(4,5) * t302;
t82 = -t112 * t179 + t154 * t183;
t22 = qJD(4) * t82 + t125 * t179 + t183 * t58;
t322 = t22 / 0.2e1;
t83 = t112 * t183 + t154 * t179;
t23 = -qJD(4) * t83 + t125 * t183 - t179 * t58;
t321 = t23 / 0.2e1;
t55 = qJDD(4) - t59;
t314 = t55 / 0.2e1;
t84 = t358 * t181 + t357 * t185;
t69 = pkin(7) * t166 + t84;
t77 = -pkin(1) * t248 + pkin(2) * t134 - pkin(7) * t135;
t14 = -t102 * t255 + t106 * t254 + t180 * t77 + t184 * t69;
t11 = pkin(8) * t125 + t14;
t85 = -t357 * t181 + t358 * t185;
t70 = -t166 * pkin(2) - t85;
t13 = -t59 * pkin(3) - t58 * pkin(8) + t70;
t1 = qJD(4) * t16 + t11 * t183 + t13 * t179;
t350 = t1 * mrSges(5,2);
t2 = -qJD(4) * t17 - t11 * t179 + t13 * t183;
t349 = t2 * mrSges(5,1);
t150 = -pkin(3) * t184 - pkin(8) * t180 - pkin(2);
t253 = qJD(4) * t179;
t346 = -t150 * t253 + t362 * t179 + t361 * t183;
t251 = qJD(4) * t183;
t345 = t150 * t251 + t361 * t179 - t362 * t183;
t38 = -mrSges(5,1) * t82 + mrSges(5,2) * t83;
t290 = mrSges(4,3) * t112;
t87 = mrSges(4,1) * t154 - t290;
t344 = -t87 + t38;
t107 = Ifges(4,4) * t111;
t343 = Ifges(4,5) * t154;
t340 = mrSges(4,1) + t190;
t228 = mrSges(4,2) - t243;
t227 = mrSges(3,3) * t237;
t339 = -mrSges(3,1) * t167 - mrSges(4,1) * t111 + mrSges(4,2) * t112 + t227;
t338 = t180 * t251 + t356;
t105 = (t179 * t181 + t183 * t263) * t259;
t252 = qJD(4) * t180;
t337 = t179 * t252 - t183 * t254 + t105;
t121 = pkin(7) * t178 + t261;
t271 = t177 * t181;
t262 = pkin(2) * t268 + pkin(7) * t271;
t299 = pkin(1) * t177;
t122 = -t262 - t299;
t336 = t184 * t121 + t180 * t122;
t168 = pkin(6) * t271;
t298 = pkin(1) * t185;
t147 = t178 * t298 - t168;
t335 = t1 * t183 - t179 * t2;
t108 = qJD(4) - t111;
t289 = Ifges(3,4) * t181;
t334 = -t181 * (Ifges(3,1) * t185 - t289) / 0.2e1 + pkin(1) * (mrSges(3,1) * t181 + mrSges(3,2) * t185);
t15 = -t57 * qJD(3) - t180 * t69 + t184 * t77;
t331 = t15 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,5) * t58 + Ifges(4,6) * t59 + Ifges(4,3) * t125;
t130 = qJD(2) * t198;
t132 = t147 * qJD(2);
t37 = -qJD(3) * t336 + t130 * t184 - t132 * t180;
t330 = t85 * mrSges(3,1) - t84 * mrSges(3,2) + Ifges(3,5) * t135 - Ifges(3,6) * t134 + Ifges(3,3) * t166;
t329 = t349 - t350;
t327 = -t351 * pkin(7) + mrSges(3,2) + t360;
t29 = Ifges(5,5) * t83 + Ifges(5,6) * t82 + Ifges(5,3) * t108;
t326 = -t57 * mrSges(4,3) + t29 / 0.2e1 + t354;
t4 = Ifges(5,4) * t22 + Ifges(5,2) * t23 + Ifges(5,6) * t55;
t325 = t4 / 0.2e1;
t324 = Ifges(5,1) * t322 + Ifges(5,4) * t321 + Ifges(5,5) * t314;
t300 = Ifges(5,4) * t83;
t30 = Ifges(5,2) * t82 + Ifges(5,6) * t108 + t300;
t319 = -t30 / 0.2e1;
t318 = t30 / 0.2e1;
t80 = Ifges(5,4) * t82;
t31 = Ifges(5,1) * t83 + Ifges(5,5) * t108 + t80;
t317 = -t31 / 0.2e1;
t316 = t31 / 0.2e1;
t311 = -t82 / 0.2e1;
t310 = t82 / 0.2e1;
t309 = -t83 / 0.2e1;
t308 = t83 / 0.2e1;
t307 = -t108 / 0.2e1;
t306 = t108 / 0.2e1;
t303 = t112 / 0.2e1;
t291 = mrSges(4,3) * t111;
t287 = Ifges(4,4) * t180;
t286 = Ifges(4,4) * t184;
t285 = Ifges(5,4) * t179;
t284 = Ifges(5,4) * t183;
t283 = Ifges(3,6) * t167;
t282 = t111 * Ifges(4,6);
t281 = t112 * Ifges(4,5);
t12 = -pkin(3) * t125 - t15;
t280 = t12 * t180;
t279 = t14 * t184;
t278 = t154 * Ifges(4,3);
t277 = t167 * Ifges(3,5);
t276 = t180 * mrSges(4,3);
t275 = t111 * t179;
t274 = t111 * t183;
t270 = t177 * t182;
t269 = t177 * t184;
t267 = t179 * t180;
t266 = t180 * t183;
t260 = t301 * pkin(1) + pkin(6) * t270;
t256 = qJD(2) * t177;
t3 = Ifges(5,5) * t22 + Ifges(5,6) * t23 + Ifges(5,3) * t55;
t235 = t181 * t256;
t234 = t185 * t256;
t231 = t254 / 0.2e1;
t230 = -t252 / 0.2e1;
t229 = -pkin(1) * t182 + pkin(6) * t240;
t92 = -t143 * t180 - t184 * t240;
t226 = mrSges(3,3) * t236;
t141 = t178 * t180 + t181 * t269;
t197 = -t141 * t183 + t179 * t268;
t90 = -t141 * t179 - t183 * t268;
t219 = mrSges(5,1) * t90 + mrSges(5,2) * t197;
t140 = -t178 * t184 + t180 * t271;
t217 = t140 * mrSges(4,1) + t141 * mrSges(4,2);
t213 = Ifges(4,1) * t184 - t287;
t212 = Ifges(5,1) * t183 - t285;
t211 = Ifges(5,1) * t179 + t284;
t210 = -Ifges(4,2) * t180 + t286;
t209 = -Ifges(5,2) * t179 + t284;
t208 = Ifges(5,2) * t183 + t285;
t207 = Ifges(4,5) * t184 - Ifges(4,6) * t180;
t206 = Ifges(5,5) * t183 - Ifges(5,6) * t179;
t205 = Ifges(5,5) * t179 + Ifges(5,6) * t183;
t120 = t168 + (-pkin(2) - t298) * t178;
t60 = t140 * pkin(3) - t141 * pkin(8) + t120;
t62 = -pkin(8) * t268 + t336;
t25 = -t179 * t62 + t183 * t60;
t26 = t179 * t60 + t183 * t62;
t56 = -t102 * t180 + t106 * t184;
t75 = -t180 * t121 + t184 * t122;
t78 = -t128 * t180 + t129 * t184;
t144 = t178 * t264 + t239;
t145 = -t178 * t265 + t238;
t200 = t145 * pkin(2) + pkin(7) * t144 + t260;
t192 = -t143 * pkin(2) - pkin(7) * t142 + t229;
t36 = -t121 * t255 + t122 * t254 + t180 * t130 + t184 * t132;
t133 = t261 * qJD(2);
t164 = Ifges(3,4) * t236;
t146 = (-mrSges(3,1) * t185 + mrSges(3,2) * t181) * t177;
t127 = -t167 * mrSges(3,2) + t226;
t119 = t150 * t179 + t183 * t296;
t118 = t150 * t183 - t179 * t296;
t100 = Ifges(3,1) * t237 + t164 + t277;
t99 = t283 + (t185 * Ifges(3,2) + t289) * t259;
t97 = t145 * t184 + t180 * t270;
t96 = t145 * t180 - t182 * t269;
t89 = -qJD(3) * t140 + t184 * t234;
t88 = qJD(3) * t141 + t180 * t234;
t86 = -mrSges(4,2) * t154 + t291;
t74 = pkin(3) * t112 - pkin(8) * t111;
t66 = t144 * t179 + t183 * t97;
t65 = t144 * t183 - t179 * t97;
t63 = -pkin(3) * t237 - t78;
t61 = pkin(3) * t268 - t75;
t52 = Ifges(4,1) * t112 + t107 + t343;
t50 = t278 + t281 + t282;
t48 = mrSges(5,1) * t108 - mrSges(5,3) * t83;
t47 = -mrSges(5,2) * t108 + mrSges(5,3) * t82;
t45 = -pkin(3) * t154 - t56;
t43 = qJD(4) * t90 + t179 * t235 + t89 * t183;
t42 = qJD(4) * t197 - t89 * t179 + t183 * t235;
t41 = t88 * pkin(3) - t89 * pkin(8) + t133;
t40 = -mrSges(4,2) * t125 + mrSges(4,3) * t59;
t39 = mrSges(4,1) * t125 - mrSges(4,3) * t58;
t33 = -pkin(3) * t235 - t37;
t32 = pkin(8) * t235 + t36;
t28 = t179 * t74 + t183 * t56;
t27 = -t179 * t56 + t183 * t74;
t24 = -mrSges(4,1) * t59 + mrSges(4,2) * t58;
t18 = t58 * Ifges(4,4) + t59 * Ifges(4,2) + t125 * Ifges(4,6);
t10 = -mrSges(5,2) * t55 + mrSges(5,3) * t23;
t9 = mrSges(5,1) * t55 - mrSges(5,3) * t22;
t8 = -mrSges(5,1) * t23 + mrSges(5,2) * t22;
t7 = -qJD(4) * t26 - t179 * t32 + t183 * t41;
t6 = qJD(4) * t25 + t179 * t41 + t183 * t32;
t5 = [(-t342 / 0.2e1 - t341 / 0.2e1 + Ifges(5,6) * t310 - Ifges(4,4) * t303 + Ifges(5,3) * t306 + Ifges(5,5) * t308 + t101 * mrSges(4,1) + t326) * t88 + (t107 / 0.2e1 + t343 / 0.2e1 + Ifges(4,1) * t303 - t56 * mrSges(4,3) + t52 / 0.2e1 + t101 * mrSges(4,2)) * t89 + t339 * t133 + ((-g(1) * t301 - g(2) * t182) * mrSges(3,3) + (-mrSges(3,1) * t134 - mrSges(3,2) * t135 + (m(3) * t299 - t146) * qJDD(1)) * pkin(1) + (-t85 * mrSges(3,3) + Ifges(3,1) * t135 - Ifges(3,4) * t134 + Ifges(3,5) * t166) * t181 + (t84 * mrSges(3,3) + Ifges(3,4) * t135 - Ifges(3,2) * t134 + Ifges(3,6) * t166 - t331) * t185 + ((t277 / 0.2e1 - t128 * mrSges(3,3) + t100 / 0.2e1) * t185 + (-t283 / 0.2e1 - t131 * mrSges(3,3) - t57 * mrSges(4,2) + t56 * mrSges(4,1) + t282 / 0.2e1 + t281 / 0.2e1 + t278 / 0.2e1 - t99 / 0.2e1 + t50 / 0.2e1) * t181 + (t185 * (Ifges(3,4) * t185 - Ifges(3,2) * t181) / 0.2e1 - t334) * t259) * qJD(2)) * t177 + (t182 * mrSges(2,1) + t301 * mrSges(2,2) - m(4) * t192 + t93 * mrSges(4,1) - m(3) * t229 + t143 * mrSges(3,1) - t142 * mrSges(3,2) - m(5) * (-pkin(3) * t93 + t192) - t352 * mrSges(5,1) - t353 * mrSges(5,2) + t228 * t92) * g(1) + (-t15 * mrSges(4,3) + 0.2e1 * t323) * t141 + (-Ifges(4,2) * t312 - Ifges(4,4) * t313 + Ifges(5,3) * t314 + Ifges(5,6) * t321 + Ifges(5,5) * t322 - Ifges(4,6) * t302 - t14 * mrSges(4,3) + t3 / 0.2e1 - t18 / 0.2e1 + t329) * t140 + t330 * t178 + m(4) * (t101 * t133 + t120 * t70 + t14 * t336 + t15 * t75 + t36 * t57 + t37 * t56) + t336 * t40 - t197 * t324 + (t1 * t90 - t16 * t43 + t17 * t42 + t197 * t2) * mrSges(5,3) + (-Ifges(5,1) * t197 + Ifges(5,4) * t90) * t322 + (-Ifges(5,4) * t197 + Ifges(5,2) * t90) * t321 + (-Ifges(5,5) * t197 + Ifges(5,6) * t90) * t314 + m(5) * (t1 * t26 + t12 * t61 + t16 * t7 + t17 * t6 + t2 * t25 + t33 * t45) + (-t301 * mrSges(2,1) + t182 * mrSges(2,2) - m(4) * t200 - t97 * mrSges(4,1) - m(3) * t260 - t145 * mrSges(3,1) + t144 * mrSges(3,2) - m(5) * (pkin(3) * t97 + t200) - t66 * mrSges(5,1) - t65 * mrSges(5,2) + t228 * t96) * g(2) + (Ifges(5,1) * t43 + Ifges(5,4) * t42) * t308 + (Ifges(5,4) * t43 + Ifges(5,2) * t42) * t310 + (g(1) * t142 - g(2) * t144) * mrSges(4,3) + (Ifges(5,5) * t43 + Ifges(5,6) * t42) * t306 + t90 * t325 + t43 * t316 + t42 * t318 + m(3) * (-t128 * t133 + t131 * t132 + t147 * t85 + t261 * t84) + t261 * (-mrSges(3,2) * t166 - mrSges(3,3) * t134) + t70 * t217 + t25 * t9 + t26 * t10 - t12 * t219 + t33 * t38 + t45 * (-mrSges(5,1) * t42 + mrSges(5,2) * t43) + t6 * t47 + t7 * t48 + t61 * t8 + t75 * t39 + Ifges(2,3) * qJDD(1) + t36 * t86 + t37 * t87 + t120 * t24 + t132 * t127 + t147 * (mrSges(3,1) * t166 - mrSges(3,3) * t135); -t184 * t349 + t345 * t47 + t346 * t48 + (-m(4) * t101 + t227 - t339) * t131 + (mrSges(5,1) * t338 - mrSges(5,2) * t337) * t45 + (-t1 * t267 + t16 * t337 - t17 * t338 - t2 * t266) * mrSges(5,3) + t330 + (Ifges(5,5) * t309 + Ifges(5,6) * t311 + Ifges(5,3) * t307 - t354) * t180 * t236 + t334 * qJD(1) ^ 2 * t177 ^ 2 + t154 * t101 * (mrSges(4,1) * t180 + mrSges(4,2) * t184) - ((-Ifges(3,2) * t237 + t180 * t29 + t184 * t52 + t100 + t164) * t185 + t167 * (Ifges(3,5) * t185 - Ifges(3,6) * t181) + t181 * t50 + t154 * (Ifges(4,3) * t181 + t185 * t207) + t112 * (Ifges(4,5) * t181 + t185 * t213) + t111 * (Ifges(4,6) * t181 + t185 * t210)) * t259 / 0.2e1 + t184 * t350 + t99 * t237 / 0.2e1 + (Ifges(5,1) * t105 + Ifges(5,4) * t104) * t309 + (t146 - t351 * t262 + (t360 * t181 + t359 * t185) * t177) * g(3) + (-t254 * t56 + t279) * mrSges(4,3) + (Ifges(5,4) * t105 + Ifges(5,2) * t104) * t311 - t15 * t276 + (t230 * t30 + t231 * t31) * t183 - t4 * t267 / 0.2e1 + (-t208 * t252 + (Ifges(5,6) * t180 + t184 * t209) * qJD(3)) * t310 + (Ifges(4,2) * t184 + t287) * t312 + (Ifges(4,1) * t180 + t286) * t313 + (-Ifges(5,3) * t184 + t180 * t206) * t314 + t105 * t317 + (-Ifges(5,6) * t184 + t180 * t209) * t321 + (-Ifges(5,5) * t184 + t180 * t212) * t322 + t180 * t323 + t266 * t324 + t40 * t296 + (Ifges(4,5) * t180 + Ifges(4,6) * t184) * t302 + (-t205 * t252 + (Ifges(5,3) * t180 + t184 * t206) * qJD(3)) * t306 + (-t211 * t252 + (Ifges(5,5) * t180 + t184 * t212) * qJD(3)) * t308 + (Ifges(5,5) * t105 + Ifges(5,6) * t104) * t307 + (t355 * t142 + t143 * t327) * g(2) + (t355 * t144 + t145 * t327) * g(1) + t356 * t319 + (-t56 * (mrSges(4,1) * t181 - mrSges(4,3) * t263) - t57 * (-mrSges(4,2) * t181 - t185 * t276)) * t259 + t179 * t31 * t230 + (t226 - t127) * t128 + (t111 * t210 + t112 * t213 + t154 * t207) * qJD(3) / 0.2e1 - t70 * t216 - pkin(2) * t24 + t52 * t231 + t214 * t280 - t63 * t38 - t79 * t86 - t78 * t87 + t118 * t9 + t119 * t10 + (t1 * t119 + t118 * t2 + t346 * t16 + t345 * t17 - t45 * t63) * m(5) + t326 * t255 + (-pkin(2) * t70 - t56 * t78 - t57 * t79) * m(4) + (m(5) * (t254 * t45 + t280) - t86 * t255 + t344 * t254 + (-t39 + t8) * t180 + m(4) * (t279 - t15 * t180 + (-t57 * t180 - t56 * t184) * qJD(3))) * pkin(7) - t184 * t3 / 0.2e1 + t184 * t18 / 0.2e1; -t112 * t348 + (-m(5) * t45 + t290 - t344) * t57 + (t228 * t97 + t340 * t96) * g(1) + (t228 * t93 - t340 * t92) * g(2) + ((-t253 + t275) * t17 + (-t251 + t274) * t16 + t335) * mrSges(5,3) + (m(5) * ((-t16 * t183 - t17 * t179) * qJD(4) + t335) - t48 * t251 - t47 * t253 - t179 * t9 + t183 * t10) * pkin(8) + t331 + t108 * t45 * t214 + t112 * t347 + (-pkin(3) * t12 - t16 * t27 - t17 * t28) * m(5) + t179 * t324 + t183 * t325 + (Ifges(5,5) * t112 + t111 * t212) * t309 + (Ifges(5,6) * t112 + t111 * t209) * t311 + t205 * t314 + t251 * t316 + t274 * t317 + t275 * t318 + t253 * t319 + t208 * t321 + t211 * t322 + t51 * t303 + (Ifges(5,3) * t112 + t111 * t206) * t307 + (t140 * t190 - t141 * t243 + t217) * g(3) + (t291 - t86) * t56 + (t108 * t206 + t209 * t82 + t212 * t83) * qJD(4) / 0.2e1 - (Ifges(4,1) * t111 - t288 + t29) * t112 / 0.2e1 - (-Ifges(4,2) * t112 + t107 + t52) * t111 / 0.2e1 - pkin(3) * t8 + t12 * t215 - t28 * t47 - t27 * t48 - t101 * (mrSges(4,1) * t112 + mrSges(4,2) * t111) - t154 * (Ifges(4,5) * t111 - Ifges(4,6) * t112) / 0.2e1; -t45 * (mrSges(5,1) * t83 + mrSges(5,2) * t82) + (Ifges(5,1) * t82 - t300) * t309 + t30 * t308 + (Ifges(5,5) * t82 - Ifges(5,6) * t83) * t307 - t16 * t47 + t17 * t48 - g(1) * (mrSges(5,1) * t65 - mrSges(5,2) * t66) - g(2) * (-t353 * mrSges(5,1) + t352 * mrSges(5,2)) - g(3) * t219 + (t16 * t82 + t17 * t83) * mrSges(5,3) + t3 + (-Ifges(5,2) * t83 + t31 + t80) * t311 + t329;];
tau = t5;
