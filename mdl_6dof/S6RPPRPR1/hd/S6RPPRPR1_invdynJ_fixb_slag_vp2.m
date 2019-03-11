% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:08
% EndTime: 2019-03-09 01:39:28
% DurationCPUTime: 12.93s
% Computational Cost: add. (7596->580), mult. (17083->762), div. (0->0), fcn. (12630->18), ass. (0->253)
t198 = sin(pkin(11));
t201 = cos(pkin(11));
t206 = sin(qJ(6));
t209 = cos(qJ(6));
t222 = t198 * t206 - t201 * t209;
t155 = t222 * qJD(6);
t199 = sin(pkin(10));
t207 = sin(qJ(4));
t202 = cos(pkin(10));
t292 = cos(qJ(4));
t246 = t292 * t202;
t220 = -t207 * t199 + t246;
t153 = t220 * qJD(1);
t97 = t222 * t153;
t337 = t155 - t97;
t165 = t198 * t209 + t201 * t206;
t156 = t165 * qJD(6);
t96 = t165 * t153;
t336 = t156 - t96;
t150 = qJD(6) - t153;
t166 = t199 * t292 + t207 * t202;
t154 = t166 * qJD(1);
t137 = qJD(4) * t198 + t154 * t201;
t239 = t201 * qJD(4) - t154 * t198;
t77 = t137 * t209 + t206 * t239;
t291 = Ifges(7,4) * t77;
t351 = -t137 * t206 + t209 * t239;
t31 = Ifges(7,2) * t351 + t150 * Ifges(7,6) + t291;
t313 = t31 / 0.2e1;
t73 = Ifges(7,4) * t351;
t32 = t77 * Ifges(7,1) + t150 * Ifges(7,5) + t73;
t312 = t32 / 0.2e1;
t200 = sin(pkin(9));
t175 = pkin(1) * t200 + qJ(3);
t161 = qJD(1) * qJD(3) + qJDD(1) * t175;
t353 = t199 ^ 2 + t202 ^ 2;
t197 = qJ(1) + pkin(9);
t188 = sin(t197);
t191 = cos(t197);
t352 = g(1) * t191 + g(2) * t188;
t157 = t220 * qJD(4);
t118 = qJD(1) * t157 + qJDD(1) * t166;
t100 = qJDD(4) * t201 - t118 * t198;
t101 = qJDD(4) * t198 + t118 * t201;
t25 = qJD(6) * t351 + t100 * t206 + t101 * t209;
t315 = t25 / 0.2e1;
t26 = -qJD(6) * t77 + t100 * t209 - t101 * t206;
t314 = t26 / 0.2e1;
t350 = -m(4) - m(3);
t349 = -m(5) - m(7);
t172 = t175 * qJD(1);
t185 = t202 * qJD(2);
t278 = pkin(7) * qJD(1);
t139 = t185 + (-t172 - t278) * t199;
t146 = t199 * qJD(2) + t202 * t172;
t140 = t202 * t278 + t146;
t80 = t207 * t139 + t140 * t292;
t72 = qJD(4) * qJ(5) + t80;
t179 = pkin(3) * t202 + pkin(2);
t203 = cos(pkin(9));
t290 = pkin(1) * t203;
t169 = -t179 - t290;
t152 = qJD(1) * t169 + qJD(3);
t89 = -pkin(4) * t153 - qJ(5) * t154 + t152;
t42 = -t198 * t72 + t201 * t89;
t21 = -pkin(5) * t153 - pkin(8) * t137 + t42;
t43 = t198 * t89 + t201 * t72;
t28 = pkin(8) * t239 + t43;
t8 = -t206 * t28 + t209 * t21;
t348 = t8 * mrSges(7,1);
t9 = t206 * t21 + t209 * t28;
t347 = t9 * mrSges(7,2);
t305 = t100 / 0.2e1;
t304 = t101 / 0.2e1;
t158 = t166 * qJD(4);
t254 = qJDD(1) * t199;
t119 = qJD(1) * t158 - qJDD(1) * t246 + t207 * t254;
t117 = qJDD(6) + t119;
t303 = t117 / 0.2e1;
t302 = t119 / 0.2e1;
t346 = t42 * mrSges(6,1);
t345 = t43 * mrSges(6,2);
t340 = t239 * Ifges(6,6);
t341 = t137 * Ifges(6,5);
t344 = t77 * Ifges(7,5) + Ifges(7,6) * t351 - t153 * Ifges(6,3) + t150 * Ifges(7,3) + t340 + t341;
t284 = pkin(8) + qJ(5);
t170 = t284 * t198;
t171 = t284 * t201;
t129 = -t170 * t209 - t171 * t206;
t266 = t153 * t201;
t115 = pkin(4) * t154 - qJ(5) * t153;
t128 = t207 * t140;
t79 = t139 * t292 - t128;
t52 = t201 * t115 - t198 * t79;
t29 = pkin(5) * t154 - pkin(8) * t266 + t52;
t267 = t153 * t198;
t53 = t198 * t115 + t201 * t79;
t41 = -pkin(8) * t267 + t53;
t343 = -qJD(5) * t222 + qJD(6) * t129 - t206 * t29 - t209 * t41;
t130 = -t170 * t206 + t171 * t209;
t342 = -qJD(5) * t165 - qJD(6) * t130 + t206 * t41 - t209 * t29;
t231 = mrSges(6,1) * t198 + mrSges(6,2) * t201;
t70 = -qJD(4) * pkin(4) + qJD(5) - t79;
t339 = t70 * t231;
t282 = mrSges(5,3) * t154;
t338 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t239 - mrSges(6,2) * t137 - t282;
t54 = -t100 * mrSges(6,1) + t101 * mrSges(6,2);
t335 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t118 + t54;
t334 = Ifges(5,5) * qJD(4);
t333 = Ifges(5,6) * qJD(4);
t196 = pkin(10) + qJ(4);
t187 = sin(t196);
t332 = t187 * t352;
t190 = cos(t196);
t232 = -mrSges(6,1) * t201 + mrSges(6,2) * t198;
t218 = m(6) * pkin(4) - t232;
t329 = t218 * t190;
t285 = pkin(7) + t175;
t159 = t285 * t199;
t160 = t285 * t202;
t328 = -t292 * t159 - t160 * t207;
t182 = t202 * qJDD(2);
t141 = -t161 * t199 + t182;
t142 = t199 * qJDD(2) + t202 * t161;
t326 = -t141 * t199 + t142 * t202;
t94 = mrSges(6,2) * t153 + mrSges(6,3) * t239;
t95 = -mrSges(6,1) * t153 - mrSges(6,3) * t137;
t325 = -t198 * t95 + t201 * t94;
t60 = -mrSges(6,2) * t119 + mrSges(6,3) * t100;
t61 = mrSges(6,1) * t119 - mrSges(6,3) * t101;
t324 = -t198 * t61 + t201 * t60;
t126 = t182 + (-pkin(7) * qJDD(1) - t161) * t199;
t253 = qJDD(1) * t202;
t127 = pkin(7) * t253 + t142;
t241 = qJD(4) * t292;
t247 = t207 * t126 + t292 * t127 + t139 * t241;
t35 = qJDD(4) * qJ(5) + (qJD(5) - t128) * qJD(4) + t247;
t149 = qJDD(1) * t169 + qJDD(3);
t50 = pkin(4) * t119 - qJ(5) * t118 - qJD(5) * t154 + t149;
t14 = -t198 * t35 + t201 * t50;
t15 = t198 * t50 + t201 * t35;
t226 = -t14 * t198 + t15 * t201;
t234 = mrSges(5,1) * t190 - mrSges(5,2) * t187;
t235 = -mrSges(4,1) * t202 + mrSges(4,2) * t199;
t321 = m(4) * pkin(2) + t187 * mrSges(7,3) + mrSges(3,1) + t234 - t235;
t320 = -m(6) * t70 + t338;
t11 = pkin(8) * t100 + t15;
t7 = pkin(5) * t119 - pkin(8) * t101 + t14;
t1 = qJD(6) * t8 + t11 * t209 + t206 * t7;
t2 = -qJD(6) * t9 - t11 * t206 + t209 * t7;
t319 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t205 = -pkin(7) - qJ(3);
t318 = -m(7) * pkin(5) * t198 - m(4) * qJ(3) + m(6) * t205 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - t231;
t317 = Ifges(7,4) * t315 + Ifges(7,2) * t314 + Ifges(7,6) * t303;
t316 = Ifges(7,1) * t315 + Ifges(7,4) * t314 + Ifges(7,5) * t303;
t311 = Ifges(6,1) * t304 + Ifges(6,4) * t305 + Ifges(6,5) * t302;
t310 = -t351 / 0.2e1;
t309 = t351 / 0.2e1;
t308 = -t77 / 0.2e1;
t307 = t77 / 0.2e1;
t306 = m(6) + m(7);
t301 = -t150 / 0.2e1;
t300 = t150 / 0.2e1;
t299 = t153 / 0.2e1;
t298 = -t153 / 0.2e1;
t296 = t154 / 0.2e1;
t293 = t201 / 0.2e1;
t208 = sin(qJ(1));
t289 = pkin(1) * t208;
t210 = cos(qJ(1));
t192 = t210 * pkin(1);
t81 = t220 * qJD(3) + qJD(4) * t328;
t90 = pkin(4) * t158 - qJ(5) * t157 - qJD(5) * t166;
t47 = t198 * t90 + t201 * t81;
t283 = mrSges(5,3) * t153;
t281 = Ifges(5,4) * t154;
t280 = Ifges(6,4) * t198;
t279 = Ifges(6,4) * t201;
t103 = -pkin(4) * t220 - qJ(5) * t166 + t169;
t113 = -t207 * t159 + t160 * t292;
t59 = t198 * t103 + t201 * t113;
t265 = t157 * t198;
t264 = t157 * t201;
t262 = t166 * t198;
t261 = t166 * t201;
t260 = t188 * t190;
t259 = t190 * t191;
t257 = t191 * t179 + t192;
t256 = qJD(4) * t207;
t44 = -mrSges(7,1) * t351 + mrSges(7,2) * t77;
t252 = -t44 + t338;
t251 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t117;
t250 = m(4) + m(5) + t306;
t249 = -t198 * (t137 * Ifges(6,4) + Ifges(6,2) * t239 - Ifges(6,6) * t153) / 0.2e1;
t248 = (t137 * Ifges(6,1) + Ifges(6,4) * t239 - Ifges(6,5) * t153) * t293;
t180 = -pkin(2) - t290;
t243 = m(6) * qJ(5) + mrSges(6,3);
t242 = m(7) * t284 + mrSges(7,3);
t10 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t46 = -t198 * t81 + t201 * t90;
t240 = t119 * mrSges(5,1) + t118 * mrSges(5,2);
t58 = t201 * t103 - t113 * t198;
t236 = -mrSges(4,1) * t253 + mrSges(4,2) * t254;
t230 = Ifges(6,1) * t201 - t280;
t229 = -Ifges(6,2) * t198 + t279;
t228 = Ifges(6,5) * t201 - Ifges(6,6) * t198;
t225 = t198 * t42 - t201 * t43;
t45 = -pkin(5) * t220 - pkin(8) * t261 + t58;
t51 = -pkin(8) * t262 + t59;
t16 = -t206 * t51 + t209 * t45;
t17 = t206 * t45 + t209 * t51;
t224 = -(-t172 * t199 + t185) * t199 + t146 * t202;
t178 = pkin(5) * t201 + pkin(4);
t223 = t178 * t190 + t187 * t284;
t143 = -qJD(4) * mrSges(5,2) + t283;
t221 = t143 + t325;
t40 = t126 * t292 - t207 * t127 - t139 * t256 - t140 * t241;
t195 = pkin(11) + qJ(6);
t186 = sin(t195);
t189 = cos(t195);
t216 = m(7) * t178 + mrSges(7,1) * t189 - mrSges(7,2) * t186;
t36 = -qJDD(4) * pkin(4) + qJDD(5) - t40;
t215 = t187 * t243 + t329;
t82 = qJD(3) * t166 + qJD(4) * t113;
t168 = qJDD(1) * t180 + qJDD(3);
t148 = Ifges(5,4) * t153;
t134 = t186 * t188 + t189 * t259;
t133 = -t186 * t259 + t188 * t189;
t132 = t186 * t191 - t189 * t260;
t131 = t186 * t260 + t189 * t191;
t111 = t154 * Ifges(5,1) + t148 + t334;
t110 = t153 * Ifges(5,2) + t281 + t333;
t109 = t222 * t166;
t108 = t165 * t166;
t105 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t119;
t85 = pkin(5) * t262 - t328;
t65 = pkin(5) * t265 + t82;
t64 = pkin(5) * t267 + t80;
t63 = mrSges(7,1) * t150 - mrSges(7,3) * t77;
t62 = -mrSges(7,2) * t150 + mrSges(7,3) * t351;
t57 = t155 * t166 - t157 * t165;
t56 = -t156 * t166 - t157 * t222;
t55 = -pkin(5) * t239 + t70;
t39 = -t140 * t256 + t247;
t37 = t101 * Ifges(6,4) + t100 * Ifges(6,2) + t119 * Ifges(6,6);
t34 = -pkin(8) * t265 + t47;
t27 = pkin(5) * t158 - pkin(8) * t264 + t46;
t20 = -t100 * pkin(5) + t36;
t19 = -mrSges(7,2) * t117 + mrSges(7,3) * t26;
t18 = mrSges(7,1) * t117 - mrSges(7,3) * t25;
t4 = -qJD(6) * t17 - t206 * t34 + t209 * t27;
t3 = qJD(6) * t16 + t206 * t27 + t209 * t34;
t5 = [(-m(5) * t79 - t320) * t82 + (mrSges(2,1) * t208 - t132 * mrSges(7,1) + mrSges(2,2) * t210 - t131 * mrSges(7,2) + t349 * (-t191 * t205 - t289) + (m(6) - t350) * t289 + t318 * t191 + (m(5) * t179 - m(7) * (-t179 - t223) - m(6) * (-qJ(5) * t187 - t179) + t187 * mrSges(6,3) + t329 + t321) * t188) * g(1) + (-m(6) * t257 - mrSges(2,1) * t210 - t134 * mrSges(7,1) + mrSges(2,2) * t208 - t133 * mrSges(7,2) + t349 * (-t188 * t205 + t257) + t350 * t192 + t318 * t188 + (-m(7) * t223 - t215 - t321) * t191) * g(2) + m(6) * (t14 * t58 + t15 * t59 + t42 * t46 + t43 * t47) + (-t1 * t108 + t109 * t2 - t56 * t8 + t57 * t9) * mrSges(7,3) + (-Ifges(7,1) * t109 - Ifges(7,4) * t108) * t315 + (-Ifges(7,5) * t109 - Ifges(7,6) * t108) * t303 + (-Ifges(7,4) * t109 - Ifges(7,2) * t108) * t314 + t20 * (mrSges(7,1) * t108 - mrSges(7,2) * t109) + (-t14 * t261 - t15 * t262 - t264 * t42 - t265 * t43) * mrSges(6,3) + (Ifges(7,1) * t56 + Ifges(7,4) * t57) * t307 + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t203 * mrSges(3,1) - 0.2e1 * t200 * mrSges(3,2) + m(3) * (t200 ^ 2 + t203 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (t149 * mrSges(5,2) - t40 * mrSges(5,3) + Ifges(5,1) * t118 - Ifges(5,4) * t119 + Ifges(5,5) * qJDD(4) + t228 * t302 + t229 * t305 + t230 * t304 + t231 * t36) * t166 + (Ifges(4,4) * t199 + Ifges(4,2) * t202) * t253 + (Ifges(4,1) * t199 + Ifges(4,4) * t202) * t254 - (-m(5) * t40 + m(6) * t36 + t335) * t328 - (t149 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2) - t39 * mrSges(5,3) - Ifges(5,4) * t118 + Ifges(6,5) * t304 + Ifges(7,5) * t315 + Ifges(5,2) * t119 - Ifges(5,6) * qJDD(4) + Ifges(6,6) * t305 + Ifges(7,6) * t314 + Ifges(6,3) * t302 + Ifges(7,3) * t303 + t319) * t220 - (Ifges(6,5) * t101 + Ifges(6,6) * t100 + Ifges(6,3) * t119 + t251) * t220 / 0.2e1 + m(4) * (t224 * qJD(3) + t168 * t180 + t175 * t326) + (Ifges(7,5) * t56 + Ifges(7,6) * t57) * t300 + (Ifges(7,4) * t56 + Ifges(7,2) * t57) * t309 + m(7) * (t1 * t17 + t16 * t2 + t20 * t85 + t3 * t9 + t4 * t8 + t55 * t65) + m(5) * (t113 * t39 + t149 * t169 + t80 * t81) + t81 * t143 + t113 * t105 + t47 * t94 + t46 * t95 + t85 * t10 + t261 * t311 + t56 * t312 + t57 * t313 - t109 * t316 + t168 * t235 + t180 * t236 - t108 * t317 + t169 * t240 + (t161 * t353 + t326) * mrSges(4,3) - t37 * t262 / 0.2e1 + (-t79 * mrSges(5,3) + t339 + t248 + t249 + t228 * t298 + Ifges(5,4) * t299 + Ifges(5,1) * t296 + t152 * mrSges(5,2) + t334 / 0.2e1 + t111 / 0.2e1 + t239 * t229 / 0.2e1 + t137 * t230 / 0.2e1) * t157 + (Ifges(6,3) * t298 - Ifges(5,2) * t299 + Ifges(7,3) * t300 + Ifges(7,5) * t307 + Ifges(7,6) * t309 - Ifges(5,4) * t296 + t348 - t347 + t152 * mrSges(5,1) - t333 / 0.2e1 - t110 / 0.2e1 + t340 / 0.2e1 + t341 / 0.2e1 + t346 - t345 + t344 / 0.2e1 - t80 * mrSges(5,3)) * t158 + t16 * t18 + t17 * t19 + t55 * (-mrSges(7,1) * t57 + mrSges(7,2) * t56) + t59 * t60 + t58 * t61 + t3 * t62 + t4 * t63 + t65 * t44; m(3) * qJDD(2) - t108 * t18 - t109 * t19 + t56 * t62 + t57 * t63 + (t105 + t324) * t166 - (t10 + t335) * t220 - t252 * t158 + t221 * t157 + (-m(3) - t250) * g(3) + m(7) * (-t1 * t109 - t108 * t2 + t158 * t55 - t20 * t220 + t56 * t9 + t57 * t8) + m(5) * (t157 * t80 - t158 * t79 + t166 * t39 + t220 * t40) + m(6) * (-t157 * t225 + t158 * t70 + t166 * t226 - t220 * t36) + m(4) * (t141 * t202 + t142 * t199); t240 + t252 * t154 + t236 - t336 * t63 - t337 * t62 - t221 * t153 + t201 * t61 + t198 * t60 - t222 * t18 + t165 * t19 - t353 * qJD(1) ^ 2 * mrSges(4,3) + (-g(1) * t188 + g(2) * t191) * t250 + (t1 * t165 - t154 * t55 - t2 * t222 - t336 * t8 - t337 * t9) * m(7) + (t14 * t201 + t15 * t198 + t153 * t225 - t154 * t70) * m(6) + (-t153 * t80 + t154 * t79 + t149) * m(5) + (-qJD(1) * t224 + t168) * m(4); (t282 + t320) * t80 - t154 * t348 - t154 * t346 + t342 * t63 + (t1 * t130 + t129 * t2 - t178 * t20 + t342 * t8 + t343 * t9 - t55 * t64) * m(7) + t343 * t62 + (t352 * (mrSges(5,1) + t216 + t218) - t242 * g(3)) * t187 + (-Ifges(7,5) * t155 - Ifges(7,6) * t156) * t300 + (-Ifges(7,1) * t155 - Ifges(7,4) * t156) * t307 + (-Ifges(7,4) * t155 - Ifges(7,2) * t156) * t309 - (Ifges(5,1) * t153 - t281 + t344) * t154 / 0.2e1 - t239 * (Ifges(6,6) * t154 + t153 * t229) / 0.2e1 + (Ifges(6,3) * t154 + t153 * t228) * t299 - qJD(4) * (Ifges(5,5) * t153 - Ifges(5,6) * t154) / 0.2e1 - t152 * (mrSges(5,1) * t154 + mrSges(5,2) * t153) - t137 * (Ifges(6,5) * t154 + t153 * t230) / 0.2e1 + (t266 * t42 + t267 * t43 + t226) * mrSges(6,3) + (-t143 + t283) * t79 + (t352 * (mrSges(5,2) - t242 - t243) - t216 * g(3)) * t190 - t336 * t313 - t337 * t312 + (mrSges(7,1) * t336 - mrSges(7,2) * t337) * t55 + t324 * qJ(5) - t153 * t248 - t153 * t249 + t154 * t347 - t153 * t339 + t154 * t345 + (-t215 - t234) * g(3) + (-Ifges(7,5) * t97 - Ifges(7,6) * t96 + Ifges(7,3) * t154) * t301 + (-Ifges(7,1) * t97 - Ifges(7,4) * t96 + Ifges(7,5) * t154) * t308 + (-Ifges(7,4) * t97 - Ifges(7,2) * t96 + Ifges(7,6) * t154) * t310 + (-Ifges(5,2) * t154 + t111 + t148) * t298 + (Ifges(7,5) * t165 - Ifges(7,6) * t222) * t303 + t20 * (mrSges(7,1) * t222 + mrSges(7,2) * t165) + (Ifges(7,4) * t165 - Ifges(7,2) * t222) * t314 + (Ifges(7,1) * t165 - Ifges(7,4) * t222) * t315 + (-t1 * t222 - t165 * t2 - t336 * t9 + t337 * t8) * mrSges(7,3) - t222 * t317 + t325 * qJD(5) + (-pkin(4) * t36 + qJ(5) * t226 - qJD(5) * t225 - t42 * t52 - t43 * t53) * m(6) + (Ifges(6,5) * t198 + Ifges(6,6) * t201) * t302 + (Ifges(6,1) * t198 + t279) * t304 + (Ifges(6,2) * t201 + t280) * t305 + t37 * t293 + t110 * t296 - t178 * t10 + t129 * t18 + t130 * t19 + Ifges(5,5) * t118 - Ifges(5,6) * t119 - t53 * t94 - t52 * t95 + t198 * t311 + t165 * t316 + t36 * t232 + Ifges(5,3) * qJDD(4) - t39 * mrSges(5,2) + t40 * mrSges(5,1) - pkin(4) * t54 - t64 * t44; t306 * t190 * g(3) + t137 * t95 - t239 * t94 - t351 * t62 + t77 * t63 + t10 + t54 + (-t351 * t9 + t77 * t8 + t20 - t332) * m(7) + (t137 * t42 - t239 * t43 - t332 + t36) * m(6); -t55 * (mrSges(7,1) * t77 + mrSges(7,2) * t351) + (Ifges(7,1) * t351 - t291) * t308 + t31 * t307 + (Ifges(7,5) * t351 - Ifges(7,6) * t77) * t301 - t8 * t62 + t9 * t63 - g(1) * (mrSges(7,1) * t133 - mrSges(7,2) * t134) - g(2) * (-mrSges(7,1) * t131 + mrSges(7,2) * t132) - g(3) * (-mrSges(7,1) * t186 - mrSges(7,2) * t189) * t187 + (t351 * t8 + t77 * t9) * mrSges(7,3) + t251 + (-Ifges(7,2) * t77 + t32 + t73) * t310 + t319;];
tau  = t5;
