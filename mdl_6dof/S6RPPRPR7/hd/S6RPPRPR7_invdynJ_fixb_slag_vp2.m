% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:54
% EndTime: 2019-03-09 01:53:15
% DurationCPUTime: 14.52s
% Computational Cost: add. (7563->581), mult. (16025->758), div. (0->0), fcn. (11516->14), ass. (0->256)
t305 = pkin(8) + qJ(5);
t383 = -m(6) * qJ(5) - m(7) * t305 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t204 = pkin(9) + qJ(4);
t194 = sin(t204);
t196 = cos(t204);
t382 = mrSges(5,1) * t194 + t196 * t383;
t209 = sin(pkin(9));
t211 = cos(pkin(9));
t271 = t209 ^ 2 + t211 ^ 2;
t250 = t271 * mrSges(4,3);
t208 = sin(pkin(10));
t210 = cos(pkin(10));
t215 = sin(qJ(6));
t218 = cos(qJ(6));
t163 = t208 * t218 + t210 * t215;
t155 = t163 * qJD(6);
t216 = sin(qJ(4));
t311 = cos(qJ(4));
t226 = -t209 * t311 - t211 * t216;
t151 = t226 * qJD(1);
t96 = t163 * t151;
t381 = t96 - t155;
t230 = t208 * t215 - t210 * t218;
t154 = t230 * qJD(6);
t97 = t230 * t151;
t380 = t97 - t154;
t217 = sin(qJ(1));
t219 = cos(qJ(1));
t379 = g(1) * t217 - g(2) * t219;
t214 = -pkin(1) - qJ(3);
t378 = -qJD(1) * qJD(3) + qJDD(1) * t214;
t186 = pkin(5) * t210 + pkin(4);
t237 = -mrSges(6,1) * t210 + mrSges(6,2) * t208;
t377 = -m(6) * pkin(4) - m(7) * t186 + t237;
t149 = qJD(6) - t151;
t258 = t311 * t211;
t270 = qJD(1) * t209;
t153 = qJD(1) * t258 - t216 * t270;
t132 = qJD(4) * t208 + t153 * t210;
t243 = qJD(4) * t210 - t153 * t208;
t74 = t132 * t218 + t215 * t243;
t310 = Ifges(7,4) * t74;
t375 = -t132 * t215 + t218 * t243;
t29 = Ifges(7,2) * t375 + Ifges(7,6) * t149 + t310;
t331 = t29 / 0.2e1;
t71 = Ifges(7,4) * t375;
t30 = Ifges(7,1) * t74 + Ifges(7,5) * t149 + t71;
t330 = t30 / 0.2e1;
t225 = -t209 * t216 + t258;
t120 = qJD(4) * t151 + qJDD(1) * t225;
t100 = qJDD(4) * t210 - t120 * t208;
t324 = t100 / 0.2e1;
t101 = qJDD(4) * t208 + t120 * t210;
t323 = t101 / 0.2e1;
t121 = qJD(4) * t153 - qJDD(1) * t226;
t321 = t121 / 0.2e1;
t23 = qJD(6) * t375 + t100 * t215 + t101 * t218;
t333 = t23 / 0.2e1;
t24 = -qJD(6) * t74 + t100 * t218 - t101 * t215;
t332 = t24 / 0.2e1;
t373 = Ifges(6,1) * t323 + Ifges(6,4) * t324 + Ifges(6,5) * t321;
t372 = -m(7) - m(5);
t173 = qJD(1) * t214 + qJD(2);
t248 = -pkin(7) * qJD(1) + t173;
t146 = t248 * t209;
t147 = t248 * t211;
t99 = t146 * t311 + t147 * t216;
t89 = qJD(4) * qJ(5) + t99;
t190 = qJD(1) * qJ(2) + qJD(3);
t168 = pkin(3) * t270 + t190;
t90 = -pkin(4) * t151 - qJ(5) * t153 + t168;
t47 = -t208 * t89 + t210 * t90;
t27 = -pkin(5) * t151 - pkin(8) * t132 + t47;
t48 = t208 * t90 + t210 * t89;
t32 = pkin(8) * t243 + t48;
t9 = -t215 * t32 + t218 * t27;
t371 = t9 * mrSges(7,1);
t119 = qJDD(6) + t121;
t322 = t119 / 0.2e1;
t10 = t215 * t27 + t218 * t32;
t370 = t10 * mrSges(7,2);
t369 = t47 * mrSges(6,1);
t368 = t48 * mrSges(6,2);
t363 = t243 * Ifges(6,6);
t364 = t132 * Ifges(6,5);
t367 = t74 * Ifges(7,5) + Ifges(7,6) * t375 - t151 * Ifges(6,3) + t149 * Ifges(7,3) + t363 + t364;
t169 = t305 * t208;
t170 = t305 * t210;
t125 = -t169 * t218 - t170 * t215;
t284 = t151 * t210;
t116 = pkin(4) * t153 - qJ(5) * t151;
t143 = t216 * t146;
t98 = t147 * t311 - t143;
t53 = t116 * t210 - t208 * t98;
t35 = pkin(5) * t153 - pkin(8) * t284 + t53;
t285 = t151 * t208;
t54 = t116 * t208 + t210 * t98;
t42 = -pkin(8) * t285 + t54;
t366 = -qJD(5) * t230 + qJD(6) * t125 - t215 * t35 - t218 * t42;
t126 = -t169 * t215 + t170 * t218;
t365 = -qJD(5) * t163 - qJD(6) * t126 + t215 * t42 - t218 * t35;
t303 = mrSges(5,3) * t153;
t362 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t243 - mrSges(6,2) * t132 - t303;
t252 = qJD(4) * t311;
t269 = qJD(4) * t216;
t157 = -t209 * t269 + t211 * t252;
t361 = qJD(1) * t230 - t154 * t226 - t157 * t163;
t360 = -qJD(1) * t163 + t155 * t226 - t157 * t230;
t52 = -mrSges(6,1) * t100 + mrSges(6,2) * t101;
t357 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t120 + t52;
t356 = Ifges(5,5) * qJD(4);
t355 = Ifges(5,6) * qJD(4);
t354 = t379 * t196;
t306 = -pkin(7) + t214;
t166 = t306 * t209;
t167 = t306 * t211;
t351 = -t166 * t216 + t167 * t311;
t197 = t209 * pkin(3);
t213 = -pkin(7) - qJ(3);
t350 = t197 * t219 + t213 * t217;
t206 = qJD(1) * qJD(2);
t174 = -qJDD(1) * qJ(2) - t206;
t203 = pkin(10) + qJ(6);
t193 = sin(t203);
t195 = cos(t203);
t349 = mrSges(7,1) * t195 - mrSges(7,2) * t193 - t377;
t94 = mrSges(6,2) * t151 + mrSges(6,3) * t243;
t95 = -mrSges(6,1) * t151 - mrSges(6,3) * t132;
t347 = -t208 * t95 + t210 * t94;
t59 = -mrSges(6,2) * t121 + mrSges(6,3) * t100;
t60 = mrSges(6,1) * t121 - mrSges(6,3) * t101;
t346 = -t208 * t60 + t210 * t59;
t156 = -t209 * t252 - t211 * t269;
t165 = qJDD(2) + t378;
t244 = -pkin(7) * qJDD(1) + t165;
t137 = t244 * t209;
t138 = t244 * t211;
t50 = -t137 * t216 + t138 * t311 - t146 * t252 - t147 * t269;
t46 = -qJDD(4) * pkin(4) + qJDD(5) - t50;
t84 = -qJD(4) * pkin(4) + qJD(5) - t98;
t345 = t156 * t84 + t225 * t46;
t260 = t137 * t311 + t138 * t216 + t147 * t252;
t44 = qJDD(4) * qJ(5) + (qJD(5) - t143) * qJD(4) + t260;
t171 = qJDD(3) - t174;
t267 = qJDD(1) * t209;
t158 = pkin(3) * t267 + t171;
t45 = pkin(4) * t121 - qJ(5) * t120 - qJD(5) * t153 + t158;
t14 = -t208 * t44 + t210 * t45;
t15 = t208 * t45 + t210 * t44;
t232 = -t14 * t208 + t15 * t210;
t49 = -t146 * t269 + t260;
t341 = -t156 * t98 - t157 * t99 - t225 * t50 + t226 * t49;
t340 = -m(6) * t84 + t362;
t11 = pkin(8) * t100 + t15;
t8 = pkin(5) * t121 - pkin(8) * t101 + t14;
t1 = qJD(6) * t9 + t11 * t218 + t215 * t8;
t2 = -qJD(6) * t10 - t11 * t215 + t218 * t8;
t339 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t240 = mrSges(4,1) * t209 + mrSges(4,2) * t211;
t338 = t194 * t377 + mrSges(2,2) - mrSges(3,3) - t240 - t382;
t337 = m(4) * t190 + m(5) * t168 - mrSges(5,1) * t151 + mrSges(5,2) * t153 + qJD(1) * t240;
t236 = mrSges(6,1) * t208 + mrSges(6,2) * t210;
t336 = m(7) * pkin(5) * t208 + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + t236;
t335 = Ifges(7,4) * t333 + Ifges(7,2) * t332 + Ifges(7,6) * t322;
t334 = Ifges(7,1) * t333 + Ifges(7,4) * t332 + Ifges(7,5) * t322;
t329 = -t375 / 0.2e1;
t328 = t375 / 0.2e1;
t327 = -t74 / 0.2e1;
t326 = t74 / 0.2e1;
t325 = -m(6) - m(7);
t320 = -t149 / 0.2e1;
t319 = t149 / 0.2e1;
t318 = t151 / 0.2e1;
t317 = -t151 / 0.2e1;
t315 = t153 / 0.2e1;
t312 = t210 / 0.2e1;
t79 = pkin(4) * t157 - qJ(5) * t156 - qJD(5) * t225 + qJD(2);
t85 = qJD(3) * t226 + qJD(4) * t351;
t39 = t208 * t79 + t210 * t85;
t304 = mrSges(5,3) * t151;
t302 = Ifges(5,4) * t153;
t301 = Ifges(6,4) * t208;
t300 = Ifges(6,4) * t210;
t286 = qJDD(1) * pkin(1);
t283 = t156 * t208;
t282 = t156 * t210;
t281 = t225 * t208;
t280 = t225 * t210;
t278 = t193 * t217;
t277 = t193 * t219;
t276 = t195 * t217;
t275 = t195 * t219;
t182 = qJ(2) + t197;
t115 = -pkin(4) * t226 - qJ(5) * t225 + t182;
t124 = t166 * t311 + t167 * t216;
t65 = t115 * t208 + t124 * t210;
t266 = qJDD(1) * t211;
t273 = mrSges(4,1) * t267 + mrSges(4,2) * t266;
t272 = pkin(1) * t219 + qJ(2) * t217;
t36 = -mrSges(7,1) * t375 + mrSges(7,2) * t74;
t265 = -t36 + t362;
t264 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t119;
t263 = -m(4) - m(5) + t325;
t262 = -t208 * (Ifges(6,4) * t132 + Ifges(6,2) * t243 - Ifges(6,6) * t151) / 0.2e1;
t261 = (Ifges(6,1) * t132 + Ifges(6,4) * t243 - Ifges(6,5) * t151) * t312;
t259 = t197 * t217 + t272;
t7 = -mrSges(7,1) * t24 + mrSges(7,2) * t23;
t199 = t219 * qJ(2);
t251 = -pkin(1) * t217 + t199;
t38 = -t208 * t85 + t210 * t79;
t249 = t121 * mrSges(5,1) + mrSges(5,2) * t120;
t247 = t271 * t173;
t246 = t271 * t165;
t64 = t115 * t210 - t124 * t208;
t235 = Ifges(6,1) * t210 - t301;
t234 = -Ifges(6,2) * t208 + t300;
t233 = Ifges(6,5) * t210 - Ifges(6,6) * t208;
t231 = t208 * t47 - t210 * t48;
t43 = -pkin(5) * t226 - pkin(8) * t280 + t64;
t51 = -pkin(8) * t281 + t65;
t16 = -t215 * t51 + t218 * t43;
t17 = t215 * t43 + t218 * t51;
t144 = -qJD(4) * mrSges(5,2) + t304;
t227 = t144 + t347;
t109 = t163 * t225;
t86 = qJD(3) * t225 + qJD(4) * t124;
t220 = qJD(1) ^ 2;
t191 = qJDD(2) - t286;
t148 = Ifges(5,4) * t151;
t142 = t194 * t275 - t278;
t141 = t194 * t277 + t276;
t140 = t194 * t276 + t277;
t139 = -t194 * t278 + t275;
t113 = t153 * Ifges(5,1) + t148 + t356;
t112 = t151 * Ifges(5,2) + t302 + t355;
t111 = t230 * t225;
t110 = t230 * t226;
t108 = t163 * t226;
t103 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t121;
t87 = pkin(5) * t281 - t351;
t70 = pkin(5) * t285 + t99;
t69 = pkin(5) * t283 + t86;
t63 = -pkin(5) * t243 + t84;
t62 = mrSges(7,1) * t149 - mrSges(7,3) * t74;
t61 = -mrSges(7,2) * t149 + mrSges(7,3) * t375;
t58 = t154 * t225 - t156 * t163;
t56 = -qJD(6) * t109 - t156 * t230;
t33 = Ifges(6,4) * t101 + Ifges(6,2) * t100 + Ifges(6,6) * t121;
t31 = -pkin(8) * t283 + t39;
t26 = pkin(5) * t157 - pkin(8) * t282 + t38;
t25 = -pkin(5) * t100 + t46;
t19 = -mrSges(7,2) * t119 + mrSges(7,3) * t24;
t18 = mrSges(7,1) * t119 - mrSges(7,3) * t23;
t4 = -qJD(6) * t17 - t215 * t31 + t218 * t26;
t3 = qJD(6) * t16 + t215 * t26 + t218 * t31;
t5 = [(Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (t191 - t286) * mrSges(3,2) + m(5) * (t124 * t49 + t158 * t182 + t85 * t99) + m(6) * (t14 * t64 + t15 * t65 + t38 * t47 + t39 * t48) + (t168 * mrSges(5,1) - t370 + t371 - t355 / 0.2e1 - t112 / 0.2e1 + t369 - t368 + t363 / 0.2e1 + t364 / 0.2e1 - Ifges(5,4) * t315 + Ifges(6,3) * t317 - Ifges(5,2) * t318 + Ifges(7,3) * t319 + Ifges(7,5) * t326 + Ifges(7,6) * t328 + t367 / 0.2e1) * t157 + t337 * qJD(2) + (-t1 * t109 + t10 * t58 + t111 * t2 - t56 * t9) * mrSges(7,3) + (-Ifges(7,1) * t111 - Ifges(7,4) * t109) * t333 + t25 * (mrSges(7,1) * t109 - mrSges(7,2) * t111) + (-Ifges(7,5) * t111 - Ifges(7,6) * t109) * t322 + (-Ifges(7,4) * t111 - Ifges(7,2) * t109) * t332 + m(7) * (t1 * t17 + t10 * t3 + t16 * t2 + t25 * t87 + t4 * t9 + t63 * t69) + (Ifges(7,1) * t56 + Ifges(7,4) * t58) * t326 + m(4) * (qJ(2) * t171 - qJD(3) * t247 + t214 * t246) + (-m(6) * t259 - t140 * mrSges(7,1) - t139 * mrSges(7,2) + (-m(4) - m(3)) * t272 + t372 * (-t213 * t219 + t259) + (-m(4) * qJ(3) + m(6) * t213 - t336) * t219 + t338 * t217) * g(2) + (-t142 * mrSges(7,1) + t141 * mrSges(7,2) - m(4) * t199 - m(6) * (t199 + t350) - m(3) * t251 + t372 * (t251 + t350) + (-m(4) * t214 + m(6) * pkin(1) + t336) * t217 + t338 * t219) * g(1) + (-t165 - t378) * t250 + t280 * t373 + (t262 + t261 + t168 * mrSges(5,2) + t356 / 0.2e1 + t113 / 0.2e1 + t243 * t234 / 0.2e1 + t132 * t235 / 0.2e1 + Ifges(5,1) * t315 + t233 * t317 + Ifges(5,4) * t318) * t156 - 0.2e1 * t174 * mrSges(3,3) - (-m(5) * t50 + m(6) * t46 + t357) * t351 + t345 * t236 + (-m(5) * t98 - t340) * t86 + t341 * mrSges(5,3) - (-t158 * mrSges(5,2) - Ifges(5,1) * t120 + Ifges(5,4) * t121 - Ifges(5,5) * qJDD(4) - t233 * t321 - t234 * t324 - t235 * t323) * t225 - (Ifges(6,5) * t101 + Ifges(6,6) * t100 + Ifges(6,3) * t121 + t264) * t226 / 0.2e1 - (t158 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2) - Ifges(5,4) * t120 + Ifges(6,5) * t323 + Ifges(7,5) * t333 + Ifges(5,2) * t121 - Ifges(5,6) * qJDD(4) + Ifges(6,6) * t324 + Ifges(7,6) * t332 + Ifges(6,3) * t321 + Ifges(7,3) * t322 + t339) * t226 + (-t14 * t280 - t15 * t281 - t282 * t47 - t283 * t48) * mrSges(6,3) + t85 * t144 + t124 * t103 + t39 * t94 + t38 * t95 + t87 * t7 + t4 * t62 + t63 * (-mrSges(7,1) * t58 + mrSges(7,2) * t56) + t64 * t60 + t65 * t59 + t69 * t36 + t3 * t61 + t16 * t18 + t17 * t19 - t33 * t281 / 0.2e1 + t171 * t240 + (Ifges(7,5) * t56 + Ifges(7,6) * t58) * t319 + (Ifges(4,1) * t211 - Ifges(4,4) * t209) * t266 + t56 * t330 + t58 * t331 - t111 * t334 - t109 * t335 + t182 * t249 + (Ifges(7,4) * t56 + Ifges(7,2) * t58) * t328 - (Ifges(4,4) * t211 - Ifges(4,2) * t209) * t267 + m(3) * (-pkin(1) * t191 + (-t174 + t206) * qJ(2)) + qJ(2) * t273; t108 * t18 + t110 * t19 + t361 * t62 + t360 * t61 + (-m(3) * qJ(2) - mrSges(3,3)) * t220 - (t103 + t346) * t226 - (t7 + t357) * t225 + t227 * t157 + t265 * t156 + (mrSges(3,2) - t250) * qJDD(1) + (-t208 * t94 - t210 * t95 - t337) * qJD(1) - m(5) * t341 + m(3) * t191 + m(4) * t246 - t379 * (m(3) - t263) + (t1 * t110 + t10 * t360 + t108 * t2 - t156 * t63 - t225 * t25 + t361 * t9) * m(7) + (-t231 * t157 - t232 * t226 - (t208 * t48 + t210 * t47) * qJD(1) - t345) * m(6); t273 + t163 * t19 - t230 * t18 + t381 * t62 + t380 * t61 + t249 + t210 * t60 + t208 * t59 - t227 * t151 - t220 * t250 + t265 * t153 + (g(1) * t219 + g(2) * t217) * t263 + (t1 * t163 + t10 * t380 - t153 * t63 - t2 * t230 + t381 * t9) * m(7) + (t14 * t210 + t15 * t208 + t151 * t231 - t153 * t84) * m(6) + (-t151 * t99 + t153 * t98 + t158) * m(5) + (qJD(1) * t247 + t171) * m(4); (-mrSges(7,1) * t381 + mrSges(7,2) * t380) * t63 + (-t1 * t230 + t10 * t381 - t163 * t2 - t380 * t9) * mrSges(7,3) - t151 * t262 + t153 * t368 + t153 * t370 + (t194 * t349 + t382) * g(3) + (-Ifges(7,4) * t97 - Ifges(7,2) * t96 + Ifges(7,6) * t153) * t329 + (-Ifges(5,2) * t153 + t113 + t148) * t317 - t151 * t261 + (-Ifges(7,5) * t154 - Ifges(7,6) * t155) * t319 + (-Ifges(7,1) * t154 - Ifges(7,4) * t155) * t326 + (-Ifges(7,4) * t154 - Ifges(7,2) * t155) * t328 - t243 * (Ifges(6,6) * t153 + t151 * t234) / 0.2e1 - t168 * (mrSges(5,1) * t153 + mrSges(5,2) * t151) - qJD(4) * (Ifges(5,5) * t151 - Ifges(5,6) * t153) / 0.2e1 - t132 * (Ifges(6,5) * t153 + t151 * t235) / 0.2e1 + (Ifges(6,3) * t153 + t151 * t233) * t318 + (-pkin(4) * t46 + qJ(5) * t232 - qJD(5) * t231 - t47 * t53 - t48 * t54) * m(6) - t153 * t371 - (Ifges(5,1) * t151 - t302 + t367) * t153 / 0.2e1 - t153 * t369 + t365 * t62 + (t1 * t126 + t10 * t366 + t125 * t2 - t186 * t25 + t365 * t9 - t63 * t70) * m(7) + t366 * t61 + t208 * t373 + t379 * ((-mrSges(5,1) - t349) * t196 + t383 * t194) - t186 * t7 - t84 * t236 * t151 + (t284 * t47 + t285 * t48 + t232) * mrSges(6,3) + t346 * qJ(5) + t347 * qJD(5) + (t303 + t340) * t99 + (Ifges(7,1) * t163 - Ifges(7,4) * t230) * t333 + t25 * (mrSges(7,1) * t230 + mrSges(7,2) * t163) + (Ifges(7,5) * t163 - Ifges(7,6) * t230) * t322 + (Ifges(7,4) * t163 - Ifges(7,2) * t230) * t332 - t230 * t335 + (-Ifges(7,5) * t97 - Ifges(7,6) * t96 + Ifges(7,3) * t153) * t320 + (-Ifges(7,1) * t97 - Ifges(7,4) * t96 + Ifges(7,5) * t153) * t327 + (-t144 + t304) * t98 + t125 * t18 + t126 * t19 + Ifges(5,5) * t120 - Ifges(5,6) * t121 - t54 * t94 - t53 * t95 - t70 * t36 - pkin(4) * t52 - t49 * mrSges(5,2) + t50 * mrSges(5,1) + Ifges(5,3) * qJDD(4) + t46 * t237 + t380 * t330 + t381 * t331 + t33 * t312 + t112 * t315 + (Ifges(6,5) * t208 + Ifges(6,6) * t210) * t321 + (Ifges(6,1) * t208 + t300) * t323 + (Ifges(6,2) * t210 + t301) * t324 + t163 * t334; t325 * t194 * g(3) + t132 * t95 - t243 * t94 - t375 * t61 + t74 * t62 + t52 + t7 + (-t10 * t375 + t74 * t9 + t25 + t354) * m(7) + (t132 * t47 - t243 * t48 + t354 + t46) * m(6); -t63 * (mrSges(7,1) * t74 + mrSges(7,2) * t375) + (Ifges(7,1) * t375 - t310) * t327 + t29 * t326 + (Ifges(7,5) * t375 - Ifges(7,6) * t74) * t320 - t9 * t61 + t10 * t62 - g(1) * (mrSges(7,1) * t139 - mrSges(7,2) * t140) - g(2) * (mrSges(7,1) * t141 + mrSges(7,2) * t142) - g(3) * (-mrSges(7,1) * t193 - mrSges(7,2) * t195) * t196 + (t10 * t74 + t375 * t9) * mrSges(7,3) + t264 + (-Ifges(7,2) * t74 + t30 + t71) * t329 + t339;];
tau  = t5;
