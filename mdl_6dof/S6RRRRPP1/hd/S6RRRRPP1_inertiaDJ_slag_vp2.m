% Calculate time derivative of joint inertia matrix for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:44:05
% EndTime: 2019-03-09 20:44:17
% DurationCPUTime: 5.25s
% Computational Cost: add. (7662->468), mult. (17512->643), div. (0->0), fcn. (15994->8), ass. (0->217)
t229 = sin(qJ(3));
t230 = sin(qJ(2));
t232 = cos(qJ(3));
t233 = cos(qJ(2));
t198 = t229 * t233 + t232 * t230;
t335 = qJD(2) + qJD(3);
t161 = t335 * t198;
t350 = Ifges(7,4) + Ifges(6,5);
t355 = t161 * t350;
t354 = Ifges(6,1) + Ifges(7,1);
t353 = Ifges(6,4) - Ifges(7,5);
t349 = Ifges(6,6) - Ifges(7,6);
t228 = sin(qJ(4));
t278 = qJD(4) * t228;
t269 = t198 * t278;
t197 = t229 * t230 - t232 * t233;
t160 = t335 * t197;
t231 = cos(qJ(4));
t294 = t160 * t231;
t247 = t269 + t294;
t277 = qJD(4) * t231;
t268 = t198 * t277;
t295 = t160 * t228;
t248 = t268 - t295;
t61 = mrSges(5,1) * t248 - mrSges(5,2) * t247;
t321 = -pkin(8) - pkin(7);
t210 = t321 * t230;
t211 = t321 * t233;
t171 = t210 * t229 - t211 * t232;
t271 = qJD(2) * t321;
t204 = t230 * t271;
t256 = t233 * t271;
t97 = qJD(3) * t171 + t204 * t229 - t232 * t256;
t352 = -m(5) * t97 - t61;
t336 = (t228 ^ 2 + t231 ^ 2) * t232;
t340 = (mrSges(6,3) + mrSges(7,2));
t351 = 2 * t340;
t227 = sin(pkin(10));
t296 = cos(pkin(10));
t262 = t296 * t228;
t195 = t227 * t231 + t262;
t261 = t296 * t231;
t279 = qJD(4) * t198;
t54 = t160 * t195 + t227 * t269 - t261 * t279;
t286 = t227 * t228;
t246 = t261 - t286;
t55 = -t160 * t246 - t195 * t279;
t348 = t353 * t54 + t354 * t55 + t355;
t125 = t195 * t198;
t126 = t246 * t198;
t347 = -t353 * t125 + t126 * t354 + t350 * t197;
t188 = t246 * qJD(4);
t346 = mrSges(6,3) * t188;
t187 = t195 * qJD(4);
t345 = -t353 * t187 + t188 * t354;
t344 = t195 * t354 + t353 * t246;
t219 = -pkin(2) * t233 - pkin(1);
t144 = t197 * pkin(3) - t198 * pkin(9) + t219;
t84 = pkin(2) * qJD(2) * t230 + pkin(3) * t161 + pkin(9) * t160;
t338 = t232 * t210 + t211 * t229;
t96 = qJD(3) * t338 + t232 * t204 + t229 * t256;
t275 = t144 * t277 + t228 * t84 + t231 * t96;
t27 = -t171 * t278 + t275;
t264 = -t228 * t96 + t231 * t84;
t164 = t231 * t171;
t90 = t228 * t144 + t164;
t28 = -t90 * qJD(4) + t264;
t343 = -t228 * t28 + t231 * t27;
t342 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t341 = Ifges(5,5) * t277 - Ifges(5,6) * t278 - t187 * t349 + t188 * t350;
t339 = m(7) * qJD(6);
t337 = -mrSges(5,1) * t277 + mrSges(5,2) * t278;
t311 = pkin(4) * t227;
t212 = qJ(6) + t311;
t334 = m(7) * t212 + mrSges(7,3);
t270 = t296 * pkin(4);
t214 = -t270 - pkin(5);
t323 = m(6) * pkin(4);
t333 = m(7) * t214 - t296 * t323 - mrSges(6,1) - mrSges(7,1);
t332 = t227 * t323 - mrSges(6,2) + t334;
t331 = 0.2e1 * m(5);
t330 = 0.2e1 * m(6);
t329 = 0.2e1 * m(7);
t328 = 0.2e1 * t97;
t132 = t187 * mrSges(7,1) - t188 * mrSges(7,3);
t327 = 0.2e1 * t132;
t133 = t187 * mrSges(6,1) + t188 * mrSges(6,2);
t326 = 0.2e1 * t133;
t153 = -mrSges(7,1) * t246 - mrSges(7,3) * t195;
t325 = 0.2e1 * t153;
t324 = 0.2e1 * t219;
t318 = t161 / 0.2e1;
t249 = qJ(5) * t160 - qJD(5) * t198;
t10 = pkin(4) * t161 + t249 * t231 + (-t164 + (qJ(5) * t198 - t144) * t228) * qJD(4) + t264;
t16 = -qJ(5) * t268 + (-qJD(4) * t171 + t249) * t228 + t275;
t6 = t227 * t10 + t296 * t16;
t313 = pkin(2) * t232;
t154 = -mrSges(6,1) * t246 + mrSges(6,2) * t195;
t312 = pkin(4) * t154;
t310 = -qJ(5) - pkin(9);
t33 = -mrSges(6,2) * t161 + mrSges(6,3) * t54;
t36 = mrSges(7,2) * t54 + mrSges(7,3) * t161;
t309 = t33 + t36;
t34 = mrSges(6,1) * t161 - mrSges(6,3) * t55;
t35 = -t161 * mrSges(7,1) + t55 * mrSges(7,2);
t308 = t35 - t34;
t224 = t231 * qJ(5);
t89 = t231 * t144 - t171 * t228;
t64 = pkin(4) * t197 - t198 * t224 + t89;
t292 = t198 * t228;
t73 = -qJ(5) * t292 + t90;
t31 = t227 * t64 + t296 * t73;
t307 = mrSges(5,1) * t228;
t306 = mrSges(6,3) * t187;
t305 = Ifges(5,4) * t228;
t304 = Ifges(5,4) * t231;
t303 = pkin(2) * qJD(3);
t302 = t338 * t97;
t179 = t188 * mrSges(7,2);
t300 = t229 * mrSges(4,1);
t298 = t232 * mrSges(4,2);
t102 = -mrSges(7,2) * t125 + mrSges(7,3) * t197;
t99 = -mrSges(6,2) * t197 - mrSges(6,3) * t125;
t297 = t99 + t102;
t291 = t198 * t231;
t221 = pkin(4) * t278;
t222 = t229 * t303;
t203 = t222 + t221;
t290 = t203 * t154;
t206 = -mrSges(5,1) * t231 + mrSges(5,2) * t228;
t284 = t229 * t206;
t216 = pkin(2) * t229 + pkin(9);
t283 = -qJ(5) - t216;
t100 = mrSges(6,1) * t197 - mrSges(6,3) * t126;
t101 = -mrSges(7,1) * t197 + mrSges(7,2) * t126;
t282 = t101 - t100;
t276 = 0.2e1 * t233;
t274 = t232 * t303;
t218 = -pkin(4) * t231 - pkin(3);
t22 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t21 = -t54 * mrSges(7,1) - t55 * mrSges(7,3);
t267 = -t278 / 0.2e1;
t266 = -Ifges(5,6) * t228 - (2 * Ifges(4,4));
t192 = t216 * t231 + t224;
t259 = t283 * t228;
t130 = t192 * t227 - t259 * t296;
t131 = t192 * t296 + t227 * t259;
t223 = t231 * qJD(5);
t257 = t231 * t274;
t258 = qJD(4) * t283;
t159 = t228 * t258 + t223 + t257;
t240 = (-qJD(5) - t274) * t228 + t231 * t258;
t87 = t159 * t227 - t240 * t296;
t88 = t159 * t296 + t227 * t240;
t265 = t130 * t87 + t131 * t88;
t263 = qJD(4) * t310;
t185 = t228 * t263 + t223;
t242 = -qJD(5) * t228 + t231 * t263;
t120 = t185 * t227 - t242 * t296;
t121 = t185 * t296 + t227 * t242;
t207 = pkin(9) * t231 + t224;
t167 = t207 * t227 - t262 * t310;
t168 = t207 * t296 + t286 * t310;
t260 = t120 * t167 + t168 * t121;
t255 = mrSges(5,3) * t336;
t119 = pkin(4) * t292 - t338;
t253 = mrSges(5,2) * t231 + t307;
t252 = Ifges(5,1) * t231 - t305;
t251 = -Ifges(5,2) * t228 + t304;
t250 = Ifges(5,5) * t228 + Ifges(5,6) * t231;
t5 = t10 * t296 - t227 * t16;
t30 = -t227 * t73 + t296 * t64;
t245 = t120 * t130 + t121 * t131 + t167 * t87 + t168 * t88;
t244 = t132 + t133;
t243 = -Ifges(5,5) * t294 + t161 * t342 + t349 * t54 + t350 * t55;
t140 = -pkin(5) * t246 - qJ(6) * t195 + t218;
t98 = pkin(5) * t187 - qJ(6) * t188 - qJD(6) * t195 + t221;
t134 = Ifges(7,5) * t188 + Ifges(7,3) * t187;
t135 = Ifges(6,4) * t188 - Ifges(6,2) * t187;
t155 = Ifges(7,5) * t195 - Ifges(7,3) * t246;
t156 = Ifges(6,4) * t195 + Ifges(6,2) * t246;
t201 = t251 * qJD(4);
t202 = t252 * qJD(4);
t209 = Ifges(5,1) * t228 + t304;
t241 = t231 * t201 + t228 * t202 + t209 * t277 + t345 * t195 + (-t134 + t135) * t246 + t344 * t188 + (t155 - t156) * t187;
t60 = pkin(4) * t248 + t97;
t236 = -t270 * t346 + t214 * t179 - t306 * t311 + (qJD(6) * t246 - t187 * t212) * mrSges(7,2) + t341;
t146 = -mrSges(5,2) * t197 - mrSges(5,3) * t292;
t147 = mrSges(5,1) * t197 - mrSges(5,3) * t291;
t71 = mrSges(5,1) * t161 + mrSges(5,3) * t247;
t72 = -mrSges(5,2) * t161 - mrSges(5,3) * t248;
t235 = -t147 * t277 - t146 * t278 + m(5) * (-t277 * t89 - t278 * t90 + t343) + t231 * t72 - t228 * t71;
t1 = qJ(6) * t161 + qJD(6) * t197 + t6;
t111 = Ifges(5,6) * t197 + t198 * t251;
t112 = Ifges(5,5) * t197 + t198 * t252;
t17 = Ifges(7,5) * t55 + Ifges(7,6) * t161 - Ifges(7,3) * t54;
t18 = Ifges(6,4) * t55 + Ifges(6,2) * t54 + Ifges(6,6) * t161;
t200 = t253 * qJD(4);
t208 = Ifges(5,2) * t231 + t305;
t26 = qJ(6) * t197 + t31;
t29 = -t197 * pkin(5) - t30;
t3 = -t161 * pkin(5) - t5;
t44 = -Ifges(5,4) * t247 - Ifges(5,2) * t248 + Ifges(5,6) * t161;
t45 = -Ifges(5,1) * t247 - Ifges(5,4) * t248 + Ifges(5,5) * t161;
t46 = pkin(5) * t125 - qJ(6) * t126 + t119;
t66 = Ifges(7,5) * t126 + Ifges(7,6) * t197 + Ifges(7,3) * t125;
t67 = Ifges(6,4) * t126 - Ifges(6,2) * t125 + Ifges(6,6) * t197;
t8 = -pkin(5) * t54 - qJ(6) * t55 - qJD(6) * t126 + t60;
t234 = (-mrSges(4,1) + t206) * t97 + t250 * t318 + (-t155 / 0.2e1 + t156 / 0.2e1) * t54 - t338 * t200 + t29 * t179 - t31 * t306 + t111 * t267 + (t348 / 0.2e1 + t3 * mrSges(7,2) - t5 * mrSges(6,3) + t318 * t350) * t195 + (-t26 * mrSges(7,2) + t66 / 0.2e1 - t67 / 0.2e1) * t187 + (t134 / 0.2e1 - t135 / 0.2e1) * t125 + (t198 * t267 - t294 / 0.2e1) * t209 + t112 * t277 / 0.2e1 + t202 * t291 / 0.2e1 - t201 * t292 / 0.2e1 - t96 * mrSges(4,2) + t341 * t197 / 0.2e1 + t46 * t132 + t119 * t133 + ((-t228 * t90 - t231 * t89) * qJD(4) + t343) * mrSges(5,3) - t248 * t208 / 0.2e1 + t344 * t55 / 0.2e1 + t345 * t126 / 0.2e1 + t8 * t153 + t60 * t154 - Ifges(4,5) * t160 - Ifges(4,6) * t161 - t30 * t346 + t347 * t188 / 0.2e1 + (t1 * mrSges(7,2) + t6 * mrSges(6,3) - t17 / 0.2e1 + t18 / 0.2e1 + t349 * t318) * t246 + t228 * t45 / 0.2e1 + t231 * t44 / 0.2e1;
t217 = -pkin(3) - t313;
t205 = t218 - t313;
t141 = t253 * t198;
t123 = t140 - t313;
t91 = t222 + t98;
t75 = mrSges(6,1) * t125 + mrSges(6,2) * t126;
t74 = mrSges(7,1) * t125 - mrSges(7,3) * t126;
t2 = [t141 * t328 + (t1 * t26 + t29 * t3 + t46 * t8) * t329 + (t119 * t60 + t30 * t5 + t31 * t6) * t330 + 0.2e1 * (t160 * t338 - t161 * t171) * mrSges(4,3) - 0.2e1 * t338 * t61 - t112 * t294 + (mrSges(4,1) * t161 - mrSges(4,2) * t160) * t324 + (mrSges(4,3) * t328 - 0.2e1 * Ifges(4,1) * t160 - t228 * t44 + t231 * t45 + (Ifges(5,5) * t231 + t266) * t161 + (-t231 * t111 - t228 * t112 - t197 * t250) * qJD(4)) * t198 + (t27 * t90 + t28 * t89 - t302) * t331 + 0.2e1 * m(4) * (t171 * t96 - t302) + t111 * t295 + 0.2e1 * t31 * t33 + 0.2e1 * t30 * t34 + 0.2e1 * t29 * t35 + 0.2e1 * t26 * t36 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t233) * t276 + (0.2e1 * pkin(2) * (mrSges(4,1) * t197 + mrSges(4,2) * t198) + m(4) * pkin(2) * t324 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t230 + (Ifges(3,1) - Ifges(3,2)) * t276) * t230) * qJD(2) + (-t66 + t67) * t54 + 0.2e1 * t46 * t21 + 0.2e1 * t8 * t74 + 0.2e1 * t60 * t75 + 0.2e1 * t89 * t71 + 0.2e1 * t90 * t72 + 0.2e1 * t6 * t99 + 0.2e1 * t5 * t100 + 0.2e1 * t3 * t101 + 0.2e1 * t1 * t102 + 0.2e1 * t119 * t22 + (-0.2e1 * t96 * mrSges(4,3) - t266 * t160 + ((2 * Ifges(4,2)) + t342) * t161 + t243) * t197 + 0.2e1 * t27 * t146 + 0.2e1 * t28 * t147 + t347 * t55 + (-t161 * t349 + t17 - t18) * t125 + (t348 + t355) * t126; (Ifges(3,5) * t233 - Ifges(3,6) * t230 + (-mrSges(3,1) * t233 + mrSges(3,2) * t230) * pkin(7)) * qJD(2) + t297 * t88 + t282 * t87 + t309 * t131 + t308 * t130 + t234 + (m(4) * (t229 * t96 - t232 * t97) + (t160 * t232 - t161 * t229) * mrSges(4,3) + ((-t197 * mrSges(4,3) + t231 * t146 - t228 * t147 + m(5) * (-t228 * t89 + t231 * t90) + m(4) * t171) * t232 + (t198 * mrSges(4,3) + t141 - (m(5) + m(4)) * t338) * t229) * qJD(3)) * pkin(2) + t235 * t216 + m(6) * (t119 * t203 - t130 * t5 + t131 * t6 + t205 * t60 - t30 * t87 + t31 * t88) + m(7) * (t1 * t131 + t123 * t8 + t130 * t3 + t26 * t88 + t29 * t87 + t46 * t91) + t91 * t74 + t123 * t21 + t203 * t75 + t205 * t22 - t352 * t217; (t203 * t205 + t265) * t330 + (t123 * t91 + t265) * t329 + t241 + (-0.2e1 * t298 + 0.2e1 * t284 - 0.2e1 * t300 + (t216 * t336 + t217 * t229) * t331 + 0.2e1 * t255) * t303 - t208 * t278 + t123 * t327 + t91 * t325 + 0.2e1 * t290 + t205 * t326 + 0.2e1 * t217 * t200 + (t130 * t188 - t131 * t187 + t87 * t195 + t246 * t88) * t351; t297 * t121 + t282 * t120 + t75 * t221 + t308 * t167 + t309 * t168 + m(6) * (t119 * t221 - t120 * t30 + t121 * t31 - t167 * t5 + t168 * t6 + t218 * t60) + m(7) * (t1 * t168 + t120 * t29 + t121 * t26 + t140 * t8 + t167 * t3 + t46 * t98) + t234 + t235 * pkin(9) + t98 * t74 + t140 * t21 + t218 * t22 + t352 * pkin(3); t241 + (-t298 + t284 - t300 + m(5) * (-pkin(3) * t229 + pkin(9) * t336) + t255) * t303 + m(6) * (t203 * t218 + t205 * t221 + t245) + m(7) * (t123 * t98 + t140 * t91 + t245) + (-t208 + t312) * t278 + (t123 + t140) * t132 + (t205 + t218) * t133 + (t98 + t91) * t153 + (-pkin(3) + t217) * t200 + t290 + t340 * ((t120 + t87) * t195 - (-t121 - t88) * t246 + (t130 + t167) * t188 + (-t131 - t168) * t187); (-t208 + 0.2e1 * t312) * t278 + t241 + (t218 * t221 + t260) * t330 + (t140 * t98 + t260) * t329 + t140 * t327 + t98 * t325 - 0.2e1 * pkin(3) * t200 + t218 * t326 + (t120 * t195 + t121 * t246 + t167 * t188 - t168 * t187) * t351; -Ifges(5,5) * t269 + (t227 * t6 + t296 * t5) * t323 + t33 * t311 + m(7) * (qJD(6) * t26 + t1 * t212 + t214 * t3) + t34 * t270 + t243 + t1 * mrSges(7,3) - t3 * mrSges(7,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t27 * mrSges(5,2) + t28 * mrSges(5,1) + qJD(6) * t102 + t212 * t36 + t214 * t35 - t248 * Ifges(5,6); -mrSges(5,2) * t257 + t131 * t339 + t216 * t337 - t274 * t307 + t332 * t88 + t333 * t87 + t236; pkin(9) * t337 + t120 * t333 + t121 * t332 + t168 * t339 + t236; 0.2e1 * t334 * qJD(6); m(6) * t60 + m(7) * t8 + t21 + t22; m(6) * t203 + m(7) * t91 + t244; m(6) * t221 + m(7) * t98 + t244; 0; 0; m(7) * t3 + t35; m(7) * t87 + t179; m(7) * t120 + t179; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
