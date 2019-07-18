% Calculate vector of inverse dynamics joint torques for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:45
% EndTime: 2019-07-18 17:21:16
% DurationCPUTime: 10.06s
% Computational Cost: add. (3823->468), mult. (8348->645), div. (0->0), fcn. (5290->10), ass. (0->232)
t345 = Ifges(3,2) + Ifges(4,2);
t237 = m(4) * pkin(1) + mrSges(4,1);
t362 = mrSges(3,1) + t237;
t347 = Ifges(3,4) + Ifges(4,4);
t361 = Ifges(3,1) + Ifges(4,1);
t346 = Ifges(3,5) + Ifges(4,5);
t344 = Ifges(4,6) + Ifges(3,6);
t171 = cos(qJ(2));
t167 = sin(qJ(2));
t290 = Ifges(4,4) * t167;
t292 = Ifges(3,4) * t167;
t360 = t345 * t171 + t290 + t292;
t277 = t167 * mrSges(4,2);
t163 = qJ(2) + qJ(4);
t158 = sin(t163);
t159 = cos(t163);
t334 = -t159 * mrSges(5,1) + mrSges(5,2) * t158;
t359 = mrSges(3,2) * t167 - t171 * t362 + t277 + t334;
t242 = qJD(1) * qJD(2);
t131 = qJDD(1) * t171 - t167 * t242;
t358 = t131 / 0.2e1;
t251 = qJD(1) * t171;
t357 = t347 * t251;
t170 = cos(qJ(4));
t299 = pkin(3) + qJ(3);
t138 = t299 * t167;
t126 = qJD(1) * t138;
t173 = pkin(2) + pkin(1);
t181 = qJD(2) * t173 - t126;
t139 = t299 * t171;
t128 = qJD(1) * t139;
t166 = sin(qJ(4));
t266 = t128 * t166;
t356 = t170 * t181 - t266;
t132 = qJDD(1) * t167 + t171 * t242;
t125 = t166 * t171 + t167 * t170;
t188 = t125 * qJD(4);
t62 = -qJD(1) * t188 + t131 * t170 - t132 * t166;
t165 = sin(qJ(5));
t282 = t165 * mrSges(6,2);
t355 = -mrSges(6,3) * t159 - t158 * t282;
t169 = cos(qJ(5));
t298 = mrSges(6,1) * t169;
t354 = -t282 + t298;
t168 = sin(qJ(1));
t172 = cos(qJ(1));
t353 = g(1) * t172 + g(2) * t168;
t340 = t171 * t173;
t352 = -m(5) * t340 - m(6) * (pkin(4) * t158 + t340) + t359;
t124 = t166 * t167 - t170 * t171;
t118 = t124 * qJD(1);
t111 = qJD(5) + t118;
t160 = qJDD(2) + qJDD(4);
t187 = t124 * qJD(4);
t61 = -qJD(1) * t187 + t131 * t166 + t132 * t170;
t119 = t125 * qJD(1);
t161 = qJD(2) + qJD(4);
t92 = -t119 * t165 + t161 * t169;
t25 = qJD(5) * t92 + t160 * t165 + t169 * t61;
t325 = t25 / 0.2e1;
t94 = t119 * t169 + t161 * t165;
t26 = -qJD(5) * t94 + t160 * t169 - t165 * t61;
t324 = t26 / 0.2e1;
t309 = Ifges(6,4) * t94;
t36 = t92 * Ifges(6,2) + t111 * Ifges(6,6) + t309;
t351 = -t36 / 0.2e1;
t60 = qJDD(5) - t62;
t322 = t60 / 0.2e1;
t350 = m(5) + m(6);
t349 = -mrSges(5,1) * t160 - mrSges(6,1) * t26 + mrSges(6,2) * t25 + mrSges(5,3) * t61;
t286 = t119 * mrSges(5,3);
t343 = -mrSges(5,1) * t161 - mrSges(6,1) * t92 + mrSges(6,2) * t94 + t286;
t208 = mrSges(6,1) * t165 + mrSges(6,2) * t169;
t342 = t208 * t356;
t63 = -mrSges(6,2) * t111 + mrSges(6,3) * t92;
t64 = mrSges(6,1) * t111 - mrSges(6,3) * t94;
t198 = -t165 * t64 + t169 * t63;
t294 = mrSges(5,3) * t118;
t98 = -mrSges(5,2) * t161 - t294;
t341 = -t198 - t98;
t339 = qJD(1) * t360 + qJD(2) * t344;
t252 = qJD(1) * t167;
t338 = t346 * qJD(2) + t252 * t361 + t357;
t337 = -t167 * t344 + t171 * t346;
t244 = qJD(5) * t169;
t87 = -qJD(2) * t124 - t187;
t193 = t125 * t244 + t165 * t87;
t11 = mrSges(6,1) * t60 - mrSges(6,3) * t25;
t12 = -mrSges(6,2) * t60 + mrSges(6,3) * t26;
t336 = -t165 * t11 + t169 * t12;
t335 = t347 * t171;
t245 = qJD(5) * t165;
t241 = qJD(1) * qJD(3);
t225 = t167 * t241;
t177 = qJDD(2) * t173 - t132 * t299 - t225;
t102 = t131 * qJ(3) + t171 * t241;
t82 = pkin(3) * t131 + t102;
t16 = qJD(4) * t356 + t166 * t177 + t170 * t82;
t14 = pkin(4) * t160 + t16;
t257 = t170 * t128;
t76 = t166 * t181 + t257;
t68 = t161 * pkin(4) + t76;
t122 = -qJD(1) * t340 + qJD(3);
t86 = -pkin(4) * t119 + t122;
t32 = t165 * t86 + t169 * t68;
t90 = -t131 * t173 + qJDD(3);
t33 = -pkin(4) * t61 + t90;
t3 = -qJD(5) * t32 - t14 * t165 + t169 * t33;
t302 = t3 * t165;
t31 = -t165 * t68 + t169 * t86;
t331 = -t31 * t244 - t32 * t245 - t302;
t330 = -m(4) * qJ(3) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t2 = qJD(5) * t31 + t14 * t169 + t165 * t33;
t329 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t174 = qJD(1) ^ 2;
t327 = Ifges(6,1) * t325 + Ifges(6,4) * t324 + Ifges(6,5) * t322;
t321 = -t92 / 0.2e1;
t320 = -t94 / 0.2e1;
t319 = t94 / 0.2e1;
t318 = -t111 / 0.2e1;
t315 = t119 / 0.2e1;
t311 = t169 / 0.2e1;
t145 = -pkin(1) * t251 + qJD(3);
t308 = pkin(1) * t145;
t307 = pkin(1) * t171;
t304 = g(3) * t158;
t303 = t169 * t2;
t83 = -t126 * t166 + t257;
t301 = t356 * t83;
t297 = mrSges(4,2) * t171;
t295 = mrSges(5,2) * t159;
t288 = Ifges(6,4) * t165;
t287 = Ifges(6,4) * t169;
t285 = t119 * Ifges(5,4);
t284 = t119 * t356;
t17 = qJD(4) * t76 + t166 * t82 - t170 * t177;
t283 = t125 * t17;
t147 = t158 * mrSges(6,3);
t278 = t166 * t356;
t96 = qJDD(2) * pkin(1) - qJ(3) * t132 - t225;
t276 = t167 * t96;
t273 = t17 * t170;
t271 = qJ(3) * t171;
t270 = t118 * t165;
t269 = t118 * t169;
t268 = t125 * t165;
t267 = t125 * t169;
t265 = t158 * t172;
t264 = t165 * t168;
t263 = t165 * t172;
t262 = t166 * t173;
t261 = t167 * t171;
t260 = t167 * t173;
t259 = t168 * t169;
t258 = t169 * t172;
t254 = t355 * t168;
t253 = t355 * t172;
t154 = pkin(1) * t252;
t129 = pkin(2) * t252 + t154;
t250 = qJD(2) * t167;
t156 = pkin(1) * t250;
t130 = pkin(2) * t250 + t156;
t249 = qJD(2) * t171;
t248 = qJD(3) * t167;
t247 = qJD(3) * t171;
t239 = Ifges(6,5) * t25 + Ifges(6,6) * t26 + Ifges(6,3) * t60;
t136 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t252;
t236 = t136 * t271;
t235 = mrSges(4,3) * t251;
t91 = Ifges(6,4) * t92;
t37 = t94 * Ifges(6,1) + t111 * Ifges(6,5) + t91;
t232 = t37 * t311;
t231 = qJ(3) * t252;
t226 = -t62 * mrSges(5,1) + t61 * mrSges(5,2);
t220 = -t245 / 0.2e1;
t219 = qJD(2) * t299;
t216 = -t131 * mrSges(4,1) + t132 * mrSges(4,2);
t215 = t168 * t299 + t172 * t340;
t214 = qJ(3) * mrSges(4,3) * t261;
t196 = -t170 * t138 - t139 * t166;
t110 = -t167 * t219 + t247;
t184 = -t171 * t219 - t248;
t95 = -t138 * t166 + t139 * t170;
t39 = qJD(4) * t95 + t110 * t166 - t170 * t184;
t212 = -t17 * t196 - t356 * t39;
t209 = -mrSges(4,1) * t171 + t277;
t207 = Ifges(6,1) * t169 - t288;
t204 = -Ifges(6,2) * t165 + t287;
t201 = Ifges(6,5) * t169 - Ifges(6,6) * t165;
t199 = -t165 * t31 + t169 * t32;
t97 = -pkin(4) * t125 - t340;
t48 = t165 * t97 + t169 * t95;
t47 = -t165 * t95 + t169 * t97;
t195 = t159 * t354 + t147;
t192 = t125 * t245 - t169 * t87;
t191 = t145 * (mrSges(4,1) * t167 + t297);
t190 = t167 * (Ifges(3,1) * t171 - t292);
t189 = t167 * (Ifges(4,1) * t171 - t290);
t183 = t295 + (mrSges(5,1) + t298) * t158;
t180 = -t302 + (-t165 * t32 - t169 * t31) * qJD(5);
t179 = m(6) * (pkin(4) * t159 - t260) - t158 * t298;
t109 = Ifges(5,4) * t118;
t35 = t94 * Ifges(6,5) + t92 * Ifges(6,6) + t111 * Ifges(6,3);
t6 = t25 * Ifges(6,4) + t26 * Ifges(6,2) + t60 * Ifges(6,6);
t72 = -t118 * Ifges(5,2) + t161 * Ifges(5,6) + t285;
t73 = t119 * Ifges(5,1) + t161 * Ifges(5,5) - t109;
t176 = (Ifges(6,5) * t165 + Ifges(6,6) * t169) * t322 + (Ifges(6,2) * t169 + t288) * t324 + (Ifges(6,1) * t165 + t287) * t325 + t165 * t327 + t6 * t311 + t72 * t315 + t76 * t286 + mrSges(6,3) * t303 + t36 * t220 - t31 * (mrSges(6,1) * t119 + mrSges(6,3) * t269) - t32 * (-mrSges(6,2) * t119 + mrSges(6,3) * t270) + (t232 - t342) * qJD(5) + (Ifges(6,3) * t119 - t118 * t201) * t318 + (Ifges(6,5) * t119 - t118 * t207) * t320 + (Ifges(6,6) * t119 - t118 * t204) * t321 - t122 * (mrSges(5,1) * t119 - mrSges(5,2) * t118) - t161 * (-Ifges(5,5) * t118 - Ifges(5,6) * t119) / 0.2e1 - t356 * t294 + (t111 * t201 + t204 * t92 + t207 * t94) * qJD(5) / 0.2e1 - (-Ifges(5,1) * t118 - t285 + t35) * t119 / 0.2e1 + (-Ifges(5,2) * t119 - t109 + t73) * t118 / 0.2e1 + (-mrSges(5,1) - t354) * t17 + Ifges(5,3) * t160 + Ifges(5,5) * t61 + Ifges(5,6) * t62 + t37 * t269 / 0.2e1 - t16 * mrSges(5,2) - t118 * t342 + t270 * t351;
t175 = qJ(3) ^ 2;
t162 = t171 ^ 2;
t137 = -qJD(2) * mrSges(4,2) + t235;
t135 = qJD(2) * pkin(1) - t231;
t127 = t209 * qJD(1);
t112 = -pkin(1) * t131 + qJDD(3);
t108 = t159 * t258 + t264;
t107 = -t159 * t263 + t259;
t106 = -t159 * t259 + t263;
t105 = t159 * t264 + t258;
t104 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t132;
t89 = pkin(4) * t118 + t129;
t88 = qJD(2) * t125 + t188;
t84 = -t126 * t170 - t266;
t81 = mrSges(5,1) * t118 + mrSges(5,2) * t119;
t65 = -pkin(4) * t87 + t130;
t51 = -mrSges(5,2) * t160 + mrSges(5,3) * t62;
t43 = pkin(4) * t270 + t169 * t356;
t42 = pkin(4) * t269 - t165 * t356;
t41 = t165 * t89 + t169 * t84;
t40 = -t165 * t84 + t169 * t89;
t38 = qJD(4) * t196 + t170 * t110 + t166 * t184;
t10 = -qJD(5) * t48 - t165 * t38 + t169 * t65;
t9 = qJD(5) * t47 + t165 * t65 + t169 * t38;
t1 = [t335 * t132 / 0.2e1 + (-t137 * t250 + m(4) * (t162 * t241 - t276 + (-qJD(2) * t135 + t102) * t171)) * qJ(3) + (-Ifges(6,1) * t192 - Ifges(6,4) * t193 + Ifges(6,5) * t88) * t319 + t267 * t327 + (Ifges(5,1) * t87 - Ifges(5,4) * t88) * t315 + t208 * t283 + t137 * t247 + qJD(2) * t191 + (t102 * t171 + t131 * t271 - t135 * t249 - t276) * mrSges(4,3) + (-m(4) * t307 + t209) * t112 + t87 * t232 + (-t356 * t87 - t76 * t88 + t283) * mrSges(5,3) - t356 * (mrSges(6,1) * t193 - mrSges(6,2) * t192) + ((-t167 * t345 + t335) * t171 + t190 + t189) * t242 / 0.2e1 - t349 * t196 + m(5) * (t122 * t130 + t16 * t95 - t340 * t90 + t38 * t76 + t212) - t340 * t226 + (-mrSges(4,2) * t271 + t167 * t346 + t171 * t344) * qJDD(2) + t343 * t39 + t337 * qJD(2) ^ 2 / 0.2e1 + t338 * t249 / 0.2e1 - t339 * t250 / 0.2e1 + (t192 * t31 - t193 * t32 - t2 * t268 - t267 * t3) * mrSges(6,3) + (t131 * t345 + t132 * t347) * t171 / 0.2e1 + t31 * mrSges(6,1) * t88 - t32 * mrSges(6,2) * t88 + t127 * t156 + (-m(5) * t215 - m(6) * (pkin(4) * t265 + t215) - t108 * mrSges(6,1) - t107 * mrSges(6,2) - mrSges(6,3) * t265 + (-mrSges(2,1) + t359) * t172 + t330 * t168) * g(2) + t360 * t358 + (Ifges(6,3) * t322 + Ifges(6,6) * t324 + Ifges(6,5) * t325 - Ifges(5,6) * t160 - Ifges(5,4) * t61 - Ifges(5,2) * t62 + t90 * mrSges(5,1) - t16 * mrSges(5,3) + t239 / 0.2e1 + t329) * t124 + (t90 * mrSges(5,2) + Ifges(5,1) * t61 + Ifges(5,4) * t62 + Ifges(5,5) * t160 + t201 * t322 + t204 * t324 + t207 * t325 + t220 * t37) * t125 + (t347 * t358 + t361 * t132 + m(4) * (-qJD(3) * t135 + (-t175 * t251 + t308) * qJD(2)) - qJ(3) * t104) * t167 + t161 * (Ifges(5,5) * t87 - Ifges(5,6) * t88) / 0.2e1 + t130 * t81 + t122 * (mrSges(5,1) * t88 + mrSges(5,2) * t87) - t216 * t307 - t118 * (Ifges(5,4) * t87 - Ifges(5,2) * t88) / 0.2e1 + t95 * t51 + t38 * t98 + t87 * t73 / 0.2e1 + t88 * t35 / 0.2e1 - t88 * t72 / 0.2e1 + t9 * t63 + t10 * t64 - t6 * t268 / 0.2e1 + t47 * t11 + t48 * t12 - t136 * t248 - t214 * t242 + t193 * t351 - qJD(2) * t236 + (-t106 * mrSges(6,1) - t105 * mrSges(6,2) + (-t299 * t350 + t330) * t172 + (mrSges(2,1) + t147 - t352) * t168) * g(1) + t111 * (-Ifges(6,5) * t192 - Ifges(6,6) * t193 + Ifges(6,3) * t88) / 0.2e1 + t92 * (-Ifges(6,4) * t192 - Ifges(6,2) * t193 + Ifges(6,6) * t88) / 0.2e1 + m(6) * (t10 * t31 + t2 * t48 + t3 * t47 + t32 * t9 + t212) + Ifges(2,3) * qJDD(1); t51 * t262 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (-t189 / 0.2e1 - t190 / 0.2e1 + m(4) * t175 * t261 + t214) * t174 + ((t166 * t343 - t170 * t341) * qJD(4) - t349 * t170 + (-t273 + (t170 * t199 - t278) * qJD(4)) * m(6) + (t16 * t166 - t273 + (t170 * t76 - t278) * qJD(4)) * m(5)) * t173 + (-t31 * t40 - t32 * t41 + t301) * m(6) + (-t122 * t129 - t76 * t84 + t301) * m(5) + (t236 - m(4) * (-t135 * t271 + t167 * t308) - t191) * qJD(1) + t176 - t343 * t83 + t344 * t131 + t346 * t132 + t237 * t96 - t337 * t242 / 0.2e1 + t339 * t252 / 0.2e1 + t353 * (m(5) * t260 + mrSges(5,1) * t158 + mrSges(3,2) * t171 + t167 * t362 + t295 + t297) + t137 * t231 + t135 * t235 + (-t64 * t244 - t63 * t245 + m(6) * (t180 + t303) + t336) * (pkin(4) + t262) + t331 * mrSges(6,3) - (-t252 * t345 + t338 + t357) * t251 / 0.2e1 - t129 * t81 + pkin(1) * t104 - t102 * mrSges(4,2) - t84 * t98 - t41 * t63 - t40 * t64 - g(1) * (t172 * t179 - t253) - g(2) * (t168 * t179 - t254) - t127 * t154 + (-t195 + t352) * g(3); t169 * t11 + t165 * t12 - t343 * t119 + t198 * qJD(5) + (t167 * t136 - t171 * t137) * qJD(1) - t341 * t118 + t216 + t226 + (-g(1) * t168 + g(2) * t172) * (m(4) + t350) + (t111 * t199 + t165 * t2 + t169 * t3 + t284) * m(6) + (t118 * t76 + t284 + t90) * m(5) + (-qJ(3) * t162 * t174 + t135 * t252 + t112) * m(4); ((-t165 * t63 - t169 * t64) * qJD(5) + (-t159 * t353 + t303 - t304 + t331) * m(6) + t336) * pkin(4) - m(6) * (t31 * t42 + t32 * t43) + t176 + t180 * mrSges(6,3) + (t172 * t183 + t253) * g(1) + (t168 * t183 + t254) * g(2) + (-t195 + t334) * g(3) - t356 * t98 - t43 * t63 - t42 * t64 + (m(6) * t356 - t343) * t76; t356 * (mrSges(6,1) * t94 + mrSges(6,2) * t92) + (Ifges(6,1) * t92 - t309) * t320 + t36 * t319 + (Ifges(6,5) * t92 - Ifges(6,6) * t94) * t318 - t31 * t63 + t32 * t64 - g(1) * (mrSges(6,1) * t107 - mrSges(6,2) * t108) - g(2) * (-mrSges(6,1) * t105 + mrSges(6,2) * t106) + t208 * t304 + (t31 * t92 + t32 * t94) * mrSges(6,3) + t239 + (-Ifges(6,2) * t94 + t37 + t91) * t321 + t329;];
tau  = t1;
