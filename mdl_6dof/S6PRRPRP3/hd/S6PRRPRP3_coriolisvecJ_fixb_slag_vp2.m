% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:12:29
% EndTime: 2018-11-23 15:12:40
% DurationCPUTime: 10.20s
% Computational Cost: add. (5398->572), mult. (14108->780), div. (0->0), fcn. (10131->10), ass. (0->262)
t343 = Ifges(6,1) + Ifges(7,1);
t342 = -Ifges(6,4) + Ifges(7,5);
t341 = Ifges(7,4) + Ifges(6,5);
t207 = sin(pkin(11));
t209 = cos(pkin(11));
t211 = sin(qJ(5));
t214 = cos(qJ(5));
t175 = t207 * t214 + t209 * t211;
t215 = cos(qJ(3));
t224 = t175 * t215;
t146 = qJD(2) * t224;
t159 = t175 * qJD(5);
t278 = t146 - t159;
t287 = t207 * t211;
t174 = -t214 * t209 + t287;
t223 = t174 * t215;
t147 = qJD(2) * t223;
t158 = t174 * qJD(5);
t277 = -t147 + t158;
t340 = -Ifges(6,6) + Ifges(7,6);
t355 = Ifges(6,3) + Ifges(7,2);
t354 = -Ifges(4,6) / 0.2e1;
t263 = t207 * qJD(3);
t212 = sin(qJ(3));
t271 = qJD(2) * t212;
t172 = t209 * t271 + t263;
t353 = t172 / 0.2e1;
t270 = qJD(2) * t215;
t352 = -t270 / 0.2e1;
t262 = qJD(2) * qJD(3);
t244 = t212 * t262;
t269 = qJD(3) * t209;
t225 = t207 * t271 - t269;
t108 = t211 * t172 + t214 * t225;
t220 = qJD(3) * t223;
t65 = -qJD(2) * t220 - qJD(5) * t108;
t219 = t214 * t172 - t211 * t225;
t221 = qJD(3) * t224;
t66 = qJD(2) * t221 + qJD(5) * t219;
t351 = t341 * t244 + t342 * t66 + t343 * t65;
t106 = Ifges(6,4) * t108;
t200 = qJD(5) - t270;
t294 = Ifges(7,5) * t108;
t339 = t341 * t200 + t343 * t219 - t106 + t294;
t213 = sin(qJ(2));
t208 = sin(pkin(6));
t274 = qJD(1) * t208;
t252 = t213 * t274;
t180 = qJD(2) * pkin(8) + t252;
t210 = cos(pkin(6));
t273 = qJD(1) * t212;
t140 = t215 * t180 + t210 * t273;
t250 = t207 * t270;
t116 = pkin(4) * t250 + t140;
t350 = -t278 * pkin(5) + t277 * qJ(6) - qJD(6) * t175 - t116;
t349 = -Ifges(4,1) / 0.2e1;
t348 = -Ifges(6,6) / 0.2e1;
t347 = Ifges(7,6) / 0.2e1;
t327 = t65 / 0.2e1;
t325 = t66 / 0.2e1;
t319 = -t108 / 0.2e1;
t318 = t108 / 0.2e1;
t312 = t200 / 0.2e1;
t346 = Ifges(4,4) * t352;
t316 = t219 / 0.2e1;
t345 = qJD(3) * t354;
t344 = Ifges(5,5) * t353;
t136 = qJD(3) * qJ(4) + t140;
t182 = -pkin(3) * t215 - qJ(4) * t212 - pkin(2);
t216 = cos(qJ(2));
t251 = t216 * t274;
t141 = qJD(2) * t182 - t251;
t68 = -t136 * t207 + t209 * t141;
t41 = -pkin(4) * t270 - pkin(9) * t172 + t68;
t69 = t209 * t136 + t207 * t141;
t48 = -pkin(9) * t225 + t69;
t14 = -t211 * t48 + t214 * t41;
t338 = qJD(6) - t14;
t171 = t209 * t182;
t282 = t209 * t212;
t113 = -pkin(9) * t282 + t171 + (-pkin(8) * t207 - pkin(4)) * t215;
t281 = t209 * t215;
t143 = pkin(8) * t281 + t207 * t182;
t286 = t207 * t212;
t126 = -pkin(9) * t286 + t143;
t337 = t211 * t113 + t214 * t126;
t233 = pkin(3) * t212 - qJ(4) * t215;
t155 = qJD(3) * t233 - qJD(4) * t212;
t268 = qJD(3) * t212;
t258 = pkin(8) * t268;
t120 = t209 * t155 + t207 * t258;
t279 = t215 * t216;
t137 = (-t207 * t279 + t209 * t213) * t274;
t336 = -t137 + t120;
t138 = (t207 * t213 + t209 * t279) * t274;
t149 = t207 * t155;
t335 = -t209 * t258 - t138 + t149;
t334 = t355 * t244 + t340 * t66 + t341 * t65;
t243 = t215 * t262;
t239 = t207 * t243;
t302 = mrSges(5,2) * t209;
t148 = mrSges(5,1) * t239 + t243 * t302;
t20 = t66 * mrSges(7,1) - t65 * mrSges(7,3);
t21 = t66 * mrSges(6,1) + t65 * mrSges(6,2);
t333 = t148 + t20 + t21;
t15 = t211 * t41 + t214 * t48;
t227 = pkin(4) * t212 - pkin(9) * t281;
t222 = t227 * qJD(3);
t118 = (t155 + t252) * qJD(2);
t168 = t212 * t180;
t272 = qJD(2) * t208;
t248 = t216 * t272;
t240 = t215 * t248;
t280 = t210 * t215;
t247 = qJD(1) * t280;
t276 = qJD(1) * t240 + qJD(3) * t247;
t95 = (qJD(4) - t168) * qJD(3) + t276;
t39 = t209 * t118 - t207 * t95;
t29 = qJD(2) * t222 + t39;
t40 = t207 * t118 + t209 * t95;
t31 = -pkin(9) * t239 + t40;
t4 = -qJD(5) * t15 - t211 * t31 + t214 * t29;
t86 = t222 + t120;
t285 = t207 * t215;
t96 = t149 + (-pkin(8) * t282 - pkin(9) * t285) * qJD(3);
t13 = -qJD(5) * t337 - t211 * t96 + t214 * t86;
t264 = qJD(5) * t214;
t266 = qJD(5) * t211;
t3 = t211 * t29 + t214 * t31 + t41 * t264 - t266 * t48;
t1 = qJ(6) * t244 + qJD(6) * t200 + t3;
t2 = -pkin(5) * t244 - t4;
t332 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t181 = -qJD(2) * pkin(2) - t251;
t255 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t256 = t347 + t348;
t257 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t7 = -pkin(5) * t200 + t338;
t9 = qJ(6) * t200 + t15;
t331 = t256 * t108 + t255 * t200 + t257 * t219 + t14 * mrSges(6,1) + t181 * mrSges(4,1) + t68 * mrSges(5,1) + t9 * mrSges(7,3) - Ifges(5,6) * t225 / 0.2e1 + Ifges(5,3) * t352 + t344 + t345 - (Ifges(4,4) * t212 + t215 * Ifges(4,2)) * qJD(2) / 0.2e1 + Ifges(6,6) * t319 + Ifges(7,6) * t318 - t140 * mrSges(4,3) - t15 * mrSges(6,2) - t69 * mrSges(5,2) - t7 * mrSges(7,1) + t341 * t316 + t355 * t312;
t330 = Ifges(5,1) / 0.2e1;
t329 = Ifges(7,5) * t327 + Ifges(7,3) * t325 + t244 * t347;
t328 = -t65 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t325 + t244 * t348;
t326 = -t66 / 0.2e1;
t317 = -t219 / 0.2e1;
t313 = -t200 / 0.2e1;
t311 = -t207 / 0.2e1;
t310 = t209 / 0.2e1;
t206 = t212 * pkin(8);
t307 = pkin(9) + qJ(4);
t49 = -mrSges(7,2) * t66 + mrSges(7,3) * t244;
t52 = -mrSges(6,2) * t244 - mrSges(6,3) * t66;
t306 = t49 + t52;
t50 = mrSges(6,1) * t244 - mrSges(6,3) * t65;
t51 = -mrSges(7,1) * t244 + t65 * mrSges(7,2);
t305 = t51 - t50;
t139 = -t168 + t247;
t177 = t233 * qJD(2);
t89 = -t139 * t207 + t209 * t177;
t67 = qJD(2) * t227 + t89;
t90 = t209 * t139 + t207 * t177;
t74 = -pkin(9) * t250 + t90;
t26 = t211 * t67 + t214 * t74;
t300 = mrSges(6,3) * t219;
t82 = mrSges(6,1) * t200 - t300;
t83 = -mrSges(7,1) * t200 + mrSges(7,2) * t219;
t304 = t82 - t83;
t301 = mrSges(6,3) * t108;
t81 = -mrSges(6,2) * t200 - t301;
t84 = -mrSges(7,2) * t108 + mrSges(7,3) * t200;
t303 = t84 + t81;
t299 = Ifges(5,1) * t172;
t298 = Ifges(5,4) * t209;
t297 = Ifges(6,4) * t219;
t295 = Ifges(5,5) * t209;
t293 = Ifges(5,2) * t207;
t292 = Ifges(5,6) * t207;
t284 = t208 * t213;
t160 = t212 * t284 - t280;
t267 = qJD(3) * t215;
t99 = t180 * t267 + (qJD(3) * t210 + t248) * t273;
t291 = t160 * t99;
t98 = -t180 * t268 + t276;
t290 = t215 * t98;
t289 = Ifges(4,5) * qJD(3);
t283 = t208 * t216;
t275 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t225 - t172 * mrSges(5,2) - mrSges(4,3) * t271;
t167 = t215 * pkin(4) * t263 + pkin(8) * t267;
t178 = pkin(4) * t286 + t206;
t265 = qJD(5) * t212;
t261 = t99 * t206;
t254 = Ifges(5,5) * t270;
t253 = Ifges(5,6) * t270;
t203 = -pkin(4) * t209 - pkin(3);
t249 = t213 * t272;
t46 = mrSges(7,1) * t108 - mrSges(7,3) * t219;
t47 = mrSges(6,1) * t108 + mrSges(6,2) * t219;
t242 = t46 + t47 - t275;
t241 = t212 * t251;
t238 = t263 / 0.2e1 - t172 / 0.2e1;
t237 = -mrSges(4,1) * t215 + mrSges(4,2) * t212;
t236 = mrSges(5,1) * t207 + t302;
t235 = Ifges(5,1) * t209 - Ifges(5,4) * t207;
t234 = -t293 + t298;
t25 = -t211 * t74 + t214 * t67;
t53 = t113 * t214 - t126 * t211;
t161 = t210 * t212 + t215 * t284;
t122 = -t161 * t207 - t209 * t283;
t123 = t161 * t209 - t207 * t283;
t229 = t214 * t122 - t123 * t211;
t56 = t122 * t211 + t123 * t214;
t184 = t307 * t207;
t185 = t307 * t209;
t228 = -t214 * t184 - t185 * t211;
t131 = -t184 * t211 + t185 * t214;
t12 = t113 * t264 - t126 * t266 + t211 * t86 + t214 * t96;
t226 = t139 * mrSges(4,3) + t271 * t349 + t346 - t289 / 0.2e1 - t181 * mrSges(4,2);
t132 = -qJD(3) * pkin(3) + qJD(4) - t139;
t76 = pkin(4) * t239 + t99;
t91 = pkin(4) * t225 + t132;
t217 = qJD(2) ^ 2;
t190 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t270;
t176 = t237 * qJD(2);
t166 = (mrSges(4,1) * t212 + mrSges(4,2) * t215) * t262;
t154 = (mrSges(5,1) * t212 - mrSges(5,3) * t281) * t262;
t153 = (-mrSges(5,2) * t212 - mrSges(5,3) * t285) * t262;
t152 = t174 * t212;
t151 = t175 * t212;
t145 = -mrSges(5,1) * t270 - mrSges(5,3) * t172;
t144 = mrSges(5,2) * t270 - mrSges(5,3) * t225;
t142 = -pkin(8) * t285 + t171;
t128 = (Ifges(5,5) * t212 + t215 * t235) * t262;
t127 = (Ifges(5,6) * t212 + t215 * t234) * t262;
t125 = qJD(3) * t161 + t212 * t248;
t124 = -qJD(3) * t160 + t240;
t105 = Ifges(7,5) * t219;
t104 = pkin(5) * t174 - qJ(6) * t175 + t203;
t102 = -Ifges(5,4) * t225 - t254 + t299;
t101 = Ifges(5,4) * t172 - Ifges(5,2) * t225 - t253;
t93 = t264 * t282 - t265 * t287 + t221;
t92 = -t175 * t265 - t220;
t88 = t124 * t209 + t207 * t249;
t87 = -t124 * t207 + t209 * t249;
t80 = qJD(4) * t175 + qJD(5) * t131;
t79 = -qJD(4) * t174 + qJD(5) * t228;
t73 = pkin(5) * t151 + qJ(6) * t152 + t178;
t71 = t137 * t211 + t138 * t214;
t70 = -t214 * t137 + t138 * t211;
t45 = pkin(5) * t219 + qJ(6) * t108;
t44 = pkin(5) * t215 - t53;
t43 = -qJ(6) * t215 + t337;
t36 = -Ifges(6,2) * t108 + Ifges(6,6) * t200 + t297;
t33 = Ifges(7,6) * t200 + Ifges(7,3) * t108 + t105;
t27 = t108 * pkin(5) - qJ(6) * t219 + t91;
t24 = -pkin(5) * t271 - t25;
t23 = qJ(6) * t271 + t26;
t22 = pkin(5) * t93 - qJ(6) * t92 + qJD(6) * t152 + t167;
t11 = qJD(5) * t56 + t211 * t88 - t214 * t87;
t10 = qJD(5) * t229 + t211 * t87 + t214 * t88;
t8 = -pkin(5) * t268 - t13;
t6 = qJ(6) * t268 - qJD(6) * t215 + t12;
t5 = pkin(5) * t66 - qJ(6) * t65 - qJD(6) * t219 + t76;
t16 = [-t161 * mrSges(4,3) * t244 + t122 * t154 + t123 * t153 + t124 * t190 + t88 * t144 + t87 * t145 + t306 * t56 - t305 * t229 - t304 * t11 + t303 * t10 + ((-mrSges(3,2) * t217 - t166) * t216 + (-mrSges(3,1) * t217 + qJD(2) * t176) * t213) * t208 + (mrSges(4,3) * t243 + t333) * t160 + t242 * t125 + m(7) * (t1 * t56 + t10 * t9 + t11 * t7 + t125 * t27 + t160 * t5 - t2 * t229) + m(6) * (t10 * t15 - t11 * t14 + t125 * t91 + t160 * t76 + t229 * t4 + t3 * t56) + m(5) * (t122 * t39 + t123 * t40 + t125 * t132 + t68 * t87 + t69 * t88 + t291) + m(4) * (t124 * t140 - t125 * t139 + t291 + t161 * t98 + (t181 - t251) * t249); (-t39 * mrSges(5,1) + t40 * mrSges(5,2) - Ifges(6,6) * t326 - Ifges(7,6) * t325 - t341 * t327 + (t235 * t353 + t132 * t236 + (Ifges(4,5) / 0.2e1 + t234 * t310) * qJD(3) + t101 * t311 + t102 * t310 + (-t207 * t69 - t209 * t68) * mrSges(5,3) + (-m(4) * t139 + m(5) * t132 - t275) * pkin(8) - t226 + ((-0.3e1 / 0.2e1 * t295 + 0.3e1 / 0.2e1 * t292 + 0.3e1 / 0.2e1 * Ifges(4,4)) * t215 + (t209 ^ 2 * t330 - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(5,3) + (-0.3e1 / 0.2e1 * t298 + t293) * t207 - t255) * t212) * qJD(2)) * qJD(3) + t332) * t215 + ((-m(4) * t140 - t190) * pkin(8) + (Ifges(5,6) * t310 + t354) * qJD(3) + t344 + t331) * t268 - t351 * t152 / 0.2e1 + (-t15 * mrSges(6,3) - t9 * mrSges(7,2) + t91 * mrSges(6,1) - t36 / 0.2e1 + t27 * mrSges(7,1) + t33 / 0.2e1 + Ifges(7,3) * t318 - Ifges(6,2) * t319 + t340 * t312 + t342 * t316) * t93 + (mrSges(6,2) * t91 + mrSges(7,2) * t7 - mrSges(6,3) * t14 - mrSges(7,3) * t27 + Ifges(6,4) * t319 + Ifges(7,5) * t318 + t339 / 0.2e1 + t341 * t312 + t343 * t316) * t92 + (-t176 * t213 - t190 * t279 - m(4) * (t181 * t213 + (-t139 * t212 + t140 * t215) * t216)) * t274 + (-t1 * t151 - t152 * t2) * mrSges(7,2) + (-t151 * t3 + t152 * t4) * mrSges(6,3) + (-Ifges(7,5) * t152 + Ifges(7,3) * t151) * t325 + (-Ifges(6,4) * t152 - Ifges(6,2) * t151) * t326 + t5 * (mrSges(7,1) * t151 + mrSges(7,3) * t152) + t76 * (mrSges(6,1) * t151 - mrSges(6,2) * t152) + (t1 * t43 + t2 * t44 + t5 * t73 + (t6 - t71) * t9 + (-t70 + t8) * t7 + (t22 - t241) * t27) * m(7) + (t178 * t76 + t3 * t337 + t4 * t53 + (t167 - t241) * t91 + (t12 - t71) * t15 + (t13 + t70) * t14) * m(6) + t337 * t52 + (pkin(8) * t148 + t127 * t311 + t128 * t310 + (mrSges(4,3) + t236) * t99 + (-t207 * t40 - t209 * t39) * mrSges(5,3) - t242 * t251) * t212 - t303 * t71 + t304 * t70 + m(4) * (pkin(8) * t290 + t261) + t178 * t21 - pkin(2) * t166 + t167 * t47 + t143 * t153 + t142 * t154 + t13 * t82 + t8 * t83 + t6 * t84 + t12 * t81 + t73 * t20 + ((-m(4) * pkin(2) + t237) * t252 + ((-0.3e1 / 0.2e1 * Ifges(4,4) - t292 + t295 / 0.2e1) * t212 - t257 * t152 + t256 * t151) * t268) * qJD(2) - t334 * t215 / 0.2e1 + t335 * t144 + t336 * t145 + (-t132 * t241 + t142 * t39 + t143 * t40 + t335 * t69 + t336 * t68 + t261) * m(5) + mrSges(4,3) * t290 + t151 * t328 + t151 * t329 + (t151 * t342 - t152 * t343) * t327 + t22 * t46 + t43 * t49 + t44 * t51 + t53 * t50; (-t132 * t140 - t68 * t89 - t69 * t90 - pkin(3) * t99 + (-t207 * t68 + t209 * t69) * qJD(4) + (-t207 * t39 + t209 * t40) * qJ(4)) * m(5) + (t340 * t146 - t341 * t147) * t313 + (-t341 * t158 + t340 * t159) * t312 + (t342 * t146 - t343 * t147) * t317 + (t342 * t174 + t343 * t175) * t327 + (-t343 * t158 + t342 * t159) * t316 + (t1 * t131 + t104 * t5 - t228 * t2 + (-t23 + t79) * t9 + (-t24 + t80) * t7 + t350 * t27) * m(7) + t350 * t46 + t339 * (-t158 / 0.2e1 + t147 / 0.2e1) + t351 * t175 / 0.2e1 + (-t116 * t91 + t228 * t4 + t131 * t3 + t203 * t76 + (-t26 + t79) * t15 + (-t25 - t80) * t14) * m(6) - t305 * t228 + (-Ifges(6,4) * t158 - Ifges(7,5) * t147 - Ifges(6,2) * t159 + Ifges(7,3) * t146) * t319 + (-Ifges(6,4) * t147 - Ifges(7,5) * t158 - Ifges(6,2) * t146 + Ifges(7,3) * t159) * t318 + (t99 * mrSges(5,2) + t128 / 0.2e1 - t39 * mrSges(5,3) - qJD(4) * t145 - qJ(4) * t154) * t207 + (-t99 * mrSges(5,1) + t127 / 0.2e1 + t40 * mrSges(5,3) + qJD(4) * t144 + qJ(4) * t153) * t209 + (t33 - t36) * (t159 / 0.2e1 - t146 / 0.2e1) + t303 * t79 - t304 * t80 + t306 * t131 + (-mrSges(7,1) * t278 + mrSges(7,3) * t277) * t27 + (-mrSges(6,1) * t278 - mrSges(6,2) * t277) * t91 + t275 * t140 + (t14 * t277 + t15 * t278 - t174 * t3 - t175 * t4) * mrSges(6,3) + (-t1 * t174 + t175 * t2 - t277 * t7 + t278 * t9) * mrSges(7,2) + t203 * t21 - t139 * t190 + t5 * (mrSges(7,1) * t174 - mrSges(7,3) * t175) + t76 * (mrSges(6,1) * t174 + mrSges(6,2) * t175) - pkin(3) * t148 - t90 * t144 - t89 * t145 - t116 * t47 + t104 * t20 - t98 * mrSges(4,2) - t99 * mrSges(4,1) - t24 * t83 - t23 * t84 - t26 * t81 - t25 * t82 + (Ifges(7,5) * t175 + Ifges(7,3) * t174) * t325 + (Ifges(6,4) * t175 - Ifges(6,2) * t174) * t326 + t174 * t328 + t174 * t329 + ((t346 + t289 / 0.2e1 + (-t299 / 0.2e1 + t68 * mrSges(5,3) - t132 * mrSges(5,2) - t102 / 0.2e1 + t254 / 0.2e1) * t209 + (t69 * mrSges(5,3) - t132 * mrSges(5,1) + t101 / 0.2e1 + t269 * t330 - t253 / 0.2e1 - t238 * Ifges(5,4)) * t207 + t226) * t215 + (t238 * Ifges(5,5) + (t292 / 0.2e1 + Ifges(4,4) / 0.2e1) * t271 + (t207 * t234 / 0.2e1 + t349 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t270 + t345 + (t174 * t340 + t175 * t341) * qJD(3) / 0.2e1 - t331) * t212) * qJD(2); t304 * t219 + t303 * t108 + t225 * t144 + t172 * t145 + (t108 * t9 - t219 * t7 + t5) * m(7) + (t108 * t15 + t14 * t219 + t76) * m(6) + (t68 * t172 + t225 * t69 + t99) * m(5) + t333; (-t108 * t343 + t105 - t297 + t33) * t317 + (t108 * t7 + t219 * t9) * mrSges(7,2) - t27 * (mrSges(7,1) * t219 + mrSges(7,3) * t108) - t91 * (mrSges(6,1) * t219 - mrSges(6,2) * t108) + (Ifges(7,3) * t219 - t294) * t319 + (-t108 * t341 + t219 * t340) * t313 + (-Ifges(6,2) * t219 - t106 + t339) * t318 - t332 + (-pkin(5) * t2 + qJ(6) * t1 - t15 * t7 - t27 * t45 + t338 * t9) * m(7) + (-t301 - t303) * t14 + (t300 + t304) * t15 + qJD(6) * t84 + t36 * t316 - t45 * t46 + qJ(6) * t49 - pkin(5) * t51 + t334; t219 * t46 - t200 * t84 + 0.2e1 * (t2 / 0.2e1 + t27 * t316 + t9 * t313) * m(7) + t51;];
tauc  = t16(:);
