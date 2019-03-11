% Calculate time derivative of joint inertia matrix for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:18
% EndTime: 2019-03-10 03:33:27
% DurationCPUTime: 5.76s
% Computational Cost: add. (17901->518), mult. (38429->768), div. (0->0), fcn. (38751->10), ass. (0->247)
t226 = sin(qJ(3));
t227 = sin(qJ(2));
t231 = cos(qJ(3));
t232 = cos(qJ(2));
t279 = t231 * t232;
t187 = -t226 * t227 + t279;
t189 = t226 * t232 + t231 * t227;
t225 = sin(qJ(4));
t230 = cos(qJ(4));
t145 = t187 * t225 + t189 * t230;
t224 = sin(qJ(5));
t273 = qJD(5) * t224;
t260 = t145 * t273;
t229 = cos(qJ(5));
t332 = qJD(2) + qJD(3);
t153 = t332 * t187;
t154 = t332 * t189;
t245 = t230 * t187 - t189 * t225;
t76 = qJD(4) * t245 + t153 * t230 - t154 * t225;
t293 = t229 * t76;
t240 = t260 - t293;
t272 = qJD(5) * t229;
t295 = t224 * t76;
t241 = t145 * t272 + t295;
t32 = mrSges(6,1) * t241 - mrSges(6,2) * t240;
t321 = -pkin(8) - pkin(7);
t206 = t321 * t227;
t190 = t226 * t206;
t208 = t321 * t232;
t163 = -t231 * t208 + t190;
t237 = (t279 * t321 - t190) * qJD(2);
t118 = -qJD(3) * t163 + t237;
t313 = t153 * pkin(9);
t276 = qJD(3) * t231;
t277 = qJD(3) * t226;
t117 = qJD(2) * t189 * t321 + t206 * t276 + t208 * t277;
t88 = -pkin(9) * t154 + t117;
t161 = t231 * t206 + t208 * t226;
t133 = -pkin(9) * t189 + t161;
t134 = pkin(9) * t187 + t163;
t92 = t133 * t225 + t134 * t230;
t37 = qJD(4) * t92 + t225 * t88 - t230 * (t118 - t313);
t342 = m(6) * t37 + t32;
t221 = t224 ^ 2;
t222 = t229 ^ 2;
t341 = t221 + t222;
t223 = sin(qJ(6));
t228 = cos(qJ(6));
t244 = t223 * t224 - t228 * t229;
t331 = qJD(5) + qJD(6);
t151 = t331 * t244;
t188 = t223 * t229 + t224 * t228;
t152 = t331 * t188;
t278 = -Ifges(7,5) * t151 - Ifges(7,6) * t152;
t340 = Ifges(6,5) * t272 + t278;
t334 = t230 * t133 - t134 * t225;
t36 = t230 * t88 + t334 * qJD(4) + (-t206 * t277 + t208 * t276 + t237 - t313) * t225;
t135 = qJD(2) * t227 * pkin(2) + pkin(3) * t154;
t77 = qJD(4) * t145 + t153 * t225 + t230 * t154;
t43 = pkin(4) * t77 - pkin(10) * t76 + t135;
t216 = -pkin(2) * t232 - pkin(1);
t167 = -t187 * pkin(3) + t216;
t90 = -pkin(4) * t245 - t145 * pkin(10) + t167;
t13 = t224 * t43 + t229 * t36 + t90 * t272 - t273 * t92;
t256 = -t224 * t36 + t229 * t43;
t89 = t229 * t92;
t55 = t224 * t90 + t89;
t14 = -qJD(5) * t55 + t256;
t339 = t13 * t229 - t14 * t224;
t338 = Ifges(6,5) * t293 + Ifges(6,3) * t77;
t214 = pkin(2) * t231 + pkin(3);
t283 = t226 * t230;
t174 = pkin(2) * t283 + t225 * t214;
t172 = pkin(10) + t174;
t309 = -pkin(11) - t172;
t255 = qJD(5) * t309;
t274 = qJD(4) * t230;
t275 = qJD(4) * t225;
t284 = t225 * t226;
t137 = t214 * t274 + (-t226 * t275 + (t230 * t231 - t284) * qJD(3)) * pkin(2);
t281 = t229 * t137;
t105 = t224 * t255 + t281;
t285 = t224 * t137;
t106 = t229 * t255 - t285;
t158 = t309 * t224;
t220 = t229 * pkin(11);
t159 = t172 * t229 + t220;
t119 = t158 * t228 - t159 * t223;
t49 = qJD(6) * t119 + t105 * t228 + t106 * t223;
t120 = t158 * t223 + t159 * t228;
t50 = -qJD(6) * t120 - t105 * t223 + t106 * t228;
t337 = t50 * mrSges(7,1) - t49 * mrSges(7,2);
t212 = pkin(3) * t225 + pkin(10);
t308 = -pkin(11) - t212;
t183 = t308 * t224;
t184 = t212 * t229 + t220;
t139 = t183 * t228 - t184 * t223;
t254 = qJD(5) * t308;
t264 = pkin(3) * t274;
t164 = t224 * t254 + t229 * t264;
t165 = -t224 * t264 + t229 * t254;
t81 = qJD(6) * t139 + t164 * t228 + t165 * t223;
t140 = t183 * t223 + t184 * t228;
t82 = -qJD(6) * t140 - t164 * t223 + t165 * t228;
t336 = t82 * mrSges(7,1) - t81 * mrSges(7,2);
t95 = t244 * t145;
t320 = -pkin(11) - pkin(10);
t205 = t320 * t224;
t207 = pkin(10) * t229 + t220;
t160 = t205 * t228 - t207 * t223;
t261 = qJD(5) * t320;
t196 = t224 * t261;
t197 = t229 * t261;
t115 = qJD(6) * t160 + t196 * t228 + t197 * t223;
t162 = t205 * t223 + t207 * t228;
t116 = -qJD(6) * t162 - t196 * t223 + t197 * t228;
t335 = t116 * mrSges(7,1) - t115 * mrSges(7,2);
t333 = t341 * t230;
t306 = mrSges(7,3) * t244;
t100 = t115 * t306;
t108 = mrSges(7,1) * t152 - mrSges(7,2) * t151;
t314 = pkin(5) * t229;
t215 = -pkin(4) - t314;
t104 = t215 * t108;
t307 = mrSges(7,3) * t152;
t121 = t162 * t307;
t155 = mrSges(7,1) * t244 + mrSges(7,2) * t188;
t218 = pkin(5) * t273;
t142 = t155 * t218;
t250 = mrSges(6,1) * t224 + mrSges(6,2) * t229;
t193 = t250 * qJD(5);
t315 = pkin(4) * t193;
t330 = -t100 + t104 - t121 + t142 - t315;
t107 = t140 * t307;
t265 = pkin(3) * t275;
t198 = t218 + t265;
t129 = t198 * t155;
t316 = pkin(3) * t230;
t213 = -pkin(4) - t316;
t170 = t213 * t193;
t200 = -mrSges(6,1) * t229 + mrSges(6,2) * t224;
t180 = t200 * t265;
t252 = mrSges(6,3) * t264;
t203 = t221 * t252;
t204 = t222 * t252;
t68 = t81 * t306;
t199 = t215 - t316;
t98 = t199 * t108;
t329 = -t107 + t129 + t170 + t180 + t203 + t204 - t68 + t98;
t328 = 2 * m(5);
t327 = 0.2e1 * m(6);
t326 = 2 * m(7);
t325 = 0.2e1 * t37;
t324 = 0.2e1 * t135;
t323 = 0.2e1 * t216;
t312 = t37 * t334;
t305 = Ifges(6,4) * t224;
t304 = Ifges(6,4) * t229;
t303 = Ifges(6,6) * t224;
t302 = pkin(3) * qJD(4);
t301 = pkin(5) * qJD(6);
t298 = t137 * mrSges(5,2);
t138 = t214 * t275 + (t226 * t274 + (t225 * t231 + t283) * qJD(3)) * pkin(2);
t297 = t138 * t334;
t291 = t137 * t221;
t290 = t137 * t222;
t289 = t145 * t224;
t288 = t145 * t229;
t271 = qJD(6) * t223;
t270 = qJD(6) * t228;
t269 = 0.2e1 * mrSges(7,3);
t268 = 0.2e1 * t232;
t30 = -t145 * t152 - t244 * t76;
t31 = -t188 * t76 + t331 * t95;
t267 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t77;
t266 = mrSges(7,3) * t301;
t263 = t228 * t151 * mrSges(7,3);
t258 = -t273 / 0.2e1;
t257 = -(2 * Ifges(5,4)) - t303;
t54 = -t224 * t92 + t229 * t90;
t253 = t341 * t137;
t173 = -pkin(2) * t284 + t214 * t230;
t171 = -pkin(4) - t173;
t251 = -t225 * mrSges(5,1) - t230 * mrSges(5,2);
t249 = Ifges(6,1) * t229 - t305;
t248 = -Ifges(6,2) * t224 + t304;
t247 = Ifges(6,5) * t224 + Ifges(6,6) * t229;
t38 = -pkin(5) * t245 - pkin(11) * t288 + t54;
t46 = -pkin(11) * t289 + t55;
t18 = -t223 * t46 + t228 * t38;
t19 = t223 * t38 + t228 * t46;
t101 = mrSges(6,2) * t245 - mrSges(6,3) * t289;
t102 = -mrSges(6,1) * t245 - mrSges(6,3) * t288;
t246 = t229 * t101 - t224 * t102;
t10 = -pkin(11) * t241 + t13;
t5 = -pkin(11) * t293 + pkin(5) * t77 + (-t89 + (pkin(11) * t145 - t90) * t224) * qJD(5) + t256;
t3 = qJD(6) * t18 + t10 * t228 + t223 * t5;
t4 = -qJD(6) * t19 - t10 * t223 + t228 * t5;
t243 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t267;
t242 = -t228 * t244 * t266 + (-pkin(5) * t307 + t188 * t266) * t223 + t340;
t239 = (-mrSges(4,1) * t226 - mrSges(4,2) * t231) * qJD(3) * pkin(2);
t109 = -Ifges(7,4) * t151 - Ifges(7,2) * t152;
t110 = -Ifges(7,1) * t151 - Ifges(7,4) * t152;
t156 = Ifges(7,4) * t188 - Ifges(7,2) * t244;
t157 = Ifges(7,1) * t188 - Ifges(7,4) * t244;
t194 = t248 * qJD(5);
t195 = t249 * qJD(5);
t201 = Ifges(6,2) * t229 + t305;
t202 = Ifges(6,1) * t224 + t304;
t238 = -t109 * t244 + t188 * t110 - t151 * t157 - t152 * t156 + t229 * t194 + t224 * t195 - t201 * t273 + t202 * t272;
t130 = t218 + t138;
t103 = t130 * t155;
t127 = t138 * t200;
t131 = mrSges(6,3) * t291;
t132 = mrSges(6,3) * t290;
t136 = t138 * mrSges(5,1);
t150 = t171 * t193;
t47 = t49 * t306;
t87 = t120 * t307;
t166 = t171 - t314;
t93 = t166 * t108;
t236 = t103 + t127 + t131 + t132 - t136 + t150 + t238 - t47 - t87 + t93 - t298;
t39 = mrSges(6,1) * t77 + mrSges(6,3) * t240;
t40 = -mrSges(6,2) * t77 - mrSges(6,3) * t241;
t235 = -t101 * t273 + m(6) * (-t272 * t54 - t273 * t55 + t339) + t229 * t40 - t224 * t39 - t102 * t272;
t21 = pkin(5) * t241 + t37;
t26 = -Ifges(6,4) * t240 - Ifges(6,2) * t241 + Ifges(6,6) * t77;
t27 = -Ifges(6,1) * t240 - Ifges(6,4) * t241 + Ifges(6,5) * t77;
t94 = t188 * t145;
t51 = -Ifges(7,4) * t95 - Ifges(7,2) * t94 - Ifges(7,6) * t245;
t52 = -Ifges(7,1) * t95 - Ifges(7,4) * t94 - Ifges(7,5) * t245;
t60 = pkin(5) * t289 - t334;
t69 = -Ifges(6,6) * t245 + t145 * t248;
t70 = -Ifges(6,5) * t245 + t145 * t249;
t8 = Ifges(7,4) * t30 + Ifges(7,2) * t31 + Ifges(7,6) * t77;
t9 = Ifges(7,1) * t30 + Ifges(7,4) * t31 + Ifges(7,5) * t77;
t234 = -t194 * t289 / 0.2e1 + t195 * t288 / 0.2e1 + t70 * t272 / 0.2e1 - t19 * t307 - t3 * t306 + ((-t224 * t55 - t229 * t54) * qJD(5) + t339) * mrSges(6,3) - t241 * t201 / 0.2e1 - (-Ifges(6,6) * t273 + t340) * t245 / 0.2e1 + t69 * t258 - t334 * t193 + (Ifges(7,5) * t188 - Ifges(7,6) * t244 + t247) * t77 / 0.2e1 - t244 * t8 / 0.2e1 + (t151 * t18 - t188 * t4) * mrSges(7,3) + (t200 - mrSges(5,1)) * t37 + (t145 * t258 + t293 / 0.2e1) * t202 + t229 * t26 / 0.2e1 + t224 * t27 / 0.2e1 + t188 * t9 / 0.2e1 - t151 * t52 / 0.2e1 - t152 * t51 / 0.2e1 + t21 * t155 + t31 * t156 / 0.2e1 + t30 * t157 / 0.2e1 + t60 * t108 - t94 * t109 / 0.2e1 - t95 * t110 / 0.2e1 - Ifges(5,6) * t77 + Ifges(5,5) * t76 - t36 * mrSges(5,2);
t233 = t118 * mrSges(4,1) - t117 * mrSges(4,2) + Ifges(4,5) * t153 - Ifges(4,6) * t154 + t234;
t182 = (-mrSges(7,1) * t223 - mrSges(7,2) * t228) * t301;
t99 = t250 * t145;
t66 = -mrSges(7,1) * t245 + mrSges(7,3) * t95;
t65 = mrSges(7,2) * t245 - mrSges(7,3) * t94;
t56 = mrSges(7,1) * t94 - mrSges(7,2) * t95;
t17 = -mrSges(7,2) * t77 + mrSges(7,3) * t31;
t16 = mrSges(7,1) * t77 - mrSges(7,3) * t30;
t11 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t1 = [(t18 * t4 + t19 * t3 + t21 * t60) * t326 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t232) * t268 + (m(4) * pkin(2) * t323 - 0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (-mrSges(4,1) * t187 + mrSges(4,2) * t189) - 0.2e1 * Ifges(3,4) * t227 + (Ifges(3,1) - Ifges(3,2)) * t268) * t227) * qJD(2) - t69 * t295 + 0.2e1 * m(4) * (t117 * t163 + t118 * t161) + t99 * t325 + t70 * t293 - 0.2e1 * t187 * Ifges(4,2) * t154 + 0.2e1 * t153 * t189 * Ifges(4,1) + 0.2e1 * (-t334 * t76 - t77 * t92) * mrSges(5,3) - 0.2e1 * t334 * t32 + (mrSges(4,1) * t154 + mrSges(4,2) * t153) * t323 + 0.2e1 * (t153 * t187 - t154 * t189) * Ifges(4,4) + 0.2e1 * (t117 * t187 - t118 * t189 - t153 * t161 - t154 * t163) * mrSges(4,3) + (t13 * t55 + t14 * t54 - t312) * t327 + (t135 * t167 + t36 * t92 - t312) * t328 - (mrSges(5,1) * t324 - 0.2e1 * t36 * mrSges(5,3) + t257 * t76 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t77 + t267 + t338) * t245 + (mrSges(5,2) * t324 + mrSges(5,3) * t325 + 0.2e1 * Ifges(5,1) * t76 - t224 * t26 + t229 * t27 + (Ifges(6,5) * t229 + t257) * t77 + (-t224 * t70 - t229 * t69 + t245 * t247) * qJD(5)) * t145 + 0.2e1 * t167 * (mrSges(5,1) * t77 + mrSges(5,2) * t76) - t94 * t8 - t95 * t9 + 0.2e1 * t13 * t101 + 0.2e1 * t14 * t102 + 0.2e1 * t3 * t65 + 0.2e1 * t4 * t66 + 0.2e1 * t54 * t39 + 0.2e1 * t55 * t40 + 0.2e1 * t21 * t56 + 0.2e1 * t60 * t11 + t31 * t51 + t30 * t52 + 0.2e1 * t18 * t16 + 0.2e1 * t19 * t17 + t77 * (-Ifges(7,5) * t95 - Ifges(7,6) * t94); (Ifges(3,5) * t232 - Ifges(3,6) * t227 + (-mrSges(3,1) * t232 + mrSges(3,2) * t227) * pkin(7)) * qJD(2) + t246 * t137 + t233 + m(6) * (t171 * t37 + t281 * t55 - t285 * t54 - t297) + m(5) * (t137 * t92 - t173 * t37 + t174 * t36 - t297) + m(7) * (t119 * t4 + t120 * t3 + t130 * t60 + t166 * t21 + t18 * t50 + t19 * t49) + t166 * t11 + t171 * t32 + t138 * t99 + t130 * t56 + t119 * t16 + t120 * t17 + t49 * t65 + t50 * t66 + t235 * t172 + (m(4) * (t117 * t226 + t118 * t231 + (-t161 * t226 + t163 * t231) * qJD(3)) + (-t231 * t153 - t226 * t154 + (t187 * t231 + t189 * t226) * qJD(3)) * mrSges(4,3)) * pkin(2) + (t137 * t245 + t138 * t145 - t173 * t76 - t174 * t77) * mrSges(5,3); 0.2e1 * t150 - 0.2e1 * t136 + 0.2e1 * t239 + 0.2e1 * t132 + 0.2e1 * t131 + 0.2e1 * t127 - 0.2e1 * t87 + 0.2e1 * t103 + (t138 * t171 + t172 * t253) * t327 + (t137 * t174 - t138 * t173) * t328 + (t119 * t50 + t120 * t49 + t130 * t166) * t326 + 0.2e1 * t93 + t238 + (t119 * t151 - t50 * t188) * t269 - 0.2e1 * t298 - 0.2e1 * t47; t233 + m(7) * (t139 * t4 + t140 * t3 + t18 * t82 + t19 * t81 + t198 * t60 + t199 * t21) + t199 * t11 + t198 * t56 + t139 * t16 + t140 * t17 + t81 * t65 + t82 * t66 + t235 * t212 + (m(5) * (t225 * t36 - t230 * t37) + (-t225 * t77 - t230 * t76) * mrSges(5,3) + ((t245 * mrSges(5,3) + t246 + m(6) * (-t224 * t54 + t229 * t55) + m(5) * t92) * t230 + (t145 * mrSges(5,3) + t99 - (m(6) + m(5)) * t334) * t225) * qJD(4)) * pkin(3) + t342 * t213; ((-t50 - t82) * t188 - (-t119 - t139) * t151) * mrSges(7,3) + m(6) * (t138 * t213 + (t290 + t291) * t212) + m(7) * (t119 * t82 + t120 * t81 + t130 * t199 + t139 * t50 + t140 * t49 + t166 * t198) + t239 + (m(5) * (t137 * t225 - t138 * t230) + (m(6) * (t171 * t225 + t172 * t333) + m(5) * (-t173 * t225 + t174 * t230) + t251) * qJD(4)) * pkin(3) + t236 + t329; 0.2e1 * t170 - 0.2e1 * t107 + 0.2e1 * (m(6) * (t212 * t333 + t213 * t225) + t251) * t302 + 0.2e1 * t129 + 0.2e1 * t98 + t238 + 0.2e1 * t203 + 0.2e1 * t204 + (t139 * t151 - t188 * t82) * t269 - 0.2e1 * t68 + 0.2e1 * t180 + (t139 * t82 + t140 * t81 + t198 * t199) * t326; t234 + t56 * t218 + t215 * t11 + t160 * t16 + t162 * t17 + t115 * t65 + t116 * t66 + t235 * pkin(10) + m(7) * (t115 * t19 + t116 * t18 + t160 * t4 + t162 * t3 + t21 * t215 + t218 * t60) - t342 * pkin(4); m(6) * (-pkin(4) * t138 + pkin(10) * t253) + ((-t116 - t50) * t188 - (-t119 - t160) * t151) * mrSges(7,3) + m(7) * (t115 * t120 + t116 * t119 + t130 * t215 + t160 * t50 + t162 * t49 + t166 * t218) + t236 + t330; m(7) * (t115 * t140 + t116 * t139 + t160 * t82 + t162 * t81 + t198 * t215 + t199 * t218) + (m(6) * (-pkin(4) * t225 + pkin(10) * t333) + t251) * t302 + ((-t116 - t82) * t188 - (-t139 - t160) * t151) * mrSges(7,3) + t238 + t329 + t330; 0.2e1 * t142 + 0.2e1 * t104 - 0.2e1 * t315 + (t115 * t162 + t116 * t160 + t215 * t218) * t326 - 0.2e1 * t100 - 0.2e1 * t121 + (-t116 * t188 + t151 * t160) * t269 + t238; -Ifges(6,5) * t260 + t14 * mrSges(6,1) - t13 * mrSges(6,2) - t241 * Ifges(6,6) + (m(7) * (-t18 * t271 + t19 * t270 + t223 * t3 + t228 * t4) + t65 * t270 + t223 * t17 - t66 * t271 + t228 * t16) * pkin(5) + t243 + t338; -t250 * t137 + (t172 * t200 - t303) * qJD(5) + (m(7) * (-t119 * t271 + t120 * t270 + t223 * t49 + t228 * t50) + t263) * pkin(5) + t242 + t337; -t250 * t264 + (t200 * t212 - t303) * qJD(5) + (m(7) * (-t139 * t271 + t140 * t270 + t223 * t81 + t228 * t82) + t263) * pkin(5) + t242 + t336; (pkin(10) * t200 - t303) * qJD(5) + (m(7) * (t115 * t223 + t116 * t228 - t160 * t271 + t162 * t270) + t263) * pkin(5) + t242 + t335; 0.2e1 * t182; t243; t278 + t337; t278 + t336; t278 + t335; t182; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
