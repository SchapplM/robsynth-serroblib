% Calculate time derivative of joint inertia matrix for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:40
% EndTime: 2019-03-09 10:55:54
% DurationCPUTime: 6.07s
% Computational Cost: add. (11735->614), mult. (30786->926), div. (0->0), fcn. (31601->12), ass. (0->247)
t246 = sin(pkin(6));
t317 = 0.2e1 * t246;
t244 = sin(pkin(12));
t247 = cos(pkin(12));
t250 = sin(qJ(6));
t253 = cos(qJ(6));
t256 = t244 * t250 - t247 * t253;
t297 = -t256 / 0.2e1;
t217 = t244 * t253 + t247 * t250;
t296 = t217 / 0.2e1;
t316 = t244 / 0.2e1;
t293 = t247 / 0.2e1;
t252 = sin(qJ(2));
t255 = cos(qJ(2));
t184 = (-qJD(3) * t252 + (pkin(2) * t252 - qJ(3) * t255) * qJD(2)) * t246;
t270 = qJD(2) * t252;
t265 = t246 * t270;
t249 = cos(pkin(6));
t291 = pkin(1) * t255;
t267 = t249 * t291;
t198 = -pkin(8) * t265 + qJD(2) * t267;
t189 = qJD(3) * t249 + t198;
t245 = sin(pkin(11));
t248 = cos(pkin(11));
t124 = t248 * t184 - t245 * t189;
t271 = qJD(2) * t246;
t274 = t248 * t255;
t100 = (pkin(3) * t252 - pkin(9) * t274) * t271 + t124;
t275 = t246 * t255;
t213 = t249 * t252 * pkin(1) + pkin(8) * t275;
t194 = qJ(3) * t249 + t213;
t195 = (-pkin(2) * t255 - qJ(3) * t252 - pkin(1)) * t246;
t137 = -t245 * t194 + t248 * t195;
t276 = t246 * t252;
t207 = t245 * t249 + t248 * t276;
t102 = -pkin(3) * t275 - t207 * pkin(9) + t137;
t125 = t245 * t184 + t248 * t189;
t264 = t255 * t271;
t261 = t245 * t264;
t116 = -pkin(9) * t261 + t125;
t138 = t248 * t194 + t245 * t195;
t206 = -t245 * t276 + t248 * t249;
t118 = pkin(9) * t206 + t138;
t251 = sin(qJ(4));
t254 = cos(qJ(4));
t268 = qJD(4) * t254;
t269 = qJD(4) * t251;
t24 = t100 * t254 - t102 * t269 - t251 * t116 - t118 * t268;
t22 = -pkin(4) * t265 - t24;
t216 = t245 * t251 - t254 * t248;
t257 = t254 * t206 - t207 * t251;
t112 = qJD(4) * t257 - t216 * t264;
t92 = -t112 * t244 + t247 * t265;
t93 = t112 * t247 + t244 * t265;
t49 = -t92 * mrSges(6,1) + t93 * mrSges(6,2);
t315 = -m(6) * t22 - t49;
t290 = pkin(9) + qJ(3);
t224 = t290 * t245;
t226 = t290 * t248;
t314 = -t254 * t224 - t226 * t251;
t208 = t256 * qJD(6);
t313 = 2 * m(4);
t312 = 2 * m(5);
t311 = 0.2e1 * m(6);
t310 = 2 * m(7);
t309 = -2 * mrSges(3,3);
t308 = t92 / 0.2e1;
t210 = t216 * qJD(4);
t218 = t245 * t254 + t248 * t251;
t211 = t218 * qJD(4);
t285 = Ifges(6,4) * t247;
t259 = -Ifges(6,2) * t244 + t285;
t120 = Ifges(6,6) * t211 - t210 * t259;
t306 = t120 / 0.2e1;
t286 = Ifges(6,4) * t244;
t260 = Ifges(6,1) * t247 - t286;
t121 = Ifges(6,5) * t211 - t210 * t260;
t305 = t121 / 0.2e1;
t209 = t217 * qJD(6);
t157 = -Ifges(7,5) * t208 - Ifges(7,6) * t209;
t304 = t157 / 0.2e1;
t158 = -Ifges(7,4) * t208 - Ifges(7,2) * t209;
t303 = t158 / 0.2e1;
t160 = -Ifges(7,1) * t208 - Ifges(7,4) * t209;
t302 = t160 / 0.2e1;
t170 = Ifges(7,4) * t217 - Ifges(7,2) * t256;
t301 = t170 / 0.2e1;
t172 = Ifges(7,1) * t217 - Ifges(7,4) * t256;
t300 = t172 / 0.2e1;
t299 = -t208 / 0.2e1;
t298 = -t209 / 0.2e1;
t295 = Ifges(6,1) * t316 + t285 / 0.2e1;
t294 = -t244 / 0.2e1;
t292 = t248 / 0.2e1;
t289 = pkin(10) + qJ(5);
t23 = t251 * t100 + t102 * t268 + t254 * t116 - t118 * t269;
t21 = (qJ(5) * t270 - qJD(5) * t255) * t246 + t23;
t146 = t206 * t251 + t207 * t254;
t113 = qJD(4) * t146 + t218 * t264;
t199 = t213 * qJD(2);
t176 = pkin(3) * t261 + t199;
t41 = t113 * pkin(4) - t112 * qJ(5) - t146 * qJD(5) + t176;
t12 = t247 * t21 + t244 * t41;
t63 = t251 * t102 + t254 * t118;
t58 = -qJ(5) * t275 + t63;
t233 = pkin(8) * t276;
t197 = t233 + (-pkin(2) - t291) * t249;
t153 = -t206 * pkin(3) + t197;
t69 = -pkin(4) * t257 - t146 * qJ(5) + t153;
t31 = t244 * t69 + t247 * t58;
t288 = Ifges(4,4) * t245;
t287 = Ifges(4,4) * t248;
t284 = t176 * mrSges(5,1);
t283 = t176 * mrSges(5,2);
t180 = -t224 * t251 + t226 * t254;
t142 = qJD(3) * t218 + qJD(4) * t180;
t282 = t142 * t314;
t281 = t210 * t244;
t280 = t210 * t247;
t279 = t218 * t244;
t278 = t218 * t247;
t133 = pkin(4) * t211 + qJ(5) * t210 - qJD(5) * t218;
t141 = -t216 * qJD(3) + qJD(4) * t314;
t76 = t244 * t133 + t247 * t141;
t239 = -pkin(3) * t248 - pkin(2);
t162 = pkin(4) * t216 - qJ(5) * t218 + t239;
t104 = t244 * t162 + t247 * t180;
t144 = -mrSges(6,1) * t281 - mrSges(6,2) * t280;
t273 = -Ifges(5,5) * t210 - Ifges(5,6) * t211;
t185 = t248 * mrSges(4,2) * t264 + mrSges(4,1) * t261;
t128 = -t244 * t146 - t247 * t275;
t129 = t247 * t146 - t244 * t275;
t70 = t128 * t253 - t129 * t250;
t28 = qJD(6) * t70 + t250 * t92 + t253 * t93;
t71 = t128 * t250 + t129 * t253;
t29 = -qJD(6) * t71 - t250 * t93 + t253 * t92;
t6 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t113;
t89 = -t209 * t218 + t210 * t256;
t90 = t208 * t218 + t210 * t217;
t45 = Ifges(7,5) * t89 + Ifges(7,6) * t90 + Ifges(7,3) * t211;
t266 = Ifges(5,5) * t112 - Ifges(5,6) * t113 + Ifges(5,3) * t265;
t263 = Ifges(6,5) * t316 + Ifges(7,5) * t296 + Ifges(6,6) * t293 + Ifges(7,6) * t297;
t10 = -t29 * mrSges(7,1) + t28 * mrSges(7,2);
t48 = -t90 * mrSges(7,1) + t89 * mrSges(7,2);
t11 = -t21 * t244 + t247 * t41;
t30 = -t244 * t58 + t247 * t69;
t64 = t113 * mrSges(5,1) + t112 * mrSges(5,2);
t156 = t211 * mrSges(5,1) - t210 * mrSges(5,2);
t62 = t102 * t254 - t251 * t118;
t75 = t247 * t133 - t141 * t244;
t103 = t247 * t162 - t180 * t244;
t59 = pkin(4) * t275 - t62;
t258 = Ifges(6,5) * t247 - Ifges(6,6) * t244;
t16 = -pkin(5) * t257 - pkin(10) * t129 + t30;
t19 = pkin(10) * t128 + t31;
t3 = t16 * t253 - t19 * t250;
t4 = t16 * t250 + t19 * t253;
t77 = pkin(5) * t216 - pkin(10) * t278 + t103;
t85 = -pkin(10) * t279 + t104;
t43 = -t250 * t85 + t253 * t77;
t44 = t250 * t77 + t253 * t85;
t222 = t289 * t244;
t225 = t289 * t247;
t177 = -t222 * t253 - t225 * t250;
t179 = -t222 * t250 + t225 * t253;
t238 = -pkin(5) * t247 - pkin(4);
t231 = Ifges(3,5) * t264;
t228 = Ifges(6,2) * t247 + t286;
t223 = -mrSges(6,1) * t247 + mrSges(6,2) * t244;
t212 = -t233 + t267;
t191 = (mrSges(4,1) * t252 - mrSges(4,3) * t274) * t271;
t190 = (-mrSges(4,3) * t245 * t255 - mrSges(4,2) * t252) * t271;
t183 = -mrSges(4,1) * t275 - t207 * mrSges(4,3);
t182 = mrSges(4,2) * t275 + t206 * mrSges(4,3);
t173 = Ifges(5,1) * t218 - Ifges(5,4) * t216;
t171 = Ifges(5,4) * t218 - Ifges(5,2) * t216;
t168 = mrSges(7,1) * t256 + mrSges(7,2) * t217;
t166 = (t252 * Ifges(4,5) + (t248 * Ifges(4,1) - t288) * t255) * t271;
t165 = (t252 * Ifges(4,6) + (-t245 * Ifges(4,2) + t287) * t255) * t271;
t164 = mrSges(6,1) * t216 - mrSges(6,3) * t278;
t163 = -mrSges(6,2) * t216 - mrSges(6,3) * t279;
t161 = -Ifges(5,1) * t210 - Ifges(5,4) * t211;
t159 = -Ifges(5,4) * t210 - Ifges(5,2) * t211;
t155 = mrSges(7,1) * t209 - mrSges(7,2) * t208;
t154 = (mrSges(6,1) * t244 + mrSges(6,2) * t247) * t218;
t152 = mrSges(6,1) * t211 + mrSges(6,3) * t280;
t151 = -mrSges(6,2) * t211 + mrSges(6,3) * t281;
t148 = t256 * t218;
t147 = t217 * t218;
t143 = pkin(5) * t279 - t314;
t140 = -qJD(5) * t217 - qJD(6) * t179;
t139 = -qJD(5) * t256 + qJD(6) * t177;
t135 = -mrSges(5,1) * t275 - t146 * mrSges(5,3);
t134 = mrSges(5,2) * t275 + mrSges(5,3) * t257;
t132 = Ifges(6,5) * t216 + t218 * t260;
t131 = Ifges(6,6) * t216 + t218 * t259;
t130 = Ifges(6,3) * t216 + t218 * t258;
t123 = mrSges(7,1) * t216 + mrSges(7,3) * t148;
t122 = -mrSges(7,2) * t216 - mrSges(7,3) * t147;
t119 = Ifges(6,3) * t211 - t210 * t258;
t117 = -pkin(5) * t281 + t142;
t96 = -mrSges(5,2) * t265 - mrSges(5,3) * t113;
t95 = mrSges(5,1) * t265 - mrSges(5,3) * t112;
t94 = mrSges(7,1) * t147 - mrSges(7,2) * t148;
t84 = Ifges(5,1) * t146 + Ifges(5,4) * t257 - Ifges(5,5) * t275;
t83 = Ifges(5,4) * t146 + Ifges(5,2) * t257 - Ifges(5,6) * t275;
t82 = -mrSges(6,1) * t257 - mrSges(6,3) * t129;
t81 = mrSges(6,2) * t257 + mrSges(6,3) * t128;
t80 = -Ifges(7,1) * t148 - Ifges(7,4) * t147 + Ifges(7,5) * t216;
t79 = -Ifges(7,4) * t148 - Ifges(7,2) * t147 + Ifges(7,6) * t216;
t78 = -Ifges(7,5) * t148 - Ifges(7,6) * t147 + Ifges(7,3) * t216;
t74 = -mrSges(7,2) * t211 + mrSges(7,3) * t90;
t73 = mrSges(7,1) * t211 - mrSges(7,3) * t89;
t72 = -mrSges(6,1) * t128 + mrSges(6,2) * t129;
t66 = pkin(10) * t281 + t76;
t65 = pkin(5) * t211 + pkin(10) * t280 + t75;
t61 = Ifges(5,1) * t112 - Ifges(5,4) * t113 + Ifges(5,5) * t265;
t60 = Ifges(5,4) * t112 - Ifges(5,2) * t113 + Ifges(5,6) * t265;
t57 = Ifges(6,1) * t129 + Ifges(6,4) * t128 - Ifges(6,5) * t257;
t56 = Ifges(6,4) * t129 + Ifges(6,2) * t128 - Ifges(6,6) * t257;
t55 = Ifges(6,5) * t129 + Ifges(6,6) * t128 - Ifges(6,3) * t257;
t53 = mrSges(6,1) * t113 - mrSges(6,3) * t93;
t52 = -mrSges(6,2) * t113 + mrSges(6,3) * t92;
t51 = -mrSges(7,1) * t257 - mrSges(7,3) * t71;
t50 = mrSges(7,2) * t257 + mrSges(7,3) * t70;
t47 = Ifges(7,1) * t89 + Ifges(7,4) * t90 + Ifges(7,5) * t211;
t46 = Ifges(7,4) * t89 + Ifges(7,2) * t90 + Ifges(7,6) * t211;
t42 = -pkin(5) * t128 + t59;
t38 = -mrSges(7,1) * t70 + mrSges(7,2) * t71;
t37 = Ifges(6,1) * t93 + Ifges(6,4) * t92 + Ifges(6,5) * t113;
t36 = Ifges(6,4) * t93 + Ifges(6,2) * t92 + Ifges(6,6) * t113;
t35 = Ifges(6,5) * t93 + Ifges(6,6) * t92 + Ifges(6,3) * t113;
t34 = Ifges(7,1) * t71 + Ifges(7,4) * t70 - Ifges(7,5) * t257;
t33 = Ifges(7,4) * t71 + Ifges(7,2) * t70 - Ifges(7,6) * t257;
t32 = Ifges(7,5) * t71 + Ifges(7,6) * t70 - Ifges(7,3) * t257;
t18 = -mrSges(7,2) * t113 + mrSges(7,3) * t29;
t17 = mrSges(7,1) * t113 - mrSges(7,3) * t28;
t15 = -pkin(5) * t92 + t22;
t14 = -qJD(6) * t44 - t250 * t66 + t253 * t65;
t13 = qJD(6) * t43 + t250 * t65 + t253 * t66;
t9 = pkin(10) * t92 + t12;
t8 = Ifges(7,1) * t28 + Ifges(7,4) * t29 + Ifges(7,5) * t113;
t7 = Ifges(7,4) * t28 + Ifges(7,2) * t29 + Ifges(7,6) * t113;
t5 = pkin(5) * t113 - pkin(10) * t93 + t11;
t2 = -qJD(6) * t4 - t250 * t9 + t253 * t5;
t1 = qJD(6) * t3 + t250 * t5 + t253 * t9;
t20 = [-(t35 + t6 - t60 + 0.2e1 * t284) * t257 + (t1 * t4 + t15 * t42 + t2 * t3) * t310 + (t11 * t30 + t12 * t31 + t22 * t59) * t311 + (t153 * t176 + t23 * t63 + t24 * t62) * t312 + (t124 * t137 + t125 * t138 + t197 * t199) * t313 + 0.2e1 * m(3) * (t198 * t213 - t199 * t212) + (t55 + t32 - t83) * t113 + (t61 + 0.2e1 * t283) * t146 - 0.2e1 * t199 * (mrSges(3,1) * t249 - mrSges(3,3) * t276) - t266 * t275 + 0.2e1 * t198 * (-t249 * mrSges(3,2) + mrSges(3,3) * t275) + t249 * t231 + t206 * t165 + 0.2e1 * t199 * (-mrSges(4,1) * t206 + mrSges(4,2) * t207) + t207 * t166 + 0.2e1 * t197 * t185 + 0.2e1 * t125 * t182 + 0.2e1 * t124 * t183 + 0.2e1 * t138 * t190 + 0.2e1 * t137 * t191 + 0.2e1 * t153 * t64 + t129 * t37 + 0.2e1 * t23 * t134 + 0.2e1 * t24 * t135 + t128 * t36 + t112 * t84 + t92 * t56 + t93 * t57 + 0.2e1 * t62 * t95 + 0.2e1 * t63 * t96 + 0.2e1 * t11 * t82 + 0.2e1 * t12 * t81 + t70 * t7 + t71 * t8 + 0.2e1 * t22 * t72 + 0.2e1 * t59 * t49 + 0.2e1 * t1 * t50 + 0.2e1 * t2 * t51 + 0.2e1 * t31 * t52 + 0.2e1 * t30 * t53 + 0.2e1 * t15 * t38 + 0.2e1 * t42 * t10 + t29 * t33 + t28 * t34 + 0.2e1 * t3 * t17 + 0.2e1 * t4 * t18 + ((t213 * t309 + Ifges(4,5) * t207 + Ifges(5,5) * t146 - 0.2e1 * Ifges(3,6) * t249 + Ifges(4,6) * t206 + Ifges(5,6) * t257 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t252) * t317) * t252 + (Ifges(3,5) * t249 - t245 * (Ifges(4,4) * t207 + Ifges(4,2) * t206) + t212 * t309 + t248 * (Ifges(4,1) * t207 + Ifges(4,4) * t206) + (-pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t248 + Ifges(4,6) * t245 + Ifges(3,4)) * t255) * t317 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3)) * t276) * t255) * t271; -(-t159 / 0.2e1 + t119 / 0.2e1 + t45 / 0.2e1) * t257 + (t199 * mrSges(4,2) + t166 / 0.2e1 - qJD(3) * t183 - qJ(3) * t191 - t124 * mrSges(4,3)) * t245 + t231 + (t72 - t135) * t142 + (t210 * t62 - t211 * t63 - t216 * t23 - t218 * t24) * mrSges(5,3) + m(5) * (t141 * t63 - t142 * t62 + t176 * t239 + t180 * t23 + t24 * t314) + m(6) * (t103 * t11 + t104 * t12 + t142 * t59 - t22 * t314 + t30 * t75 + t31 * t76) - (t49 - t95) * t314 + t131 * t308 + t129 * t305 + t128 * t306 + m(4) * (-pkin(2) * t199 + (-t137 * t245 + t138 * t248) * qJD(3) + (-t124 * t245 + t125 * t248) * qJ(3)) + m(7) * (t1 * t44 + t117 * t42 + t13 * t4 + t14 * t3 + t143 * t15 + t2 * t43) + (-t199 * mrSges(4,1) + t165 / 0.2e1 + qJ(3) * t190 + t125 * mrSges(4,3) + qJD(3) * t182) * t248 + (-t255 * t273 / 0.2e1 + ((Ifges(5,5) * t218 / 0.2e1 - Ifges(5,6) * t216 / 0.2e1 - Ifges(3,6) + Ifges(4,5) * t245 / 0.2e1 + Ifges(4,6) * t292) * t252 + (-t245 * (Ifges(4,2) * t248 + t288) / 0.2e1 + (Ifges(4,1) * t245 + t287) * t292) * t255) * qJD(2)) * t246 - (t84 / 0.2e1 + t57 * t293 + t56 * t294) * t210 + (t283 + t61 / 0.2e1 + t37 * t293 + t36 * t294) * t218 + (t284 + t35 / 0.2e1 + t6 / 0.2e1 - t60 / 0.2e1) * t216 + t239 * t64 - t198 * mrSges(3,2) - t199 * mrSges(3,1) + t180 * t96 - pkin(2) * t185 + t12 * t163 + t11 * t164 + t112 * t173 / 0.2e1 + t31 * t151 + t30 * t152 + t22 * t154 + t153 * t156 + t146 * t161 / 0.2e1 - t147 * t7 / 0.2e1 - t148 * t8 / 0.2e1 + t141 * t134 + t143 * t10 + t59 * t144 + t93 * t132 / 0.2e1 + t1 * t122 + t2 * t123 + t117 * t38 + t103 * t53 + t104 * t52 + t15 * t94 + t75 * t82 + t89 * t34 / 0.2e1 + t90 * t33 / 0.2e1 + t4 * t74 + (-t171 / 0.2e1 + t130 / 0.2e1 + t78 / 0.2e1) * t113 + t29 * t79 / 0.2e1 + t28 * t80 / 0.2e1 + t76 * t81 + t70 * t46 / 0.2e1 + t71 * t47 / 0.2e1 + t3 * t73 + t42 * t48 + t13 * t50 + t14 * t51 + t43 * t17 + t44 * t18 + (t55 / 0.2e1 + t32 / 0.2e1 - t83 / 0.2e1) * t211; (t103 * t75 + t104 * t76 - t282) * t311 + (t141 * t180 - t282) * t312 + 0.2e1 * (-t141 * t216 + t142 * t218 - t180 * t211 + t210 * t314) * mrSges(5,3) - 0.2e1 * t314 * t144 + (t117 * t143 + t13 * t44 + t14 * t43) * t310 + (qJ(3) * t313 + 0.2e1 * mrSges(4,3)) * (t245 ^ 2 + t248 ^ 2) * qJD(3) + (-t120 * t244 + t121 * t247 + t161) * t218 + (t119 + t45 - t159) * t216 + (t130 + t78 - t171) * t211 - (-t131 * t244 + t132 * t247 + t173) * t210 + 0.2e1 * t239 * t156 + 0.2e1 * t76 * t163 + 0.2e1 * t75 * t164 + 0.2e1 * t104 * t151 + 0.2e1 * t103 * t152 + 0.2e1 * t142 * t154 - t147 * t46 - t148 * t47 + 0.2e1 * t143 * t48 + 0.2e1 * t13 * t122 + 0.2e1 * t14 * t123 + 0.2e1 * t117 * t94 + t89 * t80 + t90 * t79 + 0.2e1 * t44 * t74 + 0.2e1 * t43 * t73; -t256 * t17 + t217 * t18 - t208 * t50 - t209 * t51 + t244 * t52 + t247 * t53 + m(7) * (t1 * t217 - t2 * t256 - t208 * t4 - t209 * t3) + m(6) * (t11 * t247 + t12 * t244) + m(5) * t176 + m(4) * t199 + t64 + t185; -t208 * t122 - t209 * t123 + t244 * t151 + t247 * t152 - t256 * t73 + t217 * t74 + m(7) * (t13 * t217 - t14 * t256 - t208 * t44 - t209 * t43) + m(6) * (t244 * t76 + t247 * t75) + t156; (-t208 * t217 + t209 * t256) * t310; -t257 * t304 + (-t1 * t256 - t2 * t217 + t208 * t3 - t209 * t4) * mrSges(7,3) + t228 * t308 + t93 * t295 + t8 * t296 + t7 * t297 + t33 * t298 + t34 * t299 + t28 * t300 + t29 * t301 + t71 * t302 + t70 * t303 + (t37 / 0.2e1 + m(6) * (-qJ(5) * t11 - qJD(5) * t30) - qJD(5) * t82 - t11 * mrSges(6,3) - qJ(5) * t53) * t244 + (t36 / 0.2e1 + m(6) * (qJ(5) * t12 + qJD(5) * t31) + qJD(5) * t81 + t12 * mrSges(6,3) + qJ(5) * t52) * t247 + t266 + t263 * t113 + t238 * t10 + t22 * t223 + t179 * t18 + t15 * t168 + t177 * t17 + t42 * t155 + t140 * t51 + t139 * t50 - t23 * mrSges(5,2) + t24 * mrSges(5,1) + t315 * pkin(4) + m(7) * (t1 * t179 + t139 * t4 + t140 * t3 + t15 * t238 + t177 * t2); (-t13 * t256 - t14 * t217 + t208 * t43 - t209 * t44) * mrSges(7,3) + (t305 + m(6) * (-qJ(5) * t75 - qJD(5) * t103) + t210 * t228 / 0.2e1 - qJ(5) * t152 - qJD(5) * t164 - t75 * mrSges(6,3)) * t244 + (t306 + m(6) * (qJ(5) * t76 + qJD(5) * t104) - t210 * t295 + qJ(5) * t151 + qJD(5) * t163 + t76 * mrSges(6,3)) * t247 + t47 * t296 + t46 * t297 + t79 * t298 + t80 * t299 + t89 * t300 + t90 * t301 - t148 * t302 - t147 * t303 + t216 * t304 + (-m(6) * pkin(4) - mrSges(5,1) + t223) * t142 + m(7) * (t117 * t238 + t13 * t179 + t139 * t44 + t14 * t177 + t140 * t43) + t263 * t211 + t238 * t48 + t179 * t74 + t117 * t168 + t177 * t73 + t143 * t155 + t140 * t123 - t141 * mrSges(5,2) - pkin(4) * t144 + t139 * t122 + t273; m(7) * (t139 * t217 - t140 * t256 - t177 * t209 - t179 * t208); (t139 * t179 + t140 * t177) * t310 - t208 * t172 + t217 * t160 - t209 * t170 - t256 * t158 + 0.2e1 * t238 * t155 + 0.2e1 * (-t139 * t256 - t140 * t217 + t177 * t208 - t179 * t209) * mrSges(7,3) + (qJ(5) * t311 + 0.2e1 * mrSges(6,3)) * qJD(5) * (t244 ^ 2 + t247 ^ 2); m(7) * t15 + t10 - t315; m(6) * t142 + m(7) * t117 + t144 + t48; 0; t155; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t6; mrSges(7,1) * t14 - mrSges(7,2) * t13 + t45; -t155; mrSges(7,1) * t140 - mrSges(7,2) * t139 + t157; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
