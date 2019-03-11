% Calculate time derivative of joint inertia matrix for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:27
% EndTime: 2019-03-09 18:15:36
% DurationCPUTime: 4.38s
% Computational Cost: add. (13052->504), mult. (28642->723), div. (0->0), fcn. (28412->10), ass. (0->214)
t234 = sin(pkin(11));
t235 = cos(pkin(11));
t320 = -mrSges(5,1) * t235 + mrSges(5,2) * t234 - mrSges(4,1);
t237 = sin(qJ(5));
t241 = cos(qJ(5));
t271 = t235 * t241;
t248 = t234 * t237 - t271;
t200 = t248 * qJD(5);
t213 = t234 * t241 + t235 * t237;
t201 = t213 * qJD(5);
t236 = sin(qJ(6));
t240 = cos(qJ(6));
t170 = -t213 * t236 - t240 * t248;
t117 = qJD(6) * t170 - t200 * t240 - t201 * t236;
t171 = t213 * t240 - t236 * t248;
t118 = -qJD(6) * t171 + t200 * t236 - t201 * t240;
t270 = Ifges(7,5) * t117 + Ifges(7,6) * t118;
t319 = -Ifges(6,5) * t200 - Ifges(6,6) * t201 + t270;
t119 = -mrSges(7,1) * t170 + mrSges(7,2) * t171;
t318 = 0.2e1 * t119;
t238 = sin(qJ(3));
t226 = pkin(2) * t238 + qJ(4);
t203 = (-pkin(9) - t226) * t234;
t231 = t235 * pkin(9);
t272 = t226 * t235;
t204 = t231 + t272;
t242 = cos(qJ(3));
t288 = pkin(2) * qJD(3);
t224 = t242 * t288 + qJD(4);
t266 = qJD(5) * t241;
t273 = t224 * t234;
t124 = t203 * t266 + t224 * t271 + (-qJD(5) * t204 - t273) * t237;
t198 = t201 * pkin(10);
t102 = -t198 + t124;
t161 = t237 * t203 + t241 * t204;
t125 = -qJD(5) * t161 - t213 * t224;
t297 = pkin(10) * t200;
t103 = t125 + t297;
t160 = t241 * t203 - t204 * t237;
t296 = pkin(10) * t213;
t139 = t160 - t296;
t206 = t248 * pkin(10);
t140 = -t206 + t161;
t81 = t139 * t240 - t140 * t236;
t28 = qJD(6) * t81 + t102 * t240 + t103 * t236;
t82 = t139 * t236 + t140 * t240;
t29 = -qJD(6) * t82 - t102 * t236 + t103 * t240;
t317 = t29 * mrSges(7,1) - t28 * mrSges(7,2);
t219 = (-pkin(9) - qJ(4)) * t234;
t221 = qJ(4) * t235 + t231;
t146 = t219 * t266 + qJD(4) * t271 + (-qJD(4) * t234 - qJD(5) * t221) * t237;
t130 = -t198 + t146;
t185 = t237 * t219 + t241 * t221;
t147 = -t213 * qJD(4) - qJD(5) * t185;
t131 = t147 + t297;
t184 = t241 * t219 - t221 * t237;
t151 = t184 - t296;
t152 = -t206 + t185;
t97 = t151 * t240 - t152 * t236;
t44 = qJD(6) * t97 + t130 * t240 + t131 * t236;
t98 = t151 * t236 + t152 * t240;
t45 = -qJD(6) * t98 - t130 * t236 + t131 * t240;
t316 = t45 * mrSges(7,1) - t44 * mrSges(7,2);
t239 = sin(qJ(2));
t302 = -pkin(8) - pkin(7);
t222 = t302 * t239;
t243 = cos(qJ(2));
t223 = t302 * t243;
t315 = t242 * t222 + t223 * t238;
t314 = qJD(2) + qJD(3);
t215 = t238 * t243 + t242 * t239;
t182 = t314 * t215;
t214 = t238 * t239 - t242 * t243;
t181 = t314 * t214;
t77 = t181 * t248 - t201 * t215;
t78 = t181 * t213 + t200 * t215;
t313 = Ifges(6,5) * t77 + Ifges(6,6) * t78 + Ifges(6,3) * t182;
t178 = mrSges(6,1) * t248 + mrSges(6,2) * t213;
t312 = (-mrSges(4,2) * t242 + (t178 + t320) * t238) * t288;
t311 = 2 * m(5);
t310 = 2 * m(6);
t309 = 2 * m(7);
t308 = 0.2e1 * pkin(2);
t62 = -t118 * mrSges(7,1) + t117 * mrSges(7,2);
t307 = 0.2e1 * t62;
t233 = t235 ^ 2;
t187 = t222 * t238 - t223 * t242;
t256 = qJD(2) * t302;
t217 = t239 * t256;
t253 = t243 * t256;
t136 = qJD(3) * t187 + t217 * t238 - t242 * t253;
t306 = 0.2e1 * t136;
t195 = t200 * mrSges(6,2);
t285 = t201 * mrSges(6,1);
t162 = -t195 + t285;
t305 = 0.2e1 * t162;
t230 = -pkin(2) * t243 - pkin(1);
t304 = 0.2e1 * t230;
t303 = m(4) / 0.2e1;
t299 = pkin(2) * t242;
t298 = pkin(5) * t201;
t105 = pkin(2) * qJD(2) * t239 + pkin(3) * t182 + qJ(4) * t181 - qJD(4) * t215;
t135 = qJD(3) * t315 + t242 * t217 + t238 * t253;
t65 = t235 * t105 - t135 * t234;
t293 = t65 * mrSges(5,3);
t169 = t214 * pkin(3) - t215 * qJ(4) + t230;
t127 = t234 * t169 + t235 * t187;
t276 = t215 * t234;
t104 = -pkin(9) * t276 + t127;
t126 = t235 * t169 - t187 * t234;
t275 = t215 * t235;
t93 = pkin(4) * t214 - pkin(9) * t275 + t126;
t52 = t241 * t104 + t237 * t93;
t291 = Ifges(5,4) * t234;
t290 = Ifges(5,4) * t235;
t289 = Ifges(5,2) * t234;
t287 = pkin(5) * qJD(6);
t286 = t118 * mrSges(7,3);
t284 = t234 * t65;
t66 = t234 * t105 + t235 * t135;
t283 = t235 * t66;
t282 = t104 * t237;
t281 = t127 * t235;
t280 = t136 * t315;
t279 = t181 * t234;
t278 = t181 * t235;
t277 = t315 * t238;
t122 = -mrSges(5,1) * t279 - mrSges(5,2) * t278;
t268 = t234 ^ 2 + t233;
t267 = qJD(3) * t238;
t265 = qJD(6) * t236;
t264 = qJD(6) * t240;
t263 = 2 * mrSges(6,3);
t262 = 0.2e1 * mrSges(7,3);
t261 = 0.2e1 * t243;
t260 = mrSges(7,3) * t287;
t156 = t213 * t215;
t157 = t248 * t215;
t106 = -t156 * t240 + t157 * t236;
t25 = qJD(6) * t106 + t236 * t78 + t240 * t77;
t107 = -t156 * t236 - t157 * t240;
t26 = -qJD(6) * t107 - t236 * t77 + t240 * t78;
t259 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t182;
t258 = pkin(2) * t267;
t257 = t240 * t117 * mrSges(7,3);
t227 = -pkin(4) * t235 - pkin(3);
t41 = -t78 * mrSges(6,1) + t77 * mrSges(6,2);
t11 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t51 = t241 * t93 - t282;
t255 = t268 * t224;
t254 = t268 * qJD(4);
t251 = -t195 + t62;
t150 = pkin(4) * t276 - t315;
t250 = Ifges(5,5) * t235 - Ifges(5,6) * t234;
t37 = pkin(5) * t214 + pkin(10) * t157 + t51;
t42 = -pkin(10) * t156 + t52;
t16 = -t236 * t42 + t240 * t37;
t17 = t236 * t37 + t240 * t42;
t249 = 0.2e1 * t268 * mrSges(5,3);
t190 = pkin(5) * t248 + t227;
t48 = pkin(4) * t182 + pkin(9) * t278 + t65;
t58 = pkin(9) * t279 + t66;
t14 = -qJD(5) * t52 - t237 * t58 + t241 * t48;
t5 = pkin(5) * t182 - pkin(10) * t77 + t14;
t13 = -qJD(5) * t282 + t237 * t48 + t241 * t58 + t93 * t266;
t6 = pkin(10) * t78 + t13;
t3 = qJD(6) * t16 + t236 * t5 + t240 * t6;
t4 = -qJD(6) * t17 - t236 * t6 + t240 * t5;
t247 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t259;
t246 = t240 * t170 * t260 + (pkin(5) * t286 + t171 * t260) * t236 + t319;
t120 = Ifges(7,4) * t171 + Ifges(7,2) * t170;
t121 = Ifges(7,1) * t171 + Ifges(7,4) * t170;
t163 = -Ifges(6,4) * t200 - Ifges(6,2) * t201;
t164 = -Ifges(6,1) * t200 - Ifges(6,4) * t201;
t179 = Ifges(6,4) * t213 - Ifges(6,2) * t248;
t180 = Ifges(6,1) * t213 - Ifges(6,4) * t248;
t63 = Ifges(7,4) * t117 + Ifges(7,2) * t118;
t64 = Ifges(7,1) * t117 + Ifges(7,4) * t118;
t245 = t117 * t121 + t118 * t120 - t163 * t248 + t213 * t164 + t170 * t63 + t171 * t64 - t201 * t179 - t200 * t180;
t96 = -pkin(4) * t279 + t136;
t10 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t182;
t108 = pkin(5) * t156 + t150;
t35 = Ifges(6,4) * t77 + Ifges(6,2) * t78 + Ifges(6,6) * t182;
t36 = Ifges(6,1) * t77 + Ifges(6,4) * t78 + Ifges(6,5) * t182;
t50 = -pkin(5) * t78 + t96;
t56 = Ifges(7,4) * t107 + Ifges(7,2) * t106 + Ifges(7,6) * t214;
t57 = Ifges(7,1) * t107 + Ifges(7,4) * t106 + Ifges(7,5) * t214;
t89 = t182 * Ifges(5,6) - (-t289 + t290) * t181;
t9 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t182;
t90 = t182 * Ifges(5,5) - (Ifges(5,1) * t235 - t291) * t181;
t94 = -Ifges(6,4) * t157 - Ifges(6,2) * t156 + Ifges(6,6) * t214;
t95 = -Ifges(6,1) * t157 - Ifges(6,4) * t156 + Ifges(6,5) * t214;
t244 = (-t117 * t16 + t170 * t3 - t171 * t4) * mrSges(7,3) - (Ifges(5,1) * t234 + t290) * t278 / 0.2e1 + (Ifges(5,2) * t235 + t291) * t279 / 0.2e1 + (-t13 * t248 - t14 * t213 + t200 * t51 - t201 * t52) * mrSges(6,3) + (Ifges(5,5) * t234 + Ifges(6,5) * t213 + Ifges(7,5) * t171 + Ifges(5,6) * t235 - Ifges(6,6) * t248 + Ifges(7,6) * t170) * t182 / 0.2e1 - t248 * t35 / 0.2e1 + t320 * t136 + t319 * t214 / 0.2e1 + t17 * t286 + mrSges(5,3) * t283 + t235 * t89 / 0.2e1 + t234 * t90 / 0.2e1 + t213 * t36 / 0.2e1 - t200 * t95 / 0.2e1 - t201 * t94 / 0.2e1 - Ifges(4,6) * t182 + t96 * t178 + t78 * t179 / 0.2e1 + t77 * t180 / 0.2e1 - Ifges(4,5) * t181 + t150 * t162 - t156 * t163 / 0.2e1 - t157 * t164 / 0.2e1 + t170 * t9 / 0.2e1 + t171 * t10 / 0.2e1 - t135 * mrSges(4,2) + t26 * t120 / 0.2e1 + t25 * t121 / 0.2e1 + t117 * t57 / 0.2e1 + t118 * t56 / 0.2e1 + t50 * t119 + t106 * t63 / 0.2e1 + t107 * t64 / 0.2e1 + t108 * t62;
t229 = -pkin(3) - t299;
t218 = t227 - t299;
t207 = (-mrSges(7,1) * t236 - mrSges(7,2) * t240) * t287;
t189 = t258 + t298;
t188 = t190 - t299;
t173 = mrSges(5,1) * t214 - mrSges(5,3) * t275;
t172 = -mrSges(5,2) * t214 - mrSges(5,3) * t276;
t165 = (mrSges(5,1) * t234 + mrSges(5,2) * t235) * t215;
t138 = mrSges(6,1) * t214 + mrSges(6,3) * t157;
t137 = -mrSges(6,2) * t214 - mrSges(6,3) * t156;
t129 = mrSges(5,1) * t182 + mrSges(5,3) * t278;
t128 = -mrSges(5,2) * t182 + mrSges(5,3) * t279;
t109 = mrSges(6,1) * t156 - mrSges(6,2) * t157;
t92 = mrSges(7,1) * t214 - mrSges(7,3) * t107;
t91 = -mrSges(7,2) * t214 + mrSges(7,3) * t106;
t68 = -mrSges(6,2) * t182 + mrSges(6,3) * t78;
t67 = mrSges(6,1) * t182 - mrSges(6,3) * t77;
t59 = -mrSges(7,1) * t106 + mrSges(7,2) * t107;
t21 = -mrSges(7,2) * t182 + mrSges(7,3) * t26;
t20 = mrSges(7,1) * t182 - mrSges(7,3) * t25;
t1 = [(-0.2e1 * mrSges(4,3) * t135 + t259 + t313) * t214 + 0.2e1 * (t181 * t315 - t182 * t187) * mrSges(4,3) - 0.2e1 * t315 * t122 + (mrSges(4,1) * t304 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(6,3) + Ifges(7,3)) * t214 - Ifges(6,5) * t157 - Ifges(6,6) * t156 + Ifges(7,5) * t107 + Ifges(7,6) * t106 + (-(2 * Ifges(4,4)) + t250) * t215) * t182 + 0.2e1 * m(4) * (t135 * t187 - t280) + (t126 * t65 + t127 * t66 - t280) * t311 + (mrSges(4,3) * t306 - t234 * t89 + t235 * t90 - (Ifges(5,1) * t233 + (2 * Ifges(4,1)) + (t289 - 0.2e1 * t290) * t234) * t181) * t215 + t165 * t306 + (t108 * t50 + t16 * t4 + t17 * t3) * t309 + (t13 * t52 + t14 * t51 + t150 * t96) * t310 - 0.2e1 * (-Ifges(4,4) + t250) * t181 * t214 + 0.2e1 * t66 * t172 + 0.2e1 * t65 * t173 - t156 * t35 - t157 * t36 + 0.2e1 * t150 * t41 + 0.2e1 * t13 * t137 + 0.2e1 * t14 * t138 + 0.2e1 * t127 * t128 + 0.2e1 * t126 * t129 + t106 * t9 + t107 * t10 + 0.2e1 * t108 * t11 + 0.2e1 * t96 * t109 + 0.2e1 * t3 * t91 + 0.2e1 * t4 * t92 + t78 * t94 + t77 * t95 - mrSges(4,2) * t181 * t304 + 0.2e1 * t51 * t67 + 0.2e1 * t52 * t68 + 0.2e1 * t50 * t59 + t26 * t56 + t25 * t57 + 0.2e1 * t16 * t20 + 0.2e1 * t17 * t21 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t243) * t261 + (-0.2e1 * pkin(1) * mrSges(3,1) + (mrSges(4,1) * t214 + mrSges(4,2) * t215) * t308 + m(4) * pkin(2) * t304 - 0.2e1 * Ifges(3,4) * t239 + (Ifges(3,1) - Ifges(3,2)) * t261) * t239) * qJD(2); t244 + ((t135 * t238 - t136 * t242) * t303 + (m(6) * t150 * t238 / 0.2e1 - m(5) * t277 / 0.2e1 + (t187 * t242 - t277) * t303) * qJD(3)) * t308 + m(6) * (t124 * t52 + t125 * t51 + t13 * t161 + t14 * t160 + t218 * t96) + m(5) * (-t126 * t273 + t136 * t229 + t224 * t281 - t226 * t284 + t66 * t272) + ((t109 + t165) * t267 + (t181 * t242 - t182 * t238 + (-t214 * t242 + t215 * t238) * qJD(3)) * mrSges(4,3)) * pkin(2) + (Ifges(3,5) * t243 - Ifges(3,6) * t239 + (-mrSges(3,1) * t243 + mrSges(3,2) * t239) * pkin(7)) * qJD(2) + (t226 * t128 + t224 * t172) * t235 + (-t226 * t129 - t224 * t173 - t293) * t234 + m(7) * (t108 * t189 + t16 * t29 + t17 * t28 + t188 * t50 + t3 * t82 + t4 * t81) + t229 * t122 + t218 * t41 + t188 * t11 + t189 * t59 + t160 * t67 + t161 * t68 + t124 * t137 + t125 * t138 + t28 * t91 + t29 * t92 + t81 * t20 + t82 * t21; t245 + (t226 * t255 + t229 * t258) * t311 + (t188 * t189 + t28 * t82 + t29 * t81) * t309 + (t124 * t161 + t125 * t160 + t218 * t258) * t310 + t224 * t249 + 0.2e1 * t312 + (-t117 * t81 + t118 * t82 + t170 * t28 - t171 * t29) * t262 + (-t124 * t248 - t125 * t213 + t160 * t200 - t161 * t201) * t263 + t218 * t305 + t188 * t307 + t189 * t318; t244 + m(5) * (-pkin(3) * t136 + (-t126 * t234 + t281) * qJD(4) + (t283 - t284) * qJ(4)) + t227 * t41 + t184 * t67 + t185 * t68 + t190 * t11 + t146 * t137 + t147 * t138 - pkin(3) * t122 + t97 * t20 + t98 * t21 + t44 * t91 + t45 * t92 + (qJ(4) * t128 + qJD(4) * t172) * t235 + (-qJ(4) * t129 - qJD(4) * t173 - t293) * t234 + t59 * t298 + m(7) * (t108 * t298 + t16 * t45 + t17 * t44 + t190 * t50 + t3 * t98 + t4 * t97) + m(6) * (t13 * t185 + t14 * t184 + t146 * t52 + t147 * t51 + t227 * t96); t245 + m(5) * (-pkin(3) * t258 + qJ(4) * t255 + t226 * t254) + m(7) * (t188 * t298 + t189 * t190 + t28 * t98 + t29 * t97 + t44 * t82 + t45 * t81) + m(6) * (t124 * t185 + t125 * t184 + t146 * t161 + t147 * t160 + t227 * t258) + (t188 + t190) * t62 + (t227 + t218) * t162 + (t189 + t298) * t119 + t312 + (t255 + t254) * mrSges(5,3) + ((-t29 - t45) * t171 + (t28 + t44) * t170 + (t82 + t98) * t118 + (-t81 - t97) * t117) * mrSges(7,3) + ((-t125 - t147) * t213 - (t124 + t146) * t248 - (t161 + t185) * t201 - (-t160 - t184) * t200) * mrSges(6,3); t245 + (qJ(4) * t268 * t311 + t249) * qJD(4) + (t190 * t298 + t44 * t98 + t45 * t97) * t309 + (t146 * t185 + t147 * t184) * t310 + (-t117 * t97 + t118 * t98 + t170 * t44 - t171 * t45) * t262 + (-t146 * t248 - t147 * t213 + t184 * t200 - t185 * t201) * t263 + t227 * t305 + t190 * t307 + t298 * t318; m(5) * t136 + m(6) * t96 + m(7) * t50 + t11 + t122 + t41; m(7) * t189 + t285 + (m(5) + m(6)) * t258 + t251; -(-m(7) * pkin(5) - mrSges(6,1)) * t201 + t251; 0; t14 * mrSges(6,1) - t13 * mrSges(6,2) + (t91 * t264 + t236 * t21 - t92 * t265 + t240 * t20 + m(7) * (-t16 * t265 + t17 * t264 + t236 * t3 + t240 * t4)) * pkin(5) + t247 + t313; t125 * mrSges(6,1) - t124 * mrSges(6,2) + (m(7) * (t236 * t28 + t240 * t29 + t264 * t82 - t265 * t81) - t257) * pkin(5) + t246 + t317; t147 * mrSges(6,1) - t146 * mrSges(6,2) + (m(7) * (t236 * t44 + t240 * t45 + t264 * t98 - t265 * t97) - t257) * pkin(5) + t246 + t316; 0; 0.2e1 * t207; t247; t270 + t317; t270 + t316; 0; t207; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
