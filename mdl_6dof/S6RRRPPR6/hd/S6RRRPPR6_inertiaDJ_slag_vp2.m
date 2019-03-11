% Calculate time derivative of joint inertia matrix for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:49
% EndTime: 2019-03-09 15:48:01
% DurationCPUTime: 6.27s
% Computational Cost: add. (7634->587), mult. (20053->841), div. (0->0), fcn. (19286->10), ass. (0->247)
t312 = Ifges(5,5) - Ifges(6,4);
t311 = -Ifges(5,6) + Ifges(6,5);
t214 = cos(qJ(3));
t270 = cos(pkin(11));
t236 = t270 * t214;
t207 = sin(pkin(11));
t211 = sin(qJ(3));
t265 = t207 * t211;
t172 = -t236 + t265;
t213 = cos(qJ(6));
t253 = qJD(6) * t213;
t237 = t270 * t211;
t173 = t207 * t214 + t237;
t166 = t173 * qJD(3);
t210 = sin(qJ(6));
t262 = t210 * t166;
t225 = t172 * t253 + t262;
t254 = qJD(6) * t210;
t261 = t213 * t166;
t224 = t172 * t254 - t261;
t310 = Ifges(6,1) + Ifges(4,3) + Ifges(5,3);
t309 = 2 * mrSges(6,1) + 2 * mrSges(5,3);
t290 = -t210 / 0.2e1;
t288 = t213 / 0.2e1;
t209 = cos(pkin(6));
t208 = sin(pkin(6));
t212 = sin(qJ(2));
t264 = t208 * t212;
t169 = t209 * t211 + t214 * t264;
t215 = cos(qJ(2));
t258 = qJD(2) * t208;
t246 = t215 * t258;
t134 = -qJD(3) * t169 - t211 * t246;
t263 = t208 * t215;
t171 = t209 * t212 * pkin(1) + pkin(8) * t263;
t157 = t171 * qJD(2);
t101 = -t134 * pkin(3) + t157;
t168 = t209 * t214 - t211 * t264;
t113 = t207 * t168 + t169 * t270;
t135 = qJD(3) * t168 + t214 * t246;
t83 = t207 * t134 + t135 * t270;
t216 = -t83 * qJ(5) - t113 * qJD(5) + t101;
t296 = pkin(4) + pkin(10);
t82 = -t134 * t270 + t135 * t207;
t14 = t296 * t82 + t216;
t257 = qJD(2) * t212;
t247 = t208 * t257;
t151 = pkin(9) * t209 + t171;
t152 = (-pkin(2) * t215 - pkin(9) * t212 - pkin(1)) * t208;
t103 = t214 * t151 + t211 * t152;
t155 = (pkin(2) * t212 - pkin(9) * t215) * t258;
t195 = pkin(8) * t264;
t286 = pkin(1) * t215;
t170 = t209 * t286 - t195;
t156 = t170 * qJD(2);
t56 = -t103 * qJD(3) + t214 * t155 - t156 * t211;
t30 = pkin(3) * t247 - qJ(4) * t135 - qJD(4) * t169 + t56;
t255 = qJD(3) * t214;
t256 = qJD(3) * t211;
t55 = -t151 * t256 + t152 * t255 + t211 * t155 + t214 * t156;
t39 = qJ(4) * t134 + qJD(4) * t168 + t55;
t232 = t207 * t39 - t270 * t30;
t3 = t83 * pkin(5) - t247 * t296 + t232;
t102 = -t211 * t151 + t214 * t152;
t70 = -pkin(3) * t263 - t169 * qJ(4) + t102;
t84 = qJ(4) * t168 + t103;
t37 = -t207 * t84 + t270 * t70;
t29 = pkin(4) * t263 - t37;
t21 = t113 * pkin(5) + pkin(10) * t263 + t29;
t112 = -t168 * t270 + t169 * t207;
t150 = t195 + (-pkin(2) - t286) * t209;
t114 = -t168 * pkin(3) + t150;
t219 = -t113 * qJ(5) + t114;
t23 = t112 * t296 + t219;
t5 = t21 * t213 - t210 * t23;
t1 = qJD(6) * t5 + t14 * t213 + t210 * t3;
t6 = t21 * t210 + t213 * t23;
t2 = -qJD(6) * t6 - t14 * t210 + t213 * t3;
t306 = -t1 * t210 - t2 * t213;
t308 = m(7) * ((-t210 * t5 + t213 * t6) * qJD(6) - t306);
t280 = -qJ(4) - pkin(9);
t183 = t280 * t214;
t136 = -t183 * t207 - t280 * t237;
t108 = pkin(5) * t173 + t136;
t203 = -pkin(3) * t214 - pkin(2);
t227 = -qJ(5) * t173 + t203;
t96 = t172 * t296 + t227;
t48 = t108 * t213 - t210 * t96;
t167 = qJD(3) * t236 - t207 * t256;
t206 = pkin(3) * t256;
t223 = -qJ(5) * t167 - qJD(5) * t173 + t206;
t69 = t166 * t296 + t223;
t238 = qJD(3) * t280;
t165 = qJD(4) * t214 + t211 * t238;
t221 = -qJD(4) * t211 + t214 * t238;
t106 = t165 * t207 - t270 * t221;
t86 = pkin(5) * t167 + t106;
t16 = qJD(6) * t48 + t210 * t86 + t213 * t69;
t49 = t108 * t210 + t213 * t96;
t17 = -qJD(6) * t49 - t210 * t69 + t213 * t86;
t307 = -t16 * t210 - t17 * t213;
t305 = Ifges(4,5) * t255 - Ifges(4,6) * t256 + t311 * t166 + t312 * t167;
t248 = t270 * pkin(3);
t202 = -t248 - pkin(4);
t304 = m(6) * t202 + mrSges(6,2);
t13 = t207 * t30 + t270 * t39;
t10 = -t208 * (qJ(5) * t257 - qJD(5) * t215) - t13;
t303 = 2 * m(5);
t302 = 0.2e1 * m(6);
t301 = 0.2e1 * m(7);
t300 = -2 * mrSges(3,3);
t299 = 0.2e1 * t157;
t298 = m(5) * pkin(3);
t91 = t213 * t112 + t210 * t263;
t297 = t91 / 0.2e1;
t229 = Ifges(7,5) * t210 + Ifges(7,6) * t213;
t295 = -t229 * qJD(6) / 0.2e1;
t276 = Ifges(7,4) * t213;
t231 = Ifges(7,1) * t210 + t276;
t180 = t231 * qJD(6);
t294 = -t180 / 0.2e1;
t293 = Ifges(7,5) * t288 + Ifges(7,6) * t290;
t185 = -Ifges(7,2) * t210 + t276;
t292 = t185 / 0.2e1;
t277 = Ifges(7,4) * t210;
t187 = Ifges(7,1) * t213 - t277;
t291 = t187 / 0.2e1;
t289 = t210 / 0.2e1;
t285 = pkin(3) * t207;
t281 = -mrSges(6,2) + mrSges(5,1);
t38 = t207 * t70 + t270 * t84;
t279 = Ifges(4,4) * t211;
t278 = Ifges(4,4) * t214;
t275 = t156 * mrSges(3,2);
t274 = t157 * mrSges(3,1);
t272 = t167 * mrSges(6,1);
t269 = t172 * t210;
t268 = t172 * t213;
t199 = -pkin(10) + t202;
t267 = t199 * t210;
t266 = t199 * t213;
t252 = 0.2e1 * t208;
t226 = -t210 * t112 + t213 * t263;
t44 = qJD(6) * t226 - t210 * t247 + t213 * t82;
t45 = qJD(6) * t91 + t210 * t82 + t213 * t247;
t7 = Ifges(7,5) * t45 + Ifges(7,6) * t44 + Ifges(7,3) * t83;
t250 = mrSges(7,3) * t254;
t249 = mrSges(7,3) * t253;
t245 = t199 * t254;
t244 = t199 * t253;
t240 = -t254 / 0.2e1;
t239 = -t253 / 0.2e1;
t65 = t83 * mrSges(6,1) + mrSges(6,2) * t247;
t182 = mrSges(7,1) * t210 + mrSges(7,2) * t213;
t230 = Ifges(7,2) * t213 + t277;
t107 = t165 * t270 + t207 * t221;
t137 = -t183 * t270 + t265 * t280;
t228 = t106 * t136 + t107 * t137;
t25 = qJ(5) * t263 - t38;
t220 = Ifges(4,5) * t135 + Ifges(4,6) * t134 + t310 * t247 + t311 * t82 + t312 * t83;
t50 = Ifges(7,5) * t225 - Ifges(7,6) * t224 + Ifges(7,3) * t167;
t217 = (-t210 * t48 + t213 * t49) * qJD(6) - t307;
t200 = qJ(5) + t285;
t192 = Ifges(3,5) * t246;
t188 = Ifges(4,1) * t211 + t278;
t186 = Ifges(4,2) * t214 + t279;
t181 = (Ifges(4,1) * t214 - t279) * qJD(3);
t179 = (-Ifges(4,2) * t211 + t278) * qJD(3);
t178 = t230 * qJD(6);
t176 = (mrSges(4,1) * t211 + mrSges(4,2) * t214) * qJD(3);
t175 = -mrSges(7,1) * t253 + mrSges(7,2) * t254;
t159 = t167 * mrSges(6,3);
t158 = t167 * mrSges(5,2);
t139 = -mrSges(4,1) * t263 - t169 * mrSges(4,3);
t138 = mrSges(4,2) * t263 + t168 * mrSges(4,3);
t131 = Ifges(5,1) * t173 - Ifges(5,4) * t172;
t130 = Ifges(5,4) * t173 - Ifges(5,2) * t172;
t129 = -Ifges(6,2) * t173 + Ifges(6,6) * t172;
t128 = -Ifges(6,6) * t173 + Ifges(6,3) * t172;
t127 = -mrSges(6,2) * t172 - mrSges(6,3) * t173;
t126 = mrSges(5,1) * t172 + mrSges(5,2) * t173;
t124 = -mrSges(7,2) * t173 + mrSges(7,3) * t268;
t123 = mrSges(7,1) * t173 - mrSges(7,3) * t269;
t122 = pkin(4) * t172 + t227;
t121 = (-mrSges(7,1) * t213 + mrSges(7,2) * t210) * t172;
t120 = Ifges(5,1) * t167 - Ifges(5,4) * t166;
t119 = Ifges(5,4) * t167 - Ifges(5,2) * t166;
t118 = -Ifges(6,2) * t167 + Ifges(6,6) * t166;
t117 = -Ifges(6,6) * t167 + Ifges(6,3) * t166;
t116 = -t166 * mrSges(6,2) - t159;
t115 = t166 * mrSges(5,1) + t158;
t111 = mrSges(4,1) * t247 - mrSges(4,3) * t135;
t110 = -mrSges(4,2) * t247 + mrSges(4,3) * t134;
t109 = -t172 * pkin(5) + t137;
t105 = Ifges(4,1) * t169 + Ifges(4,4) * t168 - Ifges(4,5) * t263;
t104 = Ifges(4,4) * t169 + Ifges(4,2) * t168 - Ifges(4,6) * t263;
t100 = -mrSges(5,1) * t263 - t113 * mrSges(5,3);
t99 = mrSges(5,2) * t263 - t112 * mrSges(5,3);
t98 = t113 * mrSges(6,1) - mrSges(6,2) * t263;
t97 = t112 * mrSges(6,1) + mrSges(6,3) * t263;
t95 = Ifges(7,5) * t173 + t172 * t231;
t94 = Ifges(7,6) * t173 + t172 * t230;
t93 = Ifges(7,3) * t173 + t172 * t229;
t90 = pkin(4) * t166 + t223;
t89 = mrSges(7,1) * t167 - mrSges(7,3) * t225;
t88 = -mrSges(7,2) * t167 - mrSges(7,3) * t224;
t87 = -t166 * pkin(5) + t107;
t85 = -mrSges(4,1) * t134 + mrSges(4,2) * t135;
t76 = t83 * mrSges(6,3);
t75 = t83 * mrSges(5,2);
t72 = Ifges(4,1) * t135 + Ifges(4,4) * t134 + Ifges(4,5) * t247;
t71 = Ifges(4,4) * t135 + Ifges(4,2) * t134 + Ifges(4,6) * t247;
t67 = mrSges(5,1) * t247 - mrSges(5,3) * t83;
t66 = -mrSges(5,2) * t247 - mrSges(5,3) * t82;
t64 = mrSges(6,1) * t82 - mrSges(6,3) * t247;
t63 = mrSges(7,1) * t224 + mrSges(7,2) * t225;
t62 = -mrSges(6,2) * t112 - mrSges(6,3) * t113;
t61 = mrSges(5,1) * t112 + mrSges(5,2) * t113;
t60 = Ifges(5,1) * t113 - Ifges(5,4) * t112 - Ifges(5,5) * t263;
t59 = Ifges(5,4) * t113 - Ifges(5,2) * t112 - Ifges(5,6) * t263;
t58 = -Ifges(6,4) * t263 - Ifges(6,2) * t113 + Ifges(6,6) * t112;
t57 = -Ifges(6,5) * t263 - Ifges(6,6) * t113 + Ifges(6,3) * t112;
t54 = mrSges(7,1) * t113 + mrSges(7,3) * t226;
t53 = -mrSges(7,2) * t113 + mrSges(7,3) * t91;
t52 = Ifges(7,1) * t225 - Ifges(7,4) * t224 + Ifges(7,5) * t167;
t51 = Ifges(7,4) * t225 - Ifges(7,2) * t224 + Ifges(7,6) * t167;
t47 = -mrSges(7,1) * t91 - mrSges(7,2) * t226;
t46 = t112 * pkin(4) + t219;
t43 = -t82 * mrSges(6,2) - t76;
t42 = t82 * mrSges(5,1) + t75;
t34 = Ifges(5,1) * t83 - Ifges(5,4) * t82 + Ifges(5,5) * t247;
t33 = Ifges(5,4) * t83 - Ifges(5,2) * t82 + Ifges(5,6) * t247;
t32 = Ifges(6,4) * t247 - Ifges(6,2) * t83 + Ifges(6,6) * t82;
t31 = Ifges(6,5) * t247 - Ifges(6,6) * t83 + Ifges(6,3) * t82;
t28 = -Ifges(7,1) * t226 + Ifges(7,4) * t91 + Ifges(7,5) * t113;
t27 = -Ifges(7,4) * t226 + Ifges(7,2) * t91 + Ifges(7,6) * t113;
t26 = -Ifges(7,5) * t226 + Ifges(7,6) * t91 + Ifges(7,3) * t113;
t22 = -t112 * pkin(5) - t25;
t20 = mrSges(7,1) * t83 - mrSges(7,3) * t45;
t19 = -mrSges(7,2) * t83 + mrSges(7,3) * t44;
t18 = t82 * pkin(4) + t216;
t15 = -mrSges(7,1) * t44 + mrSges(7,2) * t45;
t11 = -pkin(4) * t247 + t232;
t9 = Ifges(7,1) * t45 + Ifges(7,4) * t44 + Ifges(7,5) * t83;
t8 = Ifges(7,4) * t45 + Ifges(7,2) * t44 + Ifges(7,6) * t83;
t4 = -t82 * pkin(5) - t10;
t12 = [(t192 - 0.2e1 * t274 - 0.2e1 * t275) * t209 + (t7 + t34 - t32) * t113 + (t31 - t33) * t112 + 0.2e1 * m(3) * (t156 * t171 - t157 * t170) + 0.2e1 * m(4) * (t102 * t56 + t103 * t55 + t150 * t157) + (t26 + t60 - t58) * t83 + (t57 - t59) * t82 + (mrSges(3,3) * t212 * t299 + (0.2e1 * t156 * mrSges(3,3) - t220) * t215 + ((t170 * t300 + Ifges(3,5) * t209 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t215) * t252) * t215 + (t171 * t300 + Ifges(4,5) * t169 - 0.2e1 * Ifges(3,6) * t209 + Ifges(4,6) * t168 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t212) * t252 + t312 * t113 + t311 * t112 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t310) * t263) * t212) * qJD(2)) * t208 + (t101 * t114 + t13 * t38 - t232 * t37) * t303 - 0.2e1 * t232 * t100 - t226 * t9 + (-mrSges(4,1) * t168 + mrSges(4,2) * t169) * t299 + (t1 * t6 + t2 * t5 + t22 * t4) * t301 + (t10 * t25 + t11 * t29 + t18 * t46) * t302 + 0.2e1 * t6 * t19 + 0.2e1 * t5 * t20 + 0.2e1 * t22 * t15 + t44 * t27 + t45 * t28 + 0.2e1 * t46 * t43 + 0.2e1 * t4 * t47 + 0.2e1 * t1 * t53 + 0.2e1 * t2 * t54 + 0.2e1 * t18 * t62 + 0.2e1 * t25 * t64 + 0.2e1 * t29 * t65 + 0.2e1 * t38 * t66 + 0.2e1 * t37 * t67 + t91 * t8 + 0.2e1 * t10 * t97 + 0.2e1 * t11 * t98 + 0.2e1 * t13 * t99 + 0.2e1 * t101 * t61 + 0.2e1 * t103 * t110 + 0.2e1 * t102 * t111 + 0.2e1 * t114 * t42 + t134 * t104 + t135 * t105 + 0.2e1 * t55 * t138 + 0.2e1 * t56 * t139 + 0.2e1 * t150 * t85 + t168 * t71 + t169 * t72; ((-t211 * t56 + t214 * t55 + (-t102 * t214 - t103 * t211) * qJD(3)) * pkin(9) - pkin(2) * t157) * m(4) + t192 + (-t37 * mrSges(5,3) + t29 * mrSges(6,1) + t26 / 0.2e1 + t60 / 0.2e1 - t58 / 0.2e1) * t167 - t274 - t275 + (t50 / 0.2e1 - t118 / 0.2e1 + t120 / 0.2e1) * t113 + (t117 / 0.2e1 - t119 / 0.2e1) * t112 + (pkin(9) * t110 + t55 * mrSges(4,3) + t71 / 0.2e1 - t157 * mrSges(4,1)) * t214 + (-pkin(9) * t111 - t56 * mrSges(4,3) + t72 / 0.2e1 + t157 * mrSges(4,2)) * t211 + (t66 - t64) * t137 + (t65 - t67) * t136 + (t99 - t97) * t107 + (t98 - t100) * t106 + m(5) * (t101 * t203 - t106 * t37 + t107 * t38 + t114 * t206 + t13 * t137 + t136 * t232) + (t232 * mrSges(5,3) + t11 * mrSges(6,1) + t7 / 0.2e1 + t34 / 0.2e1 - t32 / 0.2e1) * t173 - t226 * t52 / 0.2e1 + m(6) * (-t10 * t137 + t106 * t29 - t107 * t25 + t11 * t136 + t122 * t18 + t46 * t90) + m(7) * (t1 * t49 + t109 * t4 + t16 * t6 + t17 * t5 + t2 * t48 + t22 * t87) + (t93 / 0.2e1 - t129 / 0.2e1 + t131 / 0.2e1) * t83 + (t128 / 0.2e1 - t130 / 0.2e1) * t82 + ((-Ifges(3,6) + Ifges(4,5) * t211 / 0.2e1 + Ifges(4,6) * t214 / 0.2e1 + (Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t173 + (-Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1) * t172) * t257 - t305 * t215 / 0.2e1) * t208 + (t27 * t288 + t28 * t289 + t25 * mrSges(6,1) - t38 * mrSges(5,3) + t57 / 0.2e1 - t59 / 0.2e1) * t166 + (t8 * t288 + t9 * t289 + t10 * mrSges(6,1) - t13 * mrSges(5,3) + t31 / 0.2e1 - t33 / 0.2e1 + (t27 * t290 + t28 * t288) * qJD(6)) * t172 + t51 * t297 + ((-pkin(9) * t139 - t102 * mrSges(4,3) + t105 / 0.2e1) * t214 + (-pkin(9) * t138 - t103 * mrSges(4,3) + pkin(3) * t61 - t104 / 0.2e1) * t211) * qJD(3) + t48 * t20 + t49 * t19 + t16 * t53 + t17 * t54 + t22 * t63 - pkin(2) * t85 + t87 * t47 + t6 * t88 + t5 * t89 + t90 * t62 + t44 * t94 / 0.2e1 + t45 * t95 / 0.2e1 + t109 * t15 + t114 * t115 + t46 * t116 + t4 * t121 + t122 * t43 + t2 * t123 + t1 * t124 + t101 * t126 + t18 * t127 + t150 * t176 + t168 * t179 / 0.2e1 + t169 * t181 / 0.2e1 + t134 * t186 / 0.2e1 + t135 * t188 / 0.2e1 + t203 * t42; -0.2e1 * pkin(2) * t176 + 0.2e1 * t109 * t63 + 0.2e1 * t203 * t115 + 0.2e1 * t122 * t116 + 0.2e1 * t87 * t121 + 0.2e1 * t17 * t123 + 0.2e1 * t16 * t124 + 0.2e1 * t90 * t127 + t214 * t179 + t211 * t181 + 0.2e1 * t48 * t89 + 0.2e1 * t49 * t88 + (t214 * t188 + (0.2e1 * pkin(3) * t126 - t186) * t211) * qJD(3) + (t203 * t206 + t228) * t303 + (t122 * t90 + t228) * t302 + (t109 * t87 + t16 * t49 + t17 * t48) * t301 + (t106 * t309 - t118 + t120 + t50) * t173 + (t136 * t309 - t129 + t131 + t93) * t167 + (-t137 * t309 + t210 * t95 + t213 * t94 + t128 - t130) * t166 + (t210 * t52 + t213 * t51 + t117 - t119 - t107 * t309 + (-t210 * t94 + t213 * t95) * qJD(6)) * t172; -t6 * t249 + t28 * t240 + t27 * t239 + t5 * t250 + t53 * t244 + t199 * t308 + (t13 * t207 - t232 * t270) * t298 - t232 * mrSges(5,1) - t226 * t294 + t306 * mrSges(7,3) + t304 * t11 + (-m(6) * t10 + m(7) * t4 + t15 - t64) * t200 + (-m(6) * t25 + m(7) * t22 + t47 - t97) * qJD(5) - t54 * t245 + t67 * t248 + t20 * t266 + t19 * t267 + t66 * t285 + t9 * t288 + t8 * t290 + t45 * t291 + t44 * t292 + t83 * t293 + t113 * t295 - t178 * t297 + t220 - t10 * mrSges(6,3) - t13 * mrSges(5,2) - t55 * mrSges(4,2) + t56 * mrSges(4,1) - t22 * t175 + t4 * t182 + t202 * t65; t305 - t49 * t249 + t94 * t239 + (qJD(6) * t187 - t178) * t268 / 0.2e1 + t48 * t250 + t124 * t244 + (m(6) * t137 + m(7) * t109 + t121) * qJD(5) + (t172 * t185 + t95) * t240 + (m(6) * t200 + t207 * t298 - mrSges(5,2) + mrSges(6,3)) * t107 + m(7) * (t199 * t217 + t200 * t87) + (-mrSges(4,1) * t255 + mrSges(4,2) * t256) * pkin(9) + (-t166 * t285 - t167 * t248) * mrSges(5,3) + (-qJD(5) * t172 - t166 * t200) * mrSges(6,1) + t307 * mrSges(7,3) + (-t270 * t298 - mrSges(5,1) + t304) * t106 - t123 * t245 + t89 * t266 + t88 * t267 + t202 * t272 + t52 * t288 + t51 * t290 + t262 * t291 + t261 * t292 + t167 * t293 + t269 * t294 + t173 * t295 - t109 * t175 + t87 * t182 + t200 * t63; -0.2e1 * t175 * t200 + t178 * t210 - t180 * t213 + (-t185 * t213 - t187 * t210) * qJD(6) + 0.2e1 * (mrSges(6,3) + t182 + (m(6) + m(7)) * t200) * qJD(5); t213 * t19 - t210 * t20 + t75 - t76 + t281 * t82 + (-t210 * t53 - t213 * t54) * qJD(6) + m(7) * (t1 * t213 - t2 * t210 + (-t210 * t6 - t213 * t5) * qJD(6)) + m(6) * t18 + m(5) * t101; m(5) * t206 - t210 * t89 + t213 * t88 + t158 - t159 + t281 * t166 + (-t213 * t123 - t210 * t124) * qJD(6) + m(7) * (t16 * t213 - t17 * t210 + (-t210 * t49 - t213 * t48) * qJD(6)) + m(6) * t90; 0; 0; t210 * t19 + t213 * t20 + (-t210 * t54 + t213 * t53) * qJD(6) + t308 + m(6) * t11 + t65; t272 + t210 * t88 + t213 * t89 + (-t123 * t210 + t124 * t213) * qJD(6) + m(7) * t217 + m(6) * t106; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t17 - mrSges(7,2) * t16 + t50; (-t182 * t199 - t229) * qJD(6); t175; -t182 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
