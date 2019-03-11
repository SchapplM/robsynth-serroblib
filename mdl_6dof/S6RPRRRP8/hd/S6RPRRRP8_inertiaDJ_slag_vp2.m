% Calculate time derivative of joint inertia matrix for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:58
% EndTime: 2019-03-09 06:23:07
% DurationCPUTime: 4.34s
% Computational Cost: add. (3901->349), mult. (8014->483), div. (0->0), fcn. (6857->6), ass. (0->168)
t318 = Ifges(6,1) + Ifges(7,1);
t298 = Ifges(7,4) + Ifges(6,5);
t156 = sin(qJ(5));
t157 = cos(qJ(5));
t293 = t156 ^ 2 + t157 ^ 2;
t255 = Ifges(7,5) * t156;
t257 = Ifges(6,4) * t156;
t317 = t318 * t157 + t255 - t257;
t316 = mrSges(7,2) + mrSges(6,3);
t315 = Ifges(7,2) + Ifges(6,3);
t254 = Ifges(7,5) * t157;
t256 = Ifges(6,4) * t157;
t314 = t318 * t156 - t254 + t256;
t268 = sin(qJ(4));
t269 = sin(qJ(3));
t270 = cos(qJ(4));
t271 = cos(qJ(3));
t117 = t268 * t271 + t270 * t269;
t233 = qJD(5) * t157;
t234 = qJD(5) * t156;
t206 = qJD(3) * t271;
t135 = pkin(3) * t206 + qJD(2);
t300 = -qJD(3) - qJD(4);
t83 = t300 * t117;
t116 = t268 * t269 - t270 * t271;
t84 = t300 * t116;
t31 = pkin(4) * t84 - pkin(9) * t83 + t135;
t158 = -pkin(1) - pkin(7);
t297 = -pkin(8) + t158;
t125 = t297 * t271;
t115 = t125 * qJD(3);
t124 = t297 * t269;
t167 = qJD(3) * t124;
t203 = qJD(4) * t268;
t204 = qJD(4) * t270;
t36 = t115 * t270 - t124 * t203 + t125 * t204 - t167 * t268;
t142 = t269 * pkin(3) + qJ(2);
t62 = pkin(4) * t117 + pkin(9) * t116 + t142;
t90 = t124 * t270 + t125 * t268;
t6 = t156 * t31 + t157 * t36 + t62 * t233 - t234 * t90;
t2 = qJ(6) * t84 + qJD(6) * t117 + t6;
t258 = t156 * t62 + t157 * t90;
t7 = -qJD(5) * t258 - t156 * t36 + t157 * t31;
t4 = -pkin(5) * t84 - t7;
t188 = t156 * t4 + t157 * t2;
t26 = qJ(6) * t117 + t258;
t32 = -t156 * t90 + t157 * t62;
t27 = -pkin(5) * t117 - t32;
t313 = t27 * t233 - t26 * t234 + t188;
t312 = Ifges(6,6) * t157 + t156 * t298;
t211 = t116 * t234;
t244 = t157 * t83;
t169 = t211 + t244;
t245 = t156 * t83;
t303 = -t116 * t233 + t245;
t311 = t298 * t84 + (-Ifges(6,4) + Ifges(7,5)) * t303 + t318 * t169;
t310 = -t317 * t116 + t298 * t117;
t309 = t317 * qJD(5);
t127 = -t157 * mrSges(7,1) - t156 * mrSges(7,3);
t128 = -t157 * mrSges(6,1) + t156 * mrSges(6,2);
t308 = -t127 - t128;
t198 = pkin(3) * t204;
t307 = t293 * t198;
t144 = pkin(3) * t268 + pkin(9);
t290 = -t144 * t233 - t156 * t198;
t305 = t144 * t234 - t157 * t198;
t302 = t293 * t84;
t22 = t84 * mrSges(6,1) - mrSges(6,3) * t169;
t23 = -t84 * mrSges(7,1) + t169 * mrSges(7,2);
t24 = -t84 * mrSges(6,2) - mrSges(6,3) * t303;
t25 = -mrSges(7,2) * t303 + t84 * mrSges(7,3);
t301 = (t24 + t25) * t157 + (-t22 + t23) * t156;
t299 = -Ifges(4,1) + Ifges(4,2);
t241 = t116 * t156;
t73 = -mrSges(6,2) * t117 + mrSges(6,3) * t241;
t76 = mrSges(7,2) * t241 + mrSges(7,3) * t117;
t260 = t73 + t76;
t240 = t116 * t157;
t74 = mrSges(6,1) * t117 + mrSges(6,3) * t240;
t75 = -mrSges(7,1) * t117 - mrSges(7,2) * t240;
t259 = -t74 + t75;
t232 = qJD(6) * t157;
t235 = Ifges(7,4) * t233 + Ifges(7,6) * t234;
t294 = mrSges(7,2) * t232 + t235;
t205 = qJD(3) * t269;
t292 = -mrSges(4,1) * t205 - mrSges(4,2) * t206;
t291 = t303 * Ifges(7,6) + t298 * t244 + t315 * t84;
t288 = -t268 * t84 - t270 * t83;
t266 = t156 * t7;
t287 = -t32 * t233 - t258 * t234 - t266;
t286 = t316 * t293 * t270;
t283 = 2 * m(5);
t282 = 2 * m(6);
t281 = 2 * m(7);
t280 = -2 * mrSges(5,3);
t185 = t156 * mrSges(7,1) - t157 * mrSges(7,3);
t118 = t185 * qJD(5);
t279 = 0.2e1 * t118;
t227 = t270 * pkin(3);
t145 = -t227 - pkin(4);
t278 = 0.2e1 * t145;
t277 = 0.2e1 * qJD(2);
t276 = m(5) * pkin(3);
t275 = t84 / 0.2e1;
t264 = t157 * t6;
t37 = qJD(4) * t90 + t115 * t268 + t270 * t167;
t89 = t124 * t268 - t270 * t125;
t263 = t37 * t89;
t261 = t302 * pkin(9);
t253 = Ifges(6,6) * t156;
t251 = pkin(3) * qJD(4);
t250 = t116 * t37;
t249 = t116 * t83;
t248 = t117 * t84;
t243 = t157 * t84;
t242 = t84 * t156;
t237 = t307 * t144;
t236 = t307 * pkin(9);
t231 = 0.2e1 * mrSges(5,3) * t83;
t230 = t84 * t280;
t228 = m(6) * t268;
t224 = t270 * mrSges(5,2);
t222 = t268 * mrSges(5,1);
t220 = t268 * t89;
t147 = mrSges(7,2) * t233;
t214 = t156 * t270;
t213 = t157 * t270;
t212 = t268 * t128;
t202 = t234 / 0.2e1;
t173 = t117 * t198;
t199 = t302 * t144 + t293 * t173;
t197 = pkin(3) * t203;
t149 = Ifges(6,5) * t233;
t189 = -Ifges(6,6) * t234 + t149;
t187 = -t83 * t89 + t250;
t186 = t156 * mrSges(6,1) + t157 * mrSges(6,2);
t182 = -Ifges(6,2) * t156 + t256;
t181 = Ifges(7,3) * t156 + t254;
t180 = -t157 * pkin(5) - t156 * qJ(6);
t179 = pkin(5) * t156 - qJ(6) * t157;
t172 = t116 * t197;
t126 = -pkin(4) + t180;
t103 = pkin(5) * t234 - qJ(6) * t233 - t156 * qJD(6);
t166 = -t266 + (-t156 * t258 - t157 * t32) * qJD(5);
t164 = qJD(5) * t180 + t232;
t119 = t186 * qJD(5);
t163 = -t84 * mrSges(5,2) + (t118 + t119) * t116 + (mrSges(5,1) + t308) * t83 + t316 * t302;
t120 = t181 * qJD(5);
t121 = t182 * qJD(5);
t129 = -Ifges(7,3) * t157 + t255;
t130 = Ifges(6,2) * t157 + t257;
t162 = (t129 - t130) * t234 + t314 * t233 + (-t120 + t121) * t157 + t309 * t156;
t161 = m(7) * t232 + (m(7) * t180 - t308) * qJD(5);
t14 = Ifges(7,5) * t169 + t84 * Ifges(7,6) + Ifges(7,3) * t303;
t15 = Ifges(6,4) * t169 - Ifges(6,2) * t303 + t84 * Ifges(6,6);
t38 = -t116 * t179 + t89;
t43 = Ifges(7,6) * t117 - t116 * t181;
t44 = t117 * Ifges(6,6) - t116 * t182;
t9 = t116 * t164 + t179 * t83 + t37;
t160 = t43 * t202 + mrSges(6,3) * t264 - t44 * t234 / 0.2e1 - t36 * mrSges(5,2) + Ifges(5,5) * t83 - Ifges(5,6) * t84 + t38 * t118 + t89 * t119 + t9 * t127 + (-mrSges(5,1) + t128) * t37 + t312 * t275 + (t235 + t189) * t117 / 0.2e1 + t311 * t156 / 0.2e1 + (-t130 / 0.2e1 + t129 / 0.2e1) * t245 + (t121 / 0.2e1 - t120 / 0.2e1) * t241 + (-Ifges(7,6) * t275 + t15 / 0.2e1 - t14 / 0.2e1) * t157 - (qJD(5) * t129 + t309) * t240 / 0.2e1 + (t116 * t130 + t310) * t233 / 0.2e1 + t313 * mrSges(7,2) + t314 * (t116 * t202 + t244 / 0.2e1);
t159 = (-t156 * t260 + t157 * t259) * qJD(5) + m(6) * (t264 + t287) + m(7) * t313 + t301;
t108 = -t227 + t126;
t99 = t197 + t103;
t61 = t186 * t116;
t60 = t185 * t116;
t20 = mrSges(6,1) * t303 + mrSges(6,2) * t169;
t19 = mrSges(7,1) * t303 - mrSges(7,3) * t169;
t1 = [(mrSges(4,1) * t269 + mrSges(4,2) * t271 + mrSges(3,3)) * t277 + (-t14 + t15) * t241 - t311 * t240 + t169 * t310 + (t135 * t142 + t36 * t90 + t263) * t283 + t250 * t280 + (t2 * t26 + t27 * t4 + t38 * t9) * t281 + ((m(4) + m(3)) * t277 + 0.2e1 * (mrSges(4,1) * t271 - mrSges(4,2) * t269) * qJD(3)) * qJ(2) - t303 * (t44 - t43) + (t258 * t6 + t32 * t7 + t263) * t282 + 0.2e1 * t258 * t24 + 0.2e1 * (t116 * t84 - t117 * t83) * Ifges(5,4) + (-0.2e1 * t135 * mrSges(5,2) - 0.2e1 * Ifges(5,1) * t83 + (-Ifges(7,6) * t156 - t157 * t298 + t253) * t84) * t116 + (0.2e1 * Ifges(4,4) * t269 + t271 * t299) * t205 + (-0.2e1 * Ifges(4,4) * t271 + t269 * t299) * t206 + (0.2e1 * mrSges(5,1) * t135 - Ifges(6,6) * t303 + t280 * t36 + t298 * t211 + ((2 * Ifges(5,2)) + t315) * t84 + t291) * t117 + t90 * t230 + (t231 + 0.2e1 * t20) * t89 + 0.2e1 * t26 * t25 + 0.2e1 * t27 * t23 + 0.2e1 * t32 * t22 + 0.2e1 * t38 * t19 - 0.2e1 * t9 * t60 - 0.2e1 * t37 * t61 + 0.2e1 * t6 * t73 + 0.2e1 * t7 * t74 + 0.2e1 * t4 * t75 + 0.2e1 * t2 * t76 + 0.2e1 * t142 * (mrSges(5,1) * t84 + mrSges(5,2) * t83); (t60 + t61) * t83 + (t19 + t20 + t231) * t116 + (t156 * t259 + t157 * t260) * t84 + m(7) * (t116 * t9 + t242 * t27 + t243 * t26 - t38 * t83) + m(6) * (-t242 * t32 + t243 * t258 + t187) + m(5) * (t84 * t90 + t187) + (m(5) * t36 + t159 + t230) * t117; (t248 - t249) * t283 + 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (t248 * t293 - t249); t160 + (t268 * t36 - t270 * t37 + (t270 * t90 + t220) * qJD(4)) * t276 - t61 * t197 + m(6) * (t145 * t37 + (t213 * t258 - t214 * t32 + t220) * t251) + m(7) * (t108 * t9 + t99 * t38 + (t213 * t26 + t214 * t27) * t251) - t99 * t60 + t108 * t19 + t145 * t20 - Ifges(4,5) * t205 - Ifges(4,6) * t206 + t292 * t158 + t287 * mrSges(6,3) - t305 * t260 - t290 * t259 + (m(6) * (t166 + t264) + m(7) * ((-t156 * t26 + t157 * t27) * qJD(5) + t188) + t301) * t144 + (pkin(3) * t288 - t172 - t173) * mrSges(5,3); t163 + ((t116 * t268 + t117 * t270) * qJD(4) - t288) * t276 + m(6) * (-t145 * t83 + t172 + t199) + m(7) * (-t108 * t83 + t116 * t99 + t199) + t292; t108 * t279 + t119 * t278 + 0.2e1 * t99 * t127 + (t108 * t99 + t237) * t281 + t237 * t282 + (t228 * t278 + 0.2e1 * t212 - 0.2e1 * t222 - 0.2e1 * t224 + 0.2e1 * t286) * t251 + t162; t160 + t159 * pkin(9) + t166 * mrSges(6,3) + m(7) * (t103 * t38 + t126 * t9) - t103 * t60 + t126 * t19 + (-m(6) * t37 - t20) * pkin(4); m(7) * (t103 * t116 - t126 * t83 + t261) + m(6) * (pkin(4) * t83 + t261) + t163; (t103 + t99) * t127 + (t145 - pkin(4)) * t119 + (t126 + t108) * t118 + m(7) * (t103 * t108 + t126 * t99 + t236) + m(6) * t236 + (-pkin(4) * t228 + t212 - t222 - t224 + t286) * t251 + t162; -0.2e1 * pkin(4) * t119 + t126 * t279 + 0.2e1 * (m(7) * t126 + t127) * t103 + t162; -Ifges(6,6) * t245 - pkin(5) * t23 + m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t26) + qJD(6) * t76 + qJ(6) * t25 + t2 * mrSges(7,3) - t6 * mrSges(6,2) + t7 * mrSges(6,1) - t4 * mrSges(7,1) + t312 * t116 * qJD(5) + t291; (-m(7) * t179 - t185 - t186) * t84 + t161 * t117; m(7) * ((-pkin(5) * t214 + qJ(6) * t213) * t251 + t164 * t144) - pkin(5) * t147 - qJ(6) * mrSges(7,2) * t234 + t189 + t294 + (mrSges(7,1) + mrSges(6,1)) * t290 + (-mrSges(7,3) + mrSges(6,2)) * t305; t149 + (mrSges(7,2) * t180 - t253) * qJD(5) + t161 * pkin(9) + t294; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t4 + t23; (t117 * t233 + t242) * m(7); -m(7) * t290 + t147; m(7) * pkin(9) * t233 + t147; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
