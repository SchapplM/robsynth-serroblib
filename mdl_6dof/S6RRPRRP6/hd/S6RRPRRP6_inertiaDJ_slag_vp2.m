% Calculate time derivative of joint inertia matrix for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:50
% EndTime: 2019-03-09 12:08:05
% DurationCPUTime: 6.32s
% Computational Cost: add. (8542->634), mult. (23703->875), div. (0->0), fcn. (22934->10), ass. (0->261)
t207 = sin(qJ(5));
t331 = (Ifges(7,4) + Ifges(6,5)) * t207;
t204 = sin(pkin(6));
t330 = 0.2e1 * t204;
t210 = cos(qJ(5));
t270 = t207 ^ 2 + t210 ^ 2;
t208 = sin(qJ(4));
t211 = cos(qJ(4));
t265 = qJD(4) * t211;
t247 = t207 * t265;
t262 = qJD(5) * t210;
t218 = t208 * t262 + t247;
t205 = cos(pkin(11));
t194 = -pkin(2) * t205 - pkin(3);
t149 = -pkin(4) * t211 - pkin(10) * t208 + t194;
t203 = sin(pkin(11));
t193 = pkin(2) * t203 + pkin(9);
t285 = t193 * t211;
t325 = t207 * t149 + t210 * t285;
t328 = qJD(5) * t325;
t209 = sin(qJ(2));
t212 = cos(qJ(2));
t132 = (t203 * t212 + t205 * t209) * t204;
t129 = qJD(2) * t132;
t206 = cos(pkin(6));
t111 = t132 * t208 - t206 * t211;
t269 = qJD(2) * t204;
t282 = t205 * t212;
t130 = (-t203 * t209 + t282) * t269;
t85 = -qJD(4) * t111 + t130 * t211;
t112 = t132 * t211 + t206 * t208;
t284 = t204 * t209;
t131 = t203 * t284 - t204 * t282;
t91 = t112 * t210 + t131 * t207;
t42 = qJD(5) * t91 - t129 * t210 + t207 * t85;
t90 = t112 * t207 - t131 * t210;
t43 = -qJD(5) * t90 + t129 * t207 + t210 * t85;
t84 = qJD(4) * t112 + t130 * t208;
t7 = Ifges(6,5) * t43 - Ifges(6,6) * t42 + Ifges(6,3) * t84;
t8 = Ifges(7,4) * t43 + Ifges(7,2) * t84 + Ifges(7,6) * t42;
t327 = t7 + t8;
t311 = pkin(1) * t206;
t190 = t209 * t311;
t283 = t204 * t212;
t306 = -pkin(8) - qJ(3);
t324 = (t283 * t306 - t190) * qJD(2) - qJD(3) * t284;
t15 = mrSges(6,1) * t42 + mrSges(6,2) * t43;
t267 = qJD(4) * t208;
t191 = t212 * t311;
t181 = qJD(2) * t191;
t236 = t306 * t209;
t105 = t181 + (qJD(2) * t236 + qJD(3) * t212) * t204;
t67 = t205 * t105 + t203 * t324;
t114 = pkin(2) * t206 + t204 * t236 + t191;
t146 = pkin(8) * t283 + t190;
t127 = qJ(3) * t283 + t146;
t87 = t203 * t114 + t205 * t127;
t74 = pkin(9) * t206 + t87;
t252 = t209 * t269;
t233 = pkin(2) * t252;
t89 = pkin(3) * t129 - pkin(9) * t130 + t233;
t150 = (-pkin(2) * t212 - pkin(1)) * t204;
t92 = t131 * pkin(3) - t132 * pkin(9) + t150;
t22 = -t208 * t67 + t211 * t89 - t74 * t265 - t92 * t267;
t52 = t208 * t92 + t211 * t74;
t63 = mrSges(5,1) * t129 - mrSges(5,3) * t85;
t323 = m(5) * (qJD(4) * t52 + t22) - t15 + t63;
t167 = (pkin(4) * t208 - pkin(10) * t211) * qJD(4);
t322 = t167 * t210 - t328;
t21 = t208 * t89 + t211 * t67 + t92 * t265 - t267 * t74;
t17 = pkin(10) * t129 + t21;
t66 = t105 * t203 - t205 * t324;
t29 = pkin(4) * t84 - pkin(10) * t85 + t66;
t45 = pkin(10) * t131 + t52;
t86 = t114 * t205 - t203 * t127;
t73 = -pkin(3) * t206 - t86;
t50 = pkin(4) * t111 - pkin(10) * t112 + t73;
t304 = t207 * t50 + t210 * t45;
t4 = -qJD(5) * t304 - t17 * t207 + t210 * t29;
t225 = pkin(5) * t210 + qJ(6) * t207;
t261 = qJD(6) * t210;
t321 = qJD(5) * t225 - t261;
t320 = 2 * m(6);
t319 = 2 * m(7);
t318 = -2 * mrSges(3,3);
t317 = -2 * mrSges(4,3);
t316 = -2 * Ifges(4,4);
t315 = 0.2e1 * t150;
t314 = 0.2e1 * t193;
t313 = m(4) * pkin(2);
t18 = -pkin(4) * t129 - t22;
t312 = m(6) * t18;
t310 = t66 * mrSges(4,1);
t309 = t66 * mrSges(5,1);
t308 = t66 * mrSges(5,2);
t307 = t67 * mrSges(4,2);
t55 = mrSges(6,1) * t90 + mrSges(6,2) * t91;
t94 = mrSges(5,1) * t131 - mrSges(5,3) * t112;
t303 = t55 - t94;
t56 = -mrSges(7,2) * t90 + mrSges(7,3) * t111;
t57 = -mrSges(6,2) * t111 - mrSges(6,3) * t90;
t302 = t56 + t57;
t58 = mrSges(6,1) * t111 - mrSges(6,3) * t91;
t59 = -mrSges(7,1) * t111 + mrSges(7,2) * t91;
t301 = -t58 + t59;
t300 = Ifges(5,4) * t208;
t299 = Ifges(5,4) * t211;
t298 = Ifges(6,4) * t207;
t297 = Ifges(6,4) * t210;
t296 = Ifges(7,5) * t207;
t295 = Ifges(7,5) * t210;
t294 = Ifges(6,6) * t210;
t293 = t129 * Ifges(5,5);
t292 = t129 * Ifges(5,6);
t291 = t131 * Ifges(5,6);
t140 = -pkin(8) * t252 + t181;
t290 = t140 * mrSges(3,2);
t141 = t146 * qJD(2);
t289 = t141 * mrSges(3,1);
t170 = -t210 * mrSges(6,1) + t207 * mrSges(6,2);
t288 = t170 - mrSges(5,1);
t287 = t149 * t210;
t281 = t207 * t208;
t280 = t208 * t210;
t246 = t210 * t265;
t263 = qJD(5) * t208;
t250 = t207 * t263;
t219 = t246 - t250;
t123 = mrSges(6,1) * t267 - mrSges(6,3) * t219;
t264 = qJD(5) * t207;
t124 = mrSges(7,2) * t246 + (-mrSges(7,1) * qJD(4) - mrSges(7,2) * t264) * t208;
t279 = -t123 + t124;
t125 = -mrSges(6,2) * t267 - mrSges(6,3) * t218;
t126 = -mrSges(7,2) * t218 + mrSges(7,3) * t267;
t278 = t125 + t126;
t277 = t149 * t262 + t207 * t167;
t226 = Ifges(7,3) * t207 + t295;
t134 = -Ifges(7,6) * t211 + t208 * t226;
t227 = -Ifges(6,2) * t207 + t297;
t137 = -Ifges(6,6) * t211 + t208 * t227;
t276 = t134 - t137;
t228 = Ifges(7,1) * t210 + t296;
t138 = -Ifges(7,4) * t211 + t208 * t228;
t229 = Ifges(6,1) * t210 - t298;
t139 = -Ifges(6,5) * t211 + t208 * t229;
t275 = t138 + t139;
t162 = mrSges(6,2) * t211 - mrSges(6,3) * t281;
t165 = -mrSges(7,2) * t281 - mrSges(7,3) * t211;
t274 = t162 + t165;
t163 = -mrSges(6,1) * t211 - mrSges(6,3) * t280;
t164 = mrSges(7,1) * t211 + mrSges(7,2) * t280;
t273 = -t163 + t164;
t272 = Ifges(6,5) * t246 + Ifges(6,3) * t267;
t271 = t270 * pkin(10) * t265;
t156 = Ifges(7,4) * t262 + Ifges(7,6) * t264;
t268 = qJD(4) * t207;
t266 = qJD(4) * t210;
t6 = Ifges(7,5) * t43 + Ifges(7,6) * t84 + Ifges(7,3) * t42;
t9 = Ifges(6,4) * t43 - Ifges(6,2) * t42 + Ifges(6,6) * t84;
t260 = t6 / 0.2e1 - t9 / 0.2e1;
t259 = Ifges(5,5) * t85 - Ifges(5,6) * t84 + Ifges(5,3) * t129;
t10 = Ifges(7,1) * t43 + Ifges(7,4) * t84 + Ifges(7,5) * t42;
t11 = Ifges(6,1) * t43 - Ifges(6,4) * t42 + Ifges(6,5) * t84;
t257 = t10 / 0.2e1 + t11 / 0.2e1;
t30 = Ifges(7,5) * t91 + Ifges(7,6) * t111 + Ifges(7,3) * t90;
t33 = Ifges(6,4) * t91 - Ifges(6,2) * t90 + Ifges(6,6) * t111;
t256 = t30 / 0.2e1 - t33 / 0.2e1;
t34 = Ifges(7,1) * t91 + Ifges(7,4) * t111 + Ifges(7,5) * t90;
t35 = Ifges(6,1) * t91 - Ifges(6,4) * t90 + Ifges(6,5) * t111;
t255 = t34 / 0.2e1 + t35 / 0.2e1;
t171 = -Ifges(7,3) * t210 + t296;
t96 = -t171 * t263 + (Ifges(7,6) * t208 + t211 * t226) * qJD(4);
t174 = Ifges(6,2) * t210 + t298;
t99 = -t174 * t263 + (Ifges(6,6) * t208 + t211 * t227) * qJD(4);
t254 = -t99 / 0.2e1 + t96 / 0.2e1;
t253 = Ifges(3,5) * t212 * t269 + Ifges(4,5) * t130 - Ifges(4,6) * t129;
t251 = t193 * t267;
t176 = Ifges(7,1) * t207 - t295;
t100 = -t176 * t263 + (Ifges(7,4) * t208 + t211 * t228) * qJD(4);
t177 = Ifges(6,1) * t207 + t297;
t101 = -t177 * t263 + (Ifges(6,5) * t208 + t211 * t229) * qJD(4);
t245 = t100 / 0.2e1 + t101 / 0.2e1;
t244 = t134 / 0.2e1 - t137 / 0.2e1;
t243 = t138 / 0.2e1 + t139 / 0.2e1;
t154 = t226 * qJD(5);
t157 = t227 * qJD(5);
t242 = t154 / 0.2e1 - t157 / 0.2e1;
t155 = Ifges(6,5) * t262 - Ifges(6,6) * t264;
t241 = t155 / 0.2e1 + t156 / 0.2e1;
t159 = t228 * qJD(5);
t160 = t229 * qJD(5);
t240 = t159 / 0.2e1 + t160 / 0.2e1;
t239 = t171 / 0.2e1 - t174 / 0.2e1;
t238 = t294 / 0.2e1 - Ifges(7,6) * t210 / 0.2e1 + t331 / 0.2e1;
t237 = t176 / 0.2e1 + t177 / 0.2e1;
t26 = -t84 * mrSges(7,1) + t43 * mrSges(7,2);
t235 = t193 * t207 + pkin(5);
t51 = -t208 * t74 + t211 * t92;
t232 = Ifges(7,4) * t246 + Ifges(7,2) * t267 + Ifges(7,6) * t218;
t231 = t207 * mrSges(6,1) + t210 * mrSges(6,2);
t169 = -t210 * mrSges(7,1) - t207 * mrSges(7,3);
t230 = t207 * mrSges(7,1) - t210 * mrSges(7,3);
t224 = pkin(5) * t207 - qJ(6) * t210;
t19 = -t207 * t45 + t210 * t50;
t44 = -pkin(4) * t131 - t51;
t221 = t193 + t224;
t3 = t210 * t17 + t207 * t29 + t50 * t262 - t264 * t45;
t1 = qJ(6) * t84 + qJD(6) * t111 + t3;
t12 = qJ(6) * t111 + t304;
t13 = -pkin(5) * t111 - t19;
t2 = -pkin(5) * t84 - t4;
t216 = t1 * t210 - t12 * t264 + t13 * t262 + t2 * t207;
t215 = -t19 * t262 - t207 * t4 + t210 * t3 - t264 * t304;
t24 = -mrSges(6,2) * t84 - mrSges(6,3) * t42;
t25 = mrSges(6,1) * t84 - mrSges(6,3) * t43;
t27 = -mrSges(7,2) * t42 + mrSges(7,3) * t84;
t214 = (t24 + t27) * t210 + (-t25 + t26) * t207 + (-t207 * t302 + t210 * t301) * qJD(5);
t213 = m(7) * t261 + (-m(7) * t225 + t169 + t170) * qJD(5);
t199 = Ifges(5,5) * t265;
t178 = Ifges(5,1) * t208 + t299;
t175 = Ifges(5,2) * t211 + t300;
t168 = -pkin(4) - t225;
t161 = (Ifges(5,1) * t211 - t300) * qJD(4);
t158 = (-Ifges(5,2) * t208 + t299) * qJD(4);
t153 = (mrSges(5,1) * t208 + mrSges(5,2) * t211) * qJD(4);
t152 = t231 * qJD(5);
t151 = t230 * qJD(5);
t148 = t231 * t208;
t147 = t230 * t208;
t145 = -pkin(8) * t284 + t191;
t143 = qJD(5) * t224 - qJD(6) * t207;
t136 = -Ifges(7,2) * t211 + (Ifges(7,4) * t210 + Ifges(7,6) * t207) * t208;
t135 = -Ifges(6,3) * t211 + (Ifges(6,5) * t210 - Ifges(6,6) * t207) * t208;
t119 = t130 * mrSges(4,2);
t118 = t221 * t208;
t109 = -t207 * t285 + t287;
t107 = mrSges(6,1) * t218 + mrSges(6,2) * t219;
t106 = mrSges(7,1) * t218 - mrSges(7,3) * t219;
t104 = t211 * t235 - t287;
t103 = -qJ(6) * t211 + t325;
t98 = -Ifges(7,4) * t250 + t232;
t97 = -Ifges(6,5) * t250 - Ifges(6,6) * t218 + t272;
t95 = t208 * t321 + t221 * t265;
t93 = -mrSges(5,2) * t131 - mrSges(5,3) * t111;
t76 = t207 * t251 + t322;
t75 = (-t208 * t266 - t211 * t264) * t193 + t277;
t69 = -t235 * t267 - t322;
t68 = (-t193 * t264 - qJD(6)) * t211 + (-t193 * t210 + qJ(6)) * t267 + t277;
t62 = -mrSges(5,2) * t129 - mrSges(5,3) * t84;
t61 = Ifges(5,1) * t112 - Ifges(5,4) * t111 + Ifges(5,5) * t131;
t60 = Ifges(5,4) * t112 - Ifges(5,2) * t111 + t291;
t54 = mrSges(7,1) * t90 - mrSges(7,3) * t91;
t53 = mrSges(5,1) * t84 + mrSges(5,2) * t85;
t47 = Ifges(5,1) * t85 - Ifges(5,4) * t84 + t293;
t46 = Ifges(5,4) * t85 - Ifges(5,2) * t84 + t292;
t32 = Ifges(7,4) * t91 + Ifges(7,2) * t111 + Ifges(7,6) * t90;
t31 = Ifges(6,5) * t91 - Ifges(6,6) * t90 + Ifges(6,3) * t111;
t23 = pkin(5) * t90 - qJ(6) * t91 + t44;
t14 = mrSges(7,1) * t42 - mrSges(7,3) * t43;
t5 = pkin(5) * t42 - qJ(6) * t43 - qJD(6) * t91 + t18;
t16 = [(0.2e1 * (t140 * t212 + t141 * t209) * mrSges(3,3) + ((t145 * t318 + Ifges(3,5) * t206 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t212) * t330) * t212 + (0.2e1 * pkin(2) * (mrSges(4,1) * t131 + mrSges(4,2) * t132) + t146 * t318 + t313 * t315 - 0.2e1 * Ifges(3,6) * t206 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t209 + (Ifges(3,1) - Ifges(3,2)) * t212) * t330) * t209) * qJD(2)) * t204 + (t30 - t33) * t42 + 0.2e1 * (-t131 * t67 + t132 * t66) * mrSges(4,3) + 0.2e1 * m(4) * (-t66 * t86 + t67 * t87) + 0.2e1 * t21 * t93 + 0.2e1 * t22 * t94 + t85 * t61 + 0.2e1 * t73 * t53 + 0.2e1 * t5 * t54 + 0.2e1 * t18 * t55 + 0.2e1 * t1 * t56 + 0.2e1 * t3 * t57 + 0.2e1 * t4 * t58 + 0.2e1 * t2 * t59 + 0.2e1 * t52 * t62 + 0.2e1 * t51 * t63 + 0.2e1 * t44 * t15 + 0.2e1 * t19 * t25 + 0.2e1 * t13 * t26 + 0.2e1 * t12 * t27 + 0.2e1 * t23 * t14 + (-t46 + 0.2e1 * t309 + t327) * t111 + t119 * t315 + (t1 * t12 + t13 * t2 + t23 * t5) * t319 + 0.2e1 * t304 * t24 + (t18 * t44 + t19 * t4 + t3 * t304) * t320 + 0.2e1 * m(3) * (t140 * t146 - t141 * t145) + (t253 - 0.2e1 * t289 - 0.2e1 * t290 - 0.2e1 * t307 - 0.2e1 * t310) * t206 + (t34 + t35) * t43 + t131 * t259 + (-t60 + t31 + t32) * t84 + (t47 + 0.2e1 * t308) * t112 + (mrSges(4,1) * t315 + t87 * t317 + t132 * t316 + Ifges(5,5) * t112 - Ifges(4,6) * t206 - Ifges(5,6) * t111 + ((2 * Ifges(4,2)) + Ifges(5,3)) * t131) * t129 + (0.2e1 * Ifges(4,1) * t132 + Ifges(4,5) * t206 + t131 * t316 + t317 * t86) * t130 + (t6 - t9) * t90 + (t11 + t10) * t91 + 0.2e1 * m(5) * (t21 * t52 + t22 * t51 + t66 * t73); (-t22 * mrSges(5,3) + t293 / 0.2e1 + t308 + t47 / 0.2e1 + t257 * t210 + t260 * t207 + (-t207 * t255 + t210 * t256) * qJD(5) + (-t60 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1 - t52 * mrSges(5,3) - t291 / 0.2e1) * qJD(4) + (-qJD(4) * t93 + t312 - t323) * t193) * t208 + t253 + (-t158 / 0.2e1 + t97 / 0.2e1 + t98 / 0.2e1) * t111 - t289 - t290 - t307 + t131 * t199 / 0.2e1 + (-t129 * t203 - t130 * t205) * pkin(2) * mrSges(4,3) + m(7) * (t1 * t103 + t104 * t2 + t118 * t5 + t12 * t68 + t13 * t69 + t23 * t95) + t85 * t178 / 0.2e1 + t73 * t153 + t112 * t161 / 0.2e1 + t3 * t162 + t4 * t163 + t2 * t164 + t1 * t165 + t5 * t147 + t18 * t148 + t19 * t123 + t13 * t124 + t12 * t126 + t118 * t14 + t109 * t25 + t103 * t27 + t104 * t26 + t23 * t106 + t44 * t107 + t95 * t54 + t69 * t59 + t75 * t57 + t76 * t58 + t68 * t56 + (-t175 / 0.2e1 + t135 / 0.2e1 + t136 / 0.2e1) * t84 - t310 + (t203 * t67 - t205 * t66) * t313 + t243 * t43 + t244 * t42 + t245 * t91 + t304 * t125 - Ifges(3,6) * t252 + (m(5) * t66 + t53) * t194 + t254 * t90 + (t21 * mrSges(5,3) + t292 / 0.2e1 - t309 + t46 / 0.2e1 - t7 / 0.2e1 - t8 / 0.2e1 + (m(5) * t21 + t62) * t193 + (t61 / 0.2e1 - t51 * mrSges(5,3) + t255 * t210 + t256 * t207 + (-m(5) * t51 + m(6) * t44 + t303) * t193) * qJD(4)) * t211 + t325 * t24 + m(6) * (t109 * t4 + t19 * t76 + t3 * t325 + t304 * t75); 0.2e1 * t103 * t126 + 0.2e1 * t104 * t124 + 0.2e1 * t118 * t106 + 0.2e1 * t109 * t123 + 0.2e1 * t325 * t125 + 0.2e1 * t95 * t147 + 0.2e1 * t194 * t153 + 0.2e1 * t75 * t162 + 0.2e1 * t76 * t163 + 0.2e1 * t69 * t164 + 0.2e1 * t68 * t165 + (t109 * t76 + t325 * t75) * t320 + (t103 * t68 + t104 * t69 + t118 * t95) * t319 + (t158 - t97 - t98 + (t148 * t314 + t207 * t276 + t210 * t275 + t178) * qJD(4)) * t211 + (t107 * t314 + t161 + (t100 + t101) * t210 + (t96 - t99) * t207 + (-t207 * t275 + t210 * t276) * qJD(5) + (t193 ^ 2 * t211 * t320 + t135 + t136 - t175) * qJD(4)) * t208; m(4) * t233 + t129 * mrSges(4,1) + t119 + (-t14 + (t207 * t301 + t210 * t302 + t93) * qJD(4) + m(6) * (-t19 * t268 + t266 * t304 - t18) + m(7) * (t12 * t266 + t13 * t268 - t5) + t323) * t211 + (t62 + (t54 + t303) * qJD(4) + m(6) * (qJD(4) * t44 + t215) + m(7) * (qJD(4) * t23 + t216) + m(5) * (-qJD(4) * t51 + t21) + t214) * t208; (-t107 - m(7) * t95 - t106 + (t274 * t210 + t273 * t207 + m(7) * (t103 * t210 + t104 * t207) + (-t109 * t207 + t210 * t325 - t285) * m(6)) * qJD(4)) * t211 + (t278 * t210 + t279 * t207 + (t147 + t148) * qJD(4) + (-t207 * t274 + t210 * t273) * qJD(5) + m(7) * (qJD(4) * t118 - t103 * t264 + t104 * t262 + t207 * t69 + t210 * t68) + (-t109 * t262 - t207 * t76 + t210 * t75 - t264 * t325 + t251) * m(6)) * t208; 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (-0.1e1 + t270) * t208 * t265; (m(6) * t215 + m(7) * t216 + t214) * pkin(10) + t259 + (t1 * mrSges(7,2) + t3 * mrSges(6,3) - t260) * t210 + t44 * t152 + t168 * t14 + t5 * t169 + t18 * t170 + t143 * t54 + t23 * t151 + m(7) * (t143 * t23 + t168 * t5) - t21 * mrSges(5,2) + t22 * mrSges(5,1) + t237 * t43 + t238 * t84 + t239 * t42 + t240 * t91 + t241 * t111 + t242 * t90 + ((t13 * mrSges(7,2) - t19 * mrSges(6,3) + t255) * t210 + (-t12 * mrSges(7,2) - mrSges(6,3) * t304 + t256) * t207) * qJD(5) + (t2 * mrSges(7,2) - t4 * mrSges(6,3) + t257) * t207 + (-t312 - t15) * pkin(4); t208 * t193 * t152 + t199 + m(7) * (t118 * t143 + t168 * t95) + t168 * t106 + t95 * t169 + t143 * t147 + t118 * t151 - pkin(4) * t107 - t241 * t211 + ((-Ifges(5,6) + t238) * t208 + (t208 * mrSges(5,2) + (-m(6) * pkin(4) + t288) * t211) * t193) * qJD(4) + (t75 * mrSges(6,3) + t68 * mrSges(7,2) + t240 * t208 + t237 * t265 + (t104 * mrSges(7,2) - t109 * mrSges(6,3) + t208 * t239 + t243) * qJD(5) + (t273 * qJD(5) + m(6) * (-qJD(5) * t109 + t75) + m(7) * (qJD(5) * t104 + t68) + t278) * pkin(10) - t254) * t210 + (-t76 * mrSges(6,3) + t69 * mrSges(7,2) + t242 * t208 + t239 * t265 + (-t103 * mrSges(7,2) - mrSges(6,3) * t325 - t208 * t237 + t244) * qJD(5) + (-t274 * qJD(5) + m(6) * (-t76 - t328) + m(7) * (-qJD(5) * t103 + t69) + t279) * pkin(10) + t245) * t207; (t169 + t288) * t267 + m(6) * (-pkin(4) * t267 + t271) + m(7) * (t168 * t267 + t271) + (-t152 - t151 - m(7) * t143 + (-mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t270) * qJD(4)) * t211; -0.2e1 * pkin(4) * t152 + 0.2e1 * t151 * t168 + (-t154 + t157) * t210 + (t159 + t160) * t207 + 0.2e1 * (m(7) * t168 + t169) * t143 + ((t176 + t177) * t210 + (t171 - t174) * t207) * qJD(5); m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) + t1 * mrSges(7,3) + qJD(6) * t56 + qJ(6) * t27 - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - pkin(5) * t26 + t327; -Ifges(6,6) * t247 - pkin(5) * t124 + m(7) * (-pkin(5) * t69 + qJ(6) * t68 + qJD(6) * t103) + qJD(6) * t165 + qJ(6) * t126 + t68 * mrSges(7,3) + t76 * mrSges(6,1) - t69 * mrSges(7,1) - t75 * mrSges(6,2) + (-t294 - t331) * t263 + t232 + t272; (-m(7) * t224 - t230 - t231) * t265 + t213 * t208; -mrSges(7,2) * t321 + t213 * pkin(10) + t155 + t156; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t26; m(7) * t69 + t124; t218 * m(7); (m(7) * pkin(10) + mrSges(7,2)) * t262; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
