% Calculate time derivative of joint inertia matrix for
% S6RRPRRP5
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:52
% EndTime: 2019-03-09 11:58:10
% DurationCPUTime: 7.56s
% Computational Cost: add. (8508->662), mult. (23681->910), div. (0->0), fcn. (22987->10), ass. (0->267)
t215 = cos(qJ(5));
t330 = Ifges(6,6) + Ifges(7,6);
t334 = t330 * t215;
t212 = sin(qJ(5));
t333 = (Ifges(6,5) + Ifges(7,5)) * t212;
t209 = sin(pkin(6));
t332 = 0.2e1 * t209;
t216 = cos(qJ(4));
t260 = qJD(4) * t216;
t240 = t212 * t260;
t213 = sin(qJ(4));
t257 = qJD(5) * t215;
t242 = t213 * t257;
t219 = t240 + t242;
t208 = sin(pkin(11));
t210 = cos(pkin(11));
t214 = sin(qJ(2));
t217 = cos(qJ(2));
t131 = (t208 * t217 + t210 * t214) * t209;
t128 = qJD(2) * t131;
t211 = cos(pkin(6));
t112 = t131 * t213 - t211 * t216;
t264 = qJD(2) * t209;
t278 = t210 * t217;
t129 = (-t208 * t214 + t278) * t264;
t87 = -qJD(4) * t112 + t129 * t216;
t113 = t131 * t216 + t211 * t213;
t280 = t209 * t214;
t130 = t208 * t280 - t209 * t278;
t93 = t113 * t215 + t130 * t212;
t43 = -qJD(5) * t93 + t128 * t215 - t212 * t87;
t92 = -t113 * t212 + t130 * t215;
t44 = qJD(5) * t92 + t128 * t212 + t215 * t87;
t86 = qJD(4) * t113 + t129 * t213;
t7 = Ifges(7,5) * t44 + Ifges(7,6) * t43 + Ifges(7,3) * t86;
t8 = Ifges(6,5) * t44 + Ifges(6,6) * t43 + Ifges(6,3) * t86;
t329 = t7 + t8;
t325 = -m(6) * pkin(10) - mrSges(6,3);
t258 = qJD(5) * t213;
t244 = t212 * t258;
t245 = t215 * t260;
t248 = -t219 * mrSges(7,1) - mrSges(7,2) * t245;
t107 = -mrSges(7,2) * t244 - t248;
t196 = pkin(2) * t208 + pkin(9);
t115 = pkin(5) * t219 + t196 * t260;
t324 = m(7) * t115 + t107;
t197 = -pkin(2) * t210 - pkin(3);
t152 = -pkin(4) * t216 - pkin(10) * t213 + t197;
t274 = t215 * t216;
t170 = t196 * t274;
t111 = t212 * t152 + t170;
t265 = t212 ^ 2 + t215 ^ 2;
t307 = pkin(1) * t211;
t193 = t214 * t307;
t279 = t209 * t217;
t299 = -pkin(8) - qJ(3);
t323 = (t279 * t299 - t193) * qJD(2) - qJD(3) * t280;
t15 = -mrSges(6,1) * t43 + mrSges(6,2) * t44;
t262 = qJD(4) * t213;
t194 = t217 * t307;
t186 = qJD(2) * t194;
t229 = t214 * t299;
t106 = t186 + (qJD(2) * t229 + qJD(3) * t217) * t209;
t71 = t210 * t106 + t208 * t323;
t116 = pkin(2) * t211 + t209 * t229 + t194;
t149 = pkin(8) * t279 + t193;
t127 = qJ(3) * t279 + t149;
t89 = t208 * t116 + t210 * t127;
t76 = pkin(9) * t211 + t89;
t247 = t214 * t264;
t226 = pkin(2) * t247;
t91 = pkin(3) * t128 - pkin(9) * t129 + t226;
t154 = (-pkin(2) * t217 - pkin(1)) * t209;
t94 = t130 * pkin(3) - t131 * pkin(9) + t154;
t22 = -t213 * t71 + t216 * t91 - t76 * t260 - t94 * t262;
t54 = t213 * t94 + t216 * t76;
t67 = mrSges(5,1) * t128 - mrSges(5,3) * t87;
t322 = m(5) * (qJD(4) * t54 + t22) - t15 + t67;
t321 = 2 * m(6);
t320 = 0.2e1 * m(7);
t319 = -2 * mrSges(3,3);
t318 = -2 * mrSges(4,3);
t317 = -2 * mrSges(7,3);
t316 = -2 * Ifges(4,4);
t315 = 0.2e1 * t154;
t314 = 0.2e1 * t196;
t313 = m(4) * pkin(2);
t311 = m(7) * pkin(5);
t18 = -pkin(4) * t128 - t22;
t310 = m(6) * t18;
t309 = m(6) * t216;
t306 = pkin(4) * t213;
t305 = pkin(10) * t216;
t70 = t106 * t208 - t210 * t323;
t304 = t70 * mrSges(4,1);
t303 = t70 * mrSges(5,1);
t302 = t70 * mrSges(5,2);
t301 = t71 * mrSges(4,2);
t298 = -qJ(6) - pkin(10);
t46 = pkin(10) * t130 + t54;
t88 = t116 * t210 - t208 * t127;
t75 = -pkin(3) * t211 - t88;
t52 = pkin(4) * t112 - pkin(10) * t113 + t75;
t20 = t212 * t52 + t215 * t46;
t57 = -mrSges(6,1) * t92 + mrSges(6,2) * t93;
t96 = mrSges(5,1) * t130 - mrSges(5,3) * t113;
t296 = t57 - t96;
t59 = -mrSges(7,2) * t112 + mrSges(7,3) * t92;
t60 = -mrSges(6,2) * t112 + mrSges(6,3) * t92;
t295 = t59 + t60;
t61 = mrSges(7,1) * t112 - mrSges(7,3) * t93;
t62 = mrSges(6,1) * t112 - mrSges(6,3) * t93;
t294 = -t61 - t62;
t293 = mrSges(6,2) * t215;
t292 = Ifges(5,4) * t213;
t291 = Ifges(5,4) * t216;
t290 = Ifges(6,4) * t212;
t289 = Ifges(6,4) * t215;
t288 = Ifges(7,4) * t212;
t287 = Ifges(7,4) * t215;
t286 = t128 * Ifges(5,5);
t285 = t128 * Ifges(5,6);
t284 = t130 * Ifges(5,6);
t139 = -pkin(8) * t247 + t186;
t283 = t139 * mrSges(3,2);
t140 = t149 * qJD(2);
t282 = t140 * mrSges(3,1);
t174 = -mrSges(6,1) * t215 + mrSges(6,2) * t212;
t281 = t174 - mrSges(5,1);
t277 = t212 * t213;
t276 = t212 * t216;
t275 = t213 * t215;
t171 = (-t305 + t306) * qJD(4);
t273 = t152 * t257 + t212 * t171;
t221 = -Ifges(7,2) * t212 + t287;
t135 = -Ifges(7,6) * t216 + t213 * t221;
t222 = -Ifges(6,2) * t212 + t289;
t136 = -Ifges(6,6) * t216 + t213 * t222;
t272 = -t135 - t136;
t223 = Ifges(7,1) * t215 - t288;
t137 = -Ifges(7,5) * t216 + t213 * t223;
t224 = Ifges(6,1) * t215 - t290;
t138 = -Ifges(6,5) * t216 + t213 * t224;
t271 = t137 + t138;
t246 = t196 * t262;
t270 = t215 * t171 + t212 * t246;
t166 = mrSges(7,2) * t216 - mrSges(7,3) * t277;
t167 = mrSges(6,2) * t216 - mrSges(6,3) * t277;
t269 = t166 + t167;
t168 = -mrSges(7,1) * t216 - mrSges(7,3) * t275;
t169 = -mrSges(6,1) * t216 - mrSges(6,3) * t275;
t268 = -t168 - t169;
t267 = Ifges(7,5) * t245 + Ifges(7,3) * t262;
t266 = Ifges(6,5) * t245 + Ifges(6,3) * t262;
t259 = qJD(5) * t212;
t155 = mrSges(7,1) * t259 + mrSges(7,2) * t257;
t263 = qJD(4) * t212;
t261 = qJD(4) * t215;
t256 = qJD(6) * t215;
t10 = Ifges(6,4) * t44 + Ifges(6,2) * t43 + Ifges(6,6) * t86;
t9 = Ifges(7,4) * t44 + Ifges(7,2) * t43 + Ifges(7,6) * t86;
t255 = -t9 / 0.2e1 - t10 / 0.2e1;
t254 = Ifges(5,5) * t87 - Ifges(5,6) * t86 + Ifges(5,3) * t128;
t11 = Ifges(7,1) * t44 + Ifges(7,4) * t43 + Ifges(7,5) * t86;
t12 = Ifges(6,1) * t44 + Ifges(6,4) * t43 + Ifges(6,5) * t86;
t253 = t11 / 0.2e1 + t12 / 0.2e1;
t33 = Ifges(7,4) * t93 + Ifges(7,2) * t92 + Ifges(7,6) * t112;
t34 = Ifges(6,4) * t93 + Ifges(6,2) * t92 + Ifges(6,6) * t112;
t252 = -t33 / 0.2e1 - t34 / 0.2e1;
t35 = Ifges(7,1) * t93 + Ifges(7,4) * t92 + Ifges(7,5) * t112;
t36 = Ifges(6,1) * t93 + Ifges(6,4) * t92 + Ifges(6,5) * t112;
t251 = t35 / 0.2e1 + t36 / 0.2e1;
t250 = -mrSges(6,1) - t311;
t249 = Ifges(3,5) * t217 * t264 + Ifges(4,5) * t129 - Ifges(4,6) * t128;
t243 = t216 * t259;
t178 = Ifges(7,2) * t215 + t288;
t100 = -t178 * t258 + (Ifges(7,6) * t213 + t216 * t221) * qJD(4);
t179 = Ifges(6,2) * t215 + t290;
t101 = -t179 * t258 + (Ifges(6,6) * t213 + t216 * t222) * qJD(4);
t239 = t100 / 0.2e1 + t101 / 0.2e1;
t181 = Ifges(7,1) * t212 + t287;
t102 = -t181 * t258 + (Ifges(7,5) * t213 + t216 * t223) * qJD(4);
t182 = Ifges(6,1) * t212 + t289;
t103 = -t182 * t258 + (Ifges(6,5) * t213 + t216 * t224) * qJD(4);
t238 = t102 / 0.2e1 + t103 / 0.2e1;
t237 = t135 / 0.2e1 + t136 / 0.2e1;
t236 = t137 / 0.2e1 + t138 / 0.2e1;
t203 = Ifges(7,5) * t257;
t204 = Ifges(6,5) * t257;
t235 = t203 / 0.2e1 + t204 / 0.2e1 - t330 * t259 / 0.2e1;
t160 = t221 * qJD(5);
t161 = t222 * qJD(5);
t234 = t160 / 0.2e1 + t161 / 0.2e1;
t163 = t223 * qJD(5);
t164 = t224 * qJD(5);
t233 = t163 / 0.2e1 + t164 / 0.2e1;
t232 = t333 / 0.2e1 + t334 / 0.2e1;
t231 = t178 / 0.2e1 + t179 / 0.2e1;
t230 = t181 / 0.2e1 + t182 / 0.2e1;
t14 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t19 = -t212 * t46 + t215 * t52;
t53 = -t213 * t76 + t216 * t94;
t228 = qJD(5) * t298;
t225 = mrSges(6,1) * t212 + t293;
t45 = -pkin(4) * t130 - t53;
t21 = t213 * t91 + t216 * t71 + t94 * t260 - t262 * t76;
t17 = pkin(10) * t128 + t21;
t30 = pkin(4) * t86 - pkin(10) * t87 + t70;
t3 = t215 * t17 + t212 * t30 + t52 * t257 - t259 * t46;
t220 = -t244 + t245;
t4 = -qJD(5) * t20 - t17 * t212 + t215 * t30;
t205 = Ifges(5,5) * t260;
t198 = -pkin(5) * t215 - pkin(4);
t183 = Ifges(5,1) * t213 + t291;
t180 = Ifges(5,2) * t216 + t292;
t175 = t298 * t215;
t173 = -mrSges(7,1) * t215 + mrSges(7,2) * t212;
t172 = t298 * t212;
t165 = (Ifges(5,1) * t216 - t292) * qJD(4);
t162 = (-Ifges(5,2) * t213 + t291) * qJD(4);
t157 = (t213 * mrSges(5,1) + t216 * mrSges(5,2)) * qJD(4);
t156 = t225 * qJD(5);
t151 = t225 * t213;
t150 = (mrSges(7,1) * t212 + mrSges(7,2) * t215) * t213;
t148 = -pkin(8) * t280 + t194;
t145 = (pkin(5) * t212 + t196) * t213;
t144 = -qJD(6) * t212 + t215 * t228;
t143 = t212 * t228 + t256;
t142 = t215 * t152;
t134 = -Ifges(6,3) * t216 + (Ifges(6,5) * t215 - Ifges(6,6) * t212) * t213;
t133 = -Ifges(7,3) * t216 + (Ifges(7,5) * t215 - Ifges(7,6) * t212) * t213;
t126 = -mrSges(6,2) * t262 - mrSges(6,3) * t219;
t125 = -mrSges(7,2) * t262 - mrSges(7,3) * t219;
t124 = mrSges(6,1) * t262 - mrSges(6,3) * t220;
t123 = mrSges(7,1) * t262 - mrSges(7,3) * t220;
t119 = t129 * mrSges(4,2);
t110 = -t196 * t276 + t142;
t108 = mrSges(6,1) * t219 + mrSges(6,2) * t220;
t105 = -qJ(6) * t277 + t111;
t99 = -Ifges(6,5) * t244 - Ifges(6,6) * t219 + t266;
t98 = -Ifges(7,5) * t244 - Ifges(7,6) * t219 + t267;
t97 = -qJ(6) * t275 + t142 + (-t196 * t212 - pkin(5)) * t216;
t95 = -mrSges(5,2) * t130 - mrSges(5,3) * t112;
t78 = -qJD(5) * t111 + t270;
t77 = (-t213 * t261 - t243) * t196 + t273;
t66 = -mrSges(5,2) * t128 - mrSges(5,3) * t86;
t65 = Ifges(5,1) * t113 - Ifges(5,4) * t112 + Ifges(5,5) * t130;
t64 = Ifges(5,4) * t113 - Ifges(5,2) * t112 + t284;
t63 = (-qJ(6) * qJD(5) - qJD(4) * t196) * t275 + (-qJD(6) * t213 + (-qJ(6) * qJD(4) - qJD(5) * t196) * t216) * t212 + t273;
t58 = -t213 * t256 + (pkin(5) * t213 - qJ(6) * t274) * qJD(4) + (-t170 + (qJ(6) * t213 - t152) * t212) * qJD(5) + t270;
t56 = -mrSges(7,1) * t92 + mrSges(7,2) * t93;
t55 = mrSges(5,1) * t86 + mrSges(5,2) * t87;
t48 = Ifges(5,1) * t87 - Ifges(5,4) * t86 + t286;
t47 = Ifges(5,4) * t87 - Ifges(5,2) * t86 + t285;
t32 = Ifges(6,5) * t93 + Ifges(6,6) * t92 + Ifges(6,3) * t112;
t31 = Ifges(7,5) * t93 + Ifges(7,6) * t92 + Ifges(7,3) * t112;
t27 = -pkin(5) * t92 + t45;
t26 = mrSges(6,1) * t86 - mrSges(6,3) * t44;
t25 = mrSges(7,1) * t86 - mrSges(7,3) * t44;
t24 = -mrSges(6,2) * t86 + mrSges(6,3) * t43;
t23 = -mrSges(7,2) * t86 + mrSges(7,3) * t43;
t13 = qJ(6) * t92 + t20;
t6 = pkin(5) * t112 - qJ(6) * t93 + t19;
t5 = -pkin(5) * t43 + t18;
t2 = qJ(6) * t43 + qJD(6) * t92 + t3;
t1 = pkin(5) * t86 - qJ(6) * t44 - qJD(6) * t93 + t4;
t16 = [(t249 - 0.2e1 * t282 - 0.2e1 * t283 - 0.2e1 * t301 - 0.2e1 * t304) * t211 + 0.2e1 * m(3) * (t139 * t149 - t140 * t148) + 0.2e1 * (-t130 * t71 + t131 * t70) * mrSges(4,3) + (t33 + t34) * t43 + 0.2e1 * m(4) * (-t70 * t88 + t71 * t89) + (t35 + t36) * t44 + (mrSges(4,1) * t315 + t89 * t318 + t131 * t316 + Ifges(5,5) * t113 - Ifges(4,6) * t211 - Ifges(5,6) * t112 + ((2 * Ifges(4,2)) + Ifges(5,3)) * t130) * t128 + (0.2e1 * Ifges(4,1) * t131 + Ifges(4,5) * t211 + t130 * t316 + t318 * t88) * t129 + (t48 + 0.2e1 * t302) * t113 + t130 * t254 + 0.2e1 * t21 * t95 + 0.2e1 * t22 * t96 + t87 * t65 + 0.2e1 * t75 * t55 + 0.2e1 * t54 * t66 + 0.2e1 * t53 * t67 + 0.2e1 * t2 * t59 + 0.2e1 * t3 * t60 + 0.2e1 * t1 * t61 + 0.2e1 * t4 * t62 + 0.2e1 * t5 * t56 + 0.2e1 * t18 * t57 + 0.2e1 * t45 * t15 + 0.2e1 * t27 * t14 + 0.2e1 * t13 * t23 + 0.2e1 * t20 * t24 + 0.2e1 * t6 * t25 + 0.2e1 * t19 * t26 + 0.2e1 * m(5) * (t21 * t54 + t22 * t53 + t70 * t75) + (-t47 + 0.2e1 * t303 + t329) * t112 + (0.2e1 * (t139 * t217 + t140 * t214) * mrSges(3,3) + ((t148 * t319 + Ifges(3,5) * t211 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t217) * t332) * t217 + (t313 * t315 - 0.2e1 * Ifges(3,6) * t211 + 0.2e1 * pkin(2) * (mrSges(4,1) * t130 + mrSges(4,2) * t131) + t149 * t319 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t214 + (Ifges(3,1) - Ifges(3,2)) * t217) * t332) * t214) * qJD(2)) * t209 + t119 * t315 + (t1 * t6 + t13 * t2 + t27 * t5) * t320 + (t18 * t45 + t19 * t4 + t20 * t3) * t321 + (t31 + t32 - t64) * t86 + (t9 + t10) * t92 + (t11 + t12) * t93; (m(5) * t70 + t55) * t197 + (-t22 * mrSges(5,3) + t286 / 0.2e1 + t302 + t48 / 0.2e1 + t253 * t215 + t255 * t212 + (-t212 * t251 + t215 * t252) * qJD(5) + (-t64 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1 - t284 / 0.2e1 - t54 * mrSges(5,3)) * qJD(4) + (-qJD(4) * t95 + t310 - t322) * t196) * t213 + t249 + m(6) * (t110 * t4 + t111 * t3 + t19 * t78 + t20 * t77) + (-t162 / 0.2e1 + t98 / 0.2e1 + t99 / 0.2e1) * t112 + (-t128 * t208 - t129 * t210) * pkin(2) * mrSges(4,3) - t304 + (t21 * mrSges(5,3) + t285 / 0.2e1 - t303 + t47 / 0.2e1 - t8 / 0.2e1 - t7 / 0.2e1 + (m(5) * t21 + t66) * t196 + (t65 / 0.2e1 - t53 * mrSges(5,3) + t251 * t215 + t252 * t212 + (-m(5) * t53 + m(6) * t45 + t296) * t196) * qJD(4)) * t216 - t301 + t130 * t205 / 0.2e1 + t87 * t183 / 0.2e1 + m(7) * (t1 * t97 + t105 * t2 + t115 * t27 + t13 * t63 + t145 * t5 + t58 * t6) + t113 * t165 / 0.2e1 + t2 * t166 + t3 * t167 + t1 * t168 + t4 * t169 + t5 * t150 + t18 * t151 + t75 * t157 - Ifges(3,6) * t247 + t145 * t14 + t6 * t123 + t19 * t124 + t13 * t125 + t20 * t126 + t115 * t56 + t110 * t26 + t111 * t24 + t97 * t25 + t236 * t44 + t237 * t43 + t238 * t93 + t105 * t23 + t239 * t92 + t27 * t107 + t45 * t108 + t77 * t60 + t78 * t62 + t58 * t61 + t63 * t59 + (-t180 / 0.2e1 + t133 / 0.2e1 + t134 / 0.2e1) * t86 - t282 - t283 + (t208 * t71 - t210 * t70) * t313; 0.2e1 * t105 * t125 + 0.2e1 * t145 * t107 + 0.2e1 * t110 * t124 + 0.2e1 * t111 * t126 + 0.2e1 * t115 * t150 + 0.2e1 * t97 * t123 + 0.2e1 * t197 * t157 + 0.2e1 * t63 * t166 + 0.2e1 * t77 * t167 + 0.2e1 * t58 * t168 + 0.2e1 * t78 * t169 + (t110 * t78 + t111 * t77) * t321 + (t105 * t63 + t115 * t145 + t58 * t97) * t320 + (t162 - t98 - t99 + (t151 * t314 + t212 * t272 + t215 * t271 + t183) * qJD(4)) * t216 + (t108 * t314 + t165 + (t102 + t103) * t215 + (-t100 - t101) * t212 + (-t212 * t271 + t215 * t272) * qJD(5) + (0.2e1 * t196 ^ 2 * t309 + t133 + t134 - t180) * qJD(4)) * t213; m(4) * t226 + t128 * mrSges(4,1) + t119 + (-t14 + (t212 * t294 + t215 * t295 + t95) * qJD(4) + m(6) * (-t19 * t263 + t20 * t261 - t18) + m(7) * (t13 * t261 - t263 * t6 - t5) + t322) * t216 + (t66 + (t23 + t24) * t215 + (-t25 - t26) * t212 + (t56 + t296) * qJD(4) + (-t212 * t295 + t215 * t294) * qJD(5) + m(6) * (qJD(4) * t45 - t19 * t257 - t20 * t259 - t212 * t4 + t215 * t3) + m(7) * (qJD(4) * t27 - t1 * t212 - t13 * t259 + t2 * t215 - t257 * t6) + m(5) * (-qJD(4) * t53 + t21)) * t213; (-t108 + (t269 * t215 + t268 * t212 + m(7) * (t105 * t215 - t212 * t97) + m(6) * (-t110 * t212 + t111 * t215) - t196 * t309) * qJD(4) - t324) * t216 + ((t125 + t126) * t215 + (-t123 - t124) * t212 + (t150 + t151) * qJD(4) + (-t212 * t269 + t215 * t268) * qJD(5) + m(7) * (qJD(4) * t145 - t105 * t259 - t212 * t58 + t215 * t63 - t257 * t97) + (-t110 * t257 - t111 * t259 - t212 * t78 + t215 * t77 + t246) * m(6)) * t213; 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (-0.1e1 + t265) * t213 * t260; t254 + m(7) * (t1 * t172 + t13 * t143 + t144 * t6 - t175 * t2 + t198 * t5) + (-t4 * mrSges(6,3) - t1 * mrSges(7,3) + (-m(6) * t4 - t26) * pkin(10) + (-t13 * mrSges(7,3) + pkin(5) * t56 - pkin(10) * t60 + t20 * t325 + t27 * t311 + t252) * qJD(5) + t253) * t212 + (t2 * mrSges(7,3) + t3 * mrSges(6,3) + (-t19 * mrSges(6,3) - t6 * mrSges(7,3) + t251) * qJD(5) + (m(6) * (-qJD(5) * t19 + t3) + t24 - qJD(5) * t62) * pkin(10) - t255) * t215 + t198 * t14 + t172 * t25 + t5 * t173 + t18 * t174 - t175 * t23 + t27 * t155 + t45 * t156 + t143 * t59 + t144 * t61 + t232 * t86 + t233 * t93 + t234 * t92 + t235 * t112 + t230 * t44 + t231 * t43 + t22 * mrSges(5,1) - t21 * mrSges(5,2) + (-t310 - t15) * pkin(4); t205 + t213 * t196 * t156 + m(7) * (t105 * t143 + t115 * t198 + t144 * t97 + t172 * t58 - t175 * t63) + t198 * t107 + t143 * t166 + t144 * t168 + t172 * t123 + t115 * t173 - t175 * t125 + t145 * t155 - pkin(4) * t108 - t235 * t216 + ((-Ifges(5,6) + t232) * t213 + (t213 * mrSges(5,2) + (-m(6) * pkin(4) + t281) * t216) * t196) * qJD(4) + (-t78 * mrSges(6,3) - t58 * mrSges(7,3) - t234 * t213 - t231 * t260 + (-m(6) * t78 - t124) * pkin(10) + (-t105 * mrSges(7,3) + pkin(5) * t150 - pkin(10) * t167 + t111 * t325 + t145 * t311 - t230 * t213 - t237) * qJD(5) + t238) * t212 + (t77 * mrSges(6,3) + t63 * mrSges(7,3) + t233 * t213 + t230 * t260 + (m(6) * t77 + t126) * pkin(10) + (-t97 * mrSges(7,3) - t110 * mrSges(6,3) - t231 * t213 + (-m(6) * t110 - t169) * pkin(10) + t236) * qJD(5) + t239) * t215; m(7) * (-pkin(5) * t243 + t143 * t275 - t144 * t277 - t172 * t242 + t175 * t244) + ((t173 + t281) * t213 + m(7) * (-t172 * t276 - t175 * t274 + t198 * t213) + m(6) * (t265 * t305 - t306)) * qJD(4) + (-t156 - t155 + (-mrSges(5,2) + (mrSges(6,3) + mrSges(7,3)) * t265) * qJD(4)) * t216; 0.2e1 * t198 * t155 + (-t143 * t175 + t144 * t172) * t320 - 0.2e1 * pkin(4) * t156 + (t144 * t317 + t163 + t164 + (-t175 * t317 - t178 - t179 + 0.2e1 * (m(7) * t198 + t173) * pkin(5)) * qJD(5)) * t212 + (0.2e1 * t143 * mrSges(7,3) + t160 + t161 + (t172 * t317 + t181 + t182) * qJD(5)) * t215; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 + (m(7) * t1 + t25) * pkin(5) + t329; mrSges(6,1) * t78 + mrSges(7,1) * t58 - mrSges(6,2) * t77 - mrSges(7,2) * t63 - t330 * t240 + (m(7) * t58 + t123) * pkin(5) + (-t333 - t334) * t258 + t266 + t267; (t212 * t250 - t293) * t260 + (t250 * t215 + (mrSges(6,2) + mrSges(7,2)) * t212) * t258 + t248; -mrSges(7,2) * t143 + t203 + t204 + (mrSges(7,1) + t311) * t144 + ((-mrSges(6,1) * pkin(10) - mrSges(7,3) * pkin(5)) * t215 + (mrSges(6,2) * pkin(10) - t330) * t212) * qJD(5); 0; m(7) * t5 + t14; t324; m(7) * t262; t259 * t311 + t155; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
