% Calculate time derivative of joint inertia matrix for
% S6RRPRRP9
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:44
% EndTime: 2019-03-09 12:28:00
% DurationCPUTime: 7.73s
% Computational Cost: add. (9391->648), mult. (24625->910), div. (0->0), fcn. (24677->10), ass. (0->255)
t238 = cos(qJ(5));
t330 = Ifges(6,6) + Ifges(7,6);
t334 = t330 * t238;
t235 = sin(qJ(5));
t333 = (Ifges(6,5) + Ifges(7,5)) * t235;
t232 = sin(pkin(6));
t332 = 0.2e1 * t232;
t231 = sin(pkin(11));
t233 = cos(pkin(11));
t236 = sin(qJ(4));
t239 = cos(qJ(4));
t191 = t231 * t239 + t233 * t236;
t278 = qJD(5) * t238;
t265 = t191 * t278;
t190 = t231 * t236 - t239 * t233;
t185 = t190 * qJD(4);
t289 = t235 * t185;
t242 = t265 - t289;
t279 = qJD(5) * t235;
t266 = t191 * t279;
t288 = t238 * t185;
t241 = t266 + t288;
t234 = cos(pkin(6));
t237 = sin(qJ(2));
t292 = t232 * t237;
t183 = -t231 * t292 + t233 * t234;
t184 = t231 * t234 + t233 * t292;
t126 = t183 * t236 + t184 * t239;
t240 = cos(qJ(2));
t291 = t232 * t240;
t109 = -t235 * t126 - t238 * t291;
t282 = qJD(2) * t232;
t268 = t237 * t282;
t245 = t239 * t183 - t184 * t236;
t267 = t240 * t282;
t97 = qJD(4) * t245 - t190 * t267;
t54 = qJD(5) * t109 + t235 * t268 + t238 * t97;
t243 = -t238 * t126 + t235 * t291;
t55 = qJD(5) * t243 - t235 * t97 + t238 * t268;
t98 = qJD(4) * t126 + t191 * t267;
t7 = Ifges(7,5) * t54 + Ifges(7,6) * t55 + Ifges(7,3) * t98;
t8 = Ifges(6,5) * t54 + Ifges(6,6) * t55 + Ifges(6,3) * t98;
t329 = t7 + t8;
t325 = -m(6) * pkin(10) - mrSges(6,3);
t223 = -pkin(3) * t233 - pkin(2);
t139 = pkin(4) * t190 - pkin(10) * t191 + t223;
t307 = pkin(9) + qJ(3);
t203 = t307 * t231;
t204 = t307 * t233;
t153 = -t203 * t236 + t204 * t239;
t146 = t238 * t153;
t95 = t235 * t139 + t146;
t324 = -t239 * t203 - t204 * t236;
t323 = -m(6) * pkin(4) - mrSges(6,1) * t238 + mrSges(6,2) * t235;
t322 = 2 * m(4);
t321 = 2 * m(5);
t320 = 0.2e1 * m(6);
t319 = 2 * m(7);
t318 = -2 * mrSges(3,3);
t317 = -2 * mrSges(5,3);
t316 = -2 * mrSges(7,3);
t123 = qJD(3) * t191 + qJD(4) * t153;
t315 = 0.2e1 * t123;
t314 = -0.2e1 * t324;
t311 = m(7) * pkin(5);
t310 = t233 / 0.2e1;
t309 = pkin(1) * t240;
t306 = -qJ(6) - pkin(10);
t188 = t234 * t237 * pkin(1) + pkin(8) * t291;
t169 = qJ(3) * t234 + t188;
t170 = (-pkin(2) * t240 - qJ(3) * t237 - pkin(1)) * t232;
t121 = t233 * t169 + t231 * t170;
t102 = pkin(9) * t183 + t121;
t120 = -t231 * t169 + t233 * t170;
t86 = -pkin(3) * t291 - t184 * pkin(9) + t120;
t47 = t239 * t102 + t236 * t86;
t37 = -pkin(10) * t291 + t47;
t218 = pkin(8) * t292;
t172 = t218 + (-pkin(2) - t309) * t234;
t132 = -t183 * pkin(3) + t172;
t59 = -pkin(4) * t245 - t126 * pkin(10) + t132;
t22 = t235 * t59 + t238 * t37;
t305 = mrSges(6,2) * t238;
t304 = Ifges(4,4) * t231;
t303 = Ifges(4,4) * t233;
t302 = Ifges(6,4) * t235;
t301 = Ifges(6,4) * t238;
t300 = Ifges(7,4) * t235;
t299 = Ifges(7,4) * t238;
t175 = t188 * qJD(2);
t251 = t231 * t267;
t151 = pkin(3) * t251 + t175;
t298 = t151 * mrSges(5,1);
t297 = t151 * mrSges(5,2);
t296 = t123 * t324;
t295 = t191 * t235;
t294 = t191 * t238;
t290 = t233 * t240;
t246 = -Ifges(7,2) * t235 + t299;
t113 = Ifges(7,6) * t190 + t191 * t246;
t247 = -Ifges(6,2) * t235 + t301;
t114 = Ifges(6,6) * t190 + t191 * t247;
t287 = -t113 - t114;
t248 = Ifges(7,1) * t238 - t300;
t115 = Ifges(7,5) * t190 + t191 * t248;
t249 = Ifges(6,1) * t238 - t302;
t116 = Ifges(6,5) * t190 + t191 * t249;
t286 = t115 + t116;
t157 = (-qJD(3) * t237 + (pkin(2) * t237 - qJ(3) * t240) * qJD(2)) * t232;
t277 = t234 * t309;
t174 = -pkin(8) * t268 + qJD(2) * t277;
t162 = qJD(3) * t234 + t174;
t108 = t231 * t157 + t233 * t162;
t186 = t191 * qJD(4);
t285 = -Ifges(7,5) * t288 + Ifges(7,3) * t186;
t284 = -Ifges(6,5) * t288 + Ifges(6,3) * t186;
t283 = -Ifges(5,5) * t185 - Ifges(5,6) * t186;
t158 = t233 * mrSges(4,2) * t267 + mrSges(4,1) * t251;
t192 = mrSges(7,1) * t279 + mrSges(7,2) * t278;
t281 = qJD(4) * t236;
t280 = qJD(4) * t239;
t10 = Ifges(6,4) * t54 + Ifges(6,2) * t55 + Ifges(6,6) * t98;
t9 = Ifges(7,4) * t54 + Ifges(7,2) * t55 + Ifges(7,6) * t98;
t276 = -t9 / 0.2e1 - t10 / 0.2e1;
t275 = Ifges(5,5) * t97 - Ifges(5,6) * t98 + Ifges(5,3) * t268;
t11 = Ifges(7,1) * t54 + Ifges(7,4) * t55 + Ifges(7,5) * t98;
t12 = Ifges(6,1) * t54 + Ifges(6,4) * t55 + Ifges(6,5) * t98;
t274 = t11 / 0.2e1 + t12 / 0.2e1;
t40 = -Ifges(7,4) * t243 + Ifges(7,2) * t109 - Ifges(7,6) * t245;
t41 = -Ifges(6,4) * t243 + Ifges(6,2) * t109 - Ifges(6,6) * t245;
t273 = -t40 / 0.2e1 - t41 / 0.2e1;
t42 = -Ifges(7,1) * t243 + Ifges(7,4) * t109 - Ifges(7,5) * t245;
t43 = -Ifges(6,1) * t243 + Ifges(6,4) * t109 - Ifges(6,5) * t245;
t272 = t42 / 0.2e1 + t43 / 0.2e1;
t69 = -Ifges(7,4) * t241 - Ifges(7,2) * t242 + Ifges(7,6) * t186;
t70 = -Ifges(6,4) * t241 - Ifges(6,2) * t242 + Ifges(6,6) * t186;
t271 = t69 / 0.2e1 + t70 / 0.2e1;
t71 = -Ifges(7,1) * t241 - Ifges(7,4) * t242 + Ifges(7,5) * t186;
t72 = -Ifges(6,1) * t241 - Ifges(6,4) * t242 + Ifges(6,5) * t186;
t270 = t72 / 0.2e1 + t71 / 0.2e1;
t122 = -t190 * qJD(3) + qJD(4) * t324;
t138 = pkin(4) * t186 + pkin(10) * t185;
t269 = t238 * t122 + t235 * t138 + t139 * t278;
t264 = t113 / 0.2e1 + t114 / 0.2e1;
t263 = t115 / 0.2e1 + t116 / 0.2e1;
t227 = Ifges(7,5) * t278;
t228 = Ifges(6,5) * t278;
t262 = t227 / 0.2e1 + t228 / 0.2e1 - t330 * t279 / 0.2e1;
t196 = t246 * qJD(5);
t197 = t247 * qJD(5);
t261 = -t197 / 0.2e1 - t196 / 0.2e1;
t198 = t248 * qJD(5);
t199 = t249 * qJD(5);
t260 = t198 / 0.2e1 + t199 / 0.2e1;
t259 = t333 / 0.2e1 + t334 / 0.2e1;
t211 = Ifges(7,2) * t238 + t300;
t212 = Ifges(6,2) * t238 + t302;
t258 = t211 / 0.2e1 + t212 / 0.2e1;
t213 = Ifges(7,1) * t235 + t299;
t214 = Ifges(6,1) * t235 + t301;
t257 = t213 / 0.2e1 + t214 / 0.2e1;
t48 = t98 * mrSges(5,1) + t97 * mrSges(5,2);
t19 = -t55 * mrSges(7,1) + t54 * mrSges(7,2);
t21 = -t235 * t37 + t238 * t59;
t46 = -t236 * t102 + t239 * t86;
t256 = qJD(5) * t306;
t133 = t186 * mrSges(5,1) - t185 * mrSges(5,2);
t255 = -t122 * t235 + t238 * t138;
t94 = t238 * t139 - t153 * t235;
t107 = t233 * t157 - t231 * t162;
t253 = Ifges(5,5) * t268;
t252 = Ifges(5,6) * t268;
t36 = pkin(4) * t291 - t46;
t250 = mrSges(6,1) * t235 + t305;
t101 = -pkin(9) * t251 + t108;
t84 = (pkin(3) * t237 - pkin(9) * t290) * t282 + t107;
t18 = -t236 * t101 - t102 * t280 + t239 * t84 - t86 * t281;
t244 = qJ(6) * t185 - qJD(6) * t191;
t17 = t239 * t101 - t102 * t281 + t236 * t84 + t86 * t280;
t15 = pkin(10) * t268 + t17;
t32 = t98 * pkin(4) - t97 * pkin(10) + t151;
t3 = t238 * t15 + t235 * t32 + t59 * t278 - t279 * t37;
t80 = mrSges(7,1) * t242 - mrSges(7,2) * t241;
t16 = -pkin(4) * t268 - t18;
t4 = -t22 * qJD(5) - t15 * t235 + t238 * t32;
t224 = -pkin(5) * t238 - pkin(4);
t216 = Ifges(3,5) * t267;
t208 = t306 * t238;
t206 = -mrSges(7,1) * t238 + mrSges(7,2) * t235;
t205 = t306 * t235;
t193 = t250 * qJD(5);
t187 = -t218 + t277;
t182 = -qJD(6) * t235 + t238 * t256;
t181 = qJD(6) * t238 + t235 * t256;
t164 = (mrSges(4,1) * t237 - mrSges(4,3) * t290) * t282;
t163 = (-mrSges(4,3) * t231 * t240 - mrSges(4,2) * t237) * t282;
t156 = -mrSges(4,1) * t291 - t184 * mrSges(4,3);
t155 = mrSges(4,2) * t291 + t183 * mrSges(4,3);
t148 = Ifges(5,1) * t191 - Ifges(5,4) * t190;
t147 = Ifges(5,4) * t191 - Ifges(5,2) * t190;
t145 = (t237 * Ifges(4,5) + (t233 * Ifges(4,1) - t304) * t240) * t282;
t144 = (t237 * Ifges(4,6) + (-t231 * Ifges(4,2) + t303) * t240) * t282;
t143 = mrSges(6,1) * t190 - mrSges(6,3) * t294;
t142 = mrSges(7,1) * t190 - mrSges(7,3) * t294;
t141 = -mrSges(6,2) * t190 - mrSges(6,3) * t295;
t140 = -mrSges(7,2) * t190 - mrSges(7,3) * t295;
t137 = t250 * t191;
t136 = (mrSges(7,1) * t235 + mrSges(7,2) * t238) * t191;
t135 = -Ifges(5,1) * t185 - Ifges(5,4) * t186;
t134 = -Ifges(5,4) * t185 - Ifges(5,2) * t186;
t124 = pkin(5) * t295 - t324;
t118 = -mrSges(5,1) * t291 - t126 * mrSges(5,3);
t117 = mrSges(5,2) * t291 + mrSges(5,3) * t245;
t112 = Ifges(6,3) * t190 + (Ifges(6,5) * t238 - Ifges(6,6) * t235) * t191;
t111 = Ifges(7,3) * t190 + (Ifges(7,5) * t238 - Ifges(7,6) * t235) * t191;
t106 = -mrSges(6,2) * t186 - mrSges(6,3) * t242;
t105 = -mrSges(7,2) * t186 - mrSges(7,3) * t242;
t104 = mrSges(6,1) * t186 + mrSges(6,3) * t241;
t103 = mrSges(7,1) * t186 + mrSges(7,3) * t241;
t81 = mrSges(6,1) * t242 - mrSges(6,2) * t241;
t78 = -mrSges(5,2) * t268 - mrSges(5,3) * t98;
t77 = mrSges(5,1) * t268 - mrSges(5,3) * t97;
t76 = pkin(5) * t242 + t123;
t75 = -qJ(6) * t295 + t95;
t74 = Ifges(5,1) * t126 + Ifges(5,4) * t245 - Ifges(5,5) * t291;
t73 = Ifges(5,4) * t126 + Ifges(5,2) * t245 - Ifges(5,6) * t291;
t68 = -Ifges(6,5) * t266 - Ifges(6,6) * t242 + t284;
t67 = -Ifges(7,5) * t266 - Ifges(7,6) * t242 + t285;
t66 = -mrSges(6,1) * t245 + mrSges(6,3) * t243;
t65 = -mrSges(7,1) * t245 + mrSges(7,3) * t243;
t64 = mrSges(6,2) * t245 + mrSges(6,3) * t109;
t63 = mrSges(7,2) * t245 + mrSges(7,3) * t109;
t62 = pkin(5) * t190 - qJ(6) * t294 + t94;
t61 = -mrSges(6,1) * t109 - mrSges(6,2) * t243;
t60 = -mrSges(7,1) * t109 - mrSges(7,2) * t243;
t45 = Ifges(5,1) * t97 - Ifges(5,4) * t98 + t253;
t44 = Ifges(5,4) * t97 - Ifges(5,2) * t98 + t252;
t39 = -Ifges(6,5) * t243 + Ifges(6,6) * t109 - Ifges(6,3) * t245;
t38 = -Ifges(7,5) * t243 + Ifges(7,6) * t109 - Ifges(7,3) * t245;
t34 = -qJD(5) * t95 + t255;
t33 = -t153 * t279 + t269;
t29 = -pkin(5) * t109 + t36;
t28 = -mrSges(6,2) * t98 + mrSges(6,3) * t55;
t27 = -mrSges(7,2) * t98 + mrSges(7,3) * t55;
t26 = mrSges(6,1) * t98 - mrSges(6,3) * t54;
t25 = mrSges(7,1) * t98 - mrSges(7,3) * t54;
t24 = -qJ(6) * t265 + (-qJD(5) * t153 + t244) * t235 + t269;
t23 = pkin(5) * t186 + t244 * t238 + (-t146 + (qJ(6) * t191 - t139) * t235) * qJD(5) + t255;
t20 = -mrSges(6,1) * t55 + mrSges(6,2) * t54;
t13 = qJ(6) * t109 + t22;
t6 = -pkin(5) * t245 + qJ(6) * t243 + t21;
t5 = -pkin(5) * t55 + t16;
t2 = qJ(6) * t55 + qJD(6) * t109 + t3;
t1 = pkin(5) * t98 - qJ(6) * t54 + qJD(6) * t243 + t4;
t14 = [t234 * t216 + (t9 + t10) * t109 + ((t188 * t318 + Ifges(4,5) * t184 + Ifges(5,5) * t126 - 0.2e1 * Ifges(3,6) * t234 + Ifges(4,6) * t183 + Ifges(5,6) * t245 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t237) * t332) * t237 + (t233 * (Ifges(4,1) * t184 + Ifges(4,4) * t183) - t231 * (Ifges(4,4) * t184 + Ifges(4,2) * t183) + t187 * t318 + Ifges(3,5) * t234 + (-pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t233 + Ifges(4,6) * t231 + Ifges(3,4)) * t240) * t332 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3)) * t292) * t240) * t282 + t183 * t144 + 0.2e1 * t175 * (-mrSges(4,1) * t183 + mrSges(4,2) * t184) + t184 * t145 + 0.2e1 * t121 * t163 + 0.2e1 * t120 * t164 + 0.2e1 * t172 * t158 + 0.2e1 * t108 * t155 + 0.2e1 * t107 * t156 + 0.2e1 * t132 * t48 + 0.2e1 * t17 * t117 + 0.2e1 * t18 * t118 + t97 * t74 + 0.2e1 * t46 * t77 + 0.2e1 * t47 * t78 + 0.2e1 * t2 * t63 + 0.2e1 * t3 * t64 + 0.2e1 * t1 * t65 + 0.2e1 * t4 * t66 + 0.2e1 * t5 * t60 + 0.2e1 * t16 * t61 + 0.2e1 * t36 * t20 + 0.2e1 * t13 * t27 + 0.2e1 * t22 * t28 + 0.2e1 * t29 * t19 + 0.2e1 * t6 * t25 + 0.2e1 * t21 * t26 + 0.2e1 * m(3) * (t174 * t188 - t175 * t187) + (t38 + t39 - t73) * t98 + (t40 + t41) * t55 - (t11 + t12) * t243 - (-t44 + 0.2e1 * t298 + t329) * t245 + (t42 + t43) * t54 + (t1 * t6 + t13 * t2 + t29 * t5) * t319 + (t16 * t36 + t21 * t4 + t22 * t3) * t320 + (t132 * t151 + t17 * t47 + t18 * t46) * t321 + (t107 * t120 + t108 * t121 + t172 * t175) * t322 - t275 * t291 + 0.2e1 * t174 * (-t234 * mrSges(3,2) + mrSges(3,3) * t291) - 0.2e1 * t175 * (mrSges(3,1) * t234 - mrSges(3,3) * t292) + (t45 + 0.2e1 * t297) * t126; (qJD(3) * t155 + qJ(3) * t163 + t108 * mrSges(4,3) - t175 * mrSges(4,1) + t144 / 0.2e1) * t233 + t216 + (t61 - t118) * t123 + (-qJD(3) * t156 - qJ(3) * t164 - t107 * mrSges(4,3) + t145 / 0.2e1 + t175 * mrSges(4,2)) * t231 + m(7) * (t1 * t62 + t124 * t5 + t13 * t24 + t2 * t75 + t23 * t6 + t29 * t76) + (-t147 / 0.2e1 + t111 / 0.2e1 + t112 / 0.2e1) * t98 + t223 * t48 - t174 * mrSges(3,2) - t175 * mrSges(3,1) + t153 * t78 - pkin(2) * t158 + t2 * t140 + t3 * t141 + t1 * t142 + t4 * t143 + t97 * t148 / 0.2e1 + t132 * t133 + t126 * t135 / 0.2e1 + t5 * t136 + t16 * t137 + t124 * t19 + t122 * t117 + t6 * t103 + t21 * t104 + m(4) * (-pkin(2) * t175 + (-t120 * t231 + t121 * t233) * qJD(3) + (-t107 * t231 + t108 * t233) * qJ(3)) + t13 * t105 + t22 * t106 + t95 * t28 + t94 * t26 + t76 * t60 + t29 * t80 + t36 * t81 + t24 * t63 + t33 * t64 + t23 * t65 + t34 * t66 + t75 * t27 + t62 * t25 - t270 * t243 - (-t134 / 0.2e1 + t67 / 0.2e1 + t68 / 0.2e1) * t245 + t263 * t54 + t264 * t55 + t271 * t109 - (-t46 * mrSges(5,3) + t74 / 0.2e1 + t272 * t238 + t273 * t235) * t185 - (t20 - t77) * t324 + m(5) * (t122 * t47 - t123 * t46 + t151 * t223 + t153 * t17 + t18 * t324) + m(6) * (t123 * t36 - t16 * t324 + t21 * t34 + t22 * t33 + t3 * t95 + t4 * t94) + (-t18 * mrSges(5,3) + t253 / 0.2e1 + t297 + t45 / 0.2e1 + t274 * t238 + t276 * t235 + (-t235 * t272 + t238 * t273) * qJD(5)) * t191 + (-t17 * mrSges(5,3) - t252 / 0.2e1 - t44 / 0.2e1 + t298 + t7 / 0.2e1 + t8 / 0.2e1) * t190 + (-t240 * t283 / 0.2e1 + ((-Ifges(3,6) + Ifges(4,5) * t231 / 0.2e1 + Ifges(4,6) * t310) * t237 + (-t231 * (Ifges(4,2) * t233 + t304) / 0.2e1 + (Ifges(4,1) * t231 + t303) * t310) * t240) * qJD(2)) * t232 + (-t47 * mrSges(5,3) - t73 / 0.2e1 + t38 / 0.2e1 + t39 / 0.2e1) * t186; 0.2e1 * t62 * t103 + 0.2e1 * t94 * t104 + 0.2e1 * t75 * t105 + 0.2e1 * t95 * t106 + t137 * t315 + 0.2e1 * t124 * t80 + 0.2e1 * t223 * t133 + 0.2e1 * t76 * t136 + 0.2e1 * t24 * t140 + 0.2e1 * t33 * t141 + 0.2e1 * t23 * t142 + 0.2e1 * t34 * t143 + t81 * t314 + (t122 * t317 - t134 + t67 + t68) * t190 + (t153 * t317 + t111 + t112 - t147) * t186 + (t122 * t153 - t296) * t321 + (t33 * t95 + t34 * t94 - t296) * t320 + (t124 * t76 + t23 * t62 + t24 * t75) * t319 - (mrSges(5,3) * t314 + t235 * t287 + t238 * t286 + t148) * t185 + (mrSges(5,3) * t315 + t135 + (t71 + t72) * t238 + (-t69 - t70) * t235 + (-t235 * t286 + t238 * t287) * qJD(5)) * t191 + (qJ(3) * t322 + 0.2e1 * mrSges(4,3)) * (t231 ^ 2 + t233 ^ 2) * qJD(3); (t25 + t26) * t238 + (t27 + t28) * t235 + ((t63 + t64) * t238 + (-t65 - t66) * t235) * qJD(5) + m(7) * (t1 * t238 + t2 * t235 + (t13 * t238 - t235 * t6) * qJD(5)) + m(6) * (t235 * t3 + t238 * t4 + (-t21 * t235 + t22 * t238) * qJD(5)) + m(5) * t151 + m(4) * t175 + t48 + t158; (t103 + t104) * t238 + (t105 + t106) * t235 + ((t140 + t141) * t238 + (-t142 - t143) * t235) * qJD(5) + m(6) * (t235 * t33 + t238 * t34 + (-t235 * t94 + t238 * t95) * qJD(5)) + m(7) * (t23 * t238 + t235 * t24 + (-t235 * t62 + t238 * t75) * qJD(5)) + t133; 0; m(7) * (t1 * t205 + t13 * t181 + t182 * t6 - t2 * t208 + t224 * t5) + t224 * t19 + t5 * t206 - t208 * t27 + t29 * t192 + t36 * t193 + t205 * t25 + t182 * t65 + t181 * t63 - pkin(4) * t20 - t17 * mrSges(5,2) + t18 * mrSges(5,1) + t257 * t54 + t258 * t55 + t259 * t98 - t260 * t243 - t261 * t109 - t262 * t245 + (t3 * mrSges(6,3) + t2 * mrSges(7,3) + (-t21 * mrSges(6,3) - t6 * mrSges(7,3) + t272) * qJD(5) + (m(6) * (-qJD(5) * t21 + t3) + t28 - qJD(5) * t66) * pkin(10) - t276) * t238 + (-t4 * mrSges(6,3) - t1 * mrSges(7,3) + (-m(6) * t4 - t26) * pkin(10) + (-t13 * mrSges(7,3) + pkin(5) * t60 - pkin(10) * t64 + t22 * t325 + t29 * t311 + t273) * qJD(5) + t274) * t235 + t275 + t323 * t16; m(7) * (t181 * t75 + t182 * t62 + t205 * t23 - t208 * t24 + t224 * t76) + t224 * t80 + t76 * t206 - t208 * t105 + t124 * t192 - t324 * t193 + t205 * t103 + t182 * t142 + t181 * t140 - t122 * mrSges(5,2) - pkin(4) * t81 + t262 * t190 + t259 * t186 + (-mrSges(5,1) + t323) * t123 + (-t34 * mrSges(6,3) - t23 * mrSges(7,3) + t261 * t191 + t258 * t185 + (-m(6) * t34 - t104) * pkin(10) + (-t75 * mrSges(7,3) + pkin(5) * t136 - pkin(10) * t141 + t124 * t311 - t191 * t257 + t325 * t95 - t264) * qJD(5) + t270) * t235 + (t33 * mrSges(6,3) + t24 * mrSges(7,3) + t260 * t191 - t257 * t185 + (m(6) * t33 + t106) * pkin(10) + (-t62 * mrSges(7,3) - t94 * mrSges(6,3) - t258 * t191 + (-m(6) * t94 - t143) * pkin(10) + t263) * qJD(5) + t271) * t238 + t283; m(7) * (t181 * t235 + t182 * t238 + (-t205 * t235 - t208 * t238) * qJD(5)); 0.2e1 * t224 * t192 - 0.2e1 * pkin(4) * t193 + (-t181 * t208 + t182 * t205) * t319 + (t182 * t316 + t198 + t199 + (-t208 * t316 - t211 - t212 + 0.2e1 * (m(7) * t224 + t206) * pkin(5)) * qJD(5)) * t235 + (0.2e1 * t181 * mrSges(7,3) + t196 + t197 + (t205 * t316 + t213 + t214) * qJD(5)) * t238; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 + (m(7) * t1 + t25) * pkin(5) + t329; mrSges(6,1) * t34 + mrSges(7,1) * t23 - mrSges(6,2) * t33 - mrSges(7,2) * t24 + t330 * t289 + (m(7) * t23 + t103) * pkin(5) + (-t333 - t334) * t191 * qJD(5) + t284 + t285; (-t305 + (-mrSges(6,1) - t311) * t235) * qJD(5) - t192; -mrSges(7,2) * t181 + t227 + t228 + (mrSges(7,1) + t311) * t182 + ((-mrSges(6,1) * pkin(10) - (mrSges(7,3) * pkin(5))) * t238 + (mrSges(6,2) * pkin(10) - t330) * t235) * qJD(5); 0; m(7) * t5 + t19; m(7) * t76 + t80; 0; t279 * t311 + t192; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
