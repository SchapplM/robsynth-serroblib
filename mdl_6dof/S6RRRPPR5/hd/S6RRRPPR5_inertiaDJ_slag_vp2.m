% Calculate time derivative of joint inertia matrix for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:34
% EndTime: 2019-03-09 15:38:49
% DurationCPUTime: 6.48s
% Computational Cost: add. (12147->652), mult. (32042->978), div. (0->0), fcn. (31774->12), ass. (0->274)
t342 = Ifges(4,3) + Ifges(5,3);
t254 = sin(pkin(12));
t257 = cos(pkin(12));
t259 = sin(qJ(6));
t262 = cos(qJ(6));
t268 = t254 * t259 - t257 * t262;
t323 = -t268 / 0.2e1;
t227 = t254 * t262 + t257 * t259;
t322 = t227 / 0.2e1;
t318 = t254 / 0.2e1;
t317 = t257 / 0.2e1;
t258 = cos(pkin(6));
t261 = sin(qJ(2));
t256 = sin(pkin(6));
t264 = cos(qJ(2));
t288 = t256 * t264;
t223 = t258 * t261 * pkin(1) + pkin(8) * t288;
t201 = pkin(9) * t258 + t223;
t202 = (-pkin(2) * t264 - pkin(9) * t261 - pkin(1)) * t256;
t260 = sin(qJ(3));
t263 = cos(qJ(3));
t141 = t263 * t201 + t260 * t202;
t289 = t256 * t261;
t243 = pkin(8) * t289;
t316 = pkin(1) * t264;
t222 = t258 * t316 - t243;
t255 = sin(pkin(11));
t298 = cos(pkin(11));
t272 = t298 * t260;
t226 = t255 * t263 + t272;
t214 = t226 * qJD(3);
t290 = t255 * t260;
t266 = t263 * t298 - t290;
t215 = t266 * qJD(3);
t284 = qJD(3) * t260;
t278 = pkin(3) * t284;
t127 = pkin(4) * t214 - qJ(5) * t215 - qJD(5) * t226 + t278;
t313 = -qJ(4) - pkin(9);
t273 = qJD(3) * t313;
t213 = qJD(4) * t263 + t260 * t273;
t265 = -qJD(4) * t260 + t263 * t273;
t148 = t213 * t298 + t255 * t265;
t80 = t257 * t127 - t148 * t254;
t81 = t254 * t127 + t257 * t148;
t341 = -t254 * t80 + t257 * t81;
t218 = t258 * t263 - t260 * t289;
t286 = qJD(2) * t256;
t274 = t264 * t286;
t185 = qJD(3) * t218 + t263 * t274;
t219 = t258 * t260 + t263 * t289;
t285 = qJD(2) * t261;
t275 = t256 * t285;
t204 = (pkin(2) * t261 - pkin(9) * t264) * t286;
t205 = t222 * qJD(2);
t89 = -qJD(3) * t141 + t263 * t204 - t205 * t260;
t61 = pkin(3) * t275 - qJ(4) * t185 - qJD(4) * t219 + t89;
t184 = -qJD(3) * t219 - t260 * t274;
t283 = qJD(3) * t263;
t88 = -t201 * t284 + t202 * t283 + t260 * t204 + t263 * t205;
t69 = qJ(4) * t184 + qJD(4) * t218 + t88;
t24 = t255 * t61 + t298 * t69;
t21 = (qJ(5) * t285 - qJD(5) * t264) * t256 + t24;
t117 = -t184 * t298 + t185 * t255;
t118 = t255 * t184 + t185 * t298;
t206 = t223 * qJD(2);
t139 = -t184 * pkin(3) + t206;
t154 = t255 * t218 + t219 * t298;
t41 = t117 * pkin(4) - t118 * qJ(5) - t154 * qJD(5) + t139;
t11 = -t21 * t254 + t257 * t41;
t12 = t257 * t21 + t254 * t41;
t340 = -t11 * t254 + t12 * t257;
t339 = Ifges(6,5) * t318 + Ifges(7,5) * t322 + Ifges(6,6) * t317 + Ifges(7,6) * t323;
t216 = t268 * qJD(6);
t338 = 2 * m(5);
t337 = 2 * m(6);
t336 = 2 * m(7);
t335 = -2 * mrSges(3,3);
t334 = 0.2e1 * t206;
t333 = m(5) * pkin(3);
t101 = -t118 * t254 + t257 * t275;
t332 = t101 / 0.2e1;
t217 = t227 * qJD(6);
t165 = -Ifges(7,5) * t216 - Ifges(7,6) * t217;
t331 = t165 / 0.2e1;
t166 = -Ifges(7,4) * t216 - Ifges(7,2) * t217;
t330 = t166 / 0.2e1;
t168 = -Ifges(7,1) * t216 - Ifges(7,4) * t217;
t329 = t168 / 0.2e1;
t177 = Ifges(7,4) * t227 - Ifges(7,2) * t268;
t327 = t177 / 0.2e1;
t179 = Ifges(7,1) * t227 - Ifges(7,4) * t268;
t326 = t179 / 0.2e1;
t325 = -t216 / 0.2e1;
t324 = -t217 / 0.2e1;
t308 = Ifges(6,4) * t257;
t320 = Ifges(6,1) * t318 + t308 / 0.2e1;
t319 = -t254 / 0.2e1;
t315 = pkin(3) * t255;
t247 = qJ(5) + t315;
t314 = pkin(10) + t247;
t140 = -t260 * t201 + t263 * t202;
t106 = -pkin(3) * t288 - t219 * qJ(4) + t140;
t119 = qJ(4) * t218 + t141;
t68 = t255 * t106 + t298 * t119;
t57 = -qJ(5) * t288 + t68;
t153 = -t218 * t298 + t219 * t255;
t200 = t243 + (-pkin(2) - t316) * t258;
t161 = -t218 * pkin(3) + t200;
t74 = t153 * pkin(4) - t154 * qJ(5) + t161;
t31 = t254 * t74 + t257 * t57;
t312 = mrSges(7,3) * t227;
t311 = Ifges(4,4) * t260;
t310 = Ifges(4,4) * t263;
t309 = Ifges(6,4) * t254;
t305 = t205 * mrSges(3,2);
t304 = t206 * mrSges(3,1);
t303 = t216 * mrSges(7,3);
t302 = t217 * mrSges(7,3);
t301 = t268 * mrSges(7,3);
t147 = t213 * t255 - t298 * t265;
t237 = t313 * t263;
t186 = -t237 * t255 - t313 * t272;
t297 = t147 * t186;
t296 = t215 * t254;
t295 = t215 * t257;
t294 = t226 * t254;
t293 = t226 * t257;
t292 = t247 * t254;
t291 = t247 * t257;
t250 = -pkin(3) * t263 - pkin(2);
t170 = -pkin(4) * t266 - qJ(5) * t226 + t250;
t187 = -t237 * t298 + t290 * t313;
t110 = t254 * t170 + t257 * t187;
t146 = mrSges(6,1) * t296 + mrSges(6,2) * t295;
t282 = qJD(5) * t254;
t281 = qJD(5) * t257;
t280 = 0.2e1 * t256;
t102 = t118 * t257 + t254 * t275;
t132 = -t254 * t154 - t257 * t288;
t133 = t257 * t154 - t254 * t288;
t75 = t132 * t262 - t133 * t259;
t28 = qJD(6) * t75 + t101 * t259 + t102 * t262;
t76 = t132 * t259 + t133 * t262;
t29 = -qJD(6) * t76 + t101 * t262 - t102 * t259;
t6 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t117;
t96 = -t215 * t268 - t217 * t226;
t97 = -t215 * t227 + t216 * t226;
t45 = Ifges(7,5) * t96 + Ifges(7,6) * t97 + Ifges(7,3) * t214;
t277 = Ifges(4,6) * t288;
t276 = t298 * pkin(3);
t10 = -t29 * mrSges(7,1) + t28 * mrSges(7,2);
t48 = -t97 * mrSges(7,1) + t96 * mrSges(7,2);
t30 = -t254 * t57 + t257 * t74;
t70 = t117 * mrSges(5,1) + t118 * mrSges(5,2);
t163 = t214 * mrSges(5,1) + t215 * mrSges(5,2);
t49 = -t101 * mrSges(6,1) + t102 * mrSges(6,2);
t109 = t257 * t170 - t187 * t254;
t249 = -t276 - pkin(4);
t23 = -t255 * t69 + t298 * t61;
t271 = Ifges(6,1) * t257 - t309;
t270 = -Ifges(6,2) * t254 + t308;
t269 = Ifges(6,5) * t257 - Ifges(6,6) * t254;
t16 = pkin(5) * t153 - pkin(10) * t133 + t30;
t19 = pkin(10) * t132 + t31;
t3 = t16 * t262 - t19 * t259;
t4 = t16 * t259 + t19 * t262;
t82 = -pkin(5) * t266 - pkin(10) * t293 + t109;
t95 = -pkin(10) * t294 + t110;
t43 = -t259 * t95 + t262 * t82;
t44 = t259 * t82 + t262 * t95;
t67 = t106 * t298 - t255 * t119;
t220 = t314 * t254;
t221 = t314 * t257;
t159 = -t220 * t262 - t221 * t259;
t160 = -t220 * t259 + t221 * t262;
t60 = pkin(4) * t288 - t67;
t267 = Ifges(4,5) * t185 + Ifges(5,5) * t118 + Ifges(4,6) * t184 - Ifges(5,6) * t117 + t275 * t342;
t22 = -pkin(4) * t275 - t23;
t251 = Ifges(4,5) * t283;
t242 = Ifges(3,5) * t274;
t239 = Ifges(4,1) * t260 + t310;
t238 = Ifges(4,2) * t263 + t311;
t235 = Ifges(6,2) * t257 + t309;
t233 = -mrSges(6,1) * t257 + mrSges(6,2) * t254;
t232 = -t257 * pkin(5) + t249;
t231 = (Ifges(4,1) * t263 - t311) * qJD(3);
t230 = (-Ifges(4,2) * t260 + t310) * qJD(3);
t229 = (mrSges(4,1) * t260 + mrSges(4,2) * t263) * qJD(3);
t212 = Ifges(5,5) * t215;
t210 = Ifges(5,6) * t214;
t189 = -mrSges(4,1) * t288 - t219 * mrSges(4,3);
t188 = mrSges(4,2) * t288 + t218 * mrSges(4,3);
t180 = Ifges(5,1) * t226 + Ifges(5,4) * t266;
t178 = Ifges(5,4) * t226 + Ifges(5,2) * t266;
t175 = mrSges(7,1) * t268 + mrSges(7,2) * t227;
t174 = -mrSges(5,1) * t266 + mrSges(5,2) * t226;
t172 = -mrSges(6,1) * t266 - mrSges(6,3) * t293;
t171 = mrSges(6,2) * t266 - mrSges(6,3) * t294;
t169 = Ifges(5,1) * t215 - Ifges(5,4) * t214;
t167 = Ifges(5,4) * t215 - Ifges(5,2) * t214;
t164 = mrSges(7,1) * t217 - mrSges(7,2) * t216;
t162 = (mrSges(6,1) * t254 + mrSges(6,2) * t257) * t226;
t158 = mrSges(6,1) * t214 - mrSges(6,3) * t295;
t157 = -mrSges(6,2) * t214 - mrSges(6,3) * t296;
t152 = t268 * t226;
t151 = t227 * t226;
t150 = mrSges(4,1) * t275 - mrSges(4,3) * t185;
t149 = -mrSges(4,2) * t275 + mrSges(4,3) * t184;
t145 = Ifges(4,1) * t219 + Ifges(4,4) * t218 - Ifges(4,5) * t288;
t144 = Ifges(4,4) * t219 + Ifges(4,2) * t218 - t277;
t143 = pkin(5) * t294 + t186;
t138 = -mrSges(5,1) * t288 - t154 * mrSges(5,3);
t137 = mrSges(5,2) * t288 - t153 * mrSges(5,3);
t136 = -Ifges(6,5) * t266 + t226 * t271;
t135 = -Ifges(6,6) * t266 + t226 * t270;
t134 = -Ifges(6,3) * t266 + t226 * t269;
t131 = -qJD(5) * t227 - qJD(6) * t160;
t130 = -qJD(5) * t268 + qJD(6) * t159;
t129 = -mrSges(7,1) * t266 + mrSges(7,3) * t152;
t128 = mrSges(7,2) * t266 - mrSges(7,3) * t151;
t126 = Ifges(6,5) * t214 + t215 * t271;
t125 = Ifges(6,6) * t214 + t215 * t270;
t124 = Ifges(6,3) * t214 + t215 * t269;
t121 = -mrSges(4,1) * t184 + mrSges(4,2) * t185;
t120 = pkin(5) * t296 + t147;
t108 = Ifges(4,1) * t185 + Ifges(4,4) * t184 + Ifges(4,5) * t275;
t107 = Ifges(4,4) * t185 + Ifges(4,2) * t184 + Ifges(4,6) * t275;
t104 = mrSges(5,1) * t275 - mrSges(5,3) * t118;
t103 = -mrSges(5,2) * t275 - mrSges(5,3) * t117;
t99 = mrSges(5,1) * t153 + mrSges(5,2) * t154;
t98 = mrSges(7,1) * t151 - mrSges(7,2) * t152;
t91 = Ifges(5,1) * t154 - Ifges(5,4) * t153 - Ifges(5,5) * t288;
t90 = Ifges(5,4) * t154 - Ifges(5,2) * t153 - Ifges(5,6) * t288;
t87 = mrSges(6,1) * t153 - mrSges(6,3) * t133;
t86 = -mrSges(6,2) * t153 + mrSges(6,3) * t132;
t85 = -Ifges(7,1) * t152 - Ifges(7,4) * t151 - Ifges(7,5) * t266;
t84 = -Ifges(7,4) * t152 - Ifges(7,2) * t151 - Ifges(7,6) * t266;
t83 = -Ifges(7,5) * t152 - Ifges(7,6) * t151 - Ifges(7,3) * t266;
t79 = -mrSges(7,2) * t214 + mrSges(7,3) * t97;
t78 = mrSges(7,1) * t214 - mrSges(7,3) * t96;
t77 = -mrSges(6,1) * t132 + mrSges(6,2) * t133;
t71 = -pkin(10) * t296 + t81;
t64 = Ifges(5,1) * t118 - Ifges(5,4) * t117 + Ifges(5,5) * t275;
t63 = Ifges(5,4) * t118 - Ifges(5,2) * t117 + Ifges(5,6) * t275;
t62 = pkin(5) * t214 - pkin(10) * t295 + t80;
t59 = mrSges(6,1) * t117 - mrSges(6,3) * t102;
t58 = -mrSges(6,2) * t117 + mrSges(6,3) * t101;
t55 = Ifges(6,1) * t133 + Ifges(6,4) * t132 + Ifges(6,5) * t153;
t54 = Ifges(6,4) * t133 + Ifges(6,2) * t132 + Ifges(6,6) * t153;
t53 = Ifges(6,5) * t133 + Ifges(6,6) * t132 + Ifges(6,3) * t153;
t51 = mrSges(7,1) * t153 - mrSges(7,3) * t76;
t50 = -mrSges(7,2) * t153 + mrSges(7,3) * t75;
t47 = Ifges(7,1) * t96 + Ifges(7,4) * t97 + Ifges(7,5) * t214;
t46 = Ifges(7,4) * t96 + Ifges(7,2) * t97 + Ifges(7,6) * t214;
t42 = -t132 * pkin(5) + t60;
t38 = -mrSges(7,1) * t75 + mrSges(7,2) * t76;
t37 = Ifges(6,1) * t102 + Ifges(6,4) * t101 + Ifges(6,5) * t117;
t36 = Ifges(6,4) * t102 + Ifges(6,2) * t101 + Ifges(6,6) * t117;
t35 = Ifges(6,5) * t102 + Ifges(6,6) * t101 + Ifges(6,3) * t117;
t34 = Ifges(7,1) * t76 + Ifges(7,4) * t75 + Ifges(7,5) * t153;
t33 = Ifges(7,4) * t76 + Ifges(7,2) * t75 + Ifges(7,6) * t153;
t32 = Ifges(7,5) * t76 + Ifges(7,6) * t75 + Ifges(7,3) * t153;
t18 = -mrSges(7,2) * t117 + mrSges(7,3) * t29;
t17 = mrSges(7,1) * t117 - mrSges(7,3) * t28;
t15 = -t101 * pkin(5) + t22;
t14 = -qJD(6) * t44 - t259 * t71 + t262 * t62;
t13 = qJD(6) * t43 + t259 * t62 + t262 * t71;
t9 = pkin(10) * t101 + t12;
t8 = Ifges(7,1) * t28 + Ifges(7,4) * t29 + Ifges(7,5) * t117;
t7 = Ifges(7,4) * t28 + Ifges(7,2) * t29 + Ifges(7,6) * t117;
t5 = pkin(5) * t117 - pkin(10) * t102 + t11;
t2 = -qJD(6) * t4 - t259 * t9 + t262 * t5;
t1 = qJD(6) * t3 + t259 * t5 + t262 * t9;
t20 = [(t139 * t161 + t23 * t67 + t24 * t68) * t338 + (-mrSges(4,1) * t218 + mrSges(4,2) * t219) * t334 + (t1 * t4 + t15 * t42 + t2 * t3) * t336 + (t11 * t30 + t12 * t31 + t22 * t60) * t337 + 0.2e1 * m(3) * (t205 * t223 - t206 * t222) + 0.2e1 * m(4) * (t140 * t89 + t141 * t88 + t200 * t206) + t218 * t107 + t219 * t108 + 0.2e1 * t88 * t188 + 0.2e1 * t89 * t189 + 0.2e1 * t200 * t121 + t184 * t144 + t185 * t145 + t154 * t64 + 0.2e1 * t161 * t70 + 0.2e1 * t141 * t149 + 0.2e1 * t140 * t150 + t132 * t36 + t133 * t37 + 0.2e1 * t24 * t137 + 0.2e1 * t23 * t138 + 0.2e1 * t139 * t99 + t118 * t91 + t101 * t54 + t102 * t55 + 0.2e1 * t68 * t103 + 0.2e1 * t67 * t104 + 0.2e1 * t12 * t86 + 0.2e1 * t11 * t87 + t76 * t8 + 0.2e1 * t22 * t77 + t75 * t7 + 0.2e1 * t2 * t51 + 0.2e1 * t31 * t58 + 0.2e1 * t30 * t59 + 0.2e1 * t60 * t49 + 0.2e1 * t1 * t50 + 0.2e1 * t42 * t10 + t28 * t34 + 0.2e1 * t15 * t38 + t29 * t33 + 0.2e1 * t3 * t17 + 0.2e1 * t4 * t18 + (mrSges(3,3) * t261 * t334 + (0.2e1 * t205 * mrSges(3,3) - t267) * t264 + ((t222 * t335 + Ifges(3,5) * t258 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t264) * t280) * t264 + (t223 * t335 + Ifges(4,5) * t219 + Ifges(5,5) * t154 - 0.2e1 * Ifges(3,6) * t258 + Ifges(4,6) * t218 - Ifges(5,6) * t153 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t261) * t280 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t342) * t288) * t261) * qJD(2)) * t256 + (-t63 + t35 + t6) * t153 + (t242 - 0.2e1 * t304 - 0.2e1 * t305) * t258 + (t53 + t32 - t90) * t117; -t305 - t304 + (-t104 + t49) * t186 + m(6) * (t109 * t11 + t110 * t12 + t147 * t60 + t186 * t22 + t30 * t80 + t31 * t81) + m(7) * (t1 * t44 + t120 * t42 + t13 * t4 + t14 * t3 + t143 * t15 + t2 * t43) + t242 + (t77 - t138) * t147 + (-t214 * t68 - t215 * t67 - t226 * t23 + t24 * t266) * mrSges(5,3) + ((-t251 / 0.2e1 - t212 / 0.2e1 + t210 / 0.2e1) * t264 + (Ifges(5,5) * t226 / 0.2e1 + Ifges(5,6) * t266 / 0.2e1 - Ifges(3,6) + Ifges(4,5) * t260 / 0.2e1 + Ifges(4,6) * t263 / 0.2e1) * t285) * t256 - (t35 / 0.2e1 + t6 / 0.2e1 - t63 / 0.2e1) * t266 + (-t167 / 0.2e1 + t124 / 0.2e1 + t45 / 0.2e1) * t153 + t135 * t332 + (-t178 / 0.2e1 + t134 / 0.2e1 + t83 / 0.2e1) * t117 + (-pkin(2) * t206 + (-t89 * t260 + t88 * t263 + (-t140 * t263 - t141 * t260) * qJD(3)) * pkin(9)) * m(4) + t250 * t70 + t200 * t229 + t218 * t230 / 0.2e1 + t219 * t231 / 0.2e1 + t184 * t238 / 0.2e1 + t185 * t239 / 0.2e1 + (t53 / 0.2e1 + t32 / 0.2e1 - t90 / 0.2e1) * t214 + t187 * t103 + t139 * t174 + t118 * t180 / 0.2e1 + t154 * t169 / 0.2e1 + t12 * t171 + t11 * t172 + t31 * t157 + t30 * t158 + t22 * t162 + t161 * t163 - t151 * t7 / 0.2e1 - t152 * t8 / 0.2e1 + t143 * t10 + t60 * t146 + t148 * t137 + t132 * t125 / 0.2e1 + t133 * t126 / 0.2e1 + t102 * t136 / 0.2e1 + t1 * t128 + t2 * t129 + t120 * t38 - pkin(2) * t121 + t109 * t59 + t110 * t58 + t96 * t34 / 0.2e1 + t97 * t33 / 0.2e1 + t15 * t98 + t29 * t84 / 0.2e1 + t28 * t85 / 0.2e1 + t81 * t86 + t80 * t87 + t76 * t47 / 0.2e1 + t3 * t78 + t4 * t79 + t75 * t46 / 0.2e1 + t42 * t48 + t13 * t50 + t14 * t51 + t43 * t17 + t44 * t18 + (-t206 * mrSges(4,1) + t107 / 0.2e1 + pkin(9) * t149 + t88 * mrSges(4,3)) * t263 + (t206 * mrSges(4,2) + t108 / 0.2e1 - pkin(9) * t150 - t89 * mrSges(4,3)) * t260 + ((-pkin(9) * t189 - t140 * mrSges(4,3) + t145 / 0.2e1) * t263 + (-pkin(9) * t188 - t141 * mrSges(4,3) + pkin(3) * t99 + t277 / 0.2e1 - t144 / 0.2e1) * t260) * qJD(3) + m(5) * (t139 * t250 - t147 * t67 + t148 * t68 + t161 * t278 - t186 * t23 + t187 * t24) + (t91 / 0.2e1 + t55 * t317 + t54 * t319) * t215 + (t64 / 0.2e1 + t37 * t317 + t36 * t319) * t226; -(t124 + t45 - t167) * t266 + 0.2e1 * (t147 * t226 + t148 * t266 + t186 * t215 - t187 * t214) * mrSges(5,3) + (t263 * t239 + (0.2e1 * pkin(3) * t174 - t238) * t260) * qJD(3) + (t148 * t187 + t250 * t278 + t297) * t338 + (t120 * t143 + t13 * t44 + t14 * t43) * t336 + (t109 * t80 + t110 * t81 + t297) * t337 + (-t125 * t254 + t126 * t257 + t169) * t226 + (-t135 * t254 + t136 * t257 + t180) * t215 + (t134 + t83 - t178) * t214 + t263 * t230 + t260 * t231 + 0.2e1 * t250 * t163 - 0.2e1 * pkin(2) * t229 + 0.2e1 * t186 * t146 + 0.2e1 * t81 * t171 + 0.2e1 * t80 * t172 + 0.2e1 * t110 * t157 + 0.2e1 * t109 * t158 + 0.2e1 * t147 * t162 - t151 * t46 - t152 * t47 + 0.2e1 * t143 * t48 + 0.2e1 * t13 * t128 + 0.2e1 * t14 * t129 + 0.2e1 * t120 * t98 + t96 * t85 + t97 * t84 + 0.2e1 * t43 * t78 + 0.2e1 * t44 * t79; t267 + t339 * t117 + t340 * mrSges(6,3) + m(6) * (t22 * t249 + t340 * t247 + (-t254 * t30 + t257 * t31) * qJD(5)) + t29 * t327 + t76 * t329 + t75 * t330 + t153 * t331 + t235 * t332 + (t23 * t298 + t24 * t255) * t333 + t36 * t317 + t37 * t318 + t102 * t320 + t8 * t322 + t7 * t323 + t33 * t324 + t34 * t325 + t28 * t326 - t4 * t302 + t103 * t315 + m(7) * (t1 * t160 + t130 * t4 + t131 * t3 + t15 * t232 + t159 * t2) + t249 * t49 + t232 * t10 + t22 * t233 + t15 * t175 + t42 * t164 + t159 * t17 + t160 * t18 + t130 * t50 + t131 * t51 - t88 * mrSges(4,2) + t89 * mrSges(4,1) - t24 * mrSges(5,2) + t23 * mrSges(5,1) + t104 * t276 + t86 * t281 + t58 * t291 - t1 * t301 - t87 * t282 - t59 * t292 + t3 * t303 - t2 * t312; m(7) * (t120 * t232 + t13 * t160 + t130 * t44 + t131 * t43 + t14 * t159) + t251 + (-mrSges(5,3) * t315 + t339) * t214 + m(6) * (t341 * t247 + (-t109 * t254 + t110 * t257) * qJD(5)) + t341 * mrSges(6,3) - t266 * t331 + (t255 * t333 - mrSges(5,2)) * t148 + (m(6) * t249 - t298 * t333 - mrSges(5,1) + t233) * t147 + t97 * t327 - t152 * t329 - t151 * t330 + (-mrSges(4,1) * t283 + mrSges(4,2) * t284) * pkin(9) - t215 * mrSges(5,3) * t276 + t125 * t317 + t126 * t318 + t295 * t320 + t47 * t322 + t46 * t323 + t84 * t324 + t85 * t325 + t96 * t326 - t210 + t212 + t249 * t146 + t232 * t48 + t120 * t175 + t143 * t164 + t159 * t78 + t160 * t79 + t130 * t128 + t131 * t129 + t171 * t281 + t157 * t291 - t13 * t301 - t44 * t302 - t172 * t282 - Ifges(4,6) * t284 - t158 * t292 - t235 * t296 / 0.2e1 + t43 * t303 - t14 * t312; (t130 * t160 + t131 * t159) * t336 - t216 * t179 + t227 * t168 - t217 * t177 - t268 * t166 + 0.2e1 * t232 * t164 + 0.2e1 * (-t130 * t268 - t131 * t227 + t159 * t216 - t160 * t217) * mrSges(7,3) + (t247 * t337 + 0.2e1 * mrSges(6,3)) * qJD(5) * (t254 ^ 2 + t257 ^ 2); -t268 * t17 + t227 * t18 - t216 * t50 - t217 * t51 + t254 * t58 + t257 * t59 + m(7) * (t1 * t227 - t2 * t268 - t216 * t4 - t217 * t3) + m(6) * (t11 * t257 + t12 * t254) + m(5) * t139 + t70; m(5) * t278 - t216 * t128 - t217 * t129 + t254 * t157 + t257 * t158 - t268 * t78 + t227 * t79 + m(7) * (t13 * t227 - t14 * t268 - t216 * t44 - t217 * t43) + m(6) * (t254 * t81 + t257 * t80) + t163; m(7) * (t130 * t227 - t131 * t268 - t159 * t217 - t160 * t216); (-t216 * t227 + t217 * t268) * t336; m(6) * t22 + m(7) * t15 + t10 + t49; m(6) * t147 + m(7) * t120 + t146 + t48; t164; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t6; mrSges(7,1) * t14 - mrSges(7,2) * t13 + t45; mrSges(7,1) * t131 - mrSges(7,2) * t130 + t165; -t164; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
