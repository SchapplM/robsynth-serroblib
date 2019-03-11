% Calculate time derivative of joint inertia matrix for
% S6RRPRRP10
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:36
% EndTime: 2019-03-09 12:36:53
% DurationCPUTime: 7.09s
% Computational Cost: add. (9444->637), mult. (24705->902), div. (0->0), fcn. (24714->10), ass. (0->256)
t230 = sin(qJ(5));
t335 = (Ifges(6,5) + Ifges(7,4)) * t230;
t227 = sin(pkin(6));
t334 = 0.2e1 * t227;
t226 = sin(pkin(11));
t228 = cos(pkin(11));
t231 = sin(qJ(4));
t234 = cos(qJ(4));
t186 = t226 * t234 + t228 * t231;
t233 = cos(qJ(5));
t278 = qJD(5) * t233;
t185 = t226 * t231 - t234 * t228;
t180 = t185 * qJD(4);
t292 = t230 * t180;
t237 = t186 * t278 - t292;
t279 = qJD(5) * t230;
t265 = t186 * t279;
t291 = t233 * t180;
t236 = t265 + t291;
t219 = -pkin(3) * t228 - pkin(2);
t135 = pkin(4) * t185 - pkin(10) * t186 + t219;
t314 = pkin(9) + qJ(3);
t201 = t314 * t226;
t202 = t314 * t228;
t149 = -t201 * t231 + t202 * t234;
t329 = t230 * t135 + t233 * t149;
t332 = qJD(5) * t329;
t229 = cos(pkin(6));
t232 = sin(qJ(2));
t295 = t227 * t232;
t178 = -t226 * t295 + t228 * t229;
t179 = t226 * t229 + t228 * t295;
t124 = t178 * t231 + t179 * t234;
t235 = cos(qJ(2));
t294 = t227 * t235;
t108 = t230 * t124 + t233 * t294;
t282 = qJD(2) * t227;
t267 = t232 * t282;
t238 = t234 * t178 - t179 * t231;
t266 = t235 * t282;
t95 = qJD(4) * t238 - t185 * t266;
t54 = -qJD(5) * t108 + t230 * t267 + t233 * t95;
t268 = t230 * t294;
t55 = -qJD(5) * t268 + t124 * t278 + t230 * t95 - t233 * t267;
t96 = qJD(4) * t124 + t186 * t266;
t7 = Ifges(6,5) * t54 - Ifges(6,6) * t55 + Ifges(6,3) * t96;
t8 = Ifges(7,4) * t54 + Ifges(7,2) * t96 + Ifges(7,6) * t55;
t331 = t7 + t8;
t328 = -t234 * t201 - t202 * t231;
t204 = -t233 * mrSges(6,1) + t230 * mrSges(6,2);
t327 = -m(6) * pkin(4) + t204;
t121 = -t185 * qJD(3) + qJD(4) * t328;
t181 = t186 * qJD(4);
t134 = pkin(4) * t181 + pkin(10) * t180;
t33 = -t121 * t230 + t134 * t233 - t332;
t183 = t229 * t232 * pkin(1) + pkin(8) * t294;
t165 = qJ(3) * t229 + t183;
t166 = (-pkin(2) * t235 - qJ(3) * t232 - pkin(1)) * t227;
t120 = t228 * t165 + t226 * t166;
t100 = pkin(9) * t178 + t120;
t280 = qJD(4) * t234;
t281 = qJD(4) * t231;
t153 = (-qJD(3) * t232 + (pkin(2) * t232 - qJ(3) * t235) * qJD(2)) * t227;
t315 = pkin(1) * t235;
t275 = t229 * t315;
t170 = -pkin(8) * t267 + qJD(2) * t275;
t158 = qJD(3) * t229 + t170;
t106 = t228 * t153 - t226 * t158;
t293 = t228 * t235;
t82 = (pkin(3) * t232 - pkin(9) * t293) * t282 + t106;
t119 = -t226 * t165 + t228 * t166;
t84 = -pkin(3) * t294 - t179 * pkin(9) + t119;
t107 = t226 * t153 + t228 * t158;
t251 = t226 * t266;
t99 = -pkin(9) * t251 + t107;
t17 = -t100 * t281 + t231 * t82 + t234 * t99 + t84 * t280;
t15 = pkin(10) * t267 + t17;
t171 = t183 * qJD(2);
t147 = pkin(3) * t251 + t171;
t31 = t96 * pkin(4) - t95 * pkin(10) + t147;
t46 = t234 * t100 + t231 * t84;
t36 = -pkin(10) * t294 + t46;
t214 = pkin(8) * t295;
t168 = t214 + (-pkin(2) - t315) * t229;
t128 = -t178 * pkin(3) + t168;
t58 = -pkin(4) * t238 - t124 * pkin(10) + t128;
t311 = t230 * t58 + t233 * t36;
t4 = -qJD(5) * t311 - t15 * t230 + t233 * t31;
t244 = pkin(5) * t233 + qJ(6) * t230;
t277 = qJD(6) * t233;
t326 = qJD(5) * t244 - t277;
t325 = 2 * m(4);
t324 = 2 * m(5);
t323 = 0.2e1 * m(6);
t322 = 2 * m(7);
t321 = -2 * mrSges(3,3);
t320 = -2 * mrSges(5,3);
t122 = qJD(3) * t186 + qJD(4) * t149;
t319 = 0.2e1 * t122;
t318 = -0.2e1 * t328;
t316 = t228 / 0.2e1;
t24 = mrSges(6,1) * t96 - mrSges(6,3) * t54;
t25 = -t96 * mrSges(7,1) + t54 * mrSges(7,2);
t313 = t24 - t25;
t26 = -mrSges(6,2) * t96 - mrSges(6,3) * t55;
t27 = -mrSges(7,2) * t55 + mrSges(7,3) * t96;
t312 = t26 + t27;
t61 = mrSges(6,2) * t238 - mrSges(6,3) * t108;
t62 = -mrSges(7,2) * t108 - mrSges(7,3) * t238;
t310 = t61 + t62;
t109 = t233 * t124 - t268;
t63 = -mrSges(6,1) * t238 - mrSges(6,3) * t109;
t64 = mrSges(7,1) * t238 + mrSges(7,2) * t109;
t309 = -t63 + t64;
t308 = Ifges(4,4) * t226;
t307 = Ifges(4,4) * t228;
t306 = Ifges(6,4) * t230;
t305 = Ifges(6,4) * t233;
t304 = Ifges(7,5) * t230;
t303 = Ifges(7,5) * t233;
t302 = Ifges(6,6) * t233;
t301 = t147 * mrSges(5,1);
t300 = t147 * mrSges(5,2);
t299 = t122 * t328;
t298 = t186 * t230;
t297 = t186 * t233;
t102 = mrSges(6,1) * t181 + mrSges(6,3) * t236;
t103 = -t181 * mrSges(7,1) - mrSges(7,2) * t236;
t290 = t102 - t103;
t104 = -mrSges(6,2) * t181 - mrSges(6,3) * t237;
t105 = -mrSges(7,2) * t237 + mrSges(7,3) * t181;
t289 = t104 + t105;
t245 = Ifges(7,3) * t230 + t303;
t110 = Ifges(7,6) * t185 + t186 * t245;
t246 = -Ifges(6,2) * t230 + t305;
t113 = Ifges(6,6) * t185 + t186 * t246;
t288 = t110 - t113;
t247 = Ifges(7,1) * t233 + t304;
t114 = Ifges(7,4) * t185 + t186 * t247;
t248 = Ifges(6,1) * t233 - t306;
t115 = Ifges(6,5) * t185 + t186 * t248;
t287 = t114 + t115;
t136 = -mrSges(6,2) * t185 - mrSges(6,3) * t298;
t139 = -mrSges(7,2) * t298 + mrSges(7,3) * t185;
t286 = t136 + t139;
t137 = mrSges(6,1) * t185 - mrSges(6,3) * t297;
t138 = -mrSges(7,1) * t185 + mrSges(7,2) * t297;
t285 = -t137 + t138;
t284 = -Ifges(6,5) * t291 + Ifges(6,3) * t181;
t283 = -Ifges(5,5) * t180 - Ifges(5,6) * t181;
t154 = t228 * mrSges(4,2) * t266 + mrSges(4,1) * t251;
t191 = Ifges(7,4) * t278 + Ifges(7,6) * t279;
t6 = Ifges(7,5) * t54 + Ifges(7,6) * t96 + Ifges(7,3) * t55;
t9 = Ifges(6,4) * t54 - Ifges(6,2) * t55 + Ifges(6,6) * t96;
t276 = t6 / 0.2e1 - t9 / 0.2e1;
t274 = Ifges(5,5) * t95 - Ifges(5,6) * t96 + Ifges(5,3) * t267;
t10 = Ifges(7,1) * t54 + Ifges(7,4) * t96 + Ifges(7,5) * t55;
t11 = Ifges(6,1) * t54 - Ifges(6,4) * t55 + Ifges(6,5) * t96;
t273 = t10 / 0.2e1 + t11 / 0.2e1;
t37 = Ifges(7,5) * t109 - Ifges(7,6) * t238 + Ifges(7,3) * t108;
t40 = Ifges(6,4) * t109 - Ifges(6,2) * t108 - Ifges(6,6) * t238;
t272 = t37 / 0.2e1 - t40 / 0.2e1;
t41 = Ifges(7,1) * t109 - Ifges(7,4) * t238 + Ifges(7,5) * t108;
t42 = Ifges(6,1) * t109 - Ifges(6,4) * t108 - Ifges(6,5) * t238;
t271 = -t42 / 0.2e1 - t41 / 0.2e1;
t65 = -Ifges(7,5) * t236 + Ifges(7,6) * t181 + Ifges(7,3) * t237;
t68 = -Ifges(6,4) * t236 - Ifges(6,2) * t237 + Ifges(6,6) * t181;
t270 = t65 / 0.2e1 - t68 / 0.2e1;
t69 = -Ifges(7,1) * t236 + Ifges(7,4) * t181 + Ifges(7,5) * t237;
t70 = -Ifges(6,1) * t236 - Ifges(6,4) * t237 + Ifges(6,5) * t181;
t269 = t69 / 0.2e1 + t70 / 0.2e1;
t263 = t110 / 0.2e1 - t113 / 0.2e1;
t262 = t114 / 0.2e1 + t115 / 0.2e1;
t190 = Ifges(6,5) * t278 - Ifges(6,6) * t279;
t261 = t190 / 0.2e1 + t191 / 0.2e1;
t189 = t245 * qJD(5);
t192 = t246 * qJD(5);
t260 = -t192 / 0.2e1 + t189 / 0.2e1;
t193 = t247 * qJD(5);
t194 = t248 * qJD(5);
t259 = t193 / 0.2e1 + t194 / 0.2e1;
t205 = -Ifges(7,3) * t233 + t304;
t208 = Ifges(6,2) * t233 + t306;
t258 = t205 / 0.2e1 - t208 / 0.2e1;
t257 = t302 / 0.2e1 - Ifges(7,6) * t233 / 0.2e1 + t335 / 0.2e1;
t209 = Ifges(7,1) * t230 - t303;
t210 = Ifges(6,1) * t230 + t305;
t256 = t209 / 0.2e1 + t210 / 0.2e1;
t47 = t96 * mrSges(5,1) + t95 * mrSges(5,2);
t45 = -t231 * t100 + t234 * t84;
t129 = t181 * mrSges(5,1) - t180 * mrSges(5,2);
t254 = Ifges(5,5) * t267;
t253 = Ifges(5,6) * t267;
t252 = -Ifges(7,4) * t291 + Ifges(7,2) * t181 + Ifges(7,6) * t237;
t35 = pkin(4) * t294 - t45;
t250 = mrSges(6,1) * t230 + mrSges(6,2) * t233;
t203 = -t233 * mrSges(7,1) - t230 * mrSges(7,3);
t249 = mrSges(7,1) * t230 - mrSges(7,3) * t233;
t243 = pkin(5) * t230 - qJ(6) * t233;
t21 = -t230 * t36 + t233 * t58;
t18 = -t100 * t280 - t231 * t99 + t234 * t82 - t84 * t281;
t92 = t135 * t233 - t149 * t230;
t3 = t233 * t15 + t230 * t31 + t58 * t278 - t279 * t36;
t32 = t233 * t121 + t230 * t134 + t135 * t278 - t149 * t279;
t16 = -pkin(4) * t267 - t18;
t212 = Ifges(3,5) * t266;
t198 = -pkin(4) - t244;
t188 = t250 * qJD(5);
t187 = t249 * qJD(5);
t182 = -t214 + t275;
t177 = -pkin(5) * t279 + qJ(6) * t278 + qJD(6) * t230;
t160 = (mrSges(4,1) * t232 - mrSges(4,3) * t293) * t282;
t159 = (-mrSges(4,3) * t226 * t235 - mrSges(4,2) * t232) * t282;
t152 = -mrSges(4,1) * t294 - t179 * mrSges(4,3);
t151 = mrSges(4,2) * t294 + t178 * mrSges(4,3);
t144 = Ifges(5,1) * t186 - Ifges(5,4) * t185;
t143 = Ifges(5,4) * t186 - Ifges(5,2) * t185;
t141 = (t232 * Ifges(4,5) + (t228 * Ifges(4,1) - t308) * t235) * t282;
t140 = (t232 * Ifges(4,6) + (-t226 * Ifges(4,2) + t307) * t235) * t282;
t133 = t250 * t186;
t132 = t249 * t186;
t131 = -Ifges(5,1) * t180 - Ifges(5,4) * t181;
t130 = -Ifges(5,4) * t180 - Ifges(5,2) * t181;
t117 = -mrSges(5,1) * t294 - t124 * mrSges(5,3);
t116 = mrSges(5,2) * t294 + mrSges(5,3) * t238;
t112 = Ifges(7,2) * t185 + (Ifges(7,4) * t233 + Ifges(7,6) * t230) * t186;
t111 = Ifges(6,3) * t185 + (Ifges(6,5) * t233 - Ifges(6,6) * t230) * t186;
t101 = t186 * t243 - t328;
t79 = mrSges(6,1) * t237 - mrSges(6,2) * t236;
t78 = mrSges(7,1) * t237 + mrSges(7,3) * t236;
t76 = -mrSges(5,2) * t267 - mrSges(5,3) * t96;
t75 = mrSges(5,1) * t267 - mrSges(5,3) * t95;
t74 = -pkin(5) * t185 - t92;
t73 = qJ(6) * t185 + t329;
t72 = Ifges(5,1) * t124 + Ifges(5,4) * t238 - Ifges(5,5) * t294;
t71 = Ifges(5,4) * t124 + Ifges(5,2) * t238 - Ifges(5,6) * t294;
t67 = -Ifges(7,4) * t265 + t252;
t66 = -Ifges(6,5) * t265 - Ifges(6,6) * t237 + t284;
t60 = mrSges(6,1) * t108 + mrSges(6,2) * t109;
t59 = mrSges(7,1) * t108 - mrSges(7,3) * t109;
t48 = -t243 * t180 + t186 * t326 + t122;
t44 = Ifges(5,1) * t95 - Ifges(5,4) * t96 + t254;
t43 = Ifges(5,4) * t95 - Ifges(5,2) * t96 + t253;
t39 = Ifges(7,4) * t109 - Ifges(7,2) * t238 + Ifges(7,6) * t108;
t38 = Ifges(6,5) * t109 - Ifges(6,6) * t108 - Ifges(6,3) * t238;
t29 = -pkin(5) * t181 - t33;
t28 = qJ(6) * t181 + qJD(6) * t185 + t32;
t23 = pkin(5) * t108 - qJ(6) * t109 + t35;
t20 = mrSges(6,1) * t55 + mrSges(6,2) * t54;
t19 = mrSges(7,1) * t55 - mrSges(7,3) * t54;
t13 = pkin(5) * t238 - t21;
t12 = -qJ(6) * t238 + t311;
t5 = pkin(5) * t55 - qJ(6) * t54 - qJD(6) * t109 + t16;
t2 = -pkin(5) * t96 - t4;
t1 = qJ(6) * t96 - qJD(6) * t238 + t3;
t14 = [(t10 + t11) * t109 + (t6 - t9) * t108 + ((t183 * t321 + Ifges(4,5) * t179 + Ifges(5,5) * t124 - 0.2e1 * Ifges(3,6) * t229 + Ifges(4,6) * t178 + Ifges(5,6) * t238 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t232) * t334) * t232 + (Ifges(3,5) * t229 + t182 * t321 - t226 * (Ifges(4,4) * t179 + Ifges(4,2) * t178) + t228 * (Ifges(4,1) * t179 + Ifges(4,4) * t178) + (-pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t228 + Ifges(4,6) * t226 + Ifges(3,4)) * t235) * t334 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3)) * t295) * t235) * t282 + 0.2e1 * t311 * t26 + (t16 * t35 + t21 * t4 + t3 * t311) * t323 + t229 * t212 + t178 * t140 + 0.2e1 * t171 * (-mrSges(4,1) * t178 + mrSges(4,2) * t179) + t179 * t141 + 0.2e1 * t168 * t154 + 0.2e1 * t120 * t159 + 0.2e1 * t119 * t160 + 0.2e1 * t107 * t151 + 0.2e1 * t106 * t152 + 0.2e1 * t128 * t47 + 0.2e1 * t17 * t116 + 0.2e1 * t18 * t117 + t95 * t72 + 0.2e1 * t45 * t75 + 0.2e1 * t46 * t76 + 0.2e1 * t4 * t63 + 0.2e1 * t2 * t64 + 0.2e1 * t5 * t59 + 0.2e1 * t16 * t60 + 0.2e1 * t3 * t61 + 0.2e1 * t1 * t62 + 0.2e1 * t35 * t20 + 0.2e1 * t13 * t25 + 0.2e1 * t12 * t27 + 0.2e1 * t23 * t19 + 0.2e1 * t21 * t24 + 0.2e1 * m(3) * (t170 * t183 - t171 * t182) + (-t71 + t38 + t39) * t96 + (t37 - t40) * t55 + (t41 + t42) * t54 + (t1 * t12 + t13 * t2 + t23 * t5) * t322 + (t128 * t147 + t17 * t46 + t18 * t45) * t324 + (t106 * t119 + t107 * t120 + t168 * t171) * t325 - (-t43 + 0.2e1 * t301 + t331) * t238 - t274 * t294 + 0.2e1 * t170 * (-t229 * mrSges(3,2) + mrSges(3,3) * t294) - 0.2e1 * t171 * (mrSges(3,1) * t229 - mrSges(3,3) * t295) + (t44 + 0.2e1 * t300) * t124; t212 + m(5) * (t121 * t46 - t122 * t45 + t147 * t219 + t149 * t17 + t18 * t328) + m(6) * (t122 * t35 - t16 * t328 + t21 * t33 + t3 * t329 + t311 * t32 + t4 * t92) - (-t75 + t20) * t328 + t311 * t104 + m(7) * (t1 * t73 + t101 * t5 + t12 * t28 + t13 * t29 + t2 * t74 + t23 * t48) + (-t143 / 0.2e1 + t111 / 0.2e1 + t112 / 0.2e1) * t96 + (qJD(3) * t151 + qJ(3) * t159 + t107 * mrSges(4,3) - t171 * mrSges(4,1) + t140 / 0.2e1) * t228 + t219 * t47 - t170 * mrSges(3,2) - t171 * mrSges(3,1) + t149 * t76 - pkin(2) * t154 + t3 * t136 + t4 * t137 + t2 * t138 + t1 * t139 + t95 * t144 / 0.2e1 + (-qJD(3) * t152 - t106 * mrSges(4,3) - qJ(3) * t160 + t171 * mrSges(4,2) + t141 / 0.2e1) * t226 + t128 * t129 + t124 * t131 / 0.2e1 + t5 * t132 + t16 * t133 + t121 * t116 + t101 * t19 + t21 * t102 + t13 * t103 + t12 * t105 + t92 * t24 + t35 * t79 + t73 * t27 + t74 * t25 + t23 * t78 + t33 * t63 + t29 * t64 + t48 * t59 + t32 * t61 + t28 * t62 + m(4) * (-pkin(2) * t171 + (-t119 * t226 + t120 * t228) * qJD(3) + (-t106 * t226 + t107 * t228) * qJ(3)) + (t60 - t117) * t122 + (-t46 * mrSges(5,3) + t38 / 0.2e1 + t39 / 0.2e1 - t71 / 0.2e1) * t181 + t262 * t54 + t263 * t55 + t269 * t109 + t270 * t108 - (-t45 * mrSges(5,3) + t72 / 0.2e1 - t271 * t233 + t272 * t230) * t180 + t329 * t26 - (-t130 / 0.2e1 + t66 / 0.2e1 + t67 / 0.2e1) * t238 + (-t18 * mrSges(5,3) + t254 / 0.2e1 + t300 + t44 / 0.2e1 + t273 * t233 + t276 * t230 + (t230 * t271 + t233 * t272) * qJD(5)) * t186 + (-t17 * mrSges(5,3) - t253 / 0.2e1 - t43 / 0.2e1 + t301 + t8 / 0.2e1 + t7 / 0.2e1) * t185 + (-t235 * t283 / 0.2e1 + ((-Ifges(3,6) + Ifges(4,5) * t226 / 0.2e1 + Ifges(4,6) * t316) * t232 + (-t226 * (Ifges(4,2) * t228 + t308) / 0.2e1 + (Ifges(4,1) * t226 + t307) * t316) * t235) * qJD(2)) * t227; 0.2e1 * t101 * t78 + 0.2e1 * t92 * t102 + 0.2e1 * t74 * t103 + 0.2e1 * t329 * t104 + 0.2e1 * t73 * t105 + t133 * t319 + 0.2e1 * t219 * t129 + 0.2e1 * t48 * t132 + 0.2e1 * t32 * t136 + 0.2e1 * t33 * t137 + 0.2e1 * t29 * t138 + 0.2e1 * t28 * t139 + t79 * t318 + (t121 * t320 - t130 + t66 + t67) * t185 + (t149 * t320 + t111 + t112 - t143) * t181 + (t32 * t329 + t33 * t92 - t299) * t323 + (t121 * t149 - t299) * t324 + (t101 * t48 + t28 * t73 + t29 * t74) * t322 - (mrSges(5,3) * t318 + t230 * t288 + t233 * t287 + t144) * t180 + (mrSges(5,3) * t319 + t131 + (t69 + t70) * t233 + (t65 - t68) * t230 + (-t230 * t287 + t233 * t288) * qJD(5)) * t186 + (qJ(3) * t325 + 0.2e1 * mrSges(4,3)) * qJD(3) * (t226 ^ 2 + t228 ^ 2); t313 * t233 + t312 * t230 + (t230 * t309 + t233 * t310) * qJD(5) + m(6) * (t230 * t3 + t233 * t4 + (-t21 * t230 + t233 * t311) * qJD(5)) + m(7) * (t1 * t230 - t2 * t233 + (t12 * t233 + t13 * t230) * qJD(5)) + m(5) * t147 + m(4) * t171 + t47 + t154; t290 * t233 + t289 * t230 + (t285 * t230 + t286 * t233) * qJD(5) + m(7) * (t230 * t28 - t233 * t29 + (t230 * t74 + t233 * t73) * qJD(5)) + m(6) * (t230 * t32 + t233 * t33 + (-t230 * t92 + t233 * t329) * qJD(5)) + t129; 0; t274 + (t312 * t233 - t313 * t230 + (-t230 * t310 + t233 * t309) * qJD(5) + m(7) * (t1 * t233 - t12 * t279 + t13 * t278 + t2 * t230) + m(6) * (-t21 * t278 - t230 * t4 + t233 * t3 - t279 * t311)) * pkin(10) + t5 * t203 + t198 * t19 + t23 * t187 + t35 * t188 - t177 * t59 - pkin(4) * t20 - t17 * mrSges(5,2) + t18 * mrSges(5,1) + (t1 * mrSges(7,2) + t3 * mrSges(6,3) - t276) * t233 + m(7) * (-t177 * t23 + t198 * t5) + t256 * t54 + t257 * t96 + t258 * t55 + t259 * t109 + t260 * t108 - t261 * t238 + ((t13 * mrSges(7,2) - t21 * mrSges(6,3) - t271) * t233 + (-t12 * mrSges(7,2) - mrSges(6,3) * t311 + t272) * t230) * qJD(5) + (t2 * mrSges(7,2) - t4 * mrSges(6,3) + t273) * t230 + t327 * t16; t48 * t203 + t198 * t78 + t101 * t187 - t328 * t188 - t177 * t132 - t121 * mrSges(5,2) - pkin(4) * t79 + m(7) * (-t101 * t177 + t198 * t48) + t261 * t185 + t257 * t181 + (-mrSges(5,1) + t327) * t122 + (t32 * mrSges(6,3) + t28 * mrSges(7,2) + t259 * t186 - t256 * t180 + (t74 * mrSges(7,2) - t92 * mrSges(6,3) + t258 * t186 + t262) * qJD(5) + (t285 * qJD(5) + m(7) * (qJD(5) * t74 + t28) + m(6) * (-qJD(5) * t92 + t32) + t289) * pkin(10) - t270) * t233 + (t29 * mrSges(7,2) - t33 * mrSges(6,3) + t260 * t186 - t258 * t180 + (-t73 * mrSges(7,2) - mrSges(6,3) * t329 - t256 * t186 + t263) * qJD(5) + (-t286 * qJD(5) + m(7) * (-qJD(5) * t73 + t29) + m(6) * (-t33 - t332) - t290) * pkin(10) + t269) * t230 + t283; 0; -0.2e1 * pkin(4) * t188 + 0.2e1 * t187 * t198 + (-t189 + t192) * t233 + (t193 + t194) * t230 + 0.2e1 * (-m(7) * t198 - t203) * t177 + ((t209 + t210) * t233 + (t205 - t208) * t230) * qJD(5); t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) - pkin(5) * t25 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) + qJD(6) * t62 + qJ(6) * t27 + t1 * mrSges(7,3) + t331; Ifges(6,6) * t292 - pkin(5) * t103 + m(7) * (-pkin(5) * t29 + qJ(6) * t28 + qJD(6) * t73) + qJD(6) * t139 + qJ(6) * t105 + t28 * mrSges(7,3) + t33 * mrSges(6,1) - t29 * mrSges(7,1) - t32 * mrSges(6,2) + (-t302 - t335) * t186 * qJD(5) + t252 + t284; m(7) * t177 + ((-mrSges(6,2) + mrSges(7,3)) * t233 + (-mrSges(6,1) - mrSges(7,1)) * t230) * qJD(5); -t326 * mrSges(7,2) + (m(7) * t277 + (-m(7) * t244 + t203 + t204) * qJD(5)) * pkin(10) + t190 + t191; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t25; m(7) * t29 + t103; m(7) * t279; (m(7) * pkin(10) + mrSges(7,2)) * t278; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
