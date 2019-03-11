% Calculate time derivative of joint inertia matrix for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:22
% EndTime: 2019-03-09 21:30:39
% DurationCPUTime: 7.07s
% Computational Cost: add. (6325->760), mult. (17255->1018), div. (0->0), fcn. (15290->8), ass. (0->292)
t359 = Ifges(5,5) + Ifges(6,4);
t250 = cos(qJ(4));
t251 = cos(qJ(3));
t312 = qJD(3) * t251;
t291 = t250 * t312;
t247 = sin(qJ(4));
t248 = sin(qJ(3));
t310 = qJD(4) * t248;
t294 = t247 * t310;
t257 = t291 - t294;
t309 = qJD(4) * t250;
t256 = t247 * t312 + t248 * t309;
t358 = Ifges(5,6) / 0.2e1;
t245 = sin(pkin(6));
t249 = sin(qJ(2));
t323 = t245 * t249;
t229 = pkin(8) * t323;
t246 = cos(pkin(6));
t252 = cos(qJ(2));
t346 = pkin(1) * t252;
t168 = t246 * t346 - t229;
t313 = qJD(3) * t248;
t356 = qJ(5) * t313 - qJD(5) * t251;
t166 = t246 * t248 + t251 * t323;
t315 = qJD(2) * t245;
t295 = t252 * t315;
t115 = qJD(3) * t166 + t248 * t295;
t165 = -t246 * t251 + t248 * t323;
t116 = -qJD(3) * t165 + t251 * t295;
t314 = qJD(2) * t249;
t296 = t245 * t314;
t322 = t245 * t252;
t301 = t247 * t322;
t60 = -qJD(4) * t301 + t116 * t247 + t166 * t309 - t250 * t296;
t117 = t166 * t247 + t250 * t322;
t61 = -qJD(4) * t117 + t116 * t250 + t247 * t296;
t11 = Ifges(5,5) * t61 - Ifges(5,6) * t60 + Ifges(5,3) * t115;
t13 = Ifges(6,4) * t61 + Ifges(6,2) * t115 + Ifges(6,6) * t60;
t9 = Ifges(7,5) * t61 + Ifges(7,6) * t60 - Ifges(7,3) * t115;
t355 = t11 + t13 - t9;
t325 = qJ(5) * t247;
t347 = pkin(4) + pkin(5);
t354 = -t250 * t347 - t325;
t353 = 2 * m(5);
t352 = 2 * m(6);
t351 = 2 * m(7);
t350 = 2 * pkin(9);
t349 = -2 * mrSges(3,3);
t348 = -2 * mrSges(7,3);
t345 = pkin(9) * t247;
t344 = pkin(9) * t251;
t343 = pkin(10) * t247;
t342 = pkin(10) * t250;
t341 = -mrSges(6,2) + mrSges(7,3);
t340 = -Ifges(5,6) - Ifges(7,6);
t339 = pkin(10) - qJ(6);
t144 = t229 + (-pkin(2) - t346) * t246;
t76 = pkin(3) * t165 - pkin(10) * t166 + t144;
t169 = t246 * t249 * pkin(1) + pkin(8) * t322;
t145 = pkin(9) * t246 + t169;
t146 = (-pkin(2) * t252 - pkin(9) * t249 - pkin(1)) * t245;
t90 = t251 * t145 + t248 * t146;
t78 = -pkin(10) * t322 + t90;
t27 = t247 * t76 + t250 * t78;
t338 = Ifges(4,4) * t248;
t337 = Ifges(4,4) * t251;
t336 = Ifges(5,4) * t247;
t335 = Ifges(5,4) * t250;
t334 = Ifges(7,4) * t247;
t333 = Ifges(7,4) * t250;
t332 = Ifges(6,5) * t247;
t331 = Ifges(6,5) * t250;
t330 = Ifges(7,5) * t250;
t159 = t168 * qJD(2);
t329 = t159 * mrSges(3,2);
t160 = t169 * qJD(2);
t328 = t160 * mrSges(3,1);
t327 = t160 * mrSges(4,1);
t326 = t160 * mrSges(4,2);
t324 = qJ(5) * t250;
t321 = t247 * t248;
t320 = t248 * t250;
t319 = t250 * t251;
t89 = -t248 * t145 + t251 * t146;
t198 = (pkin(3) * t248 - pkin(10) * t251) * qJD(3);
t202 = -pkin(3) * t251 - pkin(10) * t248 - pkin(2);
t318 = t247 * t198 + t202 * t309;
t234 = pkin(9) * t319;
t311 = qJD(4) * t247;
t317 = qJD(4) * t234 + t202 * t311;
t316 = Ifges(5,5) * t291 + Ifges(5,3) * t313;
t143 = t247 * t202 + t234;
t185 = Ifges(6,4) * t309 + Ifges(6,6) * t311;
t308 = qJD(5) * t247;
t307 = qJD(5) * t250;
t305 = qJD(6) * t250;
t24 = t165 * qJ(5) + t27;
t303 = Ifges(4,6) * t322;
t302 = mrSges(7,3) * t319;
t300 = Ifges(4,5) * t116 - Ifges(4,6) * t115 + Ifges(4,3) * t296;
t77 = pkin(3) * t322 - t89;
t266 = Ifges(6,3) * t247 + t331;
t148 = -Ifges(6,6) * t251 + t248 * t266;
t267 = Ifges(7,2) * t247 + t333;
t150 = Ifges(7,6) * t251 + t248 * t267;
t268 = -Ifges(5,2) * t247 + t335;
t152 = -Ifges(5,6) * t251 + t248 * t268;
t299 = t148 + t150 - t152;
t269 = Ifges(7,1) * t250 + t334;
t153 = Ifges(7,5) * t251 + t248 * t269;
t270 = Ifges(6,1) * t250 + t332;
t154 = -Ifges(6,4) * t251 + t248 * t270;
t271 = Ifges(5,1) * t250 - t336;
t155 = -Ifges(5,5) * t251 + t248 * t271;
t298 = t153 + t154 + t155;
t297 = -pkin(4) - t345;
t21 = -t60 * mrSges(7,1) + t61 * mrSges(7,2);
t207 = t339 * t250;
t33 = -t115 * mrSges(6,1) + t61 * mrSges(6,2);
t31 = -t115 * mrSges(7,1) - t61 * mrSges(7,3);
t26 = -t247 * t78 + t250 * t76;
t233 = t247 * t344;
t142 = t202 * t250 - t233;
t158 = (pkin(2) * t249 - pkin(9) * t252) * t315;
t40 = -t145 * t312 - t146 * t313 + t251 * t158 - t248 * t159;
t290 = Ifges(6,4) * t291 + Ifges(6,2) * t313 + Ifges(6,6) * t256;
t289 = t296 / 0.2e1;
t10 = Ifges(6,5) * t61 + Ifges(6,6) * t115 + Ifges(6,3) * t60;
t12 = Ifges(7,4) * t61 + Ifges(7,2) * t60 - Ifges(7,6) * t115;
t14 = Ifges(5,4) * t61 - Ifges(5,2) * t60 + Ifges(5,6) * t115;
t288 = t10 / 0.2e1 + t12 / 0.2e1 - t14 / 0.2e1;
t15 = Ifges(7,1) * t61 + Ifges(7,4) * t60 - Ifges(7,5) * t115;
t16 = Ifges(6,1) * t61 + Ifges(6,4) * t115 + Ifges(6,5) * t60;
t17 = Ifges(5,1) * t61 - Ifges(5,4) * t60 + Ifges(5,5) * t115;
t287 = t15 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1;
t118 = t166 * t250 - t301;
t45 = Ifges(6,5) * t118 + Ifges(6,6) * t165 + Ifges(6,3) * t117;
t47 = Ifges(7,4) * t118 + Ifges(7,2) * t117 - Ifges(7,6) * t165;
t49 = Ifges(5,4) * t118 - Ifges(5,2) * t117 + Ifges(5,6) * t165;
t286 = t47 / 0.2e1 - t49 / 0.2e1 + t45 / 0.2e1;
t50 = Ifges(7,1) * t118 + Ifges(7,4) * t117 - Ifges(7,5) * t165;
t51 = Ifges(6,1) * t118 + Ifges(6,4) * t165 + Ifges(6,5) * t117;
t52 = Ifges(5,1) * t118 - Ifges(5,4) * t117 + Ifges(5,5) * t165;
t285 = t50 / 0.2e1 + t51 / 0.2e1 + t52 / 0.2e1;
t213 = Ifges(5,2) * t250 + t336;
t100 = -t213 * t310 + (Ifges(5,6) * t248 + t251 * t268) * qJD(3);
t209 = -Ifges(6,3) * t250 + t332;
t96 = -t209 * t310 + (Ifges(6,6) * t248 + t251 * t266) * qJD(3);
t211 = -Ifges(7,2) * t250 + t334;
t98 = -t211 * t310 + (-Ifges(7,6) * t248 + t251 * t267) * qJD(3);
t284 = t96 / 0.2e1 + t98 / 0.2e1 - t100 / 0.2e1;
t178 = -mrSges(7,1) * t311 + mrSges(7,2) * t309;
t122 = -qJ(5) * t251 + t143;
t283 = -t198 * t250 + t317;
t215 = Ifges(7,1) * t247 - t333;
t101 = -t215 * t310 + (-Ifges(7,5) * t248 + t251 * t269) * qJD(3);
t216 = Ifges(6,1) * t247 - t331;
t102 = -t216 * t310 + (Ifges(6,4) * t248 + t251 * t270) * qJD(3);
t217 = Ifges(5,1) * t247 + t335;
t103 = -t217 * t310 + (Ifges(5,5) * t248 + t251 * t271) * qJD(3);
t282 = t101 / 0.2e1 + t102 / 0.2e1 + t103 / 0.2e1;
t281 = t148 / 0.2e1 + t150 / 0.2e1 - t152 / 0.2e1;
t280 = t153 / 0.2e1 + t154 / 0.2e1 + t155 / 0.2e1;
t241 = Ifges(5,5) * t309;
t265 = Ifges(7,6) * t247 + t330;
t279 = t265 * qJD(4) / 0.2e1 + t311 * t358 - t241 / 0.2e1 - t185 / 0.2e1;
t182 = t266 * qJD(4);
t184 = t267 * qJD(4);
t186 = t268 * qJD(4);
t278 = t182 / 0.2e1 + t184 / 0.2e1 - t186 / 0.2e1;
t188 = t269 * qJD(4);
t189 = t270 * qJD(4);
t190 = t271 * qJD(4);
t277 = t188 / 0.2e1 + t189 / 0.2e1 + t190 / 0.2e1;
t208 = Ifges(7,5) * t247 - Ifges(7,6) * t250;
t276 = -t208 / 0.2e1 + t359 * t247 / 0.2e1 + (t358 - Ifges(6,6) / 0.2e1) * t250;
t275 = t209 / 0.2e1 + t211 / 0.2e1 - t213 / 0.2e1;
t274 = t215 / 0.2e1 + t216 / 0.2e1 + t217 / 0.2e1;
t206 = -t250 * mrSges(5,1) + t247 * mrSges(5,2);
t273 = mrSges(5,1) * t247 + mrSges(5,2) * t250;
t204 = -t250 * mrSges(6,1) - t247 * mrSges(6,3);
t272 = mrSges(6,1) * t247 - mrSges(6,3) * t250;
t264 = pkin(4) * t250 + t325;
t263 = pkin(4) * t247 - t324;
t39 = -t145 * t313 + t146 * t312 + t248 * t158 + t251 * t159;
t37 = pkin(10) * t296 + t39;
t53 = pkin(3) * t115 - pkin(10) * t116 + t160;
t7 = -t247 * t37 + t250 * t53 - t78 * t309 - t76 * t311;
t262 = t250 * (m(6) * pkin(10) - t341);
t261 = pkin(9) + t263;
t260 = qJ(5) * t118 - t77;
t259 = -t247 * t347 + t324;
t6 = t247 * t53 + t250 * t37 + t76 * t309 - t311 * t78;
t258 = -pkin(9) + t259;
t4 = t115 * qJ(5) + t165 * qJD(5) + t6;
t255 = pkin(3) * t296 + t40;
t87 = (-t250 * t313 - t251 * t311) * pkin(9) + t318;
t254 = qJ(5) * t61 + qJD(5) * t118 + t255;
t106 = -mrSges(7,1) * t256 + mrSges(7,2) * t257;
t244 = t251 * pkin(4);
t242 = Ifges(4,5) * t312;
t223 = mrSges(7,3) * t294;
t222 = mrSges(6,2) * t291;
t220 = Ifges(3,5) * t295;
t218 = Ifges(4,1) * t248 + t337;
t214 = Ifges(4,2) * t251 + t338;
t205 = mrSges(7,1) * t250 + mrSges(7,2) * t247;
t203 = t339 * t247;
t199 = -pkin(3) - t264;
t197 = -mrSges(6,2) * t321 - mrSges(6,3) * t251;
t196 = mrSges(6,1) * t251 + mrSges(6,2) * t320;
t195 = -mrSges(5,1) * t251 - mrSges(5,3) * t320;
t194 = mrSges(7,1) * t251 - mrSges(7,3) * t320;
t193 = mrSges(5,2) * t251 - mrSges(5,3) * t321;
t192 = -mrSges(7,2) * t251 + mrSges(7,3) * t321;
t191 = (Ifges(4,1) * t251 - t338) * qJD(3);
t187 = (-Ifges(4,2) * t248 + t337) * qJD(3);
t180 = (mrSges(4,1) * t248 + mrSges(4,2) * t251) * qJD(3);
t179 = t273 * qJD(4);
t177 = t272 * qJD(4);
t175 = pkin(3) - t354;
t172 = t273 * t248;
t171 = (-mrSges(7,1) * t247 + mrSges(7,2) * t250) * t248;
t170 = t272 * t248;
t164 = qJD(4) * t207 - qJD(6) * t247;
t163 = qJD(4) * t263 - t308;
t162 = -t311 * t339 - t305;
t156 = t261 * t248;
t151 = -Ifges(6,2) * t251 + (Ifges(6,4) * t250 + Ifges(6,6) * t247) * t248;
t149 = -Ifges(5,3) * t251 + (Ifges(5,5) * t250 - Ifges(5,6) * t247) * t248;
t147 = Ifges(7,3) * t251 + t248 * t265;
t137 = qJD(4) * t259 + t308;
t132 = -mrSges(6,2) * t256 + mrSges(6,3) * t313;
t131 = -mrSges(5,2) * t313 - mrSges(5,3) * t256;
t130 = mrSges(7,2) * t313 + mrSges(7,3) * t256;
t129 = t222 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t311) * t248;
t128 = mrSges(5,1) * t313 - mrSges(5,3) * t257;
t127 = t223 + (-mrSges(7,1) * t248 - t302) * qJD(3);
t123 = -t142 + t244;
t121 = -mrSges(4,1) * t322 - mrSges(4,3) * t166;
t120 = mrSges(4,2) * t322 - mrSges(4,3) * t165;
t119 = t258 * t248;
t108 = qJ(6) * t321 + t122;
t107 = mrSges(5,1) * t256 + mrSges(5,2) * t257;
t105 = mrSges(6,1) * t256 - mrSges(6,3) * t257;
t104 = pkin(5) * t251 + t233 + t244 + (-qJ(6) * t248 - t202) * t250;
t99 = -Ifges(6,4) * t294 + t290;
t97 = -Ifges(5,5) * t294 - Ifges(5,6) * t256 + t316;
t95 = -t208 * t310 + (-Ifges(7,3) * t248 + t251 * t265) * qJD(3);
t94 = mrSges(4,1) * t296 - mrSges(4,3) * t116;
t93 = -mrSges(4,2) * t296 - mrSges(4,3) * t115;
t92 = Ifges(4,1) * t166 - Ifges(4,4) * t165 - Ifges(4,5) * t322;
t91 = Ifges(4,4) * t166 - Ifges(4,2) * t165 - t303;
t88 = t313 * t345 - t283;
t86 = (qJD(4) * t264 - t307) * t248 + t261 * t312;
t85 = -mrSges(6,1) * t165 + mrSges(6,2) * t118;
t84 = mrSges(5,1) * t165 - mrSges(5,3) * t118;
t83 = -mrSges(7,1) * t165 - mrSges(7,3) * t118;
t82 = -mrSges(5,2) * t165 - mrSges(5,3) * t117;
t81 = mrSges(7,2) * t165 + mrSges(7,3) * t117;
t80 = -mrSges(6,2) * t117 + mrSges(6,3) * t165;
t79 = t297 * t313 + t283;
t73 = t87 + t356;
t68 = mrSges(5,1) * t117 + mrSges(5,2) * t118;
t67 = -mrSges(7,1) * t117 + mrSges(7,2) * t118;
t66 = mrSges(6,1) * t117 - mrSges(6,3) * t118;
t65 = mrSges(4,1) * t115 + mrSges(4,2) * t116;
t64 = (qJD(4) * t354 + t307) * t248 + t258 * t312;
t63 = Ifges(4,1) * t116 - Ifges(4,4) * t115 + Ifges(4,5) * t296;
t62 = Ifges(4,4) * t116 - Ifges(4,2) * t115 + Ifges(4,6) * t296;
t48 = Ifges(6,4) * t118 + Ifges(6,2) * t165 + Ifges(6,6) * t117;
t46 = Ifges(5,5) * t118 - Ifges(5,6) * t117 + Ifges(5,3) * t165;
t44 = Ifges(7,5) * t118 + Ifges(7,6) * t117 - Ifges(7,3) * t165;
t43 = (-pkin(9) * qJD(3) + qJ(6) * qJD(4)) * t320 + (qJD(6) * t248 + (-pkin(9) * qJD(4) + qJ(6) * qJD(3)) * t251) * t247 + t318 + t356;
t41 = (-qJ(6) * t312 - t198) * t250 + (qJ(6) * t311 - t305 + (-pkin(5) + t297) * qJD(3)) * t248 + t317;
t34 = -mrSges(6,2) * t60 + mrSges(6,3) * t115;
t32 = mrSges(5,1) * t115 - mrSges(5,3) * t61;
t30 = -mrSges(5,2) * t115 - mrSges(5,3) * t60;
t29 = mrSges(7,2) * t115 + mrSges(7,3) * t60;
t28 = pkin(4) * t117 - t260;
t25 = -pkin(4) * t165 - t26;
t23 = -t117 * t347 + t260;
t22 = mrSges(5,1) * t60 + mrSges(5,2) * t61;
t20 = mrSges(6,1) * t60 - mrSges(6,3) * t61;
t19 = qJ(6) * t117 + t24;
t18 = -qJ(6) * t118 - t165 * t347 - t26;
t8 = pkin(4) * t60 - t254;
t5 = -pkin(4) * t115 - t7;
t3 = -t347 * t60 + t254;
t2 = qJ(6) * t60 + qJD(6) * t117 + t4;
t1 = -qJ(6) * t61 - qJD(6) * t118 - t115 * t347 - t7;
t35 = [(-t255 * t77 + t26 * t7 + t27 * t6) * t353 - 0.2e1 * t255 * t68 + (-t62 + 0.2e1 * t327 + t355) * t165 + (t63 + 0.2e1 * t326) * t166 + (t220 - 0.2e1 * t328 - 0.2e1 * t329) * t246 + (-t252 * t300 + 0.2e1 * (t159 * t252 + t160 * t249) * mrSges(3,3) + ((t168 * t349 + Ifges(3,5) * t246 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t252) * t245) * t252 + (t169 * t349 + Ifges(4,5) * t166 - 0.2e1 * Ifges(3,6) * t246 - Ifges(4,6) * t165 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t249 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t252) * t245) * t249) * qJD(2)) * t245 + 0.2e1 * m(3) * (t159 * t169 - t160 * t168) + 0.2e1 * m(4) * (t144 * t160 + t39 * t90 + t40 * t89) + (t50 + t51 + t52) * t61 + (t45 + t47 - t49) * t60 + (t15 + t16 + t17) * t118 + (t10 + t12 - t14) * t117 + (t46 + t48 - t44 - t91) * t115 + (t1 * t18 + t19 * t2 + t23 * t3) * t351 + (t24 * t4 + t25 * t5 + t28 * t8) * t352 + 0.2e1 * t23 * t21 + 0.2e1 * t28 * t20 + 0.2e1 * t19 * t29 + 0.2e1 * t27 * t30 + 0.2e1 * t18 * t31 + 0.2e1 * t26 * t32 + 0.2e1 * t25 * t33 + 0.2e1 * t24 * t34 + 0.2e1 * t8 * t66 + 0.2e1 * t3 * t67 + 0.2e1 * t77 * t22 + 0.2e1 * t4 * t80 + 0.2e1 * t2 * t81 + 0.2e1 * t6 * t82 + 0.2e1 * t1 * t83 + 0.2e1 * t7 * t84 + 0.2e1 * t5 * t85 + 0.2e1 * t90 * t93 + 0.2e1 * t89 * t94 + t116 * t92 + 0.2e1 * t39 * t120 + 0.2e1 * t40 * t121 + 0.2e1 * t144 * t65; m(6) * (t122 * t4 + t123 * t5 + t156 * t8 + t24 * t73 + t25 * t79 + t28 * t86) + m(7) * (t1 * t104 + t108 * t2 + t119 * t3 + t18 * t41 + t19 * t43 + t23 * t64) + (-t40 * mrSges(4,3) + t326 + t63 / 0.2e1 + Ifges(4,5) * t289 + t287 * t250 + t288 * t247 + (-t247 * t285 + t250 * t286) * qJD(4) + (-t91 / 0.2e1 - t44 / 0.2e1 + t46 / 0.2e1 + t48 / 0.2e1 - t90 * mrSges(4,3) + t303 / 0.2e1) * qJD(3) + (-qJD(3) * t120 + t22 - t94 + m(4) * (-qJD(3) * t90 - t40) - m(5) * t255) * pkin(9)) * t248 - t255 * t172 + (-m(4) * t160 - t65) * pkin(2) + t220 + (Ifges(4,6) * t289 + t39 * mrSges(4,3) + t9 / 0.2e1 - t11 / 0.2e1 - t13 / 0.2e1 + t62 / 0.2e1 - t327 + (m(4) * t39 + t93) * pkin(9) + (-t89 * mrSges(4,3) + t92 / 0.2e1 + t285 * t250 + t286 * t247 + (-m(4) * t89 + m(5) * t77 - t121 + t68) * pkin(9)) * qJD(3)) * t251 + (-t252 * t242 / 0.2e1 - Ifges(3,6) * t314) * t245 - t329 - t328 + m(5) * (t142 * t7 + t143 * t6 + t26 * t88 + t27 * t87) + (t97 / 0.2e1 + t99 / 0.2e1 - t95 / 0.2e1 - t187 / 0.2e1) * t165 + t284 * t117 + t280 * t61 + t281 * t60 + t282 * t118 + (-t147 / 0.2e1 + t149 / 0.2e1 + t151 / 0.2e1 - t214 / 0.2e1) * t115 + t64 * t67 + t73 * t80 + t43 * t81 + t41 * t83 + t79 * t85 + t86 * t66 + t87 * t82 + t88 * t84 + t104 * t31 + t28 * t105 + t23 * t106 + t77 * t107 + t108 * t29 + t119 * t21 + t122 * t34 + t123 * t33 + t18 * t127 + t26 * t128 + t25 * t129 + t19 * t130 + t27 * t131 + t24 * t132 + t142 * t32 + t143 * t30 + t156 * t20 + t8 * t170 + t3 * t171 + t144 * t180 + t166 * t191 / 0.2e1 + t2 * t192 + t6 * t193 + t1 * t194 + t7 * t195 + t5 * t196 + t4 * t197 + t116 * t218 / 0.2e1; 0.2e1 * t119 * t106 + 0.2e1 * t104 * t127 + 0.2e1 * t123 * t129 + 0.2e1 * t108 * t130 + 0.2e1 * t122 * t132 + 0.2e1 * t142 * t128 + 0.2e1 * t143 * t131 + 0.2e1 * t156 * t105 + 0.2e1 * t86 * t170 + 0.2e1 * t64 * t171 - 0.2e1 * pkin(2) * t180 + 0.2e1 * t43 * t192 + 0.2e1 * t87 * t193 + 0.2e1 * t41 * t194 + 0.2e1 * t88 * t195 + 0.2e1 * t79 * t196 + 0.2e1 * t73 * t197 + (t142 * t88 + t143 * t87) * t353 + (t122 * t73 + t123 * t79 + t156 * t86) * t352 + (t104 * t41 + t108 * t43 + t119 * t64) * t351 + (t187 + t95 - t97 - t99 + (t172 * t350 + t247 * t299 + t250 * t298 + t218) * qJD(3)) * t251 + (t107 * t350 + t191 + (t101 + t102 + t103) * t250 + (t96 + t98 - t100) * t247 + ((pkin(9) ^ 2) * t251 * t353 - t147 + t149 + t151 - t214) * qJD(3) + (-t247 * t298 + t250 * t299) * qJD(4)) * t248; m(5) * (pkin(3) * t255 + t342 * t6 - t343 * t7) - t255 * t206 + m(6) * (t163 * t28 + t199 * t8 + t342 * t4 + t343 * t5) + t300 + (-t2 * mrSges(7,3) + t4 * mrSges(6,2) + t6 * mrSges(5,3) + (t30 + t34) * pkin(10) - t288) * t250 + ((t25 * mrSges(6,2) - t26 * mrSges(5,3) - t18 * mrSges(7,3) + t285) * t250 + (-t24 * mrSges(6,2) - t27 * mrSges(5,3) + t19 * mrSges(7,3) + t286) * t247 + ((-t84 + t85) * t250 + (-t80 - t82) * t247 + m(6) * (-t24 * t247 + t25 * t250) + m(5) * (-t247 * t27 - t250 * t26)) * pkin(10)) * qJD(4) + (t5 * mrSges(6,2) - t7 * mrSges(5,3) - t1 * mrSges(7,3) + (-t32 + t33) * pkin(10) + t287) * t247 + t274 * t61 + t275 * t60 + t276 * t115 + t277 * t118 + t278 * t117 - t279 * t165 + m(7) * (t1 * t203 + t137 * t23 + t162 * t19 + t164 * t18 + t175 * t3 + t2 * t207) - pkin(3) * t22 - t39 * mrSges(4,2) + t40 * mrSges(4,1) + t137 * t67 + t162 * t81 + t163 * t66 + t164 * t83 + t175 * t21 + t28 * t177 + t23 * t178 + t77 * t179 + t199 * t20 + t203 * t31 + t8 * t204 + t3 * t205 + t207 * t29; t248 * pkin(9) * t179 - pkin(3) * t107 + t199 * t105 + t175 * t106 + t119 * t178 + t203 * t127 + t207 * t130 + t137 * t171 + t156 * t177 + t162 * t192 + t163 * t170 + t164 * t194 + t86 * t204 + t64 * t205 + t242 + m(7) * (t104 * t164 + t108 * t162 + t119 * t137 + t175 * t64 + t203 * t41 + t207 * t43) + m(6) * (t156 * t163 + t199 * t86) + t279 * t251 + ((-m(5) * pkin(3) - mrSges(4,1) + t206) * t344 + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t276) * t248) * qJD(3) + (t87 * mrSges(5,3) + t73 * mrSges(6,2) - t43 * mrSges(7,3) + t277 * t248 + t274 * t312 + (t123 * mrSges(6,2) - t142 * mrSges(5,3) - t104 * mrSges(7,3) + t275 * t248 + t280) * qJD(4) + (t131 + t132 + (-t195 + t196) * qJD(4) + m(5) * (-qJD(4) * t142 + t87) + m(6) * (qJD(4) * t123 + t73)) * pkin(10) - t284) * t250 + (-t88 * mrSges(5,3) + t79 * mrSges(6,2) - t41 * mrSges(7,3) + t278 * t248 + t275 * t312 + (-t122 * mrSges(6,2) - t143 * mrSges(5,3) + t108 * mrSges(7,3) - t274 * t248 + t281) * qJD(4) + (-t128 + t129 + (-t193 - t197) * qJD(4) + m(5) * (-qJD(4) * t143 - t88) + m(6) * (-qJD(4) * t122 + t79)) * pkin(10) + t282) * t247; -0.2e1 * pkin(3) * t179 + 0.2e1 * t199 * t177 + 0.2e1 * t137 * t205 + 0.2e1 * t175 * t178 + (t137 * t175 + t162 * t207 + t164 * t203) * t351 + 0.2e1 * (m(6) * t199 + t204) * t163 + (t162 * t348 - t182 - t184 + t186) * t250 + (t164 * t348 + t188 + t189 + t190) * t247 + ((t203 * t348 + t215 + t216 + t217) * t250 + (0.2e1 * mrSges(7,3) * t207 + t209 + t211 - t213) * t247) * qJD(4); (t80 + t81) * qJD(5) + (t34 + t29) * qJ(5) + m(6) * (-pkin(4) * t5 + qJ(5) * t4 + qJD(5) * t24) + m(7) * (qJ(5) * t2 + qJD(5) * t19 - t1 * t347) - t1 * mrSges(7,1) + t2 * mrSges(7,2) + t4 * mrSges(6,3) - t5 * mrSges(6,1) - t6 * mrSges(5,2) + t7 * mrSges(5,1) - pkin(4) * t33 - t347 * t31 + t355; t290 + (t340 * t247 - t330) * t312 + (Ifges(7,3) * qJD(3) + (t340 * t250 + (Ifges(7,5) - t359) * t247) * qJD(4)) * t248 + (t130 + t132) * qJ(5) + (t192 + t197) * qJD(5) + m(7) * (qJ(5) * t43 + qJD(5) * t108 - t347 * t41) + m(6) * (-pkin(4) * t79 + qJ(5) * t73 + qJD(5) * t122) - t41 * mrSges(7,1) + t43 * mrSges(7,2) + t73 * mrSges(6,3) - t79 * mrSges(6,1) - t87 * mrSges(5,2) + t88 * mrSges(5,1) - pkin(4) * t129 - t347 * t127 + t316; t241 + m(7) * (qJ(5) * t162 + qJD(5) * t207 - t164 * t347) - t164 * mrSges(7,1) + t162 * mrSges(7,2) + qJD(5) * t262 + ((-mrSges(6,2) * pkin(4) + mrSges(7,3) * t347 - Ifges(7,5)) * t250 + (qJ(5) * t341 + t340) * t247 + (-m(6) * t264 + t204 + t206) * pkin(10)) * qJD(4) + t185; 0.2e1 * (mrSges(7,2) + mrSges(6,3) + (m(6) + m(7)) * qJ(5)) * qJD(5); m(6) * t5 + m(7) * t1 + t31 + t33; -mrSges(6,2) * t294 + t222 + t223 + m(6) * t79 + m(7) * t41 + (-t302 + (-mrSges(6,1) - mrSges(7,1)) * t248) * qJD(3); m(7) * t164 + qJD(4) * t262; 0; 0; m(7) * t3 + t21; m(7) * t64 + t106; m(7) * t137 + t178; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t35(1) t35(2) t35(4) t35(7) t35(11) t35(16); t35(2) t35(3) t35(5) t35(8) t35(12) t35(17); t35(4) t35(5) t35(6) t35(9) t35(13) t35(18); t35(7) t35(8) t35(9) t35(10) t35(14) t35(19); t35(11) t35(12) t35(13) t35(14) t35(15) t35(20); t35(16) t35(17) t35(18) t35(19) t35(20) t35(21);];
Mq  = res;
