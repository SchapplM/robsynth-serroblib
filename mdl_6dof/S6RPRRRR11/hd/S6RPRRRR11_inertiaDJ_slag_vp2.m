% Calculate time derivative of joint inertia matrix for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:25
% EndTime: 2019-03-09 07:38:45
% DurationCPUTime: 9.09s
% Computational Cost: add. (23549->729), mult. (67567->1072), div. (0->0), fcn. (72550->14), ass. (0->291)
t244 = sin(pkin(7));
t247 = cos(pkin(7));
t248 = cos(pkin(6));
t243 = sin(pkin(13));
t245 = sin(pkin(6));
t246 = cos(pkin(13));
t307 = t245 * t246;
t328 = pkin(1) * t248;
t299 = qJ(2) * t307 + t243 * t328;
t158 = (t244 * t248 + t247 * t307) * pkin(9) + t299;
t252 = sin(qJ(3));
t256 = cos(qJ(3));
t235 = t246 * t328;
t311 = t243 * t245;
t164 = pkin(2) * t248 + t235 + (-pkin(9) * t247 - qJ(2)) * t311;
t184 = (-pkin(9) * t243 * t244 - pkin(2) * t246 - pkin(1)) * t245;
t267 = t164 * t247 + t184 * t244;
t98 = -t252 * t158 + t267 * t256;
t250 = sin(qJ(5));
t254 = cos(qJ(5));
t251 = sin(qJ(4));
t292 = qJD(5) * t251;
t255 = cos(qJ(4));
t294 = qJD(4) * t255;
t259 = -t250 * t292 + t254 * t294;
t220 = -pkin(4) * t255 - pkin(11) * t251 - pkin(3);
t302 = t254 * t255;
t236 = pkin(10) * t302;
t187 = t250 * t220 + t236;
t361 = qJD(5) * t187;
t305 = t247 * t256;
t308 = t244 * t256;
t360 = (-t243 * t252 + t246 * t305) * t245 + t248 * t308;
t359 = -2 * Ifges(4,4);
t249 = sin(qJ(6));
t253 = cos(qJ(6));
t266 = t249 * t250 - t253 * t254;
t334 = -t266 / 0.2e1;
t205 = t249 * t254 + t250 * t253;
t333 = t205 / 0.2e1;
t358 = t250 / 0.2e1;
t329 = t254 / 0.2e1;
t193 = t266 * t251;
t221 = -mrSges(6,1) * t254 + mrSges(6,2) * t250;
t356 = -m(6) * pkin(4) - mrSges(5,1) + t221;
t355 = qJD(5) + qJD(6);
t354 = 0.2e1 * m(6);
t353 = 2 * m(7);
t352 = 0.2e1 * pkin(10);
t351 = -2 * mrSges(4,3);
t350 = m(5) * pkin(3);
t306 = t247 * t252;
t309 = t244 * t252;
t163 = t248 * t309 + (t243 * t256 + t246 * t306) * t245;
t194 = -t244 * t307 + t247 * t248;
t135 = t163 * t251 - t194 * t255;
t155 = t360 * qJD(3);
t105 = -t135 * qJD(4) + t155 * t255;
t136 = t163 * t255 + t194 * t251;
t107 = t136 * t254 - t250 * t360;
t156 = t163 * qJD(3);
t60 = -t107 * qJD(5) - t105 * t250 + t156 * t254;
t348 = t60 / 0.2e1;
t106 = -t136 * t250 - t254 * t360;
t71 = t106 * t253 - t107 * t249;
t347 = t71 / 0.2e1;
t72 = t106 * t249 + t107 * t253;
t346 = t72 / 0.2e1;
t345 = -pkin(12) - pkin(11);
t298 = qJD(2) * t245;
t281 = t243 * t298;
t272 = t244 * t281;
t111 = pkin(3) * t156 - pkin(10) * t155 + t272;
t295 = qJD(4) * t251;
t82 = (-t243 * t306 + t246 * t256) * t298 + t98 * qJD(3);
t127 = -t164 * t244 + t247 * t184;
t90 = -pkin(3) * t360 - pkin(10) * t163 + t127;
t151 = t256 * t158;
t99 = t164 * t306 + t184 * t309 + t151;
t96 = pkin(10) * t194 + t99;
t32 = t111 * t255 - t251 * t82 - t96 * t294 - t90 * t295;
t30 = -pkin(4) * t156 - t32;
t344 = m(6) * t30;
t323 = Ifges(6,4) * t250;
t224 = Ifges(6,2) * t254 + t323;
t322 = Ifges(6,4) * t254;
t268 = -Ifges(6,2) * t250 + t322;
t146 = -t224 * t292 + (Ifges(6,6) * t251 + t268 * t255) * qJD(4);
t343 = t146 / 0.2e1;
t226 = Ifges(6,1) * t250 + t322;
t269 = Ifges(6,1) * t254 - t323;
t147 = -t226 * t292 + (Ifges(6,5) * t251 + t269 * t255) * qJD(4);
t342 = t147 / 0.2e1;
t165 = t355 * t266;
t341 = -t165 / 0.2e1;
t166 = t355 * t205;
t340 = -t166 / 0.2e1;
t170 = Ifges(7,4) * t205 - Ifges(7,2) * t266;
t339 = t170 / 0.2e1;
t171 = Ifges(7,1) * t205 - Ifges(7,4) * t266;
t338 = t171 / 0.2e1;
t191 = -Ifges(6,5) * t255 + t269 * t251;
t337 = t191 / 0.2e1;
t192 = t205 * t251;
t336 = -t192 / 0.2e1;
t335 = -t193 / 0.2e1;
t332 = t226 / 0.2e1;
t331 = -t250 / 0.2e1;
t330 = -t254 / 0.2e1;
t327 = pkin(10) * t250;
t61 = t106 * qJD(5) + t105 * t254 + t156 * t250;
t33 = -mrSges(6,1) * t60 + mrSges(6,2) * t61;
t79 = mrSges(5,1) * t156 - mrSges(5,3) * t105;
t326 = t33 - t79;
t48 = t251 * t90 + t255 * t96;
t46 = -pkin(11) * t360 + t48;
t95 = -t194 * pkin(3) - t98;
t59 = t135 * pkin(4) - t136 * pkin(11) + t95;
t27 = t250 * t59 + t254 * t46;
t325 = Ifges(5,4) * t251;
t324 = Ifges(5,4) * t255;
t321 = Ifges(6,6) * t250;
t320 = t156 * Ifges(5,5);
t319 = t156 * Ifges(5,6);
t318 = t360 * Ifges(5,6);
t83 = (t243 * t305 + t246 * t252) * t298 + (t267 * t252 + t151) * qJD(3);
t317 = t256 * t83;
t316 = -mrSges(5,1) * t255 + mrSges(5,2) * t251 - mrSges(4,1);
t109 = -mrSges(5,1) * t360 - mrSges(5,3) * t136;
t73 = -mrSges(6,1) * t106 + mrSges(6,2) * t107;
t315 = -t109 + t73;
t195 = -t255 * t247 + t251 * t309;
t296 = qJD(3) * t256;
t279 = t244 * t296;
t172 = -t195 * qJD(4) + t255 * t279;
t313 = t172 * t255;
t196 = t247 * t251 + t255 * t309;
t173 = t196 * qJD(4) + t251 * t279;
t312 = t173 * t251;
t141 = t195 * t173;
t304 = t250 * t251;
t303 = t251 * t254;
t301 = Ifges(4,5) * t155 - Ifges(4,6) * t156;
t123 = -Ifges(7,5) * t165 - Ifges(7,6) * t166;
t218 = (pkin(4) * t251 - pkin(11) * t255) * qJD(4);
t300 = t254 * t218 + t295 * t327;
t297 = qJD(3) * t252;
t293 = qJD(5) * t250;
t291 = qJD(5) * t254;
t290 = qJD(6) * t249;
t289 = qJD(6) * t253;
t104 = t136 * qJD(4) + t155 * t251;
t20 = t71 * qJD(6) + t249 * t60 + t253 * t61;
t21 = -t72 * qJD(6) - t249 * t61 + t253 * t60;
t6 = Ifges(7,5) * t20 + Ifges(7,6) * t21 + Ifges(7,3) * t104;
t23 = Ifges(6,5) * t61 + Ifges(6,6) * t60 + Ifges(6,3) * t104;
t287 = pkin(5) * t293;
t53 = Ifges(6,1) * t107 + Ifges(6,4) * t106 + Ifges(6,5) * t135;
t284 = t53 * t329;
t283 = Ifges(5,5) * t105 - Ifges(5,6) * t104 + Ifges(5,3) * t156;
t130 = -t166 * t251 - t266 * t294;
t131 = t355 * t193 - t205 * t294;
t87 = Ifges(7,5) * t130 + Ifges(7,6) * t131 + Ifges(7,3) * t295;
t282 = qJD(5) * t345;
t280 = t244 * t297;
t241 = Ifges(6,5) * t291;
t276 = -Ifges(6,6) * t293 / 0.2e1 + t241 / 0.2e1 + t123 / 0.2e1;
t275 = Ifges(6,5) * t358 + Ifges(7,5) * t333 + Ifges(6,6) * t329 + Ifges(7,6) * t334;
t174 = -t250 * t196 - t254 * t308;
t113 = t174 * qJD(5) + t254 * t172 + t250 * t280;
t263 = -t254 * t196 + t250 * t308;
t114 = t263 * qJD(5) - t250 * t172 + t254 * t280;
t120 = t174 * t253 + t249 * t263;
t65 = t120 * qJD(6) + t113 * t253 + t114 * t249;
t121 = t174 * t249 - t253 * t263;
t66 = -t121 * qJD(6) - t113 * t249 + t114 * t253;
t274 = t66 * mrSges(7,1) - t65 * mrSges(7,2);
t26 = -t250 * t46 + t254 * t59;
t47 = -t251 * t96 + t255 * t90;
t139 = t250 * t218 + t220 * t291 + (-t254 * t295 - t255 * t293) * pkin(10);
t203 = t254 * t220;
t186 = -t255 * t327 + t203;
t273 = -qJD(5) * t186 + t139;
t31 = t251 * t111 + t255 * t82 + t90 * t294 - t96 * t295;
t29 = pkin(11) * t156 + t31;
t43 = t104 * pkin(4) - t105 * pkin(11) + t83;
t10 = -t27 * qJD(5) - t250 * t29 + t254 * t43;
t9 = t250 * t43 + t254 * t29 + t59 * t291 - t46 * t293;
t271 = -t10 * t250 + t9 * t254;
t270 = mrSges(6,1) * t250 + mrSges(6,2) * t254;
t17 = pkin(5) * t135 - pkin(12) * t107 + t26;
t22 = pkin(12) * t106 + t27;
t12 = t17 * t253 - t22 * t249;
t13 = t17 * t249 + t22 * t253;
t159 = -pkin(12) * t303 + t203 + (-pkin(5) - t327) * t255;
t176 = -pkin(12) * t304 + t187;
t118 = t159 * t253 - t176 * t249;
t119 = t159 * t249 + t176 * t253;
t229 = t345 * t250;
t230 = t345 * t254;
t178 = t229 * t253 + t230 * t249;
t179 = t229 * t249 - t230 * t253;
t216 = t250 * t282;
t217 = t254 * t282;
t133 = t178 * qJD(6) + t216 * t253 + t217 * t249;
t134 = -t179 * qJD(6) - t216 * t249 + t217 * t253;
t265 = t134 * mrSges(7,1) - t133 * mrSges(7,2) + t123;
t45 = pkin(4) * t360 - t47;
t4 = pkin(5) * t104 - pkin(12) * t61 + t10;
t5 = pkin(12) * t60 + t9;
t2 = t12 * qJD(6) + t249 * t4 + t253 * t5;
t3 = -t13 * qJD(6) - t249 * t5 + t253 * t4;
t264 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t6;
t115 = (pkin(5) * t251 - pkin(12) * t302) * qJD(4) + (-t236 + (pkin(12) * t251 - t220) * t250) * qJD(5) + t300;
t258 = t250 * t294 + t251 * t291;
t126 = -t258 * pkin(12) + t139;
t68 = t118 * qJD(6) + t115 * t249 + t126 * t253;
t69 = -t119 * qJD(6) + t115 * t253 - t126 * t249;
t262 = t69 * mrSges(7,1) - t68 * mrSges(7,2) + t87;
t261 = t195 * t294 + t312;
t145 = Ifges(6,5) * t259 - t258 * Ifges(6,6) + Ifges(6,3) * t295;
t242 = Ifges(5,5) * t294;
t238 = -pkin(5) * t254 - pkin(4);
t227 = Ifges(5,1) * t251 + t324;
t225 = Ifges(5,2) * t255 + t325;
t219 = (pkin(5) * t250 + pkin(10)) * t251;
t215 = -mrSges(6,1) * t255 - mrSges(6,3) * t303;
t214 = mrSges(6,2) * t255 - mrSges(6,3) * t304;
t212 = (Ifges(5,1) * t255 - t325) * qJD(4);
t211 = t269 * qJD(5);
t210 = (-Ifges(5,2) * t251 + t324) * qJD(4);
t209 = t268 * qJD(5);
t207 = (mrSges(5,1) * t251 + mrSges(5,2) * t255) * qJD(4);
t206 = t270 * qJD(5);
t201 = (-mrSges(7,1) * t249 - mrSges(7,2) * t253) * qJD(6) * pkin(5);
t199 = t270 * t251;
t190 = -Ifges(6,6) * t255 + t268 * t251;
t189 = -Ifges(6,3) * t255 + (Ifges(6,5) * t254 - t321) * t251;
t185 = t258 * pkin(5) + pkin(10) * t294;
t183 = -mrSges(6,2) * t295 - t258 * mrSges(6,3);
t182 = mrSges(6,1) * t295 - t259 * mrSges(6,3);
t181 = -mrSges(7,1) * t255 + mrSges(7,3) * t193;
t180 = mrSges(7,2) * t255 - mrSges(7,3) * t192;
t168 = mrSges(7,1) * t266 + mrSges(7,2) * t205;
t157 = t258 * mrSges(6,1) + t259 * mrSges(6,2);
t148 = mrSges(7,1) * t192 - mrSges(7,2) * t193;
t144 = -Ifges(7,1) * t193 - Ifges(7,4) * t192 - Ifges(7,5) * t255;
t143 = -Ifges(7,4) * t193 - Ifges(7,2) * t192 - Ifges(7,6) * t255;
t142 = -Ifges(7,5) * t193 - Ifges(7,6) * t192 - Ifges(7,3) * t255;
t140 = t300 - t361;
t138 = mrSges(4,1) * t194 - mrSges(4,3) * t163;
t137 = -mrSges(4,2) * t194 + mrSges(4,3) * t360;
t125 = -Ifges(7,1) * t165 - Ifges(7,4) * t166;
t124 = -Ifges(7,4) * t165 - Ifges(7,2) * t166;
t122 = mrSges(7,1) * t166 - mrSges(7,2) * t165;
t117 = -mrSges(7,2) * t295 + mrSges(7,3) * t131;
t116 = mrSges(7,1) * t295 - mrSges(7,3) * t130;
t112 = mrSges(4,1) * t156 + mrSges(4,2) * t155;
t108 = mrSges(5,2) * t360 - mrSges(5,3) * t135;
t97 = mrSges(5,1) * t135 + mrSges(5,2) * t136;
t94 = -mrSges(7,1) * t131 + mrSges(7,2) * t130;
t89 = Ifges(7,1) * t130 + Ifges(7,4) * t131 + Ifges(7,5) * t295;
t88 = Ifges(7,4) * t130 + Ifges(7,2) * t131 + Ifges(7,6) * t295;
t78 = -mrSges(5,2) * t156 - mrSges(5,3) * t104;
t77 = Ifges(5,1) * t136 - Ifges(5,4) * t135 - Ifges(5,5) * t360;
t76 = Ifges(5,4) * t136 - Ifges(5,2) * t135 - t318;
t75 = mrSges(6,1) * t135 - mrSges(6,3) * t107;
t74 = -mrSges(6,2) * t135 + mrSges(6,3) * t106;
t70 = mrSges(5,1) * t104 + mrSges(5,2) * t105;
t63 = Ifges(5,1) * t105 - Ifges(5,4) * t104 + t320;
t62 = Ifges(5,4) * t105 - Ifges(5,2) * t104 + t319;
t52 = Ifges(6,4) * t107 + Ifges(6,2) * t106 + Ifges(6,6) * t135;
t51 = Ifges(6,5) * t107 + Ifges(6,6) * t106 + Ifges(6,3) * t135;
t50 = mrSges(7,1) * t135 - mrSges(7,3) * t72;
t49 = -mrSges(7,2) * t135 + mrSges(7,3) * t71;
t40 = mrSges(6,1) * t104 - mrSges(6,3) * t61;
t39 = -mrSges(6,2) * t104 + mrSges(6,3) * t60;
t38 = -pkin(5) * t106 + t45;
t37 = -mrSges(7,1) * t71 + mrSges(7,2) * t72;
t36 = Ifges(7,1) * t72 + Ifges(7,4) * t71 + Ifges(7,5) * t135;
t35 = Ifges(7,4) * t72 + Ifges(7,2) * t71 + Ifges(7,6) * t135;
t34 = Ifges(7,5) * t72 + Ifges(7,6) * t71 + Ifges(7,3) * t135;
t25 = Ifges(6,1) * t61 + Ifges(6,4) * t60 + Ifges(6,5) * t104;
t24 = Ifges(6,4) * t61 + Ifges(6,2) * t60 + Ifges(6,6) * t104;
t16 = -pkin(5) * t60 + t30;
t15 = -mrSges(7,2) * t104 + mrSges(7,3) * t21;
t14 = mrSges(7,1) * t104 - mrSges(7,3) * t20;
t11 = -mrSges(7,1) * t21 + mrSges(7,2) * t20;
t8 = Ifges(7,1) * t20 + Ifges(7,4) * t21 + Ifges(7,5) * t104;
t7 = Ifges(7,4) * t20 + Ifges(7,2) * t21 + Ifges(7,6) * t104;
t1 = [0.2e1 * (t246 * (-mrSges(3,2) * t248 + mrSges(3,3) * t307) + m(3) * (t299 * t246 + (qJ(2) * t311 - t235) * t243)) * t298 + (t23 + t6 - t62) * t135 + (t12 * t3 + t13 * t2 + t16 * t38) * t353 + (t10 * t26 + t27 * t9 + t30 * t45) * t354 + (t163 * t359 + Ifges(5,5) * t136 - Ifges(4,6) * t194 - Ifges(5,6) * t135 + t99 * t351 - ((2 * Ifges(4,2)) + Ifges(5,3)) * t360) * t156 + (0.2e1 * Ifges(4,1) * t163 + Ifges(4,5) * t194 + t98 * t351 - t359 * t360) * t155 + 0.2e1 * (-mrSges(4,1) * t360 + mrSges(4,2) * t163) * t272 - t360 * t283 + 0.2e1 * m(5) * (t31 * t48 + t32 * t47 + t83 * t95) + t136 * t63 + 0.2e1 * t82 * t137 - 0.2e1 * t83 * t138 + 0.2e1 * t127 * t112 + 0.2e1 * t31 * t108 + 0.2e1 * t32 * t109 + t106 * t24 + t107 * t25 + t105 * t77 + 0.2e1 * t83 * t97 + 0.2e1 * t95 * t70 + 0.2e1 * t48 * t78 + 0.2e1 * t47 * t79 + 0.2e1 * t30 * t73 + 0.2e1 * t9 * t74 + 0.2e1 * t10 * t75 + t71 * t7 + t72 * t8 + t60 * t52 + t61 * t53 + 0.2e1 * t2 * t49 + 0.2e1 * t3 * t50 + 0.2e1 * t45 * t33 + 0.2e1 * t26 * t40 + 0.2e1 * t27 * t39 + 0.2e1 * t16 * t37 + 0.2e1 * t38 * t11 + t21 * t35 + t20 * t36 + 0.2e1 * t12 * t14 + 0.2e1 * t13 * t15 + (-t76 + t51 + t34) * t104 + 0.2e1 * m(4) * (t127 * t272 + t82 * t99 - t83 * t98) + t194 * t301 - 0.2e1 * (mrSges(3,1) * t248 - mrSges(3,3) * t311) * t281; t172 * t108 + t247 * t112 + t113 * t74 + t114 * t75 + t120 * t14 + t121 * t15 + t174 * t40 - t263 * t39 + t196 * t78 + t65 * t49 + t66 * t50 + (t11 + t326) * t195 + (t37 + t315) * t173 + m(5) * (t172 * t48 - t173 * t47 - t195 * t32 + t196 * t31) + m(6) * (t10 * t174 + t113 * t27 + t114 * t26 + t173 * t45 + t195 * t30 - t263 * t9) + m(7) * (t12 * t66 + t120 * t3 + t121 * t2 + t13 * t65 + t16 * t195 + t173 * t38) + (-t256 * t70 + (-t155 * t256 - t156 * t252) * mrSges(4,3) + (t256 * t137 + (-t138 + t97) * t252) * qJD(3) + m(4) * (t247 * t281 + t252 * t82 + t99 * t296 - t98 * t297 - t317) + m(5) * (t95 * t297 - t317)) * t244; 0.2e1 * m(7) * (t120 * t66 + t121 * t65 + t141) + 0.2e1 * m(6) * (-t113 * t263 + t114 * t174 + t141) + 0.2e1 * m(5) * (-t244 ^ 2 * t252 * t296 + t196 * t172 + t141); t301 + (t320 / 0.2e1 + t63 / 0.2e1 + t25 * t329 + t24 * t331 - t32 * mrSges(5,3) + (t52 * t330 + t53 * t331) * qJD(5) + (t318 / 0.2e1 - t48 * mrSges(5,3) - t76 / 0.2e1 + t34 / 0.2e1 + t51 / 0.2e1) * qJD(4) + (-qJD(4) * t108 + t344 + m(5) * (-qJD(4) * t48 - t32) + t326) * pkin(10)) * t251 + (-t350 + t316) * t83 + (-t210 / 0.2e1 + t145 / 0.2e1 + t87 / 0.2e1) * t135 + (-t225 / 0.2e1 + t189 / 0.2e1 + t142 / 0.2e1) * t104 + t107 * t342 + t106 * t343 + t89 * t346 + t88 * t347 + t190 * t348 + t8 * t335 + t7 * t336 + t61 * t337 + m(7) * (t118 * t3 + t119 * t2 + t12 * t69 + t13 * t68 + t16 * t219 + t185 * t38) - t360 * t242 / 0.2e1 + t105 * t227 / 0.2e1 + t136 * t212 / 0.2e1 + t9 * t214 + t10 * t215 + t219 * t11 + t30 * t199 + t95 * t207 + t2 * t180 + t3 * t181 + t26 * t182 + t27 * t183 + t185 * t37 + t186 * t40 + t187 * t39 + t45 * t157 + t21 * t143 / 0.2e1 + t20 * t144 / 0.2e1 + t16 * t148 + t139 * t74 + t140 * t75 + t131 * t35 / 0.2e1 + t130 * t36 / 0.2e1 + t119 * t15 + t12 * t116 + t13 * t117 + t118 * t14 + t38 * t94 - t82 * mrSges(4,2) - pkin(3) * t70 + t68 * t49 + t69 * t50 + m(6) * (t10 * t186 + t139 * t27 + t140 * t26 + t187 * t9) + (t319 / 0.2e1 + t62 / 0.2e1 - t23 / 0.2e1 - t6 / 0.2e1 + t31 * mrSges(5,3) + (t284 + t52 * t331 - t47 * mrSges(5,3) + t77 / 0.2e1) * qJD(4) + (m(5) * t31 + t78 + (-m(5) * t47 + m(6) * t45 + t315) * qJD(4)) * pkin(10)) * t255; t113 * t214 + t114 * t215 + t120 * t116 + t121 * t117 + t174 * t182 - t263 * t183 + t65 * t180 + t66 * t181 + (t157 + t94) * t195 + (t148 + t199) * t173 + (-t256 * t207 + (-mrSges(4,2) * t256 + t316 * t252) * qJD(3)) * t244 + m(7) * (t118 * t66 + t119 * t65 + t120 * t69 + t121 * t68 + t173 * t219 + t185 * t195) + m(6) * (t113 * t187 + t114 * t186 - t139 * t263 + t140 * t174) - t280 * t350 + (m(6) * t261 / 0.2e1 + m(5) * (-t196 * t295 + t261 + t313) / 0.2e1) * t352 + (t313 + t312 + (t195 * t255 - t196 * t251) * qJD(4)) * mrSges(5,3); -0.2e1 * pkin(3) * t207 + 0.2e1 * t118 * t116 + 0.2e1 * t119 * t117 + t130 * t144 + t131 * t143 + 0.2e1 * t139 * t214 + 0.2e1 * t140 * t215 + 0.2e1 * t185 * t148 + 0.2e1 * t68 * t180 + 0.2e1 * t69 * t181 + 0.2e1 * t186 * t182 + 0.2e1 * t187 * t183 - t192 * t88 - t193 * t89 + 0.2e1 * t219 * t94 + (t139 * t187 + t140 * t186) * t354 + (t118 * t69 + t119 * t68 + t185 * t219) * t353 + (-t145 + t210 - t87 + (-t190 * t250 + t191 * t254 + t199 * t352 + t227) * qJD(4)) * t255 + (t157 * t352 - t250 * t146 + t254 * t147 + t212 + (-t190 * t254 - t191 * t250) * qJD(5) + (pkin(10) ^ 2 * t255 * t354 + t142 + t189 - t225) * qJD(4)) * t251; (t12 * t165 - t13 * t166 - t2 * t266 - t205 * t3) * mrSges(7,3) + (-t33 - t344) * pkin(4) + t35 * t340 + t36 * t341 + t125 * t346 + t124 * t347 + t224 * t348 + t24 * t329 + t61 * t332 + t8 * t333 + t7 * t334 + t20 * t338 + t21 * t339 + t25 * t358 + t283 + t238 * t11 + t107 * t211 / 0.2e1 + t30 * t221 + t45 * t206 + t106 * t209 / 0.2e1 + t178 * t14 + t179 * t15 + t16 * t168 + t133 * t49 + t134 * t50 + t38 * t122 + t32 * mrSges(5,1) - t31 * mrSges(5,2) + ((-t250 * t27 - t254 * t26) * qJD(5) + t271) * mrSges(6,3) + t275 * t104 + t276 * t135 + (t284 + (-t52 / 0.2e1 + pkin(5) * t37) * t250) * qJD(5) + m(7) * (t12 * t134 + t13 * t133 + t16 * t238 + t178 * t3 + t179 * t2 + t38 * t287) + (-t75 * t291 + m(6) * (-t26 * t291 - t27 * t293 + t271) + t254 * t39 - t250 * t40 - t74 * t293) * pkin(11); -t172 * mrSges(5,2) + (t122 + t206) * t195 + m(7) * (t120 * t134 + t121 * t133 + t178 * t66 + t179 * t65 + t195 * t287) + (t120 * t165 - t121 * t166 - t205 * t66 - t266 * t65) * mrSges(7,3) + (m(6) * pkin(11) + mrSges(6,3)) * (t113 * t254 - t114 * t250 + (-t174 * t254 + t250 * t263) * qJD(5)) + (m(7) * t238 + t168 + t356) * t173; t238 * t94 + t219 * t122 + t88 * t334 + t89 * t333 + t124 * t336 + t125 * t335 + t179 * t117 + t133 * t180 + t134 * t181 + t185 * t168 + t131 * t339 + t130 * t338 + t178 * t116 + t144 * t341 + t143 * t340 - pkin(4) * t157 + m(7) * (t118 * t134 + t119 * t133 + t178 * t69 + t179 * t68 + t185 * t238) + t242 + (t356 * qJD(4) * pkin(10) - t276) * t255 + (t342 - t224 * t294 / 0.2e1 - t140 * mrSges(6,3) + (-t190 / 0.2e1 - t187 * mrSges(6,3) + (m(7) * t219 + t148) * pkin(5)) * qJD(5) + (m(6) * (-t140 - t361) - t182 - qJD(5) * t214) * pkin(11)) * t250 + (t118 * t165 - t119 * t166 - t205 * t69 - t266 * t68) * mrSges(7,3) + (t343 + t294 * t332 + qJD(5) * t337 + t273 * mrSges(6,3) + (m(6) * t273 - qJD(5) * t215 + t183) * pkin(11)) * t254 + (t211 * t329 + t209 * t331 + pkin(10) * t206 + (t224 * t330 + t226 * t331) * qJD(5) + (pkin(10) * mrSges(5,2) - Ifges(5,6) + t275) * qJD(4)) * t251; 0.2e1 * t168 * t287 + 0.2e1 * t238 * t122 - t165 * t171 + t205 * t125 - t166 * t170 - t266 * t124 + (t133 * t179 + t134 * t178 + t238 * t287) * t353 + t250 * t211 - t224 * t293 - 0.2e1 * pkin(4) * t206 + (qJD(5) * t226 + t209) * t254 + 0.2e1 * (-t133 * t266 - t134 * t205 + t165 * t178 - t166 * t179) * mrSges(7,3); t10 * mrSges(6,1) - t9 * mrSges(6,2) + (m(7) * (-t12 * t290 + t13 * t289 + t2 * t249 + t253 * t3) - t50 * t290 + t253 * t14 + t49 * t289 + t249 * t15) * pkin(5) + t264 + t23; t114 * mrSges(6,1) - t113 * mrSges(6,2) + m(7) * (t249 * t65 + t253 * t66 + (-t120 * t249 + t121 * t253) * qJD(6)) * pkin(5) + t274; t140 * mrSges(6,1) - t139 * mrSges(6,2) + (m(7) * (-t118 * t290 + t119 * t289 + t249 * t68 + t253 * t69) - t181 * t290 + t253 * t116 + t180 * t289 + t249 * t117) * pkin(5) + t145 + t262; t241 + (pkin(11) * t221 - t321) * qJD(5) + (m(7) * (t133 * t249 + t134 * t253 + (-t178 * t249 + t179 * t253) * qJD(6)) + (t253 * t165 - t249 * t166 + (t205 * t249 - t253 * t266) * qJD(6)) * mrSges(7,3)) * pkin(5) + t265; 0.2e1 * t201; t264; t274; t262; t265; t201; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
