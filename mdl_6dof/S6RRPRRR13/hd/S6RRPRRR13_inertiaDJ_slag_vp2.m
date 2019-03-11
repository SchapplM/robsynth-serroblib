% Calculate time derivative of joint inertia matrix for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:45:01
% EndTime: 2019-03-09 14:45:17
% DurationCPUTime: 6.89s
% Computational Cost: add. (9809->707), mult. (24575->1029), div. (0->0), fcn. (22861->10), ass. (0->283)
t231 = sin(qJ(6));
t232 = sin(qJ(5));
t235 = cos(qJ(6));
t236 = cos(qJ(5));
t247 = t231 * t232 - t235 * t236;
t341 = qJD(5) + qJD(6);
t344 = t341 * t247;
t237 = cos(qJ(4));
t277 = qJD(5) * t237;
t233 = sin(qJ(4));
t283 = qJD(4) * t233;
t242 = t232 * t283 - t236 * t277;
t308 = pkin(10) * t237;
t309 = pkin(4) * t233;
t197 = qJ(3) - t308 + t309;
t239 = -pkin(2) - pkin(9);
t291 = t233 * t239;
t216 = t236 * t291;
t157 = t232 * t197 + t216;
t343 = qJD(5) * t157;
t319 = -t247 / 0.2e1;
t185 = t231 * t236 + t232 * t235;
t318 = t185 / 0.2e1;
t315 = t232 / 0.2e1;
t313 = t236 / 0.2e1;
t230 = cos(pkin(6));
t229 = sin(pkin(6));
t238 = cos(qJ(2));
t294 = t229 * t238;
t272 = t233 * t294;
t173 = t230 * t237 - t272;
t234 = sin(qJ(2));
t295 = t229 * t234;
t140 = -t173 * t232 + t236 * t295;
t141 = t173 * t236 + t232 * t295;
t172 = t230 * t233 + t237 * t294;
t76 = t140 * t235 - t141 * t231;
t77 = t140 * t231 + t141 * t235;
t342 = Ifges(6,5) * t141 + Ifges(7,5) * t77 + Ifges(6,6) * t140 + Ifges(7,6) * t76 + (Ifges(6,3) + Ifges(7,3)) * t172;
t287 = t232 ^ 2 + t236 ^ 2;
t278 = qJD(5) * t236;
t279 = qJD(5) * t232;
t285 = qJD(2) * t238;
t266 = t229 * t285;
t286 = qJD(2) * t229;
t267 = t234 * t286;
t213 = pkin(2) * t267;
t284 = qJD(3) * t234;
t116 = t213 + (-t284 + (pkin(9) * t234 - qJ(3) * t238) * qJD(2)) * t229;
t217 = pkin(8) * t295;
t310 = pkin(1) * t238;
t269 = -pkin(2) - t310;
t123 = pkin(3) * t295 + t217 + (-pkin(9) + t269) * t230;
t297 = qJ(3) * t234;
t143 = (t238 * t239 - pkin(1) - t297) * t229;
t220 = t230 * t234 * pkin(1);
t330 = pkin(3) + pkin(8);
t144 = (t294 * t330 + t220) * qJD(2);
t281 = qJD(4) * t237;
t38 = t237 * t116 + t123 * t281 - t143 * t283 + t233 * t144;
t34 = pkin(10) * t266 + t38;
t274 = t230 * t310;
t214 = qJD(2) * t274;
t225 = t230 * qJD(3);
t122 = -t267 * t330 + t214 + t225;
t138 = -qJD(4) * t172 + t233 * t267;
t139 = -qJD(4) * t272 + t230 * t281 - t237 * t267;
t50 = pkin(4) * t139 - pkin(10) * t138 + t122;
t73 = t233 * t123 + t237 * t143;
t63 = pkin(10) * t295 + t73;
t175 = pkin(8) * t294 + t220;
t158 = -t230 * qJ(3) - t175;
t142 = pkin(3) * t294 - t158;
t86 = pkin(4) * t172 - pkin(10) * t173 + t142;
t10 = t232 * t50 + t236 * t34 + t86 * t278 - t279 * t63;
t37 = t232 * t86 + t236 * t63;
t11 = -qJD(5) * t37 - t232 * t34 + t236 * t50;
t250 = t10 * t236 - t11 * t232;
t36 = -t232 * t63 + t236 * t86;
t59 = -qJD(5) * t141 - t138 * t232 + t236 * t266;
t45 = -mrSges(6,2) * t139 + mrSges(6,3) * t59;
t60 = qJD(5) * t140 + t138 * t236 + t232 * t266;
t46 = mrSges(6,1) * t139 - mrSges(6,3) * t60;
t340 = m(6) * (-t278 * t36 - t279 * t37 + t250) - t232 * t46 + t236 * t45;
t339 = 0.2e1 * m(6);
t338 = 2 * m(7);
t337 = -0.2e1 * pkin(1);
t336 = 2 * mrSges(4,1);
t335 = -2 * mrSges(3,3);
t251 = -pkin(2) * t238 - t297;
t159 = (-pkin(1) + t251) * t229;
t334 = -0.2e1 * t159;
t333 = m(6) * pkin(4);
t332 = t76 / 0.2e1;
t331 = t77 / 0.2e1;
t329 = -pkin(11) - pkin(10);
t328 = -t344 / 0.2e1;
t131 = t341 * t185;
t327 = -t131 / 0.2e1;
t136 = Ifges(7,4) * t185 - Ifges(7,2) * t247;
t326 = t136 / 0.2e1;
t137 = Ifges(7,1) * t185 - Ifges(7,4) * t247;
t325 = t137 / 0.2e1;
t324 = t140 / 0.2e1;
t323 = t141 / 0.2e1;
t304 = Ifges(6,4) * t232;
t253 = Ifges(6,1) * t236 - t304;
t163 = Ifges(6,5) * t233 + t237 * t253;
t322 = t163 / 0.2e1;
t167 = t185 * t237;
t321 = -t167 / 0.2e1;
t169 = t247 * t237;
t320 = -t169 / 0.2e1;
t204 = Ifges(6,2) * t236 + t304;
t317 = t204 / 0.2e1;
t316 = -t232 / 0.2e1;
t314 = -t236 / 0.2e1;
t311 = m(6) * t233;
t307 = Ifges(3,4) + Ifges(4,6);
t306 = Ifges(5,4) * t233;
t305 = Ifges(5,4) * t237;
t303 = Ifges(6,4) * t236;
t302 = Ifges(6,6) * t232;
t164 = -pkin(8) * t267 + t214;
t301 = t164 * mrSges(3,2);
t104 = mrSges(5,1) * t266 - mrSges(5,3) * t138;
t29 = -mrSges(6,1) * t59 + mrSges(6,2) * t60;
t300 = t104 - t29;
t147 = mrSges(5,1) * t295 - mrSges(5,3) * t173;
t87 = -mrSges(6,1) * t140 + mrSges(6,2) * t141;
t299 = -t147 + t87;
t201 = -mrSges(6,1) * t236 + mrSges(6,2) * t232;
t298 = t201 - mrSges(5,1);
t165 = t175 * qJD(2);
t296 = t165 * t234;
t293 = t232 * t237;
t292 = t233 * t236;
t290 = t236 * t237;
t289 = t237 * t239;
t79 = -Ifges(7,5) * t344 - Ifges(7,6) * t131;
t288 = Ifges(3,5) * t266 + Ifges(4,5) * t267;
t282 = qJD(4) * t236;
t280 = qJD(4) * t239;
t276 = qJD(6) * t231;
t275 = qJD(6) * t235;
t19 = qJD(6) * t76 + t231 * t59 + t235 * t60;
t20 = -qJD(6) * t77 - t231 * t60 + t235 * t59;
t5 = Ifges(7,5) * t19 + Ifges(7,6) * t20 + Ifges(7,3) * t139;
t25 = Ifges(6,5) * t60 + Ifges(6,6) * t59 + Ifges(6,3) * t139;
t89 = -t131 * t237 + t247 * t283;
t91 = t185 * t283 + t237 * t344;
t41 = Ifges(7,5) * t89 + Ifges(7,6) * t91 + Ifges(7,3) * t281;
t273 = pkin(5) * t279;
t224 = Ifges(6,5) * t278;
t271 = -Ifges(6,6) * t279 / 0.2e1 + t224 / 0.2e1 + t79 / 0.2e1;
t270 = Ifges(5,5) * t138 - Ifges(5,6) * t139 + Ifges(5,3) * t266;
t268 = qJD(5) * t329;
t265 = t237 * t280;
t264 = t233 * t279;
t263 = t232 * t277;
t260 = Ifges(6,5) * t315 + Ifges(7,5) * t318 + Ifges(6,6) * t313 + Ifges(7,6) * t319;
t88 = -qJD(4) * t169 - t131 * t233;
t90 = -qJD(4) * t167 + t233 * t344;
t259 = t90 * mrSges(7,1) - t88 * mrSges(7,2);
t258 = (mrSges(4,2) - mrSges(3,1)) * t165;
t109 = mrSges(7,1) * t167 - mrSges(7,2) * t169;
t182 = (pkin(5) * t232 - t239) * t237;
t257 = m(7) * t182 + t109;
t256 = -t232 * t239 + pkin(5);
t180 = t236 * t197;
t156 = -t232 * t291 + t180;
t181 = qJD(3) + (pkin(4) * t237 + pkin(10) * t233) * qJD(4);
t95 = t232 * t181 + t197 * t278 + t236 * t265 - t239 * t264;
t255 = -qJD(5) * t156 + t95;
t72 = t123 * t237 - t233 * t143;
t254 = mrSges(6,1) * t232 + mrSges(6,2) * t236;
t206 = Ifges(6,1) * t232 + t303;
t252 = -Ifges(6,2) * t232 + t303;
t24 = pkin(5) * t172 - pkin(11) * t141 + t36;
t28 = pkin(11) * t140 + t37;
t12 = -t231 * t28 + t235 * t24;
t13 = t231 * t24 + t235 * t28;
t248 = t233 * t72 - t237 * t73;
t118 = -pkin(11) * t290 + t233 * t256 + t180;
t133 = -pkin(11) * t293 + t157;
t66 = t118 * t235 - t133 * t231;
t67 = t118 * t231 + t133 * t235;
t208 = t329 * t232;
t209 = t329 * t236;
t148 = t208 * t235 + t209 * t231;
t149 = t208 * t231 - t209 * t235;
t195 = t232 * t268;
t196 = t236 * t268;
t93 = qJD(6) * t148 + t195 * t235 + t196 * t231;
t94 = -qJD(6) * t149 - t195 * t231 + t196 * t235;
t246 = t94 * mrSges(7,1) - t93 * mrSges(7,2) + t79;
t39 = -t233 * t116 - t123 * t283 - t143 * t281 + t144 * t237;
t4 = pkin(5) * t139 - pkin(11) * t60 + t11;
t8 = pkin(11) * t59 + t10;
t2 = qJD(6) * t12 + t231 * t4 + t235 * t8;
t3 = -qJD(6) * t13 - t231 * t8 + t235 * t4;
t245 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t5;
t171 = t236 * t181;
t61 = t171 + (-t216 + (pkin(11) * t237 - t197) * t232) * qJD(5) + (pkin(11) * t292 + t237 * t256) * qJD(4);
t70 = pkin(11) * t242 + t95;
t21 = qJD(6) * t66 + t231 * t61 + t235 * t70;
t22 = -qJD(6) * t67 - t231 * t70 + t235 * t61;
t244 = t22 * mrSges(7,1) - t21 * mrSges(7,2) + t41;
t62 = -pkin(4) * t295 - t72;
t243 = t233 * t282 + t263;
t35 = -pkin(4) * t266 - t39;
t106 = -Ifges(6,5) * t243 + Ifges(6,6) * t242 + Ifges(6,3) * t281;
t221 = -pkin(5) * t236 - pkin(4);
t207 = Ifges(5,1) * t237 - t306;
t205 = -Ifges(5,2) * t233 + t305;
t202 = t233 * mrSges(5,1) + t237 * mrSges(5,2);
t194 = mrSges(6,1) * t233 - mrSges(6,3) * t290;
t193 = -mrSges(6,2) * t233 - mrSges(6,3) * t293;
t192 = (-Ifges(5,1) * t233 - t305) * qJD(4);
t191 = t253 * qJD(5);
t190 = (-Ifges(5,2) * t237 - t306) * qJD(4);
t189 = t252 * qJD(5);
t187 = (mrSges(5,1) * t237 - mrSges(5,2) * t233) * qJD(4);
t186 = t254 * qJD(5);
t183 = -mrSges(4,1) * t294 - mrSges(4,3) * t230;
t178 = (-mrSges(7,1) * t231 - mrSges(7,2) * t235) * qJD(6) * pkin(5);
t177 = t254 * t237;
t174 = -t217 + t274;
t168 = t247 * t233;
t166 = t185 * t233;
t162 = Ifges(6,6) * t233 + t237 * t252;
t161 = Ifges(6,3) * t233 + (Ifges(6,5) * t236 - t302) * t237;
t160 = t230 * t269 + t217;
t155 = -t164 - t225;
t154 = -mrSges(6,2) * t281 + mrSges(6,3) * t242;
t153 = mrSges(6,1) * t281 + mrSges(6,3) * t243;
t152 = -pkin(5) * t242 + t233 * t280;
t151 = mrSges(7,1) * t233 + mrSges(7,3) * t169;
t150 = -mrSges(7,2) * t233 - mrSges(7,3) * t167;
t146 = -mrSges(5,2) * t295 - mrSges(5,3) * t172;
t145 = t213 + (-qJ(3) * t285 - t284) * t229;
t134 = mrSges(7,1) * t247 + mrSges(7,2) * t185;
t117 = -mrSges(6,1) * t242 - mrSges(6,2) * t243;
t112 = mrSges(5,1) * t172 + mrSges(5,2) * t173;
t108 = -t206 * t277 + (Ifges(6,5) * t237 - t233 * t253) * qJD(4);
t107 = -t204 * t277 + (Ifges(6,6) * t237 - t233 * t252) * qJD(4);
t105 = -mrSges(5,2) * t266 - mrSges(5,3) * t139;
t103 = Ifges(5,1) * t173 - Ifges(5,4) * t172 + Ifges(5,5) * t295;
t102 = Ifges(5,4) * t173 - Ifges(5,2) * t172 + Ifges(5,6) * t295;
t101 = -Ifges(7,1) * t169 - Ifges(7,4) * t167 + Ifges(7,5) * t233;
t100 = -Ifges(7,4) * t169 - Ifges(7,2) * t167 + Ifges(7,6) * t233;
t99 = -Ifges(7,5) * t169 - Ifges(7,6) * t167 + Ifges(7,3) * t233;
t98 = mrSges(6,1) * t172 - mrSges(6,3) * t141;
t97 = -mrSges(6,2) * t172 + mrSges(6,3) * t140;
t96 = -t232 * t265 + t171 - t343;
t82 = mrSges(5,1) * t139 + mrSges(5,2) * t138;
t81 = -Ifges(7,1) * t344 - Ifges(7,4) * t131;
t80 = -Ifges(7,4) * t344 - Ifges(7,2) * t131;
t78 = mrSges(7,1) * t131 - mrSges(7,2) * t344;
t69 = -mrSges(7,2) * t281 + mrSges(7,3) * t91;
t68 = mrSges(7,1) * t281 - mrSges(7,3) * t89;
t65 = Ifges(5,1) * t138 - Ifges(5,4) * t139 + Ifges(5,5) * t266;
t64 = Ifges(5,4) * t138 - Ifges(5,2) * t139 + Ifges(5,6) * t266;
t55 = Ifges(6,1) * t141 + Ifges(6,4) * t140 + Ifges(6,5) * t172;
t54 = Ifges(6,4) * t141 + Ifges(6,2) * t140 + Ifges(6,6) * t172;
t52 = mrSges(7,1) * t172 - mrSges(7,3) * t77;
t51 = -mrSges(7,2) * t172 + mrSges(7,3) * t76;
t47 = -pkin(5) * t140 + t62;
t44 = -mrSges(7,1) * t91 + mrSges(7,2) * t89;
t43 = Ifges(7,1) * t89 + Ifges(7,4) * t91 + Ifges(7,5) * t281;
t42 = Ifges(7,4) * t89 + Ifges(7,2) * t91 + Ifges(7,6) * t281;
t40 = -mrSges(7,1) * t76 + mrSges(7,2) * t77;
t32 = Ifges(7,1) * t77 + Ifges(7,4) * t76 + Ifges(7,5) * t172;
t31 = Ifges(7,4) * t77 + Ifges(7,2) * t76 + Ifges(7,6) * t172;
t27 = Ifges(6,1) * t60 + Ifges(6,4) * t59 + Ifges(6,5) * t139;
t26 = Ifges(6,4) * t60 + Ifges(6,2) * t59 + Ifges(6,6) * t139;
t23 = -pkin(5) * t59 + t35;
t15 = -mrSges(7,2) * t139 + mrSges(7,3) * t20;
t14 = mrSges(7,1) * t139 - mrSges(7,3) * t19;
t9 = -mrSges(7,1) * t20 + mrSges(7,2) * t19;
t7 = Ifges(7,1) * t19 + Ifges(7,4) * t20 + Ifges(7,5) * t139;
t6 = Ifges(7,4) * t19 + Ifges(7,2) * t20 + Ifges(7,6) * t139;
t1 = [(t12 * t3 + t13 * t2 + t23 * t47) * t338 + (t10 * t37 + t11 * t36 + t35 * t62) * t339 + (t5 - t64 + t25) * t172 + 0.2e1 * m(3) * (t164 * t175 - t165 * t174) + 0.2e1 * m(4) * (t145 * t159 + t155 * t158 + t160 * t165) + 0.2e1 * m(5) * (t122 * t142 + t38 * t73 + t39 * t72) + 0.2e1 * t155 * t183 + t173 * t65 + 0.2e1 * t38 * t146 + 0.2e1 * t39 * t147 + t140 * t26 + t141 * t27 + 0.2e1 * t142 * t82 + t138 * t103 + 0.2e1 * t122 * t112 + 0.2e1 * t10 * t97 + 0.2e1 * t11 * t98 + 0.2e1 * t72 * t104 + 0.2e1 * t73 * t105 + 0.2e1 * t35 * t87 + t76 * t6 + t77 * t7 + t59 * t54 + t60 * t55 + 0.2e1 * t62 * t29 + 0.2e1 * t2 * t51 + 0.2e1 * t3 * t52 + 0.2e1 * t37 * t45 + 0.2e1 * t36 * t46 + 0.2e1 * t47 * t9 + 0.2e1 * t23 * t40 + t20 * t31 + t19 * t32 + 0.2e1 * t12 * t14 + 0.2e1 * t13 * t15 + (0.2e1 * t258 + t288 - 0.2e1 * t301) * t230 + (t234 * t270 + 0.2e1 * t145 * (mrSges(4,2) * t238 - mrSges(4,3) * t234) + t296 * t336 + 0.2e1 * (t164 * t238 + t296) * mrSges(3,3) + ((t158 * t336 + mrSges(4,2) * t334 + t175 * t335 + (Ifges(4,5) - (2 * Ifges(3,6))) * t230 + (mrSges(3,1) * t337 - 0.2e1 * t234 * t307) * t229) * t234 + (t160 * t336 + t174 * t335 + mrSges(4,3) * t334 + Ifges(5,5) * t173 - Ifges(5,6) * t172 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t230 + (mrSges(3,2) * t337 + 0.2e1 * t238 * t307) * t229 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3)) * t295) * t238) * qJD(2)) * t229 + (-t102 + t342) * t139; t7 * t320 + t6 * t321 + t60 * t322 + t108 * t323 + t107 * t324 + t43 * t331 + t42 * t332 + (t238 * (Ifges(5,5) * t237 - Ifges(5,6) * t233) / 0.2e1 - Ifges(4,4) * t238 - Ifges(3,6) * t234 + t251 * mrSges(4,1)) * t286 - t301 + (t239 * t105 - t38 * mrSges(5,3) - t64 / 0.2e1 + t25 / 0.2e1 + t5 / 0.2e1) * t233 + (-t183 + t112) * qJD(3) + t258 + t122 * t202 + t138 * t207 / 0.2e1 + t10 * t193 + t11 * t194 + t142 * t187 + t173 * t192 / 0.2e1 + t35 * t177 + t182 * t9 + t59 * t162 / 0.2e1 + t152 * t40 + t36 * t153 + t37 * t154 - t155 * mrSges(4,3) + t156 * t46 + t157 * t45 + t2 * t150 + t3 * t151 + t62 * t117 + t95 * t97 + t96 * t98 + t20 * t100 / 0.2e1 + t19 * t101 / 0.2e1 + t23 * t109 + t91 * t31 / 0.2e1 + qJ(3) * t82 + t89 * t32 / 0.2e1 + t66 * t14 + t67 * t15 + t12 * t68 + t13 * t69 + t21 * t51 + t22 * t52 + t47 * t44 + m(7) * (t12 * t22 + t13 * t21 + t152 * t47 + t182 * t23 + t2 * t67 + t3 * t66) + m(4) * (-pkin(2) * t165 - qJ(3) * t155 - qJD(3) * t158) + (-t205 / 0.2e1 + t161 / 0.2e1 + t99 / 0.2e1) * t139 + t288 + ((-Ifges(5,5) * t233 - Ifges(5,6) * t237) * t295 / 0.2e1 - t237 * t102 / 0.2e1 - t233 * t103 / 0.2e1 - t55 * t292 / 0.2e1 + t233 * t54 * t315 + t248 * mrSges(5,3) + (-m(5) * t248 + t237 * t146 + t233 * t299 + t311 * t62) * t239 + t342 * t237 / 0.2e1) * qJD(4) + m(6) * (t157 * t10 + t156 * t11 - t289 * t35 + t96 * t36 + t95 * t37) + m(5) * (qJ(3) * t122 + qJD(3) * t142 + t289 * t39 + t291 * t38) + (-t190 / 0.2e1 + t106 / 0.2e1 + t41 / 0.2e1) * t172 + (t27 * t313 + t26 * t316 - t39 * mrSges(5,3) + t65 / 0.2e1 + t300 * t239 + (t314 * t54 + t316 * t55) * qJD(5)) * t237; 0.2e1 * qJ(3) * t187 + t91 * t100 + t89 * t101 + 0.2e1 * t152 * t109 + 0.2e1 * t21 * t150 + 0.2e1 * t22 * t151 + 0.2e1 * t156 * t153 + 0.2e1 * t157 * t154 - t167 * t42 - t169 * t43 + 0.2e1 * t182 * t44 + 0.2e1 * t95 * t193 + 0.2e1 * t96 * t194 + 0.2e1 * t66 * t68 + 0.2e1 * t67 * t69 + (t152 * t182 + t21 * t67 + t22 * t66) * t338 + (t156 * t96 + t157 * t95) * t339 + 0.2e1 * (mrSges(4,3) + t202 + (m(4) + m(5)) * qJ(3)) * qJD(3) + (t106 - t190 + t41 + (t162 * t232 - t163 * t236 + 0.2e1 * t177 * t239 - t207) * qJD(4)) * t233 + (-t232 * t107 + t236 * t108 - 0.2e1 * t239 * t117 + t192 + (-t162 * t236 - t163 * t232) * qJD(5) + (-0.2e1 * t239 ^ 2 * t311 + t161 - t205 + t99) * qJD(4)) * t237; mrSges(4,1) * t266 - t166 * t14 - t168 * t15 + t88 * t51 + t90 * t52 + m(7) * (t12 * t90 + t13 * t88 - t166 * t3 - t168 * t2) + m(4) * t165 + (-t9 + (-t232 * t98 + t236 * t97 + t146) * qJD(4) - m(7) * t23 + m(6) * (-qJD(4) * t232 * t36 + t282 * t37 - t35) + m(5) * (qJD(4) * t73 + t39) + t300) * t237 + (t105 + (-t232 * t97 - t236 * t98) * qJD(5) + m(5) * t38 + (-m(5) * t72 + m(6) * t62 + m(7) * t47 + t299 + t40) * qJD(4) + t340) * t233; m(7) * (-t166 * t22 - t168 * t21 + t66 * t90 + t67 * t88) + t88 * t150 - t168 * t69 + t90 * t151 - t166 * t68 + (-t232 * t194 + m(6) * (-t156 * t232 + t157 * t236) + t236 * t193) * t281 + (-t194 * t278 - t232 * t153 + m(6) * (-t156 * t278 - t157 * t279 - t232 * t96 + t236 * t95) - t193 * t279 + t236 * t154 + (-0.2e1 * m(6) * t289 + t177 + t257) * qJD(4)) * t233 + (-m(7) * t152 - t117 - t44) * t237; (-t166 * t90 - t168 * t88) * t338 + 0.4e1 * (m(6) * (-0.1e1 + t287) / 0.2e1 - m(7) / 0.2e1) * t233 * t281; t27 * t315 + t59 * t317 + t7 * t318 + t6 * t319 + t191 * t323 + t189 * t324 + t19 * t325 + t20 * t326 + t31 * t327 + t32 * t328 + t81 * t331 + t80 * t332 + t26 * t313 + (-t278 * t98 - t279 * t97 + t340) * pkin(10) + (t201 - t333) * t35 + ((-t232 * t37 - t236 * t36) * qJD(5) + t250) * mrSges(6,3) + t221 * t9 + t60 * t206 / 0.2e1 + t62 * t186 + t148 * t14 + t149 * t15 + t23 * t134 + t93 * t51 + t94 * t52 + t47 * t78 - t38 * mrSges(5,2) + t39 * mrSges(5,1) - pkin(4) * t29 + t270 + t260 * t139 + t271 * t172 + m(7) * (t12 * t94 + t13 * t93 + t148 * t3 + t149 * t2 + t221 * t23 + t273 * t47) + (t12 * t344 - t13 * t131 - t185 * t3 - t2 * t247) * mrSges(7,3) + (t55 * t313 + (-t54 / 0.2e1 + pkin(5) * t40) * t232) * qJD(5); m(7) * (t148 * t22 + t149 * t21 + t152 * t221 + t66 * t94 + t67 * t93) + t221 * t44 + t42 * t319 + t43 * t318 + t182 * t78 + t80 * t321 + t81 * t320 + t152 * t134 + t148 * t68 + t149 * t69 + t93 * t150 + t94 * t151 + t100 * t327 + t91 * t326 + t89 * t325 + t101 * t328 - pkin(4) * t117 + ((-Ifges(5,5) + (t298 - t333) * t239) * qJD(4) + t271) * t233 + (-t96 * mrSges(6,3) + t283 * t317 + t108 / 0.2e1 + (-t162 / 0.2e1 - t157 * mrSges(6,3) + t257 * pkin(5)) * qJD(5) + (-t153 + m(6) * (-t96 - t343) - qJD(5) * t193) * pkin(10)) * t232 + (-t131 * t67 - t185 * t22 - t21 * t247 + t344 * t66) * mrSges(7,3) + (qJD(5) * t322 - t206 * t283 / 0.2e1 + t107 / 0.2e1 + t255 * mrSges(6,3) + (m(6) * t255 - qJD(5) * t194 + t154) * pkin(10)) * t236 + (-t239 * t186 + t191 * t313 + t189 * t316 + (t204 * t314 + t206 * t316) * qJD(5) + (-t239 * mrSges(5,2) - Ifges(5,6) + t260) * qJD(4)) * t237; m(7) * (-pkin(5) * t263 + t148 * t90 + t149 * t88 - t166 * t94 - t168 * t93) + (t131 * t168 - t166 * t344 - t185 * t90 - t247 * t88) * mrSges(7,3) + (m(6) * (t287 * t308 - t309) + (m(7) * t221 + t134 + t298) * t233) * qJD(4) + (-t186 - t78 + (mrSges(6,3) * t287 - mrSges(5,2)) * qJD(4)) * t237; -t131 * t136 - t247 * t80 - t344 * t137 + t185 * t81 + (t148 * t94 + t149 * t93 + t221 * t273) * t338 + 0.2e1 * t134 * t273 + 0.2e1 * t221 * t78 - 0.2e1 * pkin(4) * t186 + t232 * t191 - t204 * t279 + (qJD(5) * t206 + t189) * t236 + 0.2e1 * (-t131 * t149 + t148 * t344 - t185 * t94 - t247 * t93) * mrSges(7,3); t11 * mrSges(6,1) - t10 * mrSges(6,2) + (t51 * t275 + t231 * t15 - t52 * t276 + t235 * t14 + m(7) * (-t12 * t276 + t13 * t275 + t2 * t231 + t235 * t3)) * pkin(5) + t245 + t25; t96 * mrSges(6,1) - t95 * mrSges(6,2) + (m(7) * (t21 * t231 + t22 * t235 + t275 * t67 - t276 * t66) + t150 * t275 + t231 * t69 - t151 * t276 + t235 * t68) * pkin(5) + t106 + t244; (-t236 * t281 + t264) * mrSges(6,2) + (-t232 * t281 - t233 * t278) * mrSges(6,1) + m(7) * (t231 * t88 + t235 * t90 + (t166 * t231 - t168 * t235) * qJD(6)) * pkin(5) + t259; t224 + (pkin(10) * t201 - t302) * qJD(5) + (m(7) * (t231 * t93 + t235 * t94 + (-t148 * t231 + t149 * t235) * qJD(6)) + (t235 * t344 - t231 * t131 + (t185 * t231 - t235 * t247) * qJD(6)) * mrSges(7,3)) * pkin(5) + t246; 0.2e1 * t178; t245; t244; t259; t246; t178; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
