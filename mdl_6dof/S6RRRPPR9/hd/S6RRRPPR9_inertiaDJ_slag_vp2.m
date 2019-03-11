% Calculate time derivative of joint inertia matrix for
% S6RRRPPR9
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:41
% EndTime: 2019-03-09 16:11:56
% DurationCPUTime: 6.23s
% Computational Cost: add. (7051->699), mult. (19055->1001), div. (0->0), fcn. (17699->10), ass. (0->279)
t248 = sin(pkin(11));
t250 = cos(pkin(11));
t346 = qJD(4) * (t248 ^ 2 + t250 ^ 2);
t345 = 2 * pkin(9);
t344 = t248 / 0.2e1;
t343 = -t250 / 0.2e1;
t249 = sin(pkin(6));
t254 = sin(qJ(2));
t293 = qJD(2) * t254;
t281 = t249 * t293;
t301 = t249 * t254;
t235 = pkin(8) * t301;
t251 = cos(pkin(6));
t257 = cos(qJ(2));
t320 = pkin(1) * t257;
t198 = t251 * t320 - t235;
t181 = t198 * qJD(2);
t253 = sin(qJ(3));
t256 = cos(qJ(3));
t294 = qJD(2) * t249;
t258 = (pkin(2) * t254 - pkin(9) * t257) * t294;
t300 = t249 * t257;
t199 = t251 * t254 * pkin(1) + pkin(8) * t300;
t173 = pkin(9) * t251 + t199;
t175 = (-pkin(2) * t257 - pkin(9) * t254 - pkin(1)) * t249;
t340 = t256 * t173 + t253 * t175;
t55 = -qJD(3) * t340 - t181 * t253 + t256 * t258;
t47 = -pkin(3) * t281 - t55;
t189 = -t251 * t256 + t253 * t301;
t280 = t257 * t294;
t139 = -qJD(3) * t189 + t256 * t280;
t104 = t139 * t248 - t250 * t281;
t105 = t139 * t250 + t248 * t281;
t53 = t104 * mrSges(5,1) + t105 * mrSges(5,2);
t342 = -m(5) * t47 - t53;
t331 = pkin(4) + pkin(5);
t341 = t248 * t331;
t292 = qJD(3) * t253;
t339 = qJ(5) * t292 - qJD(5) * t256;
t273 = qJ(5) * t248 + pkin(3);
t210 = -pkin(4) * t250 - t273;
t214 = -mrSges(6,1) * t250 - mrSges(6,3) * t248;
t338 = m(6) * t210 + t214;
t337 = 0.2e1 * m(5);
t336 = 0.2e1 * m(6);
t335 = 2 * m(7);
t334 = -2 * mrSges(3,3);
t190 = t251 * t253 + t256 * t301;
t136 = t190 * t248 + t250 * t300;
t137 = t190 * t250 - t248 * t300;
t252 = sin(qJ(6));
t255 = cos(qJ(6));
t73 = t136 * t255 - t137 * t252;
t333 = t73 / 0.2e1;
t74 = t136 * t252 + t137 * t255;
t332 = t74 / 0.2e1;
t202 = t248 * t255 - t250 * t252;
t260 = t248 * t252 + t250 * t255;
t129 = Ifges(7,4) * t202 - Ifges(7,2) * t260;
t329 = t129 / 0.2e1;
t130 = Ifges(7,1) * t202 - Ifges(7,4) * t260;
t328 = t130 / 0.2e1;
t176 = t202 * t253;
t327 = t176 / 0.2e1;
t177 = t260 * t253;
t326 = t177 / 0.2e1;
t187 = t260 * qJD(6);
t325 = -t187 / 0.2e1;
t188 = t202 * qJD(6);
t324 = -t188 / 0.2e1;
t323 = -t260 / 0.2e1;
t322 = t202 / 0.2e1;
t319 = pkin(4) * t248;
t318 = pkin(9) * t253;
t317 = -pkin(10) + qJ(4);
t291 = qJD(3) * t256;
t54 = -t173 * t292 + t175 * t291 + t256 * t181 + t253 * t258;
t45 = (qJ(4) * t293 - qJD(4) * t257) * t249 + t54;
t138 = qJD(3) * t190 + t253 * t280;
t182 = t199 * qJD(2);
t48 = t138 * pkin(3) - t139 * qJ(4) - t190 * qJD(4) + t182;
t16 = t248 * t48 + t250 * t45;
t172 = t235 + (-pkin(2) - t320) * t251;
t81 = t189 * pkin(3) - t190 * qJ(4) + t172;
t82 = -qJ(4) * t300 + t340;
t41 = t248 * t81 + t250 * t82;
t316 = m(7) * qJD(5);
t315 = mrSges(6,3) * t250;
t314 = Ifges(4,4) * t253;
t313 = Ifges(4,4) * t256;
t312 = Ifges(5,4) * t248;
t311 = Ifges(5,4) * t250;
t310 = Ifges(6,5) * t248;
t309 = Ifges(6,5) * t250;
t308 = t181 * mrSges(3,2);
t307 = t182 * mrSges(3,1);
t306 = t182 * mrSges(4,1);
t305 = t182 * mrSges(4,2);
t186 = -qJD(4) * t253 + (pkin(3) * t253 - qJ(4) * t256) * qJD(3);
t304 = t186 * t250;
t303 = t248 * t253;
t302 = t248 * t256;
t299 = t250 * t253;
t298 = t250 * t256;
t92 = -t253 * t173 + t256 * t175;
t119 = -Ifges(7,5) * t187 - Ifges(7,6) * t188;
t278 = t250 * t291;
t297 = -qJ(5) * t278 - qJD(5) * t299;
t279 = t248 * t291;
t179 = mrSges(5,1) * t279 + mrSges(5,2) * t278;
t212 = -pkin(3) * t256 - qJ(4) * t253 - pkin(2);
t161 = pkin(9) * t298 + t248 * t212;
t296 = qJ(4) * t346;
t290 = qJD(5) * t248;
t25 = -qJD(6) * t74 + t104 * t255 - t105 * t252;
t26 = qJD(6) * t73 + t104 * t252 + t105 * t255;
t3 = Ifges(7,5) * t26 + Ifges(7,6) * t25 - Ifges(7,3) * t138;
t32 = t189 * qJ(5) + t41;
t287 = pkin(9) * t292;
t286 = Ifges(4,6) * t300;
t34 = Ifges(6,5) * t105 + Ifges(6,6) * t138 + Ifges(6,3) * t104;
t37 = Ifges(5,4) * t105 - Ifges(5,2) * t104 + Ifges(5,6) * t138;
t285 = -t37 / 0.2e1 + t34 / 0.2e1;
t38 = Ifges(6,1) * t105 + Ifges(6,4) * t138 + Ifges(6,5) * t104;
t39 = Ifges(5,1) * t105 - Ifges(5,4) * t104 + Ifges(5,5) * t138;
t284 = t38 / 0.2e1 + t39 / 0.2e1;
t283 = Ifges(4,5) * t139 - Ifges(4,6) * t138 + Ifges(4,3) * t281;
t83 = pkin(3) * t300 - t92;
t282 = -pkin(9) * t248 - pkin(4);
t262 = Ifges(6,3) * t248 + t309;
t150 = (Ifges(6,6) * t253 + t256 * t262) * qJD(3);
t265 = -Ifges(5,2) * t248 + t311;
t153 = (Ifges(5,6) * t253 + t256 * t265) * qJD(3);
t277 = t150 / 0.2e1 - t153 / 0.2e1;
t266 = Ifges(6,1) * t250 + t310;
t154 = (Ifges(6,4) * t253 + t256 * t266) * qJD(3);
t267 = Ifges(5,1) * t250 - t312;
t155 = (Ifges(5,5) * t253 + t256 * t267) * qJD(3);
t276 = t154 / 0.2e1 + t155 / 0.2e1;
t275 = t310 / 0.2e1 - t312 / 0.2e1 + (Ifges(6,3) + Ifges(5,2)) * t343;
t274 = -t309 / 0.2e1 + t311 / 0.2e1 + (Ifges(6,1) + Ifges(5,1)) * t344;
t70 = -t138 * mrSges(6,1) + t105 * mrSges(6,2);
t52 = t104 * mrSges(6,1) - t105 * mrSges(6,3);
t272 = qJ(5) * t250 - pkin(9);
t15 = -t248 * t45 + t250 * t48;
t40 = -t248 * t82 + t250 * t81;
t234 = pkin(9) * t302;
t160 = t212 * t250 - t234;
t12 = t138 * qJ(5) + t189 * qJD(5) + t16;
t269 = t281 / 0.2e1;
t196 = -mrSges(6,1) * t292 + mrSges(6,2) * t278;
t8 = -t25 * mrSges(7,1) + t26 * mrSges(7,2);
t143 = -qJ(5) * t256 + t161;
t268 = -Ifges(7,5) * t202 / 0.2e1 + Ifges(7,6) * t260 / 0.2e1 + Ifges(5,6) * t250 / 0.2e1 + Ifges(6,6) * t343 + (Ifges(5,5) + Ifges(6,4)) * t344;
t110 = qJD(6) * t176 + t260 * t291;
t111 = -t187 * t253 + t202 * t291;
t56 = -t111 * mrSges(7,1) + t110 * mrSges(7,2);
t118 = mrSges(7,1) * t188 - mrSges(7,2) * t187;
t264 = Ifges(6,4) * t250 + Ifges(6,6) * t248;
t263 = Ifges(5,5) * t250 - Ifges(5,6) * t248;
t21 = -pkin(10) * t137 - t189 * t331 - t40;
t22 = pkin(10) * t136 + t32;
t6 = t21 * t255 - t22 * t252;
t7 = t21 * t252 + t22 * t255;
t245 = t256 * pkin(4);
t117 = pkin(5) * t256 + t234 + t245 + (-pkin(10) * t253 - t212) * t250;
t123 = pkin(10) * t303 + t143;
t65 = t117 * t255 - t123 * t252;
t66 = t117 * t252 + t123 * t255;
t213 = t317 * t248;
t216 = t317 * t250;
t140 = t213 * t255 - t216 * t252;
t141 = t213 * t252 + t216 * t255;
t165 = t248 * t186;
t135 = -t250 * t287 + t165;
t49 = Ifges(7,5) * t110 + Ifges(7,6) * t111 - Ifges(7,3) * t292;
t259 = qJ(5) * t137 - t83;
t14 = t104 * pkin(4) - t105 * qJ(5) - t137 * qJD(5) + t47;
t244 = Ifges(4,5) * t291;
t231 = Ifges(3,5) * t280;
t226 = mrSges(6,1) * t279;
t224 = Ifges(4,1) * t253 + t313;
t223 = Ifges(4,2) * t256 + t314;
t215 = -mrSges(5,1) * t250 + mrSges(5,2) * t248;
t209 = (Ifges(4,1) * t256 - t314) * qJD(3);
t208 = (-Ifges(4,2) * t253 + t313) * qJD(3);
t207 = (mrSges(4,1) * t253 + mrSges(4,2) * t256) * qJD(3);
t206 = -mrSges(6,2) * t303 - mrSges(6,3) * t256;
t205 = mrSges(6,1) * t256 + mrSges(6,2) * t299;
t204 = -mrSges(5,1) * t256 - mrSges(5,3) * t299;
t203 = mrSges(5,2) * t256 - mrSges(5,3) * t303;
t197 = (-mrSges(6,2) * t302 + mrSges(6,3) * t253) * qJD(3);
t195 = (mrSges(5,1) * t253 - mrSges(5,3) * t298) * qJD(3);
t194 = (-mrSges(5,2) * t253 - mrSges(5,3) * t302) * qJD(3);
t193 = (mrSges(5,1) * t248 + mrSges(5,2) * t250) * t253;
t192 = (mrSges(6,1) * t248 - t315) * t253;
t191 = t250 * t331 + t273;
t178 = -mrSges(6,3) * t278 + t226;
t174 = (-t272 + t319) * t253;
t171 = -Ifges(5,5) * t256 + t253 * t267;
t170 = -Ifges(6,4) * t256 + t253 * t266;
t169 = -Ifges(5,6) * t256 + t253 * t265;
t168 = -Ifges(6,2) * t256 + t253 * t264;
t167 = -Ifges(5,3) * t256 + t253 * t263;
t166 = -Ifges(6,6) * t256 + t253 * t262;
t152 = (Ifges(6,2) * t253 + t256 * t264) * qJD(3);
t151 = (Ifges(5,3) * t253 + t256 * t263) * qJD(3);
t148 = mrSges(7,1) * t256 - mrSges(7,3) * t177;
t147 = -mrSges(7,2) * t256 + mrSges(7,3) * t176;
t146 = -mrSges(4,1) * t300 - t190 * mrSges(4,3);
t145 = mrSges(4,2) * t300 - t189 * mrSges(4,3);
t144 = -t160 + t245;
t142 = (t272 - t341) * t253;
t134 = t248 * t287 + t304;
t127 = mrSges(7,1) * t260 + mrSges(7,2) * t202;
t125 = (pkin(9) + t319) * t291 + t297;
t122 = t282 * t292 - t304;
t121 = -Ifges(7,1) * t187 - Ifges(7,4) * t188;
t120 = -Ifges(7,4) * t187 - Ifges(7,2) * t188;
t116 = (pkin(9) + t341) * t291 + t297;
t115 = mrSges(4,1) * t281 - mrSges(4,3) * t139;
t114 = -mrSges(4,2) * t281 - mrSges(4,3) * t138;
t113 = t135 + t339;
t112 = -mrSges(7,1) * t176 + mrSges(7,2) * t177;
t109 = Ifges(4,1) * t190 - Ifges(4,4) * t189 - Ifges(4,5) * t300;
t108 = Ifges(4,4) * t190 - Ifges(4,2) * t189 - t286;
t102 = Ifges(7,1) * t177 + Ifges(7,4) * t176 + Ifges(7,5) * t256;
t101 = Ifges(7,4) * t177 + Ifges(7,2) * t176 + Ifges(7,6) * t256;
t100 = Ifges(7,5) * t177 + Ifges(7,6) * t176 + Ifges(7,3) * t256;
t96 = qJD(4) * t202 - qJD(6) * t141;
t95 = qJD(4) * t260 + qJD(6) * t140;
t91 = t165 + (-pkin(9) * t299 + pkin(10) * t302) * qJD(3) + t339;
t90 = -mrSges(6,1) * t189 + mrSges(6,2) * t137;
t89 = mrSges(5,1) * t189 - mrSges(5,3) * t137;
t88 = -mrSges(5,2) * t189 - mrSges(5,3) * t136;
t87 = -mrSges(6,2) * t136 + mrSges(6,3) * t189;
t86 = -t304 + (-pkin(10) * t298 + (-pkin(5) + t282) * t253) * qJD(3);
t85 = mrSges(7,2) * t292 + mrSges(7,3) * t111;
t84 = -mrSges(7,1) * t292 - mrSges(7,3) * t110;
t77 = mrSges(4,1) * t138 + mrSges(4,2) * t139;
t76 = mrSges(5,1) * t136 + mrSges(5,2) * t137;
t75 = mrSges(6,1) * t136 - mrSges(6,3) * t137;
t72 = Ifges(4,1) * t139 - Ifges(4,4) * t138 + Ifges(4,5) * t281;
t71 = Ifges(4,4) * t139 - Ifges(4,2) * t138 + Ifges(4,6) * t281;
t69 = mrSges(5,1) * t138 - mrSges(5,3) * t105;
t68 = -mrSges(5,2) * t138 - mrSges(5,3) * t104;
t67 = -mrSges(6,2) * t104 + mrSges(6,3) * t138;
t64 = Ifges(5,1) * t137 - Ifges(5,4) * t136 + Ifges(5,5) * t189;
t63 = Ifges(6,1) * t137 + Ifges(6,4) * t189 + Ifges(6,5) * t136;
t62 = Ifges(5,4) * t137 - Ifges(5,2) * t136 + Ifges(5,6) * t189;
t61 = Ifges(6,4) * t137 + Ifges(6,2) * t189 + Ifges(6,6) * t136;
t60 = Ifges(5,5) * t137 - Ifges(5,6) * t136 + Ifges(5,3) * t189;
t59 = Ifges(6,5) * t137 + Ifges(6,6) * t189 + Ifges(6,3) * t136;
t58 = -mrSges(7,1) * t189 - mrSges(7,3) * t74;
t57 = mrSges(7,2) * t189 + mrSges(7,3) * t73;
t51 = Ifges(7,1) * t110 + Ifges(7,4) * t111 - Ifges(7,5) * t292;
t50 = Ifges(7,4) * t110 + Ifges(7,2) * t111 - Ifges(7,6) * t292;
t42 = pkin(4) * t136 - t259;
t36 = Ifges(6,4) * t105 + Ifges(6,2) * t138 + Ifges(6,6) * t104;
t35 = Ifges(5,5) * t105 - Ifges(5,6) * t104 + Ifges(5,3) * t138;
t33 = -pkin(4) * t189 - t40;
t31 = -mrSges(7,1) * t73 + mrSges(7,2) * t74;
t30 = -t136 * t331 + t259;
t29 = Ifges(7,1) * t74 + Ifges(7,4) * t73 - Ifges(7,5) * t189;
t28 = Ifges(7,4) * t74 + Ifges(7,2) * t73 - Ifges(7,6) * t189;
t27 = Ifges(7,5) * t74 + Ifges(7,6) * t73 - Ifges(7,3) * t189;
t20 = -qJD(6) * t66 - t252 * t91 + t255 * t86;
t19 = qJD(6) * t65 + t252 * t86 + t255 * t91;
t18 = -mrSges(7,1) * t138 - mrSges(7,3) * t26;
t17 = mrSges(7,2) * t138 + mrSges(7,3) * t25;
t13 = -pkin(4) * t138 - t15;
t11 = pkin(5) * t104 + t14;
t10 = pkin(10) * t104 + t12;
t9 = -pkin(10) * t105 - t138 * t331 - t15;
t5 = Ifges(7,1) * t26 + Ifges(7,4) * t25 - Ifges(7,5) * t138;
t4 = Ifges(7,4) * t26 + Ifges(7,2) * t25 - Ifges(7,6) * t138;
t2 = -qJD(6) * t7 - t10 * t252 + t255 * t9;
t1 = qJD(6) * t6 + t10 * t255 + t252 * t9;
t23 = [0.2e1 * m(4) * (t172 * t182 + t340 * t54 + t55 * t92) + 0.2e1 * t340 * t114 + 0.2e1 * m(3) * (t181 * t199 - t182 * t198) + (t72 + 0.2e1 * t305) * t190 + (-t3 + t35 + t36 - t71 + 0.2e1 * t306) * t189 + (t231 - 0.2e1 * t307 - 0.2e1 * t308) * t251 + (t1 * t7 - t11 * t30 + t2 * t6) * t335 + (t12 * t32 + t13 * t33 + t14 * t42) * t336 + (t15 * t40 + t16 * t41 + t47 * t83) * t337 + (t60 + t61 - t27 - t108) * t138 + (t38 + t39) * t137 + (t34 - t37) * t136 + (-t257 * t283 + 0.2e1 * (t181 * t257 + t182 * t254) * mrSges(3,3) + ((t198 * t334 + Ifges(3,5) * t251 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t257) * t249) * t257 + (t199 * t334 + Ifges(4,5) * t190 - 0.2e1 * Ifges(3,6) * t251 - Ifges(4,6) * t189 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t254 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t257) * t249) * t254) * qJD(2)) * t249 + (t63 + t64) * t105 + (t59 - t62) * t104 + 0.2e1 * t7 * t17 + 0.2e1 * t6 * t18 + t25 * t28 + t26 * t29 + 0.2e1 * t30 * t8 - 0.2e1 * t11 * t31 + 0.2e1 * t42 * t52 + 0.2e1 * t1 * t57 + 0.2e1 * t2 * t58 + 0.2e1 * t32 * t67 + 0.2e1 * t41 * t68 + 0.2e1 * t40 * t69 + 0.2e1 * t33 * t70 + t73 * t4 + t74 * t5 + 0.2e1 * t14 * t75 + 0.2e1 * t47 * t76 + 0.2e1 * t83 * t53 + 0.2e1 * t12 * t87 + 0.2e1 * t16 * t88 + 0.2e1 * t15 * t89 + 0.2e1 * t13 * t90 + 0.2e1 * t92 * t115 + t139 * t109 + 0.2e1 * t54 * t145 + 0.2e1 * t55 * t146 + 0.2e1 * t172 * t77; t231 + ((t286 / 0.2e1 - t108 / 0.2e1 + t60 / 0.2e1 + t61 / 0.2e1 - t27 / 0.2e1 - t340 * mrSges(4,3) + (-m(4) * t340 - t145) * pkin(9)) * t253 + (t109 / 0.2e1 - t92 * mrSges(4,3) + (t64 / 0.2e1 + t63 / 0.2e1) * t250 + (t59 / 0.2e1 - t62 / 0.2e1) * t248 + (-m(4) * t92 + m(5) * t83 - t146 + t76) * pkin(9)) * t256) * qJD(3) - t307 + m(4) * (pkin(9) * t256 * t54 - pkin(2) * t182 - t318 * t55) + m(5) * (t134 * t40 + t135 * t41 + t15 * t160 + t16 * t161 + t318 * t47) + (-t55 * mrSges(4,3) + Ifges(4,5) * t269 + t305 + t72 / 0.2e1 + t284 * t250 + t285 * t248 + (-t115 + t53) * pkin(9)) * t253 + (pkin(9) * t114 + t54 * mrSges(4,3) + Ifges(4,6) * t269 + t3 / 0.2e1 - t35 / 0.2e1 - t36 / 0.2e1 + t71 / 0.2e1 - t306) * t256 + (-t257 * t244 / 0.2e1 - Ifges(3,6) * t293) * t249 + t276 * t137 + t277 * t136 + m(6) * (t113 * t32 + t12 * t143 + t122 * t33 + t125 * t42 + t13 * t144 + t14 * t174) + m(7) * (t1 * t66 - t11 * t142 - t116 * t30 + t19 * t7 + t2 * t65 + t20 * t6) + t51 * t332 + t50 * t333 + t5 * t326 + t4 * t327 - t308 + (t151 / 0.2e1 + t152 / 0.2e1 - t49 / 0.2e1 - t208 / 0.2e1) * t189 + (-t100 / 0.2e1 + t167 / 0.2e1 + t168 / 0.2e1 - t223 / 0.2e1) * t138 + (t166 / 0.2e1 - t169 / 0.2e1) * t104 + (t170 / 0.2e1 + t171 / 0.2e1) * t105 + t30 * t56 + t19 * t57 + t20 * t58 + t65 * t18 + t66 * t17 - pkin(2) * t77 + t6 * t84 + t7 * t85 + t25 * t101 / 0.2e1 + t26 * t102 / 0.2e1 + t110 * t29 / 0.2e1 + t111 * t28 / 0.2e1 - t11 * t112 + t113 * t87 - t116 * t31 + t122 * t90 + t125 * t75 + t134 * t89 + t135 * t88 + t142 * t8 + t143 * t67 + t144 * t70 + t1 * t147 + t2 * t148 + t160 * t69 + t161 * t68 + t174 * t52 + t42 * t178 + t83 * t179 + t14 * t192 + t47 * t193 + t41 * t194 + t40 * t195 + t33 * t196 + t32 * t197 + t16 * t203 + t15 * t204 + t13 * t205 + t12 * t206 + t172 * t207 + t190 * t209 / 0.2e1 + t139 * t224 / 0.2e1; (-t116 * t142 + t19 * t66 + t20 * t65) * t335 + (t113 * t143 + t122 * t144 + t125 * t174) * t336 + (t134 * t160 + t135 * t161) * t337 + (t49 - t151 - t152 + t208) * t256 + ((-t100 + t167 + t168 - t223) * t253 + (t224 + (t170 + t171) * t250 + (t166 - t169) * t248 + (m(5) * t318 + t193) * t345) * t256) * qJD(3) + (t179 * t345 + t209 + (t154 + t155) * t250 + (t150 - t153) * t248) * t253 + 0.2e1 * t65 * t84 + 0.2e1 * t66 * t85 + t110 * t102 + t111 * t101 - 0.2e1 * t116 * t112 + 0.2e1 * t142 * t56 + 0.2e1 * t19 * t147 + 0.2e1 * t20 * t148 + t176 * t50 + t177 * t51 + 0.2e1 * t174 * t178 + 0.2e1 * t125 * t192 + 0.2e1 * t161 * t194 + 0.2e1 * t160 * t195 + 0.2e1 * t144 * t196 + 0.2e1 * t143 * t197 + 0.2e1 * t135 * t203 + 0.2e1 * t134 * t204 + 0.2e1 * t122 * t205 + 0.2e1 * t113 * t206 - 0.2e1 * pkin(2) * t207; t283 + t338 * t14 + (-t1 * t260 + t187 * t6 - t188 * t7 - t2 * t202) * mrSges(7,3) + t342 * pkin(3) + (t13 * mrSges(6,2) - t15 * mrSges(5,3) + (-t75 + t31) * qJD(5) + (-t89 + t90) * qJD(4) + (-t69 + t70) * qJ(4) + m(6) * (qJ(4) * t13 + qJD(4) * t33 - qJD(5) * t42) + m(5) * (-qJ(4) * t15 - qJD(4) * t40) + t30 * t316 + t284) * t248 + (t16 * mrSges(5,3) + t12 * mrSges(6,2) + (t87 + t88) * qJD(4) + (t67 + t68) * qJ(4) + m(6) * (qJ(4) * t12 + qJD(4) * t32) + m(5) * (qJ(4) * t16 + qJD(4) * t41) - t285) * t250 + t274 * t105 + t275 * t104 + t268 * t138 + t120 * t333 + t5 * t322 + t4 * t323 + t28 * t324 + t29 * t325 + t26 * t328 + t25 * t329 + t121 * t332 + m(7) * (t1 * t141 - t11 * t191 + t140 * t2 + t6 * t96 + t7 * t95) - t54 * mrSges(4,2) + t55 * mrSges(4,1) + t95 * t57 + t96 * t58 + t30 * t118 - t11 * t127 + t140 * t18 + t141 * t17 - t189 * t119 / 0.2e1 + t191 * t8 + t210 * t52 + t47 * t215; t244 + m(7) * (-t116 * t191 + t140 * t20 + t141 * t19 + t65 * t96 + t66 * t95) + t338 * t125 + (t187 * t65 - t188 * t66 - t19 * t260 - t20 * t202) * mrSges(7,3) + (t122 * mrSges(6,2) - t134 * mrSges(5,3) + (t112 - t192) * qJD(5) + (-t204 + t205) * qJD(4) + (-t195 + t196) * qJ(4) + m(5) * (-qJ(4) * t134 - qJD(4) * t160) + m(6) * (qJ(4) * t122 + qJD(4) * t144 - qJD(5) * t174) + t142 * t316 + t276) * t248 + (t135 * mrSges(5,3) + t113 * mrSges(6,2) + (t206 + t203) * qJD(4) + (t194 + t197) * qJ(4) + m(5) * (qJ(4) * t135 + qJD(4) * t161) + m(6) * (qJ(4) * t113 + qJD(4) * t143) - t277) * t250 + ((pkin(9) * mrSges(4,2) - Ifges(4,6) + t268) * t253 + (t274 * t250 + t275 * t248 + (-m(5) * pkin(3) - mrSges(4,1) + t215) * pkin(9)) * t256) * qJD(3) + t51 * t322 + t50 * t323 + t101 * t324 + t102 * t325 + t121 * t326 + t120 * t327 + t110 * t328 + t111 * t329 + t256 * t119 / 0.2e1 - t116 * t127 + t140 * t84 + t141 * t85 + t142 * t118 + t95 * t147 + t96 * t148 - pkin(3) * t179 + t191 * t56 + t210 * t178; 0.2e1 * t191 * t118 - t260 * t120 + t202 * t121 - t188 * t129 - t187 * t130 + (t140 * t96 + t141 * t95 + t191 * t290) * t335 + (-t210 * t290 + t296) * t336 + t296 * t337 + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * t346 + 0.2e1 * (t127 - t214) * t290 + 0.2e1 * (t140 * t187 - t141 * t188 - t202 * t96 - t260 * t95) * mrSges(7,3); m(6) * t14 + m(7) * t11 - t342 + t52 - t8; t226 + (m(5) * pkin(9) - t315) * t291 + m(6) * t125 + m(7) * t116 - t56 + t179; (-m(6) - m(7)) * t290 - t118; 0; t252 * t17 + t255 * t18 + (-t252 * t58 + t255 * t57) * qJD(6) + m(7) * (t1 * t252 + t2 * t255 + (-t252 * t6 + t255 * t7) * qJD(6)) + m(6) * t13 + t70; t252 * t85 + t255 * t84 + (t147 * t255 - t148 * t252) * qJD(6) + m(7) * (t19 * t252 + t20 * t255 + (-t252 * t65 + t255 * t66) * qJD(6)) + m(6) * t122 + t196; m(7) * (t252 * t95 + t255 * t96 + (-t140 * t252 + t141 * t255) * qJD(6)) + m(6) * t248 * qJD(4) + (t255 * t187 - t252 * t188 + (t202 * t252 - t255 * t260) * qJD(6)) * mrSges(7,3); 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t3; mrSges(7,1) * t20 - mrSges(7,2) * t19 + t49; mrSges(7,1) * t96 - mrSges(7,2) * t95 + t119; 0; (-mrSges(7,1) * t252 - mrSges(7,2) * t255) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;
