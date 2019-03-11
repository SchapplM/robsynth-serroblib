% Calculate time derivative of joint inertia matrix for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:31
% EndTime: 2019-03-09 16:21:45
% DurationCPUTime: 5.66s
% Computational Cost: add. (7084->645), mult. (18508->923), div. (0->0), fcn. (17192->10), ass. (0->264)
t326 = Ifges(5,1) + Ifges(4,3);
t325 = Ifges(4,5) - Ifges(5,4);
t324 = -Ifges(4,6) + Ifges(5,5);
t237 = sin(pkin(11));
t239 = cos(pkin(11));
t242 = sin(qJ(6));
t245 = cos(qJ(6));
t255 = t245 * t237 + t242 * t239;
t303 = -t255 / 0.2e1;
t322 = -t237 * t242 + t239 * t245;
t302 = t322 / 0.2e1;
t323 = -t237 / 0.2e1;
t300 = t239 / 0.2e1;
t238 = sin(pkin(6));
t244 = sin(qJ(2));
t282 = t238 * t244;
t221 = pkin(8) * t282;
t240 = cos(pkin(6));
t247 = cos(qJ(2));
t299 = pkin(1) * t247;
t182 = t240 * t299 - t221;
t174 = t255 * qJD(6);
t264 = (t237 ^ 2 + t239 ^ 2) * qJD(5);
t321 = 2 * m(6);
t320 = 2 * m(7);
t319 = -2 * mrSges(3,3);
t318 = 2 * mrSges(3,3);
t246 = cos(qJ(3));
t277 = qJD(2) * t238;
t268 = t247 * t277;
t243 = sin(qJ(3));
t270 = t243 * t282;
t134 = -qJD(3) * t270 + (qJD(3) * t240 + t268) * t246;
t176 = t240 * t243 + t246 * t282;
t133 = qJD(3) * t176 + t243 * t268;
t276 = qJD(2) * t244;
t269 = t238 * t276;
t94 = t133 * t239 - t237 * t269;
t95 = t133 * t237 + t239 * t269;
t38 = Ifges(6,4) * t95 + Ifges(6,2) * t94 + Ifges(6,6) * t134;
t317 = -t38 / 0.2e1;
t175 = -t240 * t246 + t270;
t281 = t238 * t247;
t131 = t175 * t239 + t237 * t281;
t132 = t175 * t237 - t239 * t281;
t70 = t131 * t245 - t132 * t242;
t316 = t70 / 0.2e1;
t71 = t131 * t242 + t132 * t245;
t315 = t71 / 0.2e1;
t314 = t94 / 0.2e1;
t313 = t95 / 0.2e1;
t312 = pkin(4) + pkin(9);
t173 = t322 * qJD(6);
t113 = -Ifges(7,5) * t174 - Ifges(7,6) * t173;
t311 = t113 / 0.2e1;
t124 = Ifges(7,4) * t322 - Ifges(7,2) * t255;
t310 = t124 / 0.2e1;
t125 = Ifges(7,1) * t322 - Ifges(7,4) * t255;
t309 = t125 / 0.2e1;
t292 = Ifges(6,4) * t239;
t262 = Ifges(6,1) * t237 + t292;
t147 = (Ifges(6,5) * t246 + t243 * t262) * qJD(3);
t308 = t147 / 0.2e1;
t161 = t322 * t246;
t307 = -t161 / 0.2e1;
t162 = t255 * t246;
t306 = -t162 / 0.2e1;
t305 = -t173 / 0.2e1;
t304 = -t174 / 0.2e1;
t301 = t237 / 0.2e1;
t233 = t246 * pkin(9);
t298 = pkin(3) + qJ(5);
t297 = -pkin(10) - t298;
t183 = t240 * t244 * pkin(1) + pkin(8) * t281;
t159 = pkin(9) * t240 + t183;
t160 = (-pkin(2) * t247 - pkin(9) * t244 - pkin(1)) * t238;
t165 = (pkin(2) * t244 - pkin(9) * t247) * t277;
t166 = t182 * qJD(2);
t274 = qJD(3) * t246;
t275 = qJD(3) * t243;
t49 = -t159 * t274 - t160 * t275 + t246 * t165 - t243 * t166;
t33 = t134 * pkin(4) + (qJD(5) * t247 - t276 * t298) * t238 - t49;
t167 = t183 * qJD(2);
t248 = -t134 * qJ(4) - t176 * qJD(4) + t167;
t34 = t175 * qJD(5) + t133 * t298 + t248;
t12 = t237 * t33 + t239 * t34;
t88 = -t243 * t159 + t246 * t160;
t79 = pkin(3) * t281 - t88;
t60 = t176 * pkin(4) + qJ(5) * t281 + t79;
t158 = t221 + (-pkin(2) - t299) * t240;
t252 = -t176 * qJ(4) + t158;
t64 = t175 * t298 + t252;
t26 = t237 * t60 + t239 * t64;
t296 = mrSges(6,2) * t237;
t295 = Ifges(4,4) * t243;
t294 = Ifges(4,4) * t246;
t293 = Ifges(6,4) * t237;
t291 = Ifges(5,6) * t243;
t290 = Ifges(5,6) * t246;
t289 = t166 * mrSges(3,2);
t288 = t167 * mrSges(3,1);
t287 = t167 * mrSges(4,1);
t286 = t167 * mrSges(4,2);
t285 = qJ(4) * t243;
t283 = t237 * t243;
t279 = t239 * t246;
t263 = pkin(3) * t275 - qJD(4) * t243;
t135 = -qJD(5) * t246 + (-qJ(4) * t246 + qJ(5) * t243) * qJD(3) + t263;
t198 = t312 * t274;
t87 = t239 * t135 + t237 * t198;
t89 = t246 * t159 + t243 * t160;
t184 = -t246 * t298 - pkin(2) - t285;
t212 = t312 * t243;
t118 = t239 * t184 + t237 * t212;
t213 = t246 * pkin(4) + t233;
t273 = 2 * mrSges(7,3);
t272 = 0.2e1 * t238;
t23 = qJD(6) * t70 + t242 * t94 + t245 * t95;
t24 = -qJD(6) * t71 - t242 * t95 + t245 * t94;
t6 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t134;
t104 = -qJD(6) * t161 + t255 * t275;
t105 = t174 * t246 + t275 * t322;
t44 = Ifges(7,5) * t104 + Ifges(7,6) * t105 + Ifges(7,3) * t274;
t267 = t239 * t275;
t266 = t281 / 0.2e1;
t265 = Ifges(6,5) * t300 + Ifges(7,5) * t302 + Ifges(6,6) * t323 + Ifges(7,6) * t303;
t47 = -t94 * mrSges(6,1) + t95 * mrSges(6,2);
t10 = -t24 * mrSges(7,1) + t23 * mrSges(7,2);
t52 = -t105 * mrSges(7,1) + t104 * mrSges(7,2);
t11 = -t237 * t34 + t239 * t33;
t25 = -t237 * t64 + t239 * t60;
t108 = t134 * mrSges(5,1) + mrSges(5,2) * t269;
t112 = mrSges(7,1) * t173 - t174 * mrSges(7,2);
t86 = -t135 * t237 + t239 * t198;
t207 = t246 * mrSges(5,2) - t243 * mrSges(5,3);
t261 = Ifges(6,2) * t239 + t293;
t260 = Ifges(6,5) * t237 + Ifges(6,6) * t239;
t259 = -pkin(3) * t246 - t285;
t258 = t11 * t239 + t12 * t237;
t15 = pkin(5) * t176 - pkin(10) * t132 + t25;
t18 = pkin(10) * t131 + t26;
t3 = t15 * t245 - t18 * t242;
t4 = t15 * t242 + t18 * t245;
t257 = t237 * t87 + t239 * t86;
t111 = -pkin(10) * t279 + t118;
t188 = t239 * t212;
t99 = pkin(5) * t243 + t188 + (pkin(10) * t246 - t184) * t237;
t51 = t111 * t245 + t242 * t99;
t50 = -t111 * t242 + t245 * t99;
t256 = -t173 * t255 + t174 * t322;
t199 = t297 * t237;
t200 = t297 * t239;
t121 = t199 * t245 + t200 * t242;
t120 = -t199 * t242 + t200 * t245;
t164 = -mrSges(6,1) * t267 + t275 * t296;
t78 = qJ(4) * t281 - t89;
t253 = t324 * t133 + t325 * t134 + t269 * t326;
t48 = -t159 * t275 + t160 * t274 + t243 * t165 + t246 * t166;
t65 = -t175 * pkin(4) - t78;
t5 = pkin(5) * t134 - pkin(10) * t95 + t11;
t9 = pkin(10) * t94 + t12;
t1 = qJD(6) * t3 + t242 * t5 + t245 * t9;
t2 = -qJD(6) * t4 - t242 * t9 + t245 * t5;
t251 = -t1 * t255 - t173 * t4 + t174 * t3 - t2 * t322;
t75 = (pkin(5) * t246 - pkin(10) * t283) * qJD(3) + t86;
t76 = pkin(10) * t267 + t87;
t13 = qJD(6) * t50 + t242 * t75 + t245 * t76;
t14 = -qJD(6) * t51 - t242 * t76 + t245 * t75;
t250 = -t13 * t255 - t14 * t322 - t173 * t51 + t174 * t50;
t84 = -qJD(5) * t255 + qJD(6) * t120;
t85 = -qJD(5) * t322 - qJD(6) * t121;
t249 = -t120 * t174 + t121 * t173 + t255 * t84 + t322 * t85;
t40 = -qJ(4) * t269 + qJD(4) * t281 - t48;
t35 = -t133 * pkin(4) - t40;
t231 = Ifges(4,5) * t274;
t230 = Ifges(5,5) * t275;
t225 = pkin(5) * t237 + qJ(4);
t218 = Ifges(3,5) * t268;
t211 = Ifges(4,1) * t243 + t294;
t210 = Ifges(4,2) * t246 + t295;
t209 = -Ifges(5,2) * t243 - t290;
t208 = -Ifges(5,3) * t246 - t291;
t206 = Ifges(6,1) * t239 - t293;
t205 = -Ifges(6,2) * t237 + t292;
t203 = mrSges(6,1) * t237 + mrSges(6,2) * t239;
t201 = -pkin(2) + t259;
t197 = t312 * t275;
t196 = (Ifges(4,1) * t246 - t295) * qJD(3);
t195 = (-Ifges(4,2) * t243 + t294) * qJD(3);
t194 = (-Ifges(5,2) * t246 + t291) * qJD(3);
t193 = (Ifges(5,3) * t243 - t290) * qJD(3);
t192 = (mrSges(4,1) * t243 + mrSges(4,2) * t246) * qJD(3);
t191 = (-mrSges(5,2) * t243 - mrSges(5,3) * t246) * qJD(3);
t190 = -mrSges(6,2) * t243 - mrSges(6,3) * t279;
t189 = mrSges(6,3) * t237 * t246 + mrSges(6,1) * t243;
t181 = (mrSges(6,3) * t239 * t243 - mrSges(6,2) * t246) * qJD(3);
t180 = (mrSges(6,1) * t246 - mrSges(6,3) * t283) * qJD(3);
t179 = (mrSges(6,1) * t239 - t296) * t246;
t172 = pkin(5) * t279 + t213;
t171 = -qJ(4) * t274 + t263;
t157 = Ifges(6,5) * t243 - t246 * t262;
t156 = Ifges(6,6) * t243 - t246 * t261;
t155 = Ifges(6,3) * t243 - t246 * t260;
t154 = (-pkin(5) * t239 - t312) * t275;
t146 = (Ifges(6,6) * t246 + t243 * t261) * qJD(3);
t145 = (Ifges(6,3) * t246 + t243 * t260) * qJD(3);
t141 = mrSges(7,1) * t243 + mrSges(7,3) * t162;
t140 = -mrSges(7,2) * t243 - mrSges(7,3) * t161;
t139 = -mrSges(4,1) * t281 - t176 * mrSges(4,3);
t138 = mrSges(4,2) * t281 - t175 * mrSges(4,3);
t137 = t176 * mrSges(5,1) - mrSges(5,2) * t281;
t136 = t175 * mrSges(5,1) + mrSges(5,3) * t281;
t122 = mrSges(7,1) * t255 + mrSges(7,2) * t322;
t117 = -t184 * t237 + t188;
t116 = -mrSges(5,2) * t175 - mrSges(5,3) * t176;
t115 = -Ifges(7,1) * t174 - Ifges(7,4) * t173;
t114 = -Ifges(7,4) * t174 - Ifges(7,2) * t173;
t110 = mrSges(4,1) * t269 - mrSges(4,3) * t134;
t109 = -mrSges(4,2) * t269 - mrSges(4,3) * t133;
t107 = mrSges(5,1) * t133 - mrSges(5,3) * t269;
t106 = mrSges(7,1) * t161 - mrSges(7,2) * t162;
t103 = Ifges(4,1) * t176 - Ifges(4,4) * t175 - Ifges(4,5) * t281;
t102 = Ifges(4,4) * t176 - Ifges(4,2) * t175 - Ifges(4,6) * t281;
t101 = -Ifges(5,4) * t281 - Ifges(5,2) * t176 + Ifges(5,6) * t175;
t100 = -Ifges(5,5) * t281 - Ifges(5,6) * t176 + Ifges(5,3) * t175;
t93 = -Ifges(7,1) * t162 - Ifges(7,4) * t161 + Ifges(7,5) * t243;
t92 = -Ifges(7,4) * t162 - Ifges(7,2) * t161 + Ifges(7,6) * t243;
t91 = -Ifges(7,5) * t162 - Ifges(7,6) * t161 + Ifges(7,3) * t243;
t83 = mrSges(6,1) * t176 - mrSges(6,3) * t132;
t82 = -mrSges(6,2) * t176 + mrSges(6,3) * t131;
t81 = -mrSges(7,2) * t274 + mrSges(7,3) * t105;
t80 = mrSges(7,1) * t274 - mrSges(7,3) * t104;
t77 = t175 * pkin(3) + t252;
t74 = -mrSges(5,2) * t133 - mrSges(5,3) * t134;
t73 = mrSges(4,1) * t133 + mrSges(4,2) * t134;
t72 = -mrSges(6,1) * t131 + mrSges(6,2) * t132;
t69 = Ifges(4,1) * t134 - Ifges(4,4) * t133 + Ifges(4,5) * t269;
t68 = Ifges(4,4) * t134 - Ifges(4,2) * t133 + Ifges(4,6) * t269;
t67 = Ifges(5,4) * t269 - Ifges(5,2) * t134 + Ifges(5,6) * t133;
t66 = Ifges(5,5) * t269 - Ifges(5,6) * t134 + Ifges(5,3) * t133;
t62 = mrSges(6,1) * t134 - mrSges(6,3) * t95;
t61 = -mrSges(6,2) * t134 + mrSges(6,3) * t94;
t59 = Ifges(6,1) * t132 + Ifges(6,4) * t131 + Ifges(6,5) * t176;
t58 = Ifges(6,4) * t132 + Ifges(6,2) * t131 + Ifges(6,6) * t176;
t57 = Ifges(6,5) * t132 + Ifges(6,6) * t131 + Ifges(6,3) * t176;
t54 = mrSges(7,1) * t176 - mrSges(7,3) * t71;
t53 = -mrSges(7,2) * t176 + mrSges(7,3) * t70;
t46 = Ifges(7,1) * t104 + Ifges(7,4) * t105 + Ifges(7,5) * t274;
t45 = Ifges(7,4) * t104 + Ifges(7,2) * t105 + Ifges(7,6) * t274;
t43 = t133 * pkin(3) + t248;
t42 = -pkin(3) * t269 - t49;
t41 = -t131 * pkin(5) + t65;
t39 = Ifges(6,1) * t95 + Ifges(6,4) * t94 + Ifges(6,5) * t134;
t37 = Ifges(6,5) * t95 + Ifges(6,6) * t94 + Ifges(6,3) * t134;
t36 = -mrSges(7,1) * t70 + mrSges(7,2) * t71;
t29 = Ifges(7,1) * t71 + Ifges(7,4) * t70 + Ifges(7,5) * t176;
t28 = Ifges(7,4) * t71 + Ifges(7,2) * t70 + Ifges(7,6) * t176;
t27 = Ifges(7,5) * t71 + Ifges(7,6) * t70 + Ifges(7,3) * t176;
t19 = -t94 * pkin(5) + t35;
t17 = -mrSges(7,2) * t134 + mrSges(7,3) * t24;
t16 = mrSges(7,1) * t134 - mrSges(7,3) * t23;
t8 = Ifges(7,1) * t23 + Ifges(7,4) * t24 + Ifges(7,5) * t134;
t7 = Ifges(7,4) * t23 + Ifges(7,2) * t24 + Ifges(7,6) * t134;
t20 = [(t1 * t4 + t19 * t41 + t2 * t3) * t320 + (t11 * t25 + t12 * t26 + t35 * t65) * t321 + (t167 * t244 * t318 + (t166 * t318 - t253) * t247 + ((t182 * t319 + Ifges(3,5) * t240 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t247) * t272) * t247 + (t183 * t319 - 0.2e1 * Ifges(3,6) * t240 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t244) * t272 + t325 * t176 + t324 * t175 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t326) * t281) * t244) * qJD(2)) * t238 + 0.2e1 * m(3) * (t166 * t183 - t167 * t182) + 0.2e1 * m(4) * (t158 * t167 + t48 * t89 + t49 * t88) + 0.2e1 * m(5) * (t40 * t78 + t42 * t79 + t43 * t77) + (t27 + t57 + t103 - t101) * t134 + (t100 - t102) * t133 + 0.2e1 * t3 * t16 + 0.2e1 * t4 * t17 + t24 * t28 + t23 * t29 + 0.2e1 * t19 * t36 + 0.2e1 * t41 * t10 + 0.2e1 * t1 * t53 + 0.2e1 * t2 * t54 + 0.2e1 * t26 * t61 + 0.2e1 * t25 * t62 + 0.2e1 * t65 * t47 + t70 * t7 + t71 * t8 + 0.2e1 * t35 * t72 + 0.2e1 * t77 * t74 + 0.2e1 * t12 * t82 + 0.2e1 * t11 * t83 + (t37 + t6 - t67 + t69 + 0.2e1 * t286) * t176 + (t66 - t68 + 0.2e1 * t287) * t175 + (t218 - 0.2e1 * t288 - 0.2e1 * t289) * t240 + t94 * t58 + t95 * t59 + 0.2e1 * t78 * t107 + 0.2e1 * t79 * t108 + 0.2e1 * t89 * t109 + 0.2e1 * t88 * t110 + 0.2e1 * t43 * t116 + t131 * t38 + t132 * t39 + 0.2e1 * t40 * t136 + 0.2e1 * t42 * t137 + 0.2e1 * t48 * t138 + 0.2e1 * t49 * t139 + 0.2e1 * t158 * t73; m(6) * (t11 * t117 + t118 * t12 - t197 * t65 + t213 * t35 + t25 * t86 + t26 * t87) + (-m(4) * t167 - t73) * pkin(2) - t288 + m(7) * (t1 * t51 + t13 * t4 + t14 * t3 + t154 * t41 + t172 * t19 + t2 * t50) + t218 + (t44 / 0.2e1 + t145 / 0.2e1 - t194 / 0.2e1 + t196 / 0.2e1) * t176 + (t193 / 0.2e1 - t195 / 0.2e1) * t175 + (t91 / 0.2e1 + t155 / 0.2e1 - t209 / 0.2e1 + t211 / 0.2e1) * t134 + (t208 / 0.2e1 - t210 / 0.2e1) * t133 + t45 * t316 + t8 * t306 + t7 * t307 + t132 * t308 + t157 * t313 + t156 * t314 + t46 * t315 + m(5) * (t171 * t77 + t201 * t43) + (t239 * t317 + t39 * t323 + t48 * mrSges(4,3) - t40 * mrSges(5,1) - t66 / 0.2e1 + t68 / 0.2e1 - t287) * t246 - t289 + t50 * t16 + t51 * t17 + t41 * t52 + ((-t107 + t109) * t246 + (t108 - t110) * t243 + ((t137 - t139) * t246 + (t136 - t138) * t243) * qJD(3) + m(5) * (t243 * t42 - t246 * t40 + t274 * t79 + t275 * t78) + m(4) * (-t243 * t49 + t246 * t48 - t274 * t88 - t275 * t89)) * pkin(9) + ((-t231 / 0.2e1 - t230 / 0.2e1) * t247 + (-Ifges(3,6) + (-Ifges(5,5) / 0.2e1 + Ifges(4,6) / 0.2e1) * t246 + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t243) * t276) * t238 + t13 * t53 + t14 * t54 + t3 * t80 + t4 * t81 + t86 * t83 + t87 * t82 + (t42 * mrSges(5,1) - t49 * mrSges(4,3) + t6 / 0.2e1 + t37 / 0.2e1 - t67 / 0.2e1 + t69 / 0.2e1 + t286) * t243 + t24 * t92 / 0.2e1 + t23 * t93 / 0.2e1 + t104 * t29 / 0.2e1 + t105 * t28 / 0.2e1 + t19 * t106 + t117 * t62 + t118 * t61 + t1 * t140 + t2 * t141 + t131 * t146 / 0.2e1 + t154 * t36 + ((-t101 / 0.2e1 + t103 / 0.2e1 + t57 / 0.2e1 + t27 / 0.2e1 + Ifges(5,4) * t266 + t79 * mrSges(5,1) - t88 * mrSges(4,3)) * t246 + (Ifges(4,6) * t266 + t100 / 0.2e1 - t102 / 0.2e1 + t78 * mrSges(5,1) - t89 * mrSges(4,3) + t59 * t301 + t58 * t300) * t243) * qJD(3) + t65 * t164 + t171 * t116 + t172 * t10 + t35 * t179 + t25 * t180 + t26 * t181 + t11 * t189 + t12 * t190 + t77 * t191 + t158 * t192 - t197 * t72 + t201 * t74 + t43 * t207 + t213 * t47; 0.2e1 * t50 * t80 + 0.2e1 * t51 * t81 + t104 * t93 + t105 * t92 + 0.2e1 * t13 * t140 + 0.2e1 * t14 * t141 + 0.2e1 * t154 * t106 - t161 * t45 - t162 * t46 + 0.2e1 * t172 * t52 + 0.2e1 * t117 * t180 + 0.2e1 * t118 * t181 + 0.2e1 * t86 * t189 + 0.2e1 * t87 * t190 - 0.2e1 * pkin(2) * t192 - 0.2e1 * t197 * t179 + 0.2e1 * t201 * t191 + 0.2e1 * t213 * t164 + 0.2e1 * (m(5) * t201 + t207) * t171 + (t117 * t86 + t118 * t87 - t197 * t213) * t321 + (t13 * t51 + t14 * t50 + t154 * t172) * t320 + (-t146 * t239 - t147 * t237 - t193 + t195) * t246 + (t44 + t145 - t194 + t196) * t243 + ((t155 - t209 + t211 + t91) * t246 + (t156 * t239 + t157 * t237 + t208 - t210) * t243) * qJD(3); (-qJD(5) * t83 - t11 * mrSges(6,3) - t298 * t62 + t39 / 0.2e1) * t239 + m(6) * (qJ(4) * t35 + qJD(4) * t65 - t258 * t298 + (-t237 * t26 - t239 * t25) * qJD(5)) + (-t12 * mrSges(6,3) - qJD(5) * t82 - t298 * t61 + t317) * t237 + t253 + (t36 + t72 - t136) * qJD(4) + (t47 - t107) * qJ(4) + t8 * t302 + t7 * t303 + t29 * t304 + t28 * t305 + t23 * t309 + t24 * t310 + t176 * t311 + t206 * t313 + t205 * t314 + t115 * t315 + t114 * t316 + m(7) * (qJD(4) * t41 + t1 * t121 + t120 * t2 + t19 * t225 + t3 * t85 + t4 * t84) + m(5) * (-pkin(3) * t42 - qJ(4) * t40 - qJD(4) * t78) + t251 * mrSges(7,3) + t265 * t134 - t40 * mrSges(5,3) + t42 * mrSges(5,2) - t48 * mrSges(4,2) + t49 * mrSges(4,1) + t84 * t53 + t85 * t54 - pkin(3) * t108 + t41 * t112 + t120 * t16 + t121 * t17 + t19 * t122 + t35 * t203 + t225 * t10; t250 * mrSges(7,3) + m(6) * (-qJ(4) * t197 - t257 * t298 + (-t117 * t239 - t118 * t237) * qJD(5)) + (-t298 * t181 - t87 * mrSges(6,3) - qJD(5) * t190 - t146 / 0.2e1) * t237 + (-t86 * mrSges(6,3) - qJD(5) * t189 - t180 * t298 + t308) * t239 + t230 + t231 + t46 * t302 + t45 * t303 + t93 * t304 + t92 * t305 + t115 * t306 + t114 * t307 + t104 * t309 + t105 * t310 + t243 * t311 + m(7) * (t120 * t14 + t121 * t13 + t154 * t225 + t50 * t85 + t51 * t84) + (m(5) * t233 + m(6) * t213 + m(7) * t172 + t246 * mrSges(5,1) + t106 + t179) * qJD(4) + t120 * t80 + t121 * t81 + t84 * t140 + t85 * t141 + t154 * t122 + ((-pkin(3) * mrSges(5,1) - Ifges(5,4) + t265) * t246 + (-qJ(4) * mrSges(5,1) + t205 * t300 + t206 * t301 - Ifges(4,6)) * t243 + (m(5) * t259 - t246 * mrSges(4,1) + t243 * mrSges(4,2) + t207) * pkin(9)) * qJD(3) + qJ(4) * t164 + t172 * t112 - t197 * t203 + t225 * t52; 0.2e1 * t225 * t112 - t255 * t114 + t322 * t115 - t173 * t124 - t174 * t125 + (qJD(4) * t225 + t120 * t85 + t121 * t84) * t320 + (qJ(4) * qJD(4) + t264 * t298) * t321 - t249 * t273 + 0.2e1 * mrSges(6,3) * t264 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3) + t122 + t203) * qJD(4); m(5) * t42 + m(6) * t258 - m(7) * t251 + t16 * t322 + t17 * t255 + t173 * t53 - t174 * t54 + t237 * t61 + t239 * t62 + t108; t173 * t140 - t174 * t141 + t239 * t180 + t237 * t181 + t255 * t81 + t322 * t80 + (m(5) * pkin(9) + mrSges(5,1)) * t274 - m(7) * t250 + m(6) * t257; -m(6) * t264 + m(7) * t249 + t256 * t273; -0.2e1 * m(7) * t256; m(6) * t35 + m(7) * t19 + t10 + t47; -m(6) * t197 + m(7) * t154 + t164 + t52; (m(6) + m(7)) * qJD(4) + t112; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t6; mrSges(7,1) * t14 - mrSges(7,2) * t13 + t44; mrSges(7,1) * t85 - mrSges(7,2) * t84 + t113; -mrSges(7,1) * t174 - t173 * mrSges(7,2); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
