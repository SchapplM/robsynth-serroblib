% Calculate time derivative of joint inertia matrix for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:39:11
% EndTime: 2019-03-08 23:39:23
% DurationCPUTime: 5.17s
% Computational Cost: add. (7307->633), mult. (21038->952), div. (0->0), fcn. (20751->14), ass. (0->258)
t223 = sin(pkin(13));
t226 = cos(pkin(13));
t316 = mrSges(6,1) * t223 + mrSges(6,2) * t226;
t310 = 2 * pkin(10);
t229 = sin(qJ(6));
t233 = cos(qJ(6));
t239 = t223 * t229 - t226 * t233;
t292 = -t239 / 0.2e1;
t185 = t223 * t233 + t226 * t229;
t291 = t185 / 0.2e1;
t289 = t226 / 0.2e1;
t224 = sin(pkin(7));
t231 = sin(qJ(3));
t257 = qJD(3) * t231;
t247 = t224 * t257;
t227 = cos(pkin(7));
t235 = cos(qJ(3));
t268 = t224 * t235;
t180 = t227 * t231 * pkin(2) + pkin(9) * t268;
t159 = pkin(10) * t227 + t180;
t160 = (-pkin(3) * t235 - pkin(10) * t231 - pkin(2)) * t224;
t258 = qJD(3) * t224;
t165 = (pkin(3) * t231 - pkin(10) * t235) * t258;
t269 = t224 * t231;
t212 = pkin(9) * t269;
t288 = pkin(2) * t235;
t179 = t227 * t288 - t212;
t166 = t179 * qJD(3);
t230 = sin(qJ(4));
t234 = cos(qJ(4));
t255 = qJD(4) * t234;
t256 = qJD(4) * t230;
t56 = -t159 * t255 - t160 * t256 + t165 * t234 - t230 * t166;
t47 = -pkin(4) * t247 - t56;
t174 = -t234 * t227 + t230 * t269;
t246 = t235 * t258;
t129 = -t174 * qJD(4) + t234 * t246;
t97 = -t129 * t223 + t226 * t247;
t98 = t129 * t226 + t223 * t247;
t54 = -t97 * mrSges(6,1) + t98 * mrSges(6,2);
t315 = -m(6) * t47 - t54;
t236 = cos(qJ(2));
t262 = t235 * t236;
t232 = sin(qJ(2));
t265 = t231 * t232;
t314 = t227 * t262 - t265;
t313 = -m(5) * pkin(3) - t234 * mrSges(5,1) + t230 * mrSges(5,2) - mrSges(4,1);
t172 = t239 * qJD(6);
t312 = 0.2e1 * m(6);
t311 = 2 * m(7);
t309 = -2 * mrSges(4,3);
t308 = m(6) / 0.2e1;
t306 = m(6) * pkin(10);
t175 = t227 * t230 + t234 * t269;
t127 = -t175 * t223 - t226 * t268;
t128 = t175 * t226 - t223 * t268;
t70 = t127 * t233 - t128 * t229;
t305 = t70 / 0.2e1;
t71 = t127 * t229 + t128 * t233;
t304 = t71 / 0.2e1;
t303 = t97 / 0.2e1;
t302 = t98 / 0.2e1;
t120 = Ifges(7,4) * t185 - Ifges(7,2) * t239;
t300 = t120 / 0.2e1;
t121 = Ifges(7,1) * t185 - Ifges(7,4) * t239;
t299 = t121 / 0.2e1;
t281 = Ifges(6,4) * t226;
t242 = -Ifges(6,2) * t223 + t281;
t142 = (Ifges(6,6) * t230 + t242 * t234) * qJD(4);
t298 = t142 / 0.2e1;
t282 = Ifges(6,4) * t223;
t243 = Ifges(6,1) * t226 - t282;
t143 = (Ifges(6,5) * t230 + t243 * t234) * qJD(4);
t297 = t143 / 0.2e1;
t161 = t185 * t230;
t296 = -t161 / 0.2e1;
t162 = t239 * t230;
t295 = -t162 / 0.2e1;
t294 = -t172 / 0.2e1;
t173 = t185 * qJD(6);
t293 = -t173 / 0.2e1;
t290 = -t223 / 0.2e1;
t225 = sin(pkin(6));
t259 = qJD(2) * t225;
t248 = t232 * t259;
t244 = t224 * t248;
t228 = cos(pkin(6));
t75 = t228 * t246 + (t314 * qJD(3) + (-t227 * t265 + t262) * qJD(2)) * t225;
t263 = t232 * t235;
t264 = t231 * t236;
t237 = t227 * t264 + t263;
t117 = t237 * t225 + t228 * t269;
t171 = -t224 * t225 * t236 + t228 * t227;
t81 = t117 * t234 + t171 * t230;
t36 = t81 * qJD(4) + t230 * t75 - t234 * t244;
t80 = t117 * t230 - t171 * t234;
t21 = t80 * t36;
t287 = pkin(11) + qJ(5);
t55 = -t159 * t256 + t160 * t255 + t230 * t165 + t234 * t166;
t44 = (qJ(5) * t257 - qJD(5) * t235) * t224 + t55;
t130 = t175 * qJD(4) + t230 * t246;
t167 = t180 * qJD(3);
t48 = pkin(4) * t130 - qJ(5) * t129 - qJD(5) * t175 + t167;
t14 = t223 * t48 + t226 * t44;
t158 = t212 + (-pkin(3) - t288) * t227;
t79 = pkin(4) * t174 - qJ(5) * t175 + t158;
t89 = t234 * t159 + t230 * t160;
t82 = -qJ(5) * t268 + t89;
t42 = t223 * t79 + t226 * t82;
t284 = Ifges(5,4) * t230;
t283 = Ifges(5,4) * t234;
t116 = -t225 * t314 - t228 * t268;
t74 = t228 * t247 + (t237 * qJD(3) + (t227 * t263 + t264) * qJD(2)) * t225;
t280 = t116 * t74;
t279 = t36 * t230;
t37 = -t80 * qJD(4) + t230 * t244 + t234 * t75;
t278 = t37 * t234;
t107 = mrSges(5,1) * t247 - mrSges(5,3) * t129;
t276 = -t107 + t54;
t135 = -mrSges(5,1) * t268 - mrSges(5,3) * t175;
t72 = -mrSges(6,1) * t127 + mrSges(6,2) * t128;
t275 = -t135 + t72;
t164 = t316 * t255;
t104 = -t230 * t173 - t239 * t255;
t105 = t172 * t230 - t185 * t255;
t57 = -t105 * mrSges(7,1) + t104 * mrSges(7,2);
t274 = t164 + t57;
t198 = -mrSges(6,1) * t226 + mrSges(6,2) * t223;
t273 = t198 - mrSges(5,1);
t272 = t116 * t167;
t271 = t223 * t230;
t270 = t223 * t234;
t267 = t226 * t230;
t266 = t226 * t234;
t111 = -Ifges(7,5) * t172 - Ifges(7,6) * t173;
t261 = -mrSges(4,1) * t227 + mrSges(5,1) * t174 + mrSges(5,2) * t175 + mrSges(4,3) * t269;
t170 = -qJD(5) * t230 + (pkin(4) * t230 - qJ(5) * t234) * qJD(4);
t253 = pkin(10) * t256;
t125 = t226 * t170 + t223 * t253;
t196 = -pkin(4) * t234 - qJ(5) * t230 - pkin(3);
t148 = pkin(10) * t266 + t223 * t196;
t26 = t70 * qJD(6) + t229 * t97 + t233 * t98;
t27 = -t71 * qJD(6) - t229 * t98 + t233 * t97;
t5 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t130;
t252 = Ifges(5,6) * t268;
t49 = Ifges(7,5) * t104 + Ifges(7,6) * t105 + Ifges(7,3) * t256;
t250 = Ifges(5,5) * t129 - Ifges(5,6) * t130 + Ifges(5,3) * t247;
t249 = pkin(5) * t223 + pkin(10);
t245 = Ifges(6,5) * t223 / 0.2e1 + Ifges(6,6) * t289 + Ifges(7,5) * t291 + Ifges(7,6) * t292;
t8 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t13 = -t223 * t44 + t226 * t48;
t41 = -t223 * t82 + t226 * t79;
t88 = -t230 * t159 + t160 * t234;
t83 = pkin(4) * t268 - t88;
t241 = Ifges(6,5) * t226 - Ifges(6,6) * t223;
t15 = -t223 * t37 + t226 * t74;
t16 = t223 * t74 + t226 * t37;
t240 = -t15 * t223 + t16 * t226;
t25 = pkin(5) * t174 - pkin(11) * t128 + t41;
t33 = pkin(11) * t127 + t42;
t9 = -t229 * t33 + t233 * t25;
t10 = t229 * t25 + t233 * t33;
t52 = t116 * t226 - t223 * t81;
t53 = t116 * t223 + t226 * t81;
t19 = -t229 * t53 + t233 * t52;
t20 = t229 * t52 + t233 * t53;
t183 = t226 * t196;
t115 = -pkin(11) * t267 + t183 + (-pkin(10) * t223 - pkin(5)) * t234;
t131 = -pkin(11) * t271 + t148;
t68 = t115 * t233 - t131 * t229;
t69 = t115 * t229 + t131 * t233;
t197 = t287 * t223;
t199 = t287 * t226;
t132 = -t197 * t233 - t199 * t229;
t133 = -t197 * t229 + t199 * t233;
t238 = t80 * t255 + t279;
t220 = Ifges(5,5) * t255;
t218 = -pkin(5) * t226 - pkin(4);
t210 = Ifges(4,5) * t246;
t205 = Ifges(5,1) * t230 + t283;
t204 = Ifges(5,2) * t234 + t284;
t202 = Ifges(6,1) * t223 + t281;
t201 = Ifges(6,2) * t226 + t282;
t194 = t249 * t230;
t193 = (Ifges(5,1) * t234 - t284) * qJD(4);
t192 = (-Ifges(5,2) * t230 + t283) * qJD(4);
t191 = (mrSges(5,1) * t230 + mrSges(5,2) * t234) * qJD(4);
t190 = -mrSges(6,1) * t234 - mrSges(6,3) * t267;
t189 = mrSges(6,2) * t234 - mrSges(6,3) * t271;
t188 = -mrSges(4,2) * t227 + mrSges(4,3) * t268;
t181 = t249 * t255;
t178 = (mrSges(6,1) * t230 - mrSges(6,3) * t266) * qJD(4);
t177 = (-mrSges(6,2) * t230 - mrSges(6,3) * t270) * qJD(4);
t176 = t316 * t230;
t163 = (mrSges(4,1) * t231 + mrSges(4,2) * t235) * t258;
t157 = -Ifges(6,5) * t234 + t243 * t230;
t156 = -Ifges(6,6) * t234 + t242 * t230;
t155 = -Ifges(6,3) * t234 + t241 * t230;
t152 = t223 * t170;
t147 = -pkin(10) * t270 + t183;
t141 = (Ifges(6,3) * t230 + t241 * t234) * qJD(4);
t137 = -mrSges(7,1) * t234 + mrSges(7,3) * t162;
t136 = mrSges(7,2) * t234 - mrSges(7,3) * t161;
t134 = mrSges(5,2) * t268 - mrSges(5,3) * t174;
t126 = -t226 * t253 + t152;
t118 = mrSges(7,1) * t239 + mrSges(7,2) * t185;
t113 = -Ifges(7,1) * t172 - Ifges(7,4) * t173;
t112 = -Ifges(7,4) * t172 - Ifges(7,2) * t173;
t110 = mrSges(7,1) * t173 - mrSges(7,2) * t172;
t109 = t152 + (-pkin(10) * t267 - pkin(11) * t270) * qJD(4);
t108 = -mrSges(5,2) * t247 - mrSges(5,3) * t130;
t106 = mrSges(7,1) * t161 - mrSges(7,2) * t162;
t103 = Ifges(5,1) * t175 - Ifges(5,4) * t174 - Ifges(5,5) * t268;
t102 = Ifges(5,4) * t175 - Ifges(5,2) * t174 - t252;
t96 = (pkin(5) * t230 - pkin(11) * t266) * qJD(4) + t125;
t95 = -Ifges(7,1) * t162 - Ifges(7,4) * t161 - Ifges(7,5) * t234;
t94 = -Ifges(7,4) * t162 - Ifges(7,2) * t161 - Ifges(7,6) * t234;
t93 = -Ifges(7,5) * t162 - Ifges(7,6) * t161 - Ifges(7,3) * t234;
t91 = -t185 * qJD(5) - t133 * qJD(6);
t90 = -t239 * qJD(5) + t132 * qJD(6);
t87 = mrSges(6,1) * t174 - mrSges(6,3) * t128;
t86 = -mrSges(6,2) * t174 + mrSges(6,3) * t127;
t85 = -mrSges(7,2) * t256 + mrSges(7,3) * t105;
t84 = mrSges(7,1) * t256 - mrSges(7,3) * t104;
t73 = mrSges(5,1) * t130 + mrSges(5,2) * t129;
t67 = Ifges(5,1) * t129 - Ifges(5,4) * t130 + Ifges(5,5) * t247;
t66 = Ifges(5,4) * t129 - Ifges(5,2) * t130 + Ifges(5,6) * t247;
t65 = mrSges(6,1) * t130 - mrSges(6,3) * t98;
t64 = -mrSges(6,2) * t130 + mrSges(6,3) * t97;
t63 = Ifges(6,1) * t128 + Ifges(6,4) * t127 + Ifges(6,5) * t174;
t62 = Ifges(6,4) * t128 + Ifges(6,2) * t127 + Ifges(6,6) * t174;
t61 = Ifges(6,5) * t128 + Ifges(6,6) * t127 + Ifges(6,3) * t174;
t60 = -pkin(5) * t127 + t83;
t59 = mrSges(7,1) * t174 - mrSges(7,3) * t71;
t58 = -mrSges(7,2) * t174 + mrSges(7,3) * t70;
t51 = Ifges(7,1) * t104 + Ifges(7,4) * t105 + Ifges(7,5) * t256;
t50 = Ifges(7,4) * t104 + Ifges(7,2) * t105 + Ifges(7,6) * t256;
t40 = Ifges(6,1) * t98 + Ifges(6,4) * t97 + Ifges(6,5) * t130;
t39 = Ifges(6,4) * t98 + Ifges(6,2) * t97 + Ifges(6,6) * t130;
t38 = Ifges(6,5) * t98 + Ifges(6,6) * t97 + Ifges(6,3) * t130;
t35 = -mrSges(7,1) * t70 + mrSges(7,2) * t71;
t34 = -pkin(5) * t97 + t47;
t32 = Ifges(7,1) * t71 + Ifges(7,4) * t70 + Ifges(7,5) * t174;
t31 = Ifges(7,4) * t71 + Ifges(7,2) * t70 + Ifges(7,6) * t174;
t30 = Ifges(7,5) * t71 + Ifges(7,6) * t70 + Ifges(7,3) * t174;
t29 = -t69 * qJD(6) - t109 * t229 + t233 * t96;
t28 = t68 * qJD(6) + t109 * t233 + t229 * t96;
t18 = -mrSges(7,2) * t130 + mrSges(7,3) * t27;
t17 = mrSges(7,1) * t130 - mrSges(7,3) * t26;
t12 = pkin(11) * t97 + t14;
t11 = pkin(5) * t130 - pkin(11) * t98 + t13;
t7 = Ifges(7,1) * t26 + Ifges(7,4) * t27 + Ifges(7,5) * t130;
t6 = Ifges(7,4) * t26 + Ifges(7,2) * t27 + Ifges(7,6) * t130;
t4 = -t20 * qJD(6) + t15 * t233 - t16 * t229;
t3 = t19 * qJD(6) + t15 * t229 + t16 * t233;
t2 = -t10 * qJD(6) + t11 * t233 - t12 * t229;
t1 = t9 * qJD(6) + t11 * t229 + t12 * t233;
t22 = [0.2e1 * m(7) * (t19 * t4 + t20 * t3 + t21) + 0.2e1 * m(6) * (t15 * t52 + t16 * t53 + t21) + 0.2e1 * m(5) * (t37 * t81 + t21 + t280) + 0.2e1 * m(4) * (t117 * t75 + t171 * t244 + t280); t81 * t108 + t116 * t73 + t37 * t134 + t15 * t87 + t16 * t86 + t171 * t163 + t19 * t17 + t20 * t18 + t75 * t188 + t3 * t58 + t4 * t59 + t52 * t65 + t53 * t64 + t261 * t74 + (-mrSges(3,1) * t232 - mrSges(3,2) * t236) * t259 + (t8 + t276) * t80 + (t35 + t275) * t36 + ((-mrSges(4,1) * t235 + mrSges(4,2) * t231) * t244 + (t116 * t235 - t117 * t231) * qJD(3) * mrSges(4,3)) * t224 + m(4) * (-pkin(2) * t224 ^ 2 * t248 + t117 * t166 - t179 * t74 + t180 * t75 + t272) + m(5) * (t158 * t74 - t36 * t88 + t37 * t89 + t55 * t81 - t56 * t80 + t272) + m(7) * (t1 * t20 + t10 * t3 + t19 * t2 + t34 * t80 + t36 * t60 + t4 * t9) + m(6) * (t13 * t52 + t14 * t53 + t15 * t41 + t16 * t42 + t36 * t83 + t47 * t80); 0.2e1 * t261 * t167 + (t1 * t10 + t2 * t9 + t34 * t60) * t311 + (t13 * t41 + t14 * t42 + t47 * t83) * t312 + (-t235 * t250 - 0.2e1 * pkin(2) * t163 + ((0.2e1 * Ifges(4,4) * t268 + Ifges(4,5) * t227 + t179 * t309) * t235 + (-0.2e1 * Ifges(4,4) * t269 + t180 * t309 + Ifges(5,5) * t175 - 0.2e1 * Ifges(4,6) * t227 - Ifges(5,6) * t174 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3)) * t268) * t231) * qJD(3)) * t224 + (t5 + t38 - t66) * t174 + 0.2e1 * m(5) * (t158 * t167 + t55 * t89 + t56 * t88) + 0.2e1 * m(4) * (t166 * t180 - t167 * t179) + (t30 + t61 - t102) * t130 + t227 * t210 + 0.2e1 * t166 * t188 + t175 * t67 + 0.2e1 * t9 * t17 + 0.2e1 * t10 * t18 + t27 * t31 + t26 * t32 + 0.2e1 * t34 * t35 + 0.2e1 * t1 * t58 + 0.2e1 * t2 * t59 + 0.2e1 * t60 * t8 + 0.2e1 * t42 * t64 + 0.2e1 * t41 * t65 + t70 * t6 + t71 * t7 + 0.2e1 * t47 * t72 + 0.2e1 * t83 * t54 + 0.2e1 * t14 * t86 + 0.2e1 * t13 * t87 + t97 * t62 + t98 * t63 + 0.2e1 * t88 * t107 + 0.2e1 * t89 * t108 + t127 * t39 + t128 * t40 + t129 * t103 + 0.2e1 * t55 * t134 + 0.2e1 * t56 * t135 + 0.2e1 * t158 * t73; -t75 * mrSges(4,2) + t116 * t191 + t3 * t136 + t4 * t137 + t15 * t190 + t16 * t189 + t53 * t177 + t52 * t178 + t19 * t84 + t20 * t85 + t274 * t80 + (t106 + t176) * t36 + m(7) * (t181 * t80 + t19 * t29 + t194 * t36 + t20 * t28 + t3 * t69 + t4 * t68) + m(6) * (t125 * t52 + t126 * t53 + t147 * t15 + t148 * t16) + (t238 * t308 + m(5) * (-t81 * t256 + t238 + t278) / 0.2e1) * t310 + (t279 + t278 + (-t230 * t81 + t234 * t80) * qJD(4)) * mrSges(5,3) + t313 * t74; t313 * t167 + t210 + m(7) * (t1 * t69 + t10 * t28 + t181 * t60 + t194 * t34 + t2 * t68 + t29 * t9) + t7 * t295 + t6 * t296 + t128 * t297 + t127 * t298 + t157 * t302 + t156 * t303 + t51 * t304 + t50 * t305 + (-t192 / 0.2e1 + t49 / 0.2e1 + t141 / 0.2e1) * t174 + ((t63 * t289 + t62 * t290 - t88 * mrSges(5,3) + t103 / 0.2e1) * t234 + (-t89 * mrSges(5,3) + t252 / 0.2e1 - t102 / 0.2e1 + t61 / 0.2e1 + t30 / 0.2e1) * t230) * qJD(4) + (t40 * t289 + t39 * t290 - t56 * mrSges(5,3) + t67 / 0.2e1) * t230 + (t234 * t108 + t276 * t230 + (-t230 * t134 + t275 * t234) * qJD(4) + m(5) * (-t56 * t230 + t55 * t234 - t88 * t255 - t89 * t256) + m(6) * (t230 * t47 + t83 * t255)) * pkin(10) + (-t235 * t220 / 0.2e1 + (Ifges(5,5) * t230 / 0.2e1 + Ifges(5,6) * t234 / 0.2e1 - Ifges(4,6)) * t257) * t224 + (-t204 / 0.2e1 + t93 / 0.2e1 + t155 / 0.2e1) * t130 + (t55 * mrSges(5,3) - t5 / 0.2e1 - t38 / 0.2e1 + t66 / 0.2e1) * t234 + t129 * t205 / 0.2e1 + t14 * t189 + t13 * t190 + t158 * t191 + t175 * t193 / 0.2e1 + t194 * t8 + t47 * t176 + t42 * t177 + t41 * t178 + t181 * t35 + m(6) * (t125 * t41 + t126 * t42 + t13 * t147 + t14 * t148) + t28 * t58 + t29 * t59 + t60 * t57 + t68 * t17 + t69 * t18 - pkin(3) * t73 + t9 * t84 + t10 * t85 + t27 * t94 / 0.2e1 + t26 * t95 / 0.2e1 + t104 * t32 / 0.2e1 + t105 * t31 / 0.2e1 + t34 * t106 + t125 * t87 + t126 * t86 + t1 * t136 + t2 * t137 + t147 * t65 + t148 * t64 + t83 * t164 - t166 * mrSges(4,2); -0.2e1 * pkin(3) * t191 + t104 * t95 + t105 * t94 + 0.2e1 * t181 * t106 + 0.2e1 * t125 * t190 + 0.2e1 * t126 * t189 + 0.2e1 * t28 * t136 + 0.2e1 * t29 * t137 + 0.2e1 * t147 * t178 + 0.2e1 * t148 * t177 - t161 * t50 - t162 * t51 + 0.2e1 * t194 * t57 + 0.2e1 * t68 * t84 + 0.2e1 * t69 * t85 + (t125 * t147 + t126 * t148) * t312 + (t181 * t194 + t28 * t69 + t29 * t68) * t311 + (t192 - t49 - t141) * t234 + (-t223 * t142 + t226 * t143 + t164 * t310 + t193) * t230 + ((t155 - t204 + t93) * t230 + (-t223 * t156 + t226 * t157 + t205 + (t230 * t306 + t176) * t310) * t234) * qJD(4); -t37 * mrSges(5,2) + t80 * t110 + t240 * mrSges(6,3) + (t118 + t273) * t36 + m(7) * (t132 * t4 + t133 * t3 + t19 * t91 + t20 * t90 + t218 * t36) + m(6) * (-pkin(4) * t36 + (-t223 * t52 + t226 * t53) * qJD(5) + t240 * qJ(5)) + (t172 * t19 - t173 * t20 - t185 * t4 - t239 * t3) * mrSges(7,3); t315 * pkin(4) + t250 + t26 * t299 + t27 * t300 + t202 * t302 + t201 * t303 + t113 * t304 + t112 * t305 + t7 * t291 + t6 * t292 + t31 * t293 + t32 * t294 + (m(6) * (-qJ(5) * t13 - qJD(5) * t41) - qJD(5) * t87 - qJ(5) * t65 - t13 * mrSges(6,3) + t40 / 0.2e1) * t223 + (m(6) * (qJ(5) * t14 + qJD(5) * t42) + qJD(5) * t86 + qJ(5) * t64 + t14 * mrSges(6,3) + t39 / 0.2e1) * t226 + t245 * t130 + (-t1 * t239 - t10 * t173 + t172 * t9 - t185 * t2) * mrSges(7,3) + t218 * t8 + t47 * t198 + t174 * t111 / 0.2e1 + m(7) * (t1 * t133 + t10 * t90 + t132 * t2 + t218 * t34 + t9 * t91) - t55 * mrSges(5,2) + t56 * mrSges(5,1) + t90 * t58 + t91 * t59 + t60 * t110 + t34 * t118 + t132 * t17 + t133 * t18; t220 + m(7) * (t132 * t29 + t133 * t28 + t181 * t218 + t68 * t91 + t69 * t90) - t234 * t111 / 0.2e1 + t218 * t57 + t194 * t110 + t181 * t118 + t50 * t292 + t51 * t291 + t95 * t294 + t94 * t293 + t105 * t300 + t104 * t299 + t132 * t84 + t133 * t85 + t90 * t136 + t91 * t137 + t112 * t296 + t113 * t295 - pkin(4) * t164 + (m(6) * (qJ(5) * t126 + qJD(5) * t148) + qJD(5) * t189 + qJ(5) * t177 + t126 * mrSges(6,3) + t298) * t226 + (m(6) * (-qJ(5) * t125 - qJD(5) * t147) - qJD(5) * t190 - qJ(5) * t178 - t125 * mrSges(6,3) + t297) * t223 + (t172 * t68 - t173 * t69 - t185 * t29 - t239 * t28) * mrSges(7,3) + ((pkin(10) * mrSges(5,2) - Ifges(5,6) + t245) * t230 + (t202 * t289 + t201 * t290 + (-m(6) * pkin(4) + t273) * pkin(10)) * t234) * qJD(4); 0.2e1 * t218 * t110 - t172 * t121 + t185 * t113 - t173 * t120 - t239 * t112 + (t132 * t91 + t133 * t90) * t311 + 0.2e1 * (t132 * t172 - t133 * t173 - t185 * t91 - t239 * t90) * mrSges(7,3) + (qJ(5) * t312 + 0.2e1 * mrSges(6,3)) * qJD(5) * (t223 ^ 2 + t226 ^ 2); 0.2e1 * (m(7) / 0.2e1 + t308) * t36; m(7) * t34 - t315 + t8; m(7) * t181 + t255 * t306 + t274; t110; 0; mrSges(7,1) * t4 - mrSges(7,2) * t3; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t29 - mrSges(7,2) * t28 + t49; mrSges(7,1) * t91 - mrSges(7,2) * t90 + t111; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
