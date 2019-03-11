% Calculate time derivative of joint inertia matrix for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:29:03
% EndTime: 2019-03-09 10:29:16
% DurationCPUTime: 5.68s
% Computational Cost: add. (10507->619), mult. (29286->920), div. (0->0), fcn. (29162->12), ass. (0->258)
t219 = sin(pkin(6));
t312 = 0.2e1 * t219;
t217 = sin(pkin(12));
t220 = cos(pkin(12));
t311 = mrSges(6,1) * t217 + mrSges(6,2) * t220;
t218 = sin(pkin(11));
t210 = pkin(2) * t218 + pkin(9);
t310 = 0.2e1 * t210;
t223 = sin(qJ(6));
t226 = cos(qJ(6));
t230 = t217 * t223 - t220 * t226;
t286 = -t230 / 0.2e1;
t185 = t217 * t226 + t220 * t223;
t285 = t185 / 0.2e1;
t283 = t220 / 0.2e1;
t221 = cos(pkin(11));
t225 = sin(qJ(2));
t228 = cos(qJ(2));
t163 = (t218 * t228 + t221 * t225) * t219;
t154 = qJD(2) * t163;
t248 = qJD(2) * t219;
t250 = t221 * t228;
t155 = (-t218 * t225 + t250) * t248;
t242 = t225 * t248;
t234 = pkin(2) * t242;
t102 = pkin(3) * t154 - pkin(9) * t155 + t234;
t254 = t219 * t225;
t162 = t218 * t254 - t219 * t250;
t188 = (-pkin(2) * t228 - pkin(1)) * t219;
t103 = t162 * pkin(3) - t163 * pkin(9) + t188;
t224 = sin(qJ(4));
t227 = cos(qJ(4));
t246 = qJD(4) * t227;
t247 = qJD(4) * t224;
t222 = cos(pkin(6));
t282 = pkin(1) * t222;
t208 = t228 * t282;
t205 = qJD(2) * t208;
t276 = -pkin(8) - qJ(3);
t237 = t276 * t225;
t129 = t205 + (qJD(2) * t237 + qJD(3) * t228) * t219;
t207 = t225 * t282;
t253 = t219 * t228;
t307 = (t253 * t276 - t207) * qJD(2) - qJD(3) * t254;
t83 = t221 * t129 + t218 * t307;
t141 = pkin(2) * t222 + t219 * t237 + t208;
t181 = pkin(8) * t253 + t207;
t153 = qJ(3) * t253 + t181;
t98 = t218 * t141 + t221 * t153;
t88 = pkin(9) * t222 + t98;
t32 = t102 * t227 - t103 * t247 - t224 * t83 - t88 * t246;
t25 = -pkin(4) * t154 - t32;
t137 = t163 * t224 - t222 * t227;
t96 = -qJD(4) * t137 + t155 * t227;
t74 = t154 * t220 - t217 * t96;
t75 = t154 * t217 + t220 * t96;
t43 = -t74 * mrSges(6,1) + t75 * mrSges(6,2);
t308 = -m(6) * t25 - t43;
t249 = t217 ^ 2 + t220 ^ 2;
t175 = t230 * qJD(6);
t306 = 0.2e1 * m(6);
t305 = 2 * m(7);
t304 = -2 * mrSges(3,3);
t303 = -2 * mrSges(4,3);
t302 = -2 * Ifges(4,4);
t301 = 0.2e1 * t188;
t300 = m(4) * pkin(2);
t138 = t163 * t227 + t222 * t224;
t100 = -t138 * t217 + t162 * t220;
t101 = t138 * t220 + t162 * t217;
t61 = t100 * t226 - t101 * t223;
t299 = t61 / 0.2e1;
t62 = t100 * t223 + t101 * t226;
t298 = t62 / 0.2e1;
t297 = t74 / 0.2e1;
t296 = t75 / 0.2e1;
t134 = Ifges(7,4) * t185 - Ifges(7,2) * t230;
t294 = t134 / 0.2e1;
t135 = Ifges(7,1) * t185 - Ifges(7,4) * t230;
t293 = t135 / 0.2e1;
t268 = Ifges(6,4) * t220;
t232 = -Ifges(6,2) * t217 + t268;
t147 = (Ifges(6,6) * t224 + t227 * t232) * qJD(4);
t292 = t147 / 0.2e1;
t269 = Ifges(6,4) * t217;
t233 = Ifges(6,1) * t220 - t269;
t148 = (Ifges(6,5) * t224 + t227 * t233) * qJD(4);
t291 = t148 / 0.2e1;
t164 = t185 * t224;
t290 = -t164 / 0.2e1;
t165 = t230 * t224;
t289 = -t165 / 0.2e1;
t288 = -t175 / 0.2e1;
t176 = t185 * qJD(6);
t287 = -t176 / 0.2e1;
t284 = -t217 / 0.2e1;
t281 = pkin(4) * t224;
t82 = t129 * t218 - t221 * t307;
t280 = t82 * mrSges(4,1);
t279 = t82 * mrSges(5,1);
t278 = t82 * mrSges(5,2);
t277 = t83 * mrSges(4,2);
t275 = pkin(10) + qJ(5);
t31 = t224 * t102 + t103 * t246 + t227 * t83 - t247 * t88;
t19 = qJ(5) * t154 + qJD(5) * t162 + t31;
t95 = qJD(4) * t138 + t155 * t224;
t39 = pkin(4) * t95 - qJ(5) * t96 - qJD(5) * t138 + t82;
t12 = t220 * t19 + t217 * t39;
t60 = t224 * t103 + t227 * t88;
t52 = qJ(5) * t162 + t60;
t97 = t141 * t221 - t218 * t153;
t87 = -pkin(3) * t222 - t97;
t58 = pkin(4) * t137 - qJ(5) * t138 + t87;
t27 = t217 * t58 + t220 * t52;
t79 = mrSges(5,1) * t154 - mrSges(5,3) * t96;
t274 = -t79 + t43;
t271 = Ifges(5,4) * t224;
t270 = Ifges(5,4) * t227;
t267 = t154 * Ifges(5,5);
t266 = t154 * Ifges(5,6);
t265 = t162 * Ifges(5,6);
t169 = -pkin(8) * t242 + t205;
t264 = t169 * mrSges(3,2);
t170 = t181 * qJD(2);
t263 = t170 * mrSges(3,1);
t105 = mrSges(5,1) * t162 - mrSges(5,3) * t138;
t64 = -mrSges(6,1) * t100 + mrSges(6,2) * t101;
t262 = -t105 + t64;
t168 = t311 * t246;
t118 = -t176 * t224 - t230 * t246;
t119 = t175 * t224 - t185 * t246;
t73 = -t119 * mrSges(7,1) + t118 * mrSges(7,2);
t261 = -t168 - t73;
t194 = -mrSges(6,1) * t220 + mrSges(6,2) * t217;
t260 = t194 - mrSges(5,1);
t259 = qJ(5) * t227;
t258 = t210 * t224;
t257 = t210 * t227;
t256 = t217 * t224;
t255 = t217 * t227;
t252 = t220 * t224;
t251 = t220 * t227;
t124 = -Ifges(7,5) * t175 - Ifges(7,6) * t176;
t245 = qJD(5) * t224;
t174 = -t245 + (-t259 + t281) * qJD(4);
t241 = t210 * t247;
t127 = t220 * t174 + t217 * t241;
t211 = -pkin(2) * t221 - pkin(3);
t182 = -pkin(4) * t227 - qJ(5) * t224 + t211;
t131 = t217 * t182 + t210 * t251;
t23 = qJD(6) * t61 + t223 * t74 + t226 * t75;
t24 = -qJD(6) * t62 - t223 * t75 + t226 * t74;
t6 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t95;
t244 = Ifges(5,5) * t96 - Ifges(5,6) * t95 + Ifges(5,3) * t154;
t65 = Ifges(7,5) * t118 + Ifges(7,6) * t119 + Ifges(7,3) * t247;
t243 = Ifges(3,5) * t228 * t248 + Ifges(4,5) * t155 - Ifges(4,6) * t154;
t239 = Ifges(6,5) * t217 / 0.2e1 + Ifges(6,6) * t283 + Ifges(7,5) * t285 + Ifges(7,6) * t286;
t10 = -t24 * mrSges(7,1) + t23 * mrSges(7,2);
t238 = m(6) * t249;
t236 = pkin(5) * t217 + t210;
t11 = -t19 * t217 + t220 * t39;
t26 = -t217 * t52 + t220 * t58;
t235 = t249 * mrSges(6,3);
t59 = t103 * t227 - t224 * t88;
t231 = Ifges(6,5) * t220 - Ifges(6,6) * t217;
t16 = pkin(5) * t137 - pkin(10) * t101 + t26;
t17 = pkin(10) * t100 + t27;
t3 = t16 * t226 - t17 * t223;
t4 = t16 * t223 + t17 * t226;
t167 = t220 * t182;
t117 = -pkin(10) * t252 + t167 + (-t210 * t217 - pkin(5)) * t227;
t121 = -pkin(10) * t256 + t131;
t71 = t117 * t226 - t121 * t223;
t72 = t117 * t223 + t121 * t226;
t193 = t275 * t217;
t195 = t275 * t220;
t139 = -t193 * t226 - t195 * t223;
t140 = -t193 * t223 + t195 * t226;
t53 = -pkin(4) * t162 - t59;
t214 = Ifges(5,5) * t246;
t212 = -pkin(5) * t220 - pkin(4);
t200 = Ifges(5,1) * t224 + t270;
t199 = Ifges(5,2) * t227 + t271;
t198 = Ifges(6,1) * t217 + t268;
t197 = Ifges(6,2) * t220 + t269;
t192 = (Ifges(5,1) * t227 - t271) * qJD(4);
t191 = (-Ifges(5,2) * t224 + t270) * qJD(4);
t190 = (t224 * mrSges(5,1) + t227 * mrSges(5,2)) * qJD(4);
t187 = -mrSges(6,1) * t227 - mrSges(6,3) * t252;
t186 = mrSges(6,2) * t227 - mrSges(6,3) * t256;
t180 = -pkin(8) * t254 + t208;
t179 = (mrSges(6,1) * t224 - mrSges(6,3) * t251) * qJD(4);
t178 = (-mrSges(6,2) * t224 - mrSges(6,3) * t255) * qJD(4);
t177 = t311 * t224;
t173 = t236 * t224;
t161 = -Ifges(6,5) * t227 + t224 * t233;
t160 = -Ifges(6,6) * t227 + t224 * t232;
t159 = -Ifges(6,3) * t227 + t224 * t231;
t158 = t236 * t246;
t156 = t217 * t174;
t149 = t155 * mrSges(4,2);
t146 = (Ifges(6,3) * t224 + t227 * t231) * qJD(4);
t145 = -mrSges(7,1) * t227 + mrSges(7,3) * t165;
t144 = mrSges(7,2) * t227 - mrSges(7,3) * t164;
t132 = mrSges(7,1) * t230 + mrSges(7,2) * t185;
t130 = -t210 * t255 + t167;
t128 = -t220 * t241 + t156;
t126 = -Ifges(7,1) * t175 - Ifges(7,4) * t176;
t125 = -Ifges(7,4) * t175 - Ifges(7,2) * t176;
t123 = mrSges(7,1) * t176 - mrSges(7,2) * t175;
t120 = mrSges(7,1) * t164 - mrSges(7,2) * t165;
t114 = t156 + (-pkin(10) * t255 - t210 * t252) * qJD(4);
t113 = -Ifges(7,1) * t165 - Ifges(7,4) * t164 - Ifges(7,5) * t227;
t112 = -Ifges(7,4) * t165 - Ifges(7,2) * t164 - Ifges(7,6) * t227;
t111 = -Ifges(7,5) * t165 - Ifges(7,6) * t164 - Ifges(7,3) * t227;
t110 = -qJD(5) * t185 - qJD(6) * t140;
t109 = -qJD(5) * t230 + qJD(6) * t139;
t108 = (pkin(5) * t224 - pkin(10) * t251) * qJD(4) + t127;
t107 = -mrSges(7,2) * t247 + mrSges(7,3) * t119;
t106 = mrSges(7,1) * t247 - mrSges(7,3) * t118;
t104 = -mrSges(5,2) * t162 - mrSges(5,3) * t137;
t78 = -mrSges(5,2) * t154 - mrSges(5,3) * t95;
t77 = Ifges(5,1) * t138 - Ifges(5,4) * t137 + Ifges(5,5) * t162;
t76 = Ifges(5,4) * t138 - Ifges(5,2) * t137 + t265;
t69 = mrSges(6,1) * t137 - mrSges(6,3) * t101;
t68 = -mrSges(6,2) * t137 + mrSges(6,3) * t100;
t67 = Ifges(7,1) * t118 + Ifges(7,4) * t119 + Ifges(7,5) * t247;
t66 = Ifges(7,4) * t118 + Ifges(7,2) * t119 + Ifges(7,6) * t247;
t63 = mrSges(5,1) * t95 + mrSges(5,2) * t96;
t55 = Ifges(5,1) * t96 - Ifges(5,4) * t95 + t267;
t54 = Ifges(5,4) * t96 - Ifges(5,2) * t95 + t266;
t50 = Ifges(6,1) * t101 + Ifges(6,4) * t100 + Ifges(6,5) * t137;
t49 = Ifges(6,4) * t101 + Ifges(6,2) * t100 + Ifges(6,6) * t137;
t48 = Ifges(6,5) * t101 + Ifges(6,6) * t100 + Ifges(6,3) * t137;
t47 = mrSges(6,1) * t95 - mrSges(6,3) * t75;
t46 = -mrSges(6,2) * t95 + mrSges(6,3) * t74;
t45 = mrSges(7,1) * t137 - mrSges(7,3) * t62;
t44 = -mrSges(7,2) * t137 + mrSges(7,3) * t61;
t42 = -qJD(6) * t72 + t108 * t226 - t114 * t223;
t41 = qJD(6) * t71 + t108 * t223 + t114 * t226;
t40 = -pkin(5) * t100 + t53;
t36 = Ifges(6,1) * t75 + Ifges(6,4) * t74 + Ifges(6,5) * t95;
t35 = Ifges(6,4) * t75 + Ifges(6,2) * t74 + Ifges(6,6) * t95;
t34 = Ifges(6,5) * t75 + Ifges(6,6) * t74 + Ifges(6,3) * t95;
t33 = -mrSges(7,1) * t61 + mrSges(7,2) * t62;
t30 = Ifges(7,1) * t62 + Ifges(7,4) * t61 + Ifges(7,5) * t137;
t29 = Ifges(7,4) * t62 + Ifges(7,2) * t61 + Ifges(7,6) * t137;
t28 = Ifges(7,5) * t62 + Ifges(7,6) * t61 + Ifges(7,3) * t137;
t15 = -mrSges(7,2) * t95 + mrSges(7,3) * t24;
t14 = mrSges(7,1) * t95 - mrSges(7,3) * t23;
t13 = -pkin(5) * t74 + t25;
t9 = pkin(10) * t74 + t12;
t8 = Ifges(7,1) * t23 + Ifges(7,4) * t24 + Ifges(7,5) * t95;
t7 = Ifges(7,4) * t23 + Ifges(7,2) * t24 + Ifges(7,6) * t95;
t5 = pkin(5) * t95 - pkin(10) * t75 + t11;
t2 = -qJD(6) * t4 - t223 * t9 + t226 * t5;
t1 = qJD(6) * t3 + t223 * t5 + t226 * t9;
t18 = [0.2e1 * m(3) * (t169 * t181 - t170 * t180) + (t243 - 0.2e1 * t263 - 0.2e1 * t264 - 0.2e1 * t277 - 0.2e1 * t280) * t222 + 0.2e1 * m(4) * (-t82 * t97 + t83 * t98) + t100 * t35 + t101 * t36 + 0.2e1 * t31 * t104 + 0.2e1 * t32 * t105 + t96 * t77 + 0.2e1 * t87 * t63 + t74 * t49 + t75 * t50 + 0.2e1 * t60 * t78 + 0.2e1 * t59 * t79 + 0.2e1 * t12 * t68 + 0.2e1 * t11 * t69 + t61 * t7 + t62 * t8 + 0.2e1 * t25 * t64 + 0.2e1 * t53 * t43 + 0.2e1 * t1 * t44 + 0.2e1 * t2 * t45 + 0.2e1 * t27 * t46 + 0.2e1 * t26 * t47 + 0.2e1 * t40 * t10 + 0.2e1 * t13 * t33 + t24 * t29 + t23 * t30 + 0.2e1 * t3 * t14 + 0.2e1 * t4 * t15 + t149 * t301 + (t1 * t4 + t13 * t40 + t2 * t3) * t305 + (t11 * t26 + t12 * t27 + t25 * t53) * t306 + 0.2e1 * m(5) * (t31 * t60 + t32 * t59 + t82 * t87) + (t48 + t28 - t76) * t95 + (mrSges(4,1) * t301 + t98 * t303 + t163 * t302 + Ifges(5,5) * t138 - Ifges(4,6) * t222 - Ifges(5,6) * t137 + ((2 * Ifges(4,2)) + Ifges(5,3)) * t162) * t154 + (0.2e1 * Ifges(4,1) * t163 + Ifges(4,5) * t222 + t162 * t302 + t303 * t97) * t155 + (t55 + 0.2e1 * t278) * t138 + (t34 - t54 + t6 + 0.2e1 * t279) * t137 + t162 * t244 + (0.2e1 * (t169 * t228 + t170 * t225) * mrSges(3,3) + ((t180 * t304 + Ifges(3,5) * t222 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t228) * t312) * t228 + (t300 * t301 - 0.2e1 * Ifges(3,6) * t222 + 0.2e1 * pkin(2) * (mrSges(4,1) * t162 + mrSges(4,2) * t163) + t181 * t304 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t225 + (Ifges(3,1) - Ifges(3,2)) * t228) * t312) * t225) * qJD(2)) * t219 + 0.2e1 * (-t162 * t83 + t163 * t82) * mrSges(4,3); t243 + t211 * t63 + t96 * t200 / 0.2e1 + t12 * t186 + t11 * t187 + t87 * t190 + t138 * t192 / 0.2e1 + t25 * t177 + t27 * t178 + t26 * t179 + t53 * t168 + t173 * t10 + t158 * t33 + t1 * t144 + t2 * t145 + t130 * t47 + t131 * t46 + t127 * t69 + t128 * t68 + t118 * t30 / 0.2e1 + t119 * t29 / 0.2e1 + t13 * t120 + t3 * t106 + t4 * t107 + t24 * t112 / 0.2e1 + t23 * t113 / 0.2e1 + t72 * t15 + t40 * t73 + t71 * t14 + t41 * t44 + t42 * t45 + t162 * t214 / 0.2e1 + t161 * t296 + t160 * t297 + t67 * t298 + t66 * t299 + (t218 * t83 - t221 * t82) * t300 + t8 * t289 + t7 * t290 + t101 * t291 + t100 * t292 - t280 + (-t191 / 0.2e1 + t146 / 0.2e1 + t65 / 0.2e1) * t137 + (-t154 * t218 - t155 * t221) * pkin(2) * mrSges(4,3) + ((-t76 / 0.2e1 + t48 / 0.2e1 + t28 / 0.2e1 - t265 / 0.2e1 - t60 * mrSges(5,3) + (-m(5) * t60 - t104) * t210) * t224 + (t50 * t283 + t49 * t284 - t59 * mrSges(5,3) + t77 / 0.2e1 + (-m(5) * t59 + m(6) * t53 + t262) * t210) * t227) * qJD(4) + m(7) * (t1 * t72 + t13 * t173 + t158 * t40 + t2 * t71 + t3 * t42 + t4 * t41) + (-t199 / 0.2e1 + t159 / 0.2e1 + t111 / 0.2e1) * t95 + (-t32 * mrSges(5,3) + t36 * t283 + t35 * t284 + t267 / 0.2e1 + t278 + t55 / 0.2e1 + t274 * t210) * t224 + (t266 / 0.2e1 - t279 + t54 / 0.2e1 - t34 / 0.2e1 - t6 / 0.2e1 + t210 * t78 + t31 * mrSges(5,3)) * t227 + m(5) * (t211 * t82 + t257 * t31 - t258 * t32) + m(6) * (t11 * t130 + t12 * t131 + t127 * t26 + t128 * t27 + t25 * t258) - t263 - t264 - t277 - Ifges(3,6) * t242; 0.2e1 * t71 * t106 + 0.2e1 * t72 * t107 + t119 * t112 + t118 * t113 + 0.2e1 * t158 * t120 + 0.2e1 * t127 * t187 + 0.2e1 * t128 * t186 + 0.2e1 * t130 * t179 + 0.2e1 * t131 * t178 + 0.2e1 * t41 * t144 + 0.2e1 * t42 * t145 - t164 * t66 - t165 * t67 + 0.2e1 * t173 * t73 + 0.2e1 * t211 * t190 + (t158 * t173 + t41 * t72 + t42 * t71) * t305 + (t127 * t130 + t128 * t131) * t306 + (t191 - t146 - t65) * t227 + (-t217 * t147 + t220 * t148 + t168 * t310 + t192) * t224 + ((t111 + t159 - t199) * t224 + (-t217 * t160 + t220 * t161 + t200 + (m(6) * t258 + t177) * t310) * t227) * qJD(4); m(4) * t234 + t154 * mrSges(4,1) + t118 * t44 + t119 * t45 - t164 * t14 - t165 * t15 + t149 + (-t10 - t274) * t227 + (-t217 * t47 + t220 * t46 + t78) * t224 + ((-t217 * t69 + t220 * t68 + t104) * t227 + (t33 + t262) * t224) * qJD(4) + m(7) * (-t1 * t165 + t118 * t4 + t119 * t3 - t13 * t227 - t164 * t2 + t247 * t40) + m(6) * (-t227 * t25 + (-t11 * t217 + t12 * t220) * t224 + (t224 * t53 + (-t217 * t26 + t220 * t27) * t227) * qJD(4)) + m(5) * (t224 * t31 + t227 * t32 + (-t224 * t59 + t227 * t60) * qJD(4)); -t164 * t106 - t165 * t107 + t118 * t144 + t119 * t145 + t261 * t227 + (t220 * t178 - t217 * t179) * t224 + ((t186 * t220 - t187 * t217) * t227 + (t120 + t177) * t224) * qJD(4) + m(7) * (t118 * t72 + t119 * t71 - t158 * t227 - t164 * t42 - t165 * t41 + t173 * t247) + m(6) * ((-t127 * t217 + t128 * t220) * t224 + (t210 * t224 ^ 2 + (-t130 * t217 + t131 * t220 - t257) * t227) * qJD(4)); (-t118 * t165 - t119 * t164) * t305 + 0.4e1 * (m(6) * (-0.1e1 + t249) / 0.2e1 - m(7) / 0.2e1) * t224 * t246; (m(6) * (qJ(5) * t12 + qJD(5) * t27) + qJD(5) * t68 + qJ(5) * t46 + t12 * mrSges(6,3) + t35 / 0.2e1) * t220 + (m(6) * (-qJ(5) * t11 - qJD(5) * t26) - qJD(5) * t69 - qJ(5) * t47 - t11 * mrSges(6,3) + t36 / 0.2e1) * t217 + (-t1 * t230 + t175 * t3 - t176 * t4 - t185 * t2) * mrSges(7,3) + t212 * t10 + t25 * t194 + t137 * t124 / 0.2e1 + t139 * t14 + t140 * t15 + t13 * t132 + t40 * t123 + t109 * t44 + t110 * t45 - t31 * mrSges(5,2) + t32 * mrSges(5,1) + t244 + t198 * t296 + t197 * t297 + t126 * t298 + t125 * t299 + t8 * t285 + t7 * t286 + t29 * t287 + t30 * t288 + t23 * t293 + t24 * t294 + m(7) * (t1 * t140 + t109 * t4 + t110 * t3 + t13 * t212 + t139 * t2) + t308 * pkin(4) + t239 * t95; m(7) * (t109 * t72 + t110 * t71 + t139 * t42 + t140 * t41 + t158 * t212) - t227 * t124 / 0.2e1 + t212 * t73 + t66 * t286 + t67 * t285 + t113 * t288 + t112 * t287 + t126 * t289 - pkin(4) * t168 + t173 * t123 + t158 * t132 + t125 * t290 + t109 * t144 + t110 * t145 + t139 * t106 + t140 * t107 + t119 * t294 + t118 * t293 + t214 + (qJD(5) * t186 + qJ(5) * t178 + t128 * mrSges(6,3) + t292 + m(6) * (qJ(5) * t128 + qJD(5) * t131)) * t220 + (-qJD(5) * t187 - qJ(5) * t179 - t127 * mrSges(6,3) + t291 + m(6) * (-qJ(5) * t127 - qJD(5) * t130)) * t217 + (t175 * t71 - t176 * t72 - t185 * t42 - t230 * t41) * mrSges(7,3) + ((t210 * mrSges(5,2) - Ifges(5,6) + t239) * t224 + (t198 * t283 + t197 * t284 + (-m(6) * pkin(4) + t260) * t210) * t227) * qJD(4); -t227 * t123 + m(7) * (-t109 * t165 - t110 * t164 + t118 * t140 + t119 * t139) + t238 * t245 + ((-mrSges(5,2) + t235) * t227 + m(6) * (t249 * t259 - t281) + (m(7) * t212 + t132 + t260) * t224) * qJD(4) + (-t118 * t230 - t119 * t185 - t164 * t175 + t165 * t176) * mrSges(7,3); (t109 * t140 + t110 * t139) * t305 - t176 * t134 - t230 * t125 + 0.2e1 * t212 * t123 - t175 * t135 + t185 * t126 + 0.2e1 * (-t109 * t230 - t110 * t185 + t139 * t175 - t140 * t176) * mrSges(7,3) + 0.2e1 * (qJ(5) * t238 + t235) * qJD(5); m(7) * t13 + t10 - t308; m(6) * t210 * t246 + m(7) * t158 - t261; (m(6) + m(7)) * t247; t123; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t6; mrSges(7,1) * t42 - mrSges(7,2) * t41 + t65; -t73; mrSges(7,1) * t110 - mrSges(7,2) * t109 + t124; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
