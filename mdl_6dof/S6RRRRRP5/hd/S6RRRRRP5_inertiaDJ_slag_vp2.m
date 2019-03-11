% Calculate time derivative of joint inertia matrix for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:13
% EndTime: 2019-03-10 01:18:22
% DurationCPUTime: 5.49s
% Computational Cost: add. (10898->566), mult. (26139->814), div. (0->0), fcn. (23694->8), ass. (0->234)
t294 = Ifges(6,5) + Ifges(7,5);
t293 = Ifges(6,6) + Ifges(7,6);
t296 = Ifges(6,3) + Ifges(7,3);
t218 = sin(qJ(3));
t219 = sin(qJ(2));
t222 = cos(qJ(3));
t255 = qJD(3) * t222;
t223 = cos(qJ(2));
t258 = qJD(2) * t223;
t230 = t218 * t258 + t219 * t255;
t216 = sin(qJ(5));
t295 = pkin(4) * t216;
t220 = cos(qJ(5));
t280 = pkin(4) * t220;
t277 = -mrSges(6,1) - mrSges(7,1);
t217 = sin(qJ(4));
t221 = cos(qJ(4));
t235 = t217 * t218 - t221 * t222;
t166 = t235 * t219;
t282 = -pkin(9) - pkin(8);
t196 = t282 * t218;
t197 = t282 * t222;
t152 = t221 * t196 + t197 * t217;
t180 = t217 * t222 + t218 * t221;
t128 = -pkin(10) * t180 + t152;
t153 = t217 * t196 - t221 * t197;
t129 = -pkin(10) * t235 + t153;
t79 = t216 * t128 + t220 * t129;
t193 = -pkin(2) * t223 - t219 * pkin(8) - pkin(1);
t178 = t222 * t193;
t266 = t219 * t222;
t279 = pkin(7) * t218;
t141 = -pkin(9) * t266 + t178 + (-pkin(3) - t279) * t223;
t264 = t222 * t223;
t204 = pkin(7) * t264;
t161 = t218 * t193 + t204;
t267 = t218 * t219;
t151 = -pkin(9) * t267 + t161;
t95 = t217 * t141 + t221 * t151;
t245 = t222 * t258;
t259 = qJD(2) * t219;
t292 = -Ifges(4,5) * t245 - Ifges(4,3) * t259;
t291 = qJD(3) + qJD(4);
t147 = t291 * t180;
t103 = -t147 * t219 - t235 * t258;
t104 = t166 * t291 - t180 * t258;
t290 = -Ifges(5,5) * t103 - Ifges(5,6) * t104 - Ifges(5,3) * t259;
t289 = 2 * m(4);
t288 = 2 * m(5);
t287 = 2 * m(6);
t286 = 2 * m(7);
t285 = -0.2e1 * pkin(1);
t284 = 0.2e1 * pkin(7);
t283 = m(7) * pkin(5);
t281 = -t218 / 0.2e1;
t276 = -mrSges(7,2) - mrSges(6,2);
t165 = t180 * t219;
t120 = -t165 * t216 - t166 * t220;
t48 = -qJD(5) * t120 - t103 * t216 + t104 * t220;
t40 = -mrSges(7,2) * t259 + mrSges(7,3) * t48;
t41 = -mrSges(6,2) * t259 + mrSges(6,3) * t48;
t275 = t40 + t41;
t94 = t221 * t141 - t217 * t151;
t74 = -pkin(4) * t223 + t166 * pkin(10) + t94;
t80 = -pkin(10) * t165 + t95;
t37 = t216 * t74 + t220 * t80;
t274 = Ifges(4,4) * t218;
t273 = Ifges(4,4) * t222;
t272 = Ifges(4,6) * t218;
t271 = pkin(4) * qJD(5);
t270 = t223 * Ifges(4,6);
t269 = t216 * t217;
t268 = t217 * t220;
t195 = Ifges(4,1) * t218 + t273;
t265 = t222 * t195;
t119 = -t165 * t220 + t166 * t216;
t110 = mrSges(7,2) * t223 + t119 * mrSges(7,3);
t111 = mrSges(6,2) * t223 + t119 * mrSges(6,3);
t263 = t110 + t111;
t112 = -mrSges(7,1) * t223 - t120 * mrSges(7,3);
t113 = -mrSges(6,1) * t223 - t120 * mrSges(6,3);
t262 = t112 + t113;
t206 = pkin(3) * t221 + pkin(4);
t251 = qJD(5) * t220;
t252 = qJD(5) * t216;
t134 = t206 * t251 + (-t217 * t252 + (t220 * t221 - t269) * qJD(4)) * pkin(3);
t169 = pkin(3) * t268 + t206 * t216;
t261 = pkin(4) * t169 * t251 + t134 * t295;
t191 = (pkin(2) * t219 - pkin(8) * t223) * qJD(2);
t260 = t222 * t191 + t259 * t279;
t192 = pkin(3) * t267 + t219 * pkin(7);
t257 = qJD(3) * t218;
t256 = qJD(3) * t219;
t254 = qJD(4) * t217;
t253 = qJD(4) * t221;
t214 = pkin(7) * t258;
t213 = pkin(3) * t257;
t47 = qJD(5) * t119 + t103 * t220 + t104 * t216;
t114 = t218 * t191 + t193 * t255 + (-t222 * t259 - t223 * t257) * pkin(7);
t100 = -pkin(9) * t230 + t114;
t91 = (pkin(3) * t219 - pkin(9) * t264) * qJD(2) + (-t204 + (pkin(9) * t219 - t193) * t218) * qJD(3) + t260;
t29 = -qJD(4) * t95 - t100 * t217 + t221 * t91;
t21 = pkin(4) * t259 - pkin(10) * t103 + t29;
t28 = t221 * t100 + t141 * t253 - t151 * t254 + t217 * t91;
t23 = pkin(10) * t104 + t28;
t6 = -qJD(5) * t37 + t220 * t21 - t216 * t23;
t2 = pkin(5) * t259 - qJ(6) * t47 - qJD(6) * t120 + t6;
t38 = mrSges(7,1) * t259 - mrSges(7,3) * t47;
t250 = m(7) * t2 + t38;
t158 = pkin(3) * t230 + t214;
t207 = -pkin(3) * t222 - pkin(2);
t249 = qJD(3) * t282;
t247 = t218 * t256;
t14 = -t48 * mrSges(7,1) + t47 * mrSges(7,2);
t139 = -t180 * t216 - t220 * t235;
t146 = t291 * t235;
t60 = qJD(5) * t139 - t146 * t220 - t147 * t216;
t140 = t180 * t220 - t216 * t235;
t61 = -qJD(5) * t140 + t146 * t216 - t147 * t220;
t30 = -t61 * mrSges(7,1) + t60 * mrSges(7,2);
t244 = t276 * t220;
t243 = (2 * Ifges(3,4)) + t272;
t130 = pkin(4) * t147 + t213;
t36 = -t216 * t80 + t220 * t74;
t78 = t220 * t128 - t129 * t216;
t242 = 0.2e1 * t276;
t189 = t218 * t249;
t190 = t222 * t249;
t106 = t221 * t189 + t217 * t190 + t196 * t253 + t197 * t254;
t72 = -pkin(10) * t147 + t106;
t107 = -qJD(4) * t153 - t189 * t217 + t221 * t190;
t73 = pkin(10) * t146 + t107;
t20 = -qJD(5) * t79 - t216 * t72 + t220 * t73;
t9 = -qJ(6) * t60 - qJD(6) * t140 + t20;
t241 = m(7) * t9 - t60 * mrSges(7,3);
t168 = -pkin(3) * t269 + t220 * t206;
t144 = pkin(4) * t165 + t192;
t240 = -mrSges(4,1) * t222 + mrSges(4,2) * t218;
t239 = mrSges(4,1) * t218 + mrSges(4,2) * t222;
t238 = Ifges(4,1) * t222 - t274;
t237 = -Ifges(4,2) * t218 + t273;
t194 = Ifges(4,2) * t222 + t274;
t236 = Ifges(4,5) * t218 + Ifges(4,6) * t222;
t162 = pkin(4) * t235 + t207;
t81 = -pkin(4) * t104 + t158;
t234 = t296 * t259 + t293 * t48 + t294 * t47;
t135 = -t206 * t252 + (-t217 * t251 + (-t216 * t221 - t268) * qJD(4)) * pkin(3);
t131 = t135 * mrSges(7,1);
t132 = t135 * mrSges(6,1);
t233 = t134 * t276 + t131 + t132;
t5 = t216 * t21 + t220 * t23 + t74 * t251 - t252 * t80;
t19 = t128 * t251 - t129 * t252 + t216 * t73 + t220 * t72;
t232 = (-mrSges(5,1) * t217 - mrSges(5,2) * t221) * qJD(4) * pkin(3);
t231 = t245 - t247;
t229 = t134 * t139 - t135 * t140 + t169 * t61;
t56 = Ifges(7,6) * t61;
t57 = Ifges(6,6) * t61;
t58 = Ifges(7,5) * t60;
t59 = Ifges(6,5) * t60;
t8 = qJ(6) * t61 + qJD(6) * t139 + t19;
t228 = t20 * mrSges(6,1) + t9 * mrSges(7,1) - t19 * mrSges(6,2) - t8 * mrSges(7,2) + t56 + t57 + t58 + t59;
t227 = t216 * t61 + (t139 * t220 + t140 * t216) * qJD(5);
t3 = qJ(6) * t48 + qJD(6) * t119 + t5;
t226 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) + t234;
t142 = Ifges(5,6) * t147;
t143 = Ifges(5,5) * t146;
t225 = t107 * mrSges(5,1) - t106 * mrSges(5,2) - t142 - t143 + t228;
t224 = t29 * mrSges(5,1) - t28 * mrSges(5,2) + t226 - t290;
t212 = Ifges(4,5) * t255;
t205 = pkin(5) + t280;
t188 = -mrSges(4,1) * t223 - mrSges(4,3) * t266;
t187 = mrSges(4,2) * t223 - mrSges(4,3) * t267;
t186 = t238 * qJD(3);
t185 = t237 * qJD(3);
t184 = t239 * qJD(3);
t167 = pkin(5) + t168;
t164 = -Ifges(4,5) * t223 + t219 * t238;
t163 = t219 * t237 - t270;
t160 = -t223 * t279 + t178;
t157 = -mrSges(4,2) * t259 - mrSges(4,3) * t230;
t156 = mrSges(4,1) * t259 - mrSges(4,3) * t231;
t155 = -mrSges(5,1) * t223 + t166 * mrSges(5,3);
t154 = mrSges(5,2) * t223 - t165 * mrSges(5,3);
t150 = Ifges(5,1) * t180 - Ifges(5,4) * t235;
t149 = Ifges(5,4) * t180 - Ifges(5,2) * t235;
t148 = mrSges(5,1) * t235 + mrSges(5,2) * t180;
t138 = mrSges(4,1) * t230 + mrSges(4,2) * t231;
t127 = mrSges(5,1) * t165 - mrSges(5,2) * t166;
t125 = -t195 * t256 + (Ifges(4,5) * t219 + t223 * t238) * qJD(2);
t124 = -t194 * t256 + (Ifges(4,6) * t219 + t223 * t237) * qJD(2);
t118 = -Ifges(5,1) * t166 - Ifges(5,4) * t165 - Ifges(5,5) * t223;
t117 = -Ifges(5,4) * t166 - Ifges(5,2) * t165 - Ifges(5,6) * t223;
t115 = -qJD(3) * t161 + t260;
t109 = t169 * t134;
t108 = -pkin(5) * t139 + t162;
t99 = -Ifges(5,1) * t146 - Ifges(5,4) * t147;
t98 = -Ifges(5,4) * t146 - Ifges(5,2) * t147;
t97 = mrSges(5,1) * t147 - mrSges(5,2) * t146;
t93 = -mrSges(5,2) * t259 + mrSges(5,3) * t104;
t92 = mrSges(5,1) * t259 - mrSges(5,3) * t103;
t90 = Ifges(6,1) * t140 + Ifges(6,4) * t139;
t89 = Ifges(7,1) * t140 + Ifges(7,4) * t139;
t88 = Ifges(6,4) * t140 + Ifges(6,2) * t139;
t87 = Ifges(7,4) * t140 + Ifges(7,2) * t139;
t86 = -mrSges(6,1) * t139 + mrSges(6,2) * t140;
t85 = -mrSges(7,1) * t139 + mrSges(7,2) * t140;
t82 = -pkin(5) * t119 + t144;
t76 = -mrSges(6,1) * t119 + mrSges(6,2) * t120;
t75 = -mrSges(7,1) * t119 + mrSges(7,2) * t120;
t66 = Ifges(6,1) * t120 + Ifges(6,4) * t119 - Ifges(6,5) * t223;
t65 = Ifges(7,1) * t120 + Ifges(7,4) * t119 - Ifges(7,5) * t223;
t64 = Ifges(6,4) * t120 + Ifges(6,2) * t119 - Ifges(6,6) * t223;
t63 = Ifges(7,4) * t120 + Ifges(7,2) * t119 - Ifges(7,6) * t223;
t54 = qJ(6) * t139 + t79;
t53 = -qJ(6) * t140 + t78;
t52 = -mrSges(5,1) * t104 + mrSges(5,2) * t103;
t51 = Ifges(5,1) * t103 + Ifges(5,4) * t104 + Ifges(5,5) * t259;
t50 = Ifges(5,4) * t103 + Ifges(5,2) * t104 + Ifges(5,6) * t259;
t49 = -pkin(5) * t61 + t130;
t39 = mrSges(6,1) * t259 - mrSges(6,3) * t47;
t35 = Ifges(6,1) * t60 + Ifges(6,4) * t61;
t34 = Ifges(7,1) * t60 + Ifges(7,4) * t61;
t33 = Ifges(6,4) * t60 + Ifges(6,2) * t61;
t32 = Ifges(7,4) * t60 + Ifges(7,2) * t61;
t31 = -mrSges(6,1) * t61 + mrSges(6,2) * t60;
t26 = qJ(6) * t119 + t37;
t25 = -pkin(5) * t223 - t120 * qJ(6) + t36;
t24 = -pkin(5) * t48 + t81;
t15 = -mrSges(6,1) * t48 + mrSges(6,2) * t47;
t13 = Ifges(6,1) * t47 + Ifges(6,4) * t48 + Ifges(6,5) * t259;
t12 = Ifges(7,1) * t47 + Ifges(7,4) * t48 + Ifges(7,5) * t259;
t11 = Ifges(6,4) * t47 + Ifges(6,2) * t48 + Ifges(6,6) * t259;
t10 = Ifges(7,4) * t47 + Ifges(7,2) * t48 + Ifges(7,6) * t259;
t1 = [(t2 * t25 + t24 * t82 + t26 * t3) * t286 + (t144 * t81 + t36 * t6 + t37 * t5) * t287 + (t158 * t192 + t28 * t95 + t29 * t94) * t288 + (t161 * t114 + t160 * t115) * t289 + 0.2e1 * t114 * t187 + 0.2e1 * t115 * t188 + 0.2e1 * t192 * t52 - t165 * t50 - t166 * t51 + 0.2e1 * t28 * t154 + 0.2e1 * t29 * t155 + 0.2e1 * t158 * t127 + 0.2e1 * t160 * t156 + 0.2e1 * t161 * t157 + 0.2e1 * t144 * t15 + t104 * t117 + t103 * t118 + 0.2e1 * t3 * t110 + 0.2e1 * t5 * t111 + 0.2e1 * t2 * t112 + 0.2e1 * t6 * t113 + 0.2e1 * t94 * t92 + 0.2e1 * t95 * t93 + 0.2e1 * t24 * t75 + 0.2e1 * t81 * t76 + 0.2e1 * t82 * t14 + 0.2e1 * t25 * t38 + 0.2e1 * t36 * t39 + 0.2e1 * t26 * t40 + 0.2e1 * t37 * t41 + (t138 * t284 - t218 * t124 + t222 * t125 + (-t222 * t163 - t218 * t164 + t223 * t236) * qJD(3) + (-Ifges(5,5) * t166 - Ifges(5,6) * t165 + mrSges(3,1) * t285 + (Ifges(4,5) * t222 - t243) * t219 + t294 * t120 + t293 * t119 + (pkin(7) ^ 2 * t289 + t239 * t284 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3) - t296) * t223) * qJD(2)) * t219 + (t65 + t66) * t47 + (t64 + t63) * t48 + (t12 + t13) * t120 + ((mrSges(3,2) * t285 - t218 * t163 + t222 * t164 + t223 * t243) * qJD(2) - t234 + t290 + t292) * t223 + (t10 + t11) * t119; (t87 / 0.2e1 + t88 / 0.2e1) * t48 + (pkin(8) * t157 + t114 * mrSges(4,3) + t124 / 0.2e1) * t222 + (-pkin(8) * t156 - t115 * mrSges(4,3) + t125 / 0.2e1) * t218 + (-pkin(2) * t214 + (t114 * t222 - t115 * t218 + (-t160 * t222 - t161 * t218) * qJD(3)) * pkin(8)) * m(4) + t207 * t52 + t180 * t51 / 0.2e1 + t192 * t97 - t165 * t98 / 0.2e1 - t166 * t99 / 0.2e1 + t153 * t93 + t106 * t154 + t107 * t155 + t158 * t148 + t162 * t15 + t144 * t31 - t146 * t118 / 0.2e1 - t147 * t117 / 0.2e1 + t104 * t149 / 0.2e1 + t103 * t150 / 0.2e1 + t152 * t92 - pkin(2) * t138 + t130 * t76 + t108 * t14 + t8 * t110 + t19 * t111 + t9 * t112 + t20 * t113 + t24 * t85 + t81 * t86 + t49 * t75 + t78 * t39 + t79 * t41 + t82 * t30 + t53 * t38 + t54 * t40 + (t12 / 0.2e1 + t13 / 0.2e1) * t140 + (t10 / 0.2e1 + t11 / 0.2e1) * t139 + (t34 / 0.2e1 + t35 / 0.2e1) * t120 + (t33 / 0.2e1 + t32 / 0.2e1) * t119 + m(6) * (t130 * t144 + t162 * t81 + t19 * t37 + t20 * t36 + t5 * t79 + t6 * t78) + m(7) * (t108 * t24 + t2 * t53 + t25 * t9 + t26 * t8 + t3 * t54 + t49 * t82) + (t64 / 0.2e1 + t63 / 0.2e1) * t61 + (t65 / 0.2e1 + t66 / 0.2e1) * t60 + (t139 * t3 - t140 * t2 - t25 * t60 + t26 * t61) * mrSges(7,3) + (t139 * t5 - t140 * t6 - t36 * t60 + t37 * t61) * mrSges(6,3) + (-t58 / 0.2e1 - t56 / 0.2e1 + t143 / 0.2e1 + t142 / 0.2e1 - t59 / 0.2e1 - t57 / 0.2e1 - t212 / 0.2e1 + (Ifges(3,5) + t265 / 0.2e1 + t194 * t281 + (-mrSges(3,1) + t240) * pkin(7)) * qJD(2)) * t223 + m(5) * (t106 * t95 + t107 * t94 + t152 * t29 + t153 * t28 + t158 * t207 + t192 * t213) + ((t164 / 0.2e1 - pkin(8) * t188 - t160 * mrSges(4,3)) * t222 + (-t163 / 0.2e1 + t270 / 0.2e1 - pkin(8) * t187 - t161 * mrSges(4,3) + pkin(3) * t127) * t218) * qJD(3) + (t89 / 0.2e1 + t90 / 0.2e1) * t47 - t235 * t50 / 0.2e1 + (t222 * t186 / 0.2e1 + t185 * t281 - Ifges(3,6) * qJD(2) + (-t222 * t194 / 0.2e1 + t195 * t281) * qJD(3) + (qJD(2) * mrSges(3,2) + t184) * pkin(7) + (Ifges(5,5) * t180 - Ifges(5,6) * t235 + t293 * t139 + t294 * t140 + t236) * qJD(2) / 0.2e1) * t219 + (t146 * t94 - t147 * t95 - t180 * t29 - t235 * t28) * mrSges(5,3); -0.2e1 * pkin(2) * t184 + 0.2e1 * t108 * t30 + 0.2e1 * t130 * t86 - t146 * t150 - t147 * t149 + 0.2e1 * t162 * t31 - t235 * t98 + t180 * t99 + t222 * t185 + t218 * t186 + 0.2e1 * t207 * t97 + 0.2e1 * t49 * t85 + (t87 + t88) * t61 + (t89 + t90) * t60 + (t34 + t35) * t140 + (t33 + t32) * t139 + (t265 + (0.2e1 * pkin(3) * t148 - t194) * t218) * qJD(3) + (t106 * t153 + t107 * t152 + t207 * t213) * t288 + (t130 * t162 + t19 * t79 + t20 * t78) * t287 + (t108 * t49 + t53 * t9 + t54 * t8) * t286 + 0.2e1 * (t139 * t8 - t140 * t9 - t53 * t60 + t54 * t61) * mrSges(7,3) + 0.2e1 * (t139 * t19 - t140 * t20 - t60 * t78 + t61 * t79) * mrSges(6,3) + 0.2e1 * (-t106 * t235 - t107 * t180 + t146 * t152 - t147 * t153) * mrSges(5,3); t224 + t167 * t38 + t168 * t39 - t114 * mrSges(4,2) + t115 * mrSges(4,1) - Ifges(4,5) * t247 + (m(5) * (t217 * t28 + t221 * t29 + t253 * t95 - t254 * t94) + t221 * t92 + t217 * t93 + t154 * t253 - t155 * t254) * pkin(3) + t275 * t169 + t262 * t135 + t263 * t134 + m(7) * (t134 * t26 + t135 * t25 + t167 * t2 + t169 * t3) + m(6) * (t134 * t37 + t135 * t36 + t168 * t6 + t169 * t5) - t230 * Ifges(4,6) - t292; t212 + (m(5) * (t106 * t217 + t107 * t221 + (-t152 * t217 + t153 * t221) * qJD(4)) + (t221 * t146 - t217 * t147 + (t180 * t217 - t221 * t235) * qJD(4)) * mrSges(5,3)) * pkin(3) + t225 + (pkin(8) * t240 - t272) * qJD(3) + m(7) * (t134 * t54 + t135 * t53 + t167 * t9 + t169 * t8) + m(6) * (t134 * t79 + t135 * t78 + t168 * t20 + t169 * t19) + (-t168 * t60 + t229) * mrSges(6,3) + (-t167 * t60 + t229) * mrSges(7,3); 0.2e1 * t131 + 0.2e1 * t132 + t134 * t242 + 0.2e1 * t232 + (t135 * t167 + t109) * t286 + (t135 * t168 + t109) * t287; t224 + t250 * t205 + (t220 * t39 + t275 * t216 + (-t216 * t262 + t220 * t263) * qJD(5) + m(7) * (t216 * t3 - t25 * t252 + t251 * t26) + m(6) * (t216 * t5 + t220 * t6 + t251 * t37 - t252 * t36)) * pkin(4); t241 * t205 + (m(7) * (t216 * t8 + t251 * t54 - t252 * t53) + m(6) * (t19 * t216 + t20 * t220 + t251 * t79 - t252 * t78) + t227 * mrSges(7,3) + (-t220 * t60 + t227) * mrSges(6,3)) * pkin(4) + t225; t232 + (t216 * t277 + t244) * t271 + m(7) * (-pkin(4) * t167 * t252 + t135 * t205 + t261) + m(6) * ((t135 * t220 - t168 * t252) * pkin(4) + t261) + t233; (t242 * t280 + 0.2e1 * ((-t205 + t280) * m(7) + t277) * t295) * qJD(5); pkin(5) * t250 + t226; pkin(5) * t241 + t228; t135 * t283 + t233; (t244 + (t277 - t283) * t216) * t271; 0; m(7) * t24 + t14; m(7) * t49 + t30; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
