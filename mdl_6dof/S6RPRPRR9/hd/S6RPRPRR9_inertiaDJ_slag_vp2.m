% Calculate time derivative of joint inertia matrix for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:22
% EndTime: 2019-03-09 04:01:34
% DurationCPUTime: 5.86s
% Computational Cost: add. (16451->536), mult. (48988->824), div. (0->0), fcn. (52794->14), ass. (0->239)
t181 = sin(qJ(6));
t184 = cos(qJ(6));
t182 = sin(qJ(5));
t226 = qJD(6) * t182;
t185 = cos(qJ(5));
t230 = qJD(5) * t185;
t296 = t181 * t226 - t184 * t230;
t299 = -2 * Ifges(5,4);
t175 = sin(pkin(12));
t177 = sin(pkin(6));
t178 = cos(pkin(12));
t180 = cos(pkin(6));
t183 = sin(qJ(3));
t179 = cos(pkin(7));
t186 = cos(qJ(3));
t246 = t179 * t186;
t176 = sin(pkin(7));
t249 = t176 * t186;
t114 = t180 * t249 + (-t175 * t183 + t178 * t246) * t177;
t298 = 0.2e1 * t114;
t247 = t179 * t183;
t250 = t176 * t183;
t115 = t177 * (t175 * t186 + t178 * t247) + t180 * t250;
t297 = 0.2e1 * t115;
t276 = t184 / 0.2e1;
t248 = t177 * t178;
t136 = -t176 * t248 + t179 * t180;
t174 = sin(pkin(13));
t254 = cos(pkin(13));
t88 = t174 * t114 + t115 * t254;
t193 = t136 * t185 - t182 * t88;
t109 = t114 * qJD(3);
t110 = t115 * qJD(3);
t85 = t109 * t254 - t174 * t110;
t57 = qJD(5) * t193 + t185 * t85;
t77 = t136 * t182 + t185 * t88;
t87 = -t114 * t254 + t115 * t174;
t59 = -t181 * t77 + t184 * t87;
t84 = t109 * t174 + t110 * t254;
t28 = qJD(6) * t59 + t181 * t84 + t184 * t57;
t60 = t181 * t87 + t184 * t77;
t29 = -qJD(6) * t60 - t181 * t57 + t184 * t84;
t12 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t44 = mrSges(6,1) * t84 - mrSges(6,3) * t57;
t267 = t12 - t44;
t34 = -mrSges(7,1) * t59 + mrSges(7,2) * t60;
t62 = mrSges(6,1) * t87 - mrSges(6,3) * t77;
t266 = t34 - t62;
t239 = t181 ^ 2 + t184 ^ 2;
t225 = qJD(6) * t184;
t190 = -t181 * t230 - t182 * t225;
t295 = Ifges(4,5) * t109 + Ifges(5,5) * t85 - Ifges(4,6) * t110 - Ifges(5,6) * t84;
t206 = t254 * t186;
t129 = (-t174 * t183 + t206) * t176 * qJD(3);
t132 = (t174 * t186 + t183 * t254) * t176;
t121 = t132 * t185 + t179 * t182;
t234 = qJD(5) * t121;
t93 = t129 * t182 + t234;
t255 = t93 * t182;
t120 = t132 * t182 - t185 * t179;
t235 = qJD(5) * t120;
t92 = t129 * t185 - t235;
t294 = t92 * t185 + t255;
t232 = qJD(5) * t182;
t275 = pkin(1) * t180;
t240 = qJ(2) * t248 + t175 * t275;
t113 = (t176 * t180 + t179 * t248) * pkin(9) + t240;
t106 = t186 * t113;
t165 = t178 * t275;
t251 = t175 * t177;
t116 = pkin(2) * t180 + t165 + (-pkin(9) * t179 - qJ(2)) * t251;
t127 = (-pkin(9) * t175 * t176 - pkin(2) * t178 - pkin(1)) * t177;
t188 = (-t106 + (-t116 * t179 - t127 * t176) * t183) * qJD(3);
t192 = -t109 * qJ(4) - t115 * qJD(4);
t238 = qJD(2) * t177;
t219 = t175 * t238;
t203 = t179 * t219;
t218 = t178 * t238;
t105 = t116 * t246;
t236 = qJD(3) * t186;
t217 = t176 * t236;
t71 = qJD(3) * t105 + t127 * t217 + t186 * t218 + (-qJD(3) * t113 - t203) * t183;
t52 = -qJ(4) * t110 + qJD(4) * t114 + t71;
t32 = t254 * t52 + (-t183 * t218 - t186 * t203 + t188 + t192) * t174;
t73 = -t113 * t183 + t127 * t249 + t105;
t65 = pkin(3) * t136 - qJ(4) * t115 + t73;
t74 = t116 * t247 + t127 * t250 + t106;
t70 = qJ(4) * t114 + t74;
t39 = t174 * t65 + t254 * t70;
t37 = pkin(10) * t136 + t39;
t89 = -t116 * t176 + t179 * t127;
t75 = -pkin(3) * t114 + t89;
t48 = pkin(4) * t87 - pkin(10) * t88 + t75;
t149 = t176 * t219;
t98 = pkin(3) * t110 + t149;
t50 = pkin(4) * t84 - pkin(10) * t85 + t98;
t5 = t182 * t50 + t185 * t32 + t48 * t230 - t232 * t37;
t265 = t182 * t48 + t185 * t37;
t292 = qJD(5) * t265;
t6 = -t182 * t32 + t185 * t50 - t292;
t293 = -t182 * t6 + t185 * t5;
t285 = m(5) * pkin(3);
t291 = t174 * t285 - mrSges(5,2);
t151 = -mrSges(7,1) * t184 + mrSges(7,2) * t181;
t290 = -m(7) * pkin(5) - mrSges(6,1) + t151;
t220 = t254 * pkin(3);
t168 = -t220 - pkin(4);
t289 = m(6) * t168 - mrSges(6,1) * t185 + mrSges(6,2) * t182 - t254 * t285 - mrSges(5,1);
t288 = 0.2e1 * m(7);
t287 = m(6) / 0.2e1;
t286 = m(7) / 0.2e1;
t58 = qJD(5) * t77 + t182 * t85;
t11 = Ifges(7,1) * t28 + Ifges(7,4) * t29 + Ifges(7,5) * t58;
t284 = t11 / 0.2e1;
t283 = t29 / 0.2e1;
t282 = -t193 / 0.2e1;
t261 = Ifges(7,4) * t181;
t154 = Ifges(7,2) * t184 + t261;
t260 = Ifges(7,4) * t184;
t198 = -Ifges(7,2) * t181 + t260;
t100 = -t154 * t226 + (Ifges(7,6) * t182 + t185 * t198) * qJD(5);
t281 = t100 / 0.2e1;
t156 = Ifges(7,1) * t181 + t260;
t199 = Ifges(7,1) * t184 - t261;
t101 = -t156 * t226 + (Ifges(7,5) * t182 + t185 * t199) * qJD(5);
t280 = t101 / 0.2e1;
t135 = -Ifges(7,5) * t185 + t182 * t199;
t279 = t135 / 0.2e1;
t278 = Ifges(7,5) * t181 / 0.2e1 + Ifges(7,6) * t276;
t277 = -t181 / 0.2e1;
t274 = pkin(3) * t174;
t273 = pkin(5) * t182;
t272 = pkin(11) * t185;
t269 = t84 * mrSges(5,3);
t268 = t85 * mrSges(5,3);
t263 = Ifges(6,4) * t182;
t262 = Ifges(6,4) * t185;
t259 = Ifges(7,6) * t181;
t258 = t120 * t93;
t131 = t174 * t250 - t176 * t206;
t72 = (-t175 * t246 - t178 * t183) * t238 + t188;
t31 = t174 * t52 - t254 * (t192 + t72);
t257 = t131 * t31;
t128 = qJD(3) * t132;
t253 = t128 * t131;
t245 = t181 * t182;
t244 = t181 * t185;
t243 = t182 * t184;
t242 = t184 * t185;
t237 = qJD(3) * t183;
t233 = qJD(5) * t181;
t231 = qJD(5) * t184;
t138 = -t185 * pkin(5) - t182 * pkin(11) + t168;
t167 = pkin(10) + t274;
t118 = t138 * t184 - t167 * t244;
t229 = qJD(6) * t118;
t119 = t138 * t181 + t167 * t242;
t228 = qJD(6) * t119;
t227 = qJD(6) * t181;
t224 = t185 * t288;
t9 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t58;
t223 = Ifges(6,5) * t57 - Ifges(6,6) * t58 + Ifges(6,3) * t84;
t222 = mrSges(6,3) * t232;
t216 = t167 * t232;
t210 = t120 * t230;
t63 = t84 * mrSges(5,1) + t85 * mrSges(5,2);
t209 = -t230 / 0.2e1;
t208 = t230 / 0.2e1;
t207 = -t226 / 0.2e1;
t148 = (-t272 + t273) * qJD(5);
t90 = t148 * t181 - t184 * t216 + t229;
t205 = t90 - t229;
t91 = t148 * t184 + t181 * t216 - t228;
t204 = -t91 - t228;
t13 = pkin(5) * t58 - pkin(11) * t57 + t31;
t3 = pkin(11) * t84 + t5;
t15 = pkin(11) * t87 + t265;
t38 = -t174 * t70 + t254 * t65;
t36 = -t136 * pkin(4) - t38;
t20 = -pkin(5) * t193 - t77 * pkin(11) + t36;
t7 = -t15 * t181 + t184 * t20;
t1 = qJD(6) * t7 + t13 * t181 + t184 * t3;
t8 = t15 * t184 + t181 * t20;
t2 = -qJD(6) * t8 + t13 * t184 - t181 * t3;
t202 = t1 * t184 - t181 * t2;
t201 = t182 * mrSges(6,1) + t185 * mrSges(6,2);
t200 = mrSges(7,1) * t181 + mrSges(7,2) * t184;
t16 = mrSges(7,1) * t58 - mrSges(7,3) * t28;
t17 = -mrSges(7,2) * t58 + mrSges(7,3) * t29;
t197 = -t181 * t16 + t184 * t17;
t94 = -t121 * t181 + t131 * t184;
t66 = qJD(6) * t94 + t128 * t181 + t184 * t92;
t95 = t121 * t184 + t131 * t181;
t67 = -qJD(6) * t95 + t128 * t184 - t181 * t92;
t196 = -t181 * t67 + t184 * t66;
t18 = -t182 * t37 + t185 * t48;
t189 = -t225 * t7 - t227 * t8 + t202;
t99 = -Ifges(7,5) * t296 + Ifges(7,6) * t190 + Ifges(7,3) * t232;
t171 = Ifges(6,5) * t230;
t170 = Ifges(7,5) * t225;
t157 = Ifges(6,1) * t182 + t262;
t155 = Ifges(6,2) * t185 + t263;
t147 = -mrSges(7,1) * t185 - mrSges(7,3) * t243;
t146 = mrSges(7,2) * t185 - mrSges(7,3) * t245;
t145 = (Ifges(6,1) * t185 - t263) * qJD(5);
t144 = t199 * qJD(6);
t143 = (-Ifges(6,2) * t182 + t262) * qJD(5);
t142 = t198 * qJD(6);
t141 = -Ifges(7,6) * t227 + t170;
t140 = t201 * qJD(5);
t139 = t200 * qJD(6);
t137 = t200 * t182;
t134 = -Ifges(7,6) * t185 + t182 * t198;
t133 = -Ifges(7,3) * t185 + (Ifges(7,5) * t184 - t259) * t182;
t126 = -mrSges(7,2) * t232 + mrSges(7,3) * t190;
t125 = mrSges(7,1) * t232 + mrSges(7,3) * t296;
t112 = mrSges(7,1) * t190 + mrSges(7,2) * t296;
t97 = mrSges(4,1) * t136 - mrSges(4,3) * t115;
t96 = -mrSges(4,2) * t136 + mrSges(4,3) * t114;
t86 = mrSges(4,1) * t110 + mrSges(4,2) * t109;
t79 = mrSges(5,1) * t136 - mrSges(5,3) * t88;
t78 = -mrSges(5,2) * t136 - mrSges(5,3) * t87;
t61 = -mrSges(6,2) * t87 + mrSges(6,3) * t193;
t53 = -mrSges(6,1) * t193 + mrSges(6,2) * t77;
t45 = -mrSges(6,2) * t84 - mrSges(6,3) * t58;
t43 = Ifges(6,1) * t77 + Ifges(6,4) * t193 + Ifges(6,5) * t87;
t42 = Ifges(6,4) * t77 + Ifges(6,2) * t193 + t87 * Ifges(6,6);
t41 = -mrSges(7,1) * t193 - mrSges(7,3) * t60;
t40 = mrSges(7,2) * t193 + mrSges(7,3) * t59;
t33 = mrSges(6,1) * t58 + mrSges(6,2) * t57;
t25 = Ifges(6,1) * t57 - Ifges(6,4) * t58 + t84 * Ifges(6,5);
t24 = Ifges(6,4) * t57 - Ifges(6,2) * t58 + t84 * Ifges(6,6);
t23 = Ifges(7,1) * t60 + Ifges(7,4) * t59 - Ifges(7,5) * t193;
t22 = Ifges(7,4) * t60 + Ifges(7,2) * t59 - Ifges(7,6) * t193;
t21 = Ifges(7,5) * t60 + Ifges(7,6) * t59 - Ifges(7,3) * t193;
t14 = -pkin(5) * t87 - t18;
t10 = Ifges(7,4) * t28 + Ifges(7,2) * t29 + Ifges(7,6) * t58;
t4 = -pkin(5) * t84 - t6;
t19 = [(t21 - t42) * t58 - (t9 - t24) * t193 + (t88 * t299 + Ifges(6,5) * t77 - Ifges(5,6) * t136 + Ifges(6,6) * t193 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t87) * t84 + (t1 * t8 + t14 * t4 + t2 * t7) * t288 - (0.2e1 * mrSges(4,3) * t74 + Ifges(4,4) * t297 + Ifges(4,2) * t298 + Ifges(4,6) * t136) * t110 + (-0.2e1 * mrSges(4,3) * t73 + Ifges(4,1) * t297 + Ifges(4,4) * t298 + Ifges(4,5) * t136) * t109 + (0.2e1 * Ifges(5,1) * t88 + Ifges(5,5) * t136 + t299 * t87) * t85 + t295 * t136 + 0.2e1 * m(5) * (-t31 * t38 + t32 * t39 + t75 * t98) + 0.2e1 * t71 * t96 + 0.2e1 * t72 * t97 + 0.2e1 * t98 * (mrSges(5,1) * t87 + mrSges(5,2) * t88) + 0.2e1 * t89 * t86 + 0.2e1 * t75 * t63 + t77 * t25 + 0.2e1 * t32 * t78 - 0.2e1 * t31 * t79 + t60 * t11 + 0.2e1 * t5 * t61 + 0.2e1 * t6 * t62 + t57 * t43 + t59 * t10 + 0.2e1 * t31 * t53 + 0.2e1 * t1 * t40 + 0.2e1 * t2 * t41 - 0.2e1 * t38 * t268 - 0.2e1 * t39 * t269 + 0.2e1 * t18 * t44 + 0.2e1 * m(3) * (t240 * t178 + (qJ(2) * t251 - t165) * t175) * t238 + 0.2e1 * (-mrSges(4,1) * t114 + mrSges(4,2) * t115) * t149 + 0.2e1 * m(4) * (t149 * t89 + t71 * t74 + t72 * t73) + 0.2e1 * m(6) * (t18 * t6 + t265 * t5 + t31 * t36) + 0.2e1 * t265 * t45 + t87 * t223 + 0.2e1 * t14 * t12 + 0.2e1 * t7 * t16 + 0.2e1 * t8 * t17 + 0.2e1 * (-mrSges(3,2) * t180 + mrSges(3,3) * t248) * t218 - 0.2e1 * (mrSges(3,1) * t180 - mrSges(3,3) * t251) * t219 + t28 * t23 + t29 * t22 + 0.2e1 * t4 * t34 + 0.2e1 * t36 * t33; t121 * t45 + t129 * t78 + t131 * t33 + t94 * t16 + t95 * t17 + t66 * t40 + t67 * t41 + t92 * t61 + t266 * t93 + (t86 + t63) * t179 + (t53 - t79) * t128 + t267 * t120 + (t131 * t85 - t132 * t84) * mrSges(5,3) + m(7) * (t1 * t95 + t120 * t4 + t14 * t93 + t2 * t94 + t66 * t8 + t67 * t7) + m(6) * (-t120 * t6 + t121 * t5 + t128 * t36 - t18 * t93 + t265 * t92 + t257) + m(5) * (-t128 * t38 + t129 * t39 + t132 * t32 + t179 * t98 + t257) + (t96 * t236 - t97 * t237 + m(4) * (t183 * t71 + t186 * t72 + t236 * t74 - t237 * t73 + t203) + (-t109 * t186 - t110 * t183) * mrSges(4,3)) * t176; 0.2e1 * m(7) * (t66 * t95 + t67 * t94 + t258) + 0.2e1 * m(6) * (t121 * t92 + t253 + t258) + 0.2e1 * m(5) * (t129 * t132 + t253); t43 * t208 + (-t155 / 0.2e1 + t133 / 0.2e1) * t58 + t295 + t193 * t143 / 0.2e1 + (t266 * t230 + t267 * t182 + m(6) * ((-t18 * t185 - t182 * t265) * qJD(5) + t293) + t185 * t45 + m(7) * (t14 * t230 + t182 * t4)) * t167 + t243 * t284 + m(7) * (t1 * t119 + t118 * t2 + t7 * t91 + t8 * t90) + t28 * t279 + t60 * t280 + t59 * t281 + t99 * t282 + t134 * t283 + (-t18 * t230 + t293) * mrSges(6,3) + t291 * t32 + t289 * t31 + (t207 * t22 + t208 * t23) * t184 + (t207 * t23 + t209 * t22) * t181 + t84 * (Ifges(6,5) * t182 + Ifges(6,6) * t185) / 0.2e1 - t185 * t9 / 0.2e1 + t185 * t24 / 0.2e1 + t182 * t25 / 0.2e1 + t168 * t33 + t77 * t145 / 0.2e1 + t1 * t146 + t2 * t147 + t57 * t157 / 0.2e1 + t4 * t137 + t36 * t140 + t7 * t125 + t8 * t126 + t118 * t16 + t119 * t17 - t14 * t112 + t90 * t40 + t91 * t41 - t71 * mrSges(4,2) + t72 * mrSges(4,1) - t269 * t274 - t220 * t268 - t265 * t222 - t61 * t216 - t42 * t232 / 0.2e1 + t21 * t232 / 0.2e1 + t87 * (-Ifges(6,6) * t232 + t171) / 0.2e1 - t10 * t245 / 0.2e1; m(7) * (t118 * t67 + t119 * t66 + t90 * t95 + t91 * t94 + (t210 + t255) * t167) + t93 * t137 - t120 * t112 + t66 * t146 + t95 * t126 + t67 * t147 + t94 * t125 - t121 * t222 + t131 * t140 + m(6) * ((t120 * t185 - t121 * t182) * qJD(5) + t294) * t167 - mrSges(4,2) * t217 - t176 * mrSges(4,1) * t237 + t291 * t129 + t289 * t128 + (t210 + t294) * mrSges(6,3); 0.2e1 * t91 * t147 + 0.2e1 * t118 * t125 + 0.2e1 * t90 * t146 + 0.2e1 * t119 * t126 + (t118 * t91 + t119 * t90) * t288 + 0.2e1 * t168 * t140 + (t143 - t99 + (-t134 * t181 + t135 * t184 + 0.2e1 * t137 * t167 + t157) * qJD(5)) * t185 + (-t181 * t100 + t184 * t101 - 0.2e1 * t167 * t112 + t145 + (-t134 * t184 - t135 * t181) * qJD(6) + (t167 ^ 2 * t224 + t133 - t155) * qJD(5)) * t182; m(5) * t98 + ((-t181 * t41 + t184 * t40 + t61) * qJD(5) + m(7) * (t231 * t8 - t233 * t7 - t4) + m(6) * (t6 + t292) - t267) * t185 + (t45 + (-t181 * t40 - t184 * t41) * qJD(6) + t266 * qJD(5) + m(7) * (qJD(5) * t14 + t189) + m(6) * (-t18 * qJD(5) + t5) + t197) * t182 + t63; 0.2e1 * ((t231 * t95 - t233 * t94 - t93) * t286 + (-t93 + t234) * t287) * t185 + 0.2e1 * ((-t225 * t94 - t227 * t95 + t196 + t235) * t286 + (t92 + t235) * t287) * t182; t185 * t112 + (-t147 * t225 - t181 * t125 + m(7) * (-t118 * t225 - t119 * t227 - t181 * t91 + t184 * t90) - t146 * t227 + t184 * t126) * t182 + (-t147 * t244 + t182 * t137 + m(7) * (-t118 * t244 + t119 * t242 + (t182 ^ 2 - t185 ^ 2) * t167) + t146 * t242) * qJD(5); (-0.1e1 + t239) * t224 * t232; t10 * t276 + t181 * t284 + t59 * t142 / 0.2e1 + t60 * t144 / 0.2e1 + t4 * t151 + t58 * t278 + t154 * t283 + t28 * t156 / 0.2e1 + t14 * t139 + t141 * t282 - t5 * mrSges(6,2) + t6 * mrSges(6,1) + (t22 * t277 + t23 * t276) * qJD(6) + (-m(7) * t4 - t12) * pkin(5) + ((-t181 * t8 - t184 * t7) * qJD(6) + t202) * mrSges(7,3) + (m(7) * t189 - t225 * t41 - t227 * t40 + t197) * pkin(11) + t223; -t92 * mrSges(6,2) + t120 * t139 + (m(7) * pkin(11) + mrSges(7,3)) * ((-t181 * t95 - t184 * t94) * qJD(6) + t196) + t290 * t93; pkin(5) * t112 + t171 + (-t141 / 0.2e1 + t290 * t167 * qJD(5)) * t185 + (t156 * t208 + t281 + qJD(6) * t279 + t205 * mrSges(7,3) + (m(7) * t205 - qJD(6) * t147 + t126) * pkin(11)) * t184 + (t154 * t209 + t280 - qJD(6) * t134 / 0.2e1 + t204 * mrSges(7,3) + (m(7) * t204 - qJD(6) * t146 - t125) * pkin(11)) * t181 + (t167 * t139 + t142 * t277 + t144 * t276 + (-t184 * t154 / 0.2e1 + t156 * t277) * qJD(6) + (t167 * mrSges(6,2) - Ifges(6,6) + t278) * qJD(5)) * t182; -t185 * t139 + (t182 * t151 + m(7) * (t239 * t272 - t273) + t239 * t185 * mrSges(7,3) - t201) * qJD(5); -0.2e1 * pkin(5) * t139 + t142 * t184 + t144 * t181 + (-t154 * t181 + t156 * t184) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t9; mrSges(7,1) * t67 - mrSges(7,2) * t66; mrSges(7,1) * t91 - mrSges(7,2) * t90 + t99; t112; t170 + (pkin(11) * t151 - t259) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
