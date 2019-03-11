% Calculate time derivative of joint inertia matrix for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:07:05
% EndTime: 2019-03-09 06:07:15
% DurationCPUTime: 4.92s
% Computational Cost: add. (6717->362), mult. (14672->512), div. (0->0), fcn. (15047->8), ass. (0->159)
t280 = Ifges(6,4) + Ifges(7,4);
t279 = Ifges(6,1) + Ifges(7,1);
t278 = Ifges(6,2) + Ifges(7,2);
t154 = sin(pkin(10));
t158 = sin(qJ(3));
t155 = cos(pkin(10));
t161 = cos(qJ(3));
t208 = t155 * t161;
t120 = -t158 * t154 + t208;
t121 = t154 * t161 + t158 * t155;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t168 = t160 * t120 - t121 * t157;
t233 = Ifges(6,5) + Ifges(7,5);
t277 = t168 * t233;
t266 = Ifges(6,6) + Ifges(7,6);
t276 = t168 * t266;
t159 = cos(qJ(5));
t275 = t280 * t159;
t156 = sin(qJ(5));
t274 = t280 * t156;
t272 = -t278 * t156 + t275;
t271 = t279 * t159 - t274;
t270 = t233 * t156 + t266 * t159;
t268 = t278 * t159 + t274;
t267 = t279 * t156 + t275;
t255 = (t156 ^ 2 + t159 ^ 2) * t160;
t265 = Ifges(6,3) + Ifges(7,3);
t204 = qJD(5) * t156;
t106 = t120 * qJD(3);
t107 = t121 * qJD(3);
t70 = qJD(4) * t168 + t106 * t160 - t107 * t157;
t214 = t159 * t70;
t95 = t120 * t157 + t121 * t160;
t165 = t95 * t204 - t214;
t203 = qJD(5) * t159;
t195 = t95 * t203;
t220 = t156 * t70;
t166 = t195 + t220;
t71 = qJD(4) * t95 + t106 * t157 + t160 * t107;
t264 = -t165 * t280 - t278 * t166 + t266 * t71;
t263 = -t279 * t165 - t166 * t280 + t233 * t71;
t262 = t272 * t95 - t276;
t228 = t271 * t95 - t277;
t261 = t272 * qJD(5);
t260 = t271 * qJD(5);
t205 = t233 * t203;
t231 = pkin(7) + qJ(2);
t129 = t231 * t154;
t118 = t161 * t129;
t130 = t231 * t155;
t98 = -t130 * t158 - t118;
t91 = -pkin(8) * t121 + t98;
t99 = -t158 * t129 + t161 * t130;
t92 = pkin(8) * t120 + t99;
t256 = -t157 * t92 + t160 * t91;
t51 = t157 * t91 + t160 * t92;
t190 = -pkin(2) * t155 - pkin(1);
t100 = -pkin(3) * t120 + t190;
t52 = -pkin(4) * t168 - pkin(9) * t95 + t100;
t30 = -t156 * t51 + t159 * t52;
t47 = t159 * t51;
t31 = t156 * t52 + t47;
t254 = -t156 * t30 + t159 * t31;
t164 = t260 * t156 + t261 * t159 + t203 * t267 - t204 * t268;
t253 = 2 * m(6);
t252 = 2 * m(7);
t251 = -2 * mrSges(5,3);
t89 = -t121 * qJD(2) - qJD(3) * t99;
t163 = -t106 * pkin(8) + t89;
t88 = -qJD(3) * t118 + qJD(2) * t208 + (-qJD(2) * t154 - qJD(3) * t130) * t158;
t77 = -pkin(8) * t107 + t88;
t25 = qJD(4) * t51 + t157 * t77 - t160 * t163;
t250 = 0.2e1 * t25;
t249 = -0.2e1 * t256;
t122 = mrSges(7,1) * t204 + mrSges(7,2) * t203;
t248 = 0.2e1 * t122;
t133 = -mrSges(7,1) * t159 + mrSges(7,2) * t156;
t247 = 0.2e1 * t133;
t246 = -0.2e1 * t156;
t245 = 0.2e1 * t159;
t244 = m(7) * pkin(5);
t241 = mrSges(7,3) * pkin(5);
t238 = pkin(3) * t107;
t237 = pkin(3) * t160;
t24 = qJD(4) * t256 + t157 * t163 + t160 * t77;
t34 = pkin(4) * t71 - pkin(9) * t70 + t238;
t201 = t156 * t34 + t159 * t24 + t52 * t203;
t5 = -t204 * t51 + t201;
t236 = t159 * t5;
t235 = t25 * t256;
t185 = -t156 * t24 + t159 * t34;
t6 = -qJD(5) * t31 + t185;
t234 = t6 * t156;
t230 = -qJ(6) - pkin(9);
t227 = mrSges(6,2) * t159;
t222 = pkin(3) * qJD(4);
t219 = t156 * t95;
t218 = t157 * mrSges(5,1);
t213 = t159 * t95;
t212 = t160 * mrSges(5,2);
t141 = pkin(3) * t157 + pkin(9);
t211 = t141 * t159;
t134 = -mrSges(6,1) * t159 + mrSges(6,2) * t156;
t207 = t157 * t134;
t149 = t159 * qJ(6);
t206 = -qJ(6) - t141;
t202 = 0.2e1 * t238;
t200 = 0.2e1 * qJD(5);
t199 = t160 * t222;
t198 = pkin(5) * t204;
t197 = mrSges(7,1) + t244;
t143 = -pkin(5) * t159 - pkin(4);
t189 = t71 * mrSges(5,1) + t70 * mrSges(5,2);
t188 = t266 * t156;
t187 = -t204 / 0.2e1;
t184 = qJD(5) * t230;
t183 = t233 * t214 + t265 * t71;
t181 = qJD(5) * t206;
t178 = mrSges(6,3) * t255;
t177 = -(2 * Ifges(5,4)) - t188;
t176 = mrSges(6,1) * t156 + t227;
t17 = -pkin(5) * t168 - t149 * t95 + t30;
t23 = -qJ(6) * t219 + t31;
t171 = -t23 * t156 - t17 * t159;
t170 = -t31 * t156 - t30 * t159;
t169 = -qJ(6) * t70 - qJD(6) * t95;
t18 = t166 * mrSges(7,1) - t165 * mrSges(7,2);
t123 = t176 * qJD(5);
t3 = -qJ(6) * t195 + (-qJD(5) * t51 + t169) * t156 + t201;
t36 = pkin(5) * t219 - t256;
t8 = pkin(5) * t166 + t25;
t162 = -t24 * mrSges(5,2) + mrSges(6,3) * t236 + Ifges(5,5) * t70 + t36 * t122 - t256 * t123 + t8 * t133 + (-mrSges(5,1) + t134) * t25 - (-t266 * t204 + t205) * t168 / 0.2e1 + t263 * t156 / 0.2e1 - t261 * t219 / 0.2e1 + t260 * t213 / 0.2e1 + t262 * t187 + t228 * t203 / 0.2e1 + t267 * (t95 * t187 + t214 / 0.2e1) + t268 * (-t220 / 0.2e1 - t195 / 0.2e1) + (-Ifges(5,6) + t270 / 0.2e1) * t71 + (t3 * mrSges(7,3) + t264 / 0.2e1) * t159;
t148 = t159 * qJD(6);
t142 = -pkin(4) - t237;
t135 = pkin(9) * t159 + t149;
t132 = t230 * t156;
t131 = t143 - t237;
t128 = t157 * t222 + t198;
t114 = t149 + t211;
t113 = t206 * t156;
t104 = -qJD(6) * t156 + t159 * t184;
t103 = t156 * t184 + t148;
t102 = t106 * mrSges(4,2);
t97 = (-qJD(6) - t199) * t156 + t159 * t181;
t96 = t156 * t181 + t159 * t199 + t148;
t75 = -mrSges(6,1) * t168 - mrSges(6,3) * t213;
t74 = -mrSges(7,1) * t168 - mrSges(7,3) * t213;
t73 = mrSges(6,2) * t168 - mrSges(6,3) * t219;
t72 = mrSges(7,2) * t168 - mrSges(7,3) * t219;
t64 = t176 * t95;
t63 = (mrSges(7,1) * t156 + mrSges(7,2) * t159) * t95;
t29 = -mrSges(6,2) * t71 - mrSges(6,3) * t166;
t28 = -mrSges(7,2) * t71 - mrSges(7,3) * t166;
t27 = mrSges(6,1) * t71 + mrSges(6,3) * t165;
t26 = mrSges(7,1) * t71 + mrSges(7,3) * t165;
t19 = mrSges(6,1) * t166 - mrSges(6,2) * t165;
t1 = pkin(5) * t71 + t169 * t159 + (-t47 + (qJ(6) * t95 - t52) * t156) * qJD(5) + t185;
t2 = [0.2e1 * t100 * t189 - (mrSges(5,1) * t202 + t24 * t251 + ((2 * Ifges(5,2)) + t265) * t71 + t177 * t70 + t183) * t168 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) * (t154 ^ 2 + t155 ^ 2) + t51 * t71 * t251 + (mrSges(5,3) * t249 - t156 * t262 + t159 * t228) * t70 + 0.2e1 * m(4) * (t88 * t99 + t89 * t98) - 0.2e1 * t120 * Ifges(4,2) * t107 + 0.2e1 * t121 * t106 * Ifges(4,1) + t19 * t249 + t64 * t250 + (t1 * t17 + t23 * t3 + t36 * t8) * t252 + (mrSges(5,2) * t202 + mrSges(5,3) * t250 + 0.2e1 * Ifges(5,1) * t70 + t263 * t159 - t264 * t156 + (t159 * t233 + t177) * t71 + ((-t262 + t276) * t159 + (-t228 + t277) * t156) * qJD(5)) * t95 + (t30 * t6 + t31 * t5 - t235) * t253 + 0.2e1 * m(5) * (t100 * t238 + t24 * t51 - t235) + 0.2e1 * t17 * t26 + 0.2e1 * t23 * t28 + 0.2e1 * t30 * t27 + 0.2e1 * t31 * t29 + 0.2e1 * t36 * t18 + 0.2e1 * t8 * t63 + 0.2e1 * t3 * t72 + 0.2e1 * t5 * t73 + 0.2e1 * t1 * t74 + 0.2e1 * t6 * t75 + 0.2e1 * (t106 * t120 - t107 * t121) * Ifges(4,4) + 0.2e1 * (-t106 * t98 - t107 * t99 + t120 * t88 - t121 * t89) * mrSges(4,3) + 0.2e1 * t190 * (t107 * mrSges(4,1) + t102); t102 + (t26 + t27) * t159 + (t28 + t29) * t156 - (-m(5) * pkin(3) - mrSges(4,1)) * t107 + ((t72 + t73) * t159 + (-t74 - t75) * t156) * qJD(5) + m(6) * (qJD(5) * t254 + t156 * t5 + t159 * t6) + m(7) * (t1 * t159 + t156 * t3 + (-t156 * t17 + t159 * t23) * qJD(5)) + t189; 0; t162 + (m(5) * (t157 * t24 - t160 * t25) + (-t157 * t71 - t160 * t70) * mrSges(5,3) + ((m(5) * t51 + m(6) * t254 + mrSges(5,3) * t168 - t156 * t75 + t159 * t73) * t160 + (t95 * mrSges(5,3) + t64 - (m(6) + m(5)) * t256) * t157) * qJD(4)) * pkin(3) + (t171 * mrSges(7,3) + t170 * mrSges(6,3) + (m(6) * t170 - t156 * t73 - t159 * t75) * t141) * qJD(5) + t29 * t211 + (-t6 * mrSges(6,3) - t1 * mrSges(7,3) - t141 * t27) * t156 + m(6) * (-t141 * t234 + t142 * t25 + t5 * t211) - t88 * mrSges(4,2) + t89 * mrSges(4,1) + t96 * t72 + t97 * t74 + Ifges(4,5) * t106 - Ifges(4,6) * t107 + t113 * t26 + t114 * t28 + t128 * t63 + t131 * t18 + t142 * t19 + m(7) * (t1 * t113 + t114 * t3 + t128 * t36 + t131 * t8 + t17 * t97 + t23 * t96); m(7) * (t156 * t96 + t159 * t97 + (-t113 * t156 + t114 * t159) * qJD(5)); (t113 * t97 + t114 * t96 + t128 * t131) * t252 + t128 * t247 + t131 * t248 + 0.2e1 * t142 * t123 + (t97 * t246 + t96 * t245 + (-t113 * t159 - t114 * t156) * t200) * mrSges(7,3) + (-0.2e1 * t212 - 0.2e1 * t218 + (t141 * t255 + t142 * t157) * t253 + 0.2e1 * t207 + 0.2e1 * t178) * t222 + t164; t162 + (qJD(5) * t170 - t234) * mrSges(6,3) + (qJD(5) * t171 - t1 * t156) * mrSges(7,3) + (m(6) * (-t203 * t30 - t204 * t31 - t234 + t236) + t159 * t29 - t156 * t27 - t73 * t204 - t75 * t203) * pkin(9) + m(7) * (t1 * t132 + t103 * t23 + t104 * t17 + t135 * t3 + t143 * t8 + t198 * t36) + t63 * t198 + t103 * t72 + t104 * t74 + t132 * t26 + t135 * t28 + t143 * t18 + (-m(6) * t25 - t19) * pkin(4); m(7) * (t103 * t156 + t104 * t159 + (-t132 * t156 + t135 * t159) * qJD(5)); m(7) * (t103 * t114 + t104 * t113 + t128 * t143 + t131 * t198 + t132 * t97 + t135 * t96) + (t128 + t198) * t133 + (-pkin(4) + t142) * t123 + (t131 + t143) * t122 + (m(6) * (-pkin(4) * t157 + pkin(9) * t255) - t218 + t207 - t212 + t178) * t222 + ((t103 + t96) * t159 + (-t104 - t97) * t156 + ((-t113 - t132) * t159 + (-t114 - t135) * t156) * qJD(5)) * mrSges(7,3) + t164; -0.2e1 * pkin(4) * t123 + t198 * t247 + t143 * t248 + (t103 * t135 + t104 * t132 + t143 * t198) * t252 + (t103 * t245 + t104 * t246 + (-t132 * t159 - t135 * t156) * t200) * mrSges(7,3) + t164; mrSges(6,1) * t6 + mrSges(7,1) * t1 - mrSges(6,2) * t5 - mrSges(7,2) * t3 - t70 * t188 + (m(7) * t1 + t26) * pkin(5) - t270 * t95 * qJD(5) + t183; (-t227 + (-mrSges(6,1) - t244) * t156) * qJD(5) - t122; -mrSges(7,2) * t96 + t197 * t97 - t176 * t199 + ((-mrSges(6,1) * t141 - t241) * t159 + (mrSges(6,2) * t141 - t266) * t156) * qJD(5) + t205; -mrSges(7,2) * t103 + t197 * t104 + ((-mrSges(6,1) * pkin(9) - t241) * t159 + (mrSges(6,2) * pkin(9) - t266) * t156) * qJD(5) + t205; 0; m(7) * t8 + t18; 0; m(7) * t128 + t122; m(7) * t198 + t122; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
