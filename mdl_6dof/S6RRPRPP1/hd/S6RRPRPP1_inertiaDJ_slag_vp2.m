% Calculate time derivative of joint inertia matrix for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:10
% EndTime: 2019-03-09 09:45:21
% DurationCPUTime: 4.46s
% Computational Cost: add. (4652->403), mult. (10658->575), div. (0->0), fcn. (9897->8), ass. (0->172)
t229 = Ifges(7,4) + Ifges(6,5);
t228 = -Ifges(7,6) + Ifges(6,6);
t149 = sin(pkin(9));
t150 = cos(pkin(9));
t152 = sin(qJ(2));
t154 = cos(qJ(2));
t130 = t149 * t154 + t150 * t152;
t123 = t130 * qJD(2);
t128 = t149 * t152 - t150 * t154;
t125 = t128 * qJD(2);
t148 = sin(pkin(10));
t151 = sin(qJ(4));
t187 = qJD(4) * t151;
t179 = t130 * t187;
t153 = cos(qJ(4));
t192 = t153 * t125;
t162 = t179 + t192;
t197 = cos(pkin(10));
t172 = t197 * t153;
t173 = t197 * t151;
t188 = qJD(4) * t130;
t45 = t125 * t173 + t148 * t162 - t172 * t188;
t129 = t148 * t153 + t173;
t161 = -t148 * t151 + t172;
t46 = -t125 * t161 - t129 * t188;
t227 = (Ifges(6,1) + Ifges(7,1)) * t46 + (Ifges(6,4) - Ifges(7,5)) * t45 + t229 * t123;
t186 = qJD(4) * t153;
t178 = t130 * t186;
t163 = -t151 * t125 + t178;
t226 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t122 = t129 * qJD(4);
t124 = t161 * qJD(4);
t225 = Ifges(5,5) * t186 - t228 * t122 + t229 * t124;
t224 = m(6) + m(7);
t210 = pkin(4) * t148;
t140 = qJ(6) + t210;
t223 = m(7) * t140 + mrSges(7,3);
t222 = 2 * m(6);
t221 = 0.2e1 * m(7);
t220 = -2 * mrSges(4,3);
t208 = -qJ(3) - pkin(7);
t136 = t208 * t152;
t137 = t208 * t154;
t98 = -t150 * t136 - t137 * t149;
t218 = 0.2e1 * t98;
t146 = -pkin(2) * t154 - pkin(1);
t217 = 0.2e1 * t146;
t216 = m(6) * pkin(4);
t164 = qJ(5) * t125 - qJD(5) * t130;
t174 = qJD(2) * t208;
t121 = qJD(3) * t154 + t152 * t174;
t157 = -t152 * qJD(3) + t154 * t174;
t71 = t150 * t121 + t149 * t157;
t182 = pkin(2) * qJD(2) * t152;
t72 = pkin(3) * t123 + pkin(8) * t125 + t182;
t175 = -t151 * t71 + t153 * t72;
t87 = t128 * pkin(3) - t130 * pkin(8) + t146;
t99 = t136 * t149 - t137 * t150;
t96 = t153 * t99;
t7 = pkin(4) * t123 + t164 * t153 + (-t96 + (qJ(5) * t130 - t87) * t151) * qJD(4) + t175;
t184 = t151 * t72 + t153 * t71 + t87 * t186;
t9 = -qJ(5) * t178 + (-qJD(4) * t99 + t164) * t151 + t184;
t4 = t148 * t7 + t197 * t9;
t215 = t123 / 0.2e1;
t212 = -t130 / 0.2e1;
t203 = Ifges(5,4) * t151;
t138 = Ifges(5,2) * t153 + t203;
t211 = -t138 / 0.2e1;
t70 = t121 * t149 - t150 * t157;
t209 = t70 * t98;
t23 = -mrSges(6,2) * t123 + mrSges(6,3) * t45;
t26 = mrSges(7,2) * t45 + mrSges(7,3) * t123;
t207 = t23 + t26;
t24 = mrSges(6,1) * t123 - mrSges(6,3) * t46;
t25 = -t123 * mrSges(7,1) + t46 * mrSges(7,2);
t206 = -t24 + t25;
t195 = t130 * t153;
t51 = -t151 * t99 + t153 * t87;
t28 = pkin(4) * t128 - qJ(5) * t195 + t51;
t196 = t130 * t151;
t52 = t151 * t87 + t96;
t44 = -qJ(5) * t196 + t52;
t13 = t148 * t28 + t197 * t44;
t75 = t129 * t130;
t57 = -mrSges(6,2) * t128 - mrSges(6,3) * t75;
t60 = -mrSges(7,2) * t75 + mrSges(7,3) * t128;
t205 = t57 + t60;
t76 = t161 * t130;
t58 = mrSges(6,1) * t128 - mrSges(6,3) * t76;
t59 = -mrSges(7,1) * t128 + mrSges(7,2) * t76;
t204 = -t58 + t59;
t202 = Ifges(5,4) * t153;
t201 = t123 * Ifges(5,5);
t200 = t123 * Ifges(5,6);
t199 = t124 * mrSges(7,2);
t198 = t128 * Ifges(5,6);
t143 = pkin(2) * t149 + pkin(8);
t191 = qJ(5) + t143;
t185 = qJD(6) * t129;
t181 = pkin(4) * t187;
t145 = -pkin(2) * t150 - pkin(3);
t180 = t197 * pkin(4);
t19 = -t45 * mrSges(6,1) + t46 * mrSges(6,2);
t18 = -t45 * mrSges(7,1) - t46 * mrSges(7,3);
t168 = qJD(4) * t191;
t103 = qJD(5) * t153 - t151 * t168;
t156 = -qJD(5) * t151 - t153 * t168;
t62 = t103 * t148 - t156 * t197;
t63 = t103 * t197 + t148 * t156;
t126 = t191 * t153;
t169 = t191 * t151;
t79 = t126 * t148 + t169 * t197;
t80 = t126 * t197 - t148 * t169;
t177 = t62 * t79 + t80 * t63;
t176 = -Ifges(5,6) * t151 - (2 * Ifges(4,4));
t171 = 0.2e1 * t182;
t82 = t122 * mrSges(6,1) + t124 * mrSges(6,2);
t81 = t122 * mrSges(7,1) - t124 * mrSges(7,3);
t69 = pkin(4) * t196 + t98;
t167 = mrSges(5,1) * t151 + mrSges(5,2) * t153;
t166 = Ifges(5,1) * t153 - t203;
t165 = -Ifges(5,2) * t151 + t202;
t135 = -pkin(4) * t153 + t145;
t49 = pkin(4) * t163 + t70;
t3 = -t148 * t9 + t197 * t7;
t12 = -t148 * t44 + t197 * t28;
t159 = -t81 - t82;
t158 = -Ifges(5,5) * t192 + t226 * t123 + t228 * t45 + t229 * t46;
t144 = -t180 - pkin(5);
t139 = Ifges(5,1) * t151 + t202;
t134 = t166 * qJD(4);
t133 = t165 * qJD(4);
t132 = t167 * qJD(4);
t112 = t125 * mrSges(4,2);
t95 = Ifges(6,1) * t129 + Ifges(6,4) * t161;
t94 = Ifges(7,1) * t129 - Ifges(7,5) * t161;
t93 = Ifges(6,4) * t129 + Ifges(6,2) * t161;
t92 = Ifges(7,5) * t129 - Ifges(7,3) * t161;
t91 = -mrSges(6,1) * t161 + mrSges(6,2) * t129;
t90 = -mrSges(7,1) * t161 - mrSges(7,3) * t129;
t89 = mrSges(5,1) * t128 - mrSges(5,3) * t195;
t88 = -mrSges(5,2) * t128 - mrSges(5,3) * t196;
t86 = Ifges(6,1) * t124 - Ifges(6,4) * t122;
t85 = Ifges(7,1) * t124 + Ifges(7,5) * t122;
t84 = Ifges(6,4) * t124 - Ifges(6,2) * t122;
t83 = Ifges(7,5) * t124 + Ifges(7,3) * t122;
t73 = -pkin(5) * t161 - qJ(6) * t129 + t135;
t65 = Ifges(5,5) * t128 + t130 * t166;
t64 = t130 * t165 + t198;
t56 = pkin(5) * t122 - qJ(6) * t124 + t181 - t185;
t55 = -mrSges(5,2) * t123 - mrSges(5,3) * t163;
t54 = mrSges(5,1) * t123 + mrSges(5,3) * t162;
t50 = mrSges(5,1) * t163 - mrSges(5,2) * t162;
t48 = mrSges(6,1) * t75 + mrSges(6,2) * t76;
t47 = mrSges(7,1) * t75 - mrSges(7,3) * t76;
t34 = -Ifges(5,1) * t162 - Ifges(5,4) * t163 + t201;
t33 = -Ifges(5,4) * t162 - Ifges(5,2) * t163 + t200;
t32 = Ifges(6,1) * t76 - Ifges(6,4) * t75 + Ifges(6,5) * t128;
t31 = Ifges(7,1) * t76 + Ifges(7,4) * t128 + Ifges(7,5) * t75;
t30 = Ifges(6,4) * t76 - Ifges(6,2) * t75 + Ifges(6,6) * t128;
t29 = Ifges(7,5) * t76 + Ifges(7,6) * t128 + Ifges(7,3) * t75;
t22 = pkin(5) * t75 - qJ(6) * t76 + t69;
t21 = -t52 * qJD(4) + t175;
t20 = -t187 * t99 + t184;
t15 = Ifges(6,4) * t46 + Ifges(6,2) * t45 + t123 * Ifges(6,6);
t14 = Ifges(7,5) * t46 + t123 * Ifges(7,6) - Ifges(7,3) * t45;
t11 = -t128 * pkin(5) - t12;
t10 = qJ(6) * t128 + t13;
t5 = -pkin(5) * t45 - qJ(6) * t46 - qJD(6) * t76 + t49;
t2 = -t123 * pkin(5) - t3;
t1 = qJ(6) * t123 + qJD(6) * t128 + t4;
t6 = [-t112 * t217 + t50 * t218 + (t1 * t10 + t11 * t2 + t22 * t5) * t221 + (t12 * t3 + t13 * t4 + t49 * t69) * t222 + 0.2e1 * ((-t152 ^ 2 + t154 ^ 2) * Ifges(3,4) - pkin(1) * (mrSges(3,1) * t152 + mrSges(3,2) * t154) + (-Ifges(3,2) + Ifges(3,1)) * t152 * t154) * qJD(2) + t227 * t76 + (mrSges(4,1) * t217 + t99 * t220 - t228 * t75 + t229 * t76) * t123 + (t30 - t29) * t45 + (mrSges(4,2) * t171 - 0.2e1 * Ifges(4,1) * t125 - t151 * t33 + t153 * t34 + (Ifges(5,5) * t153 + t176) * t123 + (-t151 * t65 - t153 * t64 + t128 * (-Ifges(5,5) * t151 - Ifges(5,6) * t153)) * qJD(4) + 0.2e1 * (t167 + mrSges(4,3)) * t70) * t130 + (mrSges(4,1) * t171 + t71 * t220 - t176 * t125 + ((2 * Ifges(4,2)) + t226) * t123 + t158) * t128 + 0.2e1 * t20 * t88 + 0.2e1 * t21 * t89 + 0.2e1 * t69 * t19 + 0.2e1 * t4 * t57 + 0.2e1 * t3 * t58 + 0.2e1 * t2 * t59 + 0.2e1 * t1 * t60 + 0.2e1 * t51 * t54 + 0.2e1 * t52 * t55 + 0.2e1 * t5 * t47 + 0.2e1 * t49 * t48 + 0.2e1 * t10 * t26 + 0.2e1 * t13 * t23 + 0.2e1 * t12 * t24 + 0.2e1 * t11 * t25 + 0.2e1 * t22 * t18 - (mrSges(4,3) * t218 - t151 * t64 + t153 * t65) * t125 + 0.2e1 * m(5) * (t20 * t52 + t21 * t51 + t209) + 0.2e1 * m(4) * (t146 * t182 + t71 * t99 + t209) + (t14 - t15) * t75 + (t31 + t32) * t46; ((-t123 * t149 + t125 * t150) * mrSges(4,3) + m(4) * (t149 * t71 - t150 * t70)) * pkin(2) + (Ifges(3,5) * t154 - Ifges(3,6) * t152 + (-mrSges(3,1) * t154 + mrSges(3,2) * t152) * pkin(7)) * qJD(2) + m(7) * (t1 * t80 + t10 * t63 + t11 * t62 + t2 * t79 + t22 * t56 + t5 * t73) + (t4 * mrSges(6,3) + t1 * mrSges(7,2) - t14 / 0.2e1 + t15 / 0.2e1 + t228 * t215) * t161 + (t31 / 0.2e1 + t32 / 0.2e1) * t124 + (-t30 / 0.2e1 + t29 / 0.2e1) * t122 + (t201 / 0.2e1 + t70 * mrSges(5,2) + t34 / 0.2e1 - t21 * mrSges(5,3) + t133 * t212 - t125 * t211 + (t139 * t212 - t52 * mrSges(5,3) + pkin(4) * t48 - t198 / 0.2e1 - t64 / 0.2e1 + t69 * t216) * qJD(4) + (-m(5) * t21 - t54 + (-m(5) * t52 - t88) * qJD(4)) * t143) * t151 + (t200 / 0.2e1 - t70 * mrSges(5,1) + t33 / 0.2e1 + t20 * mrSges(5,3) + t130 * t134 / 0.2e1 - t125 * t139 / 0.2e1 + (t130 * t211 - t51 * mrSges(5,3) + t65 / 0.2e1) * qJD(4) + (-qJD(4) * t89 + t55 + m(5) * (-qJD(4) * t51 + t20)) * t143) * t153 + (-t92 / 0.2e1 + t93 / 0.2e1) * t45 + t225 * t128 / 0.2e1 + (t85 / 0.2e1 + t86 / 0.2e1) * t76 + (t83 / 0.2e1 - t84 / 0.2e1) * t75 + (-t12 * t124 - t122 * t13) * mrSges(6,3) + (-t10 * t122 + t11 * t124) * mrSges(7,2) + (m(5) * t70 + t50) * t145 + (t94 / 0.2e1 + t95 / 0.2e1) * t46 + t135 * t19 + t98 * t132 + (t227 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3) + t229 * t215) * t129 - Ifges(4,5) * t125 - Ifges(4,6) * t123 + t22 * t81 + t69 * t82 + t5 * t90 + t49 * t91 + t73 * t18 - t70 * mrSges(4,1) - t71 * mrSges(4,2) + t56 * t47 + m(6) * (-t12 * t62 + t13 * t63 + t135 * t49 - t3 * t79 + t4 * t80) + t206 * t79 + t207 * t80 + t204 * t62 + t205 * t63; 0.2e1 * t145 * t132 + t153 * t133 + t151 * t134 + 0.2e1 * t135 * t82 + 0.2e1 * t56 * t90 + 0.2e1 * t73 * t81 + (t85 + t86) * t129 - (t83 - t84) * t161 + (t94 + t95) * t124 + (-t93 + t92) * t122 + (t153 * t139 + (0.2e1 * pkin(4) * t91 - t138) * t151) * qJD(4) + (t56 * t73 + t177) * t221 + (t135 * t181 + t177) * t222 + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t122 * t80 + t124 * t79 + t129 * t62 + t161 * t63); m(4) * t182 + t123 * mrSges(4,1) + t151 * t55 + t153 * t54 - t112 + t207 * t129 - t206 * t161 + t205 * t124 + t204 * t122 + (-t151 * t89 + t153 * t88) * qJD(4) + m(6) * (-t12 * t122 + t124 * t13 + t129 * t4 + t161 * t3) + m(7) * (t1 * t129 + t10 * t124 + t11 * t122 - t161 * t2) + m(5) * (t151 * t20 + t153 * t21 + (-t151 * t51 + t153 * t52) * qJD(4)); t224 * (t122 * t79 + t124 * t80 + t129 * t63 - t161 * t62); 0.2e1 * t224 * (-t122 * t161 + t129 * t124); (t148 * t4 + t197 * t3) * t216 + t158 - Ifges(5,5) * t179 + t144 * t25 + t140 * t26 + qJD(6) * t60 - t20 * mrSges(5,2) + t21 * mrSges(5,1) - t4 * mrSges(6,2) + t3 * mrSges(6,1) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + m(7) * (qJD(6) * t10 + t1 * t140 + t144 * t2) + t24 * t180 + t23 * t210 - t163 * Ifges(5,6); m(7) * qJD(6) * t80 - Ifges(5,6) * t187 + t144 * t199 + (t148 * t216 - mrSges(6,2) + t223) * t63 + (m(7) * t144 - t197 * t216 - mrSges(6,1) - mrSges(7,1)) * t62 + (-mrSges(5,1) * t186 + mrSges(5,2) * t187) * t143 + (-t122 * t210 - t124 * t180) * mrSges(6,3) + (qJD(6) * t161 - t122 * t140) * mrSges(7,2) + t225; -mrSges(5,2) * t186 - mrSges(5,1) * t187 + (-t122 * t197 + t124 * t148) * t216 + m(7) * (t122 * t144 + t124 * t140 + t185) + t159; 0.2e1 * t223 * qJD(6); m(6) * t49 + m(7) * t5 + t18 + t19; m(6) * t181 + m(7) * t56 - t159; 0; 0; 0; m(7) * t2 + t25; m(7) * t62 + t199; m(7) * t122; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
