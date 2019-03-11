% Calculate time derivative of joint inertia matrix for
% S6RPRRRP2
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:54
% EndTime: 2019-03-09 05:59:01
% DurationCPUTime: 3.79s
% Computational Cost: add. (3954->443), mult. (9310->619), div. (0->0), fcn. (7811->8), ass. (0->179)
t151 = sin(qJ(3));
t149 = sin(qJ(5));
t150 = sin(qJ(4));
t152 = cos(qJ(5));
t153 = cos(qJ(4));
t162 = t149 * t150 - t152 * t153;
t154 = cos(qJ(3));
t188 = qJD(3) * t154;
t115 = t149 * t153 + t150 * t152;
t219 = qJD(4) + qJD(5);
t79 = t219 * t115;
t47 = -t151 * t79 - t162 * t188;
t101 = t162 * t151;
t48 = t101 * t219 - t115 * t188;
t13 = -t48 * mrSges(6,1) + t47 * mrSges(6,2);
t206 = t48 * mrSges(7,1) - t47 * mrSges(7,2);
t163 = -t13 + t206;
t185 = qJD(4) * t153;
t159 = t150 * t188 + t151 * t185;
t172 = t153 * t188;
t186 = qJD(4) * t151;
t176 = t150 * t186;
t160 = t172 - t176;
t70 = t159 * mrSges(5,1) + t160 * mrSges(5,2);
t228 = -t70 + t163;
t227 = Ifges(6,5) + Ifges(7,5);
t226 = Ifges(6,6) + Ifges(7,6);
t225 = -Ifges(6,3) - Ifges(7,3);
t180 = -cos(pkin(10)) * pkin(1) - pkin(2);
t224 = 0.2e1 * t180;
t223 = pkin(4) * t152;
t222 = -mrSges(6,1) - mrSges(7,1);
t221 = pkin(4) * qJD(5);
t214 = -pkin(9) - pkin(8);
t133 = t214 * t150;
t134 = t214 * t153;
t88 = t149 * t133 - t152 * t134;
t111 = -pkin(3) * t154 - t151 * pkin(8) + t180;
t139 = sin(pkin(10)) * pkin(1) + pkin(7);
t192 = t153 * t154;
t124 = t139 * t192;
t77 = t150 * t111 + t124;
t189 = qJD(3) * t151;
t220 = -Ifges(5,5) * t172 - Ifges(5,3) * t189;
t190 = t150 ^ 2 + t153 ^ 2;
t218 = 2 * m(5);
t217 = 2 * m(6);
t216 = 2 * m(7);
t215 = 0.2e1 * t139;
t213 = pkin(5) * t48;
t203 = Ifges(5,4) * t150;
t130 = Ifges(5,2) * t153 + t203;
t212 = -t130 / 0.2e1;
t211 = -t150 / 0.2e1;
t210 = pkin(3) * t151;
t209 = pkin(8) * t154;
t29 = -mrSges(7,2) * t189 + mrSges(7,3) * t48;
t30 = -mrSges(6,2) * t189 + mrSges(6,3) * t48;
t208 = t29 + t30;
t183 = qJD(5) * t152;
t207 = (-t101 * t183 + t149 * t47) * pkin(4);
t103 = t153 * t111;
t193 = t151 * t153;
t59 = -pkin(9) * t193 + t103 + (-t139 * t150 - pkin(4)) * t154;
t194 = t150 * t151;
t68 = -pkin(9) * t194 + t77;
t23 = t149 * t59 + t152 * t68;
t100 = t115 * t151;
t89 = mrSges(7,2) * t154 - t100 * mrSges(7,3);
t90 = mrSges(6,2) * t154 - t100 * mrSges(6,3);
t205 = t89 + t90;
t91 = -mrSges(7,1) * t154 + t101 * mrSges(7,3);
t92 = -mrSges(6,1) * t154 + t101 * mrSges(6,3);
t204 = t91 + t92;
t202 = Ifges(5,4) * t153;
t201 = Ifges(5,5) * t150;
t200 = Ifges(5,6) * t150;
t199 = Ifges(5,6) * t153;
t127 = (-t209 + t210) * qJD(3);
t178 = t139 * t189;
t191 = t153 * t127 + t150 * t178;
t46 = -qJD(4) * t77 + t191;
t198 = t150 * t46;
t197 = t154 * Ifges(5,6);
t129 = -mrSges(5,1) * t153 + mrSges(5,2) * t150;
t196 = -mrSges(4,1) + t129;
t195 = t139 * t154;
t104 = pkin(4) * t194 + t151 * t139;
t187 = qJD(4) * t150;
t184 = qJD(5) * t149;
t182 = pkin(4) * t187;
t24 = (pkin(4) * t151 - pkin(9) * t192) * qJD(3) + (-t124 + (pkin(9) * t151 - t111) * t150) * qJD(4) + t191;
t171 = t154 * t187;
t45 = t150 * t127 + t111 * t185 + (-t153 * t189 - t171) * t139;
t26 = -pkin(9) * t159 + t45;
t6 = -qJD(5) * t23 - t149 * t26 + t152 * t24;
t2 = pkin(5) * t189 - qJ(6) * t47 + qJD(6) * t101 + t6;
t27 = mrSges(7,1) * t189 - mrSges(7,3) * t47;
t181 = m(7) * t2 + t27;
t128 = t139 * t188;
t86 = pkin(4) * t159 + t128;
t141 = -pkin(4) * t153 - pkin(3);
t179 = qJD(4) * t214;
t177 = t151 * t188;
t174 = t100 * t184;
t78 = t219 * t162;
t32 = t79 * mrSges(7,1) - t78 * mrSges(7,2);
t170 = (-mrSges(6,2) - mrSges(7,2)) * t152;
t169 = (2 * Ifges(4,4)) + t200;
t22 = -t149 * t68 + t152 * t59;
t76 = -t150 * t195 + t103;
t168 = -t76 * qJD(4) + t45;
t87 = t152 * t133 + t134 * t149;
t125 = t150 * t179;
t126 = t153 * t179;
t51 = -qJD(5) * t88 - t125 * t149 + t152 * t126;
t18 = qJ(6) * t78 - qJD(6) * t115 + t51;
t167 = m(7) * t18 + t78 * mrSges(7,3);
t166 = mrSges(5,1) * t150 + mrSges(5,2) * t153;
t165 = Ifges(5,1) * t153 - t203;
t131 = Ifges(5,1) * t150 + t202;
t164 = -Ifges(5,2) * t150 + t202;
t161 = t225 * t189 - t226 * t48 - t227 * t47;
t5 = t149 * t24 + t152 * t26 + t59 * t183 - t184 * t68;
t50 = t152 * t125 + t149 * t126 + t133 * t183 + t134 * t184;
t17 = -qJ(6) * t79 - qJD(6) * t162 + t50;
t72 = Ifges(7,6) * t79;
t73 = Ifges(6,6) * t79;
t74 = Ifges(7,5) * t78;
t75 = Ifges(6,5) * t78;
t157 = t51 * mrSges(6,1) + t18 * mrSges(7,1) - t50 * mrSges(6,2) - t17 * mrSges(7,2) - t72 - t73 - t74 - t75;
t156 = -t149 * t79 + (t115 * t149 - t152 * t162) * qJD(5);
t3 = qJ(6) * t48 - qJD(6) * t100 + t5;
t155 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) - t161;
t145 = Ifges(5,5) * t185;
t140 = pkin(5) + t223;
t123 = -mrSges(5,1) * t154 - mrSges(5,3) * t193;
t122 = mrSges(5,2) * t154 - mrSges(5,3) * t194;
t121 = t165 * qJD(4);
t120 = t164 * qJD(4);
t119 = t166 * qJD(4);
t110 = t166 * t151;
t99 = -Ifges(5,5) * t154 + t151 * t165;
t98 = t151 * t164 - t197;
t96 = pkin(5) * t162 + t141;
t95 = -mrSges(5,2) * t189 - mrSges(5,3) * t159;
t94 = mrSges(5,1) * t189 - mrSges(5,3) * t160;
t85 = Ifges(6,1) * t115 - Ifges(6,4) * t162;
t84 = Ifges(7,1) * t115 - Ifges(7,4) * t162;
t83 = Ifges(6,4) * t115 - Ifges(6,2) * t162;
t82 = Ifges(7,4) * t115 - Ifges(7,2) * t162;
t81 = mrSges(6,1) * t162 + mrSges(6,2) * t115;
t80 = mrSges(7,1) * t162 + mrSges(7,2) * t115;
t69 = pkin(5) * t100 + t104;
t67 = pkin(5) * t79 + t182;
t66 = mrSges(6,1) * t100 - mrSges(6,2) * t101;
t65 = mrSges(7,1) * t100 - mrSges(7,2) * t101;
t64 = -qJ(6) * t162 + t88;
t63 = -qJ(6) * t115 + t87;
t62 = -t131 * t186 + (Ifges(5,5) * t151 + t154 * t165) * qJD(3);
t61 = -t130 * t186 + (Ifges(5,6) * t151 + t154 * t164) * qJD(3);
t58 = -Ifges(6,1) * t101 - Ifges(6,4) * t100 - Ifges(6,5) * t154;
t57 = -Ifges(7,1) * t101 - Ifges(7,4) * t100 - Ifges(7,5) * t154;
t56 = -Ifges(6,4) * t101 - Ifges(6,2) * t100 - Ifges(6,6) * t154;
t55 = -Ifges(7,4) * t101 - Ifges(7,2) * t100 - Ifges(7,6) * t154;
t37 = -Ifges(6,1) * t78 - Ifges(6,4) * t79;
t36 = -Ifges(7,1) * t78 - Ifges(7,4) * t79;
t35 = -Ifges(6,4) * t78 - Ifges(6,2) * t79;
t34 = -Ifges(7,4) * t78 - Ifges(7,2) * t79;
t33 = mrSges(6,1) * t79 - mrSges(6,2) * t78;
t28 = mrSges(6,1) * t189 - mrSges(6,3) * t47;
t19 = t86 - t213;
t16 = -qJ(6) * t100 + t23;
t14 = -pkin(5) * t154 + t101 * qJ(6) + t22;
t11 = Ifges(6,1) * t47 + Ifges(6,4) * t48 + Ifges(6,5) * t189;
t10 = Ifges(7,1) * t47 + Ifges(7,4) * t48 + Ifges(7,5) * t189;
t9 = Ifges(6,4) * t47 + Ifges(6,2) * t48 + Ifges(6,6) * t189;
t8 = Ifges(7,4) * t47 + Ifges(7,2) * t48 + Ifges(7,6) * t189;
t1 = [-(t10 + t11) * t101 - (t8 + t9) * t100 + (t55 + t56) * t48 + (t57 + t58) * t47 + (t70 * t215 - t150 * t61 + t153 * t62 + (-t150 * t99 - t154 * (-t199 - t201) - t153 * t98) * qJD(4) + (mrSges(4,1) * t224 + (Ifges(5,5) * t153 - t169) * t151 - t227 * t101 - t226 * t100 + (t139 ^ 2 * t218 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) + t225) * t154) * qJD(3)) * t151 - 0.2e1 * t69 * t206 + (t14 * t2 + t16 * t3 + t19 * t69) * t216 + (t104 * t86 + t22 * t6 + t23 * t5) * t217 + (t77 * t45 + t76 * t46) * t218 + ((mrSges(4,2) * t224 + t110 * t215 - t150 * t98 + t153 * t99 + t154 * t169) * qJD(3) + t161 + t220) * t154 + 0.2e1 * t14 * t27 + 0.2e1 * t22 * t28 + 0.2e1 * t16 * t29 + 0.2e1 * t23 * t30 + 0.2e1 * t19 * t65 + 0.2e1 * t86 * t66 + 0.2e1 * t3 * t89 + 0.2e1 * t5 * t90 + 0.2e1 * t2 * t91 + 0.2e1 * t6 * t92 + 0.2e1 * t76 * t94 + 0.2e1 * t77 * t95 + 0.2e1 * t104 * t13 + 0.2e1 * t45 * t122 + 0.2e1 * t46 * t123; t204 * t48 + t205 * t47 - t208 * t101 - (t27 + t28) * t100 + ((t122 * t153 - t123 * t150) * qJD(3) + t228) * t154 + (-t150 * t94 + t153 * t95 + (-t122 * t150 - t123 * t153) * qJD(4) + (t110 + t65 + t66) * qJD(3)) * t151 + m(7) * (-t2 * t100 - t3 * t101 + t14 * t48 - t154 * t19 + t16 * t47 + t189 * t69) + m(6) * (-t6 * t100 - t5 * t101 + t104 * t189 - t154 * t86 + t22 * t48 + t23 * t47) + m(5) * ((-t150 * t76 + t153 * t77 - t195) * t188 + (t178 - t198 + t153 * t45 + (-t150 * t77 - t153 * t76) * qJD(4)) * t151); 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (-t100 * t48 - t101 * t47 - t177) + (-0.1e1 + t190) * t177 * t218; m(7) * (t14 * t18 + t16 * t17 + t19 * t96 + t2 * t63 + t3 * t64 + t67 * t69) - (t34 / 0.2e1 + t35 / 0.2e1) * t100 + (-t46 * mrSges(5,3) - pkin(8) * t94 + t62 / 0.2e1 + t188 * t212 + (pkin(4) * t66 - t77 * mrSges(5,3) - pkin(8) * t122 - t98 / 0.2e1 + t197 / 0.2e1) * qJD(4)) * t150 + (qJD(4) * t99 / 0.2e1 + t61 / 0.2e1 + t131 * t188 / 0.2e1 + t168 * mrSges(5,3) + (m(5) * t168 - qJD(4) * t123 + t95) * pkin(8)) * t153 + (t10 / 0.2e1 + t11 / 0.2e1) * t115 + m(6) * (t104 * t182 + t141 * t86 + t22 * t51 + t23 * t50 + t5 * t88 + t6 * t87) - (t36 / 0.2e1 + t37 / 0.2e1) * t101 + (-t145 / 0.2e1 + t74 / 0.2e1 + t72 / 0.2e1 + t75 / 0.2e1 + t73 / 0.2e1 + (t139 * t196 + Ifges(4,5)) * qJD(3)) * t154 - (t55 / 0.2e1 + t56 / 0.2e1) * t79 - (t57 / 0.2e1 + t58 / 0.2e1) * t78 + (-t115 * t6 - t162 * t5 + t22 * t78 - t23 * t79) * mrSges(6,3) + (-t115 * t2 + t14 * t78 - t16 * t79 - t162 * t3) * mrSges(7,3) - t96 * t206 - (t8 / 0.2e1 + t9 / 0.2e1) * t162 + (t139 * t119 + t120 * t211 + t153 * t121 / 0.2e1 + (t131 * t211 + t153 * t212) * qJD(4) + (t139 * mrSges(4,2) - Ifges(4,6) + t201 / 0.2e1 + t199 / 0.2e1 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t115 - (Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t162) * qJD(3)) * t151 + t63 * t27 + t64 * t29 + t67 * t65 + t69 * t32 - pkin(3) * t70 + t19 * t80 + t86 * t81 + t87 * t28 + t88 * t30 + t17 * t89 + t50 * t90 + t18 * t91 + t51 * t92 + m(5) * (-pkin(3) * t128 + (-t187 * t77 - t198) * pkin(8)) + t104 * t33 + (t84 / 0.2e1 + t85 / 0.2e1) * t47 + (t82 / 0.2e1 + t83 / 0.2e1) * t48 + t141 * t13; (-t119 - t32 - t33) * t154 + m(7) * (-t18 * t100 - t17 * t101 - t154 * t67 + t64 * t47 + t63 * t48) + m(6) * (-pkin(4) * t171 - t51 * t100 - t50 * t101 + t88 * t47 + t87 * t48) + ((mrSges(5,3) * t190 - mrSges(4,2)) * t154 + m(5) * (t190 * t209 - t210) + (m(6) * t141 + m(7) * t96 + t196 + t80 + t81) * t151) * qJD(3) + (mrSges(7,3) + mrSges(6,3)) * (-t100 * t78 + t101 * t79 - t115 * t48 - t162 * t47); -0.2e1 * pkin(3) * t119 + t153 * t120 + t150 * t121 + 0.2e1 * t141 * t33 + 0.2e1 * t96 * t32 + 0.2e1 * t67 * t80 - (t82 + t83) * t79 - (t84 + t85) * t78 + (t36 + t37) * t115 - (t34 + t35) * t162 + (t153 * t131 + (0.2e1 * pkin(4) * t81 - t130) * t150) * qJD(4) + (t17 * t64 + t18 * t63 + t67 * t96) * t216 + (t141 * t182 + t50 * t88 + t51 * t87) * t217 + 0.2e1 * (-t115 * t18 - t162 * t17 + t63 * t78 - t64 * t79) * mrSges(7,3) + 0.2e1 * (-t115 * t51 - t162 * t50 + t78 * t87 - t79 * t88) * mrSges(6,3); -t159 * Ifges(5,6) + t155 + t181 * t140 + (t152 * t28 + t208 * t149 + (-t149 * t204 + t152 * t205) * qJD(5) + m(6) * (t149 * t5 + t152 * t6 + t183 * t23 - t184 * t22) + m(7) * (-t14 * t184 + t149 * t3 + t16 * t183)) * pkin(4) - Ifges(5,5) * t176 - t45 * mrSges(5,2) + t46 * mrSges(5,1) - t220; m(7) * (pkin(4) * t174 + t140 * t48 + t207) + m(6) * ((t152 * t48 + t174) * pkin(4) + t207) + t228; t145 + t167 * t140 + (pkin(8) * t129 - t200) * qJD(4) + (m(7) * (t149 * t17 + t183 * t64 - t184 * t63) + m(6) * (t149 * t50 + t152 * t51 + t183 * t88 - t184 * t87) + t156 * mrSges(7,3) + (t152 * t78 + t156) * mrSges(6,3)) * pkin(4) + t157; 0.2e1 * (t170 + ((-t140 + t223) * m(7) + t222) * t149) * t221; pkin(5) * t181 + t155; m(7) * t213 + t163; pkin(5) * t167 + t157; (t170 + (-m(7) * pkin(5) + t222) * t149) * t221; 0; m(7) * t19 - t206; m(7) * t189; m(7) * t67 + t32; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
