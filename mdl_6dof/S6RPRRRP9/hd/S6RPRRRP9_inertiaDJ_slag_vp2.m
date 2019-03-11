% Calculate time derivative of joint inertia matrix for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:11
% EndTime: 2019-03-09 06:26:19
% DurationCPUTime: 4.00s
% Computational Cost: add. (3830->453), mult. (8800->641), div. (0->0), fcn. (7295->6), ass. (0->182)
t214 = 2 * qJD(2);
t227 = mrSges(6,1) + mrSges(7,1);
t226 = Ifges(7,5) + Ifges(6,5);
t225 = Ifges(7,6) + Ifges(6,6);
t224 = Ifges(6,3) + Ifges(7,3);
t149 = sin(qJ(5));
t150 = sin(qJ(4));
t152 = cos(qJ(5));
t153 = cos(qJ(4));
t165 = t149 * t150 - t152 * t153;
t218 = qJD(4) + qJD(5);
t223 = t218 * t165;
t151 = sin(qJ(3));
t191 = qJD(3) * t151;
t174 = t150 * t191;
t154 = cos(qJ(3));
t186 = qJD(4) * t154;
t176 = t153 * t186;
t161 = t174 - t176;
t222 = pkin(4) * t152;
t220 = pkin(4) * qJD(5);
t212 = -pkin(9) - pkin(8);
t134 = t212 * t150;
t135 = t212 * t153;
t87 = t149 * t134 - t152 * t135;
t190 = qJD(3) * t154;
t219 = Ifges(5,6) * t174 + Ifges(5,3) * t190;
t209 = pkin(8) * t154;
t210 = pkin(3) * t151;
t129 = qJ(2) - t209 + t210;
t155 = -pkin(1) - pkin(7);
t194 = t151 * t155;
t137 = t153 * t194;
t97 = t150 * t129 + t137;
t192 = t150 ^ 2 + t153 ^ 2;
t217 = 2 * m(5);
t216 = 2 * m(6);
t215 = 2 * m(7);
t213 = m(7) * pkin(5);
t211 = -t150 / 0.2e1;
t208 = -mrSges(7,2) - mrSges(6,2);
t118 = t149 * t153 + t150 * t152;
t48 = t118 * t191 + t154 * t223;
t28 = -mrSges(7,2) * t190 + t48 * mrSges(7,3);
t29 = -mrSges(6,2) * t190 + t48 * mrSges(6,3);
t207 = t28 + t29;
t103 = t165 * t151;
t184 = qJD(5) * t152;
t104 = t165 * t154;
t78 = t218 * t118;
t45 = -qJD(3) * t104 - t151 * t78;
t206 = (-t103 * t184 + t149 * t45) * pkin(4);
t114 = t153 * t129;
t171 = -t150 * t155 + pkin(4);
t193 = t153 * t154;
t69 = -pkin(9) * t193 + t151 * t171 + t114;
t195 = t150 * t154;
t79 = -pkin(9) * t195 + t97;
t25 = t149 * t69 + t152 * t79;
t102 = t118 * t154;
t88 = -mrSges(7,2) * t151 - mrSges(7,3) * t102;
t89 = -mrSges(6,2) * t151 - mrSges(6,3) * t102;
t205 = t88 + t89;
t90 = mrSges(7,1) * t151 + mrSges(7,3) * t104;
t91 = mrSges(6,1) * t151 + mrSges(6,3) * t104;
t204 = t90 + t91;
t203 = Ifges(5,4) * t150;
t202 = Ifges(5,4) * t153;
t201 = Ifges(5,5) * t150;
t200 = Ifges(5,6) * t150;
t199 = Ifges(5,6) * t153;
t115 = qJD(2) + (pkin(3) * t154 + pkin(8) * t151) * qJD(3);
t106 = t153 * t115;
t189 = qJD(3) * t155;
t179 = t154 * t189;
t53 = -t97 * qJD(4) - t150 * t179 + t106;
t198 = t150 * t53;
t197 = t151 * Ifges(5,6);
t131 = -mrSges(5,1) * t153 + mrSges(5,2) * t150;
t196 = -mrSges(4,1) + t131;
t188 = qJD(4) * t150;
t187 = qJD(4) * t153;
t185 = qJD(5) * t149;
t183 = pkin(4) * t188;
t46 = -t154 * t78 + t165 * t191;
t22 = t106 + (-t137 + (pkin(9) * t154 - t129) * t150) * qJD(4) + (pkin(9) * t151 * t153 + t154 * t171) * qJD(3);
t178 = t151 * t188;
t52 = t150 * t115 + t129 * t187 + t153 * t179 - t155 * t178;
t30 = pkin(9) * t161 + t52;
t6 = -qJD(5) * t25 - t149 * t30 + t152 * t22;
t2 = pkin(5) * t190 - t46 * qJ(6) + qJD(6) * t104 + t6;
t26 = mrSges(7,1) * t190 - t46 * mrSges(7,3);
t182 = m(7) * t2 + t26;
t142 = -pkin(4) * t153 - pkin(3);
t181 = qJD(4) * t212;
t180 = t151 * t190;
t139 = t151 * t189;
t177 = t150 * t186;
t101 = t118 * t151;
t175 = t101 * t185;
t12 = -t48 * mrSges(7,1) + t46 * mrSges(7,2);
t32 = t78 * mrSges(7,1) - mrSges(7,2) * t223;
t173 = t208 * t152;
t172 = -Ifges(5,5) * t153 + (2 * Ifges(4,4));
t24 = -t149 * t79 + t152 * t69;
t96 = -t150 * t194 + t114;
t170 = -t96 * qJD(4) + t52;
t86 = t152 * t134 + t135 * t149;
t116 = pkin(4) * t195 - t154 * t155;
t127 = t150 * t181;
t128 = t153 * t181;
t51 = -qJD(5) * t87 - t127 * t149 + t152 * t128;
t16 = qJ(6) * t223 - qJD(6) * t118 + t51;
t169 = m(7) * t16 + mrSges(7,3) * t223;
t168 = mrSges(5,1) * t150 + mrSges(5,2) * t153;
t167 = Ifges(5,1) * t153 - t203;
t133 = Ifges(5,1) * t150 + t202;
t166 = -Ifges(5,2) * t150 + t202;
t132 = Ifges(5,2) * t153 + t203;
t47 = -qJD(3) * t102 + t151 * t223;
t164 = t208 * t45 + t227 * t47;
t163 = t190 * t224 + t225 * t48 + t226 * t46;
t5 = t149 * t22 + t152 * t30 + t69 * t184 - t185 * t79;
t92 = -pkin(4) * t161 + t139;
t50 = t152 * t127 + t149 * t128 + t134 * t184 + t135 * t185;
t162 = t153 * t191 + t177;
t15 = -qJ(6) * t78 - qJD(6) * t165 + t50;
t73 = Ifges(7,6) * t78;
t74 = Ifges(6,6) * t78;
t75 = Ifges(7,5) * t223;
t76 = Ifges(6,5) * t223;
t159 = t51 * mrSges(6,1) + t16 * mrSges(7,1) - t50 * mrSges(6,2) - t15 * mrSges(7,2) - t73 - t74 - t75 - t76;
t158 = -t149 * t78 + (t118 * t149 - t152 * t165) * qJD(5);
t3 = t48 * qJ(6) - qJD(6) * t102 + t5;
t157 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) + t163;
t146 = Ifges(5,5) * t187;
t141 = pkin(5) + t222;
t126 = mrSges(5,1) * t151 - mrSges(5,3) * t193;
t125 = -mrSges(5,2) * t151 - mrSges(5,3) * t195;
t124 = t167 * qJD(4);
t123 = t166 * qJD(4);
t122 = t168 * qJD(4);
t111 = t168 * t154;
t100 = Ifges(5,5) * t151 + t154 * t167;
t99 = t154 * t166 + t197;
t98 = pkin(5) * t165 + t142;
t95 = -mrSges(5,2) * t190 + mrSges(5,3) * t161;
t94 = mrSges(5,1) * t190 + mrSges(5,3) * t162;
t85 = Ifges(6,1) * t118 - Ifges(6,4) * t165;
t84 = Ifges(7,1) * t118 - Ifges(7,4) * t165;
t83 = Ifges(6,4) * t118 - Ifges(6,2) * t165;
t82 = Ifges(7,4) * t118 - Ifges(7,2) * t165;
t81 = mrSges(6,1) * t165 + mrSges(6,2) * t118;
t80 = mrSges(7,1) * t165 + mrSges(7,2) * t118;
t71 = t102 * pkin(5) + t116;
t68 = -mrSges(5,1) * t161 - mrSges(5,2) * t162;
t65 = pkin(5) * t78 + t183;
t63 = mrSges(6,1) * t102 - mrSges(6,2) * t104;
t62 = mrSges(7,1) * t102 - mrSges(7,2) * t104;
t61 = -qJ(6) * t165 + t87;
t60 = -qJ(6) * t118 + t86;
t59 = -t133 * t186 + (Ifges(5,5) * t154 - t151 * t167) * qJD(3);
t58 = -t132 * t186 + (Ifges(5,6) * t154 - t151 * t166) * qJD(3);
t57 = -Ifges(6,1) * t104 - Ifges(6,4) * t102 + Ifges(6,5) * t151;
t56 = -Ifges(7,1) * t104 - Ifges(7,4) * t102 + Ifges(7,5) * t151;
t55 = -Ifges(6,4) * t104 - Ifges(6,2) * t102 + Ifges(6,6) * t151;
t54 = -Ifges(7,4) * t104 - Ifges(7,2) * t102 + Ifges(7,6) * t151;
t37 = -Ifges(6,1) * t223 - Ifges(6,4) * t78;
t36 = -Ifges(7,1) * t223 - Ifges(7,4) * t78;
t35 = -Ifges(6,4) * t223 - Ifges(6,2) * t78;
t34 = -Ifges(7,4) * t223 - Ifges(7,2) * t78;
t33 = mrSges(6,1) * t78 - mrSges(6,2) * t223;
t27 = mrSges(6,1) * t190 - t46 * mrSges(6,3);
t19 = -t48 * pkin(5) + t92;
t18 = -qJ(6) * t102 + t25;
t17 = pkin(5) * t151 + qJ(6) * t104 + t24;
t13 = -mrSges(6,1) * t48 + mrSges(6,2) * t46;
t11 = Ifges(6,1) * t46 + Ifges(6,4) * t48 + Ifges(6,5) * t190;
t10 = Ifges(7,1) * t46 + Ifges(7,4) * t48 + Ifges(7,5) * t190;
t9 = Ifges(6,4) * t46 + Ifges(6,2) * t48 + Ifges(6,6) * t190;
t8 = Ifges(7,4) * t46 + Ifges(7,2) * t48 + Ifges(7,6) * t190;
t1 = [((mrSges(4,2) * t214) - t150 * t58 + t153 * t59 - 0.2e1 * t155 * t68 + (-t153 * t99 - t150 * t100 + t151 * (-t199 - t201)) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) + (-t172 - t200) * t154 - t226 * t104 - t225 * t102 + (-0.2e1 * m(5) * t155 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + Ifges(5,3) + t224) * t151) * qJD(3)) * t154 - (t10 + t11) * t104 - (t8 + t9) * t102 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t214 + (t54 + t55) * t48 + (t56 + t57) * t46 + (mrSges(4,1) * t214 + (-0.2e1 * qJ(2) * mrSges(4,2) - t153 * t100 + 0.2e1 * t155 * t111 + t150 * t99 + t151 * t172) * qJD(3) + t163 + t219) * t151 + (t17 * t2 + t18 * t3 + t19 * t71) * t215 + (t116 * t92 + t24 * t6 + t25 * t5) * t216 + (t97 * t52 + t96 * t53) * t217 + 0.2e1 * t17 * t26 + 0.2e1 * t24 * t27 + 0.2e1 * t18 * t28 + 0.2e1 * t25 * t29 + 0.2e1 * t19 * t62 + 0.2e1 * t71 * t12 + 0.2e1 * t3 * t88 + 0.2e1 * t5 * t89 + 0.2e1 * t2 * t90 + 0.2e1 * t6 * t91 + 0.2e1 * t92 * t63 + 0.2e1 * t96 * t94 + 0.2e1 * t97 * t95 + 0.2e1 * t116 * t13 + 0.2e1 * t52 * t125 + 0.2e1 * t53 * t126; t204 * t47 + t205 * t45 - t207 * t103 - (t26 + t27) * t101 + (-t12 - t13 - t68 + (t125 * t153 - t126 * t150) * qJD(3)) * t154 + (-t150 * t94 + t153 * t95 + (-t125 * t150 - t126 * t153) * qJD(4) + (t111 + t62 + t63) * qJD(3)) * t151 + m(6) * (-t101 * t6 - t103 * t5 + t116 * t191 - t154 * t92 + t47 * t24 + t45 * t25) + m(5) * ((-t150 * t96 + t153 * t97) * t190 + (-0.2e1 * t179 - t198 + t153 * t52 + (-t150 * t97 - t153 * t96) * qJD(4)) * t151) + m(7) * (-t101 * t2 - t103 * t3 - t154 * t19 + t47 * t17 + t45 * t18 + t191 * t71); 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (-t101 * t47 - t103 * t45 - t180) + (-0.1e1 + t192) * t180 * t217; -(t34 / 0.2e1 + t35 / 0.2e1) * t102 - (t54 / 0.2e1 + t55 / 0.2e1) * t78 - (t56 / 0.2e1 + t57 / 0.2e1) * t223 + (-t118 * t6 - t165 * t5 + t223 * t24 - t25 * t78) * mrSges(6,3) + (-t118 * t2 - t165 * t3 + t17 * t223 - t18 * t78) * mrSges(7,3) + (t146 / 0.2e1 - t75 / 0.2e1 - t73 / 0.2e1 - t76 / 0.2e1 - t74 / 0.2e1 + (t155 * t196 - Ifges(4,5)) * qJD(3)) * t151 + (t123 * t211 + t153 * t124 / 0.2e1 - t155 * t122 + (-t153 * t132 / 0.2e1 + t133 * t211) * qJD(4) + (-t155 * mrSges(4,2) + t201 / 0.2e1 + t199 / 0.2e1 - Ifges(4,6) + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t118 - (Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t165) * qJD(3)) * t154 - (t8 / 0.2e1 + t9 / 0.2e1) * t165 + (t82 / 0.2e1 + t83 / 0.2e1) * t48 + (t84 / 0.2e1 + t85 / 0.2e1) * t46 + m(5) * (-pkin(3) * t139 + (-t188 * t97 - t198) * pkin(8)) + m(6) * (t116 * t183 + t142 * t92 + t51 * t24 + t50 * t25 + t5 * t87 + t6 * t86) + (t10 / 0.2e1 + t11 / 0.2e1) * t118 + (-t133 * t191 / 0.2e1 + qJD(4) * t100 / 0.2e1 + t58 / 0.2e1 + t170 * mrSges(5,3) + (m(5) * t170 - qJD(4) * t126 + t95) * pkin(8)) * t153 + (t132 * t191 / 0.2e1 - t53 * mrSges(5,3) - pkin(8) * t94 + t59 / 0.2e1 + (-pkin(8) * t125 - t97 * mrSges(5,3) + pkin(4) * t63 - t99 / 0.2e1 - t197 / 0.2e1) * qJD(4)) * t150 + m(7) * (t15 * t18 + t16 * t17 + t19 * t98 + t60 * t2 + t61 * t3 + t65 * t71) + t60 * t26 + t61 * t28 + t65 * t62 - pkin(3) * t68 + t71 * t32 + t19 * t80 + t86 * t27 + t87 * t29 + t15 * t88 + t50 * t89 + t16 * t90 + t51 * t91 + t92 * t81 + t98 * t12 + t116 * t33 + t142 * t13 - (t36 / 0.2e1 + t37 / 0.2e1) * t104; (-t122 - t32 - t33) * t154 + m(6) * (-pkin(4) * t177 - t101 * t51 - t103 * t50 + t45 * t87 + t47 * t86) + m(7) * (-t101 * t16 - t103 * t15 - t154 * t65 + t61 * t45 + t60 * t47) + ((mrSges(5,3) * t192 - mrSges(4,2)) * t154 + m(5) * (t192 * t209 - t210) + (m(6) * t142 + m(7) * t98 + t196 + t80 + t81) * t151) * qJD(3) + (mrSges(7,3) + mrSges(6,3)) * (-t101 * t223 + t103 * t78 - t118 * t47 - t165 * t45); -0.2e1 * pkin(3) * t122 + t153 * t123 + t150 * t124 + 0.2e1 * t142 * t33 + 0.2e1 * t98 * t32 + 0.2e1 * t65 * t80 - (t82 + t83) * t78 - (t84 + t85) * t223 + (t36 + t37) * t118 - (t34 + t35) * t165 + (t153 * t133 + (0.2e1 * pkin(4) * t81 - t132) * t150) * qJD(4) + (t61 * t15 + t60 * t16 + t65 * t98) * t215 + (t142 * t183 + t50 * t87 + t51 * t86) * t216 + 0.2e1 * (-t118 * t16 - t15 * t165 + t223 * t60 - t61 * t78) * mrSges(7,3) + 0.2e1 * (-t118 * t51 - t165 * t50 + t223 * t86 - t78 * t87) * mrSges(6,3); -Ifges(5,6) * t176 - t162 * Ifges(5,5) + t182 * t141 + (t152 * t27 + t207 * t149 + (-t149 * t204 + t152 * t205) * qJD(5) + m(6) * (t149 * t5 + t152 * t6 + t184 * t25 - t185 * t24) + m(7) * (t149 * t3 - t17 * t185 + t18 * t184)) * pkin(4) + t157 - t52 * mrSges(5,2) + t53 * mrSges(5,1) + t219; (-t153 * t190 + t178) * mrSges(5,2) + (-t150 * t190 - t151 * t187) * mrSges(5,1) + m(6) * ((t152 * t47 + t175) * pkin(4) + t206) + m(7) * (pkin(4) * t175 + t141 * t47 + t206) + t164; t146 + t169 * t141 + (pkin(8) * t131 - t200) * qJD(4) + (m(6) * (t149 * t50 + t152 * t51 + t184 * t87 - t185 * t86) + m(7) * (t149 * t15 + t184 * t61 - t185 * t60) + t158 * mrSges(7,3) + (t152 * t223 + t158) * mrSges(6,3)) * pkin(4) + t159; 0.2e1 * (t173 + ((-t141 + t222) * m(7) - t227) * t149) * t220; pkin(5) * t182 + t157; t213 * t47 + t164; pkin(5) * t169 + t159; (t173 + (-t213 - t227) * t149) * t220; 0; m(7) * t19 + t12; m(7) * t191; m(7) * t65 + t32; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
