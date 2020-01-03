% Calculate time derivative of joint inertia matrix for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:41
% EndTime: 2019-12-31 21:59:50
% DurationCPUTime: 3.35s
% Computational Cost: add. (2758->357), mult. (6865->516), div. (0->0), fcn. (5613->6), ass. (0->157)
t201 = Ifges(5,5) + Ifges(6,5);
t200 = Ifges(5,6) + Ifges(6,6);
t205 = -Ifges(5,3) - Ifges(6,3);
t139 = sin(qJ(3));
t140 = sin(qJ(2));
t142 = cos(qJ(3));
t169 = qJD(3) * t142;
t143 = cos(qJ(2));
t172 = qJD(2) * t143;
t147 = t139 * t172 + t140 * t169;
t138 = sin(qJ(4));
t141 = cos(qJ(4));
t150 = t138 * t139 - t141 * t142;
t197 = qJD(3) + qJD(4);
t71 = t197 * t150;
t106 = t138 * t142 + t139 * t141;
t72 = t197 * t106;
t204 = -t200 * t72 - t201 * t71;
t203 = pkin(3) * t141;
t202 = -mrSges(5,1) - mrSges(6,1);
t199 = pkin(3) * qJD(4);
t95 = t150 * t140;
t191 = -pkin(8) - pkin(7);
t122 = t191 * t139;
t123 = t191 * t142;
t81 = t138 * t122 - t141 * t123;
t160 = t142 * t172;
t173 = qJD(2) * t140;
t198 = -Ifges(4,5) * t160 - Ifges(4,3) * t173;
t119 = -pkin(2) * t143 - t140 * pkin(7) - pkin(1);
t175 = t142 * t143;
t129 = pkin(6) * t175;
t90 = t139 * t119 + t129;
t196 = 2 * m(4);
t195 = 2 * m(5);
t194 = 2 * m(6);
t193 = -0.2e1 * pkin(1);
t192 = 0.2e1 * pkin(6);
t190 = -t139 / 0.2e1;
t189 = t142 / 0.2e1;
t187 = pkin(6) * t139;
t104 = t142 * t119;
t177 = t140 * t142;
t63 = -pkin(8) * t177 + t104 + (-pkin(3) - t187) * t143;
t178 = t139 * t140;
t79 = -pkin(8) * t178 + t90;
t27 = t138 * t63 + t141 * t79;
t183 = Ifges(4,4) * t139;
t182 = Ifges(4,4) * t142;
t181 = Ifges(4,6) * t139;
t180 = Ifges(4,6) * t143;
t154 = Ifges(4,1) * t142 - t183;
t93 = -Ifges(4,5) * t143 + t140 * t154;
t179 = t142 * t93;
t121 = Ifges(4,1) * t139 + t182;
t176 = t142 * t121;
t117 = (pkin(2) * t140 - pkin(7) * t143) * qJD(2);
t174 = t142 * t117 + t173 * t187;
t118 = pkin(3) * t178 + t140 * pkin(6);
t171 = qJD(3) * t139;
t170 = qJD(3) * t140;
t168 = qJD(4) * t138;
t167 = qJD(4) * t141;
t166 = pkin(3) * t171;
t41 = -t140 * t72 - t150 * t172;
t21 = (pkin(3) * t140 - pkin(8) * t175) * qJD(2) + (-t129 + (pkin(8) * t140 - t119) * t139) * qJD(3) + t174;
t46 = t119 * t169 + t139 * t117 + (-t142 * t173 - t143 * t171) * pkin(6);
t35 = -pkin(8) * t147 + t46;
t6 = -qJD(4) * t27 - t138 * t35 + t141 * t21;
t2 = pkin(4) * t173 - qJ(5) * t41 + qJD(5) * t95 + t6;
t22 = mrSges(6,1) * t173 - mrSges(6,3) * t41;
t165 = m(6) * t2 + t22;
t88 = pkin(3) * t147 + pkin(6) * t172;
t131 = -pkin(3) * t142 - pkin(2);
t164 = qJD(3) * t191;
t163 = t139 * t170;
t42 = -t106 * t172 + t197 * t95;
t11 = -t42 * mrSges(6,1) + t41 * mrSges(6,2);
t29 = t72 * mrSges(6,1) - t71 * mrSges(6,2);
t159 = (-mrSges(5,2) - mrSges(6,2)) * t141;
t158 = (2 * Ifges(3,4)) + t181;
t26 = -t138 * t79 + t141 * t63;
t80 = t141 * t122 + t123 * t138;
t115 = t139 * t164;
t116 = t142 * t164;
t45 = -t81 * qJD(4) - t115 * t138 + t141 * t116;
t15 = qJ(5) * t71 - qJD(5) * t106 + t45;
t157 = m(6) * t15 + t71 * mrSges(6,3);
t156 = -mrSges(4,1) * t142 + mrSges(4,2) * t139;
t155 = mrSges(4,1) * t139 + mrSges(4,2) * t142;
t153 = -Ifges(4,2) * t139 + t182;
t120 = Ifges(4,2) * t142 + t183;
t152 = Ifges(4,5) * t139 + Ifges(4,6) * t142;
t47 = -t90 * qJD(3) + t174;
t151 = -t139 * t47 + t142 * t46;
t149 = t205 * t173 - t200 * t42 - t201 * t41;
t5 = t138 * t21 + t141 * t35 + t63 * t167 - t168 * t79;
t44 = t141 * t115 + t138 * t116 + t122 * t167 + t123 * t168;
t148 = t160 - t163;
t14 = -qJ(5) * t72 - qJD(5) * t150 + t44;
t146 = t45 * mrSges(5,1) + t15 * mrSges(6,1) - t44 * mrSges(5,2) - t14 * mrSges(6,2) + t204;
t145 = -t138 * t72 + (t106 * t138 - t141 * t150) * qJD(4);
t94 = t106 * t140;
t3 = qJ(5) * t42 - qJD(5) * t94 + t5;
t144 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) - t149;
t135 = Ifges(4,5) * t169;
t130 = pkin(4) + t203;
t114 = -mrSges(4,1) * t143 - mrSges(4,3) * t177;
t113 = mrSges(4,2) * t143 - mrSges(4,3) * t178;
t112 = t154 * qJD(3);
t111 = t153 * qJD(3);
t110 = t155 * qJD(3);
t92 = t140 * t153 - t180;
t91 = pkin(4) * t150 + t131;
t89 = -t143 * t187 + t104;
t87 = -mrSges(4,2) * t173 - mrSges(4,3) * t147;
t86 = mrSges(4,1) * t173 - mrSges(4,3) * t148;
t85 = -mrSges(5,1) * t143 + t95 * mrSges(5,3);
t84 = -mrSges(6,1) * t143 + t95 * mrSges(6,3);
t83 = mrSges(5,2) * t143 - t94 * mrSges(5,3);
t82 = mrSges(6,2) * t143 - t94 * mrSges(6,3);
t78 = Ifges(5,1) * t106 - Ifges(5,4) * t150;
t77 = Ifges(6,1) * t106 - Ifges(6,4) * t150;
t76 = Ifges(5,4) * t106 - Ifges(5,2) * t150;
t75 = Ifges(6,4) * t106 - Ifges(6,2) * t150;
t74 = mrSges(5,1) * t150 + mrSges(5,2) * t106;
t73 = mrSges(6,1) * t150 + mrSges(6,2) * t106;
t69 = pkin(4) * t94 + t118;
t62 = mrSges(4,1) * t147 + mrSges(4,2) * t148;
t58 = pkin(4) * t72 + t166;
t57 = mrSges(5,1) * t94 - mrSges(5,2) * t95;
t56 = mrSges(6,1) * t94 - mrSges(6,2) * t95;
t55 = -qJ(5) * t150 + t81;
t54 = -qJ(5) * t106 + t80;
t53 = -t121 * t170 + (Ifges(4,5) * t140 + t143 * t154) * qJD(2);
t52 = -t120 * t170 + (Ifges(4,6) * t140 + t143 * t153) * qJD(2);
t51 = -Ifges(5,1) * t95 - Ifges(5,4) * t94 - Ifges(5,5) * t143;
t50 = -Ifges(6,1) * t95 - Ifges(6,4) * t94 - Ifges(6,5) * t143;
t49 = -Ifges(5,4) * t95 - Ifges(5,2) * t94 - Ifges(5,6) * t143;
t48 = -Ifges(6,4) * t95 - Ifges(6,2) * t94 - Ifges(6,6) * t143;
t34 = -Ifges(5,1) * t71 - Ifges(5,4) * t72;
t33 = -Ifges(6,1) * t71 - Ifges(6,4) * t72;
t32 = -Ifges(5,4) * t71 - Ifges(5,2) * t72;
t31 = -Ifges(6,4) * t71 - Ifges(6,2) * t72;
t30 = mrSges(5,1) * t72 - mrSges(5,2) * t71;
t25 = -mrSges(5,2) * t173 + mrSges(5,3) * t42;
t24 = -mrSges(6,2) * t173 + mrSges(6,3) * t42;
t23 = mrSges(5,1) * t173 - mrSges(5,3) * t41;
t18 = -pkin(4) * t42 + t88;
t17 = -qJ(5) * t94 + t27;
t16 = -pkin(4) * t143 + t95 * qJ(5) + t26;
t12 = -mrSges(5,1) * t42 + mrSges(5,2) * t41;
t10 = Ifges(5,1) * t41 + Ifges(5,4) * t42 + Ifges(5,5) * t173;
t9 = Ifges(6,1) * t41 + Ifges(6,4) * t42 + Ifges(6,5) * t173;
t8 = Ifges(5,4) * t41 + Ifges(5,2) * t42 + Ifges(5,6) * t173;
t7 = Ifges(6,4) * t41 + Ifges(6,2) * t42 + Ifges(6,6) * t173;
t1 = [-(t9 + t10) * t95 - (t7 + t8) * t94 + (t118 * t88 + t26 * t6 + t27 * t5) * t195 + (t90 * t46 + t89 * t47) * t196 + (t16 * t2 + t17 * t3 + t18 * t69) * t194 + (t48 + t49) * t42 + (t50 + t51) * t41 + (t62 * t192 - t139 * t52 + t142 * t53 + (-t139 * t93 - t142 * t92 + t143 * t152) * qJD(3) + (mrSges(3,1) * t193 + (Ifges(4,5) * t142 - t158) * t140 - t201 * t95 - t200 * t94 + (pkin(6) ^ 2 * t196 + t155 * t192 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) + t205) * t143) * qJD(2)) * t140 + ((mrSges(3,2) * t193 - t139 * t92 + t143 * t158 + t179) * qJD(2) + t149 + t198) * t143 + 0.2e1 * t16 * t22 + 0.2e1 * t17 * t24 + 0.2e1 * t26 * t23 + 0.2e1 * t27 * t25 + 0.2e1 * t18 * t56 + 0.2e1 * t69 * t11 + 0.2e1 * t3 * t82 + 0.2e1 * t5 * t83 + 0.2e1 * t2 * t84 + 0.2e1 * t6 * t85 + 0.2e1 * t88 * t57 + 0.2e1 * t89 * t86 + 0.2e1 * t90 * t87 + 0.2e1 * t46 * t113 + 0.2e1 * t47 * t114 + 0.2e1 * t118 * t12; (t112 * t189 + t111 * t190 - Ifges(3,6) * qJD(2) + (-t142 * t120 / 0.2e1 + t121 * t190) * qJD(3) + (qJD(2) * mrSges(3,2) + t110) * pkin(6) + (t201 * t106 - t150 * t200 + t152) * qJD(2) / 0.2e1) * t140 - (t7 / 0.2e1 + t8 / 0.2e1) * t150 + (-t2 * t106 - t150 * t3 + t16 * t71 - t17 * t72) * mrSges(6,3) + (-t106 * t6 - t150 * t5 + t26 * t71 - t27 * t72) * mrSges(5,3) + t52 * t189 + (t179 / 0.2e1 + (-t92 / 0.2e1 + pkin(3) * t57 + t180 / 0.2e1) * t139) * qJD(3) + (t142 * t87 + m(4) * (-t169 * t89 - t171 * t90 + t151) - t139 * t86 - t114 * t169 - t113 * t171) * pkin(7) + m(5) * (t118 * t166 + t131 * t88 + t26 * t45 + t27 * t44 + t5 * t81 + t6 * t80) - (t135 + t204) * t143 / 0.2e1 + ((-t139 * t90 - t142 * t89) * qJD(3) + t151) * mrSges(4,3) - (t48 / 0.2e1 + t49 / 0.2e1) * t72 - (t50 / 0.2e1 + t51 / 0.2e1) * t71 + (t9 / 0.2e1 + t10 / 0.2e1) * t106 + m(6) * (t14 * t17 + t15 * t16 + t18 * t91 + t2 * t54 + t3 * t55 + t58 * t69) - (t33 / 0.2e1 + t34 / 0.2e1) * t95 - (t31 / 0.2e1 + t32 / 0.2e1) * t94 + (Ifges(3,5) + t176 / 0.2e1 + t120 * t190 + (-m(4) * pkin(2) - mrSges(3,1) + t156) * pkin(6)) * t172 + t54 * t22 + t55 * t24 + t58 * t56 - pkin(2) * t62 + t69 * t29 + t18 * t73 + t80 * t23 + t81 * t25 + t14 * t82 + t44 * t83 + t15 * t84 + t45 * t85 + t88 * t74 + t91 * t11 + t118 * t30 + t131 * t12 + t139 * t53 / 0.2e1 + (t77 / 0.2e1 + t78 / 0.2e1) * t41 + (t75 / 0.2e1 + t76 / 0.2e1) * t42; -0.2e1 * pkin(2) * t110 + t142 * t111 + t139 * t112 + 0.2e1 * t131 * t30 + 0.2e1 * t91 * t29 + 0.2e1 * t58 * t73 - (t75 + t76) * t72 - (t77 + t78) * t71 + (t33 + t34) * t106 - (t31 + t32) * t150 + (t176 + (0.2e1 * pkin(3) * t74 - t120) * t139) * qJD(3) + (t131 * t166 + t44 * t81 + t45 * t80) * t195 + (t14 * t55 + t15 * t54 + t58 * t91) * t194 + 0.2e1 * (-t106 * t15 - t14 * t150 + t54 * t71 - t55 * t72) * mrSges(6,3) + 0.2e1 * (-t106 * t45 - t150 * t44 + t71 * t80 - t72 * t81) * mrSges(5,3); -Ifges(4,5) * t163 + t165 * t130 + (t141 * t23 + (t24 + t25) * t138 + ((t82 + t83) * t141 + (-t84 - t85) * t138) * qJD(4) + m(6) * (t138 * t3 - t16 * t168 + t167 * t17) + m(5) * (t138 * t5 + t141 * t6 + t167 * t27 - t168 * t26)) * pkin(3) + t144 - t147 * Ifges(4,6) - t46 * mrSges(4,2) + t47 * mrSges(4,1) - t198; t135 + t157 * t130 + (pkin(7) * t156 - t181) * qJD(3) + (m(6) * (t138 * t14 + t167 * t55 - t168 * t54) + m(5) * (t138 * t44 + t141 * t45 + t167 * t81 - t168 * t80) + t145 * mrSges(6,3) + (t141 * t71 + t145) * mrSges(5,3)) * pkin(3) + t146; 0.2e1 * (t159 + ((-t130 + t203) * m(6) + t202) * t138) * t199; pkin(4) * t165 + t144; pkin(4) * t157 + t146; (t159 + (-m(6) * pkin(4) + t202) * t138) * t199; 0; m(6) * t18 + t11; m(6) * t58 + t29; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
