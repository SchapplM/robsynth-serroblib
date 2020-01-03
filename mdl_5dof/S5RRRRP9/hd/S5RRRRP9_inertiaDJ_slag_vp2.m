% Calculate time derivative of joint inertia matrix for
% S5RRRRP9
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:49
% EndTime: 2019-12-31 22:03:58
% DurationCPUTime: 3.44s
% Computational Cost: add. (2705->366), mult. (6654->529), div. (0->0), fcn. (5307->6), ass. (0->157)
t200 = Ifges(6,4) + Ifges(5,5);
t199 = Ifges(6,6) - Ifges(5,6);
t203 = -mrSges(5,1) - mrSges(6,1);
t202 = -Ifges(6,2) - Ifges(5,3);
t133 = sin(qJ(2));
t132 = sin(qJ(3));
t136 = cos(qJ(2));
t170 = qJD(2) * t136;
t156 = t132 * t170;
t135 = cos(qJ(3));
t167 = qJD(3) * t135;
t141 = t133 * t167 + t156;
t197 = qJD(3) + qJD(4);
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t96 = t131 * t132 - t134 * t135;
t66 = t197 * t96;
t97 = t131 * t135 + t132 * t134;
t67 = t197 * t97;
t201 = t199 * t67 - t200 * t66;
t88 = t97 * t133;
t155 = t135 * t170;
t171 = qJD(2) * t133;
t198 = -Ifges(4,5) * t155 - Ifges(4,3) * t171;
t168 = qJD(3) * t133;
t158 = t132 * t168;
t142 = t155 - t158;
t174 = t133 * t135;
t187 = pkin(6) * t132;
t107 = -pkin(2) * t136 - t133 * pkin(7) - pkin(1);
t95 = t135 * t107;
t60 = -pkin(8) * t174 + t95 + (-pkin(3) - t187) * t136;
t175 = t132 * t133;
t172 = t135 * t136;
t119 = pkin(6) * t172;
t85 = t132 * t107 + t119;
t74 = -pkin(8) * t175 + t85;
t183 = t131 * t60 + t134 * t74;
t105 = (pkin(2) * t133 - pkin(7) * t136) * qJD(2);
t176 = t135 * t105 + t171 * t187;
t19 = (pkin(3) * t133 - pkin(8) * t172) * qJD(2) + (-t119 + (pkin(8) * t133 - t107) * t132) * qJD(3) + t176;
t169 = qJD(3) * t132;
t46 = t132 * t105 + t107 * t167 + (-t135 * t171 - t136 * t169) * pkin(6);
t33 = -pkin(8) * t141 + t46;
t6 = -qJD(4) * t183 - t131 * t33 + t134 * t19;
t196 = 2 * m(4);
t195 = 2 * m(5);
t194 = 2 * m(6);
t193 = -0.2e1 * pkin(1);
t192 = 0.2e1 * pkin(6);
t191 = -pkin(8) - pkin(7);
t190 = -t132 / 0.2e1;
t189 = t135 / 0.2e1;
t89 = t96 * t133;
t79 = -mrSges(5,1) * t136 + t89 * mrSges(5,3);
t80 = mrSges(6,1) * t136 - t89 * mrSges(6,2);
t182 = -t79 + t80;
t181 = Ifges(4,4) * t132;
t180 = Ifges(4,4) * t135;
t179 = Ifges(4,6) * t132;
t178 = Ifges(4,6) * t136;
t149 = Ifges(4,1) * t135 - t181;
t87 = -Ifges(4,5) * t136 + t133 * t149;
t177 = t135 * t87;
t111 = Ifges(4,1) * t132 + t180;
t173 = t135 * t111;
t106 = pkin(3) * t175 + t133 * pkin(6);
t166 = qJD(4) * t131;
t165 = qJD(4) * t134;
t164 = pkin(3) * t169;
t163 = pkin(3) * t166;
t162 = pkin(3) * t165;
t161 = t191 * t132;
t112 = t191 * t135;
t75 = -t131 * t112 - t134 * t161;
t160 = t75 * t166;
t83 = t141 * pkin(3) + pkin(6) * t170;
t123 = -pkin(3) * t135 - pkin(2);
t159 = qJD(3) * t191;
t104 = t132 * t159;
t152 = t135 * t159;
t43 = -qJD(4) * t75 + t134 * t104 + t131 * t152;
t76 = -t134 * t112 + t131 * t161;
t44 = qJD(4) * t76 + t131 * t104 - t134 * t152;
t154 = t76 * t43 + t44 * t75;
t153 = (2 * Ifges(3,4)) + t179;
t39 = -t96 * t170 - t197 * t88;
t22 = -mrSges(6,1) * t171 + t39 * mrSges(6,2);
t151 = -mrSges(4,1) * t135 + mrSges(4,2) * t132;
t150 = mrSges(4,1) * t132 + mrSges(4,2) * t135;
t148 = -Ifges(4,2) * t132 + t180;
t110 = Ifges(4,2) * t135 + t181;
t147 = Ifges(4,5) * t132 + Ifges(4,6) * t135;
t24 = -t131 * t74 + t134 * t60;
t47 = -t85 * qJD(3) + t176;
t144 = -t132 * t47 + t135 * t46;
t40 = -t166 * t175 + (t174 * t197 + t156) * t134 + t142 * t131;
t143 = t202 * t171 - t199 * t40 - t200 * t39;
t5 = t131 * t19 + t134 * t33 + t60 * t165 - t166 * t74;
t140 = t203 * t44 + (-mrSges(5,2) + mrSges(6,3)) * t43 + t201;
t2 = qJ(5) * t171 - qJD(5) * t136 + t5;
t3 = -pkin(4) * t171 - t6;
t139 = t6 * mrSges(5,1) - t3 * mrSges(6,1) - t5 * mrSges(5,2) + t2 * mrSges(6,3) - t143;
t117 = qJD(5) + t162;
t137 = -mrSges(5,2) * t162 + t117 * mrSges(6,3) + t203 * t163;
t130 = qJD(5) * mrSges(6,3);
t127 = Ifges(4,5) * t167;
t122 = -pkin(3) * t134 - pkin(4);
t120 = pkin(3) * t131 + qJ(5);
t103 = -mrSges(4,1) * t136 - mrSges(4,3) * t174;
t102 = mrSges(4,2) * t136 - mrSges(4,3) * t175;
t101 = t149 * qJD(3);
t100 = t148 * qJD(3);
t99 = t150 * qJD(3);
t86 = t133 * t148 - t178;
t84 = -t136 * t187 + t95;
t82 = -mrSges(4,2) * t171 - mrSges(4,3) * t141;
t81 = mrSges(4,1) * t171 - mrSges(4,3) * t142;
t78 = mrSges(5,2) * t136 - t88 * mrSges(5,3);
t77 = -t88 * mrSges(6,2) - mrSges(6,3) * t136;
t73 = Ifges(5,1) * t97 - Ifges(5,4) * t96;
t72 = Ifges(6,1) * t97 + Ifges(6,5) * t96;
t71 = Ifges(5,4) * t97 - Ifges(5,2) * t96;
t70 = Ifges(6,5) * t97 + Ifges(6,3) * t96;
t69 = mrSges(5,1) * t96 + mrSges(5,2) * t97;
t68 = mrSges(6,1) * t96 - mrSges(6,3) * t97;
t59 = mrSges(4,1) * t141 + mrSges(4,2) * t142;
t58 = pkin(4) * t96 - qJ(5) * t97 + t123;
t55 = mrSges(5,1) * t88 - mrSges(5,2) * t89;
t54 = mrSges(6,1) * t88 + mrSges(6,3) * t89;
t53 = -t111 * t168 + (Ifges(4,5) * t133 + t136 * t149) * qJD(2);
t52 = -t110 * t168 + (Ifges(4,6) * t133 + t136 * t148) * qJD(2);
t51 = -Ifges(5,1) * t89 - Ifges(5,4) * t88 - Ifges(5,5) * t136;
t50 = -Ifges(6,1) * t89 - Ifges(6,4) * t136 + Ifges(6,5) * t88;
t49 = -Ifges(5,4) * t89 - Ifges(5,2) * t88 - Ifges(5,6) * t136;
t48 = -Ifges(6,5) * t89 - Ifges(6,6) * t136 + Ifges(6,3) * t88;
t45 = pkin(4) * t88 + qJ(5) * t89 + t106;
t32 = -Ifges(5,1) * t66 - Ifges(5,4) * t67;
t31 = -Ifges(6,1) * t66 + Ifges(6,5) * t67;
t30 = -Ifges(5,4) * t66 - Ifges(5,2) * t67;
t29 = -Ifges(6,5) * t66 + Ifges(6,3) * t67;
t28 = mrSges(5,1) * t67 - mrSges(5,2) * t66;
t27 = mrSges(6,1) * t67 + mrSges(6,3) * t66;
t23 = -mrSges(5,2) * t171 - mrSges(5,3) * t40;
t21 = mrSges(5,1) * t171 - mrSges(5,3) * t39;
t20 = -mrSges(6,2) * t40 + mrSges(6,3) * t171;
t18 = pkin(4) * t136 - t24;
t17 = -qJ(5) * t136 + t183;
t15 = pkin(4) * t67 + qJ(5) * t66 - qJD(5) * t97 + t164;
t13 = mrSges(5,1) * t40 + mrSges(5,2) * t39;
t12 = mrSges(6,1) * t40 - mrSges(6,3) * t39;
t11 = Ifges(5,1) * t39 - Ifges(5,4) * t40 + Ifges(5,5) * t171;
t10 = Ifges(6,1) * t39 + Ifges(6,4) * t171 + Ifges(6,5) * t40;
t9 = Ifges(5,4) * t39 - Ifges(5,2) * t40 + Ifges(5,6) * t171;
t8 = Ifges(6,5) * t39 + Ifges(6,6) * t171 + Ifges(6,3) * t40;
t7 = pkin(4) * t40 - qJ(5) * t39 + qJD(5) * t89 + t83;
t1 = [(t59 * t192 - t132 * t52 + t135 * t53 + (-t132 * t87 - t135 * t86 + t136 * t147) * qJD(3) + (mrSges(3,1) * t193 + (Ifges(4,5) * t135 - t153) * t133 - t200 * t89 + t199 * t88 + (pkin(6) ^ 2 * t196 + t150 * t192 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) + t202) * t136) * qJD(2)) * t133 + (t48 - t49) * t40 + (t50 + t51) * t39 - (t10 + t11) * t89 + (t8 - t9) * t88 + (t106 * t83 + t183 * t5 + t24 * t6) * t195 + 0.2e1 * t183 * t23 + (t17 * t2 + t18 * t3 + t45 * t7) * t194 + (t85 * t46 + t84 * t47) * t196 + ((mrSges(3,2) * t193 - t132 * t86 + t136 * t153 + t177) * qJD(2) + t143 + t198) * t136 + 0.2e1 * t47 * t103 + 0.2e1 * t106 * t13 + 0.2e1 * t46 * t102 + 0.2e1 * t2 * t77 + 0.2e1 * t5 * t78 + 0.2e1 * t6 * t79 + 0.2e1 * t3 * t80 + 0.2e1 * t83 * t55 + 0.2e1 * t84 * t81 + 0.2e1 * t85 * t82 + 0.2e1 * t7 * t54 + 0.2e1 * t45 * t12 + 0.2e1 * t18 * t22 + 0.2e1 * t24 * t21 + 0.2e1 * t17 * t20; -(t127 + t201) * t136 / 0.2e1 + (t23 + t20) * t76 + m(5) * (t106 * t164 + t123 * t83 + t183 * t43 - t24 * t44 + t5 * t76 - t6 * t75) + (t22 - t21) * t75 + t52 * t189 + (t77 + t78) * t43 + (t10 / 0.2e1 + t11 / 0.2e1) * t97 + (t8 / 0.2e1 - t9 / 0.2e1) * t96 + m(6) * (t15 * t45 + t17 * t43 + t18 * t44 + t2 * t76 + t3 * t75 + t58 * t7) + (t48 / 0.2e1 - t49 / 0.2e1) * t67 - (t50 / 0.2e1 + t51 / 0.2e1) * t66 - (t31 / 0.2e1 + t32 / 0.2e1) * t89 + (t29 / 0.2e1 - t30 / 0.2e1) * t88 + (-t17 * t67 - t18 * t66 - t2 * t96 + t3 * t97) * mrSges(6,2) + (-t183 * t67 + t24 * t66 - t5 * t96 - t6 * t97) * mrSges(5,3) + ((-t132 * t85 - t135 * t84) * qJD(3) + t144) * mrSges(4,3) + (t100 * t190 - Ifges(3,6) * qJD(2) + t101 * t189 + (t111 * t190 - t135 * t110 / 0.2e1) * qJD(3) + (qJD(2) * mrSges(3,2) + t99) * pkin(6) + (t199 * t96 + t200 * t97 + t147) * qJD(2) / 0.2e1) * t133 + (t73 / 0.2e1 + t72 / 0.2e1) * t39 + (-t132 * t81 + t135 * t82 + m(4) * (-t167 * t84 - t169 * t85 + t144) - t102 * t169 - t103 * t167) * pkin(7) + (t177 / 0.2e1 + (-t86 / 0.2e1 + t178 / 0.2e1 + pkin(3) * t55) * t132) * qJD(3) + (t70 / 0.2e1 - t71 / 0.2e1) * t40 + t182 * t44 + (Ifges(3,5) + t110 * t190 + t173 / 0.2e1 + (-m(4) * pkin(2) - mrSges(3,1) + t151) * pkin(6)) * t170 + t132 * t53 / 0.2e1 + t123 * t13 + t106 * t28 + t83 * t69 + t7 * t68 + t15 * t54 + t58 * t12 - pkin(2) * t59 + t45 * t27; -0.2e1 * pkin(2) * t99 + t135 * t100 + t132 * t101 + 0.2e1 * t123 * t28 + 0.2e1 * t15 * t68 + 0.2e1 * t58 * t27 + (t31 + t32) * t97 + (t29 - t30) * t96 + (t70 - t71) * t67 - (t72 + t73) * t66 + (t173 + (0.2e1 * pkin(3) * t69 - t110) * t132) * qJD(3) + (t15 * t58 + t154) * t194 + (t123 * t164 + t154) * t195 + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * (-t43 * t96 + t44 * t97 - t66 * t75 - t67 * t76); m(6) * (t117 * t17 + t120 * t2 + t122 * t3) + ((t21 + m(5) * t6 + (m(5) * t183 + t78) * qJD(4)) * t134 + (t23 + m(5) * t5 + (-m(5) * t24 + m(6) * t18 + t182) * qJD(4)) * t131) * pkin(3) - t141 * Ifges(4,6) + t139 + t117 * t77 + t120 * t20 + t122 * t22 - t46 * mrSges(4,2) + t47 * mrSges(4,1) - Ifges(4,5) * t158 - t198; m(6) * (t117 * t76 + t120 * t43 + t122 * t44) + t127 + (pkin(7) * t151 - t179) * qJD(3) + (-t117 * t96 - t120 * t67 - t122 * t66) * mrSges(6,2) + (t97 * mrSges(6,2) * t166 + m(6) * t160 + m(5) * (t131 * t43 - t134 * t44 + t165 * t76 + t160) + (-t131 * t67 + t134 * t66 + (t131 * t97 - t134 * t96) * qJD(4)) * mrSges(5,3)) * pkin(3) + t140; 0.2e1 * m(6) * (t117 * t120 + t122 * t163) + 0.2e1 * t137; m(6) * (-pkin(4) * t3 + qJ(5) * t2 + qJD(5) * t17) + t139 + qJD(5) * t77 - pkin(4) * t22 + qJ(5) * t20; m(6) * (-pkin(4) * t44 + qJ(5) * t43 + qJD(5) * t76) + (pkin(4) * t66 - qJ(5) * t67 - qJD(5) * t96) * mrSges(6,2) + t140; m(6) * (-pkin(4) * t163 + qJ(5) * t117 + qJD(5) * t120) + t130 + t137; 0.2e1 * m(6) * qJ(5) * qJD(5) + 0.2e1 * t130; m(6) * t3 + t22; m(6) * t44 - t66 * mrSges(6,2); m(6) * t163; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
