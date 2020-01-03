% Calculate time derivative of joint inertia matrix for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:27
% EndTime: 2019-12-31 21:32:36
% DurationCPUTime: 3.78s
% Computational Cost: add. (1989->373), mult. (4739->547), div. (0->0), fcn. (3554->6), ass. (0->166)
t175 = Ifges(5,4) + Ifges(4,5);
t120 = cos(qJ(3));
t152 = qJD(3) * t120;
t117 = sin(qJ(3));
t154 = qJD(3) * t117;
t197 = Ifges(5,6) * t154 + t175 * t152;
t196 = -Ifges(5,2) - Ifges(4,3);
t118 = sin(qJ(2));
t121 = cos(qJ(2));
t155 = qJD(2) * t121;
t144 = t117 * t155;
t124 = t118 * t152 + t144;
t166 = Ifges(5,5) * t117;
t136 = Ifges(5,1) * t120 + t166;
t168 = Ifges(4,4) * t117;
t137 = Ifges(4,1) * t120 - t168;
t195 = (t136 + t137) * qJD(3);
t116 = sin(qJ(5));
t119 = cos(qJ(5));
t131 = t116 * t120 - t117 * t119;
t66 = t131 * t118;
t194 = -mrSges(6,1) * t116 - mrSges(6,2) * t119;
t179 = pkin(6) * t117;
t107 = t121 * t179;
t115 = t121 * pkin(3);
t93 = -pkin(2) * t121 - pkin(7) * t118 - pkin(1);
t34 = pkin(4) * t121 + t107 + t115 + (-pkin(8) * t118 - t93) * t120;
t160 = t117 * t118;
t158 = t120 * t121;
t108 = pkin(6) * t158;
t58 = t117 * t93 + t108;
t48 = -qJ(4) * t121 + t58;
t37 = pkin(8) * t160 + t48;
t6 = -t116 * t37 + t119 * t34;
t89 = (pkin(2) * t118 - pkin(7) * t121) * qJD(2);
t141 = qJD(3) * t108 - t120 * t89 + t93 * t154;
t153 = qJD(3) * t118;
t146 = t117 * t153;
t147 = -pkin(3) - t179;
t8 = pkin(8) * t146 + (-pkin(8) * t158 + (-pkin(4) + t147) * t118) * qJD(2) + t141;
t156 = qJD(2) * t118;
t109 = qJ(4) * t156;
t159 = t118 * t120;
t171 = t117 * t89 + t93 * t152;
t9 = t109 + (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t159 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t117) * t121 + t171;
t1 = qJD(5) * t6 + t116 * t8 + t119 * t9;
t7 = t116 * t34 + t119 * t37;
t2 = -qJD(5) * t7 - t116 * t9 + t119 * t8;
t193 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t192 = qJD(3) - qJD(5);
t162 = qJ(4) * t117;
t185 = pkin(3) + pkin(4);
t191 = -t120 * t185 - t162;
t133 = pkin(3) * t120 + t162;
t150 = qJD(4) * t120;
t190 = qJD(3) * t133 - t150;
t189 = 2 * m(4);
t188 = 2 * m(6);
t187 = -0.2e1 * pkin(1);
t186 = 0.2e1 * pkin(6);
t184 = pkin(7) - pkin(8);
t90 = -t116 * qJ(4) - t119 * t185;
t64 = t119 * qJD(4) + qJD(5) * t90;
t178 = t64 * mrSges(6,2);
t91 = t119 * qJ(4) - t116 * t185;
t65 = -t116 * qJD(4) - qJD(5) * t91;
t177 = t65 * mrSges(6,1);
t176 = qJD(2) / 0.2e1;
t130 = t116 * t117 + t119 * t120;
t40 = t192 * t130;
t18 = t118 * t40 - t131 * t155;
t19 = t130 * t155 + t192 * t66;
t174 = Ifges(6,5) * t19 + Ifges(6,6) * t18;
t41 = t192 * t131;
t173 = Ifges(6,5) * t40 - Ifges(6,6) * t41;
t61 = -Ifges(5,4) * t121 + t118 * t136;
t62 = -Ifges(4,5) * t121 + t118 * t137;
t172 = t61 + t62;
t167 = Ifges(4,4) * t120;
t165 = Ifges(5,5) * t120;
t164 = Ifges(4,6) * t120;
t163 = t121 * Ifges(4,6);
t161 = qJ(4) * t120;
t151 = qJD(4) * t117;
t95 = -Ifges(5,3) * t120 + t166;
t96 = Ifges(4,2) * t120 + t168;
t149 = t95 / 0.2e1 - t96 / 0.2e1;
t97 = Ifges(5,1) * t117 - t165;
t98 = Ifges(4,1) * t117 + t167;
t148 = t97 / 0.2e1 + t98 / 0.2e1;
t100 = t184 * t120;
t143 = t120 * t155;
t57 = t120 * t93 - t107;
t134 = Ifges(5,3) * t117 + t165;
t59 = -Ifges(5,6) * t121 + t118 * t134;
t135 = -Ifges(4,2) * t117 + t167;
t60 = t118 * t135 - t163;
t142 = t59 - t60 + t163;
t140 = -t120 * mrSges(4,1) + t117 * mrSges(4,2);
t139 = mrSges(4,1) * t117 + mrSges(4,2) * t120;
t94 = -t120 * mrSges(5,1) - t117 * mrSges(5,3);
t138 = mrSges(5,1) * t117 - mrSges(5,3) * t120;
t132 = pkin(3) * t117 - t161;
t99 = t184 * t117;
t47 = t100 * t119 + t116 * t99;
t46 = -t100 * t116 + t119 * t99;
t129 = pkin(6) + t132;
t128 = -t117 * t185 + t161;
t87 = t184 * t154;
t88 = qJD(3) * t100;
t20 = qJD(5) * t46 + t116 * t88 - t119 * t87;
t21 = -qJD(5) * t47 + t116 * t87 + t119 * t88;
t127 = t21 * mrSges(6,1) - t20 * mrSges(6,2) + t173;
t126 = -pkin(6) + t128;
t125 = t143 - t146;
t123 = -t124 * Ifges(5,6) - t175 * t143 + t196 * t156 + t174;
t53 = mrSges(5,2) * t143 + (-mrSges(5,1) * qJD(2) - mrSges(5,2) * t154) * t118;
t25 = (-t120 * t156 - t121 * t154) * pkin(6) + t171;
t92 = -pkin(2) - t133;
t86 = -mrSges(5,2) * t160 - mrSges(5,3) * t121;
t85 = mrSges(5,1) * t121 + mrSges(5,2) * t159;
t84 = -mrSges(4,1) * t121 - mrSges(4,3) * t159;
t83 = mrSges(4,2) * t121 - mrSges(4,3) * t160;
t80 = t135 * qJD(3);
t79 = t134 * qJD(3);
t78 = t139 * qJD(3);
t77 = t138 * qJD(3);
t73 = pkin(2) - t191;
t70 = t138 * t118;
t68 = qJD(3) * t132 - t151;
t67 = t130 * t118;
t63 = t129 * t118;
t56 = qJD(3) * t128 + t151;
t55 = -mrSges(5,2) * t124 + mrSges(5,3) * t156;
t54 = -mrSges(4,2) * t156 - mrSges(4,3) * t124;
t52 = mrSges(4,1) * t156 - mrSges(4,3) * t125;
t51 = mrSges(6,1) * t121 - mrSges(6,3) * t67;
t50 = -mrSges(6,2) * t121 - mrSges(6,3) * t66;
t49 = t115 - t57;
t45 = t126 * t118;
t44 = -Ifges(6,1) * t131 - Ifges(6,4) * t130;
t43 = -Ifges(6,4) * t131 - Ifges(6,2) * t130;
t42 = mrSges(6,1) * t130 - mrSges(6,2) * t131;
t36 = mrSges(4,1) * t124 + mrSges(4,2) * t125;
t35 = mrSges(5,1) * t124 - mrSges(5,3) * t125;
t33 = mrSges(6,1) * t66 + mrSges(6,2) * t67;
t32 = -t98 * t153 + (Ifges(4,5) * t118 + t121 * t137) * qJD(2);
t31 = -t97 * t153 + (Ifges(5,4) * t118 + t121 * t136) * qJD(2);
t30 = -t96 * t153 + (Ifges(4,6) * t118 + t121 * t135) * qJD(2);
t29 = -t95 * t153 + (Ifges(5,6) * t118 + t121 * t134) * qJD(2);
t28 = Ifges(6,1) * t67 - Ifges(6,4) * t66 + Ifges(6,5) * t121;
t27 = Ifges(6,4) * t67 - Ifges(6,2) * t66 + Ifges(6,6) * t121;
t26 = t156 * t179 - t141;
t24 = t118 * t190 + t129 * t155;
t23 = t147 * t156 + t141;
t22 = -qJD(4) * t121 + t109 + t25;
t15 = (qJD(3) * t191 + t150) * t118 + t126 * t155;
t14 = Ifges(6,1) * t40 - Ifges(6,4) * t41;
t13 = Ifges(6,4) * t40 - Ifges(6,2) * t41;
t12 = mrSges(6,1) * t41 + mrSges(6,2) * t40;
t11 = -mrSges(6,1) * t156 - mrSges(6,3) * t19;
t10 = mrSges(6,2) * t156 + mrSges(6,3) * t18;
t5 = -mrSges(6,1) * t18 + mrSges(6,2) * t19;
t4 = Ifges(6,1) * t19 + Ifges(6,4) * t18 - Ifges(6,5) * t156;
t3 = Ifges(6,4) * t19 + Ifges(6,2) * t18 - Ifges(6,6) * t156;
t16 = [0.2e1 * t25 * t83 + 0.2e1 * t26 * t84 + 0.2e1 * t23 * t85 + 0.2e1 * t22 * t86 - t66 * t3 + t67 * t4 + 0.2e1 * t24 * t70 + 0.2e1 * t1 * t50 + 0.2e1 * t2 * t51 + 0.2e1 * t49 * t53 + 0.2e1 * t48 * t55 + 0.2e1 * t57 * t52 + 0.2e1 * t58 * t54 + 0.2e1 * t63 * t35 + 0.2e1 * t45 * t5 + t18 * t27 + t19 * t28 + 0.2e1 * t15 * t33 + 0.2e1 * t7 * t10 + 0.2e1 * t6 * t11 + (t36 * t186 + (t31 + t32) * t120 + (t29 - t30) * t117 + (t142 * t120 + (t121 * t175 - t172) * t117) * qJD(3) + (mrSges(3,1) * t187 - Ifges(6,5) * t67 + Ifges(6,6) * t66 + (-(2 * Ifges(3,4)) + t175 * t120 + (-Ifges(4,6) + Ifges(5,6)) * t117) * t118 + (pkin(6) ^ 2 * t189 + t139 * t186 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(6,3)) + t196) * t121) * qJD(2)) * t118 + 0.2e1 * m(5) * (t22 * t48 + t23 * t49 + t24 * t63) + (t25 * t58 + t26 * t57) * t189 + (t1 * t7 + t15 * t45 + t2 * t6) * t188 + ((mrSges(3,2) * t187 + 0.2e1 * Ifges(3,4) * t121 + t117 * t142 + t120 * t172) * qJD(2) + t123) * t121; m(6) * (t1 * t47 + t15 * t73 + t2 * t46 + t20 * t7 + t21 * t6 + t45 * t56) + t92 * t35 + t24 * t94 + t63 * t77 - t66 * t13 / 0.2e1 + t67 * t14 / 0.2e1 + t68 * t70 + t73 * t5 + t46 * t11 + t47 * t10 + t20 * t50 + t21 * t51 + t56 * t33 - pkin(2) * t36 + t40 * t28 / 0.2e1 - t41 * t27 / 0.2e1 + t15 * t42 + t18 * t43 / 0.2e1 + t19 * t44 / 0.2e1 + t45 * t12 + (-t29 / 0.2e1 + t30 / 0.2e1 + t25 * mrSges(4,3) + t22 * mrSges(5,2) + (t49 * mrSges(5,2) - t57 * mrSges(4,3) + t61 / 0.2e1 + t62 / 0.2e1) * qJD(3)) * t120 + (-t197 / 0.2e1 + t173 / 0.2e1) * t121 + m(5) * (t24 * t92 + t63 * t68) + (Ifges(3,5) + t148 * t120 + t149 * t117 + (-m(4) * pkin(2) - mrSges(3,1) + t140) * pkin(6)) * t155 + (t164 * t176 - qJD(2) * (-Ifges(6,5) * t131 - Ifges(6,6) * t130) / 0.2e1 - Ifges(3,6) * qJD(2) + (mrSges(3,2) * qJD(2) + t78) * pkin(6) + (t79 / 0.2e1 - t80 / 0.2e1 - t148 * qJD(3) + t175 * t176) * t117 + (t195 / 0.2e1 - Ifges(5,6) * t176 + t149 * qJD(3)) * t120) * t118 + (-t1 * t130 + t131 * t2 - t40 * t6 - t41 * t7) * mrSges(6,3) - t130 * t3 / 0.2e1 - t131 * t4 / 0.2e1 + ((t54 + t55) * t120 + (-t52 + t53) * t117 + ((-t84 + t85) * t120 + (-t83 - t86) * t117) * qJD(3) + m(4) * (-t26 * t117 + t25 * t120 - t152 * t57 - t154 * t58) + m(5) * (t23 * t117 + t22 * t120 + t152 * t49 - t154 * t48)) * pkin(7) + (t31 / 0.2e1 + t32 / 0.2e1 - t26 * mrSges(4,3) + t23 * mrSges(5,2) + (-t48 * mrSges(5,2) - t58 * mrSges(4,3) + t163 / 0.2e1 + t59 / 0.2e1 - t60 / 0.2e1) * qJD(3)) * t117; 0.2e1 * t56 * t42 + 0.2e1 * t73 * t12 + (t20 * t47 + t21 * t46 + t56 * t73) * t188 + t40 * t44 - t131 * t14 - t41 * t43 - t130 * t13 + 0.2e1 * t92 * t77 - 0.2e1 * pkin(2) * t78 + (-t79 + t80) * t120 + t195 * t117 + ((t97 + t98) * t120 + (t95 - t96) * t117) * qJD(3) + 0.2e1 * (m(5) * t92 + t94) * t68 + 0.2e1 * (-t130 * t20 + t131 * t21 - t40 * t46 - t41 * t47) * mrSges(6,3); -Ifges(4,6) * t144 + t90 * t11 + t91 * t10 + qJD(4) * t86 + t64 * t50 + t65 * t51 - pkin(3) * t53 + qJ(4) * t55 - t23 * mrSges(5,1) - t25 * mrSges(4,2) + t26 * mrSges(4,1) + t22 * mrSges(5,3) - t123 + (Ifges(6,3) * qJD(2) + (-t117 * t175 - t164) * qJD(3)) * t118 + m(6) * (t1 * t91 + t2 * t90 + t6 * t65 + t64 * t7) + m(5) * (-pkin(3) * t23 + qJ(4) * t22 + qJD(4) * t48) + t193; m(6) * (t20 * t91 + t21 * t90 + t46 * t65 + t47 * t64) - Ifges(4,6) * t154 - t190 * mrSges(5,2) + (-t130 * t64 + t131 * t65 - t40 * t90 - t41 * t91) * mrSges(6,3) + (m(5) * t150 + (-m(5) * t133 + t140 + t94) * qJD(3)) * pkin(7) - t127 + t197; (t64 * t91 + t65 * t90) * t188 + 0.2e1 * t178 - 0.2e1 * t177 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); t116 * t10 + t119 * t11 + (-t116 * t51 + t119 * t50) * qJD(5) + m(6) * (t1 * t116 + t119 * t2 + (-t116 * t6 + t119 * t7) * qJD(5)) + m(5) * t23 + t53; m(6) * (t116 * t20 + t119 * t21 + (-t116 * t46 + t119 * t47) * qJD(5)) + (m(5) * pkin(7) + mrSges(5,2)) * t152 + (-t116 * t41 - t119 * t40 + (-t116 * t131 - t119 * t130) * qJD(5)) * mrSges(6,3); m(6) * (t116 * t64 + t119 * t65) + (m(6) * (-t116 * t90 + t119 * t91) - t194) * qJD(5); 0; -Ifges(6,3) * t156 + t174 - t193; t127; t177 - t178; t194 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
