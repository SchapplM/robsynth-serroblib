% Calculate time derivative of joint inertia matrix for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:27
% EndTime: 2019-03-09 08:21:32
% DurationCPUTime: 2.62s
% Computational Cost: add. (1725->378), mult. (4073->547), div. (0->0), fcn. (3121->6), ass. (0->154)
t192 = pkin(4) + qJ(3);
t125 = sin(pkin(9));
t126 = cos(pkin(9));
t191 = qJD(3) * (t125 ^ 2 + t126 ^ 2);
t150 = qJD(4) * t125;
t92 = qJD(5) * t126 + t150;
t190 = 0.2e1 * t92;
t188 = Ifges(4,1) + Ifges(6,3);
t187 = -Ifges(6,1) - Ifges(5,3);
t172 = pkin(3) + qJ(5);
t186 = t172 * t126;
t127 = sin(qJ(6));
t129 = cos(qJ(6));
t83 = t125 * t129 + t126 * t127;
t69 = t83 * qJD(6);
t130 = cos(qJ(2));
t128 = sin(qJ(2));
t157 = t126 * t128;
t185 = pkin(4) * t157 + t130 * qJ(5);
t152 = qJD(2) * t128;
t184 = -qJ(4) * t152 + qJD(4) * t130;
t183 = 2 * m(4);
t182 = 2 * m(5);
t181 = 2 * m(6);
t180 = 2 * m(7);
t179 = -2 * pkin(1);
t178 = 2 * pkin(7);
t177 = t83 / 0.2e1;
t84 = t125 * t127 - t126 * t129;
t176 = t84 / 0.2e1;
t151 = qJD(2) * t130;
t140 = t126 * t151;
t161 = -qJ(4) * t140 - qJD(4) * t157;
t38 = (pkin(3) * t125 + pkin(7)) * t151 + t161;
t175 = m(5) * t38;
t171 = -pkin(5) - qJ(4);
t70 = t84 * qJD(6);
t170 = Ifges(7,5) * t69 - Ifges(7,6) * t70;
t169 = m(6) * qJD(3);
t168 = Ifges(4,4) * t125;
t167 = Ifges(4,4) * t126;
t166 = Ifges(6,5) * t125;
t165 = Ifges(6,5) * t126;
t164 = Ifges(5,6) * t125;
t163 = Ifges(5,6) * t126;
t68 = -t128 * qJD(3) + (pkin(2) * t128 - qJ(3) * t130) * qJD(2);
t162 = t126 * t68;
t156 = t126 * t130;
t93 = -pkin(2) * t130 - t128 * qJ(3) - pkin(1);
t56 = pkin(7) * t156 + t125 * t93;
t141 = t125 * t151;
t75 = mrSges(6,1) * t152 + mrSges(6,2) * t141;
t160 = qJ(5) * t125;
t159 = t125 * t128;
t158 = t125 * t130;
t64 = mrSges(4,1) * t141 + mrSges(4,2) * t140;
t155 = pkin(3) * t159 + t128 * pkin(7);
t154 = qJ(3) * t191;
t79 = mrSges(5,1) * t140 + mrSges(5,2) * t152;
t94 = t192 * t125;
t96 = t192 * t126;
t147 = Ifges(5,4) - Ifges(4,5) + Ifges(6,6);
t146 = Ifges(6,4) - Ifges(5,5) + Ifges(4,6);
t25 = -t128 * t69 - t84 * t151;
t26 = -t128 * t70 + t83 * t151;
t145 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t152;
t144 = pkin(7) * t152;
t143 = t125 * (-pkin(4) - pkin(8));
t142 = -pkin(7) * t125 - pkin(3);
t7 = -t25 * mrSges(7,1) + t26 * mrSges(7,2);
t32 = t70 * mrSges(7,1) + t69 * mrSges(7,2);
t139 = qJ(4) * t125 + pkin(2);
t106 = pkin(7) * t158;
t55 = t126 * t93 - t106;
t45 = qJ(4) * t130 - t56;
t138 = pkin(4) * t140 + t130 * qJD(5) - t162;
t122 = t130 * pkin(3);
t46 = t122 - t55;
t137 = t172 * t125 + pkin(7);
t57 = t125 * t68;
t136 = t57 - t184;
t134 = mrSges(6,1) * t126 - mrSges(6,3) * t125;
t133 = -mrSges(5,2) * t125 - mrSges(5,3) * t126;
t20 = t106 + t122 + (pkin(8) * t128 - t93) * t126 + t185;
t21 = t128 * t143 + t171 * t130 + t56;
t4 = t127 * t21 + t129 * t20;
t3 = -t127 * t20 + t129 * t21;
t80 = pkin(8) * t125 + t94;
t81 = pkin(8) * t126 + t96;
t31 = t127 * t81 + t129 * t80;
t30 = -t127 * t80 + t129 * t81;
t43 = -t126 * t144 + t57;
t132 = (-qJ(5) + t142) * t128;
t131 = qJD(5) * t159 + t161;
t77 = (-mrSges(6,2) * t156 - mrSges(6,3) * t128) * qJD(2);
t97 = mrSges(6,1) * t125 + mrSges(6,3) * t126;
t95 = mrSges(5,2) * t126 - mrSges(5,3) * t125;
t91 = -pkin(3) * t126 - t139;
t90 = mrSges(5,1) * t157 - mrSges(5,2) * t130;
t89 = mrSges(5,1) * t159 + mrSges(5,3) * t130;
t88 = -mrSges(6,2) * t157 + mrSges(6,3) * t130;
t87 = -mrSges(4,1) * t130 - mrSges(4,3) * t157;
t86 = -mrSges(6,1) * t130 + mrSges(6,2) * t159;
t85 = mrSges(4,2) * t130 - mrSges(4,3) * t159;
t78 = (mrSges(5,1) * t158 - mrSges(5,3) * t128) * qJD(2);
t76 = (mrSges(4,1) * t128 - mrSges(4,3) * t156) * qJD(2);
t74 = (-mrSges(4,2) * t128 - mrSges(4,3) * t158) * qJD(2);
t73 = t133 * t128;
t72 = t134 * t128;
t71 = t139 + t186;
t63 = t133 * t151;
t62 = t134 * t151;
t61 = t83 * t128;
t60 = t84 * t128;
t59 = -qJ(4) * t157 + t155;
t58 = t171 * t125 - pkin(2) - t186;
t54 = (t128 * Ifges(5,4) + (-t126 * Ifges(5,2) + t164) * t130) * qJD(2);
t53 = (t128 * Ifges(5,5) + (t125 * Ifges(5,3) - t163) * t130) * qJD(2);
t52 = (t128 * Ifges(4,5) + (t126 * Ifges(4,1) - t168) * t130) * qJD(2);
t51 = (-t128 * Ifges(6,4) + (t125 * Ifges(6,1) + t165) * t130) * qJD(2);
t50 = (t128 * Ifges(4,6) + (-t125 * Ifges(4,2) + t167) * t130) * qJD(2);
t49 = (-t128 * Ifges(6,6) + (t126 * Ifges(6,3) + t166) * t130) * qJD(2);
t48 = -mrSges(7,1) * t130 - t61 * mrSges(7,3);
t47 = mrSges(7,2) * t130 - t60 * mrSges(7,3);
t44 = (qJ(4) * t126 - t160) * t128 - t155;
t42 = t125 * t144 + t162;
t41 = Ifges(7,1) * t84 + Ifges(7,4) * t83;
t40 = Ifges(7,4) * t84 + Ifges(7,2) * t83;
t39 = -mrSges(7,1) * t83 + mrSges(7,2) * t84;
t37 = -pkin(4) * t159 - t45;
t36 = (t171 * t126 + t160) * t128 + t155;
t35 = t142 * t152 - t162;
t34 = Ifges(7,1) * t69 - Ifges(7,4) * t70;
t33 = Ifges(7,4) * t69 - Ifges(7,2) * t70;
t29 = t46 + t185;
t28 = -t43 + t184;
t27 = mrSges(7,1) * t60 + mrSges(7,2) * t61;
t19 = Ifges(7,1) * t61 - Ifges(7,4) * t60 - Ifges(7,5) * t130;
t18 = Ifges(7,4) * t61 - Ifges(7,2) * t60 - Ifges(7,6) * t130;
t17 = t137 * t151 + t131;
t16 = (-pkin(4) * t158 - pkin(7) * t157) * qJD(2) + t136;
t15 = mrSges(7,1) * t152 - mrSges(7,3) * t26;
t14 = -mrSges(7,2) * t152 + mrSges(7,3) * t25;
t13 = qJD(2) * t132 + t138;
t12 = -t84 * qJD(3) - qJD(6) * t31;
t11 = t83 * qJD(3) + qJD(6) * t30;
t10 = (-pkin(5) * t126 + t137) * t151 + t131;
t9 = ((-pkin(7) * t126 + pkin(5)) * t128 + t130 * t143) * qJD(2) + t136;
t8 = (pkin(8) * t156 + t132) * qJD(2) + t138;
t6 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t152;
t5 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t152;
t2 = -qJD(6) * t4 - t127 * t8 + t129 * t9;
t1 = qJD(6) * t3 + t127 * t9 + t129 * t8;
t22 = [(t64 * t178 + (t49 + t52 - t54) * t126 + (-t50 + t51 + t53) * t125) * t128 + (((mrSges(3,1) * t179) - 0.2e1 * Ifges(3,4) * t128 + Ifges(7,5) * t61 - Ifges(7,6) * t60 - t146 * t159 - t147 * t157) * t128 + ((mrSges(3,2) * t179) + (-(2 * Ifges(6,2)) - (2 * Ifges(5,1)) - (2 * Ifges(4,3)) - (2 * Ifges(3,2)) + (2 * Ifges(3,1)) + (pkin(7) ^ 2 * t183) - Ifges(7,3) + (mrSges(4,2) * t178 + (Ifges(5,2) + t188) * t126) * t126 + (mrSges(4,1) * t178 + (Ifges(4,2) - t187) * t125 + 0.2e1 * (-Ifges(4,4) + Ifges(6,5) - Ifges(5,6)) * t126) * t125) * t128 + 0.2e1 * (t146 * t125 + t147 * t126 + Ifges(3,4)) * t130) * t130) * qJD(2) - t130 * t145 + 0.2e1 * t43 * t85 + 0.2e1 * t16 * t86 + 0.2e1 * t42 * t87 + 0.2e1 * t13 * t88 + 0.2e1 * t28 * t89 + 0.2e1 * t35 * t90 - 0.2e1 * t17 * t72 + 0.2e1 * t38 * t73 + 0.2e1 * t56 * t74 + 0.2e1 * t37 * t75 + 0.2e1 * t55 * t76 + 0.2e1 * t29 * t77 + 0.2e1 * t45 * t78 + 0.2e1 * t46 * t79 - t60 * t5 + t61 * t6 + 0.2e1 * t44 * t62 + 0.2e1 * t59 * t63 + 0.2e1 * t36 * t7 + 0.2e1 * t1 * t47 + 0.2e1 * t2 * t48 + t25 * t18 + t26 * t19 + 0.2e1 * t10 * t27 + 0.2e1 * t4 * t14 + 0.2e1 * t3 * t15 + (t1 * t4 + t10 * t36 + t2 * t3) * t180 + (t13 * t29 + t16 * t37 - t17 * t44) * t181 + (t28 * t45 + t35 * t46 + t38 * t59) * t182 + (t55 * t42 + t56 * t43) * t183; (t1 * t83 - t2 * t84 - t3 * t69 - t4 * t70) * mrSges(7,3) + ((Ifges(7,5) * t176 + Ifges(7,6) * t177 - Ifges(3,6) + (pkin(7) * mrSges(3,2)) + (Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1 + Ifges(4,6) / 0.2e1) * t126 + (-Ifges(6,6) / 0.2e1 - Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t125) * t128 + (Ifges(3,5) - t125 * (Ifges(4,2) * t126 + t168) / 0.2e1 - t126 * (-Ifges(5,2) * t125 - t163) / 0.2e1 + (-m(4) * pkin(2) - mrSges(4,1) * t126 + mrSges(4,2) * t125 - mrSges(3,1)) * pkin(7) + (t187 * t126 - t164 + t166) * t125 / 0.2e1 + (t188 * t125 - t165 + t167) * t126 / 0.2e1) * t130) * qJD(2) + m(7) * (t1 * t31 + t10 * t58 + t11 * t4 + t12 * t3 + t2 * t30 - t36 * t92) + (t43 * mrSges(4,3) - t28 * mrSges(5,1) - t16 * mrSges(6,2) + t50 / 0.2e1 - t51 / 0.2e1 - t53 / 0.2e1 + (t74 - t78) * qJ(3) + (t85 + t86 - t89) * qJD(3) + m(5) * (-qJ(3) * t28 - qJD(3) * t45) + t37 * t169 + m(4) * (qJ(3) * t43 + qJD(3) * t56)) * t126 + (-qJD(4) * t73 - t42 * mrSges(4,3) + t35 * mrSges(5,1) - t13 * mrSges(6,2) + t52 / 0.2e1 - t54 / 0.2e1 + t49 / 0.2e1 + (-t76 + t79) * qJ(3) + (-t87 + t88 + t90) * qJD(3) + m(5) * (qJ(3) * t35 + qJD(3) * t46 - qJD(4) * t59) + t29 * t169 + m(4) * (-qJ(3) * t42 - qJD(3) * t55)) * t125 - t130 * t170 / 0.2e1 + (t63 + t175) * t91 + t94 * t77 + t38 * t95 + t96 * t75 - t17 * t97 + t69 * t19 / 0.2e1 - t70 * t18 / 0.2e1 + t71 * t62 + t58 * t7 - t60 * t33 / 0.2e1 + t61 * t34 / 0.2e1 - pkin(2) * t64 + t36 * t32 + t10 * t39 + t25 * t40 / 0.2e1 + t26 * t41 / 0.2e1 + t11 * t47 + t12 * t48 + t30 * t15 + t31 * t14 + (t72 - t27) * t92 + t5 * t177 + t6 * t176 + m(6) * (t13 * t94 + t16 * t96 - t17 * t71 + t44 * t92); -0.2e1 * t95 * t150 + 0.2e1 * t58 * t32 + t83 * t33 + t84 * t34 - t70 * t40 + t69 * t41 + (t11 * t31 + t12 * t30 - t58 * t92) * t180 + (t71 * t92 + (t125 * t94 + t126 * t96) * qJD(3)) * t181 + (-t91 * t150 + t154) * t182 + t154 * t183 + 0.2e1 * (mrSges(5,1) - mrSges(6,2) + mrSges(4,3)) * t191 + (-t39 + t97) * t190 + 0.2e1 * (t11 * t83 - t12 * t84 - t30 * t69 - t31 * t70) * mrSges(7,3); t175 + m(6) * t17 + m(7) * t10 + ((m(4) * pkin(7)) + (-mrSges(6,1) - mrSges(5,3)) * t126 + (-mrSges(5,2) + mrSges(6,3)) * t125) * t151 + t7 + t64; -m(5) * t150 + (-m(7) / 0.2e1 - m(6) / 0.2e1) * t190 + t32; 0; -t127 * t15 + t129 * t14 + (-t127 * t47 - t129 * t48) * qJD(6) + t77 + m(7) * (t1 * t129 - t127 * t2 + (-t127 * t4 - t129 * t3) * qJD(6)) + m(6) * t13 + m(5) * t35 + t79; m(7) * (t129 * t11 - t127 * t12 + (-t127 * t31 - t129 * t30) * qJD(6)) + (m(5) + m(6)) * t125 * qJD(3) + (t127 * t69 - t129 * t70 + (-t127 * t83 + t129 * t84) * qJD(6)) * mrSges(7,3); 0; 0; t127 * t14 + t129 * t15 + (-t127 * t48 + t129 * t47) * qJD(6) + m(7) * (t1 * t127 + t129 * t2 + (-t127 * t3 + t129 * t4) * qJD(6)) + m(6) * t16 + t75; m(7) * (t11 * t127 + t12 * t129 + (-t127 * t30 + t129 * t31) * qJD(6)) + t126 * t169 + (-t127 * t70 - t129 * t69 + (t127 * t84 + t129 * t83) * qJD(6)) * mrSges(7,3); 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t145; mrSges(7,1) * t12 - mrSges(7,2) * t11 + t170; 0; (-mrSges(7,1) * t129 + mrSges(7,2) * t127) * qJD(6); (-mrSges(7,1) * t127 - mrSges(7,2) * t129) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
