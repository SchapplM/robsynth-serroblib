% Calculate time derivative of joint inertia matrix for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:26
% EndTime: 2019-03-08 22:10:32
% DurationCPUTime: 2.74s
% Computational Cost: add. (2782->360), mult. (6672->512), div. (0->0), fcn. (5963->10), ass. (0->156)
t89 = sin(qJ(6));
t93 = cos(qJ(6));
t146 = t89 ^ 2 + t93 ^ 2;
t123 = t146 * mrSges(7,3);
t94 = cos(qJ(5));
t196 = (mrSges(6,2) - t123) * t94;
t64 = -t93 * mrSges(7,1) + mrSges(7,2) * t89;
t150 = mrSges(6,1) - t64;
t184 = -m(7) * pkin(5) - t150;
t87 = sin(pkin(6));
t92 = sin(qJ(2));
t156 = t87 * t92;
t88 = cos(pkin(6));
t91 = sin(qJ(3));
t95 = cos(qJ(3));
t45 = t91 * t156 - t88 * t95;
t46 = t95 * t156 + t88 * t91;
t90 = sin(qJ(5));
t108 = t45 * t94 - t46 * t90;
t140 = qJD(5) * t108;
t193 = m(7) * pkin(10);
t97 = -pkin(3) - pkin(4);
t61 = -t90 * qJ(4) + t94 * t97;
t42 = t94 * qJD(4) + qJD(5) * t61;
t192 = t42 * mrSges(6,2);
t191 = -t95 * pkin(3) - t91 * qJ(4);
t62 = t94 * qJ(4) + t90 * t97;
t190 = m(5) * pkin(8) + mrSges(5,2);
t65 = -t95 * mrSges(5,1) - t91 * mrSges(5,3);
t189 = -mrSges(4,1) * t95 + mrSges(4,2) * t91 + t65;
t188 = t146 * t94;
t187 = mrSges(7,3) + t193;
t129 = qJD(6) * t93;
t171 = pkin(8) - pkin(9);
t69 = t171 * t95;
t117 = qJD(3) * t69;
t127 = t171 * t91;
t38 = -t127 * t94 + t90 * t69;
t138 = qJD(5) * t38;
t143 = qJD(3) * t91;
t60 = t171 * t143;
t16 = t117 * t90 - t94 * t60 - t138;
t49 = t90 * t91 + t94 * t95;
t33 = (qJD(3) - qJD(5)) * t49;
t134 = qJD(5) * t94;
t142 = qJD(3) * t95;
t154 = t90 * t95;
t34 = -qJD(5) * t154 + t134 * t91 + t142 * t90 - t143 * t94;
t110 = mrSges(7,1) * t89 + mrSges(7,2) * t93;
t52 = t110 * qJD(6);
t66 = Ifges(7,5) * t89 + Ifges(7,6) * t93;
t130 = qJD(6) * t89;
t74 = Ifges(7,6) * t130;
t186 = -(Ifges(6,6) - t66 / 0.2e1) * t34 + Ifges(6,5) * t33 + t38 * t52 - t16 * mrSges(6,2) + t49 * (Ifges(7,5) * t129 - t74) / 0.2e1;
t56 = Ifges(7,4) * t129 - Ifges(7,2) * t130;
t57 = Ifges(7,1) * t129 - Ifges(7,4) * t130;
t168 = Ifges(7,4) * t89;
t67 = Ifges(7,2) * t93 + t168;
t167 = Ifges(7,4) * t93;
t68 = Ifges(7,1) * t89 + t167;
t98 = -(t89 * t67 - t93 * t68) * qJD(6) + t93 * t56 + t89 * t57;
t145 = qJD(2) * t92;
t126 = t87 * t145;
t96 = cos(qJ(2));
t155 = t87 * t96;
t25 = t45 * t90 + t46 * t94;
t19 = t89 * t155 + t25 * t93;
t144 = qJD(2) * t96;
t125 = t87 * t144;
t36 = qJD(3) * t46 + t125 * t91;
t37 = -qJD(3) * t45 + t125 * t95;
t8 = t36 * t90 + t37 * t94 + t140;
t3 = -qJD(6) * t19 - t126 * t93 - t8 * t89;
t18 = t93 * t155 - t25 * t89;
t4 = qJD(6) * t18 - t126 * t89 + t8 * t93;
t114 = t3 * t89 - t4 * t93;
t185 = (t18 * t93 + t19 * t89) * qJD(6) + t114;
t183 = 2 * m(6);
t182 = 0.2e1 * m(7);
t181 = -2 * mrSges(6,3);
t39 = t127 * t90 + t94 * t69;
t137 = qJD(5) * t39;
t17 = -t117 * t94 - t90 * t60 + t137;
t180 = 0.2e1 * t17;
t179 = 0.2e1 * t38;
t178 = -0.2e1 * t52;
t177 = 0.2e1 * t90;
t176 = m(6) / 0.2e1;
t175 = m(7) / 0.2e1;
t50 = t91 * t94 - t154;
t173 = -t50 / 0.2e1;
t172 = -t67 / 0.2e1;
t139 = qJD(5) * t25;
t7 = -t36 * t94 + t37 * t90 + t139;
t170 = t108 * t7;
t169 = mrSges(7,3) * t50;
t166 = Ifges(7,5) * t93;
t165 = t17 * t38;
t164 = t36 * t91;
t163 = t37 * t95;
t162 = t42 * t90;
t43 = t90 * qJD(4) + t62 * qJD(5);
t161 = t43 * t108;
t160 = t43 * t38;
t159 = t43 * t94;
t153 = t93 * t33;
t152 = t94 * t52;
t151 = qJD(6) / 0.2e1;
t149 = Ifges(4,4) - Ifges(5,5);
t148 = Ifges(7,5) * t153 + Ifges(7,3) * t34;
t147 = qJ(4) * t142 + t91 * qJD(4);
t136 = qJD(5) * t89;
t135 = qJD(5) * t93;
t63 = -pkin(2) + t191;
t48 = t95 * pkin(4) - t63;
t23 = t49 * pkin(5) - t50 * pkin(10) + t48;
t13 = t23 * t93 - t39 * t89;
t133 = qJD(6) * t13;
t14 = t23 * t89 + t39 * t93;
t132 = qJD(6) * t14;
t59 = -pkin(10) + t62;
t131 = qJD(6) * t59;
t124 = t50 * t130;
t122 = t146 * t42;
t121 = -Ifges(7,6) * t89 - (2 * Ifges(6,4));
t40 = t97 * t143 + t147;
t10 = t34 * pkin(5) - t33 * pkin(10) + t40;
t1 = t10 * t89 + t16 * t93 + t133;
t120 = -t1 + t133;
t2 = t10 * t93 - t16 * t89 - t132;
t119 = t2 + t132;
t118 = t87 ^ 2 * t92 * t144;
t104 = t124 - t153;
t105 = t129 * t50 + t89 * t33;
t5 = -Ifges(7,4) * t104 - Ifges(7,2) * t105 + Ifges(7,6) * t34;
t116 = t5 / 0.2e1 + t33 * t68 / 0.2e1;
t6 = -Ifges(7,1) * t104 - Ifges(7,4) * t105 + Ifges(7,5) * t34;
t115 = t6 / 0.2e1 + t33 * t172;
t113 = -t8 * mrSges(6,2) - t108 * t52;
t112 = -t108 * t17 + t38 * t7;
t58 = pkin(5) - t61;
t54 = (mrSges(4,1) * t91 + mrSges(4,2) * t95) * qJD(3);
t53 = (mrSges(5,1) * t91 - mrSges(5,3) * t95) * qJD(3);
t44 = pkin(3) * t143 - t147;
t35 = mrSges(6,1) * t49 + mrSges(6,2) * t50;
t31 = t49 * mrSges(7,1) - t93 * t169;
t30 = -t49 * mrSges(7,2) - t89 * t169;
t28 = pkin(8) * t163;
t26 = t110 * t50;
t21 = Ifges(7,5) * t49 + (Ifges(7,1) * t93 - t168) * t50;
t20 = Ifges(7,6) * t49 + (-Ifges(7,2) * t89 + t167) * t50;
t15 = mrSges(6,1) * t34 + mrSges(6,2) * t33;
t12 = -t34 * mrSges(7,2) - mrSges(7,3) * t105;
t11 = t34 * mrSges(7,1) + mrSges(7,3) * t104;
t9 = mrSges(7,1) * t105 - mrSges(7,2) * t104;
t22 = [0.2e1 * m(7) * (t18 * t3 + t19 * t4 - t170) + 0.2e1 * m(6) * (t25 * t8 - t118 - t170) + 0.2e1 * (m(5) + m(4)) * (t45 * t36 + t46 * t37 - t118); t18 * t11 + t19 * t12 - t108 * t9 + t7 * t26 + t3 * t31 + t4 * t30 + (-t108 * t33 - t25 * t34 - t49 * t8 + t50 * t7) * mrSges(6,3) + ((t15 - t53 - t54) * t96 + (-t96 * mrSges(3,2) + (-mrSges(3,1) - t35 + t189) * t92) * qJD(2)) * t87 + m(7) * (t1 * t19 + t13 * t3 + t14 * t4 + t18 * t2 + t112) + m(4) * (-pkin(2) * t126 + t28) + m(6) * (t16 * t25 + t39 * t8 + (-t145 * t48 + t40 * t96) * t87 + t112) + m(5) * (t126 * t63 - t44 * t155 + t28) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t142 * t45 - t143 * t46 + t164) * pkin(8) + (mrSges(4,3) + mrSges(5,2)) * (t164 + t163 + (t45 * t95 - t46 * t91) * qJD(3)); t39 * t34 * t181 - 0.2e1 * pkin(2) * t54 + 0.2e1 * t1 * t30 + 0.2e1 * t13 * t11 + 0.2e1 * t14 * t12 + 0.2e1 * t48 * t15 + t26 * t180 + 0.2e1 * t2 * t31 + 0.2e1 * t40 * t35 + t9 * t179 + 0.2e1 * t63 * t53 + 0.2e1 * (m(5) * t63 + t65) * t44 + (t1 * t14 + t13 * t2 + t165) * t182 + (t16 * t39 + t40 * t48 + t165) * t183 + 0.2e1 * t149 * qJD(3) * t95 ^ 2 + (mrSges(6,3) * t179 - t89 * t20 + t93 * t21) * t33 + (t16 * t181 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t34 + t121 * t33 + t148) * t49 + (mrSges(6,3) * t180 + 0.2e1 * Ifges(6,1) * t33 - t89 * t5 + t93 * t6 + (t121 + t166) * t34 + (-t93 * t20 - t89 * t21 - t49 * t66) * qJD(6)) * t50 + 0.2e1 * (-t149 * t91 + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t95) * t143; t150 * t7 + (-mrSges(4,2) + mrSges(5,3)) * t37 + (-mrSges(4,1) - mrSges(5,1)) * t36 + m(6) * (t25 * t42 - t61 * t7 + t62 * t8 - t161) + m(5) * (-pkin(3) * t36 + qJ(4) * t37 + qJD(4) * t46) + t185 * mrSges(7,3) - t113 + (-t185 * t59 + t58 * t7 - t161 + (-t18 * t89 + t19 * t93) * t42) * m(7); t43 * t26 + t58 * t9 + t150 * t17 + m(7) * (t58 * t17 + t160) + m(6) * (t16 * t62 - t17 * t61 + t39 * t42 + t160) + (-t33 * t61 - t34 * t62 - t42 * t49 + t43 * t50) * mrSges(6,3) + ((-mrSges(5,2) * pkin(3) + Ifges(5,4) + Ifges(4,5)) * t95 + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t91 + (m(5) * t191 + t189) * pkin(8)) * qJD(3) + (-t31 * t131 - qJD(6) * t21 / 0.2e1 + t42 * t30 + t59 * t12 + m(7) * (t1 * t59 - t13 * t131 + t14 * t42) + (t67 * t151 - t57 / 0.2e1) * t50 + t120 * mrSges(7,3) - t116) * t93 + (-t30 * t131 + t20 * t151 - t42 * t31 - t59 * t11 + m(7) * (-t13 * t42 - t131 * t14 - t2 * t59) + (t68 * t151 + t56 / 0.2e1) * t50 + t119 * mrSges(7,3) - t115) * t89 + t190 * qJD(4) * t95 - t186; t58 * t178 + 0.2e1 * t192 + (t122 * t59 + t58 * t43) * t182 + (t42 * t62 - t43 * t61) * t183 + t98 + 0.2e1 * t150 * t43 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4) - 0.2e1 * t123 * t42; m(5) * t36 + 0.2e1 * ((t135 * t19 - t136 * t18 - t7) * t175 + (-t7 + t139) * t176) * t94 + ((-t18 * t129 - t130 * t19 - t114 - t140) * t175 + (t8 - t140) * t176) * t177; t190 * t142 + (-t33 * mrSges(6,3) - t9 + (-t49 * mrSges(6,3) + t93 * t30 - t89 * t31) * qJD(5) + m(7) * (-t13 * t136 + t135 * t14 - t17) + m(6) * (-t17 + t137)) * t94 + (qJD(5) * t26 - t89 * t11 + t93 * t12 + (-t89 * t30 - t93 * t31) * qJD(6) + (qJD(5) * t50 - t34) * mrSges(6,3) + m(7) * (t1 * t93 - t129 * t13 - t130 * t14 - t2 * t89 + t138) + m(6) * (t16 + t138)) * t90; t152 + m(7) * (t146 * t162 - t159) + m(6) * (-t159 + t162) + (t150 * t90 + t196 + m(7) * (t188 * t59 + t58 * t90) + m(6) * (-t61 * t90 + t62 * t94)) * qJD(5); m(7) * (-0.1e1 + t146) * t134 * t177; t184 * t7 - t185 * t187 + t113; -pkin(5) * t9 + t184 * t17 + (t1 * mrSges(7,3) + t50 * t57 / 0.2e1 + (-t13 * mrSges(7,3) + t50 * t172 + t21 / 0.2e1) * qJD(6) + (-m(7) * t120 - qJD(6) * t31 + t12) * pkin(10) + t116) * t93 + (-t2 * mrSges(7,3) + t56 * t173 + (-t14 * mrSges(7,3) + t68 * t173 - t20 / 0.2e1) * qJD(6) + (-m(7) * t119 - qJD(6) * t30 - t11) * pkin(10) + t115) * t89 + t186; -t192 + (pkin(5) + t58) * t52 + t187 * t122 + t184 * t43 - t98; -t152 + (t184 * t90 + t188 * t193 - t196) * qJD(5); pkin(5) * t178 + t98; mrSges(7,1) * t3 - mrSges(7,2) * t4; t2 * mrSges(7,1) - t1 * mrSges(7,2) - Ifges(7,5) * t124 - Ifges(7,6) * t105 + t148; t74 - t110 * t42 + (t59 * t64 - t166) * qJD(6); (t130 * t90 - t134 * t93) * mrSges(7,2) + (-t129 * t90 - t134 * t89) * mrSges(7,1); -t74 + (pkin(10) * t64 + t166) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
