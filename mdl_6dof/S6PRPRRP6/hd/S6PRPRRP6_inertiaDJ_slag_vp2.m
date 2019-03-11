% Calculate time derivative of joint inertia matrix for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:16
% EndTime: 2019-03-08 20:18:20
% DurationCPUTime: 1.98s
% Computational Cost: add. (1353->318), mult. (3469->453), div. (0->0), fcn. (2665->8), ass. (0->152)
t176 = m(6) / 0.2e1;
t184 = m(7) / 0.2e1 + t176;
t194 = 0.2e1 * t184;
t193 = Ifges(7,2) + Ifges(6,3);
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t146 = t88 ^ 2 + t91 ^ 2;
t89 = sin(qJ(4));
t93 = -pkin(2) - pkin(8);
t162 = t89 * t93;
t92 = cos(qJ(4));
t58 = pkin(4) * t89 - pkin(9) * t92 + qJ(3);
t187 = t91 * t162 + t88 * t58;
t192 = qJD(5) * t187;
t191 = -0.2e1 * m(6);
t190 = m(6) + m(7);
t189 = mrSges(7,2) + mrSges(6,3);
t188 = t89 * mrSges(5,1) + mrSges(5,2) * t92 + mrSges(4,3);
t158 = Ifges(7,4) + Ifges(6,5);
t186 = t158 * t91;
t172 = cos(qJ(2));
t86 = sin(pkin(6));
t123 = t86 * t172;
t87 = cos(pkin(6));
t102 = -t123 * t92 - t87 * t89;
t145 = qJD(4) * t89;
t118 = t89 * t123;
t90 = sin(qJ(2));
t165 = t86 * t90;
t131 = qJD(2) * t165;
t144 = qJD(4) * t92;
t20 = -qJD(4) * t118 - t131 * t92 + t144 * t87;
t166 = t20 * t92;
t185 = -t102 * t145 - t166;
t121 = qJD(2) * t172;
t183 = m(7) * qJ(6) + mrSges(7,3);
t19 = qJD(4) * t102 + t131 * t89;
t41 = t87 * t92 - t118;
t182 = qJD(4) * (-t102 * t89 + t41 * t92) + t19 * t89 - t166;
t46 = qJD(3) + (pkin(4) * t92 + pkin(9) * t89) * qJD(4);
t181 = -t91 * t46 + t192;
t110 = pkin(5) * t91 + qJ(6) * t88;
t139 = qJD(6) * t91;
t180 = qJD(5) * t110 - t139;
t159 = t91 * t92;
t55 = mrSges(6,1) * t89 - mrSges(6,3) * t159;
t56 = -mrSges(7,1) * t89 + mrSges(7,2) * t159;
t150 = -t55 + t56;
t140 = qJD(5) * t92;
t127 = t91 * t140;
t130 = t88 * t145;
t100 = -t127 + t130;
t28 = -mrSges(6,2) * t144 + mrSges(6,3) * t100;
t29 = mrSges(7,2) * t100 + mrSges(7,3) * t144;
t155 = t28 + t29;
t122 = t88 * t93 - pkin(5);
t160 = t91 * t58;
t25 = t122 * t89 - t160;
t30 = -t162 * t88 + t160;
t142 = qJD(5) * t88;
t128 = t93 * t142;
t129 = t93 * t144;
t141 = qJD(5) * t91;
t133 = t91 * t129 + t58 * t141 + t88 * t46;
t6 = qJ(6) * t144 + (qJD(6) - t128) * t89 + t133;
t8 = -t128 * t89 + t133;
t179 = t150 * qJD(5) + m(7) * (t25 * qJD(5) + t6) + m(6) * (-t30 * qJD(5) + t8) + t155;
t164 = t88 * t92;
t54 = -mrSges(6,2) * t89 - mrSges(6,3) * t164;
t57 = -mrSges(7,2) * t164 + mrSges(7,3) * t89;
t151 = t54 + t57;
t101 = t88 * t140 + t145 * t91;
t26 = mrSges(6,1) * t144 + mrSges(6,3) * t101;
t27 = -mrSges(7,1) * t144 - mrSges(7,2) * t101;
t156 = -t26 + t27;
t24 = qJ(6) * t89 + t187;
t7 = t122 * t144 + t181;
t9 = -t129 * t88 - t181;
t178 = -t151 * qJD(5) + m(7) * (-t24 * qJD(5) + t7) + m(6) * (-t9 - t192) + t156;
t177 = 0.2e1 * t93;
t117 = t86 * t121;
t22 = t165 * t88 + t41 * t91;
t143 = qJD(5) * t22;
t4 = -t117 * t91 + t19 * t88 + t143;
t174 = t4 * t88;
t105 = t165 * t91 - t41 * t88;
t5 = qJD(5) * t105 + t117 * t88 + t19 * t91;
t173 = t5 * t91;
t171 = Ifges(6,4) * t88;
t170 = Ifges(6,4) * t91;
t169 = Ifges(7,5) * t88;
t168 = Ifges(7,5) * t91;
t167 = Ifges(7,6) * t89;
t11 = t102 * t20;
t163 = t89 * Ifges(6,6);
t17 = -mrSges(7,1) * t100 + mrSges(7,3) * t101;
t18 = -mrSges(6,1) * t100 - mrSges(6,2) * t101;
t157 = t17 + t18;
t111 = Ifges(7,3) * t88 + t168;
t34 = t111 * t92 + t167;
t112 = -Ifges(6,2) * t88 + t170;
t35 = t112 * t92 + t163;
t154 = t34 - t35;
t115 = t88 * mrSges(7,1) - t91 * mrSges(7,3);
t43 = t115 * t92;
t116 = t88 * mrSges(6,1) + t91 * mrSges(6,2);
t44 = t116 * t92;
t153 = t43 + t44;
t47 = t115 * qJD(5);
t48 = t116 * qJD(5);
t152 = t47 + t48;
t65 = -t91 * mrSges(6,1) + t88 * mrSges(6,2);
t149 = t65 - mrSges(5,1);
t148 = qJ(3) * t117 + qJD(3) * t165;
t147 = t146 * pkin(9) * t144;
t64 = -t91 * mrSges(7,1) - t88 * mrSges(7,3);
t134 = t64 + t149;
t126 = t105 * t141;
t67 = -Ifges(7,3) * t91 + t169;
t68 = Ifges(6,2) * t91 + t171;
t125 = t67 / 0.2e1 - t68 / 0.2e1;
t69 = Ifges(7,1) * t88 - t168;
t70 = Ifges(6,1) * t88 + t170;
t124 = -t69 / 0.2e1 - t70 / 0.2e1;
t120 = Ifges(6,6) * t130 + Ifges(7,6) * t127 + t193 * t144;
t114 = Ifges(6,1) * t91 - t171;
t113 = Ifges(7,1) * t91 + t169;
t109 = pkin(5) * t88 - qJ(6) * t91;
t36 = Ifges(7,4) * t89 + t113 * t92;
t37 = Ifges(6,5) * t89 + t114 * t92;
t106 = -t158 * t89 - t36 - t37;
t103 = t109 - t93;
t95 = m(5) * t182;
t94 = m(7) * t139 + (-m(7) * t110 + t64 + t65) * qJD(5);
t83 = Ifges(7,4) * t141;
t82 = Ifges(6,5) * t141;
t80 = Ifges(7,6) * t142;
t59 = -pkin(4) - t110;
t53 = t114 * qJD(5);
t52 = t113 * qJD(5);
t51 = t112 * qJD(5);
t50 = t111 * qJD(5);
t49 = (mrSges(5,1) * t92 - mrSges(5,2) * t89) * qJD(4);
t38 = qJD(5) * t109 - qJD(6) * t88;
t33 = t103 * t92;
t15 = -t70 * t140 + (Ifges(6,5) * t92 - t114 * t89) * qJD(4);
t14 = -t69 * t140 + (Ifges(7,4) * t92 - t113 * t89) * qJD(4);
t13 = -t68 * t140 + (Ifges(6,6) * t92 - t112 * t89) * qJD(4);
t12 = -t67 * t140 + (Ifges(7,6) * t92 - t111 * t89) * qJD(4);
t10 = -t103 * t145 + t180 * t92;
t3 = pkin(9) * t173;
t1 = [0.2e1 * m(5) * (t121 * t86 ^ 2 * t90 + t41 * t19 - t11) + 0.2e1 * t190 * (-t105 * t4 + t22 * t5 - t11); t151 * t5 - t157 * t102 + t150 * t4 + t155 * t22 - t156 * t105 + t153 * t20 - t182 * mrSges(5,3) + ((t49 + (-mrSges(3,1) + mrSges(4,2)) * qJD(2)) * t90 + (-mrSges(3,2) + t188) * t121) * t86 + m(5) * t148 + m(6) * (t105 * t9 + t187 * t5 + t8 * t22 - t30 * t4) + m(4) * (-pkin(2) * t131 + t148) + m(7) * (-t10 * t102 - t105 * t7 + t20 * t33 + t22 * t6 + t24 * t5 + t25 * t4) + (t95 / 0.2e1 + t185 * t176) * t177; 0.2e1 * qJ(3) * t49 + 0.2e1 * t10 * t43 + 0.2e1 * t33 * t17 + 0.2e1 * t24 * t29 + 0.2e1 * t25 * t27 + 0.2e1 * t30 * t26 + 0.2e1 * t187 * t28 + 0.2e1 * t8 * t54 + 0.2e1 * t9 * t55 + 0.2e1 * t7 * t56 + 0.2e1 * t6 * t57 + 0.2e1 * m(7) * (t10 * t33 + t24 * t6 + t25 * t7) + 0.2e1 * m(6) * (t187 * t8 + t30 * t9) + 0.2e1 * ((m(4) + m(5)) * qJ(3) + t188) * qJD(3) + ((0.2e1 * Ifges(5,4) * t89 + t44 * t177 + (-t154 - t167) * t88 + t106 * t91) * qJD(4) + t120) * t89 + (-0.2e1 * t93 * t18 + (t14 + t15) * t91 + (t12 - t13) * t88 + ((t154 - t163) * t91 + t106 * t88) * qJD(5) + ((-0.2e1 * Ifges(5,4) + t186 + (-Ifges(6,6) + Ifges(7,6)) * t88) * t92 + (t93 ^ 2 * t191 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + t193) * t89) * qJD(4)) * t92; m(4) * t131 + t95 + t88 * (-t144 * t105 + (-t143 + t4) * t89) * t194 + t190 * (t91 * t22 * t144 + t185 + (-t126 + t173) * t89); (-m(7) * t10 + (t151 * t91 + t150 * t88 + m(7) * (t24 * t91 + t25 * t88) + m(6) * (t187 * t91 - t30 * t88)) * qJD(4) - t157) * t92 + (t129 * t191 + (m(7) * t33 + t153) * qJD(4) + t179 * t91 + t178 * t88) * t89; 0.4e1 * t184 * (-0.1e1 + t146) * t89 * t144; -t19 * mrSges(5,2) - t152 * t102 + t134 * t20 + m(6) * (-pkin(4) * t20 + t3) + m(7) * (-t102 * t38 + t20 * t59 + t3) + pkin(9) * (-t142 * t22 - t126 + t174) * t194 + t189 * (t174 + t173 + (-t105 * t91 - t22 * t88) * qJD(5)); t59 * t17 + t10 * t64 + t38 * t43 + t33 * t47 - pkin(4) * t18 + m(7) * (t10 * t59 + t33 * t38) + (t82 / 0.2e1 + t83 / 0.2e1 + t80 / 0.2e1 + (-Ifges(5,5) + (-m(6) * pkin(4) + t149) * t93) * qJD(4)) * t89 + (t7 * mrSges(7,2) - t9 * mrSges(6,3) + t14 / 0.2e1 + t15 / 0.2e1 - t125 * t145 + (t34 / 0.2e1 - t35 / 0.2e1 - t163 / 0.2e1 - t187 * mrSges(6,3) - t24 * mrSges(7,2)) * qJD(5) + t178 * pkin(9)) * t88 + (t6 * mrSges(7,2) + t8 * mrSges(6,3) - t12 / 0.2e1 + t13 / 0.2e1 + t124 * t145 + (t36 / 0.2e1 + t37 / 0.2e1 + t25 * mrSges(7,2) - t30 * mrSges(6,3)) * qJD(5) + t179 * pkin(9)) * t91 + (-t93 * t48 + (t52 / 0.2e1 + t53 / 0.2e1) * t91 + (t50 / 0.2e1 - t51 / 0.2e1) * t88 + (-t93 * mrSges(5,2) - Ifges(5,6) + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t91 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t88) * qJD(4) + (t124 * t88 + t125 * t91) * qJD(5)) * t92; t134 * t145 + m(7) * (t145 * t59 + t147) + m(6) * (-pkin(4) * t145 + t147) + (-m(7) * t38 + (t189 * t146 - mrSges(5,2)) * qJD(4) - t152) * t92; -0.2e1 * pkin(4) * t48 + 0.2e1 * t47 * t59 + (-t50 + t51) * t91 + (t52 + t53) * t88 + 0.2e1 * (m(7) * t59 + t64) * t38 + ((t69 + t70) * t91 + (t67 - t68) * t88) * qJD(5); m(7) * qJD(6) * t22 + (-mrSges(6,2) + t183) * t5 + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t4; -pkin(5) * t27 + m(7) * (-pkin(5) * t7 + qJ(6) * t6 + qJD(6) * t24) + qJD(6) * t57 + qJ(6) * t29 + t6 * mrSges(7,3) - t7 * mrSges(7,1) - t8 * mrSges(6,2) + t9 * mrSges(6,1) + (-Ifges(6,6) * t91 - t158 * t88) * t140 + (-Ifges(7,6) * t88 - t186) * t145 + t120; (-m(7) * t109 - t115 - t116) * t144 + t94 * t89; -t180 * mrSges(7,2) - Ifges(6,6) * t142 + t94 * pkin(9) + t80 + t82 + t83; 0.2e1 * t183 * qJD(6); m(7) * t4; m(7) * t7 + t27; m(7) * (t141 * t89 + t144 * t88); (m(7) * pkin(9) + mrSges(7,2)) * t141; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
