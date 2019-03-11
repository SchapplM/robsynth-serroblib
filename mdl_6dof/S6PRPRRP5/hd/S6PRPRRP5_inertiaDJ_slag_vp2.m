% Calculate time derivative of joint inertia matrix for
% S6PRPRRP5
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:05
% EndTime: 2019-03-08 20:14:10
% DurationCPUTime: 2.31s
% Computational Cost: add. (1352->328), mult. (3486->478), div. (0->0), fcn. (2682->8), ass. (0->159)
t189 = Ifges(6,6) + Ifges(7,6);
t188 = Ifges(6,3) + Ifges(7,3);
t90 = cos(qJ(5));
t132 = qJD(5) * t90;
t87 = sin(qJ(5));
t134 = qJD(5) * t87;
t46 = mrSges(7,1) * t134 + mrSges(7,2) * t132;
t106 = mrSges(6,1) * t87 + mrSges(6,2) * t90;
t47 = t106 * qJD(5);
t187 = t46 + t47;
t62 = -mrSges(7,1) * t90 + mrSges(7,2) * t87;
t76 = -pkin(5) * t90 - pkin(4);
t186 = m(7) * t76 + t62;
t88 = sin(qJ(4));
t138 = qJD(4) * t88;
t116 = t87 * t138;
t91 = cos(qJ(4));
t131 = qJD(5) * t91;
t117 = t90 * t131;
t97 = t116 - t117;
t173 = m(6) / 0.2e1;
t127 = t173 + m(7) / 0.2e1;
t185 = 0.4e1 * t127;
t93 = -pkin(2) - pkin(8);
t135 = qJD(4) * t93;
t120 = t91 * t135;
t159 = t88 * t93;
t167 = pkin(9) * t91;
t168 = pkin(4) * t88;
t57 = qJ(3) - t167 + t168;
t28 = t90 * t159 + t87 * t57;
t44 = qJD(3) + (pkin(4) * t91 + pkin(9) * t88) * qJD(4);
t96 = -t28 * qJD(5) + t90 * t44;
t8 = -t120 * t87 + t96;
t184 = m(6) * t8;
t111 = -t87 * t93 + pkin(5);
t137 = qJD(4) * t90;
t121 = t88 * t137;
t130 = qJD(6) * t90;
t3 = qJ(6) * t121 + (qJ(6) * t134 + qJD(4) * t111 - t130) * t91 + t96;
t183 = m(7) * t3;
t181 = t88 * mrSges(5,1) + t91 * mrSges(5,2) + mrSges(4,3);
t142 = t87 ^ 2 + t90 ^ 2;
t180 = -m(6) * pkin(9) - mrSges(6,3);
t43 = t90 * t57;
t27 = -t159 * t87 + t43;
t133 = qJD(5) * t88;
t119 = t87 * t133;
t125 = t90 * t120 + t57 * t132 + t87 * t44;
t7 = -t119 * t93 + t125;
t179 = m(6) * (-t27 * qJD(5) + t7);
t144 = -mrSges(6,1) * t90 + mrSges(6,2) * t87 - mrSges(5,1);
t178 = t144 + t186;
t85 = sin(pkin(6));
t89 = sin(qJ(2));
t161 = t85 * t89;
t123 = qJD(2) * t161;
t92 = cos(qJ(2));
t160 = t85 * t92;
t128 = t88 * t160;
t136 = qJD(4) * t91;
t86 = cos(pkin(6));
t18 = -qJD(4) * t128 - t123 * t91 + t136 * t86;
t156 = t91 * t18;
t37 = t160 * t91 + t86 * t88;
t139 = qJD(4) * t37;
t17 = t123 * t88 - t139;
t38 = t86 * t91 - t128;
t177 = qJD(4) * (t37 * t88 + t38 * t91) + t88 * t17 - t156;
t176 = 0.2e1 * m(7);
t175 = -2 * mrSges(7,3);
t174 = 0.2e1 * t93;
t172 = m(6) * pkin(4);
t170 = m(7) * pkin(5);
t22 = -pkin(5) * t97 + t135 * t88;
t169 = m(7) * t22;
t166 = Ifges(6,4) * t87;
t165 = Ifges(6,4) * t90;
t164 = Ifges(7,4) * t87;
t163 = Ifges(7,4) * t90;
t162 = t18 * t37;
t158 = t91 * mrSges(6,3);
t157 = t91 * mrSges(7,3);
t155 = mrSges(6,2) + mrSges(7,2);
t154 = Ifges(6,5) + Ifges(7,5);
t152 = -qJ(6) - pkin(9);
t118 = t87 * t131;
t98 = t118 + t121;
t13 = -t97 * mrSges(7,1) - mrSges(7,2) * t98;
t14 = -mrSges(6,1) * t97 - mrSges(6,2) * t98;
t151 = t13 + t14;
t23 = mrSges(7,1) * t136 + mrSges(7,3) * t98;
t24 = mrSges(6,1) * t136 + mrSges(6,3) * t98;
t150 = t23 + t24;
t25 = -mrSges(7,2) * t136 + mrSges(7,3) * t97;
t26 = -mrSges(6,2) * t136 + mrSges(6,3) * t97;
t149 = t25 + t26;
t102 = -Ifges(7,2) * t87 + t163;
t29 = t88 * Ifges(7,6) + t102 * t91;
t103 = -Ifges(6,2) * t87 + t165;
t30 = t88 * Ifges(6,6) + t103 * t91;
t148 = t29 + t30;
t40 = (mrSges(7,1) * t87 + mrSges(7,2) * t90) * t91;
t41 = t106 * t91;
t147 = t40 + t41;
t53 = -mrSges(7,2) * t88 - t157 * t87;
t54 = -mrSges(6,2) * t88 - t158 * t87;
t146 = t53 + t54;
t55 = mrSges(7,1) * t88 - t157 * t90;
t56 = mrSges(6,1) * t88 - t158 * t90;
t145 = t56 + t55;
t140 = qJD(2) * t92;
t122 = t85 * t140;
t143 = qJ(3) * t122 + qJD(3) * t161;
t141 = qJ(6) * t91;
t124 = pkin(5) * t134;
t115 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t66 = Ifges(7,2) * t90 + t164;
t67 = Ifges(6,2) * t90 + t166;
t114 = t66 / 0.2e1 + t67 / 0.2e1;
t68 = Ifges(7,1) * t87 + t163;
t69 = Ifges(6,1) * t87 + t165;
t113 = -t68 / 0.2e1 - t69 / 0.2e1;
t112 = mrSges(7,1) + t170;
t110 = t154 * t88;
t109 = qJD(5) * t152;
t108 = t189 * t116 + t188 * t136;
t107 = mrSges(6,1) + t112;
t105 = Ifges(6,1) * t90 - t166;
t104 = Ifges(7,1) * t90 - t164;
t31 = Ifges(7,5) * t88 + t104 * t91;
t32 = Ifges(6,5) * t88 + t105 * t91;
t99 = -t31 - t32 - t110;
t19 = t161 * t90 - t38 * t87;
t20 = t161 * t87 + t38 * t90;
t5 = -t20 * qJD(5) + t122 * t90 - t17 * t87;
t6 = t19 * qJD(5) + t122 * t87 + t17 * t90;
t95 = -t5 * t87 + t6 * t90 + (-t19 * t90 - t20 * t87) * qJD(5);
t94 = m(5) * t177;
t82 = Ifges(6,5) * t132;
t81 = Ifges(7,5) * t132;
t64 = t152 * t90;
t61 = t152 * t87;
t52 = t105 * qJD(5);
t51 = t104 * qJD(5);
t50 = t103 * qJD(5);
t49 = t102 * qJD(5);
t48 = (mrSges(5,1) * t91 - mrSges(5,2) * t88) * qJD(4);
t45 = (pkin(5) * t87 - t93) * t91;
t34 = -qJD(6) * t87 + t109 * t90;
t33 = t109 * t87 + t130;
t16 = -t141 * t87 + t28;
t15 = t111 * t88 - t141 * t90 + t43;
t12 = -t69 * t131 + (Ifges(6,5) * t91 - t105 * t88) * qJD(4);
t11 = -t68 * t131 + (Ifges(7,5) * t91 - t104 * t88) * qJD(4);
t10 = -t67 * t131 + (Ifges(6,6) * t91 - t103 * t88) * qJD(4);
t9 = -t66 * t131 + (Ifges(7,6) * t91 - t102 * t88) * qJD(4);
t4 = -qJ(6) * t117 + (-qJD(6) * t91 + (qJ(6) * qJD(4) - qJD(5) * t93) * t88) * t87 + t125;
t1 = [0.2e1 * m(5) * (t85 ^ 2 * t89 * t140 + t17 * t38 + t162) + (t19 * t5 + t20 * t6 + t162) * t185; t146 * t6 + t145 * t5 + t151 * t37 + t149 * t20 + t150 * t19 + t147 * t18 - t177 * mrSges(5,3) + (t89 * t48 + ((-mrSges(3,1) + mrSges(4,2)) * t89 + (-mrSges(3,2) + t181) * t92) * qJD(2)) * t85 + m(5) * t143 + m(6) * (t8 * t19 + t7 * t20 + t27 * t5 + t28 * t6) + m(4) * (-pkin(2) * t123 + t143) + m(7) * (t15 * t5 + t16 * t6 + t18 * t45 + t19 * t3 + t20 * t4 + t22 * t37) + (t94 / 0.2e1 + (t138 * t37 - t156) * t173) * t174; 0.2e1 * qJ(3) * t48 + 0.2e1 * t45 * t13 + 0.2e1 * t15 * t23 + 0.2e1 * t16 * t25 + 0.2e1 * t22 * t40 + 0.2e1 * t27 * t24 + 0.2e1 * t28 * t26 + 0.2e1 * t3 * t55 + 0.2e1 * t4 * t53 + 0.2e1 * t7 * t54 + 0.2e1 * t8 * t56 + (t15 * t3 + t16 * t4 + t22 * t45) * t176 + 0.2e1 * m(6) * (t27 * t8 + t28 * t7) + 0.2e1 * ((m(4) + m(5)) * qJ(3) + t181) * qJD(3) + ((0.2e1 * Ifges(5,4) * t88 + t148 * t87 + t174 * t41 + t90 * t99) * qJD(4) + t108) * t88 + (-0.2e1 * t93 * t14 + (t11 + t12) * t90 + (-t10 - t9) * t87 + ((-t189 * t88 - t148) * t90 + t99 * t87) * qJD(5) + ((t154 * t90 - t189 * t87 - 0.2e1 * Ifges(5,4)) * t91 + (-0.2e1 * m(6) * t93 ^ 2 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + t188) * t88) * qJD(4)) * t91; m(4) * t123 + t94 + 0.2e1 * t127 * ((-t18 + (-t19 * t87 + t20 * t90) * qJD(4)) * t91 + (t95 + t139) * t88); (-t169 + (t146 * t90 - t145 * t87 + m(7) * (-t15 * t87 + t16 * t90) + m(6) * (-t27 * t87 + t28 * t90)) * qJD(4) - t151) * t91 + (t147 * qJD(4) + m(7) * (qJD(4) * t45 - t132 * t15 - t134 * t16) + m(6) * (-t134 * t28 - 0.2e1 * t120) + (-t146 * qJD(5) - t150 - t183 - t184) * t87 + (m(7) * t4 - t145 * qJD(5) + t149 + t179) * t90) * t88; (-0.1e1 + t142) * t88 * t136 * t185; -t17 * mrSges(5,2) + t187 * t37 + m(7) * (t124 * t37 + t19 * t34 + t20 * t33 + t5 * t61 - t6 * t64) + (-t172 + t178) * t18 + (mrSges(7,3) - t180) * t95; m(7) * (t15 * t34 + t16 * t33 + t22 * t76 + t3 * t61 - t4 * t64) + t76 * t13 + t61 * t23 + t22 * t62 - t64 * t25 + t45 * t46 + t33 * t53 + t34 * t55 - pkin(4) * t14 + (t81 / 0.2e1 + t82 / 0.2e1 + (-Ifges(5,5) + (t144 - t172) * t93) * qJD(4)) * t88 + (-t8 * mrSges(6,3) - t3 * mrSges(7,3) + t11 / 0.2e1 + t12 / 0.2e1 + t114 * t138 + (-t24 - t184) * pkin(9) + (pkin(5) * t40 - pkin(9) * t54 - t16 * mrSges(7,3) - t29 / 0.2e1 - t30 / 0.2e1 - t115 * t88 + t45 * t170 + t180 * t28) * qJD(5)) * t87 + (t7 * mrSges(6,3) + t4 * mrSges(7,3) + t9 / 0.2e1 + t10 / 0.2e1 + t113 * t138 + (-t15 * mrSges(7,3) - t27 * mrSges(6,3) + t31 / 0.2e1 + t32 / 0.2e1) * qJD(5) + (-qJD(5) * t56 + t179 + t26) * pkin(9)) * t90 + (-t93 * t47 + (t51 / 0.2e1 + t52 / 0.2e1) * t90 + (-t49 / 0.2e1 - t50 / 0.2e1) * t87 + (-t93 * mrSges(5,2) - Ifges(5,6) + t115 * t90 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t87) * qJD(4) + (t113 * t87 - t114 * t90) * qJD(5)) * t91; m(7) * (-pkin(5) * t118 + t119 * t64 + (-t132 * t61 + t33 * t90 - t34 * t87) * t88) + (m(6) * (t142 * t167 - t168) + t178 * t88) * qJD(4) + ((-mrSges(5,2) + m(7) * (-t61 * t87 - t64 * t90) + (mrSges(6,3) + mrSges(7,3)) * t142) * qJD(4) - t187) * t91; 0.2e1 * t76 * t46 - 0.2e1 * pkin(4) * t47 + (-t33 * t64 + t34 * t61) * t176 + (t34 * t175 + t51 + t52 + (0.2e1 * t186 * pkin(5) - t64 * t175 - t66 - t67) * qJD(5)) * t87 + (0.2e1 * t33 * mrSges(7,3) + t49 + t50 + (t175 * t61 + t68 + t69) * qJD(5)) * t90; t107 * t5 - t155 * t6; mrSges(6,1) * t8 + mrSges(7,1) * t3 - mrSges(6,2) * t7 - mrSges(7,2) * t4 - t110 * t137 + (t23 + t183) * pkin(5) + (-t154 * t87 - t189 * t90) * t131 + t108; (-t107 * t90 + t155 * t87) * t133 + (-t107 * t87 - t155 * t90) * t136; -mrSges(7,2) * t33 + t81 + t82 + t112 * t34 + ((-mrSges(6,1) * pkin(9) - mrSges(7,3) * pkin(5)) * t90 + (mrSges(6,2) * pkin(9) - t189) * t87) * qJD(5); 0; m(7) * t18; t13 + t169; m(7) * t138; m(7) * t124 + t46; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
