% Calculate time derivative of joint inertia matrix for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:31
% EndTime: 2019-12-31 21:07:38
% DurationCPUTime: 2.44s
% Computational Cost: add. (1037->340), mult. (2647->465), div. (0->0), fcn. (1677->4), ass. (0->149)
t184 = -Ifges(6,4) - Ifges(5,5);
t170 = Ifges(4,5) + Ifges(6,5);
t187 = -Ifges(5,4) + t170;
t107 = sin(qJ(3));
t108 = sin(qJ(2));
t109 = cos(qJ(3));
t146 = qJD(3) * t109;
t110 = cos(qJ(2));
t149 = qJD(2) * t110;
t112 = t107 * t149 + t108 * t146;
t186 = -Ifges(5,1) - Ifges(6,1) - Ifges(4,3);
t148 = qJD(3) * t107;
t150 = qJD(2) * t108;
t185 = pkin(6) * (t109 * t150 + t110 * t148);
t163 = Ifges(5,6) * t109;
t122 = Ifges(5,3) * t107 - t163;
t161 = Ifges(6,6) * t109;
t123 = Ifges(6,2) * t107 + t161;
t183 = (t122 + t123) * qJD(3);
t162 = Ifges(6,6) * t107;
t119 = Ifges(6,3) * t109 + t162;
t167 = Ifges(4,4) * t107;
t128 = Ifges(4,1) * t109 - t167;
t182 = (t119 + t128) * qJD(3);
t181 = 2 * m(4);
t180 = 2 * m(6);
t179 = -2 * pkin(1);
t178 = 0.2e1 * pkin(6);
t177 = 2 * mrSges(6,1);
t176 = pkin(4) + pkin(7);
t173 = pkin(6) * t107;
t171 = -mrSges(5,1) - mrSges(6,1);
t106 = -pkin(3) - qJ(5);
t64 = (pkin(2) * t108 - pkin(7) * t110) * qJD(2);
t66 = -pkin(2) * t110 - t108 * pkin(7) - pkin(1);
t169 = t107 * t64 + t66 * t146;
t138 = t109 * t149;
t168 = mrSges(5,1) * t138 + mrSges(5,2) * t150;
t152 = t109 * t110;
t92 = pkin(6) * t152;
t31 = t107 * t66 + t92;
t166 = Ifges(4,4) * t109;
t165 = Ifges(4,6) * t109;
t164 = Ifges(5,6) * t107;
t160 = t110 * Ifges(5,4);
t159 = t110 * Ifges(4,6);
t155 = t107 * t108;
t158 = pkin(3) * t155 + t108 * pkin(6);
t157 = qJ(4) * t107;
t156 = qJD(2) * mrSges(6,3);
t154 = t107 * t110;
t153 = t108 * t109;
t151 = qJ(4) * qJD(4);
t147 = qJD(3) * t108;
t145 = qJD(5) * t107;
t33 = -Ifges(4,5) * t110 + t108 * t128;
t34 = -Ifges(6,5) * t110 + t108 * t119;
t126 = -Ifges(5,2) * t109 + t164;
t37 = t108 * t126 - t160;
t144 = t37 - t33 - t34;
t91 = pkin(6) * t154;
t76 = t176 * t109;
t143 = mrSges(6,1) * t148;
t142 = -pkin(3) - t173;
t141 = t107 * t147;
t30 = t109 * t66 - t91;
t137 = pkin(3) * t112 + pkin(6) * t149 + qJ(4) * t141;
t136 = pkin(3) * t148 - qJD(4) * t107;
t135 = qJD(3) * t92 - t109 * t64 + t66 * t148;
t120 = -Ifges(6,3) * t107 + t161;
t125 = Ifges(5,2) * t107 + t163;
t74 = Ifges(4,1) * t107 + t166;
t134 = t120 / 0.2e1 - t125 / 0.2e1 - t74 / 0.2e1;
t121 = Ifges(5,3) * t109 + t164;
t124 = Ifges(6,2) * t109 - t162;
t73 = Ifges(4,2) * t109 + t167;
t133 = -t124 / 0.2e1 - t73 / 0.2e1 - t121 / 0.2e1;
t22 = qJ(4) * t110 - t31;
t132 = -qJD(4) * t110 + t169;
t131 = mrSges(4,1) * t107 + mrSges(4,2) * t109;
t67 = t109 * mrSges(5,2) - t107 * mrSges(5,3);
t130 = -mrSges(5,2) * t107 - mrSges(5,3) * t109;
t129 = -mrSges(6,2) * t109 + mrSges(6,3) * t107;
t127 = -Ifges(4,2) * t107 + t166;
t118 = -pkin(3) * t109 - t157;
t32 = t108 * t127 - t159;
t35 = -Ifges(5,5) * t110 + t108 * t122;
t36 = -Ifges(6,4) * t110 + t108 * t123;
t117 = -t32 + t35 + t36 + t159;
t116 = -qJ(4) * t109 + qJ(5) * t107;
t115 = (m(5) * pkin(7) - t171) * t109;
t114 = t138 - t141;
t111 = -Ifges(5,4) * t141 + t184 * t112 - t138 * t170 + t186 * t150;
t28 = -mrSges(6,1) * t112 + mrSges(6,2) * t150;
t105 = t110 * pkin(3);
t99 = Ifges(6,4) * t148;
t98 = Ifges(4,5) * t146;
t97 = Ifges(5,5) * t148;
t96 = Ifges(6,5) * t146;
t78 = mrSges(6,1) * t138;
t75 = t176 * t107;
t68 = -mrSges(6,2) * t107 - mrSges(6,3) * t109;
t65 = -pkin(2) + t118;
t63 = qJD(3) * t76;
t62 = t176 * t148;
t61 = mrSges(5,1) * t153 - mrSges(5,2) * t110;
t60 = -mrSges(6,1) * t155 - mrSges(6,2) * t110;
t59 = mrSges(5,1) * t155 + mrSges(5,3) * t110;
t58 = mrSges(6,1) * t153 + mrSges(6,3) * t110;
t57 = -mrSges(4,1) * t110 - mrSges(4,3) * t153;
t56 = mrSges(4,2) * t110 - mrSges(4,3) * t155;
t54 = t127 * qJD(3);
t53 = t126 * qJD(3);
t49 = t131 * qJD(3);
t48 = t130 * qJD(3);
t47 = t129 * qJD(3);
t43 = t130 * t108;
t42 = t129 * t108;
t41 = t106 * t109 - pkin(2) - t157;
t39 = -qJ(4) * t146 + t136;
t38 = -qJ(4) * t153 + t158;
t29 = -mrSges(5,1) * t141 + t168;
t27 = mrSges(5,1) * t112 - mrSges(5,3) * t150;
t26 = t78 + (-t143 - t156) * t108;
t25 = -mrSges(4,2) * t150 - mrSges(4,3) * t112;
t24 = mrSges(4,1) * t150 - mrSges(4,3) * t114;
t23 = t105 - t30;
t21 = t108 * t116 + t158;
t20 = qJD(3) * t116 - qJD(5) * t109 + t136;
t19 = -pkin(4) * t155 - t22;
t18 = -mrSges(6,2) * t114 + mrSges(6,3) * t112;
t17 = mrSges(4,1) * t112 + mrSges(4,2) * t114;
t16 = -mrSges(5,2) * t112 - mrSges(5,3) * t114;
t15 = qJ(5) * t110 + t105 + t91 + (pkin(4) * t108 - t66) * t109;
t14 = t125 * t147 + (Ifges(5,4) * t108 + t110 * t126) * qJD(2);
t13 = t124 * t147 + (Ifges(6,4) * t108 + t110 * t123) * qJD(2);
t12 = t121 * t147 + (Ifges(5,5) * t108 + t110 * t122) * qJD(2);
t11 = t120 * t147 + (Ifges(6,5) * t108 + t110 * t119) * qJD(2);
t10 = -t74 * t147 + (Ifges(4,5) * t108 + t110 * t128) * qJD(2);
t9 = -t73 * t147 + (Ifges(4,6) * t108 + t110 * t127) * qJD(2);
t8 = t150 * t173 - t135;
t7 = t169 - t185;
t6 = (-qJ(4) * t149 - qJD(4) * t108) * t109 + t137;
t5 = t142 * t150 + t135;
t4 = -qJ(4) * t150 - t132 + t185;
t3 = t116 * t149 + (t145 + (qJ(5) * qJD(3) - qJD(4)) * t109) * t108 + t137;
t2 = (-pkin(4) * t153 - t91) * qJD(3) + (-pkin(4) * t154 + (-pkin(6) * t109 + qJ(4)) * t108) * qJD(2) + t132;
t1 = -pkin(4) * t141 + qJD(5) * t110 + (pkin(4) * t152 + (-qJ(5) + t142) * t108) * qJD(2) + t135;
t40 = [0.2e1 * t7 * t56 + 0.2e1 * t8 * t57 + 0.2e1 * t1 * t58 + 0.2e1 * t4 * t59 + 0.2e1 * t2 * t60 + 0.2e1 * t5 * t61 + 0.2e1 * t31 * t25 + 0.2e1 * t38 * t16 + 0.2e1 * t3 * t42 + 0.2e1 * t6 * t43 + 0.2e1 * t21 * t18 + 0.2e1 * t15 * t26 + 0.2e1 * t22 * t27 + 0.2e1 * t19 * t28 + 0.2e1 * t23 * t29 + 0.2e1 * t30 * t24 + (t30 * t8 + t31 * t7) * t181 + 0.2e1 * m(5) * (t22 * t4 + t23 * t5 + t38 * t6) + (t1 * t15 + t19 * t2 + t21 * t3) * t180 + (((mrSges(3,2) * t179) + 0.2e1 * Ifges(3,4) * t110 + (-t144 + t160) * t109 + t117 * t107) * qJD(2) + t111) * t110 + (t17 * t178 + (t10 + t11 - t14) * t109 + (t12 + t13 - t9) * t107 + (t117 * t109 + (t110 * t170 + t144) * t107) * qJD(3) + ((mrSges(3,1) * t179) - 0.2e1 * Ifges(3,4) * t108 + t187 * t153 + (-Ifges(4,6) - t184) * t155 + (pkin(6) ^ 2 * t181 + t131 * t178 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + t186) * t110) * qJD(2)) * t108; -pkin(2) * t17 + t65 * t16 + t41 * t18 + t20 * t42 + t21 * t47 + t75 * t26 + t76 * t28 + t3 * t68 + t38 * t48 + t39 * t43 + t63 * t58 + t6 * t67 - t62 * t60 + m(6) * (t1 * t75 + t15 * t63 - t19 * t62 + t2 * t76 + t20 * t21 + t3 * t41) + m(5) * (t38 * t39 + t6 * t65) + (-t97 / 0.2e1 - t99 / 0.2e1 - t96 / 0.2e1 - t98 / 0.2e1 + (Ifges(3,5) + (-m(4) * pkin(2) - mrSges(3,1)) * pkin(6)) * qJD(2)) * t110 + (t7 * mrSges(4,3) + t2 * mrSges(6,1) - t4 * mrSges(5,1) + t9 / 0.2e1 - t12 / 0.2e1 - t13 / 0.2e1 + (-pkin(6) * mrSges(4,1) - t134) * t149 + (t33 / 0.2e1 + t34 / 0.2e1 - t37 / 0.2e1 + t160 / 0.2e1 + t23 * mrSges(5,1) - t30 * mrSges(4,3) + t15 * mrSges(6,1)) * qJD(3) + (t25 - t27 + (-t57 + t61) * qJD(3) + m(5) * (t23 * qJD(3) - t4) + m(4) * (-t30 * qJD(3) + t7)) * pkin(7)) * t109 + (-t8 * mrSges(4,3) + t1 * mrSges(6,1) + t5 * mrSges(5,1) + t10 / 0.2e1 + t11 / 0.2e1 - t14 / 0.2e1 + (pkin(6) * mrSges(4,2) + t133) * t149 + (-t32 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1 + t159 / 0.2e1 + t22 * mrSges(5,1) - t31 * mrSges(4,3) - t19 * mrSges(6,1)) * qJD(3) + (-t24 + t29 + (-t56 + t59) * qJD(3) + m(5) * (t22 * qJD(3) + t5) + m(4) * (-t31 * qJD(3) - t8)) * pkin(7)) * t107 + (-t109 * t53 / 0.2e1 - Ifges(3,6) * qJD(2) - t107 * t54 / 0.2e1 + (qJD(2) * mrSges(3,2) + t49) * pkin(6) + (t107 * t134 + t109 * t133) * qJD(3) + t183 * t107 / 0.2e1 + t182 * t109 / 0.2e1 + (t187 * t107 + t184 * t109 + t165) * qJD(2) / 0.2e1) * t108; 0.2e1 * t65 * t48 + 0.2e1 * t20 * t68 + 0.2e1 * t41 * t47 + (t20 * t41 - t62 * t76 + t63 * t75) * t180 - 0.2e1 * pkin(2) * t49 + 0.2e1 * (m(5) * t65 + t67) * t39 + (-t62 * t177 - t183 + t54) * t109 + (t63 * t177 + t182 - t53) * t107 + ((t177 * t75 - t120 + t125 + t74) * t109 + (-0.2e1 * mrSges(6,1) * t76 - t121 - t124 - t73) * t107) * qJD(3); (-t59 + t60) * qJD(4) + (-t27 + t28) * qJ(4) + m(6) * (qJ(4) * t2 + qJD(4) * t19 - qJD(5) * t15 + t1 * t106) + m(5) * (-pkin(3) * t5 - qJ(4) * t4 - qJD(4) * t22) + (-Ifges(5,4) * t109 - Ifges(4,6) * t107) * t149 + (-t107 * t170 - t165) * t147 + t106 * t26 - qJD(5) * t58 - pkin(3) * t29 - t7 * mrSges(4,2) + t8 * mrSges(4,1) - t1 * mrSges(6,3) + t2 * mrSges(6,2) - t4 * mrSges(5,3) + t5 * mrSges(5,2) - t111; -mrSges(6,1) * t145 + m(6) * (-qJ(4) * t62 + qJD(4) * t76 - qJD(5) * t75 + t106 * t63) + t99 + t96 + t97 + t98 - t62 * mrSges(6,2) - t63 * mrSges(6,3) + qJD(4) * t115 + ((-pkin(3) * mrSges(5,1) + t106 * mrSges(6,1) - Ifges(5,4)) * t109 + (qJ(4) * t171 - Ifges(4,6)) * t107 + (m(5) * t118 - t109 * mrSges(4,1) + t107 * mrSges(4,2) + t67) * pkin(7)) * qJD(3); 0.2e1 * m(5) * t151 + 0.2e1 * qJD(5) * mrSges(6,3) + 0.2e1 * m(6) * (-qJD(5) * t106 + t151) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJD(4); t78 + m(5) * t5 + m(6) * t1 + (t148 * t171 - t156) * t108 + t168; m(6) * t63 + qJD(3) * t115; -m(6) * qJD(5); 0; m(6) * t2 + t28; -m(6) * t62 - t143; m(6) * qJD(4); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t40(1), t40(2), t40(4), t40(7), t40(11); t40(2), t40(3), t40(5), t40(8), t40(12); t40(4), t40(5), t40(6), t40(9), t40(13); t40(7), t40(8), t40(9), t40(10), t40(14); t40(11), t40(12), t40(13), t40(14), t40(15);];
Mq = res;
