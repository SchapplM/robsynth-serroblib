% Calculate time derivative of joint inertia matrix for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:06:06
% EndTime: 2019-03-09 08:06:14
% DurationCPUTime: 3.61s
% Computational Cost: add. (3051->326), mult. (6916->499), div. (0->0), fcn. (6410->8), ass. (0->140)
t187 = Ifges(5,1) + Ifges(6,1);
t182 = Ifges(6,4) + Ifges(5,5);
t181 = Ifges(5,6) - Ifges(6,6);
t116 = cos(pkin(10));
t114 = sin(pkin(10));
t157 = Ifges(6,5) * t114;
t159 = Ifges(5,4) * t114;
t186 = t187 * t116 + t157 - t159;
t185 = qJD(4) * (t114 ^ 2 + t116 ^ 2);
t115 = sin(pkin(9));
t118 = sin(qJ(2));
t151 = cos(pkin(9));
t167 = cos(qJ(2));
t100 = t115 * t167 + t118 * t151;
t184 = 0.2e1 * t100;
t183 = Ifges(3,1) - Ifges(3,2);
t180 = -qJ(3) - pkin(7);
t89 = t100 * qJD(2);
t121 = -t115 * t118 + t151 * t167;
t90 = t121 * qJD(2);
t179 = t182 * t89 + t186 * t90;
t154 = t114 * t90;
t52 = -mrSges(5,2) * t89 - mrSges(5,3) * t154;
t55 = -mrSges(6,2) * t154 + mrSges(6,3) * t89;
t178 = t52 + t55;
t152 = t116 * t90;
t53 = mrSges(5,1) * t89 - mrSges(5,3) * t152;
t54 = -t89 * mrSges(6,1) + mrSges(6,2) * t152;
t177 = t54 - t53;
t144 = qJD(2) * t118;
t137 = pkin(2) * t144;
t38 = pkin(3) * t89 - qJ(4) * t90 - qJD(4) * t100 + t137;
t107 = t180 * t167;
t120 = qJD(2) * t107 - t118 * qJD(3);
t106 = t180 * t118;
t88 = qJD(2) * t106 + qJD(3) * t167;
t48 = t115 * t120 + t151 * t88;
t43 = t114 * t48;
t16 = t116 * t38 - t43;
t17 = t114 * t38 + t116 * t48;
t176 = -t114 * t16 + t116 * t17;
t175 = -t181 * t114 + t182 * t116 - (2 * Ifges(4,4));
t174 = 2 * m(5);
t173 = 2 * m(6);
t172 = 2 * m(7);
t169 = pkin(4) + pkin(5);
t166 = pkin(2) * t115;
t47 = t115 * t88 - t151 * t120;
t72 = -t151 * t106 - t107 * t115;
t165 = t47 * t72;
t164 = t89 * mrSges(4,3);
t163 = t90 * mrSges(4,3);
t108 = qJ(4) + t166;
t162 = -pkin(8) + t108;
t111 = -pkin(2) * t167 - pkin(1);
t62 = -pkin(3) * t121 - t100 * qJ(4) + t111;
t73 = t115 * t106 - t107 * t151;
t31 = t114 * t62 + t116 * t73;
t46 = mrSges(5,1) * t154 + mrSges(5,2) * t152;
t117 = sin(qJ(6));
t119 = cos(qJ(6));
t123 = t117 * t114 + t116 * t119;
t91 = t123 * qJD(6);
t101 = t114 * t119 - t117 * t116;
t92 = t101 * qJD(6);
t161 = -Ifges(7,5) * t91 - Ifges(7,6) * t92;
t158 = Ifges(5,4) * t116;
t156 = Ifges(6,5) * t116;
t150 = qJ(5) * t116;
t149 = t100 * t114;
t148 = t100 * t116;
t145 = t108 * t185;
t142 = qJD(4) * t114;
t140 = qJD(5) * t114;
t139 = t47 * t184;
t49 = t101 * t100;
t24 = qJD(6) * t49 + t123 * t90;
t25 = -t100 * t91 + t101 * t90;
t138 = Ifges(7,5) * t24 + Ifges(7,6) * t25 - Ifges(7,3) * t89;
t26 = -qJ(5) * t121 + t31;
t134 = t151 * pkin(2);
t133 = t89 * mrSges(4,1) + t90 * mrSges(4,2);
t59 = t92 * mrSges(7,1) - mrSges(7,2) * t91;
t132 = qJD(2) * t167;
t70 = t114 * t73;
t30 = t116 * t62 - t70;
t10 = t89 * qJ(5) - qJD(5) * t121 + t17;
t45 = mrSges(6,1) * t154 - mrSges(6,3) * t152;
t110 = -t134 - pkin(3);
t7 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t128 = -t114 * Ifges(5,2) + t158;
t125 = t114 * Ifges(6,3) + t156;
t13 = t70 + t169 * t121 + (-pkin(8) * t100 - t62) * t116;
t18 = pkin(8) * t149 + t26;
t3 = -t117 * t18 + t119 * t13;
t4 = t117 * t13 + t119 * t18;
t94 = t162 * t114;
t95 = t162 * t116;
t56 = -t117 * t95 + t119 * t94;
t57 = t117 * t94 + t119 * t95;
t124 = -qJD(5) * t148 - t90 * t150 + t47;
t122 = -t114 * qJ(5) + t110;
t105 = -mrSges(6,1) * t116 - mrSges(6,3) * t114;
t93 = -t116 * pkin(4) + t122;
t80 = t116 * t169 - t122;
t69 = Ifges(7,1) * t101 - Ifges(7,4) * t123;
t68 = Ifges(7,4) * t101 - Ifges(7,2) * t123;
t67 = mrSges(7,1) * t123 + mrSges(7,2) * t101;
t66 = -mrSges(6,2) * t149 - mrSges(6,3) * t121;
t65 = mrSges(6,1) * t121 + mrSges(6,2) * t148;
t64 = -mrSges(5,1) * t121 - mrSges(5,3) * t148;
t63 = mrSges(5,2) * t121 - mrSges(5,3) * t149;
t61 = -Ifges(7,1) * t91 - Ifges(7,4) * t92;
t60 = -Ifges(7,4) * t91 - Ifges(7,2) * t92;
t58 = (mrSges(6,1) * t114 - mrSges(6,3) * t116) * t100;
t50 = t123 * t100;
t42 = qJD(4) * t101 - qJD(6) * t57;
t41 = qJD(4) * t123 + qJD(6) * t56;
t40 = mrSges(7,1) * t121 - mrSges(7,3) * t50;
t39 = -mrSges(7,2) * t121 + mrSges(7,3) * t49;
t35 = t89 * Ifges(5,6) + t128 * t90;
t34 = t89 * Ifges(6,6) + t125 * t90;
t32 = (pkin(4) * t114 - t150) * t100 + t72;
t29 = (-t114 * t169 + t150) * t100 - t72;
t28 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t27 = pkin(4) * t121 - t30;
t21 = Ifges(7,1) * t50 + Ifges(7,4) * t49 + Ifges(7,5) * t121;
t20 = Ifges(7,4) * t50 + Ifges(7,2) * t49 + Ifges(7,6) * t121;
t19 = pkin(4) * t154 + t124;
t15 = mrSges(7,2) * t89 + mrSges(7,3) * t25;
t14 = -mrSges(7,1) * t89 - mrSges(7,3) * t24;
t12 = t154 * t169 + t124;
t11 = -pkin(4) * t89 - t16;
t9 = pkin(8) * t154 + t10;
t8 = t43 - t169 * t89 + (-pkin(8) * t90 - t38) * t116;
t6 = Ifges(7,1) * t24 + Ifges(7,4) * t25 - t89 * Ifges(7,5);
t5 = Ifges(7,4) * t24 + Ifges(7,2) * t25 - t89 * Ifges(7,6);
t2 = -qJD(6) * t4 - t117 * t9 + t119 * t8;
t1 = qJD(6) * t3 + t117 * t8 + t119 * t9;
t22 = [0.2e1 * (t163 + t46) * t72 + (0.2e1 * Ifges(3,4) * t167 + t183 * t118) * t132 + (-0.2e1 * Ifges(3,4) * t118 + t183 * t167) * t144 + t179 * t148 + (t186 * t100 - t182 * t121) * t152 + (t10 * t26 + t11 * t27 + t19 * t32) * t173 + (t16 * t30 + t17 * t31 + t165) * t174 + (t1 * t4 - t12 * t29 + t2 * t3) * t172 + (mrSges(5,1) * t114 + mrSges(5,2) * t116) * t139 + (-Ifges(7,5) * t50 - Ifges(7,6) * t49 + t175 * t100 - ((2 * Ifges(4,2)) + (2 * Ifges(6,2)) + (2 * Ifges(5,3)) + Ifges(7,3)) * t121) * t89 + 0.2e1 * (-mrSges(4,1) * t121 + mrSges(4,2) * t100) * t137 + (0.2e1 * t121 * t48 + t139) * mrSges(4,3) + (Ifges(4,1) * t184 - t121 * t175) * t90 + (t181 * t121 + (t125 - t128) * t100) * t154 + t121 * t138 + (t34 - t35) * t149 + 0.2e1 * t111 * t133 - 0.2e1 * pkin(1) * (mrSges(3,1) * t118 + mrSges(3,2) * t167) * qJD(2) + 0.2e1 * t3 * t14 + 0.2e1 * t4 * t15 + t24 * t21 + t25 * t20 - 0.2e1 * t12 * t28 + 0.2e1 * t29 * t7 - 0.2e1 * t73 * t164 + 0.2e1 * t1 * t39 + 0.2e1 * m(4) * (t111 * t137 + t48 * t73 + t165) + 0.2e1 * t2 * t40 + 0.2e1 * t32 * t45 + t49 * t5 + t50 * t6 + 0.2e1 * t31 * t52 + 0.2e1 * t30 * t53 + 0.2e1 * t27 * t54 + 0.2e1 * t26 * t55 + 0.2e1 * t19 * t58 + 0.2e1 * t17 * t63 + 0.2e1 * t16 * t64 + 0.2e1 * t11 * t65 + 0.2e1 * t10 * t66; (t187 * t114 - t156 + t158) * t152 / 0.2e1 + (t10 * t116 + t11 * t114) * mrSges(6,2) + t177 * t108 * t114 + t178 * t108 * t116 + t179 * t114 / 0.2e1 + (t182 * t114 + t181 * t116) * t89 / 0.2e1 + t176 * mrSges(5,3) + m(5) * (t110 * t47 + t176 * t108 + (-t114 * t30 + t116 * t31) * qJD(4)) + (t66 + t63) * qJD(4) * t116 + m(6) * (t19 * t93 + (qJD(4) * t26 + t10 * t108) * t116 + (qJD(4) * t27 - qJD(5) * t32 + t108 * t11) * t114) + m(7) * (t1 * t57 - t12 * t80 + t2 * t56 + t3 * t42 + t4 * t41) + Ifges(3,5) * t132 + (-t1 * t123 - t101 * t2 + t3 * t91 - t4 * t92) * mrSges(7,3) - t89 * (Ifges(7,5) * t101 - Ifges(7,6) * t123) / 0.2e1 - t123 * t5 / 0.2e1 + t121 * t161 / 0.2e1 + (-mrSges(3,1) * t132 + mrSges(3,2) * t144) * pkin(7) - Ifges(3,6) * t144 + m(4) * (t115 * t48 - t151 * t47) * pkin(2) + (-Ifges(6,3) * t116 + t157) * t154 / 0.2e1 - (Ifges(5,2) * t116 + t159) * t154 / 0.2e1 - t134 * t163 + t41 * t39 + t42 * t40 - t164 * t166 - t47 * mrSges(4,1) - t48 * mrSges(4,2) + t56 * t14 + t57 * t15 + t29 * t59 + t49 * t60 / 0.2e1 + t50 * t61 / 0.2e1 - t12 * t67 + t25 * t68 / 0.2e1 + t24 * t69 / 0.2e1 + t80 * t7 - Ifges(4,6) * t89 + Ifges(4,5) * t90 - t91 * t21 / 0.2e1 - t92 * t20 / 0.2e1 + t93 * t45 + t101 * t6 / 0.2e1 + t19 * t105 + t110 * t46 - t116 * t34 / 0.2e1 + t116 * t35 / 0.2e1 + t47 * (-mrSges(5,1) * t116 + mrSges(5,2) * t114) + (m(7) * t29 + t28 - t58) * t140 + (t65 - t64) * t142; t101 * t61 + 0.2e1 * t80 * t59 - t123 * t60 - t92 * t68 - t91 * t69 + (t140 * t80 + t41 * t57 + t42 * t56) * t172 + t145 * t174 + (-t140 * t93 + t145) * t173 + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * t185 + 0.2e1 * (-t105 + t67) * t140 + 0.2e1 * (-t101 * t42 - t123 * t41 + t56 * t91 - t57 * t92) * mrSges(7,3); m(4) * t137 + t101 * t15 - t123 * t14 - t91 * t39 - t92 * t40 - t177 * t116 + t178 * t114 + m(7) * (t1 * t101 - t123 * t2 - t3 * t92 - t4 * t91) + m(6) * (t10 * t114 - t11 * t116) + m(5) * (t114 * t17 + t116 * t16) + t133; m(7) * (t101 * t41 - t123 * t42 - t56 * t92 - t57 * t91); (-t101 * t91 + t123 * t92) * t172; m(5) * t47 + m(6) * t19 + m(7) * t12 + t45 + t46 - t7; (-m(6) - m(7)) * t140 - t59; 0; 0; t117 * t15 + t119 * t14 + (-t117 * t40 + t119 * t39) * qJD(6) + m(7) * (t117 * t1 + t119 * t2 + (-t117 * t3 + t119 * t4) * qJD(6)) + m(6) * t11 + t54; m(7) * (t117 * t41 + t119 * t42 + (-t117 * t56 + t119 * t57) * qJD(6)) + m(6) * t142 + (-t117 * t92 + t119 * t91 + (t101 * t117 - t119 * t123) * qJD(6)) * mrSges(7,3); m(7) * (-t117 * t91 - t119 * t92 + (t101 * t119 + t117 * t123) * qJD(6)); 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t138; mrSges(7,1) * t42 - mrSges(7,2) * t41 + t161; -t59; 0; (-mrSges(7,1) * t117 - mrSges(7,2) * t119) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
