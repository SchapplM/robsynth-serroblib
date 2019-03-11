% Calculate time derivative of joint inertia matrix for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:23
% EndTime: 2019-03-09 02:15:25
% DurationCPUTime: 1.70s
% Computational Cost: add. (2273->287), mult. (4691->401), div. (0->0), fcn. (4198->6), ass. (0->124)
t140 = Ifges(7,4) + Ifges(6,5);
t170 = Ifges(7,2) + Ifges(6,3);
t92 = sin(qJ(5));
t93 = cos(qJ(5));
t130 = t92 ^ 2 + t93 ^ 2;
t156 = sin(qJ(4));
t157 = cos(qJ(4));
t89 = sin(pkin(9));
t90 = cos(pkin(9));
t57 = t156 * t89 - t157 * t90;
t128 = qJD(5) * t93;
t112 = qJD(4) * t156;
t113 = qJD(4) * t157;
t54 = -t112 * t90 - t113 * t89;
t144 = t92 * t54;
t98 = t57 * t128 - t144;
t129 = qJD(5) * t92;
t143 = t93 * t54;
t97 = t57 * t129 + t143;
t58 = t156 * t90 + t157 * t89;
t78 = t89 * pkin(3) + qJ(2);
t33 = pkin(4) * t58 + pkin(8) * t57 + t78;
t91 = -pkin(1) - qJ(3);
t158 = -pkin(7) + t91;
t65 = t158 * t89;
t66 = t158 * t90;
t42 = t156 * t66 + t157 * t65;
t167 = t92 * t33 + t93 * t42;
t169 = qJD(5) * t167;
t168 = 2 * qJD(2);
t41 = t156 * t65 - t157 * t66;
t166 = Ifges(7,6) * t129 + t140 * t128;
t111 = (t89 ^ 2 + t90 ^ 2) * qJD(3);
t27 = -t58 * qJD(3) - qJD(4) * t41;
t55 = -t112 * t89 + t113 * t90;
t32 = pkin(4) * t55 - pkin(8) * t54 + qJD(2);
t4 = -t92 * t27 + t32 * t93 - t169;
t3 = t33 * t128 - t129 * t42 + t93 * t27 + t92 * t32;
t1 = qJ(6) * t55 + qJD(6) * t58 + t3;
t15 = t33 * t93 - t92 * t42;
t12 = -t58 * pkin(5) - t15;
t145 = t57 * t93;
t38 = t58 * mrSges(6,1) + mrSges(6,3) * t145;
t39 = -t58 * mrSges(7,1) - mrSges(7,2) * t145;
t134 = -t38 + t39;
t20 = -t55 * mrSges(6,2) + mrSges(6,3) * t98;
t21 = mrSges(7,2) * t98 + t55 * mrSges(7,3);
t138 = t20 + t21;
t165 = t134 * qJD(5) + m(7) * (t12 * qJD(5) + t1) + m(6) * (-t15 * qJD(5) + t3) + t138;
t11 = qJ(6) * t58 + t167;
t146 = t57 * t92;
t37 = -mrSges(6,2) * t58 + mrSges(6,3) * t146;
t40 = mrSges(7,2) * t146 + mrSges(7,3) * t58;
t135 = t37 + t40;
t18 = t55 * mrSges(6,1) - mrSges(6,3) * t97;
t19 = -t55 * mrSges(7,1) + mrSges(7,2) * t97;
t139 = t18 - t19;
t2 = -t55 * pkin(5) - t4;
t164 = -t135 * qJD(5) + m(7) * (-t11 * qJD(5) + t2) + m(6) * (-t4 - t169) - t139;
t163 = 2 * m(5);
t162 = -2 * mrSges(5,3);
t161 = 0.2e1 * t41;
t155 = Ifges(6,4) * t92;
t154 = Ifges(6,4) * t93;
t153 = Ifges(7,5) * t92;
t152 = Ifges(7,5) * t93;
t151 = Ifges(6,6) * t58;
t28 = -qJD(3) * t57 + qJD(4) * t42;
t150 = t28 * t41;
t149 = t54 * t57;
t148 = t55 * t58;
t147 = t55 * t92;
t142 = t93 * t55;
t69 = -t93 * mrSges(6,1) + t92 * mrSges(6,2);
t141 = mrSges(5,1) - t69;
t105 = Ifges(7,3) * t92 + t152;
t23 = Ifges(7,6) * t58 - t105 * t57;
t106 = -Ifges(6,2) * t92 + t154;
t24 = -t106 * t57 + t151;
t137 = t23 - t24;
t107 = Ifges(7,1) * t93 + t153;
t25 = Ifges(7,4) * t58 - t107 * t57;
t108 = Ifges(6,1) * t93 - t155;
t26 = Ifges(6,5) * t58 - t108 * t57;
t136 = t25 + t26;
t133 = t130 * pkin(8) * t55;
t127 = qJD(6) * t93;
t122 = t55 * t162;
t70 = -Ifges(7,3) * t93 + t153;
t71 = Ifges(6,2) * t93 + t155;
t119 = t70 / 0.2e1 - t71 / 0.2e1;
t72 = Ifges(7,1) * t92 - t152;
t73 = Ifges(6,1) * t92 + t154;
t118 = t72 / 0.2e1 + t73 / 0.2e1;
t114 = t55 * mrSges(5,1) + t54 * mrSges(5,2);
t110 = mrSges(6,1) * t92 + t93 * mrSges(6,2);
t68 = -t93 * mrSges(7,1) - t92 * mrSges(7,3);
t109 = mrSges(7,1) * t92 - t93 * mrSges(7,3);
t104 = -pkin(5) * t93 - t92 * qJ(6);
t103 = pkin(5) * t92 - qJ(6) * t93;
t101 = t57 * t28 - t54 * t41;
t99 = -t98 * Ifges(7,6) + t140 * t143 + t170 * t55;
t96 = t134 * t92 + t135 * t93;
t95 = qJD(5) * t104 + t127;
t94 = m(7) * t127 + (m(7) * t104 + t68 + t69) * qJD(5);
t67 = -pkin(4) + t104;
t64 = t108 * qJD(5);
t63 = t107 * qJD(5);
t62 = t106 * qJD(5);
t61 = t105 * qJD(5);
t60 = t110 * qJD(5);
t59 = t109 * qJD(5);
t53 = -pkin(5) * t129 + qJ(6) * t128 + t92 * qJD(6);
t35 = t110 * t57;
t34 = t109 * t57;
t17 = -t103 * t57 + t41;
t14 = -mrSges(6,1) * t98 + mrSges(6,2) * t97;
t13 = -mrSges(7,1) * t98 - mrSges(7,3) * t97;
t9 = Ifges(6,1) * t97 + Ifges(6,4) * t98 + Ifges(6,5) * t55;
t8 = Ifges(7,1) * t97 + Ifges(7,4) * t55 - Ifges(7,5) * t98;
t7 = Ifges(6,4) * t97 + Ifges(6,2) * t98 + Ifges(6,6) * t55;
t6 = Ifges(7,5) * t97 + Ifges(7,6) * t55 - Ifges(7,3) * t98;
t5 = t103 * t54 + t57 * t95 + t28;
t10 = [0.2e1 * t78 * t114 + t42 * t122 - 0.2e1 * t5 * t34 - 0.2e1 * t28 * t35 + 0.2e1 * t3 * t37 + 0.2e1 * t4 * t38 + 0.2e1 * t2 * t39 + 0.2e1 * t1 * t40 + t14 * t161 + 0.2e1 * t15 * t18 + 0.2e1 * t12 * t19 + 0.2e1 * t167 * t20 + 0.2e1 * t11 * t21 + 0.2e1 * t17 * t13 + (mrSges(5,3) * t161 + t136 * t93 + t137 * t92) * t54 + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t111 * t91) + (qJD(2) * t78 + t27 * t42 + t150) * t163 + 0.2e1 * m(6) * (t15 * t4 + t167 * t3 + t150) + 0.2e1 * m(7) * (t1 * t11 + t12 * t2 + t17 * t5) + (mrSges(5,1) * t168 + t27 * t162 + (-Ifges(6,6) * t92 - (2 * Ifges(5,4))) * t54 + ((2 * Ifges(5,2)) + t170) * t55 + t99) * t58 + (-0.2e1 * qJD(2) * mrSges(5,2) + t28 * t162 - 0.2e1 * Ifges(5,1) * t54 + (-t8 - t9) * t93 + (-t6 + t7) * t92 + ((2 * Ifges(5,4)) - t140 * t93 + (Ifges(6,6) - Ifges(7,6)) * t92) * t55 + ((-t137 + t151) * t93 + (t140 * t58 + t136) * t92) * qJD(5)) * t57 + (m(3) * qJ(2) + mrSges(4,1) * t89 + mrSges(4,2) * t90 + mrSges(3,3)) * t168 + 0.2e1 * mrSges(4,3) * t111; (t13 + t14) * t57 + (0.2e1 * t57 * mrSges(5,3) + t34 + t35) * t54 + t96 * t55 + m(7) * (t11 * t142 + t12 * t147 - t54 * t17 + t57 * t5) + m(6) * (t142 * t167 - t15 * t147 + t101) + m(5) * (t42 * t55 + t101) - m(4) * t111 + (m(5) * t27 + t164 * t92 + t165 * t93 + t122) * t58; (t148 - t149) * t163 + 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (t130 * t148 - t149); t139 * t93 + t138 * t92 + (m(5) + m(4)) * qJD(2) + t96 * qJD(5) + m(7) * (t92 * t1 - t2 * t93 + (t11 * t93 + t12 * t92) * qJD(5)) + m(6) * (t92 * t3 + t4 * t93 + (-t15 * t92 + t167 * t93) * qJD(5)) + t114; 0; 0; t17 * t59 + t41 * t60 + t67 * t13 + t5 * t68 + t53 * t34 + Ifges(5,5) * t54 - Ifges(5,6) * t55 - t27 * mrSges(5,2) - pkin(4) * t14 + m(7) * (-t53 * t17 + t67 * t5) + (-m(6) * pkin(4) - t141) * t28 + (t8 / 0.2e1 + t9 / 0.2e1 - t4 * mrSges(6,3) + t2 * mrSges(7,2) + (-t61 / 0.2e1 + t62 / 0.2e1) * t57 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t55 + t119 * t54 + (-t151 / 0.2e1 - t24 / 0.2e1 + t23 / 0.2e1 - t11 * mrSges(7,2) - t167 * mrSges(6,3) + t118 * t57) * qJD(5) + t164 * pkin(8)) * t92 + (-t6 / 0.2e1 + t7 / 0.2e1 + t3 * mrSges(6,3) + t1 * mrSges(7,2) + (-t63 / 0.2e1 - t64 / 0.2e1) * t57 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t55 + t118 * t54 + (t25 / 0.2e1 + t26 / 0.2e1 + t12 * mrSges(7,2) - t15 * mrSges(6,3) - t119 * t57) * qJD(5) + t165 * pkin(8)) * t93 + t166 * t58 / 0.2e1; (t59 + t60) * t57 + (-t68 + t141) * t54 + m(7) * (-t53 * t57 - t54 * t67 + t133) + m(6) * (pkin(4) * t54 + t133) + (-mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t130) * t55; 0; -0.2e1 * pkin(4) * t60 + 0.2e1 * t67 * t59 + (-t61 + t62) * t93 + (t63 + t64) * t92 + 0.2e1 * (-m(7) * t67 - t68) * t53 + ((t72 + t73) * t93 + (t70 - t71) * t92) * qJD(5); -Ifges(6,6) * t144 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t11) + qJD(6) * t40 + qJ(6) * t21 + t1 * mrSges(7,3) - pkin(5) * t19 + t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + (Ifges(6,6) * t93 + t140 * t92) * t57 * qJD(5) + t99; (-m(7) * t103 - t109 - t110) * t55 + t94 * t58; m(7) * t53 + ((-mrSges(6,2) + mrSges(7,3)) * t93 + (-mrSges(6,1) - mrSges(7,1)) * t92) * qJD(5); mrSges(7,2) * t95 - Ifges(6,6) * t129 + pkin(8) * t94 + t166; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t19; (t128 * t58 + t147) * m(7); m(7) * t129; (m(7) * pkin(8) + mrSges(7,2)) * t128; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
