% Calculate time derivative of joint inertia matrix for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:13
% EndTime: 2019-03-09 02:00:17
% DurationCPUTime: 1.70s
% Computational Cost: add. (2279->278), mult. (4945->387), div. (0->0), fcn. (4472->8), ass. (0->124)
t128 = Ifges(7,4) + Ifges(6,5);
t159 = Ifges(7,2) + Ifges(6,3);
t87 = sin(qJ(5));
t89 = cos(qJ(5));
t118 = t87 ^ 2 + t89 ^ 2;
t116 = qJD(5) * t89;
t84 = sin(pkin(10));
t85 = cos(pkin(10));
t88 = sin(qJ(4));
t90 = cos(qJ(4));
t59 = t84 * t88 - t90 * t85;
t55 = t59 * qJD(4);
t130 = t87 * t55;
t60 = t84 * t90 + t88 * t85;
t94 = t60 * t116 - t130;
t96 = -cos(pkin(9)) * pkin(1) - pkin(3) * t85 - pkin(2);
t32 = pkin(4) * t59 - pkin(8) * t60 + t96;
t74 = sin(pkin(9)) * pkin(1) + qJ(3);
t142 = pkin(7) + t74;
t57 = t142 * t84;
t58 = t142 * t85;
t35 = -t88 * t57 + t58 * t90;
t156 = t87 * t32 + t89 * t35;
t158 = qJD(5) * t156;
t155 = -t90 * t57 - t58 * t88;
t23 = -t59 * qJD(3) + t155 * qJD(4);
t143 = pkin(8) * t55;
t56 = t60 * qJD(4);
t144 = pkin(4) * t56;
t38 = t143 + t144;
t4 = -t23 * t87 + t38 * t89 - t158;
t2 = -pkin(5) * t56 - t4;
t157 = m(7) * t2;
t117 = qJD(5) * t87;
t154 = Ifges(7,6) * t117 + t128 * t116;
t101 = pkin(5) * t89 + qJ(6) * t87;
t115 = qJD(6) * t89;
t153 = t101 * qJD(5) - t115;
t131 = t60 * t87;
t39 = -mrSges(6,2) * t59 - mrSges(6,3) * t131;
t42 = -mrSges(7,2) * t131 + mrSges(7,3) * t59;
t123 = t39 + t42;
t113 = t60 * t117;
t135 = t55 * t89;
t93 = t113 + t135;
t19 = mrSges(6,1) * t56 + t93 * mrSges(6,3);
t129 = t89 * mrSges(7,2);
t20 = -t56 * mrSges(7,1) - mrSges(7,2) * t113 - t55 * t129;
t127 = t19 - t20;
t152 = -t123 * qJD(5) - t127;
t40 = -mrSges(6,3) * t60 * t89 + mrSges(6,1) * t59;
t41 = -mrSges(7,1) * t59 + t60 * t129;
t122 = -t40 + t41;
t21 = -mrSges(6,2) * t56 - t94 * mrSges(6,3);
t22 = -t94 * mrSges(7,2) + mrSges(7,3) * t56;
t126 = t21 + t22;
t13 = t32 * t89 - t35 * t87;
t3 = t32 * t116 - t35 * t117 + t89 * t23 + t87 * t38;
t151 = t122 * qJD(5) + m(6) * (-t13 * qJD(5) + t3) + t126;
t150 = 2 * m(5);
t149 = -2 * mrSges(5,3);
t148 = -2 * Ifges(5,4);
t24 = t60 * qJD(3) + t35 * qJD(4);
t147 = 0.2e1 * t24;
t146 = -0.2e1 * t155;
t141 = Ifges(6,4) * t87;
t140 = Ifges(6,4) * t89;
t139 = Ifges(7,5) * t87;
t138 = Ifges(7,5) * t89;
t137 = t24 * t155;
t136 = t55 * t60;
t134 = t56 * t59;
t132 = t59 * Ifges(6,6);
t102 = Ifges(7,3) * t87 + t138;
t25 = Ifges(7,6) * t59 + t102 * t60;
t103 = -Ifges(6,2) * t87 + t140;
t26 = t103 * t60 + t132;
t125 = t25 - t26;
t104 = Ifges(7,1) * t89 + t139;
t27 = Ifges(7,4) * t59 + t104 * t60;
t105 = Ifges(6,1) * t89 - t141;
t28 = Ifges(6,5) * t59 + t105 * t60;
t124 = t27 + t28;
t121 = t118 * t143;
t69 = -t89 * mrSges(6,1) + t87 * mrSges(6,2);
t120 = t69 - mrSges(5,1);
t70 = -Ifges(7,3) * t89 + t139;
t71 = Ifges(6,2) * t89 + t141;
t111 = t70 / 0.2e1 - t71 / 0.2e1;
t72 = Ifges(7,1) * t87 - t138;
t73 = Ifges(6,1) * t87 + t140;
t110 = t73 / 0.2e1 + t72 / 0.2e1;
t51 = t55 * mrSges(5,2);
t109 = t56 * mrSges(5,1) - t51;
t107 = t87 * mrSges(6,1) + t89 * mrSges(6,2);
t68 = -t89 * mrSges(7,1) - t87 * mrSges(7,3);
t106 = t87 * mrSges(7,1) - t89 * mrSges(7,3);
t100 = pkin(5) * t87 - qJ(6) * t89;
t98 = -t155 * t56 + t24 * t59;
t95 = t94 * Ifges(7,6) - t128 * t135 + t159 * t56;
t92 = t122 * t87 + t123 * t89;
t91 = m(7) * t115 + (-m(7) * t101 + t68 + t69) * qJD(5);
t67 = -pkin(4) - t101;
t66 = t105 * qJD(5);
t65 = t104 * qJD(5);
t64 = t103 * qJD(5);
t63 = t102 * qJD(5);
t62 = t107 * qJD(5);
t61 = t106 * qJD(5);
t54 = -pkin(5) * t117 + qJ(6) * t116 + qJD(6) * t87;
t37 = t107 * t60;
t36 = t106 * t60;
t17 = t100 * t60 - t155;
t16 = t94 * mrSges(6,1) - t93 * mrSges(6,2);
t15 = t94 * mrSges(7,1) + t93 * mrSges(7,3);
t11 = -t93 * Ifges(6,1) - t94 * Ifges(6,4) + t56 * Ifges(6,5);
t10 = -t93 * Ifges(7,1) + t56 * Ifges(7,4) + t94 * Ifges(7,5);
t9 = -t93 * Ifges(6,4) - t94 * Ifges(6,2) + t56 * Ifges(6,6);
t8 = -t93 * Ifges(7,5) + t56 * Ifges(7,6) + t94 * Ifges(7,3);
t7 = -pkin(5) * t59 - t13;
t6 = qJ(6) * t59 + t156;
t5 = -t100 * t55 + t153 * t60 + t24;
t1 = qJ(6) * t56 + qJD(6) * t59 + t3;
t12 = [0.2e1 * t96 * t109 + t16 * t146 + 0.2e1 * t5 * t36 + t37 * t147 + 0.2e1 * t3 * t39 + 0.2e1 * t4 * t40 + 0.2e1 * t2 * t41 + 0.2e1 * t1 * t42 + 0.2e1 * t13 * t19 + 0.2e1 * t7 * t20 + 0.2e1 * t156 * t21 + 0.2e1 * t6 * t22 + 0.2e1 * t17 * t15 + t35 * t56 * t149 - (mrSges(5,3) * t146 + t124 * t89 + t125 * t87) * t55 + 0.2e1 * m(6) * (t13 * t4 + t156 * t3 - t137) + (t23 * t35 - t137) * t150 + 0.2e1 * m(7) * (t1 * t6 + t17 * t5 + t2 * t7) + (t23 * t149 - (-Ifges(6,6) * t87 + t148) * t55 + ((2 * Ifges(5,2)) + t159) * t56 + t95) * t59 + (mrSges(5,3) * t147 - 0.2e1 * Ifges(5,1) * t55 + (t10 + t11) * t89 + (t8 - t9) * t87 + (t148 + t128 * t89 + (-Ifges(6,6) + Ifges(7,6)) * t87) * t56 + ((t125 - t132) * t89 + (-t128 * t59 - t124) * t87) * qJD(5)) * t60 + 0.2e1 * (m(4) * t74 + mrSges(4,3)) * qJD(3) * (t84 ^ 2 + t85 ^ 2); (t15 + t16) * t59 + (t36 + t37) * t56 - t92 * t55 + m(6) * (t13 * t130 - t135 * t156 + t98) + m(7) * (-t7 * t130 - t6 * t135 + t17 * t56 + t5 * t59) + m(5) * (-t35 * t55 + t98) + (-m(6) * t156 * t117 + m(7) * (t7 * t116 - t6 * t117) + m(5) * t23 + (-m(6) * t4 + t152 + t157) * t87 + (m(7) * t1 + t151) * t89) * t60; (t134 - t136) * t150 + 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (-t118 * t136 + t134); t127 * t89 + t126 * t87 + t92 * qJD(5) + m(7) * (t1 * t87 - t2 * t89 + (t6 * t89 + t7 * t87) * qJD(5)) + m(6) * (t3 * t87 + t4 * t89 + (-t13 * t87 + t156 * t89) * qJD(5)) + t109; 0; 0; t5 * t68 + t17 * t61 - t155 * t62 + t67 * t15 - Ifges(5,6) * t56 - t54 * t36 - Ifges(5,5) * t55 - t23 * mrSges(5,2) - pkin(4) * t16 + m(7) * (-t17 * t54 + t5 * t67) + (-m(6) * pkin(4) + t120) * t24 + (-t4 * mrSges(6,3) + t2 * mrSges(7,2) + t10 / 0.2e1 + t11 / 0.2e1 + (t63 / 0.2e1 - t64 / 0.2e1) * t60 + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t56 - t111 * t55 + (t25 / 0.2e1 - t26 / 0.2e1 - t132 / 0.2e1 - t6 * mrSges(7,2) - t156 * mrSges(6,3) - t110 * t60) * qJD(5) + (m(6) * (-t4 - t158) + m(7) * (-qJD(5) * t6 + t2) + t152) * pkin(8)) * t87 + (t1 * mrSges(7,2) + t3 * mrSges(6,3) - t8 / 0.2e1 + t9 / 0.2e1 + (t65 / 0.2e1 + t66 / 0.2e1) * t60 + (Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t56 - t110 * t55 + (t27 / 0.2e1 + t28 / 0.2e1 + t7 * mrSges(7,2) - t13 * mrSges(6,3) + t111 * t60) * qJD(5) + (m(7) * (t7 * qJD(5) + t1) + t151) * pkin(8)) * t89 + t154 * t59 / 0.2e1; t51 + (t61 + t62) * t59 + (t68 + t120) * t56 + m(6) * (-t121 - t144) + m(7) * (-t54 * t59 + t56 * t67 - t121) - (mrSges(7,2) + mrSges(6,3)) * t55 * t118; 0; -0.2e1 * pkin(4) * t62 + 0.2e1 * t61 * t67 + (-t63 + t64) * t89 + (t65 + t66) * t87 + 0.2e1 * (-m(7) * t67 - t68) * t54 + ((t72 + t73) * t89 + (t70 - t71) * t87) * qJD(5); Ifges(6,6) * t130 + qJD(6) * t42 + qJ(6) * t22 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6) + t1 * mrSges(7,3) - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - pkin(5) * t20 + (-Ifges(6,6) * t89 - t128 * t87) * t60 * qJD(5) + t95; -(-m(7) * t100 - t106 - t107) * t55 + t91 * t60; m(7) * t54 + ((-mrSges(6,2) + mrSges(7,3)) * t89 + (-mrSges(6,1) - mrSges(7,1)) * t87) * qJD(5); -mrSges(7,2) * t153 - Ifges(6,6) * t117 + t91 * pkin(8) + t154; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); t20 + t157; t94 * m(7); m(7) * t117; (m(7) * pkin(8) + mrSges(7,2)) * t116; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
