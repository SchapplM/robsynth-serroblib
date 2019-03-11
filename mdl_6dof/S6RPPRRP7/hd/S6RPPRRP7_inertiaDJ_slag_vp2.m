% Calculate time derivative of joint inertia matrix for
% S6RPPRRP7
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:43
% EndTime: 2019-03-09 02:12:46
% DurationCPUTime: 1.76s
% Computational Cost: add. (2275->314), mult. (4753->440), div. (0->0), fcn. (4253->6), ass. (0->137)
t142 = Ifges(6,5) + Ifges(7,5);
t173 = Ifges(6,3) + Ifges(7,3);
t141 = Ifges(6,6) + Ifges(7,6);
t95 = sin(qJ(5));
t117 = t141 * t95;
t172 = -(2 * Ifges(5,4)) - t117;
t160 = sin(qJ(4));
t161 = cos(qJ(4));
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t58 = t160 * t92 - t161 * t93;
t96 = cos(qJ(5));
t129 = qJD(5) * t96;
t125 = t58 * t129;
t112 = qJD(4) * t160;
t113 = qJD(4) * t161;
t55 = -t112 * t93 - t113 * t92;
t100 = -t95 * t55 + t125;
t130 = qJD(5) * t95;
t126 = t58 * t130;
t145 = t96 * t55;
t99 = t126 + t145;
t171 = 2 * qJD(2);
t59 = t160 * t93 + t161 * t92;
t80 = t92 * pkin(3) + qJ(2);
t35 = pkin(4) * t59 + pkin(8) * t58 + t80;
t94 = -pkin(1) - qJ(3);
t162 = -pkin(7) + t94;
t66 = t162 * t92;
t67 = t162 * t93;
t44 = t160 * t67 + t161 * t66;
t38 = t96 * t44;
t16 = t95 * t35 + t38;
t170 = -t160 * t66 + t161 * t67;
t140 = -qJ(6) - pkin(8);
t68 = t140 * t95;
t71 = t140 * t96;
t169 = -t68 * t95 - t71 * t96;
t110 = (t92 ^ 2 + t93 ^ 2) * qJD(3);
t168 = 2 * m(5);
t167 = 2 * m(7);
t166 = -2 * mrSges(5,3);
t165 = -2 * mrSges(7,3);
t164 = -0.2e1 * t170;
t163 = m(7) * pkin(5);
t159 = mrSges(6,2) * t96;
t158 = Ifges(6,4) * t95;
t157 = Ifges(6,4) * t96;
t156 = Ifges(7,4) * t95;
t155 = Ifges(7,4) * t96;
t27 = -t58 * qJD(3) + qJD(4) * t44;
t154 = t27 * t170;
t153 = t55 * t58;
t56 = -t112 * t92 + t113 * t93;
t152 = t56 * t95;
t151 = t56 * t96;
t150 = t58 * t95;
t149 = t58 * t96;
t144 = mrSges(6,1) * t96 - t95 * mrSges(6,2) + mrSges(5,1);
t143 = mrSges(6,2) + mrSges(7,2);
t17 = t56 * mrSges(7,1) - mrSges(7,3) * t99;
t18 = t56 * mrSges(6,1) - mrSges(6,3) * t99;
t139 = t17 + t18;
t19 = -t56 * mrSges(7,2) + mrSges(7,3) * t100;
t20 = -t56 * mrSges(6,2) + mrSges(6,3) * t100;
t138 = t19 + t20;
t103 = -Ifges(7,2) * t95 + t155;
t22 = Ifges(7,6) * t59 - t103 * t58;
t104 = -Ifges(6,2) * t95 + t157;
t23 = Ifges(6,6) * t59 - t104 * t58;
t137 = t22 + t23;
t105 = Ifges(7,1) * t96 - t156;
t24 = Ifges(7,5) * t59 - t105 * t58;
t106 = Ifges(6,1) * t96 - t158;
t25 = Ifges(6,5) * t59 - t106 * t58;
t136 = t24 + t25;
t39 = -mrSges(7,2) * t59 + mrSges(7,3) * t150;
t40 = -mrSges(6,2) * t59 + mrSges(6,3) * t150;
t135 = t39 + t40;
t41 = t59 * mrSges(7,1) + mrSges(7,3) * t149;
t42 = t59 * mrSges(6,1) + mrSges(6,3) * t149;
t134 = -t41 - t42;
t60 = mrSges(7,1) * t130 + mrSges(7,2) * t129;
t132 = t95 ^ 2 + t96 ^ 2;
t131 = qJ(6) * t58;
t128 = t56 * t166;
t26 = -t59 * qJD(3) + qJD(4) * t170;
t34 = pkin(4) * t56 - pkin(8) * t55 + qJD(2);
t127 = t35 * t129 + t96 * t26 + t95 * t34;
t124 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t72 = Ifges(7,2) * t96 + t156;
t73 = Ifges(6,2) * t96 + t158;
t123 = t72 / 0.2e1 + t73 / 0.2e1;
t74 = Ifges(7,1) * t95 + t155;
t75 = Ifges(6,1) * t95 + t157;
t122 = t74 / 0.2e1 + t75 / 0.2e1;
t121 = mrSges(7,1) + t163;
t116 = t132 * t56;
t115 = t56 * mrSges(5,1) + t55 * mrSges(5,2);
t114 = -t95 * t26 + t96 * t34;
t15 = t96 * t35 - t95 * t44;
t111 = qJD(5) * t140;
t109 = t142 * t145 + t173 * t56;
t108 = -mrSges(6,1) - t121;
t107 = mrSges(6,1) * t95 + t159;
t102 = t170 * t55 + t58 * t27;
t101 = -qJ(6) * t55 + qJD(6) * t58;
t13 = -t100 * mrSges(7,1) + t99 * mrSges(7,2);
t98 = t141 * t96 + t142 * t95;
t97 = t134 * t95 + t135 * t96;
t86 = Ifges(6,5) * t129;
t85 = Ifges(7,5) * t129;
t82 = -pkin(5) * t96 - pkin(4);
t69 = -mrSges(7,1) * t96 + t95 * mrSges(7,2);
t65 = t106 * qJD(5);
t64 = t105 * qJD(5);
t63 = t104 * qJD(5);
t62 = t103 * qJD(5);
t61 = t107 * qJD(5);
t54 = -t95 * qJD(6) + t111 * t96;
t53 = qJD(6) * t96 + t111 * t95;
t37 = t107 * t58;
t36 = (-mrSges(7,1) * t95 - mrSges(7,2) * t96) * t58;
t28 = -pkin(5) * t150 - t170;
t14 = -mrSges(6,1) * t100 + mrSges(6,2) * t99;
t12 = -pkin(5) * t100 + t27;
t11 = t131 * t95 + t16;
t9 = Ifges(6,1) * t99 + Ifges(6,4) * t100 + Ifges(6,5) * t56;
t8 = Ifges(7,1) * t99 + Ifges(7,4) * t100 + Ifges(7,5) * t56;
t7 = Ifges(6,4) * t99 + Ifges(6,2) * t100 + Ifges(6,6) * t56;
t6 = Ifges(7,4) * t99 + Ifges(7,2) * t100 + Ifges(7,6) * t56;
t5 = t59 * pkin(5) + t131 * t96 + t15;
t4 = -qJD(5) * t16 + t114;
t3 = -t130 * t44 + t127;
t2 = qJ(6) * t125 + (-qJD(5) * t44 + t101) * t95 + t127;
t1 = t56 * pkin(5) + t101 * t96 + (-t38 + (-t35 - t131) * t95) * qJD(5) + t114;
t10 = [t44 * t128 + 0.2e1 * t80 * t115 + 0.2e1 * t2 * t39 + 0.2e1 * t3 * t40 + 0.2e1 * t1 * t41 + 0.2e1 * t4 * t42 + t14 * t164 + 0.2e1 * t28 * t13 + 0.2e1 * t12 * t36 - 0.2e1 * t27 * t37 + 0.2e1 * t5 * t17 + 0.2e1 * t15 * t18 + 0.2e1 * t11 * t19 + 0.2e1 * t16 * t20 + (mrSges(5,3) * t164 + t136 * t96 - t137 * t95) * t55 + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t110 * t94) + (qJD(2) * t80 + t26 * t44 - t154) * t168 + 0.2e1 * m(6) * (t15 * t4 + t16 * t3 - t154) + (t1 * t5 + t11 * t2 + t12 * t28) * t167 + (mrSges(5,1) * t171 + t26 * t166 + ((2 * Ifges(5,2)) + t173) * t56 + t172 * t55 + t109) * t59 + (-0.2e1 * qJD(2) * mrSges(5,2) + t27 * t166 - 0.2e1 * Ifges(5,1) * t55 + (-t8 - t9) * t96 + (t6 + t7) * t95 + (-t142 * t96 - t172) * t56 + (t136 * t95 + t137 * t96 + t59 * t98) * qJD(5)) * t58 + (m(3) * qJ(2) + mrSges(4,1) * t92 + mrSges(4,2) * t93 + mrSges(3,3)) * t171 + 0.2e1 * mrSges(4,3) * t110; (t13 + t14) * t58 + (0.2e1 * t58 * mrSges(5,3) - t36 + t37) * t55 + t97 * t56 + m(6) * (-t15 * t152 + t16 * t151 + t102) + m(7) * (t11 * t151 + t58 * t12 - t5 * t152 - t55 * t28) + m(5) * (t44 * t56 + t102) - m(4) * t110 + (t128 + t138 * t96 - t139 * t95 + (t134 * t96 - t135 * t95) * qJD(5) + m(6) * (-t129 * t15 - t130 * t16 + t3 * t96 - t4 * t95) + m(7) * (-t1 * t95 - t11 * t130 - t129 * t5 + t2 * t96) + m(5) * t26) * t59; (t56 * t59 - t153) * t168 + 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (t116 * t59 - t153); t139 * t96 + t138 * t95 + (m(5) + m(4)) * qJD(2) + t97 * qJD(5) + m(6) * (t95 * t3 + t4 * t96 + (-t15 * t95 + t16 * t96) * qJD(5)) + m(7) * (t1 * t96 + t95 * t2 + (t11 * t96 - t5 * t95) * qJD(5)) + t115; 0; 0; m(7) * (t1 * t68 + t11 * t53 + t12 * t82 - t2 * t71 + t5 * t54) + t68 * t17 + t12 * t69 - t71 * t19 + t82 * t13 + t28 * t60 - t170 * t61 - Ifges(5,6) * t56 + t53 * t39 + t54 * t41 + Ifges(5,5) * t55 - t26 * mrSges(5,2) - pkin(4) * t14 + (t85 / 0.2e1 + t86 / 0.2e1) * t59 + (-m(6) * pkin(4) - t144) * t27 + (-t4 * mrSges(6,3) - t1 * mrSges(7,3) + t8 / 0.2e1 + t9 / 0.2e1 + (t62 / 0.2e1 + t63 / 0.2e1) * t58 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t56 - t123 * t55 + (pkin(5) * t36 - t16 * mrSges(6,3) - t11 * mrSges(7,3) - t22 / 0.2e1 - t23 / 0.2e1 - t124 * t59 + t122 * t58 + t28 * t163) * qJD(5) + (-m(6) * t4 - t18 + (-m(6) * t16 - t40) * qJD(5)) * pkin(8)) * t95 + (t3 * mrSges(6,3) + t2 * mrSges(7,3) + t6 / 0.2e1 + t7 / 0.2e1 + (-t64 / 0.2e1 - t65 / 0.2e1) * t58 + t124 * t56 + t122 * t55 + (m(6) * t3 + t20) * pkin(8) + (-t15 * mrSges(6,3) - t5 * mrSges(7,3) + t24 / 0.2e1 + t25 / 0.2e1 + t123 * t58 + (-m(6) * t15 - t42) * pkin(8)) * qJD(5)) * t96; (t60 + t61) * t58 + (-t69 + t144) * t55 + (-mrSges(5,2) + (mrSges(6,3) + mrSges(7,3)) * t132) * t56 + m(6) * (pkin(4) * t55 + pkin(8) * t116) + (pkin(5) * t126 - t82 * t55 + (t53 * t96 - t54 * t95 + (-t68 * t96 + t71 * t95) * qJD(5)) * t59 + t169 * t56) * m(7); m(7) * (qJD(5) * t169 + t53 * t95 + t54 * t96); -0.2e1 * pkin(4) * t61 + 0.2e1 * t82 * t60 + (-t53 * t71 + t54 * t68) * t167 + (t54 * t165 + t64 + t65 + (-t71 * t165 - t72 - t73 + 0.2e1 * (m(7) * t82 + t69) * pkin(5)) * qJD(5)) * t95 + (0.2e1 * t53 * mrSges(7,3) + t62 + t63 + (t165 * t68 + t74 + t75) * qJD(5)) * t96; t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2) - t55 * t117 + (m(7) * t1 + t17) * pkin(5) + t98 * t58 * qJD(5) + t109; (t108 * t95 - t143 * t96) * t56 + (t108 * t96 + t143 * t95) * t59 * qJD(5); (-t159 + (-mrSges(6,1) - t163) * t95) * qJD(5) - t60; -t53 * mrSges(7,2) + t85 + t86 + t121 * t54 + ((-mrSges(6,1) * pkin(8) - (mrSges(7,3) * pkin(5))) * t96 + (mrSges(6,2) * pkin(8) - t141) * t95) * qJD(5); 0; m(7) * t12 + t13; -m(7) * t55; 0; t130 * t163 + t60; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
