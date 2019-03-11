% Calculate time derivative of joint inertia matrix for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:36
% EndTime: 2019-03-08 20:58:39
% DurationCPUTime: 2.17s
% Computational Cost: add. (3234->330), mult. (8029->518), div. (0->0), fcn. (7824->12), ass. (0->148)
t178 = 2 * mrSges(5,3);
t113 = sin(pkin(12));
t116 = cos(pkin(12));
t119 = sin(qJ(3));
t140 = qJD(3) * t119;
t136 = pkin(3) * t140;
t114 = sin(pkin(11));
t122 = cos(qJ(3));
t147 = cos(pkin(11));
t131 = t147 * t119;
t99 = t114 * t122 + t131;
t89 = t99 * qJD(3);
t146 = t114 * t119;
t126 = t122 * t147 - t146;
t90 = t126 * qJD(3);
t40 = pkin(4) * t89 - qJ(5) * t90 - qJD(5) * t99 + t136;
t161 = -qJ(4) - pkin(8);
t132 = qJD(3) * t161;
t125 = -qJD(4) * t119 + t122 * t132;
t88 = qJD(4) * t122 + t119 * t132;
t51 = t114 * t125 + t147 * t88;
t16 = -t113 * t51 + t116 * t40;
t17 = t113 * t40 + t116 * t51;
t177 = -t16 * t113 + t17 * t116;
t115 = sin(pkin(6));
t120 = sin(qJ(2));
t142 = qJD(2) * t120;
t134 = t115 * t142;
t123 = cos(qJ(2));
t141 = qJD(2) * t123;
t133 = t115 * t141;
t117 = cos(pkin(6));
t145 = t115 * t120;
t94 = t117 * t119 + t122 * t145;
t76 = -qJD(3) * t94 - t119 * t133;
t93 = t117 * t122 - t119 * t145;
t77 = qJD(3) * t93 + t122 * t133;
t34 = t114 * t76 + t147 * t77;
t29 = -t113 * t34 + t116 * t134;
t30 = t113 * t134 + t116 * t34;
t176 = -t113 * t29 + t116 * t30;
t118 = sin(qJ(6));
t121 = cos(qJ(6));
t127 = t113 * t118 - t116 * t121;
t91 = t127 * qJD(6);
t170 = m(5) * pkin(3);
t175 = t114 * t170 - mrSges(5,2);
t135 = t147 * pkin(3);
t108 = -t135 - pkin(4);
t174 = m(6) * t108 - mrSges(6,1) * t116 + mrSges(6,2) * t113 - t147 * t170 - mrSges(5,1);
t173 = 0.2e1 * m(6);
t172 = 2 * m(7);
t50 = t114 * t88 - t147 * t125;
t171 = 0.2e1 * t50;
t112 = t116 ^ 2;
t151 = t116 * t90;
t154 = t113 * t90;
t49 = mrSges(6,1) * t154 + mrSges(6,2) * t151;
t100 = t113 * t121 + t116 * t118;
t92 = t100 * qJD(6);
t26 = -t127 * t90 - t92 * t99;
t27 = -t100 * t90 + t91 * t99;
t9 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t168 = t49 + t9;
t167 = pkin(3) * t114;
t105 = t161 * t122;
t78 = -t105 * t114 - t161 * t131;
t166 = t50 * t78;
t33 = t114 * t77 - t147 * t76;
t54 = t114 * t94 - t147 * t93;
t18 = t54 * t33;
t165 = t91 * mrSges(7,3);
t164 = t92 * mrSges(7,3);
t163 = t127 * mrSges(7,3);
t106 = qJ(5) + t167;
t162 = pkin(9) + t106;
t109 = -pkin(3) * t122 - pkin(2);
t67 = -pkin(4) * t126 - qJ(5) * t99 + t109;
t79 = -t105 * t147 + t146 * t161;
t32 = t113 * t67 + t116 * t79;
t160 = -Ifges(7,5) * t91 - Ifges(7,6) * t92;
t159 = Ifges(6,4) * t113;
t158 = Ifges(6,4) * t116;
t157 = t100 * mrSges(7,3);
t156 = t113 * Ifges(6,2);
t153 = t113 * t99;
t150 = t116 * t99;
t144 = t115 * t123;
t139 = qJD(3) * t122;
t138 = 0.2e1 * t119;
t137 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t89;
t63 = t89 * mrSges(5,1) + t90 * mrSges(5,2);
t31 = -t113 * t79 + t116 * t67;
t130 = t115 ^ 2 * t120 * t141;
t129 = t78 * t33 + t50 * t54;
t128 = Ifges(6,5) * t116 - Ifges(6,6) * t113;
t19 = -pkin(5) * t126 - pkin(9) * t150 + t31;
t25 = -pkin(9) * t153 + t32;
t5 = -t118 * t25 + t121 * t19;
t6 = t118 * t19 + t121 * t25;
t55 = t114 * t93 + t147 * t94;
t45 = -t113 * t55 - t116 * t144;
t46 = -t113 * t144 + t116 * t55;
t12 = -t118 * t46 + t121 * t45;
t13 = t118 * t45 + t121 * t46;
t95 = t162 * t113;
t96 = t162 * t116;
t60 = -t118 * t96 - t121 * t95;
t61 = -t118 * t95 + t121 * t96;
t124 = -t76 * t119 + t77 * t122 + (-t119 * t94 - t122 * t93) * qJD(3);
t103 = -t116 * pkin(5) + t108;
t102 = (mrSges(4,1) * t119 + mrSges(4,2) * t122) * qJD(3);
t74 = Ifges(7,1) * t100 - Ifges(7,4) * t127;
t73 = Ifges(7,4) * t100 - Ifges(7,2) * t127;
t72 = mrSges(7,1) * t127 + mrSges(7,2) * t100;
t71 = -mrSges(5,1) * t126 + mrSges(5,2) * t99;
t69 = -mrSges(6,1) * t126 - mrSges(6,3) * t150;
t68 = mrSges(6,2) * t126 - mrSges(6,3) * t153;
t66 = -Ifges(7,1) * t91 - Ifges(7,4) * t92;
t65 = -Ifges(7,4) * t91 - Ifges(7,2) * t92;
t64 = mrSges(7,1) * t92 - mrSges(7,2) * t91;
t62 = (mrSges(6,1) * t113 + mrSges(6,2) * t116) * t99;
t59 = mrSges(6,1) * t89 - mrSges(6,3) * t151;
t58 = -mrSges(6,2) * t89 - mrSges(6,3) * t154;
t53 = t127 * t99;
t52 = t100 * t99;
t48 = pkin(5) * t153 + t78;
t44 = -qJD(5) * t100 - qJD(6) * t61;
t43 = -qJD(5) * t127 + qJD(6) * t60;
t42 = -mrSges(7,1) * t126 + mrSges(7,3) * t53;
t41 = mrSges(7,2) * t126 - mrSges(7,3) * t52;
t39 = t89 * Ifges(6,5) + (t116 * Ifges(6,1) - t159) * t90;
t38 = t89 * Ifges(6,6) + (-t156 + t158) * t90;
t35 = pkin(5) * t154 + t50;
t28 = mrSges(7,1) * t52 - mrSges(7,2) * t53;
t21 = -Ifges(7,1) * t53 - Ifges(7,4) * t52 - Ifges(7,5) * t126;
t20 = -Ifges(7,4) * t53 - Ifges(7,2) * t52 - Ifges(7,6) * t126;
t15 = -mrSges(7,2) * t89 + mrSges(7,3) * t27;
t14 = mrSges(7,1) * t89 - mrSges(7,3) * t26;
t11 = -pkin(9) * t154 + t17;
t10 = pkin(5) * t89 - pkin(9) * t151 + t16;
t8 = Ifges(7,1) * t26 + Ifges(7,4) * t27 + t89 * Ifges(7,5);
t7 = Ifges(7,4) * t26 + Ifges(7,2) * t27 + t89 * Ifges(7,6);
t4 = -qJD(6) * t13 - t118 * t30 + t121 * t29;
t3 = qJD(6) * t12 + t118 * t29 + t121 * t30;
t2 = -qJD(6) * t6 + t10 * t121 - t11 * t118;
t1 = qJD(6) * t5 + t10 * t118 + t11 * t121;
t22 = [0.2e1 * m(7) * (t12 * t4 + t13 * t3 + t18) + 0.2e1 * m(5) * (t55 * t34 - t130 + t18) + 0.2e1 * m(6) * (t29 * t45 + t30 * t46 + t18) + 0.2e1 * m(4) * (t93 * t76 + t94 * t77 - t130); t12 * t14 + t13 * t15 + t29 * t69 + t3 * t41 + t30 * t68 + t4 * t42 + t45 * t59 + t46 * t58 + t168 * t54 + (t28 + t62) * t33 + (t126 * t34 + t33 * t99 + t54 * t90 - t55 * t89) * mrSges(5,3) + t124 * mrSges(4,3) + ((-t102 - t63) * t123 + (-t123 * mrSges(3,2) + (-mrSges(4,1) * t122 + mrSges(4,2) * t119 - mrSges(3,1) + t71) * t120) * qJD(2)) * t115 + m(6) * (t16 * t45 + t17 * t46 + t29 * t31 + t30 * t32 + t129) + m(7) * (t1 * t13 + t12 * t2 + t3 * t6 + t33 * t48 + t35 * t54 + t4 * t5) + m(5) * (t79 * t34 + t51 * t55 + (t109 * t142 - t123 * t136) * t115 + t129) + (-pkin(2) * t134 + pkin(8) * t124) * m(4); (-Ifges(4,4) * t119 + pkin(3) * t71) * qJD(3) * t138 + (t78 * t90 - t79 * t89) * t178 + (t51 * t178 - ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t89 - t137 - 0.2e1 * (-Ifges(5,4) + t128) * t90) * t126 + t89 * (-Ifges(7,5) * t53 - Ifges(7,6) * t52) + (0.2e1 * Ifges(4,4) * t122 + (Ifges(4,1) - Ifges(4,2)) * t138) * t139 + t62 * t171 + (t1 * t6 + t2 * t5 + t35 * t48) * t172 + (t31 * t16 + t32 * t17 + t166) * t173 + 0.2e1 * t5 * t14 + 0.2e1 * t6 * t15 + 0.2e1 * m(5) * (t109 * t136 + t51 * t79 + t166) + t26 * t21 + t27 * t20 + (mrSges(5,3) * t171 - t113 * t38 + t116 * t39 + (-0.2e1 * Ifges(5,4) + t128) * t89 + (Ifges(6,1) * t112 + (2 * Ifges(5,1)) + (t156 - 0.2e1 * t158) * t113) * t90) * t99 + 0.2e1 * t35 * t28 + 0.2e1 * t1 * t41 + 0.2e1 * t2 * t42 + 0.2e1 * t48 * t9 - t52 * t7 - t53 * t8 + 0.2e1 * t32 * t58 + 0.2e1 * t31 * t59 + 0.2e1 * t17 * t68 + 0.2e1 * t16 * t69 + 0.2e1 * t78 * t49 - 0.2e1 * pkin(2) * t102 + 0.2e1 * t109 * t63; m(7) * (t12 * t44 + t13 * t43 + t3 * t61 + t4 * t60) - t3 * t163 - t13 * t164 - t4 * t157 + t12 * t165 + t54 * t64 + m(6) * (t176 * t106 + (-t113 * t45 + t116 * t46) * qJD(5)) - t77 * mrSges(4,2) + t76 * mrSges(4,1) + t175 * t34 + t176 * mrSges(6,3) + (m(7) * t103 + t174 + t72) * t33; m(7) * (t1 * t61 + t103 * t35 + t2 * t60 + t43 * t6 + t44 * t5) + (-t113 * t69 + t116 * t68 + m(6) * (-t113 * t31 + t116 * t32)) * qJD(5) + (-mrSges(4,1) * t139 + mrSges(4,2) * t140) * pkin(8) + (-t135 * t90 - t167 * t89) * mrSges(5,3) + (m(6) * t177 - t113 * t59 + t116 * t58) * t106 + t177 * mrSges(6,3) + t175 * t51 + t174 * t50 + Ifges(4,5) * t139 - t1 * t163 - t6 * t164 - t126 * t160 / 0.2e1 + (Ifges(6,5) * t113 + Ifges(7,5) * t100 + Ifges(6,6) * t116 - Ifges(7,6) * t127) * t89 / 0.2e1 - t127 * t7 / 0.2e1 - Ifges(4,6) * t140 - t2 * t157 + (Ifges(6,1) * t113 + t158) * t151 / 0.2e1 - (Ifges(6,2) * t116 + t159) * t154 / 0.2e1 + t5 * t165 + t43 * t41 + t44 * t42 + t60 * t14 + t61 * t15 + t48 * t64 - t52 * t65 / 0.2e1 - t53 * t66 / 0.2e1 + t35 * t72 + t27 * t73 / 0.2e1 + t26 * t74 / 0.2e1 - Ifges(5,6) * t89 + Ifges(5,5) * t90 - t91 * t21 / 0.2e1 - t92 * t20 / 0.2e1 + t100 * t8 / 0.2e1 + t103 * t9 + t108 * t49 + t113 * t39 / 0.2e1 + t116 * t38 / 0.2e1; -t91 * t74 + t100 * t66 - t92 * t73 - t127 * t65 + 0.2e1 * t103 * t64 + (t43 * t61 + t44 * t60) * t172 + 0.2e1 * (-t100 * t44 - t127 * t43 + t60 * t91 - t61 * t92) * mrSges(7,3) + (t106 * t173 + 0.2e1 * mrSges(6,3)) * qJD(5) * (t113 ^ 2 + t112); m(5) * t134 + m(7) * (t100 * t3 - t12 * t92 - t127 * t4 - t13 * t91) + m(6) * (t113 * t30 + t116 * t29); m(5) * t136 + t100 * t15 + t113 * t58 + t116 * t59 - t127 * t14 - t91 * t41 - t92 * t42 + m(7) * (t1 * t100 - t127 * t2 - t5 * t92 - t6 * t91) + m(6) * (t113 * t17 + t116 * t16) + t63; m(7) * (t100 * t43 - t127 * t44 - t60 * t92 - t61 * t91); (-t100 * t91 + t127 * t92) * t172; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t33; m(6) * t50 + m(7) * t35 + t168; t64; 0; 0; mrSges(7,1) * t4 - mrSges(7,2) * t3; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t137; mrSges(7,1) * t44 - mrSges(7,2) * t43 + t160; -t64; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
