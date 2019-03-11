% Calculate time derivative of joint inertia matrix for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:09
% EndTime: 2019-03-08 18:49:13
% DurationCPUTime: 1.68s
% Computational Cost: add. (1602->262), mult. (4772->399), div. (0->0), fcn. (4631->12), ass. (0->132)
t175 = m(5) + m(6);
t80 = cos(qJ(4));
t129 = qJD(6) * t80;
t77 = sin(qJ(4));
t131 = qJD(4) * t77;
t76 = sin(qJ(6));
t79 = cos(qJ(6));
t94 = t76 * t129 + t79 * t131;
t134 = cos(pkin(7));
t74 = cos(pkin(12));
t116 = t74 * t134;
t72 = sin(pkin(7));
t81 = cos(qJ(3));
t147 = t72 * t81;
t71 = sin(pkin(12));
t73 = sin(pkin(6));
t75 = cos(pkin(6));
t78 = sin(qJ(3));
t174 = t75 * t147 + (t81 * t116 - t71 * t78) * t73;
t135 = qJ(5) * t77;
t104 = -pkin(4) * t80 - t135;
t53 = -pkin(3) + t104;
t173 = m(6) * t53;
t172 = mrSges(5,3) + mrSges(6,1);
t57 = mrSges(7,1) * t76 + mrSges(7,2) * t79;
t171 = -mrSges(6,3) - t57;
t55 = t80 * mrSges(6,2) - t77 * mrSges(6,3);
t168 = -t80 * mrSges(5,1) + t77 * mrSges(5,2) + t55;
t148 = t72 * t78;
t39 = -t80 * t134 + t77 * t148;
t159 = m(6) / 0.2e1;
t165 = 0.2e1 * pkin(9) * (t159 + m(5) / 0.2e1);
t133 = qJD(3) * t78;
t124 = t72 * t133;
t132 = qJD(3) * t81;
t123 = t72 * t132;
t40 = t134 * t77 + t148 * t80;
t29 = qJD(4) * t40 + t123 * t77;
t30 = t147 * t76 + t39 * t79;
t12 = qJD(6) * t30 + t124 * t79 + t29 * t76;
t96 = t147 * t79 - t39 * t76;
t13 = qJD(6) * t96 - t124 * t76 + t29 * t79;
t164 = qJD(6) * (t30 * t76 + t79 * t96) - t12 * t76 - t13 * t79;
t85 = t75 * t148 + (t116 * t78 + t71 * t81) * t73;
t23 = t85 * qJD(3);
t92 = -t73 * t74 * t72 + t134 * t75;
t16 = t77 * t92 + t80 * t85;
t22 = t174 * qJD(3);
t5 = qJD(4) * t16 + t22 * t77;
t83 = t85 * t77;
t15 = -t80 * t92 + t83;
t8 = t15 * t76 - t174 * t79;
t1 = -qJD(6) * t8 - t23 * t76 + t5 * t79;
t7 = t15 * t79 + t174 * t76;
t2 = qJD(6) * t7 + t23 * t79 + t5 * t76;
t163 = -t1 * t79 - t2 * t76 + (t7 * t76 - t79 * t8) * qJD(6);
t82 = -pkin(4) - pkin(10);
t42 = t80 * t82 - pkin(3) - t135;
t154 = pkin(5) + pkin(9);
t61 = t154 * t77;
t26 = -t42 * t76 + t61 * t79;
t27 = t42 * t79 + t61 * t76;
t101 = t26 * t76 - t27 * t79;
t152 = mrSges(7,3) * t80;
t49 = mrSges(7,1) * t77 + t152 * t76;
t50 = -mrSges(7,2) * t77 - t152 * t79;
t162 = -m(7) * t101 - t76 * t49 + t79 * t50;
t158 = m(5) * pkin(3);
t157 = -t76 / 0.2e1;
t156 = -t79 / 0.2e1;
t155 = t79 / 0.2e1;
t28 = t39 * qJD(4) - t80 * t123;
t6 = -qJD(4) * t83 + (qJD(4) * t92 + t22) * t80;
t153 = -t28 * t16 + t40 * t6;
t3 = t16 * t6;
t151 = Ifges(7,4) * t76;
t150 = Ifges(7,4) * t79;
t14 = t174 * t23;
t17 = t40 * t28;
t107 = Ifges(7,1) * t76 + t150;
t142 = t77 * Ifges(7,5);
t37 = -t107 * t80 + t142;
t146 = t76 * t37;
t59 = Ifges(7,1) * t79 - t151;
t144 = t76 * t59;
t143 = t76 * t82;
t106 = Ifges(7,2) * t79 + t151;
t36 = t77 * Ifges(7,6) - t106 * t80;
t141 = t79 * t36;
t58 = -Ifges(7,2) * t76 + t150;
t139 = t79 * t58;
t138 = t79 * t82;
t137 = -mrSges(5,1) + mrSges(6,2);
t44 = (-mrSges(6,2) * t77 - mrSges(6,3) * t80) * qJD(4);
t45 = (mrSges(5,1) * t77 + mrSges(5,2) * t80) * qJD(4);
t136 = t44 + t45;
t130 = qJD(4) * t80;
t126 = -mrSges(4,1) + t168;
t125 = mrSges(5,2) + t171;
t62 = t154 * t80;
t121 = t79 * t129;
t120 = t76 * t131;
t118 = m(6) * pkin(9) + mrSges(6,1);
t114 = Ifges(7,5) * t120 + t94 * Ifges(7,6) + Ifges(7,3) * t130;
t113 = pkin(4) * t131 - qJD(5) * t77;
t20 = t174 * t124;
t110 = t5 * t77 + t6 * t80;
t108 = mrSges(7,1) * t79 - mrSges(7,2) * t76;
t105 = -Ifges(7,5) * t76 - Ifges(7,6) * t79;
t34 = (pkin(10) * t77 - qJ(5) * t80) * qJD(4) + t113;
t52 = qJD(4) * t62;
t10 = qJD(6) * t26 + t34 * t79 + t52 * t76;
t11 = -qJD(6) * t27 - t34 * t76 + t52 * t79;
t103 = -t10 * t76 - t11 * t79;
t100 = -t28 * t80 + t29 * t77;
t98 = t6 * qJ(5) + t16 * qJD(5);
t97 = -qJ(5) * t28 + qJD(5) * t40;
t93 = t120 - t121;
t86 = m(7) * t163;
t84 = m(7) * t164;
t51 = t154 * t131;
t47 = t107 * qJD(6);
t46 = t106 * qJD(6);
t43 = t108 * qJD(6);
t41 = t108 * t80;
t38 = -qJ(5) * t130 + t113;
t33 = mrSges(7,1) * t130 - mrSges(7,3) * t93;
t32 = -mrSges(7,2) * t130 + mrSges(7,3) * t94;
t24 = -mrSges(7,1) * t94 + mrSges(7,2) * t93;
t19 = -t59 * t129 + (t80 * Ifges(7,5) + t107 * t77) * qJD(4);
t18 = -t58 * t129 + (t80 * Ifges(7,6) + t106 * t77) * qJD(4);
t4 = [0.2e1 * m(7) * (t1 * t7 + t2 * t8 + t3) + 0.2e1 * m(4) * (t22 * t85 - t14) + 0.2e1 * t175 * (t15 * t5 - t14 + t3); m(7) * (t1 * t30 + t12 * t8 + t13 * t7 - t2 * t96 + t153) + m(4) * (t148 * t22 - t20) + t175 * (-t147 * t23 + t15 * t29 + t39 * t5 + t153 - t20); 0.2e1 * m(7) * (-t12 * t96 + t13 * t30 - t17) + 0.2e1 * t175 * (-t132 * t72 ^ 2 * t78 + t29 * t39 - t17); -t22 * mrSges(4,2) + t1 * t49 + t16 * t24 + t2 * t50 + t8 * t32 + t7 * t33 + t6 * t41 + m(7) * (t1 * t26 + t10 * t8 + t11 * t7 - t16 * t51 + t2 * t27 + t6 * t62) + t172 * ((t15 * t80 - t16 * t77) * qJD(4) + t110) + (t130 * t15 - t131 * t16 + t110) * t165 - (m(6) * t38 + t136) * t174 + (t126 - t158 + t173) * t23; t12 * t50 + t13 * t49 + t40 * t24 - t28 * t41 + t30 * t33 - t96 * t32 + m(7) * (-t10 * t96 + t11 * t30 + t12 * t27 + t13 * t26 - t28 * t62 - t40 * t51) + t172 * ((t39 * t80 - t40 * t77) * qJD(4) + t100) + (t130 * t39 - t131 * t40 + t100) * t165 + (-t136 * t81 + (-mrSges(4,2) * t81 + t126 * t78) * qJD(3) - t133 * t158 + 0.2e1 * (t133 * t53 - t38 * t81) * t159) * t72; 0.2e1 * m(7) * (t10 * t27 + t11 * t26 - t51 * t62) + 0.2e1 * t11 * t49 + 0.2e1 * t10 * t50 - 0.2e1 * t51 * t41 + 0.2e1 * t53 * t44 + 0.2e1 * t62 * t24 - 0.2e1 * pkin(3) * t45 + 0.2e1 * t27 * t32 + 0.2e1 * t26 * t33 + 0.2e1 * (t55 + t173) * t38 + ((t141 + t146 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t77) * qJD(4) + t114) * t77 + (-t79 * t18 - t76 * t19 + (t76 * t36 + (-t37 - t142) * t79) * qJD(6) + ((0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) + t105) * t80 + ((2 * Ifges(5,1)) - (2 * Ifges(5,2)) + (2 * Ifges(6,2)) - (2 * Ifges(6,3)) + Ifges(7,3)) * t77) * qJD(4)) * t80; t16 * t43 + t137 * t5 - t125 * t6 + m(6) * (-pkin(4) * t5 + t98) + m(7) * t98 - t82 * t86 + t163 * mrSges(7,3); t40 * t43 + t137 * t29 + t125 * t28 + m(7) * t97 + m(6) * (-pkin(4) * t29 + t97) - t82 * t84 + t164 * mrSges(7,3); t19 * t155 + t18 * t157 + m(7) * (-qJ(5) * t51 + qJD(5) * t62 + t10 * t143 + t11 * t138) + t33 * t138 + t32 * t143 - t51 * t57 + t62 * t43 + qJD(5) * t41 + qJ(5) * t24 + t103 * mrSges(7,3) + (t118 * qJD(5) - t46 * t156 - t47 * t157) * t80 + (t77 * t105 / 0.2e1 - t141 / 0.2e1 - t146 / 0.2e1 + (t59 * t156 + t76 * t58 / 0.2e1) * t80 + t101 * mrSges(7,3) + t162 * t82) * qJD(6) + ((-pkin(4) * mrSges(6,1) + Ifges(7,5) * t155 + Ifges(7,6) * t157 - Ifges(6,4) + Ifges(5,5)) * t80 + (t139 / 0.2e1 + t144 / 0.2e1 - qJ(5) * mrSges(6,1) - Ifges(5,6) + Ifges(6,5)) * t77 + (m(6) * t104 + t168) * pkin(9)) * qJD(4); 0.2e1 * qJ(5) * t43 + t46 * t76 - t47 * t79 + (-t139 - t144) * qJD(6) + 0.2e1 * ((m(6) + m(7)) * qJ(5) - t171) * qJD(5); m(6) * t5 - t86; m(6) * t29 - t84; -m(7) * t103 + t162 * qJD(6) + t118 * t130 + t76 * t32 + t79 * t33; 0; 0; mrSges(7,1) * t1 - mrSges(7,2) * t2; mrSges(7,1) * t13 - mrSges(7,2) * t12; mrSges(7,1) * t11 - mrSges(7,2) * t10 - Ifges(7,5) * t121 + t114; ((-mrSges(7,2) * t82 - Ifges(7,6)) * t79 + (-mrSges(7,1) * t82 - Ifges(7,5)) * t76) * qJD(6); -t57 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t4(1) t4(2) t4(4) t4(7) t4(11) t4(16); t4(2) t4(3) t4(5) t4(8) t4(12) t4(17); t4(4) t4(5) t4(6) t4(9) t4(13) t4(18); t4(7) t4(8) t4(9) t4(10) t4(14) t4(19); t4(11) t4(12) t4(13) t4(14) t4(15) t4(20); t4(16) t4(17) t4(18) t4(19) t4(20) t4(21);];
Mq  = res;
