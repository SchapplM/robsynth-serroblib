% Calculate time derivative of joint inertia matrix for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:35
% EndTime: 2019-03-08 21:45:39
% DurationCPUTime: 2.22s
% Computational Cost: add. (1437->326), mult. (3623->463), div. (0->0), fcn. (2773->8), ass. (0->145)
t171 = pkin(4) + pkin(8);
t181 = Ifges(7,4) + Ifges(6,5);
t180 = Ifges(7,2) + Ifges(6,3);
t96 = cos(qJ(5));
t144 = qJD(5) * t96;
t93 = sin(qJ(5));
t145 = qJD(5) * t93;
t95 = sin(qJ(2));
t149 = qJD(2) * t95;
t91 = sin(pkin(6));
t139 = t91 * t149;
t98 = cos(qJ(2));
t163 = t91 * t98;
t140 = t96 * t163;
t148 = qJD(2) * t98;
t138 = t91 * t148;
t164 = t91 * t95;
t92 = cos(pkin(6));
t94 = sin(qJ(3));
t97 = cos(qJ(3));
t44 = t164 * t97 + t92 * t94;
t25 = qJD(3) * t44 + t138 * t94;
t141 = t94 * t164;
t43 = -t92 * t97 + t141;
t8 = -t25 * t96 - qJD(5) * t140 + (qJD(5) * t43 + t139) * t93;
t170 = t8 * t96;
t27 = t163 * t93 + t43 * t96;
t28 = t43 * t93 - t140;
t179 = t28 * t144 - t145 * t27 - t170;
t147 = qJD(3) * t94;
t137 = t96 * t147;
t143 = qJD(5) * t97;
t108 = t93 * t143 + t137;
t178 = 2 * qJ(4);
t177 = m(6) + m(7);
t130 = -qJ(4) * t94 - pkin(2);
t99 = -pkin(3) - pkin(9);
t48 = t97 * t99 + t130;
t76 = t171 * t94;
t176 = t96 * t48 + t93 * t76;
t70 = t93 * mrSges(7,1) - t96 * mrSges(7,3);
t71 = t93 * mrSges(6,1) + t96 * mrSges(6,2);
t175 = t70 + t71;
t173 = m(7) * qJ(6) + mrSges(7,3);
t129 = pkin(3) * t147 - qJD(4) * t94;
t35 = (pkin(9) * t94 - qJ(4) * t97) * qJD(3) + t129;
t146 = qJD(3) * t97;
t64 = t171 * t146;
t7 = -qJD(5) * t176 - t35 * t93 + t64 * t96;
t172 = m(5) / 0.2e1;
t9 = qJD(5) * t27 + t139 * t96 + t25 * t93;
t5 = t93 * t9;
t169 = Ifges(6,4) * t93;
t168 = Ifges(6,4) * t96;
t167 = Ifges(7,5) * t93;
t166 = Ifges(7,5) * t96;
t165 = Ifges(7,6) * t94;
t26 = -qJD(3) * t141 + (qJD(3) * t92 + t138) * t97;
t11 = t44 * t26;
t162 = t93 * t97;
t161 = t93 * t99;
t160 = t94 * Ifges(6,6);
t159 = t96 * t97;
t158 = t96 * t99;
t157 = -mrSges(4,1) + mrSges(5,2);
t156 = mrSges(5,3) - mrSges(4,2);
t30 = -mrSges(6,2) * t146 + mrSges(6,3) * t108;
t33 = mrSges(7,2) * t108 + mrSges(7,3) * t146;
t154 = t30 + t33;
t135 = t93 * t147;
t107 = -t143 * t96 + t135;
t31 = mrSges(6,1) * t146 - mrSges(6,3) * t107;
t32 = mrSges(7,2) * t135 + (-mrSges(7,1) * qJD(3) - mrSges(7,2) * t144) * t97;
t153 = t31 - t32;
t122 = Ifges(7,1) * t93 - t166;
t38 = t94 * Ifges(7,4) - t122 * t97;
t123 = Ifges(6,1) * t93 + t168;
t39 = t94 * Ifges(6,5) - t123 * t97;
t152 = t38 + t39;
t59 = mrSges(6,1) * t94 + mrSges(6,3) * t162;
t60 = -mrSges(7,1) * t94 - mrSges(7,2) * t162;
t151 = -t59 + t60;
t61 = -mrSges(6,2) * t94 - mrSges(6,3) * t159;
t62 = -mrSges(7,2) * t159 + mrSges(7,3) * t94;
t150 = t61 + t62;
t77 = t171 * t97;
t142 = qJD(6) * t93;
t134 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t72 = Ifges(7,3) * t93 + t166;
t73 = -Ifges(6,2) * t93 + t168;
t133 = -t72 / 0.2e1 + t73 / 0.2e1;
t74 = Ifges(7,1) * t96 + t167;
t75 = Ifges(6,1) * t96 - t169;
t132 = t74 / 0.2e1 + t75 / 0.2e1;
t131 = m(5) * pkin(8) + mrSges(5,1);
t120 = -Ifges(7,3) * t96 + t167;
t36 = -t120 * t97 + t165;
t121 = Ifges(6,2) * t96 + t169;
t37 = -t121 * t97 + t160;
t127 = -t36 + t37 - t165;
t126 = (m(7) * t99 - mrSges(7,2)) * t93;
t125 = mrSges(6,1) * t96 - mrSges(6,2) * t93;
t124 = mrSges(7,1) * t96 + mrSges(7,3) * t93;
t119 = pkin(5) * t96 + qJ(6) * t93;
t118 = -pkin(5) * t93 + qJ(6) * t96;
t17 = qJ(6) * t94 + t176;
t23 = -t48 * t93 + t76 * t96;
t18 = -pkin(5) * t94 - t23;
t117 = t17 * t96 + t18 * t93;
t116 = t176 * t96 - t23 * t93;
t115 = t25 * t94 + t26 * t97;
t112 = t26 * qJ(4) + t44 * qJD(4);
t111 = t108 * Ifges(6,6) + t181 * t135 + t180 * t146;
t6 = t76 * t144 - t145 * t48 + t96 * t35 + t93 * t64;
t105 = t150 * t96 + t151 * t93;
t104 = t9 * t161 + t179 * t99;
t100 = m(7) * t118 - t175;
t85 = Ifges(7,6) * t144;
t69 = t97 * mrSges(5,2) - t94 * mrSges(5,3);
t66 = -pkin(3) * t97 + t130;
t65 = qJ(4) - t118;
t63 = t171 * t147;
t58 = t123 * qJD(5);
t57 = t122 * qJD(5);
t56 = t121 * qJD(5);
t55 = t120 * qJD(5);
t54 = (mrSges(4,1) * t94 + mrSges(4,2) * t97) * qJD(3);
t53 = (-mrSges(5,2) * t94 - mrSges(5,3) * t97) * qJD(3);
t52 = t125 * qJD(5);
t51 = t124 * qJD(5);
t47 = t125 * t97;
t46 = t124 * t97;
t42 = -qJ(4) * t146 + t129;
t40 = qJD(5) * t119 - qJD(6) * t96 + qJD(4);
t34 = t119 * t97 + t77;
t20 = -mrSges(6,1) * t108 + mrSges(6,2) * t107;
t19 = -mrSges(7,1) * t108 - mrSges(7,3) * t107;
t15 = -t75 * t143 + (t97 * Ifges(6,5) + t123 * t94) * qJD(3);
t14 = -t74 * t143 + (t97 * Ifges(7,4) + t122 * t94) * qJD(3);
t13 = -t73 * t143 + (t97 * Ifges(6,6) + t121 * t94) * qJD(3);
t12 = -t72 * t143 + (t97 * Ifges(7,6) + t120 * t94) * qJD(3);
t10 = (qJD(5) * t118 + t142) * t97 + (-t119 - t171) * t147;
t3 = -pkin(5) * t146 - t7;
t2 = qJ(6) * t146 + qJD(6) * t94 + t6;
t1 = [0.2e1 * t177 * (-t27 * t8 + t28 * t9 + t11) + 0.2e1 * (m(5) + m(4)) * (-t148 * t91 ^ 2 * t95 + t25 * t43 + t11); t150 * t9 + t151 * t8 + (t19 + t20) * t44 + t154 * t28 + t153 * t27 + (t46 + t47) * t26 + m(6) * (t176 * t9 - t23 * t8 + t26 * t77 + t27 * t7 + t28 * t6 - t44 * t63) + m(7) * (t10 * t44 + t17 * t9 + t18 * t8 + t2 * t28 + t26 * t34 - t27 * t3) + 0.2e1 * (m(4) / 0.2e1 + t172) * (t146 * t43 - t147 * t44 + t115) * pkin(8) + (mrSges(4,3) + mrSges(5,1)) * ((t43 * t97 - t44 * t94) * qJD(3) + t115) + ((-t53 - t54) * t98 + (-t98 * mrSges(3,2) + (-t97 * mrSges(4,1) + t94 * mrSges(4,2) - mrSges(3,1) + t69) * t95) * qJD(2) - m(4) * pkin(2) * t149 + 0.2e1 * (t66 * t149 - t42 * t98) * t172) * t91; -0.2e1 * pkin(2) * t54 + 0.2e1 * t10 * t46 + 0.2e1 * t17 * t33 + 0.2e1 * t18 * t32 + 0.2e1 * t34 * t19 + 0.2e1 * t2 * t62 + 0.2e1 * t77 * t20 + 0.2e1 * t23 * t31 + 0.2e1 * t176 * t30 + 0.2e1 * t3 * t60 - 0.2e1 * t63 * t47 + 0.2e1 * t66 * t53 + 0.2e1 * t7 * t59 + 0.2e1 * t6 * t61 + 0.2e1 * (m(5) * t66 + t69) * t42 + 0.2e1 * m(6) * (t176 * t6 + t23 * t7 - t63 * t77) + 0.2e1 * m(7) * (t10 * t34 + t17 * t2 + t18 * t3) + ((0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t94 + t127 * t96 + t152 * t93) * qJD(3) + t111) * t94 + ((t12 - t13) * t96 + (-t14 - t15) * t93 + (t127 * t93 + (-t181 * t94 - t152) * t96) * qJD(5) + ((0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + (-Ifges(6,6) + Ifges(7,6)) * t96 - t181 * t93) * t97 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) - (2 * Ifges(5,3)) + (2 * Ifges(5,2)) + t180) * t94) * qJD(3)) * t97; (t51 + t52) * t44 + t157 * t25 + (t156 + t175) * t26 + m(6) * (t104 + t112) + m(7) * (t65 * t26 + t40 * t44 + t104) + m(5) * (-pkin(3) * t25 + t112) + (mrSges(6,3) + mrSges(7,2)) * (t170 - t5 + (t27 * t93 - t28 * t96) * qJD(5)); t77 * t52 + t65 * t19 + t10 * t70 - t63 * t71 + t40 * t46 + qJD(4) * t47 + t34 * t51 + qJ(4) * t20 + t94 * t85 / 0.2e1 + (t14 / 0.2e1 + t15 / 0.2e1 + t3 * mrSges(7,2) - t7 * mrSges(6,3) + t153 * t99) * t96 + (t12 / 0.2e1 - t13 / 0.2e1 - t2 * mrSges(7,2) - t6 * mrSges(6,3) + t154 * t99) * t93 + m(7) * (t65 * t10 - t158 * t3 + t161 * t2 + t40 * t34) + m(6) * (-qJ(4) * t63 + qJD(4) * t77 + t158 * t7 + t161 * t6) + ((-t55 / 0.2e1 + t56 / 0.2e1) * t96 + (t57 / 0.2e1 + t58 / 0.2e1) * t93 + t131 * qJD(4)) * t97 + ((-t160 / 0.2e1 - t176 * mrSges(6,3) - t17 * mrSges(7,2) + t36 / 0.2e1 - t37 / 0.2e1 - t132 * t97) * t96 + (t23 * mrSges(6,3) - t18 * mrSges(7,2) - t38 / 0.2e1 - t39 / 0.2e1 + t133 * t97 + t134 * t94) * t93 + (m(6) * t116 + m(7) * t117 + t105) * t99) * qJD(5) + ((-pkin(3) * mrSges(5,1) - Ifges(5,4) + Ifges(4,5) - t134 * t96 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t93 + (-m(5) * pkin(3) + t157) * pkin(8)) * t97 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + t133 * t96 + t132 * t93 + (-m(5) * qJ(4) - t156) * pkin(8)) * t94) * qJD(3); t52 * t178 + 0.2e1 * t51 * t65 + (-t57 - t58) * t96 + (-t55 + t56) * t93 + ((t72 - t73) * t96 + (-t74 - t75) * t93) * qJD(5) + 0.2e1 * (m(7) * t65 + t70) * t40 + (0.2e1 * mrSges(5,3) + 0.2e1 * t71 + (m(5) + m(6)) * t178) * qJD(4); m(5) * t25 + t177 * (t5 + t179); t153 * t96 + t154 * t93 + t131 * t146 + t105 * qJD(5) + m(7) * (qJD(5) * t117 + t2 * t93 - t3 * t96) + m(6) * (qJD(5) * t116 + t6 * t93 + t7 * t96); 0; 0; m(7) * qJD(6) * t28 + (-mrSges(6,2) + t173) * t9 + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t8; -Ifges(7,6) * t137 - pkin(5) * t32 + m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t17) + qJD(6) * t62 + qJ(6) * t33 + t2 * mrSges(7,3) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t7 * mrSges(6,1) + (-Ifges(7,6) * t93 - t181 * t96) * t143 + t111; t85 + qJD(6) * t126 + ((-qJ(6) * mrSges(7,2) - Ifges(6,6)) * t96 + (mrSges(7,2) * pkin(5) - t181) * t93 + t100 * t99) * qJD(5); m(7) * t142 + qJD(5) * t100; 0.2e1 * t173 * qJD(6); m(7) * t8; m(7) * t3 + t32; qJD(5) * t126; m(7) * t145; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
