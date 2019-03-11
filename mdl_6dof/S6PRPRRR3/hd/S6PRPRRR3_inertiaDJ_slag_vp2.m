% Calculate time derivative of joint inertia matrix for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:59
% EndTime: 2019-03-08 20:32:03
% DurationCPUTime: 2.34s
% Computational Cost: add. (5008->302), mult. (11894->470), div. (0->0), fcn. (12685->12), ass. (0->140)
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t119 = cos(qJ(5));
t189 = (t114 ^ 2 + t118 ^ 2) * t119;
t97 = -mrSges(7,1) * t118 + mrSges(7,2) * t114;
t191 = t97 - mrSges(6,1);
t148 = qJD(6) * t118;
t115 = sin(qJ(5));
t110 = sin(pkin(12));
t116 = sin(qJ(4));
t112 = cos(pkin(12));
t120 = cos(qJ(4));
t155 = t112 * t120;
t90 = -t110 * t116 + t155;
t91 = t110 * t120 + t112 * t116;
t71 = t115 * t91 - t119 * t90;
t81 = t90 * qJD(4);
t82 = t91 * qJD(4);
t48 = -qJD(5) * t71 - t115 * t82 + t119 * t81;
t72 = t115 * t90 + t119 * t91;
t130 = t114 * t48 + t72 * t148;
t188 = 2 * m(7);
t187 = -2 * mrSges(6,3);
t176 = pkin(8) + qJ(3);
t95 = t176 * t110;
t96 = t176 * t112;
t74 = -t116 * t95 + t120 * t96;
t64 = -qJD(3) * t91 - qJD(4) * t74;
t125 = -pkin(9) * t81 + t64;
t88 = t120 * t95;
t73 = -t116 * t96 - t88;
t131 = -pkin(9) * t91 + t73;
t67 = pkin(9) * t90 + t74;
t35 = t115 * t131 + t119 * t67;
t63 = qJD(3) * t155 - qJD(4) * t88 + (-qJD(3) * t110 - qJD(4) * t96) * t116;
t57 = -pkin(9) * t82 + t63;
t19 = qJD(5) * t35 + t115 * t57 - t119 * t125;
t186 = 0.2e1 * t19;
t34 = t115 * t67 - t119 * t131;
t185 = 0.2e1 * t34;
t184 = m(6) / 0.2e1;
t182 = pkin(4) * t82;
t111 = sin(pkin(6));
t117 = sin(qJ(2));
t151 = qJD(2) * t117;
t144 = t111 * t151;
t113 = cos(pkin(6));
t157 = t111 * t117;
t78 = -t110 * t157 + t112 * t113;
t79 = t110 * t113 + t112 * t157;
t68 = -t116 * t79 + t120 * t78;
t69 = t116 * t78 + t120 * t79;
t39 = t115 * t69 - t119 * t68;
t121 = cos(qJ(2));
t150 = qJD(2) * t121;
t143 = t111 * t150;
t54 = qJD(4) * t68 + t143 * t90;
t55 = -qJD(4) * t69 - t143 * t91;
t15 = -qJD(5) * t39 + t115 * t55 + t119 * t54;
t156 = t111 * t121;
t40 = t115 * t68 + t119 * t69;
t31 = -t114 * t40 - t118 * t156;
t5 = qJD(6) * t31 + t114 * t144 + t118 * t15;
t181 = t118 * t5;
t180 = t19 * t34;
t18 = -qJD(5) * t34 + t115 * t125 + t119 * t57;
t101 = -pkin(3) * t112 - pkin(2);
t75 = -pkin(4) * t90 + t101;
t36 = pkin(5) * t71 - pkin(10) * t72 + t75;
t23 = t114 * t36 + t118 * t35;
t49 = qJD(5) * t72 + t115 * t81 + t119 * t82;
t24 = pkin(5) * t49 - pkin(10) * t48 + t182;
t3 = -qJD(6) * t23 - t114 * t18 + t118 * t24;
t179 = t3 * t114;
t16 = qJD(5) * t40 + t115 * t54 - t119 * t55;
t178 = t39 * t16;
t128 = t114 * t156 - t118 * t40;
t6 = qJD(6) * t128 - t114 * t15 + t118 * t144;
t177 = t6 * t114;
t163 = t118 * t48;
t175 = Ifges(7,5) * t163 + Ifges(7,3) * t49;
t174 = mrSges(7,3) * t118;
t173 = Ifges(7,4) * t114;
t172 = Ifges(7,4) * t118;
t171 = Ifges(7,6) * t114;
t170 = pkin(4) * qJD(5);
t168 = t114 * t72;
t167 = t115 * mrSges(6,1);
t166 = t115 * t34;
t165 = t115 * t39;
t164 = t115 * t97;
t162 = t118 * t72;
t161 = t119 * mrSges(6,2);
t160 = t111 ^ 2 * t117;
t154 = t114 * t119;
t153 = t118 * t119;
t152 = t110 ^ 2 + t112 ^ 2;
t149 = qJD(6) * t114;
t147 = t72 * t149;
t129 = t147 - t163;
t12 = mrSges(7,1) * t130 - mrSges(7,2) * t129;
t145 = m(7) * t19 + t12;
t25 = t49 * mrSges(6,1) + t48 * mrSges(6,2);
t141 = -(2 * Ifges(6,4)) - t171;
t139 = t150 * t160;
t138 = mrSges(7,3) * t189;
t137 = t34 * t16 + t19 * t39;
t136 = mrSges(7,1) * t114 + mrSges(7,2) * t118;
t135 = Ifges(7,1) * t118 - t173;
t134 = -Ifges(7,2) * t114 + t172;
t133 = Ifges(7,5) * t114 + Ifges(7,6) * t118;
t132 = -t110 * t78 + t112 * t79;
t22 = -t114 * t35 + t118 * t36;
t93 = t134 * qJD(6);
t94 = t135 * qJD(6);
t98 = Ifges(7,2) * t118 + t173;
t99 = Ifges(7,1) * t114 + t172;
t127 = t114 * t94 + t118 * t93 + t99 * t148 - t149 * t98;
t126 = -t177 + (t114 * t128 - t118 * t31) * qJD(6);
t92 = t136 * qJD(6);
t124 = -t15 * mrSges(6,2) + mrSges(7,3) * t126 + t191 * t16 + t5 * t174 + t39 * t92;
t2 = qJD(6) * t22 + t114 * t24 + t118 * t18;
t20 = mrSges(7,1) * t49 + mrSges(7,3) * t129;
t21 = -mrSges(7,2) * t49 - mrSges(7,3) * t130;
t50 = -mrSges(7,2) * t71 - mrSges(7,3) * t168;
t51 = mrSges(7,1) * t71 - mrSges(7,3) * t162;
t123 = -t51 * t148 - t50 * t149 + m(7) * (t118 * t2 - t148 * t22 - t149 * t23 - t179) + t118 * t21 - t114 * t20;
t10 = -Ifges(7,1) * t129 - Ifges(7,4) * t130 + t49 * Ifges(7,5);
t104 = Ifges(7,5) * t148;
t28 = t71 * Ifges(7,6) + t134 * t72;
t29 = t71 * Ifges(7,5) + t135 * t72;
t9 = -Ifges(7,4) * t129 - Ifges(7,2) * t130 + t49 * Ifges(7,6);
t122 = t2 * t174 + t29 * t148 / 0.2e1 + t34 * t92 + t99 * t163 / 0.2e1 + Ifges(6,5) * t48 - t93 * t168 / 0.2e1 + t94 * t162 / 0.2e1 + t71 * (-Ifges(7,6) * t149 + t104) / 0.2e1 + t114 * t10 / 0.2e1 + t118 * t9 / 0.2e1 + (-t179 + (-t114 * t23 - t118 * t22) * qJD(6)) * mrSges(7,3) - t18 * mrSges(6,2) + (t133 / 0.2e1 - Ifges(6,6)) * t49 + t191 * t19 - t130 * t98 / 0.2e1 - (t72 * t99 + t28) * t149 / 0.2e1;
t103 = -pkin(4) * t119 - pkin(5);
t102 = pkin(4) * t115 + pkin(10);
t77 = t81 * mrSges(5,2);
t70 = t82 * mrSges(5,1) + t77;
t53 = mrSges(6,1) * t71 + mrSges(6,2) * t72;
t43 = t136 * t72;
t1 = [0.2e1 * m(7) * (-t128 * t5 + t31 * t6 + t178) + 0.2e1 * m(6) * (t40 * t15 - t139 + t178) + 0.2e1 * m(5) * (t69 * t54 + t68 * t55 - t139) + 0.2e1 * m(4) * (t111 * t132 - t160) * t150; t39 * t12 + t16 * t43 + t31 * t20 - t128 * t21 + t5 * t50 + t6 * t51 + (-t15 * t71 + t16 * t72 + t39 * t48 - t40 * t49) * mrSges(6,3) + (t54 * t90 - t55 * t91 - t68 * t81 - t69 * t82) * mrSges(5,3) + ((-t25 - t70) * t121 + ((mrSges(4,3) * t152 - mrSges(3,2)) * t121 + (-mrSges(4,1) * t112 - mrSges(5,1) * t90 + mrSges(4,2) * t110 + mrSges(5,2) * t91 - mrSges(3,1) + t53) * t117) * qJD(2)) * t111 + m(6) * (t35 * t15 + t18 * t40 + (-t121 * t182 + t151 * t75) * t111 + t137) + m(4) * (t132 * qJD(3) + (qJ(3) * t121 * t152 - pkin(2) * t117) * t111 * qJD(2)) + m(5) * (t101 * t144 + t54 * t74 + t55 * t73 + t63 * t69 + t64 * t68) + m(7) * (-t128 * t2 + t22 * t6 + t23 * t5 + t3 * t31 + t137); t35 * t49 * t187 + 0.2e1 * t91 * t81 * Ifges(5,1) + 0.2e1 * t101 * t70 + t12 * t185 + t43 * t186 + 0.2e1 * t2 * t50 + 0.2e1 * t22 * t20 + 0.2e1 * t23 * t21 + 0.2e1 * t75 * t25 + 0.2e1 * t3 * t51 + (mrSges(6,3) * t185 - t114 * t28 + t118 * t29) * t48 + 0.2e1 * m(5) * (t63 * t74 + t64 * t73) + 0.2e1 * m(6) * (t18 * t35 + t182 * t75 + t180) + (t2 * t23 + t22 * t3 + t180) * t188 + (t18 * t187 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t49 + t141 * t48 + t175) * t71 + (mrSges(6,3) * t186 + 0.2e1 * Ifges(6,1) * t48 + t118 * t10 - t114 * t9 + (Ifges(7,5) * t118 + t141) * t49 + (-t114 * t29 - t118 * t28 - t133 * t71) * qJD(6)) * t72 - 0.2e1 * (Ifges(5,2) * t90 - pkin(4) * t53) * t82 + 0.2e1 * (t90 * t81 - t91 * t82) * Ifges(5,4) + 0.2e1 * (t63 * t90 - t64 * t91 - t73 * t81 - t74 * t82) * mrSges(5,3) + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * t152 * qJD(3); m(7) * (t114 * t5 + t118 * t6 + (-t114 * t31 - t118 * t128) * qJD(6)) + (m(4) + m(5) + m(6)) * t144; -t51 * t149 + t118 * t20 + m(7) * (t114 * t2 + t118 * t3 + (-t114 * t22 + t118 * t23) * qJD(6)) + t50 * t148 + t114 * t21 + t77 - (-m(6) * pkin(4) - mrSges(5,1)) * t82 + t25; 0; m(7) * t103 * t16 - t54 * mrSges(5,2) + t55 * mrSges(5,1) + 0.2e1 * ((t115 * t15 - t119 * t16) * t184 + (m(7) * (-t128 * t153 - t31 * t154 + t165) / 0.2e1 + (t119 * t40 + t165) * t184) * qJD(5)) * pkin(4) + t124 + m(7) * (t128 * t149 - t148 * t31 - t177 + t181) * t102; t122 - Ifges(5,6) * t82 + Ifges(5,5) * t81 - t63 * mrSges(5,2) + t64 * mrSges(5,1) + (m(6) * (t115 * t18 - t119 * t19) + (-t115 * t49 - t119 * t48) * mrSges(6,3) + ((t72 * mrSges(6,3) + t43) * t115 + (-t71 * mrSges(6,3) - t114 * t51 + t118 * t50) * t119 + m(7) * (t153 * t23 - t154 * t22 + t166) + m(6) * (t119 * t35 + t166)) * qJD(5)) * pkin(4) + t145 * t103 + t123 * t102; 0; 0.2e1 * t103 * t92 + (0.2e1 * t164 + (t189 * t102 + t103 * t115) * t188 - 0.2e1 * t161 - 0.2e1 * t167 + 0.2e1 * t138) * t170 + t127; m(7) * (-pkin(5) * t16 + (t126 + t181) * pkin(10)) + t124; -pkin(5) * t145 + pkin(10) * t123 + t122; 0; (-pkin(5) + t103) * t92 + (m(7) * (-pkin(5) * t115 + t189 * pkin(10)) + t164 - t161 - t167 + t138) * t170 + t127; -0.2e1 * pkin(5) * t92 + t127; mrSges(7,1) * t6 - mrSges(7,2) * t5; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t147 - Ifges(7,6) * t130 + t175; -t92; t104 - t136 * t119 * t170 + (t102 * t97 - t171) * qJD(6); t104 + (pkin(10) * t97 - t171) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
