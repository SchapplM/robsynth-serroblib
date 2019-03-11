% Calculate time derivative of joint inertia matrix for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:35
% EndTime: 2019-03-08 21:12:40
% DurationCPUTime: 2.35s
% Computational Cost: add. (1844->374), mult. (4902->561), div. (0->0), fcn. (4209->10), ass. (0->158)
t115 = cos(pkin(11));
t112 = t115 ^ 2;
t113 = sin(pkin(11));
t187 = qJD(4) * (t113 ^ 2 + t112);
t120 = cos(qJ(3));
t139 = qJD(3) * t120;
t131 = t115 * t139;
t132 = t113 * t139;
t68 = mrSges(5,1) * t132 + mrSges(5,2) * t131;
t116 = sin(qJ(6));
t119 = cos(qJ(6));
t124 = t113 * t116 + t115 * t119;
t117 = sin(qJ(3));
t85 = t113 * t119 - t115 * t116;
t65 = t85 * t117;
t28 = qJD(6) * t65 + t124 * t139;
t72 = t124 * qJD(6);
t29 = -t117 * t72 + t139 * t85;
t7 = -t29 * mrSges(7,1) + t28 * mrSges(7,2);
t186 = t68 - t7;
t185 = mrSges(5,3) + mrSges(6,2);
t184 = Ifges(5,1) + Ifges(6,1);
t176 = pkin(4) + pkin(5);
t183 = t176 * t113;
t140 = qJD(3) * t117;
t182 = qJ(5) * t140 - qJD(5) * t120;
t181 = 2 * m(5);
t180 = 2 * m(6);
t179 = 2 * m(7);
t178 = 2 * pkin(8);
t177 = m(5) / 0.2e1;
t145 = t115 * t117;
t151 = -qJ(5) * t131 - qJD(5) * t145;
t173 = pkin(4) * t113;
t40 = (pkin(8) + t173) * t139 + t151;
t175 = m(6) * t40;
t114 = sin(pkin(6));
t121 = cos(qJ(2));
t141 = qJD(2) * t121;
t133 = t114 * t141;
t118 = sin(qJ(2));
t147 = t114 * t118;
t150 = cos(pkin(6));
t75 = t117 * t150 + t120 * t147;
t48 = qJD(3) * t75 + t117 * t133;
t74 = t117 * t147 - t120 * t150;
t172 = t74 * t48;
t171 = -pkin(9) + qJ(4);
t134 = qJD(2) * t147;
t49 = -qJD(3) * t74 + t120 * t133;
t25 = t113 * t134 + t115 * t49;
t155 = t115 * t25;
t146 = t114 * t121;
t47 = -t113 * t146 + t75 * t115;
t170 = t115 * qJD(4) * t47 + qJ(4) * t155;
t169 = Ifges(7,5) * t28 + Ifges(7,6) * t29;
t41 = mrSges(7,1) * t124 + mrSges(7,2) * t85;
t95 = -mrSges(6,1) * t115 - mrSges(6,3) * t113;
t168 = t41 - t95;
t73 = t85 * qJD(6);
t167 = -Ifges(7,5) * t72 - Ifges(7,6) * t73;
t66 = t124 * t117;
t30 = -mrSges(7,1) * t65 + mrSges(7,2) * t66;
t160 = mrSges(6,3) * t115;
t77 = (mrSges(6,1) * t113 - t160) * t117;
t166 = -t77 + t30;
t148 = t113 * t120;
t79 = (-mrSges(5,2) * t117 - mrSges(5,3) * t148) * qJD(3);
t82 = (-mrSges(6,2) * t148 + mrSges(6,3) * t117) * qJD(3);
t165 = t79 + t82;
t144 = t115 * t120;
t80 = (mrSges(5,1) * t117 - mrSges(5,3) * t144) * qJD(3);
t81 = -mrSges(6,1) * t140 + mrSges(6,2) * t131;
t164 = -t80 + t81;
t149 = t113 * t117;
t86 = mrSges(5,2) * t120 - mrSges(5,3) * t149;
t89 = -mrSges(6,2) * t149 - mrSges(6,3) * t120;
t163 = t86 + t89;
t87 = -mrSges(5,1) * t120 - mrSges(5,3) * t145;
t88 = mrSges(6,1) * t120 + mrSges(6,2) * t145;
t162 = -t87 + t88;
t161 = -mrSges(5,1) * t115 + mrSges(5,2) * t113 - mrSges(4,1);
t159 = Ifges(5,4) * t113;
t158 = Ifges(5,4) * t115;
t157 = Ifges(6,5) * t113;
t156 = Ifges(6,5) * t115;
t71 = -qJD(4) * t117 + (pkin(3) * t117 - qJ(4) * t120) * qJD(3);
t154 = t115 * t71;
t153 = t48 * t117;
t152 = t49 * t120;
t93 = -pkin(3) * t120 - qJ(4) * t117 - pkin(2);
t62 = pkin(8) * t144 + t113 * t93;
t143 = qJ(4) * t187;
t138 = qJD(5) * t113;
t136 = pkin(8) * t140;
t135 = -pkin(8) * t113 - pkin(4);
t130 = qJ(5) * t113 + pkin(3);
t129 = qJ(5) * t115 - pkin(8);
t104 = pkin(8) * t148;
t61 = t115 * t93 - t104;
t53 = -qJ(5) * t120 + t62;
t34 = mrSges(7,1) * t73 - mrSges(7,2) * t72;
t110 = t120 * pkin(4);
t33 = pkin(5) * t120 + t104 + t110 + (-pkin(9) * t117 - t93) * t115;
t39 = pkin(9) * t149 + t53;
t8 = -t116 * t39 + t119 * t33;
t9 = t116 * t33 + t119 * t39;
t46 = t75 * t113 + t115 * t146;
t11 = -t116 * t47 + t119 * t46;
t12 = t116 * t46 + t119 * t47;
t94 = t171 * t113;
t97 = t171 * t115;
t50 = -t116 * t97 + t119 * t94;
t51 = t116 * t94 + t119 * t97;
t24 = t113 * t49 - t115 * t134;
t125 = qJ(4) * t24 + qJD(4) * t46;
t63 = t113 * t71;
t45 = -t115 * t136 + t63;
t123 = t139 * t74 + t153;
t122 = (-Ifges(6,4) - Ifges(5,5)) * t115 + (Ifges(5,6) - Ifges(6,6)) * t113;
t99 = mrSges(6,1) * t132;
t91 = -pkin(4) * t115 - t130;
t90 = (mrSges(4,1) * t117 + mrSges(4,2) * t120) * qJD(3);
t78 = (mrSges(5,1) * t113 + mrSges(5,2) * t115) * t117;
t76 = t176 * t115 + t130;
t67 = -mrSges(6,3) * t131 + t99;
t64 = (-t129 + t173) * t117;
t60 = (t117 * Ifges(5,5) + (t115 * Ifges(5,1) - t159) * t120) * qJD(3);
t59 = (t117 * Ifges(6,4) + (t115 * Ifges(6,1) + t157) * t120) * qJD(3);
t58 = (t117 * Ifges(5,6) + (-t113 * Ifges(5,2) + t158) * t120) * qJD(3);
t57 = (t117 * Ifges(6,6) + (t113 * Ifges(6,3) + t156) * t120) * qJD(3);
t56 = mrSges(7,1) * t120 - mrSges(7,3) * t66;
t55 = -mrSges(7,2) * t120 + mrSges(7,3) * t65;
t54 = t110 - t61;
t52 = (t129 - t183) * t117;
t44 = t113 * t136 + t154;
t43 = Ifges(7,1) * t85 - Ifges(7,4) * t124;
t42 = Ifges(7,4) * t85 - Ifges(7,2) * t124;
t37 = t135 * t140 - t154;
t36 = -Ifges(7,1) * t72 - Ifges(7,4) * t73;
t35 = -Ifges(7,4) * t72 - Ifges(7,2) * t73;
t32 = (pkin(8) + t183) * t139 + t151;
t31 = t45 + t182;
t21 = Ifges(7,1) * t66 + Ifges(7,4) * t65 + Ifges(7,5) * t120;
t20 = Ifges(7,4) * t66 + Ifges(7,2) * t65 + Ifges(7,6) * t120;
t19 = qJD(4) * t85 - qJD(6) * t51;
t18 = qJD(4) * t124 + qJD(6) * t50;
t16 = t63 + (-pkin(8) * t145 + pkin(9) * t148) * qJD(3) + t182;
t15 = -t154 + (-pkin(9) * t144 + (-pkin(5) + t135) * t117) * qJD(3);
t14 = mrSges(7,2) * t140 + mrSges(7,3) * t29;
t13 = -mrSges(7,1) * t140 - mrSges(7,3) * t28;
t6 = Ifges(7,1) * t28 + Ifges(7,4) * t29 - Ifges(7,5) * t140;
t5 = Ifges(7,4) * t28 + Ifges(7,2) * t29 - Ifges(7,6) * t140;
t4 = qJD(6) * t11 + t116 * t24 + t119 * t25;
t3 = -qJD(6) * t12 - t116 * t25 + t119 * t24;
t2 = -qJD(6) * t9 - t116 * t16 + t119 * t15;
t1 = qJD(6) * t8 + t116 * t15 + t119 * t16;
t10 = [0.2e1 * m(7) * (t11 * t3 + t12 * t4 + t172) + 0.2e1 * m(4) * (-t114 ^ 2 * t118 * t141 + t75 * t49 + t172) + 0.2e1 * (m(6) + m(5)) * (t24 * t46 + t47 * t25 + t172); t11 * t13 + t12 * t14 + t3 * t56 + t4 * t55 + t165 * t47 + t164 * t46 + t163 * t25 + t162 * t24 + (t67 + t186) * t74 + (t78 - t166) * t48 + (-t121 * t90 + (-t121 * mrSges(3,2) + (-t120 * mrSges(4,1) + t117 * mrSges(4,2) - mrSges(3,1)) * t118) * qJD(2)) * t114 + (t153 + t152 + (-t117 * t75 + t120 * t74) * qJD(3)) * mrSges(4,3) + m(6) * (t24 * t54 + t25 * t53 + t31 * t47 + t37 * t46 + t40 * t74 + t48 * t64) + m(7) * (t1 * t12 + t11 * t2 + t3 * t8 + t32 * t74 + t4 * t9 - t48 * t52) + m(5) * (-t24 * t61 + t25 * t62 - t44 * t46 + t45 * t47) - m(4) * pkin(2) * t134 + (t123 * t177 + m(4) * (-t140 * t75 + t123 + t152) / 0.2e1) * t178; 0.2e1 * t40 * t77 + 0.2e1 * t62 * t79 + 0.2e1 * t61 * t80 + 0.2e1 * t54 * t81 + 0.2e1 * t53 * t82 + 0.2e1 * t45 * t86 + 0.2e1 * t44 * t87 + 0.2e1 * t37 * t88 + 0.2e1 * t31 * t89 - 0.2e1 * pkin(2) * t90 + t66 * t6 + 0.2e1 * t64 * t67 + 0.2e1 * t52 * t7 + 0.2e1 * t1 * t55 + 0.2e1 * t2 * t56 + t65 * t5 - 0.2e1 * t32 * t30 + t28 * t21 + t29 * t20 + 0.2e1 * t8 * t13 + 0.2e1 * t9 * t14 + (t1 * t9 + t2 * t8 - t32 * t52) * t179 + (t31 * t53 + t37 * t54 + t40 * t64) * t180 + (t44 * t61 + t45 * t62) * t181 + ((t78 * t178 + 0.2e1 * (Ifges(4,4) + t122) * t120 + (-(2 * Ifges(6,2)) - (2 * Ifges(5,3)) - (2 * Ifges(7,3)) - (2 * Ifges(4,2)) + (2 * Ifges(4,1)) + (pkin(8) ^ 2 * t181) + t184 * t112 + ((Ifges(5,2) + Ifges(6,3)) * t113 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t115) * t113) * t117) * qJD(3) + t169) * t120 + (t68 * t178 + (t59 + t60) * t115 + (t57 - t58) * t113 + (-Ifges(7,5) * t66 - Ifges(7,6) * t65 + (-0.2e1 * Ifges(4,4) - t122) * t117) * qJD(3)) * t117; -t49 * mrSges(4,2) - t74 * t34 + (t161 - t168) * t48 + m(6) * (t48 * t91 + (-qJD(5) * t74 + t125) * t113 + t170) + m(7) * (t11 * t19 + t12 * t18 - t138 * t74 + t3 * t50 + t4 * t51 - t48 * t76) + m(5) * (-pkin(3) * t48 + t113 * t125 + t170) + (t11 * t72 - t12 * t73 - t124 * t4 - t3 * t85) * mrSges(7,3) + t185 * (t113 * t24 + t155); t40 * t95 + t76 * t7 - t124 * t5 / 0.2e1 + t85 * t6 / 0.2e1 + t66 * t36 / 0.2e1 - pkin(3) * t68 - t72 * t21 / 0.2e1 - t73 * t20 / 0.2e1 + t52 * t34 + t18 * t55 + t19 * t56 + t65 * t35 / 0.2e1 + t29 * t42 / 0.2e1 + t28 * t43 / 0.2e1 + t50 * t13 + t51 * t14 - t32 * t41 + (-t57 / 0.2e1 + t58 / 0.2e1 + t45 * mrSges(5,3) + t31 * mrSges(6,2) + t163 * qJD(4) + t165 * qJ(4) + m(5) * (qJ(4) * t45 + qJD(4) * t62) + m(6) * (qJ(4) * t31 + qJD(4) * t53)) * t115 + t120 * t167 / 0.2e1 + ((pkin(8) * mrSges(4,2) - Ifges(7,5) * t85 / 0.2e1 + Ifges(7,6) * t124 / 0.2e1 - Ifges(4,6) + (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t115 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t113) * t117 + (-t113 * (Ifges(5,2) * t115 + t159) / 0.2e1 + t113 * (-Ifges(6,3) * t115 + t157) / 0.2e1 + Ifges(4,5) + (-m(5) * pkin(3) + t161) * pkin(8) + (t184 * t113 - t156 + t158) * t115 / 0.2e1) * t120) * qJD(3) + m(7) * (t1 * t51 + t18 * t9 + t19 * t8 + t2 * t50 - t32 * t76) + (-t1 * t124 - t2 * t85 + t8 * t72 - t9 * t73) * mrSges(7,3) + (t59 / 0.2e1 + t60 / 0.2e1 + t37 * mrSges(6,2) - t44 * mrSges(5,3) + t162 * qJD(4) + t164 * qJ(4) + m(5) * (-qJ(4) * t44 - qJD(4) * t61) + m(6) * (qJ(4) * t37 + qJD(4) * t54) + (-m(6) * t64 + m(7) * t52 + t166) * qJD(5)) * t113 + (t67 + t175) * t91; 0.2e1 * t76 * t34 - t124 * t35 + t85 * t36 - t73 * t42 - t72 * t43 + (t138 * t76 + t18 * t51 + t19 * t50) * t179 + t143 * t181 + (-t138 * t91 + t143) * t180 + 0.2e1 * t185 * t187 + 0.2e1 * t168 * t138 + 0.2e1 * (-t124 * t18 - t19 * t85 + t50 * t72 - t51 * t73) * mrSges(7,3); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1 + t177) * t48; t99 + ((m(5) * pkin(8)) - t160) * t139 + t175 + m(7) * t32 + t186; (-m(6) - m(7)) * t138 - t34; 0; m(6) * t24 + m(7) * (t116 * t4 + t119 * t3 + (-t11 * t116 + t119 * t12) * qJD(6)); t116 * t14 + t119 * t13 + (-t116 * t56 + t119 * t55) * qJD(6) + m(7) * (t1 * t116 + t119 * t2 + (-t116 * t8 + t119 * t9) * qJD(6)) + m(6) * t37 + t81; m(7) * (t116 * t18 + t119 * t19 + (-t116 * t50 + t119 * t51) * qJD(6)) + m(6) * t113 * qJD(4) + (-t116 * t73 + t119 * t72 + (t116 * t85 - t119 * t124) * qJD(6)) * mrSges(7,3); 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,3) * t140 + t169; mrSges(7,1) * t19 - mrSges(7,2) * t18 + t167; 0; (-mrSges(7,1) * t116 - mrSges(7,2) * t119) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
