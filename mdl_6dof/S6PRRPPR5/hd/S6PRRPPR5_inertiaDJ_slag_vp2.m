% Calculate time derivative of joint inertia matrix for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:13
% EndTime: 2019-03-08 21:18:17
% DurationCPUTime: 1.87s
% Computational Cost: add. (1884->316), mult. (4642->484), div. (0->0), fcn. (3994->10), ass. (0->138)
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t89 = t114 * mrSges(5,2) - t111 * mrSges(5,3);
t170 = -t114 * mrSges(4,1) + t111 * mrSges(4,2) + t89;
t108 = cos(pkin(6));
t106 = sin(pkin(6));
t115 = cos(qJ(2));
t142 = qJD(2) * t115;
t134 = t106 * t142;
t112 = sin(qJ(2));
t148 = t106 * t112;
t136 = t111 * t148;
t47 = -qJD(3) * t136 + (qJD(3) * t108 + t134) * t114;
t67 = t108 * t111 + t114 * t148;
t58 = t67 * qJD(4);
t169 = qJ(4) * t47 + t58;
t105 = sin(pkin(11));
t107 = cos(pkin(11));
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t168 = -t105 * t110 + t107 * t113;
t123 = t113 * t105 + t110 * t107;
t65 = t123 * qJD(6);
t104 = t107 ^ 2;
t131 = (t105 ^ 2 + t104) * qJD(5);
t167 = 2 * m(6);
t166 = 2 * m(7);
t165 = m(5) / 0.2e1;
t164 = -t123 / 0.2e1;
t163 = t168 / 0.2e1;
t162 = pkin(4) + pkin(8);
t161 = t107 / 0.2e1;
t141 = qJD(3) * t111;
t133 = t107 * t141;
t156 = mrSges(6,2) * t105;
t57 = -mrSges(6,1) * t133 + t141 * t156;
t54 = t168 * t114;
t29 = -qJD(6) * t54 + t123 * t141;
t30 = t114 * t65 + t141 * t168;
t9 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t160 = t57 + t9;
t100 = t114 * pkin(8);
t22 = t67 * t47;
t109 = -pkin(3) - qJ(5);
t159 = -pkin(9) + t109;
t55 = t123 * t114;
t31 = mrSges(7,1) * t54 - mrSges(7,2) * t55;
t70 = (mrSges(6,1) * t107 - t156) * t114;
t158 = t31 + t70;
t132 = pkin(3) * t141 - qJD(4) * t111;
t48 = -qJD(5) * t114 + (-qJ(4) * t114 + qJ(5) * t111) * qJD(3) + t132;
t140 = qJD(3) * t114;
t83 = t162 * t140;
t19 = t105 * t83 + t107 * t48;
t151 = qJ(4) * t111;
t73 = t109 * t114 - pkin(2) - t151;
t90 = t162 * t111;
t37 = t105 * t90 + t107 * t73;
t64 = t168 * qJD(6);
t157 = -Ifges(7,5) * t65 - Ifges(7,6) * t64;
t155 = Ifges(6,4) * t105;
t154 = Ifges(6,4) * t107;
t152 = t105 * Ifges(6,1);
t149 = t105 * t111;
t147 = t106 * t115;
t145 = t107 * t114;
t91 = t114 * pkin(4) + t100;
t143 = qJD(2) * t112;
t139 = 2 * mrSges(7,3);
t41 = mrSges(7,1) * t123 + mrSges(7,2) * t168;
t88 = mrSges(6,1) * t105 + mrSges(6,2) * t107;
t138 = mrSges(5,3) + t41 + t88;
t137 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t140;
t135 = t106 * t143;
t33 = mrSges(7,1) * t64 - t65 * mrSges(7,2);
t18 = -t105 * t48 + t107 * t83;
t130 = -t123 * t64 + t168 * t65;
t128 = -Ifges(6,5) * t105 - Ifges(6,6) * t107;
t127 = -pkin(3) * t114 - t151;
t126 = t19 * t105 + t18 * t107;
t46 = qJD(3) * t67 + t111 * t134;
t23 = -t105 * t135 + t107 * t46;
t24 = t105 * t46 + t107 * t135;
t125 = t105 * t24 + t107 * t23;
t77 = t107 * t90;
t28 = pkin(5) * t111 + t77 + (pkin(9) * t114 - t73) * t105;
t32 = -pkin(9) * t145 + t37;
t7 = -t110 * t32 + t113 * t28;
t8 = t110 * t28 + t113 * t32;
t66 = -t108 * t114 + t136;
t44 = t105 * t147 + t66 * t107;
t45 = t66 * t105 - t107 * t147;
t10 = -t110 * t45 + t113 * t44;
t11 = t110 * t44 + t113 * t45;
t84 = t159 * t105;
t85 = t159 * t107;
t40 = t110 * t85 + t113 * t84;
t39 = -t110 * t84 + t113 * t85;
t124 = t46 * t111 + t47 * t114;
t12 = (pkin(5) * t114 - pkin(9) * t149) * qJD(3) + t18;
t13 = pkin(9) * t133 + t19;
t1 = qJD(6) * t7 + t110 * t12 + t113 * t13;
t2 = -qJD(6) * t8 - t110 * t13 + t113 * t12;
t120 = -t1 * t123 - t168 * t2 - t64 * t8 + t65 * t7;
t3 = qJD(6) * t10 + t110 * t23 + t113 * t24;
t4 = -qJD(6) * t11 - t110 * t24 + t113 * t23;
t119 = -t10 * t65 + t11 * t64 + t123 * t3 + t168 * t4;
t16 = -qJD(5) * t123 + qJD(6) * t39;
t17 = -qJD(5) * t168 - qJD(6) * t40;
t118 = -t123 * t16 - t168 * t17 + t39 * t65 - t40 * t64;
t94 = pkin(5) * t105 + qJ(4);
t86 = -pkin(2) + t127;
t82 = t162 * t141;
t81 = (mrSges(4,1) * t111 + mrSges(4,2) * t114) * qJD(3);
t80 = (-mrSges(5,2) * t111 - mrSges(5,3) * t114) * qJD(3);
t79 = -mrSges(6,2) * t111 - mrSges(6,3) * t145;
t78 = mrSges(6,3) * t105 * t114 + mrSges(6,1) * t111;
t72 = (mrSges(6,3) * t107 * t111 - mrSges(6,2) * t114) * qJD(3);
t71 = (mrSges(6,1) * t114 - mrSges(6,3) * t149) * qJD(3);
t63 = pkin(5) * t145 + t91;
t62 = -qJ(4) * t140 + t132;
t53 = (-pkin(5) * t107 - t162) * t141;
t52 = (t114 * Ifges(6,5) + (t152 + t154) * t111) * qJD(3);
t51 = (t114 * Ifges(6,6) + (t107 * Ifges(6,2) + t155) * t111) * qJD(3);
t50 = mrSges(7,1) * t111 + mrSges(7,3) * t55;
t49 = -mrSges(7,2) * t111 - mrSges(7,3) * t54;
t43 = Ifges(7,1) * t168 - Ifges(7,4) * t123;
t42 = Ifges(7,4) * t168 - Ifges(7,2) * t123;
t36 = -t105 * t73 + t77;
t35 = -Ifges(7,1) * t65 - Ifges(7,4) * t64;
t34 = -Ifges(7,4) * t65 - Ifges(7,2) * t64;
t21 = -Ifges(7,1) * t55 - Ifges(7,4) * t54 + Ifges(7,5) * t111;
t20 = -Ifges(7,4) * t55 - Ifges(7,2) * t54 + Ifges(7,6) * t111;
t15 = -mrSges(7,2) * t140 + mrSges(7,3) * t30;
t14 = mrSges(7,1) * t140 - mrSges(7,3) * t29;
t6 = Ifges(7,1) * t29 + Ifges(7,4) * t30 + Ifges(7,5) * t140;
t5 = Ifges(7,4) * t29 + Ifges(7,2) * t30 + Ifges(7,6) * t140;
t25 = [0.2e1 * m(7) * (t10 * t4 + t11 * t3 + t22) + 0.2e1 * m(6) * (t23 * t44 + t24 * t45 + t22) + 0.2e1 * (m(5) + m(4)) * (-t106 ^ 2 * t112 * t142 + t66 * t46 + t22); t10 * t14 + t11 * t15 + t23 * t78 + t24 * t79 + t3 * t49 + t4 * t50 + t44 * t71 + t45 * t72 + t160 * t67 + t158 * t47 + m(6) * (t18 * t44 + t19 * t45 + t23 * t36 + t24 * t37 + t47 * t91 - t67 * t82) + m(7) * (t1 * t11 + t10 * t2 + t3 * t8 + t4 * t7 + t47 * t63 + t53 * t67) + 0.2e1 * (m(4) / 0.2e1 + t165) * (t140 * t66 - t141 * t67 + t124) * pkin(8) + (mrSges(4,3) + mrSges(5,1)) * ((-t111 * t67 + t114 * t66) * qJD(3) + t124) + ((-t80 - t81) * t115 + (-t115 * mrSges(3,2) + (-mrSges(3,1) + t170) * t112) * qJD(2) - m(4) * pkin(2) * t143 + 0.2e1 * (-t115 * t62 + t143 * t86) * t165) * t106; 0.2e1 * t91 * t57 + 0.2e1 * t18 * t78 + 0.2e1 * t19 * t79 - 0.2e1 * pkin(2) * t81 - 0.2e1 * t82 * t70 + 0.2e1 * t86 * t80 + 0.2e1 * t36 * t71 + 0.2e1 * t37 * t72 + 0.2e1 * t63 * t9 + 0.2e1 * t1 * t49 + 0.2e1 * t2 * t50 + 0.2e1 * t53 * t31 - t54 * t5 - t55 * t6 + t29 * t21 + t30 * t20 + 0.2e1 * t7 * t14 + 0.2e1 * t8 * t15 + 0.2e1 * (m(5) * t86 + t89) * t62 + (t18 * t36 + t19 * t37 - t82 * t91) * t167 + (t1 * t8 + t2 * t7 + t53 * t63) * t166 + (-t105 * t52 - t107 * t51 + (-Ifges(7,5) * t55 - Ifges(7,6) * t54 + ((2 * Ifges(4,4)) + (2 * Ifges(5,6)) + t128) * t114) * qJD(3)) * t114 + ((0.2e1 * (-Ifges(4,4) - Ifges(5,6) - t128) * t111 + (-t104 * Ifges(6,2) + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + (2 * Ifges(5,2)) - (2 * Ifges(5,3)) + (2 * Ifges(6,3)) + Ifges(7,3) + (-t152 - 0.2e1 * t154) * t105) * t114) * qJD(3) + t137) * t111; t67 * t33 + (-mrSges(4,1) + mrSges(5,2)) * t46 - t125 * mrSges(6,3) + (-mrSges(4,2) + t138) * t47 + m(7) * (t10 * t17 + t11 * t16 + t3 * t40 + t39 * t4 + t47 * t94 + t58) + m(5) * (-pkin(3) * t46 + t169) + m(6) * (t125 * t109 + (-t105 * t45 - t107 * t44) * qJD(5) + t169) - t119 * mrSges(7,3); t94 * t9 - t82 * t88 + t6 * t163 + t5 * t164 + qJ(4) * t57 + t63 * t33 - t64 * t20 / 0.2e1 - t65 * t21 / 0.2e1 + t16 * t49 + t17 * t50 + t53 * t41 - t54 * t34 / 0.2e1 - t55 * t35 / 0.2e1 + t39 * t14 + t40 * t15 + t30 * t42 / 0.2e1 + t29 * t43 / 0.2e1 + m(6) * (-qJ(4) * t82 + t126 * t109 + (-t105 * t37 - t107 * t36) * qJD(5)) + (m(5) * t100 + m(6) * t91 + m(7) * t63 + t114 * mrSges(5,1) + t158) * qJD(4) + ((-pkin(3) * mrSges(5,1) + Ifges(6,5) * t161 - Ifges(6,6) * t105 / 0.2e1 + Ifges(7,5) * t163 + Ifges(7,6) * t164 + Ifges(4,5) - Ifges(5,4)) * t114 + (t105 * (Ifges(6,1) * t107 - t155) / 0.2e1 + (-Ifges(6,2) * t105 + t154) * t161 - qJ(4) * mrSges(5,1) - Ifges(4,6) + Ifges(5,5)) * t111 + (m(5) * t127 + t170) * pkin(8)) * qJD(3) + t111 * t157 / 0.2e1 + m(7) * (t1 * t40 + t16 * t8 + t17 * t7 + t2 * t39 + t53 * t94) + t120 * mrSges(7,3) + (-t51 / 0.2e1 - qJD(5) * t79 - t19 * mrSges(6,3) + t109 * t72) * t105 + (t52 / 0.2e1 - qJD(5) * t78 - t18 * mrSges(6,3) + t109 * t71) * t107; 0.2e1 * t94 * t33 - t123 * t34 + t168 * t35 - t64 * t42 - t65 * t43 + (qJD(4) * t94 + t16 * t40 + t17 * t39) * t166 + (qJ(4) * qJD(4) - t109 * t131) * t167 + t118 * t139 + 0.2e1 * mrSges(6,3) * t131 + 0.2e1 * (m(5) * qJ(4) + t138) * qJD(4); m(5) * t46 + m(6) * t125 + m(7) * t119; t105 * t72 + t107 * t71 + t168 * t14 + t123 * t15 + t64 * t49 - t65 * t50 + (m(5) * pkin(8) + mrSges(5,1)) * t140 - m(7) * t120 + m(6) * t126; -m(6) * t131 - m(7) * t118 + t130 * t139; -0.2e1 * m(7) * t130; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t47; -m(6) * t82 + m(7) * t53 + t160; (m(6) + m(7)) * qJD(4) + t33; 0; 0; mrSges(7,1) * t4 - mrSges(7,2) * t3; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t137; mrSges(7,1) * t17 - mrSges(7,2) * t16 + t157; -mrSges(7,1) * t65 - t64 * mrSges(7,2); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
