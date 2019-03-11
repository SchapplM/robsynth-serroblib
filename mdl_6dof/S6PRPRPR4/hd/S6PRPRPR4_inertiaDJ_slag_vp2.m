% Calculate time derivative of joint inertia matrix for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:39:07
% EndTime: 2019-03-08 19:39:10
% DurationCPUTime: 1.80s
% Computational Cost: add. (3058->306), mult. (7696->481), div. (0->0), fcn. (7798->12), ass. (0->127)
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t108 = sin(pkin(11));
t146 = pkin(8) + qJ(3);
t96 = t146 * t108;
t111 = cos(pkin(11));
t98 = t146 * t111;
t156 = -t114 * t98 - t117 * t96;
t112 = cos(pkin(6));
t109 = sin(pkin(6));
t115 = sin(qJ(2));
t132 = t109 * t115;
t83 = -t108 * t132 + t111 * t112;
t84 = t108 * t112 + t111 * t132;
t155 = -t108 * t83 + t111 * t84;
t107 = sin(pkin(12));
t110 = cos(pkin(12));
t113 = sin(qJ(6));
t116 = cos(qJ(6));
t119 = t107 * t113 - t110 * t116;
t85 = t119 * qJD(6);
t154 = 2 * m(6);
t153 = 2 * m(7);
t75 = -t114 * t96 + t117 * t98;
t93 = t108 * t117 + t111 * t114;
t49 = qJD(3) * t93 + qJD(4) * t75;
t152 = 0.2e1 * t49;
t105 = t110 ^ 2;
t91 = t108 * t114 - t117 * t111;
t87 = t91 * qJD(4);
t137 = t110 * t87;
t140 = t107 * t87;
t51 = -mrSges(6,1) * t140 - mrSges(6,2) * t137;
t92 = t107 * t116 + t110 * t113;
t86 = t92 * qJD(6);
t26 = t119 * t87 - t86 * t93;
t27 = t93 * t85 + t87 * t92;
t9 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t151 = t51 + t9;
t150 = t49 * t156;
t118 = cos(qJ(2));
t127 = qJD(2) * t118;
t124 = t109 * t127;
t53 = t114 * t83 + t117 * t84;
t34 = qJD(4) * t53 + t124 * t93;
t52 = t114 * t84 - t117 * t83;
t14 = t52 * t34;
t88 = t93 * qJD(4);
t149 = t88 * Ifges(6,5);
t148 = t88 * Ifges(6,6);
t147 = -mrSges(6,1) * t110 + mrSges(6,2) * t107 - mrSges(5,1);
t145 = pkin(9) + qJ(5);
t44 = pkin(4) * t88 + qJ(5) * t87 - qJD(5) * t93;
t48 = -t91 * qJD(3) + t156 * qJD(4);
t18 = t107 * t44 + t110 * t48;
t101 = -pkin(3) * t111 - pkin(2);
t65 = pkin(4) * t91 - qJ(5) * t93 + t101;
t32 = t107 * t65 + t110 * t75;
t144 = -Ifges(7,5) * t85 - Ifges(7,6) * t86;
t143 = Ifges(6,4) * t107;
t142 = Ifges(6,4) * t110;
t141 = Ifges(6,2) * t107;
t139 = t107 * t93;
t136 = t110 * t93;
t133 = t109 ^ 2 * t115;
t131 = t109 * t118;
t129 = t108 ^ 2 + t111 ^ 2;
t128 = qJD(2) * t109;
t126 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t88;
t125 = t115 * t128;
t62 = t88 * mrSges(5,1) - t87 * mrSges(5,2);
t17 = -t107 * t48 + t110 * t44;
t31 = -t107 * t75 + t110 * t65;
t122 = -t156 * t34 + t49 * t52;
t121 = Ifges(6,5) * t110 - Ifges(6,6) * t107;
t33 = -qJD(4) * t52 - t124 * t91;
t28 = -t107 * t33 + t110 * t125;
t29 = t107 * t125 + t110 * t33;
t120 = -t107 * t28 + t110 * t29;
t19 = pkin(5) * t91 - pkin(9) * t136 + t31;
t22 = -pkin(9) * t139 + t32;
t5 = -t113 * t22 + t116 * t19;
t6 = t113 * t19 + t116 * t22;
t42 = -t107 * t53 - t110 * t131;
t43 = -t107 * t131 + t110 * t53;
t12 = -t113 * t43 + t116 * t42;
t13 = t113 * t42 + t116 * t43;
t94 = t145 * t107;
t97 = t145 * t110;
t72 = -t113 * t97 - t116 * t94;
t74 = -t113 * t94 + t116 * t97;
t100 = -pkin(5) * t110 - pkin(4);
t71 = Ifges(7,1) * t92 - Ifges(7,4) * t119;
t70 = Ifges(7,4) * t92 - Ifges(7,2) * t119;
t69 = mrSges(7,1) * t119 + mrSges(7,2) * t92;
t67 = mrSges(6,1) * t91 - mrSges(6,3) * t136;
t66 = -mrSges(6,2) * t91 - mrSges(6,3) * t139;
t64 = -Ifges(7,1) * t85 - Ifges(7,4) * t86;
t63 = -Ifges(7,4) * t85 - Ifges(7,2) * t86;
t61 = mrSges(7,1) * t86 - mrSges(7,2) * t85;
t60 = (mrSges(6,1) * t107 + mrSges(6,2) * t110) * t93;
t59 = mrSges(6,1) * t88 + mrSges(6,3) * t137;
t58 = -mrSges(6,2) * t88 + mrSges(6,3) * t140;
t55 = t119 * t93;
t54 = t92 * t93;
t50 = pkin(5) * t139 - t156;
t47 = -qJD(5) * t92 - qJD(6) * t74;
t46 = -qJD(5) * t119 + qJD(6) * t72;
t39 = mrSges(7,1) * t91 + t55 * mrSges(7,3);
t38 = -mrSges(7,2) * t91 - t54 * mrSges(7,3);
t37 = t149 - (Ifges(6,1) * t110 - t143) * t87;
t36 = t148 - (-t141 + t142) * t87;
t35 = -pkin(5) * t140 + t49;
t30 = mrSges(7,1) * t54 - mrSges(7,2) * t55;
t21 = -Ifges(7,1) * t55 - Ifges(7,4) * t54 + Ifges(7,5) * t91;
t20 = -Ifges(7,4) * t55 - Ifges(7,2) * t54 + Ifges(7,6) * t91;
t16 = -mrSges(7,2) * t88 + t27 * mrSges(7,3);
t15 = mrSges(7,1) * t88 - t26 * mrSges(7,3);
t11 = pkin(9) * t140 + t18;
t10 = pkin(5) * t88 + pkin(9) * t137 + t17;
t8 = Ifges(7,1) * t26 + Ifges(7,4) * t27 + t88 * Ifges(7,5);
t7 = Ifges(7,4) * t26 + Ifges(7,2) * t27 + t88 * Ifges(7,6);
t4 = -qJD(6) * t13 - t113 * t29 + t116 * t28;
t3 = qJD(6) * t12 + t113 * t28 + t116 * t29;
t2 = -qJD(6) * t6 + t10 * t116 - t11 * t113;
t1 = qJD(6) * t5 + t10 * t113 + t11 * t116;
t23 = [0.2e1 * m(6) * (t28 * t42 + t29 * t43 + t14) + 0.2e1 * m(7) * (t12 * t4 + t13 * t3 + t14) + 0.2e1 * m(4) * (t155 * t109 - t133) * t127 + 0.2e1 * (-t127 * t133 + t53 * t33 + t14) * m(5); t12 * t15 + t13 * t16 + t28 * t67 + t29 * t66 + t3 * t38 + t4 * t39 + t42 * t59 + t43 * t58 + t151 * t52 + (t30 + t60) * t34 + (-t33 * t91 + t34 * t93 - t52 * t87 - t53 * t88) * mrSges(5,3) + (-t118 * t62 + ((mrSges(4,3) * t129 - mrSges(3,2)) * t118 + (-mrSges(4,1) * t111 + mrSges(5,1) * t91 + mrSges(4,2) * t108 + mrSges(5,2) * t93 - mrSges(3,1)) * t115) * qJD(2)) * t109 + m(7) * (t1 * t13 + t12 * t2 + t3 * t6 + t34 * t50 + t35 * t52 + t4 * t5) + m(6) * (t17 * t42 + t18 * t43 + t31 * t28 + t32 * t29 + t122) + m(5) * (t101 * t125 + t33 * t75 + t48 * t53 + t122) + m(4) * (t155 * qJD(3) + (qJ(3) * t118 * t129 - pkin(2) * t115) * t128); 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * t129 * qJD(3) + (t31 * t17 + t32 * t18 - t150) * t154 + 0.2e1 * m(5) * (t48 * t75 - t150) + t60 * t152 + (t1 * t6 + t2 * t5 + t35 * t50) * t153 + 0.2e1 * (t156 * t87 - t75 * t88) * mrSges(5,3) - 0.2e1 * t156 * t51 + (-0.2e1 * mrSges(5,3) * t48 + ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t88 + t126 - 0.2e1 * (-Ifges(5,4) + t121) * t87) * t91 + t88 * (-Ifges(7,5) * t55 - Ifges(7,6) * t54) + 0.2e1 * t5 * t15 + 0.2e1 * t6 * t16 + t26 * t21 + t27 * t20 + 0.2e1 * t35 * t30 + 0.2e1 * t1 * t38 + 0.2e1 * t2 * t39 + 0.2e1 * t50 * t9 - t54 * t7 - t55 * t8 + 0.2e1 * t32 * t58 + 0.2e1 * t31 * t59 + 0.2e1 * t18 * t66 + 0.2e1 * t17 * t67 + 0.2e1 * t101 * t62 + (mrSges(5,3) * t152 - t107 * t36 + t110 * t37 + (-0.2e1 * Ifges(5,4) + t121) * t88 - (Ifges(6,1) * t105 + (2 * Ifges(5,1)) + (t141 - 0.2e1 * t142) * t107) * t87) * t93; (m(4) + m(5)) * t125 + m(6) * (t107 * t29 + t110 * t28) + m(7) * (-t119 * t4 - t12 * t86 - t13 * t85 + t3 * t92); t107 * t58 + t110 * t59 - t119 * t15 + t92 * t16 - t85 * t38 - t86 * t39 + m(7) * (t1 * t92 - t119 * t2 - t5 * t86 - t6 * t85) + m(6) * (t107 * t18 + t110 * t17) + t62; (t119 * t86 - t85 * t92) * t153; -t33 * mrSges(5,2) + t52 * t61 + t120 * mrSges(6,3) + (t69 + t147) * t34 + m(6) * (-pkin(4) * t34 + (-t107 * t42 + t110 * t43) * qJD(5) + t120 * qJ(5)) + m(7) * (t100 * t34 + t47 * t12 + t46 * t13 + t3 * t74 + t4 * t72) + (-t119 * t3 + t12 * t85 - t13 * t86 - t4 * t92) * mrSges(7,3); m(6) * (-pkin(4) * t49 + (-t107 * t31 + t110 * t32) * qJD(5) + (-t17 * t107 + t18 * t110) * qJ(5)) - t119 * t7 / 0.2e1 + (-t1 * t119 - t2 * t92 + t5 * t85 - t6 * t86) * mrSges(7,3) + t88 * (Ifges(7,5) * t92 - Ifges(7,6) * t119) / 0.2e1 + m(7) * (t1 * t74 + t100 * t35 + t2 * t72 + t46 * t6 + t47 * t5) + (t18 * mrSges(6,3) + qJ(5) * t58 + qJD(5) * t66 + t36 / 0.2e1 + t148 / 0.2e1) * t110 + (-t17 * mrSges(6,3) - qJ(5) * t59 - qJD(5) * t67 + t37 / 0.2e1 + t149 / 0.2e1) * t107 - (t110 * (Ifges(6,1) * t107 + t142) / 0.2e1 - t107 * (Ifges(6,2) * t110 + t143) / 0.2e1 + Ifges(5,5)) * t87 + t91 * t144 / 0.2e1 + t147 * t49 + t46 * t38 + t47 * t39 - t48 * mrSges(5,2) - pkin(4) * t51 + t50 * t61 - t54 * t63 / 0.2e1 - t55 * t64 / 0.2e1 + t35 * t69 + t27 * t70 / 0.2e1 + t26 * t71 / 0.2e1 + t72 * t15 + t74 * t16 - t85 * t21 / 0.2e1 - t86 * t20 / 0.2e1 - Ifges(5,6) * t88 + t92 * t8 / 0.2e1 + t100 * t9; m(7) * (-t119 * t47 + t46 * t92 - t72 * t86 - t74 * t85); -t85 * t71 + t92 * t64 - t86 * t70 - t119 * t63 + 0.2e1 * t100 * t61 + (t46 * t74 + t47 * t72) * t153 + 0.2e1 * (-t119 * t46 - t47 * t92 + t72 * t85 - t74 * t86) * mrSges(7,3) + (qJ(5) * t154 + 0.2e1 * mrSges(6,3)) * qJD(5) * (t107 ^ 2 + t105); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t34; m(6) * t49 + m(7) * t35 + t151; 0; t61; 0; mrSges(7,1) * t4 - mrSges(7,2) * t3; t2 * mrSges(7,1) - t1 * mrSges(7,2) + t126; -t61; t47 * mrSges(7,1) - t46 * mrSges(7,2) + t144; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;
