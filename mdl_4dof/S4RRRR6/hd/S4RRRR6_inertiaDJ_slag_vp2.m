% Calculate time derivative of joint inertia matrix for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:23
% EndTime: 2019-12-31 17:29:29
% DurationCPUTime: 1.96s
% Computational Cost: add. (1787->312), mult. (5038->494), div. (0->0), fcn. (4358->8), ass. (0->140)
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t106 = sin(qJ(3));
t132 = qJD(4) * t106;
t109 = cos(qJ(3));
t134 = qJD(3) * t109;
t112 = -t105 * t132 + t108 * t134;
t163 = t105 / 0.2e1;
t151 = t108 / 0.2e1;
t104 = cos(pkin(4));
t110 = cos(qJ(2));
t150 = pkin(1) * t110;
t103 = sin(pkin(4));
t107 = sin(qJ(2));
t139 = t103 * t107;
t96 = pkin(6) * t139;
t73 = t104 * t150 - t96;
t138 = t103 * t110;
t74 = t104 * t107 * pkin(1) + pkin(6) * t138;
t63 = pkin(7) * t104 + t74;
t64 = (-pkin(2) * t110 - pkin(7) * t107 - pkin(1)) * t103;
t146 = t106 * t64 + t109 * t63;
t137 = qJD(2) * t103;
t68 = (pkin(2) * t107 - pkin(7) * t110) * t137;
t69 = t73 * qJD(2);
t14 = -qJD(3) * t146 - t106 * t69 + t109 * t68;
t162 = 2 * m(5);
t161 = 0.2e1 * pkin(7);
t160 = -2 * mrSges(3,3);
t72 = t104 * t106 + t109 * t139;
t113 = t105 * t138 - t108 * t72;
t136 = qJD(2) * t107;
t126 = t103 * t136;
t125 = t110 * t137;
t71 = -t104 * t109 + t106 * t139;
t48 = -qJD(3) * t71 + t109 * t125;
t21 = qJD(4) * t113 - t48 * t105 + t108 * t126;
t159 = t21 / 0.2e1;
t49 = -t105 * t72 - t108 * t138;
t158 = t49 / 0.2e1;
t157 = -t113 / 0.2e1;
t142 = Ifges(5,4) * t105;
t118 = Ifges(5,1) * t108 - t142;
t67 = -Ifges(5,5) * t109 + t106 * t118;
t156 = t67 / 0.2e1;
t155 = Ifges(5,5) * t163 + Ifges(5,6) * t151;
t141 = Ifges(5,4) * t108;
t91 = Ifges(5,1) * t105 + t141;
t154 = t91 / 0.2e1;
t153 = -t105 / 0.2e1;
t152 = -t108 / 0.2e1;
t149 = pkin(7) * t109;
t148 = t69 * mrSges(3,2);
t70 = t74 * qJD(2);
t147 = t70 * mrSges(3,1);
t145 = mrSges(5,3) * t106;
t144 = Ifges(4,4) * t106;
t143 = Ifges(4,4) * t109;
t140 = Ifges(5,6) * t105;
t135 = qJD(3) * t106;
t133 = qJD(4) * t105;
t131 = qJD(4) * t108;
t130 = qJD(4) * t109;
t22 = qJD(4) * t49 + t105 * t126 + t48 * t108;
t47 = qJD(3) * t72 + t106 * t125;
t3 = Ifges(5,5) * t22 + Ifges(5,6) * t21 + Ifges(5,3) * t47;
t129 = Ifges(4,5) * t48 - Ifges(4,6) * t47 + Ifges(4,3) * t126;
t127 = Ifges(4,6) * t138;
t85 = (pkin(3) * t106 - pkin(8) * t109) * qJD(3);
t86 = -pkin(3) * t109 - pkin(8) * t106 - pkin(2);
t32 = t86 * t131 + t105 * t85 + (-t105 * t130 - t108 * t135) * pkin(7);
t60 = -t105 * t149 + t108 * t86;
t122 = -qJD(4) * t60 + t32;
t33 = -t86 * t133 + t108 * t85 + (t105 * t135 - t108 * t130) * pkin(7);
t61 = t105 * t86 + t108 * t149;
t121 = -qJD(4) * t61 - t33;
t13 = t106 * t68 + t109 * t69 + t134 * t64 - t135 * t63;
t11 = pkin(8) * t126 + t13;
t18 = t47 * pkin(3) - t48 * pkin(8) + t70;
t62 = t96 + (-pkin(2) - t150) * t104;
t27 = t71 * pkin(3) - t72 * pkin(8) + t62;
t29 = -pkin(8) * t138 + t146;
t7 = -t105 * t29 + t108 * t27;
t1 = qJD(4) * t7 + t105 * t18 + t108 * t11;
t8 = t105 * t27 + t108 * t29;
t2 = -qJD(4) * t8 - t105 * t11 + t108 * t18;
t120 = t1 * t108 - t2 * t105;
t87 = -mrSges(5,1) * t108 + mrSges(5,2) * t105;
t119 = mrSges(5,1) * t105 + mrSges(5,2) * t108;
t117 = -Ifges(5,2) * t105 + t141;
t89 = Ifges(5,2) * t108 + t142;
t34 = -t106 * t63 + t109 * t64;
t16 = -Ifges(5,4) * t113 + Ifges(5,2) * t49 + Ifges(5,6) * t71;
t17 = -Ifges(5,1) * t113 + Ifges(5,4) * t49 + Ifges(5,5) * t71;
t114 = t151 * t17 + t153 * t16;
t111 = t105 * t134 + t106 * t131;
t40 = Ifges(5,5) * t112 - Ifges(5,6) * t111 + Ifges(5,3) * t135;
t102 = Ifges(4,5) * t134;
t101 = Ifges(5,5) * t131;
t94 = Ifges(3,5) * t125;
t92 = Ifges(4,1) * t106 + t143;
t90 = Ifges(4,2) * t109 + t144;
t84 = -mrSges(5,1) * t109 - t108 * t145;
t83 = mrSges(5,2) * t109 - t105 * t145;
t82 = (Ifges(4,1) * t109 - t144) * qJD(3);
t81 = t118 * qJD(4);
t80 = (-Ifges(4,2) * t106 + t143) * qJD(3);
t79 = t117 * qJD(4);
t78 = -Ifges(5,6) * t133 + t101;
t77 = (mrSges(4,1) * t106 + mrSges(4,2) * t109) * qJD(3);
t76 = t119 * qJD(4);
t75 = t119 * t106;
t66 = -Ifges(5,6) * t109 + t106 * t117;
t65 = -Ifges(5,3) * t109 + (Ifges(5,5) * t108 - t140) * t106;
t55 = -mrSges(5,2) * t135 - mrSges(5,3) * t111;
t54 = mrSges(5,1) * t135 - mrSges(5,3) * t112;
t52 = -mrSges(4,1) * t138 - mrSges(4,3) * t72;
t51 = mrSges(4,2) * t138 - mrSges(4,3) * t71;
t43 = mrSges(5,1) * t111 + mrSges(5,2) * t112;
t42 = -t91 * t132 + (Ifges(5,5) * t106 + t109 * t118) * qJD(3);
t41 = -t89 * t132 + (Ifges(5,6) * t106 + t109 * t117) * qJD(3);
t39 = mrSges(4,1) * t126 - mrSges(4,3) * t48;
t38 = -mrSges(4,2) * t126 - mrSges(4,3) * t47;
t37 = Ifges(4,1) * t72 - Ifges(4,4) * t71 - Ifges(4,5) * t138;
t36 = Ifges(4,4) * t72 - Ifges(4,2) * t71 - t127;
t31 = mrSges(5,1) * t71 + mrSges(5,3) * t113;
t30 = -mrSges(5,2) * t71 + mrSges(5,3) * t49;
t28 = pkin(3) * t138 - t34;
t26 = -mrSges(5,1) * t49 - mrSges(5,2) * t113;
t25 = mrSges(4,1) * t47 + mrSges(4,2) * t48;
t24 = Ifges(4,1) * t48 - Ifges(4,4) * t47 + Ifges(4,5) * t126;
t23 = Ifges(4,4) * t48 - Ifges(4,2) * t47 + Ifges(4,6) * t126;
t15 = -Ifges(5,5) * t113 + Ifges(5,6) * t49 + Ifges(5,3) * t71;
t12 = -pkin(3) * t126 - t14;
t10 = mrSges(5,1) * t47 - mrSges(5,3) * t22;
t9 = -mrSges(5,2) * t47 + mrSges(5,3) * t21;
t6 = -mrSges(5,1) * t21 + mrSges(5,2) * t22;
t5 = Ifges(5,1) * t22 + Ifges(5,4) * t21 + Ifges(5,5) * t47;
t4 = Ifges(5,4) * t22 + Ifges(5,2) * t21 + Ifges(5,6) * t47;
t19 = [t71 * t3 - t71 * t23 + 0.2e1 * t70 * (mrSges(4,1) * t71 + mrSges(4,2) * t72) + t72 * t24 + 0.2e1 * t14 * t52 + 0.2e1 * t62 * t25 + t48 * t37 + t49 * t4 + 0.2e1 * t13 * t51 + 0.2e1 * t34 * t39 + 0.2e1 * t12 * t26 + 0.2e1 * t28 * t6 + 0.2e1 * t1 * t30 + 0.2e1 * t2 * t31 + t21 * t16 + t22 * t17 + 0.2e1 * t8 * t9 + 0.2e1 * t7 * t10 + (-0.2e1 * t147 + t94 - 0.2e1 * t148) * t104 - t113 * t5 + 0.2e1 * m(3) * (t69 * t74 - t70 * t73) + (t1 * t8 + t12 * t28 + t2 * t7) * t162 + (t15 - t36) * t47 + (-t110 * t129 + 0.2e1 * (t107 * t70 + t110 * t69) * mrSges(3,3) + ((t73 * t160 + Ifges(3,5) * t104 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t110) * t103) * t110 + (t74 * t160 + Ifges(4,5) * t72 - 0.2e1 * Ifges(3,6) * t104 - Ifges(4,6) * t71 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t107 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t110) * t103) * t107) * qJD(2)) * t103 + 0.2e1 * t146 * t38 + 0.2e1 * m(4) * (t13 * t146 + t14 * t34 + t62 * t70); (-t90 / 0.2e1 + t65 / 0.2e1) * t47 + m(5) * (t1 * t61 + t2 * t60 + t32 * t8 + t33 * t7) + t94 + (-m(4) * t70 - t25) * pkin(2) - t147 - t148 + t12 * t75 + t62 * t77 + t72 * t82 / 0.2e1 + t1 * t83 + t2 * t84 + t48 * t92 / 0.2e1 + t7 * t54 + t8 * t55 + t60 * t10 + t61 * t9 + t28 * t43 + t32 * t30 + t33 * t31 + (t70 * mrSges(4,2) + t24 / 0.2e1 + t5 * t151 + t4 * t153 - t14 * mrSges(4,3) + (t152 * t16 + t153 * t17) * qJD(4)) * t106 + (-t110 * t102 / 0.2e1 + (Ifges(4,5) * t106 / 0.2e1 + Ifges(4,6) * t109 / 0.2e1 - Ifges(3,6)) * t136) * t103 + t41 * t158 + t66 * t159 + t22 * t156 + t42 * t157 + (t109 * t38 + (-t39 + t6) * t106 + (-t106 * t51 + (t26 - t52) * t109) * qJD(3) + m(4) * (-t14 * t106 + t13 * t109 - t134 * t34 - t135 * t146) + m(5) * (t106 * t12 + t134 * t28)) * pkin(7) + ((-t34 * mrSges(4,3) + t37 / 0.2e1 + t114) * t109 + (-t146 * mrSges(4,3) + t127 / 0.2e1 - t36 / 0.2e1 + t15 / 0.2e1) * t106) * qJD(3) + (-t70 * mrSges(4,1) + t23 / 0.2e1 - t3 / 0.2e1 + t13 * mrSges(4,3)) * t109 + (-t80 / 0.2e1 + t40 / 0.2e1) * t71; 0.2e1 * t32 * t83 + 0.2e1 * t61 * t55 + 0.2e1 * t33 * t84 + 0.2e1 * t60 * t54 + (t32 * t61 + t33 * t60) * t162 - 0.2e1 * pkin(2) * t77 + (-t40 + t80 + (-t105 * t66 + t108 * t67 + t161 * t75 + t92) * qJD(3)) * t109 + (t43 * t161 - t105 * t41 + t108 * t42 + t82 + (-t105 * t67 - t108 * t66) * qJD(4) + (pkin(7) ^ 2 * t109 * t162 + t65 - t90) * qJD(3)) * t106; t4 * t151 + t5 * t163 + t28 * t76 + t71 * t78 / 0.2e1 + t79 * t158 + t81 * t157 + t12 * t87 + t47 * t155 + t89 * t159 + t22 * t154 - t13 * mrSges(4,2) + t14 * mrSges(4,1) + t114 * qJD(4) + (-m(5) * t12 - t6) * pkin(3) + ((-t105 * t8 - t108 * t7) * qJD(4) + t120) * mrSges(5,3) + (-t31 * t131 - t30 * t133 + m(5) * (-t131 * t7 - t133 * t8 + t120) + t108 * t9 - t105 * t10) * pkin(8) + t129; -pkin(3) * t43 + t102 + (-t78 / 0.2e1 + (-m(5) * pkin(3) - mrSges(4,1) + t87) * qJD(3) * pkin(7)) * t109 + (t134 * t154 + t41 / 0.2e1 + qJD(4) * t156 + t122 * mrSges(5,3) + (m(5) * t122 - qJD(4) * t84 + t55) * pkin(8)) * t108 + (-t89 * t134 / 0.2e1 + t42 / 0.2e1 - qJD(4) * t66 / 0.2e1 + t121 * mrSges(5,3) + (m(5) * t121 - qJD(4) * t83 - t54) * pkin(8)) * t105 + (t81 * t151 + t79 * t153 + pkin(7) * t76 + (t152 * t89 + t153 * t91) * qJD(4) + (mrSges(4,2) * pkin(7) - Ifges(4,6) + t155) * qJD(3)) * t106; -0.2e1 * pkin(3) * t76 + t105 * t81 + t108 * t79 + (-t105 * t89 + t108 * t91) * qJD(4); mrSges(5,1) * t2 - mrSges(5,2) * t1 + t3; mrSges(5,1) * t33 - mrSges(5,2) * t32 + t40; t101 + (pkin(8) * t87 - t140) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t19(1), t19(2), t19(4), t19(7); t19(2), t19(3), t19(5), t19(8); t19(4), t19(5), t19(6), t19(9); t19(7), t19(8), t19(9), t19(10);];
Mq = res;
