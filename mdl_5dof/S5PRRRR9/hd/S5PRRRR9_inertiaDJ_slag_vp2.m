% Calculate time derivative of joint inertia matrix for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:06
% EndTime: 2019-12-05 17:19:14
% DurationCPUTime: 2.22s
% Computational Cost: add. (2049->320), mult. (5460->501), div. (0->0), fcn. (4764->10), ass. (0->146)
t107 = cos(qJ(4));
t108 = cos(qJ(3));
t137 = qJD(3) * t108;
t124 = t107 * t137;
t104 = sin(qJ(3));
t138 = qJD(3) * t104;
t165 = -Ifges(5,5) * t124 - Ifges(5,3) * t138;
t102 = sin(qJ(5));
t103 = sin(qJ(4));
t106 = cos(qJ(5));
t116 = t102 * t103 - t106 * t107;
t64 = t116 * t104;
t86 = -pkin(3) * t108 - pkin(8) * t104 - pkin(2);
t140 = t107 * t108;
t94 = pkin(7) * t140;
t60 = t103 * t86 + t94;
t164 = t60 * qJD(4);
t87 = -mrSges(5,1) * t107 + mrSges(5,2) * t103;
t163 = -m(5) * pkin(3) - mrSges(4,1) + t87;
t162 = qJD(4) + qJD(5);
t161 = 0.2e1 * m(5);
t160 = 2 * m(6);
t159 = 0.2e1 * pkin(7);
t150 = Ifges(5,4) * t103;
t88 = Ifges(5,2) * t107 + t150;
t157 = -t88 / 0.2e1;
t156 = -pkin(9) - pkin(8);
t155 = -t103 / 0.2e1;
t154 = pkin(7) * t103;
t100 = sin(pkin(5));
t109 = cos(qJ(2));
t139 = qJD(2) * t109;
t126 = t100 * t139;
t101 = cos(pkin(5));
t105 = sin(qJ(2));
t144 = t100 * t105;
t66 = t101 * t104 + t108 * t144;
t47 = qJD(3) * t66 + t104 * t126;
t65 = -t101 * t108 + t104 * t144;
t32 = t65 * t47;
t84 = (pkin(3) * t104 - pkin(8) * t108) * qJD(3);
t152 = t107 * t84 + t138 * t154;
t149 = Ifges(5,4) * t107;
t148 = Ifges(5,6) * t103;
t147 = t108 * Ifges(5,6);
t146 = t47 * t104;
t48 = -qJD(3) * t65 + t108 * t126;
t145 = t48 * t108;
t143 = t100 * t109;
t142 = t103 * t104;
t141 = t104 * t107;
t136 = qJD(4) * t103;
t135 = qJD(4) * t104;
t134 = qJD(4) * t107;
t133 = qJD(5) * t102;
t132 = qJD(5) * t106;
t75 = t102 * t107 + t103 * t106;
t43 = t162 * t75;
t25 = -t104 * t43 - t116 * t137;
t26 = -t75 * t137 + t162 * t64;
t131 = -Ifges(6,5) * t25 - Ifges(6,6) * t26 - Ifges(6,3) * t138;
t130 = pkin(4) * t136;
t113 = t103 * t143 - t66 * t107;
t127 = qJD(2) * t144;
t10 = qJD(4) * t113 - t48 * t103 + t107 * t127;
t49 = -t66 * t103 - t107 * t143;
t11 = qJD(4) * t49 + t103 * t127 + t48 * t107;
t17 = t102 * t113 + t106 * t49;
t2 = qJD(5) * t17 + t10 * t102 + t106 * t11;
t18 = t102 * t49 - t106 * t113;
t3 = -qJD(5) * t18 + t10 * t106 - t102 * t11;
t129 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t128 = qJD(4) * t156;
t125 = t103 * t135;
t123 = (2 * Ifges(4,4)) + t148;
t30 = t103 * t84 + t86 * t134 + (-t107 * t138 - t108 * t136) * pkin(7);
t73 = t107 * t86;
t59 = -t108 * t154 + t73;
t122 = -t59 * qJD(4) + t30;
t121 = mrSges(5,1) * t103 + mrSges(5,2) * t107;
t120 = Ifges(5,1) * t107 - t150;
t89 = Ifges(5,1) * t103 + t149;
t119 = -Ifges(5,2) * t103 + t149;
t118 = Ifges(5,5) * t103 + Ifges(5,6) * t107;
t39 = -pkin(9) * t141 + t73 + (-pkin(4) - t154) * t108;
t51 = -pkin(9) * t142 + t60;
t15 = -t102 * t51 + t106 * t39;
t16 = t102 * t39 + t106 * t51;
t90 = t156 * t103;
t91 = t156 * t107;
t52 = t102 * t91 + t106 * t90;
t53 = t102 * t90 - t106 * t91;
t82 = t103 * t128;
t83 = t107 * t128;
t28 = qJD(5) * t52 + t102 * t83 + t106 * t82;
t29 = -qJD(5) * t53 - t102 * t82 + t106 * t83;
t40 = Ifges(6,6) * t43;
t42 = t162 * t116;
t41 = Ifges(6,5) * t42;
t117 = t29 * mrSges(6,1) - t28 * mrSges(6,2) - t40 - t41;
t12 = (pkin(4) * t104 - pkin(9) * t140) * qJD(3) + (-t94 + (pkin(9) * t104 - t86) * t103) * qJD(4) + t152;
t111 = t103 * t137 + t104 * t134;
t22 = -pkin(9) * t111 + t30;
t5 = qJD(5) * t15 + t102 * t12 + t106 * t22;
t6 = -qJD(5) * t16 - t102 * t22 + t106 * t12;
t115 = t6 * mrSges(6,1) - t5 * mrSges(6,2) - t131;
t114 = t137 * t65 + t146;
t112 = t124 - t125;
t99 = Ifges(5,5) * t134;
t96 = -pkin(4) * t107 - pkin(3);
t85 = (pkin(4) * t103 + pkin(7)) * t104;
t81 = -mrSges(5,1) * t108 - mrSges(5,3) * t141;
t80 = mrSges(5,2) * t108 - mrSges(5,3) * t142;
t79 = t120 * qJD(4);
t78 = t119 * qJD(4);
t77 = (mrSges(4,1) * t104 + mrSges(4,2) * t108) * qJD(3);
t76 = t121 * qJD(4);
t71 = (-mrSges(6,1) * t102 - mrSges(6,2) * t106) * qJD(5) * pkin(4);
t69 = t121 * t104;
t63 = t75 * t104;
t62 = -Ifges(5,5) * t108 + t104 * t120;
t61 = t104 * t119 - t147;
t58 = pkin(4) * t111 + pkin(7) * t137;
t57 = -mrSges(5,2) * t138 - mrSges(5,3) * t111;
t56 = mrSges(5,1) * t138 - mrSges(5,3) * t112;
t55 = -mrSges(6,1) * t108 + mrSges(6,3) * t64;
t54 = mrSges(6,2) * t108 - mrSges(6,3) * t63;
t46 = Ifges(6,1) * t75 - Ifges(6,4) * t116;
t45 = Ifges(6,4) * t75 - Ifges(6,2) * t116;
t44 = mrSges(6,1) * t116 + mrSges(6,2) * t75;
t38 = mrSges(5,1) * t111 + mrSges(5,2) * t112;
t37 = mrSges(6,1) * t63 - mrSges(6,2) * t64;
t36 = -t89 * t135 + (Ifges(5,5) * t104 + t108 * t120) * qJD(3);
t35 = -t88 * t135 + (Ifges(5,6) * t104 + t108 * t119) * qJD(3);
t34 = -Ifges(6,1) * t64 - Ifges(6,4) * t63 - Ifges(6,5) * t108;
t33 = -Ifges(6,4) * t64 - Ifges(6,2) * t63 - Ifges(6,6) * t108;
t31 = t152 - t164;
t21 = -Ifges(6,1) * t42 - Ifges(6,4) * t43;
t20 = -Ifges(6,4) * t42 - Ifges(6,2) * t43;
t19 = mrSges(6,1) * t43 - mrSges(6,2) * t42;
t14 = -mrSges(6,2) * t138 + mrSges(6,3) * t26;
t13 = mrSges(6,1) * t138 - mrSges(6,3) * t25;
t9 = -mrSges(6,1) * t26 + mrSges(6,2) * t25;
t8 = Ifges(6,1) * t25 + Ifges(6,4) * t26 + Ifges(6,5) * t138;
t7 = Ifges(6,4) * t25 + Ifges(6,2) * t26 + Ifges(6,6) * t138;
t1 = [0.2e1 * m(6) * (t17 * t3 + t18 * t2 + t32) + 0.2e1 * m(5) * (t10 * t49 - t11 * t113 + t32) + 0.2e1 * m(4) * (-t100 ^ 2 * t105 * t139 + t66 * t48 + t32); t10 * t81 + t11 * t80 + t17 * t13 + t18 * t14 + t2 * t54 + t3 * t55 + t49 * t56 - t113 * t57 + (t9 + t38) * t65 + (t37 + t69) * t47 + (-t109 * t77 + (-t109 * mrSges(3,2) + (-mrSges(4,1) * t108 + mrSges(4,2) * t104 - mrSges(3,1)) * t105) * qJD(2)) * t100 + m(6) * (t15 * t3 + t16 * t2 + t17 * t6 + t18 * t5 + t47 * t85 + t58 * t65) + m(5) * (t10 * t59 + t11 * t60 - t113 * t30 + t31 * t49) - m(4) * pkin(2) * t127 + (m(5) * t114 / 0.2e1 + m(4) * (-t66 * t138 + t114 + t145) / 0.2e1) * t159 + (t146 + t145 + (-t104 * t66 + t108 * t65) * qJD(3)) * mrSges(4,3); -0.2e1 * pkin(2) * t77 + 0.2e1 * t15 * t13 + 0.2e1 * t16 * t14 + t25 * t34 + t26 * t33 + 0.2e1 * t30 * t80 + 0.2e1 * t31 * t81 + 0.2e1 * t58 * t37 + 0.2e1 * t5 * t54 + 0.2e1 * t6 * t55 + 0.2e1 * t59 * t56 + 0.2e1 * t60 * t57 - t63 * t7 - t64 * t8 + 0.2e1 * t85 * t9 + (t30 * t60 + t31 * t59) * t161 + (t15 * t6 + t16 * t5 + t58 * t85) * t160 + ((-t103 * t61 + t107 * t62 + t108 * t123 + t159 * t69) * qJD(3) + t131 + t165) * t108 + (t38 * t159 - t103 * t35 + t107 * t36 + (-t103 * t62 - t107 * t61 + t108 * t118) * qJD(4) + (-Ifges(6,5) * t64 - Ifges(6,6) * t63 + (Ifges(5,5) * t107 - t123) * t104 + (pkin(7) ^ 2 * t161 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) - Ifges(6,3)) * t108) * qJD(3)) * t104; -t48 * mrSges(4,2) + (t19 + t76) * t65 + m(6) * (t130 * t65 + t17 * t29 + t18 * t28 + t2 * t53 + t3 * t52) + (-t116 * t2 + t17 * t42 - t18 * t43 - t3 * t75) * mrSges(6,3) + (m(6) * t96 + t163 + t44) * t47 + (m(5) * pkin(8) + mrSges(5,3)) * (-t10 * t103 + t107 * t11 + (t103 * t113 - t107 * t49) * qJD(4)); t96 * t9 + m(6) * (t15 * t29 + t16 * t28 + t5 * t53 + t52 * t6 + t58 * t96) + t85 * t19 - t116 * t7 / 0.2e1 + t75 * t8 / 0.2e1 - t63 * t20 / 0.2e1 - t64 * t21 / 0.2e1 + t26 * t45 / 0.2e1 + t25 * t46 / 0.2e1 + t52 * t13 + t53 * t14 + t28 * t54 + t29 * t55 + t58 * t44 - pkin(3) * t38 - t42 * t34 / 0.2e1 - t43 * t33 / 0.2e1 + (-t99 / 0.2e1 + t41 / 0.2e1 + t40 / 0.2e1 + (pkin(7) * t163 + Ifges(4,5)) * qJD(3)) * t108 + (t36 / 0.2e1 - t31 * mrSges(5,3) + t137 * t157 + (t147 / 0.2e1 - t61 / 0.2e1 - t60 * mrSges(5,3) + (m(6) * t85 + t37) * pkin(4)) * qJD(4) + (m(5) * (-t31 - t164) - t56 - qJD(4) * t80) * pkin(8)) * t103 + (-t116 * t5 + t15 * t42 - t16 * t43 - t6 * t75) * mrSges(6,3) + (t35 / 0.2e1 + qJD(4) * t62 / 0.2e1 + t89 * t137 / 0.2e1 + t122 * mrSges(5,3) + (m(5) * t122 - qJD(4) * t81 + t57) * pkin(8)) * t107 + (t107 * t79 / 0.2e1 + t78 * t155 - Ifges(4,6) * qJD(3) + (t107 * t157 + t155 * t89) * qJD(4) + (qJD(3) * mrSges(4,2) + t76) * pkin(7) + (Ifges(6,5) * t75 - Ifges(6,6) * t116 + t118) * qJD(3) / 0.2e1) * t104; -t42 * t46 + t75 * t21 - t43 * t45 - t116 * t20 + (t130 * t96 + t28 * t53 + t29 * t52) * t160 + 0.2e1 * t44 * t130 + 0.2e1 * t96 * t19 - 0.2e1 * pkin(3) * t76 + t103 * t79 - t88 * t136 + (qJD(4) * t89 + t78) * t107 + 0.2e1 * (-t116 * t28 - t29 * t75 + t42 * t52 - t43 * t53) * mrSges(6,3); t10 * mrSges(5,1) - t11 * mrSges(5,2) + m(6) * (t102 * t2 + t106 * t3 + (-t102 * t17 + t106 * t18) * qJD(5)) * pkin(4) + t129; -Ifges(5,5) * t125 + t31 * mrSges(5,1) - t30 * mrSges(5,2) - t111 * Ifges(5,6) + (m(6) * (t102 * t5 + t106 * t6 + t132 * t16 - t133 * t15) + t54 * t132 + t102 * t14 - t55 * t133 + t106 * t13) * pkin(4) + t115 - t165; t99 + (pkin(8) * t87 - t148) * qJD(4) + (m(6) * (t102 * t28 + t106 * t29 + (-t102 * t52 + t106 * t53) * qJD(5)) + (-t102 * t43 + t106 * t42 + (t102 * t75 - t106 * t116) * qJD(5)) * mrSges(6,3)) * pkin(4) + t117; 0.2e1 * t71; t129; t115; t117; t71; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
