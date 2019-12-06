% Calculate time derivative of joint inertia matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:09
% EndTime: 2019-12-05 16:54:14
% DurationCPUTime: 1.65s
% Computational Cost: add. (961->274), mult. (2718->403), div. (0->0), fcn. (2097->8), ass. (0->126)
t126 = Ifges(5,5) + Ifges(6,5);
t157 = -Ifges(5,3) - Ifges(6,3);
t85 = cos(qJ(3));
t115 = qJD(3) * t85;
t84 = cos(qJ(4));
t104 = t84 * t115;
t82 = sin(qJ(3));
t113 = qJD(4) * t82;
t81 = sin(qJ(4));
t89 = -t81 * t113 + t104;
t112 = qJD(4) * t84;
t88 = t82 * t112 + t81 * t115;
t125 = Ifges(5,6) + Ifges(6,6);
t156 = t125 * t84 + t126 * t81;
t135 = Ifges(6,4) * t84;
t93 = -Ifges(6,2) * t81 + t135;
t137 = Ifges(5,4) * t84;
t94 = -Ifges(5,2) * t81 + t137;
t155 = (t93 + t94) * qJD(4);
t136 = Ifges(6,4) * t81;
t95 = Ifges(6,1) * t84 - t136;
t138 = Ifges(5,4) * t81;
t96 = Ifges(5,1) * t84 - t138;
t154 = (t95 + t96) * qJD(4);
t55 = -pkin(3) * t85 - pkin(8) * t82 - pkin(2);
t139 = pkin(7) * t85;
t70 = t84 * t139;
t26 = t81 * t55 + t70;
t153 = -m(5) * pkin(8) - mrSges(5,3);
t152 = -m(5) * pkin(3) - mrSges(5,1) * t84 + mrSges(5,2) * t81 - mrSges(4,1);
t57 = -mrSges(6,1) * t84 + mrSges(6,2) * t81;
t72 = -pkin(4) * t84 - pkin(3);
t151 = m(6) * t72 + t57;
t150 = 0.2e1 * m(5);
t149 = 0.2e1 * m(6);
t148 = 0.2e1 * pkin(7);
t147 = -2 * mrSges(6,3);
t146 = m(5) / 0.2e1;
t143 = m(6) * pkin(4);
t140 = pkin(7) * t81;
t86 = cos(qJ(2));
t117 = qJD(2) * t86;
t79 = sin(pkin(5));
t108 = t79 * t117;
t83 = sin(qJ(2));
t131 = t79 * t83;
t80 = cos(pkin(5));
t34 = t131 * t85 + t80 * t82;
t15 = qJD(3) * t34 + t108 * t82;
t134 = t15 * t82;
t33 = t131 * t82 - t80 * t85;
t16 = -qJD(3) * t33 + t108 * t85;
t133 = t16 * t85;
t132 = t33 * t15;
t130 = t79 * t86;
t129 = t81 * t82;
t128 = t82 * t84;
t124 = -qJ(5) - pkin(8);
t29 = -t85 * Ifges(6,5) + t82 * t95;
t30 = -t85 * Ifges(5,5) + t82 * t96;
t123 = t29 + t30;
t53 = (pkin(3) * t82 - pkin(8) * t85) * qJD(3);
t122 = t112 * t55 + t53 * t81;
t116 = qJD(3) * t82;
t121 = t116 * t140 + t53 * t84;
t114 = qJD(4) * t81;
t42 = mrSges(6,1) * t114 + mrSges(6,2) * t112;
t119 = qJ(5) * t82;
t118 = qJ(5) * t84;
t111 = qJD(5) * t84;
t110 = pkin(4) * t114;
t109 = qJD(2) * t131;
t60 = Ifges(6,2) * t84 + t136;
t61 = Ifges(5,2) * t84 + t138;
t103 = -t60 / 0.2e1 - t61 / 0.2e1;
t62 = Ifges(6,1) * t81 + t135;
t63 = Ifges(5,1) * t81 + t137;
t102 = -t63 / 0.2e1 - t62 / 0.2e1;
t101 = mrSges(6,1) + t143;
t100 = t125 * t81;
t99 = qJD(4) * t124;
t98 = -t104 * t126 + t116 * t157;
t97 = mrSges(5,1) * t81 + mrSges(5,2) * t84;
t27 = -t85 * Ifges(6,6) + t82 * t93;
t28 = -t85 * Ifges(5,6) + t82 * t94;
t92 = t125 * t85 - t27 - t28;
t17 = -t130 * t84 - t34 * t81;
t91 = t130 * t81 - t34 * t84;
t90 = t115 * t33 + t134;
t12 = mrSges(6,1) * t88 + mrSges(6,2) * t89;
t78 = Ifges(5,5) * t112;
t77 = Ifges(6,5) * t112;
t59 = t124 * t84;
t56 = t124 * t81;
t54 = (pkin(4) * t81 + pkin(7)) * t82;
t52 = -mrSges(5,1) * t85 - mrSges(5,3) * t128;
t51 = -mrSges(6,1) * t85 - mrSges(6,3) * t128;
t50 = mrSges(5,2) * t85 - mrSges(5,3) * t129;
t49 = mrSges(6,2) * t85 - mrSges(6,3) * t129;
t44 = (mrSges(4,1) * t82 + mrSges(4,2) * t85) * qJD(3);
t43 = t97 * qJD(4);
t41 = t84 * t55;
t38 = t97 * t82;
t37 = (mrSges(6,1) * t81 + mrSges(6,2) * t84) * t82;
t32 = -qJD(5) * t81 + t84 * t99;
t31 = t81 * t99 + t111;
t25 = -t139 * t81 + t41;
t24 = pkin(4) * t88 + pkin(7) * t115;
t23 = -mrSges(5,2) * t116 - mrSges(5,3) * t88;
t22 = -mrSges(6,2) * t116 - mrSges(6,3) * t88;
t21 = mrSges(5,1) * t116 - mrSges(5,3) * t89;
t20 = mrSges(6,1) * t116 - mrSges(6,3) * t89;
t19 = -t119 * t81 + t26;
t14 = -t82 * t118 + t41 + (-pkin(4) - t140) * t85;
t13 = mrSges(5,1) * t88 + mrSges(5,2) * t89;
t11 = -t63 * t113 + (Ifges(5,5) * t82 + t85 * t96) * qJD(3);
t10 = -t62 * t113 + (Ifges(6,5) * t82 + t85 * t95) * qJD(3);
t9 = -t61 * t113 + (Ifges(5,6) * t82 + t85 * t94) * qJD(3);
t8 = -t60 * t113 + (Ifges(6,6) * t82 + t85 * t93) * qJD(3);
t7 = -qJD(4) * t26 + t121;
t6 = (-t114 * t85 - t116 * t84) * pkin(7) + t122;
t5 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t128 + (-qJD(5) * t82 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t85) * t81 + t122;
t4 = qJD(4) * t17 + t109 * t81 + t16 * t84;
t3 = qJD(4) * t91 + t109 * t84 - t16 * t81;
t2 = -t82 * t111 + (pkin(4) * t82 - t118 * t85) * qJD(3) + (-t70 + (-t55 + t119) * t81) * qJD(4) + t121;
t1 = [0.2e1 * m(4) * (-t117 * t79 ^ 2 * t83 + t16 * t34 + t132) + 0.4e1 * (t146 + m(6) / 0.2e1) * (t17 * t3 - t4 * t91 + t132); (t50 + t49) * t4 + (t13 + t12) * t33 + (t51 + t52) * t3 - (t23 + t22) * t91 + (t20 + t21) * t17 + (t38 + t37) * t15 + (-t86 * t44 + (-t86 * mrSges(3,2) + (-mrSges(4,1) * t85 + mrSges(4,2) * t82 - mrSges(3,1)) * t83) * qJD(2)) * t79 + m(5) * (t17 * t7 + t25 * t3 + t26 * t4 - t6 * t91) + m(6) * (t14 * t3 + t15 * t54 + t17 * t2 + t19 * t4 + t24 * t33 - t5 * t91) - m(4) * pkin(2) * t109 + (t90 * t146 + m(4) * (-t34 * t116 + t133 + t90) / 0.2e1) * t148 + (t134 + t133 + (t33 * t85 - t34 * t82) * qJD(3)) * mrSges(4,3); -0.2e1 * pkin(2) * t44 + 0.2e1 * t54 * t12 + 0.2e1 * t14 * t20 + 0.2e1 * t19 * t22 + 0.2e1 * t2 * t51 + 0.2e1 * t25 * t21 + 0.2e1 * t26 * t23 + 0.2e1 * t24 * t37 + 0.2e1 * t5 * t49 + 0.2e1 * t6 * t50 + 0.2e1 * t7 * t52 + (t25 * t7 + t26 * t6) * t150 + (t14 * t2 + t19 * t5 + t24 * t54) * t149 + ((0.2e1 * Ifges(4,4) * t85 + t123 * t84 + t148 * t38 + t81 * t92) * qJD(3) + t98) * t85 + (t13 * t148 + (t10 + t11) * t84 + (-t8 - t9) * t81 + (t92 * t84 + (t126 * t85 - t123) * t81) * qJD(4) + ((t126 * t84 - 0.2e1 * Ifges(4,4) - t100) * t82 + (pkin(7) ^ 2 * t150 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t157) * t85) * qJD(3)) * t82; -t16 * mrSges(4,2) + (t42 + t43) * t33 + m(6) * (t110 * t33 + t17 * t32 + t3 * t56 - t31 * t91 - t4 * t59) + (t151 + t152) * t15 + (mrSges(6,3) - t153) * (-t3 * t81 + t4 * t84 + (-t17 * t84 + t81 * t91) * qJD(4)); m(6) * (t14 * t32 + t19 * t31 + t2 * t56 + t24 * t72 - t5 * t59) + t72 * t12 + t32 * t51 + t54 * t42 + t56 * t20 + t24 * t57 - t59 * t22 + t31 * t49 - pkin(3) * t13 + (-t77 / 0.2e1 - t78 / 0.2e1 + (pkin(7) * t152 + Ifges(4,5)) * qJD(3)) * t85 + (-t7 * mrSges(5,3) - t2 * mrSges(6,3) + t10 / 0.2e1 + t11 / 0.2e1 + t103 * t115 + (-m(5) * t7 - t21) * pkin(8) + (-t27 / 0.2e1 - t28 / 0.2e1 + pkin(4) * t37 - pkin(8) * t50 - t19 * mrSges(6,3) + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t85 + t54 * t143 + t153 * t26) * qJD(4)) * t81 + (t5 * mrSges(6,3) + t6 * mrSges(5,3) + t8 / 0.2e1 + t9 / 0.2e1 - t102 * t115 + (t29 / 0.2e1 + t30 / 0.2e1 - t14 * mrSges(6,3) - t25 * mrSges(5,3)) * qJD(4) + (t23 + m(5) * (-qJD(4) * t25 + t6) - qJD(4) * t52) * pkin(8)) * t84 + (pkin(7) * t43 + (t102 * t81 + t103 * t84) * qJD(4) - t155 * t81 / 0.2e1 + t154 * t84 / 0.2e1 + (-Ifges(4,6) + pkin(7) * mrSges(4,2) + t156 / 0.2e1) * qJD(3)) * t82; 0.2e1 * t72 * t42 + (-t31 * t59 + t32 * t56) * t149 - 0.2e1 * pkin(3) * t43 + (t32 * t147 + (0.2e1 * pkin(4) * t151 - t59 * t147 - t60 - t61) * qJD(4) + t154) * t81 + (0.2e1 * t31 * mrSges(6,3) + (t147 * t56 + t62 + t63) * qJD(4) + t155) * t84; (-mrSges(5,2) - mrSges(6,2)) * t4 + (mrSges(5,1) + t101) * t3; mrSges(5,1) * t7 + mrSges(6,1) * t2 - mrSges(5,2) * t6 - mrSges(6,2) * t5 - t100 * t115 + (m(6) * t2 + t20) * pkin(4) - t156 * t113 - t98; -mrSges(6,2) * t31 + t77 + t78 + t101 * t32 + ((-mrSges(5,1) * pkin(8) - mrSges(6,3) * pkin(4)) * t84 + (mrSges(5,2) * pkin(8) - t125) * t81) * qJD(4); 0; m(6) * t15; m(6) * t24 + t12; m(6) * t110 + t42; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
