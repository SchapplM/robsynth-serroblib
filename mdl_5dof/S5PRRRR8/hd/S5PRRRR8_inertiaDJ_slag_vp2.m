% Calculate time derivative of joint inertia matrix for
% S5PRRRR8
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:52
% EndTime: 2019-12-05 17:15:05
% DurationCPUTime: 1.90s
% Computational Cost: add. (2025->244), mult. (5151->388), div. (0->0), fcn. (4725->10), ass. (0->124)
t89 = sin(qJ(5));
t93 = cos(qJ(5));
t94 = cos(qJ(4));
t165 = (t89 ^ 2 + t93 ^ 2) * t94;
t75 = -mrSges(6,1) * t93 + mrSges(6,2) * t89;
t166 = t75 - mrSges(5,1);
t127 = qJD(5) * t93;
t164 = qJD(3) + qJD(4);
t90 = sin(qJ(4));
t91 = sin(qJ(3));
t95 = cos(qJ(3));
t67 = t90 * t91 - t94 * t95;
t48 = t164 * t67;
t68 = t90 * t95 + t91 * t94;
t104 = t68 * t127 - t48 * t89;
t163 = 2 * m(6);
t162 = -2 * mrSges(5,3);
t157 = -pkin(8) - pkin(7);
t118 = qJD(3) * t157;
t114 = t95 * t118;
t124 = t157 * t91;
t78 = t157 * t95;
t59 = t124 * t90 - t94 * t78;
t74 = t91 * t118;
t25 = qJD(4) * t59 - t114 * t94 + t90 * t74;
t161 = 0.2e1 * t25;
t58 = -t124 * t94 - t78 * t90;
t160 = 0.2e1 * t58;
t159 = m(5) / 0.2e1;
t125 = pkin(3) * qJD(3) * t91;
t49 = t164 * t68;
t18 = pkin(4) * t49 + pkin(9) * t48 + t125;
t82 = -pkin(3) * t95 - pkin(2);
t41 = pkin(4) * t67 - pkin(9) * t68 + t82;
t21 = t41 * t89 + t59 * t93;
t24 = -qJD(4) * t58 + t114 * t90 + t94 * t74;
t3 = -qJD(5) * t21 + t18 * t93 - t24 * t89;
t156 = t3 * t89;
t92 = sin(qJ(2));
t130 = qJD(2) * t92;
t87 = sin(pkin(5));
t123 = t87 * t130;
t139 = t87 * t92;
t88 = cos(pkin(5));
t61 = -t139 * t91 + t88 * t95;
t62 = t139 * t95 + t88 * t91;
t106 = t61 * t94 - t62 * t90;
t96 = cos(qJ(2));
t129 = qJD(2) * t96;
t122 = t87 * t129;
t53 = -qJD(3) * t62 - t122 * t91;
t54 = qJD(3) * t61 + t122 * t95;
t13 = qJD(4) * t106 + t53 * t90 + t54 * t94;
t138 = t87 * t96;
t37 = t61 * t90 + t62 * t94;
t28 = -t138 * t93 - t37 * t89;
t5 = qJD(5) * t28 + t123 * t89 + t93 * t13;
t155 = t5 * t93;
t105 = t138 * t89 - t37 * t93;
t6 = qJD(5) * t105 + t123 * t93 - t89 * t13;
t154 = t6 * t89;
t153 = mrSges(6,3) * t93;
t152 = Ifges(6,4) * t89;
t151 = Ifges(6,4) * t93;
t150 = Ifges(6,6) * t89;
t149 = t25 * t58;
t14 = qJD(4) * t37 - t53 * t94 + t54 * t90;
t148 = t106 * t14;
t147 = t106 * t90;
t145 = t48 * t93;
t144 = t58 * t90;
t143 = t68 * t89;
t142 = t68 * t93;
t137 = t89 * t94;
t136 = t90 * mrSges(5,1);
t135 = t90 * t75;
t134 = t93 * t94;
t133 = t94 * mrSges(5,2);
t132 = -Ifges(6,5) * t145 + Ifges(6,3) * t49;
t131 = pkin(3) * qJD(4);
t128 = qJD(5) * t89;
t126 = 0.2e1 * t91;
t121 = t68 * t128;
t103 = t121 + t145;
t15 = mrSges(6,1) * t104 - mrSges(6,2) * t103;
t119 = m(6) * t25 + t15;
t116 = -(2 * Ifges(5,4)) - t150;
t115 = t87 ^ 2 * t92 * t129;
t113 = mrSges(6,3) * t165;
t112 = -mrSges(4,1) * t95 + mrSges(4,2) * t91;
t111 = mrSges(6,1) * t89 + mrSges(6,2) * t93;
t110 = Ifges(6,1) * t93 - t152;
t109 = -Ifges(6,2) * t89 + t151;
t108 = Ifges(6,5) * t89 + Ifges(6,6) * t93;
t107 = -t106 * t25 + t58 * t14;
t20 = t41 * t93 - t59 * t89;
t72 = t109 * qJD(5);
t73 = t110 * qJD(5);
t76 = Ifges(6,2) * t93 + t152;
t77 = Ifges(6,1) * t89 + t151;
t102 = t127 * t77 - t128 * t76 + t72 * t93 + t73 * t89;
t101 = -t154 + (t105 * t89 - t28 * t93) * qJD(5);
t100 = -t53 * t91 + t54 * t95 + (-t61 * t95 - t62 * t91) * qJD(3);
t70 = t111 * qJD(5);
t99 = -t13 * mrSges(5,2) + t101 * mrSges(6,3) - t106 * t70 + t14 * t166 + t5 * t153;
t16 = mrSges(6,1) * t49 + mrSges(6,3) * t103;
t17 = -mrSges(6,2) * t49 - mrSges(6,3) * t104;
t2 = qJD(5) * t20 + t18 * t89 + t24 * t93;
t43 = -mrSges(6,2) * t67 - mrSges(6,3) * t143;
t44 = mrSges(6,1) * t67 - mrSges(6,3) * t142;
t98 = m(6) * (-t127 * t20 - t128 * t21 + t2 * t93 - t156) + t93 * t17 - t89 * t16 - t44 * t127 - t43 * t128;
t10 = -Ifges(6,1) * t103 - Ifges(6,4) * t104 + t49 * Ifges(6,5);
t31 = t67 * Ifges(6,6) + t109 * t68;
t32 = Ifges(6,5) * t67 + t110 * t68;
t83 = Ifges(6,5) * t127;
t9 = -Ifges(6,4) * t103 - Ifges(6,2) * t104 + t49 * Ifges(6,6);
t97 = t2 * t153 + t32 * t127 / 0.2e1 - t77 * t145 / 0.2e1 + t58 * t70 - Ifges(5,5) * t48 - t72 * t143 / 0.2e1 + t73 * t142 / 0.2e1 + t67 * (-Ifges(6,6) * t128 + t83) / 0.2e1 + t89 * t10 / 0.2e1 + t93 * t9 / 0.2e1 + (-t156 + (-t20 * t93 - t21 * t89) * qJD(5)) * mrSges(6,3) - t24 * mrSges(5,2) + (t108 / 0.2e1 - Ifges(5,6)) * t49 + t166 * t25 - t104 * t76 / 0.2e1 - (t68 * t77 + t31) * t128 / 0.2e1;
t81 = -pkin(3) * t94 - pkin(4);
t80 = pkin(3) * t90 + pkin(9);
t71 = (mrSges(4,1) * t91 + mrSges(4,2) * t95) * qJD(3);
t52 = mrSges(5,1) * t67 + mrSges(5,2) * t68;
t39 = t111 * t68;
t22 = mrSges(5,1) * t49 - mrSges(5,2) * t48;
t1 = [0.2e1 * m(6) * (-t105 * t5 + t28 * t6 - t148) + 0.2e1 * m(5) * (t13 * t37 - t115 - t148) + 0.2e1 * m(4) * (t53 * t61 + t54 * t62 - t115); t14 * t39 - t106 * t15 + t28 * t16 - t105 * t17 + t5 * t43 + t6 * t44 + ((-t22 - t71) * t96 + (-t96 * mrSges(3,2) + (-mrSges(3,1) + t112 + t52) * t92) * qJD(2)) * t87 + m(6) * (-t105 * t2 + t20 * t6 + t21 * t5 + t28 * t3 + t107) + m(5) * (t59 * t13 + t24 * t37 + (-t125 * t96 + t130 * t82) * t87 + t107) + (t106 * t48 - t13 * t67 + t14 * t68 - t37 * t49) * mrSges(5,3) + t100 * mrSges(4,3) + (-pkin(2) * t123 + pkin(7) * t100) * m(4); t59 * t49 * t162 - 0.2e1 * pkin(2) * t71 + t15 * t160 + 0.2e1 * t20 * t16 + 0.2e1 * t21 * t17 + 0.2e1 * t2 * t43 + 0.2e1 * t82 * t22 + t39 * t161 + 0.2e1 * t3 * t44 + 0.2e1 * m(5) * (t125 * t82 + t24 * t59 + t149) + (t2 * t21 + t20 * t3 + t149) * t163 - (mrSges(5,3) * t160 - t31 * t89 + t32 * t93) * t48 + (t24 * t162 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t49 - t116 * t48 + t132) * t67 + (mrSges(5,3) * t161 - 0.2e1 * Ifges(5,1) * t48 + t93 * t10 - t89 * t9 + (Ifges(6,5) * t93 + t116) * t49 + (-t108 * t67 - t93 * t31 - t89 * t32) * qJD(5)) * t68 + ((-Ifges(4,4) * t91 + pkin(3) * t52) * t126 + (0.2e1 * Ifges(4,4) * t95 + (Ifges(4,1) - Ifges(4,2)) * t126) * t95) * qJD(3); m(6) * t14 * t81 - t54 * mrSges(4,2) + t53 * mrSges(4,1) + 0.2e1 * ((t13 * t90 - t14 * t94) * t159 + (m(6) * (-t105 * t134 - t28 * t137 - t147) / 0.2e1 + (t37 * t94 - t147) * t159) * qJD(4)) * pkin(3) + t99 + m(6) * (t105 * t128 - t127 * t28 - t154 + t155) * t80; (Ifges(4,5) * t95 - Ifges(4,6) * t91 + t112 * pkin(7)) * qJD(3) + t97 + t119 * t81 + t98 * t80 + (m(5) * (t24 * t90 - t25 * t94) + (t48 * t94 - t49 * t90) * mrSges(5,3) + ((t68 * mrSges(5,3) + t39) * t90 + (-t67 * mrSges(5,3) + t93 * t43 - t89 * t44) * t94 + m(6) * (t134 * t21 - t137 * t20 + t144) + m(5) * (t59 * t94 + t144)) * qJD(4)) * pkin(3); 0.2e1 * t81 * t70 + (-0.2e1 * t133 - 0.2e1 * t136 + 0.2e1 * t135 + (t165 * t80 + t81 * t90) * t163 + 0.2e1 * t113) * t131 + t102; m(6) * (-pkin(4) * t14 + (t101 + t155) * pkin(9)) + t99; -pkin(4) * t119 + pkin(9) * t98 + t97; (-pkin(4) + t81) * t70 + (-t133 - t136 + m(6) * (-pkin(4) * t90 + pkin(9) * t165) + t135 + t113) * t131 + t102; -0.2e1 * pkin(4) * t70 + t102; mrSges(6,1) * t6 - mrSges(6,2) * t5; mrSges(6,1) * t3 - mrSges(6,2) * t2 - Ifges(6,5) * t121 - Ifges(6,6) * t104 + t132; t83 - t111 * t94 * t131 + (t75 * t80 - t150) * qJD(5); t83 + (pkin(9) * t75 - t150) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
