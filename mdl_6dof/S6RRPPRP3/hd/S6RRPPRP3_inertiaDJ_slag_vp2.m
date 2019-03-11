% Calculate time derivative of joint inertia matrix for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:46
% EndTime: 2019-03-09 08:33:48
% DurationCPUTime: 1.98s
% Computational Cost: add. (1352->309), mult. (2755->424), div. (0->0), fcn. (1794->4), ass. (0->131)
t162 = Ifges(6,5) + Ifges(7,5);
t159 = Ifges(6,6) + Ifges(7,6);
t161 = Ifges(6,3) + Ifges(7,3);
t138 = pkin(7) - qJ(4);
t93 = cos(qJ(2));
t123 = qJD(5) * t93;
t91 = sin(qJ(2));
t127 = qJD(2) * t91;
t90 = sin(qJ(5));
t92 = cos(qJ(5));
t96 = t90 * t123 + t92 * t127;
t119 = t92 * t123;
t95 = t90 * t127 - t119;
t55 = -t93 * pkin(2) - t91 * qJ(3) - pkin(1);
t41 = t93 * pkin(3) - t55;
t26 = pkin(4) * t91 + pkin(8) * t93 + t41;
t58 = t138 * t91;
t42 = t92 * t58;
t12 = t90 * t26 + t42;
t158 = 2 * m(7);
t157 = -0.2e1 * pkin(1);
t156 = 2 * mrSges(5,1);
t155 = -2 * mrSges(5,3);
t154 = 2 * mrSges(7,3);
t126 = qJD(2) * t93;
t131 = qJ(3) * t126 + t91 * qJD(3);
t36 = pkin(2) * t127 - t131;
t153 = -0.2e1 * t36;
t152 = 0.2e1 * t55;
t149 = pkin(7) * t93;
t85 = t93 * qJ(4);
t59 = -t85 + t149;
t151 = 0.2e1 * t59;
t89 = qJ(3) + pkin(4);
t150 = 0.2e1 * t89;
t94 = -pkin(2) - pkin(3);
t148 = mrSges(6,2) * t90;
t147 = Ifges(6,4) * t90;
t146 = Ifges(6,4) * t92;
t145 = Ifges(7,4) * t90;
t144 = Ifges(7,4) * t92;
t35 = qJD(4) * t93 + t127 * t138;
t143 = t35 * t59;
t142 = t90 * t93;
t21 = t92 * t26;
t141 = t92 * t93;
t140 = mrSges(6,2) + mrSges(7,2);
t139 = mrSges(5,3) - mrSges(4,2);
t22 = mrSges(7,1) * t126 - mrSges(7,3) * t96;
t23 = mrSges(6,1) * t126 - mrSges(6,3) * t96;
t137 = t22 + t23;
t24 = -mrSges(7,2) * t126 - mrSges(7,3) * t95;
t25 = -mrSges(6,2) * t126 - mrSges(6,3) * t95;
t136 = t24 + t25;
t99 = -Ifges(7,2) * t90 + t144;
t31 = Ifges(7,6) * t91 - t93 * t99;
t101 = -Ifges(6,2) * t90 + t146;
t32 = Ifges(6,6) * t91 - t101 * t93;
t135 = t31 + t32;
t103 = Ifges(7,1) * t92 - t145;
t33 = t91 * Ifges(7,5) - t103 * t93;
t105 = Ifges(6,1) * t92 - t147;
t34 = t91 * Ifges(6,5) - t105 * t93;
t134 = t33 + t34;
t49 = -mrSges(7,2) * t91 + mrSges(7,3) * t142;
t50 = -mrSges(6,2) * t91 + mrSges(6,3) * t142;
t133 = t49 + t50;
t51 = mrSges(7,1) * t91 + mrSges(7,3) * t141;
t52 = mrSges(6,1) * t91 + mrSges(6,3) * t141;
t132 = -t51 - t52;
t130 = mrSges(5,1) * t126 + mrSges(5,2) * t127;
t129 = qJ(6) * t93;
t88 = -pkin(8) + t94;
t128 = qJ(6) - t88;
t125 = qJD(3) * t59;
t124 = qJD(5) * t90;
t122 = qJD(6) * t93;
t17 = (pkin(4) * t93 + t88 * t91) * qJD(2) + t131;
t37 = -qJD(4) * t91 + t126 * t138;
t121 = qJD(5) * t21 + t90 * t17 + t92 * t37;
t116 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t100 = Ifges(6,2) * t92 + t147;
t98 = Ifges(7,2) * t92 + t145;
t115 = -t98 / 0.2e1 - t100 / 0.2e1;
t102 = Ifges(7,1) * t90 + t144;
t104 = Ifges(6,1) * t90 + t146;
t114 = -t102 / 0.2e1 - t104 / 0.2e1;
t113 = (m(7) * pkin(5)) + mrSges(7,1);
t112 = m(6) * t88 - mrSges(6,3);
t111 = t159 * t91;
t110 = t92 * t17 - t37 * t90;
t11 = -t58 * t90 + t21;
t109 = qJD(5) * t128;
t108 = -mrSges(6,1) - t113;
t107 = -mrSges(6,1) * t90 - mrSges(6,2) * t92;
t106 = -mrSges(7,1) * t90 - mrSges(7,2) * t92;
t43 = t106 * qJD(5);
t97 = t119 * t159 + t161 * t126 + t162 * t96;
t15 = mrSges(7,1) * t95 + mrSges(7,2) * t96;
t82 = Ifges(6,6) * t124;
t81 = Ifges(7,6) * t124;
t74 = pkin(5) * t92 + t89;
t73 = -pkin(5) * t124 + qJD(3);
t57 = mrSges(6,1) * t92 - t148;
t56 = mrSges(7,1) * t92 - mrSges(7,2) * t90;
t54 = t128 * t92;
t53 = t128 * t90;
t48 = t105 * qJD(5);
t47 = t103 * qJD(5);
t46 = t101 * qJD(5);
t45 = t99 * qJD(5);
t44 = t107 * qJD(5);
t40 = t107 * t93;
t39 = t106 * t93;
t38 = -t85 + (-pkin(5) * t90 + pkin(7)) * t93;
t30 = qJD(6) * t90 + t109 * t92;
t29 = -qJD(6) * t92 + t109 * t90;
t27 = t127 * t94 + t131;
t18 = pkin(5) * t95 - t35;
t16 = mrSges(6,1) * t95 + mrSges(6,2) * t96;
t10 = t104 * t123 + (t93 * Ifges(6,5) + t105 * t91) * qJD(2);
t9 = t102 * t123 + (t93 * Ifges(7,5) + t103 * t91) * qJD(2);
t8 = t100 * t123 + (t93 * Ifges(6,6) + t101 * t91) * qJD(2);
t7 = t98 * t123 + (t93 * Ifges(7,6) + t91 * t99) * qJD(2);
t6 = t129 * t90 + t12;
t5 = pkin(5) * t91 + t129 * t92 + t11;
t4 = -t12 * qJD(5) + t110;
t3 = -t124 * t58 + t121;
t2 = qJ(6) * t119 + (-qJ(6) * t127 - qJD(5) * t58 + t122) * t90 + t121;
t1 = t92 * t122 + (-qJ(6) * t91 * t92 + pkin(5) * t93) * qJD(2) + (-t42 + (-t26 - t129) * t90) * qJD(5) + t110;
t13 = [0.2e1 * t41 * t130 + 0.2e1 * t2 * t49 + 0.2e1 * t3 * t50 + 0.2e1 * t1 * t51 + 0.2e1 * t4 * t52 + t16 * t151 + 0.2e1 * t38 * t15 + 0.2e1 * t18 * t39 - 0.2e1 * t35 * t40 + 0.2e1 * t5 * t22 + 0.2e1 * t11 * t23 + 0.2e1 * t6 * t24 + 0.2e1 * t12 * t25 + m(4) * t36 * t152 + (t1 * t5 + t18 * t38 + t2 * t6) * t158 + 0.2e1 * m(6) * (t11 * t4 + t12 * t3 - t143) + 0.2e1 * m(5) * (t27 * t41 + t37 * t58 - t143) + (t27 * t156 + mrSges(4,3) * t153 + t37 * t155 + (mrSges(3,1) * t157 + mrSges(4,1) * t152 + mrSges(5,3) * t151 + 0.2e1 * (-Ifges(3,4) - Ifges(5,4) + Ifges(4,5)) * t91 + t134 * t92 + (-t111 - t135) * t90) * qJD(2) + t97) * t91 + (mrSges(4,1) * t153 - 0.2e1 * t27 * mrSges(5,2) + 0.2e1 * t35 * mrSges(5,3) + (-t10 - t9) * t92 + (t7 + t8) * t90 + (t134 * t90 + t135 * t92) * qJD(5) + (-0.2e1 * t55 * mrSges(4,3) + mrSges(3,2) * t157 + t58 * t155 + (t159 * t90 - t162 * t92 + 0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(5,4) - 0.2e1 * Ifges(4,5)) * t93 + (-(2 * Ifges(4,3)) - (2 * Ifges(3,2)) - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + (2 * Ifges(4,1)) + (2 * Ifges(3,1)) + t161) * t91) * qJD(2)) * t93; t37 * mrSges(5,2) + t74 * t15 + t89 * t16 + t18 * t56 + t53 * t22 - t54 * t24 + t29 * t49 + t30 * t51 + t38 * t43 + t73 * t39 + t59 * t44 + (t81 / 0.2e1 + t82 / 0.2e1) * t91 + (-t57 - mrSges(5,1)) * t35 + m(5) * (-qJ(3) * t35 + t37 * t94 + t125) + m(7) * (t1 * t53 + t18 * t74 - t2 * t54 + t29 * t6 + t30 * t5 + t38 * t73) + m(6) * (-t35 * t89 + t125) + (t88 * t25 - t2 * mrSges(7,3) - t7 / 0.2e1 - t8 / 0.2e1 + (t47 / 0.2e1 + t48 / 0.2e1) * t93 + t112 * t3 + (-t33 / 0.2e1 - t34 / 0.2e1 - t88 * t52 + t5 * mrSges(7,3) + t115 * t93 + t116 * t91 - t112 * t11) * qJD(5)) * t92 + (-t88 * t23 + t1 * mrSges(7,3) - t9 / 0.2e1 - t10 / 0.2e1 + (-t45 / 0.2e1 - t46 / 0.2e1) * t93 - t112 * t4 + (t31 / 0.2e1 + t32 / 0.2e1 - t88 * t50 + t6 * mrSges(7,3) + t114 * t93 - t112 * t12) * qJD(5)) * t90 + ((-pkin(2) * mrSges(4,2) - t94 * mrSges(5,3) + Ifges(4,4) + Ifges(3,5) + Ifges(5,6) + (-Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t92 + t116 * t90 + (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1)) * pkin(7)) * t93 + (-Ifges(5,5) + Ifges(4,6) - Ifges(3,6) + t114 * t92 - t115 * t90 + t139 * qJ(3) + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t91) * qJD(2) + (m(4) * t149 - t139 * t93 + t40) * qJD(3); (-t29 * t54 + t30 * t53 + t73 * t74) * t158 + 0.2e1 * t73 * t56 + 0.2e1 * t74 * t43 + t44 * t150 + (-0.2e1 * mrSges(7,3) * t29 + t45 + t46) * t92 + (t30 * t154 + t47 + t48) * t90 + ((t154 * t53 + t102 + t104) * t92 + (-t154 * t54 - t100 - t98) * t90) * qJD(5) + (m(6) * t150 + t156 + 0.2e1 * mrSges(4,3) + 0.2e1 * t57 + 0.2e1 * (m(4) + m(5)) * qJ(3)) * qJD(3); t136 * t92 - t137 * t90 + (m(4) * pkin(7) - t139) * t126 + (t132 * t92 - t133 * t90) * qJD(5) + m(7) * (-t1 * t90 + t2 * t92 + (-t5 * t92 - t6 * t90) * qJD(5)) + m(6) * (t3 * t92 - t4 * t90 + (-t11 * t92 - t12 * t90) * qJD(5)) + m(5) * t37; m(7) * (t29 * t92 - t30 * t90 + (-t53 * t92 + t54 * t90) * qJD(5)); 0; t137 * t92 + t136 * t90 + (t132 * t90 + t133 * t92) * qJD(5) + m(7) * (t1 * t92 + t2 * t90 + (-t5 * t90 + t6 * t92) * qJD(5)) + m(6) * (t3 * t90 + t4 * t92 + (-t11 * t90 + t12 * t92) * qJD(5)) + m(5) * t27 + t130; m(7) * (t29 * t90 + t30 * t92 + (-t53 * t90 - t54 * t92) * qJD(5)); 0; 0; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 - t90 * qJD(2) * t111 + (m(7) * t1 + t22) * pkin(5) + t97; -mrSges(7,2) * t29 + t81 + t82 + t113 * t30 + (t88 * t148 + (-mrSges(6,1) * t88 + (mrSges(7,3) * pkin(5)) - t162) * t92) * qJD(5); (t108 * t92 + t140 * t90) * qJD(5); (t108 * t90 - t140 * t92) * qJD(5); 0; m(7) * t18 + t15; m(7) * t73 + t43; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
