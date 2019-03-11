% Calculate time derivative of joint inertia matrix for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:42
% EndTime: 2019-03-09 02:07:46
% DurationCPUTime: 2.02s
% Computational Cost: add. (1077->298), mult. (2349->427), div. (0->0), fcn. (1525->4), ass. (0->131)
t158 = Ifges(6,6) + Ifges(7,6);
t157 = Ifges(6,3) + Ifges(7,3);
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t49 = -mrSges(7,1) * t81 + mrSges(7,2) * t79;
t63 = -pkin(5) * t81 - pkin(4);
t156 = m(7) * t63 + t49;
t80 = sin(qJ(4));
t117 = qJD(4) * t80;
t105 = t79 * t117;
t82 = cos(qJ(4));
t111 = qJD(5) * t82;
t106 = t81 * t111;
t84 = t105 - t106;
t115 = qJD(4) * t82;
t77 = qJ(2) - pkin(7);
t109 = t77 * t115;
t155 = qJD(2) * t80 + t109;
t153 = Ifges(6,5) + Ifges(7,5);
t122 = t79 ^ 2 + t81 ^ 2;
t108 = t79 * t111;
t116 = qJD(4) * t81;
t85 = t116 * t80 + t108;
t10 = -t84 * mrSges(7,1) - mrSges(7,2) * t85;
t13 = -pkin(5) * t84 - qJD(2) * t82 + t117 * t77;
t152 = m(7) * t13 + t10;
t129 = -qJ(6) - pkin(8);
t48 = t129 * t79;
t51 = t129 * t81;
t151 = t48 * t79 + t51 * t81;
t133 = t82 * mrSges(7,3);
t44 = t80 * mrSges(7,1) - t133 * t81;
t134 = t82 * mrSges(6,3);
t45 = t80 * mrSges(6,1) - t134 * t81;
t124 = t44 + t45;
t42 = -t80 * mrSges(7,2) - t133 * t79;
t43 = -t80 * mrSges(6,2) - t134 * t79;
t125 = t42 + t43;
t150 = t124 * t79 - t125 * t81;
t149 = 0.2e1 * m(7);
t148 = -0.2e1 * mrSges(7,3);
t147 = 2 * qJD(3);
t146 = m(7) * pkin(5);
t144 = pkin(8) * t82;
t143 = t80 * pkin(4);
t142 = Ifges(6,4) * t79;
t141 = Ifges(6,4) * t81;
t140 = Ifges(7,4) * t79;
t139 = Ifges(7,4) * t81;
t136 = t79 * t80;
t135 = t80 * t81;
t132 = mrSges(6,2) + mrSges(7,2);
t78 = pkin(1) + qJ(3);
t17 = mrSges(7,1) * t115 + mrSges(7,3) * t85;
t18 = mrSges(6,1) * t115 + mrSges(6,3) * t85;
t128 = -t17 - t18;
t19 = -mrSges(7,2) * t115 + mrSges(7,3) * t84;
t20 = -mrSges(6,2) * t115 + mrSges(6,3) * t84;
t127 = t19 + t20;
t89 = -Ifges(7,2) * t79 + t139;
t21 = Ifges(7,6) * t80 + t82 * t89;
t90 = -Ifges(6,2) * t79 + t141;
t22 = Ifges(6,6) * t80 + t82 * t90;
t126 = t21 + t22;
t123 = -mrSges(6,1) * t81 + mrSges(6,2) * t79 - mrSges(5,1);
t46 = t143 + t78 - t144;
t16 = t77 * t135 + t79 * t46;
t112 = qJD(5) * t81;
t114 = qJD(5) * t79;
t36 = mrSges(7,1) * t114 + mrSges(7,2) * t112;
t37 = mrSges(6,1) * t114 + mrSges(6,2) * t112;
t121 = qJ(6) * t82;
t74 = t80 ^ 2;
t120 = qJD(2) * t74;
t118 = qJD(3) * t78;
t113 = qJD(5) * t80;
t76 = t82 ^ 2;
t72 = t76 * qJD(2);
t110 = t114 * t146 + t36;
t107 = t80 * t112;
t104 = t79 * t113;
t103 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t52 = Ifges(7,2) * t81 + t140;
t53 = Ifges(6,2) * t81 + t142;
t102 = t52 / 0.2e1 + t53 / 0.2e1;
t54 = Ifges(7,1) * t79 + t139;
t55 = Ifges(6,1) * t79 + t141;
t101 = -t54 / 0.2e1 - t55 / 0.2e1;
t100 = mrSges(7,1) + t146;
t99 = t153 * t80;
t98 = qJD(5) * t129;
t35 = qJD(3) + (pkin(4) * t82 + pkin(8) * t80) * qJD(4);
t97 = t46 * t112 + t155 * t81 + t79 * t35;
t96 = t158 * t105 + t157 * t115;
t95 = m(6) * pkin(4) - t123;
t94 = -mrSges(6,1) - t100;
t12 = -t121 * t79 + t16;
t33 = t81 * t46;
t9 = -t81 * t121 + t33 + (-t77 * t79 + pkin(5)) * t80;
t93 = t12 * t81 - t79 * t9;
t92 = Ifges(6,1) * t81 - t142;
t91 = Ifges(7,1) * t81 - t140;
t15 = -t136 * t77 + t33;
t88 = t15 * t79 - t16 * t81;
t23 = Ifges(7,5) * t80 + t82 * t91;
t24 = Ifges(6,5) * t80 + t82 * t92;
t87 = -t23 - t24 - t99;
t83 = -qJD(6) * t82 + (qJ(6) * qJD(4) - qJD(5) * t77) * t80;
t71 = Ifges(6,5) * t112;
t70 = Ifges(7,5) * t112;
t60 = t77 * t72;
t41 = t92 * qJD(5);
t40 = t91 * qJD(5);
t39 = t90 * qJD(5);
t38 = t89 * qJD(5);
t34 = (pkin(5) * t79 - t77) * t82;
t31 = (mrSges(6,1) * t79 + mrSges(6,2) * t81) * t82;
t30 = (mrSges(7,1) * t79 + mrSges(7,2) * t81) * t82;
t28 = t81 * t35;
t26 = -qJD(6) * t79 + t81 * t98;
t25 = qJD(6) * t81 + t79 * t98;
t11 = -mrSges(6,1) * t84 - mrSges(6,2) * t85;
t8 = -t55 * t111 + (Ifges(6,5) * t82 - t80 * t92) * qJD(4);
t7 = -t54 * t111 + (Ifges(7,5) * t82 - t80 * t91) * qJD(4);
t6 = -t53 * t111 + (Ifges(6,6) * t82 - t80 * t90) * qJD(4);
t5 = -t52 * t111 + (Ifges(7,6) * t82 - t80 * t89) * qJD(4);
t4 = -t77 * t107 + t28 + (-qJD(5) * t46 - t155) * t79;
t3 = -t104 * t77 + t97;
t2 = -qJ(6) * t106 + t79 * t83 + t97;
t1 = pkin(5) * t115 + t28 + t83 * t81 + ((-t46 + t121) * qJD(5) - t155) * t79;
t14 = [(mrSges(4,3) * t147) + 0.2e1 * t1 * t44 + 0.2e1 * t34 * t10 + 0.2e1 * t12 * t19 + 0.2e1 * t13 * t30 + 0.2e1 * t15 * t18 + 0.2e1 * t16 * t20 + 0.2e1 * t9 * t17 + 0.2e1 * t2 * t42 + 0.2e1 * t3 * t43 + 0.2e1 * t4 * t45 + 0.2e1 * (m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3) + (-t74 - t76) * mrSges(5,3)) * qJD(2) + 0.2e1 * m(6) * (t15 * t4 + t16 * t3 + t60) + 0.2e1 * m(5) * (t120 * t77 + t118 + t60) + 0.2e1 * m(4) * (qJ(2) * qJD(2) + t118) + (t1 * t9 + t12 * t2 + t13 * t34) * t149 + (mrSges(5,1) * t147 + (-0.2e1 * t78 * mrSges(5,2) + 0.2e1 * Ifges(5,4) * t80 + t126 * t79 + 0.2e1 * t77 * t31 + t81 * t87) * qJD(4) + t96) * t80 + (mrSges(5,2) * t147 - 0.2e1 * qJD(2) * t31 - 0.2e1 * t77 * t11 + (t7 + t8) * t81 + (-t5 - t6) * t79 + ((-t158 * t80 - t126) * t81 + t87 * t79) * qJD(5) + (0.2e1 * t78 * mrSges(5,1) + (t153 * t81 - t158 * t79 - 0.2e1 * Ifges(5,4)) * t82 + (-0.2e1 * m(6) * t77 ^ 2 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + t157) * t80) * qJD(4)) * t82; t128 * t81 - t127 * t79 + (-t82 * mrSges(5,1) + t80 * mrSges(5,2)) * qJD(4) + (-m(5) - m(4)) * qJD(3) + t150 * qJD(5) + m(7) * (-qJD(5) * t93 - t1 * t81 - t2 * t79) + m(6) * (qJD(5) * t88 - t3 * t79 - t4 * t81); 0; m(4) * qJD(2) + m(6) * t72 + m(5) * (t72 + t120) + (-t11 + (-m(6) * t88 + m(7) * t93 - t150) * qJD(4) - t152) * t82 + (t127 * t81 + t128 * t79 + (t30 + t31) * qJD(4) + (-t124 * t81 - t125 * t79) * qJD(5) + m(7) * (qJD(4) * t34 - t1 * t79 - t112 * t9 - t114 * t12 + t2 * t81) + m(6) * (-t112 * t15 - t114 * t16 + t3 * t81 - t4 * t79 - 0.2e1 * t109)) * t80; 0; 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (-0.1e1 + t122) * t80 * t115; m(7) * (t1 * t48 + t12 * t25 + t13 * t63 - t2 * t51 + t26 * t9) + t63 * t10 + t48 * t17 + t13 * t49 - t51 * t19 + t34 * t36 + t25 * t42 + t26 * t44 - pkin(4) * t11 + (-qJD(2) * mrSges(5,2) + t70 / 0.2e1 + t71 / 0.2e1 + (-t77 * t95 - Ifges(5,5)) * qJD(4)) * t80 + (-t1 * mrSges(7,3) - t4 * mrSges(6,3) + t7 / 0.2e1 + t8 / 0.2e1 + t102 * t117 + (-t21 / 0.2e1 - t22 / 0.2e1 + pkin(5) * t30 - t16 * mrSges(6,3) - t12 * mrSges(7,3) - t103 * t80 + t34 * t146) * qJD(5) + (-m(6) * t4 - t18 + (-m(6) * t16 - t43) * qJD(5)) * pkin(8)) * t79 + (t2 * mrSges(7,3) + t3 * mrSges(6,3) + t5 / 0.2e1 + t6 / 0.2e1 + t101 * t117 + (t23 / 0.2e1 + t24 / 0.2e1 - t9 * mrSges(7,3) - t15 * mrSges(6,3)) * qJD(5) + (t20 - qJD(5) * t45 + m(6) * (-qJD(5) * t15 + t3)) * pkin(8)) * t81 + (-t77 * t37 + (t40 / 0.2e1 + t41 / 0.2e1) * t81 + (-t38 / 0.2e1 - t39 / 0.2e1) * t79 + (-t77 * mrSges(5,2) - Ifges(5,6) + t103 * t81 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t79) * qJD(4) + t95 * qJD(2) + (t101 * t79 - t102 * t81) * qJD(5)) * t82; m(7) * (qJD(5) * t151 - t25 * t79 - t26 * t81); m(7) * (-pkin(5) * t108 + t51 * t104 - t48 * t107 + t25 * t135 - t26 * t136) + (m(6) * (t122 * t144 - t143) + (t123 + t156) * t80) * qJD(4) + (-t37 - t36 + (-m(7) * t151 - mrSges(5,2) + (mrSges(6,3) + mrSges(7,3)) * t122) * qJD(4)) * t82; -0.2e1 * pkin(4) * t37 + 0.2e1 * t63 * t36 + (-t25 * t51 + t26 * t48) * t149 + (t26 * t148 + t40 + t41 + (0.2e1 * t156 * pkin(5) - t51 * t148 - t52 - t53) * qJD(5)) * t79 + (0.2e1 * t25 * mrSges(7,3) + t38 + t39 + (t148 * t48 + t54 + t55) * qJD(5)) * t81; t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2) - t99 * t116 + (m(7) * t1 + t17) * pkin(5) + (-t153 * t79 - t158 * t81) * t111 + t96; t110 + t37; (t132 * t79 + t81 * t94) * t113 + (-t132 * t81 + t79 * t94) * t115; -mrSges(7,2) * t25 + t70 + t71 + t100 * t26 + ((-mrSges(6,1) * pkin(8) - mrSges(7,3) * pkin(5)) * t81 + (mrSges(6,2) * pkin(8) - t158) * t79) * qJD(5); 0; t152; 0; m(7) * t117; t110; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
