% Calculate time derivative of joint inertia matrix for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:48
% EndTime: 2019-03-09 02:02:52
% DurationCPUTime: 1.74s
% Computational Cost: add. (1223->274), mult. (2661->390), div. (0->0), fcn. (1759->6), ass. (0->125)
t150 = Ifges(7,4) + Ifges(6,5);
t71 = sin(qJ(5));
t73 = cos(qJ(5));
t158 = -Ifges(7,6) * t71 - t150 * t73;
t141 = 2 * qJD(3);
t157 = Ifges(7,2) + Ifges(6,3);
t118 = t71 ^ 2 + t73 ^ 2;
t156 = (m(6) / 0.2e1 + m(7) / 0.2e1) * (t118 - 0.1e1);
t59 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t72 = sin(qJ(4));
t130 = t59 * t72;
t60 = sin(pkin(9)) * pkin(1) + qJ(3);
t74 = cos(qJ(4));
t31 = t72 * pkin(4) - pkin(8) * t74 + t60;
t149 = t73 * t130 + t71 * t31;
t155 = qJD(5) * t149;
t154 = -mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t118;
t11 = qJ(6) * t72 + t149;
t131 = t31 * t73;
t97 = t59 * t71 - pkin(5);
t12 = t72 * t97 - t131;
t127 = t73 * t74;
t44 = t72 * mrSges(6,1) - mrSges(6,3) * t127;
t45 = -t72 * mrSges(7,1) + mrSges(7,2) * t127;
t121 = t44 - t45;
t129 = t71 * t74;
t43 = -t72 * mrSges(6,2) - mrSges(6,3) * t129;
t46 = -mrSges(7,2) * t129 + t72 * mrSges(7,3);
t122 = t43 + t46;
t112 = qJD(5) * t74;
t101 = t73 * t112;
t116 = qJD(4) * t72;
t104 = t71 * t116;
t78 = -t101 + t104;
t79 = t112 * t71 + t116 * t73;
t13 = -mrSges(7,1) * t78 + mrSges(7,3) * t79;
t14 = -mrSges(6,1) * t78 - mrSges(6,2) * t79;
t15 = -t130 * t71 + t131;
t111 = qJD(6) * t73;
t87 = pkin(5) * t73 + qJ(6) * t71;
t145 = qJD(5) * t87 - t111;
t86 = pkin(5) * t71 - qJ(6) * t73;
t81 = -t59 + t86;
t5 = -t81 * t116 + t145 * t74;
t153 = m(7) * t5 + (m(6) * (-t149 * t73 + t15 * t71) - m(7) * (t11 * t73 + t12 * t71) + t121 * t71 - t122 * t73) * qJD(4) + t13 + t14;
t115 = qJD(4) * t74;
t148 = t72 * t115 * t156;
t147 = qJD(4) * (t72 ^ 2 - t74 ^ 2);
t146 = m(7) * qJ(6) + mrSges(7,3);
t139 = pkin(8) * t72;
t140 = pkin(4) * t74;
t35 = qJD(3) + (t139 + t140) * qJD(4);
t144 = t73 * t35 - t155;
t143 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t138 = Ifges(6,4) * t71;
t137 = Ifges(6,4) * t73;
t136 = Ifges(7,5) * t71;
t135 = Ifges(7,5) * t73;
t134 = Ifges(6,6) * t72;
t132 = Ifges(7,6) * t72;
t19 = mrSges(6,1) * t115 + mrSges(6,3) * t79;
t20 = -mrSges(7,1) * t115 - mrSges(7,2) * t79;
t125 = -t19 + t20;
t21 = -mrSges(6,2) * t115 + mrSges(6,3) * t78;
t22 = mrSges(7,2) * t78 + mrSges(7,3) * t115;
t124 = t21 + t22;
t88 = Ifges(7,3) * t71 + t135;
t24 = t74 * t88 + t132;
t89 = -Ifges(6,2) * t71 + t137;
t25 = t74 * t89 + t134;
t123 = t24 - t25;
t49 = -t73 * mrSges(6,1) + t71 * mrSges(6,2);
t120 = t49 - mrSges(5,1);
t119 = t118 * pkin(8) * t115;
t114 = qJD(5) * t71;
t113 = qJD(5) * t73;
t110 = qJ(6) * qJD(4);
t48 = -t73 * mrSges(7,1) - t71 * mrSges(7,3);
t108 = t48 + t120;
t105 = t59 * t115;
t107 = t73 * t105 + t31 * t113 + t71 * t35;
t103 = t71 * t115;
t102 = t59 * t114;
t50 = -Ifges(7,3) * t73 + t136;
t51 = Ifges(6,2) * t73 + t138;
t99 = -t51 / 0.2e1 + t50 / 0.2e1;
t52 = Ifges(7,1) * t71 - t135;
t53 = Ifges(6,1) * t71 + t137;
t98 = -t52 / 0.2e1 - t53 / 0.2e1;
t96 = Ifges(6,6) * t104 + Ifges(7,6) * t101 + t157 * t115;
t29 = qJD(5) * t86 - qJD(6) * t71;
t92 = t71 * mrSges(7,1) - t73 * mrSges(7,3);
t37 = t92 * qJD(5);
t93 = t71 * mrSges(6,1) + t73 * mrSges(6,2);
t38 = t93 * qJD(5);
t94 = m(7) * t29 + t37 + t38;
t91 = Ifges(6,1) * t73 - t138;
t90 = Ifges(7,1) * t73 + t136;
t26 = Ifges(7,4) * t72 + t74 * t90;
t27 = Ifges(6,5) * t72 + t74 * t91;
t83 = -t150 * t72 - t26 - t27;
t80 = t118 * t139;
t3 = -t102 * t72 + t107;
t4 = -t103 * t59 + t144;
t77 = -t113 * t15 - t114 * t149 + t3 * t73 - t4 * t71;
t76 = m(7) * t111 + (-m(7) * t87 + t48 + t49) * qJD(5);
t1 = t74 * t110 + (qJD(6) - t102) * t72 + t107;
t18 = t81 * t74;
t2 = t115 * t97 - t144;
t32 = t92 * t74;
t33 = t93 * t74;
t75 = m(7) * (qJD(4) * t18 + t1 * t73 - t11 * t114 + t113 * t12 + t2 * t71) + t125 * t71 + t124 * t73 + (t32 + t33) * qJD(4) + (-t121 * t73 - t122 * t71) * qJD(5);
t65 = Ifges(7,4) * t113;
t64 = Ifges(6,5) * t113;
t62 = Ifges(7,6) * t114;
t47 = -pkin(4) - t87;
t42 = t91 * qJD(5);
t41 = t90 * qJD(5);
t40 = t89 * qJD(5);
t39 = t88 * qJD(5);
t9 = -t53 * t112 + (Ifges(6,5) * t74 - t72 * t91) * qJD(4);
t8 = -t52 * t112 + (Ifges(7,4) * t74 - t72 * t90) * qJD(4);
t7 = -t51 * t112 + (Ifges(6,6) * t74 - t72 * t89) * qJD(4);
t6 = -t50 * t112 + (Ifges(7,6) * t74 - t72 * t88) * qJD(4);
t10 = [0.2e1 * t1 * t46 + 0.2e1 * t11 * t22 + 0.2e1 * t12 * t20 + 0.2e1 * t18 * t13 + 0.2e1 * t15 * t19 + 0.2e1 * t149 * t21 + 0.2e1 * t2 * t45 + 0.2e1 * t3 * t43 + 0.2e1 * t5 * t32 + 0.2e1 * t4 * t44 + 0.2e1 * m(7) * (t1 * t11 + t12 * t2 + t18 * t5) + 0.2e1 * m(6) * (t149 * t3 + t15 * t4) + (mrSges(4,3) + (m(4) + m(5)) * t60) * t141 + (mrSges(5,1) * t141 + (-0.2e1 * t60 * mrSges(5,2) + 0.2e1 * Ifges(5,4) * t72 + 0.2e1 * t59 * t33 + (-t123 - t132) * t71 + t83 * t73) * qJD(4) + t96) * t72 + (mrSges(5,2) * t141 - 0.2e1 * t59 * t14 + (t8 + t9) * t73 + (t6 - t7) * t71 + ((t123 - t134) * t73 + t83 * t71) * qJD(5) + (0.2e1 * t60 * mrSges(5,1) + (-Ifges(6,6) * t71 - 0.2e1 * Ifges(5,4) - t158) * t74 + (-0.2e1 * m(6) * t59 ^ 2 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + t157) * t72) * qJD(4)) * t74; m(6) * t59 * t147 + t153 * t72 + (m(6) * t77 + t75) * t74; -0.4e1 * t148; -t153 * t74 + (m(6) * (t77 - 0.2e1 * t105) + t75) * t72; -0.2e1 * t147 * t156; 0.4e1 * t148; t29 * t32 + t18 * t37 + t47 * t13 + t5 * t48 - pkin(4) * t14 + m(7) * (t18 * t29 + t47 * t5) + (t65 / 0.2e1 + t62 / 0.2e1 + t64 / 0.2e1 + (-Ifges(5,5) + (-m(6) * pkin(4) + t120) * t59) * qJD(4)) * t72 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(7,2) - t4 * mrSges(6,3) - t99 * t116 + (-t134 / 0.2e1 + t24 / 0.2e1 - t25 / 0.2e1 - t11 * mrSges(7,2) - t149 * mrSges(6,3)) * qJD(5) + (-t122 * qJD(5) + m(7) * (-qJD(5) * t11 + t2) + m(6) * (-t4 - t155) + t125) * pkin(8)) * t71 + (-t6 / 0.2e1 + t7 / 0.2e1 + t3 * mrSges(6,3) + t1 * mrSges(7,2) + t98 * t116 + (t26 / 0.2e1 + t27 / 0.2e1 + t12 * mrSges(7,2) - t15 * mrSges(6,3)) * qJD(5) + (-t121 * qJD(5) + m(7) * (qJD(5) * t12 + t1) + m(6) * (-qJD(5) * t15 + t3) + t124) * pkin(8)) * t73 + (-t59 * t38 + (t41 / 0.2e1 + t42 / 0.2e1) * t73 + (t39 / 0.2e1 - t40 / 0.2e1) * t71 + (-t59 * mrSges(5,2) - Ifges(5,6) + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t73 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t71) * qJD(4) + (t71 * t98 + t73 * t99) * qJD(5)) * t74; t94 * t72 + (t108 * t74 - t154 * t72 + m(6) * (-t80 - t140) + m(7) * (t47 * t74 - t80)) * qJD(4); t108 * t116 + m(7) * (t116 * t47 + t119) + m(6) * (-pkin(4) * t116 + t119) + (t154 * qJD(4) - t94) * t74; -0.2e1 * pkin(4) * t38 + 0.2e1 * t37 * t47 + (-t39 + t40) * t73 + (t41 + t42) * t71 + 0.2e1 * (m(7) * t47 + t48) * t29 + ((t52 + t53) * t73 + (t50 - t51) * t71) * qJD(5); m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t11) + qJD(6) * t46 + qJ(6) * t22 + t1 * mrSges(7,3) - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - pkin(5) * t20 + (-Ifges(6,6) * t73 - t150 * t71) * t112 + t158 * t116 + t96; m(7) * (qJD(6) * t74 - t110 * t72) * t73 + (t143 * t71 + (-mrSges(7,3) + mrSges(6,2)) * t73) * t116 + ((mrSges(6,2) - t146) * t71 - t143 * t73) * t112; (-m(7) * t86 - t92 - t93) * t115 + t76 * t72; -t145 * mrSges(7,2) - Ifges(6,6) * t114 + t76 * pkin(8) + t62 + t64 + t65; 0.2e1 * t146 * qJD(6); m(7) * t2 + t20; -t78 * m(7); (t113 * t72 + t103) * m(7); (m(7) * pkin(8) + mrSges(7,2)) * t113; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
