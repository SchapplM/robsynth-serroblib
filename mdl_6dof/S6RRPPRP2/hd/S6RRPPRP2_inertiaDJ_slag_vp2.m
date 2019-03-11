% Calculate time derivative of joint inertia matrix for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:52
% EndTime: 2019-03-09 08:29:57
% DurationCPUTime: 2.35s
% Computational Cost: add. (2465->313), mult. (5368->430), div. (0->0), fcn. (4758->6), ass. (0->142)
t152 = Ifges(6,6) + Ifges(7,6);
t100 = cos(pkin(9));
t102 = sin(qJ(2));
t104 = cos(qJ(2));
t99 = sin(pkin(9));
t73 = t100 * t102 + t104 * t99;
t177 = t152 * t73;
t153 = Ifges(6,5) + Ifges(7,5);
t176 = Ifges(6,3) + Ifges(7,3);
t101 = sin(qJ(5));
t103 = cos(qJ(5));
t141 = Ifges(7,4) * t103;
t112 = Ifges(7,1) * t101 + t141;
t143 = Ifges(6,4) * t103;
t113 = Ifges(6,1) * t101 + t143;
t135 = t100 * t104;
t72 = t102 * t99 - t135;
t147 = t153 * t73 + (t112 + t113) * t72;
t175 = t101 * t147;
t142 = Ifges(7,4) * t101;
t110 = Ifges(7,2) * t103 + t142;
t144 = Ifges(6,4) * t101;
t111 = Ifges(6,2) * t103 + t144;
t148 = t177 + (t110 + t111) * t72;
t174 = t103 * t148;
t131 = qJD(5) * t103;
t68 = t73 * qJD(2);
t107 = t101 * t68 + t72 * t131;
t132 = qJD(5) * t101;
t124 = t72 * t132;
t138 = t103 * t68;
t106 = t124 - t138;
t173 = -t101 * t153 - t152 * t103;
t169 = -2 * mrSges(7,3);
t172 = 2 * mrSges(5,1) + 2 * mrSges(4,3);
t133 = qJD(2) * t102;
t69 = qJD(2) * t135 - t133 * t99;
t98 = pkin(2) * t133;
t108 = -qJ(4) * t69 - qJD(4) * t73 + t98;
t163 = pkin(3) + pkin(8);
t17 = t163 * t68 + t108;
t151 = -qJ(3) - pkin(7);
t121 = qJD(2) * t151;
t66 = qJD(3) * t104 + t102 * t121;
t67 = -t102 * qJD(3) + t104 * t121;
t37 = -t100 * t67 + t66 * t99;
t20 = pkin(4) * t69 + t37;
t96 = -pkin(2) * t104 - pkin(1);
t109 = -t73 * qJ(4) + t96;
t33 = t163 * t72 + t109;
t82 = t151 * t102;
t85 = t151 * t104;
t48 = -t100 * t82 - t85 * t99;
t39 = pkin(4) * t73 + t48;
t3 = t101 * t20 + t103 * t17 + t39 * t131 - t132 * t33;
t19 = t103 * t20;
t8 = t101 * t39 + t103 * t33;
t4 = -qJD(5) * t8 - t101 * t17 + t19;
t171 = t101 * t3 + t103 * t4;
t170 = 2 * m(7);
t27 = pkin(3) * t68 + t108;
t168 = -0.2e1 * t27;
t43 = t72 * pkin(3) + t109;
t167 = -0.2e1 * t43;
t166 = 0.2e1 * t96;
t165 = m(7) * pkin(5);
t162 = pkin(2) * t99;
t159 = pkin(2) * t100;
t155 = mrSges(5,2) - mrSges(4,1);
t154 = mrSges(6,2) + mrSges(7,2);
t23 = -mrSges(7,2) * t69 - mrSges(7,3) * t106;
t24 = -mrSges(6,2) * t69 - mrSges(6,3) * t106;
t150 = t23 + t24;
t25 = mrSges(7,1) * t69 - mrSges(7,3) * t107;
t26 = mrSges(6,1) * t69 - mrSges(6,3) * t107;
t149 = t25 + t26;
t139 = t101 * t72;
t44 = mrSges(7,1) * t73 - mrSges(7,3) * t139;
t45 = mrSges(6,1) * t73 - mrSges(6,3) * t139;
t146 = -t44 - t45;
t137 = t103 * t72;
t46 = -mrSges(7,2) * t73 + mrSges(7,3) * t137;
t47 = -mrSges(6,2) * t73 + mrSges(6,3) * t137;
t145 = t46 + t47;
t95 = -pkin(3) - t159;
t92 = -pkin(8) + t95;
t136 = qJ(6) - t92;
t130 = qJD(6) * t103;
t129 = 0.2e1 * t68;
t86 = -Ifges(7,2) * t101 + t141;
t87 = -Ifges(6,2) * t101 + t143;
t127 = t86 / 0.2e1 + t87 / 0.2e1;
t88 = Ifges(7,1) * t103 - t142;
t89 = Ifges(6,1) * t103 - t144;
t126 = t88 / 0.2e1 + t89 / 0.2e1;
t125 = mrSges(7,1) + t165;
t93 = qJ(4) + t162;
t122 = -qJ(6) * t72 - t33;
t71 = t136 * t103;
t120 = 0.2e1 * t98;
t97 = mrSges(7,1) * t131;
t75 = -mrSges(7,2) * t132 + t97;
t36 = t103 * t39;
t5 = pkin(5) * t73 + t101 * t122 + t36;
t6 = qJ(6) * t137 + t8;
t117 = t5 * t101 - t6 * t103;
t7 = -t101 * t33 + t36;
t116 = t101 * t7 - t103 * t8;
t38 = t100 * t66 + t67 * t99;
t49 = -t100 * t85 + t82 * t99;
t115 = t37 * t48 + t38 * t49;
t114 = mrSges(6,1) * t103 - mrSges(6,2) * t101;
t14 = t106 * mrSges(7,1) + t107 * mrSges(7,2);
t105 = t153 * t107 + t152 * t138 + t176 * t69;
t21 = -pkin(4) * t68 + t38;
t90 = pkin(5) * t131 + qJD(4);
t84 = mrSges(6,1) * t101 + mrSges(6,2) * t103;
t83 = mrSges(7,1) * t101 + mrSges(7,2) * t103;
t81 = pkin(5) * t101 + t93;
t80 = t113 * qJD(5);
t79 = t112 * qJD(5);
t78 = t111 * qJD(5);
t77 = t110 * qJD(5);
t76 = t114 * qJD(5);
t70 = t136 * t101;
t63 = t69 * mrSges(5,3);
t62 = t69 * mrSges(4,2);
t51 = -qJD(5) * t71 - qJD(6) * t101;
t50 = t132 * t136 - t130;
t42 = t114 * t72;
t41 = (-mrSges(7,1) * t103 + mrSges(7,2) * t101) * t72;
t40 = -pkin(4) * t72 + t49;
t22 = (-pkin(5) * t103 - pkin(4)) * t72 + t49;
t15 = mrSges(6,1) * t106 + mrSges(6,2) * t107;
t13 = Ifges(6,1) * t107 - Ifges(6,4) * t106 + t69 * Ifges(6,5);
t12 = Ifges(7,1) * t107 - Ifges(7,4) * t106 + t69 * Ifges(7,5);
t11 = Ifges(6,4) * t107 - Ifges(6,2) * t106 + t69 * Ifges(6,6);
t10 = Ifges(7,4) * t107 - Ifges(7,2) * t106 + t69 * Ifges(7,6);
t9 = pkin(5) * t106 + t21;
t2 = -qJ(6) * t106 + t130 * t72 + t3;
t1 = pkin(5) * t69 + t19 + t122 * t131 + (-qJ(6) * t68 - qJD(5) * t39 - qJD(6) * t72 - t17) * t101;
t16 = [0.2e1 * t1 * t44 + 0.2e1 * t22 * t14 + 0.2e1 * t40 * t15 + 0.2e1 * t2 * t46 - 0.2e1 * t21 * t42 + 0.2e1 * t6 * t23 + 0.2e1 * t8 * t24 + 0.2e1 * t5 * t25 + 0.2e1 * t7 * t26 + 0.2e1 * t3 * t47 + 0.2e1 * t4 * t45 + 0.2e1 * t9 * t41 + t63 * t167 + t62 * t166 + t69 * t48 * t172 + (mrSges(4,1) * t166 + mrSges(5,2) * t167 - t172 * t49 + t174 + t175) * t68 + 0.2e1 * m(4) * (t96 * t98 + t115) + 0.2e1 * m(5) * (t27 * t43 + t115) + 0.2e1 * m(6) * (t21 * t40 + t3 * t8 + t4 * t7) + (t1 * t5 + t2 * t6 + t22 * t9) * t170 + (mrSges(4,2) * t120 + mrSges(5,3) * t168 + (-Ifges(4,4) - Ifges(5,6)) * t129 + t37 * t172 + ((2 * Ifges(4,1)) + (2 * Ifges(5,2)) + t176) * t69 + t105) * t73 + (mrSges(4,1) * t120 + mrSges(5,2) * t168 + (Ifges(5,3) + Ifges(4,2)) * t129 - t38 * t172 + (t10 + t11) * t103 + (t12 + t13) * t101 + (-0.2e1 * Ifges(4,4) - 0.2e1 * Ifges(5,6) - t173) * t69 + (t147 * t103 + (-t148 - t177) * t101) * qJD(5)) * t72 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t102 + mrSges(3,2) * t104) + (-Ifges(3,2) + Ifges(3,1)) * t102 * t104 + (-t102 ^ 2 + t104 ^ 2) * Ifges(3,4)) * qJD(2); t81 * t14 + t93 * t15 + t21 * t84 + t22 * t75 - t70 * t23 - t71 * t25 + t40 * t76 + t90 * t41 + t50 * t44 + t51 * t46 + t9 * t83 + (-mrSges(4,2) + mrSges(5,3)) * t38 + t155 * t37 + (-t72 * mrSges(5,1) - t42) * qJD(4) + (t95 * mrSges(5,1) - mrSges(4,3) * t159 - Ifges(5,4) + Ifges(4,5)) * t69 + (Ifges(3,5) * t104 - Ifges(3,6) * t102 + (-mrSges(3,1) * t104 + mrSges(3,2) * t102) * pkin(7)) * qJD(2) + (t92 * t26 - t1 * mrSges(7,3) - t4 * mrSges(6,3) + t12 / 0.2e1 + t13 / 0.2e1 + (-t77 / 0.2e1 - t78 / 0.2e1) * t72 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t69) * t103 + (t92 * t24 - t2 * mrSges(7,3) - t3 * mrSges(6,3) - t10 / 0.2e1 - t11 / 0.2e1 + (-t79 / 0.2e1 - t80 / 0.2e1) * t72 + (-Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t69) * t101 + m(5) * (qJD(4) * t49 + t37 * t95 + t38 * t93) + m(7) * (-t1 * t71 - t2 * t70 + t22 * t90 + t5 * t50 + t51 * t6 + t81 * t9) + m(6) * (qJD(4) * t40 + t171 * t92 + t21 * t93) + m(4) * (-t100 * t37 + t38 * t99) * pkin(2) + (-t93 * mrSges(5,1) - mrSges(4,3) * t162 + t101 * t126 + t103 * t127 + Ifges(5,5) - Ifges(4,6)) * t68 + (t117 * mrSges(7,3) + t116 * mrSges(6,3) + (-m(6) * t116 - t101 * t45 + t103 * t47) * t92 + (-t101 * t127 + t103 * t126) * t72 + t173 * t73 / 0.2e1 - t175 / 0.2e1 - t174 / 0.2e1) * qJD(5); 0.2e1 * t90 * t83 + 0.2e1 * t81 * t75 + (-t50 * t71 - t51 * t70 + t81 * t90) * t170 + 0.2e1 * t93 * t76 + (t50 * t169 - t79 - t80) * t103 + (t51 * t169 + t77 + t78) * t101 + 0.2e1 * (mrSges(5,3) + t84 + (m(5) + m(6)) * t93) * qJD(4) + ((-t169 * t70 - t86 - t87) * t103 + (t71 * t169 - t88 - t89) * t101) * qJD(5); m(4) * t98 + t62 - t63 - t155 * t68 + t150 * t103 - t149 * t101 + (-t145 * t101 + t146 * t103) * qJD(5) + m(7) * (-t1 * t101 + t103 * t2 + (-t101 * t6 - t103 * t5) * qJD(5)) + m(6) * (-t101 * t4 + t103 * t3 + (-t101 * t8 - t103 * t7) * qJD(5)) + m(5) * t27; m(7) * (-t101 * t50 + t103 * t51 + (t101 * t70 + t103 * t71) * qJD(5)); 0; t69 * mrSges(5,1) + t149 * t103 + t150 * t101 + (t101 * t146 + t103 * t145) * qJD(5) + m(7) * (-qJD(5) * t117 + t1 * t103 + t101 * t2) + m(6) * (-qJD(5) * t116 + t171) + m(5) * t37; m(7) * (t101 * t51 + t103 * t50 + (t101 * t71 - t103 * t70) * qJD(5)); 0; 0; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 - t152 * t124 + (m(7) * t1 + t25) * pkin(5) + t105; -mrSges(7,2) * t51 + t125 * t50 + ((-mrSges(6,2) * t92 - t152) * t103 + (-mrSges(6,1) * t92 + (mrSges(7,3) * pkin(5)) - t153) * t101) * qJD(5); -t97 + ((-mrSges(6,1) - t165) * t103 + t154 * t101) * qJD(5); (-t154 * t103 + (-mrSges(6,1) - t125) * t101) * qJD(5); 0; m(7) * t9 + t14; m(7) * t90 + t75; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
