% Calculate time derivative of joint inertia matrix for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:14
% EndTime: 2019-03-09 02:33:19
% DurationCPUTime: 2.28s
% Computational Cost: add. (4725->250), mult. (9230->381), div. (0->0), fcn. (9441->8), ass. (0->120)
t105 = sin(qJ(5));
t107 = cos(qJ(5));
t101 = sin(pkin(10));
t102 = cos(pkin(10));
t178 = sin(qJ(4));
t179 = cos(qJ(4));
t114 = t178 * t101 - t179 * t102;
t79 = -t179 * t101 - t178 * t102;
t60 = t105 * t79 - t107 * t114;
t136 = qJD(4) * t178;
t137 = qJD(4) * t179;
t71 = -t101 * t136 + t102 * t137;
t72 = t79 * qJD(4);
t111 = t60 * qJD(5) + t105 * t72 + t107 * t71;
t104 = sin(qJ(6));
t106 = cos(qJ(6));
t147 = t104 ^ 2 + t106 ^ 2;
t195 = t147 * t111;
t186 = t147 * t107;
t120 = -t105 * t114 - t107 * t79;
t41 = t120 * qJD(5) + t105 * t71 - t107 * t72;
t173 = t41 * t60;
t194 = t111 * t120 - t173;
t103 = -pkin(1) - qJ(3);
t169 = -pkin(7) + t103;
t84 = t169 * t101;
t85 = t169 * t102;
t62 = t178 * t85 + t179 * t84;
t52 = t114 * qJD(3) - t62 * qJD(4);
t110 = -t72 * pkin(8) + t52;
t61 = -t178 * t84 + t179 * t85;
t116 = pkin(8) * t114 + t61;
t56 = pkin(8) * t79 + t62;
t23 = t105 * t56 - t107 * t116;
t51 = t79 * qJD(3) - t84 * t136 + t85 * t137;
t46 = -t71 * pkin(8) + t51;
t10 = -qJD(5) * t23 + t105 * t110 + t107 * t46;
t24 = t105 * t116 + t107 * t56;
t90 = t101 * pkin(3) + qJ(2);
t63 = -pkin(4) * t79 + t90;
t25 = pkin(5) * t120 - pkin(9) * t60 + t63;
t13 = -t104 * t24 + t106 * t25;
t64 = pkin(4) * t71 + qJD(2);
t17 = pkin(5) * t111 + pkin(9) * t41 + t64;
t2 = qJD(6) * t13 + t10 * t106 + t104 * t17;
t14 = t104 * t25 + t106 * t24;
t3 = -qJD(6) * t14 - t10 * t104 + t106 * t17;
t193 = -t104 * t3 + t106 * t2;
t192 = -t105 * t111 + t107 * t41;
t166 = mrSges(7,1) * t106;
t86 = mrSges(7,2) * t104 - t166;
t191 = t86 - mrSges(6,1);
t190 = (-t104 * t14 - t106 * t13) * qJD(6) + t193;
t144 = qJD(6) * t106;
t119 = -t104 * t41 + t60 * t144;
t123 = -t104 * t13 + t106 * t14;
t131 = (t101 ^ 2 + t102 ^ 2) * qJD(3);
t185 = 2 * m(7);
t11 = qJD(5) * t24 + t105 * t46 - t107 * t110;
t184 = 0.2e1 * t11;
t183 = 0.2e1 * t23;
t182 = 0.2e1 * t64;
t181 = m(6) / 0.2e1;
t175 = t11 * t23;
t172 = t111 * mrSges(6,3);
t171 = t71 * t79;
t170 = t72 * t114;
t151 = t106 * t41;
t168 = -Ifges(7,5) * t151 + Ifges(7,3) * t111;
t165 = Ifges(7,4) * t104;
t164 = Ifges(7,4) * t106;
t163 = Ifges(7,6) * t104;
t162 = pkin(4) * qJD(5);
t158 = t104 * t60;
t157 = t105 * mrSges(6,1);
t155 = t60 * t105;
t153 = t105 * t86;
t150 = t106 * t60;
t149 = t107 * mrSges(6,2);
t145 = qJD(6) * t104;
t143 = 2 * mrSges(5,3);
t142 = t60 * t145;
t118 = t142 + t151;
t12 = t119 * mrSges(7,1) - t118 * mrSges(7,2);
t140 = m(7) * t11 + t12;
t139 = t71 * mrSges(5,1) + t72 * mrSges(5,2);
t138 = mrSges(6,1) * t111 - mrSges(6,2) * t41;
t133 = -(2 * Ifges(6,4)) - t163;
t130 = mrSges(7,3) * t186;
t129 = -t11 * t60 + t23 * t41;
t128 = t170 + t171;
t127 = mrSges(7,1) * t104 + mrSges(7,2) * t106;
t126 = Ifges(7,1) * t106 - t165;
t125 = -Ifges(7,2) * t104 + t164;
t124 = Ifges(7,5) * t104 + Ifges(7,6) * t106;
t15 = mrSges(7,1) * t111 + t118 * mrSges(7,3);
t16 = -mrSges(7,2) * t111 - t119 * mrSges(7,3);
t122 = -t104 * t15 + t106 * t16;
t43 = -mrSges(7,2) * t120 - mrSges(7,3) * t158;
t44 = mrSges(7,1) * t120 - mrSges(7,3) * t150;
t121 = -t104 * t44 + t106 * t43;
t82 = t125 * qJD(6);
t83 = t126 * qJD(6);
t87 = Ifges(7,2) * t106 + t165;
t88 = Ifges(7,1) * t104 + t164;
t117 = t104 * t83 + t106 * t82 + t88 * t144 - t87 * t145;
t81 = t127 * qJD(6);
t115 = -mrSges(6,2) * t111 + mrSges(7,3) * t195 + t191 * t41 - t60 * t81;
t113 = t114 * t52 + t51 * t79 - t61 * t72 - t62 * t71;
t109 = -t44 * t144 - t43 * t145 + m(7) * (-t13 * t144 - t14 * t145 + t193) + t122;
t20 = Ifges(7,6) * t120 + t125 * t60;
t21 = Ifges(7,5) * t120 + t126 * t60;
t6 = -t118 * Ifges(7,4) - t119 * Ifges(7,2) + Ifges(7,6) * t111;
t7 = -t118 * Ifges(7,1) - t119 * Ifges(7,4) + Ifges(7,5) * t111;
t95 = Ifges(7,5) * t144;
t108 = -t10 * mrSges(6,2) + t21 * t144 / 0.2e1 + t23 * t81 - t88 * t151 / 0.2e1 - Ifges(6,5) * t41 + t104 * t7 / 0.2e1 - t82 * t158 / 0.2e1 + t83 * t150 / 0.2e1 + t106 * t6 / 0.2e1 + t120 * (-Ifges(7,6) * t145 + t95) / 0.2e1 - t119 * t87 / 0.2e1 - (t60 * t88 + t20) * t145 / 0.2e1 + (t124 / 0.2e1 - Ifges(6,6)) * t111 + t191 * t11 + t190 * mrSges(7,3);
t94 = -pkin(4) * t107 - pkin(5);
t93 = pkin(4) * t105 + pkin(9);
t33 = t127 * t60;
t1 = [-0.2e1 * Ifges(5,1) * t170 - 0.2e1 * Ifges(5,2) * t171 + 0.2e1 * t90 * t139 + 0.2e1 * t63 * t138 + 0.2e1 * t3 * t44 + 0.2e1 * t2 * t43 + t33 * t184 + t12 * t183 + 0.2e1 * t13 * t15 + 0.2e1 * t14 * t16 - 0.2e1 * t24 * t172 - (mrSges(6,3) * t183 - t104 * t20 + t106 * t21) * t41 + t113 * t143 + (t13 * t3 + t14 * t2 + t175) * t185 + 0.2e1 * m(6) * (t10 * t24 + t63 * t64 + t175) + 0.2e1 * m(5) * (qJD(2) * t90 + t51 * t62 + t52 * t61) + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t103 * t131) + (mrSges(6,1) * t182 - 0.2e1 * t10 * mrSges(6,3) + ((2 * Ifges(6,2)) + Ifges(7,3)) * t111 - t133 * t41 + t168) * t120 + (mrSges(6,2) * t182 + mrSges(6,3) * t184 - 0.2e1 * Ifges(6,1) * t41 - t104 * t6 + t106 * t7 + (Ifges(7,5) * t106 + t133) * t111 + (-t104 * t21 - t106 * t20 - t120 * t124) * qJD(6)) * t60 + 0.2e1 * (m(3) * qJ(2) + mrSges(4,1) * t101 - mrSges(5,1) * t79 + mrSges(4,2) * t102 - mrSges(5,2) * t114 + mrSges(3,3)) * qJD(2) + 0.2e1 * mrSges(4,3) * t131 + 0.2e1 * (t114 * t71 + t72 * t79) * Ifges(5,4); -t60 * t12 + t41 * t33 + t121 * t111 + t128 * t143 + (-t194 + t173) * mrSges(6,3) + (-t172 + (-t104 * t43 - t106 * t44) * qJD(6) + t122) * t120 + m(7) * (t123 * t111 + t190 * t120 + t129) + m(6) * (t10 * t120 + t111 * t24 + t129) - m(5) * t113 - m(4) * t131; 0.2e1 * m(7) * (t120 * t195 - t173) + 0.2e1 * m(6) * t194 - 0.2e1 * m(5) * t128; t104 * t16 + t106 * t15 + t121 * qJD(6) + (m(5) + m(4)) * qJD(2) + m(7) * (t123 * qJD(6) + t104 * t2 + t106 * t3) + m(6) * t64 + t138 + t139; 0; 0; t108 + (m(6) * (t10 * t105 - t107 * t11) + t192 * mrSges(6,3) + ((m(6) * t24 + m(7) * t123 - t120 * mrSges(6,3) + t121) * t107 + (t60 * mrSges(6,3) + t33 + (m(6) + m(7)) * t23) * t105) * qJD(5)) * pkin(4) + t109 * t93 + t140 * t94 + Ifges(5,5) * t72 - Ifges(5,6) * t71 - t51 * mrSges(5,2) + t52 * mrSges(5,1); m(7) * (t195 * t93 + t94 * t41) - t71 * mrSges(5,2) + t72 * mrSges(5,1) + 0.2e1 * (-t192 * t181 + (m(7) * (t186 * t120 - t155) / 0.2e1 + (t107 * t120 - t155) * t181) * qJD(5)) * pkin(4) + t115; 0; 0.2e1 * t94 * t81 + (-0.2e1 * t149 - 0.2e1 * t157 + (t105 * t94 + t186 * t93) * t185 + 0.2e1 * t153 + 0.2e1 * t130) * t162 + t117; -t140 * pkin(5) + t109 * pkin(9) + t108; m(7) * (-pkin(5) * t41 + pkin(9) * t195) + t115; 0; (t94 - pkin(5)) * t81 + (-t149 - t157 + m(7) * (-pkin(5) * t105 + t186 * pkin(9)) + t153 + t130) * t162 + t117; -0.2e1 * pkin(5) * t81 + t117; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t142 - t119 * Ifges(7,6) + t168; (-t106 * t111 + t120 * t145) * mrSges(7,2) + (-t104 * t111 - t120 * t144) * mrSges(7,1); -t81; t95 - t127 * t107 * t162 + (-t93 * t166 + (mrSges(7,2) * t93 - Ifges(7,6)) * t104) * qJD(6); t95 + (t86 * pkin(9) - t163) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
