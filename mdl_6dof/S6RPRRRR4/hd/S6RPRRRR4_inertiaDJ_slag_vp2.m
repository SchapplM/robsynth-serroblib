% Calculate time derivative of joint inertia matrix for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:00
% EndTime: 2019-03-09 07:05:07
% DurationCPUTime: 3.37s
% Computational Cost: add. (12754->316), mult. (27058->483), div. (0->0), fcn. (29610->10), ass. (0->154)
t125 = sin(qJ(6));
t121 = t125 ^ 2;
t129 = cos(qJ(6));
t122 = t129 ^ 2;
t214 = t121 + t122;
t127 = sin(qJ(4));
t131 = cos(qJ(4));
t123 = sin(pkin(11));
t124 = cos(pkin(11));
t128 = sin(qJ(3));
t132 = cos(qJ(3));
t102 = t123 * t132 + t128 * t124;
t192 = pkin(7) + qJ(2);
t106 = t192 * t123;
t107 = t192 * t124;
t212 = -t132 * t106 - t107 * t128;
t72 = -pkin(8) * t102 + t212;
t101 = -t128 * t123 + t124 * t132;
t84 = -t128 * t106 + t132 * t107;
t73 = pkin(8) * t101 + t84;
t47 = -t127 * t73 + t131 * t72;
t164 = qJD(6) * t129;
t126 = sin(qJ(5));
t130 = cos(qJ(5));
t92 = t101 * qJD(3);
t93 = t102 * qJD(3);
t145 = t127 * t93 - t131 * t92;
t80 = t101 * t131 - t102 * t127;
t58 = t80 * qJD(4) - t145;
t81 = t101 * t127 + t102 * t131;
t59 = -t81 * qJD(4) - t127 * t92 - t131 * t93;
t60 = t126 * t81 - t130 * t80;
t35 = -t60 * qJD(5) + t126 * t59 + t130 * t58;
t61 = t126 * t80 + t130 * t81;
t141 = t125 * t35 + t61 * t164;
t211 = t214 * t130;
t142 = -pkin(9) * t81 + t47;
t48 = t127 * t72 + t131 * t73;
t45 = pkin(9) * t80 + t48;
t22 = t126 * t142 + t130 * t45;
t159 = -pkin(2) * t124 - pkin(1);
t85 = -pkin(3) * t101 + t159;
t64 = -pkin(4) * t80 + t85;
t40 = pkin(5) * t60 - pkin(10) * t61 + t64;
t16 = -t125 * t22 + t129 * t40;
t17 = t125 * t40 + t129 * t22;
t210 = -t125 * t16 + t129 * t17;
t67 = t101 * qJD(2) + qJD(3) * t212;
t68 = -t102 * qJD(2) - qJD(3) * t84;
t38 = t127 * (-t92 * pkin(8) + t68) + t131 * (-t93 * pkin(8) + t67) + t47 * qJD(4);
t39 = -t48 * qJD(4) + t145 * pkin(8) + (-t127 * t212 - t131 * t84) * qJD(3) - t81 * qJD(2);
t209 = t39 * mrSges(5,1) - t38 * mrSges(5,2) + Ifges(5,5) * t58 + Ifges(5,6) * t59;
t166 = qJD(5) * t130;
t163 = pkin(4) * t166;
t152 = mrSges(7,3) * t163;
t111 = t121 * t152;
t112 = t122 * t152;
t150 = mrSges(7,1) * t125 + mrSges(7,2) * t129;
t103 = t150 * qJD(6);
t116 = -pkin(4) * t130 - pkin(5);
t86 = t116 * t103;
t190 = mrSges(7,1) * t129;
t108 = mrSges(7,2) * t125 - t190;
t167 = qJD(5) * t126;
t96 = pkin(4) * t108 * t167;
t208 = t111 + t112 + t86 + t96;
t207 = 2 * m(6);
t206 = 2 * m(7);
t133 = -t58 * pkin(9) + t39;
t19 = pkin(9) * t59 + t38;
t9 = t22 * qJD(5) + t126 * t19 - t130 * t133;
t205 = 0.2e1 * t9;
t204 = -2 * mrSges(6,3);
t21 = t126 * t45 - t130 * t142;
t203 = 0.2e1 * t21;
t200 = t93 * pkin(3);
t46 = -pkin(4) * t59 + t200;
t202 = 0.2e1 * t46;
t201 = t21 * t9;
t198 = pkin(5) * t103;
t36 = t61 * qJD(5) + t126 * t58 - t130 * t59;
t13 = pkin(5) * t36 - pkin(10) * t35 + t46;
t175 = qJD(6) * t16;
t8 = -t21 * qJD(5) + t126 * t133 + t130 * t19;
t2 = t125 * t13 + t129 * t8 + t175;
t197 = t129 * t2;
t117 = pkin(3) * t131 + pkin(4);
t168 = t127 * t130;
t79 = t117 * t167 + (t127 * t166 + (t126 * t131 + t168) * qJD(4)) * pkin(3);
t196 = t21 * t79;
t174 = qJD(6) * t17;
t3 = -t125 * t8 + t129 * t13 - t174;
t195 = t3 * t125;
t169 = t126 * t127;
t78 = t117 * t166 + (-t127 * t167 + (t130 * t131 - t169) * qJD(4)) * pkin(3);
t193 = t78 * mrSges(6,2);
t177 = t129 * t35;
t191 = Ifges(7,5) * t177 + Ifges(7,3) * t36;
t189 = Ifges(7,4) * t125;
t188 = Ifges(7,4) * t129;
t187 = Ifges(7,6) * t125;
t186 = pkin(4) * qJD(5);
t185 = t121 * t78;
t184 = t122 * t78;
t181 = t125 * t61;
t176 = t129 * t61;
t91 = pkin(3) * t168 + t126 * t117;
t89 = pkin(10) + t91;
t173 = qJD(6) * t89;
t165 = qJD(6) * t125;
t161 = t61 * t165;
t140 = t161 - t177;
t12 = t141 * mrSges(7,1) - t140 * mrSges(7,2);
t162 = m(7) * t9 + t12;
t158 = -t59 * mrSges(5,1) + t58 * mrSges(5,2);
t157 = t36 * mrSges(6,1) + t35 * mrSges(6,2);
t156 = -t165 / 0.2e1;
t155 = -(2 * Ifges(6,4)) - t187;
t154 = t214 * t78;
t151 = -t126 * mrSges(6,1) - t130 * mrSges(6,2);
t149 = Ifges(7,1) * t129 - t189;
t148 = -Ifges(7,2) * t125 + t188;
t147 = Ifges(7,5) * t125 + Ifges(7,6) * t129;
t42 = -mrSges(7,2) * t60 - mrSges(7,3) * t181;
t43 = mrSges(7,1) * t60 - mrSges(7,3) * t176;
t146 = -t125 * t43 + t129 * t42;
t90 = -pkin(3) * t169 + t117 * t130;
t104 = t148 * qJD(6);
t105 = t149 * qJD(6);
t109 = Ifges(7,2) * t129 + t189;
t110 = Ifges(7,1) * t125 + t188;
t139 = t129 * t104 + t125 * t105 - t109 * t165 + t110 * t164;
t138 = (-mrSges(5,1) * t127 - mrSges(5,2) * t131) * qJD(4) * pkin(3);
t66 = t79 * t108;
t74 = mrSges(7,3) * t185;
t75 = mrSges(7,3) * t184;
t76 = t79 * mrSges(6,1);
t88 = -pkin(5) - t90;
t82 = t88 * t103;
t137 = t139 + t66 + t74 + t75 - t76 + t82 - t193;
t10 = -t140 * Ifges(7,4) - t141 * Ifges(7,2) + t36 * Ifges(7,6);
t11 = -t140 * Ifges(7,1) - t141 * Ifges(7,4) + t36 * Ifges(7,5);
t118 = Ifges(7,5) * t164;
t33 = t60 * Ifges(7,6) + t148 * t61;
t34 = t60 * Ifges(7,5) + t149 * t61;
t136 = -t8 * mrSges(6,2) + mrSges(7,3) * t197 + t21 * t103 + t33 * t156 + t34 * t164 / 0.2e1 + Ifges(6,5) * t35 - t104 * t181 / 0.2e1 + t105 * t176 / 0.2e1 + t60 * (-Ifges(7,6) * t165 + t118) / 0.2e1 + t125 * t11 / 0.2e1 + t129 * t10 / 0.2e1 + (-mrSges(6,1) + t108) * t9 + (t147 / 0.2e1 - Ifges(6,6)) * t36 - t141 * t109 / 0.2e1 + (t177 / 0.2e1 + t61 * t156) * t110;
t14 = mrSges(7,1) * t36 + t140 * mrSges(7,3);
t15 = -mrSges(7,2) * t36 - t141 * mrSges(7,3);
t135 = -t43 * t164 - t42 * t165 + m(7) * (-t16 * t164 - t17 * t165 - t195 + t197) + t129 * t15 - t125 * t14;
t134 = t136 + (-t195 + (-t125 * t17 - t129 * t16) * qJD(6)) * mrSges(7,3);
t115 = pkin(4) * t126 + pkin(10);
t87 = t92 * mrSges(4,2);
t41 = t150 * t61;
t1 = [0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t123 ^ 2 + t124 ^ 2) * qJD(2) + 0.2e1 * t64 * t157 + 0.2e1 * t85 * t158 + 0.2e1 * t159 * (t93 * mrSges(4,1) + t87) + 0.2e1 * (t101 * t92 - t102 * t93) * Ifges(4,4) + 0.2e1 * m(5) * (t85 * t200 + t38 * t48 + t39 * t47) + t12 * t203 + t41 * t205 + (t16 * t3 + t17 * t2 + t201) * t206 + (t22 * t8 + t46 * t64 + t201) * t207 + 0.2e1 * (t80 * t58 + t81 * t59) * Ifges(5,4) + 0.2e1 * (t38 * t80 - t39 * t81 - t47 * t58 + t48 * t59) * mrSges(5,3) + (mrSges(6,1) * t202 + t8 * t204 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t36 + t155 * t35 + t191) * t60 + (mrSges(6,2) * t202 + mrSges(6,3) * t205 + 0.2e1 * Ifges(6,1) * t35 - t125 * t10 + t129 * t11 + (Ifges(7,5) * t129 + t155) * t36 + (-t125 * t34 - t129 * t33 - t60 * t147) * qJD(6)) * t61 + (mrSges(6,3) * t203 - t125 * t33 + t129 * t34) * t35 + 0.2e1 * (-mrSges(5,1) * t80 + mrSges(5,2) * t81) * t200 + 0.2e1 * t80 * Ifges(5,2) * t59 - 0.2e1 * t101 * Ifges(4,2) * t93 + 0.2e1 * t81 * t58 * Ifges(5,1) + 0.2e1 * t102 * t92 * Ifges(4,1) + 0.2e1 * (t101 * t67 - t102 * t68 - t212 * t92 - t84 * t93) * mrSges(4,3) + 0.2e1 * m(4) * (t212 * t68 + t67 * t84) + t22 * t36 * t204 + 0.2e1 * t16 * t14 + 0.2e1 * t17 * t15 + 0.2e1 * t2 * t42 + 0.2e1 * t3 * t43; t125 * t15 + t129 * t14 + t87 - (-m(5) * pkin(3) - mrSges(4,1)) * t93 + t146 * qJD(6) + m(7) * (qJD(6) * t210 + t125 * t2 + t129 * t3) + m(6) * t46 + t157 + t158; 0; m(7) * (t88 * t9 + t196) + t136 + (-t42 * t173 + m(7) * (-t16 * t78 - t17 * t173 - t3 * t89) - t89 * t14 - t78 * t43 + (-t3 - t174) * mrSges(7,3)) * t125 + (-mrSges(7,3) * t175 - t43 * t173 + t78 * t42 + m(7) * (-t16 * t173 + t17 * t78 + t2 * t89) + t89 * t15) * t129 + (m(5) * (t127 * t38 + t131 * t39 + (-t127 * t47 + t131 * t48) * qJD(4)) + (t127 * t59 - t131 * t58 + (t127 * t81 + t131 * t80) * qJD(4)) * mrSges(5,3)) * pkin(3) + (-t35 * t90 - t36 * t91 - t60 * t78 + t61 * t79) * mrSges(6,3) + m(6) * (t22 * t78 + t8 * t91 - t9 * t90 + t196) - t67 * mrSges(4,2) + t68 * mrSges(4,1) + t79 * t41 + t88 * t12 + Ifges(4,5) * t92 - Ifges(4,6) * t93 + t209; 0; -0.2e1 * t193 + 0.2e1 * t66 + 0.2e1 * t74 + 0.2e1 * t75 - 0.2e1 * t76 + 0.2e1 * t82 + 0.2e1 * t138 + (t89 * t154 + t79 * t88) * t206 + (t78 * t91 - t79 * t90) * t207 + t139; t134 + (m(6) * (t126 * t8 - t130 * t9) + (-t126 * t36 - t130 * t35) * mrSges(6,3) + ((m(6) * t22 + m(7) * t210 - t60 * mrSges(6,3) + t146) * t130 + (t61 * mrSges(6,3) + t41 + (m(6) + m(7)) * t21) * t126) * qJD(5)) * pkin(4) + t162 * t116 + t135 * t115 + t209; 0; t137 + (m(6) * (t126 * t78 - t130 * t79) + (m(7) * (t126 * t88 + t211 * t89) + m(6) * (-t126 * t90 + t130 * t91) + t151) * qJD(5)) * pkin(4) + m(7) * (t116 * t79 + (t184 + t185) * t115) + t138 + t208; 0.2e1 * t111 + 0.2e1 * t112 + 0.2e1 * t86 + 0.2e1 * t96 + 0.2e1 * (m(7) * (t115 * t211 + t116 * t126) + t151) * t186 + t139; -t162 * pkin(5) + t135 * pkin(10) + t134; 0; -t198 + m(7) * (-pkin(5) * t79 + pkin(10) * t154) + t137; -t198 + (m(7) * (-pkin(5) * t126 + pkin(10) * t211) + t151) * t186 + t139 + t208; t139 - 0.2e1 * t198; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t161 - t141 * Ifges(7,6) + t191; -t103; t118 - t150 * t78 + (-t89 * t190 + (mrSges(7,2) * t89 - Ifges(7,6)) * t125) * qJD(6); t118 - t150 * t163 + (t108 * t115 - t187) * qJD(6); t118 + (t108 * pkin(10) - t187) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
