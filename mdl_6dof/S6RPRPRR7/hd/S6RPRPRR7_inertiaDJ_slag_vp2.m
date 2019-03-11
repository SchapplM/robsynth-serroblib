% Calculate time derivative of joint inertia matrix for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:59
% EndTime: 2019-03-09 03:55:04
% DurationCPUTime: 2.18s
% Computational Cost: add. (4910->271), mult. (9592->405), div. (0->0), fcn. (9546->8), ass. (0->133)
t111 = sin(qJ(5));
t192 = cos(qJ(5));
t109 = sin(pkin(10));
t112 = sin(qJ(3));
t114 = cos(qJ(3));
t165 = cos(pkin(10));
t86 = -t109 * t114 - t112 * t165;
t140 = t165 * t114;
t87 = -t109 * t112 + t140;
t128 = -t111 * t87 + t192 * t86;
t161 = qJD(3) * t112;
t78 = -qJD(3) * t140 + t109 * t161;
t79 = t86 * qJD(3);
t40 = qJD(5) * t128 + t111 * t78 + t192 * t79;
t60 = t111 * t86 + t192 * t87;
t185 = t40 * t60;
t102 = pkin(3) * t165 + pkin(4);
t190 = pkin(3) * t109;
t74 = t192 * t102 - t111 * t190;
t67 = t74 * qJD(5);
t182 = t128 * t67;
t75 = t111 * t102 + t192 * t190;
t68 = t75 * qJD(5);
t183 = t68 * t60;
t208 = -t74 * t40 + t182 + t183;
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t93 = -mrSges(7,1) * t113 + mrSges(7,2) * t110;
t206 = -mrSges(6,1) + t93;
t118 = t60 * qJD(5) + t111 * t79 - t192 * t78;
t205 = t118 * t128;
t162 = t110 ^ 2 + t113 ^ 2;
t141 = t162 * t118;
t158 = qJD(6) * t113;
t127 = t110 * t40 + t60 * t158;
t160 = qJD(3) * t114;
t207 = t67 * mrSges(6,2);
t142 = t67 * t162;
t204 = -t109 * t78 + t165 * t79;
t159 = qJD(6) * t110;
t171 = t110 * t60;
t43 = mrSges(7,2) * t128 - mrSges(7,3) * t171;
t167 = t113 * t60;
t44 = -mrSges(7,1) * t128 - mrSges(7,3) * t167;
t203 = -t44 * t158 - t43 * t159;
t202 = -mrSges(4,1) * t161 - mrSges(4,2) * t160;
t130 = -t110 * t44 + t113 * t43;
t115 = -pkin(1) - pkin(7);
t163 = qJ(4) - t115;
t91 = t163 * t112;
t92 = t163 * t114;
t61 = t109 * t91 - t165 * t92;
t129 = -pkin(8) * t87 + t61;
t62 = -t109 * t92 - t165 * t91;
t56 = pkin(8) * t86 + t62;
t24 = t111 * t129 + t192 * t56;
t103 = t112 * pkin(3) + qJ(2);
t66 = -pkin(4) * t86 + t103;
t25 = -pkin(5) * t128 - pkin(9) * t60 + t66;
t15 = -t110 * t24 + t113 * t25;
t16 = t110 * t25 + t113 * t24;
t122 = -t114 * qJD(4) + t161 * t163;
t123 = -qJD(3) * t92 - t112 * qJD(4);
t54 = -t109 * t123 + t122 * t165;
t116 = -t79 * pkin(8) + t54;
t23 = t111 * t56 - t192 * t129;
t55 = t109 * t122 + t165 * t123;
t46 = pkin(8) * t78 + t55;
t11 = -t23 * qJD(5) + t111 * t116 + t192 * t46;
t97 = pkin(3) * t160 + qJD(2);
t63 = -pkin(4) * t78 + t97;
t17 = pkin(5) * t118 - pkin(9) * t40 + t63;
t3 = -qJD(6) * t16 - t11 * t110 + t113 * t17;
t189 = t110 * t3;
t201 = -t15 * t158 - t16 * t159 - t189;
t200 = 2 * m(6);
t199 = 2 * m(7);
t12 = t24 * qJD(5) + t111 * t46 - t192 * t116;
t198 = 0.2e1 * t12;
t197 = 0.2e1 * t23;
t196 = 0.2e1 * t63;
t195 = m(5) * pkin(3);
t2 = qJD(6) * t15 + t11 * t113 + t110 * t17;
t188 = t113 * t2;
t187 = t12 * t23;
t186 = t23 * t68;
t184 = t118 * mrSges(6,3);
t181 = t78 * t86;
t180 = t79 * t87;
t169 = t113 * t40;
t179 = Ifges(7,5) * t169 + Ifges(7,3) * t118;
t178 = Ifges(7,4) * t110;
t177 = Ifges(7,4) * t113;
t176 = Ifges(7,6) * t110;
t153 = t60 * t159;
t126 = t153 - t169;
t13 = mrSges(7,1) * t118 + mrSges(7,3) * t126;
t174 = t110 * t13;
t72 = pkin(9) + t75;
t166 = t113 * t72;
t157 = 2 * mrSges(5,3);
t148 = -t78 * mrSges(5,1) + t79 * mrSges(5,2);
t147 = mrSges(6,1) * t118 + t40 * mrSges(6,2);
t144 = -(2 * Ifges(6,4)) - t176;
t139 = -t12 * t60 - t23 * t40;
t137 = t180 + t181;
t136 = mrSges(7,1) * t110 + mrSges(7,2) * t113;
t135 = Ifges(7,1) * t113 - t178;
t134 = -Ifges(7,2) * t110 + t177;
t133 = Ifges(7,5) * t110 + Ifges(7,6) * t113;
t14 = -mrSges(7,2) * t118 - mrSges(7,3) * t127;
t132 = t113 * t14 - t174;
t131 = -t110 * t15 + t113 * t16;
t89 = t134 * qJD(6);
t90 = t135 * qJD(6);
t94 = Ifges(7,2) * t113 + t178;
t95 = Ifges(7,1) * t110 + t177;
t125 = t110 * t90 + t113 * t89 + t95 * t158 - t159 * t94;
t88 = t136 * qJD(6);
t124 = -mrSges(6,2) * t118 + mrSges(7,3) * t141 - t206 * t40 - t60 * t88;
t121 = t54 * t87 - t55 * t86 + t61 * t79 - t62 * t78;
t120 = -t189 + (-t110 * t16 - t113 * t15) * qJD(6);
t119 = t120 + t188;
t104 = Ifges(7,5) * t158;
t20 = -Ifges(7,6) * t128 + t134 * t60;
t21 = -Ifges(7,5) * t128 + t135 * t60;
t6 = -Ifges(7,4) * t126 - Ifges(7,2) * t127 + Ifges(7,6) * t118;
t7 = -Ifges(7,1) * t126 - Ifges(7,4) * t127 + Ifges(7,5) * t118;
t117 = -t11 * mrSges(6,2) + mrSges(7,3) * t188 + t21 * t158 / 0.2e1 + t23 * t88 + t95 * t169 / 0.2e1 + Ifges(6,5) * t40 + t110 * t7 / 0.2e1 - t89 * t171 / 0.2e1 + t90 * t167 / 0.2e1 + t113 * t6 / 0.2e1 - t128 * (-Ifges(7,6) * t159 + t104) / 0.2e1 - t127 * t94 / 0.2e1 - (t60 * t95 + t20) * t159 / 0.2e1 + t206 * t12 + (t133 / 0.2e1 - Ifges(6,6)) * t118;
t71 = -pkin(5) - t74;
t33 = t136 * t60;
t9 = mrSges(7,1) * t127 - mrSges(7,2) * t126;
t1 = [0.2e1 * t97 * (-mrSges(5,1) * t86 + mrSges(5,2) * t87) + 0.2e1 * t103 * t148 + 0.2e1 * Ifges(5,2) * t181 + 0.2e1 * Ifges(5,1) * t180 + 0.2e1 * t66 * t147 + 0.2e1 * t2 * t43 + 0.2e1 * t3 * t44 + t33 * t198 + t9 * t197 + 0.2e1 * t15 * t13 + 0.2e1 * t16 * t14 - 0.2e1 * t24 * t184 + (mrSges(6,3) * t197 - t110 * t20 + t113 * t21) * t40 + (t15 * t3 + t16 * t2 + t187) * t199 + (t11 * t24 + t63 * t66 + t187) * t200 + 0.2e1 * m(5) * (t103 * t97 + t54 * t61 + t55 * t62) + 0.2e1 * (t78 * t87 + t79 * t86) * Ifges(5,4) - t121 * t157 - (mrSges(6,1) * t196 - 0.2e1 * mrSges(6,3) * t11 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t118 + t144 * t40 + t179) * t128 + (mrSges(6,2) * t196 + mrSges(6,3) * t198 + 0.2e1 * Ifges(6,1) * t40 - t110 * t6 + t113 * t7 + (Ifges(7,5) * t113 + t144) * t118 + (-t110 * t21 - t113 * t20 + t128 * t133) * qJD(6)) * t60 + 0.2e1 * (qJ(2) * (mrSges(4,1) * t114 - mrSges(4,2) * t112) + (t112 ^ 2 - t114 ^ 2) * Ifges(4,4)) * qJD(3) + 0.2e1 * (mrSges(4,1) * t112 + mrSges(4,2) * t114 + mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * qJD(2) + 0.2e1 * (Ifges(4,2) - Ifges(4,1)) * t112 * t160; -t40 * t33 - t60 * t9 + t130 * t118 - t137 * t157 - ((-t110 * t43 - t113 * t44) * qJD(6) + t132) * t128 + m(7) * (t118 * t131 - t119 * t128 + t139) + m(6) * (-t11 * t128 + t118 * t24 + t139) + m(5) * t121 + (-0.2e1 * t185 + 0.2e1 * t205) * mrSges(6,3); 0.2e1 * m(6) * (t185 - t205) + 0.2e1 * m(7) * (-t128 * t141 + t185) + 0.2e1 * m(5) * t137; t117 + m(7) * (t12 * t71 + t186) + m(6) * (t11 * t75 - t12 * t74 + t186) - t75 * t184 - Ifges(4,6) * t160 - Ifges(4,5) * t161 + (t109 * t55 + t165 * t54) * t195 + t14 * t166 + Ifges(5,6) * t78 + Ifges(5,5) * t79 + t68 * t33 + t71 * t9 + t54 * mrSges(5,1) - t55 * mrSges(5,2) + (m(6) * t24 + m(7) * t131 + t130) * t67 - t204 * mrSges(5,3) * pkin(3) + t202 * t115 + (m(7) * t119 - t174 + t203) * t72 + t201 * mrSges(7,3) + t208 * mrSges(6,3); m(6) * (t118 * t75 - t208) + m(7) * (-t40 * t71 - t183 + t162 * (t118 * t72 - t182)) + t204 * t195 + t78 * mrSges(5,2) + t79 * mrSges(5,1) + t124 + t202; 0.2e1 * t71 * t88 - 0.2e1 * t207 + (t142 * t72 + t68 * t71) * t199 + (t67 * t75 - t68 * t74) * t200 + t125 + 0.2e1 * t206 * t68 + 0.2e1 * mrSges(7,3) * t142; t110 * t14 + t113 * t13 + t130 * qJD(6) + m(7) * (qJD(6) * t131 + t110 * t2 + t113 * t3) + m(6) * t63 + m(5) * t97 + t147 + t148; 0; 0; 0; t117 + (-m(7) * t12 - t9) * pkin(5) + (m(7) * (t188 + t201) + t132 + t203) * pkin(9) + t120 * mrSges(7,3); m(7) * (pkin(5) * t40 + pkin(9) * t141) + t124; -t207 + (t71 - pkin(5)) * t88 + t125 + (m(7) * pkin(9) + mrSges(7,3)) * t142 + (-m(7) * pkin(5) + t206) * t68; 0; -0.2e1 * pkin(5) * t88 + t125; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t153 - Ifges(7,6) * t127 + t179; (-t113 * t118 - t128 * t159) * mrSges(7,2) + (-t110 * t118 + t128 * t158) * mrSges(7,1); t104 - t136 * t67 + (-mrSges(7,1) * t166 + (mrSges(7,2) * t72 - Ifges(7,6)) * t110) * qJD(6); -t88; t104 + (pkin(9) * t93 - t176) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
