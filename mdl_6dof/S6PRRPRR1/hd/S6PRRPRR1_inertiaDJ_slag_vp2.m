% Calculate time derivative of joint inertia matrix for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:51:17
% EndTime: 2019-03-08 21:51:23
% DurationCPUTime: 2.43s
% Computational Cost: add. (5231->335), mult. (12377->510), div. (0->0), fcn. (12737->12), ass. (0->149)
t125 = sin(qJ(6));
t129 = cos(qJ(6));
t184 = mrSges(7,1) * t129;
t107 = mrSges(7,2) * t125 - t184;
t202 = -mrSges(6,1) + t107;
t162 = qJD(6) * t129;
t126 = sin(qJ(5));
t193 = cos(qJ(5));
t121 = sin(pkin(12));
t127 = sin(qJ(3));
t123 = cos(pkin(12));
t130 = cos(qJ(3));
t167 = t123 * t130;
t98 = -t121 * t127 + t167;
t99 = t121 * t130 + t123 * t127;
t143 = -t126 * t99 + t193 * t98;
t91 = t99 * qJD(3);
t92 = t98 * qJD(3);
t48 = qJD(5) * t143 - t126 * t91 + t193 * t92;
t72 = t126 * t98 + t193 * t99;
t142 = t125 * t48 + t72 * t162;
t114 = pkin(3) * t123 + pkin(4);
t192 = pkin(3) * t121;
t88 = t114 * t193 - t126 * t192;
t82 = t88 * qJD(5);
t204 = t82 * mrSges(6,2);
t49 = qJD(5) * t72 + t126 * t92 + t193 * t91;
t25 = t49 * mrSges(6,1) + t48 * mrSges(6,2);
t70 = t91 * mrSges(5,1) + t92 * mrSges(5,2);
t203 = -t25 - t70;
t153 = (t125 ^ 2 + t129 ^ 2) * t82;
t89 = t126 * t114 + t193 * t192;
t186 = -qJ(4) - pkin(8);
t108 = t186 * t127;
t109 = t186 * t130;
t76 = t123 * t108 + t109 * t121;
t145 = -pkin(9) * t99 + t76;
t100 = t121 * t108;
t77 = -t123 * t109 + t100;
t67 = pkin(9) * t98 + t77;
t201 = -t126 * t67 + t193 * t145;
t122 = sin(pkin(6));
t131 = cos(qJ(2));
t168 = t122 * t131;
t124 = cos(pkin(6));
t128 = sin(qJ(2));
t169 = t122 * t128;
t93 = t124 * t130 - t127 * t169;
t94 = t124 * t127 + t130 * t169;
t68 = -t121 * t94 + t123 * t93;
t69 = t121 * t93 + t123 * t94;
t40 = t126 * t68 + t193 * t69;
t140 = t125 * t168 - t129 * t40;
t31 = -t125 * t40 - t129 * t168;
t200 = -t125 * t31 - t129 * t140;
t199 = 2 * m(6);
t198 = 2 * m(7);
t197 = -2 * mrSges(6,3);
t65 = -t99 * qJD(4) + (t167 * t186 - t100) * qJD(3);
t132 = -t92 * pkin(9) + t65;
t36 = t126 * t145 + t193 * t67;
t154 = qJD(3) * t186;
t66 = t123 * (qJD(4) * t130 + t127 * t154) + t121 * (-qJD(4) * t127 + t130 * t154);
t57 = -pkin(9) * t91 + t66;
t19 = qJD(5) * t36 + t126 * t57 - t132 * t193;
t196 = 0.2e1 * t19;
t195 = -0.2e1 * t201;
t191 = t19 * t201;
t115 = -pkin(3) * t130 - pkin(2);
t79 = -pkin(4) * t98 + t115;
t34 = -pkin(5) * t143 - pkin(10) * t72 + t79;
t23 = t125 * t34 + t129 * t36;
t171 = qJD(6) * t23;
t18 = qJD(5) * t201 + t126 * t132 + t193 * t57;
t117 = qJD(3) * t127 * pkin(3);
t78 = pkin(4) * t91 + t117;
t24 = pkin(5) * t49 - pkin(10) * t48 + t78;
t3 = -t125 * t18 + t129 * t24 - t171;
t190 = t3 * t125;
t83 = t89 * qJD(5);
t189 = t201 * t83;
t144 = -t126 * t69 + t193 * t68;
t164 = qJD(2) * t131;
t157 = t122 * t164;
t74 = -qJD(3) * t94 - t127 * t157;
t75 = qJD(3) * t93 + t130 * t157;
t54 = -t121 * t75 + t123 * t74;
t55 = t121 * t74 + t123 * t75;
t16 = qJD(5) * t40 + t126 * t55 - t193 * t54;
t188 = t144 * t16;
t187 = t144 * t83;
t174 = t129 * t48;
t185 = Ifges(7,5) * t174 + Ifges(7,3) * t49;
t183 = mrSges(7,3) * t129;
t182 = Ifges(7,4) * t125;
t181 = Ifges(7,4) * t129;
t180 = Ifges(7,6) * t125;
t177 = t125 * t72;
t173 = t129 * t72;
t22 = -t125 * t36 + t129 * t34;
t172 = qJD(6) * t22;
t86 = pkin(10) + t89;
t170 = qJD(6) * t86;
t165 = qJD(2) * t128;
t163 = qJD(6) * t125;
t161 = 0.2e1 * t127;
t160 = t72 * t163;
t158 = t122 * t165;
t156 = -t163 / 0.2e1;
t155 = -(2 * Ifges(6,4)) - t180;
t152 = t122 ^ 2 * t128 * t164;
t151 = -t144 * t19 - t16 * t201;
t150 = -mrSges(4,1) * t130 + mrSges(4,2) * t127;
t149 = mrSges(7,1) * t125 + mrSges(7,2) * t129;
t148 = Ifges(7,1) * t129 - t182;
t147 = -Ifges(7,2) * t125 + t181;
t146 = Ifges(7,5) * t125 + Ifges(7,6) * t129;
t141 = t160 - t174;
t105 = t147 * qJD(6);
t106 = t148 * qJD(6);
t110 = Ifges(7,2) * t129 + t182;
t111 = Ifges(7,1) * t125 + t181;
t139 = t129 * t105 + t125 * t106 - t110 * t163 + t111 * t162;
t15 = qJD(5) * t144 + t126 * t54 + t193 * t55;
t6 = qJD(6) * t140 - t125 * t15 + t129 * t158;
t137 = -t6 * t125 + (t125 * t140 - t129 * t31) * qJD(6);
t5 = qJD(6) * t31 + t125 * t158 + t129 * t15;
t136 = t129 * t5 + t137;
t135 = -t74 * t127 + t75 * t130 + (-t127 * t94 - t130 * t93) * qJD(3);
t103 = t149 * qJD(6);
t134 = -t15 * mrSges(6,2) + mrSges(7,3) * t137 - t144 * t103 + t202 * t16 + t5 * t183;
t10 = -Ifges(7,1) * t141 - Ifges(7,4) * t142 + t49 * Ifges(7,5);
t116 = Ifges(7,5) * t162;
t2 = t125 * t24 + t129 * t18 + t172;
t28 = -Ifges(7,6) * t143 + t147 * t72;
t29 = -Ifges(7,5) * t143 + t148 * t72;
t9 = -Ifges(7,4) * t141 - Ifges(7,2) * t142 + t49 * Ifges(7,6);
t133 = -t18 * mrSges(6,2) + t2 * t183 + t28 * t156 + t29 * t162 / 0.2e1 - t201 * t103 + Ifges(6,5) * t48 - t105 * t177 / 0.2e1 + t106 * t173 / 0.2e1 - t143 * (-Ifges(7,6) * t163 + t116) / 0.2e1 + t125 * t10 / 0.2e1 + t129 * t9 / 0.2e1 + (t146 / 0.2e1 - Ifges(6,6)) * t49 - t142 * t110 / 0.2e1 + t202 * t19 + (t174 / 0.2e1 + t72 * t156) * t111;
t104 = (mrSges(4,1) * t127 + mrSges(4,2) * t130) * qJD(3);
t85 = -pkin(5) - t88;
t73 = -mrSges(5,1) * t98 + mrSges(5,2) * t99;
t52 = -mrSges(6,1) * t143 + mrSges(6,2) * t72;
t51 = -mrSges(7,1) * t143 - mrSges(7,3) * t173;
t50 = mrSges(7,2) * t143 - mrSges(7,3) * t177;
t43 = t149 * t72;
t21 = -mrSges(7,2) * t49 - mrSges(7,3) * t142;
t20 = mrSges(7,1) * t49 + mrSges(7,3) * t141;
t11 = mrSges(7,1) * t142 - mrSges(7,2) * t141;
t1 = [0.2e1 * m(7) * (-t140 * t5 + t31 * t6 - t188) + 0.2e1 * m(6) * (t40 * t15 - t152 - t188) + 0.2e1 * m(5) * (t68 * t54 + t69 * t55 - t152) + 0.2e1 * m(4) * (t93 * t74 + t94 * t75 - t152); -t144 * t11 + t16 * t43 + t31 * t20 - t140 * t21 + t5 * t50 + t6 * t51 + (t143 * t15 - t144 * t48 + t16 * t72 - t40 * t49) * mrSges(6,3) + (-t54 * t99 + t55 * t98 - t68 * t92 - t69 * t91) * mrSges(5,3) + t135 * mrSges(4,3) + ((-t104 + t203) * t131 + (-t131 * mrSges(3,2) + (-mrSges(3,1) + t150 + t52 + t73) * t128) * qJD(2)) * t122 + m(5) * (t76 * t54 + t77 * t55 + t65 * t68 + t66 * t69 + (t115 * t165 - t117 * t131) * t122) + m(6) * (t36 * t15 + t18 * t40 + (-t131 * t78 + t165 * t79) * t122 + t151) + m(7) * (-t140 * t2 + t22 * t6 + t23 * t5 + t3 * t31 + t151) + (-pkin(2) * t158 + pkin(8) * t135) * m(4); t36 * t49 * t197 + 0.2e1 * t92 * t99 * Ifges(5,1) - 0.2e1 * t98 * Ifges(5,2) * t91 - 0.2e1 * pkin(2) * t104 + t11 * t195 + 0.2e1 * t115 * t70 + t43 * t196 + 0.2e1 * t2 * t50 + 0.2e1 * t22 * t20 + 0.2e1 * t23 * t21 + 0.2e1 * t79 * t25 + 0.2e1 * t3 * t51 + 0.2e1 * t78 * t52 + (mrSges(6,3) * t195 - t125 * t28 + t129 * t29) * t48 + 0.2e1 * m(5) * (t115 * t117 + t65 * t76 + t66 * t77) + (t18 * t36 + t78 * t79 - t191) * t199 + (t2 * t23 + t22 * t3 - t191) * t198 - (t18 * t197 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t49 + t155 * t48 + t185) * t143 + (mrSges(6,3) * t196 + 0.2e1 * Ifges(6,1) * t48 + t129 * t10 - t125 * t9 + (Ifges(7,5) * t129 + t155) * t49 + (-t125 * t29 - t129 * t28 + t143 * t146) * qJD(6)) * t72 + ((-Ifges(4,4) * t127 + pkin(3) * t73) * t161 + (0.2e1 * Ifges(4,4) * t130 + (Ifges(4,1) - Ifges(4,2)) * t161) * t130) * qJD(3) + 0.2e1 * (-t91 * t99 + t92 * t98) * Ifges(5,4) + 0.2e1 * (-t65 * t99 + t66 * t98 - t76 * t92 - t77 * t91) * mrSges(5,3); t74 * mrSges(4,1) + t54 * mrSges(5,1) - t75 * mrSges(4,2) - t55 * mrSges(5,2) + m(6) * (t15 * t89 - t16 * t88 + t40 * t82 - t187) + m(5) * (t121 * t55 + t123 * t54) * pkin(3) + t134 + (t136 * t86 + t16 * t85 + t200 * t82 - t187) * m(7); m(7) * (t19 * t85 - t189) + (Ifges(4,5) * t130 - Ifges(4,6) * t127 + pkin(8) * t150) * qJD(3) + m(6) * (t18 * t89 - t19 * t88 + t36 * t82 - t189) + t85 * t11 - Ifges(5,6) * t91 + Ifges(5,5) * t92 + t83 * t43 + t65 * mrSges(5,1) - t66 * mrSges(5,2) + (t143 * t82 - t48 * t88 - t49 * t89 + t72 * t83) * mrSges(6,3) + t133 + (m(7) * (-t170 * t23 - t22 * t82 - t3 * t86) - t86 * t20 - t82 * t51 - t50 * t170 + (-t3 - t171) * mrSges(7,3)) * t125 + (m(5) * (t121 * t66 + t123 * t65) + (-t121 * t91 - t123 * t92) * mrSges(5,3)) * pkin(3) + (m(7) * (-t170 * t22 + t2 * t86 + t23 * t82) + t86 * t21 + t82 * t50 - mrSges(7,3) * t172 - t51 * t170) * t129; 0.2e1 * t85 * t103 - 0.2e1 * t204 + (t153 * t86 + t83 * t85) * t198 + (t82 * t89 - t83 * t88) * t199 + t139 + 0.2e1 * t202 * t83 + 0.2e1 * mrSges(7,3) * t153; m(7) * (qJD(6) * t200 + t125 * t5 + t129 * t6) + (m(5) + m(6)) * t158; m(5) * t117 + t125 * t21 + t129 * t20 + (-t125 * t51 + t129 * t50) * qJD(6) + m(7) * (t125 * t2 + t129 * t3 + (-t125 * t22 + t129 * t23) * qJD(6)) + m(6) * t78 - t203; 0; 0; m(7) * (-pkin(5) * t16 + pkin(10) * t136) + t134; t133 + (-t190 + (-t125 * t23 - t129 * t22) * qJD(6)) * mrSges(7,3) + (-m(7) * t19 - t11) * pkin(5) + (m(7) * (t129 * t2 - t162 * t22 - t163 * t23 - t190) + t129 * t21 - t125 * t20 - t50 * t163 - t51 * t162) * pkin(10); -t204 + (-pkin(5) + t85) * t103 + t139 + (m(7) * pkin(10) + mrSges(7,3)) * t153 + (-m(7) * pkin(5) + t202) * t83; 0; -0.2e1 * pkin(5) * t103 + t139; mrSges(7,1) * t6 - mrSges(7,2) * t5; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t160 - Ifges(7,6) * t142 + t185; t116 - t149 * t82 + (-t86 * t184 + (mrSges(7,2) * t86 - Ifges(7,6)) * t125) * qJD(6); -t103; t116 + (pkin(10) * t107 - t180) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
