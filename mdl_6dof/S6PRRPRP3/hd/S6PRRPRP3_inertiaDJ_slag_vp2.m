% Calculate time derivative of joint inertia matrix for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:16
% EndTime: 2019-03-08 21:34:21
% DurationCPUTime: 2.71s
% Computational Cost: add. (2596->384), mult. (6919->569), div. (0->0), fcn. (6209->10), ass. (0->163)
t213 = Ifges(7,4) + Ifges(6,5);
t211 = Ifges(6,6) - Ifges(7,6);
t212 = -Ifges(7,2) - Ifges(6,3);
t210 = mrSges(6,3) + mrSges(7,2);
t209 = m(7) * qJD(6);
t140 = cos(pkin(11));
t145 = cos(qJ(5));
t138 = sin(pkin(11));
t142 = sin(qJ(5));
t180 = t138 * t142;
t113 = -t145 * t140 + t180;
t102 = t113 * qJD(5);
t114 = t138 * t145 + t140 * t142;
t103 = t114 * qJD(5);
t208 = -t213 * t102 - t211 * t103;
t143 = sin(qJ(3));
t168 = qJD(5) * t143;
t146 = cos(qJ(3));
t170 = qJD(3) * t146;
t49 = -t113 * t170 - t114 * t168;
t167 = qJD(5) * t145;
t175 = t140 * t143;
t50 = t114 * t170 + t167 * t175 - t168 * t180;
t15 = t50 * mrSges(7,1) - t49 * mrSges(7,3);
t16 = t50 * mrSges(6,1) + t49 * mrSges(6,2);
t161 = t138 * t170;
t188 = mrSges(5,2) * t140;
t94 = mrSges(5,1) * t161 + t170 * t188;
t207 = t15 + t16 + t94;
t206 = m(7) * qJ(6) + mrSges(7,3);
t119 = -pkin(3) * t146 - qJ(4) * t143 - pkin(2);
t111 = t140 * t119;
t64 = -pkin(9) * t175 + t111 + (-pkin(8) * t138 - pkin(4)) * t146;
t179 = t138 * t143;
t174 = t140 * t146;
t89 = pkin(8) * t174 + t138 * t119;
t79 = -pkin(9) * t179 + t89;
t192 = t142 * t64 + t145 * t79;
t101 = -qJD(4) * t143 + (pkin(3) * t143 - qJ(4) * t146) * qJD(3);
t171 = qJD(3) * t143;
t164 = pkin(8) * t171;
t73 = t140 * t101 + t138 * t164;
t39 = (pkin(4) * t143 - pkin(9) * t174) * qJD(3) + t73;
t178 = t138 * t146;
t90 = t138 * t101;
t54 = t90 + (-pkin(8) * t175 - pkin(9) * t178) * qJD(3);
t8 = -qJD(5) * t192 - t142 * t54 + t145 * t39;
t205 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t204 = -mrSges(6,2) + t206;
t203 = 2 * m(5);
t202 = 2 * m(6);
t201 = 0.2e1 * m(7);
t200 = 0.2e1 * pkin(8);
t137 = t140 ^ 2;
t199 = m(5) / 0.2e1;
t198 = t140 / 0.2e1;
t196 = pkin(9) + qJ(4);
t26 = -mrSges(7,2) * t50 + mrSges(7,3) * t171;
t29 = -mrSges(6,2) * t171 - mrSges(6,3) * t50;
t195 = t26 + t29;
t27 = mrSges(6,1) * t171 - mrSges(6,3) * t49;
t28 = -mrSges(7,1) * t171 + t49 * mrSges(7,2);
t194 = -t27 + t28;
t57 = t103 * mrSges(7,1) + t102 * mrSges(7,3);
t58 = t103 * mrSges(6,1) - t102 * mrSges(6,2);
t193 = t57 + t58;
t92 = t114 * t143;
t82 = -mrSges(7,2) * t92 - mrSges(7,3) * t146;
t83 = mrSges(6,2) * t146 - mrSges(6,3) * t92;
t191 = t82 + t83;
t93 = t113 * t143;
t84 = -mrSges(6,1) * t146 + mrSges(6,3) * t93;
t85 = mrSges(7,1) * t146 - mrSges(7,2) * t93;
t190 = -t84 + t85;
t187 = Ifges(5,4) * t138;
t186 = Ifges(5,4) * t140;
t141 = cos(pkin(6));
t139 = sin(pkin(6));
t144 = sin(qJ(2));
t177 = t139 * t144;
t104 = -t141 * t146 + t143 * t177;
t105 = t141 * t143 + t146 * t177;
t147 = cos(qJ(2));
t172 = qJD(2) * t147;
t162 = t139 * t172;
t77 = qJD(3) * t105 + t143 * t162;
t38 = t104 * t77;
t185 = t138 * Ifges(5,2);
t184 = t77 * t143;
t78 = -qJD(3) * t104 + t146 * t162;
t183 = t78 * t146;
t181 = -mrSges(5,1) * t140 + mrSges(5,2) * t138 - mrSges(4,1);
t176 = t139 * t147;
t134 = pkin(8) * t170;
t109 = pkin(4) * t161 + t134;
t118 = pkin(4) * t179 + t143 * pkin(8);
t169 = qJD(5) * t142;
t166 = t142 * qJD(4);
t165 = t145 * qJD(4);
t131 = -pkin(4) * t140 - pkin(3);
t163 = qJD(2) * t177;
t122 = t196 * t140;
t158 = qJD(5) * t196;
t32 = t140 * t165 - t122 * t169 + (-t145 * t158 - t166) * t138;
t33 = t140 * t166 + t122 * t167 + (-t142 * t158 + t165) * t138;
t159 = t196 * t138;
t80 = t142 * t122 + t145 * t159;
t81 = t145 * t122 - t142 * t159;
t160 = t81 * t32 + t33 * t80;
t156 = -Ifges(5,5) * t140 + Ifges(5,6) * t138;
t40 = -t138 * t78 + t140 * t163;
t41 = t138 * t163 + t140 * t78;
t155 = -t138 * t40 + t140 * t41;
t20 = -t142 * t79 + t145 * t64;
t75 = -t105 * t138 - t140 * t176;
t76 = t105 * t140 - t138 * t176;
t22 = t142 * t76 - t145 * t75;
t23 = t142 * t75 + t145 * t76;
t152 = t212 * t171 + t211 * t50 - t213 * t49;
t5 = -qJD(5) * t22 + t142 * t40 + t145 * t41;
t6 = qJD(5) * t23 + t142 * t41 - t145 * t40;
t151 = t22 * t33 + t32 * t23 + t81 * t5 + t6 * t80;
t150 = t104 * t170 + t184;
t7 = t142 * t39 + t145 * t54 + t64 * t167 - t169 * t79;
t117 = (mrSges(4,1) * t143 + mrSges(4,2) * t146) * qJD(3);
t116 = -mrSges(5,1) * t146 - mrSges(5,3) * t175;
t115 = mrSges(5,2) * t146 - mrSges(5,3) * t179;
t108 = (mrSges(5,1) * t143 - mrSges(5,3) * t174) * qJD(3);
t107 = (-mrSges(5,2) * t143 - mrSges(5,3) * t178) * qJD(3);
t106 = (mrSges(5,1) * t138 + t188) * t143;
t88 = -pkin(8) * t178 + t111;
t87 = (t143 * Ifges(5,5) + (t140 * Ifges(5,1) - t187) * t146) * qJD(3);
t86 = (t143 * Ifges(5,6) + (-t185 + t186) * t146) * qJD(3);
t74 = -t140 * t164 + t90;
t72 = Ifges(6,1) * t114 - Ifges(6,4) * t113;
t71 = Ifges(7,1) * t114 + Ifges(7,5) * t113;
t70 = Ifges(6,4) * t114 - Ifges(6,2) * t113;
t69 = Ifges(7,5) * t114 + Ifges(7,3) * t113;
t68 = mrSges(6,1) * t113 + mrSges(6,2) * t114;
t67 = mrSges(7,1) * t113 - mrSges(7,3) * t114;
t63 = pkin(5) * t113 - qJ(6) * t114 + t131;
t62 = -Ifges(6,1) * t102 - Ifges(6,4) * t103;
t61 = -Ifges(7,1) * t102 + Ifges(7,5) * t103;
t60 = -Ifges(6,4) * t102 - Ifges(6,2) * t103;
t59 = -Ifges(7,5) * t102 + Ifges(7,3) * t103;
t53 = mrSges(6,1) * t92 - mrSges(6,2) * t93;
t52 = mrSges(7,1) * t92 + mrSges(7,3) * t93;
t37 = -Ifges(6,1) * t93 - Ifges(6,4) * t92 - Ifges(6,5) * t146;
t36 = -Ifges(7,1) * t93 - Ifges(7,4) * t146 + Ifges(7,5) * t92;
t35 = -Ifges(6,4) * t93 - Ifges(6,2) * t92 - Ifges(6,6) * t146;
t34 = -Ifges(7,5) * t93 - Ifges(7,6) * t146 + Ifges(7,3) * t92;
t25 = pkin(5) * t103 + qJ(6) * t102 - qJD(6) * t114;
t24 = pkin(5) * t92 + qJ(6) * t93 + t118;
t19 = pkin(5) * t146 - t20;
t17 = -qJ(6) * t146 + t192;
t14 = Ifges(6,1) * t49 - Ifges(6,4) * t50 + Ifges(6,5) * t171;
t13 = Ifges(7,1) * t49 + Ifges(7,4) * t171 + Ifges(7,5) * t50;
t12 = Ifges(6,4) * t49 - Ifges(6,2) * t50 + Ifges(6,6) * t171;
t11 = Ifges(7,5) * t49 + Ifges(7,6) * t171 + Ifges(7,3) * t50;
t9 = pkin(5) * t50 - qJ(6) * t49 + qJD(6) * t93 + t109;
t4 = -pkin(5) * t171 - t8;
t3 = qJ(6) * t171 - qJD(6) * t146 + t7;
t1 = [0.2e1 * m(5) * (t40 * t75 + t41 * t76 + t38) + 0.2e1 * m(4) * (-t139 ^ 2 * t144 * t172 + t105 * t78 + t38) + 0.2e1 * (m(7) + m(6)) * (t22 * t6 + t23 * t5 + t38); t76 * t107 + t75 * t108 + t41 * t115 + t40 * t116 + t190 * t6 + t191 * t5 + t195 * t23 + t194 * t22 + (t52 + t53 + t106) * t77 + t207 * t104 + (-t147 * t117 + (-t147 * mrSges(3,2) + (-t146 * mrSges(4,1) + t143 * mrSges(4,2) - mrSges(3,1)) * t144) * qJD(2)) * t139 + (t184 + t183 + (t104 * t146 - t105 * t143) * qJD(3)) * mrSges(4,3) + m(7) * (t104 * t9 + t17 * t5 + t19 * t6 + t22 * t4 + t23 * t3 + t24 * t77) + m(6) * (t104 * t109 + t118 * t77 + t192 * t5 - t20 * t6 - t22 * t8 + t23 * t7) + m(5) * (t40 * t88 + t41 * t89 + t73 * t75 + t74 * t76) - m(4) * pkin(2) * t163 + (t150 * t199 + m(4) * (-t105 * t171 + t150 + t183) / 0.2e1) * t200; (t36 + t37) * t49 + (-t138 * t86 + t140 * t87 + t94 * t200) * t143 + t152 * t146 + (t17 * t3 + t19 * t4 + t24 * t9) * t201 + (t73 * t88 + t74 * t89) * t203 + (((-(2 * Ifges(4,4)) - t156) * t143 - t213 * t93 - t211 * t92) * t143 + (t106 * t200 + 0.2e1 * (Ifges(4,4) + t156) * t146 + (Ifges(5,1) * t137 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + pkin(8) ^ 2 * t203 - (2 * Ifges(5,3)) + (t185 - 0.2e1 * t186) * t138 + t212) * t143) * t146) * qJD(3) - (t13 + t14) * t93 + (t109 * t118 + t192 * t7 + t20 * t8) * t202 + 0.2e1 * t192 * t29 + (t11 - t12) * t92 + (t34 - t35) * t50 + 0.2e1 * t24 * t15 + 0.2e1 * t17 * t26 + 0.2e1 * t20 * t27 + 0.2e1 * t19 * t28 + 0.2e1 * t9 * t52 + 0.2e1 * t3 * t82 + 0.2e1 * t7 * t83 + 0.2e1 * t8 * t84 + 0.2e1 * t4 * t85 + 0.2e1 * t89 * t107 + 0.2e1 * t88 * t108 + 0.2e1 * t109 * t53 + 0.2e1 * t74 * t115 + 0.2e1 * t73 * t116 - 0.2e1 * pkin(2) * t117 + 0.2e1 * t118 * t16; -t78 * mrSges(4,2) + t193 * t104 + t155 * mrSges(5,3) + (t67 + t68 + t181) * t77 + m(7) * (t104 * t25 + t63 * t77 + t151) + m(6) * (t131 * t77 + t151) + m(5) * (-pkin(3) * t77 + (-t138 * t75 + t140 * t76) * qJD(4) + t155 * qJ(4)) + t210 * (-t102 * t22 - t103 * t23 - t113 * t5 + t114 * t6); (t13 / 0.2e1 + t14 / 0.2e1) * t114 + (t11 / 0.2e1 - t12 / 0.2e1) * t113 + (t74 * mrSges(5,3) + qJ(4) * t107 + qJD(4) * t115 + t86 / 0.2e1) * t140 + (-t73 * mrSges(5,3) - qJ(4) * t108 - qJD(4) * t116 + t87 / 0.2e1) * t138 + (t34 / 0.2e1 - t35 / 0.2e1) * t103 - (t36 / 0.2e1 + t37 / 0.2e1) * t102 + ((pkin(8) * mrSges(4,2) + Ifges(5,5) * t138 / 0.2e1 + Ifges(5,6) * t198 - Ifges(4,6) + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t114 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t113) * t143 + ((Ifges(5,1) * t138 + t186) * t198 - t138 * (Ifges(5,2) * t140 + t187) / 0.2e1 + Ifges(4,5) + (-m(5) * pkin(3) + t181) * pkin(8)) * t146) * qJD(3) + t190 * t33 + t191 * t32 + t194 * t80 + t195 * t81 + m(7) * (t17 * t32 + t19 * t33 + t24 * t25 + t3 * t81 + t4 * t80 + t63 * t9) - (t61 / 0.2e1 + t62 / 0.2e1) * t93 - t208 * t146 / 0.2e1 + (t59 / 0.2e1 - t60 / 0.2e1) * t92 + (t69 / 0.2e1 - t70 / 0.2e1) * t50 + (t71 / 0.2e1 + t72 / 0.2e1) * t49 + (-t102 * t19 - t103 * t17 - t113 * t3 + t114 * t4) * mrSges(7,2) + (t102 * t20 - t103 * t192 - t113 * t7 - t114 * t8) * mrSges(6,3) + m(6) * (t109 * t131 + t192 * t32 - t20 * t33 + t7 * t81 - t8 * t80) + m(5) * ((-t138 * t88 + t140 * t89) * qJD(4) + (-t138 * t73 + t140 * t74) * qJ(4)) + t25 * t52 + t24 * t57 + t63 * t15 + t9 * t67 - pkin(3) * t94 + t109 * t68 + t118 * t58 + t131 * t16; 0.2e1 * t131 * t58 + 0.2e1 * t25 * t67 + 0.2e1 * t63 * t57 + (t61 + t62) * t114 + (-t60 + t59) * t113 + (-t70 + t69) * t103 - (t71 + t72) * t102 + t160 * t202 + (t25 * t63 + t160) * t201 + (qJ(4) * t203 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t138 ^ 2 + t137) + 0.2e1 * t210 * (-t102 * t80 - t103 * t81 - t113 * t32 + t114 * t33); 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1 + t199) * t77; m(5) * t134 + m(6) * t109 + m(7) * t9 + t207; m(7) * t25 + t193; 0; t204 * t5 + t205 * t6 + t23 * t209; qJD(6) * t82 + qJ(6) * t26 + m(7) * (-pkin(5) * t4 + qJ(6) * t3 + qJD(6) * t17) + t3 * mrSges(7,3) - t7 * mrSges(6,2) - pkin(5) * t28 - t4 * mrSges(7,1) + t8 * mrSges(6,1) - t152; t81 * t209 + (pkin(5) * t102 - qJ(6) * t103 - qJD(6) * t113) * mrSges(7,2) + t205 * t33 + t204 * t32 + t208; 0; 0.2e1 * t206 * qJD(6); m(7) * t6; m(7) * t4 + t28; m(7) * t33 - t102 * mrSges(7,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
