% Calculate time derivative of joint inertia matrix for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:17
% EndTime: 2019-03-08 23:00:24
% DurationCPUTime: 3.16s
% Computational Cost: add. (5536->356), mult. (13353->542), div. (0->0), fcn. (13260->12), ass. (0->161)
t127 = sin(qJ(6));
t131 = cos(qJ(6));
t199 = mrSges(7,1) * t131;
t110 = mrSges(7,2) * t127 - t199;
t220 = -mrSges(6,1) + t110;
t172 = qJD(6) * t131;
t124 = sin(pkin(12));
t184 = cos(pkin(12));
t128 = sin(qJ(4));
t129 = sin(qJ(3));
t132 = cos(qJ(4));
t133 = cos(qJ(3));
t102 = -t128 * t129 + t132 * t133;
t217 = qJD(3) + qJD(4);
t79 = t217 * t102;
t103 = t128 * t133 + t129 * t132;
t80 = t217 * t103;
t58 = -t124 * t80 + t184 * t79;
t76 = t102 * t124 + t103 * t184;
t145 = t127 * t58 + t172 * t76;
t210 = -pkin(9) - pkin(8);
t169 = t210 * t129;
t219 = t210 * t133;
t222 = -t128 * t169 + t132 * t219;
t84 = t128 * t219 + t132 * t169;
t75 = -t102 * t184 + t103 * t124;
t118 = -pkin(3) * t133 - pkin(2);
t86 = -pkin(4) * t102 + t118;
t40 = pkin(5) * t75 - pkin(10) * t76 + t86;
t138 = -t103 * qJ(5) + t84;
t70 = qJ(5) * t102 - t222;
t42 = t124 * t138 + t184 * t70;
t21 = t127 * t40 + t131 * t42;
t177 = t21 * qJD(6);
t139 = t222 * qJD(3);
t146 = -t79 * qJ(5) - t103 * qJD(5);
t154 = qJD(4) * t219;
t155 = qJD(4) * t169;
t61 = qJD(3) * t84 + t128 * t154 + t132 * t155;
t35 = -t80 * qJ(5) + t102 * qJD(5) + t61;
t19 = t184 * t35 + (-t128 * t155 + t132 * t154 + t139 + t146) * t124;
t57 = t124 * t79 + t184 * t80;
t120 = qJD(3) * t129 * pkin(3);
t73 = pkin(4) * t80 + t120;
t24 = pkin(5) * t57 - pkin(10) * t58 + t73;
t3 = -t127 * t19 + t131 * t24 - t177;
t221 = -t3 - t177;
t180 = t124 * t128;
t193 = pkin(3) * qJD(4);
t90 = (t132 * t184 - t180) * t193;
t158 = (t127 ^ 2 + t131 ^ 2) * t90;
t125 = sin(pkin(6));
t134 = cos(qJ(2));
t178 = t125 * t134;
t126 = cos(pkin(6));
t130 = sin(qJ(2));
t179 = t125 * t130;
t94 = t126 * t133 - t129 * t179;
t95 = t126 * t129 + t133 * t179;
t71 = -t128 * t95 + t132 * t94;
t72 = t128 * t94 + t132 * t95;
t44 = t124 * t71 + t184 * t72;
t143 = t127 * t178 - t131 * t44;
t37 = -t127 * t44 - t131 * t178;
t218 = -t127 * t37 - t131 * t143;
t216 = (-mrSges(5,1) * t128 - mrSges(5,2) * t132) * t193;
t215 = 2 * m(6);
t214 = 2 * m(7);
t62 = qJD(4) * t222 + t139;
t18 = t124 * t35 - t184 * (t146 + t62);
t213 = 0.2e1 * t18;
t41 = t124 * t70 - t138 * t184;
t212 = 0.2e1 * t41;
t211 = m(6) * pkin(4);
t208 = pkin(4) * t124;
t207 = t18 * t41;
t157 = t184 * t128;
t89 = (t124 * t132 + t157) * t193;
t206 = t41 * t89;
t174 = qJD(2) * t134;
t161 = t125 * t174;
t82 = -qJD(3) * t95 - t129 * t161;
t83 = qJD(3) * t94 + t133 * t161;
t33 = qJD(4) * t71 + t128 * t82 + t132 * t83;
t34 = -qJD(4) * t72 - t128 * t83 + t132 * t82;
t15 = t124 * t33 - t184 * t34;
t43 = t124 * t72 - t184 * t71;
t205 = t43 * t15;
t204 = t43 * t89;
t203 = t57 * mrSges(6,3);
t187 = t131 * t58;
t201 = Ifges(7,5) * t187 + Ifges(7,3) * t57;
t197 = mrSges(7,3) * t131;
t196 = Ifges(7,4) * t127;
t195 = Ifges(7,4) * t131;
t194 = Ifges(7,6) * t127;
t192 = t127 * mrSges(7,3);
t189 = t127 * t76;
t186 = t131 * t76;
t150 = mrSges(7,1) * t127 + mrSges(7,2) * t131;
t106 = t150 * qJD(6);
t117 = pkin(3) * t132 + pkin(4);
t91 = -pkin(3) * t180 + t117 * t184;
t87 = -pkin(5) - t91;
t185 = t87 * t106;
t20 = -t127 * t42 + t131 * t40;
t183 = qJD(6) * t20;
t92 = pkin(3) * t157 + t117 * t124;
t88 = pkin(10) + t92;
t182 = qJD(6) * t88;
t163 = t184 * pkin(4);
t116 = -t163 - pkin(5);
t181 = t116 * t106;
t175 = qJD(2) * t130;
t173 = qJD(6) * t127;
t171 = 0.2e1 * t129;
t167 = mrSges(7,3) * t173;
t166 = t76 * t173;
t165 = mrSges(7,3) * t172;
t162 = t125 * t175;
t25 = mrSges(6,1) * t57 + mrSges(6,2) * t58;
t160 = -t173 / 0.2e1;
t159 = -(2 * Ifges(6,4)) - t194;
t153 = t125 ^ 2 * t130 * t174;
t152 = t41 * t15 + t18 * t43;
t151 = -mrSges(4,1) * t133 + mrSges(4,2) * t129;
t149 = Ifges(7,1) * t131 - t196;
t148 = -Ifges(7,2) * t127 + t195;
t147 = Ifges(7,5) * t127 + Ifges(7,6) * t131;
t144 = t166 - t187;
t108 = t148 * qJD(6);
t109 = t149 * qJD(6);
t111 = Ifges(7,2) * t131 + t196;
t112 = Ifges(7,1) * t127 + t195;
t142 = t108 * t131 + t109 * t127 - t111 * t173 + t112 * t172;
t16 = t124 * t34 + t184 * t33;
t5 = qJD(6) * t37 + t127 * t162 + t131 * t16;
t141 = t34 * mrSges(5,1) - t33 * mrSges(5,2) - t16 * mrSges(6,2) + t43 * t106 + t15 * t220 + t5 * t197;
t6 = qJD(6) * t143 - t127 * t16 + t131 * t162;
t140 = -t6 * t127 + (t127 * t143 - t131 * t37) * qJD(6);
t137 = t131 * t5 + t140;
t136 = -t82 * t129 + t83 * t133 + (-t129 * t95 - t133 * t94) * qJD(3);
t10 = -Ifges(7,1) * t144 - Ifges(7,4) * t145 + Ifges(7,5) * t57;
t119 = Ifges(7,5) * t172;
t2 = t127 * t24 + t131 * t19 + t183;
t28 = Ifges(7,6) * t75 + t148 * t76;
t29 = Ifges(7,5) * t75 + t149 * t76;
t9 = -Ifges(7,4) * t144 - Ifges(7,2) * t145 + Ifges(7,6) * t57;
t135 = -t61 * mrSges(5,2) - t19 * mrSges(6,2) + t2 * t197 + t28 * t160 + t29 * t172 / 0.2e1 + t41 * t106 + Ifges(6,5) * t58 + t62 * mrSges(5,1) - t108 * t189 / 0.2e1 + t109 * t186 / 0.2e1 + t75 * (-Ifges(7,6) * t173 + t119) / 0.2e1 + t127 * t10 / 0.2e1 - Ifges(5,6) * t80 + Ifges(5,5) * t79 + t131 * t9 / 0.2e1 + (t147 / 0.2e1 - Ifges(6,6)) * t57 - t145 * t111 / 0.2e1 + t220 * t18 + (t187 / 0.2e1 + t76 * t160) * t112;
t115 = pkin(10) + t208;
t107 = (mrSges(4,1) * t129 + mrSges(4,2) * t133) * qJD(3);
t81 = -mrSges(5,1) * t102 + mrSges(5,2) * t103;
t59 = mrSges(5,1) * t80 + mrSges(5,2) * t79;
t52 = mrSges(6,1) * t75 + mrSges(6,2) * t76;
t50 = mrSges(7,1) * t75 - mrSges(7,3) * t186;
t49 = -mrSges(7,2) * t75 - mrSges(7,3) * t189;
t48 = t150 * t76;
t23 = -mrSges(7,2) * t57 - mrSges(7,3) * t145;
t22 = mrSges(7,1) * t57 + mrSges(7,3) * t144;
t13 = mrSges(7,1) * t145 - mrSges(7,2) * t144;
t1 = [0.2e1 * m(7) * (-t143 * t5 + t37 * t6 + t205) + 0.2e1 * m(6) * (t16 * t44 - t153 + t205) + 0.2e1 * m(5) * (t33 * t72 + t34 * t71 - t153) + 0.2e1 * m(4) * (t82 * t94 + t83 * t95 - t153); t43 * t13 + t15 * t48 + t37 * t22 - t143 * t23 + t5 * t49 + t6 * t50 + (t15 * t76 - t16 * t75 + t43 * t58 - t44 * t57) * mrSges(6,3) + (t102 * t33 - t103 * t34 - t71 * t79 - t72 * t80) * mrSges(5,3) + t136 * mrSges(4,3) + ((-t107 - t25 - t59) * t134 + (-t134 * mrSges(3,2) + (-mrSges(3,1) + t151 + t52 + t81) * t130) * qJD(2)) * t125 + m(5) * (-t222 * t33 + t84 * t34 + t61 * t72 + t62 * t71 + (t118 * t175 - t120 * t134) * t125) + m(6) * (t42 * t16 + t19 * t44 + (-t134 * t73 + t175 * t86) * t125 + t152) + m(7) * (-t143 * t2 + t20 * t6 + t21 * t5 + t3 * t37 + t152) + (-pkin(2) * t162 + pkin(8) * t136) * m(4); -0.2e1 * t42 * t203 + 0.2e1 * t79 * t103 * Ifges(5,1) - 0.2e1 * t102 * Ifges(5,2) * t80 - 0.2e1 * pkin(2) * t107 + 0.2e1 * t118 * t59 + t13 * t212 + t48 * t213 + 0.2e1 * t2 * t49 + 0.2e1 * t20 * t22 + 0.2e1 * t21 * t23 + 0.2e1 * t86 * t25 + 0.2e1 * t3 * t50 + 0.2e1 * t73 * t52 + (mrSges(6,3) * t212 - t127 * t28 + t131 * t29) * t58 + 0.2e1 * m(5) * (t118 * t120 - t222 * t61 + t62 * t84) + (t19 * t42 + t73 * t86 + t207) * t215 + (t2 * t21 + t20 * t3 + t207) * t214 + (-0.2e1 * t19 * mrSges(6,3) + t159 * t58 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t57 + t201) * t75 + (mrSges(6,3) * t213 + 0.2e1 * Ifges(6,1) * t58 + t131 * t10 - t127 * t9 + (Ifges(7,5) * t131 + t159) * t57 + (-t127 * t29 - t131 * t28 - t147 * t75) * qJD(6)) * t76 + ((-Ifges(4,4) * t129 + pkin(3) * t81) * t171 + (0.2e1 * Ifges(4,4) * t133 + (Ifges(4,1) - Ifges(4,2)) * t171) * t133) * qJD(3) + 0.2e1 * (t102 * t79 - t103 * t80) * Ifges(5,4) + 0.2e1 * (t102 * t61 - t103 * t62 + t222 * t80 - t79 * t84) * mrSges(5,3); t82 * mrSges(4,1) - t83 * mrSges(4,2) + m(6) * (-t15 * t91 + t16 * t92 + t44 * t90 + t204) + t140 * mrSges(7,3) + m(5) * (t128 * t33 + t132 * t34 + (-t128 * t71 + t132 * t72) * qJD(4)) * pkin(3) + t141 + (t137 * t88 + t15 * t87 + t218 * t90 + t204) * m(7); (-t57 * t92 - t58 * t91 - t75 * t90 + t76 * t89) * mrSges(6,3) + t87 * t13 + t89 * t48 + t135 + m(6) * (-t18 * t91 + t19 * t92 + t42 * t90 + t206) + (m(7) * (-t182 * t20 + t2 * t88 + t21 * t90) + t88 * t23 + t90 * t49 - mrSges(7,3) * t183 - t50 * t182) * t131 + (-t49 * t182 + t221 * mrSges(7,3) + (-m(7) * t20 - t50) * t90 + (m(7) * t221 - t22) * t88) * t127 + (m(5) * (t128 * t61 + t132 * t62 + (-t128 * t84 - t132 * t222) * qJD(4)) + (-t128 * t80 - t132 * t79 + (t102 * t132 + t103 * t128) * qJD(4)) * mrSges(5,3)) * pkin(3) + m(7) * (t18 * t87 + t206) + (Ifges(4,5) * t133 - Ifges(4,6) * t129 + pkin(8) * t151) * qJD(3); 0.2e1 * t185 - 0.2e1 * t90 * mrSges(6,2) + (t158 * t88 + t87 * t89) * t214 + (-t89 * t91 + t90 * t92) * t215 + t142 + 0.2e1 * t220 * t89 + 0.2e1 * t216 + 0.2e1 * mrSges(7,3) * t158; m(7) * (t115 * t137 + t116 * t15) + t143 * t167 - t6 * t192 - t37 * t165 + (t124 * t16 - t15 * t184) * t211 + t141; (t124 * t19 - t18 * t184) * t211 + t135 - t58 * mrSges(6,3) * t163 - t3 * t192 - t203 * t208 - t21 * t167 - t20 * t165 + (m(7) * t18 + t13) * t116 + (t131 * t23 - t127 * t22 + m(7) * (-t127 * t3 + t131 * t2 + (-t127 * t21 - t131 * t20) * qJD(6)) - t49 * t173 - t50 * t172) * t115; t142 + t181 + t185 + (t124 * t211 - mrSges(6,2)) * t90 + (m(7) * t116 - t184 * t211 + t220) * t89 + t216 + (m(7) * t115 + mrSges(7,3)) * t158; t142 + 0.2e1 * t181; m(7) * (qJD(6) * t218 + t127 * t5 + t131 * t6) + m(6) * t162; t127 * t23 + t131 * t22 + (-t127 * t50 + t131 * t49) * qJD(6) + m(7) * (t127 * t2 + t131 * t3 + (-t127 * t20 + t131 * t21) * qJD(6)) + m(6) * t73 + t25; 0; 0; 0; mrSges(7,1) * t6 - mrSges(7,2) * t5; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t166 - Ifges(7,6) * t145 + t201; t119 - t150 * t90 + (-t88 * t199 + (mrSges(7,2) * t88 - Ifges(7,6)) * t127) * qJD(6); t119 + (t110 * t115 - t194) * qJD(6); -t106; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
