% Calculate time derivative of joint inertia matrix for
% S6PRRPRR5
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:26
% EndTime: 2019-03-08 22:16:31
% DurationCPUTime: 2.84s
% Computational Cost: add. (5202->444), mult. (13196->680), div. (0->0), fcn. (12788->12), ass. (0->185)
t163 = sin(pkin(12));
t212 = pkin(9) + qJ(4);
t145 = t212 * t163;
t165 = cos(pkin(12));
t147 = t212 * t165;
t168 = sin(qJ(5));
t172 = cos(qJ(5));
t106 = -t168 * t145 + t172 * t147;
t177 = t163 * t168 - t165 * t172;
t123 = t177 * qJD(5);
t173 = cos(qJ(3));
t192 = qJD(3) * t173;
t181 = t163 * t192;
t210 = mrSges(5,2) * t165;
t118 = mrSges(5,1) * t181 + t192 * t210;
t139 = t163 * t172 + t165 * t168;
t124 = t139 * qJD(5);
t169 = sin(qJ(3));
t79 = -t124 * t169 - t177 * t192;
t80 = t123 * t169 - t139 * t192;
t42 = -t80 * mrSges(6,1) + t79 * mrSges(6,2);
t167 = sin(qJ(6));
t171 = cos(qJ(6));
t116 = t139 * t169;
t117 = t177 * t169;
t69 = -t116 * t171 + t117 * t167;
t31 = qJD(6) * t69 + t167 * t80 + t171 * t79;
t70 = -t116 * t167 - t117 * t171;
t32 = -qJD(6) * t70 - t167 * t79 + t171 * t80;
t9 = -t32 * mrSges(7,1) + t31 * mrSges(7,2);
t227 = t118 + t42 + t9;
t193 = qJD(3) * t169;
t226 = -Ifges(6,5) * t79 - Ifges(6,6) * t80 - Ifges(6,3) * t193;
t225 = 2 * m(5);
t224 = 2 * m(6);
t223 = 2 * m(7);
t222 = 2 * pkin(8);
t162 = t165 ^ 2;
t221 = m(5) / 0.2e1;
t220 = m(7) * pkin(5);
t90 = -t139 * t167 - t171 * t177;
t219 = t90 / 0.2e1;
t91 = t139 * t171 - t167 * t177;
t218 = t91 / 0.2e1;
t217 = -t177 / 0.2e1;
t216 = t139 / 0.2e1;
t215 = t165 / 0.2e1;
t213 = pkin(5) * t124;
t46 = qJD(6) * t90 - t123 * t171 - t124 * t167;
t47 = -qJD(6) * t91 + t123 * t167 - t124 * t171;
t211 = Ifges(7,5) * t46 + Ifges(7,6) * t47;
t144 = -pkin(3) * t173 - qJ(4) * t169 - pkin(2);
t197 = t165 * t173;
t112 = pkin(8) * t197 + t163 * t144;
t202 = t163 * t169;
t104 = -pkin(9) * t202 + t112;
t134 = t165 * t144;
t198 = t165 * t169;
t92 = -pkin(9) * t198 + t134 + (-pkin(8) * t163 - pkin(4)) * t173;
t52 = t172 * t104 + t168 * t92;
t209 = Ifges(5,4) * t163;
t208 = Ifges(5,4) * t165;
t207 = t163 * Ifges(5,2);
t206 = -mrSges(5,1) * t165 + mrSges(5,2) * t163 - mrSges(4,1);
t166 = cos(pkin(6));
t164 = sin(pkin(6));
t170 = sin(qJ(2));
t200 = t164 * t170;
t126 = t166 * t169 + t173 * t200;
t174 = cos(qJ(2));
t194 = qJD(2) * t174;
t182 = t164 * t194;
t102 = qJD(3) * t126 + t169 * t182;
t205 = t102 * t169;
t125 = -t166 * t173 + t169 * t200;
t103 = -qJD(3) * t125 + t173 * t182;
t204 = t103 * t173;
t203 = t104 * t168;
t68 = t125 * t102;
t201 = t163 * t173;
t199 = t164 * t174;
t196 = -Ifges(6,5) * t123 - Ifges(6,6) * t124;
t122 = -qJD(4) * t169 + (pkin(3) * t169 - qJ(4) * t173) * qJD(3);
t185 = pkin(8) * t193;
t98 = t165 * t122 + t163 * t185;
t159 = pkin(8) * t192;
t130 = pkin(4) * t181 + t159;
t143 = pkin(4) * t202 + t169 * pkin(8);
t191 = qJD(4) * t163;
t190 = qJD(4) * t165;
t189 = qJD(5) * t172;
t188 = qJD(6) * t167;
t187 = qJD(6) * t171;
t186 = -Ifges(7,5) * t31 - Ifges(7,6) * t32 - Ifges(7,3) * t193;
t100 = -t126 * t163 - t165 * t199;
t101 = t126 * t165 - t163 * t199;
t53 = t100 * t172 - t101 * t168;
t183 = qJD(2) * t200;
t72 = -t103 * t163 + t165 * t183;
t73 = t103 * t165 + t163 * t183;
t17 = qJD(5) * t53 + t168 * t72 + t172 * t73;
t54 = t100 * t168 + t101 * t172;
t18 = -qJD(5) * t54 - t168 * t73 + t172 * t72;
t26 = -t167 * t54 + t171 * t53;
t5 = qJD(6) * t26 + t167 * t18 + t17 * t171;
t27 = t167 * t53 + t171 * t54;
t6 = -qJD(6) * t27 - t167 * t17 + t171 * t18;
t184 = t6 * mrSges(7,1) - t5 * mrSges(7,2);
t156 = -pkin(4) * t165 - pkin(3);
t21 = -mrSges(7,1) * t47 + t46 * mrSges(7,2);
t51 = t172 * t92 - t203;
t105 = -t172 * t145 - t147 * t168;
t180 = -Ifges(5,5) * t165 + Ifges(5,6) * t163;
t179 = -t163 * t72 + t165 * t73;
t35 = -pkin(5) * t173 + pkin(10) * t117 + t51;
t41 = -pkin(10) * t116 + t52;
t15 = -t167 * t41 + t171 * t35;
t16 = t167 * t35 + t171 * t41;
t77 = -pkin(10) * t139 + t105;
t78 = -pkin(10) * t177 + t106;
t39 = -t167 * t78 + t171 * t77;
t40 = t167 * t77 + t171 * t78;
t64 = -t145 * t189 + t172 * t190 + (-qJD(5) * t147 - t191) * t168;
t55 = -pkin(10) * t124 + t64;
t65 = -t139 * qJD(4) - qJD(5) * t106;
t56 = pkin(10) * t123 + t65;
t11 = qJD(6) * t39 + t167 * t56 + t171 * t55;
t12 = -qJD(6) * t40 - t167 * t55 + t171 * t56;
t178 = t12 * mrSges(7,1) - t11 * mrSges(7,2) + t211;
t71 = (pkin(4) * t169 - pkin(9) * t197) * qJD(3) + t98;
t114 = t163 * t122;
t83 = t114 + (-pkin(8) * t198 - pkin(9) * t201) * qJD(3);
t20 = -qJD(5) * t52 - t168 * t83 + t172 * t71;
t13 = pkin(5) * t193 - pkin(10) * t79 + t20;
t19 = -qJD(5) * t203 + t168 * t71 + t172 * t83 + t92 * t189;
t14 = pkin(10) * t80 + t19;
t2 = qJD(6) * t15 + t13 * t167 + t14 * t171;
t3 = -qJD(6) * t16 + t13 * t171 - t14 * t167;
t176 = t3 * mrSges(7,1) - t2 * mrSges(7,2) - t186;
t175 = t125 * t192 + t205;
t142 = (mrSges(4,1) * t169 + mrSges(4,2) * t173) * qJD(3);
t141 = -mrSges(5,1) * t173 - mrSges(5,3) * t198;
t140 = mrSges(5,2) * t173 - mrSges(5,3) * t202;
t132 = (-mrSges(7,1) * t167 - mrSges(7,2) * t171) * qJD(6) * pkin(5);
t129 = (mrSges(5,1) * t169 - mrSges(5,3) * t197) * qJD(3);
t128 = (-mrSges(5,2) * t169 - mrSges(5,3) * t201) * qJD(3);
t127 = (mrSges(5,1) * t163 + t210) * t169;
t119 = t123 * mrSges(6,2);
t113 = pkin(5) * t177 + t156;
t111 = -pkin(8) * t201 + t134;
t110 = (t169 * Ifges(5,5) + (t165 * Ifges(5,1) - t209) * t173) * qJD(3);
t109 = (t169 * Ifges(5,6) + (-t207 + t208) * t173) * qJD(3);
t108 = -mrSges(6,1) * t173 + mrSges(6,3) * t117;
t107 = mrSges(6,2) * t173 - mrSges(6,3) * t116;
t99 = -t165 * t185 + t114;
t97 = Ifges(6,1) * t139 - Ifges(6,4) * t177;
t96 = Ifges(6,4) * t139 - Ifges(6,2) * t177;
t95 = mrSges(6,1) * t177 + mrSges(6,2) * t139;
t94 = pkin(5) * t116 + t143;
t89 = -Ifges(6,1) * t123 - Ifges(6,4) * t124;
t88 = -Ifges(6,4) * t123 - Ifges(6,2) * t124;
t87 = mrSges(6,1) * t124 - t119;
t82 = mrSges(6,1) * t116 - mrSges(6,2) * t117;
t67 = -Ifges(6,1) * t117 - Ifges(6,4) * t116 - Ifges(6,5) * t173;
t66 = -Ifges(6,4) * t117 - Ifges(6,2) * t116 - Ifges(6,6) * t173;
t61 = -mrSges(6,2) * t193 + mrSges(6,3) * t80;
t60 = mrSges(6,1) * t193 - mrSges(6,3) * t79;
t59 = -mrSges(7,1) * t173 - mrSges(7,3) * t70;
t58 = mrSges(7,2) * t173 + mrSges(7,3) * t69;
t57 = -pkin(5) * t80 + t130;
t50 = Ifges(7,1) * t91 + Ifges(7,4) * t90;
t49 = Ifges(7,4) * t91 + Ifges(7,2) * t90;
t48 = -mrSges(7,1) * t90 + mrSges(7,2) * t91;
t38 = -mrSges(7,1) * t69 + mrSges(7,2) * t70;
t37 = Ifges(6,1) * t79 + Ifges(6,4) * t80 + Ifges(6,5) * t193;
t36 = Ifges(6,4) * t79 + Ifges(6,2) * t80 + Ifges(6,6) * t193;
t34 = Ifges(7,1) * t70 + Ifges(7,4) * t69 - Ifges(7,5) * t173;
t33 = Ifges(7,4) * t70 + Ifges(7,2) * t69 - Ifges(7,6) * t173;
t25 = -mrSges(7,2) * t193 + mrSges(7,3) * t32;
t24 = mrSges(7,1) * t193 - mrSges(7,3) * t31;
t23 = Ifges(7,1) * t46 + Ifges(7,4) * t47;
t22 = Ifges(7,4) * t46 + Ifges(7,2) * t47;
t8 = Ifges(7,1) * t31 + Ifges(7,4) * t32 + Ifges(7,5) * t193;
t7 = Ifges(7,4) * t31 + Ifges(7,2) * t32 + Ifges(7,6) * t193;
t1 = [0.2e1 * m(7) * (t26 * t6 + t27 * t5 + t68) + 0.2e1 * m(6) * (t17 * t54 + t18 * t53 + t68) + 0.2e1 * m(5) * (t100 * t72 + t101 * t73 + t68) + 0.2e1 * m(4) * (-t164 ^ 2 * t170 * t194 + t126 * t103 + t68); t100 * t129 + t101 * t128 + t17 * t107 + t18 * t108 + t73 * t140 + t72 * t141 + t26 * t24 + t27 * t25 + t5 * t58 + t53 * t60 + t54 * t61 + t6 * t59 + t227 * t125 + (t127 + t82 + t38) * t102 + (-t174 * t142 + (-t174 * mrSges(3,2) + (-t173 * mrSges(4,1) + t169 * mrSges(4,2) - mrSges(3,1)) * t170) * qJD(2)) * t164 + (t205 + t204 + (t125 * t173 - t126 * t169) * qJD(3)) * mrSges(4,3) + m(5) * (t100 * t98 + t101 * t99 + t111 * t72 + t112 * t73) - m(4) * pkin(2) * t183 + m(7) * (t102 * t94 + t125 * t57 + t15 * t6 + t16 * t5 + t2 * t27 + t26 * t3) + m(6) * (t102 * t143 + t125 * t130 + t17 * t52 + t18 * t51 + t19 * t54 + t20 * t53) + (t175 * t221 + m(4) * (-t126 * t193 + t175 + t204) / 0.2e1) * t222; (t130 * t143 + t19 * t52 + t20 * t51) * t224 + (t111 * t98 + t112 * t99) * t225 + (t186 + t226) * t173 + ((-Ifges(6,5) * t117 + Ifges(7,5) * t70 - Ifges(6,6) * t116 + Ifges(7,6) * t69 + (-(2 * Ifges(4,4)) - t180) * t169) * t169 + (t127 * t222 + 0.2e1 * (Ifges(4,4) + t180) * t173 + (Ifges(5,1) * t162 - (2 * Ifges(5,3)) - (2 * Ifges(4,2)) + (2 * Ifges(4,1)) + (pkin(8) ^ 2 * t225) - Ifges(6,3) - Ifges(7,3) + (t207 - 0.2e1 * t208) * t163) * t169) * t173) * qJD(3) + (t15 * t3 + t16 * t2 + t57 * t94) * t223 + 0.2e1 * t99 * t140 + 0.2e1 * t98 * t141 - 0.2e1 * pkin(2) * t142 + 0.2e1 * t143 * t42 + 0.2e1 * t112 * t128 + 0.2e1 * t111 * t129 + 0.2e1 * t130 * t82 - t116 * t36 - t117 * t37 + 0.2e1 * t19 * t107 + 0.2e1 * t20 * t108 + 0.2e1 * t94 * t9 + t79 * t67 + t80 * t66 + t69 * t7 + t70 * t8 + 0.2e1 * t57 * t38 + 0.2e1 * t2 * t58 + 0.2e1 * t3 * t59 + 0.2e1 * t51 * t60 + 0.2e1 * t52 * t61 + t32 * t33 + t31 * t34 + 0.2e1 * t15 * t24 + 0.2e1 * t16 * t25 + (-t109 * t163 + t110 * t165 + t118 * t222) * t169; -t103 * mrSges(4,2) + (t21 + t87) * t125 + t179 * mrSges(5,3) + (t48 + t95 + t206) * t102 + m(7) * (t102 * t113 + t11 * t27 + t12 * t26 + t125 * t213 + t39 * t6 + t40 * t5) + m(6) * (t102 * t156 + t105 * t18 + t106 * t17 + t53 * t65 + t54 * t64) + m(5) * (-pkin(3) * t102 + (-t100 * t163 + t101 * t165) * qJD(4) + t179 * qJ(4)) + (-t26 * t46 + t27 * t47 + t5 * t90 - t6 * t91) * mrSges(7,3) + (t123 * t53 - t124 * t54 - t139 * t18 - t17 * t177) * mrSges(6,3); -(t211 + t196) * t173 / 0.2e1 + m(5) * (-t111 * t191 + t112 * t190 + (-t163 * t98 + t165 * t99) * qJ(4)) + (-t15 * t46 + t16 * t47 + t2 * t90 - t3 * t91) * mrSges(7,3) + m(7) * (t11 * t16 + t113 * t57 + t12 * t15 + t2 * t40 + t213 * t94 + t3 * t39) + t37 * t216 + t36 * t217 + t8 * t218 + t7 * t219 + (t109 / 0.2e1 + qJ(4) * t128 + qJD(4) * t140 + t99 * mrSges(5,3)) * t165 + (t110 / 0.2e1 - qJ(4) * t129 - qJD(4) * t141 - t98 * mrSges(5,3)) * t163 - (t66 / 0.2e1 - pkin(5) * t38) * t124 + (t123 * t51 - t124 * t52 - t139 * t20 - t177 * t19) * mrSges(6,3) + t156 * t42 + t143 * t87 + t130 * t95 - t116 * t88 / 0.2e1 - t117 * t89 / 0.2e1 - pkin(3) * t118 - t123 * t67 / 0.2e1 + t106 * t61 + t64 * t107 + t65 * t108 + t113 * t9 + t94 * t21 + t80 * t96 / 0.2e1 + t79 * t97 / 0.2e1 + t105 * t60 + t69 * t22 / 0.2e1 + t70 * t23 / 0.2e1 + t57 * t48 + t11 * t58 + t12 * t59 + t47 * t33 / 0.2e1 + t32 * t49 / 0.2e1 + t31 * t50 / 0.2e1 + t39 * t24 + t40 * t25 + t46 * t34 / 0.2e1 + ((pkin(8) * mrSges(4,2) + Ifges(5,5) * t163 / 0.2e1 + Ifges(5,6) * t215 + Ifges(6,5) * t216 + Ifges(6,6) * t217 + Ifges(7,5) * t218 + Ifges(7,6) * t219 - Ifges(4,6)) * t169 + (-t163 * (Ifges(5,2) * t165 + t209) / 0.2e1 + (Ifges(5,1) * t163 + t208) * t215 + Ifges(4,5) + (-m(5) * pkin(3) + t206) * pkin(8)) * t173) * qJD(3) + m(6) * (t105 * t20 + t106 * t19 + t130 * t156 + t51 * t65 + t52 * t64); 0.2e1 * t113 * t21 - t123 * t97 - t177 * t88 + t139 * t89 + 0.2e1 * t156 * t87 + t90 * t22 + t91 * t23 + t46 * t50 + t47 * t49 - (-0.2e1 * pkin(5) * t48 + t96) * t124 + (t11 * t40 + t113 * t213 + t12 * t39) * t223 + (t105 * t65 + t106 * t64) * t224 + 0.2e1 * (t11 * t90 - t12 * t91 - t39 * t46 + t40 * t47) * mrSges(7,3) + 0.2e1 * (t105 * t123 - t106 * t124 - t139 * t65 - t177 * t64) * mrSges(6,3) + (qJ(4) * t225 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t163 ^ 2 + t162); 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1 + t221) * t102; m(5) * t159 + m(6) * t130 + m(7) * t57 + t227; -t119 - (-mrSges(6,1) - t220) * t124 + t21; 0; t18 * mrSges(6,1) - t17 * mrSges(6,2) + (t167 * t5 + t171 * t6 + (-t167 * t26 + t171 * t27) * qJD(6)) * t220 + t184; t20 * mrSges(6,1) - t19 * mrSges(6,2) + (t58 * t187 + t167 * t25 + m(7) * (-t15 * t188 + t16 * t187 + t167 * t2 + t171 * t3) - t59 * t188 + t171 * t24) * pkin(5) + t176 - t226; t65 * mrSges(6,1) - t64 * mrSges(6,2) + (m(7) * (t11 * t167 + t12 * t171 + (-t167 * t39 + t171 * t40) * qJD(6)) + (t167 * t47 - t171 * t46 + (t167 * t91 + t171 * t90) * qJD(6)) * mrSges(7,3)) * pkin(5) + t178 + t196; 0; 0.2e1 * t132; t184; t176; t178; 0; t132; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
