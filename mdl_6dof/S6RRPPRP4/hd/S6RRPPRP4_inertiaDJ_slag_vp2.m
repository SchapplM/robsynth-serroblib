% Calculate time derivative of joint inertia matrix for
% S6RRPPRP4
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:37:03
% EndTime: 2019-03-09 08:37:12
% DurationCPUTime: 3.61s
% Computational Cost: add. (2309->402), mult. (5540->566), div. (0->0), fcn. (4452->6), ass. (0->160)
t211 = Ifges(7,4) + Ifges(6,5);
t210 = Ifges(6,6) - Ifges(7,6);
t137 = sin(pkin(9));
t138 = cos(pkin(9));
t139 = sin(qJ(5));
t141 = cos(qJ(5));
t104 = t137 * t139 + t138 * t141;
t94 = t104 * qJD(5);
t159 = qJD(5) * t141;
t160 = qJD(5) * t139;
t95 = t137 * t159 - t138 * t160;
t209 = -t210 * t95 - t211 * t94;
t208 = qJD(3) * (t137 ^ 2 + t138 ^ 2);
t189 = -mrSges(7,1) - mrSges(6,1);
t207 = mrSges(6,3) + mrSges(7,2);
t206 = Ifges(5,1) + Ifges(4,1);
t205 = m(7) * qJD(6);
t142 = cos(qJ(2));
t140 = sin(qJ(2));
t165 = qJD(2) * t140;
t204 = qJ(4) * t165 - qJD(4) * t142;
t203 = m(7) * qJ(6) + mrSges(7,3);
t202 = pkin(7) + (pkin(3) + pkin(4)) * t137;
t111 = -pkin(2) * t142 - t140 * qJ(3) - pkin(1);
t171 = t137 * t142;
t125 = pkin(7) * t171;
t134 = t142 * pkin(3);
t48 = pkin(4) * t142 + t125 + t134 + (-pkin(8) * t140 - t111) * t138;
t173 = t137 * t140;
t169 = t138 * t142;
t82 = pkin(7) * t169 + t137 * t111;
t71 = -qJ(4) * t142 + t82;
t57 = pkin(8) * t173 + t71;
t184 = t139 * t48 + t141 * t57;
t156 = -pkin(7) * t137 - pkin(3);
t93 = -t140 * qJD(3) + (pkin(2) * t140 - qJ(3) * t142) * qJD(2);
t174 = t138 * t93;
t25 = -t174 + (-pkin(8) * t169 + (-pkin(4) + t156) * t140) * qJD(2);
t170 = t138 * t140;
t83 = t137 * t93;
t26 = t83 + (-pkin(7) * t170 + pkin(8) * t171) * qJD(2) + t204;
t4 = -qJD(5) * t184 - t139 * t26 + t141 * t25;
t201 = -m(7) * pkin(5) + t189;
t200 = -mrSges(6,2) + t203;
t199 = 2 * m(4);
t198 = 2 * m(5);
t197 = 2 * m(6);
t196 = 2 * m(7);
t195 = -0.2e1 * pkin(1);
t194 = 0.2e1 * pkin(7);
t155 = pkin(3) * t137 + pkin(7);
t164 = qJD(2) * t142;
t153 = t138 * t164;
t168 = -qJ(4) * t153 - qJD(4) * t170;
t59 = t155 * t164 + t168;
t192 = m(5) * t59;
t188 = mrSges(7,3) - mrSges(6,2);
t187 = -pkin(8) + qJ(3);
t154 = t137 * t164;
t86 = t104 * t140;
t41 = qJD(5) * t86 + t139 * t153 - t141 * t154;
t21 = -mrSges(7,2) * t41 - mrSges(7,3) * t165;
t22 = mrSges(6,2) * t165 - mrSges(6,3) * t41;
t186 = t21 + t22;
t172 = t137 * t141;
t105 = -t138 * t139 + t172;
t42 = qJD(5) * t105 * t140 + t104 * t164;
t23 = -mrSges(6,1) * t165 - mrSges(6,3) * t42;
t24 = mrSges(7,1) * t165 + t42 * mrSges(7,2);
t185 = t24 - t23;
t85 = t139 * t170 - t140 * t172;
t73 = -t85 * mrSges(7,2) + mrSges(7,3) * t142;
t74 = -mrSges(6,2) * t142 - t85 * mrSges(6,3);
t183 = t73 + t74;
t75 = mrSges(6,1) * t142 - t86 * mrSges(6,3);
t76 = -mrSges(7,1) * t142 + t86 * mrSges(7,2);
t182 = -t75 + t76;
t179 = mrSges(5,3) * t138;
t178 = Ifges(4,4) * t137;
t177 = Ifges(4,4) * t138;
t176 = Ifges(5,5) * t137;
t175 = Ifges(5,5) * t138;
t88 = mrSges(4,1) * t154 + mrSges(4,2) * t153;
t167 = qJ(3) * t208;
t163 = qJD(3) * t141;
t162 = qJD(4) * t137;
t158 = t139 * qJD(3);
t110 = -t138 * pkin(3) - t137 * qJ(4) - pkin(2);
t157 = pkin(7) * t165;
t114 = t187 * t138;
t28 = t138 * t163 - t114 * t160 + (t159 * t187 + t158) * t137;
t151 = t187 * t137;
t69 = t141 * t114 + t139 * t151;
t29 = qJD(5) * t69 - t137 * t163 + t138 * t158;
t68 = t114 * t139 - t141 * t151;
t152 = t69 * t28 + t29 * t68;
t150 = -t210 * t41 + t211 * t42;
t81 = t138 * t111 - t125;
t96 = t138 * pkin(4) - t110;
t100 = -mrSges(5,1) * t165 + mrSges(5,2) * t153;
t14 = -t139 * t57 + t141 * t48;
t67 = -t138 * t157 + t83;
t3 = t139 * t25 + t141 * t26 + t48 * t159 - t160 * t57;
t124 = qJ(4) * t170;
t70 = -t140 * t202 + t124;
t47 = t164 * t202 + t168;
t145 = (-Ifges(5,4) - Ifges(4,5)) * t138 + (Ifges(4,6) - Ifges(5,6)) * t137;
t117 = mrSges(5,1) * t154;
t113 = -mrSges(5,1) * t138 - mrSges(5,3) * t137;
t109 = -mrSges(5,2) * t173 - mrSges(5,3) * t142;
t108 = mrSges(5,1) * t142 + mrSges(5,2) * t170;
t107 = -mrSges(4,1) * t142 - mrSges(4,3) * t170;
t106 = mrSges(4,2) * t142 - mrSges(4,3) * t173;
t101 = (-mrSges(5,2) * t171 + mrSges(5,3) * t140) * qJD(2);
t99 = (mrSges(4,1) * t140 - mrSges(4,3) * t169) * qJD(2);
t98 = (-mrSges(4,2) * t140 - mrSges(4,3) * t171) * qJD(2);
t97 = (mrSges(5,1) * t137 - t179) * t140;
t87 = -mrSges(5,3) * t153 + t117;
t84 = t140 * t155 - t124;
t80 = (t140 * Ifges(4,5) + (t138 * Ifges(4,1) - t178) * t142) * qJD(2);
t79 = (t140 * Ifges(5,4) + (t138 * Ifges(5,1) + t176) * t142) * qJD(2);
t78 = (t140 * Ifges(4,6) + (-t137 * Ifges(4,2) + t177) * t142) * qJD(2);
t77 = (t140 * Ifges(5,6) + (t137 * Ifges(5,3) + t175) * t142) * qJD(2);
t72 = t134 - t81;
t66 = t137 * t157 + t174;
t65 = Ifges(6,1) * t105 - Ifges(6,4) * t104;
t64 = Ifges(7,1) * t105 + Ifges(7,5) * t104;
t63 = Ifges(6,4) * t105 - Ifges(6,2) * t104;
t62 = Ifges(7,5) * t105 + Ifges(7,3) * t104;
t61 = mrSges(6,1) * t104 + mrSges(6,2) * t105;
t60 = mrSges(7,1) * t104 - mrSges(7,3) * t105;
t56 = t156 * t165 - t174;
t55 = -Ifges(6,1) * t94 - Ifges(6,4) * t95;
t54 = -Ifges(7,1) * t94 + Ifges(7,5) * t95;
t53 = -Ifges(6,4) * t94 - Ifges(6,2) * t95;
t52 = -Ifges(7,5) * t94 + Ifges(7,3) * t95;
t51 = t95 * mrSges(6,1) - t94 * mrSges(6,2);
t50 = t95 * mrSges(7,1) + t94 * mrSges(7,3);
t46 = t67 + t204;
t44 = mrSges(6,1) * t85 + mrSges(6,2) * t86;
t43 = mrSges(7,1) * t85 - mrSges(7,3) * t86;
t34 = Ifges(6,1) * t86 - Ifges(6,4) * t85 + Ifges(6,5) * t142;
t33 = Ifges(7,1) * t86 + Ifges(7,4) * t142 + Ifges(7,5) * t85;
t32 = Ifges(6,4) * t86 - Ifges(6,2) * t85 + Ifges(6,6) * t142;
t31 = Ifges(7,5) * t86 + Ifges(7,6) * t142 + Ifges(7,3) * t85;
t30 = pkin(5) * t104 - qJ(6) * t105 + t96;
t18 = pkin(5) * t95 + qJ(6) * t94 - qJD(6) * t105 + t162;
t17 = pkin(5) * t85 - qJ(6) * t86 + t70;
t13 = -pkin(5) * t142 - t14;
t12 = qJ(6) * t142 + t184;
t11 = t41 * mrSges(6,1) + t42 * mrSges(6,2);
t10 = t41 * mrSges(7,1) - t42 * mrSges(7,3);
t9 = Ifges(6,1) * t42 - Ifges(6,4) * t41 - Ifges(6,5) * t165;
t8 = Ifges(7,1) * t42 - Ifges(7,4) * t165 + Ifges(7,5) * t41;
t7 = Ifges(6,4) * t42 - Ifges(6,2) * t41 - Ifges(6,6) * t165;
t6 = Ifges(7,5) * t42 - Ifges(7,6) * t165 + Ifges(7,3) * t41;
t5 = -t41 * pkin(5) + t42 * qJ(6) + t86 * qJD(6) + t47;
t2 = pkin(5) * t165 - t4;
t1 = -qJ(6) * t165 + qJD(6) * t142 + t3;
t15 = [0.2e1 * t14 * t23 + 0.2e1 * t13 * t24 + 0.2e1 * t17 * t10 + 0.2e1 * t12 * t21 + (t8 + t9) * t86 + (t6 - t7) * t85 + (t33 + t34) * t42 + (t31 - t32) * t41 + t150 * t142 + (t88 * t194 + (t79 + t80) * t138 + (t77 - t78) * t137) * t140 + (t1 * t12 + t13 * t2 - t17 * t5) * t196 + (t46 * t71 + t56 * t72 + t59 * t84) * t198 + (t81 * t66 + t82 * t67) * t199 + 0.2e1 * t184 * t22 + (t14 * t4 + t184 * t3 - t47 * t70) * t197 + ((mrSges(3,1) * t195 + (-(2 * Ifges(3,4)) - t145) * t140 - t211 * t86 + t210 * t85) * t140 + (mrSges(3,2) * t195 + 0.2e1 * (Ifges(3,4) + t145) * t142 + (-(2 * Ifges(6,3)) - (2 * Ifges(7,2)) + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(5,2)) - (2 * Ifges(4,3)) + pkin(7) ^ 2 * t199 + (mrSges(4,2) * t194 + t206 * t138) * t138 + (mrSges(4,1) * t194 + (Ifges(5,3) + Ifges(4,2)) * t137 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t138) * t137) * t140) * t142) * qJD(2) - 0.2e1 * t5 * t43 - 0.2e1 * t47 * t44 + 0.2e1 * t70 * t11 + 0.2e1 * t1 * t73 + 0.2e1 * t3 * t74 + 0.2e1 * t4 * t75 + 0.2e1 * t2 * t76 + 0.2e1 * t84 * t87 + 0.2e1 * t59 * t97 + 0.2e1 * t82 * t98 + 0.2e1 * t81 * t99 + 0.2e1 * t72 * t100 + 0.2e1 * t71 * t101 + 0.2e1 * t67 * t106 + 0.2e1 * t66 * t107 + 0.2e1 * t56 * t108 + 0.2e1 * t46 * t109; t30 * t10 + (t79 / 0.2e1 + t80 / 0.2e1 + t56 * mrSges(5,2) - t66 * mrSges(4,3) + (-t107 + t108) * qJD(3) + (-t99 + t100) * qJ(3) + m(4) * (-qJ(3) * t66 - qJD(3) * t81) + m(5) * (qJ(3) * t56 + qJD(3) * t72) + (-m(5) * t84 + m(6) * t70 + t44 - t97) * qJD(4)) * t137 + (t64 / 0.2e1 + t65 / 0.2e1) * t42 + m(7) * (t1 * t69 + t12 * t28 + t13 * t29 + t17 * t18 + t2 * t68 - t30 * t5) + (t62 / 0.2e1 - t63 / 0.2e1) * t41 + (t54 / 0.2e1 + t55 / 0.2e1) * t86 + (t52 / 0.2e1 - t53 / 0.2e1) * t85 + t182 * t29 + t183 * t28 + t185 * t68 + t186 * t69 + (t31 / 0.2e1 - t32 / 0.2e1) * t95 - (t33 / 0.2e1 + t34 / 0.2e1) * t94 + (t8 / 0.2e1 + t9 / 0.2e1) * t105 + (t6 / 0.2e1 - t7 / 0.2e1) * t104 + m(6) * (-t14 * t29 + t184 * t28 + t3 * t69 - t4 * t68 - t47 * t96) + (-t104 * t3 - t105 * t4 + t14 * t94 - t184 * t95) * mrSges(6,3) + (-t1 * t104 + t105 * t2 - t12 * t95 - t13 * t94) * mrSges(7,2) + t209 * t142 / 0.2e1 + (t192 + t87) * t110 + t18 * t43 + (-t77 / 0.2e1 + t78 / 0.2e1 + t46 * mrSges(5,2) + t67 * mrSges(4,3) + (t106 + t109) * qJD(3) + (t98 + t101) * qJ(3) + m(4) * (qJ(3) * t67 + qJD(3) * t82) + m(5) * (qJ(3) * t46 + qJD(3) * t71)) * t138 + t17 * t50 + ((pkin(7) * mrSges(3,2) - Ifges(3,6) + (Ifges(4,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t138 + (Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1) * t137 + (-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t105 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t104) * t140 + (Ifges(3,5) + t137 * (-Ifges(5,3) * t138 + t176) / 0.2e1 - t137 * (Ifges(4,2) * t138 + t178) / 0.2e1 + (-m(4) * pkin(2) - mrSges(4,1) * t138 + mrSges(4,2) * t137 - mrSges(3,1)) * pkin(7) + (t206 * t137 - t175 + t177) * t138 / 0.2e1) * t142) * qJD(2) - t5 * t60 - t47 * t61 + t70 * t51 - pkin(2) * t88 + t96 * t11 + t59 * t113; 0.2e1 * t18 * t60 + 0.2e1 * t30 * t50 + 0.2e1 * t96 * t51 + (t62 - t63) * t95 - (t64 + t65) * t94 + (t54 + t55) * t105 + (t52 - t53) * t104 + 0.2e1 * (-t113 + t61) * t162 + (t18 * t30 + t152) * t196 + t167 * t199 + (t162 * t96 + t152) * t197 + (-t110 * t162 + t167) * t198 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * t208 + 0.2e1 * t207 * (-t104 * t28 + t105 * t29 - t68 * t94 - t69 * t95); t117 + t188 * t42 + t189 * t41 + (m(4) * pkin(7) - t179) * t164 + t192 + m(6) * t47 + m(7) * t5 + t88; -m(7) * t18 + t189 * t95 - t188 * t94 + (-m(5) - m(6)) * t162; 0; -t185 * t141 + t186 * t139 + (t139 * t182 + t141 * t183) * qJD(5) + m(7) * (t1 * t139 - t141 * t2 + (t12 * t141 + t13 * t139) * qJD(5)) + m(6) * (t139 * t3 + t141 * t4 + (-t139 * t14 + t141 * t184) * qJD(5)) + m(5) * t56 + t100; m(5) * t137 * qJD(3) + (m(6) + m(7)) * (t139 * t28 - t141 * t29 + t69 * t159 + t160 * t68) + t207 * (-t139 * t95 + t141 * t94 + (-t104 * t141 + t105 * t139) * qJD(5)); 0; 0; t4 * mrSges(6,1) - t3 * mrSges(6,2) - pkin(5) * t24 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) + qJD(6) * t73 + qJ(6) * t21 + t1 * mrSges(7,3) - t2 * mrSges(7,1) + (-Ifges(7,2) - Ifges(6,3)) * t165 + t150; t69 * t205 + (pkin(5) * t94 - qJ(6) * t95 - qJD(6) * t104) * mrSges(7,2) + t201 * t29 + t200 * t28 + t209; 0; t200 * t159 + (qJD(5) * t201 + t205) * t139; 0.2e1 * t203 * qJD(6); m(7) * t2 + t24; m(7) * t29 - t94 * mrSges(7,2); 0; m(7) * t160; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
