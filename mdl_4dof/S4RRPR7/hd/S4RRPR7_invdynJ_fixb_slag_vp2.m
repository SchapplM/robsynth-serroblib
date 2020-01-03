% Calculate vector of inverse dynamics joint torques for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:57
% EndTime: 2019-12-31 17:06:10
% DurationCPUTime: 5.88s
% Computational Cost: add. (2286->385), mult. (5338->533), div. (0->0), fcn. (3533->10), ass. (0->186)
t235 = -m(5) - m(4);
t238 = mrSges(4,2) - mrSges(5,3);
t108 = qJ(2) + pkin(7);
t106 = sin(t108);
t107 = cos(t108);
t112 = sin(qJ(2));
t115 = cos(qJ(2));
t95 = -mrSges(3,1) * t115 + mrSges(3,2) * t112;
t237 = -mrSges(4,1) * t107 + t106 * t238 + t95;
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t109 = sin(pkin(7));
t173 = cos(pkin(7));
t159 = qJD(1) * qJD(2);
t146 = t112 * t159;
t158 = qJDD(1) * t115;
t90 = -t146 + t158;
t91 = qJDD(1) * t112 + t115 * t159;
t56 = t109 * t90 + t173 * t91;
t142 = t173 * t112;
t165 = qJD(1) * t115;
t75 = -qJD(1) * t142 - t109 * t165;
t60 = qJD(2) * t114 + t111 * t75;
t26 = qJD(4) * t60 + qJDD(2) * t111 + t114 * t56;
t213 = t26 / 0.2e1;
t61 = qJD(2) * t111 - t114 * t75;
t27 = -qJD(4) * t61 + qJDD(2) * t114 - t111 * t56;
t212 = t27 / 0.2e1;
t55 = -t109 * t91 + t173 * t90;
t52 = qJDD(4) - t55;
t211 = t52 / 0.2e1;
t236 = t90 / 0.2e1;
t122 = -t109 * t112 + t115 * t173;
t74 = t122 * qJD(1);
t67 = Ifges(4,4) * t74;
t234 = t60 * Ifges(5,6);
t227 = qJD(4) - t74;
t233 = t227 * Ifges(5,3);
t232 = t74 * Ifges(4,2);
t231 = qJD(2) / 0.2e1;
t230 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t60 - mrSges(5,2) * t61 + mrSges(4,3) * t75;
t174 = qJDD(2) / 0.2e1;
t229 = Ifges(4,5) * qJD(2);
t228 = Ifges(4,6) * qJD(2);
t160 = qJD(4) * t114;
t77 = t122 * qJD(2);
t86 = t109 * t115 + t142;
t127 = t111 * t77 + t86 * t160;
t104 = pkin(5) * t158;
t83 = -pkin(5) * t146 + t104;
t84 = t91 * pkin(5);
t225 = t112 * t84 + t115 * t83;
t198 = pkin(2) * t115;
t103 = pkin(1) + t198;
t92 = -qJD(1) * t103 + qJD(3);
t31 = -pkin(3) * t74 + pkin(6) * t75 + t92;
t110 = -qJ(3) - pkin(5);
t96 = t110 * t115;
t89 = qJD(1) * t96;
t78 = t173 * t89;
t147 = t110 * t112;
t88 = qJD(1) * t147;
t82 = qJD(2) * pkin(2) + t88;
t48 = t109 * t82 - t78;
t40 = qJD(2) * pkin(6) + t48;
t10 = -t111 * t40 + t114 * t31;
t172 = qJDD(1) * pkin(1);
t66 = -pkin(2) * t90 + qJDD(3) - t172;
t12 = -pkin(3) * t55 - pkin(6) * t56 + t66;
t163 = qJD(3) * t112;
t44 = qJDD(2) * pkin(2) - qJ(3) * t91 - qJD(1) * t163 - t84;
t164 = qJD(2) * t112;
t153 = pkin(5) * t164;
t162 = qJD(3) * t115;
t50 = qJ(3) * t90 + t104 + (-t153 + t162) * qJD(1);
t18 = t109 * t44 + t173 * t50;
t16 = qJDD(2) * pkin(6) + t18;
t1 = qJD(4) * t10 + t111 * t12 + t114 * t16;
t11 = t111 * t31 + t114 * t40;
t2 = -qJD(4) * t11 - t111 * t16 + t114 * t12;
t224 = t1 * t114 - t111 * t2;
t223 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t222 = 0.2e1 * t174;
t182 = t109 * t89;
t47 = t173 * t82 + t182;
t220 = t92 * mrSges(4,2) - t47 * mrSges(4,3);
t219 = m(3) * pkin(1) + mrSges(2,1) - t237;
t218 = t2 * mrSges(5,1) - t1 * mrSges(5,2);
t217 = t92 * mrSges(4,1) + t10 * mrSges(5,1) - t11 * mrSges(5,2) - t48 * mrSges(4,3);
t215 = Ifges(5,1) * t213 + Ifges(5,4) * t212 + Ifges(5,5) * t211;
t210 = -t60 / 0.2e1;
t209 = -t61 / 0.2e1;
t208 = t61 / 0.2e1;
t207 = -t227 / 0.2e1;
t205 = -t75 / 0.2e1;
t203 = t114 / 0.2e1;
t202 = Ifges(4,4) * t75;
t201 = Ifges(5,4) * t61;
t200 = pkin(2) * t109;
t197 = pkin(5) * t115;
t196 = pkin(6) * t106;
t187 = Ifges(3,4) * t112;
t186 = Ifges(3,4) * t115;
t185 = Ifges(5,4) * t111;
t184 = Ifges(5,4) * t114;
t181 = t111 * t74;
t179 = t111 * t86;
t177 = t114 * t74;
t176 = t114 * t86;
t116 = cos(qJ(1));
t171 = t111 * t116;
t113 = sin(qJ(1));
t170 = t113 * t111;
t169 = t113 * t114;
t168 = t114 * t116;
t166 = qJD(1) * t112;
t161 = qJD(4) * t111;
t157 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t166) * t197;
t156 = Ifges(5,5) * t26 + Ifges(5,6) * t27 + Ifges(5,3) * t52;
t155 = pkin(2) * t166;
t154 = pkin(2) * t164;
t59 = Ifges(5,4) * t60;
t21 = Ifges(5,1) * t61 + Ifges(5,5) * t227 + t59;
t150 = t21 * t203;
t149 = t173 * pkin(2);
t148 = -t55 * mrSges(4,1) + t56 * mrSges(4,2);
t144 = -t161 / 0.2e1;
t143 = qJD(2) * t110;
t139 = pkin(3) * t107 + t196;
t138 = -g(1) * t113 + g(2) * t116;
t137 = mrSges(3,1) * t112 + mrSges(3,2) * t115;
t135 = -mrSges(5,1) * t114 + mrSges(5,2) * t111;
t134 = mrSges(5,1) * t111 + mrSges(5,2) * t114;
t133 = Ifges(5,1) * t114 - t185;
t132 = t115 * Ifges(3,2) + t187;
t131 = -Ifges(5,2) * t111 + t184;
t130 = Ifges(3,5) * t115 - Ifges(3,6) * t112;
t129 = Ifges(5,5) * t114 - Ifges(5,6) * t111;
t46 = -pkin(3) * t122 - pkin(6) * t86 - t103;
t58 = t109 * t147 - t173 * t96;
t22 = -t111 * t58 + t114 * t46;
t23 = t111 * t46 + t114 * t58;
t126 = -t114 * t77 + t161 * t86;
t125 = pkin(1) * t137;
t39 = -qJD(2) * pkin(3) - t47;
t124 = t39 * t134;
t17 = -t109 * t50 + t173 * t44;
t123 = t112 * (Ifges(3,1) * t115 - t187);
t120 = m(5) * pkin(3) - t135;
t119 = t115 * t143 - t163;
t105 = Ifges(3,4) * t165;
t102 = -t149 - pkin(3);
t94 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t165;
t81 = Ifges(3,1) * t166 + Ifges(3,5) * qJD(2) + t105;
t80 = Ifges(3,6) * qJD(2) + qJD(1) * t132;
t76 = t86 * qJD(2);
t73 = t112 * t143 + t162;
t72 = t107 * t168 + t170;
t71 = -t107 * t171 + t169;
t70 = -t107 * t169 + t171;
t69 = t107 * t170 + t168;
t62 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t74;
t54 = t173 * t88 + t182;
t53 = t109 * t88 - t78;
t45 = -mrSges(4,1) * t74 - mrSges(4,2) * t75;
t42 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t56;
t41 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t55;
t38 = -t75 * Ifges(4,1) + t229 + t67;
t37 = -t202 + t228 + t232;
t35 = pkin(3) * t76 - pkin(6) * t77 + t154;
t34 = -pkin(3) * t75 - pkin(6) * t74 + t155;
t33 = t109 * t119 + t173 * t73;
t30 = mrSges(5,1) * t227 - mrSges(5,3) * t61;
t29 = -mrSges(5,2) * t227 + mrSges(5,3) * t60;
t20 = Ifges(5,2) * t60 + Ifges(5,6) * t227 + t201;
t19 = t61 * Ifges(5,5) + t233 + t234;
t15 = -qJDD(2) * pkin(3) - t17;
t14 = t111 * t34 + t114 * t54;
t13 = -t111 * t54 + t114 * t34;
t9 = -mrSges(5,2) * t52 + mrSges(5,3) * t27;
t8 = mrSges(5,1) * t52 - mrSges(5,3) * t26;
t7 = -mrSges(5,1) * t27 + mrSges(5,2) * t26;
t6 = -qJD(4) * t23 - t111 * t33 + t114 * t35;
t5 = qJD(4) * t22 + t111 * t35 + t114 * t33;
t3 = t26 * Ifges(5,4) + t27 * Ifges(5,2) + t52 * Ifges(5,6);
t4 = [(-t72 * mrSges(5,1) - t71 * mrSges(5,2) + t235 * (t116 * t103 - t113 * t110) + t223 * t113 + (-m(5) * t139 - t219) * t116) * g(2) + (-t70 * mrSges(5,1) - t69 * mrSges(5,2) + (-t110 * t235 + t223) * t116 + (-m(5) * (-t103 - t139) + m(4) * t103 + t219) * t113) * g(1) + (Ifges(3,1) * t91 + Ifges(3,4) * t236 - pkin(5) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t91) + t222 * Ifges(3,5)) * t112 + (Ifges(4,1) * t205 + t67 / 0.2e1 + t38 / 0.2e1 + t150 + t220) * t77 - t103 * t148 - t127 * t20 / 0.2e1 + t176 * t215 + t91 * t186 / 0.2e1 - t95 * t172 + (-m(4) * t47 + m(5) * t39 - t230) * (t109 * t73 - t119 * t173) + Ifges(3,6) * t115 * t174 + t39 * (mrSges(5,1) * t127 - mrSges(5,2) * t126) + (t197 * t90 + t225) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t225) + (t66 * mrSges(4,2) - t17 * mrSges(4,3) + Ifges(4,1) * t56 + Ifges(4,4) * t55 + Ifges(4,5) * t222 + t129 * t211 + t131 * t212 + t133 * t213 + t15 * t134 + t21 * t144) * t86 + t227 * (-Ifges(5,5) * t126 - Ifges(5,6) * t127) / 0.2e1 - (Ifges(5,3) * t211 + Ifges(5,6) * t212 + Ifges(5,5) * t213 - Ifges(4,4) * t56 - Ifges(4,2) * t55 + t66 * mrSges(4,1) + t156 / 0.2e1 - t18 * mrSges(4,3) - t222 * Ifges(4,6) + t218) * t122 + t115 * (Ifges(3,4) * t91 + Ifges(3,2) * t90 + Ifges(3,6) * qJDD(2)) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t90 + mrSges(3,2) * t91) + t33 * t62 + t58 * t41 - t3 * t179 / 0.2e1 + (-t1 * t179 + t10 * t126 - t11 * t127 - t176 * t2) * mrSges(5,3) - t94 * t153 + (t115 * (-Ifges(3,2) * t112 + t186) + t123) * t159 / 0.2e1 + (-m(4) * t17 + m(5) * t15 - t42 + t7) * (-t109 * t96 - t110 * t142) + m(4) * (-t103 * t66 + t154 * t92 + t18 * t58 + t33 * t48) + m(5) * (t1 * t23 + t10 * t6 + t11 * t5 + t2 * t22) + t45 * t154 - qJDD(2) * mrSges(3,2) * t197 + t5 * t29 + t6 * t30 + t22 * t8 + t23 * t9 - t125 * t159 + t60 * (-Ifges(5,4) * t126 - Ifges(5,2) * t127) / 0.2e1 + (-Ifges(5,1) * t126 - Ifges(5,4) * t127) * t208 + (t130 * t231 - t157) * qJD(2) + (-Ifges(4,4) * t205 + Ifges(5,5) * t208 + t233 / 0.2e1 + t234 / 0.2e1 + t19 / 0.2e1 - t37 / 0.2e1 - t232 / 0.2e1 + t217) * t76 - t80 * t164 / 0.2e1 + Ifges(2,3) * qJDD(1) + (Ifges(4,5) * t77 - Ifges(4,6) * t76 + t115 * t81) * t231 + t132 * t236; -t21 * t177 / 0.2e1 + (t80 / 0.2e1 + pkin(5) * t94) * t166 + t15 * t135 + (-Ifges(5,3) * t207 - Ifges(5,5) * t209 - Ifges(5,6) * t210 - t228 / 0.2e1 + t217) * t75 + (t129 * t207 + t133 * t209 + t131 * t210 - t124 - t229 / 0.2e1 - t220) * t74 + t41 * t200 + t3 * t203 + t37 * t205 + (Ifges(5,5) * t111 + Ifges(5,6) * t114) * t211 + (Ifges(5,2) * t114 + t185) * t212 + (Ifges(5,1) * t111 + t184) * t213 + t230 * t53 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + ((-t161 + t181) * t11 + (-t160 + t177) * t10 + t224) * mrSges(5,3) + (-t30 * t160 - t29 * t161 + t114 * t9 - t111 * t8 + m(5) * ((-t10 * t114 - t11 * t111) * qJD(4) + t224)) * (pkin(6) + t200) + (t129 * t227 + t131 * t60 + t133 * t61) * qJD(4) / 0.2e1 + (-m(5) * (t196 + t198) - t120 * t107 - m(4) * t198 + t237) * g(3) + Ifges(3,5) * t91 + t102 * t7 + Ifges(3,6) * t90 - t83 * mrSges(3,2) - t84 * mrSges(3,1) - t54 * t62 + Ifges(4,6) * t55 + Ifges(4,5) * t56 + t111 * t215 - t45 * t155 - (-Ifges(3,2) * t166 + t105 + t81) * t165 / 0.2e1 + (Ifges(4,1) * t74 + t19 + t202) * t75 / 0.2e1 - (Ifges(4,2) * t75 + t38 + t67) * t74 / 0.2e1 - t14 * t29 - t13 * t30 + t17 * mrSges(4,1) - t18 * mrSges(4,2) + (t150 + t124) * qJD(4) - t130 * t159 / 0.2e1 + (t157 + (-t123 / 0.2e1 + t125) * qJD(1)) * qJD(1) + (t181 / 0.2e1 + t144) * t20 + (g(1) * t116 + g(2) * t113) * (t137 - t235 * pkin(2) * t112 + (-m(5) * pkin(6) + t238) * t107 + (mrSges(4,1) + t120) * t106) + ((t109 * t18 + t17 * t173) * pkin(2) - t92 * t155 + t47 * t53 - t48 * t54) * m(4) + t42 * t149 + (-t10 * t13 + t102 * t15 - t11 * t14 - t39 * t53) * m(5); -t74 * t62 - t230 * t75 + (t227 * t29 + t8) * t114 + (-t227 * t30 + t9) * t111 + t148 + (t1 * t111 + t2 * t114 + t39 * t75 + t138 + t227 * (-t10 * t111 + t11 * t114)) * m(5) + (-t47 * t75 - t48 * t74 + t138 + t66) * m(4); -t39 * (mrSges(5,1) * t61 + mrSges(5,2) * t60) + (Ifges(5,1) * t60 - t201) * t209 + t20 * t208 + (Ifges(5,5) * t60 - Ifges(5,6) * t61) * t207 - t10 * t29 + t11 * t30 - g(1) * (mrSges(5,1) * t71 - mrSges(5,2) * t72) - g(2) * (-mrSges(5,1) * t69 + mrSges(5,2) * t70) + g(3) * t134 * t106 + (t10 * t60 + t11 * t61) * mrSges(5,3) + t156 + (-Ifges(5,2) * t61 + t21 + t59) * t210 + t218;];
tau = t4;
