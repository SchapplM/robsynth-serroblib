% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:18
% EndTime: 2019-03-08 18:49:24
% DurationCPUTime: 3.51s
% Computational Cost: add. (3284->396), mult. (8787->547), div. (0->0), fcn. (6835->12), ass. (0->203)
t125 = cos(qJ(6));
t122 = sin(qJ(6));
t203 = Ifges(7,4) * t122;
t148 = Ifges(7,2) * t125 + t203;
t245 = -t148 / 0.2e1;
t235 = qJD(3) / 0.2e1;
t153 = mrSges(7,1) * t125 - mrSges(7,2) * t122;
t128 = -pkin(4) - pkin(10);
t123 = sin(qJ(4));
t121 = cos(pkin(6));
t107 = qJD(1) * t121 + qJD(2);
t116 = sin(pkin(12));
t118 = sin(pkin(6));
t124 = sin(qJ(3));
t127 = cos(qJ(3));
t119 = cos(pkin(12));
t120 = cos(pkin(7));
t186 = t119 * t120;
t134 = (t116 * t127 + t124 * t186) * t118;
t117 = sin(pkin(7));
t188 = t117 * t124;
t47 = qJD(1) * t134 + t107 * t188;
t45 = qJD(3) * pkin(9) + t47;
t159 = pkin(5) * qJD(3) + t45;
t154 = t159 * t123;
t126 = cos(qJ(4));
t184 = qJD(1) * t118;
t164 = t119 * t184;
t74 = t107 * t120 - t117 * t164;
t68 = t126 * t74;
t22 = t68 - t154;
t240 = qJD(5) - t22;
t20 = qJD(4) * t128 + t240;
t189 = t116 * t124;
t46 = t127 * (t107 * t117 + t120 * t164) - t184 * t189;
t160 = -qJ(5) * t123 - pkin(3);
t87 = t126 * t128 + t160;
t36 = qJD(3) * t87 - t46;
t5 = -t122 * t36 + t125 * t20;
t6 = t122 * t20 + t125 * t36;
t155 = t122 * t5 - t125 * t6;
t181 = qJD(3) * t126;
t114 = pkin(5) * t181;
t193 = t123 * t74;
t31 = t126 * t45 + t193;
t27 = -qJD(4) * qJ(5) - t31;
t21 = t114 - t27;
t215 = -t125 / 0.2e1;
t216 = -t122 / 0.2e1;
t183 = qJD(3) * t123;
t110 = qJD(6) + t183;
t178 = qJD(4) * t125;
t90 = -t122 * t181 + t178;
t214 = Ifges(7,4) * t90;
t89 = -qJD(4) * t122 - t125 * t181;
t49 = t89 * Ifges(7,2) + t110 * Ifges(7,6) + t214;
t88 = Ifges(7,4) * t89;
t50 = t90 * Ifges(7,1) + t110 * Ifges(7,5) + t88;
t244 = t155 * mrSges(7,3) + t21 * t153 + t215 * t49 + t216 * t50;
t243 = -t181 / 0.2e1;
t226 = 2 * m(5);
t242 = -t226 / 0.2e1;
t241 = mrSges(6,2) - mrSges(5,1);
t232 = mrSges(5,3) + mrSges(6,1);
t135 = (t127 * t186 - t189) * t118;
t187 = t117 * t127;
t239 = t121 * t187 + t135;
t198 = Ifges(7,6) * t125;
t201 = Ifges(7,5) * t122;
t147 = t198 + t201;
t202 = Ifges(7,4) * t125;
t150 = Ifges(7,1) * t122 + t202;
t217 = t110 / 0.2e1;
t222 = t90 / 0.2e1;
t233 = qJD(4) / 0.2e1;
t234 = -qJD(4) / 0.2e1;
t97 = -pkin(4) * t126 + t160;
t39 = qJD(3) * t97 - t46;
t44 = -qJD(3) * pkin(3) - t46;
t238 = -t44 * mrSges(5,1) - t27 * mrSges(6,1) + t39 * mrSges(6,2) + t31 * mrSges(5,3) - Ifges(6,5) * t233 - Ifges(5,6) * t234 - t147 * t217 - t150 * t222 + t89 * t245 + t244 + ((Ifges(5,2) + Ifges(6,3)) * t126 + (Ifges(5,4) + Ifges(6,6)) * t123) * t235;
t237 = -Ifges(5,1) / 0.2e1;
t236 = Ifges(5,4) * t243;
t195 = t123 * t45;
t30 = -t68 + t195;
t231 = -qJD(5) - t30;
t113 = pkin(4) * t183;
t109 = qJD(4) * t113;
t176 = qJD(5) * t123;
t132 = t47 - t176;
t146 = pkin(10) * t123 - qJ(5) * t126;
t137 = t146 * qJD(4);
t32 = t109 + (t137 + t132) * qJD(3);
t42 = (qJD(1) * t135 + t107 * t187) * qJD(3);
t196 = t123 * t42;
t8 = t196 + (t126 * t159 + t193) * qJD(4);
t1 = qJD(6) * t5 + t122 * t8 + t125 * t32;
t2 = -qJD(6) * t6 - t122 * t32 + t125 * t8;
t230 = t1 * t122 + t125 * t2;
t229 = -m(5) * t31 + m(6) * t27;
t173 = qJD(4) * qJD(6);
t175 = qJD(6) * t126;
t179 = qJD(4) * t123;
t69 = -t122 * t173 + (t122 * t179 - t125 * t175) * qJD(3);
t70 = -t125 * t173 + (t122 * t175 + t123 * t178) * qJD(3);
t228 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t69 + Ifges(7,6) * t70;
t225 = -t46 / 0.2e1;
t224 = -t89 / 0.2e1;
t223 = -t90 / 0.2e1;
t221 = pkin(5) + pkin(9);
t218 = -t110 / 0.2e1;
t11 = qJD(4) * t31 + t196;
t57 = t121 * t188 + t134;
t78 = -t117 * t118 * t119 + t120 * t121;
t37 = t123 * t57 - t78 * t126;
t212 = t11 * t37;
t166 = t123 * t188;
t79 = -t126 * t120 + t166;
t211 = t11 * t79;
t43 = t47 * qJD(3);
t209 = t43 * t239;
t177 = qJD(4) * t126;
t206 = t126 * t42 + t74 * t177;
t174 = qJD(3) * qJD(4);
t85 = (-mrSges(6,2) * t123 - mrSges(6,3) * t126) * t174;
t86 = (mrSges(5,1) * t123 + mrSges(5,2) * t126) * t174;
t205 = t85 + t86;
t91 = (mrSges(6,2) * t126 - mrSges(6,3) * t123) * qJD(3);
t204 = t91 + (-mrSges(5,1) * t126 + mrSges(5,2) * t123) * qJD(3);
t200 = Ifges(7,5) * t125;
t199 = Ifges(7,6) * t122;
t194 = t123 * t46;
t192 = t127 * t43;
t101 = -mrSges(6,1) * t181 - qJD(4) * mrSges(6,3);
t60 = -mrSges(7,1) * t89 + mrSges(7,2) * t90;
t191 = t101 - t60;
t190 = -qJD(4) * t241 - t183 * t232;
t185 = qJD(4) * mrSges(5,2) - mrSges(5,3) * t181 + t101;
t182 = qJD(3) * t124;
t180 = qJD(3) * t127;
t172 = pkin(9) * t11 / 0.2e1;
t171 = -Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t170 = Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1;
t169 = -0.3e1 / 0.2e1 * Ifges(6,6) - 0.3e1 / 0.2e1 * Ifges(5,4);
t168 = t60 - t185;
t104 = t221 * t126;
t165 = qJ(5) * t177;
t163 = t117 * t182;
t162 = t117 * t180;
t161 = t126 * t174;
t158 = pkin(4) * t179 - t176;
t152 = mrSges(7,1) * t122 + mrSges(7,2) * t125;
t151 = Ifges(7,1) * t125 - t203;
t149 = -Ifges(7,2) * t122 + t202;
t16 = t122 * t239 + t125 * t37;
t17 = t122 * t37 - t125 * t239;
t54 = mrSges(7,1) * t161 - mrSges(7,3) * t69;
t55 = -mrSges(7,2) * t161 + mrSges(7,3) * t70;
t145 = t122 * t55 + t125 * t54;
t71 = -mrSges(7,2) * t110 + mrSges(7,3) * t89;
t72 = mrSges(7,1) * t110 - mrSges(7,3) * t90;
t144 = -t122 * t72 + t125 * t71;
t38 = t123 * t78 + t126 * t57;
t103 = t221 * t123;
t58 = t103 * t125 - t122 * t87;
t59 = t103 * t122 + t125 * t87;
t139 = -t122 * t79 + t125 * t187;
t63 = t122 * t187 + t125 * t79;
t80 = t120 * t123 + t126 * t188;
t24 = -qJD(4) * pkin(4) - t231;
t131 = t39 * mrSges(6,3) + t6 * mrSges(7,2) - t110 * Ifges(7,3) - t90 * Ifges(7,5) - t89 * Ifges(7,6) + t183 * t237 + Ifges(5,5) * t234 + t236 + Ifges(6,4) * t233 + (-t123 * Ifges(6,2) - Ifges(6,6) * t126) * t235 - t24 * mrSges(6,1) - t30 * mrSges(5,3) - t44 * mrSges(5,2) - t5 * mrSges(7,1);
t129 = qJD(3) ^ 2;
t108 = Ifges(7,3) * t161;
t95 = qJD(4) * t104;
t94 = t221 * t179;
t93 = -qJ(5) * t181 + t113;
t77 = t158 - t165;
t76 = qJD(3) * t146 + t113;
t73 = t137 + t158;
t62 = qJD(4) * t80 + t123 * t162;
t61 = qJD(4) * t166 - t120 * t177 - t126 * t162;
t53 = t57 * qJD(3);
t52 = t239 * qJD(3);
t40 = -mrSges(7,1) * t70 + mrSges(7,2) * t69;
t35 = t69 * Ifges(7,1) + t70 * Ifges(7,4) + Ifges(7,5) * t161;
t34 = t69 * Ifges(7,4) + t70 * Ifges(7,2) + Ifges(7,6) * t161;
t33 = t109 + (t132 - t165) * qJD(3);
t29 = qJD(6) * t139 - t122 * t163 + t125 * t62;
t28 = qJD(6) * t63 + t122 * t62 + t125 * t163;
t26 = -qJD(6) * t59 - t122 * t73 + t125 * t95;
t25 = qJD(6) * t58 + t122 * t95 + t125 * t73;
t23 = t114 + t31;
t19 = t122 * t194 + t125 * t47;
t18 = -t122 * t47 + t125 * t194;
t15 = -t57 * t179 + (qJD(4) * t78 + t52) * t126;
t14 = qJD(4) * t38 + t123 * t52;
t13 = t122 * t23 + t125 * t76;
t12 = -t122 * t76 + t125 * t23;
t10 = -t179 * t45 + t206;
t9 = (-qJD(5) + t195) * qJD(4) - t206;
t7 = (qJD(5) - t154) * qJD(4) + t206;
t4 = qJD(6) * t16 + t122 * t14 + t125 * t53;
t3 = -qJD(6) * t17 - t122 * t53 + t125 * t14;
t41 = [t16 * t54 + t17 * t55 + t3 * t72 + t38 * t40 + t4 * t71 - t205 * t239 + t204 * t53 - t190 * t14 + t168 * t15 + m(7) * (t1 * t17 + t15 * t21 + t16 * t2 + t3 * t5 + t38 * t7 + t4 * t6) + m(4) * (t42 * t57 - t46 * t53 + t47 * t52 - t209) + m(5) * (t10 * t38 + t14 * t30 + t15 * t31 + t44 * t53 - t209 + t212) + m(6) * (t14 * t24 - t15 * t27 - t239 * t33 - t38 * t9 + t39 * t53 + t212) + (-t53 * mrSges(4,1) - t52 * mrSges(4,2) + t232 * qJD(4) * (-t123 * t38 + t126 * t37)) * qJD(3); t28 * t71 + t29 * t72 + t80 * t40 + t63 * t54 - t139 * t55 - t190 * t62 - t168 * t61 + m(5) * (t10 * t80 + t30 * t62 - t31 * t61 + t211) + m(7) * (-t1 * t139 + t2 * t63 - t21 * t61 + t28 * t6 + t29 * t5 + t7 * t80) + m(6) * (t24 * t62 + t27 * t61 - t80 * t9 + t211) + ((-mrSges(4,2) * t129 - t205) * t127 + (-t129 * mrSges(4,1) + qJD(3) * t204) * t124 + m(4) * (t124 * t42 + t180 * t47 - t182 * t46 - t192) + m(5) * (t182 * t44 - t192) + m(6) * (-t127 * t33 + t182 * t39)) * t117 + t232 * t174 * (-t123 * t80 + t126 * t79); -t43 * mrSges(4,1) + t104 * t40 + t58 * t54 + t59 * t55 - t94 * t60 + t77 * t91 + t97 * t85 + (t26 - t18) * t72 + (t25 - t19) * t71 + (qJD(3) * t46 - t42) * mrSges(4,2) + m(7) * (t1 * t59 + t104 * t7 + t2 * t58 - t21 * t94 + t25 * t6 + t26 * t5) - m(7) * (t18 * t5 + t19 * t6) + m(6) * (t33 * t97 + t39 * t77) + (t43 * mrSges(5,2) - t33 * mrSges(6,3) + t108 / 0.2e1 + t190 * t46 + t232 * t11 + 0.2e1 * (t225 * t24 + t172) * m(6) + (t225 * t30 + t172) * t226 + (t169 * t183 + t170 * qJD(4) + (t185 + t229) * pkin(9) - t238) * qJD(4) + t228) * t123 + (t35 * t216 - t9 * mrSges(6,1) + t10 * mrSges(5,3) + t34 * t215 + t7 * t153 - t69 * t150 / 0.2e1 + t70 * t245 - t43 * mrSges(5,1) + t33 * mrSges(6,2) + (-t1 * t125 + t2 * t122) * mrSges(7,3) + ((t199 - t200) * t217 + t149 * t224 + t151 * t223 - t21 * t152 + t50 * t215 + t122 * t49 / 0.2e1 + (t6 * t122 + t5 * t125) * mrSges(7,3)) * qJD(6) + (-t131 + ((-t201 / 0.2e1 - t198 / 0.2e1 - t169) * t126 + (0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(6,3) + 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) + Ifges(7,3) / 0.2e1) * t123) * qJD(3) + t171 * qJD(4)) * qJD(4) + (m(5) * t10 - m(6) * t9 + (m(5) * t30 + m(6) * t24 - t190) * qJD(4)) * pkin(9) + (-m(7) * t21 - t168 + t229) * t46) * t126 + (t242 * t43 - t86) * pkin(3) + (-m(6) * t39 + qJD(3) * mrSges(4,1) + t242 * t44 - t204) * t47; t7 * t152 + t69 * t151 / 0.2e1 + t70 * t149 / 0.2e1 + t125 * t35 / 0.2e1 + t34 * t216 - t93 * t91 - t13 * t71 - t12 * t72 - t22 * t60 + qJ(5) * t40 - t9 * mrSges(6,3) - t10 * mrSges(5,2) + t190 * t31 - t185 * t30 + t241 * t11 - t191 * qJD(5) - t230 * mrSges(7,3) + (t147 * t218 + t148 * t224 + t150 * t223 + t244) * qJD(6) + ((t131 + Ifges(6,6) * t243 + (t200 / 0.2e1 - t199 / 0.2e1 - pkin(4) * mrSges(6,1) + t171) * qJD(4) + t236) * t126 + ((Ifges(5,2) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + t237) * t181 + (-qJ(5) * mrSges(6,1) + t170) * qJD(4) + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t183 + t238) * t123) * qJD(3) + (qJ(5) * t7 - t12 * t5 - t13 * t6 + t240 * t21) * m(7) + (t145 + m(7) * t230 + (-m(7) * t155 + t144) * qJD(6)) * t128 + (-pkin(4) * t11 - qJ(5) * t9 + t231 * t27 - t24 * t31 - t39 * t93) * m(6); t144 * qJD(6) + t191 * qJD(4) + (mrSges(6,1) * t177 + (t144 + t91) * t123) * qJD(3) + t145 + (-qJD(4) * t21 - t110 * t155 + t230) * m(7) + (qJD(4) * t27 + t183 * t39 + t11) * m(6); t108 - t21 * (mrSges(7,1) * t90 + mrSges(7,2) * t89) + (Ifges(7,1) * t89 - t214) * t223 + t49 * t222 + (Ifges(7,5) * t89 - Ifges(7,6) * t90) * t218 - t5 * t71 + t6 * t72 + (t5 * t89 + t6 * t90) * mrSges(7,3) + (-Ifges(7,2) * t90 + t50 + t88) * t224 + t228;];
tauc  = t41(:);
