% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:11
% EndTime: 2019-03-08 21:08:20
% DurationCPUTime: 3.72s
% Computational Cost: add. (2428->421), mult. (5921->549), div. (0->0), fcn. (3371->8), ass. (0->195)
t130 = sin(qJ(6));
t133 = cos(qJ(6));
t155 = mrSges(7,1) * t130 + mrSges(7,2) * t133;
t131 = sin(qJ(3));
t134 = cos(qJ(3));
t161 = pkin(5) * t131 + pkin(9) * t134;
t187 = qJD(2) * t134;
t189 = qJD(2) * t131;
t128 = qJD(2) * pkin(2);
t135 = cos(qJ(2));
t126 = sin(pkin(6));
t191 = qJD(1) * t126;
t99 = t135 * t191;
t86 = -t99 - t128;
t57 = -pkin(3) * t187 - qJ(4) * t189 + t86;
t35 = pkin(4) * t187 + qJD(5) - t57;
t17 = qJD(2) * t161 + t35;
t136 = -pkin(3) - pkin(4);
t124 = -pkin(9) + t136;
t109 = qJ(5) * t189;
t127 = cos(pkin(6));
t196 = t127 * t134;
t100 = qJD(1) * t196;
t132 = sin(qJ(2));
t170 = t132 * t191;
t85 = qJD(2) * pkin(8) + t170;
t198 = t131 * t85 - t100;
t174 = qJD(4) + t198;
t167 = -t109 + t174;
t21 = qJD(3) * t124 + t167;
t5 = -t130 * t21 + t133 * t17;
t6 = t130 * t17 + t133 * t21;
t158 = t130 * t6 + t133 * t5;
t222 = -t133 / 0.2e1;
t223 = t130 / 0.2e1;
t106 = qJD(6) + t189;
t79 = qJD(3) * t130 + t133 * t187;
t221 = Ifges(7,4) * t79;
t185 = qJD(3) * t133;
t78 = t130 * t187 - t185;
t23 = Ifges(7,2) * t78 + Ifges(7,6) * t106 - t221;
t77 = Ifges(7,4) * t78;
t24 = -Ifges(7,1) * t79 + Ifges(7,5) * t106 + t77;
t125 = qJD(3) * qJ(4);
t190 = qJD(1) * t131;
t55 = t127 * t190 + t134 * t85;
t41 = -qJ(5) * t187 + t55;
t33 = -t125 - t41;
t29 = qJD(3) * pkin(5) - t33;
t249 = -t158 * mrSges(7,3) + t29 * t155 - t222 * t24 - t223 * t23;
t195 = t130 * t135;
t87 = -t134 * pkin(3) - t131 * qJ(4) - pkin(2);
t76 = t134 * pkin(4) - t87;
t56 = t161 + t76;
t210 = pkin(8) - qJ(5);
t95 = t210 * t131;
t20 = t130 * t56 + t133 * t95;
t141 = pkin(5) * t134 + t124 * t131;
t140 = t141 * qJD(3);
t118 = t131 * qJD(4);
t184 = qJD(3) * t134;
t192 = qJ(4) * t184 + t118;
t32 = t140 + t192;
t62 = -qJD(5) * t131 + t184 * t210;
t248 = -qJD(6) * t20 - t130 * t62 + t133 * t32 - (-t131 * t195 - t132 * t133) * t191;
t19 = -t130 * t95 + t133 * t56;
t194 = t133 * t135;
t247 = qJD(6) * t19 + t130 * t32 + t133 * t62 - (-t130 * t132 + t131 * t194) * t191;
t246 = -mrSges(5,1) - mrSges(4,1);
t245 = mrSges(5,2) + mrSges(4,3);
t114 = Ifges(5,5) * t189;
t201 = Ifges(7,6) * t130;
t202 = Ifges(7,5) * t133;
t150 = -t201 + t202;
t204 = Ifges(7,4) * t133;
t152 = -Ifges(7,2) * t130 + t204;
t205 = Ifges(7,4) * t130;
t154 = Ifges(7,1) * t133 - t205;
t224 = t106 / 0.2e1;
t227 = -t79 / 0.2e1;
t228 = t78 / 0.2e1;
t238 = qJD(3) / 0.2e1;
t240 = qJD(2) / 0.2e1;
t241 = -qJD(2) / 0.2e1;
t45 = t125 + t55;
t244 = t35 * mrSges(6,2) + t57 * mrSges(5,1) + t86 * mrSges(4,1) + Ifges(5,6) * t238 - Ifges(5,3) * t187 / 0.2e1 + t114 / 0.2e1 + (Ifges(4,4) * t131 + t134 * Ifges(4,2)) * t241 + (-t134 * Ifges(6,1) - Ifges(6,4) * t131) * t240 + t150 * t224 - t33 * mrSges(6,3) - t45 * mrSges(5,2) - t55 * mrSges(4,3) + t152 * t228 + t154 * t227 - (Ifges(4,6) + Ifges(6,5)) * qJD(3) / 0.2e1 + t249;
t180 = qJD(3) * qJD(6);
t182 = qJD(6) * t134;
t48 = -t133 * t180 + (t130 * t182 + t131 * t185) * qJD(2);
t243 = -t48 / 0.2e1;
t186 = qJD(3) * t131;
t49 = t130 * t180 + (-t130 * t186 + t133 * t182) * qJD(2);
t242 = -t49 / 0.2e1;
t40 = -t109 + t198;
t237 = qJD(4) + t40;
t181 = qJD(2) * qJD(3);
t169 = t134 * t181;
t193 = qJ(4) * t169 + qJD(2) * t118;
t15 = (t140 - t170) * qJD(2) + t193;
t16 = t55 * qJD(3) + (-qJ(5) * t184 + (-qJD(5) + t99) * t131) * qJD(2);
t1 = qJD(6) * t5 + t130 * t15 + t133 * t16;
t2 = -qJD(6) * t6 - t130 * t16 + t133 * t15;
t160 = -t1 * t133 + t2 * t130;
t157 = -t130 * t5 + t133 * t6;
t236 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t48 + Ifges(7,6) * t49;
t173 = t136 * qJD(3);
t28 = t173 + t167;
t89 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t189;
t235 = m(6) * t28 + t89;
t207 = -qJD(3) * mrSges(6,1) + mrSges(7,1) * t78 + mrSges(7,2) * t79 + mrSges(6,3) * t187;
t233 = -m(6) * t33 + m(7) * t29 - t207;
t232 = 0.2e1 * pkin(8);
t231 = m(4) / 0.2e1;
t230 = m(5) / 0.2e1;
t229 = -t78 / 0.2e1;
t226 = t79 / 0.2e1;
t225 = -t106 / 0.2e1;
t168 = t131 * t181;
t25 = qJD(3) * t100 - t186 * t85 + t99 * t187;
t18 = qJD(3) * qJD(4) + t25;
t183 = qJD(5) * t134;
t13 = -qJ(5) * t168 + qJD(2) * t183 - t18;
t197 = t126 * t132;
t64 = t127 * t131 + t134 * t197;
t219 = t13 * t64;
t96 = t210 * t134;
t218 = t13 * t96;
t171 = qJD(2) * t126 * t135;
t142 = qJD(3) * t127 + t171;
t26 = t142 * t190 + t184 * t85;
t175 = t131 * t197;
t63 = t175 - t196;
t214 = t26 * t63;
t211 = mrSges(6,3) - mrSges(5,2);
t80 = (-mrSges(5,1) * t134 - mrSges(5,3) * t131) * qJD(2);
t83 = (mrSges(6,1) * t131 - mrSges(6,2) * t134) * qJD(2);
t209 = t83 - t80;
t208 = qJD(3) * t246 + t189 * t245;
t93 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t187;
t94 = mrSges(5,2) * t187 + qJD(3) * mrSges(5,3);
t206 = t93 + t94;
t203 = Ifges(7,5) * t130;
t200 = Ifges(7,6) * t133;
t71 = mrSges(6,1) * t169 + mrSges(6,2) * t168;
t188 = qJD(2) * t132;
t178 = (-mrSges(4,1) * t134 + mrSges(4,2) * t131) * qJD(2) - t209;
t177 = t94 - t207;
t176 = pkin(3) * t186;
t172 = t126 * t188;
t166 = t131 * t173;
t165 = -0.3e1 / 0.2e1 * Ifges(6,4) + 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4);
t164 = -Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t163 = Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t159 = t1 * t130 + t2 * t133;
t156 = -mrSges(7,1) * t133 + mrSges(7,2) * t130;
t153 = Ifges(7,1) * t130 + t204;
t151 = Ifges(7,2) * t133 + t205;
t30 = mrSges(7,1) * t169 - mrSges(7,3) * t48;
t31 = -mrSges(7,2) * t169 + mrSges(7,3) * t49;
t149 = -t130 * t30 + t133 * t31;
t52 = -mrSges(7,2) * t106 + mrSges(7,3) * t78;
t53 = mrSges(7,1) * t106 + mrSges(7,3) * t79;
t148 = -t130 * t53 + t133 * t52;
t147 = -t130 * t52 - t133 * t53;
t42 = t126 * t194 - t130 * t63;
t43 = t126 * t195 + t133 * t63;
t44 = -qJD(3) * pkin(3) + t174;
t144 = m(4) * t198 + m(5) * t44 + t208;
t143 = -m(4) * t55 - m(5) * t45 - t206;
t115 = Ifges(4,4) * t187;
t139 = t35 * mrSges(6,1) + t44 * mrSges(5,2) + t5 * mrSges(7,1) + t86 * mrSges(4,2) + t106 * Ifges(7,3) - t79 * Ifges(7,5) + t78 * Ifges(7,6) + (-Ifges(6,4) * t134 - t131 * Ifges(6,2)) * t241 + (t131 * Ifges(5,1) - Ifges(5,5) * t134) * t240 + Ifges(4,1) * t189 / 0.2e1 + t115 / 0.2e1 - t28 * mrSges(6,3) + t198 * mrSges(4,3) - t57 * mrSges(5,3) - t6 * mrSges(7,2) + (Ifges(6,6) + Ifges(5,4) + Ifges(4,5)) * t238;
t137 = qJD(2) ^ 2;
t129 = qJ(4) + pkin(5);
t111 = qJ(4) * t187;
t105 = Ifges(7,3) * t169;
t82 = pkin(3) * t189 - t111;
t73 = (mrSges(4,1) * t131 + mrSges(4,2) * t134) * t181;
t72 = (mrSges(5,1) * t131 - mrSges(5,3) * t134) * t181;
t61 = t176 - t192;
t60 = t210 * t186 + t183;
t59 = t136 * t189 + t111;
t58 = t166 + t192;
t39 = -qJD(3) * t175 + t134 * t142;
t38 = qJD(3) * t64 + t131 * t171;
t37 = qJD(2) * t141 + t111;
t34 = (t170 + t176) * qJD(2) - t193;
t27 = (t166 - t170) * qJD(2) + t193;
t14 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t12 = t48 * Ifges(7,1) + t49 * Ifges(7,4) + Ifges(7,5) * t169;
t11 = t48 * Ifges(7,4) + t49 * Ifges(7,2) + Ifges(7,6) * t169;
t10 = t130 * t37 + t133 * t41;
t9 = -t130 * t41 + t133 * t37;
t8 = qJD(6) * t42 - t130 * t172 + t133 * t38;
t7 = -qJD(6) * t43 - t130 * t38 - t133 * t172;
t3 = [t64 * t14 + t42 * t30 + t43 * t31 + t8 * t52 + t7 * t53 + (t89 + t208) * t38 + (t93 + t177) * t39 + ((-mrSges(3,2) * t137 + t71 - t72 - t73) * t135 + (-mrSges(3,1) * t137 + qJD(2) * t178) * t132) * t126 + m(4) * (t25 * t64 + t214 + t38 * t198 + t39 * t55 + (t86 - t99) * t172) + m(5) * (t18 * t64 + t214 + t38 * t44 + t39 * t45 + (-t135 * t34 + t188 * t57) * t126) + m(6) * (-t219 + t16 * t63 + t28 * t38 - t33 * t39 + (t135 * t27 - t188 * t35) * t126) + m(7) * (t1 * t43 + t2 * t42 + t29 * t39 + t5 * t7 + t6 * t8 - t219) + (-t131 * t64 + t134 * t63) * t181 * (mrSges(4,3) - t211); -pkin(2) * t73 + t96 * t14 + t19 * t30 + t20 * t31 + t58 * t83 + t61 * t80 + t62 * t89 + t76 * t71 + t87 * t72 + t207 * t60 + t248 * t53 + t247 * t52 + m(6) * (t16 * t95 + t27 * t76 + t28 * t62 + t33 * t60 + t35 * t58 - t218) + m(5) * (t34 * t87 + t57 * t61) + (t105 / 0.2e1 + t27 * mrSges(6,1) - t34 * mrSges(5,3) - t16 * mrSges(6,3) + ((t231 + t230) * t232 + t245) * t26 + (mrSges(4,2) * t188 + (-t144 - t235) * t135) * t191 + ((t96 * mrSges(6,3) + t131 * t165) * qJD(2) + t143 * pkin(8) + t164 * qJD(3) + t244) * qJD(3) + t236) * t131 + (-t27 * mrSges(6,2) - t34 * mrSges(5,1) + t12 * t222 + t11 * t223 + t25 * mrSges(4,3) + t18 * mrSges(5,2) + t154 * t243 + t152 * t242 + (mrSges(6,3) + t155) * t13 + t159 * mrSges(7,3) + (t18 * t230 + t231 * t25) * t232 + (-mrSges(4,1) * t188 + (t143 - t233) * t135) * t191 + (t133 * t23 / 0.2e1 + t24 * t223 + (t200 + t203) * t224 + t151 * t228 + t153 * t227 + t29 * t156 + t157 * mrSges(7,3)) * qJD(6) + (((-t202 / 0.2e1 + t201 / 0.2e1 - t165) * t134 - t95 * mrSges(6,3) + (-0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1) + Ifges(7,3) / 0.2e1) * t131) * qJD(2) + t144 * pkin(8) + t163 * qJD(3) + t139) * qJD(3)) * t134 + (-t178 - m(5) * t57 + m(6) * t35 + 0.2e1 * (-t128 / 0.2e1 - t86 / 0.2e1) * m(4)) * t170 + (t1 * t20 + t19 * t2 + t247 * t6 + t248 * t5 - t29 * t60 - t218) * m(7); t246 * t26 + (-pkin(3) * t26 + qJ(4) * t18 + t174 * t45 - t44 * t55 - t57 * t82) * m(5) + (((-Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t187 - t115 / 0.2e1 - t139 + (-t203 / 0.2e1 - t200 / 0.2e1 - t136 * mrSges(6,3) - pkin(3) * mrSges(5,2) + t163) * qJD(3)) * t134 + ((-Ifges(5,1) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1) * t187 + (Ifges(6,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t189 - t114 / 0.2e1 + (qJ(4) * t211 + t164) * qJD(3) - t244) * t131) * qJD(2) + t206 * t198 + (t150 * t225 + t152 * t229 + t154 * t226 - t249) * qJD(6) + t160 * mrSges(7,3) + t177 * qJD(4) - t130 * t12 / 0.2e1 + t129 * t14 - t41 * t89 - t82 * t80 - t59 * t83 - t10 * t52 - t9 * t53 + t18 * mrSges(5,3) - t25 * mrSges(4,2) + t16 * mrSges(6,2) + (-mrSges(6,1) + t156) * t13 + (-t10 * t6 - t129 * t13 + t237 * t29 - t5 * t9) * m(7) + (-qJ(4) * t13 + t136 * t16 - t237 * t33 - t28 * t41 - t35 * t59) * m(6) - t207 * t40 - t208 * t55 + (t149 - m(7) * t160 + (-m(7) * t158 + t147) * qJD(6)) * t124 + t11 * t222 + t151 * t242 + t153 * t243; t147 * qJD(6) - t177 * qJD(3) + (-t211 * t184 + (t147 - t209) * t131) * qJD(2) + t149 + (-t29 * qJD(3) - t106 * t158 - t160) * m(7) + (qJD(3) * t33 - t189 * t35 + t16) * m(6) + (-qJD(3) * t45 + t189 * t57 + t26) * m(5); t130 * t31 + t133 * t30 + t148 * qJD(6) + m(6) * t27 + m(7) * (t157 * qJD(6) + t159) + (t233 * t134 + (m(7) * t157 + t148 + t235) * t131) * qJD(2) + t71; t105 - t29 * (-mrSges(7,1) * t79 + mrSges(7,2) * t78) + (Ifges(7,1) * t78 + t221) * t226 + t23 * t227 + (Ifges(7,5) * t78 + Ifges(7,6) * t79) * t225 - t5 * t52 + t6 * t53 + (t5 * t78 - t6 * t79) * mrSges(7,3) + (Ifges(7,2) * t79 + t24 + t77) * t229 + t236;];
tauc  = t3(:);
