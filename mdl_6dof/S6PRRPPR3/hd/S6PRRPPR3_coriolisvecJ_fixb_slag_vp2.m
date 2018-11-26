% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:09:29
% EndTime: 2018-11-23 15:09:32
% DurationCPUTime: 3.59s
% Computational Cost: add. (2428->422), mult. (5921->550), div. (0->0), fcn. (3371->8), ass. (0->195)
t131 = sin(qJ(6));
t134 = cos(qJ(6));
t155 = mrSges(7,1) * t131 + mrSges(7,2) * t134;
t132 = sin(qJ(3));
t135 = cos(qJ(3));
t161 = pkin(5) * t132 + pkin(9) * t135;
t188 = qJD(2) * t135;
t190 = qJD(2) * t132;
t129 = qJD(2) * pkin(2);
t136 = cos(qJ(2));
t127 = sin(pkin(6));
t192 = qJD(1) * t127;
t99 = t136 * t192;
t86 = -t99 - t129;
t57 = -pkin(3) * t188 - qJ(4) * t190 + t86;
t35 = pkin(4) * t188 + qJD(5) - t57;
t17 = qJD(2) * t161 + t35;
t137 = -pkin(3) - pkin(4);
t125 = -pkin(9) + t137;
t110 = qJ(5) * t190;
t128 = cos(pkin(6));
t197 = t128 * t135;
t100 = qJD(1) * t197;
t133 = sin(qJ(2));
t173 = t133 * t192;
t85 = qJD(2) * pkin(8) + t173;
t199 = t132 * t85 - t100;
t175 = qJD(4) + t199;
t168 = -t110 + t175;
t21 = qJD(3) * t125 + t168;
t5 = -t131 * t21 + t134 * t17;
t6 = t131 * t17 + t134 * t21;
t158 = t131 * t6 + t134 * t5;
t223 = -t134 / 0.2e1;
t224 = t131 / 0.2e1;
t107 = qJD(6) + t190;
t79 = qJD(3) * t131 + t134 * t188;
t222 = Ifges(7,4) * t79;
t186 = qJD(3) * t134;
t78 = t131 * t188 - t186;
t23 = Ifges(7,2) * t78 + Ifges(7,6) * t107 - t222;
t77 = Ifges(7,4) * t78;
t24 = -Ifges(7,1) * t79 + Ifges(7,5) * t107 + t77;
t126 = qJD(3) * qJ(4);
t191 = qJD(1) * t132;
t55 = t128 * t191 + t135 * t85;
t41 = -qJ(5) * t188 + t55;
t33 = -t126 - t41;
t29 = qJD(3) * pkin(5) - t33;
t250 = -t158 * mrSges(7,3) + t29 * t155 - t223 * t24 - t224 * t23;
t196 = t131 * t136;
t87 = -t135 * pkin(3) - t132 * qJ(4) - pkin(2);
t76 = t135 * pkin(4) - t87;
t56 = t161 + t76;
t211 = pkin(8) - qJ(5);
t95 = t211 * t132;
t20 = t131 * t56 + t134 * t95;
t142 = pkin(5) * t135 + t125 * t132;
t141 = t142 * qJD(3);
t119 = t132 * qJD(4);
t185 = qJD(3) * t135;
t193 = qJ(4) * t185 + t119;
t32 = t141 + t193;
t62 = -qJD(5) * t132 + t185 * t211;
t249 = -qJD(6) * t20 - t131 * t62 + t134 * t32 - (-t132 * t196 - t133 * t134) * t192;
t19 = -t131 * t95 + t134 * t56;
t195 = t134 * t136;
t248 = qJD(6) * t19 + t131 * t32 + t134 * t62 - (-t131 * t133 + t132 * t195) * t192;
t247 = -mrSges(5,1) - mrSges(4,1);
t246 = mrSges(5,2) + mrSges(4,3);
t115 = Ifges(5,5) * t190;
t202 = Ifges(7,6) * t131;
t203 = Ifges(7,5) * t134;
t150 = -t202 + t203;
t205 = Ifges(7,4) * t134;
t152 = -Ifges(7,2) * t131 + t205;
t206 = Ifges(7,4) * t131;
t154 = Ifges(7,1) * t134 - t206;
t225 = t107 / 0.2e1;
t228 = -t79 / 0.2e1;
t229 = t78 / 0.2e1;
t239 = qJD(3) / 0.2e1;
t241 = qJD(2) / 0.2e1;
t242 = -qJD(2) / 0.2e1;
t45 = t126 + t55;
t245 = t35 * mrSges(6,2) + t57 * mrSges(5,1) + t86 * mrSges(4,1) + Ifges(5,6) * t239 - Ifges(5,3) * t188 / 0.2e1 + t115 / 0.2e1 + (Ifges(4,4) * t132 + t135 * Ifges(4,2)) * t242 + (-t135 * Ifges(6,1) - Ifges(6,4) * t132) * t241 + t150 * t225 - t33 * mrSges(6,3) - t45 * mrSges(5,2) - t55 * mrSges(4,3) + t152 * t229 + t154 * t228 - (Ifges(4,6) + Ifges(6,5)) * qJD(3) / 0.2e1 + t250;
t181 = qJD(3) * qJD(6);
t183 = qJD(6) * t135;
t48 = -t134 * t181 + (t131 * t183 + t132 * t186) * qJD(2);
t244 = -t48 / 0.2e1;
t187 = qJD(3) * t132;
t49 = t131 * t181 + (-t131 * t187 + t134 * t183) * qJD(2);
t243 = -t49 / 0.2e1;
t157 = -t131 * t5 + t134 * t6;
t182 = qJD(2) * qJD(3);
t170 = t135 * t182;
t194 = qJ(4) * t170 + qJD(2) * t119;
t15 = (t141 - t173) * qJD(2) + t194;
t16 = t55 * qJD(3) + (-qJ(5) * t185 + (-qJD(5) + t99) * t132) * qJD(2);
t1 = qJD(6) * t5 + t131 * t15 + t134 * t16;
t2 = -qJD(6) * t6 - t131 * t16 + t134 * t15;
t160 = -t1 * t134 + t131 * t2;
t40 = -t110 + t199;
t238 = -t40 - qJD(4);
t237 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t48 + Ifges(7,6) * t49;
t174 = t137 * qJD(3);
t28 = t174 + t168;
t89 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t190;
t236 = m(6) * t28 + t89;
t208 = -qJD(3) * mrSges(6,1) + mrSges(7,1) * t78 + mrSges(7,2) * t79 + mrSges(6,3) * t188;
t234 = -m(6) * t33 + m(7) * t29 - t208;
t233 = 0.2e1 * pkin(8);
t232 = m(4) / 0.2e1;
t231 = m(5) / 0.2e1;
t230 = -t78 / 0.2e1;
t227 = t79 / 0.2e1;
t226 = -t107 / 0.2e1;
t169 = t132 * t182;
t171 = qJD(2) * t127 * t136;
t166 = t135 * t171;
t25 = qJD(1) * t166 + qJD(3) * t100 - t187 * t85;
t18 = qJD(3) * qJD(4) + t25;
t184 = qJD(5) * t135;
t13 = -qJ(5) * t169 + qJD(2) * t184 - t18;
t198 = t127 * t133;
t64 = t128 * t132 + t135 * t198;
t220 = t13 * t64;
t96 = t211 * t135;
t219 = t13 * t96;
t26 = t85 * t185 + (qJD(3) * t128 + t171) * t191;
t176 = t132 * t198;
t63 = t176 - t197;
t215 = t26 * t63;
t212 = mrSges(5,2) - mrSges(6,3);
t80 = (-mrSges(5,1) * t135 - mrSges(5,3) * t132) * qJD(2);
t83 = (mrSges(6,1) * t132 - mrSges(6,2) * t135) * qJD(2);
t210 = t83 - t80;
t209 = qJD(3) * t247 + t190 * t246;
t93 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t188;
t94 = mrSges(5,2) * t188 + qJD(3) * mrSges(5,3);
t207 = -t93 - t94;
t204 = Ifges(7,5) * t131;
t201 = Ifges(7,6) * t134;
t71 = mrSges(6,1) * t170 + mrSges(6,2) * t169;
t189 = qJD(2) * t133;
t179 = (-mrSges(4,1) * t135 + mrSges(4,2) * t132) * qJD(2) - t210;
t178 = t94 - t208;
t177 = pkin(3) * t187;
t172 = t127 * t189;
t167 = t132 * t174;
t165 = -0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * Ifges(6,4) + 0.3e1 / 0.2e1 * Ifges(5,5);
t164 = Ifges(4,5) / 0.2e1 + Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1;
t163 = -Ifges(4,6) / 0.2e1 - Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1;
t159 = t1 * t131 + t2 * t134;
t156 = -mrSges(7,1) * t134 + mrSges(7,2) * t131;
t153 = Ifges(7,1) * t131 + t205;
t151 = Ifges(7,2) * t134 + t206;
t30 = mrSges(7,1) * t170 - mrSges(7,3) * t48;
t31 = -mrSges(7,2) * t170 + mrSges(7,3) * t49;
t149 = -t131 * t30 + t134 * t31;
t52 = -mrSges(7,2) * t107 + mrSges(7,3) * t78;
t53 = mrSges(7,1) * t107 + mrSges(7,3) * t79;
t148 = -t131 * t53 + t134 * t52;
t147 = -t131 * t52 - t134 * t53;
t42 = t127 * t195 - t131 * t63;
t43 = t127 * t196 + t134 * t63;
t44 = -qJD(3) * pkin(3) + t175;
t144 = m(4) * t199 + m(5) * t44 + t209;
t143 = -m(4) * t55 - m(5) * t45 + t207;
t116 = Ifges(4,4) * t188;
t140 = t35 * mrSges(6,1) + t44 * mrSges(5,2) + t5 * mrSges(7,1) + t86 * mrSges(4,2) + t107 * Ifges(7,3) - t79 * Ifges(7,5) + t78 * Ifges(7,6) + (-Ifges(6,4) * t135 - t132 * Ifges(6,2)) * t242 + (t132 * Ifges(5,1) - Ifges(5,5) * t135) * t241 + Ifges(4,1) * t190 / 0.2e1 + t116 / 0.2e1 - t28 * mrSges(6,3) + t199 * mrSges(4,3) - t57 * mrSges(5,3) - t6 * mrSges(7,2) + (Ifges(6,6) + Ifges(5,4) + Ifges(4,5)) * t239;
t138 = qJD(2) ^ 2;
t130 = qJ(4) + pkin(5);
t112 = qJ(4) * t188;
t106 = Ifges(7,3) * t170;
t82 = pkin(3) * t190 - t112;
t73 = (mrSges(4,1) * t132 + mrSges(4,2) * t135) * t182;
t72 = (mrSges(5,1) * t132 - mrSges(5,3) * t135) * t182;
t61 = t177 - t193;
t60 = t187 * t211 + t184;
t59 = t137 * t190 + t112;
t58 = t167 + t193;
t39 = qJD(3) * t64 + t132 * t171;
t38 = qJD(3) * t176 - t128 * t185 - t166;
t37 = qJD(2) * t142 + t112;
t34 = (t173 + t177) * qJD(2) - t194;
t27 = (t167 - t173) * qJD(2) + t194;
t14 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t12 = t48 * Ifges(7,1) + t49 * Ifges(7,4) + Ifges(7,5) * t170;
t11 = t48 * Ifges(7,4) + t49 * Ifges(7,2) + Ifges(7,6) * t170;
t10 = t131 * t37 + t134 * t41;
t9 = -t131 * t41 + t134 * t37;
t8 = -qJD(6) * t43 - t131 * t39 - t134 * t172;
t7 = qJD(6) * t42 - t131 * t172 + t134 * t39;
t3 = [t64 * t14 + t42 * t30 + t43 * t31 + t7 * t52 + t8 * t53 + (t89 + t209) * t39 + (-t93 - t178) * t38 + ((-mrSges(3,2) * t138 + t71 - t72 - t73) * t136 + (-mrSges(3,1) * t138 + qJD(2) * t179) * t133) * t127 + m(4) * (t25 * t64 + t215 - t38 * t55 + t39 * t199 + (t86 - t99) * t172) + m(5) * (t18 * t64 + t215 - t38 * t45 + t39 * t44 + (-t136 * t34 + t189 * t57) * t127) + m(6) * (-t220 + t16 * t63 + t28 * t39 + t33 * t38 + (t136 * t27 - t189 * t35) * t127) + m(7) * (t1 * t43 + t2 * t42 - t29 * t38 + t5 * t8 + t6 * t7 - t220) + (-t132 * t64 + t135 * t63) * t182 * (mrSges(4,3) + t212); -pkin(2) * t73 + t96 * t14 + t19 * t30 + t20 * t31 + t58 * t83 + t61 * t80 + t62 * t89 + t76 * t71 + t87 * t72 + t208 * t60 + t249 * t53 + t248 * t52 + m(6) * (t16 * t95 + t27 * t76 + t28 * t62 + t33 * t60 + t35 * t58 - t219) + m(5) * (t34 * t87 + t57 * t61) + (t106 / 0.2e1 - t16 * mrSges(6,3) - t34 * mrSges(5,3) + t27 * mrSges(6,1) + ((t232 + t231) * t233 + t246) * t26 + (mrSges(4,2) * t189 + (-t144 - t236) * t136) * t192 + (t143 * pkin(8) + (t96 * mrSges(6,3) + t132 * t165) * qJD(2) + t163 * qJD(3) + t245) * qJD(3) + t237) * t132 + (-t34 * mrSges(5,1) - t27 * mrSges(6,2) + t11 * t224 + t18 * mrSges(5,2) + t25 * mrSges(4,3) + t12 * t223 + t154 * t244 + t152 * t243 + (mrSges(6,3) + t155) * t13 + t159 * mrSges(7,3) + (t18 * t231 + t232 * t25) * t233 + (-mrSges(4,1) * t189 + (t143 - t234) * t136) * t192 + ((t201 + t204) * t225 + t151 * t229 + t153 * t228 + t29 * t156 + t134 * t23 / 0.2e1 + t24 * t224 + t157 * mrSges(7,3)) * qJD(6) + (t144 * pkin(8) + ((-t203 / 0.2e1 + t202 / 0.2e1 - t165) * t135 - t95 * mrSges(6,3) + (-0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + Ifges(7,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t132) * qJD(2) + t164 * qJD(3) + t140) * qJD(3)) * t135 + (-t179 - m(5) * t57 + m(6) * t35 + 0.2e1 * (-t129 / 0.2e1 - t86 / 0.2e1) * m(4)) * t173 + (t1 * t20 + t19 * t2 + t248 * t6 + t249 * t5 - t29 * t60 - t219) * m(7); -t207 * t199 + t151 * t243 + t153 * t244 + t130 * t14 - t131 * t12 / 0.2e1 - t82 * t80 - t59 * t83 - t41 * t89 - t9 * t53 - t10 * t52 - t25 * mrSges(4,2) + t16 * mrSges(6,2) + t18 * mrSges(5,3) + t247 * t26 + (((Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t188 - t116 / 0.2e1 + (-t204 / 0.2e1 - t201 / 0.2e1 - t137 * mrSges(6,3) - pkin(3) * mrSges(5,2) + t164) * qJD(3) - t140) * t135 + ((-Ifges(5,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1) * t188 + (Ifges(4,4) / 0.2e1 + Ifges(6,4) / 0.2e1) * t190 - t115 / 0.2e1 + (-qJ(4) * t212 + t163) * qJD(3) - t245) * t132) * qJD(2) + (-pkin(3) * t26 + qJ(4) * t18 + t175 * t45 - t44 * t55 - t57 * t82) * m(5) + t160 * mrSges(7,3) + (-mrSges(6,1) + t156) * t13 + (-qJ(4) * t13 + t137 * t16 + t238 * t33 - t28 * t41 - t35 * t59) * m(6) + (-t10 * t6 - t13 * t130 - t238 * t29 - t5 * t9) * m(7) + (t149 - m(7) * t160 + (-m(7) * t158 + t147) * qJD(6)) * t125 + t178 * qJD(4) + (t150 * t226 + t152 * t230 + t154 * t227 - t250) * qJD(6) - t208 * t40 - t209 * t55 + t11 * t223; t147 * qJD(6) - t178 * qJD(3) + (t212 * t185 + (t147 - t210) * t132) * qJD(2) + t149 + (-qJD(3) * t29 - t107 * t158 - t160) * m(7) + (qJD(3) * t33 - t190 * t35 + t16) * m(6) + (-qJD(3) * t45 + t190 * t57 + t26) * m(5); t131 * t31 + t134 * t30 + t148 * qJD(6) + m(6) * t27 + m(7) * (qJD(6) * t157 + t159) + (t234 * t135 + (m(7) * t157 + t148 + t236) * t132) * qJD(2) + t71; t106 - t29 * (-mrSges(7,1) * t79 + mrSges(7,2) * t78) + (Ifges(7,1) * t78 + t222) * t227 + t23 * t228 + (Ifges(7,5) * t78 + Ifges(7,6) * t79) * t226 - t5 * t52 + t6 * t53 + (t5 * t78 - t6 * t79) * mrSges(7,3) + (Ifges(7,2) * t79 + t24 + t77) * t230 + t237;];
tauc  = t3(:);
