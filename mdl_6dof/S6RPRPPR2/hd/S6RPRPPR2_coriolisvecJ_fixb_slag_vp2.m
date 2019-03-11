% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:16
% EndTime: 2019-03-09 02:41:27
% DurationCPUTime: 5.48s
% Computational Cost: add. (3850->440), mult. (9388->573), div. (0->0), fcn. (6087->8), ass. (0->214)
t248 = mrSges(6,2) - mrSges(5,1);
t212 = mrSges(6,1) + mrSges(5,3);
t247 = Ifges(5,5) - Ifges(6,4);
t246 = -Ifges(5,6) + Ifges(6,5);
t129 = sin(pkin(10));
t133 = sin(qJ(3));
t135 = cos(qJ(3));
t191 = cos(pkin(10));
t109 = t129 * t135 + t133 * t191;
t102 = t109 * qJD(1);
t164 = t191 * t135;
t159 = qJD(1) * t164;
t184 = qJD(1) * t133;
t100 = t129 * t184 - t159;
t132 = sin(qJ(6));
t134 = cos(qJ(6));
t76 = -qJD(3) * t132 + t100 * t134;
t77 = qJD(3) * t134 + t100 * t132;
t97 = Ifges(5,4) * t100;
t98 = qJD(6) + t102;
t245 = Ifges(5,1) * t102 + Ifges(5,5) * qJD(3) + Ifges(7,5) * t77 + Ifges(7,6) * t76 + Ifges(7,3) * t98 - t97;
t197 = t100 * mrSges(5,3);
t198 = t100 * mrSges(6,1);
t84 = -qJD(3) * mrSges(6,3) + t198;
t210 = -qJD(3) * mrSges(5,2) - t197 - t84;
t195 = t102 * mrSges(5,3);
t196 = t102 * mrSges(6,1);
t244 = -t248 * qJD(3) - t195 - t196;
t43 = -mrSges(7,1) * t76 + mrSges(7,2) * t77;
t243 = -t84 + t43;
t176 = qJD(1) * qJD(3);
t166 = t133 * t176;
t101 = t109 * qJD(3);
t91 = qJD(1) * t101;
t51 = qJD(6) * t76 + t132 * t91;
t92 = qJD(3) * t159 - t129 * t166;
t28 = mrSges(7,1) * t92 - mrSges(7,3) * t51;
t52 = -qJD(6) * t77 + t134 * t91;
t29 = -mrSges(7,2) * t92 + mrSges(7,3) * t52;
t242 = t132 * t29 + t134 * t28;
t177 = t133 * qJD(4);
t122 = sin(pkin(9)) * pkin(1) + pkin(7);
t111 = t122 * qJD(1);
t161 = qJ(4) * qJD(1) + t111;
t178 = t133 * qJD(2);
t79 = t135 * t161 + t178;
t136 = -qJD(1) * t177 - qJD(3) * t79;
t181 = qJD(4) * t135;
t182 = qJD(3) * t133;
t128 = t135 * qJD(2);
t86 = qJD(3) * t128 - t111 * t182;
t68 = (-qJ(4) * t182 + t181) * qJD(1) + t86;
t20 = t129 * t68 - t191 * t136;
t15 = pkin(5) * t92 + t20;
t119 = pkin(3) * t166;
t143 = -qJ(5) * t92 - qJD(5) * t102 + t119;
t229 = pkin(4) + pkin(8);
t16 = t229 * t91 + t143;
t70 = t129 * t79;
t78 = -t133 * t161 + t128;
t73 = qJD(3) * pkin(3) + t78;
t37 = t191 * t73 - t70;
t144 = qJD(5) - t37;
t218 = t102 * pkin(5);
t17 = -qJD(3) * t229 + t144 + t218;
t170 = -cos(pkin(9)) * pkin(1) - pkin(2);
t110 = -pkin(3) * t135 + t170;
t185 = qJD(1) * t110;
t99 = qJD(4) + t185;
t139 = -t102 * qJ(5) + t99;
t32 = t100 * t229 + t139;
t5 = -t132 * t32 + t134 * t17;
t1 = qJD(6) * t5 + t132 * t15 + t134 * t16;
t6 = t132 * t17 + t134 * t32;
t2 = -qJD(6) * t6 - t132 * t16 + t134 * t15;
t241 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t51 + Ifges(7,6) * t52;
t154 = mrSges(7,1) * t134 - mrSges(7,2) * t132;
t219 = pkin(5) * t100;
t167 = t191 * t79;
t38 = t129 * t73 + t167;
t35 = -qJD(3) * qJ(5) - t38;
t19 = -t35 - t219;
t223 = -t132 / 0.2e1;
t221 = Ifges(7,4) * t77;
t26 = Ifges(7,2) * t76 + Ifges(7,6) * t98 + t221;
t74 = Ifges(7,4) * t76;
t27 = Ifges(7,1) * t77 + Ifges(7,5) * t98 + t74;
t240 = t154 * t19 - t134 * t26 / 0.2e1 + t27 * t223;
t46 = t100 * pkin(4) + t139;
t239 = t99 * mrSges(5,1) - t46 * mrSges(6,2);
t238 = t5 * mrSges(7,1) + t99 * mrSges(5,2) - t6 * mrSges(7,2) - t46 * mrSges(6,3);
t237 = t51 / 0.2e1;
t236 = t52 / 0.2e1;
t235 = -t76 / 0.2e1;
t234 = t76 / 0.2e1;
t233 = -t77 / 0.2e1;
t232 = t77 / 0.2e1;
t231 = -t98 / 0.2e1;
t230 = t98 / 0.2e1;
t228 = -t100 / 0.2e1;
t227 = t100 / 0.2e1;
t226 = -t102 / 0.2e1;
t225 = t102 / 0.2e1;
t222 = t134 / 0.2e1;
t220 = pkin(3) * t129;
t186 = qJ(4) + t122;
t106 = t186 * t133;
t107 = t186 * t135;
t64 = t191 * t106 + t107 * t129;
t217 = t20 * t64;
t216 = t92 * mrSges(6,1);
t215 = -qJD(3) / 0.2e1;
t214 = qJD(3) / 0.2e1;
t211 = -Ifges(5,4) - Ifges(6,6);
t21 = t129 * t136 + t191 * t68;
t208 = mrSges(7,3) * t132;
t207 = mrSges(7,3) * t134;
t206 = Ifges(4,4) * t133;
t205 = Ifges(5,4) * t102;
t204 = Ifges(7,4) * t132;
t203 = Ifges(7,4) * t134;
t201 = Ifges(7,5) * t132;
t200 = Ifges(6,6) * t102;
t199 = Ifges(7,6) * t134;
t108 = t129 * t133 - t164;
t194 = t108 * t20;
t190 = Ifges(4,5) * qJD(3);
t189 = Ifges(4,6) * qJD(3);
t188 = t101 * t132;
t187 = t134 * t101;
t183 = qJD(1) * t135;
t180 = qJD(6) * t132;
t179 = qJD(6) * t134;
t175 = t43 + t210;
t126 = pkin(3) * t184;
t174 = mrSges(4,3) * t184;
t173 = mrSges(4,3) * t183;
t125 = Ifges(4,4) * t183;
t169 = t191 * pkin(3);
t168 = m(4) * t122 + mrSges(4,3);
t40 = t129 * t78 + t167;
t162 = qJD(3) * t186;
t80 = -t133 * t162 + t181;
t81 = -t135 * t162 - t177;
t44 = t129 * t80 - t191 * t81;
t163 = qJ(5) * t100 + t126;
t123 = -t169 - pkin(4);
t158 = t1 * t134 - t132 * t2;
t157 = t132 * t6 + t134 * t5;
t156 = t5 * t132 - t6 * t134;
t18 = -qJD(3) * qJD(5) - t21;
t113 = t170 * qJD(1);
t153 = mrSges(7,1) * t132 + mrSges(7,2) * t134;
t152 = Ifges(7,1) * t134 - t204;
t151 = Ifges(7,1) * t132 + t203;
t150 = -Ifges(7,2) * t132 + t203;
t149 = Ifges(7,2) * t134 + t204;
t148 = Ifges(7,5) * t134 - Ifges(7,6) * t132;
t147 = t199 + t201;
t140 = -t109 * qJ(5) + t110;
t42 = t108 * t229 + t140;
t47 = pkin(5) * t109 + t64;
t12 = t132 * t47 + t134 * t42;
t11 = -t132 * t42 + t134 * t47;
t53 = -mrSges(7,2) * t98 + mrSges(7,3) * t76;
t54 = mrSges(7,1) * t98 - mrSges(7,3) * t77;
t146 = -t132 * t54 + t134 * t53;
t145 = -t132 * t53 - t134 * t54;
t94 = t111 * t135 + t178;
t41 = t191 * t78 - t70;
t45 = t129 * t81 + t191 * t80;
t103 = qJD(3) * t164 - t129 * t182;
t142 = pkin(3) * t182 - qJ(5) * t103 - qJD(5) * t109;
t65 = -t129 * t106 + t107 * t191;
t141 = -t145 - t244;
t138 = -qJD(6) * t156 + t1 * t132 + t134 * t2;
t137 = qJD(6) * t146 + t242;
t120 = qJ(5) + t220;
t114 = -qJD(3) * mrSges(4,2) + t173;
t112 = qJD(3) * mrSges(4,1) - t174;
t105 = Ifges(4,1) * t184 + t125 + t190;
t104 = t189 + (Ifges(4,2) * t135 + t206) * qJD(1);
t96 = Ifges(6,6) * t100;
t93 = -t111 * t133 + t128;
t90 = Ifges(7,3) * t92;
t89 = t92 * mrSges(6,3);
t88 = t92 * mrSges(5,2);
t87 = t94 * qJD(3);
t67 = -mrSges(6,2) * t100 - mrSges(6,3) * t102;
t66 = mrSges(5,1) * t100 + mrSges(5,2) * t102;
t59 = -Ifges(5,2) * t100 + Ifges(5,6) * qJD(3) + t205;
t58 = Ifges(6,4) * qJD(3) - Ifges(6,2) * t102 + t96;
t57 = Ifges(6,5) * qJD(3) + Ifges(6,3) * t100 - t200;
t56 = t108 * pkin(4) + t140;
t55 = pkin(4) * t102 + t163;
t48 = -t108 * pkin(5) + t65;
t39 = pkin(4) * t101 + t142;
t36 = t102 * t229 + t163;
t34 = -qJD(3) * pkin(4) + t144;
t33 = pkin(4) * t91 + t143;
t31 = -t101 * pkin(5) + t45;
t30 = pkin(5) * t103 + t44;
t24 = t41 - t218;
t23 = t40 - t219;
t22 = t101 * t229 + t142;
t14 = -pkin(5) * t91 - t18;
t13 = -mrSges(7,1) * t52 + mrSges(7,2) * t51;
t10 = t51 * Ifges(7,1) + t52 * Ifges(7,4) + t92 * Ifges(7,5);
t9 = t51 * Ifges(7,4) + t52 * Ifges(7,2) + t92 * Ifges(7,6);
t8 = t132 * t23 + t134 * t36;
t7 = -t132 * t36 + t134 * t23;
t4 = -qJD(6) * t12 - t132 * t22 + t134 * t30;
t3 = qJD(6) * t11 + t132 * t30 + t134 * t22;
t25 = [-t244 * t44 + (Ifges(7,4) * t188 + Ifges(7,2) * t187) * t234 + t39 * t67 + t3 * t53 + t4 * t54 + t31 * t43 + t48 * t13 + t11 * t28 + t12 * t29 + t56 * (-t91 * mrSges(6,2) - t89) + (Ifges(7,1) * t188 + Ifges(7,4) * t187) * t232 + (Ifges(7,5) * t188 + Ifges(7,6) * t187) * t230 + (t57 / 0.2e1 - t59 / 0.2e1 - t38 * mrSges(5,3) + t35 * mrSges(6,1) - Ifges(5,4) * t225 + Ifges(6,6) * t226 + Ifges(6,3) * t227 - Ifges(5,2) * t228 + t246 * t214 + t239) * t101 + (t18 * mrSges(6,1) - t21 * mrSges(5,3) - t33 * mrSges(6,2) - t14 * t154 + t151 * t237 + t149 * t236 + t132 * t10 / 0.2e1 + t9 * t222 + (t201 / 0.2e1 + t199 / 0.2e1 + t211) * t92 + (Ifges(6,3) + Ifges(5,2)) * t91 + (t148 * t230 + t150 * t234 + t152 * t232 + t153 * t19 + t222 * t27 + t223 * t26) * qJD(6)) * t108 + ((-qJD(6) * t157 + t158) * t108 + t6 * t187 - t5 * t188) * mrSges(7,3) + t212 * (t64 * t92 - t65 * t91) + (t90 / 0.2e1 - t33 * mrSges(6,3) + t211 * t91 + t212 * t20 + (Ifges(5,1) + Ifges(6,2) + Ifges(7,3) / 0.2e1) * t92 + t241) * t109 + m(7) * (t1 * t12 + t11 * t2 + t14 * t48 + t19 * t31 + t3 * t6 + t4 * t5) + (t168 * t87 + (-t94 * mrSges(4,3) - t104 / 0.2e1 + t113 * mrSges(4,1) - t189 / 0.2e1 + (-m(4) * t94 - t114) * t122 + (t170 * mrSges(4,1) - 0.3e1 / 0.2e1 * t206 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t135) * qJD(1) + (t66 + qJD(1) * (mrSges(5,1) * t108 + mrSges(5,2) * t109) + m(5) * (t99 + t185)) * pkin(3)) * qJD(3)) * t133 + (t245 / 0.2e1 - t58 / 0.2e1 - t37 * mrSges(5,3) + t34 * mrSges(6,1) + Ifges(5,1) * t225 - Ifges(6,2) * t226 - Ifges(6,6) * t227 + Ifges(5,4) * t228 + Ifges(7,3) * t230 + Ifges(7,5) * t232 + Ifges(7,6) * t234 + t247 * t214 + t238) * t103 + t110 * (t91 * mrSges(5,1) + t88) + t26 * t187 / 0.2e1 + t27 * t188 / 0.2e1 + t19 * (-mrSges(7,1) * t187 + mrSges(7,2) * t188) + (t168 * t86 + (-t93 * mrSges(4,3) + t105 / 0.2e1 + 0.3e1 / 0.2e1 * t125 + t190 / 0.2e1 + (-m(4) * t93 - t112) * t122 + 0.2e1 * t113 * mrSges(4,2)) * qJD(3)) * t135 + t210 * t45 + m(6) * (-t18 * t65 + t33 * t56 + t34 * t44 - t35 * t45 + t39 * t46 + t217) + m(5) * (t21 * t65 - t37 * t44 + t38 * t45 + t217); (-t212 * t91 + t13) * t109 + t175 * t103 + t141 * t101 + (-t133 * t112 + t135 * t114 + (-t133 ^ 2 - t135 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (t212 * t92 + t137) * t108 + m(4) * (t86 * t133 - t135 * t87 + (-t133 * t93 + t135 * t94) * qJD(3)) + m(5) * (-t101 * t37 + t103 * t38 + t109 * t21 + t194) + m(6) * (t101 * t34 - t103 * t35 - t109 * t18 + t194) + m(7) * (t101 * t157 + t103 * t19 + t108 * t138 + t109 * t14); t14 * t153 - t86 * mrSges(4,2) - t87 * mrSges(4,1) - t55 * t67 - t8 * t53 - t7 * t54 - t24 * t43 - t18 * mrSges(6,3) - t21 * mrSges(5,2) + t240 * qJD(6) + t176 * Ifges(4,5) * t135 / 0.2e1 + (-t120 * t18 + t123 * t20 - t34 * t40 - t46 * t55 + (t41 - qJD(5)) * t35) * m(6) + t243 * qJD(5) + t244 * t40 - t210 * t41 + (-t97 + t245) * t227 + (-mrSges(6,1) * t120 - mrSges(5,3) * t220 + t246) * t91 + (-Ifges(5,2) * t227 + Ifges(6,3) * t228 + t147 * t231 + t149 * t235 + t151 * t233 - t207 * t6 + t208 * t5 + t215 * t246 - t239 + t240) * t102 + (t148 / 0.2e1 - mrSges(5,3) * t169 + t247) * t92 + (t96 + t58) * t228 + (-Ifges(5,1) * t226 - Ifges(7,5) * t233 + Ifges(6,2) * t225 - Ifges(7,6) * t235 - Ifges(7,3) * t231 - t215 * t247 + t238) * t100 + t248 * t20 - (-Ifges(4,2) * t184 + t105 + t125) * t183 / 0.2e1 - (t147 * t98 + t149 * t76 + t151 * t77) * qJD(6) / 0.2e1 + (t174 + t112) * t94 + (-t113 * (mrSges(4,1) * t133 + mrSges(4,2) * t135) - (Ifges(4,1) * t135 - t206) * t184 / 0.2e1) * qJD(1) + (m(7) * t138 + t179 * t53 - t180 * t54 + t242) * (-pkin(8) + t123) + (t120 * t14 - t5 * t7 - t6 * t8 + (-t24 + qJD(5)) * t19) * m(7) + (t173 - t114) * t93 + t123 * t216 + t10 * t222 + t9 * t223 + t38 * t195 + t34 * t198 + t150 * t236 + t152 * t237 + (-t179 * t6 + t180 * t5) * mrSges(7,3) + (-t205 + t57) * t226 + (t200 + t59) * t225 + t120 * t13 + ((t129 * t21 - t191 * t20) * pkin(3) - t126 * t99 + t37 * t40 - t38 * t41) * m(5) - Ifges(4,6) * t166 / 0.2e1 - t66 * t126 + t104 * t184 / 0.2e1 - t35 * t196 - t37 * t197 - t2 * t207 - t1 * t208; -t132 * t28 + t134 * t29 + t88 - t89 - t248 * t91 + t145 * qJD(6) + t175 * t100 - t141 * t102 + (t100 * t19 - t157 * t98 + t158) * m(7) + (-t100 * t35 - t102 * t34 + t33) * m(6) + (t100 * t38 + t102 * t37 + t119) * m(5); t216 - t243 * qJD(3) + (t146 + t67) * t102 + t137 + (-qJD(3) * t19 - t102 * t156 + t138) * m(7) + (qJD(3) * t35 + t102 * t46 + t20) * m(6); t90 - t19 * (mrSges(7,1) * t77 + mrSges(7,2) * t76) + (Ifges(7,1) * t76 - t221) * t233 + t26 * t232 + (Ifges(7,5) * t76 - Ifges(7,6) * t77) * t231 - t5 * t53 + t6 * t54 + (t5 * t76 + t6 * t77) * mrSges(7,3) + (-Ifges(7,2) * t77 + t27 + t74) * t235 + t241;];
tauc  = t25(:);
