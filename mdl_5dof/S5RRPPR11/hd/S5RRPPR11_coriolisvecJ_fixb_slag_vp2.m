% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:26
% EndTime: 2019-12-31 19:46:36
% DurationCPUTime: 4.18s
% Computational Cost: add. (2725->399), mult. (6628->565), div. (0->0), fcn. (3919->6), ass. (0->191)
t160 = cos(qJ(2));
t233 = pkin(3) + pkin(6);
t258 = t233 * t160;
t201 = qJD(1) * t160;
t257 = -Ifges(3,4) * t201 / 0.2e1;
t251 = -qJD(1) / 0.2e1;
t256 = -mrSges(3,1) + mrSges(4,2);
t157 = sin(qJ(5));
t159 = cos(qJ(5));
t154 = sin(pkin(8));
t158 = sin(qJ(2));
t205 = t154 * t158;
t170 = pkin(4) * t160 - pkin(7) * t205;
t166 = t170 * qJD(2);
t155 = cos(pkin(8));
t196 = qJD(1) * qJD(2);
t187 = t158 * t196;
t138 = pkin(2) * t187;
t173 = -qJ(3) * t160 + qJ(4) * t158;
t197 = t158 * qJD(3);
t163 = qJD(2) * t173 - qJD(4) * t160 - t197;
t54 = qJD(1) * t163 + t138;
t90 = (qJD(1) * t258 - qJD(4)) * qJD(2);
t28 = -t154 * t54 + t155 * t90;
t17 = qJD(1) * t166 + t28;
t183 = t155 * t187;
t29 = t154 * t90 + t155 * t54;
t19 = pkin(7) * t183 + t29;
t198 = t155 * qJD(2);
t110 = t154 * t201 - t198;
t202 = qJD(1) * t158;
t156 = -pkin(2) - qJ(4);
t185 = -t158 * qJ(3) - pkin(1);
t109 = t156 * t160 + t185;
t80 = t109 * qJD(1);
t144 = pkin(6) * t202;
t118 = -pkin(3) * t202 - t144;
t245 = -t118 + qJD(3);
t85 = qJD(2) * t156 + t245;
t36 = -t154 * t80 + t155 * t85;
t22 = pkin(4) * t202 + pkin(7) * t110 + t36;
t111 = -t154 * qJD(2) - t155 * t201;
t37 = t154 * t85 + t155 * t80;
t23 = pkin(7) * t111 + t37;
t5 = -t157 * t23 + t159 * t22;
t1 = qJD(5) * t5 + t157 * t17 + t159 * t19;
t172 = t159 * t154 + t157 * t155;
t6 = t157 * t22 + t159 * t23;
t2 = -qJD(5) * t6 - t157 * t19 + t159 * t17;
t167 = t172 * t158;
t79 = qJD(1) * t167;
t99 = t172 * qJD(5);
t218 = -t99 - t79;
t188 = t155 * t202;
t206 = t154 * t157;
t78 = t159 * t188 - t202 * t206;
t246 = t155 * t159 - t206;
t98 = t246 * qJD(5);
t219 = t98 + t78;
t255 = -t1 * t172 - t2 * t246 - t218 * t5 - t219 * t6;
t250 = -qJD(2) / 0.2e1;
t249 = qJD(2) / 0.2e1;
t146 = pkin(6) * t201;
t119 = pkin(3) * t201 + t146;
t145 = pkin(2) * t202;
t91 = qJD(1) * t173 + t145;
t48 = t155 * t119 - t154 * t91;
t33 = qJD(1) * t170 + t48;
t49 = t154 * t119 + t155 * t91;
t40 = pkin(7) * t188 + t49;
t220 = -pkin(7) + t156;
t122 = t220 * t154;
t123 = t220 * t155;
t64 = -t122 * t157 + t123 * t159;
t248 = -qJD(4) * t172 + qJD(5) * t64 - t157 * t33 - t159 * t40;
t65 = t122 * t159 + t123 * t157;
t247 = -qJD(4) * t246 - qJD(5) * t65 + t157 * t40 - t159 * t33;
t184 = t110 * t157 + t159 * t111;
t165 = qJD(2) * t167;
t30 = qJD(1) * t165 + qJD(5) * t184;
t200 = qJD(2) * t158;
t164 = t246 * t200;
t55 = t110 * t159 - t111 * t157;
t31 = qJD(1) * t164 + qJD(5) * t55;
t244 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t30 + Ifges(6,6) * t31;
t126 = -pkin(2) * t160 + t185;
t108 = t126 * qJD(1);
t125 = -qJD(2) * pkin(2) + qJD(3) + t144;
t139 = qJD(5) + t202;
t212 = Ifges(4,6) * t160;
t243 = -(m(4) * t125 + (mrSges(4,1) + mrSges(3,3)) * t202 + t256 * qJD(2)) * pkin(6) - t125 * mrSges(4,1) - t36 * mrSges(5,1) - t5 * mrSges(6,1) - Ifges(3,5) * t249 + t257 - Ifges(4,4) * t250 - (-Ifges(4,2) * t158 - t212) * t251 - t139 * Ifges(6,3) + t55 * Ifges(6,5) - t184 * Ifges(6,6) - t111 * Ifges(5,6) + t110 * Ifges(5,5) + t108 * mrSges(4,3) + t37 * mrSges(5,2) + t6 * mrSges(6,2) - (Ifges(3,1) + Ifges(5,3)) * t202 / 0.2e1;
t153 = qJD(2) * qJ(3);
t130 = -t146 - t153;
t132 = -mrSges(4,1) * t201 - qJD(2) * mrSges(4,3);
t174 = t36 * t154 - t37 * t155;
t215 = Ifges(5,4) * t154;
t176 = Ifges(5,2) * t155 + t215;
t214 = Ifges(5,4) * t155;
t216 = Ifges(5,1) * t154;
t177 = t214 + t216;
t217 = mrSges(5,2) * t154;
t178 = mrSges(5,1) * t155 - t217;
t225 = -t155 / 0.2e1;
t226 = -t154 / 0.2e1;
t96 = qJD(4) + t153 + t119;
t242 = -t174 * mrSges(5,3) + (m(4) * t130 + qJD(2) * mrSges(3,2) - mrSges(3,3) * t201 + t132) * pkin(6) - (-t110 * Ifges(5,1) + t111 * Ifges(5,4) + Ifges(5,5) * t202) * t226 - (-t110 * Ifges(5,4) + t111 * Ifges(5,2) + Ifges(5,6) * t202) * t225 + t130 * mrSges(4,1) - Ifges(4,5) * t250 - Ifges(3,6) * t249 - t110 * t177 / 0.2e1 + t111 * t176 / 0.2e1 - t108 * mrSges(4,2) - t96 * t178 + ((Ifges(3,2) + Ifges(4,3)) * t160 + (Ifges(3,4) + Ifges(4,6)) * t158) * t251;
t241 = t30 / 0.2e1;
t240 = t31 / 0.2e1;
t239 = -t184 / 0.2e1;
t238 = t184 / 0.2e1;
t237 = t55 / 0.2e1;
t236 = -t55 / 0.2e1;
t86 = t246 * t160;
t235 = -t86 / 0.2e1;
t87 = t172 * t160;
t234 = -t87 / 0.2e1;
t232 = pkin(1) * mrSges(3,1);
t231 = pkin(1) * mrSges(3,2);
t230 = -t172 / 0.2e1;
t229 = t246 / 0.2e1;
t228 = -t139 / 0.2e1;
t227 = t139 / 0.2e1;
t224 = t155 / 0.2e1;
t223 = Ifges(6,4) * t55;
t213 = Ifges(5,5) * t154;
t211 = Ifges(5,6) * t155;
t121 = qJD(2) * t258;
t148 = pkin(2) * t200;
t69 = t148 + t163;
t39 = t154 * t121 + t155 * t69;
t134 = t233 * t158;
t60 = t155 * t109 + t154 * t134;
t203 = t155 * t160;
t199 = qJD(2) * t160;
t18 = -mrSges(6,1) * t184 - mrSges(6,2) * t55;
t63 = -mrSges(5,1) * t111 - mrSges(5,2) * t110;
t195 = t132 - t18 - t63;
t194 = -0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(4,6);
t193 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t192 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t191 = m(4) * pkin(6) + mrSges(4,1);
t189 = -pkin(4) * t155 - pkin(3);
t9 = -t31 * mrSges(6,1) + t30 * mrSges(6,2);
t186 = t160 * t196;
t38 = t155 * t121 - t154 * t69;
t120 = t233 * t200;
t175 = t154 * t29 + t155 * t28;
t115 = t155 * t134;
t43 = t158 * pkin(4) + t115 + (pkin(7) * t160 - t109) * t154;
t46 = -pkin(7) * t203 + t60;
t12 = -t157 * t46 + t159 * t43;
t13 = t157 * t43 + t159 * t46;
t169 = -t213 / 0.2e1 - t211 / 0.2e1;
t84 = (-pkin(6) + t189) * t200;
t168 = -qJ(3) * t199 - t197;
t82 = -mrSges(5,1) * t183 + t187 * t217;
t152 = qJD(2) * qJD(3);
t140 = pkin(4) * t154 + qJ(3);
t137 = Ifges(6,3) * t186;
t124 = pkin(6) * t187 - t152;
t116 = (mrSges(4,2) * t160 - t158 * mrSges(4,3)) * qJD(1);
t97 = pkin(4) * t203 + t258;
t94 = t148 + t168;
t93 = (mrSges(5,3) * t155 * t158 - mrSges(5,2) * t160) * t196;
t92 = (mrSges(5,1) * t160 - mrSges(5,3) * t205) * t196;
t89 = -qJD(1) * t120 + t152;
t83 = t189 * t202 - t144;
t81 = qJD(1) * t168 + t138;
t77 = mrSges(5,1) * t202 + mrSges(5,3) * t110;
t76 = -mrSges(5,2) * t202 + mrSges(5,3) * t111;
t68 = (Ifges(5,5) * t160 + t158 * t177) * t196;
t67 = (Ifges(5,6) * t160 + t158 * t176) * t196;
t66 = qJD(1) * t84 + t152;
t61 = -pkin(4) * t111 + t96;
t59 = -t154 * t109 + t115;
t53 = Ifges(6,4) * t184;
t45 = t160 * t99 + t164;
t44 = -qJD(5) * t86 + t165;
t42 = mrSges(6,1) * t139 + mrSges(6,3) * t55;
t41 = -mrSges(6,2) * t139 + mrSges(6,3) * t184;
t32 = pkin(7) * t158 * t198 + t39;
t27 = t166 + t38;
t21 = -mrSges(6,2) * t186 + t31 * mrSges(6,3);
t20 = mrSges(6,1) * t186 - t30 * mrSges(6,3);
t16 = -Ifges(6,1) * t55 + Ifges(6,5) * t139 + t53;
t15 = Ifges(6,2) * t184 + Ifges(6,6) * t139 - t223;
t8 = Ifges(6,1) * t30 + Ifges(6,4) * t31 + Ifges(6,5) * t186;
t7 = Ifges(6,4) * t30 + Ifges(6,2) * t31 + Ifges(6,6) * t186;
t4 = -qJD(5) * t13 - t157 * t32 + t159 * t27;
t3 = qJD(5) * t12 + t157 * t27 + t159 * t32;
t10 = [(t137 / 0.2e1 + t28 * mrSges(5,1) - t81 * mrSges(4,3) - t29 * mrSges(5,2) + (t192 * qJD(2) + ((0.3e1 / 0.2e1 * t213 + 0.3e1 / 0.2e1 * t211 + t194) * t158 - t126 * mrSges(4,2) - 0.2e1 * t232 + (0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + Ifges(6,3) / 0.2e1 - t155 ^ 2 * Ifges(5,2) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + t191 * pkin(6) + (-t214 - t216 / 0.2e1) * t154) * t160) * qJD(1) + t242) * qJD(2) + t244) * t158 + (t81 * mrSges(4,2) + t89 * t178 + t67 * t225 + t68 * t226 - t191 * t124 + (t154 * t28 - t155 * t29) * mrSges(5,3) + (((t169 - t194) * t160 - t126 * mrSges(4,3) - 0.2e1 * t231 + Ifges(6,5) * t234 + Ifges(6,6) * t235) * qJD(1) + t193 * qJD(2) - t243) * qJD(2)) * t160 + t94 * t116 - t120 * t63 + t97 * t9 + t59 * t92 + t60 * t93 + t84 * t18 + t39 * t76 + t38 * t77 + t61 * (-mrSges(6,1) * t45 + mrSges(6,2) * t44) + m(4) * (t108 * t94 + t81 * t126) + m(6) * (t1 * t13 + t12 * t2 + t3 * t6 + t4 * t5 + t61 * t84 + t66 * t97) + t3 * t41 + t4 * t42 + t44 * t16 / 0.2e1 + t45 * t15 / 0.2e1 + t12 * t20 + t13 * t21 + t258 * t82 + m(5) * (-t120 * t96 + t258 * t89 + t28 * t59 + t29 * t60 + t36 * t38 + t37 * t39) + t66 * (mrSges(6,1) * t86 - mrSges(6,2) * t87) + (-Ifges(6,4) * t87 - Ifges(6,2) * t86) * t240 + (-Ifges(6,1) * t87 - Ifges(6,4) * t86) * t241 + (-t1 * t86 + t2 * t87 - t44 * t5 + t45 * t6) * mrSges(6,3) + (Ifges(6,5) * t44 + Ifges(6,6) * t45) * t227 + t8 * t234 + t7 * t235 + (Ifges(6,1) * t44 + Ifges(6,4) * t45) * t236 + (Ifges(6,4) * t44 + Ifges(6,2) * t45) * t238; t255 * mrSges(6,3) + (-m(4) * t108 - t116) * (-qJ(3) * t201 + t145) + (-t99 / 0.2e1 - t79 / 0.2e1) * t16 + (-Ifges(6,5) * t99 - Ifges(6,6) * t98) * t227 + (-Ifges(6,1) * t99 - Ifges(6,4) * t98) * t236 + (-Ifges(6,4) * t99 - Ifges(6,2) * t98) * t238 + ((t257 + (t231 - t212 / 0.2e1) * qJD(1) + t243) * t160 + ((t232 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1 + t169) * t158) * qJD(1) + (-Ifges(5,3) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t201 - t242) * t158 + ((-pkin(2) * mrSges(4,1) + Ifges(5,5) * t224 + Ifges(5,6) * t226 + Ifges(6,5) * t229 + Ifges(6,6) * t230 + (-m(4) * pkin(2) + t256) * pkin(6) + t193) * t160 + (t154 * (Ifges(5,1) * t155 - t215) / 0.2e1 + (-Ifges(5,2) * t154 + t214) * t224 + pkin(6) * mrSges(3,2) - qJ(3) * mrSges(4,1) + t192) * t158) * qJD(2)) * qJD(1) + (qJ(3) * t89 + t175 * t156 + (-t154 * t37 - t155 * t36) * qJD(4) - t36 * t48 - t37 * t49 + t245 * t96) * m(5) + (mrSges(6,1) * t219 + mrSges(6,2) * t218) * t61 + t140 * t9 - t118 * t63 - t124 * mrSges(4,3) + qJ(3) * t82 - t83 * t18 + t65 * t21 - t49 * t76 - t48 * t77 + t64 * t20 + m(4) * (-t124 * qJ(3) - t130 * qJD(3)) + (-t98 / 0.2e1 - t78 / 0.2e1) * t15 + (t156 * t92 - qJD(4) * t77 - t28 * mrSges(5,3) + t89 * mrSges(5,2) + t68 / 0.2e1) * t155 + (t156 * t93 - qJD(4) * t76 - t29 * mrSges(5,3) + t89 * mrSges(5,1) - t67 / 0.2e1) * t154 + t66 * (mrSges(6,1) * t172 + mrSges(6,2) * t246) + (Ifges(6,4) * t246 - Ifges(6,2) * t172) * t240 + (Ifges(6,1) * t246 - Ifges(6,4) * t172) * t241 + t247 * t42 + t248 * t41 + (t1 * t65 + t140 * t66 + t2 * t64 + (-t83 + qJD(3)) * t61 + t248 * t6 + t247 * t5) * m(6) - t195 * qJD(3) + (Ifges(6,5) * t79 + Ifges(6,6) * t78) * t228 + t8 * t229 + t7 * t230 + (Ifges(6,1) * t79 + Ifges(6,4) * t78) * t237 + (Ifges(6,4) * t79 + Ifges(6,2) * t78) * t239; t172 * t21 + t246 * t20 + t154 * t93 + t155 * t92 + t218 * t42 + t219 * t41 + t195 * qJD(2) + (t191 * t199 + (-t154 * t77 + t155 * t76 + t116) * t158) * qJD(1) - m(4) * (-qJD(2) * t130 - t108 * t202) + (-qJD(2) * t61 - t255) * m(6) + (-qJD(2) * t96 - t174 * t202 + t175) * m(5); -t110 * t77 - t111 * t76 - t184 * t41 - t55 * t42 + t82 + t9 + (-t184 * t6 - t5 * t55 + t66) * m(6) + (-t110 * t36 - t111 * t37 + t89) * m(5); t137 - t61 * (-mrSges(6,1) * t55 + mrSges(6,2) * t184) + (Ifges(6,1) * t184 + t223) * t237 + t15 * t236 + (Ifges(6,5) * t184 + Ifges(6,6) * t55) * t228 - t5 * t41 + t6 * t42 + (t184 * t5 - t55 * t6) * mrSges(6,3) + (Ifges(6,2) * t55 + t16 + t53) * t239 + t244;];
tauc = t10(:);
