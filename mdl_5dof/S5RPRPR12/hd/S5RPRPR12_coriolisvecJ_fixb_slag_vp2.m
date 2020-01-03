% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:27
% EndTime: 2019-12-31 18:29:36
% DurationCPUTime: 4.15s
% Computational Cost: add. (4239->379), mult. (11509->536), div. (0->0), fcn. (8610->8), ass. (0->175)
t244 = Ifges(4,1) / 0.2e1;
t155 = cos(pkin(8));
t214 = cos(qJ(3));
t176 = t214 * t155;
t146 = qJD(1) * t176;
t153 = sin(pkin(8));
t157 = sin(qJ(3));
t183 = t157 * t153;
t127 = qJD(1) * t183 - t146;
t243 = -t127 / 0.2e1;
t242 = Ifges(5,3) + Ifges(4,2);
t137 = t153 * t214 + t157 * t155;
t128 = t137 * qJD(1);
t152 = sin(pkin(9));
t154 = cos(pkin(9));
t113 = qJD(3) * t152 + t128 * t154;
t156 = sin(qJ(5));
t158 = cos(qJ(5));
t173 = t154 * qJD(3) - t128 * t152;
t241 = -t113 * t156 + t158 * t173;
t64 = t113 * t158 + t156 * t173;
t163 = t176 - t183;
t131 = t163 * qJD(3);
t121 = qJD(1) * t131;
t164 = t152 * t156 - t154 * t158;
t32 = qJD(5) * t241 - t121 * t164;
t227 = t32 / 0.2e1;
t136 = t152 * t158 + t154 * t156;
t33 = -qJD(5) * t64 - t121 * t136;
t226 = t33 / 0.2e1;
t132 = t137 * qJD(3);
t122 = qJD(1) * t132;
t219 = t122 / 0.2e1;
t240 = Ifges(4,4) * t243;
t239 = t127 / 0.2e1;
t238 = t128 * t244;
t208 = pkin(7) + qJ(4);
t141 = t208 * t152;
t143 = t208 * t154;
t109 = -t141 * t156 + t143 * t158;
t210 = pkin(7) * t154;
t102 = pkin(3) * t128 + qJ(4) * t127;
t209 = pkin(6) + qJ(2);
t144 = t209 * t155;
t139 = qJD(1) * t144;
t126 = t157 * t139;
t142 = t209 * t153;
t138 = qJD(1) * t142;
t177 = t214 * t138;
t104 = -t126 - t177;
t46 = t154 * t102 - t104 * t152;
t28 = pkin(4) * t128 + t127 * t210 + t46;
t187 = t127 * t152;
t47 = t152 * t102 + t154 * t104;
t37 = pkin(7) * t187 + t47;
t237 = -qJD(4) * t136 - qJD(5) * t109 + t156 * t37 - t158 * t28;
t107 = -t141 * t158 - t143 * t156;
t236 = -qJD(4) * t164 + qJD(5) * t107 - t156 * t28 - t158 * t37;
t130 = t136 * qJD(5);
t86 = t136 * t127;
t193 = -t86 - t130;
t129 = t164 * qJD(5);
t87 = t164 * t127;
t192 = -t87 - t129;
t233 = -t214 * t142 - t157 * t144;
t188 = t121 * t154;
t60 = pkin(3) * t122 - qJ(4) * t121 - qJD(4) * t128;
t180 = qJD(1) * qJD(2);
t175 = t153 * t180;
t182 = qJD(2) * t146 - qJD(3) * t177;
t66 = -t157 * t175 + (qJD(4) - t126) * qJD(3) + t182;
t24 = -t152 * t66 + t154 * t60;
t14 = pkin(4) * t122 - pkin(7) * t188 + t24;
t189 = t121 * t152;
t25 = t152 * t60 + t154 * t66;
t15 = -pkin(7) * t189 + t25;
t105 = -t157 * t138 + t139 * t214;
t101 = qJD(3) * qJ(4) + t105;
t178 = -pkin(2) * t155 - pkin(1);
t140 = qJD(1) * t178 + qJD(2);
t78 = pkin(3) * t127 - qJ(4) * t128 + t140;
t38 = -t101 * t152 + t154 * t78;
t20 = pkin(4) * t127 - pkin(7) * t113 + t38;
t39 = t154 * t101 + t152 * t78;
t27 = pkin(7) * t173 + t39;
t5 = -t156 * t27 + t158 * t20;
t1 = qJD(5) * t5 + t14 * t156 + t15 * t158;
t6 = t156 * t20 + t158 * t27;
t2 = -qJD(5) * t6 + t14 * t158 - t15 * t156;
t232 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t32 + Ifges(6,6) * t33;
t231 = (m(3) * qJ(2) + mrSges(3,3)) * (t153 ^ 2 + t155 ^ 2);
t201 = Ifges(5,2) * t152;
t204 = Ifges(5,4) * t154;
t168 = -t201 + t204;
t205 = Ifges(5,4) * t152;
t169 = Ifges(5,1) * t154 - t205;
t170 = mrSges(5,1) * t152 + mrSges(5,2) * t154;
t215 = t154 / 0.2e1;
t216 = -t152 / 0.2e1;
t99 = -qJD(3) * pkin(3) + qJD(4) - t104;
t230 = t168 * t173 / 0.2e1 + t169 * t113 / 0.2e1 + t99 * t170 + t140 * mrSges(4,2) + (-t39 * t152 - t38 * t154) * mrSges(5,3) + t240 + Ifges(4,5) * qJD(3) + t238 + (t113 * Ifges(5,4) + Ifges(5,2) * t173 + t127 * Ifges(5,6)) * t216 + (t113 * Ifges(5,1) + Ifges(5,4) * t173 + t127 * Ifges(5,5)) * t215;
t229 = Ifges(6,4) * t227 + Ifges(6,2) * t226 + Ifges(6,6) * t219;
t228 = Ifges(6,1) * t227 + Ifges(6,4) * t226 + Ifges(6,5) * t219;
t225 = -t241 / 0.2e1;
t224 = t241 / 0.2e1;
t223 = -t64 / 0.2e1;
t222 = t64 / 0.2e1;
t125 = qJD(5) + t127;
t218 = -t125 / 0.2e1;
t217 = t125 / 0.2e1;
t213 = Ifges(6,4) * t64;
t73 = pkin(3) * t132 - qJ(4) * t131 - qJD(4) * t137;
t81 = t163 * qJD(2) + t233 * qJD(3);
t35 = t152 * t73 + t154 * t81;
t206 = Ifges(4,4) * t128;
t202 = Ifges(5,5) * t154;
t199 = Ifges(5,6) * t152;
t162 = t137 * qJD(2);
t71 = qJD(1) * t162 + qJD(3) * t105;
t197 = t233 * t71;
t196 = t122 * Ifges(5,5);
t195 = t122 * Ifges(5,6);
t103 = -pkin(3) * t163 - qJ(4) * t137 + t178;
t110 = -t157 * t142 + t144 * t214;
t52 = t152 * t103 + t154 * t110;
t194 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t173 + mrSges(5,2) * t113 + t128 * mrSges(4,3);
t186 = t131 * t152;
t185 = t137 * t152;
t77 = mrSges(5,1) * t189 + mrSges(5,2) * t188;
t179 = Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t9 = -t33 * mrSges(6,1) + t32 * mrSges(6,2);
t34 = -t152 * t81 + t154 * t73;
t174 = t122 * mrSges(4,1) + t121 * mrSges(4,2);
t51 = t154 * t103 - t110 * t152;
t167 = -t199 + t202;
t166 = t25 * t152 + t24 * t154;
t165 = t152 * t38 - t154 * t39;
t36 = -pkin(4) * t163 - t137 * t210 + t51;
t40 = -pkin(7) * t185 + t52;
t12 = -t156 * t40 + t158 * t36;
t13 = t156 * t36 + t158 * t40;
t82 = qJD(3) * t110 + t162;
t160 = t140 * mrSges(4,1) + t38 * mrSges(5,1) + t5 * mrSges(6,1) + Ifges(6,3) * t125 + Ifges(6,6) * t241 + Ifges(6,5) * t64 + Ifges(5,6) * t173 + Ifges(5,5) * t113 - Ifges(4,6) * qJD(3) - t206 / 0.2e1 - t39 * mrSges(5,2) - t6 * mrSges(6,2) + t242 * t239;
t148 = -pkin(4) * t154 - pkin(3);
t120 = Ifges(6,3) * t122;
t117 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t127;
t94 = t164 * t137;
t93 = t136 * t137;
t89 = mrSges(5,1) * t122 - mrSges(5,3) * t188;
t88 = -mrSges(5,2) * t122 - mrSges(5,3) * t189;
t85 = pkin(4) * t185 - t233;
t84 = mrSges(5,1) * t127 - mrSges(5,3) * t113;
t83 = -mrSges(5,2) * t127 + mrSges(5,3) * t173;
t72 = -pkin(4) * t187 + t105;
t70 = (-qJD(3) * t139 - t175) * t157 + t182;
t61 = Ifges(6,4) * t241;
t56 = -pkin(4) * t173 + t99;
t55 = pkin(4) * t186 + t82;
t54 = t121 * t169 + t196;
t53 = t121 * t168 + t195;
t45 = pkin(4) * t189 + t71;
t44 = mrSges(6,1) * t125 - mrSges(6,3) * t64;
t43 = -mrSges(6,2) * t125 + mrSges(6,3) * t241;
t42 = t137 * t129 - t131 * t136;
t41 = -t130 * t137 - t131 * t164;
t26 = -mrSges(6,1) * t241 + mrSges(6,2) * t64;
t23 = -pkin(7) * t186 + t35;
t22 = -mrSges(6,2) * t122 + mrSges(6,3) * t33;
t21 = mrSges(6,1) * t122 - mrSges(6,3) * t32;
t19 = Ifges(6,1) * t64 + Ifges(6,5) * t125 + t61;
t18 = Ifges(6,2) * t241 + Ifges(6,6) * t125 + t213;
t16 = pkin(4) * t132 - t131 * t210 + t34;
t4 = -qJD(5) * t13 - t156 * t23 + t158 * t16;
t3 = qJD(5) * t12 + t156 * t16 + t158 * t23;
t7 = [m(5) * (t24 * t51 + t25 * t52 + t34 * t38 + t35 * t39 + t82 * t99 - t197) + m(4) * (-t104 * t82 + t105 * t81 + t110 * t70 - t197) + 0.2e1 * t231 * t180 + t178 * t174 + (t179 * t127 + t160) * t132 + (-t104 * t131 - t105 * t132 - t110 * t122 - t121 * t233 + t137 * t71 + t163 * t70) * mrSges(4,3) - t233 * t77 + (Ifges(6,5) * t41 + Ifges(6,6) * t42) * t217 + (Ifges(6,1) * t41 + Ifges(6,4) * t42) * t222 + (Ifges(6,4) * t41 + Ifges(6,2) * t42) * t224 - t94 * t228 - t93 * t229 + (t167 * t239 + t230 + t238) * t131 + m(6) * (t1 * t13 + t12 * t2 + t3 * t6 + t4 * t5 + t45 * t85 + t55 * t56) + (-Ifges(6,5) * t94 - Ifges(6,6) * t93) * t219 + (-Ifges(6,4) * t94 - Ifges(6,2) * t93) * t226 + (-Ifges(6,1) * t94 - Ifges(6,4) * t93) * t227 + (-t1 * t93 + t2 * t94 - t41 * t5 + t42 * t6) * mrSges(6,3) + t45 * (mrSges(6,1) * t93 - mrSges(6,2) * t94) + (t131 * t243 - t128 * t132 / 0.2e1 + t163 * t121 - t137 * t122) * Ifges(4,4) - (t120 / 0.2e1 - t25 * mrSges(5,2) + t24 * mrSges(5,1) + t167 * t121 + (Ifges(6,3) / 0.2e1 + t242) * t122 + t232) * t163 + t194 * t82 + (t53 * t216 + t54 * t215 + t71 * t170 + t167 * t219 - t166 * mrSges(5,3) + (Ifges(4,1) + Ifges(5,1) * t154 ^ 2 / 0.2e1 + (-t204 + t201 / 0.2e1) * t152) * t121) * t137 + t12 * t21 + t13 * t22 + t41 * t19 / 0.2e1 + t42 * t18 / 0.2e1 + t3 * t43 + t4 * t44 + t55 * t26 + t56 * (-mrSges(6,1) * t42 + mrSges(6,2) * t41) + t35 * t83 + t34 * t84 + t85 * t9 + t52 * t88 + t51 * t89 + t81 * t117; -t164 * t21 + t136 * t22 + t152 * t88 + t154 * t89 + t193 * t44 + t192 * t43 + (-t26 - t194) * t128 + (-t152 * t84 + t154 * t83 + t117) * t127 - m(4) * (-t104 * t128 - t105 * t127) + t174 - t231 * qJD(1) ^ 2 + (t1 * t136 - t128 * t56 - t164 * t2 + t192 * t6 + t193 * t5) * m(6) + (-t127 * t165 - t99 * t128 + t166) * m(5); -t164 * t229 + (Ifges(6,5) * t136 - Ifges(6,6) * t164) * t219 + (Ifges(6,4) * t136 - Ifges(6,2) * t164) * t226 + (Ifges(6,1) * t136 - Ifges(6,4) * t164) * t227 + (-t1 * t164 - t136 * t2 - t192 * t5 + t193 * t6) * mrSges(6,3) + t45 * (mrSges(6,1) * t164 + mrSges(6,2) * t136) + (Ifges(6,5) * t87 + Ifges(6,6) * t86) * t218 + (Ifges(6,1) * t87 + Ifges(6,4) * t86) * t223 + (Ifges(6,4) * t87 + Ifges(6,2) * t86) * t225 + t136 * t228 + (-t104 * mrSges(4,3) + t240 + (t202 / 0.2e1 - t199 / 0.2e1) * t127 + (t244 - t179) * t128 + t230) * t127 + t236 * t43 + (t1 * t109 + t107 * t2 + t148 * t45 + t236 * t6 + t237 * t5 - t56 * t72) * m(6) + t237 * t44 + (-pkin(3) * t71 - t165 * qJD(4) + (-t152 * t24 + t154 * t25) * qJ(4) - t105 * t99 - t38 * t46 - t39 * t47) * m(5) + (-t87 / 0.2e1 - t129 / 0.2e1) * t19 + (-t86 / 0.2e1 - t130 / 0.2e1) * t18 + (-Ifges(6,5) * t129 - Ifges(6,6) * t130) * t217 + (-Ifges(6,1) * t129 - Ifges(6,4) * t130) * t222 + (-Ifges(6,4) * t129 - Ifges(6,2) * t130) * t224 + (-mrSges(6,1) * t193 + mrSges(6,2) * t192) * t56 - t194 * t105 + (t25 * mrSges(5,3) + qJD(4) * t83 + qJ(4) * t88 + t53 / 0.2e1 - t71 * mrSges(5,1) + t195 / 0.2e1) * t154 + (-t24 * mrSges(5,3) - qJD(4) * t84 - qJ(4) * t89 + t54 / 0.2e1 + t71 * mrSges(5,2) + t196 / 0.2e1) * t152 + (-t160 + t105 * mrSges(4,3) + t206 / 0.2e1) * t128 + ((Ifges(5,1) * t152 + t204) * t215 + (Ifges(5,2) * t154 + t205) * t216 + Ifges(4,5)) * t121 - t70 * mrSges(4,2) - t71 * mrSges(4,1) - t72 * t26 - pkin(3) * t77 - t47 * t83 - t46 * t84 + t107 * t21 + t109 * t22 - t104 * t117 - Ifges(4,6) * t122 + t148 * t9; t113 * t84 - t173 * t83 - t241 * t43 + t64 * t44 + t77 + t9 + (-t241 * t6 + t5 * t64 + t45) * m(6) + (t113 * t38 - t173 * t39 + t71) * m(5); t120 - t56 * (mrSges(6,1) * t64 + mrSges(6,2) * t241) + (Ifges(6,1) * t241 - t213) * t223 + t18 * t222 + (Ifges(6,5) * t241 - Ifges(6,6) * t64) * t218 - t5 * t43 + t6 * t44 + (t241 * t5 + t6 * t64) * mrSges(6,3) + (-Ifges(6,2) * t64 + t19 + t61) * t225 + t232;];
tauc = t7(:);
