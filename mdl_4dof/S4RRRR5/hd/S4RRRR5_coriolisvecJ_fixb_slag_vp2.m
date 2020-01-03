% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:23
% EndTime: 2019-12-31 17:27:31
% DurationCPUTime: 3.90s
% Computational Cost: add. (2710->363), mult. (6985->536), div. (0->0), fcn. (4368->6), ass. (0->186)
t242 = qJD(2) / 0.2e1;
t132 = sin(qJ(3));
t214 = -pkin(7) - pkin(6);
t164 = qJD(3) * t214;
t133 = sin(qJ(2));
t135 = cos(qJ(3));
t174 = t133 * t135;
t136 = cos(qJ(2));
t175 = t132 * t136;
t155 = pkin(2) * t133 - pkin(6) * t136;
t107 = t155 * qJD(1);
t92 = t132 * t107;
t241 = t132 * t164 - t92 - (-pkin(5) * t174 - pkin(7) * t175) * qJD(1);
t173 = t135 * t136;
t144 = pkin(3) * t133 - pkin(7) * t173;
t172 = qJD(1) * t133;
t163 = t132 * t172;
t73 = pkin(5) * t163 + t135 * t107;
t240 = -qJD(1) * t144 + t135 * t164 - t73;
t161 = Ifges(3,5) * t242;
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t170 = qJD(2) * t135;
t103 = -t163 + t170;
t166 = t136 * qJD(1);
t130 = pkin(5) * t166;
t117 = qJD(2) * pkin(6) + t130;
t112 = -pkin(2) * t136 - t133 * pkin(6) - pkin(1);
t96 = t112 * qJD(1);
t62 = t117 * t135 + t132 * t96;
t44 = pkin(7) * t103 + t62;
t181 = t134 * t44;
t125 = qJD(3) - t166;
t162 = t135 * t172;
t104 = qJD(2) * t132 + t162;
t61 = -t117 * t132 + t135 * t96;
t43 = -pkin(7) * t104 + t61;
t39 = pkin(3) * t125 + t43;
t10 = t131 * t39 + t181;
t171 = qJD(2) * t133;
t159 = qJD(1) * t171;
t123 = Ifges(5,3) * t159;
t158 = t134 * t103 - t104 * t131;
t58 = t103 * t131 + t104 * t134;
t201 = Ifges(5,4) * t58;
t120 = qJD(4) + t125;
t206 = -t120 / 0.2e1;
t220 = -t58 / 0.2e1;
t222 = -t158 / 0.2e1;
t51 = Ifges(5,4) * t158;
t23 = Ifges(5,1) * t58 + Ifges(5,5) * t120 + t51;
t116 = -qJD(2) * pkin(2) + pkin(5) * t172;
t77 = -pkin(3) * t103 + t116;
t182 = t131 * t44;
t9 = t134 * t39 - t182;
t239 = t123 + (Ifges(5,5) * t158 - Ifges(5,6) * t58) * t206 + (t10 * t58 + t158 * t9) * mrSges(5,3) + (-Ifges(5,2) * t58 + t23 + t51) * t222 - t77 * (mrSges(5,1) * t58 + mrSges(5,2) * t158) + (Ifges(5,1) * t158 - t201) * t220;
t118 = t214 * t132;
t119 = t214 * t135;
t69 = t118 * t134 + t119 * t131;
t238 = qJD(4) * t69 + t240 * t131 + t241 * t134;
t70 = t118 * t131 - t119 * t134;
t237 = -qJD(4) * t70 - t241 * t131 + t240 * t134;
t165 = qJD(2) * qJD(3);
t168 = qJD(3) * t132;
t169 = qJD(2) * t136;
t71 = t135 * t165 + (-t133 * t168 + t135 * t169) * qJD(1);
t167 = qJD(3) * t135;
t140 = t132 * t169 + t133 * t167;
t72 = -qJD(1) * t140 - t132 * t165;
t18 = qJD(4) * t158 + t131 * t72 + t134 * t71;
t19 = -qJD(4) * t58 - t131 * t71 + t134 * t72;
t157 = pkin(5) * t159;
t110 = t155 * qJD(2);
t97 = qJD(1) * t110;
t32 = -qJD(3) * t62 + t132 * t157 + t135 * t97;
t15 = pkin(3) * t159 - pkin(7) * t71 + t32;
t31 = -t117 * t168 + t132 * t97 - t135 * t157 + t96 * t167;
t20 = pkin(7) * t72 + t31;
t2 = qJD(4) * t9 + t131 * t15 + t134 * t20;
t3 = -qJD(4) * t10 - t131 * t20 + t134 * t15;
t236 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t18 + Ifges(5,6) * t19;
t128 = Ifges(3,4) * t166;
t146 = t132 * t62 + t135 * t61;
t189 = Ifges(4,4) * t135;
t150 = -Ifges(4,2) * t132 + t189;
t190 = Ifges(4,4) * t132;
t152 = Ifges(4,1) * t135 - t190;
t153 = mrSges(4,1) * t132 + mrSges(4,2) * t135;
t187 = Ifges(4,6) * t132;
t188 = Ifges(4,5) * t135;
t202 = t135 / 0.2e1;
t203 = -t132 / 0.2e1;
t207 = t104 / 0.2e1;
t191 = Ifges(4,4) * t104;
t48 = Ifges(4,2) * t103 + Ifges(4,6) * t125 + t191;
t100 = Ifges(4,4) * t103;
t49 = Ifges(4,1) * t104 + Ifges(4,5) * t125 + t100;
t137 = -t146 * mrSges(4,3) + t49 * t202 + t48 * t203 + t125 * (-t187 + t188) / 0.2e1 + t116 * t153 + t103 * t150 / 0.2e1 + t152 * t207;
t235 = t137 + Ifges(3,1) * t172 / 0.2e1 + t128 / 0.2e1 + t161;
t22 = Ifges(5,2) * t158 + Ifges(5,6) * t120 + t201;
t233 = t22 / 0.2e1;
t160 = -Ifges(3,6) * qJD(2) / 0.2e1;
t26 = -mrSges(5,1) * t158 + mrSges(5,2) * t58;
t228 = m(5) * t77 + t26;
t145 = t131 * t132 - t134 * t135;
t88 = t145 * t133;
t126 = pkin(5) * t173;
t80 = t132 * t112 + t126;
t227 = qJD(3) + qJD(4);
t225 = -t32 * mrSges(4,1) + t31 * mrSges(4,2) - Ifges(4,5) * t71 - Ifges(4,6) * t72 - t236;
t224 = t18 / 0.2e1;
t223 = t19 / 0.2e1;
t221 = t158 / 0.2e1;
t219 = t58 / 0.2e1;
t218 = t71 / 0.2e1;
t217 = t72 / 0.2e1;
t106 = t131 * t135 + t132 * t134;
t87 = t106 * t133;
t216 = -t87 / 0.2e1;
t215 = -t88 / 0.2e1;
t212 = pkin(1) * mrSges(3,1);
t211 = pkin(1) * mrSges(3,2);
t209 = -t103 / 0.2e1;
t208 = -t104 / 0.2e1;
t205 = t120 / 0.2e1;
t204 = -t125 / 0.2e1;
t200 = pkin(3) * t132;
t199 = pkin(5) * t132;
t198 = pkin(7) * t133;
t63 = t227 * t145;
t142 = t145 * t136;
t82 = qJD(1) * t142;
t194 = -t63 + t82;
t64 = t227 * t106;
t143 = t106 * t136;
t81 = qJD(1) * t143;
t193 = -t64 + t81;
t192 = Ifges(3,4) * t133;
t180 = t135 * t110 + t171 * t199;
t177 = qJD(2) * mrSges(3,2);
t156 = m(4) * t116 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t103 + mrSges(4,2) * t104 + mrSges(3,3) * t172;
t154 = mrSges(4,1) * t135 - mrSges(4,2) * t132;
t151 = Ifges(4,1) * t132 + t189;
t149 = Ifges(4,2) * t135 + t190;
t148 = Ifges(4,5) * t132 + Ifges(4,6) * t135;
t102 = t135 * t112;
t60 = -pkin(7) * t174 + t102 + (-pkin(3) - t199) * t136;
t66 = -t132 * t198 + t80;
t28 = -t131 * t66 + t134 * t60;
t29 = t131 * t60 + t134 * t66;
t147 = -t132 * t32 + t135 * t31;
t41 = t132 * t110 + t112 * t167 + (-t133 * t170 - t136 * t168) * pkin(5);
t139 = t61 * mrSges(4,1) + t9 * mrSges(5,1) + t120 * Ifges(5,3) + t58 * Ifges(5,5) + t158 * Ifges(5,6) + t125 * Ifges(4,3) + t104 * Ifges(4,5) + t103 * Ifges(4,6) + t160 - (Ifges(3,2) * t136 + t192) * qJD(1) / 0.2e1 - t10 * mrSges(5,2) - t62 * mrSges(4,2);
t127 = -pkin(3) * t135 - pkin(2);
t124 = Ifges(4,3) * t159;
t114 = mrSges(3,3) * t166 - t177;
t111 = (pkin(5) + t200) * t133;
t99 = t166 * t200 + t130;
t79 = -pkin(5) * t175 + t102;
t78 = pkin(3) * t140 + pkin(5) * t169;
t76 = mrSges(4,1) * t125 - mrSges(4,3) * t104;
t75 = -mrSges(4,2) * t125 + mrSges(4,3) * t103;
t74 = -pkin(5) * t162 + t92;
t54 = -t72 * pkin(3) + qJD(2) * t130;
t53 = -mrSges(4,2) * t159 + mrSges(4,3) * t72;
t52 = mrSges(4,1) * t159 - mrSges(4,3) * t71;
t46 = mrSges(5,1) * t120 - mrSges(5,3) * t58;
t45 = -mrSges(5,2) * t120 + mrSges(5,3) * t158;
t42 = -qJD(3) * t80 + t180;
t40 = -mrSges(4,1) * t72 + mrSges(4,2) * t71;
t36 = t71 * Ifges(4,1) + t72 * Ifges(4,4) + Ifges(4,5) * t159;
t35 = t71 * Ifges(4,4) + t72 * Ifges(4,2) + Ifges(4,6) * t159;
t34 = -qJD(2) * t143 + t227 * t88;
t33 = -qJD(2) * t142 - t133 * t64;
t30 = -pkin(7) * t140 + t41;
t27 = t144 * qJD(2) + (-t126 + (-t112 + t198) * t132) * qJD(3) + t180;
t14 = t134 * t43 - t182;
t13 = -t131 * t43 - t181;
t12 = -mrSges(5,2) * t159 + mrSges(5,3) * t19;
t11 = mrSges(5,1) * t159 - mrSges(5,3) * t18;
t8 = -mrSges(5,1) * t19 + mrSges(5,2) * t18;
t7 = t18 * Ifges(5,1) + t19 * Ifges(5,4) + Ifges(5,5) * t159;
t6 = t18 * Ifges(5,4) + t19 * Ifges(5,2) + Ifges(5,6) * t159;
t5 = -qJD(4) * t29 - t131 * t30 + t134 * t27;
t4 = qJD(4) * t28 + t131 * t27 + t134 * t30;
t1 = [(pkin(5) * t40 + t35 * t203 + t150 * t217 + t152 * t218 + t36 * t202 + (-t132 * t31 - t135 * t32) * mrSges(4,3) + (-t135 * t48 / 0.2e1 + t49 * t203 + t148 * t204 + t116 * t154 + t149 * t209 + t151 * t208 + (t132 * t61 - t135 * t62) * mrSges(4,3)) * qJD(3) + (-pkin(5) * t114 + t139 + (Ifges(5,5) * t215 + Ifges(5,6) * t216 + (t188 / 0.2e1 - t187 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,4)) * t133 - 0.2e1 * t212 + (-Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - Ifges(4,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(5) + t153) * pkin(5)) * t136) * qJD(1) + t160) * qJD(2)) * t133 + (t161 + (-0.2e1 * t211 + 0.3e1 / 0.2e1 * Ifges(3,4) * t136) * qJD(1) + t156 * pkin(5) + t235) * t169 + (t10 * t34 - t2 * t87 + t3 * t88 - t33 * t9) * mrSges(5,3) + t54 * (mrSges(5,1) * t87 - mrSges(5,2) * t88) + (-Ifges(5,4) * t88 - Ifges(5,2) * t87) * t223 + (-Ifges(5,1) * t88 - Ifges(5,4) * t87) * t224 + t111 * t8 + t41 * t75 + t42 * t76 + t77 * (-mrSges(5,1) * t34 + mrSges(5,2) * t33) + t78 * t26 + t79 * t52 + t80 * t53 + t4 * t45 + t5 * t46 + t28 * t11 + t29 * t12 + t33 * t23 / 0.2e1 + (-t123 / 0.2e1 - t124 / 0.2e1 + t225) * t136 + t34 * t233 + m(5) * (t10 * t4 + t111 * t54 + t2 * t29 + t28 * t3 + t5 * t9 + t77 * t78) + (Ifges(5,5) * t33 + Ifges(5,6) * t34) * t205 + t7 * t215 + t6 * t216 + (Ifges(5,1) * t33 + Ifges(5,4) * t34) * t219 + (Ifges(5,4) * t33 + Ifges(5,2) * t34) * t221 + m(4) * (t31 * t80 + t32 * t79 + t62 * t41 + t61 * t42); t147 * mrSges(4,3) + (-mrSges(5,1) * t193 + mrSges(5,2) * t194) * t77 + t237 * t46 + (t10 * t238 + t127 * t54 + t2 * t70 + t237 * t9 + t3 * t69 - t77 * t99) * m(5) + t238 * t45 + t132 * t36 / 0.2e1 + t127 * t8 + (t228 * t200 + t137) * qJD(3) + t106 * t7 / 0.2e1 - t99 * t26 - t74 * t75 - t73 * t76 + t69 * t11 + t70 * t12 - pkin(2) * t40 - m(4) * (t61 * t73 + t62 * t74) + ((-m(4) * t146 - t132 * t75 - t135 * t76) * qJD(3) - t132 * t52 + t135 * t53 + m(4) * t147) * pkin(6) + t54 * (mrSges(5,1) * t145 + mrSges(5,2) * t106) + (Ifges(5,4) * t106 - Ifges(5,2) * t145) * t223 + (Ifges(5,1) * t106 - Ifges(5,4) * t145) * t224 + (t10 * t193 - t106 * t3 - t145 * t2 - t194 * t9) * mrSges(5,3) - t145 * t6 / 0.2e1 + (-Ifges(5,5) * t63 - Ifges(5,6) * t64) * t205 + (-Ifges(5,1) * t63 - Ifges(5,4) * t64) * t219 + (-Ifges(5,4) * t63 - Ifges(5,2) * t64) * t221 + (t81 / 0.2e1 - t64 / 0.2e1) * t22 + (t82 / 0.2e1 - t63 / 0.2e1) * t23 + (-Ifges(5,5) * t82 - Ifges(5,6) * t81) * t206 + (-Ifges(5,1) * t82 - Ifges(5,4) * t81) * t220 + (-Ifges(5,4) * t82 - Ifges(5,2) * t81) * t222 + (((t212 + t192 / 0.2e1) * qJD(1) + (t114 + t177) * pkin(5) - t139 + t160 + (Ifges(5,5) * t106 - Ifges(5,6) * t145 + t148) * t242) * t133 + (-t128 / 0.2e1 + t161 + (t211 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t133) * qJD(1) + ((-m(4) * pkin(2) - mrSges(3,1) - t154) * qJD(2) - t156) * pkin(5) - t235) * t136) * qJD(1) + t35 * t202 + t149 * t217 + t151 * t218; (-Ifges(4,2) * t104 + t100 + t49) * t209 + (t103 * t61 + t104 * t62) * mrSges(4,3) + t124 + (t134 * t11 + t131 * t12 + m(5) * (t131 * t2 + t134 * t3) - t228 * t104 + (-t131 * t46 + t134 * t45 + m(5) * (t10 * t134 - t131 * t9)) * qJD(4)) * pkin(3) - t116 * (mrSges(4,1) * t104 + mrSges(4,2) * t103) - t61 * t75 + t62 * t76 - t14 * t45 - t13 * t46 - t225 + t58 * t233 - m(5) * (t10 * t14 + t13 * t9) + (Ifges(4,5) * t103 - Ifges(4,6) * t104) * t204 + t48 * t207 + (Ifges(4,1) * t103 - t191) * t208 + t239; t10 * t46 + t22 * t219 - t9 * t45 + t236 + t239;];
tauc = t1(:);
