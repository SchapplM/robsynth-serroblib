% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:57
% EndTime: 2019-12-31 19:35:09
% DurationCPUTime: 4.87s
% Computational Cost: add. (2764->389), mult. (7219->513), div. (0->0), fcn. (4737->6), ass. (0->194)
t123 = cos(qJ(5));
t121 = sin(qJ(5));
t141 = mrSges(6,1) * t123 - mrSges(6,2) * t121;
t120 = sin(pkin(8));
t124 = cos(qJ(2));
t168 = cos(pkin(8));
t149 = t168 * t124;
t146 = qJD(1) * t149;
t122 = sin(qJ(2));
t164 = qJD(1) * t122;
t88 = t120 * t164 - t146;
t67 = qJD(2) * t123 + t121 * t88;
t198 = Ifges(6,4) * t67;
t66 = -qJD(2) * t121 + t123 * t88;
t99 = t120 * t124 + t122 * t168;
t90 = t99 * qJD(1);
t85 = qJD(5) + t90;
t18 = Ifges(6,2) * t66 + Ifges(6,6) * t85 + t198;
t65 = Ifges(6,4) * t66;
t19 = Ifges(6,1) * t67 + Ifges(6,5) * t85 + t65;
t200 = -t121 / 0.2e1;
t202 = pkin(4) * t88;
t180 = -qJ(3) - pkin(6);
t107 = t180 * t124;
t102 = qJD(1) * t107;
t150 = t168 * t102;
t106 = t180 * t122;
t101 = qJD(1) * t106;
t97 = qJD(2) * pkin(2) + t101;
t60 = t120 * t97 - t150;
t55 = -qJD(2) * qJ(4) - t60;
t30 = -t55 - t202;
t229 = t30 * t141 - t123 * t18 / 0.2e1 + t19 * t200;
t227 = mrSges(5,2) - mrSges(4,1);
t116 = -pkin(2) * t124 - pkin(1);
t165 = qJD(1) * t116;
t103 = qJD(3) + t165;
t127 = -t90 * qJ(4) + t103;
t44 = t88 * pkin(3) + t127;
t228 = t103 * mrSges(4,1) - t44 * mrSges(5,2) - t229;
t226 = Ifges(5,4) - Ifges(4,5);
t225 = Ifges(5,5) - Ifges(4,6);
t188 = t88 * mrSges(4,3);
t189 = t88 * mrSges(5,1);
t73 = -qJD(2) * mrSges(5,3) + t189;
t179 = -qJD(2) * mrSges(4,2) - t188 - t73;
t186 = t90 * mrSges(4,3);
t187 = t90 * mrSges(5,1);
t224 = -qJD(2) * t227 - t186 - t187;
t32 = -mrSges(6,1) * t66 + mrSges(6,2) * t67;
t223 = -t73 + t32;
t163 = qJD(1) * t124;
t166 = Ifges(3,6) * qJD(2);
t175 = Ifges(3,4) * t122;
t221 = t166 / 0.2e1 + (Ifges(3,2) * t124 + t175) * qJD(1) / 0.2e1 + pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t163);
t159 = qJD(1) * qJD(2);
t154 = t122 * t159;
t112 = pkin(2) * t154;
t81 = qJD(2) * t146 - t120 * t154;
t130 = -qJ(4) * t81 - qJD(4) * t90 + t112;
t203 = pkin(3) + pkin(7);
t89 = t99 * qJD(2);
t80 = qJD(1) * t89;
t14 = t203 * t80 + t130;
t151 = qJD(2) * t180;
t87 = -t122 * qJD(3) + t124 * t151;
t126 = t87 * qJD(1);
t86 = qJD(3) * t124 + t122 * t151;
t75 = t86 * qJD(1);
t36 = t120 * t75 - t168 * t126;
t22 = pkin(4) * t81 + t36;
t24 = t203 * t88 + t127;
t92 = t120 * t102;
t59 = t168 * t97 + t92;
t131 = qJD(4) - t59;
t201 = t90 * pkin(4);
t25 = -qJD(2) * t203 + t131 + t201;
t5 = -t121 * t24 + t123 * t25;
t1 = qJD(5) * t5 + t121 * t22 + t123 * t14;
t6 = t121 * t25 + t123 * t24;
t2 = -qJD(5) * t6 - t121 * t14 + t123 * t22;
t40 = qJD(5) * t66 + t121 * t80;
t41 = -qJD(5) * t67 + t123 * t80;
t220 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t40 + Ifges(6,6) * t41;
t218 = t5 * mrSges(6,1) + t103 * mrSges(4,2) - t6 * mrSges(6,2) - t44 * mrSges(5,3);
t217 = -0.2e1 * pkin(1);
t215 = t40 / 0.2e1;
t214 = t41 / 0.2e1;
t213 = -t66 / 0.2e1;
t212 = t66 / 0.2e1;
t211 = -t67 / 0.2e1;
t210 = t67 / 0.2e1;
t209 = -t85 / 0.2e1;
t208 = t85 / 0.2e1;
t207 = -t88 / 0.2e1;
t205 = -t90 / 0.2e1;
t204 = t90 / 0.2e1;
t199 = t123 / 0.2e1;
t197 = pkin(2) * t120;
t196 = pkin(6) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t164);
t63 = -t168 * t106 - t107 * t120;
t194 = t36 * t63;
t193 = t66 * Ifges(6,6);
t192 = t67 * Ifges(6,5);
t191 = t81 * mrSges(5,1);
t190 = t85 * Ifges(6,3);
t185 = t90 * Ifges(4,4);
t184 = t90 * Ifges(5,6);
t183 = -qJD(2) / 0.2e1;
t181 = -Ifges(4,4) - Ifges(5,6);
t37 = t120 * t126 + t168 * t75;
t177 = mrSges(6,3) * t121;
t176 = mrSges(6,3) * t123;
t174 = Ifges(6,4) * t121;
t173 = Ifges(6,4) * t123;
t172 = Ifges(6,5) * t121;
t171 = Ifges(6,6) * t123;
t167 = Ifges(3,5) * qJD(2);
t162 = qJD(2) * t122;
t161 = qJD(5) * t121;
t160 = qJD(5) * t123;
t118 = pkin(2) * t164;
t158 = -Ifges(4,4) / 0.2e1 - Ifges(5,6) / 0.2e1;
t155 = t168 * pkin(2);
t153 = t124 * t159;
t46 = t120 * t86 - t168 * t87;
t148 = qJ(4) * t88 + t118;
t61 = t101 * t120 - t150;
t115 = -t155 - pkin(3);
t145 = t1 * t123 - t121 * t2;
t144 = t121 * t6 + t123 * t5;
t143 = t121 * t5 - t123 * t6;
t64 = t120 * t106 - t107 * t168;
t142 = t63 * t81 - t64 * t80;
t31 = -qJD(2) * qJD(4) - t37;
t140 = mrSges(6,1) * t121 + mrSges(6,2) * t123;
t139 = Ifges(6,1) * t123 - t174;
t138 = Ifges(6,1) * t121 + t173;
t137 = -Ifges(6,2) * t121 + t173;
t136 = Ifges(6,2) * t123 + t174;
t135 = Ifges(6,5) * t123 - Ifges(6,6) * t121;
t134 = t171 + t172;
t132 = -t99 * qJ(4) + t116;
t98 = t120 * t122 - t149;
t35 = t203 * t98 + t132;
t48 = pkin(4) * t99 + t63;
t13 = t121 * t48 + t123 * t35;
t12 = -t121 * t35 + t123 * t48;
t42 = -mrSges(6,2) * t85 + mrSges(6,3) * t66;
t43 = mrSges(6,1) * t85 - mrSges(6,3) * t67;
t133 = -t121 * t42 - t123 * t43;
t91 = qJD(2) * t149 - t120 * t162;
t129 = pkin(2) * t162 - qJ(4) * t91 - qJD(4) * t99;
t47 = t120 * t87 + t168 * t86;
t62 = t101 * t168 + t92;
t125 = -qJD(5) * t143 + t1 * t121 + t123 * t2;
t117 = Ifges(3,4) * t163;
t113 = qJ(4) + t197;
t96 = Ifges(3,1) * t164 + t117 + t167;
t84 = Ifges(4,4) * t88;
t83 = Ifges(5,6) * t88;
t79 = Ifges(6,3) * t81;
t78 = t81 * mrSges(5,3);
t77 = t81 * mrSges(4,2);
t58 = t98 * pkin(3) + t132;
t57 = -mrSges(5,2) * t88 - mrSges(5,3) * t90;
t56 = mrSges(4,1) * t88 + mrSges(4,2) * t90;
t54 = t90 * Ifges(4,1) + Ifges(4,5) * qJD(2) - t84;
t53 = -t88 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t185;
t52 = Ifges(5,4) * qJD(2) - t90 * Ifges(5,2) + t83;
t51 = Ifges(5,5) * qJD(2) + t88 * Ifges(5,3) - t184;
t50 = -qJD(2) * pkin(3) + t131;
t49 = -t98 * pkin(4) + t64;
t45 = pkin(3) * t90 + t148;
t39 = t62 - t201;
t38 = t61 - t202;
t29 = pkin(3) * t89 + t129;
t28 = -t89 * pkin(4) + t47;
t27 = pkin(4) * t91 + t46;
t26 = t203 * t90 + t148;
t23 = pkin(3) * t80 + t130;
t21 = -mrSges(6,2) * t81 + mrSges(6,3) * t41;
t20 = mrSges(6,1) * t81 - mrSges(6,3) * t40;
t17 = t190 + t192 + t193;
t16 = -pkin(4) * t80 - t31;
t15 = t203 * t89 + t129;
t11 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t10 = t121 * t38 + t123 * t26;
t9 = -t121 * t26 + t123 * t38;
t8 = t40 * Ifges(6,1) + t41 * Ifges(6,4) + t81 * Ifges(6,5);
t7 = t40 * Ifges(6,4) + t41 * Ifges(6,2) + t81 * Ifges(6,6);
t4 = -qJD(5) * t13 - t121 * t15 + t123 * t27;
t3 = qJD(5) * t12 + t121 * t27 + t123 * t15;
t33 = [t116 * (t80 * mrSges(4,1) + t77) + t58 * (-t80 * mrSges(5,2) - t78) + t29 * t57 + t3 * t42 + t4 * t43 + t49 * t11 + t28 * t32 + t12 * t20 + t13 * t21 + (-t52 / 0.2e1 + t54 / 0.2e1 + t17 / 0.2e1 + t190 / 0.2e1 + t192 / 0.2e1 + t193 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t90 + t158 * t88 + t218) * t91 + (t51 / 0.2e1 - t53 / 0.2e1 + t134 * t208 + t138 * t210 + t136 * t212 + t158 * t90 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t88 - t143 * mrSges(6,3) + t228) * t89 + t179 * t47 - t224 * t46 + m(4) * (t37 * t64 - t46 * t59 + t47 * t60 + t194) + m(6) * (t1 * t13 + t12 * t2 + t16 * t49 + t28 * t30 + t3 * t6 + t4 * t5) + m(5) * (t23 * t58 + t29 * t44 - t31 * t64 + t46 * t50 - t47 * t55 + t194) + (-t59 * t91 - t60 * t89 + t142) * mrSges(4,3) + (t50 * t91 + t55 * t89 + t142) * mrSges(5,1) + ((-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t91 + (Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1) * t89 + (-t196 + t96 / 0.2e1 + t167 / 0.2e1 + (mrSges(3,2) * t217 + 0.3e1 / 0.2e1 * Ifges(3,4) * t124) * qJD(1)) * t124) * qJD(2) + (t79 / 0.2e1 - t23 * mrSges(5,3) + t181 * t80 + (mrSges(5,1) + mrSges(4,3)) * t36 + (Ifges(6,3) / 0.2e1 + Ifges(4,1) + Ifges(5,2)) * t81 + t220) * t99 + (t138 * t215 + t136 * t214 - t16 * t141 - t23 * mrSges(5,2) + t7 * t199 + t121 * t8 / 0.2e1 + t31 * mrSges(5,1) - t37 * mrSges(4,3) + (t172 / 0.2e1 + t171 / 0.2e1 + t181) * t81 + (Ifges(5,3) + Ifges(4,2)) * t80 + t145 * mrSges(6,3) + (-mrSges(6,3) * t144 + t135 * t208 + t137 * t212 + t139 * t210 + t140 * t30 + t18 * t200 + t19 * t199) * qJD(5)) * t98 + (-t166 / 0.2e1 + (mrSges(3,1) * t217 - 0.3e1 / 0.2e1 * t175 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t124) * qJD(1) + (m(4) * (t103 + t165) + qJD(1) * (mrSges(4,1) * t98 + mrSges(4,2) * t99) + t56) * pkin(2) - t221) * t162; (-t160 * t6 + t161 * t5) * mrSges(6,3) + t113 * t11 - (Ifges(3,5) * t124 - Ifges(3,6) * t122) * t159 / 0.2e1 - t56 * t118 + t137 * t214 + t139 * t215 + t50 * t189 + t115 * t191 + t163 * t196 + t8 * t199 + t7 * t200 + t60 * t186 - t45 * t57 - t39 * t32 - t10 * t42 - t9 * t43 - t31 * mrSges(5,3) - t37 * mrSges(4,2) + (t51 - t185) * t205 + (-t10 * t6 + t113 * t16 - t5 * t9 + (-t39 + qJD(4)) * t30) * m(6) - Ifges(3,6) * t154 - t55 * t187 - t59 * t188 - t2 * t176 - t1 * t177 + (t52 + t83) * t207 + t223 * qJD(4) + t224 * t61 - t179 * t62 + (-mrSges(5,1) * t113 - mrSges(4,3) * t197 + t225) * t80 + (-t113 * t31 + t115 * t36 - t44 * t45 - t50 * t61 + (t62 - qJD(4)) * t55) * m(5) + (t53 + t184) * t204 + (-t103 * t118 + t59 * t61 - t60 * t62 + (t120 * t37 - t168 * t36) * pkin(2)) * m(4) + t221 * t164 + t229 * qJD(5) + (-mrSges(4,3) * t155 + t135 / 0.2e1 - t226) * t81 + (-Ifges(4,1) * t205 - Ifges(6,5) * t211 + Ifges(5,2) * t204 - Ifges(6,6) * t213 - Ifges(6,3) * t209 + t183 * t226 + t218) * t88 + t227 * t36 + t16 * t140 + Ifges(3,5) * t153 + (Ifges(5,3) * t207 + t134 * t209 + t136 * t213 + t138 * t211 - t176 * t6 + t177 * t5 + t183 * t225 - t228) * t90 + (m(6) * t125 + t121 * t21 + t123 * t20 + t160 * t42 - t161 * t43) * (-pkin(7) + t115) - (-Ifges(3,2) * t164 + t117 + t96) * t163 / 0.2e1 - (t134 * t85 + t136 * t66 + t138 * t67) * qJD(5) / 0.2e1 + (-Ifges(4,2) * t90 + t17 + t54 - t84) * t88 / 0.2e1 + (pkin(1) * (mrSges(3,1) * t122 + mrSges(3,2) * t124) - t122 * (Ifges(3,1) * t124 - t175) / 0.2e1) * qJD(1) ^ 2 + (-mrSges(3,1) * t153 + mrSges(3,2) * t154) * pkin(6); -t121 * t20 + t123 * t21 + t77 - t78 - t227 * t80 + t133 * qJD(5) + (t32 + t179) * t88 + (t133 + t224) * t90 + (-t144 * t85 + t30 * t88 + t145) * m(6) + (-t50 * t90 - t55 * t88 + t23) * m(5) + (t59 * t90 + t60 * t88 + t112) * m(4); t191 + t90 * t57 - t223 * qJD(2) + (t42 * t85 + t20) * t123 + (-t43 * t85 + t21) * t121 + (-qJD(2) * t30 - t143 * t90 + t125) * m(6) + (qJD(2) * t55 + t44 * t90 + t36) * m(5); t79 - t30 * (mrSges(6,1) * t67 + mrSges(6,2) * t66) + (Ifges(6,1) * t66 - t198) * t211 + t18 * t210 + (Ifges(6,5) * t66 - Ifges(6,6) * t67) * t209 - t5 * t42 + t6 * t43 + (t5 * t66 + t6 * t67) * mrSges(6,3) + (-Ifges(6,2) * t67 + t19 + t65) * t213 + t220;];
tauc = t33(:);
