% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:02
% EndTime: 2019-12-31 19:54:11
% DurationCPUTime: 4.09s
% Computational Cost: add. (3576->348), mult. (9633->455), div. (0->0), fcn. (6811->6), ass. (0->162)
t213 = Ifges(5,1) + Ifges(6,1);
t212 = Ifges(6,4) + Ifges(5,5);
t215 = -mrSges(5,1) - mrSges(6,1);
t211 = -Ifges(5,6) + Ifges(6,6);
t147 = sin(qJ(2));
t148 = cos(qJ(2));
t178 = -qJ(3) - pkin(6);
t156 = qJD(2) * t178;
t106 = qJD(3) * t148 + t147 * t156;
t107 = -t147 * qJD(3) + t148 * t156;
t144 = sin(pkin(8));
t145 = cos(pkin(8));
t57 = t145 * t106 + t144 * t107;
t167 = t145 * t148;
t120 = -t144 * t147 + t167;
t109 = t120 * qJD(1);
t163 = qJD(1) * t148;
t164 = qJD(1) * t147;
t110 = -t144 * t163 - t145 * t164;
t146 = sin(qJ(4));
t188 = cos(qJ(4));
t152 = t146 * t109 - t110 * t188;
t66 = t188 * t109 + t146 * t110;
t33 = pkin(4) * t152 - qJ(5) * t66;
t179 = -Ifges(5,4) + Ifges(6,5);
t143 = qJD(2) + qJD(4);
t184 = Ifges(6,5) * t66;
t61 = Ifges(5,4) * t66;
t210 = t212 * t143 + t213 * t152 - t184 + t61;
t138 = pkin(2) * t145 + pkin(3);
t183 = pkin(2) * t144;
t105 = t146 * t138 + t188 * t183;
t182 = pkin(7) * t109;
t131 = t178 * t147;
t125 = qJD(1) * t131;
t132 = t178 * t148;
t126 = qJD(1) * t132;
t168 = t145 * t126;
t73 = -t125 * t144 + t168;
t154 = t73 - t182;
t181 = pkin(7) * t110;
t113 = t144 * t126;
t74 = t145 * t125 + t113;
t50 = t74 + t181;
t209 = t105 * qJD(4) - t146 * t50 + t154 * t188;
t119 = qJD(2) * pkin(2) + t125;
t70 = t144 * t119 - t168;
t47 = t70 + t182;
t173 = t146 * t47;
t69 = t145 * t119 + t113;
t45 = qJD(2) * pkin(3) + t181 + t69;
t12 = t188 * t45 - t173;
t208 = -t12 + qJD(5);
t207 = t211 * t152;
t10 = -t143 * pkin(4) + t208;
t139 = -pkin(2) * t148 - pkin(1);
t165 = qJD(1) * t139;
t128 = qJD(3) + t165;
t75 = -t109 * pkin(3) + t128;
t14 = -pkin(4) * t66 - qJ(5) * t152 + t75;
t206 = t75 * mrSges(5,2) + t10 * mrSges(6,2) - t12 * mrSges(5,3) - t14 * mrSges(6,3);
t13 = t146 * t45 + t188 * t47;
t11 = t143 * qJ(5) + t13;
t60 = Ifges(6,5) * t152;
t21 = Ifges(6,6) * t143 - Ifges(6,3) * t66 + t60;
t185 = Ifges(5,4) * t152;
t22 = Ifges(5,2) * t66 + Ifges(5,6) * t143 + t185;
t205 = t11 * mrSges(6,2) + t13 * mrSges(5,3) - t21 / 0.2e1 + t22 / 0.2e1 - t14 * mrSges(6,1) - t75 * mrSges(5,1);
t203 = t66 / 0.2e1;
t202 = -t66 / 0.2e1;
t200 = -t152 / 0.2e1;
t199 = t152 / 0.2e1;
t198 = pkin(1) * mrSges(3,1);
t197 = pkin(1) * mrSges(3,2);
t121 = t144 * t148 + t145 * t147;
t77 = t145 * t131 + t132 * t144;
t58 = -pkin(7) * t121 + t77;
t122 = t144 * t131;
t78 = -t145 * t132 + t122;
t59 = pkin(7) * t120 + t78;
t153 = -t146 * t59 + t188 * t58;
t48 = (-t121 * qJD(3) + (t167 * t178 - t122) * qJD(2)) * qJD(1);
t112 = t120 * qJD(2);
t98 = qJD(1) * t112;
t149 = -t98 * pkin(7) + t48;
t49 = t57 * qJD(1);
t111 = t121 * qJD(2);
t97 = qJD(1) * t111;
t41 = -pkin(7) * t97 + t49;
t5 = qJD(4) * t13 + t146 * t41 - t149 * t188;
t196 = t153 * t5;
t72 = t146 * t120 + t121 * t188;
t195 = t5 * t72;
t193 = -t110 / 0.2e1;
t192 = -t111 / 0.2e1;
t191 = t112 / 0.2e1;
t190 = -t143 / 0.2e1;
t187 = mrSges(5,3) * t66;
t186 = mrSges(5,3) * t152;
t180 = mrSges(6,2) + mrSges(5,3);
t177 = mrSges(6,2) * t152 + t143 * t215 + t186;
t51 = -mrSges(5,2) * t143 + t187;
t54 = mrSges(6,2) * t66 + mrSges(6,3) * t143;
t176 = t54 + t51;
t175 = Ifges(3,4) * t147;
t174 = Ifges(4,4) * t110;
t172 = Ifges(3,5) * qJD(2);
t171 = Ifges(3,6) * qJD(2);
t170 = qJD(2) * mrSges(3,1);
t169 = qJD(2) * mrSges(3,2);
t162 = qJD(2) * t147;
t161 = t146 * t183;
t141 = pkin(2) * t164;
t160 = t172 / 0.2e1;
t159 = -t171 / 0.2e1;
t158 = t97 * mrSges(4,1) + t98 * mrSges(4,2);
t137 = qJD(2) * t141;
t76 = pkin(3) * t97 + t137;
t157 = qJD(4) * t188;
t81 = -pkin(3) * t110 + t141;
t82 = pkin(2) * t162 + pkin(3) * t111;
t56 = -t106 * t144 + t145 * t107;
t87 = -t120 * pkin(3) + t139;
t92 = -qJD(4) * t161 + t138 * t157;
t155 = -pkin(7) * t112 + t56;
t20 = t146 * t58 + t188 * t59;
t4 = -qJD(4) * t173 + t146 * t149 + t45 * t157 + t188 * t41;
t151 = t120 * t188 - t146 * t121;
t104 = t138 * t188 - t161;
t2 = qJD(5) * t143 + t4;
t31 = qJD(4) * t66 - t146 * t97 + t188 * t98;
t32 = qJD(4) * t152 + t146 * t98 + t188 * t97;
t150 = -t4 * mrSges(5,2) + t2 * mrSges(6,3) + t211 * t32 + t212 * t31 + t215 * t5;
t140 = Ifges(3,4) * t163;
t130 = mrSges(3,3) * t163 - t169;
t129 = -mrSges(3,3) * t164 + t170;
t118 = Ifges(3,1) * t164 + t140 + t172;
t117 = t171 + (t148 * Ifges(3,2) + t175) * qJD(1);
t103 = Ifges(4,4) * t109;
t102 = -pkin(4) - t104;
t101 = qJ(5) + t105;
t86 = qJD(5) + t92;
t85 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t110;
t84 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t109;
t68 = -mrSges(4,1) * t109 - mrSges(4,2) * t110;
t63 = -Ifges(4,1) * t110 + Ifges(4,5) * qJD(2) + t103;
t62 = Ifges(4,2) * t109 + Ifges(4,6) * qJD(2) - t174;
t44 = -pkin(7) * t111 + t57;
t37 = qJD(4) * t72 + t111 * t188 + t146 * t112;
t36 = qJD(4) * t151 - t146 * t111 + t112 * t188;
t35 = -mrSges(5,1) * t66 + mrSges(5,2) * t152;
t34 = -mrSges(6,1) * t66 - mrSges(6,3) * t152;
t26 = t32 * mrSges(6,1);
t25 = t31 * mrSges(5,2);
t18 = -pkin(4) * t151 - t72 * qJ(5) + t87;
t17 = t33 + t81;
t16 = t146 * t154 + t188 * t50;
t9 = pkin(4) * t37 - qJ(5) * t36 - qJD(5) * t72 + t82;
t8 = qJD(4) * t20 + t146 * t44 - t155 * t188;
t7 = qJD(4) * t153 + t146 * t155 + t188 * t44;
t6 = pkin(4) * t32 - qJ(5) * t31 - qJD(5) * t152 + t76;
t1 = [t128 * (mrSges(4,1) * t111 + mrSges(4,2) * t112) + m(5) * (-t12 * t8 + t13 * t7 + t20 * t4 + t75 * t82 + t76 * t87 - t196) + m(6) * (t10 * t8 + t11 * t7 + t14 * t9 + t18 * t6 + t2 * t20 - t196) + (t151 * t4 + t195) * mrSges(5,3) + (t151 * t2 + t195) * mrSges(6,2) + t6 * (-mrSges(6,1) * t151 - mrSges(6,3) * t72) + t76 * (-mrSges(5,1) * t151 + mrSges(5,2) * t72) + (mrSges(5,1) * t87 + t179 * t72 - (Ifges(5,2) + Ifges(6,3)) * t151 - t180 * t20) * t32 + (-mrSges(6,3) * t18 - t151 * t179 - t153 * t180 + t213 * t72) * t31 + (t109 * t191 - t111 * t193 + t120 * t98 - t121 * t97) * Ifges(4,4) + (-t111 * t70 - t112 * t69 + t120 * t49 - t121 * t48 - t77 * t98 - t78 * t97) * mrSges(4,3) + (t109 * t192 - t120 * t97) * Ifges(4,2) + t210 * t36 / 0.2e1 + (t211 * t37 + t212 * t36) * t143 / 0.2e1 + (t179 * t37 + t213 * t36) * t199 + t82 * t35 + t57 * t84 + t56 * t85 + t9 * t34 + t87 * t25 + (Ifges(4,5) * t191 + Ifges(4,6) * t192 + (t118 / 0.2e1 - pkin(6) * t129 + t160 + (-0.2e1 * t197 + 0.3e1 / 0.2e1 * Ifges(3,4) * t148) * qJD(1)) * t148) * qJD(2) + (-t117 / 0.2e1 - pkin(6) * t130 + t159 + (-0.2e1 * t198 - 0.3e1 / 0.2e1 * t175 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t148) * qJD(1) + (m(4) * (t128 + t165) + qJD(1) * (-mrSges(4,1) * t120 + mrSges(4,2) * t121) + t68) * pkin(2)) * t162 + (t112 * t193 + t121 * t98) * Ifges(4,1) + (-Ifges(5,2) * t203 + Ifges(6,3) * t202 - t205) * t37 + (Ifges(5,4) * t203 + Ifges(6,5) * t202 + t206) * t36 + t18 * t26 + m(4) * (t48 * t77 + t49 * t78 + t56 * t69 + t57 * t70) + t176 * t7 + t177 * t8 + t139 * t158 + t63 * t191 + t62 * t192; -(Ifges(4,2) * t110 + t103 + t63) * t109 / 0.2e1 + (-t101 * t32 + t102 * t31) * mrSges(6,2) + (-t104 * t31 - t105 * t32) * mrSges(5,3) + (t69 * t109 - t70 * t110 + (-t144 * t97 - t145 * t98) * pkin(2)) * mrSges(4,3) + m(4) * (t144 * t49 + t145 * t48) * pkin(2) + (-Ifges(5,2) * t202 + Ifges(6,3) * t203 + t179 * t200 + t205) * t152 + (-t104 * t5 + t105 * t4 - t75 * t81 + (-t16 + t92) * t13 - t209 * t12) * m(5) + (t101 * t2 + t102 * t5 - t14 * t17 + (-t16 + t86) * t11 + t209 * t10) * m(6) + t209 * t177 - t81 * t35 - t74 * t84 - t73 * t85 + t48 * mrSges(4,1) - t49 * mrSges(4,2) - t17 * t34 + ((-t140 / 0.2e1 - t118 / 0.2e1 + t160 + qJD(1) * t197 + (t129 - t170) * pkin(6)) * t148 + (t117 / 0.2e1 + t159 + (t198 + t175 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t148) * qJD(1) + (t130 + t169) * pkin(6) + (-m(4) * t128 - t68) * pkin(2)) * t147) * qJD(1) + (-t210 / 0.2e1 + t212 * t190 + t213 * t200 + Ifges(5,4) * t202 + Ifges(6,5) * t203 - t206) * t66 - m(4) * (t69 * t73 + t70 * t74) - t176 * t16 + t110 * (Ifges(4,1) * t109 + t174) / 0.2e1 + t150 + t86 * t54 + t92 * t51 - Ifges(4,6) * t97 + Ifges(4,5) * t98 - qJD(2) * (Ifges(4,5) * t109 + Ifges(4,6) * t110) / 0.2e1 - t128 * (-t110 * mrSges(4,1) + t109 * mrSges(4,2)) + t207 * t190 + t62 * t193; t32 * mrSges(5,1) - t31 * mrSges(6,3) - t109 * t84 - t110 * t85 - t176 * t66 - t177 * t152 + t158 + t25 + t26 + (-t10 * t152 - t11 * t66 + t6) * m(6) + (t12 * t152 - t13 * t66 + t76) * m(5) + (-t109 * t70 - t110 * t69 + t137) * m(4); -t75 * (mrSges(5,1) * t152 + mrSges(5,2) * t66) + t22 * t199 + (Ifges(6,3) * t152 + t184) * t203 - t14 * (mrSges(6,1) * t152 - mrSges(6,3) * t66) + qJD(5) * t54 - t33 * t34 + t150 + (-t177 + t186) * t13 + (-t176 + t187) * t12 + (-pkin(4) * t31 - qJ(5) * t32 - t10 * t66 + t11 * t152) * mrSges(6,2) + (t212 * t66 + t207) * t190 + (-pkin(4) * t5 + qJ(5) * t2 - t10 * t13 + t11 * t208 - t14 * t33) * m(6) + (-Ifges(5,2) * t152 + t210 + t61) * t202 + (t213 * t66 - t185 + t21 + t60) * t200; t31 * mrSges(6,2) - t143 * t54 + t152 * t34 + 0.2e1 * (t5 / 0.2e1 + t11 * t190 + t14 * t199) * m(6);];
tauc = t1(:);
