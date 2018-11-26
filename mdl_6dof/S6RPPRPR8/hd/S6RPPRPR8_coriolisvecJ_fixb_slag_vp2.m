% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2018-11-23 15:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:43:27
% EndTime: 2018-11-23 15:43:31
% DurationCPUTime: 3.54s
% Computational Cost: add. (3633->390), mult. (8207->508), div. (0->0), fcn. (5382->6), ass. (0->178)
t109 = sin(pkin(9));
t187 = sin(qJ(4));
t153 = t187 * t109;
t110 = cos(pkin(9));
t114 = cos(qJ(4));
t158 = t110 * t114;
t89 = -t153 + t158;
t156 = qJD(4) * t114;
t88 = t114 * t109 + t110 * t187;
t207 = t88 * qJD(3);
t111 = -pkin(1) - qJ(3);
t96 = qJD(1) * t111 + qJD(2);
t149 = -pkin(7) * qJD(1) + t96;
t78 = t149 * t110;
t119 = -qJD(1) * t207 + t78 * t156;
t77 = t149 * t109;
t154 = t187 * t77;
t24 = -t119 + qJD(4) * (t154 - qJD(5));
t176 = mrSges(6,1) + mrSges(5,3);
t175 = mrSges(6,2) - mrSges(5,1);
t210 = t89 * qJD(3);
t150 = qJD(4) * t187;
t173 = -pkin(7) + t111;
t91 = t173 * t109;
t92 = t173 * t110;
t43 = t150 * t91 - t92 * t156 + t207;
t209 = Ifges(6,4) - Ifges(5,5);
t208 = Ifges(6,5) - Ifges(5,6);
t98 = qJD(1) * t153;
t84 = qJD(1) * t158 - t98;
t81 = qJD(6) + t84;
t112 = sin(qJ(6));
t113 = cos(qJ(6));
t52 = t114 * t77 + t187 * t78;
t28 = t210 * qJD(1) + qJD(4) * t52;
t85 = t88 * qJD(4);
t75 = qJD(1) * t85;
t15 = -t75 * pkin(5) + t28;
t107 = qJD(1) * qJD(2);
t124 = qJ(5) * t75 - qJD(5) * t84 + t107;
t191 = pkin(4) + pkin(8);
t152 = t110 * t156;
t76 = qJD(1) * t152 - qJD(4) * t98;
t16 = t191 * t76 + t124;
t51 = -t114 * t78 + t154;
t32 = -t84 * pkin(5) - t51;
t205 = -t32 + qJD(5);
t23 = -qJD(4) * t191 + t205;
t108 = qJD(1) * qJ(2);
t103 = qJD(3) + t108;
t104 = t109 * pkin(3);
t93 = qJD(1) * t104 + t103;
t125 = -qJ(5) * t84 + t93;
t83 = t88 * qJD(1);
t30 = t191 * t83 + t125;
t5 = -t112 * t30 + t113 * t23;
t1 = qJD(6) * t5 + t112 * t15 + t113 * t16;
t6 = t112 * t23 + t113 * t30;
t2 = -qJD(6) * t6 - t112 * t16 + t113 * t15;
t206 = t1 * t112 + t113 * t2;
t204 = -t51 - qJD(5);
t157 = t109 ^ 2 + t110 ^ 2;
t148 = qJD(1) * t157;
t64 = -qJD(4) * t112 + t113 * t83;
t40 = qJD(6) * t64 + t112 * t76;
t65 = qJD(4) * t113 + t112 * t83;
t41 = -qJD(6) * t65 + t113 * t76;
t202 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t40 + Ifges(7,6) * t41;
t167 = Ifges(7,4) * t112;
t131 = Ifges(7,2) * t113 + t167;
t166 = Ifges(7,4) * t113;
t133 = Ifges(7,1) * t112 + t166;
t136 = mrSges(7,1) * t113 - mrSges(7,2) * t112;
t142 = t112 * t5 - t113 * t6;
t162 = Ifges(7,6) * t113;
t165 = Ifges(7,5) * t112;
t185 = Ifges(7,4) * t65;
t18 = t64 * Ifges(7,2) + t81 * Ifges(7,6) + t185;
t189 = -t112 / 0.2e1;
t63 = Ifges(7,4) * t64;
t19 = t65 * Ifges(7,1) + t81 * Ifges(7,5) + t63;
t194 = -t81 / 0.2e1;
t196 = -t65 / 0.2e1;
t198 = -t64 / 0.2e1;
t190 = t83 * pkin(5);
t45 = -qJD(4) * qJ(5) - t52;
t29 = -t45 - t190;
t201 = (t162 + t165) * t194 + t131 * t198 + t133 * t196 + t29 * t136 + t142 * mrSges(7,3) + t19 * t189 - t113 * t18 / 0.2e1;
t200 = t40 / 0.2e1;
t199 = t41 / 0.2e1;
t197 = t64 / 0.2e1;
t195 = t65 / 0.2e1;
t193 = t81 / 0.2e1;
t188 = t113 / 0.2e1;
t186 = m(3) * qJ(2);
t61 = -t114 * t92 + t187 * t91;
t182 = t28 * t61;
t181 = t28 * t89;
t180 = t64 * Ifges(7,6);
t179 = t65 * Ifges(7,5);
t178 = t81 * Ifges(7,3);
t174 = -Ifges(6,6) - Ifges(5,4);
t68 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t83;
t70 = mrSges(6,1) * t83 - qJD(4) * mrSges(6,3);
t172 = t68 - t70;
t171 = -t175 * qJD(4) - t176 * t84;
t35 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t170 = -t70 + t35;
t137 = mrSges(4,1) * t109 + mrSges(4,2) * t110;
t169 = mrSges(5,1) * t83 + mrSges(5,2) * t84 + qJD(1) * t137;
t168 = m(4) * qJD(3);
t164 = Ifges(7,5) * t113;
t163 = Ifges(7,6) * t112;
t161 = qJ(5) * t83;
t86 = -t109 * t150 + t152;
t160 = t112 * t86;
t159 = t113 * t86;
t100 = qJ(2) + t104;
t155 = t68 + t170;
t146 = -qJ(5) * t89 + t100;
t145 = t1 * t113 - t112 * t2;
t143 = t112 * t6 + t113 * t5;
t42 = -qJD(4) * pkin(4) - t204;
t141 = t42 * t85 - t45 * t86;
t140 = t51 * t85 + t52 * t86;
t62 = t114 * t91 + t187 * t92;
t139 = -t61 * t75 - t62 * t76;
t135 = mrSges(7,1) * t112 + mrSges(7,2) * t113;
t134 = Ifges(7,1) * t113 - t167;
t132 = -Ifges(7,2) * t112 + t166;
t20 = -mrSges(7,1) * t75 - mrSges(7,3) * t40;
t21 = mrSges(7,2) * t75 + mrSges(7,3) * t41;
t129 = t112 * t21 + t113 * t20;
t39 = t191 * t88 + t146;
t49 = t89 * pkin(5) + t61;
t12 = t112 * t49 + t113 * t39;
t11 = -t112 * t39 + t113 * t49;
t47 = -mrSges(7,2) * t81 + mrSges(7,3) * t64;
t48 = mrSges(7,1) * t81 - mrSges(7,3) * t65;
t128 = t112 * t48 - t113 * t47;
t127 = -t112 * t47 - t113 * t48;
t123 = qJ(5) * t85 - qJD(5) * t89 + qJD(2);
t122 = -t127 - t171;
t118 = -m(7) * t142 - t128;
t44 = qJD(4) * t62 + t210;
t116 = qJD(1) ^ 2;
t80 = Ifges(5,4) * t83;
t79 = Ifges(6,6) * t83;
t74 = Ifges(7,3) * t75;
t73 = t75 * mrSges(5,2);
t72 = t75 * mrSges(6,3);
t60 = -mrSges(6,2) * t83 - mrSges(6,3) * t84;
t58 = pkin(4) * t84 + t161;
t57 = pkin(4) * t88 + t146;
t56 = t84 * Ifges(5,1) + Ifges(5,5) * qJD(4) - t80;
t55 = t84 * Ifges(5,4) - t83 * Ifges(5,2) + Ifges(5,6) * qJD(4);
t54 = Ifges(6,4) * qJD(4) - t84 * Ifges(6,2) + t79;
t53 = Ifges(6,5) * qJD(4) - t84 * Ifges(6,6) + t83 * Ifges(6,3);
t50 = -t88 * pkin(5) + t62;
t46 = pkin(4) * t83 + t125;
t36 = pkin(4) * t86 + t123;
t34 = t191 * t84 + t161;
t33 = t52 - t190;
t31 = pkin(4) * t76 + t124;
t27 = -t150 * t77 + t119;
t26 = -t85 * pkin(5) + t44;
t25 = -t86 * pkin(5) - t43;
t22 = t191 * t86 + t123;
t17 = t178 + t179 + t180;
t14 = -t76 * pkin(5) - t24;
t13 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t10 = t112 * t33 + t113 * t34;
t9 = -t112 * t34 + t113 * t33;
t8 = t40 * Ifges(7,1) + t41 * Ifges(7,4) - t75 * Ifges(7,5);
t7 = t40 * Ifges(7,4) + t41 * Ifges(7,2) - t75 * Ifges(7,6);
t4 = -qJD(6) * t12 - t112 * t22 + t113 * t26;
t3 = qJD(6) * t11 + t112 * t26 + t113 * t22;
t37 = [(t208 * t86 + t209 * t85) * qJD(4) / 0.2e1 + (-t31 * mrSges(6,3) - t74 / 0.2e1 + mrSges(5,2) * t107 + t174 * t76 + t176 * t28 + (-Ifges(7,3) / 0.2e1 - Ifges(6,2) - Ifges(5,1)) * t75 + t202) * t89 + (((2 * mrSges(3,3)) + t137 + 0.2e1 * t186) * qJD(1) + t169 + m(4) * (t103 + t108) + m(5) * (qJD(1) * t100 + t93)) * qJD(2) + 0.2e1 * qJD(3) * mrSges(4,3) * t148 + (t139 - t140) * mrSges(5,3) + (t139 - t141) * mrSges(6,1) + m(7) * (t1 * t12 + t11 * t2 + t14 * t50 + t25 * t29 + t3 * t6 + t4 * t5) + (-t31 * mrSges(6,2) + t112 * t8 / 0.2e1 + t7 * t188 - t27 * mrSges(5,3) + t24 * mrSges(6,1) - t14 * t136 + t131 * t199 + t133 * t200 + mrSges(5,1) * t107 + (Ifges(6,3) + Ifges(5,2)) * t76 + (-t165 / 0.2e1 - t162 / 0.2e1 - t174) * t75 + t145 * mrSges(7,3) + ((-t163 + t164) * t193 + t132 * t197 + t134 * t195 + t29 * t135 + t18 * t189 + t19 * t188 - t143 * mrSges(7,3)) * qJD(6)) * t88 + m(5) * (t27 * t62 - t43 * t52 + t44 * t51 + t182) + m(6) * (-t24 * t62 + t31 * t57 + t36 * t46 + t42 * t44 + t43 * t45 + t182) - t171 * t44 - t172 * t43 + t19 * t160 / 0.2e1 + t5 * (-mrSges(7,1) * t85 - mrSges(7,3) * t160) + t29 * (-mrSges(7,1) * t159 + mrSges(7,2) * t160) + t18 * t159 / 0.2e1 + t6 * (mrSges(7,2) * t85 + mrSges(7,3) * t159) + t93 * (mrSges(5,1) * t86 - mrSges(5,2) * t85) - t84 * (Ifges(6,2) * t85 + Ifges(6,6) * t86) / 0.2e1 + t83 * (Ifges(6,6) * t85 + Ifges(6,3) * t86) / 0.2e1 + t46 * (-mrSges(6,2) * t86 + mrSges(6,3) * t85) + t86 * t53 / 0.2e1 + t84 * (-Ifges(5,1) * t85 - Ifges(5,4) * t86) / 0.2e1 - t83 * (-Ifges(5,4) * t85 - Ifges(5,2) * t86) / 0.2e1 - t86 * t55 / 0.2e1 + t85 * t54 / 0.2e1 + t57 * (-t76 * mrSges(6,2) + t72) + t36 * t60 + t3 * t47 + t4 * t48 + t50 * t13 + t25 * t35 + t11 * t20 + t12 * t21 - (t56 + t17) * t85 / 0.2e1 + (-t111 * t148 - t157 * t96) * t168 + t100 * (t76 * mrSges(5,1) - t73) + (Ifges(7,5) * t160 + Ifges(7,6) * t159 - Ifges(7,3) * t85) * t193 + (Ifges(7,1) * t160 + Ifges(7,4) * t159 - Ifges(7,5) * t85) * t195 + (Ifges(7,4) * t160 + Ifges(7,2) * t159 - Ifges(7,6) * t85) * t197; (-mrSges(3,3) - t186) * t116 + (-t176 * t76 + t13) * t88 + t155 * t86 + t122 * t85 + (t128 * qJD(6) + t176 * t75 - t129) * t89 + m(7) * (t14 * t88 + t29 * t86 + t143 * t85 + (qJD(6) * t142 - t206) * t89) + m(5) * (t27 * t88 + t140 - t181) + m(6) * (-t24 * t88 + t141 - t181) + (-m(4) * t103 - m(5) * t93 - m(6) * t46 - t157 * t168 - t118 - t169 - t60) * qJD(1); -t112 * t20 + t113 * t21 + t72 - t73 - t175 * t76 + t127 * qJD(6) + (m(4) + m(5)) * t107 - t157 * t116 * mrSges(4,3) + t155 * t83 - t122 * t84 - m(5) * (t51 * t84 - t52 * t83) + m(4) * t96 * t148 + (-t81 * t143 + t29 * t83 + t145) * m(7) + (-t42 * t84 - t45 * t83 + t31) * m(6); (-mrSges(6,1) * qJ(5) + t208) * t76 + (-t164 / 0.2e1 + t163 / 0.2e1 + pkin(4) * mrSges(6,1) + t209) * t75 + (-pkin(4) * t28 - qJ(5) * t24 + t204 * t45 - t42 * t52 - t46 * t58) * m(6) + (-t93 * mrSges(5,1) - t53 / 0.2e1 + t55 / 0.2e1 + t46 * mrSges(6,2) - t45 * mrSges(6,1) + t52 * mrSges(5,3) + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t84 + (Ifges(5,6) / 0.2e1 - Ifges(6,5) / 0.2e1) * qJD(4) + (Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t83 + t201) * t84 - t206 * mrSges(7,3) + t14 * t135 + (-t80 / 0.2e1 - t79 / 0.2e1 + t178 / 0.2e1 + t180 / 0.2e1 + t179 / 0.2e1 - t54 / 0.2e1 + t56 / 0.2e1 + t17 / 0.2e1 + t93 * mrSges(5,2) + t42 * mrSges(6,1) + t51 * mrSges(5,3) - t46 * mrSges(6,3) - t6 * mrSges(7,2) + t5 * mrSges(7,1) + (-Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * qJD(4)) * t83 + t170 * qJD(5) + t171 * t52 + t172 * t51 + t175 * t28 - t58 * t60 - t10 * t47 - t9 * t48 - t32 * t35 - t27 * mrSges(5,2) - t24 * mrSges(6,3) + qJ(5) * t13 + (qJ(5) * t14 - t10 * t6 - t191 * t206 + t205 * t29 - t5 * t9) * m(7) + (-t118 * t191 + t201) * qJD(6) - t129 * t191 + t132 * t199 + t134 * t200 + t8 * t188 + t7 * t189; -t75 * mrSges(6,1) + t84 * t60 - t170 * qJD(4) + (t47 * t81 + t20) * t113 + (-t48 * t81 + t21) * t112 + (-qJD(4) * t29 - t81 * t142 + t206) * m(7) + (qJD(4) * t45 + t46 * t84 + t28) * m(6); -t74 - t29 * (mrSges(7,1) * t65 + mrSges(7,2) * t64) + (Ifges(7,1) * t64 - t185) * t196 + t18 * t195 + (Ifges(7,5) * t64 - Ifges(7,6) * t65) * t194 - t5 * t47 + t6 * t48 + (t5 * t64 + t6 * t65) * mrSges(7,3) + (-Ifges(7,2) * t65 + t19 + t63) * t198 + t202;];
tauc  = t37(:);
