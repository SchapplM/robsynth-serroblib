% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:35
% EndTime: 2019-03-08 19:34:43
% DurationCPUTime: 3.62s
% Computational Cost: add. (2610->369), mult. (6579->504), div. (0->0), fcn. (4396->10), ass. (0->189)
t116 = sin(qJ(4));
t119 = cos(qJ(4));
t135 = pkin(9) * t116 - qJ(5) * t119;
t127 = t135 * qJD(4);
t164 = qJD(5) * t116;
t167 = qJD(4) * t116;
t147 = pkin(4) * t167 - t164;
t111 = sin(pkin(11));
t112 = sin(pkin(6));
t113 = cos(pkin(11));
t117 = sin(qJ(2));
t120 = cos(qJ(2));
t67 = (t111 * t120 + t113 * t117) * t112;
t62 = qJD(1) * t67;
t234 = t127 + t147 - t62;
t115 = sin(qJ(6));
t118 = cos(qJ(6));
t142 = mrSges(7,1) * t118 - mrSges(7,2) * t115;
t121 = -pkin(4) - pkin(9);
t170 = qJD(1) * t112;
t154 = t117 * t170;
t153 = t120 * t170;
t90 = qJD(2) * pkin(2) + t153;
t51 = t111 * t90 + t113 * t154;
t48 = qJD(2) * pkin(8) + t51;
t148 = pkin(5) * qJD(2) + t48;
t143 = t148 * t116;
t114 = cos(pkin(6));
t98 = qJD(1) * t114 + qJD(3);
t92 = t119 * t98;
t28 = t92 - t143;
t226 = qJD(5) - t28;
t20 = qJD(4) * t121 + t226;
t149 = -qJ(5) * t116 - pkin(3);
t126 = t119 * t121 + t149;
t91 = t111 * t154;
t50 = t113 * t90 - t91;
t32 = qJD(2) * t126 - t50;
t5 = -t115 * t32 + t118 * t20;
t6 = t115 * t20 + t118 * t32;
t144 = t115 * t5 - t118 * t6;
t199 = -t118 / 0.2e1;
t200 = -t115 / 0.2e1;
t168 = qJD(2) * t119;
t107 = pkin(5) * t168;
t175 = t116 * t98;
t35 = t119 * t48 + t175;
t31 = -qJD(4) * qJ(5) - t35;
t25 = t107 - t31;
t169 = qJD(2) * t116;
t101 = qJD(6) + t169;
t166 = qJD(4) * t118;
t85 = -t115 * t168 + t166;
t198 = Ifges(7,4) * t85;
t84 = -qJD(4) * t115 - t118 * t168;
t40 = Ifges(7,2) * t84 + Ifges(7,6) * t101 + t198;
t81 = Ifges(7,4) * t84;
t41 = Ifges(7,1) * t85 + Ifges(7,5) * t101 + t81;
t233 = -t144 * mrSges(7,3) - t25 * t142 - t199 * t40 - t200 * t41;
t95 = -mrSges(6,1) * t168 - qJD(4) * mrSges(6,3);
t186 = qJD(4) * mrSges(5,2) - mrSges(5,3) * t168 + t95;
t44 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t160 = -t44 + t186;
t216 = -m(5) * t35 + m(6) * t31;
t232 = -m(7) * t25 + t160 + t216;
t231 = -t168 / 0.2e1;
t172 = t116 * t118;
t197 = pkin(2) * t113;
t68 = t126 - t197;
t103 = pkin(2) * t111 + pkin(8);
t190 = pkin(5) + t103;
t82 = t190 * t116;
t38 = t115 * t82 + t118 * t68;
t65 = t113 * t153 - t91;
t83 = t190 * t119;
t71 = qJD(4) * t83;
t230 = -qJD(6) * t38 - t115 * t234 + t118 * t71 - t172 * t65;
t173 = t115 * t116;
t37 = -t115 * t68 + t118 * t82;
t229 = qJD(6) * t37 + t115 * t71 + t118 * t234 - t173 * t65;
t66 = (t111 * t117 - t113 * t120) * t112;
t64 = qJD(2) * t66;
t55 = qJD(1) * t64;
t176 = t116 * t55;
t10 = -t176 + (t119 * t148 + t175) * qJD(4);
t106 = pkin(4) * t169;
t100 = qJD(4) * t106;
t125 = t62 - t164;
t27 = t100 + (t127 + t125) * qJD(2);
t1 = qJD(6) * t5 + t10 * t115 + t118 * t27;
t2 = -qJD(6) * t6 + t10 * t118 - t115 * t27;
t217 = t1 * t115 + t118 * t2;
t228 = m(7) * t217;
t219 = mrSges(6,1) + mrSges(5,3);
t227 = mrSges(6,2) - mrSges(5,1);
t179 = Ifges(7,6) * t118;
t182 = Ifges(7,5) * t115;
t136 = t179 + t182;
t184 = Ifges(7,4) * t115;
t137 = Ifges(7,2) * t118 + t184;
t183 = Ifges(7,4) * t118;
t139 = Ifges(7,1) * t115 + t183;
t201 = t101 / 0.2e1;
t205 = t85 / 0.2e1;
t220 = qJD(4) / 0.2e1;
t221 = -qJD(4) / 0.2e1;
t222 = qJD(2) / 0.2e1;
t130 = -pkin(4) * t119 + t149;
t36 = qJD(2) * t130 - t50;
t47 = -qJD(2) * pkin(3) - t50;
t225 = t31 * mrSges(6,1) + t47 * mrSges(5,1) + Ifges(5,6) * t221 - (Ifges(5,4) * t116 + t119 * Ifges(5,2)) * qJD(2) / 0.2e1 + Ifges(6,5) * t220 + (-Ifges(6,6) * t116 - t119 * Ifges(6,3)) * t222 + t136 * t201 - t35 * mrSges(5,3) - t36 * mrSges(6,2) + t84 * t137 / 0.2e1 + t139 * t205 + t233;
t224 = -Ifges(5,1) / 0.2e1;
t223 = Ifges(5,4) * t231;
t177 = t116 * t48;
t34 = -t92 + t177;
t218 = -qJD(5) - t34;
t163 = qJD(6) * t119;
t151 = t118 * t163;
t161 = qJD(4) * qJD(6);
t56 = -t115 * t161 + (t115 * t167 - t151) * qJD(2);
t152 = t115 * t163;
t57 = -t118 * t161 + (t116 * t166 + t152) * qJD(2);
t215 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t56 + Ifges(7,6) * t57;
t59 = -mrSges(7,2) * t101 + mrSges(7,3) * t84;
t60 = mrSges(7,1) * t101 - mrSges(7,3) * t85;
t133 = t115 * t60 - t118 * t59;
t162 = qJD(4) * qJD(2);
t150 = t119 * t162;
t42 = mrSges(7,1) * t150 - mrSges(7,3) * t56;
t43 = -mrSges(7,2) * t150 + mrSges(7,3) * t57;
t134 = t115 * t43 + t118 * t42;
t212 = t133 * qJD(6) - t134;
t165 = qJD(4) * t119;
t189 = -t119 * t55 + t98 * t165;
t11 = (-qJD(5) + t177) * qJD(4) - t189;
t12 = -t167 * t48 + t189;
t211 = m(5) * t12 - m(6) * t11;
t185 = -t227 * qJD(4) - t219 * t169;
t30 = -qJD(4) * pkin(4) - t218;
t210 = m(5) * t34 + m(6) * t30 - t185;
t209 = 0.2e1 * m(5);
t208 = -t65 / 0.2e1;
t207 = -t84 / 0.2e1;
t206 = -t85 / 0.2e1;
t202 = -t101 / 0.2e1;
t13 = qJD(4) * t35 - t176;
t45 = -t114 * t119 + t116 * t67;
t194 = t13 * t45;
t63 = qJD(2) * t67;
t54 = qJD(1) * t63;
t193 = t54 * t66;
t86 = (mrSges(6,2) * t119 - mrSges(6,3) * t116) * qJD(2);
t188 = t86 + (-mrSges(5,1) * t119 + mrSges(5,2) * t116) * qJD(2);
t187 = -t95 + t44;
t181 = Ifges(7,5) * t118;
t180 = Ifges(7,6) * t115;
t159 = -Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t158 = Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1;
t157 = -0.3e1 / 0.2e1 * Ifges(6,6) - 0.3e1 / 0.2e1 * Ifges(5,4);
t156 = t103 * t13 / 0.2e1;
t155 = qJ(5) * t165;
t141 = mrSges(7,1) * t115 + mrSges(7,2) * t118;
t140 = Ifges(7,1) * t118 - t184;
t138 = -Ifges(7,2) * t115 + t183;
t21 = -t115 * t66 + t118 * t45;
t22 = t115 * t45 + t118 * t66;
t46 = t114 * t116 + t119 * t67;
t124 = t36 * mrSges(6,3) + t6 * mrSges(7,2) - t101 * Ifges(7,3) - t85 * Ifges(7,5) - t84 * Ifges(7,6) + t169 * t224 + Ifges(5,5) * t221 + t223 + Ifges(6,4) * t220 + (-t116 * Ifges(6,2) - Ifges(6,6) * t119) * t222 - t30 * mrSges(6,1) - t34 * mrSges(5,3) - t47 * mrSges(5,2) - t5 * mrSges(7,1);
t99 = Ifges(7,3) * t150;
t88 = -qJ(5) * t168 + t106;
t80 = t130 - t197;
t78 = (mrSges(5,1) * t116 + mrSges(5,2) * t119) * t162;
t77 = (-mrSges(6,2) * t116 - mrSges(6,3) * t119) * t162;
t72 = t147 - t155;
t70 = t190 * t167;
t69 = qJD(2) * t135 + t106;
t33 = t100 + (t125 - t155) * qJD(2);
t29 = t107 + t35;
t26 = -mrSges(7,1) * t57 + mrSges(7,2) * t56;
t19 = -t67 * t167 + (qJD(4) * t114 - t64) * t119;
t18 = qJD(4) * t46 - t116 * t64;
t17 = t56 * Ifges(7,1) + t57 * Ifges(7,4) + Ifges(7,5) * t150;
t16 = t56 * Ifges(7,4) + t57 * Ifges(7,2) + Ifges(7,6) * t150;
t15 = t115 * t29 + t118 * t69;
t14 = -t115 * t69 + t118 * t29;
t9 = (qJD(5) - t143) * qJD(4) + t189;
t4 = qJD(6) * t21 + t115 * t18 + t118 * t63;
t3 = -qJD(6) * t22 - t115 * t63 + t118 * t18;
t7 = [t21 * t42 + t22 * t43 + t46 * t26 + t3 * t60 + t4 * t59 + (t77 + t78) * t66 + t188 * t63 - t185 * t18 - t160 * t19 + m(4) * (-t50 * t63 - t51 * t64 - t55 * t67 + t193) + m(5) * (t12 * t46 + t18 * t34 + t19 * t35 + t47 * t63 + t193 + t194) + m(6) * (-t11 * t46 + t18 * t30 - t19 * t31 + t33 * t66 + t36 * t63 + t194) + m(7) * (t1 * t22 + t19 * t25 + t2 * t21 + t3 * t5 + t4 * t6 + t46 * t9) + (-mrSges(4,1) * t63 + mrSges(4,2) * t64 + t219 * qJD(4) * (-t116 * t46 + t119 * t45) + (-mrSges(3,1) * t117 - mrSges(3,2) * t120) * t112 * qJD(2)) * qJD(2); -t54 * mrSges(4,1) + t83 * t26 + t37 * t42 + t38 * t43 - t70 * t44 + t72 * t86 + t80 * t77 + t230 * t60 + t229 * t59 + (qJD(2) * t65 + t55) * mrSges(4,2) + m(6) * (t33 * t80 + t36 * t72) + (t99 / 0.2e1 + t54 * mrSges(5,2) - t33 * mrSges(6,3) + t185 * t65 + t219 * t13 + 0.2e1 * (t208 * t30 + t156) * m(6) + (t208 * t34 + t156) * t209 + (t158 * qJD(4) + (t186 + t216) * t103 + t157 * t169 + t225) * qJD(4) + t215) * t116 + (-t11 * mrSges(6,1) + t12 * mrSges(5,3) - t56 * t139 / 0.2e1 - t57 * t137 / 0.2e1 + t9 * t142 - t54 * mrSges(5,1) + t33 * mrSges(6,2) + t16 * t199 + t17 * t200 + (-t1 * t118 + t2 * t115) * mrSges(7,3) + (t41 * t199 + t115 * t40 / 0.2e1 + (t180 - t181) * t201 + t140 * t206 + t138 * t207 - t25 * t141 + (t6 * t115 + t5 * t118) * mrSges(7,3)) * qJD(6) + (-t124 + ((-t182 / 0.2e1 - t179 / 0.2e1 - t157) * t119 + (-0.3e1 / 0.2e1 * Ifges(6,3) - 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + Ifges(7,3) / 0.2e1) * t116) * qJD(2) + t159 * qJD(4)) * qJD(4) + (qJD(4) * t210 + t211) * t103 + t232 * t65) * t119 + (t78 + t54 * t209 / 0.2e1) * (-pkin(3) - t197) + (qJD(2) * mrSges(4,1) - t188 - m(6) * t36 - t47 * t209 / 0.2e1) * t62 + (t1 * t38 + t2 * t37 + t229 * t6 + t230 * t5 - t25 * t70 + t83 * t9) * m(7) + (t50 * t62 - t51 * t65 + (-t111 * t55 - t113 * t54) * pkin(2)) * m(4); m(7) * (-t6 * t151 + t5 * t152) + (m(7) * (t172 * t5 + t173 * t6) + t219 * qJD(2) * (-t116 ^ 2 - t119 ^ 2)) * qJD(4) + (-m(6) - m(5)) * t119 * t13 + (t26 + m(7) * t9 + (t115 * t59 + t118 * t60 + t210) * qJD(4) + t211) * t116 + (-t232 * qJD(4) + t212 - t228) * t119; t56 * t140 / 0.2e1 + t57 * t138 / 0.2e1 + t9 * t141 + t118 * t17 / 0.2e1 + t16 * t200 - t88 * t86 - t15 * t59 - t14 * t60 - t28 * t44 + qJ(5) * t26 - t11 * mrSges(6,3) - t12 * mrSges(5,2) + t185 * t35 - t186 * t34 + t227 * t13 + t187 * qJD(5) - t217 * mrSges(7,3) + (t136 * t202 + t137 * t207 + t139 * t206 - t233) * qJD(6) + ((t223 + (t181 / 0.2e1 - t180 / 0.2e1 - pkin(4) * mrSges(6,1) + t159) * qJD(4) + t124 + Ifges(6,6) * t231) * t119 + ((-qJ(5) * mrSges(6,1) + t158) * qJD(4) + (Ifges(5,4) / 0.2e1 + Ifges(6,6) / 0.2e1) * t169 + (t224 + Ifges(5,2) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t168 - t225) * t116) * qJD(2) + (qJ(5) * t9 - t14 * t5 - t15 * t6 + t226 * t25) * m(7) + (t134 + t228 + (-m(7) * t144 - t133) * qJD(6)) * t121 + (-pkin(4) * t13 - qJ(5) * t11 + t218 * t31 - t30 * t35 - t36 * t88) * m(6); -t187 * qJD(4) + (mrSges(6,1) * t165 + (-t133 + t86) * t116) * qJD(2) + (-qJD(4) * t25 - t101 * t144 + t217) * m(7) + (qJD(4) * t31 + t169 * t36 + t13) * m(6) - t212; t99 - t25 * (mrSges(7,1) * t85 + mrSges(7,2) * t84) + (Ifges(7,1) * t84 - t198) * t206 + t40 * t205 + (Ifges(7,5) * t84 - Ifges(7,6) * t85) * t202 - t5 * t59 + t6 * t60 + (t5 * t84 + t6 * t85) * mrSges(7,3) + (-Ifges(7,2) * t85 + t41 + t81) * t207 + t215;];
tauc  = t7(:);
