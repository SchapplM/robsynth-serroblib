% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2018-11-23 14:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:58:52
% EndTime: 2018-11-23 14:58:55
% DurationCPUTime: 3.43s
% Computational Cost: add. (2241->373), mult. (5044->497), div. (0->0), fcn. (2875->8), ass. (0->185)
t105 = cos(qJ(6));
t198 = t105 / 0.2e1;
t221 = -qJD(2) / 0.2e1;
t102 = sin(qJ(6));
t183 = Ifges(7,4) * t102;
t128 = Ifges(7,2) * t105 + t183;
t182 = Ifges(7,4) * t105;
t130 = Ifges(7,1) * t102 + t182;
t133 = mrSges(7,1) * t105 - mrSges(7,2) * t102;
t108 = -pkin(4) - pkin(9);
t103 = sin(qJ(4));
t101 = cos(pkin(6));
t169 = qJD(1) * t101;
t146 = t103 * t169;
t106 = cos(qJ(4));
t157 = t106 * qJD(2);
t115 = pkin(5) * t157 + t146;
t109 = -pkin(2) - pkin(8);
t107 = cos(qJ(2));
t100 = sin(pkin(6));
t170 = qJD(1) * t100;
t147 = t107 * t170;
t126 = qJD(3) - t147;
t58 = qJD(2) * t109 + t126;
t53 = t106 * t58;
t114 = t115 - t53;
t222 = qJD(5) + t114;
t20 = qJD(4) * t108 + t222;
t140 = -qJ(5) * t106 + qJ(3);
t120 = pkin(9) * t103 + t140;
t104 = sin(qJ(2));
t148 = t104 * t170;
t168 = qJD(2) * t103;
t136 = pkin(4) * t168 + t148;
t37 = qJD(2) * t120 + t136;
t5 = -t102 * t37 + t105 * t20;
t6 = t102 * t20 + t105 * t37;
t134 = t5 * t102 - t6 * t105;
t178 = Ifges(7,6) * t105;
t181 = Ifges(7,5) * t102;
t199 = -t102 / 0.2e1;
t92 = qJD(6) + t157;
t203 = -t92 / 0.2e1;
t162 = qJD(4) * t105;
t72 = t102 * t168 + t162;
t205 = -t72 / 0.2e1;
t164 = qJD(4) * t102;
t71 = t105 * t168 - t164;
t206 = -t71 / 0.2e1;
t197 = Ifges(7,4) * t72;
t22 = Ifges(7,2) * t71 + Ifges(7,6) * t92 + t197;
t66 = Ifges(7,4) * t71;
t23 = Ifges(7,1) * t72 + Ifges(7,5) * t92 + t66;
t165 = qJD(4) * qJ(5);
t145 = t106 * t169;
t29 = t145 + (-pkin(5) * qJD(2) + t58) * t103;
t24 = t29 + t165;
t225 = -t134 * mrSges(7,3) - (t178 + t181) * t203 - t128 * t206 - t130 * t205 - t24 * t133 - t23 * t199 + t22 * t198;
t218 = mrSges(6,1) + mrSges(5,3);
t224 = mrSges(6,2) - mrSges(5,1);
t219 = qJD(4) / 0.2e1;
t220 = -qJD(4) / 0.2e1;
t41 = t103 * t58 + t145;
t33 = -t41 - t165;
t51 = qJD(2) * t140 + t136;
t76 = qJD(2) * qJ(3) + t148;
t223 = t76 * mrSges(5,1) + t33 * mrSges(6,1) - t51 * mrSges(6,2) - t41 * mrSges(5,3) - Ifges(6,5) * t220 - Ifges(5,6) * t219 + t225 + ((Ifges(5,4) + Ifges(6,6)) * t106 + (-Ifges(5,2) - Ifges(6,3)) * t103) * t221;
t161 = qJD(4) * t106;
t163 = qJD(4) * t103;
t217 = pkin(4) * t161 + qJ(5) * t163;
t176 = qJD(4) * t41;
t167 = qJD(2) * t104;
t144 = t100 * t167;
t139 = qJD(1) * t144;
t78 = t106 * t139;
t18 = -t78 + t176;
t216 = -t33 * qJD(4) - t18;
t40 = -t53 + t146;
t215 = -qJD(5) - t40;
t15 = qJD(4) * t29 - t78;
t116 = qJD(3) + (qJD(4) * pkin(9) - qJD(5)) * t106;
t156 = qJD(2) * qJD(4);
t142 = t103 * t156;
t97 = pkin(4) * t157;
t184 = qJ(5) * t142 + qJD(4) * t97;
t19 = (t116 + t147) * qJD(2) + t184;
t1 = qJD(6) * t5 + t102 * t15 + t105 * t19;
t2 = -qJD(6) * t6 - t102 * t19 + t105 * t15;
t135 = t1 * t102 + t2 * t105;
t214 = -t24 * qJD(4) + t135;
t213 = m(5) * t41 - m(6) * t33;
t155 = qJD(4) * qJD(6);
t158 = qJD(6) * t105;
t44 = -t102 * t155 + (t102 * t161 + t103 * t158) * qJD(2);
t159 = qJD(6) * t102;
t45 = -t105 * t155 + (-t103 * t159 + t105 * t161) * qJD(2);
t212 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t44 + Ifges(7,6) * t45;
t49 = -mrSges(7,2) * t92 + mrSges(7,3) * t71;
t50 = mrSges(7,1) * t92 - mrSges(7,3) * t72;
t123 = t102 * t50 - t105 * t49;
t25 = -mrSges(7,1) * t142 - mrSges(7,3) * t44;
t26 = mrSges(7,2) * t142 + mrSges(7,3) * t45;
t124 = t102 * t26 + t105 * t25;
t211 = -t123 * qJD(6) + t124;
t208 = t44 / 0.2e1;
t207 = t45 / 0.2e1;
t204 = t72 / 0.2e1;
t200 = m(6) * t51;
t174 = t100 * t107;
t56 = t101 * t103 + t106 * t174;
t195 = t18 * t56;
t190 = pkin(5) - t109;
t189 = t103 * t139 + t58 * t161;
t74 = (-mrSges(6,2) * t103 - mrSges(6,3) * t106) * qJD(2);
t75 = (mrSges(5,1) * t103 + mrSges(5,2) * t106) * qJD(2);
t188 = t74 + t75;
t187 = -t224 * qJD(4) - t157 * t218;
t34 = -mrSges(7,1) * t71 + mrSges(7,2) * t72;
t87 = mrSges(6,1) * t168 - qJD(4) * mrSges(6,3);
t186 = -t87 + t34;
t185 = qJD(4) * mrSges(5,2) + mrSges(5,3) * t168 + t87;
t73 = qJ(5) * t168 + t97;
t180 = Ifges(7,5) * t105;
t179 = Ifges(7,6) * t102;
t177 = t107 * t76;
t175 = t100 * t104;
t173 = t104 * t106;
t166 = qJD(2) * t107;
t160 = qJD(5) * t106;
t154 = t34 - t185;
t153 = -0.3e1 / 0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * Ifges(6,6);
t152 = Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t151 = -Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1;
t150 = t103 * t174;
t149 = -t109 * t18 / 0.2e1;
t143 = t100 * t166;
t141 = qJD(4) * t190;
t138 = t148 / 0.2e1;
t132 = mrSges(7,1) * t102 + mrSges(7,2) * t105;
t131 = Ifges(7,1) * t105 - t183;
t129 = -Ifges(7,2) * t102 + t182;
t125 = qJD(3) + t147;
t99 = t103 * pkin(4);
t65 = t120 + t99;
t81 = t190 * t106;
t31 = t102 * t81 + t105 * t65;
t30 = -t102 * t65 + t105 * t81;
t69 = t125 * qJD(2);
t117 = t104 * t69 + t166 * t76;
t39 = t102 * t56 + t105 * t175;
t38 = -t102 * t175 + t105 * t56;
t113 = -m(7) * t134 - t123;
t32 = -qJD(4) * pkin(4) - t215;
t95 = Ifges(6,6) * t168;
t112 = t51 * mrSges(6,3) + t6 * mrSges(7,2) - t92 * Ifges(7,3) - t72 * Ifges(7,5) - t71 * Ifges(7,6) + Ifges(5,5) * t220 + (t106 * Ifges(5,1) - Ifges(5,4) * t103) * t221 + Ifges(6,4) * t219 - Ifges(6,2) * t157 / 0.2e1 + t95 / 0.2e1 - t32 * mrSges(6,1) - t40 * mrSges(5,3) - t5 * mrSges(7,1) - t76 * mrSges(5,2);
t110 = qJD(2) ^ 2;
t80 = t190 * t103;
t79 = t140 + t99;
t70 = -qJD(2) * pkin(2) + t126;
t68 = t106 * t141;
t67 = t103 * t141;
t64 = (mrSges(5,1) * t106 - mrSges(5,2) * t103) * t156;
t63 = (-mrSges(6,2) * t106 + mrSges(6,3) * t103) * t156;
t57 = t101 * t106 - t150;
t55 = pkin(9) * t157 + t73;
t54 = qJD(3) - t160 + t217;
t48 = (-t102 * t173 + t105 * t107) * t170;
t47 = (-t102 * t107 - t105 * t173) * t170;
t46 = t116 + t217;
t36 = -qJD(4) * t150 + t101 * t161 - t106 * t144;
t35 = qJD(4) * t56 - t103 * t144;
t27 = (t125 - t160) * qJD(2) + t184;
t17 = -qJD(4) * t146 + t189;
t16 = (-qJD(5) + t146) * qJD(4) - t189;
t14 = (qJD(5) - t115) * qJD(4) + t189;
t13 = -mrSges(7,1) * t45 + mrSges(7,2) * t44;
t12 = t102 * t29 + t105 * t55;
t11 = -t102 * t55 + t105 * t29;
t10 = t44 * Ifges(7,1) + t45 * Ifges(7,4) - Ifges(7,5) * t142;
t9 = t44 * Ifges(7,4) + t45 * Ifges(7,2) - Ifges(7,6) * t142;
t8 = -qJD(6) * t39 - t102 * t143 + t105 * t36;
t7 = qJD(6) * t38 + t102 * t36 + t105 * t143;
t4 = -qJD(6) * t31 - t102 * t46 - t105 * t67;
t3 = qJD(6) * t30 - t102 * t67 + t105 * t46;
t21 = [t57 * t13 + t38 * t25 + t39 * t26 + t7 * t49 + t8 * t50 - t187 * t36 - t154 * t35 + m(5) * (t17 * t57 - t35 * t41 + t36 * t40 + t195) + m(6) * (-t16 * t57 + t32 * t36 + t33 * t35 + t195) + m(7) * (t1 * t39 + t14 * t57 + t2 * t38 - t24 * t35 + t5 * t8 + t6 * t7) + t218 * t156 * (-t103 * t56 - t106 * t57) + (((-mrSges(3,2) + mrSges(4,3)) * t110 + t188 * qJD(2)) * t107 + (t63 + t64 + (-mrSges(3,1) + mrSges(4,2)) * t110) * t104 + m(5) * t117 + m(6) * (t104 * t27 + t166 * t51) + (-t107 * t139 + t167 * t70 + t117) * m(4)) * t100; qJ(3) * t64 + qJD(3) * t75 - t80 * t13 + t30 * t25 + t31 * t26 - t68 * t34 + t54 * t74 + t79 * t63 + (-t47 + t4) * t50 + (-t48 + t3) * t49 - t188 * t147 + (qJD(2) * t126 + t69) * mrSges(4,3) - m(7) * (t47 * t5 + t48 * t6) + m(7) * (t1 * t31 - t14 * t80 + t2 * t30 - t24 * t68 + t3 * t6 + t4 * t5) + m(6) * (t27 * t79 + t51 * t54) + 0.2e1 * (-m(5) * t177 / 0.2e1 - t107 * t200 / 0.2e1 - (pkin(2) * t167 + t104 * t70 + t177) * m(4) / 0.2e1) * t170 + (t69 * mrSges(5,2) - t27 * mrSges(6,3) + t218 * t18 - t187 * t148 + 0.2e1 * (t138 * t32 + t149) * m(6) + 0.2e1 * (t138 * t40 + t149) * m(5) + ((-t185 + t213) * t109 + t151 * qJD(4) + t153 * t157 + t223) * qJD(4) + t212) * t106 + (-t27 * mrSges(6,2) + t9 * t198 + t102 * t10 / 0.2e1 - t17 * mrSges(5,3) + t16 * mrSges(6,1) - t14 * t133 + t130 * t208 + t128 * t207 + t69 * mrSges(5,1) + (t1 * t105 - t2 * t102) * mrSges(7,3) + (t92 * (-t179 + t180) / 0.2e1 + t131 * t204 + t71 * t129 / 0.2e1 + t24 * t132 + t23 * t198 + t22 * t199 + (-t6 * t102 - t5 * t105) * mrSges(7,3)) * qJD(6) + (t152 * qJD(4) + ((-t181 / 0.2e1 - t178 / 0.2e1 - t153) * t103 + (0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(6,3) - Ifges(7,3) - 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(6,2)) * t106) * qJD(2) + t112) * qJD(4) + (m(5) * t17 - m(6) * t16 + (m(5) * t40 + m(6) * t32 - t187) * qJD(4)) * t109 + (-m(7) * t24 - t154 - t213) * t148) * t103 + (m(5) + m(4)) * (qJ(3) * t69 + qJD(3) * t76); -t110 * mrSges(4,3) + (t13 + (t102 * t49 + t105 * t50 - t187) * qJD(4) + m(7) * (t162 * t5 + t164 * t6 + t14) + m(5) * (t40 * qJD(4) + t17) + m(6) * (t32 * qJD(4) - t16)) * t103 + (t154 * qJD(4) + m(7) * (-t158 * t6 + t159 * t5 - t214) + m(5) * (-t18 + t176) + m(6) * t216 - t211) * t106 + (-m(5) * t76 - t200 + (-t76 + t148) * m(4) - t113 - t188) * qJD(2); t131 * t208 + t129 * t207 + t14 * t132 + t10 * t198 + t9 * t199 - t73 * t74 - t12 * t49 - t11 * t50 + t114 * t34 - t16 * mrSges(6,3) - t17 * mrSges(5,2) + qJ(5) * t13 + t187 * t41 - t185 * t40 + t224 * t18 + t124 * t108 + t186 * qJD(5) - t135 * mrSges(7,3) + (t108 * t113 - t225) * qJD(6) + ((-t95 / 0.2e1 + (-t180 / 0.2e1 + t179 / 0.2e1 + pkin(4) * mrSges(6,1) + t152) * qJD(4) - t112 - Ifges(5,4) * t168 / 0.2e1) * t103 + ((-qJ(5) * mrSges(6,1) + t151) * qJD(4) + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t157 + (Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t168 - t223) * t106) * qJD(2) + (qJ(5) * t14 + t135 * t108 - t11 * t5 - t12 * t6 + t222 * t24) * m(7) + (-pkin(4) * t18 - qJ(5) * t16 + t215 * t33 - t32 * t41 - t51 * t73) * m(6); -t186 * qJD(4) + (-mrSges(6,1) * t163 + (-t123 + t74) * t106) * qJD(2) + (-t134 * t92 + t214) * m(7) + (t157 * t51 - t216) * m(6) + t211; -Ifges(7,3) * t142 - t24 * (mrSges(7,1) * t72 + mrSges(7,2) * t71) + (Ifges(7,1) * t71 - t197) * t205 + t22 * t204 + (Ifges(7,5) * t71 - Ifges(7,6) * t72) * t203 - t5 * t49 + t6 * t50 + (t5 * t71 + t6 * t72) * mrSges(7,3) + (-Ifges(7,2) * t72 + t23 + t66) * t206 + t212;];
tauc  = t21(:);
