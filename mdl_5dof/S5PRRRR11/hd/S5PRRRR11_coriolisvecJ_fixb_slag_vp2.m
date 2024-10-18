% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:14
% EndTime: 2024-09-27 21:45:18
% DurationCPUTime: 2.20s
% Computational Cost: add. (4306->328), mult. (13928->485), div. (0->0), fcn. (10558->8), ass. (0->166)
t155 = cos(pkin(5));
t146 = t155 * qJD(2) + qJD(3);
t142 = qJD(4) + t146;
t138 = qJD(5) + t142;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t185 = qJD(4) * t160;
t186 = qJD(4) * t157;
t158 = sin(qJ(3));
t154 = sin(pkin(5));
t219 = pkin(7) + pkin(8);
t181 = t154 * t219;
t172 = t158 * t181;
t161 = cos(qJ(3));
t212 = pkin(2) * t155;
t149 = t161 * t212;
t189 = qJD(1) * t154;
t191 = qJD(2) * t149 + t161 * t189;
t95 = -qJD(2) * t172 + t191;
t85 = pkin(3) * t146 + t95;
t169 = qJD(3) * t172;
t192 = t191 * qJD(3);
t88 = -qJD(2) * t169 + t192;
t148 = t158 * t212;
t164 = -t161 * t181 - t148;
t177 = t158 * t189;
t89 = (t164 * qJD(2) - t177) * qJD(3);
t195 = t154 * t161;
t96 = t177 + (t219 * t195 + t148) * qJD(2);
t17 = t157 * t89 + t160 * t88 + t85 * t185 - t96 * t186;
t126 = (t157 * t161 + t158 * t160) * t154;
t226 = qJD(3) + qJD(4);
t87 = t226 * t126;
t74 = qJD(2) * t87;
t10 = -pkin(9) * t74 + t17;
t92 = t160 * t96;
t48 = t157 * t85 + t92;
t18 = -t48 * qJD(4) - t157 * t88 + t160 * t89;
t125 = (-t157 * t158 + t160 * t161) * t154;
t86 = t226 * t125;
t73 = qJD(2) * t86;
t11 = -pkin(9) * t73 + t18;
t156 = sin(qJ(5));
t159 = cos(qJ(5));
t121 = qJD(2) * t125;
t210 = pkin(9) * t121;
t41 = t48 + t210;
t198 = t156 * t41;
t122 = qJD(2) * t126;
t119 = t122 * pkin(9);
t90 = t157 * t96;
t47 = t160 * t85 - t90;
t40 = -t119 + t47;
t35 = pkin(4) * t142 + t40;
t8 = t159 * t35 - t198;
t2 = t8 * qJD(5) + t10 * t159 + t11 * t156;
t173 = t159 * t121 - t122 * t156;
t29 = t173 * qJD(5) - t156 * t74 + t159 * t73;
t197 = t159 * t41;
t9 = t156 * t35 + t197;
t3 = -t9 * qJD(5) - t10 * t156 + t11 * t159;
t70 = t121 * t156 + t122 * t159;
t30 = -t70 * qJD(5) - t156 * t73 - t159 * t74;
t170 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t29 + Ifges(6,6) * t30;
t213 = Ifges(6,4) * t70;
t66 = Ifges(6,4) * t173;
t38 = Ifges(6,1) * t70 + Ifges(6,5) * t138 + t66;
t135 = (-pkin(3) * t161 - pkin(2)) * t154;
t152 = t155 * qJD(1);
t124 = qJD(2) * t135 + t152;
t81 = -t121 * pkin(4) + t124;
t238 = t170 - (Ifges(6,5) * t173 - Ifges(6,6) * t70) * t138 / 0.2e1 - (-Ifges(6,2) * t70 + t38 + t66) * t173 / 0.2e1 - t81 * (mrSges(6,1) * t70 + mrSges(6,2) * t173) - (Ifges(6,1) * t173 - t213) * t70 / 0.2e1;
t217 = t70 * t9;
t37 = Ifges(6,2) * t173 + Ifges(6,6) * t138 + t213;
t236 = t37 / 0.2e1;
t235 = -Ifges(4,6) * t146 / 0.2e1;
t150 = pkin(3) * t160 + pkin(4);
t183 = qJD(5) * t159;
t184 = qJD(5) * t156;
t193 = t157 * t159;
t53 = -t157 * t95 - t92;
t45 = t53 - t210;
t54 = t160 * t95 - t90;
t46 = -t119 + t54;
t234 = t156 * t46 - t159 * t45 - t150 * t184 + (-t157 * t183 + (-t156 * t160 - t193) * qJD(4)) * pkin(3);
t194 = t156 * t157;
t233 = -t156 * t45 - t159 * t46 + t150 * t183 + (-t157 * t184 + (t159 * t160 - t194) * qJD(4)) * pkin(3);
t231 = mrSges(6,3) * t173;
t211 = pkin(3) * t155;
t111 = t149 - t172 + t211;
t190 = pkin(7) * t195 + t148;
t120 = pkin(8) * t195 + t190;
t64 = t157 * t111 + t160 * t120;
t182 = qJD(2) * qJD(3);
t174 = t158 * t182;
t171 = t154 * t174;
t104 = -pkin(7) * t171 + t192;
t115 = qJD(2) * t190 + t177;
t227 = qJD(3) * t115;
t228 = -mrSges(4,1) * t227 - t104 * mrSges(4,2);
t188 = qJD(2) * t154;
t136 = -pkin(2) * t188 + t152;
t205 = Ifges(4,4) * t158;
t225 = (m(5) * t124 - mrSges(5,1) * t121 + mrSges(5,2) * t122) * pkin(3) + t136 * mrSges(4,1) + t235 - (t161 * Ifges(4,2) + t205) * t188 / 0.2e1 - t115 * mrSges(4,3);
t224 = t18 * mrSges(5,1) - t17 * mrSges(5,2) + Ifges(5,5) * t73 - Ifges(5,6) * t74;
t222 = t173 / 0.2e1;
t220 = t70 / 0.2e1;
t215 = t122 / 0.2e1;
t79 = t125 * t159 - t126 * t156;
t208 = t79 * t29;
t80 = t125 * t156 + t126 * t159;
t207 = t80 * t30;
t206 = mrSges(5,3) * t121;
t200 = t122 * mrSges(5,3);
t199 = t122 * Ifges(5,4);
t187 = qJD(3) * t158;
t180 = t158 * t188;
t179 = t161 * t188;
t178 = t154 * t187;
t175 = t154 * t182;
t63 = t160 * t111 - t120 * t157;
t49 = pkin(4) * t155 - pkin(9) * t126 + t63;
t52 = pkin(9) * t125 + t64;
t21 = -t156 * t52 + t159 * t49;
t22 = t156 * t49 + t159 * t52;
t144 = qJD(3) * t149;
t112 = t144 - t169;
t113 = t164 * qJD(3);
t33 = t111 * t185 + t160 * t112 + t157 * t113 - t120 * t186;
t114 = -pkin(7) * t180 + t191;
t141 = Ifges(4,4) * t179;
t163 = Ifges(4,1) * t180 / 0.2e1 + t141 / 0.2e1 + Ifges(4,5) * t146 + t136 * mrSges(4,2) - t114 * mrSges(4,3);
t34 = -t64 * qJD(4) - t112 * t157 + t160 * t113;
t118 = Ifges(5,4) * t121;
t61 = t121 * Ifges(5,2) + t142 * Ifges(5,6) + t199;
t62 = t122 * Ifges(5,1) + t142 * Ifges(5,5) + t118;
t162 = t8 * t231 - t124 * (mrSges(5,1) * t122 + mrSges(5,2) * t121) + t70 * t236 + t47 * t206 + t61 * t215 - t122 * (Ifges(5,1) * t121 - t199) / 0.2e1 - t142 * (Ifges(5,5) * t121 - Ifges(5,6) * t122) / 0.2e1 + t224 - (-Ifges(5,2) * t122 + t118 + t62) * t121 / 0.2e1 + t238;
t137 = Ifges(4,5) * t161 * t175;
t133 = -pkin(7) * t154 * t158 + t149;
t132 = pkin(3) * t193 + t150 * t156;
t131 = -pkin(3) * t194 + t150 * t159;
t130 = t190 * qJD(3);
t129 = -pkin(7) * t178 + t144;
t128 = -t146 * mrSges(4,2) + mrSges(4,3) * t179;
t127 = mrSges(4,1) * t146 - mrSges(4,3) * t180;
t123 = (mrSges(4,1) * t158 + mrSges(4,2) * t161) * t175;
t100 = pkin(3) * t180 + pkin(4) * t122;
t99 = -t125 * pkin(4) + t135;
t98 = mrSges(5,1) * t142 - t200;
t97 = -mrSges(5,2) * t142 + t206;
t65 = pkin(3) * t178 + pkin(4) * t87;
t59 = pkin(3) * t171 + pkin(4) * t74;
t58 = mrSges(6,1) * t138 - mrSges(6,3) * t70;
t57 = -mrSges(6,2) * t138 + t231;
t44 = mrSges(5,1) * t74 + mrSges(5,2) * t73;
t43 = -mrSges(6,1) * t173 + mrSges(6,2) * t70;
t32 = -t80 * qJD(5) - t156 * t86 - t159 * t87;
t31 = t79 * qJD(5) - t156 * t87 + t159 * t86;
t24 = -pkin(9) * t86 + t34;
t23 = -pkin(9) * t87 + t33;
t13 = t159 * t40 - t198;
t12 = -t156 * t40 - t197;
t6 = -mrSges(6,1) * t30 + mrSges(6,2) * t29;
t5 = -t22 * qJD(5) - t156 * t23 + t159 * t24;
t4 = t21 * qJD(5) + t156 * t24 + t159 * t23;
t1 = [t31 * t57 + t32 * t58 + t86 * t97 - t87 * t98 + (t207 - t208) * mrSges(6,3) + (-t125 * t73 - t126 * t74) * mrSges(5,3) + (t6 + t44 + t123) * t155 + m(5) * (t125 * t18 + t126 * t17 - t47 * t87 + t48 * t86) + m(6) * (t155 * t59 + t2 * t80 + t3 * t79 + t31 * t9 + t32 * t8) + ((-t127 * t158 + t128 * t161 + (-t158 ^ 2 - t161 ^ 2) * mrSges(4,3) * t188) * qJD(3) + m(4) * (t104 * t158 - t114 * t187) + m(5) * t174 * t211) * t154; t31 * t38 / 0.2e1 + t4 * t57 + t5 * t58 + t65 * t43 + t59 * (-mrSges(6,1) * t79 + mrSges(6,2) * t80) + t81 * (-mrSges(6,1) * t32 + mrSges(6,2) * t31) + t86 * t62 / 0.2e1 - t87 * t61 / 0.2e1 + t33 * t97 + t34 * t98 + t99 * t6 + t129 * t128 - t130 * t127 + t135 * t44 + t138 * (Ifges(6,5) * t31 + Ifges(6,6) * t32) / 0.2e1 + t32 * t236 + (t2 * t79 - t21 * t29 + t22 * t30 - t3 * t80 - t8 * t31 + t9 * t32) * mrSges(6,3) + (-pkin(2) * t123 + (t104 * t161 + t158 * t227) * mrSges(4,3) + (t163 * t161 + (t235 + t225) * t158) * qJD(3)) * t154 + m(4) * (t104 * t190 - t114 * t130 + t115 * t129 - t133 * t227) + (-t47 * t86 - t48 * t87 - t63 * t73 - t64 * t74) * mrSges(5,3) + (-mrSges(5,3) * t18 + Ifges(5,1) * t73 - Ifges(5,4) * t74) * t126 + (mrSges(5,3) * t17 + Ifges(5,4) * t73 - Ifges(5,2) * t74) * t125 + (Ifges(5,1) * t86 - Ifges(5,4) * t87) * t215 + t121 * (Ifges(5,4) * t86 - Ifges(5,2) * t87) / 0.2e1 + t124 * (mrSges(5,1) * t87 + mrSges(5,2) * t86) + t142 * (Ifges(5,5) * t86 - Ifges(5,6) * t87) / 0.2e1 + m(6) * (t2 * t22 + t21 * t3 + t4 * t9 + t5 * t8 + t59 * t99 + t65 * t81) + ((Ifges(4,5) * t155 / 0.2e1 - t133 * mrSges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,4) * t195) * t161 + (-t190 * mrSges(4,3) - Ifges(4,6) * t155 + (-0.3e1 / 0.2e1 * t205 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t161) * t154 + (m(5) * t135 - mrSges(5,1) * t125 + mrSges(5,2) * t126) * pkin(3)) * t158) * t175 + m(5) * (t17 * t64 + t18 * t63 + t33 * t48 + t34 * t47) + (t137 / 0.2e1 + t224 + t170 + t228) * t155 + (t31 * t220 + t80 * t29) * Ifges(6,1) + (t32 * t222 + t79 * t30) * Ifges(6,2) + (t32 * t220 + t31 * t222 + t207 + t208) * Ifges(6,4); t48 * t200 - t54 * t97 - t53 * t98 - t100 * t43 + t115 * t127 - t114 * t128 + t234 * t58 + t233 * t57 - m(5) * (t47 * t53 + t48 * t54) + (-t131 * t29 + t132 * t30 + t217) * mrSges(6,3) + ((-t141 / 0.2e1 - t163) * t161 + ((t205 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t161) * t188 + (t146 / 0.2e1 - qJD(3)) * Ifges(4,6) - t225) * t158) * t188 + t162 + t137 + (-t98 * t186 + m(5) * (t157 * t17 + t160 * t18 + t48 * t185 - t47 * t186) + t97 * t185 + (-t157 * t74 - t160 * t73) * mrSges(5,3)) * pkin(3) + (-t100 * t81 + t131 * t3 + t132 * t2 + t233 * t9 + t234 * t8) * m(6) + t228; mrSges(6,3) * t217 - t13 * t57 - t12 * t58 - t47 * t97 + (t98 + t200) * t48 + t162 - m(6) * (t12 * t8 + t13 * t9) + (-t122 * t43 + (-t156 * t58 + t159 * t57) * qJD(5) + (t156 * t30 - t159 * t29) * mrSges(6,3) + (-t122 * t81 + t156 * t2 + t159 * t3 + t9 * t183 - t8 * t184) * m(6)) * pkin(4); t37 * t220 - t8 * t57 + t9 * t58 + (t173 * t8 + t217) * mrSges(6,3) + t238;];
tauc = t1(:);
