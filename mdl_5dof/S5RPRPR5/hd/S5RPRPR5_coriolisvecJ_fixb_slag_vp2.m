% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR5
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:26
% EndTime: 2020-01-03 11:41:43
% DurationCPUTime: 4.27s
% Computational Cost: add. (3787->343), mult. (10388->503), div. (0->0), fcn. (7441->8), ass. (0->180)
t157 = sin(pkin(9));
t159 = cos(pkin(9));
t162 = sin(qJ(3));
t164 = cos(qJ(3));
t135 = t157 * t164 + t159 * t162;
t128 = t135 * qJD(3);
t160 = cos(pkin(8));
t171 = qJD(1) * t135;
t223 = t160 * t171 - t128;
t172 = t157 * t162 - t159 * t164;
t129 = t172 * qJD(3);
t190 = qJD(1) * t160;
t222 = t172 * t190 - t129;
t158 = sin(pkin(8));
t113 = t158 * t128;
t104 = qJD(1) * t113;
t199 = qJ(4) * t158;
t184 = t164 * t199;
t201 = qJ(2) * t162;
t186 = t160 * t201;
t170 = -t184 - t186;
t187 = qJD(4) * t158;
t167 = qJD(3) * t170 - t162 * t187;
t137 = -pkin(2) * t160 - pkin(6) * t158 - pkin(1);
t126 = qJD(1) * t137 + qJD(2);
t195 = t160 * t164;
t145 = qJD(2) * t195;
t188 = qJD(3) * t164;
t194 = qJD(1) * t145 + t126 * t188;
t49 = qJD(1) * t167 + t194;
t146 = qJ(2) * t195;
t189 = qJD(2) * t162;
t180 = t160 * t189;
t169 = -t164 * t187 - t180;
t197 = t126 * t162;
t179 = qJD(3) * t197;
t196 = t158 * t162;
t185 = qJ(4) * t196;
t50 = -t179 + ((-t146 + t185) * qJD(3) + t169) * qJD(1);
t14 = -t157 * t49 + t159 * t50;
t12 = pkin(7) * t104 + t14;
t111 = t158 * t129;
t103 = qJD(1) * t111;
t15 = t157 * t50 + t159 * t49;
t13 = pkin(7) * t103 + t15;
t161 = sin(qJ(5));
t163 = cos(qJ(5));
t147 = qJD(3) - t190;
t191 = qJD(1) * t158;
t181 = t164 * t191;
t182 = t162 * t191;
t108 = t157 * t182 - t159 * t181;
t207 = pkin(7) * t108;
t115 = t164 * t126;
t79 = qJD(1) * t170 + t115;
t68 = pkin(3) * t147 + t79;
t183 = qJ(2) * t190;
t95 = t164 * t183 + t197;
t80 = -qJ(4) * t182 + t95;
t72 = t157 * t80;
t35 = t159 * t68 - t72;
t18 = pkin(4) * t147 + t207 + t35;
t109 = t158 * t171;
t206 = pkin(7) * t109;
t202 = t159 * t80;
t36 = t157 * t68 + t202;
t19 = t36 - t206;
t7 = t161 * t18 + t163 * t19;
t3 = -qJD(5) * t7 + t12 * t163 - t13 * t161;
t1 = t3 * mrSges(6,1);
t142 = qJD(5) + t147;
t178 = t108 * t161 - t163 * t109;
t57 = t108 * t163 + t109 * t161;
t209 = Ifges(6,4) * t57;
t20 = Ifges(6,2) * t178 + Ifges(6,6) * t142 - t209;
t53 = Ifges(6,4) * t178;
t21 = -Ifges(6,1) * t57 + Ifges(6,5) * t142 + t53;
t217 = -t57 / 0.2e1;
t26 = qJD(5) * t57 + t103 * t163 + t104 * t161;
t23 = Ifges(6,6) * t26;
t25 = qJD(5) * t178 + t103 * t161 - t104 * t163;
t24 = Ifges(6,5) * t25;
t127 = pkin(3) * t182 + qJ(2) * t191 + qJD(4);
t78 = pkin(4) * t109 + t127;
t235 = t1 + t23 + t24 - (Ifges(6,5) * t178 + Ifges(6,6) * t57) * t142 / 0.2e1 - (Ifges(6,2) * t57 + t21 + t53) * t178 / 0.2e1 - t78 * (-t57 * mrSges(6,1) + mrSges(6,2) * t178) + t20 * t217 + (Ifges(6,1) * t178 + t209) * t57 / 0.2e1;
t153 = t158 * qJD(2);
t131 = t158 * pkin(3) * t188 + t153;
t6 = -t161 * t19 + t163 * t18;
t234 = t178 * t6 - t57 * t7;
t86 = t135 * t163 - t161 * t172;
t231 = -qJD(5) * t86 - t222 * t161 + t223 * t163;
t85 = -t135 * t161 - t163 * t172;
t230 = qJD(5) * t85 + t223 * t161 + t222 * t163;
t152 = pkin(3) * t159 + pkin(4);
t208 = pkin(3) * t157;
t124 = t152 * t161 + t163 * t208;
t38 = -t157 * t79 - t202;
t27 = t38 + t206;
t39 = t159 * t79 - t72;
t28 = t39 + t207;
t229 = -t124 * qJD(5) + t161 * t28 - t163 * t27;
t123 = t152 * t163 - t161 * t208;
t228 = t123 * qJD(5) - t161 * t27 - t163 * t28;
t155 = t158 ^ 2;
t156 = t160 ^ 2;
t192 = t155 + t156;
t224 = t192 * mrSges(3,3);
t107 = t162 * t137 + t146;
t94 = -t162 * t183 + t115;
t174 = t94 * t162 - t95 * t164;
t176 = -Ifges(4,5) * t162 - Ifges(4,6) * t164;
t204 = Ifges(4,4) * t164;
t205 = Ifges(4,4) * t162;
t210 = -t164 / 0.2e1;
t212 = t147 / 0.2e1;
t220 = t174 * mrSges(4,3) + (t147 * Ifges(4,6) + (-t162 * Ifges(4,2) + t204) * t191) * t210 - t162 * (t147 * Ifges(4,5) + (t164 * Ifges(4,1) - t205) * t191) / 0.2e1 + t176 * t212;
t2 = qJD(5) * t6 + t12 * t161 + t13 * t163;
t216 = t2 * mrSges(6,2);
t215 = -t108 / 0.2e1;
t211 = -t160 / 0.2e1;
t193 = t137 * t188 + t145;
t65 = t167 + t193;
t66 = (-t146 + (-t137 + t199) * t162) * qJD(3) + t169;
t31 = t157 * t66 + t159 * t65;
t133 = t164 * t137;
t84 = -t184 + t133 + (-pkin(3) - t201) * t160;
t96 = -t185 + t107;
t43 = t157 * t84 + t159 * t96;
t203 = Ifges(5,4) * t108;
t165 = qJD(1) ^ 2;
t200 = qJ(2) * t165;
t125 = t131 * qJD(1);
t136 = pkin(3) * t196 + t158 * qJ(2);
t30 = -t157 * t65 + t159 * t66;
t42 = -t157 * t96 + t159 * t84;
t177 = qJD(3) * t186;
t119 = t172 * t158;
t34 = -pkin(4) * t160 + pkin(7) * t119 + t42;
t118 = t135 * t158;
t37 = -pkin(7) * t118 + t43;
t10 = -t161 * t37 + t163 * t34;
t11 = t161 * t34 + t163 * t37;
t76 = -qJD(1) * t177 + t194;
t77 = -t179 + (-qJ(2) * t188 - t189) * t190;
t175 = -t76 * t162 - t77 * t164;
t69 = -t118 * t163 + t119 * t161;
t70 = -t118 * t161 - t119 * t163;
t120 = -mrSges(4,2) * t147 - mrSges(4,3) * t182;
t121 = mrSges(4,1) * t147 - mrSges(4,3) * t181;
t173 = t120 * t164 - t121 * t162;
t168 = -t77 * mrSges(4,1) - t14 * mrSges(5,1) + t76 * mrSges(4,2) + t15 * mrSges(5,2) + t216;
t150 = t155 * t200;
t122 = (mrSges(4,1) * t162 + mrSges(4,2) * t164) * t191;
t106 = t133 - t186;
t105 = Ifges(5,4) * t109;
t101 = Ifges(5,5) * t104;
t100 = Ifges(5,6) * t103;
t99 = t104 * mrSges(5,2);
t93 = -qJD(3) * t107 - t180;
t92 = -t177 + t193;
t91 = pkin(3) * t181 - pkin(4) * t108;
t89 = mrSges(5,1) * t147 + mrSges(5,3) * t108;
t88 = -mrSges(5,2) * t147 - t109 * mrSges(5,3);
t87 = pkin(4) * t118 + t136;
t81 = -pkin(4) * t111 + t131;
t71 = -pkin(4) * t103 + t125;
t67 = mrSges(5,1) * t109 - mrSges(5,2) * t108;
t52 = -Ifges(5,1) * t108 + Ifges(5,5) * t147 - t105;
t51 = -Ifges(5,2) * t109 + Ifges(5,6) * t147 - t203;
t48 = mrSges(6,1) * t142 + mrSges(6,3) * t57;
t47 = -mrSges(6,2) * t142 + mrSges(6,3) * t178;
t33 = -qJD(5) * t70 + t111 * t163 + t113 * t161;
t32 = qJD(5) * t69 + t111 * t161 - t113 * t163;
t29 = -mrSges(6,1) * t178 - mrSges(6,2) * t57;
t22 = t25 * mrSges(6,2);
t17 = pkin(7) * t111 + t31;
t16 = pkin(7) * t113 + t30;
t5 = -qJD(5) * t11 + t16 * t163 - t161 * t17;
t4 = qJD(5) * t10 + t16 * t161 + t163 * t17;
t8 = [(-t1 - t24 / 0.2e1 - t23 / 0.2e1 + t101 / 0.2e1 - t100 / 0.2e1 + t168) * t160 + t78 * (-mrSges(6,1) * t33 + mrSges(6,2) * t32) + t81 * t29 + t71 * (-mrSges(6,1) * t69 + mrSges(6,2) * t70) - t136 * t99 + t4 * t47 + t5 * t48 + t32 * t21 / 0.2e1 + t33 * t20 / 0.2e1 + (((mrSges(4,2) * t153 + (-t107 * mrSges(4,3) + Ifges(4,6) * t160 + (0.2e1 * qJ(2) * mrSges(4,1) - 0.3e1 / 0.2e1 * t204) * t158) * qJD(3)) * t164 + (mrSges(4,1) * t153 + (t106 * mrSges(4,3) + Ifges(4,5) * t160 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.3e1 / 0.2e1 * t205 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t164) * t158) * qJD(3)) * t162) * t158 + 0.2e1 * (t224 + (m(3) * t192 + m(4) * t155) * qJ(2)) * qJD(2)) * qJD(1) + (t175 * mrSges(4,3) + qJD(2) * t122 + t220 * qJD(3)) * t158 + (-t136 * mrSges(5,1) + t43 * mrSges(5,3) - Ifges(5,4) * t119 - Ifges(5,2) * t118 + Ifges(5,6) * t211) * t103 - (-t42 * mrSges(5,3) - Ifges(5,1) * t119 - Ifges(5,4) * t118 + Ifges(5,5) * t211) * t104 + t125 * (mrSges(5,1) * t118 - mrSges(5,2) * t119) + t178 * (Ifges(6,4) * t32 + Ifges(6,2) * t33) / 0.2e1 + m(4) * (t106 * t77 + t107 * t76 + t92 * t95 + t93 * t94) + (-t87 * mrSges(6,1) + t11 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,2) * t69 + Ifges(6,6) * t211) * t26 + (-t10 * mrSges(6,3) + Ifges(6,1) * t70 + Ifges(6,4) * t69 + Ifges(6,5) * t211) * t25 + (t2 * t69 - t3 * t70 - t6 * t32 + t7 * t33) * mrSges(6,3) + m(5) * (t125 * t136 + t127 * t131 + t14 * t42 + t15 * t43 + t30 * t35 + t31 * t36) + m(6) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + t71 * t87 + t78 * t81) + t87 * t22 + (t36 * t111 + t35 * t113 - t15 * t118 + t14 * t119) * mrSges(5,3) - t109 * (-Ifges(5,4) * t113 + Ifges(5,2) * t111) / 0.2e1 + t127 * (-mrSges(5,1) * t111 - mrSges(5,2) * t113) + (-Ifges(5,5) * t113 + Ifges(5,6) * t111) * t212 + (-Ifges(5,1) * t113 + Ifges(5,4) * t111) * t215 + t31 * t88 + t30 * t89 + t111 * t51 / 0.2e1 - t113 * t52 / 0.2e1 + t92 * t120 + t93 * t121 + t131 * t67 + t142 * (Ifges(6,5) * t32 + Ifges(6,6) * t33) / 0.2e1 + (Ifges(6,1) * t32 + Ifges(6,4) * t33) * t217; t223 * t89 + t222 * t88 + t231 * t48 + t230 * t47 + t173 * qJD(3) + (-t25 * t85 + t26 * t86) * mrSges(6,3) + (t103 * t135 - t104 * t172) * mrSges(5,3) - t165 * t224 + (-t173 * t160 + (-t122 - t29 - t67) * t158) * qJD(1) - m(3) * (t156 * t200 + t150) + (-t78 * t191 + t2 * t86 + t230 * t7 + t231 * t6 + t3 * t85) * m(6) + (-t127 * t191 + t15 * t135 - t14 * t172 + t222 * t36 + t223 * t35) * m(5) + (-t147 * t174 - t150 - t175) * m(4); -m(5) * (t35 * t38 + t36 * t39) + m(5) * (t14 * t159 + t15 * t157) * pkin(3) + ((t162 * (-Ifges(4,2) * t164 - t205) / 0.2e1 + (-Ifges(4,1) * t162 - t204) * t210 - qJ(2) * (mrSges(4,1) * t164 - mrSges(4,2) * t162)) * t191 + t176 * qJD(3) + (-m(5) * t127 - t67) * t164 * pkin(3) - t220) * t191 + (-t123 * t25 + t124 * t26 + t234) * mrSges(6,3) + t235 - t168 + t228 * t47 + t229 * t48 + (t123 * t3 + t124 * t2 + t228 * t7 + t229 * t6 - t78 * t91) * m(6) + (-t108 * t36 - t109 * t35 + (t103 * t157 + t104 * t159) * pkin(3)) * mrSges(5,3) + (Ifges(5,2) * t108 - t105 + t52) * t109 / 0.2e1 + t108 * (-Ifges(5,1) * t109 + t203) / 0.2e1 - t127 * (-mrSges(5,1) * t108 - mrSges(5,2) * t109) - t147 * (-Ifges(5,5) * t109 + Ifges(5,6) * t108) / 0.2e1 - t101 + t100 - t39 * t88 - t38 * t89 - t91 * t29 - t94 * t120 + t95 * t121 + t51 * t215; -t103 * mrSges(5,1) - t26 * mrSges(6,1) - t108 * t89 + t109 * t88 - t178 * t47 - t57 * t48 + t22 - t99 + (-t178 * t7 - t6 * t57 + t71) * m(6) + (-t35 * t108 + t36 * t109 + t125) * m(5); t234 * mrSges(6,3) - t6 * t47 + t7 * t48 - t216 + t235;];
tauc = t8(:);
