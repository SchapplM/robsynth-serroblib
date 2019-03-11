% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:07
% EndTime: 2019-03-09 02:44:12
% DurationCPUTime: 2.81s
% Computational Cost: add. (2221->370), mult. (4949->475), div. (0->0), fcn. (2365->6), ass. (0->170)
t117 = sin(qJ(6));
t119 = cos(qJ(6));
t139 = mrSges(7,1) * t117 + mrSges(7,2) * t119;
t118 = sin(qJ(3));
t120 = cos(qJ(3));
t145 = pkin(5) * t118 + pkin(8) * t120;
t167 = qJD(1) * t120;
t168 = qJD(1) * t118;
t185 = -cos(pkin(9)) * pkin(1) - pkin(2);
t77 = t185 * qJD(1);
t47 = -pkin(3) * t167 - qJ(4) * t168 + t77;
t29 = pkin(4) * t167 + qJD(5) - t47;
t15 = qJD(1) * t145 + t29;
t121 = -pkin(3) - pkin(4);
t109 = -pkin(8) + t121;
t102 = t120 * qJD(2);
t87 = sin(pkin(9)) * pkin(1) + pkin(7);
t73 = t87 * qJD(1);
t169 = t118 * t73 - t102;
t157 = qJD(4) + t169;
t90 = qJ(5) * t168;
t151 = -t90 + t157;
t20 = qJD(3) * t109 + t151;
t5 = -t117 * t20 + t119 * t15;
t6 = t117 * t15 + t119 * t20;
t142 = t6 * t117 + t5 * t119;
t69 = qJD(3) * t117 + t119 * t167;
t195 = Ifges(7,4) * t69;
t165 = qJD(3) * t119;
t68 = t117 * t167 - t165;
t86 = qJD(6) + t168;
t17 = Ifges(7,2) * t68 + Ifges(7,6) * t86 - t195;
t65 = Ifges(7,4) * t68;
t18 = -Ifges(7,1) * t69 + Ifges(7,5) * t86 + t65;
t196 = -t119 / 0.2e1;
t197 = t117 / 0.2e1;
t110 = qJD(3) * qJ(4);
t101 = t118 * qJD(2);
t53 = t120 * t73 + t101;
t39 = -qJ(5) * t167 + t53;
t31 = -t110 - t39;
t26 = qJD(3) * pkin(5) - t31;
t224 = -t142 * mrSges(7,3) + t26 * t139 - t17 * t197 - t18 * t196;
t223 = -t167 / 0.2e1;
t222 = -mrSges(5,1) - mrSges(4,1);
t221 = mrSges(5,2) + mrSges(4,3);
t175 = Ifges(7,6) * t117;
t176 = Ifges(7,5) * t119;
t134 = -t175 + t176;
t178 = Ifges(7,4) * t119;
t136 = -Ifges(7,2) * t117 + t178;
t179 = Ifges(7,4) * t117;
t138 = Ifges(7,1) * t119 - t179;
t198 = t86 / 0.2e1;
t201 = -t69 / 0.2e1;
t202 = t68 / 0.2e1;
t212 = -qJD(3) / 0.2e1;
t213 = qJD(1) / 0.2e1;
t214 = -qJD(1) / 0.2e1;
t218 = Ifges(5,6) / 0.2e1;
t46 = t110 + t53;
t96 = Ifges(5,5) * t168;
t220 = t29 * mrSges(6,2) + t47 * mrSges(5,1) + t77 * mrSges(4,1) + qJD(3) * t218 + Ifges(5,3) * t223 + t96 / 0.2e1 + (Ifges(4,4) * t118 + t120 * Ifges(4,2)) * t214 + (-t120 * Ifges(6,1) - Ifges(6,4) * t118) * t213 - t31 * mrSges(6,3) - t46 * mrSges(5,2) - t53 * mrSges(4,3) + t136 * t202 + t138 * t201 + t134 * t198 + (Ifges(4,6) + Ifges(6,5)) * t212 + t224;
t219 = -Ifges(4,1) / 0.2e1;
t163 = qJD(6) * t120;
t155 = t117 * t163;
t160 = qJD(3) * qJD(6);
t36 = -t119 * t160 + (t118 * t165 + t155) * qJD(1);
t217 = -t36 / 0.2e1;
t154 = t119 * t163;
t166 = qJD(3) * t118;
t37 = t117 * t160 + (-t117 * t166 + t154) * qJD(1);
t216 = -t37 / 0.2e1;
t215 = Ifges(4,4) * t223;
t171 = -qJ(5) + t87;
t126 = pkin(5) * t120 + t109 * t118;
t125 = t126 * qJD(3);
t100 = t118 * qJD(4);
t162 = qJD(1) * qJD(3);
t152 = t120 * t162;
t181 = qJ(4) * t152 + qJD(1) * t100;
t14 = qJD(1) * t125 + t181;
t161 = qJD(1) * qJD(5);
t21 = -t118 * t161 + (t101 + (-qJ(5) * qJD(1) + t73) * t120) * qJD(3);
t1 = qJD(6) * t5 + t117 * t14 + t119 * t21;
t2 = -qJD(6) * t6 - t117 * t21 + t119 * t14;
t144 = -t1 * t119 + t2 * t117;
t141 = -t117 * t5 + t119 * t6;
t38 = -t90 + t169;
t211 = t38 + qJD(4);
t210 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t36 + Ifges(7,6) * t37;
t40 = -mrSges(7,2) * t86 + mrSges(7,3) * t68;
t41 = mrSges(7,1) * t86 + mrSges(7,3) * t69;
t131 = t117 * t40 + t119 * t41;
t22 = mrSges(7,1) * t152 - mrSges(7,3) * t36;
t23 = -mrSges(7,2) * t152 + mrSges(7,3) * t37;
t133 = t117 * t22 - t119 * t23;
t208 = t131 * qJD(6) + t133;
t156 = t121 * qJD(3);
t25 = t156 + t151;
t207 = m(6) * (t118 * t25 - t120 * t31) - m(7) * (-t141 * t118 - t120 * t26);
t206 = 0.2e1 * t87;
t205 = m(4) / 0.2e1;
t204 = m(5) / 0.2e1;
t203 = -t68 / 0.2e1;
t200 = t69 / 0.2e1;
t199 = -t86 / 0.2e1;
t153 = t118 * t162;
t49 = qJD(3) * t102 - t166 * t73;
t35 = qJD(3) * qJD(4) + t49;
t19 = -qJ(5) * t153 + t120 * t161 - t35;
t64 = t171 * t120;
t191 = t19 * t64;
t186 = mrSges(5,2) - mrSges(6,3);
t184 = qJD(3) * t222 + t168 * t221;
t183 = -qJD(3) * mrSges(6,1) + mrSges(7,1) * t68 + mrSges(7,2) * t69 + mrSges(6,3) * t167;
t79 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t167;
t80 = mrSges(5,2) * t167 + qJD(3) * mrSges(5,3);
t182 = t79 + t80;
t180 = mrSges(6,1) * t152 + mrSges(6,2) * t153;
t177 = Ifges(7,5) * t117;
t174 = Ifges(7,6) * t119;
t173 = t118 * t19;
t50 = t53 * qJD(3);
t172 = t120 * t50;
t164 = qJD(3) * t120;
t170 = qJ(4) * t164 + t100;
t158 = t80 - t183;
t62 = -t120 * pkin(3) - t118 * qJ(4) + t185;
t150 = t118 * t156;
t149 = -0.3e1 / 0.2e1 * Ifges(6,4) + 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4);
t148 = -Ifges(6,5) / 0.2e1 + t218 - Ifges(4,6) / 0.2e1;
t147 = Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t143 = t1 * t117 + t119 * t2;
t51 = t120 * pkin(4) - t62;
t140 = -mrSges(7,1) * t119 + mrSges(7,2) * t117;
t137 = Ifges(7,1) * t117 + t178;
t135 = Ifges(7,2) * t119 + t179;
t32 = t145 + t51;
t63 = t171 * t118;
t12 = -t117 * t63 + t119 * t32;
t13 = t117 * t32 + t119 * t63;
t132 = -t117 * t41 + t119 * t40;
t74 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t168;
t129 = t132 + t74;
t42 = -qJD(3) * pkin(3) + t157;
t123 = t25 * mrSges(6,3) + t47 * mrSges(5,3) - t169 * mrSges(4,3) + t6 * mrSges(7,2) - t86 * Ifges(7,3) + t69 * Ifges(7,5) - t68 * Ifges(7,6) + (-Ifges(6,4) * t120 - Ifges(6,2) * t118) * t213 + (Ifges(5,1) * t118 - Ifges(5,5) * t120) * t214 + t168 * t219 + t215 - t29 * mrSges(6,1) - t42 * mrSges(5,2) - t5 * mrSges(7,1) - t77 * mrSges(4,2) + (Ifges(6,6) + Ifges(5,4) + Ifges(4,5)) * t212;
t116 = qJ(4) + pkin(5);
t93 = qJ(4) * t167;
t85 = Ifges(7,3) * t152;
t72 = (mrSges(6,1) * t118 - mrSges(6,2) * t120) * qJD(1);
t71 = pkin(3) * t168 - t93;
t70 = (-mrSges(5,1) * t120 - mrSges(5,3) * t118) * qJD(1);
t55 = pkin(3) * t166 - t170;
t54 = t121 * t168 + t93;
t48 = pkin(3) * t153 - t181;
t45 = -qJD(5) * t118 + t164 * t171;
t44 = qJD(5) * t120 + t166 * t171;
t43 = t150 + t170;
t30 = qJD(1) * t150 + t181;
t28 = qJD(1) * t126 + t93;
t24 = t125 + t170;
t11 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t10 = t36 * Ifges(7,1) + t37 * Ifges(7,4) + Ifges(7,5) * t152;
t9 = t36 * Ifges(7,4) + t37 * Ifges(7,2) + Ifges(7,6) * t152;
t8 = t117 * t28 + t119 * t39;
t7 = -t117 * t39 + t119 * t28;
t4 = -qJD(6) * t13 - t117 * t45 + t119 * t24;
t3 = qJD(6) * t12 + t117 * t24 + t119 * t45;
t16 = [t51 * t180 + t45 * t74 + t55 * t70 + t43 * t72 + t64 * t11 + t3 * t40 + t4 * t41 + t13 * t23 + t12 * t22 + t183 * t44 + m(5) * (t47 * t55 + t48 * t62) + m(7) * (t1 * t13 + t12 * t2 - t26 * t44 + t3 * t6 + t4 * t5 - t191) + m(6) * (t21 * t63 + t25 * t45 + t29 * t43 + t30 * t51 + t31 * t44 - t191) + (t85 / 0.2e1 - t21 * mrSges(6,3) - t48 * mrSges(5,3) + t30 * mrSges(6,1) + ((t204 + t205) * t206 + t221) * t50 + (t148 * qJD(3) + (-m(4) * t53 - m(5) * t46 - t182) * t87 + (mrSges(4,1) * t185 + t62 * mrSges(5,1) + t64 * mrSges(6,3) + t118 * t149) * qJD(1) + t220) * qJD(3) + t210) * t118 + (t9 * t197 + t49 * mrSges(4,3) + t35 * mrSges(5,2) - t48 * mrSges(5,1) - t30 * mrSges(6,2) + t138 * t217 + t136 * t216 + t10 * t196 + (mrSges(6,3) + t139) * t19 + t143 * mrSges(7,3) + (t204 * t35 + t205 * t49) * t206 + (t18 * t197 + t119 * t17 / 0.2e1 + (t174 + t177) * t198 + t135 * t202 + t137 * t201 + t26 * t140 + t141 * mrSges(7,3)) * qJD(6) + ((m(4) * t169 + m(5) * t42 + t184) * t87 + t147 * qJD(3) - t123 + (t185 * mrSges(4,2) - t62 * mrSges(5,3) - t63 * mrSges(6,3) + (-t176 / 0.2e1 + t175 / 0.2e1 - t149) * t120 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(5,1)) * t118) * qJD(1)) * qJD(3)) * t120; t118 * t11 + t208 * t120 + m(4) * (t118 * t49 - t172) + m(5) * (t118 * t35 - t172) + m(6) * (-t120 * t21 - t173) + m(7) * (t120 * t144 + t154 * t5 + t155 * t6 - t173) + ((t79 + t158) * t120 + (t129 + t184) * t118 + m(4) * (t118 * t169 + t120 * t53) + m(5) * (t118 * t42 + t120 * t46) + (t118 ^ 2 + t120 ^ 2) * qJD(1) * (-mrSges(4,3) - t186) + t207) * qJD(3); t182 * t169 - t183 * t38 - t184 * t53 + (t134 * t199 + t136 * t203 + t138 * t200 - t224) * qJD(6) - t117 * t10 / 0.2e1 + t116 * t11 - t39 * t74 - t71 * t70 - t54 * t72 - t49 * mrSges(4,2) + t35 * mrSges(5,3) - t8 * t40 - t7 * t41 + t21 * mrSges(6,2) + t158 * qJD(4) + t9 * t196 + t144 * mrSges(7,3) + (-mrSges(6,1) + t140) * t19 + (-t116 * t19 + t211 * t26 - t5 * t7 - t6 * t8) * m(7) + (-t19 * qJ(4) + t121 * t21 - t211 * t31 - t25 * t39 - t29 * t54) * m(6) + ((-m(7) * t142 - t131) * qJD(6) - m(7) * t144 - t133) * t109 + t135 * t216 + t137 * t217 + t222 * t50 + (((-pkin(3) * mrSges(5,2) - t121 * mrSges(6,3) - t177 / 0.2e1 - t174 / 0.2e1 + t147) * qJD(3) + t215 + (-Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t167 + t123) * t120 + (-t96 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t219) * t167 + (-t186 * qJ(4) + t148) * qJD(3) + (Ifges(6,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t168 - t220) * t118) * qJD(1) + (-pkin(3) * t50 + qJ(4) * t35 + t157 * t46 - t42 * t53 - t47 * t71) * m(5); -t158 * qJD(3) + (t186 * t164 + (-t131 + t70 - t72) * t118) * qJD(1) + (-qJD(3) * t26 - t142 * t86 - t144) * m(7) + (qJD(3) * t31 - t168 * t29 + t21) * m(6) + (-qJD(3) * t46 + t168 * t47 + t50) * m(5) - t208; t117 * t23 + t119 * t22 + t132 * qJD(6) + m(6) * t30 + m(7) * (qJD(6) * t141 + t143) + (t118 * t129 - t120 * t183 + t207) * qJD(1) + t180; t85 - t26 * (-mrSges(7,1) * t69 + mrSges(7,2) * t68) + (Ifges(7,1) * t68 + t195) * t200 + t17 * t201 + (Ifges(7,5) * t68 + Ifges(7,6) * t69) * t199 - t5 * t40 + t6 * t41 + (t5 * t68 - t6 * t69) * mrSges(7,3) + (Ifges(7,2) * t69 + t18 + t65) * t203 + t210;];
tauc  = t16(:);
