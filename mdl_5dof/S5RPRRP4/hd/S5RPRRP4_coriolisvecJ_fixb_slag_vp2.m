% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:16
% EndTime: 2020-01-03 11:49:31
% DurationCPUTime: 4.88s
% Computational Cost: add. (2894->302), mult. (7858->419), div. (0->0), fcn. (5288->6), ass. (0->157)
t196 = Ifges(5,4) + Ifges(6,4);
t222 = Ifges(5,1) + Ifges(6,1);
t221 = Ifges(5,5) + Ifges(6,5);
t220 = Ifges(5,2) + Ifges(6,2);
t219 = Ifges(5,6) + Ifges(6,6);
t142 = sin(pkin(8));
t144 = sin(qJ(4));
t145 = sin(qJ(3));
t146 = cos(qJ(4));
t147 = cos(qJ(3));
t121 = t144 * t147 + t145 * t146;
t153 = qJD(1) * t121;
t102 = t142 * t153;
t224 = t196 * t102;
t176 = qJD(1) * t142;
t165 = t147 * t176;
t166 = t145 * t176;
t104 = -t144 * t166 + t146 * t165;
t223 = t196 * t104;
t143 = cos(pkin(8));
t175 = qJD(1) * t143;
t132 = qJD(3) - t175;
t126 = qJD(4) + t132;
t218 = -t220 * t102 + t219 * t126 + t223;
t217 = t222 * t104 + t221 * t126 - t224;
t138 = t142 * qJD(2);
t173 = qJD(3) * t147;
t116 = t142 * pkin(3) * t173 + t138;
t122 = -pkin(2) * t143 - pkin(6) * t142 - pkin(1);
t113 = qJD(1) * t122 + qJD(2);
t106 = t147 * t113;
t185 = qJ(2) * t145;
t168 = t143 * t185;
t199 = pkin(7) * t142;
t169 = t147 * t199;
t152 = -t168 - t169;
t75 = qJD(1) * t152 + t106;
t56 = pkin(3) * t132 + t75;
t167 = qJ(2) * t175;
t182 = t113 * t145;
t92 = t147 * t167 + t182;
t76 = -pkin(7) * t166 + t92;
t68 = t144 * t76;
t19 = t146 * t56 - t68;
t96 = t104 * qJ(5);
t13 = t19 - t96;
t140 = t142 ^ 2;
t141 = t143 ^ 2;
t177 = t140 + t141;
t216 = t177 * mrSges(3,3);
t212 = qJD(3) + qJD(4);
t90 = t212 * t121;
t215 = t143 * t153 - t90;
t154 = t144 * t145 - t146 * t147;
t214 = (t175 - t212) * t154;
t108 = t154 * t142;
t180 = t143 * t147;
t131 = qJ(2) * t180;
t100 = t145 * t122 + t131;
t161 = qJD(3) * t168;
t130 = qJD(2) * t180;
t179 = qJD(1) * t130 + t113 * t173;
t66 = -qJD(1) * t161 + t179;
t163 = qJD(3) * t182;
t174 = qJD(2) * t145;
t67 = -t163 + (-qJ(2) * t173 - t174) * t175;
t213 = -t67 * mrSges(4,1) + t66 * mrSges(4,2);
t65 = t212 * t108;
t55 = qJD(1) * t65;
t171 = qJD(4) * t146;
t172 = qJD(4) * t144;
t151 = t152 * qJD(3);
t46 = qJD(1) * t151 + t179;
t164 = t143 * t174;
t181 = t142 * t145;
t170 = pkin(7) * t181;
t47 = -t163 + (-t164 + (-t131 + t170) * qJD(3)) * qJD(1);
t7 = t144 * t47 + t146 * t46 + t56 * t171 - t172 * t76;
t2 = qJ(5) * t55 - qJD(5) * t102 + t7;
t64 = t90 * t142;
t54 = qJD(1) * t64;
t70 = t146 * t76;
t20 = t144 * t56 + t70;
t8 = -qJD(4) * t20 - t144 * t46 + t146 * t47;
t3 = qJ(5) * t54 - qJD(5) * t104 + t8;
t211 = -t8 * mrSges(5,1) - t3 * mrSges(6,1) + t7 * mrSges(5,2) + t2 * mrSges(6,2);
t91 = -t145 * t167 + t106;
t156 = t91 * t145 - t92 * t147;
t159 = -Ifges(4,5) * t145 - Ifges(4,6) * t147;
t190 = Ifges(4,4) * t147;
t191 = Ifges(4,4) * t145;
t201 = -t147 / 0.2e1;
t210 = t156 * mrSges(4,3) + t132 * t159 / 0.2e1 + (t132 * Ifges(4,6) + (-t145 * Ifges(4,2) + t190) * t176) * t201 - t145 * (t132 * Ifges(4,5) + (t147 * Ifges(4,1) - t191) * t176) / 0.2e1;
t207 = t102 / 0.2e1;
t204 = t104 / 0.2e1;
t200 = pkin(4) * t104;
t23 = t146 * t75 - t68;
t192 = mrSges(6,3) * t102;
t79 = -mrSges(6,2) * t126 - t192;
t193 = mrSges(5,3) * t102;
t80 = -mrSges(5,2) * t126 - t193;
t195 = t79 + t80;
t81 = mrSges(6,1) * t126 - t104 * mrSges(6,3);
t189 = t104 * mrSges(5,3);
t82 = mrSges(5,1) * t126 - t189;
t194 = t81 + t82;
t118 = t147 * t122;
t83 = -t169 + t118 + (-pkin(3) - t185) * t143;
t93 = -t170 + t100;
t29 = t144 * t83 + t146 * t93;
t186 = t144 * t55;
t148 = qJD(1) ^ 2;
t184 = qJ(2) * t148;
t183 = qJ(5) * t102;
t178 = t122 * t173 + t130;
t112 = t116 * qJD(1);
t114 = pkin(3) * t166 + qJ(2) * t176;
t119 = pkin(3) * t181 + t142 * qJ(2);
t162 = -t55 * mrSges(6,1) - t54 * mrSges(6,2);
t22 = -t144 * t75 - t70;
t28 = -t144 * t93 + t146 * t83;
t157 = -t66 * t145 - t67 * t147;
t109 = -mrSges(4,2) * t132 - mrSges(4,3) * t166;
t110 = mrSges(4,1) * t132 - mrSges(4,3) * t165;
t155 = t147 * t109 - t145 * t110;
t71 = t151 + t178;
t72 = -t164 + (-t131 + (-t122 + t199) * t145) * qJD(3);
t9 = t144 * t72 + t146 * t71 + t83 * t171 - t172 * t93;
t10 = -qJD(4) * t29 - t144 * t71 + t146 * t72;
t12 = pkin(4) * t126 + t13;
t50 = Ifges(6,6) * t55;
t51 = Ifges(5,6) * t55;
t52 = Ifges(6,5) * t54;
t53 = Ifges(5,5) * t54;
t73 = pkin(4) * t102 + qJD(5) + t114;
t149 = -t53 - t52 + t51 + t50 - t12 * t192 - t19 * t193 - t73 * (mrSges(6,1) * t104 - mrSges(6,2) * t102) - t114 * (mrSges(5,1) * t104 - mrSges(5,2) * t102) - t211 - (-t222 * t102 - t223) * t104 / 0.2e1 + t218 * t204 - (-t221 * t102 - t219 * t104) * t126 / 0.2e1 + (-t220 * t104 + t217 - t224) * t207;
t137 = pkin(3) * t146 + pkin(4);
t135 = t140 * t184;
t111 = (mrSges(4,1) * t145 + mrSges(4,2) * t147) * t176;
t107 = t121 * t142;
t99 = t118 - t168;
t88 = pkin(3) * t165 + t200;
t86 = -t100 * qJD(3) - t164;
t85 = -t161 + t178;
t84 = pkin(4) * t107 + t119;
t60 = mrSges(5,1) * t102 + mrSges(5,2) * t104;
t59 = mrSges(6,1) * t102 + mrSges(6,2) * t104;
t35 = -pkin(4) * t65 + t116;
t34 = -pkin(4) * t55 + t112;
t21 = -qJ(5) * t107 + t29;
t18 = -pkin(4) * t143 + qJ(5) * t108 + t28;
t17 = -t96 + t23;
t16 = t22 + t183;
t14 = t20 - t183;
t5 = qJ(5) * t64 + qJD(5) * t108 + t10;
t4 = qJ(5) * t65 - qJD(5) * t107 + t9;
t1 = [t114 * (-mrSges(5,1) * t65 - mrSges(5,2) * t64) + t73 * (-mrSges(6,1) * t65 - mrSges(6,2) * t64) + (t219 * t65 - t221 * t64) * t126 / 0.2e1 - (-t196 * t64 + t220 * t65) * t102 / 0.2e1 + (t196 * t65 - t222 * t64) * t204 + t85 * t109 + t86 * t110 + t116 * t60 + t4 * t79 + t9 * t80 + t5 * t81 + t10 * t82 + t35 * t59 + t84 * t162 + m(5) * (t10 * t19 + t112 * t119 + t114 * t116 + t20 * t9 + t28 * t8 + t29 * t7) + m(6) * (t12 * t5 + t14 * t4 + t18 * t3 + t2 * t21 + t34 * t84 + t35 * t73) + (t12 * t64 + t14 * t65 + t18 * t54 + t21 * t55) * mrSges(6,3) + (t19 * t64 + t20 * t65 + t28 * t54 + t29 * t55) * mrSges(5,3) + t119 * (-mrSges(5,1) * t55 - mrSges(5,2) * t54) - (-mrSges(5,1) * t112 - mrSges(6,1) * t34 + mrSges(5,3) * t7 + mrSges(6,3) * t2 - t196 * t54 + t220 * t55) * t107 - (mrSges(5,2) * t112 + mrSges(6,2) * t34 - mrSges(5,3) * t8 - mrSges(6,3) * t3 + t196 * t55 - t222 * t54) * t108 + (((mrSges(4,2) * t138 + (-t100 * mrSges(4,3) + Ifges(4,6) * t143 + (0.2e1 * qJ(2) * mrSges(4,1) - 0.3e1 / 0.2e1 * t190) * t142) * qJD(3)) * t147 + (mrSges(4,1) * t138 + (t99 * mrSges(4,3) + Ifges(4,5) * t143 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.3e1 / 0.2e1 * t191 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t147) * t142) * qJD(3)) * t145) * t142 + 0.2e1 * (t216 + (m(3) * t177 + m(4) * t140) * qJ(2)) * qJD(2)) * qJD(1) + (t157 * mrSges(4,3) + qJD(2) * t111 + t210 * qJD(3)) * t142 + (t52 / 0.2e1 - t50 / 0.2e1 + t53 / 0.2e1 - t51 / 0.2e1 + (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t55 - (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * t54 + t211 + t213) * t143 + m(4) * (t100 * t66 + t67 * t99 + t85 * t92 + t86 * t91) - t217 * t64 / 0.2e1 + t218 * t65 / 0.2e1; t155 * qJD(3) - t148 * t216 + (-t155 * t143 + (-t111 - t59 - t60) * t142) * qJD(1) - m(3) * (t141 * t184 + t135) + t214 * t195 + t215 * t194 + (mrSges(6,3) + mrSges(5,3)) * (t121 * t55 - t154 * t54) + (t12 * t215 + t121 * t2 + t14 * t214 - t154 * t3 - t176 * t73) * m(6) + (-t114 * t176 + t121 * t7 - t154 * t8 + t19 * t215 + t20 * t214) * m(5) + (-t132 * t156 - t135 - t157) * m(4); ((-qJ(2) * (mrSges(4,1) * t147 - mrSges(4,2) * t145) + (-Ifges(4,1) * t145 - t190) * t201 + t145 * (-Ifges(4,2) * t147 - t191) / 0.2e1) * t176 + t159 * qJD(3) + (-m(5) * t114 - t60) * t147 * pkin(3) - t210) * t176 + (mrSges(6,3) * t186 + (t146 * t54 + t186) * mrSges(5,3) + (-t144 * t194 + t146 * t195) * qJD(4) + m(5) * (t144 * t7 + t146 * t8 + t171 * t20 - t172 * t19)) * pkin(3) - t91 * t109 + t92 * t110 + (t104 * t14 + t137 * t54) * mrSges(6,3) - t17 * t79 - t23 * t80 - t16 * t81 - t22 * t82 - t88 * t59 - m(5) * (t19 * t22 + t20 * t23) + t20 * t189 + t149 + ((-t12 * t172 + t14 * t171 + t144 * t2) * pkin(3) + t137 * t3 - t12 * t16 - t14 * t17 - t73 * t88) * m(6) - t213; (t20 * mrSges(5,3) + t14 * mrSges(6,3) - pkin(4) * t59) * t104 + (-t73 * t200 - (-t12 + t13) * t14 + t3 * pkin(4)) * m(6) - t13 * t79 - t19 * t80 + t14 * t81 + t20 * t82 + pkin(4) * t54 * mrSges(6,3) + t149; t102 * t79 + t104 * t81 + 0.2e1 * (t34 / 0.2e1 + t14 * t207 + t12 * t204) * m(6) + t162;];
tauc = t1(:);
