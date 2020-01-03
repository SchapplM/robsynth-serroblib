% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:19
% EndTime: 2019-12-31 19:07:26
% DurationCPUTime: 3.46s
% Computational Cost: add. (6482->353), mult. (17539->495), div. (0->0), fcn. (13412->8), ass. (0->166)
t147 = cos(qJ(5));
t142 = sin(pkin(9));
t146 = sin(qJ(3));
t143 = cos(pkin(9));
t149 = cos(qJ(3));
t190 = t143 * t149;
t128 = -t142 * t146 + t190;
t122 = t128 * qJD(3);
t114 = qJD(1) * t122;
t129 = t142 * t149 + t143 * t146;
t123 = t129 * qJD(3);
t115 = qJD(1) * t123;
t145 = sin(qJ(4));
t148 = cos(qJ(4));
t120 = t128 * qJD(1);
t121 = t129 * qJD(1);
t179 = t148 * t120 - t121 * t145;
t62 = t179 * qJD(4) + t114 * t148 - t115 * t145;
t141 = qJD(3) + qJD(4);
t144 = sin(qJ(5));
t165 = t120 * t145 + t148 * t121;
t82 = t141 * t147 - t144 * t165;
t38 = qJD(5) * t82 + t147 * t62;
t83 = t141 * t144 + t147 * t165;
t39 = -qJD(5) * t83 - t144 * t62;
t10 = -mrSges(6,1) * t39 + mrSges(6,2) * t38;
t203 = pkin(6) + qJ(2);
t133 = t203 * t142;
t130 = qJD(1) * t133;
t134 = t203 * t143;
t131 = qJD(1) * t134;
t102 = -t130 * t146 + t131 * t149;
t162 = t129 * qJD(2);
t76 = -qJD(1) * t162 - qJD(3) * t102;
t156 = -pkin(7) * t114 + t76;
t81 = pkin(7) * t120 + t102;
t192 = t148 * t81;
t101 = -t149 * t130 - t131 * t146;
t80 = -pkin(7) * t121 + t101;
t78 = qJD(3) * pkin(3) + t80;
t44 = t145 * t78 + t192;
t136 = qJD(2) * t190;
t184 = qJD(1) * qJD(2);
t187 = qJD(3) * t149;
t75 = -t130 * t187 + qJD(1) * t136 + (-qJD(3) * t131 - t142 * t184) * t146;
t69 = -pkin(7) * t115 + t75;
t14 = t44 * qJD(4) + t145 * t69 - t148 * t156;
t239 = -m(6) * t14 - t10;
t169 = Ifges(6,5) * t147 - Ifges(6,6) * t144;
t200 = Ifges(6,4) * t147;
t171 = -Ifges(6,2) * t144 + t200;
t201 = Ifges(6,4) * t144;
t173 = Ifges(6,1) * t147 - t201;
t174 = mrSges(6,1) * t144 + mrSges(6,2) * t147;
t218 = t147 / 0.2e1;
t219 = -t144 / 0.2e1;
t226 = t83 / 0.2e1;
t235 = qJD(5) - t179;
t216 = Ifges(6,4) * t83;
t34 = Ifges(6,2) * t82 + Ifges(6,6) * t235 + t216;
t79 = Ifges(6,4) * t82;
t35 = Ifges(6,1) * t83 + Ifges(6,5) * t235 + t79;
t195 = t145 * t81;
t43 = t148 * t78 - t195;
t41 = -pkin(4) * t141 - t43;
t238 = t34 * t219 + t35 * t218 + t41 * t174 + t82 * t171 / 0.2e1 + t173 * t226 + t235 * t169 / 0.2e1;
t42 = pkin(8) * t141 + t44;
t181 = -pkin(2) * t143 - pkin(1);
t132 = qJD(1) * t181 + qJD(2);
t103 = -pkin(3) * t120 + t132;
t45 = -pkin(4) * t179 - pkin(8) * t165 + t103;
t16 = t144 * t45 + t147 * t42;
t15 = -t144 * t42 + t147 * t45;
t194 = t147 * t15;
t167 = t144 * t16 + t194;
t237 = t167 * mrSges(6,3);
t105 = -t146 * t133 + t149 * t134;
t236 = -t144 * t15 + t147 * t16;
t13 = t43 * qJD(4) + t145 * t156 + t148 * t69;
t215 = pkin(3) * t115;
t63 = qJD(4) * t165 + t114 * t145 + t148 * t115;
t25 = pkin(4) * t63 - pkin(8) * t62 + t215;
t2 = qJD(5) * t15 + t13 * t147 + t144 * t25;
t3 = -qJD(5) * t16 - t13 * t144 + t147 * t25;
t234 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t38 + Ifges(6,6) * t39;
t91 = Ifges(5,4) * t179;
t233 = t165 * Ifges(5,1) / 0.2e1 + t91 / 0.2e1;
t232 = (m(3) * qJ(2) + mrSges(3,3)) * (t142 ^ 2 + t143 ^ 2);
t65 = pkin(4) * t165 - pkin(8) * t179;
t231 = t38 / 0.2e1;
t230 = t39 / 0.2e1;
t229 = t63 / 0.2e1;
t228 = -t82 / 0.2e1;
t227 = -t83 / 0.2e1;
t225 = -t235 / 0.2e1;
t222 = -t121 / 0.2e1;
t221 = t122 / 0.2e1;
t220 = -t123 / 0.2e1;
t214 = pkin(3) * t123;
t104 = -t149 * t133 - t134 * t146;
t89 = -pkin(7) * t129 + t104;
t90 = pkin(7) * t128 + t105;
t53 = t145 * t90 - t148 * t89;
t213 = t14 * t53;
t212 = t147 * t2;
t211 = t3 * t144;
t210 = t43 * mrSges(5,3);
t209 = t44 * mrSges(5,3);
t205 = t179 * Ifges(5,2);
t202 = -mrSges(5,1) * t141 - mrSges(6,1) * t82 + mrSges(6,2) * t83 + mrSges(5,3) * t165;
t199 = t121 * Ifges(4,4);
t186 = qJD(5) * t144;
t185 = qJD(5) * t147;
t180 = t63 * mrSges(5,1) + t62 * mrSges(5,2);
t176 = -t144 * t2 - t147 * t3;
t175 = mrSges(6,1) * t147 - mrSges(6,2) * t144;
t172 = Ifges(6,1) * t144 + t200;
t170 = Ifges(6,2) * t147 + t201;
t168 = Ifges(6,5) * t144 + Ifges(6,6) * t147;
t54 = t145 * t89 + t148 * t90;
t100 = t128 * t145 + t129 * t148;
t109 = -pkin(3) * t128 + t181;
t164 = t148 * t128 - t129 * t145;
t55 = -pkin(4) * t164 - pkin(8) * t100 + t109;
t27 = t144 * t55 + t147 * t54;
t26 = -t144 * t54 + t147 * t55;
t86 = -t133 * t187 + t136 + (-qJD(2) * t142 - qJD(3) * t134) * t146;
t87 = -t105 * qJD(3) - t162;
t8 = t38 * Ifges(6,4) + t39 * Ifges(6,2) + t63 * Ifges(6,6);
t9 = t38 * Ifges(6,1) + t39 * Ifges(6,4) + t63 * Ifges(6,5);
t158 = -t13 * mrSges(5,2) + mrSges(6,3) * t212 + t144 * t9 / 0.2e1 + t172 * t231 + t170 * t230 + t168 * t229 - Ifges(5,6) * t63 + Ifges(5,5) * t62 + t8 * t218 + (-t175 - mrSges(5,1)) * t14 + t238 * qJD(5);
t157 = -pkin(7) * t122 + t87;
t155 = t16 * mrSges(6,2) - t235 * Ifges(6,3) - t83 * Ifges(6,5) - t82 * Ifges(6,6) + t141 * Ifges(5,6) + t205 / 0.2e1 + Ifges(5,4) * t165 - t103 * mrSges(5,1) - t15 * mrSges(6,1);
t17 = mrSges(6,1) * t63 - mrSges(6,3) * t38;
t18 = -mrSges(6,2) * t63 + mrSges(6,3) * t39;
t51 = -mrSges(6,2) * t235 + mrSges(6,3) * t82;
t52 = mrSges(6,1) * t235 - mrSges(6,3) * t83;
t154 = -t144 * t17 - t51 * t186 - t52 * t185 + t147 * t18 + m(6) * (-t15 * t185 - t16 * t186 - t211 + t212);
t153 = -t205 / 0.2e1 - t155;
t152 = t103 * mrSges(5,2) + t141 * Ifges(5,5) + t233 + t238;
t151 = t152 + t233;
t116 = Ifges(4,4) * t120;
t110 = t114 * mrSges(4,2);
t108 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t121;
t107 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t120;
t94 = t121 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t116;
t93 = t120 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t199;
t84 = -mrSges(5,2) * t141 + mrSges(5,3) * t179;
t73 = -pkin(7) * t123 + t86;
t68 = qJD(4) * t100 + t122 * t145 + t148 * t123;
t67 = qJD(4) * t164 + t122 * t148 - t123 * t145;
t64 = -mrSges(5,1) * t179 + mrSges(5,2) * t165;
t59 = Ifges(6,3) * t63;
t49 = pkin(3) * t121 + t65;
t47 = t148 * t80 - t195;
t46 = t145 * t80 + t192;
t28 = pkin(4) * t68 - pkin(8) * t67 + t214;
t24 = t144 * t65 + t147 * t43;
t23 = -t144 * t43 + t147 * t65;
t22 = t144 * t49 + t147 * t47;
t21 = -t144 * t47 + t147 * t49;
t20 = qJD(4) * t54 + t145 * t73 - t148 * t157;
t19 = -qJD(4) * t53 + t145 * t157 + t148 * t73;
t5 = -qJD(5) * t27 - t144 * t19 + t147 * t28;
t4 = qJD(5) * t26 + t144 * t28 + t147 * t19;
t1 = [(t151 - t237) * t67 + (t169 * t229 + t173 * t231 + t171 * t230 + Ifges(5,1) * t62 - Ifges(5,4) * t63 + mrSges(5,2) * t215 + t8 * t219 + t9 * t218 + (mrSges(5,3) + t174) * t14 + t176 * mrSges(6,3) + (t168 * t225 + t41 * t175 + t170 * t228 + t172 * t227 + t35 * t219 - t147 * t34 / 0.2e1 - t236 * mrSges(6,3)) * qJD(5)) * t100 + m(4) * (t101 * t87 + t102 * t86 + t104 * t76 + t105 * t75) + t132 * (mrSges(4,1) * t123 + mrSges(4,2) * t122) + qJD(3) * (Ifges(4,5) * t122 - Ifges(4,6) * t123) / 0.2e1 + t86 * t107 + t87 * t108 + (-t43 * t67 - t44 * t68 + t53 * t62 - t54 * t63) * mrSges(5,3) + t4 * t51 + t5 * t52 + t53 * t10 + t26 * t17 + t27 * t18 + t109 * t180 + (-t101 * t122 - t102 * t123 - t104 * t114 - t105 * t115 + t128 * t75 - t129 * t76) * mrSges(4,3) + m(5) * (t13 * t54 + t213 + t19 * t44 - t20 * t43 + (t103 * t123 + t109 * t115) * pkin(3)) + t181 * (t115 * mrSges(4,1) + t110) + (-t115 * t128 + t120 * t220) * Ifges(4,2) + (t114 * t128 - t115 * t129 + t120 * t221 + t121 * t220) * Ifges(4,4) - (t59 / 0.2e1 - Ifges(5,4) * t62 - t13 * mrSges(5,3) + mrSges(5,1) * t215 + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t63 + t234) * t164 + t19 * t84 + t202 * t20 + m(6) * (t15 * t5 + t16 * t4 + t2 * t27 + t20 * t41 + t26 * t3 + t213) + t64 * t214 + t93 * t220 + t94 * t221 + (t114 * t129 + t121 * t221) * Ifges(4,1) + t153 * t68 + 0.2e1 * t232 * t184; -t120 * t107 + t121 * t108 - t179 * t84 + t110 - t202 * t165 - (-m(5) * pkin(3) - mrSges(4,1)) * t115 + (t235 * t51 + t17) * t147 + (-t235 * t52 + t18) * t144 - m(4) * (-t101 * t121 + t102 * t120) - m(5) * (-t165 * t43 + t179 * t44) + t180 - t232 * qJD(1) ^ 2 + (-t41 * t165 + t235 * t236 - t176) * m(6); t158 + (t101 * t120 + t102 * t121) * mrSges(4,3) + (-t121 * t64 + (-t145 * t63 - t148 * t62) * mrSges(5,3) + (t202 * t145 + (-t144 * t52 + t147 * t51 + t84) * t148 + m(6) * (t145 * t41 + t236 * t148)) * qJD(4) + (0.2e1 * t103 * t222 + t13 * t145 - t14 * t148 + (-t145 * t43 + t148 * t44) * qJD(4)) * m(5)) * pkin(3) + (-t153 + t209) * t165 + (-t151 + t210) * t179 - m(6) * (t15 * t21 + t16 * t22 + t41 * t46) - m(5) * (-t43 * t46 + t44 * t47) - t132 * (mrSges(4,1) * t121 + mrSges(4,2) * t120) - qJD(3) * (Ifges(4,5) * t120 - Ifges(4,6) * t121) / 0.2e1 + t121 * t93 / 0.2e1 - t101 * t107 + t102 * t108 + Ifges(4,5) * t114 - Ifges(4,6) * t115 - t22 * t51 - t21 * t52 - (-Ifges(4,2) * t121 + t116 + t94) * t120 / 0.2e1 + (-t235 * t194 + (-t16 * t235 - t3) * t144) * mrSges(6,3) - t75 * mrSges(4,2) + t76 * mrSges(4,1) - t47 * t84 - t202 * t46 + (Ifges(4,1) * t120 - t199) * t222 + t154 * (pkin(3) * t145 + pkin(8)) - t239 * (-pkin(3) * t148 - pkin(4)); t158 + (-qJD(5) * t167 - t211) * mrSges(6,3) - t202 * t44 - m(6) * (t15 * t23 + t16 * t24 + t41 * t44) - t24 * t51 - t23 * t52 - t43 * t84 + (-t91 / 0.2e1 + t210 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t165 + t237 - t152) * t179 + t154 * pkin(8) + (t155 + t209) * t165 + t239 * pkin(4); t59 - t41 * (mrSges(6,1) * t83 + mrSges(6,2) * t82) + (Ifges(6,1) * t82 - t216) * t227 + t34 * t226 + (Ifges(6,5) * t82 - Ifges(6,6) * t83) * t225 - t15 * t51 + t16 * t52 + (t15 * t82 + t16 * t83) * mrSges(6,3) + (-Ifges(6,2) * t83 + t35 + t79) * t228 + t234;];
tauc = t1(:);
