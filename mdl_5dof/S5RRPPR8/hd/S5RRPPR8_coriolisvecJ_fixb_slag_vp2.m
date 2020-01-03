% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:38:17
% DurationCPUTime: 3.72s
% Computational Cost: add. (2508->344), mult. (6210->471), div. (0->0), fcn. (3935->6), ass. (0->158)
t231 = -mrSges(3,1) - mrSges(4,1);
t156 = sin(qJ(5));
t158 = cos(qJ(5));
t153 = sin(pkin(8));
t154 = cos(pkin(8));
t159 = cos(qJ(2));
t181 = qJD(1) * t159;
t157 = sin(qJ(2));
t182 = qJD(1) * t157;
t89 = -t153 * t182 - t154 * t181;
t90 = t153 * t181 - t154 * t182;
t166 = t156 * t89 - t158 * t90;
t206 = t156 * t90 + t158 * t89;
t40 = Ifges(6,4) * t206;
t230 = Ifges(6,2) * t166 - t40;
t229 = qJD(2) - qJD(5);
t151 = qJD(2) * qJD(3);
t180 = qJD(2) * t157;
t188 = pkin(6) - qJ(4);
t86 = -qJD(4) * t159 - t180 * t188;
t66 = qJD(1) * t86 + t151;
t130 = t188 * t159;
t88 = qJD(2) * t130 - qJD(4) * t157;
t74 = t88 * qJD(1);
t25 = -t153 * t66 + t154 * t74;
t165 = t153 * t157 + t154 * t159;
t96 = t165 * qJD(2);
t83 = qJD(1) * t96;
t19 = -pkin(7) * t83 + t25;
t26 = t153 * t74 + t154 * t66;
t109 = -t153 * t159 + t154 * t157;
t95 = t109 * qJD(2);
t82 = qJD(1) * t95;
t20 = pkin(7) * t82 + t26;
t192 = pkin(7) * t90;
t143 = pkin(6) * t181;
t118 = -qJ(4) * t181 + t143;
t152 = qJD(2) * qJ(3);
t100 = t118 + t152;
t142 = pkin(6) * t182;
t115 = qJ(4) * t182 - t142;
t160 = -pkin(2) - pkin(3);
t174 = t160 * qJD(2);
t79 = qJD(3) + t174 - t115;
t34 = -t100 * t153 + t154 * t79;
t21 = -qJD(2) * pkin(4) + t192 + t34;
t193 = pkin(7) * t89;
t35 = t100 * t154 + t153 * t79;
t22 = t35 + t193;
t5 = -t156 * t22 + t158 * t21;
t1 = qJD(5) * t5 + t156 * t19 + t158 * t20;
t14 = qJD(5) * t206 + t156 * t82 + t158 * t83;
t15 = -qJD(5) * t166 - t156 * t83 + t158 * t82;
t6 = t156 * t21 + t158 * t22;
t2 = -qJD(5) * t6 - t156 * t20 + t158 * t19;
t105 = -qJD(1) * pkin(1) - pkin(2) * t181 - qJ(3) * t182;
t70 = pkin(3) * t181 + qJD(4) - t105;
t42 = -pkin(4) * t89 + t70;
t228 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t14 + Ifges(6,6) * t15 + t229 * (Ifges(6,5) * t206 - Ifges(6,6) * t166) / 0.2e1 - t42 * (mrSges(6,1) * t166 + mrSges(6,2) * t206);
t107 = -t153 * t156 + t154 * t158;
t55 = -t115 * t153 + t118 * t154;
t27 = t55 + t193;
t56 = t115 * t154 + t118 * t153;
t28 = t56 - t192;
t119 = -qJ(3) * t153 + t154 * t160;
t114 = -pkin(4) + t119;
t120 = qJ(3) * t154 + t153 * t160;
t57 = t114 * t158 - t120 * t156;
t227 = qJD(3) * t107 + qJD(5) * t57 - t156 * t27 - t158 * t28;
t110 = t153 * t158 + t154 * t156;
t58 = t114 * t156 + t120 * t158;
t226 = -qJD(3) * t110 - qJD(5) * t58 + t156 * t28 - t158 * t27;
t126 = t143 + t152;
t128 = mrSges(4,2) * t181 + qJD(2) * mrSges(4,3);
t187 = Ifges(3,4) * t157;
t212 = qJD(2) / 0.2e1;
t213 = -qJD(2) / 0.2e1;
t216 = -Ifges(4,5) * t182 / 0.2e1;
t218 = Ifges(4,3) / 0.2e1;
t225 = -(m(4) * t126 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t181 + t128) * pkin(6) - t126 * mrSges(4,2) - Ifges(4,6) * t213 - t181 * t218 - t216 - Ifges(3,6) * t212 - (Ifges(3,2) * t159 + t187) * qJD(1) / 0.2e1 + t105 * mrSges(4,1);
t122 = -qJD(2) * pkin(2) + qJD(3) + t142;
t186 = Ifges(4,5) * t159;
t215 = -Ifges(3,4) * t181 / 0.2e1;
t219 = -Ifges(3,1) / 0.2e1;
t224 = (m(4) * t122 + (mrSges(4,2) + mrSges(3,3)) * t182 + t231 * qJD(2)) * pkin(6) - t105 * mrSges(4,3) + (t157 * Ifges(4,1) - t186) * qJD(1) / 0.2e1 - t182 * t219 - t215 + t122 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t213;
t223 = t166 * t6 + t206 * t5;
t189 = Ifges(6,4) * t166;
t221 = Ifges(6,1) * t206 - t189;
t12 = Ifges(6,1) * t166 - Ifges(6,5) * t229 + t40;
t217 = t12 / 0.2e1;
t201 = -t166 / 0.2e1;
t211 = t229 * t107;
t210 = t229 * t110;
t203 = -t206 / 0.2e1;
t202 = t206 / 0.2e1;
t200 = t166 / 0.2e1;
t199 = -t90 / 0.2e1;
t197 = t95 / 0.2e1;
t196 = t96 / 0.2e1;
t195 = pkin(1) * mrSges(3,1);
t194 = pkin(1) * mrSges(3,2);
t190 = Ifges(5,4) * t89;
t37 = t153 * t88 + t154 * t86;
t129 = t188 * t157;
t63 = t129 * t153 + t130 * t154;
t138 = qJ(3) * t181;
t146 = t157 * qJD(3);
t184 = qJD(1) * t146 + qJD(2) * t138;
t179 = qJD(2) * t159;
t183 = qJ(3) * t179 + t146;
t123 = -pkin(2) * t159 - qJ(3) * t157 - pkin(1);
t178 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t177 = 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4);
t176 = Ifges(4,6) / 0.2e1 - Ifges(3,6) / 0.2e1;
t175 = m(4) * pkin(6) + mrSges(4,2);
t173 = -t82 * mrSges(5,1) + mrSges(5,2) * t83;
t172 = -t15 * mrSges(6,1) + mrSges(6,2) * t14;
t171 = qJD(1) * t180;
t36 = -t153 * t86 + t154 * t88;
t62 = t129 * t154 - t130 * t153;
t106 = pkin(3) * t159 - t123;
t169 = t157 * t174;
t168 = -t153 * t34 + t154 * t35;
t38 = -pkin(7) * t109 + t62;
t39 = -pkin(7) * t165 + t63;
t9 = -t156 * t39 + t158 * t38;
t10 = t156 * t38 + t158 * t39;
t53 = -t109 * t156 - t158 * t165;
t54 = t109 * t158 - t156 * t165;
t84 = t160 * t182 + t138;
t71 = qJD(2) * mrSges(5,2) + mrSges(5,3) * t89;
t72 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t90;
t164 = t153 * t72 - t154 * t71 - t128;
t69 = t169 + t183;
t61 = qJD(1) * t169 + t184;
t121 = -pkin(6) * t171 + t151;
t116 = (-mrSges(4,1) * t159 - mrSges(4,3) * t157) * qJD(1);
t87 = pkin(2) * t180 - t183;
t85 = Ifges(5,4) * t90;
t73 = pkin(2) * t171 - t184;
t59 = pkin(4) * t165 + t106;
t52 = pkin(4) * t90 + t84;
t51 = -mrSges(5,1) * t89 - mrSges(5,2) * t90;
t44 = -t90 * Ifges(5,1) - Ifges(5,5) * qJD(2) + t190;
t43 = t89 * Ifges(5,2) - Ifges(5,6) * qJD(2) - t85;
t41 = -pkin(4) * t95 + t69;
t33 = -mrSges(6,1) * t229 - mrSges(6,3) * t166;
t32 = mrSges(6,2) * t229 + mrSges(6,3) * t206;
t31 = -pkin(4) * t82 + t61;
t24 = pkin(7) * t95 + t37;
t23 = -pkin(7) * t96 + t36;
t18 = -qJD(5) * t54 - t156 * t96 + t158 * t95;
t17 = qJD(5) * t53 + t156 * t95 + t158 * t96;
t16 = -mrSges(6,1) * t206 + mrSges(6,2) * t166;
t11 = Ifges(6,2) * t206 - Ifges(6,6) * t229 + t189;
t4 = -qJD(5) * t10 - t156 * t24 + t158 * t23;
t3 = qJD(5) * t9 + t156 * t23 + t158 * t24;
t7 = [m(4) * (t105 * t87 + t123 * t73) + t44 * t196 + t59 * t172 + t106 * t173 + t61 * (mrSges(5,1) * t165 + mrSges(5,2) * t109) + (-t165 * t82 + t197 * t89) * Ifges(5,2) + (t109 * t82 - t165 * t83 + t196 * t89 + t199 * t95) * Ifges(5,4) + (-t109 * t25 - t165 * t26 - t34 * t96 + t35 * t95 - t62 * t83 + t63 * t82) * mrSges(5,3) + t87 * t116 + t70 * (-mrSges(5,1) * t95 + mrSges(5,2) * t96) + t69 * t51 + t37 * t71 + t36 * t72 + t31 * (-mrSges(6,1) * t53 + mrSges(6,2) * t54) + (t15 * t53 + t18 * t202) * Ifges(6,2) + (t14 * t53 + t15 * t54 + t17 * t202 + t18 * t200) * Ifges(6,4) + (t109 * t83 + t199 * t96) * Ifges(5,1) + (t14 * t54 + t17 * t200) * Ifges(6,1) - t229 * (Ifges(6,5) * t17 + Ifges(6,6) * t18) / 0.2e1 + t4 * t33 + t41 * t16 + t42 * (-mrSges(6,1) * t18 + mrSges(6,2) * t17) + t3 * t32 + t18 * t11 / 0.2e1 + m(5) * (t106 * t61 + t25 * t62 + t26 * t63 + t34 * t36 + t35 * t37 + t69 * t70) + m(6) * (t1 * t10 + t2 * t9 + t3 * t6 + t31 * t59 + t4 * t5 + t41 * t42) + (-t73 * mrSges(4,1) + t175 * t121 + (t178 * qJD(2) + (-t123 * mrSges(4,3) - 0.2e1 * t194 - t177 * t159 + (0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + t175 * pkin(6)) * t157) * qJD(1) + t224) * qJD(2)) * t159 + (-t73 * mrSges(4,3) + (t176 * qJD(2) + (mrSges(4,1) * t123 + t157 * t177 - 0.2e1 * t195) * qJD(1) + t225) * qJD(2)) * t157 + (t1 * t53 + t10 * t15 - t14 * t9 - t17 * t5 + t18 * t6 - t2 * t54) * mrSges(6,3) + t43 * t197 + (Ifges(5,5) * t96 + Ifges(5,6) * t95) * t213 + t17 * t217; -t164 * qJD(3) + (-m(4) * t105 - t116) * (pkin(2) * t182 - t138) + (-Ifges(5,1) * t89 + t43 - t85) * t90 / 0.2e1 + t206 * t217 + t121 * mrSges(4,3) + t89 * t44 / 0.2e1 - t70 * (mrSges(5,1) * t90 - mrSges(5,2) * t89) - t84 * t51 - Ifges(5,6) * t82 - Ifges(5,5) * t83 - t56 * t71 - t55 * t72 - t52 * t16 + (t11 - t221) * t201 - t89 * (-Ifges(5,2) * t90 - t190) / 0.2e1 - t25 * mrSges(5,1) + t26 * mrSges(5,2) - t228 + ((t215 + (t194 + t186 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t231) * pkin(6) + t178) * qJD(2) - t224) * t159 + (t216 + (t195 + t187 / 0.2e1) * qJD(1) + (-Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + t218 + t219) * t181 + (mrSges(3,2) * pkin(6) - mrSges(4,2) * qJ(3) + t176) * qJD(2) - t225) * t157) * qJD(1) + m(4) * (qJ(3) * t121 + qJD(3) * t126) + (qJD(3) * t168 + t119 * t25 + t120 * t26 - t34 * t55 - t35 * t56 - t70 * t84) * m(5) + t226 * t33 + (t1 * t58 + t2 * t57 + t226 * t5 + t227 * t6 - t42 * t52) * m(6) + t227 * t32 + (-t14 * t57 + t15 * t58 - t223) * mrSges(6,3) + (-Ifges(5,5) * t89 - Ifges(5,6) * t90) * t212 + (-t119 * t83 + t120 * t82 - t34 * t89 + t35 * t90) * mrSges(5,3) + t230 * t203; t210 * t33 - t211 * t32 + (-t107 * t14 + t110 * t15) * mrSges(6,3) + (t153 * t82 - t154 * t83) * mrSges(5,3) + t164 * qJD(2) + (t175 * t179 + (t116 - t16 - t51) * t157) * qJD(1) - m(4) * (qJD(2) * t126 - t105 * t182) + (t1 * t110 + t2 * t107 - t182 * t42 + t210 * t5 - t211 * t6) * m(6) + (-qJD(2) * t168 + t153 * t26 + t154 * t25 - t182 * t70) * m(5); -t206 * t32 + t166 * t33 - t89 * t71 - t90 * t72 + t172 + t173 + (t166 * t5 - t206 * t6 + t31) * m(6) + (-t34 * t90 - t35 * t89 + t61) * m(5); t221 * t201 + t11 * t200 - t5 * t32 + t6 * t33 + t223 * mrSges(6,3) + (t12 - t230) * t203 + t228;];
tauc = t7(:);
