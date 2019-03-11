% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:43
% EndTime: 2019-03-09 02:02:51
% DurationCPUTime: 3.87s
% Computational Cost: add. (2771->363), mult. (5900->463), div. (0->0), fcn. (3052->6), ass. (0->173)
t231 = -qJD(4) / 0.2e1;
t229 = Ifges(6,1) + Ifges(7,1);
t223 = Ifges(7,4) + Ifges(6,5);
t146 = Ifges(5,5) * t231;
t100 = sin(qJ(5));
t102 = cos(qJ(5));
t101 = sin(qJ(4));
t103 = cos(qJ(4));
t163 = qJD(2) * t103;
t93 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t79 = qJD(1) * t93 + qJD(3);
t64 = t101 * t79 + t163;
t53 = qJD(4) * pkin(8) + t64;
t95 = sin(pkin(9)) * pkin(1) + qJ(3);
t75 = t101 * pkin(4) - pkin(8) * t103 + t95;
t66 = t75 * qJD(1);
t17 = -t100 * t53 + t102 * t66;
t18 = t100 * t66 + t102 * t53;
t118 = t18 * t100 + t17 * t102;
t217 = qJD(6) - t17;
t165 = qJD(1) * t101;
t94 = qJD(5) + t165;
t10 = -pkin(5) * t94 + t217;
t11 = qJ(6) * t94 + t18;
t119 = t10 * t102 - t11 * t100;
t175 = Ifges(7,5) * t102;
t123 = Ifges(7,3) * t100 + t175;
t179 = Ifges(6,4) * t102;
t129 = -Ifges(6,2) * t100 + t179;
t134 = mrSges(7,1) * t100 - mrSges(7,3) * t102;
t136 = mrSges(6,1) * t100 + mrSges(6,2) * t102;
t63 = -t101 * qJD(2) + t103 * t79;
t52 = -qJD(4) * pkin(4) - t63;
t158 = t102 * qJD(4);
t164 = qJD(1) * t103;
t82 = t100 * t164 - t158;
t159 = t100 * qJD(4);
t83 = t102 * t164 + t159;
t16 = t82 * pkin(5) - t83 * qJ(6) + t52;
t173 = Ifges(7,6) * t100;
t174 = Ifges(6,6) * t100;
t177 = Ifges(6,5) * t102;
t178 = Ifges(7,4) * t102;
t196 = t102 / 0.2e1;
t198 = t100 / 0.2e1;
t199 = -t100 / 0.2e1;
t200 = t94 / 0.2e1;
t202 = t83 / 0.2e1;
t204 = t82 / 0.2e1;
t205 = -t82 / 0.2e1;
t193 = Ifges(7,5) * t82;
t77 = Ifges(6,4) * t82;
t218 = t223 * t94 + t229 * t83 + t193 - t77;
t176 = Ifges(7,5) * t100;
t180 = Ifges(6,4) * t100;
t225 = t229 * t102 + t176 - t180;
t76 = Ifges(7,5) * t83;
t27 = t94 * Ifges(7,6) + t82 * Ifges(7,3) + t76;
t194 = Ifges(6,4) * t83;
t30 = -t82 * Ifges(6,2) + t94 * Ifges(6,6) + t194;
t230 = -t196 * t218 - t198 * t27 - t199 * t30 - t123 * t204 - t129 * t205 - t16 * t134 - t52 * t136 - (t173 + t178 - t174 + t177) * t200 - t225 * t202 - t119 * mrSges(7,2) + t118 * mrSges(6,3);
t228 = Ifges(6,6) - Ifges(7,6);
t157 = qJD(1) * qJD(4);
t143 = t103 * t157;
t161 = qJD(5) * t100;
t57 = qJD(5) * t158 + (-t101 * t158 - t103 * t161) * qJD(1);
t166 = t100 * t101;
t58 = qJD(5) * t83 - t157 * t166;
t227 = (-Ifges(6,4) + Ifges(7,5)) * t58 + t229 * t57 + t223 * t143;
t181 = Ifges(5,4) * t101;
t224 = qJD(1) / 0.2e1;
t226 = t146 - (t103 * Ifges(5,1) - t181) * t224 + t63 * mrSges(5,3) + t230;
t120 = pkin(5) * t100 - qJ(6) * t102;
t222 = qJD(5) * t120 - qJD(6) * t100 - t163 - (-qJD(1) * t120 + t79) * t101;
t54 = qJD(4) * t63;
t221 = t229 * t100 - t175 + t179;
t220 = qJD(4) * t101;
t219 = -m(5) - m(4);
t145 = Ifges(5,6) * t231;
t215 = t223 * t100 + t228 * t102;
t160 = qJD(5) * t102;
t140 = pkin(4) * t103 + pkin(8) * t101;
t80 = qJD(4) * t140 + qJD(3);
t68 = t80 * qJD(1);
t3 = t100 * t68 + t102 * t54 + t66 * t160 - t161 * t53;
t4 = -qJD(5) * t18 - t100 * t54 + t102 * t68;
t138 = -t4 * t100 + t3 * t102;
t1 = qJ(6) * t143 + t94 * qJD(6) + t3;
t2 = -pkin(5) * t143 - t4;
t139 = t1 * t102 + t2 * t100;
t183 = t102 * t101 * t93 + t100 * t75;
t214 = -qJD(5) * t183 + t102 * t80;
t20 = mrSges(7,1) * t58 - mrSges(7,3) * t57;
t21 = mrSges(6,1) * t58 + mrSges(6,2) * t57;
t55 = qJD(4) * t64;
t5 = t58 * pkin(5) - t57 * qJ(6) - t83 * qJD(6) + t55;
t213 = -m(6) * (-t158 * t18 + t159 * t17 + t55) + m(7) * (t10 * t159 + t11 * t158 - t5) - t20 - t21;
t195 = mrSges(6,3) * t83;
t60 = mrSges(6,1) * t94 - t195;
t61 = -mrSges(7,1) * t94 + t83 * mrSges(7,2);
t184 = t60 - t61;
t188 = t82 * mrSges(6,3);
t59 = -mrSges(6,2) * t94 - t188;
t62 = -mrSges(7,2) * t82 + mrSges(7,3) * t94;
t185 = t59 + t62;
t111 = -t100 * t185 - t102 * t184;
t212 = -m(6) * t118 + m(7) * t119 + t111;
t149 = mrSges(5,3) * t165;
t88 = -qJD(4) * mrSges(5,2) - t149;
t211 = -t100 * t184 + t102 * t185 + t88;
t152 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t153 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t154 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t87 = qJD(1) * t95;
t210 = -t152 * t94 - t153 * t82 + t154 * t83 - t11 * mrSges(7,3) - t17 * mrSges(6,1) - t87 * mrSges(5,1) - Ifges(6,6) * t205 - Ifges(7,6) * t204 - t145 + (Ifges(5,4) * t103 - Ifges(5,2) * t101) * t224 + t10 * mrSges(7,1) + t18 * mrSges(6,2) - t223 * t202 - (Ifges(6,3) + Ifges(7,2)) * t200;
t208 = t57 / 0.2e1;
t207 = -t58 / 0.2e1;
t206 = t58 / 0.2e1;
t203 = -t83 / 0.2e1;
t201 = -t94 / 0.2e1;
t197 = -t102 / 0.2e1;
t85 = t140 * qJD(1);
t25 = t100 * t85 + t102 * t63;
t148 = mrSges(5,3) * t164;
t182 = qJD(4) * mrSges(5,1) - mrSges(6,1) * t82 - mrSges(6,2) * t83 - t148;
t172 = t102 * t75;
t169 = mrSges(4,3) * qJD(1);
t162 = qJD(4) * t103;
t42 = mrSges(7,1) * t82 - mrSges(7,3) * t83;
t156 = t42 - t182;
t151 = t93 * t162;
t155 = t100 * t80 + t102 * t151 + t75 * t160;
t150 = t93 * t161;
t147 = m(5) * t93 - mrSges(5,3);
t144 = t100 * t93 - pkin(5);
t142 = 0.2e1 * t87;
t137 = mrSges(6,1) * t102 - mrSges(6,2) * t100;
t135 = mrSges(7,1) * t102 + mrSges(7,3) * t100;
t128 = Ifges(6,2) * t102 + t180;
t122 = -Ifges(7,3) * t102 + t176;
t121 = pkin(5) * t102 + qJ(6) * t100;
t24 = -t100 * t63 + t102 * t85;
t37 = -mrSges(7,1) * t143 + t57 * mrSges(7,2);
t116 = t120 - t93;
t35 = -t58 * mrSges(7,2) + mrSges(7,3) * t143;
t36 = mrSges(6,1) * t143 - t57 * mrSges(6,3);
t38 = -mrSges(6,2) * t143 - t58 * mrSges(6,3);
t112 = (t35 + t38) * t102 + (-t36 + t37) * t100;
t108 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t105 = m(6) * (qJD(4) * t52 - t160 * t17 - t161 * t18 + t138) + m(7) * (qJD(4) * t16 + t10 * t160 - t11 * t161 + t139) + t111 * qJD(5) + t112;
t92 = Ifges(7,2) * t143;
t91 = Ifges(6,3) * t143;
t86 = -pkin(4) - t121;
t84 = (t101 * mrSges(5,1) + mrSges(5,2) * t103) * qJD(1);
t56 = t116 * t103;
t50 = Ifges(7,4) * t57;
t49 = Ifges(6,5) * t57;
t48 = Ifges(6,6) * t58;
t47 = Ifges(7,6) * t58;
t41 = pkin(5) * t83 + qJ(6) * t82;
t39 = -t166 * t93 + t172;
t34 = t101 * t144 - t172;
t33 = qJ(6) * t101 + t183;
t23 = -pkin(5) * t164 - t24;
t22 = qJ(6) * t164 + t25;
t19 = (qJD(5) * t121 - qJD(6) * t102) * t103 - t116 * t220;
t13 = Ifges(6,4) * t57 - Ifges(6,2) * t58 + Ifges(6,6) * t143;
t12 = Ifges(7,5) * t57 + Ifges(7,6) * t143 + Ifges(7,3) * t58;
t9 = -t100 * t151 + t214;
t8 = -t101 * t150 + t155;
t7 = t144 * t162 - t214;
t6 = qJ(6) * t162 + (qJD(6) - t150) * t101 + t155;
t14 = [t19 * t42 + t56 * t20 + t33 * t35 + t34 * t37 + t39 * t36 + t183 * t38 + t8 * t59 + t6 * t62 + t9 * t60 + t7 * t61 + m(7) * (t1 * t33 + t10 * t7 + t11 * t6 + t16 * t19 + t2 * t34 + t5 * t56) + m(6) * (t17 * t9 + t18 * t8 + t183 * t3 + t4 * t39) + (-t142 * t219 + 0.2e1 * t169 + t84) * qJD(3) + (t91 / 0.2e1 + t92 / 0.2e1 - t48 / 0.2e1 + t49 / 0.2e1 + t50 / 0.2e1 + t47 / 0.2e1 + qJD(1) * qJD(3) * mrSges(5,1) + t153 * t58 - t154 * t57 + t147 * t54 + (-t142 * mrSges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,4) * t165 + t146 + t226) * qJD(4) + t108) * t101 + (t12 * t198 + t13 * t199 + t5 * t134 + t123 * t206 + t129 * t207 + (-t100 * t3 - t102 * t4) * mrSges(6,3) + (-t1 * t100 + t102 * t2) * mrSges(7,2) + (mrSges(5,3) + t136) * t55 + (t147 * t64 + t145 - t210) * qJD(4) + (t128 * t204 + t122 * t205 + t52 * t137 + t16 * t135 + t30 * t197 + (t100 * t17 - t102 * t18) * mrSges(6,3) + (-t10 * t100 - t102 * t11) * mrSges(7,2) + t221 * t203 + t215 * t201 + t218 * t199) * qJD(5) + (qJD(3) * mrSges(5,2) + (t95 * mrSges(5,1) + (-0.3e1 / 0.2e1 * Ifges(5,4) + t178 / 0.2e1 + t173 / 0.2e1 + t177 / 0.2e1 - t174 / 0.2e1) * t103 + (0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,1) + t152) * t101) * qJD(4)) * qJD(1) + t225 * t208 + (t27 * qJD(5) + t227) * t196) * t103 + ((-m(5) * t63 + m(6) * t52 - t182) * t220 + (-t21 + (-m(6) - m(5)) * t55 + qJD(4) * t88) * t103) * t93; ((-t149 - t211) * qJD(4) - t213) * t101 + ((-t148 + t156) * qJD(4) + t105) * t103; (t219 * t87 - t169 + t212 - t84) * qJD(1) + (qJD(4) * t211 + t213) * t103 + (qJD(4) * t156 + t105) * t101; t112 * pkin(8) + (pkin(8) * t212 - t230) * qJD(5) + t139 * mrSges(7,2) + (pkin(8) * t139 - t10 * t23 - t11 * t22 + t222 * t16 + t5 * t86) * m(7) + t222 * t42 + t221 * t208 - t5 * t135 + (-mrSges(5,1) - t137) * t55 + t138 * mrSges(6,3) + (-pkin(4) * t55 + t138 * pkin(8) - t17 * t24 - t18 * t25 - t52 * t64) * m(6) - t63 * t88 + t86 * t20 - t25 * t59 - t24 * t60 - t23 * t61 - t22 * t62 - t54 * mrSges(5,2) - pkin(4) * t21 + t182 * t64 + ((Ifges(5,4) * t164 / 0.2e1 + t64 * mrSges(5,3) + t145 + t215 * qJD(4) / 0.2e1 + t210) * t103 + (t87 * mrSges(5,2) + (-t181 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t103) * qJD(1) + t146 - t226) * t101) * qJD(1) + t227 * t198 + t13 * t196 + t12 * t197 + t122 * t206 + t128 * t207; t108 + t91 + t92 - t48 + t49 + t50 + t47 - t16 * (mrSges(7,1) * t83 + mrSges(7,3) * t82) - t52 * (mrSges(6,1) * t83 - mrSges(6,2) * t82) + qJD(6) * t62 - pkin(5) * t37 - t41 * t42 + qJ(6) * t35 + (t10 * t82 + t11 * t83) * mrSges(7,2) + (t184 + t195) * t18 + (-t185 - t188) * t17 + t30 * t202 + (Ifges(7,3) * t83 - t193) * t205 + (-t223 * t82 - t228 * t83) * t201 + (-pkin(5) * t2 + qJ(6) * t1 - t10 * t18 + t11 * t217 - t16 * t41) * m(7) + (-Ifges(6,2) * t83 + t218 - t77) * t204 + (-t229 * t82 - t194 + t27 + t76) * t203; t83 * t42 - t94 * t62 + 0.2e1 * (t2 / 0.2e1 + t11 * t201 + t16 * t202) * m(7) + t37;];
tauc  = t14(:);
