% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:45
% EndTime: 2019-07-18 17:20:57
% DurationCPUTime: 3.07s
% Computational Cost: add. (3023->336), mult. (7551->468), div. (0->0), fcn. (4654->6), ass. (0->170)
t136 = qJD(2) + qJD(4);
t260 = t136 * Ifges(5,6) / 0.2e1;
t259 = Ifges(3,4) + Ifges(4,4);
t143 = cos(qJ(2));
t144 = pkin(2) + pkin(1);
t121 = t144 * t143;
t107 = -qJD(1) * t121 + qJD(3);
t210 = t136 * Ifges(5,5);
t258 = t107 * mrSges(5,2) + t210 / 0.2e1;
t183 = t144 * qJD(2);
t139 = sin(qJ(4));
t140 = sin(qJ(2));
t142 = cos(qJ(4));
t159 = t139 * t140 - t142 * t143;
t102 = t159 * qJD(1);
t110 = t139 * t143 + t140 * t142;
t103 = t110 * qJD(1);
t211 = t103 * Ifges(5,4);
t257 = t260 + t211 / 0.2e1 - t102 * Ifges(5,2) / 0.2e1;
t254 = -qJD(2) / 0.2e1;
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t221 = pkin(3) + qJ(3);
t178 = qJD(2) * t221;
t175 = t143 * t178;
t193 = qJD(4) * t142;
t120 = t221 * t143;
t113 = qJD(1) * t120;
t200 = t113 * t139;
t189 = qJD(1) * qJD(3);
t128 = t143 * t189;
t176 = t140 * t178;
t92 = -qJD(1) * t176 + t128;
t18 = -qJD(4) * t200 + t142 * t92 + t183 * t193 + (-t139 * t175 + (-t139 * qJD(3) - t221 * t193) * t140) * qJD(1);
t119 = t221 * t140;
t111 = qJD(1) * t119;
t150 = t183 - t111;
t199 = t142 * t113;
t64 = t139 * t150 + t199;
t56 = t136 * pkin(4) + t64;
t78 = -pkin(4) * t103 + t107;
t24 = -t138 * t56 + t141 * t78;
t190 = qJD(1) * qJD(2);
t181 = t140 * t190;
t104 = t144 * t181;
t79 = t136 * t159;
t71 = t79 * qJD(1);
t50 = pkin(4) * t71 + t104;
t2 = t24 * qJD(5) + t138 * t50 + t141 * t18;
t25 = t138 * t78 + t141 * t56;
t3 = -t25 * qJD(5) - t138 * t18 + t141 * t50;
t252 = -t3 * t138 + t141 * t2;
t96 = qJD(5) + t102;
t251 = -t107 * mrSges(5,1) - t24 * mrSges(6,1) + t25 * mrSges(6,2) + t257;
t170 = mrSges(6,1) * t138 + mrSges(6,2) * t141;
t63 = -t142 * t150 + t200;
t157 = t63 * t170;
t165 = Ifges(6,5) * t141 - Ifges(6,6) * t138;
t217 = Ifges(6,4) * t141;
t167 = -Ifges(6,2) * t138 + t217;
t218 = Ifges(6,4) * t138;
t169 = Ifges(6,1) * t141 - t218;
t229 = t141 / 0.2e1;
t230 = -t138 / 0.2e1;
t85 = t103 * t141 + t136 * t138;
t234 = t85 / 0.2e1;
t228 = Ifges(6,4) * t85;
t83 = -t103 * t138 + t136 * t141;
t30 = Ifges(6,2) * t83 + Ifges(6,6) * t96 + t228;
t82 = Ifges(6,4) * t83;
t31 = Ifges(6,1) * t85 + Ifges(6,5) * t96 + t82;
t250 = (-t138 * t25 - t141 * t24) * mrSges(6,3) + t96 * t165 / 0.2e1 + t169 * t234 + t83 * t167 / 0.2e1 + t157 + t31 * t229 + t30 * t230;
t249 = -t30 / 0.2e1;
t213 = t103 * mrSges(5,3);
t220 = -mrSges(5,1) * t136 - mrSges(6,1) * t83 + mrSges(6,2) * t85 + t213;
t244 = -t138 * t24 + t141 * t25;
t36 = t83 * qJD(5) - t141 * t71;
t37 = -t85 * qJD(5) + t138 * t71;
t243 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t36 + Ifges(6,6) * t37;
t196 = qJD(1) * t143;
t242 = -t259 * t196 / 0.2e1;
t118 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t196;
t122 = -pkin(1) * t196 + qJD(3);
t241 = (m(4) * t122 + (-mrSges(4,1) * t143 + mrSges(4,2) * t140) * qJD(1)) * pkin(1) + t122 * mrSges(4,1) - qJ(3) * t118 - ((Ifges(3,2) + Ifges(4,2)) * t143 + t259 * t140) * qJD(1) / 0.2e1 + (Ifges(4,6) + Ifges(3,6)) * t254;
t197 = qJD(1) * t140;
t215 = qJD(2) * pkin(1);
t116 = -qJ(3) * t197 + t215;
t117 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t197;
t240 = (m(4) * t116 + t117) * qJ(3) - t122 * mrSges(4,2) + t116 * mrSges(4,3) + t242 - (Ifges(3,1) + Ifges(4,1)) * t197 / 0.2e1 + (Ifges(3,5) + Ifges(4,5)) * t254;
t239 = t36 / 0.2e1;
t238 = t37 / 0.2e1;
t80 = t136 * t110;
t72 = t80 * qJD(1);
t237 = t72 / 0.2e1;
t236 = -t83 / 0.2e1;
t235 = -t85 / 0.2e1;
t233 = -t96 / 0.2e1;
t232 = t102 / 0.2e1;
t231 = -t103 / 0.2e1;
t75 = -t111 * t139 + t199;
t225 = t63 * t75;
t224 = t83 * Ifges(6,6);
t223 = t85 * Ifges(6,5);
t222 = t96 * Ifges(6,3);
t219 = mrSges(5,3) * t102;
t94 = Ifges(5,4) * t102;
t216 = qJ(3) * mrSges(4,3);
t212 = t103 * Ifges(5,1);
t194 = qJD(3) * t140;
t153 = -t175 - t194;
t149 = t142 * t153;
t19 = -qJD(1) * t149 + t64 * qJD(4) + t139 * t92;
t205 = t142 * t19;
t204 = t63 * t103;
t203 = qJ(3) * t143 ^ 2;
t202 = t102 * t138;
t201 = t102 * t141;
t198 = t143 * mrSges(4,2) * t190 + mrSges(4,1) * t181;
t114 = t144 * t197;
t115 = t140 * t183;
t192 = qJD(5) * t138;
t191 = qJD(5) * t141;
t188 = 0.3e1 / 0.2e1 * Ifges(4,4) + 0.3e1 / 0.2e1 * Ifges(3,4);
t187 = Ifges(3,5) / 0.2e1 + Ifges(4,5) / 0.2e1;
t182 = t72 * mrSges(5,1) - t71 * mrSges(5,2);
t174 = (-Ifges(3,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(2);
t173 = -t2 * t138 - t3 * t141;
t160 = -t142 * t119 - t120 * t139;
t86 = -t119 * t139 + t120 * t142;
t95 = qJD(3) * t143 - t176;
t35 = t86 * qJD(4) + t139 * t95 - t149;
t172 = -t160 * t19 + t35 * t63;
t171 = mrSges(6,1) * t141 - mrSges(6,2) * t138;
t168 = Ifges(6,1) * t138 + t217;
t166 = Ifges(6,2) * t141 + t218;
t164 = Ifges(6,5) * t138 + Ifges(6,6) * t141;
t52 = -mrSges(6,2) * t96 + mrSges(6,3) * t83;
t53 = mrSges(6,1) * t96 - mrSges(6,3) * t85;
t161 = -t138 * t53 + t141 * t52;
t88 = -pkin(4) * t110 - t121;
t45 = t138 * t88 + t141 * t86;
t44 = -t138 * t86 + t141 * t88;
t89 = -mrSges(5,2) * t136 - t219;
t158 = -t161 - t89;
t11 = mrSges(6,1) * t72 - mrSges(6,3) * t36;
t12 = -mrSges(6,2) * t72 + mrSges(6,3) * t37;
t148 = m(6) * (-t24 * t191 - t25 * t192 + t252) + t141 * t12 - t138 * t11 - t53 * t191 - t52 * t192;
t29 = t222 + t223 + t224;
t61 = t210 - t94 + t212;
t8 = t36 * Ifges(6,4) + t37 * Ifges(6,2) + t72 * Ifges(6,6);
t9 = t36 * Ifges(6,1) + t37 * Ifges(6,4) + t72 * Ifges(6,5);
t147 = t31 * t201 / 0.2e1 + t138 * t9 / 0.2e1 - Ifges(5,5) * t71 - Ifges(5,6) * t72 - t18 * mrSges(5,2) + t202 * t249 + t64 * t213 + t63 * t219 + t8 * t229 + t164 * t237 + t166 * t238 + t168 * t239 + (-t94 + t61) * t232 + (-t211 + t29) * t231 + (-mrSges(5,1) - t171) * t19 + (-Ifges(5,1) * t231 - t165 * t233 - t167 * t236 - t169 * t235 + t157 + t258) * t102 + (Ifges(6,5) * t235 - Ifges(5,2) * t232 + Ifges(6,6) * t236 + Ifges(6,3) * t233 + t251 + t260) * t103 + (-t24 * t201 - t25 * t202 + t252) * mrSges(6,3) + t250 * qJD(5);
t146 = qJ(3) ^ 2;
t145 = qJD(1) ^ 2;
t106 = (-qJ(3) * qJD(2) * t143 - t194) * qJD(1);
t105 = -qJ(3) * t181 + t128;
t81 = pkin(4) * t102 + t114;
t76 = -t111 * t142 - t200;
t74 = mrSges(5,1) * t102 + mrSges(5,2) * t103;
t68 = Ifges(6,3) * t72;
t55 = pkin(4) * t79 + t115;
t41 = pkin(4) * t202 - t141 * t63;
t40 = pkin(4) * t201 + t138 * t63;
t39 = t138 * t81 + t141 * t76;
t38 = -t138 * t76 + t141 * t81;
t34 = t160 * qJD(4) + t139 * t153 + t142 * t95;
t10 = -mrSges(6,1) * t37 + mrSges(6,2) * t36;
t7 = -t45 * qJD(5) - t138 * t34 + t141 * t55;
t6 = t44 * qJD(5) + t138 * t55 + t141 * t34;
t1 = [t115 * t74 - t121 * t182 + t34 * t89 - t160 * t10 + t6 * t52 + t7 * t53 + t44 * t11 + t45 * t12 + (t29 / 0.2e1 + t222 / 0.2e1 + t224 / 0.2e1 + t223 / 0.2e1 - t251 - t257) * t80 - (-t94 / 0.2e1 + t212 / 0.2e1 + t61 / 0.2e1 + t250 + t258) * t79 + t220 * t35 + m(5) * (-t104 * t121 + t107 * t115 + t18 * t86 + t34 * t64 + t172) + m(6) * (t2 * t45 + t24 * t7 + t25 * t6 + t3 * t44 + t172) + m(4) * t189 * t203 + (t160 * t71 - t63 * t79 - t64 * t80 - t72 * t86) * mrSges(5,3) + (-pkin(1) * t198 + qJD(3) * t118 + (m(4) * qJ(3) + mrSges(4,3)) * t105 + (t187 * qJD(2) + t188 * t196 - t240) * qJD(2)) * t143 + (-t18 * mrSges(5,3) + t68 / 0.2e1 + Ifges(5,4) * t71 + t104 * mrSges(5,1) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t72 + t243) * t159 + (t165 * t237 + t169 * t239 + t167 * t238 + t9 * t229 + t8 * t230 + t104 * mrSges(5,2) - Ifges(5,1) * t71 - Ifges(5,4) * t72 + (mrSges(5,3) + t170) * t19 + t173 * mrSges(6,3) + (-mrSges(6,3) * t244 + t141 * t249 + t164 * t233 + t166 * t236 + t168 * t235 + t63 * t171 + t31 * t230) * qJD(5)) * t110 + (m(4) * (-qJ(3) * t106 - qJD(3) * t116) - qJD(3) * t117 - t106 * mrSges(4,3) + (t174 + ((pkin(1) * mrSges(4,2) - t188) * t140 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + m(4) * (-pkin(1) ^ 2 - t146) - t216 - pkin(1) * mrSges(4,1)) * t143) * qJD(1) + t241) * qJD(2)) * t140; -t220 * t75 - t114 * t74 - t105 * mrSges(4,2) - t76 * t89 - t39 * t52 - t38 * t53 + (m(4) * pkin(1) + mrSges(4,1)) * t106 - m(5) * (t107 * t114 + t64 * t76 + t225) - m(6) * (t24 * t38 + t25 * t39 + t225) + t147 + (((-pkin(1) * mrSges(4,3) + t187) * qJD(2) + t240 + t242) * t143 + ((Ifges(4,4) / 0.2e1 + Ifges(3,4) / 0.2e1) * t197 + t174 + (Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(4,1) / 0.2e1 - Ifges(3,1) / 0.2e1) * t196 - t241) * t140) * qJD(1) + t148 * (t139 * t144 + pkin(4)) + (-t142 * t10 + (-t139 * t72 + t142 * t71) * mrSges(5,3) - m(6) * t205 + m(5) * (t139 * t18 - t205) + ((m(5) * t64 + m(6) * t244 - t158) * t142 + (t220 + (m(6) + m(5)) * t63) * t139) * qJD(4)) * t144 + (m(4) * t146 + t216) * t140 * t143 * t145; t141 * t11 + t138 * t12 - t220 * t103 + t161 * qJD(5) - t158 * t102 + (-t143 * t118 + (m(4) * t215 + t117) * t140) * qJD(1) - m(4) * (-t116 * t197 + t145 * t203) + t182 + t198 + (t244 * t96 - t173 - t204) * m(6) + (t102 * t64 + t104 - t204) * m(5); t63 * t89 - t41 * t52 - t40 * t53 + t147 - m(6) * (t24 * t40 + t25 * t41) + t148 * pkin(4) + (-m(6) * t63 - t220) * t64; t68 - t63 * (mrSges(6,1) * t85 + mrSges(6,2) * t83) + (Ifges(6,1) * t83 - t228) * t235 + t30 * t234 + (Ifges(6,5) * t83 - Ifges(6,6) * t85) * t233 - t24 * t52 + t25 * t53 + (t24 * t83 + t25 * t85) * mrSges(6,3) + (-Ifges(6,2) * t85 + t31 + t82) * t236 + t243;];
tauc  = t1(:);
