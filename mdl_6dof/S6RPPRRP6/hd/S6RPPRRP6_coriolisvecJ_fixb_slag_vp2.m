% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2018-11-23 15:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:46:22
% EndTime: 2018-11-23 15:46:26
% DurationCPUTime: 3.60s
% Computational Cost: add. (2469->396), mult. (5133->502), div. (0->0), fcn. (2511->4), ass. (0->182)
t236 = Ifges(6,1) + Ifges(7,1);
t233 = Ifges(7,4) + Ifges(6,5);
t238 = Ifges(6,6) - Ifges(7,6);
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t101 = pkin(1) + qJ(3);
t103 = sin(qJ(4));
t105 = cos(qJ(4));
t81 = pkin(4) * t103 - pkin(8) * t105 + t101;
t61 = qJD(1) * t81 - qJD(2);
t95 = qJD(1) * qJ(2) + qJD(3);
t88 = -qJD(1) * pkin(7) + t95;
t80 = t103 * t88;
t68 = qJD(4) * pkin(8) + t80;
t20 = -t102 * t68 + t104 * t61;
t21 = t102 * t61 + t104 * t68;
t120 = t102 * t21 + t104 * t20;
t225 = qJD(6) - t20;
t169 = qJD(1) * t103;
t93 = qJD(5) + t169;
t17 = -pkin(5) * t93 + t225;
t18 = qJ(6) * t93 + t21;
t122 = t102 * t18 - t104 * t17;
t184 = Ifges(7,5) * t104;
t126 = Ifges(7,3) * t102 + t184;
t188 = Ifges(6,4) * t104;
t132 = -Ifges(6,2) * t102 + t188;
t137 = mrSges(7,1) * t102 - mrSges(7,3) * t104;
t139 = mrSges(6,1) * t102 + mrSges(6,2) * t104;
t182 = Ifges(7,6) * t102;
t183 = Ifges(6,6) * t102;
t186 = Ifges(6,5) * t104;
t187 = Ifges(7,4) * t104;
t178 = t105 * t88;
t69 = -qJD(4) * pkin(4) - t178;
t161 = t104 * qJD(4);
t168 = qJD(1) * t105;
t76 = t102 * t168 - t161;
t166 = qJD(4) * t102;
t77 = t104 * t168 + t166;
t19 = pkin(5) * t76 - qJ(6) * t77 + t69;
t205 = t104 / 0.2e1;
t207 = t102 / 0.2e1;
t208 = -t102 / 0.2e1;
t209 = t93 / 0.2e1;
t211 = t77 / 0.2e1;
t213 = t76 / 0.2e1;
t214 = -t76 / 0.2e1;
t72 = Ifges(7,5) * t77;
t22 = t93 * Ifges(7,6) + t76 * Ifges(7,3) + t72;
t185 = Ifges(7,5) * t102;
t189 = Ifges(6,4) * t102;
t223 = t104 * t236 + t185 - t189;
t201 = Ifges(7,5) * t76;
t73 = Ifges(6,4) * t76;
t226 = t233 * t93 + t236 * t77 + t201 - t73;
t202 = Ifges(6,4) * t77;
t25 = -t76 * Ifges(6,2) + t93 * Ifges(6,6) + t202;
t237 = t205 * t226 + t207 * t22 + t208 * t25 + t126 * t213 + t132 * t214 + t19 * t137 + t69 * t139 + (t182 + t187 - t183 + t186) * t209 + t223 * t211 - t122 * mrSges(7,2) - t120 * mrSges(6,3);
t177 = Ifges(5,5) * qJD(4);
t190 = Ifges(5,4) * t103;
t234 = qJD(1) / 0.2e1;
t235 = t177 / 0.2e1 + (t105 * Ifges(5,1) - t190) * t234 + t237;
t158 = qJD(1) * qJD(4);
t149 = t105 * t158;
t163 = qJD(5) * t102;
t47 = qJD(5) * t161 + (-t103 * t161 - t105 * t163) * qJD(1);
t171 = t102 * t103;
t48 = qJD(5) * t77 - t158 * t171;
t232 = (-Ifges(6,4) + Ifges(7,5)) * t48 + t236 * t47 + t233 * t149;
t123 = pkin(5) * t102 - qJ(6) * t104;
t231 = -qJD(6) * t102 + t93 * t123 - t80;
t203 = mrSges(6,3) * t77;
t53 = mrSges(6,1) * t93 - t203;
t54 = -mrSges(7,1) * t93 + mrSges(7,2) * t77;
t192 = t53 - t54;
t204 = mrSges(6,3) * t76;
t52 = -mrSges(6,2) * t93 - t204;
t55 = -mrSges(7,2) * t76 + mrSges(7,3) * t93;
t193 = t52 + t55;
t113 = t102 * t192 - t104 * t193;
t172 = qJD(4) * mrSges(5,2);
t84 = -mrSges(5,3) * t169 - t172;
t230 = -t113 + t84;
t229 = t233 * t102 + t238 * t104;
t228 = t236 * t102 - t184 + t188;
t100 = qJ(2) - pkin(7);
t164 = qJD(4) * t105;
t227 = qJD(2) * t103 + t100 * t164;
t150 = -Ifges(5,6) * qJD(4) / 0.2e1;
t162 = qJD(5) * t104;
t160 = qJD(1) * qJD(2);
t60 = t103 * t160 + t164 * t88;
t146 = pkin(4) * t105 + pkin(8) * t103;
t75 = qJD(4) * t146 + qJD(3);
t62 = t75 * qJD(1);
t3 = t102 * t62 + t104 * t60 + t61 * t162 - t163 * t68;
t4 = -qJD(5) * t21 - t102 * t60 + t104 * t62;
t142 = -t4 * t102 + t3 * t104;
t1 = qJ(6) * t149 + qJD(6) * t93 + t3;
t2 = -pkin(5) * t149 - t4;
t144 = t1 * t104 + t2 * t102;
t222 = qJD(1) * t101;
t221 = (t103 ^ 2 + t105 ^ 2) * t88;
t152 = qJD(5) * t100 * t103;
t9 = -t102 * (qJD(5) * t81 + t227) - t104 * (-t75 + t152);
t112 = -t102 * t193 - t104 * t192;
t220 = -m(6) * t120 - m(7) * t122 + t112;
t154 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t155 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t156 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t89 = -qJD(2) + t222;
t219 = -t154 * t76 - t155 * t93 - t156 * t77 - t18 * mrSges(7,3) - t20 * mrSges(6,1) - t89 * mrSges(5,1) - Ifges(6,6) * t214 - Ifges(7,6) * t213 - t150 + (Ifges(5,4) * t105 - Ifges(5,2) * t103) * t234 + t17 * mrSges(7,1) + t21 * mrSges(6,2) - t233 * t211 - (Ifges(6,3) + Ifges(7,2)) * t209;
t217 = t47 / 0.2e1;
t216 = -t48 / 0.2e1;
t215 = t48 / 0.2e1;
t212 = -t77 / 0.2e1;
t210 = -t93 / 0.2e1;
t206 = -t104 / 0.2e1;
t30 = -mrSges(7,2) * t48 + mrSges(7,3) * t149;
t33 = -mrSges(6,2) * t149 - mrSges(6,3) * t48;
t195 = t30 + t33;
t31 = mrSges(6,1) * t149 - mrSges(6,3) * t47;
t32 = -mrSges(7,1) * t149 + t47 * mrSges(7,2);
t194 = -t31 + t32;
t79 = t146 * qJD(1);
t35 = t102 * t79 + t104 * t178;
t173 = qJD(4) * mrSges(5,1);
t191 = -mrSges(6,1) * t76 - mrSges(6,2) * t77 - mrSges(5,3) * t168 + t173;
t170 = t103 * t104;
t51 = t100 * t170 + t102 * t81;
t165 = qJD(4) * t103;
t59 = -t105 * t160 + t165 * t88;
t181 = t100 * t59;
t180 = t104 * t79;
t179 = t104 * t81;
t174 = qJD(3) * t89;
t159 = qJD(1) * qJD(3);
t37 = mrSges(7,1) * t76 - mrSges(7,3) * t77;
t157 = t37 - t191;
t151 = -t177 / 0.2e1;
t148 = t102 * t75 + t227 * t104 + t81 * t162;
t147 = m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3);
t143 = -t1 * t102 + t104 * t2;
t141 = -t102 * t3 - t104 * t4;
t140 = mrSges(6,1) * t104 - mrSges(6,2) * t102;
t138 = mrSges(7,1) * t104 + mrSges(7,3) * t102;
t131 = Ifges(6,2) * t104 + t189;
t125 = -Ifges(7,3) * t104 + t185;
t124 = pkin(5) * t104 + qJ(6) * t102;
t121 = -t102 * t17 - t104 * t18;
t119 = t102 * t20 - t104 * t21;
t117 = -t100 + t123;
t114 = t102 * t194 + t104 * t195;
t111 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t106 = qJD(1) ^ 2;
t92 = Ifges(7,2) * t149;
t91 = Ifges(6,3) * t149;
t82 = -pkin(4) - t124;
t78 = (mrSges(5,1) * t103 + mrSges(5,2) * t105) * qJD(1);
t58 = t117 * t105;
t50 = -t100 * t171 + t179;
t46 = Ifges(7,4) * t47;
t45 = Ifges(6,5) * t47;
t44 = Ifges(6,6) * t48;
t43 = Ifges(7,6) * t48;
t41 = -t179 + (t100 * t102 - pkin(5)) * t103;
t40 = qJ(6) * t103 + t51;
t36 = pkin(5) * t77 + qJ(6) * t76;
t34 = -t102 * t178 + t180;
t29 = -t180 + (-pkin(5) * qJD(1) + t102 * t88) * t105;
t28 = qJ(6) * t168 + t35;
t16 = mrSges(6,1) * t48 + mrSges(6,2) * t47;
t15 = mrSges(7,1) * t48 - mrSges(7,3) * t47;
t14 = -t117 * t165 + (qJD(5) * t124 - qJD(6) * t104 - qJD(2)) * t105;
t11 = t47 * Ifges(6,4) - t48 * Ifges(6,2) + Ifges(6,6) * t149;
t10 = t47 * Ifges(7,5) + Ifges(7,6) * t149 + t48 * Ifges(7,3);
t8 = -t102 * t152 + t148;
t7 = -pkin(5) * t164 - t9;
t6 = qJ(6) * t164 + (-t100 * t163 + qJD(6)) * t103 + t148;
t5 = pkin(5) * t48 - qJ(6) * t47 - qJD(6) * t77 + t59;
t12 = [qJD(3) * t78 + t14 * t37 + t58 * t15 + t40 * t30 + t50 * t31 + t41 * t32 + t51 * t33 + t8 * t52 + t9 * t53 + t7 * t54 + t6 * t55 + 0.2e1 * (qJD(3) * mrSges(4,3) + qJD(2) * t147) * qJD(1) + m(4) * (qJD(2) * t95 + t174 + (qJ(2) * qJD(2) + qJD(3) * t101) * qJD(1)) + m(6) * (t20 * t9 + t21 * t8 + t3 * t51 + t4 * t50) + m(5) * (qJD(2) * t221 + t101 * t159 + t174) + m(7) * (t1 * t40 + t14 * t19 + t17 * t7 + t18 * t6 + t2 * t41 + t5 * t58) + (t91 / 0.2e1 + t92 / 0.2e1 - t44 / 0.2e1 + t45 / 0.2e1 + t46 / 0.2e1 + t43 / 0.2e1 + mrSges(5,1) * t159 + qJD(2) * t84 + (m(5) * t100 - mrSges(5,3)) * t60 + t154 * t48 + t156 * t47 + (0.3e1 / 0.2e1 * Ifges(5,4) * t169 + (m(6) * t69 - t191) * t100 + (-t89 - t222) * mrSges(5,2) + t151 - t235) * qJD(4) + t111) * t103 + (t126 * t215 + t132 * t216 + t5 * t137 + t10 * t207 - t100 * t16 + (mrSges(5,3) + t139) * t59 + t191 * qJD(2) + t141 * mrSges(6,3) + t143 * mrSges(7,2) - m(5) * t181 + m(6) * (-qJD(2) * t69 - t181) + (t100 * t84 + t150 - t219) * qJD(4) + (mrSges(7,2) * t121 + mrSges(6,3) * t119 + t125 * t214 + t131 * t213 + t138 * t19 + t140 * t69 + t206 * t25 + t229 * t210 + t228 * t212) * qJD(5) + (qJD(3) * mrSges(5,2) + (t101 * mrSges(5,1) + (t187 / 0.2e1 + t182 / 0.2e1 + t186 / 0.2e1 - t183 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4)) * t105 + (0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,1) + t155) * t103) * qJD(4)) * qJD(1) + t223 * t217 + (qJD(5) * t226 + t11) * t208 + (qJD(5) * t22 + t232) * t205) * t105; t194 * t104 - t195 * t102 + t113 * qJD(5) + m(6) * (t119 * qJD(5) + t141) + m(7) * (t121 * qJD(5) + t143) - t147 * t106 + ((-t95 - qJD(3)) * m(4) + (t157 - t173) * t105 + (t172 - t230) * t103 - m(6) * (-t105 * t69 + t170 * t21 - t171 * t20) - m(7) * (-t105 * t19 + t17 * t171 + t170 * t18) + (-qJD(3) - t221) * m(5)) * qJD(1); -t106 * mrSges(4,3) + (-t15 - t16 + t230 * qJD(4) + m(6) * (t161 * t21 - t166 * t20 - t59) + m(7) * (t161 * t18 + t166 * t17 - t5) - m(5) * t59) * t105 + (t157 * qJD(4) + t112 * qJD(5) + m(6) * (qJD(4) * t69 - t162 * t20 - t163 * t21 + t142) + m(7) * (qJD(4) * t19 + t162 * t17 - t163 * t18 + t144) + m(5) * t60 + t114) * t103 + (-m(5) * t89 - t78 + (-t89 + qJD(2)) * m(4) + t220) * qJD(1); ((t150 + Ifges(5,4) * t168 / 0.2e1 + t229 * qJD(4) / 0.2e1 + t219) * t105 + (t89 * mrSges(5,2) + (-t190 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t105) * qJD(1) + t151 + t235) * t103) * qJD(1) + (pkin(8) * t220 + t237) * qJD(5) + (t103 * t191 - t105 * t84) * t88 + t82 * t15 - t60 * mrSges(5,2) - t35 * t52 - t34 * t53 - t29 * t54 - t28 * t55 - pkin(4) * t16 + t231 * t37 - t5 * t138 + (-mrSges(5,1) - t140) * t59 + t142 * mrSges(6,3) + t144 * mrSges(7,2) + t114 * pkin(8) + t11 * t205 + t10 * t206 + t125 * t215 + t131 * t216 + t228 * t217 + t232 * t207 + (pkin(8) * t144 - t17 * t29 - t18 * t28 + t231 * t19 + t5 * t82) * m(7) + (-pkin(4) * t59 + pkin(8) * t142 - t20 * t34 - t21 * t35 - t69 * t80) * m(6); (t192 + t203) * t21 + (-t193 - t204) * t20 + t91 + t92 - t44 + t45 + t46 + t43 - t19 * (mrSges(7,1) * t77 + mrSges(7,3) * t76) - t69 * (mrSges(6,1) * t77 - mrSges(6,2) * t76) + qJD(6) * t55 - t36 * t37 + qJ(6) * t30 - pkin(5) * t32 + (t17 * t76 + t18 * t77) * mrSges(7,2) + t111 + t25 * t211 + (Ifges(7,3) * t77 - t201) * t214 + (-t233 * t76 - t238 * t77) * t210 + (-pkin(5) * t2 + qJ(6) * t1 - t17 * t21 + t18 * t225 - t19 * t36) * m(7) + (-Ifges(6,2) * t77 + t226 - t73) * t213 + (-t236 * t76 - t202 + t22 + t72) * t212; t77 * t37 - t93 * t55 + 0.2e1 * (t2 / 0.2e1 + t18 * t210 + t19 * t211) * m(7) + t32;];
tauc  = t12(:);
