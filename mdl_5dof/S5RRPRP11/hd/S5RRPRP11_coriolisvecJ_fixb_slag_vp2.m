% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:23
% EndTime: 2019-12-31 20:12:33
% DurationCPUTime: 4.52s
% Computational Cost: add. (2147->398), mult. (5138->513), div. (0->0), fcn. (2579->4), ass. (0->188)
t237 = Ifges(5,1) + Ifges(6,1);
t234 = Ifges(6,4) + Ifges(5,5);
t121 = sin(qJ(4));
t123 = cos(qJ(4));
t124 = cos(qJ(2));
t125 = -pkin(2) - pkin(7);
t122 = sin(qJ(2));
t163 = -qJ(3) * t122 - pkin(1);
t80 = t125 * t124 + t163;
t63 = t80 * qJD(1);
t182 = qJD(1) * t122;
t112 = pkin(6) * t182;
t89 = -pkin(3) * t182 - t112;
t64 = t125 * qJD(2) + qJD(3) - t89;
t19 = -t121 * t63 + t123 * t64;
t20 = t121 * t64 + t123 * t63;
t140 = t19 * t121 - t20 * t123;
t110 = qJD(4) + t182;
t224 = qJD(5) - t19;
t14 = -pkin(4) * t110 + t224;
t15 = qJ(5) * t110 + t20;
t141 = t14 * t121 + t15 * t123;
t156 = mrSges(6,1) * t123 + mrSges(6,3) * t121;
t158 = mrSges(5,1) * t123 - mrSges(5,2) * t121;
t203 = t123 / 0.2e1;
t204 = -t123 / 0.2e1;
t206 = -t121 / 0.2e1;
t120 = qJD(2) * qJ(3);
t181 = qJD(1) * t124;
t114 = pkin(6) * t181;
t90 = pkin(3) * t181 + t114;
t75 = t120 + t90;
t84 = qJD(2) * t121 + t123 * t181;
t166 = t121 * t181;
t85 = qJD(2) * t123 - t166;
t21 = pkin(4) * t84 - qJ(5) * t85 + t75;
t200 = Ifges(6,5) * t84;
t82 = Ifges(5,4) * t84;
t225 = t234 * t110 + t237 * t85 + t200 - t82;
t81 = Ifges(6,5) * t85;
t24 = t110 * Ifges(6,6) + t84 * Ifges(6,3) + t81;
t201 = Ifges(5,4) * t85;
t27 = -t84 * Ifges(5,2) + t110 * Ifges(5,6) + t201;
t238 = -t141 * mrSges(6,2) + t140 * mrSges(5,3) + t21 * t156 + t75 * t158 + t24 * t203 + t27 * t204 + t225 * t206;
t191 = Ifges(6,5) * t121;
t145 = -Ifges(6,3) * t123 + t191;
t195 = Ifges(5,4) * t121;
t149 = Ifges(5,2) * t123 + t195;
t207 = t110 / 0.2e1;
t212 = t85 / 0.2e1;
t214 = t84 / 0.2e1;
t215 = -t84 / 0.2e1;
t190 = Ifges(6,5) * t123;
t194 = Ifges(5,4) * t123;
t221 = t237 * t121 - t190 + t194;
t187 = Ifges(6,6) * t123;
t188 = Ifges(5,6) * t123;
t192 = Ifges(5,5) * t121;
t193 = Ifges(6,4) * t121;
t222 = -t187 + t193 + t188 + t192;
t228 = qJD(2) / 0.2e1;
t229 = -qJD(2) / 0.2e1;
t230 = -qJD(1) / 0.2e1;
t96 = -pkin(2) * t124 + t163;
t76 = t96 * qJD(1);
t98 = -t114 - t120;
t236 = t98 * mrSges(4,1) + Ifges(3,6) * t229 + (Ifges(3,4) * t122 + t124 * Ifges(3,2)) * t230 + Ifges(4,5) * t228 + (-Ifges(4,6) * t122 - t124 * Ifges(4,3)) * qJD(1) / 0.2e1 - t76 * mrSges(4,2) + t149 * t215 + t145 * t214 + t221 * t212 + t222 * t207 - t238;
t211 = pkin(3) + pkin(6);
t235 = -mrSges(3,1) + mrSges(4,2);
t175 = qJD(1) * qJD(2);
t164 = t124 * t175;
t174 = qJD(2) * qJD(4);
t176 = qJD(4) * t123;
t180 = qJD(2) * t122;
t52 = -t121 * t174 + (t121 * t180 - t124 * t176) * qJD(1);
t165 = t122 * t175;
t53 = -qJD(4) * t166 + (-t165 + t174) * t123;
t233 = (-Ifges(5,4) + Ifges(6,5)) * t53 + t237 * t52 + t234 * t164;
t232 = t237 * t123 + t191 - t195;
t226 = Ifges(5,6) - Ifges(6,6);
t109 = pkin(2) * t165;
t142 = pkin(7) * t122 - qJ(3) * t124;
t178 = qJD(3) * t122;
t129 = t142 * qJD(2) - t178;
t45 = t129 * qJD(1) + t109;
t179 = qJD(2) * t124;
t92 = t211 * t179;
t79 = qJD(1) * t92;
t4 = -t20 * qJD(4) - t121 * t45 + t123 * t79;
t102 = t211 * t122;
t196 = t121 * t102 + t123 * t80;
t116 = pkin(2) * t180;
t60 = t116 + t129;
t9 = -qJD(4) * t196 - t121 * t60 + t123 * t92;
t111 = Ifges(3,4) * t181;
t168 = Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t169 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t170 = -Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t189 = Ifges(4,6) * t124;
t94 = -qJD(2) * pkin(2) + qJD(3) + t112;
t219 = t168 * t110 + t169 * t84 - t170 * t85 + t15 * mrSges(6,3) + t19 * mrSges(5,1) + t94 * mrSges(4,1) + Ifges(5,6) * t215 + Ifges(6,6) * t214 + Ifges(3,1) * t182 / 0.2e1 + Ifges(3,5) * t228 + t111 / 0.2e1 + Ifges(4,4) * t229 + (-t122 * Ifges(4,2) - t189) * t230 - t14 * mrSges(6,1) - t20 * mrSges(5,2) - t76 * mrSges(4,3) + t234 * t212 + (Ifges(5,3) + Ifges(6,2)) * t207;
t217 = -t53 / 0.2e1;
t216 = t53 / 0.2e1;
t213 = -t85 / 0.2e1;
t210 = pkin(1) * mrSges(3,1);
t209 = pkin(1) * mrSges(3,2);
t208 = -t110 / 0.2e1;
t205 = t121 / 0.2e1;
t202 = mrSges(5,3) * t84;
t199 = t85 * mrSges(5,3);
t55 = -mrSges(5,2) * t110 - t202;
t56 = -t84 * mrSges(6,2) + mrSges(6,3) * t110;
t198 = t55 + t56;
t57 = mrSges(5,1) * t110 - t199;
t58 = -mrSges(6,1) * t110 + t85 * mrSges(6,2);
t197 = -t57 + t58;
t113 = pkin(2) * t182;
t68 = t142 * qJD(1) + t113;
t31 = t121 * t90 + t123 * t68;
t100 = -mrSges(4,1) * t181 - qJD(2) * mrSges(4,3);
t43 = mrSges(5,1) * t84 + mrSges(5,2) * t85;
t186 = -t100 + t43;
t185 = qJD(2) * mrSges(3,2);
t184 = t121 * t125;
t183 = t123 * t125;
t103 = t211 * t124;
t177 = qJD(4) * t121;
t173 = -0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(4,6);
t172 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t171 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t167 = m(4) * pkin(6) + mrSges(4,1);
t162 = m(4) * t98 - mrSges(3,3) * t181 + t100 + t185;
t161 = m(4) * t94 + (mrSges(4,1) + mrSges(3,3)) * t182 + t235 * qJD(2);
t91 = t211 * t180;
t3 = t121 * t79 + t123 * t45 + t64 * t176 - t63 * t177;
t1 = qJ(5) * t164 + qJD(5) * t110 + t3;
t2 = -pkin(4) * t164 - t4;
t160 = -t1 * t121 + t123 * t2;
t159 = -t121 * t3 - t123 * t4;
t157 = mrSges(5,1) * t121 + mrSges(5,2) * t123;
t155 = mrSges(6,1) * t121 - mrSges(6,3) * t123;
t150 = -Ifges(5,2) * t121 + t194;
t146 = Ifges(6,3) * t121 + t190;
t144 = pkin(4) * t123 + qJ(5) * t121;
t143 = -pkin(4) * t121 + qJ(5) * t123;
t30 = -t121 * t68 + t123 * t90;
t39 = t102 * t123 - t121 * t80;
t36 = -mrSges(6,1) * t164 + t52 * mrSges(6,2);
t136 = -pkin(3) - t144;
t8 = t102 * t176 + t121 * t92 + t123 * t60 - t80 * t177;
t133 = -qJ(3) * t179 - t178;
t119 = qJD(2) * qJD(3);
t66 = -qJD(1) * t91 + t119;
t34 = -mrSges(6,2) * t53 + mrSges(6,3) * t164;
t35 = mrSges(5,1) * t164 - mrSges(5,3) * t52;
t37 = -mrSges(5,2) * t164 - mrSges(5,3) * t53;
t132 = (t35 - t36) * t123 + (t34 + t37) * t121;
t131 = t197 * t121 + t198 * t123;
t130 = t4 * mrSges(5,1) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + t1 * mrSges(6,3);
t108 = Ifges(6,2) * t164;
t107 = Ifges(5,3) * t164;
t95 = qJ(3) - t143;
t93 = pkin(6) * t165 - t119;
t87 = (mrSges(4,2) * t124 - mrSges(4,3) * t122) * qJD(1);
t70 = t116 + t133;
t65 = t144 * qJD(4) - qJD(5) * t123 + qJD(3);
t61 = t133 * qJD(1) + t109;
t59 = t144 * t124 + t103;
t50 = Ifges(6,4) * t52;
t49 = Ifges(5,5) * t52;
t48 = Ifges(5,6) * t53;
t47 = Ifges(6,6) * t53;
t44 = t136 * t182 - t112;
t42 = mrSges(6,1) * t84 - mrSges(6,3) * t85;
t41 = pkin(4) * t85 + qJ(5) * t84;
t33 = -pkin(4) * t122 - t39;
t32 = qJ(5) * t122 + t196;
t23 = -pkin(4) * t181 - t30;
t22 = qJ(5) * t181 + t31;
t18 = mrSges(5,1) * t53 + mrSges(5,2) * t52;
t17 = mrSges(6,1) * t53 - mrSges(6,3) * t52;
t16 = (t143 * qJD(4) + qJD(5) * t121) * t124 + (-pkin(6) + t136) * t180;
t11 = t52 * Ifges(5,4) - t53 * Ifges(5,2) + Ifges(5,6) * t164;
t10 = t52 * Ifges(6,5) + Ifges(6,6) * t164 + t53 * Ifges(6,3);
t7 = -pkin(4) * t179 - t9;
t6 = pkin(4) * t53 - qJ(5) * t52 - qJD(5) * t85 + t66;
t5 = qJ(5) * t179 + qJD(5) * t122 + t8;
t12 = [t103 * t18 + t16 * t42 + t59 * t17 + t32 * t34 + t33 * t36 + t39 * t35 + t196 * t37 - t91 * t43 + t5 * t56 + t8 * t55 + t9 * t57 + t7 * t58 + t70 * t87 + m(5) * (t103 * t66 + t19 * t9 + t196 * t3 + t20 * t8 + t39 * t4 - t75 * t91) + m(6) * (t1 * t32 + t14 * t7 + t15 * t5 + t16 * t21 + t2 * t33 + t59 * t6) + m(4) * (t61 * t96 + t70 * t76) + (-t61 * mrSges(4,3) + t107 / 0.2e1 + t108 / 0.2e1 - t48 / 0.2e1 + t49 / 0.2e1 + t50 / 0.2e1 + t47 / 0.2e1 + t169 * t53 - t170 * t52 + (t162 * pkin(6) + (-t96 * mrSges(4,2) + t173 * t122 - 0.2e1 * t210) * qJD(1) + t171 * qJD(2) + t236) * qJD(2) + t130) * t122 + (t10 * t203 + t11 * t204 + t149 * t216 + t145 * t217 + t61 * mrSges(4,2) + t66 * t158 + t6 * t156 - t167 * t93 + (t121 * t4 - t123 * t3) * mrSges(5,3) + (-t1 * t123 - t121 * t2) * mrSges(6,2) + (t150 * t214 + t146 * t215 - t75 * t157 - t21 * t155 + t27 * t205 + (t121 * t20 + t123 * t19) * mrSges(5,3) + (t121 * t15 - t123 * t14) * mrSges(6,2) + t232 * t213 + (t226 * t121 - t123 * t234) * t207 + t225 * t204) * qJD(4) + (t161 * pkin(6) + ((-t192 / 0.2e1 - t188 / 0.2e1 - t193 / 0.2e1 + t187 / 0.2e1 - t173) * t124 - t96 * mrSges(4,3) - 0.2e1 * t209 + (0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + t167 * pkin(6) + t168) * t122) * qJD(1) + t172 * qJD(2) + t219) * qJD(2) - t221 * t52 / 0.2e1 + (t24 * qJD(4) + t233) * t206) * t124; t132 * t125 + ((((-m(4) * pkin(2) + t235) * qJD(2) - t161) * pkin(6) + (-pkin(2) * mrSges(4,1) + t169 * t121 - t170 * t123 + t172) * qJD(2) + (t209 - t189 / 0.2e1) * qJD(1) - t111 / 0.2e1 - t219) * t124 + ((-t162 + t185) * pkin(6) + (-qJ(3) * mrSges(4,1) + t171) * qJD(2) + (Ifges(3,2) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1) * t181 + (t210 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t122) * qJD(1) - t236) * t122) * qJD(1) - m(5) * (t19 * t30 + t20 * t31 + t75 * t89) - m(6) * (t14 * t23 + t15 * t22 + t21 * t44) + (-m(4) * t76 - t87) * (-qJ(3) * t181 + t113) + m(4) * (-qJ(3) * t93 - qJD(3) * t98) - t89 * t43 - t93 * mrSges(4,3) + t95 * t17 - t31 * t55 - t22 * t56 - t30 * t57 - t23 * t58 + t159 * mrSges(5,3) + t160 * mrSges(6,2) + t6 * t155 + t66 * t157 + (t149 * t214 + t145 * t215 + (-m(5) * t140 + m(6) * t141 + t131) * t125 + t221 * t213 + t222 * t208 + t238) * qJD(4) + qJ(3) * t18 + m(5) * (t66 * qJ(3) + t75 * qJD(3) + t4 * t183 + t3 * t184) + m(6) * (t1 * t184 - t2 * t183 + t21 * t65 + t6 * t95) + t186 * qJD(3) + (t65 - t44) * t42 + t233 * t203 + t232 * t52 / 0.2e1 + t10 * t205 + t11 * t206 + t146 * t216 + t150 * t217; (-t42 - t186) * qJD(2) + t131 * qJD(4) + (t167 * t179 + (t87 + t131) * t122) * qJD(1) - m(4) * (-qJD(2) * t98 - t76 * t182) + t132 + (-qJD(2) * t21 + t110 * t141 - t160) * m(6) + (-qJD(2) * t75 - t110 * t140 - t159) * m(5); t130 + t107 + t108 - t21 * (mrSges(6,1) * t85 + mrSges(6,3) * t84) - t75 * (mrSges(5,1) * t85 - mrSges(5,2) * t84) + qJD(5) * t56 - t41 * t42 - t48 + t49 + t50 + t47 + qJ(5) * t34 - pkin(4) * t36 + (-t198 - t202) * t19 + (-t197 + t199) * t20 + (t14 * t84 + t15 * t85) * mrSges(6,2) + t27 * t212 + (Ifges(6,3) * t85 - t200) * t215 + (-t226 * t85 - t234 * t84) * t208 + (-pkin(4) * t2 + qJ(5) * t1 - t14 * t20 + t224 * t15 - t21 * t41) * m(6) + (-Ifges(5,2) * t85 + t225 - t82) * t214 + (-t237 * t84 - t201 + t24 + t81) * t213; -t110 * t56 + t85 * t42 + 0.2e1 * (t2 / 0.2e1 + t15 * t208 + t21 * t212) * m(6) + t36;];
tauc = t12(:);
