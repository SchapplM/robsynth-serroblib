% Calculate vector of centrifugal and coriolis load on the joints for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PPRRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:50:17
% EndTime: 2018-11-23 14:50:21
% DurationCPUTime: 4.04s
% Computational Cost: add. (3284->380), mult. (8787->527), div. (0->0), fcn. (6835->12), ass. (0->201)
t128 = qJD(3) ^ 2;
t122 = sin(qJ(4));
t125 = cos(qJ(4));
t179 = qJD(3) * qJD(4);
t84 = (-mrSges(6,2) * t122 - mrSges(6,3) * t125) * t179;
t85 = (mrSges(5,1) * t122 + mrSges(5,2) * t125) * t179;
t262 = -mrSges(4,2) * t128 - t84 - t85;
t148 = pkin(10) * t122 - qJ(5) * t125;
t139 = t148 * qJD(4);
t181 = qJD(5) * t122;
t184 = qJD(4) * t122;
t162 = pkin(4) * t184 - t181;
t120 = cos(pkin(6));
t106 = qJD(1) * t120 + qJD(2);
t115 = sin(pkin(12));
t117 = sin(pkin(6));
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t118 = cos(pkin(12));
t119 = cos(pkin(7));
t190 = t118 * t119;
t135 = (t115 * t126 + t123 * t190) * t117;
t116 = sin(pkin(7));
t192 = t116 * t123;
t48 = qJD(1) * t135 + t106 * t192;
t261 = t48 - t139 - t162;
t225 = pkin(5) + pkin(9);
t103 = t225 * t125;
t189 = qJD(1) * t117;
t169 = t118 * t189;
t193 = t115 * t123;
t47 = t126 * (t106 * t116 + t119 * t169) - t189 * t193;
t260 = qJD(4) * t103 - t122 * t47;
t243 = qJD(3) / 0.2e1;
t240 = mrSges(6,1) + mrSges(5,3);
t259 = mrSges(6,2) - mrSges(5,1);
t124 = cos(qJ(6));
t121 = sin(qJ(6));
t209 = Ifges(7,4) * t121;
t150 = Ifges(7,2) * t124 + t209;
t208 = Ifges(7,4) * t124;
t152 = Ifges(7,1) * t121 + t208;
t155 = mrSges(7,1) * t124 - mrSges(7,2) * t121;
t127 = -pkin(4) - pkin(10);
t45 = qJD(3) * pkin(9) + t48;
t163 = pkin(5) * qJD(3) + t45;
t156 = t163 * t122;
t73 = t106 * t119 - t116 * t169;
t67 = t125 * t73;
t22 = t67 - t156;
t252 = qJD(5) - t22;
t20 = qJD(4) * t127 + t252;
t164 = -qJ(5) * t122 - pkin(3);
t86 = t125 * t127 + t164;
t36 = qJD(3) * t86 - t47;
t5 = -t121 * t36 + t124 * t20;
t6 = t121 * t20 + t124 * t36;
t157 = t5 * t121 - t6 * t124;
t204 = Ifges(7,6) * t124;
t207 = Ifges(7,5) * t121;
t186 = qJD(3) * t125;
t113 = pkin(5) * t186;
t199 = t122 * t73;
t31 = t125 * t45 + t199;
t27 = -qJD(4) * qJ(5) - t31;
t21 = t113 - t27;
t220 = -t124 / 0.2e1;
t221 = -t121 / 0.2e1;
t188 = qJD(3) * t122;
t109 = qJD(6) + t188;
t222 = -t109 / 0.2e1;
t183 = qJD(4) * t124;
t89 = -t121 * t186 + t183;
t226 = -t89 / 0.2e1;
t88 = -qJD(4) * t121 - t124 * t186;
t227 = -t88 / 0.2e1;
t219 = Ifges(7,4) * t89;
t50 = t88 * Ifges(7,2) + t109 * Ifges(7,6) + t219;
t87 = Ifges(7,4) * t88;
t51 = t89 * Ifges(7,1) + t109 * Ifges(7,5) + t87;
t258 = t157 * mrSges(7,3) + (t204 + t207) * t222 + t150 * t227 + t152 * t226 + t21 * t155 + t220 * t50 + t221 * t51;
t100 = -mrSges(6,1) * t186 - qJD(4) * mrSges(6,3);
t195 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t186 - t100;
t59 = -mrSges(7,1) * t88 + mrSges(7,2) * t89;
t176 = -t59 - t195;
t234 = -m(5) * t31 + m(6) * t27;
t257 = -m(7) * t21 + t176 + t234;
t256 = -t186 / 0.2e1;
t230 = 0.2e1 * m(5);
t255 = -t230 / 0.2e1;
t102 = t225 * t122;
t57 = t102 * t124 - t121 * t86;
t254 = qJD(6) * t57 + t260 * t121 - t261 * t124;
t58 = t102 * t121 + t124 * t86;
t253 = -qJD(6) * t58 + t261 * t121 + t260 * t124;
t241 = qJD(4) / 0.2e1;
t242 = -qJD(4) / 0.2e1;
t96 = -pkin(4) * t125 + t164;
t39 = qJD(3) * t96 - t47;
t44 = -qJD(3) * pkin(3) - t47;
t251 = t44 * mrSges(5,1) + t27 * mrSges(6,1) - t39 * mrSges(6,2) - t31 * mrSges(5,3) - Ifges(6,5) * t242 - Ifges(5,6) * t241 - t258 + ((-Ifges(5,2) - Ifges(6,3)) * t125 + (-Ifges(5,4) - Ifges(6,6)) * t122) * t243;
t136 = (t126 * t190 - t193) * t117;
t191 = t116 * t126;
t133 = t120 * t191 + t136;
t196 = -t259 * qJD(4) - t240 * t188;
t201 = t122 * t45;
t30 = -t67 + t201;
t237 = -qJD(5) - t30;
t24 = -qJD(4) * pkin(4) - t237;
t249 = m(5) * t30 + m(6) * t24 - t196;
t90 = (mrSges(6,2) * t125 - mrSges(6,3) * t122) * qJD(3);
t210 = t90 + (-mrSges(5,1) * t125 + mrSges(5,2) * t122) * qJD(3);
t248 = -m(6) * t39 + qJD(3) * mrSges(4,1) - t210;
t247 = -Ifges(5,1) / 0.2e1;
t244 = Ifges(5,4) * t256;
t165 = t122 * t179;
t108 = pkin(4) * t165;
t131 = t48 - t181;
t32 = t108 + (t139 + t131) * qJD(3);
t42 = (qJD(1) * t136 + t106 * t191) * qJD(3);
t202 = t122 * t42;
t8 = t202 + (t125 * t163 + t199) * qJD(4);
t1 = qJD(6) * t5 + t121 * t8 + t124 * t32;
t2 = -qJD(6) * t6 - t121 * t32 + t124 * t8;
t235 = t1 * t121 + t124 * t2;
t178 = qJD(4) * qJD(6);
t180 = qJD(6) * t125;
t68 = -t121 * t178 + (t121 * t184 - t124 * t180) * qJD(3);
t69 = -t124 * t178 + (t121 * t180 + t122 * t183) * qJD(3);
t233 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t68 + Ifges(7,6) * t69;
t229 = -t47 / 0.2e1;
t228 = t50 / 0.2e1;
t11 = qJD(4) * t31 + t202;
t140 = -t116 * t117 * t118 + t119 * t120;
t137 = t125 * t140;
t56 = t120 * t192 + t135;
t37 = t122 * t56 - t137;
t217 = t11 * t37;
t172 = t122 * t192;
t78 = -t125 * t119 + t172;
t216 = t11 * t78;
t43 = qJD(3) * t48;
t214 = t43 * t133;
t182 = qJD(4) * t125;
t211 = t125 * t42 + t73 * t182;
t206 = Ifges(7,5) * t124;
t205 = Ifges(7,6) * t121;
t198 = t126 * t43;
t197 = t100 - t59;
t187 = qJD(3) * t123;
t185 = qJD(3) * t126;
t177 = pkin(9) * t11 / 0.2e1;
t175 = -0.3e1 / 0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * Ifges(6,6);
t174 = -Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t173 = Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1;
t170 = qJ(5) * t182;
t168 = t116 * t187;
t167 = t116 * t185;
t166 = t125 * t179;
t154 = mrSges(7,1) * t121 + mrSges(7,2) * t124;
t153 = Ifges(7,1) * t124 - t209;
t151 = -Ifges(7,2) * t121 + t208;
t16 = t121 * t133 + t124 * t37;
t17 = t121 * t37 - t124 * t133;
t53 = mrSges(7,1) * t166 - mrSges(7,3) * t68;
t54 = -mrSges(7,2) * t166 + mrSges(7,3) * t69;
t147 = t121 * t54 + t124 * t53;
t70 = -mrSges(7,2) * t109 + mrSges(7,3) * t88;
t71 = mrSges(7,1) * t109 - mrSges(7,3) * t89;
t146 = -t121 * t71 + t124 * t70;
t142 = -t121 * t78 + t124 * t191;
t62 = t121 * t191 + t124 * t78;
t79 = t119 * t122 + t125 * t192;
t38 = t122 * t140 + t125 * t56;
t132 = qJD(3) * t133;
t130 = t39 * mrSges(6,3) + t6 * mrSges(7,2) - t109 * Ifges(7,3) - t89 * Ifges(7,5) - t88 * Ifges(7,6) + t188 * t247 + Ifges(5,5) * t242 + t244 + Ifges(6,4) * t241 + (-t122 * Ifges(6,2) - Ifges(6,6) * t125) * t243 - t24 * mrSges(6,1) - t30 * mrSges(5,3) - t44 * mrSges(5,2) - t5 * mrSges(7,1);
t112 = pkin(4) * t188;
t107 = Ifges(7,3) * t166;
t93 = t225 * t184;
t92 = -qJ(5) * t186 + t112;
t77 = t162 - t170;
t76 = qJD(3) * t148 + t112;
t61 = qJD(4) * t79 + t122 * t167;
t60 = qJD(4) * t172 - t119 * t182 - t125 * t167;
t52 = t56 * qJD(3);
t40 = -mrSges(7,1) * t69 + mrSges(7,2) * t68;
t35 = t68 * Ifges(7,1) + t69 * Ifges(7,4) + Ifges(7,5) * t166;
t34 = t68 * Ifges(7,4) + t69 * Ifges(7,2) + Ifges(7,6) * t166;
t33 = t108 + (t131 - t170) * qJD(3);
t29 = qJD(6) * t142 - t121 * t168 + t124 * t61;
t28 = qJD(6) * t62 + t121 * t61 + t124 * t168;
t23 = t113 + t31;
t15 = qJD(4) * t38 + t122 * t132;
t13 = t121 * t23 + t124 * t76;
t12 = -t121 * t76 + t124 * t23;
t10 = -t184 * t45 + t211;
t9 = (-qJD(5) + t201) * qJD(4) - t211;
t7 = (qJD(5) - t156) * qJD(4) + t211;
t4 = -qJD(6) * t17 - t121 * t52 + t124 * t15;
t3 = qJD(6) * t16 + t121 * t15 + t124 * t52;
t14 = [m(6) * (-t38 * t9 + t217) + t38 * t40 + t17 * t54 + t16 * t53 + m(5) * (t10 * t38 - t214 + t217) + t4 * t71 + m(7) * (t1 * t17 + t16 * t2 + t3 * t6 + t38 * t7 + t4 * t5) + m(4) * (t132 * t48 + t42 * t56 - t214) + t3 * t70 + t249 * t15 + (-m(4) * t47 + m(5) * t44 - t248) * t52 + t257 * (-qJD(4) * t137 - t125 * t132 + t184 * t56) + t240 * (-t38 * t165 + t37 * t166) + (-m(6) * t33 + t262) * t133; t28 * t70 + t29 * t71 + t79 * t40 + t62 * t53 - t142 * t54 - t196 * t61 + t176 * t60 + m(5) * (t10 * t79 + t30 * t61 - t31 * t60 + t216) + m(7) * (-t1 * t142 + t2 * t62 - t21 * t60 + t28 * t6 + t29 * t5 + t7 * t79) + m(6) * (t24 * t61 + t27 * t60 - t79 * t9 + t216) + (t262 * t126 + (-t128 * mrSges(4,1) + qJD(3) * t210) * t123 + m(4) * (t123 * t42 + t185 * t48 - t187 * t47 - t198) + m(5) * (t187 * t44 - t198) + m(6) * (-t126 * t33 + t187 * t39)) * t116 + t240 * t179 * (-t122 * t79 + t125 * t78); -t43 * mrSges(4,1) + t103 * t40 + t57 * t53 + t58 * t54 - t93 * t59 + t77 * t90 + t96 * t84 + t253 * t71 + t254 * t70 + (qJD(3) * t47 - t42) * mrSges(4,2) + m(6) * (t33 * t96 + t39 * t77) + (t107 / 0.2e1 + t43 * mrSges(5,2) - t33 * mrSges(6,3) + t196 * t47 + t240 * t11 + 0.2e1 * (t229 * t24 + t177) * m(6) + (t229 * t30 + t177) * t230 + (t173 * qJD(4) + (-t195 + t234) * pkin(9) + t175 * t188 + t251) * qJD(4) + t233) * t122 + (-t68 * t152 / 0.2e1 + t7 * t155 - t69 * t150 / 0.2e1 - t43 * mrSges(5,1) + t33 * mrSges(6,2) + t10 * mrSges(5,3) - t9 * mrSges(6,1) + t35 * t221 + t34 * t220 + (-t1 * t124 + t2 * t121) * mrSges(7,3) + (t51 * t220 + t121 * t228 - t21 * t154 + t109 * (t205 - t206) / 0.2e1 + t151 * t227 + t153 * t226 + (t6 * t121 + t5 * t124) * mrSges(7,3)) * qJD(6) + (-t130 + t174 * qJD(4) + ((-t207 / 0.2e1 - t204 / 0.2e1 - t175) * t125 + (0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(6,2) + Ifges(7,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(6,3)) * t122) * qJD(3)) * qJD(4) + (m(5) * t10 - m(6) * t9 + qJD(4) * t249) * pkin(9) + t257 * t47) * t125 + (t255 * t43 - t85) * pkin(3) + (t255 * t44 + t248) * t48 + (t1 * t58 + t103 * t7 + t2 * t57 - t21 * t93 + t253 * t5 + t254 * t6) * m(7); t68 * t153 / 0.2e1 + t7 * t154 + t69 * t151 / 0.2e1 + t124 * t35 / 0.2e1 + t34 * t221 - t92 * t90 - t13 * t70 - t12 * t71 - t22 * t59 + qJ(5) * t40 - t9 * mrSges(6,3) - t10 * mrSges(5,2) + t196 * t31 + t195 * t30 + t259 * t11 - t197 * qJD(5) - t235 * mrSges(7,3) + t258 * qJD(6) + ((t244 + t130 + Ifges(6,6) * t256 + (t206 / 0.2e1 - t205 / 0.2e1 - pkin(4) * mrSges(6,1) + t174) * qJD(4)) * t125 + ((Ifges(5,2) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + t247) * t186 + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t188 + (-qJ(5) * mrSges(6,1) + t173) * qJD(4) - t251) * t122) * qJD(3) + (qJ(5) * t7 - t12 * t5 - t13 * t6 + t252 * t21) * m(7) + (t147 + m(7) * t235 + (-m(7) * t157 + t146) * qJD(6)) * t127 + (-pkin(4) * t11 - qJ(5) * t9 + t237 * t27 - t24 * t31 - t39 * t92) * m(6); t146 * qJD(6) + t197 * qJD(4) + (mrSges(6,1) * t182 + (t146 + t90) * t122) * qJD(3) + t147 + (-t21 * qJD(4) - t109 * t157 + t235) * m(7) + (qJD(4) * t27 + t188 * t39 + t11) * m(6); t107 - t21 * (mrSges(7,1) * t89 + mrSges(7,2) * t88) + (Ifges(7,1) * t88 - t219) * t226 + t89 * t228 + (Ifges(7,5) * t88 - Ifges(7,6) * t89) * t222 - t5 * t70 + t6 * t71 + (t5 * t88 + t6 * t89) * mrSges(7,3) + (-Ifges(7,2) * t89 + t51 + t87) * t227 + t233;];
tauc  = t14(:);
