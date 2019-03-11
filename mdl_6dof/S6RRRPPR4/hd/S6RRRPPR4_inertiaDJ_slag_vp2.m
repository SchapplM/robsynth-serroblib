% Calculate time derivative of joint inertia matrix for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:53
% EndTime: 2019-03-09 15:32:04
% DurationCPUTime: 5.12s
% Computational Cost: add. (5255->523), mult. (12654->749), div. (0->0), fcn. (10939->8), ass. (0->213)
t267 = Ifges(6,4) + Ifges(5,5);
t266 = Ifges(6,6) - Ifges(5,6);
t268 = -Ifges(6,2) - Ifges(4,3) - Ifges(5,3);
t196 = sin(qJ(2));
t199 = cos(qJ(2));
t173 = -pkin(2) * t199 - t196 * pkin(8) - pkin(1);
t198 = cos(qJ(3));
t234 = t198 * t199;
t181 = pkin(7) * t234;
t195 = sin(qJ(3));
t130 = t195 * t173 + t181;
t227 = qJD(3) * t198;
t230 = qJD(2) * t199;
t203 = t195 * t230 + t196 * t227;
t228 = qJD(3) * t196;
t220 = t195 * t228;
t221 = t198 * t230;
t265 = t220 - t221;
t194 = sin(qJ(6));
t197 = cos(qJ(6));
t264 = -mrSges(7,1) * t194 - mrSges(7,2) * t197;
t193 = sin(pkin(10));
t252 = pkin(3) * t193;
t182 = qJ(5) + t252;
t263 = m(6) * t182 + mrSges(6,3);
t231 = qJD(2) * t196;
t239 = cos(pkin(10));
t214 = t239 * t198;
t215 = t239 * t195;
t87 = t193 * t265 - t214 * t228 - t215 * t230;
t157 = t193 * t198 + t215;
t238 = t193 * t195;
t205 = t214 - t238;
t88 = -t157 * t228 + t205 * t230;
t262 = -Ifges(4,5) * t221 + t268 * t231 + t266 * t87 - t267 * t88;
t261 = 2 * m(4);
t260 = 2 * m(5);
t259 = 0.2e1 * m(6);
t258 = 2 * m(7);
t257 = -0.2e1 * pkin(1);
t256 = 0.2e1 * pkin(7);
t255 = m(5) * pkin(3);
t254 = -pkin(4) - pkin(5);
t253 = -t195 / 0.2e1;
t251 = pkin(7) * t195;
t249 = -qJ(4) - pkin(8);
t139 = t157 * t196;
t140 = t205 * t196;
t76 = t139 * t197 - t140 * t194;
t24 = qJD(6) * t76 - t194 * t87 + t197 * t88;
t77 = t139 * t194 + t140 * t197;
t25 = -qJD(6) * t77 - t194 * t88 - t197 * t87;
t248 = Ifges(7,5) * t24 + Ifges(7,6) * t25;
t226 = qJD(4) * t198;
t168 = (pkin(2) * t196 - pkin(8) * t199) * qJD(2);
t232 = t198 * t168 + t231 * t251;
t43 = -t196 * t226 + (pkin(3) * t196 - qJ(4) * t234) * qJD(2) + (-t181 + (qJ(4) * t196 - t173) * t195) * qJD(3) + t232;
t233 = t195 * t168 + t173 * t227;
t236 = t196 * t198;
t55 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t236 + (-qJD(4) * t196 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t199) * t195 + t233;
t16 = t193 * t43 + t239 * t55;
t245 = Ifges(4,4) * t195;
t244 = Ifges(4,4) * t198;
t243 = Ifges(4,6) * t199;
t223 = t239 * pkin(3);
t184 = -t223 - pkin(4);
t179 = -pkin(5) + t184;
t131 = t179 * t197 - t182 * t194;
t117 = qJD(5) * t197 + qJD(6) * t131;
t242 = t117 * mrSges(7,2);
t132 = t179 * t194 + t182 * t197;
t118 = -qJD(5) * t194 - qJD(6) * t132;
t241 = t118 * mrSges(7,1);
t151 = t205 * qJD(3);
t240 = t151 * mrSges(6,2);
t237 = t195 * t196;
t176 = Ifges(4,1) * t195 + t244;
t235 = t198 * t176;
t159 = t198 * t173;
t108 = -qJ(4) * t236 + t159 + (-pkin(3) - t251) * t199;
t119 = -qJ(4) * t237 + t130;
t58 = t193 * t108 + t239 * t119;
t216 = qJD(3) * t249;
t149 = t195 * t216 + t226;
t202 = -qJD(4) * t195 + t198 * t216;
t90 = t239 * t149 + t193 * t202;
t174 = t249 * t198;
t121 = -t239 * t174 + t249 * t238;
t169 = pkin(3) * t237 + t196 * pkin(7);
t229 = qJD(3) * t195;
t225 = pkin(3) * t229;
t224 = pkin(7) * t230;
t186 = -pkin(3) * t198 - pkin(2);
t39 = -t87 * mrSges(5,1) + t88 * mrSges(5,2);
t38 = -t87 * mrSges(6,1) - t88 * mrSges(6,3);
t218 = Ifges(4,6) * t195 + (2 * Ifges(3,4));
t120 = -t174 * t193 - t249 * t215;
t89 = t149 * t193 - t239 * t202;
t217 = t120 * t89 + t121 * t90;
t150 = t157 * qJD(3);
t98 = t150 * mrSges(5,1) + t151 * mrSges(5,2);
t97 = t150 * mrSges(6,1) - t151 * mrSges(6,3);
t69 = -mrSges(6,1) * t231 + t88 * mrSges(6,2);
t5 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t106 = -t157 * t194 - t197 * t205;
t44 = qJD(6) * t106 + t150 * t194 + t151 * t197;
t107 = t157 * t197 - t194 * t205;
t45 = -qJD(6) * t107 + t150 * t197 - t151 * t194;
t17 = -t45 * mrSges(7,1) + t44 * mrSges(7,2);
t49 = -qJ(5) * t199 + t58;
t213 = qJ(5) * t140 - t169;
t212 = t193 * t55 - t239 * t43;
t211 = mrSges(4,1) * t195 + mrSges(4,2) * t198;
t210 = Ifges(4,1) * t198 - t245;
t209 = -Ifges(4,2) * t195 + t244;
t175 = Ifges(4,2) * t198 + t245;
t208 = Ifges(4,5) * t195 + Ifges(4,6) * t198;
t57 = t108 * t239 - t193 * t119;
t51 = t199 * pkin(4) - t57;
t29 = t199 * pkin(5) - t140 * pkin(9) + t51;
t35 = t139 * pkin(9) + t49;
t10 = -t194 * t35 + t197 * t29;
t11 = t194 * t29 + t197 * t35;
t91 = -pkin(9) * t157 + t120;
t92 = -pkin(9) * t205 + t121;
t36 = -t194 * t92 + t197 * t91;
t37 = t194 * t91 + t197 * t92;
t13 = qJ(5) * t231 - qJD(5) * t199 + t16;
t207 = qJ(5) * t157 - t186;
t41 = Ifges(7,6) * t45;
t42 = Ifges(7,5) * t44;
t59 = -pkin(9) * t151 + t89;
t60 = pkin(9) * t150 + t90;
t6 = qJD(6) * t36 + t194 * t59 + t197 * t60;
t7 = -qJD(6) * t37 - t194 * t60 + t197 * t59;
t206 = t7 * mrSges(7,1) - t6 * mrSges(7,2) + t41 + t42;
t61 = t150 * pkin(4) - t151 * qJ(5) - t157 * qJD(5) + t225;
t8 = -t88 * pkin(9) + t231 * t254 + t212;
t9 = -t87 * pkin(9) + t13;
t1 = qJD(6) * t10 + t194 * t8 + t197 * t9;
t2 = -qJD(6) * t11 - t194 * t9 + t197 * t8;
t201 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) + Ifges(7,3) * t231 - t248;
t128 = pkin(3) * t203 + t224;
t26 = -t87 * pkin(4) - t88 * qJ(5) - t140 * qJD(5) + t128;
t190 = Ifges(4,5) * t227;
t167 = -mrSges(4,1) * t199 - mrSges(4,3) * t236;
t166 = mrSges(4,2) * t199 - mrSges(4,3) * t237;
t165 = t210 * qJD(3);
t164 = t209 * qJD(3);
t163 = t211 * qJD(3);
t147 = Ifges(6,4) * t151;
t146 = Ifges(5,5) * t151;
t145 = Ifges(5,6) * t150;
t144 = Ifges(6,6) * t150;
t138 = -Ifges(4,5) * t199 + t196 * t210;
t137 = t196 * t209 - t243;
t129 = -t199 * t251 + t159;
t127 = -mrSges(4,2) * t231 - mrSges(4,3) * t203;
t126 = mrSges(4,1) * t231 + mrSges(4,3) * t265;
t125 = mrSges(6,1) * t199 + t140 * mrSges(6,2);
t124 = -mrSges(5,1) * t199 - t140 * mrSges(5,3);
t123 = mrSges(5,2) * t199 - t139 * mrSges(5,3);
t122 = -t139 * mrSges(6,2) - mrSges(6,3) * t199;
t116 = Ifges(5,1) * t157 + Ifges(5,4) * t205;
t115 = Ifges(6,1) * t157 - Ifges(6,5) * t205;
t114 = Ifges(5,4) * t157 + Ifges(5,2) * t205;
t113 = Ifges(6,5) * t157 - Ifges(6,3) * t205;
t112 = -mrSges(5,1) * t205 + mrSges(5,2) * t157;
t111 = -mrSges(6,1) * t205 - mrSges(6,3) * t157;
t105 = mrSges(4,1) * t203 - mrSges(4,2) * t265;
t103 = -pkin(4) * t205 - t207;
t102 = Ifges(5,1) * t151 - Ifges(5,4) * t150;
t101 = Ifges(6,1) * t151 + Ifges(6,5) * t150;
t100 = Ifges(5,4) * t151 - Ifges(5,2) * t150;
t99 = Ifges(6,5) * t151 + Ifges(6,3) * t150;
t96 = -t176 * t228 + (Ifges(4,5) * t196 + t199 * t210) * qJD(2);
t95 = -t175 * t228 + (Ifges(4,6) * t196 + t199 * t209) * qJD(2);
t94 = mrSges(5,1) * t139 + mrSges(5,2) * t140;
t93 = mrSges(6,1) * t139 - mrSges(6,3) * t140;
t75 = Ifges(5,1) * t140 - Ifges(5,4) * t139 - Ifges(5,5) * t199;
t74 = Ifges(6,1) * t140 - Ifges(6,4) * t199 + Ifges(6,5) * t139;
t73 = Ifges(5,4) * t140 - Ifges(5,2) * t139 - Ifges(5,6) * t199;
t72 = Ifges(6,5) * t140 - Ifges(6,6) * t199 + Ifges(6,3) * t139;
t71 = -qJD(3) * t130 + t232;
t70 = (-t198 * t231 - t199 * t229) * pkin(7) + t233;
t68 = mrSges(5,1) * t231 - mrSges(5,3) * t88;
t67 = -mrSges(5,2) * t231 + mrSges(5,3) * t87;
t66 = mrSges(6,2) * t87 + mrSges(6,3) * t231;
t65 = -t205 * t254 + t207;
t64 = mrSges(7,1) * t199 - t77 * mrSges(7,3);
t63 = -mrSges(7,2) * t199 + t76 * mrSges(7,3);
t62 = pkin(4) * t139 - t213;
t54 = Ifges(7,1) * t107 + Ifges(7,4) * t106;
t53 = Ifges(7,4) * t107 + Ifges(7,2) * t106;
t52 = -mrSges(7,1) * t106 + mrSges(7,2) * t107;
t50 = t139 * t254 + t213;
t48 = pkin(5) * t150 + t61;
t34 = -mrSges(7,1) * t76 + mrSges(7,2) * t77;
t33 = Ifges(5,1) * t88 + Ifges(5,4) * t87 + Ifges(5,5) * t231;
t32 = Ifges(6,1) * t88 + Ifges(6,4) * t231 - Ifges(6,5) * t87;
t31 = Ifges(5,4) * t88 + Ifges(5,2) * t87 + Ifges(5,6) * t231;
t30 = Ifges(6,5) * t88 + Ifges(6,6) * t231 - Ifges(6,3) * t87;
t28 = Ifges(7,1) * t77 + Ifges(7,4) * t76 + Ifges(7,5) * t199;
t27 = Ifges(7,4) * t77 + Ifges(7,2) * t76 + Ifges(7,6) * t199;
t21 = mrSges(7,2) * t231 + mrSges(7,3) * t25;
t20 = -mrSges(7,1) * t231 - mrSges(7,3) * t24;
t19 = Ifges(7,1) * t44 + Ifges(7,4) * t45;
t18 = Ifges(7,4) * t44 + Ifges(7,2) * t45;
t14 = -pkin(4) * t231 + t212;
t12 = -t87 * pkin(5) + t26;
t4 = Ifges(7,1) * t24 + Ifges(7,4) * t25 - Ifges(7,5) * t231;
t3 = Ifges(7,4) * t24 + Ifges(7,2) * t25 - Ifges(7,6) * t231;
t15 = [0.2e1 * t70 * t166 + 0.2e1 * t71 * t167 + 0.2e1 * t169 * t39 + 0.2e1 * t129 * t126 + 0.2e1 * t130 * t127 + 0.2e1 * t13 * t122 + 0.2e1 * t16 * t123 + 0.2e1 * t14 * t125 + 0.2e1 * t128 * t94 + 0.2e1 * t26 * t93 + t77 * t4 + 0.2e1 * t51 * t69 + t76 * t3 + 0.2e1 * t62 * t38 + 0.2e1 * t1 * t63 + 0.2e1 * t2 * t64 + 0.2e1 * t49 * t66 + 0.2e1 * t58 * t67 + 0.2e1 * t57 * t68 + 0.2e1 * t50 * t5 - 0.2e1 * t12 * t34 + t25 * t27 + t24 * t28 + 0.2e1 * t10 * t20 + 0.2e1 * t11 * t21 + (t1 * t11 + t10 * t2 - t12 * t50) * t258 + (t13 * t49 + t14 * t51 + t26 * t62) * t259 + (t129 * t71 + t130 * t70) * t261 + (t105 * t256 - t195 * t95 + t198 * t96 + (-t198 * t137 - t195 * t138 + t199 * t208) * qJD(3) + (mrSges(3,1) * t257 - Ifges(7,5) * t77 - Ifges(7,6) * t76 + (Ifges(4,5) * t198 - t218) * t196 + t267 * t140 + t266 * t139 + (pkin(7) ^ 2 * t261 + t211 * t256 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - 0.2e1 * Ifges(7,3) + t268) * t199) * qJD(2)) * t196 + (t32 + t33) * t140 + (-t31 + t30) * t139 + (t74 + t75) * t88 + (t73 - t72) * t87 + ((mrSges(3,2) * t257 - t195 * t137 + t198 * t138 + t199 * t218) * qJD(2) + t248 + t262) * t199 + (t128 * t169 + t16 * t58 - t212 * t57) * t260 - 0.2e1 * t212 * t124; (t101 / 0.2e1 + t102 / 0.2e1) * t140 + (-t100 / 0.2e1 + t99 / 0.2e1) * t139 + (t69 - t68) * t120 + (t95 / 0.2e1 + pkin(8) * t127 + t70 * mrSges(4,3)) * t198 + (t96 / 0.2e1 - pkin(8) * t126 - t71 * mrSges(4,3)) * t195 + (t122 + t123) * t90 + (-t124 + t125) * t89 + (t1 * t106 - t10 * t44 - t107 * t2 + t11 * t45) * mrSges(7,3) + t186 * t39 + t169 * t98 + ((-t71 * t195 + t70 * t198 + (-t129 * t198 - t130 * t195) * qJD(3)) * pkin(8) - pkin(2) * t224) * m(4) + t128 * t112 + t103 * t38 - pkin(2) * t105 + t106 * t3 / 0.2e1 + t107 * t4 / 0.2e1 + t26 * t111 + t61 * t93 + t62 * t97 + t77 * t19 / 0.2e1 + t76 * t18 / 0.2e1 + t6 * t63 + t7 * t64 + t65 * t5 + t50 * t17 - t12 * t52 + t25 * t53 / 0.2e1 + t24 * t54 / 0.2e1 + t44 * t28 / 0.2e1 + t45 * t27 / 0.2e1 - t48 * t34 + t36 * t20 + t37 * t21 + (t32 / 0.2e1 + t33 / 0.2e1) * t157 + (t74 / 0.2e1 + t75 / 0.2e1) * t151 + (t72 / 0.2e1 - t73 / 0.2e1) * t150 + (t66 + t67) * t121 + (-t113 / 0.2e1 + t114 / 0.2e1) * t87 + (t115 / 0.2e1 + t116 / 0.2e1) * t88 + ((-pkin(8) * t167 - t129 * mrSges(4,3) + t138 / 0.2e1) * t198 + (-pkin(8) * t166 - t130 * mrSges(4,3) + pkin(3) * t94 + t243 / 0.2e1 - t137 / 0.2e1) * t195) * qJD(3) + (-t190 / 0.2e1 - t146 / 0.2e1 + t145 / 0.2e1 - t147 / 0.2e1 - t144 / 0.2e1 + t42 / 0.2e1 + t41 / 0.2e1 + (t235 / 0.2e1 + t175 * t253 + Ifges(3,5) + (-mrSges(4,1) * t198 + mrSges(4,2) * t195 - mrSges(3,1)) * pkin(7)) * qJD(2)) * t199 + m(6) * (t103 * t26 + t120 * t14 + t121 * t13 + t49 * t90 + t51 * t89 + t61 * t62) + m(7) * (t1 * t37 + t10 * t7 + t11 * t6 - t12 * t65 + t2 * t36 - t48 * t50) + m(5) * (t120 * t212 + t121 * t16 + t128 * t186 + t169 * t225 - t57 * t89 + t58 * t90) + (t13 * t205 + t14 * t157 - t150 * t49 + t151 * t51) * mrSges(6,2) - (t30 / 0.2e1 - t31 / 0.2e1) * t205 + (t198 * t165 / 0.2e1 + t164 * t253 - qJD(2) * (Ifges(7,5) * t107 + Ifges(7,6) * t106) / 0.2e1 - Ifges(3,6) * qJD(2) + (-t198 * t175 / 0.2e1 + t176 * t253) * qJD(3) + (mrSges(3,2) * qJD(2) + t163) * pkin(7) + (t267 * t157 - t205 * t266 + t208) * qJD(2) / 0.2e1) * t196 + (-t150 * t58 - t151 * t57 + t157 * t212 + t16 * t205) * mrSges(5,3); -0.2e1 * pkin(2) * t163 + 0.2e1 * t103 * t97 + t106 * t18 + t107 * t19 + 0.2e1 * t61 * t111 + t198 * t164 + t195 * t165 + 0.2e1 * t65 * t17 + 0.2e1 * t186 * t98 + t44 * t54 + t45 * t53 - 0.2e1 * t48 * t52 + (t101 + t102) * t157 - (t99 - t100) * t205 + (t115 + t116) * t151 + (t113 - t114) * t150 + (t235 + (0.2e1 * pkin(3) * t112 - t175) * t195) * qJD(3) + (t186 * t225 + t217) * t260 + (t103 * t61 + t217) * t259 + (t36 * t7 + t37 * t6 - t48 * t65) * t258 + 0.2e1 * (t106 * t6 - t107 * t7 - t36 * t44 + t37 * t45) * mrSges(7,3) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * (t120 * t151 - t121 * t150 + t157 * t89 + t205 * t90); (t16 * t193 - t212 * t239) * t255 + m(6) * (qJD(5) * t49 + t13 * t182 + t14 * t184) + m(7) * (t1 * t132 + t10 * t118 + t11 * t117 + t131 * t2) + t182 * t66 + t184 * t69 + t131 * t20 + t132 * t21 + t117 * t63 + t118 * t64 + qJD(5) * t122 - t70 * mrSges(4,2) + t71 * mrSges(4,1) - t14 * mrSges(6,1) - t212 * mrSges(5,1) - t16 * mrSges(5,2) + t13 * mrSges(6,3) - t203 * Ifges(4,6) + t67 * t252 + t68 * t223 - t262 + t201 - Ifges(4,5) * t220; -t206 + t190 + m(6) * qJD(5) * t121 + m(7) * (t117 * t37 + t118 * t36 + t131 * t7 + t132 * t6) - t145 + t146 + t147 + t144 + t184 * t240 - Ifges(4,6) * t229 + (t193 * t255 - mrSges(5,2) + t263) * t90 + (m(6) * t184 - t239 * t255 - mrSges(5,1) - mrSges(6,1)) * t89 + (-mrSges(4,1) * t227 + mrSges(4,2) * t229) * pkin(8) + (-t150 * t252 - t151 * t223) * mrSges(5,3) + (qJD(5) * t205 - t150 * t182) * mrSges(6,2) + (t106 * t117 - t107 * t118 - t131 * t44 + t132 * t45) * mrSges(7,3); (t117 * t132 + t118 * t131) * t258 + 0.2e1 * t242 - 0.2e1 * t241 + 0.2e1 * t263 * qJD(5); m(5) * t128 + m(6) * t26 + m(7) * t12 + t38 + t39 - t5; m(5) * t225 + m(6) * t61 + m(7) * t48 - t17 + t97 + t98; 0; 0; t194 * t21 + t197 * t20 + (-t194 * t64 + t197 * t63) * qJD(6) + m(7) * (t1 * t194 + t197 * t2 + (-t10 * t194 + t11 * t197) * qJD(6)) + m(6) * t14 + t69; t240 + m(7) * (t194 * t6 + t197 * t7 + (-t194 * t36 + t197 * t37) * qJD(6)) + m(6) * t89 + (t194 * t45 - t197 * t44 + (t106 * t197 + t107 * t194) * qJD(6)) * mrSges(7,3); m(7) * (t117 * t194 + t118 * t197) + (m(7) * (-t131 * t194 + t132 * t197) - t264) * qJD(6); 0; 0; -t201; t206; t241 - t242; 0; t264 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
