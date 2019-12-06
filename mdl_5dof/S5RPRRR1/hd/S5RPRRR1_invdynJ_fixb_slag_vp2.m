% Calculate vector of inverse dynamics joint torques for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:38
% EndTime: 2019-12-05 18:08:58
% DurationCPUTime: 7.52s
% Computational Cost: add. (2182->494), mult. (5376->702), div. (0->0), fcn. (3642->8), ass. (0->216)
t94 = cos(qJ(4));
t166 = qJD(3) * t94;
t91 = sin(qJ(3));
t168 = qJD(1) * t91;
t90 = sin(qJ(4));
t71 = -t168 * t90 + t166;
t67 = qJD(5) - t71;
t202 = t67 * Ifges(6,3);
t167 = qJD(3) * t90;
t72 = t168 * t94 + t167;
t95 = cos(qJ(3));
t160 = t95 * qJD(1);
t80 = qJD(4) - t160;
t89 = sin(qJ(5));
t93 = cos(qJ(5));
t43 = t72 * t93 + t80 * t89;
t205 = t43 * Ifges(6,5);
t42 = -t72 * t89 + t80 * t93;
t206 = t42 * Ifges(6,6);
t15 = t202 + t205 + t206;
t159 = qJ(2) * qJD(1);
t149 = t91 * t159;
t148 = t95 * t159;
t69 = qJD(2) * t90 + t148 * t94;
t45 = t149 * t89 + t69 * t93;
t203 = t45 * mrSges(6,2);
t44 = t149 * t93 - t69 * t89;
t204 = t44 * mrSges(6,1);
t244 = -t203 + t204;
t193 = t80 * Ifges(5,6);
t198 = t71 * Ifges(5,2);
t219 = Ifges(5,4) * t72;
t31 = t193 + t198 + t219;
t264 = t244 - t31 / 0.2e1 + t15 / 0.2e1;
t157 = qJD(1) * qJD(3);
t255 = t95 * t157;
t74 = qJDD(1) * t91 + t255;
t36 = qJD(4) * t71 + qJDD(3) * t90 + t74 * t94;
t263 = t36 / 0.2e1;
t37 = -qJD(4) * t72 + qJDD(3) * t94 - t74 * t90;
t262 = t37 / 0.2e1;
t73 = qJDD(1) * t95 - t91 * t157;
t70 = qJDD(4) - t73;
t261 = t70 / 0.2e1;
t162 = qJD(4) * t94;
t142 = qJ(2) * t157;
t158 = qJD(1) * qJD(2);
t78 = qJDD(1) * qJ(2) + t158;
t50 = -t142 * t91 + t78 * t95;
t28 = qJD(4) * t69 - t94 * qJDD(2) + t50 * t90;
t211 = t28 * t90;
t68 = -t94 * qJD(2) + t148 * t90;
t107 = t68 * t162 + t211;
t164 = qJD(4) * t90;
t165 = qJD(4) * t68;
t27 = qJDD(2) * t90 + t50 * t94 - t165;
t260 = -t69 * t164 + t27 * t94 + t107;
t13 = qJD(5) * t42 + t36 * t93 + t70 * t89;
t240 = t13 / 0.2e1;
t14 = -qJD(5) * t43 - t36 * t89 + t70 * t93;
t239 = t14 / 0.2e1;
t35 = qJDD(5) - t37;
t235 = t35 / 0.2e1;
t51 = t142 * t95 + t78 * t91;
t7 = qJD(5) * t44 + t27 * t93 + t51 * t89;
t259 = t7 * mrSges(6,2);
t8 = -qJD(5) * t45 - t27 * t89 + t51 * t93;
t258 = t8 * mrSges(6,1);
t97 = qJD(1) ^ 2;
t257 = t95 * t97;
t221 = mrSges(5,3) * t72;
t173 = -mrSges(5,1) * t80 - mrSges(6,1) * t42 + mrSges(6,2) * t43 + t221;
t172 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t71 + mrSges(5,2) * t72 + mrSges(4,3) * t168;
t256 = t91 * (mrSges(5,1) * t90 + mrSges(5,2) * t94);
t161 = qJD(5) * t90;
t180 = t94 * t95;
t186 = t91 * t93;
t113 = -t180 * t89 + t186;
t52 = t113 * qJD(1);
t254 = t161 * t93 + t162 * t89 + t52;
t189 = t89 * t91;
t114 = t180 * t93 + t189;
t53 = t114 * qJD(1);
t253 = -t161 * t89 + t162 * t93 - t53;
t39 = Ifges(6,4) * t42;
t17 = Ifges(6,1) * t43 + Ifges(6,5) * t67 + t39;
t183 = t93 * t17;
t194 = t80 * Ifges(5,5);
t196 = t72 * Ifges(5,1);
t66 = Ifges(5,4) * t71;
t32 = t194 + t66 + t196;
t252 = t183 + t32;
t46 = -mrSges(5,2) * t80 + t71 * mrSges(5,3);
t76 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t160;
t143 = -t94 * t46 - t76;
t251 = qJD(1) * t69;
t92 = sin(qJ(1));
t96 = cos(qJ(1));
t249 = -g(1) * t96 - g(2) * t92;
t139 = -qJD(5) + t166;
t247 = t139 * t91 + t95 * t164;
t246 = -t173 * t90 + t143;
t245 = t27 * mrSges(5,2) - Ifges(5,5) * t36 - Ifges(5,6) * t37 - Ifges(5,3) * t70;
t241 = Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * t261;
t2 = t13 * Ifges(6,4) + t14 * Ifges(6,2) + t35 * Ifges(6,6);
t243 = t2 / 0.2e1;
t242 = Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t235;
t236 = t31 / 0.2e1;
t234 = -t42 / 0.2e1;
t233 = t42 / 0.2e1;
t232 = -t43 / 0.2e1;
t231 = t43 / 0.2e1;
t230 = -t67 / 0.2e1;
t229 = t67 / 0.2e1;
t226 = g(1) * t92;
t223 = g(2) * t96;
t222 = -mrSges(5,1) * t70 - mrSges(6,1) * t14 + mrSges(6,2) * t13 + mrSges(5,3) * t36;
t220 = Ifges(4,4) * t91;
t218 = Ifges(5,4) * t90;
t217 = Ifges(5,4) * t94;
t216 = Ifges(6,4) * t43;
t215 = Ifges(6,4) * t89;
t214 = Ifges(6,4) * t93;
t12 = Ifges(6,5) * t13;
t11 = Ifges(6,6) * t14;
t29 = Ifges(6,3) * t35;
t209 = t36 * Ifges(5,4);
t207 = t37 * Ifges(5,2);
t201 = t68 * t90;
t199 = t70 * Ifges(5,6);
t197 = t71 * Ifges(5,6);
t195 = t72 * Ifges(5,5);
t192 = t80 * Ifges(5,3);
t16 = Ifges(6,2) * t42 + Ifges(6,6) * t67 + t216;
t191 = t89 * t16;
t190 = t89 * t90;
t188 = t90 * t93;
t187 = t90 * t95;
t185 = t91 * t94;
t184 = t91 * t96;
t182 = t93 * t95;
t179 = t94 * t96;
t178 = t96 * t90;
t175 = mrSges(2,1) + mrSges(3,1);
t174 = mrSges(5,2) - mrSges(6,3);
t171 = qJ(2) * t97;
t170 = Ifges(4,5) * qJD(3);
t169 = Ifges(4,6) * qJD(3);
t163 = qJD(4) * t91;
t156 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t1 = t12 + t11 + t29;
t155 = t90 * t160;
t154 = t68 * t167;
t150 = -t191 / 0.2e1;
t147 = -t168 / 0.2e1;
t141 = t78 + t158;
t140 = -qJD(5) * t94 + qJD(3);
t138 = t90 * t149;
t137 = t258 - t259;
t135 = -t223 + t226;
t134 = t7 * t93 - t8 * t89;
t132 = mrSges(5,1) * t94 - mrSges(5,2) * t90;
t112 = t185 * t89 + t182;
t59 = t185 * t93 - t89 * t95;
t130 = -mrSges(6,1) * t112 - mrSges(6,2) * t59;
t129 = mrSges(6,1) * t89 + mrSges(6,2) * t93;
t128 = Ifges(5,1) * t94 - t218;
t127 = Ifges(6,1) * t93 - t215;
t126 = Ifges(6,1) * t89 + t214;
t125 = -Ifges(5,2) * t90 + t217;
t124 = -Ifges(6,2) * t89 + t214;
t123 = Ifges(6,2) * t93 + t215;
t122 = Ifges(5,5) * t94 - Ifges(5,6) * t90;
t121 = Ifges(6,5) * t93 - Ifges(6,6) * t89;
t120 = Ifges(6,5) * t89 + Ifges(6,6) * t93;
t119 = t44 * t93 + t45 * t89;
t118 = -t44 * t89 + t45 * t93;
t117 = t69 * t94 + t201;
t115 = t140 * t95;
t111 = mrSges(6,1) * t93 - mrSges(6,2) * t89 + mrSges(5,1);
t110 = -mrSges(5,3) - t129;
t25 = -mrSges(6,2) * t67 + mrSges(6,3) * t42;
t26 = mrSges(6,1) * t67 - mrSges(6,3) * t43;
t109 = t25 * t93 - t26 * t89 + t46;
t108 = qJDD(2) - t135;
t103 = qJD(1) * t256;
t101 = t68 * mrSges(5,3) + t66 / 0.2e1 + t194 / 0.2e1 + t196 / 0.2e1 + t32 / 0.2e1;
t99 = t202 / 0.2e1 + t206 / 0.2e1 + t205 / 0.2e1 - t193 / 0.2e1 - t219 / 0.2e1 - t198 / 0.2e1 - t69 * mrSges(5,3) + t264;
t98 = qJ(2) ^ 2;
t88 = t95 ^ 2;
t87 = t91 ^ 2;
t84 = Ifges(4,4) * t160;
t81 = t87 * t171;
t77 = -mrSges(4,1) * t95 + mrSges(4,2) * t91;
t63 = t179 * t95 + t90 * t92;
t62 = t178 * t95 - t92 * t94;
t61 = t180 * t92 - t178;
t60 = t187 * t92 + t179;
t57 = Ifges(4,1) * t168 + t170 + t84;
t56 = t169 + (Ifges(4,2) * t95 + t220) * qJD(1);
t55 = t114 * qJ(2);
t54 = t113 * qJ(2);
t49 = t59 * t159;
t48 = t112 * t159;
t41 = t184 * t89 + t63 * t93;
t40 = t184 * t93 - t63 * t89;
t30 = t192 + t195 + t197;
t24 = t139 * t182 + (t140 * t89 - t164 * t93) * t91;
t23 = t140 * t186 + (-t139 * t95 + t163 * t90) * t89;
t22 = -mrSges(5,2) * t70 + mrSges(5,3) * t37;
t19 = t113 * qJD(2) + (t93 * t115 + t247 * t89) * qJ(2);
t18 = t114 * qJD(2) + (t89 * t115 - t247 * t93) * qJ(2);
t9 = t199 + t207 + t209;
t6 = -mrSges(6,2) * t35 + mrSges(6,3) * t14;
t5 = mrSges(6,1) * t35 - mrSges(6,3) * t13;
t3 = [(-t112 * t7 + t23 * t45 - t24 * t44 - t59 * t8) * mrSges(6,3) + t51 * t256 - t112 * t243 + (Ifges(6,5) * t59 - Ifges(6,6) * t112) * t235 + (Ifges(6,4) * t59 - Ifges(6,2) * t112) * t239 + (Ifges(6,1) * t59 - Ifges(6,4) * t112) * t240 + ((-t68 * mrSges(5,1) - t69 * mrSges(5,2) + t197 / 0.2e1 + t192 / 0.2e1 + t195 / 0.2e1 + qJD(2) * mrSges(4,1) - t56 / 0.2e1 + t30 / 0.2e1 - t169 / 0.2e1 + Ifges(4,4) * t147) * t91 + (qJD(2) * mrSges(4,2) + t57 / 0.2e1 + t170 / 0.2e1 + (Ifges(4,4) * t95 / 0.2e1 + (m(5) * t98 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t91) * qJD(1) + t101 * t94 + t99 * t90) * t95) * qJD(3) + (t111 * t61 + t156 * t96 - t174 * t60 + t175 * t92) * g(1) + (-mrSges(5,1) * t63 - mrSges(6,1) * t41 - mrSges(6,2) * t40 + t156 * t92 + t174 * t62 - t175 * t96) * g(2) + (Ifges(6,1) * t24 + Ifges(6,4) * t23) * t231 + (Ifges(6,4) * t24 + Ifges(6,2) * t23) * t233 + t59 * t242 + (Ifges(6,5) * t24 + Ifges(6,6) * t23) * t229 - t28 * t130 + (m(5) * t87 * t158 + qJDD(1) * mrSges(3,3) + t141 * m(3) + (-qJDD(3) * mrSges(4,1) - mrSges(5,1) * t37 + mrSges(5,2) * t36 + mrSges(4,3) * t74 + qJD(1) * t132 * t163 + t246 * qJD(3) + m(5) * (-t166 * t69 - t154 + t51) - m(6) * t154) * t91 + (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t73 + t94 * t22 + t222 * t90 + (t173 * t94 - t90 * t46) * qJD(4) + m(5) * t260 + m(6) * t107 + (t103 + t172) * qJD(3)) * t95 + (m(6) + m(5) + m(3)) * t249 + (t50 * t95 + t51 * t91 + t249 + (t87 + t88) * t158) * m(4)) * qJ(2) + t68 * (-mrSges(6,1) * t23 + mrSges(6,2) * t24) + t54 * t5 + t55 * t6 + t24 * t17 / 0.2e1 + t18 * t25 + t19 * t26 + t23 * t16 / 0.2e1 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1) + (t28 * mrSges(5,1) + Ifges(4,2) * t73 + Ifges(4,4) * t74 + Ifges(4,6) * qJDD(3) + t50 * mrSges(4,3) + t135 * mrSges(4,1) + (m(5) * t117 + m(6) * t201 - t246) * qJD(2) + t245) * t95 + m(6) * (t18 * t45 + t19 * t44 + t54 * t8 + t55 * t7) + (t51 * mrSges(4,3) + Ifges(4,1) * t74 + Ifges(4,4) * t73 + Ifges(4,5) * qJDD(3) + t172 * qJD(2) + (mrSges(4,2) - mrSges(5,3)) * t223 + (-mrSges(4,2) - t110) * t226 + (t28 * mrSges(5,3) + 0.2e1 * t241) * t94 + (-t9 / 0.2e1 + t1 / 0.2e1 - t199 / 0.2e1 - t209 / 0.2e1 - t207 / 0.2e1 - t27 * mrSges(5,3) + t29 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1 + t137) * t90 + (-t101 * t90 + t94 * t99) * qJD(4)) * t91 + t141 * mrSges(3,3) + (t77 - mrSges(3,1)) * qJDD(2); -qJDD(1) * mrSges(3,1) - t73 * mrSges(4,1) + t74 * mrSges(4,2) - t97 * mrSges(3,3) - t53 * t25 - t52 * t26 + (-t172 * t91 - t95 * t76) * qJD(1) + (t108 - t171) * m(3) + (qJD(4) * t109 - t160 * t46 - t222) * t94 + (-t89 * t5 + t93 * t6 + t22 + (-t25 * t89 - t26 * t93) * qJD(5) + t80 * t173) * t90 + (-t171 * t88 + t108 - t81) * m(4) + (-t135 - t155 * t68 - t44 * t52 - t45 * t53 + (qJD(4) * t118 - t28) * t94 + (-qJD(5) * t119 + t134 + t165) * t90) * m(6) + (t117 * t80 + t27 * t90 - t28 * t94 - t135 - t81) * m(5); (qJD(4) * t103 - t256 * t257) * qJ(2) + (t150 + t252 / 0.2e1) * t162 - t2 * t190 / 0.2e1 + (-t132 - mrSges(4,1)) * t51 + t56 * t168 / 0.2e1 + ((mrSges(5,1) * t91 - mrSges(5,3) * t180) * qJD(1) + m(6) * t138 + t253 * mrSges(6,2) + t254 * mrSges(6,1)) * t68 + (-t188 * t8 - t190 * t7 - t253 * t44 - t254 * t45) * mrSges(6,3) + (t74 - t255 / 0.2e1) * Ifges(4,5) + t94 * t9 / 0.2e1 - t94 * t1 / 0.2e1 + (-t126 * t161 + (Ifges(6,5) * t90 + t127 * t94) * qJD(4)) * t231 + (-t123 * t161 + (Ifges(6,6) * t90 + t124 * t94) * qJD(4)) * t233 + (-Ifges(6,3) * t94 + t121 * t90) * t235 + (-Ifges(6,6) * t94 + t124 * t90) * t239 + (-Ifges(6,5) * t94 + t127 * t90) * t240 + t90 * t241 + t188 * t242 + (-t120 * t161 + (Ifges(6,3) * t90 + t121 * t94) * qJD(4)) * t229 + t129 * t211 + (t187 * t251 + t260) * mrSges(5,3) + t264 * t164 + Ifges(4,6) * t73 + (-mrSges(6,1) * t114 - mrSges(6,2) * t113 - mrSges(5,3) * t91 - mrSges(6,3) * t187 - t132 * t95 + t77) * g(3) - t50 * mrSges(4,2) - t52 * t16 / 0.2e1 - t53 * t17 / 0.2e1 - t48 * t26 + t49 * t25 + (Ifges(6,5) * t53 + Ifges(6,6) * t52) * t230 - t94 * t258 + (-t158 - t249) * (mrSges(4,1) * t91 + mrSges(4,2) * t95) + t249 * (-t91 * t90 * mrSges(6,3) - mrSges(6,1) * t59 + mrSges(6,2) * t112) + t249 * (mrSges(5,3) * t95 - t132 * t91) - t172 * t148 + t173 * t138 + (-t97 * (Ifges(4,1) * t95 - t220) / 0.2e1 + mrSges(5,2) * t251 - m(5) * (-t117 * t159 + t257 * t98) + Ifges(4,6) * t157 / 0.2e1) * t91 - (-Ifges(4,2) * t168 + t90 * t15 + t94 * t32 + t57 + t84) * t160 / 0.2e1 - (t93 * t16 + t89 * t17) * t161 / 0.2e1 + (t80 * t122 + t71 * t125 + t72 * t128) * qJD(4) / 0.2e1 - (t80 * (Ifges(5,3) * t91 + t122 * t95) + t72 * (Ifges(5,5) * t91 + t128 * t95) + t71 * (Ifges(5,6) * t91 + t125 * t95)) * qJD(1) / 0.2e1 - m(6) * (t44 * t48 - t45 * t49) + (Ifges(6,4) * t53 + Ifges(6,2) * t52) * t234 - t143 * t149 + (Ifges(6,5) * t232 + Ifges(6,6) * t234 + Ifges(6,3) * t230 + t236 - t244) * t155 + Ifges(4,3) * qJDD(3) + t30 * t147 + (Ifges(5,5) * t90 + Ifges(5,6) * t94) * t261 + (Ifges(5,2) * t94 + t218) * t262 + (Ifges(5,1) * t90 + t217) * t263 + (Ifges(6,1) * t53 + Ifges(6,4) * t52) * t232 + t94 * t259; g(3) * t256 + (mrSges(5,2) * t63 + t111 * t62) * g(1) - t111 * t28 + (mrSges(5,2) * t61 + t111 * t60) * g(2) + (t68 * t129 + t121 * t229 + t124 * t233 + t127 * t231 + t183 / 0.2e1 + t150 - t119 * mrSges(6,3)) * qJD(5) + (-t173 + t221) * t69 - t72 * t204 + (-g(3) * (-mrSges(6,1) * t188 + mrSges(6,2) * t190 + mrSges(6,3) * t94) - (mrSges(5,1) * t72 + mrSges(5,2) * t71) * t159) * t91 + t71 * t191 / 0.2e1 - t80 * (Ifges(5,5) * t71 - Ifges(5,6) * t72) / 0.2e1 + (Ifges(6,3) * t72 + t121 * t71) * t230 + (Ifges(6,5) * t72 + t127 * t71) * t232 + (Ifges(6,6) * t72 + t124 * t71) * t234 + t120 * t235 + t72 * t236 + t123 * t239 + t126 * t240 + t89 * t242 + t93 * t243 + t72 * t203 + (-g(1) * t63 - g(2) * t61 + t119 * t71 + t134) * mrSges(6,3) + (-m(6) * (-t118 + t69) + t110 * t71 + t109) * t68 - (Ifges(5,1) * t71 + t15 - t219) * t72 / 0.2e1 - (-Ifges(5,2) * t72 + t252 + t66) * t71 / 0.2e1 - t245; -t68 * (mrSges(6,1) * t43 + mrSges(6,2) * t42) + (Ifges(6,1) * t42 - t216) * t232 + t16 * t231 + (Ifges(6,5) * t42 - Ifges(6,6) * t43) * t230 - t44 * t25 + t45 * t26 - g(1) * (mrSges(6,1) * t40 - mrSges(6,2) * t41) - g(2) * ((t186 * t92 - t61 * t89) * mrSges(6,1) + (-t189 * t92 - t61 * t93) * mrSges(6,2)) - g(3) * t130 + (t42 * t44 + t43 * t45) * mrSges(6,3) + t137 + t1 + (-Ifges(6,2) * t43 + t17 + t39) * t234;];
tau = t3;
