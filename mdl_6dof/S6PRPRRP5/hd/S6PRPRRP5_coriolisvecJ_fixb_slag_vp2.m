% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:04
% EndTime: 2019-03-08 20:14:15
% DurationCPUTime: 5.99s
% Computational Cost: add. (2958->423), mult. (6920->558), div. (0->0), fcn. (4236->8), ass. (0->200)
t275 = Ifges(6,4) + Ifges(7,4);
t276 = Ifges(6,1) + Ifges(7,1);
t265 = Ifges(6,5) + Ifges(7,5);
t274 = Ifges(6,2) + Ifges(7,2);
t264 = Ifges(7,6) + Ifges(6,6);
t121 = sin(qJ(4));
t197 = qJD(2) * t121;
t178 = mrSges(5,3) * t197;
t108 = -qJD(4) * mrSges(5,2) - t178;
t119 = cos(pkin(6));
t124 = cos(qJ(4));
t203 = t119 * t124;
t173 = qJD(1) * t203;
t126 = -pkin(2) - pkin(8);
t125 = cos(qJ(2));
t118 = sin(pkin(6));
t199 = qJD(1) * t118;
t174 = t125 * t199;
t140 = qJD(3) - t174;
t83 = qJD(2) * t126 + t140;
t57 = t121 * t83 + t173;
t273 = m(5) * t57 + t108;
t120 = sin(qJ(5));
t123 = cos(qJ(5));
t47 = qJD(4) * pkin(9) + t57;
t104 = pkin(4) * t121 - pkin(9) * t124 + qJ(3);
t122 = sin(qJ(2));
t175 = t122 * t199;
t74 = qJD(2) * t104 + t175;
t14 = -t120 * t47 + t123 * t74;
t15 = t120 * t74 + t123 * t47;
t139 = t120 * t15 + t123 * t14;
t153 = mrSges(7,1) * t120 + mrSges(7,2) * t123;
t155 = mrSges(6,1) * t120 + mrSges(6,2) * t123;
t116 = qJD(5) + t197;
t193 = qJD(4) * t120;
t195 = qJD(2) * t124;
t99 = t123 * t195 + t193;
t8 = -qJ(6) * t99 + t14;
t5 = pkin(5) * t116 + t8;
t192 = qJD(4) * t123;
t98 = -t120 * t195 + t192;
t9 = qJ(6) * t98 + t15;
t157 = t120 * t9 + t123 * t5;
t229 = t123 / 0.2e1;
t232 = -t120 / 0.2e1;
t233 = t116 / 0.2e1;
t238 = t99 / 0.2e1;
t240 = t98 / 0.2e1;
t271 = t275 * t98;
t250 = t265 * t116 + t276 * t99 + t271;
t267 = t275 * t99;
t251 = t264 * t116 + t274 * t98 + t267;
t268 = t275 * t120;
t256 = t276 * t123 - t268;
t269 = t275 * t123;
t258 = -t274 * t120 + t269;
t198 = qJD(1) * t121;
t56 = -t119 * t198 + t124 * t83;
t46 = -qJD(4) * pkin(4) - t56;
t27 = -pkin(5) * t98 + qJD(6) + t46;
t272 = t229 * t250 + t232 * t251 + t27 * t153 + t46 * t155 + (-t264 * t120 + t265 * t123) * t233 + t258 * t240 + t256 * t238 - t139 * mrSges(6,3) - t157 * mrSges(7,3);
t102 = qJD(2) * qJ(3) + t175;
t209 = Ifges(5,5) * qJD(4);
t215 = Ifges(5,4) * t121;
t266 = qJD(2) / 0.2e1;
t270 = t102 * mrSges(5,2) + t209 / 0.2e1 + (t124 * Ifges(5,1) - t215) * t266 - t56 * mrSges(5,3) + t272;
t185 = qJD(2) * qJD(4);
t165 = t124 * t185;
t168 = t121 * t192;
t184 = qJD(4) * qJD(5);
t187 = qJD(5) * t124;
t65 = t123 * t184 + (-t120 * t187 - t168) * qJD(2);
t169 = t123 * t187;
t134 = t121 * t193 - t169;
t66 = qJD(2) * t134 - t120 * t184;
t263 = t264 * t165 + t274 * t66 + t275 * t65;
t262 = t265 * t165 + t275 * t66 + t276 * t65;
t221 = -qJ(6) - pkin(9);
t163 = qJD(5) * t221;
t186 = qJD(6) * t123;
t160 = pkin(4) * t124 + pkin(9) * t121;
t101 = t160 * qJD(2);
t26 = t120 * t101 + t123 * t56;
t261 = t186 - t26 + (-qJ(6) * t197 + t163) * t120;
t25 = t123 * t101 - t120 * t56;
t260 = -qJD(6) * t120 + t123 * t163 - (qJ(6) * t121 * t123 + pkin(5) * t124) * qJD(2) - t25;
t259 = t274 * t123 + t268;
t257 = t276 * t120 + t269;
t255 = -m(5) * t56 + m(6) * t46;
t166 = -Ifges(5,6) * qJD(4) / 0.2e1;
t252 = m(6) * t139;
t248 = t265 * t120 + t264 * t123;
t188 = qJD(5) * t123;
t189 = qJD(5) * t120;
t196 = qJD(2) * t122;
t172 = t118 * t196;
t191 = qJD(4) * t124;
t29 = t83 * t191 + (-qJD(4) * t119 + t172) * t198;
t94 = qJD(4) * t160 + qJD(3);
t58 = (t94 + t174) * qJD(2);
t3 = t120 * t58 + t123 * t29 + t74 * t188 - t189 * t47;
t4 = -qJD(5) * t15 - t120 * t29 + t123 * t58;
t158 = -t120 * t4 + t123 * t3;
t179 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t180 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t181 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t247 = -t179 * t116 - t180 * t98 + t181 * t99 - t102 * mrSges(5,1) - t14 * mrSges(6,1) - t5 * mrSges(7,1) - t166 + (Ifges(5,4) * t124 - t121 * Ifges(5,2)) * t266 + t15 * mrSges(6,2) + t57 * mrSges(5,3) + t9 * mrSges(7,2) - t264 * t240 - t265 * t238 - (Ifges(7,3) + Ifges(6,3)) * t233;
t244 = -m(5) / 0.2e1;
t243 = t65 / 0.2e1;
t242 = t66 / 0.2e1;
t241 = -t98 / 0.2e1;
t239 = -t99 / 0.2e1;
t235 = m(7) * t27;
t234 = -t116 / 0.2e1;
t228 = pkin(5) * t120;
t162 = t124 * t172;
t206 = qJD(4) * t57;
t30 = -qJD(1) * t162 + t206;
t204 = t118 * t125;
t81 = t119 * t121 + t124 * t204;
t225 = t30 * t81;
t21 = -t66 * mrSges(7,1) + t65 * mrSges(7,2);
t22 = -mrSges(6,1) * t66 + mrSges(6,2) * t65;
t220 = -t21 - t22;
t39 = mrSges(7,1) * t165 - mrSges(7,3) * t65;
t40 = mrSges(6,1) * t165 - mrSges(6,3) * t65;
t219 = t39 + t40;
t41 = -mrSges(7,2) * t165 + mrSges(7,3) * t66;
t42 = -mrSges(6,2) * t165 + mrSges(6,3) * t66;
t218 = t41 + t42;
t70 = -mrSges(7,2) * t116 + mrSges(7,3) * t98;
t71 = -mrSges(6,2) * t116 + mrSges(6,3) * t98;
t217 = t70 + t71;
t72 = mrSges(7,1) * t116 - mrSges(7,3) * t99;
t73 = mrSges(6,1) * t116 - mrSges(6,3) * t99;
t216 = t72 + t73;
t201 = t121 * t126;
t76 = t120 * t104 + t123 * t201;
t210 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t98 + mrSges(6,2) * t99 + mrSges(5,3) * t195;
t207 = qJ(6) * t124;
t205 = t102 * t125;
t202 = t120 * t122;
t200 = t122 * t123;
t194 = qJD(2) * t125;
t190 = qJD(4) * t126;
t170 = t124 * t190;
t183 = t104 * t188 + t120 * t94 + t123 * t170;
t49 = -mrSges(7,1) * t98 + mrSges(7,2) * t99;
t182 = t49 + t210;
t177 = t121 * t204;
t176 = t120 * t201;
t171 = t118 * t194;
t167 = -t209 / 0.2e1;
t164 = -t120 * t126 + pkin(5);
t1 = pkin(5) * t165 - qJ(6) * t65 - qJD(6) * t99 + t4;
t2 = qJ(6) * t66 + qJD(6) * t98 + t3;
t159 = -t1 * t120 + t123 * t2;
t156 = mrSges(6,1) * t123 - mrSges(6,2) * t120;
t154 = mrSges(7,1) * t123 - mrSges(7,2) * t120;
t82 = -t177 + t203;
t53 = t118 * t200 - t120 * t82;
t54 = t118 * t202 + t123 * t82;
t95 = (qJD(3) + t174) * qJD(2);
t136 = t102 * t194 + t122 * t95;
t133 = -t120 * t217 - t123 * t216;
t132 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t131 = -qJD(5) * t76 + t123 * t94;
t127 = qJD(2) ^ 2;
t117 = -pkin(5) * t123 - pkin(4);
t115 = Ifges(6,3) * t165;
t114 = Ifges(7,3) * t165;
t111 = t221 * t123;
t110 = t221 * t120;
t100 = (mrSges(5,1) * t121 + mrSges(5,2) * t124) * qJD(2);
t97 = -qJD(2) * pkin(2) + t140;
t96 = (-t126 + t228) * t124;
t93 = t123 * t104;
t88 = (mrSges(5,1) * t124 - mrSges(5,2) * t121) * t185;
t75 = t93 - t176;
t69 = (t120 * t125 + t121 * t200) * t199;
t68 = (-t121 * t202 + t123 * t125) * t199;
t67 = -pkin(5) * t134 + t121 * t190;
t64 = Ifges(6,5) * t65;
t63 = Ifges(7,5) * t65;
t62 = Ifges(6,6) * t66;
t61 = Ifges(7,6) * t66;
t52 = -qJD(4) * t177 + t119 * t191 - t162;
t51 = -qJD(4) * t81 + t121 * t172;
t48 = -t120 * t207 + t76;
t38 = t121 * t164 - t123 * t207 + t93;
t37 = t173 + (-qJD(2) * t228 + t83) * t121;
t24 = -t120 * t170 + t131;
t23 = -qJD(5) * t176 + t183;
t12 = -pkin(5) * t66 + t30;
t11 = qJD(5) * t53 + t120 * t171 + t123 * t51;
t10 = -qJD(5) * t54 - t120 * t51 + t123 * t171;
t7 = -qJ(6) * t169 + (-qJD(6) * t124 + (qJ(6) * qJD(4) - qJD(5) * t126) * t121) * t120 + t183;
t6 = qJ(6) * t168 + (qJ(6) * t189 + qJD(4) * t164 - t186) * t124 + t131;
t13 = [-t82 * mrSges(5,3) * t165 + t51 * t108 + t218 * t54 + t219 * t53 + t217 * t11 + t216 * t10 + (-qJD(4) * t178 - t220) * t81 + t182 * t52 + m(5) * (t29 * t82 + t51 * t57 - t52 * t56 + t225) + m(6) * (t10 * t14 + t11 * t15 + t3 * t54 + t4 * t53 + t46 * t52 + t225) + m(7) * (t1 * t53 + t10 * t5 + t11 * t9 + t12 * t81 + t2 * t54 + t27 * t52) + (t100 * t194 + t122 * t88 + ((-mrSges(3,2) + mrSges(4,3)) * t125 + (-mrSges(3,1) + mrSges(4,2)) * t122) * t127 + m(5) * t136 + (-t175 * t194 + t196 * t97 + t136) * m(4)) * t118; qJ(3) * t88 + t96 * t21 + t23 * t71 + t24 * t73 + t38 * t39 + t75 * t40 + t48 * t41 + t76 * t42 + t67 * t49 + t6 * t72 + t7 * t70 - t217 * t69 - t216 * t68 + t140 * t100 + (qJD(2) * t140 + t95) * mrSges(4,3) + m(6) * (t14 * t24 + t15 * t23 + t3 * t76 + t4 * t75) - m(6) * (t14 * t68 + t15 * t69) + 0.2e1 * (t205 * t244 - (pkin(2) * t196 + t122 * t97 + t205) * m(4) / 0.2e1) * t199 + (t114 / 0.2e1 + t115 / 0.2e1 + t62 / 0.2e1 + t63 / 0.2e1 + t64 / 0.2e1 + t61 / 0.2e1 + t95 * mrSges(5,1) + t180 * t66 - t181 * t65 + (0.3e1 / 0.2e1 * Ifges(5,4) * t197 + (t210 + t255) * t126 + t167 - t270) * qJD(4) + t132 + (m(5) * t126 - mrSges(5,3)) * t29 - t273 * t175) * t121 + (t12 * t153 + t95 * mrSges(5,2) + (-t1 * t123 - t120 * t2) * mrSges(7,3) + (-t120 * t3 - t123 * t4) * mrSges(6,3) + (((-0.3e1 / 0.2e1 * Ifges(5,4) - t181 * t123 - t180 * t120) * t124 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2) + t179) * t121) * qJD(2) + t166 - t247) * qJD(4) + (t46 * t156 + t27 * t154 + (t120 * t5 - t123 * t9) * mrSges(7,3) + (t120 * t14 - t123 * t15) * mrSges(6,3) + t259 * t241 + t257 * t239 + t248 * t234 - t251 * t123 / 0.2e1) * qJD(5) + t256 * t243 + t258 * t242 + (qJD(5) * t250 + t263) * t232 + t262 * t229 + (t273 * qJD(4) - t22) * t126 + (t182 + t235 + t255) * t175 + (mrSges(5,3) + t155 + (-m(6) + 0.2e1 * t244) * t126) * t30) * t124 + (m(5) + m(4)) * (qJ(3) * t95 + qJD(3) * t102) + (t1 * t38 + t12 * t96 + t2 * t48 + t27 * t67 + (-t69 + t7) * t9 + (t6 - t68) * t5) * m(7); -t127 * mrSges(4,3) + (-m(5) * t102 - t100 + (-t102 + t175) * m(4) - t252 - m(7) * t157 + t133) * qJD(2) + ((-t120 * t216 + t123 * t217 + t108) * qJD(4) + m(6) * (-t14 * t193 + t15 * t192 - t30) + m(7) * (t192 * t9 - t193 * t5 - t12) + m(5) * (-t30 + t206) + t220) * t124 + (t218 * t123 - t219 * t120 + t182 * qJD(4) + t133 * qJD(5) + m(6) * (qJD(4) * t46 - t14 * t188 - t15 * t189 + t158) + m(7) * (qJD(4) * t27 - t188 * t5 - t189 * t9 + t159) + m(5) * (-qJD(4) * t56 + t29)) * t121; ((Ifges(5,4) * t195 / 0.2e1 + t166 + t248 * qJD(4) / 0.2e1 + t247) * t124 + ((-t215 / 0.2e1 + (-Ifges(5,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t124) * qJD(2) + t167 + t270) * t121) * qJD(2) + t117 * t21 - t56 * t108 + t110 * t39 - t111 * t41 - t26 * t71 - t25 * t73 - t37 * t49 - t29 * mrSges(5,2) - pkin(4) * t22 + (-pkin(4) * t30 + pkin(9) * t158 - t14 * t25 - t15 * t26 - t46 * t57) * m(6) + t257 * t243 + t259 * t242 + t260 * t72 + (t1 * t110 - t111 * t2 + t117 * t12 + t260 * t5 + t261 * t9 - t27 * t37) * m(7) + t261 * t70 + t262 * t120 / 0.2e1 + t263 * t229 - t210 * t57 + ((t49 + t235) * t228 + (-t120 * t71 - t123 * t73 - t252) * pkin(9) + t272) * qJD(5) - t12 * t154 + (-mrSges(5,1) - t156) * t30 + t158 * mrSges(6,3) + t159 * mrSges(7,3) + (-t120 * t40 + t123 * t42) * pkin(9); t114 + t115 + t62 + t63 + t64 + t61 - t46 * (mrSges(6,1) * t99 + mrSges(6,2) * t98) - t27 * (mrSges(7,1) * t99 + mrSges(7,2) * t98) - t8 * t70 - t14 * t71 + t9 * t72 + t15 * t73 + t132 + (-t49 * t99 + t39) * pkin(5) + (t5 * t98 + t9 * t99) * mrSges(7,3) + (t14 * t98 + t15 * t99) * mrSges(6,3) + (-(-t5 + t8) * t9 + (-t27 * t99 + t1) * pkin(5)) * m(7) + (t276 * t98 - t267) * t239 + t251 * t238 + (-t264 * t99 + t265 * t98) * t234 + (-t274 * t99 + t250 + t271) * t241; -t98 * t70 + t99 * t72 + 0.2e1 * (t12 / 0.2e1 + t5 * t238 + t9 * t241) * m(7) + t21;];
tauc  = t13(:);
