% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRRP6
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:02:35
% EndTime: 2018-11-23 15:02:41
% DurationCPUTime: 5.65s
% Computational Cost: add. (2969->431), mult. (6877->562), div. (0->0), fcn. (4172->8), ass. (0->201)
t263 = Ifges(6,1) + Ifges(7,1);
t259 = Ifges(6,5) + Ifges(7,4);
t118 = sin(qJ(4));
t192 = qJD(2) * t118;
t105 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t192;
t116 = cos(pkin(6));
t121 = cos(qJ(4));
t197 = t116 * t121;
t174 = qJD(1) * t197;
t123 = -pkin(2) - pkin(8);
t122 = cos(qJ(2));
t115 = sin(pkin(6));
t194 = qJD(1) * t115;
t175 = t122 * t194;
t140 = qJD(3) - t175;
t80 = qJD(2) * t123 + t140;
t55 = t118 * t80 + t174;
t265 = m(5) * t55 + t105;
t117 = sin(qJ(5));
t120 = cos(qJ(5));
t45 = qJD(4) * pkin(9) + t55;
t119 = sin(qJ(2));
t176 = t119 * t194;
t98 = pkin(4) * t118 - pkin(9) * t121 + qJ(3);
t72 = qJD(2) * t98 + t176;
t11 = -t117 * t45 + t120 * t72;
t12 = t117 * t72 + t120 * t45;
t139 = t11 * t120 + t117 * t12;
t206 = Ifges(7,5) * t120;
t144 = Ifges(7,3) * t117 + t206;
t208 = Ifges(6,4) * t120;
t150 = -Ifges(6,2) * t117 + t208;
t155 = mrSges(7,1) * t117 - mrSges(7,3) * t120;
t157 = mrSges(6,1) * t117 + mrSges(6,2) * t120;
t113 = qJD(5) + t192;
t249 = qJD(6) - t11;
t8 = -pkin(5) * t113 + t249;
t9 = qJ(6) * t113 + t12;
t159 = t117 * t9 - t120 * t8;
t193 = qJD(1) * t118;
t54 = -t116 * t193 + t121 * t80;
t44 = -qJD(4) * pkin(4) - t54;
t184 = t120 * qJD(4);
t190 = qJD(2) * t121;
t92 = t117 * t190 - t184;
t188 = qJD(4) * t117;
t93 = t120 * t190 + t188;
t17 = pkin(5) * t92 - qJ(6) * t93 + t44;
t226 = t120 / 0.2e1;
t228 = t117 / 0.2e1;
t229 = -t117 / 0.2e1;
t230 = t113 / 0.2e1;
t234 = t93 / 0.2e1;
t236 = t92 / 0.2e1;
t237 = -t92 / 0.2e1;
t222 = Ifges(7,5) * t92;
t87 = Ifges(6,4) * t92;
t251 = t259 * t113 + t263 * t93 + t222 - t87;
t207 = Ifges(7,5) * t117;
t209 = Ifges(6,4) * t117;
t255 = t120 * t263 + t207 - t209;
t262 = Ifges(6,6) - Ifges(7,6);
t86 = Ifges(7,5) * t93;
t32 = t113 * Ifges(7,6) + t92 * Ifges(7,3) + t86;
t223 = Ifges(6,4) * t93;
t35 = -t92 * Ifges(6,2) + t113 * Ifges(6,6) + t223;
t264 = t226 * t251 + t228 * t32 + t229 * t35 + t144 * t236 + t150 * t237 + t17 * t155 + t44 * t157 + (-t117 * t262 + t120 * t259) * t230 + t255 * t234 - t159 * mrSges(7,2) - t139 * mrSges(6,3);
t201 = Ifges(5,5) * qJD(4);
t210 = Ifges(5,4) * t118;
t260 = qJD(2) / 0.2e1;
t96 = qJD(2) * qJ(3) + t176;
t261 = t96 * mrSges(5,2) + t201 / 0.2e1 + (t121 * Ifges(5,1) - t210) * t260 - t54 * mrSges(5,3) + t264;
t183 = qJD(2) * qJD(4);
t167 = t121 * t183;
t186 = qJD(5) * t117;
t64 = qJD(5) * t184 + (-t118 * t184 - t121 * t186) * qJD(2);
t166 = t118 * t183;
t65 = qJD(5) * t93 - t117 * t166;
t258 = (-Ifges(6,4) + Ifges(7,5)) * t65 + t263 * t64 + t259 * t167;
t141 = pkin(5) * t117 - qJ(6) * t120;
t257 = qJD(5) * t141 - qJD(6) * t117 - t174 - (-qJD(2) * t141 + t80) * t118;
t256 = t117 * t263 - t206 + t208;
t254 = -m(5) * t54 + m(6) * t44;
t168 = -Ifges(5,6) * qJD(4) / 0.2e1;
t196 = t118 * t123;
t250 = t117 * t98 + t120 * t196;
t247 = t117 * t259 + t120 * t262;
t185 = qJD(5) * t120;
t191 = qJD(2) * t119;
t173 = t115 * t191;
t187 = qJD(4) * t121;
t30 = t80 * t187 + (-qJD(4) * t116 + t173) * t193;
t162 = pkin(4) * t121 + pkin(9) * t118;
t89 = qJD(4) * t162 + qJD(3);
t56 = (t89 + t175) * qJD(2);
t3 = t117 * t56 + t120 * t30 + t72 * t185 - t186 * t45;
t4 = -qJD(5) * t12 - t117 * t30 + t120 * t56;
t160 = -t117 * t4 + t120 * t3;
t1 = qJ(6) * t167 + qJD(6) * t113 + t3;
t2 = -pkin(5) * t167 - t4;
t161 = t1 * t120 + t117 * t2;
t246 = -qJD(5) * t250 + t120 * t89;
t225 = mrSges(6,3) * t92;
t68 = -mrSges(6,2) * t113 - t225;
t71 = -mrSges(7,2) * t92 + mrSges(7,3) * t113;
t211 = t71 + t68;
t224 = mrSges(6,3) * t93;
t69 = mrSges(6,1) * t113 - t224;
t70 = -mrSges(7,1) * t113 + mrSges(7,2) * t93;
t212 = -t70 + t69;
t129 = -t117 * t211 - t120 * t212;
t245 = -m(6) * t139 - m(7) * t159 + t129;
t178 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t179 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t180 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t244 = -t178 * t113 - t179 * t92 - t180 * t93 - t11 * mrSges(6,1) - t9 * mrSges(7,3) - t96 * mrSges(5,1) - Ifges(6,6) * t237 - Ifges(7,6) * t236 - t168 + (Ifges(5,4) * t121 - t118 * Ifges(5,2)) * t260 + t12 * mrSges(6,2) + t55 * mrSges(5,3) + t8 * mrSges(7,1) - t259 * t234 - (Ifges(6,3) + Ifges(7,2)) * t230;
t241 = -m(5) / 0.2e1;
t240 = t64 / 0.2e1;
t239 = -t65 / 0.2e1;
t238 = t65 / 0.2e1;
t235 = -t93 / 0.2e1;
t231 = -t113 / 0.2e1;
t227 = -t120 / 0.2e1;
t164 = t121 * t173;
t199 = qJD(4) * t55;
t31 = -qJD(1) * t164 + t199;
t198 = t115 * t122;
t78 = t116 * t118 + t121 * t198;
t217 = t31 * t78;
t19 = mrSges(7,1) * t65 - mrSges(7,3) * t64;
t20 = mrSges(6,1) * t65 + mrSges(6,2) * t64;
t215 = -t19 - t20;
t38 = -mrSges(7,2) * t65 + mrSges(7,3) * t167;
t41 = -mrSges(6,2) * t167 - mrSges(6,3) * t65;
t214 = t38 + t41;
t39 = mrSges(6,1) * t167 - mrSges(6,3) * t64;
t40 = -mrSges(7,1) * t167 + t64 * mrSges(7,2);
t213 = t40 - t39;
t95 = t162 * qJD(2);
t27 = t117 * t95 + t120 * t54;
t204 = t120 * t98;
t203 = t122 * t96;
t202 = -qJD(4) * mrSges(5,1) + mrSges(6,1) * t92 + mrSges(6,2) * t93 + mrSges(5,3) * t190;
t195 = t119 * t120;
t189 = qJD(2) * t122;
t171 = t123 * t187;
t182 = t117 * t89 + t120 * t171 + t98 * t185;
t47 = mrSges(7,1) * t92 - mrSges(7,3) * t93;
t181 = t47 + t202;
t177 = t118 * t198;
t172 = t115 * t189;
t170 = t123 * t186;
t169 = -t201 / 0.2e1;
t165 = t117 * t123 - pkin(5);
t158 = mrSges(6,1) * t120 - mrSges(6,2) * t117;
t156 = mrSges(7,1) * t120 + mrSges(7,3) * t117;
t149 = Ifges(6,2) * t120 + t209;
t143 = -Ifges(7,3) * t120 + t207;
t142 = pkin(5) * t120 + qJ(6) * t117;
t26 = -t117 * t54 + t120 * t95;
t136 = -t123 + t141;
t90 = (qJD(3) + t175) * qJD(2);
t134 = t119 * t90 + t189 * t96;
t79 = -t177 + t197;
t133 = t115 * t195 - t117 * t79;
t52 = t115 * t117 * t119 + t120 * t79;
t130 = t117 * t213 + t120 * t214;
t128 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t124 = qJD(2) ^ 2;
t112 = Ifges(7,2) * t167;
t111 = Ifges(6,3) * t167;
t99 = -pkin(4) - t142;
t94 = (mrSges(5,1) * t118 + mrSges(5,2) * t121) * qJD(2);
t91 = -qJD(2) * pkin(2) + t140;
t84 = (mrSges(5,1) * t121 - mrSges(5,2) * t118) * t183;
t75 = t136 * t121;
t73 = -t117 * t196 + t204;
t67 = (t117 * t122 + t118 * t195) * t194;
t66 = t117 * t118 * t176 - t120 * t175;
t63 = Ifges(7,4) * t64;
t62 = Ifges(6,5) * t64;
t61 = Ifges(6,6) * t65;
t60 = Ifges(7,6) * t65;
t58 = t118 * t165 - t204;
t57 = qJ(6) * t118 + t250;
t50 = -qJD(4) * t177 + t116 * t187 - t164;
t49 = -qJD(4) * t78 + t118 * t173;
t46 = pkin(5) * t93 + qJ(6) * t92;
t25 = (qJD(5) * t142 - qJD(6) * t120) * t121 - t136 * t118 * qJD(4);
t24 = -t117 * t171 + t246;
t23 = -t118 * t170 + t182;
t22 = -pkin(5) * t190 - t26;
t21 = qJ(6) * t190 + t27;
t18 = t165 * t187 - t246;
t14 = t64 * Ifges(6,4) - t65 * Ifges(6,2) + Ifges(6,6) * t167;
t13 = t64 * Ifges(7,5) + Ifges(7,6) * t167 + t65 * Ifges(7,3);
t10 = qJ(6) * t187 + (qJD(6) - t170) * t118 + t182;
t7 = qJD(5) * t52 + t117 * t49 - t120 * t172;
t6 = qJD(5) * t133 + t117 * t172 + t120 * t49;
t5 = pkin(5) * t65 - qJ(6) * t64 - qJD(6) * t93 + t31;
t15 = [-t79 * mrSges(5,3) * t167 + t49 * t105 - t212 * t7 + t211 * t6 + t214 * t52 - t213 * t133 + (-mrSges(5,3) * t166 - t215) * t78 + t181 * t50 + m(5) * (t30 * t79 + t49 * t55 - t50 * t54 + t217) + m(6) * (-t11 * t7 + t12 * t6 + t133 * t4 + t3 * t52 + t44 * t50 + t217) + m(7) * (t1 * t52 - t133 * t2 + t17 * t50 + t5 * t78 + t6 * t9 + t7 * t8) + (t94 * t189 + t119 * t84 + ((-mrSges(3,2) + mrSges(4,3)) * t122 + (-mrSges(3,1) + mrSges(4,2)) * t119) * t124 + m(5) * t134 + (-t176 * t189 + t91 * t191 + t134) * m(4)) * t115; qJ(3) * t84 + t10 * t71 + t18 * t70 + t75 * t19 + t23 * t68 + t24 * t69 + t25 * t47 + t57 * t38 + t73 * t39 + t58 * t40 + t250 * t41 + t140 * t94 - t211 * t67 + t212 * t66 + (qJD(2) * t140 + t90) * mrSges(4,3) - m(6) * (-t11 * t66 + t12 * t67) - m(7) * (t66 * t8 + t67 * t9) + m(7) * (t1 * t57 + t10 * t9 + t17 * t25 + t18 * t8 + t2 * t58 + t5 * t75) + m(6) * (t11 * t24 + t12 * t23 + t250 * t3 + t4 * t73) + 0.2e1 * (t203 * t241 - (pkin(2) * t191 + t119 * t91 + t203) * m(4) / 0.2e1) * t194 + (t111 / 0.2e1 + t112 / 0.2e1 + t63 / 0.2e1 + t60 / 0.2e1 - t61 / 0.2e1 + t62 / 0.2e1 + t90 * mrSges(5,1) + t179 * t65 + t180 * t64 + (0.3e1 / 0.2e1 * Ifges(5,4) * t192 + (t202 + t254) * t123 + t169 - t261) * qJD(4) + t128 + (m(5) * t123 - mrSges(5,3)) * t30 - t265 * t176) * t118 + (t144 * t238 + t150 * t239 + t90 * mrSges(5,2) + t13 * t228 + t14 * t229 + t5 * t155 + (-t117 * t3 - t120 * t4) * mrSges(6,3) + (-t1 * t117 + t120 * t2) * mrSges(7,2) + (((-0.3e1 / 0.2e1 * Ifges(5,4) + t180 * t120 + t179 * t117) * t121 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2) + t178) * t118) * qJD(2) + t168 - t244) * qJD(4) + (t143 * t237 + t149 * t236 + t44 * t158 + t17 * t156 + t35 * t227 + (t11 * t117 - t12 * t120) * mrSges(6,3) + (-t117 * t8 - t120 * t9) * mrSges(7,2) + t256 * t235 + t247 * t231 + t251 * t229) * qJD(5) + t255 * t240 + (t32 * qJD(5) + t258) * t226 + (t265 * qJD(4) - t20) * t123 + (m(7) * t17 + t181 + t254) * t176 + (mrSges(5,3) + t157 + (-m(6) + 0.2e1 * t241) * t123) * t31) * t121 + (m(5) + m(4)) * (qJ(3) * t90 + qJD(3) * t96); -t124 * mrSges(4,3) + (-m(5) * t96 - t94 + (-t96 + t176) * m(4) + t245) * qJD(2) + ((-t117 * t212 + t120 * t211 + t105) * qJD(4) + m(6) * (-t11 * t188 + t12 * t184 - t31) + m(7) * (t184 * t9 + t188 * t8 - t5) + m(5) * (-t31 + t199) + t215) * t121 + (t181 * qJD(4) + t129 * qJD(5) + m(6) * (qJD(4) * t44 - t11 * t185 - t12 * t186 + t160) + m(7) * (qJD(4) * t17 + t185 * t8 - t186 * t9 + t161) + m(5) * (-qJD(4) * t54 + t30) + t130) * t118; t256 * t240 + t99 * t19 - t54 * t105 - t22 * t70 - t21 * t71 - t27 * t68 - t26 * t69 - t30 * mrSges(5,2) - pkin(4) * t20 + (-pkin(4) * t31 + pkin(9) * t160 - t11 * t26 - t12 * t27 - t44 * t55) * m(6) + (pkin(9) * t161 + t17 * t257 - t21 * t9 - t22 * t8 + t5 * t99) * m(7) + t257 * t47 + t258 * t228 + ((Ifges(5,4) * t190 / 0.2e1 + t168 + t247 * qJD(4) / 0.2e1 + t244) * t121 + ((-t210 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t121) * qJD(2) + t169 + t261) * t118) * qJD(2) + t160 * mrSges(6,3) + t161 * mrSges(7,2) - t5 * t156 + (-mrSges(5,1) - t158) * t31 + t130 * pkin(9) + t143 * t238 + t149 * t239 + t14 * t226 + t13 * t227 + (pkin(9) * t245 + t264) * qJD(5) - t202 * t55; t111 + t112 + t63 + t60 - t61 + t62 - t44 * (mrSges(6,1) * t93 - mrSges(6,2) * t92) - t17 * (mrSges(7,1) * t93 + mrSges(7,3) * t92) + qJD(6) * t71 + qJ(6) * t38 - pkin(5) * t40 - t46 * t47 + (t8 * t92 + t9 * t93) * mrSges(7,2) + t35 * t234 + (Ifges(7,3) * t93 - t222) * t237 + t128 + (t212 + t224) * t12 + (-t211 - t225) * t11 + (-t259 * t92 - t262 * t93) * t231 + (-pkin(5) * t2 + qJ(6) * t1 - t12 * t8 - t17 * t46 + t249 * t9) * m(7) + (-Ifges(6,2) * t93 + t251 - t87) * t236 + (-t263 * t92 - t223 + t32 + t86) * t235; -t113 * t71 + t93 * t47 + 0.2e1 * (t2 / 0.2e1 + t9 * t231 + t17 * t234) * m(7) + t40;];
tauc  = t15(:);
