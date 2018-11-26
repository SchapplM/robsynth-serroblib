% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 15:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:40:59
% EndTime: 2018-11-23 15:41:04
% DurationCPUTime: 4.83s
% Computational Cost: add. (4193->401), mult. (8860->594), div. (0->0), fcn. (5434->8), ass. (0->193)
t129 = sin(pkin(9));
t131 = cos(pkin(9));
t128 = sin(pkin(10));
t130 = cos(pkin(10));
t134 = sin(qJ(4));
t136 = cos(qJ(4));
t100 = t128 * t134 - t130 * t136;
t242 = t100 * qJD(1);
t101 = t128 * t136 + t130 * t134;
t94 = t101 * qJD(4);
t254 = -t129 * t94 + t131 * t242;
t133 = sin(qJ(6));
t135 = cos(qJ(6));
t93 = t101 * qJD(1);
t70 = qJD(4) * t133 - t135 * t93;
t226 = Ifges(7,4) * t70;
t253 = qJD(6) - t242;
t69 = qJD(4) * t135 + t133 * t93;
t25 = Ifges(7,2) * t69 + Ifges(7,6) * t253 + t226;
t258 = -t25 / 0.2e1;
t216 = -qJD(4) * mrSges(6,1) - mrSges(7,1) * t69 + mrSges(7,2) * t70 - t93 * mrSges(6,3);
t84 = t100 * t129;
t154 = t131 * t133 + t135 * t84;
t183 = t129 * qJD(1);
t257 = qJD(6) * t154 - t133 * t254 - t135 * t183;
t65 = -t131 * t135 + t133 * t84;
t256 = qJD(6) * t65 - t133 * t183 + t135 * t254;
t185 = qJD(4) * t136;
t186 = qJD(4) * t134;
t245 = t128 * t186 - t130 * t185;
t255 = t245 * t129 + t131 * t93;
t252 = t93 / 0.2e1;
t251 = -t94 / 0.2e1;
t250 = t94 / 0.2e1;
t126 = t136 * qJD(3);
t137 = -pkin(1) - pkin(2);
t111 = qJD(1) * t137 + qJD(2);
t125 = t131 * qJ(2);
t88 = qJD(1) * t125 + t129 * t111;
t77 = -qJD(1) * pkin(7) + t88;
t63 = -t134 * t77 + t126;
t187 = qJD(3) * t134;
t64 = t136 * t77 + t187;
t156 = -t134 * t64 - t136 * t63;
t182 = qJD(1) * qJD(2);
t175 = t131 * t182;
t47 = qJD(4) * t126 + t136 * t175 - t186 * t77;
t144 = t64 * qJD(4);
t48 = -t134 * t175 - t144;
t243 = -t134 * t48 + t136 * t47;
t249 = m(5) * (qJD(4) * t156 + t243);
t248 = m(3) * qJ(2) + mrSges(3,3);
t86 = qJD(4) * t242;
t40 = qJD(6) * t69 + t135 * t86;
t41 = -qJD(6) * t70 - t133 * t86;
t15 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t247 = t86 * mrSges(6,3) + t15;
t188 = qJD(2) * t131;
t168 = -qJD(5) + t188;
t152 = t168 * t134;
t246 = -t144 + (qJ(5) * t185 - t152) * qJD(1);
t190 = qJD(1) * t134;
t108 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t190;
t189 = qJD(1) * t136;
t109 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t189;
t244 = -t134 * t108 + t136 * t109;
t155 = -t134 * t63 + t136 * t64;
t33 = (qJ(5) * t186 - qJD(5) * t136) * qJD(1) + t47;
t12 = t128 * t246 + t130 * t33;
t85 = qJD(1) * t94;
t113 = t129 * t182;
t179 = pkin(4) * t186;
t99 = -qJD(1) * t179 + t113;
t34 = -pkin(5) * t85 - pkin(8) * t86 + t99;
t59 = t187 + (-qJ(5) * qJD(1) + t77) * t136;
t52 = t130 * t59;
t58 = qJ(5) * t190 + t63;
t56 = qJD(4) * pkin(4) + t58;
t19 = t128 * t56 + t52;
t17 = qJD(4) * pkin(8) + t19;
t87 = -qJ(2) * t183 + t111 * t131;
t76 = qJD(1) * pkin(3) - t87;
t67 = pkin(4) * t189 + qJD(5) + t76;
t29 = -pkin(5) * t242 + pkin(8) * t93 + t67;
t5 = -t133 * t17 + t135 * t29;
t1 = qJD(6) * t5 + t12 * t135 + t133 * t34;
t6 = t133 * t29 + t135 * t17;
t2 = -qJD(6) * t6 - t12 * t133 + t135 * t34;
t167 = t1 * t135 - t133 * t2;
t241 = -Ifges(6,5) * qJD(4) / 0.2e1 - Ifges(6,4) * t242 / 0.2e1;
t206 = t128 * t59;
t18 = t130 * t56 - t206;
t16 = -qJD(4) * pkin(5) - t18;
t160 = Ifges(7,5) * t135 - Ifges(7,6) * t133;
t211 = Ifges(7,4) * t135;
t161 = -Ifges(7,2) * t133 + t211;
t212 = Ifges(7,4) * t133;
t162 = Ifges(7,1) * t135 - t212;
t163 = mrSges(7,1) * t133 + mrSges(7,2) * t135;
t166 = t6 * t133 + t5 * t135;
t229 = t135 / 0.2e1;
t231 = t253 / 0.2e1;
t233 = t70 / 0.2e1;
t234 = t69 / 0.2e1;
t68 = Ifges(7,4) * t69;
t26 = Ifges(7,1) * t70 + Ifges(7,5) * t253 + t68;
t240 = -t166 * mrSges(7,3) + t133 * t258 + t16 * t163 + t160 * t231 + t161 * t234 + t162 * t233 + t26 * t229;
t239 = -m(5) / 0.2e1;
t238 = t40 / 0.2e1;
t237 = t41 / 0.2e1;
t236 = t252 * Ifges(6,1) + t241;
t232 = -t85 / 0.2e1;
t230 = m(6) * t67;
t227 = Ifges(6,4) * t93;
t11 = t128 * t33 - t130 * t246;
t191 = t129 * t137 + t125;
t103 = -pkin(7) + t191;
t193 = qJ(5) - t103;
t171 = t193 * t134;
t80 = t193 * t136;
t38 = -t128 * t80 - t130 * t171;
t224 = t11 * t38;
t83 = t101 * t129;
t223 = t11 * t83;
t221 = t69 * Ifges(7,6);
t220 = t70 * Ifges(7,5);
t219 = t85 * mrSges(6,3);
t217 = t253 * Ifges(7,3);
t214 = Ifges(5,4) * t134;
t213 = Ifges(5,4) * t136;
t210 = Ifges(5,5) * t136;
t209 = Ifges(5,6) * t134;
t208 = t100 * t11;
t207 = t101 * t11;
t205 = t133 * t245;
t202 = t135 * t245;
t198 = Ifges(6,6) * qJD(4);
t197 = t101 * t133;
t196 = t101 * t135;
t184 = qJD(6) * t101;
t123 = t129 * qJD(2);
t180 = Ifges(7,5) * t40 + Ifges(7,6) * t41 - Ifges(7,3) * t85;
t49 = -t85 * mrSges(6,1) + t86 * mrSges(6,2);
t170 = -t129 * qJ(2) + t131 * t137;
t169 = qJD(4) * t193;
t102 = pkin(3) - t170;
t106 = t123 - t179;
t165 = -t5 * t133 + t6 * t135;
t164 = -mrSges(5,1) * t134 - mrSges(5,2) * t136;
t159 = -t87 * t129 + t88 * t131;
t27 = -mrSges(7,1) * t85 - mrSges(7,3) * t40;
t28 = mrSges(7,2) * t85 + mrSges(7,3) * t41;
t158 = -t133 * t27 + t135 * t28;
t39 = t128 * t171 - t130 * t80;
t91 = t136 * pkin(4) + t102;
t43 = -pkin(5) * t100 + pkin(8) * t101 + t91;
t14 = t133 * t43 + t135 * t39;
t13 = -t133 * t39 + t135 * t43;
t44 = -mrSges(7,2) * t253 + mrSges(7,3) * t69;
t45 = mrSges(7,1) * t253 - mrSges(7,3) * t70;
t157 = -t133 * t44 - t135 * t45;
t151 = t135 * t184 - t205;
t150 = t133 * t184 + t202;
t149 = t134 * (-Ifges(5,1) * t136 + t214);
t148 = t136 * (Ifges(5,2) * t134 - t213);
t104 = (mrSges(5,1) * t136 - mrSges(5,2) * t134) * qJD(1);
t147 = t164 * qJD(4);
t146 = (-t134 * Ifges(5,1) - t213) * qJD(1);
t145 = (-t136 * Ifges(5,2) - t214) * qJD(1);
t142 = t136 * t169 - t152;
t138 = qJD(1) ^ 2;
t118 = -pkin(4) * t130 - pkin(5);
t117 = pkin(4) * t128 + pkin(8);
t98 = qJD(1) * t147;
t97 = Ifges(5,5) * qJD(4) + t146;
t96 = Ifges(5,6) * qJD(4) + t145;
t78 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t242;
t60 = -mrSges(6,1) * t242 - mrSges(6,2) * t93;
t57 = t134 * t169 + t136 * t168;
t54 = Ifges(6,2) * t242 + t198 - t227;
t50 = -pkin(4) * t190 - pkin(5) * t93 - pkin(8) * t242;
t42 = -pkin(5) * t94 - pkin(8) * t245 + t106;
t24 = t217 + t220 + t221;
t23 = t130 * t58 - t206;
t22 = t128 * t142 + t130 * t57;
t21 = t128 * t58 + t52;
t10 = Ifges(7,1) * t40 + Ifges(7,4) * t41 - Ifges(7,5) * t85;
t9 = Ifges(7,4) * t40 + Ifges(7,2) * t41 - Ifges(7,6) * t85;
t8 = t133 * t50 + t135 * t23;
t7 = -t133 * t23 + t135 * t50;
t4 = -qJD(6) * t14 - t133 * t22 + t135 * t42;
t3 = qJD(6) * t13 + t133 * t42 + t135 * t22;
t20 = [t91 * t49 + t22 * t78 + t4 * t45 + t3 * t44 + t14 * t28 + t13 * t27 + (m(4) * ((-t129 * t170 + t131 * t191) * qJD(1) + t159) + m(5) * (t155 * t131 + (qJD(1) * t102 + t76) * t129)) * qJD(2) + (-t149 - t148) * qJD(1) * qJD(4) + (-t101 * t86 - t245 * t252) * Ifges(6,1) + (t100 * t12 - t18 * t245 + t19 * t94 - t207) * mrSges(6,3) + (t242 * t245 / 0.2e1 + t93 * t251 + t100 * t86 - t101 * t85) * Ifges(6,4) + t67 * (-mrSges(6,1) * t94 + mrSges(6,2) * t245) + qJD(4) * (Ifges(6,5) * t245 + Ifges(6,6) * t94) / 0.2e1 - t245 * t236 + t205 * t258 + (t133 * t26 + t135 * t25) * t184 / 0.2e1 - (t146 + t97) * t185 / 0.2e1 + (t145 + t96) * t186 / 0.2e1 + (t100 * t85 + t242 * t250) * Ifges(6,2) + (-t108 * t185 - t109 * t186 + t249) * t103 + (-m(6) * t18 + m(7) * t16 + t216) * (t128 * t57 - t130 * t142) + m(7) * (t1 * t14 + t13 * t2 + t3 * t6 + t4 * t5 + t224) + t247 * t38 + 0.2e1 * t248 * t182 + (t185 * t63 + t186 * t64 - t243) * mrSges(5,3) + t244 * t188 + t5 * (-mrSges(7,1) * t94 - mrSges(7,3) * t150) + t6 * (mrSges(7,2) * t94 + mrSges(7,3) * t151) + t16 * (-mrSges(7,1) * t151 + mrSges(7,2) * t150) + t99 * (-mrSges(6,1) * t100 - mrSges(6,2) * t101) + t102 * t98 + t106 * t60 + t54 * t250 + t24 * t251 + t76 * t147 + t39 * t219 + (Ifges(7,5) * t150 + Ifges(7,6) * t151 - Ifges(7,3) * t94) * t231 + (-Ifges(7,3) * t100 - t101 * t160) * t232 + (Ifges(7,1) * t150 + Ifges(7,4) * t151 - Ifges(7,5) * t94) * t233 + (Ifges(7,4) * t150 + Ifges(7,2) * t151 - Ifges(7,6) * t94) * t234 + (-Ifges(7,6) * t100 - t101 * t161) * t237 + (-Ifges(7,5) * t100 - t101 * t162) * t238 + 0.2e1 * mrSges(4,1) * t113 + 0.2e1 * mrSges(4,2) * t175 - t100 * t180 / 0.2e1 + 0.2e1 * t104 * t123 - t10 * t196 / 0.2e1 + t2 * (-mrSges(7,1) * t100 + mrSges(7,3) * t196) + t9 * t197 / 0.2e1 + t1 * (mrSges(7,2) * t100 + mrSges(7,3) * t197) + t26 * t202 / 0.2e1 - t163 * t207 + qJD(4) ^ 2 * (t209 - t210) / 0.2e1 + m(6) * (t106 * t67 + t12 * t39 + t19 * t22 + t91 * t99 + t224); t83 * t15 + t65 * t27 - t154 * t28 + t254 * t78 + t257 * t45 + t256 * t44 - t248 * t138 + (t83 * t86 - t84 * t85) * mrSges(6,3) + (-t138 * mrSges(4,2) - qJD(1) * t244 - t49 - t98) * t131 + (-t138 * mrSges(4,1) + (-t108 * t136 - t109 * t134) * qJD(4) + (-t104 - t60) * qJD(1)) * t129 + t129 * t249 + 0.2e1 * (-m(4) * t159 / 0.2e1 + (-t230 / 0.2e1 + t76 * t239) * t129 + (t123 + t155) * t239 * t131) * qJD(1) - t216 * t255 + (-t1 * t154 - t16 * t255 + t2 * t65 + t256 * t6 + t257 * t5 + t223) * m(7) + (-t12 * t84 - t99 * t131 + t18 * t255 + t19 * t254 + t223) * m(6); t216 * t94 + t247 * t100 - (-t133 * t45 + t135 * t44 + t78) * t245 + ((t134 ^ 2 + t136 ^ 2) * qJD(1) * mrSges(5,3) + t244) * qJD(4) + (qJD(6) * t157 + t158 + t219) * t101 + m(5) * (qJD(4) * t155 + t134 * t47 + t136 * t48) + m(6) * (t101 * t12 - t18 * t94 - t19 * t245 + t208) + m(7) * (t208 + t16 * t94 - t165 * t245 + (-qJD(6) * t166 + t167) * t101); (-t135 * mrSges(7,1) + t133 * mrSges(7,2) - mrSges(6,1)) * t11 + Ifges(6,6) * t85 + Ifges(6,5) * t86 - t23 * t78 - t7 * t45 - t47 * mrSges(5,2) + t48 * mrSges(5,1) - t8 * t44 - t12 * mrSges(6,2) + (t236 - t67 * mrSges(6,2) + (-Ifges(6,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t93 - t240 + t241) * t242 + (t18 * t242 - t19 * t93 + (t128 * t85 - t130 * t86) * pkin(4)) * mrSges(6,3) - m(7) * (t16 * t21 + t5 * t7 + t6 * t8) + m(7) * (t11 * t118 + t117 * t167) + ((-m(7) * t166 + t157) * t117 + t240) * qJD(6) + (t136 * t97 / 0.2e1 - t76 * t164 + (t149 / 0.2e1 + t148 / 0.2e1) * qJD(1) + t156 * mrSges(5,3) + (t209 / 0.2e1 - t210 / 0.2e1) * qJD(4) + (-t96 / 0.2e1 + (t60 + t230) * pkin(4)) * t134) * qJD(1) + t64 * t108 - t63 * t109 + t118 * t15 + t158 * t117 + t133 * t10 / 0.2e1 + t167 * mrSges(7,3) + t9 * t229 + (Ifges(7,5) * t133 + Ifges(7,6) * t135) * t232 + (Ifges(7,2) * t135 + t212) * t237 + (Ifges(7,1) * t133 + t211) * t238 + (t18 * t21 - t19 * t23 + (-t11 * t130 + t12 * t128) * pkin(4)) * m(6) - t216 * t21 + (-t198 / 0.2e1 + t5 * mrSges(7,1) - t6 * mrSges(7,2) + t67 * mrSges(6,1) + t217 / 0.2e1 + t221 / 0.2e1 + t220 / 0.2e1 - t54 / 0.2e1 + t24 / 0.2e1 + t227 / 0.2e1) * t93; -t242 * t78 + t216 * t93 + (t253 * t44 + t27) * t135 + (-t253 * t45 + t28) * t133 + t49 + (t1 * t133 + t2 * t135 + t16 * t93 + t165 * t253) * m(7) + (-t18 * t93 - t19 * t242 + t99) * m(6); -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t16 * (mrSges(7,1) * t70 + mrSges(7,2) * t69) - t70 * (Ifges(7,1) * t69 - t226) / 0.2e1 + t25 * t233 - t253 * (Ifges(7,5) * t69 - Ifges(7,6) * t70) / 0.2e1 - t5 * t44 + t6 * t45 + (t5 * t69 + t6 * t70) * mrSges(7,3) + t180 - (-Ifges(7,2) * t70 + t26 + t68) * t69 / 0.2e1;];
tauc  = t20(:);
