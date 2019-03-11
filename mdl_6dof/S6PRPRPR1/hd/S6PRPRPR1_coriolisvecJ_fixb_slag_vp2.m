% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:44
% EndTime: 2019-03-08 19:24:56
% DurationCPUTime: 5.49s
% Computational Cost: add. (4301->414), mult. (11088->581), div. (0->0), fcn. (8312->12), ass. (0->210)
t141 = sin(qJ(4));
t192 = qJD(4) * t141;
t136 = sin(pkin(11));
t137 = sin(pkin(6));
t138 = cos(pkin(11));
t142 = sin(qJ(2));
t145 = cos(qJ(2));
t100 = (t136 * t145 + t138 * t142) * t137;
t91 = qJD(1) * t100;
t263 = pkin(4) * t192 - t91;
t135 = sin(pkin(12));
t144 = cos(qJ(4));
t131 = pkin(2) * t136 + pkin(8);
t196 = qJ(5) + t131;
t172 = qJD(4) * t196;
t149 = -qJD(5) * t141 - t144 * t172;
t204 = cos(pkin(12));
t153 = -t135 * t141 + t144 * t204;
t191 = qJD(5) * t144;
t90 = -t141 * t172 + t191;
t195 = qJD(1) * t137;
t179 = t142 * t195;
t120 = t136 * t179;
t178 = t145 * t195;
t94 = t138 * t178 - t120;
t262 = -t135 * t149 + t153 * t94 - t204 * t90;
t175 = t204 * t141;
t116 = t135 * t144 + t175;
t108 = t116 * qJD(4);
t109 = t153 * qJD(4);
t261 = pkin(5) * t108 - pkin(9) * t109 + t263;
t106 = t153 * qJD(2);
t235 = -t106 / 0.2e1;
t260 = Ifges(6,2) * t235 - Ifges(6,6) * qJD(4) / 0.2e1;
t140 = sin(qJ(6));
t143 = cos(qJ(6));
t181 = -pkin(4) * t144 - pkin(3);
t229 = pkin(2) * t138;
t122 = t181 - t229;
t68 = -pkin(5) * t153 - pkin(9) * t116 + t122;
t114 = t196 * t144;
t174 = t196 * t141;
t72 = t114 * t204 - t135 * t174;
t30 = -t140 * t72 + t143 * t68;
t259 = qJD(6) * t30 + t140 * t261 - t143 * t262;
t31 = t140 * t68 + t143 * t72;
t258 = -qJD(6) * t31 + t140 * t262 + t143 * t261;
t257 = t116 * t94 - t135 * t90 + t149 * t204;
t105 = qJD(6) - t106;
t202 = Ifges(6,5) * qJD(4);
t119 = qJD(2) * pkin(2) + t178;
t83 = t119 * t138 - t120;
t73 = qJD(2) * t181 + qJD(5) - t83;
t255 = t73 * mrSges(6,2) + t202 / 0.2e1;
t193 = qJD(2) * t144;
t107 = -qJD(2) * t175 - t135 * t193;
t220 = mrSges(6,3) * t107;
t88 = qJD(4) * t143 + t107 * t140;
t89 = qJD(4) * t140 - t107 * t143;
t222 = -qJD(4) * mrSges(6,1) - mrSges(7,1) * t88 + mrSges(7,2) * t89 - t220;
t188 = qJD(2) * qJD(4);
t177 = t141 * t188;
t99 = (t136 * t142 - t138 * t145) * t137;
t84 = t119 * t136 + t138 * t179;
t80 = qJD(2) * pkin(8) + t84;
t173 = qJ(5) * qJD(2) + t80;
t139 = cos(pkin(6));
t127 = qJD(1) * t139 + qJD(3);
t197 = t141 * t127;
t60 = t144 * t173 + t197;
t93 = qJD(2) * t99;
t87 = qJD(1) * t93;
t252 = -t60 * qJD(4) + (-qJD(2) * qJD(5) + t87) * t141;
t101 = qJD(2) * t108;
t102 = qJD(2) * t109;
t57 = qJD(6) * t88 + t102 * t143;
t37 = mrSges(7,1) * t101 - mrSges(7,3) * t57;
t58 = -qJD(6) * t89 - t102 * t140;
t38 = -mrSges(7,2) * t101 + mrSges(7,3) * t58;
t251 = -t140 * t37 + t143 * t38;
t121 = t144 * t127;
t32 = qJD(4) * t121 - t144 * t87 - t192 * t80;
t27 = (-qJ(5) * t192 + t191) * qJD(2) + t32;
t10 = t135 * t252 + t204 * t27;
t92 = qJD(2) * t100;
t86 = qJD(1) * t92;
t76 = pkin(4) * t177 + t86;
t39 = pkin(5) * t101 - pkin(9) * t102 + t76;
t50 = t204 * t60;
t59 = -t141 * t173 + t121;
t52 = qJD(4) * pkin(4) + t59;
t20 = t135 * t52 + t50;
t18 = qJD(4) * pkin(9) + t20;
t40 = -pkin(5) * t106 + pkin(9) * t107 + t73;
t5 = -t140 * t18 + t143 * t40;
t1 = qJD(6) * t5 + t10 * t143 + t140 * t39;
t6 = t140 * t40 + t143 * t18;
t2 = -qJD(6) * t6 - t10 * t140 + t143 * t39;
t250 = t1 * t143 - t140 * t2;
t249 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t57 + Ifges(7,6) * t58;
t217 = Ifges(6,4) * t107;
t243 = t217 / 0.2e1 + t260;
t166 = mrSges(7,1) * t140 + mrSges(7,2) * t143;
t209 = t135 * t60;
t19 = t204 * t52 - t209;
t17 = -qJD(4) * pkin(5) - t19;
t155 = t17 * t166;
t161 = Ifges(7,5) * t143 - Ifges(7,6) * t140;
t215 = Ifges(7,4) * t143;
t163 = -Ifges(7,2) * t140 + t215;
t216 = Ifges(7,4) * t140;
t165 = Ifges(7,1) * t143 - t216;
t232 = t143 / 0.2e1;
t233 = -t140 / 0.2e1;
t240 = t89 / 0.2e1;
t223 = t89 * Ifges(7,4);
t35 = Ifges(7,2) * t88 + Ifges(7,6) * t105 + t223;
t82 = Ifges(7,4) * t88;
t36 = Ifges(7,1) * t89 + Ifges(7,5) * t105 + t82;
t248 = t155 + t105 * t161 / 0.2e1 + t165 * t240 + t88 * t163 / 0.2e1 + t36 * t232 + t35 * t233;
t247 = t73 * mrSges(6,1) + t5 * mrSges(7,1) - t6 * mrSges(7,2) + t243;
t246 = qJD(2) ^ 2;
t245 = t57 / 0.2e1;
t244 = t58 / 0.2e1;
t242 = -t88 / 0.2e1;
t241 = -t89 / 0.2e1;
t77 = -t100 * t141 + t139 * t144;
t78 = t100 * t144 + t139 * t141;
t41 = t135 * t78 - t204 * t77;
t9 = t135 * t27 - t204 * t252;
t239 = t41 * t9;
t71 = t114 * t135 + t174 * t204;
t238 = t71 * t9;
t237 = t101 / 0.2e1;
t236 = -t105 / 0.2e1;
t234 = t107 / 0.2e1;
t231 = Ifges(7,5) * t89;
t230 = Ifges(7,6) * t88;
t228 = pkin(4) * t135;
t226 = t153 * t9;
t224 = t86 * t99;
t221 = mrSges(6,3) * t106;
t219 = Ifges(6,1) * t107;
t218 = Ifges(5,4) * t141;
t104 = Ifges(6,4) * t106;
t212 = Ifges(7,3) * t105;
t211 = t101 * mrSges(6,3);
t210 = t102 * mrSges(6,3);
t63 = -t141 * t80 + t121;
t207 = t141 * t63;
t74 = -mrSges(6,1) * t106 - mrSges(6,2) * t107;
t205 = t74 + (-mrSges(5,1) * t144 + mrSges(5,2) * t141) * qJD(2);
t203 = Ifges(5,5) * qJD(4);
t201 = Ifges(5,6) * qJD(4);
t199 = t106 * t140;
t198 = t106 * t143;
t194 = qJD(2) * t141;
t190 = qJD(6) * t140;
t189 = qJD(6) * t143;
t187 = pkin(4) * t194;
t185 = mrSges(5,3) * t194;
t184 = mrSges(5,3) * t193;
t180 = t204 * pkin(4);
t65 = t101 * mrSges(6,1) + mrSges(6,2) * t102;
t170 = -t1 * t140 - t143 * t2;
t169 = -t140 * t6 - t143 * t5;
t168 = t140 * t5 - t143 * t6;
t167 = mrSges(7,1) * t143 - mrSges(7,2) * t140;
t164 = Ifges(7,1) * t140 + t215;
t162 = Ifges(7,2) * t143 + t216;
t160 = Ifges(7,5) * t140 + Ifges(7,6) * t143;
t42 = t135 * t77 + t204 * t78;
t25 = t140 * t99 + t143 * t42;
t24 = -t140 * t42 + t143 * t99;
t61 = -mrSges(7,2) * t105 + mrSges(7,3) * t88;
t62 = mrSges(7,1) * t105 - mrSges(7,3) * t89;
t159 = -t140 * t62 + t143 * t61;
t64 = t144 * t80 + t197;
t95 = -qJD(4) * mrSges(6,2) + t221;
t156 = -t159 - t95;
t154 = (mrSges(5,1) * t141 + mrSges(5,2) * t144) * qJD(2);
t147 = qJD(6) * t169 + t250;
t134 = Ifges(5,4) * t193;
t133 = -pkin(3) - t229;
t132 = -t180 - pkin(5);
t124 = -qJD(4) * mrSges(5,2) + t184;
t123 = qJD(4) * mrSges(5,1) - t185;
t112 = qJD(4) * t154;
t111 = Ifges(5,1) * t194 + t134 + t203;
t110 = t201 + (t144 * Ifges(5,2) + t218) * qJD(2);
t98 = Ifges(7,3) * t101;
t79 = -qJD(2) * pkin(3) - t83;
t70 = t104 + t202 - t219;
t66 = -pkin(5) * t107 - pkin(9) * t106 + t187;
t45 = qJD(4) * t77 - t144 * t93;
t44 = -qJD(4) * t78 + t141 * t93;
t34 = t212 + t230 + t231;
t33 = -qJD(4) * t64 + t141 * t87;
t23 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t22 = t204 * t59 - t209;
t21 = t135 * t59 + t50;
t16 = t135 * t44 + t204 * t45;
t15 = t135 * t45 - t204 * t44;
t14 = Ifges(7,1) * t57 + Ifges(7,4) * t58 + Ifges(7,5) * t101;
t13 = Ifges(7,4) * t57 + Ifges(7,2) * t58 + Ifges(7,6) * t101;
t12 = t140 * t66 + t143 * t22;
t11 = -t140 * t22 + t143 * t66;
t4 = -qJD(6) * t25 - t140 * t16 + t143 * t92;
t3 = qJD(6) * t24 + t140 * t92 + t143 * t16;
t7 = [t44 * t123 + t45 * t124 + t16 * t95 + t41 * t23 + t24 * t37 + t25 * t38 + t3 * t61 + t4 * t62 + (t65 + t112) * t99 + t205 * t92 + t222 * t15 + (-mrSges(3,1) * t142 - mrSges(3,2) * t145) * t246 * t137 + (-t101 * t42 + t102 * t41) * mrSges(6,3) + (-t92 * mrSges(4,1) + t93 * mrSges(4,2) + (-t141 * t78 - t144 * t77) * qJD(4) * mrSges(5,3)) * qJD(2) + m(5) * (t32 * t78 + t33 * t77 + t44 * t63 + t45 * t64 + t79 * t92 + t224) + m(6) * (t10 * t42 - t15 * t19 + t16 * t20 + t73 * t92 + t76 * t99 + t239) + m(4) * (-t100 * t87 - t83 * t92 - t84 * t93 + t224) + m(7) * (t1 * t25 + t15 * t17 + t2 * t24 + t3 * t6 + t4 * t5 + t239); (t70 / 0.2e1 + t104 / 0.2e1 - t219 / 0.2e1 + t169 * mrSges(7,3) + t248 + t255) * t109 - (-t10 * mrSges(6,3) + t98 / 0.2e1 - Ifges(6,4) * t102 + t76 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t101 + t249) * t153 + (t230 / 0.2e1 + t231 / 0.2e1 + t212 / 0.2e1 + t34 / 0.2e1 + t247 + t243) * t108 + m(5) * (t133 * t86 + (-t141 * t33 - t192 * t64) * t131) + (qJD(2) * t94 + t87) * mrSges(4,2) - t262 * t95 + (-t101 * t72 + t102 * t71 - t108 * t20 - t109 * t19) * mrSges(6,3) + t258 * t62 + t259 * t61 + (-t86 * mrSges(5,1) + (t111 / 0.2e1 - t63 * mrSges(5,3) + t79 * mrSges(5,2) + t203 / 0.2e1 + (-m(5) * t63 - t123) * t131 + (0.3e1 / 0.2e1 * Ifges(5,4) * t144 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t141) * qJD(2)) * qJD(4) + (-m(5) * t64 - t124) * t94 + (m(5) * t131 + mrSges(5,3)) * t32) * t144 + (t13 * t233 + t14 * t232 + Ifges(6,1) * t102 - Ifges(6,4) * t101 + t76 * mrSges(6,2) + t165 * t245 + t163 * t244 + t161 * t237 + (mrSges(6,3) + t166) * t9 + t170 * mrSges(7,3) + (t17 * t167 + t162 * t242 + t164 * t241 + t160 * t236 + t36 * t233 - t143 * t35 / 0.2e1 + t168 * mrSges(7,3)) * qJD(6)) * t116 - m(5) * (-t94 * t207 + t79 * t91) + (mrSges(4,1) * qJD(2) - t205) * t91 + (t86 * mrSges(5,2) - t33 * mrSges(5,3) + t94 * t123 + (-t110 / 0.2e1 + t79 * mrSges(5,1) - t64 * mrSges(5,3) + pkin(4) * t74 - t131 * t124 - 0.3e1 / 0.2e1 * Ifges(5,4) * t194 - t201 / 0.2e1) * qJD(4)) * t141 + t30 * t37 + t31 * t38 + t71 * t23 - t86 * mrSges(4,1) + t122 * t65 + t133 * t112 - t257 * t222 + (t1 * t31 - t17 * t257 + t2 * t30 + t258 * t5 + t259 * t6 + t238) * m(7) + (t10 * t72 + t122 * t76 + t257 * t19 - t20 * t262 + t263 * t73 + t238) * m(6) + (t83 * t91 - t84 * t94 + (-t136 * t87 - t138 * t86) * pkin(2)) * m(4); -(t23 + t210) * t153 + t222 * t108 - t156 * t109 + (-t141 * t123 + t144 * t124 + (-t141 ^ 2 - t144 ^ 2) * qJD(2) * mrSges(5,3)) * qJD(4) + (-t211 + (-t140 * t61 - t143 * t62) * qJD(6) + t251) * t116 + m(5) * (t141 * t32 + t144 * t33 + (t144 * t64 - t207) * qJD(4)) + m(6) * (t10 * t116 - t108 * t19 + t109 * t20 - t226) + m(7) * (t108 * t17 - t109 * t168 + t116 * t147 - t226); (t185 + t123) * t64 + (t184 - t124) * t63 - m(7) * (t11 * t5 + t12 * t6) + (-t187 * t73 + t19 * t21 - t20 * t22 + (t10 * t135 - t204 * t9) * pkin(4)) * m(6) + (Ifges(6,1) * t234 + t161 * t236 + t163 * t242 + t165 * t241 - t155 - t255) * t106 + t248 * qJD(6) + t188 * Ifges(5,5) * t144 / 0.2e1 + (-m(7) * t17 - t222) * t21 + ((-t190 + t199) * t6 + (-t189 + t198) * t5 + t250) * mrSges(7,3) + (m(7) * t147 - t189 * t62 - t190 * t61 + t251) * (pkin(9) + t228) + (m(7) * t132 - mrSges(6,1) - t167) * t9 - t10 * mrSges(6,2) - (-Ifges(5,2) * t194 + t111 + t134) * t193 / 0.2e1 + t160 * t237 + t162 * t244 + t164 * t245 + t19 * t221 + t13 * t232 + (-Ifges(7,5) * t241 - Ifges(7,6) * t242 - Ifges(7,3) * t236 + t247 + t260) * t107 + (t217 + t34) * t234 + (t104 + t70) * t235 - t211 * t228 - t20 * t220 - t180 * t210 - t36 * t198 / 0.2e1 + t35 * t199 / 0.2e1 + t110 * t194 / 0.2e1 - t74 * t187 - Ifges(5,6) * t177 / 0.2e1 - t32 * mrSges(5,2) + t33 * mrSges(5,1) - t12 * t61 - t11 * t62 - t22 * t95 - Ifges(6,6) * t101 + Ifges(6,5) * t102 - t141 * t246 * (Ifges(5,1) * t144 - t218) / 0.2e1 + t132 * t23 - t79 * t154 + t140 * t14 / 0.2e1; t159 * qJD(6) + t156 * t106 + t222 * t107 + t140 * t38 + t143 * t37 + t65 + (-t105 * t168 + t107 * t17 - t170) * m(7) + (-t106 * t20 - t107 * t19 + t76) * m(6); t98 - t17 * (mrSges(7,1) * t89 + mrSges(7,2) * t88) + (Ifges(7,1) * t88 - t223) * t241 + t35 * t240 + (Ifges(7,5) * t88 - Ifges(7,6) * t89) * t236 + t6 * t62 - t5 * t61 + (t5 * t88 + t6 * t89) * mrSges(7,3) + (-Ifges(7,2) * t89 + t36 + t82) * t242 + t249;];
tauc  = t7(:);
