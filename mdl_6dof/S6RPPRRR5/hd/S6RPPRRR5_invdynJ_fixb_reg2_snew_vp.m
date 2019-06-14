% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:50:37
% EndTime: 2019-05-05 15:50:45
% DurationCPUTime: 2.80s
% Computational Cost: add. (8645->300), mult. (17056->382), div. (0->0), fcn. (10689->8), ass. (0->201)
t166 = sin(qJ(5));
t160 = qJDD(4) + qJDD(5);
t167 = sin(qJ(4));
t170 = cos(qJ(5));
t171 = cos(qJ(4));
t216 = t171 * t166;
t134 = (t167 * t170 + t216) * qJD(1);
t213 = qJD(1) * t171;
t136 = -t166 * t167 * qJD(1) + t170 * t213;
t225 = t136 * t134;
t245 = t160 - t225;
t247 = t166 * t245;
t246 = t170 * t245;
t161 = qJD(4) + qJD(5);
t224 = t161 * t134;
t210 = qJD(1) * qJD(4);
t198 = t171 * t210;
t208 = t167 * qJDD(1);
t141 = -t198 - t208;
t199 = t167 * t210;
t207 = t171 * qJDD(1);
t142 = -t199 + t207;
t99 = -t134 * qJD(5) + t166 * t141 + t170 * t142;
t195 = -t99 + t224;
t165 = sin(qJ(6));
t169 = cos(qJ(6));
t117 = t165 * t136 - t169 * t161;
t119 = t169 * t136 + t165 * t161;
t96 = t119 * t117;
t192 = -t170 * t141 + t166 * t142;
t98 = -t136 * qJD(5) - t192;
t97 = qJDD(6) - t98;
t241 = -t96 + t97;
t244 = t165 * t241;
t243 = t169 * t241;
t110 = t134 * pkin(5) - t136 * pkin(9);
t239 = t161 ^ 2;
t174 = qJD(1) ^ 2;
t168 = sin(qJ(1));
t172 = cos(qJ(1));
t190 = t172 * g(1) + t168 * g(2);
t185 = 0.2e1 * qJD(2) * qJD(1) - t190;
t184 = qJDD(3) + t185;
t235 = pkin(1) + qJ(3);
t180 = -t235 * t174 + t184;
t234 = -pkin(7) + qJ(2);
t178 = t234 * qJDD(1) + t180;
t109 = t171 * g(3) - t167 * t178;
t147 = qJD(4) * pkin(4) - pkin(8) * t213;
t162 = t167 ^ 2;
t221 = t162 * t174;
t92 = -pkin(4) * t221 + t141 * pkin(8) - qJD(4) * t147 - t109;
t230 = t170 * t92;
t177 = t171 * t178;
t176 = -t142 * pkin(8) + t177;
t215 = t171 * t174;
t240 = qJDD(4) * pkin(4) + t176 + (-pkin(4) * t215 - pkin(8) * t210 + g(3)) * t167;
t57 = t240 * t166 + t230;
t39 = -t239 * pkin(5) + t160 * pkin(9) - t134 * t110 + t57;
t164 = t174 * pkin(7);
t197 = t168 * g(1) - t172 * g(2);
t189 = -qJDD(2) + t197;
t183 = t174 * qJ(2) + t189;
t193 = t235 * qJDD(1);
t179 = t193 + t183;
t238 = 2 * qJD(3);
t94 = -t141 * pkin(4) - pkin(8) * t221 - t164 + (t147 * t171 + t238) * qJD(1) + t179;
t40 = t195 * pkin(9) + (t161 * t136 - t98) * pkin(5) + t94;
t16 = t165 * t39 - t169 * t40;
t17 = t165 * t40 + t169 * t39;
t7 = t165 * t16 + t169 * t17;
t242 = t234 - pkin(8);
t131 = qJD(6) + t134;
t194 = -t169 * t160 + t165 * t99;
t64 = (qJD(6) - t131) * t119 + t194;
t115 = t117 ^ 2;
t116 = t119 ^ 2;
t129 = t131 ^ 2;
t132 = t134 ^ 2;
t133 = t136 ^ 2;
t56 = t166 * t92 - t170 * t240;
t38 = -t160 * pkin(5) - t239 * pkin(9) + t136 * t110 + t56;
t237 = -pkin(5) * t38 + pkin(9) * t7;
t236 = t167 * g(3);
t34 = t165 * t38;
t71 = t96 + t97;
t233 = t165 * t71;
t232 = t166 * t94;
t35 = t169 * t38;
t231 = t169 * t71;
t229 = t170 * t94;
t228 = qJDD(1) * pkin(1);
t227 = t131 * t165;
t226 = t131 * t169;
t223 = t161 * t166;
t222 = t161 * t170;
t163 = t171 ^ 2;
t220 = t163 * t174;
t106 = t225 + t160;
t219 = t166 * t106;
t202 = t167 * t215;
t218 = t167 * (qJDD(4) + t202);
t217 = t170 * t106;
t214 = t162 + t163;
t211 = qJD(6) + t131;
t209 = qJDD(1) * qJ(2);
t93 = -t116 - t129;
t46 = -t165 * t93 - t231;
t187 = -t165 * t160 - t169 * t99;
t69 = t211 * t117 + t187;
t206 = pkin(5) * t69 + pkin(9) * t46 + t34;
t82 = -t129 - t115;
t43 = t169 * t82 - t244;
t65 = -t211 * t119 - t194;
t205 = pkin(5) * t65 + pkin(9) * t43 - t35;
t204 = qJD(1) * t238;
t203 = t166 * t96;
t201 = t170 * t96;
t200 = -pkin(5) * t170 - pkin(4);
t28 = t166 * t56 + t170 * t57;
t103 = t131 * t117;
t76 = -t117 * qJD(6) - t187;
t68 = t103 + t76;
t33 = t165 * t68 - t169 * t64;
t78 = t115 + t116;
t196 = pkin(5) * t78 + pkin(9) * t33 + t7;
t191 = t167 * pkin(4) + t235;
t3 = t166 * t7 - t170 * t38;
t1 = t167 * (t166 * t38 + t170 * t7) + t171 * t3;
t188 = qJDD(4) - t202;
t6 = -t169 * t16 + t165 * t17;
t27 = t166 * t57 - t170 * t56;
t13 = t167 * t28 + t171 * t27;
t108 = t177 + t236;
t81 = t171 * t108 - t167 * t109;
t182 = (-qJD(5) + t161) * t136 - t192;
t123 = t179 + t204;
t173 = qJD(4) ^ 2;
t156 = 0.2e1 * t209;
t145 = t214 * t174;
t144 = t214 * qJDD(1);
t143 = -0.2e1 * t199 + t207;
t140 = 0.2e1 * t198 + t208;
t139 = t171 * t188;
t130 = t183 + t228;
t126 = -t133 + t239;
t125 = t132 - t239;
t124 = t180 + t209;
t122 = -t133 - t239;
t121 = -t218 + t171 * (-t173 - t220);
t120 = t167 * (-t173 - t221) + t139;
t114 = -t164 + t123;
t111 = t133 - t132;
t104 = -t239 - t132;
t102 = -t116 + t129;
t101 = t115 - t129;
t100 = -t132 - t133;
t95 = t116 - t115;
t91 = -t166 * t122 - t217;
t90 = t170 * t122 - t219;
t89 = t224 + t99;
t84 = (qJD(5) + t161) * t136 + t192;
t80 = t170 * t104 - t247;
t79 = t166 * t104 + t246;
t75 = -t119 * qJD(6) - t194;
t74 = (-t117 * t169 + t119 * t165) * t131;
t73 = (-t117 * t165 - t119 * t169) * t131;
t67 = -t103 + t76;
t61 = -t119 * t227 + t169 * t76;
t60 = t119 * t226 + t165 * t76;
t59 = t117 * t226 - t165 * t75;
t58 = t117 * t227 + t169 * t75;
t55 = t167 * t91 + t171 * t90;
t54 = t166 * t89 + t170 * t182;
t53 = t166 * t182 - t170 * t89;
t51 = t169 * t101 - t233;
t50 = -t165 * t102 + t243;
t49 = t165 * t101 + t231;
t48 = t169 * t102 + t244;
t47 = t167 * t80 + t171 * t79;
t45 = t169 * t93 - t233;
t42 = t165 * t82 + t243;
t32 = -t165 * t67 + t169 * t65;
t31 = -t165 * t64 - t169 * t68;
t30 = t165 * t65 + t169 * t67;
t26 = t167 * t54 + t171 * t53;
t25 = -t166 * t69 + t170 * t46;
t24 = t166 * t46 + t170 * t69;
t23 = -t166 * t65 + t170 * t43;
t22 = t166 * t43 + t170 * t65;
t21 = -t166 * t78 + t170 * t33;
t20 = t166 * t33 + t170 * t78;
t19 = -pkin(9) * t45 + t35;
t18 = -pkin(9) * t42 + t34;
t12 = -pkin(5) * t45 + t17;
t11 = -pkin(5) * t42 + t16;
t10 = t167 * t25 + t171 * t24;
t9 = t167 * t23 + t171 * t22;
t8 = t167 * t21 + t171 * t20;
t2 = -pkin(9) * t31 - t6;
t4 = [0, 0, 0, 0, 0, qJDD(1), t197, t190, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t189 - 0.2e1 * t228, t156 + t185, pkin(1) * t130 + qJ(2) * (-t174 * pkin(1) + t185 + t209), qJDD(1), 0, 0, 0, 0, 0, 0, t156 + t184, t189 + 0.2e1 * t193 + t204, qJ(2) * t124 + t235 * t123, (t142 - t199) * t171, -t171 * t140 - t167 * t143, t139 - t167 * (t173 - t220), (-t141 + t198) * t167, t171 * (-t173 + t221) - t218, 0, t167 * t114 + t234 * t120 + t235 * t140, t171 * t114 + t234 * t121 + t235 * t143, -t234 * t144 - t235 * t145 - t81, t235 * t114 + t234 * t81, t171 * (-t136 * t223 + t170 * t99) - t167 * (t136 * t222 + t166 * t99), t171 * (t166 * t195 - t170 * t84) - t167 * (-t166 * t84 - t170 * t195), t171 * (-t166 * t126 + t246) - t167 * (t170 * t126 + t247), t171 * (t134 * t222 - t166 * t98) - t167 * (t134 * t223 + t170 * t98), t171 * (t170 * t125 - t219) - t167 * (t166 * t125 + t217), (t171 * (-t134 * t170 + t136 * t166) - t167 * (-t134 * t166 - t136 * t170)) * t161, t171 * (-pkin(8) * t79 + t232) - t167 * (pkin(8) * t80 - t229) + t191 * t84 + t234 * t47, t171 * (-pkin(8) * t90 + t229) - t167 * (pkin(8) * t91 + t232) - t191 * t195 + t234 * t55, t171 * (-pkin(8) * t53 - t27) - t167 * (pkin(8) * t54 + t28) + t234 * t26 + t191 * t100, t242 * t13 + t191 * t94, t171 * (t170 * t61 + t203) - t167 * (t166 * t61 - t201), t171 * (t166 * t95 + t170 * t32) - t167 * (t166 * t32 - t170 * t95), t171 * (t166 * t68 + t170 * t50) - t167 * (t166 * t50 - t170 * t68), t171 * (t170 * t59 - t203) - t167 * (t166 * t59 + t201), t171 * (-t166 * t64 + t170 * t51) - t167 * (t166 * t51 + t170 * t64), t171 * (t166 * t97 + t170 * t74) - t167 * (t166 * t74 - t170 * t97), t171 * (-pkin(8) * t22 - t166 * t11 + t170 * t18) - t167 * (pkin(8) * t23 + t170 * t11 + t166 * t18) + t234 * t9 + t191 * t42, t171 * (-pkin(8) * t24 - t166 * t12 + t170 * t19) - t167 * (pkin(8) * t25 + t170 * t12 + t166 * t19) + t191 * t45 + t234 * t10, t171 * (-pkin(8) * t20 + t170 * t2) - t167 * (pkin(8) * t21 + t166 * t2) + t234 * t8 + (pkin(5) * t216 - t167 * t200 + t235) * t31, (t171 * (pkin(5) * t166 - pkin(9) * t170) - t167 * (-pkin(9) * t166 + t200) + t235) * t6 + t242 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t174, -t130, 0, 0, 0, 0, 0, 0, 0, -t174, -qJDD(1), -t123, 0, 0, 0, 0, 0, 0, -t140, -t143, t145, -t114, 0, 0, 0, 0, 0, 0, -t84, t195, -t100, -t94, 0, 0, 0, 0, 0, 0, -t42, -t45, -t31, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t174, t124, 0, 0, 0, 0, 0, 0, t120, t121, -t144, t81, 0, 0, 0, 0, 0, 0, t47, t55, t26, t13, 0, 0, 0, 0, 0, 0, t9, t10, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, (-t162 + t163) * t174, t207, -t202, -t208, qJDD(4), t108, t109, 0, 0, t225, t111, t89, -t225, t182, t160, pkin(4) * t79 - t56, -t230 - t166 * (-pkin(8) * t199 + t176 + t236) + (-t166 * t188 + t90) * pkin(4), pkin(4) * t53, pkin(4) * t27, t60, t30, t48, t58, t49, t73, pkin(4) * t22 + t205, pkin(4) * t24 + t206, pkin(4) * t20 + t196, pkin(4) * t3 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, t111, t89, -t225, t182, t160, -t56, -t57, 0, 0, t60, t30, t48, t58, t49, t73, t205, t206, t196, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t95, t68, -t96, -t64, t97, -t16, -t17, 0, 0;];
tauJ_reg  = t4;
