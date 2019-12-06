% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tau_reg [5x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:46
% EndTime: 2019-12-05 18:53:52
% DurationCPUTime: 2.01s
% Computational Cost: add. (2458->262), mult. (4542->394), div. (0->0), fcn. (3452->14), ass. (0->179)
t135 = cos(qJ(4));
t132 = sin(qJ(2));
t191 = qJDD(1) * t132;
t137 = cos(qJ(2));
t197 = qJD(2) * t137;
t136 = cos(qJ(3));
t234 = pkin(1) * t136;
t147 = (qJD(1) * t197 + t191) * t234;
t131 = sin(qJ(3));
t198 = qJD(1) * t132;
t187 = pkin(1) * t198;
t171 = t131 * t187;
t81 = qJD(3) * pkin(2) - t171;
t218 = qJD(4) * t81;
t243 = (t147 + t218) * t135;
t124 = qJD(1) + qJD(2);
t174 = qJD(1) * (-qJD(2) + t124);
t128 = qJ(1) + qJ(2);
t117 = sin(t128);
t119 = cos(t128);
t200 = g(1) * t119 + g(2) * t117;
t142 = (t137 * t174 - t191) * pkin(1) + t200;
t130 = sin(qJ(4));
t122 = qJDD(1) + qJDD(2);
t202 = t136 * t122;
t205 = t131 * t122;
t123 = qJD(3) + qJD(4);
t204 = t135 * t131;
t77 = t130 * t136 + t204;
t241 = t123 * t77;
t30 = t124 * t241 + t130 * t205 - t135 * t202;
t129 = sin(qJ(5));
t134 = cos(qJ(5));
t192 = qJD(5) * t134;
t163 = qJD(3) * t171;
t194 = qJD(3) * t136;
t45 = qJDD(3) * pkin(2) + (-t131 * t191 + (-t131 * t197 - t132 * t194) * qJD(1)) * pkin(1);
t226 = -t130 * t163 - t135 * t45;
t238 = t130 * t218 + ((qJD(4) * t132 * t135 + t130 * t197) * qJD(1) + t130 * t191) * t234;
t20 = t226 + t238;
t172 = t136 * t187;
t55 = t130 * t172 - t135 * t81;
t242 = t20 * t129 + t55 * t192;
t112 = g(1) * t117;
t190 = qJDD(1) * t137;
t115 = pkin(1) * t190;
t237 = t112 + t115;
t127 = qJ(3) + qJ(4);
t118 = cos(t127);
t116 = sin(t127);
t216 = t116 * t119;
t217 = t116 * t117;
t239 = g(1) * t216 + g(2) * t217 - g(3) * t118;
t121 = qJDD(3) + qJDD(4);
t203 = t135 * t136;
t100 = t124 * t203;
t178 = t124 * t194;
t207 = t130 * t131;
t180 = t124 * t207;
t29 = qJD(4) * t100 + t122 * t204 - t123 * t180 + t130 * t202 + t135 * t178;
t70 = t77 * t124;
t51 = t129 * t123 + t134 * t70;
t17 = t51 * qJD(5) - t134 * t121 + t129 * t29;
t235 = pkin(1) * t132;
t233 = g(2) * t119;
t110 = g(3) * t116;
t231 = t136 * pkin(2);
t49 = -t134 * t123 + t129 * t70;
t68 = -t100 + t180;
t61 = qJD(5) + t68;
t230 = t49 * t61;
t229 = t51 * t61;
t228 = t61 * t70;
t227 = t70 * t68;
t225 = pkin(1) * qJD(1);
t28 = qJDD(5) + t30;
t224 = t129 * t28;
t223 = t129 * t55;
t222 = t134 * t28;
t193 = qJD(5) * t129;
t16 = t129 * t121 + t123 * t192 + t134 * t29 - t70 * t193;
t221 = t16 * t129;
t186 = qJD(2) * t235;
t170 = qJD(1) * t186;
t220 = (t170 + t233) * t131;
t219 = qJD(4) * t61;
t215 = t117 * t118;
t214 = t117 * t134;
t213 = t118 * t119;
t212 = t118 * t129;
t211 = t119 * t129;
t210 = t119 * t134;
t209 = t122 * t137;
t208 = t124 * t131;
t201 = t237 * t136;
t125 = t131 ^ 2;
t199 = -t136 ^ 2 + t125;
t196 = qJD(3) * t124;
t195 = qJD(3) * t131;
t189 = qJDD(3) * t136;
t188 = -qJDD(1) - t122;
t185 = t61 * t208;
t183 = t77 * t193;
t182 = t61 * t192;
t181 = qJD(5) * t136 * t61;
t46 = t55 * t193;
t179 = (-t61 + t68) * t55;
t175 = t134 * t61;
t173 = (-qJD(1) - t124) * qJD(2);
t162 = qJD(4) * t172;
t169 = -t135 * t163 + (-t162 + t45) * t130;
t76 = -t203 + t207;
t52 = t123 * t76;
t166 = t28 * t77 - t52 * t61;
t59 = t77 * t187;
t165 = t55 * t68 - t59 * t61;
t83 = t135 * t172;
t56 = t130 * t81 + t83;
t73 = -t124 * t231 - t137 * t225;
t37 = t129 * t73 + t134 * t56;
t36 = -t129 * t56 + t134 * t73;
t164 = -0.2e1 * t196;
t101 = -t137 * pkin(1) - t231;
t67 = t76 * t235;
t161 = t134 * t101 + t129 * t67;
t160 = t129 * t101 - t134 * t67;
t19 = t169 + t243;
t159 = -qJD(5) * t73 + t110 - t19;
t48 = t170 - t115 + (t124 * t195 - t202) * pkin(2);
t158 = -g(1) * t217 + g(2) * t216 + t48 * t77 - t73 * t52;
t157 = g(1) * t215 - g(2) * t213 + t241 * t73 + t48 * t76;
t156 = -t136 * t28 + t61 * t195;
t155 = t77 * t137;
t154 = t76 * t137;
t62 = t117 * t212 + t210;
t64 = -t118 * t211 + t214;
t151 = -g(1) * t62 - g(2) * t64 - t37 * t241 - (t36 * qJD(5) + t129 * t48) * t76 + (-t19 * t76 + t20 * t77 - t55 * t52) * t134;
t150 = -t73 * t70 - t226 + t239;
t149 = t174 * t235 - t233;
t44 = t134 * t48;
t63 = -t118 * t214 + t211;
t65 = t117 * t129 + t118 * t210;
t148 = -g(1) * t63 - g(2) * t65 + t36 * t241 - t52 * t223 + (-t37 * qJD(5) - t129 * t19 + t44) * t76 + t242 * t77;
t146 = g(1) * t213 + g(2) * t215 + t73 * t68 + t110 - t169;
t144 = -t36 * t70 + t46 + (-t20 + t239) * t134;
t143 = -t116 * t129 * t200 + g(3) * t212 + t37 * t70 + t242;
t140 = (-pkin(2) * t123 - t81) * qJD(4) - t147;
t139 = qJD(3) ^ 2;
t138 = cos(qJ(1));
t133 = sin(qJ(1));
t120 = t124 ^ 2;
t93 = -t139 * t131 + t189;
t92 = qJDD(3) * t131 + t139 * t136;
t82 = pkin(2) * t195 + t186;
t71 = t125 * t122 + 0.2e1 * t131 * t178;
t66 = t77 * t235;
t60 = t154 * t225;
t58 = t155 * t225;
t57 = t130 * t171 - t83;
t54 = 0.2e1 * t131 * t202 - 0.2e1 * t199 * t196;
t35 = -t76 * t121 - t123 * t241;
t34 = t77 * t121 - t52 * t123;
t31 = -t68 ^ 2 + t70 ^ 2;
t27 = (qJD(2) * t155 - t52 * t132) * pkin(1);
t26 = (-qJD(2) * t154 - t132 * t241) * pkin(1);
t22 = t70 * t123 - t30;
t21 = t68 * t123 + t29;
t13 = t29 * t77 - t70 * t52;
t12 = t241 * t61 + t28 * t76;
t11 = t61 * t175 - t51 * t70 + t224;
t10 = -t61 ^ 2 * t129 + t49 * t70 + t222;
t9 = t51 * t175 + t221;
t6 = -t51 * t183 + (t16 * t77 - t51 * t52) * t134;
t5 = -t241 * t70 - t29 * t76 - t77 * t30 + t52 * t68;
t4 = -t166 * t129 - t17 * t76 - t77 * t182 - t241 * t49;
t3 = t166 * t134 + t16 * t76 - t61 * t183 + t241 * t51;
t2 = (t16 - t230) * t134 + (-t17 - t229) * t129;
t1 = -(-t129 * t51 - t134 * t49) * t52 + (-t221 - t134 * t17 + (t129 * t49 - t134 * t51) * qJD(5)) * t77;
t7 = [qJDD(1), g(1) * t133 - g(2) * t138, g(1) * t138 + g(2) * t133, t122, -t233 + (t132 * t173 + t209) * pkin(1) + t237, (t132 * t188 + t137 * t173) * pkin(1) + t200, t71, t54, t92, t93, 0, -t136 * t233 + ((t209 + (-t139 + t173) * t132) * t136 + (-qJDD(3) * t132 + t137 * t164) * t131) * pkin(1) + t201, -t131 * t112 + ((-t189 + (qJD(2) * t124 + t139) * t131) * t132 + (t131 * t188 + t136 * t164) * t137) * pkin(1) + t220, t13, t5, t34, t35, 0, t101 * t30 - t66 * t121 - t27 * t123 + t82 * t68 + t157, t101 * t29 + t67 * t121 - t26 * t123 + t82 * t70 + t158, t6, t1, t3, t4, t12, (-qJD(5) * t160 - t129 * t26 + t134 * t82) * t61 + t161 * t28 + t27 * t49 + t66 * t17 + t148, -(t129 * t82 + t134 * t26) * t61 - t160 * t28 + t27 * t51 + t66 * t16 + (-t161 * t61 - t77 * t223) * qJD(5) + t151; 0, 0, 0, t122, t149 + t237, t142, t71, t54, t92, t93, 0, t136 * t149 + t201, (-t112 + (-t124 * t198 - t190) * pkin(1)) * t131 + t220, t13, t5, t34, t35, 0, -t68 * t187 + t58 * t123 + (-t136 * t30 + t68 * t195) * pkin(2) + t157, -t70 * t187 - t60 * t123 + (-t136 * t29 + t70 * t195) * pkin(2) + t158, t6, t1, t3, t4, t12, -(t129 * t60 + t134 * t187) * t61 - t58 * t49 + (t129 * t181 + t134 * t156) * pkin(2) + t148, -t77 * t46 + (t129 * t187 - t134 * t60) * t61 - t58 * t51 + (-t129 * t156 + t134 * t181) * pkin(2) + t151; 0, 0, 0, 0, 0, 0, -t131 * t120 * t136, t199 * t120, t205, t202, qJDD(3), -g(3) * t136 + t131 * t142, g(3) * t131 + t136 * t142, t227, t31, t21, t22, t121, -t135 * t162 - t57 * t123 + (t121 * t135 - t68 * t208) * pkin(2) + t140 * t130 + t150, -t59 * t123 + (-t121 * t130 - t70 * t208) * pkin(2) + t140 * t135 + t146, t9, t2, t11, t10, -t228, t57 * t49 + t165 * t129 + (-t134 * t185 + (-t129 * t219 - t17) * t135 + (qJD(4) * t49 - t182 - t224) * t130) * pkin(2) + t144, t57 * t51 + t165 * t134 + (t129 * t185 + (-t134 * t219 - t16) * t135 + (qJD(4) * t51 + t61 * t193 - t222) * t130) * pkin(2) + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t31, t21, t22, t121, t56 * t123 + t150 - t238, -t55 * t123 + t146 - t243, t9, t2, t11, t10, -t228, t129 * t179 - t56 * t49 + t144, t134 * t179 - t56 * t51 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t49, -t49 ^ 2 + t51 ^ 2, t16 + t230, -t17 + t229, t28, -g(1) * t64 + g(2) * t62 + t129 * t159 - t56 * t192 + t37 * t61 - t55 * t51 + t44, g(1) * t65 - g(2) * t63 + t36 * t61 + t55 * t49 + (qJD(5) * t56 - t48) * t129 + t159 * t134;];
tau_reg = t7;
