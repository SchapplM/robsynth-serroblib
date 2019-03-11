% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:41
% EndTime: 2019-03-09 02:42:47
% DurationCPUTime: 2.38s
% Computational Cost: add. (2697->327), mult. (5737->411), div. (0->0), fcn. (3998->14), ass. (0->180)
t126 = sin(pkin(10));
t128 = cos(pkin(10));
t132 = sin(qJ(3));
t135 = cos(qJ(3));
t88 = t126 * t135 + t128 * t132;
t82 = t88 * qJD(1);
t234 = qJD(6) + t82;
t131 = sin(qJ(6));
t235 = t131 * t234;
t134 = cos(qJ(6));
t189 = qJD(1) * qJD(3);
t183 = t135 * t189;
t184 = t132 * t189;
t54 = qJDD(1) * t88 - t126 * t184 + t128 * t183;
t52 = qJDD(6) + t54;
t46 = t134 * t52;
t159 = -t234 * t235 + t46;
t190 = t131 * qJD(3);
t193 = qJD(1) * t132;
t200 = t128 * t135;
t79 = -qJD(1) * t200 + t126 * t193;
t61 = -t134 * t79 + t190;
t236 = t234 * t61;
t127 = sin(pkin(9));
t106 = t127 * pkin(1) + pkin(7);
t199 = qJ(4) + t106;
t77 = t82 ^ 2;
t233 = -t79 ^ 2 - t77;
t178 = t199 * qJD(1);
t65 = t132 * qJD(2) + t135 * t178;
t56 = t126 * t65;
t64 = t135 * qJD(2) - t132 * t178;
t35 = t128 * t64 - t56;
t198 = -qJD(5) + t35;
t232 = -qJD(6) + t234;
t122 = qJ(1) + pkin(9);
t113 = sin(t122);
t115 = cos(t122);
t171 = g(1) * t113 - g(2) * t115;
t129 = cos(pkin(9));
t108 = -t129 * pkin(1) - pkin(2);
t119 = t135 * pkin(3);
t231 = t108 - t119;
t121 = qJ(3) + pkin(10);
t112 = sin(t121);
t114 = cos(t121);
t230 = t114 * pkin(4) + t112 * qJ(5);
t172 = g(1) * t115 + g(2) * t113;
t107 = -t128 * pkin(3) - pkin(4);
t100 = -pkin(8) + t107;
t224 = t79 * pkin(5);
t213 = t128 * t65;
t59 = qJD(3) * pkin(3) + t64;
t31 = t126 * t59 + t213;
t25 = -qJD(3) * qJ(5) - t31;
t14 = -t25 - t224;
t34 = t126 * t64 + t213;
t229 = t100 * t52 + (t14 - t34 + t224) * t234;
t191 = qJD(6) * t134;
t87 = t126 * t132 - t200;
t185 = t87 * t191;
t81 = t88 * qJD(3);
t228 = -(t234 * t81 + t52 * t87) * t131 - t234 * t185;
t226 = pkin(4) + pkin(8);
t187 = t135 * qJDD(1);
t188 = t132 * qJDD(1);
t166 = -t126 * t188 + t128 * t187;
t53 = qJD(1) * t81 - t166;
t225 = t53 * pkin(4);
t223 = t82 * pkin(5);
t103 = g(3) * t114;
t219 = g(3) * t135;
t156 = -t88 * qJ(5) + t231;
t36 = t226 * t87 + t156;
t217 = t36 * t52;
t63 = t134 * qJD(3) + t131 * t79;
t216 = t63 * t79;
t215 = t79 * t61;
t18 = -qJD(6) * t190 + t134 * qJDD(3) + t131 * t53 + t79 * t191;
t192 = qJD(3) * t132;
t84 = qJD(3) * t200 - t126 * t192;
t214 = t18 * t88 + t63 * t84;
t117 = t135 * qJDD(2);
t92 = t106 * qJDD(1);
t148 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(3) * qJD(2) + t92;
t162 = t178 * qJD(3);
t28 = qJDD(3) * pkin(3) - t132 * t148 - t135 * t162 + t117;
t32 = (qJDD(2) - t162) * t132 + t148 * t135;
t10 = t126 * t28 + t128 * t32;
t211 = t134 * t234;
t210 = t18 * t134;
t209 = t54 * qJ(5);
t109 = t119 + pkin(2);
t136 = cos(qJ(1));
t208 = t136 * pkin(1) + t115 * t109;
t207 = qJD(6) * t87;
t206 = qJDD(3) * pkin(4);
t204 = t113 * t131;
t203 = t113 * t134;
t202 = t115 * t131;
t201 = t115 * t134;
t197 = t223 - t198;
t195 = qJDD(2) - g(3);
t124 = t132 ^ 2;
t194 = -t135 ^ 2 + t124;
t95 = qJD(1) * t108;
t186 = qJDD(3) * qJ(5) + t10;
t111 = pkin(3) * t192;
t30 = t128 * t59 - t56;
t173 = qJD(5) - t30;
t12 = -t226 * qJD(3) + t173 + t223;
t147 = pkin(3) * t184 + qJDD(1) * t231 + qJDD(4);
t141 = -t82 * qJD(5) + t147 - t209;
t6 = t226 * t53 + t141;
t182 = qJD(6) * t12 + t6;
t9 = -t126 * t32 + t128 * t28;
t167 = qJDD(5) - t9;
t2 = t54 * pkin(5) - t226 * qJDD(3) + t167;
t78 = qJD(1) * t231 + qJD(4);
t145 = -t82 * qJ(5) + t78;
t22 = t226 * t79 + t145;
t181 = -qJD(6) * t22 + t2;
t177 = qJD(3) * t199;
t66 = t135 * qJD(4) - t132 * t177;
t67 = -t132 * qJD(4) - t135 * t177;
t38 = t126 * t66 - t128 * t67;
t85 = t199 * t132;
t86 = t199 * t135;
t49 = t126 * t86 + t128 * t85;
t180 = pkin(3) * t193 + t79 * qJ(5);
t176 = t131 * qJDD(3) - t134 * t53;
t133 = sin(qJ(1));
t170 = g(1) * t133 - g(2) * t136;
t19 = qJD(6) * t63 + t176;
t169 = -t88 * t19 - t84 * t61;
t130 = -qJ(4) - pkin(7);
t164 = -t133 * pkin(1) - t115 * t130;
t5 = t131 * t12 + t134 * t22;
t39 = t126 * t67 + t128 * t66;
t50 = -t126 * t85 + t128 * t86;
t7 = -qJD(3) * qJD(5) - t186;
t158 = -t84 * qJ(5) - t88 * qJD(5) + t111;
t157 = -t207 * t235 + t81 * t211 + t87 * t46;
t3 = -t53 * pkin(5) - t7;
t41 = t88 * pkin(5) + t49;
t155 = t14 * t81 + t3 * t87 - t41 * t52;
t153 = -t131 * t52 - t211 * t234;
t152 = -t95 * qJD(1) + t172 - t92;
t151 = 0.2e1 * qJD(3) * t95 - qJDD(3) * t106;
t150 = -g(3) * t112 - t114 * t172;
t149 = -t88 * t53 + t87 * t54 - t84 * t79 + t81 * t82;
t137 = qJD(3) ^ 2;
t146 = -0.2e1 * qJDD(1) * t108 - t106 * t137 + t171;
t144 = t147 - t171;
t143 = t38 * t82 - t39 * t79 + t49 * t54 - t50 * t53 - t172;
t40 = t79 * pkin(4) + t145;
t142 = -t112 * t172 + t40 * t82 + t103 + t167;
t140 = t3 + (-qJD(6) * t100 + t226 * t82 + t180) * t234 + t150;
t138 = qJD(1) ^ 2;
t102 = t126 * pkin(3) + qJ(5);
t91 = qJDD(3) * t135 - t137 * t132;
t90 = qJDD(3) * t132 + t137 * t135;
t75 = qJD(3) * t79;
t71 = -t112 * t204 + t201;
t70 = t112 * t203 + t202;
t69 = t112 * t202 + t203;
t68 = t112 * t201 - t204;
t45 = t87 * pkin(4) + t156;
t43 = t82 * pkin(4) + t180;
t42 = -t87 * pkin(5) + t50;
t33 = t81 * pkin(4) + t158;
t24 = -qJD(3) * pkin(4) + t173;
t21 = -t81 * pkin(5) + t39;
t20 = t84 * pkin(5) + t38;
t15 = t226 * t81 + t158;
t11 = t141 + t225;
t8 = t167 - t206;
t4 = t134 * t12 - t131 * t22;
t1 = t134 * t2;
t13 = [qJDD(1), t170, g(1) * t136 + g(2) * t133 (t170 + (t127 ^ 2 + t129 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t124 * qJDD(1) + 0.2e1 * t132 * t183, 0.2e1 * t132 * t187 - 0.2e1 * t189 * t194, t90, t91, 0, t132 * t151 + t135 * t146, -t132 * t146 + t135 * t151, -t10 * t87 - t30 * t84 - t31 * t81 - t9 * t88 + t143, t10 * t50 + t31 * t39 - t9 * t49 - t30 * t38 + t147 * t231 + t78 * t111 - g(1) * (-t113 * t109 + t164) - g(2) * (-t113 * t130 + t208) t24 * t84 + t25 * t81 + t7 * t87 + t8 * t88 + t143, t38 * qJD(3) + t49 * qJDD(3) - t11 * t87 - t114 * t171 - t33 * t79 - t40 * t81 - t45 * t53, t39 * qJD(3) + t50 * qJDD(3) - t11 * t88 + t112 * t171 - t33 * t82 - t40 * t84 - t45 * t54, t11 * t45 + t40 * t33 - t7 * t50 - t25 * t39 + t8 * t49 + t24 * t38 - g(1) * t164 - g(2) * (t230 * t115 + t208) + (-g(1) * (-t109 - t230) + g(2) * t130) * t113, t63 * t185 + (t18 * t87 + t63 * t81) * t131 (-t131 * t61 + t134 * t63) * t81 + (-t131 * t19 + t210 + (-t131 * t63 - t134 * t61) * qJD(6)) * t87, t214 - t228, t157 + t169, t234 * t84 + t52 * t88, -g(1) * t71 - g(2) * t69 + t1 * t88 + t42 * t19 + t21 * t61 + t4 * t84 + (-t15 * t234 - t6 * t88 - t217) * t131 + (t20 * t234 - t155) * t134 + ((-t131 * t41 - t134 * t36) * t234 - t5 * t88 + t14 * t131 * t87) * qJD(6), g(1) * t70 - g(2) * t68 + t42 * t18 + t21 * t63 - t5 * t84 + (-(qJD(6) * t41 + t15) * t234 - t217 - t182 * t88 + t14 * t207) * t134 + (-(-qJD(6) * t36 + t20) * t234 - t181 * t88 + t155) * t131; 0, 0, 0, t195, 0, 0, 0, 0, 0, t91, -t90, t149, t10 * t88 - t30 * t81 + t31 * t84 - t9 * t87 - g(3), t149, t81 * qJD(3) + t87 * qJDD(3), t84 * qJD(3) + t88 * qJDD(3), t24 * t81 - t25 * t84 - t7 * t88 + t8 * t87 - g(3), 0, 0, 0, 0, 0, t157 - t169, t214 + t228; 0, 0, 0, 0, -t132 * t138 * t135, t194 * t138, t188, t187, qJDD(3), t132 * t152 + t117 - t219, -t132 * t195 + t135 * t152 (t31 - t34) * t82 + (-t30 + t35) * t79 + (-t126 * t53 - t128 * t54) * pkin(3), t30 * t34 - t31 * t35 + (-t219 + t10 * t126 + t128 * t9 + (-qJD(1) * t78 + t172) * t132) * pkin(3), -t102 * t53 + t107 * t54 + (-t25 - t34) * t82 + (t24 + t198) * t79, -t34 * qJD(3) + t43 * t79 + (-pkin(4) + t107) * qJDD(3) + t142, t102 * qJDD(3) - t40 * t79 + t43 * t82 + (0.2e1 * qJD(5) - t35) * qJD(3) + t150 + t186, -t7 * t102 + t8 * t107 - t40 * t43 - t24 * t34 - g(3) * (t119 + t230) + t198 * t25 + t172 * (pkin(3) * t132 + pkin(4) * t112 - qJ(5) * t114) -t235 * t63 + t210 (-t234 * t63 - t19) * t134 + (-t18 + t236) * t131, t159 + t216, t153 - t215, t234 * t79, t102 * t19 + t140 * t131 + t229 * t134 + t197 * t61 + t4 * t79, t102 * t18 - t229 * t131 + t140 * t134 + t197 * t63 - t5 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t30 * t82 + t31 * t79 + t144, t233, -0.2e1 * t82 * qJD(3) + t166, -t54 + t75, t225 - t209 - t25 * t79 + (-qJD(5) - t24) * t82 + t144, 0, 0, 0, 0, 0, t153 + t215, -t159 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 + t54, -t82 * t79 + qJDD(3), -t77 - t137, t25 * qJD(3) + t142 - t206, 0, 0, 0, 0, 0, -qJD(3) * t61 + t159, -qJD(3) * t63 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t61, -t61 ^ 2 + t63 ^ 2, t18 + t236, t232 * t63 - t176, t52, -g(1) * t68 - g(2) * t70 + t134 * t103 - t131 * t6 - t14 * t63 + t232 * t5 + t1, g(1) * t69 - g(2) * t71 + t14 * t61 + t4 * t234 - t182 * t134 + (-t181 - t103) * t131;];
tau_reg  = t13;
