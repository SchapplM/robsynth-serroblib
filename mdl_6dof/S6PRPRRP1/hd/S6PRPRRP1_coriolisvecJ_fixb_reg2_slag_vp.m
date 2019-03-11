% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:29
% EndTime: 2019-03-08 19:58:37
% DurationCPUTime: 3.05s
% Computational Cost: add. (3845->352), mult. (9776->488), div. (0->0), fcn. (7380->10), ass. (0->203)
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t158 = pkin(4) * t131 - pkin(9) * t134;
t109 = t158 * qJD(4);
t127 = sin(pkin(6));
t126 = sin(pkin(11));
t128 = cos(pkin(11));
t132 = sin(qJ(2));
t135 = cos(qJ(2));
t150 = t126 * t135 + t128 * t132;
t87 = t150 * t127;
t80 = qJD(1) * t87;
t249 = t109 - t80;
t183 = t134 * qJD(2);
t116 = -qJD(5) + t183;
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t187 = qJD(5) * t133;
t188 = qJD(5) * t130;
t129 = cos(pkin(6));
t115 = t129 * qJD(1) + qJD(3);
t193 = qJD(1) * t127;
t172 = t135 * t193;
t110 = qJD(2) * pkin(2) + t172;
t179 = t132 * t193;
t75 = t126 * t110 + t128 * t179;
t72 = qJD(2) * pkin(8) + t75;
t244 = t134 * t115 - t131 * t72;
t86 = (t126 * t132 - t128 * t135) * t127;
t82 = qJD(2) * t86;
t77 = qJD(1) * t82;
t27 = qJD(4) * t244 - t134 * t77;
t200 = t131 * t115;
t56 = t134 * t72 + t200;
t51 = qJD(4) * pkin(9) + t56;
t147 = -t134 * pkin(4) - t131 * pkin(9) - pkin(3);
t111 = t126 * t179;
t74 = t128 * t110 - t111;
t61 = t147 * qJD(2) - t74;
t63 = (t109 + t80) * qJD(2);
t167 = -t130 * t63 - t133 * t27 - t61 * t187 + t51 * t188;
t17 = -t130 * t51 + t133 * t61;
t248 = t17 * t116 - t167;
t18 = t130 * t61 + t133 * t51;
t7 = -qJD(5) * t18 - t130 * t27 + t133 * t63;
t247 = t18 * t116 - t7;
t118 = t126 * pkin(2) + pkin(8);
t186 = t130 * qJD(4);
t201 = t130 * t134;
t83 = t128 * t172 - t111;
t246 = t131 * t118 * t186 + t249 * t133 + t83 * t201;
t237 = t128 * pkin(2);
t100 = t147 - t237;
t198 = t133 * t134;
t245 = t100 * t187 + t249 * t130 - t83 * t198;
t192 = qJD(2) * t131;
t171 = t130 * t192;
t185 = t133 * qJD(4);
t103 = t171 - t185;
t209 = t103 * t116;
t177 = t134 * t185;
t181 = qJD(4) * qJD(5);
t78 = -qJD(2) * t177 + qJD(5) * t171 - t133 * t181;
t243 = -t78 + t209;
t105 = t133 * t192 + t186;
t206 = t105 * t116;
t175 = t131 * t187;
t141 = t134 * t186 + t175;
t79 = t141 * qJD(2) + t130 * t181;
t242 = t79 - t206;
t124 = t131 ^ 2;
t148 = qJD(2) * t124 - t116 * t134;
t176 = t131 * t188;
t241 = -t116 * t176 - t148 * t185;
t240 = t105 ^ 2;
t239 = pkin(5) * t130;
t238 = t103 * pkin(5);
t151 = t129 * t134 - t87 * t131;
t28 = t56 * qJD(4) - t131 * t77;
t236 = t28 * t151;
t81 = qJD(2) * t87;
t76 = qJD(1) * t81;
t235 = t76 * t86;
t234 = -qJ(6) - pkin(9);
t14 = -t105 * qJ(6) + t17;
t13 = -t116 * pkin(5) + t14;
t233 = t13 - t14;
t107 = t118 * t198;
t146 = pkin(5) * t131 - qJ(6) * t198;
t184 = t133 * qJD(6);
t232 = -t131 * t184 + t146 * qJD(4) + (-t107 + (qJ(6) * t131 - t100) * t130) * qJD(5) + t246;
t199 = t131 * t133;
t231 = (-qJ(6) * qJD(5) - qJD(4) * t118) * t199 + (-qJD(6) * t131 + (-qJ(6) * qJD(4) - qJD(5) * t118) * t134) * t130 + t245;
t168 = qJD(5) * t234;
t108 = t158 * qJD(2);
t34 = t133 * t108 - t130 * t244;
t230 = t146 * qJD(2) + t130 * qJD(6) - t133 * t168 + t34;
t35 = t130 * t108 + t133 * t244;
t229 = -t184 + t35 + (-qJ(6) * t183 - t168) * t130;
t228 = (-t131 * t185 - t134 * t188) * t118 + t245;
t67 = t130 * t100 + t107;
t227 = -t67 * qJD(5) + t246;
t226 = -t103 * t177 - t79 * t199;
t224 = t103 * t83;
t223 = t105 * t83;
t50 = -qJD(4) * pkin(4) - t244;
t222 = t130 * t50;
t221 = t131 * t83;
t220 = t133 * t50;
t219 = t134 * t79;
t15 = -t103 * qJ(6) + t18;
t218 = t15 * t116;
t16 = t79 * pkin(5) + t28;
t217 = t16 * t130;
t216 = t16 * t133;
t213 = t28 * t130;
t212 = t28 * t131;
t211 = t28 * t133;
t208 = t103 * t131;
t207 = t105 * t103;
t205 = t116 * t130;
t204 = t116 * t133;
t137 = qJD(2) ^ 2;
t203 = t127 * t137;
t202 = t130 * t131;
t136 = qJD(4) ^ 2;
t196 = t136 * t131;
t195 = t136 * t134;
t125 = t134 ^ 2;
t194 = t124 - t125;
t191 = qJD(4) * t131;
t190 = qJD(4) * t134;
t189 = qJD(5) * t103;
t182 = qJD(2) * qJD(4);
t180 = t131 * t137 * t134;
t178 = t105 * t190;
t174 = t116 * t187;
t173 = t116 * t192;
t120 = t131 * t182;
t169 = t105 * t191 + t78 * t134;
t166 = -qJD(6) - t238;
t71 = -qJD(2) * pkin(3) - t74;
t165 = -qJD(2) * t71 + t77;
t164 = -t78 + t189;
t163 = pkin(5) * t120;
t161 = t131 * t174;
t160 = t105 * t175;
t159 = t134 * t120;
t157 = -t13 * t133 - t130 * t15;
t156 = t13 * t130 - t133 * t15;
t155 = -t130 * t18 - t133 * t17;
t154 = t130 * t17 - t133 * t18;
t70 = t129 * t131 + t87 * t134;
t42 = t86 * t130 + t70 * t133;
t41 = -t70 * t130 + t86 * t133;
t153 = t131 * t244 - t134 * t56;
t145 = t79 * qJ(6) + t167;
t144 = t80 * qJD(2) - t118 * t136 - t76;
t119 = -pkin(3) - t237;
t143 = qJD(4) * (qJD(2) * t119 + t71 + t83);
t142 = t148 * t130;
t140 = t155 * qJD(5) - t7 * t130 - t133 * t167;
t139 = t212 + t27 * t134 + (-t131 * t56 - t134 * t244) * qJD(4);
t138 = t78 * qJ(6) + t7;
t122 = -t133 * pkin(5) - pkin(4);
t113 = t234 * t133;
t112 = t234 * t130;
t102 = t103 ^ 2;
t94 = (t118 + t239) * t131;
t90 = t133 * t100;
t84 = (-t116 - t183) * t191;
t73 = t141 * pkin(5) + t118 * t190;
t66 = -t118 * t201 + t90;
t64 = -qJ(6) * t202 + t67;
t62 = -t102 + t240;
t60 = -qJ(6) * t199 + t90 + (-t118 * t130 - pkin(5)) * t134;
t54 = -t206 - t79;
t53 = -t78 - t209;
t47 = t200 + (qJD(2) * t239 + t72) * t134;
t46 = -t174 + (t116 * t198 + (-t105 + t186) * t131) * qJD(2);
t45 = t116 * t188 + (-t116 * t201 + (t103 + t185) * t131) * qJD(2);
t40 = t151 * qJD(4) - t82 * t134;
t39 = t70 * qJD(4) - t82 * t131;
t36 = -t166 + t50;
t33 = -t103 * t205 - t79 * t133;
t32 = -t105 * t204 - t78 * t130;
t31 = t141 * t103 + t79 * t202;
t30 = -t78 * t199 + (-t176 + t177) * t105;
t24 = t161 + t219 + (-t142 - t208) * qJD(4);
t23 = t161 - t219 + (-t142 + t208) * qJD(4);
t22 = t169 - t241;
t21 = t169 + t241;
t12 = -t130 * t242 + t133 * t243;
t11 = -t160 + (-t178 + (t78 + t189) * t131) * t130 + t226;
t10 = t160 + (t164 * t131 + t178) * t130 + t226;
t9 = t41 * qJD(5) + t81 * t130 + t40 * t133;
t8 = -t42 * qJD(5) - t40 * t130 + t81 * t133;
t5 = -t103 * qJD(6) - t145;
t4 = -t105 * qJD(6) + t138 + t163;
t3 = t39 * t105 + t9 * t116 - t42 * t120 + t151 * t78;
t2 = t39 * t103 - t8 * t116 + t41 * t120 - t151 * t79;
t1 = -t9 * t103 - t8 * t105 + t41 * t78 - t42 * t79;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132 * t203, -t135 * t203, 0, 0, 0, 0, 0, 0, 0, 0, -t81 * qJD(2), t82 * qJD(2), 0, -t74 * t81 - t75 * t82 - t77 * t87 + t235, 0, 0, 0, 0, 0, 0, -t39 * qJD(4) + (-t134 * t81 + t86 * t191) * qJD(2), -t40 * qJD(4) + (t131 * t81 + t86 * t190) * qJD(2) (t131 * t39 + t134 * t40 + (-t131 * t70 - t134 * t151) * qJD(4)) * qJD(2), -t244 * t39 + t27 * t70 + t56 * t40 + t71 * t81 + t235 - t236, 0, 0, 0, 0, 0, 0, t2, t3, t1, -t167 * t42 + t17 * t8 + t18 * t9 + t50 * t39 + t7 * t41 - t236, 0, 0, 0, 0, 0, 0, t2, t3, t1, t13 * t8 + t15 * t9 - t151 * t16 + t36 * t39 + t4 * t41 + t5 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t150 * t193 + t80) * qJD(2) (qJD(1) * t86 + t83) * qJD(2), 0, t74 * t80 - t75 * t83 + (-t126 * t77 - t128 * t76) * pkin(2), 0.2e1 * t159, -0.2e1 * t194 * t182, t195, -0.2e1 * t159, -t196, 0, t131 * t143 + t144 * t134, -t144 * t131 + t134 * t143 (-t124 - t125) * t83 * qJD(2) + t139, t118 * t139 + t76 * t119 + t153 * t83 - t71 * t80, t30, t11, t22, t31, t24, t84, -t227 * t116 + (-t7 + (t103 * t118 + t222) * qJD(4)) * t134 + (t50 * t187 - t224 + t118 * t79 + t213 + (qJD(2) * t66 + t17) * qJD(4)) * t131, t228 * t116 + (-t167 + (t105 * t118 + t220) * qJD(4)) * t134 + (-t50 * t188 - t223 - t118 * t78 + t211 + (-qJD(2) * t67 - t18) * qJD(4)) * t131, t66 * t78 - t67 * t79 - t227 * t105 - t228 * t103 + t155 * t190 + (qJD(5) * t154 + t130 * t167 - t133 * t7) * t131, -t50 * t221 - t167 * t67 + t7 * t66 + t228 * t18 + t227 * t17 + (t190 * t50 + t212) * t118, t30, t11, t22, t31, t24, t84, t73 * t103 + t94 * t79 + (t186 * t36 - t4) * t134 - t232 * t116 + (t36 * t187 - t224 + t217 + (qJD(2) * t60 + t13) * qJD(4)) * t131, t73 * t105 - t94 * t78 + (t185 * t36 + t5) * t134 + t231 * t116 + (-t36 * t188 - t223 + t216 + (-qJD(2) * t64 - t15) * qJD(4)) * t131, t60 * t78 - t64 * t79 - t232 * t105 - t231 * t103 + t157 * t190 + (qJD(5) * t156 - t130 * t5 - t133 * t4) * t131, t16 * t94 + t4 * t60 + t5 * t64 + (t73 - t221) * t36 + t231 * t15 + t232 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, -t195, 0, -qJD(4) * t153 + t27 * t131 - t28 * t134, 0, 0, 0, 0, 0, 0, t23, t21, t10 (-qJD(4) * t154 - t28) * t134 + (qJD(4) * t50 + t140) * t131, 0, 0, 0, 0, 0, 0, t23, t21, t10 (-qJD(4) * t156 - t16) * t134 + (qJD(4) * t36 + qJD(5) * t157 - t4 * t130 + t5 * t133) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t194 * t137, 0, t180, 0, 0, t165 * t131, t165 * t134, 0, 0, t32, t12, t46, t33, t45, t173, -pkin(4) * t79 - t56 * t103 + t34 * t116 - t211 + (pkin(9) * t204 + t222) * qJD(5) + (-t131 * t17 + (-pkin(9) * t191 - t134 * t50) * t130) * qJD(2), pkin(4) * t78 - t56 * t105 - t35 * t116 + t213 + (-pkin(9) * t205 + t220) * qJD(5) + (-t50 * t198 + (-pkin(9) * t185 + t18) * t131) * qJD(2), t35 * t103 + t34 * t105 + ((qJD(5) * t105 - t79) * pkin(9) + t248) * t133 + (pkin(9) * t164 + t247) * t130, -t28 * pkin(4) + pkin(9) * t140 - t17 * t34 - t18 * t35 - t50 * t56, t32, t12, t46, t33, t45, t173, -t47 * t103 + t122 * t79 - t216 + t230 * t116 + (t36 + t238) * t188 + (-t36 * t201 + (qJD(4) * t112 - t13) * t131) * qJD(2), -t47 * t105 - t122 * t78 + t217 - t229 * t116 + (t105 * t239 + t133 * t36) * qJD(5) + (-t36 * t198 + (qJD(4) * t113 + t15) * t131) * qJD(2), t112 * t78 + t113 * t79 + t230 * t105 + t229 * t103 + (t116 * t13 + t5) * t133 + (-t4 + t218) * t130, t4 * t112 - t5 * t113 + t16 * t122 + (pkin(5) * t188 - t47) * t36 - t229 * t15 - t230 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t62, t53, -t207, t54, t120, -t50 * t105 - t247, t50 * t103 - t248, 0, 0, t207, t62, t53, -t207, t54, t120, 0.2e1 * t163 - t218 + (t166 - t36) * t105 + t138, -t240 * pkin(5) - t14 * t116 + (qJD(6) + t36) * t103 + t145, t78 * pkin(5) - t233 * t103, t233 * t15 + (-t36 * t105 + t4) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t243, -t102 - t240, t15 * t103 + t13 * t105 + t16;];
tauc_reg  = t6;
