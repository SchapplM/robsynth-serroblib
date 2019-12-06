% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:54:05
% EndTime: 2019-12-05 17:54:10
% DurationCPUTime: 2.28s
% Computational Cost: add. (3456->315), mult. (7570->396), div. (0->0), fcn. (5361->16), ass. (0->179)
t159 = sin(qJ(5));
t154 = sin(pkin(9));
t156 = cos(pkin(9));
t162 = cos(qJ(3));
t214 = t156 * t162;
t197 = qJD(1) * t214;
t160 = sin(qJ(3));
t208 = qJD(1) * t160;
t102 = t154 * t208 - t197;
t235 = cos(qJ(5));
t110 = t154 * t162 + t156 * t160;
t238 = t110 * qJD(1);
t180 = t159 * t102 - t235 * t238;
t200 = t162 * qJDD(1);
t201 = t160 * qJDD(1);
t183 = t154 * t201 - t156 * t200;
t243 = t238 * qJD(3);
t57 = t183 + t243;
t204 = qJD(1) * qJD(3);
t195 = t160 * t204;
t175 = qJDD(1) * t110 - t154 * t195;
t194 = t162 * t204;
t58 = t156 * t194 + t175;
t171 = qJD(5) * t180 - t159 * t58 - t235 * t57;
t149 = qJD(3) + qJD(5);
t218 = t180 * t149;
t246 = t171 - t218;
t196 = qJD(5) * t235;
t206 = qJD(5) * t159;
t179 = -t102 * t196 - t159 * t57 - t206 * t238 + t235 * t58;
t48 = -t235 * t102 - t159 * t238;
t217 = t48 * t149;
t245 = t179 - t217;
t225 = t180 ^ 2;
t226 = t48 ^ 2;
t244 = t225 - t226;
t224 = t48 * t180;
t155 = sin(pkin(8));
t133 = t155 * pkin(1) + pkin(6);
t119 = t133 * qJDD(1);
t189 = -qJD(2) * qJD(3) - t119;
t151 = qJ(1) + pkin(8);
t139 = sin(t151);
t141 = cos(t151);
t186 = g(2) * t139 - g(3) * t141;
t230 = t238 * pkin(7);
t205 = t160 * qJD(2);
t121 = t133 * qJD(1);
t190 = qJ(4) * qJD(1) + t121;
t240 = t162 * t190;
t78 = t205 + t240;
t67 = t154 * t78;
t145 = t162 * qJD(2);
t77 = -t190 * t160 + t145;
t71 = qJD(3) * pkin(3) + t77;
t35 = t156 * t71 - t67;
t22 = qJD(3) * pkin(4) - t230 + t35;
t231 = t102 * pkin(7);
t219 = t156 * t78;
t36 = t154 * t71 + t219;
t23 = t36 - t231;
t143 = t162 * qJDD(2);
t203 = qJD(1) * qJD(4);
t33 = qJDD(3) * pkin(3) + t143 - qJD(3) * t240 + (-qJ(4) * qJDD(1) + t189 - t203) * t160;
t198 = -t160 * qJDD(2) + t189 * t162;
t207 = qJD(3) * t160;
t50 = -t121 * t207 - t198;
t37 = t162 * t203 + (-t195 + t200) * qJ(4) + t50;
t13 = -t154 * t37 + t156 * t33;
t6 = qJDD(3) * pkin(4) - t58 * pkin(7) + t13;
t14 = t154 * t33 + t156 * t37;
t7 = -t57 * pkin(7) + t14;
t1 = t159 * t6 + t22 * t196 - t23 * t206 + t235 * t7;
t150 = qJ(3) + pkin(9);
t144 = qJ(5) + t150;
t131 = sin(t144);
t132 = cos(t144);
t157 = cos(pkin(8));
t135 = -t157 * pkin(1) - pkin(2);
t146 = t162 * pkin(3);
t239 = t135 - t146;
t100 = qJD(1) * t239 + qJD(4);
t61 = t102 * pkin(4) + t100;
t242 = g(1) * t131 - t186 * t132 - t61 * t48 - t1;
t9 = t159 * t22 + t235 * t23;
t2 = -qJD(5) * t9 - t159 * t7 + t235 * t6;
t241 = -g(1) * t132 - t186 * t131 + t180 * t61 + t2;
t216 = pkin(1) * qJDD(1);
t237 = t238 ^ 2;
t164 = qJD(3) ^ 2;
t236 = t57 * pkin(4);
t234 = pkin(3) * t154;
t233 = g(1) * t162;
t229 = t160 * pkin(3);
t161 = sin(qJ(1));
t228 = t161 * pkin(1);
t163 = cos(qJ(1));
t227 = t163 * pkin(1);
t158 = -qJ(4) - pkin(6);
t109 = t154 * t160 - t214;
t106 = t109 * qJD(3);
t178 = t110 * qJD(3);
t24 = t235 * t106 + t109 * t196 + t110 * t206 + t159 * t178;
t63 = -t159 * t109 + t235 * t110;
t223 = t171 * t63 - t24 * t48;
t222 = t106 * t102 - t110 * t57;
t39 = t156 * t77 - t67;
t211 = qJ(4) + t133;
t191 = qJD(3) * t211;
t81 = t162 * qJD(4) - t160 * t191;
t82 = -t160 * qJD(4) - t162 * t191;
t41 = t154 * t82 + t156 * t81;
t38 = -t154 * t77 - t219;
t26 = t38 + t231;
t27 = t39 - t230;
t134 = t156 * pkin(3) + pkin(4);
t97 = t235 * t134 - t159 * t234;
t221 = t97 * qJD(5) - t159 * t26 - t235 * t27;
t98 = t159 * t134 + t235 * t234;
t220 = -t98 * qJD(5) + t159 * t27 - t235 * t26;
t107 = t211 * t160;
t108 = t211 * t162;
t56 = -t154 * t107 + t156 * t108;
t215 = t238 * t102;
t213 = t160 * t121;
t212 = t162 * t121;
t140 = cos(t150);
t210 = pkin(4) * t140 + t146;
t152 = t160 ^ 2;
t153 = t162 ^ 2;
t209 = t152 - t153;
t122 = qJD(1) * t135;
t120 = qJDD(1) * t135;
t165 = qJD(1) ^ 2;
t199 = t160 * t165 * t162;
t40 = -t154 * t81 + t156 * t82;
t55 = -t156 * t107 - t154 * t108;
t188 = t160 * t194;
t187 = g(2) * t141 + g(3) * t139;
t185 = g(2) * t163 + g(3) * t161;
t25 = t63 * qJD(5) - t159 * t106 + t235 * t178;
t62 = t235 * t109 + t159 * t110;
t184 = t179 * t62 - t180 * t25;
t147 = qJDD(3) + qJDD(5);
t182 = t63 * t147 - t24 * t149;
t92 = t205 + t212;
t42 = -t110 * pkin(7) + t55;
t43 = -t109 * pkin(7) + t56;
t19 = -t159 * t43 + t235 * t42;
t20 = t159 * t42 + t235 * t43;
t176 = -qJD(1) * t122 - t186;
t174 = 0.2e1 * t122 * qJD(3) - qJDD(3) * t133;
t76 = pkin(3) * t195 + qJDD(1) * t239 + qJDD(4);
t173 = -t58 * t109 - t178 * t238;
t172 = -t133 * t164 - 0.2e1 * t120 + t187;
t51 = -t92 * qJD(3) - t160 * t119 + t143;
t91 = t145 - t213;
t170 = -t51 * t160 + t50 * t162 + (-t160 * t92 - t162 * t91) * qJD(3);
t169 = -t187 + t76;
t148 = -pkin(7) + t158;
t138 = sin(t150);
t137 = t146 + pkin(2);
t118 = qJDD(3) * t162 - t164 * t160;
t117 = qJDD(3) * t160 + t164 * t162;
t114 = pkin(2) + t210;
t99 = t102 ^ 2;
t80 = (pkin(4) * t110 + t229) * qJD(3);
t79 = pkin(3) * t208 + pkin(4) * t238;
t75 = t109 * pkin(4) + t239;
t60 = -t106 * qJD(3) + t110 * qJDD(3);
t59 = -t109 * qJDD(3) - t164 * t110;
t34 = t76 + t236;
t29 = -pkin(7) * t178 + t41;
t28 = t106 * pkin(7) + t40;
t15 = -t62 * t147 - t25 * t149;
t8 = -t159 * t23 + t235 * t22;
t4 = -t20 * qJD(5) - t159 * t29 + t235 * t28;
t3 = t19 * qJD(5) + t159 * t28 + t235 * t29;
t5 = [0, 0, 0, 0, 0, qJDD(1), t185, -g(2) * t161 + g(3) * t163, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t157 * t216 + t187, -0.2e1 * t155 * t216 - t186, 0, (t185 + (t155 ^ 2 + t157 ^ 2) * t216) * pkin(1), t152 * qJDD(1) + 0.2e1 * t188, 0.2e1 * t160 * t200 - 0.2e1 * t209 * t204, t117, t153 * qJDD(1) - 0.2e1 * t188, t118, 0, t160 * t174 + t162 * t172, -t160 * t172 + t162 * t174, (t152 + t153) * t119 + t170 + t186, t120 * t135 - g(2) * (-t141 * pkin(2) - t139 * pkin(6) - t227) - g(3) * (-t139 * pkin(2) + t141 * pkin(6) - t228) + t170 * t133, -t106 * t238 + t58 * t110, t173 + t222, t60, t102 * t178 + t57 * t109, t59, 0, t55 * qJDD(3) + t76 * t109 + t239 * t57 + t187 * t140 + (t100 * t110 + t102 * t229 + t40) * qJD(3), -t56 * qJDD(3) - t100 * t106 + t76 * t110 + t239 * t58 - t187 * t138 + (t229 * t238 - t41) * qJD(3), -t41 * t102 + t35 * t106 - t14 * t109 - t13 * t110 - t178 * t36 - t238 * t40 - t55 * t58 - t56 * t57 + t186, t14 * t56 + t36 * t41 + t13 * t55 + t35 * t40 + t76 * t239 + t100 * pkin(3) * t207 - g(2) * (-t141 * t137 + t139 * t158 - t227) - g(3) * (-t139 * t137 - t141 * t158 - t228), t179 * t63 + t180 * t24, -t184 + t223, t182, -t171 * t62 - t25 * t48, t15, 0, t132 * t187 + t19 * t147 + t4 * t149 - t171 * t75 + t61 * t25 + t34 * t62 - t48 * t80, -t131 * t187 - t20 * t147 - t3 * t149 + t179 * t75 - t180 * t80 - t61 * t24 + t34 * t63, -t1 * t62 + t171 * t20 - t179 * t19 + t180 * t4 - t2 * t63 + t8 * t24 - t9 * t25 + t3 * t48 + t186, t1 * t20 + t9 * t3 + t2 * t19 + t8 * t4 + t34 * t75 + t61 * t80 - g(2) * (-t141 * t114 + t139 * t148 - t227) - g(3) * (-t139 * t114 - t141 * t148 - t228); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, t118, -t117, 0, t50 * t160 + t51 * t162 - g(1) + (-t160 * t91 + t162 * t92) * qJD(3), 0, 0, 0, 0, 0, 0, t59, -t60, -t173 + t222, -t36 * t106 - t13 * t109 + t14 * t110 - t178 * t35 - g(1), 0, 0, 0, 0, 0, 0, t15, -t182, t184 + t223, t1 * t63 - t2 * t62 - t9 * t24 - t8 * t25 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, t209 * t165, t201, t199, t200, qJDD(3), -t233 + t143 + (t92 - t212) * qJD(3) + (t176 + t189) * t160, g(1) * t160 + (t91 + t213) * qJD(3) + t176 * t162 + t198, 0, 0, t215, -t99 + t237, (t102 + t197) * qJD(3) + t175, -t215, -t183, qJDD(3), -g(1) * t140 - t38 * qJD(3) - t100 * t238 - t186 * t138 + (qJDD(3) * t156 - t102 * t208) * pkin(3) + t13, g(1) * t138 + t39 * qJD(3) + t100 * t102 - t186 * t140 + (-qJDD(3) * t154 - t208 * t238) * pkin(3) - t14, (t36 + t38) * t238 + (-t35 + t39) * t102 + (-t154 * t57 - t156 * t58) * pkin(3), -t35 * t38 - t36 * t39 + (-t233 + t13 * t156 + t14 * t154 + (-qJD(1) * t100 - t186) * t160) * pkin(3), t224, t244, t245, -t224, t246, t147, t97 * t147 + t220 * t149 + t48 * t79 + t241, -t98 * t147 - t221 * t149 + t180 * t79 + t242, t98 * t171 - t179 * t97 + (t221 + t8) * t48 + (t220 - t9) * t180, t1 * t98 + t2 * t97 - t61 * t79 - g(1) * t210 + t221 * t9 + t220 * t8 + t186 * (-pkin(4) * t138 - t229); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183 + 0.2e1 * t243, (-t102 + t197) * qJD(3) + t175, -t99 - t237, t36 * t102 + t238 * t35 + t169, 0, 0, 0, 0, 0, 0, -t171 - t218, t179 + t217, -t225 - t226, -t180 * t8 - t9 * t48 + t169 + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t244, t245, -t224, t246, t147, t9 * t149 + t241, t8 * t149 + t242, 0, 0;];
tau_reg = t5;
