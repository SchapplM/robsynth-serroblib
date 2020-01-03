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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:39:27
% EndTime: 2020-01-03 11:39:33
% DurationCPUTime: 2.20s
% Computational Cost: add. (3456->315), mult. (7570->396), div. (0->0), fcn. (5361->16), ass. (0->177)
t163 = sin(qJ(5));
t158 = sin(pkin(9));
t160 = cos(pkin(9));
t166 = cos(qJ(3));
t218 = t160 * t166;
t200 = qJD(1) * t218;
t164 = sin(qJ(3));
t211 = qJD(1) * t164;
t102 = t158 * t211 - t200;
t237 = cos(qJ(5));
t110 = t158 * t166 + t160 * t164;
t239 = t110 * qJD(1);
t183 = t163 * t102 - t237 * t239;
t203 = t166 * qJDD(1);
t204 = t164 * qJDD(1);
t186 = t158 * t204 - t160 * t203;
t246 = t239 * qJD(3);
t57 = t186 + t246;
t207 = qJD(1) * qJD(3);
t198 = t164 * t207;
t178 = t110 * qJDD(1) - t158 * t198;
t197 = t166 * t207;
t58 = t160 * t197 + t178;
t174 = t183 * qJD(5) - t163 * t58 - t237 * t57;
t153 = qJD(3) + qJD(5);
t222 = t183 * t153;
t249 = t174 - t222;
t199 = qJD(5) * t237;
t209 = qJD(5) * t163;
t182 = -t102 * t199 - t163 * t57 - t209 * t239 + t237 * t58;
t48 = -t237 * t102 - t163 * t239;
t221 = t48 * t153;
t248 = t182 - t221;
t229 = t183 ^ 2;
t230 = t48 ^ 2;
t247 = t229 - t230;
t228 = t48 * t183;
t159 = sin(pkin(8));
t135 = t159 * pkin(1) + pkin(6);
t119 = t135 * qJDD(1);
t192 = -qJD(2) * qJD(3) - t119;
t155 = qJ(1) + pkin(8);
t141 = sin(t155);
t143 = cos(t155);
t241 = g(2) * t141 - g(3) * t143;
t232 = t239 * pkin(7);
t208 = t164 * qJD(2);
t121 = t135 * qJD(1);
t193 = qJ(4) * qJD(1) + t121;
t243 = t193 * t166;
t78 = t208 + t243;
t67 = t158 * t78;
t147 = t166 * qJD(2);
t77 = -t193 * t164 + t147;
t71 = qJD(3) * pkin(3) + t77;
t35 = t160 * t71 - t67;
t22 = qJD(3) * pkin(4) - t232 + t35;
t233 = t102 * pkin(7);
t223 = t160 * t78;
t36 = t158 * t71 + t223;
t23 = t36 - t233;
t145 = t166 * qJDD(2);
t206 = qJD(1) * qJD(4);
t33 = qJDD(3) * pkin(3) + t145 - qJD(3) * t243 + (-qJ(4) * qJDD(1) + t192 - t206) * t164;
t201 = -t164 * qJDD(2) + t192 * t166;
t210 = qJD(3) * t164;
t50 = -t121 * t210 - t201;
t37 = t166 * t206 + (-t198 + t203) * qJ(4) + t50;
t13 = -t158 * t37 + t160 * t33;
t6 = qJDD(3) * pkin(4) - t58 * pkin(7) + t13;
t14 = t158 * t33 + t160 * t37;
t7 = -t57 * pkin(7) + t14;
t1 = t163 * t6 + t22 * t199 - t23 * t209 + t237 * t7;
t154 = qJ(3) + pkin(9);
t146 = qJ(5) + t154;
t133 = sin(t146);
t134 = cos(t146);
t161 = cos(pkin(8));
t137 = -t161 * pkin(1) - pkin(2);
t149 = t166 * pkin(3);
t240 = t137 - t149;
t100 = qJD(1) * t240 + qJD(4);
t61 = t102 * pkin(4) + t100;
t245 = g(1) * t133 + t241 * t134 - t61 * t48 - t1;
t9 = t163 * t22 + t237 * t23;
t2 = -t9 * qJD(5) - t163 * t7 + t237 * t6;
t244 = -g(1) * t134 + t241 * t133 + t183 * t61 + t2;
t220 = pkin(1) * qJDD(1);
t242 = g(2) * t143 + g(3) * t141;
t238 = t239 ^ 2;
t168 = qJD(3) ^ 2;
t236 = pkin(3) * t158;
t235 = g(1) * t166;
t231 = t164 * pkin(3);
t162 = -qJ(4) - pkin(6);
t109 = t158 * t164 - t218;
t106 = t109 * qJD(3);
t181 = t110 * qJD(3);
t24 = t237 * t106 + t109 * t199 + t110 * t209 + t163 * t181;
t63 = -t163 * t109 + t237 * t110;
t227 = t174 * t63 - t24 * t48;
t226 = t106 * t102 - t110 * t57;
t39 = t160 * t77 - t67;
t215 = qJ(4) + t135;
t194 = qJD(3) * t215;
t81 = t166 * qJD(4) - t164 * t194;
t82 = -t164 * qJD(4) - t166 * t194;
t41 = t158 * t82 + t160 * t81;
t38 = -t158 * t77 - t223;
t26 = t38 + t233;
t27 = t39 - t232;
t136 = t160 * pkin(3) + pkin(4);
t97 = t237 * t136 - t163 * t236;
t225 = t97 * qJD(5) - t163 * t26 - t237 * t27;
t98 = t163 * t136 + t237 * t236;
t224 = -t98 * qJD(5) + t163 * t27 - t237 * t26;
t107 = t215 * t164;
t108 = t215 * t166;
t56 = -t158 * t107 + t160 * t108;
t219 = t239 * t102;
t217 = t164 * t121;
t216 = t166 * t121;
t142 = cos(t154);
t213 = pkin(4) * t142 + t149;
t156 = t164 ^ 2;
t157 = t166 ^ 2;
t212 = t156 - t157;
t122 = qJD(1) * t137;
t120 = qJDD(1) * t137;
t169 = qJD(1) ^ 2;
t202 = t164 * t169 * t166;
t40 = -t158 * t81 + t160 * t82;
t55 = -t160 * t107 - t158 * t108;
t191 = t164 * t197;
t165 = sin(qJ(1));
t167 = cos(qJ(1));
t188 = -g(2) * t167 - g(3) * t165;
t25 = t63 * qJD(5) - t163 * t106 + t237 * t181;
t62 = t237 * t109 + t163 * t110;
t187 = t182 * t62 - t183 * t25;
t151 = qJDD(3) + qJDD(5);
t185 = t63 * t151 - t24 * t153;
t92 = t208 + t216;
t42 = -t110 * pkin(7) + t55;
t43 = -t109 * pkin(7) + t56;
t19 = -t163 * t43 + t237 * t42;
t20 = t163 * t42 + t237 * t43;
t179 = -qJD(1) * t122 + t241;
t177 = 0.2e1 * t122 * qJD(3) - qJDD(3) * t135;
t76 = pkin(3) * t198 + qJDD(1) * t240 + qJDD(4);
t176 = -t58 * t109 - t181 * t239;
t175 = t135 * t168 + 0.2e1 * t120 + t242;
t51 = -t92 * qJD(3) - t164 * t119 + t145;
t91 = t147 - t217;
t173 = -t51 * t164 + t50 * t166 + (-t164 * t92 - t166 * t91) * qJD(3);
t34 = t57 * pkin(4) + t76;
t152 = -pkin(7) + t162;
t150 = t167 * pkin(1);
t148 = t165 * pkin(1);
t140 = sin(t154);
t139 = t149 + pkin(2);
t118 = qJDD(3) * t166 - t168 * t164;
t117 = qJDD(3) * t164 + t168 * t166;
t114 = pkin(2) + t213;
t99 = t102 ^ 2;
t80 = (t110 * pkin(4) + t231) * qJD(3);
t79 = pkin(3) * t211 + pkin(4) * t239;
t75 = t109 * pkin(4) + t240;
t60 = -t106 * qJD(3) + t110 * qJDD(3);
t59 = -t109 * qJDD(3) - t110 * t168;
t29 = -pkin(7) * t181 + t41;
t28 = t106 * pkin(7) + t40;
t15 = -t62 * t151 - t25 * t153;
t8 = -t163 * t23 + t237 * t22;
t4 = -t20 * qJD(5) - t163 * t29 + t237 * t28;
t3 = t19 * qJD(5) + t163 * t28 + t237 * t29;
t5 = [0, 0, 0, 0, 0, qJDD(1), t188, g(2) * t165 - g(3) * t167, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t161 * t220 - t242, -0.2e1 * t159 * t220 + t241, 0, (t188 + (t159 ^ 2 + t161 ^ 2) * t220) * pkin(1), t156 * qJDD(1) + 0.2e1 * t191, 0.2e1 * t164 * t203 - 0.2e1 * t212 * t207, t117, t157 * qJDD(1) - 0.2e1 * t191, t118, 0, t164 * t177 - t166 * t175, t164 * t175 + t166 * t177, (t156 + t157) * t119 + t173 - t241, t120 * t137 - g(2) * (t143 * pkin(2) + t141 * pkin(6) + t150) - g(3) * (t141 * pkin(2) - t143 * pkin(6) + t148) + t173 * t135, -t106 * t239 + t58 * t110, t176 + t226, t60, t102 * t181 + t57 * t109, t59, 0, t55 * qJDD(3) + t76 * t109 + t240 * t57 - t242 * t142 + (t100 * t110 + t102 * t231 + t40) * qJD(3), -t56 * qJDD(3) - t100 * t106 + t76 * t110 + t240 * t58 + t242 * t140 + (t231 * t239 - t41) * qJD(3), -t41 * t102 + t35 * t106 - t14 * t109 - t13 * t110 - t181 * t36 - t239 * t40 - t55 * t58 - t56 * t57 - t241, t14 * t56 + t36 * t41 + t13 * t55 + t35 * t40 + t76 * t240 + t100 * pkin(3) * t210 - g(2) * (t143 * t139 - t141 * t162 + t150) - g(3) * (t141 * t139 + t143 * t162 + t148), t182 * t63 + t183 * t24, -t187 + t227, t185, -t174 * t62 - t25 * t48, t15, 0, -t134 * t242 + t19 * t151 + t4 * t153 - t174 * t75 + t61 * t25 + t34 * t62 - t48 * t80, t133 * t242 - t20 * t151 - t3 * t153 + t182 * t75 - t183 * t80 - t61 * t24 + t34 * t63, -t1 * t62 + t174 * t20 - t182 * t19 + t183 * t4 - t2 * t63 + t8 * t24 - t9 * t25 + t3 * t48 - t241, t1 * t20 + t9 * t3 + t2 * t19 + t8 * t4 + t34 * t75 + t61 * t80 - g(2) * (t143 * t114 - t141 * t152 + t150) - g(3) * (t141 * t114 + t143 * t152 + t148); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, t118, -t117, 0, t50 * t164 + t51 * t166 - g(1) + (-t164 * t91 + t166 * t92) * qJD(3), 0, 0, 0, 0, 0, 0, t59, -t60, -t176 + t226, -t36 * t106 - t13 * t109 + t14 * t110 - t181 * t35 - g(1), 0, 0, 0, 0, 0, 0, t15, -t185, t187 + t227, t1 * t63 - t2 * t62 - t9 * t24 - t8 * t25 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t212 * t169, t204, t202, t203, qJDD(3), -t235 + t145 + (t92 - t216) * qJD(3) + (t179 + t192) * t164, g(1) * t164 + (t91 + t217) * qJD(3) + t179 * t166 + t201, 0, 0, t219, -t99 + t238, (t102 + t200) * qJD(3) + t178, -t219, -t186, qJDD(3), -g(1) * t142 - t38 * qJD(3) - t100 * t239 + t241 * t140 + (qJDD(3) * t160 - t102 * t211) * pkin(3) + t13, g(1) * t140 + t39 * qJD(3) + t100 * t102 + t241 * t142 + (-qJDD(3) * t158 - t211 * t239) * pkin(3) - t14, (t36 + t38) * t239 + (-t35 + t39) * t102 + (-t158 * t57 - t160 * t58) * pkin(3), -t35 * t38 - t36 * t39 + (-t235 + t13 * t160 + t14 * t158 + (-qJD(1) * t100 + t241) * t164) * pkin(3), t228, t247, t248, -t228, t249, t151, t97 * t151 + t153 * t224 + t48 * t79 + t244, -t98 * t151 - t153 * t225 + t183 * t79 + t245, t98 * t174 - t182 * t97 + (t225 + t8) * t48 + (t224 - t9) * t183, t1 * t98 + t2 * t97 - t61 * t79 - g(1) * t213 + t225 * t9 + t224 * t8 - t241 * (-pkin(4) * t140 - t231); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186 + 0.2e1 * t246, (-t102 + t200) * qJD(3) + t178, -t99 - t238, t36 * t102 + t239 * t35 + t242 + t76, 0, 0, 0, 0, 0, 0, -t174 - t222, t182 + t221, -t229 - t230, -t183 * t8 - t9 * t48 + t242 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t247, t248, -t228, t249, t151, t9 * t153 + t244, t8 * t153 + t245, 0, 0;];
tau_reg = t5;
