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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:23:10
% EndTime: 2022-01-23 09:23:16
% DurationCPUTime: 2.14s
% Computational Cost: add. (3456->315), mult. (7570->396), div. (0->0), fcn. (5361->16), ass. (0->179)
t162 = sin(qJ(5));
t157 = sin(pkin(9));
t159 = cos(pkin(9));
t165 = cos(qJ(3));
t217 = t159 * t165;
t200 = qJD(1) * t217;
t163 = sin(qJ(3));
t211 = qJD(1) * t163;
t102 = t157 * t211 - t200;
t239 = cos(qJ(5));
t110 = t157 * t165 + t159 * t163;
t242 = t110 * qJD(1);
t183 = t162 * t102 - t239 * t242;
t203 = t165 * qJDD(1);
t204 = t163 * qJDD(1);
t186 = t157 * t204 - t159 * t203;
t247 = t242 * qJD(3);
t57 = t186 + t247;
t207 = qJD(1) * qJD(3);
t198 = t163 * t207;
t178 = t110 * qJDD(1) - t157 * t198;
t197 = t165 * t207;
t58 = t159 * t197 + t178;
t174 = t183 * qJD(5) - t162 * t58 - t239 * t57;
t152 = qJD(3) + qJD(5);
t221 = t183 * t152;
t250 = t174 - t221;
t199 = qJD(5) * t239;
t209 = qJD(5) * t162;
t182 = -t102 * t199 - t162 * t57 - t209 * t242 + t239 * t58;
t48 = -t239 * t102 - t162 * t242;
t220 = t48 * t152;
t249 = t182 - t220;
t228 = t48 ^ 2;
t229 = t183 ^ 2;
t248 = -t228 + t229;
t227 = t48 * t183;
t158 = sin(pkin(8));
t135 = t158 * pkin(1) + pkin(6);
t121 = t135 * qJDD(1);
t192 = -qJD(2) * qJD(3) - t121;
t154 = qJ(1) + pkin(8);
t141 = sin(t154);
t143 = cos(t154);
t190 = g(1) * t143 + g(2) * t141;
t232 = t242 * pkin(7);
t208 = t163 * qJD(2);
t123 = t135 * qJD(1);
t193 = qJ(4) * qJD(1) + t123;
t244 = t193 * t165;
t78 = t208 + t244;
t67 = t157 * t78;
t147 = t165 * qJD(2);
t77 = -t193 * t163 + t147;
t71 = qJD(3) * pkin(3) + t77;
t35 = t159 * t71 - t67;
t22 = qJD(3) * pkin(4) - t232 + t35;
t233 = t102 * pkin(7);
t222 = t159 * t78;
t36 = t157 * t71 + t222;
t23 = t36 - t233;
t145 = t165 * qJDD(2);
t206 = qJD(1) * qJD(4);
t33 = qJDD(3) * pkin(3) + t145 - qJD(3) * t244 + (-qJ(4) * qJDD(1) + t192 - t206) * t163;
t201 = -t163 * qJDD(2) + t192 * t165;
t210 = qJD(3) * t163;
t50 = -t123 * t210 - t201;
t37 = t165 * t206 + (-t198 + t203) * qJ(4) + t50;
t13 = -t157 * t37 + t159 * t33;
t6 = qJDD(3) * pkin(4) - t58 * pkin(7) + t13;
t14 = t157 * t33 + t159 * t37;
t7 = -pkin(7) * t57 + t14;
t1 = t162 * t6 + t22 * t199 - t23 * t209 + t239 * t7;
t153 = qJ(3) + pkin(9);
t146 = qJ(5) + t153;
t133 = sin(t146);
t134 = cos(t146);
t160 = cos(pkin(8));
t137 = -t160 * pkin(1) - pkin(2);
t148 = t165 * pkin(3);
t243 = t137 - t148;
t100 = qJD(1) * t243 + qJD(4);
t61 = t102 * pkin(4) + t100;
t246 = g(3) * t133 + t190 * t134 - t61 * t48 - t1;
t9 = t162 * t22 + t239 * t23;
t2 = -t9 * qJD(5) - t162 * t7 + t239 * t6;
t245 = -g(3) * t134 + t190 * t133 + t183 * t61 + t2;
t219 = pkin(1) * qJDD(1);
t196 = -g(1) * t141 + g(2) * t143;
t241 = t242 ^ 2;
t167 = qJD(3) ^ 2;
t240 = t57 * pkin(4);
t238 = pkin(3) * t157;
t234 = g(3) * t165;
t231 = t163 * pkin(3);
t164 = sin(qJ(1));
t230 = t164 * pkin(1);
t161 = -qJ(4) - pkin(6);
t109 = t157 * t163 - t217;
t106 = t109 * qJD(3);
t181 = t110 * qJD(3);
t24 = t239 * t106 + t109 * t199 + t110 * t209 + t162 * t181;
t63 = -t162 * t109 + t239 * t110;
t226 = t174 * t63 - t24 * t48;
t225 = t106 * t102 - t110 * t57;
t39 = t159 * t77 - t67;
t214 = qJ(4) + t135;
t194 = qJD(3) * t214;
t81 = t165 * qJD(4) - t163 * t194;
t82 = -t163 * qJD(4) - t165 * t194;
t41 = t157 * t82 + t159 * t81;
t38 = -t157 * t77 - t222;
t26 = t38 + t233;
t27 = t39 - t232;
t136 = t159 * pkin(3) + pkin(4);
t97 = t239 * t136 - t162 * t238;
t224 = t97 * qJD(5) - t162 * t26 - t239 * t27;
t98 = t162 * t136 + t239 * t238;
t223 = -t98 * qJD(5) + t162 * t27 - t239 * t26;
t107 = t214 * t163;
t108 = t214 * t165;
t56 = -t157 * t107 + t159 * t108;
t218 = t242 * t102;
t216 = t163 * t123;
t215 = t165 * t123;
t142 = cos(t153);
t213 = pkin(4) * t142 + t148;
t155 = t163 ^ 2;
t156 = t165 ^ 2;
t212 = t155 - t156;
t124 = qJD(1) * t137;
t122 = qJDD(1) * t137;
t168 = qJD(1) ^ 2;
t202 = t163 * t168 * t165;
t40 = -t157 * t81 + t159 * t82;
t55 = -t159 * t107 - t157 * t108;
t191 = t163 * t197;
t166 = cos(qJ(1));
t188 = g(1) * t164 - g(2) * t166;
t25 = t63 * qJD(5) - t162 * t106 + t239 * t181;
t62 = t239 * t109 + t162 * t110;
t187 = t182 * t62 - t183 * t25;
t150 = qJDD(3) + qJDD(5);
t185 = t63 * t150 - t24 * t152;
t92 = t208 + t215;
t42 = -t110 * pkin(7) + t55;
t43 = -t109 * pkin(7) + t56;
t19 = -t162 * t43 + t239 * t42;
t20 = t162 * t42 + t239 * t43;
t179 = -qJD(1) * t124 + t190;
t177 = 0.2e1 * t124 * qJD(3) - qJDD(3) * t135;
t76 = pkin(3) * t198 + qJDD(1) * t243 + qJDD(4);
t176 = -t58 * t109 - t181 * t242;
t175 = -t135 * t167 - 0.2e1 * t122 - t196;
t51 = -t92 * qJD(3) - t163 * t121 + t145;
t91 = t147 - t216;
t173 = -t51 * t163 + t50 * t165 + (-t163 * t92 - t165 * t91) * qJD(3);
t172 = t196 + t76;
t151 = pkin(7) - t161;
t149 = t166 * pkin(1);
t140 = sin(t153);
t139 = t148 + pkin(2);
t120 = qJDD(3) * t165 - t167 * t163;
t119 = qJDD(3) * t163 + t167 * t165;
t116 = pkin(2) + t213;
t99 = t102 ^ 2;
t80 = (t110 * pkin(4) + t231) * qJD(3);
t79 = pkin(3) * t211 + pkin(4) * t242;
t75 = t109 * pkin(4) + t243;
t60 = -t106 * qJD(3) + t110 * qJDD(3);
t59 = -t109 * qJDD(3) - t110 * t167;
t34 = t76 + t240;
t29 = -pkin(7) * t181 + t41;
t28 = t106 * pkin(7) + t40;
t15 = -t62 * t150 - t25 * t152;
t8 = -t162 * t23 + t239 * t22;
t4 = -t20 * qJD(5) - t162 * t29 + t239 * t28;
t3 = t19 * qJD(5) + t162 * t28 + t239 * t29;
t5 = [0, 0, 0, 0, 0, qJDD(1), t188, g(1) * t166 + g(2) * t164, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t160 * t219 - t196, -0.2e1 * t158 * t219 + t190, 0, (t188 + (t158 ^ 2 + t160 ^ 2) * t219) * pkin(1), t155 * qJDD(1) + 0.2e1 * t191, 0.2e1 * t163 * t203 - 0.2e1 * t212 * t207, t119, t156 * qJDD(1) - 0.2e1 * t191, t120, 0, t163 * t177 + t165 * t175, -t163 * t175 + t165 * t177, (t155 + t156) * t121 + t173 - t190, t122 * t137 - g(1) * (-t141 * pkin(2) + t143 * pkin(6) - t230) - g(2) * (t143 * pkin(2) + t141 * pkin(6) + t149) + t173 * t135, -t106 * t242 + t58 * t110, t176 + t225, t60, t102 * t181 + t57 * t109, t59, 0, t55 * qJDD(3) + t76 * t109 + t243 * t57 - t196 * t142 + (t100 * t110 + t102 * t231 + t40) * qJD(3), -t56 * qJDD(3) - t100 * t106 + t76 * t110 + t243 * t58 + t196 * t140 + (t231 * t242 - t41) * qJD(3), -t41 * t102 + t35 * t106 - t14 * t109 - t13 * t110 - t181 * t36 - t242 * t40 - t55 * t58 - t56 * t57 - t190, t14 * t56 + t36 * t41 + t13 * t55 + t35 * t40 + t76 * t243 + t100 * pkin(3) * t210 - g(1) * (-t141 * t139 - t143 * t161 - t230) - g(2) * (t143 * t139 - t141 * t161 + t149), t182 * t63 + t183 * t24, -t187 + t226, t185, -t174 * t62 - t25 * t48, t15, 0, -t134 * t196 + t19 * t150 + t4 * t152 - t174 * t75 + t61 * t25 + t34 * t62 - t48 * t80, t133 * t196 - t20 * t150 - t3 * t152 + t182 * t75 - t183 * t80 - t61 * t24 + t34 * t63, -t1 * t62 + t174 * t20 - t182 * t19 + t183 * t4 - t2 * t63 + t8 * t24 - t9 * t25 + t3 * t48 - t190, t1 * t20 + t9 * t3 + t2 * t19 + t8 * t4 + t34 * t75 + t61 * t80 - g(1) * (-t141 * t116 + t151 * t143 - t230) - g(2) * (t143 * t116 + t141 * t151 + t149); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t120, -t119, 0, t50 * t163 + t51 * t165 - g(3) + (-t163 * t91 + t165 * t92) * qJD(3), 0, 0, 0, 0, 0, 0, t59, -t60, -t176 + t225, -t36 * t106 - t13 * t109 + t14 * t110 - t181 * t35 - g(3), 0, 0, 0, 0, 0, 0, t15, -t185, t187 + t226, t1 * t63 - t2 * t62 - t24 * t9 - t25 * t8 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t212 * t168, t204, t202, t203, qJDD(3), -t234 + t145 + (t92 - t215) * qJD(3) + (t179 + t192) * t163, g(3) * t163 + (t91 + t216) * qJD(3) + t179 * t165 + t201, 0, 0, t218, -t99 + t241, (t102 + t200) * qJD(3) + t178, -t218, -t186, qJDD(3), -g(3) * t142 - t38 * qJD(3) - t100 * t242 + t190 * t140 + (qJDD(3) * t159 - t102 * t211) * pkin(3) + t13, g(3) * t140 + t39 * qJD(3) + t100 * t102 + t190 * t142 + (-qJDD(3) * t157 - t211 * t242) * pkin(3) - t14, (t36 + t38) * t242 + (-t35 + t39) * t102 + (-t157 * t57 - t159 * t58) * pkin(3), -t35 * t38 - t36 * t39 + (-t234 + t13 * t159 + t14 * t157 + (-qJD(1) * t100 + t190) * t163) * pkin(3), t227, t248, t249, -t227, t250, t150, t97 * t150 + t152 * t223 + t48 * t79 + t245, -t98 * t150 - t152 * t224 + t183 * t79 + t246, t174 * t98 - t182 * t97 + (t224 + t8) * t48 + (t223 - t9) * t183, t1 * t98 + t2 * t97 - t61 * t79 - g(3) * t213 + t224 * t9 + t223 * t8 - t190 * (-pkin(4) * t140 - t231); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186 + 0.2e1 * t247, (-t102 + t200) * qJD(3) + t178, -t99 - t241, t36 * t102 + t242 * t35 + t172, 0, 0, 0, 0, 0, 0, -t174 - t221, t182 + t220, -t228 - t229, -t183 * t8 - t9 * t48 + t172 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t248, t249, -t227, t250, t150, t9 * t152 + t245, t8 * t152 + t246, 0, 0;];
tau_reg = t5;
