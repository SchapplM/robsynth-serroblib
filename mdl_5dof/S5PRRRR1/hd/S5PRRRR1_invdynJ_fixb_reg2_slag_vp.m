% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR1
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:55
% EndTime: 2019-07-18 13:29:01
% DurationCPUTime: 2.80s
% Computational Cost: add. (2277->346), mult. (6135->492), div. (0->0), fcn. (5001->10), ass. (0->196)
t105 = sin(qJ(4));
t109 = cos(qJ(3));
t110 = cos(qJ(2));
t191 = qJD(1) * qJD(2);
t161 = t110 * t191;
t107 = sin(qJ(2));
t186 = t107 * qJDD(1);
t124 = (t161 + t186) * t109;
t106 = sin(qJ(3));
t199 = qJD(1) * t107;
t167 = t106 * t199;
t132 = -qJD(3) * pkin(2) + t167;
t251 = (qJD(4) * t132 - t124) * t105;
t241 = cos(qJ(4));
t159 = qJDD(2) * t241;
t187 = t106 * qJDD(2);
t183 = qJD(3) + qJD(4);
t211 = t105 * t109;
t64 = t241 * t106 + t211;
t250 = t183 * t64;
t25 = qJD(2) * t250 + t105 * t187 - t109 * t159;
t45 = t241 * t132 + t199 * t211;
t249 = qJD(4) * t45;
t103 = qJ(3) + qJ(4);
t98 = sin(t103);
t218 = t110 * t98;
t236 = g(3) * t107;
t99 = cos(t103);
t94 = g(2) * t99;
t247 = g(1) * t218 + t98 * t236 + t94;
t100 = t106 ^ 2;
t102 = t109 ^ 2;
t201 = t100 + t102;
t246 = t107 * (-0.1e1 + t201);
t163 = t241 * qJD(4);
t245 = t241 * qJD(3) + t163;
t112 = qJD(2) ^ 2;
t244 = qJDD(2) * t107 + t112 * t110;
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t182 = qJDD(3) + qJDD(4);
t212 = t105 * t106;
t143 = t183 * t212;
t185 = t109 * qJDD(2);
t165 = qJD(2) * t241;
t87 = t109 * t165;
t158 = -t105 * t185 - t106 * t159 - t183 * t87;
t24 = qJD(2) * t143 + t158;
t61 = t64 * qJD(2);
t44 = t104 * t183 + t108 * t61;
t15 = t44 * qJD(5) - t104 * t24 - t108 * t182;
t148 = t241 * t186;
t151 = t110 * t165;
t196 = qJD(2) * t110;
t207 = t107 * t109;
t38 = -t106 * t186 + qJDD(3) * pkin(2) + (-qJD(3) * t207 - t106 * t196) * qJD(1);
t243 = (qJD(1) * t151 + t148) * t109 + t105 * t38;
t203 = qJDD(1) - g(3);
t239 = g(1) * t107;
t130 = t203 * t110 + t239;
t93 = g(2) * t98;
t242 = t110 ^ 2 * qJDD(1) - g(3);
t240 = pkin(2) * t106;
t238 = g(1) * t110;
t237 = g(2) * t109;
t235 = g(3) * t110;
t150 = qJD(3) * t167;
t118 = -t241 * t150 - t249;
t12 = t118 + t243;
t169 = t241 * t109;
t152 = qJD(1) * t169;
t136 = t107 * t152;
t46 = -t105 * t132 + t136;
t192 = t110 * qJD(1);
t70 = -t109 * qJD(2) * pkin(2) - t192;
t27 = -t104 * t46 + t108 * t70;
t214 = qJD(5) * t27;
t190 = qJD(2) * qJD(3);
t162 = t106 * t190;
t47 = t107 * t191 - t110 * qJDD(1) + (t162 - t185) * pkin(2);
t3 = t104 * t47 + t108 * t12 + t214;
t2 = t3 * t108;
t156 = t108 * t183;
t42 = t104 * t61 - t156;
t198 = qJD(2) * t106;
t168 = t105 * t198;
t59 = -t87 + t168;
t58 = qJD(5) + t59;
t234 = t42 * t58;
t233 = t44 * t42;
t232 = t44 * t58;
t231 = t45 * t42;
t230 = t45 * t44;
t229 = t58 * t61;
t228 = t61 * t59;
t23 = qJDD(5) + t25;
t227 = t104 * t23;
t226 = t104 * t42;
t225 = t104 * t44;
t224 = t104 * t59;
t222 = t108 * t23;
t221 = t108 * t42;
t220 = t108 * t44;
t219 = t108 * t59;
t217 = t110 * t99;
t194 = qJD(5) * t104;
t14 = -qJD(5) * t156 - t104 * t182 + t108 * t24 + t61 * t194;
t216 = t14 * t104;
t215 = t15 * t108;
t28 = t104 * t70 + t108 * t46;
t213 = qJD(5) * t28;
t210 = t106 * t107;
t209 = t107 * t104;
t208 = t107 * t108;
t206 = t110 * t104;
t205 = t110 * t108;
t202 = t100 - t102;
t111 = qJD(3) ^ 2;
t200 = t111 + t112;
t197 = qJD(2) * t107;
t195 = qJD(3) * t106;
t193 = qJD(5) * t108;
t188 = qJDD(3) * t106;
t184 = t110 * qJDD(2);
t180 = -g(1) * t217 - t99 * t236 + t93;
t179 = pkin(2) * t198;
t123 = qJD(4) * t136 - t105 * t150 - t241 * t38;
t13 = t123 - t251;
t178 = t241 * t13;
t177 = t64 * t194;
t176 = t64 * t193;
t175 = qJD(5) * t109 * t58;
t173 = t106 * t112 * t109;
t36 = t45 * t194;
t37 = t45 * t193;
t172 = t104 * t241;
t171 = t108 * t241;
t170 = (-t58 + t59) * t45;
t49 = t64 * t192;
t166 = g(1) * pkin(2) * t207 - t45 * t49;
t157 = t108 * t58;
t155 = t201 * qJDD(1);
t154 = t58 * t163;
t149 = t109 * t162;
t147 = -t236 - t238;
t34 = -t245 * t109 + t143;
t146 = t13 * t64 - t45 * t34;
t145 = t23 * t64 - t34 * t58;
t142 = t104 * t28 + t108 * t27;
t141 = t104 * t27 - t108 * t28;
t140 = t220 + t226;
t131 = t169 - t212;
t48 = t131 * t199;
t139 = t238 * t240 - t45 * t48 + (g(3) * t210 + t237) * pkin(2);
t138 = -t27 * t219 - t28 * t224 + t180 + t2;
t53 = t131 * t107;
t39 = -t104 * t53 - t205;
t137 = -t108 * t53 + t206;
t21 = t106 * t151 + (t105 * t196 + t245 * t107) * t109 - t183 * t105 * t210;
t52 = t64 * t107;
t135 = t13 * t52 + t45 * t21 - g(3);
t134 = -qJD(5) * t70 - t12 - t93;
t133 = -t109 * t23 + t58 * t195;
t129 = -t203 * t107 + t238;
t126 = t110 * t131;
t125 = -t27 * t61 + t36 + (-t13 + t247) * t108;
t41 = t108 * t47;
t4 = -t104 * t12 - t213 + t41;
t121 = -t142 * qJD(5) - t4 * t104;
t119 = t141 * qJD(5) - t104 * t3 - t108 * t4;
t117 = -t70 * t61 - t123 + t247;
t116 = t28 * t61 + t37 + (t147 * t98 + t13 - t94) * t104;
t114 = t70 * t59 - t118 - t180;
t101 = t107 ^ 2;
t57 = t99 * t205 + t209;
t56 = -t99 * t206 + t208;
t55 = -t99 * t208 + t206;
t54 = t99 * t209 + t205;
t51 = qJD(1) * t126;
t50 = t64 * t199;
t33 = t104 * t199 + t108 * t51;
t32 = -t104 * t51 + t108 * t199;
t30 = t104 * t179 - t108 * t50;
t29 = t104 * t50 + t108 * t179;
t26 = -t59 ^ 2 + t61 ^ 2;
t20 = qJD(2) * t126 - t107 * t250;
t17 = t61 * t183 - t25;
t16 = -t158 + (-t168 + t59) * t183;
t10 = t137 * qJD(5) - t104 * t20 + t108 * t197;
t9 = t39 * qJD(5) + t104 * t197 + t108 * t20;
t8 = t58 * t157 - t44 * t61 + t227;
t7 = -t58 ^ 2 * t104 + t42 * t61 + t222;
t6 = t58 * t226 - t215;
t5 = t44 * t157 - t216;
t1 = (-t14 - t234) * t108 + (-t15 - t232) * t104;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t203, 0, 0, 0, 0, 0, 0, -t112 * t107 + t184, -t244, 0, t101 * qJDD(1) + t242, 0, 0, 0, 0, 0, 0, (-0.2e1 * t162 + t185) * t110 + (-t200 * t109 - t188) * t107, (-qJDD(3) * t107 - 0.2e1 * t110 * t190) * t109 + (t200 * t107 - t184) * t106, t244 * t201, t101 * t155 + 0.2e1 * t161 * t246 + t242, 0, 0, 0, 0, 0, 0, -t110 * t25 - t52 * t182 - t21 * t183 + t59 * t197, t110 * t24 - t53 * t182 - t20 * t183 + t61 * t197, -t20 * t59 + t21 * t61 - t52 * t24 - t53 * t25, -t47 * t110 + t12 * t53 + t70 * t197 + t46 * t20 + t135, 0, 0, 0, 0, 0, 0, t10 * t58 + t52 * t15 + t21 * t42 + t39 * t23, t137 * t23 - t52 * t14 + t21 * t44 - t9 * t58, -t10 * t44 + t137 * t15 + t39 * t14 - t9 * t42, t27 * t10 - t137 * t3 + t28 * t9 + t4 * t39 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t130, t129, 0, 0, t100 * qJDD(2) + 0.2e1 * t149, 0.2e1 * t106 * t185 - 0.2e1 * t202 * t190, t111 * t109 + t188, t102 * qJDD(2) - 0.2e1 * t149, qJDD(3) * t109 - t111 * t106, 0, t130 * t109, -t130 * t106, -t238 + (-g(3) + t155) * t107, -qJD(1) ^ 2 * t110 * t246, -t24 * t64 - t61 * t34, -t131 * t24 - t64 * t25 - t250 * t61 + t34 * t59, t64 * t182 - t34 * t183, -t131 * t25 + t250 * t59, t131 * t182 - t183 * t250, 0, -t47 * t131 + t70 * t250 + t49 * t183 - g(3) * t217 + (g(1) * t99 - qJD(1) * t59) * t107 + (-t109 * t25 + t59 * t195) * pkin(2), t47 * t64 - t70 * t34 + t51 * t183 + g(3) * t218 + (-g(1) * t98 - qJD(1) * t61) * t107 + (t109 * t24 + t61 * t195) * pkin(2), t12 * t131 - t250 * t46 - t49 * t61 + t51 * t59 + t146 + t147, -t70 * t199 - t46 * t51 + (t70 * t195 + (-t47 - t235) * t109) * pkin(2) + t166, -t44 * t177 + (-t14 * t64 - t34 * t44) * t108, (t221 + t225) * t34 + (t216 - t215 + (-t220 + t226) * qJD(5)) * t64, t108 * t145 + t131 * t14 - t177 * t58 + t250 * t44, t42 * t176 + (t15 * t64 - t34 * t42) * t104, -t104 * t145 + t131 * t15 - t176 * t58 - t250 * t42, -t131 * t23 + t250 * t58, t64 * t37 - g(1) * t55 - g(3) * t57 + t27 * t250 - t32 * t58 - t4 * t131 - t49 * t42 + t146 * t104 + (t104 * t175 + t108 * t133) * pkin(2), -t64 * t36 - g(1) * t54 - g(3) * t56 - t28 * t250 + t3 * t131 + t33 * t58 - t49 * t44 + t146 * t108 + (-t104 * t133 + t108 * t175) * pkin(2), t32 * t44 + t33 * t42 + (-t235 + t239) * t98 + t142 * t34 + t119 * t64 + (-t140 * t195 + (t104 * t15 - t14 * t108 + (t221 - t225) * qJD(5)) * t109) * pkin(2), -t27 * t32 - t28 * t33 + (t142 * t195 + (t119 - t235) * t109) * pkin(2) + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t202 * t112, t187, t173, t185, qJDD(3), t106 * t129 + t237, -g(2) * t106 + t109 * t129, 0, 0, t228, t26, t16, -t228, t17, t182, t48 * t183 + (t241 * t182 - t59 * t198) * pkin(2) + ((t167 + (-0.2e1 * qJD(3) - qJD(4)) * pkin(2)) * qJD(4) - t124) * t105 + t117, -t109 * t148 + (-pkin(2) * t182 - t38) * t105 + (-t110 * t152 - t61 * t240) * qJD(2) + t114 + (-pkin(2) * t163 - t50) * t183, (t46 - t48) * t61 + (t45 - t50) * t59 + (t241 * t24 - t105 * t25 + (t105 * t61 - t241 * t59) * qJD(4)) * pkin(2), t46 * t50 + (-t70 * t198 - t178 + t105 * t12 + (t105 * t45 + t241 * t46) * qJD(4)) * pkin(2) + t139, t5, t1, t8, t6, t7, -t229, t45 * t224 - t29 * t58 - t48 * t42 + (-t104 * t154 - t241 * t15 + (qJD(4) * t42 - t193 * t58 - t227) * t105) * pkin(2) + t125, t45 * t219 + t30 * t58 - t48 * t44 + (-t108 * t154 + t241 * t14 + (qJD(4) * t44 + t194 * t58 - t222) * t105) * pkin(2) + t116, t29 * t44 + t30 * t42 + ((-t171 * t42 + t172 * t44) * qJD(4) + (qJD(5) * t140 - t215 - t216) * t105) * pkin(2) + t121 + t138, -t27 * t29 - t28 * t30 + (-t178 + (t171 * t28 - t172 * t27) * qJD(4) + (t121 + t2 + t249) * t105) * pkin(2) + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t26, t16, -t228, t17, t182, t46 * t183 + t117 + t251, -t45 * t183 + t114 - t243, 0, 0, t5, t1, t8, t6, t7, -t229, t104 * t170 - t46 * t42 + t125, t108 * t170 - t46 * t44 + t116, (-t214 - t231) * t108 + (-t213 - t4 + t230) * t104 + t138, (-t141 - t46) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, -t42 ^ 2 + t44 ^ 2, -t14 + t234, -t233, -t15 + t232, t23, -g(1) * t56 + g(3) * t54 + t104 * t134 - t193 * t46 + t28 * t58 - t230 + t41, g(1) * t57 - g(3) * t55 + t27 * t58 + t231 + (qJD(5) * t46 - t47) * t104 + t134 * t108, 0, 0;];
tau_reg  = t11;
