% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP13
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP13_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:42
% EndTime: 2019-12-31 18:59:48
% DurationCPUTime: 2.69s
% Computational Cost: add. (2713->384), mult. (5218->466), div. (0->0), fcn. (3046->6), ass. (0->196)
t108 = sin(qJ(3));
t199 = t108 * qJD(1);
t80 = qJD(4) + t199;
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t111 = cos(qJ(3));
t204 = qJD(3) * t108;
t194 = t111 * qJDD(1);
t197 = qJD(1) * qJD(3);
t225 = t107 * t108;
t200 = t107 * qJD(3);
t205 = qJD(1) * t111;
t66 = t110 * t205 + t200;
t24 = qJD(4) * t66 - t110 * qJDD(3) + t107 * t194 - t197 * t225;
t231 = t24 * t110;
t202 = qJD(4) * t107;
t180 = t111 * t202;
t198 = t110 * qJD(3);
t23 = qJD(1) * (t108 * t198 + t180) - qJD(4) * t198 - t107 * qJDD(3) - t110 * t194;
t232 = t23 * t107;
t236 = t110 * t66;
t64 = t107 * t205 - t198;
t267 = (qJD(4) * (t107 * t64 - t236) - t231 + t232) * t111 + (t107 * t66 + t110 * t64) * t204;
t147 = t66 * t80;
t148 = t64 * t80;
t266 = (t24 + t147) * t107 + (t23 + t148) * t110;
t105 = t108 ^ 2;
t106 = t111 ^ 2;
t207 = t105 + t106;
t113 = -pkin(1) - pkin(6);
t76 = t113 * qJDD(1) + qJDD(2);
t174 = t207 * t76;
t99 = t111 * pkin(7);
t265 = t108 * pkin(3) + qJ(2) - t99;
t263 = t24 - t147;
t261 = t111 * t24 - t64 * t204;
t112 = cos(qJ(1));
t102 = g(2) * t112;
t109 = sin(qJ(1));
t103 = g(1) * t109;
t259 = t103 - t102;
t78 = t113 * qJD(1) + qJD(2);
t58 = -qJD(3) * pkin(3) - t111 * t78;
t18 = t64 * pkin(4) - t66 * qJ(5) + t58;
t177 = t111 * t197;
t195 = t108 * qJDD(1);
t62 = qJDD(4) + t177 + t195;
t249 = pkin(7) * t62;
t258 = t18 * t80 - t249;
t191 = pkin(7) * qJD(4) * t80;
t246 = g(3) * t108;
t257 = -t111 * t259 - t191 + t246;
t136 = -t110 * t62 + t80 * t202;
t256 = qJD(1) * (-t111 * t64 + t80 * t225) + t136;
t220 = t108 * t113;
t239 = t107 * t265 + t110 * t220;
t162 = pkin(3) * t111 + pkin(7) * t108;
t63 = t162 * qJD(3) + qJD(2);
t254 = -qJD(4) * t239 + t110 * t63;
t201 = qJD(4) * t110;
t137 = t107 * t62 + t80 * t201;
t253 = t108 * (t80 * t200 - t24) - (qJD(3) * t64 + t137) * t111;
t185 = t80 * t198;
t226 = qJD(3) * t66;
t252 = t108 * (t136 + t226) - t111 * (-t23 + t185) + t107 * qJD(1) * t80;
t251 = t66 ^ 2;
t250 = 0.2e1 * qJ(2);
t248 = pkin(7) * t66;
t247 = t62 * pkin(4);
t245 = g(3) * t111;
t51 = t265 * qJD(1);
t70 = t108 * t78;
t57 = qJD(3) * pkin(7) + t70;
t22 = t107 * t51 + t110 * t57;
t17 = t80 * qJ(5) + t22;
t244 = t17 * t80;
t243 = t22 * t80;
t242 = t66 * t64;
t156 = pkin(4) * t107 - qJ(5) * t110;
t241 = -t107 * qJD(5) + t80 * t156 - t70;
t216 = t110 * t111;
t69 = t162 * qJD(1);
t33 = t107 * t69 + t78 * t216;
t215 = t111 * t112;
t190 = g(2) * t215;
t222 = t108 * t110;
t240 = g(3) * t222 + t110 * t190;
t217 = t109 * t111;
t238 = g(1) * t215 + g(2) * t217;
t235 = t110 * t69;
t234 = t110 * t265;
t230 = t62 * qJ(5);
t229 = t112 * pkin(1) + t109 * qJ(2);
t188 = 0.2e1 * qJD(1) * qJD(2);
t228 = (qJDD(1) * qJ(2) + t188) * qJ(2);
t227 = pkin(1) * qJDD(1);
t224 = t107 * t111;
t223 = t108 * t109;
t221 = t108 * t112;
t219 = t109 * t107;
t218 = t109 * t110;
t214 = t112 * t107;
t213 = t112 * t110;
t114 = qJD(3) ^ 2;
t212 = t113 * t114;
t115 = qJD(1) ^ 2;
t211 = t115 * qJ(2);
t21 = -t107 * t57 + t110 * t51;
t210 = qJD(5) - t21;
t208 = t105 - t106;
t206 = -t114 - t115;
t203 = qJD(3) * t111;
t196 = qJDD(3) * t108;
t181 = t113 * t203;
t193 = t107 * t63 + t110 * t181 + t201 * t265;
t88 = g(2) * t221;
t189 = t64 * t216;
t187 = t112 * pkin(6) + t229;
t186 = t64 ^ 2 - t251;
t183 = t111 * t115 * t108;
t182 = t80 * t205;
t179 = t113 * t202;
t178 = -t88 + t245;
t175 = t107 * t113 - pkin(4);
t29 = qJD(1) * t63 + qJDD(1) * t265;
t38 = qJDD(3) * pkin(7) + t108 * t76 + t78 * t203;
t173 = t107 * t38 - t110 * t29 + t57 * t201 + t51 * t202;
t172 = qJD(4) * t64 - t23;
t171 = t207 * qJDD(1);
t170 = qJDD(2) - t227;
t169 = qJD(4) * t108 + qJD(1);
t168 = t108 * t177;
t52 = t108 * t219 - t213;
t54 = t108 * t214 + t218;
t167 = g(1) * t54 + g(2) * t52;
t53 = t108 * t218 + t214;
t55 = t108 * t213 - t219;
t166 = -g(1) * t55 - g(2) * t53;
t98 = t112 * qJ(2);
t165 = t113 * t109 + t98;
t161 = g(1) * t112 + g(2) * t109;
t159 = -t259 - t211;
t158 = t172 * pkin(7);
t157 = t110 * pkin(4) + t107 * qJ(5);
t16 = -t80 * pkin(4) + t210;
t155 = t107 * t17 - t110 * t16;
t154 = t107 * t16 + t110 * t17;
t153 = t107 * t22 + t110 * t21;
t152 = t107 * t21 - t110 * t22;
t145 = pkin(3) + t157;
t142 = -g(1) * (pkin(3) * t217 + pkin(7) * t223) - g(3) * t99;
t141 = qJDD(1) * t250 + t188;
t140 = t211 - t76 + t103;
t37 = -qJDD(3) * pkin(3) - t111 * t76 + t78 * t204;
t138 = -t113 + t156;
t4 = t107 * t29 + t110 * t38 + t51 * t201 - t57 * t202;
t134 = -g(1) * t217 - t191;
t133 = qJDD(3) * t113 + t197 * t250;
t132 = pkin(3) * t223 - pkin(7) * t217 + t187;
t131 = -pkin(7) * t231 - g(1) * t223 - t178;
t129 = t80 * t58 - t249;
t128 = pkin(3) * t221 - pkin(7) * t215 + t165;
t125 = g(1) * t52 - g(2) * t54 + g(3) * t224 - t173;
t124 = t107 * t148 - t231;
t123 = t141 - t161;
t1 = t80 * qJD(5) + t230 + t4;
t2 = qJDD(5) + t173 - t247;
t122 = -qJD(4) * t155 + t1 * t110 + t2 * t107;
t121 = -qJD(4) * t153 + t107 * t173 + t4 * t110;
t120 = t18 * t66 + qJDD(5) - t125;
t119 = qJD(4) * t189 + t261 * t107;
t118 = -g(1) * t53 + g(2) * t55 - g(3) * t216 + t4;
t117 = -t62 * t225 + (-t169 * t110 - t111 * t200) * t80 - t261;
t116 = -t24 * t222 - qJD(3) * t189 + t169 * t236 + (qJD(1) * t64 + t172 * t108 + t66 * t203) * t107;
t95 = qJDD(3) * t111;
t43 = t138 * t111;
t41 = -t107 * t220 + t234;
t40 = t175 * t108 - t234;
t39 = t108 * qJ(5) + t239;
t36 = t62 * t108 + t80 * t203;
t34 = t66 * pkin(4) + t64 * qJ(5);
t32 = -t78 * t224 + t235;
t28 = -t235 + (-pkin(4) * qJD(1) + t107 * t78) * t111;
t27 = qJ(5) * t205 + t33;
t15 = (qJD(4) * t157 - qJD(5) * t110) * t111 - t138 * t204;
t14 = -t107 * t181 + t254;
t13 = -t108 * t179 + t193;
t12 = t148 - t23;
t11 = t175 * t203 - t254;
t10 = (-t111 * t66 + t80 * t222) * qJD(1) + t137;
t9 = qJ(5) * t203 + (qJD(5) - t179) * t108 + t193;
t8 = t110 * t147 - t232;
t7 = -t66 * t180 + (-t111 * t23 - t66 * t204) * t110;
t6 = (-t23 - t185) * t108 + (-t136 + t226) * t111;
t3 = t24 * pkin(4) + t23 * qJ(5) - t66 * qJD(5) + t37;
t5 = [0, 0, 0, 0, 0, qJDD(1), t259, t161, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t259 - 0.2e1 * t227, t123, -t170 * pkin(1) - g(1) * (-t109 * pkin(1) + t98) - g(2) * t229 + t228, t106 * qJDD(1) - 0.2e1 * t168, -0.2e1 * t108 * t194 + 0.2e1 * t208 * t197, -t114 * t108 + t95, t105 * qJDD(1) + 0.2e1 * t168, -t114 * t111 - t196, 0, t133 * t111 + (t123 - t212) * t108, -t133 * t108 + (t141 - t212) * t111 - t238, -t113 * t171 - t174 + t259, -g(1) * t165 - g(2) * t187 + t113 * t174 + t228, t7, t267, t6, t119, t253, t36, t14 * t80 + t41 * t62 + (-t173 + (-t107 * t58 + t113 * t64) * qJD(3)) * t108 + (qJD(3) * t21 + t37 * t107 - t113 * t24 + t201 * t58) * t111 + t166, -t13 * t80 - t239 * t62 + (-t4 + (-t110 * t58 + t113 * t66) * qJD(3)) * t108 + (-qJD(3) * t22 + t37 * t110 + t113 * t23 - t202 * t58) * t111 + t167, -t13 * t64 - t14 * t66 + t41 * t23 - t239 * t24 + t153 * t204 + (qJD(4) * t152 - t107 * t4 + t110 * t173) * t111 + t238, t4 * t239 + t22 * t13 - t173 * t41 + t21 * t14 - g(1) * t128 - g(2) * t132 + (-t37 * t111 + t204 * t58) * t113, t7, t6, -t267, t36, -t253, t119, -t11 * t80 + t15 * t64 + t43 * t24 - t40 * t62 + (-t18 * t200 - t2) * t108 + (-qJD(3) * t16 + t3 * t107 + t18 * t201) * t111 + t166, t11 * t66 - t40 * t23 - t39 * t24 - t9 * t64 + t155 * t204 + (-qJD(4) * t154 - t1 * t107 + t110 * t2) * t111 + t238, -t15 * t66 + t43 * t23 + t39 * t62 + t9 * t80 + (t18 * t198 + t1) * t108 + (qJD(3) * t17 - t3 * t110 + t18 * t202) * t111 - t167, t1 * t39 + t17 * t9 + t3 * t43 + t18 * t15 + t2 * t40 + t16 * t11 - g(1) * (t55 * pkin(4) + t54 * qJ(5) + t128) - g(2) * (t53 * pkin(4) + t52 * qJ(5) + t132); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t115, t159 + t170, 0, 0, 0, 0, 0, 0, t206 * t108 + t95, t206 * t111 - t196, -t171, t174 + t159, 0, 0, 0, 0, 0, 0, t117, t252, t116, -t153 * qJD(1) + (-qJD(3) * t152 - t37) * t111 + (qJD(3) * t58 + t121) * t108 - t259, 0, 0, 0, 0, 0, 0, t117, t116, -t252, -t155 * qJD(1) + (qJD(3) * t154 - t3) * t111 + (qJD(3) * t18 + t122) * t108 - t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, -t208 * t115, t194, -t183, -t195, qJDD(3), t246 + (-t140 + t102) * t111, t108 * t140 + t178, 0, 0, t8, -t266, t10, t124, -t256, -t182, -t21 * t205 - t64 * t70 - pkin(3) * t24 - t32 * t80 + (t134 - t37) * t110 + t129 * t107 + t240, t22 * t205 - t66 * t70 + pkin(3) * t23 + t33 * t80 + t129 * t110 + (-t257 + t37) * t107, t32 * t66 + t33 * t64 + (-t21 * t199 + t4 + (-t21 + t248) * qJD(4)) * t110 + (t158 + t173 - t243) * t107 + t131, -t22 * t33 - t21 * t32 - t58 * t70 + (t190 - t37 + t246) * pkin(3) + (t121 + t88) * pkin(7) + t142, t8, t10, t266, -t182, t256, t124, t16 * t205 - t145 * t24 + t28 * t80 + t241 * t64 + (t134 - t3) * t110 + t258 * t107 + t240, t27 * t64 - t28 * t66 + (t16 * t199 + t1 + (t16 + t248) * qJD(4)) * t110 + (t158 + t2 - t244) * t107 + t131, -t17 * t205 - t145 * t23 - t27 * t80 - t241 * t66 - t258 * t110 + (t257 - t3) * t107, -t17 * t27 - t16 * t28 + t241 * t18 - t157 * t103 * t111 + (t122 + t88) * pkin(7) + t142 + (t102 * t111 + t246 - t3) * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, -t186, t12, -t242, -t263, t62, -t58 * t66 + t125 + t243, t21 * t80 + t58 * t64 - t118, 0, 0, t242, t12, t186, t62, t263, -t242, -t34 * t64 - t120 + t243 + 0.2e1 * t247, pkin(4) * t23 - t24 * qJ(5) + (t17 - t22) * t66 + (t16 - t210) * t64, 0.2e1 * t230 - t18 * t64 + t34 * t66 + (0.2e1 * qJD(5) - t21) * t80 + t118, t1 * qJ(5) - t2 * pkin(4) - t18 * t34 - t16 * t22 - g(1) * (-pkin(4) * t52 + qJ(5) * t53) - g(2) * (pkin(4) * t54 - qJ(5) * t55) + t210 * t17 + t156 * t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 + t242, t12, -t80 ^ 2 - t251, t120 - t244 - t247;];
tau_reg = t5;
