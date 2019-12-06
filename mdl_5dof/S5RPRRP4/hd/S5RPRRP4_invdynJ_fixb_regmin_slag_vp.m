% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:55
% EndTime: 2019-12-05 18:07:04
% DurationCPUTime: 2.22s
% Computational Cost: add. (2117->290), mult. (5124->388), div. (0->0), fcn. (3628->10), ass. (0->180)
t125 = sin(pkin(8));
t127 = sin(qJ(4));
t181 = t125 * qJDD(1);
t128 = sin(qJ(3));
t183 = qJDD(1) * t128;
t130 = cos(qJ(4));
t131 = cos(qJ(3));
t201 = t130 * t131;
t179 = qJD(3) + qJD(4);
t79 = t127 * t131 + t130 * t128;
t238 = t179 * t79;
t19 = (qJD(1) * t238 + t127 * t183) * t125 - t181 * t201;
t120 = t125 ^ 2;
t233 = 0.2e1 * t120;
t132 = cos(qJ(1));
t225 = g(3) * t132;
t126 = cos(pkin(8));
t188 = t126 * qJD(1);
t105 = -qJD(3) + t188;
t207 = t125 * t131;
t178 = pkin(7) * t207;
t213 = qJ(2) * t128;
t145 = -t126 * t213 - t178;
t85 = -t126 * pkin(2) - t125 * pkin(6) - pkin(1);
t67 = t85 * qJD(1) + qJD(2);
t59 = t131 * t67;
t38 = t145 * qJD(1) + t59;
t24 = -t105 * pkin(3) + t38;
t212 = qJ(2) * t131;
t104 = t126 * t212;
t194 = qJD(1) * t125;
t175 = t128 * t194;
t39 = -pkin(7) * t175 + qJD(1) * t104 + t128 * t67;
t31 = t127 * t39;
t167 = t130 * t24 - t31;
t156 = t127 * t175;
t193 = qJD(1) * t131;
t174 = t125 * t193;
t54 = t130 * t174 - t156;
t48 = t54 * qJ(5);
t237 = t48 - t167;
t184 = qJDD(1) * qJ(2);
t186 = qJD(1) * qJD(2);
t149 = t184 + t186;
t236 = t149 * t126;
t124 = qJ(3) + qJ(4);
t115 = sin(t124);
t228 = g(1) * t125;
t116 = cos(t124);
t199 = t132 * t116;
t129 = sin(qJ(1));
t205 = t129 * t115;
t60 = t126 * t205 + t199;
t200 = t132 * t115;
t204 = t129 * t116;
t62 = t126 * t200 - t204;
t235 = -g(2) * t60 + g(3) * t62 + t115 * t228;
t234 = t54 ^ 2;
t121 = t126 ^ 2;
t97 = -qJD(4) + t105;
t5 = -t97 * pkin(4) - t237;
t232 = t5 + t237;
t229 = pkin(7) * t125;
t227 = g(2) * t129;
t226 = qJ(2) * t225;
t146 = qJD(1) * t79;
t51 = t125 * t146;
t68 = pkin(3) * t175 + qJ(2) * t194;
t36 = t51 * pkin(4) + qJD(5) + t68;
t224 = t36 * t54;
t223 = t54 * t51;
t222 = t130 * t38 - t31;
t76 = t131 * t85;
t42 = -t178 + t76 + (-pkin(3) - t213) * t126;
t208 = t125 * t128;
t214 = t128 * t85 + t104;
t47 = -pkin(7) * t208 + t214;
t221 = t127 * t42 + t130 * t47;
t206 = t127 * t128;
t78 = t201 - t206;
t220 = (t179 - t188) * t78;
t219 = t126 * t146 - t238;
t218 = t179 * t156;
t83 = t128 * pkin(3) + pkin(4) * t115;
t217 = t126 * t83;
t33 = t130 * t39;
t216 = t51 * qJ(5);
t191 = qJD(3) * t131;
t192 = qJD(2) * t126;
t215 = t131 * t192 + t85 * t191;
t211 = qJD(3) * t67;
t210 = qJDD(1) * pkin(1);
t133 = qJD(1) ^ 2;
t209 = t120 * t133;
t203 = t129 * t128;
t202 = t129 * t131;
t198 = t132 * t128;
t197 = t132 * t131;
t84 = t131 * pkin(3) + pkin(4) * t116;
t100 = t125 * pkin(3) * t191;
t74 = t125 * qJD(2) + t100;
t106 = pkin(3) * t208;
t77 = t125 * qJ(2) + t106;
t196 = t120 + t121;
t123 = t131 ^ 2;
t195 = t128 ^ 2 - t123;
t190 = qJD(4) * t127;
t189 = qJD(4) * t130;
t187 = qJD(3) + t105;
t185 = qJD(1) * qJD(3);
t182 = qJDD(1) * t131;
t180 = t126 * qJDD(1);
t176 = qJ(2) * qJD(3) * t126;
t173 = qJ(2) * t180;
t172 = t131 * t186;
t171 = t131 * t185;
t169 = t127 * t182;
t103 = -qJDD(3) + t180;
t155 = qJD(1) * t176;
t66 = t85 * qJDD(1) + qJDD(2);
t58 = t131 * t66;
t12 = -t103 * pkin(3) + t58 + (-pkin(7) * t181 - t155) * t131 + (-t173 - t211 + (qJD(3) * t229 - t192) * qJD(1)) * t128;
t162 = t126 * t172 + t128 * t66 + t131 * t173 + t67 * t191;
t139 = -t128 * t155 + t162;
t18 = (-t171 - t183) * t229 + t139;
t168 = t130 * t12 - t127 * t18;
t34 = t145 * qJD(3) + t215;
t35 = -t128 * t192 + (-t104 + (-t85 + t229) * t128) * qJD(3);
t166 = -t127 * t34 + t130 * t35;
t165 = -t127 * t38 - t33;
t164 = -t127 * t47 + t130 * t42;
t163 = -t127 * t12 - t130 * t18 - t24 * t189 + t39 * t190;
t161 = t196 * t133;
t44 = qJ(2) * t181 + qJD(1) * t100 + qJDD(1) * t106 + t125 * t186;
t160 = qJD(1) * t187;
t159 = t179 * t131;
t158 = t103 + t180;
t157 = 0.2e1 * t196;
t154 = g(2) * t132 + g(3) * t129;
t153 = -t225 + t227;
t152 = -t127 * t24 - t33;
t150 = qJD(3) * (t105 + t188);
t148 = (-qJ(5) - pkin(7) - pkin(6)) * t125 - t126 * (pkin(2) + t84) - pkin(1);
t147 = t127 * t35 + t130 * t34 + t42 * t189 - t47 * t190;
t144 = -t154 - t210;
t20 = (t169 + (qJD(1) * t159 + t183) * t130) * t125 - t218;
t143 = t20 * pkin(4) + qJDD(5) + t44;
t112 = qJDD(2) - t210;
t142 = -t112 - t144;
t141 = -t105 ^ 2 - t209;
t140 = t157 * t186 + t227;
t138 = t152 * qJD(4) + t168;
t61 = t126 * t204 - t200;
t63 = -t126 * t199 - t205;
t137 = -g(2) * t61 - g(3) * t63 + t116 * t228 + t68 * t51 + t163;
t135 = -t68 * t54 + t138 + t235;
t111 = t130 * pkin(3) + pkin(4);
t91 = -qJDD(4) + t103;
t72 = -t126 * t197 - t203;
t71 = t126 * t198 - t202;
t70 = t126 * t202 - t198;
t69 = t126 * t203 + t197;
t65 = t78 * t125;
t64 = t79 * t125;
t49 = t51 ^ 2;
t30 = (t130 * t159 - t179 * t206) * t125;
t29 = t238 * t125;
t21 = -t49 + t234;
t17 = -qJ(5) * t64 + t221;
t15 = -t126 * pkin(4) - t65 * qJ(5) + t164;
t14 = -t54 * t97 + (-t169 + (-t179 * t193 - t183) * t130) * t125 + t218;
t13 = -t51 * t97 - t19;
t9 = -t48 + t222;
t8 = t165 + t216;
t7 = -t152 - t216;
t4 = t29 * qJ(5) - qJD(4) * t221 - t65 * qJD(5) + t166;
t3 = -t30 * qJ(5) - t64 * qJD(5) + t147;
t2 = -qJ(5) * t20 - qJD(5) * t51 - t163;
t1 = -t91 * pkin(4) + t19 * qJ(5) - t54 * qJD(5) + t138;
t6 = [qJDD(1), t154, -t153, t142 * t126, -t142 * t125, t157 * t184 + t140 - t225, -t226 + (-t112 + t154) * pkin(1) + (t196 * t184 + t140) * qJ(2), (qJDD(1) * t123 - 0.2e1 * t128 * t171) * t120, (-t128 * t182 + t195 * t185) * t233, (t128 * t150 - t158 * t131) * t125, (t158 * t128 + t131 * t150) * t125, t103 * t126, -g(2) * t72 + g(3) * t70 - t76 * t103 - t58 * t126 + (t105 * t126 + (t233 + t121) * qJD(1)) * qJ(2) * t191 + (qJD(3) * t85 * t105 + t149 * t233 + (qJ(2) * t103 + qJD(2) * t105 + t211 + t236) * t126) * t128, (-t128 * t176 + t215) * t105 + t214 * t103 + t139 * t126 - g(2) * t71 - g(3) * t69 + (t172 + (-t128 * t185 + t182) * qJ(2)) * t233, -t19 * t65 - t29 * t54, t19 * t64 - t20 * t65 + t29 * t51 - t30 * t54, t19 * t126 + t29 * t97 - t65 * t91, t20 * t126 + t30 * t97 + t64 * t91, t91 * t126, -t166 * t97 - t164 * t91 - t168 * t126 + t74 * t51 + t77 * t20 + t44 * t64 + t68 * t30 - g(2) * t63 + g(3) * t61 + (-t126 * t152 + t221 * t97) * qJD(4), -g(2) * t62 - g(3) * t60 - t163 * t126 + t147 * t97 - t77 * t19 + t221 * t91 - t68 * t29 + t44 * t65 + t74 * t54, -t1 * t65 + t125 * t154 + t15 * t19 - t17 * t20 - t2 * t64 + t5 * t29 - t3 * t51 - t7 * t30 - t4 * t54, t2 * t17 + t7 * t3 + t1 * t15 + t5 * t4 + t143 * (t64 * pkin(4) + t77) + t36 * (t30 * pkin(4) + t74) - t226 + (-g(2) * t148 - g(3) * t83) * t132 + (-g(2) * (-qJ(2) - t83) - g(3) * t148) * t129; 0, 0, 0, -t180, t181, -t161, -qJ(2) * t161 + qJDD(2) + t144, 0, 0, 0, 0, 0, -t131 * t103 + t141 * t128, t128 * t103 + t141 * t131, 0, 0, 0, 0, 0, -t51 * t194 - t219 * t97 - t78 * t91, -t54 * t194 + t220 * t97 + t79 * t91, t19 * t78 - t20 * t79 - t219 * t54 - t220 * t51, t1 * t78 - t36 * t194 + t2 * t79 + t219 * t5 + t220 * t7 - t154; 0, 0, 0, 0, 0, 0, 0, t131 * t128 * t209, -t195 * t209, (-t128 * t160 + t182) * t125, (-t131 * t160 - t183) * t125, -t103, -g(2) * t69 + g(3) * t71 + t58 + (-t126 * t160 - t209) * t212 + (-t187 * t67 + t228 - t236) * t128, g(1) * t207 - g(2) * t70 - g(3) * t72 - t59 * t105 + (t187 * t188 + t209) * t213 - t162, t223, t21, t13, t14, -t91, t165 * t97 + (-t130 * t91 - t174 * t51 + t97 * t190) * pkin(3) + t135, -t222 * t97 + (t127 * t91 - t174 * t54 + t97 * t189) * pkin(3) + t137, t111 * t19 + (t7 + t8) * t54 + (-t5 + t9) * t51 + (-t127 * t20 + (t127 * t54 - t130 * t51) * qJD(4)) * pkin(3), t1 * t111 - t7 * t9 - t5 * t8 - pkin(4) * t224 + t83 * t228 - g(2) * (t129 * t217 + t132 * t84) - g(3) * (t129 * t84 - t132 * t217) + (-t36 * t174 + t2 * t127 + (-t127 * t5 + t130 * t7) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, t21, t13, t14, -t91, t152 * t97 + t135, -t167 * t97 + t137, pkin(4) * t19 - t232 * t51, t232 * t7 + (t1 - t224 + t235) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49 - t234, g(1) * t126 + t125 * t153 + t5 * t54 + t7 * t51 + t143;];
tau_reg = t6;
