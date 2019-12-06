% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:35
% EndTime: 2019-12-05 16:04:43
% DurationCPUTime: 2.45s
% Computational Cost: add. (2105->332), mult. (4478->463), div. (0->0), fcn. (3294->10), ass. (0->177)
t89 = cos(pkin(5));
t179 = qJD(1) * t89;
t87 = sin(pkin(5));
t180 = qJD(1) * t87;
t95 = cos(qJ(2));
t142 = t95 * t180;
t123 = qJD(3) - t142;
t96 = -pkin(2) - pkin(7);
t53 = t96 * qJD(2) + t123;
t91 = sin(qJ(4));
t94 = cos(qJ(4));
t32 = t94 * t179 + t53 * t91;
t20 = qJD(4) * pkin(8) + t32;
t92 = sin(qJ(2));
t147 = t92 * t180;
t131 = pkin(4) * t91 - pkin(8) * t94;
t62 = qJ(3) + t131;
t37 = t62 * qJD(2) + t147;
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t122 = t20 * t90 - t37 * t93;
t160 = qJDD(1) * t87;
t141 = t92 * t160;
t132 = pkin(4) * t94 + pkin(8) * t91;
t55 = t132 * qJD(4) + qJD(3);
t14 = t141 + t62 * qJDD(2) + (t55 + t142) * qJD(2);
t146 = t91 * t179;
t159 = qJDD(1) * t89;
t170 = qJD(4) * t94;
t140 = t95 * t160;
t197 = t87 * t92;
t145 = qJD(2) * t197;
t67 = qJD(1) * t145;
t113 = qJDD(3) + t67 - t140;
t33 = t96 * qJDD(2) + t113;
t150 = -t94 * t159 - t53 * t170 - t91 * t33;
t7 = -qJD(4) * t146 - t150;
t5 = qJDD(4) * pkin(8) + t7;
t1 = -t122 * qJD(5) + t90 * t14 + t93 * t5;
t177 = qJD(2) * t91;
t76 = qJD(5) + t177;
t224 = t122 * t76 + t1;
t31 = t53 * t94 - t146;
t127 = t91 * t159 - t94 * t33;
t8 = -qJD(4) * t32 - t127;
t100 = -(t31 * t91 - t32 * t94) * qJD(4) + t7 * t91 + t8 * t94;
t195 = t87 * t95;
t193 = t89 * t95;
t86 = sin(pkin(9));
t88 = cos(pkin(9));
t47 = -t88 * t193 + t86 * t92;
t49 = t86 * t193 + t88 * t92;
t220 = g(1) * t49 + g(2) * t47 - g(3) * t195;
t223 = t100 - t220;
t154 = qJD(2) * qJD(4);
t155 = t94 * qJDD(2);
t222 = -t91 * t154 + t155;
t10 = t20 * t93 + t37 * t90;
t2 = -qJD(5) * t10 + t93 * t14 - t90 * t5;
t221 = -t10 * t76 - t2;
t172 = qJD(4) * t90;
t176 = qJD(2) * t94;
t59 = t93 * t176 + t172;
t168 = qJD(5) * t59;
t18 = -t93 * qJDD(4) + t222 * t90 + t168;
t194 = t89 * t92;
t48 = t88 * t194 + t86 * t95;
t50 = -t86 * t194 + t88 * t95;
t219 = -g(1) * t50 - g(2) * t48;
t162 = qJDD(1) - g(3);
t218 = -t162 * t197 - t219;
t163 = t93 * qJD(4);
t165 = qJD(5) * t94;
t110 = t91 * t163 + t90 * t165;
t17 = t110 * qJD(2) - qJD(5) * t163 - t90 * qJDD(4) - t93 * t155;
t161 = qJD(2) * qJ(3);
t61 = t147 + t161;
t217 = qJD(4) * (t147 - t61 - t161) - qJDD(4) * t96;
t216 = pkin(7) * t47;
t215 = pkin(7) * t49;
t211 = t89 ^ 2 * qJDD(1) - g(3);
t209 = t17 * t90;
t208 = t18 * t93;
t153 = qJDD(2) * qJ(3);
t34 = t141 + t153 + (qJD(3) + t142) * qJD(2);
t207 = t34 * t92;
t138 = t94 * t154;
t156 = t91 * qJDD(2);
t54 = qJDD(5) + t138 + t156;
t206 = t54 * t90;
t205 = t54 * t93;
t57 = t90 * t176 - t163;
t204 = t57 * t76;
t203 = t59 * t57;
t202 = t59 * t76;
t201 = t61 * t95;
t200 = t76 * t90;
t199 = t76 * t93;
t198 = t87 * t91;
t196 = t87 * t94;
t192 = t90 * t92;
t191 = t91 * t96;
t190 = t92 * t93;
t189 = t94 * t17;
t188 = t94 * t18;
t98 = qJD(2) ^ 2;
t187 = t95 * t98;
t169 = qJD(4) * t96;
t143 = t94 * t169;
t38 = -t90 * t191 + t62 * t93;
t186 = t38 * qJD(5) + t93 * t143 + t55 * t90 - (t91 * t190 + t90 * t95) * t180;
t39 = t93 * t191 + t62 * t90;
t185 = -t39 * qJD(5) - t90 * t143 + t55 * t93 - (-t91 * t192 + t93 * t95) * t180;
t184 = pkin(2) * t195 + qJ(3) * t197;
t84 = t91 ^ 2;
t85 = t94 ^ 2;
t183 = t84 - t85;
t182 = t84 + t85;
t97 = qJD(4) ^ 2;
t181 = -t97 - t98;
t178 = qJD(2) * t61;
t175 = qJD(2) * t95;
t174 = qJD(4) * t57;
t173 = qJD(4) * t59;
t171 = qJD(4) * t91;
t167 = qJD(5) * t90;
t166 = qJD(5) * t93;
t164 = qJDD(2) * pkin(2);
t158 = qJDD(4) * t91;
t152 = t91 * t195;
t151 = t94 * t98 * t91;
t148 = pkin(7) * t195 + t184;
t144 = t87 * t175;
t43 = t47 * pkin(2);
t137 = qJ(3) * t48 - t43;
t44 = t49 * pkin(2);
t136 = qJ(3) * t50 - t44;
t135 = t182 * qJDD(2);
t134 = qJD(5) * t91 + qJD(2);
t133 = t91 * t138;
t126 = t10 * t93 + t122 * t90;
t125 = t10 * t90 - t122 * t93;
t120 = (-qJD(2) * pkin(2) + t123) * t92 + t201;
t119 = t34 * qJ(3) + t61 * qJD(3);
t118 = qJDD(2) * t92 + t187;
t52 = t89 * t94 - t152;
t27 = t87 * t190 - t52 * t90;
t28 = t87 * t192 + t52 * t93;
t51 = t94 * t195 + t89 * t91;
t115 = t76 * t166 + t206;
t114 = -t76 * t167 + t205;
t112 = -g(1) * (-t86 * t198 + t49 * t94) - g(2) * (t88 * t198 + t47 * t94) + g(3) * t51;
t22 = t86 * t196 + t49 * t91;
t24 = t88 * t196 - t47 * t91;
t111 = -g(1) * t22 + g(2) * t24 - g(3) * t52;
t6 = -qJDD(4) * pkin(4) - t8;
t108 = t112 - t6;
t106 = g(3) * t197 - t219;
t19 = -qJD(4) * pkin(4) - t31;
t105 = -pkin(8) * t54 + t76 * t19;
t104 = t220 + t140;
t103 = qJDD(3) - t104;
t102 = pkin(8) * qJD(5) * t76 - t108;
t101 = -t125 * qJD(5) + t1 * t93 - t2 * t90;
t99 = t123 * qJD(2) - t96 * t97 - t106 + t153 + t34;
t81 = qJDD(4) * t94;
t60 = t132 * qJD(2);
t46 = t118 * t87;
t45 = (-qJDD(2) * t95 + t92 * t98) * t87;
t40 = t113 - t164;
t26 = -qJD(4) * t152 - t94 * t145 + t89 * t170;
t25 = -t51 * qJD(4) + t91 * t145;
t16 = t31 * t93 + t60 * t90;
t15 = -t31 * t90 + t60 * t93;
t4 = t27 * qJD(5) + t90 * t144 + t25 * t93;
t3 = -t28 * qJD(5) + t93 * t144 - t25 * t90;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t162, 0, 0, 0, 0, 0, 0, -t45, -t46, 0, (t92 ^ 2 + t95 ^ 2) * t87 ^ 2 * qJDD(1) + t211, 0, 0, 0, 0, 0, 0, 0, t45, t46, (t120 * qJD(2) - t40 * t95 + t207) * t87 + t211, 0, 0, 0, 0, 0, 0, -qJD(4) * t26 - qJDD(4) * t51 + (t118 * t91 + t138 * t92) * t87, -qJD(4) * t25 - qJDD(4) * t52 + (t94 * t187 + t222 * t92) * t87, (t51 * t94 - t52 * t91) * qJDD(2) + (-t25 * t91 + t26 * t94 + (-t51 * t91 - t52 * t94) * qJD(4)) * qJD(2), t25 * t32 - t26 * t31 - t51 * t8 + t52 * t7 - g(3) + (t61 * t175 + t207) * t87, 0, 0, 0, 0, 0, 0, t18 * t51 + t26 * t57 + t27 * t54 + t3 * t76, -t17 * t51 + t26 * t59 - t28 * t54 - t4 * t76, t17 * t27 - t18 * t28 - t3 * t59 - t4 * t57, t1 * t28 + t10 * t4 - t122 * t3 + t19 * t26 + t2 * t27 + t51 * t6 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t104, t218, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t103 - 0.2e1 * t164, 0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t153 - t218, -t40 * pkin(2) - g(1) * t136 - g(2) * t137 - g(3) * t184 - t120 * t180 + t119, qJDD(2) * t85 - 0.2e1 * t133, 0.2e1 * t183 * t154 - 0.2e1 * t91 * t155, -t91 * t97 + t81, qJDD(2) * t84 + 0.2e1 * t133, -t94 * t97 - t158, 0, -t217 * t94 + t99 * t91, t217 * t91 + t99 * t94, -t96 * t135 + t182 * t67 - t223, -g(1) * (t136 - t215) - g(2) * (t137 - t216) - g(3) * t148 + (-t201 + (-t31 * t94 - t32 * t91) * t92) * t180 + t100 * t96 + t119, -t110 * t59 - t93 * t189, (t57 * t93 + t59 * t90) * t171 + (t209 - t208 + (t57 * t90 - t59 * t93) * qJD(5)) * t94, (-t76 * t163 - t17) * t91 + (t114 + t173) * t94, t90 * t188 + (t93 * t165 - t90 * t171) * t57, (t76 * t172 - t18) * t91 + (-t115 - t174) * t94, t76 * t170 + t54 * t91, t38 * t54 + t185 * t76 + t220 * t90 + (t2 + (-t19 * t90 + t57 * t96) * qJD(4) - t106 * t93) * t91 + (-qJD(4) * t122 + t147 * t57 + t19 * t166 - t96 * t18 + t6 * t90) * t94, -t39 * t54 - t186 * t76 + t220 * t93 + (-t1 + (-t19 * t93 + t59 * t96) * qJD(4) + t106 * t90) * t91 + (-t10 * qJD(4) + t147 * t59 - t19 * t167 + t96 * t17 + t6 * t93) * t94, t17 * t38 - t18 * t39 - t185 * t59 - t186 * t57 + t125 * t171 + (-qJD(5) * t126 - t1 * t90 - t2 * t93 + t106) * t94, t1 * t39 + t2 * t38 - t6 * t94 * t96 - g(1) * (-t44 - t215) - g(2) * (-t43 - t216) - g(3) * (t131 * t197 + t148) - t185 * t122 + (t147 * t94 + t91 * t169) * t19 + t186 * t10 + t219 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t98, t103 - t164 + t67 - t178, 0, 0, 0, 0, 0, 0, t181 * t91 + t81, t181 * t94 - t158, -t135, -t178 + t223, 0, 0, 0, 0, 0, 0, -t188 + (t174 - t206) * t91 + (-t134 * t93 - t90 * t170) * t76, t189 + (t173 - t205) * t91 + (t134 * t90 - t94 * t163) * t76, (t134 * t59 - t57 * t170 - t18 * t91) * t93 + (t134 * t57 - t17 * t91 + t59 * t170) * t90, -t125 * qJD(2) + (qJD(4) * t126 - t6) * t94 + (qJD(4) * t19 + t101) * t91 - t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, -t183 * t98, t155, -t151, -t156, qJDD(4), -t61 * t176 + t112 - t127, t61 * t177 + (t31 + t146) * qJD(4) - t111 + t150, 0, 0, t59 * t199 - t209, (-t17 - t204) * t93 + (-t18 - t202) * t90, (t91 * t199 - t59 * t94) * qJD(2) + t115, t57 * t200 - t208, (-t91 * t200 + t57 * t94) * qJD(2) + t114, -t76 * t176, -pkin(4) * t18 - t102 * t93 + t105 * t90 + t122 * t176 - t15 * t76 - t32 * t57, pkin(4) * t17 + t10 * t176 + t102 * t90 + t105 * t93 + t16 * t76 - t32 * t59, t15 * t59 + t16 * t57 + ((-t18 + t168) * pkin(8) + t224) * t93 + ((qJD(5) * t57 - t17) * pkin(8) + t221) * t90 + t111, -t10 * t16 + t122 * t15 - t19 * t32 + t108 * pkin(4) + (t101 + t111) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, -t57 ^ 2 + t59 ^ 2, -t17 + t204, -t203, t202 - t18, t54, -t19 * t59 - g(1) * (-t22 * t90 + t50 * t93) - g(2) * (t24 * t90 + t48 * t93) - g(3) * t27 - t221, t19 * t57 - g(1) * (-t22 * t93 - t50 * t90) - g(2) * (t24 * t93 - t48 * t90) + g(3) * t28 - t224, 0, 0;];
tau_reg = t9;
