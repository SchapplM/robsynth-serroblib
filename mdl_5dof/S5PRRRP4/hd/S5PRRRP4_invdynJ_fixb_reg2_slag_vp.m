% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:26
% EndTime: 2019-12-05 16:46:29
% DurationCPUTime: 1.55s
% Computational Cost: add. (1777->267), mult. (3097->316), div. (0->0), fcn. (2101->10), ass. (0->172)
t237 = 2 * qJD(4);
t107 = qJDD(2) + qJDD(3);
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t178 = qJD(3) * t118;
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t171 = qJD(1) * qJD(2);
t231 = -qJDD(1) * t116 - t119 * t171;
t101 = t119 * qJDD(1);
t181 = qJD(1) * t116;
t234 = -qJDD(2) * pkin(2) + qJD(3) * t181 + t116 * t171 - t101;
t180 = qJD(1) * t119;
t85 = qJD(2) * pkin(2) + t180;
t149 = t234 * t115 + t231 * t118 - t85 * t178;
t14 = pkin(7) * t107 - t149;
t114 = sin(qJ(4));
t109 = t114 ^ 2;
t117 = cos(qJ(4));
t110 = t117 ^ 2;
t182 = t109 + t110;
t236 = t182 * t14;
t111 = qJ(2) + qJ(3);
t103 = sin(t111);
t113 = cos(pkin(8));
t197 = t103 * t113;
t112 = sin(pkin(8));
t198 = t103 * t112;
t235 = g(1) * t197 + g(2) * t198;
t108 = qJD(2) + qJD(3);
t217 = pkin(4) * t117;
t136 = qJ(5) * t114 + pkin(3) + t217;
t233 = t108 * t136;
t232 = t136 * t107;
t192 = t107 * t110;
t193 = t107 * t109;
t230 = t192 + t193;
t12 = t117 * t14;
t169 = qJDD(4) * qJ(5);
t39 = t115 * t85 + t118 * t181;
t33 = pkin(7) * t108 + t39;
t204 = t114 * t33;
t7 = t169 + t12 + (qJD(5) - t204) * qJD(4);
t10 = t114 * t14;
t174 = qJD(4) * t117;
t199 = qJDD(4) * pkin(4);
t225 = t33 * t174 - t199;
t8 = qJDD(5) + t10 + t225;
t229 = t8 * t114 + t7 * t117;
t59 = t115 * t116 - t118 * t119;
t227 = t108 * t59;
t226 = t182 * t33;
t165 = pkin(2) * t178;
t90 = t115 * t181;
t51 = t118 * t180 - t90;
t147 = -t51 + t165;
t60 = t115 * t119 + t116 * t118;
t224 = t107 * t60 - t108 * t227;
t144 = g(1) * t113 + g(2) * t112;
t175 = qJD(4) * t114;
t153 = qJD(5) + t204;
t208 = qJD(4) * pkin(4);
t21 = t153 - t208;
t177 = qJD(4) * qJ(5);
t203 = t117 * t33;
t22 = t177 + t203;
t223 = t21 * t174 - t22 * t175 + t229;
t176 = qJD(4) * t108;
t183 = -t109 + t110;
t185 = t117 * t107;
t222 = 0.2e1 * t114 * t185 + 0.2e1 * t183 * t176;
t221 = pkin(2) * t116;
t220 = pkin(2) * t118;
t219 = pkin(3) * t107;
t218 = pkin(3) * t108;
t120 = qJD(4) ^ 2;
t216 = pkin(7) * t120;
t96 = g(3) * t103;
t104 = cos(t111);
t213 = g(3) * t104;
t179 = qJD(3) * t115;
t166 = pkin(2) * t179;
t173 = qJD(5) * t114;
t48 = pkin(4) * t175 - qJ(5) * t174 - t173;
t35 = t48 + t166;
t50 = t60 * qJD(1);
t212 = t35 - t50;
t211 = t48 - t39;
t210 = t235 * t117;
t209 = t104 * pkin(3) + t103 * pkin(7);
t205 = t108 * t39;
t98 = pkin(2) * t115 + pkin(7);
t202 = t120 * t98;
t201 = pkin(7) * qJDD(4);
t200 = qJD(4) * t22;
t196 = t104 * t112;
t195 = t104 * t113;
t194 = t104 * t114;
t191 = t108 * t114;
t190 = t108 * t117;
t189 = t112 * t114;
t188 = t112 * t117;
t187 = t113 * t114;
t186 = t113 * t117;
t184 = qJDD(1) - g(3);
t172 = qJDD(4) * t98;
t168 = -g(3) * t194 + t235 * t114;
t167 = g(1) * t195 + g(2) * t196 + t96;
t163 = t108 * t179;
t162 = t108 * t178;
t17 = t231 * t115 - t234 * t118 - t85 * t179;
t15 = -t17 - t219;
t161 = -t15 - t213;
t38 = t118 * t85 - t90;
t157 = t38 * t175 + t39 * t190 + t210;
t156 = t51 * t175 + t50 * t190 + t210;
t155 = qJ(5) * t194 + t104 * t217 + t209;
t154 = t108 * t182;
t151 = t174 * t191;
t32 = -t38 - t218;
t150 = t15 * t114 + t32 * t174 - t168;
t148 = -t50 + t166;
t146 = -pkin(3) * t103 - t221;
t145 = t216 - t219;
t142 = pkin(4) * t114 - qJ(5) * t117;
t24 = t108 * t60;
t141 = -t107 * t59 - t108 * t24;
t99 = -pkin(3) - t220;
t140 = t107 * t99 + t202;
t139 = t114 * t21 + t117 * t22;
t45 = t104 * t188 - t187;
t47 = t104 * t186 + t189;
t138 = g(1) * t47 + g(2) * t45 - t12;
t5 = (t142 * qJD(4) - t173) * t108 - t232 - t17;
t137 = -t216 - t5 + t232;
t53 = -t136 - t220;
t135 = -t107 * t53 - t202 - t5;
t134 = t108 * t53 - t165;
t133 = t108 * t99 - t165;
t44 = t104 * t189 + t186;
t46 = t104 * t187 - t188;
t132 = g(1) * t46 + g(2) * t44 + t114 * t96 - t10;
t131 = t149 + t167;
t130 = t120 * t60 - t141;
t129 = -qJDD(5) + t132;
t128 = -qJDD(4) * t60 + t227 * t237;
t127 = t17 - t213 + t235;
t126 = -g(3) * t119 + t116 * t144;
t125 = t230 * pkin(7) - t38 * t154 - t167;
t124 = t182 * pkin(2) * t162 - t51 * t154 + t230 * t98 - t167;
t123 = (-t114 * t22 + t117 * t21) * qJD(4) + t229;
t122 = t144 * t136 * t103;
t121 = qJD(2) ^ 2;
t106 = t108 ^ 2;
t105 = t119 * pkin(2);
t91 = t114 * t107;
t73 = pkin(7) * t195;
t71 = pkin(7) * t196;
t70 = t114 * t106 * t117;
t69 = qJDD(4) * t117 - t120 * t114;
t68 = qJDD(4) * t114 + t117 * t120;
t54 = t183 * t106;
t49 = t142 * t108;
t41 = -0.2e1 * t151 + t192;
t40 = 0.2e1 * t151 + t193;
t26 = t32 * t175;
t19 = -t38 - t233;
t18 = t19 * t175;
t3 = t224 * t182;
t2 = t114 * t128 - t117 * t130;
t1 = t114 * t130 + t117 * t128;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t184, 0, 0, 0, 0, 0, 0, qJDD(2) * t119 - t116 * t121, -qJDD(2) * t116 - t119 * t121, 0, -g(3) + (t116 ^ 2 + t119 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t141, -t224, 0, -t149 * t60 - t17 * t59 - t227 * t39 - t24 * t38 - g(3), 0, 0, 0, 0, 0, 0, t2, t1, t3, t15 * t59 + t24 * t32 - g(3) + t182 * (t14 * t60 - t227 * t33), 0, 0, 0, 0, 0, 0, t2, t3, -t1, t123 * t60 - t139 * t227 + t19 * t24 + t5 * t59 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t101 + t126, -t184 * t116 + t144 * t119, 0, 0, 0, 0, 0, 0, 0, t107, t108 * t50 + (t107 * t118 - t163) * pkin(2) + t127, t108 * t51 + (-t107 * t115 - t162) * pkin(2) + t131, 0, t38 * t50 - t39 * t51 + (-t115 * t149 + t118 * t17 + (-t115 * t38 + t118 * t39) * qJD(3) + t126) * pkin(2), t40, t222, t68, t41, t69, 0, t26 + (t133 * qJD(4) - t172) * t114 + (-pkin(2) * t163 - t140 + t161) * t117 + t156, (-t172 + (t133 + t51) * qJD(4)) * t117 + (t108 * t148 + t140) * t114 + t150, t124 + t236, t15 * t99 - g(1) * (t113 * t146 + t73) - g(2) * (t112 * t146 + t71) - g(3) * (t105 + t209) + t148 * t32 + t98 * t236 + t147 * t226, t40, t68, -t222, 0, -t69, t41, t18 + (qJD(4) * t134 - t172) * t114 + (-t108 * t35 + t135 - t213) * t117 + t156, t124 + t223, (t172 + (-t134 - t19 - t51) * qJD(4)) * t117 + (-t108 * t212 + t135) * t114 + t168, t5 * t53 - g(1) * (-t113 * t221 + t73) - g(2) * (-t112 * t221 + t71) - g(3) * (t105 + t155) + t212 * t19 + ((qJD(4) * t21 + t7) * t98 + t147 * t22) * t117 + ((t8 - t200) * t98 + t147 * t21) * t114 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t127 + t205, t108 * t38 + t131, 0, 0, t40, t222, t68, t41, t69, 0, t26 + (-pkin(3) * t176 - t201) * t114 + (-t145 + t161) * t117 + t157, (-t201 + (t38 - t218) * qJD(4)) * t117 + (t145 - t205) * t114 + t150, t125 + t236, -t15 * pkin(3) - t32 * t39 - g(1) * (-pkin(3) * t197 + t73) - g(2) * (-pkin(3) * t198 + t71) - g(3) * t209 - t38 * t226 + pkin(7) * t236, t40, t68, -t222, 0, -t69, t41, t18 + (-t136 * t176 - t201) * t114 + (-t108 * t48 + t137 - t213) * t117 + t157, t125 + t223, (t201 + (-t19 - t38 + t233) * qJD(4)) * t117 + (-t108 * t211 + t137) * t114 + t168, pkin(7) * t123 - g(1) * t73 - g(2) * t71 - g(3) * t155 - t136 * t5 - t139 * t38 + t19 * t211 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t54, t91, t70, t185, qJDD(4), -t32 * t191 + t132, (-t108 * t32 + t96) * t117 + t138, 0, 0, -t70, t91, t54, qJDD(4), -t185, t70, 0.2e1 * t199 + (-t114 * t19 + t117 * t49) * t108 + t129, -t142 * t107 + ((t22 - t177) * t114 + (qJD(5) - t21 - t208) * t117) * t108, -t117 * t96 + 0.2e1 * t169 + qJD(5) * t237 + (t114 * t49 + t117 * t19) * t108 - t138, t7 * qJ(5) - t8 * pkin(4) - t19 * t49 - t21 * t203 - g(1) * (-pkin(4) * t46 + qJ(5) * t47) - g(2) * (-pkin(4) * t44 + qJ(5) * t45) + t153 * t22 + t142 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t70, t91, -t106 * t109 - t120, t19 * t191 - t129 - t200 + t225;];
tau_reg = t4;
