% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:45
% EndTime: 2019-12-31 19:10:54
% DurationCPUTime: 3.30s
% Computational Cost: add. (6818->358), mult. (17598->487), div. (0->0), fcn. (13256->8), ass. (0->175)
t155 = cos(pkin(9));
t237 = cos(qJ(3));
t192 = t237 * t155;
t141 = qJD(1) * t192;
t154 = sin(pkin(9));
t158 = sin(qJ(3));
t205 = t158 * t154;
t188 = qJD(1) * t205;
t120 = -t141 + t188;
t156 = sin(qJ(5));
t159 = cos(qJ(4));
t157 = sin(qJ(4));
t236 = cos(qJ(5));
t191 = t236 * t157;
t132 = t156 * t159 + t191;
t244 = qJD(4) + qJD(5);
t96 = t244 * t132;
t230 = -t132 * t120 - t96;
t190 = t236 * t159;
t210 = t156 * t157;
t168 = t190 - t210;
t186 = t236 * qJD(5);
t245 = t236 * qJD(4) + t186;
t229 = -t168 * t120 - t245 * t159 + t244 * t210;
t128 = t237 * t154 + t158 * t155;
t122 = t128 * qJD(1);
t105 = -t159 * qJD(3) + t157 * t122;
t107 = t157 * qJD(3) + t159 * t122;
t169 = -t156 * t105 + t236 * t107;
t49 = t236 * t105 + t156 * t107;
t235 = t49 * t169;
t238 = -pkin(8) - pkin(7);
t194 = qJD(4) * t238;
t208 = t157 * t120;
t234 = pkin(6) + qJ(2);
t136 = t234 * t154;
t129 = qJD(1) * t136;
t137 = t234 * t155;
t130 = qJD(1) * t137;
t247 = -t237 * t129 - t158 * t130;
t87 = t122 * pkin(3) + t120 * pkin(7);
t41 = t157 * t87 + t159 * t247;
t257 = -pkin(8) * t208 + t157 * t194 - t41;
t40 = -t157 * t247 + t159 * t87;
t256 = t122 * pkin(4) + t40 + (pkin(8) * t120 - t194) * t159;
t140 = qJD(3) * t141;
t172 = qJD(3) * t188 - t140;
t255 = -qJD(3) * qJD(4) + t172;
t202 = qJD(4) * t157;
t254 = t202 + t208;
t253 = t169 ^ 2 - t49 ^ 2;
t117 = qJD(4) + t120;
t115 = qJD(5) + t117;
t201 = qJD(4) * t159;
t195 = t122 * t201 - t255 * t157;
t200 = qJD(5) * t156;
t61 = t122 * t202 + t255 * t159;
t21 = t105 * t186 + t107 * t200 + t156 * t195 + t236 * t61;
t252 = t49 * t115 - t21;
t146 = -t155 * pkin(2) - pkin(1);
t135 = t146 * qJD(1) + qJD(2);
t65 = t120 * pkin(3) - t122 * pkin(7) + t135;
t91 = -t158 * t129 + t237 * t130;
t86 = qJD(3) * pkin(7) + t91;
t35 = -t157 * t86 + t159 * t65;
t29 = -t107 * pkin(8) + t35;
t25 = t117 * pkin(4) + t29;
t36 = t157 * t65 + t159 * t86;
t30 = -t105 * pkin(8) + t36;
t125 = t128 * qJD(3);
t116 = qJD(1) * t125;
t170 = t192 - t205;
t163 = t170 * qJD(2);
t53 = qJD(1) * t163 + qJD(3) * t247;
t73 = t116 * pkin(3) + t172 * pkin(7);
t19 = -qJD(4) * t36 - t157 * t53 + t159 * t73;
t6 = t116 * pkin(4) + t61 * pkin(8) + t19;
t18 = t157 * t73 + t159 * t53 + t65 * t201 - t86 * t202;
t9 = -t195 * pkin(8) + t18;
t162 = -t156 * t6 - t25 * t186 + t30 * t200 - t236 * t9;
t85 = -qJD(3) * pkin(3) - t247;
t44 = t105 * pkin(4) + t85;
t251 = t44 * t49 + t162;
t175 = -t35 * t117 + t18;
t249 = t36 * t117 + t19;
t182 = t117 * t157;
t248 = t107 * t182;
t246 = -t237 * t136 - t158 * t137;
t196 = t236 * t30;
t8 = t156 * t25 + t196;
t2 = -t8 * qJD(5) - t156 * t9 + t236 * t6;
t243 = -t44 * t169 + t2;
t22 = t169 * qJD(5) - t156 * t61 + t236 * t195;
t242 = t115 * t169 - t22;
t241 = -t229 * t115 + t132 * t116;
t240 = t168 * t21 - t169 * t230;
t239 = t122 ^ 2;
t138 = t238 * t157;
t139 = t238 * t159;
t104 = t156 * t138 - t236 * t139;
t233 = t104 * qJD(5) + t257 * t156 + t256 * t236;
t103 = t236 * t138 + t156 * t139;
t232 = -t103 * qJD(5) + t256 * t156 - t257 * t236;
t199 = t105 * qJD(4);
t57 = t157 * t195;
t231 = -t159 * t199 - t57;
t89 = -pkin(3) * t170 - t128 * pkin(7) + t146;
t102 = -t158 * t136 + t237 * t137;
t94 = t159 * t102;
t43 = t157 * t89 + t94;
t228 = t122 * t49;
t226 = t156 * t30;
t223 = t169 * t122;
t164 = t128 * qJD(2);
t54 = qJD(1) * t164 + t91 * qJD(3);
t222 = t54 * t246;
t221 = t54 * t157;
t220 = t54 * t159;
t219 = t61 * t157;
t218 = t105 * t120;
t217 = t107 * t105;
t216 = t107 * t122;
t92 = t116 * t170;
t215 = t122 * t105;
t214 = t122 * t120;
t213 = t128 * t157;
t212 = t128 * t159;
t209 = t157 * t116;
t124 = t170 * qJD(3);
t207 = t157 * t124;
t110 = t159 * t116;
t204 = t159 * t124;
t203 = t154 ^ 2 + t155 ^ 2;
t198 = qJD(1) * qJD(2);
t189 = t128 * t202;
t66 = t246 * qJD(3) + t163;
t88 = t125 * pkin(3) - t124 * pkin(7);
t184 = -t157 * t66 + t159 * t88;
t42 = -t157 * t102 + t159 * t89;
t183 = t203 * qJD(1) ^ 2;
t181 = t117 * t159;
t180 = -t132 * t22 + t229 * t49;
t179 = t230 * t115 + t168 * t116;
t178 = t254 * pkin(4) - t91;
t177 = t195 * t159;
t176 = t36 * t157 + t35 * t159;
t174 = 0.2e1 * t203 * t198;
t173 = -t254 * t117 + t110;
t33 = -pkin(4) * t170 - pkin(8) * t212 + t42;
t37 = -pkin(8) * t213 + t43;
t16 = -t156 * t37 + t236 * t33;
t17 = t156 * t33 + t236 * t37;
t167 = t128 * t201 + t207;
t166 = -t189 + t204;
t23 = -t102 * t202 + t157 * t88 + t159 * t66 + t89 * t201;
t165 = -pkin(7) * t116 + t117 * t85;
t67 = t102 * qJD(3) + t164;
t150 = -t159 * pkin(4) - pkin(3);
t118 = t120 ^ 2;
t80 = t168 * t128;
t79 = t132 * t128;
t68 = pkin(4) * t213 - t246;
t39 = t167 * pkin(4) + t67;
t32 = t195 * pkin(4) + t54;
t27 = t124 * t191 - t156 * t189 - t200 * t213 + (t124 * t156 + t245 * t128) * t159;
t26 = -t124 * t190 + t96 * t128 + t156 * t207;
t24 = -t43 * qJD(4) + t184;
t15 = -t167 * pkin(8) + t23;
t12 = -pkin(8) * t204 + t125 * pkin(4) + (-t94 + (pkin(8) * t128 - t89) * t157) * qJD(4) + t184;
t11 = t236 * t29 - t226;
t10 = -t156 * t29 - t196;
t7 = t236 * t25 - t226;
t4 = -t17 * qJD(5) + t236 * t12 - t156 * t15;
t3 = t16 * qJD(5) + t156 * t12 + t236 * t15;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, qJ(2) * t174, t122 * t124 - t172 * t128, -t128 * t116 - t124 * t120 - t122 * t125 - t170 * t172, t124 * qJD(3), t120 * t125 - t92, -t125 * qJD(3), 0, -t67 * qJD(3) + t146 * t116 + t135 * t125, -t66 * qJD(3) + t135 * t124 - t146 * t172, -t102 * t116 - t66 * t120 + t67 * t122 - t124 * t247 - t91 * t125 + t54 * t128 + t170 * t53 + t172 * t246, t53 * t102 - t247 * t67 + t91 * t66 - t222, t166 * t107 - t61 * t212, -(t159 * t105 + t107 * t157) * t124 + (-t177 + t219 + (t105 * t157 - t107 * t159) * qJD(4)) * t128, t107 * t125 + t128 * t110 + t166 * t117 + t170 * t61, t167 * t105 + t128 * t57, -t105 * t125 - t167 * t117 - t128 * t209 + t170 * t195, t117 * t125 - t92, t24 * t117 + t42 * t116 - t19 * t170 + t35 * t125 + t67 * t105 - t246 * t195 + t85 * t207 + (t201 * t85 + t221) * t128, t85 * t204 + t246 * t61 + t67 * t107 - t43 * t116 - t23 * t117 - t36 * t125 + t18 * t170 + (-t202 * t85 + t220) * t128, -t23 * t105 - t43 * t195 - t24 * t107 + t42 * t61 - t176 * t124 + (-t18 * t157 - t19 * t159 + (t157 * t35 - t159 * t36) * qJD(4)) * t128, t18 * t43 + t19 * t42 + t36 * t23 + t35 * t24 + t85 * t67 - t222, -t169 * t26 - t21 * t80, -t169 * t27 + t21 * t79 - t80 * t22 + t26 * t49, -t26 * t115 + t80 * t116 + t125 * t169 + t170 * t21, t22 * t79 + t49 * t27, -t27 * t115 - t79 * t116 - t49 * t125 + t170 * t22, t115 * t125 - t92, t4 * t115 + t16 * t116 + t7 * t125 - t170 * t2 + t68 * t22 + t44 * t27 + t32 * t79 + t39 * t49, -t3 * t115 - t17 * t116 - t8 * t125 - t162 * t170 + t169 * t39 - t68 * t21 - t44 * t26 + t32 * t80, t16 * t21 + t162 * t79 - t169 * t4 - t17 * t22 - t2 * t80 + t7 * t26 - t8 * t27 - t3 * t49, t2 * t16 - t162 * t17 + t8 * t3 + t32 * t68 + t44 * t39 + t7 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, -qJ(2) * t183, 0, 0, 0, 0, 0, 0, 0.2e1 * t122 * qJD(3), t140 + (-t120 - t188) * qJD(3), -t118 - t239, t91 * t120 + t122 * t247, 0, 0, 0, 0, 0, 0, t173 - t215, -t117 ^ 2 * t159 - t209 - t216, (t61 - t218) * t159 + t248 + t231, -t85 * t122 + t175 * t157 + t249 * t159, 0, 0, 0, 0, 0, 0, t179 - t228, -t223 - t241, t180 + t240, -t44 * t122 - t132 * t162 + t168 * t2 - t229 * t8 + t230 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, -t118 + t239, t140 + (t120 - t188) * qJD(3), -t214, 0, 0, -(qJD(2) + t135) * t122, t135 * t120 - t170 * t198, 0, 0, t107 * t181 - t219, (-t61 - t218) * t159 - t248 + t231, t117 * t181 + t209 - t216, t105 * t182 - t177, t173 + t215, -t117 * t122, -pkin(3) * t195 - t220 - t35 * t122 - t91 * t105 + (-pkin(7) * t201 - t40) * t117 + t165 * t157, pkin(3) * t61 - t91 * t107 + t36 * t122 + t221 + (pkin(7) * t202 + t41) * t117 + t165 * t159, t41 * t105 + t40 * t107 + ((t107 * qJD(4) - t195) * pkin(7) + t175) * t159 + ((-t61 + t199) * pkin(7) - t249) * t157, -t54 * pkin(3) - t35 * t40 - t36 * t41 - t85 * t91 + (-qJD(4) * t176 - t19 * t157 + t18 * t159) * pkin(7), -t21 * t132 - t169 * t229, t180 - t240, -t223 + t241, -t168 * t22 - t230 * t49, t179 + t228, -t115 * t122, t103 * t116 - t233 * t115 - t7 * t122 + t150 * t22 - t168 * t32 + t178 * t49 - t230 * t44, -t104 * t116 + t232 * t115 + t8 * t122 + t32 * t132 - t150 * t21 + t169 * t178 - t229 * t44, t103 * t21 - t104 * t22 - t2 * t132 - t162 * t168 + t169 * t233 + t229 * t7 + t230 * t8 + t232 * t49, t2 * t103 - t104 * t162 + t32 * t150 + t178 * t44 - t232 * t8 - t233 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, -t105 ^ 2 + t107 ^ 2, t105 * t117 - t61, -t217, t107 * t117 - t195, t116, -t85 * t107 + t249, t85 * t105 - t175, 0, 0, t235, t253, t252, -t235, t242, t116, -t10 * t115 + (-t107 * t49 - t115 * t200 + t236 * t116) * pkin(4) + t243, t11 * t115 + (-t107 * t169 - t115 * t186 - t116 * t156) * pkin(4) + t251, t10 * t169 + t11 * t49 + t8 * t169 - t7 * t49 + (t236 * t21 - t156 * t22 + (t156 * t169 - t236 * t49) * qJD(5)) * pkin(4), -t7 * t10 - t8 * t11 + (t236 * t2 - t162 * t156 - t107 * t44 + (-t156 * t7 + t236 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, t253, t252, -t235, t242, t116, t8 * t115 + t243, t7 * t115 + t251, 0, 0;];
tauc_reg = t1;
