% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:49
% EndTime: 2019-03-09 01:53:57
% DurationCPUTime: 3.66s
% Computational Cost: add. (6561->360), mult. (14770->481), div. (0->0), fcn. (10903->8), ass. (0->187)
t160 = sin(qJ(6));
t161 = cos(qJ(6));
t156 = sin(pkin(9));
t239 = sin(qJ(4));
t198 = t239 * t156;
t186 = qJD(1) * t198;
t158 = cos(pkin(9));
t162 = cos(qJ(4));
t207 = t158 * t162;
t197 = qJD(1) * t207;
t119 = -t186 + t197;
t155 = sin(pkin(10));
t107 = t155 * t119;
t157 = cos(pkin(10));
t251 = qJD(4) * t157 - t107;
t176 = t161 * t251;
t98 = t155 * qJD(4) + t119 * t157;
t55 = -t160 * t98 + t176;
t260 = t55 ^ 2;
t127 = t162 * t156 + t239 * t158;
t170 = qJD(1) * t127;
t113 = qJD(6) + t170;
t259 = t113 * t55;
t54 = t160 * t251 + t161 * t98;
t258 = t54 ^ 2;
t129 = -t198 + t207;
t206 = t161 * t157;
t249 = -t155 * t160 + t206;
t250 = t249 * qJD(6);
t222 = -t170 * t249 - t250;
t128 = t161 * t155 + t160 * t157;
t122 = t128 * qJD(6);
t221 = t128 * t170 + t122;
t256 = t129 * qJD(3);
t123 = t127 * qJD(4);
t108 = qJD(1) * t123;
t210 = t155 * t108;
t255 = t170 * t98 - t210;
t208 = t157 * t108;
t254 = t170 * t251 - t208;
t192 = qJD(4) * t239;
t202 = qJD(4) * t162;
t196 = t158 * t202;
t124 = -t156 * t192 + t196;
t136 = qJD(4) * t186;
t109 = qJD(1) * t196 - t136;
t89 = t109 * t127;
t179 = t124 * t170 + t89;
t253 = t108 * t128;
t252 = t119 * t251;
t82 = t249 * t127;
t159 = -pkin(1) - qJ(3);
t235 = -pkin(7) + t159;
t130 = t235 * t156;
t131 = t235 * t158;
t91 = t239 * t130 - t162 * t131;
t248 = t159 * qJD(1);
t205 = t156 ^ 2 + t158 ^ 2;
t247 = t205 * qJD(3);
t246 = -qJD(6) + t113;
t22 = t160 * (qJD(6) * t98 - t210) - qJD(6) * t176 + t108 * t206;
t245 = -t22 * t249 - t221 * t54;
t244 = t109 * t128 - t222 * t113;
t114 = t170 ^ 2;
t243 = -t109 * t155 - t114 * t157;
t242 = t119 ^ 2;
t153 = qJD(1) * qJD(2);
t200 = 0.2e1 * t153;
t238 = pkin(8) * t157;
t137 = qJD(2) + t248;
t191 = -pkin(7) * qJD(1) + t137;
t110 = t191 * t156;
t105 = t239 * t110;
t111 = t191 * t158;
t74 = t162 * t111 - t105;
t85 = pkin(4) * t119 + qJ(5) * t170;
t32 = -t155 * t74 + t157 * t85;
t20 = pkin(5) * t119 + t170 * t238 + t32;
t217 = t170 * t155;
t33 = t155 * t85 + t157 * t74;
t26 = pkin(8) * t217 + t33;
t234 = pkin(8) + qJ(5);
t133 = t234 * t155;
t134 = t234 * t157;
t96 = -t133 * t160 + t134 * t161;
t241 = qJD(5) * t128 + qJD(6) * t96 - t160 * t26 + t161 * t20;
t95 = -t133 * t161 - t134 * t160;
t240 = qJD(5) * t249 + qJD(6) * t95 - t160 * t20 - t161 * t26;
t148 = t156 * pkin(3);
t175 = t256 * qJD(1);
t75 = t162 * t110 + t239 * t111;
t48 = t75 * qJD(4) + t175;
t237 = t48 * t91;
t236 = t54 * t55;
t102 = t111 * t202;
t171 = t127 * qJD(3);
t164 = -qJD(1) * t171 + t102;
t44 = (-t105 + qJD(5)) * qJD(4) + t164;
t50 = pkin(4) * t109 + qJ(5) * t108 - qJD(5) * t119 + t153;
t17 = t155 * t50 + t157 * t44;
t58 = pkin(4) * t124 + qJ(5) * t123 - qJD(5) * t129 + qJD(2);
t64 = -t91 * qJD(4) - t171;
t25 = t155 * t58 + t157 * t64;
t68 = qJD(4) * qJ(5) + t75;
t154 = qJD(1) * qJ(2);
t147 = qJD(3) + t154;
t132 = qJD(1) * t148 + t147;
t69 = pkin(4) * t170 - qJ(5) * t119 + t132;
t29 = t155 * t69 + t157 * t68;
t142 = qJ(2) + t148;
t84 = pkin(4) * t127 - qJ(5) * t129 + t142;
t92 = t162 * t130 + t239 * t131;
t40 = t155 * t84 + t157 * t92;
t233 = t119 * t55;
t232 = t119 * t54;
t231 = t119 * t98;
t229 = t48 * t129;
t228 = t91 * t108;
t227 = t98 * t155;
t225 = t157 * t109 - t155 * t114;
t224 = qJD(1) * t249 + qJD(6) * t82 + t124 * t128;
t80 = t128 * t127;
t223 = t128 * qJD(1) + qJD(6) * t80 - t124 * t249;
t218 = t170 * t119;
t216 = t119 * t123;
t215 = t123 * t155;
t214 = t123 * t157;
t212 = t129 * t108;
t211 = t129 * t155;
t204 = qJD(1) * t170;
t201 = t124 * qJD(4);
t16 = -t155 * t44 + t157 * t50;
t24 = -t155 * t64 + t157 * t58;
t28 = -t155 * t68 + t157 * t69;
t39 = -t155 * t92 + t157 * t84;
t190 = qJD(1) * t205;
t23 = qJD(6) * t54 - t253;
t188 = -t128 * t23 - t222 * t55;
t187 = t109 * t249 - t221 * t113;
t12 = pkin(5) * t109 + pkin(8) * t208 + t16;
t13 = pkin(8) * t210 + t17;
t185 = t12 * t160 + t13 * t161;
t63 = -qJD(4) * pkin(4) + qJD(5) - t74;
t184 = t63 * t123 - t229;
t183 = t28 * t124 + t16 * t127;
t182 = -t124 * t29 - t127 * t17;
t15 = pkin(5) * t170 - pkin(8) * t98 + t28;
t19 = pkin(8) * t251 + t29;
t5 = t15 * t161 - t160 * t19;
t6 = t15 * t160 + t161 * t19;
t181 = -t155 * t28 + t157 * t29;
t27 = pkin(5) * t127 - t129 * t238 + t39;
t30 = -pkin(8) * t211 + t40;
t9 = -t160 * t30 + t161 * t27;
t10 = t160 * t27 + t161 * t30;
t178 = t155 * t251;
t177 = t157 * t251;
t174 = -t184 - t228;
t173 = t108 * t127 - t129 * t109 + t123 * t170;
t172 = -t179 + t212;
t169 = t177 - t227;
t168 = t177 + t227;
t167 = pkin(4) * t108 - qJ(5) * t109 + (-qJD(5) + t63) * t170;
t47 = -t110 * t192 + t164;
t165 = -t123 * t74 + t124 * t75 + t127 * t47 - t229;
t2 = -qJD(6) * t6 + t161 * t12 - t13 * t160;
t65 = qJD(4) * t92 + t256;
t31 = -pkin(5) * t210 + t48;
t163 = qJD(1) ^ 2;
t151 = t157 ^ 2;
t149 = t155 ^ 2;
t145 = -t157 * pkin(5) - pkin(4);
t112 = t123 * qJD(4);
t83 = t249 * t129;
t81 = t128 * t129;
t66 = pkin(5) * t211 + t91;
t49 = -pkin(5) * t217 + t75;
t43 = -pkin(5) * t215 + t65;
t38 = -pkin(5) * t251 + t63;
t37 = -t123 * t128 + t250 * t129;
t35 = t122 * t129 + t249 * t123;
t18 = pkin(8) * t215 + t25;
t14 = pkin(5) * t124 + pkin(8) * t214 + t24;
t4 = -qJD(6) * t10 + t14 * t161 - t160 * t18;
t3 = qJD(6) * t9 + t14 * t160 + t161 * t18;
t1 = qJD(6) * t5 + t185;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, qJ(2) * t200, 0, 0, 0, 0, 0, 0, t156 * t200, t158 * t200, 0.2e1 * qJD(3) * t190 (t147 + t154) * qJD(2) + (-t137 - t248) * t247, -t212 - t216, -t119 * t124 + t173, -t112, t179, -t201, 0, 0.2e1 * t170 * qJD(2) - qJD(4) * t65 + t109 * t142 + t124 * t132, -qJD(4) * t64 - t108 * t142 - t123 * t132 + (qJD(1) * t129 + t119) * qJD(2), -t109 * t92 + t119 * t65 - t170 * t64 - t165 - t228, t47 * t92 + t237 + t64 * t75 - t65 * t74 + (qJD(1) * t142 + t132) * qJD(2), -t151 * t212 - t98 * t214, -t123 * t169 + 0.2e1 * t208 * t211, t124 * t98 - t157 * t173, t123 * t178 - t149 * t212, t124 * t251 + t155 * t173, t179, t39 * t109 + t155 * t174 + t170 * t24 - t251 * t65 + t183, -t109 * t40 + t157 * t174 - t170 * t25 + t65 * t98 + t182, -t25 * t107 - t24 * t98 + (t25 * qJD(4) + t39 * t108 + t28 * t123 - t16 * t129) * t157 + (t40 * t108 + t29 * t123 - t17 * t129) * t155, t16 * t39 + t17 * t40 + t24 * t28 + t25 * t29 + t63 * t65 + t237, -t22 * t83 - t35 * t54, t22 * t81 - t23 * t83 - t35 * t55 - t37 * t54, t109 * t83 - t113 * t35 + t124 * t54 - t127 * t22, t23 * t81 - t37 * t55, -t109 * t81 - t113 * t37 + t124 * t55 - t127 * t23, t113 * t124 + t89, t109 * t9 + t113 * t4 + t124 * t5 + t127 * t2 + t23 * t66 + t31 * t81 + t37 * t38 - t43 * t55, -t1 * t127 - t10 * t109 - t113 * t3 - t124 * t6 - t22 * t66 + t31 * t83 - t35 * t38 + t43 * t54, -t1 * t81 - t10 * t23 - t2 * t83 + t22 * t9 + t3 * t55 + t35 * t5 - t37 * t6 - t4 * t54, t1 * t10 + t2 * t9 + t3 * t6 + t31 * t66 + t38 * t43 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t163 * qJ(2), 0, 0, 0, 0, 0, 0, -t163 * t156, -t163 * t158, 0 (-t147 - t247) * qJD(1), 0, 0, 0, 0, 0, 0, -t112 - t204, -qJD(1) * t119 - t201, t172 + t216, -qJD(1) * t132 + t165, 0, 0, 0, 0, 0, 0, -t123 * t251 + t155 * t172 - t157 * t204, t98 * t123 + t155 * t204 + t157 * t172, t168 * t124 + (t157 * t98 - t178) * qJD(1) (-qJD(1) * t28 - t182) * t157 + (-qJD(1) * t29 - t183) * t155 + t184, 0, 0, 0, 0, 0, 0, -t109 * t80 - t224 * t113 - t123 * t55 - t129 * t23, -t109 * t82 + t223 * t113 + t123 * t54 + t129 * t22, -t22 * t80 - t223 * t55 + t224 * t54 - t23 * t82, t1 * t82 + t123 * t38 - t129 * t31 - t2 * t80 - t223 * t6 - t224 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205 * t163, t137 * t190 + t153, 0, 0, 0, 0, 0, 0, -t136 + (t119 + t197) * qJD(4), -0.2e1 * t170 * qJD(4), -t114 - t242, t119 * t74 + t170 * t75 + t153, 0, 0, 0, 0, 0, 0, t225 + t252, -t231 + t243, t168 * t170 + (t149 + t151) * t108, -t119 * t63 + t17 * t155 + t16 * t157 + t170 * t181, 0, 0, 0, 0, 0, 0, t187 + t233, -t232 - t244, t188 - t245, t1 * t128 - t119 * t38 + t2 * t249 - t221 * t5 - t222 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t114 + t242, 0, -t218, t136 + (t119 - t197) * qJD(4), 0, -t132 * t119 - t175, -t102 + (t105 + t74) * qJD(4) + (t132 + qJD(3)) * t170, 0, 0, t255 * t157, t169 * t170 + (t149 - t151) * t108, -t231 - t243, -t254 * t155, t225 - t252, -t218, -t28 * t119 + t155 * t167 - t48 * t157 - t170 * t32 + t251 * t75, t119 * t29 + t155 * t48 + t157 * t167 + t170 * t33 - t75 * t98, t33 * t107 + t32 * t98 + (-qJD(5) * t107 - t28 * t170 + t17 + (t157 * qJD(5) - t33) * qJD(4)) * t157 + (qJD(5) * t98 - t170 * t29 - t16) * t155, -pkin(4) * t48 - t28 * t32 - t29 * t33 - t63 * t75 + t181 * qJD(5) + (-t16 * t155 + t17 * t157) * qJ(5), -t128 * t22 - t222 * t54, t188 + t245, -t232 + t244, -t221 * t55 - t23 * t249, t187 - t233, -t113 * t119, t109 * t95 - t241 * t113 - t119 * t5 + t145 * t23 + t221 * t38 - t249 * t31 + t49 * t55, -t109 * t96 - t240 * t113 + t119 * t6 + t128 * t31 - t145 * t22 - t222 * t38 - t49 * t54, t1 * t249 - t128 * t2 + t22 * t95 - t221 * t6 + t222 * t5 - t23 * t96 + t240 * t55 + t241 * t54, t1 * t96 + t145 * t31 + t2 * t95 + t240 * t6 - t241 * t5 - t38 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255, t254, -t251 ^ 2 - t98 ^ 2, -t251 * t29 + t28 * t98 + t48, 0, 0, 0, 0, 0, 0, t54 * t113 + t23, -t22 + t259, -t258 - t260, t5 * t54 - t6 * t55 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, t258 - t260, -t22 - t259, t236, t246 * t54 + t253, t109, t113 * t6 - t38 * t54 + t2, t246 * t5 - t38 * t55 - t185, 0, 0;];
tauc_reg  = t7;
