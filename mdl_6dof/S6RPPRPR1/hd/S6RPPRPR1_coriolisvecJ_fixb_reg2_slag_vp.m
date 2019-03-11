% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:03
% EndTime: 2019-03-09 01:40:11
% DurationCPUTime: 3.82s
% Computational Cost: add. (6476->339), mult. (15783->450), div. (0->0), fcn. (11964->10), ass. (0->180)
t161 = sin(pkin(10));
t164 = cos(pkin(10));
t167 = sin(qJ(4));
t243 = cos(qJ(4));
t140 = t243 * t161 + t167 * t164;
t128 = t140 * qJD(1);
t160 = sin(pkin(11));
t163 = cos(pkin(11));
t107 = qJD(4) * t160 + t128 * t163;
t166 = sin(qJ(6));
t168 = cos(qJ(6));
t210 = t160 * t128;
t262 = t163 * qJD(4) - t210;
t177 = t168 * t262;
t64 = -t107 * t166 + t177;
t266 = t64 ^ 2;
t199 = t243 * t164;
t147 = qJD(1) * t199;
t202 = qJD(1) * t161;
t198 = t167 * t202;
t126 = -t147 + t198;
t123 = qJD(6) + t126;
t265 = t123 * t64;
t63 = t168 * t107 + t166 * t262;
t264 = t63 ^ 2;
t146 = qJD(4) * t147;
t208 = t161 * t167;
t197 = qJD(4) * t208;
t119 = qJD(1) * t197 - t146;
t207 = t163 * t119;
t263 = t126 * t262 - t207;
t206 = t168 * t163;
t251 = -t160 * t166 + t206;
t252 = t251 * qJD(6);
t226 = -t126 * t251 - t252;
t139 = t168 * t160 + t166 * t163;
t131 = t139 * qJD(6);
t225 = t139 * t126 + t131;
t211 = t160 * t119;
t259 = t107 * t126 - t211;
t133 = t140 * qJD(4);
t120 = qJD(1) * t133;
t175 = t199 - t208;
t96 = t120 * t175;
t258 = t126 * t133 - t96;
t257 = t119 * t139;
t256 = t128 * t262;
t194 = qJD(4) * t243;
t132 = -t164 * t194 + t197;
t233 = -t140 * t120 + t132 * t126;
t255 = t163 * t233;
t254 = t262 * t133;
t150 = sin(pkin(9)) * pkin(1) + qJ(3);
t239 = pkin(7) + t150;
t134 = t239 * t161;
t135 = t239 * t164;
t253 = -t243 * t134 - t135 * t167;
t250 = -qJD(6) + t123;
t31 = t166 * (qJD(6) * t107 - t211) - qJD(6) * t177 + t119 * t206;
t249 = -t225 * t63 - t251 * t31;
t248 = t120 * t139 - t226 * t123;
t125 = t126 ^ 2;
t247 = -t120 * t160 - t125 * t163;
t246 = t128 ^ 2;
t238 = pkin(8) + qJ(5);
t142 = t238 * t160;
t143 = t238 * t163;
t105 = -t142 * t166 + t143 * t168;
t242 = pkin(8) * t163;
t144 = t150 * qJD(1);
t118 = t161 * qJD(2) + t164 * t144;
t232 = pkin(7) * qJD(1);
t110 = t164 * t232 + t118;
t103 = t167 * t110;
t155 = t164 * qJD(2);
t109 = t155 + (-t144 - t232) * t161;
t65 = t243 * t109 - t103;
t92 = pkin(4) * t128 + qJ(5) * t126;
t36 = -t160 * t65 + t163 * t92;
t21 = pkin(5) * t128 + t126 * t242 + t36;
t217 = t126 * t160;
t37 = t160 * t92 + t163 * t65;
t25 = pkin(8) * t217 + t37;
t245 = qJD(5) * t139 + qJD(6) * t105 - t166 * t25 + t168 * t21;
t104 = -t142 * t168 - t143 * t166;
t244 = qJD(5) * t251 + qJD(6) * t104 - t166 * t21 - t168 * t25;
t172 = t140 * qJD(3);
t66 = t167 * t109 + t243 * t110;
t46 = qJD(1) * t172 + t66 * qJD(4);
t241 = t46 * t253;
t240 = t63 * t64;
t32 = qJD(6) * t63 - t257;
t39 = t131 * t140 + t251 * t132;
t89 = t251 * t140;
t237 = -t89 * t32 - t39 * t64;
t40 = -t132 * t139 + t252 * t140;
t88 = t139 * t140;
t236 = -t88 * t120 - t40 * t123;
t193 = qJD(3) * t202;
t205 = qJD(3) * t147 + t109 * t194;
t174 = -t167 * t193 + t205;
t44 = (qJD(5) - t103) * qJD(4) + t174;
t57 = pkin(4) * t120 + qJ(5) * t119 - qJD(5) * t128;
t18 = t160 * t57 + t163 * t44;
t235 = t63 * t133 + t175 * t31;
t56 = qJD(4) * qJ(5) + t66;
t141 = -cos(pkin(9)) * pkin(1) - pkin(3) * t164 - pkin(2);
t124 = qJD(1) * t141 + qJD(3);
t74 = pkin(4) * t126 - qJ(5) * t128 + t124;
t27 = t160 * t74 + t163 * t56;
t67 = t175 * qJD(3) + t253 * qJD(4);
t75 = pkin(4) * t133 + qJ(5) * t132 - qJD(5) * t140;
t30 = t160 * t75 + t163 * t67;
t234 = t107 * t133 + t175 * t207;
t85 = -pkin(4) * t175 - qJ(5) * t140 + t141;
t91 = -t167 * t134 + t243 * t135;
t42 = t160 * t85 + t163 * t91;
t231 = t128 * t64;
t230 = t128 * t63;
t228 = t175 * t46;
t227 = t163 * t120 - t160 * t125;
t224 = t107 * t128;
t223 = t107 * t160;
t222 = t119 * t140;
t219 = t126 * t128;
t216 = t132 * t160;
t215 = t132 * t163;
t213 = t175 * t119;
t212 = t140 * t160;
t203 = t161 ^ 2 + t164 ^ 2;
t201 = t132 * qJD(4);
t17 = -t160 * t44 + t163 * t57;
t26 = -t160 * t56 + t163 * t74;
t29 = -t160 * t67 + t163 * t75;
t41 = -t160 * t91 + t163 * t85;
t192 = qJD(1) * t203;
t190 = -t139 * t32 - t226 * t64;
t189 = t120 * t251 - t225 * t123;
t188 = -t31 * t88 + t40 * t63;
t187 = t119 * t253 + t46 * t140;
t12 = pkin(5) * t120 + pkin(8) * t207 + t17;
t13 = pkin(8) * t211 + t18;
t186 = t12 * t166 + t13 * t168;
t185 = -t120 * t89 + t123 * t39;
t184 = t133 * t64 + t175 * t32;
t14 = pkin(5) * t126 - pkin(8) * t107 + t26;
t19 = pkin(8) * t262 + t27;
t5 = t14 * t168 - t166 * t19;
t6 = t14 * t166 + t168 * t19;
t183 = -t17 * t160 + t18 * t163;
t182 = t160 * t26 - t163 * t27;
t28 = -pkin(5) * t175 - t140 * t242 + t41;
t34 = -pkin(8) * t212 + t42;
t9 = -t166 * t34 + t168 * t28;
t10 = t166 * t28 + t168 * t34;
t181 = (-t144 * t161 + t155) * t161 - t118 * t164;
t180 = t128 * t133 + t213;
t178 = t163 * t262;
t54 = -qJD(4) * pkin(4) + qJD(5) - t65;
t173 = -t54 * t132 + t187;
t171 = pkin(4) * t119 - qJ(5) * t120 + (-qJD(5) + t54) * t126;
t2 = -qJD(6) * t6 + t168 * t12 - t13 * t166;
t68 = qJD(4) * t91 + t172;
t35 = -pkin(5) * t211 + t46;
t158 = t163 ^ 2;
t156 = t160 ^ 2;
t153 = -t163 * pkin(5) - pkin(4);
t122 = t133 * qJD(4);
t80 = t132 * t178;
t71 = pkin(5) * t212 - t253;
t48 = -pkin(5) * t216 + t68;
t47 = -pkin(5) * t217 + t66;
t45 = (-qJD(4) * t110 - t193) * t167 + t205;
t38 = -pkin(5) * t262 + t54;
t22 = pkin(8) * t216 + t30;
t15 = pkin(5) * t133 + pkin(8) * t215 + t29;
t4 = -qJD(6) * t10 + t15 * t168 - t166 * t22;
t3 = qJD(6) * t9 + t15 * t166 + t168 * t22;
t1 = qJD(6) * t5 + t186;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t192 (t150 * t192 - t181) * qJD(3), -t128 * t132 - t222, -t180 + t233, -t201, t258, -t122, 0, -qJD(4) * t68 + t120 * t141 + t124 * t133, -qJD(4) * t67 - t119 * t141 - t124 * t132, -t120 * t91 - t126 * t67 + t128 * t68 + t132 * t65 - t133 * t66 + t175 * t45 + t187, t45 * t91 - t65 * t68 + t66 * t67 - t241, -t107 * t215 - t158 * t222, -t80 + (t107 * t132 + 0.2e1 * t140 * t207) * t160, t234 - t255, -t156 * t222 + t216 * t262, t254 + (t233 - t213) * t160, t258, t41 * t120 + t29 * t126 + t26 * t133 + t160 * t173 - t17 * t175 - t262 * t68, t107 * t68 - t120 * t42 - t126 * t30 - t133 * t27 + t163 * t173 + t175 * t18, -t29 * t107 - t30 * t210 + (t30 * qJD(4) + t41 * t119 + t26 * t132 - t17 * t140) * t163 + (t42 * t119 + t27 * t132 - t18 * t140) * t160, t17 * t41 + t18 * t42 + t26 * t29 + t27 * t30 + t54 * t68 - t241, -t31 * t89 - t39 * t63, -t188 + t237, -t185 + t235, t32 * t88 - t40 * t64, t184 + t236, t123 * t133 - t96, t120 * t9 + t123 * t4 + t133 * t5 - t175 * t2 + t32 * t71 + t35 * t88 + t38 * t40 - t48 * t64, t1 * t175 - t10 * t120 - t123 * t3 - t133 * t6 - t31 * t71 + t35 * t89 - t38 * t39 + t48 * t63, -t1 * t88 - t10 * t32 - t2 * t89 + t3 * t64 + t31 * t9 + t39 * t5 - t4 * t63 - t40 * t6, t1 * t10 + t2 * t9 + t3 * t6 + t35 * t71 + t38 * t48 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t201, t180 + t233, -t132 * t66 - t133 * t65 + t140 * t45 - t228, 0, 0, 0, 0, 0, 0, -t254 + (t233 + t213) * t160, t234 + t255, -t107 * t216 - t80, t132 * t182 + t133 * t54 + t140 * t183 - t228, 0, 0, 0, 0, 0, 0, -t184 + t236, t185 + t235, t188 + t237, t1 * t89 + t133 * t38 - t175 * t35 - t2 * t88 - t39 * t6 - t40 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203 * qJD(1) ^ 2, t181 * qJD(1), 0, 0, 0, 0, 0, 0, 0.2e1 * t128 * qJD(4), t146 + (-t126 - t198) * qJD(4), -t125 - t246, t126 * t66 + t128 * t65, 0, 0, 0, 0, 0, 0, t227 + t256, -t224 + t247 (t178 + t223) * t126 + (t156 + t158) * t119, -t126 * t182 - t128 * t54 + t18 * t160 + t17 * t163, 0, 0, 0, 0, 0, 0, t189 + t231, -t230 - t248, t190 - t249, t1 * t139 - t128 * t38 + t2 * t251 - t225 * t5 - t226 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, -t125 + t246, t146 + (t126 - t198) * qJD(4), -t219, 0, 0 -(qJD(3) + t124) * t128, t124 * t126 + (t65 + t103) * qJD(4) - t174, 0, 0, t259 * t163 (t178 - t223) * t126 + (t156 - t158) * t119, -t224 - t247, -t263 * t160, t227 - t256, -t219, -t36 * t126 - t26 * t128 + t160 * t171 - t46 * t163 + t262 * t66, -t107 * t66 + t126 * t37 + t128 * t27 + t160 * t46 + t163 * t171, t36 * t107 + t37 * t210 + (-qJD(5) * t210 - t26 * t126 + t18 + (t163 * qJD(5) - t37) * qJD(4)) * t163 + (qJD(5) * t107 - t27 * t126 - t17) * t160, -pkin(4) * t46 + qJ(5) * t183 - qJD(5) * t182 - t26 * t36 - t27 * t37 - t54 * t66, -t139 * t31 - t226 * t63, t190 + t249, -t230 + t248, -t225 * t64 - t251 * t32, t189 - t231, -t123 * t128, t104 * t120 - t245 * t123 - t128 * t5 + t153 * t32 + t225 * t38 - t251 * t35 + t47 * t64, -t105 * t120 - t244 * t123 + t128 * t6 + t139 * t35 - t153 * t31 - t226 * t38 - t47 * t63, t1 * t251 + t104 * t31 - t105 * t32 - t139 * t2 - t225 * t6 + t226 * t5 + t244 * t64 + t245 * t63, t1 * t105 + t104 * t2 + t153 * t35 + t244 * t6 - t245 * t5 - t38 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, t263, -t107 ^ 2 - t262 ^ 2, t107 * t26 - t262 * t27 + t46, 0, 0, 0, 0, 0, 0, t63 * t123 + t32, -t31 + t265, -t264 - t266, t5 * t63 - t6 * t64 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, t264 - t266, -t31 - t265, t240, t250 * t63 + t257, t120, t123 * t6 - t38 * t63 + t2, t250 * t5 - t38 * t64 - t186, 0, 0;];
tauc_reg  = t7;
