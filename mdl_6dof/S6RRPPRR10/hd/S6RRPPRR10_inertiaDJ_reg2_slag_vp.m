% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:01
% EndTime: 2019-03-09 09:37:12
% DurationCPUTime: 3.70s
% Computational Cost: add. (5441->284), mult. (11339->514), div. (0->0), fcn. (10412->8), ass. (0->159)
t102 = sin(qJ(6));
t192 = cos(qJ(6));
t100 = sin(pkin(10));
t101 = cos(pkin(10));
t191 = sin(qJ(5));
t193 = cos(qJ(5));
t73 = t100 * t193 + t101 * t191;
t165 = t191 * t100;
t166 = t193 * t101;
t74 = -t165 + t166;
t141 = t102 * t74 + t192 * t73;
t140 = t102 * t73 - t192 * t74;
t142 = t74 * qJD(5);
t205 = t73 * qJD(5);
t21 = qJD(6) * t140 + t102 * t205 - t142 * t192;
t208 = t141 * t21;
t128 = t73 * t142;
t190 = t74 * t205;
t229 = 0.2e1 * t128 - 0.2e1 * t190;
t104 = cos(qJ(2));
t103 = sin(qJ(2));
t175 = t103 * qJD(3);
t144 = -qJD(4) * t104 - t175;
t130 = t100 * t144;
t194 = pkin(3) + pkin(7);
t139 = t100 * qJ(3) + t101 * t194;
t133 = pkin(4) + t139;
t189 = pkin(2) + qJ(4);
t173 = pkin(8) + t189;
t182 = t100 * t103;
t112 = -t130 + (t104 * t133 - t173 * t182) * qJD(2);
t129 = t101 * t144;
t168 = t100 * t194;
t183 = qJ(3) * t104;
t113 = t129 + (t104 * t168 + (t103 * t173 - t183) * t101) * qJD(2);
t202 = t104 * t173 + pkin(1);
t121 = t100 * t202 + t103 * t133;
t180 = t103 * qJ(3);
t123 = t103 * t168 + (-t180 - t202) * t101;
t20 = t121 * t191 + t123 * t193;
t11 = -qJD(5) * t20 + t193 * t112 - t191 * t113;
t92 = qJD(2) * t104;
t172 = pkin(5) * t92;
t178 = qJD(2) * t103;
t44 = -t104 * t142 + t73 * t178;
t105 = -t44 * pkin(9) + t11 + t172;
t19 = t121 * t193 - t123 * t191;
t64 = t73 * t104;
t108 = t103 * pkin(5) + t64 * pkin(9) + t19;
t106 = t192 * t108;
t10 = -qJD(5) * t19 - t191 * t112 - t193 * t113;
t219 = t104 * t205;
t164 = t100 * t178;
t79 = t191 * t164;
t124 = -t79 + t219;
t162 = t101 * t178;
t146 = t193 * t162;
t117 = t146 + t124;
t109 = pkin(9) * t117 - t10;
t82 = t104 * t166;
t143 = t104 * t165 - t82;
t15 = pkin(9) * t143 + t20;
t176 = qJD(6) * t102;
t1 = -qJD(6) * t106 - t102 * t105 - t109 * t192 + t15 * t176;
t107 = t102 * t108;
t158 = qJD(6) * t192;
t2 = -qJD(6) * t107 - t102 * t109 + t105 * t192 - t15 * t158;
t9 = t15 * t192 + t107;
t228 = -t1 * t141 - t140 * t2 - t21 * t9;
t147 = t173 * t191;
t148 = t173 * t193;
t53 = t100 * t147 - t101 * t148;
t122 = -t74 * pkin(9) + t53;
t120 = t102 * t122;
t54 = -t100 * t148 - t101 * t147;
t37 = -pkin(9) * t73 + t54;
t18 = t192 * t37 + t120;
t204 = qJD(4) * t193 - qJD(5) * t147;
t217 = qJD(4) * t191 + qJD(5) * t148;
t32 = t100 * t217 - t101 * t204;
t110 = pkin(9) * t205 + t32;
t31 = t204 * t100 + t101 * t217;
t111 = pkin(9) * t142 + t31;
t116 = t192 * t122;
t6 = -qJD(6) * t116 - t102 * t110 + t111 * t192 + t176 * t37;
t7 = -qJD(6) * t120 + t102 * t111 + t110 * t192 - t158 * t37;
t227 = t140 * t7 + t141 * t6 + t18 * t21;
t42 = t102 * t143 - t192 * t64;
t13 = qJD(6) * t42 + t102 * t44 - t117 * t192;
t127 = t192 * t143;
t41 = -t102 * t64 - t127;
t226 = t13 * t141 - t21 * t41;
t225 = t103 * t21 - t141 * t92;
t224 = -(t102 * t140 + t141 * t192) * qJD(6) + t102 * t21;
t12 = -qJD(6) * t127 - t102 * t117 - t176 * t64 - t192 * t44;
t223 = t12 * t140;
t201 = -t102 * t142 - t192 * t205;
t22 = t158 * t73 + t176 * t74 - t201;
t222 = t140 * t22;
t221 = t140 * t92;
t181 = t101 * t104;
t81 = t194 * t104;
t69 = pkin(4) * t181 + t81;
t220 = t205 * t69;
t23 = -qJD(6) * t141 + t201;
t214 = -t140 * t23 - t208;
t156 = t189 * t103;
t34 = -t130 + (-t100 * t156 + t104 * t139) * qJD(2);
t138 = -t101 * qJ(3) + t168;
t35 = t129 + (t101 * t156 + t104 * t138) * qJD(2);
t16 = t100 * t35 + t101 * t34;
t203 = qJD(2) * (t103 ^ 2 - t104 ^ 2);
t96 = t100 ^ 2;
t97 = t101 ^ 2;
t77 = (t96 + t97) * qJD(4);
t200 = t142 * t54 - t205 * t53 - t31 * t73 + t32 * t74;
t199 = -t10 * t73 + t11 * t74 + t142 * t20 - t19 * t205;
t198 = -t103 * t142 - t73 * t92;
t197 = -t117 * t73 - t142 * t143;
t157 = t189 * t104;
t177 = qJD(3) * t104;
t196 = qJD(2) * (-t157 - t180) - qJD(4) * t103 + t177;
t195 = 0.2e1 * qJD(3);
t76 = t194 * t178;
t184 = t76 * t100;
t88 = pkin(4) * t100 + qJ(3);
t179 = qJ(3) * qJD(3);
t174 = -0.2e1 * pkin(1) * qJD(2);
t171 = pkin(5) * t176;
t170 = pkin(7) * t178;
t169 = pkin(7) * t92;
t163 = t100 * t92;
t161 = t103 * t92;
t84 = -0.2e1 * t161;
t155 = pkin(5) * t158;
t154 = 0.2e1 * t203;
t153 = t100 * t162;
t151 = -pkin(1) - t157;
t150 = t205 * t64 + t44 * t74;
t149 = -pkin(2) * t104 - t180;
t134 = -t146 + t79;
t125 = qJD(2) * t149 + t177;
t118 = t134 - t219;
t85 = t101 * t92;
t83 = 0.2e1 * t161;
t78 = -pkin(1) + t149;
t75 = -0.2e1 * t203;
t66 = -t175 + (pkin(2) * t103 - t183) * qJD(2);
t59 = pkin(5) * t142 + qJD(3);
t58 = (-pkin(4) * t101 - t194) * t178;
t57 = pkin(5) * t73 + t88;
t52 = t101 * t151 + t103 * t138;
t51 = -t100 * t151 + t103 * t139;
t46 = -pkin(5) * t143 + t69;
t45 = -t103 * t205 + t74 * t92;
t28 = -t124 * pkin(5) + ((-pkin(5) * t193 - pkin(4)) * t101 - t194) * t178;
t17 = -t102 * t37 + t116;
t8 = -t102 * t15 + t106;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t75, 0, t84, 0, 0, t103 * t174, t104 * t174, 0, 0, 0, 0, 0, t83, t75, t84, 0, 0.2e1 * t104 * t66 - 0.2e1 * t178 * t78, -0.2e1 * t103 * t66 - 0.2e1 * t78 * t92, 0.2e1 * t78 * t66, t96 * t84, -0.4e1 * t104 * t153, t100 * t154, t97 * t84, t101 * t154, t83, -0.2e1 * t76 * t181 + 0.2e1 * t34 * t103 + 0.2e1 * (-t101 * t103 * t81 + t104 * t51) * qJD(2), 0.2e1 * t104 * t184 - 0.2e1 * t35 * t103 + 0.2e1 * (-t104 * t52 + t182 * t81) * qJD(2), 0.2e1 * (t100 * t34 - t101 * t35) * t104 + 0.2e1 * (-t100 * t51 + t101 * t52) * t178, 0.2e1 * t34 * t51 + 0.2e1 * t35 * t52 - 0.2e1 * t76 * t81, -0.2e1 * t64 * t44, -0.2e1 * t117 * t64 + 0.2e1 * t143 * t44, 0.2e1 * t103 * t44 - 0.2e1 * t64 * t92, 0.2e1 * t143 * t117, 0.2e1 * t103 * t117 + 0.2e1 * t143 * t92, t83, 0.2e1 * t11 * t103 + 0.2e1 * t58 * t82 + 0.2e1 * t69 * t134 + 0.2e1 * (t19 * qJD(2) - t165 * t58 - t220) * t104, 0.2e1 * t10 * t103 - 0.2e1 * t20 * t92 + 0.2e1 * t44 * t69 - 0.2e1 * t58 * t64, -0.2e1 * t10 * t143 + 0.2e1 * t11 * t64 + 0.2e1 * t117 * t20 - 0.2e1 * t19 * t44, -0.2e1 * t10 * t20 + 0.2e1 * t11 * t19 + 0.2e1 * t58 * t69, -0.2e1 * t42 * t12, 0.2e1 * t12 * t41 - 0.2e1 * t13 * t42, -0.2e1 * t103 * t12 + 0.2e1 * t42 * t92, 0.2e1 * t41 * t13, -0.2e1 * t103 * t13 - 0.2e1 * t41 * t92, t83, 0.2e1 * t103 * t2 + 0.2e1 * t13 * t46 + 0.2e1 * t28 * t41 + 0.2e1 * t8 * t92, 0.2e1 * t1 * t103 - 0.2e1 * t12 * t46 + 0.2e1 * t28 * t42 - 0.2e1 * t9 * t92, 0.2e1 * t1 * t41 + 0.2e1 * t12 * t8 - 0.2e1 * t13 * t9 - 0.2e1 * t2 * t42, -0.2e1 * t1 * t9 + 0.2e1 * t2 * t8 + 0.2e1 * t28 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, -t178, 0, -t169, t170, 0, 0, 0, -t92, t178, 0, 0, 0, t125, t169, -t170, t125 * pkin(7), t153 (-t96 + t97) * t178, t85, -t153, -t163, 0, t101 * t196 - t184, -t100 * t196 - t76 * t101, -t16, -t76 * qJ(3) + t81 * qJD(3) + (-t100 * t52 - t101 * t51) * qJD(4) - t16 * t189, t150, t117 * t74 + t142 * t64 - t143 * t205 - t44 * t73, t45, t197, t198, 0, -qJD(3) * t143 + t32 * t103 + t118 * t88 + t142 * t69 + t53 * t92 + t58 * t73, -qJD(3) * t64 + t103 * t31 + t44 * t88 - t54 * t92 + t58 * t74 - t220, t117 * t54 - t143 * t31 + t32 * t64 - t53 * t44 - t199, qJD(3) * t69 - t10 * t54 + t11 * t53 + t19 * t32 - t20 * t31 + t58 * t88, -t22 * t42 + t223, t12 * t141 + t13 * t140 + t21 * t42 + t22 * t41, -t103 * t22 - t221, t226, t225, 0, t103 * t7 + t13 * t57 + t141 * t28 + t17 * t92 - t21 * t46 + t41 * t59, t103 * t6 - t12 * t57 - t140 * t28 - t18 * t92 - t22 * t46 + t42 * t59, t12 * t17 - t13 * t18 + t22 * t8 + t41 * t6 - t42 * t7 - t228, -t1 * t18 + t17 * t2 + t28 * t57 + t46 * t59 - t6 * t9 + t7 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, 0.2e1 * t179, 0, 0, 0, 0, 0, 0, t100 * t195, t101 * t195, 0.2e1 * t77, 0.2e1 * t189 * t77 + 0.2e1 * t179, -0.2e1 * t190, -0.2e1 * t142 * t74 + 0.2e1 * t205 * t73, 0, 0.2e1 * t128, 0, 0, 0.2e1 * qJD(3) * t73 + 0.2e1 * t142 * t88, 0.2e1 * qJD(3) * t74 - 0.2e1 * t205 * t88, -0.2e1 * t200, 0.2e1 * qJD(3) * t88 - 0.2e1 * t31 * t54 + 0.2e1 * t32 * t53, 0.2e1 * t222, -0.2e1 * t140 * t21 + 0.2e1 * t141 * t22, 0, -0.2e1 * t208, 0, 0, 0.2e1 * t141 * t59 - 0.2e1 * t21 * t57, -0.2e1 * t140 * t59 - 0.2e1 * t22 * t57, 0.2e1 * t17 * t22 + 0.2e1 * t227, 0.2e1 * t17 * t7 - 0.2e1 * t18 * t6 + 0.2e1 * t57 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, t169, 0, 0, 0, 0, 0, 0, t85, -t163, 0, t16, 0, 0, 0, 0, 0, 0, t45, t198, -t150 - t197, t199, 0, 0, 0, 0, 0, 0, t103 * t23 - t221, t225, -t23 * t42 - t223 - t226, t23 * t8 + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, 0, 0, 0, 0, 0, 0, 0, 0, -t229, t200, 0, 0, 0, 0, 0, 0, 0, 0, t208 - t214 - t222, t17 * t23 - t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t164, 0, -t76, 0, 0, 0, 0, 0, 0, t118, t44, 0, t58, 0, 0, 0, 0, 0, 0, t13, -t12, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, 0, t142, -t205, 0, qJD(3), 0, 0, 0, 0, 0, 0, -t21, -t22, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t117, t92, t11, t10, 0, 0, 0, 0, -t12, 0, -t13, t92, -t103 * t171 + t172 * t192 + t2 (-t102 * t92 - t103 * t158) * pkin(5) + t1 (t192 * t12 - t102 * t13 + (t102 * t42 - t192 * t41) * qJD(6)) * pkin(5) (t192 * t2 - t1 * t102 + (-t102 * t8 + t192 * t9) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, 0, -t142, 0, t32, t31, 0, 0, 0, 0, -t22, 0, t21, 0, t7, t6 (t192 * t22 + t224) * pkin(5) (t192 * t7 - t102 * t6 + (-t102 * t17 + t18 * t192) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, -t142, 0, 0, 0, 0, 0, 0, 0, 0, t23, t21, 0 (t192 * t23 - t224) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t171, -0.2e1 * t155, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t13, t92, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, t21, 0, t7, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t21, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, -t155, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
