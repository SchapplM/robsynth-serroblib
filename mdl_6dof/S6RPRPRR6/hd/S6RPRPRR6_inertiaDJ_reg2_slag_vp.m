% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:53
% EndTime: 2019-03-09 03:53:06
% DurationCPUTime: 4.52s
% Computational Cost: add. (8481->296), mult. (18785->551), div. (0->0), fcn. (19706->10), ass. (0->142)
t101 = sin(pkin(11));
t103 = cos(pkin(11));
t199 = sin(qJ(5));
t202 = cos(qJ(5));
t89 = t101 * t202 + t103 * t199;
t102 = sin(pkin(10));
t104 = cos(pkin(10));
t200 = sin(qJ(3));
t203 = cos(qJ(3));
t158 = t102 * t203 + t104 * t200;
t151 = t103 * t158;
t157 = -t102 * t200 + t104 * t203;
t96 = -t104 * pkin(2) - pkin(1);
t128 = -pkin(3) * t157 - qJ(4) * t158 + t96;
t195 = pkin(7) + qJ(2);
t206 = t195 * t203;
t207 = t195 * t200;
t132 = -t102 * t207 + t104 * t206;
t43 = -t101 * t132 + t103 * t128;
t115 = -pkin(4) * t157 - pkin(8) * t151 + t43;
t152 = t101 * t158;
t44 = t101 * t128 + t103 * t132;
t116 = -pkin(8) * t152 + t44;
t17 = t115 * t202 - t116 * t199;
t144 = t158 * t199;
t149 = t157 * qJD(3);
t212 = qJD(5) * t144 - t149 * t202;
t145 = t158 * t202;
t211 = qJD(5) * t145 + t149 * t199;
t194 = pkin(8) + qJ(4);
t210 = t89 * t194;
t209 = qJD(2) * t200 + qJD(3) * t206;
t208 = -qJD(2) * t203 + qJD(3) * t207;
t18 = t115 * t199 + t116 * t202;
t65 = t102 * t206 + t104 * t207;
t119 = -t102 * t209 - t104 * t208;
t205 = 0.2e1 * t96;
t83 = t158 * qJD(3);
t204 = t83 * pkin(5);
t201 = cos(qJ(6));
t105 = sin(qJ(6));
t182 = qJD(5) * t202;
t181 = qJD(5) * t199;
t91 = t101 * t181;
t160 = -t103 * t182 + t91;
t190 = t101 * t182 + t103 * t181;
t156 = t101 * t199 - t103 * t202;
t60 = -t105 * t156 + t201 * t89;
t42 = qJD(6) * t60 - t105 * t160 + t190 * t201;
t143 = t201 * t156;
t59 = t105 * t89 + t143;
t198 = t59 * t42;
t187 = qJD(6) * t105;
t41 = qJD(6) * t143 + t105 * t190 + t160 * t201 + t187 * t89;
t197 = t60 * t41;
t53 = -t102 * t208 + t104 * t209;
t196 = t65 * t53;
t127 = -t101 * t212 + t103 * t211;
t37 = t101 * t211 + t103 * t212;
t120 = t101 * t145 + t103 * t144;
t58 = -t101 * t144 + t103 * t145;
t39 = -t105 * t120 + t201 * t58;
t12 = qJD(6) * t39 - t105 * t37 + t127 * t201;
t118 = t201 * t120;
t38 = t105 * t58 + t118;
t193 = -t12 * t60 + t38 * t41;
t191 = t120 * t160 - t127 * t89;
t167 = t202 * t194;
t168 = t194 * t199;
t66 = -t101 * t168 + t103 * t167;
t189 = t101 * t83;
t75 = t103 * t83;
t188 = qJD(4) * t157;
t61 = -0.2e1 * t157 * t83;
t186 = pkin(5) * t187;
t185 = t190 * pkin(5);
t95 = -t103 * pkin(4) - pkin(3);
t180 = qJD(6) * t201;
t178 = t202 * qJD(4);
t176 = t199 * qJD(4);
t175 = pkin(5) * t180;
t97 = t101 ^ 2;
t99 = t103 ^ 2;
t174 = 0.2e1 * (t97 + t99) * qJD(4);
t173 = 0.2e1 * (t102 ^ 2 + t104 ^ 2) * qJD(2);
t11 = qJD(6) * t118 + t105 * t127 + t187 * t58 + t201 * t37;
t166 = -t11 * t59 + t39 * t42;
t165 = -t157 * t41 - t60 * t83;
t133 = t83 * pkin(3) - qJD(4) * t158;
t150 = qJD(2) * t157;
t153 = qJ(4) * t157;
t28 = t103 * t133 - t101 * t150 + (t101 * t65 - t103 * t153) * qJD(3);
t29 = t101 * t133 + t103 * t150 + (-t101 * t153 - t103 * t65) * qJD(3);
t164 = -t101 * t28 + t103 * t29;
t121 = -qJ(4) * t149 + t133;
t141 = t103 * t149;
t111 = t83 * pkin(4) - pkin(8) * t141 - t101 * t119 + t103 * t121;
t142 = t101 * t149;
t112 = -pkin(8) * t142 + t101 * t121 + t103 * t119;
t9 = -qJD(5) * t18 + t202 * t111 - t199 * t112;
t106 = t37 * pkin(9) + t204 + t9;
t109 = -pkin(5) * t157 - t58 * pkin(9) + t17;
t107 = t201 * t109;
t8 = -qJD(5) * t17 - t199 * t111 - t202 * t112;
t110 = -pkin(9) * t127 - t8;
t14 = -pkin(9) * t120 + t18;
t1 = -qJD(6) * t107 - t105 * t106 - t110 * t201 + t14 * t187;
t159 = t89 * t160;
t54 = pkin(4) * t152 + t65;
t155 = t97 * t157;
t154 = t99 * t157;
t148 = t157 * t160 + t83 * t89;
t140 = t156 * t190;
t139 = t101 * t141;
t134 = (qJD(5) * t168 - t178) * t101;
t130 = -t156 * t37 + t190 * t58;
t51 = qJD(5) * t210 + t101 * t176 - t103 * t178;
t129 = -t89 * pkin(9) - t210;
t126 = t105 * t129;
t125 = t201 * t129;
t124 = t158 * t149;
t123 = -pkin(9) * t190 - t51;
t122 = 0.2e1 * t124;
t117 = t91 * pkin(9) + (-t176 + (-pkin(9) * t202 - t167) * qJD(5)) * t103 + t134;
t108 = t105 * t109;
t2 = -qJD(6) * t108 - t105 * t110 + t106 * t201 - t14 * t180;
t74 = pkin(4) * t142;
t70 = pkin(5) * t156 + t95;
t57 = -pkin(9) * t156 + t66;
t52 = (-qJD(5) * t167 - t176) * t103 + t134;
t48 = -t156 * t83 + t157 * t190;
t47 = t74 + t53;
t40 = pkin(5) * t120 + t54;
t33 = t201 * t57 + t126;
t32 = -t105 * t57 + t125;
t22 = t157 * t42 - t59 * t83;
t20 = pkin(5) * t127 + qJD(2) * t158 + qJD(3) * t132 + t74;
t16 = -qJD(6) * t126 - t105 * t123 + t117 * t201 - t180 * t57;
t15 = -qJD(6) * t125 - t105 * t117 - t123 * t201 + t187 * t57;
t7 = t14 * t201 + t108;
t6 = -t105 * t14 + t107;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, qJ(2) * t173, t122, 0.2e1 * t149 * t157 - 0.2e1 * t158 * t83, 0, t61, 0, 0, t83 * t205, t149 * t205, 0.2e1 * t119 * t157 - 0.2e1 * t132 * t83 + 0.2e1 * t149 * t65 + 0.2e1 * t158 * t53, 0.2e1 * t119 * t132 + 0.2e1 * t196, t99 * t122, -0.4e1 * t103 * t101 * t124, -0.2e1 * t141 * t157 + 0.2e1 * t151 * t83, t97 * t122, 0.2e1 * t142 * t157 - 0.2e1 * t152 * t83, t61, 0.2e1 * t142 * t65 + 0.2e1 * t152 * t53 - 0.2e1 * t157 * t28 + 0.2e1 * t43 * t83, 0.2e1 * t141 * t65 + 0.2e1 * t151 * t53 + 0.2e1 * t157 * t29 - 0.2e1 * t44 * t83, -0.2e1 * t141 * t43 - 0.2e1 * t142 * t44 - 0.2e1 * t151 * t28 - 0.2e1 * t152 * t29, 0.2e1 * t28 * t43 + 0.2e1 * t29 * t44 + 0.2e1 * t196, -0.2e1 * t58 * t37, 0.2e1 * t120 * t37 - 0.2e1 * t127 * t58, 0.2e1 * t157 * t37 + 0.2e1 * t58 * t83, 0.2e1 * t120 * t127, -0.2e1 * t120 * t83 + 0.2e1 * t127 * t157, t61, 0.2e1 * t120 * t47 + 0.2e1 * t127 * t54 - 0.2e1 * t157 * t9 + 0.2e1 * t17 * t83, -0.2e1 * t157 * t8 - 0.2e1 * t18 * t83 - 0.2e1 * t37 * t54 + 0.2e1 * t47 * t58, 0.2e1 * t120 * t8 - 0.2e1 * t127 * t18 + 0.2e1 * t17 * t37 - 0.2e1 * t9 * t58, 0.2e1 * t17 * t9 - 0.2e1 * t18 * t8 + 0.2e1 * t47 * t54, -0.2e1 * t39 * t11, 0.2e1 * t11 * t38 - 0.2e1 * t12 * t39, 0.2e1 * t11 * t157 + 0.2e1 * t39 * t83, 0.2e1 * t38 * t12, 0.2e1 * t12 * t157 - 0.2e1 * t38 * t83, t61, 0.2e1 * t12 * t40 - 0.2e1 * t157 * t2 + 0.2e1 * t20 * t38 + 0.2e1 * t6 * t83, -0.2e1 * t1 * t157 - 0.2e1 * t11 * t40 + 0.2e1 * t20 * t39 - 0.2e1 * t7 * t83, 0.2e1 * t1 * t38 + 0.2e1 * t11 * t6 - 0.2e1 * t12 * t7 - 0.2e1 * t2 * t39, -0.2e1 * t1 * t7 + 0.2e1 * t2 * t6 + 0.2e1 * t20 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t149, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t189 (-t155 - t154) * qJD(3), t101 * t29 + t103 * t28, 0, 0, 0, 0, 0, 0, t48, -t148, t130 + t191, -t156 * t9 - t160 * t18 - t17 * t190 - t8 * t89, 0, 0, 0, 0, 0, 0, t22, t165, t166 + t193, -t1 * t60 - t2 * t59 - t41 * t7 - t42 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t159 + 0.2e1 * t140, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t197 + 0.2e1 * t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, -t83, 0, -t53, -t119, 0, 0, t139 (-t155 + t154) * qJD(3), t189, -t139, t75, 0, -pkin(3) * t142 - qJ(4) * t189 + t101 * t188 - t53 * t103, -pkin(3) * t141 - qJ(4) * t75 + t53 * t101 + t103 * t188, t164, -pkin(3) * t53 + (-t101 * t43 + t103 * t44) * qJD(4) + t164 * qJ(4), -t160 * t58 - t37 * t89, -t130 + t191, t148, t120 * t190 + t127 * t156, t48, 0, t127 * t95 + t156 * t47 - t157 * t52 + t190 * t54 - t210 * t83, -t157 * t51 - t160 * t54 - t37 * t95 + t47 * t89 - t66 * t83, t120 * t51 - t127 * t66 + t156 * t8 + t160 * t17 - t18 * t190 - t210 * t37 - t52 * t58 - t9 * t89, t17 * t52 - t18 * t51 - t210 * t9 + t47 * t95 - t66 * t8, -t11 * t60 - t39 * t41, -t166 + t193, -t165, t12 * t59 + t38 * t42, t22, 0, t12 * t70 - t157 * t16 + t185 * t38 + t20 * t59 + t32 * t83 + t40 * t42, -t11 * t70 - t15 * t157 + t185 * t39 + t20 * t60 - t33 * t83 - t40 * t41, t1 * t59 + t11 * t32 - t12 * t33 + t15 * t38 - t16 * t39 - t2 * t60 + t41 * t6 - t42 * t7, -t1 * t33 - t15 * t7 + t16 * t6 + t185 * t40 + t2 * t32 + t20 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156 * t52 - t160 * t66 + t190 * t210 - t89 * t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t60 - t16 * t59 - t32 * t42 - t33 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, qJ(4) * t174, -0.2e1 * t159, 0.2e1 * t156 * t160 - 0.2e1 * t190 * t89, 0, 0.2e1 * t140, 0, 0, 0.2e1 * t95 * t190, -0.2e1 * t95 * t160, 0.2e1 * t156 * t51 - 0.2e1 * t160 * t210 - 0.2e1 * t190 * t66 - 0.2e1 * t52 * t89, -0.2e1 * t210 * t52 - 0.2e1 * t51 * t66, -0.2e1 * t197, 0.2e1 * t41 * t59 - 0.2e1 * t42 * t60, 0, 0.2e1 * t198, 0, 0, 0.2e1 * t185 * t59 + 0.2e1 * t42 * t70, 0.2e1 * t185 * t60 - 0.2e1 * t41 * t70, 0.2e1 * t15 * t59 - 0.2e1 * t16 * t60 + 0.2e1 * t32 * t41 - 0.2e1 * t33 * t42, -0.2e1 * t15 * t33 + 0.2e1 * t16 * t32 + 0.2e1 * t185 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, t141, 0, t53, 0, 0, 0, 0, 0, 0, t127, -t37, 0, t47, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, -t160, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, 0, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, -t127, t83, t9, t8, 0, 0, 0, 0, -t11, 0, -t12, t83, t157 * t186 + t201 * t204 + t2 (-t105 * t83 + t157 * t180) * pkin(5) + t1 (t201 * t11 - t105 * t12 + (t105 * t39 - t201 * t38) * qJD(6)) * pkin(5) (t201 * t2 - t1 * t105 + (-t105 * t6 + t201 * t7) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, t160, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t41, 0 (-t201 * t42 - t105 * t41 + (t105 * t59 + t201 * t60) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, 0, -t190, 0, t52, t51, 0, 0, 0, 0, -t41, 0, -t42, 0, t16, t15 (t201 * t41 - t105 * t42 + (t105 * t60 - t201 * t59) * qJD(6)) * pkin(5) (t201 * t16 - t105 * t15 + (-t105 * t32 + t201 * t33) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t186, -0.2e1 * t175, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, t83, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, -t42, 0, t16, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t186, -t175, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
