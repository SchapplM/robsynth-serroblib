% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:34
% EndTime: 2019-03-09 06:24:41
% DurationCPUTime: 3.20s
% Computational Cost: add. (3830->257), mult. (7612->396), div. (0->0), fcn. (6798->6), ass. (0->158)
t108 = sin(qJ(5));
t106 = t108 ^ 2;
t110 = cos(qJ(5));
t107 = t110 ^ 2;
t174 = t106 + t107;
t109 = sin(qJ(3));
t203 = sin(qJ(4));
t153 = t203 * qJD(4);
t156 = t203 * t109;
t111 = cos(qJ(4));
t112 = cos(qJ(3));
t177 = t111 * t112;
t211 = qJD(3) + qJD(4);
t51 = -qJD(3) * t156 - t109 * t153 + t211 * t177;
t222 = t174 * t51;
t172 = qJD(4) * t111;
t164 = pkin(3) * t172;
t221 = t174 * t164;
t81 = t111 * t109 + t203 * t112;
t82 = -t156 + t177;
t121 = (-t111 * t81 + t203 * t82) * qJD(4);
t50 = t211 * t81;
t220 = (t111 * t50 - t203 * t51 + t121) * pkin(3);
t212 = (t106 - t107) * qJD(5);
t102 = qJD(5) * t108;
t197 = t82 * t50;
t199 = t81 * t51;
t215 = 0.2e1 * t197 - 0.2e1 * t199;
t77 = t81 ^ 2;
t78 = t82 ^ 2;
t217 = (t78 + t77) * t102 + t110 * t215;
t216 = 0.2e1 * t82;
t170 = t112 * qJD(3);
t94 = pkin(3) * t170 + qJD(2);
t118 = t51 * pkin(4) + t50 * pkin(9) + t94;
t97 = t109 * pkin(3) + qJ(2);
t128 = t81 * pkin(4) - t82 * pkin(9) + t97;
t113 = -pkin(1) - pkin(7);
t196 = pkin(8) - t113;
t83 = t196 * t109;
t84 = t196 * t112;
t54 = -t111 * t83 - t203 * t84;
t193 = t108 * t128 + t110 * t54;
t171 = t109 * qJD(3);
t145 = t196 * t171;
t76 = qJD(3) * t84;
t26 = t111 * t76 - t203 * t145 - t83 * t153 + t84 * t172;
t6 = -qJD(5) * t193 + t108 * t26 + t110 * t118;
t210 = -0.2e1 * t212;
t209 = 0.2e1 * qJD(2);
t208 = 2 * qJD(6);
t207 = pkin(9) * t51;
t206 = pkin(9) * t81;
t205 = t50 * pkin(4);
t204 = t51 * pkin(5);
t202 = t111 * pkin(3);
t99 = t203 * pkin(3) + pkin(9);
t201 = t51 * t99;
t27 = t54 * qJD(4) - t111 * t145 - t203 * t76;
t53 = t111 * t84 - t203 * t83;
t200 = t53 * t27;
t198 = t81 * t99;
t103 = qJD(5) * t110;
t46 = t53 * t103;
t195 = t27 * t108 + t46;
t194 = t222 * pkin(9);
t147 = pkin(3) * t153;
t66 = pkin(5) * t102 - qJ(6) * t103 - t108 * qJD(6);
t55 = t147 + t66;
t192 = -t55 - t66;
t191 = t221 * t99;
t190 = t221 * pkin(9);
t100 = -pkin(4) - t202;
t189 = t100 * t103 + t108 * t147;
t188 = t106 * t50;
t186 = t107 * t50;
t184 = t108 * t51;
t183 = t50 * t100;
t182 = t50 * t108;
t181 = t50 * t110;
t180 = t51 * qJ(6);
t28 = (pkin(5) * t108 - qJ(6) * t110) * t82 + t53;
t179 = qJD(5) * t28;
t176 = t81 * qJD(6);
t173 = qJ(6) * qJD(5);
t37 = 0.2e1 * t199;
t169 = -0.2e1 * t184;
t168 = -t78 * t103 + t182 * t216;
t166 = pkin(4) * t102;
t165 = pkin(4) * t103;
t163 = pkin(9) * t102;
t162 = pkin(9) * t103;
t161 = t108 * t181;
t160 = t203 * t53;
t159 = t82 * t102;
t158 = t81 * t103;
t157 = t82 * t103;
t25 = t28 * t102;
t45 = t53 * t102;
t155 = t108 * t103;
t154 = t109 * t170;
t151 = t221 * t81 + t222 * t99;
t150 = pkin(5) * qJD(5) - qJD(6);
t146 = t78 * t155;
t119 = qJ(6) * t50 + t150 * t82;
t132 = -pkin(5) * t50 + t82 * t173;
t7 = t132 * t108 + t119 * t110 + t27;
t144 = t28 * t50 - t7 * t82;
t143 = t27 * t82 - t53 * t50;
t141 = t50 * t81 - t82 * t51;
t137 = -t110 * pkin(5) - t108 * qJ(6);
t87 = -pkin(4) + t137;
t72 = t87 - t202;
t140 = t50 * t72 - t82 * t55;
t139 = t50 * t87 - t82 * t66;
t138 = t100 * t82 - t198;
t12 = t81 * qJ(6) + t193;
t120 = t110 * t128;
t18 = -t108 * t54 + t120;
t13 = -t81 * pkin(5) - t18;
t136 = t108 * t13 + t110 * t12;
t135 = t108 * t12 - t110 * t13;
t134 = t108 * t193 + t110 * t18;
t133 = t108 * t18 - t110 * t193;
t32 = t157 - t182;
t129 = t159 + t181;
t29 = t81 * t102 - t110 * t51;
t5 = -qJD(5) * t120 + t54 * t102 - t108 * t118 + t110 * t26;
t127 = t100 * t102 - t110 * t147;
t124 = -t139 - t207;
t123 = -t7 + (t82 * t87 - t206) * qJD(5);
t122 = -t7 + (t72 * t82 - t198) * qJD(5);
t117 = t26 * t81 - t54 * t51 + t143;
t65 = t137 * qJD(5) + t110 * qJD(6);
t115 = -t81 * t164 - t140 - t201;
t3 = t176 - t5 + t180;
t4 = -t204 - t6;
t1 = -t135 * qJD(5) + t4 * t108 + t3 * t110;
t2 = -t134 * qJD(5) - t6 * t108 - t5 * t110;
t114 = pkin(3) * t121 - t183 - t201;
t105 = qJ(2) * t209;
t93 = -0.2e1 * t155;
t92 = 0.2e1 * t155;
t71 = t87 * t102;
t62 = t72 * t102;
t61 = t99 * t103 + t108 * t164;
t60 = t99 * t102 - t110 * t164;
t59 = 0.2e1 * t221;
t31 = t158 + t184;
t21 = -0.2e1 * t82 * t186 - 0.2e1 * t146;
t20 = -0.2e1 * t82 * t188 + 0.2e1 * t146;
t17 = t212 * t82 + t161;
t14 = t161 * t216 + t78 * t212;
t11 = 0.4e1 * t82 * t155 + t186 - t188;
t10 = -t141 * t108 + t81 * t157;
t9 = -0.2e1 * t141 * t110 - 0.2e1 * t81 * t159;
t8 = 0.2e1 * t222 * t81 - 0.2e1 * t197;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, t105, -0.2e1 * t154, 0.2e1 * (t109 ^ 2 - t112 ^ 2) * qJD(3), 0, 0.2e1 * t154, 0, 0, 0.2e1 * qJ(2) * t170 + 0.2e1 * qJD(2) * t109, -0.2e1 * qJ(2) * t171 + 0.2e1 * qJD(2) * t112, 0, t105, -0.2e1 * t197, 0.2e1 * t141, 0, t37, 0, 0, 0.2e1 * t97 * t51 + 0.2e1 * t94 * t81, -0.2e1 * t97 * t50 + 0.2e1 * t94 * t82, 0.2e1 * t117, -0.2e1 * t54 * t26 + 0.2e1 * t97 * t94 + 0.2e1 * t200, t21, 0.2e1 * t14, t9, t20, -0.2e1 * t10, t37, 0.2e1 * t143 * t108 + 0.2e1 * t18 * t51 + 0.2e1 * t82 * t46 + 0.2e1 * t6 * t81, 0.2e1 * t143 * t110 - 0.2e1 * t193 * t51 - 0.2e1 * t82 * t45 + 0.2e1 * t5 * t81, 0.2e1 * t134 * t50 + 0.2e1 * (qJD(5) * t133 + t108 * t5 - t110 * t6) * t82, 0.2e1 * t18 * t6 - 0.2e1 * t193 * t5 + 0.2e1 * t200, t21, t9, -0.2e1 * t14, t37, 0.2e1 * t10, t20, -0.2e1 * t108 * t144 - 0.2e1 * t13 * t51 + 0.2e1 * t28 * t157 - 0.2e1 * t4 * t81, 0.2e1 * t135 * t50 + 0.2e1 * (-qJD(5) * t136 - t108 * t3 + t110 * t4) * t82, 0.2e1 * t110 * t144 + 0.2e1 * t12 * t51 + 0.2e1 * t82 * t25 + 0.2e1 * t3 * t81, 0.2e1 * t12 * t3 + 0.2e1 * t13 * t4 + 0.2e1 * t28 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, -t117, 0, 0, 0, 0, 0, 0, -t77 * t103 + t81 * t169 + t168, t217, 0, -t133 * t51 + t2 * t81 - t143, 0, 0, 0, 0, 0, 0 (-t158 + t169) * t81 + t168, 0, -t217, t1 * t81 + t136 * t51 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, 0, -t170, 0, -t113 * t171, -t113 * t170, 0, 0, 0, 0, -t50, 0, -t51, 0, -t27, t26, t220 (-t203 * t26 - t111 * t27 + (t111 * t54 + t160) * qJD(4)) * pkin(3), -t17, -t11, t31, t17, -t29, 0, t45 + (t138 * qJD(5) - t27) * t110 + t114 * t108, -t138 * t102 + t110 * t114 + t195, t2, t27 * t100 + (-t133 * t111 + t160) * qJD(4) * pkin(3) + t2 * t99, -t17, t31, t11, 0, t29, t17, t108 * t115 + t110 * t122 + t25, t1, t122 * t108 + (-t115 - t179) * t110, t1 * t99 + t136 * t164 + t28 * t55 + t7 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, -t170, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t51, 0, -t220, 0, 0, 0, 0, 0, 0, -t129, -t32, t222, -t147 * t82 + t151 + t183, 0, 0, 0, 0, 0, 0, -t129, t222, t32, t140 + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t147, -0.2e1 * t164, 0, 0, t92, t210, 0, t93, 0, 0, 0.2e1 * t127, 0.2e1 * t189, t59, 0.2e1 * t100 * t147 + 0.2e1 * t191, t92, 0, -t210, 0, 0, t93, -0.2e1 * t55 * t110 + 0.2e1 * t62, t59, -0.2e1 * t72 * t103 - 0.2e1 * t55 * t108, 0.2e1 * t72 * t55 + 0.2e1 * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, 0, -t51, 0, -t27, t26, 0, 0, -t17, -t11, t31, t17, -t29, 0, t45 + (t205 - t207) * t108 + (-t27 + (-pkin(4) * t82 - t206) * qJD(5)) * t110, pkin(4) * t129 + pkin(9) * t29 + t195, t2, -t27 * pkin(4) + pkin(9) * t2, -t17, t31, t11, 0, t29, t17, t108 * t124 + t110 * t123 + t25, t1, t123 * t108 + (-t124 - t179) * t110, pkin(9) * t1 + t28 * t66 + t7 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t51, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t32, t222, t194 - t205, 0, 0, 0, 0, 0, 0, -t129, t222, t32, t139 + t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, -t164, 0, 0, t92, t210, 0, t93, 0, 0, t127 - t166, -t165 + t189, t221, -pkin(4) * t147 + t190, t92, 0, -t210, 0, 0, t93, t192 * t110 + t62 + t71, t221, t192 * t108 + (-t72 - t87) * t103, t55 * t87 + t72 * t66 + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t210, 0, t93, 0, 0, -0.2e1 * t166, -0.2e1 * t165, 0, 0, t92, 0, -t210, 0, 0, t93, -0.2e1 * t66 * t110 + 0.2e1 * t71, 0, -0.2e1 * t87 * t103 - 0.2e1 * t66 * t108, 0.2e1 * t87 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, 0, -t32, t51, t6, t5, 0, 0, 0, -t129, 0, t51, t32, 0, t6 + 0.2e1 * t204, t108 * t119 - t110 * t132, 0.2e1 * t176 - t5 + 0.2e1 * t180, -t4 * pkin(5) + t3 * qJ(6) + t12 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t29, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t29 (-t173 * t81 - t204) * t108 + (-t150 * t81 + t180) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, -t102, 0, -t61, t60, 0, 0, 0, t103, 0, 0, t102, 0, -t61, t65, -t60 (-pkin(5) * t164 - t173 * t99) * t108 + (qJ(6) * t164 - t150 * t99) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, -t102, 0, -t162, t163, 0, 0, 0, t103, 0, 0, t102, 0, -t162, t65, -t163, t65 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, qJ(6) * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t129, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
