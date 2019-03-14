% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaDJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:22
% EndTime: 2019-03-08 22:00:30
% DurationCPUTime: 3.52s
% Computational Cost: add. (5219->310), mult. (12533->573), div. (0->0), fcn. (12733->12), ass. (0->152)
t102 = cos(qJ(3));
t184 = -qJ(4) - pkin(8);
t146 = qJD(3) * t184;
t99 = sin(qJ(3));
t117 = -qJD(4) * t99 + t102 * t146;
t118 = qJD(4) * t102 + t99 * t146;
t175 = cos(pkin(12));
t94 = sin(pkin(12));
t109 = t94 * t117 + t175 * t118;
t144 = t175 * t102;
t76 = t94 * t99 - t144;
t77 = t94 * t102 + t175 * t99;
t90 = -pkin(3) * t102 - pkin(2);
t119 = pkin(4) * t76 - pkin(9) * t77 + t90;
t202 = -qJD(5) * t119 - t109;
t114 = t76 * t184;
t172 = t99 * qJD(3);
t163 = pkin(3) * t172;
t70 = t77 * qJD(3);
t71 = qJD(3) * t144 - t94 * t172;
t201 = -pkin(4) * t70 + pkin(9) * t71 + qJD(5) * t114 - t163;
t193 = cos(qJ(6));
t149 = t193 * qJD(6);
t200 = t193 * qJD(5) + t149;
t101 = cos(qJ(5));
t169 = qJD(5) * t101;
t98 = sin(qJ(5));
t185 = t98 * t71;
t120 = t77 * t169 + t185;
t174 = qJD(5) * t98;
t160 = t77 * t174;
t199 = -t101 * t71 + t160;
t92 = t98 ^ 2;
t93 = t101 ^ 2;
t181 = t92 - t93;
t145 = qJD(5) * t181;
t198 = qJD(5) + qJD(6);
t103 = cos(qJ(2));
t95 = sin(pkin(6));
t177 = t103 * t95;
t100 = sin(qJ(2));
t180 = t100 * t95;
t96 = cos(pkin(6));
t72 = t102 * t96 - t99 * t180;
t73 = t102 * t180 + t96 * t99;
t46 = t175 * t73 + t94 * t72;
t122 = t101 * t177 + t98 * t46;
t123 = -t101 * t46 + t98 * t177;
t171 = qJD(2) * t100;
t155 = t95 * t171;
t170 = qJD(2) * t103;
t154 = t95 * t170;
t58 = t72 * qJD(3) + t102 * t154;
t59 = -t73 * qJD(3) - t99 * t154;
t33 = t175 * t58 + t94 * t59;
t15 = t122 * qJD(5) - t101 * t33 - t98 * t155;
t16 = t123 * qJD(5) + t101 * t155 - t98 * t33;
t197 = (-t101 * t123 + t122 * t98) * qJD(5) + t101 * t16 - t15 * t98;
t10 = t202 * t101 + t201 * t98;
t11 = -t201 * t101 + t202 * t98;
t29 = t101 * t119 - t98 * t114;
t30 = t101 * t114 + t98 * t119;
t196 = -t10 * t98 + t101 * t11 + (t101 * t30 - t29 * t98) * qJD(5);
t195 = t70 * pkin(5);
t87 = pkin(3) * t94 + pkin(9);
t194 = pkin(10) + t87;
t32 = -t175 * t59 + t58 * t94;
t45 = -t175 * t72 + t73 * t94;
t22 = t45 * t32;
t43 = -t175 * t117 + t94 * t118;
t60 = t77 * t184;
t192 = t60 * t43;
t191 = t77 * t71;
t190 = t77 * t98;
t158 = t193 * t98;
t97 = sin(qJ(6));
t80 = t97 * t101 + t158;
t57 = t198 * t80;
t151 = t193 * t101;
t186 = t97 * t98;
t79 = -t151 + t186;
t189 = t79 * t57;
t56 = -t200 * t101 + t198 * t186;
t188 = t80 * t56;
t166 = t77 * t186;
t18 = t71 * t158 - t97 * t160 - qJD(6) * t166 + (t200 * t77 + t71 * t97) * t101;
t47 = t80 * t77;
t183 = -t80 * t18 + t56 * t47;
t178 = t101 * t98;
t173 = qJD(6) * t97;
t168 = t102 * qJD(3);
t54 = 0.2e1 * t76 * t70;
t167 = -0.2e1 * pkin(2) * qJD(3);
t88 = -t175 * pkin(3) - pkin(4);
t165 = 0.2e1 * qJD(5) * t88;
t164 = pkin(5) * t174;
t162 = pkin(5) * t173;
t161 = t45 * t174;
t159 = t97 * t194;
t157 = t103 * t172;
t153 = t98 * t169;
t152 = t99 * t168;
t147 = -0.4e1 * t77 * t178;
t143 = t98 * t159;
t75 = t77 ^ 2;
t142 = t75 * t153;
t141 = pkin(5) * t149;
t140 = t100 * t95 ^ 2 * t170;
t139 = t194 * t193;
t17 = -t71 * t151 + t97 * t185 + t57 * t77;
t48 = t77 * t151 - t166;
t138 = -t17 * t79 + t48 * t57;
t137 = -t32 * t60 + t45 * t43;
t136 = t32 * t77 + t45 * t71;
t135 = t43 * t77 - t60 * t71;
t134 = t56 * t76 - t70 * t80;
t133 = t70 * t77 + t71 * t76;
t132 = -t70 * t87 + t71 * t88;
t131 = t76 * t87 - t77 * t88;
t127 = -t101 * t29 - t30 * t98;
t125 = t101 * t122 + t123 * t98;
t124 = t98 * t139;
t104 = t199 * pkin(10) + t11 + t195;
t107 = -t101 * t77 * pkin(10) + t76 * pkin(5) + t29;
t105 = t193 * t107;
t108 = -t120 * pkin(10) - t10;
t23 = -pkin(10) * t190 + t30;
t1 = -qJD(6) * t105 - t97 * t104 - t193 * t108 + t23 * t173;
t20 = -t122 * t97 - t123 * t193;
t121 = t76 * t169 + t70 * t98;
t74 = t194 * t101;
t52 = t193 * t74 - t143;
t113 = t127 * qJD(5) - t10 * t101 - t11 * t98;
t112 = t125 * qJD(5) - t101 * t15 - t16 * t98;
t111 = t102 * t58 - t59 * t99 + (-t102 * t72 - t73 * t99) * qJD(3);
t106 = t97 * t107;
t2 = -qJD(6) * t106 + t193 * t104 - t97 * t108 - t23 * t149;
t81 = -t101 * pkin(5) + t88;
t51 = -t97 * t74 - t124;
t50 = t101 * t70 - t76 * t174;
t40 = pkin(5) * t190 - t60;
t34 = t77 * t145 - t71 * t178;
t27 = -t57 * t76 - t70 * t79;
t26 = -t52 * qJD(6) + (-t101 * t139 + t143) * qJD(5);
t25 = t198 * t124 + t159 * t169 + t74 * t173;
t24 = t120 * pkin(5) + t43;
t19 = -t122 * t193 + t123 * t97;
t9 = t193 * t23 + t106;
t8 = -t97 * t23 + t105;
t4 = -t20 * qJD(6) + t97 * t15 + t193 * t16;
t3 = t122 * t149 - t123 * t173 + t193 * t15 - t97 * t16;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t73 * t58 + 0.2e1 * t72 * t59 - 0.2e1 * t140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t46 * t33 - 0.2e1 * t140 + 0.2e1 * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t122 * t16 + 0.2e1 * t123 * t15 + 0.2e1 * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t19 * t4 - 0.2e1 * t20 * t3 + 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t154, 0, 0, 0, 0, 0, 0, 0, 0 (-t102 * t171 - t157) * t95 (-t103 * t168 + t99 * t171) * t95, t111, -pkin(2) * t155 + t111 * pkin(8), 0, 0, 0, 0, 0, 0 (-t103 * t70 + t76 * t171) * t95 (-t103 * t71 + t77 * t171) * t95, -t33 * t76 - t46 * t70 + t136, -t95 * pkin(3) * t157 + t46 * t109 + t33 * t114 + t90 * t155 + t137, 0, 0, 0, 0, 0, 0, t120 * t45 - t122 * t70 + t16 * t76 + t32 * t190, t101 * t136 + t123 * t70 + t15 * t76 - t160 * t45, t125 * t71 - t197 * t77, t10 * t123 - t11 * t122 - t15 * t30 + t16 * t29 + t137, 0, 0, 0, 0, 0, 0, t18 * t45 + t19 * t70 + t32 * t47 + t4 * t76, -t17 * t45 - t20 * t70 + t3 * t76 + t32 * t48, t17 * t19 - t18 * t20 + t3 * t47 - t4 * t48, -t1 * t20 + t19 * t2 + t24 * t45 - t3 * t9 + t32 * t40 + t4 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t152, 0.2e1 * (t102 ^ 2 - t99 ^ 2) * qJD(3), 0, -0.2e1 * t152, 0, 0, t99 * t167, t102 * t167, 0, 0, 0.2e1 * t191, -0.2e1 * t133, 0, t54, 0, 0, 0.2e1 * t76 * t163 + 0.2e1 * t70 * t90, 0.2e1 * t77 * t163 + 0.2e1 * t71 * t90, -0.2e1 * t109 * t76 - 0.2e1 * t114 * t70 + 0.2e1 * t135, 0.2e1 * t114 * t109 + 0.2e1 * t90 * t163 - 0.2e1 * t192, 0.2e1 * t93 * t191 - 0.2e1 * t142, 0.2e1 * t75 * t145 + t71 * t147, 0.2e1 * t133 * t101 - 0.2e1 * t76 * t160, 0.2e1 * t92 * t191 + 0.2e1 * t142, -0.2e1 * t120 * t76 - 0.2e1 * t70 * t190, t54, 0.2e1 * t11 * t76 - 0.2e1 * t120 * t60 + 0.2e1 * t43 * t190 + 0.2e1 * t29 * t70, 0.2e1 * t10 * t76 + 0.2e1 * t101 * t135 + 0.2e1 * t160 * t60 - 0.2e1 * t30 * t70, 0.2e1 * t127 * t71 - 0.2e1 * t196 * t77, -0.2e1 * t10 * t30 + 0.2e1 * t11 * t29 - 0.2e1 * t192, -0.2e1 * t48 * t17, 0.2e1 * t17 * t47 - 0.2e1 * t18 * t48, -0.2e1 * t17 * t76 + 0.2e1 * t48 * t70, 0.2e1 * t47 * t18, -0.2e1 * t18 * t76 - 0.2e1 * t47 * t70, t54, 0.2e1 * t18 * t40 + 0.2e1 * t2 * t76 + 0.2e1 * t24 * t47 + 0.2e1 * t70 * t8, 0.2e1 * t1 * t76 - 0.2e1 * t17 * t40 + 0.2e1 * t24 * t48 - 0.2e1 * t70 * t9, 0.2e1 * t1 * t47 + 0.2e1 * t17 * t8 - 0.2e1 * t18 * t9 - 0.2e1 * t2 * t48, -0.2e1 * t1 * t9 + 0.2e1 * t2 * t8 + 0.2e1 * t24 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t58, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0 (-t175 * t32 + t33 * t94) * pkin(3), 0, 0, 0, 0, 0, 0, -t101 * t32 + t161, t169 * t45 + t32 * t98, t112, t112 * t87 + t32 * t88, 0, 0, 0, 0, 0, 0, t32 * t79 + t45 * t57, t32 * t80 - t45 * t56, t19 * t56 - t20 * t57 + t3 * t79 - t4 * t80, pkin(5) * t161 + t19 * t26 - t20 * t25 - t3 * t52 + t32 * t81 + t4 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, 0, -t172, 0, -pkin(8) * t168, pkin(8) * t172, 0, 0, 0, 0, t71, 0, -t70, 0, -t43, -t109 (-t175 * t71 - t70 * t94) * pkin(3) (t109 * t94 - t43 * t175) * pkin(3), -t34, qJD(5) * t147 - t181 * t71, t121, t34, t50, 0, -t101 * t43 + t132 * t98 + (-t101 * t131 - t60 * t98) * qJD(5), t43 * t98 + t132 * t101 + (-t101 * t60 + t131 * t98) * qJD(5), t113, t113 * t87 + t43 * t88, -t17 * t80 - t48 * t56, -t138 + t183, -t134, t18 * t79 + t47 * t57, t27, 0, t164 * t47 + t18 * t81 + t24 * t79 + t26 * t76 + t40 * t57 + t51 * t70, t164 * t48 - t17 * t81 + t24 * t80 + t25 * t76 - t40 * t56 - t52 * t70, t1 * t79 + t17 * t51 - t18 * t52 - t2 * t80 + t25 * t47 - t26 * t48 + t56 * t8 - t57 * t9, -t1 * t52 + t164 * t40 + t2 * t51 + t24 * t81 - t25 * t9 + t26 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t153, -0.2e1 * t145, 0, -0.2e1 * t153, 0, 0, t98 * t165, t101 * t165, 0, 0, -0.2e1 * t188, 0.2e1 * t79 * t56 - 0.2e1 * t80 * t57, 0, 0.2e1 * t189, 0, 0, 0.2e1 * t164 * t79 + 0.2e1 * t57 * t81, 0.2e1 * t164 * t80 - 0.2e1 * t56 * t81, 0.2e1 * t25 * t79 - 0.2e1 * t26 * t80 + 0.2e1 * t51 * t56 - 0.2e1 * t52 * t57, 0.2e1 * t164 * t81 - 0.2e1 * t25 * t52 + 0.2e1 * t26 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t57 - t20 * t56 - t3 * t80 - t4 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t71, 0, t163, 0, 0, 0, 0, 0, 0, t50, -t121 (-t92 - t93) * t71, t196, 0, 0, 0, 0, 0, 0, t27, t134, t138 + t183, -t1 * t80 - t2 * t79 - t56 * t9 - t57 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t80 - t26 * t79 - t51 * t57 - t52 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t188 + 0.2e1 * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0 (t193 * t4 - t3 * t97 + (-t19 * t97 + t193 * t20) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, -t120, t70, t11, t10, 0, 0, 0, 0, -t17, 0, -t18, t70, -t162 * t76 + t193 * t195 + t2 (-t149 * t76 - t70 * t97) * pkin(5) + t1 (t193 * t17 - t18 * t97 + (-t193 * t47 + t48 * t97) * qJD(6)) * pkin(5) (t193 * t2 - t1 * t97 + (t193 * t9 - t8 * t97) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, 0, -t174, 0, -t87 * t169, t87 * t174, 0, 0, 0, 0, -t56, 0, -t57, 0, t26, t25 (t193 * t56 - t57 * t97 + (-t193 * t79 + t80 * t97) * qJD(6)) * pkin(5) (t193 * t26 - t25 * t97 + (t193 * t52 - t51 * t97) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, -t169, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t56, 0 (-t193 * t57 - t56 * t97 + (t193 * t80 + t79 * t97) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t162, -0.2e1 * t141, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, -t18, t70, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, 0, -t57, 0, t26, t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, -t141, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;