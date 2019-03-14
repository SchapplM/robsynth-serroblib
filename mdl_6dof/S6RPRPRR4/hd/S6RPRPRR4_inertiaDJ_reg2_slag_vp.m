% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:27
% EndTime: 2019-03-09 03:46:35
% DurationCPUTime: 2.70s
% Computational Cost: add. (2819->233), mult. (5714->407), div. (0->0), fcn. (4802->8), ass. (0->136)
t141 = qJD(5) + qJD(6);
t82 = sin(qJ(5));
t84 = cos(qJ(6));
t164 = t84 * t82;
t81 = sin(qJ(6));
t85 = cos(qJ(5));
t56 = t81 * t85 + t164;
t29 = t141 * t56;
t149 = qJD(6) * t81;
t152 = qJD(5) * t82;
t163 = t84 * t85;
t30 = t141 * t163 - t82 * t149 - t81 * t152;
t57 = -t81 * t82 + t163;
t182 = (-t29 * t84 + t30 * t81 + (t56 * t84 - t57 * t81) * qJD(6)) * pkin(5);
t145 = t85 * qJD(3);
t83 = sin(qJ(3));
t129 = t83 * t145;
t86 = cos(qJ(3));
t150 = qJD(5) * t86;
t132 = t82 * t150;
t49 = t129 + t132;
t120 = sin(pkin(10)) * pkin(1) + pkin(7);
t107 = pkin(4) + t120;
t103 = t107 * t83;
t165 = t82 * t86;
t154 = t83 * qJ(4);
t121 = cos(pkin(10)) * pkin(1) + pkin(2);
t170 = pkin(3) + pkin(8);
t133 = t170 * t86;
t98 = -t133 - t121;
t93 = -t98 + t154;
t89 = t83 * pkin(5) + pkin(9) * t165 + t85 * t103 + t82 * t93;
t146 = t83 * qJD(4);
t155 = qJ(4) * t86;
t177 = t170 * t83;
t174 = -t155 + t177;
t179 = -t174 * qJD(3) - qJD(5) * t103 + t146;
t74 = t86 * qJD(3);
t180 = qJD(5) * t93 + t107 * t74;
t9 = t179 * t85 - t180 * t82;
t181 = -t49 * pkin(9) - qJD(6) * t89 + t9;
t162 = t85 * t86;
t140 = t84 * t162;
t73 = t83 * qJD(3);
t130 = t82 * t73;
t131 = t85 * t150;
t175 = t130 - t131;
t18 = -qJD(6) * t140 + (t141 * t165 + t129) * t81 + t175 * t84;
t41 = t56 * t86;
t159 = -t18 * t56 + t41 * t30;
t19 = t84 * t129 - t81 * t130 + t29 * t86;
t40 = t81 * t165 - t140;
t160 = t57 * t19 - t40 * t29;
t178 = t160 - t159;
t68 = t83 * t74;
t64 = 0.2e1 * t68;
t77 = t82 ^ 2;
t79 = t85 ^ 2;
t157 = t77 + t79;
t142 = pkin(9) + t170;
t124 = t142 * t85;
t176 = t141 * t124;
t80 = t86 ^ 2;
t126 = qJD(3) * (t83 ^ 2 - t80);
t158 = t77 - t79;
t125 = t158 * qJD(5);
t144 = t86 * qJD(4);
t62 = t86 * t120;
t54 = t86 * pkin(4) + t62;
t147 = t54 * qJD(5);
t173 = (-t133 - t154) * qJD(3) + t144 + t147;
t10 = t179 * t82 + t180 * t85;
t94 = t82 * qJ(4) + t85 * t107;
t26 = -t82 * t98 + t94 * t83;
t27 = t82 * t103 - t85 * t93;
t114 = t26 * t82 - t27 * t85;
t3 = -t114 * qJD(5) + t10 * t85 - t9 * t82;
t171 = 0.2e1 * qJD(4);
t169 = t40 * t19;
t168 = t41 * t18;
t167 = t56 * t30;
t166 = t57 * t29;
t153 = qJD(3) * t54;
t151 = qJD(5) * t85;
t148 = qJD(6) * t84;
t143 = qJ(4) * qJD(5);
t139 = -0.2e1 * qJD(3) * t121;
t138 = pkin(5) * t74;
t137 = pkin(5) * t149;
t136 = pkin(5) * t148;
t135 = t82 * t170;
t134 = t85 * t170;
t128 = t82 * t151;
t127 = qJD(5) * t170;
t123 = t82 * t129;
t122 = t80 * t128;
t118 = -t86 * pkin(3) - t154;
t117 = t18 * t57 + t29 * t41;
t116 = -t56 * t19 - t30 * t40;
t115 = t26 * t85 + t27 * t82;
t113 = t166 - t167;
t108 = t142 * t152;
t21 = t83 * t30 + t56 * t74;
t25 = -pkin(9) * t162 + t27;
t87 = -t175 * pkin(9) + t10 + t138;
t1 = t25 * t149 + t181 * t84 - t81 * t87;
t104 = t120 * qJD(3);
t102 = t83 * t104;
t101 = t86 * t104;
t2 = -t25 * t148 + t181 * t81 + t84 * t87;
t7 = -t81 * t25 + (t142 * t86 + t121) * t164 + t84 * (pkin(5) + t94) * t83;
t8 = t84 * t25 + t81 * t89;
t97 = -t1 * t56 + t2 * t57 - t7 * t29 + t8 * t30;
t43 = t107 * t73;
t96 = t174 * qJD(5) - t43;
t58 = t142 * t82;
t15 = -t81 * t108 - t58 * t149 + t176 * t84;
t16 = t84 * t108 + t58 * t148 + t176 * t81;
t31 = -t84 * t124 + t81 * t58;
t32 = -t81 * t124 - t84 * t58;
t95 = t15 * t56 - t16 * t57 + t31 * t29 - t32 * t30;
t76 = qJ(4) * t171;
t71 = t82 * pkin(5) + qJ(4);
t66 = pkin(5) * t151 + qJD(4);
t65 = -0.2e1 * t68;
t55 = -0.2e1 * t126;
t53 = t118 - t121;
t50 = t83 * t151 + t82 * t74;
t48 = -t83 * t152 + t85 * t74;
t47 = t157 * t73;
t44 = t146 + (-pkin(3) * t83 + t155) * qJD(3);
t37 = pkin(5) * t162 + t54;
t34 = -t86 * t125 - t123;
t28 = -pkin(5) * t132 + (-t85 * pkin(5) - t107) * t73;
t22 = -t29 * t83 + t57 * t74;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t55, 0, t65, 0, 0, t83 * t139, t86 * t139, 0, 0, 0, 0, 0, t64, t55, t65, 0, -0.2e1 * t44 * t86 - 0.2e1 * t53 * t73, 0.2e1 * t44 * t83 - 0.2e1 * t53 * t74, -0.2e1 * t53 * t44, -0.2e1 * t77 * t68 + 0.2e1 * t122, -0.4e1 * t86 * t123 - 0.2e1 * t80 * t125, 0.2e1 * t82 * t126 - 0.2e1 * t83 * t131, -0.2e1 * t79 * t68 - 0.2e1 * t122, 0.2e1 * t85 * t126 + 0.2e1 * t83 * t132, t64, 0.2e1 * (-t54 * t145 + t10) * t83 + 0.2e1 * (qJD(3) * t26 - t82 * t147 - t43 * t85) * t86, 0.2e1 * (t82 * t153 + t9) * t83 + 0.2e1 * (-qJD(3) * t27 - t85 * t147 + t43 * t82) * t86, -0.2e1 * t114 * t73 + 0.2e1 * (qJD(5) * t115 + t10 * t82 + t85 * t9) * t86, 0.2e1 * t26 * t10 - 0.2e1 * t27 * t9 - 0.2e1 * t54 * t43, -0.2e1 * t168, 0.2e1 * t18 * t40 - 0.2e1 * t19 * t41, 0.2e1 * t18 * t83 - 0.2e1 * t41 * t74, 0.2e1 * t169, 0.2e1 * t19 * t83 + 0.2e1 * t40 * t74, t64, -0.2e1 * t37 * t19 + 0.2e1 * t2 * t83 - 0.2e1 * t28 * t40 + 0.2e1 * t7 * t74, 0.2e1 * t1 * t83 + 0.2e1 * t37 * t18 - 0.2e1 * t28 * t41 - 0.2e1 * t74 * t8, -0.2e1 * t1 * t40 - 0.2e1 * t7 * t18 + 0.2e1 * t8 * t19 + 0.2e1 * t2 * t41, -0.2e1 * t8 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t37 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (qJD(3) * t115 - t43) * t83 + (t153 - t3) * t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t41 + t8 * t18 + t7 * t19 + t2 * t40 + t28 * t83 + t37 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t157 + 0.1e1) * t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t168 + 0.2e1 * t68 + 0.2e1 * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, -t73, 0, -t101, t102, 0, 0, 0, -t74, t73, 0, 0, 0, t118 * qJD(3) + t144, t101, -t102, qJD(4) * t62 + (-pkin(3) * t62 - t120 * t154) * qJD(3), -t34, 0.4e1 * t86 * t128 - t158 * t73, t48, t34, -t50, 0, t173 * t85 + t96 * t82, -t173 * t82 + t96 * t85, -t3, t9 * t135 - t10 * t134 - t43 * qJ(4) + t54 * qJD(4) + (-t134 * t27 + t135 * t26) * qJD(5), t117, t159 + t160, t22, t116, -t21, 0, t16 * t83 - t71 * t19 + t28 * t56 + t37 * t30 + t31 * t74 - t66 * t40, t15 * t83 + t71 * t18 + t28 * t57 - t37 * t29 - t32 * t74 - t66 * t41, -t15 * t40 + t16 * t41 - t31 * t18 + t32 * t19 - t97, -t1 * t32 - t8 * t15 + t7 * t16 + t2 * t31 + t28 * t71 + t37 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t74, t44, 0, 0, 0, 0, 0, 0, t50, t48, -t47, t146 + (-t157 * t177 + t155) * qJD(3), 0, 0, 0, 0, 0, 0, t21, t22, -t178, t41 * t15 + t40 * t16 + t18 * t32 + t19 * t31 + t83 * t66 + t71 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t76, -0.2e1 * t128, 0.2e1 * t125, 0, 0.2e1 * t128, 0, 0, 0.2e1 * qJD(4) * t82 + 0.2e1 * t85 * t143, 0.2e1 * qJD(4) * t85 - 0.2e1 * t82 * t143, 0, t76, -0.2e1 * t166, 0.2e1 * t29 * t56 - 0.2e1 * t57 * t30, 0, 0.2e1 * t167, 0, 0, 0.2e1 * t71 * t30 + 0.2e1 * t66 * t56, -0.2e1 * t71 * t29 + 0.2e1 * t66 * t57, 0.2e1 * t95, -0.2e1 * t32 * t15 + 0.2e1 * t31 * t16 + 0.2e1 * t71 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, t101, 0, 0, 0, 0, 0, 0, t48, -t50, 0, t3, 0, 0, 0, 0, 0, 0, t22, -t21, -t116 - t117, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t113, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, 0, t49, t74, t10, t9, 0, 0, 0, 0, t18, 0, t19, t74, -t137 * t83 + t138 * t84 + t2 (-t148 * t83 - t74 * t81) * pkin(5) + t1 (-t18 * t84 + t19 * t81 + (t40 * t84 - t41 * t81) * qJD(6)) * pkin(5) (-t1 * t81 + t2 * t84 + (-t7 * t81 + t8 * t84) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t175, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0 (t18 * t81 + t19 * t84 + (-t40 * t81 - t41 * t84) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, 0, -t151, 0, t82 * t127, t85 * t127, 0, 0, 0, 0, -t29, 0, -t30, 0, t16, t15, -t182 (-t15 * t81 + t16 * t84 + (-t31 * t81 + t32 * t84) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t151, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30, 0, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t137, -0.2e1 * t136, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t19, t74, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, -t30, 0, t16, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, -t136, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;