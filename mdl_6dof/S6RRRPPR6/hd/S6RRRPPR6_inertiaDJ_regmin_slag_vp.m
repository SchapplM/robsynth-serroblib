% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:53:55
% EndTime: 2019-03-09 15:54:02
% DurationCPUTime: 2.31s
% Computational Cost: add. (4161->301), mult. (11208->561), div. (0->0), fcn. (10769->10), ass. (0->148)
t117 = sin(qJ(6));
t111 = t117 ^ 2;
t120 = cos(qJ(6));
t172 = -t120 ^ 2 + t111;
t148 = t172 * qJD(6);
t114 = sin(pkin(6));
t184 = 0.2e1 * t114;
t183 = 2 * qJD(5);
t182 = pkin(4) + pkin(10);
t181 = pkin(9) * t114;
t180 = -qJ(4) - pkin(9);
t113 = sin(pkin(11));
t115 = cos(pkin(11));
t119 = sin(qJ(2));
t171 = qJD(2) * t119;
t100 = t114 * t171;
t118 = sin(qJ(3));
t121 = cos(qJ(3));
t116 = cos(pkin(6));
t122 = cos(qJ(2));
t173 = t114 * t122;
t78 = pkin(8) * t173 + (pkin(1) * t119 + pkin(9)) * t116;
t79 = (-pkin(2) * t122 - pkin(9) * t119 - pkin(1)) * t114;
t179 = t118 * t79 + t121 * t78;
t81 = (pkin(2) * t119 - pkin(9) * t122) * t114 * qJD(2);
t170 = qJD(2) * t122;
t82 = -t116 * pkin(1) * t170 + pkin(8) * t100;
t36 = -t179 * qJD(3) + t118 * t82 + t121 * t81;
t153 = t114 * t170;
t174 = t114 * t119;
t87 = -t116 * t121 + t118 * t174;
t66 = -t87 * qJD(3) + t121 * t153;
t88 = t116 * t118 + t121 * t174;
t22 = pkin(3) * t100 - t66 * qJ(4) - t88 * qJD(4) + t36;
t168 = qJD(3) * t121;
t169 = qJD(3) * t118;
t35 = -t118 * t81 + t121 * t82 - t79 * t168 + t78 * t169;
t65 = t88 * qJD(3) + t118 * t153;
t27 = -t65 * qJ(4) - t87 * qJD(4) - t35;
t10 = t113 * t22 + t115 * t27;
t150 = -t118 * t78 + t121 * t79;
t39 = -pkin(3) * t173 - t88 * qJ(4) + t150;
t45 = -t87 * qJ(4) + t179;
t26 = t113 * t39 + t115 * t45;
t161 = pkin(1) * t171;
t83 = pkin(8) * t153 + t116 * t161;
t91 = t113 * t121 + t115 * t118;
t85 = t91 * qJD(3);
t178 = t120 * t85;
t43 = t113 * t66 + t115 * t65;
t58 = t113 * t88 + t115 * t87;
t49 = t117 * t173 + t120 * t58;
t29 = t49 * qJD(6) + t120 * t100 + t117 * t43;
t177 = t29 * t120;
t107 = -t115 * pkin(3) - pkin(4);
t103 = -pkin(10) + t107;
t176 = t103 * t117;
t175 = t103 * t120;
t167 = qJD(3) * t122;
t166 = qJD(5) * t122;
t165 = qJD(6) * t117;
t164 = qJD(6) * t120;
t163 = -0.2e1 * pkin(2) * qJD(3);
t162 = qJ(5) * t100 + t10;
t109 = pkin(3) * t169;
t160 = t117 * t178;
t90 = t113 * t118 - t115 * t121;
t159 = t90 * t165;
t158 = t91 * t165;
t157 = t90 * t164;
t156 = t91 * t164;
t155 = t120 * t173;
t108 = -t121 * pkin(3) - pkin(2);
t110 = t114 ^ 2;
t154 = t110 * t170;
t152 = t117 * t164;
t151 = t180 * t118;
t149 = qJD(3) * t180;
t132 = -t118 * qJD(4) + t121 * t149;
t84 = t121 * qJD(4) + t118 * t149;
t54 = t113 * t84 - t115 * t132;
t94 = t180 * t121;
t67 = -t113 * t94 - t115 * t151;
t9 = -t113 * t27 + t115 * t22;
t25 = -t113 * t45 + t115 * t39;
t52 = t65 * pkin(3) + t83;
t147 = t119 * t154;
t21 = pkin(4) * t173 - t25;
t20 = qJ(5) * t173 - t26;
t17 = -t58 * pkin(5) - t20;
t7 = t114 * t166 - t162;
t4 = -t43 * pkin(5) - t7;
t146 = t17 * t85 + t4 * t90;
t44 = -t113 * t65 + t115 * t66;
t59 = -t113 * t87 + t115 * t88;
t145 = t44 * t90 + t59 * t85;
t55 = t113 * t132 + t115 * t84;
t47 = -t85 * pkin(5) + t55;
t68 = t113 * t151 - t115 * t94;
t57 = -t90 * pkin(5) + t68;
t144 = t47 * t90 + t57 * t85;
t143 = t67 * t54 + t68 * t55;
t86 = -t113 * t169 + t115 * t168;
t142 = t85 * t91 + t86 * t90;
t16 = t59 * pkin(5) + pkin(10) * t173 + t21;
t77 = pkin(8) * t174 + (-pkin(1) * t122 - pkin(2)) * t116;
t62 = t87 * pkin(3) + t77;
t126 = -t59 * qJ(5) + t62;
t18 = t182 * t58 + t126;
t6 = t117 * t16 + t120 * t18;
t50 = t117 * t58 - t155;
t141 = t117 * t49 + t120 * t50;
t139 = -t91 * qJ(5) + t108;
t51 = t182 * t90 + t139;
t56 = t91 * pkin(5) + t67;
t34 = t117 * t56 + t120 * t51;
t104 = t113 * pkin(3) + qJ(5);
t140 = -qJD(5) * t90 - t104 * t85;
t138 = t86 * pkin(5) + t54;
t32 = -t117 * t44 - t59 * t164;
t137 = -t120 * t44 + t59 * t165;
t61 = -t117 * t86 - t156;
t136 = -t120 * t86 + t158;
t135 = -t86 * qJ(5) - t91 * qJD(5) + t109;
t134 = t118 * t167 + t121 * t171;
t133 = t118 * t171 - t121 * t167;
t131 = t47 + (-t103 * t91 + t104 * t90) * qJD(6);
t130 = -t44 * qJ(5) - t59 * qJD(5) + t52;
t129 = -t68 * t43 + t67 * t44 + t54 * t59 - t55 * t58;
t128 = -t182 * t85 - t135;
t127 = -qJD(6) * t57 - t103 * t86 - t140;
t125 = t44 * pkin(5) - t182 * t100 - t9;
t124 = 0.2e1 * t54 * t91 - 0.2e1 * t55 * t90 + 0.2e1 * t67 * t86 - 0.2e1 * t68 * t85;
t123 = t182 * t43 + t130;
t89 = t90 ^ 2;
t63 = t90 * pkin(4) + t139;
t48 = t85 * pkin(4) + t135;
t33 = -t117 * t51 + t120 * t56;
t30 = t58 * pkin(4) + t126;
t28 = -t120 * t43 - qJD(6) * t155 + (qJD(6) * t58 + t100) * t117;
t14 = t43 * pkin(4) + t130;
t13 = -t34 * qJD(6) + t117 * t128 + t120 * t138;
t12 = -t117 * t138 + t120 * t128 - t56 * t164 + t51 * t165;
t8 = -pkin(4) * t100 - t9;
t5 = -t117 * t18 + t120 * t16;
t2 = -t6 * qJD(6) - t117 * t123 + t120 * t125;
t1 = -t117 * t125 - t120 * t123 - t16 * t164 + t18 * t165;
t3 = [0, 0, 0, 0.2e1 * t147, 0.2e1 * (-t119 ^ 2 + t122 ^ 2) * t110 * qJD(2), 0.2e1 * t116 * t153, -0.2e1 * t116 * t100, 0, -0.2e1 * t110 * t161 - 0.2e1 * t83 * t116, -0.2e1 * pkin(1) * t154 + 0.2e1 * t82 * t116, 0.2e1 * t88 * t66, -0.2e1 * t88 * t65 - 0.2e1 * t66 * t87 (-t122 * t66 + t88 * t171) * t184 (t122 * t65 - t171 * t87) * t184, -0.2e1 * t147, 0.2e1 * t77 * t65 + 0.2e1 * t83 * t87 + 0.2e1 * (-t36 * t122 + t150 * t171) * t114, 0.2e1 * t77 * t66 + 0.2e1 * t83 * t88 + 0.2e1 * (-t35 * t122 - t179 * t171) * t114, -0.2e1 * t10 * t58 - 0.2e1 * t25 * t44 - 0.2e1 * t26 * t43 - 0.2e1 * t9 * t59, 0.2e1 * t26 * t10 + 0.2e1 * t25 * t9 + 0.2e1 * t62 * t52, 0.2e1 * t20 * t43 + 0.2e1 * t21 * t44 + 0.2e1 * t7 * t58 + 0.2e1 * t8 * t59, -0.2e1 * t14 * t58 - 0.2e1 * t30 * t43 + 0.2e1 * (-t122 * t8 + t21 * t171) * t114, -0.2e1 * t14 * t59 - 0.2e1 * t30 * t44 + 0.2e1 * (t122 * t7 - t20 * t171) * t114, 0.2e1 * t30 * t14 + 0.2e1 * t20 * t7 + 0.2e1 * t21 * t8, 0.2e1 * t50 * t29, -0.2e1 * t50 * t28 + 0.2e1 * t29 * t49, 0.2e1 * t29 * t59 + 0.2e1 * t50 * t44, -0.2e1 * t28 * t59 + 0.2e1 * t49 * t44, 0.2e1 * t59 * t44, 0.2e1 * t17 * t28 + 0.2e1 * t2 * t59 - 0.2e1 * t4 * t49 + 0.2e1 * t5 * t44, 0.2e1 * t1 * t59 + 0.2e1 * t17 * t29 + 0.2e1 * t4 * t50 - 0.2e1 * t6 * t44; 0, 0, 0, 0, 0, t153, -t100, 0, -t83, t82, t66 * t118 + t88 * t168, -t118 * t65 + t66 * t121 + (-t118 * t88 - t121 * t87) * qJD(3), t133 * t114, t134 * t114, 0, -pkin(2) * t65 - t83 * t121 - t133 * t181 + t77 * t169, -pkin(2) * t66 + t83 * t118 - t134 * t181 + t77 * t168, -t10 * t90 - t25 * t86 - t26 * t85 - t9 * t91 + t129, t10 * t68 + t52 * t108 + t62 * t109 - t25 * t54 + t26 * t55 - t9 * t67, t20 * t85 + t21 * t86 + t7 * t90 + t8 * t91 + t129, -t14 * t90 - t30 * t85 - t63 * t43 - t48 * t58 + (-t122 * t54 + t67 * t171) * t114, -t14 * t91 - t30 * t86 - t63 * t44 - t48 * t59 + (-t122 * t55 + t68 * t171) * t114, t14 * t63 - t20 * t55 + t21 * t54 + t30 * t48 + t8 * t67 - t7 * t68, t50 * t157 + (t29 * t90 + t50 * t85) * t117, t141 * t85 + (-t117 * t28 + t177 + (-t117 * t50 + t120 * t49) * qJD(6)) * t90, t145 * t117 + t59 * t157 + t29 * t91 + t50 * t86, t145 * t120 - t59 * t159 - t28 * t91 + t49 * t86, t44 * t91 + t59 * t86, -t120 * t146 + t13 * t59 + t17 * t159 + t2 * t91 + t57 * t28 + t33 * t44 - t47 * t49 + t5 * t86, t1 * t91 + t117 * t146 + t12 * t59 + t17 * t157 + t57 * t29 - t34 * t44 + t47 * t50 - t6 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t118 * t168, 0.2e1 * (-t118 ^ 2 + t121 ^ 2) * qJD(3), 0, 0, 0, t118 * t163, t121 * t163, t124, 0.2e1 * t108 * t109 + 0.2e1 * t143, t124, -0.2e1 * t48 * t90 - 0.2e1 * t63 * t85, -0.2e1 * t48 * t91 - 0.2e1 * t63 * t86, 0.2e1 * t63 * t48 + 0.2e1 * t143, 0.2e1 * t111 * t90 * t85 + 0.2e1 * t89 * t152, -0.2e1 * t89 * t148 + 0.4e1 * t90 * t160, 0.2e1 * t142 * t117 + 0.2e1 * t90 * t156, 0.2e1 * t142 * t120 - 0.2e1 * t90 * t158, 0.2e1 * t91 * t86, -0.2e1 * t120 * t144 + 0.2e1 * t13 * t91 + 0.2e1 * t57 * t159 + 0.2e1 * t33 * t86, 0.2e1 * t117 * t144 + 0.2e1 * t12 * t91 + 0.2e1 * t57 * t157 - 0.2e1 * t34 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t65, t100, t36, t35 (-t113 * t43 - t115 * t44) * pkin(3) (t10 * t113 + t115 * t9) * pkin(3), -qJD(5) * t58 - t104 * t43 + t107 * t44 (-pkin(4) + t107) * t100 - t9 (t104 * t171 - 0.2e1 * t166) * t114 + t162, -t20 * qJD(5) - t7 * t104 + t8 * t107, -t50 * t165 + t177, -t141 * qJD(6) - t29 * t117 - t120 * t28, -t137, t32, 0, t44 * t175 - qJD(5) * t49 + t104 * t28 + t4 * t117 + (t120 * t17 - t59 * t176) * qJD(6), -t44 * t176 + qJD(5) * t50 + t104 * t29 + t4 * t120 + (-t117 * t17 - t59 * t175) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, -t169, 0, -pkin(9) * t168, pkin(9) * t169 (-t113 * t85 - t115 * t86) * pkin(3) (t113 * t55 - t115 * t54) * pkin(3), t107 * t86 + t140, t54, t55, t68 * qJD(5) + t55 * t104 + t54 * t107, -t90 * t148 + t160, -0.4e1 * t90 * t152 - t172 * t85, -t136, t61, 0, t117 * t131 - t120 * t127, t117 * t127 + t120 * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t104 * t183, -0.2e1 * t152, 0.2e1 * t148, 0, 0, 0, 0.2e1 * qJD(5) * t117 + 0.2e1 * t104 * t164, 0.2e1 * qJD(5) * t120 - 0.2e1 * t104 * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, -t43, -t44, t14, 0, 0, 0, 0, 0, t32, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, -t85, -t86, t48, 0, 0, 0, 0, 0, t61, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t100, 0, t8, 0, 0, 0, 0, 0, -t137, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, t54, 0, 0, 0, 0, 0, -t136, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, t44, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t85 + t157, -t159 + t178, t86, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, -t164, 0, -t103 * t165, -t103 * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, -t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
