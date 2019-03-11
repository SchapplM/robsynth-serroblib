% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:06
% EndTime: 2019-03-08 20:25:12
% DurationCPUTime: 2.03s
% Computational Cost: add. (2459->239), mult. (6094->356), div. (0->0), fcn. (4856->12), ass. (0->155)
t106 = qJD(4) + qJD(5);
t114 = sin(qJ(5));
t118 = cos(qJ(4));
t196 = cos(qJ(5));
t152 = qJD(2) * t196;
t115 = sin(qJ(4));
t168 = qJD(2) * t115;
t201 = -t114 * t168 + t118 * t152;
t51 = t201 * t106;
t78 = qJD(6) - t201;
t200 = qJD(6) - t78;
t104 = pkin(4) * t114 + pkin(10);
t162 = pkin(4) * t168;
t173 = t114 * t118;
t81 = -qJD(2) * t173 - t115 * t152;
t54 = -pkin(5) * t81 - pkin(10) * t201;
t199 = (qJD(6) * t104 + t162 + t54) * t78;
t167 = qJD(4) * t115;
t109 = sin(pkin(12));
t119 = cos(qJ(2));
t110 = sin(pkin(6));
t169 = qJD(1) * t110;
t154 = t119 * t169;
t111 = cos(pkin(12));
t116 = sin(qJ(2));
t155 = t116 * t169;
t89 = t111 * t155;
t70 = t109 * t154 + t89;
t140 = pkin(4) * t167 - t70;
t87 = qJD(2) * pkin(2) + t154;
t64 = t109 * t87 + t89;
t149 = t64 + (pkin(8) + pkin(9)) * qJD(2);
t112 = cos(pkin(6));
t98 = qJD(1) * t112 + qJD(3);
t43 = -t149 * t115 + t118 * t98;
t74 = (t109 * t116 - t111 * t119) * t110;
t126 = -t114 * t115 + t118 * t196;
t44 = t115 * t98 + t118 * t149;
t178 = t114 * t44;
t42 = qJD(4) * pkin(4) + t43;
t17 = t196 * t42 - t178;
t14 = -t106 * pkin(5) - t17;
t101 = pkin(2) * t109 + pkin(8);
t186 = pkin(9) + t101;
t83 = t186 * t115;
t84 = t186 * t118;
t129 = -t114 * t84 - t196 * t83;
t88 = t109 * t155;
t73 = t111 * t154 - t88;
t146 = qJD(4) * t186;
t76 = t115 * t146;
t77 = t118 * t146;
t184 = -t129 * qJD(5) + t114 * t77 + t126 * t73 + t196 * t76;
t157 = -pkin(4) * t118 - pkin(3);
t63 = t111 * t87 - t88;
t53 = qJD(2) * t157 - t63;
t29 = -pkin(5) * t201 + pkin(10) * t81 + t53;
t151 = qJD(5) * t196;
t166 = qJD(5) * t114;
t72 = qJD(2) * t74;
t69 = qJD(1) * t72;
t23 = qJD(4) * t43 - t118 * t69;
t24 = -qJD(4) * t44 + t115 * t69;
t3 = t114 * t24 + t42 * t151 - t44 * t166 + t196 * t23;
t147 = t114 * t23 - t196 * t24;
t160 = t196 * t44;
t18 = t114 * t42 + t160;
t4 = qJD(5) * t18 + t147;
t86 = t115 * t196 + t173;
t195 = pkin(2) * t111;
t91 = t157 - t195;
t48 = -pkin(5) * t126 - pkin(10) * t86 + t91;
t50 = -t114 * t83 + t196 * t84;
t58 = t106 * t86;
t52 = t58 * qJD(2);
t57 = t106 * t126;
t198 = (qJD(6) * t29 + t3) * t126 + t14 * t57 + t4 * t86 + (-qJD(6) * t48 + t184) * t78 - t50 * t52;
t113 = sin(qJ(6));
t117 = cos(qJ(6));
t67 = t106 * t113 - t117 * t81;
t31 = t67 * qJD(6) + t113 * t51;
t191 = t52 * t86;
t137 = t57 * t78 + t191;
t165 = qJD(6) * t113;
t159 = t86 * t165;
t197 = -t117 * t137 + t78 * t159;
t194 = t14 * t201;
t193 = t14 * t86;
t192 = t48 * t52;
t65 = -t117 * t106 - t113 * t81;
t190 = t65 * t78;
t189 = t67 * t78;
t188 = t78 * t81;
t187 = t81 * t201;
t164 = qJD(6) * t117;
t30 = t106 * t164 + t117 * t51 + t165 * t81;
t185 = -t126 * t30 + t67 * t58;
t183 = t50 * qJD(5) - t114 * t76 + t196 * t77 - t86 * t73;
t182 = pkin(5) * t58 - pkin(10) * t57 + t140;
t181 = t113 * t30;
t179 = t113 * t52;
t176 = t117 * t52;
t175 = t57 * t106;
t121 = qJD(2) ^ 2;
t174 = t110 * t121;
t120 = qJD(4) ^ 2;
t172 = t120 * t115;
t171 = t120 * t118;
t170 = t115 ^ 2 - t118 ^ 2;
t163 = qJD(2) * qJD(4);
t158 = t78 * t164;
t15 = t106 * pkin(10) + t18;
t135 = t113 * t15 - t117 * t29;
t156 = -t135 * t81 + t14 * t165;
t150 = t115 * t163;
t61 = -qJD(2) * pkin(3) - t63;
t145 = -qJD(2) * t61 + t69;
t144 = t117 * t78;
t8 = t113 * t29 + t117 * t15;
t142 = t4 * t113 + t14 * t164 - t8 * t81;
t20 = t114 * t43 + t160;
t139 = pkin(4) * t166 - t20;
t138 = t126 * t31 - t58 * t65;
t75 = (t109 * t119 + t111 * t116) * t110;
t59 = t112 * t118 - t115 * t75;
t60 = t112 * t115 + t118 * t75;
t33 = t114 * t59 + t196 * t60;
t134 = t113 * t74 + t117 * t33;
t133 = -t113 * t33 + t117 * t74;
t131 = t53 * t81 - t147;
t130 = -t114 * t60 + t196 * t59;
t71 = qJD(2) * t75;
t68 = qJD(1) * t71;
t128 = qJD(2) * t70 - t101 * t120 - t68;
t127 = qJD(4) * (qJD(2) * (-pkin(3) - t195) + t61 + t73);
t56 = pkin(4) * t150 + t68;
t125 = -t113 * t137 - t158 * t86;
t21 = t196 * t43 - t178;
t124 = -t104 * t52 - t194 + (-pkin(4) * t151 + t21) * t78;
t122 = -t201 * t53 - t3;
t105 = -pkin(4) * t196 - pkin(5);
t55 = t58 * t106;
t45 = -t201 ^ 2 + t81 ^ 2;
t40 = (-qJD(2) * t86 - t81) * t106;
t35 = qJD(4) * t59 - t118 * t72;
t34 = -qJD(4) * t60 + t115 * t72;
t19 = pkin(5) * t52 - pkin(10) * t51 + t56;
t16 = t117 * t19;
t11 = t144 * t78 + t67 * t81 + t179;
t10 = -t113 * t78 ^ 2 - t65 * t81 + t176;
t9 = t144 * t67 + t181;
t6 = t33 * qJD(5) + t114 * t35 - t196 * t34;
t5 = t130 * qJD(5) + t114 * t34 + t196 * t35;
t1 = (t30 - t190) * t117 + (-t31 - t189) * t113;
t2 = [0, 0, -t116 * t174, -t119 * t174, -t63 * t71 - t64 * t72 + t68 * t74 - t69 * t75, 0, 0, 0, 0, 0, qJD(4) * t34 + (-t118 * t71 + t167 * t74) * qJD(2), -qJD(4) * t35 + (qJD(4) * t118 * t74 + t115 * t71) * qJD(2), 0, 0, 0, 0, 0, -t106 * t6 - t201 * t71 + t74 * t52, -t106 * t5 + t51 * t74 - t71 * t81, 0, 0, 0, 0, 0 (-qJD(6) * t134 - t113 * t5 + t117 * t71) * t78 + t133 * t52 + t6 * t65 - t130 * t31 -(qJD(6) * t133 + t113 * t71 + t117 * t5) * t78 - t134 * t52 + t6 * t67 - t130 * t30; 0, 0, 0, 0, t63 * t70 - t64 * t73 + (-t109 * t69 - t111 * t68) * pkin(2), 0.2e1 * t118 * t150, -0.2e1 * t170 * t163, t171, -t172, 0, t115 * t127 + t118 * t128, -t115 * t128 + t118 * t127, t51 * t86 - t57 * t81, t126 * t51 + t201 * t57 + t58 * t81 - t191, t175, -t55, 0, -t106 * t183 - t126 * t56 - t140 * t201 + t52 * t91 + t53 * t58, t106 * t184 - t140 * t81 + t51 * t91 + t53 * t57 + t56 * t86, -t67 * t159 + (t30 * t86 + t57 * t67) * t117 (-t113 * t67 - t117 * t65) * t57 + (-t181 - t117 * t31 + (t113 * t65 - t117 * t67) * qJD(6)) * t86, t185 - t197, t125 + t138, -t126 * t52 + t58 * t78, -t16 * t126 - t129 * t31 - t135 * t58 + t183 * t65 + (t192 + t182 * t78 + (t126 * t15 - t50 * t78 + t193) * qJD(6)) * t117 + t198 * t113, -t129 * t30 - t8 * t58 + t183 * t67 + (-t192 + (-qJD(6) * t15 + t19) * t126 - qJD(6) * t193 + (qJD(6) * t50 - t182) * t78) * t113 + t198 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -t171, 0, 0, 0, 0, 0, -t55, -t175, 0, 0, 0, 0, 0, t125 - t138, t185 + t197; 0, 0, 0, 0, 0, -t115 * t121 * t118, t170 * t121, 0, 0, 0, t145 * t115, t145 * t118, t187, t45, 0, t40, 0, t201 * t162 + t20 * t106 + (-t160 + (-pkin(4) * t106 - t42) * t114) * qJD(5) + t131, t21 * t106 + (-t106 * t151 + t168 * t81) * pkin(4) + t122, t9, t1, t11, t10, t188, t105 * t31 + t139 * t65 + (-t4 - t199) * t117 + t124 * t113 + t156, t105 * t30 + t113 * t199 + t117 * t124 + t139 * t67 + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, t45, 0, t40, 0, t131 + (-qJD(5) + t106) * t18, t17 * t106 + t122, t9, t1, t11, t10, t188, -pkin(5) * t31 - t4 * t117 - (-t113 * t17 + t117 * t54) * t78 - t18 * t65 - t113 * t194 + (-t158 - t179) * pkin(10) + t156, -pkin(5) * t30 + (t113 * t54 + t117 * t17) * t78 - t18 * t67 - t117 * t194 + (t165 * t78 - t176) * pkin(10) + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t65, -t65 ^ 2 + t67 ^ 2, t30 + t190, t189 - t31, t52, -t113 * t3 - t14 * t67 - t200 * t8 + t16, -t113 * t19 - t117 * t3 + t200 * t135 + t14 * t65;];
tauc_reg  = t2;
