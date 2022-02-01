% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:19
% EndTime: 2022-01-23 09:30:23
% DurationCPUTime: 1.16s
% Computational Cost: add. (1928->198), mult. (4670->252), div. (0->0), fcn. (3045->6), ass. (0->127)
t104 = sin(qJ(4));
t105 = sin(qJ(3));
t106 = cos(qJ(3));
t161 = cos(qJ(4));
t82 = t104 * t106 + t161 * t105;
t147 = qJD(1) * t82;
t148 = t147 * qJ(5);
t93 = sin(pkin(8)) * pkin(1) + pkin(6);
t84 = t93 * qJD(1);
t125 = pkin(7) * qJD(1) + t84;
t138 = t105 * qJD(2);
t56 = t125 * t106 + t138;
t50 = t104 * t56;
t118 = t125 * t105;
t98 = t106 * qJD(2);
t55 = t98 - t118;
t53 = qJD(3) * pkin(3) + t55;
t18 = t161 * t53 - t50;
t14 = t18 - t148;
t140 = qJD(3) * t105;
t134 = pkin(3) * t140;
t164 = 0.2e1 * t134;
t99 = qJD(3) + qJD(4);
t163 = t147 ^ 2;
t162 = pkin(7) + t93;
t146 = t104 * t105;
t117 = t99 * t146;
t128 = t161 * qJD(4);
t131 = t161 * t106;
t45 = -qJD(3) * t131 - t106 * t128 + t117;
t160 = t45 * t99;
t121 = qJD(1) * t131;
t141 = qJD(1) * t105;
t130 = t104 * t141;
t74 = -t121 + t130;
t159 = t74 * t99;
t158 = t147 * t74;
t94 = -cos(pkin(8)) * pkin(1) - pkin(2);
t83 = -t106 * pkin(3) + t94;
t78 = qJD(1) * t83;
t156 = t78 * t147;
t11 = t99 * pkin(4) + t14;
t155 = t11 - t14;
t46 = t99 * t82;
t35 = t46 * qJD(1);
t154 = -t82 * t35 + t45 * t74;
t122 = pkin(3) * t128;
t153 = -t104 * pkin(3) * t35 - t74 * t122;
t23 = t161 * t55 - t50;
t79 = t162 * t105;
t80 = t162 * t106;
t33 = -t104 * t79 + t161 * t80;
t152 = t99 * t121;
t151 = t105 * t84;
t34 = qJD(1) * t117 - t152;
t150 = t34 * qJ(5);
t149 = t74 * qJ(5);
t85 = qJD(1) * t94;
t107 = qJD(3) ^ 2;
t145 = t107 * t105;
t144 = t107 * t106;
t126 = t74 * pkin(4) + qJD(5);
t39 = t126 + t78;
t143 = qJD(5) + t39;
t142 = t105 ^ 2 - t106 ^ 2;
t139 = qJD(4) * t104;
t137 = qJD(1) * qJD(3);
t136 = t161 * pkin(3);
t135 = pkin(3) * t139;
t133 = pkin(3) * t141;
t52 = t161 * t56;
t108 = qJD(1) ^ 2;
t132 = t105 * t108 * t106;
t127 = t105 * t137;
t29 = pkin(3) * t127 + t35 * pkin(4);
t129 = qJD(3) * t162;
t95 = qJD(3) * t98;
t48 = -qJD(3) * t118 + t95;
t49 = t56 * qJD(3);
t124 = -t104 * t48 - t161 * t49;
t22 = -t104 * t55 - t52;
t32 = -t104 * t80 - t161 * t79;
t123 = t104 * t49 - t53 * t128 + t56 * t139 - t161 * t48;
t120 = t106 * t127;
t81 = -t131 + t146;
t119 = t147 * t46 - t81 * t34;
t116 = 0.2e1 * qJD(3) * t85;
t66 = t106 * t84 + t138;
t115 = t78 * t74 + t123;
t114 = t35 * qJ(5) + t123;
t19 = t104 * t53 + t52;
t71 = t105 * t129;
t72 = t106 * t129;
t12 = -t104 * t72 - t79 * t128 - t80 * t139 - t161 * t71;
t113 = -t99 * t130 + t152;
t112 = t143 * t74 + t114;
t8 = -t19 * qJD(4) + t124;
t13 = -t33 * qJD(4) + t104 * t71 - t161 * t72;
t58 = -t84 * t140 + t95;
t59 = t66 * qJD(3);
t65 = t98 - t151;
t111 = t59 * t105 + t58 * t106 + (-t105 * t66 - t106 * t65) * qJD(3);
t110 = t8 + t150;
t109 = (-t52 + (-pkin(3) * t99 - t53) * t104) * qJD(4) + t124;
t97 = t136 + pkin(4);
t86 = t99 * t122;
t73 = t74 ^ 2;
t57 = pkin(4) * t147 + t133;
t54 = t81 * pkin(4) + t83;
t38 = t46 * t99;
t31 = t46 * pkin(4) + t134;
t27 = -t73 + t163;
t25 = -t81 * qJ(5) + t33;
t24 = -t82 * qJ(5) + t32;
t20 = t113 + t159;
t17 = -t148 + t23;
t16 = t22 + t149;
t15 = t19 - t149;
t10 = t35 * t81 + t74 * t46;
t9 = -t147 * t45 - t34 * t82;
t6 = t45 * qJ(5) - t82 * qJD(5) + t13;
t5 = -t46 * qJ(5) - t81 * qJD(5) + t12;
t4 = -qJD(5) * t147 + t110;
t3 = -t74 * qJD(5) - t114;
t2 = t119 + t154;
t1 = -t119 + t154;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t120, -0.2e1 * t142 * t137, t144, -0.2e1 * t120, -t145, 0, t105 * t116 - t93 * t144, t106 * t116 + t93 * t145, t111, t111 * t93, t9, t1, -t160, t10, -t38, 0, t13 * t99 + t83 * t35 + t78 * t46 + (qJD(1) * t81 + t74) * t134, -t12 * t99 + t147 * t164 - t83 * t34 - t78 * t45, -t12 * t74 + t123 * t81 - t13 * t147 + t18 * t45 - t19 * t46 + t32 * t34 - t33 * t35 - t8 * t82, t19 * t12 - t123 * t33 + t18 * t13 + t78 * t164 + t8 * t32, t9, t1, -t160, t10, -t38, 0, t29 * t81 + t31 * t74 + t54 * t35 + t39 * t46 + t6 * t99, t147 * t31 + t29 * t82 - t54 * t34 - t39 * t45 - t5 * t99, t11 * t45 - t147 * t6 - t15 * t46 + t24 * t34 - t25 * t35 - t3 * t81 - t4 * t82 - t5 * t74, t11 * t6 + t15 * t5 + t4 * t24 + t3 * t25 + t29 * t54 + t39 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t144, 0, t58 * t105 - t59 * t106 + (-t105 * t65 + t106 * t66) * qJD(3), 0, 0, 0, 0, 0, 0, -t38, t160, t2, -t123 * t82 - t18 * t46 - t19 * t45 - t8 * t81, 0, 0, 0, 0, 0, 0, -t38, t160, t2, -t11 * t46 - t15 * t45 + t3 * t82 - t4 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, t142 * t108, 0, t132, 0, 0, -t85 * t141, -t85 * t106 * qJD(1) - t95 + (t65 + t151) * qJD(3), 0, 0, t158, t27, t20, -t158, 0, 0, -t74 * t133 - t22 * t99 + t109 - t156, -t133 * t147 + t23 * t99 + t115 - t86, t34 * t136 + (-t18 + t23) * t74 + (t19 + t22 + t135) * t147 + t153, -t18 * t22 - t19 * t23 + (-t78 * t141 + t161 * t8 - t104 * t123 + (-t104 * t18 + t161 * t19) * qJD(4)) * pkin(3), t158, t27, t20, -t158, 0, 0, -t143 * t147 - t16 * t99 - t57 * t74 + t109 + t150, -t147 * t57 + t17 * t99 + t112 - t86, t97 * t34 + (-t11 + t17) * t74 + (t15 + t16 + t135) * t147 + t153, -t11 * t16 - t15 * t17 - t39 * t57 + t4 * t97 + (t104 * t3 + (-t104 * t11 + t161 * t15) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t27, t20, -t158, 0, 0, t19 * t99 - t156 + t8, t18 * t99 + t115, 0, 0, t158, t27, t20, -t158, 0, 0, t15 * t99 + (-t126 - t39) * t147 + t110, -t163 * pkin(4) + t14 * t99 + t112, t34 * pkin(4) - t155 * t74, t155 * t15 + (-t147 * t39 + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99 * t147 + t35, t113 - t159, -t73 - t163, t11 * t147 + t15 * t74 + t29;];
tauc_reg = t7;
