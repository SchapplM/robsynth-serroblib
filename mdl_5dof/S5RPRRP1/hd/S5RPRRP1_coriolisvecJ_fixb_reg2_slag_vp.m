% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:58
% EndTime: 2019-12-05 18:00:03
% DurationCPUTime: 1.05s
% Computational Cost: add. (1869->191), mult. (3971->242), div. (0->0), fcn. (2458->4), ass. (0->122)
t97 = qJD(3) + qJD(4);
t104 = cos(qJ(4));
t102 = sin(qJ(4));
t103 = sin(qJ(3));
t137 = qJD(1) * t103;
t106 = -pkin(1) - pkin(6);
t83 = qJD(1) * t106 + qJD(2);
t56 = -pkin(7) * t137 + t103 * t83;
t49 = t102 * t56;
t105 = cos(qJ(3));
t136 = qJD(1) * t105;
t57 = -pkin(7) * t136 + t105 * t83;
t52 = qJD(3) * pkin(3) + t57;
t24 = t104 * t52 - t49;
t124 = t104 * t136;
t125 = t102 * t137;
t65 = t124 - t125;
t58 = t65 * qJ(5);
t14 = t24 - t58;
t73 = t102 * t105 + t103 * t104;
t63 = t73 * qJD(1);
t158 = t65 ^ 2;
t98 = qJD(1) * qJD(2);
t157 = 0.2e1 * t98;
t156 = t104 * pkin(3);
t145 = t104 * t105;
t117 = t97 * t145;
t133 = qJD(4) * t102;
t135 = qJD(3) * t103;
t37 = -t102 * t135 - t103 * t133 + t117;
t35 = t37 * t97;
t154 = t65 * t63;
t153 = t65 * t97;
t79 = pkin(3) * t137 + qJD(1) * qJ(2);
t152 = t79 * t65;
t151 = pkin(7) - t106;
t11 = pkin(4) * t97 + t14;
t150 = t11 - t14;
t132 = qJD(4) * t104;
t128 = pkin(3) * t132;
t148 = t97 * t125;
t33 = qJD(1) * t117 - t148;
t149 = -pkin(3) * t102 * t33 - t128 * t63;
t29 = t104 * t57 - t49;
t77 = t151 * t103;
t78 = t151 * t105;
t39 = -t102 * t78 - t104 * t77;
t131 = qJD(1) * qJD(3);
t123 = t105 * t131;
t76 = pkin(3) * t123 + t98;
t50 = t104 * t56;
t36 = t97 * t73;
t32 = t36 * qJD(1);
t147 = t32 * qJ(5);
t146 = t63 * qJ(5);
t89 = pkin(3) * t103 + qJ(2);
t107 = qJD(3) ^ 2;
t144 = t107 * t103;
t143 = t107 * t105;
t108 = qJD(1) ^ 2;
t142 = t108 * qJ(2);
t141 = t108 * t105;
t134 = qJD(3) * t105;
t84 = pkin(3) * t134 + qJD(2);
t122 = -pkin(4) * t63 - qJD(5);
t40 = -t122 + t79;
t140 = qJD(5) + t40;
t139 = t103 ^ 2 - t105 ^ 2;
t138 = -t107 - t108;
t130 = 0.2e1 * qJD(1);
t129 = pkin(3) * t133;
t127 = pkin(3) * t136;
t126 = t103 * t141;
t121 = pkin(7) * qJD(1) - t83;
t53 = t121 * t135;
t54 = t121 * t134;
t120 = t102 * t54 + t104 * t53;
t28 = -t102 * t57 - t50;
t38 = t102 * t77 - t104 * t78;
t119 = -t102 * t53 + t104 * t54 - t132 * t52 + t133 * t56;
t22 = pkin(4) * t33 + t76;
t118 = t103 * t123;
t74 = -t102 * t103 + t145;
t9 = -t32 * t74 - t36 * t65;
t10 = t33 * t73 + t37 * t63;
t25 = t102 * t52 + t50;
t116 = t63 * t79 + t119;
t115 = t33 * qJ(5) + t119;
t71 = t151 * t135;
t72 = qJD(3) * t78;
t12 = t102 * t71 - t104 * t72 - t132 * t78 + t133 * t77;
t15 = t25 - t146;
t3 = -qJD(5) * t63 - t115;
t8 = -qJD(4) * t25 + t120;
t111 = t8 + t147;
t4 = -t65 * qJD(5) + t111;
t114 = -t11 * t36 + t15 * t37 + t3 * t73 + t4 * t74;
t113 = -t119 * t73 - t24 * t36 + t25 * t37 + t74 * t8;
t112 = t140 * t63 + t115;
t13 = -qJD(4) * t39 + t102 * t72 + t104 * t71;
t110 = (-t50 + (-pkin(3) * t97 - t52) * t102) * qJD(4) + t120;
t95 = qJ(2) * t157;
t91 = pkin(4) + t156;
t80 = t97 * t128;
t62 = t63 ^ 2;
t55 = pkin(4) * t73 + t89;
t46 = pkin(4) * t65 + t127;
t34 = t36 * t97;
t31 = pkin(4) * t37 + t84;
t27 = -qJ(5) * t73 + t39;
t26 = -qJ(5) * t74 + t38;
t23 = -t62 + t158;
t21 = -qJD(1) * t65 - t35;
t20 = -qJD(1) * t63 - t34;
t19 = -t124 * t97 + t148 + t153;
t17 = -t58 + t29;
t16 = t28 + t146;
t6 = t36 * qJ(5) - t74 * qJD(5) + t13;
t5 = -qJ(5) * t37 - qJD(5) * t73 + t12;
t2 = -t10 - t9;
t1 = t32 * t73 - t33 * t74 + t36 * t63 - t37 * t65;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t95, -0.2e1 * t118, 0.2e1 * t139 * t131, -t144, 0.2e1 * t118, -t143, 0, -t106 * t144 + (qJ(2) * t134 + qJD(2) * t103) * t130, -t106 * t143 + (-qJ(2) * t135 + qJD(2) * t105) * t130, 0, t95, t9, t1, -t34, t10, -t35, 0, t13 * t97 + t33 * t89 + t37 * t79 + t63 * t84 + t73 * t76, -t12 * t97 - t32 * t89 - t36 * t79 + t65 * t84 + t74 * t76, -t12 * t63 - t13 * t65 + t32 * t38 - t33 * t39 - t113, -t119 * t39 + t12 * t25 + t13 * t24 + t38 * t8 + t76 * t89 + t79 * t84, t9, t1, -t34, t10, -t35, 0, t22 * t73 + t31 * t63 + t33 * t55 + t37 * t40 + t6 * t97, t22 * t74 + t31 * t65 - t32 * t55 - t36 * t40 - t5 * t97, t26 * t32 - t27 * t33 - t5 * t63 - t6 * t65 - t114, t11 * t6 + t15 * t5 + t22 * t55 + t26 * t4 + t27 * t3 + t31 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t142, 0, 0, 0, 0, 0, 0, t138 * t103, t138 * t105, 0, -t142, 0, 0, 0, 0, 0, 0, t20, t21, t2, -qJD(1) * t79 + t113, 0, 0, 0, 0, 0, 0, t20, t21, t2, -qJD(1) * t40 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t139 * t108, 0, -t126, 0, 0, -qJ(2) * t141, t103 * t142, 0, 0, t154, t23, 0, -t154, t19, 0, -t127 * t63 - t28 * t97 + t110 - t152, -t127 * t65 + t29 * t97 + t116 - t80, t32 * t156 + (-t24 + t29) * t63 + (t25 + t28 + t129) * t65 + t149, -t24 * t28 - t25 * t29 + (-t79 * t136 - t102 * t119 + t104 * t8 + (-t102 * t24 + t104 * t25) * qJD(4)) * pkin(3), t154, t23, 0, -t154, t19, 0, -t140 * t65 - t16 * t97 - t46 * t63 + t110 + t147, t17 * t97 - t46 * t65 + t112 - t80, t91 * t32 + (-t11 + t17) * t63 + (t15 + t16 + t129) * t65 + t149, -t11 * t16 - t15 * t17 + t4 * t91 - t40 * t46 + (t102 * t3 + (-t102 * t11 + t104 * t15) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t23, 0, -t154, t19, 0, t25 * t97 - t152 + t8, t24 * t97 + t116, 0, 0, t154, t23, 0, -t154, t19, 0, t15 * t97 + (t122 - t40) * t65 + t111, -pkin(4) * t158 + t14 * t97 + t112, t32 * pkin(4) - t150 * t63, t150 * t15 + (-t40 * t65 + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 + t153, -0.2e1 * t63 * t97, -t62 - t158, t11 * t65 + t15 * t63 + t22;];
tauc_reg = t7;
