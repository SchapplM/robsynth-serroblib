% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:09
% EndTime: 2019-12-05 16:10:14
% DurationCPUTime: 0.94s
% Computational Cost: add. (1121->186), mult. (2998->252), div. (0->0), fcn. (2061->6), ass. (0->119)
t128 = cos(pkin(8));
t136 = -qJ(4) - pkin(6);
t110 = qJD(3) * t136;
t80 = sin(qJ(3));
t82 = cos(qJ(3));
t47 = t82 * qJD(4) + t80 * t110;
t79 = sin(pkin(8));
t83 = cos(qJ(2));
t112 = t128 * t82;
t141 = t79 * t80;
t97 = t112 - t141;
t92 = t97 * t83;
t94 = -t80 * qJD(4) + t82 * t110;
t134 = -qJD(1) * t92 + t128 * t47 + t79 * t94;
t119 = (qJD(2) * qJD(3));
t154 = -2 * t119;
t81 = sin(qJ(2));
t84 = qJD(3) ^ 2;
t85 = qJD(2) ^ 2;
t153 = (t84 + t85) * t81;
t113 = t128 * t80;
t140 = t79 * t82;
t60 = t113 + t140;
t150 = t60 * qJD(2);
t148 = t150 ^ 2;
t108 = qJD(2) * t128;
t102 = t82 * t108;
t127 = qJD(2) * t80;
t51 = t79 * t127 - t102;
t48 = t51 ^ 2;
t152 = -t48 - t148;
t151 = -t48 + t148;
t53 = t60 * qJD(3);
t124 = qJD(3) * t80;
t56 = qJD(3) * t112 - t79 * t124;
t123 = t81 * qJD(1);
t149 = pkin(3) * t124 - t123;
t147 = pkin(3) * t80;
t122 = t83 * qJD(1);
t105 = qJD(4) + t122;
t120 = qJ(4) * qJD(3);
t129 = qJD(2) * pkin(6);
t67 = t123 + t129;
t28 = -t67 * t124 + (t105 * t82 - t80 * t120) * qJD(2);
t86 = -t82 * t67 * qJD(3) + (-t105 * t80 - t82 * t120) * qJD(2);
t2 = -t128 * t86 + t79 * t28;
t65 = t136 * t82;
t30 = -t136 * t113 - t79 * t65;
t146 = t2 * t30;
t42 = t60 * t81;
t145 = t2 * t42;
t76 = -t82 * pkin(3) - pkin(2);
t57 = t76 * qJD(2) + qJD(4) - t122;
t10 = t51 * pkin(4) - qJ(5) * t150 + t57;
t144 = t10 * t150;
t143 = t150 * t51;
t107 = qJ(4) * qJD(2) + t67;
t45 = t107 * t82;
t142 = t79 * t45;
t139 = t84 * t80;
t138 = t84 * t82;
t137 = t85 * t83;
t135 = -t60 * t122 - t128 * t94 + t79 * t47;
t3 = t128 * t28 + t79 * t86;
t33 = t128 * t45;
t44 = t107 * t80;
t35 = qJD(3) * pkin(3) - t44;
t13 = t79 * t35 + t33;
t114 = t80 * t119;
t58 = pkin(3) * t114 + qJD(2) * t123;
t77 = t80 ^ 2;
t78 = t82 ^ 2;
t133 = t77 - t78;
t132 = t77 + t78;
t130 = qJD(2) * pkin(2);
t126 = qJD(2) * t81;
t125 = qJD(3) * t53;
t16 = -t128 * t44 - t142;
t121 = qJD(5) - t16;
t118 = t80 * t85 * t82;
t116 = pkin(3) * t127;
t111 = -t53 * pkin(4) + t56 * qJ(5) + t60 * qJD(5) - t149;
t68 = -t122 - t130;
t109 = -t68 - t122;
t106 = t83 * t154;
t104 = t82 * t114;
t40 = qJD(2) * t53;
t101 = -t40 * t97 + t51 * t53;
t100 = qJD(2) * t109;
t15 = -t79 * t44 + t33;
t99 = t15 * qJD(3) - t2;
t17 = (-qJD(2) * t140 - t108 * t80) * t83 - t56 * t81;
t98 = t17 * qJD(3) + t51 * t126 - t83 * t40;
t12 = t128 * t35 - t142;
t18 = qJD(2) * t92 - t81 * t53;
t66 = t79 * t114;
t41 = qJD(3) * t102 - t66;
t43 = t97 * t81;
t96 = -t150 * t17 - t18 * t51 - t43 * t40 + t42 * t41;
t95 = t40 * pkin(4) - t41 * qJ(5) + t58;
t93 = qJD(3) * (-t109 - t130);
t90 = t18 * qJD(3) - t126 * t150 + t83 * t41;
t89 = t150 * t53 + t60 * t40 - t41 * t97 + t56 * t51;
t88 = 0.2e1 * t150 * qJD(3);
t31 = -t128 * t65 + t136 * t141;
t87 = -t134 * t51 + t135 * t150 + t2 * t60 + t30 * t41 - t31 * t40;
t74 = -t128 * pkin(3) - pkin(4);
t72 = t79 * pkin(3) + qJ(5);
t46 = t56 * qJD(3);
t29 = -pkin(4) * t97 - t60 * qJ(5) + t76;
t27 = -t66 + (t102 + t51) * qJD(3);
t26 = -t66 + (t102 - t51) * qJD(3);
t19 = pkin(4) * t150 + t51 * qJ(5) + t116;
t11 = qJD(3) * qJ(5) + t13;
t9 = -qJD(3) * pkin(4) + qJD(5) - t12;
t5 = t150 * t56 + t41 * t60;
t4 = -qJD(5) * t150 + t95;
t1 = qJD(3) * qJD(5) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85 * t81, -t137, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t106 - t82 * t153, t82 * t106 + t80 * t153, t132 * t137, (t68 * t81 + (-t123 + (t67 + t123) * t132) * t83) * qJD(2), 0, 0, 0, 0, 0, 0, t98, -t90, t96, t12 * t17 + t57 * t126 + t13 * t18 + t3 * t43 - t58 * t83 + t145, 0, 0, 0, 0, 0, 0, t98, t96, t90, t1 * t43 + t10 * t126 + t11 * t18 - t9 * t17 - t4 * t83 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t104, t133 * t154, t138, -0.2e1 * t104, -t139, 0, -pkin(6) * t138 + t80 * t93, pkin(6) * t139 + t82 * t93, 0, ((-t68 - t130) * t81 + (t129 - t67) * t83 * t132) * qJD(1), t5, -t89, t46, t101, -t125, 0, -t51 * t123 + t76 * t40 + t57 * t53 - t58 * t97 + (t51 * t147 - t135) * qJD(3), -t150 * t123 + t76 * t41 + t57 * t56 + t58 * t60 + (t147 * t150 - t134) * qJD(3), -t12 * t56 - t13 * t53 + t3 * t97 + t87, -t135 * t12 + t134 * t13 + t149 * t57 + t3 * t31 + t58 * t76 + t146, t5, t46, t89, 0, t125, t101, -t135 * qJD(3) + t10 * t53 - t111 * t51 + t29 * t40 - t4 * t97, t1 * t97 - t11 * t53 + t9 * t56 + t87, t134 * qJD(3) - t10 * t56 + t111 * t150 - t29 * t41 - t4 * t60, t1 * t31 - t111 * t10 + t134 * t11 + t135 * t9 + t4 * t29 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t133 * t85, 0, t118, 0, 0, t80 * t100, t82 * t100, 0, 0, t143, t151, t27, -t143, 0, 0, -t51 * t116 - t150 * t57 + t99, t16 * qJD(3) - t116 * t150 + t57 * t51 - t3, (t13 - t15) * t150 + (-t12 + t16) * t51 + (-t128 * t41 - t40 * t79) * pkin(3), t12 * t15 - t13 * t16 + (-t57 * t127 - t128 * t2 + t3 * t79) * pkin(3), t143, t27, -t151, 0, 0, -t143, -t19 * t51 - t144 + t99, -t72 * t40 + t74 * t41 + (t11 - t15) * t150 + (t9 - t121) * t51, -t10 * t51 + t19 * t150 + (0.2e1 * qJD(5) - t16) * qJD(3) + t3, t1 * t72 - t10 * t19 + t121 * t11 - t9 * t15 + t2 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t26, t152, t12 * t150 + t13 * t51 + t58, 0, 0, 0, 0, 0, 0, t88, t152, -t26, t11 * t51 + (-qJD(5) - t9) * t150 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t27, -t148 - t84, -t11 * qJD(3) + t144 + t2;];
tauc_reg = t6;
