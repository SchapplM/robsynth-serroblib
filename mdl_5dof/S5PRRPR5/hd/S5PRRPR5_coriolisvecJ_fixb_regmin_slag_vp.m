% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:54
% EndTime: 2019-12-05 16:28:00
% DurationCPUTime: 1.01s
% Computational Cost: add. (1294->194), mult. (3525->302), div. (0->0), fcn. (2714->10), ass. (0->124)
t79 = sin(pkin(5));
t132 = qJD(1) * t79;
t87 = cos(qJ(2));
t118 = t87 * t132;
t102 = qJD(4) + t118;
t81 = cos(pkin(5));
t131 = qJD(1) * t81;
t83 = sin(qJ(3));
t117 = t83 * t131;
t125 = qJ(4) * qJD(3);
t84 = sin(qJ(2));
t119 = t84 * t132;
t63 = qJD(2) * pkin(7) + t119;
t86 = cos(qJ(3));
t158 = (-t86 * t63 - t117) * qJD(3) + (-t102 * t83 - t86 * t125) * qJD(2);
t130 = qJD(2) * t83;
t80 = cos(pkin(10));
t143 = t80 * t86;
t78 = sin(pkin(10));
t54 = qJD(2) * t143 - t78 * t130;
t51 = qJD(5) - t54;
t157 = -qJD(5) + t51;
t133 = qJD(3) * pkin(3);
t156 = t83 * t133 - t119;
t69 = t86 * t131;
t17 = (-t83 * t63 + t69) * qJD(3) + (t102 * t86 - t83 * t125) * qJD(2);
t3 = -t80 * t158 + t78 * t17;
t61 = t78 * t86 + t80 * t83;
t56 = t61 * qJD(2);
t72 = t78 * pkin(3) + pkin(8);
t155 = (pkin(3) * t130 + t56 * pkin(4) - t54 * pkin(8) + qJD(5) * t72) * t51 + t3;
t138 = -qJ(4) - pkin(7);
t115 = t138 * t83;
t65 = t138 * t86;
t40 = t78 * t115 - t80 * t65;
t55 = t61 * qJD(3);
t48 = qJD(2) * t55;
t103 = t3 * t61 - t40 * t48;
t111 = qJD(3) * t138;
t52 = t86 * qJD(4) + t83 * t111;
t60 = t78 * t83 - t143;
t94 = -t83 * qJD(4) + t86 * t111;
t136 = t60 * t118 + t80 * t52 + t78 * t94;
t116 = -t86 * pkin(3) - pkin(2);
t50 = t116 * qJD(2) + qJD(4) - t118;
t16 = -t54 * pkin(4) - t56 * pkin(8) + t50;
t28 = t60 * pkin(4) - t61 * pkin(8) + t116;
t4 = t158 * t78 + t80 * t17;
t108 = qJ(4) * qJD(2) + t63;
t36 = t108 * t86 + t117;
t147 = t78 * t36;
t35 = -t108 * t83 + t69;
t32 = t35 + t133;
t7 = t80 * t32 - t147;
t5 = -qJD(3) * pkin(4) - t7;
t58 = t60 * qJD(3);
t154 = -(qJD(5) * t16 + t4) * t60 - t5 * t58 + (-qJD(5) * t28 - t136) * t51 + t103;
t82 = sin(qJ(5));
t85 = cos(qJ(5));
t43 = t82 * qJD(3) + t85 * t56;
t124 = qJD(2) * qJD(3);
t113 = t86 * t124;
t114 = t83 * t124;
t49 = t80 * t113 - t78 * t114;
t19 = t43 * qJD(5) + t82 * t49;
t126 = t85 * qJD(3);
t127 = qJD(5) * t82;
t18 = qJD(5) * t126 - t56 * t127 + t85 * t49;
t153 = t18 * t82;
t152 = t28 * t48;
t41 = t82 * t56 - t126;
t151 = t41 * t51;
t150 = t43 * t51;
t149 = t43 * t56;
t148 = t56 * t41;
t146 = t79 * t84;
t145 = t79 * t87;
t89 = qJD(2) ^ 2;
t144 = t79 * t89;
t30 = t80 * t36;
t142 = t82 * t48;
t44 = t85 * t48;
t88 = qJD(3) ^ 2;
t140 = t88 * t83;
t139 = t88 * t86;
t137 = -t61 * t118 + t78 * t52 - t80 * t94;
t8 = t78 * t32 + t30;
t53 = pkin(3) * t114 + qJD(2) * t119;
t135 = t83 ^ 2 - t86 ^ 2;
t134 = qJD(2) * pkin(2);
t129 = qJD(2) * t84;
t128 = qJD(5) * t61;
t123 = t84 * t144;
t121 = t79 * t129;
t120 = qJD(2) * t145;
t110 = t51 * t85;
t107 = t83 * t120;
t106 = t86 * t120;
t105 = t55 * pkin(4) + t58 * pkin(8) + t156;
t6 = qJD(3) * pkin(8) + t8;
t1 = t85 * t16 - t82 * t6;
t2 = t82 * t16 + t85 * t6;
t101 = t44 + (t54 * t82 - t127) * t51;
t59 = t86 * t146 + t81 * t83;
t98 = -t83 * t146 + t81 * t86;
t25 = t80 * t59 + t78 * t98;
t100 = -t85 * t145 - t82 * t25;
t99 = t82 * t145 - t85 * t25;
t97 = -t61 * t127 - t85 * t58;
t96 = t134 * qJD(2);
t93 = -0.2e1 * qJD(3) * t134;
t12 = t80 * t35 - t147;
t92 = -t72 * t48 + (t12 + t5) * t51;
t73 = -t80 * pkin(3) - pkin(4);
t39 = -t80 * t115 - t78 * t65;
t34 = -t59 * qJD(3) - t107;
t33 = t98 * qJD(3) + t106;
t24 = t78 * t59 - t80 * t98;
t14 = t48 * pkin(4) - t49 * pkin(8) + t53;
t13 = t85 * t14;
t11 = t80 * t33 + t78 * t34;
t10 = t78 * t35 + t30;
t9 = t78 * t33 - t80 * t34;
t15 = [0, 0, -t123, -t87 * t144, 0, 0, 0, 0, 0, -t86 * t123 + (t34 - t107) * qJD(3), t83 * t123 + (-t33 - t106) * qJD(3), t11 * t54 + t24 * t49 - t25 * t48 + t9 * t56, t8 * t11 + t3 * t24 + t4 * t25 - t7 * t9 + (t129 * t50 - t53 * t87) * t79, 0, 0, 0, 0, 0, (qJD(5) * t99 - t82 * t11 + t121 * t85) * t51 + t100 * t48 + t9 * t41 + t24 * t19, -(qJD(5) * t100 + t85 * t11 + t121 * t82) * t51 + t99 * t48 + t9 * t43 + t24 * t18; 0, 0, 0, 0, 0.2e1 * t83 * t113, -0.2e1 * t135 * t124, t139, -t140, 0, -pkin(7) * t139 + t83 * t93, pkin(7) * t140 + t86 * t93, t136 * t54 + t137 * t56 + t39 * t49 - t4 * t60 - t8 * t55 + t7 * t58 + t103, t53 * t116 + t136 * t8 - t137 * t7 + t156 * t50 + t3 * t39 + t4 * t40, t18 * t85 * t61 + t43 * t97, -(-t41 * t85 - t43 * t82) * t58 + (-t153 - t19 * t85 + (t41 * t82 - t43 * t85) * qJD(5)) * t61, t18 * t60 + t43 * t55 + t44 * t61 + t51 * t97, -t61 * t142 - t19 * t60 - t41 * t55 + (-t128 * t85 + t82 * t58) * t51, t48 * t60 + t51 * t55, t1 * t55 + t13 * t60 + t39 * t19 + t137 * t41 + (t152 + t105 * t51 + (-t40 * t51 + t5 * t61 - t6 * t60) * qJD(5)) * t85 + t154 * t82, t39 * t18 - t2 * t55 + t137 * t43 + (-t152 - (-qJD(5) * t6 + t14) * t60 - t5 * t128 + (qJD(5) * t40 - t105) * t51) * t82 + t154 * t85; 0, 0, 0, 0, -t83 * t89 * t86, t135 * t89, 0, 0, 0, t83 * t96, t86 * t96, (-t10 + t8) * t56 + (-t12 + t7) * t54 + (-t48 * t78 - t49 * t80) * pkin(3), t7 * t10 - t8 * t12 + (-t130 * t50 - t3 * t80 + t4 * t78) * pkin(3), t110 * t43 + t153, (t18 - t151) * t85 + (-t19 - t150) * t82, t110 * t51 + t142 - t149, t101 + t148, -t51 * t56, -t1 * t56 - t10 * t41 - t155 * t85 + t73 * t19 + t92 * t82, -t10 * t43 + t155 * t82 + t73 * t18 + t2 * t56 + t92 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 ^ 2 - t56 ^ 2, -t8 * t54 + t7 * t56 + t53, 0, 0, 0, 0, 0, t101 - t148, -t51 ^ 2 * t85 - t142 - t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t41, -t41 ^ 2 + t43 ^ 2, t18 + t151, t150 - t19, t48, t157 * t2 - t82 * t4 - t5 * t43 + t13, t157 * t1 - t82 * t14 - t85 * t4 + t5 * t41;];
tauc_reg = t15;
