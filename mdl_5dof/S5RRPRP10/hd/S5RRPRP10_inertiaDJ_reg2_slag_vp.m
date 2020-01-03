% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:05
% EndTime: 2019-12-31 20:11:10
% DurationCPUTime: 1.21s
% Computational Cost: add. (857->152), mult. (1939->275), div. (0->0), fcn. (1353->4), ass. (0->104)
t116 = pkin(3) + pkin(6);
t117 = pkin(2) + pkin(7);
t107 = qJ(5) + t117;
t63 = sin(qJ(2));
t108 = t63 * qJ(3);
t65 = cos(qJ(2));
t122 = t117 * t65 + t108;
t110 = qJ(3) * t65;
t121 = t117 * t63 - t110;
t62 = sin(qJ(4));
t58 = t62 ^ 2;
t64 = cos(qJ(4));
t60 = t64 ^ 2;
t113 = t58 - t60;
t61 = t65 ^ 2;
t82 = qJD(2) * (t63 ^ 2 - t61);
t120 = t113 * qJD(4);
t42 = t116 * t65;
t96 = t65 * qJD(3);
t119 = t122 * qJD(2) - qJD(4) * t42 - t96;
t118 = 0.2e1 * qJD(3);
t115 = pkin(4) * t64;
t114 = t62 * t63;
t31 = -pkin(1) - t122;
t90 = t116 * t63;
t10 = t64 * t31 + t62 * t90;
t111 = qJ(3) * t62;
t109 = qJ(5) * t65;
t106 = qJD(2) * t62;
t105 = qJD(2) * t64;
t104 = qJD(4) * t62;
t103 = qJD(4) * t64;
t102 = qJD(4) * t65;
t101 = qJD(4) * t117;
t100 = t62 * qJD(5);
t99 = t63 * qJD(2);
t98 = t63 * qJD(3);
t97 = t64 * qJD(5);
t54 = t65 * qJD(2);
t95 = qJ(3) * qJD(4);
t94 = -0.2e1 * pkin(1) * qJD(2);
t93 = pkin(4) * t104;
t92 = pkin(6) * t99;
t91 = pkin(6) * t54;
t89 = t62 * t102;
t88 = t64 * t102;
t87 = t62 * t54;
t86 = t63 * t54;
t85 = t64 * t99;
t84 = t62 * t103;
t83 = qJ(5) * t102;
t39 = t107 * t64;
t36 = t64 * t90;
t81 = t62 * t85;
t80 = t61 * t84;
t79 = t116 * t54;
t12 = -pkin(4) * t89 + (-t115 - t116) * t99;
t53 = t62 * pkin(4) + qJ(3);
t78 = -t53 * t102 + t12;
t77 = -t31 * t103 + t64 * t79;
t7 = t63 * pkin(4) + t36 + (-t31 + t109) * t62;
t8 = -t64 * t109 + t10;
t76 = -t62 * t7 + t64 * t8;
t9 = -t62 * t31 + t36;
t75 = t10 * t64 - t62 * t9;
t74 = -t65 * pkin(2) - t108;
t72 = -qJD(4) * t116 + qJD(3);
t5 = t31 * t104 - t64 * (t121 * qJD(2) - t98) - t62 * t79 - qJD(4) * t36;
t28 = t62 * t99 - t88;
t37 = t116 * t99;
t70 = t121 * qJD(4) - t37;
t24 = t65 * t115 + t42;
t50 = pkin(4) * t103 + qJD(3);
t69 = -qJD(4) * t24 - t50 * t65 + t53 * t99;
t68 = t74 * qJD(2) + t96;
t6 = qJ(3) * t87 + (-t117 * qJD(2) + t72) * t114 + t77;
t1 = t75 * qJD(4) - t5 * t62 + t6 * t64;
t67 = t64 * t83 + t65 * t100 + (-t107 * qJD(2) + t72) * t114 + t77;
t57 = qJ(3) * t118;
t49 = -0.2e1 * t84;
t48 = 0.2e1 * t84;
t47 = -0.2e1 * t86;
t46 = 0.2e1 * t86;
t40 = -pkin(1) + t74;
t38 = t107 * t62;
t34 = -0.2e1 * t82;
t33 = 0.2e1 * t120;
t27 = -t63 * t103 - t87;
t26 = t85 + t89;
t25 = -t63 * t104 + t64 * t54;
t23 = -t98 + (pkin(2) * t63 - t110) * qJD(2);
t21 = -qJD(4) * t39 - t100;
t20 = t107 * t104 - t97;
t19 = -0.2e1 * t60 * t86 - 0.2e1 * t80;
t18 = -0.2e1 * t58 * t86 + 0.2e1 * t80;
t17 = -t113 * t102 - t81;
t16 = -t113 * t99 + 0.4e1 * t65 * t84;
t14 = 0.2e1 * t62 * t82 - 0.2e1 * t63 * t88;
t13 = 0.2e1 * t63 * t89 + 0.2e1 * t64 * t82;
t11 = -0.2e1 * t61 * t120 - 0.4e1 * t65 * t81;
t4 = t20 * t64 + t21 * t62 + (-t38 * t64 + t39 * t62) * qJD(4);
t3 = -qJ(5) * t85 - t62 * t83 + t65 * t97 + t5;
t2 = (pkin(4) + t111) * t54 + t67;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t34, 0, t47, 0, 0, t63 * t94, t65 * t94, 0, 0, 0, 0, 0, t46, t34, t47, 0, 0.2e1 * t23 * t65 - 0.2e1 * t40 * t99, -0.2e1 * t23 * t63 - 0.2e1 * t40 * t54, 0.2e1 * t40 * t23, t18, t11, t14, t19, t13, t46, 0.2e1 * (-t42 * t105 + t6) * t63 + 0.2e1 * (qJD(2) * t9 - t42 * t104 - t37 * t64) * t65, 0.2e1 * (t42 * t106 + t5) * t63 + 0.2e1 * (-qJD(2) * t10 - t42 * t103 + t37 * t62) * t65, 0.2e1 * t75 * t99 + 0.2e1 * (t5 * t64 + t6 * t62 + (t10 * t62 + t64 * t9) * qJD(4)) * t65, -0.2e1 * t10 * t5 - 0.2e1 * t42 * t37 + 0.2e1 * t9 * t6, t18, t11, t14, t19, t13, t46, 0.2e1 * (-t24 * t105 + t2) * t63 + 0.2e1 * (qJD(2) * t7 - t24 * t104 + t12 * t64) * t65, 0.2e1 * (t24 * t106 + t3) * t63 + 0.2e1 * (-qJD(2) * t8 - t24 * t103 - t12 * t62) * t65, 0.2e1 * t76 * t99 + 0.2e1 * (t2 * t62 + t3 * t64 + (t62 * t8 + t64 * t7) * qJD(4)) * t65, 0.2e1 * t24 * t12 + 0.2e1 * t7 * t2 - 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t99, 0, -t91, t92, 0, 0, 0, -t54, t99, 0, 0, 0, t68, t91, -t92, t68 * pkin(6), -t17, t16, t25, t17, t27, 0, -t119 * t64 + t70 * t62, t119 * t62 + t70 * t64, -t1, -t37 * qJ(3) + t42 * qJD(3) - t1 * t117, -t17, t16, t25, t17, t27, 0, t20 * t63 - t39 * t54 + t78 * t62 - t69 * t64, -t21 * t63 + t38 * t54 + t69 * t62 + t78 * t64, (-t38 * t99 - t21 * t65 - t2 + (-t39 * t65 - t8) * qJD(4)) * t64 + (t39 * t99 + t20 * t65 + t3 + (-t38 * t65 + t7) * qJD(4)) * t62, t12 * t53 - t2 * t39 + t7 * t20 + t8 * t21 + t24 * t50 + t3 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t57, t49, t33, 0, t48, 0, 0, 0.2e1 * qJD(3) * t62 + 0.2e1 * t64 * t95, 0.2e1 * qJD(3) * t64 - 0.2e1 * t62 * t95, 0, t57, t49, t33, 0, t48, 0, 0, 0.2e1 * t53 * t103 + 0.2e1 * t50 * t62, -0.2e1 * t53 * t104 + 0.2e1 * t50 * t64, -0.2e1 * t4, -0.2e1 * t39 * t20 - 0.2e1 * t38 * t21 + 0.2e1 * t53 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, t91, 0, 0, 0, 0, 0, 0, t25, t27, 0, t1, 0, 0, 0, 0, 0, 0, t25, t27, 0, t76 * qJD(4) + t2 * t64 - t3 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t26, t54, t6, t5, 0, 0, 0, 0, t28, 0, t26, t54, (0.2e1 * pkin(4) + t111) * t54 + t67, t3, -t28 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, 0, -t103, 0, t62 * t101, t64 * t101, 0, 0, 0, 0, -t104, 0, -t103, 0, t20, -t21, t93, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t103, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t103, 0, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t28, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, -t104, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
