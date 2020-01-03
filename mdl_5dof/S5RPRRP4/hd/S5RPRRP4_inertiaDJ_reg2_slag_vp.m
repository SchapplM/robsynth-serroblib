% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP4
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
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:25
% EndTime: 2020-01-03 11:50:28
% DurationCPUTime: 0.91s
% Computational Cost: add. (1332->132), mult. (3266->231), div. (0->0), fcn. (2890->6), ass. (0->92)
t71 = sin(qJ(4));
t72 = sin(qJ(3));
t73 = cos(qJ(4));
t74 = cos(qJ(3));
t49 = t71 * t74 + t73 * t72;
t69 = sin(pkin(8));
t40 = t49 * t69;
t111 = t69 * t72;
t70 = cos(pkin(8));
t57 = t74 * t70 * qJ(2);
t84 = t70 * pkin(2) + pkin(1);
t78 = -t69 * pkin(6) - t84;
t38 = t72 * t78 + t57;
t35 = -pkin(7) * t111 + t38;
t105 = qJ(2) * t72;
t46 = t74 * t78;
t95 = t74 * t69 * pkin(7);
t77 = -t95 + t46 + (-pkin(3) - t105) * t70;
t11 = t73 * t35 + t71 * t77;
t96 = qJD(3) + qJD(4);
t104 = qJD(2) * t74;
t106 = qJD(3) * t46 + t70 * t104;
t97 = qJ(2) * qJD(3);
t83 = t72 * t97;
t30 = t70 * t83 - t106;
t99 = t72 * qJD(2);
t85 = t70 * t99;
t31 = -t38 * qJD(3) - t85;
t90 = t70 * t105;
t37 = t46 - t90;
t117 = t30 * t72 - t31 * t74 + (t37 * t72 - t38 * t74) * qJD(3);
t116 = 0.2e1 * t70;
t115 = 0.2e1 * qJD(3);
t114 = pkin(3) * t71;
t113 = t73 * pkin(3);
t34 = t96 * t49;
t112 = t34 * t70;
t110 = t71 * t72;
t109 = t73 * t74;
t103 = qJD(3) * t72;
t89 = t69 * t103;
t93 = t69 * t109;
t94 = t69 * t110;
t23 = -qJD(4) * t94 - t71 * t89 + t96 * t93;
t100 = qJD(4) * t73;
t91 = pkin(3) * t100;
t108 = -t23 * t114 - t40 * t91;
t102 = qJD(3) * t74;
t33 = -t74 * t100 - t73 * t102 + t96 * t110;
t107 = -t33 * t114 + t49 * t91;
t88 = t69 * t102;
t44 = pkin(3) * t88 + t69 * qJD(2);
t47 = pkin(3) * t111 + t69 * qJ(2);
t101 = qJD(4) * t71;
t98 = qJ(2) * qJD(2);
t92 = pkin(3) * t101;
t41 = t93 - t94;
t87 = t41 * t101;
t48 = t109 - t110;
t86 = t48 * t101;
t26 = t73 * t77;
t10 = -t71 * t35 + t26;
t82 = t69 * t70 * t115;
t67 = t69 ^ 2;
t81 = t67 * t72 * t102;
t75 = -t85 + (-t57 + ((pkin(6) + pkin(7)) * t69 + t84) * t72) * qJD(3);
t76 = (-t90 - t95) * qJD(3) + t106;
t4 = -qJD(4) * t26 + t35 * t101 - t71 * t75 - t73 * t76;
t1 = t23 * qJ(5) + t40 * qJD(5) + t4;
t5 = -t11 * qJD(4) - t71 * t76 + t73 * t75;
t22 = t96 * t40;
t2 = t22 * qJ(5) - t41 * qJD(5) + t5;
t68 = t70 ^ 2;
t64 = pkin(4) + t113;
t63 = t67 * t98;
t62 = -0.2e1 * t91;
t61 = -0.2e1 * t92;
t54 = t70 * t91;
t53 = t70 * t92;
t29 = t40 * pkin(4) + t47;
t28 = t33 * t70;
t16 = t23 * t116;
t15 = t22 * t116;
t14 = t23 * pkin(4) + t44;
t13 = -0.2e1 * t41 * t22;
t12 = 0.2e1 * t40 * t23;
t9 = -0.2e1 * t49 * t33 - 0.2e1 * t48 * t34;
t8 = -t40 * qJ(5) + t11;
t7 = -t70 * pkin(4) - t41 * qJ(5) + t10;
t6 = 0.2e1 * t22 * t40 - 0.2e1 * t41 * t23;
t3 = t48 * t22 - t49 * t23 + t33 * t40 + t34 * t41;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t67 + t68) * qJD(2), 0.2e1 * t68 * t98 + 0.2e1 * t63, -0.2e1 * t81, (t72 ^ 2 - t74 ^ 2) * t67 * t115, t72 * t82, 0.2e1 * t81, t74 * t82, 0, -0.2e1 * t31 * t70 + 0.2e1 * (t74 * t97 + t99) * t67, -0.2e1 * t30 * t70 + 0.2e1 * (-t83 + t104) * t67, 0.2e1 * t117 * t69, -0.2e1 * t38 * t30 + 0.2e1 * t37 * t31 + 0.2e1 * t63, t13, t6, t15, t12, t16, 0, 0.2e1 * t47 * t23 + 0.2e1 * t44 * t40 - 0.2e1 * t5 * t70, -0.2e1 * t47 * t22 - 0.2e1 * t4 * t70 + 0.2e1 * t44 * t41, 0.2e1 * t10 * t22 - 0.2e1 * t11 * t23 + 0.2e1 * t4 * t40 - 0.2e1 * t5 * t41, 0.2e1 * t10 * t5 - 0.2e1 * t11 * t4 + 0.2e1 * t47 * t44, t13, t6, t15, t12, t16, 0, 0.2e1 * t14 * t40 - 0.2e1 * t2 * t70 + 0.2e1 * t29 * t23, -0.2e1 * t1 * t70 + 0.2e1 * t14 * t41 - 0.2e1 * t29 * t22, 0.2e1 * t1 * t40 - 0.2e1 * t2 * t41 + 0.2e1 * t7 * t22 - 0.2e1 * t8 * t23, -0.2e1 * t8 * t1 + 0.2e1 * t29 * t14 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t103, t70 * t102, 0, -t117, 0, 0, 0, 0, 0, 0, t112, -t28, t3, -t10 * t34 - t11 * t33 - t4 * t49 + t5 * t48, 0, 0, 0, 0, 0, 0, t112, -t28, t3, -t1 * t49 + t2 * t48 - t8 * t33 - t7 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, 0, -t88, 0, t31, t30, 0, 0, 0, 0, -t22, 0, -t23, 0, t53 + t5, t4 + t54, (t22 * t73 + t87) * pkin(3) + t108, (-t4 * t71 + t5 * t73 + (-t10 * t71 + t11 * t73) * qJD(4)) * pkin(3), 0, 0, -t22, 0, -t23, 0, t53 + t2, t1 + t54, pkin(3) * t87 + t64 * t22 + t108, t2 * t64 + (-t1 * t71 + (-t7 * t71 + t73 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t102, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, (-t34 * t73 - t86) * pkin(3) + t107, 0, 0, 0, 0, 0, 0, -t34, t33, 0, -pkin(3) * t86 - t34 * t64 + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t62, 0, 0, 0, 0, 0, 0, 0, 0, t61, t62, 0, 0.2e1 * (-t64 + t113) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, -t23, 0, t5, t4, 0, 0, 0, 0, -t22, 0, -t23, 0, t2, t1, t22 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, -t34 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t91, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t91, 0, -pkin(4) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
