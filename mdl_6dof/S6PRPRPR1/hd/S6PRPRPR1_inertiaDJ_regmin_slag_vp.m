% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:18
% EndTime: 2021-01-16 01:06:22
% DurationCPUTime: 0.80s
% Computational Cost: add. (896->136), mult. (2433->266), div. (0->0), fcn. (2456->12), ass. (0->96)
t50 = sin(pkin(6));
t53 = sin(qJ(4));
t49 = sin(pkin(11));
t51 = cos(pkin(11));
t54 = sin(qJ(2));
t57 = cos(qJ(2));
t71 = t49 * t54 - t51 * t57;
t64 = t71 * qJD(2);
t72 = t49 * t57 + t51 * t54;
t98 = cos(pkin(6));
t82 = t98 * qJD(4);
t56 = cos(qJ(4));
t92 = t56 * qJD(4);
t58 = -t53 * t82 + (t53 * t64 - t72 * t92) * t50;
t110 = 2 * qJD(6);
t86 = pkin(2) * t49 + pkin(8);
t81 = qJ(5) + t86;
t67 = qJD(4) * t81;
t27 = t56 * qJD(5) - t53 * t67;
t48 = sin(pkin(12));
t61 = -t53 * qJD(5) - t56 * t67;
t97 = cos(pkin(12));
t11 = t48 * t27 - t61 * t97;
t52 = sin(qJ(6));
t109 = t11 * t52;
t65 = t72 * t50;
t28 = qJD(2) * t65;
t30 = t71 * t50;
t108 = t30 * t28;
t37 = t48 * t56 + t53 * t97;
t32 = t37 * qJD(4);
t84 = t97 * t56;
t36 = t48 * t53 - t84;
t107 = t36 * t32;
t93 = t53 * qJD(4);
t33 = qJD(4) * t84 - t48 * t93;
t106 = t37 * t33;
t55 = cos(qJ(6));
t105 = t37 * t55;
t104 = t52 * t32;
t103 = t52 * t33;
t102 = t55 * t32;
t101 = t55 * t33;
t100 = t36 * t101 + t37 * t102;
t47 = t55 ^ 2;
t99 = t52 ^ 2 - t47;
t96 = qJD(2) * t50;
t95 = qJD(6) * t52;
t94 = qJD(6) * t55;
t43 = -pkin(4) * t97 - pkin(5);
t91 = t43 * t110;
t90 = 0.2e1 * t92;
t45 = pkin(4) * t93;
t89 = t37 * t95;
t88 = t30 * t93;
t87 = t52 * t94;
t44 = -pkin(2) * t51 - pkin(3);
t85 = -0.4e1 * t52 * t105;
t83 = t99 * qJD(6);
t24 = t53 * t98 + t56 * t65;
t63 = t53 * t65;
t60 = t56 * t98 - t63;
t8 = t24 * t97 + t48 * t60;
t79 = t30 * t55 - t52 * t8;
t78 = t30 * t52 + t55 * t8;
t77 = t86 * qJD(4);
t38 = -pkin(4) * t56 + t44;
t14 = pkin(5) * t36 - pkin(9) * t37 + t38;
t34 = t81 * t56;
t70 = t81 * t53;
t20 = t34 * t97 - t48 * t70;
t76 = t14 * t55 - t20 * t52;
t75 = t14 * t52 + t20 * t55;
t42 = pkin(4) * t48 + pkin(9);
t74 = -t32 * t42 + t33 * t43;
t73 = t36 * t42 - t37 * t43;
t62 = t50 * t64;
t10 = qJD(4) * t63 + (t62 - t82) * t56;
t5 = -t48 * t10 - t58 * t97;
t7 = t48 * t24 - t60 * t97;
t69 = t5 * t52 + t7 * t94;
t68 = -t5 * t55 + t7 * t95;
t17 = t36 * t94 + t104;
t66 = t37 * t94 + t103;
t16 = t89 - t101;
t35 = t37 ^ 2;
t19 = t48 * t34 + t70 * t97;
t15 = t36 * t95 - t102;
t13 = pkin(5) * t32 - pkin(9) * t33 + t45;
t12 = t27 * t97 + t48 * t61;
t6 = -t97 * t10 + t48 * t58;
t4 = -qJD(6) * t75 - t52 * t12 + t55 * t13;
t3 = -qJD(6) * t76 - t55 * t12 - t52 * t13;
t2 = -qJD(6) * t78 + t55 * t28 - t52 * t6;
t1 = -qJD(6) * t79 - t52 * t28 - t55 * t6;
t9 = [0, 0, 0, 0, -0.2e1 * t50 ^ 2 * t64 * t72 + 0.2e1 * t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t5 * t7 + 0.2e1 * t6 * t8 + 0.2e1 * t108, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t54 * t96, -t57 * t96, (-t28 * t51 - t49 * t62) * pkin(2), 0, 0, 0, 0, 0, -t28 * t56 + t88, t28 * t53 + t30 * t92, t28 * t36 + t30 * t32, t28 * t37 + t30 * t33, -t32 * t8 + t33 * t7 - t36 * t6 + t37 * t5, pkin(4) * t88 + t11 * t7 + t12 * t8 + t19 * t5 + t20 * t6 + t28 * t38, 0, 0, 0, 0, 0, t103 * t7 + t2 * t36 + t32 * t79 + t37 * t69, t1 * t36 + t101 * t7 - t32 * t78 - t37 * t68; 0, 0, 0, 0, 0, t53 * t90, 0.2e1 * (-t53 ^ 2 + t56 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * t44 * t93, t44 * t90, 0.2e1 * t32 * t38 + 0.2e1 * t36 * t45, 0.2e1 * t33 * t38 + 0.2e1 * t37 * t45, 0.2e1 * t11 * t37 - 0.2e1 * t12 * t36 + 0.2e1 * t19 * t33 - 0.2e1 * t20 * t32, 0.2e1 * t11 * t19 + 0.2e1 * t12 * t20 + 0.2e1 * t38 * t45, 0.2e1 * t106 * t47 - 0.2e1 * t35 * t87, t110 * t35 * t99 + t33 * t85, -0.2e1 * t36 * t89 + 0.2e1 * t100, -0.2e1 * t104 * t37 - 0.2e1 * t36 * t66, 0.2e1 * t107, 0.2e1 * t109 * t37 + 0.2e1 * t19 * t66 + 0.2e1 * t32 * t76 + 0.2e1 * t4 * t36, 0.2e1 * t105 * t11 - 0.2e1 * t16 * t19 + 0.2e1 * t3 * t36 - 0.2e1 * t32 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t7 + t33 * t8 + t36 * t5 + t37 * t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t36 + t12 * t37 + t19 * t32 + t20 * t33, 0, 0, 0, 0, 0, 0, (-t32 * t37 - t33 * t36) * t55 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t106 + 0.2e1 * t107, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t10, -t5, -t6, 0, (t48 * t6 - t5 * t97) * pkin(4), 0, 0, 0, 0, 0, t68, t69; 0, 0, 0, 0, 0, 0, 0, t92, -t93, 0, -t56 * t77, t53 * t77, -t11, -t12, (-t32 * t48 - t33 * t97) * pkin(4), (-t11 * t97 + t12 * t48) * pkin(4), t101 * t52 - t37 * t83, qJD(6) * t85 - t33 * t99, t17, -t15, 0, -t11 * t55 + t74 * t52 + (t19 * t52 - t55 * t73) * qJD(6), t109 + t74 * t55 + (t19 * t55 + t52 * t73) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t92, -t32, -t33, 0, (-t32 * t97 + t33 * t48) * pkin(4), 0, 0, 0, 0, 0, t15, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t87, -0.2e1 * t83, 0, 0, 0, t52 * t91, t55 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, 0, t45, 0, 0, 0, 0, 0, -t15, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t66, t32, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t95, 0, -t42 * t94, t42 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
