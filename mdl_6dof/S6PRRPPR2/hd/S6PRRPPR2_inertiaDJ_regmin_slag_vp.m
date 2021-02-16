% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:18
% EndTime: 2021-01-16 02:21:22
% DurationCPUTime: 0.86s
% Computational Cost: add. (894->142), mult. (2371->268), div. (0->0), fcn. (2248->10), ass. (0->94)
t48 = sin(pkin(11));
t50 = cos(pkin(11));
t55 = cos(qJ(3));
t87 = t55 * qJD(3);
t52 = sin(qJ(3));
t88 = t52 * qJD(3);
t29 = -t48 * t88 + t50 * t87;
t33 = t48 * t55 + t50 * t52;
t49 = sin(pkin(6));
t56 = cos(qJ(2));
t53 = sin(qJ(2));
t92 = qJD(2) * t53;
t105 = (-t29 * t56 + t33 * t92) * t49;
t28 = t33 * qJD(3);
t32 = t48 * t52 - t50 * t55;
t104 = (-t28 * t56 + t32 * t92) * t49;
t51 = sin(qJ(6));
t46 = t51 ^ 2;
t54 = cos(qJ(6));
t94 = -t54 ^ 2 + t46;
t80 = t94 * qJD(6);
t103 = 2 * qJD(5);
t102 = pkin(4) + pkin(9);
t101 = t32 * t51;
t100 = t32 * t54;
t99 = t49 * t53;
t98 = t49 * t56;
t97 = t51 * t28;
t96 = t54 * t28;
t95 = qJ(4) + pkin(8);
t93 = cos(pkin(6));
t91 = qJD(2) * t56;
t90 = qJD(6) * t51;
t89 = qJD(6) * t54;
t86 = -0.2e1 * pkin(2) * qJD(3);
t85 = t51 * t96;
t44 = pkin(3) * t88;
t84 = t56 * t88;
t38 = t49 * t92;
t83 = t49 * t91;
t82 = t51 * t89;
t43 = -t55 * pkin(3) - pkin(2);
t42 = -pkin(3) * t50 - pkin(4);
t81 = qJD(3) * t95;
t27 = t55 * qJD(4) - t52 * t81;
t64 = -t52 * qJD(4) - t55 * t81;
t12 = t27 * t48 - t50 * t64;
t35 = t95 * t52;
t36 = t95 * t55;
t23 = t50 * t35 + t36 * t48;
t75 = -t33 * qJ(5) + t43;
t11 = t102 * t32 + t75;
t14 = pkin(5) * t33 + t23;
t79 = t11 * t54 + t14 * t51;
t78 = t11 * t51 - t14 * t54;
t13 = t50 * t27 + t48 * t64;
t24 = -t35 * t48 + t36 * t50;
t77 = t12 * t23 + t13 * t24;
t40 = pkin(3) * t48 + qJ(5);
t76 = -qJD(5) * t32 - t28 * t40;
t30 = t52 * t93 + t55 * t99;
t65 = -t52 * t99 + t55 * t93;
t16 = t48 * t30 - t50 * t65;
t74 = -t16 * t51 + t54 * t98;
t73 = t16 * t54 + t51 * t98;
t70 = t32 * t89 + t97;
t69 = t32 * t90 - t96;
t68 = t29 * t51 + t33 * t89;
t67 = -t29 * t54 + t33 * t90;
t66 = -t29 * qJ(5) - t33 * qJD(5) + t44;
t39 = -pkin(9) + t42;
t9 = -t28 * pkin(5) + t13;
t63 = t9 + (t32 * t40 - t33 * t39) * qJD(6);
t17 = t50 * t30 + t48 * t65;
t22 = qJD(3) * t65 + t55 * t83;
t57 = -qJD(3) * t30 - t52 * t83;
t6 = t48 * t22 - t50 * t57;
t7 = t50 * t22 + t48 * t57;
t62 = t16 * t12 + t17 * t13 + t6 * t23 + t7 * t24;
t61 = t16 * t29 - t17 * t28 - t32 * t7 + t33 * t6;
t15 = -pkin(5) * t32 + t24;
t60 = -qJD(6) * t15 - t29 * t39 - t76;
t59 = -0.2e1 * t49 ^ 2 * t53 * t91 + 0.2e1 * t16 * t6 + 0.2e1 * t17 * t7;
t58 = 0.2e1 * t12 * t33 - 0.2e1 * t13 * t32 + 0.2e1 * t23 * t29 - 0.2e1 * t24 * t28;
t31 = t32 ^ 2;
t20 = pkin(4) * t32 + t75;
t10 = pkin(4) * t28 + t66;
t8 = pkin(5) * t29 + t12;
t5 = t102 * t28 + t66;
t4 = qJD(6) * t73 + t38 * t54 + t51 * t6;
t3 = qJD(6) * t74 - t38 * t51 + t54 * t6;
t2 = -qJD(6) * t79 - t51 * t5 + t54 * t8;
t1 = qJD(6) * t78 - t54 * t5 - t51 * t8;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t38, -t83, 0, 0, 0, 0, 0, (-t55 * t92 - t84) * t49, (t52 * t92 - t56 * t87) * t49, t104, t105, t61, (-pkin(3) * t84 + t43 * t92) * t49 + t62, t61, -t104, -t105, (-t10 * t56 + t20 * t92) * t49 + t62, 0, 0, 0, 0, 0, -t100 * t7 + t17 * t69 + t29 * t73 + t3 * t33, t101 * t7 + t17 * t70 + t29 * t74 - t4 * t33; 0, 0, 0, 0, 0.2e1 * t52 * t87, 0.2e1 * (-t52 ^ 2 + t55 ^ 2) * qJD(3), 0, 0, 0, t52 * t86, t55 * t86, 0.2e1 * t28 * t43 + 0.2e1 * t32 * t44, 0.2e1 * t29 * t43 + 0.2e1 * t33 * t44, t58, 0.2e1 * t43 * t44 + 0.2e1 * t77, t58, -0.2e1 * t10 * t32 - 0.2e1 * t20 * t28, -0.2e1 * t10 * t33 - 0.2e1 * t20 * t29, 0.2e1 * t10 * t20 + 0.2e1 * t77, 0.2e1 * t28 * t32 * t46 + 0.2e1 * t31 * t82, -0.2e1 * t31 * t80 + 0.4e1 * t32 * t85, 0.2e1 * t32 * t68 + 0.2e1 * t33 * t97, -0.2e1 * t32 * t67 + 0.2e1 * t33 * t96, 0.2e1 * t33 * t29, -0.2e1 * t100 * t9 + 0.2e1 * t15 * t69 + 0.2e1 * t2 * t33 - 0.2e1 * t29 * t78, 0.2e1 * t1 * t33 + 0.2e1 * t101 * t9 + 0.2e1 * t15 * t70 - 0.2e1 * t29 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t22, -t6, -t7, 0, (t48 * t7 - t50 * t6) * pkin(3), 0, t6, t7, qJD(5) * t17 + t40 * t7 + t42 * t6, 0, 0, 0, 0, 0, t17 * t89 + t51 * t7, -t17 * t90 + t54 * t7; 0, 0, 0, 0, 0, 0, t87, -t88, 0, -pkin(8) * t87, pkin(8) * t88, -t12, -t13, (-t28 * t48 - t29 * t50) * pkin(3), (-t12 * t50 + t13 * t48) * pkin(3), t29 * t42 + t76, t12, t13, qJD(5) * t24 + t12 * t42 + t13 * t40, -t32 * t80 + t85, -t28 * t94 - 0.4e1 * t32 * t82, -t67, -t68, 0, t51 * t63 - t54 * t60, t51 * t60 + t54 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t40 * t103, -0.2e1 * t82, 0.2e1 * t80, 0, 0, 0, 0.2e1 * qJD(5) * t51 + 0.2e1 * t40 * t89, 0.2e1 * qJD(5) * t54 - 0.2e1 * t40 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t29, 0, t44, 0, -t28, -t29, t10, 0, 0, 0, 0, 0, -t68, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, t12, 0, 0, 0, 0, 0, -t67, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t69, t29, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t89, 0, -t39 * t90, -t39 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t18;
