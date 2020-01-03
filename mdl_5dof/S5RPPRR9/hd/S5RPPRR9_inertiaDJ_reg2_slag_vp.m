% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:51
% DurationCPUTime: 0.86s
% Computational Cost: add. (675->116), mult. (1334->233), div. (0->0), fcn. (1087->6), ass. (0->88)
t33 = sin(qJ(5));
t27 = t33 ^ 2;
t35 = cos(qJ(5));
t29 = t35 ^ 2;
t57 = qJD(5) * (t27 - t29);
t34 = sin(qJ(4));
t28 = t34 ^ 2;
t36 = cos(qJ(4));
t30 = t36 ^ 2;
t87 = t28 - t30;
t56 = t87 * qJD(4);
t95 = 2 * qJD(2);
t94 = -pkin(1) - pkin(2);
t93 = pkin(4) * t34;
t92 = pkin(7) * t36;
t91 = t33 * t36;
t90 = t35 * t36;
t31 = sin(pkin(8));
t32 = cos(pkin(8));
t89 = t32 * qJ(2) + t31 * t94;
t86 = qJD(5) * t33;
t85 = qJD(5) * t34;
t84 = qJD(5) * t35;
t83 = qJD(5) * t36;
t82 = t31 * qJD(2);
t81 = t32 * qJD(2);
t80 = t34 * qJD(4);
t79 = t36 * qJD(2);
t78 = t36 * qJD(4);
t77 = -0.2e1 * pkin(4) * qJD(5);
t76 = t31 * t91;
t75 = t31 * t90;
t18 = -pkin(6) + t89;
t74 = t18 * t91;
t73 = t18 * t90;
t72 = 0.2e1 * t81;
t71 = pkin(7) * t84;
t70 = t31 * t81;
t69 = t35 * t81;
t68 = t29 * t78;
t67 = t35 * t78;
t66 = qJD(5) * t28 * t31;
t65 = t33 * t83;
t64 = t35 * t83;
t63 = t33 * t84;
t62 = t34 * t78;
t61 = t31 * t80;
t60 = t35 * t80;
t59 = t31 * t78;
t58 = t32 * t79;
t55 = 0.2e1 * t62;
t54 = t33 * t67;
t53 = t28 * t63;
t52 = t36 * pkin(4) + t34 * pkin(7);
t51 = t92 - t93;
t45 = -t31 * qJ(2) + t32 * t94;
t17 = pkin(3) - t45;
t39 = t17 + t52;
t3 = t35 * t39 - t74;
t4 = t33 * t39 + t73;
t50 = t3 * t35 + t33 * t4;
t49 = -t3 * t33 + t35 * t4;
t12 = -t33 * t32 + t75;
t46 = t35 * t32 + t76;
t48 = t12 * t33 - t35 * t46;
t47 = t12 * t35 + t33 * t46;
t44 = t18 * t80 - t58;
t43 = -t18 * t78 - t34 * t81;
t13 = t33 * t85 - t67;
t42 = t60 + t65;
t15 = t33 * t78 + t34 * t84;
t41 = t51 * qJD(4) + t82;
t37 = qJD(5) * t39;
t1 = t42 * t18 - t33 * t41 + (-t37 - t58) * t35;
t2 = (-t18 * t83 + t41) * t35 + (-t37 + t44) * t33;
t40 = -t50 * qJD(5) - t1 * t35 - t2 * t33;
t5 = t46 * qJD(5) + t31 * t60;
t6 = -qJD(5) * t12 + t33 * t61;
t38 = -t48 * qJD(5) - t33 * t6 - t35 * t5;
t24 = t27 * t78;
t23 = -0.2e1 * t62;
t22 = t29 * t62;
t21 = t27 * t62;
t19 = t28 * t70;
t16 = -t33 * t80 + t64;
t9 = t28 * t18 * t81;
t8 = t34 * t57 - t54;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, qJ(2) * t95, 0, 0, 0, 0, 0, 0, 0.2e1 * t82, t72, 0, (-t45 * t31 + t89 * t32) * t95, t55, -0.2e1 * t56, 0, t23, 0, 0, -0.2e1 * t17 * t80 + 0.2e1 * t31 * t79, -0.2e1 * t17 * t78 - 0.2e1 * t34 * t82, (-t28 - t30) * t72, 0.2e1 * t9 + 0.2e1 * (t18 * t30 * t32 + t17 * t31) * qJD(2), 0.2e1 * t22 - 0.2e1 * t53, 0.2e1 * t28 * t57 - 0.4e1 * t34 * t54, 0.2e1 * t34 * t65 + 0.2e1 * t35 * t56, 0.2e1 * t21 + 0.2e1 * t53, -0.2e1 * t33 * t56 + 0.2e1 * t34 * t64, t23, 0.2e1 * t2 * t36 + 0.2e1 * (-t18 * t84 - t33 * t81) * t28 + 0.2e1 * (-t3 - 0.2e1 * t74) * t80, 0.2e1 * t1 * t36 + 0.2e1 * (t18 * t86 - t69) * t28 + 0.2e1 * (t4 - 0.2e1 * t73) * t80, 0.2e1 * t50 * t78 + 0.2e1 * (t49 * qJD(5) - t1 * t33 + t2 * t35) * t34, 0.2e1 * t18 ^ 2 * t62 - 0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t80, t32 * t78, 0, t19 + (t30 - 0.1e1) * t70, 0, 0, 0, 0, 0, 0, -t35 * t66 + t36 * t6 + (t46 - 0.2e1 * t76) * t80, t33 * t66 + t36 * t5 + (t12 - 0.2e1 * t75) * t80, t48 * t78 + (t47 * qJD(5) - t33 * t5 + t35 * t6) * t34, t18 * t31 * t55 - t1 * t12 - t2 * t46 + t3 * t6 - t4 * t5 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t31 ^ 2 * t62 - 0.2e1 * t12 * t5 - 0.2e1 * t46 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t87 * t18 + t49 * t36) * qJD(4) + (t40 - t58) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t34 + (t87 * t31 + t47 * t36) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t21 + 0.2e1 * t22 - 0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, t80, 0, t43, t44, 0, 0, t8, 0.4e1 * t34 * t63 + t24 - t68, t16, -t8, -t42, 0, (-t71 + (pkin(4) * t33 - t18 * t35) * qJD(4)) * t36 + (pkin(7) * qJD(4) * t33 - t69 + (pkin(4) * t35 + t18 * t33) * qJD(5)) * t34, (t52 * qJD(4) + t18 * t85) * t35 + (t51 * qJD(5) - t43) * t33, t40, t43 * pkin(4) + t40 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t61, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t31, t15 * t31, t38, -pkin(4) * t59 + t38 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t78, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t16, t24 + t68, (-t93 + (t27 + t29) * t92) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t63, -0.2e1 * t57, 0, -0.2e1 * t63, 0, 0, t33 * t77, t35 * t77, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t15, -t80, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, -t86, 0, -t71, pkin(7) * t86, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
