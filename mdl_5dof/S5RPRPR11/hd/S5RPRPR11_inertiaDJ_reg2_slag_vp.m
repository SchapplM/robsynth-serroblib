% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:55
% EndTime: 2019-12-31 18:27:57
% DurationCPUTime: 0.73s
% Computational Cost: add. (901->95), mult. (2024->170), div. (0->0), fcn. (1904->6), ass. (0->58)
t51 = sin(pkin(8));
t52 = cos(pkin(8));
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t40 = t81 * t51 + t79 * t52;
t78 = pkin(6) + qJ(2);
t24 = t40 * t78;
t39 = t79 * t51 - t81 * t52;
t33 = t39 * qJD(3);
t85 = -0.2e1 * t33;
t34 = t40 * qJD(3);
t84 = 0.2e1 * t34;
t83 = 2 * qJD(4);
t82 = -pkin(3) - pkin(4);
t80 = cos(qJ(5));
t64 = t78 * t81;
t65 = t78 * t79;
t25 = -t51 * t65 + t52 * t64;
t53 = sin(qJ(5));
t77 = qJD(5) * t53;
t76 = t39 * t84;
t47 = -t52 * pkin(2) - pkin(1);
t68 = t79 * qJD(2);
t69 = t81 * qJD(2);
t16 = t24 * qJD(3) + t51 * t68 - t52 * t69;
t17 = (qJD(3) * t64 + t68) * t52 + (-qJD(3) * t65 + t69) * t51;
t71 = -t25 * t16 + t24 * t17;
t70 = qJD(5) * t80;
t67 = 0.2e1 * (t51 ^ 2 + t52 ^ 2) * qJD(2);
t66 = t80 * t82;
t63 = t33 * t39 - t40 * t34;
t62 = -t33 * qJ(4) + t40 * qJD(4);
t60 = t40 * qJ(4) - t47;
t22 = t53 * t39 + t80 * t40;
t42 = t80 * qJ(4) + t53 * t82;
t59 = 0.2e1 * t16 * t39 + 0.2e1 * t17 * t40 - 0.2e1 * t24 * t33 - 0.2e1 * t25 * t34;
t58 = -t40 * pkin(7) + t24;
t57 = t53 * t58;
t56 = t80 * t58;
t55 = t34 * pkin(7) - t16;
t54 = t33 * pkin(7) + t17;
t41 = -t53 * qJ(4) + t66;
t30 = t53 * qJD(4) + t42 * qJD(5);
t29 = qJ(4) * t77 - t80 * qJD(4) - qJD(5) * t66;
t23 = t40 * t85;
t21 = -t80 * t39 + t53 * t40;
t20 = t39 * pkin(3) - t60;
t19 = t39 * pkin(7) + t25;
t14 = t82 * t39 + t60;
t13 = t34 * pkin(3) - t62;
t8 = t82 * t34 + t62;
t6 = qJD(5) * t22 - t53 * t33 - t80 * t34;
t5 = t80 * t33 - t53 * t34 - t39 * t70 + t40 * t77;
t4 = t80 * t19 + t57;
t3 = -t53 * t19 + t56;
t2 = qJD(5) * t57 + t19 * t70 + t53 * t55 - t80 * t54;
t1 = -qJD(5) * t56 + t19 * t77 - t53 * t54 - t80 * t55;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, qJ(2) * t67, t23, 0.2e1 * t63, 0, t76, 0, 0, t47 * t84, t47 * t85, t59, 0.2e1 * t71, t23, 0, -0.2e1 * t63, 0, 0, t76, 0.2e1 * t13 * t39 + 0.2e1 * t20 * t34, t59, -0.2e1 * t13 * t40 + 0.2e1 * t20 * t33, 0.2e1 * t20 * t13 + 0.2e1 * t71, -0.2e1 * t22 * t5, 0.2e1 * t5 * t21 - 0.2e1 * t22 * t6, 0, 0.2e1 * t21 * t6, 0, 0, 0.2e1 * t14 * t6 + 0.2e1 * t8 * t21, -0.2e1 * t14 * t5 + 0.2e1 * t8 * t22, 0.2e1 * t1 * t21 + 0.2e1 * t2 * t22 + 0.2e1 * t3 * t5 - 0.2e1 * t4 * t6, -0.2e1 * t4 * t1 + 0.2e1 * t14 * t8 - 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t33, t13, 0, 0, 0, 0, 0, 0, -t6, t5, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, -t34, 0, -t17, t16, 0, 0, 0, -t33, 0, 0, t34, 0, -t17, pkin(3) * t33 - t34 * qJ(4) - t39 * qJD(4), -t16, -t17 * pkin(3) - t16 * qJ(4) + t25 * qJD(4), 0, 0, t5, 0, t6, 0, t2, -t1, t29 * t21 + t30 * t22 + t41 * t5 - t42 * t6, -t1 * t42 - t2 * t41 - t4 * t29 - t3 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, qJ(4) * t83, 0, 0, 0, 0, 0, 0, 0.2e1 * t30, -0.2e1 * t29, 0, -0.2e1 * t42 * t29 - 0.2e1 * t41 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t5 - t53 * t6 + (-t80 * t21 + t22 * t53) * qJD(5), -t2 * t80 - t1 * t53 + (-t3 * t53 + t80 * t4) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t70, 0, -t30 * t80 - t29 * t53 + (-t41 * t53 + t80 * t42) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, -t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t70, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
