% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:17
% EndTime: 2019-12-31 18:04:20
% DurationCPUTime: 0.60s
% Computational Cost: add. (783->95), mult. (1794->181), div. (0->0), fcn. (1716->6), ass. (0->59)
t50 = sin(pkin(8));
t48 = t50 ^ 2;
t51 = cos(pkin(8));
t81 = (t51 ^ 2 + t48) * qJD(2);
t53 = sin(qJ(4));
t54 = cos(qJ(4));
t34 = t50 * t54 - t51 * t53;
t58 = qJD(4) * t34;
t80 = qJD(4) + qJD(5);
t79 = t34 * pkin(7);
t78 = cos(qJ(5));
t52 = sin(qJ(5));
t77 = t52 * t53;
t75 = -pkin(6) + qJ(2);
t37 = t75 * t51;
t76 = t53 * t37;
t62 = t75 * t50;
t19 = t54 * t37 + t53 * t62;
t74 = qJ(2) * t81;
t73 = qJD(5) * t52;
t72 = t50 * qJD(3);
t71 = t53 * qJD(2);
t70 = t53 * qJD(4);
t69 = t54 * qJD(2);
t68 = t54 * qJD(4);
t66 = -t51 * pkin(2) - t50 * qJ(3) - pkin(1);
t30 = t54 * t62;
t65 = qJD(4) * t30 + t50 * t71 + t51 * t69;
t64 = pkin(4) * t73;
t63 = t78 * t54;
t18 = t30 - t76;
t61 = t78 * qJD(5);
t27 = t51 * pkin(3) - t66;
t60 = pkin(4) * t61;
t33 = t50 * t53 + t51 * t54;
t59 = -t18 + t79;
t14 = -t52 * t33 + t34 * t78;
t36 = t52 * t54 + t53 * t78;
t57 = t78 * t59;
t12 = -t33 * pkin(7) + t19;
t4 = t12 * t78 - t52 * t59;
t56 = (-t76 - t79) * qJD(4) + t65;
t11 = -t51 * t71 - t37 * t68 + (-t75 * t70 + t69) * t50;
t26 = qJD(4) * t33;
t55 = t26 * pkin(7) + t11;
t35 = t63 - t77;
t32 = 0.2e1 * t81;
t20 = pkin(4) * t58 + t72;
t17 = t80 * t36;
t16 = -qJD(4) * t63 - t54 * t61 + t80 * t77;
t15 = t33 * pkin(4) + t27;
t13 = t33 * t78 + t52 * t34;
t10 = t37 * t70 - t65;
t6 = t14 * qJD(5) - t52 * t26 + t58 * t78;
t5 = t26 * t78 + t33 * t61 + t34 * t73 + t52 * t58;
t3 = -t52 * t12 - t57;
t2 = -t4 * qJD(5) - t52 * t56 + t55 * t78;
t1 = qJD(5) * t57 + t12 * t73 - t52 * t55 - t56 * t78;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0.2e1 * t74, 0, 0, 0, 0, 0, 0, 0.2e1 * t51 * t72, t32, 0.2e1 * t48 * qJD(3), -0.2e1 * t66 * t72 + 0.2e1 * t74, -0.2e1 * t34 * t26, 0.2e1 * t26 * t33 - 0.2e1 * t34 * t58, 0, 0.2e1 * t33 * t58, 0, 0, 0.2e1 * t27 * t58 + 0.2e1 * t33 * t72, -0.2e1 * t27 * t26 + 0.2e1 * t34 * t72, 0.2e1 * t10 * t33 - 0.2e1 * t11 * t34 + 0.2e1 * t18 * t26 - 0.2e1 * t19 * t58, -0.2e1 * t19 * t10 + 0.2e1 * t18 * t11 + 0.2e1 * t27 * t72, -0.2e1 * t14 * t5, 0.2e1 * t5 * t13 - 0.2e1 * t14 * t6, 0, 0.2e1 * t13 * t6, 0, 0, 0.2e1 * t20 * t13 + 0.2e1 * t15 * t6, 0.2e1 * t20 * t14 - 0.2e1 * t15 * t5, 0.2e1 * t1 * t13 - 0.2e1 * t2 * t14 + 0.2e1 * t3 * t5 - 0.2e1 * t4 * t6, -0.2e1 * t4 * t1 + 0.2e1 * t15 * t20 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, 0, 0, 0, 0, 0, 0, -t58, t26, 0, -t72, 0, 0, 0, 0, 0, 0, -t6, t5, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, t54 * t26 - t33 * t68, -t10 * t53 + t11 * t54 + (-t18 * t53 + t19 * t54) * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, t16 * t13 + t17 * t14 + t35 * t5 - t36 * t6, -t1 * t36 - t4 * t16 - t3 * t17 + t2 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t36 * t16 - 0.2e1 * t35 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t58, 0, t11, t10, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, (t78 * t5 - t52 * t6 + (-t13 * t78 + t14 * t52) * qJD(5)) * pkin(4), (t78 * t2 - t1 * t52 + (-t3 * t52 + t4 * t78) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t68, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, (-t78 * t17 - t16 * t52 + (-t35 * t52 + t78 * t36) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t64, -0.2e1 * t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
