% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:51
% EndTime: 2019-12-31 17:47:52
% DurationCPUTime: 0.50s
% Computational Cost: add. (446->68), mult. (994->153), div. (0->0), fcn. (907->6), ass. (0->55)
t31 = sin(pkin(7));
t27 = t31 ^ 2;
t33 = cos(pkin(7));
t29 = t33 ^ 2;
t66 = (t27 + t29) * qJD(2);
t65 = -0.2e1 * t33;
t55 = t29 * qJD(2);
t30 = sin(pkin(8));
t53 = t31 * qJD(3);
t19 = -t33 * qJD(4) - t53;
t32 = cos(pkin(8));
t54 = t31 * qJD(2);
t9 = t30 * t19 - t32 * t54;
t64 = t9 * t30;
t63 = t9 * t32;
t62 = t30 * t33;
t61 = t32 * t33;
t60 = pkin(3) + qJ(2);
t44 = -t31 * qJ(3) - pkin(1);
t16 = (-pkin(2) - qJ(4)) * t33 + t44;
t20 = t60 * t31;
t59 = t32 * t16 + t30 * t20;
t58 = qJ(2) * t66;
t34 = sin(qJ(5));
t57 = qJD(5) * t34;
t35 = cos(qJ(5));
t56 = qJD(5) * t35;
t52 = t33 * qJD(2);
t50 = 0.2e1 * t61;
t49 = qJD(5) * t32 ^ 2 * t33;
t48 = t30 * t57;
t47 = t32 * t57;
t46 = t30 * t56;
t45 = t32 * t56;
t43 = t33 * t46;
t10 = t32 * t19 + t30 * t54;
t42 = t10 * t32 + t64;
t41 = -t30 * t16 + t32 * t20;
t14 = t31 * t35 + t34 * t62;
t40 = pkin(4) * t32 + pkin(6) * t30 + t60;
t39 = t40 * t33;
t38 = t35 * t39;
t6 = t31 * pkin(6) + t59;
t1 = -qJD(5) * t38 - t35 * t10 - t34 * t52 + t6 * t57;
t2 = -t6 * t56 - t34 * t10 + (t35 * qJD(2) - t40 * t57) * t33;
t3 = -t34 * t6 + t38;
t4 = t34 * t39 + t35 * t6;
t37 = -t1 * t35 - t2 * t34 + (-t3 * t35 - t34 * t4) * qJD(5);
t12 = qJD(5) * t14;
t13 = -t31 * t57 + t43;
t15 = -t31 * t34 + t35 * t62;
t36 = t12 * t34 + t13 * t35 + (-t14 * t34 - t15 * t35) * qJD(5);
t18 = 0.2e1 * t66;
t5 = -t31 * pkin(4) - t41;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0.2e1 * t58, 0, 0, 0, 0, 0, 0, t18, t53 * t65, 0.2e1 * t27 * qJD(3), -0.2e1 * (-t33 * pkin(2) + t44) * t53 + 0.2e1 * t58, 0, 0, 0, 0, 0, 0, -0.2e1 * t9 * t31 + 0.2e1 * t32 * t55, -0.2e1 * t10 * t31 - 0.2e1 * t30 * t55, t42 * t65, 0.2e1 * t59 * t10 - 0.2e1 * t41 * t9 + 0.2e1 * t60 * t55, -0.2e1 * t15 * t12, 0.2e1 * t12 * t14 - 0.2e1 * t15 * t13, t12 * t50, 0.2e1 * t14 * t13, t13 * t50, 0, -0.2e1 * t5 * t13 - 0.2e1 * t9 * t14 + 0.2e1 * t2 * t61, 0.2e1 * t1 * t61 + 0.2e1 * t5 * t12 - 0.2e1 * t9 * t15, -0.2e1 * t1 * t14 - 0.2e1 * t3 * t12 + 0.2e1 * t4 * t13 + 0.2e1 * t2 * t15, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, -t30 * t13 - t35 * t49, t30 * t12 + t34 * t49, t36 * t32, t37 * t32 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t30 - t63, 0, 0, 0, 0, 0, 0, (t13 - t43) * t32, (t33 * t48 - t12) * t32, t36 * t30, t37 * t30 - t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, -t33 * t47, -t33 * t45, -t35 * t12 + t34 * t13 + (t14 * t35 - t15 * t34) * qJD(5), -t1 * t34 + t2 * t35 + (-t3 * t34 + t35 * t4) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t13, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
