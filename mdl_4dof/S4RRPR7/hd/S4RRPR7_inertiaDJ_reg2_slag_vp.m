% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:44
% EndTime: 2019-12-31 17:06:46
% DurationCPUTime: 0.50s
% Computational Cost: add. (697->85), mult. (1643->189), div. (0->0), fcn. (1428->6), ass. (0->65)
t31 = sin(qJ(4));
t27 = t31 ^ 2;
t33 = cos(qJ(4));
t28 = t33 ^ 2;
t64 = t27 - t28;
t50 = qJD(4) * t64;
t34 = cos(qJ(2));
t65 = -qJ(3) - pkin(5);
t21 = t65 * t34;
t29 = sin(pkin(7));
t30 = cos(pkin(7));
t32 = sin(qJ(2));
t53 = t65 * t32;
t12 = -t30 * t21 + t29 * t53;
t51 = qJD(2) * t65;
t14 = t34 * qJD(3) + t32 * t51;
t39 = t32 * qJD(3) - t34 * t51;
t35 = t30 * t14 - t29 * t39;
t18 = t29 * t32 - t30 * t34;
t19 = t29 * t34 + t30 * t32;
t26 = -t34 * pkin(2) - pkin(1);
t40 = -t18 * pkin(3) + t19 * pkin(6) - t26;
t37 = t33 * t40;
t15 = t19 * qJD(2);
t60 = t34 * qJD(2);
t61 = t32 * qJD(2);
t16 = -t29 * t61 + t30 * t60;
t56 = pkin(2) * t61;
t38 = t15 * pkin(3) - t16 * pkin(6) + t56;
t63 = qJD(4) * t31;
t1 = qJD(4) * t37 + t12 * t63 - t31 * t38 - t33 * t35;
t4 = t33 * t12 - t31 * t40;
t2 = -qJD(4) * t4 - t31 * t35 + t33 * t38;
t3 = -t31 * t12 - t37;
t73 = t1 * t31 - t2 * t33 + (t3 * t31 - t33 * t4) * qJD(4);
t11 = -t29 * t21 - t30 * t53;
t8 = t29 * t14 + t30 * t39;
t72 = t11 * t8;
t71 = t8 * t31;
t70 = t19 * t16;
t69 = t19 * t33;
t68 = t31 * t15;
t67 = t33 * t15;
t66 = t33 * t16;
t62 = qJD(4) * t33;
t59 = 0.2e1 * t18 * t15;
t58 = -0.2e1 * pkin(1) * qJD(2);
t25 = -t30 * pkin(2) - pkin(3);
t57 = 0.2e1 * qJD(4) * t25;
t55 = t31 * t62;
t54 = t32 * t60;
t52 = -0.4e1 * t31 * t69;
t17 = t19 ^ 2;
t49 = t17 * t55;
t47 = -t3 * t33 - t31 * t4;
t24 = t29 * pkin(2) + pkin(6);
t45 = -t15 * t24 + t16 * t25;
t44 = t18 * t24 - t19 * t25;
t43 = t18 * t62 + t68;
t42 = t31 * t16 + t19 * t62;
t41 = -t19 * t63 + t66;
t36 = t47 * qJD(4) - t1 * t33 - t2 * t31;
t10 = -t18 * t63 + t67;
t5 = t19 * t50 - t31 * t66;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t54, 0.2e1 * (-t32 ^ 2 + t34 ^ 2) * qJD(2), 0, -0.2e1 * t54, 0, 0, t32 * t58, t34 * t58, 0, 0, 0.2e1 * t70, -0.2e1 * t19 * t15 - 0.2e1 * t16 * t18, 0, t59, 0, 0, 0.2e1 * t26 * t15 + 0.2e1 * t18 * t56, 0.2e1 * t26 * t16 + 0.2e1 * t19 * t56, 0.2e1 * t11 * t16 - 0.2e1 * t12 * t15 - 0.2e1 * t35 * t18 + 0.2e1 * t8 * t19, 0.2e1 * t12 * t35 + 0.2e1 * t26 * t56 + 0.2e1 * t72, 0.2e1 * t28 * t70 - 0.2e1 * t49, t16 * t52 + 0.2e1 * t17 * t50, 0.2e1 * t41 * t18 + 0.2e1 * t19 * t67, 0.2e1 * t27 * t70 + 0.2e1 * t49, -0.2e1 * t42 * t18 - 0.2e1 * t19 * t68, t59, 0.2e1 * t42 * t11 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t18 + 0.2e1 * t19 * t71, 0.2e1 * t1 * t18 + 0.2e1 * t41 * t11 - 0.2e1 * t4 * t15 + 0.2e1 * t8 * t69, 0.2e1 * t47 * t16 + 0.2e1 * t73 * t19, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, -t61, 0, -pkin(5) * t60, pkin(5) * t61, 0, 0, 0, 0, t16, 0, -t15, 0, -t8, -t35, (-t15 * t29 - t16 * t30) * pkin(2), (t35 * t29 - t8 * t30) * pkin(2), -t5, qJD(4) * t52 - t64 * t16, t43, t5, t10, 0, -t8 * t33 + t45 * t31 + (t11 * t31 - t44 * t33) * qJD(4), t71 + t45 * t33 + (t11 * t33 + t44 * t31) * qJD(4), t36, t36 * t24 + t8 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t55, -0.2e1 * t50, 0, -0.2e1 * t55, 0, 0, t31 * t57, t33 * t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t16, 0, t56, 0, 0, 0, 0, 0, 0, t10, -t43, (-t27 - t28) * t16, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t42, t15, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t63, 0, -t24 * t62, t24 * t63, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t62, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
