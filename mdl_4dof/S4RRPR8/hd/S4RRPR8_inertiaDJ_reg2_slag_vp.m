% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:23
% EndTime: 2019-12-31 17:08:23
% DurationCPUTime: 0.34s
% Computational Cost: add. (288->68), mult. (686->132), div. (0->0), fcn. (506->4), ass. (0->48)
t33 = sin(qJ(4));
t34 = sin(qJ(2));
t54 = t34 * t33;
t56 = pkin(5) - pkin(6);
t60 = t56 * t54;
t35 = cos(qJ(2));
t32 = t35 * qJD(2);
t55 = cos(qJ(4));
t43 = t34 * t55;
t59 = qJD(2) * t43 - t33 * t32;
t58 = 2 * qJD(3);
t57 = pkin(2) + pkin(3);
t53 = qJ(3) * t35;
t52 = t34 * qJ(3);
t51 = qJD(4) * t33;
t50 = t34 * qJD(2);
t49 = t34 * qJD(3);
t48 = -0.2e1 * pkin(1) * qJD(2);
t47 = pkin(5) * t50;
t46 = pkin(5) * t32;
t44 = t34 * t32;
t42 = qJD(4) * t55;
t41 = t56 * t55;
t40 = t55 * t57;
t38 = -t35 * pkin(2) - t52;
t37 = t34 * t41;
t16 = t35 * t55 + t54;
t19 = t55 * qJ(3) - t33 * t57;
t36 = t38 * qJD(2) + t35 * qJD(3);
t24 = -0.2e1 * t44;
t23 = 0.2e1 * t44;
t22 = t56 * t35;
t21 = (-t34 ^ 2 + t35 ^ 2) * qJD(2);
t20 = -pkin(1) + t38;
t18 = -t33 * qJ(3) - t40;
t17 = -t35 * t33 + t43;
t14 = t57 * t35 + pkin(1) + t52;
t10 = -t49 + (pkin(2) * t34 - t53) * qJD(2);
t9 = t33 * qJD(3) + t19 * qJD(4);
t8 = qJ(3) * t51 - t55 * qJD(3) + qJD(4) * t40;
t7 = t49 + (-t57 * t34 + t53) * qJD(2);
t6 = t55 * t22 + t60;
t5 = -t33 * t22 + t37;
t4 = t16 * qJD(2) - t34 * t51 - t35 * t42;
t3 = t34 * t42 - t35 * t51 - t59;
t2 = t22 * t42 - t41 * t32 + (-qJD(2) + qJD(4)) * t60;
t1 = -qJD(4) * t37 + t22 * t51 + t59 * t56;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0.2e1 * t21, 0, t24, 0, 0, t34 * t48, t35 * t48, 0, 0, t23, 0, -0.2e1 * t21, 0, 0, t24, -0.2e1 * t10 * t35 + 0.2e1 * t20 * t50, 0, -0.2e1 * t10 * t34 - 0.2e1 * t20 * t32, 0.2e1 * t20 * t10, 0.2e1 * t17 * t4, -0.2e1 * t4 * t16 - 0.2e1 * t17 * t3, 0, 0.2e1 * t16 * t3, 0, 0, 0.2e1 * t14 * t3 + 0.2e1 * t7 * t16, 0.2e1 * t14 * t4 + 0.2e1 * t7 * t17, 0.2e1 * t1 * t16 + 0.2e1 * t2 * t17 - 0.2e1 * t6 * t3 - 0.2e1 * t5 * t4, -0.2e1 * t6 * t1 + 0.2e1 * t14 * t7 - 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t50, 0, -t46, t47, 0, 0, 0, t32, 0, 0, t50, 0, -t46, t36, -t47, t36 * pkin(5), 0, 0, -t4, 0, t3, 0, t2, -t1, t8 * t16 + t9 * t17 - t18 * t4 - t19 * t3, -t1 * t19 - t2 * t18 - t5 * t9 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, qJ(3) * t58, 0, 0, 0, 0, 0, 0, 0.2e1 * t9, -0.2e1 * t8, 0, -0.2e1 * t18 * t9 - 0.2e1 * t19 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * t4 - t33 * t3 + (-t55 * t16 + t17 * t33) * qJD(4), -t2 * t55 - t1 * t33 + (-t33 * t5 + t55 * t6) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t42, 0, -t9 * t55 - t8 * t33 + (-t18 * t33 + t55 * t19) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, -t3, 0, -t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
