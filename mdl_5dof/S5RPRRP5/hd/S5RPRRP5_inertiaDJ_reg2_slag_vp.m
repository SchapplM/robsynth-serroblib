% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:59
% DurationCPUTime: 0.48s
% Computational Cost: add. (270->52), mult. (703->87), div. (0->0), fcn. (451->6), ass. (0->43)
t39 = sin(qJ(3));
t44 = cos(pkin(8)) * pkin(1) + pkin(2);
t56 = cos(qJ(3));
t58 = sin(pkin(8)) * pkin(1);
t62 = t39 * t58 - t56 * t44;
t15 = t62 * qJD(3);
t38 = sin(qJ(4));
t35 = t38 ^ 2;
t40 = cos(qJ(4));
t36 = t40 ^ 2;
t64 = (-t35 - t36) * t15;
t61 = 0.2e1 * (-t35 + t36) * qJD(4);
t60 = 2 * qJD(5);
t51 = t39 * t44 + t56 * t58;
t18 = pkin(7) + t51;
t59 = t64 * t18;
t33 = t38 * qJD(4);
t34 = t40 * qJD(4);
t20 = -pkin(4) * t33 + qJ(5) * t34 + t38 * qJD(5);
t16 = t51 * qJD(3);
t7 = t16 - t20;
t57 = t20 - t7;
t53 = t64 * pkin(7);
t17 = -pkin(3) + t62;
t52 = t16 * t38 + t17 * t34;
t50 = pkin(3) * t33;
t49 = pkin(3) * t34;
t48 = pkin(7) * t33;
t47 = pkin(7) * t34;
t46 = t38 * t34;
t45 = -t16 * t40 + t17 * t33;
t43 = -t40 * pkin(4) - t38 * qJ(5);
t19 = t43 * qJD(4) + t40 * qJD(5);
t29 = -0.2e1 * t46;
t28 = 0.2e1 * t46;
t25 = -pkin(3) + t43;
t21 = t25 * t33;
t9 = t17 + t43;
t8 = t9 * t33;
t4 = -t38 * t15 + t18 * t34;
t3 = t40 * t15 + t18 * t33;
t1 = 0.2e1 * t64;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t16, 0.2e1 * t15, 0, -0.2e1 * t51 * t15 + 0.2e1 * t16 * t62, t28, t61, 0, t29, 0, 0, 0.2e1 * t45, 0.2e1 * t52, t1, 0.2e1 * t17 * t16 + 0.2e1 * t59, t28, 0, -t61, 0, 0, t29, -0.2e1 * t7 * t40 + 0.2e1 * t8, t1, -0.2e1 * t9 * t34 - 0.2e1 * t7 * t38, 0.2e1 * t9 * t7 + 0.2e1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, 0, 0, t28, t61, 0, t29, 0, 0, t45 - t50, -t49 + t52, t64, -t16 * pkin(3) + t53, t28, 0, -t61, 0, 0, t29, t57 * t40 + t21 + t8, t64, t57 * t38 + (-t25 - t9) * t34, -t9 * t20 + t7 * t25 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t61, 0, t29, 0, 0, -0.2e1 * t50, -0.2e1 * t49, 0, 0, t28, 0, -t61, 0, 0, t29, 0.2e1 * t20 * t40 + 0.2e1 * t21, 0, 0.2e1 * t20 * t38 - 0.2e1 * t25 * t34, -0.2e1 * t25 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t33, 0, -t4, t3, 0, 0, 0, t34, 0, 0, t33, 0, -t4, t19, -t3, (pkin(4) * t38 - qJ(5) * t40) * t15 + t19 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t34, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, t34, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t33, 0, -t47, t48, 0, 0, 0, t34, 0, 0, t33, 0, -t47, t19, -t48, t19 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, qJ(5) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
