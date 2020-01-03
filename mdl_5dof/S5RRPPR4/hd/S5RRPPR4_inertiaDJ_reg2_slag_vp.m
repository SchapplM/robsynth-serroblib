% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:49
% EndTime: 2019-12-31 19:27:50
% DurationCPUTime: 0.43s
% Computational Cost: add. (239->59), mult. (528->111), div. (0->0), fcn. (343->6), ass. (0->48)
t35 = sin(qJ(5));
t31 = t35 ^ 2;
t37 = cos(qJ(5));
t32 = t37 ^ 2;
t53 = t31 + t32;
t38 = cos(qJ(2));
t52 = pkin(1) * qJD(2);
t28 = t38 * t52;
t20 = t28 + qJD(3);
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t36 = sin(qJ(2));
t48 = t36 * t52;
t7 = t34 * t20 + t33 * t48;
t55 = t53 * t7;
t40 = 0.2e1 * qJD(3);
t17 = t34 * t48;
t6 = t33 * t20 - t17;
t54 = t6 * t34;
t46 = -t38 * pkin(1) - pkin(2);
t24 = -pkin(3) + t46;
t26 = t36 * pkin(1) + qJ(3);
t5 = t33 * t24 + t34 * t26;
t39 = -pkin(2) - pkin(3);
t12 = t34 * qJ(3) + t33 * t39;
t51 = t33 * qJD(3);
t50 = t34 * qJD(3);
t30 = t35 * qJD(5);
t49 = t37 * qJD(5);
t47 = t35 * t49;
t4 = t34 * t24 - t33 * t26;
t2 = pkin(4) - t4;
t11 = -t33 * qJ(3) + t34 * t39;
t9 = pkin(4) - t11;
t43 = qJD(5) * (-t2 - t9);
t42 = t53 * t34;
t41 = t53 * t50;
t25 = -0.2e1 * t48;
t23 = t34 * t49;
t22 = t34 * t30;
t21 = t37 * t51;
t19 = -0.2e1 * t47;
t18 = 0.2e1 * t47;
t10 = -pkin(7) + t12;
t8 = 0.2e1 * (-t31 + t32) * qJD(5);
t3 = -pkin(7) + t5;
t1 = t6 * t37;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -0.2e1 * t28, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0.2e1 * t20, 0.2e1 * t26 * t20 + 0.2e1 * t46 * t48, 0, 0, 0, 0, 0, 0, 0.2e1 * t6, 0.2e1 * t7, 0, -0.2e1 * t4 * t6 + 0.2e1 * t5 * t7, t18, t8, 0, t19, 0, 0, -0.2e1 * t2 * t30 + 0.2e1 * t1, -0.2e1 * t2 * t49 - 0.2e1 * t6 * t35, -0.2e1 * t55, 0.2e1 * t2 * t6 + 0.2e1 * t3 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t28, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, t40 + t28, -pkin(2) * t48 + t20 * qJ(3) + t26 * qJD(3), 0, 0, 0, 0, 0, 0, -t17 + (qJD(3) + t20) * t33, t50 + t7, 0, -t6 * t11 + t7 * t12 + (-t33 * t4 + t34 * t5) * qJD(3), t18, t8, 0, t19, 0, 0, t35 * t43 + t1 + t21, (-t6 - t51) * t35 + t37 * t43, -t55 - t41, t6 * t9 + t10 * t55 + (t2 * t33 + t3 * t42) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, qJ(3) * t40, 0, 0, 0, 0, 0, 0, 0.2e1 * t51, 0.2e1 * t50, 0, (-t11 * t33 + t12 * t34) * t40, t18, t8, 0, t19, 0, 0, -0.2e1 * t9 * t30 + 0.2e1 * t21, -0.2e1 * t35 * t51 - 0.2e1 * t9 * t49, -0.2e1 * t41, (t10 * t42 + t33 * t9) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t33 - t54, 0, 0, 0, 0, 0, 0, t22, t23, 0, t33 * t55 - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, (-0.1e1 + t53) * t33 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, t30, 0, -t3 * t49 - t35 * t7, t3 * t30 - t37 * t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, t30, 0, -t10 * t49 - t35 * t50, t10 * t30 - t37 * t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33 * t49, t33 * t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
