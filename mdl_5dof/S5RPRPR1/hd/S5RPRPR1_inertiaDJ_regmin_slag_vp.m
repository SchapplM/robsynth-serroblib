% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:43
% EndTime: 2019-12-05 17:47:44
% DurationCPUTime: 0.33s
% Computational Cost: add. (422->61), mult. (877->117), div. (0->0), fcn. (797->6), ass. (0->47)
t47 = sin(pkin(8));
t51 = cos(qJ(3));
t62 = cos(pkin(8));
t58 = qJD(3) * t62;
t49 = sin(qJ(3));
t61 = t49 * qJD(3);
t28 = t47 * t61 - t51 * t58;
t60 = t51 * qJD(3);
t29 = -t47 * t60 - t49 * t58;
t67 = (-t28 * t47 + t62 * t29) * pkin(3);
t33 = -t47 * t49 + t62 * t51;
t34 = -t47 * t51 - t62 * t49;
t48 = sin(qJ(5));
t50 = cos(qJ(5));
t14 = t50 * t33 + t48 * t34;
t5 = t14 * qJD(5) - t50 * t28 + t48 * t29;
t66 = 2 * qJD(2);
t65 = pkin(3) * t47;
t52 = -pkin(1) - pkin(6);
t63 = qJ(4) - t52;
t26 = -t51 * qJD(4) + t63 * t61;
t36 = t63 * t51;
t27 = -qJD(3) * t36 - t49 * qJD(4);
t10 = t47 * t26 + t62 * t27;
t35 = t63 * t49;
t16 = -t62 * t35 - t47 * t36;
t64 = t49 * pkin(3) + qJ(2);
t37 = pkin(3) * t60 + qJD(2);
t59 = qJ(2) * qJD(3);
t9 = t62 * t26 - t47 * t27;
t15 = t47 * t35 - t62 * t36;
t57 = t28 * t34 + t33 * t29;
t13 = t48 * t33 - t50 * t34;
t54 = t10 * t34 - t15 * t29 + t16 * t28 - t9 * t33;
t53 = -t13 * qJD(5) + t48 * t28 + t50 * t29;
t44 = t62 * pkin(3) + pkin(4);
t23 = (-t44 * t48 - t50 * t65) * qJD(5);
t22 = (-t44 * t50 + t48 * t65) * qJD(5);
t21 = -t34 * pkin(4) + t64;
t17 = -t28 * pkin(4) + t37;
t12 = t34 * pkin(7) + t16;
t11 = -t33 * pkin(7) + t15;
t8 = t28 * pkin(7) + t10;
t7 = -t29 * pkin(7) + t9;
t2 = -t48 * t8 + t50 * t7 + (-t11 * t48 - t12 * t50) * qJD(5);
t1 = -t48 * t7 - t50 * t8 + (-t11 * t50 + t12 * t48) * qJD(5);
t3 = [0, 0, 0, 0, t66, qJ(2) * t66, -0.2e1 * t49 * t60, 0.2e1 * (t49 ^ 2 - t51 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t49 + 0.2e1 * t51 * t59, 0.2e1 * qJD(2) * t51 - 0.2e1 * t49 * t59, 0.2e1 * t54, 0.2e1 * t16 * t10 + 0.2e1 * t15 * t9 + 0.2e1 * t64 * t37, 0.2e1 * t14 * t53, -0.2e1 * t13 * t53 - 0.2e1 * t14 * t5, 0, 0, 0, 0.2e1 * t17 * t13 + 0.2e1 * t21 * t5, 0.2e1 * t17 * t14 + 0.2e1 * t21 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t57, -t54, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t60, 0, -t52 * t61, -t52 * t60, -t67, (t10 * t47 + t62 * t9) * pkin(3), 0, 0, t53, -t5, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t60, 0, t67, 0, 0, 0, 0, 0, t53, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23, 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, t5, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t5, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
