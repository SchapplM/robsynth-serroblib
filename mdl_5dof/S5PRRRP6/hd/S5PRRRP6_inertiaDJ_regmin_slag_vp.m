% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:17
% EndTime: 2019-12-05 16:52:20
% DurationCPUTime: 0.41s
% Computational Cost: add. (486->100), mult. (1317->161), div. (0->0), fcn. (1134->6), ass. (0->52)
t60 = cos(qJ(4));
t45 = t60 * qJD(4);
t63 = t60 * qJD(3) + t45;
t62 = qJD(3) + qJD(4);
t38 = 2 * qJD(5);
t61 = -pkin(7) - pkin(6);
t33 = sin(qJ(4));
t34 = sin(qJ(3));
t59 = t33 * t34;
t36 = cos(qJ(3));
t58 = t33 * t36;
t57 = qJD(4) * t33;
t56 = t34 * qJD(3);
t35 = sin(qJ(2));
t55 = t35 * qJD(2);
t54 = t36 * qJD(3);
t37 = cos(qJ(2));
t53 = t37 * qJD(2);
t52 = -0.2e1 * pkin(2) * qJD(3);
t51 = pkin(3) * t56;
t50 = pkin(3) * t57;
t49 = t35 * t56;
t48 = t34 * t53;
t31 = -t36 * pkin(3) - pkin(2);
t47 = t61 * qJD(3);
t44 = t61 * t60;
t43 = t60 * t53;
t42 = t34 * t44;
t41 = qJD(3) * t44;
t19 = t60 * t34 + t58;
t8 = -t63 * t36 + t62 * t59;
t40 = t19 * t55 + t37 * t8;
t39 = t60 * t36 - t59;
t9 = t62 * t19;
t32 = pkin(3) * t45;
t30 = -t60 * pkin(3) - pkin(4);
t28 = t33 * pkin(3) + qJ(5);
t27 = -0.2e1 * t50;
t24 = t32 + qJD(5);
t23 = t61 * t36;
t13 = t39 * t35;
t12 = t19 * t35;
t11 = -t60 * t23 + t61 * t59;
t10 = -t33 * t23 - t42;
t7 = -pkin(4) * t39 - t19 * qJ(5) + t31;
t6 = -t37 * t9 - t39 * t55;
t5 = -t23 * t45 - t36 * t41 + (qJD(4) * t61 + t47) * t59;
t4 = -qJD(4) * t42 - t23 * t57 - t34 * t41 - t47 * t58;
t3 = -t33 * t49 + (t33 * t53 + t63 * t35) * t36 + (-t35 * t57 + t43) * t34;
t2 = t33 * t48 + t9 * t35 - t36 * t43;
t1 = t9 * pkin(4) + t8 * qJ(5) - t19 * qJD(5) + t51;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12 * t3 - 0.2e1 * t13 * t2 - 0.2e1 * t35 * t53; 0, 0, -t55, -t53, 0, 0, 0, 0, 0, -t36 * t55 - t37 * t56, t34 * t55 - t37 * t54, 0, 0, 0, 0, 0, t6, t40, t6, -t12 * t8 - t13 * t9 + t3 * t19 - t2 * t39, -t40, -t37 * t1 + t3 * t10 - t2 * t11 + t12 * t5 - t13 * t4 + t7 * t55; 0, 0, 0, 0, 0.2e1 * t34 * t54, 0.2e1 * (-t34 ^ 2 + t36 ^ 2) * qJD(3), 0, 0, 0, t34 * t52, t36 * t52, -0.2e1 * t19 * t8, -0.2e1 * t19 * t9 - 0.2e1 * t39 * t8, 0, 0, 0, 0.2e1 * t31 * t9 - 0.2e1 * t39 * t51, 0.2e1 * t19 * t51 - 0.2e1 * t31 * t8, -0.2e1 * t1 * t39 + 0.2e1 * t7 * t9, -0.2e1 * t10 * t8 - 0.2e1 * t11 * t9 + 0.2e1 * t5 * t19 - 0.2e1 * t39 * t4, -0.2e1 * t1 * t19 + 0.2e1 * t7 * t8, 0.2e1 * t7 * t1 + 0.2e1 * t10 * t5 - 0.2e1 * t11 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * t54 - t48, -t36 * t53 + t49, 0, 0, 0, 0, 0, -t3, t2, -t3, 0, -t2, t12 * t50 + t13 * t24 - t2 * t28 + t3 * t30; 0, 0, 0, 0, 0, 0, t54, -t56, 0, -pkin(6) * t54, pkin(6) * t56, 0, 0, -t8, -t9, 0, -t5, t4, -t5, t19 * t50 + t24 * t39 - t28 * t9 - t30 * t8, -t4, t10 * t50 + t11 * t24 - t4 * t28 + t5 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -0.2e1 * t32, t27, 0, 0.2e1 * t24, 0.2e1 * t28 * t24 + 0.2e1 * t30 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, -t3, 0, -t2, -t3 * pkin(4) - t2 * qJ(5) + t13 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, -t5, t4, -t5, pkin(4) * t8 - t9 * qJ(5) + qJD(5) * t39, -t4, -t5 * pkin(4) - t4 * qJ(5) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t32, -t50, 0, t38 + t32, -pkin(4) * t50 + t24 * qJ(5) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, qJ(5) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
