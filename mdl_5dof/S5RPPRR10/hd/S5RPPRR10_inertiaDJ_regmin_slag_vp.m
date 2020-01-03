% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:17
% EndTime: 2019-12-31 18:04:18
% DurationCPUTime: 0.30s
% Computational Cost: add. (274->54), mult. (675->104), div. (0->0), fcn. (640->6), ass. (0->46)
t34 = sin(pkin(8));
t32 = t34 ^ 2;
t35 = cos(pkin(8));
t60 = (t35 ^ 2 + t32) * qJD(2);
t39 = cos(qJ(4));
t37 = sin(qJ(4));
t56 = t35 * t37;
t23 = t34 * t39 - t56;
t55 = -pkin(6) + qJ(2);
t24 = t55 * t34;
t25 = t55 * t35;
t45 = t39 * t24 - t37 * t25;
t7 = -t23 * pkin(7) + t45;
t22 = t34 * t37 + t35 * t39;
t59 = t22 * qJD(2);
t41 = t23 * qJD(4);
t58 = qJD(4) + qJD(5);
t54 = qJ(2) * t60;
t36 = sin(qJ(5));
t53 = qJD(5) * t36;
t38 = cos(qJ(5));
t52 = qJD(5) * t38;
t51 = t34 * qJD(2);
t50 = t34 * qJD(3);
t48 = -t35 * pkin(2) - t34 * qJ(3) - pkin(1);
t47 = pkin(4) * t53;
t46 = pkin(4) * t52;
t20 = t35 * pkin(3) - t48;
t10 = -t36 * t22 + t38 * t23;
t44 = -t37 * t24 - t39 * t25;
t40 = -qJD(2) * t56 + t44 * qJD(4) + t39 * t51;
t21 = 0.2e1 * t60;
t19 = t22 * qJD(4);
t14 = pkin(4) * t41 + t50;
t13 = t58 * (-t36 * t39 - t37 * t38);
t12 = t58 * (t36 * t37 - t38 * t39);
t11 = t22 * pkin(4) + t20;
t9 = t38 * t22 + t36 * t23;
t8 = -t22 * pkin(7) - t44;
t6 = t19 * pkin(7) + t40;
t5 = t7 * qJD(4) + t59;
t4 = t10 * qJD(5) - t36 * t19 + t38 * t41;
t3 = t38 * t19 + t22 * t52 + t23 * t53 + t36 * t41;
t2 = -t36 * t5 + t38 * t6 + (-t36 * t7 - t38 * t8) * qJD(5);
t1 = -t36 * t6 - t38 * t5 + (t36 * t8 - t38 * t7) * qJD(5);
t15 = [0, 0, 0, 0, 0, t21, 0.2e1 * t54, 0.2e1 * t35 * t50, t21, 0.2e1 * t32 * qJD(3), -0.2e1 * t48 * t50 + 0.2e1 * t54, -0.2e1 * t23 * t19, 0.2e1 * t19 * t22 - 0.2e1 * t23 * t41, 0, 0, 0, 0.2e1 * t20 * t41 + 0.2e1 * t22 * t50, -0.2e1 * t20 * t19 + 0.2e1 * t23 * t50, -0.2e1 * t10 * t3, -0.2e1 * t10 * t4 + 0.2e1 * t3 * t9, 0, 0, 0, 0.2e1 * t11 * t4 + 0.2e1 * t14 * t9, 0.2e1 * t14 * t10 - 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, 0, 0, 0, 0, 0, -t41, t19, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t41, 0, t40, -t45 * qJD(4) - t59, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * qJD(4), -t39 * qJD(4), 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t47, -0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
