% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:50
% EndTime: 2019-12-31 17:47:51
% DurationCPUTime: 0.26s
% Computational Cost: add. (151->47), mult. (409->111), div. (0->0), fcn. (352->6), ass. (0->48)
t57 = pkin(3) + qJ(2);
t28 = sin(pkin(7));
t24 = t28 ^ 2;
t30 = cos(pkin(7));
t26 = t30 ^ 2;
t56 = (t24 + t26) * qJD(2);
t55 = -0.2e1 * t30;
t36 = -t28 * qJ(3) - pkin(1);
t13 = (-pkin(2) - qJ(4)) * t30 + t36;
t17 = t57 * t28;
t27 = sin(pkin(8));
t29 = cos(pkin(8));
t54 = t29 * t13 + t27 * t17;
t53 = t27 * t30;
t52 = t29 * t30;
t51 = qJ(2) * t56;
t50 = t57 * t30;
t31 = sin(qJ(5));
t49 = qJD(5) * t31;
t32 = cos(qJ(5));
t48 = qJD(5) * t32;
t47 = t26 * qJD(2);
t46 = t28 * qJD(2);
t45 = t28 * qJD(3);
t44 = t30 * qJD(2);
t42 = 0.2e1 * t52;
t41 = qJD(5) * t29 ^ 2 * t30;
t40 = t27 * t49;
t39 = t29 * t49;
t38 = t27 * t48;
t37 = t29 * t48;
t35 = t30 * t38;
t16 = -t30 * qJD(4) - t45;
t6 = t27 * t16 - t29 * t46;
t7 = t29 * t16 + t27 * t46;
t34 = t6 * t27 + t7 * t29;
t33 = -t27 * t13 + t29 * t17;
t11 = t28 * t32 + t31 * t53;
t15 = 0.2e1 * t56;
t12 = t28 * t31 - t32 * t53;
t10 = -t28 * t49 + t35;
t9 = t11 * qJD(5);
t5 = (pkin(4) * t29 + pkin(6) * t27) * t30 + t50;
t4 = t28 * pkin(6) + t54;
t3 = -t28 * pkin(4) - t33;
t2 = t32 * t44 - t31 * t7 + (-t31 * t5 - t32 * t4) * qJD(5);
t1 = -t31 * t44 - t32 * t7 + (t31 * t4 - t32 * t5) * qJD(5);
t8 = [0, 0, 0, 0, 0, t15, 0.2e1 * t51, t15, t45 * t55, 0.2e1 * t24 * qJD(3), -0.2e1 * (-t30 * pkin(2) + t36) * t45 + 0.2e1 * t51, -0.2e1 * t6 * t28 + 0.2e1 * t29 * t47, -0.2e1 * t27 * t47 - 0.2e1 * t7 * t28, t34 * t55, -0.2e1 * t33 * t6 + 0.2e1 * t50 * t44 + 0.2e1 * t54 * t7, 0.2e1 * t12 * t9, 0.2e1 * t12 * t10 + 0.2e1 * t9 * t11, t9 * t42, t10 * t42, 0, -0.2e1 * t3 * t10 - 0.2e1 * t6 * t11 + 0.2e1 * t2 * t52, 0.2e1 * t1 * t52 + 0.2e1 * t6 * t12 + 0.2e1 * t3 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, 0, 0, t34, 0, 0, 0, 0, 0, -t27 * t10 - t32 * t41, t27 * t9 + t31 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, t7 * t27 - t6 * t29, 0, 0, 0, 0, 0, (t10 - t35) * t29, (t30 * t40 - t9) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, -t30 * t39, -t30 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
