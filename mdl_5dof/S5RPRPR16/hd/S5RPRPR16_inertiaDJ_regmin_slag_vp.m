% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR16_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:33
% DurationCPUTime: 0.36s
% Computational Cost: add. (188->73), mult. (441->138), div. (0->0), fcn. (290->4), ass. (0->59)
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t59 = t31 * qJ(4);
t66 = -t29 * pkin(3) + t59;
t25 = t29 ^ 2;
t27 = t31 ^ 2;
t42 = (t25 - t27) * qJD(3);
t28 = sin(qJ(5));
t24 = t28 ^ 2;
t30 = cos(qJ(5));
t61 = -t30 ^ 2 + t24;
t41 = t61 * qJD(5);
t32 = -pkin(3) - pkin(7);
t52 = t29 * qJD(4);
t33 = -pkin(1) - pkin(6);
t62 = pkin(4) - t33;
t14 = t62 * t29;
t53 = t14 * qJD(5);
t65 = (t29 * t32 + t59) * qJD(3) + t52 + t53;
t64 = 2 * qJD(2);
t63 = 0.2e1 * qJD(4);
t58 = qJD(3) * t14;
t57 = qJD(5) * t28;
t56 = qJD(5) * t30;
t55 = qJD(5) * t31;
t54 = qJD(5) * t32;
t22 = t29 * qJD(3);
t51 = t31 * qJD(3);
t50 = qJ(2) * qJD(3);
t49 = qJ(4) * qJD(5);
t48 = t28 * t55;
t47 = t30 * t55;
t46 = t30 * t51;
t45 = t28 * t56;
t44 = t29 * t51;
t43 = pkin(3) * t51 + qJ(4) * t22 + qJD(2);
t40 = qJD(5) * (t25 + t27);
t39 = t28 * t46;
t13 = qJ(2) - t66;
t10 = pkin(7) * t29 + t13;
t15 = t62 * t31;
t38 = t10 * t30 + t15 * t28;
t37 = t10 * t28 - t15 * t30;
t19 = t33 * t51;
t12 = -pkin(4) * t51 + t19;
t34 = t12 + (qJ(4) * t29 - t31 * t32) * qJD(5);
t5 = qJD(3) * t66 + t52;
t18 = t33 * t22;
t16 = -0.2e1 * t44;
t11 = -pkin(4) * t22 + t18;
t9 = -t22 * t28 + t47;
t8 = t28 * t51 + t29 * t56;
t7 = t22 * t30 + t48;
t6 = -t29 * t57 + t46;
t4 = -qJD(4) * t31 + t43;
t3 = (qJD(3) * pkin(7) - qJD(4)) * t31 + t43;
t2 = -qJD(5) * t38 + t30 * t11 - t28 * t3;
t1 = qJD(5) * t37 - t28 * t11 - t30 * t3;
t17 = [0, 0, 0, 0, t64, qJ(2) * t64, t16, 0.2e1 * t42, 0, 0, 0, 0.2e1 * qJD(2) * t29 + 0.2e1 * t31 * t50, 0.2e1 * qJD(2) * t31 - 0.2e1 * t29 * t50, 0, -0.2e1 * t13 * t51 - 0.2e1 * t29 * t4, 0.2e1 * t13 * t22 - 0.2e1 * t31 * t4, 0.2e1 * t13 * t4, 0.2e1 * t24 * t44 + 0.2e1 * t25 * t45, -0.2e1 * t25 * t41 + 0.4e1 * t29 * t39, -0.2e1 * t28 * t42 + 0.2e1 * t29 * t47, -0.2e1 * t29 * t48 - 0.2e1 * t30 * t42, t16, 0.2e1 * (t30 * t58 + t2) * t31 + 0.2e1 * (qJD(3) * t37 - t12 * t30 - t28 * t53) * t29, 0.2e1 * (-t28 * t58 + t1) * t31 + 0.2e1 * (qJD(3) * t38 + t12 * t28 - t30 * t53) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t40, t30 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t51, 0, -t18, -t19, -t5, t18, t19, t5 * t33, -t29 * t41 + t39, -0.4e1 * t29 * t45 - t51 * t61, -t7, -t9, 0, t34 * t28 - t30 * t65, t28 * t65 + t34 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t51, 0, t22, t51, t5, 0, 0, 0, 0, 0, t8, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, qJ(4) * t63, -0.2e1 * t45, 0.2e1 * t41, 0, 0, 0, 0.2e1 * qJD(4) * t28 + 0.2e1 * t30 * t49, 0.2e1 * qJD(4) * t30 - 0.2e1 * t28 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, t18, 0, 0, 0, 0, 0, -t7, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t6, -t22, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56, 0, -t28 * t54, -t30 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
