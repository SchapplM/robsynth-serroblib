% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:37
% EndTime: 2019-12-31 18:24:38
% DurationCPUTime: 0.33s
% Computational Cost: add. (185->65), mult. (461->128), div. (0->0), fcn. (319->6), ass. (0->60)
t30 = cos(qJ(3));
t31 = -pkin(3) - pkin(7);
t66 = t30 * t31;
t25 = t30 ^ 2;
t28 = sin(qJ(3));
t42 = qJD(3) * (t28 ^ 2 - t25);
t27 = sin(qJ(5));
t22 = t27 ^ 2;
t29 = cos(qJ(5));
t61 = -t29 ^ 2 + t22;
t41 = t61 * qJD(5);
t51 = t30 * qJD(4);
t17 = sin(pkin(8)) * pkin(1) + pkin(6);
t62 = pkin(4) + t17;
t14 = t62 * t30;
t52 = t14 * qJD(5);
t58 = t28 * qJ(4);
t65 = (t58 - t66) * qJD(3) - t51 - t52;
t64 = 0.2e1 * qJD(4);
t63 = t30 * pkin(3);
t59 = qJ(4) * t30;
t57 = qJD(3) * t14;
t56 = qJD(5) * t27;
t55 = qJD(5) * t29;
t54 = qJD(5) * t30;
t53 = qJD(5) * t31;
t20 = t28 * qJD(3);
t21 = t30 * qJD(3);
t50 = qJ(4) * qJD(5);
t18 = -cos(pkin(8)) * pkin(1) - pkin(2);
t49 = 0.2e1 * qJD(3) * t18;
t48 = t27 * t54;
t47 = t29 * t54;
t46 = t28 * t21;
t45 = t17 * t20;
t44 = t29 * t20;
t43 = t27 * t55;
t40 = pkin(3) * t20 - t28 * qJD(4);
t39 = t27 * t44;
t13 = t62 * t28;
t34 = t18 - t58;
t4 = t34 + t66;
t38 = t29 * t13 - t27 * t4;
t37 = t27 * t13 + t29 * t4;
t5 = t62 * t20;
t33 = -t5 + (-t28 * t31 - t59) * qJD(5);
t32 = t51 + (-t58 - t63) * qJD(3);
t16 = 0.2e1 * t46;
t15 = t17 * t21;
t12 = t34 - t63;
t11 = -t27 * t20 + t47;
t10 = t27 * t21 + t28 * t55;
t9 = t44 + t48;
t8 = t29 * t21 - t28 * t56;
t7 = qJ(4) * t21 - t40;
t6 = pkin(4) * t21 + t15;
t3 = (pkin(7) * t28 - t59) * qJD(3) + t40;
t2 = -t37 * qJD(5) - t27 * t3 + t29 * t6;
t1 = -t38 * qJD(5) - t27 * t6 - t29 * t3;
t19 = [0, 0, 0, 0, t16, -0.2e1 * t42, 0, 0, 0, t28 * t49, t30 * t49, 0, -0.2e1 * t12 * t20 - 0.2e1 * t7 * t30, -0.2e1 * t12 * t21 + 0.2e1 * t7 * t28, -0.2e1 * t12 * t7, -0.2e1 * t22 * t46 + 0.2e1 * t25 * t43, -0.2e1 * t25 * t41 - 0.4e1 * t30 * t39, 0.2e1 * t27 * t42 - 0.2e1 * t28 * t47, 0.2e1 * t28 * t48 + 0.2e1 * t29 * t42, t16, 0.2e1 * (-t29 * t57 + t2) * t28 + 0.2e1 * (t38 * qJD(3) - t27 * t52 - t5 * t29) * t30, 0.2e1 * (t27 * t57 + t1) * t28 + 0.2e1 * (-t37 * qJD(3) + t5 * t27 - t29 * t52) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t15, t45, t32, t15, -t45, t32 * t17, t30 * t41 + t39, -t61 * t20 + 0.4e1 * t30 * t43, t8, -t10, 0, t33 * t27 - t65 * t29, t65 * t27 + t33 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, t20, t21, t7, 0, 0, 0, 0, 0, t10, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, qJ(4) * t64, -0.2e1 * t43, 0.2e1 * t41, 0, 0, 0, 0.2e1 * qJD(4) * t27 + 0.2e1 * t29 * t50, 0.2e1 * qJD(4) * t29 - 0.2e1 * t27 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, t15, 0, 0, 0, 0, 0, t8, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t9, t21, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t55, 0, -t27 * t53, -t29 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t19;
