% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:41
% EndTime: 2019-12-05 15:38:44
% DurationCPUTime: 0.49s
% Computational Cost: add. (364->64), mult. (1039->118), div. (0->0), fcn. (954->6), ass. (0->50)
t37 = sin(pkin(8));
t38 = cos(pkin(8));
t68 = cos(qJ(4));
t54 = qJD(4) * t68;
t39 = sin(qJ(4));
t63 = qJD(4) * t39;
t72 = -t37 * t63 + t38 * t54;
t46 = -t39 * t37 + t68 * t38;
t71 = 0.2e1 * t72;
t66 = t39 * t38;
t26 = t68 * t37 + t66;
t22 = t26 * qJD(4);
t70 = 0.2e1 * t22;
t69 = 2 * qJD(5);
t65 = pkin(6) + qJ(3);
t64 = t37 ^ 2 + t38 ^ 2;
t40 = sin(qJ(2));
t34 = t40 * qJD(2);
t41 = cos(qJ(2));
t62 = t41 * qJD(2);
t61 = t46 * t70;
t59 = t40 * t62;
t33 = -t38 * pkin(3) - pkin(2);
t29 = t65 * t38;
t56 = t65 * t37;
t48 = t68 * t56;
t17 = t39 * t29 + t48;
t18 = t68 * t29 - t39 * t56;
t53 = t68 * qJD(3);
t8 = qJD(4) * t48 - t38 * t53 + (qJD(3) * t37 + qJD(4) * t29) * t39;
t9 = t29 * t54 + qJD(3) * t66 + (-t65 * t63 + t53) * t37;
t58 = t17 * t9 - t18 * t8;
t55 = t64 * t41;
t52 = t64 * qJD(3);
t51 = 0.2e1 * t52;
t49 = -t26 * t22 + t46 * t72;
t47 = t26 * t34 - t41 * t72;
t11 = t40 * t22 - t46 * t62;
t12 = t26 * t62 + t72 * t40;
t19 = t26 * t40;
t20 = t46 * t40;
t45 = -t11 * t18 + t12 * t17 + t19 * t9 - t20 * t8;
t44 = -t11 * t46 + t12 * t26 + t19 * t72 - t20 * t22;
t43 = 0.2e1 * t17 * t72 - 0.2e1 * t18 * t22 + 0.2e1 * t9 * t26 - 0.2e1 * t46 * t8;
t42 = -0.2e1 * t20 * t11 + 0.2e1 * t19 * t12 - 0.2e1 * t59;
t16 = t26 * t71;
t14 = -pkin(4) * t46 - t26 * qJ(5) + t33;
t13 = -t41 * t22 - t34 * t46;
t7 = t22 * pkin(4) - qJ(5) * t72 - t26 * qJD(5);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t64) * t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t62, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t34, t37 * t34, qJD(2) * t55, t40 * t52 + (-pkin(2) * t40 + qJ(3) * t55) * qJD(2), 0, 0, 0, 0, 0, 0, t13, t47, t44, t33 * t34 + t45, 0, 0, 0, 0, 0, 0, t13, t44, -t47, t14 * t34 - t41 * t7 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, qJ(3) * t51, t16, 0.2e1 * t49, 0, -t61, 0, 0, t33 * t70, t33 * t71, t43, 0.2e1 * t58, t16, 0, -0.2e1 * t49, 0, 0, -t61, 0.2e1 * t14 * t22 - 0.2e1 * t46 * t7, t43, -0.2e1 * t14 * t72 - 0.2e1 * t7 * t26, 0.2e1 * t14 * t7 + 0.2e1 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t72, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t72, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t11, -t12 * pkin(4) - t11 * qJ(5) + t20 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, -t22, 0, -t9, t8, 0, 0, 0, t72, 0, 0, t22, 0, -t9, -pkin(4) * t72 - t22 * qJ(5) + qJD(5) * t46, -t8, -t9 * pkin(4) - t8 * qJ(5) + t18 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, qJ(5) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
