% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:55
% EndTime: 2020-01-03 11:22:58
% DurationCPUTime: 0.62s
% Computational Cost: add. (798->100), mult. (1882->215), div. (0->0), fcn. (1905->8), ass. (0->63)
t39 = sin(pkin(9));
t42 = cos(pkin(9));
t44 = cos(pkin(7));
t41 = sin(pkin(7));
t43 = cos(pkin(8));
t67 = t41 * t43;
t24 = t39 * t67 + t44 * t42;
t74 = -0.2e1 * t24;
t40 = sin(pkin(8));
t61 = qJD(3) * t41;
t62 = qJD(2) * t44;
t27 = -t40 * t61 + t43 * t62;
t18 = -t44 * qJD(4) + t27;
t51 = (-qJD(4) * t43 + qJD(2)) * t41;
t9 = t39 * t18 - t42 * t51;
t73 = t9 * t42;
t26 = t40 * t62 + t43 * t61;
t72 = t26 * t43;
t71 = t39 * t40;
t70 = t40 * t41;
t45 = sin(qJ(5));
t69 = t40 * t45;
t46 = cos(qJ(5));
t68 = t40 * t46;
t31 = -t44 * pkin(2) - t41 * qJ(3) - pkin(1);
t66 = t43 * t31;
t63 = qJ(2) * t44;
t64 = t40 * t31 + t43 * t63;
t16 = -t44 * qJ(4) + t64;
t19 = (pkin(3) * t40 - qJ(4) * t43 + qJ(2)) * t41;
t65 = t42 * t16 + t39 * t19;
t60 = qJD(5) * t45;
t59 = qJD(5) * t46;
t37 = t41 ^ 2;
t58 = t37 * qJD(2);
t57 = qJ(2) * qJD(2);
t56 = t41 * t68;
t55 = t39 * t60;
t54 = t39 * t59;
t53 = -t39 * t16 + t42 * t19;
t52 = -t27 * t40 + t72;
t25 = -t44 * t39 + t42 * t67;
t15 = t25 * t46 + t41 * t69;
t29 = t42 * t68 - t45 * t43;
t50 = t42 * t69 + t46 * t43;
t49 = -t66 + (qJ(2) * t40 + pkin(3)) * t44;
t48 = t24 * pkin(4) - t25 * pkin(6) + t49;
t47 = t46 * t48;
t7 = pkin(6) * t70 + t65;
t4 = t45 * t48 + t46 * t7;
t38 = t44 ^ 2;
t36 = t37 * t57;
t22 = qJD(5) * t29;
t21 = t50 * qJD(5);
t14 = t25 * t45 - t56;
t12 = qJD(5) * t15;
t11 = -qJD(5) * t56 + t25 * t60;
t10 = t42 * t18 + t39 * t51;
t6 = -pkin(4) * t70 - t53;
t3 = -t45 * t7 + t47;
t2 = -qJD(5) * t4 - t45 * t10 + t46 * t26;
t1 = -qJD(5) * t47 - t46 * t10 - t45 * t26 + t60 * t7;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t37 + t38) * qJD(2), 0.2e1 * t38 * t57 + 0.2e1 * t36, 0, 0, 0, 0, 0, 0, 0.2e1 * t26 * t44 + 0.2e1 * t40 * t58, 0.2e1 * t27 * t44 + 0.2e1 * t43 * t58, 0.2e1 * t52 * t41, 0.2e1 * t64 * t27 - 0.2e1 * (-t40 * t63 + t66) * t26 + 0.2e1 * t36, 0, 0, 0, 0, 0, 0, 0.2e1 * t26 * t24 - 0.2e1 * t70 * t9, -0.2e1 * t10 * t70 + 0.2e1 * t26 * t25, -0.2e1 * t10 * t24 + 0.2e1 * t9 * t25, 0.2e1 * t10 * t65 + 0.2e1 * t26 * t49 - 0.2e1 * t53 * t9, -0.2e1 * t15 * t11, 0.2e1 * t11 * t14 - 0.2e1 * t15 * t12, t11 * t74, 0.2e1 * t14 * t12, t12 * t74, 0, 0.2e1 * t6 * t12 + 0.2e1 * t9 * t14 + 0.2e1 * t2 * t24, 0.2e1 * t1 * t24 - 0.2e1 * t6 * t11 + 0.2e1 * t9 * t15, 0.2e1 * t1 * t14 + 0.2e1 * t3 * t11 - 0.2e1 * t4 * t12 - 0.2e1 * t2 * t15, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 + (t10 * t42 + t39 * t9) * t40, 0, 0, 0, 0, 0, 0, t12 * t71 - t22 * t24, -t11 * t71 + t21 * t24, -t11 * t50 - t29 * t12 + t21 * t14 + t22 * t15, -t1 * t29 - t2 * t50 - t4 * t21 - t3 * t22 + t71 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t29 * t21 + 0.2e1 * t22 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t39 - t73, 0, 0, 0, 0, 0, 0, -t42 * t12 - t24 * t54, t42 * t11 + t24 * t55, (-t11 * t45 - t12 * t46 + (t14 * t45 + t15 * t46) * qJD(5)) * t39, -t73 + (-t1 * t46 - t2 * t45 + (-t3 * t46 - t4 * t45) * qJD(5)) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t21 * t46 + t22 * t45 + (-t29 * t45 + t46 * t50) * qJD(5)) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, -t24 * t60, -t24 * t59, t46 * t11 - t45 * t12 + (-t14 * t46 + t15 * t45) * qJD(5), -t1 * t45 + t2 * t46 + (-t3 * t45 + t4 * t46) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t45 - t22 * t46 + (t29 * t46 + t45 * t50) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t55, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
