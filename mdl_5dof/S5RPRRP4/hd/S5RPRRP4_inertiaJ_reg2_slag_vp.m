% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:24
% EndTime: 2020-01-03 11:50:27
% DurationCPUTime: 0.66s
% Computational Cost: add. (502->87), mult. (1083->146), div. (0->0), fcn. (1143->6), ass. (0->62)
t42 = sin(pkin(8));
t66 = -0.2e1 * t42;
t43 = cos(pkin(8));
t65 = -0.2e1 * t43;
t64 = 0.2e1 * t43;
t63 = t43 * pkin(4);
t44 = sin(qJ(4));
t62 = t44 * pkin(3);
t46 = cos(qJ(4));
t37 = t46 * pkin(3);
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t26 = -t44 * t45 + t46 * t47;
t61 = t26 * t43;
t60 = t45 * t42;
t59 = t47 * t42;
t58 = t47 * t43;
t25 = pkin(3) * t60 + t42 * qJ(2);
t40 = t45 ^ 2;
t41 = t47 ^ 2;
t57 = t40 + t41;
t56 = qJ(2) * t45;
t38 = t42 ^ 2;
t55 = t38 * qJ(2);
t54 = t42 * t64;
t53 = qJ(2) * t58;
t28 = -pkin(2) * t43 - pkin(6) * t42 - pkin(1);
t10 = t53 + (-pkin(7) * t42 + t28) * t45;
t24 = t47 * t28;
t8 = -pkin(7) * t59 + t24 + (-pkin(3) - t56) * t43;
t3 = -t44 * t10 + t46 * t8;
t4 = t10 * t46 + t44 * t8;
t15 = -t43 * t56 + t24;
t16 = t28 * t45 + t53;
t52 = t15 * t47 + t16 * t45;
t27 = t44 * t47 + t45 * t46;
t21 = -t44 * t60 + t46 * t59;
t51 = -t21 * qJ(5) + t3;
t19 = t27 * t42;
t2 = -qJ(5) * t19 + t4;
t50 = pkin(3) ^ 2;
t49 = 0.2e1 * pkin(4);
t48 = qJ(2) ^ 2;
t39 = t43 ^ 2;
t36 = t44 ^ 2 * t50;
t34 = t38 * t48;
t33 = -0.2e1 * t62;
t32 = t37 + pkin(4);
t30 = t43 * t62;
t23 = t27 * t43;
t22 = t27 * t62;
t18 = t21 ^ 2;
t17 = t19 ^ 2;
t14 = t19 * t62;
t13 = t21 * t65;
t12 = t19 * t65;
t11 = t26 ^ 2 + t27 ^ 2;
t9 = pkin(4) * t19 + t25;
t6 = -0.2e1 * t21 * t19;
t5 = -t19 * t27 - t21 * t26;
t1 = t51 - t63;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t38, t54, 0, t39, 0, 0, pkin(1) * t64, pkin(1) * t66, 0.2e1 * (t38 + t39) * qJ(2), pkin(1) ^ 2 + t39 * t48 + t34, t41 * t38, -0.2e1 * t47 * t38 * t45, t58 * t66, t40 * t38, t45 * t54, t39, -0.2e1 * t15 * t43 + 0.2e1 * t45 * t55, 0.2e1 * t16 * t43 + 0.2e1 * t47 * t55, t52 * t66, t15 ^ 2 + t16 ^ 2 + t34, t18, t6, t13, t17, -t12, t39, 0.2e1 * t19 * t25 - 0.2e1 * t3 * t43, 0.2e1 * t21 * t25 + 0.2e1 * t4 * t43, -0.2e1 * t19 * t4 - 0.2e1 * t21 * t3, t25 ^ 2 + t3 ^ 2 + t4 ^ 2, t18, t6, t13, t17, -t12, t39, -0.2e1 * t1 * t43 + 0.2e1 * t19 * t9, 0.2e1 * t2 * t43 + 0.2e1 * t21 * t9, -0.2e1 * t1 * t21 - 0.2e1 * t19 * t2, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t42, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t58, t45 * t43, -t57 * t42, t52, 0, 0, 0, 0, 0, 0, -t61, t23, t5, t26 * t3 + t27 * t4, 0, 0, 0, 0, 0, 0, -t61, t23, t5, t1 * t26 + t2 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, -t60, -t43, t15, -t16, 0, 0, 0, 0, t21, 0, -t19, -t43, -t37 * t43 + t3, t30 - t4, -t21 * t37 - t14, (t3 * t46 + t4 * t44) * pkin(3), 0, 0, t21, 0, -t19, -t43, (-pkin(4) - t32) * t43 + t51, -t2 + t30, -t21 * t32 - t14, t1 * t32 + t2 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t45, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t27, 0, t26 * t37 + t22, 0, 0, 0, 0, 0, 0, t26, -t27, 0, t26 * t32 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, t33, 0, t46 ^ 2 * t50 + t36, 0, 0, 0, 0, 0, 1, 0.2e1 * t32, t33, 0, t32 ^ 2 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, -t43, t3, -t4, 0, 0, 0, 0, t21, 0, -t19, -t43, t51 - 0.2e1 * t63, -t2, -t21 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t27, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t27, 0, t26 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t62, 0, 0, 0, 0, 0, 0, 0, 1, t49 + t37, -t62, 0, t32 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t49, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t21, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
