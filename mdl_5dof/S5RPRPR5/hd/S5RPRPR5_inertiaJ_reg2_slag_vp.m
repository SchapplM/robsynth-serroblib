% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:00
% EndTime: 2022-01-23 09:26:02
% DurationCPUTime: 0.65s
% Computational Cost: add. (740->87), mult. (1569->176), div. (0->0), fcn. (1709->8), ass. (0->56)
t42 = sin(pkin(8));
t66 = -0.2e1 * t42;
t44 = cos(pkin(8));
t65 = -0.2e1 * t44;
t64 = 0.2e1 * t44;
t41 = sin(pkin(9));
t63 = t41 * pkin(3);
t43 = cos(pkin(9));
t62 = t43 * pkin(3);
t61 = t44 * pkin(3);
t46 = sin(qJ(3));
t60 = t46 * t42;
t48 = cos(qJ(3));
t59 = t48 * t42;
t58 = t48 * t44;
t31 = -pkin(2) * t44 - pkin(6) * t42 - pkin(1);
t27 = t48 * t31;
t55 = qJ(4) * t42;
t56 = qJ(2) * t46;
t13 = -t48 * t55 + t27 + (-pkin(3) - t56) * t44;
t52 = qJ(2) * t58;
t18 = t52 + (t31 - t55) * t46;
t7 = t41 * t13 + t43 * t18;
t30 = pkin(3) * t60 + t42 * qJ(2);
t39 = t46 ^ 2;
t40 = t48 ^ 2;
t57 = t39 + t40;
t37 = t42 ^ 2;
t54 = t37 * qJ(2);
t53 = t42 * t64;
t24 = -t41 * t60 + t43 * t59;
t6 = t43 * t13 - t18 * t41;
t4 = -pkin(4) * t44 - pkin(7) * t24 + t6;
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t29 = t41 * t48 + t43 * t46;
t22 = t29 * t42;
t5 = -pkin(7) * t22 + t7;
t1 = t47 * t4 - t45 * t5;
t2 = t4 * t45 + t47 * t5;
t20 = -t44 * t56 + t27;
t21 = t31 * t46 + t52;
t51 = t20 * t48 + t21 * t46;
t49 = qJ(2) ^ 2;
t38 = t44 ^ 2;
t35 = t37 * t49;
t34 = pkin(4) + t62;
t28 = -t41 * t46 + t43 * t48;
t26 = t34 * t45 + t47 * t63;
t25 = t34 * t47 - t45 * t63;
t16 = pkin(4) * t22 + t30;
t15 = t28 * t45 + t29 * t47;
t14 = t28 * t47 - t29 * t45;
t10 = -t22 * t45 + t24 * t47;
t8 = t47 * t22 + t24 * t45;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t53, 0, t38, 0, 0, pkin(1) * t64, pkin(1) * t66, 0.2e1 * (t37 + t38) * qJ(2), pkin(1) ^ 2 + t38 * t49 + t35, t40 * t37, -0.2e1 * t48 * t37 * t46, t58 * t66, t39 * t37, t46 * t53, t38, -0.2e1 * t20 * t44 + 0.2e1 * t46 * t54, 0.2e1 * t21 * t44 + 0.2e1 * t48 * t54, t51 * t66, t20 ^ 2 + t21 ^ 2 + t35, t24 ^ 2, -0.2e1 * t24 * t22, t24 * t65, t22 ^ 2, -t22 * t65, t38, 0.2e1 * t22 * t30 - 0.2e1 * t44 * t6, 0.2e1 * t24 * t30 + 0.2e1 * t44 * t7, -0.2e1 * t22 * t7 - 0.2e1 * t24 * t6, t30 ^ 2 + t6 ^ 2 + t7 ^ 2, t10 ^ 2, -0.2e1 * t10 * t8, t10 * t65, t8 ^ 2, t8 * t64, t38, -0.2e1 * t1 * t44 + 0.2e1 * t16 * t8, 0.2e1 * t10 * t16 + 0.2e1 * t2 * t44, -0.2e1 * t1 * t10 - 0.2e1 * t2 * t8, t1 ^ 2 + t16 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t42, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t58, t46 * t44, -t57 * t42, t51, 0, 0, 0, 0, 0, 0, -t28 * t44, t29 * t44, -t22 * t29 - t24 * t28, t28 * t6 + t29 * t7, 0, 0, 0, 0, 0, 0, -t14 * t44, t15 * t44, -t10 * t14 - t15 * t8, t1 * t14 + t15 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 ^ 2 + t29 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, -t60, -t44, t20, -t21, 0, 0, 0, 0, t24, 0, -t22, -t44, -t43 * t61 + t6, t41 * t61 - t7, (-t22 * t41 - t24 * t43) * pkin(3), (t41 * t7 + t43 * t6) * pkin(3), 0, 0, t10, 0, -t8, -t44, -t25 * t44 + t1, t26 * t44 - t2, -t10 * t25 - t26 * t8, t1 * t25 + t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t46, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, (t28 * t43 + t29 * t41) * pkin(3), 0, 0, 0, 0, 0, 0, t14, -t15, 0, t14 * t25 + t15 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t62, -0.2e1 * t63, 0, (t41 ^ 2 + t43 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t26, 0, t25 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t24, 0, t30, 0, 0, 0, 0, 0, 0, t8, t10, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t8, -t44, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
