% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR14_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:16
% EndTime: 2019-12-31 18:35:18
% DurationCPUTime: 0.60s
% Computational Cost: add. (428->63), mult. (734->118), div. (0->0), fcn. (815->6), ass. (0->53)
t29 = sin(pkin(8));
t30 = cos(pkin(8));
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t13 = -t29 * t32 + t30 * t34;
t63 = -0.2e1 * t13;
t11 = t29 * t34 + t30 * t32;
t62 = (t11 * t29 + t13 * t30) * pkin(3);
t10 = t13 ^ 2;
t9 = t11 ^ 2;
t61 = t9 + t10;
t35 = -pkin(1) - pkin(6);
t48 = -qJ(4) + t35;
t16 = t48 * t32;
t43 = t48 * t34;
t4 = t29 * t16 - t30 * t43;
t60 = t4 ^ 2;
t22 = t32 * pkin(3) + qJ(2);
t59 = 0.2e1 * t22;
t58 = 0.2e1 * qJ(2);
t57 = t29 * pkin(3);
t56 = t30 * pkin(3);
t55 = t4 * t13;
t21 = -pkin(4) - t56;
t54 = t13 * t21;
t31 = sin(qJ(5));
t53 = t31 * t11;
t52 = t31 * t13;
t33 = cos(qJ(5));
t51 = t31 * t33;
t50 = t33 * t13;
t24 = t31 ^ 2;
t26 = t33 ^ 2;
t49 = t24 + t26;
t25 = t32 ^ 2;
t27 = t34 ^ 2;
t17 = t25 + t27;
t47 = t11 * t63;
t46 = t31 * t50;
t20 = pkin(7) + t57;
t44 = t49 * t20;
t3 = t11 * pkin(4) - t13 * pkin(7) + t22;
t6 = t30 * t16 + t29 * t43;
t1 = t33 * t3 - t31 * t6;
t2 = t31 * t3 + t33 * t6;
t42 = t1 * t33 + t2 * t31;
t41 = -t1 * t31 + t2 * t33;
t40 = t6 * t11 - t55;
t39 = -t11 * t20 + t54;
t36 = qJ(2) ^ 2;
t15 = t17 * t35;
t7 = t33 * t11;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t58, (pkin(1) ^ 2) + t36, t27, -0.2e1 * t34 * t32, 0, t25, 0, 0, t32 * t58, t34 * t58, -0.2e1 * t15, t17 * t35 ^ 2 + t36, t10, t47, 0, t9, 0, 0, t11 * t59, t13 * t59, -0.2e1 * t40, t22 ^ 2 + t6 ^ 2 + t60, t26 * t10, -0.2e1 * t10 * t51, 0.2e1 * t11 * t50, t24 * t10, t31 * t47, t9, 0.2e1 * t1 * t11 + 0.2e1 * t4 * t52, -0.2e1 * t2 * t11 + 0.2e1 * t4 * t50, t42 * t63, t1 ^ 2 + t2 ^ 2 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t17, t15, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t40, 0, 0, 0, 0, 0, 0, -t61 * t31, -t61 * t33, 0, t11 * t41 - t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t9 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t32, 0, t34 * t35, -t32 * t35, 0, 0, 0, 0, t13, 0, -t11, 0, -t4, -t6, -t62, (t29 * t6 - t30 * t4) * pkin(3), t46, (-t24 + t26) * t13, t53, -t46, t7, 0, t31 * t39 - t4 * t33, t4 * t31 + t33 * t39, t41, t20 * t41 + t4 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t11, 0, t62, 0, 0, 0, 0, 0, 0, t50, -t52, t49 * t11, t11 * t44 - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t57, 0, (t29 ^ 2 + t30 ^ 2) * pkin(3) ^ 2, t24, 0.2e1 * t51, 0, t26, 0, 0, -0.2e1 * t21 * t33, 0.2e1 * t21 * t31, 0.2e1 * t44, t49 * t20 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t13, 0, t22, 0, 0, 0, 0, 0, 0, t7, -t53, -t49 * t13, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t52, t11, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t33, 0, -t31 * t20, -t33 * t20, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
