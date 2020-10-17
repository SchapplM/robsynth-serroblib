% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR1
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:43
% EndTime: 2019-12-05 17:47:45
% DurationCPUTime: 0.48s
% Computational Cost: add. (425->54), mult. (711->90), div. (0->0), fcn. (832->6), ass. (0->46)
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t18 = t33 * t38 + t34 * t36;
t20 = -t33 * t36 + t34 * t38;
t35 = sin(qJ(5));
t37 = cos(qJ(5));
t42 = -t35 * t18 + t37 * t20;
t58 = t42 ^ 2;
t43 = t37 * t18 + t35 * t20;
t57 = t43 ^ 2;
t56 = (t18 * t33 + t20 * t34) * pkin(3);
t17 = t20 ^ 2;
t51 = t18 ^ 2;
t55 = t17 + t51;
t54 = t57 + t58;
t39 = -pkin(1) - pkin(6);
t22 = (-qJ(4) + t39) * t36;
t28 = t38 * t39;
t23 = -t38 * qJ(4) + t28;
t10 = -t33 * t22 + t34 * t23;
t3 = -t20 * pkin(7) + t10;
t11 = t34 * t22 + t33 * t23;
t4 = -t18 * pkin(7) + t11;
t1 = t37 * t3 - t35 * t4;
t2 = t35 * t3 + t37 * t4;
t53 = t1 * t42 + t2 * t43;
t46 = t34 * pkin(3);
t26 = pkin(4) + t46;
t47 = t33 * pkin(3);
t13 = t37 * t26 - t35 * t47;
t14 = t35 * t26 + t37 * t47;
t52 = t13 * t42 + t14 * t43;
t27 = t36 * pkin(3) + qJ(2);
t12 = t18 * pkin(4) + t27;
t50 = 0.2e1 * t12;
t49 = 0.2e1 * t27;
t48 = 0.2e1 * qJ(2);
t30 = t36 ^ 2;
t31 = t38 ^ 2;
t24 = t30 + t31;
t45 = t10 * t20 + t11 * t18;
t40 = qJ(2) ^ 2;
t21 = t24 * t39;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t48, (pkin(1) ^ 2) + t40, t31, -0.2e1 * t38 * t36, 0, t30, 0, 0, t36 * t48, t38 * t48, -0.2e1 * t21, t24 * t39 ^ 2 + t40, t17, -0.2e1 * t20 * t18, 0, t51, 0, 0, t18 * t49, t20 * t49, -0.2e1 * t45, t10 ^ 2 + t11 ^ 2 + t27 ^ 2, t58, -0.2e1 * t42 * t43, 0, t57, 0, 0, t43 * t50, t42 * t50, -0.2e1 * t53, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t24, t21, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t45, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t36, 0, t28, -t36 * t39, 0, 0, 0, 0, t20, 0, -t18, 0, t10, -t11, -t56, (t10 * t34 + t11 * t33) * pkin(3), 0, 0, t42, 0, -t43, 0, t1, -t2, -t52, t1 * t13 + t2 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, t56, 0, 0, 0, 0, 0, 0, t42, -t43, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t47, 0, (t33 ^ 2 + t34 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t13, -0.2e1 * t14, 0, t13 ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t20, 0, t27, 0, 0, 0, 0, 0, 0, t43, t42, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, -t43, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t43, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t13, -t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
