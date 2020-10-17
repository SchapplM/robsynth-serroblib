% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:04
% EndTime: 2019-12-31 17:21:06
% DurationCPUTime: 0.38s
% Computational Cost: add. (175->62), mult. (418->122), div. (0->0), fcn. (365->4), ass. (0->51)
t26 = sin(qJ(3));
t22 = t26 ^ 2;
t28 = cos(qJ(3));
t24 = t28 ^ 2;
t55 = t22 + t24;
t33 = -t28 * pkin(3) - t26 * qJ(4);
t10 = -pkin(2) + t33;
t54 = -0.2e1 * t10;
t27 = sin(qJ(2));
t53 = -0.2e1 * t27;
t52 = 0.2e1 * t27;
t51 = pkin(2) * t26;
t50 = pkin(2) * t28;
t49 = pkin(5) * t26;
t23 = t27 ^ 2;
t48 = t23 * pkin(5);
t47 = t26 * pkin(6);
t46 = t27 * pkin(5);
t45 = t28 * pkin(6);
t29 = cos(qJ(2));
t11 = -t29 * pkin(2) - t27 * pkin(6) - pkin(1);
t40 = t28 * t29;
t4 = pkin(5) * t40 + t26 * t11;
t44 = t26 * t27;
t43 = t26 * t28;
t42 = t26 * t29;
t41 = t27 * t29;
t18 = t28 * t27;
t39 = t55 * pkin(6) ^ 2;
t38 = t29 * qJ(4);
t37 = t26 * t41;
t36 = t23 * t43;
t1 = -t38 + t4;
t7 = t28 * t11;
t2 = -t7 + (pkin(3) + t49) * t29;
t35 = t1 * t28 + t2 * t26;
t3 = -pkin(5) * t42 + t7;
t34 = -t3 * t26 + t4 * t28;
t32 = -pkin(3) * t26 + t28 * qJ(4);
t31 = pkin(5) ^ 2;
t25 = t29 ^ 2;
t20 = t23 * t31;
t17 = t24 * t23;
t16 = t22 * t23;
t14 = pkin(6) * t42;
t13 = t26 * t18;
t12 = t40 * t53;
t9 = 0.2e1 * t55 * pkin(6);
t8 = (t22 - t24) * t27;
t5 = (pkin(5) - t32) * t27;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t23, 0.2e1 * t41, 0, t25, 0, 0, 0.2e1 * pkin(1) * t29, pkin(1) * t53, 0.2e1 * (t23 + t25) * pkin(5), pkin(1) ^ 2 + t25 * t31 + t20, t17, -0.2e1 * t36, t12, t16, 0.2e1 * t37, t25, 0.2e1 * t26 * t48 - 0.2e1 * t3 * t29, 0.2e1 * t28 * t48 + 0.2e1 * t4 * t29, (-t26 * t4 - t28 * t3) * t52, t3 ^ 2 + t4 ^ 2 + t20, t17, t12, 0.2e1 * t36, t25, -0.2e1 * t37, t16, 0.2e1 * t2 * t29 + 0.2e1 * t5 * t44, (-t1 * t26 + t2 * t28) * t52, -0.2e1 * t1 * t29 - 0.2e1 * t5 * t18, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t29, 0, -t46, -t29 * pkin(5), 0, 0, t13, -t8, -t42, -t13, -t40, 0, t14 + (-pkin(5) * t28 - t51) * t27, pkin(6) * t40 + (t49 - t50) * t27, t34, -pkin(2) * t46 + t34 * pkin(6), t13, -t42, t8, 0, t40, -t13, t10 * t44 - t5 * t28 + t14, t35, -t5 * t26 + (-pkin(6) * t29 - t10 * t27) * t28, t35 * pkin(6) + t5 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t22, 0.2e1 * t43, 0, t24, 0, 0, 0.2e1 * t50, -0.2e1 * t51, t9, pkin(2) ^ 2 + t39, t22, 0, -0.2e1 * t43, 0, 0, t24, t28 * t54, t9, t26 * t54, t10 ^ 2 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t44, -t29, t3, -t4, 0, 0, 0, t18, 0, -t29, t44, 0, t7 + (-0.2e1 * pkin(3) - t49) * t29, t33 * t27, -0.2e1 * t38 + t4, -t2 * pkin(3) + t1 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t28, 0, -t47, -t45, 0, 0, 0, t26, 0, 0, -t28, 0, -t47, t32, t45, t32 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t18, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
