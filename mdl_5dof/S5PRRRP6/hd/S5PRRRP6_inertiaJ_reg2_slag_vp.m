% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:17
% EndTime: 2019-12-05 16:52:21
% DurationCPUTime: 0.59s
% Computational Cost: add. (236->67), mult. (542->109), div. (0->0), fcn. (595->6), ass. (0->44)
t30 = sin(qJ(4));
t31 = sin(qJ(3));
t33 = cos(qJ(4));
t34 = cos(qJ(3));
t14 = t30 * t31 - t33 * t34;
t57 = t14 ^ 2;
t24 = -t34 * pkin(3) - pkin(2);
t56 = 0.2e1 * t24;
t55 = 0.2e1 * t34;
t54 = -pkin(7) - pkin(6);
t53 = t33 * pkin(3);
t16 = t30 * t34 + t33 * t31;
t52 = t16 * t14;
t32 = sin(qJ(2));
t51 = t31 * t32;
t50 = t34 * t32;
t35 = cos(qJ(2));
t49 = t35 * t14;
t48 = t35 * t16;
t26 = t31 ^ 2;
t28 = t34 ^ 2;
t47 = t26 + t28;
t18 = t54 * t34;
t45 = t54 * t31;
t6 = -t30 * t18 - t33 * t45;
t8 = -t33 * t18 + t30 * t45;
t46 = t6 ^ 2 + t8 ^ 2;
t10 = t16 * t32;
t12 = -t30 * t51 + t33 * t50;
t44 = t10 * t6 + t12 * t8;
t43 = t10 * t16 - t12 * t14;
t29 = t35 ^ 2;
t42 = t10 ^ 2 + t12 ^ 2 + t29;
t41 = t47 * t32;
t40 = -0.2e1 * t8 * t14 + 0.2e1 * t6 * t16;
t37 = 2 * pkin(4);
t36 = 2 * qJ(5);
t27 = t32 ^ 2;
t25 = t30 * pkin(3);
t22 = pkin(4) + t53;
t20 = t25 + qJ(5);
t13 = t16 ^ 2;
t3 = t14 * pkin(4) - t16 * qJ(5) + t24;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 + t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t27 + t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t34, -t35 * t31, t41, t35 * pkin(2) + pkin(6) * t41, 0, 0, 0, 0, 0, 0, -t49, -t48, t43, -t35 * t24 + t44, 0, 0, 0, 0, 0, 0, -t49, t43, t48, -t35 * t3 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t26, t31 * t55, 0, t28, 0, 0, pkin(2) * t55, -0.2e1 * pkin(2) * t31, 0.2e1 * t47 * pkin(6), t47 * pkin(6) ^ 2 + pkin(2) ^ 2, t13, -0.2e1 * t52, 0, t57, 0, 0, t14 * t56, t16 * t56, t40, t24 ^ 2 + t46, t13, 0, 0.2e1 * t52, 0, 0, t57, 0.2e1 * t3 * t14, t40, -0.2e1 * t3 * t16, t3 ^ 2 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, 0, (-t10 * t33 + t12 * t30) * pkin(3), 0, 0, 0, 0, 0, 0, -t10, 0, t12, -t10 * t22 + t12 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t34, 0, -t31 * pkin(6), -t34 * pkin(6), 0, 0, 0, 0, t16, 0, -t14, 0, -t6, -t8, (-t14 * t30 - t16 * t33) * pkin(3), (t30 * t8 - t33 * t6) * pkin(3), 0, t16, 0, 0, t14, 0, -t6, -t20 * t14 - t22 * t16, t8, t8 * t20 - t6 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t53, -0.2e1 * t25, 0, (t30 ^ 2 + t33 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t22, 0, 0.2e1 * t20, t20 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, t12, -t10 * pkin(4) + t12 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, -t6, -t8, 0, 0, 0, t16, 0, 0, t14, 0, -t6, -pkin(4) * t16 - t14 * qJ(5), t8, -t6 * pkin(4) + t8 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t53, -t25, 0, 0, 0, 0, 0, 1, 0, 0, t37 + t53, 0, t36 + t25, t22 * pkin(4) + t20 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t37, 0, t36, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
