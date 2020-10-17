% Calculate inertial parameters regressor of joint inertia matrix for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:16
% EndTime: 2019-12-05 15:11:18
% DurationCPUTime: 0.41s
% Computational Cost: add. (87->42), mult. (246->62), div. (0->0), fcn. (266->6), ass. (0->41)
t28 = cos(qJ(4));
t26 = sin(qJ(4));
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t29 = cos(qJ(3));
t41 = t29 * t26;
t6 = t24 * t41 + t25 * t28;
t47 = t6 * t26;
t15 = t29 * t28;
t8 = t24 * t15 - t25 * t26;
t35 = t8 * t28 + t47;
t27 = sin(qJ(3));
t20 = t26 ^ 2;
t22 = t28 ^ 2;
t38 = t20 + t22;
t9 = t38 * t27;
t51 = -0.2e1 * t26;
t49 = t26 * pkin(6);
t48 = t28 * pkin(6);
t46 = t26 * t27;
t45 = t26 * t28;
t44 = t27 * t24;
t43 = t28 * t27;
t42 = t29 * t24;
t40 = pkin(6) * t9;
t39 = t38 * pkin(6) ^ 2;
t37 = t26 * t44;
t36 = t24 * t43;
t18 = t24 ^ 2;
t21 = t27 ^ 2;
t14 = t21 * t18;
t34 = t6 ^ 2 + t8 ^ 2 + t14;
t33 = t35 * pkin(6);
t32 = -t26 * pkin(4) + t28 * qJ(5);
t31 = t8 * t43 + (-t42 + t47) * t27;
t23 = t29 ^ 2;
t19 = t25 ^ 2;
t11 = -t28 * pkin(4) - t26 * qJ(5) - pkin(3);
t10 = 0.2e1 * t38 * pkin(6);
t5 = t38 * t21 + t23;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 + t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t18 + t14 + t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 + t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t42, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t37, t35, -pkin(3) * t44 + t33, 0, 0, 0, 0, 0, 0, -t36, t35, -t37, t11 * t44 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t27, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t41, t9, t29 * pkin(3) + t40, 0, 0, 0, 0, 0, 0, t15, t9, t41, -t29 * t11 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t20, 0.2e1 * t45, 0, t22, 0, 0, 0.2e1 * pkin(3) * t28, pkin(3) * t51, t10, pkin(3) ^ 2 + t39, t20, 0, -0.2e1 * t45, 0, 0, t22, -0.2e1 * t11 * t28, t10, t11 * t51, t11 ^ 2 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t8, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, t8, -t6 * pkin(4) + t8 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t43, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, t43, t32 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t28, 0, -t49, -t48, 0, 0, 0, t26, 0, 0, -t28, 0, -t49, t32, t48, t32 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
