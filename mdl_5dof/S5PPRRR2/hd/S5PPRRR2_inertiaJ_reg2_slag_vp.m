% Calculate inertial parameters regressor of joint inertia matrix for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:49
% EndTime: 2019-12-05 15:14:51
% DurationCPUTime: 0.40s
% Computational Cost: add. (177->52), mult. (420->91), div. (0->0), fcn. (516->8), ass. (0->35)
t23 = sin(pkin(9));
t24 = cos(pkin(9));
t27 = sin(qJ(3));
t36 = cos(qJ(3));
t8 = t27 * t23 - t36 * t24;
t6 = t8 ^ 2;
t29 = cos(qJ(4));
t20 = -t29 * pkin(4) - pkin(3);
t41 = 0.2e1 * t20;
t40 = 0.2e1 * t29;
t39 = -pkin(7) - pkin(6);
t25 = sin(qJ(5));
t38 = t25 * pkin(4);
t28 = cos(qJ(5));
t37 = t28 * pkin(4);
t10 = t36 * t23 + t27 * t24;
t26 = sin(qJ(4));
t35 = t26 * t10;
t34 = t29 * t10;
t21 = t26 ^ 2;
t22 = t29 ^ 2;
t33 = t21 + t22;
t32 = t33 * t10;
t15 = t25 * t29 + t28 * t26;
t17 = t39 * t29;
t16 = t39 * t26;
t13 = t25 * t26 - t28 * t29;
t12 = t15 ^ 2;
t11 = t13 ^ 2;
t7 = t10 ^ 2;
t4 = t25 * t16 - t28 * t17;
t3 = t28 * t16 + t25 * t17;
t2 = -t25 * t35 + t28 * t34;
t1 = t15 * t10;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 ^ 2 + t24 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 + t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t7 + t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t13 + t2 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t10, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t29, t8 * t26, t32, -t8 * pkin(3) + pkin(6) * t32, 0, 0, 0, 0, 0, 0, t8 * t13, t8 * t15, t1 * t15 - t2 * t13, -t1 * t3 + t2 * t4 + t8 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t3 + t15 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t21, t26 * t40, 0, t22, 0, 0, pkin(3) * t40, -0.2e1 * pkin(3) * t26, 0.2e1 * t33 * pkin(6), t33 * pkin(6) ^ 2 + pkin(3) ^ 2, t12, -0.2e1 * t15 * t13, 0, t11, 0, 0, t13 * t41, t15 * t41, -0.2e1 * t4 * t13 - 0.2e1 * t3 * t15, t20 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t34, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, (-t1 * t28 + t2 * t25) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t26, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, (-t13 * t28 + t15 * t25) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t29, 0, -t26 * pkin(6), -t29 * pkin(6), 0, 0, 0, 0, t15, 0, -t13, 0, t3, -t4, (-t13 * t25 - t15 * t28) * pkin(4), (t25 * t4 + t28 * t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t38, 0, (t25 ^ 2 + t28 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t13, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
