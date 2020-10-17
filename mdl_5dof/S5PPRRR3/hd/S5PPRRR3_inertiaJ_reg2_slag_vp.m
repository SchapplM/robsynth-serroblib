% Calculate inertial parameters regressor of joint inertia matrix for
% S5PPRRR3
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:59
% EndTime: 2019-12-05 15:17:02
% DurationCPUTime: 0.45s
% Computational Cost: add. (186->67), mult. (484->116), div. (0->0), fcn. (580->8), ass. (0->43)
t30 = cos(qJ(4));
t17 = -t30 * pkin(4) - pkin(3);
t47 = 0.2e1 * t17;
t46 = 0.2e1 * t30;
t45 = -pkin(7) - pkin(6);
t26 = sin(qJ(5));
t44 = t26 * pkin(4);
t29 = cos(qJ(5));
t43 = t29 * pkin(4);
t27 = sin(qJ(4));
t28 = sin(qJ(3));
t42 = t27 * t28;
t24 = sin(pkin(9));
t41 = t28 * t24;
t40 = t30 * t28;
t31 = cos(qJ(3));
t39 = t31 * t24;
t38 = t31 * t27;
t37 = t31 * t30;
t20 = t27 ^ 2;
t22 = t30 ^ 2;
t36 = t20 + t22;
t35 = t36 * t28;
t25 = cos(pkin(9));
t7 = -t24 * t38 - t25 * t30;
t8 = t24 * t37 - t25 * t27;
t34 = -t7 * t27 + t8 * t30;
t11 = t26 * t30 + t29 * t27;
t23 = t31 ^ 2;
t21 = t28 ^ 2;
t19 = t25 ^ 2;
t18 = t24 ^ 2;
t15 = t21 * t18;
t13 = t45 * t30;
t12 = t45 * t27;
t9 = t26 * t27 - t29 * t30;
t6 = -t26 * t42 + t29 * t40;
t5 = t11 * t28;
t4 = t26 * t12 - t29 * t13;
t3 = t29 * t12 + t26 * t13;
t2 = t26 * t7 + t29 * t8;
t1 = -t26 * t8 + t29 * t7;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 + t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t18 + t15 + t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t34 - t39) * t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t5 + t2 * t6 - t28 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 + t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t21 + t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t39, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t40, t27 * t41, t34, -pkin(3) * t41 + pkin(6) * t34, 0, 0, 0, 0, 0, 0, t9 * t41, t11 * t41, -t1 * t11 - t2 * t9, t1 * t3 + t17 * t41 + t2 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t28, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t38, t35, t31 * pkin(3) + pkin(6) * t35, 0, 0, 0, 0, 0, 0, -t31 * t9, -t31 * t11, t5 * t11 - t6 * t9, -t31 * t17 - t5 * t3 + t6 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t20, t27 * t46, 0, t22, 0, 0, pkin(3) * t46, -0.2e1 * pkin(3) * t27, 0.2e1 * t36 * pkin(6), t36 * pkin(6) ^ 2 + pkin(3) ^ 2, t11 ^ 2, -0.2e1 * t11 * t9, 0, t9 ^ 2, 0, 0, t9 * t47, t11 * t47, -0.2e1 * t3 * t11 - 0.2e1 * t4 * t9, t17 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, (t1 * t29 + t2 * t26) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t40, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, (t26 * t6 - t29 * t5) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t30, 0, -t27 * pkin(6), -t30 * pkin(6), 0, 0, 0, 0, t11, 0, -t9, 0, t3, -t4, (-t11 * t29 - t26 * t9) * pkin(4), (t26 * t4 + t29 * t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t44, 0, (t26 ^ 2 + t29 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t43, -t44, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
