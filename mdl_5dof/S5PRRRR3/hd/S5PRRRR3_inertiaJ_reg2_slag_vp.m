% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:18
% EndTime: 2019-12-05 17:06:20
% DurationCPUTime: 0.46s
% Computational Cost: add. (158->41), mult. (331->71), div. (0->0), fcn. (259->6), ass. (0->32)
t25 = sin(qJ(5));
t23 = t25 ^ 2;
t28 = cos(qJ(5));
t24 = t28 ^ 2;
t37 = t23 + t24;
t30 = cos(qJ(3));
t21 = t30 * pkin(2);
t17 = t21 + pkin(3);
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t27 = sin(qJ(3));
t42 = t27 * pkin(2);
t36 = t29 * t42;
t8 = t26 * t17 + t36;
t6 = pkin(8) + t8;
t45 = t37 * t6;
t43 = t26 * pkin(3);
t15 = pkin(8) + t43;
t46 = t37 * t15;
t44 = pkin(4) * t25;
t34 = -t29 * t17 + t26 * t42;
t5 = -pkin(4) + t34;
t41 = t5 * t28;
t20 = t29 * pkin(3);
t16 = -t20 - pkin(4);
t39 = t16 * t28;
t38 = pkin(8) * t37;
t22 = pkin(4) * t28;
t13 = 0.2e1 * t25 * t28;
t12 = t16 * t25;
t3 = t5 * t25;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t42, 0, (t27 ^ 2 + t30 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t34, -0.2e1 * t8, 0, t34 ^ 2 + t8 ^ 2, t23, t13, 0, t24, 0, 0, -0.2e1 * t41, 0.2e1 * t3, 0.2e1 * t45, t37 * t6 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, -t42, 0, 0, 0, 0, 0, 0, 0, 1, t20 - t34, -t36 + (-pkin(3) - t17) * t26, 0, (t26 * t8 - t29 * t34) * pkin(3), t23, t13, 0, t24, 0, 0, (-t16 - t5) * t28, t12 + t3, t46 + t45, t5 * t16 + t46 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t43, 0, (t26 ^ 2 + t29 ^ 2) * pkin(3) ^ 2, t23, t13, 0, t24, 0, 0, -0.2e1 * t39, 0.2e1 * t12, 0.2e1 * t46, t15 ^ 2 * t37 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t34, -t8, 0, 0, t23, t13, 0, t24, 0, 0, t22 - t41, t3 - t44, t38 + t45, -t5 * pkin(4) + pkin(8) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t20, -t43, 0, 0, t23, t13, 0, t24, 0, 0, t22 - t39, t12 - t44, t38 + t46, -t16 * pkin(4) + pkin(8) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t23, t13, 0, t24, 0, 0, 0.2e1 * t22, -0.2e1 * t44, 0.2e1 * t38, pkin(8) ^ 2 * t37 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * t6, -t28 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * t15, -t28 * t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * pkin(8), -t28 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
