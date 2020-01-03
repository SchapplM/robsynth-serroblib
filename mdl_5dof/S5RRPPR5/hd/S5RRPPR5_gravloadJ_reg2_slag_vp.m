% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t16 = g(1) * t34 + g(2) * t31;
t27 = qJ(2) + pkin(8);
t23 = sin(t27);
t54 = t16 * t23;
t19 = t23 * qJ(4);
t24 = cos(t27);
t53 = -t24 * pkin(3) - t19;
t30 = sin(qJ(2));
t51 = pkin(2) * t30;
t48 = t24 * pkin(4);
t28 = -qJ(3) - pkin(6);
t47 = pkin(7) + t28;
t32 = cos(qJ(5));
t46 = t23 * t32;
t45 = t24 * t34;
t44 = t34 * t28;
t43 = qJ(4) * t24;
t33 = cos(qJ(2));
t25 = t33 * pkin(2);
t42 = t25 - t53;
t22 = t25 + pkin(1);
t18 = t34 * t22;
t41 = g(2) * (pkin(3) * t45 + t34 * t19 + t18);
t40 = -pkin(3) * t23 - t51;
t15 = g(1) * t31 - g(2) * t34;
t29 = sin(qJ(5));
t10 = -t24 * t29 + t46;
t39 = t23 * t29 + t24 * t32;
t38 = -t22 + t53;
t3 = t10 * t31;
t5 = t29 * t45 - t34 * t46;
t37 = g(1) * t5 - g(2) * t3 + g(3) * t39;
t4 = t39 * t31;
t6 = t39 * t34;
t36 = g(1) * t6 + g(2) * t4 + g(3) * t10;
t35 = -g(3) * t33 + t16 * t30;
t14 = t34 * t43;
t12 = t31 * t43;
t8 = t15 * t24;
t7 = t15 * t23;
t2 = g(3) * t23 + t16 * t24;
t1 = -g(3) * t24 + t54;
t9 = [0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t33, -t15 * t30, -t16, -g(1) * (-t31 * pkin(1) + t34 * pkin(6)) - g(2) * (t34 * pkin(1) + t31 * pkin(6)), 0, 0, 0, 0, 0, 0, t8, -t7, -t16, -g(1) * (-t31 * t22 - t44) - g(2) * (-t31 * t28 + t18), 0, 0, 0, 0, 0, 0, t8, -t16, t7, g(1) * t44 - t41 + (-g(1) * t38 + g(2) * t28) * t31, 0, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, g(1) * t3 + g(2) * t5, t16, -t41 + (g(1) * t47 - g(2) * t48) * t34 + (-g(1) * (t38 - t48) + g(2) * t47) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, g(3) * t30 + t16 * t33, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t35 * pkin(2), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t40 * t34 + t14) - g(2) * (t40 * t31 + t12) - g(3) * t42, 0, 0, 0, 0, 0, 0, -t37, -t36, 0, -g(1) * (-t34 * t51 + t14) - g(2) * (-t31 * t51 + t12) - g(3) * (t42 + t48) + (pkin(3) + pkin(4)) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t36, 0, 0;];
taug_reg = t9;
